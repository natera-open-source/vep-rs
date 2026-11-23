// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The GWAS and SpliceAI classes and what they read: the pairs of a catalog's associations
//! file, the site a row states, the records of a score file at a span, the rows each class
//! takes and the grounds on which `spliceai_multi_entry_cell` refuses the others, on made-up
//! rows, catalogs and score files; then the `score` command over the same shapes and over
//! the committed SpliceAI corpora. The expected values of the class tests are those
//! `tools/vep-formats-compare` gives on the same inputs.

use flate2::read::MultiGzDecoder;
use flate2::write::{DeflateEncoder, GzEncoder};
use flate2::{Compression, Crc};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;
use vep_plugin_pairs::registry::BuiltInRegistry;
use vep_plugin_pairs::rules::entry::{entry_kind, parse_default, Entry, Row, ABSENT};
use vep_plugin_pairs::rules::gwas::{gwas_site, GwasCatalogPairs, Pair};
use vep_plugin_pairs::rules::spliceai::{
    refusal_lines, RefusalCounts, ScoreRecord, SpliceAiTable, GROUNDS,
};
use vep_plugin_pairs::rules::{BuiltIn, Consult, RecordRows, Refusal, Rules};
use vep_plugin_pairs::score::run;

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";
const CROSS_PRODUCT: &str = "gwas_haplotype_cross_product";
const INTERACTION: &str = "gwas_interaction_first_snp_allele_whitespace";
const CHOPPED: &str = "spliceai_chopped_multiallelic";
const MULTI: &str = "spliceai_multi_entry_cell";

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-gwas-spliceai-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    fn file(&self, name: &str, content: impl AsRef<[u8]>) -> PathBuf {
        let path = self.0.join(name);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(&path, content).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn text_of(path: &Path) -> &str {
    path.to_str().unwrap()
}

fn refusal<T>(read: Result<T, Refusal>) -> String {
    match read {
        Ok(_) => panic!("the file was read"),
        Err(Refusal(message)) => message,
    }
}

// ---- the default format's rows and what a class makes of two outputs -------------------

fn default_row(name: &str, location: &str, allele: &str, feature: &str, extra: &str) -> String {
    format!("{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\t{extra}")
}

fn default_text(rows: &[String]) -> String {
    format!("{HEADER}\n{}\n", rows.join("\n"))
}

/// An entry of the VCF or JSON format, or of the tab format where `location` is given: the
/// fields a parser of that format files a row or a record under.
fn entry(
    key: &str,
    location: Option<&str>,
    allele: &str,
    fields: &[(&str, &str)],
    indel: Option<bool>,
    allele_location: Option<&str>,
) -> Row {
    Rc::new(Entry {
        key: key.to_string(),
        record: location.is_none().then_some(0),
        location: location.map(str::to_string),
        allele: allele.to_string(),
        feature: "ENST1".to_string(),
        fields: fields.iter().copied().collect(),
        indel,
        allele_location: allele_location.map(str::to_string),
        ..Entry::default()
    })
}

/// What a run hands the class under test.
#[derive(Default, Clone, Copy)]
struct Readers<'a> {
    gwas: Option<&'a GwasCatalogPairs>,
    spliceai: Option<&'a SpliceAiTable>,
    input: &'a [&'a str],
}

/// What a class makes of the rows two outputs write under one key: per field the rows on
/// which it accepts the differing value, the rows it accepts a field of, and the differing
/// values it leaves.
#[derive(Debug, Default)]
struct Compared {
    masked: BTreeMap<String, usize>,
    rows: usize,
    left: usize,
}

impl Compared {
    /// Whether no difference is left and the class took a field of `rows` rows.
    fn green(&self, rows: usize) -> bool {
        self.left == 0 && self.rows == rows
    }

    /// Whether a difference is left and the class took a field of `rows` rows.
    fn red(&self, rows: usize) -> bool {
        self.left > 0 && self.rows == rows
    }
}

fn compare_entries(expected: &[Row], actual: &[Row], klass: &str, readers: Readers) -> Compared {
    let index = RecordRows::new(expected, readers.input);
    let declared: BTreeSet<String> = [klass.to_string()].into();
    let mut compared = Compared::default();
    for e in expected.iter().filter(|e| entry_kind(e) != "record") {
        let Some(a) = actual.iter().find(|a| a.key == e.key) else {
            continue;
        };
        let record_rows = index.of(e).unwrap();
        let names: Vec<&str> = e
            .fields
            .iter()
            .map(|(name, _)| name)
            .chain(
                a.fields
                    .iter()
                    .map(|(name, _)| name)
                    .filter(|name| !e.fields.contains(name)),
            )
            .filter(|name| *name != "Extra")
            .collect();
        let mut taken = false;
        for name in names {
            let expected_value = e.fields.get(name).unwrap_or(ABSENT);
            let actual_value = a.fields.get(name).unwrap_or(ABSENT);
            if expected_value == actual_value {
                continue;
            }
            let accepted = BuiltIn
                .accepts(&Consult {
                    klass,
                    e,
                    name,
                    expected: expected_value,
                    actual: actual_value,
                    a: Some(a),
                    record_rows: record_rows.as_deref(),
                    candidate_rows: None,
                    declared: &declared,
                    assembly: None,
                    source: None,
                    gwas_pairs: readers.gwas,
                    revel_table: None,
                    utrannotator_sweep: None,
                    downstream_sequences: None,
                    splice_windows: None,
                    spliceai_table: readers.spliceai,
                    exon_spans: None,
                })
                .unwrap();
            if accepted {
                *compared.masked.entry(name.to_string()).or_default() += 1;
                taken = true;
            } else {
                compared.left += 1;
            }
        }
        compared.rows += usize::from(taken);
    }
    compared
}

fn compare(reference: &str, engine: &str, klass: &str, readers: Readers) -> Compared {
    compare_entries(
        &parse_default(reference).unwrap(),
        &parse_default(engine).unwrap(),
        klass,
        readers,
    )
}

// ---- the GWAS catalog and its two classes ----------------------------------------------

const CATALOG_HEADER: &str = "CHR_ID\tCHR_POS\tSTRONGEST SNP-RISK ALLELE\tSNPS\tSTUDY ACCESSION";
const GWAS_ROW_FIELDS: &str =
    "GWAS_associations=rs1;GWAS_risk_allele=rs1-A;GWAS_p_value=1e-8;GWAS_study=GCST1";

/// The pairs of a catalog of `rows`: the chromosomes, the positions, the risk alleles and
/// the study of each.
fn gwas_pairs(scratch: &Scratch, rows: &[[&str; 4]]) -> GwasCatalogPairs {
    let lines: Vec<String> = rows
        .iter()
        .map(|[chromosomes, positions, alleles, study]| {
            format!("{chromosomes}\t{positions}\t{alleles}\t-\t{study}")
        })
        .collect();
    let catalog = scratch.file(
        "gwas_catalog.tsv",
        format!("{CATALOG_HEADER}\n{}\n", lines.join("\n")),
    );
    GwasCatalogPairs::from_catalog(&catalog).unwrap()
}

fn pairs(held: &HashSet<Pair>) -> BTreeSet<[&str; 4]> {
    held.iter()
        .map(|(chromosome, position, allele, study)| {
            [
                chromosome.as_str(),
                position.as_str(),
                allele.as_str(),
                study.as_str(),
            ]
        })
        .collect()
}

fn gwas_row(name: &str, location: &str, allele: &str, fields: &str) -> String {
    let more = if fields.is_empty() {
        String::new()
    } else {
        format!(";{fields}")
    };
    default_row(
        name,
        location,
        allele,
        "ENST1",
        &format!("IMPACT=MODIFIER;STRAND=1{more}"),
    )
}

fn with_catalog(pairs: &GwasCatalogPairs) -> Readers<'_> {
    Readers {
        gwas: Some(pairs),
        ..Readers::default()
    }
}

/// The record entry and the row a VCF parser files for one line whose CSQ names the
/// `GWAS_study` field alone.
fn vcf_entries(
    position: &str,
    reference: &str,
    alternates: &str,
    allele: &str,
    indel: bool,
    allele_location: &str,
    study: &str,
) -> Vec<Row> {
    vec![
        entry(
            "rec0|record",
            None,
            "",
            &[
                ("CHROM", "21"),
                ("POS", position),
                ("ID", "named"),
                ("REF", reference),
                ("ALT", alternates),
            ],
            None,
            None,
        ),
        entry(
            &format!("rec0|{allele}|ENST1|Transcript"),
            None,
            allele,
            &[
                ("Allele", allele),
                ("Consequence", "intron_variant"),
                ("Feature", "ENST1"),
                ("Feature_type", "Transcript"),
                ("GWAS_study", study),
            ],
            Some(indel),
            Some(allele_location),
        ),
    ]
}

/// The record entry and the row a JSON parser files for one object: the record's
/// `seq_region_name` and span, and the row's `gwas_` members.
fn json_entries(
    span: &[(&str, &str)],
    allele: &str,
    indel: Option<bool>,
    members: &[(&str, &str)],
) -> Vec<Row> {
    let record: Vec<(&str, &str)> = [("id", "named"), ("seq_region_name", "21")]
        .into_iter()
        .chain(span.iter().copied())
        .collect();
    let row: Vec<(&str, &str)> = [
        ("variant_allele", allele),
        ("transcript_id", "ENST1"),
        ("consequence_terms", "intron_variant"),
    ]
    .into_iter()
    .chain(members.iter().copied())
    .collect();
    vec![
        entry("rec0|record", None, "", &record, None, None),
        entry(
            &format!("rec0|{allele}|ENST1|transcript"),
            None,
            allele,
            &row,
            indel,
            None,
        ),
    ]
}

#[test]
fn gwas_haplotype_cross_product_accepts_reference_filled_fields_only() {
    let scratch = Scratch::new("cross-product");
    let held = gwas_pairs(
        &scratch,
        &[
            ["21;21", "100;150", "rs1-G; rs2-A", "GCST1"],
            ["21;21", "200;250", "rs3-T; rs4-C", "GCST1"],
        ],
    );
    let reference_only = pairs(held.reference_only());
    assert!(reference_only.contains(&["21", "100", "A", "GCST1"]));
    assert!(reference_only.contains(&["21", "150", "G", "GCST1"]));
    assert!(held.engine_only().is_empty());
    let second = GWAS_ROW_FIELDS.replace("rs1", "rs2");
    let reference = default_text(&[
        gwas_row("a", "21:100", "A", GWAS_ROW_FIELDS),
        gwas_row("b", "21:200", "C", &second),
        gwas_row("c", "21:300", "G", ""),
    ]);
    let other_study = second.replace("GCST1", "GCST2").replace("1e-8", "3e-9");
    let candidate = default_text(&[
        gwas_row("a", "21:100", "A", ""),
        gwas_row("b", "21:200", "C", &other_study),
        gwas_row("c", "21:300", "G", ""),
    ]);
    let compared = compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held));
    assert!(compared.green(2), "{compared:?}");
    assert_eq!(compared.masked["GWAS_study"], 2);
    assert_eq!(compared.masked["GWAS_risk_allele"], 1);
    // A field the engine alone fills stays a difference.
    let candidate_only = candidate.replace(
        &gwas_row("c", "21:300", "G", ""),
        &gwas_row("c", "21:300", "G", GWAS_ROW_FIELDS),
    );
    assert!(
        compare(
            &reference,
            &candidate_only,
            CROSS_PRODUCT,
            with_catalog(&held)
        )
        .left
            > 0
    );
    // A field of another plugin differing on the same row is not the class's.
    let other_plugin = candidate.replacen(
        "IMPACT=MODIFIER;STRAND=1\n",
        "IMPACT=MODIFIER;STRAND=1;LoFtool=0.5\n",
        1,
    );
    assert!(
        compare(
            &reference,
            &other_plugin,
            CROSS_PRODUCT,
            with_catalog(&held)
        )
        .left
            > 0
    );
    // The G row of the line 21 99 CA CG,C starts at 21:100, one base past the line's POS:
    // the VCF entry is read there, and a catalog site at the POS itself pairs nothing.
    let expected = vcf_entries("99", "CA", "CG,C", "G", false, "21:100", "GCST1");
    let actual = vcf_entries("99", "CA", "CG,C", "G", false, "21:100", "");
    for (position, rows) in [("100", 1), ("99", 0)] {
        let site = gwas_pairs(
            &scratch,
            &[["21;21", &format!("{position};150"), "rs1-T; rs2-G", "GCST1"]],
        );
        let compared = compare_entries(&expected, &actual, CROSS_PRODUCT, with_catalog(&site));
        assert_eq!(
            (compared.rows, compared.left),
            (rows, 1 - rows),
            "{position}"
        );
    }
}

#[test]
fn gwas_classes_read_no_insertions_row_at_the_records_start_in_any_format() {
    let scratch = Scratch::new("site");
    let held = gwas_pairs(&scratch, &[["21;21", "100;150", "rs1-G; rs2-A", "GCST1"]]);
    assert!(pairs(held.reference_only()).contains(&["21", "100", "A", "GCST1"]));
    let site = Some(("21".to_string(), "100".to_string()));
    for (name, location, taken) in [("sub", "21:100", true), ("ins", "21:99-100", false)] {
        let reference = default_text(&[gwas_row(name, location, "A", "GWAS_study=GCST1")]);
        let candidate = default_text(&[gwas_row(name, location, "A", "")]);
        let compared = compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held));
        assert_eq!(
            (compared.rows, compared.left),
            (usize::from(taken), usize::from(!taken)),
            "{name}"
        );
    }
    // The VCF format: a substitution at the POS, and an insertion whose span is empty.
    for (position, reference, alternate, indel, allele_location, taken) in [
        ("100", "C", "A", false, "21:100", true),
        ("99", "C", "CA", true, "21:99-100", false),
    ] {
        let entries = |study| {
            vcf_entries(
                position,
                reference,
                alternate,
                "A",
                indel,
                allele_location,
                study,
            )
        };
        let (expected, actual) = (entries("GCST1"), entries(""));
        let compared = compare_entries(&expected, &actual, CROSS_PRODUCT, with_catalog(&held));
        assert_eq!(compared.rows, usize::from(taken), "{position}");
        assert_eq!(
            gwas_site(&expected[1], Some(&expected)),
            site.clone().filter(|_| taken),
            "{position}"
        );
    }
    // The JSON format: the record object's start and end.
    for (span, indel, taken) in [
        (
            [("start", "100"), ("end", "100"), ("allele_string", "C/A")],
            false,
            true,
        ),
        (
            [("start", "100"), ("end", "99"), ("allele_string", "-/A")],
            true,
            false,
        ),
    ] {
        let expected = json_entries(&span, "A", Some(indel), &[("gwas_study", "GCST1")]);
        let actual = json_entries(&span, "A", Some(indel), &[]);
        let compared = compare_entries(&expected, &actual, CROSS_PRODUCT, with_catalog(&held));
        assert_eq!(compared.rows, usize::from(taken), "{span:?}");
        assert_eq!(
            gwas_site(&expected[1], Some(&expected)),
            site.clone().filter(|_| taken),
            "{span:?}"
        );
    }
}

#[test]
fn gwas_haplotype_cross_product_refuses_a_row_the_catalog_does_not_pair() {
    let scratch = Scratch::new("unpaired");
    let held = gwas_pairs(
        &scratch,
        &[
            ["21;21", "100;150", "rs1-G; rs2-A", "GCST1"],
            ["21", "400", "rs9-A", "GCST1"],
        ],
    );
    let reference = default_text(&[gwas_row("a", "21:100", "A", GWAS_ROW_FIELDS)]);
    let candidate = default_text(&[gwas_row("a", "21:100", "A", "")]);
    assert!(compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held)).green(1));
    assert!(compare(&reference, &candidate, CROSS_PRODUCT, Readers::default()).red(0));
    for (from, to) in [
        ("21:100", "21:400"),
        ("\tA\t", "\tG\t"),
        ("GCST1", "GCST7"),
        ("21:100", "21:100-101"),
    ] {
        let compared = compare(
            &reference.replace(from, to),
            &candidate.replace(from, to),
            CROSS_PRODUCT,
            with_catalog(&held),
        );
        assert!(compared.red(0), "{from} {compared:?}");
    }
    let ordinary = default_text(&[gwas_row("a", "21:400", "A", GWAS_ROW_FIELDS)]);
    for planted in [
        default_text(&[gwas_row("a", "21:400", "A", "")]),
        ordinary.replace("GCST1", "GCST2"),
    ] {
        let compared = compare(&ordinary, &planted, CROSS_PRODUCT, with_catalog(&held));
        assert!(compared.red(0), "{compared:?}");
    }
    // A catalog row whose chromosomes and positions differ in count pairs nothing,
    // whichever list is the longer.
    for lists in [
        ["21;21;21", "100;150", "rs1-G; rs2-A", "GCST1"],
        ["21;21", "100;150;175", "rs1-G; rs2-A", "GCST1"],
    ] {
        let uneven = gwas_pairs(&scratch, &[lists]);
        assert!(uneven.reference_only().is_empty() && uneven.engine_only().is_empty());
    }
}

const GWAS_SKIPPED_ROW_FIELDS: &str = "GWAS_accessions=http://purl.obolibrary.org/obo/MONDO_0005086;GWAS_associated_gene=MYC,PVT1;GWAS_beta_coef=;GWAS_odds_ratio=1.27;GWAS_p_value=5E-11;GWAS_pmid=24220699;GWAS_risk_allele=CG;GWAS_study=GCST002273";

#[test]
fn gwas_haplotype_cross_product_pairs_a_row_whose_snps_and_locations_differ_in_count() {
    let scratch = Scratch::new("skipped-row");
    let lists = "8;8\t127877125;127877126\trs35252396-CG; rs6470588-C; rs6470589-G";
    let catalog = |snps: &str| {
        let path = scratch.file(
            "gwas_catalog.tsv",
            format!("{CATALOG_HEADER}\n{lists}\t{snps}\tGCST002273\n"),
        );
        GwasCatalogPairs::from_catalog(&path).unwrap()
    };
    let held = catalog("rs35252396; rs6470588; rs6470589");
    let expected: BTreeSet<[&str; 4]> = ["127877125", "127877126"]
        .into_iter()
        .flat_map(|position| ["CG", "C", "G"].map(|allele| ["8", position, allele, "GCST002273"]))
        .collect();
    assert_eq!(pairs(held.reference_only()), expected);
    assert!(held.engine_only().is_empty());
    let planted = |location: &str, allele: &str, fields: &str, name: &str| {
        (
            default_text(&[gwas_row(name, location, allele, fields)]),
            default_text(&[gwas_row(name, location, allele, "")]),
        )
    };
    let fields = GWAS_SKIPPED_ROW_FIELDS;
    let (reference, candidate) = planted("8:127877125-127877126", "CG", fields, "plant05");
    let compared = compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held));
    assert!(compared.green(1), "{compared:?}");
    assert_eq!(compared.masked["GWAS_study"], 1);
    assert!(compare(&reference, &candidate, CROSS_PRODUCT, Readers::default()).red(0));
    for (location, allele) in [
        ("8:127877125", "CG"),
        ("8:127877126-127877127", "CG"),
        ("8:127877125", "G"),
    ] {
        let (reference, candidate) = planted(location, allele, fields, "plant05");
        let compared = compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held));
        assert!(compared.green(1), "{location} {allele} {compared:?}");
    }
    let (reversed_candidate, reversed_reference) =
        planted("8:127877125-127877126", "CG", fields, "plant05");
    for (reference, candidate) in [
        planted("8:127877125", "T", fields, "plant05"),
        planted("8:127877127", "G", fields, "plant05"),
        planted("8:127877124-127877125", "CG", fields, "plant05"),
        planted("8:127877125-127877127", "CG", fields, "plant05"),
        planted(
            "8:127877125-127877126",
            "CG",
            &fields.replace("GCST002273", "GCST1"),
            "plant05",
        ),
        planted("8:127877125-127877126", "CG", fields, "8_127877125_A/ACG"),
        (reversed_reference, reversed_candidate),
    ] {
        let compared = compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held));
        assert!(compared.red(0), "{reference} {compared:?}");
    }
    // As many SNPs as locations: the engine pairs each with its own, the reference with both.
    let even = gwas_pairs(
        &scratch,
        &[["8;8", "127877125;127877126", "rs1-CG; rs2-C", "GCST1"]],
    );
    assert_eq!(
        pairs(even.reference_only()),
        BTreeSet::from([
            ["8", "127877125", "C", "GCST1"],
            ["8", "127877126", "CG", "GCST1"]
        ])
    );
    let two_snps = catalog("rs35252396; rs6470588");
    assert!(two_snps.reference_only().is_empty() && two_snps.engine_only().is_empty());
    // GWAS.pm:261 reads a risk-allele cell holding ? as empty, so no SNP of the row is
    // written.
    for alleles in ["rs1-CG; rs2-?", "rs1-CG; rs2-C; rs3-?", "rs1-CG x rs2-?"] {
        let unknown = gwas_pairs(
            &scratch,
            &[["8;8", "127877125;127877126", alleles, "GCST1"]],
        );
        assert!(
            unknown.reference_only().is_empty() && unknown.engine_only().is_empty(),
            "{alleles}"
        );
    }
}

/// The fields of a row at a site two catalog rows name: the accession, the beta
/// coefficient, the p-value, the publication and the study.
fn site_fields(values: [&str; 5]) -> String {
    let [accession, beta, p_value, pmid, study] = values;
    format!("GWAS_accessions={accession};GWAS_associated_gene=;GWAS_beta_coef={beta};GWAS_odds_ratio=;GWAS_p_value={p_value};GWAS_pmid={pmid};GWAS_risk_allele=A;GWAS_study={study}")
}

fn haplotype_row() -> String {
    site_fields([
        "http://purl.obolibrary.org/obo/MONDO_9990001",
        "",
        "3E-18",
        "99000001",
        "GCST990001",
    ])
}

fn own_row() -> String {
    site_fields([
        "http://www.ebi.ac.uk/efo/EFO_9990002",
        "0.0125 unit increase",
        "4E-8",
        "99000002",
        "GCST990002",
    ])
}

/// The catalog rows that name 2:5000402: the last location of a haplotype row of study
/// GCST990001 and the one location of a single-SNP row of study GCST990002, whose risk
/// allele another SNP of the haplotype row carries.
fn own_row_pairs(scratch: &Scratch) -> GwasCatalogPairs {
    gwas_pairs(
        scratch,
        &[
            [
                "2;2",
                "5000149;5000402",
                "rs9990001-A; rs9990002-C",
                "GCST990001",
            ],
            ["2", "5000402", "rs9990002-A", "GCST990002"],
        ],
    )
}

fn own_row_texts(
    reference_fields: &str,
    candidate_fields: &str,
    location: &str,
    allele: &str,
) -> (String, String) {
    (
        default_text(&[gwas_row("gwas_own_row", location, allele, reference_fields)]),
        default_text(&[gwas_row("gwas_own_row", location, allele, candidate_fields)]),
    )
}

#[test]
fn gwas_haplotype_cross_product_accepts_a_field_the_sites_own_catalog_row_fills() {
    let scratch = Scratch::new("own-row");
    let held = own_row_pairs(&scratch);
    assert_eq!(
        pairs(held.reference_only()),
        BTreeSet::from([
            ["2", "5000149", "C", "GCST990001"],
            ["2", "5000402", "A", "GCST990001"]
        ])
    );
    let (reference, candidate) = own_row_texts(&haplotype_row(), &own_row(), "2:5000402", "A");
    let compared = compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held));
    assert!(compared.green(1), "{compared:?}");
    assert_eq!(compared.masked["GWAS_beta_coef"], 1);
    assert_eq!(compared.masked["GWAS_study"], 1);
}

#[test]
fn gwas_haplotype_cross_product_refuses_a_candidate_only_field_where_the_site_has_no_pair() {
    let scratch = Scratch::new("own-row-refused");
    let held = own_row_pairs(&scratch);
    let (haplotype, own) = (haplotype_row(), own_row());
    for ((reference, candidate), rows) in [
        (own_row_texts(&haplotype, &own, "2:5000403", "A"), 0),
        (own_row_texts(&haplotype, &own, "2:5000402", "G"), 0),
        (
            own_row_texts(
                &haplotype.replace("GCST990001", "GCST1"),
                &own,
                "2:5000402",
                "A",
            ),
            0,
        ),
        (own_row_texts("", &own, "2:5000402", "A"), 0),
        (
            own_row_texts(
                &haplotype,
                &own.replace("GCST990002", "GCST990001"),
                "2:5000402",
                "A",
            ),
            1,
        ),
        (
            own_row_texts(
                &haplotype,
                "GWAS_beta_coef=0.0125 unit increase",
                "2:5000402",
                "A",
            ),
            1,
        ),
        (
            own_row_texts(&haplotype, &format!("{own};LoFtool=0.5"), "2:5000402", "A"),
            1,
        ),
    ] {
        let compared = compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held));
        assert!(compared.red(rows), "{candidate} {compared:?}");
    }
    // Where the engine reports the reference's study, or the site has no pair, the beta
    // coefficient it alone fills is left.
    for (candidate_fields, location) in [
        (own.replace("GCST990002", "GCST990001"), "2:5000402"),
        (own.clone(), "2:5000403"),
    ] {
        let (reference, candidate) = own_row_texts(&haplotype, &candidate_fields, location, "A");
        let compared = compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held));
        assert!(
            !compared.masked.contains_key("GWAS_beta_coef"),
            "{location}"
        );
    }
}

#[test]
fn gwas_haplotype_cross_product_leaves_a_planted_wrong_or_missing_value_a_difference() {
    let scratch = Scratch::new("own-row-planted");
    let held = own_row_pairs(&scratch);
    let own = own_row();
    let unpaired = default_text(&[gwas_row("gwas_own_row", "2:5000500", "A", &own)]);
    for planted in [
        own.replace("0.0125", "0.5"),
        own.replace("GWAS_beta_coef=0.0125 unit increase", "GWAS_beta_coef="),
    ] {
        let candidate = default_text(&[gwas_row("gwas_own_row", "2:5000500", "A", &planted)]);
        let compared = compare(&unpaired, &candidate, CROSS_PRODUCT, with_catalog(&held));
        assert!(
            compared.red(0) && compared.left == 1,
            "{planted} {compared:?}"
        );
    }
}

#[test]
fn gwas_haplotype_cross_product_accepts_the_annotated_rows_empty_fields() {
    let scratch = Scratch::new("empty-fields");
    let held = gwas_pairs(&scratch, &[["21;21", "100;150", "rs1-A; rs2-G", "GCST1"]]);
    let filled = "GWAS_accessions=http://purl.obolibrary.org/obo/MONDO_0004992;GWAS_associated_gene=;GWAS_beta_coef=;GWAS_odds_ratio=0.41;GWAS_p_value=5E-7;GWAS_pmid=29299148;GWAS_risk_allele=G;GWAS_study=GCST1";
    let reference = default_text(&[gwas_row("a", "21:100", "G", filled)]);
    let candidate = default_text(&[gwas_row("a", "21:100", "G", "")]);
    let compared = compare(&reference, &candidate, CROSS_PRODUCT, with_catalog(&held));
    assert!(compared.green(1), "{compared:?}");
    assert_eq!(compared.masked["GWAS_beta_coef"], 1);
    // The engine alone fills the empty field: a difference.
    let filled_alone = candidate.replace("STRAND=1\n", "STRAND=1;GWAS_beta_coef=0.5\n");
    assert!(
        compare(
            &reference,
            &filled_alone,
            CROSS_PRODUCT,
            with_catalog(&held)
        )
        .left
            > 0
    );
    // The empty string is the reference's: a key the engine alone writes empty on an
    // annotated row stays a difference.
    let annotated_both = default_text(&[gwas_row("a", "21:100", "G", "GWAS_study=GCST1")]);
    let empty_key = annotated_both.replace("GWAS_study=GCST1", "GWAS_study=GCST1;GWAS_beta_coef=");
    assert!(compare(
        &annotated_both,
        &empty_key,
        CROSS_PRODUCT,
        with_catalog(&held)
    )
    .red(0));
    // Every GWAS_ field of the reference row empty: nothing annotated the row.
    let empty_only = default_text(&[gwas_row(
        "a",
        "21:100",
        "G",
        "GWAS_beta_coef=;GWAS_odds_ratio=",
    )]);
    assert!(compare(&empty_only, &candidate, CROSS_PRODUCT, with_catalog(&held)).red(0));
    // The tab format: an empty cell against the engine's `-`.
    let tab = |values: [&str; 3]| {
        vec![entry(
            "21:100|G|ENST1|Transcript",
            Some("21:100"),
            "G",
            &[
                ("Uploaded_variation", "a"),
                ("Location", "21:100"),
                ("Allele", "G"),
                ("GWAS_beta_coef", values[0]),
                ("GWAS_odds_ratio", values[1]),
                ("GWAS_study", values[2]),
            ],
            None,
            None,
        )]
    };
    let compared = compare_entries(
        &tab(["", "0.41", "GCST1"]),
        &tab(["-", "-", "-"]),
        CROSS_PRODUCT,
        with_catalog(&held),
    );
    assert!(
        compared.green(1) && compared.masked.len() == 3,
        "{compared:?}"
    );
    // The JSON format: an empty member against no key; the record object states the site.
    let members = [
        ("gwas_beta_coef", ""),
        ("gwas_odds_ratio", "0.41"),
        ("gwas_study", "GCST1"),
    ];
    let compared = compare_entries(
        &json_entries(&[("start", "100")], "G", None, &members),
        &json_entries(&[("start", "100")], "G", None, &[]),
        CROSS_PRODUCT,
        with_catalog(&held),
    );
    assert!(
        compared.green(1) && compared.masked.len() == 3,
        "{compared:?}"
    );
}

const INTERACTION_FIELDS: &str =
    "GWAS_p_value=3E-7;GWAS_pmid=29691896;GWAS_risk_allele=C;GWAS_study=GCST006494";

/// The catalog's interaction row of an SNP of chromosome 21 and one of chromosome 11, the
/// chromosome 21 SNP written first.
fn interaction_pairs(scratch: &Scratch) -> GwasCatalogPairs {
    gwas_pairs(
        scratch,
        &[[
            "21 x 11",
            "39375657 x 92975544",
            "rs4818032-C x rs10830963-G",
            "GCST006494",
        ]],
    )
}

#[test]
fn gwas_interaction_first_snp_allele_whitespace_accepts_the_first_snps_own_allele() {
    let scratch = Scratch::new("interaction");
    let held = interaction_pairs(&scratch);
    assert_eq!(
        pairs(held.engine_only()),
        BTreeSet::from([["21", "39375657", "C", "GCST006494"]])
    );
    assert_eq!(
        pairs(held.reference_only()),
        BTreeSet::from([["21", "39375657", "G", "GCST006494"]])
    );
    let row = |fields: &str| default_text(&[gwas_row("rs4818032", "21:39375657", "C", fields)]);
    let reference = row("");
    let compared = compare(
        &reference,
        &row(INTERACTION_FIELDS),
        INTERACTION,
        with_catalog(&held),
    );
    assert!(compared.green(1), "{compared:?}");
    assert_eq!(compared.masked["GWAS_study"], 1);
    // The engine's row carries the catalog row's empty cell as an empty field: the class's
    // with the filled ones, and a difference on a row the engine fills nowhere.
    let with_empty = row(&format!("GWAS_odds_ratio=;{INTERACTION_FIELDS}"));
    let compared = compare(&reference, &with_empty, INTERACTION, with_catalog(&held));
    assert!(compared.green(1) && compared.masked["GWAS_odds_ratio"] == 1);
    let empty_alone = row("GWAS_odds_ratio=");
    assert!(compare(&reference, &empty_alone, INTERACTION, with_catalog(&held)).red(0));
    // The JSON and VCF formats state the site on their record entry.
    let members = [("gwas_risk_allele", "C"), ("gwas_study", "GCST006494")];
    let compared = compare_entries(
        &json_entries(&[("start", "39375657")], "C", None, &[]),
        &json_entries(&[("start", "39375657")], "C", None, &members),
        INTERACTION,
        with_catalog(&held),
    );
    assert!(
        compared.green(1) && compared.masked.len() == 2,
        "{compared:?}"
    );
    let vcf = |study| vcf_entries("39375657", "T", "C", "C", false, "21:39375657", study);
    let compared = compare_entries(
        &vcf(""),
        &vcf("GCST006494"),
        INTERACTION,
        with_catalog(&held),
    );
    assert!(compared.green(1), "{compared:?}");
}

#[test]
fn gwas_interaction_first_snp_allele_whitespace_refuses_any_other_candidate_only_field() {
    let scratch = Scratch::new("interaction-refused");
    let held = interaction_pairs(&scratch);
    let row = |location: &str, fields: &str| {
        default_text(&[gwas_row("rs4818032", location, "C", fields)])
    };
    let reference = row("21:39375657", "");
    let candidate = row("21:39375657", INTERACTION_FIELDS);
    assert!(compare(&reference, &candidate, INTERACTION, Readers::default()).red(0));
    for (from, to) in [
        ("21:39375657", "11:92975544"),
        ("\tC\t", "\tA\t"),
        ("GCST006494", "GCST000001"),
    ] {
        let compared = compare(
            &reference.replace(from, to),
            &candidate.replace(from, to),
            INTERACTION,
            with_catalog(&held),
        );
        assert!(compared.red(0), "{from} {compared:?}");
    }
    let annotated = row("21:39375657", "GWAS_study=GCST000028");
    assert!(compare(&annotated, &candidate, INTERACTION, with_catalog(&held)).red(0));
    let ordinary = row("21:500", INTERACTION_FIELDS);
    for planted in [
        row("21:500", ""),
        ordinary.replace("GCST006494", "GCST000001"),
    ] {
        assert!(compare(&ordinary, &planted, INTERACTION, with_catalog(&held)).red(0));
    }
}

#[test]
fn a_gwas_class_leaves_a_planted_value_on_a_one_allele_row_a_difference() {
    let row = |pair: &str| {
        let extra = if pair.is_empty() {
            "IMPACT=MODERATE;STRAND=1".to_string()
        } else {
            format!("IMPACT=MODERATE;STRAND=1;{pair}")
        };
        default_text(&[default_row("21_100_A/G", "21:100", "G", "ENST1", &extra)])
    };
    let reference = row("GWAS_study=GCST1");
    for klass in [CROSS_PRODUCT, INTERACTION] {
        for candidate in [row(""), row("GWAS_study=GCST2")] {
            let compared = compare(&reference, &candidate, klass, Readers::default());
            assert!(
                compared.red(0) && compared.left == 1,
                "{klass} {compared:?}"
            );
        }
    }
}

#[test]
fn a_file_that_is_no_associations_file_is_refused_by_name() {
    let scratch = Scratch::new("no-catalog");
    let other = scratch.file("other.tsv", "a\tb\n");
    assert_eq!(
        refusal(GwasCatalogPairs::from_catalog(&other)),
        format!("{}: not a GWAS catalog associations file (columns CHR_ID, CHR_POS, STRONGEST SNP-RISK ALLELE, STUDY ACCESSION)", other.display())
    );
    let empty = scratch.file("empty.tsv", "");
    assert!(refusal(GwasCatalogPairs::from_catalog(&empty))
        .ends_with("not a GWAS catalog associations file (columns CHR_ID, CHR_POS, STRONGEST SNP-RISK ALLELE, STUDY ACCESSION)"));
    let missing = scratch.path("nowhere.tsv");
    assert_eq!(
        refusal(GwasCatalogPairs::from_catalog(&missing)),
        format!(
            "cannot read {}: No such file or directory",
            missing.display()
        )
    );
    // The columns are found by name wherever they stand, a line short of one of them is
    // left out, and a carriage return ends a line as a line feed does.
    let reordered = scratch.file(
        "reordered.tsv",
        "STUDY ACCESSION\tDATE\tCHR_POS\tSTRONGEST SNP-RISK ALLELE\tCHR_ID\r\nGCST1\t2020\t100;150\trs1-G; rs2-A\t21;21\r\nGCST2\t2021\t200;250\trs3-T; rs4-C\r\n",
    );
    let held = GwasCatalogPairs::from_catalog(&reordered).unwrap();
    assert_eq!(
        pairs(held.reference_only()),
        BTreeSet::from([["21", "100", "A", "GCST1"], ["21", "150", "G", "GCST1"]])
    );
}

// ---- the SpliceAI score files and their two classes ------------------------------------

const SPLICEAI_HEADER: &str = "##fileformat=VCFv4.2\n##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAIv1.3 synthetic fixture, fabricated values. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
/// A value of the field's shape, not a score of any file.
const SPLICEAI_VALUE: &str = "SpliceAI_pred=GENE1|0.10|0.20|0.30|0.40|-4|-3|2|1";

fn gzip(data: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

/// One BGZF block of `data`: a gzip member whose `BC` extra field holds its size.
fn bgzf_block(data: &[u8]) -> Vec<u8> {
    let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    let body = encoder.finish().unwrap();
    let mut crc = Crc::new();
    crc.update(data);
    let mut block = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00".to_vec();
    block.extend_from_slice(&(body.len() as u16 + 25).to_le_bytes());
    block.extend_from_slice(&body);
    block.extend_from_slice(&crc.sum().to_le_bytes());
    block.extend_from_slice(&(data.len() as u32).to_le_bytes());
    block
}

/// A made-up table of `lines` under `header`, each line a BGZF block of its own, with a
/// tabix index on its second column that files every row of a chromosome in its first bin.
fn write_table(scratch: &Scratch, name: &str, header: &str, lines: &[String]) -> PathBuf {
    let mut bytes = bgzf_block(format!("{header}\n").as_bytes());
    let mut spans: Vec<(String, u64, u64)> = Vec::new();
    for line in lines {
        let at = bytes.len() as u64;
        bytes.extend_from_slice(&bgzf_block(format!("{line}\n").as_bytes()));
        let chromosome = line.split('\t').next().unwrap();
        match spans.iter_mut().find(|(name, _, _)| name == chromosome) {
            Some(span) => span.2 = bytes.len() as u64,
            None => spans.push((chromosome.to_string(), at, bytes.len() as u64)),
        }
    }
    bytes.extend_from_slice(&bgzf_block(b""));
    let names: Vec<u8> = spans
        .iter()
        .flat_map(|(name, _, _)| name.bytes().chain([0]))
        .collect();
    let mut index = b"TBI\x01".to_vec();
    for field in [
        spans.len() as i32,
        0,
        1,
        2,
        2,
        i32::from(b'#'),
        0,
        names.len() as i32,
    ] {
        index.extend_from_slice(&field.to_le_bytes());
    }
    index.extend_from_slice(&names);
    for (_, first, end) in &spans {
        index.extend_from_slice(&1i32.to_le_bytes());
        index.extend_from_slice(&0u32.to_le_bytes());
        index.extend_from_slice(&1i32.to_le_bytes());
        index.extend_from_slice(&(first << 16).to_le_bytes());
        index.extend_from_slice(&(end << 16).to_le_bytes());
        index.extend_from_slice(&0i32.to_le_bytes());
    }
    scratch.file(&format!("{name}.tbi"), gzip(&index));
    scratch.file(name, bytes)
}

fn score_file(scratch: &Scratch, name: &str, lines: &[String]) -> PathBuf {
    write_table(scratch, name, SPLICEAI_HEADER, lines)
}

fn table(files: &[PathBuf], cache: Option<&Path>) -> SpliceAiTable {
    SpliceAiTable::new(files, cache).unwrap()
}

fn with_table(table: &SpliceAiTable) -> Readers<'_> {
    Readers {
        spliceai: Some(table),
        ..Readers::default()
    }
}

/// Made-up scores of three symbols.
fn scores(symbol: &str) -> &'static str {
    match symbol {
        "GENEA" => "0.15|0.08|0.06|0.17|-6|-41|-17|-9",
        "GENEB" => "0.02|0.12|0.10|0.05|-25|-23|37|17",
        _ => "0.31|0.00|0.44|0.09|3|-2|11|-30",
    }
}

/// The entry of `symbol`: the symbol and its scores.
fn own_value(symbol: &str) -> String {
    format!("{symbol}|{}", scores(symbol))
}

/// A score record at `chromosome` `position` of one entry per value: the alternate allele
/// before a symbol and its scores.
fn record_on(
    chromosome: &str,
    position: u64,
    reference: &str,
    alternate: &str,
    values: &[&str],
) -> String {
    let entries: Vec<String> = values
        .iter()
        .map(|value| format!("{alternate}|{value}"))
        .collect();
    format!(
        "{chromosome}\t{position}\t.\t{reference}\t{alternate}\t.\t.\tSpliceAI={}",
        entries.join(",")
    )
}

fn record(position: u64, reference: &str, alternate: &str, values: &[&str]) -> String {
    record_on("21", position, reference, alternate, values)
}

/// What SpliceAI.pm:326-332 make of a record: the first symbol and everything after it.
fn record_value(record: &str) -> &str {
    let info = record.split_once("SpliceAI=").unwrap().1;
    info.split_once('|').unwrap().1
}

/// The records of `symbols` at the deletion 21 5000 GCA G, whose rows stand at
/// 21:5001-5002 with the allele written `-`: as many entries a record as `entries` says.
fn deletion_records(symbols: &[&str], entries: usize) -> Vec<String> {
    let values: Vec<String> = symbols.iter().map(|symbol| own_value(symbol)).collect();
    values
        .chunks(entries)
        .map(|group| {
            let group: Vec<&str> = group.iter().map(String::as_str).collect();
            record(5000, "GCA", "G", &group)
        })
        .collect()
}

/// A made-up cache that names the gene of ENST1, ENST2 and ENST3 on `chromosome`; no score
/// entry names the third's.
fn spliceai_cache(scratch: &Scratch, name: &str, chromosome: &str) -> PathBuf {
    scratch.file(
        &format!("{name}/transcripts/{chromosome}/1-1000000.json"),
        r#"[{"stable_id": "ENST1", "gene_symbol": "GENEA"}, {"stable_id": "ENST2", "gene_symbol": "GENEB"}, {"stable_id": "ENST3", "gene_symbol": "GENEX"}]"#,
    );
    scratch.path(name)
}

/// Default-format rows on ENST1, a GENEA transcript, ENST2, a GENEB transcript, and ENST3,
/// a transcript of a gene no entry names, as many as values are given: each value the
/// row's SpliceAI pairs of the Extra text, or nothing.
#[derive(Clone, Copy)]
struct Rows<'a> {
    symbols: bool,
    name: &'a str,
    allele: &'a str,
    site: &'a str,
}

const DELETION: Rows = Rows {
    symbols: true,
    name: "del_two_genes",
    allele: "-",
    site: "21:5001-5002",
};
/// The rows of the substitution 21 7000 A C: they stand at 21:7000 with the allele C.
const SUBSTITUTION: Rows = Rows {
    symbols: true,
    name: "snv_two_genes",
    allele: "C",
    site: "21:7000",
};

impl Rows<'_> {
    fn text(&self, values: &[&str]) -> String {
        let rows: Vec<String> = [("ENST1", "GENEA"), ("ENST2", "GENEB"), ("ENST3", "GENEX")]
            .iter()
            .zip(values)
            .map(|((feature, gene), value)| {
                let mut extra = "IMPACT=MODIFIER".to_string();
                if self.symbols {
                    extra.push_str(&format!(";SYMBOL={gene}"));
                }
                if !value.is_empty() {
                    extra.push_str(&format!(";{value}"));
                }
                default_row(self.name, self.site, self.allele, feature, &extra)
            })
            .collect();
        default_text(&rows)
    }

    fn unnamed(self) -> Self {
        Rows {
            symbols: false,
            ..self
        }
    }
}

fn pred(value: &str) -> String {
    format!("SpliceAI_pred={value}")
}

/// The nine `SpliceAI_pred_` pairs SpliceAI.pm:276-287 cut from a value, and a flag.
fn split_pairs(value: &str, cutoff: &str) -> String {
    let members = [
        "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL",
    ];
    let mut pairs: Vec<String> = members
        .iter()
        .zip(value.split('|'))
        .map(|(member, part)| format!("SpliceAI_pred_{member}={part}"))
        .collect();
    if !cutoff.is_empty() {
        pairs.push(format!("SpliceAI_cutoff={cutoff}"));
    }
    pairs.join(";")
}

fn chopped_rows(name: &str, rows: &[(&str, &str)]) -> String {
    let rows: Vec<String> = rows
        .iter()
        .map(|(allele, spliceai)| {
            let more = if spliceai.is_empty() {
                String::new()
            } else {
                format!(";{spliceai}")
            };
            default_row(
                name,
                "21:26965500",
                allele,
                "ENST00000284987",
                &format!("IMPACT=MODIFIER;STRAND=1{more}"),
            )
        })
        .collect();
    default_text(&rows)
}

#[test]
fn spliceai_chopped_multiallelic_accepts_the_substitution_scored_beside_a_deletion() {
    let spelt = "21_26965500_C/-/A";
    let reference = chopped_rows(spelt, &[("-", ""), ("A", "")]);
    let candidate = chopped_rows(spelt, &[("-", ""), ("A", SPLICEAI_VALUE)]);
    let compared = compare(&reference, &candidate, CHOPPED, Readers::default());
    assert!(compared.green(1), "{compared:?}");
    assert_eq!(compared.masked["SpliceAI_pred"], 1);
    // Under an identifier the default format's rows may be two records', and the class
    // leaves them.
    let named = "chopped_multiallelic";
    let compared = compare(
        &chopped_rows(named, &[("-", ""), ("A", "")]),
        &chopped_rows(named, &[("-", ""), ("A", SPLICEAI_VALUE)]),
        CHOPPED,
        Readers::default(),
    );
    assert!(compared.red(0), "{compared:?}");
    // The VCF line and the JSON allele string spell the record's alleles.
    let row = |kind: &str, allele: &str, indel: bool, fields: &[(&str, &str)]| {
        entry(
            &format!("rec0|{allele}|ENST1|{kind}"),
            None,
            allele,
            fields,
            Some(indel),
            None,
        )
    };
    let vcf = |score: &str| {
        vec![
            entry(
                "rec0|record",
                None,
                "",
                &[
                    ("CHROM", "21"),
                    ("POS", "26965499"),
                    ("REF", "GC"),
                    ("ALT", "G,GA"),
                ],
                None,
                None,
            ),
            row("Transcript", "-", true, &[("SpliceAI_pred_DS_AG", "")]),
            row("Transcript", "A", false, &[("SpliceAI_pred_DS_AG", score)]),
        ]
    };
    let compared = compare_entries(&vcf(""), &vcf("0.10"), CHOPPED, Readers::default());
    assert!(compared.green(1), "{compared:?}");
    let json = |spliceai: &[(&str, &str)]| {
        vec![
            entry(
                "rec0|record",
                None,
                "",
                &[("id", named), ("allele_string", "C/-/A")],
                None,
                None,
            ),
            row("transcript", "-", true, &[]),
            row("transcript", "A", false, spliceai),
        ]
    };
    let object = [("spliceai", r#"{"DS_AG":0.1,"SYMBOL":"GENE1"}"#)];
    let compared = compare_entries(&json(&[]), &json(&object), CHOPPED, Readers::default());
    assert!(
        compared.green(1) && compared.masked["spliceai"] == 1,
        "{compared:?}"
    );
}

#[test]
fn spliceai_chopped_multiallelic_refuses_a_record_without_a_deletion_or_another_direction() {
    let name = "21_26965500_C/-/A";
    let other_value = SPLICEAI_VALUE.replace("0.10", "0.11");
    for (reference, candidate) in [
        (
            chopped_rows("21_26965500_C/G/A", &[("G", ""), ("A", "")]),
            chopped_rows("21_26965500_C/G/A", &[("G", ""), ("A", SPLICEAI_VALUE)]),
        ),
        (
            chopped_rows("21_26965500_C/A", &[("A", "")]),
            chopped_rows("21_26965500_C/A", &[("A", SPLICEAI_VALUE)]),
        ),
        (
            chopped_rows(name, &[("-", ""), ("A", "")]),
            chopped_rows(name, &[("-", SPLICEAI_VALUE), ("A", "")]),
        ),
        (
            chopped_rows("21_26965500_C/-/CA", &[("-", ""), ("CA", "")]),
            chopped_rows("21_26965500_C/-/CA", &[("-", ""), ("CA", SPLICEAI_VALUE)]),
        ),
        (
            chopped_rows("chopped_multiallelic", &[("-", ""), ("CA", "")]),
            chopped_rows("chopped_multiallelic", &[("-", ""), ("CA", SPLICEAI_VALUE)]),
        ),
        (
            chopped_rows(name, &[("-", ""), ("A", SPLICEAI_VALUE)]),
            chopped_rows(name, &[("-", ""), ("A", "")]),
        ),
        (
            chopped_rows(name, &[("-", ""), ("A", SPLICEAI_VALUE)]),
            chopped_rows(name, &[("-", ""), ("A", &other_value)]),
        ),
    ] {
        let compared = compare(&reference, &candidate, CHOPPED, Readers::default());
        assert!(compared.red(0), "{candidate} {compared:?}");
    }
    // A one-allele record's row, its value missing or wrong.
    let row = |pair: &str| {
        default_text(&[default_row(
            "21_100_A/G",
            "21:100",
            "G",
            "ENST1",
            &format!("IMPACT=MODERATE;STRAND=1{pair}"),
        )])
    };
    let reference = row(";SpliceAI_pred=GENE1|0.10|0.20|0.30|0.40|-4|-3|2|1");
    for candidate in [
        row(""),
        row(";SpliceAI_pred=GENE1|0.11|0.20|0.30|0.40|-4|-3|2|1"),
    ] {
        let compared = compare(&reference, &candidate, CHOPPED, Readers::default());
        assert!(compared.red(0) && compared.left == 1, "{compared:?}");
    }
}

#[test]
fn spliceai_chopped_multiallelic_reads_a_named_records_alleles_from_the_input() {
    // The line's REF and ALT say what a row under an identifier cannot: beside the deletion
    // of GC to G,GCA the CA row at the one-base Location is an insertion's, the A row of
    // GC to G,GA a substitution's.
    for (alternates, allele, rows) in [("G,GA", "A", 1), ("G,GCA", "CA", 0)] {
        let line = format!("21\t26965499\tchopped_multiallelic\tGC\t{alternates}\t.\t.\t.");
        let readers = Readers {
            input: &[&line],
            ..Readers::default()
        };
        let compared = compare(
            &chopped_rows("chopped_multiallelic", &[("-", ""), (allele, "")]),
            &chopped_rows(
                "chopped_multiallelic",
                &[("-", ""), (allele, SPLICEAI_VALUE)],
            ),
            CHOPPED,
            readers,
        );
        assert_eq!(
            (compared.rows, compared.left),
            (rows, 1 - rows),
            "{alternates}"
        );
    }
}

#[test]
fn spliceai_multi_entry_cell_takes_both_rows_of_a_record_of_several_entries() {
    let scratch = Scratch::new("both-rows");
    let two = score_file(
        &scratch,
        "two/indel.synthetic.vcf.gz",
        &deletion_records(&["GENEA", "GENEB"], 2),
    );
    let whole_record = record(
        5000,
        "GCA",
        "G",
        &[&own_value("GENEA"), &own_value("GENEB")],
    );
    let whole = record_value(&whole_record);
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let held = table(std::slice::from_ref(&two), None);
    let reference = DELETION.text(&[&pred(whole), ""]);
    let candidate = DELETION.text(&[&pred(&first), &pred(&second)]);
    assert!(compare(&reference, &candidate, MULTI, with_table(&held)).green(2));
    // Rows that state no SYMBOL take their gene from the cache.
    let cache = spliceai_cache(&scratch, "cache", "21");
    let cached = table(std::slice::from_ref(&two), Some(&cache));
    let unnamed = DELETION.unnamed();
    let compared = compare(
        &unnamed.text(&[&pred(whole), ""]),
        &unnamed.text(&[&pred(&first), &pred(&second)]),
        MULTI,
        with_table(&cached),
    );
    assert!(compared.green(2), "{compared:?}");
    let flagged = compare(
        &DELETION.text(&[&format!("{};SpliceAI_cutoff=FAIL", pred(whole)), ""]),
        &DELETION.text(&[
            &format!("{};SpliceAI_cutoff=FAIL", pred(&first)),
            &format!("{};SpliceAI_cutoff=PASS", pred(&second)),
        ]),
        MULTI,
        with_table(&held),
    );
    assert!(flagged.green(2), "{flagged:?}");
    // The nine split fields: one differs on the first gene's row, nine on the second's.
    let split = compare(
        &DELETION.text(&[&split_pairs(whole, "FAIL"), ""]),
        &DELETION.text(&[&split_pairs(&first, "FAIL"), &split_pairs(&second, "FAIL")]),
        MULTI,
        with_table(&held),
    );
    assert!(split.green(2) && split.masked.len() == 10, "{split:?}");
    assert_eq!(split.masked["SpliceAI_pred_DP_DL"], 2);
    let three_symbols = ["GENEA", "GENEB", "GENEC"];
    let three = table(
        &[score_file(
            &scratch,
            "three/indel.synthetic.vcf.gz",
            &deletion_records(&three_symbols, 3),
        )],
        None,
    );
    let three_values: Vec<String> = three_symbols
        .iter()
        .map(|symbol| own_value(symbol))
        .collect();
    let three_values: Vec<&str> = three_values.iter().map(String::as_str).collect();
    let whole_of_three = record(5000, "GCA", "G", &three_values);
    let compared = compare(
        &DELETION.text(&[&pred(record_value(&whole_of_three)), ""]),
        &candidate,
        MULTI,
        with_table(&three),
    );
    assert!(compared.green(2), "{compared:?}");
    // Two records at the one allele, GENEA's and GENEB's entries then GENEC's alone: each
    // gene one entry names is read.
    let in_pairs = table(
        &[score_file(
            &scratch,
            "pairs/indel.synthetic.vcf.gz",
            &deletion_records(&three_symbols, 2),
        )],
        None,
    );
    assert!(compare(&reference, &candidate, MULTI, with_table(&in_pairs)).green(2));

    // The VCF spelling: the split fields of a CSQ entry, a comma written `&`, on rows that
    // state no SYMBOL field and no Location of their own.
    let vcf_row = |feature: &str, value: &str| {
        let members: Vec<&str> = if value.is_empty() {
            vec![""; 9]
        } else {
            value.split('|').collect()
        };
        let names: Vec<String> = [
            "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL",
        ]
        .iter()
        .map(|member| format!("SpliceAI_pred_{member}"))
        .collect();
        let fields: Vec<(&str, &str)> = names
            .iter()
            .map(String::as_str)
            .zip(members.iter().copied())
            .collect();
        Rc::new(Entry {
            key: format!("rec0|-|{feature}|Transcript"),
            record: Some(0),
            allele: "-".to_string(),
            feature: feature.to_string(),
            fields: fields.into_iter().collect(),
            indel: Some(true),
            allele_location: Some("21:5001-5002".to_string()),
            ..Entry::default()
        })
    };
    let cut: Vec<&str> = whole.split('|').take(9).collect();
    let cut = cut.join("|").replace(',', "&");
    let expected = vec![vcf_row("ENST1", &cut), vcf_row("ENST2", "")];
    let actual = vec![vcf_row("ENST1", &first), vcf_row("ENST2", &second)];
    let compared = compare_entries(&expected, &actual, MULTI, with_table(&cached));
    assert!(
        compared.green(2) && compared.masked.len() == 9,
        "{compared:?}"
    );
    assert!(compare_entries(&expected, &actual, MULTI, with_table(&held)).red(0));

    // The JSON object: its members compared as numbers where both are numbers.
    let json_row = |feature: &str, spliceai: &str| {
        let fields: Vec<(&str, &str)> = if spliceai.is_empty() {
            vec![("variant_allele", "-")]
        } else {
            vec![("variant_allele", "-"), ("spliceai", spliceai)]
        };
        Rc::new(Entry {
            key: format!("rec0|-|{feature}|transcript"),
            record: Some(0),
            allele: "-".to_string(),
            feature: feature.to_string(),
            fields: fields.into_iter().collect(),
            indel: Some(true),
            allele_location: Some("21:5001-5002".to_string()),
            ..Entry::default()
        })
    };
    let expected = vec![
        json_row(
            "ENST1",
            r#"{"DP_AG":-6,"DP_AL":-41,"DP_DG":-17,"DP_DL":"-9,G","DS_AG":0.15,"DS_AL":0.08,"DS_DG":0.06,"DS_DL":0.17,"SYMBOL":"GENEA"}"#,
        ),
        json_row("ENST2", ""),
    ];
    let actual = vec![
        json_row(
            "ENST1",
            r#"{"DP_AG":-6,"DP_AL":-41,"DP_DG":-17,"DP_DL":-9,"DS_AG":0.15,"DS_AL":0.08,"DS_DG":0.06,"DS_DL":0.17,"SYMBOL":"GENEA"}"#,
        ),
        json_row(
            "ENST2",
            r#"{"DP_AG":-25,"DP_AL":-23,"DP_DG":37,"DP_DL":17,"DS_AG":0.02,"DS_AL":0.12,"DS_DG":0.1,"DS_DL":0.05,"SYMBOL":"GENEB"}"#,
        ),
    ];
    let compared = compare_entries(&expected, &actual, MULTI, with_table(&cached));
    assert!(
        compared.green(2) && compared.masked["spliceai"] == 2,
        "{compared:?}"
    );
}

#[test]
fn spliceai_multi_entry_cell_takes_nothing_without_the_score_files() {
    let scratch = Scratch::new("without-files");
    let whole_record = record(
        5000,
        "GCA",
        "G",
        &[&own_value("GENEA"), &own_value("GENEB")],
    );
    let whole = record_value(&whole_record);
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let reference = DELETION.text(&[&pred(whole), ""]);
    let candidate = DELETION.text(&[&pred(&first), &pred(&second)]);
    assert!(compare(&reference, &candidate, MULTI, Readers::default()).red(0));
    let held = table(
        &[score_file(
            &scratch,
            "two/indel.synthetic.vcf.gz",
            &deletion_records(&["GENEA", "GENEB"], 2),
        )],
        None,
    );
    assert!(compare(&reference, &candidate, MULTI, with_table(&held)).green(2));
    // A file that holds no record at the row's variant.
    let compared = compare(
        &reference.replace(DELETION.site, "21:6001-6002"),
        &candidate.replace(DELETION.site, "21:6001-6002"),
        MULTI,
        with_table(&held),
    );
    assert!(compared.red(0), "{compared:?}");
    // Neither a SYMBOL nor a cache names the row's gene.
    let unnamed = DELETION.unnamed();
    let compared = compare(
        &unnamed.text(&[&pred(whole), ""]),
        &unnamed.text(&[&pred(&first), &pred(&second)]),
        MULTI,
        with_table(&held),
    );
    assert!(compared.red(0), "{compared:?}");
}

#[test]
fn spliceai_multi_entry_cell_leaves_a_row_off_the_shape_a_difference() {
    let scratch = Scratch::new("off-shape");
    let whole_record = record(
        5000,
        "GCA",
        "G",
        &[&own_value("GENEA"), &own_value("GENEB")],
    );
    let whole = pred(record_value(&whole_record));
    let (first, second) = (pred(&own_value("GENEA")), pred(&own_value("GENEB")));
    let third = pred(&own_value("GENEC"));
    let joined = table(
        &[score_file(
            &scratch,
            "two/indel.synthetic.vcf.gz",
            &deletion_records(&["GENEA", "GENEB"], 2),
        )],
        None,
    );
    let reference = DELETION.text(&[&whole, ""]);
    let candidate = DELETION.text(&[&first, &second]);
    let single = table(
        &[score_file(
            &scratch,
            "single/indel.synthetic.vcf.gz",
            &deletion_records(&["GENEA", "GENEB"], 1),
        )],
        None,
    );
    assert!(compare(&reference, &candidate, MULTI, with_table(&single)).red(0));
    let off = |reference: &str, candidate: &str, rows: usize| {
        let compared = compare(reference, candidate, MULTI, with_table(&joined));
        assert!(compared.red(rows), "{candidate} {compared:?}");
    };
    // On the first entry's gene's row the reference must write the record's whole value.
    off(&DELETION.text(&["", ""]), &candidate, 1);
    let longer = format!("{whole},G|GENEC|{}", scores("GENEC"));
    off(&DELETION.text(&[&longer, ""]), &candidate, 1);
    // The engine's value must be the entry of the row's own gene.
    off(&reference, &DELETION.text(&[&first, &first]), 1);
    off(
        &reference,
        &DELETION.text(&[&first, &second.replace("0.02", "0.03")]),
        1,
    );
    off(&reference, &DELETION.text(&[&first, &third]), 1);
    off(
        &reference.replace("SYMBOL=GENEB", "SYMBOL=GENEC"),
        &candidate.replace("SYMBOL=GENEB", "SYMBOL=GENEC"),
        1,
    );
    // The cache names each row's gene: the second gene's row writes the third entry, and
    // the row of a gene no entry names writes the second.
    let cache = spliceai_cache(&scratch, "cache", "21");
    let three_symbols = ["GENEA", "GENEB", "GENEC"];
    let three = table(
        &[score_file(
            &scratch,
            "three/indel.synthetic.vcf.gz",
            &deletion_records(&three_symbols, 3),
        )],
        Some(&cache),
    );
    let three_values: Vec<String> = three_symbols
        .iter()
        .map(|symbol| own_value(symbol))
        .collect();
    let three_values: Vec<&str> = three_values.iter().map(String::as_str).collect();
    let whole_three = pred(record_value(&record(5000, "GCA", "G", &three_values)));
    let unnamed = DELETION.unnamed();
    let compared = compare(
        &unnamed.text(&[&whole_three, ""]),
        &unnamed.text(&[&first, &third]),
        MULTI,
        with_table(&three),
    );
    assert!(compared.red(1), "{compared:?}");
    let compared = compare(
        &unnamed.text(&[&whole_three, "", ""]),
        &unnamed.text(&[&first, &second, &second]),
        MULTI,
        with_table(&three),
    );
    assert!(compared.red(2), "{compared:?}");
    // The two flags of the first gene's row agree, and the reference writes none on the
    // second's.
    off(
        &DELETION.text(&[&format!("{whole};SpliceAI_cutoff=FAIL"), ""]),
        &DELETION.text(&[&format!("{first};SpliceAI_cutoff=PASS"), &second]),
        1,
    );
    off(
        &DELETION.text(&[&whole, "SpliceAI_cutoff=FAIL"]),
        &DELETION.text(&[&first, &format!("{second};SpliceAI_cutoff=PASS")]),
        1,
    );
    // Another allele at the site.
    let other_allele = Rows {
        allele: "T",
        ..DELETION
    };
    off(
        &other_allele.text(&[&whole, ""]),
        &other_allele.text(&[&first, &second]),
        0,
    );
    // A record whose INFO column holds a second key.
    let second_key = format!("{};DB", whole_record);
    let keyed = table(
        &[score_file(
            &scratch,
            "second-key/indel.synthetic.vcf.gz",
            &[second_key],
        )],
        None,
    );
    assert!(compare(&reference, &candidate, MULTI, with_table(&keyed)).red(0));
}

#[test]
fn spliceai_multi_entry_cell_takes_the_rows_of_a_single_base_record_of_several_entries() {
    let scratch = Scratch::new("single-base");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let candidate = SUBSTITUTION.text(&[&pred(&first), &pred(&second)]);
    let two = record(7000, "A", "C", &[&first, &second]);
    let reference = SUBSTITUTION.text(&[&pred(record_value(&two)), ""]);
    let held = table(
        &[score_file(
            &scratch,
            "two/snv.synthetic.vcf.gz",
            std::slice::from_ref(&two),
        )],
        None,
    );
    assert!(compare(&reference, &candidate, MULTI, with_table(&held)).green(2));
    // A record of nine entries with GENEA's first and GENEB's fifth.
    let others: Vec<String> = "DEFGHIJ"
        .chars()
        .map(|letter| format!("GENE{letter}|{}", scores("GENEC")))
        .collect();
    let mut values: Vec<&str> = vec![&first];
    values.extend(others[..3].iter().map(String::as_str));
    values.push(&second);
    values.extend(others[3..].iter().map(String::as_str));
    let nine = record(7000, "A", "C", &values);
    assert_eq!(record_value(&nine).matches(',').count(), 8);
    let held = table(
        &[score_file(
            &scratch,
            "nine/snv.synthetic.vcf.gz",
            std::slice::from_ref(&nine),
        )],
        None,
    );
    let whole_of_nine = SUBSTITUTION.text(&[&pred(record_value(&nine)), ""]);
    assert!(compare(&whole_of_nine, &candidate, MULTI, with_table(&held)).green(2));
    // The records read from one base before the row's Location: a single-base record at
    // the position before and, in a second file, an insertion and a deletion anchored at
    // the site and at the position before it, every one of them naming both genes.
    let before = record(6999, "T", "C", &[&first, &second]);
    let indels: Vec<String> = [
        (6999, "T", "TC"),
        (6999, "TA", "T"),
        (7000, "A", "AC"),
        (7000, "AG", "A"),
    ]
    .iter()
    .map(|(position, reference, alternate)| {
        record(*position, reference, alternate, &[&first, &second])
    })
    .collect();
    let beside = table(
        &[
            score_file(&scratch, "beside/snv.synthetic.vcf.gz", &[before, two]),
            score_file(&scratch, "beside/indel.synthetic.vcf.gz", &indels),
        ],
        None,
    );
    assert_eq!(beside.records("21", 6999, 7000).unwrap().len(), 6);
    assert!(compare(&reference, &candidate, MULTI, with_table(&beside)).green(2));
}

#[test]
fn spliceai_multi_entry_cell_leaves_a_single_base_row_a_difference_off_the_shape() {
    let scratch = Scratch::new("single-base-off");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let two = record(7000, "A", "C", &[&first, &second]);
    let reference = SUBSTITUTION.text(&[&pred(record_value(&two)), ""]);
    let candidate = SUBSTITUTION.text(&[&pred(&first), &pred(&second)]);
    assert!(compare(&reference, &candidate, MULTI, Readers::default()).red(0));
    let single = table(
        &[score_file(
            &scratch,
            "single/snv.synthetic.vcf.gz",
            &[
                record(7000, "A", "C", &[&first]),
                record(7000, "A", "C", &[&second]),
            ],
        )],
        None,
    );
    assert!(compare(&reference, &candidate, MULTI, with_table(&single)).red(0));
    let joined = table(
        &[score_file(&scratch, "two/snv.synthetic.vcf.gz", &[two])],
        None,
    );
    let other_entry = SUBSTITUTION.text(&[&pred(&first), &pred(&first)]);
    assert!(compare(&reference, &other_entry, MULTI, with_table(&joined)).red(1));
    let not_held = SUBSTITUTION.text(&[&pred(&format!("{first},C|{}", own_value("GENEC"))), ""]);
    assert!(compare(&not_held, &candidate, MULTI, with_table(&joined)).red(1));
}

/// Made-up entries of GENEA with other scores than `scores` gives it: the first fails a
/// cutoff of 0.5 and the second reaches it.
const GENEA_LOW: &str = "GENEA|0.11|0.21|0.31|0.41|-1|-2|-3|-4";
const GENEA_HIGH: &str = "GENEA|0.12|0.92|0.32|0.42|9|8|7|6";

/// The one-field value of a row whose gene the `entries` of one record name: the first
/// behind its allele and each later one with the allele `C`.
fn entries_value(entries: &[&str]) -> String {
    entries.join(",C|")
}

/// The nine `SpliceAI_pred_` pairs of a row whose gene the `entries` of one record name,
/// each member its values over the entries joined by commas, and a flag.
fn entries_pairs(entries: &[&str], cutoff: &str) -> String {
    let members = [
        "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL",
    ];
    let parts: Vec<Vec<&str>> = entries
        .iter()
        .map(|entry| entry.split('|').collect())
        .collect();
    let mut pairs: Vec<String> = members
        .iter()
        .enumerate()
        .map(|(index, member)| {
            let values: Vec<&str> = parts.iter().map(|entry| entry[index]).collect();
            format!("SpliceAI_pred_{member}={}", values.join(","))
        })
        .collect();
    if !cutoff.is_empty() {
        pairs.push(format!("SpliceAI_cutoff={cutoff}"));
    }
    pairs.join(";")
}

/// The six kinds of record in which several entries name one gene: the engine's row carries
/// every entry of the record that names its gene, in the record's order, and the
/// reference's the record's whole value where the gene is the first entry's and nothing
/// where it is a later entry's, in the one-field and in the nine-field form.
#[test]
fn spliceai_multi_entry_cell_takes_the_rows_of_a_gene_several_entries_of_one_record_name() {
    let scratch = Scratch::new("one-gene-entries");
    let (first, other) = (own_value("GENEA"), own_value("GENEB"));
    let kinds: [(&str, Vec<&str>); 6] = [
        ("equal", vec![&first, &first]),
        ("differ", vec![GENEA_LOW, GENEA_HIGH]),
        ("four", vec![GENEA_LOW, &first, GENEA_HIGH, &first]),
        ("twice-then-other", vec![&first, GENEA_LOW, &other]),
        ("around-other", vec![GENEA_LOW, &other, &first]),
        ("other-then-twice", vec![&other, GENEA_LOW, GENEA_HIGH]),
    ];
    for (kind, entries) in kinds {
        let line = record(7000, "A", "C", &entries);
        let held = table(
            &[score_file(
                &scratch,
                &format!("{kind}/snv.synthetic.vcf.gz"),
                std::slice::from_ref(&line),
            )],
            None,
        );
        let of = |gene: &str| -> Vec<&str> {
            entries
                .iter()
                .copied()
                .filter(|entry| entry.starts_with(gene))
                .collect()
        };
        let (own, others) = (of("GENEA"), of("GENEB"));
        let leads = entries[0].starts_with("GENEA");
        let whole = record_value(&line);

        // The one field: a record that names GENEA alone gives both engines one value.
        let candidate_others = if others.is_empty() {
            String::new()
        } else {
            pred(&entries_value(&others))
        };
        let candidate = SUBSTITUTION.text(&[&pred(&entries_value(&own)), &candidate_others]);
        let reference = if leads {
            SUBSTITUTION.text(&[&pred(whole), ""])
        } else {
            SUBSTITUTION.text(&["", &pred(whole)])
        };
        let differing = if others.is_empty() { 0 } else { 2 };
        assert!(
            compare(&reference, &candidate, MULTI, with_table(&held)).green(differing),
            "{kind}: one field"
        );

        // The nine fields, with the flags a cutoff of 0.5 gives each engine.
        let largest = |entry: &str| -> f64 {
            entry.split('|').collect::<Vec<_>>()[1..5]
                .iter()
                .map(|score| score.parse::<f64>().unwrap())
                .fold(0.0, f64::max)
        };
        let verdict_of = |reaches: bool| if reaches { "PASS" } else { "FAIL" };
        let ours_flag = verdict_of(own.iter().any(|entry| largest(entry) >= 0.5));
        let theirs_flag = verdict_of(largest(entries[0]) >= 0.5);
        let other_flag = verdict_of(others.iter().any(|entry| largest(entry) >= 0.5));
        let candidate_others = if others.is_empty() {
            String::new()
        } else {
            entries_pairs(&others, other_flag)
        };
        let candidate = SUBSTITUTION.text(&[&entries_pairs(&own, ours_flag), &candidate_others]);
        let cut = split_pairs(whole, theirs_flag);
        let reference = if leads {
            SUBSTITUTION.text(&[&cut, ""])
        } else {
            SUBSTITUTION.text(&["", &cut])
        };
        let compared = compare(&reference, &candidate, MULTI, with_table(&held));
        assert!(
            compared.left == 0 && compared.rows >= 1,
            "{kind}: nine fields"
        );
    }
}

/// A row of a gene several entries of one record name stays a difference off the shape:
/// an engine that writes the later entry alone, the first alone or the entries in another
/// order, and a pair of flags no one cutoff gives.
#[test]
fn spliceai_multi_entry_cell_leaves_a_row_of_a_gene_several_entries_name_a_difference_off_the_shape(
) {
    let scratch = Scratch::new("one-gene-entries-refused");
    let other = own_value("GENEB");
    let line = record(7000, "A", "C", &[GENEA_LOW, &other, GENEA_HIGH]);
    let held = table(
        &[score_file(
            &scratch,
            "around/snv.synthetic.vcf.gz",
            std::slice::from_ref(&line),
        )],
        None,
    );
    let whole = record_value(&line);
    let reference = SUBSTITUTION.text(&[&pred(whole), ""]);
    let others = pred(&other);
    for (what, value) in [
        ("the later entry alone", GENEA_HIGH.to_string()),
        ("the first entry alone", GENEA_LOW.to_string()),
        (
            "the entries in the other order",
            entries_value(&[GENEA_HIGH, GENEA_LOW]),
        ),
    ] {
        let candidate = SUBSTITUTION.text(&[&pred(&value), &others]);
        assert!(
            compare(&reference, &candidate, MULTI, with_table(&held)).red(1),
            "{what}"
        );
    }

    // Flags: the first entry fails 0.5 and a later entry of the gene reaches it, so FAIL
    // beside PASS stands and PASS beside FAIL does not; over entries none of which passes
    // the first's largest score, no cutoff parts the two flags.
    let own = [GENEA_LOW, GENEA_HIGH];
    let cut = |flag: &str| SUBSTITUTION.text(&[&split_pairs(whole, flag), ""]);
    let joined = |flag: &str| {
        SUBSTITUTION.text(&[
            &entries_pairs(&own, flag),
            &entries_pairs(&[&other], "FAIL"),
        ])
    };
    assert_eq!(
        compare(&cut("FAIL"), &joined("PASS"), MULTI, with_table(&held)).left,
        0
    );
    assert!(compare(&cut("PASS"), &joined("FAIL"), MULTI, with_table(&held)).left > 0);
    let level = record(7000, "A", "C", &[GENEA_HIGH, GENEA_LOW]);
    let level_held = table(
        &[score_file(
            &scratch,
            "level/snv.synthetic.vcf.gz",
            std::slice::from_ref(&level),
        )],
        None,
    );
    let level_whole = record_value(&level);
    let level_own = [GENEA_HIGH, GENEA_LOW];
    let level_cut = |flag: &str| SUBSTITUTION.text(&[&split_pairs(level_whole, flag)]);
    let level_joined = |flag: &str| SUBSTITUTION.text(&[&entries_pairs(&level_own, flag)]);
    assert_eq!(
        compare(
            &level_cut("PASS"),
            &level_joined("PASS"),
            MULTI,
            with_table(&level_held)
        )
        .left,
        0
    );
    assert!(
        compare(
            &level_cut("FAIL"),
            &level_joined("PASS"),
            MULTI,
            with_table(&level_held)
        )
        .left
            > 0
    );
}

#[test]
fn spliceai_multi_entry_cell_reads_every_record_line_at_the_allele() {
    let scratch = Scratch::new("every-line");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let changed = format!("GENEA|{}", scores("GENEC"));
    let other = format!("GENEB|{}", scores("GENEA"));
    let third = format!("GENEA|{}", scores("GENEB"));
    let last = format!("GENEB|{}", scores("GENEC"));
    let both = [first.as_str(), second.as_str()];
    let reversed = [second.as_str(), first.as_str()];
    let own = [pred(&first), pred(&second)];
    let whole = |values: &[&str]| pred(record_value(&record(7000, "A", "C", values)));
    type Shape<'a> = (
        &'a str,
        Vec<&'a [&'a str]>,
        Vec<String>,
        Vec<String>,
        bool,
        usize,
    );
    let first_twice = [first.as_str(), first.as_str()];
    let changed_other = [changed.as_str(), other.as_str()];
    let third_last = [third.as_str(), last.as_str()];
    let other_changed = [other.as_str(), changed.as_str()];
    let first_changed = [first.as_str(), changed.as_str()];
    let shapes: Vec<Shape> = vec![
        (
            "one-line",
            vec![&both],
            vec![whole(&both), String::new()],
            own.to_vec(),
            true,
            2,
        ),
        (
            "twice",
            vec![&both, &both],
            vec![whole(&both), String::new()],
            own.to_vec(),
            true,
            2,
        ),
        (
            "three-times",
            vec![&both, &both, &both],
            vec![whole(&both), String::new()],
            own.to_vec(),
            true,
            2,
        ),
        (
            "other-order",
            vec![&both, &reversed],
            vec![whole(&both), whole(&reversed)],
            own.to_vec(),
            true,
            2,
        ),
        (
            "one-entry-last",
            vec![&both, &both[..1]],
            vec![whole(&both[..1]), String::new()],
            own.to_vec(),
            true,
            1,
        ),
        (
            "one-entry-first",
            vec![&both[..1], &both],
            vec![whole(&both), String::new()],
            own.to_vec(),
            true,
            2,
        ),
        (
            "one-gene-twice",
            vec![&first_twice],
            vec![whole(&[&first, &first]), String::new()],
            vec![own[0].clone(), String::new()],
            false,
            0,
        ),
        (
            "gene-unnamed",
            vec![&both, &both],
            vec![whole(&both), String::new(), String::new()],
            vec![own[0].clone(), own[1].clone(), own[1].clone()],
            false,
            2,
        ),
        (
            "other-scores",
            vec![&both, &changed_other],
            vec![whole(&[&changed, &other]), String::new()],
            vec![pred(&changed), pred(&other)],
            true,
            2,
        ),
        (
            "other-scores-three-times",
            vec![&both, &changed_other, &third_last],
            vec![whole(&[&third, &last]), String::new()],
            vec![pred(&third), pred(&last)],
            true,
            2,
        ),
        (
            "other-scores-other-order",
            vec![&both, &other_changed],
            vec![whole(&both), whole(&[&other, &changed])],
            vec![pred(&changed), pred(&other)],
            true,
            2,
        ),
        (
            "earlier-line",
            vec![&both, &changed_other],
            vec![whole(&[&changed, &other]), String::new()],
            vec![own[0].clone(), pred(&other)],
            false,
            1,
        ),
        (
            "later-entry-of-a-line",
            vec![&first_changed],
            vec![whole(&[&first, &changed]), String::new()],
            vec![pred(&changed), String::new()],
            false,
            0,
        ),
        (
            "earlier-entry-of-a-line",
            vec![&first_changed],
            vec![whole(&[&first, &changed]), String::new()],
            vec![own[0].clone(), String::new()],
            false,
            0,
        ),
    ];
    for (name, lines, reference, candidate, green, rows) in shapes {
        let records: Vec<String> = lines
            .iter()
            .map(|line| record(7000, "A", "C", line))
            .collect();
        let held = table(
            &[score_file(
                &scratch,
                &format!("{name}/snv.synthetic.vcf.gz"),
                &records,
            )],
            None,
        );
        let texts = |values: &[String]| {
            let values: Vec<&str> = values.iter().map(String::as_str).collect();
            SUBSTITUTION.text(&values)
        };
        let (reference, candidate) = (texts(&reference), texts(&candidate));
        let compared = compare(&reference, &candidate, MULTI, with_table(&held));
        assert_eq!((compared.left == 0, compared.rows), (green, rows), "{name}");
        assert!(
            compare(&reference, &candidate, MULTI, Readers::default()).red(0),
            "{name}"
        );
    }
    // The flags: equal on the row of the gene that leads a line, the engine's alone on the
    // other row.
    let line = record(7000, "A", "C", &both);
    let twice = table(
        &[score_file(
            &scratch,
            "flagged/snv.synthetic.vcf.gz",
            &[line.clone(), line],
        )],
        None,
    );
    let flagged = SUBSTITUTION.text(&[&format!("{};SpliceAI_cutoff=FAIL", whole(&both)), ""]);
    let flags = |first_flag: &str| {
        SUBSTITUTION.text(&[
            &format!("{};SpliceAI_cutoff={first_flag}", own[0]),
            &format!("{};SpliceAI_cutoff=PASS", own[1]),
        ])
    };
    assert!(compare(&flagged, &flags("FAIL"), MULTI, with_table(&twice)).green(2));
    assert!(compare(&flagged, &flags("PASS"), MULTI, with_table(&twice)).red(1));
}

/// The rows each ground of `GROUNDS` refuses among the made-up rows of `refused_cases`.
const REFUSED: [(&str, u64); 12] = [
    ("no_candidate_value", 1),
    ("no_span", 1),
    ("other_info", 1),
    ("several_references", 1),
    ("no_record", 2),
    ("no_gene", 1),
    ("gene_not_named", 1),
    ("lone_entry", 1),
    ("candidate_value", 1),
    ("reference_not_silent", 1),
    ("reference_value", 1),
    ("reference_flag", 1),
];

/// A made-up row and the record lines at its variant: its Location, its gene (none for a
/// row that states none), the lines, and the SpliceAI pairs the reference and the engine
/// write.
type RefusedCase = (String, &'static str, Vec<String>, String, String);

/// For each ground of the class one row that meets it first, two for the ground of no
/// record line, then the two rows of a line of GENEA's entry and GENEB's, which the class
/// takes; on `chromosome`, the positions counted from `base`. The row under the ground of
/// the engine's value writes the entry of an earlier line than the last that names its
/// gene.
fn refused_cases(chromosome: &str, base: u64) -> Vec<RefusedCase> {
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let changed = format!("GENEA|{}", scores("GENEC"));
    let whole = pred(&format!("{first},C|{second}"));
    let own = pred(&first);
    let at = |offset: u64| format!("{chromosome}:{}", base + offset);
    let line =
        |offset: u64, entries: &[&str]| record_on(chromosome, base + offset, "A", "C", entries);
    let pair = [first.as_str(), second.as_str()];
    vec![
        (
            at(410),
            "GENEA",
            vec![line(410, &pair)],
            whole.clone(),
            String::new(),
        ),
        (
            chromosome.to_string(),
            "GENEA",
            vec![],
            whole.clone(),
            own.clone(),
        ),
        (
            at(430),
            "GENEA",
            vec![format!("{};DB", line(430, &pair))],
            whole.clone(),
            own.clone(),
        ),
        (
            at(440),
            "GENEA",
            vec![
                line(440, &pair),
                record_on(chromosome, base + 440, "G", "C", &pair),
            ],
            whole.clone(),
            own.clone(),
        ),
        (at(450), "GENEA", vec![], whole.clone(), own.clone()),
        (at(455), "GENEA", vec![], whole.clone(), own.clone()),
        (
            at(460),
            "",
            vec![line(460, &pair)],
            whole.clone(),
            own.clone(),
        ),
        (
            at(470),
            "GENEX",
            vec![line(470, &pair)],
            String::new(),
            pred(&second),
        ),
        (
            at(490),
            "GENEA",
            vec![line(490, &pair[..1])],
            own.clone(),
            pred(&changed),
        ),
        (
            at(500),
            "GENEA",
            vec![line(500, &pair), line(500, &[&changed, &second])],
            pred(&format!("{changed},C|{second}")),
            own.clone(),
        ),
        (
            at(510),
            "GENEA",
            vec![line(510, &[&second, &first])],
            format!("{};SpliceAI_cutoff=FAIL", pred(&second)),
            format!("{own};SpliceAI_cutoff=PASS"),
        ),
        (
            at(520),
            "GENEA",
            vec![line(520, &pair)],
            pred(&second),
            own.clone(),
        ),
        (
            at(530),
            "GENEA",
            vec![line(530, &pair)],
            format!("{whole};SpliceAI_cutoff=FAIL"),
            format!("{own};SpliceAI_cutoff=PASS"),
        ),
        (at(550), "GENEA", vec![line(550, &pair)], whole, own),
        (at(550), "GENEB", vec![], String::new(), pred(&second)),
    ]
}

fn refused_lines(cases: &[RefusedCase]) -> Vec<String> {
    cases.iter().flat_map(|case| case.2.clone()).collect()
}

fn expected_counts() -> RefusalCounts {
    RefusalCounts {
        asked: 15,
        refused: 13,
        grounds: REFUSED
            .iter()
            .map(|(ground, rows)| (ground.to_string(), *rows))
            .collect(),
    }
}

#[test]
fn spliceai_multi_entry_cell_counts_the_rows_it_refuses_by_ground() {
    let scratch = Scratch::new("by-ground");
    let cases = refused_cases("21", 7000);
    let rows = |side: usize, chosen: &[RefusedCase]| {
        let rows: Vec<String> = chosen
            .iter()
            .map(|(location, gene, _, reference, engine)| {
                let mut extra = "IMPACT=MODIFIER".to_string();
                if !gene.is_empty() {
                    extra.push_str(&format!(";SYMBOL={gene}"));
                }
                let written = [reference, engine][side];
                if !written.is_empty() {
                    extra.push_str(&format!(";{written}"));
                }
                default_row(
                    &format!("snv_{location}_{gene}"),
                    location,
                    "C",
                    &format!("ENST_{gene}"),
                    &extra,
                )
            })
            .collect();
        default_text(&rows)
    };
    let held = table(
        &[score_file(
            &scratch,
            "scores/snv.synthetic.vcf.gz",
            &refused_lines(&cases),
        )],
        None,
    );
    let grounds: Vec<&str> = GROUNDS.iter().map(|(ground, _)| *ground).collect();
    assert_eq!(grounds, REFUSED.map(|(ground, _)| ground));
    compare(&rows(0, &cases), &rows(1, &cases), MULTI, with_table(&held));
    let counts = held.refusal_counts();
    assert_eq!(counts, expected_counts());
    let lines = refusal_lines(&counts, "");
    assert_eq!(
        lines[0],
        "spliceai_multi_entry_cell under --spliceai-table: rows asked about 15, rows refused 13"
    );
    assert_eq!(
        lines[5],
        "spliceai_multi_entry_cell: rows refused 2: no record line of the score files stands at the row's variant"
    );
    assert_eq!(lines.len(), 13);
    for word in ["GENE", "ENST", "74", "75", "0.1", "snv_"] {
        assert!(!lines.iter().any(|line| line.contains(word)), "{word}");
    }
    // Each comparison counts its own rows: the two rows the class takes, then none.
    held.clear_refusals();
    let taken = &cases[cases.len() - 2..];
    assert!(compare(&rows(0, taken), &rows(1, taken), MULTI, with_table(&held)).green(2));
    let counts = held.refusal_counts();
    assert_eq!(
        (counts.asked, counts.refused, counts.grounds.len()),
        (2, 0, 0)
    );
    assert_eq!(
        refusal_lines(&counts, ", engine vep-rs"),
        ["spliceai_multi_entry_cell under --spliceai-table, engine vep-rs: rows asked about 2, rows refused 0"]
    );
    held.clear_refusals();
    compare(&rows(0, &cases), &rows(0, &cases), MULTI, with_table(&held));
    assert_eq!(held.refusal_counts(), RefusalCounts::default());
}

#[test]
fn spliceai_table_reads_the_records_at_a_span_and_refuses_what_is_no_score_file() {
    let scratch = Scratch::new("records");
    let two = score_file(
        &scratch,
        "two/indel.synthetic.vcf.gz",
        &deletion_records(&["GENEA", "GENEB"], 2),
    );
    let held = table(std::slice::from_ref(&two), None);
    let info = format!(
        "SpliceAI=G|GENEA|{},G|GENEB|{}",
        scores("GENEA"),
        scores("GENEB")
    );
    let expected = vec![ScoreRecord {
        position: 5000,
        reference: "GCA".to_string(),
        alternate: "G".to_string(),
        info,
    }];
    assert_eq!(held.records("21", 4999, 5001).unwrap(), expected);
    assert_eq!(held.records("chr21", 5000, 5000).unwrap(), expected);
    assert!(held.records("21", 5001, 5002).unwrap().is_empty());
    assert!(held.records("22", 5000, 5000).unwrap().is_empty());
    let single = score_file(
        &scratch,
        "single/indel.synthetic.vcf.gz",
        &deletion_records(&["GENEA", "GENEB"], 1),
    );
    let both = table(&[two.clone(), single], None);
    let commas: Vec<usize> = both
        .records("21", 5000, 5000)
        .unwrap()
        .iter()
        .map(|record| record.info.matches(',').count())
        .collect();
    assert_eq!(commas, [1, 0, 0]);
    // A line that lacks one of the four cells, or whose POS is no whole number, is left
    // out.
    let partial = table(
        &[score_file(
            &scratch,
            "partial/snv.synthetic.vcf.gz",
            &[
                "21\t7000\t.\tA\tC\t.\t.\t".to_string(),
                "21\t7000\t.\tA\t\t.\t.\tSpliceAI=C|GENEA|1".to_string(),
                "21\t7000x\t.\tA\tC\t.\t.\tSpliceAI=C|GENEA|1".to_string(),
                "21\t7000\t.\tA\tC\t.\t.\tSpliceAI=C|GENEA|1".to_string(),
            ],
        )],
        None,
    );
    assert_eq!(partial.records("21", 7000, 7000).unwrap().len(), 1);
    let seven = write_table(
        &scratch,
        "seven.tsv.gz",
        "#chr\thg19_pos\tref\talt\taaref\taaalt\tREVEL",
        &["21\t5000\tA\tT\tV\tL\t0.300".to_string()],
    );
    assert_eq!(
        refusal(SpliceAiTable::new(std::slice::from_ref(&seven), None)),
        format!(
            "--spliceai-table: {} is not a VCF (a header line of eight columns)",
            seven.display()
        )
    );
    let no_cache = scratch.path("no-cache");
    assert_eq!(
        refusal(SpliceAiTable::new(
            std::slice::from_ref(&two),
            Some(&no_cache)
        )),
        format!(
            "--json-cache: {} has no transcripts/ directory",
            no_cache.display()
        )
    );
    fs::remove_file(format!("{}.tbi", two.display())).unwrap();
    assert!(
        refusal(SpliceAiTable::new(std::slice::from_ref(&two), None)).starts_with(&format!(
            "--spliceai-table: cannot read the tabix index of {two}: none of {two}.csi, ",
            two = two.display()
        ))
    );
}

#[test]
fn a_rows_gene_is_read_from_the_cache_by_its_transcript() {
    let scratch = Scratch::new("cache-genes");
    let two = score_file(
        &scratch,
        "two/indel.synthetic.vcf.gz",
        &deletion_records(&["GENEA", "GENEB"], 2),
    );
    let whole_record = record(
        5000,
        "GCA",
        "G",
        &[&own_value("GENEA"), &own_value("GENEB")],
    );
    let whole = pred(record_value(&whole_record));
    let (first, second) = (pred(&own_value("GENEA")), pred(&own_value("GENEB")));
    let unnamed = DELETION.unnamed();
    let reference = unnamed.text(&[&whole, ""]);
    let candidate = unnamed.text(&[&first, &second]);
    let cached = |name: &str, chromosome: &str, shard: &str| {
        scratch.file(
            &format!("{name}/transcripts/{chromosome}/1-1000000.json"),
            shard,
        );
        table(std::slice::from_ref(&two), Some(&scratch.path(name)))
    };
    // A Feature carrying a version the cache's stable id lacks is looked up without it,
    // and a later transcript of a stable id replaces an earlier one's symbol.
    let held = cached(
        "versions",
        "21",
        r#"[{"stable_id": "ENST1", "gene_symbol": "GENEX"}, null, {"stable_id": "ENST1", "gene_symbol": "GENEA"}, {"stable_id": "ENST2", "gene_symbol": "GENEB"}]"#,
    );
    let versioned = |text: &str| text.replace("\tENST1\t", "\tENST1.4\t");
    let compared = compare(
        &versioned(&reference),
        &versioned(&candidate),
        MULTI,
        with_table(&held),
    );
    assert!(compared.green(2), "{compared:?}");
    // A row whose Location names no directory of the cache is looked up in every shard,
    // one that is a range-keyed map included.
    let held = cached(
        "elsewhere",
        "chr21",
        r#"{"1-1000000": [{"stable_id": "ENST1", "gene_symbol": "GENEA"}, {"stable_id": "ENST2", "gene_symbol": "GENEB"}]}"#,
    );
    assert!(compare(&reference, &candidate, MULTI, with_table(&held)).green(2));
    // A transcript without a symbol, or with an empty one, names no gene; a symbol that is
    // no text names a gene no entry names.
    let held = cached(
        "unfilled",
        "21",
        r#"[{"stable_id": "ENST1"}, {"stable_id": "ENST2", "gene_symbol": ""}, {"stable_id": "ENST3", "gene_symbol": 7}]"#,
    );
    let compared = compare(
        &unnamed.text(&[&whole, "", ""]),
        &unnamed.text(&[&first, &second, &second]),
        MULTI,
        with_table(&held),
    );
    assert!(compared.red(0), "{compared:?}");
    assert_eq!(
        held.refusal_counts().grounds,
        [
            ("no_gene".to_string(), 2),
            ("gene_not_named".to_string(), 1)
        ]
    );
    // A shard that is no JSON is refused by name.
    let held = cached("damaged", "21", "[{");
    let rows = parse_default(&reference).unwrap();
    let engine = parse_default(&candidate).unwrap();
    let found = BuiltIn.accepts(&Consult {
        klass: MULTI,
        e: &rows[0],
        name: "SpliceAI_pred",
        expected: "x",
        actual: "y",
        a: Some(&engine[0]),
        record_rows: None,
        candidate_rows: None,
        declared: &BTreeSet::new(),
        assembly: None,
        source: None,
        gwas_pairs: None,
        revel_table: None,
        utrannotator_sweep: None,
        downstream_sequences: None,
        splice_windows: None,
        spliceai_table: Some(&held),
        exon_spans: None,
    });
    let shard = scratch.path("damaged/transcripts/21/1-1000000.json");
    assert!(refusal(found).starts_with(&format!("{} is not JSON", shard.display())));
}

// ---- the score command -----------------------------------------------------------------

const ERROR: &str = "ERROR: [score_plugin_tuples] ";
const MULTI_TABLE: &str =
    "spliceai_multi_entry_cell\tSpliceAI\tSpliceAI_pred,SpliceAI_pred_*,SpliceAI_cutoff\tmet\t-\n";
const SPLICEAI_NOTICE: &str = "no --spliceai-table: spliceai_multi_entry_cell sets nothing aside";
const SPLICEAI_CACHE_NOTICE: &str =
    "no --json-cache: spliceai_multi_entry_cell sets aside no row that writes no SYMBOL";
/// The tree's class table.
const COMMITTED_CLASS_TABLE: &str =
    include_str!("../../../scripts/concordance/plugin_class_table.tsv");

/// The exit status, the lines of standard output and standard error of one command line.
fn run_main(arguments: &[&str]) -> (i32, Vec<String>, String) {
    let arguments: Vec<String> = arguments.iter().map(|word| word.to_string()).collect();
    let (mut out, mut error) = (Vec::new(), Vec::new());
    let code = run(
        &arguments,
        &BuiltInRegistry,
        &vep_plugin_pairs::score::api::BuiltIn,
        &mut out,
        &mut error,
    );
    let out = String::from_utf8(out).unwrap();
    (
        code,
        out.lines().map(str::to_string).collect(),
        String::from_utf8(error).unwrap(),
    )
}

/// One scored run of two outputs under the class table `classes`: the files stand under
/// `name` in the scratch, `more` follows the options every run takes.
struct Scored {
    code: i32,
    lines: Vec<String>,
    error: String,
    report: Value,
}

fn scored(
    scratch: &Scratch,
    name: &str,
    plugin: &str,
    classes: &str,
    reference: &[String],
    engine: &[String],
    more: &[&str],
) -> Scored {
    let output = |side: &str, rows: &[String]| {
        scratch.file(
            &format!("{name}/{side}.txt"),
            format!("{HEADER}\n{}", rows.concat()),
        )
    };
    let reference = output("reference", reference);
    let engine = format!("engine={}", output("engine", engine).display());
    let table = scratch.file(&format!("{name}/classes.tsv"), classes);
    let work = scratch.path(&format!("{name}/work"));
    let report = scratch.path(&format!("{name}/report.json"));
    let mut line = vec![
        "--reference",
        text_of(&reference),
        "--engine",
        &engine,
        "--plugins",
        plugin,
        "--class-table",
        text_of(&table),
        "--work-dir",
        text_of(&work),
        "--report",
        text_of(&report),
    ];
    line.extend_from_slice(more);
    let (code, lines, error) = run_main(&line);
    let report =
        fs::read(&report).map_or(Value::Null, |bytes| serde_json::from_slice(&bytes).unwrap());
    Scored {
        code,
        lines,
        error,
        report,
    }
}

impl Scored {
    /// The report's entry of `source` and the engine in one projection.
    fn entry(&self, source: &str, projection: &str) -> &Value {
        self.report["scores"]
            .as_array()
            .unwrap()
            .iter()
            .find(|entry| entry["source"] == source && entry["projection"] == projection)
            .unwrap_or_else(|| panic!("no entry of {source}, {projection}: {}", self.error))
    }

    /// The keys under which `klass` set a tuple aside and the rows left, in the feature
    /// projection.
    fn aside(&self, source: &str, klass: &str) -> (i64, i64) {
        let entry = self.entry(source, "feature");
        (
            entry["masked"][klass].as_i64().unwrap(),
            entry["residual_rows"].as_i64().unwrap(),
        )
    }

    fn passed(&self) -> &Scored {
        assert_eq!((self.code, self.error.as_str()), (0, ""));
        self
    }
}

fn counts(entry: &Value) -> [i64; 6] {
    [
        "perl",
        "rust",
        "intersection",
        "adj_perl",
        "adj_rust",
        "adj_intersection",
    ]
    .map(|key| entry[key].as_i64().unwrap_or(-1))
}

fn row(location: &str, allele: &str, feature: &str, extra: &str, uploaded: &str) -> String {
    format!("{uploaded}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\t{extra}\n")
}

/// The scores of two symbols at the deletion 3 7000199 GCA G, whose rows stand at
/// 3:7000200-7000201 with the allele written `-`, or at the substitution 3 7000300 A C,
/// whose rows stand at 3:7000300 with the allele C.
#[derive(Clone, Copy)]
struct Variant {
    position: u64,
    reference: &'static str,
    alternate: &'static str,
    file: &'static str,
}

const DELETED: Variant = Variant {
    position: 7000199,
    reference: "GCA",
    alternate: "G",
    file: "indel.synthetic.vcf.gz",
};
const SUBSTITUTED: Variant = Variant {
    position: 7000300,
    reference: "A",
    alternate: "C",
    file: "snv.synthetic.vcf.gz",
};

impl Variant {
    /// A score file of the variant's record lines, each the entries of one line.
    fn score_file(&self, scratch: &Scratch, directory: &str, lines: &[&[&str]]) -> PathBuf {
        let records: Vec<String> = lines
            .iter()
            .map(|entries| record_on("3", self.position, self.reference, self.alternate, entries))
            .collect();
        score_file(scratch, &format!("{directory}/{}", self.file), &records)
    }

    /// A row of the variant on `feature` with `value` as its SpliceAI_pred and `flag` as
    /// its cutoff flag, or without either.
    fn row(&self, feature: &str, value: &str, flag: &str) -> String {
        let mut extra = "IMPACT=MODIFIER".to_string();
        if !value.is_empty() {
            extra.push_str(&format!(";SpliceAI_pred={value}"));
        }
        if !flag.is_empty() {
            extra.push_str(&format!(";SpliceAI_cutoff={flag}"));
        }
        if self.file == "snv.synthetic.vcf.gz" {
            row("3:7000300", "C", feature, &extra, "snv_two_genes")
        } else {
            row("3:7000200-7000201", "-", feature, &extra, "del_two_genes")
        }
    }

    /// The whole value of a line of `values`, behind its first allele.
    fn whole(&self, values: &[&str]) -> String {
        values.join(&format!(",{}|", self.alternate))
    }

    fn rows(&self, values: &[&str]) -> Vec<String> {
        ["ENST1", "ENST2", "ENST3"]
            .iter()
            .zip(values)
            .map(|(feature, value)| self.row(feature, value, ""))
            .collect()
    }
}

fn table_words<'a>(files: &[&'a Path]) -> Vec<&'a str> {
    files
        .iter()
        .flat_map(|file| ["--spliceai-table", text_of(file)])
        .collect()
}

#[test]
fn spliceai_table_is_handed_to_the_classes_of_spliceai_which_set_nothing_aside_without_it() {
    let scratch = Scratch::new("handed");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let whole = DELETED.whole(&[&first, &second]);
    let reference = [
        DELETED.row("ENST1", &whole, "FAIL"),
        DELETED.row("ENST2", "", ""),
    ];
    let engine = [
        DELETED.row("ENST1", &first, "FAIL"),
        DELETED.row("ENST2", &second, "PASS"),
    ];
    let run = |name: &str, more: &[&str]| {
        scored(
            &scratch,
            name,
            "SpliceAI",
            MULTI_TABLE,
            &reference,
            &engine,
            more,
        )
    };
    let without = run("without", &[]);
    assert!(without.passed().report["spliceai_table"].is_null());
    assert_eq!(
        counts(without.entry("SpliceAI", "feature")),
        [2, 4, 1, 2, 4, 1]
    );
    assert_eq!(without.aside("SpliceAI", MULTI), (0, 4));
    let joined = DELETED.score_file(&scratch, "joined", &[&[&first, &second]]);
    let cache = spliceai_cache(&scratch, "cache", "3");
    let unnamed = run("unnamed", &["--spliceai-table", text_of(&joined)]);
    assert_eq!(unnamed.passed().aside("SpliceAI", MULTI), (0, 4));
    let both = [
        "--spliceai-table",
        text_of(&joined),
        "--json-cache",
        text_of(&cache),
    ];
    let report = run("with", &both);
    assert_eq!(
        report.passed().report["spliceai_table"],
        json!([{"path": text_of(&joined)}])
    );
    assert_eq!(
        counts(report.entry("SpliceAI", "feature")),
        [2, 4, 1, 1, 1, 1]
    );
    assert_eq!(report.aside("SpliceAI", MULTI), (3, 0));
    let allele = report.entry("SpliceAI", "allele");
    assert_eq!(counts(allele)[3..], [1, 1, 1]);
    assert_eq!(allele["residual_rows"], 0);
    // A file that scores each gene in a record of its own.
    let single = DELETED.score_file(&scratch, "single", &[&[&first], &[&second]]);
    let apart = run(
        "apart",
        &[
            "--spliceai-table",
            text_of(&single),
            "--json-cache",
            text_of(&cache),
        ],
    );
    assert_eq!(apart.passed().aside("SpliceAI", MULTI), (0, 4));
}

#[test]
fn spliceai_rows_off_the_shape_of_a_record_of_several_entries_stay_residual() {
    let scratch = Scratch::new("residual");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let whole = DELETED.whole(&[&first, &second]);
    let joined = DELETED.score_file(&scratch, "joined", &[&[&first, &second]]);
    let cache = spliceai_cache(&scratch, "cache", "3");
    let changed = second.replace("0.02", "0.03");
    for (name, reference, engine, masked, residual) in [
        (
            "other-entry",
            vec![whole.as_str(), ""],
            vec![first.as_str(), first.as_str()],
            1,
            1,
        ),
        (
            "changed",
            vec![whole.as_str(), ""],
            vec![first.as_str(), changed.as_str()],
            1,
            1,
        ),
        (
            "silent-first",
            vec!["", ""],
            vec![first.as_str(), second.as_str()],
            1,
            1,
        ),
        (
            "unnamed-gene",
            vec![whole.as_str(), "", ""],
            vec![first.as_str(), second.as_str(), second.as_str()],
            2,
            1,
        ),
    ] {
        let report = scored(
            &scratch,
            name,
            "SpliceAI",
            MULTI_TABLE,
            &DELETED.rows(&reference),
            &DELETED.rows(&engine),
            &[
                "--spliceai-table",
                text_of(&joined),
                "--json-cache",
                text_of(&cache),
            ],
        );
        assert_eq!(
            report.passed().aside("SpliceAI", MULTI),
            (masked, residual),
            "{name}"
        );
    }
}

#[test]
fn spliceai_rows_of_a_single_base_record_of_several_entries_are_set_aside() {
    let scratch = Scratch::new("single-base-rows");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let whole = SUBSTITUTED.whole(&[&first, &second]);
    let first_twice = SUBSTITUTED.whole(&[&first, &first]);
    let joined = SUBSTITUTED.score_file(&scratch, "joined", &[&[&first, &second]]);
    let single = SUBSTITUTED.score_file(&scratch, "single", &[&[&first], &[&second]]);
    let cache = spliceai_cache(&scratch, "cache", "3");
    for (name, reference, engine, files, masked, residual) in [
        (
            "several-entries",
            [whole.as_str(), ""],
            [first.as_str(), second.as_str()],
            Some(&joined),
            2,
            0,
        ),
        (
            "no-file",
            [whole.as_str(), ""],
            [first.as_str(), second.as_str()],
            None,
            0,
            3,
        ),
        (
            "one-entry-records",
            [whole.as_str(), ""],
            [first.as_str(), second.as_str()],
            Some(&single),
            0,
            3,
        ),
        (
            "other-entry",
            [whole.as_str(), ""],
            [first.as_str(), first.as_str()],
            Some(&joined),
            1,
            1,
        ),
        (
            "first-twice",
            [first_twice.as_str(), ""],
            [first.as_str(), second.as_str()],
            Some(&joined),
            1,
            2,
        ),
    ] {
        let mut more = vec!["--json-cache", text_of(&cache)];
        if let Some(file) = files {
            more.extend(["--spliceai-table", text_of(file)]);
        }
        let report = scored(
            &scratch,
            name,
            "SpliceAI",
            MULTI_TABLE,
            &SUBSTITUTED.rows(&reference),
            &SUBSTITUTED.rows(&engine),
            &more,
        );
        assert_eq!(
            report.passed().aside("SpliceAI", MULTI),
            (masked, residual),
            "{name}"
        );
    }
}

#[test]
fn spliceai_rows_are_set_aside_under_every_record_line_at_the_allele() {
    let scratch = Scratch::new("every-line-rows");
    let cache = spliceai_cache(&scratch, "cache", "3");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let changed = "GENEA|0.31|0.00|0.44|0.09|3|-2|11|-30".to_string();
    let other = format!("GENEB|{}", scores("GENEA"));
    let third = format!("GENEA|{}", scores("GENEB"));
    let last = "GENEB|0.31|0.00|0.44|0.09|3|-2|11|-30".to_string();
    let both = [first.as_str(), second.as_str()];
    let reversed = [second.as_str(), first.as_str()];
    let whole = |values: &[&str]| SUBSTITUTED.whole(values);
    type Shape<'a> = (
        &'a str,
        Vec<&'a [&'a str]>,
        Vec<String>,
        Vec<&'a str>,
        i64,
        i64,
    );
    let changed_other = [changed.as_str(), other.as_str()];
    let third_last = [third.as_str(), last.as_str()];
    let other_changed = [other.as_str(), changed.as_str()];
    let shapes: Vec<Shape> = vec![
        (
            "one-line",
            vec![&both],
            vec![whole(&both), String::new()],
            both.to_vec(),
            2,
            0,
        ),
        (
            "twice",
            vec![&both, &both],
            vec![whole(&both), String::new()],
            both.to_vec(),
            2,
            0,
        ),
        (
            "three-times",
            vec![&both, &both, &both],
            vec![whole(&both), String::new()],
            both.to_vec(),
            2,
            0,
        ),
        (
            "other-order",
            vec![&both, &reversed],
            vec![whole(&both), whole(&reversed)],
            both.to_vec(),
            2,
            0,
        ),
        (
            "one-entry-last",
            vec![&both, &both[..1]],
            vec![first.clone(), String::new()],
            both.to_vec(),
            1,
            0,
        ),
        (
            "one-entry-first",
            vec![&both[..1], &both],
            vec![whole(&both), String::new()],
            both.to_vec(),
            2,
            0,
        ),
        (
            "gene-unnamed",
            vec![&both, &both],
            vec![whole(&both), String::new(), String::new()],
            vec![&first, &second, &second],
            2,
            1,
        ),
        (
            "other-scores",
            vec![&both, &changed_other],
            vec![whole(&[&changed, &other]), String::new()],
            vec![&changed, &other],
            2,
            0,
        ),
        (
            "other-scores-three-times",
            vec![&both, &changed_other, &third_last],
            vec![whole(&[&third, &last]), String::new()],
            vec![&third, &last],
            2,
            0,
        ),
        (
            "other-scores-other-order",
            vec![&both, &other_changed],
            vec![whole(&both), whole(&[&other, &changed])],
            vec![&changed, &other],
            2,
            0,
        ),
        (
            "earlier-line",
            vec![&both, &changed_other],
            vec![whole(&[&changed, &other]), String::new()],
            vec![&first, &other],
            1,
            2,
        ),
    ];
    for (name, lines, reference, engine, masked, residual) in shapes {
        let file = SUBSTITUTED.score_file(&scratch, &format!("{name}/scores"), &lines);
        let reference: Vec<&str> = reference.iter().map(String::as_str).collect();
        let (reference, engine) = (SUBSTITUTED.rows(&reference), SUBSTITUTED.rows(&engine));
        let cached = ["--json-cache", text_of(&cache)];
        let with = [&cached[..], &["--spliceai-table", text_of(&file)]].concat();
        let report = scored(
            &scratch,
            name,
            "SpliceAI",
            MULTI_TABLE,
            &reference,
            &engine,
            &with,
        );
        assert_eq!(
            report.passed().aside("SpliceAI", MULTI),
            (masked, residual),
            "{name}"
        );
        let report = scored(
            &scratch,
            name,
            "SpliceAI",
            MULTI_TABLE,
            &reference,
            &engine,
            &cached,
        );
        let (masked, residual) = report.passed().aside("SpliceAI", MULTI);
        assert!(masked == 0 && residual > 0, "{name}");
    }
}

#[test]
fn spliceai_class_counts_the_same_refused_rows_on_any_thread_count() {
    use vep_plugin_pairs::score::{report::dumps, score, ClassTable, Options};
    let scratch = Scratch::new("refusals-on-threads");
    let cases = refused_cases("21", 7000);
    let output = |side: usize| {
        let rows: Vec<String> = cases
            .iter()
            .map(|(location, gene, _, reference, engine)| {
                let mut extra = "IMPACT=MODIFIER".to_string();
                if !gene.is_empty() {
                    extra.push_str(&format!(";SYMBOL={gene}"));
                }
                let written = [reference, engine][side];
                if !written.is_empty() {
                    extra.push_str(&format!(";{written}"));
                }
                default_row(
                    &format!("snv_{location}_{gene}"),
                    location,
                    "C",
                    &format!("ENST_{gene}"),
                    &extra,
                )
            })
            .collect();
        scratch.file(&format!("side{side}.txt"), default_text(&rows))
    };
    let (reference, engine) = (output(0), output(1));
    let scores = score_file(
        &scratch,
        "scores/snv.synthetic.vcf.gz",
        &refused_lines(&cases),
    );
    let table = scratch.file("classes.tsv", MULTI_TABLE);
    // Ranges of one byte hold one Location each, so the rows are asked about on as many
    // threads as the run has.
    let report = |threads: i64, range_bytes: Option<u64>| {
        let options = Options {
            class_table: Some(ClassTable::File(text_of(&table).to_string())),
            spliceai_table: vec![text_of(&scores).to_string()],
            pairing_threads: Some(threads),
            settle_range_bytes: range_bytes,
            ..Options::default()
        };
        score(
            text_of(&reference),
            &[("engine".to_string(), text_of(&engine).to_string())],
            &["SpliceAI".to_string()],
            &[],
            Some(text_of(&scratch.path("work"))),
            &options,
            &BuiltInRegistry,
            &BuiltIn,
        )
        .unwrap_or_else(|stop| panic!("the run on {threads} threads was stopped: {stop:?}"))
    };
    let first = report(1, None);
    let grounds: serde_json::Map<String, Value> = REFUSED
        .iter()
        .map(|(ground, rows)| (ground.to_string(), json!(rows)))
        .collect();
    assert_eq!(
        first["spliceai_refusals"],
        json!([{"engine": "engine", "asked": 15, "refused": 13, "grounds": grounds}])
    );
    for threads in [1, 2, 5, 16] {
        assert!(
            dumps(&report(threads, Some(1))) == dumps(&first),
            "{threads} threads"
        );
    }
}

#[test]
fn summary_counts_the_rows_the_spliceai_class_refuses_by_ground_under_the_score_files() {
    let scratch = Scratch::new("summary-grounds");
    // The row of a gene the cache does not hold stands on ENST9, the row of a gene no
    // entry names on ENST3, and the last row on ENST2.
    let cases = refused_cases("3", 7_000_000);
    let features = [
        "ENST1", "ENST1", "ENST1", "ENST1", "ENST1", "ENST1", "ENST9", "ENST3", "ENST1", "ENST1",
        "ENST1", "ENST1", "ENST1", "ENST1", "ENST2",
    ];
    let rows = |side: usize| -> Vec<String> {
        cases
            .iter()
            .zip(features)
            .enumerate()
            .map(|(at, ((location, _, _, reference, engine), feature))| {
                let written = [reference, engine][side];
                let pairs = if written.is_empty() {
                    String::new()
                } else {
                    format!(";{written}")
                };
                row(
                    location,
                    "C",
                    feature,
                    &format!("IMPACT=MODIFIER{pairs}"),
                    &format!("snv_{at}"),
                )
            })
            .collect()
    };
    let (reference, engine) = (rows(0), rows(1));
    let file = score_file(
        &scratch,
        "scores/snv.synthetic.vcf.gz",
        &refused_lines(&cases),
    );
    let cache = spliceai_cache(&scratch, "cache", "3");
    let cached = ["--json-cache", text_of(&cache)];
    let with = [&cached[..], &["--spliceai-table", text_of(&file)]].concat();
    let report = scored(
        &scratch,
        "with",
        "SpliceAI",
        MULTI_TABLE,
        &reference,
        &engine,
        &with,
    );
    let mut counted = vec![
        "spliceai_multi_entry_cell under --spliceai-table, engine engine: rows asked about 15, rows refused 13"
            .to_string(),
    ];
    for ((ground, rows), (key, read)) in REFUSED.iter().zip(GROUNDS) {
        assert_eq!(*ground, key);
        counted.push(format!(
            "spliceai_multi_entry_cell, engine engine: rows refused {rows}: {read}"
        ));
    }
    let lines = &report.passed().lines;
    let (summary, last) = lines.split_at(lines.len() - counted.len());
    assert_eq!(last, counted);
    assert!(counted[1].ends_with(": rows refused 1: the candidate's row writes no SpliceAI value"));
    assert!(counted[5].ends_with(
        ": rows refused 2: no record line of the score files stands at the row's variant"
    ));
    assert!(counted[7].ends_with(
        ": rows refused 1: no entry of a record line at the row's variant names the row's gene"
    ));
    assert!(counted[9].ends_with(
        ": rows refused 1: the candidate's value is not the gene's entries of the last record line that names it"
    ));
    for word in ["GENE", "ENST", "7000", "0.1", "snv_"] {
        assert!(!counted.iter().any(|line| line.contains(word)), "{word}");
    }
    assert!(summary.contains(
        &"SpliceAI\tengine\tfeature\tclass spliceai_multi_entry_cell set aside 2 [met, -]"
            .to_string()
    ));
    assert!(!summary
        .iter()
        .any(|line| line.contains("rows refused") || line.contains("rows asked about")));
    let grounds: serde_json::Map<String, Value> = REFUSED
        .iter()
        .map(|(ground, rows)| (ground.to_string(), json!(rows)))
        .collect();
    assert_eq!(
        report.report["spliceai_refusals"],
        json!([{"engine": "engine", "asked": 15, "refused": 13, "grounds": grounds}])
    );
    // Without the score files nothing is counted or printed.
    let without = scored(
        &scratch,
        "without",
        "SpliceAI",
        MULTI_TABLE,
        &reference,
        &engine,
        &cached,
    );
    assert_eq!(without.passed().lines.last().unwrap(), SPLICEAI_NOTICE);
    assert!(without.report["spliceai_refusals"].is_null());
    assert!(!without
        .lines
        .iter()
        .any(|line| line.contains("rows refused") || line.contains("rows asked about")));
    // Rows the class takes whole print the one line.
    let taken = scored(
        &scratch,
        "taken",
        "SpliceAI",
        MULTI_TABLE,
        &reference[13..],
        &engine[13..],
        &with,
    );
    assert_eq!(
        taken.passed().lines.last().unwrap(),
        "spliceai_multi_entry_cell under --spliceai-table, engine engine: rows asked about 2, rows refused 0"
    );
    assert_eq!(
        taken.report["spliceai_refusals"],
        json!([{"engine": "engine", "asked": 2, "refused": 0, "grounds": {}}])
    );
    // An engine that writes the reference's rows asks the class about none.
    let equal = scored(
        &scratch,
        "equal",
        "SpliceAI",
        MULTI_TABLE,
        &reference,
        &reference,
        &with,
    );
    assert_eq!(equal.passed().report["spliceai_refusals"], json!([]));
    assert!(!equal
        .lines
        .iter()
        .any(|line| line.contains("rows asked about")));
}

#[test]
fn summary_says_that_the_spliceai_class_sets_nothing_aside_without_the_score_files() {
    let scratch = Scratch::new("summary-files");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let joined = DELETED.score_file(&scratch, "joined", &[&[&first, &second]]);
    let single = DELETED.score_file(&scratch, "single", &[&[&first], &[&second]]);
    let whole = DELETED.whole(&[&first, &second]);
    let reference = DELETED.rows(&[&whole, ""]);
    let engine = DELETED.rows(&[&first, &second]);
    let cache = spliceai_cache(&scratch, "cache", "3");
    let run = |name: &str, files: &[&Path]| {
        let more = [&["--json-cache", text_of(&cache)][..], &table_words(files)].concat();
        scored(
            &scratch,
            name,
            "SpliceAI",
            MULTI_TABLE,
            &reference,
            &engine,
            &more,
        )
    };
    let set_aside = |keys: &str| {
        format!(
            "SpliceAI\tengine\tfeature\tclass spliceai_multi_entry_cell set aside {keys} [met, -]"
        )
    };
    let without = run("without", &[]);
    assert_eq!(without.passed().lines.last().unwrap(), SPLICEAI_NOTICE);
    assert!(without.lines.contains(&set_aside("0 (none)")));
    // Under both files the last line that names a gene decides: the records of one entry.
    let with = run("with", &[&joined, &single]);
    assert!(!with
        .passed()
        .lines
        .iter()
        .any(|line| line == SPLICEAI_NOTICE));
    assert!(!with.lines.contains(&set_aside("2")));
    let alone = run("joined", &[&joined]);
    assert!(alone.passed().lines.contains(&set_aside("2")));
    assert_eq!(
        alone.report["spliceai_table"],
        json!([{"path": text_of(&joined)}])
    );
    let nowhere = scratch.path("nowhere.synthetic.vcf.gz");
    let missing = run("missing", &[&nowhere]);
    assert_eq!(
        (missing.code, missing.lines.len(), missing.error.as_str()),
        (
            2,
            0,
            format!("{ERROR}not a file: {}\n", nowhere.display()).as_str()
        )
    );
    fs::remove_file(format!("{}.tbi", joined.display())).unwrap();
    let unindexed = run("unindexed", &[&joined]);
    assert_eq!((unindexed.code, unindexed.lines.len()), (2, 0));
    assert!(unindexed.error.starts_with(&format!(
        "{ERROR}--spliceai-table: cannot read the tabix index of {joined}: none of {joined}.csi, ",
        joined = joined.display()
    )));
}

#[test]
fn summary_says_that_the_spliceai_class_names_no_gene_from_a_cache_without_one() {
    let scratch = Scratch::new("summary-cache");
    let (first, second) = (own_value("GENEA"), own_value("GENEB"));
    let joined = DELETED.score_file(&scratch, "joined", &[&[&first, &second]]);
    let whole = DELETED.whole(&[&first, &second]);
    let reference = DELETED.rows(&[&whole, ""]);
    let engine = DELETED.rows(&[&first, &second]);
    let cache = spliceai_cache(&scratch, "cache", "3");
    let run = |name: &str, more: &[&str]| {
        scored(
            &scratch,
            name,
            "SpliceAI",
            MULTI_TABLE,
            &reference,
            &engine,
            more,
        )
    };
    let set_aside = |keys: &str| {
        format!(
            "SpliceAI\tengine\tfeature\tclass spliceai_multi_entry_cell set aside {keys} [met, -]"
        )
    };
    // The notice stands before the counts of the rows the class refused, both for want of
    // a named gene.
    let table_alone = run("table-alone", &["--spliceai-table", text_of(&joined)]);
    let lines = &table_alone.passed().lines;
    assert_eq!(
        lines[lines.len() - 3..],
        [
            SPLICEAI_CACHE_NOTICE.to_string(),
            "spliceai_multi_entry_cell under --spliceai-table, engine engine: rows asked about 2, rows refused 2".to_string(),
            "spliceai_multi_entry_cell, engine engine: rows refused 2: neither a SYMBOL nor the candidate cache names the row's gene".to_string(),
        ]
    );
    assert!(lines.contains(&set_aside("0 (none)")));
    let both = run(
        "both",
        &[
            "--spliceai-table",
            text_of(&joined),
            "--json-cache",
            text_of(&cache),
        ],
    );
    assert!(!both
        .passed()
        .lines
        .iter()
        .any(|line| line == SPLICEAI_CACHE_NOTICE));
    assert!(both.lines.contains(&set_aside("2")));
    let neither = run("neither", &[]);
    assert_eq!(neither.passed().lines.last().unwrap(), SPLICEAI_NOTICE);
    assert!(!neither
        .lines
        .iter()
        .any(|line| line == SPLICEAI_CACHE_NOTICE));
}

/// Per SpliceAI corpus and set (`base` the manifest's own run): the keys under which the
/// class sets a tuple aside in the feature projection and in the allele projection.
const CORPUS_COUNTS: [(&str, &str, i64, i64); 8] = [
    ("GRCh37-plugin-spliceai-fixture", "base", 32, 8),
    ("GRCh37-plugin-spliceai-fixture", "cutoff", 52, 15),
    ("GRCh37-plugin-spliceai-fixture", "split-output", 246, 99),
    (
        "GRCh37-plugin-spliceai-fixture",
        "split-output-cutoff",
        266,
        106,
    ),
    ("GRCh38-plugin-spliceai-fixture", "base", 32, 8),
    ("GRCh38-plugin-spliceai-fixture", "cutoff", 48, 16),
    ("GRCh38-plugin-spliceai-fixture", "split-output", 241, 99),
    (
        "GRCh38-plugin-spliceai-fixture",
        "split-output-cutoff",
        257,
        107,
    ),
];

/// The engine's default-format output of a corpus run as the corpus's committed files
/// state it: the reference's expected rows with each field difference the run's manifest
/// entry documents under a compared field written as the engine's value.
fn corpus_engine_output(target: &Value, expected: &Path) -> String {
    let compared: Vec<&str> = target["fields_compared"]
        .as_array()
        .unwrap()
        .iter()
        .map(|field| field.as_str().unwrap())
        .collect();
    let mut documented: BTreeMap<[&str; 3], Vec<(&str, &str)>> = BTreeMap::new();
    for item in target["field_divergences"].as_array().unwrap() {
        let text = |key: &str| item[key].as_str().unwrap();
        if compared.contains(&text("field")) {
            documented
                .entry([text("location"), text("allele"), text("feature")])
                .or_default()
                .push((text("field"), text("vep_rs")));
        }
    }
    let mut text = String::new();
    MultiGzDecoder::new(fs::File::open(expected).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let mut out = String::new();
    for line in text.lines() {
        let mut cells: Vec<String> = line.split('\t').map(str::to_string).collect();
        let fields = if line.starts_with('#') {
            None
        } else {
            documented.get(&[cells[1].as_str(), cells[2].as_str(), cells[4].as_str()])
        };
        if let Some(fields) = fields {
            let mut pairs: Vec<(String, String)> = Vec::new();
            if cells[13] != "-" {
                for pair in cells[13].split(';') {
                    let (name, value) = pair.split_once('=').unwrap_or((pair, ""));
                    pairs.push((name.to_string(), value.to_string()));
                }
            }
            for (name, value) in fields {
                match pairs.iter_mut().find(|(held, _)| held == name) {
                    Some(pair) => pair.1 = value.to_string(),
                    None => pairs.push((name.to_string(), value.to_string())),
                }
            }
            let joined: Vec<String> = pairs
                .iter()
                .map(|(name, value)| format!("{name}={value}"))
                .collect();
            cells[13] = joined.join(";");
        }
        out.push_str(&cells.join("\t"));
        out.push('\n');
    }
    out
}

#[test]
fn the_spliceai_corpora_apply_the_class_of_a_record_of_several_entries_at_their_count() {
    let scratch = Scratch::new("corpora");
    let rows: Vec<&str> = COMMITTED_CLASS_TABLE
        .split_inclusive('\n')
        .filter(|line| line.split('\t').nth(1) == Some("SpliceAI"))
        .collect();
    let met = scratch.file("met.tsv", rows.concat());
    let unmet: Vec<String> = rows
        .iter()
        .map(|line| {
            if line.starts_with(MULTI) {
                line.replace("\tmet\t", "\tnot met\t")
            } else {
                line.to_string()
            }
        })
        .collect();
    let unmet = scratch.file("unmet.tsv", unmet.concat());
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116");
    for (corpus, name, feature, allele) in CORPUS_COUNTS {
        let directory = golden.join(corpus);
        let manifest: Value =
            serde_json::from_slice(&fs::read(directory.join("manifest.json")).unwrap()).unwrap();
        let target = if name == "base" {
            &manifest
        } else {
            manifest["flag_sets"]
                .as_array()
                .unwrap()
                .iter()
                .find(|entry| entry["name"] == name)
                .unwrap()
        };
        let expected_dir = if name == "base" {
            "expected"
        } else {
            target["expected_dir"].as_str().unwrap()
        };
        let reference = directory.join(expected_dir).join("default.txt.gz");
        let engine = scratch.file(
            &format!("{corpus}/{name}/vep-rs.txt"),
            corpus_engine_output(target, &reference),
        );
        let mut files: Vec<PathBuf> = Vec::new();
        for kind in ["snv", "indel"] {
            let mut of_kind: Vec<PathBuf> = fs::read_dir(directory.join("plugin_data"))
                .unwrap()
                .map(|entry| entry.unwrap().path())
                .filter(|path| {
                    let name = path.file_name().unwrap().to_str().unwrap();
                    name.starts_with(&format!("spliceai_scores.raw.{kind}."))
                        && name.ends_with(".vcf.gz")
                })
                .collect();
            of_kind.sort();
            files.extend(of_kind);
        }
        let files: Vec<&Path> = files.iter().map(PathBuf::as_path).collect();
        let cache = directory.join("json_cache");
        let scored = |table: &Path, label: &str| {
            let work = scratch.path(&format!("{corpus}/{name}/work-{label}"));
            let report = scratch.path(&format!("{corpus}/{name}/{label}.json"));
            let engine = format!("vep-rs={}", engine.display());
            let mut line = vec![
                "--reference",
                text_of(&reference),
                "--engine",
                &engine,
                "--plugins",
                "SpliceAI",
                "--class-table",
                text_of(table),
                "--json-cache",
                text_of(&cache),
                "--work-dir",
                text_of(&work),
                "--report",
                text_of(&report),
            ];
            line.extend(table_words(&files));
            let (code, _, error) = run_main(&line);
            assert_eq!((code, error.as_str()), (0, ""), "{corpus} {name}");
            let report: Value = serde_json::from_slice(&fs::read(&report).unwrap()).unwrap();
            report
        };
        let (applied, left) = (scored(&met, "met"), scored(&unmet, "unmet"));
        let entry = |report: &Value, projection: &str| -> Value {
            report["scores"]
                .as_array()
                .unwrap()
                .iter()
                .find(|entry| entry["source"] == "SpliceAI" && entry["projection"] == projection)
                .unwrap()
                .clone()
        };
        for (projection, count) in [("feature", feature), ("allele", allele)] {
            let item = entry(&applied, projection);
            assert_eq!(item["masked"][MULTI], count, "{corpus} {name} {projection}");
            assert_eq!(item["adj_perl"], item["adj_rust"]);
            assert_eq!(item["adj_perl"], item["adj_intersection"]);
            assert_eq!(item["residual_rows"], 0, "{corpus} {name} {projection}");
            let item = entry(&left, projection);
            assert_eq!(
                item["not_applied"],
                json!({MULTI: count}),
                "{corpus} {name} {projection}"
            );
            assert!(item["residual_rows"].as_i64().unwrap() > 0);
        }
    }
}

#[test]
fn summary_says_that_classes_of_gwas_set_nothing_aside_without_the_catalog() {
    let scratch = Scratch::new("summary-catalog");
    let catalog = scratch.file(
        "associations.tsv",
        format!("{CATALOG_HEADER}\n21;21\t100;150\trs1-G; rs2-A\t-\tGCST1\n"),
    );
    let classes = "gwas_haplotype_cross_product\tGWAS\tGWAS_*\tmet\t-\n";
    let reference = [
        row("21:100", "A", "ENST1", "GWAS_pmid=17;GWAS_study=GCST1", "v"),
        row("21:400", "A", "ENST1", "GWAS_pmid=18;GWAS_study=GCST1", "v"),
        row("21:500", "T", "ENST1", "IMPACT=LOW", "v"),
    ];
    let engine = [
        row("21:100", "A", "ENST1", "IMPACT=LOW", "v"),
        row("21:400", "A", "ENST1", "IMPACT=LOW", "v"),
        row("21:500", "T", "ENST1", "IMPACT=LOW", "v"),
    ];
    let run = |name: &str, more: &[&str]| {
        scored(&scratch, name, "GWAS", classes, &reference, &engine, more)
    };
    let notice = "no --gwas-catalog: a class of GWAS that reads the catalog sets nothing aside";
    let set_aside = |keys: &str| {
        format!(
            "GWAS\tengine\tfeature\tclass gwas_haplotype_cross_product set aside {keys} [met, -]"
        )
    };
    let without = run("without", &[]);
    assert_eq!(without.passed().lines.last().unwrap(), notice);
    assert!(without.lines.contains(&set_aside("0 (none)")));
    assert!(without.report["gwas_catalog"].is_null());
    assert_eq!(without.aside("GWAS", CROSS_PRODUCT), (0, 4));
    let with = run("with", &["--gwas-catalog", text_of(&catalog)]);
    assert!(!with.passed().lines.iter().any(|line| line == notice));
    assert!(with.lines.contains(&set_aside("2")));
    assert_eq!(with.aside("GWAS", CROSS_PRODUCT), (2, 2));
    let digest: String = Sha256::digest(fs::read(&catalog).unwrap())
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    assert_eq!(
        with.report["gwas_catalog"],
        json!({"path": text_of(&catalog), "sha256": digest})
    );
    let nowhere = scratch.path("nowhere.tsv");
    let missing = run("missing", &["--gwas-catalog", text_of(&nowhere)]);
    assert_eq!(
        (missing.code, missing.lines.len(), missing.error.as_str()),
        (
            2,
            0,
            format!("{ERROR}not a file: {}\n", nowhere.display()).as_str()
        )
    );
    let other = scratch.file("other.tsv", "a\tb\n");
    let refused = run("refused", &["--gwas-catalog", text_of(&other)]);
    assert_eq!(
        (refused.code, refused.lines.len(), refused.error.as_str()),
        (
            2,
            0,
            format!("{ERROR}{}: not a GWAS catalog associations file (columns CHR_ID, CHR_POS, STRONGEST SNP-RISK ALLELE, STUDY ACCESSION)\n", other.display()).as_str()
        )
    );
}
