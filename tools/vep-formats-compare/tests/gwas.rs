// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two GWAS catalog classes and the classes of a field the candidate alone fills, read
//! through the comparison of two entry lists on planted rows and a planted catalog.

mod compare_common;

use compare_common::*;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::path::Path;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Uncaught;
use vep_formats_compare::fields::entry_kind;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::field_class_accepts;
use vep_formats_compare::rules::gwas::gwas_site;
use vep_formats_compare::rules::gwas::GwasCatalogPairs;
use vep_formats_compare::rules::gwas::GwasPair;
use vep_formats_compare::rules::RuleContext;
use vep_formats_compare::rules::Side;

const CROSS_PRODUCT: &str = "gwas_haplotype_cross_product";
const INTERACTION: &str = "gwas_interaction_first_snp_allele_whitespace";

const GWAS_CATALOG_HEADER: &str =
    "CHR_ID\tCHR_POS\tSTRONGEST SNP-RISK ALLELE\tSNPS\tSTUDY ACCESSION";
const GWAS_ROW_FIELDS: &str =
    "GWAS_associations=rs1;GWAS_risk_allele=rs1-A;GWAS_p_value=1e-8;GWAS_study=GCST1";

const STUDY_CSQ_HEADER: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|GWAS_study\">\n";

/// A catalog row: its CHR_ID, CHR_POS, STRONGEST SNP-RISK ALLELE and STUDY ACCESSION cells.
type CatalogRow<'a> = (&'a str, &'a str, &'a str, &'a str);

fn pair(chromosome: &str, position: &str, allele: &str, study: &str) -> GwasPair {
    (
        chromosome.to_string(),
        position.to_string(),
        allele.to_string(),
        study.to_string(),
    )
}

fn pairs_of(items: &[CatalogRow]) -> BTreeSet<GwasPair> {
    items
        .iter()
        .map(|(chromosome, position, allele, study)| pair(chromosome, position, allele, study))
        .collect()
}

/// The pairs of a catalog of `rows`, written under `directory`.
fn gwas_pairs(directory: &Path, rows: &[CatalogRow]) -> GwasCatalogPairs {
    let mut lines = vec![GWAS_CATALOG_HEADER.to_string()];
    lines.extend(rows.iter().map(|(chromosomes, positions, alleles, study)| {
        format!("{chromosomes}\t{positions}\t{alleles}\t-\t{study}")
    }));
    let catalog = planted(directory, "gwas_catalog.tsv", &(lines.join("\n") + "\n"));
    GwasCatalogPairs::from_catalog(&catalog).expect("the catalog reads")
}

fn side_of(pairs: GwasCatalogPairs) -> Side {
    Side {
        gwas_pairs: Some(pairs),
        ..Side::default()
    }
}

/// The side tables of a comparison that reads the catalog of `rows`.
fn gwas_side(directory: &Path, rows: &[CatalogRow]) -> Side {
    side_of(gwas_pairs(directory, rows))
}

fn pairs(side: &Side) -> &GwasCatalogPairs {
    side.gwas_pairs.as_ref().expect("the side names a catalog")
}

/// The report of two entry lists under the built-in `classes`, and what was applied.
fn compare_under(
    expected: &[Entry],
    actual: &[Entry],
    side: &Side,
    classes: &[&str],
) -> (Value, Documented) {
    let mut documented = documented_with(classes);
    let report = compare_entries(
        expected,
        actual,
        &mut documented,
        &CompareOptions::new(side),
    )
    .expect("the entries compare");
    (report, documented)
}

/// The report of two default-format texts read beside the catalog of `side`.
fn compare_gwas(
    reference: &str,
    candidate: &str,
    side: &Side,
    classes: &[&str],
) -> (Value, Documented) {
    let expected = parse_default(reference).expect("the reference parses");
    let actual = parse_default(candidate).expect("the candidate parses");
    compare_under(&expected, &actual, side, classes)
}

fn gwas_row(name: &str, location: &str, allele: &str, fields: &str) -> String {
    let extra = if fields.is_empty() {
        "IMPACT=MODIFIER;STRAND=1".to_string()
    } else {
        format!("IMPACT=MODIFIER;STRAND=1;{fields}")
    };
    default_row(name, location, allele, "ENST1", "intron_variant", &extra)
}

fn applied(documented: &Documented) -> BTreeMap<String, i64> {
    counts(&documented.applied)
}

/// The classes a column's values were set aside under, each with its count.
fn masked(report: &Value, column: &str) -> BTreeMap<String, i64> {
    at(report, &["columns", column, "masked"])
        .as_dict()
        .expect("the masked counts")
        .iter()
        .map(|(klass, counted)| match counted {
            Value::Int(counted) => (klass.clone(), *counted),
            other => panic!("a count is {other:?}"),
        })
        .collect()
}

fn differing(report: &Value, column: &str) -> i64 {
    number(report, &["columns", column, "differing"])
}

fn lines_of(lines: &[&str]) -> Vec<String> {
    lines.iter().map(|line| line.to_string()).collect()
}

fn refs(entries: &[Entry]) -> Vec<&Entry> {
    entries.iter().collect()
}

/// `line` without its last character, the empty INFO column a CSQ value is written into.
fn before_info(line: &str) -> &str {
    &line[..line.len() - 1]
}

#[test]
fn the_cross_product_class_accepts_the_fields_the_reference_alone_fills() {
    // Two haplotype rows of the catalog, each two SNPs on chromosome 21: at the first SNP's
    // location the reference annotates the second SNP's risk allele too. A GWAS_ field the
    // reference fills there and the candidate leaves empty, or both fill differently, is
    // the class; a field the candidate alone fills stays a difference, and the class is
    // counted once per entry however many of its fields differ.
    let directory = scratch("gwas-reference-filled");
    let side = gwas_side(
        &directory,
        &[
            ("21;21", "100;150", "rs1-G; rs2-A", "GCST1"),
            ("21;21", "200;250", "rs3-T; rs4-C", "GCST1"),
        ],
    );
    let held = pairs(&side);
    assert!(held
        .reference_only
        .contains(&pair("21", "100", "A", "GCST1")));
    assert!(held
        .reference_only
        .contains(&pair("21", "150", "G", "GCST1")));
    assert!(held.candidate_only.is_empty());
    let second = GWAS_ROW_FIELDS.replace("rs1", "rs2");
    let reference = default_text(&[
        gwas_row("a", "21:100", "A", GWAS_ROW_FIELDS),
        gwas_row("b", "21:200", "C", &second),
        gwas_row("c", "21:300", "G", ""),
    ]);
    let candidate = default_text(&[
        gwas_row("a", "21:100", "A", ""),
        gwas_row(
            "b",
            "21:200",
            "C",
            &second.replace("GCST1", "GCST2").replace("1e-8", "3e-9"),
        ),
        gwas_row("c", "21:300", "G", ""),
    ]);
    let (report, documented) = compare_gwas(&reference, &candidate, &side, &[CROSS_PRODUCT]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(CROSS_PRODUCT, 2)]));
    assert_eq!(
        masked(&report, "GWAS_study"),
        count_map(&[(CROSS_PRODUCT, 2)])
    );
    assert_eq!(
        masked(&report, "GWAS_risk_allele"),
        count_map(&[(CROSS_PRODUCT, 1)])
    );
    assert_eq!(
        verdict(&compare_gwas(&reference, &candidate, &side, &[]).0),
        "red"
    );
    let candidate_only = candidate.replace(
        &gwas_row("c", "21:300", "G", ""),
        &gwas_row("c", "21:300", "G", GWAS_ROW_FIELDS),
    );
    assert_eq!(
        verdict(&compare_gwas(&reference, &candidate_only, &side, &[CROSS_PRODUCT]).0),
        "red"
    );
    // A field of another plugin differing on the same row is not the class's.
    let other_plugin = candidate.replacen(
        "IMPACT=MODIFIER;STRAND=1\n",
        "IMPACT=MODIFIER;STRAND=1;LoFtool=0.5\n",
        1,
    );
    assert_eq!(
        verdict(&compare_gwas(&reference, &other_plugin, &side, &[CROSS_PRODUCT]).0),
        "red"
    );
    // The G row of the line 21 99 CA CG,C starts at 21:100, one base past the line's POS,
    // where the default format's Location and the JSON object's start put it: the VCF entry
    // is read there, and a catalog site at the POS itself pairs nothing with it.
    let inputs = lines_of(&["21\t99\tshared_base\tCA\tCG,C\t.\t.\t."]);
    let vcf = format!(
        "{STUDY_CSQ_HEADER}{}CSQ=G|intron_variant|ENST1|Transcript|GCST1,-|intron_variant|ENST1|Transcript|\n",
        before_info(&inputs[0])
    );
    for (position, wanted, counted) in [
        ("100", "green", count_map(&[(CROSS_PRODUCT, 1)])),
        ("99", "red", count_map(&[])),
    ] {
        let site = gwas_side(
            &directory,
            &[("21;21", &format!("{position};150"), "rs1-T; rs2-G", "GCST1")],
        );
        let expected = parse_vcf(&vcf, &inputs).expect("the reference parses");
        let actual =
            parse_vcf(&vcf.replace("|GCST1,", "|,"), &inputs).expect("the candidate parses");
        let (report, documented) = compare_under(&expected, &actual, &site, &[CROSS_PRODUCT]);
        assert_eq!(verdict(&report), wanted, "{position}");
        assert_eq!(applied(&documented), counted, "{position}");
    }
}

#[test]
fn the_gwas_classes_read_no_row_of_an_insertion_at_its_records_start_in_the_vcf_and_json_formats() {
    // The default and tab formats read a row at a site only where its Location spans one
    // base or as many bases as its sequence allele, so an insertion's row, printed over its
    // two flanking bases, states no site there. The VCF and JSON formats hold a row to the
    // same condition on the span its record states: an insertion's row, whose span is
    // empty, stays a difference in the three formats, where a substitution's row at the
    // same site and under the same allele is the class in each.
    let directory = scratch("gwas-insertion-rows");
    let side = gwas_side(&directory, &[("21;21", "100;150", "rs1-G; rs2-A", "GCST1")]);
    assert!(pairs(&side)
        .reference_only
        .contains(&pair("21", "100", "A", "GCST1")));
    let one = count_map(&[(CROSS_PRODUCT, 1)]);
    let none = count_map(&[]);
    let study = "GWAS_study=GCST1";
    for (name, location, wanted, counted) in [
        ("sub", "21:100", "green", &one),
        ("ins", "21:99-100", "red", &none),
    ] {
        let (report, documented) = compare_gwas(
            &default_text(&[gwas_row(name, location, "A", study)]),
            &default_text(&[gwas_row(name, location, "A", "")]),
            &side,
            &[CROSS_PRODUCT],
        );
        assert_eq!(verdict(&report), wanted, "{name}");
        assert_eq!(&applied(&documented), counted, "{name}");
    }
    let site = Some(("21".to_string(), "100".to_string()));
    for (line, wanted, counted) in [
        ("21\t100\tsub\tC\tA\t.\t.\t.", "green", &one),
        ("21\t99\tins\tC\tCA\t.\t.\t.", "red", &none),
    ] {
        let vcf = format!(
            "{STUDY_CSQ_HEADER}{}CSQ=A|intron_variant|ENST1|Transcript|GCST1\n",
            before_info(line)
        );
        let inputs = lines_of(&[line]);
        let entries = parse_vcf(&vcf, &inputs).expect("the reference parses");
        let actual = parse_vcf(&vcf.replace("|GCST1", "|"), &inputs).expect("the candidate parses");
        let (report, documented) = compare_under(&entries, &actual, &side, &[CROSS_PRODUCT]);
        assert_eq!(verdict(&report), wanted, "{line}");
        assert_eq!(&applied(&documented), counted, "{line}");
        let row = entries
            .iter()
            .find(|e| entry_kind(e) != "record")
            .expect("a row");
        let stated = if wanted == "green" {
            site.clone()
        } else {
            None
        };
        assert_eq!(gwas_site(row, Some(&refs(&entries))), stated, "{line}");
    }
    let entry = r#""variant_allele": "A", "transcript_id": "ENST1", "consequence_terms": ["intron_variant"], "impact": "MODIFIER", "strand": 1"#;
    for (line, id, span, wanted, counted) in [
        (
            "21\t100\tsub\tC\tA\t.\t.\t.",
            "sub",
            r#""start": 100, "end": 100, "allele_string": "C/A""#,
            "green",
            &one,
        ),
        (
            "21\t99\tins\tC\tCA\t.\t.\t.",
            "ins",
            r#""start": 100, "end": 99, "allele_string": "-/A""#,
            "red",
            &none,
        ),
    ] {
        let record = format!(
            r#""input": {}, "id": "{id}", "seq_region_name": "21", {span}, "most_severe_consequence": "intron_variant""#,
            quoted(line)
        );
        let reference_json = format!(
            "{{{record}, \"transcript_consequences\": [{{{entry}, \"gwas_study\": \"GCST1\"}}]}}\n"
        );
        let candidate_json = format!("{{{record}, \"transcript_consequences\": [{{{entry}}}]}}\n");
        let inputs = lines_of(&[line]);
        let entries = parse_json(&reference_json, &inputs).expect("the reference parses");
        let actual = parse_json(&candidate_json, &inputs).expect("the candidate parses");
        let (report, documented) = compare_under(&entries, &actual, &side, &[CROSS_PRODUCT]);
        assert_eq!(verdict(&report), wanted, "{line}");
        assert_eq!(&applied(&documented), counted, "{line}");
        let row = entries
            .iter()
            .find(|e| entry_kind(e) != "record")
            .expect("a row");
        let stated = if wanted == "green" {
            site.clone()
        } else {
            None
        };
        assert_eq!(gwas_site(row, Some(&refs(&entries))), stated, "{line}");
    }
}

/// The JSON text of `text` as a string.
fn quoted(text: &str) -> String {
    dumped(&format!(
        "\"{}\"",
        text.replace('\\', "\\\\")
            .replace('"', "\\\"")
            .replace('\t', "\\t")
    ))
}

#[test]
fn the_cross_product_class_refuses_a_row_the_catalog_does_not_pair() {
    // The same pair of values stays a difference without the catalog, at a location no row
    // of two or more SNPs names, under an allele that is no other SNP's risk allele, and
    // under a study the catalog row does not carry; a value planted missing or wrong on
    // such a row is counted under its field.
    let directory = scratch("gwas-unpaired-row");
    let side = gwas_side(
        &directory,
        &[
            ("21;21", "100;150", "rs1-G; rs2-A", "GCST1"),
            ("21", "400", "rs9-A", "GCST1"),
        ],
    );
    let no_catalog = Side::default();
    let reference = default_text(&[gwas_row("a", "21:100", "A", GWAS_ROW_FIELDS)]);
    let candidate = default_text(&[gwas_row("a", "21:100", "A", "")]);
    assert_eq!(
        verdict(&compare_gwas(&reference, &candidate, &side, &[CROSS_PRODUCT]).0),
        "green"
    );
    let (report, documented) = compare_gwas(&reference, &candidate, &no_catalog, &[CROSS_PRODUCT]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), count_map(&[]));
    for (old, new) in [
        ("21:100", "21:400"),
        ("\tA\t", "\tG\t"),
        ("GCST1", "GCST7"),
        ("21:100", "21:100-101"),
    ] {
        let (report, documented) = compare_gwas(
            &reference.replace(old, new),
            &candidate.replace(old, new),
            &side,
            &[CROSS_PRODUCT],
        );
        assert_eq!(verdict(&report), "red", "{old} {new}");
        assert_eq!(applied(&documented), count_map(&[]), "{old} {new}");
    }
    let ordinary = default_text(&[gwas_row("a", "21:400", "A", GWAS_ROW_FIELDS)]);
    for planted in [
        default_text(&[gwas_row("a", "21:400", "A", "")]),
        ordinary.replace("GCST1", "GCST2"),
    ] {
        let (report, documented) = compare_gwas(&ordinary, &planted, &side, &[CROSS_PRODUCT]);
        assert_eq!(verdict(&report), "red");
        assert_eq!(applied(&documented), count_map(&[]));
        assert_eq!(differing(&report, "GWAS_study"), 1);
    }
    // A catalog row whose chromosomes and positions differ in count pairs nothing,
    // whichever list is the longer.
    for (chromosomes, positions, alleles) in [
        ("21;21;21", "100;150", "rs1-G; rs2-A"),
        ("21;21", "100;150;175", "rs1-G; rs2-A"),
    ] {
        let uneven = gwas_pairs(&directory, &[(chromosomes, positions, alleles, "GCST1")]);
        assert!(
            uneven.reference_only.is_empty(),
            "{chromosomes} {positions}"
        );
        assert!(
            uneven.candidate_only.is_empty(),
            "{chromosomes} {positions}"
        );
    }
}

const GWAS_SKIPPED_ROW_FIELDS: &str = "GWAS_accessions=http://purl.obolibrary.org/obo/MONDO_0005086;GWAS_associated_gene=MYC,PVT1;GWAS_beta_coef=;GWAS_odds_ratio=1.27;GWAS_p_value=5E-11;GWAS_pmid=24220699;GWAS_risk_allele=CG;GWAS_study=GCST002273";

/// A reference row that carries `fields` and the candidate's row without them.
fn planted_pair(location: &str, allele: &str, fields: &str, name: &str) -> (String, String) {
    (
        default_text(&[gwas_row(name, location, allele, fields)]),
        default_text(&[gwas_row(name, location, allele, "")]),
    )
}

#[test]
fn the_cross_product_class_pairs_a_row_whose_snps_and_locations_differ_in_count() {
    // The catalog row of study GCST002273 names three SNPs over two locations of chromosome
    // 8. The candidate skips such a row and Ensembl VEP writes every SNP at every location
    // of it, so a two-base substitution over both locations whose allele is the first SNP's
    // two-base risk allele carries the row's fields in the reference alone: the class,
    // counted once. So are that allele against one reference base at the first location and
    // from the second location on, and the third SNP's allele at the first location. An
    // allele no SNP of the row names, a site off the row, another study, an allele shorter
    // than its Location, an insertion the record's own alleles spell and a value the
    // candidate alone writes stay differences; a row whose SNPs and locations agree in
    // count pairs each allele away from its own location alone, and one whose risk alleles
    // differ from them pairs nothing. A row whose risk-allele cell holds a ? is written by
    // neither engine and pairs nothing, whatever its counts.
    let directory = scratch("gwas-uneven-snps");
    let lists = "8;8\t127877125;127877126\trs35252396-CG; rs6470588-C; rs6470589-G";
    let catalog = planted(
        &directory,
        "gwas_catalog.tsv",
        &format!("{GWAS_CATALOG_HEADER}\n{lists}\trs35252396; rs6470588; rs6470589\tGCST002273\n"),
    );
    let read = GwasCatalogPairs::from_catalog(&catalog).expect("the catalog reads");
    let mut expected = BTreeSet::new();
    for position in ["127877125", "127877126"] {
        for allele in ["CG", "C", "G"] {
            expected.insert(pair("8", position, allele, "GCST002273"));
        }
    }
    assert_eq!(read.reference_only, expected);
    assert!(read.candidate_only.is_empty());
    let side = side_of(read);
    let planted_row = |location: &str, allele: &str| {
        planted_pair(location, allele, GWAS_SKIPPED_ROW_FIELDS, "plant05")
    };
    let one = count_map(&[(CROSS_PRODUCT, 1)]);

    let (reference, candidate) = planted_row("8:127877125-127877126", "CG");
    let (report, documented) = compare_gwas(&reference, &candidate, &side, &[CROSS_PRODUCT]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    assert_eq!(masked(&report, "GWAS_study"), one);
    assert_eq!(
        verdict(&compare_gwas(&reference, &candidate, &side, &[]).0),
        "red"
    );
    for (location, allele) in [
        ("8:127877125", "CG"),
        ("8:127877126-127877127", "CG"),
        ("8:127877125", "G"),
    ] {
        let (reference, candidate) = planted_row(location, allele);
        let (report, documented) = compare_gwas(&reference, &candidate, &side, &[CROSS_PRODUCT]);
        assert_eq!(verdict(&report), "green", "{location} {allele}");
        assert_eq!(applied(&documented), one, "{location} {allele}");
    }
    let swapped = planted_row("8:127877125-127877126", "CG");
    for (reference, candidate) in [
        planted_row("8:127877125", "T"),
        planted_row("8:127877127", "G"),
        planted_row("8:127877124-127877125", "CG"),
        planted_row("8:127877125-127877127", "CG"),
        planted_pair(
            "8:127877125-127877126",
            "CG",
            &GWAS_SKIPPED_ROW_FIELDS.replace("GCST002273", "GCST1"),
            "plant05",
        ),
        planted_pair(
            "8:127877125-127877126",
            "CG",
            GWAS_SKIPPED_ROW_FIELDS,
            "8_127877125_A/ACG",
        ),
        (swapped.1, swapped.0),
    ] {
        let (report, documented) = compare_gwas(&reference, &candidate, &side, &[CROSS_PRODUCT]);
        assert_eq!(verdict(&report), "red", "{reference}");
        assert_eq!(applied(&documented), count_map(&[]), "{reference}");
    }
    assert_eq!(
        verdict(&compare_gwas(&reference, &candidate, &Side::default(), &[CROSS_PRODUCT]).0),
        "red"
    );
    // As many SNPs as locations: the candidate pairs each with its own, the reference with
    // both.
    let even = gwas_pairs(
        &directory,
        &[("8;8", "127877125;127877126", "rs1-CG; rs2-C", "GCST1")],
    );
    assert_eq!(
        even.reference_only,
        pairs_of(&[
            ("8", "127877125", "C", "GCST1"),
            ("8", "127877126", "CG", "GCST1")
        ])
    );
    let catalog = planted(
        &directory,
        "gwas_catalog.tsv",
        &format!("{GWAS_CATALOG_HEADER}\n{lists}\trs35252396; rs6470588\tGCST002273\n"),
    );
    let two_snps = GwasCatalogPairs::from_catalog(&catalog).expect("the catalog reads");
    assert!(two_snps.reference_only.is_empty() && two_snps.candidate_only.is_empty());
    // GWAS.pm:261 reads a risk-allele cell holding ? as empty, so no SNP of the row is
    // written.
    for alleles in ["rs1-CG; rs2-?", "rs1-CG; rs2-C; rs3-?", "rs1-CG x rs2-?"] {
        let unknown = gwas_pairs(
            &directory,
            &[("8;8", "127877125;127877126", alleles, "GCST1")],
        );
        assert!(unknown.reference_only.is_empty(), "{alleles}");
        assert!(unknown.candidate_only.is_empty(), "{alleles}");
    }
}

/// The `GWAS_` fields of a row at the site 2:5000402: the accession, the beta coefficient,
/// the p-value, the publication and the study.
fn site_fields(accession: &str, beta: &str, p_value: &str, pmid: &str, study: &str) -> String {
    format!("GWAS_accessions={accession};GWAS_associated_gene=;GWAS_beta_coef={beta};GWAS_odds_ratio=;GWAS_p_value={p_value};GWAS_pmid={pmid};GWAS_risk_allele=A;GWAS_study={study}")
}

fn haplotype_row() -> String {
    site_fields(
        "http://purl.obolibrary.org/obo/MONDO_9990001",
        "",
        "3E-18",
        "99000001",
        "GCST990001",
    )
}

fn own_row() -> String {
    site_fields(
        "http://www.ebi.ac.uk/efo/EFO_9990002",
        "0.0125 unit increase",
        "4E-8",
        "99000002",
        "GCST990002",
    )
}

/// The catalog rows that name 2:5000402: the last location of a haplotype row of study
/// GCST990001 and the one location of a single-SNP row of study GCST990002, whose risk
/// allele another SNP of the haplotype row carries.
fn own_row_side(directory: &Path) -> Side {
    gwas_side(
        directory,
        &[
            (
                "2;2",
                "5000149;5000402",
                "rs9990001-A; rs9990002-C",
                "GCST990001",
            ),
            ("2", "5000402", "rs9990002-A", "GCST990002"),
        ],
    )
}

/// The reference's and the candidate's text of the one row `gwas_own_row`.
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
fn the_cross_product_class_accepts_a_field_the_sites_own_catalog_row_fills() {
    // At 2:5000402 Ensembl VEP reports the haplotype row, which has no beta coefficient,
    // and the candidate the site's own single-SNP row, which has one: the GWAS_beta_coef
    // the candidate alone fills is the class's with the fields both fill differently, the
    // row counted once.
    let directory = scratch("gwas-own-row");
    let side = own_row_side(&directory);
    assert_eq!(
        pairs(&side).reference_only,
        pairs_of(&[
            ("2", "5000402", "A", "GCST990001"),
            ("2", "5000149", "C", "GCST990001")
        ])
    );
    let (reference, candidate) = own_row_texts(&haplotype_row(), &own_row(), "2:5000402", "A");
    let (report, documented) = compare_gwas(&reference, &candidate, &side, &[CROSS_PRODUCT]);
    let one = count_map(&[(CROSS_PRODUCT, 1)]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    assert_eq!(masked(&report, "GWAS_beta_coef"), one);
    assert_eq!(masked(&report, "GWAS_study"), one);
    assert_eq!(
        verdict(&compare_gwas(&reference, &candidate, &side, &[]).0),
        "red"
    );
}

#[test]
fn the_cross_product_class_refuses_a_candidate_only_field_where_the_site_has_no_pair() {
    // The same field stays a difference at a site or under an allele the catalog does not
    // pair, on a reference row annotated under another study or not annotated, where the
    // candidate reports the reference's study or none, and for a field of another plugin.
    let directory = scratch("gwas-own-row-refused");
    let side = own_row_side(&directory);
    let (haplotype, own) = (haplotype_row(), own_row());
    let one = count_map(&[(CROSS_PRODUCT, 1)]);
    let none = count_map(&[]);
    let same_study = own_row_texts(
        &haplotype,
        &own.replace("GCST990002", "GCST990001"),
        "2:5000402",
        "A",
    );
    let off_site = own_row_texts(&haplotype, &own, "2:5000403", "A");
    for ((reference, candidate), counted) in [
        (off_site.clone(), &none),
        (own_row_texts(&haplotype, &own, "2:5000402", "G"), &none),
        (
            own_row_texts(
                &haplotype.replace("GCST990001", "GCST1"),
                &own,
                "2:5000402",
                "A",
            ),
            &none,
        ),
        (own_row_texts("", &own, "2:5000402", "A"), &none),
        (same_study.clone(), &one),
        (
            own_row_texts(
                &haplotype,
                "GWAS_beta_coef=0.0125 unit increase",
                "2:5000402",
                "A",
            ),
            &one,
        ),
        (
            own_row_texts(&haplotype, &format!("{own};LoFtool=0.5"), "2:5000402", "A"),
            &one,
        ),
    ] {
        let (report, documented) = compare_gwas(&reference, &candidate, &side, &[CROSS_PRODUCT]);
        assert_eq!(verdict(&report), "red", "{candidate}");
        assert_eq!(&applied(&documented), counted, "{candidate}");
    }
    let (report, _) = compare_gwas(&same_study.0, &same_study.1, &side, &[CROSS_PRODUCT]);
    assert_eq!(differing(&report, "GWAS_beta_coef"), 1);
    let (report, _) = compare_gwas(&off_site.0, &off_site.1, &side, &[CROSS_PRODUCT]);
    assert_eq!(differing(&report, "GWAS_beta_coef"), 1);
}

#[test]
fn the_cross_product_class_leaves_a_planted_wrong_or_missing_value_a_difference() {
    // At a site the catalog does not pair, where both engines report the site's own row, a
    // GWAS_beta_coef the candidate writes wrong and one it lacks both stay differences, the
    // class applied to nothing.
    let directory = scratch("gwas-planted-value");
    let side = own_row_side(&directory);
    let own = own_row();
    let unpaired = default_text(&[gwas_row("gwas_own_row", "2:5000500", "A", &own)]);
    for planted in [
        own.replace("0.0125", "0.5"),
        own.replace("GWAS_beta_coef=0.0125 unit increase", "GWAS_beta_coef="),
    ] {
        let candidate = default_text(&[gwas_row("gwas_own_row", "2:5000500", "A", &planted)]);
        let (report, documented) = compare_gwas(&unpaired, &candidate, &side, &[CROSS_PRODUCT]);
        assert_eq!(verdict(&report), "red", "{planted}");
        assert_eq!(applied(&documented), count_map(&[]), "{planted}");
        assert_eq!(differing(&report, "GWAS_beta_coef"), 1);
    }
}

#[test]
fn the_cross_product_class_accepts_the_empty_fields_of_the_annotated_row() {
    // The plugin writes a field the catalog row leaves empty as a defined empty string
    // (`GWAS_beta_coef=`, an empty tab cell, `""` in JSON), which a row the candidate never
    // annotated has no key for: beside the row's filled `GWAS_` fields the pair is the
    // class's in every format; a field the candidate alone fills, an empty field on a row
    // the reference fills nowhere, and an empty key the candidate alone writes stay
    // differences.
    let directory = scratch("gwas-empty-fields");
    let side = gwas_side(&directory, &[("21;21", "100;150", "rs1-A; rs2-G", "GCST1")]);
    let one = count_map(&[(CROSS_PRODUCT, 1)]);
    let filled = "GWAS_accessions=http://purl.obolibrary.org/obo/MONDO_0004992;GWAS_associated_gene=;GWAS_beta_coef=;GWAS_odds_ratio=0.41;GWAS_p_value=5E-7;GWAS_pmid=29299148;GWAS_risk_allele=G;GWAS_study=GCST1";
    let reference = default_text(&[gwas_row("a", "21:100", "G", filled)]);
    let candidate = default_text(&[gwas_row("a", "21:100", "G", "")]);
    let (report, documented) = compare_gwas(&reference, &candidate, &side, &[CROSS_PRODUCT]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    assert_eq!(masked(&report, "Extra"), one);
    assert_eq!(
        verdict(&compare_gwas(&reference, &candidate, &side, &[]).0),
        "red"
    );
    // The candidate alone fills the empty field: a difference.
    let filled_by_candidate = candidate.replace("STRAND=1\n", "STRAND=1;GWAS_beta_coef=0.5\n");
    assert_eq!(
        verdict(&compare_gwas(&reference, &filled_by_candidate, &side, &[CROSS_PRODUCT]).0),
        "red"
    );
    // The empty string is the reference's: a key the candidate alone writes empty on an
    // annotated row stays a difference.
    let annotated_both = default_text(&[gwas_row("a", "21:100", "G", "GWAS_study=GCST1")]);
    let empty_key = annotated_both.replace("GWAS_study=GCST1", "GWAS_study=GCST1;GWAS_beta_coef=");
    assert_eq!(
        verdict(&compare_gwas(&annotated_both, &empty_key, &side, &[CROSS_PRODUCT]).0),
        "red"
    );
    // Every GWAS_ field of the reference row empty: nothing annotated the row, so no class.
    let empty_only = default_text(&[gwas_row(
        "a",
        "21:100",
        "G",
        "GWAS_beta_coef=;GWAS_odds_ratio=",
    )]);
    assert_eq!(
        verdict(&compare_gwas(&empty_only, &candidate, &side, &[CROSS_PRODUCT]).0),
        "red"
    );
    // Tab: an empty cell against the candidate's `-`.
    let tab_header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIMPACT\tSTRAND\tGWAS_beta_coef\tGWAS_odds_ratio\tGWAS_study\n";
    let tab_reference = format!(
        "{tab_header}a\t21:100\tG\tENSG1\tENST1\tTranscript\tintron_variant\tMODIFIER\t1\t\t0.41\tGCST1\n"
    );
    let tab_candidate = format!(
        "{tab_header}a\t21:100\tG\tENSG1\tENST1\tTranscript\tintron_variant\tMODIFIER\t1\t-\t-\t-\n"
    );
    let (report, documented) = compare_under(
        &parse_tab(&tab_reference).expect("the reference parses"),
        &parse_tab(&tab_candidate).expect("the candidate parses"),
        &side,
        &[CROSS_PRODUCT],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    // JSON: `"gwas_beta_coef": ""` against no key; the record object states the site.
    let inputs = lines_of(&["21\t100\ta\tA\tG\t.\t.\t."]);
    let entry = r#""variant_allele": "G", "transcript_id": "ENST1", "consequence_terms": ["intron_variant"], "impact": "MODIFIER", "strand": 1"#;
    let record = format!(
        r#""input": {}, "id": "a", "seq_region_name": "21", "start": 100, "most_severe_consequence": "intron_variant""#,
        quoted(&inputs[0])
    );
    let reference_json = format!(
        "{{{record}, \"transcript_consequences\": [{{{entry}, \"gwas_beta_coef\": \"\", \"gwas_odds_ratio\": \"0.41\", \"gwas_study\": \"GCST1\"}}]}}\n"
    );
    let candidate_json = format!("{{{record}, \"transcript_consequences\": [{{{entry}}}]}}\n");
    let (report, documented) = compare_under(
        &parse_json(&reference_json, &inputs).expect("the reference parses"),
        &parse_json(&candidate_json, &inputs).expect("the candidate parses"),
        &side,
        &[CROSS_PRODUCT],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
}

const INTERACTION_FIELDS: &str =
    "GWAS_p_value=3E-7;GWAS_pmid=29691896;GWAS_risk_allele=C;GWAS_study=GCST006494";

/// The catalog's interaction row of rs4818032 (chromosome 21) and rs10830963 (chromosome
/// 11), the chromosome 21 SNP written first.
fn interaction_side(directory: &Path) -> Side {
    gwas_side(
        directory,
        &[(
            "21 x 11",
            "39375657 x 92975544",
            "rs4818032-C x rs10830963-G",
            "GCST006494",
        )],
    )
}

#[test]
fn the_interaction_class_accepts_the_first_snps_own_allele() {
    // 21 39375657 T>C against the interaction row rs4818032-C x rs10830963-G: the reference
    // keeps the space after the first risk allele and annotates nothing, the candidate
    // writes the row's fields; the VCF and JSON formats state the site on their record
    // entry.
    let directory = scratch("gwas-interaction");
    let side = interaction_side(&directory);
    assert_eq!(
        pairs(&side).candidate_only,
        pairs_of(&[("21", "39375657", "C", "GCST006494")])
    );
    assert_eq!(
        pairs(&side).reference_only,
        pairs_of(&[("21", "39375657", "G", "GCST006494")])
    );
    let one = count_map(&[(INTERACTION, 1)]);
    let reference = default_text(&[gwas_row("rs4818032", "21:39375657", "C", "")]);
    let candidate = default_text(&[gwas_row(
        "rs4818032",
        "21:39375657",
        "C",
        INTERACTION_FIELDS,
    )]);
    let (report, documented) = compare_gwas(&reference, &candidate, &side, &[INTERACTION]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    assert_eq!(masked(&report, "GWAS_study"), one);
    assert_eq!(
        verdict(&compare_gwas(&reference, &candidate, &side, &[]).0),
        "red"
    );
    // The candidate's row carries the catalog row's empty cell as an empty field: the
    // class's with the filled ones, and a difference on a row the candidate fills nowhere.
    let with_empty = default_text(&[gwas_row(
        "rs4818032",
        "21:39375657",
        "C",
        &format!("GWAS_odds_ratio=;{INTERACTION_FIELDS}"),
    )]);
    let (report, documented) = compare_gwas(&reference, &with_empty, &side, &[INTERACTION]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    let empty_alone = default_text(&[gwas_row(
        "rs4818032",
        "21:39375657",
        "C",
        "GWAS_odds_ratio=",
    )]);
    let (report, documented) = compare_gwas(&reference, &empty_alone, &side, &[INTERACTION]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), count_map(&[]));
    let inputs = lines_of(&["21\t39375657\trs4818032\tT\tC\t.\t.\t."]);
    let entry = r#""variant_allele": "C", "transcript_id": "ENST1", "consequence_terms": ["intron_variant"]"#;
    let record = format!(
        r#""input": {}, "id": "rs4818032", "seq_region_name": "21", "start": 39375657, "most_severe_consequence": "intron_variant""#,
        quoted(&inputs[0])
    );
    let scored = format!(
        "{{{record}, \"transcript_consequences\": [{{{entry}, \"gwas_risk_allele\": \"C\", \"gwas_study\": \"GCST006494\"}}]}}\n"
    );
    let plain = format!("{{{record}, \"transcript_consequences\": [{{{entry}}}]}}\n");
    let (report, documented) = compare_under(
        &parse_json(&plain, &inputs).expect("the reference parses"),
        &parse_json(&scored, &inputs).expect("the candidate parses"),
        &side,
        &[INTERACTION],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    let vcf = format!(
        "{STUDY_CSQ_HEADER}{}CSQ=C|intron_variant|ENST1|Transcript|\n",
        before_info(&inputs[0])
    );
    let (report, documented) = compare_under(
        &parse_vcf(&vcf, &inputs).expect("the reference parses"),
        &parse_vcf(
            &vcf.replace("Transcript|\n", "Transcript|GCST006494\n"),
            &inputs,
        )
        .expect("the candidate parses"),
        &side,
        &[INTERACTION],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
}

#[test]
fn the_interaction_class_refuses_any_other_field_the_candidate_alone_fills() {
    // A GWAS_ field the candidate alone fills stays a difference without the catalog, at
    // the row's last SNP, under an allele that is not the SNP's risk allele, under another
    // study, and where the reference annotates the row too; a value planted missing or
    // wrong on an ordinary row is counted under its field.
    let directory = scratch("gwas-interaction-refused");
    let side = interaction_side(&directory);
    let none = count_map(&[]);
    let reference = default_text(&[gwas_row("rs4818032", "21:39375657", "C", "")]);
    let candidate = default_text(&[gwas_row(
        "rs4818032",
        "21:39375657",
        "C",
        INTERACTION_FIELDS,
    )]);
    let (report, documented) =
        compare_gwas(&reference, &candidate, &Side::default(), &[INTERACTION]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none);
    for (old, new) in [
        ("21:39375657", "11:92975544"),
        ("\tC\t", "\tA\t"),
        ("GCST006494", "GCST000001"),
    ] {
        let (report, documented) = compare_gwas(
            &reference.replace(old, new),
            &candidate.replace(old, new),
            &side,
            &[INTERACTION],
        );
        assert_eq!(verdict(&report), "red", "{old} {new}");
        assert_eq!(applied(&documented), none, "{old} {new}");
    }
    let annotated = default_text(&[gwas_row(
        "rs4818032",
        "21:39375657",
        "C",
        "GWAS_study=GCST000028",
    )]);
    let (report, documented) = compare_gwas(&annotated, &candidate, &side, &[INTERACTION]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none);
    let ordinary = default_text(&[gwas_row("a", "21:500", "C", INTERACTION_FIELDS)]);
    for planted in [
        default_text(&[gwas_row("a", "21:500", "C", "")]),
        ordinary.replace("GCST006494", "GCST000001"),
    ] {
        let (report, documented) = compare_gwas(&ordinary, &planted, &side, &[INTERACTION]);
        assert_eq!(verdict(&report), "red");
        assert_eq!(applied(&documented), none);
        assert_eq!(differing(&report, "GWAS_study"), 1);
    }
}

/// Whether the class `klass` accepts the field `name` of `e`, the reference's value
/// `expected_value` against no value of the candidate's, beside the rows given.
fn accepts(
    klass: &str,
    e: &Entry,
    name: &str,
    expected_value: &str,
    record_rows: Option<&[Entry]>,
    candidate_rows: Option<&[Entry]>,
) -> bool {
    let record_rows = record_rows.map(refs);
    let candidate_rows = candidate_rows.map(refs);
    let side = Side::default();
    let context = RuleContext {
        record_rows: record_rows.as_deref(),
        candidate_rows: candidate_rows.as_deref(),
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    field_class_accepts(klass, e, name, expected_value, "-", None, &context)
        .expect("the class is asked")
}

#[test]
fn the_classes_of_a_field_the_candidate_alone_fills_accept_its_value_alone() {
    // alphamissense_two_digit_transcript_version, cadd_structural_allele_match (on a
    // structural allele's row alone), splicevault_json_later_rows (the JSON object alone)
    // and loeuf_zero_minimum (a 0 alone) accept a field the candidate fills where the
    // reference has none; the reverse direction and any other value stay differences.
    let reference = default_text(&[
        default_row(
            "a",
            "21:100",
            "G",
            "ENST1",
            "missense_variant",
            "IMPACT=MODERATE;STRAND=1",
        ),
        default_row(
            "b",
            "21:200-900",
            "deletion",
            "ENST2",
            "coding_sequence_variant",
            "IMPACT=MODIFIER;STRAND=1",
        ),
        default_row(
            "c",
            "21:300",
            "T",
            "ENST3",
            "intron_variant",
            "IMPACT=MODIFIER;STRAND=1",
        ),
        default_row(
            "d",
            "21:400",
            "C",
            "ENST4",
            "intron_variant",
            "IMPACT=MODIFIER;STRAND=1",
        ),
    ]);
    let candidate = default_text(&[
        default_row(
            "a",
            "21:100",
            "G",
            "ENST1",
            "missense_variant",
            "IMPACT=MODERATE;STRAND=1;am_class=likely_pathogenic;am_pathogenicity=0.9",
        ),
        default_row(
            "b",
            "21:200-900",
            "deletion",
            "ENST2",
            "coding_sequence_variant",
            "IMPACT=MODIFIER;STRAND=1;CADD_PHRED=12.5;CADD_RAW=0.9",
        ),
        default_row(
            "c",
            "21:300",
            "T",
            "ENST3",
            "intron_variant",
            "IMPACT=MODIFIER;STRAND=1;CADD_PHRED=3.1",
        ),
        default_row(
            "d",
            "21:400",
            "C",
            "ENST4",
            "intron_variant",
            "IMPACT=MODIFIER;STRAND=1;LOEUF=0",
        ),
    ]);
    let classes = [
        "alphamissense_two_digit_transcript_version",
        "cadd_structural_allele_match",
        "loeuf_zero_minimum",
    ];
    let (report, documented) = compare_default(&reference, &candidate, &classes);
    // The sequence allele's CADD row is not the CADD class's.
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "CADD_PHRED"), 1);
    assert_eq!(
        applied(&documented),
        count_map(&[
            ("alphamissense_two_digit_transcript_version", 1),
            ("cadd_structural_allele_match", 1),
            ("loeuf_zero_minimum", 1)
        ])
    );
    let without_c = candidate.replace(";CADD_PHRED=3.1", "");
    assert_eq!(verdict_of(&reference, &without_c, &classes), "green");
    assert_eq!(verdict_of(&without_c, &reference, &classes), "red");
    assert_eq!(
        verdict_of(
            &reference,
            &without_c.replace("LOEUF=0", "LOEUF=0.5"),
            &classes
        ),
        "red"
    );
    // The JSON splicevault object the candidate alone carries, on a transcript entry.
    let inputs = lines_of(&["21\t100\trs1\tA\tG\t.\t.\t."]);
    let record = format!(
        r#""input": {}, "id": "rs1", "most_severe_consequence": "intron_variant""#,
        quoted(&inputs[0])
    );
    let entry = r#""variant_allele": "G", "transcript_id": "ENST1", "consequence_terms": ["intron_variant"]"#;
    let reference_json = format!("{{{record}, \"transcript_consequences\": [{{{entry}}}]}}\n");
    let candidate_json = format!(
        "{{{record}, \"transcript_consequences\": [{{{entry}, \"splicevault\": {{\"1\": {{\"site_type\": \"donor\"}}}}}}]}}\n"
    );
    let without = parse_json(&reference_json, &inputs).expect("the object parses");
    let with_object = parse_json(&candidate_json, &inputs).expect("the object parses");
    let no_side = Side::default();
    let later_rows = "splicevault_json_later_rows";
    let (report, documented) = compare_under(&without, &with_object, &no_side, &[later_rows]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(later_rows, 1)]));
    let (report, _) = compare_under(&with_object, &without, &no_side, &[later_rows]);
    assert_eq!(verdict(&report), "red");
    // The text formats' SpliceVault_ fields are the sibling-allele class's on a record of
    // two alleles: the reference writes the first allele's table row on the G row as well,
    // where the candidate leaves it empty (the table has no row for the allele) or writes
    // the allele's own row; a field the candidate alone fills stays a difference.
    let klass = "splicevault_sibling_allele_rows";
    let first = "SpliceVault_site_type=donor;SpliceVault_site_pos=5";
    let splicevault_rows = |g_fields: &str, name: &str, alleles: &str| {
        let rows: Vec<String> = [("C", first), ("G", g_fields)]
            .iter()
            .filter(|(allele, _)| alleles.contains(allele))
            .map(|(allele, fields)| {
                let extra = if fields.is_empty() {
                    "IMPACT=HIGH;STRAND=1".to_string()
                } else {
                    format!("IMPACT=HIGH;STRAND=1;{fields}")
                };
                default_row(
                    name,
                    "21:100",
                    allele,
                    "ENST1",
                    "splice_donor_variant",
                    &extra,
                )
            })
            .collect();
        default_text(&rows)
    };
    let two_alleles = |g_fields: &str| splicevault_rows(g_fields, "21_100_T/C/G", "CG");
    let own_fields = "SpliceVault_site_type=acceptor;SpliceVault_site_pos=8";
    let reference = two_alleles(first);
    let own_row = two_alleles(own_fields);
    let one = count_map(&[(klass, 1)]);
    for candidate in [two_alleles(""), own_row.clone()] {
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "green");
        assert_eq!(applied(&documented), one);
    }
    assert_eq!(verdict_of(&two_alleles(""), &reference, &[klass]), "red");
    assert_eq!(verdict_of(&reference, &own_row, &[]), "red");
    // The same pair of values on a record of one allele, on rows under an identifier, and
    // against a sibling row carrying another value is no sibling's row: a difference.
    let one_allele = |g_fields: &str| splicevault_rows(g_fields, "21_100_T/G", "G");
    for (expected_rows, candidate_rows) in [
        (one_allele(first), one_allele("")),
        (one_allele(first), one_allele(own_fields)),
        (
            splicevault_rows(first, "sibling", "CG"),
            splicevault_rows("", "sibling", "CG"),
        ),
        (
            two_alleles("SpliceVault_site_type=acceptor;SpliceVault_site_pos=6"),
            own_row.clone(),
        ),
    ] {
        let (report, documented) = compare_default(&expected_rows, &candidate_rows, &[klass]);
        assert_eq!(verdict(&report), "red", "{candidate_rows}");
        assert_eq!(applied(&documented), count_map(&[]), "{candidate_rows}");
    }
    // The VCF and JSON entries name their record: the class reads the sibling's entry
    // there.
    let inputs = lines_of(&["21\t100\tsibling\tT\tC,G\t.\t.\t."]);
    let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|Feature_type|SpliceVault_site_type\">\n";
    let vcf = format!(
        "{header}{}CSQ=C|splice_donor_variant|ENST1|Transcript|donor,G|splice_donor_variant|ENST1|Transcript|donor\n",
        before_info(&inputs[0])
    );
    let (report, documented) = compare_under(
        &parse_vcf(&vcf, &inputs).expect("the reference parses"),
        &parse_vcf(
            &vcf.replace("Transcript|donor\n", "Transcript|acceptor\n"),
            &inputs,
        )
        .expect("the candidate parses"),
        &no_side,
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one);
    // A direct call without the record's rows refuses.
    let reference_rows = parse_default(&reference).expect("the reference parses");
    let entry = &reference_rows[1];
    let field = "SpliceVault_site_type";
    assert!(!accepts(klass, entry, field, "donor", None, None));
    let kept = parse_default(&two_alleles("")).expect("the candidate parses");
    assert!(accepts(
        klass,
        entry,
        field,
        "donor",
        Some(&reference_rows),
        Some(&kept)
    ));
    assert!(!accepts(
        klass,
        entry,
        field,
        "donor",
        Some(&reference_rows),
        None
    ));
    // The rows of another record at the Location are no sibling's.
    let other_record =
        parse_default(&splicevault_rows(first, "21_100_A/C/G", "CG")).expect("the rows parse");
    assert!(!accepts(
        klass,
        entry,
        field,
        "donor",
        Some(&other_record),
        Some(&kept)
    ));
}

/// What a run says that ends on an input no reading of the comparison expects.
fn uncaught<T>(run: impl FnOnce() -> T) -> String {
    let payload = std::panic::catch_unwind(std::panic::AssertUnwindSafe(run))
        .err()
        .expect("the run ends");
    payload
        .downcast_ref::<Uncaught>()
        .expect("an ending of the run's own")
        .0
        .clone()
}

#[test]
fn a_catalog_row_is_cut_at_every_separator_and_an_allele_at_its_second_hyphen() {
    // A risk allele is the text between the first and the second `-` of its token, and a
    // list is cut at every `;`, `,` and small `x`, one inside a name included.
    let mut read = GwasCatalogPairs::new();
    read.add_row("21;21", "100;150", "rs1-G-T; rs2-A", "GCST1", "");
    assert_eq!(
        read.reference_only,
        pairs_of(&[("21", "100", "A", "GCST1"), ("21", "150", "G", "GCST1")])
    );
    let mut cut = GwasCatalogPairs::new();
    cut.add_row("21;21", "100;150", "rs1-G; rs2-A", "GCST1", "");
    cut.add_row("2x1", "100;150", "rs1-G; rs2-A", "GCST2", "");
    cut.add_row("X;X", "5;6", "rs3-C; rs4-T", "GCST3", "");
    assert_eq!(
        cut.reference_only,
        pairs_of(&[
            ("1", "150", "G", "GCST2"),
            ("2", "100", "A", "GCST2"),
            ("21", "100", "A", "GCST1"),
            ("21", "150", "G", "GCST1"),
            ("X", "5", "T", "GCST3"),
            ("X", "6", "C", "GCST3"),
        ])
    );
    assert!(cut.candidate_only.is_empty());
}

fn located(location: &str, allele: &str) -> Entry {
    entry(
        &format!("{location}|{allele}|ENST1|Transcript"),
        None,
        Some(location),
        allele,
        "ENST1",
        "intron_variant",
        &[],
    )
}

#[test]
fn a_site_is_read_from_positions_of_any_length_and_as_the_interpreter_converts_them() {
    // A text row's position is its Location's text as written; a span is measured on whole
    // numbers of any length; and a text of digits the interpreter does not convert ends the
    // run where the span is measured, in the interpreter's words.
    let site =
        |chromosome: &str, position: &str| Some((chromosome.to_string(), position.to_string()));
    assert_eq!(
        gwas_site(&located("chr21:0100", "A"), None),
        site("21", "0100")
    );
    let long = "123456789012345678901234567890";
    assert_eq!(
        gwas_site(
            &located(&format!("21:{long}-123456789012345678901234567891"), "CG"),
            None
        ),
        site("21", long)
    );
    assert_eq!(
        gwas_site(
            &located(&format!("21:{long}-123456789012345678901234567892"), "CG"),
            None
        ),
        None
    );
    assert_eq!(gwas_site(&located("21:1²", "A"), None), site("21", "1²"));
    assert_eq!(
        uncaught(|| gwas_site(&located("21:1²-3", "AC"), None)),
        "invalid literal for int() with base 10: '1²'"
    );
    // A VCF row's start is the line's POS moved by the base its alleles share, the POS a
    // number of any length.
    let record = entry(
        "3|record",
        Some(3),
        None,
        "",
        "",
        "",
        &[
            ("CHROM", "chr21"),
            ("POS", "99999999999999999999"),
            ("REF", "CA"),
            ("ALT", "CG,C"),
        ],
    );
    let row = entry(
        "3|G|ENST1|Transcript",
        Some(3),
        None,
        "G",
        "ENST1",
        "intron_variant",
        &[],
    );
    assert_eq!(
        gwas_site(&row, Some(&[&record, &row])),
        site("21", "100000000000000000000")
    );
}
