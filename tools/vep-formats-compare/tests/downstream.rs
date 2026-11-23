// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The peptide paths of a transcript under another codon table, the HGNC identifier copy,
//! the genetic codes and the two classes of the Downstream plugin's start rule, on planted
//! rows and caches.

mod downstream_common;

use downstream_common::*;
use std::path::Path;
use vep_formats_compare::classes::tables;
use vep_formats_compare::classes::HGNC_ID_COPY;
use vep_formats_compare::classes::ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE as ORGANELLE_PEPTIDE;
use vep_formats_compare::classes::START_COOCCURRENCE_SWAP;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::fields::OrganelleCodonTables;
use vep_formats_compare::genetic_codes::genetic_code;
use vep_formats_compare::genetic_codes::held_tables;
use vep_formats_compare::genetic_codes::CodonTable;
use vep_formats_compare::genetic_codes::Reading;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::downstream::codon_residues;
use vep_formats_compare::rules::downstream::codon_table_id;
use vep_formats_compare::rules::downstream::downstream_class;
use vep_formats_compare::rules::downstream::translate_until_stop;
use vep_formats_compare::rules::downstream::DownstreamSequences;
use vep_formats_compare::rules::downstream::DownstreamTranscript;
use vep_formats_compare::rules::downstream::BLOSUM62_MATRIX;
use vep_formats_compare::rules::field_class_accepts;
use vep_formats_compare::rules::RuleContext;
use vep_formats_compare::rules::Side;

fn default_entries(text: &str) -> Vec<Entry> {
    parse_default(text).expect("the text parses")
}

/// A JSON record whose consequence entry is a registry pair may summarise the pair's term
/// (start_lost in the reference, start_retained_variant in the candidate): the record's
/// most_severe_consequence is set aside as for a documented key.
#[test]
fn a_registry_pair_leaves_the_json_records_most_severe_consequence_uncompared() {
    let inputs = ["21\t100\trs1\tA\tG\t.\t.\t.".to_string()];
    let record = |terms: &str, severe: &str, impact: &str| {
        dumped_value(&object(vec![
            ("input", text(&inputs[0])),
            ("id", text("rs1")),
            ("seq_region_name", text("21")),
            ("start", whole(100)),
            ("end", whole(100)),
            ("allele_string", text("A/G")),
            ("strand", whole(1)),
            ("most_severe_consequence", text(severe)),
            (
                "transcript_consequences",
                list(vec![object(vec![
                    ("variant_allele", text("G")),
                    ("transcript_id", text("ENST1")),
                    ("gene_id", text("ENSG1")),
                    ("impact", text(impact)),
                    (
                        "consequence_terms",
                        list(terms.split(',').map(text).collect()),
                    ),
                    ("strand", whole(1)),
                ])]),
            ),
        ])) + "\n"
    };
    let reference = record("start_lost,start_retained_variant", "start_lost", "HIGH");
    let candidate = record("start_retained_variant", "start_retained_variant", "LOW");
    let expected = parse_json(&reference, &inputs).expect("the reference parses");
    let actual = parse_json(&candidate, &inputs).expect("the candidate parses");
    let side = Side::default();
    let (report, documented) = compared(&expected, &actual, &side, &[START_COOCCURRENCE_SWAP]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(
        counted(&documented.applied),
        count_map(&[(START_COOCCURRENCE_SWAP, 1)])
    );
    let mut none = Documented::new();
    let report = compared_under(&expected, &actual, &mut none, &side, None);
    assert_eq!(verdict_of(&report), "red");
}

/// Two records sharing a Location, Allele and Feature (an insertion and a substitution of
/// the same bases) each keep their own rows: a documented field value applies to the row
/// whose consequence set the item states; a manifest without the set applies to every row
/// of the key.
#[test]
fn a_field_divergence_applies_to_the_row_of_its_consequence_set_alone() {
    let rows = [
        (
            "ins",
            "frameshift_variant",
            "IMPACT=HIGH;STRAND=1;HGVSp=Q0080.1:p.Met1?",
        ),
        (
            "sub",
            "missense_variant",
            "IMPACT=MODERATE;STRAND=1;HGVSp=Q0080.1:p.Ile2Thr",
        ),
    ]
    .map(|(name, consequence, extra)| {
        default_row(name, "Mito:27667-27668", "A", "Q0080", consequence, extra)
    });
    let reference = default_text(&rows);
    let candidate = reference.replace("p.Met1?", "p.Pro2_?1");
    let item = |consequence_set: Option<&str>| {
        let mut members = vec![
            ("location", text("Mito:27667-27668")),
            ("allele", text("A")),
            ("feature", text("Q0080")),
            ("field", text("HGVSp")),
            ("vep", text("Q0080.1:p.Met1?")),
            ("vep_rs", text("Q0080.1:p.Pro2_?1")),
            ("expected_divergence", text(ORGANELLE_PEPTIDE)),
        ];
        members.extend(consequence_set.map(|terms| ("consequence_set", text(terms))));
        object(vec![("field_divergences", list(vec![object(members)]))])
    };
    let masks = [ORGANELLE_PEPTIDE.to_string()].into_iter().collect();
    let documented = |consequence_set: Option<&str>| {
        Documented::from_manifest(Some(&item(consequence_set)), &masks, None, None, None)
            .expect("the manifest reads")
    };
    let side = Side::default();
    let verdict = |reference: &str, candidate: &str, documented: &mut Documented| {
        let report = compared_under(
            &default_entries(reference),
            &default_entries(candidate),
            documented,
            &side,
            None,
        );
        verdict_of(&report).to_string()
    };
    let mut keyed = documented(Some("frameshift_variant"));
    assert_eq!(verdict(&reference, &candidate, &mut keyed), "green");
    // The item names the other row's set: the frameshift row's HGVSp is compared exactly
    // and differs.
    let mut other = documented(Some("missense_variant"));
    assert_eq!(verdict(&reference, &candidate, &mut other), "red");
    // A manifest written without the set applies to every row of the key: the
    // substitution's HGVSp is then read against the documented value and is stale; on a key
    // with one row the set-less item still applies.
    let mut unkeyed = documented(None);
    assert_eq!(verdict(&reference, &candidate, &mut unkeyed), "red");
    let single = default_text(&rows[..1]);
    assert_eq!(
        verdict(
            &single,
            &single.replace("p.Met1?", "p.Pro2_?1"),
            &mut unkeyed
        ),
        "green"
    );
}

/// A candidate JSON cache of six coding transcripts: TXMT on MT under codon table 2, the
/// dump's string, TXN on chromosome 21 under the standard table, a number, TXPT on Pt
/// under table 11 and TXSTRICT on St under table 24, which read every codon as the standard
/// table does, TXPACHY on Pa under genetic code 26, which Ensembl VEP reads under its
/// module's last table, and TXNONE on No under table 7, which no table is held for.
fn organelle_cache(root: &Path) -> std::path::PathBuf {
    for (chromosome, stable_id, table) in [
        ("MT", "TXMT", text("2")),
        ("21", "TXN", whole(1)),
        ("Pt", "TXPT", whole(11)),
        ("St", "TXSTRICT", text("24")),
        ("Pa", "TXPACHY", whole(26)),
        ("No", "TXNONE", whole(7)),
    ] {
        let transcript = object(vec![
            ("stable_id", text(stable_id)),
            ("start", whole(100)),
            ("end", whole(400)),
            ("strand", whole(1)),
            ("biotype", text("protein_coding")),
            (
                "variation_effect_feature_cache",
                object(vec![("codon_table", table)]),
            ),
        ]);
        written(
            root,
            &format!("transcripts/{chromosome}/1-1000000.json"),
            &dumped_value(&list(vec![transcript])),
        );
    }
    root.to_path_buf()
}

fn codon_tables_of(cache: &Path) -> Side {
    Side {
        codon_tables: Some(
            OrganelleCodonTables::new(cache).expect("the cache has a transcripts directory"),
        ),
        ..Side::default()
    }
}

const ORGANELLE_CSQ_HEADER: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|HGVSp\">\n";

/// A frameshift on the transcript under codon table 2 as each engine writes it: an unknown
/// stop distance against a number.
fn organelle_frameshift_texts() -> (String, String) {
    let reference = default_text(&[default_row(
        "r1",
        "MT:104",
        "-",
        "TXMT",
        "frameshift_variant",
        "IMPACT=HIGH;STRAND=1;HGVSc=TXMT:c.5del;HGVSp=PMT.1:p.Ala2GlyfsTer?",
    )]);
    let candidate = reference.replace("fsTer?", "fsTer3");
    (reference, candidate)
}

/// Ensembl VEP translates a frameshift's new frame with the standard code
/// (TranscriptVariationAllele.pm:2333, :2492), so on a transcript the candidate cache files
/// under another codon table the HGVSp differs in its stop distance, its residues or both;
/// the class accepts the pair in every format, counted once per row, and reads the table
/// from the cache: without one it accepts nothing. A lost stop's extension count and a
/// duplication against an insertion are the class's other paths.
#[test]
fn organelle_peptide_paths_accepts_the_bare_translation_paths_on_a_transcript_of_another_codon_table(
) {
    let directory = scratch_of("organelle-accepts");
    let cache = organelle_cache(&directory.join("cache"));
    let tables = codon_tables_of(&cache);
    let no_tables = Side::default();
    let applied_once = count_map(&[(ORGANELLE_PEPTIDE, 1)]);
    let (reference, candidate) = organelle_frameshift_texts();
    let (report, documented) =
        compare_downstream(&reference, &candidate, &tables, &[ORGANELLE_PEPTIDE]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), applied_once);
    assert_eq!(masked(&report, "HGVSp"), applied_once);
    assert_eq!(masked(&report, "Extra"), applied_once);
    let (report, documented) =
        compare_downstream(&reference, &candidate, &no_tables, &[ORGANELLE_PEPTIDE]);
    assert_eq!(verdict_of(&report), "red");
    assert!(counted(&documented.applied).is_empty());
    let (report, _) = compare_downstream(&reference, &candidate, &tables, &[]);
    assert_eq!(verdict_of(&report), "red");
    // Two stop distances, another residue, an immediate stop, a Feature carrying its
    // version, a shift both rows state, and a two-base insertion shifted by its own length,
    // by three bases on the forward strand and by one base on the reverse strand.
    let inserted = |text: &str, strand: &str, shift: i64| {
        text.replace("\t-\tENSG1", "\tCA\tENSG1")
            .replace("STRAND=1", &format!("STRAND={strand}"))
            .replace("c.5del;", &format!("c.5_6dup;HGVS_OFFSET={shift};"))
    };
    let mut pairs = vec![
        (reference.replace("fsTer?", "fsTer9"), candidate.clone()),
        (
            reference.replace("c.5del;", "c.5del;HGVS_OFFSET=4;"),
            candidate.replace("c.5del;", "c.5del;HGVS_OFFSET=4;"),
        ),
    ];
    for (strand, shift) in [("1", 2), ("1", 3), ("-1", 1), ("-1", 4)] {
        pairs.push((
            inserted(&reference, strand, shift),
            inserted(&candidate, strand, shift),
        ));
    }
    pairs.extend([
        (
            reference.clone(),
            reference.replace("Ala2GlyfsTer?", "Ala2SerfsTer3"),
        ),
        (
            reference.clone(),
            reference.replace("Ala2GlyfsTer?", "Ala2Ter"),
        ),
        (
            reference.replace("\tTXMT\t", "\tTXMT.7\t"),
            candidate.replace("\tTXMT\t", "\tTXMT.7\t"),
        ),
    ]);
    for (reference_text, candidate_text) in &pairs {
        let (report, documented) = compare_downstream(
            reference_text,
            candidate_text,
            &tables,
            &[ORGANELLE_PEPTIDE],
        );
        assert_eq!(verdict_of(&report), "green", "{reference_text}");
        assert_eq!(
            counted(&documented.applied),
            applied_once,
            "{reference_text}"
        );
    }
    // The extension count after a lost stop, and a duplication the reference writes as an
    // insertion.
    for (consequence, reference_value, candidate_value) in [
        ("stop_lost", "p.Ter9GlnextTer?", "p.Ter9GlnextTer4"),
        ("inframe_insertion", "p.Met3_Gln4insMet", "p.Met3dup"),
    ] {
        let row = |value: &str| {
            default_text(&[default_row(
                "r2",
                "MT:120-121",
                "T",
                "TXMT",
                consequence,
                &format!("IMPACT=HIGH;STRAND=1;HGVSp=PMT.1:{value}"),
            )])
        };
        let (report, documented) = compare_downstream(
            &row(reference_value),
            &row(candidate_value),
            &tables,
            &[ORGANELLE_PEPTIDE],
        );
        assert_eq!(verdict_of(&report), "green", "{consequence}");
        assert_eq!(counted(&documented.applied), applied_once, "{consequence}");
    }
    // The tab format; the VCF and JSON formats, whose entries carry no Location, read
    // every shard.
    let tab = |distance: &str| {
        parse_tab(&format!(
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tHGVSp\nr1\tMT:104\t-\tG1\tTXMT\tTranscript\tframeshift_variant\tPMT.1:p.Ala2GlyfsTer{distance}\n"
        ))
        .expect("the tab row parses")
    };
    let (report, documented) = compared(&tab("?"), &tab("3"), &tables, &[ORGANELLE_PEPTIDE]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), applied_once);
    let inputs = ["MT\t103\tr1\tCA\tC\t.\t.\t.".to_string()];
    let vcf = |distance: &str| {
        parse_vcf(
            &format!(
                "{ORGANELLE_CSQ_HEADER}MT\t103\tr1\tCA\tC\t.\t.\tCSQ=-|frameshift_variant|Transcript|TXMT|PMT.1:p.Ala2GlyfsTer{distance}\n"
            ),
            &inputs,
        )
        .expect("the VCF line parses")
    };
    let (report, documented) = compared(&vcf("?"), &vcf("3"), &tables, &[ORGANELLE_PEPTIDE]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), applied_once);
    assert_eq!(masked(&report, "HGVSp"), applied_once);
    let json = |distance: &str| {
        let entry = object(vec![
            ("variant_allele", text("-")),
            ("transcript_id", text("TXMT")),
            ("consequence_terms", list(vec![text("frameshift_variant")])),
            ("hgvsp", text(&format!("PMT.1:p.Ala2GlyfsTer{distance}"))),
        ]);
        let record = object(vec![
            ("input", text(&inputs[0])),
            ("id", text("r1")),
            ("most_severe_consequence", text("frameshift_variant")),
            ("transcript_consequences", list(vec![entry])),
        ]);
        parse_json(&(dumped_value(&record) + "\n"), &inputs).expect("the JSON object parses")
    };
    let (report, documented) = compared(&json("?"), &json("3"), &tables, &[ORGANELLE_PEPTIDE]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), applied_once);
}

/// The command line reads a frameshift on a transcript of another codon table: the class is
/// unused without --json-cache and applied with it, and a cache without a transcripts
/// directory is refused.
#[test]
fn organelle_peptide_paths_reads_the_cache_the_command_line_names() {
    let directory = scratch_of("organelle-command-line");
    let cache = organelle_cache(&directory.join("cache"));
    let applied_once = count_map(&[(ORGANELLE_PEPTIDE, 1)]);
    let (reference, candidate) = organelle_frameshift_texts();
    let arguments = with(
        &default_arguments(&directory, &reference, &candidate),
        &["--masks", ORGANELLE_PEPTIDE],
    );
    assert_eq!(run(&arguments).code, 1);
    assert_eq!(
        texts_of(member(&report_of(&directory), &["masks", "unused"])),
        [ORGANELLE_PEPTIDE]
    );
    assert_eq!(
        run(&with(&arguments, &["--json-cache", path_text(&cache)])).code,
        0
    );
    assert_eq!(
        counts_of(member(&report_of(&directory), &["masks", "applied"])),
        applied_once
    );
    let nowhere = directory.join("nowhere");
    let ran = run(&with(&arguments, &["--json-cache", path_text(&nowhere)]));
    assert_eq!(ran.code, 2);
    assert!(
        ran.stderr.ends_with(&format!(
            "--json-cache: {} has no transcripts/ directory\n",
            path_text(&nowhere)
        )),
        "{}",
        ran.stderr
    );
}

/// The nearest shape that stays a difference: the same frameshift pair, an unknown stop
/// distance against a number, on a transcript the cache files under the standard table,
/// where both engines translate alike. So do a transcript under table 11 or table 24, which
/// read every codon as the standard table does, a transcript under genetic code 26, which
/// Ensembl VEP reads under its module's last table and so as the standard code (the class
/// is a difference between two translations of Ensembl VEP's own), a transcript the cache
/// lacks, a substitution's HGVSp on
/// the organelle transcript, a value one side leaves empty and two protein names. On the
/// organelle transcript the class explains the HGVSp alone: a row whose consequence set,
/// positions, codons, HGVSc or HGVS_OFFSET differ too, and a two-base insertion on a
/// reverse-strand transcript shifted by three bases
/// (TranscriptVariationAllele.pm:1337-1339 leave it unrotated), keep their HGVSp a
/// difference.
#[test]
fn organelle_peptide_paths_leaves_a_stop_distance_on_a_standard_table_transcript_a_difference() {
    let directory = scratch_of("organelle-leaves");
    let tables = codon_tables_of(&organelle_cache(&directory.join("cache")));
    let organelle = default_text(&[default_row(
        "r1",
        "MT:104",
        "-",
        "TXMT",
        "frameshift_variant",
        "IMPACT=HIGH;STRAND=1;HGVSp=PMT.1:p.Ala2GlyfsTer?",
    )]);
    let shifted = organelle.replace("STRAND=1;", "STRAND=1;HGVSc=TXMT:c.5del;HGVS_OFFSET=3;");
    let inserted = shifted
        .replace("\t-\tENSG1", "\tCA\tENSG1")
        .replace("STRAND=1", "STRAND=-1");
    let standard = default_text(&[default_row(
        "r1",
        "21:104",
        "-",
        "TXN",
        "frameshift_variant",
        "IMPACT=HIGH;STRAND=1;HGVSp=PN.1:p.Ala2GlyfsTer?",
    )]);
    let substitution = default_text(&[default_row(
        "r1",
        "MT:104",
        "T",
        "TXMT",
        "missense_variant",
        "IMPACT=MODERATE;STRAND=1;HGVSp=PMT.1:p.Ala2Gly",
    )]);
    let numbered = |text: &str| text.replace("fsTer?", "fsTer3");
    let absent = standard.replace("TXN", "TXOTHER");
    let on = |chromosome: &str, transcript: &str| {
        standard
            .replace("21:104", &format!("{chromosome}:104"))
            .replace("TXN", transcript)
    };
    let (plastid, strict, pachysolen) =
        (on("Pt", "TXPT"), on("St", "TXSTRICT"), on("Pa", "TXPACHY"));
    let without_value = organelle.replace(";HGVSp=PMT.1:p.Ala2GlyfsTer?", "");
    let lost_stop = organelle.replace("frameshift_variant", "stop_lost");
    let pairs = [
        (standard.clone(), numbered(&standard)),
        (absent.clone(), numbered(&absent)),
        (plastid.clone(), numbered(&plastid)),
        (strict.clone(), numbered(&strict)),
        (pachysolen.clone(), numbered(&pachysolen)),
        (
            shifted.clone(),
            numbered(&shifted.replace("HGVS_OFFSET=3", "HGVS_OFFSET=5")),
        ),
        (
            shifted.clone(),
            numbered(&shifted.replace("c.5del", "c.7del")),
        ),
        (
            organelle.replace("frameshift_variant\t-\t-\t-", "frameshift_variant\t5\t5\t2"),
            numbered(
                &organelle.replace("frameshift_variant\t-\t-\t-", "frameshift_variant\t5\t8\t3"),
            ),
        ),
        (
            organelle.replace("\t-\t-\t-\tIMPACT", "\t-\tgCa/ga\t-\tIMPACT"),
            numbered(&organelle.replace("\t-\t-\t-\tIMPACT", "\t-\tGca/ca\t-\tIMPACT")),
        ),
        (
            organelle.replace("frameshift_variant", "frameshift_variant,start_lost"),
            numbered(&organelle),
        ),
        (inserted.clone(), numbered(&inserted)),
        (
            substitution.clone(),
            substitution.replace("Ala2Gly", "Ala2Ser"),
        ),
        (organelle.clone(), without_value.clone()),
        (without_value, organelle.clone()),
        (
            organelle.clone(),
            organelle.replace("PMT.1:p.Ala2GlyfsTer?", "PMT.2:p.Ala2GlyfsTer3"),
        ),
        (lost_stop.clone(), numbered(&lost_stop)),
    ];
    for (reference_text, candidate_text) in &pairs {
        let (report, documented) = compare_downstream(
            reference_text,
            candidate_text,
            &tables,
            &[ORGANELLE_PEPTIDE],
        );
        assert_eq!(verdict_of(&report), "red", "{reference_text}");
        assert!(counted(&documented.applied).is_empty(), "{reference_text}");
        assert_eq!(differing(&report, "HGVSp"), 1, "{reference_text}");
    }
    // The standard-table transcript in the VCF format, whose entry names no chromosome.
    let inputs = ["21\t103\tr1\tCA\tC\t.\t.\t.".to_string()];
    let vcf = |distance: &str| {
        parse_vcf(
            &format!(
                "{ORGANELLE_CSQ_HEADER}21\t103\tr1\tCA\tC\t.\t.\tCSQ=-|frameshift_variant|Transcript|TXN|PN.1:p.Ala2GlyfsTer{distance}\n"
            ),
            &inputs,
        )
        .expect("the VCF line parses")
    };
    let (report, documented) = compared(&vcf("?"), &vcf("3"), &tables, &[ORGANELLE_PEPTIDE]);
    assert_eq!(verdict_of(&report), "red");
    assert!(counted(&documented.applied).is_empty());
}

/// Ensembl VEP lends a same-symbol HGNC gene's id to the transcripts of one buffer
/// (AnnotationType/Transcript.pm:290, :305), so the HGNC_ID of a transcript whose gene is
/// named by EntrezGene or RFAM follows the buffer's other records: one side empty and the
/// other an id, in either direction, is the class on a Transcript row whose SYMBOL_SOURCE is
/// not HGNC, counted once per row in every format; two differing ids, a differing id on an
/// HGNC-named gene, a value that is not an id and a row without a SYMBOL_SOURCE stay
/// differences.
#[test]
fn hgnc_id_copy_accepts_an_id_lent_to_a_gene_named_outside_hgnc() {
    let side = Side::default();
    let applied_once = count_map(&[(HGNC_ID_COPY, 1)]);
    let extra = "IMPACT=MODIFIER;STRAND=1;SYMBOL=NMRAL2P;SYMBOL_SOURCE=EntrezGene";
    let reference = default_text(&[default_row(
        "r1",
        "3:186010000",
        "C",
        "ENST00000824369",
        "intron_variant",
        extra,
    )]);
    let candidate = reference.replace(
        "SYMBOL_SOURCE=EntrezGene",
        "SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:52332",
    );
    let verdict = |reference: &str, candidate: &str, classes: &[&str]| {
        let (report, _) = compare_downstream(reference, candidate, &side, classes);
        verdict_of(&report).to_string()
    };
    let (report, documented) = compare_downstream(&reference, &candidate, &side, &[HGNC_ID_COPY]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), applied_once);
    assert_eq!(masked(&report, "Extra"), applied_once);
    assert_eq!(masked(&report, "HGNC_ID"), applied_once);
    assert_eq!(
        mask_examples(&report, HGNC_ID_COPY),
        [
            "3:186010000|C|ENST00000824369|Transcript: Extra expected 'IMPACT=MODIFIER;STRAND=1;SYMBOL=NMRAL2P;SYMBOL_SOURCE=EntrezGene', got 'IMPACT=MODIFIER;STRAND=1;SYMBOL=NMRAL2P;SYMBOL_SOURCE=EntrezGene;HGNC_ID=HGNC:52332'",
            "3:186010000|C|ENST00000824369|Transcript: HGNC_ID expected '<absent>', got 'HGNC:52332'",
        ]
    );
    assert_eq!(verdict(&reference, &candidate, &[]), "red");
    // The reference lends the id and the candidate prints none; and a cache before the
    // prefix writes the id as digits alone.
    assert_eq!(verdict(&candidate, &reference, &[HGNC_ID_COPY]), "green");
    assert_eq!(
        verdict(
            &reference.replace("EntrezGene", "RFAM"),
            &candidate.replace("EntrezGene", "RFAM"),
            &[HGNC_ID_COPY]
        ),
        "green"
    );
    assert_eq!(
        verdict(
            &reference,
            &candidate.replace("HGNC:52332", "52332"),
            &[HGNC_ID_COPY]
        ),
        "green"
    );
    for (reference_text, candidate_text) in [
        (
            reference.replace("EntrezGene", "HGNC"),
            candidate.replace("EntrezGene", "HGNC"),
        ),
        (
            candidate.clone(),
            candidate.replace("HGNC:52332", "HGNC:10101"),
        ),
        (
            reference.clone(),
            candidate.replace("HGNC:52332", "NMRAL2P"),
        ),
        (
            reference.replace(";SYMBOL_SOURCE=EntrezGene", ""),
            candidate.replace(";SYMBOL_SOURCE=EntrezGene", ""),
        ),
    ] {
        let (red, documented) =
            compare_downstream(&reference_text, &candidate_text, &side, &[HGNC_ID_COPY]);
        assert_eq!(verdict_of(&red), "red", "{reference_text}");
        assert!(counted(&documented.applied).is_empty(), "{reference_text}");
    }
    // The VCF format carries the field whatever the flags (Config.pm:501-507), the JSON
    // format under its renamed source key; the row is counted once in each.
    let inputs = ["3\t186010000\tr1\tA\tC\t.\t.\t.".to_string()];
    let vcf = |id: &str| {
        parse_vcf(
            &format!(
                "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|SYMBOL|SYMBOL_SOURCE|HGNC_ID\">\n3\t186010000\tr1\tA\tC\t.\t.\tCSQ=C|intron_variant|Transcript|ENST00000824369|NMRAL2P|EntrezGene|{id}\n"
            ),
            &inputs,
        )
        .expect("the VCF line parses")
    };
    let (report, documented) = compared(&vcf(""), &vcf("HGNC:52332"), &side, &[HGNC_ID_COPY]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), applied_once);
    assert_eq!(masked(&report, "HGNC_ID"), applied_once);
    let (report, _) = compared(&vcf(""), &vcf("HGNC:52332"), &side, &[]);
    assert_eq!(verdict_of(&report), "red");
    let json = |source: &str, id: Option<&str>| {
        let mut entry = vec![
            ("variant_allele", text("C")),
            ("transcript_id", text("ENST00000824369")),
            ("consequence_terms", list(vec![text("intron_variant")])),
            ("gene_symbol", text("NMRAL2P")),
            ("gene_symbol_source", text(source)),
        ];
        entry.extend(id.map(|id| ("hgnc_id", text(id))));
        let record = object(vec![
            ("input", text(&inputs[0])),
            ("id", text("r1")),
            ("most_severe_consequence", text("intron_variant")),
            ("transcript_consequences", list(vec![object(entry)])),
        ]);
        parse_json(&(dumped_value(&record) + "\n"), &inputs).expect("the JSON object parses")
    };
    let (report, documented) = compared(
        &json("EntrezGene", None),
        &json("EntrezGene", Some("HGNC:52332")),
        &side,
        &[HGNC_ID_COPY],
    );
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), applied_once);
    let (report, _) = compared(
        &json("HGNC", None),
        &json("HGNC", Some("HGNC:52332")),
        &side,
        &[HGNC_ID_COPY],
    );
    assert_eq!(verdict_of(&report), "red");
}

/// The genetic codes `crates/vep-core/src/codon.rs` reads, taken from its source: the
/// tables of the module by id (`CODES`, an id without a table left out), the id of the
/// module's last table, and the code the module lacks with its id.
struct CrateCodes {
    module: Vec<(i64, String)>,
    module_last: i64,
    lacked: (i64, String),
}

fn crate_codes() -> CrateCodes {
    let source = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../crates/vep-core/src/codon.rs"),
    )
    .expect("the codon module reads");
    let start = source
        .find("const CODES:")
        .expect("the module lists its codes");
    let end = start + source[start..].find("\n];").expect("the list ends");
    let entry = regex::Regex::new(r#"(?m)^    (?:(None)|code\(\s*b"([^"]{64})")"#).unwrap();
    let module = entry
        .captures_iter(&source[start..end])
        .enumerate()
        .filter_map(|(table, found)| {
            Some((
                i64::try_from(table).ok()?,
                found.get(2)?.as_str().to_string(),
            ))
        })
        .collect();
    let number = |name: &str| -> i64 {
        let stated = regex::Regex::new(&format!(r"const {name}: u8 = (\d+);")).unwrap();
        stated.captures(&source).expect("the constant")[1]
            .parse()
            .expect("a table id")
    };
    let lacked = regex::Regex::new(
        r#"const PACHYSOLEN_NUCLEAR_CODE: GeneticCode = GeneticCode \{\s*amino_acids: b"([^"]{64})""#,
    )
    .unwrap();
    CrateCodes {
        module,
        module_last: number("MODULE_LAST_TABLE"),
        lacked: (
            number("PACHYSOLEN_NUCLEAR"),
            lacked.captures(&source).expect("the code's constant")[1].to_string(),
        ),
    }
}

/// The amino acids of the 64 codons under `table` as `reading` reads it.
fn read_under(table: i64, reading: Reading) -> String {
    let code = genetic_code(table, reading).expect("a table is held");
    String::from_utf8(code.amino_acids().to_vec()).expect("letters")
}

/// The comparison holds a table for every genetic code crates/vep-core/src/codon.rs reads
/// and for no other id: the module's tables (`CODES`) with the same amino acids in the same
/// codon order under either engine's reading, and the code the module lacks
/// (`PACHYSOLEN_NUCLEAR_CODE`), which vep-rs reads as that constant and Ensembl VEP under
/// the module's last table (`module_table`). A cache's `codon_table` of 0 or none names
/// table 1, and a value that names no held table is kept as the cache spells it.
#[test]
fn the_genetic_codes_are_those_of_the_codon_module() {
    let held = crate_codes();
    let mut ids: Vec<i64> = held.module.iter().map(|(table, _)| *table).collect();
    ids.push(held.lacked.0);
    assert_eq!(held_tables(), ids);
    for (table, amino_acids) in &held.module {
        for reading in [Reading::EnsemblVep, Reading::VepRs] {
            assert_eq!(&read_under(*table, reading), amino_acids, "table {table}");
        }
    }
    let (lacked, amino_acids) = &held.lacked;
    assert_eq!(&read_under(*lacked, Reading::VepRs), amino_acids);
    let (last, last_amino_acids) = held.module.last().expect("a table");
    assert_eq!(*last, held.module_last);
    assert_eq!(&read_under(*lacked, Reading::EnsemblVep), last_amino_acids);
    assert_ne!(amino_acids, last_amino_acids);
    let table_of =
        |values: &[Value]| -> Vec<CodonTable> { values.iter().map(codon_table_id).collect() };
    assert_eq!(
        table_of(&[
            Value::Null,
            whole(0),
            text("0"),
            whole(1),
            text("1"),
            whole(2),
            text("2"),
            whole(24),
            whole(26),
            text("26")
        ]),
        [1, 1, 1, 1, 1, 2, 2, 24, 26, 26].map(CodonTable::Held)
    );
    assert_eq!(
        table_of(&[
            whole(7),
            text("8"),
            whole(17),
            whole(20),
            whole(25),
            whole(27),
            text("x"),
            Value::Bool(true)
        ]),
        ["7", "8", "17", "20", "25", "27", "x", "True"]
            .map(|spelt| CodonTable::Lacking(spelt.to_string()))
    );
}

/// A coding transcript of the three codons `codons` under the genetic code `table`.
fn transcript_under(table: i64, codons: &str) -> DownstreamTranscript {
    DownstreamTranscript {
        coding: codons.to_string(),
        utr: String::new(),
        table: codon_table_id(&whole(table)),
        peptide_length: None,
        peptide: None,
        model: None,
    }
}

/// A transcript under each genetic code the caches use beside those the human caches name
/// is read under that code, by either engine's reading: the middle codon of each case is one
/// its table reads unlike the standard code. The residues are the rows of
/// crates/vep-core/tests/data/genetic_codes.tsv: line 277 (table 4), 401 (table 6), 497
/// (table 9), 666 (table 12) and 757 (table 13), against lines 85, 81, 113, 90 and 117 of
/// table 1. Every other codon of each table is held to the table crates/vep-core/src/codon.rs
/// lists for it (`CODES`).
#[test]
fn a_transcript_under_a_genetic_code_of_another_genome_is_read_under_that_code() {
    let held = crate_codes();
    let every_codon: String = (0..64)
        .flat_map(|codon: usize| [codon / 16, codon / 4 % 4, codon % 4])
        .map(|base| ['T', 'C', 'A', 'G'][base])
        .collect();
    let cases = [
        (4, "TGA", 'W', '*'),
        (6, "TAA", 'Q', '*'),
        (9, "AAA", 'N', 'K'),
        (12, "CTG", 'S', 'L'),
        (13, "AGA", 'G', 'R'),
    ];
    for (table, codon, residue, standard) in cases {
        let coding = format!("ATG{codon}GGG");
        for reading in [Reading::EnsemblVep, Reading::VepRs] {
            let under = |table: i64| {
                let transcript = transcript_under(table, &coding);
                assert_eq!(transcript.table, CodonTable::Held(table));
                (
                    codon_residues(&transcript, (2, 2), reading),
                    codon_residues(&transcript, (1, 3), reading),
                    translate_until_stop(&coding, table, reading),
                )
            };
            let through = |residue: char| match residue {
                '*' => "M".to_string(),
                residue => format!("M{residue}G"),
            };
            assert_eq!(
                under(table),
                (
                    Some(residue.to_string()),
                    Some(format!("M{residue}G")),
                    through(residue)
                ),
                "table {table}"
            );
            assert_eq!(
                under(1),
                (
                    Some(standard.to_string()),
                    Some(format!("M{standard}G")),
                    through(standard)
                ),
                "table 1 beside table {table}"
            );
            let whole_table = transcript_under(table, &every_codon);
            let listed = &held
                .module
                .iter()
                .find(|(id, _)| *id == table)
                .expect("the crate lists the table")
                .1;
            assert_eq!(
                codon_residues(&whole_table, (1, 64), reading).as_ref(),
                Some(listed),
                "table {table}"
            );
        }
    }
}

/// Genetic code 26 is read as each engine reads it, and under neither is it a refusal:
/// vep-rs under the code's own table (crates/vep-core/src/codon.rs,
/// `PACHYSOLEN_NUCLEAR_CODE`, the one codon that differs from the standard code in its
/// second row of sixteen), Ensembl VEP under its module's last table, which holds the
/// standard code's residue for that codon (crates/vep-core/tests/data/genetic_codes.tsv
/// line 1178, as line 90 of table 1). The Downstream plugin's value follows: the
/// reference's is the module's reading (Downstream.pm:130-139) and the candidate's the
/// code's own.
#[test]
fn genetic_code_26_is_read_as_each_engine_reads_it() {
    let transcript = transcript_under(26, "ATGAAACTGTAA");
    assert_eq!(transcript.table, CodonTable::Held(26));
    let residues = |codons, reading| codon_residues(&transcript, codons, reading);
    assert_eq!(residues((3, 3), Reading::EnsemblVep).as_deref(), Some("L"));
    assert_eq!(residues((3, 3), Reading::VepRs).as_deref(), Some("A"));
    assert_eq!(
        residues((1, 4), Reading::EnsemblVep).as_deref(),
        Some("MKL*")
    );
    assert_eq!(residues((1, 4), Reading::VepRs).as_deref(), Some("MKA*"));
    let translated = |reading| translate_until_stop(&transcript.coding, 26, reading);
    assert_eq!(translated(Reading::EnsemblVep), "MKL");
    assert_eq!(translated(Reading::VepRs), "MKA");
    // The second codon replaced: each engine's plugin translates from that codon on.
    let value = |ensembl: bool| transcript.downstream_value(4, 6, "AAG", false, ensembl);
    assert_eq!(value(true).as_deref(), Some("KL"));
    assert_eq!(value(false).as_deref(), Some("KA"));
}

/// A transcript under genetic code 26 on the toy forward transcript's model, coding base c
/// at 21:1000+c. Codon 3 is the codon the code reads unlike the standard code; codons 2
/// and 4 read alike under both.
const CODE_26: &str = "ENST00000000207";
const CODE_26_CODING: &str = "ATGACCCTGGGTGCTGTCCGTGGTTAA";

/// The classes of a later allele's row read each reference residue as the engine that
/// writes it reads genetic code 26: the reference's residue, at the pair the first allele's
/// notation left, as Ensembl VEP reads the codon (its module's last table); the candidate's,
/// at the record's own pair, as the code itself reads it; and the BLOSUM62 the reference
/// alone writes is the matrix cell of Ensembl VEP's residue. The same row with the other
/// engine's residue stays a difference. The first allele's rows are written alike on both
/// sides.
#[test]
fn a_later_alleles_row_under_genetic_code_26_is_read_as_each_engine_reads_it() {
    let directory = scratch_of("code-26-later-allele");
    let cache = directory.join("cache");
    let transcript = modelled(
        CODE_26,
        1,
        &[(1, 42, 1001, 1042)],
        Modelled {
            coding: CODE_26_CODING,
            table: 26,
            ..Modelled::default()
        },
    );
    written(
        &cache,
        "transcripts/21/1-1000000.json",
        &dumped_value(&list(vec![transcript])),
    );
    let side = transcripts_of(&cache);
    // A row's name, Location, Allele and consequence, its CDS_position, Protein_position
    // and Amino_acids, and what its Extra holds after the STRAND.
    let row =
        |name: &str, location: &str, allele: &str, terms: &str, cells: [&str; 3], extra: &str| {
            let [cds, protein, amino_acids] = cells;
            with_coding_cells(
                &default_row(
                    name,
                    location,
                    allele,
                    CODE_26,
                    terms,
                    &format!("IMPACT=HIGH;STRAND=1{extra}"),
                ),
                &format!("-\t{cds}\t{protein}\t{amino_acids}\t-\t-"),
            )
        };
    let columns = ["CDS_position", "Protein_position", "Amino_acids"];
    // The columns the class takes of the later allele's row of `record`, whose two alleles
    // the rows spell `alleles`: the first allele's row carries `first` and an offset of
    // three bases, the later allele's row `reference` and `candidate`.
    let taken = |record: &str,
                 alleles: [&str; 2],
                 first: [&str; 3],
                 reference: [&str; 3],
                 candidate: [&str; 3]| {
        let start = record
            .split('\t')
            .nth(1)
            .and_then(|start| start.parse::<i64>().ok());
        let start = start.expect("the record's POS is a number");
        let location = format!("21:{start}-{}", start + 1);
        let text = |later: [&str; 3]| {
            default_text(&[
                row(
                    "v",
                    &location,
                    alleles[0],
                    "frameshift_variant",
                    first,
                    ";HGVS_OFFSET=3",
                ),
                row("v", &location, alleles[1], "frameshift_variant", later, ""),
            ])
        };
        let (report, _) = compare_record(
            &text(reference),
            &text(candidate),
            &side,
            &[SHARED_SHIFT],
            &[record],
        );
        columns
            .iter()
            .copied()
            .filter(|column| masked(&report, column) == count_map(&[(SHARED_SHIFT, 1)]))
            .collect::<Vec<_>>()
    };
    // An insertion after coding base 4 moved into codon 3: the reference's residue there is
    // Ensembl VEP's reading of the codon.
    let moved_into = "21\t1004\tv\tA\tAC,ACC\t.\t.\t.";
    let (own, first) = (["4-5", "2", "T/TX"], ["4-5", "2", "T/TX"]);
    assert_eq!(
        taken(moved_into, ["C", "CC"], first, ["7-8", "3", "L/TX"], own),
        columns
    );
    assert_eq!(
        taken(moved_into, ["C", "CC"], first, ["7-8", "3", "A/TX"], own),
        columns[..2]
    );
    // An insertion after coding base 8, inside codon 3, moved into codon 4: the candidate's
    // residue at its own pair is the code's own reading of the codon.
    let moved_out_of = "21\t1008\tv\tT\tTG,TA\t.\t.\t.";
    let (moved, first) = (["11-12", "4", "G/LX"], ["8-9", "3", "L/LX"]);
    assert_eq!(
        taken(moved_out_of, ["G", "A"], first, moved, ["8-9", "3", "A/LX"]),
        columns
    );
    assert_eq!(
        taken(moved_out_of, ["G", "A"], first, moved, ["8-9", "3", "L/LX"]),
        columns[..2]
    );
    // Three bases inserted between codons 2 and 3, read one base on at codon 3 alone: the
    // score is the matrix cell of Ensembl VEP's residue there.
    let between = "21\t1006\tb\tC\tCC,CAAA\t.\t.\t.";
    let scored = |residue: char| {
        let score = BLOSUM62_MATRIX[&(residue, 'K')];
        let location = "21:1006-1007";
        let text = |cells: [&str; 3], extra: &str| {
            default_text(&[
                row(
                    "b",
                    location,
                    "C",
                    "inframe_insertion",
                    ["6-7", "2-3", "-/X"],
                    ";HGVS_OFFSET=1",
                ),
                row("b", location, "AAA", "inframe_insertion", cells, extra),
            ])
        };
        let pair = format!("{residue}/K");
        let (report, _) = compare_record(
            &text(["7-8", "3", pair.as_str()], &format!(";BLOSUM62={score}")),
            &text(["6-7", "2-3", "-/K"], ""),
            &side,
            &[BLOSUM62_SHARED_SHIFT],
            &[between],
        );
        masked(&report, "BLOSUM62") == count_map(&[(BLOSUM62_SHARED_SHIFT, 1)])
    };
    assert_ne!(BLOSUM62_MATRIX[&('L', 'K')], BLOSUM62_MATRIX[&('A', 'K')]);
    assert!(scored('L'));
    assert!(!scored('A'));
}

const LACKING_TABLE_7: &str = "the cache's codon_table 7 names no genetic code the comparison \
     holds a table for (it holds 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, \
     26)";

/// A row that asks for a transcript under a genetic code no table is held for ends the
/// comparison with an error that names the transcript and the table, under each class that
/// reads the table: the row is read under no other code. A row of such a transcript that no
/// class asks about is compared as any other.
#[test]
fn a_transcript_under_a_genetic_code_without_a_table_ends_the_comparison_by_name() {
    let directory = scratch_of("code-without-a-table");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let refused = format!("transcript {TOY_TABLE_7}: {LACKING_TABLE_7}");
    let coding_end = downstream_row("21:124-127", "-", TOY_TABLE_7, "24-27", "GPK", "GGPNS", 4);
    let own_codon = downstream_row(
        "21:106",
        "-",
        TOY_TABLE_7,
        "6",
        "LLLSVVKGPNS",
        "KLLLSVVKGPNS",
        5,
    );
    for (row, class) in [
        (coding_end, UTR_FRAME_PAST_CODING_END),
        (own_codon, START_PAST_OWN_CODON),
    ] {
        let (reference, candidate) = downstream_texts(&[row]);
        assert_eq!(
            refusal_of(&reference, &candidate, &side, &[class]),
            refused,
            "{class}"
        );
        // Two equal rows ask no class: the comparison is green.
        let (report, _) = compare_downstream(&reference, &reference, &side, &[class]);
        assert_eq!(verdict_of(&report), "green", "{class}");
    }
    let tables = codon_tables_of(&organelle_cache(&directory.join("organelle")));
    let (reference, candidate) = organelle_frameshift_texts();
    let on_no_table = |text: &str| text.replace("MT:104", "No:104").replace("TXMT", "TXNONE");
    assert_eq!(
        refusal_of(
            &on_no_table(&reference),
            &on_no_table(&candidate),
            &tables,
            &[ORGANELLE_PEPTIDE]
        ),
        format!("transcript TXNONE: {LACKING_TABLE_7}")
    );
}

fn coding_end_rows() -> Vec<DownstreamRow> {
    vec![
        downstream_row("21:124-127", "-", TOY, "24-27", "GPK", "GGPNS", 4),
        downstream_row(
            "21:124-127",
            "-",
            &format!("{TOY}.3"),
            "24-27",
            "GPK",
            "GGPNS",
            4,
        ),
        downstream_row("21:1024-1027", "-", TOY_FORWARD, "24-27", "GPK", "GGPNS", 4).offset(2),
        downstream_row("21:2016-2019", "-", TOY_REVERSE, "24-27", "GPK", "GGPNS", 4)
            .offset(-2)
            .strand(-1),
        downstream_row("21:224-227", "-", TOY_UTR_STOP, "24-27", "", "GDPNS", 4),
        downstream_row("21:324-327", "-", TOY_TABLE_2, "24-27", "VTKMA", "GWP", 2),
    ]
}

/// Downstream.pm:113 rounds the first deleted base's coordinate down to a codon and :115
/// reads the changed coding sequence from that offset, which lies past its end when a
/// deletion starting on a codon's third base runs to the last coding base: the reference
/// translates the 3' UTR from its first base (:118-119), the candidate the two coding bases
/// before the deletion with the UTR. The class recomputes both values from the candidate
/// cache and takes the row once in every format, under a Feature carrying a version, beside
/// an HGVS_OFFSET of 2 or -2, which moves the span past the last coding base, and on a
/// transcript whose UTR opens with a stop codon, where the reference's value is empty. Both
/// values follow the transcript's codon table. Without the cache the class accepts nothing.
#[test]
fn downstream_utr_frame_past_coding_end_accepts_a_deletion_running_to_the_last_coding_base() {
    let directory = scratch_of("coding-end-accepts");
    let cache = downstream_cache(&directory.join("cache"));
    let (reference, candidate) = downstream_texts(&coding_end_rows());
    let side = transcripts_of(&cache);
    let no_cache = Side::default();
    let all_six = count_map(&[(UTR_FRAME_PAST_CODING_END, 6)]);
    let (report, documented) =
        compare_downstream(&reference, &candidate, &side, &[UTR_FRAME_PAST_CODING_END]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), all_six);
    assert_eq!(masked(&report, "DownstreamProtein"), all_six);
    assert_eq!(masked(&report, "Extra"), all_six);
    assert_eq!(
        mask_examples(&report, UTR_FRAME_PAST_CODING_END)[0],
        format!("21:124-127|-|{TOY}|Transcript: DownstreamProtein expected 'GPK', got 'GGPNS'")
    );
    // The class without the cache and the cache without the class set nothing aside; beside
    // the other Downstream class the rows stay this one's.
    let (red, documented) = compare_downstream(
        &reference,
        &candidate,
        &no_cache,
        &[UTR_FRAME_PAST_CODING_END],
    );
    assert_eq!(verdict_of(&red), "red");
    assert!(counted(&documented.applied).is_empty());
    assert_eq!(differing(&red, "DownstreamProtein"), 6);
    let (red, _) = compare_downstream(&reference, &candidate, &side, &[]);
    assert_eq!(verdict_of(&red), "red");
    let (report, documented) = compare_downstream(
        &reference,
        &candidate,
        &side,
        &[UTR_FRAME_PAST_CODING_END, START_PAST_OWN_CODON],
    );
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), all_six);
    // The tab format carries the field as a column, the VCF as a CSQ field, the JSON under
    // its lower-cased key with the span as cds_start and cds_end; the row is counted once in
    // each.
    let record = "21\t123\tv\tGTTAA\tG\t.\t.\t.";
    let expected = downstream_formats(record, "21:124-127", "-", TOY, "24-27", Some("GPK"), 4);
    let actual = downstream_formats(record, "21:124-127", "-", TOY, "24-27", Some("GGPNS"), 4);
    assert!(every_format_reads(
        &applied_by_format(&expected, &actual, &side, &[UTR_FRAME_PAST_CODING_END]),
        "green",
        &[(UTR_FRAME_PAST_CODING_END, 1)]
    ));
    for (fmt, verdict, _) in
        applied_by_format(&expected, &actual, &no_cache, &[UTR_FRAME_PAST_CODING_END])
    {
        assert_eq!(verdict, "red", "{fmt}");
    }
    // A transcript the cache files under two chromosomes, one coding sequence under each
    // one's own mapper pairs, is read without a model by an entry that states no Location,
    // so the row is the class's in every format; one filed with two coding sequences is left
    // out for such an entry.
    let filed = |root: &Path, chromosome: &str, start: i64, coding: &str| {
        let transcript = modelled(
            TOY,
            1,
            &[(1, 42, start, start + 41)],
            Modelled {
                coding,
                ..Modelled::default()
            },
        );
        written(
            root,
            &format!("transcripts/{chromosome}/1-1000000.json"),
            &dumped_value(&list(vec![transcript])),
        );
    };
    let record = "X\t1023\tv\tGTTAA\tG\t.\t.\t.";
    let expected = downstream_formats(record, "X:1024-1027", "-", TOY, "24-27", Some("GPK"), 4);
    let actual = downstream_formats(record, "X:1024-1027", "-", TOY, "24-27", Some("GGPNS"), 4);
    let (two_models, two_sequences) = (
        directory.join("two_models"),
        directory.join("two_sequences"),
    );
    let other_stop = format!("{}TGA", &TOY_CODING[..TOY_CODING.len() - 3]);
    filed(&two_models, "X", 1001, TOY_CODING);
    filed(&two_models, "Y", 501, TOY_CODING);
    filed(&two_sequences, "X", 1001, TOY_CODING);
    filed(&two_sequences, "Y", 501, &other_stop);
    for ((fmt, entries), (_, others)) in expected.iter().zip(&actual) {
        let row = entries
            .iter()
            .find(|entry| entry.feature == TOY)
            .expect("the format has the row");
        let read = transcripts_of(&two_models);
        let transcript = sequences_of(&read)
            .transcript(row)
            .expect("the cache reads")
            .expect("the cache holds the transcript");
        assert_eq!(transcript.model.is_none(), *fmt != "tab", "{fmt}");
        let (report, documented) = compared(entries, others, &read, &[UTR_FRAME_PAST_CODING_END]);
        assert_eq!(verdict_of(&report), "green", "{fmt}");
        assert_eq!(
            counted(&documented.applied),
            count_map(&[(UTR_FRAME_PAST_CODING_END, 1)]),
            "{fmt}"
        );
        let left_out = DownstreamSequences::new(&two_sequences)
            .expect("the cache has a transcripts directory")
            .transcript(row)
            .expect("the cache reads")
            .is_none();
        assert_eq!(left_out, *fmt != "tab", "{fmt}");
    }
}

/// The command line reads the rows of a deletion running to the last coding base: the class
/// is unused without --json-cache and applied with it, and a cache without a transcripts
/// directory is refused.
#[test]
fn downstream_utr_frame_past_coding_end_reads_the_cache_the_command_line_names() {
    let directory = scratch_of("coding-end-command-line");
    let cache = downstream_cache(&directory.join("cache"));
    let (reference, candidate) = downstream_texts(&coding_end_rows());
    let all_six = count_map(&[(UTR_FRAME_PAST_CODING_END, 6)]);
    let arguments = with(
        &default_arguments(&directory, &reference, &candidate),
        &["--masks", UTR_FRAME_PAST_CODING_END],
    );
    assert_eq!(run(&arguments).code, 1);
    assert_eq!(
        texts_of(member(&report_of(&directory), &["masks", "unused"])),
        [UTR_FRAME_PAST_CODING_END]
    );
    assert_eq!(
        run(&with(&arguments, &["--json-cache", path_text(&cache)])).code,
        0
    );
    assert_eq!(
        counts_of(member(&report_of(&directory), &["masks", "applied"])),
        all_six
    );
    let nowhere = directory.join("nowhere");
    let ran = run(&with(&arguments, &["--json-cache", path_text(&nowhere)]));
    assert_eq!(ran.code, 2);
    assert_eq!(
        ran.stderr,
        format!(
            "ERROR: [compare_vep_formats] --json-cache: {} has no transcripts/ directory\n",
            path_text(&nowhere)
        )
    );
}

/// The class takes the pair the mechanism writes on the row it reaches and nothing wider.
/// The same pair on a deletion that ends before the last coding base or starts on a codon's
/// first or second base, on a deletion with inserted bases, and on a transcript without a 3'
/// UTR or absent from the cache stays a difference with the class unused; so
/// do the pair the standard code gives on a transcript under table 2, a value planted on
/// either side, a value one side lacks, a differing ProteinLengthChange or CDS_position, a
/// consequence set without frameshift_variant, with a splice term or differing between the
/// entries, and the row of a deletion that ends inside the coding sequence with that row's
/// own pair. The pair is the class's on the record's own coordinates alone: on a column
/// printed two bases before, beside the HGVS_OFFSET that moves it there or beside none, it
/// stays a difference.
#[test]
fn downstream_utr_frame_past_coding_end_leaves_a_row_the_mechanism_does_not_reach_a_difference() {
    let directory = scratch_of("coding-end-leaves");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let base = || downstream_row("21:124-127", "-", TOY, "24-27", "GPK", "GGPNS", 4);
    let texts = |row: DownstreamRow| downstream_texts(&[row]);
    let (reference, candidate) = texts(base());
    let (report, _) =
        compare_downstream(&reference, &candidate, &side, &[UTR_FRAME_PAST_CODING_END]);
    assert_eq!(verdict_of(&report), "green");
    let refused = [
        (
            "ends before the last coding base",
            texts(base().at("21:121-124", "21-24")),
        ),
        (
            "starts on a codon's first base",
            texts(base().at("21:125-127", "25-27")),
        ),
        (
            "starts on a codon's second base",
            texts(base().at("21:126-127", "26-27")),
        ),
        (
            "the column two bases before without an offset",
            texts(base().at("21:122-125", "22-25")),
        ),
        (
            "the column two bases before with the offset",
            texts(
                base()
                    .at("21:1022-1025", "22-25")
                    .feature(TOY_FORWARD)
                    .offset(2),
            ),
        ),
        ("inserted bases", texts(base().allele("GC"))),
        ("no 3' UTR", texts(base().feature(TOY_NO_UTR))),
        (
            "the standard code's pair under table 2",
            texts(base().feature(TOY_TABLE_2).values("VTKIA", "G", 0)),
        ),
        (
            "absent from the cache",
            texts(base().feature("ENST00000000199")),
        ),
        ("a planted candidate value", texts(base().vep_rs("GGPNT"))),
        ("a planted reference value", texts(base().vep("GPR"))),
        (
            "the reference lacks the field",
            (
                reference.replace("DownstreamProtein=GPK;", ""),
                candidate.clone(),
            ),
        ),
        (
            "the candidate lacks the field",
            (
                reference.clone(),
                candidate.replace("DownstreamProtein=GGPNS;", ""),
            ),
        ),
        (
            "a differing ProteinLengthChange",
            (
                reference.clone(),
                candidate.replace("ProteinLengthChange=4", "ProteinLengthChange=5"),
            ),
        ),
        (
            "a differing CDS_position",
            (
                reference.clone(),
                candidate.replace("\t24-27\t", "\t23-27\t"),
            ),
        ),
        (
            "no frameshift_variant",
            texts(base().consequence("inframe_deletion")),
        ),
        (
            "a splice term",
            texts(base().consequence("frameshift_variant,splice_region_variant")),
        ),
        (
            "differing consequence sets",
            (
                reference.clone(),
                candidate.replace("frameshift_variant", "frameshift_variant,stop_lost"),
            ),
        ),
        (
            "a deletion ending inside the coding sequence",
            texts(
                base()
                    .at("21:106-107", "6-7")
                    .values("SCCPWLRAQIA", "NSCCPWLRAQIA", 5),
            ),
        ),
    ];
    for (label, (reference_text, candidate_text)) in &refused {
        let (red, documented) = compare_downstream(
            reference_text,
            candidate_text,
            &side,
            &[UTR_FRAME_PAST_CODING_END],
        );
        assert_eq!(verdict_of(&red), "red", "{label}");
        assert!(counted(&documented.applied).is_empty(), "{label}");
    }
    // A row without the value in the other formats; the reference's side is refused alike.
    let record = "21\t123\tv\tGTTAA\tG\t.\t.\t.";
    let formats = |protein: Option<&str>| {
        downstream_formats(record, "21:124-127", "-", TOY, "24-27", protein, 4)
    };
    let without = formats(None);
    for (expected, actual) in [
        (formats(Some("GPK")), without.clone()),
        (without, formats(Some("GGPNS"))),
    ] {
        assert!(every_format_reads(
            &applied_by_format(&expected, &actual, &side, &[UTR_FRAME_PAST_CODING_END]),
            "red",
            &[]
        ));
    }
}

fn own_codon_rows() -> Vec<DownstreamRow> {
    vec![
        downstream_row("21:106", "-", TOY, "6", "LLLSVVKGPNS", "KLLLSVVKGPNS", 5),
        downstream_row(
            "21:106-107",
            "-",
            TOY,
            "6-7",
            "SCCPWLRAQIA",
            "NSCCPWLRAQIA",
            5,
        ),
        downstream_row("21:106", "GC", TOY, "6", "HSCCPWLRAQIA", "KHSCCPWLRAQIA", 6),
        downstream_row("22:206", "-", SHORT, "6", "KVK", "", -3),
        downstream_row("21:206", "-", TOY_NO_UTR, "6", "LLLSVV", "KLLLSVV", 0),
        downstream_row("21:306", "AC", TOY, "6", "YSCCPWLRAQIA", "KYSCCPWLRAQIA", 6).strand(-1),
        downstream_row(
            "21:1006",
            "-",
            TOY_FORWARD,
            "6",
            "LLLSVVKGPNS",
            "KLLLSVVKGPNS",
            5,
        )
        .offset(1),
        downstream_row(
            "21:2037",
            "-",
            TOY_REVERSE,
            "6",
            "LLLSVVKGPNS",
            "KLLLSVVKGPNS",
            5,
        )
        .offset(-1)
        .strand(-1),
        downstream_row(
            "21:406",
            "-",
            TOY_TABLE_2,
            "6",
            "LLLSVVKWP",
            "KLLLSVVKWP",
            3,
        ),
    ]
}

/// Downstream.pm:112-115 start the downstream translation at the first changed base's own
/// coordinate rounded down to a codon, the codon after the touched one when that base is a
/// codon's third; the candidate starts at the touched codon, so the reference's value is
/// the candidate's less its first residue. The class recomputes both values from the
/// candidate cache and takes the row once in every format, on a transcript without a 3'
/// UTR, on a reverse-strand row, whose Allele it reads reverse-complemented, and beside an
/// HGVS_OFFSET of 1 or -1 that moves the deletion of base 6 to base 7. Both values follow
/// the transcript's codon table. Without the cache the class accepts nothing.
#[test]
fn downstream_start_past_own_codon_accepts_a_translation_started_one_codon_late() {
    let directory = scratch_of("own-codon-accepts");
    let cache = downstream_cache(&directory.join("cache"));
    let rows = own_codon_rows();
    let (reference, candidate) = downstream_texts(&rows);
    let side = transcripts_of(&cache);
    let no_cache = Side::default();
    let all_nine = count_map(&[(START_PAST_OWN_CODON, 9)]);
    let (report, documented) =
        compare_downstream(&reference, &candidate, &side, &[START_PAST_OWN_CODON]);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(counted(&documented.applied), all_nine);
    assert_eq!(masked(&report, "DownstreamProtein"), all_nine);
    assert_eq!(masked(&report, "Extra"), all_nine);
    assert_eq!(
        mask_examples(&report, START_PAST_OWN_CODON)[0],
        format!(
            "21:106|-|{TOY}|Transcript: DownstreamProtein expected 'LLLSVVKGPNS', got 'KLLLSVVKGPNS'"
        )
    );
    let (red, documented) =
        compare_downstream(&reference, &candidate, &no_cache, &[START_PAST_OWN_CODON]);
    assert_eq!(verdict_of(&red), "red");
    assert!(counted(&documented.applied).is_empty());
    assert_eq!(differing(&red, "DownstreamProtein"), 9);
    let (red, _) = compare_downstream(&reference, &candidate, &side, &[]);
    assert_eq!(verdict_of(&red), "red");
    // The other formats count each row once: a kept residue, inserted bases, and the empty
    // value, an empty tab cell, an empty CSQ cell and an empty JSON string.
    for (record, location, allele, feature, cds, vep, vep_rs, change) in [
        (
            "21\t105\tv\tAA\tA\t.\t.\t.",
            "21:106",
            "-",
            TOY,
            "6",
            "LLLSVVKGPNS",
            "KLLLSVVKGPNS",
            5,
        ),
        (
            "21\t106\tv\tA\tGC\t.\t.\t.",
            "21:106",
            "GC",
            TOY,
            "6",
            "HSCCPWLRAQIA",
            "KHSCCPWLRAQIA",
            6,
        ),
        (
            "22\t205\tv\tAC\tA\t.\t.\t.",
            "22:206",
            "-",
            SHORT,
            "6",
            "KVK",
            "",
            -3,
        ),
    ] {
        let expected =
            downstream_formats(record, location, allele, feature, cds, Some(vep), change);
        let actual =
            downstream_formats(record, location, allele, feature, cds, Some(vep_rs), change);
        assert!(
            every_format_reads(
                &applied_by_format(&expected, &actual, &side, &[START_PAST_OWN_CODON]),
                "green",
                &[(START_PAST_OWN_CODON, 1)]
            ),
            "{feature} {allele}"
        );
        for (fmt, verdict, _) in
            applied_by_format(&expected, &actual, &no_cache, &[START_PAST_OWN_CODON])
        {
            assert_eq!(verdict, "red", "{fmt} {feature} {allele}");
        }
    }
}

/// The command line with both Downstream classes named: each row is one class's, unused
/// without --json-cache and applied with it.
#[test]
fn downstream_start_past_own_codon_reads_the_cache_the_command_line_names() {
    let directory = scratch_of("own-codon-command-line");
    let cache = downstream_cache(&directory.join("cache"));
    let mut both = own_codon_rows();
    both.push(downstream_row(
        "21:124-127",
        "-",
        TOY,
        "24-27",
        "GPK",
        "GGPNS",
        4,
    ));
    let (reference, candidate) = downstream_texts(&both);
    let masks = format!("{UTR_FRAME_PAST_CODING_END},{START_PAST_OWN_CODON}");
    let arguments = with(
        &default_arguments(&directory, &reference, &candidate),
        &["--masks", &masks],
    );
    assert_eq!(run(&arguments).code, 1);
    assert_eq!(
        texts_of(member(&report_of(&directory), &["masks", "unused"])),
        [START_PAST_OWN_CODON, UTR_FRAME_PAST_CODING_END]
    );
    assert_eq!(
        run(&with(&arguments, &["--json-cache", path_text(&cache)])).code,
        0
    );
    assert_eq!(
        counts_of(member(&report_of(&directory), &["masks", "applied"])),
        count_map(&[(START_PAST_OWN_CODON, 9), (UTR_FRAME_PAST_CODING_END, 1)])
    );
}

/// The class takes the pair the mechanism writes on the row it reaches and nothing wider.
/// The same pair on a deletion starting on a codon's first or second base and on a
/// transcript absent from the cache stays a difference with the class
/// unused; so do the pair the standard code gives on a transcript under table 2, a value
/// planted on either side, a value one side lacks, a differing ProteinLengthChange, a
/// consequence set without frameshift_variant or with a splice term, a reverse-strand row's
/// pair on a forward-strand row, and the row of a deletion running to the last coding base,
/// the other Downstream class's. An insertion prints the two bases around it as its
/// CDS_position and both engines start its translation at the codon taking the first
/// inserted base; read as a replacement of bases 6 and 7 the recomputed pair is neither
/// that value against a planted one nor the pair of the replacement of base 6. A deletion
/// of base 4 that an HGVS_OFFSET of 2 moves to base 6 carries the pair because the plugin
/// read the shifted coordinates, and stays a difference under this class.
#[test]
fn downstream_start_past_own_codon_leaves_a_row_the_mechanism_does_not_reach_a_difference() {
    let directory = scratch_of("own-codon-leaves");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let base = || downstream_row("21:106", "-", TOY, "6", "LLLSVVKGPNS", "KLLLSVVKGPNS", 5);
    let texts = |row: DownstreamRow| downstream_texts(&[row]);
    let green = |reference: &str, candidate: &str| {
        let (report, _) = compare_downstream(reference, candidate, &side, &[START_PAST_OWN_CODON]);
        verdict_of(&report) == "green"
    };
    let (reference, candidate) = texts(base());
    assert!(green(&reference, &candidate));
    let (short_reference, short_candidate) = texts(
        base()
            .at("22:206", "6")
            .feature(SHORT)
            .values("KVK", "", -3),
    );
    assert!(green(&short_reference, &short_candidate));
    let refused = [
        (
            "starts on a codon's first base",
            texts(base().at("21:107", "7")),
        ),
        (
            "starts on a codon's second base",
            texts(base().at("21:105", "5")),
        ),
        (
            "the column two bases before without an offset",
            texts(base().at("21:104", "4")),
        ),
        (
            "the column two bases before with the offset",
            texts(base().at("21:1004", "4").feature(TOY_FORWARD).offset(2)),
        ),
        (
            "an insertion carrying a replacement's pair",
            texts(base().at("21:106-107", "6-7").allele("GC").values(
                "HSCCPWLRAQIA",
                "KHSCCPWLRAQIA",
                6,
            )),
        ),
        (
            "an insertion's value against a planted one",
            texts(base().at("21:106-107", "6-7").allele("GC").values(
                "ALLLSVVKGPNS",
                "KALLLSVVKGPNS",
                6,
            )),
        ),
        (
            "a reverse-strand row's pair on a forward-strand row",
            texts(
                base()
                    .allele("AC")
                    .values("YSCCPWLRAQIA", "KYSCCPWLRAQIA", 6),
            ),
        ),
        (
            "the standard code's pair under table 2",
            texts(base().feature(TOY_TABLE_2).values("LLLSVVK", "KLLLSVVK", 1)),
        ),
        (
            "absent from the cache",
            texts(base().feature("ENST00000000199")),
        ),
        (
            "a planted candidate value",
            texts(base().vep_rs("KLLLSVVKGPNT")),
        ),
        (
            "a planted reference value",
            texts(base().vep("LLLSVVKGPNT")),
        ),
        (
            "the reference lacks the field",
            (
                reference.replace("DownstreamProtein=LLLSVVKGPNS;", ""),
                candidate.clone(),
            ),
        ),
        (
            "the candidate lacks the field",
            (
                reference.clone(),
                candidate.replace("DownstreamProtein=KLLLSVVKGPNS;", ""),
            ),
        ),
        (
            "the candidate lacks the field where its value is empty",
            (
                short_reference.clone(),
                short_candidate.replace("DownstreamProtein=;", ""),
            ),
        ),
        (
            "a differing ProteinLengthChange",
            (
                reference.clone(),
                candidate.replace("ProteinLengthChange=5", "ProteinLengthChange=4"),
            ),
        ),
        (
            "no frameshift_variant",
            texts(base().consequence("inframe_deletion")),
        ),
        (
            "a splice term",
            texts(base().consequence("frameshift_variant,splice_donor_variant")),
        ),
        (
            "a deletion running to the last coding base",
            texts(base().at("21:124-127", "24-27").values("GPK", "GGPNS", 4)),
        ),
    ];
    for (label, (reference_text, candidate_text)) in &refused {
        let (red, documented) = compare_downstream(
            reference_text,
            candidate_text,
            &side,
            &[START_PAST_OWN_CODON],
        );
        assert_eq!(verdict_of(&red), "red", "{label}");
        assert!(counted(&documented.applied).is_empty(), "{label}");
    }
    // A row without the value in the tab and JSON formats, which tell it from the empty
    // value.
    let record = "22\t205\tv\tAC\tA\t.\t.\t.";
    let expected = downstream_formats(record, "22:206", "-", SHORT, "6", Some("KVK"), -3);
    let without = downstream_formats(record, "22:206", "-", SHORT, "6", None, -3);
    for (fmt, verdict, applied) in
        applied_by_format(&expected, &without, &side, &[START_PAST_OWN_CODON])
    {
        if fmt != "vcf" {
            assert_eq!(verdict, "red", "{fmt}");
            assert!(applied.is_empty(), "{fmt}");
        }
    }
}

/// Ensembl VEP writes `start_lost` beside `start_retained_variant` on one allele and the
/// candidate writes the set without `start_lost`: the registry's `start_cooccurrence_swap`
/// sets that pair of consequence sets aside, and a frameshift row carrying it has the
/// Downstream values of any other row of its class. Each Downstream class takes its row
/// under that pair of sets, as `plugin_hgvs_shifted_cds_coordinates` takes a row the HGVS
/// shift changes, and under no other difference: a candidate that keeps `start_lost` and
/// drops `start_retained_variant`, one that drops both, a further term on one side, a
/// reference without both start terms and a splice term beside the pair stay differences.
#[test]
fn downstream_classes_take_the_start_term_difference_the_registry_sets_aside() {
    let directory = scratch_of("start-terms");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let both_start_terms = "frameshift_variant,start_lost,start_retained_variant";
    let without_start_lost = "frameshift_variant,start_retained_variant";
    let judged = |reference_terms: &str, candidate_terms: &str, row: &DownstreamRow| {
        let entry = |value: &str, terms: &str| {
            default_entries(&default_text(&[row.clone().consequence(terms).line(value)])).remove(0)
        };
        let (e, a) = (
            entry(&row.vep, reference_terms),
            entry(&row.vep_rs, candidate_terms),
        );
        let context = RuleContext {
            record_rows: None,
            candidate_rows: None,
            declared: &[],
            assembly: None,
            dbnsfp_columns: &[],
            side: &side,
        };
        let accepting: Vec<&str> = tables()
            .transcript_reading_classes()
            .iter()
            .copied()
            .filter(|class| {
                field_class_accepts(
                    class,
                    &e,
                    "DownstreamProtein",
                    &row.vep,
                    &row.vep_rs,
                    Some(&a),
                    &context,
                )
                .expect("the class reads the row")
            })
            .collect();
        let found = downstream_class(
            &e,
            &a,
            "DownstreamProtein",
            &row.vep,
            &row.vep_rs,
            sequences_of(&side),
            None,
        )
        .expect("the cache reads");
        // A row one cause changes: the class alone, none to declare beside it.
        let found = found.map(|(class, second)| {
            assert_eq!(second, None);
            class
        });
        assert_eq!(accepting, found.as_deref().into_iter().collect::<Vec<_>>());
        found
    };
    let own_codon = downstream_row("22:206", "-", SHORT, "6", "KVK", "", -3);
    let coding_end = downstream_row("21:124-127", "-", TOY, "24-27", "GPK", "GGPNS", 4);
    let shifted = downstream_row(
        "21:1005",
        "-",
        TOY_FORWARD,
        "5",
        "LLLSVVKGPNS",
        "KLLLSVVKGPNS",
        5,
    )
    .offset(2);
    let class = |name: &str| Some(name.to_string());
    assert_eq!(
        judged("frameshift_variant", "frameshift_variant", &own_codon),
        class(START_PAST_OWN_CODON)
    );
    assert_eq!(
        judged(both_start_terms, without_start_lost, &own_codon),
        class(START_PAST_OWN_CODON)
    );
    assert_eq!(
        judged(both_start_terms, without_start_lost, &coding_end),
        class(UTR_FRAME_PAST_CODING_END)
    );
    assert_eq!(
        judged("frameshift_variant", "frameshift_variant", &shifted),
        class(SHIFTED_CDS_COORDINATES)
    );
    assert_eq!(
        judged(both_start_terms, without_start_lost, &shifted),
        class(SHIFTED_CDS_COORDINATES)
    );
    let with_stop_lost = format!("{without_start_lost},stop_lost");
    let reference_with_splice = format!("{both_start_terms},splice_region_variant");
    let candidate_with_splice = format!("{without_start_lost},splice_region_variant");
    for (reference_terms, candidate_terms) in [
        (both_start_terms, "frameshift_variant,start_lost"),
        (both_start_terms, "frameshift_variant"),
        (both_start_terms, with_stop_lost.as_str()),
        ("frameshift_variant,start_lost", "frameshift_variant"),
        (without_start_lost, both_start_terms),
        (
            reference_with_splice.as_str(),
            candidate_with_splice.as_str(),
        ),
    ] {
        assert_eq!(
            judged(reference_terms, candidate_terms, &own_codon),
            None,
            "{reference_terms} {candidate_terms}"
        );
        assert_eq!(
            judged(reference_terms, candidate_terms, &shifted),
            None,
            "{reference_terms} {candidate_terms}"
        );
    }
}
