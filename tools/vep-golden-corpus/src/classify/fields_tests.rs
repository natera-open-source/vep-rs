// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the dispatch of field classes and of its three classes that read no
//! cache: the fields of a covered plugin on a record `--minimal` split, the motif score of
//! a reverse-strand motif under an allele of several bases, and the cell type of a name
//! that holds a comma.

use super::super::sequence_edges::HGVS_SHIFT_SEQUENCE_EDGE_WINDOW;
use super::super::spliceai::SPLICEAI_MULTI_ENTRY_CELL;
use super::*;

fn key(location: &str, allele: &str, feature_type: &str) -> FieldKey {
    FieldKey {
        location: location.to_string(),
        allele: allele.to_string(),
        feature: "ENSF1".to_string(),
        feature_type: feature_type.to_string(),
        consequences: "intron_variant".to_string(),
    }
}

fn question<'a>(
    name: &'a str,
    key: &'a FieldKey,
    vep: &'a str,
    vep_rs: &'a str,
) -> FieldQuestion<'a> {
    FieldQuestion {
        name,
        key,
        vep,
        vep_rs,
        strand: "-",
        flags: &[],
        perl_row: None,
        rust_row: None,
        minimal_split_record: false,
        genomic_span: None,
    }
}

fn class(question: &FieldQuestion<'_>) -> Option<&'static str> {
    let named = classify_field_divergence(question, &FieldModels::default());
    named.map(|(class, second)| {
        assert_eq!(second, None);
        class
    })
}

#[test]
fn the_dispatch_asks_twenty_seven_classes_and_its_own_three_first() {
    assert_eq!(FIELD_CLASSES.len(), 27);
    let transcript = key("21:100", "G", "Transcript");
    // A field no class reads is asked of every class and taken by none.
    let note = question("NOTE", &transcript, "x", "y");
    for asked in FIELD_CLASSES {
        assert_eq!(asked(&note, &FieldModels::default()), Verdict::Pass);
    }
    assert_eq!(class(&note), None);
    // The class of a split record is asked before any other.
    let motif = key("21:100-101", "AG", "MotifFeature");
    let mut both = question("MOTIF_SCORE_CHANGE", &motif, "-0.092", "-0.136");
    both.strand = "-1";
    assert_eq!(class(&both), Some(MOTIF_SCORE_REVERSE_STRAND_MNV));
    both.minimal_split_record = true;
    assert_eq!(class(&both), Some(MOTIF_SCORE_REVERSE_STRAND_MNV));
    both.name = "MaxEntScan_diff";
    assert_eq!(class(&both), Some(MINIMAL_PLUGIN_SEES_THE_REJOINED_RECORD));
}

#[test]
fn a_field_reaches_the_class_of_the_format_comparator_past_every_class_before_it() {
    let transcript = key("21:100", "G", "Transcript");
    let note = question("NOTE", &transcript, "x", "y");
    // Under models that hold the comparator's acceptance of the field, the class at the
    // place takes the field and no other class does.
    let mut models = FieldModels::default();
    let accepted = models.spliceai.accepted.entry(transcript.clone());
    accepted.or_default().insert("NOTE".to_string());
    for (place, asked) in FIELD_CLASSES.iter().enumerate() {
        let taken = Verdict::Class(SPLICEAI_MULTI_ENTRY_CELL, None);
        let expected = if place == ASKED_OF_THE_COMPARATOR {
            taken
        } else {
            Verdict::Pass
        };
        assert_eq!(asked(&note, &models), expected, "{place}");
    }
    assert_eq!(
        classify_field_divergence(&note, &models),
        Some((SPLICEAI_MULTI_ENTRY_CELL, None))
    );
    // A field every class before the place passes is a question, whatever the answers
    // held; one a class before it takes, or ends the dispatch on, is none.
    assert!(reaches_the_comparator(&note, &FieldModels::default()));
    assert!(reaches_the_comparator(&note, &models));
    let motif = key("21:100-101", "AG", "MotifFeature");
    let mut score = question("MOTIF_SCORE_CHANGE", &motif, "-0.092", "-0.136");
    score.strand = "-1";
    assert!(!reaches_the_comparator(&score, &models));
    let protein = question("DownstreamProtein", &transcript, "MK", "ML");
    assert_eq!(
        super::super::downstream::downstream_field_class(&protein, &models),
        Verdict::NoClass
    );
    assert!(!reaches_the_comparator(&protein, &models));
    assert_eq!(classify_field_divergence(&protein, &models), None);
}

#[test]
fn a_field_reaches_the_sequence_edge_class_where_every_other_class_passes() {
    assert_eq!(
        SEQUENCE_EDGES_ASKED_OF_THE_COMPARATOR,
        FIELD_CLASSES.len() - 1
    );
    let transcript = key("edge:500", "-", "Transcript");
    let offset = question("HGVS_OFFSET", &transcript, "-", "-2");
    // Under models that hold the comparator's acceptance of the field, the last class
    // takes the field and no other class does.
    let mut models = FieldModels::default();
    let accepted = models.sequence_edges.accepted.entry(transcript.clone());
    accepted.or_default().insert("HGVS_OFFSET".to_string());
    for (place, asked) in FIELD_CLASSES.iter().enumerate() {
        let taken = Verdict::Class(HGVS_SHIFT_SEQUENCE_EDGE_WINDOW, None);
        let expected = if place == SEQUENCE_EDGES_ASKED_OF_THE_COMPARATOR {
            taken
        } else {
            Verdict::Pass
        };
        assert_eq!(asked(&offset, &models), expected, "{place}");
    }
    assert_eq!(
        classify_field_divergence(&offset, &models),
        Some((HGVS_SHIFT_SEQUENCE_EDGE_WINDOW, None))
    );
    // A field every other class passes is a question, whatever the answers held.
    assert!(reaches_the_sequence_edge_class(
        &offset,
        &FieldModels::default()
    ));
    assert!(reaches_the_sequence_edge_class(&offset, &models));
    // A field the class of the score files takes is none: that class is asked first.
    let mut scored = FieldModels::default();
    let accepted = scored.spliceai.accepted.entry(transcript.clone());
    accepted.or_default().insert("HGVS_OFFSET".to_string());
    assert!(reaches_the_comparator(&offset, &scored));
    assert!(!reaches_the_sequence_edge_class(&offset, &scored));
    assert_eq!(
        classify_field_divergence(&offset, &scored),
        Some((SPLICEAI_MULTI_ENTRY_CELL, None))
    );
    // Nor is one a class before it takes, or ends the dispatch on.
    let motif = key("21:100-101", "AG", "MotifFeature");
    let mut score = question("MOTIF_SCORE_CHANGE", &motif, "-0.092", "-0.136");
    score.strand = "-1";
    assert!(!reaches_the_sequence_edge_class(&score, &models));
    let protein = question("DownstreamProtein", &transcript, "MK", "ML");
    assert!(!reaches_the_sequence_edge_class(&protein, &models));
}

#[test]
fn a_covered_plugins_field_on_a_split_record_is_the_rejoined_record_class() {
    for name in REJOINED_RECORD_PLUGIN_FIELDS {
        assert!(rejoined_record_plugin_field(name), "{name}");
    }
    for name in [
        "gnomAD_mean_cov",
        "GWAS_PMID",
        "AVI_x",
        "OpenTargets_l2g",
        "am_class",
        "BayesDel_score",
        "CADD_PHRED",
        "MaxEntScan_alt",
        "SpliceAI_pred",
        "SpliceAI_pred_DS_AG",
        "SpliceVault_top_events",
        "5UTR_annotation",
    ] {
        assert!(rejoined_record_plugin_field(name), "{name}");
    }
    for name in [
        "NOTE",
        "gnomAD",
        "REVEL_score",
        "revel",
        "SpliceAI_cutoff",
        "HGVSc",
        "",
    ] {
        assert!(!rejoined_record_plugin_field(name), "{name}");
    }
    let transcript = key("21:100", "CA", "Transcript");
    let mut covered = question("gnomAD_mean_cov", &transcript, "1.7938", "1.7937");
    assert_eq!(class(&covered), None);
    covered.minimal_split_record = true;
    assert_eq!(
        class(&covered),
        Some(MINIMAL_PLUGIN_SEES_THE_REJOINED_RECORD)
    );
    let mut other = question("NOTE", &transcript, "x", "y");
    other.minimal_split_record = true;
    assert_eq!(class(&other), None);
}

#[test]
fn a_motif_score_on_a_reverse_strand_row_of_several_bases_is_the_motif_class() {
    let taken = |location: &str, allele: &str, feature_type: &str, strand: &str, name: &str| {
        let key = key(location, allele, feature_type);
        let mut asked = question(name, &key, "-0.092", "-0.136");
        asked.strand = strand;
        class(&asked) == Some(MOTIF_SCORE_REVERSE_STRAND_MNV)
    };
    let score = "MOTIF_SCORE_CHANGE";
    assert!(taken("21:100-101", "AG", "MotifFeature", "-1", score));
    assert!(taken("21:100-102", "AGT", "MotifFeature", "-1", score));
    assert!(taken("100-101", "AG", "MotifFeature", "-1", score));
    // The forward strand, one base, an allele that is not as long as the span, an allele
    // that holds another letter, a span that is no pair of numbers, another feature type
    // and another field.
    assert!(!taken("chr21:5:100-101", "AG", "MotifFeature", "-1", score));
    assert!(!taken("21:100-101", "AG", "MotifFeature", "1", score));
    assert!(!taken("21:100-101", "AG", "MotifFeature", "-", score));
    assert!(!taken("21:300", "G", "MotifFeature", "-1", score));
    assert!(!taken("21:100-102", "AG", "MotifFeature", "-1", score));
    assert!(!taken("21:100-101", "AN", "MotifFeature", "-1", score));
    assert!(!taken("21:100-101", "ag", "MotifFeature", "-1", score));
    assert!(!taken("21:100", "AG", "MotifFeature", "-1", score));
    assert!(!taken("21:a-b", "AG", "MotifFeature", "-1", score));
    assert!(!taken("21:101-100", "AG", "MotifFeature", "-1", score));
    assert!(!taken("21:100-101", "AG", "RegulatoryFeature", "-1", score));
    assert!(!taken(
        "21:100-101",
        "AG",
        "MotifFeature",
        "-1",
        "MOTIF_POS"
    ));
    // Coordinates past 64 bits are still numbers.
    assert!(taken(
        "21:100000000000000000000-100000000000000000001",
        "AG",
        "MotifFeature",
        "-1",
        score
    ));
    // A span whose count of positions 128 bits do not hold is the span of no allele.
    let widest = "21:0-170141183460469231731687303715884105727";
    assert_eq!(span_positions(widest), None);
    assert!(!taken(widest, "AG", "MotifFeature", "-1", score));
    let below = "21:1-170141183460469231731687303715884105727";
    assert_eq!(span_positions(below), Some(i128::MAX));
}

#[test]
fn an_empty_cell_type_beside_an_activity_under_the_flag_is_the_comma_class() {
    let regulatory = key("21:400", "G", "RegulatoryFeature");
    let asked = |vep: &'static str, vep_rs: &'static str, flags: &[&str]| {
        let flags: Vec<String> = flags.iter().map(|flag| flag.to_string()).collect();
        let mut asked = question("CELL_TYPE", &regulatory, vep, vep_rs);
        asked.flags = &flags;
        class(&asked) == Some(CELL_TYPE_COMMA_NAMES)
    };
    let activity = "A549_(m,_58_y):ACTIVE";
    assert!(asked(
        "",
        activity,
        &["--regulatory", "--cell_type", "A549_(m,_58_y)"]
    ));
    assert!(asked("", activity, &["--cell_type=A549_(m,_58_y)"]));
    // Without the flag, with a reference value, and without a vep-rs value.
    assert!(!asked("", activity, &["--regulatory"]));
    assert!(!asked("", activity, &["--cell_types", "--cell_typ"]));
    assert!(!asked("-", activity, &["--cell_type"]));
    assert!(!asked("HepG2:ACTIVE", activity, &["--cell_type"]));
    assert!(!asked("", "-", &["--cell_type"]));
    let flags = ["--cell_type".to_string()];
    let mut other = question("MOTIF_NAME", &regulatory, "", activity);
    other.flags = &flags;
    assert_eq!(class(&other), None);
}

#[test]
fn a_codon_ensembl_vep_reads_into_the_utr_is_the_coding_end_class() {
    let transcript = key("21:5000", "A", "Transcript");
    let row = |position: &str, codons: &str| -> FieldRow {
        [("CDS_position", position), ("Codons", codons)]
            .iter()
            .map(|(name, value)| (name.to_string(), value.to_string()))
            .collect()
    };
    // A window past a coding end inside a codon: the codons, and the residues beside them.
    let (perl, rust) = (row("298", "Ta/Aag"), row("298", "Ta/Aa"));
    let asked = |name: &'static str, vep: &'static str, vep_rs: &'static str| {
        let mut asked = question(name, &transcript, vep, vep_rs);
        asked.perl_row = Some(&perl);
        asked.rust_row = Some(&rust);
        class(&asked)
    };
    assert_eq!(
        asked("Codons", "Ta/Aag", "Ta/Aa"),
        Some(ALTERNATE_CODON_PAST_THE_CODING_END)
    );
    assert_eq!(
        asked("Amino_acids", "X/K", "X"),
        Some(ALTERNATE_CODON_PAST_THE_CODING_END)
    );
    // Residues that are no translation of the row's codons, residues asked without the
    // rows' codons, the two sides exchanged and another field are no class.
    assert_eq!(asked("Amino_acids", "X/R", "X"), None);
    assert_eq!(
        class(&question("Amino_acids", &transcript, "X/K", "X")),
        None
    );
    assert_eq!(asked("Codons", "Ta/Aa", "Ta/Aag"), None);
    assert_eq!(asked("HGVSp", "Ta/Aag", "Ta/Aa"), None);
    // A deletion that leaves one coding base, where the row's coding positions open in the
    // first codon.
    let (perl, rust) = (row("2-5", "aTGTAa/gg"), row("2-5", "aTGTAa/aa"));
    let mut left = question("Codons", &transcript, "aTGTAa/gg", "aTGTAa/aa");
    left.perl_row = Some(&perl);
    left.rust_row = Some(&rust);
    assert_eq!(class(&left), Some(ALTERNATE_CODON_PAST_THE_CODING_END));
    let later = row("4-8", "aTGTAa/gg");
    left.perl_row = Some(&later);
    assert_eq!(class(&left), None);
}
