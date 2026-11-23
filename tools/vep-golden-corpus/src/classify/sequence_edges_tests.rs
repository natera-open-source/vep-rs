// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the sequence-edge class: the line a question is, the fields the class
//! takes of those the format comparator accepts, and the refusal of a corpus cache
//! without transcripts.

use super::*;
use crate::classify::fields::{FieldKey, FieldRow};
use crate::json::Object;
use crate::TempDir;
use std::path::PathBuf;

fn row(fields: &[(&str, &str)]) -> FieldRow {
    fields
        .iter()
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect()
}

fn key(feature: &str) -> FieldKey {
    FieldKey {
        location: "edge:500".to_string(),
        allele: "-".to_string(),
        feature: feature.to_string(),
        feature_type: "Transcript".to_string(),
        consequences: "non_coding_transcript_exon_variant".to_string(),
    }
}

fn question<'a>(
    name: &'a str,
    key: &'a FieldKey,
    values: (&'a str, &'a str),
    rows: (Option<&'a FieldRow>, Option<&'a FieldRow>),
) -> FieldQuestion<'a> {
    FieldQuestion {
        name,
        key,
        vep: values.0,
        vep_rs: values.1,
        strand: "-1",
        flags: &[],
        perl_row: rows.0,
        rust_row: rows.1,
        minimal_split_record: false,
        genomic_span: None,
    }
}

#[test]
fn a_question_is_a_line_of_twenty_four_columns_with_a_dash_for_a_field_a_row_lacks() {
    let key = key("TXR");
    // The run compares the notation and the offset and neither a protein notation nor a
    // shift length: each row holds the two beside its position columns, and the row of the
    // reference a strand vep-rs's row does not load.
    let perl_row = row(&[
        ("HGVSc", "TXR.2:n.401del"),
        ("HGVS_OFFSET", "-"),
        ("STRAND", "-1"),
        ("cDNA_position", "401"),
        ("CDS_position", "-"),
        ("Protein_position", "-"),
    ]);
    let rust_row = row(&[
        ("HGVSc", "TXR.2:n.403del"),
        ("HGVS_OFFSET", "-2"),
        ("cDNA_position", "401"),
        ("CDS_position", "-"),
        ("Protein_position", "-"),
    ]);
    let asked = question(
        "HGVS_OFFSET",
        &key,
        ("-", "-2"),
        (Some(&perl_row), Some(&rust_row)),
    );
    let line = question_line(&asked);
    assert_eq!(
        line,
        "edge:500\t-\tTXR\tHGVS_OFFSET\t-\t-2\
         \tTranscript\t-1\tTXR.2:n.401del\t-\t-\t-\t401\t-\t-\
         \tTranscript\t-\tTXR.2:n.403del\t-\t-2\t-\t401\t-\t-\n"
    );
    assert_eq!(
        line.split('\t').count(),
        6 + 2 * (2 + SEQUENCE_EDGE_ROW_FIELDS.len())
    );
    // A question whose rows are not at hand carries the key and the strand alone.
    let bare = question_line(&question("HGVSc", &key, ("a", "b"), (None, None)));
    assert_eq!(
        bare,
        "edge:500\t-\tTXR\tHGVSc\ta\tb\
         \tTranscript\t-1\t-\t-\t-\t-\t-\t-\t-\
         \tTranscript\t-\t-\t-\t-\t-\t-\t-\t-\n"
    );
}

#[test]
fn the_class_takes_the_fields_the_format_comparator_accepts() {
    let (first, second) = (key("TXR"), key("TXF"));
    let asked = |name: &'static str, key: &FieldKey, models: Models| {
        let asked = question(name, key, ("-", "-2"), (None, None));
        let models = FieldModels {
            sequence_edges: models,
            ..FieldModels::default()
        };
        hgvs_shift_sequence_edge_window(&asked, &models)
    };
    // A run without a reference sequence holds no answer and takes no field.
    assert_eq!(
        asked("HGVS_OFFSET", &first, Models::default()),
        Verdict::Pass
    );
    let accepting = || {
        let mut models = Models::default();
        let fields = models.accepted.entry(first.clone()).or_default();
        fields.insert("HGVS_OFFSET".to_string());
        models
    };
    assert_eq!(
        asked("HGVS_OFFSET", &first, accepting()),
        Verdict::Class(HGVS_SHIFT_SEQUENCE_EDGE_WINDOW, None)
    );
    // Another field of the key, and the field on another key, are not the comparator's.
    assert_eq!(asked("HGVSc", &first, accepting()), Verdict::Pass);
    assert_eq!(asked("HGVS_OFFSET", &second, accepting()), Verdict::Pass);
}

#[test]
fn a_corpus_cache_without_transcripts_is_refused_before_the_command_starts() {
    let scratch = TempDir::create(&std::env::temp_dir(), "golden_sequence_edges_").unwrap();
    let cache = scratch.path.join("json_cache");
    std::fs::create_dir(&cache).unwrap();
    let cache = cache.to_str().unwrap();
    let manifest = Object::new();
    let run = FieldRun {
        corpus: "c",
        cache_dir: cache,
        manifest: &manifest,
        fields: &["HGVSc".to_string()],
        spliceai_tables: &[],
        reference_fasta: Some("reference.fa"),
    };
    // The program named is no file: a run that reached it would end on that.
    let comparator = Comparator {
        program: PathBuf::from("no/such/vep-formats-compare"),
    };
    let Err(failure) = ask(&comparator, &run, "reference.fa", &[]) else {
        panic!("the run is refused");
    };
    assert_eq!(failure.code, 1);
    assert_eq!(
        failure.line,
        format!(
            "ERROR: [vep-golden-corpus] --reference-fasta: the corpus cache {cache} holds no \
             transcripts directory, which hgvs_shift_sequence_edge_window reads"
        )
    );
    std::fs::create_dir(scratch.path.join("json_cache/transcripts")).unwrap();
    let Err(failure) = ask(&comparator, &run, "reference.fa", &[]) else {
        panic!("the command does not start");
    };
    assert!(
        failure
            .line
            .starts_with("ERROR: [vep-golden-corpus] the format comparator no/such/vep-formats-compare did not start ("),
        "{}",
        failure.line
    );
}
