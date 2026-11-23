// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over fixture pairs that each hold one shape of input (an empty file,
//! repeated keys, Extra keys on one side, ragged rows, line endings, odd bytes, tied
//! mismatch shapes, large counts) and compares `fields.json` and stdout byte for byte with
//! the goldens.

mod common;

use common::{assert_clean_exit, golden, run, run_in, stdout, Scratch};
use std::fmt::Write;
use std::fs;

/// Runs the binary over `tests/fixtures/<reference>` and `<scored>` with `options` and holds
/// the report against `expected_<name>_fields.json`.
fn assert_report(name: &str, reference: &str, scored: &str, options: &[&str]) -> String {
    let scratch = Scratch::new(name);
    let (output, json) = run(
        &scratch,
        &format!("tests/fixtures/{reference}"),
        &format!("tests/fixtures/{scored}"),
        options,
    );
    assert_clean_exit(&output);
    assert_eq!(
        json,
        golden(&format!("expected_{name}_fields.json")),
        "{name}"
    );
    stdout(&output)
}

/// The pair `<name>_reference.txt` and `<name>_scored.txt` against `expected_<name>_*`.
fn assert_pair(name: &str, options: &[&str]) {
    let summary = assert_report(
        name,
        &format!("{name}_reference.txt"),
        &format!("{name}_scored.txt"),
        options,
    );
    assert_eq!(
        summary,
        golden(&format!("expected_{name}_stdout.txt")),
        "{name}"
    );
}

#[test]
fn a_zero_byte_input_is_a_side_with_no_key() {
    let summary = assert_report("empty_reference", "empty.txt", "pair_scored.txt", &[]);
    assert_eq!(summary, golden("expected_empty_reference_stdout.txt"));
    let summary = assert_report("empty_both", "empty.txt", "empty.txt", &[]);
    assert_eq!(summary, golden("expected_empty_both_stdout.txt"));
}

#[test]
fn a_zero_byte_gzip_input_is_a_side_with_no_key() {
    let summary = assert_report("empty_gz", "pair_reference.txt", "empty.txt.gz", &[]);
    assert_eq!(summary, golden("expected_empty_gz_stdout.txt"));
}

#[test]
fn a_repeated_key_counts_once_and_keeps_its_first_row_in_bytewise_order() {
    assert_pair("duplicates", &[]);
}

#[test]
fn extra_keys_compare_one_by_one_and_a_key_on_one_side_compares_with_the_absent_marker() {
    assert_pair("extra", &[]);
}

#[test]
fn rows_under_seven_columns_are_dropped_and_missing_value_columns_read_empty() {
    assert_pair("ragged", &[]);
}

#[test]
fn records_end_at_lf_at_crlf_at_a_lone_cr_and_at_the_end_of_the_file() {
    assert_pair("line_endings", &[]);
}

#[test]
fn fields_are_stripped_and_bytes_outside_utf8_replaced_before_they_are_compared() {
    assert_pair("odd_bytes", &[]);
}

#[test]
fn a_byte_order_mark_makes_the_header_line_a_row() {
    assert_pair("byte_order_mark", &[]);
}

#[test]
fn a_consequence_set_that_ends_below_tab_is_merged_where_its_line_sorts() {
    assert_pair("below_tab", &[]);
}

#[test]
fn mismatch_shapes_are_listed_by_count_then_by_first_observation() {
    assert_pair("ties", &[]);
    assert_report(
        "ties_examples2",
        "ties_reference.txt",
        "ties_scored.txt",
        &["--examples", "2"],
    );
    assert_report(
        "ties_examples0",
        "ties_reference.txt",
        "ties_scored.txt",
        &["--examples", "0"],
    );
}

#[test]
fn examples_at_or_below_zero_keeps_no_mismatch_shape() {
    for count in ["0", "-1", "-99999999999999999999999999"] {
        let scratch = Scratch::new("examples0");
        let (output, json) = run(
            &scratch,
            "tests/fixtures/reference.txt",
            "tests/fixtures/scored.txt",
            &["--examples", count],
        );
        assert_clean_exit(&output);
        assert_eq!(
            json,
            golden("expected_fields_examples0.json"),
            "--examples {count}"
        );
    }
}

#[test]
fn keep_only_reads_the_stripped_feature_type_of_rows_of_seven_columns_or_more() {
    assert_pair(
        "keep_only_edges",
        &["--keep-only", "regulatory_feature_rows"],
    );
}

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

/// One side of a pair of 20,001 shared keys: the first column equal on one key, Gene on all
/// but one, cDNA_position on two of three.
fn generated_side(scored: bool) -> String {
    let mut out = String::from(HEADER);
    for n in 1..=20_001 {
        let uploaded = match (n, scored) {
            (1, _) => "same".to_string(),
            (_, false) => format!("r{n}"),
            (_, true) => format!("scored{n}"),
        };
        let gene = if scored && n == 7 { "ENSG2" } else { "ENSG1" };
        let cdna = if scored && n % 3 == 0 { "6" } else { "5" };
        writeln!(
            out,
            "{uploaded}\t1:{n}\tT\t{gene}\tENST1\tTranscript\tintron_variant\t{cdna}\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1"
        )
        .unwrap();
    }
    out
}

/// FNV-1a over the bytes, 64 bits.
fn fnv1a(data: &[u8]) -> u64 {
    data.iter().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(*byte)).wrapping_mul(0x0100_0000_01b3)
    })
}

#[test]
fn counts_past_a_thousand_are_grouped_and_agreements_round_to_six_places() {
    let scratch = Scratch::new("generated");
    let reference = generated_side(false);
    let scored = generated_side(true);
    // The goldens are over exactly these bytes.
    assert_eq!(
        (reference.len(), fnv1a(reference.as_bytes())),
        (1_818_044, 0x8355_7a39_29f7_63d6)
    );
    assert_eq!(
        (scored.len(), fnv1a(scored.as_bytes())),
        (1_918_044, 0x2fee_c7c5_2b7f_7eb8)
    );
    fs::write(scratch.0.join("reference.txt"), reference).unwrap();
    fs::write(scratch.0.join("scored.txt"), scored).unwrap();
    let output = run_in(
        &scratch.0,
        &[
            "--perl",
            "reference.txt",
            "--rust",
            "scored.txt",
            "--out",
            "fields.json",
            "--examples",
            "2",
        ],
    );
    assert_clean_exit(&output);
    assert_eq!(
        fs::read_to_string(scratch.0.join("fields.json")).unwrap(),
        golden("expected_generated_fields.json")
    );
    let summary = stdout(&output);
    assert_eq!(summary, golden("expected_generated_stdout.txt"));
    for line in [
        "keys: reference 20,001 scored 20,001 intersection 20,001 (reference-only 0, scored-only 0)\n",
        "  Uploaded_variation       equal            1 /       20,001  agreement 5e-05\n",
        "  cDNA_position            equal       13,334 /       20,001  agreement 0.666667\n",
    ] {
        assert!(summary.contains(line), "{line}");
    }
}
