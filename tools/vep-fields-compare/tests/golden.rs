// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over the fixture pairs and compares `fields.json` and stdout byte for byte
//! with the goldens, which hold the tool's contract.

mod common;

use common::{assert_clean_exit, fixture, golden, run, run_args, run_in, stdout, Scratch, ROOT};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::ffi::OsStr;
use std::fs;
use std::io::Write;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

#[test]
fn fixture_pair_matches_its_goldens() {
    let scratch = Scratch::new("default");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/reference.txt",
        "tests/fixtures/scored.txt",
        &[],
    );
    assert_clean_exit(&output);
    assert_eq!(json, golden("expected_fields.json"));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        golden("expected_stdout.txt")
    );
    let leftovers: Vec<String> = fs::read_dir(&scratch.0)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(
        leftovers,
        ["report"],
        "the sort intermediates' directory is removed"
    );
}

#[test]
fn writes_report_and_exits_zero() {
    let scratch = Scratch::new("pair");
    let out = scratch.0.join("r.json");
    let output = run_args(&[
        "--perl",
        "tests/fixtures/pair_reference.txt",
        "--rust",
        "tests/fixtures/pair_scored.txt",
        "--out",
        out.to_str().unwrap(),
    ]);
    assert_clean_exit(&output);
    let json = fs::read_to_string(&out).unwrap();
    assert_eq!(json, golden("expected_pair_fields.json"));
    assert!(json.contains("\"intersection\": 4,"));
    let summary = stdout(&output);
    assert_eq!(summary, golden("expected_pair_stdout.txt"));
    assert!(summary.contains("Uploaded_variation") && summary.contains("Extra:DISTANCE"));
}

#[test]
fn examples_caps_the_mismatch_shapes_in_counter_order() {
    let scratch = Scratch::new("examples1");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/reference.txt",
        "tests/fixtures/scored.txt",
        &["--examples", "1"],
    );
    assert_clean_exit(&output);
    assert_eq!(json, golden("expected_fields_examples1.json"));
}

#[test]
fn disjoint_pair_reports_null_agreements() {
    let scratch = Scratch::new("disjoint");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/disjoint_reference.txt",
        "tests/fixtures/disjoint_scored.txt",
        &[],
    );
    assert_clean_exit(&output);
    assert_eq!(json, golden("expected_disjoint_fields.json"));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        golden("expected_disjoint_stdout.txt")
    );
}

#[test]
fn keep_only_compares_the_regulatory_and_motif_rows_alone() {
    let scratch = Scratch::new("keep-only");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/regulatory_reference.txt",
        "tests/fixtures/regulatory_scored.txt",
        &["--keep-only", "regulatory_feature_rows"],
    );
    assert_clean_exit(&output);
    assert_eq!(json, golden("expected_regulatory_fields.json"));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        golden("expected_regulatory_stdout.txt")
    );
}

#[test]
fn keep_only_filters_gzipped_inputs_alike() {
    let scratch = Scratch::new("keep-only-gz");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/regulatory_reference.txt.gz",
        "tests/fixtures/regulatory_scored.txt.gz",
        &["--keep-only=regulatory_feature_rows"],
    );
    assert_clean_exit(&output);
    assert_eq!(json, golden("expected_regulatory_gz_fields.json"));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        golden("expected_regulatory_gz_stdout.txt")
    );
}

#[test]
fn without_keep_only_regulatory_rows_are_compared_beside_every_other_row() {
    let scratch = Scratch::new("keep-only-absent");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/regulatory_reference.txt",
        "tests/fixtures/regulatory_scored.txt",
        &[],
    );
    assert_clean_exit(&output);
    assert_eq!(json, golden("expected_regulatory_whole_fields.json"));
    assert!(!json.contains("keep_only"));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        golden("expected_regulatory_whole_stdout.txt")
    );
}

#[test]
fn keep_only_refuses_a_class_it_does_not_name() {
    let scratch = Scratch::new("keep-only-bogus");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/regulatory_reference.txt",
        "tests/fixtures/regulatory_scored.txt",
        &["--keep-only", "transcript_rows"],
    );
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).ends_with(
        "vep-fields-compare: error: argument --keep-only: invalid choice: 'transcript_rows' (choose from 'regulatory_feature_rows')\n"
    ));
    assert!(output.stdout.is_empty());
    assert!(json.is_empty());
}

#[test]
fn gzip_input_is_read_through_its_suffix() {
    let scratch = Scratch::new("gz");
    let gz_path = scratch.0.join("reference.txt.gz");
    let mut encoder = GzEncoder::new(fs::File::create(&gz_path).unwrap(), Compression::fast());
    encoder
        .write_all(&fs::read(Path::new(ROOT).join("tests/fixtures/reference.txt")).unwrap())
        .unwrap();
    encoder.finish().unwrap();
    let gz_arg = gz_path.to_str().unwrap();
    let (output, json) = run(&scratch, gz_arg, "tests/fixtures/scored.txt", &[]);
    assert_clean_exit(&output);
    let expected = golden("expected_fields.json").replace(
        "\"reference\": \"tests/fixtures/reference.txt\"",
        &format!("\"reference\": \"{gz_arg}\""),
    );
    assert_eq!(json, expected);
}

#[test]
fn gzip_members_and_the_zero_padding_after_them_read_as_one_input() {
    let scratch = Scratch::new("gz-padded");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/pair_reference_padded.txt.gz",
        "tests/fixtures/pair_scored.txt",
        &[],
    );
    assert_clean_exit(&output);
    assert_eq!(json, golden("expected_pair_padded_gz_fields.json"));
    assert_eq!(stdout(&output), golden("expected_pair_stdout.txt"));
}

#[test]
fn an_empty_gzip_input_is_an_empty_side() {
    let scratch = Scratch::new("gz-empty");
    let gz_path = scratch.0.join("reference.txt.gz");
    fs::write(&gz_path, b"").unwrap();
    let (output, json) = run(
        &scratch,
        gz_path.to_str().unwrap(),
        "tests/fixtures/scored.txt",
        &[],
    );
    assert_clean_exit(&output);
    assert!(json.contains("\"reference_keys\": 0,\n    \"scored_keys\": 129,"));
    assert!(String::from_utf8_lossy(&output.stdout).starts_with(
        "keys: reference 0 scored 129 intersection 0 (reference-only 0, scored-only 129)\n"
    ));
}

#[test]
fn input_paths_are_normalised_as_pathlib_does() {
    let scratch = Scratch::new("normalised");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/reference.txt/",
        "./tests//fixtures/scored.txt/.",
        &[],
    );
    assert_clean_exit(&output);
    assert_eq!(json, golden("expected_fields.json"));
}

#[test]
fn input_names_outside_utf8_are_reported_as_surrogate_escapes() {
    let scratch = Scratch::new("undecodable-names");
    let reference = OsStr::from_bytes(b"ref\xff\xe2\x82.txt");
    let scored = OsStr::from_bytes(b"d\xc0\x80r/sco\xed\xa0\x80red.txt");
    if let Err(e) = fs::create_dir(scratch.0.join(OsStr::from_bytes(b"d\xc0\x80r"))) {
        // A file system that takes UTF-8 names alone cannot hold these inputs.
        eprintln!("not run: {e}");
        return;
    }
    fs::copy(fixture("pair_reference.txt"), scratch.0.join(reference)).unwrap();
    fs::copy(fixture("pair_scored.txt"), scratch.0.join(scored)).unwrap();
    let output = run_in(
        &scratch.0,
        &[
            OsStr::new("--perl"),
            reference,
            OsStr::new("--rust"),
            scored,
            OsStr::new("--out"),
            OsStr::new("fields.json"),
        ],
    );
    assert_clean_exit(&output);
    assert_eq!(
        fs::read_to_string(scratch.0.join("fields.json")).unwrap(),
        golden("expected_pair_undecodable_names_fields.json")
    );
    assert_eq!(stdout(&output), golden("expected_pair_stdout.txt"));
}

#[test]
fn a_missing_input_exits_2_with_one_line_that_names_it() {
    let scratch = Scratch::new("missing");
    let (output, json) = run(
        &scratch,
        "tests/fixtures/reference.txt",
        "./tests//fixtures/missing.txt",
        &[],
    );
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "ERROR: [compare_vep_fields] not a file: tests/fixtures/missing.txt\n"
    );
    assert!(output.stdout.is_empty());
    assert!(json.is_empty());
}
