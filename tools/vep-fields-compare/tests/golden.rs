// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over the fixture pairs and compares `fields.json` and stdout byte for byte
//! with the goldens `scripts/concordance/compare_vep_fields.py` wrote over the same files.

use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-fields-compare");
const ROOT: &str = env!("CARGO_MANIFEST_DIR");

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-fields-compare-golden-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn golden(name: &str) -> String {
    fs::read_to_string(Path::new(ROOT).join("tests/fixtures").join(name)).unwrap()
}

/// Runs the binary from the package root with the fixture paths relative to it, so the
/// report's `inputs` block carries the strings the golden holds.
fn run(scratch: &Scratch, reference: &str, scored: &str, extra: &[&str]) -> (Output, String) {
    let out = scratch.0.join("report").join("fields.json");
    let output = Command::new(BIN)
        .current_dir(ROOT)
        .args(["--perl", reference, "--rust", scored, "--out"])
        .arg(&out)
        .arg("--tmp-dir")
        .arg(&scratch.0)
        .args(extra)
        .output()
        .unwrap();
    let json = fs::read_to_string(&out).unwrap_or_default();
    (output, json)
}

fn assert_clean_exit(output: &Output) {
    assert!(
        output.status.success(),
        "exit {:?}: {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn fixture_pair_matches_the_python_goldens() {
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
fn a_missing_input_exits_2_with_the_python_message() {
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
