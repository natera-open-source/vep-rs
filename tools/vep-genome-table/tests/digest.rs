// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `digest` over the two outputs under `tests/fixtures`, which hold the same records under
//! other `##` lines.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::{gzip, scratch, tool};
use vep_genome_table::digest_file;

/// `grep -v '^##' tests/fixtures/output_one.txt | sha256sum`, and the same for
/// `output_two.txt`.
const EXPECTED: &str = "e2ee2e106c3094584b3d23477cda36c6e63026844b3ec23145b15b31737a4b03";

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name)
}

fn digest(path: &Path) -> String {
    let run = tool(&["digest", path.to_str().unwrap()]);
    assert_eq!(run.stderr, "");
    assert_eq!(run.code, 0);
    run.stdout
}

#[test]
fn the_digest_is_the_sha256_of_the_lines_that_do_not_start_with_two_hashes() {
    assert_eq!(digest(&fixture("output_one.txt")), format!("{EXPECTED}\n"));
    assert_eq!(digest_file(&fixture("output_one.txt")).unwrap(), EXPECTED);
}

#[test]
fn the_same_output_under_other_meta_lines_has_the_same_digest() {
    let one = fs::read(fixture("output_one.txt")).unwrap();
    let two = fs::read(fixture("output_two.txt")).unwrap();
    assert_ne!(one, two);
    assert_eq!(digest(&fixture("output_two.txt")), format!("{EXPECTED}\n"));
}

#[test]
fn the_gzip_of_an_output_has_the_output_s_digest() {
    let dir = scratch("digest_gzip");
    let whole = dir.join("output_one.txt.gz");
    fs::write(&whole, gzip(&fs::read(fixture("output_one.txt")).unwrap())).unwrap();
    assert_eq!(digest(&whole), format!("{EXPECTED}\n"));
}

#[test]
fn a_gzip_of_several_members_is_read_to_its_end() {
    let bytes = fs::read(fixture("output_two.txt")).unwrap();
    let (head, tail) = bytes.split_at(bytes.len() / 2);
    let mut members = gzip(head);
    members.extend(gzip(tail));
    let path = scratch("digest_members").join("output_two.txt.gz");
    fs::write(&path, members).unwrap();
    assert_eq!(digest(&path), format!("{EXPECTED}\n"));
}

#[test]
fn a_changed_record_changes_the_digest() {
    let text = fs::read_to_string(fixture("output_one.txt")).unwrap();
    let path = scratch("digest_changed").join("output.txt");
    fs::write(&path, text.replace("intron_variant", "intergenic_variant")).unwrap();
    assert_ne!(digest(&path), format!("{EXPECTED}\n"));
}

#[test]
fn an_output_that_cannot_be_read_exits_1_with_one_error_line() {
    let path = scratch("digest_missing").join("absent.txt");
    let run = tool(&["digest", path.to_str().unwrap()]);
    let opening = format!("ERROR: [vep-genome-table] {}: ", path.display());
    assert!(run.stderr.starts_with(&opening), "{}", run.stderr);
    assert_eq!(run.stderr.lines().count(), 1);
    assert_eq!(run.stdout, "");
    assert_eq!(run.code, 1);
}
