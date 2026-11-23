// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The record that names a rewritten file: the directory a path stands in, the nearest
//! `provenance.json` no higher than the named directory, and the digests recorded again in
//! the order they changed.

use super::*;
use crate::scratch::TempDir;
use std::fs;

#[test]
fn the_parent_of_a_path_is_the_one_pathlib_gives() {
    let cases = [
        ("a/b/c", "a/b"),
        ("a/b", "a"),
        ("a", "."),
        (".", "."),
        ("/a/b", "/a"),
        ("/a", "/"),
        ("/", "/"),
        ("//a", "//"),
    ];
    for (path, expected) in cases {
        assert_eq!(parent(path), expected, "{path}");
    }
}

#[test]
fn the_nearest_record_at_or_below_the_named_directory_records_a_file() {
    let scratch = TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-record-").unwrap();
    let above = scratch.path.to_str().unwrap().to_string();
    let root = format!("{above}/corpora");
    fs::create_dir_all(format!("{root}/one/expected")).unwrap();
    fs::create_dir_all(format!("{root}/two/expected")).unwrap();
    let file = format!("{root}/one/expected/vcf.vcf.gz");
    assert_eq!(recording_provenance(&file, &root), None);
    // A record above the named directory is not read.
    fs::write(format!("{above}/provenance.json"), "{}").unwrap();
    assert_eq!(recording_provenance(&file, &root), None);
    fs::write(format!("{root}/provenance.json"), "{}").unwrap();
    assert_eq!(
        recording_provenance(&file, &root),
        Some(format!("{root}/provenance.json"))
    );
    fs::write(format!("{root}/one/provenance.json"), "{}").unwrap();
    assert_eq!(
        recording_provenance(&file, &root),
        Some(format!("{root}/one/provenance.json"))
    );
    assert_eq!(
        recording_provenance(&format!("{root}/two/expected/vcf.vcf.gz"), &root),
        Some(format!("{root}/provenance.json"))
    );
    // A directory of that name is no record.
    fs::create_dir(format!("{root}/two/expected/provenance.json")).unwrap();
    assert_eq!(
        recording_provenance(&format!("{root}/two/expected/vcf.vcf.gz"), &root),
        Some(format!("{root}/provenance.json"))
    );
}

#[test]
fn a_digest_named_twice_keeps_its_first_place_and_its_last_value() {
    let mut changed = Vec::new();
    for (before, after) in [("a", "1"), ("b", "2"), ("a", "3")] {
        record_digest(&mut changed, before.to_string(), after.to_string());
    }
    let expected = [("a", "3"), ("b", "2")].map(|(a, b)| (a.to_string(), b.to_string()));
    assert_eq!(changed, expected);
}
