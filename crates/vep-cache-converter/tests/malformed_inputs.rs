// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A flat cache whose region and variation files are named for the parent
//! directory (`.._1-1000000.json`, `.._all_vars.gz` and its index): the
//! converter skips them with its unexpected-name warning and writes nothing
//! outside the output directory, whose entries are `info.json`, the
//! `transcripts/` directory and the checksum manifest written over them.

use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/parent_named_shards")
}

/// A scratch directory under the system temporary directory, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "vep-cache-converter-malformed-{}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The names of every entry directly under `dir`, sorted.
fn entries(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

#[test]
fn files_named_for_the_parent_directory_are_skipped_and_write_nothing_beside_the_output() {
    let scratch = Scratch::new();
    let out = scratch.0.join("out");
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .arg("--json-dir")
        .arg(fixture())
        .arg("--output-dir")
        .arg(&out)
        .output()
        .expect("vep-cache-converter runs");
    let err = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{err}");
    assert_eq!(
        err.matches("Skipping file with unexpected name").count(),
        3,
        "{err}"
    );
    assert_eq!(entries(&scratch.0), ["out"]);
    assert_eq!(
        entries(&out),
        ["CHECKSUMS.sha256", "info.json", "transcripts"]
    );
    assert_eq!(entries(&out.join("transcripts")), ["21"]);
    assert_eq!(entries(&out.join("transcripts/21")), ["1-1000000.json"]);
}
