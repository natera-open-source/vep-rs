// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--include-predictions` serves `homo_sapiens` on GRCh37 and GRCh38 alone.
//! A build or an installation of any other species or assembly is refused
//! with one `ERROR: [cache-builder]` line that names the limit and the two
//! routes that serve every species, an installation without the flag and the
//! conversion, before any input is resolved: nothing is downloaded, no site
//! is read and no directory is written.

use std::path::Path;
use std::process::Command;

const REFUSAL: &str = "ERROR: [cache-builder] --include-predictions fetches the SIFT and \
     PolyPhen matrices of homo_sapiens on GRCh37 or GRCh38 alone, not of mus_musculus on \
     GRCm39; --auto without the flag reads the matrices of any species from the release's \
     table dumps, and the conversion carries those of its Ensembl cache \
     (docs/cache-conversion.md)";

/// Runs the builder with `args` and asserts the refusal: exit 1, the one error
/// line, no download begun and nothing left under `scratch`.
fn assert_refused(scratch: &Path, args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(args)
        .output()
        .expect("vep-cache-builder runs");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{args:?}: {stderr}");
    let errors: Vec<&str> = stderr
        .lines()
        .filter(|line| line.starts_with("ERROR: ["))
        .collect();
    assert_eq!(errors, [REFUSAL], "{args:?}: {stderr}");
    assert!(!stderr.contains("Downloading"), "{args:?}: {stderr}");
    assert_eq!(
        std::fs::read_dir(scratch).unwrap().count(),
        0,
        "{args:?} wrote under {}",
        scratch.display()
    );
}

/// A build named one flag at a time, with no local input, would download the
/// species' GFF3 first; an installation would read the site's listings first.
/// Each is refused before either, the installation against a site directory
/// that holds no release.
#[test]
fn predictions_of_another_species_are_refused_before_any_input_is_resolved() {
    let scratch = tempfile::tempdir().unwrap();
    let cache = scratch.path().join("cache");
    assert_refused(
        scratch.path(),
        &[
            "--release",
            "116",
            "--assembly",
            "GRCm39",
            "--species",
            "mus_musculus",
            "--include-predictions",
            "--output-dir",
            cache.to_str().unwrap(),
        ],
    );

    let site = tempfile::tempdir().unwrap();
    assert_refused(
        scratch.path(),
        &[
            "--auto",
            "mus_musculus",
            "GRCm39",
            "--release",
            "116",
            "--include-predictions",
            "--dir_cache",
            cache.to_str().unwrap(),
            "--site",
            site.path().to_str().unwrap(),
        ],
    );
}
