// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--species` through the `vep` binary on the mouse corpus cache: the value
//! names the species directory under the cache root, so the Latin name finds
//! the cache and an alias names a directory the root does not hold.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn mouse_corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCm39-mouse")
}

/// The number of corpus records the runs annotate, from the top of its input.
const RECORDS: usize = 5;

/// Writes `in.vcf` under `dir`: the first [`RECORDS`] records of the corpus.
fn write_input(dir: &Path) {
    let corpus = fs::read_to_string(mouse_corpus().join("variants.vcf")).unwrap();
    let records: Vec<&str> = corpus
        .lines()
        .filter(|line| !line.starts_with('#'))
        .take(RECORDS)
        .collect();
    assert_eq!(records.len(), RECORDS, "the corpus holds the records");
    fs::write(
        dir.join("in.vcf"),
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}\n",
            records.join("\n")
        ),
    )
    .unwrap();
}

/// A cache root holding the corpus cache as `<root>/mus_musculus/116_GRCm39`
/// alone.
fn cache_root(dir: &Path) -> PathBuf {
    let root = dir.join("cache_root");
    copy_dir(
        &mouse_corpus().join("json_cache"),
        &root.join("mus_musculus").join("116_GRCm39"),
    );
    root
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

/// Runs the binary on `in.vcf` under `dir` against `root` in Ensembl VEP's
/// layout, writing `<dir>/<species>.txt`. Nothing is asserted about the exit
/// status.
fn run(dir: &Path, root: &Path, species: &str) -> Output {
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(dir.join("in.vcf"))
        .arg("-o")
        .arg(dir.join(format!("{species}.txt")))
        .args(["--cache", "--dir_cache"])
        .arg(root)
        .args(["--cache_version", "116", "--species", species])
        .arg("--force_overwrite")
        .stdin(Stdio::null());
    cmd.output().unwrap()
}

/// The Latin name annotates against the species directory, and a common name
/// or an abbreviation of the same species is the cache error naming the
/// directory the run looked for, with no output written.
///
/// Documented: "Species: names the cache directory under `--dir_cache`; the
/// Latin name only, with no alias such as human" (Cache, species and runtime,
/// Species and assembly). Ensembl VEP names the directory by the flag's value
/// as given (`CacheDir.pm:309-315`) and dies when the root lacks it (`:318`);
/// it replaces an alias by the Latin name only where it sets up a database
/// connection (`BaseRunner.pm:146-153`), which `--offline` skips (`:116`).
#[test]
fn a_species_alias_names_a_directory_the_cache_root_lacks() {
    let dir = tempfile::tempdir().unwrap();
    write_input(dir.path());
    let root = cache_root(dir.path());

    let output = run(dir.path(), &root, "mus_musculus");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "mus_musculus: {stderr}");
    let text = fs::read_to_string(dir.path().join("mus_musculus.txt")).unwrap();
    let leaf = root.join("mus_musculus").join("116_GRCm39");
    assert!(
        text.lines()
            .any(|line| line == format!("## Using cache in {}", leaf.display())),
        "mus_musculus: the header names the species directory:\n{text}"
    );
    let annotated: BTreeSet<&str> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.split('\t').next())
        .collect();
    assert_eq!(
        annotated.len(),
        RECORDS,
        "mus_musculus: one row per record at least:\n{text}"
    );

    for alias in ["mouse", "mmusculus"] {
        let output = run(dir.path(), &root, alias);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{alias}: {stderr}");
        assert!(
            stderr.lines().any(|line| {
                line == format!(
                    "ERROR: [cache] Cache directory {} not found",
                    root.join(alias).display()
                )
            }),
            "{alias}: {stderr}"
        );
        assert!(
            !dir.path().join(format!("{alias}.txt")).exists(),
            "{alias}: no output is written"
        );
    }
}
