// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A byte-order mark at the start of a line changes nothing: each plain-file
//! reader of the builder builds a file with a mark before every one of its
//! lines to the tree of the same file without its marks. The readers are the
//! GFF3 model, the GTF model, the GTF read for its tags beside a GFF3, the
//! regulatory bake and the protein FASTA.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

/// A path under this crate's directory.
fn crate_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

/// The text of a plain or gzipped file.
fn text(path: &Path) -> String {
    let bytes = std::fs::read(path).unwrap();
    if path.extension().is_some_and(|e| e == "gz") {
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(&bytes[..])
            .read_to_string(&mut text)
            .unwrap();
        text
    } else {
        String::from_utf8(bytes).unwrap()
    }
}

/// Every file under `dir`, by its path within it.
fn tree(dir: &Path) -> BTreeMap<PathBuf, Vec<u8>> {
    let mut files = BTreeMap::new();
    let mut pending = vec![dir.to_path_buf()];
    while let Some(next) = pending.pop() {
        for entry in std::fs::read_dir(&next).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                pending.push(path);
            } else {
                let name = path.strip_prefix(dir).unwrap().to_path_buf();
                files.insert(name, std::fs::read(&path).unwrap());
            }
        }
    }
    files
}

/// Builds a cache twice, `flag` naming a plain copy of `source` and then the
/// same lines with a byte-order mark before each, `others` given to both, and
/// holds the marked build's tree to the plain build's, file by file.
fn assert_marks_change_nothing(flag: &str, source: &Path, others: &[(&str, &Path)]) {
    let dir = tempfile::tempdir().unwrap();
    let plain = text(source);
    let marked: String = plain.lines().map(|l| format!("\u{feff}{l}\n")).collect();
    assert_eq!(marked.matches('\u{feff}').count(), plain.lines().count());
    let build = |name: &str, lines: &str| {
        let input = dir.path().join(format!("{name}.txt"));
        std::fs::write(&input, lines).unwrap();
        let out = dir.path().join(name);
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"));
        cmd.args(["--release", "116", "--assembly", "GRCh38"]);
        cmd.arg(flag).arg(&input).arg("--output-dir").arg(&out);
        for (other, path) in others {
            cmd.arg(other).arg(path);
        }
        let output = cmd.output().expect("vep-cache-builder runs");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{name}: {stderr}");
        tree(&out)
    };
    let (plain, marked) = (build("plain", &plain), build("marked", &marked));
    assert!(plain.keys().any(|file| file.starts_with("transcripts")));
    let names = |tree: &BTreeMap<PathBuf, Vec<u8>>| tree.keys().cloned().collect::<Vec<_>>();
    assert_eq!(names(&marked), names(&plain));
    for (file, bytes) in &plain {
        assert!(marked[file] == *bytes, "{} differs", file.display());
    }
}

#[test]
fn a_marked_gff3_model_builds_the_tree_of_the_file_without_its_marks() {
    assert_marks_change_nothing(
        "--gff3",
        &crate_path("tests/fixtures/grch38_chr21/slice.gff3"),
        &[(
            "--fasta",
            &crate_path("tests/fixtures/grch38_chr21/peptides.fa"),
        )],
    );
}

/// The model's five `#!` header lines are marked with its records, so the
/// keys `info.json` takes from them are read behind a mark.
#[test]
fn a_marked_gtf_model_builds_the_tree_of_the_file_without_its_marks() {
    let corpus = crate_path("../../tests/golden/116/GRCh38-gtf-builder");
    assert_marks_change_nothing(
        "--gtf",
        &corpus.join("chr21_genes.gtf.gz"),
        &[("--genome-fasta", &corpus.join("reference.fa.gz"))],
    );
}

/// A GTF beside a GFF3 is read for its transcript tags and its header: both
/// are read behind a mark.
#[test]
fn a_marked_gtf_read_for_its_tags_builds_the_tree_of_the_file_without_its_marks() {
    let fixtures = crate_path("tests/fixtures/grch38_chr21");
    assert_marks_change_nothing(
        "--gtf",
        &fixtures.join("slice.gtf"),
        &[
            ("--gff3", &fixtures.join("slice.gff3")),
            ("--fasta", &fixtures.join("peptides.fa")),
        ],
    );
}

/// A protein FASTA with a mark before each header and each sequence line
/// loads the proteins of the same file without its marks.
#[test]
fn a_marked_protein_fasta_builds_the_tree_of_the_file_without_its_marks() {
    let fixtures = crate_path("tests/fixtures/grch38_chr21");
    assert_marks_change_nothing(
        "--fasta",
        &fixtures.join("peptides.fa"),
        &[("--gff3", &fixtures.join("slice.gff3"))],
    );
}

#[test]
fn a_marked_regulatory_file_bakes_the_tree_of_the_file_without_its_marks() {
    let fixtures = crate_path("tests/fixtures/grch38_chr21");
    assert_marks_change_nothing(
        "--regulatory-gff",
        &crate_path("../vep-cli/tests/fixtures/regulatory-gff/regulatory_features.gff3.gz"),
        &[
            ("--gff3", &fixtures.join("slice.gff3")),
            ("--fasta", &fixtures.join("peptides.fa")),
        ],
    );
}
