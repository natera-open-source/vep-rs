// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Two readings of `--transcript_filter` the `vep` binary shares with Ensembl
//! VEP, run through the binary: the refusal of an `in` list file that exists
//! and cannot be read (`FilterSet.pm:397` at release/116.2, `open IN, $list or
//! throw("ERROR: Could not read from file $list\n")`), on a path every user is
//! refused, root included; and the one root under which every filter string and
//! the predicate `--exclude_predicted` appends are joined, each unit with its
//! own `and` or `or` and no parentheses around a string (`parse_filters`,
//! `FilterSet.pm:217-221`, every string's units pushed under the root;
//! `Config.pm:566-571` and `:870-871`, the predicate pushed after the written
//! strings).

use std::fs;
use std::os::unix::net::UnixListener;
use std::path::{Path, PathBuf};

mod common;

/// A release 116 golden corpus.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/golden/116/{name}"))
}

/// The RefSeq corpus record within reach of five models of two genes:
/// `NM_003225.3` of `TFF1`, protein coding, and four predicted `XR_` models of
/// `LOC105372815`, lncRNA.
const TFF1: &str = "21\t42362510\t21:42362510:G:T\tG\tT\t.\t.\t.";

const FIVE: [&str; 5] = [
    "NM_003225.3",
    "XR_007067875.1",
    "XR_007067876.1",
    "XR_007067877.1",
    "XR_937755.3",
];

/// A filter string with an `or` and no parentheses, true of every one of the
/// five models.
const CODING_OR_LNCRNA: &str = "biotype is protein_coding or biotype is lncRNA";

/// The predicate `--exclude_predicted` appends, written as a filter string of
/// its own.
const NOT_PREDICTED: &str = "not stable_id match X";

/// The GRCh38 RefSeq corpus cache under `--refseq`, with its reference.
fn refseq_cache() -> Vec<String> {
    let corpus = corpus("GRCh38-refseq");
    vec![
        "--cache".to_string(),
        "--dir_cache".to_string(),
        corpus.join("cache_root").display().to_string(),
        "--cache_version".to_string(),
        "116".to_string(),
        "--refseq".to_string(),
        "--fasta".to_string(),
        corpus.join("reference.fa.gz").display().to_string(),
    ]
}

/// Runs the binary on the RefSeq corpus record `TFF1` with `flags` and returns
/// the `Feature` column of its default-format rows, `-` on an intergenic row.
fn features(dir: &Path, flags: &[&str]) -> Vec<String> {
    let input = dir.join("in.vcf");
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{TFF1}\n"),
    )
    .unwrap();
    let out = dir.join("out.txt");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args([
            "--offline",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
        ])
        .args(refseq_cache())
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .args(flags)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{flags:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').nth(4).unwrap().to_string())
        .collect()
}

/// Runs the binary with `filter` on no input and returns its exit status and
/// its standard error; a filter Ensembl VEP would refuse ends the run before
/// the input is opened.
fn refusal(dir: &Path, filter: &str) -> (Option<i32>, String) {
    let json_cache = corpus("GRCh38-transcript-filter").join("json_cache");
    let output = common::isolated_vep(dir)
        .args(["--transcript_filter", filter])
        .arg("--json_cache")
        .arg(&json_cache)
        .args(["-i", "/nonexistent/file.vcf", "--quiet"])
        .output()
        .unwrap();
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stderr)
            .trim_end()
            .to_string(),
    )
}

/// An `in` list file that exists and cannot be read is refused with Ensembl
/// VEP's message for an `open` that fails, and a list path that is no file
/// with its message for a path `-e` denies (`FilterSet.pm:397` and `:415`).
/// Both fixtures refuse every user, root included, since neither failure is a
/// permission check: `open(2)` on a Unix domain socket fails with `ENXIO`, and
/// `stat(2)` or `open(2)` on a path through a regular file with `ENOTDIR`.
/// Ensembl VEP 116.2 refuses both paths with the same two messages.
#[test]
fn an_in_list_file_that_cannot_be_read_is_refused_under_every_user() {
    let dir = tempfile::tempdir().unwrap();
    let socket = dir.path().join("list.sock");
    let _listener = UnixListener::bind(&socket).unwrap();
    assert!(socket.exists());
    let (status, stderr) = refusal(dir.path(), &format!("stable_id in {}", socket.display()));
    assert_eq!(status, Some(1));
    assert_eq!(
        stderr,
        format!(
            "ERROR: [filter] Could not read from file {}",
            socket.display()
        )
    );

    let file = dir.path().join("ids.txt");
    fs::write(&file, "NM_003225.3\n").unwrap();
    let through_file = file.join("x");
    assert!(!through_file.exists());
    let (status, stderr) = refusal(
        dir.path(),
        &format!("stable_id in {}", through_file.display()),
    );
    assert_eq!(status, Some(1));
    assert_eq!(
        stderr,
        format!(
            "ERROR: [filter] Could not find/parse list {}",
            through_file.display()
        )
    );
}

/// Every `--transcript_filter` string and the predicate `--exclude_predicted`
/// appends are units under one root, each joined to the units before it by its
/// own `and` or `or`, with no parentheses around a string, and the appended
/// predicate stands last whatever the flags' order. So `--exclude_predicted`
/// with the `or` filter keeps `NM_003225.3` alone, the predicate's `and`
/// applying to the whole `or` ((coding or lncRNA) and not predicted), while the
/// same predicate written first and the `or` filter second keeps all five, the
/// second string's `or` unit joining at the root ((not predicted and coding) or
/// lncRNA). Ensembl VEP 116.2 keeps the same models under each of the five
/// command lines.
#[test]
fn exclude_predicted_and_the_filter_strings_share_one_root() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(features(dir.path(), &[]), FIVE);
    assert_eq!(
        features(dir.path(), &["--transcript_filter", CODING_OR_LNCRNA]),
        FIVE
    );
    for flags in [
        &[
            "--exclude_predicted",
            "--transcript_filter",
            CODING_OR_LNCRNA,
        ][..],
        &[
            "--transcript_filter",
            CODING_OR_LNCRNA,
            "--exclude_predicted",
        ][..],
    ] {
        assert_eq!(features(dir.path(), flags), ["NM_003225.3"], "{flags:?}");
    }
    assert_eq!(
        features(
            dir.path(),
            &[
                "--transcript_filter",
                CODING_OR_LNCRNA,
                "--transcript_filter",
                NOT_PREDICTED,
            ]
        ),
        ["NM_003225.3"]
    );
    assert_eq!(
        features(
            dir.path(),
            &[
                "--transcript_filter",
                NOT_PREDICTED,
                "--transcript_filter",
                CODING_OR_LNCRNA,
            ]
        ),
        FIVE
    );
}
