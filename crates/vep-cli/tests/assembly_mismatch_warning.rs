// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `Transcript-assembly mismatch` warning through the `vep` binary on the
//! GRCh38 RefSeq golden corpus. Ensembl VEP draws it inside the output hash of
//! each row it prints (OutputFactory.pm:1717-1718), after the row filters and
//! the pick family have dropped theirs (:477-514), and never under `--summary`
//! or `--most_severe`, which build no row hash (:387). `rs922296917`, a
//! nineteen-base coding deletion on COL18A1, marks its rows on the edited models
//! `NM_030582.4` and `NM_130444.3`, whose first `_rna_edit` deletes, while the
//! pick family keeps the unedited `NM_001379500.1`.

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The GRCh38 RefSeq golden corpus of release 116.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-refseq")
}

/// The record that marks its rows on the two edited models.
const RECORD: &str = "21\t45505371\trs922296917\tCGGCCCTCCGGGCCCCCCTG\tC\t.\t.\t.";

/// The warning the record draws.
const WARNING: &str = "WARNING: Transcript-assembly mismatch in rs922296917";

/// Runs the binary on the record under `--refseq`, `flags` and the pipeline
/// `mode`, and returns the models of the default-format rows and the lines of
/// the warnings file (none when no file was written). Each run writes its own
/// output, so a warnings file one run leaves is never read as another's.
fn annotate(dir: &Path, mode: &str, flags: &[&str]) -> (Vec<String>, Vec<String>) {
    let corpus = corpus();
    let input = dir.join("in.vcf");
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{RECORD}\n"),
    )
    .unwrap();
    let out = dir.join(format!("{mode}{}.txt", flags.concat()));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["--offline", "--cache", "--dir_cache"])
        .arg(corpus.join("cache_root"))
        .args(["--cache_version", "116", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--refseq", "--fasta"])
        .arg(corpus.join("reference.fa.gz"))
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .args(["--pipeline_mode", mode])
        .args(flags)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let models = fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').nth(4).unwrap_or_default().to_string())
        .collect();
    let warnings = fs::read_to_string(format!("{}_warnings.txt", out.display()))
        .unwrap_or_default()
        .lines()
        .map(str::to_string)
        .collect();
    (models, warnings)
}

/// The pipeline modes a run takes.
const MODES: [&str; 2] = ["batched", "streaming"];

/// The two coordinator arms: a run without a plugin annotates, filters and
/// renders each task's records in one pass, and a run with a plugin crosses the
/// pool once per phase; `TSSDistance` needs no data file and adds nothing to a
/// coding row, so it switches the arm and changes no asserted value.
const ARMS: [&[&str]; 2] = [&[], &["--plugin", "TSSDistance"]];

/// Every mode and arm with `flags` appended: the mode, a label and the run's
/// flags.
fn shapes<'a>(flags: &'a [&'a str]) -> Vec<(&'static str, String, Vec<&'a str>)> {
    let mut shapes = Vec::new();
    for mode in MODES {
        for arm in ARMS {
            let mut run: Vec<&str> = arm.to_vec();
            run.extend_from_slice(flags);
            shapes.push((mode, format!("{mode} {}", run.join(" ")), run));
        }
    }
    shapes
}

/// Without a record-level flag every row is printed, the two marked rows among
/// them, and the record warns once; `--flag_pick` and `--coding_only` keep the
/// marked rows and warn alike; `--use_given_ref` prints the rows and warns
/// nothing.
#[test]
fn a_printed_marked_row_warns_once() {
    let dir = tempfile::tempdir().unwrap();
    let all = [
        "NM_001379500.1",
        "NM_030582.4",
        "NM_130444.3",
        "XM_047440957.1",
        "XM_047440963.1",
        "XM_047440964.1",
    ];
    let coding = ["NM_001379500.1", "NM_030582.4", "NM_130444.3"];
    for (mode, label, run) in shapes(&[]) {
        let (models, warnings) = annotate(dir.path(), mode, &run);
        assert_eq!(models, all, "{label}");
        assert_eq!(warnings, [WARNING], "{label}");
    }
    for (mode, label, run) in shapes(&["--flag_pick"]) {
        let (models, warnings) = annotate(dir.path(), mode, &run);
        assert_eq!(models, all, "{label}");
        assert_eq!(warnings, [WARNING], "{label}");
    }
    for (mode, label, run) in shapes(&["--coding_only"]) {
        let (models, warnings) = annotate(dir.path(), mode, &run);
        assert_eq!(models, coding, "{label}");
        assert_eq!(warnings, [WARNING], "{label}");
    }
    for (mode, label, run) in shapes(&["--use_given_ref"]) {
        let (models, warnings) = annotate(dir.path(), mode, &run);
        assert_eq!(models, all, "{label}");
        assert!(warnings.is_empty(), "{label}: {warnings:?}");
    }
}

/// A removing pick keeps the unedited model alone, the marked rows are dropped
/// before any row hash is built, and the record warns nothing; a summary flag
/// builds no row hash and warns nothing either. Ensembl VEP writes no warnings
/// file for any of these runs.
#[test]
fn a_marked_row_the_record_step_drops_warns_nothing() {
    let dir = tempfile::tempdir().unwrap();
    for flag in ["--pick", "--pick_allele"] {
        for (mode, label, run) in shapes(&[flag]) {
            let (models, warnings) = annotate(dir.path(), mode, &run);
            assert_eq!(models, ["NM_001379500.1"], "{label}");
            assert!(warnings.is_empty(), "{label}: {warnings:?}");
        }
    }
    for flag in ["--per_gene", "--pick_allele_gene"] {
        for (mode, label, run) in shapes(&[flag]) {
            let (models, warnings) = annotate(dir.path(), mode, &run);
            assert_eq!(models.len(), 2, "{label}: {models:?}");
            assert!(
                models.contains(&"NM_001379500.1".to_string()),
                "{label}: {models:?}"
            );
            assert!(warnings.is_empty(), "{label}: {warnings:?}");
        }
    }
    for flag in ["--most_severe", "--summary"] {
        for (mode, label, run) in shapes(&[flag]) {
            let (models, warnings) = annotate(dir.path(), mode, &run);
            assert_eq!(models, ["-"], "{label}");
            assert!(warnings.is_empty(), "{label}: {warnings:?}");
        }
    }
}
