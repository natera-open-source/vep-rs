// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Flags accepted without effect, through the `vep` binary: a run with the flag
//! writes the same default and JSON output as a run without it, apart from the
//! command-line header that records the flag; and the input order is never
//! checked.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

/// Sequence records of the release-116 corpus on both of its chromosomes.
const RECORDS: &str = "\
21\t30282169\tstop_ins4_twostops_rev\tC\tCTATT\t.\t.\t.
21\t30348409\tstop_del2_lost_rev\tCTA\tC\t.\t.\t.
21\t30430756\tstop_ins5_newstop_fwd\tT\tTCCTAA\t.\t.\t.
21\t34716478\tpenult_ins3_stop_fwd\tA\tAATA\t.\t.\t.
22\t20061573\tpartial_codon_del3_fwd\tGAAG\tG\t.\t.\t.
";

fn write_input(dir: &Path, records: &str) -> PathBuf {
    let path = dir.join("in.vcf");
    std::fs::write(
        &path,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"),
    )
    .unwrap();
    path
}

/// Runs the binary on `input` in `format` (`default` or `json`) with `extra`
/// flags and returns the output text and standard error.
fn annotate(dir: &Path, input: &Path, format: &str, extra: &[&str]) -> (String, String) {
    let out = dir.join(format!(
        "out-{format}-{}.txt",
        extra.join("_").replace('/', "_")
    ));
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache())
        .args(["--assembly", "GRCh38", "--force_overwrite", "--no_stats"])
        .args(extra)
        .stdin(Stdio::null());
    if format == "json" {
        cmd.arg("--json");
    }
    let output = cmd.output().unwrap();
    assert!(
        output.status.success(),
        "{extra:?} ({format}): {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (
        std::fs::read_to_string(&out).unwrap(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// The output without the header lines that differ between any two runs: the
/// command line, which names the flag under test, and the time of the run.
fn without_volatile_headers(text: &str) -> String {
    text.lines()
        .filter(|l| !l.starts_with("## VEP command-line:"))
        .filter(|l| !l.starts_with("## Output produced at "))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Asserts the flag changes nothing in the default and JSON outputs, the
/// annotation itself included, against a run without it.
fn assert_no_effect(flags: &[&str]) {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(dir.path(), RECORDS);
    for format in ["default", "json"] {
        let (without, _) = annotate(dir.path(), &input, format, &[]);
        let (with, stderr) = annotate(dir.path(), &input, format, flags);
        assert!(
            without.lines().filter(|l| !l.starts_with('#')).count() >= 5,
            "{format}: the run annotates"
        );
        assert_eq!(
            without_volatile_headers(&with),
            without_volatile_headers(&without),
            "{flags:?} changes the {format} output"
        );
        assert!(!stderr.contains("ERROR"), "{flags:?} ({format}): {stderr}");
    }
}

#[test]
fn condel_is_accepted_without_effect() {
    assert_no_effect(&["--condel", "b"]);
}

#[test]
fn hgnc_is_accepted_without_effect() {
    assert_no_effect(&["--hgnc"]);
}

#[test]
fn keep_ann_is_accepted_without_effect() {
    assert_no_effect(&["--keep_ann"]);
}

#[test]
fn phased_is_accepted_without_effect() {
    assert_no_effect(&["--phased"]);
}

#[test]
fn sereal_is_accepted_without_effect() {
    assert_no_effect(&["--sereal"]);
}

#[test]
fn solr_is_accepted_without_effect() {
    assert_no_effect(&["--solr"]);
}

#[test]
fn tmpdir_is_accepted_without_effect() {
    assert_no_effect(&["--tmpdir", "/var/tmp"]);
}

#[test]
fn use_given_ref_is_accepted_without_effect() {
    assert_no_effect(&["--use_given_ref"]);
}

#[test]
fn no_check_variants_order_is_accepted_without_effect() {
    assert_no_effect(&["--no_check_variants_order"]);
}

/// vep-rs always runs offline: `--offline` changes nothing (in Ensembl VEP it
/// switches on `--cache` and off `--database`, Config.pm:464-469).
#[test]
fn offline_changes_no_output() {
    assert_no_effect(&["--offline"]);
    assert_no_effect(&["--offline", "--cache"]);
}

/// An IUPAC ambiguity code in an allele is annotated, not skipped: Ensembl
/// VEP's allele check (Parser.pm:556-560) is unanchored, so `A/R` passes on its
/// `A` and only a string with no `ACGTN-` run at all is skipped; vep-rs applies
/// the same test. The record gets rows under the allele `R` and no warning.
#[test]
fn iupac_allele_is_annotated_not_skipped() {
    let dir = tempfile::tempdir().unwrap();
    let input = write_input(
        dir.path(),
        "21\t30348411\tiupac_r\tA\tR\t.\t.\t.\n21\t30348411\tno_run\tX\tY\t.\t.\t.\n",
    );
    let (text, stderr) = annotate(dir.path(), &input, "default", &[]);
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(!rows.is_empty(), "the IUPAC record annotates:\n{text}");
    for row in &rows {
        let columns: Vec<&str> = row.split('\t').collect();
        assert_eq!(columns[0], "iupac_r", "{row}");
        assert_eq!(columns[2], "R", "{row}");
    }
    assert!(!stderr.contains("iupac_r"), "{stderr}");
    assert!(
        stderr.contains("Invalid allele string X/Y or possible parsing error"),
        "{stderr}"
    );
}

/// An unsorted input is annotated in input order with no order check: Ensembl
/// VEP's `InputBuffer.pm:189-258` aborts after 100 records more than 100 bp
/// behind their predecessor, vep-rs never does, with or without
/// `--no_check_variants_order`.
#[test]
fn unsorted_input_is_annotated_in_input_order() {
    let dir = tempfile::tempdir().unwrap();
    let unsorted = RECORDS.lines().rev().collect::<Vec<_>>().join("\n") + "\n";
    let input = write_input(dir.path(), &unsorted);
    for flags in [&[][..], &["--no_check_variants_order"][..]] {
        let (text, stderr) = annotate(dir.path(), &input, "default", flags);
        let mut names: Vec<&str> = text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        names.dedup();
        let expected: Vec<&str> = unsorted
            .lines()
            .map(|l| l.split('\t').nth(2).unwrap())
            .collect();
        assert_eq!(names, expected, "{flags:?}: rows follow the input order");
        assert!(
            !stderr.to_lowercase().contains("unsorted") && !stderr.contains("ERROR"),
            "{flags:?}: {stderr}"
        );
    }
}
