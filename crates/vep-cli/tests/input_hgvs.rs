// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads HGVS genomic notation: the format needs the
//! reference FASTA, each notation annotates as its VCF twin does and is named
//! by itself, a notation Ensembl VEP cannot resolve is skipped with its
//! warning, and `--max_sv_size` bounds the span a notation may describe.

use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-input-hgvs")
}

/// The binary's arguments after the input: the corpus cache, the output path
/// and `extra`, with the corpus reference under `--fasta` when `fasta` is set.
fn args(out: &Path, fasta: bool, extra: &[&str]) -> Vec<String> {
    let mut args: Vec<String> = vec![
        "-o".into(),
        out.display().to_string(),
        "--json_cache".into(),
        corpus().join("json_cache").display().to_string(),
        "--assembly".into(),
        "GRCh38".into(),
        "--force_overwrite".into(),
        "--no_stats".into(),
        "--quiet".into(),
    ];
    if fasta {
        args.push("--fasta".into());
        args.push(corpus().join("reference.fa.gz").display().to_string());
    }
    args.extend(extra.iter().map(|s| s.to_string()));
    args
}

fn vep(dir: &Path, input: &Path, out: &Path, fasta: bool, extra: &[&str]) -> Output {
    common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .args(args(out, fasta, extra))
        .stdin(Stdio::null())
        .output()
        .unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The output's lines that are not `#` headers.
fn rows(out: &Path) -> Vec<String> {
    String::from_utf8_lossy(&std::fs::read(out).unwrap())
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// The rows of a successful run on `input` with the reference and `extra` flags.
fn annotate(dir: &Path, input: &Path, out: &Path, extra: &[&str]) -> Vec<String> {
    let output = vep(dir, input, out, true, extra);
    assert!(
        output.status.success(),
        "{} {extra:?}: {}",
        input.display(),
        stderr(&output)
    );
    rows(out)
}

/// A row's `Uploaded_variation` column.
fn name(row: &str) -> &str {
    row.split_once('\t').unwrap().0
}

/// A row past its `Uploaded_variation` column.
fn after_name(row: &str) -> &str {
    row.split_once('\t').unwrap().1
}

/// HGVS input needs the reference the notations are resolved against: without
/// `--fasta` the run stops before any output, under an explicit `--format` and
/// under the detection of the format from the first line, as Ensembl VEP's own
/// parser refuses the format offline (Parser/HGVS.pm:104).
#[test]
fn hgvs_input_needs_a_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(dir.path(), "records.hgvs", "21:g.30282169C>A\n");
    for extra in [&["--format", "hgvs"][..], &[][..]] {
        let out = dir.path().join("refused.txt");
        let output = vep(dir.path(), &input, &out, false, extra);
        assert_eq!(output.status.code(), Some(1), "{extra:?}");
        assert_eq!(
            stderr(&output),
            "ERROR: [input] Cannot use HGVS format in offline mode without a FASTA file\n",
            "{extra:?}"
        );
        assert!(!out.exists(), "{extra:?}: nothing is written");
    }
}

/// Every notation of the corpus annotates as its VCF twin does, row for row
/// past the name column, under `--format hgvs` and under format detection, and
/// `Uploaded_variation` is the notation itself, as Ensembl VEP names the
/// record (Parser/HGVS.pm:237): the accession `NC_000021.9` and the `chr21`
/// spelling print their rows at `21:`, the reference FASTA's name for the
/// chromosome, and the `(p.=)` suffix stays in the name.
#[test]
fn hgvs_records_annotate_as_their_vcf_twins_and_are_named_by_the_notation() {
    let dir = tempfile::tempdir().unwrap();
    let vcf_rows = annotate(
        dir.path(),
        &corpus().join("variants.vcf"),
        &dir.path().join("vcf.txt"),
        &[],
    );
    let notations: Vec<String> = std::fs::read_to_string(corpus().join("variants.hgvs.txt"))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    assert_eq!(notations.len(), 25);
    for (label, extra) in [("explicit", &["--format", "hgvs"][..]), ("guess", &[][..])] {
        let hgvs_rows = annotate(
            dir.path(),
            &corpus().join("variants.hgvs.txt"),
            &dir.path().join(format!("{label}.txt")),
            extra,
        );
        assert_eq!(hgvs_rows.len(), vcf_rows.len(), "{label}");
        for (hgvs_row, vcf_row) in hgvs_rows.iter().zip(&vcf_rows) {
            assert_eq!(after_name(hgvs_row), after_name(vcf_row), "{label}");
        }
        let names: Vec<&str> = hgvs_rows.iter().map(|row| name(row)).collect();
        for notation in &notations {
            assert!(names.contains(&notation.as_str()), "{label}: {notation}");
        }
        for notation in ["NC_000021.9:g.30282169C>A", "chr21:g.30282169C>A"] {
            let row = hgvs_rows.iter().find(|row| name(row) == notation).unwrap();
            assert_eq!(
                row.split('\t').nth(1),
                Some("21:30282169"),
                "{label}: {row}"
            );
        }
    }
}

/// A notation Ensembl VEP cannot resolve is skipped and warned about with the
/// resolver's message (Parser/HGVS.pm:204-216), the line quoted as every skip
/// of a parser-fed record quotes it: a reference allele the genome contradicts,
/// a repeat, a notation with no type letter. A variant the resolver builds but
/// `validate_vf` refuses (Parser.pm:564-592) is skipped likewise: an insertion
/// at one position, and an upper-case `DEL`, whose listed bases (none) stay its
/// reference allele. The notations around them annotate.
#[test]
fn unresolvable_notations_are_skipped_with_the_warning() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(
        dir.path(),
        "mixed.hgvs",
        "21:g.30282169G>A\n21:g.30282169C>A\n21:g.30282169C[3]\n21:30282169C>A\n\
         21:g.30282172insG\n21:g.30282172DEL\n21:g.30282172del\n",
    );
    let out = dir.path().join("mixed.txt");
    let rows = annotate(dir.path(), &input, &out, &["--format", "hgvs"]);
    let names: Vec<&str> = rows.iter().map(|row| name(row)).collect();
    assert_eq!(names, ["21:g.30282169C>A", "21:g.30282172del"]);
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(
        warnings,
        "WARNING: line 1 skipped (21:g.30282169G>A): Unable to parse HGVS notation '21:g.30282169G>A': \
         Reference allele extracted from 21:30282169-30282169 (C) does not match reference \
         allele given by HGVS notation 21:g.30282169G>A (G)\n\
         WARNING: line 3 skipped (21:g.30282169C[3]): Unable to parse HGVS notation '21:g.30282169C[3]'\n\
         WARNING: line 4 skipped (21:30282169C>A): Unable to parse HGVS notation '21:30282169C>A': \
         HGVS type  not recognised for 21:30282169C>A\n\
         WARNING: line 5 skipped (21:g.30282172insG): Alleles look like an insertion (-/G) \
         but coordinates are not start = end + 1 (START=30282172, END=30282172) \
         for variant (21:g.30282172insG) \n\
         WARNING: line 6 skipped (21:g.30282172DEL): Length of reference allele ( length 0) \
         does not match coordinates 30282172-30282172\n"
    );
}

/// The span a genomic notation describes must be shorter than `--max_sv_size`
/// (VariationFeatureAdaptor.pm:1886): a twenty-base deletion annotates under
/// the default limit and is skipped under `--max_sv_size 10` with Ensembl
/// VEP's message and its hint.
#[test]
fn max_sv_size_bounds_the_span_of_a_notation() {
    let dir = tempfile::tempdir().unwrap();
    let input = write(dir.path(), "span.hgvs", "21:g.30282180_30282199del\n");
    let rows = annotate(
        dir.path(),
        &input,
        &dir.path().join("default.txt"),
        &["--format", "hgvs"],
    );
    assert!(!rows.is_empty());
    assert!(rows[0].starts_with("21:g.30282180_30282199del\t21:30282180-30282199\t-\t"));
    let out = dir.path().join("limited.txt");
    let rows = annotate(
        dir.path(),
        &input,
        &out,
        &["--format", "hgvs", "--max_sv_size", "10"],
    );
    assert!(rows.is_empty(), "{rows:?}");
    let warnings = std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
    assert_eq!(
        warnings,
        "WARNING: line 1 skipped (21:g.30282180_30282199del): Unable to parse HGVS notation '21:g.30282180_30282199del': \
         Region requested must be smaller than 10. This value can be adjusted with --max_sv_size.\n"
    );
}
