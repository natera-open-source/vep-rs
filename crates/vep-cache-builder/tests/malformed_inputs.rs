// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Malformed inputs the builder refuses or reads whole rather than failing on:
//! a GFF3 whose exons mix ranked and unranked lines, a seqid that would place a
//! chromosome's directory outside the cache, and a FASTA index naming more bases
//! than the file holds. Each fixture under `tests/fixtures/malformed/` is the
//! shape as a file; a failed build prints one `ERROR: [cache-builder]` line and
//! exits 1, and a refused build writes nothing outside the output directory.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use serde_json::Value;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/malformed")
        .join(name)
}

/// Runs the builder on `gff3` and `peptides` into `out` with `extra` flags.
fn build(gff3: &Path, peptides: &Path, out: &Path, extra: &[(&str, &Path)]) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"));
    cmd.args(["--release", "116", "--assembly", "GRCh37"])
        .arg("--gff3")
        .arg(gff3)
        .arg("--fasta")
        .arg(peptides)
        .arg("--output-dir")
        .arg(out);
    for (flag, path) in extra {
        cmd.arg(flag).arg(path);
    }
    cmd.output().expect("vep-cache-builder runs")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
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

/// Exit 1 with `line` as the one `ERROR: [` line on standard error.
fn assert_refused(output: &Output, line: &str) {
    let err = stderr(output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let errors: Vec<&str> = err.lines().filter(|l| l.starts_with("ERROR: [")).collect();
    assert_eq!(errors, [line], "{err}");
}

/// Sixty-four exons of one forward-strand transcript, 33 with a `rank`
/// attribute and 31 without, listed in no particular order: the build
/// completes, and the cache carries the exons in genomic order, numbered 1 to
/// 64 in that order, as the transcript order falls back to the start coordinate
/// when any exon lacks a rank and the file's ranks go with it.
#[test]
fn exons_mixing_ranked_and_unranked_lines_build_in_start_order() {
    let out = tempfile::tempdir().unwrap();
    let output = build(
        &fixture("mixed_exon_ranks.gff3"),
        &fixture("mixed_exon_ranks_peptides.fa"),
        out.path(),
        &[],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let shards = entries(&out.path().join("transcripts/21"));
    assert_eq!(shards, ["1-1000000.json"]);
    let text = std::fs::read_to_string(out.path().join("transcripts/21/1-1000000.json")).unwrap();
    let records: Vec<Value> = serde_json::from_str(&text).unwrap();
    assert_eq!(records.len(), 1);
    let exons = records[0]["exons"].as_array().unwrap();
    assert_eq!(exons.len(), 64);
    let starts: Vec<u64> = exons.iter().map(|e| e["start"].as_u64().unwrap()).collect();
    let mut sorted = starts.clone();
    sorted.sort_unstable();
    assert_eq!(starts, sorted);
    let ranks: Vec<u64> = exons.iter().map(|e| e["rank"].as_u64().unwrap()).collect();
    assert_eq!(ranks, (1..=64).collect::<Vec<u64>>());
}

/// A seqid of `../../escaped`, and one that is an absolute path, name a
/// directory outside the output: the build is refused at the GFF3, naming the
/// seqid and its line, and the output directory's parent holds nothing but the
/// output directory, which is empty.
#[test]
fn a_seqid_escaping_the_output_directory_is_refused_and_writes_nothing() {
    let parent = tempfile::tempdir().unwrap();
    let out = parent.path().join("cache").join("deep");
    std::fs::create_dir_all(&out).unwrap();
    let output = build(
        &fixture("seqid_parent_escape.gff3"),
        &fixture("seqid_escape_peptides.fa"),
        &out,
        &[],
    );
    assert_refused(
        &output,
        "ERROR: [cache-builder] Failed to parse GFF3: invalid chromosome name '../../escaped' at line 2",
    );
    assert_eq!(entries(parent.path()), ["cache"]);
    assert_eq!(entries(&parent.path().join("cache")), ["deep"]);
    assert!(entries(&out).is_empty(), "{:?}", entries(&out));

    let target = parent.path().join("elsewhere");
    let absolute = std::fs::read_to_string(fixture("seqid_parent_escape.gff3"))
        .unwrap()
        .replace("../../escaped", target.to_str().unwrap());
    let gff3 = parent.path().join("absolute.gff3");
    std::fs::write(&gff3, absolute).unwrap();
    let output = build(&gff3, &fixture("seqid_escape_peptides.fa"), &out, &[]);
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to parse GFF3: invalid chromosome name '{}' at line 2",
            target.display()
        ),
    );
    assert_eq!(entries(parent.path()), ["absolute.gff3", "cache"]);
    assert!(entries(&out).is_empty(), "{:?}", entries(&out));
}

/// A `.fai` naming a 4e18-base chromosome 21 over a 126-byte FASTA of 120
/// bases, with a GFF3 whose transcript spans 3e18 bases: the genome FASTA is
/// refused as it opens, before any sequence is read, naming the bases the file
/// lacks (the 120 bases after the 4-byte header are held, the rest are not).
#[test]
fn a_fasta_index_naming_bases_beyond_the_file_is_refused() {
    let out = tempfile::tempdir().unwrap();
    let output = build(
        &fixture("index_overstates_contig.gff3"),
        &fixture("index_overstates_contig_peptides.fa"),
        out.path(),
        &[("--genome-fasta", &fixture("index_overstates_contig.fa"))],
    );
    assert_refused(
        &output,
        &format!(
            "ERROR: [cache-builder] Failed to load genome FASTA (needs .fa + .fai index): FASTA index names 3999999999999999880 bases beyond the end of {}",
            fixture("index_overstates_contig.fa").display()
        ),
    );
    assert!(entries(out.path()).is_empty(), "{:?}", entries(out.path()));
}
