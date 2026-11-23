// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Two transcripts of a merged cache that share a `dbID`, through the `vep`
//! binary.
//!
//! A merged cache joins the region files of an Ensembl-set cache and a
//! RefSeq-set cache (`Pipeline/DumpVEP/MergeVEP.pm:94-120`), and each transcript
//! keeps the `dbID` of its own set, so one number can name a transcript of each
//! set. Ensembl VEP keeps one transcript per `dbID` among those a buffer loads
//! (`AnnotationType/Transcript.pm:283-287`) and writes no row for the other;
//! vep-rs keys a transcript on its set and its `dbID` and annotates both.
//!
//! The cache is the committed conversion of the synthetic merged cache that
//! `scripts/divergences/reproductions/inputs/merged_same_dbid_collapse/make_cache.pl`
//! builds: `XR_000001.1` (RefSeq set) and `ENSEXT00000000001` (Ensembl set)
//! overlap in one region file under `dbID` 7; `ENSEXT00000000002` (Ensembl set)
//! and `XR_000002.1` (RefSeq set) lie in two region files under `dbID` 9. Each
//! expected row is the row Ensembl VEP 116.3 writes for the record on the same
//! cache built with four distinct `dbID`s.

use std::fs;
use std::path::PathBuf;
use std::process::Stdio;

mod common;

const SHARED_LOCUS: &str = "contig\t101500\tshared_locus\tA\tG\t.\t.\t.";
const FIRST_OF_FAR_PAIR: &str = "contig\t501000\tfirst_of_far_pair\tA\tG\t.\t.\t.";
const SECOND_OF_FAR_PAIR: &str = "contig\t3201000\tsecond_of_far_pair\tA\tG\t.\t.\t.";

const ENSEMBL_ROW_OF_THE_SHARED_LOCUS: &str =
    "shared_locus\tcontig:101500\tG\tENSEXG00000000001\tENSEXT00000000001\tTranscript\t\
     non_coding_transcript_exon_variant\t500\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl";
const REFSEQ_ROW_OF_THE_SHARED_LOCUS: &str =
    "shared_locus\tcontig:101500\tG\t900001\tXR_000001.1\tTranscript\t\
     non_coding_transcript_exon_variant\t1500\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq";
const ROW_OF_THE_FIRST_OF_THE_FAR_PAIR: &str =
    "first_of_far_pair\tcontig:501000\tG\tENSEXG00000000002\tENSEXT00000000002\tTranscript\t\
     non_coding_transcript_exon_variant\t1000\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER;STRAND=1;SOURCE=Ensembl";
const ROW_OF_THE_SECOND_OF_THE_FAR_PAIR: &str =
    "second_of_far_pair\tcontig:3201000\tG\t900002\tXR_000002.1\tTranscript\t\
     non_coding_transcript_exon_variant\t1000\t-\t-\t-\t-\t-\t\
     IMPACT=MODIFIER;STRAND=1;SOURCE=RefSeq";

/// The committed conversion of the synthetic merged cache.
fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../scripts/divergences/reproductions/inputs/merged_same_dbid_collapse/json_cache")
}

/// Runs `vep --merged` on `records` with `flags` and returns the output's
/// data lines (those not opening with `#`).
fn rows(records: &[&str], flags: &[&str]) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("records.vcf");
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for record in records {
        text.push_str(record);
        text.push('\n');
    }
    fs::write(&input, text).unwrap();
    let out = dir.path().join("out.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["--offline", "--json_cache"])
        .arg(cache())
        .args([
            "--species",
            "two_set_example",
            "--assembly",
            "EXAMPLE1",
            "--merged",
            "--no_stats",
            "--force_overwrite",
        ])
        .args(flags)
        .stdin(Stdio::null())
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
        .map(str::to_string)
        .collect()
}

/// Both transcripts of one region file that share a `dbID` get their row on
/// the record they overlap.
#[test]
fn two_sets_sharing_a_dbid_in_one_region_file_are_both_annotated() {
    assert_eq!(
        rows(&[SHARED_LOCUS], &[]),
        [
            ENSEMBL_ROW_OF_THE_SHARED_LOCUS,
            REFSEQ_ROW_OF_THE_SHARED_LOCUS
        ]
    );
}

/// A transcript keeps its row when an earlier record of the input lies on the
/// transcript of the other set that shares its `dbID`.
#[test]
fn two_sets_sharing_a_dbid_in_two_region_files_are_both_annotated() {
    assert_eq!(
        rows(&[FIRST_OF_FAR_PAIR, SECOND_OF_FAR_PAIR], &[]),
        [
            ROW_OF_THE_FIRST_OF_THE_FAR_PAIR,
            ROW_OF_THE_SECOND_OF_THE_FAR_PAIR
        ]
    );
}

/// The rows of the three records are the same whatever the buffer holds: all
/// records in one buffer, and one record a buffer.
#[test]
fn the_rows_do_not_change_with_the_buffer_size() {
    let records = [SHARED_LOCUS, FIRST_OF_FAR_PAIR, SECOND_OF_FAR_PAIR];
    let expected = [
        ENSEMBL_ROW_OF_THE_SHARED_LOCUS,
        REFSEQ_ROW_OF_THE_SHARED_LOCUS,
        ROW_OF_THE_FIRST_OF_THE_FAR_PAIR,
        ROW_OF_THE_SECOND_OF_THE_FAR_PAIR,
    ];
    assert_eq!(rows(&records, &[]), expected);
    assert_eq!(rows(&records, &["--buffer_size", "1"]), expected);
}
