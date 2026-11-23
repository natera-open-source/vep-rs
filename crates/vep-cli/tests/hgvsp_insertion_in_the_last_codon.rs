// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary writes HGVSp for an insertion between two bases of
//! the last residue's codon that leaves the protein as it was.
//! The first inserted base completes that codon to the same residue and the
//! next three are a stop codon, which stands at the stop's own residue, so
//! the alternate protein is the reference protein and the notation is the
//! stop at its residue number written as unchanged, at five, six, eight and
//! nine inserted bases on either strand.
//! Ensembl VEP types the clipped codon window an insertion
//! (TranscriptVariationAllele.pm:2063-2064), finds no residue after the last
//! one to flank it (TranscriptVariationAllele.pm:2368-2380) and writes no
//! HGVSp (TranscriptVariationAllele.pm:2146), except where its terms hold
//! `frameshift_variant` and it walks the two translations to the stop
//! (TranscriptVariationAllele.pm:2348-2357).
//! A record of one alternate allele that holds a base its two alleles share
//! is read as the insertion alone; in a record of two alternate alleles the
//! allele keeps the shared base, and the notation is the same.
//! The records are made up on the chromosome 21 reference the corpus carries
//! and read on its release-116 GRCh38 cache, under `--shift_hgvs 0`.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Stdio;

mod common;

/// The identifiers corpus: a chromosome 21 cache slice and the bgzipped
/// reference that holds the bases of both transcripts' last codons.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-identifiers")
}

/// The forward-strand transcript whose protein has 365 residues: the codon of
/// residue 365 is coding bases 1093 to 1095, at 21:17565687-17565689, and the
/// stop codon follows it.
const FORWARD: &str = "ENST00000284878";
const FORWARD_STOP_UNCHANGED: &str = "ENSP00000284878.7:p.Ter366%3D";

/// The reverse-strand transcript whose protein has 905 residues: the codon of
/// residue 905 is coding bases 2713 to 2715, at 21:29537235-29537233, and the
/// stop codon follows it.
const REVERSE: &str = "ENST00000389125";
const REVERSE_STOP_UNCHANGED: &str = "ENSP00000373777.3:p.Ter906%3D";

/// The terms of an insertion that keeps the stop codon.
const KEPT_STOP: &str = "inframe_insertion,stop_retained_variant";

/// The terms both engines write for five or eight bases held with a base the
/// two alleles share, where the allele's span is that base alone.
const SHIFTED_FRAME: &str = "stop_gained,frameshift_variant";

/// Each record inserts five, six, eight or nine bases between the second and
/// third base of the last residue's codon. A `pure` record writes them after
/// the base before the gap; a `held` record replaces the base after the gap
/// by the inserted bases and that base, beside a second alternate allele that
/// keeps the parser from reading the record as the insertion alone.
const RECORDS: &[(&str, &str, &str)] = &[
    ("21\t17565688\tfwd_5_pure\tT\tTCTGAG", FORWARD, KEPT_STOP),
    ("21\t17565688\tfwd_6_pure\tT\tTCTGAGG", FORWARD, KEPT_STOP),
    ("21\t17565688\tfwd_8_pure\tT\tTCTGAGGAG", FORWARD, KEPT_STOP),
    (
        "21\t17565688\tfwd_9_pure\tT\tTCTGAGGAGT",
        FORWARD,
        KEPT_STOP,
    ),
    (
        "21\t17565689\tfwd_5_held\tA\tCTGAGA,G",
        FORWARD,
        SHIFTED_FRAME,
    ),
    ("21\t17565689\tfwd_6_held\tA\tCTGAGGA,G", FORWARD, KEPT_STOP),
    (
        "21\t17565689\tfwd_8_held\tA\tCTGAGGAGA,G",
        FORWARD,
        SHIFTED_FRAME,
    ),
    (
        "21\t17565689\tfwd_9_held\tA\tCTGAGGAGTA,G",
        FORWARD,
        KEPT_STOP,
    ),
    ("21\t29537233\trev_5_pure\tC\tCCTCAT", REVERSE, KEPT_STOP),
    ("21\t29537233\trev_6_pure\tC\tCCCTCAT", REVERSE, KEPT_STOP),
    ("21\t29537233\trev_8_pure\tC\tCCTCCTCAT", REVERSE, KEPT_STOP),
    (
        "21\t29537233\trev_9_pure\tC\tCACTCCTCAT",
        REVERSE,
        KEPT_STOP,
    ),
    (
        "21\t29537234\trev_5_held\tG\tCTCATG,A",
        REVERSE,
        SHIFTED_FRAME,
    ),
    ("21\t29537234\trev_6_held\tG\tCCTCATG,A", REVERSE, KEPT_STOP),
    (
        "21\t29537234\trev_8_held\tG\tCTCCTCATG,A",
        REVERSE,
        SHIFTED_FRAME,
    ),
    (
        "21\t29537234\trev_9_held\tG\tACTCCTCATG,C",
        REVERSE,
        KEPT_STOP,
    ),
];

/// The record's name, its third column.
fn name(record: &str) -> &str {
    record.split('\t').nth(2).expect("a record of five columns")
}

/// The Consequence and the HGVSp of the row each record's inserting allele
/// has on its transcript, by the record's name.
fn annotated(extra: &[&str]) -> BTreeMap<String, (String, String)> {
    let dir = tempfile::tempdir().expect("a scratch directory");
    let input = dir.path().join("records.vcf");
    let output = dir.path().join("out.txt");
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (record, _, _) in RECORDS {
        text += &format!("{record}\t.\t.\t.\n");
    }
    std::fs::write(&input, text).expect("the input is written");
    let corpus = corpus();
    let status = common::isolated_vep(dir.path())
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38"])
        .args(["--json_cache", corpus.join("json_cache").to_str().unwrap()])
        .args(["--fasta", corpus.join("reference.fa.gz").to_str().unwrap()])
        .args(["--hgvs", "--shift_hgvs", "0", "--tab"])
        .args([
            "--fields",
            "Uploaded_variation,Allele,Feature,Consequence,HGVSp",
        ])
        .args(["--no_stats", "--quiet", "--force_overwrite"])
        .args(extra)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("vep runs");
    assert!(status.success(), "vep exits 0");
    let transcript_of: BTreeMap<&str, &str> = RECORDS
        .iter()
        .map(|(record, transcript, _)| (name(record), *transcript))
        .collect();
    let mut rows = BTreeMap::new();
    for line in std::fs::read_to_string(&output)
        .expect("the output is read")
        .lines()
        .filter(|line| !line.starts_with('#'))
    {
        let cells: Vec<&str> = line.split('\t').collect();
        // The second allele of a `held` record is one base.
        if transcript_of.get(cells[0]) == Some(&cells[2]) && cells[1].len() > 1 {
            let row = (cells[3].to_string(), cells[4].to_string());
            assert!(
                rows.insert(cells[0].to_string(), row).is_none(),
                "one row a record: {line}"
            );
        }
    }
    rows
}

/// The stop written as unchanged for `transcript`.
fn stop_unchanged(transcript: &str) -> &'static str {
    if transcript == FORWARD {
        FORWARD_STOP_UNCHANGED
    } else {
        REVERSE_STOP_UNCHANGED
    }
}

/// Every record carries the stop at its residue number as unchanged: residue
/// 366 on the forward-strand transcript and residue 906 on the reverse-strand
/// one, at every length, for the insertion alone and for the allele that
/// holds a shared base. The terms are those of the allele's own span.
#[test]
fn insertion_in_the_last_codon_that_keeps_the_protein_is_the_stop_unchanged_at_every_length() {
    let rows = annotated(&[]);
    assert_eq!(rows.len(), RECORDS.len());
    for (record, transcript, terms) in RECORDS {
        let name = name(record);
        let (consequence, hgvsp) = &rows[name];
        assert_eq!(consequence, terms, "{name}");
        assert_eq!(hgvsp, stop_unchanged(transcript), "{name}");
    }
}

/// Under `--hgvsp_use_prediction` the same value stands inside parentheses.
#[test]
fn insertion_in_the_last_codon_that_keeps_the_protein_is_the_stop_unchanged_in_the_prediction_form()
{
    let rows = annotated(&["--hgvsp_use_prediction"]);
    assert_eq!(rows.len(), RECORDS.len());
    for (record, transcript, _) in RECORDS {
        let name = name(record);
        let predicted = stop_unchanged(transcript).replacen(":p.", ":p.(", 1) + ")";
        assert_eq!(rows[name].1, predicted, "{name}");
    }
}
