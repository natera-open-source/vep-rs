// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The direction term of a `<CNV:TR>` allele that is read symbolically, on the
//! rows the binary writes for the tandem-repeat corpus: its own run, which
//! reads no reference sequence, and its run under `--max_sv_size 10`, which
//! keeps a record symbolic when one of its alleles is over ten bases.
//!
//! A record states the length of an alternate allele in `RUS` with `RUC` or
//! `RB`, and the length of the run it replaces in `SVLEN` or `END`. Where the
//! allele is shorter than the run the row names `feature_truncation`. Ensembl
//! VEP 116.3 names `feature_elongation` on every such row, since
//! `copy_number_gain` holds for every `tandem_repeat` whatever its counts
//! (Utils/VariationEffect.pm:246-256) and `feature_truncation` asks for a
//! loss or a deletion (:350-373). Where the record states the allele to be
//! exactly as long as a run of more than one base, and no other length for it,
//! the row names neither term and carries the IMPACT of the first term left,
//! where Ensembl VEP 116.3 names `feature_elongation` as well. Every other cell
//! is the one Ensembl VEP writes, and so is every row of an allele that is
//! longer than the run or whose length the record does not state.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

mod common;

/// The alleles of the corpus that state fewer bases than their run, each as
/// its record's ID and the `ALLELE_NUM` of its rows: one unit fewer on line 4,
/// the second allele of line 8 and the second allele of line 9.
const STATED_LOSSES: [(&str, &str); 3] = [
    ("tr_loss", "ALLELE_NUM=1"),
    ("tr_multi", "ALLELE_NUM=2"),
    ("tr_two_svlen", "ALLELE_NUM=2"),
];

/// The allele of the corpus that states as many bases as its run, as its
/// record's ID and the `ALLELE_NUM` of its rows: three units on line 5.
const AS_LONG_AS_ITS_RUN: (&str, &str) = ("tr_same", "ALLELE_NUM=1");

/// The row of line 4 on the first coding transcript, whole.
const LOSS_ROW: &str = "tr_loss\t21:13620284-13620292\ttandem_repeat\tENSG00000166351\tENST00000299443\tTranscript\tfeature_truncation,coding_sequence_variant\t1496-1504\t1044-1052\t348-351\t-\t-\t-\tALLELE_NUM=1;IMPACT=HIGH;STRAND=1;VARIANT_CLASS=tandem_repeat;OverlapBP=9;OverlapPC=0.02";

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The data rows of a default-format output.
fn data_rows(text: &str) -> Vec<String> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// The rows Ensembl VEP 116.3 wrote for the corpus under `expected_dir`.
fn stored_rows(expected_dir: &str) -> Vec<String> {
    let path = corpus("GRCh38-vcf-tandem-repeats")
        .join(expected_dir)
        .join("default.txt.gz");
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    data_rows(&text)
}

/// The rows the binary writes for the corpus's input under the corpus's own
/// flags and `extra`.
fn written_rows(extra: &[&str]) -> Vec<String> {
    let dir = tempfile::tempdir().unwrap();
    let tandem = corpus("GRCh38-vcf-tandem-repeats");
    let out = dir.path().join("out.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(tandem.join("variants.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(tandem.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--buffer_size", "5000", "--force_overwrite", "--no_stats"])
        .args(["--quiet", "--uploaded_allele", "--allele_number"])
        .arg("--variant_class")
        .args(extra)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "vep ends with exit 0: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    data_rows(&fs::read_to_string(out).unwrap())
}

/// Whether `row` is a row of an allele that states a loss.
fn states_a_loss(row: &str) -> bool {
    let id = row.split('\t').next().unwrap_or("");
    let extra = row.rsplit('\t').next().unwrap_or("");
    STATED_LOSSES
        .iter()
        .any(|(record, allele)| id == *record && extra.split(';').any(|item| item == *allele))
}

/// Ensembl VEP's rows with the direction term of each symbolic allele that
/// states a loss read as a truncation, and the number of rows so read.
fn with_stated_losses_named(stored: Vec<String>) -> (Vec<String>, usize) {
    let mut named = 0;
    let rows = stored
        .into_iter()
        .map(|row| {
            let symbolic = row.split('\t').nth(2) == Some("tandem_repeat");
            if symbolic && states_a_loss(&row) && row.contains("\tfeature_elongation,") {
                named += 1;
                row.replacen("\tfeature_elongation,", "\tfeature_truncation,", 1)
            } else {
                row
            }
        })
        .collect();
    (rows, named)
}

/// `rows` with each row of the symbolic allele as long as its run read without
/// `feature_elongation` and with the IMPACT of `coding_sequence_variant`, the
/// one term left on it (Utils/Constants.pm:1061-1078), and the number of rows
/// so read.
fn with_no_length_term_on_the_allele_as_long_as_its_run(rows: Vec<String>) -> (Vec<String>, usize) {
    let (record, allele) = AS_LONG_AS_ITS_RUN;
    let mut read = 0;
    let rows = rows
        .into_iter()
        .map(|row| {
            let mut cells = row.split('\t');
            let id = cells.next() == Some(record);
            let symbolic = cells.nth(1) == Some("tandem_repeat");
            let numbered = row
                .rsplit('\t')
                .next()
                .is_some_and(|extra| extra.split(';').any(|item| item == allele));
            let gain = "\tfeature_elongation,coding_sequence_variant\t";
            if id && symbolic && numbered && row.contains(gain) {
                read += 1;
                row.replacen(gain, "\tcoding_sequence_variant\t", 1)
                    .replacen(";IMPACT=HIGH;", ";IMPACT=MODIFIER;", 1)
            } else {
                row
            }
        })
        .collect();
    (rows, read)
}

#[test]
fn an_allele_that_states_fewer_bases_than_its_run_is_named_a_truncation_without_a_reference() {
    let (expected, named) = with_stated_losses_named(stored_rows("expected"));
    assert_eq!(
        named, 6,
        "three alleles state a loss, each on two coding transcripts"
    );
    let (expected, read) = with_no_length_term_on_the_allele_as_long_as_its_run(expected);
    assert_eq!(
        read, 2,
        "one allele is as long as its run, on two coding transcripts"
    );
    let written = written_rows(&[]);
    assert!(
        written.iter().any(|row| row == LOSS_ROW),
        "the row of line 4 on ENST00000299443 is written whole"
    );
    assert_eq!(written, expected);
}

#[test]
fn an_allele_kept_symbolic_for_its_records_size_is_named_by_its_own_stated_length() {
    let reference = corpus("GRCh38-identifiers").join("reference.fa.gz");
    let (expected, named) = with_stated_losses_named(stored_rows("expected_fasta-max10"));
    assert_eq!(
        named, 4,
        "the second alleles of lines 8 and 9 stay symbolic, each on two coding transcripts"
    );
    let (expected, read) = with_no_length_term_on_the_allele_as_long_as_its_run(expected);
    assert_eq!(
        read, 0,
        "line 5 is expanded to its own reference run and writes no row"
    );
    let written = written_rows(&[
        "--fasta",
        reference.to_str().unwrap(),
        "--max_sv_size",
        "10",
    ]);
    assert_eq!(written, expected);
}
