// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field HGVSp` on the four made-up runs of `tests/fixtures/stop_codon/`
//! and on the one of `tests/fixtures/frameshift/`, one row a pair of the pair table of the
//! same name: an `HGVSp` is an entry of `field_divergences` under
//! `hgvsp_stop_codon_protein_comparison`, `hgvsp_stop_named_on_incomplete_coding_sequence`
//! or `hgvsp_frameshift_with_no_changed_residue` where the table says the class takes the
//! pair and an `unexplained_residual` where it does not, and the manifest is written to the
//! bytes of the fixture's `expected.1.json`.

mod classify_support;
mod common;
mod recorded_runs;

use common::fixture;
use recorded_runs::{call, classes, Recorded, UNEXPLAINED};
use std::collections::BTreeMap;

/// Runs the fixture `name` of the directory `family`, whose record of ordinal `i` stands at
/// 21:100(i+1) and whose row is written one position on, and holds the class of each row to
/// the last column of the pair table: `class` for 1, `unexplained_residual` for 0.
fn pairs_run_of(family: &str, name: &str, class: &str, summary: (usize, &str)) {
    let table = std::fs::read_to_string(fixture(family).join(format!("{name}.tsv"))).unwrap();
    let taken: Vec<bool> = table
        .lines()
        .map(|line| line.rsplit('\t').next() == Some("1"))
        .collect();
    assert_eq!(taken.len(), summary.0);
    let run = Recorded::new(&format!("{family}/{name}"));
    let entries = run.holds(&call(&["HGVSp"], summary, "expected.1.json"));
    let locations: Vec<String> = (0..taken.len())
        .map(|ordinal| format!("21:{}", 100 * (ordinal + 1) + 1))
        .collect();
    let expected: BTreeMap<(&str, &str), &str> = locations
        .iter()
        .zip(&taken)
        .map(|(location, taken)| {
            let class = if *taken { class } else { UNEXPLAINED };
            ((location.as_str(), "HGVSp"), class)
        })
        .collect();
    assert_eq!(
        classes(&entries, |entry| entry.location.as_str()),
        expected,
        "{name}"
    );
}

/// `pairs_run_of` on a fixture of `tests/fixtures/stop_codon/`.
fn pairs_run(name: &str, class: &str, summary: (usize, &str)) {
    pairs_run_of("stop_codon", name, class, summary);
}

#[test]
fn classify_names_the_stop_codon_protein_comparison_class() {
    // An `HGVSp` at the stop codon of one protein, Ensembl VEP's description read from the
    // codon window and vep-rs's the one the two proteins give, is
    // `hgvsp_stop_codon_protein_comparison` in each of its first nine shapes; the nearest
    // shapes stay an `unexplained_residual`.
    pairs_run(
        "protein_comparison",
        "hgvsp_stop_codon_protein_comparison",
        (
            67,
            "{'hgvsp_stop_codon_protein_comparison': 23, 'unexplained_residual': 44}",
        ),
    );
}

#[test]
fn classify_names_a_residue_and_a_stop_inserted_before_the_last_residue() {
    // An `HGVSp` in which Ensembl VEP writes one residue and a stop inserted between two
    // neighbouring residues, the residue the second neighbour's, on a row whose terms hold
    // `stop_retained_variant`, and vep-rs the stop unchanged at the residue after that
    // neighbour, is `hgvsp_stop_codon_protein_comparison`; the nearest shapes stay an
    // `unexplained_residual`.
    pairs_run(
        "residue_and_stop_inserted",
        "hgvsp_stop_codon_protein_comparison",
        (
            19,
            "{'unexplained_residual': 16, 'hgvsp_stop_codon_protein_comparison': 3}",
        ),
    );
}

#[test]
fn classify_names_the_stop_named_on_incomplete_coding_sequence_class() {
    // An `HGVSp` in which Ensembl VEP names the stop alone on a row that shows an incomplete
    // coding sequence, where vep-rs writes none, is
    // `hgvsp_stop_named_on_incomplete_coding_sequence`; the nearest shapes stay an
    // `unexplained_residual`.
    pairs_run(
        "named_on_incomplete_coding_sequence",
        "hgvsp_stop_named_on_incomplete_coding_sequence",
        (
            19,
            "{'unexplained_residual': 10, 'hgvsp_stop_named_on_incomplete_coding_sequence': 9}",
        ),
    );
}

#[test]
fn classify_reads_the_stop_named_alone_inside_prediction_parentheses() {
    // The stop named alone inside one pair of prediction parentheses, on a row that shows an
    // incomplete coding sequence where vep-rs writes none, is
    // `hgvsp_stop_named_on_incomplete_coding_sequence` as it is without them; a parenthesis
    // without its pair and a second pair around the first stay an `unexplained_residual`.
    pairs_run(
        "named_inside_prediction_parentheses",
        "hgvsp_stop_named_on_incomplete_coding_sequence",
        (
            13,
            "{'unexplained_residual': 8, 'hgvsp_stop_named_on_incomplete_coding_sequence': 5}",
        ),
    );
}

#[test]
fn classify_names_the_frameshift_with_no_changed_residue_class() {
    // An `HGVSp` in which Ensembl VEP writes a frameshift with the unknown count that names
    // one residue, no stop, on both sides of its number, alone or inside one pair of
    // prediction parentheses, on a row whose terms hold `frameshift_variant`, where vep-rs
    // writes none, is `hgvsp_frameshift_with_no_changed_residue`; the nearest shapes stay an
    // `unexplained_residual`.
    pairs_run_of(
        "frameshift",
        "no_changed_residue",
        "hgvsp_frameshift_with_no_changed_residue",
        (
            19,
            "{'unexplained_residual': 15, 'hgvsp_frameshift_with_no_changed_residue': 4}",
        ),
    );
}
