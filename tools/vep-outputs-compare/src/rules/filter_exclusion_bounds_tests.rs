// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pass over `discordant.tsv` under the 115.2 reference, on reports written for each test:
//! a side's exclusions never exceed its rows, a bucket counts a pair whether or not the
//! adjusted F1 removes it, and the open file holds the rows no rule removed, as they stand.

use super::*;
use crate::input::KeyPrep;
use crate::TempDir;

const HEADER: &str = "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set";

/// The consequence set as the comparison keys it: its terms unique, sorted and comma-joined.
fn norm(csq: &str) -> String {
    let mut set = Vec::new();
    KeyPrep::default().normalized_set(csq, &mut set);
    String::from_utf8(set).unwrap()
}

/// One row of a discordance report under the file name `f`.
fn row(source: &str, key: [&str; 4], csq: &str) -> String {
    let [location, allele, feature, feature_type] = key;
    [
        "f",
        source,
        location,
        allele,
        feature,
        feature_type,
        &norm(csq),
    ]
    .join("\t")
}

/// One row on a transcript.
fn transcript_row(source: &str, location: &str, allele: &str, feature: &str, csq: &str) -> String {
    row(source, [location, allele, feature, "Transcript"], csq)
}

struct Pass {
    exclusions: Exclusions,
    /// The open file, whole.
    open: String,
    /// The names of the files in the directory of the report after the pass, sorted.
    files: Vec<String>,
}

impl Pass {
    fn excluded(&self) -> (u64, u64) {
        (self.exclusions.excluded_perl, self.exclusions.excluded_rust)
    }

    fn bucket(&self, name: &str) -> BucketStats {
        let (_, stats) = self
            .exclusions
            .buckets
            .iter()
            .find(|(bucket, _)| bucket == name)
            .expect("bucket");
        stats.clone()
    }

    /// The pairs the buckets report as excluded, summed.
    fn attributed(&self) -> u64 {
        let stats = self.exclusions.buckets.iter().flat_map(|(_, stats)| stats);
        stats
            .filter(|(field, _)| *field == "excluded")
            .map(|(_, n)| n)
            .sum()
    }

    fn open_lines(&self) -> Vec<&str> {
        self.open.lines().collect()
    }
}

/// Runs the pass under the 115.2 reference over a report of `rows`: no record in scope, no
/// row-aware rule and no output directory.
fn pass(rows: &[String]) -> Pass {
    let scratch = TempDir::create(&std::env::temp_dir(), "filter_bounds_test_").unwrap();
    let discordant_path = scratch.path.join("discordant.tsv");
    let mut report = format!("{HEADER}\n");
    for row in rows {
        report.push_str(row);
        report.push('\n');
    }
    std::fs::write(&discordant_path, report).unwrap();
    let open_path = scratch.path.join("discordant_open.tsv");
    let scope = ScopeFeatures::default();
    let exclusions = filter_intended_divergences(
        &discordant_path,
        &open_path,
        Collector::Release115(&scope),
        &[],
        None,
        None,
    )
    .unwrap();
    let mut files: Vec<String> = std::fs::read_dir(&scratch.path)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    files.sort_unstable();
    Pass {
        exclusions,
        open: std::fs::read_to_string(&open_path).unwrap(),
        files,
    }
}

fn stats(count: u64, excluded: u64) -> BucketStats {
    vec![("count", count), ("excluded", excluded)]
}

#[test]
fn two_feature_types_on_one_location_allele_and_feature_are_two_pairs_not_four() {
    let (location, allele, feature) = ("21:100", "A", "ENST1");
    let key = |feature_type| [location, allele, feature, feature_type];
    let outcome = pass(&[
        row("missing_in_rust", key("Transcript"), "intron_variant"),
        row(
            "missing_in_rust",
            key("RegulatoryFeature"),
            "intron_variant",
        ),
        row(
            "extra_in_rust",
            key("Transcript"),
            "intron_variant,splice_region_variant",
        ),
        row(
            "extra_in_rust",
            key("RegulatoryFeature"),
            "intron_variant,splice_polypyrimidine_tract_variant",
        ),
    ]);
    assert_eq!(
        outcome.bucket("intronic_overcall_swap"),
        stats(2, 1),
        "both pairs are counted and the splice-region one alone is excluded"
    );
    assert_eq!(
        outcome.excluded(),
        (1, 1),
        "a key without the feature type would pair the four rows crosswise and exclude two"
    );
}

#[test]
fn one_excluded_pair_takes_one_row_from_each_side() {
    let outcome = pass(&[
        transcript_row("missing_in_rust", "21:200", "A", "ENST2", "intron_variant"),
        transcript_row(
            "extra_in_rust",
            "21:200",
            "A",
            "ENST2",
            "intron_variant,splice_region_variant",
        ),
    ]);
    assert_eq!(outcome.excluded(), (1, 1));
}

/// Holds the one pair of `perl_csq` and `rust_csq` on a transcript to be counted in `bucket`
/// and to be excluded or not, the excluded pairs of the buckets summing to the exclusions.
fn assert_counted(perl_csq: &str, rust_csq: &str, bucket: &str, excluded: bool) {
    let outcome = pass(&[
        transcript_row("missing_in_rust", "21:100", "A", "ENST1", perl_csq),
        transcript_row("extra_in_rust", "21:100", "A", "ENST1", rust_csq),
    ]);
    let n = u64::from(excluded);
    assert_eq!(
        outcome.bucket(bucket),
        stats(1, n),
        "{perl_csq} -> {rust_csq} in {bucket}"
    );
    assert_eq!(outcome.excluded(), (n, n));
    assert_eq!(outcome.attributed(), n);
}

#[test]
fn a_polypyrimidine_term_added_to_an_intronic_set_is_counted_and_not_excluded() {
    assert_counted(
        "intron_variant",
        "intron_variant,splice_polypyrimidine_tract_variant",
        "intronic_overcall_swap",
        false,
    );
}

#[test]
fn a_polypyrimidine_term_added_beside_a_utr_term_is_counted_and_not_excluded() {
    assert_counted(
        "5_prime_UTR_variant,intron_variant",
        "5_prime_UTR_variant,intron_variant,splice_polypyrimidine_tract_variant",
        "ppt_overcall_swap",
        false,
    );
}

#[test]
fn splice_terms_added_to_a_frameshift_set_are_counted_and_not_excluded() {
    assert_counted(
        "frameshift_variant,stop_gained",
        "frameshift_variant,intron_variant,splice_acceptor_variant,stop_gained",
        "frameshift_splice_cascade_swap",
        false,
    );
}

#[test]
fn a_splice_sub_term_added_to_a_coding_set_is_counted_and_not_excluded() {
    assert_counted(
        "missense_variant",
        "missense_variant,splice_donor_5th_base_variant",
        "splice_lastwrite_swap",
        false,
    );
}

#[test]
fn splice_region_added_to_a_coding_set_is_counted_and_excluded() {
    assert_counted(
        "missense_variant",
        "missense_variant,splice_region_variant",
        "splice_lastwrite_swap",
        true,
    );
}

#[test]
fn the_start_co_emission_is_counted_and_excluded() {
    assert_counted(
        "start_lost,start_retained_variant",
        "start_retained_variant",
        "start_cooccurrence_swap",
        true,
    );
}

/// A report of two excluded pairs, one pair that is counted and stays, and one row a side
/// that has no partner.
fn mixed_report() -> Vec<String> {
    let missing = |location, allele, feature, csq| {
        transcript_row("missing_in_rust", location, allele, feature, csq)
    };
    let extra = |location, allele, feature, csq| {
        transcript_row("extra_in_rust", location, allele, feature, csq)
    };
    vec![
        missing("21:100", "A", "ENST1", "missense_variant"),
        extra(
            "21:100",
            "A",
            "ENST1",
            "missense_variant,splice_region_variant",
        ),
        missing("21:200", "-", "ENST2", "start_lost,start_retained_variant"),
        extra("21:200", "-", "ENST2", "start_retained_variant"),
        missing("21:300", "T", "ENST3", "intron_variant"),
        extra(
            "21:300",
            "T",
            "ENST3",
            "intron_variant,splice_polypyrimidine_tract_variant",
        ),
        missing("21:400", "G", "ENST4", "downstream_gene_variant"),
        extra("21:500", "C", "ENST5", "upstream_gene_variant"),
    ]
}

#[test]
fn the_open_rows_are_the_pairs_not_excluded_and_the_rows_without_a_partner() {
    let rows = mixed_report();
    let outcome = pass(&rows);
    assert_eq!(outcome.excluded(), (2, 2));
    let mut expected = vec![HEADER];
    expected.extend(rows[4..].iter().map(String::as_str));
    assert_eq!(outcome.open_lines(), expected);
}

#[test]
fn the_open_rows_number_the_adjusted_one_sided_total() {
    let rows = mixed_report();
    let outcome = pass(&rows);
    let on_side = |source: &str| {
        let marked = rows.iter().filter(|row| row.contains(source));
        marked.count() as u64
    };
    let (perl, rust) = (on_side("\tmissing_in_rust\t"), on_side("\textra_in_rust\t"));
    let (excluded_perl, excluded_rust) = outcome.excluded();
    assert_eq!(
        outcome.open_lines().len() as u64 - 1,
        (perl - excluded_perl) + (rust - excluded_rust)
    );
}

#[test]
fn the_open_rows_keep_the_order_and_the_bytes_of_the_report() {
    let mut rows = mixed_report();
    rows.reverse();
    let outcome = pass(&rows);
    let mut report = vec![HEADER];
    report.extend(rows.iter().map(String::as_str));
    let positions: Vec<usize> = outcome
        .open_lines()
        .iter()
        .map(|line| report.iter().position(|row| row == line).expect("a row"))
        .collect();
    assert!(positions.is_sorted(), "{positions:?}");
    assert_eq!(positions.len(), 5);
}

#[test]
fn the_pass_writes_the_open_file_it_is_given_and_no_other() {
    let outcome = pass(&mixed_report());
    assert_eq!(outcome.files, ["discordant.tsv", "discordant_open.tsv"]);
    assert_eq!(outcome.excluded(), (2, 2));
    assert_eq!(outcome.attributed(), 2);
}

#[test]
fn a_report_whose_every_row_is_excluded_leaves_the_header_alone() {
    let outcome = pass(&mixed_report()[..4]);
    assert_eq!(outcome.open, format!("{HEADER}\n"));
}

#[test]
fn a_report_of_a_header_alone_gives_an_open_file_of_the_header_alone() {
    let outcome = pass(&[]);
    assert_eq!(outcome.excluded(), (0, 0));
    assert_eq!(outcome.open, format!("{HEADER}\n"));
}
