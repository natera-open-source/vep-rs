// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `Novel / existing variants` row of the `--stats_text` report on the
//! `GRCh38-freqfilter` golden corpus: a processed variant is existing when a
//! known variant of the cache is co-located with it, with `--check_frequency`
//! as with `--check_existing` alone, while Ensembl VEP's report of a
//! `--check_frequency` run counts every processed variant as existing
//! (`docs/divergences/stats-existing-variants-check-frequency.md`).

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The corpus's committed `freqfilter` flags: an allele whose `AF` in the
/// cache is above 0.01 is dropped.
const FREQUENCY_CHECK: [&str; 9] = [
    "--check_frequency",
    "--freq_pop",
    "AF",
    "--freq_freq",
    "0.01",
    "--freq_gt_lt",
    "gt",
    "--freq_filter",
    "exclude",
];

/// The position of the corpus record with no known variant: the cache's
/// known-variant file holds no row within one base of it, the window Ensembl
/// VEP reads for a variant (`VariationTabix.pm:229`).
const NOVEL: &str = "33030010";

/// The position of a corpus record with a co-located known variant,
/// `rs1254937786`, which the frequency check keeps.
const KNOWN: &str = "30107044";

/// The committed GRCh38-freqfilter corpus.
fn corpus() -> common::Corpus {
    common::corpora()
        .into_iter()
        .find(|corpus| corpus.name == "GRCh38-freqfilter")
        .expect("the GRCh38-freqfilter corpus is committed")
}

/// The corpus's input cut to its header and the records at `positions`,
/// written under `dir`.
fn records_at(dir: &Path, positions: &[&str]) -> PathBuf {
    let text = std::fs::read_to_string(corpus().dir.join("variants.vcf")).unwrap();
    let kept: Vec<&str> = text
        .lines()
        .filter(|line| {
            line.starts_with('#')
                || line
                    .split('\t')
                    .nth(1)
                    .is_some_and(|position| positions.contains(&position))
        })
        .collect();
    let records = kept.iter().filter(|line| !line.starts_with('#')).count();
    assert_eq!(records, positions.len(), "{positions:?}");
    let input = dir.join("records.vcf");
    std::fs::write(&input, format!("{}\n", kept.join("\n"))).unwrap();
    input
}

/// The rows of the `[General statistics]` section of a `--stats_text` run on
/// `input` against the corpus's cache with `flags`, in the report's order.
fn general_statistics(dir: &Path, input: &Path, flags: &[&str]) -> Vec<String> {
    let corpus = corpus();
    let out = dir.join("out.txt");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .args(corpus.cache_args())
        .args([
            "--assembly",
            &corpus.assembly,
            "--force_overwrite",
            "--stats_text",
            "--quiet",
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
    std::fs::read_to_string(dir.join("out.txt_summary.txt"))
        .unwrap()
        .lines()
        .skip_while(|line| *line != "[General statistics]")
        .skip(1)
        .take_while(|line| !line.is_empty())
        .map(str::to_string)
        .collect()
}

/// Of the two records, one has a known variant at its position and one has
/// none, so the row is `1 (50.0) / 1 (50.0)`, under the frequency check as
/// under `--check_existing` alone, which is also Ensembl VEP 116.3's row for
/// the run without the check. Under the check Ensembl VEP 116.3 prints
/// `0 (0.0) / 2 (100.0)`: its frequency routine loops over each variant's list
/// of known variants and so creates an empty list where there was none
/// (`BaseCacheVariation.pm:195`), and the report counts a variant as existing
/// when the list is there, empty or not (`Stats.pm:244`).
#[test]
fn a_frequency_check_leaves_the_count_of_existing_variants_as_it_is() {
    let dir = tempfile::tempdir().unwrap();
    let input = records_at(dir.path(), &[KNOWN, NOVEL]);
    let expected = [
        "Lines of input read\t2",
        "Variants processed\t2",
        "Variants filtered out\t0",
        "Novel / existing variants\t1 (50.0) / 1 (50.0)",
        "Overlapped genes\t3",
        "Overlapped transcripts\t20",
        "Overlapped regulatory features\t-",
    ];
    assert_eq!(
        general_statistics(dir.path(), &input, &FREQUENCY_CHECK),
        expected
    );
    assert_eq!(
        general_statistics(dir.path(), &input, &["--check_existing"]),
        expected
    );
}

/// A run whose one record has no known variant counts no existing variant,
/// and the row is `-`, as `Stats.pm:727-734` writes it for a count never made;
/// Ensembl VEP 116.3 prints `0 (0.0) / 1 (100.0)` under the check.
#[test]
fn a_frequency_check_over_a_record_with_no_known_variant_counts_none_existing() {
    let dir = tempfile::tempdir().unwrap();
    let input = records_at(dir.path(), &[NOVEL]);
    assert_eq!(
        general_statistics(dir.path(), &input, &FREQUENCY_CHECK),
        [
            "Lines of input read\t1",
            "Variants processed\t1",
            "Variants filtered out\t0",
            "Novel / existing variants\t-",
            "Overlapped genes\t3",
            "Overlapped transcripts\t20",
            "Overlapped regulatory features\t-",
        ]
    );
}

/// The corpus's fifteen records under its own flags: the check drops four,
/// ten of the eleven it keeps carry an identifier in both engines' rows and
/// the record at [`NOVEL`] carries none, so the row is `1 (9.1) / 10 (90.9)`
/// where Ensembl VEP 116.3 prints `0 (0.0) / 11 (100.0)`. Under
/// `--check_existing` alone the fifteen are fourteen existing and one novel in
/// both engines. Under an `include` filter at 0.0001 a record with no known
/// variant has no frequency and is dropped, so the twelve records kept are all
/// existing, in both engines too.
#[test]
fn the_corpus_run_counts_its_one_record_with_no_known_variant_as_novel() {
    let dir = tempfile::tempdir().unwrap();
    let input = corpus().dir.join("variants.vcf");
    assert_eq!(
        general_statistics(dir.path(), &input, &FREQUENCY_CHECK),
        [
            "Lines of input read\t15",
            "Variants processed\t11",
            "Variants filtered out\t4",
            "Novel / existing variants\t1 (9.1) / 10 (90.9)",
            "Overlapped genes\t19",
            "Overlapped transcripts\t216",
            "Overlapped regulatory features\t-",
        ]
    );
    assert_eq!(
        general_statistics(dir.path(), &input, &["--check_existing"]),
        [
            "Lines of input read\t15",
            "Variants processed\t15",
            "Variants filtered out\t0",
            "Novel / existing variants\t1 (6.7) / 14 (93.3)",
            "Overlapped genes\t23",
            "Overlapped transcripts\t254",
            "Overlapped regulatory features\t-",
        ]
    );
    let include = [
        "--check_frequency",
        "--freq_pop",
        "AF",
        "--freq_freq",
        "0.0001",
        "--freq_gt_lt",
        "gt",
        "--freq_filter",
        "include",
    ];
    assert_eq!(
        general_statistics(dir.path(), &input, &include),
        [
            "Lines of input read\t15",
            "Variants processed\t12",
            "Variants filtered out\t3",
            "Novel / existing variants\t0 (0.0) / 12 (100.0)",
            "Overlapped genes\t18",
            "Overlapped transcripts\t195",
            "Overlapped regulatory features\t-",
        ]
    );
}
