// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `[Position in protein]` chart of the `--stats_text` report on the
//! `GRCh37-hgvs` golden corpus: the chart counts the position the
//! `Protein_position` column prints, so it is the same with `--hgvs` and
//! without it, while Ensembl VEP's chart of the corpus's `--hgvs` run counts
//! the shifted positions (`docs/divergences/stats-protein-position.md`).

use std::path::{Path, PathBuf};

mod common;

/// The ten deciles `00-10%` to `90-100%` of the chart on the corpus input, the
/// chart Ensembl VEP 116.3 writes for the run without `--hgvs`.
const DECILES: [(&str, u32); 10] = [
    ("00-10%", 171),
    ("10-20%", 75),
    ("20-30%", 79),
    ("30-40%", 95),
    ("40-50%", 90),
    ("50-60%", 129),
    ("60-70%", 65),
    ("70-80%", 78),
    ("80-90%", 113),
    ("90-100%", 158),
];

/// The rows the chart counts, the sum of its ten deciles.
const ROWS: u32 = 1_053;

/// The chart of the corpus's committed Ensembl VEP 116.3 report, a `--hgvs`
/// run: the shifted positions, and the rows whose shifted start maps to no
/// coding position left out.
const ENSEMBL_HGVS_DECILES: [(&str, u32); 10] = [
    ("00-10%", 166),
    ("10-20%", 76),
    ("20-30%", 79),
    ("30-40%", 95),
    ("40-50%", 87),
    ("50-60%", 126),
    ("60-70%", 66),
    ("70-80%", 80),
    ("80-90%", 112),
    ("90-100%", 145),
];

/// The rows Ensembl VEP's `--hgvs` chart counts.
const ENSEMBL_HGVS_ROWS: u32 = 1_032;

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37-hgvs")
}

/// The `[Position in protein]` section of a statistics report as
/// `(decile, count)` rows.
fn position_in_protein(report: &str) -> Vec<(String, u32)> {
    let section = report
        .split("\n[Position in protein]\n")
        .nth(1)
        .expect("a [Position in protein] section");
    section
        .lines()
        .take_while(|line| !line.is_empty() && !line.starts_with('['))
        .map(|line| {
            let (decile, count) = line.split_once('\t').unwrap();
            (decile.to_string(), count.parse().unwrap())
        })
        .collect()
}

fn rows(chart: &[(&str, u32)]) -> Vec<(String, u32)> {
    chart
        .iter()
        .map(|(decile, count)| (decile.to_string(), *count))
        .collect()
}

fn total(chart: &[(String, u32)]) -> u32 {
    chart.iter().map(|(_, count)| count).sum()
}

/// Runs the corpus input with `--stats_text` and `extra`, and returns the
/// report's `[Position in protein]` chart.
fn vep_rs_chart(dir: &Path, extra: &[&str]) -> Vec<(String, u32)> {
    let corpus = corpus_dir();
    let output = dir.join("out.txt");
    let result = common::isolated_vep(dir)
        .arg("-i")
        .arg(corpus.join("variants.vcf"))
        .arg("-o")
        .arg(&output)
        .arg("--json_cache")
        .arg(corpus.join("json_cache"))
        .arg("--fasta")
        .arg(corpus.join("reference.fa.gz"))
        .args([
            "--offline",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh37",
        ])
        .args(["--stats_text", "--force_overwrite", "--quiet"])
        .args(extra)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{extra:?}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report = std::fs::read_to_string(dir.join("out.txt_summary.txt")).unwrap();
    position_in_protein(&report)
}

/// vep-rs's chart is the same ten rows with and without `--hgvs`, and the rows
/// it counts sum to [`ROWS`] in both runs.
#[test]
fn position_in_protein_chart_is_the_same_with_and_without_hgvs() {
    let dir = tempfile::tempdir().unwrap();
    for extra in [&[][..], &["--hgvs"][..]] {
        let chart = vep_rs_chart(dir.path(), extra);
        assert_eq!(chart, rows(&DECILES), "{extra:?}");
        assert_eq!(total(&chart), ROWS, "{extra:?}");
    }
}

/// The corpus's committed Ensembl VEP report, written with `--hgvs`, carries
/// the shifted chart, 21 rows short of vep-rs's.
#[test]
fn committed_ensembl_vep_hgvs_report_counts_the_shifted_positions() {
    let report =
        std::fs::read_to_string(corpus_dir().join("expected/default.summary.log")).unwrap();
    assert!(
        report.contains("--hgvs"),
        "the committed report is a --hgvs run"
    );
    let chart = position_in_protein(&report);
    assert_eq!(chart, rows(&ENSEMBL_HGVS_DECILES));
    assert_eq!(total(&chart), ENSEMBL_HGVS_ROWS);
    assert_eq!(ROWS - ENSEMBL_HGVS_ROWS, 21);
}
