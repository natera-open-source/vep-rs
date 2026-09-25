// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Format-parity tests: the `--tab`, `--vcf` and `--json` outputs (and the
//! Parquet output when the DuckDB CLI is installed) on each golden corpus are
//! compared with the Ensembl VEP files committed beside it. Header lines VEP
//! fills with run-specific values (time, paths, command line, Perl API
//! component versions) compare by shape; every other header line and every data
//! value compares exactly, except the consequence keys the corpus manifest
//! documents as divergent.

mod common;

use std::path::Path;
use std::process::Command;

use common::{
    compare_entries, compare_tab_headers, compare_vcf_headers, corpora, input_records, parse_json,
    parse_tab, parse_vcf, read_expected, run_vep, Corpus, Documented,
};

fn check(corpus: &Corpus, label: &str, header: String, body: String, failures: &mut Vec<String>) {
    if !header.is_empty() || !body.is_empty() {
        failures.push(format!(
            "== corpus {}/{} ({label})\n{header}{body}",
            corpus.release, corpus.name
        ));
    }
}

#[test]
fn tab_output_matches_vep_on_every_corpus() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        let expected = read_expected(&corpus.dir, "tab.txt");
        let actual = run_vep(corpus, "tab", tmp.path(), &[]);
        let documented = Documented::from_manifest(&corpus.manifest);
        check(
            corpus,
            "--tab",
            compare_tab_headers(&expected, &actual),
            compare_entries(&parse_tab(&expected), &parse_tab(&actual), &documented),
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn vcf_output_matches_vep_on_every_corpus() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        let expected = read_expected(&corpus.dir, "vcf.vcf");
        let actual = run_vep(corpus, "vcf", tmp.path(), &[]);
        let documented = Documented::from_manifest(&corpus.manifest);
        check(
            corpus,
            "--vcf",
            compare_vcf_headers(&expected, &actual),
            compare_entries(&parse_vcf(&expected), &parse_vcf(&actual), &documented),
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn json_output_matches_vep_on_every_corpus() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        let expected = read_expected(&corpus.dir, "json.jsonl");
        let actual = run_vep(corpus, "json", tmp.path(), &[]);
        let documented = Documented::from_manifest(&corpus.manifest);
        let inputs = input_records(corpus);
        check(
            corpus,
            "--json",
            String::new(),
            compare_entries(
                &parse_json(&expected, &inputs),
                &parse_json(&actual, &inputs),
                &documented,
            ),
            &mut failures,
        );
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

fn duckdb_available() -> bool {
    Command::new("duckdb")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// The Parquet output, read back through the adapter, holds every `--tab` row
/// (Parquet is sorted by variant key, so rows compare as a multiset), in both
/// shapes; every scalar column of every row group carries a dictionary, a Bloom
/// filter and statistics, and the footer carries the run's key-value metadata.
#[test]
fn parquet_round_trips_to_the_tab_output() {
    if !duckdb_available() {
        eprintln!("skipping: duckdb CLI not installed");
        return;
    }
    let adapter =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../scripts/adapters/parquet_to_vep_tab.py");
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        let tab = run_vep(corpus, "tab", tmp.path(), &[]);
        let mut expected_rows: Vec<&str> = tab.lines().filter(|l| !l.starts_with("##")).collect();
        expected_rows.sort_unstable();
        for shape in ["flat", "nested"] {
            let parquet_dir = tmp.path().join(format!(
                "{}_{}_{shape}.parquet",
                corpus.release, corpus.name
            ));
            let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
            cmd.arg("-i")
                .arg(corpus.dir.join("variants.vcf"))
                .arg("-o")
                .arg(&parquet_dir)
                .arg("--offline")
                .arg("--json_cache")
                .arg(corpus.dir.join("json_cache"))
                .args(["--species", "homo_sapiens", "--assembly", &corpus.name])
                .args([
                    "--buffer_size",
                    "5000",
                    "--force_overwrite",
                    "--no_stats",
                    "--quiet",
                ])
                .args(["--output_format", "parquet", "--parquet_shape", shape]);
            let out = cmd.output().unwrap();
            assert!(
                out.status.success(),
                "vep --output_format parquet ({shape}) failed:\n{}",
                String::from_utf8_lossy(&out.stderr)
            );
            let back = Command::new("python3")
                .arg(&adapter)
                .arg(&parquet_dir)
                .output()
                .unwrap();
            assert!(
                back.status.success(),
                "adapter failed:\n{}",
                String::from_utf8_lossy(&back.stderr)
            );
            let round_trip = String::from_utf8(back.stdout).unwrap();
            let mut actual_rows: Vec<&str> = round_trip
                .lines()
                .filter(|l| !l.starts_with("##"))
                .collect();
            actual_rows.sort_unstable();
            if expected_rows != actual_rows {
                let first = expected_rows
                    .iter()
                    .zip(actual_rows.iter())
                    .position(|(a, b)| a != b)
                    .unwrap_or(expected_rows.len().min(actual_rows.len()));
                failures.push(format!(
                    "== corpus {}/{} (parquet {shape} round trip): {} tab rows vs {} round-trip rows; first difference at sorted line {first}:\n  tab:     {:?}\n  parquet: {:?}",
                    corpus.release,
                    corpus.name,
                    expected_rows.len(),
                    actual_rows.len(),
                    expected_rows.get(first),
                    actual_rows.get(first)
                ));
            }
            if shape == "flat" {
                if let Some(report) = parquet_metadata_report(&parquet_dir) {
                    failures.push(format!(
                        "== corpus {}/{} (parquet metadata): {report}",
                        corpus.release, corpus.name
                    ));
                }
            }
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// Several row groups per partition, each with a Bloom filter and statistics on
/// every scalar column, and the Bloom filters prune: a `Feature` absent from the
/// file but inside the column's value range (so min/max statistics cannot rule it
/// out) is excluded by every row group's filter, and a `Feature` present in one
/// row is excluded by every row group but its own.
#[test]
fn parquet_row_groups_prune_on_bloom_filters() {
    if !duckdb_available() {
        eprintln!("skipping: duckdb CLI not installed");
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    let parquet_dir = tmp.path().join("rowgroups.parquet");
    let out = Command::new(env!("CARGO_BIN_EXE_vep"))
        .arg("-i")
        .arg(corpus.dir.join("variants.vcf"))
        .arg("-o")
        .arg(&parquet_dir)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.dir.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", &corpus.name])
        .args([
            "--buffer_size",
            "5000",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
        ])
        .args(["--output_format", "parquet", "--parquet_shape", "flat"])
        .args(["--parquet_row_group_size", "2048"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "vep --output_format parquet failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    if let Some(report) = parquet_metadata_report(&parquet_dir) {
        panic!("parquet metadata at 2,048-row groups: {report}");
    }
    // The partition with the most rows is the one probed.
    let glob = format!("{}/**/*.parquet", parquet_dir.display()).replace('\'', "''");
    let file = duckdb_csv(&format!(
        "SELECT file_name FROM parquet_file_metadata('{glob}') ORDER BY num_rows DESC LIMIT 1"
    ))
    .trim()
    .replace('\'', "''");
    let row_groups: usize = duckdb_csv(&format!(
        "SELECT num_row_groups FROM parquet_file_metadata('{file}')"
    ))
    .trim()
    .parse()
    .unwrap();
    assert!(
        row_groups >= 3,
        "expected at least three row groups in {file}, found {row_groups}"
    );
    let bounds = duckdb_csv(&format!(
        "SELECT min(Feature), max(Feature) FROM read_parquet('{file}')"
    ));
    let (min, max) = bounds.trim().split_once(',').unwrap();
    // Appending a character to the smallest id makes a value no row holds that
    // still sorts inside [min, max].
    let absent = format!("{min}0");
    assert!(min < absent.as_str() && absent.as_str() < max);
    let excluded_absent = bloom_excluded_row_groups(&file, &absent);
    assert_eq!(
        excluded_absent, row_groups,
        "the absent Feature {absent} was excluded by {excluded_absent} of {row_groups} row groups"
    );
    // A Feature that occurs on exactly one row sits in exactly one row group.
    let single = duckdb_csv(&format!(
        "SELECT Feature FROM read_parquet('{file}') WHERE Feature IS NOT NULL \
         GROUP BY Feature HAVING count(*) = 1 ORDER BY Feature LIMIT 1"
    ))
    .trim()
    .to_string();
    assert!(!single.is_empty(), "no single-row Feature in {file}");
    let excluded_present = bloom_excluded_row_groups(&file, &single);
    assert_eq!(
        excluded_present,
        row_groups - 1,
        "the single-row Feature {single} was excluded by {excluded_present} of {row_groups} row groups"
    );
}

fn duckdb_csv(sql: &str) -> String {
    let out = Command::new("duckdb")
        .args(["-csv", "-noheader", "-c", sql])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "duckdb failed:\n{}\nsql: {sql}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

/// Row groups of one Parquet file whose `Feature` Bloom filter rules out `value`.
fn bloom_excluded_row_groups(file: &str, value: &str) -> usize {
    duckdb_csv(&format!(
        "SELECT count(*) FILTER (WHERE bloom_filter_excludes) \
         FROM parquet_bloom_probe('{file}', 'Feature', '{}')",
        value.replace('\'', "''")
    ))
    .trim()
    .parse()
    .unwrap()
}

/// Checks the flat Parquet directory's metadata; `None` when everything holds.
fn parquet_metadata_report(dir: &Path) -> Option<String> {
    let glob = format!("{}/**/*.parquet", dir.display()).replace('\'', "''");
    let sql = format!(
        "SELECT path_in_schema, encodings, bloom_filter_offset IS NOT NULL, stats_min IS NOT NULL, \
         num_values, stats_null_count FROM parquet_metadata('{glob}')"
    );
    let out = Command::new("duckdb")
        .args(["-csv", "-noheader", "-c", &sql])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "parquet_metadata failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let mut problems = Vec::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() < 6 {
            continue;
        }
        let (column, encodings, bloom, stats, values, nulls) =
            (cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]);
        // A column that is NULL throughout a row group has no dictionary to
        // filter on; every other column must carry the full read-side kit.
        let all_null = values == "0" || values == nulls;
        let complete = encodings.contains("DICTIONARY") && bloom == "true" && stats == "true";
        if !all_null && !complete {
            problems.push(format!(
                "{column}: encodings={encodings} bloom={bloom} stats={stats}"
            ));
        }
    }
    let kv_sql = format!("SELECT key FROM parquet_kv_metadata('{glob}')");
    let kv = Command::new("duckdb")
        .args(["-csv", "-noheader", "-c", &kv_sql])
        .output()
        .unwrap();
    let keys = String::from_utf8_lossy(&kv.stdout);
    for wanted in ["vep_rs_version", "assembly", "command_line"] {
        if !keys.lines().any(|k| k == wanted) {
            problems.push(format!("footer lacks key-value metadata {wanted}"));
        }
    }
    if problems.is_empty() {
        None
    } else {
        problems.sort();
        problems.dedup();
        Some(problems.join("; "))
    }
}
