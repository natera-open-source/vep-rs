// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Format-parity tests: the `--tab`, `--vcf` and `--json` outputs (and the
//! Parquet output when the DuckDB CLI is installed) on each golden corpus, and
//! on each `flag_sets` run of its manifest that names the format, are compared
//! with the Ensembl VEP files committed beside it. Header lines VEP fills with
//! run-specific values (time, paths, command line, Perl API component versions)
//! compare by shape; every other header line and every data value compares
//! exactly, except the consequence keys the corpus manifest documents as
//! divergent. The Parquet round trip runs on the manifest's own run alone.

mod common;

use std::path::Path;
use std::process::Command;

use common::{
    compare_run, corpora, default_run, parse_tab, parse_vcf, run_output, run_vep, runs,
    vep_command, Run,
};

fn check(run: &Run, label: &str, header: String, body: String, failures: &mut Vec<String>) {
    if !header.is_empty() || !body.is_empty() {
        failures.push(format!(
            "== corpus {} ({label})\n{header}{body}",
            run.label()
        ));
    }
}

#[test]
fn tab_output_matches_vep_on_every_corpus() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        for run in runs(corpus).iter().filter(|r| r.compares("tab")) {
            let (header, body) = compare_run(run, "tab", tmp.path());
            check(run, "--tab", header, body, &mut failures);
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn vcf_output_matches_vep_on_every_corpus() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        for run in runs(corpus).iter().filter(|r| r.compares("vcf")) {
            let (header, body) = compare_run(run, "vcf", tmp.path());
            check(run, "--vcf", header, body, &mut failures);
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn json_output_matches_vep_on_every_corpus() {
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora() {
        for run in runs(corpus).iter().filter(|r| r.compares("json")) {
            let (header, body) = compare_run(run, "json", tmp.path());
            check(run, "--json", header, body, &mut failures);
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// Ensembl VEP lists `OverlapBP` and `OverlapPC` as the fields `--overlaps` adds
/// (`Constants.pm` 131), last among the flag fields of the tab columns and of the
/// VCF CSQ fields (after `FLAGS`, and after the `SYMBOL_SOURCE` and `HGNC_ID` that
/// `--vcf` switches on), and fills them on every structural variant row except a
/// breakpoint's, when the overlap is positive, with the overlap length and
/// `sprintf("%.2f")` of the percentage of the feature it covers (`OutputFactory.pm`
/// 2087-2094). On the GRCh38 corpus the 272-base deletion
/// `gnomAD-SV_v3_DEL_chr21_9b2929e7` (21:13610236-13610507) against ENST00000299443
/// carries 272 and 0.75 and against ENST00000620442 272 and 1.73, the values Ensembl
/// VEP 116.2 writes for those rows; the breakend `synth_bnd_inter_0055_A` carries
/// neither (`-` in tab, empty in VCF). Without the flag neither format names the
/// two fields.
#[test]
fn overlaps_adds_overlapbp_overlappc_columns_to_tab_and_vcf() {
    let tmp = tempfile::tempdir().unwrap();
    let corpora = corpora();
    let corpus = corpora
        .iter()
        .find(|c| c.release == "116" && c.name == "GRCh38")
        .expect("the 116/GRCh38 corpus");
    let run = default_run(corpus);
    let deletion_id = "gnomAD-SV_v3_DEL_chr21_9b2929e7";
    let breakend_id = "synth_bnd_inter_0055_A";
    let column_line = |tab: &str| {
        tab.lines()
            .find(|l| l.starts_with("#Uploaded_variation"))
            .unwrap()
            .to_string()
    };
    let csq_format = |vcf: &str| {
        vcf.lines()
            .find(|l| l.starts_with("##INFO=<ID=CSQ"))
            .unwrap()
            .to_string()
    };

    let tab = run_vep(&run, "tab", tmp.path(), &[]);
    assert!(column_line(&tab).ends_with("\tFLAGS"));
    assert!(!tab.contains("OverlapBP"));
    let vcf = run_vep(&run, "vcf", tmp.path(), &[]);
    assert!(csq_format(&vcf).ends_with("|FLAGS|SYMBOL_SOURCE|HGNC_ID\">"));
    assert!(!vcf.contains("OverlapBP"));

    let tab = run_vep(&run, "tab", tmp.path(), &["--overlaps"]);
    assert!(column_line(&tab).ends_with("\tFLAGS\tOverlapBP\tOverlapPC"));
    let rows = parse_tab(&tab);
    let tab_values = |id: &str, feature: &str| {
        let row = rows
            .iter()
            .find(|e| e.fields["Uploaded_variation"] == id && e.feature == feature)
            .unwrap_or_else(|| panic!("tab row {id} x {feature}"));
        (
            row.fields["OverlapBP"].as_str(),
            row.fields["OverlapPC"].as_str(),
        )
    };
    assert_eq!(tab_values(deletion_id, "ENST00000299443"), ("272", "0.75"));
    assert_eq!(tab_values(deletion_id, "ENST00000620442"), ("272", "1.73"));
    assert_eq!(tab_values(breakend_id, "ENST00000374832"), ("-", "-"));

    let vcf = run_vep(&run, "vcf", tmp.path(), &["--overlaps"]);
    assert!(csq_format(&vcf).ends_with("|FLAGS|SYMBOL_SOURCE|HGNC_ID|OverlapBP|OverlapPC\">"));
    let inputs = run.input_records();
    let entries = parse_vcf(&vcf, &inputs);
    let csq_values = |id: &str, feature: &str| {
        let id_column = format!("\t{id}\t");
        let record = inputs
            .iter()
            .position(|l| l.contains(&id_column))
            .unwrap_or_else(|| panic!("input record {id}"));
        let entry = entries
            .iter()
            .find(|e| e.record == Some(record) && e.feature == feature)
            .unwrap_or_else(|| panic!("CSQ entry of record {record} x {feature}"));
        (
            entry.fields["OverlapBP"].as_str(),
            entry.fields["OverlapPC"].as_str(),
        )
    };
    assert_eq!(csq_values(deletion_id, "ENST00000299443"), ("272", "0.75"));
    assert_eq!(csq_values(deletion_id, "ENST00000620442"), ("272", "1.73"));
    assert_eq!(csq_values(breakend_id, "ENST00000374832"), ("", ""));

    // The flat Parquet output takes the tab columns, so it gains the two under the flag.
    if !duckdb_available() {
        eprintln!("skipping the Parquet half: duckdb CLI not installed");
        return;
    }
    let write_parquet = |flags: &[&str], dir: &str| {
        let parquet_dir = tmp.path().join(dir);
        let mut cmd = vep_command(&run, tmp.path());
        cmd.arg("-o")
            .arg(&parquet_dir)
            .args(["--output_format", "parquet", "--parquet_shape", "flat"])
            .args(flags);
        let out = run_output(&run, &mut cmd);
        assert!(
            out.status.success(),
            "vep --output_format parquet {flags:?} failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
        format!("{}/**/*.parquet", parquet_dir.display()).replace('\'', "''")
    };
    let overlap_columns = |glob: &str| {
        duckdb_csv(&format!(
            "SELECT count(*) FROM (DESCRIBE SELECT * FROM read_parquet('{glob}')) \
             WHERE column_name IN ('OverlapBP', 'OverlapPC')"
        ))
        .trim()
        .to_string()
    };
    let plain_glob = write_parquet(&[], "plain.parquet");
    assert_eq!(overlap_columns(&plain_glob), "0");
    let overlaps_glob = write_parquet(&["--overlaps"], "overlaps.parquet");
    assert_eq!(overlap_columns(&overlaps_glob), "2");
    let values = duckdb_csv(&format!(
        "SELECT OverlapBP, OverlapPC FROM read_parquet('{overlaps_glob}') \
         WHERE Uploaded_variation = '{deletion_id}' AND Feature = 'ENST00000299443'"
    ));
    assert_eq!(values.trim(), "272,0.75");
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
        let run = default_run(corpus);
        let tab = run_vep(&run, "tab", tmp.path(), &[]);
        let mut expected_rows: Vec<&str> = tab.lines().filter(|l| !l.starts_with("##")).collect();
        expected_rows.sort_unstable();
        for shape in ["flat", "nested"] {
            let parquet_dir = tmp.path().join(format!(
                "{}_{}_{shape}.parquet",
                corpus.release, corpus.name
            ));
            let mut cmd = vep_command(&run, tmp.path());
            cmd.arg("-o").arg(&parquet_dir).args([
                "--output_format",
                "parquet",
                "--parquet_shape",
                shape,
            ]);
            let out = run_output(&run, &mut cmd);
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
    let corpora = corpora();
    let run = default_run(&corpora[0]);
    let parquet_dir = tmp.path().join("rowgroups.parquet");
    let mut cmd = vep_command(&run, tmp.path());
    cmd.arg("-o")
        .arg(&parquet_dir)
        .args(["--output_format", "parquet", "--parquet_shape", "flat"])
        .args(["--parquet_row_group_size", "2048"]);
    let out = run_output(&run, &mut cmd);
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
