// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `convert-plugin-data` holds its layout flags to the table's tabix index.
//!
//! Each test stages a committed table of made-up scores with its index and
//! runs the binary on it. The index of the single-position CADD table
//! declares one-based starts and no end column; the index of the CADD span
//! table declares zero-based starts and an end column, the third one; the
//! index of the SpliceAI VCF is of tabix's `vcf` preset, which declares
//! one-based starts and no end column.
//!
//! The flags name what the index declares, or single one-based positions.
//! Any other `--end-col` given without `--ref-col` is refused before the
//! table is read, and so is any other `--zero-based` unless it is given with
//! `--end-col` and `--ref-col` under an index that declares zero-based
//! starts: the run ends with exit code 1 on one
//! error line, which names the flags that match the index, and leaves no
//! file. Flags that name single one-based positions where the index declares
//! more draw one warning line, which names the flags that match the index,
//! and the conversion goes on to its store. Flags that match the index draw
//! neither.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// The committed single-position table, under `tests/golden/116`.
const SINGLE_POSITIONS: &str =
    "GRCh38-plugin-cadd-fixture/plugin_data/whole_genome_SNVs_synthetic.tsv.gz";

/// The committed span table, under `tests/golden/116`.
const SPANS: &str = "GRCh38-plugin-cadd-sv-fixture/plugin_data/1000G_phase3_SVs_synthetic.tsv.gz";

/// The committed VCF, under `tests/golden/116`, indexed with the `vcf` preset.
const VCF: &str =
    "GRCh38-plugin-spliceai-fixture/plugin_data/spliceai_scores.raw.snv.hg38.synthetic.vcf.gz";

/// What a refusal says after the flags of an index of single one-based
/// positions.
const AS_DECLARED: &str = "--chr-col 0 --start-col 1, the layout the index declares";

/// Copies a committed table and its `.tbi` into `dir`; returns the copy.
fn staged(dir: &Path, table: &str) -> PathBuf {
    let source = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(table);
    let name = source.file_name().unwrap().to_string_lossy().into_owned();
    let copy = dir.join(&name);
    std::fs::copy(&source, &copy).unwrap();
    std::fs::copy(
        source.with_file_name(format!("{name}.tbi")),
        dir.join(format!("{name}.tbi")),
    )
    .unwrap();
    copy
}

/// The store the runs of one test write.
fn store(dir: &Path) -> PathBuf {
    dir.join("scores.vpd")
}

/// `vep-cache-converter convert-plugin-data` on `table` with `flags`, its log
/// lines without colour codes.
fn convert(table: &Path, flags: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"))
        .env("NO_COLOR", "1")
        .arg("convert-plugin-data")
        .arg("--input")
        .arg(table)
        .arg("--output")
        .arg(store(table.parent().unwrap()))
        .args(flags)
        .output()
        .unwrap()
}

/// Standard error of `output` as text.
fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// The warning lines of a run's log, each without its timestamp and the
/// space after it.
fn warning_lines(output: &Output) -> Vec<String> {
    stderr(output)
        .lines()
        .filter_map(|line| line.split_once(' ').map(|(_, rest)| rest.to_string()))
        .filter(|rest| rest.trim_start().starts_with("WARN "))
        .collect()
}

/// The names of the files in `dir`, sorted.
fn files(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// Asserts that a run was refused before it read its table: exit code 1, the
/// one error line, which names the flags given as `named` and ends its
/// sentence on `advice`, and no file beside those `dir` held before.
fn assert_refused(table: &Path, flags: &[&str], named: &str, advice: &str) {
    let dir = table.parent().unwrap();
    let before = files(dir);
    let output = convert(table, flags);
    assert_eq!(output.status.code(), Some(1), "{flags:?}");
    assert_eq!(
        stderr(&output),
        format!(
            "ERROR: [cache-converter] {named} name another layout than the index of {} \
             declares; convert with {advice} (docs/plugins/data-files.md, \"The layout flags \
             of a store\")\n",
            table.display()
        ),
        "{flags:?}"
    );
    assert_eq!(files(dir), before, "{flags:?} left a file");
}

/// Asserts that a run converted its table: exit code 0, both files of the
/// store, and the warning lines `expected`.
fn assert_converted(table: &Path, flags: &[&str], expected: &[String]) {
    let dir = table.parent().unwrap();
    let output = convert(table, flags);
    assert_eq!(
        output.status.code(),
        Some(0),
        "{flags:?}: {}",
        stderr(&output)
    );
    assert_eq!(warning_lines(&output), expected, "{flags:?}");
    let written = store(dir);
    assert!(
        written.is_file() && written.with_extension("vpdi").is_file(),
        "{flags:?} wrote no store"
    );
}

/// On the single-position table, `--zero-based` and an `--end-col` other
/// than the start column each name what the index does not declare, alone and
/// together, and each run is refused. An `--end-col` naming the start column
/// is no end column and converts, as the default flags do.
#[test]
fn zero_based_starts_or_an_end_column_the_index_does_not_declare_are_refused() {
    let dir = tempfile::tempdir().unwrap();
    let table = staged(dir.path(), SINGLE_POSITIONS);
    assert_refused(
        &table,
        &["--zero-based"],
        "--chr-col 0 --start-col 1 --zero-based",
        AS_DECLARED,
    );
    assert_refused(
        &table,
        &["--end-col", "2"],
        "--chr-col 0 --start-col 1 --end-col 2",
        AS_DECLARED,
    );
    assert_refused(
        &table,
        &["--end-col", "2", "--zero-based"],
        "--chr-col 0 --start-col 1 --end-col 2 --zero-based",
        AS_DECLARED,
    );
    assert_converted(&table, &["--end-col", "1"], &[]);
    assert_converted(&table, &[], &[]);
}

/// On the span table every `--end-col` or `--zero-based` given without
/// `--ref-col` but the index's own pair is refused: another end column, with
/// and without `--zero-based`, the
/// index's end column alone and its zero-based starts alone. The refusal
/// names the index's own flags.
#[test]
fn an_end_column_or_starts_other_than_the_index_declares_are_refused() {
    let dir = tempfile::tempdir().unwrap();
    let table = staged(dir.path(), SPANS);
    let advice =
        "--chr-col 0 --start-col 1 --end-col 2 --zero-based, the layout the index declares";
    for (flags, named) in [
        (
            &["--end-col", "3"][..],
            "--chr-col 0 --start-col 1 --end-col 3",
        ),
        (
            &["--end-col", "3", "--zero-based"][..],
            "--chr-col 0 --start-col 1 --end-col 3 --zero-based",
        ),
        (
            &["--end-col", "2"][..],
            "--chr-col 0 --start-col 1 --end-col 2",
        ),
        (
            &["--zero-based"][..],
            "--chr-col 0 --start-col 1 --zero-based",
        ),
    ] {
        assert_refused(&table, flags, named, advice);
    }
}

/// On the span table the default flags, single one-based positions, convert
/// with one warning line, and the flags the line names convert with none.
#[test]
fn single_positions_under_an_index_that_declares_more_convert_with_one_warning() {
    let dir = tempfile::tempdir().unwrap();
    let table = staged(dir.path(), SPANS);
    let declared = "--chr-col 0 --start-col 1 --end-col 2 --zero-based";
    let warning = format!(
        " WARN vep_cache_converter: the index of {} declares {declared} and the flags name \
         --chr-col 0 --start-col 1: CADD, which takes its columns from the index, refuses this \
         store; convert with {declared} for CADD",
        table.display()
    );
    assert_converted(&table, &[], &[warning]);
    assert_converted(&table, &["--end-col", "2", "--zero-based"], &[]);
}

/// An index of the `vcf` preset declares one-based starts and no end column.
/// The VCF converts with the default flags and with `--end-col` beside
/// `--ref-col`, its `INFO` and reference allele columns, each time with no
/// warning line. An `--end-col` without `--ref-col` is refused, and the
/// refusal names both ways out, the second one without `--zero-based` where
/// the flags name it; `--zero-based` is refused, with and without the two
/// columns.
#[test]
fn under_a_vcf_index_the_flags_name_single_positions_or_the_info_and_reference_columns() {
    let dir = tempfile::tempdir().unwrap();
    let table = staged(dir.path(), VCF);
    assert_refused(
        &table,
        &["--end-col", "7"],
        "--chr-col 0 --start-col 1 --end-col 7",
        "--chr-col 0 --start-col 1, the layout the index declares, or with --ref-col beside \
         --end-col, a VCF's INFO and reference allele columns",
    );
    assert_refused(
        &table,
        &["--end-col", "7", "--zero-based"],
        "--chr-col 0 --start-col 1 --end-col 7 --zero-based",
        "--chr-col 0 --start-col 1, the layout the index declares, or without --zero-based and \
         with --ref-col beside --end-col, a VCF's INFO and reference allele columns",
    );
    assert_refused(
        &table,
        &["--zero-based"],
        "--chr-col 0 --start-col 1 --zero-based",
        AS_DECLARED,
    );
    assert_refused(
        &table,
        &["--end-col", "7", "--ref-col", "3", "--zero-based"],
        "--chr-col 0 --start-col 1 --end-col 7 --ref-col 3 --zero-based",
        "--chr-col 0 --start-col 1 --end-col 7 --ref-col 3, the starts the index declares",
    );
    assert_converted(&table, &["--end-col", "7", "--ref-col", "3"], &[]);
    assert_converted(&table, &[], &[]);
}
