// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over the fixtures of `tests/fixtures/` and holds its standard output,
//! its error line and its exit code to the committed bytes of `tests/fixtures/expected/`.
//!
//! `tests/fixtures/cases.tsv` lists the cases: a name, the exit code, how the error stream
//! is held (`exact`: byte for byte; `usage`: it opens with the usage line; `any`: not
//! read), then the arguments, each relative to the fixtures directory, where the binary
//! runs.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const THREADS_VARIABLE: &str = "VEP_PARQUET_TO_TAB_THREADS";

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn run(args: &[&str], threads: Option<&str>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-parquet-to-tab"));
    command.args(args).current_dir(fixtures());
    match threads {
        Some(threads) => command.env(THREADS_VARIABLE, threads),
        None => command.env_remove(THREADS_VARIABLE),
    };
    command.output().expect("the binary runs")
}

struct Case {
    name: String,
    code: i32,
    stderr: String,
    args: Vec<String>,
}

fn cases() -> Vec<Case> {
    let listing = fs::read_to_string(fixtures().join("cases.tsv")).expect("cases.tsv");
    listing
        .lines()
        .map(|line| {
            let mut fields = line.split('\t').map(str::to_string);
            Case {
                name: fields.next().expect("a name"),
                code: fields
                    .next()
                    .expect("an exit code")
                    .parse()
                    .expect("a number"),
                stderr: fields.next().expect("a stderr rule"),
                args: fields.collect(),
            }
        })
        .collect()
}

fn expected(name: &str, stream: &str) -> Option<Vec<u8>> {
    fs::read(fixtures().join("expected").join(format!("{name}.{stream}"))).ok()
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn hold(case: &Case, threads: Option<&str>) {
    let args: Vec<&str> = case.args.iter().map(String::as_str).collect();
    let output = run(&args, threads);
    let label = format!("{} on {threads:?} threads", case.name);
    assert_eq!(
        output.status.code(),
        Some(case.code),
        "{label}: {}",
        text(&output.stderr)
    );
    if let Some(stdout) = expected(&case.name, "stdout") {
        assert_eq!(text(&output.stdout), text(&stdout), "{label}: stdout");
        assert_eq!(output.stdout, stdout, "{label}: stdout bytes");
    }
    match case.stderr.as_str() {
        "exact" => {
            let stderr = expected(&case.name, "stderr").expect("an expected stderr");
            assert_eq!(text(&output.stderr), text(&stderr), "{label}: stderr");
            assert_eq!(output.stderr, stderr, "{label}: stderr bytes");
        }
        "usage" => assert!(
            output.stderr.starts_with(b"usage: "),
            "{label}: {}",
            text(&output.stderr)
        ),
        "any" => {}
        other => panic!("{label}: stderr rule {other}"),
    }
}

#[test]
fn every_case_writes_the_committed_bytes_on_one_thread() {
    let cases = cases();
    assert_eq!(cases.len(), 81);
    for case in &cases {
        hold(case, Some("1"));
    }
}

#[test]
fn every_case_writes_the_committed_bytes_on_three_threads() {
    for case in &cases() {
        hold(case, Some("3"));
    }
}

#[test]
fn every_case_writes_the_committed_bytes_on_the_default_threads() {
    for case in &cases() {
        hold(case, None);
    }
}

#[test]
fn every_expected_file_belongs_to_a_case() {
    let cases = cases();
    let mut files = 0;
    for entry in fs::read_dir(fixtures().join("expected")).expect("expected/") {
        let name = entry.expect("an entry").file_name();
        let name = name.to_string_lossy();
        let (case, stream) = name.rsplit_once('.').expect("name.stream");
        assert!(matches!(stream, "stdout" | "stderr"), "{name}");
        let known = cases.iter().find(|c| c.name == case);
        let known = known.unwrap_or_else(|| panic!("{name}: no such case"));
        assert!(stream == "stdout" || known.stderr == "exact", "{name}");
        files += 1;
    }
    let exact = cases.iter().filter(|c| c.stderr == "exact").count();
    let with_stdout = cases
        .iter()
        .filter(|c| expected(&c.name, "stdout").is_some())
        .count();
    assert_eq!(files, exact + with_stdout);
    assert_eq!((exact, with_stdout), (71, 78));
}

#[test]
fn a_thread_bound_that_is_no_count_is_a_usage_error() {
    for bound in ["0", "-1", "two", "", "1.5"] {
        let output = run(&["cells.parquet"], Some(bound));
        assert_eq!(output.status.code(), Some(2), "{bound:?}");
        assert!(output.stdout.is_empty(), "{bound:?}");
        let stderr = text(&output.stderr);
        assert!(stderr.starts_with("usage: "), "{stderr}");
        assert!(
            stderr.contains("VEP_PARQUET_TO_TAB_THREADS") && stderr.contains("1 or more"),
            "{stderr}"
        );
    }
}

#[test]
fn a_type_outside_the_cell_contract_is_refused_by_file_and_column() {
    use arrow_array::{ArrayRef, BinaryArray, RecordBatch, StringArray};
    use parquet::arrow::ArrowWriter;
    use std::sync::Arc;

    let dir = std::env::temp_dir().join(format!("vep-parquet-to-tab-types-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("a scratch directory");
    let path = dir.join("bytes.parquet");
    let batch = RecordBatch::try_from_iter([
        ("chrom", Arc::new(StringArray::from(vec!["21"])) as ArrayRef),
        (
            "Feature",
            Arc::new(StringArray::from(vec!["f"])) as ArrayRef,
        ),
        (
            "RAW",
            Arc::new(BinaryArray::from(vec![&b"\x00\x01"[..]])) as ArrayRef,
        ),
    ])
    .expect("a batch");
    let mut writer = ArrowWriter::try_new(
        fs::File::create(&path).expect("a file"),
        batch.schema(),
        None,
    )
    .expect("a writer");
    writer.write(&batch).expect("a write");
    writer.close().expect("a footer");
    for mode in [&[][..], &["--canonical"][..]] {
        let mut args = mode.to_vec();
        args.push(path.to_str().expect("a UTF-8 path"));
        let output = run(&args, Some("2"));
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
        assert_eq!(
            text(&output.stderr),
            format!(
                "ERROR: [parquet_to_vep_tab] {}: column RAW is stored as Binary, which the cell \
                 contract has no tab text for\n",
                path.display()
            )
        );
    }
    fs::remove_dir_all(&dir).expect("the scratch directory goes");
}
