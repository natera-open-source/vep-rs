// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A failed run of either mode ends with one line on standard error,
//! `ERROR: [cache-converter] <message>`, and exits 1. `<message>` is the
//! failure and its causes, outermost first, joined by `: `. The line follows
//! whatever the run logged before it failed, and a run refused before its
//! first log line exits 1 on a standard error that takes no bytes.

use std::path::Path;
use std::process::{Command, Output, Stdio};

/// The converter, its log lines without colour codes.
fn converter() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep-cache-converter"));
    cmd.env("NO_COLOR", "1");
    cmd
}

/// Standard error of `output` as text.
fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

/// A cache conversion refused before its first log line, a `--json-dir` that
/// does not exist, prints the error line and nothing else.
#[test]
fn a_conversion_refused_before_it_logs_prints_the_error_line_alone() {
    let scratch = tempfile::tempdir().unwrap();
    let missing = scratch.path().join("no_such_cache");
    let output = converter()
        .arg("--json-dir")
        .arg(&missing)
        .arg("--output-dir")
        .arg(scratch.path().join("out"))
        .output()
        .expect("vep-cache-converter runs");
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr(&output),
        format!(
            "ERROR: [cache-converter] JSON directory does not exist: {}\n",
            missing.display()
        )
    );
}

/// A failure with a cause, a flat cache whose `info.json` is an empty object,
/// ends with one line that carries the failure and its cause; every line
/// before it is a log line.
#[test]
fn a_failure_with_a_cause_ends_on_one_line_that_carries_both() {
    let scratch = tempfile::tempdir().unwrap();
    let flat = scratch.path().join("flat");
    std::fs::create_dir_all(flat.join("transcripts")).unwrap();
    std::fs::write(flat.join("info.json"), "{}").unwrap();
    std::fs::write(flat.join("transcripts/21_1-1000000.json"), "[]\n").unwrap();

    let output = converter()
        .arg("--json-dir")
        .arg(&flat)
        .arg("--output-dir")
        .arg(scratch.path().join("out"))
        .output()
        .expect("vep-cache-converter runs");
    let err = stderr(&output);
    assert_eq!(output.status.code(), Some(1), "{err}");
    let (logged, last) = err
        .strip_suffix('\n')
        .and_then(|text| text.rsplit_once('\n'))
        .expect("log lines precede the error line");
    assert_eq!(
        last,
        format!(
            "ERROR: [cache-converter] {} is not a cache's info: \
             info.json lacks a non-empty species string",
            flat.join("info.json").display()
        ),
        "{err}"
    );
    assert!(logged.contains("Total transcripts: 1 region, "), "{err}");
    for line in logged.lines() {
        assert!(
            line.contains(" vep_cache_converter: "),
            "not a log line: {line:?}"
        );
    }
}

/// `convert-plugin-data` ends a failed run the same way: an input that does
/// not exist prints the error line alone.
#[test]
fn a_failed_plugin_data_conversion_prints_the_error_line_alone() {
    let scratch = tempfile::tempdir().unwrap();
    let missing = scratch.path().join("no_such_table.tsv.gz");
    let output = converter()
        .arg("convert-plugin-data")
        .arg("--input")
        .arg(&missing)
        .arg("--output")
        .arg(scratch.path().join("store.vpd"))
        .output()
        .expect("vep-cache-converter runs");
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr(&output),
        format!(
            "ERROR: [cache-converter] Input file not found: {}\n",
            missing.display()
        )
    );
}

/// The device that refuses every write, where the platform has one.
const FULL_DEVICE: &str = "/dev/full";

/// On a standard error that takes no bytes, the device that refuses every
/// write or a pipe whose reader has gone, a run refused before its first log
/// line loses the error line and exits 1.
#[test]
fn a_failed_run_exits_1_when_standard_error_takes_no_bytes() {
    let scratch = tempfile::tempdir().unwrap();
    let missing = scratch.path().join("no_such_cache");
    let mut streams: Vec<(&str, Stdio)> = Vec::new();
    if Path::new(FULL_DEVICE).exists() {
        let device = std::fs::File::options()
            .write(true)
            .open(FULL_DEVICE)
            .unwrap();
        streams.push((FULL_DEVICE, device.into()));
    }
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    streams.push(("a pipe whose reader has gone", writer.into()));

    for (name, stream) in streams {
        let status = converter()
            .arg("--json-dir")
            .arg(&missing)
            .arg("--output-dir")
            .arg(scratch.path().join("out"))
            .stdout(Stdio::null())
            .stderr(stream)
            .status()
            .expect("vep-cache-converter runs");
        assert_eq!(status.code(), Some(1), "{name}: {status}");
    }
}
