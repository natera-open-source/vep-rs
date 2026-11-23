// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary's exit code and bytes for what no golden case holds: help, a count of
//! operands other than ten, paths as `pathlib.Path` reads them, a verdict that cannot be
//! written, a verdict that is written over, a name outside UTF-8, and the two bounds of
//! the tool: an integer of more than 128 bits and a threshold in digits outside ASCII.

mod common;

use common::{error_line, golden, report, run_in, Scratch, DEFAULTS};
use std::ffi::{OsStr, OsString};
use std::fs;
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

/// The ten operands of a run over `baseline` and `builder` with the reports of the
/// `equal_pass` case, the verdict operand as given.
fn operands(baseline: &Path, builder: &Path, verdict: &OsStr) -> Vec<OsString> {
    let mut args: Vec<OsString> = vec![
        baseline.into(),
        builder.into(),
        report("cache_diff_equal.json").into(),
        report("concordance_one.json").into(),
        verdict.into(),
    ];
    args.extend(DEFAULTS.iter().map(OsString::from));
    args
}

fn equal_metrics(verdict: &str) -> Vec<OsString> {
    let metrics = report("metrics_baseline.json");
    operands(&metrics, &metrics, OsStr::new(verdict))
}

/// The lines of the `equal_pass` case up to the verdict's own.
fn five_gate_lines() -> Vec<u8> {
    let whole = golden("equal_pass.stdout");
    let end = whole
        .windows(2)
        .position(|pair| pair == b"\n\n")
        .expect("an empty line before the verdict");
    whole[..=end].to_vec()
}

#[test]
fn help_prints_the_operands_and_leaves_with_zero() {
    let scratch = Scratch::new("help");
    for flag in ["-h", "--help"] {
        let output = run_in(&scratch.0, &[flag]);
        assert_eq!(output.status.code(), Some(0));
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stdout.starts_with("usage: vep-cache-output-gates BASELINE_METRICS "));
        assert!(stdout.contains("\nexit code: 0 for a PASS verdict, 1 for a FAIL verdict"));
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn any_count_of_operands_but_ten_is_refused_with_the_usage() {
    let scratch = Scratch::new("count");
    let ten = equal_metrics("verdict.json");
    for count in [0, 1, 9, 11] {
        let mut args = ten.clone();
        args.resize(count, OsString::from("0.01"));
        let output = run_in(&scratch.0, &args);
        assert_eq!(output.status.code(), Some(2), "{count} operands");
        assert!(output.stdout.is_empty());
        let stderr = String::from_utf8(output.stderr).unwrap();
        let expected = format!(
            "ERROR: [vep-cache-output-gates] ten operands are needed, {count} given\nusage: vep-cache-output-gates "
        );
        assert!(stderr.starts_with(&expected), "{stderr}");
        assert!(!scratch.0.join("verdict.json").exists());
    }
    // A help flag among ten operands is an operand: here a baseline that names nothing.
    let mut args = ten.clone();
    args[0] = "--help".into();
    let output = run_in(&scratch.0, &args);
    assert_eq!(output.status.code(), Some(2));
    assert!(error_line(&output).contains("cannot read the baseline metrics --help: "));
}

#[test]
fn a_verdict_that_cannot_be_written_leaves_the_five_gate_lines() {
    let scratch = Scratch::new("unwritable");
    let output = run_in(&scratch.0, &equal_metrics("no_such_directory/verdict.json"));
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(output.stdout, five_gate_lines());
    assert!(error_line(&output).starts_with(
        "ERROR: [vep-cache-output-gates] cannot write the verdict no_such_directory/verdict.json: "
    ));
}

#[test]
fn a_verdict_is_written_over_what_its_path_held() {
    let scratch = Scratch::new("over");
    fs::write(scratch.0.join("verdict.json"), "x".repeat(100_000)).unwrap();
    let output = run_in(&scratch.0, &equal_metrics("verdict.json"));
    assert_eq!(output.status.code(), Some(0));
    assert_eq!(output.stdout, golden("equal_pass.stdout"));
    assert_eq!(
        fs::read(scratch.0.join("verdict.json")).unwrap(),
        golden("equal_pass.verdict.json")
    );
}

#[test]
fn paths_are_read_as_pathlib_writes_them() {
    let scratch = Scratch::new("paths");
    fs::create_dir(scratch.0.join("out")).unwrap();
    fs::copy(
        report("metrics_baseline.json"),
        scratch.0.join("metrics.json"),
    )
    .unwrap();
    // The verdict is written at, and printed as, the path without its repeated and
    // trailing separators and its `.` components; a report's trailing separator is
    // dropped the same way before it is opened.
    let args = operands(
        Path::new(".//metrics.json/"),
        Path::new("metrics.json/."),
        OsStr::new("./out//./verdict.json/"),
    );
    let output = run_in(&scratch.0, &args);
    assert_eq!(output.status.code(), Some(0));
    let expected = String::from_utf8(golden("equal_pass.stdout"))
        .unwrap()
        .replace("  Report:  verdict.json\n", "  Report:  out/verdict.json\n");
    assert_eq!(String::from_utf8(output.stdout).unwrap(), expected);
    assert_eq!(
        fs::read(scratch.0.join("out/verdict.json")).unwrap(),
        golden("equal_pass.verdict.json")
    );
    assert!(output.stderr.is_empty());

    // The concordance summary is looked for at the path `pathlib` writes and opened at
    // the operand as it stands: a file named with a trailing separator is found, cannot be
    // opened and ends the run before any line.
    fs::copy(
        report("concordance_one.json"),
        scratch.0.join("concordance.json"),
    )
    .unwrap();
    let mut args = equal_metrics("verdict.json");
    args[3] = "concordance.json/".into();
    let output = run_in(&scratch.0, &args);
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert!(error_line(&output).contains("cannot read the concordance report concordance.json/: "));
    assert!(!scratch.0.join("verdict.json").exists());

    // A directory in a report's place cannot be read.
    let mut args = equal_metrics("verdict.json");
    args[2] = "out".into();
    let output = run_in(&scratch.0, &args);
    assert_eq!(output.status.code(), Some(2));
    assert!(error_line(&output).contains("cannot read the cache comparison report out: "));
}

#[test]
fn a_verdict_path_outside_utf8_is_printed_in_its_own_bytes() {
    let scratch = Scratch::new("bytes");
    let name = OsStr::from_bytes(b"verdict\xff.json");
    let output = run_in(
        &scratch.0,
        &operands(
            &report("metrics_baseline.json"),
            &report("metrics_baseline.json"),
            name,
        ),
    );
    assert_eq!(output.status.code(), Some(0));
    let mut expected = five_gate_lines();
    expected.extend_from_slice(b"\n  Verdict: PASS\n  Report:  verdict\xff.json\n");
    assert_eq!(output.stdout, expected);
    assert_eq!(
        fs::read(scratch.0.join(name)).unwrap(),
        golden("equal_pass.verdict.json")
    );
}

#[test]
fn the_first_report_that_cannot_be_read_is_the_one_named() {
    let scratch = Scratch::new("order");
    let mut args = equal_metrics("verdict.json");
    args[0] = "no_baseline.json".into();
    args[1] = report("metrics_truncated.json").into();
    args[2] = "no_cache_diff.json".into();
    args[3] = report("concordance_truncated.json").into();
    let output = run_in(&scratch.0, &args);
    assert_eq!(output.status.code(), Some(2));
    assert!(error_line(&output).contains("cannot read the baseline metrics no_baseline.json: "));
    args[0] = report("metrics_baseline.json").into();
    let output = run_in(&scratch.0, &args);
    assert!(error_line(&output).contains("metrics_truncated.json is not JSON: "));
    args[1] = report("metrics_baseline.json").into();
    let output = run_in(&scratch.0, &args);
    assert!(error_line(&output)
        .contains("cannot read the cache comparison report no_cache_diff.json: "));
    args[2] = report("cache_diff_equal.json").into();
    let output = run_in(&scratch.0, &args);
    assert!(error_line(&output).contains("concordance_truncated.json is not JSON: "));
    assert!(output.stdout.is_empty());
    // A threshold is read before any report.
    args[9] = "one percent".into();
    args[0] = "no_baseline.json".into();
    let output = run_in(&scratch.0, &args);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        error_line(&output),
        "ERROR: [vep-cache-output-gates] CONSEQUENCE_DELTA is not a number: \"one percent\"\n"
    );
}

#[test]
fn an_integer_of_more_than_128_bits_is_echoed_and_refused_in_arithmetic() {
    let scratch = Scratch::new("long");
    let long = "340282366920938463463374607431768211456";
    let write = |name: &str, total: &str| {
        let path = scratch.0.join(name);
        fs::write(&path, format!("{{\"total_annotation_lines\": {total}}}\n")).unwrap();
        path
    };
    let small = write("small.json", "7");
    let zero = write("zero.json", "0");
    let large = write("large.json", long);
    // Beside a baseline of zero the count is not computed with: it scores 1.0 and stands
    // in the line and in the report as written.
    let output = run_in(
        &scratch.0,
        &operands(&zero, &large, OsStr::new("verdict.json")),
    );
    assert_eq!(output.status.code(), Some(1));
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains(&format!(
        "\n  Total annotations delta: 1.000000 (0 vs {long}) (gate < 0.001) FAIL\n"
    )));
    let verdict = fs::read_to_string(scratch.0.join("verdict.json")).unwrap();
    assert!(verdict.contains(&format!("\"builder\": {long},\n")));
    // Where the delta divides, such a count ends the run after the first gate's line.
    for (baseline, builder) in [(&small, &large), (&large, &small)] {
        let output = run_in(
            &scratch.0,
            &operands(baseline, builder, OsStr::new("refused.json")),
        );
        assert_eq!(output.status.code(), Some(2));
        assert_eq!(
            output.stdout,
            b"  Consequence F1: 1.000000 (gate >= 0.999) PASS\n"
        );
        assert!(error_line(&output).contains("is an integer of more than 128 bits"));
        assert!(!scratch.0.join("refused.json").exists());
    }
}

#[test]
fn a_threshold_in_digits_or_blanks_outside_ascii_is_refused() {
    let scratch = Scratch::new("digits");
    for word in ["\u{661}", "\u{a0}0.01", "0.01\u{2003}"] {
        let mut args = equal_metrics("verdict.json");
        args[6] = word.into();
        let output = run_in(&scratch.0, &args);
        assert_eq!(output.status.code(), Some(2), "{word:?}");
        assert!(output.stdout.is_empty());
        assert!(error_line(&output).contains("TOTAL_DELTA is not a number: "));
    }
}
