// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over the fixtures and compares the report, standard output and the exit
//! code with the expected bytes under `tests/fixtures/`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-intermediates-compare");
const ROOT: &str = env!("CARGO_MANIFEST_DIR");
const OUTPUT: &str = "tests/fixtures/output.txt";

/// A directory of the test's own, with a `tmp` inside it that the binary is given as
/// `TMPDIR`.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-intermediates-compare-golden-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(path.join("tmp")).unwrap();
        Scratch(path)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_string()
    }

    fn write(&self, name: &str, bytes: &[u8]) -> String {
        fs::write(self.0.join(name), bytes).unwrap();
        self.path(name)
    }

    /// Runs the binary from the package root, so fixture paths are relative to it.
    fn run(&self, args: &[&str]) -> Output {
        let output = Command::new(BIN)
            .current_dir(ROOT)
            .env("TMPDIR", self.0.join("tmp"))
            .args(args)
            .output()
            .unwrap();
        assert_eq!(
            fs::read_dir(self.0.join("tmp")).unwrap().count(),
            0,
            "the sorted intermediates' directory is removed"
        );
        output
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn expected(name: &str) -> Vec<u8> {
    fs::read(
        Path::new(ROOT)
            .join("tests/fixtures")
            .join(format!("expected_{name}.txt")),
    )
    .unwrap()
}

fn trace(name: &str) -> String {
    format!("tests/fixtures/{name}.tsv")
}

fn show(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// A complete report: exit code 0, nothing on standard error.
fn assert_report(output: &Output, name: &str) {
    assert_eq!(output.status.code(), Some(0), "{}", show(&output.stderr));
    assert!(output.stderr.is_empty(), "{}", show(&output.stderr));
    assert_eq!(show(&output.stdout), show(&expected(name)));
}

/// A report cut short or refused: exit code 1, the expected bytes on standard output and
/// one line that names the tool on standard error.
fn assert_failure(output: &Output, name: &str) {
    assert_eq!(output.status.code(), Some(1), "{}", show(&output.stderr));
    assert_eq!(show(&output.stdout), show(&expected(name)));
    assert_one_error_line(output);
}

fn assert_one_error_line(output: &Output) {
    let stderr = show(&output.stderr);
    assert!(
        stderr.starts_with("ERROR: [vep-intermediates-compare] ")
            && stderr.ends_with('\n')
            && stderr.matches('\n').count() == 1,
        "{stderr}"
    );
}

#[test]
fn fixture_pair_gives_the_expected_report() {
    let scratch = Scratch::new("report");
    let output = scratch.run(&["--perl-trace", &trace("trace"), "--rust-output", OUTPUT]);
    assert_report(&output, "report");
}

#[test]
fn discordant_only_drops_the_matching_rows() {
    let scratch = Scratch::new("discordant-only");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("trace"),
        "--rust-output",
        OUTPUT,
        "--discordant-only",
    ]);
    assert_report(&output, "discordant_only");
}

#[test]
fn output_names_the_report_file_and_a_dash_is_standard_output() {
    let scratch = Scratch::new("output");
    let report = scratch.path("report.txt");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("trace"),
        "--rust-output",
        OUTPUT,
        "--output",
        &report,
    ]);
    assert_eq!(output.status.code(), Some(0), "{}", show(&output.stderr));
    assert!(output.stdout.is_empty() && output.stderr.is_empty());
    assert_eq!(show(&fs::read(&report).unwrap()), show(&expected("report")));
    let output = scratch.run(&[
        "--perl-trace",
        &trace("trace"),
        "--rust-output",
        OUTPUT,
        "--output",
        "-",
    ]);
    assert_report(&output, "report");
}

#[test]
fn options_are_read_by_prefix_and_after_an_equals_sign() {
    let scratch = Scratch::new("prefix");
    let output = scratch.run(&[
        &format!("--perl={}", trace("trace")),
        "--rust",
        OUTPUT,
        "--disc",
    ]);
    assert_report(&output, "discordant_only");
}

#[test]
fn input_paths_are_normalised_as_pathlib_does() {
    let scratch = Scratch::new("normalised");
    let output = scratch.run(&[
        "--perl-trace",
        "tests/fixtures/trace.tsv/",
        "--rust-output",
        "./tests//fixtures/output.txt/.",
    ]);
    assert_report(&output, "report");
}

#[test]
fn trace_columns_are_found_by_header_name() {
    let scratch = Scratch::new("reordered");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("reordered_trace"),
        "--rust-output",
        OUTPUT,
        "--discordant-only",
    ]);
    assert_report(&output, "reordered");
}

#[test]
fn values_a_short_trace_row_lacks_print_as_none() {
    let scratch = Scratch::new("none-values");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("none_values_trace"),
        "--rust-output",
        OUTPUT,
        "--discordant-only",
    ]);
    assert_report(&output, "none_values");
}

#[test]
fn short_trace_rows_under_keys_of_their_own_are_counted() {
    let scratch = Scratch::new("short-unmatched");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("short_unmatched_trace"),
        "--rust-output",
        OUTPUT,
        "--discordant-only",
    ]);
    assert_report(&output, "short_unmatched");
}

#[test]
fn a_trace_of_a_header_alone_and_two_empty_inputs_give_reports() {
    let scratch = Scratch::new("empty");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("header_only_trace"),
        "--rust-output",
        OUTPUT,
    ]);
    assert_report(&output, "header_only");
    let empty_trace = scratch.write("trace.tsv", b"");
    let empty_output = scratch.write("output.txt", b"");
    let output = scratch.run(&["--perl-trace", &empty_trace, "--rust-output", &empty_output]);
    assert_report(&output, "empty");
}

#[test]
fn a_compared_row_without_its_consequence_ends_the_report_before_it() {
    let scratch = Scratch::new("short-row");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("short_row_trace"),
        "--rust-output",
        OUTPUT,
    ]);
    assert_failure(&output, "short_row");
    // The report file holds the same cut-short bytes.
    let report = scratch.path("report.txt");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("short_row_trace"),
        "--rust-output",
        OUTPUT,
        "--output",
        &report,
    ]);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert_eq!(
        show(&fs::read(&report).unwrap()),
        show(&expected("short_row"))
    );
}

#[test]
fn a_reference_codon_without_an_alternative_ends_the_report_before_its_row() {
    let scratch = Scratch::new("no-alt-codon");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("no_alt_codon_trace"),
        "--rust-output",
        OUTPUT,
    ]);
    assert_failure(&output, "no_alt_codon");
}

#[test]
fn keys_without_an_order_end_the_report_after_its_counts() {
    for name in ["unorderable", "absent_location"] {
        let scratch = Scratch::new(name);
        let output = scratch.run(&[
            "--perl-trace",
            &trace(&format!("{name}_trace")),
            "--rust-output",
            OUTPUT,
        ]);
        assert_failure(&output, name);
    }
}

#[test]
fn a_refused_input_writes_nothing_and_leaves_the_report_file_alone() {
    let scratch = Scratch::new("refused");
    let not_utf8 = scratch.write("comment.txt", b"## comment \xff\n");
    let gzip = scratch.write("output.txt.gz", b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03");
    let over_limit = scratch.write(
        "long_field.tsv",
        format!(
            "location\tallele\ttranscript_id\n{}\tA\tENST001\n",
            "x".repeat(131_073)
        )
        .as_bytes(),
    );
    let missing = scratch.path("missing.tsv");
    for (perl_trace, rust_output) in [
        (trace("no_location_trace"), OUTPUT.to_string()),
        (trace("trace"), not_utf8),
        (trace("trace"), gzip),
        (over_limit, OUTPUT.to_string()),
        (missing.clone(), OUTPUT.to_string()),
        (trace("trace"), missing),
        (trace("trace"), "tests/fixtures".to_string()),
    ] {
        let report = scratch.path("report.txt");
        fs::write(&report, b"kept").unwrap();
        let output = scratch.run(&[
            "--perl-trace",
            &perl_trace,
            "--rust-output",
            &rust_output,
            "--output",
            &report,
        ]);
        assert_eq!(output.status.code(), Some(1), "{perl_trace} {rust_output}");
        assert!(output.stdout.is_empty());
        assert_one_error_line(&output);
        assert_eq!(fs::read(&report).unwrap(), b"kept");
    }
}

#[test]
fn a_field_of_the_trace_holds_the_field_limit_and_no_more() {
    let scratch = Scratch::new("field-limit");
    let at_limit = scratch.write(
        "at_limit.tsv",
        format!(
            "location\tallele\ttranscript_id\n{}\tA\tENST001\n",
            "\u{e9}".repeat(131_072)
        )
        .as_bytes(),
    );
    let output = scratch.run(&["--perl-trace", &at_limit, "--rust-output", OUTPUT]);
    assert_eq!(output.status.code(), Some(0), "{}", show(&output.stderr));
    assert!(show(&output.stdout).starts_with("Perl rows: 1\nRust rows: 60\nUnique keys: 61\n\n"));
}

#[test]
fn an_unwritable_report_path_is_a_failure_with_nothing_written() {
    let scratch = Scratch::new("unwritable");
    let report = scratch.path("no/such/directory/report.txt");
    let output = scratch.run(&[
        "--perl-trace",
        &trace("trace"),
        "--rust-output",
        OUTPUT,
        "--output",
        &report,
    ]);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert_one_error_line(&output);
}

#[test]
fn the_inputs_are_read_whole_before_the_report_file_is_opened() {
    let scratch = Scratch::new("in-place");
    let in_place = scratch.write(
        "output.txt",
        &fs::read(Path::new(ROOT).join(OUTPUT)).unwrap(),
    );
    let output = scratch.run(&[
        "--perl-trace",
        &trace("trace"),
        "--rust-output",
        &in_place,
        "--output",
        &in_place,
    ]);
    assert_eq!(output.status.code(), Some(0), "{}", show(&output.stderr));
    assert_eq!(
        show(&fs::read(&in_place).unwrap()),
        show(&expected("report"))
    );
}

#[test]
fn usage_errors_exit_2_and_help_exits_0() {
    let scratch = Scratch::new("usage");
    for args in [
        vec![],
        vec!["--perl-trace", "tests/fixtures/trace.tsv"],
        vec![
            "--perl-trace",
            "tests/fixtures/trace.tsv",
            "--rust-output",
            OUTPUT,
            "stray",
        ],
        vec![
            "--perl-trace",
            "tests/fixtures/trace.tsv",
            "--rust-output",
            OUTPUT,
            "--bogus",
        ],
        vec![
            "--perl-trace",
            "tests/fixtures/trace.tsv",
            "--rust-output",
            OUTPUT,
            "--output",
        ],
        vec![
            "--perl-trace",
            "tests/fixtures/trace.tsv",
            "--rust-output",
            OUTPUT,
            "--discordant-only=1",
        ],
        vec![
            "--perl-trace",
            "tests/fixtures/trace.tsv",
            "--rust-output",
            OUTPUT,
            "--=x",
        ],
        vec!["-=x", "-h"],
    ] {
        let output = scratch.run(&args);
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(output.stdout.is_empty(), "{args:?}");
        let stderr = show(&output.stderr);
        assert!(
            stderr.starts_with("usage: vep-intermediates-compare ")
                && stderr.contains("\nvep-intermediates-compare: error: "),
            "{stderr}"
        );
    }
    for args in [vec!["-h"], vec!["--help"], vec!["--bogus", "--he"]] {
        let output = scratch.run(&args);
        assert_eq!(output.status.code(), Some(0), "{args:?}");
        assert!(show(&output.stdout).starts_with("usage: vep-intermediates-compare "));
        assert!(output.stderr.is_empty());
    }
}
