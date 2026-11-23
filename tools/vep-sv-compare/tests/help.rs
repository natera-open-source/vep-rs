// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary's help and usage errors: the help on stdout with exit code 0 for either
//! spelling of its option, and the usage block before the error line of a usage error
//! on stderr with exit code 2. Neither reads an input or writes a file.

use std::process::Command;

use vep_sv_compare::cli::{usage, HELP};
use vep_sv_compare::TempDir;

/// The exit code, stdout and stderr of the binary run on `arguments` in an empty
/// directory, which must be empty afterwards.
fn run(arguments: &[&str]) -> (Option<i32>, String, String) {
    let directory = TempDir::new("vep_sv_compare_help_").expect("a scratch directory");
    let output = Command::new(env!("CARGO_BIN_EXE_vep-sv-compare"))
        .args(arguments)
        .current_dir(&directory.path)
        .env("COLUMNS", "40")
        .output()
        .expect("the binary runs");
    let entries = std::fs::read_dir(&directory.path)
        .expect("the directory is read")
        .count();
    assert_eq!(entries, 0, "{arguments:?}");
    (
        output.status.code(),
        String::from_utf8(output.stdout).expect("stdout is UTF-8"),
        String::from_utf8(output.stderr).expect("stderr is UTF-8"),
    )
}

#[test]
fn the_help_is_printed_on_stdout_with_exit_code_0() {
    for arguments in [
        &["--help"][..],
        &["-h"],
        &["--he"],
        &["-hh"],
        &[
            "unrecognised",
            "--assembly",
            "grch38",
            "-h",
            "--assembly",
            "bad",
        ],
    ] {
        let (code, stdout, stderr) = run(arguments);
        assert_eq!(code, Some(0), "{arguments:?}");
        assert_eq!(stdout, HELP, "{arguments:?}");
        assert_eq!(stderr, "", "{arguments:?}");
    }
}

#[test]
fn a_usage_error_prints_the_usage_block_and_the_error_on_stderr_with_exit_code_2() {
    let cases: [(&[&str], &str); 5] = [
        (
            &["--assembly", "bad", "--help"],
            "argument --assembly: invalid choice: 'bad' (choose from 'grch37', 'grch38')",
        ),
        (
            &["--help=yes"],
            "argument -h/--help: ignored explicit argument 'yes'",
        ),
        (&["--", "-h"], "unrecognized arguments: -- -h"),
        (
            &["--perl-dir", "-h"],
            "argument --perl-dir: expected one argument",
        ),
        (&["extra"], "unrecognized arguments: extra"),
    ];
    for (arguments, message) in cases {
        let (code, stdout, stderr) = run(arguments);
        assert_eq!(code, Some(2), "{arguments:?}");
        assert_eq!(stdout, "", "{arguments:?}");
        assert_eq!(
            stderr,
            format!("{}vep-sv-compare: error: {message}\n", usage()),
            "{arguments:?}"
        );
    }
    assert!(usage().starts_with("usage: vep-sv-compare [-h] "));
    assert!(HELP.starts_with(usage()));
}
