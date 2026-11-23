// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary's exit code and bytes for the ways its two files are named: the table beside
//! the input, a table written over its own input, a pipe, a device, files that cannot be
//! read or written, refused command lines, the thread bound and help.

mod common;

use common::{
    assert_clean_exit, assert_refused, assert_usage_error, bytes, fixture, run_bound, run_in, text,
    Scratch, BIN, THREADS_ENV,
};
use std::fs;
use std::io::Write;
use std::os::unix::fs::{symlink, PermissionsExt};
use std::path::Path;
use std::process::{Command, Stdio};

/// A scratch directory that holds the `tuple.tsv` fixture under `name`.
fn with_table(tag: &str, name: &str) -> Scratch {
    let scratch = Scratch::new(tag);
    let path = scratch.path(name);
    fs::create_dir_all(path.parent().unwrap()).unwrap();
    fs::copy(fixture("tuple.tsv"), path).unwrap();
    scratch
}

fn assert_tuple_goldens(table: &Path, stdout: &[u8]) {
    assert_eq!(fs::read(table).unwrap(), bytes("expected_tuple.tsv"));
    assert_eq!(stdout, bytes("expected_tuple.stdout"));
}

#[test]
fn the_table_is_written_beside_the_input_where_no_output_is_named() {
    let scratch = with_table("beside", "report/discordant.tsv");
    let elsewhere = Scratch::new("beside-elsewhere");
    let output = run_in(&elsewhere.0, &[scratch.path("report/discordant.tsv")]);
    assert_clean_exit(&output);
    assert_tuple_goldens(
        &scratch.path("report/discordant_classified.tsv"),
        &output.stdout,
    );
    assert_eq!(fs::read_dir(&elsewhere.0).unwrap().count(), 0);
}

#[test]
fn the_default_table_s_directory_is_the_input_s_as_pathlib_normalises_it() {
    let scratch = with_table("beside-relative", "report/discordant.tsv");
    fs::copy(fixture("tuple.tsv"), scratch.path("here.tsv")).unwrap();
    let output = run_in(&scratch.0, &["here.tsv"]);
    assert_clean_exit(&output);
    assert_tuple_goldens(&scratch.path("discordant_classified.tsv"), &output.stdout);

    fs::create_dir(scratch.path("other")).unwrap();
    let output = run_in(&scratch.0, &["./other//../report/./discordant.tsv"]);
    assert_clean_exit(&output);
    assert_tuple_goldens(
        &scratch.path("report/discordant_classified.tsv"),
        &output.stdout,
    );
}

#[test]
fn the_output_option_names_the_table_in_any_of_its_spellings() {
    let scratch = with_table("spellings", "in.tsv");
    for (index, words) in [
        &["in.tsv", "-o", "out.tsv"][..],
        &["-o", "out.tsv", "in.tsv"],
        &["--output=out.tsv", "in.tsv"],
        &["in.tsv", "--out", "out.tsv"],
        &["-oout.tsv", "--", "in.tsv"],
        &["-o", "first.tsv", "in.tsv", "-o=out.tsv"],
    ]
    .into_iter()
    .enumerate()
    {
        let _ = fs::remove_file(scratch.path("out.tsv"));
        let output = run_in(&scratch.0, words);
        assert_clean_exit(&output);
        assert_tuple_goldens(&scratch.path("out.tsv"), &output.stdout);
        assert!(
            !scratch.path("discordant_classified.tsv").exists(),
            "case {index}"
        );
        assert!(!scratch.path("first.tsv").exists(), "case {index}");
    }
}

#[test]
fn a_table_written_over_its_own_input_is_read_whole_first() {
    // The default table of an input that carries the default's own name.
    let scratch = with_table("in-place", "discordant_classified.tsv");
    for bound in [Some("1"), Some("4")] {
        fs::copy(
            fixture("tuple.tsv"),
            scratch.path("discordant_classified.tsv"),
        )
        .unwrap();
        let output = run_bound(&scratch.0, &["discordant_classified.tsv"], bound);
        assert_clean_exit(&output);
        assert_tuple_goldens(&scratch.path("discordant_classified.tsv"), &output.stdout);
    }
    // The input named as the output, under its own name, a symbolic link and a hard link.
    fs::copy(fixture("tuple.tsv"), scratch.path("in.tsv")).unwrap();
    symlink(scratch.path("in.tsv"), scratch.path("link.tsv")).unwrap();
    fs::hard_link(scratch.path("in.tsv"), scratch.path("hard.tsv")).unwrap();
    for name in ["in.tsv", "link.tsv", "hard.tsv"] {
        fs::write(scratch.path("in.tsv"), bytes("tuple.tsv")).unwrap();
        let output = run_in(&scratch.0, &["in.tsv", "-o", name]);
        assert_clean_exit(&output);
        assert_tuple_goldens(&scratch.path("in.tsv"), &output.stdout);
    }
}

#[test]
fn an_existing_table_is_replaced_whole() {
    let scratch = with_table("replaced", "in.tsv");
    fs::write(scratch.path("out.tsv"), "x".repeat(1 << 16)).unwrap();
    let output = run_in(&scratch.0, &["in.tsv", "-o", "out.tsv"]);
    assert_clean_exit(&output);
    assert_tuple_goldens(&scratch.path("out.tsv"), &output.stdout);
}

#[test]
fn a_refused_input_leaves_an_existing_table_as_it_was() {
    let scratch = Scratch::new("left-alone");
    for file in [
        "not_utf8.tsv",
        "no_location.tsv",
        "missing_columns.tsv",
        "empty.tsv",
    ] {
        fs::write(scratch.path("out.tsv"), "kept\n").unwrap();
        let output = run_in(
            &scratch.0,
            &[fixture(file).as_os_str(), "-o".as_ref(), "out.tsv".as_ref()],
        );
        assert_eq!(output.status.code(), Some(1), "{file}");
        assert!(output.stdout.is_empty(), "{file}");
        assert_eq!(
            fs::read(scratch.path("out.tsv")).unwrap(),
            b"kept\n",
            "{file}"
        );
    }
}

#[test]
fn a_pipe_is_read_as_a_file_is() {
    let scratch = Scratch::new("pipe");
    for bound in ["1", "3"] {
        let mut child = Command::new(BIN)
            .current_dir(&scratch.0)
            .args(["/dev/stdin", "-o", "out.tsv"])
            .env(THREADS_ENV, bound)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child
            .stdin
            .take()
            .unwrap()
            .write_all(&bytes("tuple.tsv"))
            .unwrap();
        let output = child.wait_with_output().unwrap();
        assert_clean_exit(&output);
        assert_tuple_goldens(&scratch.path("out.tsv"), &output.stdout);
    }
}

#[test]
fn a_table_sent_to_stdout_stands_before_the_summary() {
    if !Path::new("/dev/stdout").exists() {
        return;
    }
    let scratch = with_table("to-stdout", "in.tsv");
    let output = run_in(&scratch.0, &["in.tsv", "-o", "/dev/stdout"]);
    assert_clean_exit(&output);
    let expected = [bytes("expected_tuple.tsv"), bytes("expected_tuple.stdout")].concat();
    assert_eq!(output.stdout, expected);

    let output = run_in(&scratch.0, &["in.tsv", "-o", "/dev/null"]);
    assert_clean_exit(&output);
    assert_eq!(output.stdout, bytes("expected_tuple.stdout"));
}

#[test]
fn an_input_that_cannot_be_read_is_exit_code_1_and_writes_nothing() {
    let scratch = with_table("unreadable", "in.tsv");
    fs::create_dir(scratch.path("directory")).unwrap();
    let cases = [
        (
            "missing.tsv",
            "missing.tsv: No such file or directory (os error 2)",
        ),
        ("", ": No such file or directory (os error 2)"),
        ("-", "-: No such file or directory (os error 2)"),
        ("directory", "directory: is a directory"),
        ("in.tsv/", "in.tsv/: Not a directory (os error 20)"),
    ];
    for (input, detail) in cases {
        let output = run_in(&scratch.0, &[input, "-o", "out.tsv"]);
        assert_refused(&output, detail);
        assert!(output.stdout.is_empty());
        assert!(!scratch.path("out.tsv").exists(), "{input}");
    }
    let output = run_in(&scratch.0, &["directory"]);
    assert_refused(&output, "directory: is a directory");
    assert!(!scratch.path("discordant_classified.tsv").exists());

    // A file without the read permission, where the test's user is held to permissions.
    fs::set_permissions(scratch.path("in.tsv"), fs::Permissions::from_mode(0o000)).unwrap();
    if fs::File::open(scratch.path("in.tsv")).is_err() {
        let output = run_in(&scratch.0, &["in.tsv", "-o", "out.tsv"]);
        assert_refused(&output, "in.tsv: Permission denied (os error 13)");
        assert!(!scratch.path("out.tsv").exists());
    }
}

#[test]
fn a_table_that_cannot_be_written_is_exit_code_1_with_no_summary() {
    let scratch = with_table("unwritable", "in.tsv");
    fs::create_dir(scratch.path("directory")).unwrap();
    let cases = [
        (
            "missing/out.tsv",
            "missing/out.tsv: No such file or directory (os error 2)",
        ),
        ("", ": No such file or directory (os error 2)"),
        ("directory", "directory: Is a directory (os error 21)"),
        (
            "in.tsv/out.tsv",
            "in.tsv/out.tsv: Not a directory (os error 20)",
        ),
    ];
    for (table, detail) in cases {
        let output = run_in(&scratch.0, &["in.tsv", "-o", table]);
        assert_refused(&output, detail);
        assert!(output.stdout.is_empty(), "{table}");
    }
    assert_eq!(
        fs::read(scratch.path("in.tsv")).unwrap(),
        bytes("tuple.tsv")
    );
    if Path::new("/dev/full").exists() {
        for bound in [Some("1"), Some("2")] {
            let output = run_bound(&scratch.0, &["in.tsv", "-o", "/dev/full"], bound);
            assert_refused(&output, "/dev/full: No space left on device (os error 28)");
            assert!(output.stdout.is_empty());
        }
    }
}

#[test]
fn a_summary_that_stdout_does_not_take_is_exit_code_120_behind_a_whole_table() {
    if !Path::new("/dev/full").exists() {
        return;
    }
    let scratch = with_table("stdout-full", "in.tsv");
    let output = Command::new(BIN)
        .current_dir(&scratch.0)
        .args(["in.tsv", "-o", "out.tsv"])
        .stdin(Stdio::null())
        .stdout(
            fs::OpenOptions::new()
                .write(true)
                .open("/dev/full")
                .unwrap(),
        )
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(120));
    assert_eq!(
        text(&output.stderr),
        "ERROR: [vep-discordant-classes] stdout: No space left on device (os error 28)\n"
    );
    assert_eq!(
        fs::read(scratch.path("out.tsv")).unwrap(),
        bytes("expected_tuple.tsv")
    );
}

#[test]
fn a_refused_command_line_is_exit_code_2_with_the_usage_line() {
    let scratch = with_table("usage", "in.tsv");
    let cases: [(&[&str], &str); 12] = [
        (&[], "the following arguments are required: input"),
        (
            &["-o", "out.tsv"],
            "the following arguments are required: input",
        ),
        (&["--"], "the following arguments are required: input"),
        (
            &["in.tsv", "-o"],
            "argument -o/--output: expected one argument",
        ),
        (
            &["in.tsv", "-o", "--", "out.tsv"],
            "argument -o/--output: expected one argument",
        ),
        (
            &["in.tsv", "-o", "-x"],
            "argument -o/--output: expected one argument",
        ),
        (
            &["in.tsv", "extra", "more"],
            "unrecognized arguments: extra more",
        ),
        (
            &["in.tsv", "--threads", "2"],
            "unrecognized arguments: --threads 2",
        ),
        (
            &["--", "in.tsv", "-o", "out.tsv"],
            "unrecognized arguments: -o out.tsv",
        ),
        (
            &["--=x", "in.tsv"],
            "ambiguous option: --=x could match --help, --output",
        ),
        (
            &["--help=1"],
            "argument -h/--help: ignored explicit argument '1'",
        ),
        (
            &["-h-", "in.tsv"],
            "argument -h/--help: ignored explicit argument '-'",
        ),
    ];
    for (words, message) in cases {
        let output = run_in(&scratch.0, words);
        assert_usage_error(&output, message);
        assert!(!scratch.path("out.tsv").exists(), "{words:?}");
        assert!(
            !scratch.path("discordant_classified.tsv").exists(),
            "{words:?}"
        );
    }
}

#[test]
fn help_prints_the_usage_and_the_options_and_reads_no_table() {
    let scratch = Scratch::new("help");
    for words in [
        &["-h"][..],
        &["--help"],
        &["--he"],
        &["missing.tsv", "-h"],
        &["-hx"],
    ] {
        let output = run_in(&scratch.0, words);
        assert_clean_exit(&output);
        assert_eq!(
            text(&output.stdout),
            "usage: vep-discordant-classes [-h] [-o OUTPUT] input\n\
             \n\
             Classify discordant tuples\n\
             \n\
             positional arguments:\n\
             \x20 input                Path to discordant.tsv\n\
             \n\
             options:\n\
             \x20 -h, --help           show this help message and exit\n\
             \x20 -o, --output OUTPUT  Output path (default: discordant_classified.tsv)\n"
        );
    }
}

#[test]
fn a_thread_bound_that_is_no_count_is_refused_before_the_input_is_read() {
    let scratch = with_table("bound", "in.tsv");
    for value in ["0", "many", "-2", "1.5"] {
        let output = run_bound(&scratch.0, &["in.tsv", "-o", "out.tsv"], Some(value));
        assert_refused(
            &output,
            &format!("{THREADS_ENV}={value}: expected a whole number of threads, 1 or more"),
        );
        assert!(!scratch.path("out.tsv").exists());
    }
    let output = run_bound(&scratch.0, &["in.tsv", "-o", "out.tsv"], Some(""));
    assert_clean_exit(&output);
    let output = run_bound(&scratch.0, &["missing.tsv"], Some("0"));
    assert_refused(
        &output,
        &format!("{THREADS_ENV}=0: expected a whole number of threads, 1 or more"),
    );
}
