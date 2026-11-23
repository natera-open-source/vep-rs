// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary with the spellings of its options, with refused command lines and with a
//! request for help, and holds the exit code and the bytes of each.

mod common;

use common::{
    assert_clean_exit, fixture, golden, run, run_args, run_with_path, stderr, stdout, Scratch,
};
use std::ffi::OsStr;
use std::os::unix::ffi::OsStrExt;

const USAGE: &str = "usage: vep-fields-compare [-h] --perl PERL --rust RUST --out OUT [--examples EXAMPLES] [--tmp-dir TMP_DIR] [--keep-only {regulatory_feature_rows}]\n";

#[test]
fn abbreviated_options_write_the_bytes_of_their_full_names() {
    let scratch = Scratch::new("abbreviated");
    let out = scratch.report();
    let output = run_args(&[
        "--per",
        "tests/fixtures/reference.txt",
        "--ru=tests/fixtures/scored.txt",
        "--o",
        out.to_str().unwrap(),
        "--tmp",
        scratch.0.to_str().unwrap(),
        "--ex",
        "1",
    ]);
    assert_clean_exit(&output);
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        golden("expected_fields_examples1.json")
    );
    assert_eq!(stdout(&output), golden("expected_stdout.txt"));

    let output = run_args(&[
        "--p",
        "tests/fixtures/regulatory_reference.txt",
        "--r",
        "tests/fixtures/regulatory_scored.txt",
        "--o",
        out.to_str().unwrap(),
        "--t",
        scratch.0.to_str().unwrap(),
        "--k",
        "regulatory_feature_rows",
    ]);
    assert_clean_exit(&output);
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        golden("expected_regulatory_fields.json")
    );
}

#[test]
fn examples_takes_a_count_with_underscores_and_one_past_any_number_of_shapes() {
    for count in ["1_000", "99999999999999999999999999", " 20 ", "+20"] {
        let scratch = Scratch::new("count");
        let (output, json) = run(
            &scratch,
            "tests/fixtures/reference.txt",
            "tests/fixtures/scored.txt",
            &["--examples", count],
        );
        assert_clean_exit(&output);
        assert_eq!(json, golden("expected_fields.json"), "--examples {count:?}");
    }
}

#[test]
fn missing_input_is_a_single_line_error() {
    let scratch = Scratch::new("missing-both");
    let missing = scratch.0.join("nope");
    let out = scratch.0.join("r.json");
    let output = run_args(&[
        OsStr::new("--perl"),
        missing.as_os_str(),
        OsStr::new("--rust"),
        missing.as_os_str(),
        OsStr::new("--out"),
        out.as_os_str(),
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        stderr(&output),
        format!(
            "ERROR: [compare_vep_fields] not a file: {}\n",
            missing.display()
        )
    );
    assert!(output.stdout.is_empty());
    assert!(!out.exists());
}

#[test]
fn a_missing_input_named_outside_utf8_is_named_by_its_surrogate_escapes() {
    let output = run_args(&[
        OsStr::new("--perl"),
        OsStr::from_bytes(b"no\xffpe\xe2\x82/x\xc3\xa9.txt"),
        OsStr::new("--rust"),
        OsStr::new("tests/fixtures/scored.txt"),
        OsStr::new("--out"),
        OsStr::new("o.json"),
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        output.stderr,
        b"ERROR: [compare_vep_fields] not a file: no\\udcffpe\\udce2\\udc82/x\xc3\xa9.txt\n"
    );
    assert!(output.stdout.is_empty());
}

/// A run that fails after its command line is read: exit code 1, one stderr line that
/// opens with the tool's name, no summary on stdout and no directory of sort intermediates
/// left under `--tmp-dir`.
fn assert_failure(scratch: &Scratch, output: &std::process::Output) {
    assert_eq!(output.status.code(), Some(1), "{}", stderr(output));
    let message = stderr(output);
    assert!(
        message.starts_with("ERROR: [vep-fields-compare] "),
        "{message}"
    );
    assert_eq!(message.matches('\n').count(), 1, "{message}");
    assert!(output.stdout.is_empty());
    let left: Vec<String> = std::fs::read_dir(&scratch.0)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name.starts_with("vep_fields_"))
        .collect();
    assert!(left.is_empty(), "{left:?}");
}

#[test]
fn an_unwritable_out_exits_1_and_prints_no_summary() {
    let scratch = Scratch::new("unwritable-out");
    let tmp = scratch.0.to_str().unwrap();
    let reference = "tests/fixtures/pair_reference.txt";
    let scored = "tests/fixtures/pair_scored.txt";
    // `--out` names a directory.
    let output = run_args(&[
        "--perl",
        reference,
        "--rust",
        scored,
        "--out",
        tmp,
        "--tmp-dir",
        tmp,
    ]);
    assert_failure(&scratch, &output);
    // `--out` lies under a regular file.
    let under_a_file = format!("{reference}/x/fields.json");
    let output = run_args(&[
        "--perl",
        reference,
        "--rust",
        scored,
        "--out",
        &under_a_file,
        "--tmp-dir",
        tmp,
    ]);
    assert_failure(&scratch, &output);
}

#[test]
fn a_tmp_dir_that_is_no_directory_exits_1() {
    let scratch = Scratch::new("no-tmp-dir");
    let out = scratch.report();
    let absent = scratch.0.join("absent");
    for tmp in [absent.to_str().unwrap(), "tests/fixtures/pair_scored.txt"] {
        let output = run_args(&[
            "--perl",
            "tests/fixtures/pair_reference.txt",
            "--rust",
            "tests/fixtures/pair_scored.txt",
            "--out",
            out.to_str().unwrap(),
            "--tmp-dir",
            tmp,
        ]);
        assert_failure(&scratch, &output);
        assert!(!out.exists());
    }
}

#[test]
fn an_unreadable_input_exits_1() {
    use std::os::unix::fs::PermissionsExt;
    let scratch = Scratch::new("unreadable");
    let input = scratch.0.join("reference.txt");
    std::fs::copy(fixture("pair_reference.txt"), &input).unwrap();
    std::fs::set_permissions(&input, std::fs::Permissions::from_mode(0o000)).unwrap();
    if std::fs::File::open(&input).is_ok() {
        // A user who may read any file cannot meet this failure.
        return;
    }
    let output = run_args(&[
        "--perl",
        input.to_str().unwrap(),
        "--rust",
        "tests/fixtures/pair_scored.txt",
        "--out",
        scratch.report().to_str().unwrap(),
        "--tmp-dir",
        scratch.0.to_str().unwrap(),
    ]);
    assert_failure(&scratch, &output);
    assert!(stderr(&output).ends_with("Permission denied (os error 13)\n"));
    assert!(!scratch.report().exists());
}

#[test]
fn a_host_without_sort_orders_the_rows_in_memory_to_the_same_bytes() {
    let scratch = Scratch::new("no-sort");
    let no_programs = scratch.0.join("bin");
    std::fs::create_dir(&no_programs).unwrap();
    let out = scratch.report();
    for (pair, expected) in [("", "expected"), ("odd_bytes_", "expected_odd_bytes")] {
        let output = run_with_path(
            &no_programs,
            &[
                "--perl",
                &format!("tests/fixtures/{pair}reference.txt"),
                "--rust",
                &format!("tests/fixtures/{pair}scored.txt"),
                "--out",
                out.to_str().unwrap(),
                "--tmp-dir",
                scratch.0.to_str().unwrap(),
            ],
        );
        assert_clean_exit(&output);
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            golden(&format!("{expected}_fields.json"))
        );
        assert_eq!(stdout(&output), golden(&format!("{expected}_stdout.txt")));
    }
}

#[test]
fn a_sort_that_fails_exits_1_with_its_exit_code_and_its_stderr() {
    use std::os::unix::fs::PermissionsExt;
    let scratch = Scratch::new("failing-sort");
    let programs = scratch.0.join("bin");
    std::fs::create_dir(&programs).unwrap();
    let sort = programs.join("sort");
    std::fs::write(&sort, "#!/bin/sh\necho boom >&2\nexit 3\n").unwrap();
    std::fs::set_permissions(&sort, std::fs::Permissions::from_mode(0o755)).unwrap();
    let output = run_with_path(
        &programs,
        &[
            "--perl",
            "tests/fixtures/pair_reference.txt",
            "--rust",
            "tests/fixtures/pair_scored.txt",
            "--out",
            scratch.report().to_str().unwrap(),
            "--tmp-dir",
            scratch.0.to_str().unwrap(),
        ],
    );
    assert_failure(&scratch, &output);
    assert_eq!(
        stderr(&output),
        "ERROR: [vep-fields-compare] sort failed (rc=3): boom\n"
    );
    assert!(!scratch.report().exists());
}

#[test]
fn an_input_that_is_a_directory_is_not_a_file() {
    let output = run_args(&[
        "--perl",
        "tests/fixtures",
        "--rust",
        "tests/fixtures/pair_scored.txt",
        "--out",
        "o.json",
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        stderr(&output),
        "ERROR: [compare_vep_fields] not a file: tests/fixtures\n"
    );
    assert!(output.stdout.is_empty());
}

#[test]
fn a_refused_command_line_exits_2_with_the_usage_line_and_one_error_line() {
    let reference = "tests/fixtures/reference.txt";
    let scored = "tests/fixtures/scored.txt";
    let cases: [(&[&str], &str); 9] = [
        (
            &[],
            "the following arguments are required: --perl, --rust, --out",
        ),
        (
            &["--perl", reference, "--rust", scored],
            "the following arguments are required: --out",
        ),
        (
            &[
                "--perl",
                reference,
                "--rust",
                scored,
                "--out",
                "o.json",
                "--bogus=1",
                "stray",
            ],
            "unrecognized arguments: --bogus=1 stray",
        ),
        (
            &[
                "--perl",
                reference,
                "--rust",
                scored,
                "--out",
                "o.json",
                "--examples",
                "1.5",
            ],
            "argument --examples: invalid int value: '1.5'",
        ),
        (
            &[
                "--perl",
                reference,
                "--rust",
                scored,
                "--out",
                "o.json",
                "--examples",
            ],
            "argument --examples: expected one argument",
        ),
        (
            &["--perl", "-x", "--rust", scored, "--out", "o.json"],
            "argument --perl: expected one argument",
        ),
        (
            &[
                "--perl", reference, "--rust", scored, "--out", "o.json", "--",
            ],
            "unrecognized arguments: --",
        ),
        (
            &[
                "--perl",
                reference,
                "--rust",
                scored,
                "--out",
                "o.json",
                "--keep-only",
                "all",
            ],
            "argument --keep-only: invalid choice: 'all' (choose from 'regulatory_feature_rows')",
        ),
        (
            &["--help=x"],
            "argument -h/--help: ignored explicit argument 'x'",
        ),
    ];
    for (args, message) in cases {
        let output = run_args(args);
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert_eq!(
            stderr(&output),
            format!("{USAGE}vep-fields-compare: error: {message}\n"),
            "{args:?}"
        );
        assert!(output.stdout.is_empty(), "{args:?}");
    }
    assert!(!std::path::Path::new(common::ROOT).join("o.json").exists());
}

#[test]
fn help_prints_the_usage_and_the_options_and_exits_0() {
    let asked: [&[&str]; 5] = [
        &["-h"],
        &["--help"],
        &["--he"],
        &["--bogus", "-h"],
        &["-h", "--examples", "x"],
    ];
    for args in asked {
        let output = run_args(args);
        assert_clean_exit(&output);
        let text = stdout(&output);
        assert!(text.starts_with(USAGE), "{args:?}: {text}");
        for line in [
            "\nField-level agreement between two VEP default-format outputs.\n",
            "\n  -h, --help            show this help message and exit\n",
            "\n  --examples EXAMPLES   mismatch shapes kept per field (default 20)\n",
            "\n  --keep-only {regulatory_feature_rows}\n",
        ] {
            assert!(text.contains(line), "{args:?}: {line:?} in {text}");
        }
    }
}
