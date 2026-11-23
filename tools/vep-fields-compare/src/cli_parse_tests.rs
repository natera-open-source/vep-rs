// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line against what `argparse` answers for the same arguments: the values an
//! accepted line gives, and the message after `error: ` of a refused one.

use super::*;

const BASE: [&str; 6] = ["--perl", "p", "--rust", "r", "--out", "o"];

fn parse_strs(args: &[&str]) -> Result<Command, String> {
    parse(&args.iter().map(OsString::from).collect::<Vec<_>>())
}

fn run_args(args: &[&str]) -> Args {
    match parse_strs(args) {
        Ok(Command::Run(args)) => args,
        other => panic!("{args:?}: {other:?}"),
    }
}

/// `BASE` then `extra`.
fn with(extra: &[&str]) -> Result<Command, String> {
    parse_strs(&[&BASE[..], extra].concat())
}

fn error(args: &[&str]) -> String {
    match parse_strs(args) {
        Err(message) => message,
        Ok(command) => panic!("{args:?}: {command:?}"),
    }
}

fn error_with(extra: &[&str]) -> String {
    error(&[&BASE[..], extra].concat())
}

fn is_help(args: &[&str]) -> bool {
    matches!(parse_strs(args), Ok(Command::Help))
}

#[test]
fn arguments_take_both_spellings_and_default_examples_to_20() {
    let args = run_args(&["--perl", "p", "--rust=r", "--out", "o"]);
    assert_eq!(args.perl, PathBuf::from("p"));
    assert_eq!(args.rust, PathBuf::from("r"));
    assert_eq!(args.out, PathBuf::from("o"));
    assert_eq!(args.examples, 20);
    assert!(args.tmp_dir.is_none());
    assert!(args.keep_only.is_none());
    let args = run_args(&[&BASE[..], &["--examples", "1", "--tmp-dir", "t"]].concat());
    assert_eq!(args.examples, 1);
    assert_eq!(args.tmp_dir, Some(PathBuf::from("t")));
}

#[test]
fn keep_only_takes_its_one_class_in_both_spellings_and_refuses_any_other_value() {
    for spelling in [
        &["--keep-only", "regulatory_feature_rows"][..],
        &["--keep-only=regulatory_feature_rows"][..],
        &[
            "--keep-only=regulatory_feature_rows",
            "--k=regulatory_feature_rows",
        ][..],
    ] {
        let Ok(Command::Run(args)) = with(spelling) else {
            panic!("{spelling:?}");
        };
        assert_eq!(args.keep_only, Some(KeepOnly::RegulatoryFeatureRows));
    }
    let refused = |value: &str| {
        format!(
            "argument --keep-only: invalid choice: '{value}' (choose from 'regulatory_feature_rows')"
        )
    };
    assert_eq!(error_with(&["--keep-only", "bogus"]), refused("bogus"));
    assert_eq!(error_with(&["--keep-only", ""]), refused(""));
    assert_eq!(
        error_with(&["--keep-only", "Regulatory_feature_rows"]),
        refused("Regulatory_feature_rows")
    );
    assert_eq!(
        error_with(&["--keep-only", " regulatory_feature_rows"]),
        refused(" regulatory_feature_rows")
    );
    assert_eq!(
        error_with(&[
            "--keep-only",
            "regulatory_feature_rows",
            "--keep-only",
            "bogus"
        ]),
        refused("bogus")
    );
    assert_eq!(
        error_with(&["--keep-only"]),
        "argument --keep-only: expected one argument"
    );
    assert_eq!(
        error(&["--keep-only", "--perl", "p", "--rust", "r", "--out", "o"]),
        "argument --keep-only: expected one argument"
    );
}

#[test]
fn path_values_are_normalised_on_their_bytes_in_both_spellings() {
    let args = run_args(&[
        "--perl",
        "./p//q/",
        "--rust=r/.",
        "--out",
        "o",
        "--tmp-dir=t/",
    ]);
    assert_eq!(args.perl, PathBuf::from("p/q"));
    assert_eq!(args.rust, PathBuf::from("r"));
    assert_eq!(args.tmp_dir, Some(PathBuf::from("t")));
    let raw = OsStr::from_bytes(b"./ref\xff.txt/");
    let mut inline = OsString::from("--perl=");
    inline.push(raw);
    for argv in [
        vec![
            inline,
            "--rust".into(),
            "r".into(),
            "--out".into(),
            "o".into(),
        ],
        vec![
            "--perl".into(),
            raw.to_os_string(),
            "--rust".into(),
            "r".into(),
            "--out".into(),
            "o".into(),
        ],
    ] {
        let Ok(Command::Run(args)) = parse(&argv) else {
            panic!("{argv:?}");
        };
        assert_eq!(args.perl.as_os_str(), OsStr::from_bytes(b"ref\xff.txt"));
    }
    assert_eq!(
        run_args(&["--perl=", "--rust", "r", "--out="]).perl,
        PathBuf::from(".")
    );
    assert_eq!(
        run_args(&["--perl=", "--rust", "r", "--out="]).out,
        PathBuf::from(".")
    );
}

#[test]
fn an_option_is_named_by_any_prefix_that_fits_it_alone() {
    let args = run_args(&[
        "--p",
        "p",
        "--r",
        "r",
        "--o",
        "o",
        "--t",
        "t",
        "--e",
        "2",
        "--k",
        "regulatory_feature_rows",
    ]);
    assert_eq!(args.perl, PathBuf::from("p"));
    assert_eq!(args.rust, PathBuf::from("r"));
    assert_eq!(args.out, PathBuf::from("o"));
    assert_eq!(args.tmp_dir, Some(PathBuf::from("t")));
    assert_eq!(args.examples, 2);
    assert_eq!(args.keep_only, Some(KeepOnly::RegulatoryFeatureRows));
    let args = run_args(&[
        "--pe=p",
        "--rus=r",
        "--ou=o",
        "--tmp-d=t",
        "--exam=3",
        "--keep=regulatory_feature_rows",
    ]);
    assert_eq!(args.perl, PathBuf::from("p"));
    assert_eq!(args.tmp_dir, Some(PathBuf::from("t")));
    assert_eq!(args.examples, 3);
    assert_eq!(args.keep_only, Some(KeepOnly::RegulatoryFeatureRows));
    assert_eq!(
        error_with(&["--exampl", "2", "--EXAMPLES", "3"]),
        "unrecognized arguments: --EXAMPLES 3"
    );
    assert_eq!(
        error_with(&["--perl-x", "v"]),
        "unrecognized arguments: --perl-x v"
    );
    assert_eq!(
        error(&["--=x", "--perl", "p", "--rust", "r", "--out", "o"]),
        "ambiguous option: --=x could match --help, --perl, --rust, --out, --examples, --tmp-dir, --keep-only"
    );
}

#[test]
fn a_later_value_replaces_an_earlier_one() {
    assert_eq!(
        run_args(&[&BASE[..], &["--examples", "1", "--examples", "2"]].concat()).examples,
        2
    );
    assert_eq!(
        run_args(&["--perl", "p", "--perl", "q", "--rust", "r", "--out", "o"]).perl,
        PathBuf::from("q")
    );
}

#[test]
fn examples_is_read_as_int_reads_a_count() {
    let examples = |value: &str| run_args(&[&BASE[..], &["--examples", value]].concat()).examples;
    assert_eq!(examples("0"), 0);
    assert_eq!(examples("-1"), -1);
    assert_eq!(examples("+2"), 2);
    assert_eq!(examples(" 3 "), 3);
    assert_eq!(examples("007"), 7);
    assert_eq!(examples("1_000"), 1000);
    assert_eq!(examples("99999999999999999999999999"), i64::MAX);
    assert_eq!(
        run_args(&[&BASE[..], &["--examples=-2"]].concat()).examples,
        -2
    );
    for value in ["x", "1.5", "", "_1", "1__0", "1_", "0x10", "1e2", "-", "+"] {
        assert_eq!(
            error_with(&["--examples", value]),
            format!("argument --examples: invalid int value: '{value}'")
        );
    }
    assert_eq!(
        error_with(&["--examples="]),
        "argument --examples: invalid int value: ''"
    );
}

#[test]
fn a_value_that_reads_as_an_option_is_not_taken() {
    let expected_one = |name: &str| format!("argument {name}: expected one argument");
    assert_eq!(error(&["--perl"]), expected_one("--perl"));
    assert_eq!(error_with(&["--examples"]), expected_one("--examples"));
    assert_eq!(error_with(&["--out"]), expected_one("--out"));
    assert_eq!(error_with(&["--tmp-dir"]), expected_one("--tmp-dir"));
    assert_eq!(
        error_with(&["--tmp-dir", "--examples", "3"]),
        expected_one("--tmp-dir")
    );
    assert_eq!(
        error(&["--perl", "--rust", "r", "--out", "o"]),
        expected_one("--perl")
    );
    assert_eq!(
        error(&["--perl", "-x", "--rust", "r", "--out", "o"]),
        expected_one("--perl")
    );
    assert_eq!(
        error(&["--perl", "--", "p", "--rust", "r", "--out", "o"]),
        expected_one("--perl")
    );
    assert_eq!(
        error_with(&["--examples", "--3"]),
        expected_one("--examples")
    );
    for value in ["-", "-1", "-1.5", "-.5", "-1_0", "-1.", "-1e5", "-x y", ""] {
        let args = run_args(&["--perl", value, "--rust", "r", "--out", "o"]);
        let expected = if value.is_empty() { "." } else { value };
        assert_eq!(args.perl, PathBuf::from(expected), "{value:?}");
    }
}

#[test]
fn missing_options_are_named_together_before_any_argument_left_unread() {
    let required = |names: &str| format!("the following arguments are required: {names}");
    assert_eq!(error(&[]), required("--perl, --rust, --out"));
    assert_eq!(error(&["--out", "o"]), required("--perl, --rust"));
    assert_eq!(error(&["--perl", "p", "--rust", "r"]), required("--out"));
    assert_eq!(
        error(&["--perl", "p", "--bogus"]),
        required("--rust, --out")
    );
    assert_eq!(
        error(&["--", "--perl", "p", "--rust", "r", "--out", "o"]),
        required("--perl, --rust, --out")
    );
}

#[test]
fn arguments_no_option_takes_are_named_together() {
    let unrecognized = |args: &str| format!("unrecognized arguments: {args}");
    assert_eq!(error_with(&["--bogus"]), unrecognized("--bogus"));
    assert_eq!(
        error_with(&["--bogus", "--other"]),
        unrecognized("--bogus --other")
    );
    assert_eq!(
        error_with(&["--bogus=1", "stray"]),
        unrecognized("--bogus=1 stray")
    );
    assert_eq!(error_with(&["stray"]), unrecognized("stray"));
    assert_eq!(error_with(&["-x"]), unrecognized("-x"));
    assert_eq!(error_with(&["-"]), unrecognized("-"));
    assert_eq!(error_with(&["---"]), unrecognized("---"));
    assert_eq!(error_with(&["--"]), unrecognized("--"));
    assert_eq!(error_with(&["--", "x"]), unrecognized("-- x"));
    assert_eq!(
        error_with(&["--", "--examples", "1"]),
        unrecognized("-- --examples 1")
    );
    assert_eq!(
        error_with(&["-x", "--examples", "q"]),
        "argument --examples: invalid int value: 'q'"
    );
}

#[test]
fn help_is_answered_where_it_stands_on_the_line() {
    assert!(is_help(&["-h"]));
    assert!(is_help(&["--help"]));
    assert!(is_help(&["--he"]));
    assert!(is_help(&["--h"]));
    assert!(is_help(&["-hx"]));
    assert!(is_help(&["-hh"]));
    assert!(is_help(&["-h", "--examples", "x"]));
    assert!(is_help(&[&BASE[..], &["--bogus", "-h"]].concat()));
    assert!(is_help(&[
        "--p", "p", "--r", "r", "--o", "o", "--e", "1", "--h"
    ]));
    assert_eq!(
        error_with(&["--examples", "x", "-h"]),
        "argument --examples: invalid int value: 'x'"
    );
    let refused = |value: &str| format!("argument -h/--help: ignored explicit argument '{value}'");
    assert_eq!(error(&["--help=x"]), refused("x"));
    assert_eq!(error_with(&["-h=3"]), refused("3"));
    assert_eq!(error_with(&["--he=3"]), refused("3"));
    assert_eq!(error(&["-hh=3"]), refused("3"));
    assert_eq!(error(&["-h-x"]), refused("-x"));
}
