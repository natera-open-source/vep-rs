// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line under `argparse`'s reading: each case is one command line and what
//! `argparse` makes of it for this parser.

use super::*;

fn parsed(words: &[&str]) -> Result<Command, String> {
    let argv: Vec<OsString> = words.iter().map(OsString::from).collect();
    parse(&argv)
}

fn run(input: &str, output: Option<&str>) -> Result<Command, String> {
    Ok(Command::Run(Args {
        input: OsString::from(input),
        output: output.map(OsString::from),
    }))
}

fn refused(message: &str) -> Result<Command, String> {
    Err(message.to_string())
}

const ONE_ARGUMENT: &str = "argument -o/--output: expected one argument";
const NO_INPUT: &str = "the following arguments are required: input";

#[test]
fn the_input_alone_leaves_the_output_to_its_default() {
    assert_eq!(parsed(&["in.tsv"]), run("in.tsv", None));
    assert_eq!(parsed(&["./a/../in.tsv/"]), run("./a/../in.tsv/", None));
}

#[test]
fn the_output_is_named_before_or_after_the_input_whole_or_by_a_prefix() {
    for words in [
        &["in.tsv", "-o", "out.tsv"][..],
        &["-o", "out.tsv", "in.tsv"],
        &["in.tsv", "--output", "out.tsv"],
        &["in.tsv", "--output=out.tsv"],
        &["in.tsv", "--out", "out.tsv"],
        &["in.tsv", "--o", "out.tsv"],
        &["in.tsv", "--outpu=out.tsv"],
        &["in.tsv", "-oout.tsv"],
        &["in.tsv", "-o=out.tsv"],
        &["in.tsv", "-o", "a.tsv", "-o", "out.tsv"],
        &["in.tsv", "--output", "a.tsv", "--output", "out.tsv"],
        &["--ou=x", "--o=out.tsv", "in.tsv"],
        &["-o", "x", "in.tsv", "-oout.tsv"],
        &["-o", "out.tsv", "--", "in.tsv"],
    ] {
        assert_eq!(parsed(words), run("in.tsv", Some("out.tsv")), "{words:?}");
    }
}

#[test]
fn a_one_letter_option_takes_the_rest_of_its_argument_as_the_value() {
    assert_eq!(parsed(&["-oh", "in.tsv"]), run("in.tsv", Some("h")));
    assert_eq!(parsed(&["-o-h", "in.tsv"]), run("in.tsv", Some("-h")));
    assert_eq!(parsed(&["-o--", "in.tsv"]), run("in.tsv", Some("--")));
    assert_eq!(parsed(&["-o x", "in.tsv"]), run("in.tsv", Some(" x")));
    assert_eq!(parsed(&["-ox=y", "in.tsv"]), run("in.tsv", Some("x=y")));
    assert_eq!(
        parsed(&["--output=--", "in.tsv"]),
        run("in.tsv", Some("--"))
    );
}

#[test]
fn an_empty_value_is_a_value() {
    assert_eq!(parsed(&["in.tsv", "-o", ""]), run("in.tsv", Some("")));
    assert_eq!(parsed(&["in.tsv", "-o="]), run("in.tsv", Some("")));
    assert_eq!(parsed(&["in.tsv", "--output="]), run("in.tsv", Some("")));
    assert_eq!(parsed(&[""]), run("", None));
}

#[test]
fn a_value_may_start_with_a_dash_where_it_reads_as_no_option() {
    assert_eq!(parsed(&["in.tsv", "-o", "-5"]), run("in.tsv", Some("-5")));
    assert_eq!(
        parsed(&["in.tsv", "-o", "-.5x"]),
        run("in.tsv", Some("-.5x"))
    );
    assert_eq!(
        parsed(&["in.tsv", "-o", "-\u{663}"]),
        run("in.tsv", Some("-\u{663}"))
    );
    assert_eq!(parsed(&["in.tsv", "-o", "- x"]), run("in.tsv", Some("- x")));
    assert_eq!(parsed(&["in.tsv", "-o", "-"]), run("in.tsv", Some("-")));
    assert_eq!(parsed(&["-5"]), run("-5", None));
    assert_eq!(parsed(&["-"]), run("-", None));
    assert_eq!(parsed(&["- x"]), run("- x", None));
    assert_eq!(parsed(&["in.tsv", "-o", "-x"]), refused(ONE_ARGUMENT));
    assert_eq!(parsed(&["in.tsv", "-o", "-..5"]), refused(ONE_ARGUMENT));
    assert_eq!(parsed(&["in.tsv", "-o", "-\u{b2}"]), refused(ONE_ARGUMENT));
}

#[test]
fn an_option_without_its_value_is_refused() {
    for words in [
        &["in.tsv", "-o"][..],
        &["-o"],
        &["--o"],
        &["--output"],
        &["in.tsv", "-o", "--"],
        &["in.tsv", "-o", "--", "out.tsv"],
        &["-o", "out.tsv", "-o"],
        &["-o", "--output", "in.tsv"],
        &["-o", "-h", "in.tsv"],
        &["-ho"],
    ] {
        assert_eq!(parsed(words), refused(ONE_ARGUMENT), "{words:?}");
    }
}

#[test]
fn the_marker_makes_every_later_argument_a_value() {
    assert_eq!(parsed(&["--", "in.tsv"]), run("in.tsv", None));
    assert_eq!(parsed(&["in.tsv", "--"]), run("in.tsv", None));
    assert_eq!(parsed(&["--", "-o"]), run("-o", None));
    assert_eq!(
        parsed(&["--", "in.tsv", "-o", "out.tsv"]),
        refused("unrecognized arguments: -o out.tsv")
    );
    assert_eq!(
        parsed(&["in.tsv", "--", "-o", "out.tsv"]),
        refused("unrecognized arguments: -o out.tsv")
    );
    assert_eq!(
        parsed(&["--", "--", "in.tsv"]),
        refused("unrecognized arguments: in.tsv")
    );
    assert_eq!(
        parsed(&["in.tsv", "--", "--"]),
        refused("unrecognized arguments: --")
    );
    assert_eq!(parsed(&["--"]), refused(NO_INPUT));
}

#[test]
fn a_missing_input_is_named_before_the_arguments_nothing_took() {
    assert_eq!(parsed(&[]), refused(NO_INPUT));
    assert_eq!(parsed(&["-o", "out.tsv"]), refused(NO_INPUT));
    assert_eq!(parsed(&["-x"]), refused(NO_INPUT));
}

#[test]
fn arguments_nothing_took_are_listed_in_the_order_they_were_met() {
    assert_eq!(
        parsed(&["in.tsv", "extra", "more"]),
        refused("unrecognized arguments: extra more")
    );
    assert_eq!(
        parsed(&["-x", "in.tsv"]),
        refused("unrecognized arguments: -x")
    );
    assert_eq!(
        parsed(&["--unknown=1", "in.tsv"]),
        refused("unrecognized arguments: --unknown=1")
    );
    assert_eq!(
        parsed(&["in.tsv", "-O", "out.tsv"]),
        refused("unrecognized arguments: -O out.tsv")
    );
    assert_eq!(
        parsed(&["---output", "x", "in.tsv"]),
        refused("unrecognized arguments: ---output in.tsv")
    );
    assert_eq!(
        parsed(&["--output-x", "x", "in.tsv"]),
        refused("unrecognized arguments: --output-x in.tsv")
    );
    assert_eq!(
        parsed(&["--output x", "in.tsv"]),
        refused("unrecognized arguments: in.tsv")
    );
}

#[test]
fn a_prefix_of_every_option_is_ambiguous() {
    assert_eq!(
        parsed(&["-=x", "in.tsv"]),
        refused("ambiguous option: -=x could match -h, --help, -o, --output")
    );
    assert_eq!(
        parsed(&["--=x", "in.tsv"]),
        refused("ambiguous option: --=x could match --help, --output")
    );
}

#[test]
fn help_is_shown_as_soon_as_it_is_read() {
    for words in [
        &["-h"][..],
        &["--help"],
        &["--he"],
        &["--h"],
        &["-h", "in.tsv"],
        &["in.tsv", "-h"],
        &["in.tsv", "x", "-h"],
        &["--help", "--unknown"],
        &["--unknown", "--help"],
        &["-hx"],
        &["-hh"],
        &["-h 1"],
        &["-ho", "out.tsv", "in.tsv"],
        &["-hoout.tsv", "in.tsv"],
        &["-ho=", "in.tsv"],
        &["-h", "-o"],
    ] {
        assert_eq!(parsed(words), Ok(Command::Help), "{words:?}");
    }
}

#[test]
fn help_takes_no_value() {
    let ignored = |value: &str| {
        refused(&format!(
            "argument -h/--help: ignored explicit argument '{value}'"
        ))
    };
    assert_eq!(parsed(&["-h=1"]), ignored("1"));
    assert_eq!(parsed(&["-h="]), ignored(""));
    assert_eq!(parsed(&["--help=1"]), ignored("1"));
    assert_eq!(parsed(&["--help="]), ignored(""));
    assert_eq!(parsed(&["--he=1"]), ignored("1"));
    assert_eq!(parsed(&["-h-"]), ignored("-"));
    assert_eq!(parsed(&["-hh=1"]), ignored("1"));
    assert_eq!(parsed(&["-hh-x"]), ignored("-x"));
}

#[test]
fn an_error_met_before_the_help_is_reported_and_one_behind_it_is_not() {
    assert_eq!(parsed(&["-o", "-h"]), refused(ONE_ARGUMENT));
    assert_eq!(parsed(&["-h", "-o"]), Ok(Command::Help));
    assert_eq!(
        parsed(&["--=x", "-h"]),
        refused("ambiguous option: --=x could match --help, --output")
    );
}
