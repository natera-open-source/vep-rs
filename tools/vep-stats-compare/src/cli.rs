// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line, read as Python's `argparse` reads five one-value options and `-h`:
//! an option by its full name or any unambiguous prefix of it, its value as the next
//! argument or after `=`, a value that starts with `-` taken only where it reads as a
//! negative number or carries a space, the last of a repeated option kept, `--mask-row`
//! collected, help where `-h` or `--help` is met before the first fault.

use crate::pyfmt::{py_path_str, PyInt};
use crate::unicode::decimal_value;

pub const USAGE: &str = "usage: vep-stats-compare [-h] --expected EXPECTED --actual ACTUAL [--tolerance TOLERANCE] [--mask-row LABEL] [--report REPORT]";
pub const HELP: &str = "\
Section-by-section comparison of two VEP statistics reports: exit 0 on green, 1 on red, 2 on an error.

options:
  -h, --help             show this help message and exit
  --expected EXPECTED    the Ensembl VEP statistics report
  --actual ACTUAL        the candidate's statistics report
  --tolerance TOLERANCE  the count difference a chart row may show (default 0: the sets must be equal)
  --mask-row LABEL       a further row label whose value is set aside (repeatable)
  --report REPORT        where the JSON report is written";

/// The options of one run. A path is the bytes of `str(pathlib.Path(value))`.
#[derive(Debug, PartialEq, Eq)]
pub struct Args {
    pub expected: Vec<u8>,
    pub actual: Vec<u8>,
    pub tolerance: PyInt,
    pub mask_rows: Vec<Vec<u8>>,
    pub report: Option<Vec<u8>>,
}

/// What the command line asks for.
#[derive(Debug, PartialEq, Eq)]
pub enum Parsed {
    Run(Args),
    Help,
    /// A usage error, with its message.
    Usage(String),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Action {
    Help,
    Expected,
    Actual,
    Tolerance,
    MaskRow,
    Report,
}

/// The option strings in the order the parser declares them.
const OPTIONS: [(&str, Action); 7] = [
    ("-h", Action::Help),
    ("--help", Action::Help),
    ("--expected", Action::Expected),
    ("--actual", Action::Actual),
    ("--tolerance", Action::Tolerance),
    ("--mask-row", Action::MaskRow),
    ("--report", Action::Report),
];

fn action_of(option: &[u8]) -> Option<Action> {
    OPTIONS
        .iter()
        .find(|(name, _)| name.as_bytes() == option)
        .map(|(_, action)| *action)
}

/// What joins an option to a value written in the same argument.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Joint {
    /// `--option=value`.
    Equals,
    /// `-hvalue`.
    Nothing,
}

/// One reading of an argument as an option: the action (None for an option this parser
/// does not have), the option string and the value written in the same argument.
#[derive(Clone, Debug)]
struct Reading {
    action: Option<Action>,
    option: Vec<u8>,
    explicit: Option<(Joint, Vec<u8>)>,
}

fn split_at_equals(arg: &[u8]) -> (&[u8], Option<&[u8]>) {
    match arg.iter().position(|b| *b == b'=') {
        Some(at) => (&arg[..at], Some(&arg[at + 1..])),
        None => (arg, None),
    }
}

/// Every option an argument that is no option string itself may stand for: a long option
/// by a prefix of its name, `-h` by the argument's first two characters.
fn readings_by_prefix(arg: &[u8]) -> Vec<Reading> {
    let (prefix, value) = split_at_equals(arg);
    let explicit = value.map(|value| (Joint::Equals, value.to_vec()));
    let mut readings = Vec::new();
    for (name, action) in OPTIONS {
        let name = name.as_bytes();
        if arg[1] != b'-' && name == &arg[..2] {
            readings.push(Reading {
                action: Some(action),
                option: name.to_vec(),
                explicit: Some((Joint::Nothing, arg[2..].to_vec())),
            });
        } else if name.starts_with(prefix) {
            readings.push(Reading {
                action: Some(action),
                option: name.to_vec(),
                explicit: explicit.clone(),
            });
        }
    }
    readings
}

/// Whether the argument opens like a negative number: `-`, an optional `.`, a decimal
/// digit of any script.
fn opens_like_a_negative_number(arg: &[u8]) -> bool {
    let rest = arg.strip_prefix(b"-").unwrap_or(arg);
    let rest = rest.strip_prefix(b".").unwrap_or(rest);
    let first = String::from_utf8_lossy(&rest[..rest.len().min(4)])
        .chars()
        .next();
    first.is_some_and(|c| decimal_value(c).is_some())
}

/// The readings of an argument as an option, or None where it is a value.
fn readings(arg: &[u8]) -> Option<Vec<Reading>> {
    if arg.first() != Some(&b'-') {
        return None;
    }
    let exact = |option: &[u8], explicit: Option<(Joint, Vec<u8>)>| {
        action_of(option).map(|action| {
            vec![Reading {
                action: Some(action),
                option: option.to_vec(),
                explicit,
            }]
        })
    };
    if let Some(found) = exact(arg, None) {
        return Some(found);
    }
    if arg.len() == 1 {
        return None;
    }
    if let (option, Some(value)) = split_at_equals(arg) {
        if let Some(found) = exact(option, Some((Joint::Equals, value.to_vec()))) {
            return Some(found);
        }
    }
    let by_prefix = readings_by_prefix(arg);
    if !by_prefix.is_empty() {
        return Some(by_prefix);
    }
    if opens_like_a_negative_number(arg) || arg.contains(&b' ') {
        return None;
    }
    Some(vec![Reading {
        action: None,
        option: arg.to_vec(),
        explicit: None,
    }])
}

/// One argument's place in the command line.
enum Kind {
    Value,
    /// The first `--`; every argument after it is a value.
    Separator,
    Option(Vec<Reading>),
}

fn shown(arg: &[u8]) -> String {
    String::from_utf8_lossy(arg).into_owned()
}

/// Reads the arguments after the program's name.
pub fn parse(argv: &[Vec<u8>]) -> Parsed {
    let mut kinds = Vec::with_capacity(argv.len());
    let mut separated = false;
    for arg in argv {
        kinds.push(if separated {
            Kind::Value
        } else if arg == b"--" {
            separated = true;
            Kind::Separator
        } else {
            readings(arg).map_or(Kind::Value, Kind::Option)
        });
    }
    let mut expected = None;
    let mut actual = None;
    let mut tolerance = PyInt::from_u64(0);
    let mut mask_rows = Vec::new();
    let mut report = None;
    // No argument is positional, so a value no option takes is unrecognized; the fault is
    // raised after every option is read, as is a missing required option.
    let mut unrecognized: Vec<String> = Vec::new();
    let mut at = 0;
    while at < argv.len() {
        let Kind::Option(readings) = &kinds[at] else {
            unrecognized.push(shown(&argv[at]));
            at += 1;
            continue;
        };
        if readings.len() > 1 {
            let matches: Vec<String> = readings.iter().map(|r| shown(&r.option)).collect();
            return Parsed::Usage(format!(
                "ambiguous option: {} could match {}",
                shown(&argv[at]),
                matches.join(", ")
            ));
        }
        let first = &readings[0];
        let mut option: &[u8] = &first.option;
        let mut explicit: Option<(Joint, &[u8])> = first
            .explicit
            .as_ref()
            .map(|(joint, value)| (*joint, value.as_slice()));
        let mut help = false;
        let mut stored: Option<(Action, &[u8])> = None;
        let mut next = at + 1;
        match first.action {
            None => unrecognized.push(shown(&argv[at])),
            Some(Action::Help) => loop {
                match explicit {
                    None => {
                        help = true;
                        break;
                    }
                    // `-h` run together with more: each further character is a flag of
                    // its own, and no other flag exists.
                    Some((joint, tail)) if option == b"-h" && !tail.is_empty() => {
                        if joint == Joint::Equals || tail[0] == b'-' {
                            return Parsed::Usage(format!(
                                "argument -h/--help: ignored explicit argument '{}'",
                                shown(tail)
                            ));
                        }
                        help = true;
                        if tail[0] != b'h' {
                            unrecognized.push(format!("-{}", shown(tail)));
                            break;
                        }
                        option = b"-h";
                        explicit = match &tail[1..] {
                            [] => None,
                            [b'=', value @ ..] => Some((Joint::Equals, value)),
                            value => Some((Joint::Nothing, value)),
                        };
                    }
                    Some((_, value)) => {
                        return Parsed::Usage(format!(
                            "argument -h/--help: ignored explicit argument '{}'",
                            shown(value)
                        ));
                    }
                }
            },
            Some(action) => match explicit {
                Some((_, value)) => stored = Some((action, value)),
                None if matches!(kinds.get(next), Some(Kind::Value)) => {
                    stored = Some((action, &argv[next]));
                    next += 1;
                }
                None => {
                    return Parsed::Usage(format!(
                        "argument {}: expected one argument",
                        shown(option)
                    ));
                }
            },
        }
        if help {
            return Parsed::Help;
        }
        if let Some((action, value)) = stored {
            match action {
                Action::Expected => expected = Some(py_path_str(value)),
                Action::Actual => actual = Some(py_path_str(value)),
                Action::Report => report = Some(py_path_str(value)),
                Action::MaskRow => mask_rows.push(value.to_vec()),
                Action::Tolerance => {
                    let read = std::str::from_utf8(value).ok().and_then(PyInt::parse);
                    let Some(read) = read else {
                        return Parsed::Usage(format!(
                            "argument --tolerance: invalid int value: '{}'",
                            shown(value)
                        ));
                    };
                    tolerance = read;
                }
                Action::Help => unreachable!("help stores no value"),
            }
        }
        at = next;
    }
    let (expected, actual) = match (expected, actual) {
        (Some(expected), Some(actual)) => (expected, actual),
        (expected, actual) => {
            let mut missing = Vec::new();
            if expected.is_none() {
                missing.push("--expected");
            }
            if actual.is_none() {
                missing.push("--actual");
            }
            return Parsed::Usage(format!(
                "the following arguments are required: {}",
                missing.join(", ")
            ));
        }
    };
    if !unrecognized.is_empty() {
        return Parsed::Usage(format!(
            "unrecognized arguments: {}",
            unrecognized.join(" ")
        ));
    }
    Parsed::Run(Args {
        expected,
        actual,
        tolerance,
        mask_rows,
        report,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_strs<S: AsRef<str>>(args: &[S]) -> Parsed {
        let argv: Vec<Vec<u8>> = args
            .iter()
            .map(|a| a.as_ref().as_bytes().to_vec())
            .collect();
        parse(&argv)
    }

    fn run<S: AsRef<str>>(args: &[S]) -> Args {
        match parse_strs(args) {
            Parsed::Run(args) => args,
            other => panic!("not a run: {other:?}"),
        }
    }

    fn is_usage<S: AsRef<str>>(args: &[S]) -> bool {
        matches!(parse_strs(args), Parsed::Usage(_))
    }

    const E: [&str; 4] = ["--expected", "e", "--actual", "a"];

    /// The two required options, then `extra`.
    fn with(extra: &[&str]) -> Vec<String> {
        E.iter().chain(extra).map(|arg| arg.to_string()).collect()
    }

    #[test]
    fn options_take_both_spellings_any_prefix_and_their_defaults() {
        let args = run(&E);
        assert_eq!(args.expected, b"e");
        assert_eq!(args.actual, b"a");
        assert_eq!(args.tolerance, PyInt::from_u64(0));
        assert!(args.mask_rows.is_empty());
        assert_eq!(args.report, None);
        let args = run(&[
            "--exp",
            "./x//e/",
            "--a=a/.",
            "--tol",
            "2",
            "--mask",
            "Row one",
            "--mask-ro=Row=two",
            "--rep",
            "out/r.json",
        ]);
        assert_eq!(args.expected, b"x/e");
        assert_eq!(args.actual, b"a");
        assert_eq!(args.tolerance, PyInt::from_u64(2));
        assert_eq!(args.mask_rows, [b"Row one".to_vec(), b"Row=two".to_vec()]);
        assert_eq!(args.report.as_deref(), Some(&b"out/r.json"[..]));
        // The last of a repeated option is the one kept.
        let args = run(&with(&["--expected", "z", "--tolerance", "1", "--t=3"]));
        assert_eq!(args.expected, b"z");
        assert_eq!(args.tolerance, PyInt::from_u64(3));
        let raw = [
            b"--expected=./e\xff/".to_vec(),
            b"--actual".to_vec(),
            b"a".to_vec(),
            b"--mask-row".to_vec(),
            b"\xfflabel".to_vec(),
        ];
        let Parsed::Run(args) = parse(&raw) else {
            panic!("raw bytes are a path and a label");
        };
        assert_eq!(args.expected, b"e\xff");
        assert_eq!(args.mask_rows, [b"\xfflabel".to_vec()]);
    }

    #[test]
    fn a_value_that_starts_with_a_dash_is_taken_as_argparse_takes_it() {
        for taken in [
            "-5",
            "-x y",
            "-",
            "",
            "-.5x",
            "--a b",
            "--expected x",
            "-\u{663}",
        ] {
            assert_eq!(
                run(&with(&["--mask-row", taken])).mask_rows,
                [taken.as_bytes().to_vec()],
                "{taken:?}"
            );
        }
        for refused in ["-x", "--foo", "--", "--tolerance", "-h", "--he"] {
            assert!(
                !matches!(parse_strs(&with(&["--mask-row", refused])), Parsed::Run(_)),
                "{refused:?}"
            );
        }
        assert_eq!(run(&with(&["--mask-row=-x"])).mask_rows, [b"-x".to_vec()]);
        assert_eq!(run(&with(&["--mask-row="])).mask_rows, [b"".to_vec()]);
        assert_eq!(run(&with(&["--expected", "-5"])).expected, b"-5");
    }

    #[test]
    fn tolerance_is_an_integer_as_int_reads_one() {
        for (text, expected) in [
            ("1", "1"),
            (" 2 ", "2"),
            ("+2", "2"),
            ("1_0", "10"),
            ("\u{663}", "3"),
            ("-0", "0"),
            ("-1", "-1"),
            ("-1 ", "-1"),
            (" -1", "-1"),
            ("-\u{663}", "-3"),
            (
                "99999999999999999999999999999999999999999",
                "99999999999999999999999999999999999999999",
            ),
        ] {
            assert_eq!(
                run(&with(&["--tolerance", text])).tolerance.decimal(),
                expected,
                "{text:?}"
            );
        }
        assert_eq!(run(&with(&["--tolerance=-1"])).tolerance.decimal(), "-1");
        for refused in [
            "x", "", "1.0", "0x10", "1e3", "-1.5", "-1e3", "- 1", "-x", "-", "-.5", "_1", "1_",
            "1__0", "\u{1c}1", "+-1",
        ] {
            assert!(is_usage(&with(&["--tolerance", refused])), "{refused:?}");
        }
        assert!(is_usage(&with(&["--tolerance=1=2"])));
        assert!(is_usage(&with(&["--tolerance"])));
        assert!(is_usage(&with(&["--tolerance", "--mask-row", "x"])));
        let raw = [
            b"--expected=e".to_vec(),
            b"--actual=a".to_vec(),
            b"--tolerance".to_vec(),
            b"1\xff".to_vec(),
        ];
        assert!(matches!(parse(&raw), Parsed::Usage(_)));
    }

    #[test]
    fn faults_are_usage_errors() {
        assert!(is_usage::<&str>(&[]));
        assert!(is_usage(&["--expected", "e"]));
        assert!(is_usage(&["--actual", "a"]));
        assert!(is_usage(&["-e", "e", "-a", "a"]));
        assert!(is_usage(&["--expected", "--actual", "a"]));
        for extra in [
            &["extra"][..],
            &["--", "extra"],
            &["--"],
            &["--bogus"],
            &["--mask-rowx", "a"],
            &["--tolerancee", "1"],
            &["-t", "1"],
            &["--=x"],
            &["-=x"],
            &["-"],
            &["--tolerance", "1", "2"],
            &["--report"],
            &["--expected"],
            &["--help=x"],
            &["-h=x"],
            &["-h="],
            &["-h-x"],
        ] {
            assert!(is_usage(&with(extra)), "{extra:?}");
        }
    }

    #[test]
    fn help_wins_where_it_is_met_before_a_fault() {
        for args in [
            &["-h"][..],
            &["--help"],
            &["--he"],
            &["--h"],
            &["-hh"],
            &["-hx"],
            &["-h5"],
            &["--bogus", "-h"],
            &["extra", "--help"],
            &["-h", "--tolerance", "x"],
            &["-h", "--"],
            &["--expected", "e", "-h"],
        ] {
            assert_eq!(parse_strs(args), Parsed::Help, "{args:?}");
        }
        assert!(is_usage(&["--tolerance", "x", "-h"]));
        assert!(is_usage(&["--", "-h"]));
        assert!(is_usage(&["--mask-row", "-h"]));
    }
}
