// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line as CPython's `argparse` reads it: every argument is classed as an
//! option or a value before any is consumed, an option may be shortened to a prefix that
//! fits it alone, and a value follows its option or an `=`.

use std::ffi::{OsStr, OsString};
use std::os::unix::ffi::OsStrExt;

pub const USAGE: &str = "usage: vep-cache-diff [-h] --baseline-dir BASELINE_DIR --builder-dir BUILDER_DIR --output-report OUTPUT_REPORT";
pub const HELP: &str = "\
Structurally compare two JSON cache directories.

options:
  -h, --help            show this help message and exit
  --baseline-dir BASELINE_DIR
                        Path to the baseline JSON cache directory
  --builder-dir BUILDER_DIR
                        Path to the builder-generated JSON cache directory
  --output-report OUTPUT_REPORT
                        Path for the output JSON diff report

environment:
  VEP_CACHE_DIFF_THREADS       worker threads (default: the host's available cores)
  VEP_CACHE_DIFF_BUFFER_BYTES  bytes of variation records a chromosome holds in memory
                               before it writes them to a temporary file (default:
                               268435456)";

#[derive(Debug, PartialEq)]
pub struct Args {
    pub baseline_dir: OsString,
    pub builder_dir: OsString,
    pub output_report: OsString,
}

#[derive(Debug, PartialEq)]
pub enum Parsed {
    Run(Args),
    Help,
    Usage(String),
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Opt {
    Help,
    BaselineDir,
    BuilderDir,
    OutputReport,
}

const LONG_OPTIONS: [(&str, Opt); 4] = [
    ("--help", Opt::Help),
    ("--baseline-dir", Opt::BaselineDir),
    ("--builder-dir", Opt::BuilderDir),
    ("--output-report", Opt::OutputReport),
];

/// What one argument is before it is used, by argparse's rules.
#[derive(Debug, PartialEq)]
enum Token<'a> {
    /// A value: an option takes it, and on its own it is an argument too many.
    Value,
    /// An option, with the text behind its `=` when it has one.
    Option(Opt, Option<&'a [u8]>),
    /// `-h` with letters behind it, which argparse reads as further one-letter options.
    HelpCluster(&'a [u8]),
    /// A prefix of more than one long option, with the options it fits.
    Ambiguous(Vec<&'static str>),
    /// A dash-led argument that names no option.
    Unknown,
}

/// Classifies one argument. An option is one of the five names, a prefix of one long
/// name alone, or either with `=value` behind it. An argument that is empty, `-`, dash-led
/// like a negative number (`-5`, `-.5`) or holds a space is a value.
fn classify(arg: &[u8]) -> Token<'_> {
    if arg.first() != Some(&b'-') || arg.len() == 1 {
        return Token::Value;
    }
    if arg == b"-h" {
        return Token::Option(Opt::Help, None);
    }
    let (name, explicit) = match memchr::memchr(b'=', arg) {
        Some(i) => (&arg[..i], Some(&arg[i + 1..])),
        None => (arg, None),
    };
    if name == b"-h" {
        return Token::Option(Opt::Help, explicit);
    }
    if arg.starts_with(b"--") {
        if let Some((_, opt)) = LONG_OPTIONS
            .iter()
            .find(|(long, _)| long.as_bytes() == name)
        {
            return Token::Option(*opt, explicit);
        }
        let fitting: Vec<&(&str, Opt)> = LONG_OPTIONS
            .iter()
            .filter(|(long, _)| long.as_bytes().starts_with(name))
            .collect();
        match fitting[..] {
            [] => {}
            [(_, opt)] => return Token::Option(*opt, explicit),
            _ => return Token::Ambiguous(fitting.iter().map(|(long, _)| *long).collect()),
        }
    } else if arg.starts_with(b"-h") {
        return Token::HelpCluster(&arg[2..]);
    }
    let digit = |i: usize| arg.get(i).is_some_and(u8::is_ascii_digit);
    if digit(1) || (arg[1] == b'.' && digit(2)) || arg.contains(&b' ') {
        return Token::Value;
    }
    Token::Unknown
}

/// Reads the command line as argparse reads it, left to right: `--help` ends the reading
/// where it stands, a missing value is an error where it stands, and an argument too many
/// or a missing required option is an error at the end.
pub fn parse(argv: &[OsString]) -> Parsed {
    let mut baseline_dir = None;
    let mut builder_dir = None;
    let mut output_report = None;
    let mut surplus: Vec<&OsStr> = Vec::new();
    let mut rest = argv.iter().map(OsString::as_os_str);
    while let Some(arg) = rest.next() {
        if arg == "--" {
            // It and everything behind it are values, and no option is left to take one.
            surplus.push(arg);
            surplus.extend(rest.by_ref());
            break;
        }
        let (opt, explicit) = match classify(arg.as_bytes()) {
            Token::Option(opt, explicit) => (opt, explicit),
            Token::HelpCluster(letters) if letters[0] != b'-' => return Parsed::Help,
            Token::HelpCluster(letters) => {
                return Parsed::Usage(format!(
                    "argument -h/--help: ignored explicit argument '{}'",
                    String::from_utf8_lossy(letters)
                ));
            }
            Token::Ambiguous(fitting) => {
                return Parsed::Usage(format!(
                    "ambiguous option: {} could match {}",
                    arg.to_string_lossy(),
                    fitting.join(", ")
                ));
            }
            Token::Value | Token::Unknown => {
                surplus.push(arg);
                continue;
            }
        };
        let name = LONG_OPTIONS
            .iter()
            .find(|(_, o)| *o == opt)
            .map_or("", |(long, _)| long);
        if opt == Opt::Help {
            return match explicit {
                None => Parsed::Help,
                Some(text) => Parsed::Usage(format!(
                    "argument -h/--help: ignored explicit argument '{}'",
                    String::from_utf8_lossy(text)
                )),
            };
        }
        let value = match explicit {
            Some(text) => OsStr::from_bytes(text),
            None => {
                let next = rest.clone().next();
                match next {
                    Some(next) if next != "--" && classify(next.as_bytes()) == Token::Value => {
                        rest.next();
                        next
                    }
                    _ => return Parsed::Usage(format!("argument {name}: expected one argument")),
                }
            }
        };
        match opt {
            Opt::BaselineDir => baseline_dir = Some(value.to_os_string()),
            Opt::BuilderDir => builder_dir = Some(value.to_os_string()),
            Opt::OutputReport => output_report = Some(value.to_os_string()),
            Opt::Help => unreachable!("--help is handled above"),
        }
    }
    let (Some(baseline_dir), Some(builder_dir), Some(output_report)) =
        (&baseline_dir, &builder_dir, &output_report)
    else {
        let missing: Vec<&str> = [
            (&baseline_dir, "--baseline-dir"),
            (&builder_dir, "--builder-dir"),
            (&output_report, "--output-report"),
        ]
        .iter()
        .filter(|(value, _)| value.is_none())
        .map(|(_, name)| *name)
        .collect();
        return Parsed::Usage(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    };
    if !surplus.is_empty() {
        let listed: Vec<_> = surplus.iter().map(|arg| arg.to_string_lossy()).collect();
        return Parsed::Usage(format!("unrecognized arguments: {}", listed.join(" ")));
    }
    Parsed::Run(Args {
        baseline_dir: baseline_dir.clone(),
        builder_dir: builder_dir.clone(),
        output_report: output_report.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parsed(args: &[&str]) -> Parsed {
        parse(&args.iter().map(OsString::from).collect::<Vec<_>>())
    }

    fn run_with(baseline_dir: &str, builder_dir: &str, output_report: &str) -> Parsed {
        Parsed::Run(Args {
            baseline_dir: baseline_dir.into(),
            builder_dir: builder_dir.into(),
            output_report: output_report.into(),
        })
    }

    fn usage(args: &[&str]) -> String {
        match parsed(args) {
            Parsed::Usage(message) => message,
            other => panic!("{args:?}: {other:?}"),
        }
    }

    const ALL: [&str; 6] = [
        "--baseline-dir",
        "a",
        "--builder-dir",
        "b",
        "--output-report",
        "o",
    ];

    fn all_and(more: &[&'static str]) -> Vec<&'static str> {
        ALL.iter().chain(more).copied().collect()
    }

    #[test]
    fn options_take_both_spellings_and_a_later_value_replaces_an_earlier_one() {
        assert_eq!(parsed(&ALL), run_with("a", "b", "o"));
        assert_eq!(
            parsed(&[
                "--output-report=o",
                "--builder-dir",
                "b",
                "--baseline-dir=a"
            ]),
            run_with("a", "b", "o")
        );
        assert_eq!(
            parsed(&all_and(&["--baseline-dir", "later", "--output-report="])),
            run_with("later", "b", "")
        );
    }

    #[test]
    fn a_prefix_names_the_one_option_it_fits() {
        assert_eq!(
            parsed(&["--ba", "a", "--bu=b", "--o", "o"]),
            run_with("a", "b", "o")
        );
        assert_eq!(
            parsed(&["--baseline", "a", "--builder-d", "b", "--output-repor=o"]),
            run_with("a", "b", "o")
        );
        assert_eq!(parsed(&["--h"]), Parsed::Help);
        assert_eq!(
            usage(&["--b", "a", "--builder-dir", "b", "--output-report", "o"]),
            "ambiguous option: --b could match --baseline-dir, --builder-dir"
        );
        assert_eq!(
            usage(&["--b=a"]),
            "ambiguous option: --b=a could match --baseline-dir, --builder-dir"
        );
        assert_eq!(
            usage(&["--=x"]),
            "ambiguous option: --=x could match --help, --baseline-dir, --builder-dir, --output-report"
        );
    }

    #[test]
    fn values_follow_argparse_s_reading_of_dash_led_arguments() {
        assert_eq!(
            parsed(&[
                "--baseline-dir",
                "-",
                "--builder-dir",
                "-5",
                "--output-report",
                "-.5x"
            ]),
            run_with("-", "-5", "-.5x")
        );
        assert_eq!(
            parsed(&[
                "--baseline-dir",
                "-my dir",
                "--builder-dir",
                "",
                "--output-report=-x"
            ]),
            run_with("-my dir", "", "-x")
        );
        assert_eq!(
            parsed(&[
                "--baseline-dir=--builder-dir",
                "--builder-dir=b",
                "--output-report=o"
            ]),
            run_with("--builder-dir", "b", "o")
        );
        for (args, option) in [
            (
                &["--baseline-dir", "-x", "--builder-dir", "b"][..],
                "--baseline-dir",
            ),
            (&["--baseline-dir", "--builder-dir", "b"], "--baseline-dir"),
            (&["--baseline-dir", "a", "--builder-dir"], "--builder-dir"),
            (&["--output-report", "--", "o"], "--output-report"),
            (&["--output-report", "-.x"], "--output-report"),
        ] {
            assert_eq!(
                usage(args),
                format!("argument {option}: expected one argument")
            );
        }
    }

    #[test]
    fn help_ends_the_reading_where_it_stands() {
        assert_eq!(parsed(&["-h"]), Parsed::Help);
        assert_eq!(parsed(&["--help"]), Parsed::Help);
        assert_eq!(parsed(&["--bogus", "extra", "--help"]), Parsed::Help);
        assert_eq!(parsed(&["--baseline-dir", "a", "-h", "--b"]), Parsed::Help);
        assert_eq!(parsed(&["-hx"]), Parsed::Help);
        for (args, text) in [
            (&["--help=x"][..], "x"),
            (&["--help="], ""),
            (&["-h=x"], "x"),
            (&["-h-x"], "-x"),
        ] {
            assert_eq!(
                usage(args),
                format!("argument -h/--help: ignored explicit argument '{text}'")
            );
        }
        assert!(matches!(parsed(&["--b", "--help"]), Parsed::Usage(_)));
        assert!(matches!(
            parsed(&["--baseline-dir", "--help"]),
            Parsed::Usage(_)
        ));
        assert!(matches!(parsed(&["--", "--help"]), Parsed::Usage(_)));
    }

    #[test]
    fn a_missing_option_is_named_before_an_argument_too_many() {
        assert_eq!(
            usage(&[]),
            "the following arguments are required: --baseline-dir, --builder-dir, --output-report"
        );
        assert_eq!(
            usage(&["--builder-dir", "b", "extra"]),
            "the following arguments are required: --baseline-dir, --output-report"
        );
        assert_eq!(
            usage(&[
                "--baseline-dir",
                "a",
                "--",
                "--builder-dir",
                "b",
                "--output-report",
                "o"
            ]),
            "the following arguments are required: --builder-dir, --output-report"
        );
        assert_eq!(
            usage(&[
                "-baseline-dir",
                "a",
                "--builder-dir",
                "b",
                "--output-report",
                "o"
            ]),
            "the following arguments are required: --baseline-dir"
        );
        for (more, listed) in [
            (&["extra", "more"][..], "extra more"),
            (&["--bogus"], "--bogus"),
            (&["--bogus=1"], "--bogus=1"),
            (&["-x"], "-x"),
            (&["--"], "--"),
            (&["--", "x", "--help"], "-- x --help"),
            (&[""], ""),
            (&["-"], "-"),
            (&["--output-reports", "o"], "--output-reports o"),
        ] {
            assert_eq!(
                usage(&all_and(more)),
                format!("unrecognized arguments: {listed}"),
                "{more:?}"
            );
        }
        assert_eq!(
            usage(&[
                "first",
                "--baseline-dir",
                "a",
                "--builder-dir",
                "b",
                "--output-report",
                "o"
            ]),
            "unrecognized arguments: first"
        );
    }

    #[test]
    fn values_that_are_not_utf8_reach_their_option_as_given() {
        let raw = OsStr::from_bytes(b"dir\xff/cache");
        let mut inline = OsString::from("--baseline-dir=");
        inline.push(raw);
        for argv in [
            vec![inline, "--builder-dir=b".into(), "--output-report=o".into()],
            vec![
                "--baseline-dir".into(),
                raw.to_os_string(),
                "--builder-dir=b".into(),
                "--output-report=o".into(),
            ],
        ] {
            match parse(&argv) {
                Parsed::Run(args) => assert_eq!(args.baseline_dir, raw),
                other => panic!("{other:?}"),
            }
        }
    }
}
