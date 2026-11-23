// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line as CPython's `argparse` reads it: every argument is classed as an
//! option or a value before any is consumed, an option may be shortened to a prefix that
//! fits it alone, and a value follows its option or an `=`.

use memchr::memchr;
use std::ffi::{OsStr, OsString};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::PathBuf;

pub const USAGE: &str = "\
usage: vep-intermediates-compare [-h] --perl-trace PERL_TRACE --rust-output RUST_OUTPUT
                                 [--discordant-only] [--output OUTPUT]";
pub const HELP: &str = "\
Joins an Ensembl VEP intermediate-value trace with a vep-rs default-format output on
(location, allele, transcript) and reports the fields that diverge.

options:
  -h, --help                 show this help message and exit
  --perl-trace PERL_TRACE    tab-separated trace of trace_perl_intermediates.pl
  --rust-output RUST_OUTPUT  vep-rs default-format output
  --discordant-only          print only the rows with a diverging field
  --output OUTPUT            report file (default: -, standard output)";

#[derive(Debug, PartialEq)]
pub struct Args {
    pub perl_trace: PathBuf,
    pub rust_output: PathBuf,
    pub discordant_only: bool,
    /// The report path as given; `-` is standard output.
    pub output: OsString,
}

#[derive(Debug, PartialEq)]
pub enum Parsed {
    Run(Args),
    Help,
}

#[derive(Clone, Copy, PartialEq)]
enum Opt {
    Help,
    PerlTrace,
    RustOutput,
    DiscordantOnly,
    Output,
}

impl Opt {
    fn takes_value(self) -> bool {
        matches!(self, Opt::PerlTrace | Opt::RustOutput | Opt::Output)
    }
}

const SHORT_HELP: &[u8] = b"-h";
const LONG: [(&[u8], Opt); 5] = [
    (b"--help", Opt::Help),
    (b"--perl-trace", Opt::PerlTrace),
    (b"--rust-output", Opt::RustOutput),
    (b"--discordant-only", Opt::DiscordantOnly),
    (b"--output", Opt::Output),
];

/// One reading of an option argument: the option, whether it is spelled `-h`, and the
/// text attached to it with whether an `=` introduced that text.
struct Candidate {
    opt: Opt,
    short: bool,
    attached: Option<(bool, Vec<u8>)>,
}

enum Token {
    /// Not an option: a value for the option before it, or an argument nothing takes.
    Value,
    /// The `--` that makes every later argument a value.
    Separator,
    /// Shaped like an option and naming none.
    Unknown,
    Option(Vec<Candidate>),
}

fn exact(name: &[u8]) -> Option<Opt> {
    if name == SHORT_HELP {
        return Some(Opt::Help);
    }
    LONG.iter()
        .find(|(long, _)| *long == name)
        .map(|&(_, opt)| opt)
}

fn classify(arg: &[u8]) -> Token {
    if arg.first() != Some(&b'-') {
        return Token::Value;
    }
    if let Some(opt) = exact(arg) {
        return Token::Option(vec![Candidate {
            opt,
            short: arg == SHORT_HELP,
            attached: None,
        }]);
    }
    if arg.len() == 1 {
        return Token::Value;
    }
    let (name, value) = match memchr(b'=', arg) {
        Some(i) => (&arg[..i], Some(&arg[i + 1..])),
        None => (arg, None),
    };
    if let (Some(opt), Some(value)) = (exact(name), value) {
        return Token::Option(vec![Candidate {
            opt,
            short: name == SHORT_HELP,
            attached: Some((true, value.to_vec())),
        }]);
    }
    // `-=...` reads as `--=...` does: `-` is a prefix of every option.
    let candidates: Vec<Candidate> = if arg[1] == b'-' || name == b"-" {
        LONG.iter()
            .filter(|(long, _)| long.starts_with(name))
            .map(|&(_, opt)| Candidate {
                opt,
                short: false,
                attached: value.map(|v| (true, v.to_vec())),
            })
            .collect()
    } else if arg.starts_with(SHORT_HELP) {
        vec![Candidate {
            opt: Opt::Help,
            short: true,
            attached: Some((false, arg[2..].to_vec())),
        }]
    } else {
        Vec::new()
    };
    if !candidates.is_empty() {
        return Token::Option(candidates);
    }
    // argparse's negative-number test, `-\.?\d` at the start, on ASCII digits.
    let digit_at = if arg[1] == b'.' { 2 } else { 1 };
    if arg.get(digit_at).is_some_and(u8::is_ascii_digit) || arg.contains(&b' ') {
        return Token::Value;
    }
    Token::Unknown
}

/// `str(pathlib.PurePosixPath(arg))` on the argument's bytes: repeated separators and `.`
/// components dropped, a trailing separator dropped, exactly two leading slashes kept,
/// `..` left alone, `.` for an empty result.
fn py_path(arg: &[u8]) -> PathBuf {
    let leading_slashes = arg.iter().take_while(|b| **b == b'/').count();
    let root: &[u8] = match leading_slashes {
        0 => b"",
        2 => b"//",
        _ => b"/",
    };
    let parts: Vec<&[u8]> = arg
        .split(|b| *b == b'/')
        .filter(|part| !part.is_empty() && *part != b".")
        .collect();
    let joined = [root, &parts.join(&b'/')].concat();
    if joined.is_empty() {
        PathBuf::from(".")
    } else {
        PathBuf::from(OsStr::from_bytes(&joined))
    }
}

pub fn parse(argv: Vec<OsString>) -> Result<Parsed, String> {
    let argv: Vec<Vec<u8>> = argv.into_iter().map(OsString::into_vec).collect();
    let mut after_separator = false;
    let tokens: Vec<Token> = argv
        .iter()
        .map(|arg| {
            if after_separator {
                Token::Value
            } else if arg == b"--" {
                after_separator = true;
                Token::Separator
            } else {
                classify(arg)
            }
        })
        .collect();

    let shown = |arg: &[u8]| String::from_utf8_lossy(arg).into_owned();
    let mut perl_trace = None;
    let mut rust_output = None;
    let mut discordant_only = false;
    let mut output = b"-".to_vec();
    let mut unrecognized: Vec<String> = Vec::new();
    let mut i = 0;
    while i < tokens.len() {
        let Token::Option(candidates) = &tokens[i] else {
            unrecognized.push(shown(&argv[i]));
            i += 1;
            continue;
        };
        if candidates.len() > 1 {
            return Err(format!("ambiguous option: {}", shown(&argv[i])));
        }
        let first = &candidates[0];
        let opt = first.opt;
        let mut attached = first.attached.clone();
        let mut help = false;
        let mut consumed = 1;
        loop {
            match attached.take() {
                // Text run onto `-h`: each `h` is the flag again, as a cluster of short
                // flags reads; any other letter is an argument nothing takes.
                Some((with_equals, text)) if first.short && !text.is_empty() => {
                    if with_equals || text[0] == b'-' {
                        return Err(format!("ignored explicit argument '{}'", shown(&text)));
                    }
                    help = true;
                    if text[0] != b'h' {
                        unrecognized.push(format!("-{}", shown(&text)));
                        break;
                    }
                    attached = match &text[1..] {
                        [] => None,
                        [b'=', rest @ ..] => Some((true, rest.to_vec())),
                        rest => Some((false, rest.to_vec())),
                    };
                }
                Some((_, text)) => {
                    match opt {
                        Opt::PerlTrace => perl_trace = Some(text),
                        Opt::RustOutput => rust_output = Some(text),
                        Opt::Output => output = text,
                        Opt::Help | Opt::DiscordantOnly => {
                            return Err(format!("ignored explicit argument '{}'", shown(&text)));
                        }
                    }
                    break;
                }
                None => {
                    if opt.takes_value() {
                        let Some(Token::Value) = tokens.get(i + 1) else {
                            return Err(format!(
                                "argument {}: expected one argument",
                                shown(&argv[i])
                            ));
                        };
                        let value = argv[i + 1].clone();
                        match opt {
                            Opt::PerlTrace => perl_trace = Some(value),
                            Opt::RustOutput => rust_output = Some(value),
                            _ => output = value,
                        }
                        consumed = 2;
                    } else if opt == Opt::Help {
                        help = true;
                    } else {
                        discordant_only = true;
                    }
                    break;
                }
            }
        }
        if help {
            return Ok(Parsed::Help);
        }
        i += consumed;
    }

    let missing: Vec<&str> = [
        ("--perl-trace", perl_trace.is_none()),
        ("--rust-output", rust_output.is_none()),
    ]
    .iter()
    .filter(|(_, absent)| *absent)
    .map(|(name, _)| *name)
    .collect();
    let (Some(perl_trace), Some(rust_output)) = (perl_trace, rust_output) else {
        return Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ));
    };
    if !unrecognized.is_empty() {
        return Err(format!(
            "unrecognized arguments: {}",
            unrecognized.join(" ")
        ));
    }
    Ok(Parsed::Run(Args {
        perl_trace: py_path(&perl_trace),
        rust_output: py_path(&rust_output),
        discordant_only,
        output: OsString::from_vec(output),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_strs(args: &[&str]) -> Result<Parsed, String> {
        parse(args.iter().map(OsString::from).collect())
    }

    fn run(args: &[&str]) -> Args {
        match parse_strs(args) {
            Ok(Parsed::Run(args)) => args,
            other => panic!("{args:?}: {other:?}"),
        }
    }

    #[test]
    fn options_take_a_following_value_or_one_after_an_equals_sign() {
        let args = run(&["--perl-trace", "t.tsv", "--rust-output=o.txt"]);
        assert_eq!(args.perl_trace, PathBuf::from("t.tsv"));
        assert_eq!(args.rust_output, PathBuf::from("o.txt"));
        assert!(!args.discordant_only);
        assert_eq!(args.output, OsString::from("-"));
        let args = run(&[
            "--discordant-only",
            "--output",
            "report.txt",
            "--perl-trace=a",
            "--rust-output",
            "b",
            "--perl-trace",
            "c",
        ]);
        assert!(args.discordant_only);
        assert_eq!(args.output, OsString::from("report.txt"));
        assert_eq!(args.perl_trace, PathBuf::from("c"));
    }

    #[test]
    fn an_option_is_named_by_any_prefix_that_fits_it_alone() {
        let args = run(&["--p", "t", "--r=o", "--d", "--o", "x"]);
        assert_eq!(args.perl_trace, PathBuf::from("t"));
        assert_eq!(args.rust_output, PathBuf::from("o"));
        assert!(args.discordant_only);
        assert_eq!(args.output, OsString::from("x"));
        assert_eq!(parse_strs(&["--he"]), Ok(Parsed::Help));
        assert!(
            parse_strs(&["--=x", "--perl-trace", "t", "--rust-output", "o"])
                .unwrap_err()
                .starts_with("ambiguous option")
        );
        for argv in [vec!["-=x", "-h"], vec!["-=", "--help"]] {
            assert!(
                parse_strs(&argv)
                    .unwrap_err()
                    .starts_with("ambiguous option"),
                "{argv:?}"
            );
        }
        assert_eq!(parse_strs(&["-h", "-=x"]), Ok(Parsed::Help));
    }

    #[test]
    fn a_value_that_starts_with_a_dash_is_taken_only_in_three_shapes() {
        let base = ["--perl-trace", "t", "--rust-output", "o", "--output"];
        for value in ["-", "-5", "-.5", "-a b", "--no such"] {
            let mut argv = base.to_vec();
            argv.push(value);
            assert_eq!(run(&argv).output, OsString::from(value), "{value}");
        }
        for value in ["-x", "--other", "-.", "--", "--discordant-only"] {
            let mut argv = base.to_vec();
            argv.push(value);
            assert!(parse_strs(&argv).is_err(), "{value}");
        }
        assert_eq!(
            run(&["--perl-trace", "t", "--rust-output", "o", "--output=-x"]).output,
            OsString::from("-x")
        );
    }

    #[test]
    fn help_is_taken_where_it_stands_in_the_argument_order() {
        for argv in [
            vec!["-h"],
            vec!["--help"],
            vec!["stray", "-h"],
            vec!["--bogus", "--help"],
            vec!["-hh"],
            vec!["-hx"],
            vec!["--discordant-only", "-h", "--perl-trace"],
        ] {
            assert_eq!(parse_strs(&argv), Ok(Parsed::Help), "{argv:?}");
        }
        for argv in [
            vec!["--perl-trace", "-h"],
            vec!["--discordant-only=x", "-h"],
            vec!["-h=x"],
            vec!["-h="],
            vec!["-h-"],
            vec!["-hh=x"],
            vec!["--help=x"],
        ] {
            assert!(parse_strs(&argv).is_err(), "{argv:?}");
        }
    }

    #[test]
    fn missing_required_options_and_stray_arguments_are_errors() {
        assert_eq!(
            parse_strs(&[]).unwrap_err(),
            "the following arguments are required: --perl-trace, --rust-output"
        );
        assert_eq!(
            parse_strs(&["--perl-trace", "t"]).unwrap_err(),
            "the following arguments are required: --rust-output"
        );
        let base = ["--perl-trace", "t", "--rust-output", "o"];
        for stray in [vec!["stray"], vec!["--bogus"], vec!["--"], vec!["--", "x"]] {
            let argv = [base.to_vec(), stray.clone()].concat();
            assert!(
                parse_strs(&argv)
                    .unwrap_err()
                    .starts_with("unrecognized arguments"),
                "{stray:?}"
            );
        }
        assert!(parse_strs(&["--perl-trace"]).is_err());
        assert!(parse_strs(&[
            "--perl-trace",
            "t",
            "--rust-output",
            "o",
            "--discordant-only="
        ])
        .is_err());
    }

    #[test]
    fn input_paths_are_normalised_on_their_bytes_and_the_report_path_is_kept() {
        let args = run(&[
            "--perl-trace",
            "./a//b/",
            "--rust-output=o/.",
            "--output",
            "./r//x/",
        ]);
        assert_eq!(args.perl_trace, PathBuf::from("a/b"));
        assert_eq!(args.rust_output, PathBuf::from("o"));
        assert_eq!(args.output, OsString::from("./r//x/"));
        for (arg, expected) in [
            ("//x", "//x"),
            ("///x", "/x"),
            ("", "."),
            ("./", "."),
            ("x/../y", "x/../y"),
        ] {
            assert_eq!(py_path(arg.as_bytes()), PathBuf::from(expected), "{arg:?}");
        }
        let raw = OsStr::from_bytes(b"./t\xff.tsv/");
        let argv = vec![
            OsString::from("--perl-trace"),
            raw.to_os_string(),
            OsString::from("--rust-output"),
            OsString::from("o"),
        ];
        match parse(argv) {
            Ok(Parsed::Run(args)) => {
                assert_eq!(args.perl_trace.as_os_str(), OsStr::from_bytes(b"t\xff.tsv"));
            }
            other => panic!("{other:?}"),
        }
    }
}
