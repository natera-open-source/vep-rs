// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The command line, read as CPython's `argparse` reads a parser of these seven options with
//! its defaults (`allow_abbrev`, no positionals): an option is named whole or by a prefix
//! that fits it alone, takes its value from the next argument or after `=`, and a later
//! value replaces an earlier one. Errors come in argparse's order: what an option meets as
//! it is read, then the required options that are missing, then the arguments no option
//! took.

use crate::input::KeepOnly;
use crate::pyfmt::{py_int_clamped, py_path_str};
use memchr::memchr;
use std::ffi::{OsStr, OsString};
use std::os::unix::ffi::OsStrExt;
use std::path::PathBuf;

pub const USAGE: &str = "usage: vep-fields-compare [-h] --perl PERL --rust RUST --out OUT [--examples EXAMPLES] [--tmp-dir TMP_DIR] [--keep-only {regulatory_feature_rows}]";
pub const HELP: &str = "\
Field-level agreement between two VEP default-format outputs.

options:
  -h, --help            show this help message and exit
  --perl PERL           reference engine default-format output
  --rust RUST           scored engine default-format output
  --out OUT             JSON report path
  --examples EXAMPLES   mismatch shapes kept per field (default 20)
  --tmp-dir TMP_DIR     directory for the sorted intermediates
  --keep-only {regulatory_feature_rows}
                        compare only RegulatoryFeature and MotifFeature rows on both sides";

#[derive(Debug)]
pub struct Args {
    pub perl: PathBuf,
    pub rust: PathBuf,
    pub out: PathBuf,
    /// Clamped to `i64`: a count at or below zero keeps no mismatch shape and one above the
    /// number of shapes keeps them all, so the clamp changes no byte of the report.
    pub examples: i64,
    pub tmp_dir: Option<PathBuf>,
    pub keep_only: Option<KeepOnly>,
}

/// What a command line asks for.
#[derive(Debug)]
pub enum Command {
    Run(Args),
    Help,
}

/// The long options in the order a usage error lists them.
const LONG_OPTIONS: [&str; 7] = [
    "--help",
    "--perl",
    "--rust",
    "--out",
    "--examples",
    "--tmp-dir",
    "--keep-only",
];

/// How one argument reads before any is consumed.
enum Reading<'a> {
    /// A value for the option before it, or a stray positional.
    Argument,
    /// One option, with the text after `=` when the argument carries one.
    Option(&'static str, Option<&'a [u8]>),
    /// `-h` with more characters behind it.
    ShortHelp(&'a [u8]),
    /// A prefix that more than one option starts with.
    Ambiguous(Vec<&'static str>),
    /// Shaped like an option and naming none.
    Unknown,
}

/// Reads one argument other than the `--` marker. An argument that starts with `-` is an
/// option unless it is `-` alone, starts like a negative number (`-` or `-.` and a digit)
/// or holds a space; such a one is a value.
fn read(arg: &[u8]) -> Reading<'_> {
    if arg.first() != Some(&b'-') || arg.len() == 1 {
        return Reading::Argument;
    }
    let (head, inline) = match memchr(b'=', arg) {
        Some(i) => (&arg[..i], Some(&arg[i + 1..])),
        None => (arg, None),
    };
    if head == b"-h" {
        return Reading::Option("--help", inline);
    }
    if arg[1] == b'-' {
        let exact = LONG_OPTIONS.iter().find(|o| o.as_bytes() == head);
        let matches: Vec<&'static str> = match exact {
            Some(option) => vec![option],
            None => LONG_OPTIONS
                .iter()
                .copied()
                .filter(|o| o.as_bytes().starts_with(head))
                .collect(),
        };
        match matches.as_slice() {
            [] => {}
            [option] => return Reading::Option(option, inline),
            _ => return Reading::Ambiguous(matches),
        }
    } else if arg.starts_with(b"-h") {
        return Reading::ShortHelp(&arg[2..]);
    }
    let negative_number = matches!(arg, [b'-', b'0'..=b'9', ..] | [b'-', b'.', b'0'..=b'9', ..]);
    if negative_number || arg.contains(&b' ') {
        Reading::Argument
    } else {
        Reading::Unknown
    }
}

fn lossy(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// A path option as argparse's `type=Path` reads it: normalised like `str(pathlib.Path(arg))`
/// on its bytes, so `reference.txt/` and `./reference.txt/.` name the file `reference.txt`
/// whether or not the value is UTF-8.
fn path(value: &[u8]) -> PathBuf {
    PathBuf::from(OsStr::from_bytes(&py_path_str(value)))
}

/// The characters behind a `-h`. The option takes no value, so they read as more one-letter
/// options: a further `h` is help again, a letter that names no option is left unread and
/// the help is shown all the same, and a `-` or an `=` there is a value the option refuses.
fn short_help(mut rest: &[u8]) -> Result<(), String> {
    let refused = |value: &[u8]| {
        format!(
            "argument -h/--help: ignored explicit argument '{}'",
            lossy(value)
        )
    };
    loop {
        match rest {
            [b'-', ..] => return Err(refused(rest)),
            [b'h', b'=', value @ ..] => return Err(refused(value)),
            [b'h', more @ ..] if !more.is_empty() => rest = more,
            _ => return Ok(()),
        }
    }
}

pub fn parse(argv: &[OsString]) -> Result<Command, String> {
    let mut perl = None;
    let mut rust = None;
    let mut out = None;
    let mut examples = 20i64;
    let mut tmp_dir = None;
    let mut keep_only = None;
    let mut unread: Vec<String> = Vec::new();
    let mut index = 0;
    while index < argv.len() {
        let arg = argv[index].as_bytes();
        index += 1;
        if arg == b"--" {
            // The marker makes every later argument a positional, and the parser takes none.
            unread.extend(argv[index - 1..].iter().map(|a| lossy(a.as_bytes())));
            break;
        }
        let (name, inline) = match read(arg) {
            Reading::Argument | Reading::Unknown => {
                unread.push(lossy(arg));
                continue;
            }
            Reading::Ambiguous(matches) => {
                return Err(format!(
                    "ambiguous option: {} could match {}",
                    lossy(arg),
                    matches.join(", ")
                ))
            }
            Reading::ShortHelp(rest) => return short_help(rest).map(|()| Command::Help),
            Reading::Option(name, inline) => (name, inline),
        };
        if name == "--help" {
            return match inline {
                Some(value) => Err(format!(
                    "argument -h/--help: ignored explicit argument '{}'",
                    lossy(value)
                )),
                None => Ok(Command::Help),
            };
        }
        let value = match inline {
            Some(value) => value,
            None => match argv.get(index).map(|next| next.as_bytes()) {
                Some(next) if next != b"--" && matches!(read(next), Reading::Argument) => {
                    index += 1;
                    next
                }
                _ => return Err(format!("argument {name}: expected one argument")),
            },
        };
        match name {
            "--perl" => perl = Some(path(value)),
            "--rust" => rust = Some(path(value)),
            "--out" => out = Some(path(value)),
            "--tmp-dir" => tmp_dir = Some(path(value)),
            "--examples" => {
                let text = lossy(value);
                examples = py_int_clamped(&text)
                    .ok_or_else(|| format!("argument --examples: invalid int value: '{text}'"))?;
            }
            "--keep-only" => {
                let text = lossy(value);
                keep_only = Some(KeepOnly::from_name(&text).ok_or_else(|| {
                    format!(
                        "argument --keep-only: invalid choice: '{text}' (choose from '{}')",
                        KeepOnly::RegulatoryFeatureRows.name()
                    )
                })?);
            }
            _ => unreachable!("an option of LONG_OPTIONS"),
        }
    }
    let missing: Vec<&str> = [
        ("--perl", perl.is_none()),
        ("--rust", rust.is_none()),
        ("--out", out.is_none()),
    ]
    .iter()
    .filter(|(_, absent)| *absent)
    .map(|(name, _)| *name)
    .collect();
    match (perl, rust, out) {
        (Some(perl), Some(rust), Some(out)) if unread.is_empty() => Ok(Command::Run(Args {
            perl,
            rust,
            out,
            examples,
            tmp_dir,
            keep_only,
        })),
        (Some(_), Some(_), Some(_)) => Err(format!("unrecognized arguments: {}", unread.join(" "))),
        _ => Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        )),
    }
}

#[cfg(test)]
#[path = "cli_parse_tests.rs"]
mod parse_tests;
