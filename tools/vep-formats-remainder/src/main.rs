// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Front end of `scripts/concordance/compare_vep_formats.py` for the default and tab formats.
//!
//! Reads a reference and a candidate output on every core, groups each side's data lines by
//! their Location, keeps back the groups whose lines are the same text in the same order on
//! both sides, and writes the rest, the remainder, with the counts the comparator adds for
//! what was kept back. No comparison rule lives here: a group is either the same text on
//! both sides or handed to the comparator whole.

mod counters;
mod groups;
mod input;
mod layout;
mod output;
mod scan;

use std::ffi::{OsStr, OsString};
use std::os::unix::ffi::OsStrExt;
use std::path::PathBuf;
use std::process::ExitCode;

const COMPONENT: &str = "vep-formats-remainder";
const THREADS_VARIABLE: &str = "VEP_FORMATS_REMAINDER_THREADS";
const USAGE: &str = "usage: vep-formats-remainder --expected EXPECTED --actual ACTUAL --out-dir OUT_DIR [--force-locations FILE] [--force-text FILE] [--keep-values COLUMN FILE] [--unordered]";
const HELP: &str = "\
Keeps back the record groups two VEP default- or tab-format outputs write alike and writes
the rest for scripts/concordance/compare_vep_formats.py.

options:
  --expected EXPECTED       reference output (gunzipped when its name ends in .gz)
  --actual ACTUAL           candidate output (gunzipped when its name ends in .gz)
  --out-dir OUT_DIR         directory for expected.txt, actual.txt and counters.json
  --force-locations FILE    one Location per line; every group of a listed Location is
                            handed over
  --force-text FILE         one text per line; every group holding a reference line that
                            contains a listed text is handed over
  --keep-values COLUMN FILE one value per line; every group holding a reference line whose
                            COLUMN is not a listed value is handed over
  --unordered               the comparison reads no order: groups are kept back whatever
                            their place in the two files

environment:
  VEP_FORMATS_REMAINDER_THREADS  worker threads (default: the host's available cores)";

pub struct Args {
    pub expected: PathBuf,
    pub actual: PathBuf,
    pub out_dir: PathBuf,
    pub force_locations: Option<PathBuf>,
    pub force_text: Option<PathBuf>,
    pub keep_values: Option<(Vec<u8>, PathBuf)>,
    pub unordered: bool,
}

/// The arguments, or `None` where `-h` or `--help` stands as an option: an option's value
/// spelt that way is the value.
fn parse_args(argv: Vec<OsString>) -> Result<Option<Args>, String> {
    let mut expected = None;
    let mut actual = None;
    let mut out_dir = None;
    let mut force_locations = None;
    let mut force_text = None;
    let mut keep_values = None;
    let mut unordered = false;
    let mut rest = argv.into_iter();
    while let Some(arg) = rest.next() {
        let name = arg.to_string_lossy().into_owned();
        let mut value = |what: &str| {
            rest.next()
                .ok_or_else(|| format!("argument {name}: expected {what}"))
        };
        match name.as_str() {
            "--expected" => expected = Some(PathBuf::from(value("one argument")?)),
            "--actual" => actual = Some(PathBuf::from(value("one argument")?)),
            "--out-dir" => out_dir = Some(PathBuf::from(value("one argument")?)),
            "--force-locations" => force_locations = Some(PathBuf::from(value("one argument")?)),
            "--force-text" => force_text = Some(PathBuf::from(value("one argument")?)),
            "--keep-values" => {
                let column = value("two arguments")?;
                let file = value("two arguments")?;
                keep_values = Some((column.as_bytes().to_vec(), PathBuf::from(file)));
            }
            "--unordered" => unordered = true,
            "-h" | "--help" => return Ok(None),
            _ => return Err(format!("unrecognized arguments: {name}")),
        }
    }
    let required = |value: Option<PathBuf>, name: &str| {
        value.ok_or_else(|| format!("the following arguments are required: {name}"))
    };
    Ok(Some(Args {
        expected: required(expected, "--expected")?,
        actual: required(actual, "--actual")?,
        out_dir: required(out_dir, "--out-dir")?,
        force_locations,
        force_text,
        keep_values,
        unordered,
    }))
}

/// The worker count: `VEP_FORMATS_REMAINDER_THREADS` when set, else the cores the host
/// makes available to the process.
fn thread_count() -> Result<usize, String> {
    match std::env::var_os(THREADS_VARIABLE) {
        Some(value) if !value.is_empty() => parse_threads(&value),
        _ => Ok(std::thread::available_parallelism().map_or(1, usize::from)),
    }
}

fn parse_threads(value: &OsStr) -> Result<usize, String> {
    let text = value.to_string_lossy();
    match text.trim().parse::<usize>() {
        Ok(n) if n >= 1 => Ok(n),
        _ => Err(format!(
            "{THREADS_VARIABLE}={text} is not a positive integer"
        )),
    }
}

fn run(args: &Args) -> Result<(), String> {
    let threads = thread_count()?;
    let force = scan::Force::load(args)?;
    let (expected, actual) = std::thread::scope(|scope| {
        let expected = scope.spawn(|| input::Side::read(&args.expected, threads));
        let actual = scope.spawn(|| input::Side::read(&args.actual, threads));
        (
            expected.join().expect("reference reader thread"),
            actual.join().expect("candidate reader thread"),
        )
    });
    let (expected, actual) = (expected?, actual?);
    let split = groups::split(&expected, &actual, &force, args.unordered, threads)?;
    output::write(&args.out_dir, &expected, &actual, &split, threads)
}

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().skip(1).collect();
    let args = match parse_args(argv) {
        Ok(Some(args)) => args,
        Ok(None) => {
            println!("{USAGE}\n\n{HELP}");
            return ExitCode::SUCCESS;
        }
        Err(message) => {
            eprintln!("{USAGE}\n{COMPONENT}: error: {message}");
            return ExitCode::from(2);
        }
    };
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("ERROR: [{COMPONENT}] {message}");
            ExitCode::from(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_thread_count_is_a_positive_integer() {
        assert_eq!(parse_threads(OsStr::new("16")), Ok(16));
        assert_eq!(parse_threads(OsStr::new(" 2 ")), Ok(2));
        assert!(parse_threads(OsStr::new("0")).is_err());
        assert!(parse_threads(OsStr::new("many")).is_err());
    }

    #[test]
    fn keep_values_takes_a_column_and_a_file() {
        let argv = [
            "--expected",
            "e",
            "--actual",
            "a",
            "--out-dir",
            "o",
            "--keep-values",
            "Feature_type",
            "kept.txt",
            "--unordered",
        ];
        let args = parse_args(argv.iter().map(OsString::from).collect())
            .unwrap()
            .unwrap();
        let (column, file) = args.keep_values.unwrap();
        assert_eq!(column, b"Feature_type");
        assert_eq!(file, PathBuf::from("kept.txt"));
        assert!(args.unordered);
        assert!(parse_args(vec![OsString::from("--expected")]).is_err());
    }

    #[test]
    fn the_help_option_is_an_option_and_never_a_value() {
        let parse = |argv: &[&str]| parse_args(argv.iter().map(OsString::from).collect());
        assert!(parse(&["--help"]).unwrap().is_none());
        assert!(parse(&["--expected", "e", "-h"]).unwrap().is_none());
        let args = parse(&["--expected", "-h", "--actual", "--help", "--out-dir", "-h"])
            .unwrap()
            .unwrap();
        assert_eq!(args.expected, PathBuf::from("-h"));
        assert_eq!(args.actual, PathBuf::from("--help"));
        assert_eq!(args.out_dir, PathBuf::from("-h"));
    }
}
