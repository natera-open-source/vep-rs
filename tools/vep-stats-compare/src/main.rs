// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Section-by-section comparison of two VEP statistics reports.
//!
//! `--expected F --actual F [--tolerance N] [--mask-row LABEL ...] [--report F]`: one
//! Ensembl VEP statistics report (the `<output>_summary.txt` that `--stats_text` writes)
//! against a candidate's, a JSON report with a `verdict`, `stats: <verdict>`, one line
//! per reason and one per row set aside under a documented difference on stdout, exit 0
//! on `green`, 1 on `red` and 2 on an error. Either file may be gzipped (a `.gz` name).
//! The two files are read on two threads.
//!
//! The report: `mode` and `format` (both `stats`), `verdict`, `reasons` (one line per
//! section that differs), `tolerance`, `masked_rows`, `sections` (each file's section
//! titles in order and whether the shared ones agree in order), `compared` (sections and
//! rows compared), `documented` (each row set aside under a documented difference: its
//! `class`, `section`, `row`, `expected` and `actual`), `differences`, a map from section
//! title to the list of its differences (`kind` one of `value`, `count`, `missing`,
//! `unexpected`, `order`, `duplicate`, `section_missing`, `section_unexpected`; `row`,
//! `expected` and `actual` where they apply), then `expected` and `actual` (the paths).

mod cli;
mod input;
mod json;
mod pyfmt;
mod report;
mod terms;
mod unicode;

use crate::pyfmt::{backslashreplace, strerror, surrogateescape};
use std::ffi::OsStr;
use std::fs;
use std::io::{self, Write};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::Path;
use std::process::ExitCode;

/// The component every error line names.
const COMPONENT: &str = "compare_vep_stats";

/// A path or a label of the command line as Python writes it to standard error.
fn shown(arg: &[u8]) -> String {
    backslashreplace(&surrogateescape(arg))
}

/// Writes the report file, its directory made first, the text followed by one newline.
fn write_report(path: &[u8], text: &str) -> io::Result<()> {
    let path = Path::new(OsStr::from_bytes(path));
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, format!("{text}\n"))
}

/// Compares the two files and writes the report; the error is the line's text after the
/// component.
fn run(args: &cli::Args) -> Result<report::Report, String> {
    let (expected, actual) = std::thread::scope(|scope| {
        let expected = scope.spawn(|| input::read_report(&args.expected));
        let actual = input::read_report(&args.actual);
        (expected.join().expect("the reading thread"), actual)
    });
    let unreadable = |path: &[u8], error: String| format!("cannot read {}: {error}", shown(path));
    let expected = expected.map_err(|error| unreadable(&args.expected, error))?;
    let actual = actual.map_err(|error| unreadable(&args.actual, error))?;
    let report = report::compare_reports(expected, actual, &args.tolerance, &args.mask_rows)?;
    if let Some(path) = &args.report {
        let text = report.json(&args.expected, &args.actual).dumps();
        write_report(path, &text)
            .map_err(|error| format!("cannot write {}: {}", shown(path), strerror(&error)))?;
    }
    Ok(report)
}

fn main() -> ExitCode {
    let argv: Vec<Vec<u8>> = std::env::args_os()
        .skip(1)
        .map(OsStringExt::into_vec)
        .collect();
    let args = match cli::parse(&argv) {
        cli::Parsed::Run(args) => args,
        cli::Parsed::Help => {
            println!("{}\n\n{}", cli::USAGE, cli::HELP);
            return ExitCode::SUCCESS;
        }
        cli::Parsed::Usage(message) => {
            eprintln!("{}\nvep-stats-compare: error: {message}", cli::USAGE);
            return ExitCode::from(2);
        }
    };
    match run(&args) {
        Ok(report) => {
            // A reader that closed stdout early is not this run's failure.
            let _ = io::stdout().write_all(report.summary().as_bytes());
            if report.reasons.is_empty() {
                ExitCode::SUCCESS
            } else {
                ExitCode::FAILURE
            }
        }
        Err(message) => {
            eprintln!("ERROR: [{COMPONENT}] {message}");
            ExitCode::from(2)
        }
    }
}
