// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The supplementary columns of fastVEP's tab output as plugin tuples.
//!
//! `<fastvep.tab> <tuples.tsv> [--field-map <map.tsv>] [--fields <name,name>]`: a fastVEP
//! run with `--sa-dir` adds one tab column per loaded source after the built-in columns,
//! and this tool writes those columns as plugin tuples, the form in which a comparator's
//! values are scored against Ensembl VEP's plugin and `--custom` fields: tab-separated,
//! five columns, no header: the row's Location, Allele and Feature, the Ensembl VEP field
//! name, the decoded value. The field map names the Ensembl VEP field of each part of a
//! column. One summary line on stderr counts what was written and what was left out.
//!
//! Exit 0 with the tuple file and the summary line; 1 with one
//! `ERROR: [fastvep_plugin_tuples]` line for a refused run, which leaves an earlier tuple
//! file as it was; 2 for a malformed command line. The rows are written on every core, in
//! input order.

mod cells;
mod cli;
mod columns;
mod map;
mod pyfmt;
mod rows;
mod unicode;

use crate::map::Failure;
use crate::pyfmt::{py_path_str, shown, strerror, uncaught_os_error};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{self, Write};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::os::unix::fs::MetadataExt;
use std::path::Path;
use std::process::ExitCode;

/// The component every refusal names.
const COMPONENT: &str = "fastvep_plugin_tuples";

/// The field map a run without `--field-map` reads: the one of this tool, at the path it
/// was built from.
const FIELD_MAP: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/data/fastvep_field_map.tsv");

/// The error number of a directory opened where a file is read, as Python raises it.
const EISDIR: i32 = 21;

fn path(bytes: &[u8]) -> &Path {
    Path::new(OsStr::from_bytes(bytes))
}

/// Whether two paths name one file, as `pathlib.Path.samefile` decides it: each is read as
/// `pathlib` normalizes it, and a path that cannot be examined names no file.
fn same_file(one: &[u8], other: &[u8]) -> bool {
    let examined = |name: &[u8]| fs::metadata(path(&py_path_str(name)));
    match (examined(one), examined(other)) {
        (Ok(one), Ok(other)) => one.ino() == other.ino() && one.dev() == other.dev(),
        _ => false,
    }
}

/// Opens the input, refusing a directory as Python's `open` does.
fn open_input(name: &[u8]) -> io::Result<File> {
    let file = File::open(path(name))?;
    if file.metadata()?.is_dir() {
        return Err(io::Error::from_raw_os_error(EISDIR));
    }
    Ok(file)
}

/// Writes the tuple file and returns the summary line.
fn run(args: &cli::Args) -> Result<String, Failure> {
    let refused = |message: String| Failure::Refused(message);
    let uncaught = |error: io::Error| Failure::Uncaught(uncaught_os_error(&error));
    let kept = map::kept_fields(args.fields.as_deref())?;
    // An error line names the map as the command line wrote it; the file is the one
    // `pathlib` reads that name as.
    let (map_name, map_file) = match &args.field_map {
        Some(given) => (shown(given), py_path_str(given)),
        None => (FIELD_MAP.to_string(), FIELD_MAP.as_bytes().to_vec()),
    };
    let map_bytes = fs::read(path(&map_file)).map_err(|error| {
        refused(format!(
            "cannot read the field map {map_name}: {}",
            strerror(&error)
        ))
    })?;
    let field_map = map::read_field_map(&map_name, &map_bytes, kept.as_deref())?;
    let input = open_input(&args.input).map_err(|error| {
        refused(format!(
            "cannot read the input {}: {}",
            shown(&args.input),
            strerror(&error)
        ))
    })?;
    let threads = rows::worker_threads(&input).map_err(Failure::Threads)?;
    let mut input = rows::Input::new(input);
    let header = input.header().map_err(uncaught)?;
    let placed =
        columns::place_columns(header.as_deref(), &field_map, kept.is_some()).map_err(refused)?;
    if same_file(&args.tuples, &args.input) {
        return Err(refused(format!(
            "the tuple file {} is the input",
            shown(&args.tuples)
        )));
    }
    let tuples = File::create(path(&args.tuples)).map_err(|error| {
        refused(format!(
            "cannot write the tuple file {}: {}",
            shown(&args.tuples),
            strerror(&error)
        ))
    })?;
    let counts = rows::write_tuples(input, &placed, tuples, threads).map_err(uncaught)?;
    Ok(rows::summary(&placed, &counts))
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
            eprintln!(
                "{}\nvep-fastvep-plugin-tuples: error: {message}",
                cli::USAGE
            );
            return ExitCode::from(2);
        }
    };
    let (line, code) = match run(&args) {
        Ok(summary) => (summary, ExitCode::SUCCESS),
        Err(Failure::Refused(message)) => {
            (format!("ERROR: [{COMPONENT}] {message}"), ExitCode::FAILURE)
        }
        Err(Failure::Uncaught(line)) => (line, ExitCode::FAILURE),
        Err(Failure::Threads(message)) => (
            format!("ERROR: [vep-fastvep-plugin-tuples] {message}"),
            ExitCode::FAILURE,
        ),
    };
    // A reader that closed stderr early is not this run's failure.
    let _ = io::stderr().write_all(format!("{line}\n").as_bytes());
    code
}
