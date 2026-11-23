// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Classes the rows of a discordant-tuple table.
//!
//! A concordance comparison lists the tuples one engine wrote and the other did not. This
//! tool pairs the reference-only rows of such a table with its scored-only rows and names,
//! for every row, why it has no equal on the other side: the end of its Location alone
//! differs, its consequence set differs, the other engine chose another transcript, or
//! the other engine wrote nothing for the variant. It writes the table with that class in
//! a last column and a summary on stdout; both are fixed to the byte, and `README.md`
//! states the contract.
//!
//! The table is read twice and its text is never held: the first read keeps a few numbers
//! a row, the classes are derived from those, and the second read writes each row with
//! its class.

mod classes;
mod cli;
mod csv;
mod intern;
mod pyfmt;
mod report;
mod scan;
mod table;
mod text;
mod threads;

use cli::{Args, Command};
use std::ffi::{OsStr, OsString};
use std::fs::{File, OpenOptions};
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::os::unix::fs::MetadataExt;
use std::path::Path;
use std::process::ExitCode;
use text::Text;

/// The name of the table written beside the input where `-o` names none.
const DEFAULT_OUTPUT: &[u8] = b"discordant_classified.tsv";

fn named(path: &OsStr, detail: impl std::fmt::Display) -> String {
    format!("{}: {detail}", Path::new(path).display())
}

/// The input: a file that is read twice, or the whole of an input that can be read only
/// once.
enum Source {
    File(File),
    Held(Vec<u8>),
}

impl Source {
    fn open(path: &OsStr) -> io::Result<Source> {
        let mut file = File::open(path)?;
        let metadata = file.metadata()?;
        if metadata.is_dir() {
            return Err(io::ErrorKind::IsADirectory.into());
        }
        if metadata.is_file() {
            return Ok(Source::File(file));
        }
        let mut held = Vec::new();
        file.read_to_end(&mut held)?;
        Ok(Source::Held(held))
    }

    /// The text of the input from its start.
    fn text(&mut self) -> io::Result<Text<Box<dyn Read + '_>>> {
        Ok(Text::new(match self {
            Source::File(file) => {
                file.seek(SeekFrom::Start(0))?;
                Box::new(file)
            }
            Source::Held(held) => Box::new(held.as_slice()),
        }))
    }

    /// Holds the whole input where `output` is the file it is read from, so that the
    /// table can be written over it.
    fn hold_if_written_to(&mut self, output: &File) -> io::Result<()> {
        let Source::File(file) = self else {
            return Ok(());
        };
        let (read, written) = (file.metadata()?, output.metadata()?);
        if (read.dev(), read.ino()) == (written.dev(), written.ino()) {
            let mut held = Vec::new();
            file.seek(SeekFrom::Start(0))?;
            file.read_to_end(&mut held)?;
            *self = Source::Held(held);
        }
        Ok(())
    }
}

fn input_failure(input: &OsStr, failure: scan::Failure) -> String {
    match failure {
        scan::Failure::Read(error) => named(input, error),
        scan::Failure::NotUtf8(offset) => named(
            input,
            format_args!("the byte at offset {offset} is not UTF-8"),
        ),
        scan::Failure::Table(detail) => named(input, detail),
    }
}

/// The exit code of a summary that stdout does not take: the code CPython ends with where
/// it cannot flush its stdout at exit.
const STDOUT_FAILED: u8 = 120;

/// A run that ends early: its exit code and the detail of its error line.
struct Refusal {
    code: u8,
    detail: String,
}

impl From<String> for Refusal {
    fn from(detail: String) -> Refusal {
        Refusal { code: 1, detail }
    }
}

fn run(args: &Args) -> Result<(), Refusal> {
    let threads = threads::thread_count()?;
    let input = args.input.as_os_str();
    let mut source = Source::open(input).map_err(|error| named(input, error))?;
    let text = source.text().map_err(|error| named(input, error))?;
    let (header, keys) =
        scan::scan(text, threads).map_err(|failure| input_failure(input, failure))?;
    let classes = classes::classify(&keys, header.class.is_some(), threads);
    drop(keys);

    // The table is made, or emptied where it exists, only once every row has its class: a
    // run that refuses its input leaves the table as it found it.
    let output = match &args.output {
        Some(output) => output.clone(),
        None => OsString::from_vec(pyfmt::sibling(input.as_bytes(), DEFAULT_OUTPUT)),
    };
    let output = output.as_os_str();
    let sink = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(false)
        .open(output)
        .map_err(|error| named(output, error))?;
    source
        .hold_if_written_to(&sink)
        .map_err(|error| named(input, error))?;
    let is_file = sink
        .metadata()
        .map_err(|error| named(output, error))?
        .is_file();
    if is_file {
        sink.set_len(0).map_err(|error| named(output, error))?;
    }
    let text = source.text().map_err(|error| named(input, error))?;
    let tally =
        report::write_table(text, &header, &classes, sink, threads).map_err(
            |failure| match failure {
                report::Failure::Input(failure) => input_failure(input, failure),
                report::Failure::Write(error) => named(output, error),
            },
        )?;

    let (summary, unprintable) = report::summary(&tally, &header);
    let mut stdout = io::stdout().lock();
    stdout
        .write_all(&summary)
        .and_then(|()| stdout.flush())
        .map_err(|error| Refusal {
            code: STDOUT_FAILED,
            detail: format!("stdout: {error}"),
        })?;
    match unprintable {
        None => Ok(()),
        Some(column) => Err(named(
            input,
            format_args!("the summary ends at a row that ends before its {column} column"),
        )
        .into()),
    }
}

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().skip(1).collect();
    let args = match cli::parse(&argv) {
        Ok(Command::Run(args)) => args,
        Ok(Command::Help) => {
            let _ = writeln!(io::stdout(), "{}\n\n{}", cli::USAGE, cli::HELP);
            return ExitCode::SUCCESS;
        }
        Err(error) => {
            let _ = writeln!(
                io::stderr(),
                "{}\nvep-discordant-classes: error: {error}",
                cli::USAGE
            );
            return ExitCode::from(2);
        }
    };
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(Refusal { code, detail }) => {
            let _ = writeln!(io::stderr(), "ERROR: [vep-discordant-classes] {detail}");
            ExitCode::from(code)
        }
    }
}
