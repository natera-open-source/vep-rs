// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Projects a vep-rs Parquet output back to the `--tab` format.
//!
//! A Parquet output is one file, a directory of `part-*.parquet` files, or a directory of
//! `chrom=<c>/` directories holding them. Every file carries the key columns and then the
//! `--tab` columns, each stored under the type the writer's cell contract declares, and
//! [`cells`] inverts that contract cell by cell. The output is the `#`-prefixed header line
//! of the tab columns, then one line per row in file order, the files of a directory in
//! path order. A tab column whose name an earlier column carries, in any letter case, is
//! stored behind `tab_` and printed under its own name, which the file's footer states.
//!
//! `--canonical` prints the canonical form instead: the header line, then the data lines
//! sorted, every DOUBLE cell as its shortest round-trip decimal and every INT32 cell as
//! its integer ([`canonical`]), a cell stored as a string as stored. The input may then
//! also be a `--tab` output, plain or gzip, which `--types-of` reads under the types of the
//! Parquet output it is compared with.
//!
//! Row groups are read across the host's cores and written in order.

mod args;
mod canonical;
mod cells;
mod gunzip;
mod parquet_input;
mod pyfmt;
mod tab_input;

use std::fmt::Display;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;
use std::process::ExitCode;
use std::sync::mpsc::{channel, sync_channel};

use rayon::prelude::*;

use args::{Args, Parsed};
use canonical::{canonical_kind, push_canonical_line, CanonicalKind, Unconvertible};
use parquet_input::{starts_with_parquet_magic, ParquetOutput, Unit};
use tab_input::{data_lines, TabInput};

const COMPONENT: &str = "parquet_to_vep_tab";
const THREADS_VARIABLE: &str = "VEP_PARQUET_TO_TAB_THREADS";
const USAGE: &str = "usage: vep-parquet-to-tab [-h] [--input INPUT_PATH] [--output OUTPUT] \
                     [--canonical] [--types-of PARQUET] [path]";
const HELP: &str = "\
Project a vep-rs Parquet output back to the --tab format: the #-prefixed header line of the
tab columns, then one line per row in file order, the files of a directory in path order.

positional arguments:
  path                  Parquet file or directory (with --canonical, a --tab output too)

options:
  -h, --help            show this help message and exit
  --input, -i INPUT_PATH
                        the same as the positional path
  --output, -o OUTPUT   output path (default stdout)
  --canonical           the canonical form: the header, then the data lines sorted with
                        every DOUBLE cell as its shortest round-trip decimal and every INT32
                        cell as its integer; the input may be a --tab output
  --types-of PARQUET    with --canonical on a --tab output: the Parquet output it is
                        compared with; a column that output stores as strings stays as
                        written

environment:
  VEP_PARQUET_TO_TAB_THREADS
                        the threads at work (default: the host's cores)";

/// Why a run ends without its output.
#[derive(Debug)]
pub enum Failure {
    /// The line after `ERROR: [parquet_to_vep_tab] `; exit code 2.
    Reported(Vec<u8>),
    /// The reader of standard output closed it: nothing more is owed; exit code 0.
    ClosedPipe,
    /// Deflate data no decompressor accepts, in the file named; exit code 1.
    CorruptDeflate(Vec<u8>),
}

impl Failure {
    pub fn message(text: impl Into<String>) -> Self {
        Self::Reported(text.into().into_bytes())
    }

    /// `<name>: <detail>`, the name as given.
    pub fn named(name: &[u8], detail: impl Display) -> Self {
        let mut line = name.to_vec();
        line.extend_from_slice(format!(": {detail}").as_bytes());
        Self::Reported(line)
    }

    pub fn at(path: &Path, detail: impl Display) -> Self {
        Self::named(path.as_os_str().as_bytes(), detail)
    }

    pub fn io(error: &io::Error) -> Self {
        if error.kind() == io::ErrorKind::BrokenPipe {
            Self::ClosedPipe
        } else {
            Self::message(pyfmt::os_error_text(error))
        }
    }

    pub fn corrupt_deflate(name: &[u8]) -> Self {
        Self::CorruptDeflate(name.to_vec())
    }
}

/// The lines of one unit of input: a row group, or a block of a `--tab` file.
struct Lines {
    bytes: Vec<u8>,
    /// Where each line ends in `bytes`; a line holds no line end of its own to split at,
    /// since a stored string may hold one.
    ends: Vec<usize>,
}

impl Lines {
    fn iter(&self) -> impl Iterator<Item = &[u8]> {
        let mut start = 0;
        self.ends.iter().map(move |end| {
            let line = &self.bytes[start..*end];
            start = *end;
            line
        })
    }
}

fn unconvertible(token: Unconvertible) -> Failure {
    let mut line = b"could not convert string to float: ".to_vec();
    pyfmt::push_py_repr(&mut line, &token.0);
    Failure::Reported(line)
}

/// The row groups read at once at most: a single writer keeps pace with no more of them.
const WAVE: usize = 16;

/// The rows of a Parquet output in file order, a wave of row groups at a time: the row
/// groups of a wave are read across the threads while the wave before it is written, and
/// a text that is written goes back to hold a later row group.
fn write_rows(
    output: &ParquetOutput,
    units: &[Unit],
    sink: &mut (dyn Write + Send),
) -> Result<(), Failure> {
    let wave = rayon::current_num_threads().clamp(1, WAVE);
    let (sender, receiver) = sync_channel::<Vec<Vec<u8>>>(0);
    let (spare_sender, spare_receiver) = channel::<Vec<u8>>();
    std::thread::scope(|scope| {
        let writer = scope.spawn(move || -> io::Result<()> {
            for texts in receiver {
                for mut text in texts {
                    sink.write_all(&text)?;
                    text.clear();
                    // The reader may have ended; a text nobody takes is dropped.
                    let _ = spare_sender.send(text);
                }
            }
            Ok(())
        });
        let mut read = Ok(());
        for group in units.chunks(wave) {
            let spares: Vec<Vec<u8>> = group
                .iter()
                .map(|_| spare_receiver.try_recv().unwrap_or_default())
                .collect();
            let texts: Vec<Result<Vec<u8>, Failure>> = group
                .par_iter()
                .zip(spares)
                .map(|(unit, mut text)| {
                    output
                        .for_each_row(*unit, &mut text, |text, _| {
                            text.push(b'\n');
                            Ok(())
                        })
                        .map(|()| text)
                })
                .collect();
            // The row groups before one that cannot be read are still written.
            let mut whole = Vec::with_capacity(texts.len());
            for text in texts {
                match text {
                    Ok(text) => whole.push(text),
                    Err(failure) => {
                        read = Err(failure);
                        break;
                    }
                }
            }
            if sender.send(whole).is_err() || read.is_err() {
                break;
            }
        }
        drop(sender);
        let written = writer.join().expect("the writing thread");
        written.map_err(|error| Failure::io(&error)).and(read)
    })
}

/// The canonical lines of every unit, sorted, written one a line.
fn write_sorted(units: Vec<Lines>, sink: &mut dyn Write) -> Result<(), Failure> {
    let mut lines: Vec<&[u8]> = Vec::with_capacity(units.iter().map(|unit| unit.ends.len()).sum());
    for unit in &units {
        lines.extend(unit.iter());
    }
    lines.par_sort_unstable();
    let mut sink = BufWriter::with_capacity(1 << 20, sink);
    for line in lines {
        sink.write_all(line)
            .and_then(|()| sink.write_all(b"\n"))
            .map_err(|e| Failure::io(&e))?;
    }
    sink.flush().map_err(|e| Failure::io(&e))
}

fn write_header<C: AsRef<[u8]>>(sink: &mut dyn Write, columns: &[C]) -> Result<(), Failure> {
    let mut header = vec![b'#'];
    for (i, column) in columns.iter().enumerate() {
        if i > 0 {
            header.push(b'\t');
        }
        header.extend_from_slice(column.as_ref());
    }
    header.push(b'\n');
    sink.write_all(&header).map_err(|e| Failure::io(&e))
}

/// The canonical kind of each column: the kind of its name, and `Verbatim` where
/// `as_written` takes the column, by its place and its name, for one whose cells are
/// strings.
fn kinds<C: AsRef<[u8]>>(
    columns: &[C],
    as_written: impl Fn(usize, &[u8]) -> bool,
) -> Vec<CanonicalKind> {
    columns
        .iter()
        .enumerate()
        .map(|(place, column)| {
            if as_written(place, column.as_ref()) {
                CanonicalKind::Verbatim
            } else {
                canonical_kind(column.as_ref())
            }
        })
        .collect()
}

/// Whether a Parquet output stands at `path`; a path that holds neither a regular file
/// nor a directory is refused.
fn holds_parquet(path: &Path, display: &[u8]) -> Result<bool, Failure> {
    let Ok(metadata) = fs::metadata(path) else {
        return Err(Failure::named(display, "no such file or directory"));
    };
    if !(metadata.is_file() || metadata.is_dir()) {
        return Err(Failure::named(
            display,
            "neither a regular file nor a directory",
        ));
    }
    if metadata.is_dir() {
        return Ok(true);
    }
    File::open(path)
        .and_then(|mut file| starts_with_parquet_magic(&mut file))
        .map_err(|e| Failure::Reported(pyfmt::os_error_line(&e, display)))
}

/// The names of the tab columns that the Parquet output `named` stores as strings.
fn text_columns_of(named: &[u8]) -> Result<Vec<Vec<u8>>, Failure> {
    let display = pyfmt::py_path_str(named);
    let path = Path::new(std::ffi::OsStr::from_bytes(&display));
    if !holds_parquet(path, &display)? {
        return Err(Failure::named(&display, "not a Parquet file"));
    }
    let output = ParquetOutput::open(path)?;
    Ok(output
        .text_columns()
        .into_iter()
        .map(|name| name.as_bytes().to_vec())
        .collect())
}

/// The first failure in input order, or every unit's lines.
fn in_order(units: Vec<Result<Lines, Failure>>) -> Result<Vec<Lines>, Failure> {
    units.into_iter().collect()
}

fn write_tab(
    path: &Path,
    display: &[u8],
    canonical: bool,
    types_of: Option<&[u8]>,
    sink: &mut (dyn Write + Send),
) -> Result<(), Failure> {
    let parquet = holds_parquet(path, display)?;
    if !parquet && !canonical {
        return Err(Failure::named(display, "not a Parquet file"));
    }
    if !parquet {
        let input = TabInput::read(path, display)?;
        // A `--tab` output does not state its columns' types: the Parquet output it is
        // compared with does, and the columns that output stores as strings stay as written.
        let strings = types_of
            .map(text_columns_of)
            .transpose()?
            .unwrap_or_default();
        write_header(sink, &input.columns)?;
        let kinds = kinds(&input.columns, |_, name| {
            strings.iter().any(|column| column == name)
        });
        let units = input
            .into_blocks()
            .into_par_iter()
            .map(|block| {
                let mut lines = Lines {
                    bytes: Vec::with_capacity(block.len()),
                    ends: Vec::new(),
                };
                for line in data_lines(&block) {
                    push_canonical_line(&mut lines.bytes, &kinds, line).map_err(unconvertible)?;
                    lines.ends.push(lines.bytes.len());
                }
                Ok(lines)
            })
            .collect();
        return write_sorted(in_order(units)?, sink);
    }
    let output = ParquetOutput::open(path)?;
    write_header(sink, &output.columns)?;
    let units = output.units();
    if !canonical {
        return write_rows(&output, &units, sink);
    }
    let units = units
        .par_iter()
        .map(|unit| {
            // A cell printed from a field that stores strings is in its canonical form.
            let prints_text = output.prints_text(*unit);
            let kinds = kinds(&output.columns, |place, _| prints_text[place]);
            let mut ends = Vec::new();
            let mut bytes = Vec::new();
            let mut line = Vec::new();
            output.for_each_row(*unit, &mut bytes, |bytes, start| {
                // The row's tab line is split at its tabs again, so that a stored string
                // holding a tab moves the cells after it as it does in the `--tab` format.
                line.clear();
                line.extend_from_slice(&bytes[start..]);
                bytes.truncate(start);
                push_canonical_line(bytes, &kinds, &line).map_err(unconvertible)?;
                ends.push(bytes.len());
                Ok(())
            })?;
            Ok(Lines { bytes, ends })
        })
        .collect();
    write_sorted(in_order(units)?, sink)
}

fn run(args: &Args) -> Result<(), Failure> {
    let display = pyfmt::py_path_str(&args.path);
    let path = Path::new(std::ffi::OsStr::from_bytes(&display));
    let types_of = args.types_of.as_deref();
    if args.output == b"-" {
        let mut stdout = io::stdout();
        write_tab(path, &display, args.canonical, types_of, &mut stdout)?;
        stdout.flush().map_err(|e| Failure::io(&e))
    } else {
        let mut file = File::create(args::os(&args.output))
            .map_err(|e| Failure::Reported(pyfmt::os_error_line(&e, &args.output)))?;
        write_tab(path, &display, args.canonical, types_of, &mut file)
    }
}

/// Bounds the threads at work when the variable asks for it.
fn bound_threads() -> Result<(), String> {
    let Some(value) = std::env::var_os(THREADS_VARIABLE) else {
        return Ok(());
    };
    let threads = value
        .to_str()
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|threads| *threads >= 1)
        .ok_or_else(|| {
            format!(
                "{THREADS_VARIABLE}={}: not a whole number of 1 or more",
                value.to_string_lossy()
            )
        })?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .map_err(|error| format!("{THREADS_VARIABLE}: {error}"))
}

fn error_line(detail: &[u8]) {
    let mut line = format!("ERROR: [{COMPONENT}] ").into_bytes();
    pyfmt::push_error_stream_text(&mut line, detail);
    line.push(b'\n');
    let _ = io::stderr().write_all(&line);
}

fn main() -> ExitCode {
    let args = match args::parse(std::env::args_os().skip(1).collect()) {
        Parsed::Run(args) => args,
        Parsed::Help => {
            println!("{USAGE}\n\n{HELP}");
            return ExitCode::SUCCESS;
        }
        Parsed::Usage(message) => {
            eprintln!("{USAGE}\nvep-parquet-to-tab: error: {message}");
            return ExitCode::from(2);
        }
    };
    if let Err(message) = bound_threads() {
        eprintln!("{USAGE}\nvep-parquet-to-tab: error: {message}");
        return ExitCode::from(2);
    }
    match run(&args) {
        Ok(()) | Err(Failure::ClosedPipe) => ExitCode::SUCCESS,
        Err(Failure::Reported(detail)) => {
            error_line(&detail);
            ExitCode::from(2)
        }
        Err(Failure::CorruptDeflate(name)) => {
            let mut detail = name;
            detail.extend_from_slice(b": the gzip stream holds deflate data that cannot be read");
            error_line(&detail);
            ExitCode::FAILURE
        }
    }
}
