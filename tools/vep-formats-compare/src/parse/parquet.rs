// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reading a Parquet output as tab text.
//!
//! This module writes the text `tools/vep-parquet-to-tab` writes on its standard output for
//! the same output: a tab column prints the cells of its own field, whatever its name. A Parquet
//! output is one file, a directory of part files or a directory of `chrom=<c>/` directories
//! holding them: the key columns lead every file and the columns after them are the `--tab`
//! columns, a column whose name an earlier column carries in any letter case stored behind
//! `tab_` and printed under its own name, which the file's footer states, each under the
//! type of the writer's cell contract, which is inverted cell by cell:
//!
//! - a NULL cell prints `-`, a string as stored, an integer in decimal, a double as the
//!   shortest decimal that reads back to the same value, in plain notation;
//! - a boolean prints its column's spelling, `1` and `0` inside a list;
//! - a list prints its elements `,`-joined, a NULL element as `-`;
//! - a struct of positions prints `start` where `start` and `end` agree, `start-end`
//!   otherwise with `?` for a NULL member, `-` where both are NULL, then `/total`.
//!
//! The text is the `#` line of the tab columns, then one line a row, the files in the
//! order of their paths. The canonical form has the rows sorted, every DOUBLE cell as its
//! shortest decimal and every INT32 cell as its integer, and is also read off a `--tab`
//! output, whose `##` lines are dropped. A column an output stores as strings holds what a
//! plugin or a `--custom` source wrote, whatever its name: its cells are as stored, and as
//! written in the `--tab` output a comparison pairs with it.
//!
//! `CanonicalText` is the canonical form with its rows in the output's own order, written
//! a few thousand rows at a time and never held whole: the text the stream of record
//! groups reads of each output of a Parquet comparison.
//!
//! A failure is reported in the words of the adapter's standard error. Its reader raises
//! on a file that is no Parquet file with texts of its own; the ones a file's size and
//! last eight bytes decide are written here, and any other refusal of a file carries the
//! text of this reader.

use crate::argparse::{backslashreplace, path_text, repr as repr_argument};
use crate::cache;
use crate::error::{Error, Result};
use crate::input::{gunzip, universal_newlines};
use crate::pyfmt;
use crate::report::strerror;
use crate::stream::gzip::{self, Gunzip};
use crate::stream::Text;
use parquet::basic::{ConvertedType, LogicalType, Repetition, Type as PhysicalType};
use parquet::column::reader::ColumnReader;
use parquet::data_type::ByteArray;
use parquet::file::metadata::FileMetaData;
use parquet::file::reader::FileReader;
use parquet::file::serialized_reader::SerializedFileReader;
use parquet::schema::types::{SchemaDescriptor, Type};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ffi::OsString;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, Read};
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::os::unix::fs::FileExt;
use std::panic;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;

/// The component of the adapter's error line.
const COMPONENT: &str = "parquet_to_vep_tab";

/// The name and the usage line the adapter's option parser opens a refusal with.
const ADAPTER_PROGRAM: &str = "parquet_to_vep_tab.py";
const ADAPTER_USAGE: &str =
    "usage: parquet_to_vep_tab.py [-h] [--input INPUT_PATH] [--output OUTPUT] [--canonical] [path]";

const KEY_COLUMNS: &[&str] = &["contig_index", "chrom", "pos", "end", "ref", "alt"];

/// What the writer stores a tab column named as a key column behind.
const TAB_COLUMN_PREFIX: &str = "tab_";

/// Whether `name` is a key column's name behind any number of `TAB_COLUMN_PREFIX`, none
/// included.
fn names_a_key_column(name: &str) -> bool {
    let mut rest = name;
    while let Some(inner) = rest.strip_prefix(TAB_COLUMN_PREFIX) {
        rest = inner;
    }
    KEY_COLUMNS.contains(&rest)
}

/// The tab column a field after the key columns holds, by the name it is stored under, as
/// `tab_column_name` of the writer's module reads it. The writer stores a tab column named
/// as a key column, behind any number of the prefix, behind one more (`ref` as `tab_ref`,
/// `tab_ref` as `tab_tab_ref`), so that no two columns of a file share a name; every other
/// column is stored under its own where no earlier column carries that name in other
/// letters (`tab_columns`).
fn tab_column_name(stored: &str) -> &str {
    match stored.strip_prefix(TAB_COLUMN_PREFIX) {
        Some(tab_column) if names_a_key_column(tab_column) => tab_column,
        _ => stored,
    }
}

/// The spelling of a boolean inside a list and of a BOOLEAN column of no listed spelling:
/// true, then false.
const FLAG_SPELLING: (&str, &str) = ("1", "0");

/// The cache columns behind the population frequency fields; the frequency column of each
/// is `<column>_AF`, `AF` for `AF`.
const FREQUENCY_POPULATIONS: &[&str] = &[
    "AF",
    "AFR",
    "AMR",
    "ASN",
    "EAS",
    "EUR",
    "SAS",
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_MID",
    "gnomADe_NFE",
    "gnomADe_REMAINING",
    "gnomADe_SAS",
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_REMAINING",
    "gnomADg_SAS",
];
const DOUBLE_COLUMNS: &[&str] = &["MAX_AF", "MOTIF_SCORE_CHANGE"];
const INT32_COLUMNS: &[&str] = &[
    "STRAND",
    "DISTANCE",
    "ALLELE_NUM",
    "TSL",
    "HGVS_OFFSET",
    "OverlapBP",
    "MOTIF_POS",
];

const NULL_TEXT: &str = "-";
const PARQUET_MAGIC: &[u8] = b"PAR1";
const ENCRYPTED_MAGIC: &[u8] = b"PARE";
const GZIP_MAGIC: &[u8] = &[0x1f, 0x8b];

/// The magic bytes that end a Parquet file and the length of its footer before them.
const FOOTER_SIZE: u64 = 8;

/// The rows the adapter reads of a file at a time, a row group after another: it turns
/// every column of that many rows into values before it prints the first of them, so a
/// string that is not UTF-8 anywhere in them is met before a cell of any of their rows.
const BATCH_ROWS: usize = 65536;

/// The data lines a worker puts in the canonical form at a time.
const TAB_LINES: usize = 1 << 16;

/// The rows of a Parquet file a `CanonicalText` prints at a time.
const TEXT_ROWS: usize = 1 << 13;

/// The bytes of a `--tab` output a `CanonicalText` reads at a time; a longer line is read
/// whole.
const TEXT_BYTES: usize = 1 << 20;

/// What the adapter writes to its standard error before it ends with a failure.
struct Failure(String);

/// The failure of the adapter's one error line.
fn failure(message: impl std::fmt::Display) -> Failure {
    Failure(format!("ERROR: [{COMPONENT}] {message}"))
}

/// The failure of an operation of the system on `path`, in the words of `str(OSError)`.
fn system_failure(error: &std::io::Error, path: &Path) -> Failure {
    let named = repr_argument(&path_text(path));
    match error.raw_os_error() {
        Some(code) => failure(format!("[Errno {code}] {}: {named}", strerror(error))),
        None => failure(format!("{}: {named}", strerror(error))),
    }
}

/// A path as a line of the adapter spells it.
fn shown(path: &Path) -> String {
    backslashreplace(&path_text(path))
}

/// The `--tab` spelling of a BOOLEAN column's values: true, then false.
fn boolean_spelling(column: &str) -> (&'static str, &'static str) {
    match column {
        "CANONICAL" => ("YES", NULL_TEXT),
        "PICK" | "GENCODE_PRIMARY" => ("1", NULL_TEXT),
        "HIGH_INF_POS" => ("Y", "N"),
        _ => FLAG_SPELLING,
    }
}

/// What the canonical form does to a cell, by its column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Canonical {
    Verbatim,
    Double,
    DoubleList,
    Int32,
}

fn canonical_of(column: &str) -> Canonical {
    let frequency = column == "AF"
        || column.strip_suffix("_AF").is_some_and(|population| {
            population != "AF" && FREQUENCY_POPULATIONS.contains(&population)
        });
    if DOUBLE_COLUMNS.contains(&column) {
        Canonical::Double
    } else if frequency {
        Canonical::DoubleList
    } else if INT32_COLUMNS.contains(&column) {
        Canonical::Int32
    } else {
        Canonical::Verbatim
    }
}

/// The shortest decimal that reads back to `value`, in plain notation: `0.095`,
/// `0.000098` for 9.8e-05, `1` for 1.0, `-0` for -0.0. Of two shortest decimals equally
/// close to the value this is the one of larger magnitude.
fn write_shortest_decimal(out: &mut String, value: f64) {
    let _ = write!(out, "{value}");
}

/// Whether `c` reads as the letter `i` without regard to case, as the interpreter's
/// pattern matcher has it: the dotted capital and the dotless small letter do.
fn is_letter_i(c: char) -> bool {
    matches!(c, 'i' | 'I' | '\u{130}' | '\u{131}')
}

/// Whether `text` is `name` without regard to case; `name` is in lower case.
fn is_named(text: &str, name: &str) -> bool {
    let mut chars = text.chars();
    name.chars().all(|letter| {
        chars.next().is_some_and(|c| match letter {
            'i' => is_letter_i(c),
            _ => c.to_ascii_lowercase() == letter,
        })
    }) && chars.next().is_none()
}

/// Whether a token is a decimal number as the writer's cell contract parses it: an
/// optional sign, digits with an optional fraction or a fraction alone, an optional
/// exponent; or an infinity or NaN. A digit is a decimal digit of any script, and one
/// line end may follow the number.
fn follows_the_number_grammar(token: &str) -> bool {
    let body = token.strip_suffix('\n').unwrap_or(token);
    let body = body.strip_prefix(['+', '-']).unwrap_or(body);
    if ["inf", "infinity", "nan"]
        .iter()
        .any(|name| is_named(body, name))
    {
        return true;
    }
    let digits = |text: &str| -> usize {
        text.chars()
            .take_while(|c| pyfmt::char_is_decimal(*c))
            .map(char::len_utf8)
            .sum()
    };
    let whole = digits(body);
    let mut rest = &body[whole..];
    if whole > 0 {
        if let Some(fraction) = rest.strip_prefix('.') {
            rest = &fraction[digits(fraction)..];
        }
    } else {
        let Some(fraction) = rest.strip_prefix('.') else {
            return false;
        };
        let places = digits(fraction);
        if places == 0 {
            return false;
        }
        rest = &fraction[places..];
    }
    if rest.is_empty() {
        return true;
    }
    let Some(exponent) = rest.strip_prefix(['e', 'E']) else {
        return false;
    };
    let exponent = exponent.strip_prefix(['+', '-']).unwrap_or(exponent);
    let places = digits(exponent);
    places > 0 && places == exponent.len()
}

/// Appends a DOUBLE cell's token as its shortest decimal; a token outside the grammar
/// stays as written. The error is the interpreter's refusal of a token the grammar takes
/// and `float()` does not.
fn write_canonical_number(out: &mut String, token: &str) -> std::result::Result<(), String> {
    if !follows_the_number_grammar(token) {
        out.push_str(token);
        return Ok(());
    }
    match pyfmt::float_of(token) {
        Some(value) => {
            write_shortest_decimal(out, value);
            Ok(())
        }
        None => Err(pyfmt::float_error(token)),
    }
}

/// Appends an INT32 cell's token as its integer (`007` is `7`), as an INT32 column prints
/// the value. A token that is no INT32 stays as written, as the string column a plugin or a
/// `--custom` source declares under the name stores it; an INT32 column stores NULL for it,
/// which prints `-` and so differs.
fn write_canonical_int32(out: &mut String, token: &str) {
    let (negative, digits) = match token.as_bytes() {
        [b'-', digits @ ..] => (true, digits),
        [b'+', digits @ ..] => (false, digits),
        digits => (false, digits),
    };
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        out.push_str(token);
        return;
    }
    let significant = digits.iter().take_while(|digit| **digit == b'0').count();
    let digits = &digits[significant.min(digits.len() - 1)..];
    // Ten digits hold every INT32.
    if digits.len() > 10 {
        out.push_str(token);
        return;
    }
    let magnitude = digits
        .iter()
        .fold(0i64, |value, digit| value * 10 + i64::from(digit - b'0'));
    let value = if negative { -magnitude } else { magnitude };
    if i32::try_from(value).is_ok() {
        let _ = write!(out, "{value}");
    } else {
        out.push_str(token);
    }
}

/// Appends a `--tab` cell in canonical form: a DOUBLE or `LIST<DOUBLE>` cell with each
/// number as its shortest decimal, an INT32 cell as its integer, every other cell unchanged.
fn write_canonical_cell(
    out: &mut String,
    canonical: Canonical,
    text: &str,
) -> std::result::Result<(), String> {
    if text == NULL_TEXT {
        out.push_str(text);
        return Ok(());
    }
    match canonical {
        Canonical::Verbatim => out.push_str(text),
        Canonical::Double => write_canonical_number(out, text)?,
        Canonical::DoubleList => {
            for (at, token) in text.split(',').enumerate() {
                if at > 0 {
                    out.push(',');
                }
                write_canonical_number(out, token)?;
            }
        }
        Canonical::Int32 => write_canonical_int32(out, text),
    }
    Ok(())
}

/// Appends a data line in canonical form; a cell past the last header column stays as
/// written.
fn write_canonical_line(
    out: &mut String,
    columns: &[Canonical],
    line: &str,
) -> std::result::Result<(), String> {
    for (at, cell) in line.split('\t').enumerate() {
        if at > 0 {
            out.push('\t');
        }
        let canonical = columns.get(at).copied().unwrap_or(Canonical::Verbatim);
        write_canonical_cell(out, canonical, cell)?;
    }
    Ok(())
}

/// Data lines one after another, each followed by a line end. A line may hold line ends of
/// its own, so the lines are told apart by where they start.
#[derive(Default)]
struct Lines {
    text: String,
    /// The start of every line, then the length of the text.
    starts: Vec<usize>,
    /// The lines in sorted order, where the canonical form is written.
    order: Vec<u32>,
}

impl Lines {
    fn new() -> Lines {
        Lines {
            text: String::new(),
            starts: vec![0],
            order: Vec::new(),
        }
    }

    /// Leaves no line, the text's room kept.
    fn clear(&mut self) {
        self.text.clear();
        self.starts.clear();
        self.starts.push(0);
        self.order.clear();
    }

    /// Ends the line written since the last one.
    fn end_line(&mut self) {
        self.text.push('\n');
        self.starts.push(self.text.len());
    }

    fn len(&self) -> usize {
        self.starts.len() - 1
    }

    /// A line without its line end.
    fn line(&self, at: usize) -> &str {
        &self.text[self.starts[at]..self.starts[at + 1] - 1]
    }

    fn sort(&mut self) {
        let mut order: Vec<u32> = (0..self.len() as u32).collect();
        order.sort_unstable_by(|a, b| self.line(*a as usize).cmp(self.line(*b as usize)));
        self.order = order;
    }
}

/// The text of the adapter's standard output: the header line, then the lines of every
/// part, sorted across the parts in the canonical form.
fn joined(columns: &[String], parts: Vec<Lines>, canonical: bool) -> String {
    let size: usize = parts.iter().map(|part| part.text.len()).sum();
    let header: usize = columns.iter().map(|column| column.len() + 1).sum();
    let mut out = String::with_capacity(size + header + 1);
    out.push('#');
    out.push_str(&columns.join("\t"));
    out.push('\n');
    if !canonical {
        for part in parts {
            out.push_str(&part.text);
        }
        return out;
    }
    // Every part is sorted: the least line of all is the least of the parts' next lines.
    let mut next: BinaryHeap<Reverse<(&str, usize, usize)>> = parts
        .iter()
        .enumerate()
        .filter(|(_, part)| part.len() > 0)
        .map(|(at, part)| Reverse((part.line(part.order[0] as usize), at, 0)))
        .collect();
    while let Some(Reverse((line, at, place))) = next.pop() {
        out.push_str(line);
        out.push('\n');
        let part = &parts[at];
        if let Some(following) = part.order.get(place + 1) {
            next.push(Reverse((part.line(*following as usize), at, place + 1)));
        }
    }
    out
}

/// `work(state, task)` of every task from 0 to `tasks`, on the workers a reader of many
/// shards runs on, each with a state of its own; the results in the tasks' order.
fn on_workers<T, S>(
    tasks: usize,
    state: impl Fn() -> S + Sync,
    work: impl Fn(&mut S, usize) -> T + Sync,
) -> Vec<T>
where
    T: Send,
{
    let workers = cache::threads().clamp(1, tasks.max(1));
    if workers == 1 {
        let mut state = state();
        return (0..tasks).map(|task| work(&mut state, task)).collect();
    }
    let next = AtomicUsize::new(0);
    let mut results: Vec<Option<T>> = Vec::new();
    results.resize_with(tasks, || None);
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    let mut state = state();
                    let mut done = Vec::new();
                    loop {
                        let task = next.fetch_add(1, Ordering::Relaxed);
                        if task >= tasks {
                            break;
                        }
                        done.push((task, work(&mut state, task)));
                    }
                    done
                })
            })
            .collect();
        for handle in handles {
            match handle.join() {
                Ok(done) => {
                    for (task, result) in done {
                        results[task] = Some(result);
                    }
                }
                Err(payload) => panic::resume_unwind(payload),
            }
        }
    });
    results
        .into_iter()
        .map(|result| result.expect("every task was taken by a worker"))
        .collect()
}

/// The canonical form of a `--tab` output: its columns and its data lines, the `##` lines
/// dropped. A line ends at `\n` alone, a `\r` before it dropped, so a cell holding another
/// separator character stays in its row. A column named in `as_written` is a string column
/// of the Parquet output the tab output is compared with, and its cells stay as written.
fn tab_output(path: &Path, as_written: &[String]) -> std::result::Result<String, Failure> {
    let bytes = std::fs::read(path).map_err(|error| system_failure(&error, path))?;
    let bytes = if bytes.starts_with(GZIP_MAGIC) {
        gunzip(&bytes).map_err(|refusal| {
            // The adapter catches what `gzip` raises and not what its inflater raises.
            if refusal.starts_with("Error ") {
                Failure(format!("zlib.error: {refusal}"))
            } else {
                failure(refusal)
            }
        })?
    } else {
        bytes
    };
    let text = String::from_utf8(bytes).map_err(|error| {
        failure(
            pyfmt::decode_utf8(error.as_bytes())
                .err()
                .unwrap_or_default(),
        )
    })?;
    let mut header: Option<&str> = None;
    let mut data: Vec<&str> = Vec::new();
    for line in text.split('\n') {
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            continue;
        }
        if let Some(columns) = line.strip_prefix('#') {
            if header.is_some() {
                return Err(failure(format!("{}: two header lines", shown(path))));
            }
            header = Some(columns);
        } else if !line.is_empty() {
            data.push(line);
        }
    }
    let Some(header) = header else {
        return Err(failure(format!(
            "{}: no # header line; a tab output needs its header",
            shown(path)
        )));
    };
    let columns: Vec<String> = header.split('\t').map(str::to_string).collect();
    let canonical: Vec<Canonical> = columns
        .iter()
        .map(|column| {
            if as_written.contains(column) {
                Canonical::Verbatim
            } else {
                canonical_of(column)
            }
        })
        .collect();
    let parts = on_workers(
        data.len().div_ceil(TAB_LINES),
        || (),
        |(), part| {
            let mut lines = Lines::new();
            let from = part * TAB_LINES;
            for line in &data[from..data.len().min(from + TAB_LINES)] {
                write_canonical_line(&mut lines.text, &canonical, line).map_err(failure)?;
                lines.end_line();
            }
            lines.sort();
            Ok(lines)
        },
    );
    let parts = parts
        .into_iter()
        .collect::<std::result::Result<Vec<Lines>, Failure>>()?;
    Ok(joined(&columns, parts, true))
}

/// The value of one leaf column as the adapter's reader hands it back.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LeafKind {
    Bool,
    Int32,
    /// An INT32 that holds an unsigned integer.
    UInt32,
    Int64,
    /// An INT64 that holds an unsigned integer.
    UInt64,
    Float,
    Double,
    Text,
}

/// One leaf column of a file.
#[derive(Debug, Clone, Copy)]
struct Leaf {
    index: usize,
    kind: LeafKind,
    /// The definition level of a value that is not NULL.
    max_def: i16,
}

/// How the leaf columns of a tab column make its cell.
#[derive(Debug, Clone, Copy)]
enum Shape {
    /// A value a row.
    Scalar(Leaf),
    /// A list a row: NULL below the level `present`, empty at it.
    List { leaf: Leaf, present: i16 },
    /// A struct a row, of which `start`, `end` and `total` are read: NULL below `present`.
    Positions {
        members: [Option<Leaf>; 3],
        present: i16,
    },
    /// A struct with none of the three members, which prints NULL in every row.
    NoPositions,
}

/// One tab column of a file.
#[derive(Debug, Clone)]
struct Column {
    shape: Shape,
    spelling: (&'static str, &'static str),
    /// What the canonical form does to the cell as it is printed; `Verbatim` for a column
    /// that stores strings and where the printed text is in that form already.
    canonical: Canonical,
}

/// The value kind of a leaf column, or what it holds that this reader does not read.
fn leaf_kind(schema: &SchemaDescriptor, index: usize) -> std::result::Result<LeafKind, String> {
    let column = schema.column(index);
    let logical = column.logical_type_ref();
    let converted = column.converted_type();
    let unsigned = match logical {
        Some(LogicalType::Integer(integer)) => Some(!integer.is_signed),
        // A column of NULL alone is stored as an INT32 of unknown type.
        Some(LogicalType::Unknown) => Some(false),
        Some(_) => None,
        None => match converted {
            ConvertedType::NONE
            | ConvertedType::INT_8
            | ConvertedType::INT_16
            | ConvertedType::INT_32
            | ConvertedType::INT_64 => Some(false),
            ConvertedType::UINT_8
            | ConvertedType::UINT_16
            | ConvertedType::UINT_32
            | ConvertedType::UINT_64 => Some(true),
            _ => None,
        },
    };
    let annotated = logical.is_some() || converted != ConvertedType::NONE;
    let text = matches!(logical, Some(LogicalType::String))
        || (logical.is_none() && converted == ConvertedType::UTF8);
    match (column.physical_type(), unsigned) {
        (PhysicalType::BOOLEAN, _) if !annotated => Ok(LeafKind::Bool),
        (PhysicalType::INT32, Some(false)) => Ok(LeafKind::Int32),
        (PhysicalType::INT32, Some(true)) => Ok(LeafKind::UInt32),
        (PhysicalType::INT64, Some(false)) => Ok(LeafKind::Int64),
        (PhysicalType::INT64, Some(true)) => Ok(LeafKind::UInt64),
        (PhysicalType::FLOAT, _) if !annotated => Ok(LeafKind::Float),
        (PhysicalType::DOUBLE, _) if !annotated => Ok(LeafKind::Double),
        (PhysicalType::BYTE_ARRAY, _) if text => Ok(LeafKind::Text),
        (physical, _) => Err(match logical {
            Some(logical) => format!("{physical} values annotated {}", annotation_name(logical)),
            None if converted != ConvertedType::NONE => {
                format!("{physical} values annotated {converted}")
            }
            None => format!("{physical} values"),
        }),
    }
}

/// The name of an annotation, without its parameters.
fn annotation_name(logical: &LogicalType) -> String {
    let spelt = format!("{logical:?}");
    spelt
        .split(['(', ' ', '{'])
        .next()
        .unwrap_or("")
        .to_string()
}

fn leaf(schema: &SchemaDescriptor, index: usize) -> std::result::Result<Leaf, String> {
    Ok(Leaf {
        index,
        kind: leaf_kind(schema, index)?,
        max_def: schema.column(index).max_def_level(),
    })
}

fn repetition(field: &Type) -> Repetition {
    let info = field.get_basic_info();
    if info.has_repetition() {
        info.repetition()
    } else {
        Repetition::REQUIRED
    }
}

fn leaves_below(field: &Type) -> usize {
    if field.is_primitive() {
        1
    } else {
        field
            .get_fields()
            .iter()
            .map(|child| leaves_below(child))
            .sum()
    }
}

/// How the top-level field `field`, whose first leaf column is `first`, makes a cell; or
/// what it holds that this reader does not read.
fn shape_of(
    schema: &SchemaDescriptor,
    field: &Type,
    first: usize,
) -> std::result::Result<Shape, String> {
    if field.is_primitive() {
        let leaf = leaf(schema, first)?;
        return Ok(match repetition(field) {
            Repetition::REPEATED => Shape::List { leaf, present: 0 },
            _ => Shape::Scalar(leaf),
        });
    }
    let present = match repetition(field) {
        Repetition::REPEATED => return Err("a repeated group".to_string()),
        Repetition::OPTIONAL => 1,
        Repetition::REQUIRED => 0,
    };
    let info = field.get_basic_info();
    let annotation = (info.logical_type_ref(), info.converted_type());
    let fields = field.get_fields();
    match annotation {
        (Some(LogicalType::List), _) | (None, ConvertedType::LIST) => {
            let [repeated] = fields else {
                return Err("a list group without its one repeated field".to_string());
            };
            if repetition(repeated) != Repetition::REPEATED {
                return Err("a list group without its one repeated field".to_string());
            }
            // The repeated field is the element itself where it is a value, a group of
            // several fields or a group named as the older writers name an element.
            let element = if repeated.is_primitive() {
                repeated
            } else {
                let named_as_an_element = repeated.name() == "array"
                    || repeated.name() == format!("{}_tuple", field.name());
                match repeated.get_fields() {
                    [element] if !named_as_an_element => element,
                    _ => return Err("a list of structs".to_string()),
                }
            };
            if !element.is_primitive()
                || (!std::ptr::eq(element, repeated) && repetition(element) == Repetition::REPEATED)
            {
                return Err("a list of lists or of structs".to_string());
            }
            Ok(Shape::List {
                leaf: leaf(schema, first)?,
                present,
            })
        }
        (None, ConvertedType::NONE) => {
            let mut members = [None; 3];
            let mut index = first;
            for child in fields {
                let member = ["start", "end", "total"]
                    .iter()
                    .position(|name| child.name() == *name);
                if let Some(member) = member {
                    if !child.is_primitive() || repetition(child) == Repetition::REPEATED {
                        return Err(format!("a struct whose {} is no value", child.name()));
                    }
                    members[member] = Some(leaf(schema, index)?);
                }
                index += leaves_below(child);
            }
            if members.iter().all(Option::is_none) {
                Ok(Shape::NoPositions)
            } else {
                Ok(Shape::Positions { members, present })
            }
        }
        (Some(logical), _) => Err(format!("a group annotated {}", annotation_name(logical))),
        (None, converted) => Err(format!("a group annotated {converted}")),
    }
}

/// Whether the cell a column of this shape prints is in the canonical form of `canonical`
/// as it stands: a double prints as its shortest decimal and an INT32 as its integer.
fn prints_canonical(canonical: Canonical, shape: &Shape) -> bool {
    match (canonical, shape) {
        (Canonical::Verbatim, _) => true,
        (
            Canonical::Double | Canonical::DoubleList,
            Shape::Scalar(leaf) | Shape::List { leaf, .. },
        ) => {
            matches!(leaf.kind, LeafKind::Float | LeafKind::Double)
        }
        (Canonical::Int32, Shape::Scalar(leaf)) => leaf.kind == LeafKind::Int32,
        _ => false,
    }
}

/// Whether a column of this shape stores strings, one a row or the elements of a list.
fn stores_text(shape: &Shape) -> bool {
    matches!(
        shape,
        Shape::Scalar(leaf) | Shape::List { leaf, .. } if leaf.kind == LeafKind::Text
    )
}

/// One Parquet file of an output.
struct PartFile {
    path: PathBuf,
    /// The count of the key columns that lead the file's fields.
    key_columns: usize,
    /// The tab columns: every top-level field after the leading key columns.
    names: Vec<String>,
    columns: Vec<Column>,
    /// The rows of every row group.
    row_groups: Vec<usize>,
}

/// The Parquet files of an output: the file itself, or every `.parquet` file under the
/// directory (a `chrom=<c>/` directory included) in sorted path order. A directory
/// reached through a link is not searched, and one that cannot be listed holds nothing.
fn parquet_files(path: &Path, is_dir: bool) -> std::result::Result<Vec<PathBuf>, Failure> {
    if !is_dir {
        return Ok(vec![path.to_path_buf()]);
    }
    fn search(directory: &Path, above: &[Vec<u8>], found: &mut Vec<Vec<Vec<u8>>>) {
        let Ok(entries) = std::fs::read_dir(directory) else {
            return;
        };
        for entry in entries.flatten() {
            let name = entry.file_name().as_bytes().to_vec();
            let mut parts = above.to_vec();
            parts.push(name.clone());
            let below = entry.path();
            if name.ends_with(b".parquet") && below.is_file() {
                found.push(parts.clone());
            }
            if entry.file_type().is_ok_and(|kind| kind.is_dir()) {
                search(&below, &parts, found);
            }
        }
    }
    let mut found = Vec::new();
    search(path, &[], &mut found);
    if found.is_empty() {
        return Err(failure(format!(
            "{}: no .parquet file under the directory",
            shown(path)
        )));
    }
    // Paths compare part by part.
    found.sort();
    let base = path.as_os_str().as_bytes();
    Ok(found
        .into_iter()
        .map(|parts| {
            let mut raw = if base == b"." {
                Vec::new()
            } else {
                base.to_vec()
            };
            for part in parts {
                if !raw.is_empty() && !raw.ends_with(b"/") {
                    raw.push(b'/');
                }
                raw.extend_from_slice(&part);
            }
            PathBuf::from(OsString::from_vec(raw))
        })
        .collect())
}

/// Opens a Parquet file. What the size and the last eight bytes of a file that is no
/// Parquet file decide is refused in the words of the adapter's reader.
fn open(path: &Path) -> std::result::Result<SerializedFileReader<File>, Failure> {
    let file = File::open(path).map_err(|error| {
        let detail = strerror(&error);
        let code = error.raw_os_error().unwrap_or(0);
        failure(format!(
            "[Errno {code}] Failed to open local file '{}'. Detail: [errno {code}] {detail}",
            shown(path)
        ))
    })?;
    let size = file
        .metadata()
        .map_err(|error| system_failure(&error, path))?
        .len();
    if size == 0 {
        return Err(failure("Parquet file size is 0 bytes"));
    }
    if size < FOOTER_SIZE {
        return Err(failure(format!(
            "Parquet file size is {size} bytes, smaller than the minimum file footer ({FOOTER_SIZE} bytes)"
        )));
    }
    let mut footer = [0u8; FOOTER_SIZE as usize];
    file.read_exact_at(&mut footer, size - FOOTER_SIZE)
        .map_err(|error| system_failure(&error, path))?;
    let (length, magic) = footer.split_at(4);
    if magic != PARQUET_MAGIC && magic != ENCRYPTED_MAGIC {
        return Err(failure(
            "Parquet magic bytes not found in footer. Either the file is corrupted or this is not a parquet file.",
        ));
    }
    let length = u32::from_le_bytes([length[0], length[1], length[2], length[3]]);
    if u64::from(length) > size - FOOTER_SIZE {
        return Err(failure(format!(
            "Parquet file size is {size} bytes, smaller than the size reported by footer's ({length}bytes)"
        )));
    }
    if magic == ENCRYPTED_MAGIC {
        return Err(failure(
            "Could not read encrypted metadata, no decryption found in reader's properties",
        ));
    }
    SerializedFileReader::new(file).map_err(failure)
}

/// The text a file's footer holds under the key `tab_columns`: its `--tab` column names,
/// joined by tabs.
fn stated_tab_columns(footer: &FileMetaData) -> Option<&str> {
    footer
        .key_value_metadata()
        .into_iter()
        .flatten()
        .find(|entry| entry.key == "tab_columns")
        .and_then(|entry| entry.value.as_deref())
}

/// The tab columns the top-level fields after the leading key columns hold, each under the
/// name the `--tab` format gives it, and the count of those key columns. The names are the
/// ones `stated` holds, the file's own list, where it names each such field: the writer
/// stores a name an earlier column carries in other letters behind the prefix too, and such
/// a name does not read back from the stored name. A file that states no list, or a list of
/// another length, names each column by `tab_column_name` of its field.
fn tab_columns(
    schema: &SchemaDescriptor,
    stated: Option<&str>,
) -> std::result::Result<(usize, Vec<String>), Failure> {
    let fields = schema.root_schema().get_fields();
    let skip = fields
        .iter()
        .take_while(|field| KEY_COLUMNS.contains(&field.name()))
        .count();
    if skip == fields.len() {
        return Err(failure("no tab columns in the Parquet schema"));
    }
    let stated: Vec<&str> = stated.map_or_else(Vec::new, |names| names.split('\t').collect());
    let names = if stated.len() == fields.len() - skip {
        stated.into_iter().map(str::to_string).collect()
    } else {
        fields[skip..]
            .iter()
            .map(|field| tab_column_name(field.name()).to_string())
            .collect()
    };
    Ok((skip, names))
}

/// The tab columns of a file with how each is read.
fn columns_of(
    path: &Path,
    schema: &SchemaDescriptor,
    skip: usize,
    names: &[String],
) -> std::result::Result<Vec<Column>, Failure> {
    let fields = schema.root_schema().get_fields();
    let mut first: usize = fields[..skip].iter().map(|field| leaves_below(field)).sum();
    let mut columns = Vec::with_capacity(names.len());
    for (field, name) in fields[skip..].iter().zip(names) {
        let shape = shape_of(schema, field, first).map_err(|held| {
            failure(format!(
                "{}: the column {name} holds {held}, which is not read",
                shown(path)
            ))
        })?;
        let canonical = canonical_of(name);
        columns.push(Column {
            shape,
            spelling: boolean_spelling(name),
            canonical: if stores_text(&shape) || prints_canonical(canonical, &shape) {
                Canonical::Verbatim
            } else {
                canonical
            },
        });
        first += leaves_below(field);
    }
    Ok(columns)
}

/// The values of some records of a leaf column: the levels of every value and NULL, and
/// the values that are not NULL.
#[derive(Default)]
struct LeafValues {
    def: Vec<i16>,
    rep: Vec<i16>,
    bools: Vec<bool>,
    int32s: Vec<i32>,
    int64s: Vec<i64>,
    floats: Vec<f32>,
    doubles: Vec<f64>,
    texts: Vec<ByteArray>,
}

impl LeafValues {
    fn clear(&mut self) {
        self.def.clear();
        self.rep.clear();
        self.bools.clear();
        self.int32s.clear();
        self.int64s.clear();
        self.floats.clear();
        self.doubles.clear();
        self.texts.clear();
    }
}

/// Reads the next `records` records of a leaf column into `values`, which is emptied first.
fn read_leaf(
    reader: &mut ColumnReader,
    records: usize,
    values: &mut LeafValues,
) -> parquet::errors::Result<()> {
    values.clear();
    let mut read = 0;
    while read < records {
        let left = records - read;
        let (def, rep) = (Some(&mut values.def), Some(&mut values.rep));
        let (records_read, _, levels_read) = match reader {
            ColumnReader::BoolColumnReader(reader) => {
                reader.read_records(left, def, rep, &mut values.bools)?
            }
            ColumnReader::Int32ColumnReader(reader) => {
                reader.read_records(left, def, rep, &mut values.int32s)?
            }
            ColumnReader::Int64ColumnReader(reader) => {
                reader.read_records(left, def, rep, &mut values.int64s)?
            }
            ColumnReader::FloatColumnReader(reader) => {
                reader.read_records(left, def, rep, &mut values.floats)?
            }
            ColumnReader::DoubleColumnReader(reader) => {
                reader.read_records(left, def, rep, &mut values.doubles)?
            }
            ColumnReader::ByteArrayColumnReader(reader) => {
                reader.read_records(left, def, rep, &mut values.texts)?
            }
            ColumnReader::Int96ColumnReader(_) | ColumnReader::FixedLenByteArrayColumnReader(_) => {
                return Err(parquet::errors::ParquetError::General(
                    "a column of this type is not read".to_string(),
                ))
            }
        };
        if records_read == 0 && levels_read == 0 {
            return Err(parquet::errors::ParquetError::General(format!(
                "a column ends {left} rows before its row group"
            )));
        }
        read += records_read;
    }
    Ok(())
}

/// Why the cells of some rows were not printed.
enum Unread {
    /// A string is not UTF-8: the text of the interpreter's refusal to decode it.
    Text(String),
    /// The reader's refusal of the file's pages.
    File(parquet::errors::ParquetError),
}

/// A member of a positions struct as the interpreter holds it.
#[derive(Clone, Copy)]
enum Member<'a> {
    Bool(bool),
    Whole(i128),
    Real(f64),
    Text(&'a str),
}

impl Member<'_> {
    /// `a == b` of the interpreter: numbers by value whatever their type.
    fn equals(&self, other: &Member) -> bool {
        let number = |member: &Member| match *member {
            Member::Bool(flag) => Some(Ok(i128::from(flag))),
            Member::Whole(whole) => Some(Ok(whole)),
            Member::Real(real) => Some(Err(real)),
            Member::Text(_) => None,
        };
        match (number(self), number(other)) {
            (None, None) => matches!((self, other), (Member::Text(a), Member::Text(b)) if a == b),
            (Some(Ok(a)), Some(Ok(b))) => a == b,
            (Some(Err(a)), Some(Err(b))) => a == b,
            (Some(Ok(whole)), Some(Err(real))) | (Some(Err(real)), Some(Ok(whole))) => {
                // Every whole number a column holds is one a double compares exactly with.
                real.fract() == 0.0 && real.abs() < 1e30 && real as i128 == whole
            }
            _ => false,
        }
    }

    /// `str(member)` of the interpreter.
    fn write(&self, out: &mut String) {
        match *self {
            Member::Bool(flag) => out.push_str(if flag { "True" } else { "False" }),
            Member::Whole(whole) => {
                let _ = write!(out, "{whole}");
            }
            Member::Real(real) => pyfmt::write_repr_f64(out, real),
            Member::Text(text) => out.push_str(text),
        }
    }
}

/// The cells of some rows of one column, one after another.
#[derive(Default)]
struct Cells {
    text: String,
    /// Where every cell ends.
    ends: Vec<usize>,
}

impl Cells {
    fn clear(&mut self) {
        self.text.clear();
        self.ends.clear();
    }

    fn end_cell(&mut self) {
        self.ends.push(self.text.len());
    }

    fn cell(&self, row: usize) -> &str {
        let start = if row == 0 { 0 } else { self.ends[row - 1] };
        &self.text[start..self.ends[row]]
    }
}

fn decoded(bytes: &ByteArray) -> std::result::Result<&str, Unread> {
    pyfmt::decode_utf8(bytes.data()).map_err(Unread::Text)
}

/// Appends the text of the value `at` of a leaf's values that are not NULL, a boolean in
/// `spelling`.
fn write_value(
    out: &mut String,
    kind: LeafKind,
    values: &LeafValues,
    at: usize,
    spelling: (&str, &str),
) -> std::result::Result<(), Unread> {
    match kind {
        LeafKind::Bool => out.push_str(if values.bools[at] {
            spelling.0
        } else {
            spelling.1
        }),
        LeafKind::Int32 => {
            let _ = write!(out, "{}", values.int32s[at]);
        }
        LeafKind::UInt32 => {
            let _ = write!(out, "{}", values.int32s[at] as u32);
        }
        LeafKind::Int64 => {
            let _ = write!(out, "{}", values.int64s[at]);
        }
        LeafKind::UInt64 => {
            let _ = write!(out, "{}", values.int64s[at] as u64);
        }
        LeafKind::Float => write_shortest_decimal(out, f64::from(values.floats[at])),
        LeafKind::Double => write_shortest_decimal(out, values.doubles[at]),
        LeafKind::Text => out.push_str(decoded(&values.texts[at])?),
    }
    Ok(())
}

/// The value `at` of a leaf's values that are not NULL, as a member of a positions struct.
fn member_of(
    kind: LeafKind,
    values: &LeafValues,
    at: usize,
) -> std::result::Result<Member<'_>, Unread> {
    Ok(match kind {
        LeafKind::Bool => Member::Bool(values.bools[at]),
        LeafKind::Int32 => Member::Whole(i128::from(values.int32s[at])),
        LeafKind::UInt32 => Member::Whole(i128::from(values.int32s[at] as u32)),
        LeafKind::Int64 => Member::Whole(i128::from(values.int64s[at])),
        LeafKind::UInt64 => Member::Whole(i128::from(values.int64s[at] as u64)),
        LeafKind::Float => Member::Real(f64::from(values.floats[at])),
        LeafKind::Double => Member::Real(values.doubles[at]),
        LeafKind::Text => Member::Text(decoded(&values.texts[at])?),
    })
}

/// Appends the cell of a positions struct: `start` where `start` and `end` agree,
/// `start-end` otherwise with `?` for a NULL member, `-` where both are NULL, then `/total`.
fn write_positions(
    out: &mut String,
    start: Option<Member>,
    end: Option<Member>,
    total: Option<Member>,
) {
    match (start, end) {
        (None, None) => out.push_str(NULL_TEXT),
        (Some(start), Some(end)) if start.equals(&end) => start.write(out),
        _ => {
            for (at, member) in [start, end].into_iter().enumerate() {
                if at > 0 {
                    out.push('-');
                }
                match member {
                    Some(member) => member.write(out),
                    None => out.push('?'),
                }
            }
        }
    }
    if let Some(total) = total {
        out.push('/');
        total.write(out);
    }
}

/// The open leaf columns of one tab column in a row group, with the values last read.
struct ColumnReaders {
    readers: [Option<ColumnReader>; 3],
    values: [LeafValues; 3],
}

/// Reads the next `rows` rows of a column and prints their cells into `cells`.
fn read_cells(
    column: &Column,
    open: &mut ColumnReaders,
    rows: usize,
    cells: &mut Cells,
) -> std::result::Result<(), Unread> {
    cells.clear();
    for (reader, values) in open.readers.iter_mut().zip(open.values.iter_mut()) {
        if let Some(reader) = reader {
            read_leaf(reader, rows, values).map_err(Unread::File)?;
        }
    }
    let short = || {
        Unread::File(parquet::errors::ParquetError::General(
            "a column holds fewer values than its levels name".to_string(),
        ))
    };
    match column.shape {
        Shape::NoPositions => {
            for _ in 0..rows {
                cells.text.push_str(NULL_TEXT);
                cells.end_cell();
            }
        }
        Shape::Scalar(leaf) => {
            let values = &open.values[0];
            let mut at = 0;
            for row in 0..rows {
                if leaf.max_def > 0
                    && values.def.get(row).copied().ok_or_else(short)? < leaf.max_def
                {
                    cells.text.push_str(NULL_TEXT);
                } else {
                    write_value(&mut cells.text, leaf.kind, values, at, column.spelling)?;
                    at += 1;
                }
                cells.end_cell();
            }
        }
        Shape::List { leaf, present } => {
            let values = &open.values[0];
            let mut at = 0;
            let mut level = 0;
            for _ in 0..rows {
                let def = values.def.get(level).copied().ok_or_else(short)?;
                if def < present {
                    cells.text.push_str(NULL_TEXT);
                    level += 1;
                } else if def == present && leaf.max_def > present {
                    level += 1;
                } else {
                    let mut first = true;
                    while first || values.rep.get(level).is_some_and(|rep| *rep > 0) {
                        if !first {
                            cells.text.push(',');
                        }
                        first = false;
                        if values.def[level] < leaf.max_def {
                            cells.text.push_str(NULL_TEXT);
                        } else {
                            write_value(&mut cells.text, leaf.kind, values, at, FLAG_SPELLING)?;
                            at += 1;
                        }
                        level += 1;
                    }
                }
                cells.end_cell();
            }
        }
        Shape::Positions { members, present } => {
            let mut at = [0usize; 3];
            for row in 0..rows {
                let mut held: [Option<Member>; 3] = [None; 3];
                let mut null = false;
                for (place, member) in members.iter().enumerate() {
                    let Some(leaf) = member else { continue };
                    let values = &open.values[place];
                    let def = if leaf.max_def > 0 {
                        values.def.get(row).copied().ok_or_else(short)?
                    } else {
                        0
                    };
                    if def < present {
                        null = true;
                    } else if def == leaf.max_def {
                        held[place] = Some(member_of(leaf.kind, values, at[place])?);
                        at[place] += 1;
                    }
                }
                if null {
                    cells.text.push_str(NULL_TEXT);
                } else {
                    write_positions(&mut cells.text, held[0], held[1], held[2]);
                }
                cells.end_cell();
            }
        }
    }
    Ok(())
}

/// Where in the adapter's reading of a file a failure is met: the batch of rows, the cells
/// of every column of the batch read before the first of its rows is printed.
type Met = (usize, usize, usize, usize);

/// The failure of the rows from `row` on at column `column`, and where it is met.
fn unread(column: usize, row: usize, unread: Unread) -> (Met, Failure) {
    let met = (row / BATCH_ROWS, 0, column, row);
    match unread {
        Unread::Text(refusal) => (met, failure(refusal)),
        Unread::File(refusal) => (met, failure(refusal)),
    }
}

/// The leaf columns of every tab column of one row group, opened at the row group's first
/// row, `first_row` of its file.
fn open_columns(
    reader: &SerializedFileReader<File>,
    part: &PartFile,
    group: usize,
    first_row: usize,
) -> std::result::Result<Vec<ColumnReaders>, (Met, Failure)> {
    let row_group = reader
        .get_row_group(group)
        .map_err(|refusal| unread(0, first_row, Unread::File(refusal)))?;
    let mut open: Vec<ColumnReaders> = Vec::with_capacity(part.columns.len());
    for (at, column) in part.columns.iter().enumerate() {
        let leaves: [Option<Leaf>; 3] = match column.shape {
            Shape::Scalar(leaf) | Shape::List { leaf, .. } => [Some(leaf), None, None],
            Shape::Positions { members, .. } => members,
            Shape::NoPositions => [None; 3],
        };
        let mut readers: [Option<ColumnReader>; 3] = [None, None, None];
        for (place, leaf) in leaves.iter().enumerate() {
            if let Some(leaf) = leaf {
                let reader = row_group
                    .get_column_reader(leaf.index)
                    .map_err(|refusal| unread(at, first_row, Unread::File(refusal)))?;
                readers[place] = Some(reader);
            }
        }
        open.push(ColumnReaders {
            readers,
            values: Default::default(),
        });
    }
    Ok(open)
}

/// Reads the next `rows` rows of a row group's open columns, the first of them `row` of
/// its file, and appends their lines to `lines`: every column's cells are read before the
/// first row is printed.
#[allow(clippy::too_many_arguments)]
fn print_rows(
    part: &PartFile,
    open: &mut [ColumnReaders],
    row: usize,
    rows: usize,
    header: &[Canonical],
    canonical: bool,
    cells: &mut [Cells],
    scratch: &mut String,
    lines: &mut Lines,
) -> std::result::Result<(), (Met, Failure)> {
    for (at, column) in part.columns.iter().enumerate() {
        read_cells(column, &mut open[at], rows, &mut cells[at])
            .map_err(|refusal| unread(at, row, refusal))?;
    }
    // A cell that holds a tab is several cells of the canonical form's line.
    let split = canonical && cells.iter().any(|cells| cells.text.contains('\t'));
    for at in 0..rows {
        let refused = |refusal: String| ((row / BATCH_ROWS, 1, 0, row + at), failure(refusal));
        if split {
            scratch.clear();
            for (column, cells) in cells.iter().enumerate() {
                if column > 0 {
                    scratch.push('\t');
                }
                scratch.push_str(cells.cell(at));
            }
            write_canonical_line(&mut lines.text, header, scratch).map_err(refused)?;
        } else {
            for (column, cells) in cells.iter().enumerate() {
                if column > 0 {
                    lines.text.push('\t');
                }
                let form = if canonical {
                    part.columns[column].canonical
                } else {
                    Canonical::Verbatim
                };
                write_canonical_cell(&mut lines.text, form, cells.cell(at)).map_err(refused)?;
            }
        }
        lines.end_line();
    }
    Ok(())
}

/// The lines of one row group, or the failure met first in it and where.
fn row_group_lines(
    reader: &SerializedFileReader<File>,
    part: &PartFile,
    group: usize,
    first_row: usize,
    header: &[Canonical],
    canonical: bool,
) -> std::result::Result<Lines, (Met, Failure)> {
    let mut open = open_columns(reader, part, group, first_row)?;
    let mut lines = Lines::new();
    let mut cells: Vec<Cells> = Vec::new();
    cells.resize_with(part.columns.len(), Cells::default);
    let mut scratch = String::new();
    let mut row = first_row;
    let end = first_row + part.row_groups[group];
    while row < end {
        let rows = (end - row).min(BATCH_ROWS - row % BATCH_ROWS);
        print_rows(
            part,
            &mut open,
            row,
            rows,
            header,
            canonical,
            &mut cells,
            &mut scratch,
            &mut lines,
        )?;
        row += rows;
    }
    if canonical {
        lines.sort();
    }
    Ok(lines)
}

/// The Parquet files of an output with how the tab columns of each are read. Every file
/// is opened before the columns of the first are read; none is held open.
fn part_files(path: &Path, is_dir: bool) -> std::result::Result<Vec<PartFile>, Failure> {
    let files = parquet_files(path, is_dir)?;
    let mut opened = Vec::with_capacity(files.len());
    let mut stated = Vec::with_capacity(files.len());
    for file in &files {
        let reader = open(file)?;
        let metadata = reader.metadata();
        stated.push(stated_tab_columns(metadata.file_metadata()).map(str::to_string));
        let row_groups: Vec<usize> = metadata
            .row_groups()
            .iter()
            .map(|group| group.num_rows().max(0) as usize)
            .collect();
        opened.push((metadata.file_metadata().schema_descr_ptr(), row_groups));
    }
    let mut parts: Vec<PartFile> = Vec::with_capacity(files.len());
    for (file, (schema, row_groups)) in files.iter().zip(&opened) {
        let (key_columns, names) = tab_columns(schema, stated[parts.len()].as_deref())?;
        if parts.first().is_some_and(|first| first.names != names) {
            return Err(failure(format!(
                "{}: its tab columns differ from {}'s",
                shown(file),
                shown(&files[0])
            )));
        }
        parts.push(PartFile {
            path: file.clone(),
            columns: Vec::new(),
            key_columns,
            names,
            row_groups: row_groups.clone(),
        });
    }
    for (part, (schema, _)) in parts.iter_mut().zip(&opened) {
        part.columns = columns_of(&part.path, schema, part.key_columns, &part.names)?;
    }
    Ok(parts)
}

/// What the canonical form does to a cell of `part` by its place, for a line that is
/// split again.
fn forms_by_place(part: &PartFile) -> Vec<Canonical> {
    let by_place = part.names.iter().zip(&part.columns);
    by_place
        .map(|(name, column)| {
            if stores_text(&column.shape) {
                Canonical::Verbatim
            } else {
                canonical_of(name)
            }
        })
        .collect()
}

/// The text of a Parquet output: its tab columns and its rows as `--tab` lines, in file
/// order or in the canonical form.
fn parquet_output(
    path: &Path,
    is_dir: bool,
    canonical: bool,
) -> std::result::Result<String, Failure> {
    let parts = part_files(path, is_dir)?;
    let names = parts[0].names.clone();
    let headers: Vec<Vec<Canonical>> = parts.iter().map(forms_by_place).collect();
    // One task a row group: the file, the row group and the rows of the file before it.
    let mut tasks: Vec<(usize, usize, usize)> = Vec::new();
    for (file, part) in parts.iter().enumerate() {
        let mut first_row = 0;
        for (group, rows) in part.row_groups.iter().enumerate() {
            tasks.push((file, group, first_row));
            first_row += rows;
        }
    }
    let read = on_workers(
        tasks.len(),
        || None::<(usize, SerializedFileReader<File>)>,
        |opened, task| {
            let (file, group, first_row) = tasks[task];
            let part = &parts[file];
            if opened.as_ref().is_none_or(|(held, _)| *held != file) {
                let reader = open(&part.path).map_err(|refusal| ((file, (0, 0, 0, 0)), refusal))?;
                *opened = Some((file, reader));
            }
            let (_, reader) = opened.as_ref().expect("the file was opened");
            row_group_lines(reader, part, group, first_row, &headers[file], canonical)
                .map_err(|(met, refusal)| ((file, met), refusal))
        },
    );
    let mut lines = Vec::with_capacity(read.len());
    let mut first: Option<((usize, Met), Failure)> = None;
    for result in read {
        match result {
            Ok(part) => lines.push(part),
            Err((met, refusal)) => {
                if first.as_ref().is_none_or(|(earliest, _)| met < *earliest) {
                    first = Some((met, refusal));
                }
            }
        }
    }
    match first {
        Some((_, refusal)) => Err(refusal),
        None => Ok(joined(&names, lines, canonical)),
    }
}

/// Whether the adapter's option parser takes its one argument for an option it does not
/// have. An argument that opens with `-` is a path to that parser only where it is `-`
/// alone, opens as a negative number does or holds a space; one that spells an option of
/// the parser, whole or by its first letters, is read as that option, which is not
/// followed here.
fn read_as_an_unknown_option(argument: &str) -> bool {
    let Some(after) = argument.strip_prefix('-') else {
        return false;
    };
    let mut chars = after.chars();
    let opens_a_number = match chars.next() {
        None => return false,
        Some('.') => chars.next().is_some_and(pyfmt::char_is_decimal),
        Some(c) => pyfmt::char_is_decimal(c),
    };
    let spells_an_option = match after.strip_prefix('-') {
        Some(long) => {
            let name = long.split_once('=').map_or(long, |(name, _)| name);
            ["help", "input", "output", "canonical"]
                .iter()
                .any(|option| option.starts_with(name))
        }
        // Before a `=` the dash alone opens every option.
        None => after.starts_with(['h', 'i', 'o', '=']),
    };
    !(opens_a_number || argument.contains(' ') || spells_an_option)
}

/// What the adapter writes to its standard output for `path`; a `--tab` output keeps the
/// columns named in `as_written` as written.
fn adapter_output(
    path: &Path,
    canonical: bool,
    as_written: &[String],
) -> std::result::Result<String, Failure> {
    let argument = path_text(path);
    if read_as_an_unknown_option(&argument) {
        return Err(Failure(format!(
            "{ADAPTER_USAGE}\n{ADAPTER_PROGRAM}: error: unrecognized arguments: {}",
            backslashreplace(&argument)
        )));
    }
    let found = std::fs::metadata(path)
        .map_err(|_| failure(format!("{}: no such file or directory", shown(path))))?;
    if !(found.is_file() || found.is_dir()) {
        return Err(failure(format!(
            "{}: neither a regular file nor a directory",
            shown(path)
        )));
    }
    let parquet = found.is_dir() || {
        let mut head = Vec::with_capacity(PARQUET_MAGIC.len());
        File::open(path)
            .and_then(|file| file.take(PARQUET_MAGIC.len() as u64).read_to_end(&mut head))
            .map_err(|error| system_failure(&error, path))?;
        head == PARQUET_MAGIC
    };
    if parquet {
        parquet_output(path, found.is_dir(), canonical)
    } else if canonical {
        tab_output(path, as_written)
    } else {
        Err(failure(format!("{}: not a Parquet file", shown(path))))
    }
}

/// A Parquet output (a file or a directory of part files), or with `canonical` a
/// `--tab` output too, read through the adapter as tab text. The canonical form is the
/// header line, then the rows sorted, every DOUBLE cell as its shortest decimal: a Parquet
/// file stores a frequency as a number, so the `0.0010` and `1.47e-05` spellings of a
/// `--tab` output come back as `0.001` and `0.0000147`, and reading the tab side
/// through the same form compares the two as numbers and every other cell verbatim. Without
/// `canonical` the rows come in file order, which self-parity needs to pair the rows of a
/// record that share a key with the default format's rows in their order.
///
/// The text is the adapter's standard output as the comparator reads it, every `\r\n`
/// and every other `\r` a `\n`.
pub fn parquet_as_tab(path: &Path, canonical: bool) -> Result<String> {
    as_the_comparator_reads(path, adapter_output(path, canonical, &[]))
}

/// The `--tab` output at `path` in the canonical form, read beside the Parquet output
/// `parquet` it is compared with: a column that output stores as strings stays as written,
/// where [`parquet_as_tab`] reprints a column of a `--tab` output by its name alone.
pub fn tab_beside_parquet(path: &Path, parquet: &Path) -> Result<String> {
    as_the_comparator_reads(path, adapter_output(path, true, &text_columns(parquet)))
}

/// The names of the tab columns the Parquet output at `path` stores as strings: each
/// column of the name, in every file. An output that cannot be read names none, and the
/// comparison meets that failure where it reads the output itself.
pub fn text_columns(path: &Path) -> Vec<String> {
    let mut strings: Vec<String> = Vec::new();
    let mut others: Vec<String> = Vec::new();
    let Ok(found) = std::fs::metadata(path) else {
        return strings;
    };
    // A path that is neither a regular file nor a directory is refused where the output is
    // read, and opening a named pipe waits for a writer.
    if !(found.is_file() || found.is_dir()) {
        return strings;
    }
    let Ok(files) = parquet_files(path, found.is_dir()) else {
        return strings;
    };
    for file in &files {
        let Ok(reader) = open(file) else {
            return Vec::new();
        };
        let schema = reader.metadata().file_metadata().schema_descr_ptr();
        let stated = stated_tab_columns(reader.metadata().file_metadata());
        let Ok((skip, names)) = tab_columns(&schema, stated) else {
            return Vec::new();
        };
        let Ok(columns) = columns_of(file, &schema, skip, &names) else {
            return Vec::new();
        };
        for (name, column) in names.into_iter().zip(&columns) {
            let held = if stores_text(&column.shape) {
                &mut strings
            } else {
                &mut others
            };
            if !held.contains(&name) {
                held.push(name);
            }
        }
    }
    strings.retain(|name| !others.contains(name));
    strings
}

/// The adapter's standard output as the comparator reads it, or its failure as the
/// comparator's error.
fn as_the_comparator_reads(
    path: &Path,
    output: std::result::Result<String, Failure>,
) -> Result<String> {
    match output {
        Ok(text) => Ok(universal_newlines(text)),
        Err(Failure(stderr)) => Err(Error::comparator(format!(
            "the Parquet adapter failed on {}: {}",
            path_text(path),
            pyfmt::split_whitespace(&stderr).join(" ")
        ))),
    }
}

/// What the readings of a `CanonicalText` met.
#[derive(Default)]
struct Seen {
    failed: AtomicBool,
    split: AtomicBool,
    header_row: AtomicBool,
}

/// Whether a line of `lines` lines of text holds a line end of its own: a `\r`, or a `\n`
/// before its last byte.
fn holds_a_line_end(text: &str, lines: usize) -> bool {
    let bytes = text.as_bytes();
    memchr::memchr(b'\r', bytes).is_some() || memchr::memchr_iter(b'\n', bytes).count() != lines
}

/// Whether a line of `text`, whose lines each end in `\n`, opens with `#`.
fn holds_a_header_line(text: &str) -> bool {
    text.starts_with('#') || memchr::memmem::find(text.as_bytes(), b"\n#").is_some()
}

/// The rows of a Parquet output in the canonical form, the files in the order of their
/// paths and the rows of each in the file's order.
struct ParquetRows {
    parts: Vec<PartFile>,
    /// What the canonical form does to a cell of each file by its place.
    headers: Vec<Vec<Canonical>>,
    started: bool,
    /// The file in work, its next row group and that row group's first row.
    file: usize,
    group: usize,
    first_row: usize,
    reader: Option<SerializedFileReader<File>>,
    /// The open columns of the row group in work, its next row and the row after its last.
    open: Option<(Vec<ColumnReaders>, usize, usize)>,
    cells: Vec<Cells>,
    scratch: String,
    lines: Lines,
}

impl ParquetRows {
    fn open(path: &Path, is_dir: bool) -> std::result::Result<ParquetRows, Failure> {
        let parts = part_files(path, is_dir)?;
        let headers = parts.iter().map(forms_by_place).collect();
        Ok(ParquetRows {
            parts,
            headers,
            started: false,
            file: 0,
            group: 0,
            first_row: 0,
            reader: None,
            open: None,
            cells: Vec::new(),
            scratch: String::new(),
            lines: Lines::new(),
        })
    }

    /// Appends the header line or the next rows to `out`; `false` after the last row.
    fn next(&mut self, out: &mut String, seen: &Seen) -> std::result::Result<bool, Failure> {
        if !self.started {
            self.started = true;
            let from = out.len();
            out.push('#');
            out.push_str(&self.parts[0].names.join("\t"));
            out.push('\n');
            if holds_a_line_end(&out[from..], 1) {
                seen.split.store(true, Ordering::Relaxed);
            }
            return Ok(true);
        }
        loop {
            if let Some((open, row, end)) = &mut self.open {
                if *row < *end {
                    let part = &self.parts[self.file];
                    let rows = (*end - *row).min(TEXT_ROWS);
                    self.cells.resize_with(part.columns.len(), Cells::default);
                    self.lines.clear();
                    print_rows(
                        part,
                        open,
                        *row,
                        rows,
                        &self.headers[self.file],
                        true,
                        &mut self.cells,
                        &mut self.scratch,
                        &mut self.lines,
                    )
                    .map_err(|(_, refusal)| refusal)?;
                    *row += rows;
                    if holds_a_line_end(&self.lines.text, rows) {
                        seen.split.store(true, Ordering::Relaxed);
                    }
                    if holds_a_header_line(&self.lines.text) {
                        seen.header_row.store(true, Ordering::Relaxed);
                    }
                    out.push_str(&self.lines.text);
                    return Ok(true);
                }
                self.open = None;
            }
            let Some(part) = self.parts.get(self.file) else {
                return Ok(false);
            };
            if self.group == part.row_groups.len() {
                self.file += 1;
                self.group = 0;
                self.first_row = 0;
                self.reader = None;
                continue;
            }
            let reader = match &self.reader {
                Some(reader) => reader,
                None => self.reader.insert(open(&part.path)?),
            };
            let columns = open_columns(reader, part, self.group, self.first_row)
                .map_err(|(_, refusal)| refusal)?;
            let rows = part.row_groups[self.group];
            self.open = Some((columns, self.first_row, self.first_row + rows));
            self.first_row += rows;
            self.group += 1;
        }
    }
}

/// The bytes of a `--tab` output, a gzip file inflated.
enum TabBytes {
    Plain(File),
    Gzip(Box<Gunzip>),
}

/// The rows of a `--tab` output in the canonical form, in the file's order. The column
/// line stands before the first data line of the output these rows are read from; the
/// sorted form takes it wherever it stands.
struct TabRows {
    path: PathBuf,
    /// The columns whose cells stay as written, as `tab_output` takes them.
    as_written: Vec<String>,
    bytes: TabBytes,
    /// The bytes read whose lines are still to print: the start of a line, once the lines
    /// before it are printed.
    held: Vec<u8>,
    ended: bool,
    columns: Option<Vec<Canonical>>,
}

impl TabRows {
    fn open(path: &Path, as_written: Vec<String>) -> std::result::Result<TabRows, Failure> {
        let mut head = Vec::with_capacity(GZIP_MAGIC.len());
        File::open(path)
            .and_then(|file| file.take(GZIP_MAGIC.len() as u64).read_to_end(&mut head))
            .map_err(|error| system_failure(&error, path))?;
        let bytes = if head == GZIP_MAGIC {
            let named = path.to_path_buf();
            let file: gzip::Open = Box::new(move || Ok(Box::new(File::open(&named)?)));
            TabBytes::Gzip(Box::new(
                Gunzip::new(file).map_err(|error| system_failure(&error, path))?,
            ))
        } else {
            TabBytes::Plain(File::open(path).map_err(|error| system_failure(&error, path))?)
        };
        Ok(TabRows {
            path: path.to_path_buf(),
            as_written,
            bytes,
            held: Vec::new(),
            ended: false,
            columns: None,
        })
    }

    /// Reads up to `TEXT_BYTES` more bytes onto `held`.
    fn fill(&mut self) -> std::result::Result<(), Failure> {
        let from = self.held.len();
        self.held.resize(from + TEXT_BYTES, 0);
        let mut filled = from;
        let read = match &mut self.bytes {
            TabBytes::Gzip(gunzip) => match gunzip.read(&mut self.held, &mut filled) {
                Ok(more) => Ok(more),
                Err(gzip::Failure::System(error)) => Err(system_failure(&error, &self.path)),
                Err(gzip::Failure::Text(refusal)) => Err(failure(refusal)),
            },
            TabBytes::Plain(file) => loop {
                if filled == self.held.len() {
                    break Ok(true);
                }
                match file.read(&mut self.held[filled..]) {
                    Ok(0) => break Ok(false),
                    Ok(read) => filled += read,
                    Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                    Err(error) => break Err(system_failure(&error, &self.path)),
                }
            },
        };
        self.held.truncate(filled);
        self.ended = !read?;
        Ok(())
    }

    /// Appends the canonical form of one line of the output to `out`: nothing for a `##`
    /// line and for an empty line.
    fn print(
        &mut self,
        line: &[u8],
        out: &mut String,
        seen: &Seen,
    ) -> std::result::Result<(), Failure> {
        let line = std::str::from_utf8(line)
            .map_err(|_| failure(format!("{}: a line is not UTF-8", shown(&self.path))))?
            .trim_end_matches('\r');
        if line.starts_with("##") || line.is_empty() {
            return Ok(());
        }
        let from = out.len();
        if let Some(columns) = line.strip_prefix('#') {
            if self.columns.is_some() {
                return Err(failure(format!("{}: two header lines", shown(&self.path))));
            }
            let form = |column: &str| {
                if self.as_written.iter().any(|written| written == column) {
                    Canonical::Verbatim
                } else {
                    canonical_of(column)
                }
            };
            self.columns = Some(columns.split('\t').map(form).collect());
            out.push_str(line);
        } else {
            let Some(columns) = &self.columns else {
                return Err(failure(format!(
                    "{}: a data line precedes the # header line",
                    shown(&self.path)
                )));
            };
            write_canonical_line(out, columns, line).map_err(failure)?;
        }
        out.push('\n');
        if holds_a_line_end(&out[from..], 1) {
            seen.split.store(true, Ordering::Relaxed);
        }
        Ok(())
    }

    /// Appends the next lines to `out`; `false` after the last line.
    fn next(&mut self, out: &mut String, seen: &Seen) -> std::result::Result<bool, Failure> {
        loop {
            if !self.ended {
                self.fill()?;
            }
            let whole = match memchr::memrchr(b'\n', &self.held) {
                Some(last) => last + 1,
                None if self.ended => self.held.len(),
                // A line longer than the bytes read so far.
                None => continue,
            };
            let held = std::mem::take(&mut self.held);
            let mut lines = held[..whole].split(|byte| *byte == b'\n').peekable();
            while let Some(line) = lines.next() {
                // What follows the last line end is no line.
                if lines.peek().is_some() || !line.is_empty() {
                    self.print(line, out, seen)?;
                }
            }
            self.held = held;
            self.held.drain(..whole);
            if self.ended && self.held.is_empty() {
                if self.columns.is_none() {
                    return Err(failure(format!(
                        "{}: no # header line; a tab output needs its header",
                        shown(&self.path)
                    )));
                }
                return Ok(false);
            }
            return Ok(true);
        }
    }
}

enum Rows {
    Parquet(Box<ParquetRows>),
    Tab(Box<TabRows>),
}

impl Rows {
    /// The rows of the output at `path`, a Parquet output or a `--tab` output, as
    /// `adapter_output` tells the two apart. A `--tab` output keeps as written the columns
    /// that `beside`, the Parquet output it is compared with, stores as strings.
    fn open(path: &Path, beside: Option<&Path>) -> std::result::Result<Rows, Failure> {
        let argument = path_text(path);
        if read_as_an_unknown_option(&argument) {
            return Err(failure(format!(
                "{}: read as an option",
                backslashreplace(&argument)
            )));
        }
        let found = std::fs::metadata(path)
            .map_err(|_| failure(format!("{}: no such file or directory", shown(path))))?;
        if !(found.is_file() || found.is_dir()) {
            return Err(failure(format!(
                "{}: neither a regular file nor a directory",
                shown(path)
            )));
        }
        let parquet = found.is_dir() || {
            let mut head = Vec::with_capacity(PARQUET_MAGIC.len());
            File::open(path)
                .and_then(|file| file.take(PARQUET_MAGIC.len() as u64).read_to_end(&mut head))
                .map_err(|error| system_failure(&error, path))?;
            head == PARQUET_MAGIC
        };
        if parquet {
            Ok(Rows::Parquet(Box::new(ParquetRows::open(
                path,
                found.is_dir(),
            )?)))
        } else {
            let as_written = beside.map(text_columns).unwrap_or_default();
            Ok(Rows::Tab(Box::new(TabRows::open(path, as_written)?)))
        }
    }

    fn next(&mut self, out: &mut String, seen: &Seen) -> std::result::Result<bool, Failure> {
        match self {
            Rows::Parquet(rows) => rows.next(out, seen),
            Rows::Tab(rows) => rows.next(out, seen),
        }
    }
}

/// One reading of a `CanonicalText`.
struct Reading {
    rows: Rows,
    seen: Arc<Seen>,
    text: String,
    at: usize,
    ended: bool,
}

impl Read for Reading {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        while self.at == self.text.len() {
            if self.ended {
                return Ok(0);
            }
            self.text.clear();
            self.at = 0;
            match self.rows.next(&mut self.text, &self.seen) {
                Ok(more) => self.ended = !more,
                Err(Failure(refusal)) => {
                    self.seen.failed.store(true, Ordering::Relaxed);
                    self.ended = true;
                    self.text.clear();
                    return Err(io::Error::other(refusal));
                }
            }
        }
        let count = out.len().min(self.text.len() - self.at);
        out[..count].copy_from_slice(&self.text.as_bytes()[self.at..self.at + count]);
        self.at += count;
        Ok(count)
    }
}

/// The canonical form of an output with its rows in the output's own order: the header
/// line, then one line a row, a Parquet output's files in the order of their paths and a
/// `--tab` output's data lines in the file's order. Sorted after their header line, the
/// lines are the lines the sorted canonical form gives for the output.
///
/// The text is written anew at each reading, a few thousand rows at a time, and is never
/// held whole. A reading fails where the sorted form fails, and also on a `--tab` output
/// whose column line follows a data line, which the sorted form reads; its failure is
/// worded here and not as the adapter words it, so a caller that needs the adapter's
/// words reads the sorted form. `unread` says whether a reading failed, and `split`
/// whether a reading met a row that holds a line end, which a reader of lines takes for
/// two rows. `header_row` says whether a reading met a row of a Parquet output that opens
/// with `#`, which a reader of lines takes for a header line: the rows that follow it in
/// the output's order are not the rows that follow it in the sorted form.
pub struct CanonicalText {
    path: PathBuf,
    beside: Option<PathBuf>,
    seen: Arc<Seen>,
}

impl CanonicalText {
    /// The text of the output at `path`, whose sorted form is `parquet_as_tab`'s.
    pub fn of(path: &Path) -> Arc<CanonicalText> {
        Arc::new(CanonicalText {
            path: path.to_path_buf(),
            beside: None,
            seen: Arc::new(Seen::default()),
        })
    }

    /// The text of the output at `path` read beside the Parquet output `parquet` it is
    /// compared with, whose sorted form is `tab_beside_parquet`'s.
    pub fn beside(path: &Path, parquet: &Path) -> Arc<CanonicalText> {
        Arc::new(CanonicalText {
            path: path.to_path_buf(),
            beside: Some(parquet.to_path_buf()),
            seen: Arc::new(Seen::default()),
        })
    }

    /// Whether a reading failed.
    pub fn unread(&self) -> bool {
        self.seen.failed.load(Ordering::Relaxed)
    }

    /// Whether a reading met a row or a column line that holds a line end.
    pub fn split(&self) -> bool {
        self.seen.split.load(Ordering::Relaxed)
    }

    /// Whether a reading met a row of a Parquet output that opens with `#`.
    pub fn header_row(&self) -> bool {
        self.seen.header_row.load(Ordering::Relaxed)
    }
}

impl Text for CanonicalText {
    fn open(&self) -> io::Result<Box<dyn Read + Send>> {
        match Rows::open(&self.path, self.beside.as_deref()) {
            Ok(rows) => Ok(Box::new(Reading {
                rows,
                seen: Arc::clone(&self.seen),
                text: String::new(),
                at: 0,
                ended: false,
            })),
            Err(Failure(refusal)) => {
                self.seen.failed.store(true, Ordering::Relaxed);
                Err(io::Error::other(refusal))
            }
        }
    }

    fn failed(&self, detail: &str) -> Error {
        Error::comparator(format!(
            "cannot read the canonical form of {}: {detail}",
            path_text(&self.path)
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn canonical(column: &str, text: &str) -> String {
        let mut out = String::new();
        write_canonical_cell(&mut out, canonical_of(column), text).expect("the cell is rewritten");
        out
    }

    fn shortest(value: f64) -> String {
        let mut out = String::new();
        write_shortest_decimal(&mut out, value);
        out
    }

    #[test]
    fn the_shortest_decimal_is_the_writers_text() {
        let cases = [
            (0.095, "0.095"),
            (1.0, "1"),
            (9.8e-05, "0.000098"),
            (0.0, "0"),
            (-0.0, "-0"),
            (1e16, "10000000000000000"),
            (1e21, "1000000000000000000000"),
            (0.1 + 0.2, "0.30000000000000004"),
            (f64::NAN, "NaN"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            // Two shortest decimals equally close: the one of larger magnitude.
            (-126430982921026.62, "-126430982921026.63"),
            (126430982921026.62, "126430982921026.63"),
        ];
        for (value, text) in cases {
            assert_eq!(shortest(value), text);
            if value.is_finite() {
                assert_eq!(
                    text.parse::<f64>().unwrap().to_bits(),
                    value.to_bits(),
                    "{text}"
                );
            }
        }
        assert_eq!(shortest(5e-324), format!("0.{}5", "0".repeat(323)));
        assert_eq!(shortest(f64::from(0.1f32)), "0.10000000149011612");
    }

    #[test]
    fn the_number_grammar_is_the_cell_contracts() {
        for token in [
            "0.0950",
            "1",
            "-0.012",
            "9.8e-05",
            "1E5",
            ".5",
            "3.",
            "+1",
            "inf",
            "-Infinity",
            "nan",
            "+.5e-1",
            "1.e5",
            "NaN",
            "INFINITY",
            "\u{663}.\u{665}",
            "1.50\n",
        ] {
            assert!(follows_the_number_grammar(token), "{token:?}");
            assert!(pyfmt::float_of(token).is_some(), "{token:?}");
        }
        for token in [
            "n/a", "1_000", " 1", "0x10", "1e", "-", "", "1,2", ".", "+", "1e+", "e5", "1.5.2",
            "infinit", "nano", "1 ", "1\n\n", "\n", "--1",
        ] {
            assert!(!follows_the_number_grammar(token), "{token:?}");
        }
        // A dotless or a dotted `i` reads as the letter to the grammar and not to `float()`.
        for token in ["\u{131}nf", "\u{130}NF", "\u{131}nf\u{131}n\u{131}ty"] {
            assert!(follows_the_number_grammar(token), "{token:?}");
            assert_eq!(pyfmt::float_of(token), None);
        }
    }

    #[test]
    fn a_canonical_cell_reprints_a_double_and_an_int32() {
        assert_eq!(canonical("MAX_AF", "0.0950"), "0.095");
        assert_eq!(canonical("MAX_AF", "1.0"), "1");
        assert_eq!(canonical("MAX_AF", "9.8e-05"), "0.000098");
        assert_eq!(canonical("MAX_AF", NULL_TEXT), NULL_TEXT);
        assert_eq!(canonical("MAX_AF", "n/a"), "n/a");
        assert_eq!(canonical("MOTIF_SCORE_CHANGE", "0.100"), "0.1");
        assert_eq!(canonical("gnomADe_AF", "0.0950,0.1,1"), "0.095,0.1,1");
        assert_eq!(canonical("Feature", "0.0950"), "0.0950");
        // An INT32 cell prints as its integer, and a text that is no INT32 as written.
        assert_eq!(canonical("TSL", "007"), "7");
        assert_eq!(canonical("TSL", "00"), "0");
        assert_eq!(canonical("STRAND", "-1"), "-1");
        assert_eq!(canonical("TSL", "2147483648"), "2147483648");
        assert_eq!(canonical("TSL", "1.5"), "1.5");
        assert_eq!(canonical("TSL", &format!("{}7", "0".repeat(5000))), "7");
        assert_eq!(canonical("TSL", &"9".repeat(5000)), "9".repeat(5000));
        assert_eq!(canonical("TSL", " 7"), " 7");
        assert_eq!(canonical("TSL", NULL_TEXT), NULL_TEXT);
        assert_eq!(canonical("Feature", "007"), "007");
        let frequency: Vec<&str> = FREQUENCY_POPULATIONS.to_vec();
        assert_eq!(frequency.len(), 28);
        for population in frequency {
            let column = match population {
                "AF" => "AF".to_string(),
                other => format!("{other}_AF"),
            };
            assert_eq!(canonical_of(&column), Canonical::DoubleList, "{column}");
        }
        for column in ["MAX_AF_POPS", "AF_AF", "_AF", "af", "Feature"] {
            assert_eq!(canonical_of(column), Canonical::Verbatim, "{column}");
        }
    }

    #[test]
    fn a_canonical_cell_keeps_what_its_type_does_not_read() {
        assert_eq!(canonical("TSL", "+5"), "5");
        assert_eq!(canonical("TSL", "-0"), "0");
        assert_eq!(canonical("TSL", "-2147483648"), "-2147483648");
        assert_eq!(canonical("TSL", "-2147483649"), "-2147483649");
        assert_eq!(canonical("TSL", ""), "");
        assert_eq!(canonical("TSL", "+"), "+");
        assert_eq!(canonical("TSL", "7\n"), "7\n");
        assert_eq!(canonical("TSL", "\u{663}"), "\u{663}");
        // A list, as a plugin's field of the name holds one, and its integers alone.
        assert_eq!(canonical("TSL", "1,5,2"), "1,5,2");
        assert_eq!(canonical("TSL", "05"), "5");
        assert_eq!(canonical("DISTANCE", "1;5;2"), "1;5;2");
        assert_eq!(canonical("STRAND", "+1,-1"), "+1,-1");
        assert_eq!(canonical("TSL", "--1"), "--1");
        assert_eq!(canonical("TSL", "+-1"), "+-1");
        assert_eq!(canonical("MAX_AF", "1.50\n"), "1.5");
        assert_eq!(canonical("MAX_AF", "-nan"), "NaN");
        assert_eq!(canonical("MAX_AF", "INF"), "inf");
        assert_eq!(canonical("MAX_AF", "-Infinity"), "-inf");
        assert_eq!(canonical("MAX_AF", "\u{663}.\u{665}"), "3.5");
        assert_eq!(canonical("MAX_AF", "1e400"), "inf");
        assert_eq!(canonical("MAX_AF", "0.5,0.25"), "0.5,0.25");
        assert_eq!(canonical("AF", "0.50,-,x,,1e-3"), "0.5,-,x,,0.001");
        assert_eq!(canonical("AF", ""), "");
        let mut out = String::new();
        assert_eq!(
            write_canonical_cell(&mut out, Canonical::DoubleList, "0.5,\u{131}nf"),
            Err("could not convert string to float: '\u{131}nf'".to_string())
        );
        // A cell past the last column of the header stays as written.
        let mut out = String::new();
        write_canonical_line(
            &mut out,
            &[Canonical::Verbatim, Canonical::Double],
            "0.50\t0.50\t0.50",
        )
        .unwrap();
        assert_eq!(out, "0.50\t0.5\t0.50");
    }

    #[test]
    fn a_positions_struct_prints_its_span_and_its_total() {
        let text = |start: Option<Member>, end: Option<Member>, total: Option<Member>| {
            let mut out = String::new();
            write_positions(&mut out, start, end, total);
            out
        };
        let whole = |value: i128| Some(Member::Whole(value));
        assert_eq!(text(whole(3), whole(4), whole(12)), "3-4/12");
        assert_eq!(text(whole(3), whole(3), whole(12)), "3/12");
        assert_eq!(text(whole(44), whole(44), None), "44");
        assert_eq!(text(None, None, whole(1200)), "-/1200");
        assert_eq!(text(whole(146), None, None), "146-?");
        assert_eq!(text(None, whole(52), None), "?-52");
        assert_eq!(text(None, None, None), NULL_TEXT);
        // Members of other types print as the interpreter prints them and agree as its
        // numbers do.
        assert_eq!(text(whole(1), Some(Member::Real(1.0)), None), "1");
        assert_eq!(text(Some(Member::Real(1.0)), whole(1), None), "1.0");
        assert_eq!(
            text(
                Some(Member::Bool(true)),
                whole(1),
                Some(Member::Bool(false))
            ),
            "True/False"
        );
        assert_eq!(
            text(
                Some(Member::Real(f64::NAN)),
                Some(Member::Real(f64::NAN)),
                None
            ),
            "nan-nan"
        );
        assert_eq!(text(Some(Member::Text("1")), whole(1), None), "1-1");
        assert_eq!(
            text(
                Some(Member::Text("x")),
                Some(Member::Text("x")),
                Some(Member::Real(1e-5))
            ),
            "x/1e-05"
        );
    }

    /// A stored name reads back as the writer's `tab_column_name` reads it: one prefix
    /// goes where a key column's name stands behind it, and every other name stays.
    #[test]
    fn a_stored_name_reads_back_as_its_tab_column() {
        for key in KEY_COLUMNS {
            assert_eq!(tab_column_name(&format!("tab_{key}")), *key);
            assert_eq!(
                tab_column_name(&format!("tab_tab_{key}")),
                format!("tab_{key}")
            );
            assert_eq!(tab_column_name(key), *key);
        }
        for name in [
            "Allele",
            "REF_ALLELE",
            "Ref",
            "ALT",
            "reference",
            "ref_",
            "tab",
            "tab_",
            "tab_tab_",
            "tab_Allele",
            "tab_Ref",
            "TAB_ref",
            "xtab_ref",
            "",
        ] {
            assert_eq!(tab_column_name(name), name);
        }
    }

    #[test]
    fn a_boolean_prints_the_spelling_of_its_column() {
        assert_eq!(boolean_spelling("CANONICAL"), ("YES", "-"));
        assert_eq!(boolean_spelling("PICK"), ("1", "-"));
        assert_eq!(boolean_spelling("GENCODE_PRIMARY"), ("1", "-"));
        assert_eq!(boolean_spelling("HIGH_INF_POS"), ("Y", "N"));
        assert_eq!(boolean_spelling("FLAGGED"), ("1", "0"));
    }

    #[test]
    fn the_adapters_parser_takes_a_path_that_opens_with_a_dash_for_an_option() {
        for argument in [
            "-x.parquet",
            "--x",
            "-x=y",
            "--x=y",
            "-a-b",
            "--cache.parquet",
            "-\u{e9}",
        ] {
            assert!(read_as_an_unknown_option(argument), "{argument}");
        }
        // A path, an argument the parser reads as a value, and a spelling of its options.
        for argument in [
            "x.parquet",
            "",
            "-",
            "-1",
            "-.5x",
            "-1.parquet",
            "-a b",
            "--",
            "--=x",
            "-=x",
            "-h",
            "-hx",
            "--he",
            "--help",
            "-i",
            "-ix",
            "--in=x",
            "-ox",
            "--o",
            "--canonical",
            "--c",
            "-i x",
        ] {
            assert!(!read_as_an_unknown_option(argument), "{argument}");
        }
    }

    #[test]
    fn the_lines_of_several_parts_are_joined_in_order_or_sorted() {
        let part = |rows: &[&str]| {
            let mut lines = Lines::new();
            for row in rows {
                lines.text.push_str(row);
                lines.end_line();
            }
            lines.sort();
            lines
        };
        let columns = ["A".to_string(), "B".to_string()];
        let parts = || {
            vec![
                part(&["b\t1", "a\n\t2", "a\t3"]),
                part(&[]),
                part(&["a\t2", "b"]),
            ]
        };
        assert_eq!(
            joined(&columns, parts(), false),
            "#A\tB\nb\t1\na\n\t2\na\t3\na\t2\nb\n"
        );
        // A line is compared without its line end, and whole where it holds one.
        assert_eq!(
            joined(&columns, parts(), true),
            "#A\tB\na\t2\na\t3\na\n\t2\nb\nb\t1\n"
        );
        assert_eq!(joined(&columns, Vec::new(), true), "#A\tB\n");
    }

    /// Plants a Parquet file of one row group whose tab columns `names` store strings, a
    /// row an item of `rows`.
    fn plant(path: &Path, names: &[&str], rows: &[&[&str]]) {
        use parquet::data_type::ByteArrayType;
        use parquet::file::properties::WriterProperties;
        use parquet::file::writer::SerializedFileWriter;
        let fields = names
            .iter()
            .map(|name| {
                let field = Type::primitive_type_builder(name, PhysicalType::BYTE_ARRAY)
                    .with_repetition(Repetition::REQUIRED)
                    .with_logical_type(Some(LogicalType::String))
                    .build();
                Arc::new(field.expect("a string field"))
            })
            .collect();
        let schema = Type::group_type_builder("schema")
            .with_fields(fields)
            .build()
            .expect("a schema");
        let file = File::create(path).expect("the file is created");
        let properties = Arc::new(WriterProperties::builder().build());
        let mut writer =
            SerializedFileWriter::new(file, Arc::new(schema), properties).expect("a file writer");
        let mut group = writer.next_row_group().expect("a row group");
        let mut at = 0;
        while let Some(mut column) = group.next_column().expect("a column") {
            let cells: Vec<ByteArray> = rows.iter().map(|row| ByteArray::from(row[at])).collect();
            column
                .typed::<ByteArrayType>()
                .write_batch(&cells, None, None)
                .expect("the column is written");
            column.close().expect("the column is closed");
            at += 1;
        }
        group.close().expect("the row group is closed");
        writer.close().expect("the file is closed");
    }

    #[test]
    fn a_parquet_row_that_opens_with_a_number_sign_is_compared_whole() {
        use crate::compare_files::compare_files_through;
        use crate::manifest::Documented;
        use crate::remainder::{FrontEnd, STREAM};
        use crate::rules::Side;
        use indexmap::IndexMap;
        use std::collections::BTreeSet;

        assert!(holds_a_header_line("#a\tb\nc\td\n"));
        assert!(holds_a_header_line("a\tb\n#c\td\n"));
        assert!(!holds_a_header_line("a#\t#b\nc\t#\n"));
        assert!(!holds_a_header_line(""));
        let directory = std::env::temp_dir().join(format!(
            "vep-formats-compare-parquet-{}-header-row",
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).expect("a directory for the planted files");
        let names = ["Uploaded_variation", "Location"];
        // The row that opens with `#` stands first in the file and second in byte order.
        let headed: [&[&str]; 4] = [
            &["#h", "c"],
            &[" rec0", "21:1000"],
            &["rec1", "21:1001"],
            &["rec2", "21:1002"],
        ];
        let rows = directory.join("rows.parquet");
        plant(&rows, &names, &headed[1..]);
        let header_row = directory.join("header-row.parquet");
        plant(&header_row, &names, &headed);
        let compared = |path: &Path, front_end: Option<&mut FrontEnd>| {
            let report = compare_files_through(
                front_end,
                "parquet",
                path,
                path,
                None,
                &mut Documented::new(),
                &BTreeSet::new(),
                &IndexMap::new(),
                &IndexMap::new(),
                4,
                false,
                None,
                None,
                None,
                &mut Side::default(),
            );
            crate::pyjson::dumps_compact(&report.expect("the output compares with itself"))
        };
        for (path, met, note) in [
            (
                &rows,
                false,
                "3 of 3 reference rows were kept back by the stream of record groups, \
                 0 of 3 record groups handed over",
            ),
            (
                &header_row,
                true,
                "the whole outputs were compared: a row that opens with # is outside the \
                 front end",
            ),
        ] {
            let text = CanonicalText::of(path);
            let mut read = String::new();
            let mut reading = Text::open(text.as_ref()).expect("the text opens");
            reading.read_to_string(&mut read).expect("the text is read");
            assert_eq!(text.header_row(), met, "{read:?}");
            assert!(!text.unread() && !text.split(), "{read:?}");
            let mut front_end = FrontEnd::new(STREAM);
            let handed = compared(path, Some(&mut front_end));
            assert_eq!(handed, compared(path, None));
            assert_eq!(front_end.note, note);
        }
        std::fs::remove_dir_all(&directory).expect("the planted files are removed");
    }
}
