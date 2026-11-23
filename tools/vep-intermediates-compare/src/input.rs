// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two inputs as streams: records under universal newlines, strict UTF-8, the trace
//! as CPython's `csv.DictReader` reads a tab-separated file, the output's 14-column rows,
//! and the child `sort` that orders each side's rows by key.

use crate::rowcodec::{push_field, push_sequence};
use memchr::{memchr, memchr2, memchr_iter};
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};

const READ_BUF: usize = 4 << 20;
const PIPE_BUF: usize = 1 << 20;

/// `csv.field_size_limit()`: the most code points one field of the trace may hold.
const FIELD_LIMIT: usize = 131_072;

/// The trace columns that key a row, then the ones the comparison reads, by header name.
const KEY_COLUMNS: [&str; 3] = ["location", "allele", "transcript_id"];
const VALUE_COLUMNS: [&str; 7] = [
    "cds_start",
    "cds_end",
    "ref_codon",
    "alt_codon",
    "ref_peptide",
    "alt_peptide",
    "consequence_set",
];

/// Why one side stopped before its rows were sorted.
#[derive(Debug, PartialEq)]
pub enum Stop {
    /// The other side failed first.
    Cancelled,
    Failed(String),
}

impl From<io::Error> for Stop {
    fn from(e: io::Error) -> Stop {
        Stop::Failed(e.to_string())
    }
}

fn failed<T>(detail: impl Into<String>) -> Result<T, Stop> {
    Err(Stop::Failed(detail.into()))
}

/// Calls `f` with each record of `reader` without its terminator and with whether it had
/// one. `\n`, `\r\n` and a lone `\r` each end a record, as Python's universal-newlines
/// text mode reads them; a final record needs no terminator.
fn for_each_record<R: Read>(
    reader: R,
    f: impl FnMut(&[u8], bool) -> Result<(), Stop>,
) -> Result<(), Stop> {
    for_each_record_buffered(reader, READ_BUF, f)
}

fn for_each_record_buffered<R: Read>(
    mut reader: R,
    capacity: usize,
    mut f: impl FnMut(&[u8], bool) -> Result<(), Stop>,
) -> Result<(), Stop> {
    let mut buf = vec![0u8; capacity.max(2)];
    let mut start = 0;
    let mut end = 0;
    let mut eof = false;
    loop {
        while let Some(i) = memchr2(b'\n', b'\r', &buf[start..end]) {
            let pos = start + i;
            if buf[pos] == b'\r' && pos + 1 == end && !eof {
                break;
            }
            f(&buf[start..pos], true)?;
            start = pos + 1;
            if buf[pos] == b'\r' && start < end && buf[start] == b'\n' {
                start += 1;
            }
        }
        if eof {
            if start < end {
                f(&buf[start..end], false)?;
            }
            return Ok(());
        }
        if start > 0 {
            buf.copy_within(start..end, 0);
            end -= start;
            start = 0;
        }
        if end == buf.len() {
            buf.resize(buf.len() * 2, 0);
        }
        match reader.read(&mut buf[end..]) {
            Ok(0) => eof = true,
            Ok(n) => end += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e.into()),
        }
    }
}

/// The fields of one record of the trace.
#[derive(Default)]
struct Fields {
    bytes: Vec<u8>,
    ends: Vec<usize>,
}

impl Fields {
    fn len(&self) -> usize {
        self.ends.len()
    }

    fn get(&self, i: usize) -> Option<&[u8]> {
        let end = *self.ends.get(i)?;
        let start = if i == 0 { 0 } else { self.ends[i - 1] };
        Some(&self.bytes[start..end])
    }

    fn clear(&mut self) {
        self.bytes.clear();
        self.ends.clear();
    }
}

#[derive(Clone, Copy, PartialEq, Default)]
enum CsvState {
    #[default]
    StartRecord,
    StartField,
    InField,
    InQuotedField,
    QuoteInQuotedField,
    EatNewline,
}

/// The reader of CPython's `csv` module under the `excel` dialect with a tab delimiter:
/// a field that starts with `"` runs to its closing quote across tabs and line ends, `""`
/// inside it is one quote, and text after the closing quote is appended as it stands.
#[derive(Default)]
struct CsvReader {
    state: CsvState,
    fields: Fields,
    field_chars: usize,
}

impl CsvReader {
    /// Feeds one line of the file. Returns whether the line ended a record, which
    /// `self.fields` then holds until the next call.
    fn feed(&mut self, line: &[u8], terminated: bool) -> Result<bool, Stop> {
        if self.state == CsvState::StartRecord {
            self.fields.clear();
            if memchr(b'"', line).is_none() {
                return self.feed_unquoted(line);
            }
        }
        for &byte in line {
            self.step(Some(byte))?;
        }
        if terminated {
            self.step(Some(b'\n'))?;
        }
        self.step(None)?;
        Ok(self.state == CsvState::StartRecord)
    }

    /// A line without a quote is its tab-separated fields; an empty line has none.
    fn feed_unquoted(&mut self, line: &[u8]) -> Result<bool, Stop> {
        if line.is_empty() {
            return Ok(true);
        }
        let mut start = 0;
        for tab in memchr_iter(b'\t', line).chain(std::iter::once(line.len())) {
            let field = &line[start..tab];
            if field.len() > FIELD_LIMIT
                && field.iter().filter(|b| is_char_start(**b)).count() > FIELD_LIMIT
            {
                return field_too_large();
            }
            self.fields.bytes.extend_from_slice(field);
            self.fields.ends.push(self.fields.bytes.len());
            start = tab + 1;
        }
        Ok(true)
    }

    /// At the end of the file: a quoted field left open ends the last record.
    fn finish(&mut self) -> bool {
        if self.state != CsvState::InQuotedField {
            return false;
        }
        self.save_field();
        self.state = CsvState::StartRecord;
        true
    }

    fn save_field(&mut self) {
        self.fields.ends.push(self.fields.bytes.len());
        self.field_chars = 0;
    }

    fn add(&mut self, byte: u8) -> Result<(), Stop> {
        if is_char_start(byte) {
            if self.field_chars >= FIELD_LIMIT {
                return field_too_large();
            }
            self.field_chars += 1;
        }
        self.fields.bytes.push(byte);
        Ok(())
    }

    /// One step of the `csv` module's state machine; `None` is the end of a line.
    fn step(&mut self, c: Option<u8>) -> Result<(), Stop> {
        use CsvState::*;
        match (self.state, c) {
            (StartRecord, None) => {}
            (StartRecord, Some(b'\n')) => self.state = EatNewline,
            (StartRecord, Some(_)) => {
                self.state = StartField;
                return self.step(c);
            }
            (StartField | InField | QuoteInQuotedField, None) => {
                self.save_field();
                self.state = StartRecord;
            }
            (StartField | InField | QuoteInQuotedField, Some(b'\n')) => {
                self.save_field();
                self.state = EatNewline;
            }
            (StartField, Some(b'"')) => self.state = InQuotedField,
            (StartField, Some(b'\t')) => self.save_field(),
            (StartField, Some(byte)) => {
                self.add(byte)?;
                self.state = InField;
            }
            (InField | QuoteInQuotedField, Some(b'\t')) => {
                self.save_field();
                self.state = StartField;
            }
            (InField, Some(byte)) => self.add(byte)?,
            (InQuotedField, None) => {}
            (InQuotedField, Some(b'"')) => self.state = QuoteInQuotedField,
            (InQuotedField, Some(byte)) => self.add(byte)?,
            (QuoteInQuotedField, Some(b'"')) => {
                self.add(b'"')?;
                self.state = InQuotedField;
            }
            (QuoteInQuotedField, Some(byte)) => {
                self.add(byte)?;
                self.state = InField;
            }
            (EatNewline, None) => self.state = StartRecord,
            (EatNewline, Some(b'\n')) => {}
            (EatNewline, Some(_)) => {
                return failed("a line end inside an unquoted field");
            }
        }
        Ok(())
    }
}

fn is_char_start(byte: u8) -> bool {
    byte & 0xc0 != 0x80
}

fn field_too_large<T>() -> Result<T, Stop> {
    failed(format!("field larger than field limit ({FIELD_LIMIT})"))
}

/// Where each wanted column stands in the trace's header: the last column of that name,
/// as a `dict` built from the header keeps it.
struct Columns {
    key: [Option<usize>; 3],
    values: [Option<usize>; 7],
}

impl Columns {
    fn of(header: &Fields) -> Columns {
        let last = |name: &str| {
            (0..header.len())
                .rev()
                .find(|&i| header.get(i) == Some(name.as_bytes()))
        };
        Columns {
            key: KEY_COLUMNS.map(last),
            values: VALUE_COLUMNS.map(last),
        }
    }
}

/// Turns the records of the trace into rows: the first record is the header, a record
/// without a field is skipped, and each other record is one row.
#[derive(Default)]
struct TraceRows {
    columns: Option<Columns>,
    row: Vec<u8>,
    sequence: u64,
}

impl TraceRows {
    fn record(&mut self, fields: &Fields, sorter: &mut Sorter) -> Result<(), Stop> {
        let Some(columns) = &self.columns else {
            self.columns = Some(Columns::of(fields));
            return Ok(());
        };
        if fields.len() == 0 {
            return Ok(());
        }
        self.row.clear();
        for (name, column) in KEY_COLUMNS.iter().zip(columns.key) {
            // A key column the header lacks is `DictReader`'s `KeyError`; one the row is
            // too short to reach is its `None`.
            let Some(column) = column else {
                return failed(format!("the header has no {name} column"));
            };
            push_field(&mut self.row, fields.get(column));
        }
        push_sequence(&mut self.row, self.sequence);
        for column in columns.values {
            // A value column the header lacks reads as the empty string.
            let value = match column {
                Some(column) => fields.get(column),
                None => Some(&b""[..]),
            };
            push_field(&mut self.row, value);
        }
        self.sequence += 1;
        sorter.push(&self.row)?;
        Ok(())
    }
}

fn feed_trace(file: File, sorter: &mut Sorter, cancel: &AtomicBool) -> Result<(), Stop> {
    let mut csv = CsvReader::default();
    let mut rows = TraceRows::default();
    let mut line_number = 0u64;
    for_each_record(file, |line, terminated| {
        line_number += 1;
        check(line, line_number, cancel)?;
        if csv.feed(line, terminated)? {
            rows.record(&csv.fields, sorter)?;
        }
        Ok(())
    })?;
    if csv.finish() {
        rows.record(&csv.fields, sorter)?;
    }
    Ok(())
}

/// The columns of the 14-column default output the comparison reads: Location, Allele and
/// Feature key a row; Consequence, CDS_position, Amino_acids and Codons are its values.
const OUTPUT_KEY: [usize; 3] = [1, 2, 4];
const OUTPUT_VALUES: [usize; 4] = [6, 8, 10, 11];
const OUTPUT_COLUMNS: usize = 14;

fn feed_output(file: File, sorter: &mut Sorter, cancel: &AtomicBool) -> Result<(), Stop> {
    let mut row = Vec::with_capacity(1024);
    let mut sequence = 0u64;
    let mut line_number = 0u64;
    for_each_record(file, |line, _| {
        line_number += 1;
        check(line, line_number, cancel)?;
        if line.first() == Some(&b'#') {
            return Ok(());
        }
        // The tabs that close the first thirteen columns; a row with fewer is skipped.
        let mut tabs = [0usize; OUTPUT_COLUMNS - 1];
        let mut found = 0;
        for tab in memchr_iter(b'\t', line).take(tabs.len()) {
            tabs[found] = tab;
            found += 1;
        }
        if found < tabs.len() {
            return Ok(());
        }
        let column = |i: usize| {
            let start = if i == 0 { 0 } else { tabs[i - 1] + 1 };
            &line[start..tabs[i]]
        };
        row.clear();
        for i in OUTPUT_KEY {
            push_field(&mut row, Some(column(i)));
        }
        push_sequence(&mut row, sequence);
        for i in OUTPUT_VALUES {
            push_field(&mut row, Some(column(i)));
        }
        sequence += 1;
        sorter.push(&row)?;
        Ok(())
    })
}

/// Every line of an input is decoded before anything is reported, so a byte sequence that
/// is not UTF-8 anywhere in the file, a comment line included, refuses the file.
fn check(line: &[u8], line_number: u64, cancel: &AtomicBool) -> Result<(), Stop> {
    if cancel.load(Ordering::Relaxed) {
        return Err(Stop::Cancelled);
    }
    if std::str::from_utf8(line).is_err() {
        return failed(format!("line {line_number} is not UTF-8"));
    }
    Ok(())
}

/// Orders rows bytewise: a child `sort -S 1G` under `LC_ALL=C` that writes `sorted`, or
/// the rows in memory on a host without `sort`.
enum Sorter {
    Child {
        child: Child,
        stdin: BufWriter<ChildStdin>,
    },
    Memory {
        rows: Vec<Vec<u8>>,
        sorted: PathBuf,
    },
}

impl Sorter {
    fn start(sorted: &Path) -> Result<Sorter, Stop> {
        let out_file =
            File::create(sorted).map_err(|e| Stop::Failed(format!("{}: {e}", sorted.display())))?;
        let spawned = Command::new("sort")
            .args(["-S", "1G"])
            .env("LC_ALL", "C")
            .env_remove("LC_COLLATE")
            .stdin(Stdio::piped())
            .stdout(out_file)
            .stderr(Stdio::piped())
            .spawn();
        match spawned {
            Ok(mut child) => {
                let stdin = child.stdin.take().expect("piped stdin");
                Ok(Sorter::Child {
                    child,
                    stdin: BufWriter::with_capacity(PIPE_BUF, stdin),
                })
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(Sorter::Memory {
                rows: Vec::new(),
                sorted: sorted.to_path_buf(),
            }),
            Err(e) => failed(format!("sort: {e}")),
        }
    }

    fn push(&mut self, row: &[u8]) -> io::Result<()> {
        match self {
            Sorter::Child { stdin, .. } => {
                stdin.write_all(row)?;
                stdin.write_all(b"\n")
            }
            Sorter::Memory { rows, .. } => {
                rows.push(row.to_vec());
                Ok(())
            }
        }
    }

    /// Ends the sort after the rows were fed, or after `fed` says why they were not.
    fn finish(self, fed: Result<(), Stop>) -> Result<(), Stop> {
        match self {
            Sorter::Child {
                mut child,
                mut stdin,
            } => {
                let fed = fed.and_then(|()| stdin.flush().map_err(Stop::from));
                drop(stdin);
                let mut stderr = Vec::new();
                if let Some(mut pipe) = child.stderr.take() {
                    let _ = pipe.read_to_end(&mut stderr);
                }
                let status = child
                    .wait()
                    .map_err(|e| Stop::Failed(format!("sort: {e}")))?;
                if !status.success() {
                    let rc = status
                        .code()
                        .map_or_else(|| "signal".to_string(), |c| c.to_string());
                    let stderr = String::from_utf8_lossy(&stderr);
                    return failed(format!("sort failed (rc={rc}): {}", stderr.trim()));
                }
                fed
            }
            Sorter::Memory { mut rows, sorted } => {
                fed?;
                rows.sort_unstable();
                let on_path = |e: io::Error| Stop::Failed(format!("{}: {e}", sorted.display()));
                let file = File::create(&sorted).map_err(on_path)?;
                let mut out = BufWriter::with_capacity(PIPE_BUF, file);
                for row in &rows {
                    out.write_all(row)
                        .and_then(|()| out.write_all(b"\n"))
                        .map_err(on_path)?;
                }
                out.flush().map_err(on_path)
            }
        }
    }
}

fn sort_side(
    input: &Path,
    sorted: &Path,
    cancel: &AtomicBool,
    feed: fn(File, &mut Sorter, &AtomicBool) -> Result<(), Stop>,
) -> Result<(), Stop> {
    let at_input = |stop: Stop| match stop {
        Stop::Failed(detail) => Stop::Failed(format!("{}: {detail}", input.display())),
        Stop::Cancelled => Stop::Cancelled,
    };
    let result = File::open(input)
        .map_err(|e| at_input(e.into()))
        .and_then(|file| {
            let mut sorter = Sorter::start(sorted)?;
            let fed = feed(file, &mut sorter, cancel).map_err(at_input);
            sorter.finish(fed)
        });
    if result.is_err() {
        cancel.store(true, Ordering::Relaxed);
    }
    result
}

/// Writes the rows of the trace at `input` to `sorted`, one per line in key order and,
/// under one key, in the order the trace holds them.
pub fn sort_trace(input: &Path, sorted: &Path, cancel: &AtomicBool) -> Result<(), Stop> {
    sort_side(input, sorted, cancel, feed_trace)
}

/// The same for the rows of the default-format output at `input`.
pub fn sort_output(input: &Path, sorted: &Path, cancel: &AtomicBool) -> Result<(), Stop> {
    sort_side(input, sorted, cancel, feed_output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rowcodec::{key_len, take_field, SEQUENCE_DIGITS};

    fn records(bytes: &[u8], capacity: usize) -> Vec<(Vec<u8>, bool)> {
        let mut out = Vec::new();
        for_each_record_buffered(bytes, capacity, |r, terminated| {
            out.push((r.to_vec(), terminated));
            Ok(())
        })
        .unwrap();
        out
    }

    #[test]
    fn records_end_at_lf_crlf_and_lone_cr() {
        let data = b"a\nb\r\nc\rd\r\n\r\ne";
        let expected: Vec<(Vec<u8>, bool)> = [
            (&b"a"[..], true),
            (b"b", true),
            (b"c", true),
            (b"d", true),
            (b"", true),
            (b"e", false),
        ]
        .iter()
        .map(|(r, t)| (r.to_vec(), *t))
        .collect();
        for capacity in [2, 3, 4, 5, 7, 4096] {
            assert_eq!(records(data, capacity), expected, "capacity {capacity}");
        }
        assert_eq!(records(b"", 4096), Vec::new());
        assert_eq!(records(b"\r", 4096), vec![(Vec::new(), true)]);
        assert_eq!(records(b"x\n", 4096), vec![(b"x".to_vec(), true)]);
    }

    /// The records `csv.reader(file, delimiter="\t")` yields for `text` opened with
    /// universal newlines.
    fn csv_records(text: &[u8]) -> Result<Vec<Vec<String>>, Stop> {
        let mut csv = CsvReader::default();
        let mut out = Vec::new();
        let mut take = |fields: &Fields| {
            out.push(
                (0..fields.len())
                    .map(|i| String::from_utf8(fields.get(i).unwrap().to_vec()).unwrap())
                    .collect(),
            );
        };
        for_each_record_buffered(text, 4096, |line, terminated| {
            if csv.feed(line, terminated)? {
                take(&csv.fields);
            }
            Ok(())
        })?;
        if csv.finish() {
            take(&csv.fields);
        }
        Ok(out)
    }

    fn rows(records: &[&[&str]]) -> Vec<Vec<String>> {
        records
            .iter()
            .map(|r| r.iter().map(|f| f.to_string()).collect())
            .collect()
    }

    #[test]
    fn unquoted_lines_split_on_tabs_and_a_blank_line_has_no_field() {
        assert_eq!(
            csv_records(b"a\tb\n\n\t\nc\t\r\n x \n").unwrap(),
            rows(&[&["a", "b"], &[], &["", ""], &["c", ""], &[" x "]])
        );
        assert_eq!(csv_records(b"").unwrap(), rows(&[]));
        assert_eq!(csv_records(b"a\tb").unwrap(), rows(&[&["a", "b"]]));
    }

    #[test]
    fn quoted_fields_span_tabs_and_line_ends_and_double_their_quotes() {
        assert_eq!(
            csv_records(b"\"a\tb\"\tc\n\"x\r\ny\rz\"\t\"\"\n\"a \"\"b\"\" c \"\"d\"\"\"\n")
                .unwrap(),
            rows(&[&["a\tb", "c"], &["x\ny\nz", ""], &["a \"b\" c \"d\""]])
        );
        // A quote inside an unquoted field is text; text after a closing quote is kept.
        assert_eq!(
            csv_records(b"a\"b\t\"c\"d\"e\n\"\"\n").unwrap(),
            rows(&[&["a\"b", "cd\"e"], &[""]])
        );
    }

    #[test]
    fn a_quoted_field_left_open_ends_at_the_end_of_the_file() {
        assert_eq!(csv_records(b"a\t\"b\nc").unwrap(), rows(&[&["a", "b\nc"]]));
        assert_eq!(
            csv_records(b"a\t\"b\nc\n").unwrap(),
            rows(&[&["a", "b\nc\n"]])
        );
        assert_eq!(csv_records(b"\"").unwrap(), rows(&[&[""]]));
    }

    #[test]
    fn a_field_holds_at_most_the_field_limit_in_code_points() {
        let at_limit = "\u{e9}".repeat(FIELD_LIMIT);
        let over_limit = "\u{e9}".repeat(FIELD_LIMIT + 1);
        let too_large = failed::<()>("field larger than field limit (131072)").unwrap_err();
        for (prefix, suffix) in [("", ""), ("\"", "\""), ("a\t", "\tb"), ("a\t\"", "\"\tb")] {
            let line = |field: &str| format!("{prefix}{field}{suffix}\n").into_bytes();
            assert!(csv_records(&line(&at_limit)).is_ok(), "{prefix:?}");
            assert_eq!(
                csv_records(&line(&over_limit)).unwrap_err(),
                too_large,
                "{prefix:?}"
            );
        }
    }

    /// The sorted rows of one side, decoded: the key fields, then the value fields.
    fn sorted_rows(
        text: &[u8],
        sort: fn(&Path, &Path, &AtomicBool) -> Result<(), Stop>,
        tag: &str,
    ) -> Result<Vec<Vec<Option<String>>>, Stop> {
        let dir = std::env::temp_dir().join(format!(
            "vep-intermediates-compare-input-{}-{tag}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let (input, sorted) = (dir.join("input"), dir.join("sorted"));
        std::fs::write(&input, text).unwrap();
        let result = sort(&input, &sorted, &AtomicBool::new(false)).map(|()| {
            let bytes = std::fs::read(&sorted).unwrap();
            bytes
                .split(|b| *b == b'\n')
                .filter(|line| !line.is_empty())
                .map(|line| {
                    let key = key_len(line);
                    let mut fields = Vec::new();
                    for mut part in [&line[..key], &line[key + SEQUENCE_DIGITS..]] {
                        while !part.is_empty() {
                            fields.push(
                                take_field(&mut part)
                                    .map(|f| String::from_utf8(f.into_owned()).unwrap()),
                            );
                        }
                    }
                    fields
                })
                .collect()
        });
        std::fs::remove_dir_all(&dir).unwrap();
        result
    }

    fn some(fields: &[&str]) -> Vec<Option<String>> {
        fields.iter().map(|f| Some(f.to_string())).collect()
    }

    #[test]
    fn trace_rows_follow_the_header_by_name() {
        let text = b"transcript_id\tnote\tallele\tlocation\tconsequence_set\tcds_start\n\
            ENST2\tx\tA\t1:200\tmissense_variant\t7\textra\n\
            \n\
            ENST1\tx\tA\t1:100\tstop_gained\n\
            ENST1\ty\tA\t1:100\n";
        let rows = sorted_rows(text, sort_trace, "by-name").unwrap();
        assert_eq!(rows.len(), 3);
        // Key order, then file order; a column the header lacks is empty and one the row
        // does not reach is absent.
        assert_eq!(
            rows[0],
            [
                some(&["1:100", "A", "ENST1"]),
                vec![None],
                some(&["", "", "", "", ""]),
                some(&["stop_gained"])
            ]
            .concat()
        );
        assert_eq!(
            rows[1],
            [
                some(&["1:100", "A", "ENST1"]),
                vec![None],
                some(&["", "", "", "", ""]),
                vec![None]
            ]
            .concat()
        );
        assert_eq!(
            rows[2],
            [
                some(&["1:200", "A", "ENST2"]),
                some(&["7", "", "", "", "", "", "missense_variant"])
            ]
            .concat()
        );
    }

    #[test]
    fn a_repeated_header_name_reads_its_last_column() {
        let text =
            b"location\tallele\ttranscript_id\tlocation\n1:1\tA\tENST1\t1:2\n1:3\tC\tENST1\n";
        let rows = sorted_rows(text, sort_trace, "repeated").unwrap();
        assert_eq!(rows[0][..3], [None, Some("C".into()), Some("ENST1".into())]);
        assert_eq!(rows[1][..3], some(&["1:2", "A", "ENST1"]));
    }

    #[test]
    fn a_trace_without_a_key_column_fails_on_its_first_row_only() {
        let header = b"location\tallele\tconsequence_set\n";
        assert_eq!(sorted_rows(header, sort_trace, "no-row").unwrap().len(), 0);
        assert_eq!(sorted_rows(b"\n\n", sort_trace, "blank").unwrap().len(), 0);
        let with_row = [&header[..], b"\n1:1\tA\tstop_gained\n"].concat();
        let Stop::Failed(detail) = sorted_rows(&with_row, sort_trace, "row").unwrap_err() else {
            panic!("cancelled");
        };
        assert!(
            detail.ends_with("the header has no transcript_id column"),
            "{detail}"
        );
        // A blank first line is a header without a column.
        assert!(sorted_rows(
            b"\nlocation\tallele\ttranscript_id\n",
            sort_trace,
            "blank-first"
        )
        .is_err());
    }

    #[test]
    fn output_rows_need_fourteen_columns_and_keep_their_bytes() {
        let text = b"## comment\n\
            #Uploaded_variation\tLocation\n\
            v1\t1:100\tA\tG\tENST1\tTranscript\tmissense_variant\t1\t2\t3\tA/T\tgCa/gTa\t-\tIMPACT=MODERATE\n\
            v2\t1:100\tA\tG\tENST1\tTranscript\t stop_gained \t1\t5-6\t3\tA/*\tgCa/tGa\t-\tX\tfifteenth\r\n\
            v3\t1:050\tA\tG\tENST1\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\n\
            v4\t1:050\t-\tG\tENST9\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\t";
        let rows = sorted_rows(text, sort_output, "output").unwrap();
        assert_eq!(
            rows,
            [
                some(&["1:050", "-", "ENST9", "intron_variant", "-", "-", "-"]),
                some(&[
                    "1:100",
                    "A",
                    "ENST1",
                    "missense_variant",
                    "2",
                    "A/T",
                    "gCa/gTa"
                ]),
                some(&[
                    "1:100",
                    "A",
                    "ENST1",
                    " stop_gained ",
                    "5-6",
                    "A/*",
                    "gCa/tGa"
                ]),
            ]
        );
    }

    #[test]
    fn a_byte_sequence_outside_utf8_refuses_either_input() {
        for (text, sort) in [
            (
                &b"location\tallele\ttranscript_id\n1:1\tA\tENST\xff\n"[..],
                sort_trace as fn(&Path, &Path, &AtomicBool) -> Result<(), Stop>,
            ),
            (&b"## comment \xc3\n"[..], sort_output),
            (&b"ok\n\xed\xa0\x80\n"[..], sort_output),
        ] {
            let Stop::Failed(detail) = sorted_rows(text, sort, "utf8").unwrap_err() else {
                panic!("cancelled");
            };
            assert!(detail.ends_with("is not UTF-8"), "{detail}");
        }
    }

    #[test]
    fn a_missing_input_fails_and_cancels_the_other_side() {
        let cancel = AtomicBool::new(false);
        let missing = Path::new("/nonexistent/vep-intermediates-compare/input");
        let sorted = std::env::temp_dir().join("vep-intermediates-compare-unused");
        assert!(matches!(
            sort_trace(missing, &sorted, &cancel),
            Err(Stop::Failed(_))
        ));
        assert!(cancel.load(Ordering::Relaxed));
        assert_eq!(
            sorted_rows_cancelled(&cancel),
            Err(Stop::Cancelled),
            "a side that reads a line after the other failed stops"
        );
    }

    fn sorted_rows_cancelled(cancel: &AtomicBool) -> Result<(), Stop> {
        let dir = std::env::temp_dir().join(format!(
            "vep-intermediates-compare-input-{}-cancelled",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let (input, sorted) = (dir.join("input"), dir.join("sorted"));
        std::fs::write(&input, b"## comment\n").unwrap();
        let result = sort_output(&input, &sorted, cancel);
        std::fs::remove_dir_all(&dir).unwrap();
        result
    }
}
