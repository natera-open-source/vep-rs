// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The records of the table as CPython's `csv.reader` reads them in the `excel` dialect
//! with a tab delimiter, and a field as `csv.writer` writes it there under `QUOTE_MINIMAL`.
//!
//! The reader takes text whose line endings are LF. A field that opens with `"` runs to
//! its closing quote across tabs and line endings, `""` inside it is one quote, text
//! behind the closing quote is appended, and a quote left open ends with the input. A `"`
//! elsewhere is a character like any other. An empty line is a record of no field.

use memchr::{memchr, memchr_iter};

/// The longest field `csv` reads, in code points (`csv.field_size_limit()`).
pub const FIELD_LIMIT: usize = 131_072;

/// The text of the error `csv` raises for a longer field.
pub const FIELD_TOO_LARGE: &str = "field larger than field limit (131072)";

/// One record: its fields in order, each a slice of `data`.
pub struct Record<'a> {
    data: &'a [u8],
    spans: &'a [(usize, usize)],
}

impl<'a> Record<'a> {
    /// Whether the record holds no field, as an empty line does.
    pub fn is_empty(&self) -> bool {
        self.spans.is_empty()
    }

    /// The field at `index`, None behind the last one.
    pub fn get(&self, index: usize) -> Option<&'a [u8]> {
        self.spans
            .get(index)
            .map(|(from, to)| &self.data[*from..*to])
    }

    pub fn fields(&self) -> impl Iterator<Item = &'a [u8]> + '_ {
        self.spans.iter().map(|(from, to)| &self.data[*from..*to])
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum State {
    StartRecord,
    StartField,
    InField,
    InQuotedField,
    QuoteInQuotedField,
    EndOfLine,
}

/// The reader. Text is fed in runs of any length; a record is handed to the sink when its
/// last line is read.
pub struct Reader {
    /// The last line of the run before, where that run ended inside it.
    line: Vec<u8>,
    state: State,
    /// The fields of a record that holds a quote, unquoted, end to end.
    data: Vec<u8>,
    spans: Vec<(usize, usize)>,
    /// Where the field being read starts in `data`, and its code points so far.
    field_start: usize,
    field_chars: usize,
}

impl Default for Reader {
    fn default() -> Reader {
        Reader::new()
    }
}

impl Reader {
    pub fn new() -> Reader {
        Reader {
            line: Vec::new(),
            state: State::StartRecord,
            data: Vec::new(),
            spans: Vec::new(),
            field_start: 0,
            field_chars: 0,
        }
    }

    /// Reads one run of text. The sink is given every record that ends in it.
    pub fn feed<E>(
        &mut self,
        run: &[u8],
        sink: &mut impl FnMut(Record<'_>) -> Result<(), E>,
    ) -> Result<(), E>
    where
        E: From<&'static str>,
    {
        let mut rest = run;
        if !self.line.is_empty() {
            let Some(end) = memchr(b'\n', rest) else {
                self.line.extend_from_slice(rest);
                return Ok(());
            };
            self.line.extend_from_slice(&rest[..=end]);
            rest = &rest[end + 1..];
            let line = std::mem::take(&mut self.line);
            self.read_line(&line, true, sink)?;
            self.line = line;
            self.line.clear();
        }
        let quoted = memchr(b'"', rest).is_some();
        while let Some(end) = memchr(b'\n', rest) {
            self.read_line(&rest[..=end], quoted, sink)?;
            rest = &rest[end + 1..];
        }
        self.line.extend_from_slice(rest);
        Ok(())
    }

    /// Reads what is left behind the last line ending and hands on the record that the end
    /// of the input closes.
    pub fn finish<E>(&mut self, sink: &mut impl FnMut(Record<'_>) -> Result<(), E>) -> Result<(), E>
    where
        E: From<&'static str>,
    {
        if !self.line.is_empty() {
            let line = std::mem::take(&mut self.line);
            self.read_line(&line, true, sink)?;
        }
        if self.state == State::InQuotedField {
            self.save_field();
            self.state = State::StartRecord;
            self.hand_on(sink)?;
        }
        Ok(())
    }

    /// Reads one line, its LF included where it has one. `may_quote` is false where the
    /// caller knows that the line holds no quote.
    fn read_line<E>(
        &mut self,
        line: &[u8],
        may_quote: bool,
        sink: &mut impl FnMut(Record<'_>) -> Result<(), E>,
    ) -> Result<(), E>
    where
        E: From<&'static str>,
    {
        if self.state == State::StartRecord && !(may_quote && memchr(b'"', line).is_some()) {
            return self.read_plain_line(line, sink);
        }
        for byte in line {
            self.read_byte(*byte)?;
        }
        // The end of the line: it closes the record unless a quoted field is open.
        match self.state {
            State::InQuotedField => return Ok(()),
            State::StartField | State::InField | State::QuoteInQuotedField => self.save_field(),
            State::StartRecord | State::EndOfLine => {}
        }
        self.state = State::StartRecord;
        self.hand_on(sink)
    }

    /// A line without a quote, read at the start of a record: its fields lie between its
    /// tabs, and an empty line holds none.
    fn read_plain_line<E>(
        &mut self,
        line: &[u8],
        sink: &mut impl FnMut(Record<'_>) -> Result<(), E>,
    ) -> Result<(), E>
    where
        E: From<&'static str>,
    {
        let text = line.strip_suffix(b"\n").unwrap_or(line);
        self.spans.clear();
        if !text.is_empty() {
            let mut from = 0;
            for tab in memchr_iter(b'\t', text).chain(std::iter::once(text.len())) {
                if tab - from > FIELD_LIMIT && char_count(&text[from..tab]) > FIELD_LIMIT {
                    return Err(FIELD_TOO_LARGE.into());
                }
                self.spans.push((from, tab));
                from = tab + 1;
            }
        }
        let result = sink(Record {
            data: text,
            spans: &self.spans,
        });
        self.spans.clear();
        result
    }

    fn read_byte(&mut self, byte: u8) -> Result<(), &'static str> {
        match (self.state, byte) {
            (State::StartRecord, b'\n') => self.state = State::EndOfLine,
            (State::StartRecord | State::StartField, b'"') => self.state = State::InQuotedField,
            (State::StartRecord | State::StartField, b'\t') => {
                self.save_field();
                self.state = State::StartField;
            }
            (State::StartField, b'\n') => {
                self.save_field();
                self.state = State::EndOfLine;
            }
            (State::StartRecord | State::StartField, _) => {
                self.add(byte)?;
                self.state = State::InField;
            }
            (State::InField | State::QuoteInQuotedField, b'\t') => {
                self.save_field();
                self.state = State::StartField;
            }
            (State::InField | State::QuoteInQuotedField, b'\n') => {
                self.save_field();
                self.state = State::EndOfLine;
            }
            (State::InField, _) => self.add(byte)?,
            (State::InQuotedField, b'"') => self.state = State::QuoteInQuotedField,
            (State::InQuotedField, _) => self.add(byte)?,
            (State::QuoteInQuotedField, b'"') => {
                self.add(byte)?;
                self.state = State::InQuotedField;
            }
            (State::QuoteInQuotedField, _) => {
                self.add(byte)?;
                self.state = State::InField;
            }
            // A line holds one LF, its last byte.
            (State::EndOfLine, _) => {}
        }
        Ok(())
    }

    /// Appends one byte to the field being read. The limit counts code points, so it is
    /// checked where a code point starts.
    fn add(&mut self, byte: u8) -> Result<(), &'static str> {
        if byte & 0xc0 != 0x80 {
            if self.field_chars >= FIELD_LIMIT {
                return Err(FIELD_TOO_LARGE);
            }
            self.field_chars += 1;
        }
        self.data.push(byte);
        Ok(())
    }

    fn save_field(&mut self) {
        self.spans.push((self.field_start, self.data.len()));
        self.field_start = self.data.len();
        self.field_chars = 0;
    }

    fn hand_on<E>(&mut self, sink: &mut impl FnMut(Record<'_>) -> Result<(), E>) -> Result<(), E> {
        let result = sink(Record {
            data: &self.data,
            spans: &self.spans,
        });
        self.data.clear();
        self.spans.clear();
        self.field_start = 0;
        result
    }
}

fn char_count(text: &[u8]) -> usize {
    crate::pyfmt::char_count(text)
}

/// Appends one field as `csv.writer` writes it: as it stands, or between quotes with every
/// quote doubled where it holds a tab, a quote, an LF or a CR.
pub fn push_field(out: &mut Vec<u8>, field: &[u8]) {
    if !field
        .iter()
        .any(|byte| matches!(byte, b'\t' | b'"' | b'\n' | b'\r'))
    {
        out.extend_from_slice(field);
        return;
    }
    out.push(b'"');
    for byte in field {
        if *byte == b'"' {
            out.push(b'"');
        }
        out.push(*byte);
    }
    out.push(b'"');
}

#[cfg(test)]
#[path = "csv_tests.rs"]
mod tests;
