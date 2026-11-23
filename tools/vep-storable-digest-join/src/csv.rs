// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A table as CPython's `csv` module reads and writes it in the `excel` dialect: the
//! records of `csv.reader` over a file opened with `newline=""`, the rows of
//! `csv.DictReader`, and a record of `csv.writer` under `QUOTE_MINIMAL` with
//! `lineterminator="\n"`.

/// The longest field `csv` reads, in characters (`csv.field_size_limit()`).
const FIELD_LIMIT: usize = 131_072;

#[derive(Clone, Copy, PartialEq, Eq)]
enum State {
    StartRecord,
    StartField,
    InField,
    InQuotedField,
    QuoteInQuotedField,
    EatNewline,
}

/// The records of a text. A field that opens with `"` runs to its closing quote across
/// commas and line ends, `""` inside it is one quote and text behind the closing quote
/// is appended; a `"` elsewhere is a character like any other; a quote left open ends
/// with the text. A line ends at `\n`, `\r` or `\r\n`, and an empty line is a record of
/// no field. The error is `csv`'s message.
pub fn records(text: &str) -> Result<Vec<Vec<String>>, String> {
    let mut all = Vec::new();
    let mut fields: Vec<String> = Vec::new();
    let mut field = String::new();
    let mut length = 0usize;
    let mut state = State::StartRecord;
    // Whether a record is under way: a line was read since the last record ended.
    let mut open = false;
    let add = |field: &mut String, length: &mut usize, c: char| -> Result<(), String> {
        if *length >= FIELD_LIMIT {
            return Err(format!("field larger than field limit ({FIELD_LIMIT})"));
        }
        field.push(c);
        *length += 1;
        Ok(())
    };
    let save = |fields: &mut Vec<String>, field: &mut String, length: &mut usize| {
        fields.push(std::mem::take(field));
        *length = 0;
    };
    // `None` is the end of a line: the reader is told so after each line's characters.
    let step = |c: Option<char>,
                state: &mut State,
                fields: &mut Vec<String>,
                field: &mut String,
                length: &mut usize|
     -> Result<(), String> {
        let newline = matches!(c, Some('\n' | '\r'));
        if *state == State::StartRecord {
            match c {
                None => return Ok(()),
                Some(_) if newline => {
                    *state = State::EatNewline;
                    return Ok(());
                }
                Some(_) => *state = State::StartField,
            }
        }
        match *state {
            State::StartRecord => {}
            State::StartField => match c {
                None => {
                    save(fields, field, length);
                    *state = State::StartRecord;
                }
                Some(_) if newline => {
                    save(fields, field, length);
                    *state = State::EatNewline;
                }
                Some('"') => *state = State::InQuotedField,
                Some(',') => save(fields, field, length),
                Some(c) => {
                    add(field, length, c)?;
                    *state = State::InField;
                }
            },
            State::InField => match c {
                None => {
                    save(fields, field, length);
                    *state = State::StartRecord;
                }
                Some(_) if newline => {
                    save(fields, field, length);
                    *state = State::EatNewline;
                }
                Some(',') => {
                    save(fields, field, length);
                    *state = State::StartField;
                }
                Some(c) => add(field, length, c)?,
            },
            State::InQuotedField => match c {
                None => {}
                Some('"') => *state = State::QuoteInQuotedField,
                Some(c) => add(field, length, c)?,
            },
            State::QuoteInQuotedField => match c {
                Some('"') => {
                    add(field, length, '"')?;
                    *state = State::InQuotedField;
                }
                Some(',') => {
                    save(fields, field, length);
                    *state = State::StartField;
                }
                None => {
                    save(fields, field, length);
                    *state = State::StartRecord;
                }
                Some(_) if newline => {
                    save(fields, field, length);
                    *state = State::EatNewline;
                }
                Some(c) => {
                    add(field, length, c)?;
                    *state = State::InField;
                }
            },
            State::EatNewline => match c {
                None => *state = State::StartRecord,
                Some(_) if newline => {}
                Some(_) => {
                    return Err("new-line character seen in unquoted field - \
                                do you need to open the file with newline=''?"
                        .to_string())
                }
            },
        }
        Ok(())
    };
    for line in lines_keeping_ends(text) {
        open = true;
        for c in line.chars() {
            step(Some(c), &mut state, &mut fields, &mut field, &mut length)?;
        }
        step(None, &mut state, &mut fields, &mut field, &mut length)?;
        if state == State::StartRecord {
            all.push(std::mem::take(&mut fields));
            open = false;
        }
    }
    // The text ended inside a quoted field: what was read of it is the last field.
    if open && (length != 0 || state == State::InQuotedField) {
        fields.push(field);
        all.push(fields);
    }
    Ok(all)
}

/// The lines of a file opened with `newline=""`: each ends after its `\n`, `\r` or
/// `\r\n`, which it keeps.
fn lines_keeping_ends(text: &str) -> Vec<&str> {
    let bytes = text.as_bytes();
    let mut lines = Vec::new();
    let mut start = 0;
    let mut at = 0;
    while at < bytes.len() {
        let end = match bytes[at] {
            b'\n' => at + 1,
            b'\r' if bytes.get(at + 1) == Some(&b'\n') => at + 2,
            b'\r' => at + 1,
            _ => {
                at += 1;
                continue;
            }
        };
        lines.push(&text[start..end]);
        start = end;
        at = end;
    }
    if start < bytes.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// A table has no column of the name asked for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NoColumn;

/// A table as `csv.DictReader` gives it: the names of the first record, and every later
/// record that is not empty.
pub struct Table {
    pub header: Vec<String>,
    /// Each row holds one value a name of the header: `None` where the record is short.
    /// A record's values past the header are not kept.
    pub rows: Vec<Vec<Option<String>>>,
}

impl Table {
    pub fn read(text: &str) -> Result<Table, String> {
        let mut records = records(text)?.into_iter();
        let header = records.next().unwrap_or_default();
        let rows = records
            .filter(|record| !record.is_empty())
            .map(|record| {
                let mut values: Vec<Option<String>> = record.into_iter().map(Some).collect();
                values.resize(header.len(), None);
                values
            })
            .collect();
        Ok(Table { header, rows })
    }

    /// The value a row's dict holds under `name`: the last column of that name.
    /// [`NoColumn`] is the `KeyError` of a table without such a column.
    pub fn get<'a>(
        &self,
        row: &'a [Option<String>],
        name: &str,
    ) -> Result<Option<&'a str>, NoColumn> {
        self.header
            .iter()
            .rposition(|column| column == name)
            .map(|index| row[index].as_deref())
            .ok_or(NoColumn)
    }
}

/// One record as `csv.writer` writes it: a field in quotes when it holds a comma, a
/// quote, `\n` or `\r`, with its quotes doubled; `None` as an empty field; a record of
/// one empty field as `""`; `\n` after the record.
pub fn write_record(out: &mut String, fields: &[Option<&str>]) {
    let start = out.len();
    for (index, field) in fields.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let field = field.unwrap_or("");
        if field.contains([',', '"', '\n', '\r']) {
            out.push('"');
            for c in field.chars() {
                if c == '"' {
                    out.push('"');
                }
                out.push(c);
            }
            out.push('"');
        } else {
            out.push_str(field);
        }
    }
    if fields.len() == 1 && out.len() == start {
        out.push_str("\"\"");
    }
    out.push('\n');
}

#[cfg(test)]
#[path = "csv_tests.rs"]
mod tests;
