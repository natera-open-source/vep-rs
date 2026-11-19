// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The subset of Python's `csv` module the reports use: the default dialect's reader (a
//! `"`-quoted field, `""` for a quote inside one) and its minimal-quoting writer.

use crate::input::for_each_record;
use std::fs::File;
use std::io;
use std::path::Path;

/// Splits one record into fields as `csv.reader` does under the default dialect with the
/// given delimiter: a field opening with `"` runs to its closing quote, a doubled quote
/// inside it is one quote, and text after the closing quote is kept.
pub fn split_record(record: &str, delimiter: char) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = record.chars().peekable();
    while let Some(c) = chars.next() {
        if quoted {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    chars.next();
                    field.push('"');
                } else {
                    quoted = false;
                }
            } else {
                field.push(c);
            }
        } else if c == '"' && field.is_empty() {
            quoted = true;
        } else if c == delimiter {
            fields.push(std::mem::take(&mut field));
        } else {
            field.push(c);
        }
    }
    fields.push(field);
    fields
}

/// Reads a delimited file as its header and rows, blank records skipped as `csv.DictReader`
/// skips them.
pub fn read_table(path: &Path, delimiter: u8) -> io::Result<(Vec<String>, Vec<Vec<String>>)> {
    let mut header: Option<Vec<String>> = None;
    let mut rows = Vec::new();
    for_each_record(File::open(path)?, |record| {
        if record.is_empty() {
            return Ok(());
        }
        let fields = split_record(&String::from_utf8_lossy(record), char::from(delimiter));
        if header.is_none() {
            header = Some(fields);
        } else {
            rows.push(fields);
        }
        Ok(())
    })?;
    Ok((header.unwrap_or_default(), rows))
}

/// Appends one field as `csv.writer` writes it under `QUOTE_MINIMAL`: quoted when it holds
/// the delimiter, a quote, `\r` or `\n`, a quote inside it doubled.
fn push_field(out: &mut Vec<u8>, field: &str, delimiter: char) {
    let needs_quotes = field.contains([delimiter, '"', '\r', '\n']);
    if needs_quotes {
        out.push(b'"');
        out.extend_from_slice(field.replace('"', "\"\"").as_bytes());
        out.push(b'"');
    } else {
        out.extend_from_slice(field.as_bytes());
    }
}

/// Appends one row with the given delimiter and line terminator.
pub fn push_row<S: AsRef<str>>(out: &mut Vec<u8>, fields: &[S], delimiter: char, terminator: &str) {
    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            out.extend_from_slice(delimiter.to_string().as_bytes());
        }
        push_field(out, field.as_ref(), delimiter);
    }
    out.extend_from_slice(terminator.as_bytes());
}

/// `csv.DictWriter`'s default line terminator.
pub const CRLF: &str = "\r\n";

/// `repr(list_of_str)` for an error message.
pub fn py_list_repr(items: &[&str]) -> String {
    let quoted: Vec<String> = items.iter().map(|s| format!("'{s}'")).collect();
    format!("[{}]", quoted.join(", "))
}
