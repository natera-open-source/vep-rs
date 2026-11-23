// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The annotation metrics of one default-format output: its rows, its consequence terms
//! and the rows that carry each of twelve values.
//!
//! The file is UTF-8 text read under universal newlines, so a row ends at `\n`, `\r\n`
//! or a lone `\r`. A line that starts with `##` is skipped, a line that starts with `#`
//! names the columns, an empty line is skipped, and every other line is a row. The rows
//! after the file's leading lines are read across the threads; a line that names the
//! columns after the first row changes what the rows after it mean, so such a file is
//! read again on one thread, in order.

use crate::pipeline::{lines, Ends, Reading};
use memchr::memchr2;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};

/// Whether `str.isspace()` holds for the character.
fn is_py_space(c: char) -> bool {
    matches!(
        c,
        '\u{9}'..='\u{d}'
            | '\u{1c}'..='\u{20}'
            | '\u{85}'
            | '\u{a0}'
            | '\u{1680}'
            | '\u{2000}'..='\u{200a}'
            | '\u{2028}'
            | '\u{2029}'
            | '\u{202f}'
            | '\u{205f}'
            | '\u{3000}'
    )
}

/// `str.strip()`.
fn py_strip(text: &str) -> &str {
    text.trim_matches(is_py_space)
}

/// What `int()` reads of `text`, in ASCII digits: blanks around the number, a sign and
/// single underscores between digits, at any length; the number as `str(int)` writes it.
pub fn py_int(text: &str) -> Option<String> {
    let text = py_strip(text);
    let (negative, digits) = match text.as_bytes().first() {
        Some(b'-') => (true, &text[1..]),
        Some(b'+') => (false, &text[1..]),
        _ => (false, text),
    };
    let mut number = String::new();
    let mut after_digit = false;
    for byte in digits.bytes() {
        match byte {
            b'0'..=b'9' => {
                number.push(char::from(byte));
                after_digit = true;
            }
            b'_' if after_digit => after_digit = false,
            _ => return None,
        }
    }
    if !after_digit {
        return None;
    }
    let number = number.trim_start_matches('0');
    Some(match (number.is_empty(), negative) {
        (true, _) => "0".to_string(),
        (false, true) => format!("-{number}"),
        (false, false) => number.to_string(),
    })
}

/// A string as `json.dumps` writes it under `ensure_ascii`.
fn json_string(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            ' '..='\u{7e}' => out.push(c),
            other => {
                let mut units = [0u16; 2];
                for unit in other.encode_utf16(&mut units) {
                    out.push_str(&format!("\\u{unit:04x}"));
                }
            }
        }
    }
    out.push('"');
    out
}

/// The Extra keys whose rows are counted: three by name, then the frequency fields in
/// the order the report lists them.
const KEYS: [&str; 12] = [
    "CLIN_SIG",
    "HGVSc",
    "HGVSp",
    "AF",
    "gnomADe_AF",
    "gnomADg_AF",
    "MAX_AF",
    "EUR_AF",
    "AFR_AF",
    "AMR_AF",
    "EAS_AF",
    "SAS_AF",
];

/// The place of the column named `Extra` on a line that names the columns.
fn extra_column(header: &str) -> Option<usize> {
    header
        .trim_start_matches('#')
        .split('\t')
        .position(|name| py_strip(name) == "Extra")
}

#[derive(Default)]
pub struct Counts {
    rows: u64,
    term_index: HashMap<String, usize>,
    /// The consequence terms in the order first met, each with its count.
    terms: Vec<(String, u64)>,
    existing: u64,
    keys: [u64; 12],
}

impl Counts {
    fn term(&mut self, term: &str, count: u64) {
        match self.term_index.get(term) {
            Some(at) => self.terms[*at].1 += count,
            None => {
                self.term_index.insert(term.to_string(), self.terms.len());
                self.terms.push((term.to_string(), count));
            }
        }
    }

    /// Counts one row; `extra` is the place of the Extra column where a line named it.
    fn row(&mut self, line: &str, extra: Option<usize>) {
        self.rows += 1;
        let mut consequence = None;
        let mut existing = None;
        let mut named = None;
        let mut fourteenth = None;
        let mut columns = 0;
        for (at, cell) in line.split('\t').enumerate() {
            columns = at + 1;
            match at {
                6 => consequence = Some(cell),
                12 => existing = Some(cell),
                13 => fourteenth = Some(cell),
                _ => {}
            }
            if Some(at) == extra {
                named = Some(cell);
            }
        }
        if columns < 7 {
            return;
        }
        for term in py_strip(consequence.unwrap_or("")).split(',') {
            let term = py_strip(term);
            if !term.is_empty() {
                self.term(term, 1);
            }
        }
        let existing = py_strip(existing.unwrap_or(""));
        if !existing.is_empty() && existing != "-" {
            self.existing += 1;
        }
        let cell = named.or(fourteenth).unwrap_or("");
        if cell.is_empty() || cell == "-" {
            return;
        }
        // A key that stands twice keeps its last value.
        let mut filled = [false; 12];
        for pair in cell.split(';') {
            if let Some((key, value)) = pair.split_once('=') {
                if let Some(at) = KEYS.iter().position(|name| *name == key) {
                    filled[at] = !value.is_empty();
                }
            }
        }
        for (count, filled) in self.keys.iter_mut().zip(filled) {
            *count += u64::from(filled);
        }
    }

    /// Adds `other`, the counts of what follows in the file.
    fn merge(&mut self, other: Counts) {
        self.rows += other.rows;
        self.existing += other.existing;
        for (count, more) in self.keys.iter_mut().zip(other.keys) {
            *count += more;
        }
        for (term, count) in other.terms {
            self.term(&term, count);
        }
    }

    /// The report, as `json.dump(indent=2)` and a final newline write it; the terms the
    /// most frequent first, equal counts in the order first met.
    pub fn report(&self, elapsed: &str) -> String {
        let mut terms: Vec<&(String, u64)> = self.terms.iter().collect();
        terms.sort_by_key(|term| std::cmp::Reverse(term.1));
        let object = |pairs: Vec<(String, u64)>| -> String {
            if pairs.is_empty() {
                return "{}".to_string();
            }
            let body: Vec<String> = pairs
                .iter()
                .map(|(key, count)| format!("    {}: {count}", json_string(key)))
                .collect();
            format!("{{\n{}\n  }}", body.join(",\n"))
        };
        let terms = object(terms.into_iter().cloned().collect());
        let frequencies = object(
            KEYS[3..]
                .iter()
                .zip(&self.keys[3..])
                .map(|(key, count)| ((*key).to_string(), *count))
                .collect(),
        );
        format!(
            "{{\n  \"total_annotation_lines\": {},\n  \"consequence_counts\": {terms},\n  \"existing_variation_count\": {},\n  \"clin_sig_count\": {},\n  \"hgvsc_count\": {},\n  \"hgvsp_count\": {},\n  \"frequency_field_counts\": {frequencies},\n  \"elapsed_seconds\": {elapsed}\n}}\n",
            self.rows, self.existing, self.keys[0], self.keys[1], self.keys[2]
        )
    }

    /// The seven lines of the summary.
    pub fn summary(&self) -> String {
        let most = self.keys[3..].iter().max().copied().unwrap_or(0);
        format!(
            "  Total annotation lines: {}\n  Unique consequence types: {}\n  Existing_variation non-empty: {}\n  CLIN_SIG non-empty: {}\n  HGVSc non-empty: {}\n  HGVSp non-empty: {}\n  Max frequency field count: {most}\n",
            self.rows,
            self.terms.len(),
            self.existing,
            self.keys[0],
            self.keys[1],
            self.keys[2]
        )
    }
}

fn not_utf8() -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, "the output is not UTF-8 text")
}

/// The lines of `data`, cut at `\n` and `\r`; an empty line stands between a `\r` and
/// the `\n` after it, which no caller counts.
fn text_lines(data: &str) -> impl Iterator<Item = &str> {
    data.split(['\n', '\r'])
}

/// What a run of whole lines after the first row comes to.
struct Part {
    counts: Counts,
    /// A line that names the columns stands in the run.
    names_columns: bool,
    invalid: bool,
}

fn part(data: &[u8], extra: Option<usize>) -> Part {
    let mut made = Part {
        counts: Counts::default(),
        names_columns: false,
        invalid: false,
    };
    let Ok(text) = std::str::from_utf8(data) else {
        made.invalid = true;
        return made;
    };
    for line in text_lines(text) {
        if line.is_empty() || line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            made.names_columns = true;
            return made;
        }
        made.counts.row(line, extra);
    }
    made
}

/// Reads `file` line by line from its start and hands each line and the offset after it
/// to `each`, until `each` answers false or the file ends.
fn in_order(mut file: &File, mut each: impl FnMut(&str) -> io::Result<bool>) -> io::Result<u64> {
    let mut buf = vec![0u8; 1 << 20];
    let mut line: Vec<u8> = Vec::new();
    // The offset of the first byte of the line being read.
    let mut offset = 0u64;
    loop {
        let count = match file.read(&mut buf) {
            Ok(count) => count,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error),
        };
        if count == 0 {
            if !line.is_empty() {
                let text = std::str::from_utf8(&line).map_err(|_| not_utf8())?;
                each(text)?;
            }
            return Ok(offset);
        }
        let mut rest = &buf[..count];
        while let Some(end) = memchr2(b'\n', b'\r', rest) {
            line.extend_from_slice(&rest[..end]);
            let text = std::str::from_utf8(&line).map_err(|_| not_utf8())?;
            if !each(text)? {
                return Ok(offset);
            }
            offset += line.len() as u64 + 1;
            line.clear();
            rest = &rest[end + 1..];
        }
        line.extend_from_slice(rest);
    }
}

/// The counts of the output `file`.
pub fn counts(file: &File, threads: usize) -> io::Result<Counts> {
    // The leading lines: the place of the Extra column the last of them names, and where
    // the first row starts.
    let mut extra = None;
    let first_row = in_order(file, |line| {
        if line.is_empty() || line.starts_with("##") {
            return Ok(true);
        }
        if line.starts_with('#') {
            extra = extra_column(line);
            return Ok(true);
        }
        Ok(false)
    })?;
    let mut total = Counts::default();
    let mut names_columns = false;
    let mut invalid = false;
    lines(
        file,
        Reading::Plain { from: first_row },
        threads,
        Ends::NewlineOrReturn,
        |data| part(data, extra),
        |made: Part| {
            names_columns |= made.names_columns;
            invalid |= made.invalid;
            total.merge(made.counts);
            Ok(())
        },
    )?;
    if invalid {
        return Err(not_utf8());
    }
    if !names_columns {
        return Ok(total);
    }
    // A line names the columns after the first row: the file in order, on this thread.
    let mut total = Counts::default();
    let mut extra = None;
    let mut again = File::try_clone(file)?;
    io::Seek::rewind(&mut again)?;
    in_order(&again, |line| {
        if line.is_empty() || line.starts_with("##") {
            return Ok(true);
        }
        if line.starts_with('#') {
            extra = extra_column(line);
        } else {
            total.row(line, extra);
        }
        Ok(true)
    })?;
    Ok(total)
}

#[cfg(test)]
#[path = "metrics_tests.rs"]
mod tests;
