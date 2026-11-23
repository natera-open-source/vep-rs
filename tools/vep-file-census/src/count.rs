// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Counts of the lines of a file, one count a rule, in one pass.
//!
//! The first three rules are `grep -c` on a fixed text: a line ends at a newline and at a
//! zero byte, as `grep` reads a file that holds one, and a last line without either is a
//! line. The other two are the row rules of an `awk -F'\t'` pass, where a line ends at a
//! newline alone.

use memchr::{memchr, memchr2, memchr2_iter, memmem, memrchr};

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Rule {
    /// `grep -c '^TEXT'`: the lines that start with the text.
    Starting(Vec<u8>),
    /// `grep -vc '^TEXT'`: the lines that do not start with the text.
    NotStarting(Vec<u8>),
    /// `grep -c 'TEXT'`: the lines that hold the text.
    Containing(Vec<u8>),
    /// `awk -F'\t' '!/^#/ && index(";" $8, ";NAME=") {c++} END {print c + 0}'`: the lines
    /// that do not start with `#` and whose eighth column holds the key `NAME`.
    InfoField(Vec<u8>),
    /// `awk -F'\t' '!/^#/ && index(";" $NF ";", ";NAME=") {c++} END {print c + 0}'`: the
    /// lines that do not start with `#` and whose last column holds the key `NAME`.
    ExtraKey(Vec<u8>),
}

/// Whether `cell`, a `;`-separated list, holds `key=`.
fn holds_key(cell: &[u8], key_equals: &[u8], finder: &memmem::Finder<'_>) -> bool {
    cell.starts_with(key_equals) || finder.find(cell).is_some()
}

/// The lines of `data` as `grep` counts them: ended by a newline or a zero byte.
fn grep_lines(data: &[u8]) -> u64 {
    let ended = memchr2_iter(b'\n', 0, data).count() as u64;
    let open = data.last().is_some_and(|byte| !matches!(byte, b'\n' | 0));
    ended + u64::from(open)
}

/// The lines of `data`, as `grep` cuts them, that start with `text`.
fn starting(data: &[u8], text: &[u8]) -> u64 {
    if text.is_empty() {
        return grep_lines(data);
    }
    let mut count = u64::from(data.starts_with(text));
    for end in [b'\n', 0u8] {
        let mut needle = Vec::with_capacity(text.len() + 1);
        needle.push(end);
        needle.extend_from_slice(text);
        count += memmem::find_iter(data, &needle).count() as u64;
    }
    count
}

/// The lines of `data`, as `grep` cuts them, that hold `text`.
fn containing(data: &[u8], text: &[u8]) -> u64 {
    if text.is_empty() {
        return grep_lines(data);
    }
    let finder = memmem::Finder::new(text);
    let mut count = 0;
    let mut at = 0;
    while let Some(found) = finder.find(&data[at..]) {
        count += 1;
        // The rest of the line holds no further line.
        match memchr2(b'\n', 0, &data[at + found..]) {
            Some(end) => at += found + end + 1,
            None => break,
        }
    }
    count
}

/// The rows of `data`, cut at newlines, that do not start with `#` and whose cell picked
/// by `cell` holds `key=`.
fn rows_holding(data: &[u8], key: &[u8], cell: impl Fn(&[u8]) -> &[u8]) -> u64 {
    let mut key_equals = Vec::with_capacity(key.len() + 1);
    key_equals.extend_from_slice(key);
    key_equals.push(b'=');
    let mut inner = Vec::with_capacity(key.len() + 2);
    inner.push(b';');
    inner.extend_from_slice(&key_equals);
    let finder = memmem::Finder::new(&inner);
    let mut count = 0;
    let mut rest = data;
    while !rest.is_empty() {
        let (row, after) = match memchr(b'\n', rest) {
            Some(end) => (&rest[..end], &rest[end + 1..]),
            None => (rest, &rest[rest.len()..]),
        };
        rest = after;
        if row.first() != Some(&b'#') && holds_key(cell(row), &key_equals, &finder) {
            count += 1;
        }
    }
    count
}

/// The eighth tab-separated column of `row`; empty where the row has fewer.
fn eighth_column(row: &[u8]) -> &[u8] {
    let mut rest = row;
    for _ in 0..7 {
        match memchr(b'\t', rest) {
            Some(tab) => rest = &rest[tab + 1..],
            None => return &[],
        }
    }
    &rest[..memchr(b'\t', rest).unwrap_or(rest.len())]
}

/// The last tab-separated column of `row`.
fn last_column(row: &[u8]) -> &[u8] {
    &row[memrchr(b'\t', row).map_or(0, |tab| tab + 1)..]
}

impl Rule {
    /// Whether the rule is a `grep -c`, which leaves with status 1 where it counts no
    /// line; an `awk` pass leaves with 0 whatever it counts.
    pub fn is_grep(&self) -> bool {
        matches!(
            self,
            Rule::Starting(_) | Rule::NotStarting(_) | Rule::Containing(_)
        )
    }

    /// How many lines of `data`, a run of whole lines, the rule counts.
    pub fn count(&self, data: &[u8]) -> u64 {
        match self {
            Rule::Starting(text) => starting(data, text),
            Rule::NotStarting(text) => grep_lines(data) - starting(data, text),
            Rule::Containing(text) => containing(data, text),
            Rule::InfoField(name) => rows_holding(data, name, eighth_column),
            Rule::ExtraKey(name) => rows_holding(data, name, last_column),
        }
    }
}

/// One count a rule over `data`, a run of whole lines.
pub fn counts(rules: &[Rule], data: &[u8]) -> Vec<u64> {
    rules.iter().map(|rule| rule.count(data)).collect()
}

#[cfg(test)]
#[path = "count_tests.rs"]
mod tests;
