// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The field map: the rows that name the Ensembl VEP field of each part of a supplementary
//! column, the fields `--fields` keeps, and every refusal of a map that is malformed.

use crate::pyfmt::{
    cmp_decimal, py_repr_str, py_splitlines, saturating_u64, shown, surrogateescape,
    INT_MAX_STR_DIGITS,
};
use crate::unicode::decimal_value;
use std::cmp::Ordering;

/// The header line of a field map, column by column.
pub const MAP_HEADER: [&str; 8] = [
    "fastvep_column",
    "part",
    "vep_field",
    "source",
    "kind",
    "note",
    "list_separator",
    "reference_set",
];

/// Why a run ends without tuples.
#[derive(Debug, PartialEq, Eq)]
pub enum Failure {
    /// A refusal, reported as one `ERROR:` line.
    Refused(String),
    /// A failure no refusal covers, reported as the last line of the traceback Python
    /// prints for it.
    Uncaught(String),
    /// A thread bound that is no whole number, reported under the tool's own name.
    Threads(String),
}

/// What a rule reads of its part: all of it, or one half of a part written `word(score)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Half {
    Whole,
    Score,
    Word,
}

/// One row of the field map: the parts it reads and the field they are written under.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Rule {
    /// The first part read, counted from 1 after the dropped part; the largest `u64` for
    /// a number past it, which no entry has a part at.
    pub first: u64,
    /// The last part read, `first` for a row that reads one part.
    pub last: u64,
    pub half: Half,
    pub field: String,
    /// What joins the items of a list the column joins with `&`; empty for no list.
    pub separator: String,
    /// Whether the run keeps the field.
    pub selected: bool,
}

/// The field map: each fastVEP column with its rules, in the order the map first names it.
#[derive(Debug, PartialEq, Eq)]
pub struct FieldMap {
    pub columns: Vec<(String, Vec<Rule>)>,
}

impl FieldMap {
    /// The rules of one column, by the bytes of its header name.
    pub fn rules(&self, name: &[u8]) -> Option<&[Rule]> {
        self.columns
            .iter()
            .find(|(column, _)| column.as_bytes() == name)
            .map(|(_, rules)| rules.as_slice())
    }
}

/// The field names `--fields` lists, or None for a run without the option. An empty value,
/// and an empty name beside a comma, names no field and is refused.
pub fn kept_fields(value: Option<&[u8]>) -> Result<Option<Vec<Vec<u8>>>, Failure> {
    let Some(value) = value else {
        return Ok(None);
    };
    let names: Vec<Vec<u8>> = value.split(|b| *b == b',').map(<[u8]>::to_vec).collect();
    if names.iter().any(Vec::is_empty) {
        let beside = if value.is_empty() {
            String::new()
        } else {
            format!(" beside a comma: {}", shown(value))
        };
        return Err(Failure::Refused(format!("--fields names no field{beside}")));
    }
    Ok(Some(names))
}

/// A part as the map writes it: `N`, `N-M`, `N.word` or `N.score`, the digits of any script
/// rewritten in ASCII.
struct Part {
    first: String,
    last: Option<String>,
    half: Half,
}

/// The decimal digits `text` starts with, in ASCII, and what follows them.
fn leading_digits(text: &str) -> (String, &str) {
    let mut digits = String::new();
    for (at, c) in text.char_indices() {
        match decimal_value(c) {
            Some(value) => digits.push(char::from(b'0' + value)),
            None => return (digits, &text[at..]),
        }
    }
    (digits, "")
}

fn read_part(text: &str) -> Option<Part> {
    let (first, rest) = leading_digits(text);
    if first.is_empty() {
        return None;
    }
    let (last, half) = match rest {
        "" => (None, Half::Whole),
        ".word" => (None, Half::Word),
        ".score" => (None, Half::Score),
        _ => {
            let (last, after) = leading_digits(rest.strip_prefix('-')?);
            if last.is_empty() || !after.is_empty() {
                return None;
            }
            (Some(last), Half::Whole)
        }
    };
    Some(Part { first, last, half })
}

/// `int()` of a run of digits: Python refuses one past its digit limit with an error no
/// handler catches.
fn int_of(digits: &str) -> Result<(), Failure> {
    if digits.len() > INT_MAX_STR_DIGITS {
        return Err(Failure::Uncaught(format!(
            "ValueError: Exceeds the limit ({INT_MAX_STR_DIGITS} digits) for integer string \
             conversion: value has {} digits; use sys.set_int_max_str_digits() to increase the \
             limit",
            digits.len()
        )));
    }
    Ok(())
}

/// Reads the text of a field map. `path` is the map's name as the error lines write it,
/// `kept` the field names the run keeps, or None for all of them. A malformed row, a field
/// mapped twice and a kept name the map lacks are refused.
pub fn read_field_map(
    path: &str,
    bytes: &[u8],
    kept: Option<&[Vec<u8>]>,
) -> Result<FieldMap, Failure> {
    let refused = |message: String| Err(Failure::Refused(message));
    let Ok(text) = std::str::from_utf8(bytes) else {
        return refused(format!("the field map {path} is not UTF-8 text"));
    };
    let lines = py_splitlines(text);
    if lines
        .first()
        .is_none_or(|header| !header.split('\t').eq(MAP_HEADER))
    {
        return refused(format!(
            "{path}:1: the header is not {}",
            MAP_HEADER.join(" ")
        ));
    }
    let mut map = FieldMap {
        columns: Vec::new(),
    };
    let mut named: Vec<&str> = Vec::new();
    for (number, line) in lines.iter().enumerate().skip(1) {
        let place = format!("{path}:{}", number + 1);
        let cells: Vec<&str> = line.split('\t').collect();
        let Ok([column, part, field, _source, kind, _note, separator, reference_set]) =
            <[&str; 8]>::try_from(cells.as_slice())
        else {
            return refused(format!(
                "{place}: {} columns, not {}",
                cells.len(),
                MAP_HEADER.len()
            ));
        };
        let Some(read) = read_part(part) else {
            return refused(format!(
                "{place}: part {} is not N, N-M, N.word or N.score",
                py_repr_str(part)
            ));
        };
        int_of(&read.first)?;
        if let Some(last) = &read.last {
            int_of(last)?;
        }
        let falls = read
            .last
            .as_ref()
            .is_some_and(|last| cmp_decimal(last, &read.first) != Ordering::Greater);
        if cmp_decimal(&read.first, "1") == Ordering::Less || falls {
            return refused(format!(
                "{place}: part {} must count from 1, and a range must rise",
                py_repr_str(part)
            ));
        }
        if kind != "text" && kind != "number" {
            return refused(format!(
                "{place}: kind {} is not text or number",
                py_repr_str(kind)
            ));
        }
        if column.is_empty() || field.is_empty() || separator.is_empty() || reference_set.is_empty()
        {
            return refused(format!(
                "{place}: an empty fastvep_column, vep_field, list_separator or reference_set"
            ));
        }
        if named.contains(&field) {
            return refused(format!("{place}: {field} is mapped twice"));
        }
        named.push(field);
        let separator = if separator == "-" { "" } else { separator };
        if !separator.is_empty() && (read.half != Half::Whole || read.last.is_some()) {
            return refused(format!(
                "{place}: a list separator needs a part written N, not {}",
                py_repr_str(part)
            ));
        }
        let first = saturating_u64(&read.first);
        let rule = Rule {
            first,
            last: read.last.as_deref().map_or(first, saturating_u64),
            half: read.half,
            field: field.to_string(),
            separator: separator.to_string(),
            selected: kept.is_none_or(|kept| kept.iter().any(|name| name == field.as_bytes())),
        };
        match map.columns.iter_mut().find(|(name, _)| name == column) {
            Some((_, rules)) => rules.push(rule),
            None => map.columns.push((column.to_string(), vec![rule])),
        }
    }
    // The names are ordered as Python orders the `str` of each: by code point, a byte
    // outside UTF-8 standing as its lone surrogate.
    let mut unknown: Vec<Vec<u32>> = kept
        .unwrap_or_default()
        .iter()
        .filter(|name| {
            !named
                .iter()
                .any(|field| field.as_bytes() == name.as_slice())
        })
        .map(|name| surrogateescape(name))
        .collect();
    unknown.sort();
    unknown.dedup();
    if !unknown.is_empty() {
        let names: Vec<String> = unknown
            .iter()
            .map(|units| crate::pyfmt::backslashreplace(units))
            .collect();
        return refused(format!(
            "--fields names no row of {path}: {}",
            names.join(", ")
        ));
    }
    Ok(map)
}

#[cfg(test)]
#[path = "map_rows_tests.rs"]
mod rows_tests;
