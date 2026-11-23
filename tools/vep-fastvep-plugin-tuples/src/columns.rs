// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The column header of the input: its cells, the mapped columns it names with what each
//! part of them is read as, the `FV_` columns the map has no row for, and every refusal of
//! a header the run cannot place.

use crate::map::{FieldMap, Half, Rule};
use crate::pyfmt::shown;

/// The cells of one line of the input. The line feed that ends the line, and the carriage
/// return before it, belong to no cell; any other carriage return is written `%0D`, so a
/// Location, Allele or Feature cell holding one keeps its tuple on one line.
pub fn line_text(line: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    if memchr::memchr(b'\r', line).is_none() {
        return std::borrow::Cow::Borrowed(line);
    }
    let mut text = Vec::with_capacity(line.len() + 8);
    for &byte in line {
        if byte == b'\r' {
            text.extend_from_slice(b"%0D");
        } else {
            text.push(byte);
        }
    }
    std::borrow::Cow::Owned(text)
}

/// The cells of a column header line: the line without its leading `#`.
pub fn header_cells(line: &[u8]) -> Vec<Vec<u8>> {
    let leading = line.iter().take_while(|b| **b == b'#').count();
    line_text(&line[leading..])
        .split(|b| *b == b'\t')
        .map(<[u8]>::to_vec)
        .collect()
}

/// What the rules of a column read of one part.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reading {
    /// No rule reads the part.
    Unread,
    /// A rule reads all of it.
    Whole,
    /// Rules read these halves of a part written `word(score)`, and no rule all of it.
    Halves { word: bool, score: bool },
}

/// How many part readings a column holds ready; a part past them is read from the rules.
const READY: u64 = 64;

/// One mapped column of the header.
#[derive(Debug)]
pub struct Column {
    /// Its place in a row, counted from 0.
    pub index: usize,
    pub name: String,
    pub rules: Vec<Rule>,
    /// For each rule, the place of its field among the kept fields, None for a rule the
    /// run leaves out.
    pub slots: Vec<Option<usize>>,
    ready: Vec<Reading>,
}

fn reading_of(rules: &[Rule], index: u64) -> Reading {
    let mut reading = Reading::Unread;
    for rule in rules {
        if index < rule.first || index > rule.last {
            continue;
        }
        let (word, score) = match (rule.half, reading) {
            (Half::Whole, _) => return Reading::Whole,
            (_, Reading::Halves { word, score }) => (word, score),
            _ => (false, false),
        };
        reading = Reading::Halves {
            word: word || rule.half == Half::Word,
            score: score || rule.half == Half::Score,
        };
    }
    reading
}

impl Column {
    /// What the column's rules read of the part at `index`.
    pub fn reading(&self, index: u64) -> Reading {
        match self.ready.get(index as usize) {
            Some(reading) => *reading,
            None => reading_of(&self.rules, index),
        }
    }
}

/// The header's mapped columns and what the summary line says of the rest.
#[derive(Debug)]
pub struct Placed {
    pub columns: Vec<Column>,
    /// The kept fields of the mapped columns, in column order.
    pub fields: Vec<String>,
    /// The `FV_` columns the map has no row for.
    pub unmapped: Vec<Vec<u8>>,
    /// The cells of a row that reaches the last mapped column, five at the least.
    pub width: usize,
}

/// Finds the header's mapped columns. `fields_named` says the run named its fields. A
/// missing header, a header naming none of the map's columns, a mapped column named more
/// than once and the absent column of a named field are refused; two columns of one name
/// would write two cells of a row under one key.
pub fn place_columns(
    header: Option<&[Vec<u8>]>,
    map: &FieldMap,
    fields_named: bool,
) -> Result<Placed, String> {
    let Some(header) = header else {
        return Err("the input has no column header line".to_string());
    };
    let mut columns = Vec::new();
    let mut fields = Vec::new();
    for (index, name) in header.iter().enumerate() {
        let Some(rules) = map.rules(name) else {
            continue;
        };
        let slots = rules
            .iter()
            .map(|rule| {
                rule.selected.then(|| {
                    fields.push(rule.field.clone());
                    fields.len() - 1
                })
            })
            .collect();
        columns.push(Column {
            index,
            name: String::from_utf8(name.clone()).expect("a mapped name is the map's text"),
            rules: rules.to_vec(),
            slots,
            ready: (0..READY).map(|part| reading_of(rules, part)).collect(),
        });
    }
    if columns.is_empty() {
        let names: Vec<String> = header.iter().map(|name| shown(name)).collect();
        return Err(format!(
            "the field map names none of the header's columns: {}",
            names.join(", ")
        ));
    }
    let mut places: Vec<(&str, Vec<String>)> = Vec::new();
    for column in &columns {
        let place = (column.index + 1).to_string();
        match places.iter_mut().find(|(name, _)| *name == column.name) {
            Some((_, found)) => found.push(place),
            None => places.push((&column.name, vec![place])),
        }
    }
    let repeated: Vec<String> = places
        .iter()
        .filter(|(_, found)| found.len() > 1)
        .map(|(name, found)| {
            let (last, others) = found.split_last().expect("a place");
            format!("{name} (columns {} and {last})", others.join(", "))
        })
        .collect();
    if !repeated.is_empty() {
        return Err(format!(
            "the header names a mapped column more than once: {}",
            repeated.join(", ")
        ));
    }
    if fields_named {
        let mut absent: Vec<String> = map
            .columns
            .iter()
            .filter(|(column, _)| !header.iter().any(|name| name == column.as_bytes()))
            .flat_map(|(column, rules)| {
                rules
                    .iter()
                    .filter(|rule| rule.selected)
                    .map(move |rule| format!("{} ({column})", rule.field))
            })
            .collect();
        absent.sort();
        if !absent.is_empty() {
            return Err(format!(
                "the header lacks the column of a field --fields names: {}",
                absent.join(", ")
            ));
        }
    }
    let unmapped = header
        .iter()
        .filter(|name| name.starts_with(b"FV_") && map.rules(name).is_none())
        .cloned()
        .collect();
    let last = columns
        .iter()
        .map(|column| column.index)
        .max()
        .expect("a mapped column");
    Ok(Placed {
        columns,
        fields,
        unmapped,
        width: (last + 1).max(5),
    })
}

#[cfg(test)]
#[path = "columns_placing_tests.rs"]
mod placing_tests;
