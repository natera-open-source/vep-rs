// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tool's two readers read, written as JSON: the sections and rows of the index
//! page and the values of the reproduction records. Another reader of the same page or of
//! the same TOML is held to the tool's by comparing the two dumps.

use crate::page::Section;
use crate::pyjson::dumps_compact_sorted;
use crate::record::Reproduction;
use crate::value::Value;

fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

fn number(value: usize) -> Value {
    Value::Int(value.to_string())
}

pub(crate) fn section_value(section: &Section) -> Value {
    Value::Dict(vec![
        ("title".to_string(), text(&section.title)),
        ("line".to_string(), number(section.line)),
        (
            "rows".to_string(),
            Value::List(
                section
                    .rows
                    .iter()
                    .map(|row| {
                        Value::Dict(vec![
                            ("key".to_string(), text(&row.key)),
                            (
                                "cells".to_string(),
                                Value::List(row.cells.iter().map(|cell| text(cell)).collect()),
                            ),
                            ("line".to_string(), number(row.line)),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "subsections".to_string(),
            Value::List(section.subsections.iter().map(section_value).collect()),
        ),
        (
            "named".to_string(),
            Value::List(
                section
                    .named
                    .iter()
                    .map(|(key, line, linked)| {
                        Value::List(vec![text(key), number(*line), Value::Bool(*linked)])
                    })
                    .collect(),
            ),
        ),
    ])
}

/// The page's sections, one line of JSON a `##` section.
pub fn dump_page(sections: &[Section]) -> String {
    sections
        .iter()
        .map(|section| dumps_compact_sorted(&section_value(section)) + "\n")
        .collect()
}

/// A value with its type told apart where JSON alone would not: a date or a time as a
/// one-key object that names the type and holds its `repr()`.
fn typed(value: &Value) -> Value {
    match value {
        Value::List(items) => Value::List(items.iter().map(typed).collect()),
        Value::Dict(items) => Value::Dict(
            items
                .iter()
                .map(|(key, item)| (key.clone(), typed(item)))
                .collect(),
        ),
        Value::DateTime(_) | Value::Date(_) | Value::Time(_) => Value::Dict(vec![(
            format!("<{}>", value.py_type()),
            Value::Str(value.py_repr()),
        )]),
        other => other.clone(),
    }
}

/// The records' values, one line of JSON a record: its name, its keys in the order the
/// file writes them (every table's keys as a list beside its values), and its values.
pub fn dump_records(repros: &[Reproduction]) -> String {
    fn key_order(value: &Value) -> Value {
        match value {
            Value::Dict(items) => Value::List(
                items
                    .iter()
                    .map(|(key, item)| Value::List(vec![Value::Str(key.clone()), key_order(item)]))
                    .collect(),
            ),
            Value::List(items) => Value::List(items.iter().map(key_order).collect()),
            _ => Value::None,
        }
    }
    repros
        .iter()
        .map(|repro| {
            let line = Value::Dict(vec![
                ("name".to_string(), text(&repro.name)),
                ("keys".to_string(), key_order(&repro.data)),
                ("values".to_string(), typed(&repro.data)),
            ]);
            dumps_compact_sorted(&line) + "\n"
        })
        .collect()
}
