// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reproduction records: the `<entry>.toml` files of one directory, each read as TOML
//! into values, their keys in the order the file writes them.

use crate::pyfmt::decode_utf8_strict;
use crate::pypath;
use crate::value::{Date, DateTime, Time, Value};

/// One reproduction: the file's stem, its path and its values.
#[derive(Clone, Debug)]
pub struct Reproduction {
    pub name: String,
    pub path: String,
    pub data: Value,
}

impl Reproduction {
    pub fn input(&self) -> &Value {
        self.data.section("input")
    }

    pub fn cache(&self) -> &Value {
        self.data.section("cache")
    }

    pub fn run(&self) -> &Value {
        self.data.section("run")
    }

    pub fn compare(&self) -> &Value {
        self.data.section("compare")
    }

    pub fn expected(&self) -> &Value {
        self.data.section("expected")
    }

    /// The runs the record holds, its own first: the record without its `also` list, then
    /// one record for each table of that list, which is the record's own with the tables of
    /// `RUN_TABLES` that the table gives in place of the record's. Every run reads the
    /// record's cache and its `prepare` steps.
    pub fn runs(&self) -> Vec<Reproduction> {
        let own: Vec<(String, Value)> = self
            .data
            .items()
            .iter()
            .filter(|(key, _)| key != "also")
            .cloned()
            .collect();
        let run_of = |data: Vec<(String, Value)>| Reproduction {
            name: self.name.clone(),
            path: self.path.clone(),
            data: Value::Dict(data),
        };
        let mut runs = vec![run_of(own.clone())];
        if let Some(Value::List(further)) = self.data.get("also") {
            for table in further {
                let mut data = own.clone();
                for key in RUN_TABLES {
                    let Some(value) = table.get(key) else {
                        continue;
                    };
                    match data.iter_mut().find(|(name, _)| name == key) {
                        Some(slot) => slot.1 = value.clone(),
                        None => data.push((key.to_string(), value.clone())),
                    }
                }
                runs.push(run_of(data));
            }
        }
        runs
    }
}

/// The tables a further run of a record's `also` list gives in place of the record's own.
pub const RUN_TABLES: [&str; 4] = ["input", "run", "compare", "expected"];

/// Why the records of a directory could not be read.
#[derive(Debug)]
pub enum LoadError {
    /// A record that is missing or is not TOML.
    Reproduction(String),
    /// A file the system would not open or read, in `str(OSError)` form.
    Os(String),
}

impl LoadError {
    pub fn message(&self) -> &str {
        match self {
            LoadError::Reproduction(message) | LoadError::Os(message) => message,
        }
    }
}

fn convert_date(date: &toml_edit::Date) -> Date {
    Date {
        year: date.year,
        month: date.month,
        day: date.day,
    }
}

fn convert_time(time: &toml_edit::Time) -> Time {
    Time {
        hour: time.hour,
        minute: time.minute,
        second: time.second,
        microsecond: time.nanosecond / 1000,
    }
}

fn convert_value(value: &toml_edit::Value) -> Value {
    match value {
        toml_edit::Value::String(text) => Value::Str(text.value().clone()),
        toml_edit::Value::Integer(number) => Value::Int(number.value().to_string()),
        toml_edit::Value::Float(number) => Value::Float(*number.value()),
        toml_edit::Value::Boolean(flag) => Value::Bool(*flag.value()),
        toml_edit::Value::Datetime(datetime) => {
            let datetime = datetime.value();
            match (&datetime.date, &datetime.time) {
                (Some(date), Some(time)) => Value::DateTime(DateTime {
                    date: convert_date(date),
                    time: convert_time(time),
                    offset_minutes: datetime.offset.map(|offset| match offset {
                        toml_edit::Offset::Z => 0,
                        toml_edit::Offset::Custom { minutes } => i32::from(minutes),
                    }),
                }),
                (Some(date), None) => Value::Date(convert_date(date)),
                (None, Some(time)) => Value::Time(convert_time(time)),
                (None, None) => Value::None,
            }
        }
        toml_edit::Value::Array(items) => Value::List(items.iter().map(convert_value).collect()),
        toml_edit::Value::InlineTable(table) => Value::Dict(
            table
                .iter()
                .map(|(key, value)| (key.to_string(), convert_value(value)))
                .collect(),
        ),
    }
}

/// The first header, in the document's order, that names a table or a table under it.
fn first_header(table: &toml_edit::Table) -> Option<usize> {
    table
        .iter()
        .filter_map(|(_, item)| item_first_header(item))
        .chain(table.position())
        .min()
}

fn item_first_header(item: &toml_edit::Item) -> Option<usize> {
    match item {
        toml_edit::Item::Table(table) => first_header(table),
        toml_edit::Item::ArrayOfTables(tables) => tables.iter().filter_map(first_header).min(),
        _ => None,
    }
}

/// A table's keys in the order the document first names them: a table that a header of
/// one of its sub-tables names before its own header stands where that header stands. The
/// keys written under a header stand at that header, `block`, in their own order.
fn convert_table(table: &toml_edit::Table, block: Option<usize>) -> Value {
    let block = table.position().or(block);
    let mut items: Vec<((Option<usize>, usize), String, Value)> = table
        .iter()
        .enumerate()
        .map(|(index, (key, item))| {
            let (first, value) = match item {
                toml_edit::Item::Table(inner) if !inner.is_dotted() => {
                    (first_header(inner), convert_table(inner, block))
                }
                toml_edit::Item::Table(inner) => (None, convert_table(inner, block)),
                toml_edit::Item::ArrayOfTables(tables) => (
                    item_first_header(item),
                    Value::List(
                        tables
                            .iter()
                            .map(|inner| convert_table(inner, block))
                            .collect(),
                    ),
                ),
                toml_edit::Item::Value(value) => (None, convert_value(value)),
                toml_edit::Item::None => (None, Value::None),
            };
            let place = match first {
                Some(header) => (Some(header), 0),
                None => (block, index),
            };
            (place, key.to_string(), value)
        })
        .collect();
    items.sort_by_key(|(place, _, _)| *place);
    Value::Dict(
        items
            .into_iter()
            .map(|(_, key, value)| (key, value))
            .collect(),
    )
}

/// The line and column, counted from 1 in characters, of a byte offset of a text.
fn line_and_column(text: &str, offset: usize) -> (usize, usize) {
    let before = &text[..offset.min(text.len())];
    let line = before.matches('\n').count() + 1;
    let column = before
        .rsplit('\n')
        .next()
        .map_or(0, |tail| tail.chars().count())
        + 1;
    (line, column)
}

/// A TOML document as values; the error is one line that says what the reader refused and
/// where.
pub fn parse_toml(bytes: &[u8]) -> Result<Value, String> {
    let text = decode_utf8_strict(bytes)?;
    // A record written with `\r\n` line ends holds the same values as one written with `\n`.
    let text = text.replace("\r\n", "\n");
    if let Some(offset) = text.find('\r') {
        // What is left is no line end of TOML, wherever it stands.
        let (line, column) = line_and_column(&text, offset);
        return Err(format!(
            "a carriage return that no line feed follows (at line {line}, column {column})"
        ));
    }
    match toml_edit::ImDocument::parse(text.as_str()) {
        Ok(document) => Ok(convert_table(document.as_table(), None)),
        Err(error) => {
            let message: Vec<&str> = error.message().lines().map(str::trim).collect();
            let place = match error.span() {
                Some(span) if span.start < text.len() => {
                    let (line, column) = line_and_column(&text, span.start);
                    format!("at line {line}, column {column}")
                }
                _ => "at end of document".to_string(),
            };
            Err(format!("{} ({place})", message.join(", ")))
        }
    }
}

/// Every record of a directory, or the named ones in the order named. A name without a
/// file, a directory without a record and a file that is not TOML are errors.
pub fn load_reproductions(
    names: Option<&[String]>,
    directory: &str,
) -> Result<Vec<Reproduction>, LoadError> {
    let mut found: Vec<(String, String)> = Vec::new();
    if let Ok(entries) = std::fs::read_dir(directory) {
        let mut file_names: Vec<String> = entries
            .filter_map(Result::ok)
            .map(|entry| entry.file_name().to_string_lossy().into_owned())
            .filter(|name| name.ends_with(".toml"))
            .collect();
        file_names.sort();
        for file_name in file_names {
            found.push((
                pypath::stem(&file_name).to_string(),
                pypath::join(directory, &file_name),
            ));
        }
    }
    let path_of = |name: &str| found.iter().find(|(stem, _)| stem == name).map(|(_, p)| p);
    let selected: Vec<(String, String)> = match names.filter(|names| !names.is_empty()) {
        Some(names) => {
            let missing: Vec<&str> = names
                .iter()
                .map(String::as_str)
                .filter(|name| path_of(name).is_none())
                .collect();
            if !missing.is_empty() {
                return Err(LoadError::Reproduction(format!(
                    "no reproduction for {} under {directory}",
                    missing.join(", ")
                )));
            }
            let mut selected: Vec<(String, String)> = Vec::new();
            for name in names {
                if !selected.iter().any(|(stem, _)| stem == name) {
                    selected.push((name.clone(), path_of(name).expect("a found name").clone()));
                }
            }
            selected
        }
        None if found.is_empty() => {
            return Err(LoadError::Reproduction(format!(
                "no reproduction under {directory}"
            )))
        }
        None => found,
    };
    let mut reproductions = Vec::with_capacity(selected.len());
    for (name, path) in selected {
        let bytes =
            std::fs::read(&path).map_err(|error| LoadError::Os(pypath::os_error(&error, &path)))?;
        let data = parse_toml(&bytes)
            .map_err(|error| LoadError::Reproduction(format!("{path}: {error}")))?;
        reproductions.push(Reproduction { name, path, data });
    }
    Ok(reproductions)
}

#[cfg(test)]
#[path = "record_tests.rs"]
mod tests;
