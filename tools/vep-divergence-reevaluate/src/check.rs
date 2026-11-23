// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--check`: the index page held against the reproduction records, its `cache
//! transcripts` cells to one total and each record's `implementation` table to the
//! repository's files, with nothing run.

use crate::item::items;
use crate::page::{
    self, check_classes, leaves, parse_page, PagePatterns, Row, COUNT_CELL, VERIFIED_CELL,
};
use crate::pyfmt::{decode_utf8_strict, py_int_grouped, py_repr_str, universal_newlines};
use crate::pypath;
use crate::record::{load_reproductions, LoadError, Reproduction};
use crate::schema::{implementation_row, Schema, IMPLEMENTATION_LISTS, SIGHTING_FORM};
use crate::value::Value;
use std::collections::HashMap;

/// The input class of each record that names one of the three, in the records' order.
pub fn input_classes(repros: &[Reproduction]) -> Vec<(String, String)> {
    repros
        .iter()
        .filter_map(|repro| {
            let input_class = repro.data.get_str("input_class")?;
            page::INPUT_CLASSES
                .contains(&input_class)
                .then(|| (repro.name.clone(), input_class.to_string()))
        })
        .collect()
}

/// The unit of a figure counted over the transcripts of the committed caches, one quantity
/// with one total, as `scripts/divergences/check_divergence_order.py` writes it in a
/// `Tier and count` cell.
const CACHE_TRANSCRIPTS: &str = "cache transcripts";

/// One line when the `Tier and count` cells whose figure reads `<n> of <total> cache
/// transcripts` name more than one total: each total with the number of rows that name it
/// and the first of those entries, in page order.
fn cache_transcript_total_errors(rows: &[(&str, &Row)], patterns: &PagePatterns) -> Vec<String> {
    let unit = format!(" {CACHE_TRANSCRIPTS}");
    let mut named: Vec<(String, Vec<&str>)> = Vec::new();
    for &(name, row) in rows {
        let figure = patterns
            .cell
            .match_start(&row.cells[COUNT_CELL])
            .and_then(|m| m.name("figure"))
            .unwrap_or_default();
        let total = figure
            .strip_suffix(&unit)
            .and_then(|counted| counted.rsplit_once(" of "))
            .map(|(_, total)| total);
        let Some(total) = total else {
            continue;
        };
        match named.iter_mut().find(|(known, _)| known == total) {
            Some((_, names)) => names.push(name),
            None => named.push((total.to_string(), vec![name])),
        }
    }
    if named.len() < 2 {
        return Vec::new();
    }
    let totals: Vec<String> = named
        .iter()
        .map(|(total, names)| {
            let rows = if names.len() == 1 { "row" } else { "rows" };
            format!(
                "{total} in {} {rows}, the first {}",
                py_int_grouped(names.len()),
                names[0]
            )
        })
        .collect();
    vec![format!(
        "the page's `Tier and count` cells name {} totals of {CACHE_TRANSCRIPTS}: {}",
        named.len(),
        totals.join("; ")
    )]
}

/// What a file of the repository defines, as an `implementation` row names it.
enum Source {
    /// A Rust file's items by name.
    Rust(Vec<String>),
    /// The text of a Python or a shell file.
    Script(String),
    /// A file that is missing, is not UTF-8 or is neither Rust, Python nor shell, with the
    /// words that say so.
    Refused(String),
}

/// Whether a Python or a shell text defines a function or a class of the name: a line
/// that opens `def <name>(`, `class <name>` or `<name>()`.
fn script_defines(text: &str, name: &str) -> bool {
    text.lines().any(|line| {
        let line = line.trim_start();
        let line = line.strip_prefix("async ").unwrap_or(line);
        let after = |opening: &str| {
            line.strip_prefix(opening)
                .and_then(|rest| rest.strip_prefix(name))
        };
        after("def ").is_some_and(|rest| rest.trim_start().starts_with('('))
            || after("class ").is_some_and(|rest| rest.starts_with(['(', ':']))
            || after("function ")
                .is_some_and(|rest| !rest.starts_with(|ch: char| ch.is_alphanumeric() || ch == '_'))
            || after("").is_some_and(|rest| rest.trim_start().starts_with("()"))
    })
}

fn read_source(repo: &str, path: &str) -> Source {
    let full = pypath::join(repo, path);
    let rust = path.ends_with(".rs");
    if !rust && !path.ends_with(".py") && !path.ends_with(".sh") {
        return Source::Refused(format!(
            "{path} is neither a Rust, a Python nor a shell file"
        ));
    }
    if !pypath::is_file(&full) {
        return Source::Refused(format!("the repository holds no file {path}"));
    }
    let text = match std::fs::read(&full).map(String::from_utf8) {
        Ok(Ok(text)) => text,
        _ => return Source::Refused(format!("{path} cannot be read as UTF-8")),
    };
    if rust {
        Source::Rust(items(&text).into_iter().map(|(name, _, _)| name).collect())
    } else {
        Source::Script(text)
    }
}

/// What is wrong with one record's `implementation` table against the repository: no
/// table, a path the repository lacks and an item its file does not define. A row whose
/// shape the schema refuses is not read here as well.
fn implementation_errors(
    repro: &Reproduction,
    repo: &str,
    sources: &mut HashMap<String, Source>,
) -> Vec<String> {
    let Some(table) = repro.data.get("implementation") else {
        return vec!["the record has no implementation table".to_string()];
    };
    let mut errors = Vec::new();
    for key in IMPLEMENTATION_LISTS {
        let Some(Value::List(rows)) = table.get(key) else {
            continue;
        };
        for row in rows.iter().filter_map(Value::as_str) {
            let Some((path, item)) = implementation_row(row) else {
                continue;
            };
            let source = sources
                .entry(path.to_string())
                .or_insert_with(|| read_source(repo, path));
            let refusal = match source {
                Source::Refused(why) => Some(why.clone()),
                Source::Rust(names) => (!names.iter().any(|name| name == item))
                    .then(|| format!("{path} defines no item {item}")),
                Source::Script(text) => (!script_defines(text, item))
                    .then(|| format!("{path} defines no function or class {item}")),
            };
            if let Some(why) = refusal {
                errors.push(format!(
                    "implementation.{key} names {}: {why}",
                    py_repr_str(row)
                ));
            }
        }
    }
    errors
}

/// Every entry without a record, every record without an entry, every schema violation of
/// an entry's record, every entry whose record has no `implementation` table or names in
/// it a path the repository at `repo` lacks or an item the file does not define, every
/// version that is not the reference, every `observed` basis the page's cell contradicts,
/// every input class the page's tiers and sections contradict, and the page's `cache
/// transcripts` cells when they name more than one total, as one line each. The outer
/// error is what stops the check: a page or a record that cannot be read.
pub fn check(
    page_path: &str,
    reference: &str,
    directory: &str,
    repo: &str,
    schema: &Schema,
) -> Result<Vec<String>, String> {
    let mut errors = Vec::new();
    let bytes = std::fs::read(page_path).map_err(|error| pypath::os_error(&error, page_path))?;
    let text = universal_newlines(decode_utf8_strict(&bytes)?);
    let sections = parse_page(&text, &schema.page)?;
    // A key the page lists twice keeps its first place and its last row.
    let mut rows: Vec<(&str, &Row)> = Vec::new();
    for section in leaves(&sections) {
        for row in &section.rows {
            match rows.iter_mut().find(|(key, _)| *key == row.key) {
                Some(entry) => entry.1 = row,
                None => rows.push((&row.key, row)),
            }
        }
    }
    let repros = load_reproductions(None, directory).map_err(|error| match error {
        LoadError::Reproduction(message) | LoadError::Os(message) => message,
    })?;
    let mut sources: HashMap<String, Source> = HashMap::new();
    for (name, row) in &rows {
        let Some(repro) = repros.iter().find(|repro| repro.name == *name) else {
            errors.push(format!("{name}: no reproduction under {directory}"));
            continue;
        };
        for problem in schema.validate(repro) {
            errors.push(format!("{name}: {problem}"));
        }
        for problem in implementation_errors(repro, repo, &mut sources) {
            errors.push(format!("{name}: {problem}"));
        }
        let own = repro.data.get("verified_against");
        if own.and_then(Value::as_str) != Some(reference) {
            errors.push(format!(
                "{name}: verified against {}, the reference is {reference}",
                own.map_or_else(|| "None".to_string(), Value::py_str)
            ));
        }
        let verified = &row.cells[VERIFIED_CELL];
        if verified != reference {
            errors.push(format!(
                "{name}: the page says verified against {}, the reference is {reference}",
                py_repr_str(verified)
            ));
        }
        let page_cell = &row.cells[COUNT_CELL];
        let observed = repro.data.get_str("input_class") == Some("observed");
        // A basis the schema refuses is not held to the page as well.
        if observed && schema.input_class_errors(&repro.data).is_empty() {
            let basis = repro
                .data
                .get_str("input_class_basis")
                .expect("a basis the schema passed");
            let tier = schema
                .page
                .cell
                .match_start(page_cell)
                .and_then(|m| m.name("tier"));
            if tier.as_deref() == Some("3") {
                // At tier 3 ordinary data does not carry the shape, so the cell is no
                // sighting and the basis names the real data the entry's page records.
                if schema.sighting.match_start(basis).is_none() {
                    errors.push(format!(
                        "{name}: input_class_basis is {}: an observed entry of tier 3 reads `{SIGHTING_FORM}`",
                        py_repr_str(basis)
                    ));
                }
            } else if basis != page_cell {
                errors.push(format!(
                    "{name}: input_class_basis is {}, the page's `Tier and count` cell is {}",
                    py_repr_str(basis),
                    py_repr_str(page_cell)
                ));
            }
        }
    }
    for repro in &repros {
        if !rows.iter().any(|(key, _)| *key == repro.name) {
            errors.push(format!(
                "{}: a reproduction with no row on the page",
                repro.name
            ));
        }
    }
    errors.extend(check_classes(
        &sections,
        &input_classes(&repros),
        &schema.page,
    ));
    errors.extend(cache_transcript_total_errors(&rows, &schema.page));
    Ok(errors)
}

#[cfg(test)]
#[path = "check_tests.rs"]
mod tests;
