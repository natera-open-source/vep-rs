// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `pairs.json`, the counts the scorer takes in place of its own extraction and merge.
//!
//! Field names and the column header are hexadecimal, since an input can carry bytes that
//! no text holds; every object lists its keys in ascending order.

use crate::extract::Extracted;
use crate::input::Kind;
use crate::order::hex;
use crate::pair::{Counts, PROJECTIONS};
use crate::resolve::Table;
use std::fmt::Write;

pub struct InputFacts {
    pub kind: Kind,
    pub sha256: String,
    pub extracted: Extracted,
}

pub struct EngineFacts {
    pub table: Table,
    pub counts: [Counts; 2],
    /// The rows written for the reference and for the engine.
    pub rows: [u64; 2],
}

/// A JSON string of a source's name, which the driver states as text.
fn quoted(text: &str) -> String {
    let mut out = String::from("\"");
    for character in text.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            control if control < ' ' => write!(out, "\\u{:04x}", control as u32).unwrap(),
            other => out.push(other),
        }
    }
    out.push('"');
    out
}

fn object(entries: Vec<(String, String)>) -> String {
    let mut entries = entries;
    entries.sort();
    let body: Vec<String> = entries
        .into_iter()
        .map(|(key, value)| format!("{key}: {value}"))
        .collect();
    format!("{{{}}}", body.join(", "))
}

fn input(facts: &InputFacts) -> String {
    let unclaimed = facts
        .extracted
        .unclaimed
        .iter()
        .map(|(name, times)| (quoted(&hex(name)), times.to_string()))
        .collect();
    format!(
        "{{\"kind\": {}, \"sha256\": {}, \"lines\": {}, \"outside\": {}, \"header\": {}, \"unclaimed\": {}}}",
        quoted(facts.kind.name()),
        quoted(&facts.sha256),
        facts.extracted.lines,
        facts.extracted.outside,
        quoted(&hex(&facts.extracted.header)),
        object(unclaimed),
    )
}

fn projection(table: &Table, counts: &Counts) -> String {
    let mut values = Vec::new();
    for (field, name) in table.names.iter().enumerate() {
        let (reference, engine) = (counts.values[0][field], counts.values[1][field]);
        if reference > 0 || engine > 0 {
            let stated = format!("[{reference}, {engine}, {}]", counts.paired[field]);
            values.push((quoted(&hex(name)), stated));
        }
    }
    let shared = table
        .sources
        .iter()
        .zip(&counts.shared)
        .filter(|(_, &keys)| keys > 0)
        .map(|(source, keys)| (quoted(source), keys.to_string()))
        .collect();
    format!(
        "{{\"values\": {}, \"shared_keys\": {}, \"remainder\": [{}, {}]}}",
        object(values),
        object(shared),
        counts.left[0],
        counts.left[1],
    )
}

/// The hand-over's text: the inputs, the reference first, and per engine each projection's
/// counts and the rows written for the remainder.
pub fn pairs_json(inputs: &[InputFacts], engines: &[EngineFacts]) -> String {
    let mut text = String::from("{\n  \"inputs\": [\n");
    let stated: Vec<String> = inputs
        .iter()
        .map(|facts| format!("    {}", input(facts)))
        .collect();
    text.push_str(&stated.join(",\n"));
    text.push_str("\n  ],\n  \"engines\": [\n");
    let stated: Vec<String> = engines
        .iter()
        .map(|facts| {
            let mut entries: Vec<String> = PROJECTIONS
                .iter()
                .zip(&facts.counts)
                .map(|(name, counts)| {
                    format!("{}: {}", quoted(name), projection(&facts.table, counts))
                })
                .collect();
            entries.push(format!("\"rows\": [{}, {}]", facts.rows[0], facts.rows[1]));
            format!("    {{{}}}", entries.join(", "))
        })
        .collect();
    text.push_str(&stated.join(",\n"));
    text.push_str("\n  ]\n}\n");
    text
}
