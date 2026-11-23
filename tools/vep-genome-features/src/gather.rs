// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Gathers the features files under a directory into one table, a row a genome in
//! bytewise order of the genome's key.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::features::{Features, FEATURES_VERSION};

/// The name a features file carries under the gathered directory.
pub const FEATURES_FILE: &str = "features.json";

/// The table's columns, in order.
pub const COLUMNS: [&str; 20] = [
    "genome",
    "species",
    "assembly",
    "transcripts",
    "coding_transcripts",
    "genetic_codes",
    "biotypes",
    "incomplete_coding",
    "incomplete_coding_share",
    "sequences_with_features",
    "sequences_with_transcripts",
    "sequences_listed",
    "longest_sequence",
    "longest_sequence_length",
    "circular_sequences",
    "synonym_pairs",
    "variant_data",
    "sift_data",
    "regulatory_data",
    "exons_against_position",
];

/// The shape number of a features file, read ahead of its other keys: a file of another
/// shape can lack a key this build reads.
#[derive(Deserialize)]
struct Shape {
    features_version: u32,
}

/// Every features file under `dir`, at any depth, each directory's entries taken in
/// bytewise order of their names.
fn features_files(dir: &Path, found: &mut Vec<PathBuf>) -> Result<(), String> {
    let failed = |error: std::io::Error| format!("{}: {error}", dir.display());
    let mut entries: Vec<PathBuf> = Vec::new();
    for entry in fs::read_dir(dir).map_err(failed)? {
        entries.push(entry.map_err(failed)?.path());
    }
    entries.sort();
    for entry in entries {
        if entry.is_dir() {
            features_files(&entry, found)?;
        } else if entry.file_name().is_some_and(|name| name == FEATURES_FILE) {
            found.push(entry);
        }
    }
    Ok(())
}

/// A value as a table cell: `-` for none, and refused where it would break the row.
fn cell(file: &Path, column: &str, value: Option<&str>) -> Result<String, String> {
    match value {
        None | Some("") => Ok("-".to_owned()),
        Some(text) if text.contains(['\t', '\n', '\r']) => Err(format!(
            "{}: the {column} holds a tab or a line break",
            file.display()
        )),
        Some(text) => Ok(text.to_owned()),
    }
}

fn flag(set: bool) -> String {
    if set { "1" } else { "0" }.to_owned()
}

/// The row of one features file, without its genome key.
fn row(file: &Path, features: &Features) -> Result<Vec<String>, String> {
    let codes: Vec<String> = features
        .genetic_codes
        .iter()
        .map(|code| format!("{}:{}", code.code, code.coding_transcripts))
        .collect();
    let side = features.sequences.side_file.as_ref();
    let longest = side.and_then(|side| side.longest.as_ref());
    let count =
        |count: Option<u64>| count.map_or_else(|| "-".to_owned(), |count| count.to_string());
    Ok(vec![
        cell(file, "species", features.species.as_deref())?,
        cell(file, "assembly", features.assembly.as_deref())?,
        features.transcripts.to_string(),
        features.coding_transcripts.to_string(),
        if codes.is_empty() {
            "-".to_owned()
        } else {
            codes.join(",")
        },
        features.biotypes.len().to_string(),
        features.incomplete_coding.transcripts.to_string(),
        format!("{:.6}", features.incomplete_coding.share),
        features.sequences.with_features.to_string(),
        features.sequences.with_transcripts.to_string(),
        count(side.map(|side| side.listed)),
        cell(
            file,
            "longest sequence",
            longest.map(|longest| longest.name.as_str()),
        )?,
        count(longest.map(|longest| longest.length)),
        count(side.map(|side| side.circular.len() as u64)),
        features.synonyms.pairs.to_string(),
        flag(features.variant_data),
        flag(features.sift_data),
        flag(features.regulatory_data),
        features.exons_against_position.to_string(),
    ])
}

/// The table of every features file under `dir`: the header line, then a row a genome.
/// A genome's key is the one its file names, or the name of the directory that holds
/// the file where it names none.
pub fn gather(dir: &Path) -> Result<String, String> {
    let mut files = Vec::new();
    features_files(dir, &mut files)?;
    if files.is_empty() {
        return Err(format!("{}: no {FEATURES_FILE} under it", dir.display()));
    }
    let mut rows: BTreeMap<String, (PathBuf, Vec<String>)> = BTreeMap::new();
    for file in files {
        let failed = |what: String| format!("{}: {what}", file.display());
        let text = fs::read_to_string(&file).map_err(|error| failed(error.to_string()))?;
        let Shape { features_version } =
            serde_json::from_str(&text).map_err(|error| failed(error.to_string()))?;
        if features_version != FEATURES_VERSION {
            return Err(failed(format!(
                "features_version is {features_version}, and this build reads {FEATURES_VERSION}"
            )));
        }
        let features: Features =
            serde_json::from_str(&text).map_err(|error| failed(error.to_string()))?;
        let key = match features.genome.as_deref() {
            Some(key) if !key.is_empty() => key.to_owned(),
            _ => file
                .parent()
                .and_then(Path::file_name)
                .map(|name| name.to_string_lossy().into_owned())
                .unwrap_or_default(),
        };
        if key.is_empty() || key.contains(['\t', '\n', '\r']) {
            return Err(failed(
                "the genome key is empty or holds a tab or a line break".to_owned(),
            ));
        }
        let cells = row(&file, &features)?;
        if let Some((earlier, _)) = rows.get(&key) {
            return Err(failed(format!(
                "the genome {key} also stands in {}",
                earlier.display()
            )));
        }
        rows.insert(key, (file, cells));
    }
    let mut table = COLUMNS.join("\t");
    table.push('\n');
    for (key, (_, cells)) in rows {
        table.push_str(&key);
        for cell in cells {
            table.push('\t');
            table.push_str(&cell);
        }
        table.push('\n');
    }
    Ok(table)
}
