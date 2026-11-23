// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The JSON files kept under `tests/fixtures/`, a made-up corpus's `manifest.json` and the
//! `expected.<n>.json` a call writes among them: a record names its dataset under
//! `dataset`, the key `select` writes and `check` reads, and no object of a file holds
//! `suite`, `suites` or `multi_allelic_per_suite`, keys no subcommand reads.

mod common;

use common::ROOT;
use std::fs;
use std::path::{Path, PathBuf};
use vep_golden_corpus::json::{loads, Value};

/// Keys no subcommand reads: a manifest states a record's dataset under `dataset`, its
/// datasets under `datasets` and the multi-allelic sample of each under
/// `multi_allelic_per_dataset`. `classify` carries a record's keys as it finds them, so a
/// fixture that states one of these is written back with it and no call fails on it.
const UNREAD_KEYS: [&str; 3] = ["suite", "suites", "multi_allelic_per_suite"];

/// Every `.json` file under `directory`.
fn json_files(directory: &Path, found: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(directory).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            json_files(&path, found);
        } else if path.extension().is_some_and(|kind| kind == "json") {
            found.push(path);
        }
    }
}

/// Walks `value`, whose place in its file is `place`: each key of `UNREAD_KEYS` goes to
/// `unread` with its place, and the count of objects that hold `dataset` is returned.
fn walk(value: &Value, place: &str, unread: &mut Vec<String>) -> usize {
    match value {
        Value::Object(object) => {
            let mut datasets = usize::from(object.contains_key("dataset"));
            for (key, held) in object.iter() {
                let inside = format!("{place}/{key}");
                if UNREAD_KEYS.contains(&key) {
                    unread.push(inside.clone());
                }
                datasets += walk(held, &inside, unread);
            }
            datasets
        }
        Value::List(items) => items
            .iter()
            .enumerate()
            .map(|(index, item)| walk(item, &format!("{place}/{index}"), unread))
            .sum(),
        _ => 0,
    }
}

#[test]
fn no_fixture_names_a_dataset_under_a_key_no_subcommand_reads() {
    let fixtures = Path::new(ROOT).join("tests/fixtures");
    let mut files = Vec::new();
    json_files(&fixtures, &mut files);
    files.sort();
    let mut unread = Vec::new();
    let mut datasets = 0;
    for file in &files {
        let name = file.strip_prefix(&fixtures).unwrap().display().to_string();
        let text = fs::read_to_string(file).unwrap();
        let value = loads(&text).unwrap_or_else(|error| panic!("{name}: {error:?}"));
        datasets += walk(&value, &name, &mut unread);
    }
    // The walk reads the manifests: the made-up corpora name a dataset on their records.
    assert!(datasets > 0, "no fixture record names a dataset");
    assert_eq!(unread, Vec::<String>::new());
}
