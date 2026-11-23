// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison's reading of a cache record's coding ends (`coding_ends`), held to the
//! committed table `tests/fixtures/coding_ends/verdicts.tsv` on every transcript of the
//! corpus caches the table names. The table is written by vep-rs's cache loader and
//! consequence engine (`crates/vep-cli/tests/coding_ends_of_corpus_caches.rs`), which a
//! test holds to it as well, so the stop codon and the start codon the comparison reads at
//! a transcript's coding ends are the ones vep-rs writes its terms by.
//!
//! A cache is read here as the loader reads it: every sequence directory under
//! `transcripts/` whose name opens with no dot, and a record a second shard repeats read
//! once, by its `dbID` within its `source_cache`, or without a `dbID` by its stable id and
//! version. A coding transcript is one whose record holds a coding sequence that is not
//! empty.

use std::collections::{BTreeMap, HashSet};
use std::path::{Path, PathBuf};
use vep_formats_compare::cache;
use vep_formats_compare::coding_ends::CodingEnds;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::partial_codon_models::{coding_ends_of, PartialCodonModels};
use vep_formats_compare::rules::transcript_model::cache_integer;

fn repository() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The table's caches in its order, each with the count of its coding transcripts and its
/// transcript lines less their first two cells.
fn committed() -> Vec<(String, usize, Vec<String>)> {
    let path = repository().join("tests/fixtures/coding_ends/verdicts.tsv");
    let table = std::fs::read_to_string(&path).expect("the committed table");
    let mut caches: Vec<(String, usize, Vec<String>)> = Vec::new();
    for line in table.lines().filter(|line| !line.starts_with('#')) {
        let cells: Vec<&str> = line.split('\t').collect();
        match cells.as_slice() {
            ["cache", cache, count] => {
                caches.push((
                    cache.to_string(),
                    count.parse().expect("a count"),
                    Vec::new(),
                ));
            }
            ["transcript", cache, sequence, stable_id, stop, start] => {
                let (named, _, lines) = caches.last_mut().expect("a cache line first");
                assert_eq!(named, cache, "a transcript line under another cache's line");
                lines.push([*sequence, *stable_id, *stop, *start].join("\t"));
            }
            _ => panic!("the table holds a line of neither kind: {line}"),
        }
    }
    caches
}

/// The sequence directories of a cache's `transcripts/`, in the order of their names.
fn sequences(root: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(root)
        .expect("the cache's transcripts/")
        .map(|entry| entry.expect("a directory entry"))
        .filter(|entry| entry.path().is_dir())
        .map(|entry| entry.file_name().into_string().expect("a UTF-8 name"))
        .filter(|name| !name.starts_with('.'))
        .collect();
    names.sort();
    names
}

/// What the loader reads a record once by: its `source_cache` and `dbID`, or without a
/// `dbID` its `source_cache`, stable id and version.
type RecordKey = (Option<String>, Option<i64>, String, Option<i64>);

/// The records of one sequence, each once, with the stable id and the coding ends the
/// comparison reads of it; a record without a stable id that is text is left out.
fn records(root: &Path, sequence: &str) -> Vec<(Value, String, CodingEnds)> {
    let shards = cache::shards(root, Some(sequence)).expect("the shards");
    let read = cache::shards_transcripts(&shards, cache::threads()).expect("the shards read");
    let mut seen: HashSet<RecordKey> = HashSet::new();
    let mut kept = Vec::new();
    for record in read.into_iter().flatten() {
        let Some((stable_id, ends)) = coding_ends_of(&record) else {
            continue;
        };
        let source = record
            .get("source_cache")
            .and_then(Value::as_str)
            .map(str::to_string);
        let key = match record.get("dbID").and_then(cache_integer) {
            Some(id) => (source, Some(id), String::new(), None),
            None => (
                source,
                None,
                stable_id.clone(),
                record.get("version").and_then(cache_integer),
            ),
        };
        if seen.insert(key) {
            kept.push((record, stable_id, ends));
        }
    }
    kept
}

fn is_coding(record: &Value) -> bool {
    record
        .get("variation_effect_feature_cache")
        .and_then(|cached| cached.get("translateable_seq"))
        .and_then(Value::as_str)
        .is_some_and(|coding| !coding.is_empty())
}

#[test]
fn the_coding_ends_the_comparison_reads_are_the_committed_table() {
    let caches = committed();
    assert!(!caches.is_empty(), "the table names no cache");
    let mut differences: Vec<String> = Vec::new();
    for (cache, count, lines) in &caches {
        let root = repository()
            .join("tests/golden/116")
            .join(cache)
            .join("transcripts");
        let (mut coding, mut read) = (0, Vec::new());
        for sequence in sequences(&root) {
            for (record, stable_id, ends) in records(&root, &sequence) {
                if !is_coding(&record) {
                    continue;
                }
                coding += 1;
                let (stop, start) = (ends.stop.to_string(), ends.start.to_string());
                if stop != "whole" || start != "shown" {
                    read.push([sequence.as_str(), &stable_id, &stop, &start].join("\t"));
                }
            }
        }
        read.sort();
        if coding != *count {
            differences.push(format!(
                "{cache}: {coding} coding transcripts read, the table counts {count}"
            ));
        }
        let (committed, read): (HashSet<&String>, HashSet<&String>) =
            (lines.iter().collect(), read.iter().collect());
        let mut apart: Vec<String> = committed
            .difference(&read)
            .map(|line| format!("{cache}: the table alone holds {line:?}"))
            .chain(
                read.difference(&committed)
                    .map(|line| format!("{cache}: the comparison alone reads {line:?}")),
            )
            .collect();
        apart.sort();
        if lines.len() != read.len() && apart.is_empty() {
            apart.push(format!(
                "{cache}: {} lines read, the table holds {}",
                read.len(),
                lines.len()
            ));
        }
        differences.extend(apart);
    }
    assert!(
        differences.is_empty(),
        "{} differences from the committed table, the first of them:\n{}",
        differences.len(),
        differences
            .iter()
            .take(20)
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// The reader a comparison asks answers a row that lies on a record of a transcript with
/// the coding ends of that record, and with none where another record of the stable id
/// on the row's sequence holds the row's position too.
#[test]
fn a_row_of_each_corpus_transcript_is_answered_with_its_record() {
    let (mut answered, mut unanswered) = (0, 0);
    for (cache, _, _) in &committed() {
        let directory = repository().join("tests/golden/116").join(cache);
        let root = directory.join("transcripts");
        let models = PartialCodonModels::new(&directory).expect("the cache");
        for sequence in sequences(&root) {
            let mut by_id: BTreeMap<String, Vec<CodingEnds>> = BTreeMap::new();
            let shards = cache::shards(&root, Some(&sequence)).expect("the shards");
            let read = cache::shards_transcripts(&shards, cache::threads()).expect("the shards");
            for (stable_id, ends) in read.iter().flatten().filter_map(coding_ends_of) {
                let readings = by_id.entry(stable_id).or_default();
                if !readings.contains(&ends) {
                    readings.push(ends);
                }
            }
            for (stable_id, readings) in &by_id {
                let row_at = |position: i64| Entry {
                    key: "k".to_string(),
                    location: Some(format!("{sequence}:{position}")),
                    feature: stable_id.clone(),
                    ..Entry::default()
                };
                for ends in readings {
                    let (first, last) = ends.span.expect("a record states its span");
                    for position in [first, last] {
                        let holding: Vec<&CodingEnds> = readings
                            .iter()
                            .filter(|other| {
                                other.span.is_some_and(|(start, end)| {
                                    start <= position && position <= end
                                })
                            })
                            .collect();
                        let expected = match holding.as_slice() {
                            [only] => Some((*only).clone()),
                            _ => None,
                        };
                        match expected {
                            Some(_) => answered += 1,
                            None => unanswered += 1,
                        }
                        assert_eq!(
                            models
                                .coding_ends(&row_at(position))
                                .expect("the cache reads"),
                            expected,
                            "{cache} {sequence} {stable_id} at {position}"
                        );
                    }
                    assert!(
                        models
                            .coding_ends(&row_at(first - 1))
                            .expect("the cache reads")
                            .is_none_or(|other| other != *ends),
                        "{cache} {sequence} {stable_id} before its first base"
                    );
                }
                assert_eq!(
                    models.holds(&row_at(1)).expect("the cache reads"),
                    readings.iter().any(|ends| ends.tail_length != 0),
                    "{cache} {sequence} {stable_id}"
                );
            }
        }
    }
    assert!(answered > 0);
    println!("{answered} rows answered with their record, {unanswered} with none");
}
