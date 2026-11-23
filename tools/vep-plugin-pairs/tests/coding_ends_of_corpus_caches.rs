// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The plugin comparison's reading of a cache record's coding ends
//! (`rules::downstream::ends::coding_ends_of`), held to the committed table
//! `tests/fixtures/coding_ends/verdicts.tsv` on every coding transcript of the corpus caches
//! the table names. The table is written by vep-rs's cache loader and consequence engine
//! (`crates/vep-cli/tests/coding_ends_of_corpus_caches.rs`), which a test holds to it as
//! well, so the stop codon and the start codon this tool reads at a transcript's coding
//! ends are the ones vep-rs writes its terms by.
//!
//! A cache is read here as the loader reads it: every sequence directory under
//! `transcripts/` whose name opens with no dot, and a record a second shard repeats read
//! once, by its `dbID` within its `source_cache`, or without a `dbID` by its stable id and
//! version. A coding transcript is one whose record holds a coding sequence that is not
//! empty.

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use vep_plugin_pairs::pytext::Json;
use vep_plugin_pairs::rules::cache::TranscriptShards;
use vep_plugin_pairs::rules::downstream::ends::coding_ends_of;
use vep_plugin_pairs::rules::shared::cache_integer;

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

/// What the loader reads a record once by: its `source_cache` and `dbID`, or without a
/// `dbID` its `source_cache`, stable id and version.
type RecordKey = (Option<String>, Option<i128>, String, Option<i128>);

fn is_coding(record: &Json) -> bool {
    let coding = record
        .get("variation_effect_feature_cache")
        .and_then(|cached| cached.get("translateable_seq"));
    matches!(coding, Some(Json::Str(coding)) if !coding.is_empty())
}

#[test]
fn the_coding_ends_the_plugin_comparison_reads_are_the_committed_table() {
    let caches = committed();
    assert!(!caches.is_empty(), "the table names no cache");
    let mut differences: Vec<String> = Vec::new();
    for (cache, count, lines) in &caches {
        let directory = repository().join("tests/golden/116").join(cache);
        let shards = TranscriptShards::new(&directory).expect("the cache");
        let (mut coding, mut read) = (0, Vec::new());
        let sequences = shards.chromosomes().iter();
        for sequence in sequences.filter(|name| !name.starts_with('.')) {
            let mut seen: HashSet<RecordKey> = HashSet::new();
            let visit = |record: Json| {
                let Some((stable_id, ends)) = coding_ends_of(&record) else {
                    return Ok(());
                };
                let source = match record.get("source_cache") {
                    Some(Json::Str(source)) => Some(source.clone()),
                    _ => None,
                };
                let key = match cache_integer(record.get("dbID")) {
                    Some(id) => (source, Some(id), String::new(), None),
                    None => (
                        source,
                        None,
                        stable_id.to_string(),
                        cache_integer(record.get("version")),
                    ),
                };
                if !seen.insert(key) || !is_coding(&record) {
                    return Ok(());
                }
                coding += 1;
                let (stop, start) = (ends.stop.to_string(), ends.start.to_string());
                if stop != "whole" || start != "shown" {
                    read.push([sequence.as_str(), stable_id, &stop, &start].join("\t"));
                }
                Ok(())
            };
            shards
                .of_chromosome(sequence, visit)
                .expect("the shards read");
        }
        if coding != *count {
            differences.push(format!(
                "{cache}: {coding} coding transcripts read, the table counts {count}"
            ));
        }
        let (committed, read_lines): (HashSet<&String>, HashSet<&String>) =
            (lines.iter().collect(), read.iter().collect());
        let mut apart: Vec<String> = committed
            .difference(&read_lines)
            .map(|line| format!("{cache}: the table alone holds {line:?}"))
            .chain(
                read_lines
                    .difference(&committed)
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
