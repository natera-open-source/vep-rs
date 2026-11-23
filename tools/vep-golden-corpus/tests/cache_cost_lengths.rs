// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Holds the cost `select --cache` gives a transcript record to
//! `tests/fixtures/cache_cost/lengths.tsv`: for every cache the table names, the count of
//! its shards and of its records, the sum of the records' costs, and a hash of the text
//! that lists every record's cost, as CPython's
//! `len(zlib.compress(json.dumps(record, separators=(",", ":"), ensure_ascii=False).encode()))`
//! gives them. A cache of a committed corpus that the table does not name is not read.

mod common;

use common::{fixture, repository};
use std::collections::BTreeSet;
use std::fs;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use vep_golden_corpus::cache_cost::shard_transcripts;

/// The FNV-1a hash of 64 bits of `data`, carried on from `state`.
fn fnv1a64(state: u64, data: &[u8]) -> u64 {
    data.iter().fold(state, |state, byte| {
        (state ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// The shard files of a cache as `(chromosome directory, file name)`, in that order.
fn shards(cache: &Path) -> Vec<(String, String)> {
    let mut found = Vec::new();
    for chromosome in fs::read_dir(cache.join("transcripts")).unwrap() {
        let chromosome = chromosome.unwrap().path();
        for shard in fs::read_dir(&chromosome).unwrap() {
            let name = shard.unwrap().file_name().into_string().unwrap();
            if name.ends_with(".json") || name.ends_with(".json.gz") {
                let directory = chromosome.file_name().unwrap().to_str().unwrap();
                found.push((directory.to_string(), name));
            }
        }
    }
    found.sort();
    found
}

/// One row of the table for the cache at `cache`: shards, records, summed costs, hash.
fn measured(cache: &Path) -> String {
    let shards = shards(cache);
    let (mut records, mut bytes) = (0, 0);
    let mut state = 0xcbf2_9ce4_8422_2325;
    for (chromosome, name) in &shards {
        let path = cache.join("transcripts").join(chromosome).join(name);
        let transcripts = shard_transcripts(path.to_str().unwrap()).unwrap();
        let costs: Vec<String> = transcripts
            .iter()
            .map(|(_, _, cost)| cost.to_string())
            .collect();
        records += costs.len();
        bytes += transcripts.iter().map(|(_, _, cost)| cost).sum::<usize>();
        let line = format!("transcripts/{chromosome}/{name}\t{}\n", costs.join(","));
        state = fnv1a64(state, line.as_bytes());
    }
    format!("{}\t{records}\t{bytes}\t{state:016x}", shards.len())
}

#[test]
fn every_record_of_the_listed_caches_costs_what_zlib_under_cpython_gives() {
    let table = fs::read_to_string(fixture("cache_cost/lengths.tsv")).unwrap();
    let mut lines = table.lines();
    assert_eq!(lines.next(), Some("cache\tshards\trecords\tbytes\tfnv1a64"));
    let rows: Vec<(&str, &str)> = lines.map(|line| line.split_once('\t').unwrap()).collect();
    assert_eq!(rows.len(), 113);
    let root = repository().join("tests/golden/116");
    let mut named = BTreeSet::new();
    for (cache, _) in &rows {
        assert!(named.insert(*cache), "{cache} has two rows");
        assert!(root.join(cache).join("transcripts").is_dir(), "{cache}");
    }
    let next = AtomicUsize::new(0);
    let threads = std::thread::available_parallelism().map_or(1, |count| count.get().min(8));
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                while let Some((cache, expected)) = rows.get(next.fetch_add(1, Ordering::Relaxed)) {
                    assert_eq!(measured(&root.join(cache)), *expected, "{cache}");
                }
            });
        }
    });
    let records: usize = rows
        .iter()
        .map(|(_, row)| row.split('\t').nth(1).unwrap().parse::<usize>().unwrap())
        .sum();
    assert_eq!(records, 88_528);
}
