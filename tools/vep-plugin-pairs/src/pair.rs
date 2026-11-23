// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Pairs the reference's shards with an engine's, a shard a thread: under each key and
//! field the values equal as text pair one to one, and what neither side pairs is the
//! remainder, written as tuple lines in the order `sort` gives them under C collation.
//!
//! The feature projection keeps a value as often as its rows write it; the allele
//! projection drops the Feature and keeps a value once. A key under which a side writes
//! one field more than once is counted once for the field's source.

use crate::extract::Run;
use crate::input::{scratch_file, Kind};
use crate::order::take_varint;
use crate::resolve::Table;
use crate::rules::Res;
use crate::score::interrupt;
use crate::slots::Slot;
use memchr::memchr;
use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc;

pub const PROJECTIONS: [&str; 2] = ["feature", "allele"];
pub const SIDES: [&str; 2] = ["reference", "engine"];
const FEATURE: usize = 0;
const ALLELE: usize = 1;
const REFERENCE: usize = 0;
const ENGINE: usize = 1;

/// The counts of one projection: per field each side's values and the pairs equal as text,
/// per source the keys a side writes a field under more than once, and the remainder's
/// lines per side.
pub struct Counts {
    pub values: [Vec<u64>; 2],
    pub paired: Vec<u64>,
    pub shared: Vec<u64>,
    pub left: [u64; 2],
}

impl Counts {
    fn new(table: &Table) -> Counts {
        let fields = || vec![0; table.names.len()];
        Counts {
            values: [fields(), fields()],
            paired: fields(),
            shared: vec![0; table.sources.len()],
            left: [0; 2],
        }
    }

    fn add(&mut self, other: &Counts) {
        let sum = |into: &mut [u64], from: &[u64]| {
            for (total, count) in into.iter_mut().zip(from) {
                *total += count;
            }
        };
        for side in [REFERENCE, ENGINE] {
            sum(&mut self.values[side], &other.values[side]);
        }
        sum(&mut self.paired, &other.paired);
        sum(&mut self.shared, &other.shared);
        sum(&mut self.left, &other.left);
    }
}

/// One engine paired with the reference.
pub struct Paired {
    pub counts: [Counts; 2],
    /// The Locations of the feature remainder: each side's rows there are what a class can
    /// ask for.
    pub locations: HashSet<Vec<u8>>,
    /// The lines written as the rows of a tuple-file engine.
    pub tuple_rows: u64,
}

#[derive(Clone, Copy)]
struct Entry<'a> {
    rank: u32,
    field: u32,
    value: &'a [u8],
}

struct Row {
    /// Where the key starts in the shard.
    at: usize,
    key_len: u32,
}

struct Shard {
    bytes: Vec<u8>,
    rows: Vec<Row>,
}

impl Shard {
    fn read(path: &Path) -> io::Result<Shard> {
        let bytes = match std::fs::read(path) {
            Ok(bytes) => bytes,
            Err(error) if error.kind() == io::ErrorKind::NotFound => Vec::new(),
            Err(error) => return Err(error),
        };
        let mut rows = Vec::new();
        let mut at = 0;
        while at < bytes.len() {
            let length = take_varint(&bytes, &mut at) as usize;
            let end = at + length;
            let key_len = take_varint(&bytes, &mut at) as u32;
            rows.push(Row { at, key_len });
            at = end;
        }
        let key = |row: &Row| &bytes[row.at..row.at + row.key_len as usize];
        rows.sort_unstable_by(|a, b| key(a).cmp(key(b)));
        Ok(Shard { bytes, rows })
    }

    fn key(&self, row: usize) -> &[u8] {
        let row = &self.rows[row];
        &self.bytes[row.at..row.at + row.key_len as usize]
    }

    /// The end of the run of rows from `from` whose keys start with `prefix`.
    fn run(&self, from: usize, prefix: &[u8]) -> usize {
        let mut end = from;
        while end < self.rows.len() && self.key(end).starts_with(prefix) {
            end += 1;
        }
        end
    }

    fn entries<'a>(
        &'a self,
        rows: std::ops::Range<usize>,
        table: &Table,
        into: &mut Vec<Entry<'a>>,
    ) {
        into.clear();
        for row in &self.rows[rows] {
            let mut at = row.at + row.key_len as usize;
            for _ in 0..take_varint(&self.bytes, &mut at) {
                let field = take_varint(&self.bytes, &mut at) as u32;
                let length = take_varint(&self.bytes, &mut at) as usize;
                into.push(Entry {
                    rank: table.rank[field as usize],
                    field,
                    value: &self.bytes[at..at + length],
                });
                at += length;
            }
        }
        into.sort_unstable_by(|a, b| a.rank.cmp(&b.rank).then_with(|| a.value.cmp(b.value)));
    }
}

/// The key without its Feature: the Location and the Allele, a tab after each.
fn allele_key(key: &[u8]) -> &[u8] {
    let first = memchr(b'\t', key).expect("a key has three columns");
    let second = memchr(b'\t', &key[first + 1..]).expect("a key has three columns");
    &key[..first + second + 2]
}

/// The smaller of the two sides' next keys, where a side has one.
fn next_key<'a>(one: Option<&'a [u8]>, other: Option<&'a [u8]>) -> &'a [u8] {
    match (one, other) {
        (Some(one), Some(other)) => one.min(other),
        (Some(key), None) | (None, Some(key)) => key,
        (None, None) => unreachable!("a side has a row left"),
    }
}

#[derive(Default)]
struct ShardOut {
    /// Per projection and side, the remainder's lines.
    remainder: [[Vec<u8>; 2]; 2],
    tuple_rows: Vec<u8>,
    tuple_row_count: u64,
    locations: Vec<Vec<u8>>,
}

struct Pairer<'a> {
    table: &'a Table,
    counts: [Counts; 2],
    out: ShardOut,
    /// The sources already counted as shared under the key in hand.
    shared: Vec<u32>,
}

impl Pairer<'_> {
    /// Pairs two sides' sorted entries under `prefix`, a key or a key without its Feature,
    /// and says whether a value is left over.
    fn pair(&mut self, projection: usize, prefix: &[u8], sides: [&[Entry]; 2]) -> bool {
        let table = self.table;
        let counts = &mut self.counts[projection];
        let [reference, engine] = sides;
        let (mut r, mut e) = (0, 0);
        let mut left = false;
        self.shared.clear();
        while r < reference.len() || e < engine.len() {
            let rank = match (reference.get(r), engine.get(e)) {
                (Some(one), Some(other)) => one.rank.min(other.rank),
                (Some(entry), None) | (None, Some(entry)) => entry.rank,
                (None, None) => unreachable!(),
            };
            let run = |entries: &[Entry], from: usize| {
                from + entries[from..]
                    .iter()
                    .take_while(|entry| entry.rank == rank)
                    .count()
            };
            let (r_end, e_end) = (run(reference, r), run(engine, e));
            let field = if r_end > r {
                reference[r].field
            } else {
                engine[e].field
            } as usize;
            counts.values[REFERENCE][field] += (r_end - r) as u64;
            counts.values[ENGINE][field] += (e_end - e) as u64;
            if r_end - r > 1 || e_end - e > 1 {
                let source = table.source[field];
                if !self.shared.contains(&source) {
                    self.shared.push(source);
                    counts.shared[source as usize] += 1;
                }
            }
            let mut leave = |side: usize, value: &[u8]| {
                let lines = &mut self.out.remainder[projection][side];
                lines.extend_from_slice(prefix);
                lines.extend_from_slice(&table.names[field]);
                lines.push(b'\t');
                lines.extend_from_slice(value);
                lines.push(b'\n');
                counts.left[side] += 1;
                left = true;
            };
            while r < r_end && e < e_end {
                match reference[r].value.cmp(engine[e].value) {
                    std::cmp::Ordering::Equal => {
                        counts.paired[field] += 1;
                        r += 1;
                        e += 1;
                    }
                    std::cmp::Ordering::Less => {
                        leave(REFERENCE, reference[r].value);
                        r += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        leave(ENGINE, engine[e].value);
                        e += 1;
                    }
                }
            }
            for entry in &reference[r..r_end] {
                leave(REFERENCE, entry.value);
            }
            for entry in &engine[e..e_end] {
                leave(ENGINE, entry.value);
            }
            (r, e) = (r_end, e_end);
        }
        left
    }

    fn shard(&mut self, reference: &Shard, engine: &Shard, tuple_rows: bool) {
        let table = self.table;
        let (mut r, mut e) = (0, 0);
        let mut entries: [Vec<Entry>; 2] = [Vec::new(), Vec::new()];
        while r < reference.rows.len() || e < engine.rows.len() {
            let key = next_key(
                (r < reference.rows.len()).then(|| reference.key(r)),
                (e < engine.rows.len()).then(|| engine.key(e)),
            );
            let alleles = allele_key(key);
            let (r_end, e_end) = (reference.run(r, alleles), engine.run(e, alleles));
            reference.entries(r..r_end, table, &mut entries[REFERENCE]);
            engine.entries(e..e_end, table, &mut entries[ENGINE]);
            for side in &mut entries {
                side.dedup_by(|a, b| a.rank == b.rank && a.value == b.value);
            }
            self.pair(ALLELE, alleles, [&entries[REFERENCE], &entries[ENGINE]]);
            while r < r_end || e < e_end {
                let key = next_key(
                    (r < r_end).then(|| reference.key(r)),
                    (e < e_end).then(|| engine.key(e)),
                );
                let (r_next, e_next) = (reference.run(r, key), engine.run(e, key));
                reference.entries(r..r_next, table, &mut entries[REFERENCE]);
                engine.entries(e..e_next, table, &mut entries[ENGINE]);
                if self.pair(FEATURE, key, [&entries[REFERENCE], &entries[ENGINE]]) {
                    let location = &key[..memchr(b'\t', key).expect("a key has three columns")];
                    if self
                        .out
                        .locations
                        .last()
                        .is_none_or(|last| last != location)
                    {
                        self.out.locations.push(location.to_vec());
                    }
                    if tuple_rows {
                        for entry in &entries[ENGINE] {
                            let lines = &mut self.out.tuple_rows;
                            lines.extend_from_slice(key);
                            lines.extend_from_slice(&table.names[entry.field as usize]);
                            lines.push(b'\t');
                            lines.extend_from_slice(entry.value);
                            lines.push(b'\n');
                        }
                        self.out.tuple_row_count += entries[ENGINE].len() as u64;
                    }
                }
                (r, e) = (r_next, e_next);
            }
        }
    }
}

/// Pairs the engine that is the run's input `number` with the reference, input 0, and
/// writes the remainder's files, `input<number>.<projection>.<side>`, and with `rows` the
/// rows of a tuple-file engine, `input<number>.rows.engine`, into `out`. The caller's slot,
/// `own`, is lent while the pass's threads pair, and the caller takes one for each turn in
/// which it writes the shards that ended.
pub fn pair(
    run: &Run,
    table: &Table,
    number: usize,
    kind: Kind,
    rows: bool,
    out: &Path,
    own: &mut Slot,
) -> Res<Paired> {
    let tuple_rows = rows && kind == Kind::Tuples;
    let mut files = Vec::new();
    for projection in PROJECTIONS {
        for side in SIDES {
            let path = scratch_file(out, number, &format!("{projection}.{side}"));
            files.push(BufWriter::with_capacity(1 << 20, File::create(path)?));
        }
    }
    let mut tuple_file = match tuple_rows {
        true => Some(BufWriter::with_capacity(
            1 << 20,
            File::create(scratch_file(out, number, "rows.engine"))?,
        )),
        false => None,
    };
    let mut paired = Paired {
        counts: [Counts::new(table), Counts::new(table)],
        locations: HashSet::new(),
        tuple_rows: 0,
    };
    let next = AtomicUsize::new(0);
    let (send, receive) = mpsc::channel::<(usize, io::Result<ShardOut>)>();
    let shards = run.shards();
    let lent = own.lent();
    let paired_all = std::thread::scope(|scope| -> Res<()> {
        // Dropped on every way out, after which a thread's send fails and the thread ends.
        let receive = receive;
        let workers: Vec<_> = (0..run.threads.min(shards).max(1))
            .map(|_| {
                let send = send.clone();
                let next = &next;
                scope.spawn(move || {
                    let mut pairer = Pairer {
                        table,
                        counts: [Counts::new(table), Counts::new(table)],
                        out: ShardOut::default(),
                        shared: Vec::new(),
                    };
                    loop {
                        let shard = next.fetch_add(1, Ordering::Relaxed);
                        if shard >= shards || interrupt::requested() {
                            break;
                        }
                        let slot = run.slots.take();
                        let read = |input: usize| {
                            Shard::read(&scratch_file(run.scratch, input, &format!("shard{shard}")))
                        };
                        let done = read(0).and_then(|reference| {
                            let engine = read(number)?;
                            pairer.shard(&reference, &engine, tuple_rows);
                            Ok(std::mem::take(&mut pairer.out))
                        });
                        drop(slot);
                        if send.send((shard, done)).is_err() {
                            break;
                        }
                    }
                    pairer.counts
                })
            })
            .collect();
        drop(send);
        // Shards end in any order and are written in their own, which is the lines' order.
        let mut waiting: BTreeMap<usize, ShardOut> = BTreeMap::new();
        let mut due = 0;
        while let Ok(ended) = receive.recv() {
            // One turn writes every shard that ended before it.
            let _slot = run.slots.take();
            for (shard, done) in std::iter::once(ended).chain(receive.try_iter()) {
                waiting.insert(shard, done?);
                while let Some(found) = waiting.remove(&due) {
                    for (file, lines) in files.iter_mut().zip(found.remainder.iter().flatten()) {
                        file.write_all(lines)?;
                    }
                    if let Some(file) = tuple_file.as_mut() {
                        file.write_all(&found.tuple_rows)?;
                    }
                    paired.tuple_rows += found.tuple_row_count;
                    paired.locations.extend(found.locations);
                    due += 1;
                }
            }
        }
        for worker in workers {
            let counts = worker.join().expect("a pairing thread panicked");
            for (total, counted) in paired.counts.iter_mut().zip(&counts) {
                total.add(counted);
            }
        }
        Ok(())
    });
    drop(lent);
    paired_all?;
    for file in files.iter_mut().chain(tuple_file.as_mut()) {
        file.flush()?;
    }
    Ok(paired)
}
