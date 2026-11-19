// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The extraction of `score_plugin_tuples.py`, block by block on every thread: the
//! populated fields of each row that a source of the run claims, written to the shard of
//! the row's Location as one record per row.
//!
//! A record is its length, then the key (Location, Allele and Feature, a tab after each),
//! then the row's fields as their numbers and values, every length a varint.

use crate::input::{lines, scratch_file, Input, Kind};
use crate::order::{canonical_location, cmp_column, decode_extra, put_varint, strip};
use crate::resolve::{Registry, Source};
use crate::{Refusal, Res};
use memchr::{memchr, memchr_iter, memmem};
use std::collections::HashMap;
use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// A shard buffer is written out once it holds this much.
const FLUSH_BYTES: usize = 32 << 10;
pub const DEFAULT_HEADER: &[u8] = b"#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";

/// What every pass of a run shares.
pub struct Run<'a> {
    pub registry: &'a Registry,
    /// The Locations that start each shard but the first, ascending.
    pub splitters: &'a [Vec<u8>],
    pub scratch: &'a Path,
    pub threads: usize,
    pub distinct_keys: usize,
}

impl Run<'_> {
    pub fn shards(&self) -> usize {
        self.splitters.len() + 1
    }
}

/// The canonical form and the shard of the Location last set, kept while the rows that
/// follow write the same one.
pub struct Locator<'a> {
    run: &'a Run<'a>,
    written: Vec<u8>,
    pub location: Vec<u8>,
    pub shard: usize,
}

impl<'a> Locator<'a> {
    pub fn new(run: &'a Run<'a>) -> Locator<'a> {
        let mut locator = Locator {
            run,
            written: Vec::new(),
            location: Vec::new(),
            shard: 0,
        };
        locator.place();
        locator
    }

    pub fn set(&mut self, written: &[u8]) -> Res<()> {
        if written != self.written {
            self.location = match canonical_location(written) {
                Some(location) => location,
                None => self.run.registry.location(written)?,
            };
            self.written.clear();
            self.written.extend_from_slice(written);
            self.place();
        }
        Ok(())
    }

    fn place(&mut self) {
        let location = &self.location;
        self.shard = self
            .run
            .splitters
            .partition_point(|splitter| cmp_column(splitter, location).is_le());
    }
}

/// One input after extraction, as the hand-over states it.
pub struct Extracted {
    /// The data rows of a default-format input, the lines of a tuple file.
    pub lines: u64,
    pub outside: u64,
    pub unclaimed: Vec<(Vec<u8>, u64)>,
    pub header: Vec<u8>,
    /// Per block, the data rows before it.
    pub rows_before: Vec<u64>,
}

enum Fault {
    Width(usize),
    Header(usize),
    DistinctKeys,
    Told(String),
}

struct Worker<'a> {
    run: &'a Run<'a>,
    input: &'a Input,
    number: usize,
    locks: &'a [Mutex<()>],
    locator: Locator<'a>,
    sources: HashMap<Vec<u8>, Source>,
    shards: Vec<Vec<u8>>,
    key: Vec<u8>,
    fields: Vec<u8>,
    decoded: Vec<u8>,
    /// The block and the line in hand.
    at: (usize, u64),
    outside: u64,
    /// Per key no source claims: the times it is written, and the block and line met first.
    unclaimed: HashMap<Vec<u8>, (u64, (usize, u64))>,
    /// Per block read: its number, its lines and its data rows.
    counted: Vec<(usize, u64, u64)>,
    header: Option<(usize, u64, Vec<u8>)>,
    fault: Option<(usize, u64, Fault)>,
}

impl Worker<'_> {
    fn source_of(&mut self, name: &[u8]) -> Result<Source, Fault> {
        if let Some(&source) = self.sources.get(name) {
            return Ok(source);
        }
        let source = self
            .run
            .registry
            .source_of(name)
            .map_err(|refusal| Fault::Told(refusal.0))?;
        self.sources.insert(name.to_vec(), source);
        Ok(source)
    }

    fn unclaimed(&mut self, name: &[u8]) -> Result<(), Fault> {
        match self.unclaimed.get_mut(name) {
            Some((times, _)) => *times += 1,
            None => {
                self.unclaimed.insert(name.to_vec(), (1, self.at));
                if self.unclaimed.len() > self.run.distinct_keys {
                    return Err(Fault::DistinctKeys);
                }
            }
        }
        Ok(())
    }

    fn locate(&mut self, written: &[u8]) -> Result<(), Fault> {
        self.locator
            .set(written)
            .map_err(|refusal| Fault::Told(refusal.0))
    }

    fn start_key(&mut self, allele: &[u8], feature: &[u8]) {
        self.key.clear();
        for column in [
            self.locator.location.as_slice(),
            strip(allele),
            strip(feature),
        ] {
            self.key.extend_from_slice(column);
            self.key.push(b'\t');
        }
    }

    fn field(&mut self, number: u32, value: &[u8]) {
        put_varint(&mut self.fields, u64::from(number));
        put_varint(&mut self.fields, value.len() as u64);
        self.fields.extend_from_slice(value);
    }

    /// Writes the row's record, its `count` fields gathered in `self.fields`, to its shard.
    fn record(&mut self, count: u64) -> std::io::Result<()> {
        let mut head = Vec::with_capacity(20);
        put_varint(&mut head, self.key.len() as u64);
        let mut tail = Vec::with_capacity(10);
        put_varint(&mut tail, count);
        let shard = self.locator.shard;
        let buffer = &mut self.shards[shard];
        let length = head.len() + self.key.len() + tail.len() + self.fields.len();
        put_varint(buffer, length as u64);
        buffer.extend_from_slice(&head);
        buffer.extend_from_slice(&self.key);
        buffer.extend_from_slice(&tail);
        buffer.extend_from_slice(&self.fields);
        self.fields.clear();
        if buffer.len() >= FLUSH_BYTES {
            self.flush(shard)?;
        }
        Ok(())
    }

    fn flush(&mut self, shard: usize) -> std::io::Result<()> {
        if self.shards[shard].is_empty() {
            return Ok(());
        }
        let path = scratch_file(self.run.scratch, self.number, &format!("shard{shard}"));
        let _held = self.locks[shard].lock().unwrap();
        let mut file = OpenOptions::new().create(true).append(true).open(path)?;
        file.write_all(&self.shards[shard])?;
        self.shards[shard].clear();
        Ok(())
    }

    /// One default-format row: `extract` of the scorer from the split of its columns on.
    fn default_row(&mut self, row: &[u8]) -> Result<(), Fault> {
        let mut tabs = [0usize; 13];
        let mut found = 0;
        for tab in memchr_iter(b'\t', row) {
            if found < tabs.len() {
                tabs[found] = tab;
            }
            found += 1;
        }
        if found != tabs.len() {
            return Err(Fault::Width(found + 1));
        }
        self.locate(&row[tabs[0] + 1..tabs[1]])?;
        self.start_key(&row[tabs[1] + 1..tabs[2]], &row[tabs[3] + 1..tabs[4]]);
        let extra = &row[tabs[12] + 1..];
        if extra == b"-" {
            return Ok(());
        }
        let mut count = 0;
        for pair in extra.split(|&byte| byte == b';') {
            let Some(equals) = memchr(b'=', pair) else {
                continue;
            };
            let (name, value) = (&pair[..equals], &pair[equals + 1..]);
            if value.is_empty() || value == b"-" {
                continue;
            }
            match self.source_of(name)? {
                Source::Core => {}
                Source::Unclaimed => self.unclaimed(name)?,
                Source::Field(number) => {
                    if memmem::find(value, b"%3").is_some() {
                        let mut decoded = std::mem::take(&mut self.decoded);
                        decode_extra(value, &mut decoded);
                        self.field(number, &decoded);
                        self.decoded = decoded;
                    } else {
                        self.field(number, value);
                    }
                    count += 1;
                }
            }
        }
        if count > 0 {
            self.record(count)
                .map_err(|error| Fault::Told(error.to_string()))?;
        }
        Ok(())
    }

    /// One line of a tuple file: its value as written, under its field's source.
    fn tuple_line(&mut self, row: &[u8]) -> Result<(), Fault> {
        let mut columns = row.split(|&byte| byte == b'\t');
        let mut taken = [&row[..0]; 5];
        let mut found = 0;
        for column in &mut columns {
            if found < taken.len() {
                taken[found] = column;
            }
            found += 1;
        }
        if found != taken.len() {
            return Err(Fault::Width(found));
        }
        let [location, allele, feature, name, value] = taken;
        self.locate(location)?;
        if value.is_empty() || value == b"-" {
            return Ok(());
        }
        match self.source_of(name)? {
            Source::Core => self.outside += 1,
            Source::Unclaimed => {
                self.outside += 1;
                self.unclaimed(name)?;
            }
            Source::Field(number) => {
                self.start_key(allele, feature);
                self.field(number, value);
                self.record(1)
                    .map_err(|error| Fault::Told(error.to_string()))?;
            }
        }
        Ok(())
    }

    fn block(&mut self, index: usize, block: &[u8]) {
        let (mut number, mut rows) = (0u64, 0u64);
        for row in lines(block) {
            number += 1;
            if row.is_empty() {
                continue;
            }
            self.at = (index, number);
            let outcome = if self.input.kind == Kind::Default && row.starts_with(b"#") {
                let named = memchr_iter(b'\t', row).count() + 1;
                if row.starts_with(b"##") {
                    Ok(())
                } else if named != crate::input::DEFAULT_COLUMNS {
                    Err(Fault::Header(named))
                } else {
                    self.header = Some((index, number, row.to_vec()));
                    Ok(())
                }
            } else {
                rows += 1;
                match self.input.kind {
                    Kind::Default => self.default_row(row),
                    Kind::Tuples => self.tuple_line(row),
                }
            };
            if let Err(fault) = outcome {
                self.fault = Some((index, number, fault));
                return;
            }
        }
        self.counted.push((index, number, rows));
    }
}

/// Of what the threads found, each with its block and its line there, the first in the
/// file's order: the fault the scorer, reading line by line, ends on.
fn first_in_the_file<T>(found: Vec<(usize, u64, T)>) -> Option<(usize, u64, T)> {
    found
        .into_iter()
        .min_by_key(|(block, line, _)| (*block, *line))
}

/// The last in the file's order: the column header the scorer holds at the file's end.
fn last_in_the_file<T>(found: Vec<(usize, u64, T)>) -> Option<(usize, u64, T)> {
    found
        .into_iter()
        .max_by_key(|(block, line, _)| (*block, *line))
}

/// Extracts `input`, the run's input `number`, into its shards.
pub fn extract(run: &Run, input: &Input, number: usize) -> Res<Extracted> {
    let next = AtomicUsize::new(0);
    // No block after the first one found at fault is read: the first fault is the one named.
    let faulty = AtomicUsize::new(usize::MAX);
    let locks: Vec<Mutex<()>> = (0..run.shards()).map(|_| Mutex::new(())).collect();
    let workers = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..run.threads.min(input.blocks()).max(1))
            .map(|_| {
                scope.spawn(|| {
                    let mut worker = Worker {
                        run,
                        input,
                        number,
                        locks: &locks,
                        locator: Locator::new(run),
                        sources: HashMap::new(),
                        shards: vec![Vec::new(); run.shards()],
                        key: Vec::new(),
                        fields: Vec::new(),
                        decoded: Vec::new(),
                        at: (0, 0),
                        outside: 0,
                        unclaimed: HashMap::new(),
                        counted: Vec::new(),
                        header: None,
                        fault: None,
                    };
                    let mut block = Vec::new();
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        if index >= input.blocks() || index > faulty.load(Ordering::Relaxed) {
                            break;
                        }
                        if let Err(error) = input.read_block(index, &mut block) {
                            worker.fault = Some((index, 0, Fault::Told(error.to_string())));
                        } else {
                            worker.block(index, &block);
                        }
                        if worker.fault.is_some() {
                            faulty.fetch_min(index, Ordering::Relaxed);
                            break;
                        }
                    }
                    if worker.fault.is_none() {
                        for shard in 0..worker.shards.len() {
                            if let Err(error) = worker.flush(shard) {
                                worker.fault =
                                    Some((usize::MAX, 0, Fault::Told(error.to_string())));
                                break;
                            }
                        }
                    }
                    worker
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("an extraction thread panicked"))
            .collect::<Vec<_>>()
    });

    let mut lines_of = vec![0u64; input.blocks()];
    let mut rows_of = vec![0u64; input.blocks()];
    let mut headers = Vec::new();
    let mut faults = Vec::new();
    let mut outside = 0;
    let mut unclaimed: HashMap<Vec<u8>, (u64, (usize, u64))> = HashMap::new();
    for worker in workers {
        for (index, lines, rows) in worker.counted {
            lines_of[index] = lines;
            rows_of[index] = rows;
        }
        headers.extend(worker.header);
        faults.extend(worker.fault);
        outside += worker.outside;
        for (name, (times, met)) in worker.unclaimed {
            let found = unclaimed.entry(name).or_insert((0, met));
            found.0 += times;
            found.1 = found.1.min(met);
        }
    }
    let shown = &input.shown;
    let distinct_keys = || {
        Refusal(format!(
            "{shown}: more than {} distinct keys that no source claims",
            run.distinct_keys
        ))
    };
    // The line on which the scorer, reading in the file's order, meets one distinct key too
    // many, whichever threads met the keys: it reaches no fault after that line.
    let overflow = (unclaimed.len() > run.distinct_keys).then(|| {
        let mut met: Vec<(usize, u64)> = unclaimed.values().map(|&(_, met)| met).collect();
        met.sort_unstable();
        met[run.distinct_keys]
    });
    let fault = first_in_the_file(faults)
        .filter(|&(index, line, _)| overflow.is_none_or(|met| (index, line) < met));
    if let Some((index, line, fault)) = fault {
        let number = lines_of.iter().take(index).sum::<u64>() + line;
        return Err(match fault {
            Fault::Width(columns) => {
                let tuples = if input.tuple_file_allowed {
                    ", a tuple file line 5"
                } else {
                    ""
                };
                Refusal(format!(
                    "{shown}: line {number} has {columns} columns; a default-format row has 14{tuples}"
                ))
            }
            Fault::Header(named) => Refusal(format!(
                "{shown}: line {number}, the column header, names {named} columns, not 14"
            )),
            Fault::DistinctKeys => distinct_keys(),
            Fault::Told(message) => Refusal(message),
        });
    }
    if unclaimed.len() > run.distinct_keys {
        return Err(distinct_keys());
    }
    let mut unclaimed: Vec<_> = unclaimed
        .into_iter()
        .map(|(name, (times, _))| (name, times))
        .collect();
    unclaimed.sort();
    let mut rows_before = Vec::with_capacity(rows_of.len());
    let mut lines = 0;
    for rows in rows_of {
        rows_before.push(lines);
        lines += rows;
    }
    Ok(Extracted {
        lines,
        outside,
        unclaimed,
        header: last_in_the_file(headers).map_or_else(|| DEFAULT_HEADER.to_vec(), |found| found.2),
        rows_before,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_fault_and_last_header_of_the_file_are_kept_whichever_thread_found_them() {
        let found = || {
            vec![
                (4, 1, "d"),
                (2, 9, "b"),
                (7, 3, "e"),
                (2, 5, "a"),
                (4, 0, "c"),
            ]
        };
        assert_eq!(first_in_the_file(found()), Some((2, 5, "a")));
        assert_eq!(last_in_the_file(found()), Some((7, 3, "e")));
        assert_eq!(first_in_the_file(Vec::<(usize, u64, ())>::new()), None);
    }
}
