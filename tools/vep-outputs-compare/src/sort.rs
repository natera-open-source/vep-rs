// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `sort -u` under `LC_ALL=C` in process: the keys of one side ordered bytewise with
//! duplicates dropped, written one per line. Keys past the memory bound are sorted and
//! spilled as runs under the scratch directory and merged on finish.

use crate::input::push_collapsed_key;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Bytes of keys held in memory before a run is spilled.
pub const DEFAULT_MEMORY_BOUND: usize = 1 << 30;

/// The environment variable under which a debug build exchanges the last two keys of every
/// sort, so that a test of the whole tool can hand the merge keys that do not ascend. A
/// release build never reads it.
#[cfg(debug_assertions)]
const EXCHANGE_LAST_KEYS_ENV: &str = "VEP_OUTPUTS_COMPARE_DEBUG_EXCHANGE_LAST_KEYS";

#[cfg(debug_assertions)]
fn exchange_last_keys_where_asked(bounds: &mut [(usize, usize)]) {
    let asked = std::env::var_os(EXCHANGE_LAST_KEYS_ENV).is_some_and(|value| !value.is_empty());
    if let ([.., last_but_one, last], true) = (bounds, asked) {
        std::mem::swap(last_but_one, last);
    }
}

pub struct KeySorter {
    arena: Vec<u8>,
    bounds: Vec<(usize, usize)>,
    memory_bound: usize,
    runs: Vec<PathBuf>,
    scratch: PathBuf,
}

impl KeySorter {
    /// `scratch` holds the spilled runs, named after `tag`.
    pub fn new(scratch: &Path, tag: &str, memory_bound: usize) -> KeySorter {
        KeySorter {
            arena: Vec::new(),
            bounds: Vec::new(),
            memory_bound,
            runs: Vec::new(),
            scratch: scratch.join(tag),
        }
    }

    /// Adds one key, given without its newline.
    pub fn push(&mut self, key: &[u8]) -> io::Result<()> {
        let start = self.arena.len();
        self.arena.extend_from_slice(key);
        self.bounds.push((start, key.len()));
        if self.arena.len() >= self.memory_bound {
            self.spill()?;
        }
        Ok(())
    }

    fn sorted_unique(&mut self) -> Vec<(usize, usize)> {
        let arena = &self.arena;
        let mut bounds = std::mem::take(&mut self.bounds);
        bounds.sort_unstable_by(|a, b| arena[a.0..a.0 + a.1].cmp(&arena[b.0..b.0 + b.1]));
        bounds.dedup_by(|a, b| arena[a.0..a.0 + a.1] == arena[b.0..b.0 + b.1]);
        #[cfg(debug_assertions)]
        exchange_last_keys_where_asked(&mut bounds);
        bounds
    }

    fn write_keys(&self, bounds: &[(usize, usize)], out: &mut impl Write) -> io::Result<()> {
        for &(start, len) in bounds {
            out.write_all(&self.arena[start..start + len])?;
            out.write_all(b"\n")?;
        }
        Ok(())
    }

    fn spill(&mut self) -> io::Result<()> {
        let bounds = self.sorted_unique();
        let path = PathBuf::from(format!("{}.run{}", self.scratch.display(), self.runs.len()));
        let mut out = BufWriter::new(File::create(&path)?);
        self.write_keys(&bounds, &mut out)?;
        out.flush()?;
        self.runs.push(path);
        self.arena.clear();
        Ok(())
    }

    /// Writes the sorted unique keys to `out_path`, one per line.
    pub fn finish(mut self, out_path: &Path) -> io::Result<()> {
        let mut out = BufWriter::with_capacity(1 << 20, File::create(out_path)?);
        if self.runs.is_empty() {
            let bounds = self.sorted_unique();
            self.write_keys(&bounds, &mut out)?;
            return out.flush();
        }
        if !self.bounds.is_empty() {
            self.spill()?;
        }
        let runs = std::mem::take(&mut self.runs);
        let result = merge_runs(&runs, &mut out);
        for run in &runs {
            let _ = std::fs::remove_file(run);
        }
        result?;
        out.flush()
    }
}

/// Sorted unique keys held in memory.
pub struct SortedKeys {
    arena: Vec<u8>,
    bounds: Vec<(usize, usize)>,
}

/// The start and length of each line of `arena`; the last needs no newline.
fn line_bounds(arena: &[u8]) -> Vec<(usize, usize)> {
    let mut bounds = Vec::new();
    let mut start = 0;
    for end in memchr::memchr_iter(b'\n', arena) {
        bounds.push((start, end - start));
        start = end + 1;
    }
    if start < arena.len() {
        bounds.push((start, arena.len() - start));
    }
    bounds
}

impl SortedKeys {
    /// The keys of `arena`, a newline after each, ordered bytewise with duplicates dropped.
    pub fn from_lines(arena: Vec<u8>) -> SortedKeys {
        let mut bounds = line_bounds(&arena);
        bounds.sort_unstable_by(|a, b| arena[a.0..a.0 + a.1].cmp(&arena[b.0..b.0 + b.1]));
        bounds.dedup_by(|a, b| arena[a.0..a.0 + a.1] == arena[b.0..b.0 + b.1]);
        #[cfg(debug_assertions)]
        exchange_last_keys_where_asked(&mut bounds);
        SortedKeys { arena, bounds }
    }

    /// The keys of `arena`, a newline after each, in the order they are written in.
    #[cfg(test)]
    pub fn in_the_order_written(arena: Vec<u8>) -> SortedKeys {
        let bounds = line_bounds(&arena);
        SortedKeys { arena, bounds }
    }

    /// The first and the last key, None where there is no key.
    pub fn ends(&self) -> Option<(Vec<u8>, Vec<u8>)> {
        let last = self.len().checked_sub(1)?;
        Some((self.get(0).to_vec(), self.get(last).to_vec()))
    }

    /// Writes the keys in the order they are held in, one to a line.
    pub fn write_lines(&self, out: &mut impl Write) -> io::Result<()> {
        for &(start, len) in &self.bounds {
            out.write_all(&self.arena[start..start + len])?;
            out.write_all(b"\n")?;
        }
        Ok(())
    }

    /// The context-collapsed forms of these keys, sorted and unique.
    pub fn collapsed(&self) -> SortedKeys {
        let mut arena = Vec::with_capacity(self.arena.len());
        for &(start, len) in &self.bounds {
            push_collapsed_key(&self.arena[start..start + len], &mut arena);
            arena.push(b'\n');
        }
        SortedKeys::from_lines(arena)
    }

    pub fn len(&self) -> usize {
        self.bounds.len()
    }

    pub fn get(&self, index: usize) -> &[u8] {
        let (start, len) = self.bounds[index];
        &self.arena[start..start + len]
    }
}

#[derive(PartialEq, Eq)]
struct Head {
    line: Vec<u8>,
    run: usize,
}

impl PartialOrd for Head {
    fn partial_cmp(&self, other: &Head) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Head {
    /// Reversed, so the heap's maximum is the smallest line.
    fn cmp(&self, other: &Head) -> Ordering {
        other
            .line
            .cmp(&self.line)
            .then_with(|| other.run.cmp(&self.run))
    }
}

fn next_line(reader: &mut BufReader<File>) -> io::Result<Option<Vec<u8>>> {
    let mut line = Vec::new();
    let n = reader.read_until(b'\n', &mut line)?;
    if n == 0 {
        return Ok(None);
    }
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    Ok(Some(line))
}

/// Merges sorted unique runs into one sorted unique stream.
fn merge_runs(runs: &[PathBuf], out: &mut impl Write) -> io::Result<()> {
    let mut readers: Vec<BufReader<File>> = runs
        .iter()
        .map(|path| File::open(path).map(|f| BufReader::with_capacity(1 << 20, f)))
        .collect::<io::Result<_>>()?;
    let mut heap = BinaryHeap::new();
    for (run, reader) in readers.iter_mut().enumerate() {
        if let Some(line) = next_line(reader)? {
            heap.push(Head { line, run });
        }
    }
    let mut previous: Option<Vec<u8>> = None;
    while let Some(Head { line, run }) = heap.pop() {
        if previous.as_ref() != Some(&line) {
            out.write_all(&line)?;
            out.write_all(b"\n")?;
            previous = Some(line);
        }
        if let Some(next) = next_line(&mut readers[run])? {
            heap.push(Head { line: next, run });
        }
    }
    Ok(())
}
