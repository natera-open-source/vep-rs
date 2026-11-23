// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of the input as tuple lines: the input read in blocks of whole lines, each
//! block written as tuples on a worker thread, the blocks' tuples written in input order,
//! and the counts the summary line states.

use crate::cells::{decode_cell, Decoded, Scratch};
use crate::columns::{header_cells, line_text, Placed};
use crate::map::Half;
use crate::pyfmt::shown;
use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, Read, Write};
use std::num::NonZero;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Condvar, Mutex};

/// The environment variable that bounds the worker threads; unset or empty, the host's
/// available cores are the bound. An environment variable and not an option: the command
/// line names the run, and the thread count changes no byte of its outputs.
pub const THREADS_ENV: &str = "VEP_FASTVEP_PLUGIN_TUPLES_THREADS";

/// The bytes of input a block of a run holds before it is cut at its last line end.
const BLOCK: usize = 1 << 18;

/// The worker threads of a run over `input`: the bound, and no more than a thread a block
/// of a regular file, so a short input starts few threads.
pub fn worker_threads(input: &File) -> Result<usize, String> {
    let blocks = match input.metadata() {
        Ok(file) if file.is_file() => 1 + file.len() / BLOCK as u64,
        _ => u64::MAX,
    };
    let bound = thread_bound(std::env::var_os(THREADS_ENV).as_deref())?;
    Ok(bound.min(usize::try_from(blocks).unwrap_or(usize::MAX)))
}

/// The worker-thread bound a value of the environment variable states; unset or empty,
/// the host's available cores.
fn thread_bound(value: Option<&OsStr>) -> Result<usize, String> {
    let available = || std::thread::available_parallelism().map_or(1, NonZero::get);
    match value {
        None => Ok(available()),
        Some(value) if value.is_empty() => Ok(available()),
        Some(value) => value
            .to_str()
            .and_then(|text| text.trim().parse::<usize>().ok())
            .filter(|threads| *threads >= 1)
            .ok_or_else(|| {
                format!(
                    "{THREADS_ENV}={}: expected a whole number of threads, 1 or more",
                    value.to_string_lossy()
                )
            }),
    }
}

/// The input file, read in blocks that end at a line end.
pub struct Input {
    file: File,
    /// Bytes read that no block has handed out.
    carry: Vec<u8>,
    ended: bool,
    /// The bytes a block holds before it is cut at its last line end.
    block: usize,
}

impl Input {
    pub fn new(file: File) -> Input {
        Input::with_blocks_of(file, BLOCK)
    }

    /// The input read in blocks of at least `block` bytes.
    fn with_blocks_of(file: File, block: usize) -> Input {
        Input {
            file,
            carry: Vec::new(),
            ended: false,
            block,
        }
    }

    /// Appends up to one block of the file to `buffer`; the end of the file is remembered.
    fn read_more(&mut self, buffer: &mut Vec<u8>) -> io::Result<()> {
        let old = buffer.len();
        buffer.resize(old + self.block, 0);
        let read = loop {
            match self.file.read(&mut buffer[old..]) {
                Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
                other => break other,
            }
        };
        buffer.truncate(old + *read.as_ref().unwrap_or(&0));
        self.ended = read? == 0;
        Ok(())
    }

    /// Reads the lines before the first row and returns the cells of the column header:
    /// the last of them that starts with one `#`. None for an input without one.
    pub fn header(&mut self) -> io::Result<Option<Vec<Vec<u8>>>> {
        let mut header = None;
        let mut carry = std::mem::take(&mut self.carry);
        // The line being read starts at `at`; no line feed lies between it and `searched`.
        let (mut at, mut searched) = (0, 0);
        loop {
            if carry.get(at).is_some_and(|first| *first != b'#') {
                break;
            }
            let end = match memchr::memchr(b'\n', &carry[searched..]) {
                Some(found) => searched + found + 1,
                None if self.ended => carry.len(),
                None => {
                    carry.drain(..at);
                    (at, searched) = (0, carry.len());
                    self.read_more(&mut carry)?;
                    continue;
                }
            };
            if end == at {
                break;
            }
            if !carry[at..end].starts_with(b"##") {
                header = Some(header_cells(&carry[at..end]));
            }
            (at, searched) = (end, end);
        }
        carry.drain(..at);
        self.carry = carry;
        Ok(header)
    }

    /// The next block of whole lines, the last line of the file with or without its line
    /// feed. None at the end of the file.
    fn next_block(&mut self) -> io::Result<Option<Vec<u8>>> {
        let mut block = std::mem::take(&mut self.carry);
        // No line feed lies before `searched`, unless the block is still under its size.
        let mut searched = 0;
        loop {
            if self.ended {
                return Ok((!block.is_empty()).then_some(block));
            }
            if block.len() < self.block {
                self.read_more(&mut block)?;
                continue;
            }
            if let Some(found) = memchr::memrchr(b'\n', &block[searched..]) {
                self.carry = block.split_off(searched + found + 1);
                return Ok(Some(block));
            }
            searched = block.len();
            self.read_more(&mut block)?;
        }
    }
}

/// What the summary line counts.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub rows: u64,
    /// Rows that end before the last mapped column.
    pub short: u64,
    pub deselected: u64,
    pub unshaped: u64,
    pub several: u64,
    /// The tuples of each kept field, in the order of `Placed::fields`.
    pub by_field: Vec<u64>,
    /// The dropped parts and halves: the column's place among the mapped columns, the
    /// part and the half, with how often each was dropped.
    pub dropped: BTreeMap<(usize, u64, Half), u64>,
}

impl Counts {
    fn new(fields: usize) -> Counts {
        Counts {
            by_field: vec![0; fields],
            ..Counts::default()
        }
    }

    fn add(&mut self, other: Counts) {
        self.rows += other.rows;
        self.short += other.short;
        self.deselected += other.deselected;
        self.unshaped += other.unshaped;
        self.several += other.several;
        for (total, count) in self.by_field.iter_mut().zip(other.by_field) {
            *total += count;
        }
        for (part, count) in other.dropped {
            *self.dropped.entry(part).or_default() += count;
        }
    }
}

/// The last cell a worker decoded in one column: a source's cell repeats on every row of
/// its allele, so the decoding is kept.
#[derive(Default)]
struct Memo {
    cell: Vec<u8>,
    decoded: Decoded,
    filled: bool,
}

/// One worker's state: its counts and what it keeps between rows.
struct Worker<'a> {
    placed: &'a Placed,
    counts: Counts,
    memos: Vec<Memo>,
    scratch: Scratch,
    tabs: Vec<usize>,
}

impl<'a> Worker<'a> {
    fn new(placed: &'a Placed) -> Worker<'a> {
        Worker {
            placed,
            counts: Counts::new(placed.fields.len()),
            memos: placed.columns.iter().map(|_| Memo::default()).collect(),
            scratch: Scratch::default(),
            tabs: Vec::new(),
        }
    }

    /// Appends the tuples of every row of `block` to `out`. A line that starts with `#`
    /// is no row.
    fn write_block(&mut self, block: &[u8], out: &mut Vec<u8>) {
        let mut start = 0;
        while start < block.len() {
            let end =
                memchr::memchr(b'\n', &block[start..]).map_or(block.len(), |at| start + at + 1);
            if block[start] != b'#' {
                self.write_row(&block[start..end], out);
            }
            start = end;
        }
    }

    fn write_row(&mut self, line: &[u8], out: &mut Vec<u8>) {
        self.counts.rows += 1;
        let text = line_text(line);
        let text: &[u8] = &text;
        self.tabs.clear();
        self.tabs.extend(memchr::memchr_iter(b'\t', text));
        let tabs = &self.tabs;
        let cells = tabs.len() + 1;
        if cells < self.placed.width {
            self.counts.short += 1;
        }
        if cells < 5 {
            return;
        }
        let cell = |index: usize| {
            let start = if index == 0 { 0 } else { tabs[index - 1] + 1 };
            let end = tabs.get(index).copied().unwrap_or(text.len());
            &text[start..end]
        };
        for (order, column) in self.placed.columns.iter().enumerate() {
            if column.index >= cells {
                continue;
            }
            let value = cell(column.index);
            if value.is_empty() || value == b"-" {
                continue;
            }
            let memo = &mut self.memos[order];
            if !memo.filled || memo.cell != value {
                decode_cell(value, column, &mut memo.decoded, &mut self.scratch);
                memo.cell.clear();
                memo.cell.extend_from_slice(value);
                memo.filled = true;
            }
            for (slot, tail) in &memo.decoded.pairs {
                for key in [1, 2, 4] {
                    out.extend_from_slice(cell(key));
                    out.push(b'\t');
                }
                out.extend_from_slice(tail);
                self.counts.by_field[*slot] += 1;
            }
            for (part, half) in &memo.decoded.dropped {
                *self
                    .counts
                    .dropped
                    .entry((order, *part, *half))
                    .or_default() += 1;
            }
            self.counts.deselected += memo.decoded.deselected;
            self.counts.unshaped += memo.decoded.unshaped;
            self.counts.several += u64::from(memo.decoded.several);
        }
    }
}

/// The input and the number of the block it hands out next.
struct Source {
    input: Input,
    next: u64,
    done: bool,
    error: Option<io::Error>,
}

/// The tuple file, the number of the block it writes next and the blocks that wait for it.
struct Sink {
    out: File,
    next: u64,
    waiting: BTreeMap<u64, Vec<u8>>,
    error: Option<io::Error>,
}

/// Writes the rows of `input` to `out` as tuples on `threads` threads and returns the
/// counts. The tuples of a block are written when every block before it is written, and a
/// worker takes no block further than one block a thread past the last one written, so
/// the blocks in memory are bounded by the threads.
pub fn write_tuples(
    input: Input,
    placed: &Placed,
    out: File,
    threads: usize,
) -> io::Result<Counts> {
    let window = threads as u64;
    let handed = AtomicU64::new(0);
    let source = Mutex::new(Source {
        input,
        next: 0,
        done: false,
        error: None,
    });
    let sink = Mutex::new(Sink {
        out,
        next: 0,
        waiting: BTreeMap::new(),
        error: None,
    });
    let room = Condvar::new();
    let work = || {
        let mut worker = Worker::new(placed);
        loop {
            {
                let mut sink = sink.lock().expect("the sink");
                while sink.error.is_none() && handed.load(Ordering::Relaxed) >= sink.next + window {
                    sink = room.wait(sink).expect("the sink");
                }
                if sink.error.is_some() {
                    break;
                }
            }
            let (number, block) = {
                let mut source = source.lock().expect("the source");
                if source.done {
                    break;
                }
                match source.input.next_block() {
                    Ok(Some(block)) => {
                        source.next += 1;
                        handed.store(source.next, Ordering::Relaxed);
                        (source.next - 1, block)
                    }
                    Ok(None) => {
                        source.done = true;
                        break;
                    }
                    Err(error) => {
                        source.done = true;
                        source.error = Some(error);
                        break;
                    }
                }
            };
            let mut tuples = Vec::with_capacity(block.len());
            worker.write_block(&block, &mut tuples);
            drop(block);
            let mut sink = sink.lock().expect("the sink");
            sink.waiting.insert(number, tuples);
            while sink.error.is_none() {
                let next = sink.next;
                let Some(tuples) = sink.waiting.remove(&next) else {
                    break;
                };
                match sink.out.write_all(&tuples) {
                    Ok(()) => sink.next += 1,
                    Err(error) => sink.error = Some(error),
                }
            }
            room.notify_all();
        }
        worker.counts
    };
    let mut counts = Counts::new(placed.fields.len());
    std::thread::scope(|scope| {
        let others: Vec<_> = (1..threads).map(|_| scope.spawn(work)).collect();
        counts.add(work());
        for other in others {
            counts.add(
                other
                    .join()
                    .unwrap_or_else(|panic| std::panic::resume_unwind(panic)),
            );
        }
    });
    if let Some(error) = source.into_inner().expect("the source").error {
        return Err(error);
    }
    if let Some(error) = sink.into_inner().expect("the sink").error {
        return Err(error);
    }
    Ok(counts)
}

/// The summary line: the tuples written, the rows read and the rows that end before the
/// last mapped column, the tuples of each kept field, the mapped columns, the dropped
/// parts and halves by column, the values `--fields` left out, the parts a `word(score)`
/// row could not split, the cells with several entries and the `FV_` columns the map has
/// no row for.
pub fn summary(placed: &Placed, counts: &Counts) -> String {
    // A column is named once, so the order of names, then parts, then halves is the order
    // of `(name, part number, part text)`.
    let mut dropped: Vec<(&str, u64, Half, u64)> = counts
        .dropped
        .iter()
        .map(|((order, part, half), count)| {
            (placed.columns[*order].name.as_str(), *part, *half, *count)
        })
        .collect();
    dropped.sort();
    let dropped_total: u64 = dropped.iter().map(|(_, _, _, count)| count).sum();
    let dropped: Vec<String> = dropped
        .iter()
        .map(|(name, part, half, count)| {
            let half = match half {
                Half::Whole => "",
                Half::Word => ".word",
                Half::Score => ".score",
            };
            format!("{name}:{part}{half}={count}")
        })
        .collect();
    let dropped = if dropped.is_empty() {
        String::new()
    } else {
        format!(" ({})", dropped.join(", "))
    };
    let by_field: Vec<String> = placed
        .fields
        .iter()
        .zip(&counts.by_field)
        .map(|(field, count)| format!("{field}={count}"))
        .collect();
    let by_field = if by_field.is_empty() {
        "none".to_string()
    } else {
        by_field.join(", ")
    };
    let mapped: Vec<&str> = placed
        .columns
        .iter()
        .map(|column| column.name.as_str())
        .collect();
    let unmapped: Vec<String> = placed.unmapped.iter().map(|name| shown(name)).collect();
    let unmapped = if unmapped.is_empty() {
        "none".to_string()
    } else {
        unmapped.join(", ")
    };
    format!(
        "fastvep_plugin_tuples: wrote {} tuples from {} rows ({} too short); by field: {by_field}; \
         mapped columns: {}; dropped parts with no Ensembl VEP field: {dropped_total}{dropped}; \
         values outside --fields: {}; parts not written word(score): {}; \
         cells with several entries: {}; unmapped columns: {unmapped}",
        counts.by_field.iter().sum::<u64>(),
        counts.rows,
        counts.short,
        mapped.join(", "),
        counts.deselected,
        counts.unshaped,
        counts.several,
    )
}

#[cfg(test)]
#[path = "rows_blocks_tests.rs"]
mod blocks_tests;
