// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A whole reference output read on the host's cores: the file is cut into blocks of whole
//! lines, worker threads read what a pass needs of each block's rows, and the pass takes the
//! blocks' readings in file order, so that the order of the rows decides what it decides
//! under one thread. The lines and the rows are the ones `TextLines` and `iter_rows` give.

use super::text::{not_readable_text, ByteSource, ReadError};
use super::{row_view, RowKeys, RowView, Substitutions};
use crate::classes::FamilySpec;
use crate::error::Failure;
use memchr::memmem::Finder;
use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::num::NonZero;
use std::sync::mpsc::{sync_channel, Receiver};
use std::sync::Mutex;

/// The environment variable that bounds the threads that read rows; unset or empty, the
/// host's available cores are the bound. The thread count changes no byte the tool writes.
pub const THREADS_ENV: &str = "VEP_GOLDEN_CORPUS_THREADS";

/// The bytes a block is cut near.
const BLOCK: usize = 2 << 20;

/// The most threads that read the blocks of one file.
const MOST_THREADS: usize = 32;

/// The bound on the threads that read rows in this run.
pub fn thread_count() -> Result<usize, Failure> {
    thread_bound(std::env::var_os(THREADS_ENV).as_deref())
}

/// The bound a value of the environment variable states.
fn thread_bound(value: Option<&OsStr>) -> Result<usize, Failure> {
    let available = || std::thread::available_parallelism().map_or(1, NonZero::get);
    match value {
        None => Ok(available()),
        Some(value) if value.is_empty() => Ok(available()),
        Some(value) => value
            .to_str()
            .and_then(|text| text.trim().parse::<usize>().ok())
            .filter(|threads| *threads >= 1)
            .ok_or_else(|| {
                Failure::own(format!(
                    "{THREADS_ENV}={}: expected a whole number of threads, 1 or more",
                    value.to_string_lossy()
                ))
            }),
    }
}

/// The lines of a block, each without its line end: a line ends at `\n`, `\r` or `\r\n`,
/// and a last line needs no line end.
pub fn lines(block: &str) -> Lines<'_> {
    Lines { rest: block }
}

/// The iterator `lines` returns.
pub struct Lines<'a> {
    rest: &'a str,
}

impl<'a> Iterator for Lines<'a> {
    type Item = &'a str;

    fn next(&mut self) -> Option<&'a str> {
        if self.rest.is_empty() {
            return None;
        }
        let bytes = self.rest.as_bytes();
        let Some(at) = memchr::memchr2(b'\n', b'\r', bytes) else {
            return Some(std::mem::take(&mut self.rest));
        };
        let after = if bytes[at] == b'\r' && bytes.get(at + 1) == Some(&b'\n') {
            at + 2
        } else {
            at + 1
        };
        let line = &self.rest[..at];
        self.rest = &self.rest[after..];
        Some(line)
    }
}

/// What tells a row of a default-format output from its other lines, and reads one: the
/// dataset the rows belong to, the specs of `--family-columns` and the substitutions of the
/// dataset's input, as `iter_rows` takes them.
pub struct RowReader<'a> {
    dataset: String,
    keys: RowKeys,
    specs: Vec<FamilySpec>,
    substitutions: Option<&'a Substitutions>,
}

impl<'a> RowReader<'a> {
    pub fn new<S: AsRef<str>>(
        dataset: &str,
        family_specs: &[S],
        substitutions: Option<&'a Substitutions>,
    ) -> RowReader<'a> {
        RowReader {
            dataset: dataset.to_string(),
            keys: RowKeys {
                strand: Finder::new("STRAND"),
                flags: Finder::new("FLAGS"),
            },
            specs: family_specs
                .iter()
                .map(|spec| FamilySpec::read(spec.as_ref()))
                .collect(),
            substitutions: substitutions.filter(|held| !held.is_empty()),
        }
    }

    /// The row a line holds; `None` for a line that starts with `#` and for one of fewer
    /// than seven columns.
    pub fn view<'l>(&'l self, line: &'l str) -> Option<RowView<'l>> {
        let tabs = memchr::memchr_iter(b'\t', line.as_bytes()).take(6).count();
        if tabs < 6 || line.starts_with('#') {
            return None;
        }
        Some(row_view(
            line,
            &self.dataset,
            &self.keys,
            &self.specs,
            self.substitutions,
        ))
    }
}

/// A file as blocks of whole lines.
struct Blocks {
    source: ByteSource,
    /// The bytes read that no block has taken: the start of a line.
    carry: Vec<u8>,
    ended: bool,
}

impl Blocks {
    /// The next block as text, bytes outside UTF-8 read as `TextLines` reads them; `None`
    /// after the last. A block ends behind a line end, the last one where the file ends.
    fn next(&mut self) -> Result<Option<String>, ReadError> {
        let mut bytes = std::mem::take(&mut self.carry);
        let mut filled = bytes.len();
        let mut searched = 0;
        let cut = loop {
            if self.ended {
                break filled;
            }
            if filled >= BLOCK {
                if let Some(cut) = line_end_before(&bytes[..filled], searched) {
                    break cut;
                }
                searched = filled;
            }
            if bytes.len() == filled {
                bytes.resize(filled + BLOCK, 0);
            }
            let count = self.source.read(&mut bytes[filled..])?;
            filled += count;
            self.ended = count == 0;
        };
        if filled == 0 {
            return Ok(None);
        }
        self.carry = bytes[cut..filled].to_vec();
        bytes.truncate(cut);
        Ok(Some(match String::from_utf8(bytes) {
            Ok(text) => text,
            Err(error) => String::from_utf8_lossy(error.as_bytes()).into_owned(),
        }))
    }
}

/// The place behind the last line end of `bytes` at or after `from` that is known to be
/// one: behind a `\n`, or behind a `\r` that a byte other than `\n` follows.
fn line_end_before(bytes: &[u8], from: usize) -> Option<usize> {
    let from = from.saturating_sub(1);
    if let Some(at) = memchr::memrchr(b'\n', &bytes[from..]) {
        return Some(from + at + 1);
    }
    let known = &bytes[from..bytes.len() - 1];
    memchr::memrchr(b'\r', known).map(|at| from + at + 1)
}

/// Reads `path` a block of whole lines at a time: `read` takes what the pass needs of a
/// block, on one of at most `threads` threads, and `take` is given each block with its
/// reading in file order on the caller's thread. A file that cannot be opened or read is
/// refused as `iter_rows` refuses it, once the blocks before the failure are taken; a
/// failure `take` returns ends the pass.
pub fn read_blocks<T: Send>(
    path: &str,
    threads: usize,
    read: impl Fn(&str) -> T + Sync,
    mut take: impl FnMut(&str, T) -> Result<(), Failure>,
) -> Result<(), Failure> {
    let source = ByteSource::open(path).map_err(|error| not_readable_text(path, &error))?;
    let mut blocks = Blocks {
        source,
        carry: Vec::new(),
        ended: false,
    };
    if threads <= 1 {
        loop {
            match blocks.next() {
                Ok(Some(block)) => take(&block, read(&block))?,
                Ok(None) => return Ok(()),
                Err(error) => return Err(not_readable_text(path, &error)),
            }
        }
    }
    // More threads than these wait on the one that reads the file. The blocks read ahead
    // of the one `take` is given are bounded, so that a pass holds two blocks a thread
    // however fast the file is read.
    let threads = threads.min(MOST_THREADS);
    let window = 2 * threads;
    let (permit_out, permit_in) = sync_channel::<()>(window);
    let (block_out, block_in) = sync_channel::<(usize, String)>(threads);
    let (done_out, done_in) = sync_channel::<(usize, String, T)>(threads);
    let block_in: Mutex<Receiver<(usize, String)>> = Mutex::new(block_in);
    std::thread::scope(|scope| {
        let reader = scope.spawn(move || -> Result<(), ReadError> {
            let mut index = 0;
            while permit_out.send(()).is_ok() {
                let Some(block) = blocks.next()? else {
                    break;
                };
                if block_out.send((index, block)).is_err() {
                    break;
                }
                index += 1;
            }
            Ok(())
        });
        for _ in 0..threads {
            let done_out = done_out.clone();
            let (block_in, read) = (&block_in, &read);
            scope.spawn(move || loop {
                let next = block_in.lock().expect("the block queue").recv();
                let Ok((index, block)) = next else {
                    break;
                };
                let reading = read(&block);
                if done_out.send((index, block, reading)).is_err() {
                    break;
                }
            });
        }
        drop(done_out);
        let mut waiting: BTreeMap<usize, (String, T)> = BTreeMap::new();
        let mut next = 0;
        let mut outcome = Ok(());
        'blocks: for (index, block, reading) in &done_in {
            waiting.insert(index, (block, reading));
            while let Some((block, reading)) = waiting.remove(&next) {
                if let Err(failure) = take(&block, reading) {
                    outcome = Err(failure);
                    break 'blocks;
                }
                next += 1;
                let _ = permit_in.recv();
            }
        }
        // A thread that waits on one of the two ends or on the block queue finishes once
        // both ends are closed and the queue is empty.
        drop(done_in);
        drop(permit_in);
        while block_in.lock().expect("the block queue").recv().is_ok() {}
        let read_to_the_end = reader
            .join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic));
        outcome?;
        read_to_the_end.map_err(|error| not_readable_text(path, &error))
    })
}

#[cfg(test)]
#[path = "blocks_tests.rs"]
mod tests;
