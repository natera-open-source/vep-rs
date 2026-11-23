// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Records sorted on disk: a writer that spills its buffer as a sorted run
//! each time the buffer fills, and a merge that reads the runs back as one
//! stream in key order.
//!
//! A record is a key of two unsigned 64-bit numbers, compared first then
//! second, and a payload of bytes. The writer holds one buffer of records,
//! never more than its capacity in bytes unless one record alone is larger;
//! the merge holds one read buffer and one record a run. Records of equal
//! keys come back in the order they were pushed.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

/// The key records are ordered by.
pub(crate) type Key = (u64, u64);

/// The bytes a run reader buffers.
const READ_BUFFER_BYTES: usize = 256 << 10;

/// A buffered record: its key and where its payload lies in the buffer.
struct Entry {
    key: Key,
    at: u32,
    len: u32,
}

/// Takes records in any order and writes them as sorted runs under a
/// directory.
pub(crate) struct RunWriter {
    dir: PathBuf,
    label: String,
    capacity: usize,
    entries: Vec<Entry>,
    payloads: Vec<u8>,
    runs: Runs,
}

impl RunWriter {
    /// A writer whose runs are `<dir>/<label>.<n>.run`, spilling a run each
    /// time its buffered records pass `capacity` bytes.
    pub(crate) fn new(dir: &Path, label: &str, capacity: usize) -> Self {
        Self {
            dir: dir.to_path_buf(),
            label: label.to_string(),
            // A payload's place in the buffer is a 32-bit offset.
            capacity: capacity.min(u32::MAX as usize),
            entries: Vec::new(),
            payloads: Vec::new(),
            runs: Runs::default(),
        }
    }

    /// The bytes the buffered records take.
    fn buffered(&self) -> usize {
        self.entries.len() * std::mem::size_of::<Entry>() + self.payloads.len()
    }

    /// Adds one record, first spilling the buffer when the record does not
    /// fit beside what it holds.
    pub(crate) fn push(&mut self, key: Key, payload: &[u8]) -> Result<()> {
        let record = std::mem::size_of::<Entry>() + payload.len();
        if !self.entries.is_empty() && self.buffered() + record > self.capacity {
            self.spill()?;
        }
        let at = u32::try_from(self.payloads.len()).context("a run buffer past 4 GiB")?;
        let len = u32::try_from(payload.len()).context("a record past 4 GiB")?;
        self.payloads.extend_from_slice(payload);
        self.entries.push(Entry { key, at, len });
        Ok(())
    }

    /// Writes the buffered records as one run, in key order, and empties the
    /// buffer.
    fn spill(&mut self) -> Result<()> {
        if self.entries.is_empty() {
            return Ok(());
        }
        // Stable, so records of one key keep the order they were pushed in.
        self.entries.sort_by_key(|entry| entry.key);
        let path = self
            .dir
            .join(format!("{}.{}.run", self.label, self.runs.paths.len()));
        let file =
            File::create(&path).with_context(|| format!("Failed to create {}", path.display()))?;
        // Listed before it is written: a run that fails part way is removed
        // with the others.
        self.runs.paths.push(path.clone());
        let mut writer = BufWriter::with_capacity(1 << 20, file);
        let mut written = 0u64;
        for entry in &self.entries {
            let payload = &self.payloads[entry.at as usize..(entry.at + entry.len) as usize];
            writer
                .write_all(&entry.key.0.to_le_bytes())
                .and_then(|()| writer.write_all(&entry.key.1.to_le_bytes()))
                .and_then(|()| writer.write_all(&entry.len.to_le_bytes()))
                .and_then(|()| writer.write_all(payload))
                .with_context(|| format!("Failed to write {}", path.display()))?;
            written += HEADER_BYTES as u64 + u64::from(entry.len);
        }
        writer
            .flush()
            .with_context(|| format!("Failed to write {}", path.display()))?;
        self.runs.records += self.entries.len() as u64;
        self.runs.bytes += written;
        self.entries.clear();
        self.payloads.clear();
        Ok(())
    }

    /// Spills what the buffer still holds and returns the runs.
    pub(crate) fn finish(mut self) -> Result<Runs> {
        self.spill()?;
        Ok(std::mem::take(&mut self.runs))
    }
}

/// The bytes before a record's payload in a run: its key and its payload's
/// length.
const HEADER_BYTES: usize = 20;

/// The sorted runs a [`RunWriter`] wrote. The files are removed when the
/// value is dropped.
#[derive(Default)]
pub(crate) struct Runs {
    paths: Vec<PathBuf>,
    /// The records of every run.
    pub(crate) records: u64,
    /// The bytes of every run.
    pub(crate) bytes: u64,
}

impl Runs {
    /// The number of run files.
    pub(crate) fn count(&self) -> usize {
        self.paths.len()
    }

    /// Opens every run for one read in key order.
    pub(crate) fn merge(&self) -> Result<Merge> {
        let mut readers = Vec::with_capacity(self.paths.len());
        let mut heap = BinaryHeap::with_capacity(self.paths.len());
        for (run, path) in self.paths.iter().enumerate() {
            let file =
                File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
            let mut reader = RunReader {
                path: path.clone(),
                reader: BufReader::with_capacity(READ_BUFFER_BYTES, file),
                key: (0, 0),
                payload: Vec::new(),
            };
            if reader.read_next()? {
                heap.push(Reverse((reader.key, run)));
            }
            readers.push(reader);
        }
        Ok(Merge {
            readers,
            heap,
            head: None,
        })
    }
}

impl Drop for Runs {
    fn drop(&mut self) {
        for path in &self.paths {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// One run being read: its next record.
struct RunReader {
    path: PathBuf,
    reader: BufReader<File>,
    key: Key,
    payload: Vec<u8>,
}

impl RunReader {
    /// Reads the run's next record; `false` at its end.
    fn read_next(&mut self) -> Result<bool> {
        let failed = || format!("Failed to read {}", self.path.display());
        if self.reader.fill_buf().with_context(failed)?.is_empty() {
            return Ok(false);
        }
        let mut header = [0u8; HEADER_BYTES];
        self.reader.read_exact(&mut header).with_context(failed)?;
        let number = |at: usize| u64::from_le_bytes(header[at..at + 8].try_into().expect("8"));
        self.key = (number(0), number(8));
        let len = u32::from_le_bytes(header[16..].try_into().expect("4")) as usize;
        self.payload.resize(len, 0);
        self.reader
            .read_exact(&mut self.payload)
            .with_context(failed)?;
        Ok(true)
    }
}

/// The records of every run as one stream in key order, runs of equal keys
/// in the order the runs were written.
pub(crate) struct Merge {
    readers: Vec<RunReader>,
    heap: BinaryHeap<Reverse<(Key, usize)>>,
    /// The run whose record is the stream's head.
    head: Option<usize>,
}

impl Merge {
    /// The stream's next record, its key and its payload, without taking it;
    /// `None` at the end of the stream.
    pub(crate) fn head(&mut self) -> Option<(Key, &[u8])> {
        if self.head.is_none() {
            self.head = self.heap.pop().map(|Reverse((_, run))| run);
        }
        let reader = &self.readers[self.head?];
        Some((reader.key, &reader.payload))
    }

    /// Takes the head off the stream.
    pub(crate) fn advance(&mut self) -> Result<()> {
        if self.head.is_none() {
            self.head();
        }
        if let Some(run) = self.head.take() {
            if self.readers[run].read_next()? {
                self.heap.push(Reverse((self.readers[run].key, run)));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every record of the merged stream of `runs`.
    fn merged(runs: &Runs) -> Vec<(Key, String)> {
        let mut merge = runs.merge().unwrap();
        let mut records = Vec::new();
        while let Some((key, payload)) = merge.head() {
            records.push((key, String::from_utf8(payload.to_vec()).unwrap()));
            merge.advance().unwrap();
        }
        records
    }

    /// Records pushed in no order through a buffer that holds three of them
    /// come back in key order from several runs, those of one key in the
    /// order they were pushed, an empty payload among them; the counts are
    /// those of the records and of the run files' bytes, and the files go
    /// when the runs are dropped.
    #[test]
    fn records_spilled_as_several_runs_merge_in_key_order() {
        let dir = tempfile::tempdir().unwrap();
        let capacity = 3 * (std::mem::size_of::<Entry>() + 8);
        let mut writer = RunWriter::new(dir.path(), "names", capacity);
        let pushed: Vec<(Key, String)> = (0..40u64)
            .map(|n| (((n * 7) % 13, (n * 5) % 3), format!("record{n:02}")))
            .chain([((4, 1), String::new()), ((0, 0), "last".to_string())])
            .collect();
        for (key, payload) in &pushed {
            writer.push(*key, payload.as_bytes()).unwrap();
        }
        let runs = writer.finish().unwrap();
        assert_eq!(runs.count(), 14);
        assert_eq!(runs.records, 42);
        let on_disk: u64 = std::fs::read_dir(dir.path())
            .unwrap()
            .map(|entry| entry.unwrap().metadata().unwrap().len())
            .sum();
        assert_eq!(runs.bytes, on_disk);
        assert_eq!(
            runs.bytes,
            42 * HEADER_BYTES as u64 + 40 * 8 + "last".len() as u64
        );

        let mut expected = pushed.clone();
        expected.sort_by_key(|(key, _)| *key);
        assert_eq!(merged(&runs), expected);
        // A second merge reads the same stream.
        assert_eq!(merged(&runs), expected);

        drop(runs);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    /// A writer given nothing writes no run, whose merge is an empty stream;
    /// a record larger than the buffer is a run of its own.
    #[test]
    fn no_record_is_no_run_and_an_oversized_record_is_its_own_run() {
        let dir = tempfile::tempdir().unwrap();
        let runs = RunWriter::new(dir.path(), "none", 64).finish().unwrap();
        assert_eq!((runs.count(), runs.records, runs.bytes), (0, 0, 0));
        assert!(merged(&runs).is_empty());

        let mut writer = RunWriter::new(dir.path(), "large", 64);
        let large = "x".repeat(500);
        writer.push((2, 0), large.as_bytes()).unwrap();
        writer.push((1, 0), b"small").unwrap();
        writer.push((3, 0), large.as_bytes()).unwrap();
        let runs = writer.finish().unwrap();
        assert_eq!(runs.count(), 3);
        assert_eq!(
            merged(&runs),
            [
                ((1, 0), "small".to_string()),
                ((2, 0), large.clone()),
                ((3, 0), large)
            ]
        );
    }

    /// A run the device refuses is removed with its writer, so a spill that
    /// fails leaves no file under the directory. The check runs where the
    /// platform has a device that refuses every write.
    #[cfg(unix)]
    #[test]
    fn a_run_that_is_not_written_whole_is_removed_with_its_writer() {
        let full = Path::new("/dev/full");
        if !full.exists() {
            return;
        }
        let dir = tempfile::tempdir().unwrap();
        std::os::unix::fs::symlink(full, dir.path().join("refused.0.run")).unwrap();
        let mut writer = RunWriter::new(dir.path(), "refused", 64);
        writer.push((1, 0), b"record").unwrap();
        let error = match writer.finish() {
            Ok(_) => panic!("the device took the run"),
            Err(error) => error,
        };
        assert!(
            format!("{error:#}").starts_with("Failed to write "),
            "{error:#}"
        );
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }
}
