// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The second read of the table: every row written with its class behind it, the counts
//! of the summary taken on the way, and the summary's text.
//!
//! The calling thread reads and lays out the rows; where the bound allows, one thread
//! writes them and one counts, each taking its batches in row order.

use crate::classes::Class;
use crate::csv::{push_field, Reader, Record};
use crate::intern::{Batch, FoldBuild, Interner, ABSENT};
use crate::pyfmt::{percent, push_padded};
use crate::scan;
use crate::table::{start_of, Header, Shown, CLASS_COLUMN};
use crate::text::Text;
use std::collections::HashMap;
use std::hash::Hash;
use std::io::{self, Read, Write};
use std::sync::mpsc::{sync_channel, SyncSender};
use std::thread::ScopedJoinHandle;

/// The bytes of the written table handed on at a time, and the room a chunk has for the
/// row that carries it past them.
const CHUNK_BYTES: usize = 1 << 20;
const CHUNK_ROOM: usize = CHUNK_BYTES + (1 << 16);

/// The rows counted at a time.
const BATCH_ROWS: usize = 1 << 15;

/// The chunks or batches that may wait for the writing or the counting thread.
const WAITING: usize = 4;

/// The rows of the `By file + mismatch type` block.
const TOP_FILES: usize = 30;

/// What a second read finds changed under it.
const CHANGED: &str = "the table changed while it was read";

/// Why the second read ends early.
#[derive(Debug)]
pub enum Failure {
    Input(scan::Failure),
    Write(io::Error),
}

impl From<scan::Failure> for Failure {
    fn from(failure: scan::Failure) -> Failure {
        Failure::Input(failure)
    }
}

impl From<&'static str> for Failure {
    fn from(detail: &'static str) -> Failure {
        Failure::Input(detail.into())
    }
}

/// Counts under keys, the keys kept in the order they are first met.
struct Counter<K> {
    places: HashMap<K, usize, FoldBuild>,
    counts: Vec<(K, u64)>,
}

impl<K: Hash + Eq + Copy> Counter<K> {
    fn new() -> Counter<K> {
        Counter {
            places: HashMap::default(),
            counts: Vec::new(),
        }
    }

    fn add(&mut self, key: K) {
        let place = *self.places.entry(key).or_insert(self.counts.len());
        if place == self.counts.len() {
            self.counts.push((key, 0));
        }
        self.counts[place].1 += 1;
    }

    /// The counts from the largest down, equal counts in the order their keys were met.
    fn descending(&self) -> Vec<(K, u64)> {
        let mut counts = self.counts.clone();
        counts.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
        counts
    }
}

/// The counts of the summary. A class, a source and a file are counted as the text they
/// are, a value that its row ends before as no text.
pub struct Tally {
    rows: u64,
    types: Interner,
    sources: Interner,
    files: Interner,
    named: [u32; Class::NAMED.len()],
    by_type: Counter<u32>,
    by_type_and_source: Counter<[u32; 2]>,
    by_file_and_type: Counter<[u32; 2]>,
}

impl Tally {
    fn new() -> Tally {
        let mut types = Interner::new();
        let named =
            Class::NAMED.map(|class| types.intern(class.name().unwrap_or_default().as_bytes()));
        Tally {
            rows: 0,
            types,
            sources: Interner::new(),
            files: Interner::new(),
            named,
            by_type: Counter::new(),
            by_type_and_source: Counter::new(),
            by_file_and_type: Counter::new(),
        }
    }

    /// Counts the rows of one batch: their classes, and their source, file and read class
    /// in the three columns.
    fn absorb(&mut self, classes: &[Class], [sources, files, read]: &[Batch; 3]) {
        let values = sources.values().zip(files.values()).zip(read.values());
        for (class, ((source, file), read)) in classes.iter().zip(values) {
            let kind = match Class::NAMED.iter().position(|named| named == class) {
                Some(index) => self.named[index],
                None => self.types.intern_optional(read),
            };
            let source = self.sources.intern_optional(source);
            let file = self.files.intern_optional(file);
            self.rows += 1;
            self.by_type.add(kind);
            self.by_type_and_source.add([kind, source]);
            self.by_file_and_type.add([file, kind]);
        }
    }
}

fn text_of(values: &Interner, number: u32) -> Option<&[u8]> {
    (number != ABSENT).then(|| values.value(number))
}

/// Where the laid-out rows go: to the sink on the calling thread, or to the thread that
/// holds it.
enum Writing<'scope, W> {
    Here(W),
    Thread(
        SyncSender<Vec<u8>>,
        ScopedJoinHandle<'scope, io::Result<()>>,
    ),
}

impl<W: Write> Writing<'_, W> {
    /// Ends the writing and returns what became of it.
    fn finish(self) -> io::Result<()> {
        match self {
            Writing::Here(mut sink) => sink.flush(),
            Writing::Thread(sender, thread) => {
                drop(sender);
                thread.join().expect("the writing thread")
            }
        }
    }
}

/// Where the batches to count go: to the counts on the calling thread, or to the thread
/// that holds them.
enum Counting<'scope, 'a> {
    Here(Box<Tally>),
    Thread(
        SyncSender<(&'a [Class], [Batch; 3])>,
        ScopedJoinHandle<'scope, Tally>,
    ),
}

impl Counting<'_, '_> {
    /// Ends the counting and returns the counts.
    fn finish(self) -> Tally {
        match self {
            Counting::Here(tally) => *tally,
            Counting::Thread(sender, thread) => {
                drop(sender);
                thread.join().expect("the counting thread")
            }
        }
    }
}

struct Rows<'scope, 'a, W> {
    header: &'a Header,
    classes: &'a [Class],
    header_written: bool,
    /// The rows laid out, and how many of them are handed on to be counted.
    row: usize,
    counted: usize,
    out: Vec<u8>,
    batches: [Batch; 3],
    writing: Writing<'scope, W>,
    counting: Counting<'scope, 'a>,
}

impl<W: Write> Rows<'_, '_, W> {
    fn record(&mut self, record: Record<'_>) -> Result<(), Failure> {
        if !self.header_written {
            self.header_written = true;
            for name in &self.header.names {
                push_field(&mut self.out, name);
                self.out.push(b'\t');
            }
            self.out.extend_from_slice(CLASS_COLUMN);
            self.out.extend_from_slice(b"\r\n");
            return Ok(());
        }
        if record.is_empty() {
            return Ok(());
        }
        let (Some(class), Some(location)) =
            (self.classes.get(self.row), record.get(self.header.location))
        else {
            return Err(CHANGED.into());
        };
        let read = self.header.class.and_then(|column| record.get(column));
        let class_text = match class.name() {
            Some(name) => name.as_bytes(),
            None => read.unwrap_or_default(),
        };
        for shown in &self.header.shown {
            let value = match shown {
                Shown::Field(column) => record.get(*column).unwrap_or_default(),
                Shown::Start => start_of(location),
                Shown::Class => class_text,
            };
            push_field(&mut self.out, value);
            self.out.push(b'\t');
        }
        push_field(&mut self.out, class_text);
        self.out.extend_from_slice(b"\r\n");
        self.batches[0].push(record.get(self.header.source));
        self.batches[1].push(record.get(self.header.file));
        self.batches[2].push(read);
        self.row += 1;
        if self.row - self.counted == BATCH_ROWS {
            self.count();
        }
        if self.out.len() >= CHUNK_BYTES {
            self.write()?;
        }
        Ok(())
    }

    /// Hands on the rows laid out since the last call to be counted.
    fn count(&mut self) {
        let classes = &self.classes[self.counted..self.row];
        let batches = std::mem::take(&mut self.batches);
        self.counted = self.row;
        match &mut self.counting {
            Counting::Here(tally) => tally.absorb(classes, &batches),
            // The counting thread ends only when its sender is dropped.
            Counting::Thread(sender, _) => drop(sender.send((classes, batches))),
        }
    }

    /// Hands on the bytes laid out since the last call to be written.
    fn write(&mut self) -> Result<(), Failure> {
        match &mut self.writing {
            Writing::Here(sink) => {
                sink.write_all(&self.out).map_err(Failure::Write)?;
                self.out.clear();
                Ok(())
            }
            Writing::Thread(sender, _) => {
                let chunk = std::mem::replace(&mut self.out, Vec::with_capacity(CHUNK_ROOM));
                // The writing thread ends early only on an error of its own, which is the
                // one reported.
                sender.send(chunk).map_err(|_| CHANGED.into())
            }
        }
    }
}

/// Reads the table a second time and writes it to `sink` with the class of every row in a
/// last column, `mismatch_type`: rows end in CR LF, a row shorter than the header is
/// filled with empty fields and the fields behind the header's last column are dropped.
/// Returns the counts of the summary. `threads` bounds the threads at work.
pub fn write_table<R: Read, W: Write + Send>(
    mut text: Text<R>,
    header: &Header,
    classes: &[Class],
    sink: W,
    threads: usize,
) -> Result<Tally, Failure> {
    std::thread::scope(|scope| {
        let writing = if threads >= 2 {
            let (sender, receiver) = sync_channel::<Vec<u8>>(WAITING);
            let thread = scope.spawn(move || {
                let mut sink = sink;
                for chunk in receiver {
                    sink.write_all(&chunk)?;
                }
                sink.flush()
            });
            Writing::Thread(sender, thread)
        } else {
            Writing::Here(sink)
        };
        let counting = if threads >= 3 {
            let (sender, receiver) = sync_channel::<(&[Class], [Batch; 3])>(WAITING);
            let thread = scope.spawn(move || {
                let mut tally = Tally::new();
                for (classes, batches) in receiver {
                    tally.absorb(classes, &batches);
                }
                tally
            });
            Counting::Thread(sender, thread)
        } else {
            Counting::Here(Box::new(Tally::new()))
        };
        let mut rows = Rows {
            header,
            classes,
            header_written: false,
            row: 0,
            counted: 0,
            out: Vec::with_capacity(CHUNK_ROOM),
            batches: Default::default(),
            writing,
            counting,
        };
        let mut read = read_records(&mut text, &mut rows);
        if read.is_ok() && (!rows.header_written || rows.row != classes.len()) {
            read = Err(CHANGED.into());
        }
        if read.is_ok() {
            rows.count();
            read = rows.write();
        }
        let tally = rows.counting.finish();
        rows.writing.finish().map_err(Failure::Write)?;
        read?;
        Ok(tally)
    })
}

fn read_records<R: Read, W: Write>(
    text: &mut Text<R>,
    rows: &mut Rows<'_, '_, W>,
) -> Result<(), Failure> {
    let mut reader = Reader::new();
    let mut run = Vec::new();
    while text.next_run(&mut run).map_err(scan::Failure::from)? {
        reader.feed(&run, &mut |record| rows.record(record))?;
    }
    reader.finish(&mut |record| rows.record(record))
}

/// The summary: the total, the rows of each class, of each class and source, and of the
/// thirty largest pairs of file and class. A count under a class, a source or a file that
/// its row ends before has no text to print: the summary ends at the line before it, and
/// the second value then names the column.
pub fn summary(tally: &Tally, header: &Header) -> (Vec<u8>, Option<&'static str>) {
    let mut out = Vec::new();
    let rule = "=".repeat(70);
    out.extend_from_slice(
        format!(
            "{rule}\nDISCORDANT CLASSIFICATION SUMMARY\n{rule}\n\nTotal discordants: {}\n",
            tally.rows
        )
        .as_bytes(),
    );
    let class_column = std::str::from_utf8(CLASS_COLUMN).unwrap_or_default();
    let share = |count: u64| format!(" {count:>8} ({}%)\n", percent(count, tally.rows));

    out.extend_from_slice(b"\nBy mismatch type:\n");
    for (kind, count) in tally.by_type.descending() {
        let Some(kind) = text_of(&tally.types, kind) else {
            return (out, Some(class_column));
        };
        out.extend_from_slice(b"  ");
        push_padded(&mut out, kind, 25);
        out.extend_from_slice(share(count).as_bytes());
    }

    out.extend_from_slice(b"\nBy mismatch type + direction:\n");
    for ([kind, source], count) in tally.by_type_and_source.descending() {
        let (Some(kind), Some(source)) =
            (text_of(&tally.types, kind), text_of(&tally.sources, source))
        else {
            return (out, Some(header.source_name));
        };
        out.extend_from_slice(b"  ");
        push_padded(&mut out, kind, 25);
        out.push(b' ');
        push_padded(&mut out, source, 12);
        out.extend_from_slice(format!(" {count:>8}\n").as_bytes());
    }

    out.extend_from_slice(b"\nBy file + mismatch type (top 30):\n");
    for ([file, kind], count) in tally
        .by_file_and_type
        .descending()
        .into_iter()
        .take(TOP_FILES)
    {
        let (Some(file), Some(kind)) = (text_of(&tally.files, file), text_of(&tally.types, kind))
        else {
            return (out, Some(header.file_name));
        };
        out.extend_from_slice(b"  ");
        push_padded(&mut out, file, 35);
        out.push(b' ');
        push_padded(&mut out, kind, 25);
        out.extend_from_slice(share(count).as_bytes());
    }
    (out, None)
}
