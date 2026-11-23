// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The scratch of one file pair: each side's plugin records in shards by key, so that the
//! two sides of a key lie in one shard and a thread compares one shard at a time.
//!
//! A record is what one row gives the run's plugins: its length, the row's place in its
//! output (the offset of its line and, in a VCF line, the number of its entry), the key's
//! three lengths and its three texts, the count of plugins and, per plugin, its number, the
//! count of its fields and each field's number and value. Every number is a varint.

use crate::order::{put_varint, take_varint};
use std::cmp::Ordering;
use std::fs::{self, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

/// The most shards of a file pair.
pub const MOST_SHARDS: u64 = 4096;
/// What a reading thread's shard buffers hold together before they are written out, and
/// the bounds of one buffer.
const BUFFERED_BYTES: usize = 16 << 20;
const FLUSH_BYTES: std::ops::RangeInclusive<usize> = (4 << 10)..=(256 << 10);

/// One row's key as its three texts.
pub type KeyParts<'a> = (&'a str, &'a str, &'a str);

/// What a row gives one plugin: the plugin's number in the run and its fields' numbers and
/// values in the row's order.
pub type PluginFields = (u16, Vec<(u32, String)>);

/// The directory of one file pair's shards, removed with everything in it when dropped.
pub struct Scratch {
    directory: PathBuf,
    shards: usize,
    /// Per side and shard, held by the thread that appends to the shard's file.
    locks: [Vec<Mutex<()>>; 2],
}

impl Scratch {
    /// Makes `directory` for `shards` shards a side.
    pub fn new(directory: PathBuf, shards: usize) -> io::Result<Scratch> {
        fs::create_dir(&directory)?;
        let locks = || (0..shards).map(|_| Mutex::new(())).collect();
        Ok(Scratch {
            directory,
            shards,
            locks: [locks(), locks()],
        })
    }

    pub fn directory(&self) -> &Path {
        &self.directory
    }

    /// The bytes the shards' files hold.
    pub fn bytes(&self) -> u64 {
        fs::read_dir(&self.directory)
            .map(|entries| {
                entries
                    .filter_map(|entry| entry.ok()?.metadata().ok())
                    .map(|found| found.len())
                    .sum()
            })
            .unwrap_or(0)
    }

    fn path(&self, side: usize, shard: usize) -> PathBuf {
        self.directory.join(format!("side{side}.shard{shard}"))
    }

    /// The shard of a key: the same for both sides and for every thread.
    fn shard_of(&self, key: KeyParts) -> usize {
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        let mut mix =
            |byte: u8| hash = (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3);
        for part in [key.0, key.1, key.2] {
            part.bytes().for_each(&mut mix);
            mix(0xff);
        }
        (hash % self.shards as u64) as usize
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.directory);
    }
}

/// One thread's records of one side on their way to the shards.
pub struct Writer<'a> {
    scratch: &'a Scratch,
    side: usize,
    buffers: Vec<Vec<u8>>,
    flush_bytes: usize,
    record: Vec<u8>,
}

impl<'a> Writer<'a> {
    pub fn new(scratch: &'a Scratch, side: usize) -> Writer<'a> {
        let flush_bytes =
            (BUFFERED_BYTES / scratch.shards).clamp(*FLUSH_BYTES.start(), *FLUSH_BYTES.end());
        Writer {
            scratch,
            side,
            buffers: vec![Vec::new(); scratch.shards],
            flush_bytes,
            record: Vec::new(),
        }
    }

    /// Writes what the row at `place` gives each of `plugins` under `key`.
    pub fn row(
        &mut self,
        place: (u64, u64),
        key: KeyParts,
        plugins: &[PluginFields],
    ) -> io::Result<()> {
        let record = &mut self.record;
        record.clear();
        put_varint(record, place.0);
        put_varint(record, place.1);
        for part in [key.0, key.1, key.2] {
            put_varint(record, part.len() as u64);
        }
        for part in [key.0, key.1, key.2] {
            record.extend_from_slice(part.as_bytes());
        }
        put_varint(record, plugins.len() as u64);
        for (plugin, fields) in plugins {
            put_varint(record, u64::from(*plugin));
            put_varint(record, fields.len() as u64);
            for (number, value) in fields {
                put_varint(record, u64::from(*number));
                put_varint(record, value.len() as u64);
                record.extend_from_slice(value.as_bytes());
            }
        }
        let shard = self.scratch.shard_of(key);
        let buffer = &mut self.buffers[shard];
        put_varint(buffer, self.record.len() as u64);
        buffer.extend_from_slice(&self.record);
        if buffer.len() >= self.flush_bytes {
            self.flush(shard)?;
        }
        Ok(())
    }

    fn flush(&mut self, shard: usize) -> io::Result<()> {
        if self.buffers[shard].is_empty() {
            return Ok(());
        }
        let path = self.scratch.path(self.side, shard);
        let _held = self.scratch.locks[self.side][shard]
            .lock()
            .expect("no writer panics");
        let mut file = OpenOptions::new().create(true).append(true).open(path)?;
        file.write_all(&self.buffers[shard])?;
        self.buffers[shard].clear();
        Ok(())
    }

    /// Writes out what the buffers still hold.
    pub fn finish(mut self) -> io::Result<()> {
        (0..self.buffers.len()).try_for_each(|shard| self.flush(shard))
    }
}

/// Where a key's three texts lie in a shard's bytes.
#[derive(Clone, Copy)]
struct KeyAt {
    at: usize,
    lengths: [u32; 3],
}

impl KeyAt {
    fn parts<'a>(&self, bytes: &'a [u8]) -> (&'a [u8], &'a [u8], &'a [u8]) {
        let [location, allele, feature] = self.lengths.map(|length| length as usize);
        let key = &bytes[self.at..self.at + location + allele + feature];
        (
            &key[..location],
            &key[location..location + allele],
            &key[location + allele..],
        )
    }
}

/// One row's record in a shard's bytes: its place in its output and its key, after which
/// its plugins stand.
struct RowAt {
    place: (u64, u32),
    key: KeyAt,
}

impl RowAt {
    fn plugins_at(&self) -> usize {
        self.key.at
            + self
                .key
                .lengths
                .iter()
                .map(|length| *length as usize)
                .sum::<usize>()
    }
}

/// One side's shard as the comparison reads it: its rows' records in the order of their
/// keys and, under one key, of their places in the output.
pub struct Shard {
    bytes: Vec<u8>,
    rows: Vec<RowAt>,
}

/// What one output holds for one plugin under one key.
#[derive(Clone, Copy)]
pub struct Record<'a> {
    pub key: KeyParts<'a>,
    pub plugin: u16,
    fields: &'a [u8],
}

fn text(bytes: &[u8]) -> &str {
    // A writer takes text alone, so the bytes read back are text.
    std::str::from_utf8(bytes).expect("a record holds what a writer wrote")
}

/// Steps over one plugin's fields, which start at `at` with their count.
fn skip_fields(bytes: &[u8], at: &mut usize) {
    for _ in 0..take_varint(bytes, at) {
        take_varint(bytes, at);
        *at += take_varint(bytes, at) as usize;
    }
}

impl<'a> Record<'a> {
    /// The record's fields as their numbers and values, in the row's order.
    pub fn fields(&self) -> impl Iterator<Item = (u32, &'a str)> + 'a {
        let bytes = self.fields;
        let mut at = 0;
        let count = take_varint(bytes, &mut at);
        (0..count).map(move |_| {
            let number = take_varint(bytes, &mut at) as u32;
            let length = take_varint(bytes, &mut at) as usize;
            at += length;
            (number, text(&bytes[at - length..at]))
        })
    }

    /// The key and then the plugin: the order a shard's records stand in. Text orders by
    /// its code points, which is the order of its bytes.
    pub fn cmp_place(&self, other: &Record) -> Ordering {
        (self.key, self.plugin).cmp(&(other.key, other.plugin))
    }
}

impl Shard {
    /// Reads a side's shard of the scratch; a shard no record was written to is empty.
    pub fn load(scratch: &Scratch, side: usize, shard: usize) -> io::Result<Shard> {
        let bytes = match fs::read(scratch.path(side, shard)) {
            Ok(bytes) => bytes,
            Err(error) if error.kind() == io::ErrorKind::NotFound => Vec::new(),
            Err(error) => return Err(error),
        };
        let mut rows = Vec::new();
        let mut at = 0;
        while at < bytes.len() {
            let length = take_varint(&bytes, &mut at) as usize;
            let end = at + length;
            let place = (
                take_varint(&bytes, &mut at),
                take_varint(&bytes, &mut at) as u32,
            );
            let lengths = [(); 3].map(|()| take_varint(&bytes, &mut at) as u32);
            rows.push(RowAt {
                place,
                key: KeyAt { at, lengths },
            });
            at = end;
        }
        rows.sort_unstable_by(|a, b| {
            (a.key.parts(&bytes), a.place).cmp(&(b.key.parts(&bytes), b.place))
        });
        Ok(Shard { bytes, rows })
    }

    /// What the shard holds per key and plugin: the fields of the last row that gave the
    /// plugin a field, in the order of the keys and, under one key, of the plugins.
    pub fn records(&self) -> Records<'_> {
        Records {
            shard: self,
            next_row: 0,
            key: None,
            of_the_key: Vec::new(),
            given: 0,
        }
    }
}

/// A walk over a shard's records, one key's rows in hand at a time.
pub struct Records<'a> {
    shard: &'a Shard,
    next_row: usize,
    key: Option<KeyParts<'a>>,
    /// Per plugin of the key in hand, where the fields of its last row stand.
    of_the_key: Vec<(u16, usize)>,
    given: usize,
}

impl<'a> Iterator for Records<'a> {
    type Item = Record<'a>;

    fn next(&mut self) -> Option<Record<'a>> {
        let (bytes, rows) = (&self.shard.bytes, &self.shard.rows);
        if self.given == self.of_the_key.len() {
            let first = rows.get(self.next_row)?;
            let key = first.key.parts(bytes);
            self.of_the_key.clear();
            self.given = 0;
            // The rows of the key in their output's order: a later row stands for an
            // earlier one, plugin by plugin.
            while let Some(row) = rows
                .get(self.next_row)
                .filter(|row| row.key.parts(bytes) == key)
            {
                let mut at = row.plugins_at();
                for _ in 0..take_varint(bytes, &mut at) {
                    let plugin = take_varint(bytes, &mut at) as u16;
                    match self
                        .of_the_key
                        .iter_mut()
                        .find(|(known, _)| *known == plugin)
                    {
                        Some(known) => known.1 = at,
                        None => self.of_the_key.push((plugin, at)),
                    }
                    skip_fields(bytes, &mut at);
                }
                self.next_row += 1;
            }
            self.of_the_key.sort_unstable();
            self.key = Some((text(key.0), text(key.1), text(key.2)));
        }
        let (plugin, fields_at) = *self.of_the_key.get(self.given)?;
        self.given += 1;
        Some(Record {
            key: self.key?,
            plugin,
            fields: &bytes[fields_at..],
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str, shards: usize) -> Scratch {
        let directory = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-shards-{}-{name}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&directory);
        Scratch::new(directory, shards).unwrap()
    }

    fn fields(values: &[(u32, &str)]) -> Vec<(u32, String)> {
        values
            .iter()
            .map(|(number, value)| (*number, value.to_string()))
            .collect()
    }

    type Read = Vec<((String, String, String), u16, Vec<(u32, String)>)>;

    fn read(scratch: &Scratch, side: usize) -> Read {
        let mut found = Vec::new();
        for shard in 0..scratch.shards {
            let shard = Shard::load(scratch, side, shard).unwrap();
            for record in shard.records() {
                found.push((
                    (
                        record.key.0.to_string(),
                        record.key.1.to_string(),
                        record.key.2.to_string(),
                    ),
                    record.plugin,
                    record
                        .fields()
                        .map(|(number, value)| (number, value.to_string()))
                        .collect(),
                ));
            }
        }
        found.sort();
        found
    }

    #[test]
    fn the_last_row_of_a_key_stands_for_a_plugin_whatever_order_the_rows_are_written_in() {
        for shards in [1, 2, 7, 64] {
            let scratch = scratch(&format!("last-{shards}"), shards);
            let mut early = Writer::new(&scratch, 0);
            let mut late = Writer::new(&scratch, 0);
            let first = ("1:10", "T", "ENST1");
            // The later rows are written first, as a thread that reads a later range may.
            late.row((900, 0), first, &[(0, fields(&[(4, "late")]))])
                .unwrap();
            late.row((900, 1), first, &[(0, fields(&[(4, "latest")]))])
                .unwrap();
            late.row(
                (950, 0),
                ("1:10", "T", "ENST2"),
                &[(1, fields(&[(7, "only")]))],
            )
            .unwrap();
            early
                .row(
                    (100, 0),
                    first,
                    &[
                        (1, fields(&[(7, "other plugin")])),
                        (0, fields(&[(4, "early"), (5, "x")])),
                    ],
                )
                .unwrap();
            early
                .row((200, 0), ("1:100", "", "ENST1"), &[(0, fields(&[(4, "")]))])
                .unwrap();
            late.finish().unwrap();
            early.finish().unwrap();
            let key = |location: &str, allele: &str, feature: &str| {
                (
                    location.to_string(),
                    allele.to_string(),
                    feature.to_string(),
                )
            };
            assert_eq!(
                read(&scratch, 0),
                vec![
                    (key("1:10", "T", "ENST1"), 0, fields(&[(4, "latest")])),
                    (key("1:10", "T", "ENST1"), 1, fields(&[(7, "other plugin")])),
                    (key("1:10", "T", "ENST2"), 1, fields(&[(7, "only")])),
                    (key("1:100", "", "ENST1"), 0, fields(&[(4, "")])),
                ],
                "{shards} shards"
            );
            assert!(read(&scratch, 1).is_empty());
        }
    }

    #[test]
    fn a_shards_records_stand_in_the_order_of_the_key_as_three_texts_and_then_of_the_plugin() {
        let scratch = scratch("order", 1);
        let mut writer = Writer::new(&scratch, 1);
        let keys = [
            ("1:100", "T", "E"),
            ("1:10", "T", "E"),
            ("1:10", "", "E"),
            ("1:10", "T", ""),
            ("1:0", "0T", "E"),
            ("\u{e9}:1", "T", "E"),
            ("1:10", "T", "ENST2"),
            ("1:10", "TT", "E"),
        ];
        for (at, key) in keys.iter().enumerate() {
            let plugins = [(9, fields(&[(1, "v")])), (3, fields(&[(2, "w")]))];
            writer.row((at as u64, 0), *key, &plugins).unwrap();
        }
        writer.finish().unwrap();
        let shard = Shard::load(&scratch, 1, 0).unwrap();
        let read: Vec<(KeyParts, u16)> = shard
            .records()
            .map(|record| (record.key, record.plugin))
            .collect();
        let mut sorted = keys.to_vec();
        sorted.sort();
        let expected: Vec<(KeyParts, u16)> = sorted
            .into_iter()
            .flat_map(|key| [(key, 3), (key, 9)])
            .collect();
        assert_eq!(read, expected);
        assert_eq!(read[0].0, ("1:0", "0T", "E"));
    }

    #[test]
    fn the_scratch_is_removed_with_its_shards() {
        let scratch = scratch("removed", 3);
        let directory = scratch.directory().to_path_buf();
        let mut writer = Writer::new(&scratch, 0);
        writer
            .row((0, 0), ("1:5", "T", "E"), &[(0, fields(&[(1, "v")]))])
            .unwrap();
        writer.finish().unwrap();
        assert!(scratch.bytes() > 10);
        assert!(fs::read_dir(&directory).unwrap().count() > 0);
        drop(scratch);
        assert!(!directory.exists());
    }
}
