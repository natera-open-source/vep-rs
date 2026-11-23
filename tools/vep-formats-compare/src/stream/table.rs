// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the first reading keeps of a side: one row a run of data lines that share a
//! group, and, until the two sides are held against each other, one digest a group.

use super::digest;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

/// A run of data lines of one group, as a block's scan reports it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Run {
    /// The hash of the group's key.
    pub(crate) key: u64,
    pub(crate) lines: u32,
    /// The digest of the run's lines.
    pub(crate) digest: u128,
    /// Whether a line of the run hands its group over.
    pub(crate) forced: bool,
    /// Whether a line of the run hands its group over where the format's lines are read
    /// strictly.
    pub(crate) strict: bool,
}

/// A run of data lines of one group, as the tables keep it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Row {
    pub(crate) group: u32,
    pub(crate) lines: u32,
}

const CHUNK_BITS: usize = 16;
const CHUNK: usize = 1 << CHUNK_BITS;

/// A sequence that grows by whole chunks, so that growing it never copies it.
pub(crate) struct Chunked<T> {
    chunks: Vec<Vec<T>>,
    length: usize,
}

impl<T: Copy> Chunked<T> {
    pub(crate) fn new() -> Chunked<T> {
        Chunked {
            chunks: Vec::new(),
            length: 0,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.length
    }

    pub(crate) fn push(&mut self, value: T) {
        if self.length.is_multiple_of(CHUNK) {
            self.chunks.push(Vec::with_capacity(CHUNK));
        }
        self.chunks
            .last_mut()
            .expect("a chunk was pushed")
            .push(value);
        self.length += 1;
    }

    pub(crate) fn get(&self, index: usize) -> T {
        self.chunks[index >> CHUNK_BITS][index & (CHUNK - 1)]
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = T> + '_ {
        self.chunks.iter().flatten().copied()
    }

    /// The values from `from` to before `to`.
    pub(crate) fn range(&self, from: usize, to: usize) -> impl Iterator<Item = T> + '_ {
        (from..to).map(|index| self.get(index))
    }
}

/// One value a group number, grown to the highest number set; a number never set reads
/// the default.
struct Sparse<T> {
    chunks: Vec<Option<Box<[T]>>>,
}

impl<T: Copy + Default> Sparse<T> {
    fn new() -> Sparse<T> {
        Sparse { chunks: Vec::new() }
    }

    fn get(&self, index: usize) -> T {
        match self.chunks.get(index >> CHUNK_BITS) {
            Some(Some(chunk)) => chunk[index & (CHUNK - 1)],
            _ => T::default(),
        }
    }

    fn at(&mut self, index: usize) -> &mut T {
        let chunk = index >> CHUNK_BITS;
        if self.chunks.len() <= chunk {
            self.chunks.resize_with(chunk + 1, || None);
        }
        &mut self.chunks[chunk].get_or_insert_with(|| vec![T::default(); CHUNK].into_boxed_slice())
            [index & (CHUNK - 1)]
    }
}

/// A set of group numbers.
pub(crate) struct Bits {
    words: Vec<u64>,
}

impl Bits {
    pub(crate) fn new(length: usize) -> Bits {
        Bits {
            words: vec![0; length.div_ceil(64)],
        }
    }

    pub(crate) fn get(&self, index: u32) -> bool {
        self.words[index as usize / 64] >> (index % 64) & 1 == 1
    }

    pub(crate) fn set(&mut self, index: u32) {
        self.words[index as usize / 64] |= 1 << (index % 64);
    }

    pub(crate) fn clear(&mut self, index: u32) {
        self.words[index as usize / 64] &= !(1 << (index % 64));
    }

    pub(crate) fn clear_all(&mut self) {
        self.words.fill(0);
    }

    pub(crate) fn count(&self) -> u64 {
        self.words
            .iter()
            .map(|word| u64::from(word.count_ones()))
            .sum()
    }
}

/// A hasher for keys that are hashes already, written as their eight bytes.
#[derive(Default)]
struct Passed(u64);

impl Hasher for Passed {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, key: &[u8]) {
        if let Ok(key) = <[u8; 8]>::try_from(key) {
            self.0 = u64::from_ne_bytes(key);
        }
    }

    /// The length a slice is hashed with: every key is eight bytes.
    fn write_usize(&mut self, _: usize) {}
}

/// Spreads a key hash over all 64 bits, one-to-one.
fn spread(mut key: u64) -> u64 {
    key ^= key >> 33;
    key = key.wrapping_mul(0xff51_afd7_ed55_8ccd);
    key ^= key >> 33;
    key = key.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    key ^ key >> 33
}

const SHARDS: usize = 256;

/// A key and its number as bytes, which lie closer in a map than the integers.
type Numbers = HashMap<[u8; 8], [u8; 4], BuildHasherDefault<Passed>>;

/// The number of every group of the two sides, by the hash of its key. Two keys of one
/// hash are one group.
pub(crate) struct Groups {
    shards: Vec<Mutex<Numbers>>,
    next: AtomicU64,
}

impl Groups {
    pub(crate) fn new() -> Groups {
        Groups {
            shards: (0..SHARDS).map(|_| Mutex::default()).collect(),
            next: AtomicU64::new(0),
        }
    }

    /// The number of the group of `key`, a new one for a key not seen; `None` once the
    /// numbers below the highest a `u32` holds are used up, so that the count of the
    /// groups is a `u32` too.
    pub(crate) fn number(&self, key: u64) -> Option<u32> {
        let key = spread(key);
        // The map places a key by its low bits and tags it by its high bits; the shard is
        // read off bits between the two.
        let mut shard = self.shards[(key >> 40) as usize % SHARDS]
            .lock()
            .expect("group numbers");
        let bytes = key.to_ne_bytes();
        if let Some(&number) = shard.get(&bytes) {
            return Some(u32::from_ne_bytes(number));
        }
        let number = u32::try_from(self.next.fetch_add(1, Ordering::Relaxed))
            .ok()
            .filter(|&number| number != u32::MAX)?;
        shard.insert(bytes, number.to_ne_bytes());
        Some(number)
    }

    pub(crate) fn count(&self) -> usize {
        self.next.load(Ordering::Relaxed).min(u64::from(u32::MAX)) as usize
    }
}

const PRESENT: u8 = 1;
const FORCED: u8 = 2;
const STRICT: u8 = 4;

/// What the first reading keeps of one side.
pub(crate) struct Table {
    pub(crate) rows: Chunked<Row>,
    /// Per block, the number of its first row.
    pub(crate) block_rows: Vec<u64>,
    /// Per group, the digest of a line of the hash 1 followed by the group's lines, which
    /// tells the group's lines from the same lines after others of the hash 0.
    digests: Sparse<[u64; 2]>,
    flags: Sparse<u8>,
    /// Per group, its lines so far, in a format whose groups have a limit.
    lines: Sparse<u32>,
}

impl Table {
    pub(crate) fn new() -> Table {
        Table {
            rows: Chunked::new(),
            block_rows: Vec::new(),
            digests: Sparse::new(),
            flags: Sparse::new(),
            lines: Sparse::new(),
        }
    }

    /// Starts the rows of the next block.
    pub(crate) fn block(&mut self) {
        self.block_rows.push(self.rows.len() as u64);
    }

    /// Adds a run of the group `group`.
    pub(crate) fn run(&mut self, group: u32, run: &Run) {
        self.rows.push(Row {
            group,
            lines: run.lines,
        });
        let index = group as usize;
        let flags = self.flags.at(index);
        let first = *flags & PRESENT == 0;
        *flags |= PRESENT;
        if run.forced {
            *flags |= FORCED;
        }
        if run.strict {
            *flags |= STRICT;
        }
        let held = self.digests.at(index);
        let before = if first {
            1
        } else {
            u128::from(held[0]) << 64 | u128::from(held[1])
        };
        let joined = digest::join(before, run.digest, u64::from(run.lines));
        *held = [(joined >> 64) as u64, joined as u64];
    }

    /// Adds a run of the group `group`, which the comparison names by the group's key up
    /// to its `limit`th line: the lines past it are lines of the group `placed`, and
    /// both groups are handed over.
    pub(crate) fn limited_run(&mut self, group: u32, run: &Run, limit: u64, placed: u32) {
        let seen = self.lines.at(group as usize);
        let within = limit.saturating_sub(u64::from(*seen));
        *seen = seen.saturating_add(run.lines);
        if u64::from(run.lines) <= within {
            self.run(group, run);
            return;
        }
        // No more than `run.lines`.
        let within = within as u32;
        *self.flags.at(group as usize) |= FORCED;
        if within > 0 {
            *self.flags.at(group as usize) |= PRESENT;
            self.rows.push(Row {
                group,
                lines: within,
            });
        }
        *self.flags.at(placed as usize) |= PRESENT | FORCED;
        self.rows.push(Row {
            group: placed,
            lines: run.lines - within,
        });
    }

    /// The rows of block `block`, as a range of row numbers.
    pub(crate) fn rows_of(&self, block: usize) -> (usize, usize) {
        let from = self.block_rows[block] as usize;
        let to = self
            .block_rows
            .get(block + 1)
            .map_or(self.rows.len(), |&row| row as usize);
        (from, to)
    }

    pub(crate) fn present(&self, group: u32) -> bool {
        self.flags.get(group as usize) & PRESENT != 0
    }

    pub(crate) fn forced(&self, group: u32) -> bool {
        self.flags.get(group as usize) & FORCED != 0
    }

    pub(crate) fn strict(&self, group: u32) -> bool {
        self.flags.get(group as usize) & STRICT != 0
    }

    /// Whether the group's lines have the same digest on both sides.
    pub(crate) fn alike(&self, other: &Table, group: u32) -> bool {
        let index = group as usize;
        self.present(group)
            && other.present(group)
            && self.digests.get(index) == other.digests.get(index)
    }

    /// Drops what only the comparison of the two sides reads.
    pub(crate) fn forget_digests(&mut self) {
        self.digests = Sparse::new();
        self.flags = Sparse::new();
        self.lines = Sparse::new();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_chunked_sequence_reads_back_what_was_pushed() {
        let mut values = Chunked::new();
        for value in 0..(CHUNK as u32 * 2 + 5) {
            values.push(value);
        }
        assert_eq!(values.len(), CHUNK * 2 + 5);
        assert_eq!(values.get(CHUNK + 3), CHUNK as u32 + 3);
        assert!(values.iter().eq(0..(CHUNK as u32 * 2 + 5)));
        assert!(values
            .range(CHUNK - 2, CHUNK + 2)
            .eq(CHUNK as u32 - 2..CHUNK as u32 + 2));
    }

    #[test]
    fn a_set_of_numbers_holds_what_was_set() {
        let mut bits = Bits::new(200);
        for index in [0, 63, 64, 199] {
            bits.set(index);
        }
        assert!(bits.get(63) && bits.get(64) && !bits.get(65));
        assert_eq!(bits.count(), 4);
        bits.clear(63);
        assert!(!bits.get(63));
        bits.clear_all();
        assert_eq!(bits.count(), 0);
    }

    #[test]
    fn a_key_keeps_its_number_and_a_new_key_takes_the_next() {
        let groups = Groups::new();
        let numbers: Vec<u32> = [7, 9, 7, 0, 9, u64::MAX]
            .iter()
            .map(|&key| groups.number(key).unwrap())
            .collect();
        assert_eq!(numbers, [0, 1, 0, 2, 1, 3]);
        assert_eq!(groups.count(), 4);
        // One-to-one: no two of many keys spread to one value.
        let mut spread_keys: Vec<u64> = (0..100_000).map(spread).collect();
        spread_keys.sort_unstable();
        spread_keys.dedup();
        assert_eq!(spread_keys.len(), 100_000);
    }

    #[test]
    fn the_numbers_end_one_below_the_highest_and_the_count_fits_a_u32() {
        let groups = Groups::new();
        groups
            .next
            .store(u64::from(u32::MAX) - 1, Ordering::Relaxed);
        assert_eq!(groups.number(1), Some(u32::MAX - 1));
        assert_eq!(groups.number(2), None);
        assert_eq!(groups.number(3), None);
        assert_eq!(groups.number(1), Some(u32::MAX - 1));
        assert_eq!(groups.count(), u32::MAX as usize);
    }

    #[test]
    fn a_group_s_digest_is_the_same_however_its_lines_are_cut_into_runs() {
        let digest_of = |lines: &[&str]| {
            lines.iter().fold(digest::EMPTY, |digest, line| {
                digest::push(digest, digest::of_line(line.as_bytes()))
            })
        };
        let run = |lines: &[&str], forced: bool| Run {
            key: 1,
            lines: lines.len() as u32,
            digest: digest_of(lines),
            forced,
            strict: false,
        };
        let mut one = Table::new();
        one.block();
        one.run(4, &run(&["a", "b", "c"], false));
        let mut two = Table::new();
        two.block();
        two.run(4, &run(&["a"], false));
        two.run(9, &run(&["x"], true));
        two.block();
        two.run(4, &run(&["b", "c"], false));
        assert!(one.alike(&two, 4));
        assert!(!one.alike(&two, 9) && !two.alike(&one, 9));
        assert!(two.forced(9) && !two.forced(4) && !one.present(9));
        assert_eq!(two.rows_of(0), (0, 2));
        assert_eq!(two.rows_of(1), (2, 3));
        let mut other = Table::new();
        other.block();
        other.run(4, &run(&["a", "c", "b"], false));
        assert!(!one.alike(&other, 4));
        // A line whose hash is zero before the others: another count of lines, another
        // digest.
        let zero = Run {
            key: 1,
            lines: 1,
            digest: digest::push(digest::EMPTY, 0),
            forced: false,
            strict: false,
        };
        let mut longer = Table::new();
        longer.block();
        longer.run(4, &zero);
        longer.run(4, &run(&["a", "b", "c"], false));
        assert!(!one.alike(&longer, 4));
    }
}
