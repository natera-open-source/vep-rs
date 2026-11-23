// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Values as numbers: every distinct value of a column is held once and named by a
//! number, so that the rows of a table are compared and paired as tuples of numbers and
//! the memory of a run follows its rows and its distinct values, not its text.

use std::hash::{BuildHasherDefault, Hasher};

/// The number that stands for the value of a column that a row ends before.
pub const ABSENT: u32 = u32::MAX;

/// The most rows a table may hold, so that a row index and the number of a value both fit
/// 32 bits beside `ABSENT`.
pub const MOST_ROWS: usize = (u32::MAX - 1) as usize;

const SEED: u64 = 0x243f_6a88_85a3_08d3;
const MULTIPLIER: u64 = 0x9e37_79b9_7f4a_7c15;

/// The two halves of a 128-bit product, folded together.
fn fold(a: u64, b: u64) -> u64 {
    let product = u128::from(a) * u128::from(b);
    (product as u64) ^ ((product >> 64) as u64)
}

/// A multiply-and-fold hash over bytes and small integers. It orders nothing and no byte
/// written depends on it.
#[derive(Clone, Copy)]
pub struct FoldHasher(u64);

impl Default for FoldHasher {
    fn default() -> FoldHasher {
        FoldHasher(SEED)
    }
}

impl Hasher for FoldHasher {
    fn write(&mut self, bytes: &[u8]) {
        let mut words = bytes.chunks_exact(8);
        for word in &mut words {
            let word = u64::from_le_bytes(word.try_into().expect("eight bytes"));
            self.0 = fold(self.0 ^ word, MULTIPLIER);
        }
        let mut last = [0u8; 8];
        last[..words.remainder().len()].copy_from_slice(words.remainder());
        self.0 = fold(self.0 ^ u64::from_le_bytes(last), MULTIPLIER);
        self.0 = fold(self.0 ^ bytes.len() as u64, SEED);
    }

    fn write_u32(&mut self, number: u32) {
        self.0 = fold(self.0 ^ u64::from(number), MULTIPLIER);
    }

    fn write_usize(&mut self, number: usize) {
        self.0 = fold(self.0 ^ number as u64, MULTIPLIER);
    }

    fn finish(&self) -> u64 {
        self.0
    }
}

pub type FoldBuild = BuildHasherDefault<FoldHasher>;

fn hash_of(value: &[u8]) -> u64 {
    let mut hasher = FoldHasher::default();
    hasher.write(value);
    hasher.finish()
}

const EMPTY: u64 = u64::MAX;

/// The most slots a table holds: a slot is found by the high half of its value's hash.
const MOST_SLOTS: usize = 1 << 32;

/// The distinct values of one column, numbered from 0 in the order they are first met.
pub struct Interner {
    /// The values end to end, and where each ends.
    bytes: Vec<u8>,
    ends: Vec<usize>,
    /// An open-addressed table: the low half of a slot is a value's number, the high half
    /// the high half of its hash, whose low bits are the slot the value is looked for
    /// from. A table grows by its slots alone, no value hashed again.
    slots: Vec<u64>,
    /// The number of the value asked for last: the rows of one variant follow each other,
    /// so a value is most often the one before it.
    last: u32,
}

impl Default for Interner {
    fn default() -> Interner {
        Interner::new()
    }
}

impl Interner {
    pub fn new() -> Interner {
        Interner {
            bytes: Vec::new(),
            ends: Vec::new(),
            slots: vec![EMPTY; 1 << 10],
            last: ABSENT,
        }
    }

    /// The value a number stands for.
    pub fn value(&self, number: u32) -> &[u8] {
        let index = number as usize;
        let from = if index == 0 { 0 } else { self.ends[index - 1] };
        &self.bytes[from..self.ends[index]]
    }

    /// The number of `value`, which is held from this call on. The caller keeps the
    /// distinct values of a column at or under `MOST_ROWS`.
    pub fn intern(&mut self, value: &[u8]) -> u32 {
        if self.last != ABSENT && self.value(self.last) == value {
            return self.last;
        }
        self.last = self.look_up(value);
        self.last
    }

    fn look_up(&mut self, value: &[u8]) -> u32 {
        let tag = hash_of(value) >> 32;
        let mask = self.slots.len() - 1;
        let mut at = tag as usize & mask;
        loop {
            let slot = self.slots[at];
            if slot == EMPTY {
                break;
            }
            if slot >> 32 == tag && self.value(slot as u32) == value {
                return slot as u32;
            }
            at = (at + 1) & mask;
        }
        let number = self.ends.len() as u32;
        self.bytes.extend_from_slice(value);
        self.ends.push(self.bytes.len());
        self.slots[at] = (tag << 32) | u64::from(number);
        // Three slots in four at most are taken while the table can grow.
        if self.ends.len() * 4 > self.slots.len() * 3 && self.slots.len() < MOST_SLOTS {
            self.grow();
        }
        number
    }

    /// The number of an optional value: `ABSENT` for none.
    pub fn intern_optional(&mut self, value: Option<&[u8]>) -> u32 {
        value.map_or(ABSENT, |value| self.intern(value))
    }

    fn grow(&mut self) {
        let mask = self.slots.len() * 2 - 1;
        let mut slots = vec![EMPTY; mask + 1];
        for slot in self.slots.iter().filter(|slot| **slot != EMPTY) {
            let mut at = (slot >> 32) as usize & mask;
            while slots[at] != EMPTY {
                at = (at + 1) & mask;
            }
            slots[at] = *slot;
        }
        self.slots = slots;
    }
}

/// Values of one column on their way to its interner: the bytes end to end and the length
/// of each, `ABSENT` standing for a value that its row ends before.
#[derive(Default)]
pub struct Batch {
    bytes: Vec<u8>,
    lengths: Vec<u32>,
}

impl Batch {
    pub fn push(&mut self, value: Option<&[u8]>) {
        match value {
            Some(value) => {
                self.bytes.extend_from_slice(value);
                self.lengths.push(value.len() as u32);
            }
            None => self.lengths.push(ABSENT),
        }
    }

    /// The values in the order they were pushed.
    pub fn values(&self) -> impl Iterator<Item = Option<&[u8]>> + '_ {
        let mut from = 0;
        self.lengths.iter().map(move |length| {
            (*length != ABSENT).then(|| {
                let value = &self.bytes[from..from + *length as usize];
                from += *length as usize;
                value
            })
        })
    }
}

/// A set of row indices.
#[derive(Default)]
pub struct Bits {
    words: Vec<u64>,
    rows: usize,
}

impl Bits {
    /// The empty set over `rows` rows.
    pub fn new(rows: usize) -> Bits {
        Bits {
            words: vec![0; rows.div_ceil(64)],
            rows,
        }
    }

    /// Adds one row behind the last, in the set or outside it.
    pub fn push(&mut self, held: bool) {
        if self.rows.is_multiple_of(64) {
            self.words.push(0);
        }
        self.rows += 1;
        if held {
            self.set(self.rows - 1);
        }
    }

    pub fn get(&self, row: usize) -> bool {
        self.words[row / 64] & (1 << (row % 64)) != 0
    }

    /// The number of rows in the set.
    pub fn count(&self) -> usize {
        self.words
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum()
    }

    pub fn set(&mut self, row: usize) {
        self.words[row / 64] |= 1 << (row % 64);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_values_share_a_number_and_numbers_follow_first_sight() {
        let mut interner = Interner::new();
        let values: Vec<String> = (0..5000).map(|n| format!("value-{}", n % 1777)).collect();
        let numbers: Vec<u32> = values
            .iter()
            .map(|value| interner.intern(value.as_bytes()))
            .collect();
        assert_eq!(interner.ends.len(), 1777);
        for (index, value) in values.iter().enumerate() {
            assert_eq!(numbers[index] as usize, index % 1777);
            assert_eq!(interner.value(numbers[index]), value.as_bytes());
        }
    }

    #[test]
    fn a_table_that_grows_keeps_every_number() {
        let mut interner = Interner::new();
        let value = |n: u32| format!("{}:{n}-{}", n % 23, n.wrapping_mul(2_654_435_761));
        for n in 0..200_000 {
            assert_eq!(interner.intern(value(n).as_bytes()), n);
        }
        assert!(interner.slots.len() > 200_000 && interner.slots.len() <= 1 << 19);
        for n in (0..200_000).rev() {
            assert_eq!(interner.intern(value(n).as_bytes()), n);
            assert_eq!(interner.value(n), value(n).as_bytes());
        }
        assert_eq!(interner.ends.len(), 200_000);
    }

    #[test]
    fn the_empty_value_is_a_value_and_an_absent_one_is_not() {
        let mut interner = Interner::new();
        assert_eq!(interner.intern_optional(Some(b"")), 0);
        assert_eq!(interner.intern_optional(None), ABSENT);
        assert_eq!(interner.intern_optional(Some(b"\0")), 1);
        assert_eq!(interner.intern_optional(Some(b"\0\0")), 2);
        assert_eq!(interner.intern_optional(Some(b"")), 0);
        assert_eq!(interner.value(0), b"");
        assert_eq!(interner.ends.len(), 3);
    }

    #[test]
    fn a_batch_hands_back_its_values_in_order() {
        let mut batch = Batch::default();
        assert_eq!(batch.values().count(), 0);
        let values = [Some(&b"a"[..]), None, Some(b""), Some(b"bcd"), None];
        for value in values {
            batch.push(value);
        }
        assert_eq!(batch.values().collect::<Vec<_>>(), values);
    }

    #[test]
    fn bits_hold_the_rows_set_and_the_rows_pushed() {
        let mut bits = Bits::new(130);
        for row in [0, 63, 64, 129] {
            bits.set(row);
        }
        let held: Vec<usize> = (0..130).filter(|row| bits.get(*row)).collect();
        assert_eq!(held, [0, 63, 64, 129]);
        assert_eq!(bits.count(), 4);
        let mut pushed = Bits::default();
        for row in 0..200 {
            pushed.push(row % 7 == 0);
        }
        pushed.set(1);
        let held: Vec<usize> = (0..200).filter(|row| pushed.get(*row)).collect();
        let expected: Vec<usize> = (0..200).filter(|row| row % 7 == 0 || *row == 1).collect();
        assert_eq!(held, expected);
        assert_eq!(pushed.count(), expected.len());
    }
}
