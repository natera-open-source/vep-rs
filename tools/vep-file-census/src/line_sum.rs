// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The lines of a file in any order: how many there are and the sum of their SHA256
//! digests modulo 2**256.
//!
//! A line ends at `\n`, at `\r\n` and at a lone `\r`, as universal newlines end one, and
//! a last line needs no line end. An empty line is no line. The bytes of a line without
//! its line end are digested, each digest is a 256-bit integer with its first byte the
//! most significant, and the sum wraps at 2**256. A sum takes its terms in any order, so
//! the runs of a file are summed on the worker threads and added up as they arrive, and
//! two files that hold the same lines, each as often, come to the same count and sum
//! wherever the lines stand.

use memchr::memchr2;
use sha2::{Digest, Sha256};

/// The count of some lines and the sum of their digests: four 64-bit limbs, the most
/// significant first.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub struct LineSum {
    pub lines: u64,
    sum: [u64; 4],
}

impl LineSum {
    /// The lines of `data`, a run of whole lines, but those that start with `without`.
    pub fn of(data: &[u8], without: Option<&[u8]>) -> LineSum {
        let mut made = LineSum::default();
        let mut rest = data;
        while !rest.is_empty() {
            let (line, after) = match memchr2(b'\n', b'\r', rest) {
                Some(end) => (&rest[..end], &rest[end + 1..]),
                None => (rest, &rest[rest.len()..]),
            };
            rest = after;
            if line.is_empty() || without.is_some_and(|text| line.starts_with(text)) {
                continue;
            }
            let digest: [u8; 32] = Sha256::digest(line).into();
            let mut limbs = [0u64; 4];
            for (limb, bytes) in limbs.iter_mut().zip(digest.chunks_exact(8)) {
                *limb = u64::from_be_bytes(bytes.try_into().expect("a chunk of eight bytes"));
            }
            made.add(limbs);
            made.lines += 1;
        }
        made
    }

    /// Adds `limbs` to the sum; the carry out of the most significant limb is dropped.
    fn add(&mut self, limbs: [u64; 4]) {
        let mut carry = false;
        for (limb, more) in self.sum.iter_mut().zip(limbs).rev() {
            let (sum, over) = limb.overflowing_add(more);
            let (sum, over_carry) = sum.overflowing_add(u64::from(carry));
            *limb = sum;
            carry = over || over_carry;
        }
    }

    /// Adds the lines of `other`.
    pub fn merge(&mut self, other: LineSum) {
        self.lines += other.lines;
        self.add(other.sum);
    }

    /// The sum in 64 hexadecimal digits.
    pub fn hex(&self) -> String {
        self.sum.iter().map(|limb| format!("{limb:016x}")).collect()
    }
}

#[cfg(test)]
#[path = "line_sum_tests.rs"]
mod tests;
