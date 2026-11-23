// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a format supplies to the stream: which lines are header lines, how a data line
//! names its group, what hands a group over, and the shapes under which nothing is kept
//! back.

use std::collections::BTreeMap;
use std::ops::Range;

/// What a side's scan reads off one data line.
pub(crate) struct Scanned {
    /// The hash of the key of the line's group. Two keys of one hash are one group.
    pub(crate) key: u64,
    /// The hash of what the line holds: two lines of one group that the comparison reads
    /// alike have one hash, and the digest of a group is made of its lines' hashes.
    pub(crate) hash: u128,
    /// Whether the line hands its group over.
    pub(crate) forced: bool,
    /// Whether the comparison may end with an error on the line.
    pub(crate) fragile: bool,
    /// The shapes the line shows under which every line is handed over, one bit a shape,
    /// numbered as `Grouping::SHAPES`.
    pub(crate) shapes: u32,
    /// Whether the line hands its group over where `Grouping::strict` holds.
    pub(crate) strict: bool,
    /// What the line carries, as the bits `Grouping::strict` reads.
    pub(crate) carried: u64,
}

/// The counts of a side's kept-back lines, each by name; a name of no count is left out.
#[derive(Default)]
pub(crate) struct Tally {
    /// Per field name, the entries that carry a value.
    pub(crate) populated: BTreeMap<Vec<u8>, u64>,
    /// Per member name, the entries that carry the member, in a format that counts them.
    pub(crate) present: BTreeMap<Vec<u8>, u64>,
    /// Per floor, the entries counted under it, in a format that counts them.
    pub(crate) under: BTreeMap<Vec<u8>, u64>,
}

impl Tally {
    pub(crate) fn add(&mut self, other: Tally) {
        for (counts, more) in [
            (&mut self.populated, other.populated),
            (&mut self.present, other.present),
            (&mut self.under, other.under),
        ] {
            for (name, count) in more {
                *counts.entry(name).or_insert(0) += count;
            }
        }
    }
}

pub(crate) trait Grouping: Sync {
    /// What a side's header lines before its first data line tell of its data lines.
    type Header: Default + Send + 'static;
    type Side: Side + 'static;

    /// The shapes under which every line is handed over, in the order one is reported
    /// when several hold.
    const SHAPES: &'static [&'static str];

    /// The key of the group of the lines past their own group's `limit`, in a format
    /// whose groups have one. No line's scan gives this key.
    const PLACED: Option<u64> = None;

    /// Whether `line` is a header line: handed over as it stands, outside every group.
    fn is_header(line: &[u8]) -> bool;

    /// Whether `line` is no line of the output's data: left out of the remainder, and
    /// counted among the lines the comparison numbers.
    fn is_blank(_line: &[u8]) -> bool {
        false
    }

    /// Takes a header line that stands before a side's first data line.
    fn header_line(header: &mut Self::Header, line: &[u8]);

    /// The scan of a side's data lines under its `header`; `data` says whether the side
    /// has a data line.
    fn side(&self, header: Self::Header, data: bool, reference: bool) -> Self::Side;

    /// The shape the two sides' headers show together, as its number in `SHAPES`.
    fn between(&self, expected: &Self::Side, actual: &Self::Side) -> Option<u32>;

    /// The shape the two sides' headers show together where no line shows one.
    fn trailing(&self, _expected: &Self::Side, _actual: &Self::Side) -> Option<u32> {
        None
    }

    /// The most lines of the group of `key` a side holds that the comparison names by the
    /// group's key: a line past them is a line of the group of `PLACED`, and its own
    /// group is handed over. `None` where every line of a group is named by its key.
    fn limit(&self, _key: u64) -> Option<u64> {
        None
    }

    /// Whether the lines of `Scanned::strict` hand their groups over, `carried` being
    /// what the lines of the two sides carry.
    fn strict(&self, _carried: u64) -> bool {
        false
    }
}

pub(crate) trait Side: Send + Sync {
    /// What a scan carries from one data line of a block to the next.
    type Memo: Default;

    /// The shape this side's header lines show alone, as its number in `SHAPES`;
    /// `header_after_data` says whether a header line follows a data line.
    fn shape(&self, header_after_data: bool) -> Option<u32>;

    fn scan(&self, memo: &mut Self::Memo, line: &[u8]) -> Scanned;

    /// The entries the comparison reads from a data line.
    fn entries(&self, _line: &[u8]) -> u32 {
        1
    }

    /// Appends to `text` what a kept-back line is held against the other side's by: one
    /// line of text, ended by `\n`, that two lines share exactly when the comparison
    /// reads the same entries from both.
    fn kept_text(&self, line: &[u8], text: &mut Vec<u8>) {
        text.extend_from_slice(line);
        text.push(b'\n');
    }

    /// Adds to `tally` the counts of the lines among `kept`, each the place of a kept-back
    /// line in `block`, and appends to `entries` the entries of each of them, in order,
    /// where a line may be more entries than one.
    fn populated(
        &self,
        block: &[u8],
        kept: &[Range<usize>],
        tally: &mut Tally,
        entries: &mut Vec<u32>,
    );
}
