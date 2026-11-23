// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class of every row, from the numbers of its values.
//!
//! Reference-only rows are paired with scored-only rows in file order, each row once.
//! First by start, Allele, Feature and consequence set where the two Locations differ
//! (`coordinate_only`); then, among the rows left, by Location, Allele and Feature
//! (`consequence_swap`). A row left unpaired is `transcript_selection` where the other
//! engine wrote its start and Allele under other Features alone, and `true_missing`
//! (reference-only) or `true_extra` (scored-only) otherwise.

use crate::intern::{Bits, FoldBuild};
use crate::scan::Keys;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Class {
    CoordinateOnly,
    ConsequenceSwap,
    TranscriptSelection,
    TrueMissing,
    TrueExtra,
    /// The value the table read holds in its own `mismatch_type` column.
    AsRead,
}

impl Class {
    /// The classes the tool names, in the order of `index`.
    pub const NAMED: [Class; 5] = [
        Class::CoordinateOnly,
        Class::ConsequenceSwap,
        Class::TranscriptSelection,
        Class::TrueMissing,
        Class::TrueExtra,
    ];

    /// The text written for the class; None for a value read from the table.
    pub fn name(self) -> Option<&'static str> {
        match self {
            Class::CoordinateOnly => Some("coordinate_only"),
            Class::ConsequenceSwap => Some("consequence_swap"),
            Class::TranscriptSelection => Some("transcript_selection"),
            Class::TrueMissing => Some("true_missing"),
            Class::TrueExtra => Some("true_extra"),
            Class::AsRead => None,
        }
    }
}

/// The fewest rows a thread fills in with their classes.
const PART_ROWS: usize = 1 << 16;

/// The scored-only rows under each key in file order, and how far the reference-only rows
/// of a key have looked through them.
struct Candidates<const N: usize> {
    groups: HashMap<[u32; N], u32, FoldBuild>,
    /// Where the rows of each group start in `rows`, and the end of the last group.
    first: Vec<u32>,
    rows: Vec<u32>,
    /// Where in `rows` the next reference-only row of each group starts to look.
    next: Vec<u32>,
}

impl<const N: usize> Candidates<N> {
    /// Gathers the scored-only rows under their keys. `asked` holds the keys of the
    /// reference-only rows, the only keys ever looked up: a scored-only row under another
    /// key is then left out. Without it every key is held.
    fn of(
        asked: Option<impl Iterator<Item = [u32; N]>>,
        scored: impl Iterator<Item = (usize, [u32; N])>,
    ) -> Candidates<N> {
        let mut groups: HashMap<[u32; N], u32, FoldBuild> = HashMap::default();
        let mut sizes: Vec<u32> = Vec::new();
        let every_key = asked.is_none();
        for key in asked.into_iter().flatten() {
            if let Entry::Vacant(place) = groups.entry(key) {
                place.insert(sizes.len() as u32);
                sizes.push(0);
            }
        }
        let mut members: Vec<(u32, u32)> = Vec::new();
        for (row, key) in scored {
            let group = if every_key {
                let group = *groups.entry(key).or_insert(sizes.len() as u32);
                if group as usize == sizes.len() {
                    sizes.push(0);
                }
                group
            } else {
                match groups.get(&key) {
                    Some(group) => *group,
                    None => continue,
                }
            };
            sizes[group as usize] += 1;
            members.push((group, row as u32));
        }
        let mut first = Vec::with_capacity(sizes.len() + 1);
        let mut total = 0;
        for size in &sizes {
            first.push(total);
            total += size;
        }
        first.push(total);
        let mut place = first.clone();
        let mut rows = vec![0; members.len()];
        for (group, row) in members {
            rows[place[group as usize] as usize] = row;
            place[group as usize] += 1;
        }
        // A copy of `first`, where a count from zero would do: an array made of zeroes and
        // read before it is written pays a page copy at the first write to each page.
        let next = first[..sizes.len()].to_vec();
        Candidates {
            groups,
            first,
            rows,
            next,
        }
    }

    /// The group of a key, where the table holds the key.
    fn group(&self, key: &[u32; N]) -> Option<usize> {
        self.groups.get(key).map(|group| *group as usize)
    }
}

/// Whether the reference-only rows are the smaller half of the table: the smaller half is
/// held whole, and of the larger only what the smaller asks for.
fn reference_is_smaller(keys: &Keys) -> bool {
    keys.reference.count() * 2 <= keys.rows
}

/// The two pairings: the rows paired on everything but the end of their Location, then
/// the rows paired on Location, Allele and Feature. With `asked`, the scored-only rows are
/// held under the keys of the reference-only rows alone; the pairs are the same either way.
fn pair(keys: &Keys, asked: bool) -> (Bits, Bits) {
    let reference_rows = || (0..keys.rows).filter(|row| keys.reference.get(*row));
    let scored_rows = || (0..keys.rows).filter(|row| !keys.reference.get(*row));

    let mut coordinate = Bits::new(keys.rows);
    {
        let key = |row: usize| {
            [
                keys.start(row),
                keys.allele[row],
                keys.feature[row],
                keys.consequence[row],
            ]
        };
        let mut candidates = Candidates::of(
            asked.then(|| reference_rows().map(key)),
            scored_rows().map(|row| (row, key(row))),
        );
        for row in reference_rows() {
            let Some(group) = candidates.group(&key(row)) else {
                continue;
            };
            // A scored row whose Location is the reference row's own is passed for good, and
            // a reference row that finds no partner leaves none for the rows of its key
            // behind it.
            let end = candidates.first[group + 1];
            let mut at = candidates.next[group];
            while at < end {
                let other = candidates.rows[at as usize] as usize;
                at += 1;
                if !coordinate.get(other) && keys.location[other] != keys.location[row] {
                    coordinate.set(row);
                    coordinate.set(other);
                    break;
                }
            }
            candidates.next[group] = at;
        }
    }

    let mut swap = Bits::new(keys.rows);
    let key = |row: usize| [keys.location[row], keys.allele[row], keys.feature[row]];
    let unpaired = |row: &usize| !coordinate.get(*row);
    let mut candidates = Candidates::of(
        asked.then(|| reference_rows().filter(unpaired).map(key)),
        scored_rows().filter(unpaired).map(|row| (row, key(row))),
    );
    for row in reference_rows().filter(unpaired) {
        let Some(group) = candidates.group(&key(row)) else {
            continue;
        };
        let at = candidates.next[group];
        if at < candidates.first[group + 1] {
            swap.set(row);
            swap.set(candidates.rows[at as usize] as usize);
            candidates.next[group] = at + 1;
        }
    }
    (coordinate, swap)
}

/// What one engine's rows hold: every start and Allele, and every start, Allele and
/// Feature.
struct Side {
    variants: HashSet<[u32; 2], FoldBuild>,
    features: HashSet<[u32; 3], FoldBuild>,
}

/// What the rows of one engine hold; with `within`, only at the variants that side holds.
fn side(keys: &Keys, reference: bool, within: Option<&Side>) -> Side {
    let mut side = Side {
        variants: HashSet::default(),
        features: HashSet::default(),
    };
    for row in (0..keys.rows).filter(|row| keys.reference.get(*row) == reference) {
        let variant = [keys.start(row), keys.allele[row]];
        if within.is_some_and(|other| !other.variants.contains(&variant)) {
            continue;
        }
        side.variants.insert(variant);
        side.features
            .insert([variant[0], variant[1], keys.feature[row]]);
    }
    side
}

/// What each engine wrote, the reference first. A row asks the other side about its own
/// variant alone, so one side is held whole, the reference under `reference_whole`, and
/// the other only at the variants of the first; a row is answered the same either way.
fn sides(keys: &Keys, reference_whole: bool) -> (Side, Side) {
    if reference_whole {
        let reference = side(keys, true, None);
        let scored = side(keys, false, Some(&reference));
        (reference, scored)
    } else {
        let scored = side(keys, false, None);
        let reference = side(keys, true, Some(&scored));
        (reference, scored)
    }
}

/// The class of every row in file order. Under `as_read` the table read holds a
/// `mismatch_type` column of its own: a row keeps the value it holds there unless the
/// second pairing takes it. `threads` bounds the threads at work.
pub fn classify(keys: &Keys, as_read: bool, threads: usize) -> Vec<Class> {
    classify_holding(keys, as_read, threads, reference_is_smaller(keys))
}

/// `classify`, with the half of the table that is held whole named: the reference-only
/// rows under `reference_whole`. The classes do not depend on it.
fn classify_holding(
    keys: &Keys,
    as_read: bool,
    threads: usize,
    reference_whole: bool,
) -> Vec<Class> {
    let ((coordinate, swap), sides) = std::thread::scope(|scope| {
        if as_read {
            (pair(keys, reference_whole), None)
        } else if threads >= 2 {
            let pairing = scope.spawn(|| pair(keys, reference_whole));
            let sides = sides(keys, reference_whole);
            (pairing.join().expect("the pairing thread"), Some(sides))
        } else {
            let pairs = pair(keys, reference_whole);
            (pairs, Some(sides(keys, reference_whole)))
        }
    });
    let class_of = |row: usize| {
        if swap.get(row) {
            return Class::ConsequenceSwap;
        }
        let Some((reference_side, scored_side)) = &sides else {
            return Class::AsRead;
        };
        if coordinate.get(row) {
            return Class::CoordinateOnly;
        }
        let reference = keys.reference.get(row);
        let other = if reference {
            scored_side
        } else {
            reference_side
        };
        let (start, allele) = (keys.start(row), keys.allele[row]);
        let elsewhere = other.variants.contains(&[start, allele])
            && !other.features.contains(&[start, allele, keys.feature[row]]);
        match (elsewhere, reference) {
            (true, _) => Class::TranscriptSelection,
            (false, true) => Class::TrueMissing,
            (false, false) => Class::TrueExtra,
        }
    };
    let fill = |part: &mut [Class], from: usize| {
        for (offset, class) in part.iter_mut().enumerate() {
            *class = class_of(from + offset);
        }
    };
    let mut classes = vec![Class::AsRead; keys.rows];
    let part_rows = keys.rows.div_ceil(threads.max(1)).max(PART_ROWS);
    let mut parts = classes.chunks_mut(part_rows).enumerate();
    let first = parts.next();
    std::thread::scope(|scope| {
        for (index, part) in parts {
            scope.spawn(move || fill(part, index * part_rows));
        }
        if let Some((_, part)) = first {
            fill(part, 0);
        }
    });
    classes
}

#[cfg(test)]
#[path = "classes_tests.rs"]
mod tests;
