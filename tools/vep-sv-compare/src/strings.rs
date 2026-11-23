// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The annotation tuple and the sets of tuples every count and mask works on.
//!
//! A tuple is `(Location, Allele, Feature, Feature_type, Consequence)` after normalisation.
//! Each distinct column text of a run is held once in [`Strings`], and a column of a tuple
//! is the rank of its text in the bytewise order of all texts, so a [`Tuple`] is five
//! integers whose derived order is the order of the five texts: the order `sorted()` gives
//! tuples of `str` in Python. A [`TupleSet`] is the sorted, duplicate-free list of its
//! tuples, so iteration is in that order at any thread count.

use std::cmp::Ordering;

/// The rank of a column text in [`Strings`].
pub type Id = u32;

/// One annotation tuple, each column the [`Id`] of its text.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Tuple(pub [Id; 5]);

impl Tuple {
    pub fn location(self) -> Id {
        self.0[0]
    }

    pub fn allele(self) -> Id {
        self.0[1]
    }

    pub fn feature(self) -> Id {
        self.0[2]
    }

    pub fn feature_type(self) -> Id {
        self.0[3]
    }

    pub fn consequence(self) -> Id {
        self.0[4]
    }

    /// The first four columns: the key a pair of one-sided tuples shares.
    pub fn key4(self) -> [Id; 4] {
        [self.0[0], self.0[1], self.0[2], self.0[3]]
    }

    /// `(Location, Allele)`: one allele of one record as an output prints it.
    pub fn variant(self) -> (Id, Id) {
        (self.0[0], self.0[1])
    }
}

/// Every distinct column text of a run, sorted bytewise; a text's [`Id`] is its index.
#[derive(Debug, Default)]
pub struct Strings {
    sorted: Vec<Box<str>>,
}

impl Strings {
    /// The table of the given texts; duplicates collapse.
    pub fn from_texts<I, S>(texts: I) -> Strings
    where
        I: IntoIterator<Item = S>,
        S: Into<Box<str>>,
    {
        let mut sorted: Vec<Box<str>> = texts.into_iter().map(Into::into).collect();
        sorted.sort_unstable();
        sorted.dedup();
        assert!(
            sorted.len() <= Id::MAX as usize,
            "more distinct column texts than an Id counts"
        );
        Strings { sorted }
    }

    /// The table over a list that is already sorted bytewise and duplicate-free.
    pub fn from_sorted_unique(sorted: Vec<Box<str>>) -> Strings {
        debug_assert!(sorted.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(
            sorted.len() <= Id::MAX as usize,
            "more distinct column texts than an Id counts"
        );
        Strings { sorted }
    }

    pub fn len(&self) -> usize {
        self.sorted.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sorted.is_empty()
    }

    /// The id of a text, or `None` when no tuple of the run carries it in any column.
    pub fn get(&self, text: &str) -> Option<Id> {
        self.sorted
            .binary_search_by(|held| (**held).cmp(text))
            .ok()
            .map(|index| index as Id)
    }

    pub fn text(&self, id: Id) -> &str {
        &self.sorted[id as usize]
    }

    /// The tuple of five texts, or `None` when one of them is not in the table (no tuple of
    /// the run equals it then).
    pub fn tuple(&self, columns: [&str; 5]) -> Option<Tuple> {
        Some(Tuple([
            self.get(columns[0])?,
            self.get(columns[1])?,
            self.get(columns[2])?,
            self.get(columns[3])?,
            self.get(columns[4])?,
        ]))
    }

    /// The five texts of a tuple.
    pub fn resolve(&self, tuple: Tuple) -> [&str; 5] {
        tuple.0.map(|id| self.text(id))
    }

    /// Every text in id order.
    pub fn iter(&self) -> impl Iterator<Item = &str> {
        self.sorted.iter().map(|text| &**text)
    }
}

/// A set of tuples: sorted, duplicate-free, immutable.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TupleSet {
    sorted: Vec<Tuple>,
}

impl TupleSet {
    pub fn new() -> TupleSet {
        TupleSet::default()
    }

    /// The set of the given tuples in any order, duplicates collapsing.
    pub fn from_unsorted(mut tuples: Vec<Tuple>) -> TupleSet {
        tuples.sort_unstable();
        tuples.dedup();
        TupleSet { sorted: tuples }
    }

    /// The set over a list that is already sorted and duplicate-free.
    pub fn from_sorted_unique(sorted: Vec<Tuple>) -> TupleSet {
        debug_assert!(sorted.windows(2).all(|pair| pair[0] < pair[1]));
        TupleSet { sorted }
    }

    pub fn len(&self) -> usize {
        self.sorted.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sorted.is_empty()
    }

    pub fn contains(&self, tuple: Tuple) -> bool {
        self.sorted.binary_search(&tuple).is_ok()
    }

    /// The tuples in order.
    pub fn iter(&self) -> impl Iterator<Item = Tuple> + '_ {
        self.sorted.iter().copied()
    }

    pub fn as_slice(&self) -> &[Tuple] {
        &self.sorted
    }

    /// The tuples whose first columns equal `prefix` (one to five leading columns), a
    /// contiguous run of the sorted list: every tuple at a Location, at a
    /// `(Location, Allele)`, and so on.
    pub fn with_prefix(&self, prefix: &[Id]) -> &[Tuple] {
        let width = prefix.len().min(5);
        let prefix = &prefix[..width];
        let start = self
            .sorted
            .partition_point(|tuple| tuple.0[..width] < *prefix);
        let run = self.sorted[start..].partition_point(|tuple| tuple.0[..width] == *prefix);
        &self.sorted[start..start + run]
    }

    /// `self | other`.
    pub fn union(&self, other: &TupleSet) -> TupleSet {
        let (a, b) = (&self.sorted, &other.sorted);
        let mut out = Vec::with_capacity(a.len().max(b.len()));
        let (mut i, mut j) = (0, 0);
        while i < a.len() && j < b.len() {
            match a[i].cmp(&b[j]) {
                Ordering::Less => {
                    out.push(a[i]);
                    i += 1;
                }
                Ordering::Greater => {
                    out.push(b[j]);
                    j += 1;
                }
                Ordering::Equal => {
                    out.push(a[i]);
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&a[i..]);
        out.extend_from_slice(&b[j..]);
        TupleSet { sorted: out }
    }

    /// The union of any number of sets.
    pub fn union_all<'a>(sets: impl IntoIterator<Item = &'a TupleSet>) -> TupleSet {
        let mut all: Vec<Tuple> = Vec::new();
        for set in sets {
            all.extend_from_slice(&set.sorted);
        }
        TupleSet::from_unsorted(all)
    }

    /// `self - other`.
    pub fn difference(&self, other: &TupleSet) -> TupleSet {
        let (a, b) = (&self.sorted, &other.sorted);
        let mut out = Vec::new();
        let mut j = 0;
        for &tuple in a {
            while j < b.len() && b[j] < tuple {
                j += 1;
            }
            if j == b.len() || b[j] != tuple {
                out.push(tuple);
            }
        }
        TupleSet { sorted: out }
    }

    /// `self & other`.
    pub fn intersection(&self, other: &TupleSet) -> TupleSet {
        let (a, b) = (&self.sorted, &other.sorted);
        let mut out = Vec::new();
        let mut j = 0;
        for &tuple in a {
            while j < b.len() && b[j] < tuple {
                j += 1;
            }
            if j < b.len() && b[j] == tuple {
                out.push(tuple);
            }
        }
        TupleSet { sorted: out }
    }

    /// `len(self & other)` without building the set.
    pub fn intersection_len(&self, other: &TupleSet) -> usize {
        let (a, b) = (&self.sorted, &other.sorted);
        let mut count = 0;
        let mut j = 0;
        for &tuple in a {
            while j < b.len() && b[j] < tuple {
                j += 1;
            }
            if j < b.len() && b[j] == tuple {
                count += 1;
            }
        }
        count
    }

    /// `self <= other`.
    pub fn is_subset(&self, other: &TupleSet) -> bool {
        self.intersection_len(other) == self.len()
    }
}

impl FromIterator<Tuple> for TupleSet {
    fn from_iter<I: IntoIterator<Item = Tuple>>(iter: I) -> TupleSet {
        TupleSet::from_unsorted(iter.into_iter().collect())
    }
}

impl<'a> IntoIterator for &'a TupleSet {
    type Item = Tuple;
    type IntoIter = std::iter::Copied<std::slice::Iter<'a, Tuple>>;

    fn into_iter(self) -> Self::IntoIter {
        self.sorted.iter().copied()
    }
}

/// The two sides of a comparison in one table of texts: the reference's tuples (`perl`), the
/// scored engine's (`rust`) and the tuples only one side holds.
pub struct Sides<'a> {
    pub strings: &'a Strings,
    pub perl: &'a TupleSet,
    pub rust: &'a TupleSet,
    /// `perl - rust`.
    pub perl_only: TupleSet,
    /// `rust - perl`.
    pub rust_only: TupleSet,
}

impl<'a> Sides<'a> {
    pub fn new(strings: &'a Strings, perl: &'a TupleSet, rust: &'a TupleSet) -> Sides<'a> {
        Sides {
            strings,
            perl,
            rust,
            perl_only: perl.difference(rust),
            rust_only: rust.difference(perl),
        }
    }

    /// The five texts of a tuple.
    pub fn resolve(&self, tuple: Tuple) -> [&'a str; 5] {
        self.strings.resolve(tuple)
    }
}

/// Two tuple sets written as texts, in one table: what a test of a mask starts from.
pub struct OwnedSides {
    pub strings: Strings,
    pub perl: TupleSet,
    pub rust: TupleSet,
}

impl OwnedSides {
    /// The sides whose tuples are the given rows of five texts.
    pub fn from_rows(perl: &[[&str; 5]], rust: &[[&str; 5]]) -> OwnedSides {
        OwnedSides::from_rows_with(perl, rust, &[])
    }

    /// As [`OwnedSides::from_rows`], the table also holding `extra` texts no tuple carries, so
    /// that a test can name a tuple of neither side.
    pub fn from_rows_with(perl: &[[&str; 5]], rust: &[[&str; 5]], extra: &[&str]) -> OwnedSides {
        let strings = Strings::from_texts(
            perl.iter()
                .chain(rust)
                .flat_map(|row| row.iter().copied())
                .chain(extra.iter().copied()),
        );
        let set = |rows: &[[&str; 5]]| {
            rows.iter()
                .map(|row| strings.tuple(*row).expect("every text is in the table"))
                .collect::<TupleSet>()
        };
        let (perl, rust) = (set(perl), set(rust));
        OwnedSides {
            strings,
            perl,
            rust,
        }
    }

    pub fn sides(&self) -> Sides<'_> {
        Sides::new(&self.strings, &self.perl, &self.rust)
    }

    /// A set of this table's tuples from rows of texts, for an expected value; a row with a
    /// text the table lacks is a tuple of neither side and panics.
    pub fn set(&self, rows: &[[&str; 5]]) -> TupleSet {
        rows.iter()
            .map(|row| {
                self.strings
                    .tuple(*row)
                    .expect("an expected tuple names a text the table holds")
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_tuple_orders_as_its_five_texts() {
        let rows = [
            [
                "21:9-10",
                "deletion",
                "ENST2",
                "Transcript",
                "intron_variant",
            ],
            ["21:10", "deletion", "ENST1", "Transcript", "intron_variant"],
            ["21:10", "-", "ENST1", "Transcript", "intron_variant"],
        ];
        let owned = OwnedSides::from_rows(&rows, &[]);
        let mut expected: Vec<[&str; 5]> = rows.to_vec();
        expected.sort_unstable();
        let resolved: Vec<[&str; 5]> = owned
            .perl
            .iter()
            .map(|t| owned.strings.resolve(t))
            .collect();
        assert_eq!(resolved, expected);
    }

    #[test]
    fn set_algebra_matches_the_definitions() {
        let a = [
            ["1:20", "A", "f1", "Transcript", "x"],
            ["1:30", "A", "f1", "Transcript", "x"],
            ["1:40", "A", "f1", "Transcript", "x"],
        ];
        let b = [
            ["1:30", "A", "f1", "Transcript", "x"],
            ["1:50", "A", "f1", "Transcript", "x"],
        ];
        let owned = OwnedSides::from_rows(&a, &b);
        let sides = owned.sides();
        assert_eq!(sides.perl_only, owned.set(&[a[0], a[2]]));
        assert_eq!(sides.rust_only, owned.set(&[b[1]]));
        assert_eq!(owned.perl.intersection(&owned.rust), owned.set(&[a[1]]));
        assert_eq!(owned.perl.intersection_len(&owned.rust), 1);
        assert_eq!(owned.perl.union(&owned.rust).len(), 4);
        assert_eq!(
            TupleSet::union_all([&owned.perl, &owned.rust]),
            owned.perl.union(&owned.rust)
        );
        assert!(sides.perl_only.is_subset(&owned.perl));
        assert!(!owned.perl.is_subset(&owned.rust));
    }

    #[test]
    fn a_prefix_selects_a_contiguous_run() {
        let rows = [
            ["1:20", "A", "f1", "Transcript", "x"],
            ["1:20", "A", "f2", "Transcript", "x"],
            ["1:20", "C", "f1", "Transcript", "x"],
            ["1:30", "A", "f1", "Transcript", "x"],
        ];
        let owned = OwnedSides::from_rows(&rows, &[]);
        let id = |text| owned.strings.get(text).unwrap();
        assert_eq!(owned.perl.with_prefix(&[id("1:20")]).len(), 3);
        assert_eq!(owned.perl.with_prefix(&[id("1:20"), id("A")]).len(), 2);
        assert_eq!(owned.perl.with_prefix(&[id("1:30"), id("C")]).len(), 0);
        assert!(owned.strings.get("1:40").is_none());
    }
}
