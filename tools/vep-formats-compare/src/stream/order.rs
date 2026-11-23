// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two order conditions of what is kept back.
//!
//! The comparison reads the order of the keys both sides share and reports the first
//! place where the two orders part. A group is therefore handed over when the kept-back
//! lines would not stand in one order on both sides, and when it stands between two lines
//! of a handed-over group that both sides hold: the first place where the shared keys of
//! the remainder part is then that place in the whole outputs, moved by the kept-back
//! lines before it.

use super::table::{Bits, Chunked, Row};

/// Rounds of the second order condition before every line is handed over: each round
/// hands over the kept-back groups that stand inside a handed-over group's span.
pub(crate) const ORDER_ROUNDS: usize = 64;

/// A place in a side's lines: a row and the lines of it already passed.
struct Place<'a> {
    rows: &'a Chunked<Row>,
    row: usize,
    passed: u32,
}

impl<'a> Place<'a> {
    /// The first line of a side.
    fn start(rows: &'a Chunked<Row>) -> Place<'a> {
        Place {
            rows,
            row: 0,
            passed: 0,
        }
    }

    /// Moves to the next line of a kept-back group; `None` at the end of the side.
    fn kept(&mut self, kept: &Bits) -> Option<Row> {
        while self.row < self.rows.len() {
            let row = self.rows.get(self.row);
            if kept.get(row.group) {
                return Some(row);
            }
            self.row += 1;
            self.passed = 0;
        }
        None
    }

    fn pass(&mut self, row: Row, lines: u32) {
        self.passed += lines;
        if self.passed == row.lines {
            self.row += 1;
            self.passed = 0;
        }
    }
}

/// The first order condition: walks the kept-back lines of both sides together and hands
/// over the two groups met wherever the sides name different groups.
pub(crate) fn one_order(expected: &Chunked<Row>, actual: &Chunked<Row>, kept: &mut Bits) {
    let (mut e, mut a) = (Place::start(expected), Place::start(actual));
    while let (Some(row_e), Some(row_a)) = (e.kept(kept), a.kept(kept)) {
        if row_e.group == row_a.group {
            let lines = (row_e.lines - e.passed).min(row_a.lines - a.passed);
            e.pass(row_e, lines);
            a.pass(row_a, lines);
        } else {
            kept.clear(row_e.group);
            kept.clear(row_a.group);
        }
    }
    // Both places stand at the next kept-back line of their side, or at its end.
    for place in [e, a] {
        for row in place.rows.range(place.row, place.rows.len()) {
            kept.clear(row.group);
        }
    }
}

/// Whether the kept-back lines of the two sides stand in one order.
pub(crate) fn in_one_order(expected: &Chunked<Row>, actual: &Chunked<Row>, kept: &Bits) -> bool {
    let (mut e, mut a) = (Place::start(expected), Place::start(actual));
    loop {
        match (e.kept(kept), a.kept(kept)) {
            (Some(row_e), Some(row_a)) if row_e.group == row_a.group => {
                let lines = (row_e.lines - e.passed).min(row_a.lines - a.passed);
                e.pass(row_e, lines);
                a.pass(row_a, lines);
            }
            (None, None) => return true,
            _ => return false,
        }
    }
}

/// What the rounds of the second order condition share: per group, the fewest and the
/// most kept-back lines that precede one of its lines, on either side.
pub(crate) struct Spans {
    low: Vec<u64>,
    high: Vec<u64>,
}

impl Spans {
    pub(crate) fn new(groups: usize) -> Spans {
        Spans {
            low: vec![u64::MAX; groups],
            high: vec![0; groups],
        }
    }

    /// One round: hands over every kept-back group with a line between two lines of one
    /// handed-over group that both sides hold. Returns whether a handed-over group spans
    /// a kept-back line.
    pub(crate) fn close(
        &mut self,
        expected: &Chunked<Row>,
        actual: &Chunked<Row>,
        two_sided: &Bits,
        kept: &mut Bits,
    ) -> bool {
        self.low.fill(u64::MAX);
        self.high.fill(0);
        for side in [expected, actual] {
            let mut before = 0u64;
            for row in side.iter() {
                let group = row.group as usize;
                if kept.get(row.group) {
                    before += u64::from(row.lines);
                } else if two_sided.get(row.group) {
                    self.low[group] = self.low[group].min(before);
                    self.high[group] = self.high[group].max(before);
                }
            }
        }
        // The spans as changes, at a kept-back line's ordinal, in the number of spans
        // open over it.
        let mut changes: Vec<(u64, i32)> = Vec::new();
        for (&low, &high) in self.low.iter().zip(&self.high) {
            if low < high {
                changes.push((low, 1));
                changes.push((high, -1));
            }
        }
        if changes.is_empty() {
            return false;
        }
        changes.sort_unstable();
        let mut inside: Vec<u32> = Vec::new();
        let (mut open, mut next, mut line) = (0i64, 0usize, 0u64);
        for row in expected.iter() {
            if !kept.get(row.group) {
                continue;
            }
            let end = line + u64::from(row.lines);
            // The changes up to the row's first line, then those at each later line of it.
            let mut at = line;
            loop {
                while next < changes.len() && changes[next].0 <= at {
                    open += i64::from(changes[next].1);
                    next += 1;
                }
                if open > 0 {
                    inside.push(row.group);
                    break;
                }
                match changes.get(next) {
                    Some(&(place, _)) if place < end => at = place,
                    _ => break,
                }
            }
            line = end;
        }
        for group in inside {
            kept.clear(group);
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(side: &str) -> Chunked<Row> {
        let mut rows: Chunked<Row> = Chunked::new();
        let mut last: Option<Row> = None;
        for group in side.chars().map(|c| c as u32 - 'a' as u32) {
            match &mut last {
                Some(row) if row.group == group => row.lines += 1,
                _ => {
                    rows.extend_one(last);
                    last = Some(Row { group, lines: 1 });
                }
            }
        }
        rows.extend_one(last);
        rows
    }

    impl Chunked<Row> {
        fn extend_one(&mut self, row: Option<Row>) {
            if let Some(row) = row {
                self.push(row);
            }
        }
    }

    /// Groups by letter, one line each: which stay kept back under the order conditions
    /// when `handed` are handed over from the start.
    fn order(expected: &str, actual: &str, handed: &str) -> String {
        let id = |c: char| c as u32 - 'a' as u32;
        let groups = expected
            .chars()
            .chain(actual.chars())
            .map(id)
            .max()
            .map_or(0, |m| m as usize + 1);
        let count = |side: &str, g: u32| side.chars().filter(|c| id(*c) == g).count();
        let mut kept = Bits::new(groups);
        let mut two_sided = Bits::new(groups);
        for g in 0..groups as u32 {
            let c = char::from(b'a' + g as u8);
            if !handed.contains(c)
                && count(expected, g) > 0
                && count(expected, g) == count(actual, g)
            {
                kept.set(g);
            }
            if count(expected, g) > 0 && count(actual, g) > 0 {
                two_sided.set(g);
            }
        }
        let (e, a) = (rows(expected), rows(actual));
        one_order(&e, &a, &mut kept);
        assert!(in_one_order(&e, &a, &kept));
        let mut spans = Spans::new(groups);
        while spans.close(&e, &a, &two_sided, &mut kept) {}
        (0..groups as u32)
            .filter(|g| kept.get(*g))
            .map(|g| char::from(b'a' + g as u8))
            .collect()
    }

    #[test]
    fn groups_in_one_order_are_all_kept_back() {
        assert_eq!(order("abcd", "abcd", ""), "abcd");
        assert_eq!(order("abcd", "abcd", "b"), "acd");
    }

    #[test]
    fn groups_that_change_places_are_handed_over() {
        assert_eq!(order("abcd", "acbd", ""), "ad");
        // d ends the reference side and starts the candidate side: b and c stand inside its span.
        assert_eq!(order("abcd", "dabc", ""), "");
    }

    #[test]
    fn a_kept_back_group_inside_a_handed_over_group_s_span_is_handed_over() {
        // Group b stands apart on the reference side, with c between its lines.
        assert_eq!(order("abcbd", "abbcd", "b"), "ad");
        // The same on both sides: c is still between two lines of a handed-over group.
        assert_eq!(order("abcbd", "abcbd", "b"), "ad");
        assert_eq!(order("abcbd", "abcbd", ""), "abcd");
    }

    #[test]
    fn a_handed_over_group_at_two_places_hands_over_what_stands_between() {
        // b is before c on one side and after d on the other.
        assert_eq!(order("abcde", "acdbe", "b"), "ae");
        // A chain: handing c over leaves it around d on the candidate side.
        assert_eq!(order("abcbdce", "abbcdce", "b"), "ae");
    }

    #[test]
    fn a_one_sided_group_moves_nothing() {
        assert_eq!(order("axbc", "abcy", ""), "abc");
    }

    #[test]
    fn runs_of_several_lines_are_walked_line_by_line() {
        // A group's lines cut into other runs on the two sides, in one order all the same.
        assert_eq!(order("aabaa", "aabaa", ""), "ab");
        assert_eq!(order("aaab", "abaa", ""), "");
        // A one-sided group between a group's lines moves nothing.
        assert_eq!(order("aaxab", "ayaab", ""), "ab");
        // A span that opens inside a run of a kept-back group.
        assert_eq!(order("aabcc", "abacc", "b"), "c");
        assert_eq!(order("aabcc", "aabcc", "b"), "ac");
        assert!(!in_one_order(&rows("ab"), &rows("ba"), &{
            let mut kept = Bits::new(2);
            kept.set(0);
            kept.set(1);
            kept
        }));
    }
}
