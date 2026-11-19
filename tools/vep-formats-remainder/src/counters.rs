// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Per field name, the kept-back lines that carry a value, as the comparator's
//! `populated_counts` counts its entries: every column, and every `key=value` pair of the
//! `Extra` column as a field of its own, the last pair of a name standing for every pair
//! of that name before it and for the column of that name; `-` and the empty text are no
//! value.

use crate::input::Side;
use crate::layout::{Layout, EXTRA};
use memchr::{memchr, memchr_iter};
use std::collections::{BTreeMap, HashMap};

fn valued(value: &[u8]) -> bool {
    !value.is_empty() && value != b"-"
}

struct Counter<'a> {
    ids: HashMap<&'a [u8], u32>,
    names: Vec<&'a [u8]>,
    counts: Vec<u64>,
    /// Per name, the line that last carried it and its place in that line's `fields`.
    seen: Vec<(u64, u32)>,
    fields: Vec<(u32, bool)>,
    column_ids: Vec<u32>,
    serial: u64,
}

impl<'a> Counter<'a> {
    fn new(layout: &Layout<'a>) -> Counter<'a> {
        let mut counter = Counter {
            ids: HashMap::new(),
            names: Vec::new(),
            counts: Vec::new(),
            seen: Vec::new(),
            fields: Vec::new(),
            column_ids: Vec::new(),
            serial: 0,
        };
        counter.column_ids = layout.names.iter().map(|&n| counter.id(n)).collect();
        counter
    }

    fn id(&mut self, name: &'a [u8]) -> u32 {
        if let Some(&id) = self.ids.get(name) {
            return id;
        }
        let id = self.names.len() as u32;
        self.ids.insert(name, id);
        self.names.push(name);
        self.counts.push(0);
        self.seen.push((0, 0));
        id
    }

    fn field(&mut self, id: u32, valued: bool) {
        let seen = &mut self.seen[id as usize];
        if seen.0 == self.serial {
            self.fields[seen.1 as usize].1 = valued;
        } else {
            *seen = (self.serial, self.fields.len() as u32);
            self.fields.push((id, valued));
        }
    }

    fn line(&mut self, text: &'a [u8], extra: Option<usize>) {
        self.serial += 1;
        self.fields.clear();
        let mut extra_value: Option<&'a [u8]> = None;
        let mut start = 0;
        let ends = memchr_iter(b'\t', text).chain(std::iter::once(text.len()));
        for (column, end) in ends.enumerate() {
            let value = &text[start..end];
            if let Some(&id) = self.column_ids.get(column) {
                self.field(id, valued(value));
            }
            if extra == Some(column) {
                extra_value = Some(value);
            }
            start = end + 1;
        }
        // The comparator leaves an absent `Extra` (`-`) unread and, after the pairs, sets
        // the field `Extra` to the column's own text again.
        let Some(extra_value) = extra_value.filter(|v| *v != b"-") else {
            return;
        };
        for pair in extra_value.split(|&b| b == b';') {
            let (name, value) = match memchr(b'=', pair) {
                Some(i) => (&pair[..i], &pair[i + 1..]),
                None => (pair, &pair[pair.len()..]),
            };
            if name != EXTRA {
                let id = self.id(name);
                self.field(id, valued(value));
            }
        }
    }

    fn finish_line(&mut self) {
        for &(id, valued) in &self.fields {
            if valued {
                self.counts[id as usize] += 1;
            }
        }
    }
}

/// The populated counts of the kept-back lines of `side`, by field name, a name with no
/// valued line left out.
pub fn populated<'a>(
    side: &'a Side,
    kept: &[bool],
    layout: &Layout<'a>,
    threads: usize,
) -> BTreeMap<Vec<u8>, u64> {
    let per_worker = side.lines.len().div_ceil(threads.max(1)).max(1);
    let parts: Vec<Vec<(Vec<u8>, u64)>> = std::thread::scope(|scope| {
        let workers: Vec<_> = side
            .lines
            .chunks(per_worker)
            .zip(kept.chunks(per_worker))
            .map(|(lines, kept)| {
                scope.spawn(move || {
                    let mut counter = Counter::new(layout);
                    for (line, _) in lines.iter().zip(kept).filter(|(_, k)| **k) {
                        counter.line(side.text(line), layout.extra);
                        counter.finish_line();
                    }
                    counter
                        .names
                        .iter()
                        .zip(&counter.counts)
                        .filter(|(_, n)| **n > 0)
                        .map(|(name, n)| (name.to_vec(), *n))
                        .collect()
                })
            })
            .collect();
        workers
            .into_iter()
            .map(|w| w.join().expect("counter thread"))
            .collect()
    });
    let mut total = BTreeMap::new();
    for (name, n) in parts.into_iter().flatten() {
        *total.entry(name).or_insert(0) += n;
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(text: &str) -> Vec<(String, u64)> {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-remainder-counters-{}-{:x}",
            std::process::id(),
            text.len() * 31 + text.bytes().map(usize::from).sum::<usize>()
        ));
        std::fs::write(&path, text).unwrap();
        let side = Side::read(&path, 1).unwrap();
        std::fs::remove_file(&path).unwrap();
        let layout = Layout::of(&side).unwrap();
        let kept: Vec<bool> = side.lines.iter().map(|l| l.flags == 0).collect();
        populated(&side, &kept, &layout, 2)
            .into_iter()
            .map(|(name, n)| (String::from_utf8(name).unwrap(), n))
            .collect()
    }

    fn pairs(expected: &[(&str, u64)]) -> Vec<(String, u64)> {
        expected
            .iter()
            .map(|(n, c)| ((*n).to_string(), *c))
            .collect()
    }

    #[test]
    fn columns_and_extra_pairs_count_once_per_line() {
        let text =
            "#Location\tGene\tExtra\nl1\tg\tIMPACT=LOW;STRAND=1\nl2\t-\tIMPACT=HIGH\nl3\t\t-\n";
        assert_eq!(
            counts(text),
            pairs(&[
                ("Extra", 2),
                ("Gene", 1),
                ("IMPACT", 2),
                ("Location", 3),
                ("STRAND", 1)
            ])
        );
    }

    #[test]
    fn a_dash_an_empty_value_and_a_bare_name_are_no_value() {
        let text = "#Location\tExtra\nl1\tA=-;B=;C;D=x\n";
        assert_eq!(
            counts(text),
            pairs(&[("D", 1), ("Extra", 1), ("Location", 1)])
        );
    }

    #[test]
    fn the_last_pair_of_a_name_stands_for_the_pairs_before_it_and_the_column() {
        let text = "#Location\tGene\tExtra\nl1\tg\tA=x;A=-;Gene=-;B=-;B=y\n";
        assert_eq!(
            counts(text),
            pairs(&[("B", 1), ("Extra", 1), ("Location", 1)])
        );
    }

    #[test]
    fn an_empty_extra_and_a_pair_named_extra_leave_the_column_s_own_text() {
        let text = "#Location\tExtra\nl1\t\nl2\tExtra=-;A=1\n";
        assert_eq!(
            counts(text),
            pairs(&[("A", 1), ("Extra", 1), ("Location", 2)])
        );
    }
}
