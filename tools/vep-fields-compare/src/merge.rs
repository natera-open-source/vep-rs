// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The merge of the two sorted sides: duplicate-key collapse, the per-field and per-Extra-key
//! counts, and the mismatch shapes in first-observation order.

use memchr::{memchr, memchr_iter};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

pub const PLAIN_FIELDS: [&str; 8] = [
    "Uploaded_variation",
    "Gene",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
];
pub const ABSENT: &[u8] = b"<absent>";

const KEY_FIELDS: usize = 5;
const VALUE_FIELDS: usize = 9;
const READ_BUF: usize = 4 << 20;

#[derive(Default)]
pub struct FieldStats {
    pub compared: u64,
    pub equal: u64,
    /// `reference\tscored` to (count, first-observation ordinal); neither value holds a tab.
    shapes: HashMap<Box<[u8]>, (u64, u64)>,
}

pub struct Shape<'a> {
    pub reference: &'a [u8],
    pub scored: &'a [u8],
    pub count: u64,
}

impl FieldStats {
    fn observe(&mut self, reference: &[u8], scored: &[u8], scratch: &mut Vec<u8>) {
        self.compared += 1;
        if reference == scored {
            self.equal += 1;
            return;
        }
        scratch.clear();
        scratch.extend_from_slice(reference);
        scratch.push(b'\t');
        scratch.extend_from_slice(scored);
        if let Some(entry) = self.shapes.get_mut(scratch.as_slice()) {
            entry.0 += 1;
        } else {
            let ordinal = self.shapes.len() as u64;
            self.shapes
                .insert(scratch.clone().into_boxed_slice(), (1, ordinal));
        }
    }

    fn observe_equal(&mut self) {
        self.compared += 1;
        self.equal += 1;
    }

    pub fn distinct_mismatch_shapes(&self) -> usize {
        self.shapes.len()
    }

    /// The `n` most frequent shapes, equal counts in first-observation order, as
    /// `Counter.most_common(n)` orders them; `n <= 0` gives none.
    pub fn top_mismatches(&self, n: i64) -> Vec<Shape<'_>> {
        let mut all: Vec<(&[u8], u64, u64)> = self
            .shapes
            .iter()
            .map(|(k, &(count, ordinal))| (k.as_ref(), count, ordinal))
            .collect();
        all.sort_unstable_by(|(_, count_a, ordinal_a), (_, count_b, ordinal_b)| {
            count_b.cmp(count_a).then(ordinal_a.cmp(ordinal_b))
        });
        all.into_iter()
            .take(n.max(0) as usize)
            .map(|(k, count, _)| {
                let tab = memchr(b'\t', k).expect("shape separator");
                Shape {
                    reference: &k[..tab],
                    scored: &k[tab + 1..],
                    count,
                }
            })
            .collect()
    }
}

#[derive(Default)]
pub struct Totals {
    pub reference_keys: u64,
    pub scored_keys: u64,
    pub intersection: u64,
    pub reference_only: u64,
    pub scored_only: u64,
    pub reference_duplicate_keys: u64,
    pub scored_duplicate_keys: u64,
}

pub struct MergeResult {
    pub totals: Totals,
    /// In `PLAIN_FIELDS` order.
    pub fields: Vec<FieldStats>,
    /// Sorted by key.
    pub extra: Vec<(Vec<u8>, FieldStats)>,
}

/// One sorted side. `line` holds the first row of the current key's run and `lookahead` the
/// row that ended the run, so a run of equal keys collapses without re-reading;
/// `duplicate_keys` counts the rows the runs dropped.
struct Side {
    reader: BufReader<File>,
    line: Vec<u8>,
    lookahead: Vec<u8>,
    has_lookahead: bool,
    key_end: usize,
    duplicate_keys: u64,
}

fn read_row(reader: &mut BufReader<File>, buf: &mut Vec<u8>) -> io::Result<bool> {
    buf.clear();
    if reader.read_until(b'\n', buf)? == 0 {
        return Ok(false);
    }
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    Ok(true)
}

/// The end of the first five tab-separated fields, or the whole row when it has fewer.
fn key_end(row: &[u8]) -> usize {
    memchr_iter(b'\t', row)
        .nth(KEY_FIELDS - 1)
        .unwrap_or(row.len())
}

impl Side {
    fn open(path: &Path) -> io::Result<Side> {
        Ok(Side {
            reader: BufReader::with_capacity(READ_BUF, File::open(path)?),
            line: Vec::new(),
            lookahead: Vec::new(),
            has_lookahead: false,
            key_end: 0,
            duplicate_keys: 0,
        })
    }

    /// Advances to the next distinct key, counting the rows that repeat it under
    /// `duplicate_keys`; `false` at end of file.
    fn next_unique(&mut self) -> io::Result<bool> {
        if self.has_lookahead {
            std::mem::swap(&mut self.line, &mut self.lookahead);
            self.has_lookahead = false;
        } else if !read_row(&mut self.reader, &mut self.line)? {
            return Ok(false);
        }
        self.key_end = key_end(&self.line);
        while read_row(&mut self.reader, &mut self.lookahead)? {
            if key_end(&self.lookahead) == self.key_end
                && self.lookahead[..self.key_end] == self.line[..self.key_end]
            {
                self.duplicate_keys += 1;
                continue;
            }
            self.has_lookahead = true;
            break;
        }
        Ok(true)
    }

    fn key(&self) -> &[u8] {
        &self.line[..self.key_end]
    }

    fn values(&self) -> &[u8] {
        self.line.get(self.key_end + 1..).unwrap_or(b"")
    }
}

fn split_values(values: &[u8]) -> [&[u8]; VALUE_FIELDS] {
    let mut out = [b"" as &[u8]; VALUE_FIELDS];
    let mut start = 0;
    let mut n = 0;
    for tab in memchr_iter(b'\t', values) {
        out[n] = &values[start..tab];
        start = tab + 1;
        n += 1;
        if n == VALUE_FIELDS {
            return out;
        }
    }
    out[n] = &values[start..];
    out
}

/// A `KEY=VALUE` pair of the Extra column as byte ranges into the column.
#[derive(Clone, Copy)]
struct Pair {
    key: (usize, usize),
    value: (usize, usize),
}

fn slice(raw: &[u8], range: (usize, usize)) -> &[u8] {
    &raw[range.0..range.1]
}

/// `parse_extra`: `""` and `-` hold no pairs; parts split on `;`, empty parts skipped, a
/// bare token is a key with value `""`, a repeated key keeps its last value.
fn parse_extra(raw: &[u8], out: &mut Vec<Pair>) {
    out.clear();
    if raw.is_empty() || raw == b"-" {
        return;
    }
    let mut start = 0;
    for end in memchr_iter(b';', raw).chain(std::iter::once(raw.len())) {
        if end > start {
            let pair = match memchr(b'=', &raw[start..end]) {
                Some(i) => Pair {
                    key: (start, start + i),
                    value: (start + i + 1, end),
                },
                None => Pair {
                    key: (start, end),
                    value: (end, end),
                },
            };
            match out
                .iter_mut()
                .find(|p| slice(raw, p.key) == slice(raw, pair.key))
            {
                Some(existing) => existing.value = pair.value,
                None => out.push(pair),
            }
        }
        start = end + 1;
    }
}

#[derive(Default)]
struct ExtraScratch {
    reference: Vec<Pair>,
    scored: Vec<Pair>,
    /// (side, index) over both pair lists, sorted and deduplicated by key.
    keys: Vec<(u8, usize)>,
    shape: Vec<u8>,
}

fn stats_for<'a>(extra: &'a mut Vec<(Vec<u8>, FieldStats)>, key: &[u8]) -> &'a mut FieldStats {
    let index = match extra.iter().position(|(k, _)| k.as_slice() == key) {
        Some(i) => i,
        None => {
            extra.push((key.to_vec(), FieldStats::default()));
            extra.len() - 1
        }
    };
    &mut extra[index].1
}

/// Observes every key either side's Extra column holds, in code point order, a key missing
/// on one side comparing as `<absent>`.
fn compare_extra(
    reference_raw: &[u8],
    scored_raw: &[u8],
    extra: &mut Vec<(Vec<u8>, FieldStats)>,
    scratch: &mut ExtraScratch,
) {
    let ExtraScratch {
        reference,
        scored,
        keys,
        shape,
    } = scratch;
    parse_extra(reference_raw, reference);
    if reference_raw == scored_raw {
        for pair in reference.iter() {
            stats_for(extra, slice(reference_raw, pair.key)).observe_equal();
        }
        return;
    }
    parse_extra(scored_raw, scored);
    let key_of = |&(side, i): &(u8, usize)| -> &[u8] {
        if side == 0 {
            slice(reference_raw, reference[i].key)
        } else {
            slice(scored_raw, scored[i].key)
        }
    };
    keys.clear();
    keys.extend((0..reference.len()).map(|i| (0u8, i)));
    keys.extend((0..scored.len()).map(|i| (1u8, i)));
    keys.sort_unstable_by(|a, b| key_of(a).cmp(key_of(b)));
    keys.dedup_by(|a, b| key_of(a) == key_of(b));
    for entry in keys.iter() {
        let key = key_of(entry);
        let reference_value = reference
            .iter()
            .find(|p| slice(reference_raw, p.key) == key)
            .map_or(ABSENT, |p| slice(reference_raw, p.value));
        let scored_value = scored
            .iter()
            .find(|p| slice(scored_raw, p.key) == key)
            .map_or(ABSENT, |p| slice(scored_raw, p.value));
        stats_for(extra, key).observe(reference_value, scored_value, shape);
    }
}

/// Merges the two sorted files as the Python does: rows are read in file order, keys compare
/// bytewise, and a run of equal keys on one side counts its first row only.
pub fn merge(reference: &Path, scored: &Path) -> io::Result<MergeResult> {
    let mut reference_side = Side::open(reference)?;
    let mut scored_side = Side::open(scored)?;
    let mut totals = Totals::default();
    let mut fields: Vec<FieldStats> = (0..PLAIN_FIELDS.len())
        .map(|_| FieldStats::default())
        .collect();
    let mut extra: Vec<(Vec<u8>, FieldStats)> = Vec::new();
    let mut extra_scratch = ExtraScratch::default();
    let mut shape = Vec::new();

    let mut reference_live = reference_side.next_unique()?;
    let mut scored_live = scored_side.next_unique()?;
    while reference_live || scored_live {
        let order = if !scored_live {
            Ordering::Less
        } else if !reference_live {
            Ordering::Greater
        } else {
            reference_side.key().cmp(scored_side.key())
        };
        match order {
            Ordering::Less => {
                totals.reference_keys += 1;
                totals.reference_only += 1;
                reference_live = reference_side.next_unique()?;
            }
            Ordering::Greater => {
                totals.scored_keys += 1;
                totals.scored_only += 1;
                scored_live = scored_side.next_unique()?;
            }
            Ordering::Equal => {
                totals.reference_keys += 1;
                totals.scored_keys += 1;
                totals.intersection += 1;
                let reference_values = split_values(reference_side.values());
                let scored_values = split_values(scored_side.values());
                for (i, stats) in fields.iter_mut().enumerate() {
                    stats.observe(reference_values[i], scored_values[i], &mut shape);
                }
                compare_extra(
                    reference_values[VALUE_FIELDS - 1],
                    scored_values[VALUE_FIELDS - 1],
                    &mut extra,
                    &mut extra_scratch,
                );
                reference_live = reference_side.next_unique()?;
                scored_live = scored_side.next_unique()?;
            }
        }
    }
    totals.reference_duplicate_keys = reference_side.duplicate_keys;
    totals.scored_duplicate_keys = scored_side.duplicate_keys;
    extra.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(MergeResult {
        totals,
        fields,
        extra,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pairs(raw: &[u8]) -> Vec<(String, String)> {
        let mut out = Vec::new();
        parse_extra(raw, &mut out);
        out.iter()
            .map(|p| {
                (
                    String::from_utf8_lossy(slice(raw, p.key)).into_owned(),
                    String::from_utf8_lossy(slice(raw, p.value)).into_owned(),
                )
            })
            .collect()
    }

    fn owned(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect()
    }

    #[test]
    fn extra_parses_as_the_python_dict() {
        assert_eq!(pairs(b""), owned(&[]));
        assert_eq!(pairs(b"-"), owned(&[]));
        assert_eq!(
            pairs(b"IMPACT=HIGH;;CANONICAL;STRAND=1;STRAND=-1;=x;"),
            owned(&[
                ("IMPACT", "HIGH"),
                ("CANONICAL", ""),
                ("STRAND", "-1"),
                ("", "x")
            ])
        );
        assert_eq!(pairs(b"a=b=c"), owned(&[("a", "b=c")]));
    }

    #[test]
    fn top_mismatches_keep_first_observation_order_on_ties() {
        let mut stats = FieldStats::default();
        let mut scratch = Vec::new();
        for (r, s) in [("b", "x"), ("a", "x"), ("c", "x"), ("a", "x"), ("c", "x")] {
            stats.observe(r.as_bytes(), s.as_bytes(), &mut scratch);
        }
        stats.observe(b"same", b"same", &mut scratch);
        assert_eq!((stats.compared, stats.equal), (6, 1));
        assert_eq!(stats.distinct_mismatch_shapes(), 3);
        let order = |n: i64| -> Vec<(String, u64)> {
            stats
                .top_mismatches(n)
                .iter()
                .map(|s| (String::from_utf8_lossy(s.reference).into_owned(), s.count))
                .collect()
        };
        assert_eq!(order(1), vec![("a".to_string(), 2)]);
        assert_eq!(
            order(20),
            vec![
                ("a".to_string(), 2),
                ("c".to_string(), 2),
                ("b".to_string(), 1)
            ]
        );
        assert!(order(0).is_empty());
        assert!(order(-1).is_empty());
    }

    #[test]
    fn key_end_is_the_fifth_tab() {
        assert_eq!(key_end(b"a\tb\tc\td\te\tv1\tv2"), 9);
        assert_eq!(key_end(b"a\tb\tc"), 5);
        let values = split_values(b"1\t2\t3\t4\t5\t6\t7\t8\tX=1\tdropped");
        assert_eq!(values[0], b"1");
        assert_eq!(values[7], b"8");
        assert_eq!(values[8], b"X=1");
        let values = split_values(b"only");
        assert_eq!(values[0], b"only");
        assert_eq!(values[8], b"");
    }
}
