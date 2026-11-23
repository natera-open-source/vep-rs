// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The streaming merge of two sorted-unique key streams: the tuple counts, the one-sided rows
//! for `discordant.tsv` (the reference's rows first, then vep-rs's) and the per-set counters,
//! with the selected release's one-sided collector observing every key. The merge holds each
//! key against the one before it on its side and stops at the first that does not sort after
//! it.

use crate::classes::Counter;
use crate::rules::scope::{ReferenceRecordScope, ScopeFeatures};
use crate::sort::SortedKeys;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// A side of the comparison, named in a message by the option that gives its directory.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Reference,
    Scored,
}

impl Side {
    fn option(self) -> &'static str {
        match self {
            Side::Reference => "--perl-dir",
            Side::Scored => "--rust-dir",
        }
    }
}

/// The refusal of a key that does not sort after the key before it on its side. The merge
/// steps forward on the comparison of the two current keys alone, so past such a key it
/// would count a key both sides hold as missing and as extra at once. The message names the
/// side and, once [`written_at`] has placed the two keys, the output file and the line of it
/// that first writes each. It never names a key, which holds a row's Location, allele,
/// feature and consequence terms, and it names no partition and no place among the keys a
/// merge read, so one fault reads the same at any thread count.
#[derive(Debug)]
pub struct KeyOrderError {
    side: Side,
    collapsed: bool,
    /// The key that does not sort after `before`, the key before it on its side.
    key: Vec<u8>,
    before: Vec<u8>,
    written_at: Option<WrittenAt>,
}

/// Where a side's output file writes the two keys of a [`KeyOrderError`]: the file, and the
/// lines of it that first write `key` and `before`.
#[derive(Debug)]
struct WrittenAt {
    output: PathBuf,
    key_line: u64,
    before_line: u64,
}

impl KeyOrderError {
    fn of(side: Side, key: &[u8], before: &[u8]) -> KeyOrderError {
        KeyOrderError {
            side,
            collapsed: false,
            key: key.to_vec(),
            before: before.to_vec(),
            written_at: None,
        }
    }

    pub fn side(&self) -> Side {
        self.side
    }

    /// Whether the keys are context-collapsed ones.
    pub fn collapsed(&self) -> bool {
        self.collapsed
    }

    /// The key that does not sort after the key before it.
    pub fn key(&self) -> &[u8] {
        &self.key
    }

    /// The key before [`KeyOrderError::key`] on its side.
    pub fn before(&self) -> &[u8] {
        &self.before
    }
}

impl fmt::Display for KeyOrderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let keys = match self.collapsed {
            true => "context-collapsed keys",
            false => "keys",
        };
        let side = self.side.option();
        write!(f, "the sorted {keys} of the {side} side do not ascend: ")?;
        match &self.written_at {
            Some(at) => write!(
                f,
                "the key of line {} of {} does not sort after the key of line {}, which stands \
                 before it",
                at.key_line,
                at.output.display(),
                at.before_line
            )?,
            None => f.write_str("a key of the side does not sort after the key before it")?,
        }
        f.write_str("; the comparison needs each side's keys in order with every key once")
    }
}

impl std::error::Error for KeyOrderError {}

impl From<KeyOrderError> for io::Error {
    fn from(error: KeyOrderError) -> io::Error {
        io::Error::other(error)
    }
}

fn key_order_of(error: &mut io::Error) -> Option<&mut KeyOrderError> {
    error.get_mut()?.downcast_mut()
}

/// The refusal of keys that do not ascend that `error` carries, None for any other error.
pub fn key_order(error: &io::Error) -> Option<&KeyOrderError> {
    error.get_ref()?.downcast_ref()
}

/// Whether `error` is the refusal of keys that do not ascend.
pub fn is_key_order(error: &io::Error) -> bool {
    key_order(error).is_some()
}

/// `error` with the output file of its side and the lines of that file that first write its
/// key and the key before it named as the places of the two keys, where it is the refusal of
/// keys that do not ascend.
pub fn written_at(
    mut error: io::Error,
    output: &Path,
    key_line: u64,
    before_line: u64,
) -> io::Error {
    if let Some(order) = key_order_of(&mut error) {
        order.written_at = Some(WrittenAt {
            output: output.to_path_buf(),
            key_line,
            before_line,
        });
    }
    error
}

/// `error` with the context-collapsed keys named as the keys of the merge, where it is the
/// refusal of keys that do not ascend.
pub fn in_collapsed_keys(mut error: io::Error) -> io::Error {
    if let Some(order) = key_order_of(&mut error) {
        order.collapsed = true;
    }
    error
}

/// One side's keys from each partition to the next: within a partition the merge holds
/// every key against the one before it, and this holds the first key of each partition
/// against the last key of the partition before it, so that the side's keys ascend over
/// the whole of it in the order the discordant rows are written.
pub struct PartitionOrder {
    side: Side,
    last: Option<(usize, Vec<u8>)>,
    refused: Option<Vec<u8>>,
}

impl PartitionOrder {
    pub fn new(side: Side) -> PartitionOrder {
        PartitionOrder {
            side,
            last: None,
            refused: None,
        }
    }

    /// After a refusal by [`PartitionOrder::follow`]: the partition before the one refused,
    /// its last key, and the first key that does not sort after it.
    pub fn refused(&self) -> Option<(usize, &[u8], &[u8])> {
        let (before, last) = self.last.as_ref()?;
        Some((*before, last, self.refused.as_deref()?))
    }

    /// Takes the first and last keys of `partition`, the partitions given in their order;
    /// a partition without a key gives None.
    pub fn follow(&mut self, partition: usize, ends: Option<(Vec<u8>, Vec<u8>)>) -> io::Result<()> {
        let Some((first, last)) = ends else {
            return Ok(());
        };
        if let Some((_, previous)) = &self.last {
            if first <= *previous {
                let refusal = KeyOrderError::of(self.side, &first, previous);
                self.refused = Some(first);
                return Err(refusal.into());
            }
        }
        self.last = Some((partition, last));
        Ok(())
    }
}

/// The one-sided collector the merge feeds, if any.
pub enum Collector<'a> {
    None,
    Release115(&'a mut ScopeFeatures),
    Release116(&'a mut ReferenceRecordScope),
}

#[derive(Default)]
pub struct MergeCounts {
    pub perl_count: u64,
    pub rust_count: u64,
    pub intersection: u64,
    pub missing: u64,
    pub extra: u64,
    pub missing_counter: Counter,
    pub extra_counter: Counter,
}

impl MergeCounts {
    /// Adds the counts of `next`, the merge of the key range that follows this one.
    pub fn append(&mut self, next: &MergeCounts) {
        self.perl_count += next.perl_count;
        self.rust_count += next.rust_count;
        self.intersection += next.intersection;
        self.missing += next.missing;
        self.extra += next.extra;
        self.missing_counter.update(&next.missing_counter);
        self.extra_counter.update(&next.extra_counter);
    }
}

/// Sorted unique keys read one at a time.
pub trait KeyStream {
    /// The key the stream stands at, None past its last.
    fn current(&self) -> Option<&[u8]>;
    /// The key before the one the stream stands at, None at its first key and past its last.
    fn previous(&self) -> Option<&[u8]>;
    fn advance(&mut self) -> io::Result<()>;
}

struct Lines {
    reader: BufReader<File>,
    line: Vec<u8>,
    before: Vec<u8>,
    read: u64,
    done: bool,
}

impl Lines {
    fn open(path: &Path) -> io::Result<Lines> {
        let mut lines = Lines {
            reader: BufReader::with_capacity(1 << 20, File::open(path)?),
            line: Vec::new(),
            before: Vec::new(),
            read: 0,
            done: false,
        };
        lines.advance()?;
        Ok(lines)
    }
}

impl KeyStream for Lines {
    fn current(&self) -> Option<&[u8]> {
        (!self.done).then_some(self.line.as_slice())
    }

    fn previous(&self) -> Option<&[u8]> {
        (!self.done && self.read > 1).then_some(self.before.as_slice())
    }

    fn advance(&mut self) -> io::Result<()> {
        std::mem::swap(&mut self.line, &mut self.before);
        self.line.clear();
        if self.reader.read_until(b'\n', &mut self.line)? == 0 {
            self.done = true;
            return Ok(());
        }
        if self.line.last() == Some(&b'\n') {
            self.line.pop();
        }
        self.read += 1;
        Ok(())
    }
}

/// The keys of a [`SortedKeys`] as a stream.
pub struct MemoryKeys<'a> {
    keys: &'a SortedKeys,
    at: usize,
}

impl<'a> MemoryKeys<'a> {
    pub fn new(keys: &'a SortedKeys) -> MemoryKeys<'a> {
        MemoryKeys { keys, at: 0 }
    }
}

impl KeyStream for MemoryKeys<'_> {
    fn current(&self) -> Option<&[u8]> {
        (self.at < self.keys.len()).then(|| self.keys.get(self.at))
    }

    fn previous(&self) -> Option<&[u8]> {
        (0 < self.at && self.at < self.keys.len()).then(|| self.keys.get(self.at - 1))
    }

    fn advance(&mut self) -> io::Result<()> {
        self.at += 1;
        Ok(())
    }
}

/// Steps `stream` to its next key on `side`, which must sort after the key before it.
fn advance_ascending(stream: &mut impl KeyStream, side: Side) -> io::Result<()> {
    stream.advance()?;
    match (stream.previous(), stream.current()) {
        (Some(previous), Some(key)) if key <= previous => {
            Err(KeyOrderError::of(side, key, previous).into())
        }
        _ => Ok(()),
    }
}

fn consequence_set(key: &[u8]) -> &str {
    let start = memchr::memrchr(b'\t', key).map_or(0, |i| i + 1);
    std::str::from_utf8(&key[start..]).expect("utf-8 key")
}

fn observe(collector: &mut Collector, file_name: &str, key: &[u8], perl: bool, rust: bool) {
    match collector {
        Collector::None => {}
        Collector::Release115(scope) => {
            let key = std::str::from_utf8(key).expect("utf-8 key");
            if perl {
                scope.record(key, false);
            }
            if rust {
                scope.record(key, true);
            }
        }
        Collector::Release116(scope) => {
            scope.observe(
                file_name,
                std::str::from_utf8(key).expect("utf-8 key"),
                perl,
                rust,
            );
        }
    }
}

/// Merges the two files. With `discordant_out`, the reference-only rows are written as they
/// are met and the vep-rs-only rows after them; with `track_buckets`, each side's one-sided
/// consequence sets are counted. A 116.2 collector's open Location group is closed when both
/// streams are exhausted.
pub fn compare_sorted_key_files(
    perl_sorted: &Path,
    rust_sorted: &Path,
    file_name: &str,
    discordant_out: Option<&mut dyn Write>,
    track_buckets: bool,
    collector: Collector,
) -> io::Result<MergeCounts> {
    let mut perl = Lines::open(perl_sorted)?;
    let mut rust = Lines::open(rust_sorted)?;
    let Some(out) = discordant_out else {
        return merge_keys(
            &mut perl,
            &mut rust,
            file_name,
            None::<&mut Vec<u8>>,
            None,
            track_buckets,
            collector,
        );
    };
    let mut extra_rows: Vec<u8> = Vec::new();
    let counts = merge_keys(
        &mut perl,
        &mut rust,
        file_name,
        Some(&mut *out),
        Some(&mut extra_rows),
        track_buckets,
        collector,
    )?;
    out.write_all(&extra_rows)?;
    Ok(counts)
}

/// Merges two key streams. The reference-only rows go to `missing_out` as they are met and
/// the vep-rs-only rows are appended to `extra_rows`, each as a `discordant.tsv` row; with
/// `track_buckets`, each side's one-sided consequence sets are counted. A 116.2 collector's
/// open Location group is closed when both streams are exhausted. A key that does not sort
/// after the key before it on its side ends the merge with a [`KeyOrderError`].
pub fn merge_keys<W: Write + ?Sized>(
    perl: &mut impl KeyStream,
    rust: &mut impl KeyStream,
    file_name: &str,
    mut missing_out: Option<&mut W>,
    mut extra_rows: Option<&mut Vec<u8>>,
    track_buckets: bool,
    mut collector: Collector,
) -> io::Result<MergeCounts> {
    let mut counts = MergeCounts::default();
    loop {
        let (p, r) = (perl.current(), rust.current());
        match (p, r) {
            (None, None) => break,
            (Some(p), r) if r.is_none_or(|r| p < r) => {
                counts.perl_count += 1;
                counts.missing += 1;
                observe(&mut collector, file_name, p, true, false);
                if track_buckets {
                    counts.missing_counter.add(consequence_set(p), 1);
                }
                if let Some(out) = missing_out.as_deref_mut() {
                    out.write_all(file_name.as_bytes())?;
                    out.write_all(b"\tmissing_in_rust\t")?;
                    out.write_all(p)?;
                    out.write_all(b"\n")?;
                }
                advance_ascending(perl, Side::Reference)?;
            }
            (p, Some(r)) if p.is_none_or(|p| r < p) => {
                counts.rust_count += 1;
                counts.extra += 1;
                observe(&mut collector, file_name, r, false, true);
                if track_buckets {
                    counts.extra_counter.add(consequence_set(r), 1);
                }
                if let Some(extra_rows) = extra_rows.as_deref_mut() {
                    extra_rows.extend_from_slice(file_name.as_bytes());
                    extra_rows.extend_from_slice(b"\textra_in_rust\t");
                    extra_rows.extend_from_slice(r);
                    extra_rows.push(b'\n');
                }
                advance_ascending(rust, Side::Scored)?;
            }
            (Some(p), Some(_)) => {
                counts.perl_count += 1;
                counts.rust_count += 1;
                counts.intersection += 1;
                observe(&mut collector, file_name, p, true, true);
                advance_ascending(perl, Side::Reference)?;
                advance_ascending(rust, Side::Scored)?;
            }
            _ => unreachable!("one side is present"),
        }
    }
    if let Collector::Release116(scope) = &mut collector {
        scope.close();
    }
    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TempDir;

    /// Four keys in ascending order.
    const KEYS: [&str; 4] = [
        "1:100\tk\tt1\tTranscript\tintron_variant",
        "1:100\tk\tt2\tTranscript\tintron_variant",
        "1:250\t-\tt1\tTranscript\tupstream_gene_variant",
        "2:7\tdeletion\tt3\tTranscript\tintron_variant",
    ];

    const NEEDS: &str = "; the comparison needs each side's keys in order with every key once";

    fn lines(keys: &[&str]) -> Vec<u8> {
        keys.iter()
            .flat_map(|key| key.bytes().chain(*b"\n"))
            .collect()
    }

    /// The merge of two key files that hold `reference` and `scored` in the order given, with
    /// the discordant rows it wrote.
    fn merge_of_files(reference: &[&str], scored: &[&str]) -> io::Result<(MergeCounts, Vec<u8>)> {
        let scratch = TempDir::create(&std::env::temp_dir(), "merge_test_")?;
        let (perl, rust) = (scratch.path.join("perl"), scratch.path.join("rust"));
        std::fs::write(&perl, lines(reference))?;
        std::fs::write(&rust, lines(scored))?;
        let mut rows = Vec::new();
        let counts = compare_sorted_key_files(
            &perl,
            &rust,
            "pair.txt",
            Some(&mut rows),
            true,
            Collector::None,
        )?;
        Ok((counts, rows))
    }

    /// The merge of the same keys held in memory, as a partition's are.
    fn merge_in_memory(reference: &[&str], scored: &[&str]) -> io::Result<(MergeCounts, Vec<u8>)> {
        let perl = SortedKeys::in_the_order_written(lines(reference));
        let rust = SortedKeys::in_the_order_written(lines(scored));
        let (mut rows, mut extra_rows) = (Vec::new(), Vec::new());
        let counts = merge_keys(
            &mut MemoryKeys::new(&perl),
            &mut MemoryKeys::new(&rust),
            "pair.txt",
            Some(&mut rows),
            Some(&mut extra_rows),
            true,
            Collector::None,
        )?;
        rows.extend(extra_rows);
        Ok((counts, rows))
    }

    /// What both merges stop with, the same over key files and over keys in memory: the
    /// message, the key the stop is on and the key before it.
    fn stop_of(reference: &[&str], scored: &[&str]) -> (String, String, String) {
        let stops = [
            merge_of_files(reference, scored),
            merge_in_memory(reference, scored),
        ]
        .map(|merged| match merged {
            Err(error) => {
                let order = key_order(&error).expect("the refusal of keys that do not ascend");
                let text = |key: &[u8]| String::from_utf8(key.to_vec()).unwrap();
                (error.to_string(), text(order.key()), text(order.before()))
            }
            Ok((counts, _)) => panic!(
                "the merge went on and counted {} keys on both sides, {} missing and {} extra",
                counts.intersection, counts.missing, counts.extra
            ),
        });
        assert_eq!(stops[0], stops[1], "key files against keys in memory");
        stops[0].clone()
    }

    /// A stop on `key` behind the key `before` on the side `option` names, as [`stop_of`]
    /// gives it, with the message of a stop whose keys no line of an output file is named
    /// for: it names the side and neither key.
    fn stop_on(option: &str, key: &str, before: &str) -> (String, String, String) {
        let message = format!(
            "the sorted keys of the {option} side do not ascend: a key of the side does not \
             sort after the key before it{NEEDS}"
        );
        (message, key.to_string(), before.to_string())
    }

    #[test]
    fn keys_that_ascend_are_counted_and_written() {
        let reference = [KEYS[0], KEYS[1], KEYS[3]];
        let scored = [KEYS[1], KEYS[2], KEYS[3]];
        for merged in [
            merge_of_files(&reference, &scored),
            merge_in_memory(&reference, &scored),
        ] {
            let (counts, rows) = merged.unwrap();
            assert_eq!(
                [counts.perl_count, counts.rust_count, counts.intersection],
                [3, 3, 2]
            );
            assert_eq!([counts.missing, counts.extra], [1, 1]);
            assert_eq!(
                String::from_utf8(rows).unwrap(),
                format!(
                    "pair.txt\tmissing_in_rust\t{}\npair.txt\textra_in_rust\t{}\n",
                    KEYS[0], KEYS[2]
                )
            );
        }
    }

    #[test]
    fn a_reference_key_below_the_key_before_it_stops_the_merge() {
        assert_eq!(
            stop_of(&[KEYS[0], KEYS[2], KEYS[1], KEYS[3]], &KEYS),
            stop_on("--perl-dir", KEYS[1], KEYS[2])
        );
    }

    #[test]
    fn a_scored_key_below_the_key_before_it_stops_the_merge() {
        assert_eq!(
            stop_of(&KEYS, &[KEYS[0], KEYS[1], KEYS[3], KEYS[2]]),
            stop_on("--rust-dir", KEYS[2], KEYS[3])
        );
    }

    #[test]
    fn a_key_written_twice_stops_the_merge() {
        assert_eq!(
            stop_of(&[KEYS[0], KEYS[1], KEYS[1], KEYS[2]], &KEYS[..3]),
            stop_on("--perl-dir", KEYS[1], KEYS[1])
        );
        assert_eq!(
            stop_of(&KEYS[..2], &[KEYS[0], KEYS[0], KEYS[1]]),
            stop_on("--rust-dir", KEYS[0], KEYS[0])
        );
    }

    #[test]
    fn the_message_names_the_collapsed_keys_and_the_lines_that_write_the_two_keys_and_no_key() {
        let stop = || {
            merge_in_memory(&[KEYS[1], KEYS[0]], &[])
                .map(drop)
                .unwrap_err()
        };
        let refusal = stop();
        let order = key_order(&refusal).expect("the refusal of keys that do not ascend");
        assert_eq!(
            (order.side(), order.collapsed(), order.key(), order.before()),
            (
                Side::Reference,
                false,
                KEYS[0].as_bytes(),
                KEYS[1].as_bytes()
            )
        );
        let output = Path::new("reference/pair.txt");
        let placed = |keys: &str| {
            format!(
                "the sorted {keys} of the --perl-dir side do not ascend: the key of line 7 of \
                 reference/pair.txt does not sort after the key of line 3, which stands before \
                 it{NEEDS}"
            )
        };
        assert_eq!(written_at(stop(), output, 7, 3).to_string(), placed("keys"));
        let collapsed = written_at(in_collapsed_keys(stop()), output, 7, 3);
        assert!(key_order(&collapsed).is_some_and(KeyOrderError::collapsed));
        assert_eq!(collapsed.to_string(), placed("context-collapsed keys"));
        let other = written_at(
            in_collapsed_keys(io::Error::other("unreadable")),
            output,
            7,
            3,
        );
        assert!(!is_key_order(&other));
        assert_eq!(other.to_string(), "unreadable");
    }

    #[test]
    fn the_first_key_of_a_partition_sorts_after_the_last_key_before_it() {
        let ends = |first: usize, last: usize| {
            Some((
                KEYS[first].as_bytes().to_vec(),
                KEYS[last].as_bytes().to_vec(),
            ))
        };
        let mut order = PartitionOrder::new(Side::Reference);
        order.follow(0, ends(0, 1)).unwrap();
        order.follow(1, None).unwrap();
        order.follow(2, ends(2, 2)).unwrap();
        order.follow(3, ends(3, 3)).unwrap();
        for (first, side, option) in [
            (2, Side::Reference, "--perl-dir"),
            (1, Side::Scored, "--rust-dir"),
        ] {
            let mut order = PartitionOrder::new(side);
            order.follow(1, ends(0, 2)).unwrap();
            order.follow(2, None).unwrap();
            let stop = order.follow(4, ends(first, 3)).unwrap_err();
            let between = key_order(&stop).expect("the refusal of keys that do not ascend");
            assert_eq!(
                (between.key(), between.before()),
                (KEYS[first].as_bytes(), KEYS[2].as_bytes())
            );
            assert_eq!(
                stop.to_string(),
                stop_on(option, KEYS[first], KEYS[2]).0,
                "a stop between two partitions reads as a stop inside one"
            );
        }
    }
}

#[cfg(test)]
#[path = "merge_key_order_tests.rs"]
mod key_order_tests;

#[cfg(test)]
#[path = "merge_collectors_tests.rs"]
mod collectors_tests;
