// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The walk over two sorted tuple files, and the forward-only reader of one sorted file.

use super::error::{Scored, ScorerError};
use super::interrupt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Take};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};

/// A run of whole lines of a file: the offset of its first byte and the offset after its
/// last.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: u64,
    pub end: u64,
}

impl Span {
    /// Every line of a file.
    pub const WHOLE: Span = Span {
        start: 0,
        end: u64::MAX,
    };
}

/// The bytes a reader of a span holds of its file at a time, at most.
const MOST_BUFFERED: u64 = 1 << 20;
const LEAST_BUFFERED: u64 = 1 << 13;

/// The lines of a span of a sorted file of the scratch, each with its line feed. A line
/// that sorts below the one before it in the file, the span's first line included, is a
/// refusal: the merge and the cursors read forward alone, so past such a line they would
/// miss the lines it belongs with and leave their values without a pair or a class. Lines
/// are compared without their line feed, as `sort` compares them under C collation.
pub struct Ascending {
    path: PathBuf,
    reader: BufReader<Take<File>>,
    start: u64,
    /// The lines read of the span.
    number: u64,
    line: Vec<u8>,
    previous: Vec<u8>,
}

fn body(line: &[u8]) -> &[u8] {
    line.strip_suffix(b"\n").unwrap_or(line)
}

/// The line of `file` that ends where the line at `start` begins, with its line feed.
fn line_before(file: &File, start: u64) -> std::io::Result<Vec<u8>> {
    let mut window = [0u8; 4096];
    let mut line: Vec<u8> = Vec::new();
    // The byte before `start` is the line's own line feed.
    let mut reached = start.saturating_sub(1);
    while reached > 0 {
        let from = reached.saturating_sub(window.len() as u64);
        let read = &mut window[..(reached - from) as usize];
        file.read_exact_at(read, from)?;
        match memchr::memrchr(b'\n', read) {
            Some(feed) => {
                line.splice(0..0, read[feed + 1..].iter().copied());
                break;
            }
            None => line.splice(0..0, read.iter().copied()),
        };
        reached = from;
    }
    if start > 0 {
        line.push(b'\n');
    }
    Ok(line)
}

/// The line feeds of `file` below `end`.
fn lines_before(file: &File, end: u64) -> std::io::Result<u64> {
    let mut window = vec![0u8; 1 << 16];
    let (mut reached, mut lines) = (0u64, 0u64);
    while reached < end {
        let wanted = (end - reached).min(window.len() as u64) as usize;
        let read = file.read_at(&mut window[..wanted], reached)?;
        if read == 0 {
            break;
        }
        lines += memchr::memchr_iter(b'\n', &window[..read]).count() as u64;
        reached += read as u64;
    }
    Ok(lines)
}

impl Ascending {
    /// Every line of the file.
    pub fn open(path: &Path) -> Scored<Ascending> {
        Ascending::part(path, Span::WHOLE)
    }

    /// The lines of `span`, which starts and ends where a line of the file does.
    pub fn part(path: &Path, span: Span) -> Scored<Ascending> {
        let failed = |error: std::io::Error| ScorerError::os(&error, path);
        let mut file = File::open(path).map_err(failed)?;
        let mut previous = Vec::new();
        if span.start > 0 {
            previous = line_before(&file, span.start).map_err(failed)?;
            file.seek(SeekFrom::Start(span.start)).map_err(failed)?;
        }
        let length = span.end.saturating_sub(span.start);
        let buffered = length.clamp(LEAST_BUFFERED, MOST_BUFFERED) as usize;
        Ok(Ascending {
            path: path.to_path_buf(),
            reader: BufReader::with_capacity(buffered, file.take(length)),
            start: span.start,
            number: 0,
            line: Vec::new(),
            previous,
        })
    }

    /// Reads the next line; it is empty past the span's last.
    pub fn advance(&mut self) -> Scored<()> {
        if self.number > 0 {
            std::mem::swap(&mut self.line, &mut self.previous);
        }
        self.line.clear();
        let read = self
            .reader
            .read_until(b'\n', &mut self.line)
            .map_err(|error| ScorerError::os(&error, &self.path))?;
        if read == 0 {
            return Ok(());
        }
        self.number += 1;
        if body(&self.line) < body(&self.previous) {
            return Err(self.out_of_order());
        }
        Ok(())
    }

    /// The refusal of the line in hand, which sorts below the one before it: it names the
    /// two lines by their numbers in the file.
    fn out_of_order(&self) -> ScorerError {
        let before = match self.start {
            0 => Ok(0),
            start => lines_before(self.reader.get_ref().get_ref(), start),
        };
        let before = match before {
            Ok(before) => before,
            Err(error) => return ScorerError::os(&error, &self.path),
        };
        let name = self.path.file_name().unwrap_or_default().to_string_lossy();
        let number = before + self.number;
        ScorerError(format!(
            "{name}: line {number} sorts below line {}: the file is not in ascending order",
            number - 1
        ))
    }

    /// The line read last, with its line feed; empty before the first and past the last.
    pub fn line(&self) -> &[u8] {
        &self.line
    }
}

/// What the walk over two tuple files hands on.
pub trait Settling {
    /// A line both files carry: a pair under its field name.
    fn paired(&mut self, field: &[u8]) -> Scored<()>;
    /// A run of lines sharing the columns before the value, `prefix` with its last tab,
    /// and the two sides' values under it.
    fn settle(&mut self, prefix: &[u8], ref_values: &[&[u8]], eng_values: &[&[u8]]) -> Scored<()>;
}

/// A line's columns before its value, with the tab that ends them.
fn prefix_of(line: &[u8], path: &Path) -> Scored<usize> {
    match memchr::memrchr(b'\t', line) {
        Some(tab) => Ok(tab + 1),
        None => Err(ScorerError(format!(
            "{}: a tuple line has no tab",
            path.file_name().unwrap_or_default().to_string_lossy()
        ))),
    }
}

/// Lines, or the values of lines, held back to back in one buffer.
#[derive(Debug, Default)]
pub struct Lines {
    bytes: Vec<u8>,
    ends: Vec<usize>,
}

impl Lines {
    pub fn clear(&mut self) {
        self.bytes.clear();
        self.ends.clear();
    }

    pub fn push(&mut self, line: &[u8]) {
        self.bytes.extend_from_slice(line);
        self.ends.push(self.bytes.len());
    }

    pub fn len(&self) -> usize {
        self.ends.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ends.is_empty()
    }

    /// The line `at`, counted from zero.
    pub fn get(&self, at: usize) -> &[u8] {
        let start = match at {
            0 => 0,
            _ => self.ends[at - 1],
        };
        &self.bytes[start..self.ends[at]]
    }

    pub fn iter(&self) -> impl Iterator<Item = &[u8]> {
        (0..self.ends.len()).map(|at| self.get(at))
    }

    fn slices(&self) -> Vec<&[u8]> {
        self.iter().collect()
    }

    /// Takes the lines of `file` that start with `prefix`, its line in hand first.
    fn take(&mut self, file: &mut Ascending, prefix: &[u8]) -> Scored<()> {
        while file.line().starts_with(prefix) {
            self.push(body(&file.line()[prefix.len()..]));
            file.advance()?;
        }
        Ok(())
    }
}

/// Walks two sorted tuple files together. A line both carry is a pair, counted under its
/// field name, the column before the value; any other run of lines sharing the columns
/// before the value goes to the settling step as that prefix and the two sides' values. A
/// file that is not in ascending order is a refusal. A stop signal the run has taken ends the
/// walk in a refusal before the next key.
pub fn merge(reference: &Path, engine: &Path, settling: &mut dyn Settling) -> Scored<()> {
    merge_spans((reference, Span::WHOLE), (engine, Span::WHOLE), settling)
}

/// As `merge`, over one span of each file.
pub fn merge_spans(
    (reference, reference_span): (&Path, Span),
    (engine, engine_span): (&Path, Span),
    settling: &mut dyn Settling,
) -> Scored<()> {
    let mut reference_lines = Ascending::part(reference, reference_span)?;
    let mut engine_lines = Ascending::part(engine, engine_span)?;
    reference_lines.advance()?;
    engine_lines.advance()?;
    let mut prefix: Vec<u8> = Vec::new();
    let (mut ref_values, mut eng_values) = (Lines::default(), Lines::default());
    loop {
        let (r, e) = (reference_lines.line(), engine_lines.line());
        if r.is_empty() && e.is_empty() {
            return Ok(());
        }
        if interrupt::requested() {
            return Err(ScorerError(
                "a stop signal ended the walk over the tuple files".to_string(),
            ));
        }
        if r == e {
            let head = &r[..prefix_of(r, reference)? - 1];
            let field = memchr::memrchr(b'\t', head).map_or(head, |tab| &head[tab + 1..]);
            settling.paired(field)?;
            reference_lines.advance()?;
            engine_lines.advance()?;
            continue;
        }
        let ref_prefix = match r.is_empty() {
            true => None,
            false => Some(&r[..prefix_of(r, reference)?]),
        };
        let eng_prefix = match e.is_empty() {
            true => None,
            false => Some(&e[..prefix_of(e, engine)?]),
        };
        let (takes_reference, takes_engine) = match (ref_prefix, eng_prefix) {
            (Some(one), Some(other)) => (one <= other, other <= one),
            (one, other) => (one.is_some(), other.is_some()),
        };
        prefix.clear();
        prefix.extend_from_slice(match takes_engine {
            true => eng_prefix.unwrap_or_default(),
            false => ref_prefix.unwrap_or_default(),
        });
        ref_values.clear();
        eng_values.clear();
        if takes_reference {
            ref_values.take(&mut reference_lines, &prefix)?;
        }
        if takes_engine {
            eng_values.take(&mut engine_lines, &prefix)?;
        }
        settling.settle(&prefix, &ref_values.slices(), &eng_values.slices())?;
    }
}

/// A forward-only reader of one sorted file; a file that is not in ascending order is a
/// refusal.
pub struct Cursor {
    lines: Ascending,
}

impl Cursor {
    pub fn open(path: &Path) -> Scored<Cursor> {
        Cursor::part(path, Span::WHOLE)
    }

    /// A reader of the lines of `span`, which starts and ends where a line of the file
    /// does.
    pub fn part(path: &Path, span: Span) -> Scored<Cursor> {
        let mut lines = Ascending::part(path, span)?;
        lines.advance()?;
        Ok(Cursor { lines })
    }

    /// The lines starting with `prefix`, a tab-terminated run of whole columns, without
    /// their line feed; successive calls take prefixes in ascending byte order.
    pub fn starting(&mut self, prefix: &[u8]) -> Scored<Vec<Vec<u8>>> {
        let mut found = Lines::default();
        self.starting_into(prefix, &mut found)?;
        Ok(found.iter().map(<[u8]>::to_vec).collect())
    }

    /// As `starting`, the lines appended to `found`.
    pub fn starting_into(&mut self, prefix: &[u8], found: &mut Lines) -> Scored<()> {
        while !self.lines.line().is_empty() && self.lines.line() < prefix {
            self.lines.advance()?;
        }
        while self.lines.line().starts_with(prefix) {
            found.push(body(self.lines.line()));
            self.lines.advance()?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn scratch(name: &str) -> PathBuf {
        let directory = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-merge-{name}-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&directory);
        fs::create_dir_all(&directory).unwrap();
        directory
    }

    fn file(directory: &Path, name: &str, lines: &[&str]) -> PathBuf {
        let path = directory.join(name);
        let text: String = lines.iter().map(|line| format!("{line}\n")).collect();
        fs::write(&path, text).unwrap();
        path
    }

    #[derive(Default)]
    struct Recorded {
        paired: Vec<String>,
        settled: Vec<(String, Vec<String>, Vec<String>)>,
    }

    impl Settling for Recorded {
        fn paired(&mut self, field: &[u8]) -> Scored<()> {
            self.paired
                .push(String::from_utf8_lossy(field).into_owned());
            Ok(())
        }

        fn settle(
            &mut self,
            prefix: &[u8],
            ref_values: &[&[u8]],
            eng_values: &[&[u8]],
        ) -> Scored<()> {
            let shown = |values: &[&[u8]]| {
                values
                    .iter()
                    .map(|value| String::from_utf8_lossy(value).into_owned())
                    .collect()
            };
            self.settled.push((
                String::from_utf8_lossy(prefix).into_owned(),
                shown(ref_values),
                shown(eng_values),
            ));
            Ok(())
        }
    }

    fn strings(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn runs_under_one_prefix_reach_the_settling_step_with_each_sides_values() {
        let directory = scratch("runs");
        let reference = file(
            &directory,
            "reference",
            &[
                "1:100\tT\tENST1\tPLANT_label\tsame",
                "1:100\tT\tENST1\tPLANT_score\t1",
                "1:100\tT\tENST1\tPLANT_score\t2",
                "1:200\tT\tENST1\tPLANT_score\t5",
                "1:400\tT\tENST1\tPLANT_score\t",
            ],
        );
        let engine = file(
            &directory,
            "engine",
            &[
                "1:100\tT\tENST1\tPLANT_label\tsame",
                "1:100\tT\tENST1\tPLANT_score\t2",
                "1:100\tT\tENST1\tPLANT_score\t3",
                "1:300\tT\tENST1\tPLANT_score\t7",
            ],
        );
        let mut recorded = Recorded::default();
        merge(&reference, &engine, &mut recorded).unwrap();
        assert_eq!(recorded.paired, ["PLANT_label"]);
        assert_eq!(
            recorded.settled,
            [
                (
                    "1:100\tT\tENST1\tPLANT_score\t".to_string(),
                    strings(&["1", "2"]),
                    strings(&["2", "3"])
                ),
                (
                    "1:200\tT\tENST1\tPLANT_score\t".to_string(),
                    strings(&["5"]),
                    strings(&[])
                ),
                (
                    "1:300\tT\tENST1\tPLANT_score\t".to_string(),
                    strings(&[]),
                    strings(&["7"])
                ),
                (
                    "1:400\tT\tENST1\tPLANT_score\t".to_string(),
                    strings(&[""]),
                    strings(&[])
                ),
            ]
        );
        let mut empty = Recorded::default();
        let nothing = file(&directory, "nothing", &[]);
        merge(&nothing, &nothing, &mut empty).unwrap();
        assert!(empty.paired.is_empty() && empty.settled.is_empty());
        fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn sorted_file_that_is_not_ascending_is_refused_by_the_merge_and_by_a_cursor() {
        let directory = scratch("order");
        let head = "1:100\tT\tPLANT_label\t";
        // A line that opens the next one sorts first whatever byte follows it there, and a
        // line may repeat.
        let sorted_lines = [
            format!("{head}a"),
            format!("{head}a\x01b"),
            format!("{head}a\x01b"),
            format!("{head}b"),
        ];
        let sorted_lines: Vec<&str> = sorted_lines.iter().map(String::as_str).collect();
        let in_order = file(&directory, "in-order", &sorted_lines);
        let mut lines = Ascending::open(&in_order).unwrap();
        for expected in &sorted_lines {
            lines.advance().unwrap();
            assert_eq!(lines.line(), format!("{expected}\n").as_bytes());
        }
        lines.advance().unwrap();
        assert!(lines.line().is_empty());
        let descending = ["1:200\tT\tPLANT_label\ta".to_string(), format!("{head}a")];
        let opened_line_last = [
            format!("{head}a"),
            format!("{head}a\x01b"),
            format!("{head}a"),
        ];
        for (name, lines, number) in [
            ("descending", descending.as_slice(), 2),
            ("opened-line-last", opened_line_last.as_slice(), 3),
        ] {
            let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
            let path = file(&directory, name, &lines);
            let refusal = ScorerError(format!(
                "{name}: line {number} sorts below line {}: the file is not in ascending order",
                number - 1
            ));
            let mut cursor = Cursor::open(&path).unwrap();
            assert_eq!(cursor.starting(b"9:1\t").err(), Some(refusal.clone()));
            for (reference, engine) in [(&path, &in_order), (&in_order, &path)] {
                let mut recorded = Recorded::default();
                assert_eq!(
                    merge(reference, engine, &mut recorded).err(),
                    Some(refusal.clone())
                );
            }
        }
        fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn cursor_reads_forward_the_lines_under_each_prefix() {
        let directory = scratch("cursor");
        let rows = file(
            &directory,
            "rows",
            &["7:100\tA\tx", "7:100\tT\ty", "7:1000\tT\tz", "7:200\tT\tw"],
        );
        let mut cursor = Cursor::open(&rows).unwrap();
        assert!(cursor.starting(b"7:10\t").unwrap().is_empty());
        assert_eq!(
            cursor.starting(b"7:100\t").unwrap(),
            [b"7:100\tA\tx".to_vec(), b"7:100\tT\ty".to_vec()]
        );
        assert!(cursor.starting(b"7:100\t").unwrap().is_empty());
        assert!(cursor.starting(b"7:150\t").unwrap().is_empty());
        assert_eq!(
            cursor.starting(b"7:200\t").unwrap(),
            [b"7:200\tT\tw".to_vec()]
        );
        assert!(cursor.starting(b"1:300\t").unwrap().is_empty());
        let missing = Cursor::open(&directory.join("absent")).err().unwrap();
        assert!(missing
            .0
            .starts_with("[Errno 2] No such file or directory: '"));
        fs::remove_dir_all(&directory).unwrap();
    }
}
