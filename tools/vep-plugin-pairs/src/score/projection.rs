// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two projections of one engine: the settling of each key's unpaired tuples, the
//! verdict lines the feature projection leaves for the allele projection, and the
//! book-keeping after a key is judged.

use super::class_table::Classes;
use super::context::RowContext;
use super::error::{Scored, ScorerError};
use super::interrupt;
use super::merge::{Cursor, Settling};
use super::number::pair_values;
use super::residuals::{ResidualLines, Side};
use super::sources::{text, Sources, RESOLVED_LIMIT};
use super::tally::Tally;
use std::cmp::Reverse;
use std::collections::{BTreeSet, BinaryHeap, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

/// The bytes of verdict lines held before they are sorted into a run of the scratch.
const RUN_BYTES: usize = 64 << 20;

/// The verdict lines of the feature projection, which reach their file sorted as `sort`
/// sorts lines under C collation, by the bytes of the line without its line feed, each
/// line once.
pub struct Verdicts {
    path: PathBuf,
    run_bytes: usize,
    bytes: Vec<u8>,
    starts: Vec<usize>,
    runs: Vec<PathBuf>,
}

impl Verdicts {
    pub fn new(path: &Path) -> Verdicts {
        Verdicts::with_runs_of(path, RUN_BYTES)
    }

    /// The same with a run written whenever `run_bytes` of lines are held.
    pub fn with_runs_of(path: &Path, run_bytes: usize) -> Verdicts {
        Verdicts {
            path: path.to_path_buf(),
            run_bytes,
            bytes: Vec::new(),
            starts: Vec::new(),
            runs: Vec::new(),
        }
    }

    /// The path of the next run: the verdicts file's with the run's number after it.
    fn run_path(&self) -> PathBuf {
        let mut name = self.path.clone().into_os_string();
        name.push(format!(".run{}", self.runs.len()));
        PathBuf::from(name)
    }

    /// Takes one line, its parts written back to back with a line feed after them.
    pub fn write(&mut self, parts: &[&[u8]]) -> Scored<()> {
        self.starts.push(self.bytes.len());
        for part in parts {
            self.bytes.extend_from_slice(part);
        }
        self.bytes.push(b'\n');
        if self.bytes.len() >= self.run_bytes {
            let run = self.run_path();
            self.runs.push(run.clone());
            self.write_sorted(&run)?;
        }
        Ok(())
    }

    /// Writes the lines held to `path` in order, each once, and holds none.
    fn write_sorted(&mut self, path: &Path) -> Scored<()> {
        stop_at_a_signal()?;
        let mut lines: Vec<&[u8]> = Vec::with_capacity(self.starts.len());
        for (at, &start) in self.starts.iter().enumerate() {
            let end = self.starts.get(at + 1).copied().unwrap_or(self.bytes.len());
            lines.push(&self.bytes[start..end - 1]);
        }
        lines.sort_unstable();
        lines.dedup();
        let failed = |error: std::io::Error| ScorerError::os(&error, path);
        let mut file = BufWriter::with_capacity(1 << 20, File::create(path).map_err(failed)?);
        for line in lines {
            file.write_all(line).map_err(failed)?;
            file.write_all(b"\n").map_err(failed)?;
        }
        file.flush().map_err(failed)?;
        self.bytes.clear();
        self.starts.clear();
        Ok(())
    }

    /// Leaves every line taken in the verdicts file, sorted, each once, and removes the
    /// runs.
    pub fn finish(mut self) -> Scored<()> {
        let path = self.path.clone();
        if self.runs.is_empty() {
            return self.write_sorted(&path);
        }
        if !self.starts.is_empty() {
            let run = self.run_path();
            self.runs.push(run.clone());
            self.write_sorted(&run)?;
        }
        let mut readers = Vec::new();
        let mut heap = BinaryHeap::new();
        for (at, run) in self.runs.iter().enumerate() {
            let file = File::open(run).map_err(|error| ScorerError::os(&error, run))?;
            readers.push(BufReader::with_capacity(1 << 20, file));
            let mut line = Vec::new();
            if next_line(&mut readers[at], &mut line, run)? {
                heap.push(Reverse((line, at)));
            }
        }
        let failed = |error: std::io::Error| ScorerError::os(&error, &path);
        let mut file = BufWriter::with_capacity(1 << 20, File::create(&path).map_err(failed)?);
        let mut written: Option<Vec<u8>> = None;
        while let Some(Reverse((line, at))) = heap.pop() {
            stop_at_a_signal()?;
            let mut next = Vec::new();
            if next_line(&mut readers[at], &mut next, &self.runs[at])? {
                heap.push(Reverse((next, at)));
            }
            if written.as_ref() != Some(&line) {
                file.write_all(&line).map_err(failed)?;
                file.write_all(b"\n").map_err(failed)?;
                written = Some(line);
            }
        }
        file.flush().map_err(failed)
    }
}

/// A refusal where the run has taken a stop signal, which ends the sorting of the verdict
/// lines before its next step.
fn stop_at_a_signal() -> Scored<()> {
    match interrupt::requested() {
        true => Err(ScorerError(
            "a stop signal ended the sorting of the verdict lines".to_string(),
        )),
        false => Ok(()),
    }
}

impl Drop for Verdicts {
    fn drop(&mut self) {
        for run in &self.runs {
            let _ = fs::remove_file(run);
        }
    }
}

/// Reads the next line of a run into `line` without its line feed; `false` past the last.
fn next_line(reader: &mut BufReader<File>, line: &mut Vec<u8>, path: &Path) -> Scored<bool> {
    let read = reader
        .read_until(b'\n', line)
        .map_err(|error| ScorerError::os(&error, path))?;
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    Ok(read > 0)
}

/// One projection of one engine over a run of keys: its tally, its residual lines where the
/// run keeps them, and the book-keeping after the unpaired tuples of one key are judged.
struct Books<'a> {
    tally: Tally,
    sources: &'a Sources<'a>,
    classes: &'a Classes,
    residuals: Option<ResidualLines>,
    named: HashMap<Vec<u8>, Rc<str>>,
}

impl<'a> Books<'a> {
    fn new(sources: &'a Sources<'a>, classes: &'a Classes, residuals: bool) -> Books<'a> {
        Books {
            tally: Tally::default(),
            sources,
            classes,
            residuals: residuals.then(ResidualLines::default),
            named: HashMap::new(),
        }
    }

    /// Pairs the two sides' values of one field under one key and counts the pairs and what
    /// each side has left, which is returned.
    fn pair<'v>(
        &mut self,
        name: &[u8],
        ref_values: &[&'v [u8]],
        eng_values: &[&'v [u8]],
    ) -> (Vec<&'v [u8]>, Vec<&'v [u8]>) {
        let (as_text, as_numbers, ref_left, eng_left) = pair_values(ref_values, eng_values);
        self.tally.count(
            name,
            (as_text + as_numbers) as u64,
            as_numbers as u64,
            ref_left.len() as u64,
            eng_left.len() as u64,
        );
        (ref_left, eng_left)
    }

    /// The source the field belongs to, looked up once per field name.
    fn source(&mut self, name: &[u8]) -> Scored<Rc<str>> {
        if let Some(source) = self.named.get(name) {
            return Ok(Rc::clone(source));
        }
        let source: Rc<str> = Rc::from(self.sources.named(name)?);
        if self.named.len() >= RESOLVED_LIMIT {
            self.named.clear();
        }
        self.named.insert(name.to_vec(), Rc::clone(&source));
        Ok(source)
    }

    fn applied(&self, class: &str) -> bool {
        self.classes.applied().contains(class)
    }

    fn set_aside(&mut self, side: Side, source: &str) {
        let aside = match side {
            Side::Reference => &mut self.tally.ref_aside,
            Side::Engine => &mut self.tally.eng_aside,
        };
        *aside.entry(source.to_string()).or_default() += 1;
    }

    /// Counts the key once under each class that accepted a tuple of it and keeps the
    /// `(value, side)` pairs left unexplained as residual lines.
    fn record(
        &mut self,
        source: &str,
        prefix: &[u8],
        accepted: &BTreeSet<String>,
        unexplained: &mut [(&[u8], Side)],
    ) {
        for class in accepted {
            let counts = if self.applied(class) {
                &mut self.tally.masked
            } else {
                &mut self.tally.not_applied
            };
            *counts
                .entry(source.to_string())
                .or_default()
                .entry(class.clone())
                .or_default() += 1;
        }
        if let Some(residuals) = &mut self.residuals {
            if !unexplained.is_empty() {
                residuals.take(source, prefix, unexplained);
            }
        }
    }
}

fn no_column(prefix: &[u8]) -> ScorerError {
    ScorerError(format!(
        "the tuple line under {} lacks a column",
        text(prefix)
    ))
}

/// The feature projection of one engine over a run of keys. Each tuple left unpaired is
/// written to the verdicts without its Feature, with its side and the class that accepted
/// it (empty for none), the lines the allele projection reads.
pub struct FeatureScore<'a> {
    books: Books<'a>,
    context: RowContext<'a>,
    verdicts: Verdicts,
}

impl<'a> FeatureScore<'a> {
    /// The projection, which keeps its residual lines where `residuals` is set.
    pub fn new(
        context: RowContext<'a>,
        verdicts: Verdicts,
        sources: &'a Sources<'a>,
        classes: &'a Classes,
        residuals: bool,
    ) -> FeatureScore<'a> {
        FeatureScore {
            books: Books::new(sources, classes, residuals),
            context,
            verdicts,
        }
    }

    /// The projection's tally, its residual lines and its verdict lines, which are not
    /// sorted yet.
    pub fn done(self) -> (Tally, Option<ResidualLines>, Verdicts) {
        (self.books.tally, self.books.residuals, self.verdicts)
    }
}

impl Settling for FeatureScore<'_> {
    fn paired(&mut self, field: &[u8]) -> Scored<()> {
        *self.books.tally.paired.entry(field.to_vec()).or_default() += 1;
        Ok(())
    }

    fn settle(&mut self, prefix: &[u8], ref_values: &[&[u8]], eng_values: &[&[u8]]) -> Scored<()> {
        let body = &prefix[..prefix.len() - 1];
        let cut = memchr::memrchr(b'\t', body).ok_or_else(|| no_column(prefix))? + 1;
        let (key, name) = body.split_at(cut);
        let (ref_left, eng_left) = self.books.pair(name, ref_values, eng_values);
        let allele_cut =
            memchr::memrchr(b'\t', &key[..key.len() - 1]).ok_or_else(|| no_column(prefix))?;
        let allele_key = &key[..=allele_cut];
        let source = self.books.source(name)?;
        let mut accepted = BTreeSet::new();
        let mut unexplained: Vec<(&[u8], Side)> = Vec::new();
        for (ref_value, eng_value, ref_nth, eng_nth) in
            self.context.pairs(key, name, &ref_left, &eng_left)?
        {
            let class = self
                .context
                .judge(key, name, ref_value, eng_value, ref_nth, eng_nth)?;
            let tag = class.as_deref().unwrap_or("");
            let applied = class
                .as_deref()
                .is_some_and(|class| self.books.applied(class));
            for (value, side) in [(ref_value, Side::Reference), (eng_value, Side::Engine)] {
                let Some(value) = value else {
                    continue;
                };
                self.verdicts.write(&[
                    allele_key,
                    name,
                    b"\t",
                    value,
                    b"\t",
                    &[side.letter()],
                    b"\t",
                    tag.as_bytes(),
                ])?;
                if applied {
                    self.books.set_aside(side, &source);
                } else {
                    unexplained.push((value, side));
                }
            }
            accepted.extend(class);
        }
        self.books
            .record(&source, prefix, &accepted, &mut unexplained);
        Ok(())
    }
}

/// The allele projection of one engine over a run of keys, the classes of its unpaired
/// tuples read from the feature projection's sorted verdict lines of the same keys.
pub struct AlleleScore<'a> {
    books: Books<'a>,
    reference_verdicts: Cursor,
    engine_verdicts: Cursor,
}

impl<'a> AlleleScore<'a> {
    /// The projection, which keeps its residual lines where `residuals` is set.
    pub fn new(
        verdicts: &Path,
        sources: &'a Sources<'a>,
        classes: &'a Classes,
        residuals: bool,
    ) -> Scored<AlleleScore<'a>> {
        Ok(AlleleScore {
            books: Books::new(sources, classes, residuals),
            reference_verdicts: Cursor::open(verdicts)?,
            engine_verdicts: Cursor::open(verdicts)?,
        })
    }

    /// The projection's tally and its residual lines.
    pub fn done(self) -> (Tally, Option<ResidualLines>) {
        (self.books.tally, self.books.residuals)
    }

    /// Of `tags`, the class that comes first in the class table.
    fn first_in_table<'t>(&self, tags: impl Iterator<Item = &'t String>) -> Scored<String> {
        let mut first: Option<(usize, &String)> = None;
        for tag in tags {
            let order = self.books.classes.order(tag).ok_or_else(|| {
                ScorerError(format!(
                    "a verdict names {tag}, a class the class table does not declare"
                ))
            })?;
            if first.is_none_or(|(least, _)| order < least) {
                first = Some((order, tag));
            }
        }
        first
            .map(|(_, tag)| tag.clone())
            .ok_or_else(|| ScorerError("no class among the verdicts of a tuple".to_string()))
    }
}

impl Settling for AlleleScore<'_> {
    fn paired(&mut self, field: &[u8]) -> Scored<()> {
        *self.books.tally.paired.entry(field.to_vec()).or_default() += 1;
        Ok(())
    }

    fn settle(&mut self, prefix: &[u8], ref_values: &[&[u8]], eng_values: &[&[u8]]) -> Scored<()> {
        let body = &prefix[..prefix.len() - 1];
        let cut = memchr::memrchr(b'\t', body).ok_or_else(|| no_column(prefix))? + 1;
        let name = &body[cut..];
        let (ref_left, eng_left) = self.books.pair(name, ref_values, eng_values);
        let source = self.books.source(name)?;
        let mut accepted = BTreeSet::new();
        let mut unexplained: Vec<(&[u8], Side)> = Vec::new();
        for (values, side) in [(ref_left, Side::Reference), (eng_left, Side::Engine)] {
            let mut lookups: Vec<(Vec<u8>, &[u8])> = values
                .into_iter()
                .map(|value| {
                    let lookup = [prefix, value, b"\t", &[side.letter()], b"\t"].concat();
                    (lookup, value)
                })
                .collect();
            lookups.sort();
            for (lookup, value) in lookups {
                let cursor = match side {
                    Side::Reference => &mut self.reference_verdicts,
                    Side::Engine => &mut self.engine_verdicts,
                };
                let tags: Vec<String> = cursor
                    .starting(&lookup)?
                    .iter()
                    .map(|line| text(&line[lookup.len()..]))
                    .collect();
                // The verdict without a class sorts first, so a class there means no
                // feature tuple of this allele tuple was left unexplained.
                let classed = tags.first().is_some_and(|tag| !tag.is_empty());
                if classed && tags.iter().all(|tag| self.books.applied(tag)) {
                    self.books.set_aside(side, &source);
                    accepted.insert(self.first_in_table(tags.iter())?);
                } else {
                    unexplained.push((value, side));
                    if classed {
                        let held = tags.iter().filter(|tag| !self.books.applied(tag));
                        accepted.insert(self.first_in_table(held)?);
                    }
                }
            }
        }
        self.books
            .record(&source, prefix, &accepted, &mut unexplained);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let directory = std::env::temp_dir().join(format!(
            "vep-plugin-pairs-projection-{name}-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&directory);
        fs::create_dir_all(&directory).unwrap();
        directory
    }

    const LINES: [&[u8]; 9] = [
        b"1:100\tT\tPLANT_score\tb\tE\t",
        b"1:100\tT\tPLANT_score\ta\tR\tplanted",
        b"1:100\tT\tPLANT_score\ta\tR\t",
        b"1:100\tT\tPLANT_score\ta\x01\tR\t",
        b"1:100\tT\tPLANT_score\ta\tR\tplanted",
        b"1:100\tT\tPLANT_score\ta\tE\t",
        b"1:100\tT\tPLANT_score\t\xff\tE\t",
        b"1:100\tT\tPLANT_score\ta\tR\t",
        b"1:100\tT\tPLANT_label\tz\tR\t",
    ];
    const SORTED: &[u8] = b"1:100\tT\tPLANT_label\tz\tR\t\n\
        1:100\tT\tPLANT_score\ta\x01\tR\t\n\
        1:100\tT\tPLANT_score\ta\tE\t\n\
        1:100\tT\tPLANT_score\ta\tR\t\n\
        1:100\tT\tPLANT_score\ta\tR\tplanted\n\
        1:100\tT\tPLANT_score\tb\tE\t\n\
        1:100\tT\tPLANT_score\t\xff\tE\t\n";

    #[test]
    fn verdict_lines_reach_their_file_in_byte_order_each_once() {
        let directory = scratch("held");
        let path = directory.join("input1.verdicts");
        let mut verdicts = Verdicts::new(&path);
        for line in LINES {
            verdicts.write(&[&line[..5], &line[5..]]).unwrap();
        }
        verdicts.finish().unwrap();
        assert!(fs::read(&path).unwrap() == SORTED);
        let none = directory.join("input2.verdicts");
        Verdicts::new(&none).finish().unwrap();
        assert_eq!(fs::read(&none).unwrap(), b"");
        fs::remove_dir_all(&directory).unwrap();
    }

    #[test]
    fn verdict_lines_sorted_in_runs_merge_into_the_same_file() {
        let directory = scratch("runs");
        for run_bytes in [1, 40, 70, 150] {
            let path = directory.join("input1.verdicts");
            let mut verdicts = Verdicts::with_runs_of(&path, run_bytes);
            for line in LINES {
                verdicts.write(&[line]).unwrap();
            }
            verdicts.finish().unwrap();
            assert_eq!(
                fs::read(&path).unwrap(),
                SORTED,
                "runs of {run_bytes} bytes"
            );
            let left: Vec<_> = fs::read_dir(&directory)
                .unwrap()
                .map(|entry| entry.unwrap().file_name())
                .collect();
            assert_eq!(left, ["input1.verdicts"], "runs of {run_bytes} bytes");
        }
        fs::remove_dir_all(&directory).unwrap();
    }
}
