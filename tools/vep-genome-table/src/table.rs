// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The four files of a genome table directory: reading them, checking them against each
//! other and rewriting the two row files in key order.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt;
use std::fs::{self, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process;

use flate2::{Compression, GzBuilder};

use crate::reader;
use crate::row::{self, OutcomeKind, Row};

/// The plain row file: every row whose set is not bacterial.
pub const GENOMES_FILE: &str = "genomes.tsv";
/// The compressed row file: every row of a `bacteria_<N>_collection` set.
pub const BACTERIA_FILE: &str = "bacteria.tsv.gz";
/// The declarations of every `differs:` and `excluded:` value a row uses.
pub const OUTCOMES_FILE: &str = "outcomes.tsv";
/// The sets and how many ready-made caches each holds.
pub const SETS_FILE: &str = "sets.tsv";

/// The columns of [`OUTCOMES_FILE`].
pub const OUTCOMES_COLUMNS: [&str; 2] = ["outcome", "meaning"];
/// The columns of [`SETS_FILE`].
pub const SETS_COLUMNS: [&str; 4] = ["set", "release", "cache_version", "caches"];

/// The `OS` byte of a gzip header that names no operating system.
const GZIP_OS_UNKNOWN: u8 = 255;

/// One thing a file breaks. `line` counts from 1 in the file's text, decompressed for the
/// bacterial file; 0 marks a file that cannot be read or written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Failure {
    pub file: PathBuf,
    pub line: usize,
    pub what: String,
}

impl Failure {
    fn new(file: &Path, line: usize, what: impl Into<String>) -> Failure {
        Failure {
            file: file.to_path_buf(),
            line,
            what: what.into(),
        }
    }
}

impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}: {}", self.file.display(), self.line, self.what)
    }
}

/// How many rows the two row files hold and how many sets those rows name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Summary {
    pub rows: usize,
    pub sets: usize,
}

impl fmt::Display for Summary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} rows in {} sets", self.rows, self.sets)
    }
}

/// What a [`put`] changed: the table afterwards, the rows whose key was new and the rows
/// that took the place of one with the same key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Put {
    pub summary: Summary,
    pub added: usize,
    pub replaced: usize,
}

struct Placed {
    row: Row,
    file: PathBuf,
    line: usize,
}

struct SetRow {
    set: String,
    caches: u64,
    line: usize,
}

/// The lines after a header that names exactly `columns`, each with its line number, or
/// `None` with a failure when the file cannot be read or opens with another line.
fn body(
    path: &Path,
    columns: &[&str],
    failures: &mut Vec<Failure>,
) -> Option<Vec<(usize, String)>> {
    let lines = match reader::lines(path) {
        Ok(lines) => lines,
        Err(error) => {
            failures.push(Failure::new(path, 0, format!("cannot be read: {error}")));
            return None;
        }
    };
    if lines.first().map(String::as_str) != Some(columns.join("\t").as_str()) {
        failures.push(Failure::new(
            path,
            1,
            format!(
                "the first line is not the tab-separated header {}",
                columns.join(" ")
            ),
        ));
        return None;
    }
    Some(
        lines
            .into_iter()
            .enumerate()
            .skip(1)
            .map(|(index, line)| (index + 1, line))
            .collect(),
    )
}

/// The rows of one row file that pass the field rules; every other line is a failure.
fn read_rows(path: &Path, failures: &mut Vec<Failure>) -> Option<Vec<Placed>> {
    let mut rows = Vec::new();
    for (line, text) in body(path, &row::COLUMNS, failures)? {
        match Row::parse(&text) {
            Ok(row) => rows.push(Placed {
                row,
                file: path.to_path_buf(),
                line,
            }),
            Err(what) => failures.push(Failure::new(path, line, what)),
        }
    }
    Some(rows)
}

/// The rows of both row files, each held to its file and no key held twice; with
/// `ordered`, each file's rows also stand in bytewise key order. `None` when a file cannot
/// be read or lacks its header.
fn read_table(dir: &Path, ordered: bool, failures: &mut Vec<Failure>) -> Option<Vec<Placed>> {
    let files = [(GENOMES_FILE, false), (BACTERIA_FILE, true)];
    let mut table = Vec::new();
    let mut seen: HashMap<String, (PathBuf, usize)> = HashMap::new();
    let mut whole = true;
    for (name, bacterial) in files {
        let path = dir.join(name);
        let Some(rows) = read_rows(&path, failures) else {
            whole = false;
            continue;
        };
        let mut previous: Option<String> = None;
        for placed in rows {
            let key = placed.row.key();
            if placed.row.is_bacterial() != bacterial {
                let home = if bacterial {
                    GENOMES_FILE
                } else {
                    BACTERIA_FILE
                };
                failures.push(Failure::new(
                    &path,
                    placed.line,
                    format!("a row of the set {} belongs in {home}", placed.row.set),
                ));
            }
            if let Some((file, line)) = seen.get(&key) {
                failures.push(Failure::new(
                    &path,
                    placed.line,
                    format!("the key {key} is also at {}:{line}", file.display()),
                ));
            } else {
                seen.insert(key.clone(), (path.clone(), placed.line));
                if ordered && previous.as_ref().is_some_and(|before| key < *before) {
                    failures.push(Failure::new(
                        &path,
                        placed.line,
                        format!(
                            "out of key order: {key} follows {}",
                            previous.as_deref().unwrap_or_default()
                        ),
                    ));
                }
            }
            previous = Some(key);
            table.push(placed);
        }
    }
    whole.then_some(table)
}

/// The outcomes the directory declares, each `differs:<name>` or `excluded:<name>` with a
/// meaning, none twice.
fn read_outcomes(dir: &Path, failures: &mut Vec<Failure>) -> Option<HashSet<String>> {
    let path = dir.join(OUTCOMES_FILE);
    let mut declared = HashSet::new();
    for (line, text) in body(&path, &OUTCOMES_COLUMNS, failures)? {
        let fields: Vec<&str> = text.split('\t').collect();
        if fields.len() != OUTCOMES_COLUMNS.len() {
            failures.push(Failure::new(
                &path,
                line,
                format!(
                    "{} fields where the file has {} columns",
                    fields.len(),
                    OUTCOMES_COLUMNS.len()
                ),
            ));
            continue;
        }
        let (outcome, meaning) = (fields[0], fields[1]);
        if !matches!(
            row::outcome_kind(outcome),
            Some(OutcomeKind::Differs | OutcomeKind::Excluded)
        ) {
            failures.push(Failure::new(
                &path,
                line,
                format!("outcome \"{outcome}\" is not differs:<name> or excluded:<name>"),
            ));
        } else if meaning.trim().is_empty() {
            failures.push(Failure::new(
                &path,
                line,
                format!("outcome {outcome} has no meaning"),
            ));
        } else if !declared.insert(outcome.to_owned()) {
            failures.push(Failure::new(
                &path,
                line,
                format!("outcome {outcome} is declared twice"),
            ));
        }
    }
    Some(declared)
}

/// The sets the directory lists with their cache counts, none twice.
fn read_sets(dir: &Path, failures: &mut Vec<Failure>) -> Option<Vec<SetRow>> {
    let path = dir.join(SETS_FILE);
    let mut sets: Vec<SetRow> = Vec::new();
    for (line, text) in body(&path, &SETS_COLUMNS, failures)? {
        let fields: Vec<&str> = text.split('\t').collect();
        if fields.len() != SETS_COLUMNS.len() {
            failures.push(Failure::new(
                &path,
                line,
                format!(
                    "{} fields where the file has {} columns",
                    fields.len(),
                    SETS_COLUMNS.len()
                ),
            ));
            continue;
        }
        let (set, release, cache_version, caches) = (fields[0], fields[1], fields[2], fields[3]);
        let what = if !row::is_set(set) {
            Some(format!("set \"{set}\" is not {}", row::set_forms()))
        } else if !row::is_name(release) {
            Some(format!(
                "release \"{release}\" is not letters, digits, '.', '_' and '-'"
            ))
        } else if !row::is_positive(cache_version) {
            Some(format!(
                "cache_version \"{cache_version}\" is not an integer above zero"
            ))
        } else if !row::is_count(caches) {
            Some(format!("caches \"{caches}\" is not an unsigned integer"))
        } else if sets.iter().any(|listed| listed.set == set) {
            Some(format!("set {set} is listed twice"))
        } else {
            None
        };
        match what {
            Some(what) => failures.push(Failure::new(&path, line, what)),
            None => sets.push(SetRow {
                set: set.to_owned(),
                caches: caches.parse().unwrap_or(u64::MAX),
                line,
            }),
        }
    }
    Some(sets)
}

fn summary(table: &[Placed]) -> Summary {
    let sets: BTreeSet<&str> = table.iter().map(|placed| placed.row.set.as_str()).collect();
    Summary {
        rows: table.len(),
        sets: sets.len(),
    }
}

/// The failures in the order of `files`, then in line order. A failure of a file that
/// `files` leaves out comes after the others.
fn in_file_order(mut failures: Vec<Failure>, files: &[PathBuf]) -> Vec<Failure> {
    failures.sort_by_key(|failure| {
        let rank = files
            .iter()
            .position(|file| *file == failure.file)
            .unwrap_or(files.len());
        (rank, failure.line)
    });
    failures
}

/// Checks the directory's four files. Each file opens with its header and every line has
/// its columns in their grammar; no key stands twice; each row is in the file its set
/// names, in bytewise key order; every `differs:` and `excluded:` value a row uses is
/// declared. With `complete`, every listed set has exactly its `caches` rows and every row
/// names a listed set. The failures come in file order, then line order.
pub fn check(dir: &Path, complete: bool) -> Result<Summary, Vec<Failure>> {
    let mut failures = Vec::new();
    let table = read_table(dir, true, &mut failures);
    let outcomes = read_outcomes(dir, &mut failures);
    let sets = read_sets(dir, &mut failures);

    if let (Some(table), Some(outcomes)) = (&table, &outcomes) {
        for placed in table {
            let outcome = &placed.row.outcome;
            if outcome != "exact" && !outcomes.contains(outcome) {
                failures.push(Failure::new(
                    &placed.file,
                    placed.line,
                    format!("outcome {outcome} is not declared in {OUTCOMES_FILE}"),
                ));
            }
        }
    }

    if let (true, Some(table), Some(sets)) = (complete, &table, &sets) {
        let mut counts: BTreeMap<&str, (u64, &Placed)> = BTreeMap::new();
        for placed in table {
            counts
                .entry(placed.row.set.as_str())
                .or_insert((0, placed))
                .0 += 1;
        }
        let path = dir.join(SETS_FILE);
        for listed in sets {
            let rows = counts.get(listed.set.as_str()).map_or(0, |entry| entry.0);
            if rows != listed.caches {
                failures.push(Failure::new(
                    &path,
                    listed.line,
                    format!(
                        "set {} has {rows} rows where caches is {}",
                        listed.set, listed.caches
                    ),
                ));
            }
        }
        for (set, (rows, first)) in counts {
            if !sets.iter().any(|listed| listed.set == set) {
                failures.push(Failure::new(
                    &first.file,
                    first.line,
                    format!("set {set} is not in {SETS_FILE} ({rows} rows name it)"),
                ));
            }
        }
    }

    if failures.is_empty() {
        return Ok(summary(&table.unwrap_or_default()));
    }
    let files = [GENOMES_FILE, BACTERIA_FILE, OUTCOMES_FILE, SETS_FILE].map(|name| dir.join(name));
    Err(in_file_order(failures, &files))
}

/// The rows of both row files in bytewise key order. A directory whose row files fail a
/// rule [`check`] holds them to gives its failures in place of rows, in [`check`]'s order.
pub fn load(dir: &Path) -> Result<Vec<Row>, Vec<Failure>> {
    let mut failures = Vec::new();
    let table = read_table(dir, true, &mut failures);
    if !failures.is_empty() {
        let files = [GENOMES_FILE, BACTERIA_FILE].map(|name| dir.join(name));
        return Err(in_file_order(failures, &files));
    }
    let mut rows: Vec<Row> = table
        .unwrap_or_default()
        .into_iter()
        .map(|placed| placed.row)
        .collect();
    rows.sort_by_cached_key(Row::key);
    Ok(rows)
}

fn gzip(bytes: &[u8]) -> io::Result<Vec<u8>> {
    let mut encoder = GzBuilder::new()
        .mtime(0)
        .operating_system(GZIP_OS_UNKNOWN)
        .write(Vec::new(), Compression::best());
    encoder.write_all(bytes)?;
    encoder.finish()
}

/// The path a row file is written to before it is renamed over `path`. It carries the
/// process id, so two runs over one directory never write one file.
fn partial_path(path: &Path) -> PathBuf {
    let mut partial = path.as_os_str().to_owned();
    partial.push(format!(".{}.tmp", process::id()));
    PathBuf::from(partial)
}

/// Writes each of `files` beside its target, then renames each over its target. `staged`
/// gains every file this call creates, whether or not it reaches its target.
fn replace(files: &[(PathBuf, &[u8])], staged: &mut Vec<PathBuf>) -> Result<(), Failure> {
    let unwritten = |path: &Path, error: io::Error| {
        Failure::new(path, 0, format!("cannot be written: {error}"))
    };
    for (path, bytes) in files {
        let partial = partial_path(path);
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&partial)
            .map_err(|error| unwritten(path, error))?;
        staged.push(partial);
        file.write_all(bytes)
            .map_err(|error| unwritten(path, error))?;
    }
    for (path, _) in files {
        fs::rename(partial_path(path), path).map_err(|error| unwritten(path, error))?;
    }
    Ok(())
}

/// Writes both row files from `rows`, which iterate in key order. Each file is written
/// beside its target and renamed over it, so a reader never sees half a file, and a write
/// that fails leaves no file beside a target.
fn write_table(dir: &Path, rows: &BTreeMap<String, Row>) -> Result<(), Failure> {
    let mut plain = row::header() + "\n";
    let mut bacterial = plain.clone();
    for row in rows.values() {
        let text = if row.is_bacterial() {
            &mut bacterial
        } else {
            &mut plain
        };
        text.push_str(&row.line());
        text.push('\n');
    }
    let bacteria = dir.join(BACTERIA_FILE);
    let compressed = gzip(bacterial.as_bytes())
        .map_err(|error| Failure::new(&bacteria, 0, format!("cannot be written: {error}")))?;
    let files = [
        (dir.join(GENOMES_FILE), plain.as_bytes()),
        (bacteria, compressed.as_slice()),
    ];
    let mut staged = Vec::new();
    let written = replace(&files, &mut staged);
    if written.is_err() {
        for partial in &staged {
            let _ = fs::remove_file(partial);
        }
    }
    written
}

/// Puts the rows of `row_files` into the directory's table. A row file is the table's
/// header and rows, plain or gzip; a row takes the place of the row with its key, the last
/// one given winning, and both row files are rewritten in key order. The compressed file
/// carries no name and a zero time, so equal rows give equal bytes. Nothing is written
/// when a row of the table or of a row file fails the field rules, or when the table holds
/// a key twice or a row in the other file.
pub fn put(dir: &Path, row_files: &[PathBuf]) -> Result<Put, Vec<Failure>> {
    let mut failures = Vec::new();
    let table = read_table(dir, false, &mut failures);
    let mut incoming = Vec::new();
    for path in row_files {
        incoming.extend(read_rows(path, &mut failures).unwrap_or_default());
    }
    if !failures.is_empty() {
        let mut files = vec![dir.join(GENOMES_FILE), dir.join(BACTERIA_FILE)];
        files.extend(row_files.iter().cloned());
        return Err(in_file_order(failures, &files));
    }
    let mut rows: BTreeMap<String, Row> = table
        .unwrap_or_default()
        .into_iter()
        .map(|placed| (placed.row.key(), placed.row))
        .collect();
    let held: HashSet<String> = rows.keys().cloned().collect();
    let mut given: HashSet<String> = HashSet::new();
    for placed in incoming {
        let key = placed.row.key();
        given.insert(key.clone());
        rows.insert(key, placed.row);
    }
    let replaced = given.intersection(&held).count();
    write_table(dir, &rows).map_err(|failure| vec![failure])?;
    let sets: BTreeSet<&str> = rows.values().map(|row| row.set.as_str()).collect();
    Ok(Put {
        summary: Summary {
            rows: rows.len(),
            sets: sets.len(),
        },
        added: given.len() - replaced,
        replaced,
    })
}
