// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The files of a corpus read for identifiers under restricted redistribution terms: every
//! file under a directory in the order CPython sorts paths, each as text where it is text,
//! a gzip file inflated first, and the count of identifiers each holds in another form than
//! the fabricated one, the files of one directory read as one corpus.

use super::gzip::{decompress, GZIP_MAGIC};
use super::{holds_identifier_bytes, restricted_identifier_count, restricted_identifiers};
use crate::error::{failure, Failure};
use crate::pytext::{os_error_text, path_join};
use crate::rows::blocks::thread_count;
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::channel;

/// A corpus file read as text: its bytes as stored, whether they are gzip, whether that
/// gzip is block-compressed (BGZF, which an index beside the file may describe), and the
/// text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorpusText {
    pub raw: Vec<u8>,
    pub gzipped: bool,
    pub blocked: bool,
    pub text: String,
}

/// A file of a corpus as text (`read_corpus_text`), inflated when its first bytes are
/// gzip's whatever its name; `None` for a file that can hold no restricted identifier: one
/// whose bytes `holds_identifier_bytes` refuses, or that is not UTF-8 text once inflated (an
/// index, a binary track). A file that cannot be read or inflated is refused with its path.
pub fn read_corpus_text(path: &str) -> Result<Option<CorpusText>, Failure> {
    let unreadable = |reason: &str| format!("{path} is not a readable file ({reason})");
    let raw =
        std::fs::read(path).map_err(|error| failure(unreadable(&os_error_text(&error, path))))?;
    let gzipped = raw.starts_with(&GZIP_MAGIC);
    let data = if gzipped {
        decompress(&raw).map_err(|error| error.failure(unreadable))?
    } else {
        raw.clone()
    };
    if !holds_identifier_bytes(&data) || data.contains(&0) {
        return Ok(None);
    }
    let Ok(text) = String::from_utf8(data) else {
        return Ok(None);
    };
    // A BGZF block is a gzip member whose extra field opens with the subfield `BC`.
    let blocked = gzipped && raw.len() > 13 && raw[3] & 4 != 0 && raw[12..14] == *b"BC";
    Ok(Some(CorpusText {
        raw,
        gzipped,
        blocked,
        text,
    }))
}

/// The entries of `directory`, each as its path below the walk's root, onto `found` when it
/// is a file; a directory that is no link is walked in turn, and one that cannot be listed
/// holds nothing.
fn walk(directory: &Path, below: &[String], found: &mut Vec<Vec<String>>) -> Result<(), Failure> {
    let Ok(entries) = std::fs::read_dir(directory) else {
        return Ok(());
    };
    for entry in entries.flatten() {
        let name = entry.file_name().into_string().map_err(|_| {
            let held = directory.display();
            Failure::own(format!("{held} holds a file whose name is not UTF-8"))
        })?;
        let path = entry.path();
        let mut parts = below.to_vec();
        parts.push(name);
        if entry.file_type().is_ok_and(|kind| kind.is_dir()) {
            walk(&path, &parts, found)?;
        } else if std::fs::metadata(&path).is_ok_and(|metadata| metadata.is_file()) {
            found.push(parts);
        }
    }
    Ok(())
}

/// Every file under the directory `root` (`corpus_files`), in the order CPython sorts
/// paths: by their parts, one part after the other. A link to a file is a file, and a link
/// to a directory is not followed.
pub fn corpus_files(root: &str) -> Result<Vec<String>, Failure> {
    let mut found = Vec::new();
    walk(Path::new(root), &[], &mut found)?;
    found.sort();
    Ok(found
        .iter()
        .map(|parts| path_join(root, &parts.join("/")))
        .collect())
}

/// `read` of each file, in the order of `files`, on the threads of the run's bound.
fn of_each_file<T: Send>(files: &[String], read: impl Fn(&str) -> T + Sync) -> Vec<T> {
    let threads = thread_count().unwrap_or(1).min(files.len());
    if threads < 2 {
        return files.iter().map(|file| read(file)).collect();
    }
    let next = AtomicUsize::new(0);
    let mut results: Vec<Option<T>> = files.iter().map(|_| None).collect();
    std::thread::scope(|scope| {
        let (sender, receiver) = channel();
        for _ in 0..threads {
            let sender = sender.clone();
            let (next, read) = (&next, &read);
            scope.spawn(move || {
                while let Some((place, file)) = {
                    let place = next.fetch_add(1, Ordering::Relaxed);
                    files.get(place).map(|file| (place, file))
                } {
                    if sender.send((place, read(file))).is_err() {
                        break;
                    }
                }
            });
        }
        drop(sender);
        for (place, result) in receiver {
            results[place] = Some(result);
        }
    });
    let results = results.into_iter();
    results
        .map(|result| result.expect("a result a file"))
        .collect()
}

/// What the files of one corpus hold (`restricted_identifier_counts`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Counted {
    /// Each file that holds an identifier under restricted redistribution terms in another
    /// form than the fabricated one, with how many it holds, in the order of the files.
    pub counts: Vec<(String, usize)>,
    /// The identifiers read across the files.
    pub seen: BTreeSet<String>,
}

/// Per file that holds one, how many identifiers under restricted redistribution terms its
/// text, inflated, holds in another form than the fabricated one, and the identifiers read
/// across the files (`restricted_identifier_counts`). A second-shape identifier that
/// `restricted_identifiers` reads in one file is counted wherever it stands as a word in
/// every file, so the files of one directory are read as one corpus.
pub fn restricted_identifier_counts(files: &[String]) -> Result<Counted, Failure> {
    let none = BTreeSet::new();
    let first = of_each_file(files, |file| {
        let read = read_corpus_text(file)?;
        Ok(read.map(|read| {
            let count = restricted_identifier_count(&read.text, &none);
            (restricted_identifiers(&read.text), count)
        }))
    });
    let mut seen = BTreeSet::new();
    let mut counts: Vec<(String, usize)> = Vec::new();
    for (file, read) in files.iter().zip(first) {
        if let Some((identifiers, count)) = read? {
            seen.extend(identifiers);
            counts.push((file.clone(), count));
        }
    }
    if seen.iter().any(|identifier| !identifier.starts_with("COS")) {
        let read: Vec<String> = counts.iter().map(|(file, _)| file.clone()).collect();
        let again = of_each_file(&read, |file| {
            let read = read_corpus_text(file)?;
            Ok(read.map_or(0, |read| restricted_identifier_count(&read.text, &seen)))
        });
        for ((_, count), counted) in counts.iter_mut().zip(again) {
            *count = counted?;
        }
    }
    counts.retain(|(_, count)| *count > 0);
    Ok(Counted { counts, seen })
}

/// A path as it reads below the directory `above`, which holds it; a path below `.` reads
/// as it is.
fn below<'p>(path: &'p str, above: &str) -> &'p str {
    match path.strip_prefix(above) {
        _ if above == "." => path,
        Some(rest) if above.ends_with('/') => rest,
        Some(rest) => rest.strip_prefix('/').unwrap_or(path),
        None => path,
    }
}

/// One line per file that `restricted_identifier_counts` counts an identifier in
/// (`restricted_identifier_errors`): the file, below `relative_to` when given, and the
/// count, never the identifier.
pub fn restricted_identifier_errors(
    files: &[String],
    relative_to: Option<&str>,
) -> Result<Vec<String>, Failure> {
    let Counted { counts, .. } = restricted_identifier_counts(files)?;
    let lines = counts.iter().map(|(path, count)| {
        let name = relative_to.map_or(path.as_str(), |above| below(path, above));
        let plural = if *count == 1 { "" } else { "s" };
        format!(
            "{name} holds {count} identifier{plural} under restricted redistribution terms in \
             another form than the fabricated one; the identifiers subcommand replaces them"
        )
    });
    Ok(lines.collect())
}

#[cfg(test)]
#[path = "stored_tests.rs"]
mod tests;
