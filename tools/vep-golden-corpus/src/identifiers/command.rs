// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `identifiers` subcommand: every identifier under restricted redistribution terms in
//! the named files and directories replaced by the fabricated one of the same shape, and the
//! digests a directory's `provenance.json` holds of a rewritten file recorded again.

use super::fabricated_identifiers;
use super::gzip::compress;
use super::stored::{
    corpus_files, read_corpus_text, restricted_identifier_counts, restricted_identifier_errors,
    CorpusText, Counted,
};
use crate::cli::IdentifiersArgs;
use crate::error::{complain, failure, say, Failure};
use crate::pytext::{is_file, os_error_text, path_join, path_str, read_text, sha256_hex};
use std::collections::BTreeSet;

/// One named path read: the path, whether it is a directory, the files to rewrite with
/// their texts and counts, and the identifiers read across its files.
struct Held {
    root: String,
    directory: bool,
    holding: Vec<(String, CorpusText, usize)>,
    seen: BTreeSet<String>,
}

/// The directory a path stands in, as `pathlib` gives the parent of a path.
fn parent(path: &str) -> &str {
    match path.rfind('/') {
        None => ".",
        Some(slash) if path[..slash].chars().all(|c| c == '/') => &path[..=slash],
        Some(slash) => &path[..slash],
    }
}

/// The `provenance.json` that records a file under the directory `root`
/// (`recording_provenance`): the nearest one in the file's own directory or above it, no
/// higher than `root`.
fn recording_provenance(path: &str, root: &str) -> Option<String> {
    let mut directory = parent(path);
    loop {
        let record = path_join(directory, "provenance.json");
        if is_file(&record) {
            return Some(record);
        }
        let above = parent(directory);
        if directory == root || above == directory {
            return None;
        }
        directory = above;
    }
}

/// Sets the digest `before` to `after` among `changed`, a digest named twice keeping its
/// first place.
fn record_digest(changed: &mut Vec<(String, String)>, before: String, after: String) {
    match changed.iter_mut().find(|(held, _)| *held == before) {
        Some((_, held)) => *held = after,
        None => changed.push((before, after)),
    }
}

fn not_writable(path: &str, error: &std::io::Error) -> Failure {
    Failure::own(format!(
        "{path} is not writable ({})",
        os_error_text(error, path)
    ))
}

/// Reads one named path: every file under a directory, or the one file.
fn hold(name: &str, check: bool, status: &mut i32) -> Result<Option<Held>, Failure> {
    let root = path_str(name);
    let Ok(metadata) = std::fs::metadata(&root) else {
        return Err(failure(format!("{root} is neither a file nor a directory")));
    };
    let directory = metadata.is_dir();
    let files = if directory {
        corpus_files(&root)?
    } else {
        vec![root.clone()]
    };
    if check {
        let errors = restricted_identifier_errors(&files, None)?;
        for error in &errors {
            complain(error);
        }
        say(&format!(
            "{root}: {} files read, {} hold an identifier in another form than the fabricated one",
            files.len(),
            errors.len()
        ));
        if !errors.is_empty() {
            *status = 1;
        }
        return Ok(None);
    }
    let Counted { counts, seen } = restricted_identifier_counts(&files)?;
    let mut holding = Vec::new();
    for (path, count) in counts {
        let Some(read) = read_corpus_text(&path)? else {
            continue;
        };
        if read.blocked {
            let plural = if count == 1 { "" } else { "s" };
            return Err(failure(format!(
                "{path} holds {count} identifier{plural} under restricted redistribution terms \
                 and is block-compressed: decompress it, pass the text through this subcommand, \
                 then compress and index it again; nothing was written"
            )));
        }
        holding.push((path, read, count));
    }
    Ok(Some(Held {
        root,
        directory,
        holding,
        seen,
    }))
}

/// Rewrites the files of one named path and records their digests again.
fn rewrite(held: &Held) -> Result<(), Failure> {
    // A record names a file by the digest of its bytes or of its uncompressed bytes.
    let mut digests: Vec<(String, Vec<(String, String)>)> = Vec::new();
    for (path, read, _) in &held.holding {
        let text = fabricated_identifiers(&read.text, &held.seen).into_owned();
        let stored = if read.gzipped {
            compress(text.as_bytes())
        } else {
            text.clone().into_bytes()
        };
        std::fs::write(path, &stored).map_err(|error| not_writable(path, &error))?;
        let record = held
            .directory
            .then(|| recording_provenance(path, &held.root));
        let Some(record) = record.flatten().filter(|record| record != path) else {
            continue;
        };
        let place = digests.iter().position(|(held, _)| *held == record);
        let place = place.unwrap_or_else(|| {
            digests.push((record, Vec::new()));
            digests.len() - 1
        });
        let changed = &mut digests[place].1;
        record_digest(changed, sha256_hex(&read.raw), sha256_hex(&stored));
        let (before, after) = (read.text.as_bytes(), text.as_bytes());
        record_digest(changed, sha256_hex(before), sha256_hex(after));
    }
    let mut recorded = 0;
    for (record, changed) in &digests {
        let mut text = read_text(record)
            .map_err(|reason| Failure::own(format!("{record} is not readable text ({reason})")))?;
        for (before, after) in changed {
            recorded += text.matches(before.as_str()).count();
            text = text.replace(before.as_str(), after);
        }
        std::fs::write(record, text).map_err(|error| not_writable(record, &error))?;
    }
    let replaced: usize = held.holding.iter().map(|(_, _, count)| count).sum();
    say(&format!(
        "{}: {replaced} identifiers replaced in {} files, {recorded} digests recorded again",
        held.root,
        held.holding.len()
    ));
    Ok(())
}

/// Runs `identifiers`. Every file of every named path is read before one is written, so a
/// refusal leaves each tree as it was.
pub fn run(args: &IdentifiersArgs) -> Result<i32, Failure> {
    let mut status = 0;
    let mut held = Vec::new();
    for name in &args.path {
        held.extend(hold(name, args.check, &mut status)?);
    }
    for held in &held {
        rewrite(held)?;
    }
    Ok(status)
}

#[cfg(test)]
#[path = "command_tests.rs"]
mod tests;
