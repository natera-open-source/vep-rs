// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `synonym_of_several_sequences_by_hash_order` over the discordant rows: the
//! reference's rows of a record whose sequence name the cache holds no sequence by and the
//! cache's synonyms file gives to two or more sequences it holds. The file does not say
//! which of them the name stands for. Ensembl VEP reads such a record on whichever of them
//! comes first in the order `keys` returns the name's synonyms (`BaseVEP.pm:711-716`),
//! which changes between runs of one command, so its rows for the record are not
//! reproducible; vep-rs skips the record with a warning. The rows alone cannot tell such a
//! record from one vep-rs lost, so the class is applied only where the comparison's
//! collector holds the names read from the cache the reference run read.

use super::one_sided::{OneSidedStats, Row};
use super::scope::{Group, ReferenceRecordScope};
use crate::input::{for_each_record, for_each_row_at};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io;
use std::path::Path;

/// The synonyms file of an Ensembl VEP cache, beside its sequence directories.
pub const SYNONYMS_FILE: &str = "chr_synonyms.txt";

/// The characters Perl's `\s` matches in a file read as bytes.
fn is_perl_space(c: char) -> bool {
    c.is_ascii_whitespace() || c == '\u{0B}'
}

/// The names a synonyms file's `text` gives to two or more of `held`, the sequences a
/// cache holds, among the names the cache holds no sequence by. The file is read as both
/// engines read it (`BaseVEP.pm:646-656`): the first field of a line and each later field
/// name each other, one hop and in both directions, and a line that opens with whitespace
/// names the empty string first.
pub fn names_of_several_sequences(text: &str, held: &HashSet<String>) -> HashSet<String> {
    let mut held_partners: HashMap<&str, HashSet<&str>> = HashMap::new();
    for line in text.split('\n') {
        let (reference, rest) = line.split_once(is_perl_space).unwrap_or((line, ""));
        for synonym in rest.split(is_perl_space).filter(|name| !name.is_empty()) {
            for (name, partner) in [(reference, synonym), (synonym, reference)] {
                if !held.contains(name) && held.contains(partner) {
                    held_partners.entry(name).or_default().insert(partner);
                }
            }
        }
    }
    held_partners
        .into_iter()
        .filter(|(_, partners)| partners.len() >= 2)
        .map(|(name, _)| name.to_string())
        .collect()
}

/// The sequences an Ensembl VEP cache directory holds, as Ensembl VEP lists them
/// (`CacheDir.pm:431`): the names of its subdirectories that open with no dot.
fn held_sequences(directory: &Path) -> io::Result<HashSet<String>> {
    let mut held = HashSet::new();
    for entry in std::fs::read_dir(directory)? {
        let entry = entry?;
        let name = entry.file_name().to_string_lossy().into_owned();
        if !name.starts_with('.') && entry.path().is_dir() {
            held.insert(name);
        }
    }
    Ok(held)
}

/// [`names_of_several_sequences`] of the Ensembl VEP cache at `directory`, the one the
/// reference run read: its synonyms file (`CacheDir.pm:254`) against the sequences it holds.
/// A cache without the file gives no name. The names are the reference's own, so a sequence
/// vep-rs's converted cache lacks is still a held one here.
pub fn shared_synonym_names(directory: &Path) -> io::Result<HashSet<String>> {
    let with_path = |path: &Path, error: io::Error| {
        io::Error::new(error.kind(), format!("{}: {error}", path.display()))
    };
    let held = held_sequences(directory).map_err(|error| with_path(directory, error))?;
    let path = directory.join(SYNONYMS_FILE);
    match std::fs::read(&path) {
        Ok(bytes) => Ok(names_of_several_sequences(
            &String::from_utf8_lossy(&bytes),
            &held,
        )),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(HashSet::new()),
        Err(error) => Err(with_path(&path, error)),
    }
}

/// The sequence name of a Location: the text before its last colon, since a name may hold a
/// colon and the coordinates hold none.
fn sequence_name(location: &str) -> &str {
    location.rsplit_once(':').map_or(location, |(name, _)| name)
}

/// The header line in which Ensembl VEP states the command line of a run.
const COMMAND_LINE: &str = "## VEP command-line:";

/// Whether the reference output at `path` states, in its header, a command line that names
/// no `--synonyms`: only such a run read the cache's own synonyms file (`CacheDir.pm:254`
/// names it where the command gives none). An output that states no command line, as one
/// written under `--no_headers`, is not known to have read it.
fn read_the_caches_synonyms(path: &Path) -> io::Result<bool> {
    let with_path = |e: io::Error| io::Error::new(e.kind(), format!("{}: {e}", path.display()));
    let output = File::open(path).map_err(with_path)?;
    let mut answer = false;
    let mut in_header = true;
    for_each_record(output, |record| {
        if !in_header {
            return Ok(());
        }
        if record.first() != Some(&b'#') {
            in_header = false;
            return Ok(());
        }
        let line = String::from_utf8_lossy(record);
        if let Some(command) = line.strip_prefix(COMMAND_LINE) {
            let names_a_file = |word: &str| word == "--synonyms" || word.starts_with("--synonyms=");
            answer = !command.split_whitespace().any(names_a_file);
        }
        Ok(())
    })
    .map_err(with_path)?;
    Ok(answer)
}

/// The groups of `candidates` at which the output directory `rust_dir` holds a row. A row
/// both engines wrote is not among the discordant rows, so each output that holds a
/// candidate group is read once more.
fn groups_written_at<'a>(
    candidates: impl Iterator<Item = &'a Group>,
    rust_dir: &Path,
) -> io::Result<HashSet<Group>> {
    let mut asked: HashMap<&str, HashSet<&[u8]>> = HashMap::new();
    for (file_name, location) in candidates {
        asked
            .entry(file_name.as_str())
            .or_default()
            .insert(location.as_bytes());
    }
    let mut written = HashSet::new();
    for (file_name, locations) in asked {
        let path = rust_dir.join(file_name);
        let with_path = |e: io::Error| io::Error::new(e.kind(), format!("{}: {e}", path.display()));
        let output = File::open(&path).map_err(with_path)?;
        for_each_row_at(output, &locations, |key, _| {
            if let Some(location) = key.split(|b| *b == b'\t').next() {
                let location = String::from_utf8_lossy(location).into_owned();
                written.insert((file_name.to_string(), location));
            }
        })
        .map_err(with_path)?;
    }
    Ok(written)
}

/// Marks the rows of the class in `masked` and returns their statistics, or `None` for a
/// comparison whose collector holds no such name.
///
/// The rule reads the two output directories, `outputs`, and applies nothing without them.
/// The reference's rows at a (file, Location) are of the class when the Location's sequence
/// name is one the synonyms file gives to two or more held sequences, no earlier rule set
/// them aside, the reference output's header states a command line that names no
/// `--synonyms`, and vep-rs wrote no row there: none among `rows` and none that the
/// reference wrote too, which vep-rs's output is read for. A record vep-rs wrote a row for
/// is one it did not skip, and stays a difference.
pub fn apply_shared_synonym_rule(
    rows: &[Row],
    masked: &mut HashSet<usize>,
    scope: &ReferenceRecordScope,
    outputs: Option<(&Path, &Path)>,
) -> io::Result<Option<OneSidedStats>> {
    let Some(names) = scope.shared_synonyms.as_ref() else {
        return Ok(None);
    };
    let mut stats = OneSidedStats::default();
    let Some((reference_dir, rust_dir)) = outputs else {
        return Ok(Some(stats));
    };
    let mut reference_rows: HashMap<Group, Vec<usize>> = HashMap::new();
    let mut rust_wrote: HashSet<Group> = HashSet::new();
    for (idx, row) in rows.iter().enumerate() {
        if !names.contains(sequence_name(&row.location)) {
            continue;
        }
        let group = (row.file_name.clone(), row.location.clone());
        match row.source.as_str() {
            "missing_in_rust" if !masked.contains(&idx) => {
                reference_rows.entry(group).or_default().push(idx);
            }
            "extra_in_rust" => {
                rust_wrote.insert(group);
            }
            _ => {}
        }
    }
    reference_rows.retain(|group, _| !rust_wrote.contains(group));
    let mut read_the_file: HashMap<String, bool> = HashMap::new();
    for (file_name, _) in reference_rows.keys() {
        if !read_the_file.contains_key(file_name) {
            let reads = read_the_caches_synonyms(&reference_dir.join(file_name))?;
            read_the_file.insert(file_name.clone(), reads);
        }
    }
    reference_rows.retain(|(file_name, _), _| read_the_file[file_name]);
    if !reference_rows.is_empty() {
        let written = groups_written_at(reference_rows.keys(), rust_dir)?;
        reference_rows.retain(|group, _| !written.contains(group));
    }
    for of_group in reference_rows.into_values() {
        stats.records += 1;
        stats.perl += of_group.len() as u64;
        masked.extend(of_group);
    }
    stats.count = stats.perl;
    stats.excluded = stats.perl;
    Ok(Some(stats))
}

#[cfg(test)]
#[path = "shared_synonym_tests.rs"]
mod tests;
