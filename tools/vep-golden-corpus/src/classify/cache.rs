// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The corpus cache: its directory under either layout and every transcript record of its
//! shards.

use crate::constants::DEFAULT_SPECIES;
use crate::error::Failure;
use crate::json::{loads, JsonError, Object, Value};
use crate::pytext::{is_dir, path_join, py_str, read_text};
use crate::rows::{ByteSource, ReadError};

/// Whether a manifest's `flags` names `flag`, as `flag in flags` reads a list, a text or an
/// object.
fn names_flag(flags: Option<&Value>, flag: &str) -> bool {
    match flags {
        Some(Value::List(items)) => items.iter().any(|item| item.as_str() == Some(flag)),
        Some(Value::Str(text)) => text.contains(flag),
        Some(Value::Object(object)) => object.contains_key(flag),
        _ => false,
    }
}

/// The corpus's cache directory: `json_cache/`, or under the Ensembl layout
/// `cache_root/<species>/<release>_<assembly>` (`homo_sapiens` when the manifest names no
/// species; the species suffixed `_refseq` and then `_merged` for each of `--refseq` and
/// `--merged` the manifest's flags carry, as Ensembl VEP names the directory, CacheDir.pm:312
/// at release/116.3), the directory `crates/vep-cli/tests/common/mod.rs` resolves.
pub fn corpus_cache_dir(corpus: &str, manifest: &Object) -> String {
    if manifest.get("cache_layout").and_then(Value::as_str) != Some("ensembl") {
        return path_join(corpus, "json_cache");
    }
    let mut species = manifest
        .get("species")
        .map_or_else(|| DEFAULT_SPECIES.to_string(), py_str);
    let flags = manifest.get("flags");
    if names_flag(flags, "--refseq") {
        species.push_str("_refseq");
    }
    if names_flag(flags, "--merged") {
        species.push_str("_merged");
    }
    let text = |key: &str| manifest.get(key).map(py_str).unwrap_or_default();
    let version = format!("{}_{}", text("release"), text("assembly"));
    path_join(
        &path_join(&path_join(corpus, "cache_root"), &species),
        &version,
    )
}

/// Every transcript record of a JSON cache with the chromosome it is filed under
/// (`iter_cache_transcripts`).
pub struct CacheTranscripts {
    /// The chromosome directories not read yet, each with its name and its path.
    directories: std::vec::IntoIter<(String, String)>,
    /// The shards of the directory being read that are not read yet.
    shards: std::vec::IntoIter<String>,
    chromosome: String,
    /// The records of the shard being read that are not given yet.
    records: std::vec::IntoIter<Value>,
    failure: Option<Failure>,
    finished: bool,
}

/// Why the entries of a directory are not listed.
enum Unlisted {
    /// The directory cannot be read.
    Directory(std::io::Error),
    /// An entry that would be kept has a name that is not UTF-8.
    Name,
}

impl Unlisted {
    fn failure(&self, directory: &str) -> Failure {
        match self {
            Unlisted::Directory(error) => {
                Failure::own(format!("{directory} cannot be listed ({error})"))
            }
            Unlisted::Name => Failure::own(format!("{directory} holds a name that is not UTF-8")),
        }
    }
}

/// The entries of `directory` that `keep` takes by name, as `(name, path)` in the order of
/// their names.
fn sorted_entries(
    directory: &str,
    keep: impl Fn(&str) -> bool,
) -> Result<Vec<(String, String)>, Unlisted> {
    let mut entries = Vec::new();
    for entry in std::fs::read_dir(directory).map_err(Unlisted::Directory)? {
        let name = match entry
            .map_err(Unlisted::Directory)?
            .file_name()
            .into_string()
        {
            Ok(name) => name,
            Err(name) if keep(&name.to_string_lossy()) => return Err(Unlisted::Name),
            Err(_) => continue,
        };
        if keep(&name) {
            let path = path_join(directory, &name);
            entries.push((name, path));
        }
    }
    entries.sort();
    Ok(entries)
}

/// `(chromosome, transcript record)` for every transcript of a JSON cache. A shard is filed
/// as `{"chr": [record, null, ...]}` or as a list, plain (`*.json`) or gzipped (`*.json.gz`);
/// the chromosome is the name of the shard's directory under `<cache_dir>/transcripts`. The
/// directories are read in the order of their names, the shards of each in the order of
/// theirs and the records of a shard in the order written; a record is an object with a
/// `stable_id` that is not empty. A cache without a `transcripts` directory holds no record,
/// a directory that cannot be listed holds no shard, and a shard that cannot be read as
/// UTF-8 JSON is skipped.
///
/// An item is an `Err` for what the tool refuses in its own words: a `transcripts`
/// directory that cannot be listed, a name that is not UTF-8, deflate data that does not
/// inflate, a shard outside what the JSON reader holds, and a number that stands where a
/// shard keeps the records of a chromosome. No item follows an `Err`.
pub fn iter_cache_transcripts(cache_dir: &str) -> CacheTranscripts {
    let root = path_join(cache_dir, "transcripts");
    let (directories, failure) = if is_dir(&root) {
        match sorted_entries(&root, |_| true) {
            Ok(directories) => (directories, None),
            Err(unlisted) => (Vec::new(), Some(unlisted.failure(&root))),
        }
    } else {
        (Vec::new(), None)
    };
    CacheTranscripts {
        directories: directories.into_iter(),
        shards: Vec::new().into_iter(),
        chromosome: String::new(),
        records: Vec::new().into_iter(),
        failure,
        finished: false,
    }
}

/// The text of a shard, `None` for one the reader skips.
fn shard_text(path: &str) -> Result<Option<String>, Failure> {
    if !path.ends_with(".gz") {
        return Ok(read_text(path).ok());
    }
    let unreadable = |reason: &str| Failure::own(format!("{path} is not readable ({reason})"));
    match ByteSource::open_gzip(path).and_then(ByteSource::read_to_end) {
        Ok(bytes) => Ok(String::from_utf8(bytes).ok()),
        // CPython raises `zlib.error` for deflate data that does not inflate, in zlib's
        // words, and the reader of a cache does not catch it.
        Err(ReadError::Cpython(reason)) if reason.starts_with("Error -") => {
            Err(unreadable(&reason))
        }
        Err(ReadError::Tool(reason)) => Err(unreadable(&reason)),
        Err(ReadError::Cpython(_)) => Ok(None),
    }
}

/// The records of a shard in the order written, `None` for a shard the reader skips.
fn shard_records(path: &str) -> Result<Option<Vec<Value>>, Failure> {
    let Some(text) = shard_text(path)? else {
        return Ok(None);
    };
    let data = match loads(&text) {
        Ok(data) => data,
        Err(JsonError::Decode { .. }) => return Ok(None),
        // `int()` refuses the digits with a `ValueError`, which skips the shard.
        Err(JsonError::Unsupported(detail)) if detail.starts_with("an integer") => return Ok(None),
        Err(JsonError::Unsupported(detail)) => {
            return Err(Failure::own(format!("{path} holds {detail}")))
        }
    };
    match data {
        Value::List(records) => Ok(Some(records)),
        Value::Object(mut chromosomes) => {
            let mut records = Vec::new();
            for (name, held) in chromosomes.iter_mut() {
                match std::mem::replace(held, Value::Null) {
                    Value::List(items) => records.extend(items.into_iter().filter(Value::truthy)),
                    // The characters of a text and the keys of an object are no records.
                    Value::Str(_) | Value::Object(_) => {}
                    other if !other.truthy() => {}
                    _ => {
                        return Err(Failure::own(format!(
                            "{path} holds a number where the records of {name} stand"
                        )))
                    }
                }
            }
            Ok(Some(records))
        }
        _ => Ok(None),
    }
}

impl CacheTranscripts {
    fn fail(&mut self, failure: Failure) -> Option<Result<(String, Object), Failure>> {
        self.finished = true;
        Some(Err(failure))
    }
}

impl Iterator for CacheTranscripts {
    type Item = Result<(String, Object), Failure>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(failure) = self.failure.take() {
            return self.fail(failure);
        }
        loop {
            if self.finished {
                return None;
            }
            for record in self.records.by_ref() {
                let Value::Object(transcript) = record else {
                    continue;
                };
                if transcript.get("stable_id").is_some_and(Value::truthy) {
                    return Some(Ok((self.chromosome.clone(), transcript)));
                }
            }
            if let Some(shard) = self.shards.next() {
                match shard_records(&shard) {
                    Ok(Some(records)) => self.records = records.into_iter(),
                    Ok(None) => {}
                    Err(failure) => return self.fail(failure),
                }
                continue;
            }
            let Some((name, path)) = self.directories.next() else {
                self.finished = true;
                return None;
            };
            if !is_dir(&path) {
                continue;
            }
            let is_shard = |name: &str| name.ends_with(".json") || name.ends_with(".json.gz");
            let shards = match sorted_entries(&path, is_shard) {
                Ok(shards) => shards,
                Err(Unlisted::Directory(_)) => Vec::new(),
                Err(unlisted) => return self.fail(unlisted.failure(&path)),
            };
            let paths: Vec<String> = shards.into_iter().map(|(_, path)| path).collect();
            self.shards = paths.into_iter();
            self.chromosome = name;
        }
    }
}

#[cfg(test)]
#[path = "cache_tests.rs"]
mod tests;
