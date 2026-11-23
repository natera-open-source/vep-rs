// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcript shards of a JSON cache: under the cache's `transcripts/` one directory per
//! chromosome, each holding the shards `<name>.json` and `<name>.json.gz`. A shard is a bare
//! array of transcript objects or a range-keyed map of such arrays
//! (`shared::shard_transcripts`).

use super::shared::{read_shard, shard_transcripts, strerror};
use super::{Refusal, Res};
use crate::ordered::in_order;
use crate::pytext::Json;
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

/// The shards read at a time, at most: a shard is held as its text and its parsed value
/// while it is read.
const MOST_READERS: usize = 16;

/// The shards of `directory` in the order they are read: the files whose names end in
/// `.json` by name, then those ending in `.json.gz` by name. A name that opens with a dot
/// is a shard like any other; a directory that cannot be listed holds none.
pub fn shards(directory: &Path) -> Vec<PathBuf> {
    let Ok(entries) = fs::read_dir(directory) else {
        return Vec::new();
    };
    let names: Vec<_> = entries
        .filter_map(|entry| entry.ok().map(|entry| entry.file_name()))
        .collect();
    let mut found = Vec::new();
    for ending in [&b".json"[..], b".json.gz"] {
        let mut named: Vec<_> = names
            .iter()
            .filter(|name| name.as_encoded_bytes().ends_with(ending))
            .collect();
        named.sort_unstable_by_key(|name| name.as_encoded_bytes());
        found.extend(named.into_iter().map(|name| directory.join(name)));
    }
    found
}

/// The `transcripts/` directory of a JSON cache and the chromosome directories it held when
/// it was opened.
pub struct TranscriptShards {
    root: PathBuf,
    chromosomes: BTreeSet<String>,
    /// The threads a directory's shards are read on.
    readers: AtomicUsize,
}

impl TranscriptShards {
    /// The shards of the cache at `json_cache`, refused as `--json-cache: <cache> has no
    /// transcripts/ directory` where it holds none.
    pub fn new(json_cache: &Path) -> Res<TranscriptShards> {
        let root = json_cache.join("transcripts");
        if !root.is_dir() {
            return Err(Refusal(format!(
                "--json-cache: {} has no transcripts/ directory",
                json_cache.display()
            )));
        }
        let unreadable = |error| {
            Refusal(format!(
                "cannot read {}: {}",
                root.display(),
                strerror(&error)
            ))
        };
        let mut chromosomes = BTreeSet::new();
        for entry in fs::read_dir(&root).map_err(unreadable)? {
            let entry = entry.map_err(unreadable)?;
            if let Some(name) = entry.file_name().to_str() {
                if entry.path().is_dir() {
                    chromosomes.insert(name.to_string());
                }
            }
        }
        Ok(TranscriptShards {
            root,
            chromosomes,
            readers: AtomicUsize::new(1),
        })
    }

    /// Reads the shards of a directory on `threads` threads, `MOST_READERS` at most, from
    /// the next walk on. A walk hands over the same objects in the same order, and ends in
    /// the same refusal, on any number of threads.
    pub fn read_on(&self, threads: usize) {
        self.readers
            .store(threads.clamp(1, MOST_READERS), Ordering::Relaxed);
    }

    /// The cache's `transcripts/` directory.
    pub fn root(&self) -> &Path {
        &self.root
    }

    /// The names of the directories `transcripts/` held when the cache was opened, in the
    /// order of their code points.
    pub fn chromosomes(&self) -> &BTreeSet<String> {
        &self.chromosomes
    }

    /// Hands `visit` each transcript object of the shards of `transcripts/<chromosome>`, the
    /// shards in reading order (`shards`) and each shard's transcripts in the shard's. A
    /// shard that cannot be read, is not JSON or holds another shape is refused by name,
    /// and a refusal of `visit` ends the walk.
    pub fn of_chromosome(
        &self,
        chromosome: &str,
        mut visit: impl FnMut(Json) -> Res<()>,
    ) -> Res<()> {
        self.walk(chromosome, &Some, &mut visit)
    }

    /// Hands `visit` each transcript object of every chromosome directory, the directories
    /// in the order of their names and each as `of_chromosome` reads it.
    pub fn of_every_chromosome(&self, mut visit: impl FnMut(Json) -> Res<()>) -> Res<()> {
        self.chromosomes
            .iter()
            .try_for_each(|chromosome| self.walk(chromosome, &Some, &mut visit))
    }

    /// As `of_chromosome`, with `visit` handed what `read` makes of each transcript object
    /// and nothing for an object it makes `None` of. `read` runs on the threads that read
    /// the shards.
    pub fn read_of_chromosome<T: Send>(
        &self,
        chromosome: &str,
        read: impl Fn(Json) -> Option<T> + Sync,
        mut visit: impl FnMut(T) -> Res<()>,
    ) -> Res<()> {
        self.walk(chromosome, &read, &mut visit)
    }

    /// As `of_every_chromosome`, with `visit` handed what `read` makes of each transcript
    /// object, as `read_of_chromosome` hands it.
    pub fn read_of_every_chromosome<T: Send>(
        &self,
        read: impl Fn(Json) -> Option<T> + Sync,
        mut visit: impl FnMut(T) -> Res<()>,
    ) -> Res<()> {
        self.chromosomes
            .iter()
            .try_for_each(|chromosome| self.walk(chromosome, &read, &mut visit))
    }

    fn walk<T: Send>(
        &self,
        chromosome: &str,
        read: &(dyn Fn(Json) -> Option<T> + Sync),
        visit: &mut dyn FnMut(T) -> Res<()>,
    ) -> Res<()> {
        let shards = shards(&self.root.join(chromosome));
        let made = |at: usize| -> Res<Vec<T>> {
            let parsed = read_shard(&shards[at])?;
            let transcripts = shard_transcripts(&shards[at], parsed)?;
            Ok(transcripts.into_iter().filter_map(read).collect())
        };
        let mut walked = Ok(());
        in_order(
            shards.len(),
            self.readers.load(Ordering::Relaxed),
            &made,
            &mut |read: Res<Vec<T>>| {
                walked = read.and_then(|read| read.into_iter().try_for_each(&mut *visit));
                walked.is_ok()
            },
        );
        walked
    }
}
