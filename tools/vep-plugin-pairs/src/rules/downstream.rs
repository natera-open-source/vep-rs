// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The Downstream, NMD and BLOSUM62 classes and the transcript sequences they read, with
//! the coding ends of every transcript the same shards hold.

pub mod classes;
pub mod codes;
pub mod ends;
pub mod readings;
pub mod transcript;

pub use classes::{
    blosum62_shift_shared_across_alleles, reverse_strand_rotation, transcript_reading,
};
pub use transcript::{DownstreamTranscript, TranscriptModel};

use super::cache::TranscriptShards;
use super::coding_end::{location_bases, row_location};
use super::coding_ends::CodingEnds;
use super::entry::Entry;
use super::shared::without_version;
use super::{Refusal, Res};
use crate::pytext::{self, Json};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::OnceLock;

/// The transcripts of some shards by stable id.
type Transcripts = HashMap<String, DownstreamTranscript>;

/// Where the transcript objects come from: the shards of a cache, or objects handed over by
/// chromosome, in the order of the chromosomes' names.
enum Source {
    Cache(TranscriptShards),
    Objects(Vec<(String, Vec<Json>)>),
}

/// What some shards hold: their transcripts as the plugins read them, and under each
/// stable id the coding ends of its records, a reading that two records share once.
#[derive(Default)]
struct Held {
    transcripts: Transcripts,
    ends: HashMap<String, Vec<CodingEnds>>,
}

/// The transcripts read so far, in the order of the shards: a stable id the shards give one
/// reading under two models is kept without a model; one they give two different readings
/// otherwise is left out at the end. The coding ends are kept of every record.
#[derive(Default)]
struct Reading {
    held: Held,
    conflicting: HashSet<String>,
}

/// What one cache transcript object gives: its stable id, its reading for the plugins,
/// `None` for an object they read nothing of, and its coding ends.
type Read = (String, Option<DownstreamTranscript>, CodingEnds);

/// A cache transcript object as the plugins and the classes of a coding end read it;
/// `None` for an object without a stable id that is text.
fn read_transcript(transcript: &Json) -> Option<Read> {
    let (stable_id, ends) = ends::coding_ends_of(transcript)?;
    let read = DownstreamTranscript::from_cache(transcript).map(|(_, read)| read);
    Some((stable_id.to_string(), read, ends))
}

impl Reading {
    fn take(&mut self, (stable_id, read, ends): Read) {
        match self.held.ends.get_mut(&stable_id) {
            Some(records) if records.contains(&ends) => {}
            Some(records) => records.push(ends),
            None => {
                self.held.ends.insert(stable_id.clone(), vec![ends]);
            }
        }
        let Some(read) = read else {
            return;
        };
        let transcripts = &mut self.held.transcripts;
        let Some(held) = transcripts.get_mut(&stable_id) else {
            transcripts.insert(stable_id, read);
            return;
        };
        if *held == read {
            return;
        }
        if held.same_sequences(&read) {
            *held = DownstreamTranscript {
                model: None,
                ..read
            };
        } else {
            self.conflicting.insert(stable_id);
        }
    }

    fn finish(mut self) -> Held {
        for stable_id in &self.conflicting {
            self.held.transcripts.remove(stable_id);
        }
        self.held.transcripts.shrink_to_fit();
        self.held.ends.shrink_to_fit();
        self.held
    }
}

/// The transcripts of a candidate JSON cache as the Downstream and NMD plugins read them
/// (`DownstreamTranscript`), by transcript stable id, kept for a transcript with a
/// `translateable_seq`. A chromosome's shards are read on the first look-up of an entry
/// whose Location names its directory, every shard on the first look-up of an entry without
/// one. A stable id the shards read give one reading under two models (a transcript the
/// cache files under two chromosomes, each with its own mapper pairs) is kept without a
/// model; one they give two different readings otherwise is left out. The same read keeps
/// the coding ends of every record with a stable id (`ends::coding_ends_of`), a transcript
/// without a coding sequence among them.
pub struct DownstreamSequences {
    source: Source,
    /// What each chromosome directory holds, read once.
    chromosomes: HashMap<String, OnceLock<Res<Held>>>,
    /// What every directory holds together, read once.
    every: OnceLock<Res<Held>>,
}

impl DownstreamSequences {
    /// The transcripts of the cache at `json_cache`, which is refused without a
    /// `transcripts/` directory; no shard is read yet.
    pub fn new(json_cache: &Path) -> Result<DownstreamSequences, Refusal> {
        let shards = TranscriptShards::new(json_cache)?;
        let names: Vec<String> = shards.chromosomes().iter().cloned().collect();
        Ok(DownstreamSequences::over(Source::Cache(shards), names))
    }

    /// The transcripts of the objects `chromosomes` hands over, each chromosome's name with
    /// the transcript objects of its shards in reading order.
    pub fn of_transcripts(mut chromosomes: Vec<(String, Vec<Json>)>) -> DownstreamSequences {
        chromosomes.sort_by(|(one, _), (other, _)| one.cmp(other));
        let names = chromosomes.iter().map(|(name, _)| name.clone()).collect();
        DownstreamSequences::over(Source::Objects(chromosomes), names)
    }

    fn over(source: Source, names: Vec<String>) -> DownstreamSequences {
        DownstreamSequences {
            source,
            chromosomes: names
                .into_iter()
                .map(|name| (name, OnceLock::new()))
                .collect(),
            every: OnceLock::new(),
        }
    }

    /// Reads a chromosome's shards on `threads` threads from the next chromosome read on.
    /// The transcripts held and a read's refusal are the same on any number of threads.
    pub fn read_on(&self, threads: usize) {
        if let Source::Cache(shards) = &self.source {
            shards.read_on(threads);
        }
    }

    /// What the chromosome directory `chromosome` holds, what every directory holds for
    /// `None`.
    fn read(&self, chromosome: Option<&str>) -> Res<Held> {
        let mut reading = Reading::default();
        match &self.source {
            Source::Cache(shards) => {
                // The threads that parse the shards read each transcript too, and drop
                // its object; the readings are folded here in the shards' order.
                let read = |transcript: Json| read_transcript(&transcript);
                let take = |read| {
                    reading.take(read);
                    Ok(())
                };
                match chromosome {
                    Some(chromosome) => shards.read_of_chromosome(chromosome, read, take)?,
                    None => shards.read_of_every_chromosome(read, take)?,
                }
            }
            Source::Objects(chromosomes) => chromosomes
                .iter()
                .filter(|(name, _)| chromosome.is_none_or(|chromosome| chromosome == name))
                .flat_map(|(_, transcripts)| transcripts)
                .filter_map(read_transcript)
                .for_each(|read| reading.take(read)),
        }
        Ok(reading.finish())
    }

    /// The entry's transcript, or `None` for one the cache does not hold with a
    /// `translateable_seq`; a Feature carrying a version the cache's stable id lacks is
    /// looked up without it. A shard that cannot be read is refused.
    pub fn transcript(&self, e: &Entry) -> Res<Option<&DownstreamTranscript>> {
        let transcripts = &self.held(e.location.as_deref())?.transcripts;
        Ok(transcripts
            .get(&e.feature)
            .or_else(|| transcripts.get(without_version(&e.feature).as_ref())))
    }

    /// What the shards of the chromosome directory `location` names hold, or those of
    /// every directory without a Location and for one that names no directory of the
    /// cache.
    fn held(&self, location: Option<&str>) -> Res<&Held> {
        let chromosome = location.map(|location| pytext::rpartition(location, ":").0);
        let held = chromosome.and_then(|chromosome| self.chromosomes.get_key_value(chromosome));
        let held = match held {
            Some((chromosome, held)) => held.get_or_init(|| self.read(Some(chromosome))),
            None => self.every.get_or_init(|| self.read(None)),
        };
        held.as_ref().map_err(Refusal::clone)
    }

    /// The coding ends of the entry's transcript: of the record of its Feature, on the
    /// sequence its Location names, whose span shares a base with the Location, an
    /// insertion's being the two bases it stands between. An entry without a Location is
    /// looked up at the Location its record gives its allele, and a Feature carrying a
    /// version the cache's stable id lacks without it. `None` for a transcript the cache
    /// does not hold, and where no record of the stable id holds the Location or more than
    /// one does: one id stands on several records of a sequence, and a row names one of
    /// them by where it lies. An entry without a position is answered where the id names
    /// one record. A shard that cannot be read is refused.
    pub fn coding_ends(&self, e: &Entry) -> Res<Option<&CodingEnds>> {
        let location = row_location(e);
        let held = self.held(location)?;
        let Some(records) = held
            .ends
            .get(&e.feature)
            .or_else(|| held.ends.get(without_version(&e.feature).as_ref()))
        else {
            return Ok(None);
        };
        let span = location.and_then(location_bases);
        let mut holding = records.iter().filter(|ends| match (span, ends.span) {
            (None, _) => true,
            (Some((first, last)), Some((start, end))) => first <= end && start <= last,
            (Some(_), None) => false,
        });
        Ok(match (holding.next(), holding.next()) {
            (Some(ends), None) => Some(ends),
            _ => None,
        })
    }
}
