// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The exon spans of a JSON cache's transcripts and the place of a base along them.

use super::cache::TranscriptShards;
use super::shared::{cache_integer, without_version};
use super::Res;
use crate::pytext::Json;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::{Arc, Mutex};

/// The largest exon bound read: the distances `TranscriptExons::exon_base` takes between
/// such a bound and a position a rule reads from a row stay inside `i128`.
const LARGEST_BOUND: i128 = 1 << 62;

/// A transcript's exons on its chromosome and its strand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranscriptExons {
    /// Whether the transcript lies on the reverse strand.
    pub reverse: bool,
    /// The first and last base of each exon, the exons in ascending order.
    pub spans: Vec<(i128, i128)>,
}

impl TranscriptExons {
    /// The exon base that the chromosome's base `position` is numbered from, as its place
    /// along the exons counted from the transcript's first base, and the intron part from
    /// it: 0 in an exon; in an intron the distance past the last base of the exon before it
    /// on the transcript's strand or, where the exon after it is nearer, that distance below
    /// zero from that exon's first base (TranscriptVariationAllele.pm:2798-2829). A base
    /// outside the exons' span has place 0.
    pub fn exon_base(&self, position: i128) -> (i128, i128) {
        let mut place: i128 = 0;
        let mut before: Option<(i128, i128)> = None;
        for index in 0..self.spans.len() {
            let (start, end) = self.spans[if self.reverse {
                self.spans.len() - 1 - index
            } else {
                index
            }];
            if start <= position && position <= end {
                let along = if self.reverse {
                    end - position
                } else {
                    position - start
                };
                return (place + along + 1, 0);
            }
            if let Some((before_start, before_end)) = before {
                let (past, ahead) = if self.reverse {
                    (before_start - position, position - end)
                } else {
                    (position - before_end, start - position)
                };
                if past > 0 && ahead > 0 {
                    return if past <= ahead {
                        (place, past)
                    } else {
                        (place + 1, -ahead)
                    };
                }
            }
            place += end - start + 1;
            before = Some((start, end));
        }
        (0, 0)
    }
}

/// A cache transcript object's exons under its stable id; `None` for an object without a
/// stable id, without a `strand` of 1 or -1 or without exons, and for one with an exon
/// that is no object, states no whole `start` and `end` from 1 to `LARGEST_BOUND`, ends
/// before it starts or shares a base with another.
fn read_transcript(transcript: &Json) -> Option<(String, TranscriptExons)> {
    let (Some(Json::Str(stable_id)), Some(Json::List(exons))) =
        (transcript.get("stable_id"), transcript.get("exons"))
    else {
        return None;
    };
    let reverse = match cache_integer(transcript.get("strand"))? {
        1 => false,
        -1 => true,
        _ => return None,
    };
    let mut spans: Vec<(i128, i128)> = exons
        .iter()
        .map(|exon| {
            let start = cache_integer(exon.get("start"))?;
            let end = cache_integer(exon.get("end"))?;
            (1 <= start && start <= end && end <= LARGEST_BOUND).then_some((start, end))
        })
        .collect::<Option<_>>()?;
    spans.sort_unstable();
    if spans.is_empty() || spans.windows(2).any(|pair| pair[1].0 <= pair[0].1) {
        return None;
    }
    Some((stable_id.clone(), TranscriptExons { reverse, spans }))
}

/// The exon spans of a JSON cache's transcripts, by transcript stable id. A chromosome's
/// shards are read on the first look-up that names its directory. A transcript is held
/// where `read_transcript` reads it; a stable id a chromosome's shards give two different
/// readings is left out.
pub struct ExonSpans {
    shards: TranscriptShards,
    read: Mutex<HashMap<String, Arc<HashMap<String, TranscriptExons>>>>,
}

impl ExonSpans {
    /// The exon spans of the cache at `json_cache`, refused where it holds no
    /// `transcripts/` directory. Opening the cache reads no shard.
    pub fn new(json_cache: &Path) -> Res<ExonSpans> {
        Ok(ExonSpans {
            shards: TranscriptShards::new(json_cache)?,
            read: Mutex::new(HashMap::new()),
        })
    }

    /// Reads a chromosome's shards on `threads` threads from the next chromosome read on.
    pub fn read_on(&self, threads: usize) {
        self.shards.read_on(threads);
    }

    /// The cache's directory of the chromosome `chromosome`: its name as written, without
    /// or with a leading `chr`, or the mitochondrion's other name, the first of them the
    /// cache holds.
    fn directory(&self, chromosome: &str) -> Option<String> {
        if chromosome.is_empty() {
            return None;
        }
        let prefixed = chromosome
            .get(..3)
            .is_some_and(|prefix| prefix.eq_ignore_ascii_case("chr"));
        let bare = if prefixed {
            &chromosome[3..]
        } else {
            chromosome
        };
        let other = match bare {
            "M" => "MT",
            "MT" => "M",
            _ => bare,
        };
        [
            chromosome.to_string(),
            bare.to_string(),
            format!("chr{bare}"),
            other.to_string(),
            format!("chr{other}"),
        ]
        .into_iter()
        .find(|name| self.shards.chromosomes().contains(name))
    }

    fn read(&self, directory: &str) -> Res<HashMap<String, TranscriptExons>> {
        let mut transcripts: HashMap<String, TranscriptExons> = HashMap::new();
        let mut conflicting: HashSet<String> = HashSet::new();
        self.shards.read_of_chromosome(
            directory,
            |transcript| read_transcript(&transcript),
            |(stable_id, read)| {
                match transcripts.get(&stable_id) {
                    Some(held) if *held != read => {
                        conflicting.insert(stable_id);
                    }
                    Some(_) => {}
                    None => {
                        transcripts.insert(stable_id, read);
                    }
                }
                Ok(())
            },
        )?;
        for stable_id in conflicting {
            transcripts.remove(&stable_id);
        }
        Ok(transcripts)
    }

    /// The exons of the transcript `feature` on the chromosome `chromosome`; `None` for a
    /// chromosome the cache holds no directory of and for a transcript it does not hold
    /// there. A Feature carrying a version the cache's stable id lacks is looked up without
    /// it. A shard that cannot be read is refused.
    pub fn transcript(&self, chromosome: &str, feature: &str) -> Res<Option<TranscriptExons>> {
        let Some(directory) = self.directory(chromosome) else {
            return Ok(None);
        };
        let transcripts = {
            let mut read = self.read.lock().expect("lock");
            match read.get(&directory) {
                Some(transcripts) => Arc::clone(transcripts),
                None => {
                    let transcripts = Arc::new(self.read(&directory)?);
                    read.insert(directory, Arc::clone(&transcripts));
                    transcripts
                }
            }
        };
        Ok(transcripts
            .get(feature)
            .or_else(|| transcripts.get(without_version(feature).as_ref()))
            .cloned())
    }
}
