// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The MaxEntScan plugin's class and the introns it reads.
//!
//! A substitution in an exon of at most five bases lies in the acceptor window of the intron
//! before it and the donor window of the intron after it; MaxEntScan.pm:270 scores the
//! window of whichever intron Ensembl VEP's interval tree returns first.

use super::cache::TranscriptShards;
use super::entry::{as_number, field_named, Entry};
use super::shared::{cache_integer, without_version};
use super::{Refusal, Res};
use crate::pytext::{self, Json};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::{Arc, Mutex};

pub const MAXENTSCAN_WINDOW_TEXT_FIELDS: [&str; 2] = ["MaxEntScan_ref_seq", "MaxEntScan_alt_seq"];

/// A transcript's strand and its introns in its own order, each as the last base of the
/// exon before it and the first base of the exon after it on the genome.
type Introns = (i128, Vec<(i128, i128)>);

/// The introns of a JSON cache's transcripts, by transcript stable id, each transcript's in
/// its own order with its strand. A chromosome's shards are read on the first look-up of an
/// entry whose record names its directory, as written or without a leading `chr`; a record
/// on a chromosome the cache holds no directory of lies in no window. A stable id a
/// chromosome's shards give two different intron lists is left out.
pub struct SpliceWindows {
    shards: TranscriptShards,
    introns: Mutex<HashMap<String, Arc<HashMap<String, Introns>>>>,
}

impl SpliceWindows {
    /// The introns of the cache at `json_cache`, refused where it holds no `transcripts/`
    /// directory.
    pub fn new(json_cache: &Path) -> Result<SpliceWindows, Refusal> {
        Ok(SpliceWindows {
            shards: TranscriptShards::new(json_cache)?,
            introns: Mutex::new(HashMap::new()),
        })
    }

    /// Reads a chromosome's shards on `threads` threads from the next chromosome read on.
    pub fn read_on(&self, threads: usize) {
        self.shards.read_on(threads);
    }

    fn read(&self, chromosome: &str) -> Res<HashMap<String, Introns>> {
        let mut introns: HashMap<String, Introns> = HashMap::new();
        let mut conflicting: HashSet<String> = HashSet::new();
        self.shards.of_chromosome(chromosome, |transcript| {
            let (Some(Json::Str(stable_id)), Some(Json::List(exons))) =
                (transcript.get("stable_id"), transcript.get("exons"))
            else {
                return Ok(());
            };
            let Some(strand @ (1 | -1)) = cache_integer(transcript.get("strand")) else {
                return Ok(());
            };
            if !exons.iter().all(|exon| matches!(exon, Json::Object(_))) {
                return Ok(());
            }
            let bounds: Option<Vec<(i128, i128)>> = exons
                .iter()
                .map(|exon| {
                    Some((
                        cache_integer(exon.get("start"))?,
                        cache_integer(exon.get("end"))?,
                    ))
                })
                .collect();
            let Some(mut spans) = bounds else {
                return Ok(());
            };
            spans.sort_unstable();
            let mut gaps: Vec<(i128, i128)> = spans
                .windows(2)
                .map(|pair| (pair[0].1, pair[1].0))
                .collect();
            if strand == -1 {
                gaps.reverse();
            }
            let read = (strand, gaps);
            match introns.get(stable_id) {
                Some(held) if *held != read => {
                    conflicting.insert(stable_id.clone());
                }
                Some(_) => {}
                None => {
                    introns.insert(stable_id.clone(), read);
                }
            }
            Ok(())
        })?;
        for stable_id in conflicting {
            introns.remove(&stable_id);
        }
        Ok(introns)
    }

    /// The chromosome and position of an entry's record where the record is one base: a
    /// Location that names one position. `None` for any other record, and for an entry
    /// without a Location, which no input line of the run places.
    pub fn place(&self, e: &Entry) -> Option<(String, i128)> {
        let (chromosome, _, span) = pytext::rpartition(e.location.as_deref()?, ":");
        if chromosome.is_empty() || !pytext::isdigit(span) {
            return None;
        }
        Some((chromosome.to_string(), pytext::int(span)?))
    }

    /// Whether the entry's record is one base that lies in the acceptor window of one intron
    /// of the entry's transcript and in the donor window of the next intron in the
    /// transcript's order, the windows as MaxEntScan.pm:282-291 places them: the donor's
    /// from three bases before an intron to its sixth base, the acceptor's from an intron's
    /// twentieth base from its end to three bases after it. A Feature carrying a version the
    /// cache's stable id lacks is looked up without it.
    pub fn holds(&self, e: &Entry) -> Res<bool> {
        let Some((chromosome, position)) = self.place(e) else {
            return Ok(false);
        };
        let bare = chromosome.strip_prefix("chr").unwrap_or(&chromosome);
        let Some(directory) = [chromosome.as_str(), bare]
            .into_iter()
            .find(|name| self.shards.chromosomes().contains(*name))
        else {
            return Ok(false);
        };
        let introns = {
            let mut read = self.introns.lock().expect("lock");
            match read.get(directory) {
                Some(introns) => Arc::clone(introns),
                None => {
                    let introns = Arc::new(self.read(directory)?);
                    read.insert(directory.to_string(), Arc::clone(&introns));
                    introns
                }
            }
        };
        let Some((strand, gaps)) = introns
            .get(&e.feature)
            .or_else(|| introns.get(without_version(&e.feature).as_ref()))
        else {
            return Ok(false);
        };
        // Whether the position lies `low` to `high` bases past the exon base `beside` an
        // intron; a distance outside `i128` is in no window.
        let within = |beside: i128, low: i128, high: i128| {
            position
                .checked_sub(beside)
                .is_some_and(|past| low <= past && past <= high)
        };
        Ok(gaps.windows(2).any(|pair| {
            let ((upstream_before, upstream_after), (downstream_before, downstream_after)) =
                (pair[0], pair[1]);
            // The windows are counted from the intron's own bases, one base inside the two
            // exon bases held here.
            if *strand == 1 {
                within(upstream_after, -20, 2) && within(downstream_before, -2, 6)
            } else {
                within(upstream_before, -2, 20) && within(downstream_after, -6, 2)
            }
        }))
    }
}

/// Whether `text` is `count` bases and nothing else.
fn bases(text: &str, count: usize) -> bool {
    text.len() == count && text.bytes().all(|b| b"ACGT".contains(&b))
}

/// The shape of `maxentscan_overlapping_windows`: on a row whose record is one base in the
/// acceptor window of one intron of the row's transcript and the donor window of the next
/// (`SpliceWindows::holds`), a `MaxEntScan_` field both engines fill with numbers, and a
/// window text where the reference writes a donor's nine bases and the engine an acceptor's
/// twenty-three. Without the cache nothing is accepted.
pub fn overlapping_windows(
    e: &Entry,
    name: &str,
    ev: &str,
    av: &str,
    splice_windows: Option<&SpliceWindows>,
) -> Res<bool> {
    let Some(windows) = splice_windows else {
        return Ok(false);
    };
    if !windows.holds(e)? {
        return Ok(false);
    }
    if field_named(name, &MAXENTSCAN_WINDOW_TEXT_FIELDS) {
        return Ok(bases(ev, 9) && bases(av, 23));
    }
    Ok(as_number(ev).is_some() && as_number(av).is_some())
}
