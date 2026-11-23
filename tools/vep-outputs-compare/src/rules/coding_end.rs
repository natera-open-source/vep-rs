// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a comparison hands a pair rule that reads the coding ends of the pair's transcript:
//! the reading of the transcript's cache record, and the reference bases the one rule that
//! reads past a transcript's last base asks.

use super::cache_transcripts::CacheTranscripts;
use super::coding_ends::CodingEnds;
use super::reference_bases::ReferenceBases;
use super::row::transcript_feature_type;
use std::io;
use std::sync::Arc;

/// What a comparison read of the coding ends of one pair's transcript.
#[derive(Clone, Copy)]
pub enum CodingEndsRead<'a> {
    /// The comparison was given no cache: a rule that asks the coding ends takes nothing,
    /// and a rule that stood before the coding ends were read stands as it stood.
    NoCache,
    /// The comparison was given a cache, and the cache answers the pair's row with no one
    /// record: no rule that asks the coding ends takes the pair.
    NoRecord,
    /// The reading of the cache record the pair's row is answered with.
    Record(&'a CodingEnds),
}

/// What a pair rule of a coding end is asked with, beside the two rows of the pair.
#[derive(Clone, Copy)]
pub struct AtACodingEnd<'a> {
    pub ends: CodingEndsRead<'a>,
    /// The reference base at a genomic position of the pair's sequence, `None` past the
    /// sequence's ends; `None` where the comparison was given no reference sequence, or
    /// one that does not name the pair's sequence.
    pub reference_bases: Option<&'a dyn Fn(i64) -> Option<u8>>,
}

impl AtACodingEnd<'static> {
    /// A pair of a comparison that was given no cache and no reference sequence.
    pub const NO_CACHE: AtACodingEnd<'static> = AtACodingEnd {
        ends: CodingEndsRead::NoCache,
        reference_bases: None,
    };
}

/// The cache and the reference sequence a comparison was given, as its pair rules ask them.
#[derive(Clone, Copy, Default)]
pub struct CodingEndSources<'a> {
    /// The transcripts of the candidate's cache.
    pub transcripts: Option<&'a CacheTranscripts>,
    /// The reference sequence the candidate was run on.
    pub reference: Option<&'a ReferenceBases>,
}

/// What the sources hold for the pair of the transcript `feature` at `location`, kept
/// while a rule is asked with it.
pub struct PairSources<'a> {
    given: bool,
    record: Option<Arc<CodingEnds>>,
    reference: Option<&'a ReferenceBases>,
    sequence: &'a str,
}

impl<'a> CodingEndSources<'a> {
    /// Whether the comparison was given a cache.
    pub fn given(&self) -> bool {
        self.transcripts.is_some()
    }

    /// Reads what the sources hold for the pair of `feature` at `location`: the record of a
    /// Feature whose `feature_type` is a transcript's, and no record for any other. A cache
    /// that cannot be read is the failure.
    pub fn of_pair(
        &self,
        location: &'a str,
        feature: &str,
        feature_type: &str,
    ) -> io::Result<PairSources<'a>> {
        let record = match self.transcripts {
            Some(transcripts) if transcript_feature_type(feature_type) => {
                transcripts.coding_ends(location, feature)?
            }
            _ => None,
        };
        Ok(PairSources {
            given: self.given(),
            record,
            reference: self.reference,
            sequence: location
                .rsplit_once(':')
                .map_or("", |(sequence, _)| sequence),
        })
    }
}

impl PairSources<'static> {
    /// What a comparison that was given no cache and no reference sequence holds for a pair.
    #[allow(dead_code)]
    pub fn no_cache() -> PairSources<'static> {
        PairSources {
            given: false,
            record: None,
            reference: None,
            sequence: "",
        }
    }
}

impl PairSources<'_> {
    /// Asks `rule` with what was read for the pair.
    pub fn ask<T>(&self, rule: impl FnOnce(&AtACodingEnd<'_>) -> T) -> T {
        let ends = match (&self.record, self.given) {
            (Some(record), _) => CodingEndsRead::Record(record),
            (None, true) => CodingEndsRead::NoRecord,
            (None, false) => CodingEndsRead::NoCache,
        };
        // A file that does not name the pair's sequence hands a rule no bases: a base it
        // cannot answer for reads to a rule as one past the sequence's end.
        let reference = self
            .reference
            .filter(|reference| reference.holds(self.sequence));
        let base =
            |position: i64| reference.and_then(|reference| reference.base(self.sequence, position));
        let reference_bases: Option<&dyn Fn(i64) -> Option<u8>> = match reference {
            Some(_) => Some(&base),
            None => None,
        };
        rule(&AtACodingEnd {
            ends,
            reference_bases,
        })
    }
}
