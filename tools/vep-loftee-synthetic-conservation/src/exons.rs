// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a run reads of one sequence of the cache: the stable id and the exon count of
//! each transcript record the cache reader gives, whatever its biotype, and the bases
//! that lie in an exon, each once. The reader gives a record filed in several region
//! files once, and a stable id under two `dbID`s as two records.

use std::sync::atomic::{AtomicBool, Ordering};

use vep_seeded_variants::cache::{Cache, CacheError, Span};

/// A transcript record that has an exon.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// The transcript's stable id.
    pub stable_id: String,
    /// How many exons the record has: its exons are numbered from 1 to this count.
    pub exons: u64,
}

/// The transcript records and the exon bases of one sequence.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SequenceExons {
    /// The records with an exon, in the reader's order.
    pub records: Vec<Record>,
    /// How many records the reader gives no exon for: a record that places an exon
    /// outside its own span, or one stored across the origin of a circular sequence.
    pub without_exons: u64,
    /// The stable id of the first such record.
    pub first_without_exons: Option<String>,
    /// The bases that lie in an exon of a record: runs in ascending order that neither
    /// overlap nor touch.
    pub spans: Vec<Span>,
}

impl SequenceExons {
    /// Reads the transcripts of `sequence`. A raised `stopped` ends the read early with
    /// `None`.
    pub fn read(
        cache: &Cache,
        sequence: &str,
        stopped: &AtomicBool,
    ) -> Result<Option<SequenceExons>, CacheError> {
        let mut read = SequenceExons::default();
        let mut exons: Vec<Span> = Vec::new();
        for transcript in cache.transcripts(sequence)? {
            if stopped.load(Ordering::Relaxed) {
                return Ok(None);
            }
            let transcript = transcript?;
            if transcript.exons.is_empty() {
                read.without_exons += 1;
                read.first_without_exons.get_or_insert(transcript.stable_id);
                continue;
            }
            read.records.push(Record {
                exons: transcript.exons.len() as u64,
                stable_id: transcript.stable_id,
            });
            exons.extend(transcript.exons);
        }
        read.spans = union(exons);
        Ok(Some(read))
    }

    /// How many transcript records the reader gives, with an exon or without.
    pub fn transcripts(&self) -> u64 {
        self.records.len() as u64 + self.without_exons
    }

    /// How many exons the records with an exon have.
    pub fn exons(&self) -> u64 {
        self.records.iter().map(|record| record.exons).sum()
    }

    /// How many bases lie in an exon, each counted once.
    pub fn exon_bases(&self) -> u64 {
        self.spans.iter().map(Span::length).sum()
    }
}

/// The bases of `spans`, each once: runs in ascending order, two runs that overlap or
/// touch joined into one.
pub fn union(mut spans: Vec<Span>) -> Vec<Span> {
    spans.sort_unstable();
    let mut joined: Vec<Span> = Vec::new();
    for span in spans {
        match joined.last_mut() {
            Some(last) if span.start <= last.end.saturating_add(1) => {
                last.end = last.end.max(span.end);
            }
            _ => joined.push(span),
        }
    }
    joined
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spans(runs: &[(u64, u64)]) -> Vec<Span> {
        runs.iter()
            .map(|(start, end)| Span {
                start: *start,
                end: *end,
            })
            .collect()
    }

    #[test]
    fn runs_that_overlap_or_touch_are_joined() {
        let joined = union(spans(&[
            (50, 60),
            (10, 20),
            (15, 30),
            (31, 40),
            (12, 13),
            (60, 60),
            (62, 70),
        ]));
        assert_eq!(joined, spans(&[(10, 40), (50, 60), (62, 70)]));
        assert_eq!(union(Vec::new()), Vec::new());
        assert_eq!(union(spans(&[(5, 5), (5, 5)])), spans(&[(5, 5)]));
    }
}
