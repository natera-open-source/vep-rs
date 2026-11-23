// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A transcript model read from the cache.

use crate::consts::INTEGER_RE_PATTERN;
use crate::fields::python_type;
use crate::pyfmt;
use crate::pyjson::Value;
use crate::rules::whole_text;
use std::sync::LazyLock;

/// The keys of one pair of a cache transcript's `mapper`: a run of cDNA coordinates, the
/// genomic bases it maps to and their orientation.
pub const MAPPER_PAIR_KEYS: &[&str] = &["from_start", "from_end", "to_start", "to_end", "ori"];

static WHOLE_INTEGER_RE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text(INTEGER_RE_PATTERN));

/// `int(text)` of a text `int` reads, where the number fits 64 bits. A text of more digits
/// than the interpreter converts ends the run as its `ValueError` does.
pub(crate) fn int_of_text(text: &str) -> Option<i64> {
    match pyfmt::int_checked(text) {
        Ok(number) => number.parse().ok(),
        Err(error) => crate::uncaught!("{error}"),
    }
}

/// An integer of a cache transcript, a number in a native cache and the dump's string in a
/// converted one; `None` for any other value, and for a number 64 bits do not hold.
pub fn cache_integer(value: &Value) -> Option<i64> {
    match value {
        Value::Int(number) => Some(*number),
        Value::Str(text) if WHOLE_INTEGER_RE.is_match(text) => int_of_text(text),
        _ => None,
    }
}

fn wide(value: i64) -> i128 {
    i128::from(value)
}

fn narrow(value: i128) -> Option<i64> {
    i64::try_from(value).ok()
}

/// A cache transcript's coordinates: its `strand`, its `exons` in transcript order as
/// (start, end), its `cdna_coding_start` and `cdna_coding_end`, the count of `N` its
/// `translateable_seq` opens with (the start exon's phase, which the CDS coordinates count
/// in) and its `mapper` pairs as (cDNA start, cDNA end, genomic start, genomic end,
/// orientation).
///
/// A coordinate a method works out that 64 bits do not hold is read as no coordinate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranscriptModel {
    pub strand: i64,
    pub exons: Vec<(i64, i64)>,
    pub coding_start: i64,
    pub coding_end: i64,
    pub padding: i64,
    pub pairs: Vec<(i64, i64, i64, i64, i64)>,
}

impl TranscriptModel {
    /// The model of a cache transcript object, `None` where the cache states a part of
    /// it in no form this reads.
    pub fn from_cache(transcript: &Value, padding: i64) -> Option<Self> {
        let cached = match transcript.get("variation_effect_feature_cache") {
            Some(cached) => cached,
            None => crate::uncaught!("'variation_effect_feature_cache'"),
        };
        let mapper = match cached.as_dict() {
            Some(members) => members.get("mapper"),
            None => crate::uncaught!("'{}' object has no attribute 'get'", python_type(cached)),
        };
        let pairs = mapper.and_then(|mapper| mapper.get("pairs"));
        let exons = transcript.get("exons");
        let integer = |item: &Value, key: &str| item.get(key).and_then(cache_integer);
        let strand = integer(transcript, "strand");
        let coding_start = integer(transcript, "cdna_coding_start");
        let coding_end = integer(transcript, "cdna_coding_end");
        let (Some(strand @ (1 | -1)), Some(coding_start), Some(coding_end)) =
            (strand, coding_start, coding_end)
        else {
            return None;
        };
        let pairs = pairs.and_then(Value::as_list)?;
        let exons = exons.and_then(Value::as_list)?;
        if pairs.is_empty() || exons.is_empty() {
            return None;
        }
        if !pairs
            .iter()
            .chain(exons)
            .all(|item| item.as_dict().is_some())
        {
            return None;
        }
        // Every number is read before one is judged: a number `int` refuses ends the run
        // whatever stands before it.
        let read_pairs: Vec<Vec<Option<i64>>> = pairs
            .iter()
            .map(|pair| {
                MAPPER_PAIR_KEYS
                    .iter()
                    .map(|key| integer(pair, key))
                    .collect()
            })
            .collect();
        let read_exons: Vec<(Option<i64>, Option<i64>)> = exons
            .iter()
            .map(|exon| (integer(exon, "start"), integer(exon, "end")))
            .collect();
        let pairs = read_pairs
            .iter()
            .map(|pair| Some((pair[0]?, pair[1]?, pair[2]?, pair[3]?, pair[4]?)))
            .collect::<Option<Vec<_>>>()?;
        let mut exons = read_exons
            .iter()
            .map(|(start, end)| Some(((*start)?, (*end)?)))
            .collect::<Option<Vec<_>>>()?;
        exons.sort_unstable();
        if strand == -1 {
            exons.reverse();
        }
        Some(TranscriptModel {
            strand,
            exons,
            coding_start,
            coding_end,
            padding,
            pairs,
        })
    }

    fn cdna_wide(&self, genomic: i128) -> Option<i128> {
        for &(cdna_start, _, genomic_start, genomic_end, orientation) in &self.pairs {
            if wide(genomic_start) <= genomic && genomic <= wide(genomic_end) {
                let along = if orientation == -1 {
                    wide(genomic_end) - genomic
                } else {
                    genomic - wide(genomic_start)
                };
                return Some(wide(cdna_start) + along);
            }
        }
        None
    }

    /// The cDNA coordinate of a genomic base, `None` for a base outside the mapper pairs.
    pub fn cdna(&self, genomic: i64) -> Option<i64> {
        self.cdna_wide(wide(genomic)).and_then(narrow)
    }

    fn genomic_wide(&self, cds: i128) -> Option<i128> {
        let cdna = cds - wide(self.padding) + wide(self.coding_start) - 1;
        if !(wide(self.coding_start) <= cdna && cdna <= wide(self.coding_end)) {
            return None;
        }
        for &(cdna_start, cdna_end, genomic_start, genomic_end, orientation) in &self.pairs {
            if wide(cdna_start) <= cdna && cdna <= wide(cdna_end) {
                return Some(if orientation == -1 {
                    wide(genomic_end) - (cdna - wide(cdna_start))
                } else {
                    wide(genomic_start) + cdna - wide(cdna_start)
                });
            }
        }
        None
    }

    /// The genomic base of a CDS coordinate, `None` for a coordinate outside the coding
    /// bounds or the mapper pairs.
    pub fn genomic(&self, cds: i64) -> Option<i64> {
        self.genomic_wide(wide(cds)).and_then(narrow)
    }

    /// `cds_pair` of a span whose ends are sums of 64-bit numbers.
    pub(crate) fn cds_pair_wide(&self, start: i128, end: i128) -> (Option<i64>, Option<i64>) {
        let to_cds = wide(self.padding) - wide(self.coding_start) + 1;
        let (coding_start, coding_end) = (wide(self.coding_start), wide(self.coding_end));
        if start == end + 1 {
            let (five, three) = if self.strand == 1 {
                (end, start)
            } else {
                (start, end)
            };
            let (before, after) = match (self.cdna_wide(five), self.cdna_wide(three)) {
                (None, None) => return (None, None),
                (None, Some(after)) => (after - 1, after),
                (Some(before), None) => (before, before + 1),
                (Some(before), Some(after)) => (before, after),
            };
            return (
                (coding_start <= before && before < coding_end)
                    .then(|| narrow(before + 1 + to_cds))
                    .flatten(),
                (coding_start < after && after <= coding_end)
                    .then(|| narrow(after - 1 + to_cds))
                    .flatten(),
            );
        }
        let (five, three) = if self.strand == 1 {
            (start, end)
        } else {
            (end, start)
        };
        let coding = |base: i128| {
            self.cdna_wide(base)
                .filter(|cdna| coding_start <= *cdna && *cdna <= coding_end)
                .and_then(|cdna| narrow(cdna + to_cds))
        };
        (coding(five), coding(three))
    }

    /// The `cds_start` and `cds_end` of the genomic span `start` to `end`
    /// (BaseTranscriptVariation.pm:252-264 over `TranscriptMapper::genomic2cds`), `None`
    /// for an end that is no coding base: the CDS coordinates of the span's 5' and 3' bases;
    /// for an insertion (`start` is `end` + 1) `map_insert`'s result
    /// (Mapper.pm:504-558), the base after the 5' flank's coordinate and the base before the
    /// 3' flank's, a flank outside the mapper pairs read as the base beside the other.
    pub fn cds_pair(&self, start: i64, end: i64) -> (Option<i64>, Option<i64>) {
        self.cds_pair_wide(wide(start), wide(end))
    }

    /// The one genomic span (start, end; an insertion's start one past its end) whose
    /// `cds_pair` is `cds_start` and `cds_end`, or `None`: an end the mapper pairs do
    /// not hold, or an insertion's pair (`cds_start` is `cds_end` + 1) between two bases
    /// the genome does not keep adjacent, which the point after the one and the point before
    /// the other both give.
    pub fn genomic_span(&self, cds_start: i64, cds_end: i64) -> Option<(i64, i64)> {
        let first = self.genomic_wide(wide(cds_start.min(cds_end)))?;
        let last = self.genomic_wide(wide(cds_start.max(cds_end)))?;
        if cds_start > cds_end && (first - last).abs() != 1 {
            return None;
        }
        let (low, high) = (narrow(first.min(last))?, narrow(first.max(last))?);
        Some(if cds_start > cds_end {
            (high, low)
        } else {
            (low, high)
        })
    }

    /// `variant_exon_check` (NMD.pm:132-176) on a record's genomic end: it lies in the
    /// last exon, or in the second to last within 51 bases of that exon's 3' end, its end on
    /// the forward strand and its start on the reverse.
    pub fn ends_where_nmd_escapes(&self, genomic_end: i64) -> bool {
        let Some(&(start, end)) = self.exons.last() else {
            crate::uncaught!("tuple index out of range")
        };
        if start <= genomic_end && genomic_end <= end {
            return true;
        }
        if self.exons.len() < 2 {
            return false;
        }
        let (start, end) = self.exons[self.exons.len() - 2];
        let (start, end, genomic_end) = (wide(start), wide(end), wide(genomic_end));
        if self.strand == -1 {
            start <= genomic_end && genomic_end <= start + 51
        } else {
            end - 51 <= genomic_end && genomic_end <= end
        }
    }
}
