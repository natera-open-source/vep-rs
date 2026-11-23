// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Transcript edges: an allele placed on a transcript's exons.

use crate::argparse::{path_name, path_text};
use crate::cache;
use crate::entry::Entry;
use crate::error::{Error, Result};
use crate::input::{is_directory, read_text};
use crate::pyfmt;
use crate::pyjson::{self, Value};
use crate::rules::clipped::{length, parted};
use crate::rules::nearest_i64;
use crate::shared::integer_of;
use indexmap::IndexMap;
use std::cell::RefCell;
use std::path::Path;
use std::path::PathBuf;

/// A transcript's exon spans in genomic order and whether it lies on the reverse strand.
pub type ExonSpans = (Vec<(i64, i64)>, bool);

/// The transcripts of one shard by stable id, once per locus the shard files the id at.
pub type ShardExons = IndexMap<String, Vec<ExonSpans>>;

/// An entry's allele beside its transcript: the allele's first and last genomic position,
/// the transcript's exon spans in genomic order and whether it lies on the reverse strand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlacedAllele {
    pub start: i64,
    pub end: i64,
    pub exons: Vec<(i64, i64)>,
    pub reverse: bool,
}

impl PlacedAllele {
    /// The bases of the allele before the transcript's first base and past its last base,
    /// on the transcript's strand.
    pub fn cuts(&self) -> (i64, i64) {
        let first = self.exons.iter().map(|(start, _)| *start).min();
        let last = self.exons.iter().map(|(_, end)| *end).max();
        let (Some(first), Some(last)) = (first, last) else {
            crate::uncaught!("min() iterable argument is empty")
        };
        let low = nearest_i64((i128::from(first) - i128::from(self.start)).max(0));
        let high = nearest_i64((i128::from(self.end) - i128::from(last)).max(0));
        if self.reverse {
            (high, low)
        } else {
            (low, high)
        }
    }

    /// The exon base that the allele's base `offset` bases along the transcript's strand
    /// is numbered from, as its place along the exons counted from the transcript's first
    /// base, and the intron offset from it: 0 in an exon; in an intron the distance past the
    /// last base of the exon before it or, where the exon after it is nearer, that distance
    /// below zero from that exon's first base (TranscriptVariationAllele.pm:2798-2829). A
    /// base outside the exons' span has place 0.
    pub fn exon_base(&self, offset: i64) -> (i64, i64) {
        let position = if self.reverse {
            i128::from(self.end) - i128::from(offset)
        } else {
            i128::from(self.start) + i128::from(offset)
        };
        let mut place: i128 = 0;
        let mut before: Option<(i128, i128)> = None;
        for index in 0..self.exons.len() {
            let (start, end) = self.exons[if self.reverse {
                self.exons.len() - 1 - index
            } else {
                index
            }];
            let (start, end) = (i128::from(start), i128::from(end));
            if start <= position && position <= end {
                let along = if self.reverse {
                    end - position
                } else {
                    position - start
                };
                return (nearest_i64(place + along + 1), 0);
            }
            if let Some((before_start, before_end)) = before {
                let (past, ahead) = if self.reverse {
                    (before_start - position, position - end)
                } else {
                    (position - before_end, start - position)
                };
                if past > 0 && ahead > 0 {
                    return if past <= ahead {
                        (nearest_i64(place), nearest_i64(past))
                    } else {
                        (nearest_i64(place + 1), nearest_i64(-ahead))
                    };
                }
            }
            place += end - start + 1;
            before = Some((start, end));
        }
        (0, 0)
    }

    /// The intron offset of the allele's base `offset` bases along the transcript's
    /// strand (`exon_base`).
    pub fn intron_offset(&self, offset: i64) -> i64 {
        self.exon_base(offset).1
    }
}

/// Where an entry's allele lies on its transcript, from the exon spans of the candidate
/// JSON cache (`--json-cache`): spans and lengths alone, no sequence. A text entry's span is
/// its Location; a VCF or JSON entry's is its input record's POS and REF length (`inputs`),
/// less the shared first base where the parser drops it. A shard of the cache is read when an
/// entry first needs it.
#[derive(Debug)]
pub struct TranscriptEdges {
    pub cache: PathBuf,
    /// The input records, each cut after its ALT column: no later column is read.
    pub inputs: Option<Vec<String>>,
    pub _chromosomes: RefCell<IndexMap<String, Vec<(i64, PathBuf)>>>,
    pub _shards: RefCell<IndexMap<PathBuf, ShardExons>>,
}

impl TranscriptEdges {
    /// `inputs` is `None` where the caller has no input records.
    pub fn new(cache: &Path, inputs: Option<&[String]>) -> Result<Self> {
        if !is_directory(&cache::joined(cache, "transcripts"))? {
            return Err(Error::comparator(format!(
                "--json-cache: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        Ok(TranscriptEdges {
            cache: cache.to_path_buf(),
            inputs: inputs.map(|records| {
                records
                    .iter()
                    .map(|record| through_the_alternates(record).to_string())
                    .collect()
            }),
            _chromosomes: RefCell::new(IndexMap::new()),
            _shards: RefCell::new(IndexMap::new()),
        })
    }

    /// The chromosome and the first and last genomic position of the entry's allele. A
    /// record of several alternates, one of another length than REF, whose alleles all start
    /// with one base is read without that base, one position on (Parser/VCF.pm:301-321); an
    /// allele the record does not write so, as one --minimal trims, has no span here.
    pub fn allele_span(&self, e: &Entry) -> Result<Option<(String, i64, i64)>> {
        if let Some(location) = e.location.as_deref() {
            let (chrom, span) = location.rsplit_once(':').unwrap_or(("", location));
            let (start, end) = parted(span, '-');
            if !pyfmt::is_digit(start) || (!end.is_empty() && !pyfmt::is_digit(end)) {
                return Ok(None);
            }
            let first = integer_of(start, "Location position")?;
            let last = integer_of(
                if end.is_empty() { start } else { end },
                "Location position",
            )?;
            return Ok(Some((chrom.to_string(), first, last)));
        }
        let (Some(inputs), Some(record)) = (self.inputs.as_deref(), e.record) else {
            return Ok(None);
        };
        let count = inputs.len() as i128;
        if i128::from(record) >= count {
            return Ok(None);
        }
        // A record below zero counts from the last line, as a list's index does.
        let at = if record < 0 {
            i128::from(record) + count
        } else {
            i128::from(record)
        };
        if at < 0 {
            crate::uncaught!("list index out of range");
        }
        let cells: Vec<&str> = inputs[at as usize].split('\t').collect();
        if cells.len() < 5 || !pyfmt::is_digit(cells[1]) {
            return Ok(None);
        }
        let reference = pyfmt::upper(cells[3]);
        let alternate_cell = pyfmt::upper(cells[4]);
        let alternates: Vec<&str> = alternate_cell.split(',').collect();
        let chopped = alternates.len() > 1
            && alternates
                .iter()
                .any(|alternate| length(alternate) != length(&reference))
            && {
                let mut first_bases = std::iter::once(reference.as_str())
                    .chain(alternates.iter().copied())
                    .filter(|allele| !allele.contains('*'))
                    .map(|allele| pyfmt::prefix(allele, 1));
                first_bases
                    .next()
                    .is_some_and(|first| first_bases.all(|other| other == first))
            };
        let written = |alternate: &&str| {
            e.allele
                == if chopped && !alternate.contains('*') {
                    pyfmt::slice(alternate, Some(1), None)
                } else {
                    alternate
                }
        };
        if !alternates.iter().any(written)
            || length(&e.allele) != length(&reference) - i128::from(chopped)
        {
            return Ok(None);
        }
        let start = i128::from(integer_of(cells[1], "input record POS")?) + i128::from(chopped);
        Ok(Some((
            cells[0].to_string(),
            nearest_i64(start),
            nearest_i64(start + length(&e.allele) - 1),
        )))
    }

    /// Lists the shards of a chromosome on the first call that names it.
    fn list_chromosome(&self, chrom: &str) -> Result<()> {
        if self._chromosomes.borrow().contains_key(chrom) {
            return Ok(());
        }
        let bare = if pyfmt::lower(pyfmt::prefix(chrom, 3)) == "chr" {
            pyfmt::slice(chrom, Some(3), None)
        } else {
            chrom
        };
        let other = match bare {
            "M" => "MT",
            "MT" => "M",
            _ => bare,
        };
        let mut names: Vec<String> = Vec::new();
        for name in [
            chrom.to_string(),
            bare.to_string(),
            format!("chr{bare}"),
            other.to_string(),
            format!("chr{other}"),
        ] {
            if !name.is_empty() && !name.contains('/') && name != ".." && !names.contains(&name) {
                names.push(name);
            }
        }
        let root = cache::joined(&self.cache, "transcripts");
        let mut directory = None;
        for name in &names {
            if is_directory(&cache::joined(&root, name))? {
                directory = Some(name);
                break;
            }
        }
        let mut shards: Vec<(i64, PathBuf)> = Vec::new();
        if let Some(name) = directory {
            for shard in cache::shards(&root, Some(name))? {
                let file_name = path_name(&shard);
                let (low, high) = parted(parted(&file_name, '.').0, '-');
                if pyfmt::is_digit(low) && pyfmt::is_digit(high) {
                    let what =
                        format!("{}: the range start in the shard's name", path_text(&shard));
                    shards.push((integer_of(low, &what)?, shard));
                }
            }
        }
        shards.sort_by(|a, b| b.cmp(a));
        self._chromosomes
            .borrow_mut()
            .insert(chrom.to_string(), shards);
        Ok(())
    }

    /// The shards of a chromosome, each with the first position of the range its name
    /// gives, the highest first. The directory is the chromosome's name as written, without
    /// or with `chr`, or the mitochondrion's other name.
    pub fn chromosome_shards(&self, chrom: &str) -> Result<Vec<(i64, PathBuf)>> {
        self.list_chromosome(chrom)?;
        Ok(self
            ._chromosomes
            .borrow()
            .get(chrom)
            .cloned()
            .unwrap_or_default())
    }

    /// The exon spans and strand of the transcript `feature` where its exons meet
    /// `start` to `end`. A converted cache files a transcript in every shard its gene
    /// overlaps and a built one in the shard of its first base alone, so the shards are read
    /// from the allele's own toward the chromosome's start until one files the transcript. The
    /// id looked up is the Feature as written, then without its version.
    pub fn transcript(
        &self,
        chrom: &str,
        start: i64,
        end: i64,
        feature: &str,
    ) -> Result<Option<ExonSpans>> {
        self.list_chromosome(chrom)?;
        let chromosomes = self._chromosomes.borrow();
        let unversioned = parted(feature, '.').0;
        let stable_ids = if unversioned == feature {
            &[feature][..]
        } else {
            &[feature, unversioned][..]
        };
        for (low, shard) in chromosomes.get(chrom).into_iter().flatten() {
            if *low > end {
                continue;
            }
            if !self._shards.borrow().contains_key(shard) {
                let read = Self::shard_exons(shard)?;
                self._shards.borrow_mut().insert(shard.clone(), read);
            }
            let shards = self._shards.borrow();
            let Some(filed) = shards.get(shard) else {
                continue;
            };
            for stable_id in stable_ids {
                for (exons, reverse) in filed.get(*stable_id).into_iter().flatten() {
                    let first = exons.iter().map(|(first, _)| *first).min();
                    let last = exons.iter().map(|(_, last)| *last).max();
                    if first.is_some_and(|first| first <= end)
                        && last.is_some_and(|last| start <= last)
                    {
                        return Ok(Some((exons.clone(), *reverse)));
                    }
                }
            }
        }
        Ok(None)
    }

    /// Each transcript of a shard by stable id, once per locus the shard files the id at:
    /// its exon spans in genomic order and whether its strand is -1. A transcript without an
    /// id or without exons is left out. A shard that is not JSON or is nested deeper than the
    /// interpreter reads, and a transcript whose id is no string, whose exons are no array of
    /// objects or whose exon bound is no integer, end the run with the comparator's error.
    pub fn shard_exons(shard: &Path) -> Result<ShardExons> {
        let named = path_text(shard);
        let parsed = match pyjson::loads(&read_text(shard)?) {
            Ok(parsed) => parsed,
            Err(error) if error.is_recursion_error() => {
                return Err(Error::comparator(format!(
                    "{named} is nested too deep to read"
                )))
            }
            Err(error) if error.is_value_error() => {
                return Err(Error::comparator(format!("{named} is not JSON: {error}")))
            }
            Err(error) => crate::uncaught!("{named}: {error}"),
        };
        let mut transcripts = ShardExons::new();
        for transcript in cache::transcripts_of(shard, parsed)? {
            let filled = |name: &str| transcript.get(name).filter(|value| pyjson::truthy(value));
            let (Some(stable_id), Some(exons)) = (filled("stable_id"), filled("exons")) else {
                continue;
            };
            let (Some(stable_id), Some(exons)) = (stable_id.as_str(), exons.as_list()) else {
                return Err(malformed_transcript(&named));
            };
            if !exons
                .iter()
                .all(|exon| matches!(exon, Value::Null | Value::Dict(_)))
            {
                return Err(malformed_transcript(&named));
            }
            let shown = pyfmt::repr_str(pyfmt::prefix(stable_id, 40));
            let bounds: Option<Vec<String>> = exons
                .iter()
                .filter(|exon| exon.as_dict().is_some())
                .flat_map(|exon| [exon.get("start"), exon.get("end")])
                .map(bound_digits)
                .collect();
            let Some(bounds) = bounds else {
                return Err(Error::comparator(format!(
                    "{named}: transcript {shown} has an exon without an integer start and end"
                )));
            };
            let bound_of = format!("{named}: an exon bound of transcript {shown},");
            let mut spans: Vec<(i64, i64)> = Vec::with_capacity(bounds.len() / 2);
            for pair in bounds.chunks_exact(2) {
                let first = integer_of(&pair[0], &bound_of)?;
                let second = integer_of(&pair[1], &bound_of)?;
                spans.push((first.min(second), first.max(second)));
            }
            spans.sort_unstable();
            if !spans.is_empty() {
                let reverse = transcript.get("strand").is_some_and(|strand| {
                    pyjson::py_eq(strand, &Value::Int(-1)) || strand.as_str() == Some("-1")
                });
                transcripts
                    .entry(stable_id.to_string())
                    .or_default()
                    .push((spans, reverse));
            }
        }
        Ok(transcripts)
    }

    /// The entry's allele beside its transcript; `None` where the allele's span or the
    /// transcript's exons cannot be read, the allele is not as long as its span, or the two do
    /// not meet.
    pub fn place(&self, e: &Entry) -> Result<Option<PlacedAllele>> {
        let Some((chrom, start, end)) = self.allele_span(e)? else {
            return Ok(None);
        };
        if i128::from(end) - i128::from(start) + 1 != length(&e.allele) {
            return Ok(None);
        }
        let Some((exons, reverse)) = self.transcript(&chrom, start, end, &e.feature)? else {
            return Ok(None);
        };
        Ok(Some(PlacedAllele {
            start,
            end,
            exons,
            reverse,
        }))
    }
}

/// An input record up to the end of its fifth column, ALT; a shorter record whole.
fn through_the_alternates(record: &str) -> &str {
    record
        .match_indices('\t')
        .nth(4)
        .map_or(record, |(at, _)| &record[..at])
}

fn malformed_transcript(shard: &str) -> Error {
    Error::comparator(format!(
        "{shard}: a transcript has an id that is no string or exons that are no array of objects"
    ))
}

/// The digits of an exon bound: a whole number that is not below zero, or a string of
/// decimal digits.
fn bound_digits(bound: Option<&Value>) -> Option<String> {
    match bound? {
        Value::Int(whole) if *whole >= 0 => Some(whole.to_string()),
        Value::BigInt(whole) if !whole.starts_with('-') => Some(whole.clone()),
        Value::Str(text) if pyfmt::is_decimal(text) => Some(text.clone()),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache::testing::{plant, scratch};

    fn midway(reverse: bool) -> PlacedAllele {
        PlacedAllele {
            start: 1001,
            end: 1005,
            exons: vec![(990, 1000), (1006, 1020)],
            reverse,
        }
    }

    #[test]
    fn an_intron_base_midway_between_two_exons_counts_from_the_exon_before_it() {
        for reverse in [false, true] {
            let offsets: Vec<i64> = (0..5)
                .map(|offset| midway(reverse).intron_offset(offset))
                .collect();
            assert_eq!(offsets, [1, 2, 3, -2, -1]);
        }
        let abutting = PlacedAllele {
            exons: vec![(990, 1000), (1005, 1020)],
            ..midway(false)
        };
        assert_eq!(abutting.intron_offset(4), 0);
    }

    #[test]
    fn an_intron_base_is_numbered_from_the_nearer_exon_base_along_the_strand() {
        let places = |reverse: bool| -> Vec<(i64, i64)> {
            (0..5)
                .map(|offset| midway(reverse).exon_base(offset))
                .collect()
        };
        assert_eq!(
            places(false),
            [(11, 1), (11, 2), (11, 3), (12, -2), (12, -1)]
        );
        assert_eq!(
            places(true),
            [(15, 1), (15, 2), (15, 3), (16, -2), (16, -1)]
        );
        assert_eq!(midway(false).exon_base(-20), (0, 0));
        assert_eq!(midway(false).exon_base(30), (0, 0));
    }

    #[test]
    fn the_bases_outside_the_transcript_are_counted_on_its_strand() {
        let over_both = |reverse: bool| PlacedAllele {
            start: 5001,
            end: 5005,
            exons: vec![(5003, 5004)],
            reverse,
        };
        assert_eq!(over_both(false).cuts(), (2, 1));
        assert_eq!(over_both(true).cuts(), (1, 2));
        assert_eq!(midway(false).cuts(), (0, 0));
    }

    const TRANSCRIPT: &str = r#"{"stable_id": "ENST00000000010", "strand": 1, "exons": [{"start": 999990, "end": 1000003}]}"#;

    fn built_cache(name: &str) -> PathBuf {
        let root = scratch(name);
        plant(
            &root,
            "transcripts/21/1-1000000.json",
            &format!("[{TRANSCRIPT}]"),
        );
        plant(&root, "transcripts/21/1000001-2000000.json", "[]");
        root
    }

    fn record_entry() -> Entry {
        Entry {
            key: "0|ACGT|ENST00000000010|Transcript".to_string(),
            record: Some(0),
            allele: "ACGT".to_string(),
            feature: "ENST00000000010".to_string(),
            consequence: "3_prime_UTR_variant".to_string(),
            ..Entry::default()
        }
    }

    #[test]
    fn an_entry_without_a_location_is_placed_by_its_input_record() {
        let cache = built_cache("edges-record");
        let placed = |records: &[&str]| {
            let inputs: Vec<String> = records.iter().map(|record| record.to_string()).collect();
            TranscriptEdges::new(&cache, Some(&inputs))
                .expect("the cache has transcripts")
                .place(&record_entry())
                .expect("the shards read")
        };
        // Several alternates that share their first base, one of another length: the
        // parser drops that base, and the allele stands one position on.
        let chopped = placed(&["21\t1000001\tr\tTTTTT\tTACGT,TT\t.\t.\t."]).expect("placed");
        assert_eq!((chopped.start, chopped.end), (1000002, 1000005));
        assert_eq!(chopped.cuts(), (0, 2));
        assert_eq!(
            placed(&["21\t1000002\tr\tTTTT\tACGT\t.\t.\t."]),
            Some(chopped)
        );
        for record in [
            "21\t1000001\tr\tTTTTT\tAACGT\t.\t.\t.",
            "21\t1000002\tr\tTTTTT\tACGTT\t.\t.\t.",
            "21\t1000001\tr\tTTTTT\tTACGT,AT\t.\t.\t.",
        ] {
            assert_eq!(placed(&[record]), None, "{record}");
        }
        assert_eq!(placed(&[]), None);
        std::fs::remove_dir_all(cache).unwrap();
    }

    #[test]
    fn a_text_entry_is_placed_by_its_location_from_either_layout_of_the_cache() {
        let cache = built_cache("edges-location");
        let edges = TranscriptEdges::new(&cache, None).expect("the cache has transcripts");
        let entry = |location: &str, feature: &str| Entry {
            location: Some(location.to_string()),
            feature: feature.to_string(),
            ..record_entry()
        };
        let placed = |location: &str, feature: &str| {
            edges
                .place(&entry(location, feature))
                .expect("the shards read")
                .map(|placed| (placed.start, placed.end, placed.exons, placed.reverse))
        };
        let found = Some((1000002, 1000005, vec![(999990, 1000003)], false));
        assert_eq!(placed("21:1000002-1000005", "ENST00000000010"), found);
        assert_eq!(placed("chr21:1000002-1000005", "ENST00000000010.2"), found);
        assert_eq!(placed("21:1000002-1000004", "ENST00000000010"), None);
        assert_eq!(placed("21:1000004-1000007", "ENST00000000010"), None);
        assert_eq!(placed("22:1000002-1000005", "ENST00000000010"), None);
        assert_eq!(placed("21:1000002-1000005", "ENST00000000011"), None);
        let shards: Vec<i64> = edges
            .chromosome_shards("21")
            .expect("the shards list")
            .into_iter()
            .map(|(low, _)| low)
            .collect();
        assert_eq!(shards, [1000001, 1]);
        std::fs::remove_dir_all(cache).unwrap();
    }

    #[test]
    fn a_cache_without_transcripts_is_refused() {
        let root = scratch("edges-no-transcripts");
        assert_eq!(
            TranscriptEdges::new(&root, None).unwrap_err().to_string(),
            format!(
                "--json-cache: {} has no transcripts/ directory",
                root.display()
            )
        );
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_shard_files_each_transcript_under_its_id_with_its_spans_in_order() {
        let root = scratch("edges-shard");
        let shard = root.join("1-5.json.gz.json");
        std::fs::write(
            &shard,
            r#"{"1-5": [null, {"stable_id": "a", "strand": "-1", "exons": [{"start": "20", "end": "10"}, null, {"end": 7, "start": 5}]}],
                "6-9": [{"stable_id": "a", "strand": -1.0, "exons": [{"start": 30, "end": 40}]},
                        {"stable_id": "b", "exons": []}, {"stable_id": "", "exons": [{"start": 1, "end": 2}]},
                        {"stable_id": "c", "strand": 1, "exons": [{"start": 1, "end": 2}]}]}"#,
        )
        .unwrap();
        let filed = TranscriptEdges::shard_exons(&shard).expect("the shard reads");
        assert_eq!(
            filed.get("a"),
            Some(&vec![
                (vec![(5, 7), (10, 20)], true),
                (vec![(30, 40)], true)
            ])
        );
        assert_eq!(filed.get("c"), Some(&vec![(vec![(1, 2)], false)]));
        assert_eq!(filed.len(), 2);
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn a_malformed_shard_ends_with_the_comparators_error() {
        let root = scratch("edges-malformed");
        let refused = |name: &str, text: &str| {
            let shard = root.join(name);
            std::fs::write(&shard, text).unwrap();
            let error = TranscriptEdges::shard_exons(&shard)
                .unwrap_err()
                .to_string();
            error
                .strip_prefix(&shard.display().to_string())
                .unwrap_or_else(|| panic!("{error}"))
                .to_string()
        };
        let no_objects =
            ": a transcript has an id that is no string or exons that are no array of objects";
        let with =
            |member: &str, value: &str| format!("[{}]", TRANSCRIPT.replacen(member, value, 1));
        let exons = r#""exons": [{"start": 999990, "end": 1000003}]"#;
        assert_eq!(
            refused("exons.json", &with(exons, r#""exons": 7"#)),
            no_objects
        );
        assert_eq!(
            refused("exon.json", &with(exons, r#""exons": [7]"#)),
            no_objects
        );
        assert_eq!(
            refused(
                "id.json",
                &with(r#""ENST00000000010""#, r#"["ENST00000000010"]"#)
            ),
            no_objects
        );
        let no_integer =
            ": transcript 'ENST00000000010' has an exon without an integer start and end";
        for bound in [
            r#""first""#,
            "true",
            "1.5",
            "-1",
            "null",
            r#""""#,
            r#"" 5""#,
        ] {
            assert_eq!(
                refused("bound.json", &with("999990", bound)),
                no_integer,
                "{bound}"
            );
        }
        assert_eq!(
            refused("end.json", &with(r#", "end": 1000003"#, "")),
            no_integer
        );
        assert_eq!(
            refused("long.json", &with("999990", &format!("\"{}\"", "9".repeat(5000)))),
            ": an exon bound of transcript 'ENST00000000010', '999999999999999999999999...' is not an integer"
        );
        assert_eq!(
            refused("deep.json", &("[".repeat(100000) + &"]".repeat(100000))),
            " is nested too deep to read"
        );
        assert_eq!(
            refused("text.json", "{"),
            " is not JSON: Expecting property name enclosed in double quotes: line 1 column 2 (char 1)"
        );
        std::fs::remove_dir_all(root).unwrap();
    }
}
