// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One transcript of a JSON cache as the Downstream and NMD plugins read it: its coding
//! sequence and 3' UTR, its genetic code and the coordinates that map a genomic span to a
//! span of the coding sequence.

use super::codes::{genetic_code, translate_until_stop, translated_length, GeneticCode};
use crate::pytext::{self, Json};
use crate::rules::shared::cache_integer;

/// The keys of one pair of a cache transcript's `mapper`: a run of cDNA coordinates, the
/// genomic bases it maps to and their orientation.
const MAPPER_PAIR_KEYS: [&str; 5] = ["from_start", "from_end", "to_start", "to_end", "ori"];

/// The genetic code a cache transcript's `codon_table` names: a number in a native cache,
/// the dump's string in a converted one, table 1 for 0 or no key (Downstream.pm:132);
/// `None` for a table the codes lack and for any other value.
pub fn codon_table_id(value: Option<&Json>) -> Option<u8> {
    let text = match value {
        None | Some(Json::Null) => return Some(1),
        Some(Json::Int(digits)) => digits.as_str(),
        Some(Json::Str(text)) => text.as_str(),
        Some(_) => return None,
    };
    if !pytext::isdecimal(text) {
        return None;
    }
    let table = match pytext::int(text)? {
        0 => 1,
        table => u8::try_from(table).ok()?,
    };
    genetic_code(table).map(|_| table)
}

/// The byte at which the character `index` of `text` starts, the text's length for an index
/// at or past its end.
fn byte_offset(text: &str, ascii: bool, index: usize) -> usize {
    if ascii {
        return index.min(text.len());
    }
    text.char_indices()
        .nth(index)
        .map_or(text.len(), |(offset, _)| offset)
}

/// A slice bound as Python reads one: an index below zero counts from the end, and every
/// index is held inside the text.
fn bound(index: Option<i128>, length: usize, open: usize) -> usize {
    match index {
        None => open,
        Some(index) if index < 0 => (length as i128 + index).max(0) as usize,
        Some(index) => index.min(length as i128) as usize,
    }
}

/// A text indexed by character, as Python indexes a `str`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    text: Box<str>,
    /// The number of characters.
    length: usize,
}

impl Sequence {
    pub fn new(text: String) -> Sequence {
        let length = if text.is_ascii() {
            text.len()
        } else {
            text.chars().count()
        };
        Sequence {
            text: text.into_boxed_str(),
            length,
        }
    }

    /// The number of characters.
    pub fn len(&self) -> usize {
        self.length
    }

    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    pub fn as_str(&self) -> &str {
        &self.text
    }

    fn ascii(&self) -> bool {
        self.length == self.text.len()
    }

    /// `text[from:to]`: the characters from the index `from` up to the index `to`, an index
    /// below zero counted from the end and an open bound the text's own.
    pub fn slice(&self, from: Option<i128>, to: Option<i128>) -> &str {
        let (from, to) = (
            bound(from, self.length, 0),
            bound(to, self.length, self.length),
        );
        if from >= to {
            return "";
        }
        let ascii = self.ascii();
        &self.text[byte_offset(&self.text, ascii, from)..byte_offset(&self.text, ascii, to)]
    }
}

/// A cache transcript's coordinates: its `strand`, its `exons` in transcript order as
/// (start, end), its `cdna_coding_start` and `cdna_coding_end`, the count of `N` its
/// `translateable_seq` opens with (the start exon's phase, which the CDS coordinates count
/// in) and its `mapper` pairs as (cDNA start, cDNA end, genomic start, genomic end,
/// orientation).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranscriptModel {
    pub strand: i8,
    pub exons: Vec<(i128, i128)>,
    pub coding_start: i128,
    pub coding_end: i128,
    pub padding: i128,
    pub pairs: Vec<[i128; 5]>,
}

/// `x + y - z` where it lies in `i128`, whatever the size of a sum of two of the three: where
/// the result lies in `i128` one of the three orders passes through no sum outside it.
fn sum(x: i128, y: i128, z: i128) -> Option<i128> {
    x.checked_add(y)
        .and_then(|both| both.checked_sub(z))
        .or_else(|| x.checked_sub(z).and_then(|rest| rest.checked_add(y)))
        .or_else(|| y.checked_sub(z).and_then(|rest| rest.checked_add(x)))
}

impl TranscriptModel {
    /// The model of a cache transcript object, `None` where the cache states a part of it
    /// in no form this reads.
    pub fn from_cache(transcript: &Json, padding: i128) -> Option<TranscriptModel> {
        let mapper = transcript
            .get("variation_effect_feature_cache")
            .and_then(|sequences| sequences.get("mapper"));
        let pairs = match mapper {
            Some(mapper @ Json::Object(_)) => mapper.get("pairs"),
            _ => None,
        };
        let strand = match cache_integer(transcript.get("strand"))? {
            1 => 1,
            -1 => -1,
            _ => return None,
        };
        let coding_start = cache_integer(transcript.get("cdna_coding_start"))?;
        let coding_end = cache_integer(transcript.get("cdna_coding_end"))?;
        let (Some(Json::List(pairs)), Some(Json::List(exons))) = (pairs, transcript.get("exons"))
        else {
            return None;
        };
        if pairs.is_empty() || exons.is_empty() {
            return None;
        }
        let object = |item: &Json| matches!(item, Json::Object(_));
        if !pairs.iter().chain(exons).all(object) {
            return None;
        }
        // Every pair and every exon is read before one is refused, as the cache's own
        // objects are: a pair that lacks a number refuses the model.
        let read_pairs: Vec<[Option<i128>; 5]> = pairs
            .iter()
            .map(|pair| MAPPER_PAIR_KEYS.map(|key| cache_integer(pair.get(key))))
            .collect();
        let read_exons: Vec<(Option<i128>, Option<i128>)> = exons
            .iter()
            .map(|exon| {
                (
                    cache_integer(exon.get("start")),
                    cache_integer(exon.get("end")),
                )
            })
            .collect();
        let pairs = read_pairs
            .into_iter()
            .map(|pair| Some([pair[0]?, pair[1]?, pair[2]?, pair[3]?, pair[4]?]))
            .collect::<Option<Vec<_>>>()?;
        let mut exons = read_exons
            .into_iter()
            .map(|(start, end)| Some((start?, end?)))
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

    /// The cDNA coordinate of a genomic base: `None` for a base outside the mapper pairs,
    /// `Some(None)` for a coordinate outside `i128`.
    fn cdna_of(&self, genomic: i128) -> Option<Option<i128>> {
        self.pairs
            .iter()
            .find(|pair| pair[2] <= genomic && genomic <= pair[3])
            .map(
                |&[cdna_start, _, genomic_start, genomic_end, orientation]| {
                    if orientation == -1 {
                        sum(cdna_start, genomic_end, genomic)
                    } else {
                        sum(cdna_start, genomic, genomic_start)
                    }
                },
            )
    }

    /// The cDNA coordinate of a genomic base, `None` for a base outside the mapper pairs and
    /// for a coordinate outside `i128`.
    pub fn cdna(&self, genomic: i128) -> Option<i128> {
        self.cdna_of(genomic).flatten()
    }

    /// The CDS coordinate `past` bases after the one of a coding cDNA coordinate, `None`
    /// outside `i128`.
    fn cds(&self, cdna: i128, past: i128) -> Option<i128> {
        cdna.checked_sub(self.coding_start)?
            .checked_add(self.padding)?
            .checked_add(past)
    }

    /// The genomic base of a CDS coordinate, `None` for a coordinate outside the coding
    /// bounds or the mapper pairs, and for a base outside `i128`.
    pub fn genomic(&self, cds: i128) -> Option<i128> {
        // The coordinate's place in the coding region, from 0: one that leaves `i128` here
        // or as a cDNA coordinate lies outside the coding bounds.
        let place = cds.checked_sub(self.padding)?.checked_sub(1)?;
        let cdna = self.coding_start.checked_add(place)?;
        if !(self.coding_start <= cdna && cdna <= self.coding_end) {
            return None;
        }
        self.pairs
            .iter()
            .find(|pair| pair[0] <= cdna && cdna <= pair[1])
            .and_then(
                |&[cdna_start, _, genomic_start, genomic_end, orientation]| {
                    if orientation == -1 {
                        sum(genomic_end, cdna_start, cdna)
                    } else {
                        sum(genomic_start, cdna, cdna_start)
                    }
                },
            )
    }

    /// The `cds_start` and `cds_end` of the genomic span `start` to `end`
    /// (BaseTranscriptVariation.pm:252-264 over `TranscriptMapper::genomic2cds`), `None`
    /// for an end that is no coding base: the CDS coordinates of the span's 5' and 3'
    /// bases; for an insertion (`start` is `end` + 1) `map_insert`'s result
    /// (Mapper.pm:504-558), the base after the 5' flank's coordinate and the base before
    /// the 3' flank's, a flank outside the mapper pairs read as the base beside the other.
    pub fn cds_pair(&self, start: i128, end: i128) -> (Option<i128>, Option<i128>) {
        if end.checked_add(1) == Some(start) {
            return self.insertion_pair(Some(start), Some(end));
        }
        let (five, three) = if self.strand == 1 {
            (start, end)
        } else {
            (end, start)
        };
        let coding = |cdna: i128| {
            (self.coding_start <= cdna && cdna <= self.coding_end)
                .then(|| self.cds(cdna, 1))
                .flatten()
        };
        (
            self.cdna(five).and_then(coding),
            self.cdna(three).and_then(coding),
        )
    }

    /// `cds_pair` of an insertion, whose start is one past its end: a flank outside `i128`,
    /// which no mapper pair holds, is `None`.
    pub fn insertion_pair(
        &self,
        start: Option<i128>,
        end: Option<i128>,
    ) -> (Option<i128>, Option<i128>) {
        let (five, three) = if self.strand == 1 {
            (end, start)
        } else {
            (start, end)
        };
        let cdna_of = |flank: Option<i128>| flank.and_then(|flank| self.cdna_of(flank));
        // A cDNA coordinate outside `i128` is no coding base, and neither is the one
        // beside it.
        let (before, after) = match (cdna_of(five), cdna_of(three)) {
            (None, None) => return (None, None),
            (None, Some(after)) => (after.and_then(|after| after.checked_sub(1)), after),
            (Some(before), None) => (before, before.and_then(|before| before.checked_add(1))),
            (Some(before), Some(after)) => (before, after),
        };
        (
            before
                .filter(|before| self.coding_start <= *before && *before < self.coding_end)
                .and_then(|before| self.cds(before, 2)),
            after
                .filter(|after| self.coding_start < *after && *after <= self.coding_end)
                .and_then(|after| self.cds(after, 0)),
        )
    }

    /// The one genomic span (start, end; an insertion's start one past its end) whose
    /// `cds_pair` is `cds_start` and `cds_end`, or `None`: an end the mapper pairs do not
    /// hold, or an insertion's pair (`cds_start` is `cds_end` + 1) between two bases the
    /// genome does not keep adjacent, which the point after the one and the point before
    /// the other both give.
    pub fn genomic_span(&self, cds_start: i128, cds_end: i128) -> Option<(i128, i128)> {
        let first = self.genomic(cds_start.min(cds_end))?;
        let last = self.genomic(cds_start.max(cds_end))?;
        if cds_start > cds_end {
            (first.abs_diff(last) == 1).then(|| (first.max(last), first.min(last)))
        } else {
            Some((first.min(last), first.max(last)))
        }
    }

    /// `variant_exon_check` (NMD.pm:132-176) on a record's genomic end: it lies in the last
    /// exon, or in the second to last within 51 bases of that exon's 3' end, its end on the
    /// forward strand and its start on the reverse.
    pub fn ends_where_nmd_escapes(&self, genomic_end: i128) -> bool {
        let Some(&(start, end)) = self.exons.last() else {
            return false;
        };
        if start <= genomic_end && genomic_end <= end {
            return true;
        }
        if self.exons.len() < 2 {
            return false;
        }
        let (start, end) = self.exons[self.exons.len() - 2];
        // A distance outside `i128` is more than 51 bases.
        let within = |from: i128, to: i128| {
            to.checked_sub(from)
                .is_some_and(|bases| (0..=51).contains(&bases))
        };
        if self.strand == -1 {
            within(start, genomic_end)
        } else {
            within(genomic_end, end)
        }
    }
}

/// What the Downstream and NMD plugins read of one cache transcript: the
/// `translateable_seq`, the `three_prime_utr` (empty where the cache spells none), the
/// codon table (`codon_table_id`, `None` for a table the codes lack) and the length of the
/// `peptide` (`None` without one) of its `variation_effect_feature_cache`, and its
/// `TranscriptModel` (`None` where the cache does not state one).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DownstreamTranscript {
    pub coding: Sequence,
    pub utr: Sequence,
    pub table: Option<u8>,
    pub peptide_length: Option<usize>,
    pub model: Option<TranscriptModel>,
}

/// The characters of a text from its character `skip` on.
fn after(text: &str, skip: usize) -> &str {
    &text[byte_offset(text, text.is_ascii(), skip)..]
}

impl DownstreamTranscript {
    /// The stable id and the reading of a cache transcript object, `None` for one without a
    /// stable id, without a `translateable_seq` or with a 3' UTR that is no text.
    pub fn from_cache(transcript: &Json) -> Option<(&str, DownstreamTranscript)> {
        let Some(Json::Str(stable_id)) = transcript.get("stable_id") else {
            return None;
        };
        let sequences = match transcript.get("variation_effect_feature_cache") {
            Some(sequences @ Json::Object(_)) => sequences,
            _ => return None,
        };
        let coding = match sequences.get("translateable_seq") {
            Some(Json::Str(coding)) if !coding.is_empty() => coding.to_uppercase(),
            _ => return None,
        };
        let utr = match sequences.get("three_prime_utr") {
            Some(Json::Str(utr)) => utr.to_uppercase(),
            None | Some(Json::Null) => String::new(),
            Some(_) => return None,
        };
        let peptide_length = match sequences.get("peptide") {
            Some(Json::Str(peptide)) => Some(peptide.chars().count()),
            _ => None,
        };
        let padding = coding.chars().take_while(|letter| *letter == 'N').count();
        let read = DownstreamTranscript {
            coding: Sequence::new(coding),
            utr: Sequence::new(utr),
            table: codon_table_id(sequences.get("codon_table")),
            peptide_length,
            model: TranscriptModel::from_cache(transcript, padding as i128),
        };
        Some((stable_id, read))
    }

    /// Whether the two readings are one but for their models.
    pub fn same_sequences(&self, other: &DownstreamTranscript) -> bool {
        self.coding == other.coding
            && self.utr == other.utr
            && self.table == other.table
            && self.peptide_length == other.peptide_length
    }

    /// The transcript's genetic code, `None` under a table the codes lack.
    pub fn code(&self) -> Option<&'static GeneticCode> {
        self.table.and_then(genetic_code)
    }

    /// The Downstream plugin's value on the CDS pair `start` and `end` (an insertion's
    /// start one past its end) with the bases `allele` on the transcript's strand: the
    /// DownstreamProtein, or with `length_change` the ProteinLengthChange as printed. Both
    /// rules write the allele over the pair (Downstream.pm:108-110) and translate the
    /// changed coding sequence with the 3' UTR up to the first stop; the
    /// ProteinLengthChange is the whole translation's length less the peptide's (:141-148;
    /// the translated coding sequence's for a transcript the cache gives no peptide). With
    /// `ensembl` the DownstreamProtein starts at the pair's lower coordinate rounded down
    /// to a multiple of three (:112-115) and there is no value for a start past the coding
    /// sequence's end, where the assignment of :110 fails; without it, vep-rs's rule, at
    /// the codon holding the first changed base. `None` as well under a table the codes
    /// lack.
    pub fn downstream_value(
        &self,
        start: i128,
        end: i128,
        allele: &str,
        length_change: bool,
        ensembl: bool,
    ) -> Option<String> {
        let code = self.code()?;
        let length = self.coding.len() as i128;
        // One below the smallest `i128` is read as that one: both lie before the sequence.
        let before_start = start.saturating_sub(1);
        if ensembl && before_start > length {
            return None;
        }
        let head = self.coding.slice(None, Some(before_start));
        let tail = self.coding.slice(Some(end.max(before_start)), None);
        let utr = self.utr.as_str();
        if length_change {
            let reference = match self.peptide_length {
                Some(length) => length,
                None => translated_length(self.coding.as_str().chars(), code),
            };
            let changed = [head, allele, tail, utr];
            let translated =
                translated_length(changed.iter().flat_map(|piece| piece.chars()), code);
            return Some((translated as i128 - reference as i128).to_string());
        }
        let first = if ensembl {
            start.min(end)
        } else {
            before_start.min(length)
        };
        // The changed sequence from its character `first - first % 3` on, as a slice reads
        // an index below zero and one past the end.
        let pieces = [head, allele, tail];
        let lengths = pieces.map(|piece| {
            if piece.is_ascii() {
                piece.len()
            } else {
                piece.chars().count()
            }
        });
        let codon_start = first.saturating_sub(first.rem_euclid(3));
        let mut skip = bound(Some(codon_start), lengths.iter().sum(), 0);
        let mut kept = [""; 3];
        for (at, piece) in pieces.iter().enumerate() {
            let dropped = skip.min(lengths[at]);
            kept[at] = after(piece, dropped);
            skip -= dropped;
        }
        let read = kept.iter().chain([&utr]).flat_map(|piece| piece.chars());
        Some(translate_until_stop(read, code))
    }
}
