// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The Downstream field classes: the sequences and the coding models of the corpus cache's
//! transcripts, translation under a transcript's codon table as `Bio::Tools::CodonTable`
//! translates, the span Ensembl VEP's VCF parser gives a record, and the three classes of a
//! `DownstreamProtein` or `ProteinLengthChange` the plugin reads at another span or in
//! another frame than vep-rs.

use super::cache::iter_cache_transcripts;
use super::fields::{FieldModels, FieldQuestion, FieldRow, FieldRun, Verdict};
use crate::error::Failure;
use crate::json::{Object, Value};
use crate::names::{first, length, record_pos};
use crate::pyfmt::is_py_space;
use crate::pytext::{py_int, py_repr_str, py_str};
use std::collections::{BTreeSet, HashMap, HashSet};

/// The fields the Downstream plugin writes, in the order `downstream_plugin_values` and
/// `downstream_own_values` give their values.
pub const DOWNSTREAM_FIELDS: [&str; 2] = ["DownstreamProtein", "ProteinLengthChange"];

/// The class of a `DownstreamProtein` on a row of a deletion, with or without inserted
/// bases, whose first base is a codon's third: Ensembl VEP's Downstream plugin rounds that
/// base's coordinate down to a codon boundary, the end of the codon the deletion changes,
/// and translates from the codon after it (Downstream.pm:112-115), where vep-rs translates
/// from the changed codon.
pub const DOWNSTREAM_START_PAST_OWN_CODON: &str = "downstream_start_past_own_codon";

/// The class of such a row whose deletion inserts no base and runs to the last coding base:
/// the changed coding sequence ends one base before that boundary, so the plugin takes
/// nothing of it (Downstream.pm:115) and translates the 3' UTR from its first base
/// (Downstream.pm:116-119), one base out of the changed sequence's frame, where vep-rs
/// translates the two bases the changed codon keeps with the UTR, in frame.
pub const DOWNSTREAM_UTR_FRAME_PAST_CODING_END: &str = "downstream_utr_frame_past_coding_end";

/// The class of a `DownstreamProtein` or `ProteinLengthChange` on a row whose HGVS notation
/// is shifted 3': under `--hgvs` Ensembl VEP's `hgvs_protein` reads the row's CDS pair again
/// at the shifted span and leaves it there (TranscriptVariationAllele.pm:844-858,
/// :2406-2410) with the allele's bases rotated for that span (:1330-1345, :1778), and the
/// Downstream plugin, run after it (OutputFactory.pm:398), computes on both
/// (Downstream.pm:108-110), so its two values change with a flag that asks for a notation;
/// vep-rs computes on the record's own span and the allele as written, with the flag and
/// without it.
pub const PLUGIN_HGVS_SHIFTED_CDS_COORDINATES: &str = "plugin_hgvs_shifted_cds_coordinates";

/// The keys of a cached mapper pair, in the order `CodingModel::pairs` holds their values:
/// the pair's first and last cDNA base, its first and last genomic base and its orientation.
pub const MAPPER_PAIR_KEYS: [&str; 5] = ["from_start", "from_end", "to_start", "to_end", "ori"];

/// The bases of a codon in the order the genetic codes list their codons: the first base
/// moves slowest and the third fastest.
const CODON_BASES: [char; 4] = ['T', 'C', 'A', 'G'];

/// The genetic codes `Bio::Tools::CodonTable` translates with, by table id, as `CODES` of
/// crates/vep-core/src/codon.rs lists them: the amino acid of each of the 64 codons in the
/// order of `CODON_BASES`. Ids 7, 8 and 17 to 20 hold no table there and ids from 25 on
/// are unknown to it, so a transcript under one of them has no code here.
pub const CODON_TABLE_AMINO_ACIDS: [(i64, &str); 18] = [
    (
        1,
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        2,
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    ),
    (
        3,
        "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        4,
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        5,
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        6,
        "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        9,
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        10,
        "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        11,
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        12,
        "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        13,
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    ),
    (
        14,
        "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        15,
        "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        16,
        "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        21,
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    ),
    (
        22,
        "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        23,
        "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
    (
        24,
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    ),
];

/// The longest run of digits `int()` reads.
const MAX_INT_TEXT_DIGITS: usize = 4300;

/// The amino acids of the genetic code of id `table`, `None` for an id that has no code.
pub fn genetic_code(table: i64) -> Option<&'static str> {
    CODON_TABLE_AMINO_ACIDS
        .iter()
        .find(|(id, _)| *id == table)
        .map(|(_, amino_acids)| *amino_acids)
}

/// The bases an IUPAC nucleotide letter stands for, `U` read as `T`; none for any other
/// character.
fn iupac_bases(letter: char) -> &'static str {
    match letter {
        'T' | 'U' => "T",
        'C' => "C",
        'A' => "A",
        'G' => "G",
        'Y' => "TC",
        'W' => "TA",
        'K' => "TG",
        'M' => "CA",
        'S' => "CG",
        'R' => "AG",
        'H' => "TCA",
        'B' => "TCG",
        'D' => "TAG",
        'V' => "CAG",
        'N' | 'X' => "TCAG",
        _ => "",
    }
}

/// The amino acid `code` gives a codon of unambiguous bases, `None` for any other codon.
fn coded_amino_acid(codon: [char; 3], code: &str) -> Option<char> {
    let mut place = 0;
    for base in codon {
        place = place * 4 + CODON_BASES.iter().position(|held| *held == base)?;
    }
    code.as_bytes().get(place).copied().map(char::from)
}

/// `Bio::Tools::CodonTable::_translate_ambiguous_codon` under the genetic code `code`: the
/// amino acid of an upper-case codon that is no unambiguous codon, from the distinct amino
/// acids its readings code: one, that one; aspartate and asparagine `B`; glutamate and
/// glutamine `Z`; anything else `X`, as is a codon with a character that is no nucleotide
/// letter.
pub fn ambiguous_codon_amino_acid(codon: [char; 3], code: &str) -> char {
    let [firsts, seconds, thirds] = codon.map(iupac_bases);
    let mut amino_acids = BTreeSet::new();
    for first_base in firsts.chars() {
        for second_base in seconds.chars() {
            for third_base in thirds.chars() {
                amino_acids.extend(coded_amino_acid(
                    [first_base, second_base, third_base],
                    code,
                ));
            }
        }
    }
    let amino_acids: Vec<char> = amino_acids.into_iter().collect();
    match amino_acids.as_slice() {
        [one] => *one,
        ['D', 'N'] => 'B',
        ['E', 'Q'] => 'Z',
        _ => 'X',
    }
}

/// The translation of `sequence` under `code`: its whole codons in order, in either case.
fn translate(sequence: impl Iterator<Item = char>, code: &str) -> String {
    let upper: Vec<char> = sequence.flat_map(char::to_uppercase).collect();
    let amino_acid = |codon: &[char]| {
        let codon = [codon[0], codon[1], codon[2]];
        coded_amino_acid(codon, code).unwrap_or_else(|| {
            if codon == ['-'; 3] {
                '-'
            } else {
                ambiguous_codon_amino_acid(codon, code)
            }
        })
    };
    upper.chunks_exact(3).map(amino_acid).collect()
}

/// `Bio::Tools::CodonTable::translate` under the genetic code of id `table`, as
/// `Bio::Seq::translate` calls it for a sequence with no frame and no partial-codon reading
/// and as `translate_sequence` of crates/vep-core/src/codon.rs restates it: the sequence's
/// whole codons in order, in either case, a trailing one or two bases dropped; a codon of
/// unambiguous bases by the table, three hyphens as `-`, any other by
/// `ambiguous_codon_amino_acid`. A stop stays `*`. `None` for an id that has no code.
pub fn translate_codons(sequence: &str, table: i64) -> Option<String> {
    Some(translate(sequence.chars(), genetic_code(table)?))
}

/// A translation up to its first stop.
fn before_stop(peptide: &str) -> &str {
    peptide.split('*').next().unwrap_or_default()
}

/// A count as a number that takes a sign.
fn signed(count: usize) -> i128 {
    i128::try_from(count).unwrap_or(i128::MAX)
}

/// The `ProteinLengthChange` of a translation: its length up to its first stop less the
/// cached peptide's.
fn length_change(whole: &str, peptide_length: usize) -> String {
    (signed(length(before_stop(whole))) - signed(peptide_length)).to_string()
}

/// `coding` with `bases` in place of the bases after its first `kept` up to base `end`,
/// counted from 1.
fn changed_sequence(coding: &[char], kept: usize, end: usize, bases: &str) -> Vec<char> {
    let after = end.max(kept).min(coding.len());
    let mut changed = coding[..kept].to_vec();
    changed.extend(bases.chars());
    changed.extend_from_slice(&coding[after..]);
    changed
}

/// Ensembl VEP's `DownstreamProtein` and `ProteinLengthChange` (Downstream.pm:104-149), in
/// the order of `DOWNSTREAM_FIELDS`, under the codon table `table` for a change that puts
/// `bases`, the allele on the transcript's strand, in place of the coding bases `start` to
/// `end` (counted from 1 on `coding`, the transcript's `translateable_seq`; an insertion's
/// start lies one past its end), and whether the protein's translation starts past the
/// changed coding sequence; `None` where `start` lies more than one base past `coding`,
/// which ends the plugin's call (:110), and for a table that has no code. The protein is
/// the translation of the changed sequence from a codon boundary, the 3' UTR appended, up
/// to its first stop (:115-140). The boundary is the lower of `start` and `end` rounded
/// down to a multiple of three (:112-113), the end of the codon holding that base when it
/// is the codon's third, and a boundary past the changed sequence's end takes none of it
/// (:115), which leaves the UTR translated from its first base. The length change is the
/// length of the whole changed sequence's translation with the UTR, up to its first stop,
/// less `peptide_length`, the cached peptide's (:124-148).
pub fn downstream_plugin_values(
    coding: &str,
    utr: &str,
    table: i64,
    peptide_length: usize,
    start: usize,
    end: usize,
    bases: &str,
) -> Option<([String; 2], bool)> {
    let code = genetic_code(table)?;
    let coding: Vec<char> = coding.chars().collect();
    let kept = start.saturating_sub(1);
    if kept > coding.len() {
        return None;
    }
    let changed = changed_sequence(&coding, kept, end, bases);
    let low = start.min(end);
    let boundary = low - low % 3;
    let from_boundary = changed.get(boundary..).unwrap_or_default();
    let protein = translate(from_boundary.iter().copied().chain(utr.chars()), code);
    let whole = translate(changed.iter().copied().chain(utr.chars()), code);
    let values = [
        before_stop(&protein).to_string(),
        length_change(&whole, peptide_length),
    ];
    Some((values, boundary > changed.len()))
}

/// vep-rs's `DownstreamProtein` and `ProteinLengthChange` for the change
/// `downstream_plugin_values` describes (`translate_change` of
/// crates/vep-builtins/src/plugins/downstream.rs): the span is clipped to `coding`, the
/// protein is read from the first base of the codon holding the first changed base, and the
/// length change is computed as Ensembl VEP's. `None` for a table that has no code.
pub fn downstream_own_values(
    coding: &str,
    utr: &str,
    table: i64,
    peptide_length: usize,
    start: usize,
    end: usize,
    bases: &str,
) -> Option<[String; 2]> {
    let code = genetic_code(table)?;
    let coding: Vec<char> = coding.chars().collect();
    let kept = start.saturating_sub(1).min(coding.len());
    let changed = changed_sequence(&coding, kept, end, bases);
    let own_codon = kept - kept % 3;
    let from_own_codon = changed[own_codon..].iter().copied();
    let protein = translate(from_own_codon.chain(utr.chars()), code);
    let whole = translate(changed.iter().copied().chain(utr.chars()), code);
    Some([
        before_stop(&protein).to_string(),
        length_change(&whole, peptide_length),
    ])
}

/// The bases `shift_feature_seqs` leaves on a row whose notation moved `offset` bases
/// (TranscriptVariationAllele.pm:1330-1345): the allele as written, its first base moved to
/// its end once per base of the shift on a forward-strand transcript and its length less
/// that many times on a reverse-strand one (`reverse`), not at all where that count is not
/// positive.
pub fn hgvs_shifted_allele(allele: &str, offset: i128, reverse: bool) -> String {
    let mut bases: Vec<char> = allele.chars().collect();
    let held = u128::try_from(bases.len()).unwrap_or(u128::MAX);
    let shift = offset.unsigned_abs();
    let turns = if reverse {
        held.saturating_sub(shift)
    } else {
        shift
    };
    if held > 0 && turns > 0 {
        bases.rotate_left(usize::try_from(turns % held).unwrap_or(0));
    }
    bases.into_iter().collect()
}

/// The genomic start and end Ensembl VEP's VCF parser gives a record (Parser/VCF.pm:254,
/// :295-336): POS and the last REF base, the start one base later where an ALT differs from
/// REF in length and REF and the ALTs open with one base, an ALT naming `*` left out of
/// several, so an insertion's start lies one past its end. A line of fewer than five
/// columns, a POS that is no integer and a span 64 bits do not hold are refused.
pub fn record_span(vcf_line: &str) -> Result<(i64, i64), Failure> {
    let columns: Vec<&str> = vcf_line.split('\t').collect();
    let [chrom, pos, id, reference, alternates, ..] = columns.as_slice() else {
        return Err(Failure::own(format!(
            "the record line {} has fewer than five columns",
            py_repr_str(vcf_line)
        )));
    };
    let pos = record_pos(chrom, pos, id)?;
    let alternates: Vec<&str> = alternates.split(',').collect();
    let reference_length = length(reference);
    let indel = alternates
        .iter()
        .any(|alternate| length(alternate) != reference_length);
    let first_bases: HashSet<&str> = std::iter::once(reference)
        .chain(alternates.iter())
        .filter(|allele| alternates.len() == 1 || !allele.contains('*'))
        .map(|allele| first(allele))
        .collect();
    let start = if indel && first_bases.len() == 1 {
        pos.checked_add(1)
    } else {
        Some(pos)
    };
    let end = i64::try_from(reference_length)
        .ok()
        .and_then(|bases| pos.checked_add(bases - 1));
    match (start, end) {
        (Some(start), Some(end)) => Ok((start, end)),
        _ => Err(Failure::own(format!(
            "the record {chrom} {pos} {id} has a span 64 bits do not hold"
        ))),
    }
}

/// The one span among the records `vcf_lines` (`record_span`) that a row's Location prints,
/// an insertion's two flanks lower first, or `None` where none or more than one does.
pub fn row_record_span<'a>(
    location: &str,
    vcf_lines: impl IntoIterator<Item = &'a str>,
) -> Result<Option<(i64, i64)>, Failure> {
    let printed = location.rsplit_once(':').map_or(location, |(_, span)| span);
    let mut spans: HashSet<(i64, i64)> = HashSet::new();
    for vcf_line in vcf_lines {
        let span = record_span(vcf_line)?;
        let (low, high) = (span.0.min(span.1), span.0.max(span.1));
        let text = if low == high {
            low.to_string()
        } else {
            format!("{low}-{high}")
        };
        if printed == text {
            spans.insert(span);
        }
    }
    let mut spans = spans.into_iter();
    match (spans.next(), spans.next()) {
        (Some(span), None) => Ok(Some(span)),
        _ => Ok(None),
    }
}

/// The sequences Ensembl VEP's Downstream plugin translates for a transcript and the table
/// it translates them under.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CodingSequence {
    /// The transcript's `translateable_seq` (Downstream.pm:104-106).
    pub coding: String,
    /// Its `three_prime_utr`, the empty sequence when the cache spells none (the plugin
    /// reads the transcript's own, :116).
    pub utr: String,
    /// Its `codon_table`, 1 when the cache names none or 0 (:130-133).
    pub codon_table: i64,
}

/// The coding transcripts of a corpus cache by stable id.
pub type CodingSequences = HashMap<String, CodingSequence>;

/// What places a genomic base on a cached transcript's coding sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CodingModel {
    /// The mapper pairs, each the values of `MAPPER_PAIR_KEYS`.
    pub pairs: Vec<[i64; 5]>,
    /// The first cDNA base of the coding sequence.
    pub coding_start: i64,
    /// Its last cDNA base.
    pub coding_end: i64,
    /// The length of the cached peptide, `None` where the cache holds none.
    pub peptide_length: Option<usize>,
}

/// The mapped coding transcripts of a corpus cache by chromosome and stable id.
pub type CodingModels = HashMap<(String, String), CodingModel>;

impl CodingModel {
    /// The cDNA coordinate of a genomic base, or `None` for a base no pair holds.
    /// Coordinates are computed in 128 bits, which hold every sum of a cache's numbers.
    pub fn cdna(&self, genomic: i128) -> Option<i128> {
        self.pairs.iter().find_map(|pair| {
            let [from_start, _, to_start, to_end, orientation] = pair.map(i128::from);
            (to_start <= genomic && genomic <= to_end).then(|| {
                let along = if orientation == -1 {
                    to_end - genomic
                } else {
                    genomic - to_start
                };
                from_start + along
            })
        })
    }

    /// The CDS start and end Ensembl VEP gives the genomic span `start` to `end` on the
    /// transcript, `reverse` for one on the reverse strand, counted from 1 on a
    /// `translateable_seq` that opens with `padding` bases of start-exon phase
    /// (BaseTranscriptVariation.pm:252-264), or `None` where either is no coding base. A
    /// span of bases maps by its two ends. An insertion, whose start lies one past its end,
    /// takes the bases between its two flanks' coordinates (Mapper.pm:504-558): none
    /// between neighbours, which leaves the start one past the end, and the bases the model
    /// carries where two pairs touch on the genome; where no pair holds one flank it lies
    /// against the other. Each flank's coordinate is read alone, so a flank outside the
    /// coding sequence leaves no pair even where the bases between reach into it
    /// (TranscriptMapper.pm:441-443).
    pub fn cds_pair(
        &self,
        start: i128,
        end: i128,
        reverse: bool,
        padding: usize,
    ) -> Option<(i128, i128)> {
        let coding_start = i128::from(self.coding_start);
        let coding_end = i128::from(self.coding_end);
        let (first, last) = if Some(start) == end.checked_add(1) {
            let (before, after) = if reverse {
                (self.cdna(start), self.cdna(end))
            } else {
                (self.cdna(end), self.cdna(start))
            };
            if before.is_some_and(|before| before < coding_start)
                || after.is_some_and(|after| after > coding_end)
            {
                return None;
            }
            match (before, after) {
                (None, None) => return None,
                (Some(before), None) => (before + 1, before),
                (None, Some(after)) => (after, after - 1),
                (Some(before), Some(after)) => (before + 1, after - 1),
            }
        } else if reverse {
            (self.cdna(end)?, self.cdna(start)?)
        } else {
            (self.cdna(start)?, self.cdna(end)?)
        };
        let coding = coding_start..=coding_end;
        if !coding.contains(&first) || !coding.contains(&last) {
            return None;
        }
        let to_cds = signed(padding) + 1 - coding_start;
        Some((first + to_cds, last + to_cds))
    }
}

/// A transcript's stable id as a message names it.
fn stable_id_text(transcript: &Object) -> String {
    transcript.get("stable_id").map(py_str).unwrap_or_default()
}

/// The `variation_effect_feature_cache` of a transcript record, `None` for a record that
/// has none or an empty one. Anything else that is no object is refused.
pub fn effect_feature_cache(transcript: &Object) -> Result<Option<&Object>, Failure> {
    match transcript.get("variation_effect_feature_cache") {
        Some(Value::Object(cache)) => Ok(Some(cache)),
        Some(other) if other.truthy() => Err(Failure::own(format!(
            "the cache record of {} holds a variation_effect_feature_cache that is no object",
            stable_id_text(transcript)
        ))),
        _ => Ok(None),
    }
}

/// The stable id a transcript record is filed under, `None` for an id that is no text,
/// which no row names. An id that is a list or an object is refused.
pub fn filed_stable_id(transcript: &Object) -> Result<Option<&str>, Failure> {
    match transcript.get("stable_id") {
        Some(Value::Str(stable_id)) => Ok(Some(stable_id)),
        Some(Value::List(_) | Value::Object(_)) => Err(Failure::own(
            "a cache record holds a stable_id that is a list or an object",
        )),
        _ => Ok(None),
    }
}

/// The number a `codon_table` text spells in ASCII digits, `None` for a text that is empty
/// or holds another character, and the largest number 64 bits hold for one past it. A text
/// of digits that holds one outside ASCII, and one of more than 4,300 digits, are refused.
pub fn codon_table_digits(text: &str, transcript: &Object) -> Result<Option<i64>, Failure> {
    let refused = |why: &str| {
        Failure::own(format!(
            "the cache record of {} holds a codon_table {why}",
            stable_id_text(transcript)
        ))
    };
    if text.is_empty() {
        return Ok(None);
    }
    if text.bytes().all(|byte| byte.is_ascii_digit()) {
        if text.len() > MAX_INT_TEXT_DIGITS {
            return Err(refused("of more than 4,300 digits"));
        }
        return Ok(Some(text.parse().unwrap_or(i64::MAX)));
    }
    let digit = |c: char| c.is_ascii_digit() || (!c.is_ascii() && c.is_numeric());
    if text.chars().all(digit) {
        return Err(refused("with a digit outside ASCII"));
    }
    Ok(None)
}

/// Each transcript of a corpus cache that carries a `translateable_seq`, by stable id, with
/// the sequences Ensembl VEP's Downstream plugin translates and the table it translates
/// them under (`CodingSequence`); a converted cache carries the table as text, a native one
/// as a number. A transcript whose table is no number is left out.
pub fn downstream_sequences(cache_dir: &str) -> Result<CodingSequences, Failure> {
    let mut sequences = HashMap::new();
    for item in iter_cache_transcripts(cache_dir) {
        let (_, transcript) = item?;
        let Some(cache) = effect_feature_cache(&transcript)? else {
            continue;
        };
        let coding = cache.get("translateable_seq").and_then(Value::as_str);
        let Some(coding) = coding.filter(|coding| !coding.is_empty()) else {
            continue;
        };
        let codon_table = match cache.get("codon_table").filter(|table| table.truthy()) {
            None => Some(1),
            Some(Value::Int(table)) if !table.is_negative() => {
                Some(table.as_i64().unwrap_or(i64::MAX))
            }
            Some(Value::Str(text)) => codon_table_digits(text, &transcript)?,
            Some(_) => None,
        };
        let Some(codon_table) = codon_table else {
            continue;
        };
        let Some(stable_id) = filed_stable_id(&transcript)? else {
            continue;
        };
        let utr = cache.get("three_prime_utr").and_then(Value::as_str);
        sequences.insert(
            stable_id.to_string(),
            CodingSequence {
                coding: coding.to_string(),
                utr: utr.unwrap_or_default().to_string(),
                codon_table: if codon_table == 0 { 1 } else { codon_table },
            },
        );
    }
    Ok(sequences)
}

/// Why a value of the cache gives no number.
enum Unread {
    /// `int()` raises an error that leaves the transcript out.
    Left,
    /// CPython reads a number the tool does not hold, or ends on a traceback.
    Outside,
}

/// `int(text)`: blanks around the number, a sign, ASCII digits with single underscores
/// between them, within 64 bits.
fn int_of_text(text: &str) -> Result<i64, Unread> {
    if let Some(number) = py_int(text) {
        return Ok(number);
    }
    let body = text.trim_matches(|c| is_py_space(c) && !('\u{1c}'..='\u{1f}').contains(&c));
    let body = body.strip_prefix(['-', '+']).unwrap_or(body);
    if body.is_empty() || body.starts_with('_') || body.ends_with('_') || body.contains("__") {
        return Err(Unread::Left);
    }
    let digits = body.chars().filter(|c| *c != '_');
    if digits.clone().all(|c| c.is_ascii_digit()) {
        // A number past 64 bits: `int()` reads one of up to 4,300 digits.
        return Err(if digits.count() > MAX_INT_TEXT_DIGITS {
            Unread::Left
        } else {
            Unread::Outside
        });
    }
    let digit = |c: char| c.is_ascii_digit() || (!c.is_ascii() && c.is_numeric());
    Err(if digits.clone().all(digit) {
        Unread::Outside
    } else {
        Unread::Left
    })
}

/// `int(value)` of a value a record holds under a key, `None` for a key it lacks.
fn int_of(value: Option<&Value>) -> Result<i64, Unread> {
    // The doubles from this bound up, and below its negative, are past 64 bits.
    const PAST_64_BITS: f64 = 9_223_372_036_854_775_808.0;
    match value {
        Some(Value::Int(number)) => number.as_i64().ok_or(Unread::Outside),
        Some(Value::Str(text)) => int_of_text(text),
        Some(Value::Bool(flag)) => Ok(i64::from(*flag)),
        Some(Value::Float(number)) if number.is_nan() => Err(Unread::Left),
        Some(Value::Float(number)) => {
            let whole = number.trunc();
            if (-PAST_64_BITS..PAST_64_BITS).contains(&whole) {
                Ok(whole as i64)
            } else {
                Err(Unread::Outside)
            }
        }
        Some(Value::Null | Value::List(_) | Value::Object(_)) | None => Err(Unread::Left),
    }
}

/// The coding model of a transcript record with the effect cache given.
fn coding_model(transcript: &Object, cache: Option<&Object>) -> Result<CodingModel, Unread> {
    let held = cache
        .and_then(|cache| cache.get("mapper"))
        .and_then(Value::as_object)
        .and_then(|mapper| mapper.get("pairs"));
    let items: &[Value] = match held {
        Some(Value::List(items)) => items,
        // An empty text and an empty object hold no pair; the characters of a text and
        // the keys of an object are no pairs.
        Some(pairs @ (Value::Str(_) | Value::Object(_))) if !pairs.truthy() => &[],
        _ => return Err(Unread::Left),
    };
    let mut pairs = Vec::with_capacity(items.len());
    for item in items {
        let pair = item.as_object().ok_or(Unread::Left)?;
        let mut values = [0i64; 5];
        for (value, key) in values.iter_mut().zip(MAPPER_PAIR_KEYS) {
            *value = int_of(pair.get(key))?;
        }
        pairs.push(values);
    }
    let coding_start = int_of(transcript.get("cdna_coding_start"))?;
    let coding_end = int_of(transcript.get("cdna_coding_end"))?;
    let peptide = cache
        .and_then(|cache| cache.get("peptide"))
        .and_then(Value::as_str);
    Ok(CodingModel {
        pairs,
        coding_start,
        coding_end,
        peptide_length: peptide.map(length),
    })
}

/// Each transcript of a corpus cache that carries a `cdna_coding_start`, a
/// `cdna_coding_end` and mapper pairs, by chromosome and stable id, as a `CodingModel`; a
/// converted cache carries the pairs' numbers as text. The chromosome is part of the key
/// because a cache files a pseudoautosomal transcript under X and under Y, with its own
/// pairs on each. A transcript lacking one of the three, or holding a value that is no
/// number, is left out; a number 64 bits do not hold is refused.
pub fn downstream_coding_models(cache_dir: &str) -> Result<CodingModels, Failure> {
    let mut models = HashMap::new();
    for item in iter_cache_transcripts(cache_dir) {
        let (chromosome, transcript) = item?;
        let cache = effect_feature_cache(&transcript)?;
        let model = match coding_model(&transcript, cache) {
            Ok(model) => model,
            Err(Unread::Left) => continue,
            Err(Unread::Outside) => {
                return Err(Failure::own(format!(
                    "the cache record of {} holds a coordinate that is no number of 64 bits in ASCII digits",
                    stable_id_text(&transcript)
                )))
            }
        };
        if let Some(Value::Str(stable_id)) = transcript.get("stable_id") {
            models.insert((chromosome, stable_id.clone()), model);
        }
    }
    Ok(models)
}

/// The sequences and coding models of the corpus cache's transcripts, each `None` in a run
/// that compares no field of the Downstream plugin.
#[derive(Default)]
pub struct Models {
    /// `downstream_sequences` of the corpus cache.
    pub coding_sequences: Option<CodingSequences>,
    /// `downstream_coding_models` of the corpus cache.
    pub coding_models: Option<CodingModels>,
}

impl Models {
    pub fn load(run: &FieldRun<'_>) -> Result<Models, Failure> {
        let compared = |name: &String| DOWNSTREAM_FIELDS.contains(&name.as_str());
        if !run.fields.iter().any(compared) {
            return Ok(Models::default());
        }
        Ok(Models {
            coding_sequences: Some(downstream_sequences(run.cache_dir)?),
            coding_models: Some(downstream_coding_models(run.cache_dir)?),
        })
    }
}

/// A run of ASCII digits as a number, `None` for any other text and for a number a `usize`
/// does not hold.
fn decimal(text: &str) -> Option<usize> {
    if text.is_empty() || !text.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    text.parse().ok()
}

/// The one coding base or the span of them a `CDS_position` column names, neither end a
/// gap (`?`): the first number, and the second where the column holds two.
fn cds_position_span(column: &str) -> Option<(usize, usize)> {
    match column.split_once('-') {
        Some((start, end)) => Some((decimal(start)?, decimal(end)?)),
        None => decimal(column).map(|base| (base, base)),
    }
}

/// A non-zero `HGVS_OFFSET`, the only value Ensembl VEP prints for the field; `None` for
/// any other text and for a number 128 bits do not hold.
fn non_zero_integer(text: &str) -> Option<i128> {
    let digits = text.strip_prefix('-').unwrap_or(text);
    let well_formed = !digits.is_empty()
        && !digits.starts_with('0')
        && digits.bytes().all(|byte| byte.is_ascii_digit());
    if !well_formed {
        return None;
    }
    text.parse().ok()
}

/// The reverse complement of `bases`; a character that is no upper-case base stays.
fn reverse_complement(bases: &str) -> String {
    let complement = |base: char| match base {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        other => other,
    };
    bases.chars().rev().map(complement).collect()
}

/// A field of a row, `None` for a question without the row and for a field the row lacks.
fn field<'a>(row: Option<&'a FieldRow>, name: &str) -> Option<&'a str> {
    row.and_then(|row| row.get(name)).map(String::as_str)
}

/// The class of a differing `DownstreamProtein` or `ProteinLengthChange` pair with the
/// class of a second cause that changes the value, or `None`; `index` is the field's place
/// in `DOWNSTREAM_FIELDS`.
fn downstream_class(
    question: &FieldQuestion<'_>,
    models: &Models,
    index: usize,
) -> Option<(&'static str, Option<&'static str>)> {
    let key = question.key;
    let mut terms = key.consequences.split(',');
    if !terms.clone().any(|term| term == "frameshift_variant")
        || terms.any(|term| term.contains("splice"))
    {
        return None;
    }
    let sequence = models.coding_sequences.as_ref()?.get(&key.feature)?;
    genetic_code(sequence.codon_table)?;
    let column = field(question.perl_row, "CDS_position").unwrap_or_default();
    let (mut start, mut end) = cds_position_span(column)?;
    if question.vep == "-" || question.vep_rs == "-" {
        return None;
    }
    let chromosome = key
        .location
        .rsplit_once(':')
        .map_or(key.location.as_str(), |(chromosome, _)| chromosome);
    let model = models
        .coding_models
        .as_ref()
        .and_then(|held| held.get(&(chromosome.to_string(), key.feature.clone())));
    let peptide_length = model.and_then(|model| model.peptide_length);
    if field(question.rust_row, "CDS_position").unwrap_or(column) != column
        || (index == 1 && peptide_length.is_none())
    {
        return None;
    }
    let reverse = question.strand == "-1";
    let padding = sequence
        .coding
        .chars()
        .take_while(|base| matches!(base, 'N' | 'n'))
        .count();
    let written = if key.allele == "-" {
        ""
    } else {
        key.allele.as_str()
    };
    let on_strand = |bases: &str| {
        if reverse {
            reverse_complement(bases)
        } else {
            bases.to_string()
        }
    };
    let span = question
        .genomic_span
        .map(|(start, end)| (i128::from(start), i128::from(end)));
    let insertion = span.is_some_and(|(start, end)| start == end + 1);
    if insertion && start.checked_add(1) == Some(end) {
        let carried = model
            .zip(span)
            .and_then(|(model, (start, end))| model.cds_pair(start, end, reverse, padding));
        if carried != Some((signed(start), signed(end))) {
            (start, end) = (end, start);
        }
    }
    if start < 1 || end > length(&sequence.coding) || start > end + 1 {
        return None;
    }
    let own = (start, end, on_strand(written));
    let mut shifted = own.clone();
    let offset = field(question.perl_row, "HGVS_OFFSET").and_then(non_zero_integer);
    if let (Some(model), Some((genomic_start, genomic_end)), Some(offset)) = (model, span, offset) {
        let moved = genomic_start
            .checked_add(offset)
            .zip(genomic_end.checked_add(offset))
            .and_then(|(start, end)| model.cds_pair(start, end, reverse, padding));
        if let Some((moved_start, moved_end)) = moved {
            // A start a `usize` does not hold lies past the coding sequence, where the
            // plugin's call ends; an end it does not hold is any end past the sequence.
            shifted = (
                usize::try_from(moved_start).ok()?,
                usize::try_from(moved_end).unwrap_or(usize::MAX),
                on_strand(&hgvs_shifted_allele(written, offset, reverse)),
            );
        }
    }
    let cached_length = peptide_length.unwrap_or(0);
    let (coding, utr, table) = (&sequence.coding, &sequence.utr, sequence.codon_table);
    let plugin_values = |(start, end, bases): &(usize, usize, String)| {
        downstream_plugin_values(coding, utr, table, cached_length, *start, *end, bases)
    };
    let (reference_own, past_coding_end) = plugin_values(&own)?;
    let (reference, _) = plugin_values(&shifted)?;
    let own_values =
        downstream_own_values(coding, utr, table, cached_length, own.0, own.1, &own.2)?;
    if question.vep != reference[index] || question.vep_rs != own_values[index] {
        return None;
    }
    let own_class = if past_coding_end {
        DOWNSTREAM_UTR_FRAME_PAST_CODING_END
    } else {
        DOWNSTREAM_START_PAST_OWN_CODON
    };
    if reference[index] != reference_own[index] {
        let both_causes = index == 0 && question.vep_rs != reference_own[index];
        return Some((
            PLUGIN_HGVS_SHIFTED_CDS_COORDINATES,
            both_causes.then_some(own_class),
        ));
    }
    if index == 1 || question.vep == question.vep_rs {
        return None;
    }
    Some((own_class, None))
}

/// The class of a differing `DownstreamProtein` or `ProteinLengthChange` pair with the
/// class of a second cause that changes the value; any other field passes, and a field of
/// the plugin that takes no class is asked of no further class. The row carries
/// `frameshift_variant` and no term naming a splice site, the rows the plugin annotates
/// (Downstream.pm:95-99); both engines wrote a value and print one `CDS_position` with no
/// gap; its transcript is one of the cache's coding transcripts (`downstream_sequences`),
/// and the values are read under the transcript's codon table, so a transcript under a
/// table that has no code gets no class.
///
/// The row has two readings. Its own is the `CDS_position` column with the row's allele on
/// the transcript's strand, none for `-`: one base or `start-end` in order, and for an
/// insertion (the record's genomic span of `row_record_span` has its start one past its
/// end) the two flanking bases the column prints lower first, read as a start one past its
/// end unless the transcript's model (`downstream_coding_models`) carries exactly those
/// bases at the insertion's point. The shifted one is what the plugin reads under `--hgvs`
/// on a row whose `HGVS_OFFSET` in Ensembl VEP's row is a non-zero integer: the CDS pair of
/// the record's span moved by that offset along the genome (`CodingModel::cds_pair`) with
/// the allele of `hgvs_shifted_allele`, where both ends of the moved span are coding bases
/// (TranscriptVariationAllele.pm:1730-1738 returns before the move otherwise); on any other
/// row, and without a span or a model, it is the own reading.
///
/// vep-rs's value must be `downstream_own_values` on the own reading and Ensembl VEP's
/// `downstream_plugin_values` on the shifted one. Where Ensembl VEP's rule gives another
/// value on the own reading, the flag changes its value and the class is
/// `plugin_hgvs_shifted_cds_coordinates`, for either field; a `ProteinLengthChange` needs
/// the cached peptide's length and has no other class. Otherwise the two rules differ on
/// the own reading itself, a span on the coding sequence whose first base is a codon's
/// third: `downstream_utr_frame_past_coding_end` where Ensembl VEP's translation starts
/// past the changed coding sequence, else `downstream_start_past_own_codon`. A
/// `DownstreamProtein` the flag changes whose vep-rs value is not Ensembl VEP's on the own
/// reading either is changed by both causes: its second class is the one its own reading
/// has.
pub fn downstream_field_class(question: &FieldQuestion<'_>, models: &FieldModels) -> Verdict {
    let Some(index) = DOWNSTREAM_FIELDS
        .iter()
        .position(|field| *field == question.name)
    else {
        return Verdict::Pass;
    };
    match downstream_class(question, &models.downstream, index) {
        Some((class, second)) => Verdict::Class(class, second),
        None => Verdict::NoClass,
    }
}

#[cfg(test)]
#[path = "downstream_tests.rs"]
mod tests;
