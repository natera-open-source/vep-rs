// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Peptide paths, the HGNC identifier copy, the Downstream and NMD readings and the shared shifts.

use crate::argparse::{path_text, repr};
use crate::cache;
use crate::classes::DOWNSTREAM_START_PAST_OWN_CODON;
use crate::classes::DOWNSTREAM_UTR_FRAME_PAST_CODING_END;
use crate::classes::PLUGIN_HGVS_SHIFTED_CDS_COORDINATES;
use crate::classes::START_COOCCURRENCE_SWAP;
use crate::consts::AMINO_ACIDS_FIELDS;
use crate::consts::CDS_POSITION_FIELDS;
use crate::consts::HGVS_OFFSET_FIELDS;
use crate::consts::INTEGER_RE_PATTERN;
use crate::consts::NON_ZERO_INTEGER_RE_PATTERN;
use crate::consts::PROTEIN_POSITION_FIELDS;
use crate::consts::SHIFT_LENGTH_FIELDS;
use crate::consts::STRAND_FIELDS;
use crate::entry::Entry;
use crate::error::Error;
use crate::error::Result;
use crate::fields::field_named;
use crate::genetic_codes::{genetic_code, CodonTable, GeneticCode, Reading};
use crate::pyfmt;
use crate::pyjson;
use crate::pyjson::Value;
use crate::rules::coding_end_pairs::{coding_end_class, AtACodingEnd};
use crate::rules::intron_ends::VERSION_SUFFIX_RE;
use crate::rules::kept_start::start_codon_kept_without_five_prime_utr_pair;
use crate::rules::record_alleles::earlier_allele_shift;
use crate::rules::record_alleles::first_allele_shift;
use crate::rules::record_alleles::later_alternate;
use crate::rules::registry::registry_class;
use crate::rules::term_pairs::stop_retained_without_stop_codon_test_pair;
use crate::rules::transcript_model::int_of_text;
use crate::rules::transcript_model::TranscriptModel;
use crate::rules::whole_text;
use crate::self_parity::reverse_complement;
use crate::shared::blank;
use crate::shared::consequence_set;
use crate::shared::feature_kind;
use crate::shared::integer_text_of;
use crate::shared::ABSENT;
use crate::shared::FRAMESHIFT_TERM;
use indexmap::map::Entry as Held;
use indexmap::IndexMap;
use std::borrow::Cow;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::path::Path;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::LazyLock;

/// The start and the end of a CDS span, a gap `None`.
pub type CdsPair = (Option<i64>, Option<i64>);

/// One reading of a row the notation shifted: the own CDS pair, the pair a plugin reads under
/// --hgvs where the shift reaches it, and the record's genomic end.
pub type ShiftedReading = (CdsPair, Option<(i64, i64)>, i64);

static WHOLE_INTEGER_RE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text(INTEGER_RE_PATTERN));
static WHOLE_NON_ZERO_INTEGER_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(NON_ZERO_INTEGER_RE_PATTERN));
static WHOLE_HGNC_ID_RE: LazyLock<regex::Regex> = LazyLock::new(|| whole_text(HGNC_ID_RE_PATTERN));
static WHOLE_PEPTIDE_PAIR_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(PEPTIDE_PAIR_RE_PATTERN));

/// Whether `text` is one or more of `letters`.
fn spelt_with(text: &str, letters: &str) -> bool {
    !text.is_empty() && text.chars().all(|letter| letters.contains(letter))
}

/// `int(text)` of a text of digits, a number 64 bits do not hold at the nearest one they do.
fn clamped_int(text: &str) -> i64 {
    int_of_text(text).unwrap_or(if text.starts_with('-') {
        i64::MIN
    } else {
        i64::MAX
    })
}

/// The moves of a later allele's row in the order they are tried: the row's own offset,
/// then an earlier allele's where it is another; an offset of zero moves nothing.
fn distinct_moves(shift: i64, earlier: Option<i64>) -> Vec<i64> {
    let mut moves = Vec::with_capacity(2);
    if shift != 0 {
        moves.push(shift);
    }
    if let Some(earlier) = earlier.filter(|earlier| *earlier != 0 && *earlier != shift) {
        moves.push(earlier);
    }
    moves
}

/// The fields that state a row's coding change beside its consequence set; a row of
/// `organelle_peptide_paths_translation_table` carries each alike on both sides.
pub const CODING_CHANGE_FIELDS: &[&str] = &[
    "cDNA_position",
    "cdna_start",
    "cdna_end",
    "CDS_position",
    "cds_start",
    "cds_end",
    "Protein_position",
    "protein_start",
    "protein_end",
    "Amino_acids",
    "amino_acids",
    "Codons",
    "codons",
    "HGVSc",
    "hgvsc",
    "HGVS_OFFSET",
    "hgvs_offset",
];

pub const SYMBOL_SOURCE_FIELDS: &[&str] = &["SYMBOL_SOURCE", "gene_symbol_source"];

pub const HGNC_SOURCE: &str = "HGNC";

pub const HGNC_ID_RE_PATTERN: &str = r"(?:HGNC:)?\d+";
pub static HGNC_ID_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(HGNC_ID_RE_PATTERN).expect("the pattern compiles"));

pub const PROTEIN_LENGTH_CHANGE_FIELDS: &[&str] = &["ProteinLengthChange", "proteinlengthchange"];

pub const NMD_ESCAPING: &str = "NMD_escaping_variant";

pub const NMD_GATE_TERMS: &[&str] = &[
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
];

/// The bases each nucleotide letter stands for.
pub const NUCLEOTIDE_READINGS: &[(&str, &str)] = &[
    ("A", "A"),
    ("C", "C"),
    ("G", "G"),
    ("T", "T"),
    ("U", "T"),
    ("R", "AG"),
    ("Y", "CT"),
    ("S", "CG"),
    ("W", "AT"),
    ("K", "GT"),
    ("M", "AC"),
    ("B", "CGT"),
    ("D", "AGT"),
    ("H", "ACT"),
    ("V", "ACG"),
    ("N", "ACGT"),
    ("X", "ACGT"),
];

/// The letter Bio::Tools::CodonTable writes for a codon whose readings code two amino acids.
pub const AMBIGUOUS_RESIDUE_PAIRS: &[(&[&str], &str)] = &[(&["D", "N"], "B"), (&["E", "Q"], "Z")];

pub const BLOSUM62_FIELDS: &[&str] = &["BLOSUM62"];

pub const BLOSUM62_RESIDUES: &str = "ARNDCQEGHILKMFPSTWYV";

pub const BLOSUM62_ROWS: &[&str] = &[
    "4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0",
    "-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3",
    "-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3",
    "-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3",
    "0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1",
    "-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2",
    "-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2",
    "0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3",
    "-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3",
    "-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3",
    "-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1",
    "-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2",
    "-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1",
    "-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1",
    "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2",
    "1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2",
    "0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0",
    "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3",
    "-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1",
    "0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4",
];

/// The plugin's score of each pair of residues, as it prints it.
pub static BLOSUM62_MATRIX: LazyLock<BTreeMap<(char, char), &'static str>> = LazyLock::new(|| {
    let mut matrix = BTreeMap::new();
    for (first, row) in BLOSUM62_RESIDUES.chars().zip(BLOSUM62_ROWS) {
        let cells = pyfmt::split_whitespace(row);
        for (second, cell) in BLOSUM62_RESIDUES.chars().zip(cells) {
            matrix.insert((first, second), cell);
        }
    }
    matrix
});

pub const PEPTIDE_PAIR_RE_PATTERN: &str = r"([A-Z])/([A-Z])";
pub static PEPTIDE_PAIR_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(PEPTIDE_PAIR_RE_PATTERN).expect("the pattern compiles"));

/// `dup`, `ins` or `other` for the text after `p.`; a `delins` is `other`.
pub fn peptide_description_shape(description: &str) -> String {
    if description.ends_with("dup") {
        return "dup".to_string();
    }
    let inserted = description
        .match_indices("ins")
        .any(|(at, _)| !description[..at].ends_with("del"));
    if inserted { "ins" } else { "other" }.to_string()
}

/// The fields of `CODING_CHANGE_FIELDS` an entry carries a value for, under their names in
/// lower case.
fn coding_change(entry: &Entry) -> BTreeMap<String, &str> {
    entry
        .fields
        .iter()
        .filter(|(name, value)| field_named(name, CODING_CHANGE_FIELDS) && !blank(value))
        .map(|(name, value)| (pyfmt::lower(name), value.as_str()))
        .collect()
}

/// Whether the whole number `digits` writes, without a sign or a leading zero, is more
/// than `bound`.
fn decimal_exceeds(digits: &str, bound: u128) -> bool {
    digits.len() > 38 || digits.parse::<u128>().is_ok_and(|number| number > bound)
}

/// The whole number `digits` writes modulo `modulus`.
fn decimal_modulo(digits: &str, modulus: u128) -> u128 {
    digits.bytes().fold(0, |rest, digit| {
        (rest * 10 + u128::from(digit - b'0')) % modulus
    })
}

/// Whether a differing HGVSp pair names one protein and differs in a path of
/// TranscriptVariationAllele.pm that translates a sequence bare, on a row the bare translation
/// alone can explain: `_get_fs_peptides` and `_stop_loss_extra_AA` write the whole
/// description of a frameshift (its residues, an immediate stop, the stop distance and the
/// lost-start shapes), `_stop_loss_extra_AA` the count of an extension after a lost stop, and
/// `_check_for_peptide_duplication` the choice between a duplication and an insertion. A
/// substitution, a `delins` and a plain deletion read their residues through `peptide`,
/// which follows the transcript's table, and a value one side leaves empty names no protein,
/// so each stays a difference. The class changes the HGVSp alone: the reference entry `e`
/// and the candidate's `a` state one coding change, the same consequence set and the same
/// `CODING_CHANGE_FIELDS`, so a row whose terms, positions, codons, coding notation or shift
/// differ (a later allele annotated at another allele's shifted span, an indel shifted short at
/// a sequence edge) is left to the class of that difference.
/// An inserted allele of two or more bases on a reverse-strand transcript, shifted past its
/// own length by a length that is no multiple of it, stays a difference too:
/// TranscriptVariationAllele.pm:1337-1339 rotate a reverse-strand transcript's allele by its
/// length less the shift, nothing once the shift is the longer, so the HGVSp is read from
/// bases the record does not produce. A row that does not state a forward strand is read as
/// a reverse-strand one.
pub fn organelle_peptide_path(
    e: &Entry,
    a: &Entry,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    let (expected_name, _, expected_description) = pyfmt::rpartition(expected_value, ":p.");
    let (actual_name, _, actual_description) = pyfmt::rpartition(actual_value, ":p.");
    if expected_name.is_empty() || expected_name != actual_name {
        return false;
    }
    let terms = consequence_set(&e.consequence);
    if terms != consequence_set(&a.consequence) || coding_change(e) != coding_change(a) {
        return false;
    }
    let offset = field_of(Some(e), HGVS_OFFSET_FIELDS);
    if pyfmt::len(&e.allele) >= 2
        && spelt_with(&e.allele, "ACGTN")
        && WHOLE_INTEGER_RE.is_match(&offset)
    {
        let shift = match pyfmt::int_checked(&offset) {
            Ok(number) => number,
            Err(error) => crate::uncaught!("{error}"),
        };
        let shift = shift.trim_start_matches('-');
        let length = e.allele.len() as u128;
        if decimal_exceeds(shift, length)
            && decimal_modulo(shift, length) != 0
            && field_of(Some(e), STRAND_FIELDS) != "1"
        {
            return false;
        }
    }
    let terms: BTreeSet<&str> = terms.split(',').collect();
    if terms.contains("frameshift_variant") {
        return true;
    }
    if terms.contains("stop_lost") {
        return expected_description.contains("ext") && actual_description.contains("ext");
    }
    let shapes = BTreeSet::from([
        peptide_description_shape(expected_description),
        peptide_description_shape(actual_description),
    ]);
    (terms.contains("inframe_insertion") || terms.contains("protein_altering_variant"))
        && shapes == BTreeSet::from(["dup".to_string(), "ins".to_string()])
}

/// Whether the HGNC_ID `expected_value` of the reference entry `e` and `actual_value`
/// of the candidate entry `a` have the shape of `hgnc_id_copy`: a Transcript entry, one
/// value empty and the other an HGNC id, and the SYMBOL_SOURCE either side prints a source
/// other than HGNC. A gene HGNC names carries its own id on both engines, so a differing id
/// there, and two differing ids anywhere, stay differences.
pub fn hgnc_id_copy_pair(
    e: &Entry,
    a: Option<&Entry>,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    let last = e.key.rsplit('|').next().unwrap_or(&e.key);
    let kind = last.split('#').next().unwrap_or(last);
    if feature_kind(kind) != "transcript" {
        return false;
    }
    if blank(expected_value) == blank(actual_value) {
        return false;
    }
    let filled = if blank(expected_value) {
        actual_value
    } else {
        expected_value
    };
    if !WHOLE_HGNC_ID_RE.is_match(filled) {
        return false;
    }
    let sources: BTreeSet<String> = [Some(e), a]
        .into_iter()
        .flatten()
        .map(|entry| field_of(Some(entry), SYMBOL_SOURCE_FIELDS))
        .filter(|source| source != ABSENT)
        .collect();
    !sources.is_empty()
        && sources
            .iter()
            .all(|source| !blank(source) && source != HGNC_SOURCE)
}

/// The genetic code a cache transcript's `codon_table` names: a number in a native cache,
/// the dump's string in a converted one, table 1 for 0 or no key (Downstream.pm:132); a
/// number no table is held for (`genetic_codes::held_tables`) and any other value are a
/// code the comparison lacks, kept as the cache spells it.
pub fn codon_table_id(value: &Value) -> CodonTable {
    if matches!(value, Value::Null) {
        return CodonTable::Held(1);
    }
    let text = pyjson::py_str(value);
    let table = match pyfmt::is_decimal(&text) {
        true => int_of_text(&text),
        false => None,
    };
    match table {
        Some(0) => CodonTable::Held(1),
        Some(table) if genetic_code(table, Reading::VepRs).is_some() => CodonTable::Held(table),
        _ => CodonTable::Lacking(text),
    }
}

/// What the Downstream and NMD plugins read of one cache transcript: the
/// `translateable_seq`, the `three_prime_utr` (empty where the cache spells none), the
/// codon table (`codon_table_id`, a code no table is held for kept as the cache spells it)
/// and the `peptide` with its length (each `None` without one) of its
/// `variation_effect_feature_cache`, and its `TranscriptModel` (`None` where the cache
/// does not state one).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DownstreamTranscript {
    pub coding: String,
    pub utr: String,
    pub table: CodonTable,
    pub peptide_length: Option<i64>,
    pub peptide: Option<String>,
    pub model: Option<TranscriptModel>,
}

/// A count of code points or of residues as a whole number.
fn count(length: usize) -> i64 {
    i64::try_from(length).unwrap_or(i64::MAX)
}

impl DownstreamTranscript {
    /// The Downstream plugin's value on the CDS pair `start` and `end` (an insertion's
    /// start one past its end) with the bases `allele` on the transcript's strand: the
    /// DownstreamProtein, or with `length_change` the ProteinLengthChange as printed. Both
    /// rules write the allele over the pair (Downstream.pm:108-110) and translate the changed
    /// coding sequence with the 3' UTR up to the first stop; the ProteinLengthChange is the
    /// whole translation's length less the peptide's (:141-148; the translated coding
    /// sequence's for a transcript the cache gives no peptide). With `ensembl` the
    /// DownstreamProtein starts at the pair's lower coordinate rounded down to a multiple of
    /// three (:112-115) and there is no value for a start past the coding sequence's end,
    /// where the assignment of :110 fails; without it, the candidate's rule, at the codon
    /// holding the first changed base. Each value is translated as its engine reads the
    /// transcript's code (`Reading`): the plugin hands the cache's id to the module (:130-139).
    pub fn downstream_value(
        &self,
        start: i64,
        end: i64,
        allele: &str,
        length_change: bool,
        ensembl: bool,
    ) -> Option<String> {
        let coding_length = count(pyfmt::len(&self.coding));
        let before = start.saturating_sub(1);
        if ensembl && before > coding_length {
            return None;
        }
        // A transcript under a code no table is held for has no value: its callers end the
        // run on it by name before they ask.
        let table = match &self.table {
            CodonTable::Held(table) => *table,
            CodonTable::Lacking(spelt) => crate::uncaught!("{spelt}"),
        };
        let reading = match ensembl {
            true => Reading::EnsemblVep,
            false => Reading::VepRs,
        };
        let mut changed = String::with_capacity(self.coding.len() + allele.len());
        changed.push_str(pyfmt::slice(&self.coding, None, Some(before)));
        changed.push_str(allele);
        changed.push_str(pyfmt::slice(&self.coding, Some(end.max(before)), None));
        if length_change {
            let reference = match self.peptide_length {
                Some(length) => length,
                None => count(translate_until_stop(&self.coding, table, reading).len()),
            };
            changed.push_str(&self.utr);
            let translated = count(translate_until_stop(&changed, table, reading).len());
            return Some((i128::from(translated) - i128::from(reference)).to_string());
        }
        let first = if ensembl {
            start.min(end)
        } else {
            before.min(coding_length)
        };
        let from = first.saturating_sub(first.rem_euclid(3));
        let mut read = pyfmt::slice(&changed, Some(from), None).to_string();
        read.push_str(&self.utr);
        Some(translate_until_stop(&read, table, reading))
    }

    /// Whether the two readings differ in nothing but their models.
    fn alike_but_for_the_model(&self, other: &DownstreamTranscript) -> bool {
        self.coding == other.coding
            && self.utr == other.utr
            && self.table == other.table
            && self.peptide_length == other.peptide_length
            && self.peptide == other.peptide
    }
}

/// What `DownstreamSequences` keeps of one cache transcript object: its stable id and its
/// reading; `None` for a transcript without a stable id or a `translateable_seq`, or whose
/// `three_prime_utr` is no text.
fn downstream_transcript(transcript: &Value) -> Option<(String, DownstreamTranscript)> {
    let stable_id = transcript.get("stable_id")?.as_str()?;
    let cached = transcript.get("variation_effect_feature_cache")?;
    let members = cached.as_dict()?;
    let coding = members.get("translateable_seq")?.as_str()?;
    if coding.is_empty() {
        return None;
    }
    let utr = match members.get("three_prime_utr") {
        None | Some(Value::Null) => "",
        Some(Value::Str(utr)) => utr,
        Some(_) => return None,
    };
    let coding = pyfmt::upper(coding);
    let null = Value::Null;
    let table = codon_table_id(members.get("codon_table").unwrap_or(&null));
    let peptide = members.get("peptide").and_then(Value::as_str);
    let peptide_length = peptide.map(|peptide| count(pyfmt::len(peptide)));
    let padding = count(coding.chars().take_while(|base| *base == 'N').count());
    let model = TranscriptModel::from_cache(transcript, padding);
    Some((
        stable_id.to_string(),
        DownstreamTranscript {
            coding,
            utr: pyfmt::upper(utr),
            table,
            peptide_length,
            peptide: peptide.map(str::to_string),
            model,
        },
    ))
}

/// The letters of `sequence`, one byte a code point: a code point outside ASCII as a byte
/// no nucleotide letter is.
fn letters_of(sequence: &str) -> Cow<'_, [u8]> {
    if sequence.is_ascii() {
        return Cow::Borrowed(sequence.as_bytes());
    }
    Cow::Owned(
        sequence
            .chars()
            .map(|letter| {
                u8::try_from(letter)
                    .ok()
                    .filter(u8::is_ascii)
                    .unwrap_or(0x80)
            })
            .collect(),
    )
}

/// The letter of a codon `code` does not hold: the one amino acid its readings code, the
/// letter of `AMBIGUOUS_RESIDUE_PAIRS` for two, else `X`.
pub(crate) fn ambiguous_residue(code: GeneticCode, codon: &[u8]) -> u8 {
    let readings: Vec<&[u8]> = codon
        .iter()
        .map(|base| {
            NUCLEOTIDE_READINGS
                .iter()
                .find(|(letter, _)| letter.as_bytes() == [*base])
                .map_or(&b""[..], |(_, bases)| bases.as_bytes())
        })
        .collect();
    let mut coded = BTreeSet::new();
    for first in readings[0] {
        for second in readings[1] {
            for third in readings[2] {
                coded.extend(code.amino_acid(&[*first, *second, *third]));
            }
        }
    }
    if coded.len() == 1 {
        return coded.pop_first().unwrap_or(b'X');
    }
    AMBIGUOUS_RESIDUE_PAIRS
        .iter()
        .find(|(pair, _)| {
            pair.iter()
                .map(|residue| residue.as_bytes()[0])
                .collect::<BTreeSet<u8>>()
                == coded
        })
        .map_or(b'X', |(_, letter)| letter.as_bytes()[0])
}

/// The transcripts of a candidate JSON cache as the Downstream and NMD plugins read them
/// (`DownstreamTranscript`), by transcript stable id, kept for a transcript with a
/// `translateable_seq`. A chromosome's shards are read on the first lookup of an entry
/// whose Location names its directory, every shard on the first lookup of an entry without
/// one (a VCF or JSON entry). A stable id the shards read give one reading under two models
/// (a transcript the cache files under two chromosomes, each with its own mapper pairs) is
/// kept without a model; one they give two different readings otherwise is left out.
///
/// `_chromosomes` holds the directories whose names are UTF-8: a Location names no other.
#[derive(Debug)]
pub struct DownstreamSequences {
    pub root: PathBuf,
    pub _chromosomes: BTreeSet<String>,
    pub _sequences: RefCell<IndexMap<Option<String>, DownstreamReadings>>,
}

/// The readings of one chromosome's shards, or of every shard, by transcript stable id.
pub type DownstreamReadings = IndexMap<String, Arc<DownstreamTranscript>>;

/// Ends the run on a directory that cannot be listed, in the system's words for it.
fn unlisted(directory: &Path, error: &std::io::Error) -> ! {
    crate::uncaught!(
        "[Errno {}] {}: {}",
        error.raw_os_error().unwrap_or(0),
        crate::report::strerror(error),
        repr(&path_text(directory))
    )
}

impl DownstreamSequences {
    pub fn new(cache: &Path) -> Result<Self> {
        let root = cache::joined(cache, "transcripts");
        if !root.is_dir() {
            return Err(Error::comparator(format!(
                "--json-cache: {} has no transcripts/ directory",
                path_text(cache)
            )));
        }
        let entries = match std::fs::read_dir(&root) {
            Ok(entries) => entries,
            Err(error) => unlisted(&root, &error),
        };
        let mut chromosomes = BTreeSet::new();
        for entry in entries {
            let entry = match entry {
                Ok(entry) => entry,
                Err(error) => unlisted(&root, &error),
            };
            if entry.path().is_dir() {
                chromosomes.extend(entry.file_name().into_string());
            }
        }
        Ok(DownstreamSequences {
            root,
            _chromosomes: chromosomes,
            _sequences: RefCell::new(IndexMap::new()),
        })
    }

    pub fn _read(&self, chromosome: Option<&str>) -> Result<DownstreamReadings> {
        let shards = cache::shards(&self.root, chromosome)?;
        let read = cache::map_shards(&shards, cache::threads(), |_, transcripts| {
            Ok(transcripts
                .iter()
                .filter_map(downstream_transcript)
                .collect::<Vec<_>>())
        })?;
        let mut sequences = DownstreamReadings::new();
        let mut conflicting = BTreeSet::new();
        for (stable_id, read) in read.into_iter().flatten() {
            match sequences.entry(stable_id) {
                Held::Vacant(place) => {
                    place.insert(Arc::new(read));
                }
                Held::Occupied(mut held) => {
                    if **held.get() == read {
                        continue;
                    }
                    if held.get().alike_but_for_the_model(&read) {
                        held.insert(Arc::new(DownstreamTranscript {
                            model: None,
                            ..read
                        }));
                    } else {
                        conflicting.insert(held.key().clone());
                    }
                }
            }
        }
        for stable_id in &conflicting {
            sequences.shift_remove(stable_id);
        }
        Ok(sequences)
    }

    /// The entry's transcript, or `None` for one the cache does not hold with a
    /// `translateable_seq`; a Feature carrying a version the cache's stable id lacks is
    /// looked up without it.
    pub fn transcript(&self, e: &Entry) -> Result<Option<Arc<DownstreamTranscript>>> {
        let chromosome = e
            .location
            .as_deref()
            .filter(|location| !location.is_empty())
            .map(|location| pyfmt::rpartition(location, ":").0)
            .filter(|name| self._chromosomes.contains(*name))
            .map(str::to_string);
        if !self._sequences.borrow().contains_key(&chromosome) {
            let read = self._read(chromosome.as_deref())?;
            self._sequences
                .borrow_mut()
                .insert(chromosome.clone(), read);
        }
        let sequences = self._sequences.borrow();
        let sequences = &sequences[&chromosome];
        if let Some(held) = sequences.get(&e.feature) {
            return Ok(Some(Arc::clone(held)));
        }
        let unversioned = VERSION_SUFFIX_RE.replace_all(&e.feature, "");
        Ok(sequences.get(unversioned.as_ref()).cloned())
    }
}

/// The translation of `sequence` under the genetic code `table` as `reading` reads it (a
/// code a table is held for) up to its first stop, as `Bio::Tools::CodonTable::translate`
/// reads it for Downstream.pm:139-142: the whole codons in order, a trailing one or two bases
/// dropped; a codon with an ambiguous letter by the amino acids its readings code under the
/// table, one amino acid as itself, aspartate with asparagine `B`, glutamate with glutamine
/// `Z` and anything else `X`, as is a codon with a letter that is no nucleotide's.
pub fn translate_until_stop(sequence: &str, table: i64, reading: Reading) -> String {
    let Some(code) = genetic_code(table, reading) else {
        crate::uncaught!("{table}")
    };
    let letters = letters_of(sequence);
    let mut residues = String::with_capacity(letters.len() / 3);
    for codon in letters.chunks_exact(3) {
        let residue = code
            .amino_acid(codon)
            .unwrap_or_else(|| ambiguous_residue(code, codon));
        if residue == b'*' {
            break;
        }
        residues.push(char::from(residue));
    }
    residues
}

/// The two bounds a position column or its two JSON keys spell, each read up to a
/// `/total` suffix; one number of a column stands for both.
pub(crate) fn bound_texts<'a>(e: &'a Entry, column: &str, keys: [&str; 2]) -> [&'a str; 2] {
    let up_to_total = |text: &'a str| pyfmt::partition(text, "/").0;
    match e.fields.get(column) {
        None => keys.map(|key| up_to_total(e.fields.get(key).map_or(ABSENT, String::as_str))),
        Some(position) => {
            let (first, _, last) = pyfmt::partition(up_to_total(position), "-");
            [first, if last.is_empty() { first } else { last }]
        }
    }
}

/// The number a bound of digits writes, `None` for any other bound. A digit `int` does
/// not read ends the run as its `ValueError` does; a number 64 bits do not hold is read as
/// no number.
pub(crate) fn bound_number(bound: &str) -> Option<i64> {
    if pyfmt::is_digit(bound) {
        int_of_text(bound)
    } else {
        None
    }
}

/// The `cds_start` and `cds_end` an entry's coding span may spell, a gap (`?`, or a
/// JSON bound the entry lacks) `None`, and no pair for an entry without a coding
/// coordinate: the text formats' `CDS_position` or the JSON `cds_start` and `cds_end`
/// as given and, for a sequence allele at two adjacent numbers, swapped as well, the pair of
/// an insertion, whose start is one past its end and which every format spells lower first
/// (`Utils::format_coords`; the JSON entry's two keys likewise). Each bound is read up to
/// the `/total` suffix `--total_length` appends.
pub fn cds_pair_readings(e: &Entry) -> Vec<CdsPair> {
    let [first, last] = bound_texts(e, "CDS_position", ["cds_start", "cds_end"]);
    let start = bound_number(first);
    let end = bound_number(last);
    if start.is_none() && end.is_none() {
        return Vec::new();
    }
    let mut readings = vec![(start, end)];
    if let (Some(low), Some(high)) = (start, end) {
        if low.checked_add(1) == Some(high) && spelt_with(&e.allele, "ACGTNacgtn") {
            readings.push((end, start));
        }
    }
    readings
}

/// The genomic start and end a text entry's Location may spell: the span as printed and,
/// for a sequence allele at two adjacent bases, the point between them as well, an insertion's
/// span, whose start is one past its end.
pub fn record_span_readings(e: &Entry) -> Vec<(i64, i64)> {
    let location = e.location.as_deref().unwrap_or("");
    let (first, _, last) = pyfmt::partition(pyfmt::rpartition(location, ":").2, "-");
    let last = if last.is_empty() { first } else { last };
    if !(pyfmt::is_digit(first) && pyfmt::is_digit(last)) {
        return Vec::new();
    }
    let start = int_of_text(first);
    let end = int_of_text(last);
    let (Some(start), Some(end)) = (start, end) else {
        return Vec::new();
    };
    let mut readings = vec![(start, end)];
    if start.checked_add(1) == Some(end) && spelt_with(&e.allele, "ACGTNacgtn") {
        readings.push((end, start));
    }
    readings
}

/// The readings of a row whose reference entry `e` carries the HGVS_OFFSET `shift`,
/// each the own CDS pair (one of `own_pairs`), the pair a plugin reads under --hgvs where
/// the shift reaches it (the record's span moved by `shift` has a CDS coordinate at both
/// ends), else `None`, and the record's genomic end. The record's span is a reading of the
/// Location that maps to an own pair or, for an entry without a Location, the one span an own
/// pair maps back to (`TranscriptModel.genomic_span`).
pub fn shifted_pair_readings(
    e: &Entry,
    model: &TranscriptModel,
    shift: i64,
    own_pairs: &[CdsPair],
) -> Vec<ShiftedReading> {
    let has_location = e
        .location
        .as_deref()
        .is_some_and(|location| !location.is_empty());
    let spans: Vec<((i64, i64), CdsPair)> = if has_location {
        record_span_readings(e)
            .into_iter()
            .map(|span| (span, model.cds_pair(span.0, span.1)))
            .filter(|(_, own)| own_pairs.contains(own))
            .collect()
    } else {
        own_pairs
            .iter()
            .filter_map(|own| match own {
                (Some(start), Some(end)) => Some((model.genomic_span(*start, *end)?, *own)),
                _ => None,
            })
            .collect()
    };
    spans
        .into_iter()
        .map(|((start, end), own)| {
            let moved = |base: i64| i128::from(base) + i128::from(shift);
            let shifted = match model.cds_pair_wide(moved(start), moved(end)) {
                (Some(shifted_start), Some(shifted_end)) => Some((shifted_start, shifted_end)),
                _ => None,
            };
            (own, shifted, end)
        })
        .collect()
}

/// The bases a plugin reads under --hgvs for the allele as written, `allele`, on a row
/// the notation moved by `shift`: `shift_feature_seqs`
/// (TranscriptVariationAllele.pm:1330-1345) turns the allele left one base as many times as
/// the shift's length on a forward-strand transcript and as the allele's length less it on a
/// reverse-strand one, and `feature_seq` reads the result on the transcript's strand.
pub fn hgvs_rotated_allele(allele: &str, shift: i64, strand: i64) -> String {
    let length = i128::from(count(pyfmt::len(allele)));
    let moved = i128::from(shift).abs();
    let turns = if strand == 1 { moved } else { length - moved };
    let turned = if !allele.is_empty() && turns > 0 {
        let at = i64::try_from(turns % length).ok();
        Cow::Owned(format!(
            "{}{}",
            pyfmt::slice(allele, at, None),
            pyfmt::slice(allele, None, at)
        ))
    } else {
        Cow::Borrowed(allele)
    };
    if strand == -1 {
        reverse_complement(&turned)
    } else {
        turned.into_owned()
    }
}

/// Whether the two CDS pairs of a row are those of a later allele under --hgvs: `e` is
/// the reference's row of a second or later alternate allele of its record
/// (`later_alternate`), the candidate entry `a` prints the record's own pair, and the
/// reference prints the record's span moved by its HGVS_OFFSET `shift` or by the offset
/// of the last earlier allele's row that carries one (`earlier_allele_shift`): the pair an
/// earlier allele's protein notation left on the TranscriptVariation the record's alleles
/// share (OutputFactory.pm:1746 over BaseTranscriptVariation.pm:252-267), which is the
/// row's own offset only where the two alleles move alike.
pub fn later_allele_prints_shifted_pair(
    e: &Entry,
    a: &Entry,
    model: Option<&TranscriptModel>,
    shift: i64,
    record_rows: Option<&[&Entry]>,
) -> bool {
    let Some(model) = model else {
        return false;
    };
    if !later_alternate(e, record_rows) {
        return false;
    }
    distinct_moves(shift, earlier_allele_shift(e, record_rows))
        .into_iter()
        .any(|moved_by| printed_pair_is_moved(e, a, model, moved_by, None))
}

/// Whether the reference entry `e` prints the record's span moved by `move` bases
/// (`shifted_pair_readings`), the span that of a CDS pair the candidate entry `a` prints,
/// or of its pair `own` alone.
pub fn printed_pair_is_moved(
    e: &Entry,
    a: &Entry,
    model: &TranscriptModel,
    r#move: i64,
    own: Option<CdsPair>,
) -> bool {
    let printed = cds_pair_readings(e);
    let own_pairs = match own {
        None => cds_pair_readings(a),
        Some(own) => vec![own],
    };
    shifted_pair_readings(e, model, r#move, &own_pairs)
        .into_iter()
        .any(|(_, moved, _)| {
            moved.is_some_and(|(start, end)| printed.contains(&(Some(start), Some(end))))
        })
}

/// The one class that accepts the Downstream field `name` (a DownstreamProtein or a
/// ProteinLengthChange) whose values are `expected_value` on the reference entry `e` and
/// `actual_value` on the candidate entry `a`, with the class that must be declared beside
/// it (`None` for a row one cause changes), or `None`. Both are values, the entries
/// share one consequence set, or two that differ as the registry's
/// `start_cooccurrence_swap` admits, as the row rule of
/// `start_codon_kept_without_five_prime_utr` takes them
/// (`start_codon_kept_without_five_prime_utr_pair`) or as a class of a coding end takes
/// them on the transcript's record in the cache (`coding_end_class`), with
/// `frameshift_variant` and
/// no splice term and one CDS span without a gap, and the cache holds the transcript; a transcript under a
/// code no table is held for ends the run by name (`CodonTable::held`). Two CDS spans are one on the row of a later allele of `record_rows`
/// whose reference entry prints the record's span moved by its HGVS_OFFSET or by an earlier
/// allele's (`later_allele_prints_shifted_pair`): the own pair is then the candidate's. The
/// candidate's value is its rule on the own pair (`DownstreamTranscript.downstream_value`)
/// and the reference's is Downstream.pm:104-149 on the pair the plugin read: the pair and
/// the allele the notation shifted where the reference's non-zero HGVS_OFFSET reaches the
/// plugin (`shifted_pair_readings`, `hgvs_rotated_allele`), whatever earlier pair the
/// row prints; on a later allele's row without an offset, the printed pair an earlier
/// allele's notation left (`earlier_allele_shift`), with the allele as written; else the
/// own pair. Where
/// Downstream.pm gives another value on the own pair the class is
/// `plugin_hgvs_shifted_cds_coordinates`: alone where the candidate's value is that
/// own-pair value, the value Ensembl VEP writes without the flag, and beside the own pair's
/// class (`downstream_start_rule_class`) where it is not, since the start rule changes the
/// value as well. Where Downstream.pm gives one value on both pairs, a DownstreamProtein on
/// a row with one ProteinLengthChange is the own pair's class.
///
/// `record_rows` is `None` where the rows of the entry's record are not given.
pub fn downstream_class(
    e: &Entry,
    a: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
    sequences: &DownstreamSequences,
    record_rows: Option<&[&Entry]>,
) -> Result<Option<(String, Option<String>)>> {
    downstream_class_at(
        e,
        a,
        name,
        expected_value,
        actual_value,
        sequences,
        record_rows,
        &AtACodingEnd::NO_CACHE,
    )
}

/// `downstream_class` in a comparison that may have read the coding ends of the row's
/// transcript (`at`): two consequence sets that differ are also a row of the classes where
/// a class of a coding end takes the pair on that reading (`coding_end_class`).
#[allow(clippy::too_many_arguments)]
pub fn downstream_class_at(
    e: &Entry,
    a: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
    sequences: &DownstreamSequences,
    record_rows: Option<&[&Entry]>,
    at: &AtACodingEnd<'_>,
) -> Result<Option<(String, Option<String>)>> {
    if [expected_value, actual_value]
        .iter()
        .any(|value| matches!(*value, "-" | ABSENT))
    {
        return Ok(None);
    }
    let length_change = field_named(name, PROTEIN_LENGTH_CHANGE_FIELDS);
    if length_change
        && !(WHOLE_INTEGER_RE.is_match(expected_value) && WHOLE_INTEGER_RE.is_match(actual_value))
    {
        return Ok(None);
    }
    let terms = consequence_set(&e.consequence);
    let candidate_terms = consequence_set(&a.consequence);
    // One consequence set, the pair the registry sets aside as `start_cooccurrence_swap`,
    // the pair of `start_codon_kept_without_five_prime_utr`, or a pair a class of a coding
    // end takes on the transcript.
    if terms != candidate_terms
        && registry_class(
            &terms,
            &candidate_terms,
            &e.allele,
            &[START_COOCCURRENCE_SWAP],
        )
        .is_none()
        && !start_codon_kept_without_five_prime_utr_pair(e, a)?
        && coding_end_class(e, a, at).is_none()
    {
        return Ok(None);
    }
    if !terms.split(',').any(|term| term == FRAMESHIFT_TERM) || terms.contains("splice") {
        return Ok(None);
    }
    let own_pairs = cds_pair_readings(a);
    let two_pairs = own_pairs != cds_pair_readings(e);
    let Some(&(Some(start), Some(end))) = own_pairs.first() else {
        return Ok(None);
    };
    let allele = if e.allele == "-" {
        String::new()
    } else {
        pyfmt::upper(&e.allele)
    };
    if e.allele != "-" && !spelt_with(&allele, "ACGTN") {
        return Ok(None);
    }
    let offset = field_of(Some(e), HGVS_OFFSET_FIELDS);
    let shift = if WHOLE_NON_ZERO_INTEGER_RE.is_match(&offset) {
        clamped_int(&offset)
    } else {
        0
    };
    // A later allele's row that prints another pair reads the pair at its own offset, or with
    // none the pair the last earlier allele's notation left.
    let earlier = if two_pairs {
        earlier_allele_shift(e, record_rows).unwrap_or(0)
    } else {
        0
    };
    let read_at = if shift != 0 { shift } else { earlier };
    // The cache is read for a row the notation shifted or whose span starts on a codon's
    // third base.
    if read_at == 0 && (two_pairs || length_change || start.rem_euclid(3) != 0 || start > end) {
        return Ok(None);
    }
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(None);
    };
    transcript
        .table
        .held(&e.feature)
        .map_err(Error::comparator)?;
    let model = transcript.model.as_ref();
    if read_at != 0 && model.is_none() {
        return Ok(None);
    }
    if two_pairs && !later_allele_prints_shifted_pair(e, a, model, shift, record_rows) {
        return Ok(None);
    }
    let strand = match model {
        Some(model) => model.strand,
        None if field_of(Some(e), STRAND_FIELDS) == "-1" => -1,
        None => 1,
    };
    let own_allele = if strand == -1 {
        reverse_complement(&allele)
    } else {
        allele.clone()
    };
    let readings = match model {
        Some(model) if read_at != 0 => shifted_pair_readings(e, model, read_at, &own_pairs),
        _ => vec![(own_pairs[0], None, 0)],
    };
    let length_change_of = |entry: &Entry| field_of(Some(entry), PROTEIN_LENGTH_CHANGE_FIELDS);
    let one_length_change =
        !blank(&length_change_of(e)) && length_change_of(e) == length_change_of(a);
    let moves = distinct_moves(shift, Some(earlier));
    for (own, shifted, _) in readings {
        let (Some(start), Some(end)) = own else {
            continue;
        };
        // A later allele's reading is the one whose pair, moved by the row's offset or by the
        // earlier allele's, the reference prints.
        if two_pairs
            && !moves.iter().any(|moved_by| {
                model.is_some_and(|model| printed_pair_is_moved(e, a, model, *moved_by, Some(own)))
            })
        {
            continue;
        }
        let candidate = transcript.downstream_value(start, end, &own_allele, length_change, false);
        if candidate.as_deref() != Some(actual_value) {
            continue;
        }
        let unshifted = transcript.downstream_value(start, end, &own_allele, length_change, true);
        let reference = match shifted {
            Some((shifted_start, shifted_end)) => transcript.downstream_value(
                shifted_start,
                shifted_end,
                &hgvs_rotated_allele(&allele, shift, strand),
                length_change,
                true,
            ),
            None => unshifted.clone(),
        };
        if reference.as_deref() != Some(expected_value) {
            continue;
        }
        let start_rule = if length_change {
            None
        } else {
            downstream_start_rule_class(&transcript, start, end, &own_allele)
        };
        if reference != unshifted {
            if unshifted.as_deref() == Some(actual_value) {
                return Ok(Some((
                    PLUGIN_HGVS_SHIFTED_CDS_COORDINATES.to_string(),
                    None,
                )));
            }
            if start_rule.is_some() {
                return Ok(Some((
                    PLUGIN_HGVS_SHIFTED_CDS_COORDINATES.to_string(),
                    start_rule,
                )));
            }
        } else if one_length_change {
            if let Some(start_rule) = start_rule {
                return Ok(Some((start_rule, None)));
            }
        }
    }
    Ok(None)
}

/// The class of `DOWNSTREAM_CLASSES` whose shape the own pair `start` to `end` has
/// with the bases `allele` in its place, or `None`: a span that starts on a codon's third
/// base, `downstream_start_past_own_codon` where the offset of Downstream.pm:113 lies
/// inside the changed sequence or at its end, `downstream_utr_frame_past_coding_end` where
/// it lies past that end and the transcript has a 3' UTR.
pub fn downstream_start_rule_class(
    transcript: &DownstreamTranscript,
    start: i64,
    end: i64,
    allele: &str,
) -> Option<String> {
    let coding_length = i128::from(count(pyfmt::len(&transcript.coding)));
    let (start, end) = (i128::from(start), i128::from(end));
    if start.rem_euclid(3) != 0 || start > end || end > coding_length {
        return None;
    }
    let changed_length = start - 1 + i128::from(count(pyfmt::len(allele))) + coding_length - end;
    if start <= changed_length {
        return Some(DOWNSTREAM_START_PAST_OWN_CODON.to_string());
    }
    (!transcript.utr.is_empty()).then(|| DOWNSTREAM_UTR_FRAME_PAST_CODING_END.to_string())
}

/// Whether `value` is how the entry's format spells a field the row does not carry: no
/// Extra pair of a default-format row and no key of a JSON entry (`ABSENT`; `-` where
/// the Extra text is compared pair by pair), the `-` of a tab cell, and the empty cell of a
/// CSQ entry (no Location and no `variant_allele`), which holds every field of its header.
/// An empty value anywhere else is a value.
pub fn field_absent(e: &Entry, value: &str) -> bool {
    matches!(value, "-" | ABSENT)
        || (value.is_empty() && e.location.is_none() && !e.fields.contains_key("variant_allele"))
}

/// Whether the NMD `expected_value` of the reference entry `e` and `actual_value` of
/// the candidate entry `a` have the shape of `plugin_hgvs_shifted_cds_coordinates`: one
/// is `NMD_escaping_variant` and the other no value (`field_absent`), the entries share one
/// consequence set,
/// or two that differ as `stop_retained_without_stop_codon_test_pair` admits, each
/// with a term of `NMD_GATE_TERMS`, and one CDS span, or the two of a later allele of
/// `record_rows` (`later_allele_prints_shifted_pair`), whose own pair is the
/// candidate's and whose printed pair is the span moved by the row's offset or an earlier
/// allele's, the reference's non-zero HGVS_OFFSET
/// reaches the plugin (`shifted_pair_readings`), the reference carries the value exactly
/// when the rule of NMD.pm:98-109 holds on the shifted pair and the candidate exactly when it
/// holds on the own pair: both ends of the pair are coordinates and the CDS end is at most
/// 101 or the record's genomic end passes `variant_exon_check`
/// (`TranscriptModel.ends_where_nmd_escapes`).
///
/// `record_rows` is `None` where the rows of the entry's record are not given.
pub fn nmd_reads_shifted_pair(
    e: &Entry,
    a: &Entry,
    expected_value: &str,
    actual_value: &str,
    sequences: &DownstreamSequences,
    record_rows: Option<&[&Entry]>,
) -> Result<bool> {
    let reference = expected_value == NMD_ESCAPING;
    let candidate = actual_value == NMD_ESCAPING;
    let other_absent = if reference {
        field_absent(a, actual_value)
    } else {
        field_absent(e, expected_value)
    };
    if reference == candidate || !other_absent {
        return Ok(false);
    }
    let terms = consequence_set(&e.consequence);
    let candidate_terms = consequence_set(&a.consequence);
    // One consequence set, or the pair `stop_retained_without_stop_codon_test_pair` admits.
    if terms != candidate_terms && !stop_retained_without_stop_codon_test_pair(e, a)? {
        return Ok(false);
    }
    let gated = |side: &str| side.split(',').any(|term| NMD_GATE_TERMS.contains(&term));
    if !(gated(&terms) && gated(&candidate_terms)) {
        return Ok(false);
    }
    let offset = field_of(Some(e), HGVS_OFFSET_FIELDS);
    let own_pairs = cds_pair_readings(a);
    let printed = cds_pair_readings(e);
    if !WHOLE_NON_ZERO_INTEGER_RE.is_match(&offset) || own_pairs.is_empty() {
        return Ok(false);
    }
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(false);
    };
    let Some(model) = transcript.model.as_ref() else {
        return Ok(false);
    };
    let two_pairs = own_pairs != printed;
    let shift = clamped_int(&offset);
    if two_pairs && !later_allele_prints_shifted_pair(e, a, Some(model), shift, record_rows) {
        return Ok(false);
    }
    let moves = distinct_moves(shift, earlier_allele_shift(e, record_rows));
    for (own, shifted, genomic_end) in shifted_pair_readings(e, model, shift, &own_pairs) {
        let Some((_, shifted_end)) = shifted else {
            continue;
        };
        if two_pairs
            && !moves
                .iter()
                .any(|moved_by| printed_pair_is_moved(e, a, model, *moved_by, Some(own)))
        {
            continue;
        }
        let escapes = model.ends_where_nmd_escapes(genomic_end);
        let own_escapes = match own {
            (Some(_), Some(own_end)) => own_end <= 101 || escapes,
            _ => false,
        };
        if (shifted_end <= 101 || escapes) == reference && own_escapes == candidate {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The first and last codon an entry's protein position names, the lower first: the
/// `Protein_position` read up to a `/total` suffix, or the JSON `protein_start` and
/// `protein_end`; `None` where either is no number.
pub fn protein_span(e: &Entry) -> Option<(i64, i64)> {
    let [first, last] = bound_texts(e, "Protein_position", ["protein_start", "protein_end"]);
    if !(pyfmt::is_digit(first) && pyfmt::is_digit(last)) {
        return None;
    }
    let first = int_of_text(first);
    let last = int_of_text(last);
    let (first, last) = (first?, last?);
    Some((first.min(last), first.max(last)))
}

/// The first and last codon of a CDS pair, the lower first; an insertion's pair names the
/// codons of its two flanks.
pub fn codon_span(pair: (i64, i64)) -> (i64, i64) {
    let codon = |base: i64| {
        let codon = (i128::from(base) + 2).div_euclid(3);
        i64::try_from(codon).unwrap_or(i64::MAX)
    };
    (codon(pair.0.min(pair.1)), codon(pair.0.max(pair.1)))
}

/// The coding sequence's translation over its codons `codons` (first, last) under the
/// transcript's code as `reading` reads it, a stop `*`; `None` for a codon past the
/// sequence's end or one the table does not read.
pub fn codon_residues(
    transcript: &DownstreamTranscript,
    codons: (i64, i64),
    reading: Reading,
) -> Option<String> {
    let code = match &transcript.table {
        CodonTable::Held(table) => genetic_code(*table, reading),
        CodonTable::Lacking(spelt) => crate::uncaught!("{spelt}"),
    };
    let Some(code) = code else {
        crate::uncaught!("{:?}", transcript.table)
    };
    let mut residues = String::new();
    for codon in codons.0..=codons.1 {
        let end = codon.saturating_mul(3);
        let bases = pyfmt::slice(&transcript.coding, Some(end.saturating_sub(3)), Some(end));
        residues.push(char::from(code.amino_acid(bases.as_bytes())?));
    }
    Some(residues)
}

/// `hgvs_shift_shared_across_alleles` on the differing field `name` of a later allele's
/// row: the reference entry `e` is a second or later alternate's whose first allele's row
/// carries a non-zero HGVS_OFFSET (`first_allele_shift`), the candidate entry `a`
/// prints the record's own CDS pair and the reference that span moved by the first allele's
/// offset (`shifted_pair_readings`). A CDS position is the class's there. A protein
/// position is where each entry's is the codons of its pair (`codon_span`) and the two
/// differ. An `Amino_acids` is on such a row where the two values share their alternate
/// part and each reference part is the translation of its entry's codons
/// (`codon_residues`).
pub fn later_allele_shared_shift(
    e: &Entry,
    a: &Entry,
    name: &str,
    record_rows: Option<&[&Entry]>,
    sequences: &DownstreamSequences,
) -> Result<bool> {
    let shift = first_allele_shift(e, record_rows);
    let own_pairs = cds_pair_readings(a);
    let printed = cds_pair_readings(e);
    let Some(shift) = shift else {
        return Ok(false);
    };
    if !matches!(own_pairs.first(), Some((Some(_), Some(_)))) || own_pairs == printed {
        return Ok(false);
    }
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(false);
    };
    let Some(model) = transcript.model.as_ref() else {
        return Ok(false);
    };
    for (own, moved, _) in shifted_pair_readings(e, model, shift, &own_pairs) {
        let (Some(moved), (Some(own_start), Some(own_end))) = (moved, own) else {
            continue;
        };
        let own = (own_start, own_end);
        if moved == own || !printed.contains(&(Some(moved.0), Some(moved.1))) {
            continue;
        }
        if field_named(name, CDS_POSITION_FIELDS) {
            return Ok(true);
        }
        let (own_codons, moved_codons) = (codon_span(own), codon_span(moved));
        if own_codons == moved_codons
            || protein_span(a) != Some(own_codons)
            || protein_span(e) != Some(moved_codons)
        {
            continue;
        }
        if field_named(name, PROTEIN_POSITION_FIELDS) {
            return Ok(true);
        }
        transcript
            .table
            .held(&e.feature)
            .map_err(Error::comparator)?;
        let amino_acids = field_of(Some(e), AMINO_ACIDS_FIELDS);
        let (reference, _, changed) = pyfmt::partition(&amino_acids, "/");
        let own_amino_acids = field_of(Some(a), AMINO_ACIDS_FIELDS);
        let (own_reference, _, own_changed) = pyfmt::partition(&own_amino_acids, "/");
        if changed.is_empty() || changed != own_changed {
            continue;
        }
        let read = |codons, reading| codon_residues(&transcript, codons, reading);
        if read(moved_codons, Reading::EnsemblVep).as_deref() == Some(reference)
            && read(own_codons, Reading::VepRs).as_deref() == Some(own_reference)
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// `blosum62_hgvs_shift_shared_across_alleles` on a differing BLOSUM62: the reference
/// entry `e` writes a value and the candidate entry `a` none, on a later allele's row
/// whose CDS position and protein position `hgvs_shift_shared_across_alleles` takes
/// (`later_allele_shared_shift`: the candidate prints the record's own pair and the
/// reference that span moved by the first allele's offset, each protein position the codons
/// of its pair). The reference's `Amino_acids` is one residue over one (Blosum62.pm:110),
/// the first the translation of the moved pair's one codon (`codon_residues`) and the
/// second the alternate part of the candidate's own `Amino_acids`, which fails that test;
/// and the value is the matrix cell of the reference's pair (`BLOSUM62_MATRIX`).
pub fn blosum62_later_allele_shared_shift(
    e: &Entry,
    a: &Entry,
    ev: &str,
    av: &str,
    record_rows: Option<&[&Entry]>,
    sequences: &DownstreamSequences,
) -> Result<bool> {
    if blank(ev) || !blank(av) || !blank(&field_of(Some(a), BLOSUM62_FIELDS)) {
        return Ok(false);
    }
    for fields in [CDS_POSITION_FIELDS, PROTEIN_POSITION_FIELDS] {
        if !later_allele_shared_shift(e, a, fields[0], record_rows, sequences)? {
            return Ok(false);
        }
    }
    let amino_acids = field_of(Some(e), AMINO_ACIDS_FIELDS);
    let scored = WHOLE_PEPTIDE_PAIR_RE.captures(&amino_acids);
    let own = field_of(Some(a), AMINO_ACIDS_FIELDS);
    let Some(scored) = scored else {
        return Ok(false);
    };
    let residue = |group: usize| scored[group].chars().next().unwrap_or_default();
    let (reference, changed) = (residue(1), residue(2));
    if WHOLE_PEPTIDE_PAIR_RE.is_match(&own) || pyfmt::partition(&own, "/").2 != &scored[2] {
        return Ok(false);
    }
    let transcript = sequences.transcript(e)?;
    let shift = first_allele_shift(e, record_rows);
    let (Some(transcript), Some(shift)) = (transcript, shift) else {
        return Ok(false);
    };
    let Some(model) = transcript.model.as_ref() else {
        return Ok(false);
    };
    transcript
        .table
        .held(&e.feature)
        .map_err(Error::comparator)?;
    let printed = cds_pair_readings(e);
    let moved_codons: BTreeSet<(i64, i64)> =
        shifted_pair_readings(e, model, shift, &cds_pair_readings(a))
            .into_iter()
            .filter_map(|(_, moved, _)| moved)
            .filter(|moved| printed.contains(&(Some(moved.0), Some(moved.1))))
            .map(codon_span)
            .collect();
    if !moved_codons.iter().any(|codons| {
        codon_residues(&transcript, *codons, Reading::EnsemblVep)
            .is_some_and(|residues| residues == scored[1])
    }) {
        return Ok(false);
    }
    Ok(BLOSUM62_MATRIX.get(&(reference, changed)).copied() == Some(ev))
}

/// Whether the entry's `SHIFT_LENGTH` (the JSON `shift_length`) is a positive integer:
/// a row `--shift_3prime 1` moved, printed under `--shift_length`. Digits `int` does
/// not read are a malformed input (`integer_of`).
pub fn positive_shift(e: &Entry) -> Result<bool> {
    let value = SHIFT_LENGTH_FIELDS
        .iter()
        .find_map(|field| e.fields.get(*field))
        .map_or(ABSENT, String::as_str);
    Ok(pyfmt::is_digit(value) && integer_text_of(value, "SHIFT_LENGTH")? != "0")
}

/// Whether the entry's `SHIFT_LENGTH`, where the run prints one, is positive: a row
/// `--shift_3prime 1` moved; a run without the field says nothing about the shift.
pub fn shift_field_consistent(e: &Entry) -> Result<bool> {
    if SHIFT_LENGTH_FIELDS
        .iter()
        .all(|field| !e.fields.contains_key(*field))
    {
        return Ok(true);
    }
    positive_shift(e)
}

/// The value of the entry's first field named by one of `patterns` (`ABSENT` for none),
/// as `FIELD_CLASSES` name fields: a full name or a prefix, the case set aside.
pub fn field_of(e: Option<&Entry>, patterns: &[&str]) -> String {
    e.and_then(|e| {
        e.fields
            .iter()
            .find(|(name, _)| field_named(name, patterns))
    })
    .map_or_else(|| ABSENT.to_string(), |(_, value)| value.clone())
}
