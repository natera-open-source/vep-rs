// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein field classes: at a stop codon, Ensembl VEP's description read from the codon
//! window against the one the two proteins give, and a stop named alone on an incomplete
//! coding sequence; a frameshift that names one residue on both sides of its number; and,
//! at the first residue, the count of a frameshift on a transcript on which no start codon
//! stands.

use super::cache::iter_cache_transcripts;
use super::downstream::{downstream_sequences, translate_codons, CodingSequence, CodingSequences};
use super::fields::{FieldKey, FieldModels, FieldQuestion, FieldRow, FieldRun, Verdict};
use super::keys::record_reference;
use super::records::{first_five_columns, RunRecord};
use crate::error::Failure;
use crate::json::Value;
use crate::pytext::{is_dir, path_join, py_str};
use crate::rules::cache_transcripts::{location_span, CacheTranscripts as CodingEndTranscripts};
use crate::rules::coding_ends::terms::{
    count_alone_differs_at_the_first_residue, frameshift_counted_at_the_first_residue, Cells,
    RowFacts, WhatWasRead,
};
use crate::rules::row::transcript_feature_type;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::LazyLock;

/// The class of an `HGVSp` at the stop codon of one protein, Ensembl VEP's description read
/// from the codon window and vep-rs's the one the two proteins give.
pub const HGVSP_STOP_CODON_PROTEIN_COMPARISON: &str = "hgvsp_stop_codon_protein_comparison";

/// The class of an `HGVSp` in which Ensembl VEP names the stop alone on a row that shows an
/// incomplete coding sequence, where vep-rs writes none.
pub const HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE: &str =
    "hgvsp_stop_named_on_incomplete_coding_sequence";

/// The class of an `HGVSp` in which Ensembl VEP writes a frameshift with the unknown count
/// that names one residue on both sides of its number, where vep-rs writes none.
pub const HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE: &str =
    "hgvsp_frameshift_with_no_changed_residue";

/// The class of an `HGVSp` in which both engines name a frameshift at residue 1 of a
/// transcript on which no start codon stands, Ensembl VEP with the unknown count
/// (TranscriptVariationAllele.pm:2483) and vep-rs with the count of residues to the first
/// stop codon of the new reading frame.
pub const HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE: &str =
    "hgvsp_frameshift_count_at_the_first_residue";

/// A description that names the stop alone: unchanged, deleted, replaced by one residue or
/// extended.
pub const STOP_NAMED_ALONE: &str = r"Ter([1-9]\d*)(=|del|[A-Z][a-z]{2}(?:extTer(?:[1-9]\d*|\?))?)";

/// The ten descriptions the stop codon class reads: the stop unchanged, an extension named
/// at it, the stop alone replaced, residues inserted between the last residue and the stop,
/// a range that ends at the stop replaced, a frameshift, the stop deleted before an
/// extension's count, a deletion before an extension's count, one residue and a stop
/// inserted between two residues, and two or more residues and a stop inserted between two
/// residues.
pub const STOP_UNCHANGED: &str = r"Ter([1-9]\d*)=";
pub const STOP_EXTENSION: &str = r"Ter([1-9]\d*)([A-Z][a-z]{2})extTer([1-9]\d*|\?)";
pub const STOP_REPLACED: &str = r"Ter([1-9]\d*)delins((?:[A-Z][a-z]{2})+)";
pub const STOP_INSERTION: &str = r"[A-Z][a-z]{2}([1-9]\d*)_Ter([1-9]\d*)ins((?:[A-Z][a-z]{2})+)";
pub const STOP_RANGE_REPLACED: &str =
    r"([A-Z][a-z]{2})([1-9]\d*)_Ter([1-9]\d*)delins((?:[A-Z][a-z]{2})+)";
pub const STOP_FRAMESHIFT: &str = r"([A-Z][a-z]{2})([1-9]\d*)([A-Z][a-z]{2})fsTer([1-9]\d*|\?)";
pub const STOP_DELETED_EXTENSION: &str = r"Ter([1-9]\d*)delextTer([1-9]\d*|\?)";
/// A deletion of one residue or of a range of residues with an extension's count after it:
/// the deletion as written, its first residue's number, its last residue's number where it
/// names a range, and the count.
pub const DELETION_EXTENSION: &str =
    r"([A-Z][a-z]{2}([1-9]\d*)(?:_[A-Z][a-z]{2}([1-9]\d*))?del)extTer([1-9]\d*|\?)";
pub const RESIDUE_AND_STOP_INSERTED: &str =
    r"([A-Z][a-z]{2})([1-9]\d*)_([A-Z][a-z]{2})([1-9]\d*)ins([A-Z][a-z]{2})Ter";
pub const RESIDUES_AND_STOP_INSERTED: &str =
    r"([A-Z][a-z]{2})([1-9]\d*)_([A-Z][a-z]{2})([1-9]\d*)ins((?:[A-Z][a-z]{2}){2,})Ter";

/// A pattern of a notation as one that takes a whole text, its digits the ASCII ones.
pub fn whole(pattern: &str) -> Regex {
    Regex::new(&format!("^(?-u:{pattern})$")).expect("a pattern of a notation")
}

static NAMED_ALONE: LazyLock<Regex> = LazyLock::new(|| whole(STOP_NAMED_ALONE));
static UNCHANGED: LazyLock<Regex> = LazyLock::new(|| whole(STOP_UNCHANGED));
static EXTENSION: LazyLock<Regex> = LazyLock::new(|| whole(STOP_EXTENSION));
static REPLACED: LazyLock<Regex> = LazyLock::new(|| whole(STOP_REPLACED));
static INSERTION: LazyLock<Regex> = LazyLock::new(|| whole(STOP_INSERTION));
static RANGE_REPLACED: LazyLock<Regex> = LazyLock::new(|| whole(STOP_RANGE_REPLACED));
static FRAMESHIFT: LazyLock<Regex> = LazyLock::new(|| whole(STOP_FRAMESHIFT));
static DELETED_EXTENSION: LazyLock<Regex> = LazyLock::new(|| whole(STOP_DELETED_EXTENSION));
static DELETED_WITH_EXTENSION: LazyLock<Regex> = LazyLock::new(|| whole(DELETION_EXTENSION));
static AND_STOP_INSERTED: LazyLock<Regex> = LazyLock::new(|| whole(RESIDUE_AND_STOP_INSERTED));
static RESIDUES_INSERTED: LazyLock<Regex> = LazyLock::new(|| whole(RESIDUES_AND_STOP_INSERTED));

/// The transcripts of the corpus cache whose coding sequence ends in an incomplete codon,
/// the coding sequence and genetic code of each of its coding transcripts, the coding ends
/// of its transcripts, and the records of the run.
#[derive(Default)]
pub struct Models {
    pub partial_codon_models: Option<HashSet<String>>,
    /// `downstream_sequences` of the corpus cache, from which the ninth shape reads the
    /// reference protein where more than one residue stands before the inserted stop.
    pub coding_sequences: Option<CodingSequences>,
    /// `peptides_that_hold_a_stop` of the corpus cache, beside which the eleventh shape
    /// reads a deletion that takes a stop.
    pub held_stop_peptides: Option<HashMap<String, String>>,
    /// The coding ends of the corpus cache's transcripts as the tuple comparator reads a
    /// record (`crate::rules::cache_transcripts`), the transcripts of a sequence read when
    /// a row on it first asks; `None` for a cache without a `transcripts` directory. That
    /// reader hands the reading no 3' UTR, so the reading keeps the sequence of a transcript
    /// that shows no start codon only where the transcript ends at its coding end.
    pub coding_ends: Option<CodingEndTranscripts>,
    /// The run's records, each cut to its first five columns: their alleles give the count
    /// of reference bases a row's Allele stands in for (`record_reference`).
    pub records: Vec<RunRecord>,
}

impl Models {
    /// The models of a run whose records are `records`.
    pub fn load(run: &FieldRun<'_>, records: &[RunRecord]) -> Result<Models, Failure> {
        let threads = crate::rows::blocks::thread_count()?;
        let coding_ends = is_dir(&path_join(run.cache_dir, "transcripts"))
            .then(|| CodingEndTranscripts::new(Path::new(run.cache_dir), threads));
        Ok(Models {
            partial_codon_models: Some(partial_terminal_codon_models(run.cache_dir)?),
            coding_sequences: Some(downstream_sequences(run.cache_dir)?),
            held_stop_peptides: Some(peptides_that_hold_a_stop(run.cache_dir)?),
            coding_ends,
            records: records
                .iter()
                .map(|record| RunRecord {
                    vcf: first_five_columns(&record.vcf).to_string(),
                    manifest_index: record.manifest_index,
                    reference_warning: record.reference_warning.clone(),
                })
                .collect(),
        })
    }
}

/// The decimal number one above `digits`, a number in ASCII digits of any length.
pub fn successor(digits: &str) -> String {
    let mut above: Vec<u8> = digits.bytes().collect();
    for digit in above.iter_mut().rev() {
        if *digit == b'9' {
            *digit = b'0';
        } else {
            *digit += 1;
            return String::from_utf8(above).expect("ASCII digits");
        }
    }
    format!("1{}", String::from_utf8(above).expect("ASCII digits"))
}

/// Whether the number `low` is below the number `high`, both in ASCII digits without a
/// leading zero.
fn below(low: &str, high: &str) -> bool {
    (low.len(), low) < (high.len(), high)
}

/// The sum of two numbers in ASCII digits of any length.
fn sum(left: &str, right: &str) -> String {
    let (mut left, mut right) = (left.bytes().rev(), right.bytes().rev());
    let mut digits = Vec::new();
    let mut carry = 0;
    loop {
        let (a, b) = (left.next(), right.next());
        if a.is_none() && b.is_none() && carry == 0 {
            break;
        }
        let total = a.map_or(0, |a| a - b'0') + b.map_or(0, |b| b - b'0') + carry;
        digits.push(b'0' + total % 10);
        carry = total / 10;
    }
    digits.reverse();
    String::from_utf8(digits).expect("ASCII digits")
}

/// The residues of a run of three-letter names.
fn residues(names: &str) -> Vec<&str> {
    (0..names.len())
        .step_by(3)
        .map(|at| &names[at..at + 3])
        .collect()
}

/// The `FLAGS` of the reference row, empty for a row that is absent or carries none.
fn row_flags(perl_row: Option<&FieldRow>) -> &str {
    perl_row
        .and_then(|row| row.get("FLAGS"))
        .map_or("", String::as_str)
}

/// Whether a consequence set holds `term`.
fn holds_term(consequences: &str, term: &str) -> bool {
    consequences.split(',').any(|held| held == term)
}

/// The name and the description of an `HGVSp` value, cut at its first `:p.` and an escaped
/// `=` set aside; `None` for a value without the separator.
fn described(value: &str) -> Option<(&str, String)> {
    let (name, description) = value.split_once(":p.")?;
    Some((name, description.replace("%3D", "=")))
}

/// The positions a column names, up to a `/total` suffix: one position, or the two of a
/// range. `None` where a bound is no number in ASCII digits.
fn column_bounds(column: &str) -> Option<(u64, u64)> {
    let named = column.split('/').next().unwrap_or_default();
    let (first, last) = named.split_once('-').unwrap_or((named, named));
    let number = |text: &str| {
        (text.bytes().all(|digit| digit.is_ascii_digit()) && !text.starts_with('0'))
            .then(|| text.parse::<u64>().ok())
            .flatten()
    };
    Some((number(first)?, number(last)?))
}

/// The tenth shape of `hgvsp_stop_codon_protein_comparison`, on an `HGVSp` that Ensembl VEP
/// leaves empty: vep-rs writes the stop as unchanged at the residue after the one codon the
/// reference row's `Protein_position` names; the row's terms hold `stop_retained_variant`
/// and its `FLAGS` name no `cds_end_NF`; its `CDS_position` lies inside that codon; and its
/// `Amino_acids` are that codon's residue against the same residue and then a stop. By the
/// row the alternate protein is the reference protein up to that residue and then a stop,
/// at the residue where vep-rs's value has the reference protein's stop.
fn stop_unchanged_where_ensembl_vep_writes_none(
    key: &FieldKey,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
) -> bool {
    if row_flags(perl_row).contains("cds_end_NF")
        || !holds_term(&key.consequences, "stop_retained_variant")
    {
        return false;
    }
    let Some((_, candidate)) = described(vep_rs) else {
        return false;
    };
    let Some(stays) = UNCHANGED.captures(&candidate) else {
        return false;
    };
    let column = |name: &str| {
        perl_row
            .and_then(|row| row.get(name))
            .map_or("-", String::as_str)
    };
    let (Some((codon, last_codon)), Some((first_base, last_base))) = (
        column_bounds(column("Protein_position")),
        column_bounds(column("CDS_position")),
    ) else {
        return false;
    };
    let codon_end = 3 * u128::from(codon);
    let inside_the_codon = codon_end <= u128::from(first_base.min(last_base)) + 2
        && u128::from(first_base.max(last_base)) <= codon_end;
    let kept_and_stopped =
        column("Amino_acids")
            .split_once('/')
            .is_some_and(|(residue, alternate)| {
                residue.len() == 1
                    && residue
                        .bytes()
                        .all(|aa| aa.is_ascii_uppercase() && aa != b'X')
                    && alternate
                        .strip_prefix(residue)
                        .is_some_and(|after| after.starts_with('*'))
            });
    codon == last_codon
        && inside_the_codon
        && kept_and_stopped
        && stays[1] == successor(&codon.to_string())
}

/// `hgvsp_stop_codon_protein_comparison` on an `HGVSp` the engines write differently: both
/// values name one protein, the reference row's `FLAGS` (`perl_row`) name no `cds_end_NF`,
/// and with an escaped `=` set aside Ensembl VEP's description is the stop replaced by a
/// stop against vep-rs's stop unchanged at the same residue above 1; one residue and a stop
/// inserted between two neighbouring residues, the residue the second neighbour's, on a row
/// whose terms hold `stop_retained_variant`, against vep-rs's stop unchanged at the residue
/// after that neighbour; a range from a residue to the stop replaced by residues, none of
/// them a stop, on a row whose terms hold
/// `stop_lost`, against vep-rs's frameshift named at that residue with the first of them,
/// its count unknown, or above their number and not the one that leaves the stop at the
/// stop's own residue; or, against vep-rs's extension named at the stop: the stop alone
/// replaced by the first residues that follow the protein (the extension's first residue the
/// first named, its count unknown, above the residues named, or their number where no `Ter`
/// closes them), the stop replaced by every added residue and then a stop (the count their
/// number), residues inserted between the last residue and the stop (the count their
/// number), the stop unchanged on a row whose terms hold `stop_lost`, the same extension
/// with an unknown count on such a row, or the stop written as deleted before the
/// extension's own count on such a row. Where Ensembl VEP writes no value, the pair is the
/// tenth shape's (`stop_unchanged_where_ensembl_vep_writes_none`).
pub fn stop_codon_protein_comparison_pair(
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
) -> bool {
    if matches!(vep, "" | "-") {
        return stop_unchanged_where_ensembl_vep_writes_none(key, vep_rs, perl_row);
    }
    if matches!(vep_rs, "" | "-") || row_flags(perl_row).contains("cds_end_NF") {
        return false;
    }
    let Some((vep_name, reference)) = described(vep) else {
        return false;
    };
    let (vep_rs_name, candidate) = described(vep_rs).unwrap_or((vep_rs, String::new()));
    if vep_name != vep_rs_name {
        return false;
    }
    let lost = holds_term(&key.consequences, "stop_lost");
    let replaced = REPLACED.captures(&reference);
    let stays = UNCHANGED.captures(&candidate);
    // The stop replaced by a stop.
    if let (Some(replaced), Some(stays)) = (&replaced, &stays) {
        if &replaced[2] == "Ter" && replaced[1] == stays[1] && &stays[1] != "1" {
            return true;
        }
    }

    // One residue and a stop inserted between two neighbouring residues, against the stop
    // unchanged.
    if let (Some(placed), Some(stays)) = (AND_STOP_INSERTED.captures(&reference), &stays) {
        let (left, before, right, after) = (&placed[1], &placed[2], &placed[3], &placed[4]);
        return left != "Ter"
            && right != "Ter"
            && &placed[5] == right
            && successor(before) == after
            && stays[1] == successor(after)
            && holds_term(&key.consequences, "stop_retained_variant");
    }

    // A range from a residue to the stop replaced, against a frameshift named at that residue.
    if let (Some(ranged), Some(frameshift)) = (
        RANGE_REPLACED.captures(&reference),
        FRAMESHIFT.captures(&candidate),
    ) {
        let (residue, at, stop) = (&ranged[1], &ranged[2], &ranged[3]);
        let new = residues(&ranged[4]);
        let shifted = &frameshift[4];
        let counted = shifted == "?"
            || (below(&new.len().to_string(), shifted) && sum(shifted, at) != successor(stop));
        return residue != "Ter"
            && !new.contains(&"Ter")
            && below(at, stop)
            && &frameshift[1] == residue
            && &frameshift[2] == at
            && frameshift[3] == *new[0]
            && counted
            && lost;
    }

    let Some(extension) = EXTENSION.captures(&candidate) else {
        return false;
    };
    let (number, first, count) = (&extension[1], &extension[2], &extension[3]);
    if number == "1" || first == "Ter" {
        return false;
    }
    // An extension with an unknown count against the same extension with a count.
    if let Some(unknown) = EXTENSION.captures(&reference) {
        let same = &unknown[1] == number && &unknown[2] == first;
        if same && &unknown[3] == "?" && count != "?" && lost {
            return true;
        }
    }
    // The stop written as deleted with an extension count, against an extension with that count.
    if let Some(deleted) = DELETED_EXTENSION.captures(&reference) {
        if &deleted[1] == number && &deleted[2] == count && lost {
            return true;
        }
    }
    let mut named = match &replaced {
        Some(replaced) => residues(&replaced[2]),
        None => Vec::new(),
    };
    let closed = named.last() == Some(&"Ter");
    if closed {
        named.pop();
    }
    let at_the_stop = replaced
        .as_ref()
        .is_some_and(|replaced| &replaced[1] == number);
    let opens = !named.is_empty() && !named.contains(&"Ter") && at_the_stop && named[0] == first;
    let residues_named = named.len().to_string();
    // The stop alone replaced by the first residues that follow the protein.
    if opens
        && (count == "?" || below(&residues_named, count) || (count == residues_named && !closed))
    {
        return true;
    }

    // The stop replaced by every added residue and then the stop the inserted bases bring.
    if opens && closed && count == residues_named {
        return true;
    }

    // Residues inserted between the last residue and the stop.
    if let Some(inserted) = INSERTION.captures(&reference) {
        let between = residues(&inserted[3]);
        if successor(&inserted[1]) == inserted[2]
            && &inserted[2] == number
            && !between.contains(&"Ter")
            && between[0] == first
            && count == between.len().to_string()
        {
            return true;
        }
    }

    // The stop written as unchanged where the row's terms say it is lost.
    UNCHANGED
        .captures(&reference)
        .is_some_and(|unchanged| &unchanged[1] == number && lost)
}

/// The cache's transcripts whose coding sequence is not a whole number of codons, by stable
/// id: a coding sequence that ends inside a codon, which no stop codon closes. Empty when
/// the cache is absent. A `variation_effect_feature_cache` that is no object and a
/// `translateable_seq` that is a number are refused.
pub fn partial_terminal_codon_models(cache_dir: &str) -> Result<HashSet<String>, Failure> {
    let mut partial = HashSet::new();
    for item in iter_cache_transcripts(cache_dir) {
        let (_chromosome, transcript) = item?;
        // The reader gives only the records that carry a stable id.
        let stable_id = transcript.get("stable_id").map(py_str).unwrap_or_default();
        let refused = |what: &str| {
            Failure::own(format!(
                "the transcript {stable_id} of the corpus cache holds {what}"
            ))
        };
        let coding = match transcript
            .get("variation_effect_feature_cache")
            .filter(|held| held.truthy())
        {
            None => None,
            Some(Value::Object(held)) => held.get("translateable_seq").filter(|seq| seq.truthy()),
            Some(_) => {
                return Err(refused(
                    "a variation_effect_feature_cache that is no object",
                ))
            }
        };
        let length = match coding {
            None => 0,
            Some(Value::Str(sequence)) => sequence.chars().count(),
            Some(Value::List(items)) => items.len(),
            Some(Value::Object(members)) => members.len(),
            Some(_) => return Err(refused("a translateable_seq that is a number")),
        };
        if length % 3 != 0 {
            partial.insert(stable_id);
        }
    }
    Ok(partial)
}

/// The three-letter name Ensembl VEP writes for one of the twenty residues of the standard
/// alphabet, by its letter; `None` for any other letter.
fn residue_name(letter: u8) -> Option<&'static str> {
    Some(match letter {
        b'A' => "Ala",
        b'R' => "Arg",
        b'N' => "Asn",
        b'D' => "Asp",
        b'C' => "Cys",
        b'Q' => "Gln",
        b'E' => "Glu",
        b'G' => "Gly",
        b'H' => "His",
        b'I' => "Ile",
        b'L' => "Leu",
        b'K' => "Lys",
        b'M' => "Met",
        b'F' => "Phe",
        b'P' => "Pro",
        b'S' => "Ser",
        b'T' => "Thr",
        b'W' => "Trp",
        b'Y' => "Tyr",
        b'V' => "Val",
        _ => return None,
    })
}

/// The reference protein a cache transcript's coding sequence spells, one letter a
/// residue: its codons translated under the transcript's own genetic code up to the stop
/// codon that closes the sequence. `None` under a code that has no table and where no stop
/// codon closes the sequence: it is no whole number of codons, or its first stop codon is
/// not its last codon.
pub fn protein_closed_by_a_stop(sequence: &CodingSequence) -> Option<String> {
    let translation = translate_codons(&sequence.coding, sequence.codon_table)?;
    let residues = translation.strip_suffix('*')?;
    let whole_codons = sequence.coding.chars().count() == 3 * translation.len();
    (whole_codons && !residues.contains('*')).then(|| residues.to_string())
}

/// The ninth shape of `hgvsp_stop_codon_protein_comparison` with more than one residue
/// before the inserted stop, which the row alone does not decide: both values name one
/// protein, the reference row's `FLAGS` (`perl_row`) name no `cds_end_NF`, and with an
/// escaped `=` and a pair of prediction parentheses that both carry set aside, Ensembl VEP's
/// description is two or more residues and a stop inserted between two neighbouring
/// residues against vep-rs's stop unchanged. The coding sequence of the row's transcript in
/// the corpus cache (`coding_sequences`) shows the protein kept: the left neighbour is the
/// residue of its number, the right neighbour and the inserted residues are the protein's
/// from the next number to its last, and the unchanged stop stands at the number after the
/// last. The row's terms are not read. `false` without the cache's sequences.
pub fn inserted_residues_keep_the_protein_pair(
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    coding_sequences: Option<&CodingSequences>,
) -> bool {
    if matches!(vep, "" | "-")
        || matches!(vep_rs, "" | "-")
        || row_flags(perl_row).contains("cds_end_NF")
    {
        return false;
    }
    let (Some((vep_name, reference)), Some((vep_rs_name, candidate))) =
        (described(vep), described(vep_rs))
    else {
        return false;
    };
    if vep_name != vep_rs_name {
        return false;
    }
    let inside = |description: &str| -> Option<String> {
        let inside = description.strip_prefix('(')?.strip_suffix(')')?;
        Some(inside.to_string())
    };
    let (reference, candidate) = match (inside(&reference), inside(&candidate)) {
        (Some(reference), Some(candidate)) => (reference, candidate),
        _ => (reference, candidate),
    };
    let (Some(placed), Some(stays)) = (
        RESIDUES_INSERTED.captures(&reference),
        UNCHANGED.captures(&candidate),
    ) else {
        return false;
    };
    let Some(protein) = coding_sequences
        .and_then(|sequences| sequences.get(&key.feature))
        .and_then(protein_closed_by_a_stop)
    else {
        return false;
    };
    let protein = protein.as_bytes();
    let named = residues(&placed[5]);
    // The count of residues before the inserted ones, which is the left neighbour's number.
    let Some(kept) = protein.len().checked_sub(named.len()) else {
        return false;
    };
    if kept == 0
        || placed[2] != kept.to_string()
        || successor(&placed[2]) != placed[4]
        || stays[1] != (protein.len() + 1).to_string()
    {
        return false;
    }
    let is_residue = |name: &str, at: usize| residue_name(protein[at]) == Some(name);
    is_residue(&placed[1], kept - 1)
        && is_residue(&placed[3], kept)
        && named
            .iter()
            .enumerate()
            .all(|(place, name)| is_residue(name, kept + place))
}

/// The peptide of each transcript of the corpus cache whose peptide holds a stop, by stable
/// id: a peptide that goes on past a stop and names the residues after it. Empty when the
/// cache is absent.
pub fn peptides_that_hold_a_stop(cache_dir: &str) -> Result<HashMap<String, String>, Failure> {
    let mut peptides = HashMap::new();
    for item in iter_cache_transcripts(cache_dir) {
        let (_chromosome, transcript) = item?;
        let peptide = match transcript.get("variation_effect_feature_cache") {
            Some(Value::Object(held)) => held.get("peptide").and_then(Value::as_str),
            _ => None,
        };
        if let Some(peptide) = peptide.filter(|peptide| peptide.contains('*')) {
            let stable_id = transcript.get("stable_id").map(py_str).unwrap_or_default();
            peptides.insert(stable_id, peptide.to_string());
        }
    }
    Ok(peptides)
}

/// The count Ensembl VEP writes after a deletion that is no frameshift and whose
/// description opens at residue number `first_named`, on a peptide of `peptide_length`
/// residues, where `alternate` is the alternate coding sequence with the 3' UTR after it:
/// `?` for none. `_stop_loss_extra_AA` translates those bases under the standard code, takes
/// their first stop wherever it stands and subtracts one and the peptide's length
/// (TranscriptVariationAllele.pm:2492, 2503 and 2513); it makes no count where the
/// description opens at the first residue (2483) or the figure is not above nought (2519).
fn count_ensembl_vep_writes_after_a_deletion(
    alternate: &str,
    peptide_length: usize,
    first_named: usize,
) -> String {
    let first_stop = translate_codons(alternate, 1).and_then(|residues| residues.find('*'));
    match first_stop {
        Some(place) if first_named > 1 && place > peptide_length => {
            (place - peptide_length).to_string()
        }
        _ => "?".to_string(),
    }
}

/// The eleventh shape of `hgvsp_stop_codon_protein_comparison` on a deletion, which the row
/// alone does not decide: both values name one protein and, with a pair of prediction
/// parentheses that both carry set aside, are one deletion, of one residue or of a range,
/// with an extension's count after it; vep-rs's has a count and Ensembl VEP's the unknown
/// count or another count; the row's terms hold `stop_lost`; and the reference row's
/// `CDS_position` (`perl_row`) and the row's transcript in the corpus cache give both counts.
/// The transcript's peptide (`peptides`) holds a stop and is the translation of its coding
/// sequence (`coding_sequences`) at each codon, a first residue written as a methionine over
/// another first codon set aside. The coding bases the row names are three for each residue
/// the description names. With them out of the coding sequence, and the 3' UTR read
/// after it, the first residue that differs is the description's first; the first stop of
/// the reference at or after it is a stop the peptide holds or the stop codon that closes
/// the coding sequence; and the first stop of the alternate at or after it stands vep-rs's
/// count of places after that stop's. Ensembl VEP's count is the one its own arithmetic
/// gives on the alternate bases (`count_ensembl_vep_writes_after_a_deletion`). The row's
/// `FLAGS` are not read. `false` without the cache's sequences and peptides.
pub fn deletion_counted_from_a_stop_the_peptide_holds_pair(
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    coding_sequences: Option<&CodingSequences>,
    peptides: Option<&HashMap<String, String>>,
) -> bool {
    if matches!(vep, "" | "-") || matches!(vep_rs, "" | "-") {
        return false;
    }
    let (Some((vep_name, reference)), Some((vep_rs_name, candidate))) =
        (described(vep), described(vep_rs))
    else {
        return false;
    };
    if vep_name != vep_rs_name || !holds_term(&key.consequences, "stop_lost") {
        return false;
    }
    let inside = |description: &str| -> Option<String> {
        let inside = description.strip_prefix('(')?.strip_suffix(')')?;
        Some(inside.to_string())
    };
    let (reference, candidate) = match (inside(&reference), inside(&candidate)) {
        (Some(reference), Some(candidate)) => (reference, candidate),
        _ => (reference, candidate),
    };
    let (Some(theirs), Some(counted)) = (
        DELETED_WITH_EXTENSION.captures(&reference),
        DELETED_WITH_EXTENSION.captures(&candidate),
    ) else {
        return false;
    };
    if theirs[1] != counted[1] || &counted[4] == "?" || theirs[4] == counted[4] {
        return false;
    }
    let number = |digits: &str| digits.parse::<usize>().ok();
    let Some(first_named) = number(&counted[2]) else {
        return false;
    };
    let last_named = counted
        .get(3)
        .map_or(Some(first_named), |last| number(last.as_str()));
    let named = last_named.and_then(|last| last.checked_add(1)?.checked_sub(first_named));
    let bounds = perl_row
        .and_then(|row| row.get("CDS_position"))
        .and_then(|column| column_bounds(column))
        .and_then(|(from, to)| Some((usize::try_from(from).ok()?, usize::try_from(to).ok()?)));
    let sequence = coding_sequences.and_then(|sequences| sequences.get(&key.feature));
    let peptide = peptides.and_then(|peptides| peptides.get(&key.feature));
    let (Some(named), Some((from, to)), Some(sequence), Some(peptide)) =
        (named, bounds, sequence, peptide)
    else {
        return false;
    };
    if !sequence.coding.is_ascii() || !sequence.utr.is_ascii() {
        return false;
    }
    let Some(reference) = translate_codons(&sequence.coding, sequence.codon_table) else {
        return false;
    };
    let (reference, peptide) = (reference.as_bytes(), peptide.as_bytes());
    let residues = peptide.len().min(reference.len());
    let fits = (0..residues).all(|residue| {
        reference[residue] == peptide[residue] || (residue == 0 && peptide[0] == b'M')
    });
    let whole_codons = 1 <= from
        && to <= sequence.coding.len()
        && named
            .checked_mul(3)
            .and_then(|bases| from.checked_add(bases))
            == Some(to + 1);
    if !peptide.contains(&b'*') || !fits || !whole_codons {
        return false;
    }
    let bases = format!(
        "{}{}{}",
        &sequence.coding[..from - 1],
        &sequence.coding[to..],
        sequence.utr
    );
    let Some(alternate) = translate_codons(&bases, sequence.codon_table) else {
        return false;
    };
    let alternate = alternate.as_bytes();
    let Some(first) = alternate
        .iter()
        .zip(reference)
        .position(|(alternate, reference)| alternate != reference)
    else {
        return false;
    };
    let stop_from = |protein: &[u8]| {
        protein[first..]
            .iter()
            .position(|&residue| residue == b'*')
            .map(|residues_on| first + residues_on)
    };
    let (Some(lost), Some(new)) = (stop_from(reference), stop_from(alternate)) else {
        return false;
    };
    let held = lost < residues && peptide[lost] == b'*';
    let closing = lost == peptide.len() && sequence.coding.len() == 3 * (lost + 1);
    first + 1 == first_named
        && (held || closing)
        && new > lost
        && (new - lost).to_string() == counted[4]
        && count_ensembl_vep_writes_after_a_deletion(&bases, peptide.len(), first_named)
            == theirs[4]
}

/// `hgvsp_stop_named_on_incomplete_coding_sequence` on an `HGVSp` the engines write
/// differently: vep-rs writes none; Ensembl VEP's description, an escaped `=` and a pair of
/// prediction parentheses set aside, names the stop alone as unchanged, deleted, replaced by
/// one residue or extended; and the
/// row shows the incomplete coding sequence: the reference row's `FLAGS` (`perl_row`) name
/// `cds_end_NF`, its terms hold `incomplete_terminal_codon_variant`, or the description is
/// `Ter1=` and its `FLAGS` name `cds_start_NF`. Beyond the row, the corpus cache shows it:
/// the row's transcript is one of `partial_codon_models`
/// (`partial_terminal_codon_models`).
pub fn stop_named_on_incomplete_coding_sequence_pair(
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    partial_codon_models: Option<&HashSet<String>>,
) -> bool {
    if matches!(vep, "" | "-") || !matches!(vep_rs, "" | "-") {
        return false;
    }
    let Some((_, reference)) = described(vep) else {
        return false;
    };
    let reference = reference
        .strip_prefix('(')
        .and_then(|inside| inside.strip_suffix(')'))
        .unwrap_or(&reference);
    let Some(named) = NAMED_ALONE.captures(reference) else {
        return false;
    };
    let flags = row_flags(perl_row);
    flags.contains("cds_end_NF")
        || holds_term(&key.consequences, "incomplete_terminal_codon_variant")
        || (&named[1] == "1" && &named[2] == "=" && flags.contains("cds_start_NF"))
        || partial_codon_models.is_some_and(|models| models.contains(&key.feature))
}

/// `hgvsp_frameshift_with_no_changed_residue` on an `HGVSp` the engines write differently:
/// vep-rs writes none; Ensembl VEP's description, a pair of prediction parentheses set
/// aside, is a frameshift with the unknown count whose two residue names are one name and
/// no stop's, the pair its walk reads last where the alternate translation ends before a
/// residue differs (TranscriptVariationAllele.pm:2348, :2351-2352, :2361 and :2364); and
/// the row's terms hold `frameshift_variant`.
pub fn frameshift_with_no_changed_residue_pair(key: &FieldKey, vep: &str, vep_rs: &str) -> bool {
    if matches!(vep, "" | "-") || !matches!(vep_rs, "" | "-") {
        return false;
    }
    if !holds_term(&key.consequences, "frameshift_variant") {
        return false;
    }
    let Some((_, reference)) = described(vep) else {
        return false;
    };
    let reference = reference
        .strip_prefix('(')
        .and_then(|inside| inside.strip_suffix(')'))
        .unwrap_or(&reference);
    FRAMESHIFT
        .captures(reference)
        .is_some_and(|named| &named[4] == "?" && named[1] == named[3] && &named[1] != "Ter")
}

/// The class `hgvsp_stop_codon_protein_comparison` of the dispatch.
pub fn stop_codon_protein_comparison(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
) -> Verdict {
    if question.name != "HGVSp" {
        return Verdict::Pass;
    }
    let (key, vep, vep_rs, perl_row) = (
        question.key,
        question.vep,
        question.vep_rs,
        question.perl_row,
    );
    if stop_codon_protein_comparison_pair(key, vep, vep_rs, perl_row)
        || inserted_residues_keep_the_protein_pair(
            key,
            vep,
            vep_rs,
            perl_row,
            models.protein.coding_sequences.as_ref(),
        )
        || deletion_counted_from_a_stop_the_peptide_holds_pair(
            key,
            vep,
            vep_rs,
            perl_row,
            models.protein.coding_sequences.as_ref(),
            models.protein.held_stop_peptides.as_ref(),
        )
    {
        return Verdict::Class(HGVSP_STOP_CODON_PROTEIN_COMPARISON, None);
    }
    Verdict::Pass
}

/// The class `hgvsp_stop_named_on_incomplete_coding_sequence` of the dispatch.
pub fn stop_named_on_incomplete_coding_sequence(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
) -> Verdict {
    if question.name == "HGVSp"
        && stop_named_on_incomplete_coding_sequence_pair(
            question.key,
            question.vep,
            question.vep_rs,
            question.perl_row,
            models.protein.partial_codon_models.as_ref(),
        )
    {
        return Verdict::Class(HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE, None);
    }
    Verdict::Pass
}

/// The class `hgvsp_frameshift_with_no_changed_residue` of the dispatch.
pub fn frameshift_with_no_changed_residue(
    question: &FieldQuestion<'_>,
    _models: &FieldModels,
) -> Verdict {
    if question.name == "HGVSp"
        && frameshift_with_no_changed_residue_pair(question.key, question.vep, question.vep_rs)
    {
        return Verdict::Class(HGVSP_FRAMESHIFT_WITH_NO_CHANGED_RESIDUE, None);
    }
    Verdict::Pass
}

/// A cell of a row by its name, `-` for a cell the row lacks.
fn row_cell<'a>(row: &'a FieldRow, name: &str) -> &'a str {
    row.get(name).map_or("-", String::as_str)
}

/// The five cells of a row that a derivation of `crate::rules::coding_ends::terms` reads.
fn row_cells(row: &FieldRow) -> Cells<'_> {
    Cells {
        cdna_position: row_cell(row, "cDNA_position"),
        cds_position: row_cell(row, "CDS_position"),
        protein_position: row_cell(row, "Protein_position"),
        amino_acids: row_cell(row, "Amino_acids"),
        codons: row_cell(row, "Codons"),
    }
}

/// `hgvsp_frameshift_count_at_the_first_residue` on an `HGVSp` the engines write
/// differently. The rule is `frameshift_counted_at_the_first_residue` of
/// `crate::rules::coding_ends::terms`, the source the comparison tools compile, asked with
/// what `classify` reads of the row. Both values name one protein and a frameshift at
/// residue 1, Ensembl VEP's with the unknown count and vep-rs's with a number
/// (`count_alone_differs_at_the_first_residue`). The record of the row's Feature in the
/// corpus cache whose span shares a base with the row's Location (`Models::coding_ends`)
/// shows no start codon, with `cds_start_NF` or without it. The derivation reads the key's
/// Allele and its consequence set, which holds no `start_lost`; the `CDS_position` of both
/// rows; the `HGVS_OFFSET` of both rows, which they carry where the run compares that
/// field, for a deletion whose `CDS_position` is open at its start; and the count of
/// reference bases the Allele stands in for, which the run's records give by position
/// where they read the row alike (`record_reference` under no identifier). The pair is the
/// class's where the two values are the two texts derived: vep-rs's with the count of
/// residues from residue 1 to the first stop codon of the new reading frame, Ensembl VEP's
/// the same text with the unknown count. `false` without the cache's coding ends, without
/// either row, on a sequence whose shards are not readable, where no one record answers
/// the row, and where the reading keeps no sequence of the transcript, as on a transcript
/// with a 3' UTR.
pub fn frameshift_count_at_the_first_residue_pair(
    key: &FieldKey,
    vep: &str,
    vep_rs: &str,
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
    models: &Models,
) -> bool {
    if !count_alone_differs_at_the_first_residue(vep, vep_rs) {
        return false;
    }
    let (Some(transcripts), Some(theirs), Some(ours), Some(location)) = (
        models.coding_ends.as_ref(),
        perl_row,
        rust_row,
        location_span(&key.location),
    ) else {
        return false;
    };
    let Ok(Some(ends)) = transcripts.coding_ends(&key.location, &key.feature) else {
        return false;
    };
    let reference = record_reference(&models.records, &key.location, &key.allele, "");
    let row = RowFacts {
        transcript_row: transcript_feature_type(&key.feature_type),
        allele: &key.allele,
        location,
        ensembl_vep_terms: &key.consequences,
    };
    let read = WhatWasRead {
        cells: Some((row_cells(theirs), row_cells(ours))),
        reference_length: reference.map(|bases| bases.chars().count() as u64),
        reference_bases: None,
        // A row's reference bases are read from a VCF record, which is given on the forward
        // strand.
        record_on_the_reverse_strand: false,
    };
    let offsets = (
        row_cell(theirs, "HGVS_OFFSET"),
        row_cell(ours, "HGVS_OFFSET"),
    );
    let counted = frameshift_counted_at_the_first_residue(
        &row,
        &key.consequences,
        Some(ends.as_ref()),
        &read,
        offsets,
    );
    counted.is_ok_and(|counted| counted.takes(vep, vep_rs))
}

/// The class `hgvsp_frameshift_count_at_the_first_residue` of the dispatch.
pub fn frameshift_count_at_the_first_residue(
    question: &FieldQuestion<'_>,
    models: &FieldModels,
) -> Verdict {
    if question.name == "HGVSp"
        && frameshift_count_at_the_first_residue_pair(
            question.key,
            question.vep,
            question.vep_rs,
            question.perl_row,
            question.rust_row,
            &models.protein,
        )
    {
        return Verdict::Class(HGVSP_FRAMESHIFT_COUNT_AT_THE_FIRST_RESIDUE, None);
    }
    Verdict::Pass
}

#[cfg(test)]
#[path = "protein_tests.rs"]
mod tests;
