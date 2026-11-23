// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HGVSp stop-codon classes.

use crate::consts::AMINO_ACIDS_FIELDS;
use crate::entry::Entry;
use crate::error::Result;
use crate::genetic_codes::genetic_code;
use crate::genetic_codes::CodonTable;
use crate::genetic_codes::GeneticCode;
use crate::genetic_codes::Reading;
use crate::genetic_codes::STANDARD_CODE;
use crate::pyfmt;
use crate::rules::downstream::ambiguous_residue;
use crate::rules::downstream::field_of;
use crate::rules::downstream::protein_span;
use crate::rules::downstream::translate_until_stop;
use crate::rules::downstream::DownstreamSequences;
use crate::rules::downstream::DownstreamTranscript;
use crate::rules::term_pairs::cds_bounds;
use crate::rules::whole_text;
use crate::shared::blank;
use crate::shared::consequence_set;
use crate::shared::integer_text_of;
use crate::shared::STOP_RETAINED_TERM;
use std::cmp::Ordering;
use std::sync::Arc;
use std::sync::LazyLock;

pub const RESIDUE_NAME: &str = "[A-Z][a-z]{2}";

pub const STOP_UNCHANGED_RE_PATTERN: &str = r"Ter([1-9]\d*)=";
pub static STOP_UNCHANGED_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(STOP_UNCHANGED_RE_PATTERN).expect("the pattern compiles"));

pub const STOP_EXTENSION_RE_PATTERN: &str = r"Ter([1-9]\d*)([A-Z][a-z]{2})extTer([1-9]\d*|\?)";
pub static STOP_EXTENSION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(STOP_EXTENSION_RE_PATTERN).expect("the pattern compiles"));

pub const STOP_REPLACED_RE_PATTERN: &str = r"Ter([1-9]\d*)delins((?:[A-Z][a-z]{2})+)";
pub static STOP_REPLACED_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(STOP_REPLACED_RE_PATTERN).expect("the pattern compiles"));

pub const STOP_INSERTION_RE_PATTERN: &str =
    r"[A-Z][a-z]{2}([1-9]\d*)_Ter([1-9]\d*)ins((?:[A-Z][a-z]{2})+)";
pub static STOP_INSERTION_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(STOP_INSERTION_RE_PATTERN).expect("the pattern compiles"));

pub const STOP_RANGE_REPLACED_RE_PATTERN: &str =
    r"([A-Z][a-z]{2})([1-9]\d*)_Ter([1-9]\d*)delins((?:[A-Z][a-z]{2})+)";
pub static STOP_RANGE_REPLACED_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(STOP_RANGE_REPLACED_RE_PATTERN).expect("the pattern compiles")
});

pub const STOP_FRAMESHIFT_RE_PATTERN: &str =
    r"([A-Z][a-z]{2})([1-9]\d*)([A-Z][a-z]{2})fsTer([1-9]\d*|\?)";
pub static STOP_FRAMESHIFT_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(STOP_FRAMESHIFT_RE_PATTERN).expect("the pattern compiles"));

pub const STOP_DELETED_EXTENSION_RE_PATTERN: &str = r"Ter([1-9]\d*)delextTer([1-9]\d*|\?)";
pub static STOP_DELETED_EXTENSION_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(STOP_DELETED_EXTENSION_RE_PATTERN).expect("the pattern compiles")
});

/// A deletion of one residue or of a range of residues with an extension's count after it:
/// the deletion as written, its first residue's number, its last residue's number where it
/// names a range, and the count.
const DELETION_EXTENSION_RE_PATTERN: &str =
    r"([A-Z][a-z]{2}([1-9]\d*)(?:_[A-Z][a-z]{2}([1-9]\d*))?del)extTer([1-9]\d*|\?)";

pub const RESIDUE_AND_STOP_INSERTED_RE_PATTERN: &str =
    r"([A-Z][a-z]{2})([1-9]\d*)_([A-Z][a-z]{2})([1-9]\d*)ins([A-Z][a-z]{2})Ter";
pub static RESIDUE_AND_STOP_INSERTED_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(RESIDUE_AND_STOP_INSERTED_RE_PATTERN).expect("the pattern compiles")
});

pub const RESIDUES_AND_STOP_INSERTED_RE_PATTERN: &str =
    r"([A-Z][a-z]{2})([1-9]\d*)_([A-Z][a-z]{2})([1-9]\d*)ins((?:[A-Z][a-z]{2}){2,})Ter";
pub static RESIDUES_AND_STOP_INSERTED_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(RESIDUES_AND_STOP_INSERTED_RE_PATTERN).expect("the pattern compiles")
});

pub const STOP_NAMED_ALONE_RE_PATTERN: &str =
    r"Ter([1-9]\d*)(=|del|[A-Z][a-z]{2}(?:extTer(?:[1-9]\d*|\?))?)";
pub static STOP_NAMED_ALONE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(STOP_NAMED_ALONE_RE_PATTERN).expect("the pattern compiles"));

static STOP_UNCHANGED_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_UNCHANGED_RE_PATTERN));
static STOP_EXTENSION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_EXTENSION_RE_PATTERN));
static STOP_REPLACED_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_REPLACED_RE_PATTERN));
static STOP_INSERTION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_INSERTION_RE_PATTERN));
static STOP_RANGE_REPLACED_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_RANGE_REPLACED_RE_PATTERN));
static STOP_FRAMESHIFT_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_FRAMESHIFT_RE_PATTERN));
static STOP_DELETED_EXTENSION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_DELETED_EXTENSION_RE_PATTERN));
static DELETION_EXTENSION_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(DELETION_EXTENSION_RE_PATTERN));
static RESIDUE_AND_STOP_INSERTED_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(RESIDUE_AND_STOP_INSERTED_RE_PATTERN));
static RESIDUES_AND_STOP_INSERTED_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(RESIDUES_AND_STOP_INSERTED_RE_PATTERN));
static STOP_NAMED_ALONE_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(STOP_NAMED_ALONE_RE_PATTERN));

/// `int(digits)` of a run of digits a pattern matched, as its decimal text: the number at
/// any length the interpreter converts. A longer run ends the run as the interpreter's
/// `ValueError` does.
fn whole_number(digits: &str) -> String {
    match pyfmt::int_checked(digits) {
        Ok(number) => number,
        Err(error) => crate::uncaught!("{error}"),
    }
}

/// The order of two whole numbers that are not negative, each its decimal text.
fn compared(a: &str, b: &str) -> Ordering {
    a.len().cmp(&b.len()).then_with(|| a.cmp(b))
}

/// The sum of two whole numbers that are not negative, each its decimal text.
fn sum(a: &str, b: &str) -> String {
    let (mut a, mut b) = (a.bytes().rev(), b.bytes().rev());
    let mut digits = Vec::new();
    let mut carried = 0;
    loop {
        let (x, y) = (a.next(), b.next());
        if x.is_none() && y.is_none() && carried == 0 {
            break;
        }
        let total = x.map_or(0, |digit| digit - b'0') + y.map_or(0, |digit| digit - b'0') + carried;
        digits.push(char::from(b'0' + total % 10));
        carried = total / 10;
    }
    digits.iter().rev().collect()
}

/// Whether the row's terms hold `term`.
fn holds_term(consequence: &str, term: &str) -> bool {
    consequence_set(consequence)
        .split(',')
        .any(|held| held == term)
}

/// The two descriptions of a differing HGVSp whose values name one protein, the reference's
/// first, an escaped `=` and a pair of prediction parentheses that both carry set aside;
/// `None` where a value is blank, the values name two proteins or none, or the row's FLAGS
/// (`flags`) name `cds_end_NF`.
fn descriptions_of_one_protein(flags: &str, ev: &str, av: &str) -> Option<(String, String)> {
    if blank(ev) || blank(av) || flags.contains("cds_end_NF") {
        return None;
    }
    let (reference_name, separator, reference) = pyfmt::partition(ev, ":p.");
    let (candidate_name, _, candidate) = pyfmt::partition(av, ":p.");
    if separator.is_empty() || reference_name != candidate_name {
        return None;
    }
    let (reference, candidate) = (reference.replace("%3D", "="), candidate.replace("%3D", "="));
    let predicted = |description: &str| description.starts_with('(') && description.ends_with(')');
    if predicted(&reference) && predicted(&candidate) {
        return Some((
            reference[1..reference.len() - 1].to_string(),
            candidate[1..candidate.len() - 1].to_string(),
        ));
    }
    Some((reference, candidate))
}

/// `hgvsp_stop_codon_protein_comparison` on a differing HGVSp: both values name one
/// protein, the row's FLAGS (`flags`) name no `cds_end_NF`, and the two descriptions, an
/// escaped `=` and a pair of prediction parentheses set aside, are one of the class's first
/// nine shapes; `consequence` is the reference row's terms, which the last five of them read:
/// `stop_lost` in the fifth to the eighth, `stop_retained_variant` in the ninth.
pub fn stop_codon_protein_comparison(
    consequence: &str,
    flags: &str,
    ev: &str,
    av: &str,
) -> Result<bool> {
    let Some((reference, candidate)) = descriptions_of_one_protein(flags, ev, av) else {
        return Ok(false);
    };
    let (reference, candidate) = (reference.as_str(), candidate.as_str());
    Ok(stop_replaced_by_a_stop(reference, candidate)
        || stop_replaced_by_first_residues(reference, candidate)
        || stop_replaced_through_an_inserted_stop(reference, candidate)
        || residues_inserted_before_the_stop(reference, candidate)
        || lost_stop_written_as_unchanged(consequence, reference, candidate)
        || range_to_the_stop_replaced(consequence, reference, candidate)
        || unknown_count_where_a_stop_is_met(consequence, reference, candidate)
        || stop_deleted_with_an_extension_count(consequence, reference, candidate)
        || residue_and_stop_inserted_before_the_last_residue(consequence, reference, candidate)?)
}

/// The three-letter residue names a run of them holds, in order.
pub fn residue_names(text: &str) -> Vec<String> {
    text.chars()
        .collect::<Vec<char>>()
        .chunks(3)
        .map(|name| name.iter().collect())
        .collect()
}

/// The residue number, first new residue and count (`?` for an unknown one) of an
/// extension named at the stop of a protein of at least one residue; `None` for any other
/// description.
pub fn stop_extension(description: &str) -> Option<(String, String, String)> {
    let m = STOP_EXTENSION_WHOLE.captures(description)?;
    if &m[1] == "1" || &m[2] == "Ter" {
        return None;
    }
    Some((m[1].to_string(), m[2].to_string(), m[3].to_string()))
}

/// The stop replaced by a stop against the stop written as unchanged at the same residue.
pub fn stop_replaced_by_a_stop(reference: &str, candidate: &str) -> bool {
    let (Some(replaced), Some(unchanged)) = (
        STOP_REPLACED_WHOLE.captures(reference),
        STOP_UNCHANGED_WHOLE.captures(candidate),
    ) else {
        return false;
    };
    &replaced[2] == "Ter" && replaced[1] == unchanged[1] && &unchanged[1] != "1"
}

/// The stop alone replaced by the first residues that follow the protein, closed by `Ter`
/// or not, against an extension at the same residue with the first of them: its count unknown,
/// above the residues the reference names, or their number where the reference closes with no
/// `Ter`.
pub fn stop_replaced_by_first_residues(reference: &str, candidate: &str) -> bool {
    let (Some(replaced), Some((number, first, count))) = (
        STOP_REPLACED_WHOLE.captures(reference),
        stop_extension(candidate),
    ) else {
        return false;
    };
    let mut named = residue_names(&replaced[2]);
    let closed = named.last().is_some_and(|name| name == "Ter");
    if closed {
        named.pop();
    }
    if named.is_empty()
        || named.iter().any(|name| name == "Ter")
        || replaced[1] != number
        || named[0] != first
    {
        return false;
    }
    if count == "?" {
        return true;
    }
    match compared(&whole_number(&count), &named.len().to_string()) {
        Ordering::Greater => true,
        Ordering::Equal => !closed,
        Ordering::Less => false,
    }
}

/// The stop replaced by every added residue and then a stop against an extension at the
/// same residue with the first of them, its count the number of residues the reference names.
pub fn stop_replaced_through_an_inserted_stop(reference: &str, candidate: &str) -> bool {
    let (Some(replaced), Some((number, first, count))) = (
        STOP_REPLACED_WHOLE.captures(reference),
        stop_extension(candidate),
    ) else {
        return false;
    };
    let mut named = residue_names(&replaced[2]);
    if named.pop().as_deref() != Some("Ter")
        || named.is_empty()
        || named.iter().any(|name| name == "Ter")
    {
        return false;
    }
    replaced[1] == number && named[0] == first && count == named.len().to_string()
}

/// Residues inserted between the last residue and the stop against an extension at the stop
/// with the first of them, its count their number.
pub fn residues_inserted_before_the_stop(reference: &str, candidate: &str) -> bool {
    let (Some(inserted), Some((number, first, count))) = (
        STOP_INSERTION_WHOLE.captures(reference),
        stop_extension(candidate),
    ) else {
        return false;
    };
    let named = residue_names(&inserted[3]);
    let (last_residue, stop) = (whole_number(&inserted[1]), whole_number(&inserted[2]));
    sum(&last_residue, "1") == stop
        && inserted[2] == number
        && !named.iter().any(|name| name == "Ter")
        && named[0] == first
        && count == named.len().to_string()
}

/// The stop written as unchanged against an extension at the same residue, on a row whose
/// terms hold `stop_lost`.
pub fn lost_stop_written_as_unchanged(consequence: &str, reference: &str, candidate: &str) -> bool {
    let (Some(unchanged), Some((number, _, _))) = (
        STOP_UNCHANGED_WHOLE.captures(reference),
        stop_extension(candidate),
    ) else {
        return false;
    };
    unchanged[1] == number && holds_term(consequence, "stop_lost")
}

/// A range from a residue to the stop replaced by residues, none of them a stop, against a
/// frameshift named at that residue with the first of them, on a row whose terms hold
/// `stop_lost`. The frameshift's count is unknown, or above their number and not the one
/// that leaves the new stop at the stop's own residue: a protein that keeps its length is no
/// frameshift of the candidate's.
pub fn range_to_the_stop_replaced(consequence: &str, reference: &str, candidate: &str) -> bool {
    let (Some(replaced), Some(frameshift)) = (
        STOP_RANGE_REPLACED_WHOLE.captures(reference),
        STOP_FRAMESHIFT_WHOLE.captures(candidate),
    ) else {
        return false;
    };
    if !holds_term(consequence, "stop_lost") {
        return false;
    }
    let (residue, number, stop) = (&replaced[1], &replaced[2], &replaced[3]);
    let named = residue_names(&replaced[4]);
    let count = &frameshift[4];
    if residue == "Ter" || named.iter().any(|name| name == "Ter") {
        return false;
    }
    let (first_replaced, replaced_stop) = (whole_number(number), whole_number(stop));
    if compared(&first_replaced, &replaced_stop) != Ordering::Less
        || &frameshift[1] != residue
        || &frameshift[2] != number
        || frameshift[3] != named[0]
    {
        return false;
    }
    if count == "?" {
        return true;
    }
    let counted = whole_number(count);
    // The count is the residues from the range's first to the stop where its sum with the
    // first is one more than the stop.
    compared(&counted, &named.len().to_string()) == Ordering::Greater
        && sum(&counted, &first_replaced) != sum(&replaced_stop, "1")
}

/// `hgvsp_stop_named_on_incomplete_coding_sequence` on a differing HGVSp: the candidate
/// writes none; the reference's description, an escaped `=` and a pair of prediction
/// parentheses set aside, names the stop alone as unchanged, deleted, replaced by one residue
/// or extended; and the row shows the incomplete coding sequence: its FLAGS (`flags`) name
/// `cds_end_NF`, its terms (`consequence`) hold `incomplete_terminal_codon_variant`, or
/// the description is `Ter1=` and its FLAGS name `cds_start_NF`.
pub fn stop_named_on_incomplete_coding_sequence(
    consequence: &str,
    flags: &str,
    ev: &str,
    av: &str,
) -> bool {
    let Some((residue, change)) = stop_named_alone_against_no_value(ev, av) else {
        return false;
    };
    flags.contains("cds_end_NF")
        || holds_term(consequence, "incomplete_terminal_codon_variant")
        || (residue == "1" && change == "=" && flags.contains("cds_start_NF"))
}

/// The residue number and the change of a reference description that names the stop alone,
/// as unchanged, deleted, replaced by one residue or extended, an escaped `=` and a pair of
/// prediction parentheses set aside, where the candidate writes none; `None` for any other
/// two values.
fn stop_named_alone_against_no_value(ev: &str, av: &str) -> Option<(String, String)> {
    if blank(ev) || !blank(av) {
        return None;
    }
    let (_, separator, reference) = pyfmt::partition(ev, ":p.");
    let reference = reference.replace("%3D", "=");
    let mut reference = reference.as_str();
    if reference.starts_with('(') && reference.ends_with(')') {
        reference = &reference[1..reference.len() - 1];
    }
    let named = STOP_NAMED_ALONE_WHOLE.captures(reference)?;
    if separator.is_empty() {
        return None;
    }
    Some((named[1].to_string(), named[2].to_string()))
}

/// `hgvsp_stop_named_on_incomplete_coding_sequence` on a differing HGVSp of a row that does
/// not itself show the incomplete coding sequence: the candidate writes none; the
/// reference's description names the stop alone in one of the class's four ways; and the
/// candidate's JSON cache (`models`) holds the row's transcript with a coding sequence that
/// ends inside a codon. A frameshift that the 3' shift carries past the last whole codon is
/// such a row: `_get_fs_peptides` closes the reference translation with a stop whatever the
/// coding sequence's end (`TranscriptVariationAllele.pm:2338`), the row's terms are the
/// written allele's, and its FLAGS name no `cds_end_NF` where the transcript carries none.
pub fn stop_named_on_a_transcript_of_a_partial_last_codon(
    e: &crate::entry::Entry,
    ev: &str,
    av: &str,
    models: &crate::rules::partial_codon_models::PartialCodonModels,
) -> Result<bool> {
    if stop_named_alone_against_no_value(ev, av).is_none() {
        return Ok(false);
    }
    models.holds(e)
}

/// An extension with an unknown count against the same extension, at the same residue and
/// with the same first residue, with a count, on a row whose terms hold `stop_lost`.
pub fn unknown_count_where_a_stop_is_met(
    consequence: &str,
    reference: &str,
    candidate: &str,
) -> bool {
    let (Some(unknown), Some(counted)) = (stop_extension(reference), stop_extension(candidate))
    else {
        return false;
    };
    unknown.2 == "?"
        && counted.2 != "?"
        && (&unknown.0, &unknown.1) == (&counted.0, &counted.1)
        && holds_term(consequence, "stop_lost")
}

/// The stop written as deleted with an extension count and no residue in the stop's place,
/// against an extension at the same residue with the same count, on a row whose terms hold
/// `stop_lost`.
pub fn stop_deleted_with_an_extension_count(
    consequence: &str,
    reference: &str,
    candidate: &str,
) -> bool {
    let (Some(deleted), Some((number, _, count))) = (
        STOP_DELETED_EXTENSION_WHOLE.captures(reference),
        stop_extension(candidate),
    ) else {
        return false;
    };
    deleted[1] == number && deleted[2] == count && holds_term(consequence, "stop_lost")
}

/// One residue and a stop inserted between two neighbouring residues, the residue the right
/// flank's, against the stop written as unchanged at the residue after the right flank, on a
/// row whose terms hold `stop_retained_variant`.
pub fn residue_and_stop_inserted_before_the_last_residue(
    consequence: &str,
    reference: &str,
    candidate: &str,
) -> Result<bool> {
    let (Some(inserted), Some(unchanged)) = (
        RESIDUE_AND_STOP_INSERTED_WHOLE.captures(reference),
        STOP_UNCHANGED_WHOLE.captures(candidate),
    ) else {
        return Ok(false);
    };
    if !holds_term(consequence, STOP_RETAINED_TERM) {
        return Ok(false);
    }
    let (left, right, residue) = (&inserted[1], &inserted[3], &inserted[5]);
    let left_number = integer_text_of(&inserted[2], "HGVSp residue number")?;
    let right_number = integer_text_of(&inserted[4], "HGVSp residue number")?;
    let stop = integer_text_of(&unchanged[1], "HGVSp residue number")?;
    Ok(left != "Ter"
        && right != "Ter"
        && residue == right
        && sum(&left_number, "1") == right_number
        && sum(&right_number, "1") == stop)
}

/// The three-letter name `hgvs_protein` writes for one of the twenty residues of the
/// standard alphabet, by its letter; `None` for any other letter.
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
/// residue: its codons translated under the transcript's own genetic code, as the
/// candidate reads that code, up to the stop codon that closes the sequence. `None` under a
/// code no table is held for and where no stop codon closes the sequence: it is no whole
/// number of codons, or its first stop codon is not its last codon.
pub fn protein_closed_by_a_stop(transcript: &DownstreamTranscript) -> Option<String> {
    let CodonTable::Held(table) = &transcript.table else {
        return None;
    };
    let residues = translate_until_stop(&transcript.coding, *table, Reading::VepRs);
    let bases = transcript.coding.chars().count();
    (bases == 3 * (residues.len() + 1)).then_some(residues)
}

/// Two or more residues and a stop inserted between two neighbouring residues, against the
/// stop written as unchanged, where `protein`, the reference protein with one letter a
/// residue, shows the protein kept: the left flank is the residue of its number, the right
/// flank and the inserted residues are the protein's from the next number to its last, and
/// the unchanged stop stands at the number after the last. By the reference's value the
/// alternate protein is then the reference protein and a stop at its own residue.
pub fn residues_and_stop_inserted_that_keep_the_protein(
    reference: &str,
    candidate: &str,
    protein: &str,
) -> Result<bool> {
    let (Some(inserted), Some(unchanged)) = (
        RESIDUES_AND_STOP_INSERTED_WHOLE.captures(reference),
        STOP_UNCHANGED_WHOLE.captures(candidate),
    ) else {
        return Ok(false);
    };
    let left_number = integer_text_of(&inserted[2], "HGVSp residue number")?;
    let right_number = integer_text_of(&inserted[4], "HGVSp residue number")?;
    let stop = integer_text_of(&unchanged[1], "HGVSp residue number")?;
    let named = residue_names(&inserted[5]);
    let residues = protein.as_bytes();
    // The count of residues before the inserted ones, which is the left flank's number.
    let Some(kept) = residues.len().checked_sub(named.len()) else {
        return Ok(false);
    };
    if kept == 0
        || left_number != kept.to_string()
        || sum(&left_number, "1") != right_number
        || stop != (residues.len() + 1).to_string()
    {
        return Ok(false);
    }
    let is_residue = |name: &str, at: usize| residue_name(residues[at]) == Some(name);
    Ok(is_residue(&inserted[1], kept - 1)
        && is_residue(&inserted[3], kept)
        && named
            .iter()
            .enumerate()
            .all(|(place, name)| is_residue(name, kept + place)))
}

/// The ninth shape of `hgvsp_stop_codon_protein_comparison` with more than one residue
/// before the inserted stop, which the row alone does not decide: both values name one
/// protein, the row's FLAGS (`flags`) name no `cds_end_NF`, and the coding sequence of the
/// row's transcript in the candidate's cache (`sequences`) shows the protein kept
/// (`residues_and_stop_inserted_that_keep_the_protein`). The row's terms are not read. The
/// cache is read for a pair of the shape alone.
pub fn stop_unchanged_where_inserted_residues_keep_the_protein(
    e: &Entry,
    flags: &str,
    ev: &str,
    av: &str,
    sequences: &DownstreamSequences,
) -> Result<bool> {
    let Some((reference, candidate)) = descriptions_of_one_protein(flags, ev, av) else {
        return Ok(false);
    };
    if !RESIDUES_AND_STOP_INSERTED_WHOLE.is_match(&reference)
        || !STOP_UNCHANGED_WHOLE.is_match(&candidate)
    {
        return Ok(false);
    }
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(false);
    };
    let Some(protein) = protein_closed_by_a_stop(&transcript) else {
        return Ok(false);
    };
    residues_and_stop_inserted_that_keep_the_protein(&reference, &candidate, &protein)
}

/// The residue each whole codon of `sequence` spells under `code`: a codon with an
/// ambiguous letter by the amino acids its readings code.
fn residues_spelled(sequence: &str, code: GeneticCode) -> Vec<u8> {
    sequence
        .as_bytes()
        .chunks_exact(3)
        .map(|codon| {
            code.amino_acid(codon)
                .unwrap_or_else(|| ambiguous_residue(code, codon))
        })
        .collect()
}

/// The count Ensembl VEP writes for a substitution that loses the stop at place `at` of the
/// translation of `bases`, a coding sequence with its 3' UTR, counted from nought, on a
/// peptide of `peptide_length` residues: `?` for none. `_stop_loss_extra_AA` translates the
/// alternate bases under the standard code, takes their first stop wherever it stands and
/// subtracts one and the peptide's length (TranscriptVariationAllele.pm:2492, 2503 and
/// 2513); it makes no count where the lost stop is the first residue (2483) or the figure is
/// not above nought (2519).
fn count_ensembl_vep_writes(bases: &str, peptide_length: usize, at: usize) -> String {
    let first_stop = genetic_code(STANDARD_CODE, Reading::EnsemblVep).and_then(|standard| {
        residues_spelled(bases, standard)
            .iter()
            .enumerate()
            .find(|(place, residue)| **residue == b'*' && *place != at)
            .map(|(place, _)| place)
    });
    match first_stop {
        Some(place) if at > 0 && place > peptide_length => (place - peptide_length).to_string(),
        _ => "?".to_string(),
    }
}

/// Whether the sequences of `transcript` give `counted` residues from the stop of residue
/// number `lost` to the next stop of the frame, and `theirs` as the count Ensembl VEP writes
/// (`count_ensembl_vep_writes`):
///
/// - the transcript's peptide holds a stop and is the protein of its coding sequence: at
///   each codon it has the codon's residue, its first residue a methionine whatever the
///   first codon;
/// - residue `lost` is a stop the peptide holds, or the stop codon that closes the coding
///   sequence, one past the peptide's last residue;
/// - the next stop codon after it, sought through the coding sequence and the 3' UTR and no
///   farther, stands `counted` residues on.
fn sequence_counts_from_the_lost_stop(
    transcript: &DownstreamTranscript,
    code: GeneticCode,
    lost: usize,
    (theirs, counted): (&str, &str),
) -> bool {
    if !transcript.coding.is_ascii() || !transcript.utr.is_ascii() {
        return false;
    }
    let Some(peptide) = transcript.peptide.as_deref().map(str::as_bytes) else {
        return false;
    };
    let Some(at) = lost.checked_sub(1) else {
        return false;
    };
    let reference = residues_spelled(&transcript.coding, code);
    let residues = peptide.len().min(reference.len());
    let fits = (0..residues).all(|residue| {
        reference[residue] == peptide[residue] || (residue == 0 && peptide[0] == b'M')
    });
    let held = at < residues && peptide[at] == b'*';
    let closing = at == peptide.len()
        && transcript.coding.len() == 3 * (at + 1)
        && reference.get(at) == Some(&b'*');
    if !peptide.contains(&b'*') || !fits || !(held || closing) {
        return false;
    }
    let bases = format!("{}{}", transcript.coding, transcript.utr);
    count_ensembl_vep_writes(&bases, peptide.len(), at) == theirs
        && residues_spelled(&bases, code)
            .get(at + 1..)
            .and_then(|after| after.iter().position(|&residue| residue == b'*'))
            .is_some_and(|residues_on| (residues_on + 1).to_string() == counted)
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
    let first_stop = genetic_code(STANDARD_CODE, Reading::EnsemblVep).and_then(|standard| {
        residues_spelled(alternate, standard)
            .iter()
            .position(|&residue| residue == b'*')
    });
    match first_stop {
        Some(place) if first_named > 1 && place > peptide_length => {
            (place - peptide_length).to_string()
        }
        _ => "?".to_string(),
    }
}

/// Whether the sequences of `transcript` give `counted` as the count of a deletion of the
/// coding bases `from` to `to` whose description names `named` residues from residue number
/// `first_named`, and `theirs` as the count Ensembl VEP writes
/// (`count_ensembl_vep_writes_after_a_deletion`):
///
/// - the transcript's peptide holds a stop and is the protein of its coding sequence: at
///   each codon it has the codon's residue, its first residue a methionine whatever the
///   first codon;
/// - the bases lie inside the coding sequence, three for each residue the description
///   names, and with them out the first residue that differs is `first_named`;
/// - the first stop of the reference at or after that residue is a stop the peptide holds,
///   or the stop codon that closes the coding sequence, one past the peptide's last residue;
/// - the first stop the alternate bases meet at or after that residue, sought through the
///   coding sequence and the 3' UTR and no farther, stands `counted` places after the lost
///   stop's place.
fn sequence_counts_a_deletion_from_the_lost_stop(
    transcript: &DownstreamTranscript,
    code: GeneticCode,
    (from, to): (usize, usize),
    (first_named, named): (usize, usize),
    (theirs, counted): (&str, &str),
) -> bool {
    if !transcript.coding.is_ascii() || !transcript.utr.is_ascii() {
        return false;
    }
    let Some(peptide) = transcript.peptide.as_deref().map(str::as_bytes) else {
        return false;
    };
    let reference = residues_spelled(&transcript.coding, code);
    let residues = peptide.len().min(reference.len());
    let fits = (0..residues).all(|residue| {
        reference[residue] == peptide[residue] || (residue == 0 && peptide[0] == b'M')
    });
    let whole_codons = 1 <= from
        && to <= transcript.coding.len()
        && named
            .checked_mul(3)
            .and_then(|bases| from.checked_add(bases))
            == Some(to + 1);
    if !peptide.contains(&b'*') || !fits || !whole_codons {
        return false;
    }
    let bases = format!(
        "{}{}{}",
        &transcript.coding[..from - 1],
        &transcript.coding[to..],
        transcript.utr
    );
    let alternate = residues_spelled(&bases, code);
    let Some(first) = alternate
        .iter()
        .zip(&reference)
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
    let (Some(lost), Some(new)) = (stop_from(&reference), stop_from(&alternate)) else {
        return false;
    };
    let held = lost < residues && peptide[lost] == b'*';
    let closing = lost == peptide.len() && transcript.coding.len() == 3 * (lost + 1);
    first + 1 == first_named
        && (held || closing)
        && new > lost
        && (new - lost).to_string() == counted
        && count_ensembl_vep_writes_after_a_deletion(&bases, peptide.len(), first_named) == theirs
}

/// The row's transcript in the candidate's cache (`sequences`) and the genetic code vep-rs
/// reads it under; `None` where the cache holds no such transcript or names a table the
/// tool holds no code for.
fn transcript_under_its_code(
    e: &Entry,
    sequences: &DownstreamSequences,
) -> Result<Option<(Arc<DownstreamTranscript>, GeneticCode)>> {
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(None);
    };
    let CodonTable::Held(table) = &transcript.table else {
        return Ok(None);
    };
    Ok(genetic_code(*table, Reading::VepRs).map(|code| (transcript, code)))
}

/// The eleventh shape on a deletion: both descriptions are one deletion, of one residue or
/// of a range, with an extension's count after it (`theirs` the reference's captures,
/// `counted` the candidate's); the candidate's has a count and the reference's the unknown
/// count or another count; the reference row's terms hold `stop_lost`; and the row's coding
/// bases and its transcript in the candidate's cache give both counts
/// (`sequence_counts_a_deletion_from_the_lost_stop`).
fn deletion_counted_from_a_stop_the_peptide_holds(
    e: &Entry,
    theirs: &regex::Captures<'_>,
    counted: &regex::Captures<'_>,
    sequences: &DownstreamSequences,
) -> Result<bool> {
    if theirs[1] != counted[1]
        || &counted[4] == "?"
        || theirs[4] == counted[4]
        || !holds_term(&e.consequence, "stop_lost")
    {
        return Ok(false);
    }
    let number = |digits: &str| digits.parse::<usize>().ok();
    let Some(first_named) = number(&counted[2]) else {
        return Ok(false);
    };
    let last_named = counted
        .get(3)
        .map_or(Some(first_named), |last| number(last.as_str()));
    let named = last_named.and_then(|last| last.checked_add(1)?.checked_sub(first_named));
    let bounds = cds_bounds(e)
        .and_then(|(from, to)| Some((usize::try_from(from).ok()?, usize::try_from(to).ok()?)));
    let (Some(named), Some(bounds)) = (named, bounds) else {
        return Ok(false);
    };
    let Some((transcript, code)) = transcript_under_its_code(e, sequences)? else {
        return Ok(false);
    };
    Ok(sequence_counts_a_deletion_from_the_lost_stop(
        &transcript,
        code,
        bounds,
        (first_named, named),
        (&theirs[4], &counted[4]),
    ))
}

/// The eleventh shape of `hgvsp_stop_codon_protein_comparison`, which the row alone does
/// not decide: both values name one protein and are an extension at one residue number
/// with one first residue, a pair of prediction parentheses set aside; the candidate's has
/// a count and the reference's the unknown count or another count; the reference row's
/// terms hold `stop_lost`; and the row's transcript in the candidate's cache (`sequences`)
/// has a peptide that holds a stop and gives the candidate's count from the stop of that
/// residue number and the reference's count as Ensembl VEP's own arithmetic
/// (`sequence_counts_from_the_lost_stop`). Two values that are one deletion with an
/// extension's count after it are read the same way, beside the row's coding bases
/// (`deletion_counted_from_a_stop_the_peptide_holds`). The row's FLAGS are not read: such
/// a peptide names the stop it holds on a coding sequence that no stop codon closes as on
/// any other. The cache is read for a pair of the shape alone.
pub fn extension_counted_from_a_stop_the_peptide_holds(
    e: &Entry,
    ev: &str,
    av: &str,
    sequences: &DownstreamSequences,
) -> Result<bool> {
    if blank(ev) || blank(av) {
        return Ok(false);
    }
    let (reference_name, separator, reference) = pyfmt::partition(ev, ":p.");
    let (candidate_name, _, candidate) = pyfmt::partition(av, ":p.");
    if separator.is_empty() || reference_name != candidate_name {
        return Ok(false);
    }
    let predicted = |description: &str| description.starts_with('(') && description.ends_with(')');
    let (mut reference, mut candidate) = (reference, candidate);
    if predicted(reference) && predicted(candidate) {
        reference = &reference[1..reference.len() - 1];
        candidate = &candidate[1..candidate.len() - 1];
    }
    if let (Some(theirs), Some(counted)) = (
        DELETION_EXTENSION_WHOLE.captures(reference),
        DELETION_EXTENSION_WHOLE.captures(candidate),
    ) {
        return deletion_counted_from_a_stop_the_peptide_holds(e, &theirs, &counted, sequences);
    }
    let (Some(theirs), Some(counted)) = (
        STOP_EXTENSION_WHOLE.captures(reference),
        STOP_EXTENSION_WHOLE.captures(candidate),
    ) else {
        return Ok(false);
    };
    if theirs[1] != counted[1]
        || theirs[2] != counted[2]
        || &counted[2] == "Ter"
        || &counted[3] == "?"
        || theirs[3] == counted[3]
        || !holds_term(&e.consequence, "stop_lost")
    {
        return Ok(false);
    }
    let Ok(lost) = counted[1].parse::<usize>() else {
        return Ok(false);
    };
    let Some((transcript, code)) = transcript_under_its_code(e, sequences)? else {
        return Ok(false);
    };
    Ok(sequence_counts_from_the_lost_stop(
        &transcript,
        code,
        lost,
        (&theirs[3], &counted[3]),
    ))
}

/// The tenth shape of `hgvsp_stop_codon_protein_comparison`, on an HGVSp the reference
/// leaves empty: the candidate writes the stop as unchanged, alone or inside a pair of
/// prediction parentheses, at the residue after the one codon the reference row `e` names;
/// the row's terms hold `stop_retained_variant` and its FLAGS (`flags`) name no
/// `cds_end_NF`; its coding positions lie inside that codon; and its amino acids are that
/// codon's residue against the same residue and then a stop. By the row the alternate
/// protein is the reference protein up to that residue and then a stop, at the residue where
/// the candidate's value has the reference protein's stop.
pub fn stop_unchanged_where_the_reference_writes_none(
    e: &Entry,
    flags: &str,
    ev: &str,
    av: &str,
) -> bool {
    if !blank(ev) || blank(av) || flags.contains("cds_end_NF") {
        return false;
    }
    let (_, separator, candidate) = pyfmt::partition(av, ":p.");
    if separator.is_empty() || !holds_term(&e.consequence, STOP_RETAINED_TERM) {
        return false;
    }
    let candidate = candidate.replace("%3D", "=");
    let mut candidate = candidate.as_str();
    if candidate.starts_with('(') && candidate.ends_with(')') {
        candidate = &candidate[1..candidate.len() - 1];
    }
    let Some(unchanged) = STOP_UNCHANGED_WHOLE.captures(candidate) else {
        return false;
    };
    let (Some((first, last)), Some((cds_start, cds_end))) = (protein_span(e), cds_bounds(e)) else {
        return false;
    };
    let codon_end = 3 * i128::from(first);
    let inside_the_codon = codon_end - 2 <= i128::from(cds_start.min(cds_end))
        && i128::from(cds_start.max(cds_end)) <= codon_end;
    let amino_acids = field_of(Some(e), AMINO_ACIDS_FIELDS);
    let (residue, slash, alternate) = pyfmt::partition(&amino_acids, "/");
    let kept_and_stopped = slash == "/"
        && residue.len() == 1
        && residue
            .bytes()
            .all(|aa| aa.is_ascii_uppercase() && aa != b'X')
        && alternate
            .strip_prefix(residue)
            .is_some_and(|after| after.starts_with('*'));
    first == last
        && first >= 1
        && inside_the_codon
        && kept_and_stopped
        && whole_number(&unchanged[1]) == (i128::from(first) + 1).to_string()
}
