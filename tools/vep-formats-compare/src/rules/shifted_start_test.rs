// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `hgvsp_shifted_indel_written_span_start_test`: the protein notation of an
//! insertion or deletion that the 3' shift moves into the first codon, which Ensembl VEP
//! makes from the moved allele's residues and the start codon test of the allele as the
//! record writes it. The shapes are read from the cells' texts alone;
//! `tools/vep-golden-corpus` compiles this file by path.

use regex::Regex;
use std::collections::BTreeSet;
use std::sync::LazyLock;

/// The term of a row whose allele as written reaches the 5' UTR.
const UTR_TERM: &str = "5_prime_UTR_variant";

/// The term of a row whose allele as written holds coding bases and has no CDS start.
const CODING_TERM: &str = "coding_sequence_variant";

/// The terms of a start codon test: a row that holds one was tested where it is written.
const START_TERMS: [&str; 2] = ["start_lost", "start_retained_variant"];

/// `pattern` held to the whole text.
fn whole(pattern: &str) -> Regex {
    Regex::new(&format!(r"\A(?:{pattern})\z")).expect("the pattern compiles")
}

/// A frameshift at residue 1 with the unknown count: the first residue's name is group 1.
static FRAMESHIFT_AT_THE_FIRST_RESIDUE: LazyLock<Regex> =
    LazyLock::new(|| whole(r"([A-Z][a-z]{2})1[A-Z][a-z]{2}fsTer\?"));

/// A deletion of residues from the first: the first residue's name, the last residue's
/// name and the last residue's number.
static DELETED_FROM_THE_FIRST_RESIDUE: LazyLock<Regex> =
    LazyLock::new(|| whole(r"([A-Z][a-z]{2})1_([A-Z][a-z]{2})([1-9][0-9]*)del"));

/// The first residue unchanged: its name.
static FIRST_RESIDUE_UNCHANGED: LazyLock<Regex> = LazyLock::new(|| whole(r"([A-Z][a-z]{2})1="));

/// The lost start named at the first residue: its name, and the number of the last residue
/// of the span where the span ends past the first codon.
static LOST_START: LazyLock<Regex> =
    LazyLock::new(|| whole(r"([A-Z][a-z]{2})1(?:\?|_\?([1-9][0-9]*))"));

/// The lost start over a span of named residues: the names and the last residue's number.
static LOST_START_OVER_NAMED_RESIDUES: LazyLock<Regex> =
    LazyLock::new(|| whole(r"((?:[A-Z][a-z]{2})+)1_\?([1-9][0-9]*)"));

/// A `c.` deletion that ends on a coding base: that base's number.
static DELETION_ENDING_ON_A_CODING_BASE: LazyLock<Regex> =
    LazyLock::new(|| whole(r"c\.(?:\*?-?[0-9]+(?:[+-][0-9]+)?_)?([1-9][0-9]*)del"));

static NON_ZERO_INTEGER: LazyLock<Regex> = LazyLock::new(|| whole(r"-?[1-9][0-9]*"));

/// What the class reads of the two rows of one key beside the two HGVSp values, the
/// reference row's cell before the candidate row's in each pair.
#[derive(Clone, Copy, Debug)]
pub struct ShiftedRows<'a> {
    /// Whether the rows are an insertion's or a deletion's.
    pub indel: bool,
    /// The consequence terms of each row, joined by `,` or `&`.
    pub terms: [&'a str; 2],
    /// The HGVS_OFFSET of each row.
    pub offsets: [&'a str; 2],
    /// Whether the CDS position of each row names a start, and whether it names an end.
    pub coding_bounds: [(bool, bool); 2],
    /// The candidate row's HGVSc.
    pub candidate_hgvsc: &'a str,
}

/// The terms of a row as a set.
fn term_set(terms: &str) -> BTreeSet<&str> {
    terms
        .split([',', '&'])
        .filter(|term| !term.is_empty())
        .collect()
}

/// The descriptions of two protein notations of one protein identifier, the text after
/// `:p.`; a pair of prediction parentheses both carry is set aside.
fn descriptions_of_one_protein<'a>(
    reference: &'a str,
    candidate: &'a str,
) -> Option<(&'a str, &'a str)> {
    let (reference_name, reference) = reference.split_once(":p.")?;
    let (candidate_name, candidate) = candidate.split_once(":p.")?;
    if reference_name.is_empty() || reference_name != candidate_name {
        return None;
    }
    let predicted = |description: &str| {
        description.len() >= 2 && description.starts_with('(') && description.ends_with(')')
    };
    if predicted(reference) && predicted(candidate) {
        return Some((
            &reference[1..reference.len() - 1],
            &candidate[1..candidate.len() - 1],
        ));
    }
    Some((reference, candidate))
}

/// Whether the lost start `candidate` names ends on the residue its row's HGVSc gives: a
/// `c.` deletion that ends on coding base b ends on residue b / 3, rounded up, and a span
/// that ends inside the first codon is written without a last residue. Any other notation
/// says nothing of the residue, which is then 2 or more where the candidate names one.
fn lost_start_ends_where_the_deletion_ends(last_residue: Option<&str>, hgvsc: &str) -> bool {
    let body = hgvsc.split_once(':').map_or(hgvsc, |(_, body)| body);
    let Some(deletion) = DELETION_ENDING_ON_A_CODING_BASE.captures(body) else {
        return last_residue != Some("1");
    };
    let Ok(last_base) = deletion[1].parse::<u128>() else {
        return false;
    };
    match (last_base.div_ceil(3), last_residue) {
        (1, last_residue) => last_residue.is_none(),
        (residue, Some(last_residue)) => residue.to_string() == last_residue,
        (_, None) => false,
    }
}

/// `hgvsp_shifted_indel_written_span_start_test` on a differing HGVSp, `reference` Ensembl
/// VEP's value and `candidate` vep-rs's. `hgvs_protein` reads the frameshift, the deleted
/// residues and the peptides on the allele moved 3' (TranscriptVariationAllele.pm:1711,
/// 1724-1733, 1774-1793) and deletes the shift (:1796) before `_get_hgvs_peptides` asks
/// `start_lost` (:2161), whose test then reads the allele as the record writes it
/// (Utils/VariationEffect.pm:977-981). The candidate makes the start codon test on the
/// moved allele.
///
/// The rows are an insertion's or a deletion's, carry one non-zero HGVS_OFFSET and the same
/// terms, `5_prime_UTR_variant` among them and no term of a start codon test, and the two
/// values name one protein. Three pairs are taken:
///
/// - on rows with no CDS position, a frameshift at residue 1 with the unknown count against
///   the same residue unchanged (`Met1=`, the `=` as the text formats escape it too);
/// - on rows that hold `coding_sequence_variant` and whose CDS position has an end and no
///   start, the same frameshift against the lost start named at that residue (`Met1?`,
///   `Met1_?3`), ending on the residue of the last base of the candidate's `c.` deletion;
/// - on such rows, a deletion from residue 1 to residue n against the lost start over the
///   same n residues, each named, the first and the last as the reference names them.
///
/// Every other pair stays a difference: an empty value, another residue number, a count in
/// the reference's frameshift and the sides exchanged among them.
pub fn start_test_read_on_the_written_span(
    rows: &ShiftedRows<'_>,
    reference: &str,
    candidate: &str,
) -> bool {
    let [reference_offset, candidate_offset] = rows.offsets;
    if !rows.indel
        || reference_offset != candidate_offset
        || !NON_ZERO_INTEGER.is_match(reference_offset)
    {
        return false;
    }
    let terms = term_set(rows.terms[0]);
    if terms != term_set(rows.terms[1])
        || !terms.contains(UTR_TERM)
        || START_TERMS.iter().any(|term| terms.contains(term))
    {
        return false;
    }
    let Some((reference, candidate)) = descriptions_of_one_protein(reference, candidate) else {
        return false;
    };
    let candidate = candidate.replace("%3D", "=");
    let before_the_coding_start = rows.coding_bounds == [(false, false); 2];
    let across_the_first_base =
        terms.contains(CODING_TERM) && rows.coding_bounds == [(false, true); 2];
    if let Some(frameshift) = FRAMESHIFT_AT_THE_FIRST_RESIDUE.captures(reference) {
        if before_the_coding_start {
            return FIRST_RESIDUE_UNCHANGED
                .captures(&candidate)
                .is_some_and(|unchanged| unchanged[1] == frameshift[1]);
        }
        return across_the_first_base
            && LOST_START.captures(&candidate).is_some_and(|lost| {
                lost[1] == frameshift[1]
                    && lost_start_ends_where_the_deletion_ends(
                        lost.get(2).map(|last| last.as_str()),
                        rows.candidate_hgvsc,
                    )
            });
    }
    let Some(deleted) = DELETED_FROM_THE_FIRST_RESIDUE.captures(reference) else {
        return false;
    };
    let Some(lost) = LOST_START_OVER_NAMED_RESIDUES.captures(&candidate) else {
        return false;
    };
    let named: Vec<&str> = lost[1]
        .as_bytes()
        .chunks(3)
        .map(|name| std::str::from_utf8(name).unwrap_or_default())
        .collect();
    across_the_first_base
        && named.len() >= 2
        && named.len().to_string() == deleted[3]
        && lost[2] == deleted[3]
        && named.first().copied() == Some(&deleted[1])
        && named.last().copied() == Some(&deleted[2])
}
