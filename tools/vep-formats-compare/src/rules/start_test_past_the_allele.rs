// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `hgvsp_start_codon_test_past_the_shifted_allele`: the protein notation of an
//! insertion or deletion that the 3' shift moves onto the start codon and that leaves the
//! whole coding sequence standing, which Ensembl VEP reads from the coding sequence edited
//! without the 5' UTR where its start codon test edits the sequence the shift length past
//! the allele. The shapes are read from the cells' texts alone;
//! `tools/vep-golden-corpus` compiles this file by path.

use regex::Regex;
use std::collections::BTreeSet;
use std::sync::LazyLock;

/// The term of a row whose start codon stands after the allele's edit.
const START_RETAINED_TERM: &str = "start_retained_variant";

/// The term of a row whose start codon the allele's edit takes.
const START_LOST_TERM: &str = "start_lost";

/// `pattern` held to the whole text.
fn whole(pattern: &str) -> Regex {
    Regex::new(&format!(r"\A(?:{pattern})\z")).expect("the pattern compiles")
}

/// A positive integer, the SHIFT_LENGTH of a row the 3' shift moved.
static POSITIVE_INTEGER: LazyLock<Regex> = LazyLock::new(|| whole(r"[1-9][0-9]*"));

/// The first residue unchanged: its name.
static FIRST_RESIDUE_UNCHANGED: LazyLock<Regex> = LazyLock::new(|| whole(r"([A-Z][a-z]{2})1="));

/// A frameshift at residue 1 with the unknown count: the first residue's name.
static FRAMESHIFT_AT_THE_FIRST_RESIDUE: LazyLock<Regex> =
    LazyLock::new(|| whole(r"([A-Z][a-z]{2})1[A-Z][a-z]{2}fsTer\?"));

/// An insertion of residues between a residue 0 and residue 1: the name of each flank.
static INSERTION_BEFORE_THE_FIRST_RESIDUE: LazyLock<Regex> =
    LazyLock::new(|| whole(r"([A-Z][a-z]{2})0_([A-Z][a-z]{2})1ins(?:[A-Z][a-z]{2})+"));

/// What the class reads of the two rows of one key beside the two HGVSp values, the
/// reference row's cell before the candidate row's in each pair.
#[derive(Clone, Copy, Debug)]
pub struct StartTestRows<'a> {
    /// Whether the rows are an insertion's or a deletion's.
    pub indel: bool,
    /// The consequence terms of each row, joined by `,` or `&`.
    pub terms: [&'a str; 2],
    /// The SHIFT_LENGTH of each row.
    pub shift_lengths: [&'a str; 2],
    /// Whether the CDS position of each row names a start, and whether it names an end.
    pub coding_bounds: [(bool, bool); 2],
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

/// `hgvsp_start_codon_test_past_the_shifted_allele` on a differing HGVSp, `reference`
/// Ensembl VEP's value and `candidate` vep-rs's. Under `--shift_3prime 1` `_inv_start_altered`
/// adds the shift length to a cDNA span the allele serves with the shift
/// (Utils/VariationEffect.pm:920, 928-930), so the second start codon test of `start_lost`
/// (:871) edits the sequence the shift length past the allele, and `_get_hgvs_peptides`,
/// which writes the lost start only under `start_lost` (TranscriptVariationAllele.pm:2161),
/// reads the notation from the coding sequence edited without the 5' UTR (:2394-2439). The
/// candidate tests the start codon at the moved allele.
///
/// The rows are an insertion's or a deletion's, carry one positive SHIFT_LENGTH and the
/// same terms, `start_retained_variant` among them and no `start_lost`, and a CDS position
/// with a start and an end, and the two values name one protein. The candidate's value is
/// the first residue unchanged (`Met1=`, the `=` as the text formats escape it too). Two
/// references are taken:
///
/// - a frameshift at residue 1 with the unknown count whose first residue is the
///   candidate's (`Met1TrpfsTer?`);
/// - an insertion between a residue 0 and residue 1 whose two flanks carry one name
///   (`Ala0_Ala1insSerVal`), which need not be the candidate's residue.
///
/// Every other pair stays a difference: an empty value, another residue number, a count in
/// the reference's frameshift, unequal flank names and the sides exchanged among them.
pub fn start_codon_test_past_the_shifted_allele(
    rows: &StartTestRows<'_>,
    reference: &str,
    candidate: &str,
) -> bool {
    let [reference_shift, candidate_shift] = rows.shift_lengths;
    if !rows.indel
        || reference_shift != candidate_shift
        || !POSITIVE_INTEGER.is_match(reference_shift)
        || rows.coding_bounds != [(true, true); 2]
    {
        return false;
    }
    let terms = term_set(rows.terms[0]);
    if terms != term_set(rows.terms[1])
        || !terms.contains(START_RETAINED_TERM)
        || terms.contains(START_LOST_TERM)
    {
        return false;
    }
    let Some((reference, candidate)) = descriptions_of_one_protein(reference, candidate) else {
        return false;
    };
    let candidate = candidate.replace("%3D", "=");
    let Some(unchanged) = FIRST_RESIDUE_UNCHANGED.captures(&candidate) else {
        return false;
    };
    if let Some(frameshift) = FRAMESHIFT_AT_THE_FIRST_RESIDUE.captures(reference) {
        return frameshift[1] == unchanged[1];
    }
    INSERTION_BEFORE_THE_FIRST_RESIDUE
        .captures(reference)
        .is_some_and(|inserted| inserted[1] == inserted[2])
}
