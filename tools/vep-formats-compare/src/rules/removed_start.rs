// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row class of a deletion that removes whole codons from the first coding base of a
//! transcript with no 5' UTR: `start_codon_removed_without_five_prime_utr`.

use crate::consts::HGVSP_FIELDS;
use crate::entry::Entry;
use crate::error::Result;
use crate::fields::field_named;
use crate::pyfmt;
use crate::rules::first_field;
use crate::rules::following::term_pair_following_cell;
use crate::rules::settled_stop::RESIDUE_FIELDS;
use crate::rules::term_pairs::transcript_row;
use crate::rules::term_pairs::INFRAME_DELETION_TERM;
use crate::rules::whole_text;
use crate::shared::blank;
use crate::shared::consequence_set;
use crate::shared::transcript_positions;
use std::collections::BTreeSet;
use std::sync::LazyLock;

pub const START_LOST_TERM: &str = "start_lost";

pub const START_RETAINED_TERM: &str = "start_retained_variant";

/// A deletion from residue 1: the first residue, and past one residue the last with its
/// number.
pub const DELETION_FROM_RESIDUE_1_RE_PATTERN: &str =
    r"([A-Z][a-z]{2})1(?:_([A-Z][a-z]{2})([1-9][0-9]*))?del";
static DELETION_FROM_RESIDUE_1_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(DELETION_FROM_RESIDUE_1_RE_PATTERN));

/// A lost start: every residue removed, and past one residue the last residue's number.
pub const LOST_START_RE_PATTERN: &str = r"((?:[A-Z][a-z]{2})+)1(?:\?|_\?([1-9][0-9]*))";
static LOST_START_WHOLE: LazyLock<regex::Regex> =
    LazyLock::new(|| whole_text(LOST_START_RE_PATTERN));

/// The protein a notation names, its change and whether the change stands in parentheses;
/// `None` for a text without `:p.` and for a change with one parenthesis alone.
fn protein_change(notation: &str) -> Option<(&str, &str, bool)> {
    let (protein, separator, change) = pyfmt::partition(notation, ":p.");
    if separator.is_empty() || protein.is_empty() {
        return None;
    }
    match (change.strip_prefix('('), change.strip_suffix(')')) {
        (Some(_), Some(_)) => Some((protein, &change[1..change.len() - 1], true)),
        (None, None) => Some((protein, change, false)),
        _ => None,
    }
}

/// Whether `reference` is the protein notation of a deletion from residue 1 and
/// `candidate` the notation of a lost start over the same residues of the same protein.
/// A deletion is written with its first residue and, past one residue, its last
/// (`p.Met1del`, `p.Met1_Ala2del`; TranscriptVariationAllele.pm:2018-2027). A row that
/// carries `start_lost` is written with every residue the allele removes and `?` for the
/// alternate (:2161-2165): `p.Met1?` for one residue, `p.MetAla1_?2` past one
/// (:2030-2037). Both notations stand in parentheses or neither does (:1908, :2040).
pub fn lost_start_notations(reference: &str, candidate: &str) -> bool {
    let (Some((protein, deleted, predicted)), Some((named, lost, as_predicted))) =
        (protein_change(reference), protein_change(candidate))
    else {
        return false;
    };
    if protein != named || predicted != as_predicted {
        return false;
    }
    let (Some(deleted), Some(lost)) = (
        DELETION_FROM_RESIDUE_1_WHOLE.captures(deleted),
        LOST_START_WHOLE.captures(lost),
    ) else {
        return false;
    };
    let residues = &lost[1];
    if !residues.starts_with(&deleted[1]) {
        return false;
    }
    match (deleted.get(2), deleted.get(3), lost.get(2)) {
        // One residue on both sides.
        (None, None, None) => residues.len() == 3,
        // A range that ends on one residue and one number on both sides, the candidate
        // writing as many residues as the number counts.
        (Some(last), Some(end), Some(lost_end)) => {
            end.as_str() == lost_end.as_str()
                && end.as_str() != "1"
                && residues.ends_with(last.as_str())
                && end
                    .as_str()
                    .parse::<usize>()
                    .is_ok_and(|count| count.checked_mul(3) == Some(residues.len()))
        }
        _ => false,
    }
}

/// The cDNA end and the CDS end of an entry whose cDNA position and CDS position both
/// begin at 1 and end at a number: the coding sequence opens on the transcript's first
/// base and the span on the first coding base. `None` for any other entry.
fn ends_of_a_span_from_the_first_base(e: &Entry) -> Option<(String, String)> {
    let (cdna_start, cdna_end, cds_start, cds_end) = transcript_positions(e)?;
    let number = |text: &str| !text.is_empty() && text.bytes().all(|digit| digit.is_ascii_digit());
    (cdna_start == "1" && cds_start == "1" && number(&cdna_end) && number(&cds_end))
        .then_some((cdna_end, cds_end))
}

/// Whether the reference entry `e` and the candidate entry `a` of one key have the shape
/// of `start_codon_removed_without_five_prime_utr`: a Transcript entry; the candidate's
/// terms the reference's and `start_lost`, both sets holding `inframe_deletion` and
/// neither `start_retained_variant`; a cDNA position and a CDS position that begin at 1
/// and end at a number, the same on both sides; and residues, the same on both sides,
/// whose alternate half is `-`.
pub fn start_codon_removed_without_five_prime_utr_pair(e: &Entry, a: &Entry) -> Result<bool> {
    // The sets are built for the rows whose two Consequence texts differ alone.
    if !transcript_row(e) || e.consequence == a.consequence {
        return Ok(false);
    }
    let (expected_set, actual_set) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    let expected: BTreeSet<&str> = expected_set.split(',').collect();
    let actual: BTreeSet<&str> = actual_set.split(',').collect();
    if !expected.is_subset(&actual)
        || !actual.difference(&expected).eq([&START_LOST_TERM])
        || !expected.contains(INFRAME_DELETION_TERM)
        || expected.contains(START_RETAINED_TERM)
    {
        return Ok(false);
    }
    let ends = ends_of_a_span_from_the_first_base(e);
    if ends.is_none() || ends != ends_of_a_span_from_the_first_base(a) {
        return Ok(false);
    }
    let residues = first_field(e, RESIDUE_FIELDS);
    if residues != first_field(a, RESIDUE_FIELDS) {
        return Ok(false);
    }
    Ok(matches!(
        residues.split_once('/'),
        Some((removed, "-")) if !blank(removed)
    ))
}

/// Whether the differing values `ev` and `av` of the field `name` on a row
/// `start_codon_removed_without_five_prime_utr_pair` takes are the ones the class reads:
/// a protein notation `lost_start_notations` takes, and the cells that follow the two sets
/// of terms (`term_pair_following_cell`).
pub fn start_codon_removed_without_five_prime_utr_cell(
    e: &Entry,
    a: &Entry,
    name: &str,
    ev: &str,
    av: &str,
) -> bool {
    if field_named(name, HGVSP_FIELDS) {
        return lost_start_notations(ev, av);
    }
    term_pair_following_cell(e, a, name, ev, av)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A Transcript entry of a deletion in the spelling of the text formats: its terms,
    /// its two position columns and its residues.
    fn row(terms: &str, cdna: &str, cds: &str, residues: &str) -> Entry {
        let mut entry = Entry {
            key: "21:100-102|-|ENST1|Transcript".to_string(),
            allele: "-".to_string(),
            consequence: terms.to_string(),
            ..Entry::default()
        };
        for (name, value) in [
            ("cDNA_position", cdna),
            ("CDS_position", cds),
            ("Amino_acids", residues),
        ] {
            entry.fields.insert(name.to_string(), value.to_string());
        }
        entry
    }

    /// The same entry in the spelling of the JSON format: a start key and an end key a
    /// position, left out where `None`.
    fn object(terms: &str, cdna: [Option<&str>; 2], cds: [Option<&str>; 2]) -> Entry {
        let mut entry = Entry {
            key: "0|-|ENST1|transcript".to_string(),
            allele: "-".to_string(),
            consequence: terms.to_string(),
            ..Entry::default()
        };
        let keys = ["cdna_start", "cdna_end", "cds_start", "cds_end"];
        for (name, value) in keys.into_iter().zip(cdna.into_iter().chain(cds)) {
            if let Some(value) = value {
                entry.fields.insert(name.to_string(), value.to_string());
            }
        }
        entry
            .fields
            .insert("amino_acids".to_string(), "M/-".to_string());
        entry
    }

    fn pair(e: &Entry, a: &Entry) -> bool {
        start_codon_removed_without_five_prime_utr_pair(e, a).expect("the pair is read")
    }

    const KEPT: &str = "inframe_deletion";
    const LOST: &str = "start_lost,inframe_deletion";

    #[test]
    fn the_pair_is_a_whole_codon_removal_from_the_first_base_that_gains_start_lost() {
        for (cdna, cds, residues) in [
            ("1-3", "1-3", "M/-"),
            ("1-6", "1-6", "MA/-"),
            ("1-705", "1-705", "MAGS/-"),
            // Under --total_length each position carries the sequence's length.
            ("1-3/1500", "1-3/900", "M/-"),
        ] {
            let e = row(KEPT, cdna, cds, residues);
            let a = row(LOST, cdna, cds, residues);
            assert!(pair(&e, &a), "{cdna}");
            assert!(!pair(&a, &e), "{cdna}");
            assert!(!pair(&e, &e), "{cdna}");
        }
        // The VCF format joins a row's terms with `&`, and a term both sides write stays.
        let e = row(
            "splice_donor_variant&inframe_deletion&intron_variant",
            "1-12",
            "1-12",
            "MAGS/-",
        );
        let a = row(
            "splice_donor_variant&start_lost&inframe_deletion&intron_variant",
            "1-12",
            "1-12",
            "MAGS/-",
        );
        assert!(pair(&e, &a));
        let whole = [Some("1"), Some("3")];
        assert!(pair(
            &object(KEPT, whole, whole),
            &object(LOST, whole, whole)
        ));
    }

    #[test]
    fn a_row_that_breaks_one_condition_is_refused() {
        let e = row(KEPT, "1-3", "1-3", "M/-");
        let a = row(LOST, "1-3", "1-3", "M/-");
        assert!(pair(&e, &a));
        // A feature that is no transcript.
        let regulatory = |entry: &Entry| Entry {
            key: entry.key.replace("Transcript", "RegulatoryFeature"),
            ..entry.clone()
        };
        assert!(!pair(&regulatory(&e), &regulatory(&a)));
        // Terms that differ otherwise than by the candidate's start_lost.
        for (reference, candidate) in [
            (KEPT, "start_lost,inframe_deletion,splice_region_variant"),
            ("inframe_deletion,splice_region_variant", LOST),
            (KEPT, "stop_lost,inframe_deletion"),
            (KEPT, "start_lost"),
            (LOST, LOST),
            (
                "inframe_deletion,start_retained_variant",
                "start_lost,inframe_deletion,start_retained_variant",
            ),
        ] {
            assert!(
                !pair(
                    &row(reference, "1-3", "1-3", "M/-"),
                    &row(candidate, "1-3", "1-3", "M/-")
                ),
                "{reference} against {candidate}"
            );
        }
        // No inframe_deletion.
        for kept in ["coding_sequence_variant", "frameshift_variant", ""] {
            let lost = format!("start_lost,{kept}");
            assert!(
                !pair(
                    &row(kept, "1-3", "1-3", "M/-"),
                    &row(&lost, "1-3", "1-3", "M/-")
                ),
                "{kept}"
            );
        }
        // A position that does not begin at 1 and end at a number, or differs between the
        // sides.
        for (cdna, cds) in [
            ("2-4", "1-3"),
            ("1-3", "4-6"),
            ("?-3", "1-3"),
            ("1-?", "1-3"),
            ("1-3", "?-3"),
            ("1-3", "1-?"),
            ("-", "1-3"),
            ("1-3", "-"),
            ("11-13", "1-3"),
            ("01-3", "1-3"),
        ] {
            assert!(
                !pair(&row(KEPT, cdna, cds, "M/-"), &row(LOST, cdna, cds, "M/-")),
                "{cdna} {cds}"
            );
        }
        assert!(!pair(&e, &row(LOST, "1-6", "1-3", "M/-")));
        assert!(!pair(&e, &row(LOST, "1-3", "1-6", "M/-")));
        let whole = [Some("1"), Some("3")];
        for (cdna, cds) in [
            ([Some("2"), Some("4")], whole),
            (whole, [Some("4"), Some("6")]),
            ([None, Some("3")], whole),
            (whole, [None, Some("3")]),
            ([Some("1"), None], whole),
        ] {
            assert!(
                !pair(&object(KEPT, cdna, cds), &object(LOST, cdna, cds)),
                "{cdna:?} {cds:?}"
            );
        }
        assert!(!pair(
            &object(KEPT, whole, whole),
            &object(LOST, [Some("1"), Some("6")], whole)
        ));
        // Residues with an alternate half, without one, or that differ between the sides.
        for residues in ["MA/A", "M/L", "M", "-", "", "/-", "-/-", "M/-/-"] {
            assert!(
                !pair(
                    &row(KEPT, "1-3", "1-3", residues),
                    &row(LOST, "1-3", "1-3", residues)
                ),
                "{residues}"
            );
        }
        assert!(!pair(&e, &row(LOST, "1-3", "1-3", "L/-")));
    }

    #[test]
    fn a_deletion_from_residue_1_stands_against_the_lost_start_of_the_same_residues() {
        let taken = |reference: &str, candidate: &str| {
            lost_start_notations(
                &format!("ENSP1.1:p.{reference}"),
                &format!("ENSP1.1:p.{candidate}"),
            )
        };
        for (reference, candidate) in [
            ("Met1del", "Met1?"),
            ("Leu1del", "Leu1?"),
            ("Met1_Ala2del", "MetAla1_?2"),
            ("Met1_Ser4del", "MetAlaGlySer1_?4"),
            ("Met1_Ter3del", "MetAlaTer1_?3"),
            ("(Met1del)", "(Met1?)"),
            ("(Met1_Ala2del)", "(MetAla1_?2)"),
        ] {
            assert!(taken(reference, candidate), "{reference} {candidate}");
            assert!(!taken(candidate, reference), "{reference} {candidate}");
        }
        for (reference, candidate) in [
            // Another first residue, another last residue, another count of residues.
            ("Met1del", "Leu1?"),
            ("Met1_Ala2del", "MetGly1_?2"),
            ("Met1_Ala2del", "LeuAla1_?2"),
            ("Met1del", "MetAla1_?2"),
            ("Met1_Ala2del", "Met1?"),
            ("Met1_Ala2del", "MetAla1_?3"),
            ("Met1_Ala3del", "MetAla1_?3"),
            ("Met1_Ala2del", "MetAlaGly1_?2"),
            ("Met1_Met1del", "Met1?"),
            ("Met1del", "Met1_?1"),
            // A deletion that begins past residue 1, and another residue numbered first.
            ("Ala2del", "Met1?"),
            ("Met2_Ala3del", "MetAla1_?2"),
            ("Met1del", "Met2?"),
            ("Met11del", "Met11?"),
            // A notation of another kind on either side.
            ("Met1del", "Met1del"),
            ("Met1?", "Met1?"),
            ("Met1del", "Met1Leu"),
            ("Met1del", "Met1="),
            ("Met1dup", "Met1?"),
            ("Met1_Ala2delinsLeu", "MetAla1_?2"),
            ("Met1AlafsTer4", "Met1?"),
            ("Met1del", "?"),
            ("Met1del", "1?"),
            ("M1del", "M1?"),
            ("Met1del", "MetA1?"),
            // Parentheses on one side alone, and unpaired.
            ("(Met1del)", "Met1?"),
            ("Met1del", "(Met1?)"),
            ("(Met1del", "(Met1?"),
            ("Met1del", "Met1? "),
        ] {
            assert!(!taken(reference, candidate), "{reference} {candidate}");
        }
        // Two proteins, a notation without a protein, and no notation.
        assert!(!lost_start_notations(
            "ENSP1.1:p.Met1del",
            "ENSP2.1:p.Met1?"
        ));
        assert!(!lost_start_notations("p.Met1del", "p.Met1?"));
        assert!(!lost_start_notations("Met1del", "Met1?"));
        assert!(!lost_start_notations("-", "ENSP1.1:p.Met1?"));
        assert!(!lost_start_notations("ENSP1.1:p.Met1del", "-"));
        assert!(!lost_start_notations("", ""));
    }

    #[test]
    fn the_class_reads_the_notation_and_the_impact_of_each_sides_terms() {
        let e = row(KEPT, "1-3", "1-3", "M/-");
        let a = row(LOST, "1-3", "1-3", "M/-");
        let cell = |name: &str, ev: &str, av: &str| {
            start_codon_removed_without_five_prime_utr_cell(&e, &a, name, ev, av)
        };
        assert!(cell("HGVSp", "ENSP1.1:p.Met1del", "ENSP1.1:p.Met1?"));
        assert!(cell("hgvsp", "ENSP1.1:p.Met1del", "ENSP1.1:p.Met1?"));
        assert!(!cell("HGVSp", "ENSP1.1:p.Met1?", "ENSP1.1:p.Met1del"));
        assert!(!cell("HGVSp", "-", "ENSP1.1:p.Met1?"));
        // Utils/Constants.pm gives inframe_deletion the impact MODERATE and start_lost HIGH.
        assert!(cell("IMPACT", "MODERATE", "HIGH"));
        assert!(cell("impact", "MODERATE", "HIGH"));
        assert!(!cell("IMPACT", "HIGH", "MODERATE"));
        assert!(!cell("IMPACT", "LOW", "HIGH"));
        assert!(!cell("IMPACT", "MODERATE", "MODIFIER"));
        assert!(!cell("Codons", "ATG/-", "ATC/-"));
        assert!(!cell("STRAND", "1", "-1"));
    }
}
