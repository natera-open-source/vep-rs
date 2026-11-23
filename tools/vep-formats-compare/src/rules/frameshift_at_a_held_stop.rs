// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row class of an insertion or a deletion whose first base lies inside
//! the codon of a stop that the transcript's peptide holds as a stop and
//! reads on past: `frameshift_inside_a_held_stop_codon`.

use crate::consts::HGVSP_FIELDS;
use crate::entry::Entry;
use crate::error::Result;
use crate::fields::field_named;
use crate::rules::edited_model::CODONS_FIELDS;
use crate::rules::first_field;
use crate::rules::following::term_pair_following_cell;
use crate::rules::following::terms_write_field;
use crate::rules::settled_stop::RESIDUE_FIELDS;
use crate::rules::term_pairs::INFRAME_DELETION_TERM;
use crate::shared::blank;
use crate::shared::consequence_set;
use crate::shared::ABSENT;
use crate::shared::CODING_SEQUENCE_TERM;
use crate::shared::FRAMESHIFT_TERM;
use crate::shared::INFRAME_INSERTION_TERM;
use crate::shared::STOP_LOST_TERM;
use crate::shared::STOP_RETAINED_TERM;
use std::collections::BTreeSet;

/// The groups of terms of which the reference writes one whole where the
/// candidate writes `frameshift_variant`.
pub const REFERENCE_TERM_GROUPS: &[&[&str]] = &[
    &[CODING_SEQUENCE_TERM],
    &[INFRAME_DELETION_TERM],
    &[INFRAME_INSERTION_TERM, CODING_SEQUENCE_TERM],
];

/// The terms a reference row of the class does not hold. Membership only.
pub const REFUSED_REFERENCE_TERMS: &[&str] = &[STOP_RETAINED_TERM, STOP_LOST_TERM, FRAMESHIFT_TERM];

/// The bases one half of a `Codons` cell holds: none for a half written `-`.
fn bases(half: &str) -> usize {
    if half == "-" {
        0
    } else {
        half.chars().count()
    }
}

/// Whether the two halves of a `Codons` cell differ in length by a number of
/// bases that is no multiple of three. A cell without a `/` has no halves.
pub fn length_change_shifts_the_frame(codons: &str) -> bool {
    codons
        .split_once('/')
        .is_some_and(|(reference, alternate)| {
            !bases(reference)
                .abs_diff(bases(alternate))
                .is_multiple_of(3)
        })
}

/// Whether an `Amino_acids` cell opens with a stop and holds an `X` after its
/// `/`: the residues of a change that begins in a stop codon and ends inside
/// a codon.
pub fn residues_open_with_a_stop(residues: &str) -> bool {
    residues
        .split_once('/')
        .is_some_and(|(reference, alternate)| reference.starts_with('*') && alternate.contains('X'))
}

/// Whether the reference entry `e` and the candidate entry `a` of one key
/// have the shape of `frameshift_inside_a_held_stop_codon`. Ensembl VEP's
/// `frameshift` returns 0 wherever the reference peptide of the allele's
/// codons begins with a stop (Utils/VariationEffect.pm:1554), so the
/// reference writes no `frameshift_variant` on a row whose first base lies
/// inside the codon of a stop the peptide holds and reads on past.
///
/// The pair holds where all of these do: the candidate's terms hold
/// `frameshift_variant` and the reference's do not; the candidate's other
/// terms are the reference's without one group of `REFERENCE_TERM_GROUPS`,
/// which the reference holds whole; the reference holds none of
/// `REFUSED_REFERENCE_TERMS`; the `Amino_acids` cell is the same on both
/// sides, opens with `*` and holds `X` after its `/`; and the `Codons` cell
/// is the same on both sides, its halves differing in length by no multiple
/// of three. The rule reads the two entries alone.
pub fn frameshift_inside_a_held_stop_codon_pair(e: &Entry, a: &Entry) -> Result<bool> {
    let (expected_set, actual_set) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    let expected: BTreeSet<&str> = expected_set.split(',').collect();
    let actual: BTreeSet<&str> = actual_set.split(',').collect();
    if !actual.contains(FRAMESHIFT_TERM)
        || REFUSED_REFERENCE_TERMS
            .iter()
            .any(|term| expected.contains(term))
    {
        return Ok(false);
    }
    let beside_the_frameshift: Vec<&str> = actual
        .iter()
        .copied()
        .filter(|term| *term != FRAMESHIFT_TERM)
        .collect();
    let one_group_replaced = REFERENCE_TERM_GROUPS.iter().any(|group| {
        group.iter().all(|term| expected.contains(term))
            && expected
                .iter()
                .copied()
                .filter(|term| !group.contains(term))
                .eq(beside_the_frameshift.iter().copied())
    });
    if !one_group_replaced {
        return Ok(false);
    }
    let residues = first_field(e, RESIDUE_FIELDS);
    if residues != first_field(a, RESIDUE_FIELDS) || !residues_open_with_a_stop(residues) {
        return Ok(false);
    }
    let codons = first_field(e, CODONS_FIELDS);
    Ok(codons == first_field(a, CODONS_FIELDS) && length_change_shifts_the_frame(codons))
}

/// Whether the differing values `expected_value` and `actual_value` of the
/// field `name`, as the two entries carry them, are ones the class reads on
/// a row `frameshift_inside_a_held_stop_codon_pair` takes: an HGVSp of any
/// value on either side; a field of Downstream, NMD or the LoF plugin that
/// the candidate alone carries, where the terms let the plugin run on the
/// candidate's row and keep it from the reference's (`terms_write_field`);
/// and the IMPACT each side's terms give (`term_pair_following_cell`).
///
/// A plugin field written with an empty value is a field the candidate
/// carries, as Downstream writes one where the new frame ends at once; `-`
/// and an absent key are no field.
pub fn frameshift_inside_a_held_stop_codon_cell(
    e: &Entry,
    a: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    if expected_value == actual_value {
        return false;
    }
    if field_named(name, HGVSP_FIELDS) {
        return !(blank(expected_value) && blank(actual_value));
    }
    let (expected_set, actual_set) = (
        consequence_set(&e.consequence),
        consequence_set(&a.consequence),
    );
    let expected_terms: Vec<&str> = expected_set.split(',').collect();
    let actual_terms: Vec<&str> = actual_set.split(',').collect();
    if let Some(candidate_runs) = terms_write_field(name, &actual_terms) {
        let reference_runs = terms_write_field(name, &expected_terms) == Some(true);
        return candidate_runs
            && !reference_runs
            && blank(expected_value)
            && !matches!(actual_value, "-" | ABSENT);
    }
    let dashed = |value| if blank(value) { "-" } else { value };
    term_pair_following_cell(e, a, name, dashed(expected_value), dashed(actual_value))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An entry of the terms `terms` whose residues and codons stand under
    /// the names `fields` gives them.
    fn entry(terms: &str, fields: [(&str, &str); 2]) -> Entry {
        let mut entry = Entry {
            consequence: terms.to_string(),
            ..Entry::default()
        };
        for (name, value) in fields {
            entry.fields.insert(name.to_string(), value.to_string());
        }
        entry
    }

    /// Whether the pair rule takes the two rows of the terms `reference` and
    /// `candidate` with the cells `residues` and `codons` on both.
    fn taken(reference: &str, candidate: &str, residues: &str, codons: &str) -> bool {
        let cells = [("Amino_acids", residues), ("Codons", codons)];
        frameshift_inside_a_held_stop_codon_pair(&entry(reference, cells), &entry(candidate, cells))
            .expect("the pair is read")
    }

    #[test]
    fn a_length_change_of_no_whole_codons_shifts_the_frame() {
        for codons in [
            "tga/tTga",
            "tgA/tg",
            "tGAGCT/t",
            "tga/tgAATGa",
            "-/t",
            "tg/-",
        ] {
            assert!(length_change_shifts_the_frame(codons), "{codons}");
        }
        for codons in [
            "tga/tTTTga",
            "tga/tca",
            "tga/-",
            "-/TGA",
            "-/-",
            "tga",
            "-",
            "",
        ] {
            assert!(!length_change_shifts_the_frame(codons), "{codons}");
        }
    }

    #[test]
    fn the_residues_open_with_a_stop_and_end_inside_a_codon() {
        for residues in ["*/LX", "*/X", "*A/X", "*/*MX", "*/XL"] {
            assert!(residues_open_with_a_stop(residues), "{residues}");
        }
        for residues in ["A/LX", "A*/X", "*/L", "*X/L", "*X", "-/X", "-", ""] {
            assert!(!residues_open_with_a_stop(residues), "{residues}");
        }
    }

    #[test]
    fn the_reference_writes_one_group_where_the_candidate_writes_the_frameshift() {
        let shifted = "frameshift_variant";
        for (reference, candidate) in [
            ("coding_sequence_variant", shifted),
            ("inframe_deletion", shifted),
            ("inframe_insertion,coding_sequence_variant", shifted),
            // The VCF format joins a row's terms with `&`.
            ("inframe_insertion&coding_sequence_variant", shifted),
            (
                "coding_sequence_variant,NMD_transcript_variant",
                "frameshift_variant,NMD_transcript_variant",
            ),
            (
                "splice_region_variant&inframe_deletion",
                "frameshift_variant&splice_region_variant",
            ),
        ] {
            assert!(
                taken(reference, candidate, "*/LX", "tga/tTga"),
                "{reference}"
            );
            assert!(
                !taken(candidate, reference, "*/LX", "tga/tTga"),
                "{reference}"
            );
        }
        for (reference, candidate) in [
            ("protein_altering_variant", shifted),
            ("inframe_insertion", shifted),
            ("inframe_deletion,coding_sequence_variant", shifted),
            (
                "coding_sequence_variant",
                "frameshift_variant,coding_sequence_variant",
            ),
            (
                "coding_sequence_variant",
                "frameshift_variant,splice_region_variant",
            ),
            ("coding_sequence_variant,NMD_transcript_variant", shifted),
            ("stop_retained_variant", shifted),
            (
                "stop_retained_variant,coding_sequence_variant",
                "frameshift_variant,stop_retained_variant",
            ),
            (
                "stop_lost,coding_sequence_variant",
                "stop_lost,frameshift_variant",
            ),
            ("frameshift_variant,coding_sequence_variant", shifted),
            ("coding_sequence_variant", "coding_sequence_variant"),
            ("", shifted),
        ] {
            assert!(
                !taken(reference, candidate, "*/LX", "tga/tTga"),
                "{reference}"
            );
        }
    }

    #[test]
    fn the_two_rows_carry_the_same_residues_and_codons() {
        let (reference, candidate) = ("coding_sequence_variant", "frameshift_variant");
        assert!(taken(reference, candidate, "*/X", "tgA/tg"));
        assert!(!taken(reference, candidate, "A/X", "tgA/tg"));
        assert!(!taken(reference, candidate, "*/X", "tga/tca"));
        let e = entry(reference, [("Amino_acids", "*/X"), ("Codons", "tgA/tg")]);
        for cells in [
            [("Amino_acids", "*/LX"), ("Codons", "tgA/tg")],
            [("Amino_acids", "*/X"), ("Codons", "tGa/ta")],
        ] {
            let a = entry(candidate, cells);
            assert!(!frameshift_inside_a_held_stop_codon_pair(&e, &a).expect("the pair is read"));
        }
        // A JSON entry carries the two cells under lower-case keys.
        let cells = [("amino_acids", "*/X"), ("codons", "tgA/tg")];
        assert!(frameshift_inside_a_held_stop_codon_pair(
            &entry(reference, cells),
            &entry(candidate, cells)
        )
        .expect("the pair is read"));
        // A row without either cell is no row of the class.
        let bare = [("IMPACT", "HIGH"), ("STRAND", "1")];
        assert!(!frameshift_inside_a_held_stop_codon_pair(
            &entry(reference, bare),
            &entry(candidate, bare)
        )
        .expect("the pair is read"));
    }

    #[test]
    fn a_plugin_field_written_empty_is_carried_and_a_dash_is_not() {
        let cells = [("Amino_acids", "*/X"), ("Codons", "tGa/ta")];
        let e = entry("coding_sequence_variant", cells);
        let a = entry("frameshift_variant", cells);
        let cell = |name: &str, expected_value: &str, actual_value: &str| {
            frameshift_inside_a_held_stop_codon_cell(&e, &a, name, expected_value, actual_value)
        };
        assert!(cell("DownstreamProtein", ABSENT, ""));
        assert!(cell("DownstreamProtein", "-", ""));
        assert!(cell("DownstreamProtein", "", "LLP"));
        assert!(!cell("DownstreamProtein", "", ""));
        assert!(!cell("DownstreamProtein", ABSENT, "-"));
        assert!(!cell("DownstreamProtein", "LLP", ""));
        assert!(cell("NMD", ABSENT, "NMD_escaping_variant"));
        assert!(!cell("NMD", "NMD_escaping_variant", ABSENT));
        assert!(cell("HGVSp", "p.Ter10=", "p.Ala11LeufsTer12"));
        assert!(cell("HGVSp", "-", "p.Ala11LeufsTer12"));
        assert!(!cell("HGVSp", "-", ABSENT));
        assert!(cell("IMPACT", "MODIFIER", "HIGH"));
        assert!(!cell("IMPACT", "MODIFIER", "MODERATE"));
        assert!(!cell("STRAND", "1", "-1"));
    }
}
