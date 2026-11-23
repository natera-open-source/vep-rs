// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `alternate_codon_past_the_coding_end`: the alternate codons Ensembl VEP cuts
//! from a coding sequence with the 3' UTR joined on, and the residues read from them. The
//! shapes are read from the cells' texts alone; `tools/vep-golden-corpus` compiles this file
//! by path.

use crate::genetic_codes::{genetic_code, Reading, STANDARD_CODE};

/// The bases a codon string may hold of the 3' UTR: unchanged, so lower-case.
const UNCHANGED_BASES: &[u8] = b"acgtn";

/// The two halves of a Codons value, `None` for a value without its `/`.
fn halves(codons: &str) -> Option<(&str, &str)> {
    codons.split_once('/')
}

/// Whether every letter of `text` is one of `letters`.
fn spelt_in(text: &str, letters: &[u8]) -> bool {
    text.bytes().all(|letter| letters.contains(&letter))
}

/// What `peptide` of TranscriptVariationAllele.pm (lines 724-778) writes for a codon string
/// under the standard genetic code: the residues of its whole codons, `X` behind them where
/// bases are left over and the residues are not a lone stop, and `-` for no codon. `None`
/// for a string that holds another letter than the four bases.
fn peptide(codons: &str) -> Option<String> {
    if codons == "-" {
        return Some("-".to_string());
    }
    let bases = codons.to_ascii_uppercase();
    if bases.is_empty() || !spelt_in(&bases, b"ACGT") {
        return None;
    }
    let code = genetic_code(STANDARD_CODE, Reading::EnsemblVep)?;
    let whole = bases.len() / 3 * 3;
    let mut residues = String::new();
    for codon in bases.as_bytes()[..whole].chunks(3) {
        residues.push(char::from(code.amino_acid(codon)?));
    }
    if whole < bases.len() && residues != "*" {
        residues.push('X');
    }
    Some(residues)
}

/// The first shape: a codon window that runs past a coding sequence's end. The two
/// reference halves are equal and no whole number of codons, since the window is cut
/// where the coding sequence ends; the reference's alternate half is the candidate's with
/// one or two unchanged bases behind it, at most as many as the last codon lacks. A
/// candidate half of `-` holds no coding base; an empty half is refused.
fn window_past_the_coding_end(expected: &str, actual: &str) -> bool {
    let (Some((expected_reference, expected_alternate)), Some((reference, alternate))) =
        (halves(expected), halves(actual))
    else {
        return false;
    };
    if expected_reference != reference || !spelt_in(reference, b"ACGTacgt") {
        return false;
    }
    let lacking = match reference.len() % 3 {
        0 => return false,
        over => 3 - over,
    };
    let coding = if alternate == "-" { "" } else { alternate };
    if alternate.is_empty() || !spelt_in(coding, b"ACGTacgt") {
        return false;
    }
    expected_alternate
        .strip_prefix(coding)
        .is_some_and(|read| (1..=lacking).contains(&read.len()) && spelt_in(read, UNCHANGED_BASES))
}

/// The second shape: a deletion that leaves one or two bases of the coding sequence, on a
/// row whose protein positions open at `first_residue`. The two reference halves are
/// equal, whole codons from residue 1, and their lower-case bases are what the deletion
/// leaves; the candidate's alternate half is those bases, and the reference's is as many
/// unchanged bases that spell something else.
fn window_of_a_coding_sequence_under_one_codon(
    first_residue: &str,
    expected: &str,
    actual: &str,
) -> bool {
    let (Some((expected_reference, expected_alternate)), Some((reference, alternate))) =
        (halves(expected), halves(actual))
    else {
        return false;
    };
    if expected_reference != reference
        || reference.is_empty()
        || reference.len() % 3 != 0
        || !spelt_in(reference, b"ACGTacgt")
        || first_residue != "1"
    {
        return false;
    }
    let left: String = reference.chars().filter(char::is_ascii_lowercase).collect();
    (1..=2).contains(&left.len())
        && alternate == left
        && expected_alternate.len() == left.len()
        && expected_alternate != left
        && spelt_in(expected_alternate, UNCHANGED_BASES)
}

/// The alternate half of an Amino_acids value, which is its one half where the two
/// peptides are equal.
fn alternate_residues(amino_acids: &str) -> &str {
    amino_acids
        .split_once('/')
        .map_or(amino_acids, |(_, alternate)| alternate)
}

/// The reference half of an Amino_acids value.
fn reference_residues(amino_acids: &str) -> &str {
    amino_acids
        .split_once('/')
        .map_or(amino_acids, |(reference, _)| reference)
}

/// Whether the Codons value of the reference (`expected`) and of the candidate (`actual`)
/// have one of the class's two shapes, on a row whose protein positions open at
/// `first_residue`.
pub fn codons_read_past_the_coding_end(first_residue: &str, expected: &str, actual: &str) -> bool {
    window_past_the_coding_end(expected, actual)
        || window_of_a_coding_sequence_under_one_codon(first_residue, expected, actual)
}

/// Whether the Amino_acids value of the reference (`expected`) and of the candidate
/// (`actual`) have the class's shape, on a row with the two Codons values given: the
/// codons have the first shape, the two values open with the same reference residues, and
/// each alternate half is the translation of that side's alternate codons. A value of two
/// equal halves is refused: Ensembl VEP writes one half where the two peptides are equal
/// (TranscriptVariationAllele.pm:629), and so does the candidate.
pub fn residues_read_past_the_coding_end(
    expected_codons: &str,
    actual_codons: &str,
    expected: &str,
    actual: &str,
) -> bool {
    let two_equal_halves = |amino_acids: &str| {
        amino_acids
            .split_once('/')
            .is_some_and(|(reference, alternate)| reference == alternate)
    };
    if two_equal_halves(expected) || two_equal_halves(actual) {
        return false;
    }
    if !window_past_the_coding_end(expected_codons, actual_codons) {
        return false;
    }
    let (Some((_, expected_alternate)), Some((_, alternate))) =
        (halves(expected_codons), halves(actual_codons))
    else {
        return false;
    };
    reference_residues(expected) == reference_residues(actual)
        && peptide(expected_alternate).as_deref() == Some(alternate_residues(expected))
        && peptide(alternate).as_deref() == Some(alternate_residues(actual))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_codon_string_translates_as_ensembl_vep_writes_its_peptide() {
        assert_eq!(peptide("Aag").as_deref(), Some("K"));
        assert_eq!(peptide("Aa").as_deref(), Some("X"));
        assert_eq!(peptide("tGa").as_deref(), Some("*"));
        assert_eq!(peptide("tgaa").as_deref(), Some("*"));
        assert_eq!(peptide("aagtg").as_deref(), Some("KX"));
        assert_eq!(peptide("-").as_deref(), Some("-"));
        assert_eq!(peptide("aNg"), None);
        assert_eq!(peptide(""), None);
    }

    #[test]
    fn a_window_past_the_coding_end_lacks_at_most_what_the_last_codon_lacks() {
        assert!(window_past_the_coding_end("Ta/Aag", "Ta/Aa"));
        assert!(window_past_the_coding_end("T/Cag", "T/C"));
        assert!(window_past_the_coding_end("T/Ca", "T/C"));
        assert!(window_past_the_coding_end("tA/tg", "tA/t"));
        assert!(window_past_the_coding_end("T/ag", "T/-"));
        // Two bases behind a last codon that lacks one, an upper-case base, whole
        // reference codons, another reference half and a candidate that is no prefix.
        assert!(!window_past_the_coding_end("Ta/Aagg", "Ta/Aa"));
        assert!(!window_past_the_coding_end("Ta/AaG", "Ta/Aa"));
        assert!(!window_past_the_coding_end("Tag/Aagc", "Tag/Aag"));
        assert!(!window_past_the_coding_end("Ta/Aag", "Tc/Aa"));
        assert!(!window_past_the_coding_end("Ta/Aag", "Ta/Ca"));
        assert!(!window_past_the_coding_end("Ta/Aa", "Ta/Aa"));
        assert!(!window_past_the_coding_end("Ta", "Ta"));
        assert!(!window_past_the_coding_end("T/ag", "T/"));
    }

    #[test]
    fn a_deletion_that_leaves_one_or_two_coding_bases_is_read_at_residue_one() {
        assert!(codons_read_past_the_coding_end("1", "ATGTAa/g", "ATGTAa/a"));
        assert!(codons_read_past_the_coding_end(
            "1",
            "aTGTAa/gg",
            "aTGTAa/aa"
        ));
        assert!(!codons_read_past_the_coding_end(
            "2", "ATGTAa/g", "ATGTAa/a"
        ));
        assert!(!codons_read_past_the_coding_end(
            "1", "ATGTAa/a", "ATGTAa/a"
        ));
        assert!(!codons_read_past_the_coding_end(
            "1", "ATGTAa/G", "ATGTAa/a"
        ));
        assert!(!codons_read_past_the_coding_end(
            "1",
            "ATGtaa/g",
            "ATGtaa/taa"
        ));
    }

    #[test]
    fn residues_follow_each_side_of_a_window_past_the_coding_end() {
        assert!(residues_read_past_the_coding_end(
            "Ta/Aag", "Ta/Aa", "X/K", "X"
        ));
        assert!(residues_read_past_the_coding_end("T/ag", "T/-", "X", "X/-"));
        assert!(!residues_read_past_the_coding_end(
            "Ta/Aag", "Ta/Aa", "X/R", "X"
        ));
        assert!(!residues_read_past_the_coding_end(
            "Ta/Aag", "Ta/Aa", "X/K", "X/N"
        ));
        assert!(!residues_read_past_the_coding_end(
            "Ta/Aag", "Ta/Aa", "Y/K", "X"
        ));
        assert!(!residues_read_past_the_coding_end(
            "Ta/Aa", "Ta/Aa", "X/K", "X"
        ));
        assert!(!residues_read_past_the_coding_end(
            "Ta/Aag", "Ta/Aa", "X/K", "X/X"
        ));
    }
}
