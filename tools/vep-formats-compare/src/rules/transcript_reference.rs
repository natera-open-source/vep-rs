// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reference allele of a one-base substitution on a RefSeq model whose alignment
//! correction failed, and the cells and notations that follow it.

use crate::compare::CONSEQUENCE_FIELDS;
use crate::consts::BAM_EDIT_FIELDS;
use crate::consts::HGVSC_FIELDS;
use crate::consts::HGVSP_FIELDS;
use crate::consts::STRAND_FIELDS;
use crate::consts::USED_REF_FIELDS;
use crate::entry::Entry;
use crate::fields::field_named;
use crate::rules::edited_model::CODONS_FIELDS;
use crate::rules::first_field;
use crate::rules::following::terms_impact;
use crate::rules::term_pairs::transcript_row;
use crate::rules::term_pairs::IMPACT_FIELDS;
use crate::shared::blank;
use crate::shared::consequence_set;
use std::collections::BTreeSet;

/// The reference allele as the input record gives it.
pub const GIVEN_REF_FIELDS: &[&str] = &["GIVEN_REF", "given_ref"];

/// The two residues of a row, `K/T`, or the one both codons give, `S`.
pub const RESIDUE_FIELDS: &[&str] = &["Amino_acids", "amino_acids"];

/// The `BAM_EDIT` of a model whose alignment correction failed.
pub const FAILED_EDIT: &str = "FAILED";

/// The terms the two residues of a one-base substitution give, one of them on every such row
/// that is no start codon's: VariationEffect.pm reads each from the reference and the
/// alternate peptide. At the codon after the peptide's last residue an alternate stop is
/// `stop_retained_variant` beside any reference residue (:1344); beside one that is no stop
/// the class takes no Consequence there.
pub const RESIDUE_TERMS: &[&str] = &[
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "stop_retained_variant",
    "synonymous_variant",
];

/// The amino acid of each codon under the standard genetic code, a stop `*`, the codons in
/// the order `TTT TTC TTA TTG TCT ... GGG`.
const STANDARD_CODE_RESIDUES: &[u8; 64] =
    b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// The residue of one codon of the letters `A`, `C`, `G` and `T`, in either case, under the
/// standard genetic code, a stop `*`; `None` for any other text.
fn residue(codon: &str) -> Option<char> {
    let &[first, second, third] = codon.as_bytes() else {
        return None;
    };
    let place = |base: u8| match base.to_ascii_uppercase() {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    let codon = 16 * place(first)? + 4 * place(second)? + place(third)?;
    Some(char::from(STANDARD_CODE_RESIDUES[codon]))
}

/// The three-letter name HGVS gives a residue of the standard genetic code, a stop `Ter`;
/// `None` for any other letter.
fn residue_name(residue: char) -> Option<&'static str> {
    Some(match residue {
        'A' => "Ala",
        'C' => "Cys",
        'D' => "Asp",
        'E' => "Glu",
        'F' => "Phe",
        'G' => "Gly",
        'H' => "His",
        'I' => "Ile",
        'K' => "Lys",
        'L' => "Leu",
        'M' => "Met",
        'N' => "Asn",
        'P' => "Pro",
        'Q' => "Gln",
        'R' => "Arg",
        'S' => "Ser",
        'T' => "Thr",
        'V' => "Val",
        'W' => "Trp",
        'Y' => "Tyr",
        '*' => "Ter",
        _ => return None,
    })
}

/// Whether `text` is one of the four bases, in upper case.
fn one_base(text: &str) -> bool {
    matches!(text, "A" | "C" | "G" | "T")
}

/// The base that pairs with `base`, or `base` itself for another byte.
fn paired(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// Whether the entry carries `BAM_EDIT` `FAILED`.
pub(crate) fn failed_edit_row(e: &Entry) -> bool {
    BAM_EDIT_FIELDS.iter().any(|name| {
        e.fields
            .get(*name)
            .is_some_and(|value| value == FAILED_EDIT)
    })
}

/// Whether the reference entry `e` and the candidate entry `a` of one key have the row
/// shape of `use_transcript_ref_cache_cdna`: a Transcript entry of a one-base allele carrying
/// `BAM_EDIT` `FAILED` on both, the GIVEN_REF one base and the same on both, the candidate's
/// USED_REF that base and the reference's another.
pub fn failed_edit_reference_row(e: &Entry, a: &Entry) -> bool {
    if !transcript_row(e) || !one_base(&e.allele) || e.allele != a.allele {
        return false;
    }
    if !failed_edit_row(e) || !failed_edit_row(a) {
        return false;
    }
    let given = first_field(e, GIVEN_REF_FIELDS);
    let read = first_field(e, USED_REF_FIELDS);
    one_base(given)
        && given == first_field(a, GIVEN_REF_FIELDS)
        && given == first_field(a, USED_REF_FIELDS)
        && one_base(read)
        && read != given
}

/// The residues two Codons values of a row of the shape give, where the values are the
/// class's pair: the reference's reference residue, the candidate's and the alternate one
/// both share.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Residues {
    reference: char,
    candidate: char,
    alternate: char,
}

impl Residues {
    /// The Amino_acids text of a reference residue beside the alternate one: one letter
    /// where the two are equal.
    fn cell(reference: char, alternate: char) -> String {
        if reference == alternate {
            reference.to_string()
        } else {
            format!("{reference}/{alternate}")
        }
    }

    /// The term a reference residue and the alternate one give.
    fn term(reference: char, alternate: char) -> &'static str {
        match (reference == alternate, reference == '*', alternate == '*') {
            (true, true, _) => "stop_retained_variant",
            (true, false, _) => "synonymous_variant",
            (false, _, true) => "stop_gained",
            (false, true, false) => "stop_lost",
            (false, false, false) => "missense_variant",
        }
    }
}

/// The base of a strand text's transcript that the forward-strand `base` is: the base
/// itself on `1`, its pair on `-1`, `None` for any other strand.
fn on_strand(strand: &str, base: &str) -> Option<u8> {
    let base = *base.as_bytes().first()?;
    match strand {
        "1" => Some(base),
        "-1" => Some(paired(base)),
        _ => None,
    }
}

/// The residues of the Codons `expected_value` of the reference entry `e` and
/// `actual_value` of the candidate entry `a`, both of a row of the shape, where the two are
/// the class's pair: each one codon against one codon, the alternate codons equal, the
/// reference codons differing at one place alone, an upper-case base on both sides that is
/// the reference's USED_REF and the candidate's GIVEN_REF on the transcript's strand, and
/// the alternate codon the candidate's reference codon with the allele at that place.
/// `None` for any other two values and for a codon the standard genetic code does not read.
fn codon_pair(e: &Entry, a: &Entry, expected_value: &str, actual_value: &str) -> Option<Residues> {
    let (expected_reference, expected_alternate) = expected_value.split_once('/')?;
    let (actual_reference, actual_alternate) = actual_value.split_once('/')?;
    if expected_alternate != actual_alternate
        || [expected_reference, actual_reference, actual_alternate]
            .iter()
            .any(|codon| codon.len() != 3 || !codon.is_ascii())
    {
        return None;
    }
    let strand = first_field(e, STRAND_FIELDS);
    if strand != first_field(a, STRAND_FIELDS) {
        return None;
    }
    let read = on_strand(strand, first_field(e, USED_REF_FIELDS))?;
    let given = on_strand(strand, first_field(e, GIVEN_REF_FIELDS))?;
    let allele = on_strand(strand, &e.allele)?;
    let (reference, candidate, alternate) = (
        expected_reference.as_bytes(),
        actual_reference.as_bytes(),
        actual_alternate.as_bytes(),
    );
    let mut places = (0..3).filter(|place| reference[*place] != candidate[*place]);
    let place = places.next()?;
    if places.next().is_some() {
        return None;
    }
    if reference[place] != read || candidate[place] != given || alternate[place] != allele {
        return None;
    }
    // Every other base is the model's own, in lower case on both sides and in both halves.
    if (0..3).any(|other| {
        other != place
            && (!candidate[other].is_ascii_lowercase() || alternate[other] != candidate[other])
    }) {
        return None;
    }
    Some(Residues {
        reference: residue(expected_reference)?,
        candidate: residue(actual_reference)?,
        alternate: residue(actual_alternate)?,
    })
}

/// A transcript notation of one substitution as the text before its two bases, the base
/// changed and the base changed to: `NM_1.2:c.656`, `G` and `T` of `NM_1.2:c.656G>T`. `None`
/// for any other text.
fn substitution_notation(value: &str) -> Option<(&str, u8, u8)> {
    let &[.., changed, b'>', changed_to] = value.as_bytes() else {
        return None;
    };
    let bases = [changed, changed_to];
    if !bases.iter().all(|base| b"ACGT".contains(base)) {
        return None;
    }
    Some((&value[..value.len() - 3], changed, changed_to))
}

/// Whether two HGVSc values of a row of the shape differ as the two reference alleles say:
/// each one substitution, the text before the two bases the same on both sides, the base
/// changed to the allele on both, the base changed the reference's USED_REF on its side and
/// the GIVEN_REF on the candidate's, each on the transcript's strand.
/// TranscriptVariationAllele.pm:1510-1513 write the reference allele of a model that
/// carries an `_rna_edit` into the notation as the base changed.
fn transcript_notation_pair(
    e: &Entry,
    a: &Entry,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    let strand = first_field(e, STRAND_FIELDS);
    if strand != first_field(a, STRAND_FIELDS) {
        return false;
    }
    let on_this_strand = |base: &str| on_strand(strand, base);
    let (Some(read), Some(given), Some(allele)) = (
        on_this_strand(first_field(e, USED_REF_FIELDS)),
        on_this_strand(first_field(e, GIVEN_REF_FIELDS)),
        on_this_strand(&e.allele),
    ) else {
        return false;
    };
    let (Some(expected), Some(actual)) = (
        substitution_notation(expected_value),
        substitution_notation(actual_value),
    ) else {
        return false;
    };
    expected == (actual.0, read, allele) && actual.1 == given && actual.2 == allele
}

/// A protein notation of one residue as its parts: the text before `:p.`, the three letters
/// after it, the digits of the position and what follows them. `None` for a text without
/// one of the four.
fn residue_notation(value: &str) -> Option<[&str; 4]> {
    let (protein, described) = value.split_once(":p.")?;
    let residue = described.get(..3)?;
    let rest = &described[3..];
    let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
    if protein.is_empty() || digits == 0 {
        return None;
    }
    Some([protein, residue, &rest[..digits], &rest[digits..]])
}

/// Whether the parts `notation` of a protein notation describe the change of the residue
/// `reference`, which is no stop, to `alternate`: the reference's three-letter name, then
/// `=` (`%3D` in the text formats) where the two are one residue and the alternate's name
/// otherwise, a stop's `Ter`. TranscriptVariationAllele.pm:1793 takes the reference residue
/// from the reference codon, and :1911-1912 write `=` for two residues that are one.
fn describes(notation: [&str; 4], reference: char, alternate: char) -> bool {
    let [_, residue, _, change] = notation;
    if reference == '*' || residue_name(reference) != Some(residue) {
        return false;
    }
    if reference == alternate {
        return matches!(change, "=" | "%3D");
    }
    residue_name(alternate) == Some(change)
}

/// Whether two HGVSp values differ as the residues say: the two reference residues differ,
/// each value is the notation of one residue, the protein and the position the same on both
/// sides, and each side's the change of its own reference residue to the alternate one
/// (`describes`). Two reference codons that give one residue explain no difference, a sign
/// spelled `=` on one side and `%3D` on the other included.
fn protein_notation_pair(residues: Residues, expected_value: &str, actual_value: &str) -> bool {
    let (Some(expected), Some(actual)) = (
        residue_notation(expected_value),
        residue_notation(actual_value),
    ) else {
        return false;
    };
    residues.reference != residues.candidate
        && expected[0] == actual[0]
        && expected[2] == actual[2]
        && describes(expected, residues.reference, residues.alternate)
        && describes(actual, residues.candidate, residues.alternate)
}

/// A set of terms as its `RESIDUE_TERMS` and its other terms, each in the set's order.
fn residue_terms_apart<'a>(terms: &BTreeSet<&'a str>) -> (Vec<&'a str>, Vec<&'a str>) {
    terms
        .iter()
        .copied()
        .partition(|term| RESIDUE_TERMS.contains(term))
}

/// Whether two Consequence values differ as the residues say: each side holds one of
/// `RESIDUE_TERMS`, the term its own reference residue gives beside the alternate one, the
/// two terms differ and every other term is shared.
fn residue_term_pair(residues: Residues, expected_value: &str, actual_value: &str) -> bool {
    let (expected_set, actual_set) = (
        consequence_set(expected_value),
        consequence_set(actual_value),
    );
    let expected_terms: BTreeSet<&str> = expected_set.split(',').collect();
    let actual_terms: BTreeSet<&str> = actual_set.split(',').collect();
    let (expected_of_the_residues, expected_others) = residue_terms_apart(&expected_terms);
    let (actual_of_the_residues, actual_others) = residue_terms_apart(&actual_terms);
    let expected_term = Residues::term(residues.reference, residues.alternate);
    let actual_term = Residues::term(residues.candidate, residues.alternate);
    expected_term != actual_term
        && expected_of_the_residues == [expected_term]
        && actual_of_the_residues == [actual_term]
        && expected_others == actual_others
}

/// Whether the differing field `name` of the reference entry `e` and the candidate entry
/// `a` of one key, with the values `expected_value` and `actual_value`, is a cell
/// `use_transcript_ref_cache_cdna` takes. The row has the class's shape
/// (`failed_edit_reference_row`), both values are written, and the field is one of seven:
///
/// - the USED_REF, which the shape itself reads;
/// - the HGVSc that names each side's reference allele as the base changed, at one position
///   (`transcript_notation_pair`);
/// - the Codons that are the class's pair (`codon_pair`);
/// - the Amino_acids of a row whose two Codons are such a pair, each side's the text of its
///   own reference residue beside the alternate one;
/// - the Consequence of such a row whose two sets differ in the one term each side's
///   residues give (`residue_term_pair`);
/// - the IMPACT of a row whose two sets of terms are such a pair, each side's the impact of
///   its own terms (`terms_impact`);
/// - the HGVSp of a row whose two Codons are such a pair and whose two reference residues
///   differ, each side's the notation of its own reference residue changed to the alternate
///   one, neither reference residue a stop (`protein_notation_pair`).
pub fn failed_edit_reference_cell(
    e: &Entry,
    a: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
) -> bool {
    if blank(expected_value) || blank(actual_value) || !failed_edit_reference_row(e, a) {
        return false;
    }
    if field_named(name, USED_REF_FIELDS) {
        return true;
    }
    if field_named(name, HGVSC_FIELDS) {
        return transcript_notation_pair(e, a, expected_value, actual_value);
    }
    if field_named(name, CODONS_FIELDS) {
        return codon_pair(e, a, expected_value, actual_value).is_some();
    }
    let Some(residues) = codon_pair(
        e,
        a,
        first_field(e, CODONS_FIELDS),
        first_field(a, CODONS_FIELDS),
    ) else {
        return false;
    };
    if field_named(name, RESIDUE_FIELDS) {
        return expected_value == Residues::cell(residues.reference, residues.alternate)
            && actual_value == Residues::cell(residues.candidate, residues.alternate);
    }
    if field_named(name, CONSEQUENCE_FIELDS) {
        return residue_term_pair(residues, expected_value, actual_value);
    }
    if field_named(name, IMPACT_FIELDS) {
        if !residue_term_pair(residues, &e.consequence, &a.consequence) {
            return false;
        }
        let impact = |entry: &Entry| {
            let terms = consequence_set(&entry.consequence);
            terms_impact(&terms.split(',').collect::<Vec<_>>())
        };
        return impact(e).as_deref() == Some(expected_value)
            && impact(a).as_deref() == Some(actual_value);
    }
    if field_named(name, HGVSP_FIELDS) {
        return protein_notation_pair(residues, expected_value, actual_value);
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entry::Fields;

    /// A Transcript row of a failed-edit model: its terms, its Codons, its Amino_acids, its
    /// IMPACT, its strand and the reference allele it was analysed against, beside a
    /// GIVEN_REF of `A` and an allele of `C`.
    fn row(
        terms: &str,
        codons: &str,
        residues: &str,
        impact: &str,
        strand: &str,
        read: &str,
    ) -> Entry {
        let fields: Fields = [
            ("Consequence", terms),
            ("Codons", codons),
            ("Amino_acids", residues),
            ("IMPACT", impact),
            ("STRAND", strand),
            ("GIVEN_REF", "A"),
            ("USED_REF", read),
            ("BAM_EDIT", FAILED_EDIT),
        ]
        .into_iter()
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect();
        Entry {
            key: "7:1500|C|NM_990007.1|Transcript".to_string(),
            location: Some("7:1500".to_string()),
            allele: "C".to_string(),
            feature: "NM_990007.1".to_string(),
            consequence: terms.to_string(),
            fields,
            ..Entry::default()
        }
    }

    /// The reference's row, which read the base after the record's, and the candidate's.
    fn forward_pair() -> (Entry, Entry) {
        (
            row("missense_variant", "aGg/aCg", "R/T", "MODERATE", "1", "G"),
            row("missense_variant", "aAg/aCg", "K/T", "MODERATE", "1", "A"),
        )
    }

    /// Whether the class takes the cell `name` of a pair of rows.
    fn takes(pair: &(Entry, Entry), name: &str) -> bool {
        let (e, a) = pair;
        failed_edit_reference_cell(e, a, name, first_field(e, &[name]), first_field(a, &[name]))
    }

    #[test]
    fn the_reference_allele_its_codons_and_its_residues_are_taken_on_a_row_of_the_shape() {
        let pair = forward_pair();
        assert!(failed_edit_reference_row(&pair.0, &pair.1));
        for name in ["USED_REF", "Codons", "Amino_acids"] {
            assert!(takes(&pair, name), "{name}");
        }
        // The terms and the impact are equal, and a cell the class does not name is not its.
        for name in ["Consequence", "IMPACT", "STRAND", "GIVEN_REF"] {
            assert!(!takes(&pair, name), "{name}");
        }
        // The two engines swapped are no row of the shape.
        assert!(!failed_edit_reference_row(&pair.1, &pair.0));
        assert!(!takes(&(pair.1.clone(), pair.0.clone()), "USED_REF"));
    }

    #[test]
    fn a_row_without_the_shape_keeps_every_cell() {
        let changed = |name: &str, value: &str, both: bool| {
            let (mut e, mut a) = forward_pair();
            a.fields.insert(name.to_string(), value.to_string());
            if both {
                e.fields.insert(name.to_string(), value.to_string());
            }
            (e, a)
        };
        let mut cases = vec![
            // A model the alignment corrected, and one without a status.
            changed("BAM_EDIT", "OK", true),
            changed("BAM_EDIT", "OK", false),
            changed("BAM_EDIT", "-", true),
            // A candidate that read another base than the record's, and a record of two.
            changed("USED_REF", "T", false),
            changed("GIVEN_REF", "AT", true),
            changed("GIVEN_REF", "T", false),
        ];
        // An allele that is no single base, and a row that is no transcript's.
        let (mut e, mut a) = forward_pair();
        e.allele = "-".to_string();
        a.allele = "-".to_string();
        cases.push((e, a));
        let (mut e, mut a) = forward_pair();
        for entry in [&mut e, &mut a] {
            entry.key = "7:1500|C|ENSR990007|RegulatoryFeature".to_string();
        }
        cases.push((e, a));
        // A reference that read the record's own base differs nowhere in USED_REF.
        let (mut e, a) = forward_pair();
        e.fields.insert("USED_REF".to_string(), "A".to_string());
        cases.push((e, a));
        for (at, pair) in cases.iter().enumerate() {
            assert!(!failed_edit_reference_row(&pair.0, &pair.1), "case {at}");
            for name in ["USED_REF", "Codons", "Amino_acids"] {
                assert!(!takes(pair, name), "case {at} {name}");
            }
        }
    }

    #[test]
    fn the_codons_are_held_to_the_two_reference_alleles_on_the_transcripts_strand() {
        let codons = |strand: &str, read: &str, reference: &str, candidate: &str| {
            let e = row("missense_variant", reference, "-", "MODERATE", strand, read);
            let a = row("missense_variant", candidate, "-", "MODERATE", strand, "A");
            codon_pair(&e, &a, reference, candidate)
        };
        assert_eq!(
            codons("1", "G", "aGg/aCg", "aAg/aCg"),
            Some(Residues {
                reference: 'R',
                candidate: 'K',
                alternate: 'T'
            })
        );
        // On the reverse strand the codon holds the pairs of the three bases.
        assert_eq!(
            codons("-1", "G", "Cga/Gga", "Tga/Gga"),
            Some(Residues {
                reference: 'R',
                candidate: '*',
                alternate: 'G'
            })
        );
        for (strand, read, reference, candidate) in [
            // The reference's upper-case base is not the base it read.
            ("1", "T", "aGg/aCg", "aAg/aCg"),
            // The candidate's is not the record's.
            ("1", "G", "aGg/aCg", "aTg/aCg"),
            // The alternate codons differ, or do not hold the allele at the place.
            ("1", "G", "aGg/aCg", "aAg/aTg"),
            ("1", "G", "aGg/aTg", "aAg/aTg"),
            // A second base differs, or a base beside the place is in upper case.
            ("1", "G", "cGg/aCg", "aAg/aCg"),
            ("1", "G", "AGg/ACg", "AAg/ACg"),
            // The strand is read as the forward one.
            ("-1", "G", "aGg/aCg", "aAg/aCg"),
            ("0", "G", "aGg/aCg", "aAg/aCg"),
            // More than one codon, no alternate half, and a letter that is no base.
            ("1", "G", "aGgttt/aCgttt", "aAgttt/aCgttt"),
            ("1", "G", "aGg", "aAg"),
            ("1", "G", "nGg/nCg", "nAg/nCg"),
        ] {
            assert_eq!(
                codons(strand, read, reference, candidate),
                None,
                "{strand} {read} {reference} {candidate}"
            );
        }
    }

    #[test]
    fn the_residues_are_each_sides_own_codons_translated() {
        let pair = forward_pair();
        let with_residues = |reference: &str, candidate: &str| {
            let (mut e, mut a) = pair.clone();
            e.fields
                .insert("Amino_acids".to_string(), reference.to_string());
            a.fields
                .insert("Amino_acids".to_string(), candidate.to_string());
            takes(&(e, a), "Amino_acids")
        };
        assert!(with_residues("R/T", "K/T"));
        for (reference, candidate) in [("K/T", "R/T"), ("R/T", "K/S"), ("R", "K/T"), ("R/T", "K")] {
            assert!(
                !with_residues(reference, candidate),
                "{reference} {candidate}"
            );
        }
        // Without the pair of Codons beside them the residues are tied to nothing.
        let (mut e, a) = pair.clone();
        e.fields.insert("Codons".to_string(), "aTg/aCg".to_string());
        assert!(!takes(&(e, a), "Amino_acids"));
    }

    #[test]
    fn the_term_and_the_impact_are_the_ones_each_sides_residues_give() {
        // The base the reference read is the allele itself, so its two codons are one.
        let unchanged = (
            row("synonymous_variant", "aCg/aCg", "T", "LOW", "1", "C"),
            row("missense_variant", "aAg/aCg", "K/T", "MODERATE", "1", "A"),
        );
        for name in ["USED_REF", "Codons", "Amino_acids", "Consequence", "IMPACT"] {
            assert!(takes(&unchanged, name), "{name}");
        }
        // A term both sides share stands beside the pair, and the impact follows the terms.
        let spliced = (
            row(
                "splice_region_variant,synonymous_variant",
                "aCg/aCg",
                "T",
                "LOW",
                "1",
                "C",
            ),
            row(
                "missense_variant,splice_region_variant",
                "aAg/aCg",
                "K/T",
                "MODERATE",
                "1",
                "A",
            ),
        );
        assert!(takes(&spliced, "Consequence"));
        assert!(takes(&spliced, "IMPACT"));
        // The candidate's codon is a stop codon and the reference's, with the base it read, is
        // none: `stop_lost` against `missense_variant`, and the impact of each.
        let stop = (
            row("missense_variant", "tTa/tCa", "L/S", "MODERATE", "1", "T"),
            row("stop_lost", "tAa/tCa", "*/S", "HIGH", "1", "A"),
        );
        for name in ["Codons", "Amino_acids", "Consequence", "IMPACT"] {
            assert!(takes(&stop, name), "{name}");
        }
        // Two codons that give one term differ in no term the class takes.
        let both_stops = (
            row("stop_lost", "tGa/tCa", "*/S", "HIGH", "1", "G"),
            row("stop_lost", "tAa/tCa", "*/S", "HIGH", "1", "A"),
        );
        assert!(takes(&both_stops, "Codons"));
        assert!(!takes(&both_stops, "Consequence"));
        assert!(!takes(&both_stops, "IMPACT"));
        // A term that is not the residues', on either side, stays a difference.
        for (reference, candidate) in [
            ("missense_variant", "missense_variant"),
            ("stop_gained", "missense_variant"),
            ("synonymous_variant", "stop_gained"),
            ("synonymous_variant", "missense_variant,start_lost"),
            ("synonymous_variant,intron_variant", "missense_variant"),
            ("coding_sequence_variant", "missense_variant"),
        ] {
            let (mut e, mut a) = unchanged.clone();
            e.consequence = reference.to_string();
            e.fields
                .insert("Consequence".to_string(), reference.to_string());
            a.consequence = candidate.to_string();
            a.fields
                .insert("Consequence".to_string(), candidate.to_string());
            let pair = (e, a);
            assert!(!takes(&pair, "Consequence"), "{reference} {candidate}");
            assert!(!takes(&pair, "IMPACT"), "{reference} {candidate}");
        }
        // An impact that is not its own terms' stays a difference beside the pair.
        let (mut e, a) = unchanged.clone();
        e.fields
            .insert("IMPACT".to_string(), "MODIFIER".to_string());
        assert!(!takes(&(e, a), "IMPACT"));
    }

    /// A pair of rows with the HGVSc and the HGVSp each side writes under `--hgvs`.
    fn notated(
        pair: &(Entry, Entry),
        reference: (&str, &str),
        candidate: (&str, &str),
    ) -> (Entry, Entry) {
        let (mut e, mut a) = pair.clone();
        for (entry, (transcript, protein)) in [(&mut e, reference), (&mut a, candidate)] {
            entry
                .fields
                .insert("HGVSc".to_string(), transcript.to_string());
            entry
                .fields
                .insert("HGVSp".to_string(), protein.to_string());
        }
        (e, a)
    }

    /// The reverse-strand pair: the reference read the base before the record's, the pair of
    /// which stands in its codon, and the allele's pair makes a stop codon on both sides.
    fn reverse_pair() -> (Entry, Entry) {
        (
            row("stop_gained", "tCa/tGa", "S/*", "HIGH", "-1", "G"),
            row("stop_gained", "tTa/tGa", "L/*", "HIGH", "-1", "A"),
        )
    }

    #[test]
    fn the_transcript_notation_names_the_base_each_side_read_at_one_position() {
        let forward = forward_pair();
        let takes_notations = |pair: &(Entry, Entry), reference: &str, candidate: &str| {
            takes(&notated(pair, (reference, "-"), (candidate, "-")), "HGVSc")
        };
        for (reference, candidate) in [
            ("NM_990007.1:c.58G>C", "NM_990007.1:c.58A>C"),
            ("NM_990007.1:c.*12G>C", "NM_990007.1:c.*12A>C"),
            ("NM_990007.1:c.-7G>C", "NM_990007.1:c.-7A>C"),
            ("NR_990007.1:n.157G>C", "NR_990007.1:n.157A>C"),
        ] {
            assert!(
                takes_notations(&forward, reference, candidate),
                "{reference} {candidate}"
            );
        }
        // The base the reference read is the allele itself: its notation names no change.
        let unchanged = (
            row("synonymous_variant", "aCg/aCg", "T", "LOW", "1", "C"),
            forward.1.clone(),
        );
        assert!(takes_notations(
            &unchanged,
            "NM_990007.1:c.58C>C",
            "NM_990007.1:c.58A>C"
        ));
        // On the reverse strand the notation holds the pairs of the three bases.
        let reverse = reverse_pair();
        assert!(takes_notations(
            &reverse,
            "NM_990007.1:c.59C>G",
            "NM_990007.1:c.59T>G"
        ));
        assert!(!takes_notations(
            &reverse,
            "NM_990007.1:c.59G>C",
            "NM_990007.1:c.59A>C"
        ));
        for (reference, candidate) in [
            // Another position, another numbering and another model stay differences.
            ("NM_990007.1:c.*13G>C", "NM_990007.1:c.*12A>C"),
            ("NM_990007.1:n.58G>C", "NM_990007.1:c.58A>C"),
            ("NM_990007.2:c.58G>C", "NM_990007.1:c.58A>C"),
            // The reference's base is not the base it read, or the candidate's not the record's.
            ("NM_990007.1:c.58T>C", "NM_990007.1:c.58A>C"),
            ("NM_990007.1:c.58G>C", "NM_990007.1:c.58T>C"),
            // The base changed to is not the allele, on either side.
            ("NM_990007.1:c.58G>T", "NM_990007.1:c.58A>T"),
            ("NM_990007.1:c.58G>T", "NM_990007.1:c.58A>C"),
            ("NM_990007.1:c.58G>C", "NM_990007.1:c.58A>T"),
            // The two engines' notations swapped.
            ("NM_990007.1:c.58A>C", "NM_990007.1:c.58G>C"),
            // A notation that is no substitution of one base.
            ("NM_990007.1:c.58del", "NM_990007.1:c.58A>C"),
            ("NM_990007.1:c.58_59delinsCC", "NM_990007.1:c.58A>C"),
            ("NM_990007.1:c.58G>C", "NM_990007.1:c.58_59delinsCC"),
            ("G>C", "A>"),
        ] {
            assert!(
                !takes_notations(&forward, reference, candidate),
                "{reference} {candidate}"
            );
        }
        // A row without the shape has no notation of the class.
        let (mut e, a) = notated(
            &forward,
            ("NM_990007.1:c.58G>C", "-"),
            ("NM_990007.1:c.58A>C", "-"),
        );
        e.fields.insert("BAM_EDIT".to_string(), "OK".to_string());
        assert!(!takes(&(e, a), "HGVSc"));
    }

    #[test]
    fn the_protein_notation_names_the_residue_of_each_sides_own_reference_codon() {
        let forward = forward_pair();
        let takes_notations = |pair: &(Entry, Entry), reference: &str, candidate: &str| {
            takes(&notated(pair, ("-", reference), ("-", candidate)), "HGVSp")
        };
        assert!(takes_notations(
            &forward,
            "NP_990007.1:p.Arg20Thr",
            "NP_990007.1:p.Lys20Thr"
        ));
        // The reference's two codons are one: its notation is the unchanged residue's, with
        // the sign as the text formats or the JSON format write it.
        let unchanged = (
            row("synonymous_variant", "aCg/aCg", "T", "LOW", "1", "C"),
            forward.1.clone(),
        );
        for sign in ["%3D", "="] {
            assert!(
                takes_notations(
                    &unchanged,
                    &format!("NP_990007.1:p.Thr20{sign}"),
                    "NP_990007.1:p.Lys20Thr"
                ),
                "{sign}"
            );
        }
        assert!(!takes_notations(
            &unchanged,
            "NP_990007.1:p.Thr20Thr",
            "NP_990007.1:p.Lys20Thr"
        ));
        // A stop codon both alternate codons spell.
        assert!(takes_notations(
            &reverse_pair(),
            "NP_990007.1:p.Ser20Ter",
            "NP_990007.1:p.Leu20Ter"
        ));
        for (reference, candidate) in [
            // Another position and another protein stay differences.
            ("NP_990007.1:p.Arg21Thr", "NP_990007.1:p.Lys20Thr"),
            ("NP_990007.2:p.Arg20Thr", "NP_990007.1:p.Lys20Thr"),
            // A residue that is not the translation of the side's own codons.
            ("NP_990007.1:p.Gly20Thr", "NP_990007.1:p.Lys20Thr"),
            ("NP_990007.1:p.Arg20Thr", "NP_990007.1:p.Gln20Thr"),
            ("NP_990007.1:p.Arg20Ser", "NP_990007.1:p.Lys20Thr"),
            ("NP_990007.1:p.Arg20Thr", "NP_990007.1:p.Lys20Ser"),
            ("NP_990007.1:p.Arg20%3D", "NP_990007.1:p.Lys20Thr"),
            // The two engines' notations swapped.
            ("NP_990007.1:p.Lys20Thr", "NP_990007.1:p.Arg20Thr"),
            // Another form of the notation.
            ("NP_990007.1:p.R20T", "NP_990007.1:p.K20T"),
            ("NP_990007.1:p.(Arg20Thr)", "NP_990007.1:p.(Lys20Thr)"),
            ("NP_990007.1:p.Arg20_Thr21del", "NP_990007.1:p.Lys20Thr"),
            ("NP_990007.1:p.ArgThr", "NP_990007.1:p.LysThr"),
            ("p.Arg20Thr", "NP_990007.1:p.Lys20Thr"),
        ] {
            assert!(
                !takes_notations(&forward, reference, candidate),
                "{reference} {candidate}"
            );
        }
        // A reference residue that is a stop is written with an extension the class does
        // not read, on either side.
        let stop = (
            row("missense_variant", "tTa/tCa", "L/S", "MODERATE", "1", "T"),
            row("stop_lost", "tAa/tCa", "*/S", "HIGH", "1", "A"),
        );
        assert!(!takes_notations(
            &stop,
            "NP_990007.1:p.Leu20Ser",
            "NP_990007.1:p.Ter20SerextTer5"
        ));
        assert!(!takes_notations(
            &stop,
            "NP_990007.1:p.Leu20Ser",
            "NP_990007.1:p.Ter20Ser"
        ));
        // Two reference codons that give one residue differ in no notation the class takes.
        let one_residue = (
            row("synonymous_variant", "ccG/ccC", "P", "LOW", "1", "G"),
            row("synonymous_variant", "ccA/ccC", "P", "LOW", "1", "A"),
        );
        assert!(takes(&one_residue, "Codons"));
        assert!(!takes_notations(
            &one_residue,
            "NP_990007.1:p.Pro20%3D",
            "NP_990007.1:p.Pro20="
        ));
        // Without the pair of Codons beside it the notation is tied to nothing.
        let (mut e, a) = notated(
            &forward,
            ("-", "NP_990007.1:p.Arg20Thr"),
            ("-", "NP_990007.1:p.Lys20Thr"),
        );
        e.fields.insert("Codons".to_string(), "aTg/aCg".to_string());
        assert!(!takes(&(e, a), "HGVSp"));
    }

    #[test]
    fn each_residue_of_the_standard_code_has_its_three_letter_name() {
        for (residue, name) in [
            ('A', "Ala"),
            ('C', "Cys"),
            ('D', "Asp"),
            ('E', "Glu"),
            ('F', "Phe"),
            ('G', "Gly"),
            ('H', "His"),
            ('I', "Ile"),
            ('K', "Lys"),
            ('L', "Leu"),
            ('M', "Met"),
            ('N', "Asn"),
            ('P', "Pro"),
            ('Q', "Gln"),
            ('R', "Arg"),
            ('S', "Ser"),
            ('T', "Thr"),
            ('V', "Val"),
            ('W', "Trp"),
            ('Y', "Tyr"),
            ('*', "Ter"),
        ] {
            assert_eq!(residue_name(residue), Some(name), "{residue}");
        }
        let mut named: Vec<char> = STANDARD_CODE_RESIDUES
            .iter()
            .map(|residue| char::from(*residue))
            .collect();
        named.sort_unstable();
        named.dedup();
        assert_eq!(named.len(), 21);
        for residue in named {
            assert!(residue_name(residue).is_some(), "{residue}");
        }
        for unnamed in ['X', 'U', 'B', '-', 'a'] {
            assert_eq!(residue_name(unnamed), None, "{unnamed}");
        }
    }

    #[test]
    fn two_residues_give_one_term() {
        for (reference, alternate, term) in [
            ('K', 'T', "missense_variant"),
            ('T', 'T', "synonymous_variant"),
            ('K', '*', "stop_gained"),
            ('*', 'K', "stop_lost"),
            ('*', '*', "stop_retained_variant"),
        ] {
            assert_eq!(Residues::term(reference, alternate), term);
            assert!(RESIDUE_TERMS.contains(&term));
        }
        assert_eq!(Residues::cell('K', 'T'), "K/T");
        assert_eq!(Residues::cell('T', 'T'), "T");
    }

    #[test]
    fn one_codon_is_read_under_the_standard_code() {
        for (codon, read) in [
            ("ATG", 'M'),
            ("tgg", 'W'),
            ("aAg", 'K'),
            ("TAA", '*'),
            ("TAG", '*'),
            ("TGA", '*'),
            ("TTT", 'F'),
            ("GGG", 'G'),
            ("AGA", 'R'),
            ("ATA", 'I'),
        ] {
            assert_eq!(residue(codon), Some(read), "{codon}");
        }
        let stops = STANDARD_CODE_RESIDUES
            .iter()
            .filter(|read| **read == b'*')
            .count();
        assert_eq!(stops, 3);
        for unread in ["", "AT", "ATGA", "ANG", "A-G", "AUG"] {
            assert_eq!(residue(unread), None, "{unread}");
        }
    }
}
