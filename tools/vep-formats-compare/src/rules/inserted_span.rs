// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reference codons of an insertion on the bases an edited model holds and the genome lacks.

use crate::consts::AMINO_ACIDS_FIELDS;
use crate::consts::USED_REF_FIELDS;
use crate::entry::Entry;
use crate::error::Error;
use crate::error::Result;
use crate::fields::field_named;
use crate::genetic_codes::genetic_code;
use crate::genetic_codes::Reading;
use crate::pyfmt;
use crate::rules::downstream::bound_number;
use crate::rules::downstream::bound_texts;
use crate::rules::downstream::DownstreamSequences;
use crate::rules::edited_model::CODONS_FIELDS;
use crate::rules::first_field;
use crate::rules::term_pairs::cds_bounds;
use crate::rules::term_pairs::edited_model_row;
use crate::rules::term_pairs::transcript_row;
use crate::shared::blank;

/// The two readings of a row's reference codons where its coding span `cds` is bases the
/// model holds and the genome lacks, both in lower case and the reference's first.
///
/// TranscriptVariationAllele.pm:841-843 fix the window: it opens at the first base of the
/// row's first residue and is three bases a residue long. The candidate reads the window
/// from `coding`, the model's coding sequence. For the reference allele of a model carrying
/// `_rna_edit`, :872-873 first write that allele over the span, and an insertion's is empty,
/// so :877 reads the window from a sequence that has lost the span: the bases before the
/// span, then as many bases past it as the window still holds, `substr` giving what the
/// sequence has left.
///
/// `None` for a span or residues that are no range counted from 1, a span that is not
/// inside both the window and the sequence, and a sequence that is not ASCII.
pub fn reference_codon_readings(
    coding: &str,
    cds: (i64, i64),
    residues: (i64, i64),
) -> Option<(String, String)> {
    let (span_start, span_end) = cds;
    let (first_residue, last_residue) = residues;
    if !coding.is_ascii() || span_start < 1 || first_residue < 1 || last_residue < first_residue {
        return None;
    }
    let window_start = first_residue.checked_mul(3)?.checked_sub(2)?;
    let window_end = last_residue.checked_mul(3)?;
    if span_start < window_start || span_start > span_end || span_end > window_end {
        return None;
    }
    let place = |position: i64| usize::try_from(position).ok();
    let (span_start, span_end) = (place(span_start)?, place(span_end)?);
    let (window_start, window_end) = (place(window_start)?, place(window_end)?);
    if span_end > coding.len() {
        return None;
    }
    let window_length = window_end - window_start + 1;
    let window_of = |sequence: &str| {
        let from = window_start - 1;
        let to = sequence.len().min(from + window_length);
        sequence[from..to].to_ascii_lowercase()
    };
    let without_the_span = format!("{}{}", &coding[..span_start - 1], &coding[span_end..]);
    Some((window_of(&without_the_span), window_of(coding)))
}

/// The first and last residue of the entry's span: the `Protein_position` of the text
/// formats (one position or a range) or the JSON `protein_start` and `protein_end`; `None`
/// where a bound is missing or is not a number.
fn residue_bounds(e: &Entry) -> Option<(i64, i64)> {
    let [first, last] = bound_texts(e, "Protein_position", ["protein_start", "protein_end"]);
    Some((bound_number(first)?, bound_number(last)?))
}

/// Whether the allele is the bases of an insertion or a substitution: letters of a
/// sequence and no `-`.
fn allele_of_bases(allele: &str) -> bool {
    !allele.is_empty() && allele.bytes().all(|letter| b"ACGTN".contains(&letter))
}

/// The residues of `codons`, whole codons of the letters `a`, `c`, `g` and `t`, under the
/// genetic code `table` as `reading` reads it, a stop as `*`; `None` for a table no code
/// is held for, for bases that are no whole number of codons and for a codon with any
/// other letter.
fn translation(codons: &str, table: i64, reading: Reading) -> Option<String> {
    let code = genetic_code(table, reading)?;
    let bases = codons.to_ascii_uppercase();
    if !bases.len().is_multiple_of(3) {
        return None;
    }
    bases
        .as_bytes()
        .chunks_exact(3)
        .map(|codon| code.amino_acid(codon).map(char::from))
        .collect()
}

/// The residues Ensembl VEP writes for the reference codons `codons` under the genetic code
/// `table`: TranscriptVariationAllele.pm:724-736 translate the whole codons and :774-776
/// add `X` for a last base or two, unless the whole codons read one stop alone. `None`
/// where `translation` has none for the whole codons.
fn reference_residues(codons: &str, table: i64) -> Option<String> {
    let whole = codons.len() - codons.len() % 3;
    let mut residues = translation(codons.get(..whole)?, table, Reading::EnsemblVep)?;
    if whole < codons.len() && residues != "*" {
        residues.push('X');
    }
    Some(residues)
}

/// Whether the differing field `name` of the reference entry `e` and the candidate entry
/// `a` of one key, with the values `expected_value` and `actual_value`, is a cell
/// `refseq_edited_model_inserted_span_codons` takes.
///
/// The row is a Transcript entry of an allele of bases carrying `BAM_EDIT` `OK` on both
/// sides, its USED_REF empty on both (an insertion's reference allele), its coding span
/// and its residues the same numbers on both, and `sequences` holds its transcript. The
/// two readings of the reference codons are derived from that transcript's coding
/// sequence (`reference_codon_readings`) and differ. The field is one of two:
///
/// - the Codons whose alternate halves are equal and whose reference halves are the two
///   readings, the reference's the one that has lost the span;
/// - the Amino_acids of a row whose two Codons are such a pair: the alternate halves
///   equal, the reference halves different, each the translation of its own side's
///   reading under the transcript's genetic code as that engine reads it, the reference's
///   as `reference_residues` writes it. A transcript whose cache record names no code a
///   table is held for ends the comparison with the error that names it.
pub fn inserted_span_codons_cell(
    e: &Entry,
    a: &Entry,
    name: &str,
    expected_value: &str,
    actual_value: &str,
    sequences: &DownstreamSequences,
) -> Result<bool> {
    let codons_cell = field_named(name, CODONS_FIELDS);
    if !codons_cell && !field_named(name, AMINO_ACIDS_FIELDS) {
        return Ok(false);
    }
    if blank(expected_value) || blank(actual_value) {
        return Ok(false);
    }
    if !transcript_row(e) || !allele_of_bases(&e.allele) {
        return Ok(false);
    }
    if !edited_model_row(e) || !edited_model_row(a) {
        return Ok(false);
    }
    if !blank(first_field(e, USED_REF_FIELDS)) || !blank(first_field(a, USED_REF_FIELDS)) {
        return Ok(false);
    }
    let Some(cds) = cds_bounds(e).filter(|bounds| cds_bounds(a) == Some(*bounds)) else {
        return Ok(false);
    };
    let Some(residues) = residue_bounds(e).filter(|bounds| residue_bounds(a) == Some(*bounds))
    else {
        return Ok(false);
    };
    let Some(transcript) = sequences.transcript(e)? else {
        return Ok(false);
    };
    let Some((without_the_span, model)) =
        reference_codon_readings(&transcript.coding, cds, residues)
    else {
        return Ok(false);
    };
    if without_the_span == model {
        return Ok(false);
    }
    let (expected_codons, actual_codons) = if codons_cell {
        (expected_value, actual_value)
    } else {
        (first_field(e, CODONS_FIELDS), first_field(a, CODONS_FIELDS))
    };
    // TranscriptVariationAllele.pm:680 joins the two codon strings with a slash, always.
    let (expected_reference, expected_slash, expected_alternate) =
        pyfmt::partition(expected_codons, "/");
    let (actual_reference, actual_slash, actual_alternate) = pyfmt::partition(actual_codons, "/");
    if expected_slash.is_empty()
        || actual_slash.is_empty()
        || expected_alternate != actual_alternate
        || expected_reference != written(&without_the_span)
        || actual_reference != model
    {
        return Ok(false);
    }
    if codons_cell {
        return Ok(true);
    }
    let (expected_residues, expected_alternate) = residue_halves(expected_value);
    let (actual_residues, actual_alternate) = residue_halves(actual_value);
    if expected_alternate != actual_alternate || expected_residues == actual_residues {
        return Ok(false);
    }
    let table = transcript
        .table
        .held(&e.feature)
        .map_err(Error::comparator)?;
    let reads = |translated: Option<String>, residues: &str| {
        translated.is_some_and(|read| written(&read) == residues)
    };
    Ok(reads(
        reference_residues(&without_the_span, table),
        expected_residues,
    ) && reads(translation(&model, table, Reading::VepRs), actual_residues))
}

/// A reading as the row writes it: TranscriptVariationAllele.pm:878-880 write `-` for
/// codons, and for their residues, where the window holds no base.
fn written(reading: &str) -> &str {
    if reading.is_empty() {
        "-"
    } else {
        reading
    }
}

/// The reference and the alternate half of an Amino_acids value:
/// TranscriptVariationAllele.pm:629 writes two equal halves as one, without a slash.
fn residue_halves(value: &str) -> (&str, &str) {
    let (reference, slash, alternate) = pyfmt::partition(value, "/");
    if slash.is_empty() {
        (reference, reference)
    } else {
        (reference, alternate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Nine codons whose bases are all told apart by their place: no cut of one, two or
    /// three bases leaves a window of the sequence as it was.
    const CODING: &str = "ATGAAAACTCTTGCTGTCCGTGGTTAA";

    #[test]
    fn the_reference_reads_the_window_from_a_sequence_that_has_lost_the_span() {
        // Coding bases 9 and 10 in residues 3 and 4: the window is bases 7 to 12.
        assert_eq!(
            reference_codon_readings(CODING, (9, 10), (3, 4)),
            Some(("acttgc".to_string(), "actctt".to_string()))
        );
        // A span that is one whole codon: the reference reads the next codon.
        assert_eq!(
            reference_codon_readings(CODING, (10, 12), (4, 4)),
            Some(("gct".to_string(), "ctt".to_string()))
        );
        // A span of one base that opens the window.
        assert_eq!(
            reference_codon_readings(CODING, (4, 4), (2, 2)),
            Some(("aaa".to_string(), "aaa".to_string()))
        );
        assert_eq!(
            reference_codon_readings(CODING, (7, 7), (3, 3)),
            Some(("ctc".to_string(), "act".to_string()))
        );
        // The model's letters are read in either case.
        assert_eq!(
            reference_codon_readings(&CODING.to_ascii_lowercase(), (9, 10), (3, 4)),
            Some(("acttgc".to_string(), "actctt".to_string()))
        );
    }

    #[test]
    fn a_window_past_the_end_of_the_sequence_holds_what_is_left() {
        // The last two codons, less two bases of the last but one: 27 - 2 = 25 bases are
        // left, of which the window, opening at base 22, reads four.
        assert_eq!(
            reference_codon_readings(CODING, (23, 24), (8, 9)),
            Some(("gtaa".to_string(), "ggttaa".to_string()))
        );
        // The whole of the last codon: nothing is left to read.
        assert_eq!(
            reference_codon_readings(CODING, (25, 27), (9, 9)),
            Some((String::new(), "taa".to_string()))
        );
    }

    #[test]
    fn a_span_outside_its_window_or_its_sequence_has_no_readings() {
        for (cds, residues) in [
            // Before the window and after it.
            ((6, 8), (3, 3)),
            ((8, 10), (3, 3)),
            // Written backwards, and counted from zero or below.
            ((10, 9), (3, 4)),
            ((0, 2), (1, 1)),
            ((-2, 2), (1, 1)),
            ((1, 2), (0, 1)),
            ((4, 5), (2, 1)),
            // Past the end of the sequence.
            ((28, 29), (10, 10)),
            ((26, 28), (9, 10)),
            // Residues past what 64 bits count.
            ((1, 2), (1, i64::MAX)),
            ((i64::MAX - 1, i64::MAX), (i64::MAX / 3, i64::MAX / 3)),
        ] {
            assert_eq!(
                reference_codon_readings(CODING, cds, residues),
                None,
                "{cds:?} {residues:?}"
            );
        }
        assert_eq!(reference_codon_readings("ATGAAÄACT", (4, 5), (2, 2)), None);
    }

    #[test]
    fn whole_codons_are_translated_under_the_named_code_as_each_engine_reads_it() {
        for reading in [Reading::EnsemblVep, Reading::VepRs] {
            assert_eq!(translation("actctt", 1, reading).as_deref(), Some("TL"));
            assert_eq!(translation("ACTtaa", 1, reading).as_deref(), Some("T*"));
            // The vertebrate mitochondrial code reads two codons unlike the standard one.
            assert_eq!(translation("atatga", 1, reading).as_deref(), Some("I*"));
            assert_eq!(translation("atatga", 2, reading).as_deref(), Some("MW"));
            assert_eq!(translation("", 1, reading).as_deref(), Some(""));
            for no_codons in ["actct", "actcnt", "act-tt"] {
                assert_eq!(translation(no_codons, 1, reading), None, "{no_codons}");
            }
            assert_eq!(translation("actctt", 7, reading), None);
        }
        // Genetic code 26 is vep-rs's alone: Ensembl VEP reads its codons under another.
        assert_eq!(translation("ctg", 26, Reading::VepRs).as_deref(), Some("A"));
        assert_eq!(
            translation("ctg", 26, Reading::EnsemblVep).as_deref(),
            Some("L")
        );
    }

    #[test]
    fn the_reference_residues_end_in_x_for_a_last_base_or_two() {
        assert_eq!(reference_residues("acttgc", 1).as_deref(), Some("TC"));
        assert_eq!(reference_residues("gtaa", 1).as_deref(), Some("VX"));
        assert_eq!(reference_residues("gtaaa", 1).as_deref(), Some("VX"));
        assert_eq!(reference_residues("t", 1).as_deref(), Some("X"));
        assert_eq!(reference_residues("", 1).as_deref(), Some(""));
        // One stop alone takes no `X`; a stop beside another residue does.
        assert_eq!(reference_residues("taag", 1).as_deref(), Some("*"));
        assert_eq!(reference_residues("taatgac", 1).as_deref(), Some("**X"));
        assert_eq!(reference_residues("gtataag", 1).as_deref(), Some("V*X"));
        // A stop of the standard code that the vertebrate mitochondrial one reads as a residue.
        assert_eq!(reference_residues("tgag", 1).as_deref(), Some("*"));
        assert_eq!(reference_residues("tgag", 2).as_deref(), Some("WX"));
        assert_eq!(reference_residues("gtna", 1), None);
        assert_eq!(reference_residues("gtaa", 7), None);
        // The whole codons end inside a letter that is not ASCII.
        assert_eq!(reference_residues("ggÄa", 1), None);
    }

    #[test]
    fn an_inserted_allele_is_bases_alone() {
        for bases in ["A", "GG", "ACGTN"] {
            assert!(allele_of_bases(bases), "{bases}");
        }
        for other in ["", "-", "A-", "deletion", "acgt"] {
            assert!(!allele_of_bases(other), "{other}");
        }
    }

    #[test]
    fn an_empty_reading_is_written_as_a_hyphen_and_equal_residues_as_one_half() {
        assert_eq!(written(""), "-");
        assert_eq!(written("gct"), "gct");
        assert_eq!(residue_halves("TC/TX"), ("TC", "TX"));
        assert_eq!(residue_halves("A"), ("A", "A"));
        assert_eq!(residue_halves("A/"), ("A", ""));
        assert_eq!(residue_halves("-/X"), ("-", "X"));
    }
}
