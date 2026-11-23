// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row class of an insertion or deletion that overlaps the start codon of a
//! transcript with no 5' UTR and leaves a start codon at the start of the coding
//! sequence: `start_codon_kept_without_five_prime_utr`.

use crate::consts::HGVSP_FIELDS;
use crate::entry::Entry;
use crate::error::Result;
use crate::fields::field_named;
use crate::pyfmt;
use crate::rules::following::term_pair_following_cell;
use crate::rules::term_pairs::transcript_row;
use crate::rules::term_pairs::INFRAME_DELETION_TERM;
use crate::shared::consequence_set;
use crate::shared::transcript_positions;
use crate::shared::FRAMESHIFT_TERM;
use crate::shared::INFRAME_INSERTION_TERM;
use std::collections::BTreeSet;

pub const START_LOST_TERM: &str = "start_lost";

pub const START_RETAINED_TERM: &str = "start_retained_variant";

/// The terms of an insertion or deletion inside the coding sequence. Membership only.
pub const LENGTH_TERMS: &[&str] = &[
    FRAMESHIFT_TERM,
    INFRAME_INSERTION_TERM,
    INFRAME_DELETION_TERM,
];

/// Whether `reference` is the protein notation of a lost start and `candidate` the
/// notation of the change itself, both naming one protein. A row that carries `start_lost`
/// is written with `?` for its alternate (TranscriptVariationAllele.pm:2161-2165):
/// `p.Leu2?` where the notation's two positions are one and `p.Ser3_?2` where they are
/// two (:2030-2037). The candidate's row carries no `start_lost`, so its notation is the
/// frameshift's or the deletion's; one that holds a `?`, the count of a frameshift whose
/// new stop is not predicted (:2005-2006), is refused.
pub fn lost_start_against_the_change(reference: &str, candidate: &str) -> bool {
    let (protein, separator, lost) = pyfmt::partition(reference, ":p.");
    let (named, parted, change) = pyfmt::partition(candidate, ":p.");
    !separator.is_empty()
        && !parted.is_empty()
        && !protein.is_empty()
        && protein == named
        && lost.contains('?')
        && !change.is_empty()
        && !change.contains('?')
}

/// The CDS start and end of an entry whose cDNA span is its CDS span and opens on a base
/// of the start codon: the coding sequence opens on the transcript's first base, and the
/// smaller of the two bounds is 1, 2 or 3. A JSON entry can write the start key above the
/// end key, so the span opens on the smaller. `None` for any other entry, a bound that is
/// no number without a leading zero included.
fn span_on_the_start_codon(e: &Entry) -> Option<(String, String)> {
    let (cdna_start, cdna_end, cds_start, cds_end) = transcript_positions(e)?;
    let number = |text: &str| {
        !text.is_empty()
            && !text.starts_with('0')
            && text.bytes().all(|digit| digit.is_ascii_digit())
    };
    if cdna_start != cds_start || cdna_end != cds_end || !number(&cds_start) || !number(&cds_end) {
        return None;
    }
    let opens = if (cds_start.len(), cds_start.as_str()) <= (cds_end.len(), cds_end.as_str()) {
        cds_start.as_str()
    } else {
        cds_end.as_str()
    };
    matches!(opens, "1" | "2" | "3").then_some((cds_start, cds_end))
}

/// Whether the reference entry `e` and the candidate entry `a` of one key have the shape
/// of `start_codon_kept_without_five_prime_utr`: a Transcript entry; the candidate's terms
/// the reference's without `start_lost` and with `start_retained_variant`, which the
/// reference does not hold, every other term the same; one of `LENGTH_TERMS` on both
/// sides; and a cDNA span that is the CDS span and opens on a base of the start codon
/// (`span_on_the_start_codon`), the same on both sides.
pub fn start_codon_kept_without_five_prime_utr_pair(e: &Entry, a: &Entry) -> Result<bool> {
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
    if !actual.difference(&expected).eq([&START_RETAINED_TERM])
        || !expected
            .difference(&actual)
            .all(|term| *term == START_LOST_TERM)
        || actual.contains(START_LOST_TERM)
        || ![&expected, &actual]
            .iter()
            .all(|terms| LENGTH_TERMS.iter().any(|term| terms.contains(term)))
    {
        return Ok(false);
    }
    let span = span_on_the_start_codon(e);
    Ok(span.is_some() && span == span_on_the_start_codon(a))
}

/// Whether the differing values `ev` and `av` of the field `name` on a row
/// `start_codon_kept_without_five_prime_utr_pair` takes are the ones the class reads: a
/// protein notation `lost_start_against_the_change` takes on a row whose reference holds
/// `start_lost`, and the cells that follow the two sets of terms
/// (`term_pair_following_cell`).
pub fn start_codon_kept_without_five_prime_utr_cell(
    e: &Entry,
    a: &Entry,
    name: &str,
    ev: &str,
    av: &str,
) -> bool {
    if field_named(name, HGVSP_FIELDS) {
        return consequence_set(&e.consequence)
            .split(',')
            .any(|term| term == START_LOST_TERM)
            && lost_start_against_the_change(ev, av);
    }
    term_pair_following_cell(e, a, name, ev, av)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A Transcript entry in the spelling of the text formats: its terms and its two
    /// position columns.
    fn row(terms: &str, cdna: &str, cds: &str) -> Entry {
        let mut entry = Entry {
            key: "21:102|-|ENST1|Transcript".to_string(),
            allele: "-".to_string(),
            consequence: terms.to_string(),
            ..Entry::default()
        };
        for (name, value) in [("cDNA_position", cdna), ("CDS_position", cds)] {
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
    }

    fn pair(e: &Entry, a: &Entry) -> bool {
        start_codon_kept_without_five_prime_utr_pair(e, a).expect("the pair is read")
    }

    const LOST: &str = "frameshift_variant,start_lost";
    const KEPT: &str = "frameshift_variant,start_retained_variant";

    #[test]
    fn the_pair_is_a_row_on_the_start_codon_that_gains_the_kept_start_term() {
        for (reference, candidate, cdna, cds) in [
            // The reference holds start_lost: a frameshift and a deletion of whole codons.
            (LOST, KEPT, "3", "3"),
            (LOST, KEPT, "2-3", "2-3"),
            (LOST, KEPT, "1-250", "1-250"),
            (
                "start_lost,inframe_deletion",
                "inframe_deletion,start_retained_variant",
                "3-5",
                "3-5",
            ),
            // The reference holds neither start term.
            (
                "inframe_insertion",
                "inframe_insertion,start_retained_variant",
                "2-3",
                "2-3",
            ),
            (
                "inframe_deletion",
                "inframe_deletion,start_retained_variant",
                "1-3",
                "1-3",
            ),
            // A term both sides write stays, in either side's order.
            (
                "frameshift_variant,start_lost,NMD_transcript_variant",
                "frameshift_variant,NMD_transcript_variant,start_retained_variant",
                "2-3",
                "2-3",
            ),
            // Under --total_length each position carries the sequence's length.
            (LOST, KEPT, "3-5/1500", "3-5/900"),
        ] {
            let e = row(reference, cdna, cds);
            let a = row(candidate, cdna, cds);
            assert!(pair(&e, &a), "{reference} {cdna}");
            assert!(!pair(&a, &e), "{reference} {cdna}");
            assert!(!pair(&e, &e), "{reference} {cdna}");
            assert!(!pair(&a, &a), "{reference} {cdna}");
        }
        // The VCF format joins a row's terms with `&`.
        let e = row("frameshift_variant&start_lost", "3", "3");
        let a = row("frameshift_variant&start_retained_variant", "3", "3");
        assert!(pair(&e, &a));
        // A JSON entry: one key a bound, the start key below the end key or above it.
        for (start, end) in [("3", "3"), ("2", "3"), ("3", "2"), ("1", "250"), ("4", "3")] {
            let span = [Some(start), Some(end)];
            assert!(
                pair(&object(LOST, span, span), &object(KEPT, span, span)),
                "{start} {end}"
            );
        }
    }

    #[test]
    fn a_span_opens_on_one_of_the_first_three_coding_bases() {
        for (span, taken) in [
            ("1", true),
            ("2", true),
            ("3", true),
            ("4", false),
            ("1-2", true),
            ("3-4", true),
            ("3-900", true),
            ("4-5", false),
            ("4-6", false),
            ("10-12", false),
            ("30", false),
            ("13-15", false),
        ] {
            let e = row(LOST, span, span);
            let a = row(KEPT, span, span);
            assert_eq!(pair(&e, &a), taken, "{span}");
        }
        for (start, end, taken) in [
            ("4", "3", true),
            ("3", "4", true),
            ("5", "4", false),
            ("4", "5", false),
            ("10", "3", true),
            ("10", "9", false),
        ] {
            let span = [Some(start), Some(end)];
            assert_eq!(
                pair(&object(LOST, span, span), &object(KEPT, span, span)),
                taken,
                "{start} {end}"
            );
        }
    }

    #[test]
    fn a_row_that_breaks_one_condition_is_refused() {
        let e = row(LOST, "3", "3");
        let a = row(KEPT, "3", "3");
        assert!(pair(&e, &a));
        // A feature that is no transcript.
        let regulatory = |entry: &Entry| Entry {
            key: entry.key.replace("Transcript", "RegulatoryFeature"),
            ..entry.clone()
        };
        assert!(!pair(&regulatory(&e), &regulatory(&a)));
        // Terms that differ otherwise than by the two start terms.
        for (reference, candidate) in [
            (
                LOST,
                "frameshift_variant,start_retained_variant,splice_region_variant",
            ),
            ("frameshift_variant,start_lost,splice_region_variant", KEPT),
            ("frameshift_variant,stop_lost", KEPT),
            (LOST, "frameshift_variant,stop_retained_variant"),
            (LOST, "frameshift_variant"),
            (LOST, LOST),
            (KEPT, KEPT),
            // The candidate keeps start_lost.
            (LOST, "frameshift_variant,start_lost,start_retained_variant"),
            // The reference holds start_retained_variant.
            ("frameshift_variant,start_lost,start_retained_variant", KEPT),
            // A frameshift on one side alone.
            (LOST, "inframe_deletion,start_retained_variant"),
        ] {
            assert!(
                !pair(&row(reference, "3", "3"), &row(candidate, "3", "3")),
                "{reference} against {candidate}"
            );
        }
        // No frameshift_variant, inframe_insertion or inframe_deletion.
        for other in [
            "",
            "coding_sequence_variant",
            "protein_altering_variant",
            "missense_variant",
            "5_prime_UTR_variant",
        ] {
            let lost = format!("start_lost,{other}");
            let kept = format!("{other},start_retained_variant");
            assert!(
                !pair(&row(&lost, "3", "3"), &row(&kept, "3", "3")),
                "{other}"
            );
        }
        // A cDNA span that is not the CDS span, a bound that is no number, and positions
        // that differ between the sides.
        for (cdna, cds) in [
            ("103", "3"),
            ("4", "3"),
            ("3", "4"),
            ("103-105", "3-5"),
            ("3-6", "3-5"),
            ("2-5", "3-5"),
            ("?-5", "?-5"),
            ("3-?", "3-?"),
            ("?", "?"),
            ("-", "3"),
            ("3", "-"),
            ("-", "-"),
            ("", ""),
            ("0-2", "0-2"),
            ("03", "03"),
            ("3-05", "3-05"),
            ("x", "x"),
        ] {
            assert!(
                !pair(&row(LOST, cdna, cds), &row(KEPT, cdna, cds)),
                "{cdna} {cds}"
            );
        }
        assert!(!pair(&e, &row(KEPT, "2", "2")));
        assert!(!pair(&e, &row(KEPT, "3-4", "3-4")));
        let whole = [Some("3"), Some("5")];
        for (cdna, cds) in [
            ([Some("103"), Some("105")], whole),
            (whole, [Some("3"), Some("6")]),
            ([None, Some("5")], [None, Some("5")]),
            ([Some("3"), None], [Some("3"), None]),
            ([None, None], [None, None]),
        ] {
            assert!(
                !pair(&object(LOST, cdna, cds), &object(KEPT, cdna, cds)),
                "{cdna:?} {cds:?}"
            );
        }
        assert!(!pair(
            &object(LOST, whole, whole),
            &object(KEPT, [Some("3"), Some("6")], [Some("3"), Some("6")])
        ));
    }

    #[test]
    fn a_lost_start_notation_stands_against_the_notation_of_the_change() {
        let taken = |reference: &str, candidate: &str| {
            lost_start_against_the_change(
                &format!("ENSP1.1:p.{reference}"),
                &format!("ENSP1.1:p.{candidate}"),
            )
        };
        for (reference, candidate) in [
            ("Ala2?", "Ala2LeufsTer37"),
            ("Ser3_?2", "Ser3AlafsTer2"),
            ("Ala2_?84", "Ala2LysfsTer6"),
            ("Leu2?", "Leu2del"),
            ("Ala3?", "Ala3del"),
            ("(Ala3?)", "(Ala3del)"),
        ] {
            assert!(taken(reference, candidate), "{reference} {candidate}");
            assert!(!taken(candidate, reference), "{reference} {candidate}");
            assert!(!taken(reference, reference), "{reference}");
            assert!(!taken(candidate, candidate), "{candidate}");
        }
        // A candidate that holds a question mark, and a notation without a change.
        assert!(!taken("Ala2?", "Ala2LeufsTer?"));
        assert!(!taken("Ala2?", "?"));
        assert!(!taken("Ala2?", ""));
        assert!(!taken("", "Ala2del"));
        // Two proteins, a notation without a protein, and no notation.
        assert!(!lost_start_against_the_change(
            "ENSP1.1:p.Ala3?",
            "ENSP2.1:p.Ala3del"
        ));
        assert!(!lost_start_against_the_change("p.Ala3?", "p.Ala3del"));
        assert!(!lost_start_against_the_change(":p.Ala3?", ":p.Ala3del"));
        assert!(!lost_start_against_the_change("Ala3?", "Ala3del"));
        assert!(!lost_start_against_the_change("-", "ENSP1.1:p.Ala3del"));
        assert!(!lost_start_against_the_change("ENSP1.1:p.Ala3?", "-"));
        assert!(!lost_start_against_the_change("", ""));
        // A question mark in the protein's name is none of the change.
        assert!(!lost_start_against_the_change(
            "ENSP?.1:p.Ala3Val",
            "ENSP?.1:p.Ala3del"
        ));
    }

    #[test]
    fn the_class_reads_the_impact_of_each_sides_terms_and_the_notation_beside_start_lost() {
        let e = row("start_lost,inframe_deletion", "3-5", "3-5");
        let a = row("inframe_deletion,start_retained_variant", "3-5", "3-5");
        let cell = |name: &str, ev: &str, av: &str| {
            start_codon_kept_without_five_prime_utr_cell(&e, &a, name, ev, av)
        };
        assert!(cell("HGVSp", "ENSP1.1:p.Ala3?", "ENSP1.1:p.Ala3del"));
        assert!(cell("hgvsp", "ENSP1.1:p.Ala3?", "ENSP1.1:p.Ala3del"));
        assert!(!cell("HGVSp", "ENSP1.1:p.Ala3del", "ENSP1.1:p.Ala3?"));
        assert!(!cell("HGVSp", "ENSP1.1:p.Ala3?", "-"));
        // Utils/Constants.pm gives start_lost the impact HIGH and inframe_deletion MODERATE.
        assert!(cell("IMPACT", "HIGH", "MODERATE"));
        assert!(cell("impact", "HIGH", "MODERATE"));
        assert!(!cell("IMPACT", "MODERATE", "HIGH"));
        assert!(!cell("IMPACT", "HIGH", "LOW"));
        assert!(!cell("IMPACT", "LOW", "MODERATE"));
        assert!(!cell("Codons", "atGGCg/atg", "atGGCg/atc"));
        assert!(!cell("STRAND", "1", "-1"));
        // Both sets of a frameshift's row give HIGH, so no IMPACT follows the pair.
        let e = row(LOST, "3", "3");
        let a = row(KEPT, "3", "3");
        assert!(start_codon_kept_without_five_prime_utr_cell(
            &e,
            &a,
            "HGVSp",
            "ENSP1.1:p.Ala2?",
            "ENSP1.1:p.Ala2LeufsTer37"
        ));
        assert!(!start_codon_kept_without_five_prime_utr_cell(
            &e, &a, "IMPACT", "HIGH", "LOW"
        ));
        // A reference without start_lost writes no lost-start notation the class reads.
        let e = row("inframe_insertion", "2-3", "2-3");
        let a = row("inframe_insertion,start_retained_variant", "2-3", "2-3");
        assert!(!start_codon_kept_without_five_prime_utr_cell(
            &e,
            &a,
            "HGVSp",
            "ENSP1.1:p.Met1?",
            "ENSP1.1:p.Met1_Ala2insGlu"
        ));
        assert!(!start_codon_kept_without_five_prime_utr_cell(
            &e, &a, "IMPACT", "MODERATE", "LOW"
        ));
    }
}
