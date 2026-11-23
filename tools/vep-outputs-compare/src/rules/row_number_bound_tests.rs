// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The numbers a row-aware pair rule reads of an output row are read within 128 bits: the
//! position a generated name spells, the two ends of a Location and the two cDNA positions.
//! A number past that bound reads as no number, and every rule that needs the number refuses
//! the pair, whose two rows stay differences. Each pair below is taken where its numbers lie
//! within the bound and refused by every row-aware rule where one of them lies past it. Every
//! row is fabricated.

use super::super::row_pair_rules;
use super::*;

/// The last value of 128 bits and the value after it.
const LAST_WITHIN: &str = "170141183460469231731687303715884105727";
const FIRST_PAST: &str = "170141183460469231731687303715884105728";

/// A number past 128 bits that ends in `tail`: a one before forty digits.
fn past(tail: u32) -> String {
    format!("1{tail:040}")
}

/// One Transcript row under the name `name`.
fn row(name: &str, location: &str, allele: &str, terms: &str, cdna: &str, cds: &str) -> PairRow {
    pair_row(&[
        name,
        location,
        allele,
        "ENSG1",
        "ENST90000000001",
        "Transcript",
        terms,
        cdna,
        cds,
    ])
}

/// Whether any row-aware pair rule takes the reference's row `e` and vep-rs's row `a`.
fn taken(e: &PairRow, a: &PairRow) -> bool {
    row_pair_rules()
        .iter()
        .any(|rule| rule.matches_rows(std::slice::from_ref(e), std::slice::from_ref(a)))
}

#[test]
fn each_number_is_read_up_to_the_last_value_of_128_bits_and_no_further() {
    assert_eq!(LAST_WITHIN, i128::MAX.to_string());
    let name = |position: &str| format!("7_{position}_TA/GTGTG");
    assert_eq!(uploaded_position(&name(LAST_WITHIN)), Some(i128::MAX));
    assert_eq!(uploaded_position(&name(FIRST_PAST)), None);
    assert_eq!(
        location_ends(&format!("7:5-{LAST_WITHIN}")),
        Some((5, i128::MAX))
    );
    assert_eq!(location_ends(&format!("7:5-{FIRST_PAST}")), None);
    assert_eq!(location_ends(&format!("7:{FIRST_PAST}")), None);
    assert_eq!(
        positions_named(&format!("7:{}-{LAST_WITHIN}", i128::MAX - 2)),
        Some(3)
    );
    assert_eq!(
        positions_named(&format!("7:{LAST_WITHIN}-{FIRST_PAST}")),
        None
    );
    assert_eq!(
        location_span(&format!("7:{}-{LAST_WITHIN}", i128::MAX - 1)),
        Some(2)
    );
    assert_eq!(
        location_span(&format!("7:{LAST_WITHIN}-{FIRST_PAST}")),
        None
    );
    assert_eq!(location_span(&format!("7:{}-{}", past(46), past(47))), None);
    let before_last = (i128::MAX - 2).to_string();
    assert_eq!(
        sequence_width(&[&before_last, LAST_WITHIN, "-", "-"]),
        Some(3)
    );
    assert_eq!(sequence_width(&[LAST_WITHIN, FIRST_PAST, "-", "-"]), None);
    assert_eq!(sequence_width(&[&past(998), &past(1000), "-", "-"]), None);
}

/// `location_span` through `insertion_row` and `inserted_between_neighbours`: two
/// neighbouring positions under an identifier.
#[test]
fn an_insertions_pair_at_a_location_past_the_bound_stays_a_difference() {
    let pair = |location: &str| {
        let row = |term: &str| {
            let terms = format!("inframe_insertion,{term}");
            row(
                "1993671",
                location,
                "TTGA",
                &terms,
                "7054-7055",
                "6977-6978",
            )
        };
        (row("coding_sequence_variant"), row("stop_retained_variant"))
    };
    let (e, a) = pair("19:13208046-13208047");
    assert!(stop_codon_overlap_insertion_range_pair(&e, &a));
    let (e, a) = pair(&format!("19:{}-{}", past(13208046), past(13208047)));
    assert!(!insertion_row(&e) && !inserted_between_neighbours(&e));
    assert!(!taken(&e, &a));
}

/// `location_span` through `indel_row` and `length_changing_row`: an allele of one letter at
/// two positions under an identifier.
#[test]
fn an_indels_pair_at_a_location_past_the_bound_stays_a_difference() {
    let pair = |location: &str| {
        let row = |term: &str| {
            let terms = format!("splice_region_variant,{term}");
            row("1993671", location, "A", &terms, "1704-1705", "1629-1630")
        };
        (row("coding_sequence_variant"), row("stop_retained_variant"))
    };
    let (e, a) = pair("3:36993662-36993663");
    assert!(stop_codon_overlap_genomic_span_pair(&e, &a));
    let (e, a) = pair(&format!("3:{}-{}", past(36993662), past(36993663)));
    assert!(!indel_row(&e) && !length_changing_row(&e));
    assert!(!taken(&e, &a));
}

/// The reference's row and vep-rs's of a longer allele whose donor term the reference alone
/// writes, at `location` under `name` with the cDNA positions `cdna`.
fn past_the_span_pair(
    name: &str,
    location: &str,
    reference_terms: &str,
    cdna: &str,
) -> (PairRow, PairRow) {
    let row = |terms: &str| row(name, location, "GTGTG", terms, cdna, "-");
    (
        row(reference_terms),
        row("frameshift_variant,splice_region_variant"),
    )
}

/// `location_ends` through `positions_named`, and `sequence_width`, where no record is read:
/// three positions in the transcript's sequence under an identifier.
#[test]
fn a_longer_alleles_pair_without_its_record_stays_a_difference_past_the_bound() {
    let terms = "splice_donor_variant,frameshift_variant";
    let pair = |location: &str, cdna: &str| past_the_span_pair("1993671", location, terms, cdna);
    let (e, a) = pair("7:100998-101000", "998-1000");
    assert!(splice_term_past_the_record_span_pair(&e, &a));
    let far_location = format!("7:{}-{}", past(100998), past(101000));
    let (e, a) = pair(&far_location, "998-1000");
    assert_eq!(positions_named(&e.location), None);
    assert!(!taken(&e, &a));
    let far_positions = format!("{}-{}", past(998), past(1000));
    let (e, a) = pair("7:100998-101000", &far_positions);
    assert_eq!(
        transcript_positions(&e).map(|positions| sequence_width(&positions)),
        Some(None)
    );
    assert!(!taken(&e, &a));
}

const RECORD_TERMS: &str =
    "splice_donor_variant,frameshift_variant,splice_donor_region_variant,intron_variant";

/// `sequence_width` through `intron_placements`, where the row's name spells its record: the
/// two replaced bases of a generated name.
#[test]
fn a_longer_alleles_pair_beside_its_record_stays_a_difference_past_the_bound() {
    let name = "7_100999_TA/GTGTG";
    let pair = |cdna: &str| past_the_span_pair(name, "7:100999-101000", RECORD_TERMS, cdna);
    let (e, a) = pair("999-1000");
    assert_eq!(e.reference.as_deref(), Some("TA"));
    assert!(splice_term_past_the_record_span_pair(&e, &a));
    let (e, a) = pair(&format!("{}-{}", past(999), past(1000)));
    assert_eq!(e.reference.as_deref(), Some("TA"));
    assert!(!taken(&e, &a));
}

/// `uploaded_position` and `location_ends` through `pair_row`: a generated name whose position
/// lies past the bound gives its row no reference bases, as an identifier gives none, and a
/// row of two neighbouring positions without them is refused.
#[test]
fn a_generated_name_past_the_bound_gives_its_row_no_reference_bases() {
    let location = "7:100999-101000";
    let pair =
        |name: &str, location: &str| past_the_span_pair(name, location, RECORD_TERMS, "999-1000");
    let (e, a) = pair("7_100999_TA/GTGTG", location);
    assert_eq!(e.reference.as_deref(), Some("TA"));
    assert!(taken(&e, &a));
    let far_name = format!("7_{}_TA/GTGTG", past(100999));
    let (e, a) = pair(&far_name, location);
    assert_eq!((e.reference.as_ref(), e.indel), (None, Some(true)));
    assert!(!taken(&e, &a));
    let far_location = format!("7:{}-{}", past(100999), past(101000));
    let (e, a) = pair(&far_name, &far_location);
    assert_eq!(e.reference, None, "a record that lies past the bound whole");
    assert!(!taken(&e, &a));
}
