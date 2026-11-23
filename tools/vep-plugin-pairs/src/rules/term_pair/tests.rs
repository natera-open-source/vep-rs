// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pair rules on rows the default-format parser builds: the expected values are those
//! `tools/vep-outputs-compare` gives on the same rows, and for the rules of
//! `stop_codon_settled_before_unknown_bases`, `start_codon_removed_without_five_prime_utr`
//! and `start_codon_kept_without_five_prime_utr` those `tools/vep-formats-compare` gives.
//! The rule of `frameshift_inside_a_held_stop_codon` is held to the five made-up rows of
//! its class and their two controls.

use super::*;
use crate::rules::entry::parse_default;
use std::cell::RefCell;

const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

/// The entry of a row's first nine columns, Uploaded_variation to CDS_position.
fn pair_row(columns: &[&str; 9]) -> Entry {
    let text = format!("{HEADER}{}\t-\t-\t-\t-\t-\n", columns.join("\t"));
    let rows = parse_default(&text).expect("a row of fourteen columns");
    Entry::clone(&rows[0])
}

/// A row of one transcript under an identifier: the Location, the allele, the terms and
/// the two position columns.
fn row(location: &str, allele: &str, consequence: &str, cdna: &str, cds: &str) -> Entry {
    pair_row(&[
        "1993671",
        location,
        allele,
        "ENSG1",
        "ENST00000635895",
        "Transcript",
        consequence,
        cdna,
        cds,
    ])
}

/// A row of one transcript whose allele stands in for the reference bases `reference`.
fn record_row(reference: &str, location: &str, allele: &str, terms: &str, cdna: &str) -> Entry {
    Entry {
        reference: RefCell::new(Some(reference.to_string())),
        ..row(location, allele, terms, cdna, "-")
    }
}

fn stop_codon_overlap_genomic_span_pair(e: &Entry, a: &Entry) -> bool {
    super::stop_codon_overlap_genomic_span_pair(e, a).expect("no refusal")
}

fn stop_retained_without_stop_codon_test_pair(e: &Entry, a: &Entry) -> bool {
    super::stop_retained_without_stop_codon_test_pair(e, a).expect("no refusal")
}

fn stop_retained_three_base_stop_codon_test_pair(e: &Entry, a: &Entry) -> bool {
    super::stop_retained_three_base_stop_codon_test_pair(e, a).expect("no refusal")
}

fn stop_codon_removed_without_stop_codon_test_pair(e: &Entry, a: &Entry) -> bool {
    super::stop_codon_removed_without_stop_codon_test_pair(e, a).expect("no refusal")
}

fn stop_codon_overlap_insertion_range_pair(e: &Entry, a: &Entry) -> bool {
    super::stop_codon_overlap_insertion_range_pair(e, a).expect("no refusal")
}

fn splice_term_past_the_record_span_pair(e: &Entry, a: &Entry) -> bool {
    splice_term_past_the_record_span_rows(e, a).expect("no refusal")
}

/// The bounds of a cDNA_position column.
fn position_bounds(position: &str) -> (bool, bool) {
    let e = row("3:36993663", "-", "", position, "-");
    super::position_bounds(&e, "cDNA_position", "cdna_start", "cdna_end")
}

/// The start and the end a cDNA_position column writes.
fn position_parts(position: &str) -> Option<(String, String)> {
    let e = row("3:36993663", "-", "", position, "-");
    transcript_positions(&e).map(|parts| (parts[0].to_string(), parts[1].to_string()))
}

/// The places of an intron beside a span that two cDNA positions allow, in their order.
fn places(span: i64, cdna: (&str, &str)) -> Vec<(bool, i64, i64)> {
    let mut places = Vec::new();
    let visit = |forward, start, end| {
        places.push((forward, start, end));
        false
    };
    intron_placements(span, 2, &[cdna.0, cdna.1, "-", "-"], visit).expect("no refusal");
    places
}

#[test]
fn a_position_column_carries_its_start_and_its_end() {
    for (position, bounds) in [
        ("1705", (true, true)),
        ("1704-1705", (true, true)),
        ("?-7057", (false, true)),
        ("1705-?", (true, false)),
        ("-", (false, false)),
        ("", (false, false)),
        ("-/7140", (false, false)),
        ("120-7100/7140", (true, true)),
        ("?-6980/7023", (false, true)),
    ] {
        assert_eq!(position_bounds(position), bounds, "{position}");
    }
}

#[test]
fn the_terms_of_a_consequence_are_a_set() {
    assert!(is(&terms("b,a&b,,"), &["a", "b"]));
    assert!(is(&terms(""), &[""]));
    let (expected, gone, come) = (terms("a,b"), terms("a,b"), terms("a,c"));
    assert!(expected == gone && gone != come);
}

#[test]
fn the_genomic_span_class_reads_a_tested_row_at_an_exon_end() {
    let terms = |term: &str| format!("splice_region_variant,{term},NMD_transcript_variant");
    let pair = |e: &str, a: &str, cdna: &str, cds: &str| {
        stop_codon_overlap_genomic_span_pair(
            &row("3:36993663", "-", &terms(e), cdna, cds),
            &row("3:36993663", "-", &terms(a), cdna, cds),
        )
    };
    let (coding, retained) = ("coding_sequence_variant", "stop_retained_variant");
    assert!(pair(coding, retained, "1705", "1630"));
    assert!(pair(coding, retained, "1705/1800", "1630/1653"));
    assert!(!pair(retained, coding, "1705", "1630"));
    assert!(!pair(coding, retained, "1705-?", "1630-?"));
    assert!(!pair(coding, retained, "1705", "-"));
    assert!(!pair(
        coding,
        "stop_retained_variant,stop_lost",
        "1705",
        "1630"
    ));
    assert!(!stop_codon_overlap_genomic_span_pair(
        &row("3:36993663", "-", coding, "1705", "1630"),
        &row("3:36993663", "-", retained, "1705", "1630"),
    ));
    assert!(!stop_codon_overlap_genomic_span_pair(
        &row("3:36993663", "A", &terms(coding), "1705", "1630"),
        &row("3:36993663", "A", &terms(retained), "1705", "1630"),
    ));
}

#[test]
fn the_untested_class_reads_rows_whose_stop_codon_test_is_not_made() {
    let pair = stop_retained_without_stop_codon_test_pair;
    let acceptor = |term: &str| format!("splice_acceptor_variant,{term}");
    let (location, retained) = ("19:13208046-13208050", "stop_retained_variant");
    for (cdna, cds) in [
        ("?-7057", "-"),
        ("-", "-"),
        ("1705-?", "1630-?"),
        ("120-7100", "-"),
    ] {
        let reference = row(location, "-", &acceptor(retained), cdna, cds);
        for term in ["coding_sequence_variant", "frameshift_variant"] {
            assert!(pair(
                &reference,
                &row(location, "-", &acceptor(term), cdna, cds)
            ));
        }
        let nothing = row(location, "-", "splice_acceptor_variant", cdna, cds);
        assert!(pair(&reference, &nothing), "{cdna} {cds}");
    }
    let lettered = |term: &str| {
        row(
            "7:152648662-152648663",
            "ACRT",
            term,
            "7038-7039",
            "6961-6962",
        )
    };
    assert!(pair(&lettered(retained), &lettered("frameshift_variant")));
    let tested = |term: &str| row(location, "-", &acceptor(term), "7055-7059", "6978-6982");
    assert!(!pair(&tested(retained), &tested("coding_sequence_variant")));
    let untested = |term: &str| row(location, "-", &acceptor(term), "?-7057", "-");
    assert!(!pair(&untested(retained), &untested("stop_lost")));
    assert!(!pair(
        &untested(retained),
        &untested("coding_sequence_variant,frameshift_variant")
    ));
    assert!(!pair(
        &untested("coding_sequence_variant"),
        &untested(retained)
    ));
    assert!(!pair(
        &untested(retained),
        &tested("coding_sequence_variant")
    ));
    let substitution = |term: &str| row("19:13208046", "A", term, "-", "-");
    assert!(!pair(
        &substitution(retained),
        &substitution("coding_sequence_variant")
    ));
}

/// On a row whose first column spells a reference allele and two alternates, none of them
/// the row's allele once the shared first base is dropped, the Location decides whether the
/// row changes the length, as under an identifier; an alternate that is the row's allele
/// decides before it.
#[test]
fn the_three_base_class_reads_an_allele_no_alternate_names_by_its_location() {
    let taken = |location: &str, allele: &str| {
        let row = |consequence: &str| {
            pair_row(&[
                "19_13208044_ACGTT/ACA/AG",
                location,
                allele,
                "ENSG1",
                "ENST00000635895",
                "Transcript",
                consequence,
                "7054-7056",
                "6977-6979",
            ])
        };
        stop_retained_three_base_stop_codon_test_pair(
            &row("stop_retained_variant"),
            &row("coding_sequence_variant"),
        )
    };
    assert!(taken("19:13208046-13208048", "A"));
    assert!(!taken("19:13208046", "A"));
    assert!(taken("19:13208045", "G"));
}

#[test]
fn the_three_base_class_reads_rows_whose_stop_codon_test_is_made() {
    let pair = stop_retained_three_base_stop_codon_test_pair;
    let location = "2:190970701-190970707";
    let utr = |term: &str, cdna: &str| {
        row(
            location,
            "-",
            &format!("{term},3_prime_UTR_variant"),
            cdna,
            "2249-?",
        )
    };
    let (retained, coding) = ("stop_retained_variant", "coding_sequence_variant");
    assert!(pair(&utr(retained, "2558-2564"), &utr(coding, "2558-2564")));
    let coded = |term: &str| row(location, "-", term, "961-967", "264-270");
    assert!(pair(&coded(retained), &coded("frameshift_variant")));
    assert!(pair(
        &coded("inframe_deletion,stop_retained_variant"),
        &coded("frameshift_variant")
    ));
    assert!(pair(
        &coded("protein_altering_variant,stop_retained_variant"),
        &coded("stop_gained,frameshift_variant")
    ));
    assert!(!pair(&utr(retained, "2558-?"), &utr(coding, "2558-?")));
    assert!(!pair(&utr(retained, "2558-2564"), &utr(coding, "2558-?")));
    assert!(!pair(
        &utr(coding, "2558-2564"),
        &utr(retained, "2558-2564")
    ));
    assert!(!pair(
        &coded("inframe_deletion,stop_retained_variant"),
        &coded(coding)
    ));
    assert!(!pair(
        &coded("missense_variant,stop_retained_variant"),
        &coded("frameshift_variant")
    ));
    assert!(!pair(
        &utr(retained, "2558-2564"),
        &utr("stop_lost", "2558-2564")
    ));
    let lettered = |term: &str| {
        row(
            "7:152648662-152648663",
            "ACGTACGN",
            term,
            "1213-1214",
            "1100-1101",
        )
    };
    assert!(!pair(&lettered(retained), &lettered("frameshift_variant")));
}

#[test]
fn the_removed_codon_class_reads_rows_without_the_tests_coordinates() {
    let pair = stop_codon_removed_without_stop_codon_test_pair;
    let location = "19:13208046-13208050";
    let acceptor = |term: &str, cdna: &str, cds: &str| {
        row(
            location,
            "-",
            &format!("splice_acceptor_variant,{term}"),
            cdna,
            cds,
        )
    };
    let (retained, lost) = ("stop_retained_variant", "stop_lost");
    for (cdna, cds) in [
        ("?-7057", "-"),
        ("-", "-"),
        ("1705-?", "1630-?"),
        ("120-7100", "-"),
        ("21-99", "?-69"),
    ] {
        assert!(
            pair(&acceptor(retained, cdna, cds), &acceptor(lost, cdna, cds)),
            "{cdna} {cds}"
        );
    }
    for (cdna, cds) in [("7055-7059", "6978-6982"), ("7055-7059", "6978-?")] {
        assert!(
            !pair(&acceptor(retained, cdna, cds), &acceptor(lost, cdna, cds)),
            "{cdna} {cds}"
        );
    }
    assert!(!pair(
        &acceptor(retained, "?-7057", "-"),
        &acceptor("stop_lost,frameshift_variant", "?-7057", "-")
    ));
    assert!(!pair(
        &acceptor(lost, "?-7057", "-"),
        &acceptor(retained, "?-7057", "-")
    ));
    assert!(!pair(
        &acceptor(retained, "?-7057", "-"),
        &acceptor("coding_sequence_variant", "?-7057", "-")
    ));
    let substitution = |term: &str| row("19:13208046", "A", term, "-", "-");
    assert!(!pair(&substitution(retained), &substitution(lost)));
}

#[test]
fn the_insertion_range_class_reads_an_insertion_that_keeps_the_stop_codon() {
    let pair = stop_codon_overlap_insertion_range_pair;
    let inserted = |allele: &str, terms: &str, cdna: &str, cds: &str| {
        row("19:13208046-13208047", allele, terms, cdna, cds)
    };
    let whole = |terms: &str| inserted("TTGA", terms, "7054-7055", "6977-6978");
    for (reference, candidate) in [
        (
            "inframe_insertion,coding_sequence_variant",
            "inframe_insertion,stop_retained_variant",
        ),
        (
            "stop_gained,frameshift_variant",
            "inframe_insertion,stop_retained_variant",
        ),
        (
            "frameshift_variant",
            "inframe_insertion,stop_retained_variant",
        ),
        ("stop_gained,frameshift_variant", "stop_retained_variant"),
    ] {
        assert!(pair(&whole(reference), &whole(candidate)), "{reference}");
        assert!(!pair(&whole(candidate), &whole(reference)), "{reference}");
        assert!(pair(
            &inserted("TNGA", reference, "7054-7055", "6978-?"),
            &inserted("TNGA", candidate, "7054-7055", "6978-?")
        ));
        assert!(!pair(
            &inserted("TTGA", reference, "7054-?", "6977-?"),
            &inserted("TTGA", candidate, "7054-?", "6977-?")
        ));
        assert!(!pair(
            &inserted("TRGA", reference, "7054-7055", "6977-6978"),
            &inserted("TRGA", candidate, "7054-7055", "6977-6978")
        ));
        let deleted =
            |terms: &str| row("19:13208046-13208049", "-", terms, "7055-7058", "6978-6981");
        assert!(!pair(&deleted(reference), &deleted(candidate)));
    }
    let region = |term: &str| whole(&format!("inframe_insertion,splice_region_variant,{term}"));
    assert!(!pair(
        &region("coding_sequence_variant"),
        &region("stop_retained_variant")
    ));
    assert!(stop_codon_overlap_genomic_span_pair(
        &region("coding_sequence_variant"),
        &region("stop_retained_variant")
    ));
    assert!(!pair(
        &whole("frameshift_variant"),
        &whole("protein_altering_variant,stop_retained_variant")
    ));
    assert!(!pair(
        &whole("stop_gained"),
        &whole("stop_retained_variant")
    ));
}

#[test]
fn the_insertion_range_class_takes_an_insertion_between_two_positions_alone() {
    // `read` is the reference bases the lines of the input give a row that carries none.
    let taken_beside = |name: &str, location: &str, allele: &str, read: Option<&str>| {
        let row = |terms: &str| {
            let row = pair_row(&[
                name,
                location,
                allele,
                "ENSG1",
                "ENST00000635895",
                "Transcript",
                terms,
                "7054-7055",
                "6977-6978",
            ]);
            if row.reference.borrow().is_none() {
                *row.reference.borrow_mut() = read.map(str::to_string);
            }
            row
        };
        stop_codon_overlap_insertion_range_pair(
            &row("coding_sequence_variant"),
            &row("stop_retained_variant"),
        )
    };
    let taken =
        |name: &str, location: &str, allele: &str| taken_beside(name, location, allele, None);
    let two = "19:13208046-13208047";
    for (name, location, allele, shape) in [
        (
            "19_13208046_AC/GTT",
            two,
            "GTT",
            "two bases replaced by three",
        ),
        (
            "19_13208046_AC/AGT",
            two,
            "AGT",
            "the same, the first base shared",
        ),
        (
            "19_13208046_AC/A/T",
            two,
            "A",
            "a deletion written as a sequence",
        ),
        (
            "19_13208046_A/AC/T",
            "19:13208046",
            "AC",
            "an insertion kept with its first base",
        ),
        (
            "19_13208046_-/TGA",
            "19:13208046",
            "TGA",
            "a spelt insertion at one position",
        ),
        (
            "19_13208046_AC/A/T",
            two,
            "GTT",
            "an allele no alternate names, the reference of two",
        ),
        (
            "1993671",
            two,
            "CG",
            "a two-base allele under an identifier",
        ),
        ("1993671", "19:13208046-13208048", "TGA", "three positions"),
    ] {
        assert!(!taken(name, location, allele), "{shape}");
    }
    for (name, allele, shape) in [
        (
            "19_13208046_-/TGA",
            "TGA",
            "an insertion, the reference allele empty",
        ),
        (
            "19_13208046_A/ATGA",
            "TGA",
            "an insertion spelt with its shared first base",
        ),
        (
            "19_13208046_-/CG/T",
            "CG",
            "a two-base insertion, the alleles spelt",
        ),
        ("1993671", "TGA", "an insertion under an identifier"),
        (
            "1993671",
            "G",
            "under an identifier, two bases replaced by one",
        ),
        (
            "1993671",
            "CGT",
            "under an identifier, two bases replaced by three",
        ),
        (
            "19_13208045_TA/TATGA",
            "TGA",
            "an insertion after the second base its record keeps",
        ),
        (
            "19_13208046_AC/ATGAC",
            "TGA",
            "an insertion between the two bases its record keeps",
        ),
        (
            "19_13208046_A/ATGA/C",
            "TGA",
            "an insertion beside an alternate of another first base",
        ),
    ] {
        assert!(taken(name, two, allele), "{shape}");
    }
    for (name, allele, read, shape) in [
        ("1993671", "CG", "", "a two-base insertion"),
        ("1993671", "TGA", "", "an insertion"),
        (
            "19_13208046_-/TGA",
            "TGA",
            "AC",
            "replaced bases beside a first column that spells an insertion",
        ),
    ] {
        assert!(
            taken_beside(name, two, allele, Some(read)),
            "the lines of the input give {shape}"
        );
    }
    for (allele, shape) in [
        ("CGT", "two bases replaced by three"),
        ("G", "two bases replaced by one"),
        ("CG", "two bases replaced by two"),
    ] {
        assert!(
            !taken_beside("1993671", two, allele, Some("AC")),
            "under an identifier, the lines of the input give {shape}"
        );
    }
    assert!(
        !taken_beside("1993671", "19:13208046-13208048", "TGA", Some("")),
        "three positions, whatever the bases read"
    );
}

#[test]
fn the_three_base_class_takes_the_retained_term_alone_or_beside_one_read_away_term() {
    let pair = stop_retained_three_base_stop_codon_test_pair;
    let coded = |terms: &str| row("2:190970701-190970707", "-", terms, "961-967", "264-270");
    let retained = "stop_retained_variant";
    for candidate in ["frameshift_variant", "stop_gained,frameshift_variant"] {
        assert!(pair(&coded(retained), &coded(candidate)), "alone");
        for one in FRAMESHIFT_READ_TERMS {
            let beside = format!("{one},{retained}");
            assert!(pair(&coded(&beside), &coded(candidate)), "{one}");
        }
        for several in [
            "inframe_insertion,inframe_deletion",
            "inframe_insertion,protein_altering_variant",
            "inframe_deletion,protein_altering_variant",
            "inframe_insertion,inframe_deletion,protein_altering_variant",
        ] {
            let beside = format!("{several},{retained}");
            assert!(!pair(&coded(&beside), &coded(candidate)), "{several}");
        }
    }
}

#[test]
fn a_position_column_writes_its_start_and_its_end() {
    for (position, parts) in [
        ("999-1000", Some(("999", "1000"))),
        ("999-1000/2000", Some(("999", "1000"))),
        ("1705", Some(("1705", "1705"))),
        ("1705/1800", Some(("1705", "1705"))),
        ("1000-?", Some(("1000", "?"))),
        ("?-1000", Some(("?", "1000"))),
        ("?", Some(("?", "?"))),
        ("-", Some(("-", "-"))),
        ("", Some(("-", "-"))),
        ("-/7140", None),
        ("999-", None),
        ("-1000", None),
        ("999-1000-1001", None),
        ("9x9-1000", None),
        ("<absent>", None),
    ] {
        let written = position_parts(position);
        let written = written
            .as_ref()
            .map(|(start, end)| (start.as_str(), end.as_str()));
        assert_eq!(written, parts, "{position}");
    }
    assert_eq!(
        transcript_positions(&row("7:100998-101000", "GTGTG", "", "998-1000/2000", "-")),
        Some(["998", "1000", "-", "-"])
    );
    assert_eq!(
        transcript_positions(&row("7:100998-101000", "GTGTG", "", "998-1000", "9x9")),
        None
    );
}

const TWO_POSITIONS: &str = "7:100999-101000";

fn letters(sequence: &str) -> Vec<char> {
    sequence.chars().collect()
}

#[test]
fn a_record_is_cut_into_the_runs_of_indices_at_which_its_alleles_differ() {
    let regions =
        |reference: &str, allele: &str| differing_regions(&letters(reference), &letters(allele));
    assert_eq!(regions("TA", "GAGT"), [(0, 0), (2, 3)]);
    assert_eq!(regions("TA", "GTGTG"), [(0, 4)]);
    assert_eq!(regions("TA", "TA"), [(0, 1)], "no index differs");
    assert_eq!(regions("T", "GTG"), [(0, 0)], "a reference of one base");
    assert_eq!(regions("TAC", "G"), [(0, 2)], "an allele of one base");
}

/// Whether the class takes the reference's and the candidate's row of one transcript, the
/// Location, the allele and the two position columns the same on both: the reference
/// writes the donor term and the candidate `splice_region_variant` in its place.
fn past_the_span_rows(location: &str, allele: &str, cdna: &str, cds: &str) -> bool {
    let row = |terms: &str| row(location, allele, terms, cdna, cds);
    splice_term_past_the_record_span_pair(
        &row("splice_donor_variant,frameshift_variant"),
        &row("frameshift_variant,splice_region_variant"),
    )
}

/// Whether the class takes two rows of three positions replaced by five letters in the
/// transcript's sequence, the reference's terms and the candidate's as given.
fn past_the_span_terms(reference: &str, candidate: &str) -> bool {
    let row = |terms: &str| row("7:100998-101000", "GTGTG", terms, "998-1000", "898-900");
    splice_term_past_the_record_span_pair(&row(reference), &row(candidate))
}

#[test]
fn the_past_the_span_class_takes_a_donor_term_gone_and_the_region_term_come() {
    assert!(past_the_span_rows(
        "7:100998-101000",
        "GTGTG",
        "998-1000",
        "898-900"
    ));
    assert!(past_the_span_rows(
        "7:100998-101000",
        "GTGTG",
        "998-1000/2000",
        "898-900/1800"
    ));
}

#[test]
fn the_past_the_span_class_takes_terms_of_the_six_gone_and_none_come() {
    assert!(past_the_span_terms(
        "splice_donor_variant,frameshift_variant,splice_donor_region_variant",
        "frameshift_variant"
    ));
    assert!(past_the_span_terms(
        "splice_polypyrimidine_tract_variant,frameshift_variant",
        "frameshift_variant"
    ));
    assert!(past_the_span_terms(
        "splice_donor_variant,frameshift_variant,intron_variant",
        "frameshift_variant"
    ));
}

#[test]
fn the_past_the_span_class_refuses_a_span_across_an_exons_end() {
    assert!(!past_the_span_rows(
        "7:100998-101000",
        "GTGTG",
        "999-?",
        "899-?"
    ));
    assert!(!past_the_span_rows(
        "7:100998-101000",
        "GTGTG",
        "?-1000",
        "?-900"
    ));
}

#[test]
fn the_past_the_span_class_refuses_a_span_inside_an_intron() {
    assert!(!past_the_span_rows("7:100998-101000", "GTGTG", "-", "-"));
}

#[test]
fn the_past_the_span_class_refuses_an_allele_no_longer_than_the_span() {
    assert!(!past_the_span_rows(
        "7:100998-101000",
        "GTG",
        "998-1000",
        "898-900"
    ));
    assert!(!past_the_span_rows(
        "7:100998-101000",
        "GT",
        "998-1000",
        "898-900"
    ));
}

#[test]
fn the_past_the_span_class_refuses_a_location_of_one_position() {
    assert!(!past_the_span_rows("7:101000", "GTGTG", "1000", "900"));
}

#[test]
fn the_past_the_span_class_refuses_a_location_of_two_neighbouring_positions() {
    assert!(!past_the_span_rows(
        TWO_POSITIONS,
        "GTGTG",
        "999-1000",
        "899-900"
    ));
}

#[test]
fn the_first_column_tells_two_replaced_bases_from_an_insertion_between_them() {
    let spelt = |name: &str| {
        let row = |terms: &str| {
            pair_row(&[
                name,
                TWO_POSITIONS,
                "GTGTG",
                "ENSG1",
                "ENST00000635895",
                "Transcript",
                terms,
                "999-1000",
                "899-900",
            ])
        };
        splice_term_past_the_record_span_pair(
            &row("splice_donor_variant,frameshift_variant"),
            &row("frameshift_variant,splice_region_variant"),
        )
    };
    assert!(spelt("7_100999_TA/GTGTG"), "two bases replaced by five");
    assert!(!spelt("7_100999_T/TGTGTG"), "an insertion of five");
}

#[test]
fn the_past_the_span_class_refuses_a_deletion() {
    assert!(!past_the_span_rows(
        "7:100998-101000",
        "-",
        "998-1000",
        "898-900"
    ));
}

#[test]
fn the_past_the_span_class_refuses_a_candidate_that_keeps_one_of_the_six() {
    assert!(!past_the_span_terms(
        "splice_donor_variant,frameshift_variant,splice_donor_region_variant",
        "frameshift_variant,splice_donor_region_variant"
    ));
    assert!(!past_the_span_terms(
        "splice_donor_variant,frameshift_variant,intron_variant",
        "frameshift_variant,intron_variant"
    ));
}

#[test]
fn the_past_the_span_class_refuses_a_donor_region_term_that_comes() {
    assert!(!past_the_span_terms(
        "splice_donor_variant,frameshift_variant",
        "frameshift_variant,splice_donor_region_variant"
    ));
}

#[test]
fn the_past_the_span_class_refuses_the_region_term_where_the_tract_term_alone_went() {
    assert!(!past_the_span_terms(
        "splice_polypyrimidine_tract_variant,frameshift_variant",
        "frameshift_variant,splice_region_variant"
    ));
}

#[test]
fn the_past_the_span_class_refuses_a_further_differing_term() {
    assert!(!past_the_span_terms(
        "splice_donor_variant,frameshift_variant",
        "stop_gained,frameshift_variant,splice_region_variant"
    ));
    assert!(!past_the_span_terms(
        "splice_donor_variant,frameshift_variant,stop_lost",
        "frameshift_variant,splice_region_variant"
    ));
}

#[test]
fn the_past_the_span_class_refuses_rows_whose_transcript_positions_differ() {
    let row =
        |terms: &str, cdna: &str, cds: &str| row("7:100998-101000", "GTGTG", terms, cdna, cds);
    let reference = row(
        "splice_donor_variant,frameshift_variant",
        "998-1000",
        "898-900",
    );
    let candidate = "frameshift_variant,splice_region_variant";
    assert!(!splice_term_past_the_record_span_pair(
        &reference,
        &row(candidate, "997-999", "898-900")
    ));
    assert!(!splice_term_past_the_record_span_pair(
        &reference,
        &row(candidate, "998-1000", "897-899")
    ));
}

#[test]
fn the_past_the_span_class_refuses_the_terms_the_other_way_round() {
    assert!(!past_the_span_terms(
        "frameshift_variant,splice_region_variant",
        "splice_donor_variant,frameshift_variant"
    ));
}

#[test]
fn the_past_the_span_class_takes_a_row_that_differs_by_a_stop_codon_classes_terms_as_well() {
    let row = |terms: &str| row("7:100998-101000", "GTGTG", terms, "998-1000", "-");
    let reference = row("splice_donor_variant,stop_retained_variant");
    let candidate = row("stop_lost,splice_region_variant");
    assert!(splice_term_past_the_record_span_pair(
        &reference, &candidate
    ));
    assert!(
        stop_codon_removed_without_stop_codon_test_pair(
            &row("stop_retained_variant,splice_region_variant"),
            &candidate
        ),
        "the rest of the difference is that class's"
    );
    assert!(
        !splice_term_past_the_record_span_pair(
            &reference,
            &row("stop_lost,missense_variant,splice_region_variant")
        ),
        "a further term outside both classes"
    );
}

/// The terms of a reference row whose donor, donor-region and intron terms come from
/// positions past the span, and of the candidate's row with `splice_region_variant` in their
/// place.
const DONOR_SIDE_TERMS: &str =
    "splice_donor_variant,frameshift_variant,splice_donor_region_variant,intron_variant";
const REGION_TERM_BACK: &str = "frameshift_variant,splice_region_variant";

/// Whether the class takes the reference's and the candidate's row of one transcript on the
/// record that replaces `reference` at `location` by `allele`, their terms as given.
fn past_the_span_record(
    reference: &str,
    location: &str,
    allele: &str,
    cdna: &str,
    terms: (&str, &str),
) -> bool {
    let row = |terms: &str| record_row(reference, location, allele, terms, cdna);
    splice_term_past_the_record_span_pair(&row(terms.0), &row(terms.1))
}

#[test]
fn with_its_record_the_class_takes_a_span_inside_an_exon_of_two_positions_and_of_three() {
    let terms = (DONOR_SIDE_TERMS, REGION_TERM_BACK);
    assert!(past_the_span_record(
        "TA",
        TWO_POSITIONS,
        "GTGTG",
        "999-1000",
        terms
    ));
    assert!(past_the_span_record(
        "TAC",
        "7:100998-101000",
        "GTGTGT",
        "998-1000",
        terms
    ));
    assert!(
        past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGT",
            "901-902",
            (
                "splice_acceptor_variant,frameshift_variant",
                REGION_TERM_BACK
            )
        ),
        "the acceptor term, on the reverse strand"
    );
}

#[test]
fn with_its_record_the_class_reads_the_reference_bases_of_either_row() {
    let reference = row(TWO_POSITIONS, "GTGTG", DONOR_SIDE_TERMS, "999-1000", "-");
    let candidate = record_row("TA", TWO_POSITIONS, "GTGTG", REGION_TERM_BACK, "999-1000");
    assert!(splice_term_past_the_record_span_pair(
        &reference, &candidate
    ));
}

#[test]
fn with_its_record_the_class_takes_a_span_across_an_exons_end() {
    assert!(past_the_span_record(
        "TA",
        "7:101000-101001",
        "GTGTG",
        "1000-?",
        (
            "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
            "splice_donor_variant,coding_sequence_variant,intron_variant"
        )
    ));
}

#[test]
fn with_its_record_the_class_takes_a_span_in_an_intron() {
    assert!(past_the_span_record(
        "TA",
        "7:101003-101004",
        "GTGT",
        "-",
        (
            "splice_donor_5th_base_variant,intron_variant",
            "splice_donor_region_variant,intron_variant"
        )
    ));
}

#[test]
fn with_its_record_the_class_takes_the_intron_term_gone_one_to_three_bases_before_the_intron() {
    for (location, allele, cdna, intron_start) in [
        (TWO_POSITIONS, "GTGTG", "999-1000", 2),
        ("7:100998-100999", "GTGTGTG", "998-999", 3),
        ("7:100997-100998", "GTGTGTGT", "997-998", 4),
    ] {
        assert!(
            past_the_span_record(
                "TA",
                location,
                allele,
                cdna,
                (DONOR_SIDE_TERMS, REGION_TERM_BACK)
            ),
            "{cdna}"
        );
        let (start, end) = cdna.split_once('-').unwrap();
        let place = record_span_placement(
            &[start, end, "-", "-"],
            &letters("TA"),
            &letters(allele),
            DONOR | DONOR_REGION | INTRON,
            REGION,
        )
        .expect("no refusal");
        assert_eq!(
            place.map(|(forward, start, ..)| (forward, start)),
            Some((true, intron_start)),
            "{cdna}"
        );
    }
}

#[test]
fn with_its_record_the_class_takes_a_donor_region_term_gone_beside_a_kept_donor_term() {
    assert!(past_the_span_record(
        "TA",
        "7:101001-101002",
        "GTGT",
        "-",
        (
            "splice_donor_variant,splice_donor_region_variant,intron_variant",
            "splice_donor_variant,intron_variant"
        )
    ));
}

#[test]
fn with_its_record_the_class_refuses_a_pure_insertion() {
    assert!(!past_the_span_record(
        "",
        TWO_POSITIONS,
        "GTG",
        "999-1000",
        ("splice_donor_variant,frameshift_variant", REGION_TERM_BACK)
    ));
}

#[test]
fn with_its_record_the_class_refuses_an_allele_as_long_as_the_reference_bases() {
    assert!(!past_the_span_record(
        "TAC",
        "7:100998-101000",
        "GTG",
        "998-1000",
        ("splice_donor_variant,frameshift_variant", REGION_TERM_BACK)
    ));
}

#[test]
fn with_its_record_the_class_refuses_one_reference_base() {
    assert!(!past_the_span_record(
        "T",
        "7:101000",
        "GTG",
        "1000",
        ("splice_donor_variant,frameshift_variant", REGION_TERM_BACK)
    ));
}

/// The record replaces two bases across an exon's end by four letters whose second is the
/// reference's: the intron's first base lies in the span at an index where the alleles
/// agree, so the candidate keeps the donor term.
#[test]
fn with_its_record_the_class_refuses_a_term_gone_whose_window_holds_a_base_of_the_span() {
    let record = |candidate: &str| {
        past_the_span_record(
            "TA",
            "7:101000-101001",
            "GAGT",
            "1000-?",
            (
                "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
                candidate,
            ),
        )
    };
    assert!(record(
        "splice_donor_variant,coding_sequence_variant,intron_variant"
    ));
    assert!(!record("coding_sequence_variant,intron_variant"));
}

#[test]
fn with_its_record_the_class_refuses_rows_no_place_of_the_intron_gives() {
    assert!(!past_the_span_record(
        "TA",
        TWO_POSITIONS,
        "GTGTG",
        "999-1000",
        (DONOR_SIDE_TERMS, "frameshift_variant")
    ));
}

#[test]
fn with_its_record_the_class_refuses_a_further_differing_term() {
    assert!(!past_the_span_record(
        "TA",
        TWO_POSITIONS,
        "GTGTG",
        "999-1000",
        (
            DONOR_SIDE_TERMS,
            "stop_gained,frameshift_variant,splice_region_variant"
        )
    ));
}

#[test]
fn with_its_record_the_class_refuses_rows_whose_transcript_positions_differ() {
    let row = |terms: &str, cdna: &str| record_row("TA", TWO_POSITIONS, "GTGTG", terms, cdna);
    assert!(!splice_term_past_the_record_span_pair(
        &row(DONOR_SIDE_TERMS, "999-1000"),
        &row(REGION_TERM_BACK, "499-500")
    ));
}

#[test]
fn two_numbers_as_far_apart_as_the_span_put_the_intron_after_the_span_or_before_it() {
    let places = places(3, ("998", "1000"));
    let starts = |strand: bool, after: bool| -> BTreeSet<i64> {
        let on =
            |(forward, start, _): &(bool, i64, i64)| *forward == strand && (*start > 0) == after;
        places
            .iter()
            .filter(|place| on(place))
            .map(|place| place.1)
            .collect()
    };
    for strand in [true, false] {
        assert_eq!(
            starts(strand, true),
            (3..10).collect(),
            "one to seven bases on"
        );
        assert_eq!(starts(strand, false), (-16..0).collect(), "before the span");
    }
    assert!(places.contains(&(true, 3, 3)), "an intron of one base");
    assert!(
        places.contains(&(false, 9, 10)),
        "a two-base intron seven bases on"
    );
    assert!(
        places.contains(&(true, 3, FAR)),
        "an intron longer than any window"
    );
    assert!(
        places.contains(&(true, -2, -1)),
        "a two-base intron that ends beside the span"
    );
    let before = places.iter().filter(|place| place.1 < 0);
    assert!(before
        .clone()
        .all(|(_, start, end)| start <= end && *end < 0));
    assert_eq!(before.count(), 2 * (1..=16).sum::<usize>());
    assert_eq!(places.len(), 1406);
}

#[test]
fn two_numbers_closer_than_the_span_put_an_intron_of_the_missing_length_inside_it() {
    assert_eq!(
        places(5, ("998", "1000")),
        [(true, 1, 2), (true, 2, 3), (false, 1, 2), (false, 2, 3)],
        "two bases of five, the span's ends kept"
    );
    assert_eq!(
        places(3, ("999", "1000")),
        [(true, 1, 1), (false, 1, 1)],
        "the middle base of three"
    );
    assert_eq!(
        places(2, ("1000", "1000")),
        [],
        "no base inside a span of two"
    );
    assert_eq!(
        places(3, ("997", "1000")),
        [],
        "numbers further apart than the span"
    );
    assert_eq!(
        places(3, ("1000", "998")),
        [],
        "the second number below the first"
    );
}

/// Eight bases inside an intron that begins two bases before the span, on the reverse
/// strand: the span's first six bases lie in a window of `splice_region` and all of them
/// in the tract's.
#[test]
fn the_last_region_assigned_decides_the_references_region_term_and_any_region_the_candidates() {
    let rows = |regions: &[(i64, i64)]| intron_rows(regions, 8, -2, FAR, false, false);
    let two_runs = differing_regions(&letters("TACGATCA"), &letters("GACGATCAGT"));
    assert_eq!(two_runs, [(0, 0), (8, 9)]);
    let (reference, candidate) = rows(&two_runs);
    assert_eq!(
        (reference, candidate),
        (TRACT | INTRON, TRACT | INTRON | REGION),
        "the first run in a window and the last in none"
    );
    assert_eq!(candidate & TRACT, TRACT, "the span holds the tract's bases");
    assert_eq!(
        rows(&[(0, 0), (4, 9)]),
        (TRACT | INTRON | REGION, TRACT | INTRON | REGION),
        "the last run in a window"
    );
    assert_eq!(
        rows(&[(7, 9)]),
        (TRACT | INTRON, TRACT | INTRON),
        "no run in a window"
    );
    assert_eq!(
        intron_rows(&[(0, 0), (2, 2), (6, 8)], 6, 2, 101, true, false),
        (
            DONOR | DONOR_5TH_BASE | INTRON,
            DONOR | DONOR_REGION | INTRON
        ),
        "a run on the intron's first two bases between two runs in a window"
    );
}

/// Eight bases of an intron replaced by eleven letters alike at four of them, five runs in
/// all: the last run lies in a window of `splice_region` and the span holds bases of the
/// tract's window at every place that gives the reference's row.
#[test]
fn with_its_record_the_class_refuses_the_tract_term_gone_from_a_span_of_five_runs() {
    let (reference, allele) = ("TACGATCA", "GATGCTGAGTG");
    assert_eq!(
        differing_regions(&letters(reference), &letters(allele)),
        [(0, 0), (2, 2), (4, 4), (6, 6), (8, 10)]
    );
    assert!(!past_the_span_record(
        reference,
        "7:101007-101014",
        allele,
        "-",
        (
            "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
            "splice_region_variant,intron_variant"
        )
    ));
}

const FIFTH_BASE_TERMS: &str = "splice_donor_5th_base_variant,frameshift_variant";

/// The windows of an intron of two bases run on past its last base: its fifth base is the
/// first base past a span of two that follows it.
#[test]
fn with_its_record_the_class_takes_a_fifth_base_term_gone_beside_a_two_base_intron_before_the_span()
{
    assert_eq!(
        intron_rows(&[(0, 4)], 2, -2, -1, true, false),
        (DONOR_5TH_BASE, DONOR_REGION)
    );
    assert!(past_the_span_record(
        "TA",
        TWO_POSITIONS,
        "GTGTG",
        "999-1000",
        (
            FIFTH_BASE_TERMS,
            "splice_donor_region_variant,frameshift_variant"
        )
    ));
}

/// On the reverse strand the fifth base of an intron of two bases that begins six bases
/// after a span of two lies three bases past the span, and the eleven bases around the
/// intron's end hold the span's last base.
#[test]
fn with_its_record_the_class_takes_a_fifth_base_term_gone_beside_a_two_base_intron_six_bases_after_the_span(
) {
    assert_eq!(
        intron_rows(&[(0, 5)], 2, 7, 8, false, false),
        (DONOR_5TH_BASE, REGION)
    );
    assert!(past_the_span_record(
        "TA",
        TWO_POSITIONS,
        "GTGTGT",
        "999-1000",
        (FIFTH_BASE_TERMS, REGION_TERM_BACK)
    ));
}

/// Six positions whose two cDNA positions are four apart hold an intron of two bases: the
/// alleles agree on it, and its fifth base is the first base past the span.
#[test]
fn with_its_record_the_class_takes_a_fifth_base_term_gone_beside_an_intron_inside_the_span() {
    let record = |cdna: &str| {
        past_the_span_record(
            "TACGAT",
            "7:100999-101004",
            "GACGATGT",
            cdna,
            (
                "splice_donor_5th_base_variant,coding_sequence_variant",
                "splice_donor_region_variant,coding_sequence_variant",
            ),
        )
    };
    assert_eq!(
        intron_rows(&[(0, 0), (6, 7)], 6, 2, 3, true, false),
        (DONOR_5TH_BASE, DONOR_REGION)
    );
    assert!(record("999-1002"));
    assert!(!record("999-998"), "the second number below the first");
}

#[test]
fn an_allele_that_opens_with_its_whole_reference_allele_writes_an_insertion_after_the_span() {
    let after_the_span = |reference: &str, allele: &str| {
        insertion_after_the_span(
            &differing_regions(&letters(reference), &letters(allele)),
            reference.len() as i64,
        )
    };
    assert!(after_the_span("TA", "TAGTG"));
    assert!(after_the_span("TA", "TAG"), "one further base");
    assert!(after_the_span("TAC", "TACGT"), "three reference bases");
    assert!(!after_the_span("T", "TGTG"), "one reference base");
    assert!(!after_the_span("TA", "TGGTG"), "the second base differs");
    assert!(!after_the_span("TA", "GAGTG"), "the first base differs");
    assert!(!after_the_span("TA", "TA"), "no further base");
    assert!(!after_the_span("TAC", "TA"), "a shorter allele");
}

/// A span of two bases beside one intron, under an allele that holds further bases after
/// them: the reference's terms come from the bases after the span, the candidate's from the
/// insertion between the span's last base and the base after it.
#[test]
fn an_insertion_after_the_span_is_read_by_its_two_neighbours() {
    let three_further =
        |start: i64, end: i64, forward: bool| intron_rows(&[(2, 4)], 2, start, end, forward, false);
    assert_eq!(
        three_further(-15, FAR, false),
        (INTRON, INTRON | TRACT),
        "after the tract's last base, on the reverse strand"
    );
    assert_eq!(
        three_further(-16, FAR, false),
        (INTRON, INTRON),
        "one base past the tract's last base"
    );
    assert_eq!(
        three_further(-FAR, 18, true),
        (INTRON | TRACT, INTRON | TRACT),
        "before the tract's first base, on the forward strand"
    );
    assert_eq!(
        three_further(-FAR, 3, true),
        (ACCEPTOR, INTRON | TRACT | REGION),
        "between the tract's innermost base and the acceptor bases"
    );
    assert_eq!(
        three_further(-FAR, 3, false),
        (DONOR, INTRON | REGION),
        "between the third intronic base and the donor bases, on the reverse strand"
    );
    assert_eq!(
        three_further(1, FAR, true),
        (DONOR | DONOR_REGION | INTRON, DONOR),
        "between the two donor bases"
    );
    assert_eq!(
        three_further(1, FAR, false),
        (ACCEPTOR | INTRON, ACCEPTOR),
        "between the two acceptor bases, on the reverse strand"
    );
    let two_further = |start: i64| intron_rows(&[(2, 3)], 2, start, FAR, true, false);
    assert_eq!(
        two_further(2),
        (DONOR, REGION),
        "between an exon's last base and the donor bases"
    );
    assert_eq!(
        two_further(0),
        (DONOR_REGION | INTRON, INTRON | REGION),
        "after the two donor bases"
    );
    assert_eq!(
        intron_rows(&[(0, 0)], 1, -15, FAR, false, false),
        (INTRON | TRACT, INTRON | TRACT),
        "one reference base is read as a change of that base"
    );
}

/// The reference's terms and the candidate's at the three places where the insertion's own
/// neighbours give the candidate a term the bases after the span do not give the reference.
const AFTER_THE_SPAN_PLACES: [(&str, (&str, &str)); 3] = [
    (
        "after the tract's last base, on the reverse strand",
        (
            "intron_variant",
            "splice_polypyrimidine_tract_variant,intron_variant",
        ),
    ),
    (
        "between the tract's innermost base and the acceptor bases",
        (
            "splice_acceptor_variant",
            "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
        ),
    ),
    (
        "between the third intronic base and the donor bases, on the reverse strand",
        (
            "splice_donor_variant",
            "splice_region_variant,intron_variant",
        ),
    ),
];

/// Whether the class takes the two rows of an allele at the two positions of an intron
/// whose reference bases are `reference`, their terms as given.
fn after_the_span_record(reference: &str, allele: &str, terms: (&str, &str)) -> bool {
    past_the_span_record(reference, TWO_POSITIONS, allele, "-", terms)
}

#[test]
fn with_its_record_the_class_takes_an_insertion_after_the_span_at_the_places_its_terms_differ() {
    for (place, terms) in AFTER_THE_SPAN_PLACES {
        assert!(after_the_span_record("TA", "TAGTG", terms), "{place}");
        assert!(
            after_the_span_record("TA", "TAG", terms),
            "{place}, one further base"
        );
    }
    assert!(
        past_the_span_record(
            "TAC",
            "7:100998-101000",
            "TACGTG",
            "-",
            AFTER_THE_SPAN_PLACES[0].1
        ),
        "three reference bases"
    );
}

#[test]
fn with_its_record_the_class_takes_the_terms_an_insertion_between_two_site_bases_withholds() {
    assert!(
        past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "TAGTG",
            "1000-?",
            (
                "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
                "splice_donor_variant,coding_sequence_variant"
            )
        ),
        "between the two donor bases"
    );
    assert!(
        past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "TAGTG",
            "?-1000",
            (
                "splice_acceptor_variant,coding_sequence_variant,intron_variant",
                "splice_acceptor_variant,coding_sequence_variant"
            )
        ),
        "between the two acceptor bases, on the reverse strand"
    );
}

/// The geometry of line 166 of `tests/golden/116/TAIR10-arabidopsis/variants.vcf`, in
/// made-up bases: three reference bases, a first alternate allele of one base and a second
/// that adds two bases to the reference allele. The second's rows stand at the two positions
/// after the shared first base, an exon's last two bases on one transcript and the two donor
/// bases of an intron on another.
#[test]
fn with_its_record_the_class_takes_both_rows_of_a_second_alternate_that_opens_with_its_reference() {
    let row = |name: &str, terms: &str, cdna: &str, cds: &str| {
        pair_row(&[
            name,
            "4:10523300-10523301",
            "TAGT",
            "AT4G19230",
            "AT4G19230.2",
            "Transcript",
            terms,
            cdna,
            cds,
        ])
    };
    let spelt = "4_10523299_CTA/C/CTAGT";
    let exon = |name: &str, terms: &str| row(name, terms, "1444-1445", "1299-1300");
    assert!(splice_term_past_the_record_span_pair(
        &exon(spelt, "splice_donor_variant,frameshift_variant"),
        &exon(spelt, "frameshift_variant,splice_region_variant")
    ));
    assert!(
        !splice_term_past_the_record_span_pair(
            &exon(spelt, "splice_donor_variant,frameshift_variant"),
            &exon(spelt, "frameshift_variant")
        ),
        "a candidate's row without the insertion's splice_region_variant"
    );
    let intron = |terms: &str| row(spelt, terms, "-", "-");
    assert!(splice_term_past_the_record_span_pair(
        &intron("splice_donor_region_variant,intron_variant"),
        &intron("intron_variant,splice_region_variant")
    ));
    assert!(
        !splice_term_past_the_record_span_pair(
            &exon("rs1", "splice_donor_variant,frameshift_variant"),
            &exon("rs1", "frameshift_variant,splice_region_variant")
        ),
        "two positions under an identifier, without the record"
    );
}

#[test]
fn an_allele_that_opens_with_its_one_reference_base_writes_an_insertion_after_it() {
    let after_one_base = |reference: &str, allele: &str| {
        insertion_after_one_base(&letters(reference), &letters(allele))
    };
    assert!(after_one_base("T", "TGTG"));
    assert!(after_one_base("T", "TG"), "one further base");
    assert!(!after_one_base("T", "GTG"), "the base differs");
    assert!(!after_one_base("T", "T"), "no further base");
    assert!(!after_one_base("T", "-"), "a deletion");
    assert!(!after_one_base("TA", "TAG"), "two reference bases");
    assert!(!after_one_base("", "GTG"), "a pure insertion");
}

/// The reference's terms and the candidate's for an allele that opens with its one
/// reference base, at the places where the base read as changed and the insertion after
/// it give different terms, with the cDNA position of the base.
const AFTER_ONE_BASE_PLACES: [(&str, &str, (&str, &str)); 10] = [
    (
        "the base before the tract's outermost base",
        "-",
        (
            "intron_variant",
            "splice_polypyrimidine_tract_variant,intron_variant",
        ),
    ),
    (
        "the intron's last base, before the exon",
        "-",
        ("splice_acceptor_variant", "splice_region_variant"),
    ),
    (
        "the exon's third base",
        "3",
        (
            "splice_region_variant,non_coding_transcript_exon_variant",
            "non_coding_transcript_exon_variant",
        ),
    ),
    (
        "the donor's second base",
        "-",
        (
            "splice_donor_variant",
            "splice_region_variant,intron_variant",
        ),
    ),
    (
        "the intron's fifth base",
        "-",
        (
            "splice_donor_5th_base_variant,intron_variant",
            "splice_donor_region_variant,intron_variant",
        ),
    ),
    (
        "the intron's sixth base",
        "-",
        (
            "splice_donor_region_variant,intron_variant",
            "splice_region_variant,intron_variant",
        ),
    ),
    (
        "the intron's eighth base",
        "-",
        ("splice_region_variant,intron_variant", "intron_variant"),
    ),
    (
        "the acceptor's base beside the tract, on the reverse strand",
        "-",
        (
            "splice_acceptor_variant",
            "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
        ),
    ),
    (
        "the eighth base from the intron's 3' end, on the reverse strand",
        "-",
        (
            "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
            "splice_polypyrimidine_tract_variant,intron_variant",
        ),
    ),
    (
        "the donor's first base, before the exon on the reverse strand",
        "-",
        ("splice_donor_variant", "splice_region_variant"),
    ),
];

#[test]
fn with_its_record_the_class_takes_an_insertion_after_one_base_at_the_places_its_terms_differ() {
    for (place, cdna, terms) in AFTER_ONE_BASE_PLACES {
        assert!(
            past_the_span_record("T", "7:101000", "TGTG", cdna, terms),
            "{place}"
        );
        assert!(
            past_the_span_record("T", "7:101000", "TG", cdna, terms),
            "{place}, one further base"
        );
        assert!(
            !past_the_span_record("T", "7:101000", "GTGT", cdna, terms),
            "{place}, an allele that does not open with the reference base"
        );
        let row = |terms: &str| row("7:101000", "TGTG", terms, cdna, "-");
        assert!(
            !splice_term_past_the_record_span_pair(&row(terms.0), &row(terms.1)),
            "{place}, without the record"
        );
    }
}

/// A record of two alternate alleles as the first column spells it: the second opens with
/// the reference allele and adds three bases, and its row stands at the one position after
/// the shared first base, the donor's second base.
#[test]
fn with_its_record_the_class_takes_the_row_of_an_alternate_that_opens_with_one_kept_base() {
    let row = |name: &str, terms: &str| {
        pair_row(&[
            name,
            "22:20996968",
            "TGTG",
            "ENSG00000099949",
            "ENST00000646124",
            "Transcript",
            terms,
            "-",
            "-",
        ])
    };
    let spelt = "22_20996967_CT/C/CTGTG";
    assert!(splice_term_past_the_record_span_pair(
        &row(spelt, "splice_donor_variant"),
        &row(spelt, "splice_region_variant,intron_variant")
    ));
    assert!(
        !splice_term_past_the_record_span_pair(
            &row(spelt, "splice_donor_variant"),
            &row(spelt, "splice_donor_region_variant,intron_variant")
        ),
        "a candidate's row with a term no place of the insertion gives beside the donor term"
    );
    assert!(
        !splice_term_past_the_record_span_pair(
            &row("rs1", "splice_donor_variant"),
            &row("rs1", "splice_region_variant,intron_variant")
        ),
        "one position under an identifier, without the record"
    );
}

#[test]
fn with_its_record_the_class_refuses_a_term_an_insertion_after_one_base_does_not_explain() {
    let refused = [
        (
            "a donor term the candidate alone holds",
            "-",
            "intron_variant",
            "intron_variant,splice_donor_variant",
        ),
        (
            "another term the candidate alone holds",
            "-",
            "splice_acceptor_variant",
            "splice_region_variant,missense_variant",
        ),
        (
            "another term the reference alone holds",
            "-",
            "splice_donor_variant,non_coding_transcript_variant",
            "splice_region_variant,intron_variant",
        ),
        (
            "the intron term the insertion's place gives, missing after the donor's second base",
            "-",
            "splice_donor_variant,intron_variant",
            "splice_region_variant",
        ),
        (
            "splice_region_variant gone from a base inside the transcript's sequence, with the intron term",
            "3",
            "splice_region_variant,intron_variant",
            "intron_variant",
        ),
        (
            "the same terms on both rows",
            "-",
            "splice_region_variant,intron_variant",
            "splice_region_variant,intron_variant",
        ),
    ];
    for (shape, cdna, reference, candidate) in refused {
        assert!(
            !past_the_span_record("T", "7:101000", "TGTG", cdna, (reference, candidate)),
            "{shape}"
        );
    }
}

/// The kept base is an intron of one base: the insertion after it stands between the two
/// bases of that intron's first site window, and the base read as changed draws nothing.
#[test]
fn with_its_record_the_class_takes_the_site_term_of_an_insertion_after_a_one_base_intron() {
    let lone = "non_coding_transcript_variant";
    for site in ["splice_acceptor_variant", "splice_donor_variant"] {
        let candidate = format!("{lone},{site}");
        let rows = (lone, candidate.as_str());
        assert!(
            past_the_span_record("T", "7:101000", "TGT", "-", rows),
            "{site}"
        );
        assert!(
            !past_the_span_record("T", "7:101000", "TGT", "3", rows),
            "{site} on a base of the transcript's sequence"
        );
        assert!(
            !past_the_span_record("T", "7:101000", "GGT", "-", rows),
            "{site} under an allele that does not keep the base"
        );
    }
    let refused = |reference: &str, candidate: &str| {
        !past_the_span_record("T", "7:101000", "TGT", "-", (reference, candidate))
    };
    let site = "splice_acceptor_variant";
    for term in INTRON_ROW_TERMS.iter().filter(|term| **term != site) {
        let held = format!("{lone},{term}");
        assert!(
            refused(&held, &format!("{held},{site}")),
            "a site term beside the reference's {term}"
        );
        assert!(
            refused(&held, &format!("{lone},{site}")),
            "a site term in place of the reference's {term}"
        );
        assert!(
            refused(lone, &format!("{held},{site}")),
            "a site term with {term}"
        );
    }
    for beside in ["upstream_gene_variant", "downstream_gene_variant"] {
        assert!(
            refused(beside, &format!("{beside},{site}")),
            "a site term on a row of {beside}"
        );
    }
}

/// The rows of `--terms display`: both site terms print one name, the fifth base, the
/// donor region, the tract and `splice_region_variant` another, the intron term a third.
#[test]
fn under_the_display_names_the_class_reads_the_names_the_seven_terms_print() {
    assert!(
        past_the_span_record(
            "T",
            "7:101000",
            "TT",
            "3",
            ("FRAMESHIFT_CODING,SPLICE_SITE", "FRAMESHIFT_CODING")
        ),
        "the exon's third base, one base added after it"
    );
    assert!(
        past_the_span_record(
            "T",
            "7:101000",
            "TGTG",
            "-",
            ("ESSENTIAL_SPLICE_SITE", "SPLICE_SITE,INTRONIC")
        ),
        "the donor's second base, three bases added after it"
    );
    assert!(
        past_the_span_record(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            (
                "ESSENTIAL_SPLICE_SITE,FRAMESHIFT_CODING,SPLICE_SITE,INTRONIC",
                "FRAMESHIFT_CODING,SPLICE_SITE"
            )
        ),
        "two replaced bases inside an exon, the donor and intron names gone"
    );
    let refused = [
        (
            "the same names on both rows, the fifth base against the donor region",
            "-",
            "SPLICE_SITE,INTRONIC",
            "SPLICE_SITE,INTRONIC",
        ),
        (
            "another term that differs",
            "3",
            "FRAMESHIFT_CODING,SPLICE_SITE",
            "STOP_GAINED",
        ),
        (
            "names no place of the intron gives",
            "-",
            "INTRONIC",
            "ESSENTIAL_SPLICE_SITE",
        ),
        (
            "a Sequence Ontology name beside a display name",
            "-",
            "ESSENTIAL_SPLICE_SITE",
            "splice_region_variant,INTRONIC",
        ),
        (
            "the name of the donor and the acceptor term on the candidate's row alone",
            "-",
            "WITHIN_NON_CODING_GENE",
            "ESSENTIAL_SPLICE_SITE,WITHIN_NON_CODING_GENE",
        ),
    ];
    for (shape, cdna, reference, candidate) in refused {
        assert!(
            !past_the_span_record("T", "7:101000", "TGTG", cdna, (reference, candidate)),
            "{shape}"
        );
    }
    assert!(
        !past_the_span_record(
            "TAC",
            "7:100998-101000",
            "GAT",
            "-",
            ("INTRONIC", "SPLICE_SITE,INTRONIC")
        ),
        "three replaced bases under an allele as long, its first and last letter differing"
    );
    assert!(
        !past_the_span_record(
            "TAC",
            "7:100998-101000",
            "GACTG",
            "-",
            ("INTRONIC", "SPLICE_SITE,INTRONIC")
        ),
        "a longer allele that does not open with its reference, no name gone from the reference's row"
    );
    let unread = |terms: &str| row("7:101000", "TT", terms, "3", "-");
    assert!(
        !splice_term_past_the_record_span_pair(
            &unread("FRAMESHIFT_CODING,SPLICE_SITE"),
            &unread("FRAMESHIFT_CODING")
        ),
        "without the record"
    );
}

#[test]
fn with_its_record_the_class_refuses_an_allele_that_does_not_open_with_its_whole_reference() {
    for (place, terms) in AFTER_THE_SPAN_PLACES {
        assert!(
            !after_the_span_record("TA", "TGGTG", terms),
            "{place}, the second base differs"
        );
        assert!(
            !after_the_span_record("TA", "GAGTG", terms),
            "{place}, the first base differs"
        );
        let row = |terms: &str| row(TWO_POSITIONS, "TAGTG", terms, "-", "-");
        assert!(
            !splice_term_past_the_record_span_pair(&row(terms.0), &row(terms.1)),
            "{place}, without the record"
        );
    }
}

#[test]
fn with_its_record_the_class_refuses_a_term_an_insertion_after_the_span_does_not_explain() {
    let refused = [
        (
            "a donor term the candidate alone holds",
            "intron_variant",
            "splice_polypyrimidine_tract_variant,intron_variant,splice_donor_variant",
        ),
        (
            "splice_region_variant the reference alone holds",
            "splice_region_variant,intron_variant",
            "splice_polypyrimidine_tract_variant,intron_variant",
        ),
        (
            "another term the candidate alone holds",
            "splice_acceptor_variant",
            "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant,missense_variant",
        ),
        (
            "another term the reference alone holds",
            "splice_donor_variant,non_coding_transcript_variant",
            "splice_region_variant,intron_variant",
        ),
        (
            "the intron term the insertion's place gives, missing beside the acceptor bases",
            "splice_acceptor_variant",
            "splice_region_variant,splice_polypyrimidine_tract_variant",
        ),
        (
            "the intron term the insertion's place gives, missing beside the donor bases",
            "splice_donor_variant",
            "splice_region_variant",
        ),
        (
            "the same terms on both rows",
            "intron_variant",
            "intron_variant",
        ),
    ];
    for (shape, reference, candidate) in refused {
        assert!(
            !after_the_span_record("TA", "TAGTG", (reference, candidate)),
            "{shape}"
        );
    }
}

/// The pairing's two rules over two sets of terms (`registry_class`) take no pair on which the
/// candidate alone holds the tract term, and each class's rule is consulted by its name: the
/// pair of an insertion written after the span reaches the class's rule.
#[test]
fn no_rule_over_two_sets_stands_before_the_class_for_a_tract_term_the_candidate_alone_holds() {
    let (reference, candidate) = AFTER_THE_SPAN_PLACES[0].1;
    assert_eq!(
        registry_class(
            &consequence_set(reference),
            &consequence_set(candidate),
            "TAGTG",
            &[START_COOCCURRENCE_SWAP, SPLICE_FAMILY_SWAP],
        ),
        None
    );
    assert!(after_the_span_record("TA", "TAGTG", (reference, candidate)));
}

#[test]
fn with_its_record_the_class_reads_a_pure_insertion_as_no_insertion_after_the_span() {
    assert!(!insertion_after_the_span(
        &differing_regions(&letters(""), &letters("GTG")),
        0
    ));
    for (place, terms) in AFTER_THE_SPAN_PLACES {
        assert!(!after_the_span_record("", "GTG", terms), "{place}");
    }
}

/// An intron of thirteen bases or fewer is a frameshift intron
/// (BaseTranscriptVariation.pm:1032-1035): an insertion whose two neighbours it holds sets
/// nothing (BaseTranscriptVariationAllele.pm:138-141, :171-174).
#[test]
fn an_insertion_after_the_span_sets_no_term_inside_a_frameshift_intron() {
    let three_further = |start: i64, end: i64| intron_rows(&[(2, 4)], 2, start, end, true, false);
    assert_eq!(
        three_further(-3, 6),
        (0, 0),
        "both neighbours inside an intron of ten bases"
    );
    assert_eq!(
        three_further(2, 5),
        (0, REGION),
        "between an exon's last base and the first base of an intron of four bases"
    );
}

/// A transcript with a frameshift intron reads its exons twelve bases wider
/// (BaseTranscriptVariation.pm:861), the reference against the span and the candidate against
/// the insertion's two neighbours: where the intron ends twelve bases after the span, the base
/// after the span lies on the wider exon and the span's last base does not.
#[test]
fn an_insertion_after_the_span_loses_the_tract_term_on_the_wider_exons() {
    let three_further = |wide_exons: bool| intron_rows(&[(2, 4)], 2, -5, FAR, false, wide_exons);
    assert_eq!(
        three_further(false),
        (TRACT | INTRON | REGION, TRACT | INTRON | REGION)
    );
    assert_eq!(
        three_further(true),
        (INTRON | REGION, INTRON | REGION),
        "five bases after the intron's first base, on the reverse strand"
    );
    let one_further = |wide_exons: bool| intron_rows(&[(2, 2)], 2, -12, 13, true, wide_exons);
    assert_eq!(one_further(false), (TRACT | INTRON, TRACT | INTRON));
    assert_eq!(
        one_further(true),
        (TRACT | INTRON, INTRON),
        "the candidate alone loses the tract term"
    );
    let place = record_span_placement(
        &["-", "-", "-", "-"],
        &letters("TA"),
        &letters("TAG"),
        TRACT | INTRON,
        INTRON,
    )
    .expect("no refusal");
    assert_eq!(place, Some((true, -12, 13, true)));
    assert!(after_the_span_record(
        "TA",
        "TAG",
        (
            "splice_polypyrimidine_tract_variant,intron_variant",
            "intron_variant"
        )
    ));
}

const UNKNOWN_BESIDE_REGION: &str = "splice_region_variant,coding_sequence_variant";
const KEPT_BESIDE_REGION: &str = "splice_region_variant,stop_retained_variant";

/// A row of three bases, the last unknown, between an intron's last base and the exon base
/// that ends a stop codon: the first column, the Location, the allele and the two position
/// columns, under `terms`.
fn unknown_base_row(columns: [&str; 5], terms: &str) -> Entry {
    let [name, location, allele, cdna, cds] = columns;
    pair_row(&[
        name,
        location,
        allele,
        "ENSG1",
        "ENST00000399326",
        "Transcript",
        terms,
        cdna,
        cds,
    ])
}

const UNKNOWN_BASE_ROW: [&str; 5] = [
    "21_35201927_G/GAAN",
    "21:35201927-35201928",
    "AAN",
    "3048-3049",
    "3047-3048",
];

/// Whether the rule of `stop_codon_overlap_insertion_range` takes the row of `columns` under
/// the reference's terms and the candidate's.
fn in_the_insertion_range(columns: [&str; 5], reference: &str, candidate: &str) -> bool {
    stop_codon_overlap_insertion_range_pair(
        &unknown_base_row(columns, reference),
        &unknown_base_row(columns, candidate),
    )
}

/// Beside `splice_region_variant` the first pair of the class is taken where the allele
/// holds an unknown base, wherever it stands, with the record's alleles and under an
/// identifier; the split codon's rule refuses those rows and takes the same row of known
/// bases, which this one refuses.
#[test]
fn the_insertion_range_class_takes_its_first_pair_beside_splice_region_variant_for_an_unknown_base()
{
    let (reference, candidate) = (UNKNOWN_BESIDE_REGION, KEPT_BESIDE_REGION);
    let [name, location, allele, cdna, cds] = UNKNOWN_BASE_ROW;
    for columns in [
        UNKNOWN_BASE_ROW,
        ["1993671", location, allele, cdna, cds],
        ["1993671", location, "ANA", cdna, cds],
        ["21_35201927_G/GNAAAAT", location, "NAAAAT", cdna, cds],
        [name, location, allele, "3048-3049/3600", "3047-?/3048"],
    ] {
        assert!(
            in_the_insertion_range(columns, reference, candidate),
            "{columns:?}"
        );
        assert!(
            !stop_codon_overlap_genomic_span_pair(
                &unknown_base_row(columns, reference),
                &unknown_base_row(columns, candidate)
            ),
            "{columns:?}"
        );
    }
    for columns in [
        ["21_35201927_G/GAAT", location, "AAT", cdna, cds],
        ["1993671", location, "AAA", cdna, cds],
    ] {
        assert!(
            !in_the_insertion_range(columns, reference, candidate),
            "{columns:?}"
        );
        assert!(
            stop_codon_overlap_genomic_span_pair(
                &unknown_base_row(columns, reference),
                &unknown_base_row(columns, candidate)
            ),
            "{columns:?}"
        );
    }
}

/// The sides exchanged are no row of the class: the reference's `stop_retained_variant`
/// against the candidate's `coding_sequence_variant` on an allele holding an unknown base is
/// the shape of a stop codon test that is not made.
#[test]
fn the_first_pair_beside_splice_region_variant_is_not_taken_with_the_sides_exchanged() {
    let (reference, candidate) = (UNKNOWN_BESIDE_REGION, KEPT_BESIDE_REGION);
    assert!(in_the_insertion_range(
        UNKNOWN_BASE_ROW,
        reference,
        candidate
    ));
    assert!(!in_the_insertion_range(
        UNKNOWN_BASE_ROW,
        candidate,
        reference
    ));
    assert!(stop_retained_without_stop_codon_test_pair(
        &unknown_base_row(UNKNOWN_BASE_ROW, candidate),
        &unknown_base_row(UNKNOWN_BASE_ROW, reference)
    ));
}

/// A row off the shape is not taken: an allele holding a letter that is neither a base nor
/// `N`, a deletion, a row without a cDNA end or a CDS start, a Location of one position or
/// of three, a record whose reference allele keeps a base, the rejoined row of an allele
/// `--minimal` split from a record of two alternates, and terms that are not the pair.
#[test]
fn a_row_beside_splice_region_variant_off_the_shape_is_not_taken() {
    let (reference, candidate) = (UNKNOWN_BESIDE_REGION, KEPT_BESIDE_REGION);
    assert!(in_the_insertion_range(
        UNKNOWN_BASE_ROW,
        reference,
        candidate
    ));
    let [name, location, allele, cdna, cds] = UNKNOWN_BASE_ROW;
    for columns in [
        ["1993671", location, "AAR", cdna, cds],
        ["1993671", "21:35201928-35201930", "-", cdna, cds],
        [name, location, allele, "3048-?", "3047-?"],
        [name, location, allele, cdna, "?-3048"],
        ["1993671", "21:35201928", allele, cdna, cds],
        ["1993671", "21:35201927-35201929", allele, cdna, cds],
        ["21_35201927_GA/GAAN", location, allele, cdna, cds],
        [
            "21_35201926_AG/AGAAN/AGNAG",
            "21:35201927",
            allele,
            "-",
            "-",
        ],
    ] {
        assert!(
            !in_the_insertion_range(columns, reference, candidate),
            "{columns:?}"
        );
    }
    for (reference, candidate) in [
        (UNKNOWN_BESIDE_REGION, "splice_region_variant,stop_lost"),
        (UNKNOWN_BESIDE_REGION, "stop_retained_variant"),
        (
            UNKNOWN_BESIDE_REGION,
            "inframe_insertion,splice_region_variant,stop_retained_variant",
        ),
        ("splice_region_variant,stop_gained", KEPT_BESIDE_REGION),
        ("splice_region_variant", KEPT_BESIDE_REGION),
    ] {
        assert!(
            !in_the_insertion_range(UNKNOWN_BASE_ROW, reference, candidate),
            "{reference} {candidate}"
        );
    }
}

/// Four unknown bases.
const UNKNOWN_BASES: &str = "NNNN";

/// One Transcript row of an insertion after 21:10567700 at a codon boundary: its allele,
/// its terms, its residues and the alternate half of its Codons cell.
fn settled_stop_row(allele: &str, terms: &str, residues: &str, codons: &str) -> Entry {
    let text = format!(
        "{HEADER}v\t21:10567700-10567701\t{allele}\tENSG1\tENST1\tTranscript\t{terms}\t775-776\t423-424\t141-142\t{residues}\t-/{codons}\t-\t-\n"
    );
    let rows = parse_default(&text).expect("a row of fourteen columns");
    Entry::clone(&rows[0])
}

fn settled_stop_pair(e: &Entry, a: &Entry) -> bool {
    super::stop_codon_settled_before_unknown_bases_pair(e, a).expect("no refusal")
}

/// The alternate peptide of a codons half is the one `tools/vep-formats-compare` reads from
/// the same half.
#[test]
fn a_stop_codon_before_the_first_unknown_codon_gives_the_alternate_peptide() {
    let peptide = |codons: &str| settled_alternate_peptide(codons);
    // A sense codon, a stop codon, an unknown codon and a codon whose first base is
    // unknown, with and without bases past the last whole codon.
    assert_eq!(peptide("GAATAGNNNNCA").as_deref(), Some("E*XX"));
    assert_eq!(peptide("GAATAGNNNN").as_deref(), Some("E*XX"));
    // A codon every reading of which is one residue keeps it, and the stop alone takes no
    // `X` for the bases past it.
    assert_eq!(peptide("TAACTNGGA").as_deref(), Some("*LG"));
    assert_eq!(peptide("taaN").as_deref(), Some("*"));
    assert_eq!(peptide("TAGRATSAATAR").as_deref(), Some("*BZ*"));
    // No stop codon among the known codons, a stop codon after the unknown bases, one that
    // holds the unknown base, and one of the vertebrate mitochondrial code alone.
    assert_eq!(peptide("GAAGCTNNNNCA"), None);
    assert_eq!(peptide("GAANNNTAG"), None);
    assert_eq!(peptide("GAATANGGA"), None);
    assert_eq!(peptide("GAAAGANNNN"), None);
    // Known bases alone, a letter that names no base and an empty half.
    assert_eq!(peptide("GAATAGGGA"), None);
    assert_eq!(peptide("TAGJJJ"), None);
    assert_eq!(peptide(""), None);
}

/// The rule of `stop_codon_settled_before_unknown_bases` takes the frameshift and the
/// in-frame insertion of the tests of `tools/vep-formats-compare`, the in-frame insertion
/// beside a term both sides write and a stop codon alone before one unknown codon, and
/// refuses each with the sides exchanged.
#[test]
fn a_gained_stop_settled_before_unknown_bases_is_a_pair_of_the_class() {
    let frameshift = format!("GAATAG{UNKNOWN_BASES}");
    let inframe = format!("GAATAG{UNKNOWN_BASES}CA");
    for (allele, reference, candidate, residues) in [
        (
            frameshift.as_str(),
            "frameshift_variant",
            "stop_gained,frameshift_variant",
            "-/E*XX",
        ),
        (
            inframe.as_str(),
            "coding_sequence_variant",
            "stop_gained,inframe_insertion",
            "-/E*XX",
        ),
        (
            inframe.as_str(),
            "coding_sequence_variant,NMD_transcript_variant",
            "stop_gained,NMD_transcript_variant",
            "-/E*XX",
        ),
        (
            "TAGNNN",
            "coding_sequence_variant",
            "stop_gained,protein_altering_variant",
            "-/*X",
        ),
    ] {
        let e = settled_stop_row(allele, reference, "-", allele);
        let a = settled_stop_row(allele, candidate, residues, allele);
        assert!(settled_stop_pair(&e, &a), "{allele} {reference}");
        assert!(!settled_stop_pair(&a, &e), "{allele} {reference}");
    }
}

/// The rule of `stop_codon_settled_before_unknown_bases` takes a stop codon the inserted
/// bases keep as `tools/vep-formats-compare` does: the candidate's `inframe_insertion`
/// added to terms that hold `stop_retained_variant` on both sides, or the candidate's two
/// terms against the reference's `frameshift_variant` or `coding_sequence_variant` on an
/// insertion between two neighbouring positions. It refuses the sides exchanged, a
/// reference that writes `inframe_insertion` itself, another term of the reference in the
/// place of the two, and the second pair on a row that is no such insertion.
#[test]
fn a_kept_stop_before_unknown_bases_is_a_pair_of_the_class() {
    let frameshift = format!("TAG{UNKNOWN_BASES}");
    let inframe = format!("TAG{UNKNOWN_BASES}CA");
    let kept = "inframe_insertion,stop_retained_variant";
    for (allele, reference, candidate) in [
        (frameshift.as_str(), "stop_retained_variant", kept),
        (inframe.as_str(), "stop_retained_variant", kept),
        (frameshift.as_str(), "frameshift_variant", kept),
        (inframe.as_str(), "coding_sequence_variant", kept),
        (
            frameshift.as_str(),
            "frameshift_variant,NMD_transcript_variant",
            "inframe_insertion,stop_retained_variant,NMD_transcript_variant",
        ),
    ] {
        let e = settled_stop_row(allele, reference, "-", allele);
        let a = settled_stop_row(allele, candidate, "-/*XX", allele);
        assert!(settled_stop_pair(&e, &a), "{allele} {reference}");
        assert!(!settled_stop_pair(&a, &e), "{allele} {reference}");
    }
    for (what, reference, candidate) in [
        (
            "the reference writes inframe_insertion itself",
            "inframe_insertion,stop_retained_variant",
            "inframe_insertion,stop_retained_variant,NMD_transcript_variant",
        ),
        (
            "another term in the place of the two",
            "missense_variant",
            kept,
        ),
        (
            "a gained stop in the reference",
            "stop_gained,frameshift_variant",
            kept,
        ),
        (
            "no stop_retained_variant in the candidate",
            "frameshift_variant",
            "inframe_insertion",
        ),
    ] {
        let e = settled_stop_row(&frameshift, reference, "-", &frameshift);
        let a = settled_stop_row(&frameshift, candidate, "-/*XX", &frameshift);
        assert!(!settled_stop_pair(&e, &a), "{what}");
    }
    // The second pair on a Location of three positions: no insertion between neighbours.
    let replaced = |terms: &str, residues: &str| {
        let text = format!(
            "{HEADER}v\t21:10567700-10567702\t{frameshift}\tENSG1\tENST1\tTranscript\t{terms}\t775-777\t423-425\t141-142\t{residues}\t-/{frameshift}\t-\t-\n"
        );
        Entry::clone(&parse_default(&text).expect("a row of fourteen columns")[0])
    };
    assert!(!settled_stop_pair(
        &replaced("frameshift_variant", "-"),
        &replaced(kept, "-/*XX")
    ));
    // The first pair asks no insertion: both rows hold the kept stop.
    assert!(settled_stop_pair(
        &replaced("stop_retained_variant", "-"),
        &replaced(kept, "-/*XX")
    ));
}

/// Where the reference writes `stop_retained_variant` for a stop codon test it does not
/// make, the rule takes the candidate's gained stop with `frameshift_variant`, with an
/// in-frame term or alone, beside a term both sides write as well. It refuses the sides
/// exchanged, the reference's term kept, a second term or `coding_sequence_variant` beside
/// the gained stop, a further term of the reference alone, an allele of known bases, whose
/// test is made, and a row whose allele has the length of its Location.
#[test]
fn a_gained_stop_against_an_untested_stop_retained_term_is_a_pair_of_the_class() {
    let frameshift = format!("GAATAG{UNKNOWN_BASES}");
    let inframe = format!("GAATAG{UNKNOWN_BASES}CA");
    let retained = "stop_retained_variant";
    for (allele, reference, candidate) in [
        (
            frameshift.as_str(),
            retained,
            "stop_gained,frameshift_variant",
        ),
        (inframe.as_str(), retained, "stop_gained,inframe_insertion"),
        (inframe.as_str(), retained, "stop_gained"),
        (
            frameshift.as_str(),
            "stop_retained_variant,NMD_transcript_variant",
            "stop_gained,frameshift_variant,NMD_transcript_variant",
        ),
    ] {
        let e = settled_stop_row(allele, reference, "-", allele);
        let a = settled_stop_row(allele, candidate, "-/E*XX", allele);
        assert!(settled_stop_pair(&e, &a), "{allele} {candidate}");
        assert!(!settled_stop_pair(&a, &e), "{allele} {candidate}");
    }
    let e = settled_stop_row(&frameshift, retained, "-", &frameshift);
    for (what, reference, candidate_terms) in [
        (
            "the reference's term kept",
            e.clone(),
            "stop_gained,frameshift_variant,stop_retained_variant",
        ),
        (
            "two terms beside the gained stop",
            e.clone(),
            "stop_gained,frameshift_variant,inframe_insertion",
        ),
        (
            "coding_sequence_variant beside the gained stop",
            e.clone(),
            "stop_gained,coding_sequence_variant",
        ),
        (
            "a further term of the reference alone",
            settled_stop_row(
                &frameshift,
                "splice_region_variant,stop_retained_variant",
                "-",
                &frameshift,
            ),
            "stop_gained,frameshift_variant",
        ),
    ] {
        let a = settled_stop_row(&frameshift, candidate_terms, "-/E*XX", &frameshift);
        assert!(!settled_stop_pair(&reference, &a), "{what}");
    }
    // The Allele is of known bases, so both rows show the test made: the codons alone hold
    // the unknown bases.
    let known = "GAATAGACGT";
    let tested = settled_stop_row(known, retained, "-", &frameshift);
    let gained = settled_stop_row(
        known,
        "stop_gained,frameshift_variant",
        "-/E*XX",
        &frameshift,
    );
    assert!(!settled_stop_pair(&tested, &gained));
    // Ten positions under an allele of ten letters: no insertion's and no deletion's row.
    let same_length = |terms: &str, residues: &str| {
        let text = format!(
            "{HEADER}v\t21:10567700-10567709\t{frameshift}\tENSG1\tENST1\tTranscript\t{terms}\t775-776\t423-424\t141-142\t{residues}\t-/{frameshift}\t-\t-\n"
        );
        Entry::clone(&parse_default(&text).expect("a row of fourteen columns")[0])
    };
    assert!(!settled_stop_pair(
        &same_length(retained, "-"),
        &same_length("stop_gained,frameshift_variant", "-/E*XX"),
    ));
}

/// A row off the shape is no pair of the class: residues that are not the codons'
/// translation, known bases that settle no stop codon, codons that differ between the
/// sides, any other difference of terms, residues in the reference and a row of another
/// feature type.
#[test]
fn a_row_off_the_shape_of_the_settled_stop_is_refused() {
    let allele = format!("GAATAG{UNKNOWN_BASES}CA");
    let unsettled = format!("GAAGCT{UNKNOWN_BASES}CA");
    let (unknown, gained) = ("coding_sequence_variant", "stop_gained,inframe_insertion");
    let e = settled_stop_row(&allele, unknown, "-", &allele);
    let a = settled_stop_row(&allele, gained, "-/E*XX", &allele);
    assert!(settled_stop_pair(&e, &a));
    let row = |terms: &str, residues: &str, codons: &str| {
        settled_stop_row(&allele, terms, residues, codons)
    };
    for (what, reference, candidate) in [
        (
            "residues that are not the codons' translation",
            e.clone(),
            row(gained, "-/E*XK", &allele),
        ),
        (
            "known bases that settle no stop codon before the unknown ones",
            row(unknown, "-", &unsettled),
            row(gained, "-/EAXX", &unsettled),
        ),
        (
            "codons that differ between the sides",
            e.clone(),
            row(gained, "-/E*XX", &format!("GAATAA{UNKNOWN_BASES}CA")),
        ),
        (
            "a term other than the gained stop and an in-frame term",
            e.clone(),
            row("stop_gained,missense_variant", "-/E*XX", &allele),
        ),
        (
            "two in-frame terms beside the gained stop",
            e.clone(),
            row(
                "stop_gained,inframe_insertion,protein_altering_variant",
                "-/E*XX",
                &allele,
            ),
        ),
        (
            "no gained stop",
            e.clone(),
            row("inframe_insertion", "-/E*XX", &allele),
        ),
        (
            "a frameshift the candidate does not keep",
            row("frameshift_variant", "-", &allele),
            row("stop_gained", "-/E*XX", &allele),
        ),
        (
            "residues in the reference",
            row(unknown, "-/E*XX", &allele),
            a.clone(),
        ),
        (
            "no residues in the candidate",
            e.clone(),
            row(gained, "-", &allele),
        ),
        (
            "no reference half of the residues",
            e.clone(),
            row(gained, "/E*XX", &allele),
        ),
    ] {
        assert!(!settled_stop_pair(&reference, &candidate), "{what}");
    }
    let regulatory = |entry: &Entry| Entry {
        key: entry.key.replace("Transcript", "RegulatoryFeature"),
        ..entry.clone()
    };
    assert!(!settled_stop_pair(&regulatory(&e), &regulatory(&a)));
    let without_codons = settled_stop_row(&allele, gained, "-/E*XX", &allele);
    let mut blanked = e.clone();
    blanked.fields.insert("Codons", "-");
    let mut blanked_candidate = without_codons;
    blanked_candidate.fields.insert("Codons", "-");
    assert!(!settled_stop_pair(&blanked, &blanked_candidate));
}

/// On a row of the settled stop a field of NMD or LOFTEE the candidate alone writes
/// follows the terms where the reference's hold no truncating term, and no field does where
/// `frameshift_variant` stands on both sides.
#[test]
fn a_plugin_field_follows_the_terms_of_a_settled_stop_only_where_the_reference_holds_no_truncating_term(
) {
    let inframe = format!("GAATAG{UNKNOWN_BASES}CA");
    let e = settled_stop_row(&inframe, "coding_sequence_variant", "-", &inframe);
    let a = settled_stop_row(
        &inframe,
        "stop_gained,inframe_insertion",
        "-/E*XX",
        &inframe,
    );
    for field in TRUNCATING_TERM_FIELDS {
        assert!(
            term_pair_following_cell(&e, &a, field, "-", "HC"),
            "{field}"
        );
        assert!(
            !term_pair_following_cell(&e, &a, field, "HC", "-"),
            "{field}"
        );
    }
    for field in DOWNSTREAM_FIELDS {
        assert!(
            !term_pair_following_cell(&e, &a, field, "-", "MK"),
            "{field}"
        );
    }
    let frameshift = format!("GAATAG{UNKNOWN_BASES}");
    let e = settled_stop_row(&frameshift, "frameshift_variant", "-", &frameshift);
    let a = settled_stop_row(
        &frameshift,
        "stop_gained,frameshift_variant",
        "-/E*XX",
        &frameshift,
    );
    assert!(settled_stop_pair(&e, &a));
    for field in TRUNCATING_TERM_FIELDS.iter().chain(&DOWNSTREAM_FIELDS) {
        assert!(
            !term_pair_following_cell(&e, &a, field, "-", "HC"),
            "{field}"
        );
        assert!(
            !term_pair_following_cell(&e, &a, field, "HC", "LC"),
            "{field}"
        );
    }
}

/// One Transcript row of a deletion that begins at 21:100: its terms, its two position
/// columns and its residues.
fn removed_start_row(terms: &str, cdna: &str, cds: &str, residues: &str) -> Entry {
    let text = format!(
        "{HEADER}v\t21:100-102\t-\tENSG1\tENST1\tTranscript\t{terms}\t{cdna}\t{cds}\t1\t{residues}\tATG/-\t-\t-\n"
    );
    let rows = parse_default(&text).expect("a row of fourteen columns");
    Entry::clone(&rows[0])
}

fn removed_start_pair(e: &Entry, a: &Entry) -> bool {
    super::start_codon_removed_without_five_prime_utr_pair(e, a).expect("no refusal")
}

const START_KEPT: &str = "inframe_deletion";
const START_LOST: &str = "start_lost,inframe_deletion";

/// The rule of `start_codon_removed_without_five_prime_utr` takes the rows of the tests of
/// `tools/vep-formats-compare`: a removal of whole codons from the first base of a
/// transcript whose coding sequence opens there, the candidate's terms the reference's and
/// `start_lost`, beside a term both sides write and under `--total_length`; and refuses
/// each with the sides exchanged.
#[test]
fn a_start_codon_removed_from_the_first_base_is_a_pair_of_the_class() {
    for (reference, candidate, cdna, cds, residues) in [
        (START_KEPT, START_LOST, "1-3", "1-3", "M/-"),
        (
            "inframe_deletion,NMD_transcript_variant",
            "start_lost,inframe_deletion,NMD_transcript_variant",
            "1-6",
            "1-6",
            "MA/-",
        ),
        (
            "splice_acceptor_variant,splice_donor_variant,inframe_deletion,intron_variant",
            "splice_acceptor_variant,splice_donor_variant,start_lost,inframe_deletion,intron_variant",
            "1-12",
            "1-12",
            "MAGS/-",
        ),
        (START_KEPT, START_LOST, "1-3/1500", "1-3/900", "M/-"),
    ] {
        let e = removed_start_row(reference, cdna, cds, residues);
        let a = removed_start_row(candidate, cdna, cds, residues);
        assert!(removed_start_pair(&e, &a), "{reference} {cdna}");
        assert!(!removed_start_pair(&a, &e), "{reference} {cdna}");
        assert!(!removed_start_pair(&e, &e), "{reference} {cdna}");
    }
}

/// A row off the shape is no pair of the class: a row of another feature type, terms that
/// differ otherwise than by the candidate's `start_lost`, no `inframe_deletion`, a
/// position that does not begin at 1 and end at a number or differs between the sides,
/// and residues with an alternate half, without one, or that differ between the sides.
#[test]
fn a_row_off_the_shape_of_the_removed_start_is_refused() {
    let e = removed_start_row(START_KEPT, "1-3", "1-3", "M/-");
    let a = removed_start_row(START_LOST, "1-3", "1-3", "M/-");
    assert!(removed_start_pair(&e, &a));
    let regulatory = |entry: &Entry| Entry {
        key: entry.key.replace("Transcript", "RegulatoryFeature"),
        ..entry.clone()
    };
    assert!(!removed_start_pair(&regulatory(&e), &regulatory(&a)));
    for (reference, candidate) in [
        (
            START_KEPT,
            "start_lost,inframe_deletion,splice_region_variant",
        ),
        ("inframe_deletion,splice_region_variant", START_LOST),
        (START_KEPT, "stop_lost,inframe_deletion"),
        (START_KEPT, "start_lost"),
        (START_LOST, START_LOST),
        (
            "inframe_deletion,start_retained_variant",
            "start_lost,inframe_deletion,start_retained_variant",
        ),
        (
            "coding_sequence_variant",
            "start_lost,coding_sequence_variant",
        ),
        ("frameshift_variant", "start_lost,frameshift_variant"),
    ] {
        let e = removed_start_row(reference, "1-3", "1-3", "M/-");
        let a = removed_start_row(candidate, "1-3", "1-3", "M/-");
        assert!(
            !removed_start_pair(&e, &a),
            "{reference} against {candidate}"
        );
    }
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
    ] {
        let e = removed_start_row(START_KEPT, cdna, cds, "M/-");
        let a = removed_start_row(START_LOST, cdna, cds, "M/-");
        assert!(!removed_start_pair(&e, &a), "{cdna} {cds}");
    }
    assert!(!removed_start_pair(
        &e,
        &removed_start_row(START_LOST, "1-6", "1-3", "M/-")
    ));
    assert!(!removed_start_pair(
        &e,
        &removed_start_row(START_LOST, "1-3", "1-6", "M/-")
    ));
    for residues in ["MA/A", "M/L", "M", "-", "-/-"] {
        let e = removed_start_row(START_KEPT, "1-3", "1-3", residues);
        let a = removed_start_row(START_LOST, "1-3", "1-3", residues);
        assert!(!removed_start_pair(&e, &a), "{residues}");
    }
    assert!(!removed_start_pair(
        &e,
        &removed_start_row(START_LOST, "1-3", "1-3", "L/-")
    ));
}

/// On a row of the removed start the IMPACT follows the two sets of terms, and no field of
/// Downstream, NMD or LOFTEE does on either side: `start_lost` is no term their plugins
/// run on (Downstream.pm:97-99, NMD.pm:65 and :83, LoF.pm:166-167).
#[test]
fn the_impact_follows_the_terms_of_a_removed_start_and_no_plugin_field_does() {
    let e = removed_start_row(START_KEPT, "1-3", "1-3", "M/-");
    let a = removed_start_row(START_LOST, "1-3", "1-3", "M/-");
    assert!(removed_start_pair(&e, &a));
    for field in IMPACT_FIELDS {
        assert!(
            term_pair_following_cell(&e, &a, field, "MODERATE", "HIGH"),
            "{field}"
        );
        assert!(
            !term_pair_following_cell(&e, &a, field, "HIGH", "MODERATE"),
            "{field}"
        );
        assert!(
            !term_pair_following_cell(&e, &a, field, "LOW", "HIGH"),
            "{field}"
        );
    }
    for field in TRUNCATING_TERM_FIELDS.iter().chain(&DOWNSTREAM_FIELDS) {
        assert!(
            !term_pair_following_cell(&e, &a, field, "-", "HC"),
            "{field}"
        );
        assert!(
            !term_pair_following_cell(&e, &a, field, "HC", "-"),
            "{field}"
        );
    }
}

/// The five made-up rows of `frameshift_inside_a_held_stop_codon`, an insertion or a
/// deletion inside the codon of the stop a peptide holds at residue 10 and reads on past:
/// the Location with its allele, the reference's terms, and the residues and the Codons
/// cell both sides write.
const HELD_STOP_ROWS: [([&str; 2], &str, &str, &str); 5] = [
    (
        ["1:134-135", "T"],
        "coding_sequence_variant",
        "*/LX",
        "tga/tTga",
    ),
    (["1:136", "-"], "coding_sequence_variant", "*/X", "tgA/tg"),
    (["1:135-139", "-"], "inframe_deletion", "*A/X", "tGAGCT/t"),
    (
        ["1:135-136", "AATG"],
        "inframe_insertion,coding_sequence_variant",
        "*/*MX",
        "tga/tgAATGa",
    ),
    (["1:135", "-"], "coding_sequence_variant", "*/X", "tGa/ta"),
];

/// One Transcript row at `place`, a Location with its allele: its terms, its residues and
/// its Codons cell. The three position columns are empty, since the rule reads none.
fn held_stop_row(place: [&str; 2], terms: &str, residues: &str, codons: &str) -> Entry {
    let [location, allele] = place;
    let text = format!(
        "{HEADER}v\t{location}\t{allele}\tENSG1\tENST1\tTranscript\t{terms}\t-\t-\t-\t{residues}\t{codons}\t-\t-\n"
    );
    let rows = parse_default(&text).expect("a row of fourteen columns");
    Entry::clone(&rows[0])
}

fn held_stop_pair(e: &Entry, a: &Entry) -> bool {
    super::frameshift_inside_a_held_stop_codon_pair(e, a).expect("no refusal")
}

/// The rule of `frameshift_inside_a_held_stop_codon` takes each of the five made-up rows,
/// one for each group of terms the reference writes and two whose new frame keeps the stop,
/// alone and beside a term both sides write, and refuses each with the sides exchanged.
#[test]
fn a_frameshift_inside_a_held_stop_codon_is_a_pair_of_the_class() {
    for (place, reference, residues, codons) in HELD_STOP_ROWS {
        let e = held_stop_row(place, reference, residues, codons);
        let a = held_stop_row(place, "frameshift_variant", residues, codons);
        assert!(held_stop_pair(&e, &a), "{codons}");
        assert!(!held_stop_pair(&a, &e), "{codons}");
        for kept in ["NMD_transcript_variant", "splice_region_variant"] {
            let e = held_stop_row(place, &format!("{reference},{kept}"), residues, codons);
            let a = held_stop_row(
                place,
                &format!("frameshift_variant,{kept}"),
                residues,
                codons,
            );
            assert!(held_stop_pair(&e, &a), "{codons} {kept}");
            assert!(!held_stop_pair(&a, &e), "{codons} {kept}");
        }
    }
}

/// The two made-up controls are no pair of the class: the insertion written one base before
/// the stop codon, where both sides write `frameshift_variant`, and a base lost from the
/// closing stop codon, where both write `stop_lost`.
#[test]
fn the_controls_beside_a_held_stop_codon_are_no_pair_of_the_class() {
    let before = held_stop_row(["1:133-134", "T"], "frameshift_variant", "-/X", "-/T");
    assert!(!held_stop_pair(&before, &before));
    let closing = held_stop_row(["1:178", "-"], "stop_lost", "*/X", "taG/ta");
    assert!(!held_stop_pair(&closing, &closing));
}

/// A row off the shape is no pair of the class: terms that differ in another way, a
/// reference that holds `stop_retained_variant` or `stop_lost`, residues or codons that
/// differ between the sides or lack the shape, and a row of another feature type. A part of
/// the Codons cell written `-` holds no base.
#[test]
fn a_row_off_the_shape_of_the_held_stop_is_refused() {
    let place = ["1:134-135", "T"];
    let (residues, codons) = ("*/LX", "tga/tTga");
    let (unshifted, shifted) = ("coding_sequence_variant", "frameshift_variant");
    let row =
        |terms: &str, residues: &str, codons: &str| held_stop_row(place, terms, residues, codons);
    let of = |terms: &str| row(terms, residues, codons);
    let e = of(unshifted);
    let a = of(shifted);
    assert!(held_stop_pair(&e, &a));
    for (what, reference, candidate) in [
        ("frameshift_variant on both sides", of(shifted), a.clone()),
        (
            "no frameshift_variant in the candidate",
            e.clone(),
            of("protein_altering_variant"),
        ),
        (
            "the reference's term kept",
            e.clone(),
            of("frameshift_variant,coding_sequence_variant"),
        ),
        (
            "inframe_insertion without coding_sequence_variant",
            of("inframe_insertion"),
            a.clone(),
        ),
        (
            "inframe_deletion with coding_sequence_variant",
            of("inframe_deletion,coding_sequence_variant"),
            a.clone(),
        ),
        (
            "a term of no group",
            of("protein_altering_variant"),
            a.clone(),
        ),
        (
            "a further term of the reference alone",
            of("coding_sequence_variant,splice_region_variant"),
            a.clone(),
        ),
        (
            "a further term of the candidate alone",
            e.clone(),
            of("frameshift_variant,NMD_transcript_variant"),
        ),
        (
            "stop_retained_variant on both sides",
            of("coding_sequence_variant,stop_retained_variant"),
            of("frameshift_variant,stop_retained_variant"),
        ),
        (
            "stop_lost on both sides",
            of("coding_sequence_variant,stop_lost"),
            of("frameshift_variant,stop_lost"),
        ),
        (
            "residues that differ between the sides",
            e.clone(),
            row(shifted, "*/LXX", codons),
        ),
        (
            "codons that differ between the sides",
            e.clone(),
            row(shifted, residues, "tga/tAga"),
        ),
    ] {
        assert!(!held_stop_pair(&reference, &candidate), "{what}");
    }
    for (what, residues, codons) in [
        ("a first reference residue that is no stop", "W/LX", codons),
        ("a stop after the first reference residue", "A*/LX", codons),
        ("alternate residues without a partial codon", "*/L", codons),
        ("residues of one part", "*", codons),
        ("no residues", "-", codons),
        ("a length change of three bases", residues, "tga/tTTTga"),
        ("no length change", residues, "tga/tTa"),
        ("an alternate part that holds no base", residues, "tga/-"),
        ("codons of one part", residues, "tga"),
        ("no codons", residues, "-"),
    ] {
        let e = row(unshifted, residues, codons);
        let a = row(shifted, residues, codons);
        assert!(!held_stop_pair(&e, &a), "{what}");
    }
    let regulatory = |entry: &Entry| Entry {
        key: entry.key.replace("Transcript", "RegulatoryFeature"),
        ..entry.clone()
    };
    assert!(!held_stop_pair(&regulatory(&e), &regulatory(&a)));
}

/// On a pair of the class a field of Downstream, NMD or LOFTEE follows the terms where the
/// candidate alone writes it, and the IMPACT is that of each side's most severe term. A
/// field the reference writes, a field no term decides and, beside a splice term both sides
/// write, a Downstream field are not taken.
#[test]
fn a_plugin_field_follows_the_terms_of_a_frameshift_inside_a_held_stop_codon() {
    for (place, reference, residues, codons) in HELD_STOP_ROWS {
        let e = held_stop_row(place, reference, residues, codons);
        let a = held_stop_row(place, "frameshift_variant", residues, codons);
        assert!(held_stop_pair(&e, &a), "{codons}");
        for field in DOWNSTREAM_FIELDS.iter().chain(&TRUNCATING_TERM_FIELDS) {
            assert!(
                term_pair_following_cell(&e, &a, field, "-", "HC"),
                "{codons} {field}"
            );
            assert!(
                !term_pair_following_cell(&e, &a, field, "HC", "-"),
                "{codons} {field}"
            );
            assert!(
                !term_pair_following_cell(&e, &a, field, "HC", "LC"),
                "{codons} {field}"
            );
        }
        let impact = if reference == "coding_sequence_variant" {
            "MODIFIER"
        } else {
            "MODERATE"
        };
        assert!(
            term_pair_following_cell(&e, &a, "IMPACT", impact, "HIGH"),
            "{codons}"
        );
        assert!(
            !term_pair_following_cell(&e, &a, "IMPACT", "LOW", "HIGH"),
            "{codons}"
        );
        for field in ["REVEL", "HGVS_OFFSET"] {
            assert!(
                !term_pair_following_cell(&e, &a, field, "-", "1"),
                "{codons} {field}"
            );
        }
    }
    let (place, reference, residues, codons) = HELD_STOP_ROWS[0];
    let e = held_stop_row(
        place,
        &format!("{reference},splice_region_variant"),
        residues,
        codons,
    );
    let a = held_stop_row(
        place,
        "frameshift_variant,splice_region_variant",
        residues,
        codons,
    );
    assert!(held_stop_pair(&e, &a));
    for field in DOWNSTREAM_FIELDS {
        assert!(
            !term_pair_following_cell(&e, &a, field, "-", "MK"),
            "{field}"
        );
    }
    for field in TRUNCATING_TERM_FIELDS {
        assert!(
            term_pair_following_cell(&e, &a, field, "-", "HC"),
            "{field}"
        );
    }
}

/// One Transcript row of an insertion or deletion at the start codon of a transcript whose
/// first base is 21:100: its Location, its allele, its terms and its two position columns.
fn kept_start_row(location: &str, allele: &str, terms: &str, cdna: &str, cds: &str) -> Entry {
    let text = format!(
        "{HEADER}v\t{location}\t{allele}\tENSG1\tENST1\tTranscript\t{terms}\t{cdna}\t{cds}\t1\t-\t-\t-\t-\n"
    );
    let rows = parse_default(&text).expect("a row of fourteen columns");
    Entry::clone(&rows[0])
}

/// A row of the one-base deletion at 21:102, the third base of the start codon.
fn kept_start_deletion(terms: &str, cdna: &str, cds: &str) -> Entry {
    kept_start_row("21:102", "-", terms, cdna, cds)
}

fn kept_start_pair(e: &Entry, a: &Entry) -> bool {
    super::start_codon_kept_without_five_prime_utr_pair(e, a).expect("no refusal")
}

const START_LOST_SHIFT: &str = "frameshift_variant,start_lost";
const START_KEPT_SHIFT: &str = "frameshift_variant,start_retained_variant";

/// The rule of `start_codon_kept_without_five_prime_utr` takes the rows of the tests of
/// `tools/vep-formats-compare`: a frameshift deletion and a frameshift insertion whose
/// reference holds `start_lost`, a deletion of whole codons whose reference holds it, an
/// insertion of whole codons whose reference holds neither start term, a deletion across an
/// intron with the candidate's terms in the order vep-rs writes them, and a row under
/// `--total_length`. It refuses each with the sides exchanged.
#[test]
fn a_start_codon_kept_without_a_five_prime_utr_is_a_pair_of_the_class() {
    for (location, allele, reference, candidate, cdna, cds) in [
        ("21:102", "-", START_LOST_SHIFT, START_KEPT_SHIFT, "3", "3"),
        (
            "21:101-102",
            "GA",
            "frameshift_variant,start_lost,NMD_transcript_variant",
            "frameshift_variant,start_retained_variant,NMD_transcript_variant",
            "2-3",
            "2-3",
        ),
        (
            "21:102-104",
            "-",
            "start_lost,inframe_deletion",
            "inframe_deletion,start_retained_variant",
            "3-5",
            "3-5",
        ),
        (
            "21:101-102",
            "GGA",
            "inframe_insertion",
            "inframe_insertion,start_retained_variant",
            "2-3",
            "2-3",
        ),
        (
            "21:100-420",
            "-",
            "splice_acceptor_variant,splice_donor_variant,frameshift_variant,start_lost,splice_donor_5th_base_variant,intron_variant",
            "splice_acceptor_variant,splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant,start_retained_variant,intron_variant",
            "1-10",
            "1-10",
        ),
        (
            "21:102-104",
            "-",
            "start_lost,inframe_deletion",
            "inframe_deletion,start_retained_variant",
            "3-5/1500",
            "3-5/900",
        ),
        (
            "21:100-102",
            "-",
            "inframe_deletion",
            "inframe_deletion,start_retained_variant",
            "1-3",
            "1-3",
        ),
    ] {
        let e = kept_start_row(location, allele, reference, cdna, cds);
        let a = kept_start_row(location, allele, candidate, cdna, cds);
        assert!(kept_start_pair(&e, &a), "{reference} {cdna}");
        assert!(!kept_start_pair(&a, &e), "{reference} {cdna}");
        assert!(!kept_start_pair(&e, &e), "{reference} {cdna}");
        assert!(!kept_start_pair(&a, &a), "{reference} {cdna}");
    }
}

/// A row off the shape is no pair of the class: a row of another feature type, terms that
/// differ otherwise than by the two start terms, a candidate that keeps `start_lost`, a
/// reference that holds `start_retained_variant`, a row without `frameshift_variant`,
/// `inframe_insertion` and `inframe_deletion`, a row behind a 5' UTR, a span that opens
/// past the third coding base, a position that is no number and positions that differ
/// between the sides.
#[test]
fn a_row_off_the_shape_of_the_kept_start_is_refused() {
    let e = kept_start_deletion(START_LOST_SHIFT, "3", "3");
    let a = kept_start_deletion(START_KEPT_SHIFT, "3", "3");
    assert!(kept_start_pair(&e, &a));
    let regulatory = |entry: &Entry| Entry {
        key: entry.key.replace("Transcript", "RegulatoryFeature"),
        ..entry.clone()
    };
    assert!(!kept_start_pair(&regulatory(&e), &regulatory(&a)));
    for (reference, candidate) in [
        (
            START_LOST_SHIFT,
            "frameshift_variant,start_retained_variant,splice_region_variant",
        ),
        (
            "frameshift_variant,start_lost,splice_region_variant",
            START_KEPT_SHIFT,
        ),
        ("frameshift_variant,stop_lost", START_KEPT_SHIFT),
        (START_LOST_SHIFT, "frameshift_variant,stop_retained_variant"),
        (
            START_LOST_SHIFT,
            "frameshift_variant,start_lost,start_retained_variant",
        ),
        (START_LOST_SHIFT, "frameshift_variant"),
        (
            "frameshift_variant,start_lost,start_retained_variant",
            START_KEPT_SHIFT,
        ),
        ("start_lost", "start_retained_variant"),
        (
            "start_lost,protein_altering_variant",
            "protein_altering_variant,start_retained_variant",
        ),
        (
            "start_lost,coding_sequence_variant",
            "coding_sequence_variant,start_retained_variant",
        ),
        ("", "start_retained_variant"),
    ] {
        let e = kept_start_deletion(reference, "3", "3");
        let a = kept_start_deletion(candidate, "3", "3");
        assert!(!kept_start_pair(&e, &a), "{reference} against {candidate}");
    }
    for (cdna, cds) in [
        ("103", "3"),
        ("4", "3"),
        ("3", "4"),
        ("4", "4"),
        ("4-6", "4-6"),
        ("10-12", "10-12"),
        ("103-105", "3-5"),
        ("3-5", "3-6"),
        ("?-5", "?-5"),
        ("3-?", "3-?"),
        ("-", "3"),
        ("3", "-"),
        ("-", "-"),
        ("0-2", "0-2"),
        ("03", "03"),
    ] {
        let e = kept_start_deletion(START_LOST_SHIFT, cdna, cds);
        let a = kept_start_deletion(START_KEPT_SHIFT, cdna, cds);
        assert!(!kept_start_pair(&e, &a), "{cdna} {cds}");
    }
    assert!(!kept_start_pair(
        &e,
        &kept_start_deletion(START_KEPT_SHIFT, "2", "2")
    ));
    assert!(!kept_start_pair(
        &e,
        &kept_start_deletion(START_KEPT_SHIFT, "3-4", "3-4")
    ));
    // Each of the first three coding bases opens a span of the class, and the fourth none.
    for (first, taken) in [("1", true), ("2", true), ("3", true), ("4", false)] {
        let e = kept_start_deletion(START_LOST_SHIFT, first, first);
        let a = kept_start_deletion(START_KEPT_SHIFT, first, first);
        assert_eq!(kept_start_pair(&e, &a), taken, "{first}");
    }
}

/// On a row of the kept start the IMPACT follows the two sets of terms where they give two
/// impacts, and no field of Downstream, NMD or LOFTEE follows them on either side: the two
/// sets differ by the start terms alone, which no plugin's condition reads
/// (Downstream.pm:97-99, NMD.pm:65 and :83, LoF.pm:166-167).
#[test]
fn the_impact_follows_the_terms_of_a_kept_start_and_no_plugin_field_does() {
    let e = kept_start_row(
        "21:102-104",
        "-",
        "start_lost,inframe_deletion",
        "3-5",
        "3-5",
    );
    let a = kept_start_row(
        "21:102-104",
        "-",
        "inframe_deletion,start_retained_variant",
        "3-5",
        "3-5",
    );
    assert!(kept_start_pair(&e, &a));
    for field in IMPACT_FIELDS {
        assert!(
            term_pair_following_cell(&e, &a, field, "HIGH", "MODERATE"),
            "{field}"
        );
        assert!(
            !term_pair_following_cell(&e, &a, field, "MODERATE", "HIGH"),
            "{field}"
        );
        assert!(
            !term_pair_following_cell(&e, &a, field, "HIGH", "LOW"),
            "{field}"
        );
    }
    let shifted = kept_start_deletion(START_LOST_SHIFT, "3", "3");
    let kept = kept_start_deletion(START_KEPT_SHIFT, "3", "3");
    assert!(kept_start_pair(&shifted, &kept));
    for field in IMPACT_FIELDS {
        assert!(
            !term_pair_following_cell(&shifted, &kept, field, "HIGH", "LOW"),
            "{field}"
        );
    }
    for (e, a) in [(&e, &a), (&shifted, &kept)] {
        for field in TRUNCATING_TERM_FIELDS.iter().chain(&DOWNSTREAM_FIELDS) {
            assert!(!term_pair_following_cell(e, a, field, "-", "HC"), "{field}");
            assert!(!term_pair_following_cell(e, a, field, "HC", "-"), "{field}");
        }
    }
}
