// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify --field` on the made-up runs of `tests/fixtures/intron_classes/` and
//! `tests/fixtures/maxentscan/`: an `HGVSp` Ensembl VEP writes and vep-rs leaves empty on
//! an insertion at an intron's end or start is an entry of `field_divergences` under the
//! intron class that takes it, a `MaxEntScan_` field on a record that lies in two splice
//! windows is one under `maxentscan_overlapping_windows`, every other differing field is an
//! `unexplained_residual`, and each manifest is written to the bytes of the fixture's
//! `expected.<n>.json`.

mod classify_support;
mod common;
mod recorded_runs;

use recorded_runs::{call, classes, classes_by_record, Entry, Recorded, UNEXPLAINED};
use std::collections::{BTreeMap, BTreeSet};

const END: &str = "hgvsp_intron_end_insertion";
const START: &str = "hgvsp_intron_start_insertion";
const WINDOWS: &str = "maxentscan_overlapping_windows";

/// The fields a run compares for the intron classes: the class reads the row's `HGVSc` and
/// `HGVS_OFFSET` from the compared fields.
const NOTATION_FIELDS: [&str; 3] = ["HGVSc", "HGVSp", "HGVS_OFFSET"];

/// A scratch copy of a run of `tests/fixtures/intron_classes/` whose corpus holds the
/// cache of the transcripts of `tests/fixtures/cache_shards/intron_classes.json`.
fn intron_run(name: &str) -> Recorded {
    let run = Recorded::new(&format!("intron_classes/{name}"));
    run.cache_shard("", "1", "intron_classes.json");
    run
}

/// The classes expected of the `HGVSp` of each record: the record's name and its class.
fn of_records<'a>(cases: &[(&'a str, &'a str)]) -> BTreeMap<(&'a str, &'a str), &'a str> {
    cases
        .iter()
        .map(|(name, class)| ((*name, "HGVSp"), *class))
        .collect()
}

/// The names of the records in their order.
fn names<'a>(cases: &[(&'a str, &'a str)]) -> Vec<&'a str> {
    cases.iter().map(|(name, _)| *name).collect()
}

/// The values vep-rs writes on the entries of a class.
fn candidate_values<'a>(entries: &'a [Entry], class: &str) -> BTreeSet<&'a str> {
    entries
        .iter()
        .filter(|entry| entry.class == class)
        .map(|entry| entry.vep_rs.as_str())
        .collect()
}

/// The one run of a fixture that compares the three notation fields, held to the class of
/// each record's `HGVSp`.
fn one_run(name: &str, cases: &[(&str, &str)], summary: &str) {
    let run = intron_run(name);
    let entries = run.holds(&call(
        &NOTATION_FIELDS,
        (cases.len(), summary),
        "expected.1.json",
    ));
    assert_eq!(
        classes_by_record(&entries, &names(cases)),
        of_records(cases),
        "{name}"
    );
}

#[test]
fn classify_names_the_intron_end_insertion_class() {
    // An `HGVSp` Ensembl VEP writes and vep-rs leaves empty, on an insertion written between
    // an intron's last base and the first coding base of the next exon, is
    // `hgvsp_intron_end_insertion` for each form of the row's `HGVSc`: a copy of the
    // intron's last base, a duplication of two or more intron positions ending on it,
    // numbered from either side of the intron, an insertion of one base or of bases that do
    // not end in a pair that ends an intron, and an insertion of two or more bases that end
    // in either such pair, the pair alone included, on either strand. A row both engines
    // give an `HGVS_OFFSET`, an exonic `HGVSc`, an intron's end before a 5' UTR or a 3' UTR
    // position, two `HGVSc` values that differ, the end of an intron of 13 bases and of 2
    // bases, a transcript the cache does not hold, a vep-rs notation that is filled and
    // differs, a reference without a notation and a Location whose span lies beside no gap
    // stay an `unexplained_residual`.
    let run = intron_run("intron_end");
    // The lower position beside the gap each record is written at, and whether the class
    // takes its row; the last record's row is written over three positions.
    let cases = [
        ("rec_copy", 11199, true),
        ("rec_duplication", 21199, true),
        ("rec_first_half", 231313, true),
        ("rec_one_base", 31199, true),
        ("rec_insertion", 41199, true),
        ("rec_reverse", 211299, true),
        ("rec_major_pair", 51199, true),
        ("rec_minor_pair", 61199, true),
        ("rec_pair_alone", 141199, true),
        ("rec_pair_reverse", 211099, true),
        ("rec_pair_moved", 131199, false),
        ("rec_pair_13_bases", 241112, false),
        ("rec_pair_2_bases", 241201, false),
        ("rec_moved", 71199, false),
        ("rec_exonic", 81249, false),
        ("rec_utr5", 221199, false),
        ("rec_utr3", 221599, false),
        ("rec_two_notations", 91199, false),
        ("rec_13_bases", 231112, false),
        ("rec_2_bases", 231201, false),
        ("rec_absent", 991199, false),
        ("rec_other_protein", 101199, false),
        ("rec_no_reference", 111199, false),
        ("rec_wide", 121199, false),
    ];
    let location = |name: &str, position: i32| match name {
        "rec_wide" => "1:121201-121203".to_string(),
        _ => format!("1:{position}-{}", position + 1),
    };
    let locations: Vec<String> = cases
        .iter()
        .map(|(name, position, _)| location(name, *position))
        .collect();
    let mut expected: BTreeMap<(&str, &str), &str> = cases
        .iter()
        .zip(&locations)
        .map(|((_, _, taken), location)| {
            let class = if *taken { END } else { UNEXPLAINED };
            ((location.as_str(), "HGVSp"), class)
        })
        .collect();
    expected.insert((locations[17].as_str(), "HGVSc"), UNEXPLAINED);
    let summary = "{'unexplained_residual': 15, 'hgvsp_intron_end_insertion': 10}";
    let entries = run.holds(&call(&NOTATION_FIELDS, (25, summary), "expected.1.json"));
    assert_eq!(classes(&entries, |entry| entry.location.as_str()), expected);
    assert_eq!(candidate_values(&entries, END), BTreeSet::from(["-"]));
    // The class reads the row's HGVSc and HGVS_OFFSET from the compared fields: a run that
    // compares HGVSp alone names nothing.
    let alone = call(
        &["HGVSp"],
        (24, "{'unexplained_residual': 24}"),
        "expected.2.json",
    );
    run.holds(&alone);
}

#[test]
fn classify_names_the_intron_start_insertion_class() {
    // An `HGVSp` Ensembl VEP writes and vep-rs leaves empty, on an insertion the 3' rule
    // leaves between a coding exon's last base and the first base of the intron after it,
    // is `hgvsp_intron_start_insertion`: written in that gap or inside the exon and carried
    // to it by one `HGVS_OFFSET` on both rows, on either strand, before an intron of 13
    // bases and of 2 bases, and with that offset on a row whose Location names the gap
    // itself, as a run under `--shift_genomic 1` prints it. A Location with an offset that
    // is at the gap neither as named nor moved, a notation on both rows, a reference without
    // one, an `HGVSc` of another form, two `HGVSc` values that differ, offsets that differ
    // or are zero, a Location inside the exon without an offset, at the intron's end, at an
    // exon end outside the coding region or over a span beside no gap, and a transcript the
    // cache does not hold stay an `unexplained_residual`.
    let run = intron_run("intron_start");
    let cases = [
        ("rec_gap", START),
        ("rec_carried", START),
        ("rec_reverse", START),
        ("rec_reverse_carried", START),
        ("rec_13_bases", START),
        ("rec_2_bases", START),
        ("rec_both_notations", UNEXPLAINED),
        ("rec_no_reference", UNEXPLAINED),
        ("rec_duplication", UNEXPLAINED),
        ("rec_coding_positions", UNEXPLAINED),
        ("rec_two_numbers", UNEXPLAINED),
        ("rec_utr5", UNEXPLAINED),
        ("rec_utr3", UNEXPLAINED),
        ("rec_two_notations", UNEXPLAINED),
        ("rec_offsets_differ", UNEXPLAINED),
        ("rec_zero_offset", UNEXPLAINED),
        ("rec_inside_exon", UNEXPLAINED),
        ("rec_intron_end", UNEXPLAINED),
        ("rec_utr_exon_end", UNEXPLAINED),
        ("rec_absent", UNEXPLAINED),
        ("rec_wide", UNEXPLAINED),
        ("rec_located_at_gap", START),
        ("rec_reverse_located_at_gap", START),
        ("rec_fits_neither", UNEXPLAINED),
    ];
    let mut expected = of_records(&cases);
    expected.insert(("rec_two_notations", "HGVSc"), UNEXPLAINED);
    expected.insert(("rec_offsets_differ", "HGVS_OFFSET"), UNEXPLAINED);
    let summary = "{'unexplained_residual': 18, 'hgvsp_intron_start_insertion': 8}";
    let entries = run.holds(&call(&NOTATION_FIELDS, (26, summary), "expected.1.json"));
    assert_eq!(classes_by_record(&entries, &names(&cases)), expected);
    assert_eq!(candidate_values(&entries, START), BTreeSet::from(["-"]));
    // A run that compares HGVSp alone names nothing.
    let alone = call(
        &["HGVSp"],
        (24, "{'unexplained_residual': 24}"),
        "expected.2.json",
    );
    run.holds(&alone);
}

#[test]
fn classify_names_the_intron_classes_for_an_insertion_held_with_a_shared_base() {
    // A record that carries another allele beside an insertion holds the insertion with the
    // reference base it shares, so its row's Location is that base's one position and its
    // Allele the base followed by the inserted bases. An `HGVSp` Ensembl VEP writes and
    // vep-rs leaves empty on such a row is `hgvsp_intron_start_insertion` where the clipped
    // allele lies after a coding exon's last base and `hgvsp_intron_end_insertion` where it
    // lies before an exon's first coding base, on a reverse-strand transcript included. A
    // span beside no gap of the class and an Allele whose base count does not fit the bases
    // the `HGVSc` inserts stay an `unexplained_residual`.
    one_run(
        "held_with_a_shared_base",
        &[
            ("rec_start_held", START),
            ("rec_start_beside_no_gap", UNEXPLAINED),
            ("rec_start_count", UNEXPLAINED),
            ("rec_end_held", END),
            ("rec_end_beside_no_gap", UNEXPLAINED),
            ("rec_end_count", UNEXPLAINED),
        ],
        "{'hgvsp_intron_start_insertion': 1, 'unexplained_residual': 4, 'hgvsp_intron_end_insertion': 1}",
    );
}

#[test]
fn classify_places_a_row_by_its_location_span_for_the_intron_classes() {
    // A record that carries another allele beside the insertion keeps its own Location, one
    // position for a one-base REF, and either intron class is named by a gap of that span.
    // Under `--minimal` the row's Allele is the insertion alone: the bases an `ins` notation
    // inserts, turned by the offset where the 3' rule moved the record, or as many bases as
    // a `dup` notation duplicates. A held Allele spelled on the reverse strand is read
    // reverse-complemented. A minimised Allele that is not the notation's bases, a `dup` row
    // whose Allele's base count does not fit, a duplication that starts in the intron's
    // first half, a moved Allele that holds the shared base and a span beside no gap stay an
    // `unexplained_residual`.
    one_run(
        "location_span",
        &[
            ("rec_end_minimised", END),
            ("rec_end_minimised_reverse", END),
            ("rec_end_minimised_copy", END),
            ("rec_end_minimised_duplication", END),
            ("rec_end_reverse_strand_record", END),
            ("rec_start_minimised", START),
            ("rec_start_minimised_past_the_gap", START),
            ("rec_start_minimised_moved", START),
            ("rec_start_reverse_strand_record", START),
            ("rec_other_bases", UNEXPLAINED),
            ("rec_count", UNEXPLAINED),
            ("rec_held_moved", UNEXPLAINED),
            ("rec_first_half", UNEXPLAINED),
            ("rec_below_the_gap", UNEXPLAINED),
            ("rec_above_the_gap", UNEXPLAINED),
        ],
        "{'hgvsp_intron_start_insertion': 4, 'unexplained_residual': 6, 'hgvsp_intron_end_insertion': 5}",
    );
}

#[test]
fn classify_places_a_row_of_two_adjacent_positions_by_each_gap_of_its_span() {
    // A record whose reference allele spans two bases and that carries another allele beside
    // the insertion keeps both bases under `--minimal`: the row's Location is their two
    // adjacent positions and its Allele the insertion alone. Either intron class is named
    // by a gap of that span, its own, as with an insertion written alone, the one before its
    // first position or the one after its last, moved by the offset where the 3' rule moved
    // the record. A Location none of whose three gaps is a gap of the class, and an Allele
    // that is not the insertion at a Location whose own gap is none, stay an
    // `unexplained_residual`.
    one_run(
        "two_adjacent_positions",
        &[
            ("rec_end_after_the_span", END),
            ("rec_end_before_the_span", END),
            ("rec_start_after_the_span", START),
            ("rec_start_before_the_span", START),
            ("rec_start_moved", START),
            ("rec_end_no_gap", UNEXPLAINED),
            ("rec_start_no_gap", UNEXPLAINED),
            ("rec_start_other_bases", UNEXPLAINED),
        ],
        "{'hgvsp_intron_start_insertion': 3, 'unexplained_residual': 3, 'hgvsp_intron_end_insertion': 2}",
    );
}

#[test]
fn classify_names_the_intron_start_class_for_an_allele_printed_as_the_rule_moved_it() {
    // A run under `--shift_3prime 1` prints a moved record's Allele as the 3' rule moved it,
    // the bases the row's `HGVSc` inserts as they stand, beside the offset: a minimised row
    // of a record that carries another allele is `hgvsp_intron_start_insertion` with the
    // record's own Location, which the offset carries to the gap, and with the moved
    // Location a run under `--shift_genomic 1` prints, on a forward-strand transcript with a
    // positive offset and on a reverse-strand one with a negative offset. An Allele that is
    // the notation's bases under neither reading and a span whose gaps are a gap of the
    // class neither as named nor moved stay an `unexplained_residual`.
    one_run(
        "allele_as_the_rule_moved_it",
        &[
            ("rec_forward", START),
            ("rec_forward_located_moved", START),
            ("rec_reverse", START),
            ("rec_reverse_located_moved", START),
            ("rec_neither_reading", UNEXPLAINED),
            ("rec_neither_place", UNEXPLAINED),
        ],
        "{'hgvsp_intron_start_insertion': 4, 'unexplained_residual': 2}",
    );
}

#[test]
fn classify_names_maxentscan_overlapping_windows_where_the_cache_places_the_record_in_two_windows()
{
    // `classify --field` names `maxentscan_overlapping_windows` on a `MaxEntScan_` field of
    // a row whose one-base record lies in the acceptor window of one intron of its
    // transcript and the donor window of the next, read from the corpus cache: the three
    // scores where both engines write a number, and under `verbose` the two window texts
    // where Ensembl VEP writes nine bases and vep-rs twenty-three. A row of the same shape
    // whose record lies in one window or in none, on a transcript without such an exon or
    // over several bases, stays an `unexplained_residual`, as does every such row where the
    // corpus holds no cache.
    let run = Recorded::new("maxentscan/two_windows");
    let fields = [
        "MaxEntScan_ref",
        "MaxEntScan_alt",
        "MaxEntScan_diff",
        "MaxEntScan_ref_seq",
        "MaxEntScan_alt_seq",
    ];
    let without_a_cache = "{'unexplained_residual': 30}";
    run.holds(&call(&fields, (30, without_a_cache), "expected.1.json"));
    run.cache_shard("", "1", "maxentscan.json");
    let summary = "{'unexplained_residual': 20, 'maxentscan_overlapping_windows': 10}";
    let entries = run.holds(&call(&fields, (30, summary), "expected.2.json"));
    let taken: Vec<&Entry> = entries
        .iter()
        .filter(|entry| entry.class == WINDOWS)
        .collect();
    let rows: BTreeSet<(&str, &str)> = taken
        .iter()
        .map(|entry| (entry.feature.as_str(), entry.location.as_str()))
        .collect();
    assert_eq!(
        rows,
        BTreeSet::from([("ENST1", "1:202"), ("ENST2", "1:202")])
    );
    let named: BTreeSet<&str> = taken.iter().map(|entry| entry.field.as_str()).collect();
    assert_eq!(named, BTreeSet::from(fields));
    // Under the scores alone the class names the three numbers of each of the two rows.
    let summary = "{'unexplained_residual': 12, 'maxentscan_overlapping_windows': 6}";
    run.holds(&call(&fields[..3], (18, summary), "expected.3.json"));
}
