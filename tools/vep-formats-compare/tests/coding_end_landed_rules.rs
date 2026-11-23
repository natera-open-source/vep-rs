// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The classes of a stop codon or a start codon that stands, beside the candidate's cache
//! and without it: beside the cache each takes its pair only on a transcript that shows
//! its codon, with the allele at it; without the cache each reads the two rows alone.

mod coding_end_common;

use coding_end_common::*;

const NOTHING: &str = "{}";
const REMOVED: &str = "stop_codon_removed_without_stop_codon_test";
const UNTESTED: &str = "stop_retained_without_stop_codon_test";
const THREE_BASE: &str = "stop_retained_three_base_stop_codon_test";
const INSERTION_RANGE: &str = "stop_codon_overlap_insertion_range";
const GENOMIC_SPAN: &str = "stop_codon_overlap_genomic_span";
const UNSHIFTED: &str = "stop_codon_unshifted_span";
const SWAP: &str = "start_cooccurrence_swap";
const NAMED: &str = "hgvsp_stop_named_on_incomplete_coding_sequence";

fn taken(class: &str) -> (String, String) {
    (one(class), GREEN.to_string())
}

fn refused(result: &(String, String)) -> bool {
    result.0 == NOTHING && result.1 != GREEN
}

/// A deletion from the last residue's codon of `closed` past the transcript's last base:
/// it removes the whole stop codon, and no test of the codon is made.
fn stop_codon_removed() -> Record {
    let (theirs, ours) = pair(
        "stop_retained_variant,3_prime_UTR_variant",
        "stop_lost,3_prime_UTR_variant",
        ["14-?", "8-?", "3-?", "-", "-"],
    );
    Record::replacing("removed", &closed(), 14, 30, "-", theirs, ours)
}

/// A deletion from the last base of the stop codon of `closed` past the transcript's last
/// base: no test of the codon is made.
fn stop_codon_untested() -> Record {
    let (theirs, ours) = pair(
        "stop_retained_variant,3_prime_UTR_variant",
        "coding_sequence_variant,3_prime_UTR_variant",
        ["18-?", "12-?", "4-?", "-", "-"],
    );
    Record::replacing("untested", &closed(), 18, 30, "-", theirs, ours)
}

/// A deletion of the middle base of the stop codon of `closed`: the test is made.
fn stop_codon_tested() -> Record {
    let (theirs, ours) = pair(
        "stop_retained_variant",
        "frameshift_variant",
        ["17", "11", "4", "*", "tAa/ta"],
    );
    Record::replacing("tested", &closed(), 17, 17, "-", theirs, ours)
}

/// An insertion before the stop codon of `closed` that keeps it.
fn insertion_that_keeps_the_stop_codon() -> Record {
    let (theirs, ours) = pair(
        "frameshift_variant",
        "stop_retained_variant",
        ["15-16", "9-10", "3-4", "-", "-"],
    );
    Record::inserting("kept-insertion", &closed(), 15, "C", theirs, ours)
}

/// A deletion of the second base of the stop codon of `on`, which ends an exon of
/// `split_stop`.
fn stop_base_at_an_exon_end(on: &Transcript) -> Record {
    let (theirs, ours) = pair(
        "coding_sequence_variant,splice_region_variant",
        "stop_retained_variant,splice_region_variant",
        ["17", "11", "4", "*", "tAa/ta"],
    );
    Record::replacing("exon-end", on, 17, 17, "-", theirs, ours)
}

/// A deletion the 3' shift carries onto the stop codon of `closed`.
fn deletion_shifted_onto_the_stop_codon() -> Record {
    let five = ["15", "9", "3", "-", "-"];
    let theirs =
        cells("coding_sequence_variant,splice_region_variant", five).with("SHIFT_LENGTH", "1");
    let ours = cells("stop_retained_variant,splice_region_variant", five).with("SHIFT_LENGTH", "1");
    Record::replacing("shifted", &closed(), 15, 15, "-", theirs, ours)
}

/// A substitution in the start codon of `closed` for which Ensembl VEP writes both start
/// terms.
fn both_start_terms() -> Record {
    let (theirs, ours) = pair(
        "start_lost,start_retained_variant",
        "start_retained_variant",
        ["8", "2", "1", "M/T", "aTg/aCg"],
    );
    Record::replacing("both-starts", &closed(), 8, 8, "C", theirs, ours)
}

/// The pair of each rule on a transcript whose stop codon stands, and the transcript it is
/// taken on.
fn pairs() -> Vec<(&'static str, Record)> {
    vec![
        (REMOVED, stop_codon_removed()),
        (UNTESTED, stop_codon_untested()),
        (THREE_BASE, stop_codon_tested()),
        (INSERTION_RANGE, insertion_that_keeps_the_stop_codon()),
        (GENOMIC_SPAN, stop_base_at_an_exon_end(&split_stop())),
        (UNSHIFTED, deletion_shifted_onto_the_stop_codon()),
    ]
}

#[test]
fn a_rule_of_a_stop_codon_takes_its_pair_where_the_stop_codon_stands_under_the_allele() {
    for (class, record) in pairs() {
        assert_eq!(
            beside_the_cache(&format!("stands-{class}"), &[record], "default", class),
            taken(class),
            "{class}"
        );
    }
}

/// The same two rows on a transcript where no stop codon stands, with the flag or without,
/// are a wrong row of vep-rs: beside the cache no rule of a stop codon takes them.
#[test]
fn a_rule_of_a_stop_codon_refuses_its_pair_where_no_stop_codon_stands() {
    for (class, record) in pairs() {
        for other in ["ENST_OPEN", "ENST_FLAGGED"] {
            let left = beside_the_cache(
                &format!("none-stands-{class}-{other}"),
                &[record.clone().on(other)],
                "default",
                class,
            );
            assert!(refused(&left), "{class} on {other}: {left:?}");
        }
    }
}

/// A comparison that was given no cache reads each rule as it reads the two rows alone.
#[test]
fn a_rule_of_a_stop_codon_reads_the_rows_alone_without_a_cache() {
    for (class, record) in pairs() {
        for id in [record.rows[0].0, "ENST_OPEN", "ENST_FLAGGED", "ENST_ABSENT"] {
            assert_eq!(
                beside(
                    &format!("no-cache-{class}-{id}"),
                    &[record.clone().on(id)],
                    "default",
                    class,
                    false
                ),
                taken(class),
                "{class} on {id}"
            );
        }
    }
}

/// Beside a cache a rule takes nothing without a reading of the row's transcript: a
/// transcript the cache lacks, a record without a coding sequence, and a stable id whose
/// record at the row's Location is another one.
#[test]
fn landed_rule_beside_a_cache_without_a_verdict() {
    for (class, record) in pairs() {
        for other in ["ENST_ABSENT", "ENST_NONCODING"] {
            let left = beside_the_cache(
                &format!("no-verdict-{class}-{other}"),
                &[record.clone().on(other)],
                "default",
                class,
            );
            assert!(refused(&left), "{class} on {other}: {left:?}");
        }
    }
    // The pair of the removed stop codon on the two records of one stable id: the record
    // at the row's Location shows no stop codon, the one four thousand bases on shows
    // one.
    let [open_record, closed_record] = two_records();
    let on_the_open_record = stop_codon_removed().on("ENST_TWO");
    let left = beside_the_cache(
        "no-verdict-two-records-open",
        std::slice::from_ref(&on_the_open_record),
        "default",
        REMOVED,
    );
    assert!(refused(&left), "{left:?}");
    let on_the_closed_record = on_the_open_record.moved(closed_record.first - open_record.first);
    assert_eq!(
        beside_the_cache(
            "no-verdict-two-records-closed",
            &[on_the_closed_record],
            "default",
            REMOVED
        ),
        taken(REMOVED)
    );
    for other in ["ENST_ABSENT", "ENST_NONCODING"] {
        let left = beside_the_cache(
            &format!("no-verdict-swap-{other}"),
            &[both_start_terms().on(other)],
            "default",
            SWAP,
        );
        assert!(refused(&left), "{other}: {left:?}");
    }
}

/// Where the stop codon stands a rule takes its pair wherever the allele lies: what a rule
/// reads of the allele's place is its own reading of the two rows, beside a cache as
/// without one.
#[test]
fn a_rule_of_a_stop_codon_takes_its_pair_wherever_the_allele_lies_where_the_stop_codon_stands() {
    // The pair of the removed stop codon on a span from the codon's last base on.
    let mut part = stop_codon_removed();
    part.first += 4;
    part.reference = part.reference[4..].to_string();
    // The pair of the untested stop codon on a span that holds the whole codon.
    let whole = stop_codon_untested().moved(-4);
    // The pair of the tested stop codon on a deletion in the codon before it.
    let before = stop_codon_tested().moved(-3);
    // The pair of the stop base at an exon's end on a stop codon no intron splits.
    let unsplit = stop_base_at_an_exon_end(&closed());
    for (class, record) in [
        (REMOVED, part),
        (UNTESTED, whole),
        (THREE_BASE, before),
        (GENOMIC_SPAN, unsplit),
    ] {
        assert_eq!(
            beside_the_cache(
                &format!("allele-place-{class}"),
                &[record],
                "default",
                class
            ),
            taken(class),
            "{class}"
        );
    }
}

/// On a stop codon that polyadenylation completes two rules of a stop codon take nothing:
/// the rule of the three-base test, since the engines read two different triplets there,
/// and the rule of the untested stop codon, which stands where the genome carries the stop
/// codon whole. A span without its coordinates that removes such a stop codon is the pair
/// of the removed stop codon.
#[test]
fn three_base_test_class_on_a_shown_tail() {
    let (theirs, ours) = pair(
        "stop_retained_variant",
        "frameshift_variant",
        ["16", "10", "4", "X/-", "T/-"],
    );
    let tested = Record::replacing("tested-tail", &shown(), 16, 16, "-", theirs, ours);
    let left = beside_the_cache("shown-tail-tested", &[tested], "default", THREE_BASE);
    assert!(refused(&left), "{left:?}");
    let past_the_end = ["16-?", "10-?", "4-?", "-", "-"];
    let (theirs, ours) = pair(
        "stop_retained_variant",
        "coding_sequence_variant",
        past_the_end,
    );
    let untested = Record::replacing("untested-tail", &shown(), 16, 20, "-", theirs, ours);
    let left = beside_the_cache("shown-tail-untested", &[untested], "default", UNTESTED);
    assert!(refused(&left), "{left:?}");
    let (theirs, ours) = pair("stop_retained_variant", "stop_lost", past_the_end);
    let removed = Record::replacing("removed-tail", &shown(), 16, 20, "-", theirs, ours);
    assert_eq!(
        beside_the_cache("shown-tail-removed", &[removed], "default", REMOVED),
        taken(REMOVED)
    );
}

/// A pair of a coding end where no stop codon stands has the terms of a rule of a stop
/// codon that stands: beside the cache the class of the coding end takes it and the rule
/// does not, and without the cache the rule reads it by its terms.
#[test]
fn a_pair_of_a_coding_end_is_its_class_beside_the_cache_and_a_landed_rules_without_it() {
    for (kind, record, landed) in [
        ("past", deletion_past_the_end(&open()), UNTESTED),
        ("kept", kept_stop_deletion(&open()), THREE_BASE),
    ] {
        let both = format!("{OVERLAP},{landed}");
        assert_eq!(
            beside_the_cache(
                &format!("coding-end-{kind}"),
                std::slice::from_ref(&record),
                "default",
                &both
            ),
            taken(OVERLAP),
            "{kind}"
        );
        let left = beside_the_cache(
            &format!("coding-end-{kind}-landed"),
            std::slice::from_ref(&record),
            "default",
            landed,
        );
        assert!(refused(&left), "{kind}: {left:?}");
        assert_eq!(
            beside(
                &format!("coding-end-{kind}-no-cache"),
                &[record],
                "default",
                landed,
                false
            ),
            taken(landed),
            "{kind}"
        );
    }
}

#[test]
fn the_pair_of_both_start_terms_is_taken_only_where_a_start_codon_stands() {
    assert_eq!(
        beside_the_cache("swap-shown", &[both_start_terms()], "default", SWAP),
        taken(SWAP)
    );
    for other in ["ENST_START_OPEN", "ENST_START_FLAGGED"] {
        let left = beside_the_cache(
            &format!("swap-{other}"),
            &[both_start_terms().on(other)],
            "default",
            SWAP,
        );
        assert!(refused(&left), "{other}: {left:?}");
        assert_eq!(
            beside(
                &format!("swap-no-cache-{other}"),
                &[both_start_terms().on(other)],
                "default",
                SWAP,
                false
            ),
            taken(SWAP),
            "{other}"
        );
    }
    // Both classes named: on a transcript that shows no start codon the pair in which
    // vep-rs keeps a start term is neither's.
    let left = beside_the_cache(
        "swap-and-start",
        &[both_start_terms().on("ENST_START_OPEN")],
        "default",
        &format!("{SWAP},{START}"),
    );
    assert!(refused(&left), "{left:?}");
}

/// On a stop codon that polyadenylation completes vep-rs owes a protein notation: a row
/// on which it writes none where Ensembl VEP names the stop alone is a difference beside
/// the cache, and the class's as the row reads without it.
#[test]
fn no_protein_notation_on_a_completed_stop_codon() {
    let record = completed_stop_lost(&shown());
    let (id, theirs, ours) = record.rows[0].clone();
    let theirs = theirs.with("HGVSp", "ENSP_SHOWN.1:p.Ter4%3D");
    let mut none = ours;
    none.more.clear();
    let written = Record {
        rows: vec![(id, theirs, none)],
        ..record
    };
    let both = format!("{COMPLETED},{NAMED}");
    let (applied, reasons) = beside_the_cache(
        "named-completed",
        std::slice::from_ref(&written),
        "tab",
        &both,
    );
    assert_eq!(applied, one(COMPLETED));
    assert_eq!(reasons, "[\"field HGVSp: 1 differing entries\"]");
    let (applied, _) = beside("named-completed-no-cache", &[written], "tab", NAMED, false);
    assert_eq!(applied, one(NAMED));
}
