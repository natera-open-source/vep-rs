// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The splice and start pair rules take an addition and refuse a substitution. Of the writes
//! in `_intron_effects` only `splice_region` is a plain assignment a later region can undo
//! (BaseTranscriptVariationAllele.pm:215); the donor and acceptor flags are set once and stay,
//! so a pair in which the reference holds one of those terms and vep-rs does not is a vep-rs
//! downgrade, which no rule sets aside.

use super::*;
use crate::input::KeyPrep;

/// The consequence set as the comparison keys it: its terms unique, sorted and comma-joined.
fn norm(csq: &str) -> String {
    let mut set = Vec::new();
    KeyPrep::default().normalized_set(csq, &mut set);
    String::from_utf8(set).unwrap()
}

#[test]
fn a_donor_term_exchanged_for_splice_region_is_not_a_splice_overcall() {
    assert!(!is_splice_overcall_swap(
        &norm("missense_variant,splice_donor_variant"),
        &norm("missense_variant,splice_region_variant"),
    ));
}

#[test]
fn an_acceptor_term_exchanged_for_splice_region_is_not_a_splice_overcall() {
    assert!(!is_splice_overcall_swap(
        &norm("intron_variant,splice_acceptor_variant"),
        &norm("intron_variant,splice_region_variant"),
    ));
}

#[test]
fn an_acceptor_term_exchanged_for_splice_region_is_not_a_frameshift_cascade() {
    assert!(!is_frameshift_splice_cascade_swap(
        &norm("frameshift_variant,splice_acceptor_variant"),
        &norm("frameshift_variant,splice_region_variant"),
    ));
}

#[test]
fn intron_variant_added_alone_is_not_a_frameshift_cascade() {
    assert!(!is_frameshift_splice_cascade_swap(
        &norm("frameshift_variant"),
        &norm("frameshift_variant,intron_variant"),
    ));
}

#[test]
fn splice_terms_added_to_a_set_the_reference_lost_nothing_of_are_a_splice_overcall() {
    assert!(is_splice_overcall_swap(
        &norm("frameshift_variant,stop_gained"),
        &norm("frameshift_variant,intron_variant,splice_acceptor_variant,stop_gained"),
    ));
}

#[test]
fn a_start_term_dropped_beside_a_term_vep_rs_adds_is_not_the_start_co_emission() {
    assert!(!is_start_cooccurrence_swap(
        &norm("start_lost,start_retained_variant"),
        &norm("start_retained_variant,frameshift_variant"),
        "-",
    ));
}

#[test]
fn start_lost_dropped_alone_under_a_sequence_allele_is_the_start_co_emission() {
    assert!(is_start_cooccurrence_swap(
        &norm("start_lost,start_retained_variant"),
        &norm("start_retained_variant"),
        "-",
    ));
}
