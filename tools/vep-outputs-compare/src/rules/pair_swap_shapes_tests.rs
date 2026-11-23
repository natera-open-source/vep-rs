// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The swaps each pair rule takes and refuses, and the nesting of the narrower rules inside
//! the splice rule. Classification is first-match, so a narrower rule ahead of the splice
//! rule names the bucket of the pairs both take, and every pair a narrower rule takes has the
//! splice rule's bucket to fall to.

use super::*;
use crate::input::KeyPrep;

/// The consequence set as a tuple key carries it, through the key preparation's own
/// normaliser, so that no test reads an ordering the comparison never produces.
fn norm(csq: &str) -> String {
    let mut set = Vec::new();
    KeyPrep::default().normalized_set(csq, &mut set);
    String::from_utf8(set).expect("utf-8 consequence set")
}

fn ppt(reference: &str, scored: &str) -> bool {
    is_ppt_overcall_swap(&norm(reference), &norm(scored))
}

fn splice(reference: &str, scored: &str) -> bool {
    is_splice_overcall_swap(&norm(reference), &norm(scored))
}

fn start(reference: &str, scored: &str, allele: &str) -> bool {
    is_start_cooccurrence_swap(&norm(reference), &norm(scored), allele)
}

fn cascade(reference: &str, scored: &str) -> bool {
    is_frameshift_splice_cascade_swap(&norm(reference), &norm(scored))
}

fn intronic(reference: &str, scored: &str) -> bool {
    is_intronic_overcall_swap(&norm(reference), &norm(scored))
}

fn covered(reference: &str, scored: &str) -> bool {
    is_covered_splice_region_swap(&norm(reference), &norm(scored))
}

#[test]
fn the_ppt_rule_takes_splice_region_added_to_intron_and_tract() {
    assert!(ppt(
        "intron_variant,splice_polypyrimidine_tract_variant",
        "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
    ));
}

#[test]
fn the_ppt_rule_takes_the_tract_term_added_to_intron_and_splice_region() {
    assert!(ppt(
        "intron_variant,splice_region_variant",
        "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
    ));
}

#[test]
fn the_ppt_rule_takes_the_tract_term_added_to_a_bare_intron() {
    assert!(ppt(
        "intron_variant",
        "intron_variant,splice_polypyrimidine_tract_variant",
    ));
}

#[test]
fn the_ppt_rule_takes_splice_region_added_to_a_bare_intron() {
    assert!(ppt(
        "intron_variant",
        "intron_variant,splice_region_variant"
    ));
}

#[test]
fn the_ppt_rule_takes_both_terms_added_at_once() {
    assert!(ppt(
        "intron_variant",
        "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
    ));
}

#[test]
fn the_ppt_rule_refuses_a_reference_set_without_intron() {
    assert!(!ppt(
        "splice_region_variant",
        "intron_variant,splice_region_variant",
    ));
}

#[test]
fn the_ppt_rule_refuses_an_added_donor_term() {
    assert!(!ppt(
        "intron_variant",
        "intron_variant,splice_donor_variant"
    ));
}

#[test]
fn the_ppt_rule_refuses_a_term_the_reference_alone_holds() {
    assert!(!ppt(
        "intron_variant,splice_region_variant",
        "intron_variant"
    ));
}

#[test]
fn the_ppt_rule_refuses_identical_sets() {
    assert!(!ppt(
        "intron_variant,splice_region_variant",
        "intron_variant,splice_region_variant",
    ));
}

#[test]
fn the_ppt_rule_refuses_a_missense_to_synonymous_swap() {
    assert!(!ppt("missense_variant", "synonymous_variant"));
}

#[test]
fn the_splice_rule_takes_intron_and_acceptor_added_to_frameshift_and_stop_gained() {
    assert!(splice(
        "frameshift_variant,stop_gained",
        "frameshift_variant,intron_variant,splice_acceptor_variant,stop_gained",
    ));
}

#[test]
fn the_splice_rule_takes_the_donor_group_added_to_frameshift() {
    assert!(splice(
        "frameshift_variant",
        "frameshift_variant,intron_variant,splice_donor_5th_base_variant,splice_donor_variant",
    ));
}

#[test]
fn the_splice_rule_takes_the_donor_group_added_to_frameshift_and_stop_gained() {
    assert!(splice(
        "frameshift_variant,stop_gained",
        "frameshift_variant,intron_variant,splice_donor_5th_base_variant,\
         splice_donor_variant,stop_gained",
    ));
}

#[test]
fn the_splice_rule_takes_intron_and_fifth_base_added_to_a_donor() {
    assert!(splice(
        "splice_donor_variant",
        "intron_variant,splice_donor_5th_base_variant,splice_donor_variant",
    ));
}

/// A term downgrade is a difference of vep-rs's, never one the reference's assignment order
/// accounts for.
#[test]
fn the_splice_rule_refuses_frameshift_downgraded_to_coding_sequence() {
    assert!(!splice(
        "frameshift_variant",
        "coding_sequence_variant,intron_variant,splice_acceptor_variant",
    ));
}

#[test]
fn the_splice_rule_refuses_a_missense_to_synonymous_swap() {
    assert!(!splice("missense_variant", "synonymous_variant"));
}

/// On a sequence allele the reference's retention check has found the start codon intact,
/// so `start_lost` is the member vep-rs leaves out.
#[test]
fn on_a_sequence_allele_the_start_rule_takes_the_set_without_start_lost() {
    assert!(start(
        "frameshift_variant,start_lost,start_retained_variant",
        "frameshift_variant,start_retained_variant",
        "-",
    ));
    assert!(start(
        "5_prime_UTR_variant,start_lost,start_retained_variant",
        "5_prime_UTR_variant,start_retained_variant",
        "TGGCGGCGCG",
    ));
}

#[test]
fn on_a_sequence_allele_the_start_rule_refuses_the_set_without_start_retained() {
    assert!(!start(
        "start_lost,start_retained_variant",
        "start_lost",
        "-"
    ));
}

/// On a structural allele the retention check is skipped, so `start_retained_variant` is the
/// member vep-rs leaves out and `start_lost` the one it keeps.
#[test]
fn on_a_structural_allele_the_start_rule_takes_the_set_without_start_retained_alone() {
    assert!(start(
        "feature_truncation,frameshift_variant,start_lost,start_retained_variant",
        "feature_truncation,frameshift_variant,start_lost",
        "deletion",
    ));
    assert!(!start(
        "start_lost,start_retained_variant",
        "start_retained_variant",
        "deletion",
    ));
}

#[test]
fn the_start_rule_refuses_a_reference_set_with_one_start_term() {
    assert!(!start("start_retained_variant", "start_lost", "-"));
}

#[test]
fn an_allele_outside_the_dash_and_the_base_letters_is_structural() {
    for allele in ["-", "A", "TGGCGGCGCG", "acgt", "N"] {
        assert!(!is_structural_allele(allele), "{allele}");
    }
    for allele in [
        "deletion",
        "duplication",
        "tandem_repeat",
        "insertion",
        "inversion",
    ] {
        assert!(is_structural_allele(allele), "{allele}");
    }
}

#[test]
fn every_pair_the_ppt_rule_takes_the_splice_rule_takes_too() {
    for (reference, scored) in [
        (
            "intron_variant,splice_polypyrimidine_tract_variant",
            "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
        ),
        ("intron_variant", "intron_variant,splice_region_variant"),
        (
            "intron_variant",
            "intron_variant,splice_polypyrimidine_tract_variant",
        ),
        (
            "intron_variant,splice_region_variant",
            "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
        ),
    ] {
        assert!(ppt(reference, scored), "{reference} -> {scored}");
        assert!(splice(reference, scored), "{reference} -> {scored}");
    }
}

#[test]
fn every_pair_the_cascade_rule_takes_the_splice_rule_takes_too() {
    for (reference, scored) in [
        (
            "frameshift_variant,stop_gained",
            "frameshift_variant,intron_variant,splice_acceptor_variant,stop_gained",
        ),
        (
            "frameshift_variant",
            "frameshift_variant,intron_variant,splice_donor_5th_base_variant,splice_donor_variant",
        ),
        (
            "frameshift_variant,stop_gained",
            "frameshift_variant,intron_variant,splice_donor_5th_base_variant,\
             splice_donor_variant,stop_gained",
        ),
    ] {
        assert!(cascade(reference, scored), "{reference} -> {scored}");
        assert!(splice(reference, scored), "{reference} -> {scored}");
    }
}

#[test]
fn every_pair_the_intronic_rule_takes_the_splice_rule_takes_too() {
    for (reference, scored) in [
        ("intron_variant", "intron_variant,splice_region_variant"),
        (
            "intron_variant,splice_polypyrimidine_tract_variant",
            "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
        ),
    ] {
        assert!(intronic(reference, scored), "{reference} -> {scored}");
        assert!(splice(reference, scored), "{reference} -> {scored}");
    }
}

#[test]
fn every_pair_of_the_covered_shape_the_splice_rule_takes_too() {
    for (reference, scored) in [
        ("intron_variant", "intron_variant,splice_region_variant"),
        ("missense_variant", "missense_variant,splice_region_variant"),
    ] {
        assert!(covered(reference, scored), "{reference} -> {scored}");
        assert!(splice(reference, scored), "{reference} -> {scored}");
    }
}

#[test]
fn the_covered_shape_holds_a_bare_intron_gaining_splice_region() {
    assert!(covered(
        "intron_variant",
        "intron_variant,splice_region_variant"
    ));
}

#[test]
fn the_covered_shape_holds_an_intron_and_tract_set_gaining_splice_region_alone() {
    assert!(covered(
        "intron_variant,splice_polypyrimidine_tract_variant",
        "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
    ));
}

/// The shape is not one of intronic sets alone.
#[test]
fn the_covered_shape_holds_a_missense_set_gaining_splice_region() {
    assert!(covered(
        "missense_variant",
        "missense_variant,splice_region_variant",
    ));
}

/// A region the reference skips withholds the intron term and `splice_region_variant`
/// together.
#[test]
fn the_covered_shape_admits_intron_beside_splice_region() {
    assert!(covered(
        "missense_variant",
        "intron_variant,missense_variant,splice_region_variant",
    ));
}

/// The reference sets the tract term's flag once and nothing clears it, so its absence there
/// is not the assignment order's doing.
#[test]
fn the_covered_shape_refuses_an_added_tract_term() {
    assert!(!covered(
        "intron_variant",
        "intron_variant,splice_polypyrimidine_tract_variant",
    ));
}

#[test]
fn the_covered_shape_refuses_the_tract_term_added_beside_splice_region() {
    assert!(!covered(
        "intron_variant",
        "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
    ));
}

#[test]
fn the_covered_shape_refuses_an_added_acceptor_term() {
    assert!(!covered(
        "frameshift_variant,stop_gained",
        "frameshift_variant,intron_variant,splice_acceptor_variant,stop_gained",
    ));
}

#[test]
fn the_covered_shape_refuses_an_added_fifth_base_term() {
    assert!(!covered(
        "intron_variant",
        "intron_variant,splice_donor_5th_base_variant",
    ));
}

/// `splice_region_variant` is the mark of the mechanism; the intron term alone is not.
#[test]
fn the_covered_shape_refuses_intron_added_alone() {
    assert!(!covered(
        "missense_variant",
        "intron_variant,missense_variant"
    ));
}

/// A term the reference holds and vep-rs lacks is a downgrade of vep-rs's, whatever else
/// vep-rs adds.
#[test]
fn the_covered_shape_refuses_a_term_the_reference_alone_holds() {
    assert!(!covered(
        "missense_variant,splice_donor_variant",
        "missense_variant,splice_region_variant",
    ));
}

#[test]
fn the_covered_shape_refuses_identical_sets() {
    assert!(!covered(
        "intron_variant,splice_region_variant",
        "intron_variant,splice_region_variant",
    ));
}

#[test]
fn the_cascade_rule_takes_intron_and_acceptor_added_to_frameshift_and_stop_gained() {
    assert!(cascade(
        "frameshift_variant,stop_gained",
        "frameshift_variant,intron_variant,splice_acceptor_variant,stop_gained",
    ));
}

#[test]
fn the_cascade_rule_takes_intron_and_donor_terms_added_to_frameshift() {
    assert!(cascade(
        "frameshift_variant",
        "frameshift_variant,intron_variant,splice_donor_5th_base_variant,splice_donor_variant",
    ));
}

#[test]
fn the_cascade_rule_takes_the_donor_cascade_added_to_frameshift_and_stop_gained() {
    assert!(cascade(
        "frameshift_variant,stop_gained",
        "frameshift_variant,intron_variant,splice_donor_5th_base_variant,\
         splice_donor_variant,stop_gained",
    ));
}

#[test]
fn the_cascade_rule_takes_the_donor_cascade_added_to_an_inframe_deletion() {
    assert!(cascade(
        "inframe_deletion",
        "inframe_deletion,intron_variant,splice_donor_5th_base_variant,splice_donor_variant",
    ));
}

/// A frameshift written as `coding_sequence_variant` is a difference of vep-rs's and stays
/// one.
#[test]
fn the_cascade_rule_refuses_frameshift_downgraded_to_coding_sequence() {
    assert!(!cascade(
        "frameshift_variant,intron_variant,splice_acceptor_variant,\
         splice_donor_5th_base_variant,splice_donor_variant",
        "coding_sequence_variant,intron_variant,splice_acceptor_variant,\
         splice_donor_5th_base_variant,splice_donor_variant",
    ));
}

#[test]
fn the_cascade_rule_refuses_an_inframe_deletion_downgraded_to_coding_sequence() {
    assert!(!cascade(
        "inframe_deletion,intron_variant,splice_donor_variant",
        "coding_sequence_variant,intron_variant,splice_donor_variant",
    ));
}

/// A missense set has no frameshift or in-frame backbone; the splice rule takes the pair.
#[test]
fn the_cascade_rule_refuses_splice_region_added_to_missense() {
    assert!(!cascade(
        "missense_variant",
        "missense_variant,splice_region_variant",
    ));
}

#[test]
fn the_cascade_rule_refuses_sets_without_a_coding_backbone() {
    assert!(!cascade(
        "splice_donor_variant",
        "intron_variant,splice_donor_variant",
    ));
}

#[test]
fn the_intronic_rule_takes_splice_region_added_to_intron_and_tract() {
    assert!(intronic(
        "intron_variant,splice_polypyrimidine_tract_variant",
        "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
    ));
}

#[test]
fn the_intronic_rule_takes_the_tract_term_added_to_intron_and_splice_region() {
    assert!(intronic(
        "intron_variant,splice_region_variant",
        "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
    ));
}

#[test]
fn the_intronic_rule_takes_the_tract_term_added_to_a_bare_intron() {
    assert!(intronic(
        "intron_variant",
        "intron_variant,splice_polypyrimidine_tract_variant",
    ));
}

#[test]
fn the_intronic_rule_takes_splice_region_added_to_a_bare_intron() {
    assert!(intronic(
        "intron_variant",
        "intron_variant,splice_region_variant"
    ));
}

#[test]
fn the_intronic_rule_refuses_an_added_coding_term() {
    assert!(!intronic(
        "intron_variant",
        "intron_variant,frameshift_variant"
    ));
}

#[test]
fn the_intronic_rule_refuses_sets_that_hold_a_utr_term() {
    assert!(!intronic(
        "5_prime_UTR_variant,intron_variant",
        "5_prime_UTR_variant,intron_variant,splice_region_variant",
    ));
}

#[test]
fn the_intronic_rule_refuses_a_reference_set_without_intron() {
    assert!(!intronic(
        "splice_region_variant",
        "intron_variant,splice_region_variant,splice_polypyrimidine_tract_variant",
    ));
}

#[test]
fn the_intronic_rule_refuses_a_term_the_reference_alone_holds() {
    assert!(!intronic(
        "intron_variant,splice_region_variant",
        "intron_variant"
    ));
}
