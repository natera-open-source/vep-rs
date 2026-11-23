// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The invariants of the registry, under which classifying and excluding are independent
//! roles: a pair an excluding rule takes always has a bucket to be counted in, no pair is
//! excluded twice, and every rule names one bucket and one taxonomy class.

use super::table::DIVERGENCE_CLASSES_116_2;
use super::*;
use crate::input::KeyPrep;
use std::collections::HashSet;

/// The consequence set as a tuple key carries it, through the key preparation's own
/// normaliser.
fn norm(csq: &str) -> String {
    let mut set = Vec::new();
    KeyPrep::default().normalized_set(csq, &mut set);
    String::from_utf8(set).expect("utf-8 consequence set")
}

/// A pair excluded without a bucket would leave the per-bucket excluded counts short of the
/// excluded total.
#[test]
fn a_pair_an_excluding_rule_takes_has_a_bucket() {
    for (reference, scored) in [
        ("intron_variant", "intron_variant,splice_region_variant"),
        (
            "intron_variant,splice_polypyrimidine_tract_variant",
            "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant",
        ),
        ("missense_variant", "missense_variant,splice_region_variant"),
        (
            "missense_variant",
            "intron_variant,missense_variant,splice_region_variant",
        ),
        (
            "start_lost,start_retained_variant",
            "start_retained_variant",
        ),
        (
            "5_prime_UTR_variant,start_lost,start_retained_variant",
            "5_prime_UTR_variant,start_retained_variant",
        ),
    ] {
        let (pn, rn) = (norm(reference), norm(scored));
        for rule in EXCLUSION_REGISTRY.iter().filter(|rule| rule.excludes) {
            if !rule.matches(&pn, &rn, "-") {
                continue;
            }
            assert!(
                EXCLUSION_REGISTRY
                    .iter()
                    .filter(|rule| rule.classifies)
                    .any(|rule| rule.matches(&pn, &rn, "-")),
                "{} takes {reference} -> {scored} and no classifying rule does",
                rule.bucket
            );
        }
    }
}

/// The covered shape needs a term vep-rs adds and the start rule needs none, so the excluded
/// count holds no pair twice.
#[test]
fn no_pair_is_taken_by_two_excluding_rules() {
    for (reference, scored) in [
        (
            "start_lost,start_retained_variant",
            "start_retained_variant",
        ),
        ("intron_variant", "intron_variant,splice_region_variant"),
    ] {
        let (pn, rn) = (norm(reference), norm(scored));
        let taking: Vec<&str> = EXCLUSION_REGISTRY
            .iter()
            .filter(|rule| rule.excludes && rule.matches(&pn, &rn, "-"))
            .map(|rule| rule.bucket)
            .collect();
        assert!(taking.len() <= 1, "{reference} -> {scored}: {taking:?}");
    }
}

#[test]
fn the_buckets_of_the_registry_are_unique() {
    let buckets: HashSet<&str> = EXCLUSION_REGISTRY.iter().map(|rule| rule.bucket).collect();
    assert_eq!(buckets.len(), EXCLUSION_REGISTRY.len());
}

#[test]
fn of_the_rules_over_two_consequence_sets_the_start_rule_alone_classifies_and_excludes() {
    let both: Vec<&str> = EXCLUSION_REGISTRY
        .iter()
        .filter(|rule| !rule.row_aware() && rule.classifies && rule.excludes)
        .map(|rule| rule.bucket)
        .collect();
    assert_eq!(both, ["start_cooccurrence_swap"]);
}

/// The taxonomy class is what a report reader maps a bucket onto; an empty one would map it
/// onto nothing.
#[test]
fn every_rule_of_the_registry_has_a_taxonomy_class() {
    for rule in &EXCLUSION_REGISTRY {
        let class = DIVERGENCE_CLASSES_116_2
            .iter()
            .find(|class| class.bucket == rule.bucket)
            .unwrap_or_else(|| panic!("{} has no row in the class table", rule.bucket));
        assert!(!class.taxonomy_class.is_empty(), "{}", rule.bucket);
    }
}
