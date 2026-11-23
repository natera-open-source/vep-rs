// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What each release's collector holds after the merge of two fabricated key streams.

use super::*;
use crate::rules::breakend::{BND_MATE_32_BIT_LIMIT, BND_MATE_32_BIT_LIMIT_DIGITS};
use crate::rules::scope::Group;
use crate::TempDir;
use std::collections::{HashMap, HashSet};

/// One key of a row on a transcript.
fn key(location: &str, allele: &str, feature: &str, consequence: &str) -> String {
    format!("{location}\t{allele}\t{feature}\tTranscript\t{consequence}")
}

/// One key of an intergenic row.
fn intergenic(location: &str, allele: &str) -> String {
    format!("{location}\t{allele}\t-\t-\tintergenic_variant")
}

/// Merges the two sides' keys, sorted, as the file `f`, feeding `collector`.
fn merge(perl: &[String], rust: &[String], collector: Collector) {
    let scratch = TempDir::create(&std::env::temp_dir(), "collector_test_").unwrap();
    let written = |name: &str, keys: &[String]| {
        let mut keys = keys.to_vec();
        keys.sort_unstable();
        let path = scratch.path.join(name);
        let lines: String = keys.iter().map(|key| format!("{key}\n")).collect();
        std::fs::write(&path, lines).unwrap();
        path
    };
    let (perl, rust) = (written("perl.keys", perl), written("rust.keys", rust));
    compare_sorted_key_files(&perl, &rust, "f", None, false, collector).unwrap();
}

/// The 116.2 collector after the merge of the two sides' keys.
fn record_scope(perl: &[String], rust: &[String]) -> ReferenceRecordScope {
    let mut scope = ReferenceRecordScope::default();
    merge(perl, rust, Collector::Release116(&mut scope));
    scope
}

fn names(items: &[&str]) -> HashSet<String> {
    items.iter().map(|item| item.to_string()).collect()
}

fn groups(locations: &[&str]) -> HashSet<Group> {
    locations
        .iter()
        .map(|location| ("f".to_string(), location.to_string()))
        .collect()
}

#[test]
fn the_merge_collects_the_transcripts_each_engine_named_on_a_record_in_the_batch_scope() {
    let giant = "18:62909529-80261528";
    let small = "18:62909529-62919528";
    let perl = [
        key(giant, "deletion", "T1", "intron_variant"),
        key(giant, "deletion", "T2", "intron_variant"),
        key(small, "deletion", "S1", "intron_variant"),
        key("18:100", "A", "TA", "missense_variant"),
    ];
    let rust = [
        key(giant, "deletion", "T1", "intron_variant"),
        key(giant, "deletion", "T3", "transcript_ablation"),
        key(small, "deletion", "S1", "intron_variant"),
        key("18:200-200", "CPX", "C1", "intron_variant"),
    ];
    let mut scope = ScopeFeatures::default();
    merge(&perl, &rust, Collector::Release115(&mut scope));
    let record = |location: &str, allele: &str| (location.to_string(), allele.to_string());
    assert_eq!(
        scope.records,
        HashMap::from([
            (
                record(giant, "deletion"),
                (names(&["T1", "T2"]), names(&["T1", "T3"]))
            ),
            (record("18:200-200", "CPX"), (names(&[]), names(&["C1"]))),
        ]),
        "the giant record and the type without a term are in scope, the ten-kilobase deletion \
         and the single base are not; the transcript both engines named counts on both sides"
    );
}

#[test]
fn a_structural_record_with_no_reference_row_is_absent() {
    let rust = [
        key("21:800-900", "CPX", "C3", "intron_variant"),
        key("21:800-900", "CPX", "C4", "coding_sequence_variant"),
        key("21:30000001", "deletion", "TX_DEL", "feature_truncation"),
        key("21:1", "A", "TX_SNV", "intron_variant"),
    ];
    let perl = [key("21:1", "A", "TX_SNV", "intron_variant")];
    assert_eq!(
        record_scope(&perl, &rust).absent,
        groups(&["21:800-900", "21:30000001"])
    );
}

#[test]
fn a_reference_row_under_another_allele_keeps_the_record_present() {
    let location = "21:30165801-30166000";
    let perl = [key(
        location,
        "copy_number_variation",
        "TX_FIRST",
        "3_prime_UTR_variant",
    )];
    let rust = [
        key(
            location,
            "deletion",
            "TX_FIRST",
            "feature_truncation,3_prime_UTR_variant",
        ),
        key(
            location,
            "duplication",
            "TX_FIRST",
            "feature_elongation,3_prime_UTR_variant",
        ),
    ];
    assert_eq!(record_scope(&perl, &rust).absent, groups(&[]));
}

#[test]
fn sequence_alleles_and_tandem_repeats_never_make_a_record_absent() {
    let rust = [
        key("21:100", "A", "TX_FIRST", "missense_variant"),
        key("21:900-950", "tandem_repeat", "T7", "feature_elongation"),
    ];
    assert_eq!(record_scope(&[], &rust).absent, groups(&[]));
}

#[test]
fn a_mate_position_at_the_buffer_bound_is_recorded_from_the_reference_side() {
    assert_eq!(BND_MATE_32_BIT_LIMIT, 2_147_478_648);
    let perl = [
        intergenic("21:30425651", "N."),
        intergenic("21:30425651", "N[21:2147478648["),
        key("21:40000001", "N[21:2147478647[", "T1", "intron_variant"),
    ];
    let rust = [
        key("21:30425651", "N.", "T2", "upstream_gene_variant"),
        key("21:50000001", "N[21:2147478648[", "T3", "intron_variant"),
    ];
    assert_eq!(
        record_scope(&perl, &rust).mate_overflow,
        groups(&["21:30425651"]),
        "a mate one below the bound and a bracket the scored engine alone wrote are not the shape"
    );
}

#[test]
fn the_last_location_of_a_file_is_closed() {
    let rust = [key("21:999999-1999999", "CPX", "C1", "intron_variant")];
    let perl = [key("21:1", "A", "TX_SNV", "intron_variant")];
    assert_eq!(
        record_scope(&perl, &rust).absent,
        groups(&["21:999999-1999999"])
    );
}

#[test]
fn a_breakend_location_carries_each_sides_alleles_and_the_mate_chromosomes() {
    let perl = [
        key("21:26037668", "N[22:26037668[", "M1", "intron_variant"),
        intergenic("21:11400001", "N."),
        intergenic("21:11400001", "N[22:14000000["),
    ];
    let rust = [
        key("21:26037668", "N[22:26037668[", "M1", "intron_variant"),
        key(
            "21:26037668",
            "chromosome_breakpoint",
            "O1",
            "upstream_gene_variant",
        ),
    ];
    let scope = record_scope(&perl, &rust);
    let at = |location: &str| {
        let breakends = &scope.breakends[&("f".to_string(), location.to_string())];
        (
            breakends.perl_alleles.clone(),
            breakends.rust_alleles.clone(),
            breakends.mates.clone(),
        )
    };
    assert_eq!(
        at("21:26037668"),
        (
            names(&["N[22:26037668["]),
            names(&["N[22:26037668[", "chromosome_breakpoint"]),
            names(&["22"])
        )
    );
    assert_eq!(
        at("21:11400001"),
        (names(&["N.", "N[22:14000000["]), names(&[]), names(&["22"]))
    );
    assert_eq!(
        scope.absent,
        groups(&[]),
        "every breakend record here has a reference row"
    );
}

#[test]
fn the_collector_records_the_locations_only_the_reference_wrote() {
    let perl = [
        intergenic("21:-4-295", "deletion"),
        key("21:100", "A", "T1", "missense_variant"),
        key("21:200", "G", "T2", "intron_variant"),
    ];
    let rust = [
        key("21:100", "A", "T1", "missense_variant"),
        key("21:200", "G", "T2", "upstream_gene_variant"),
    ];
    assert_eq!(
        record_scope(&perl, &rust).reference_only,
        groups(&["21:-4-295"]),
        "a Location the scored engine wrote under another consequence is not reference-only"
    );
}

#[test]
fn a_reference_coordinate_past_the_buffer_bound_is_an_overflow() {
    let perl = [
        intergenic("21:30425651-2147483647", "deletion"),
        intergenic("21:2147483649", "N."),
        intergenic("21:30425651-2147478647", "deletion"),
    ];
    let rust = [
        key(
            "21:2147483649",
            "N[21:30425700[",
            "T1",
            "upstream_gene_variant",
        ),
        key(
            "21:40000001-2147483647",
            "deletion",
            "T2",
            "transcript_ablation",
        ),
    ];
    assert_eq!(
        record_scope(&perl, &rust).mate_overflow,
        groups(&["21:30425651-2147483647", "21:2147483649"]),
        "a Location the scored engine alone wrote and an end one below the bound are not"
    );
}

#[test]
fn a_reference_coordinate_at_the_buffer_bound_is_counted_through_the_digit_guard() {
    assert_eq!(BND_MATE_32_BIT_LIMIT.to_string().len(), 10);
    assert_eq!(BND_MATE_32_BIT_LIMIT_DIGITS, 10);
    let perl = [
        intergenic("21:2147478648", "N."),
        intergenic("21:30425651-2147478648", "deletion"),
        intergenic("2:5-2147478648", "deletion"),
        intergenic("21:2147478647", "N."),
        intergenic("21:30425651-999999999", "deletion"),
        intergenic("21:-2147483647--2147483348", "deletion"),
    ];
    assert_eq!(
        record_scope(&perl, &[]).mate_overflow,
        groups(&["21:2147478648", "21:30425651-2147478648", "2:5-2147478648"]),
        "one below the bound, a coordinate of nine digits and a negative one of ten are not \
         counted"
    );
}
