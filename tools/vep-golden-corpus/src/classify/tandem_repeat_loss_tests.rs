// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the tandem-repeat-loss module: what the class asks of a key and its two
//! sets, for an allele that states a loss and for one as long as its run, and the records
//! that write a row of both sets on a key several records share.

use super::*;

const CLASS: Option<&str> = Some(CNV_TR_SYMBOLIC_LOSS_NAMED_A_GAIN);

/// Ensembl VEP's set on a row of the class and vep-rs's, as `classify` holds a set: the
/// terms sorted.
const PERL_SET: &str = "coding_sequence_variant,feature_elongation";
const RUST_SET: &str = "coding_sequence_variant,feature_truncation";

/// vep-rs's set on the row of an allele as long as its run: Ensembl VEP's without
/// `feature_elongation`.
const SAME_SET: &str = "coding_sequence_variant";

fn row(name: &str) -> PositionRow {
    PositionRow::of_positions("10-18", "4-12", name)
}

#[test]
fn a_transcript_key_of_a_tandem_repeat_whose_sets_swap_the_two_terms_takes_the_class() {
    let class = |perl: &str, rust: &str| {
        classify_symbolic_loss_pair("Transcript", "tandem_repeat", perl, rust)
    };
    assert_eq!(class(PERL_SET, RUST_SET), CLASS);
    assert_eq!(class("feature_elongation", "feature_truncation"), CLASS);
    assert_eq!(
        class(
            "feature_elongation,intron_variant,splice_region_variant",
            "feature_truncation,intron_variant,splice_region_variant"
        ),
        CLASS
    );
}

#[test]
fn a_transcript_key_of_a_tandem_repeat_whose_vep_rs_set_lacks_the_gain_term_alone_takes_the_class()
{
    let class = |perl: &str, rust: &str| {
        classify_symbolic_loss_pair("Transcript", "tandem_repeat", perl, rust)
    };
    assert_eq!(class(PERL_SET, SAME_SET), CLASS);
    assert_eq!(
        class(
            "coding_sequence_variant,feature_elongation,intron_variant",
            "coding_sequence_variant,intron_variant"
        ),
        CLASS
    );
    assert_eq!(
        class(
            "5_prime_UTR_variant,coding_sequence_variant,feature_elongation,start_lost",
            "5_prime_UTR_variant,coding_sequence_variant,start_lost"
        ),
        CLASS
    );
    assert_eq!(
        class(
            "feature_elongation,non_coding_transcript_exon_variant",
            "non_coding_transcript_exon_variant"
        ),
        CLASS
    );
}

#[test]
fn every_other_pair_beside_an_allele_as_long_as_its_run_takes_no_class() {
    let class = classify_symbolic_loss_pair;
    // vep-rs alone holds the gain term.
    assert_eq!(
        class("Transcript", "tandem_repeat", SAME_SET, PERL_SET),
        None
    );
    // Another Allele, and a feature that is no transcript.
    for allele in ["duplication", "deletion", "-", "TANDEM_REPEAT", ""] {
        assert_eq!(
            class("Transcript", allele, PERL_SET, SAME_SET),
            None,
            "{allele}"
        );
    }
    for feature_type in ["RegulatoryFeature", "MotifFeature", "-"] {
        let named = class(feature_type, "tandem_repeat", PERL_SET, SAME_SET);
        assert_eq!(named, None, "{feature_type}");
    }
    // A set that lacks another term as well, a set with another term in its place, a set
    // without a term, a set that holds `feature_truncation` on both sides, and a set of
    // Ensembl VEP's that holds another term than the gain's.
    for (perl, rust) in [
        (
            "coding_sequence_variant,feature_elongation,intron_variant",
            SAME_SET,
        ),
        (PERL_SET, "3_prime_UTR_variant"),
        (PERL_SET, "coding_sequence_variant,intron_variant"),
        ("feature_elongation", ""),
        (PERL_SET, ""),
        (
            "coding_sequence_variant,feature_elongation,feature_truncation",
            RUST_SET,
        ),
        ("coding_sequence_variant,transcript_amplification", SAME_SET),
        (RUST_SET, SAME_SET),
        (SAME_SET, SAME_SET),
    ] {
        let named = class("Transcript", "tandem_repeat", perl, rust);
        assert_eq!(named, None, "{perl} against {rust}");
    }
}

#[test]
fn every_other_key_and_every_other_pair_of_sets_takes_no_class() {
    let class = classify_symbolic_loss_pair;
    // The two terms the other way round.
    assert_eq!(
        class("Transcript", "tandem_repeat", RUST_SET, PERL_SET),
        None
    );
    // Another Allele, and a feature that is no transcript.
    for allele in ["duplication", "deletion", "-", "TANDEM_REPEAT", ""] {
        assert_eq!(
            class("Transcript", allele, PERL_SET, RUST_SET),
            None,
            "{allele}"
        );
    }
    for feature_type in ["RegulatoryFeature", "MotifFeature", "-"] {
        let named = class(feature_type, "tandem_repeat", PERL_SET, RUST_SET);
        assert_eq!(named, None, "{feature_type}");
    }
    // A second term that differs, a term one side alone holds, and two equal sets.
    for (perl, rust) in [
        (PERL_SET, "3_prime_UTR_variant,feature_truncation"),
        (PERL_SET, "feature_truncation"),
        (
            PERL_SET,
            "coding_sequence_variant,feature_elongation,feature_truncation",
        ),
        ("coding_sequence_variant", RUST_SET),
        ("coding_sequence_variant,transcript_amplification", RUST_SET),
        (PERL_SET, PERL_SET),
        (RUST_SET, RUST_SET),
    ] {
        let named = class("Transcript", "tandem_repeat", perl, rust);
        assert_eq!(named, None, "{perl} against {rust}");
    }
}

#[test]
fn the_records_of_both_are_the_names_a_row_of_each_side_carries() {
    let reference = [row("tr_gain"), row("tr_loss"), row("tr_multi")];
    let candidate = [row("tr_multi"), row("tr_loss"), row("tr_other")];
    let named: Vec<&str> = records_of_both(&reference, &candidate)
        .into_iter()
        .map(String::as_str)
        .collect();
    assert_eq!(named, ["tr_loss", "tr_multi"]);
    assert!(records_of_both(&reference, &[]).is_empty());
    assert!(records_of_both(&[], &candidate).is_empty());
}
