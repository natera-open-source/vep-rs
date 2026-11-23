// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What this comparator adds to the registry of divergence classes: the verdict of each
//! one-sided class under a cache authority and the bucket names the masks report under.

use crate::failure::Failure;
use crate::rules::table::{DivergenceClass, DIVERGENCE_CLASSES_116_2, PAIR_RULE_COUNT};
use crate::rules::{ExclusionRule, EXCLUSION_REGISTRY};

pub use crate::rules::one_sided::{
    BND_BARE_REFERENCE_ALLELE_ROW_BUCKET, BND_MATE_32_BIT_BUFFER_BUCKET, BND_OWN_ROWS_LOST_BUCKET,
    BND_SYNONYM_MATE_ROWS_LOST_BUCKET, DOT_BREAKEND_PAIR_BUCKET,
    FASTA_NAMED_SLICE_INTERGENIC_BUCKET, GVCF_VARIANT_SITE_DROPPED_BUCKET, NON_INTEGER_END_BUCKET,
    NON_NUMERIC_POS_NUMIFIED_BUCKET, REFERENCE_SKIPPED_RECORD_BUCKET,
};

/// The bucket of a record whose ALT is the missing value and that the reference annotates.
pub const MISSING_ALT_ANNOTATED_BUCKET: &str = "missing_alt_annotated";
/// The bucket of a tandem repeat whose expansion yields no alternate allele.
pub const CNV_TR_EMPTY_EXPANSION_BUCKET: &str = "cnv_tr_empty_expansion_dropped";
/// The bucket of a record whose every ALT is `*` under a true INFO/SVTYPE.
pub const STAR_ALONE_SVTYPE_ANNOTATED_BUCKET: &str = "star_alone_svtype_annotated";
/// The bucket of a two-ALT tandem repeat that loses its `RN=0` allele.
pub const CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET: &str = "cnv_tr_rn_zero_allele_dropped";

/// The one-sided classes of the 116.2 registry whose rows this comparator, which reads the
/// input records and a cache authority, sets aside. The two it only counts are a property
/// of the run (`fasta_named_slice_intergenic`) and of the reference's buffer boundaries
/// (`bnd_mate_32bit_buffer_intergenic`).
const SET_ASIDE_WITH_CACHE: [&str; 12] = [
    REFERENCE_SKIPPED_RECORD_BUCKET,
    BND_OWN_ROWS_LOST_BUCKET,
    BND_SYNONYM_MATE_ROWS_LOST_BUCKET,
    BND_BARE_REFERENCE_ALLELE_ROW_BUCKET,
    GVCF_VARIANT_SITE_DROPPED_BUCKET,
    DOT_BREAKEND_PAIR_BUCKET,
    MISSING_ALT_ANNOTATED_BUCKET,
    NON_INTEGER_END_BUCKET,
    CNV_TR_EMPTY_EXPANSION_BUCKET,
    STAR_ALONE_SVTYPE_ANNOTATED_BUCKET,
    NON_NUMERIC_POS_NUMIFIED_BUCKET,
    CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET,
];

/// Whether this comparator sets the rows of a one-sided class aside: the verdict the report
/// states as the class's `excludes`. False for a bucket the list does not hold.
pub fn excludes_with_cache(bucket: &str) -> bool {
    SET_ASIDE_WITH_CACHE.contains(&bucket)
}

/// The one-sided classes of the 116.2 registry, in registry order.
pub fn one_sided_classes_116_2() -> &'static [DivergenceClass] {
    &DIVERGENCE_CLASSES_116_2[PAIR_RULE_COUNT..]
}

/// The pair rules that read the two consequence sets and that the adjusted F1 sets aside,
/// in registry order.
pub fn excluding_rules() -> Vec<&'static ExclusionRule> {
    let registry: &'static [ExclusionRule] = &EXCLUSION_REGISTRY;
    registry
        .iter()
        .filter(|rule| rule.excludes && !rule.row_aware())
        .collect()
}

/// Refuses a list of one-sided rules one of which has no count among `masked`, the buckets
/// this comparator has a mask for: the report would describe such a rule as applied. The
/// failure names every such bucket between single quotes, in the order of `rules`.
pub fn check_every_rule_has_a_mask(
    rules: &[DivergenceClass],
    masked: &[&str],
) -> Result<(), Failure> {
    let missing: Vec<String> = rules
        .iter()
        .filter(|rule| !masked.contains(&rule.bucket))
        .map(|rule| format!("'{}'", rule.bucket))
        .collect();
    if missing.is_empty() {
        return Ok(());
    }
    Err(Failure::assertion(format!(
        "one-sided rule(s) without a mask in this comparator: [{}]",
        missing.join(", ")
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::failure::FailureKind;

    const ONE_SIDED_BUCKETS: [&str; 14] = [
        REFERENCE_SKIPPED_RECORD_BUCKET,
        BND_OWN_ROWS_LOST_BUCKET,
        BND_SYNONYM_MATE_ROWS_LOST_BUCKET,
        FASTA_NAMED_SLICE_INTERGENIC_BUCKET,
        BND_BARE_REFERENCE_ALLELE_ROW_BUCKET,
        BND_MATE_32_BIT_BUFFER_BUCKET,
        GVCF_VARIANT_SITE_DROPPED_BUCKET,
        DOT_BREAKEND_PAIR_BUCKET,
        MISSING_ALT_ANNOTATED_BUCKET,
        NON_INTEGER_END_BUCKET,
        CNV_TR_EMPTY_EXPANSION_BUCKET,
        STAR_ALONE_SVTYPE_ANNOTATED_BUCKET,
        NON_NUMERIC_POS_NUMIFIED_BUCKET,
        CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET,
    ];

    #[test]
    fn the_named_buckets_are_the_one_sided_classes_of_the_table_in_its_order() {
        let table: Vec<&str> = one_sided_classes_116_2()
            .iter()
            .map(|class| class.bucket)
            .collect();
        assert_eq!(table, ONE_SIDED_BUCKETS);
    }

    #[test]
    fn two_classes_are_counted_and_the_others_set_aside() {
        let verdicts: Vec<bool> = one_sided_classes_116_2()
            .iter()
            .map(|class| excludes_with_cache(class.bucket))
            .collect();
        assert_eq!(
            verdicts,
            [
                true, true, true, false, true, false, true, true, true, true, true, true, true,
                true
            ]
        );
        assert!(!excludes_with_cache("a_sixth_shape"));
    }

    #[test]
    fn the_excluding_rules_are_the_two_over_consequence_sets() {
        let buckets: Vec<&str> = excluding_rules().iter().map(|rule| rule.bucket).collect();
        assert_eq!(
            buckets,
            ["start_cooccurrence_swap", "covered_splice_region_swap"]
        );
    }

    #[test]
    fn a_rule_without_a_mask_is_refused_by_name() {
        assert_eq!(
            check_every_rule_has_a_mask(one_sided_classes_116_2(), &ONE_SIDED_BUCKETS),
            Ok(())
        );
        let unmasked = |bucket: &'static str| DivergenceClass {
            bucket,
            definition: "A shape.",
            excludes: false,
            taxonomy_class: "x",
            perl_citation: "X.pm:1",
        };
        let rules = [
            unmasked(GVCF_VARIANT_SITE_DROPPED_BUCKET),
            unmasked("a_sixth_shape"),
            unmasked("a_seventh_shape"),
        ];
        let failure = check_every_rule_has_a_mask(&rules, &ONE_SIDED_BUCKETS).unwrap_err();
        assert_eq!(failure.kind, FailureKind::Assertion);
        assert_eq!(
            failure.message,
            "one-sided rule(s) without a mask in this comparator: ['a_sixth_shape', 'a_seventh_shape']"
        );
    }
}
