// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The classes of the 116.2 reference as the table states them: which they are and in which
//! order, which of them the adjusted F1 sets aside, what each definition and citation must
//! say, and that the taxonomy page names every class.

use super::*;
use crate::rules::filter::SV_TRANSCRIPT_SELECTION_BUCKET;
use crate::rules::one_sided::{
    BND_BARE_REFERENCE_ALLELE_ROW_BUCKET, BND_MATE_32_BIT_BUFFER_BUCKET, BND_OWN_ROWS_LOST_BUCKET,
    BND_SYNONYM_MATE_ROWS_LOST_BUCKET, DOT_BREAKEND_PAIR_BUCKET,
    FASTA_NAMED_SLICE_INTERGENIC_BUCKET, GVCF_VARIANT_SITE_DROPPED_BUCKET, NON_INTEGER_END_BUCKET,
    NON_NUMERIC_POS_NUMIFIED_BUCKET, REFERENCE_SKIPPED_RECORD_BUCKET,
};
use crate::rules::{row_pair_rules, EXCLUSION_REGISTRY};
use std::collections::HashSet;

/// The page that defines every class for a reader of a report.
const TAXONOMY_PAGE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../docs/concordance-provenance/116.2-divergence-classes.md"
);

/// The classes that read the two consequence sets of a pair, in classification order.
const SET_PAIR_CLASSES: [&str; 6] = [
    "intronic_overcall_swap",
    "ppt_overcall_swap",
    "frameshift_splice_cascade_swap",
    "splice_lastwrite_swap",
    "start_cooccurrence_swap",
    "covered_splice_region_swap",
];

/// The classes that read the two rows of a pair, in classification order.
const ROW_PAIR_CLASSES: [&str; 16] = [
    "structural_deletion_stop_lost_without_stop_codon",
    "stop_codon_completed_by_polyadenylation",
    "stop_codon_overlap_without_stop_codon",
    "stop_codon_behind_last_residue_without_stop_codon",
    "start_codon_overlap_without_start_codon",
    "stop_codon_overlap_genomic_span",
    "stop_retained_without_stop_codon_test",
    "stop_retained_three_base_stop_codon_test",
    "stop_codon_removed_without_stop_codon_test",
    "stop_codon_overlap_insertion_range",
    "stop_codon_settled_before_unknown_bases",
    "frameshift_inside_a_held_stop_codon",
    "splice_term_past_the_record_span",
    "pachysolen_nuclear_code_read_as_standard",
    "start_codon_removed_without_five_prime_utr",
    "start_codon_kept_without_five_prime_utr",
];

/// The one-sided classes: the skipped record, then the thirteen record shapes.
const ONE_SIDED_CLASSES: [&str; 14] = [
    "reference_skipped_record",
    "bnd_own_chromosome_rows_lost",
    "bnd_synonym_mate_rows_lost",
    "fasta_named_slice_intergenic",
    "bnd_bare_reference_allele_row",
    "bnd_mate_32bit_buffer_intergenic",
    "gvcf_variant_site_dropped",
    "dot_breakend_without_svtype_pair",
    "missing_alt_annotated",
    "non_integer_end_annotated",
    "cnv_tr_empty_expansion_dropped",
    "star_alone_svtype_annotated",
    "non_numeric_pos_numified",
    "cnv_tr_rn_zero_allele_dropped",
];

fn buckets(classes: &[DivergenceClass]) -> Vec<&str> {
    classes.iter().map(|class| class.bucket).collect()
}

fn one_sided_classes() -> &'static [DivergenceClass] {
    &DIVERGENCE_CLASSES_116_2[PAIR_RULE_COUNT..]
}

fn one_sided(bucket: &str) -> &'static DivergenceClass {
    one_sided_classes()
        .iter()
        .find(|class| class.bucket == bucket)
        .unwrap_or_else(|| panic!("no one-sided class {bucket}"))
}

#[test]
fn the_one_sided_classes_are_the_skipped_record_and_the_thirteen_record_shapes() {
    let classes = one_sided_classes();
    assert_eq!(buckets(classes), ONE_SIDED_CLASSES);
    for named_by_the_pass in [
        REFERENCE_SKIPPED_RECORD_BUCKET,
        BND_OWN_ROWS_LOST_BUCKET,
        BND_SYNONYM_MATE_ROWS_LOST_BUCKET,
        FASTA_NAMED_SLICE_INTERGENIC_BUCKET,
        BND_BARE_REFERENCE_ALLELE_ROW_BUCKET,
        BND_MATE_32_BIT_BUFFER_BUCKET,
        GVCF_VARIANT_SITE_DROPPED_BUCKET,
        DOT_BREAKEND_PAIR_BUCKET,
        NON_INTEGER_END_BUCKET,
        NON_NUMERIC_POS_NUMIFIED_BUCKET,
    ] {
        assert!(
            ONE_SIDED_CLASSES.contains(&named_by_the_pass),
            "{named_by_the_pass}"
        );
    }
    let sets_aside: Vec<bool> = classes.iter().map(|class| class.excludes).collect();
    assert_eq!(
        sets_aside,
        [
            true, false, false, false, true, false, true, true, false, true, false, false, true,
            false
        ],
        "the lost-row breakend shapes, the missing ALT, the empty expansion, the `*`-only ALT \
         and the dropped `RN=0` allele need a cache or the input record, which this comparison \
         does not read for them"
    );
    let never_set_aside: Vec<&str> = classes
        .iter()
        .filter(|class| class.definition.contains("never set aside"))
        .map(|class| class.bucket)
        .collect();
    assert_eq!(
        never_set_aside,
        [
            FASTA_NAMED_SLICE_INTERGENIC_BUCKET,
            BND_MATE_32_BIT_BUFFER_BUCKET
        ],
        "the two classes no comparison sets aside say so"
    );
    assert!(one_sided(NON_INTEGER_END_BUCKET)
        .definition
        .contains("SVLEN"));
    assert!(one_sided(BND_MATE_32_BIT_BUFFER_BUCKET)
        .definition
        .contains("span's end"));
    let overflow = &classes[5];
    assert_eq!(overflow.bucket, "bnd_mate_32bit_buffer_intergenic");
    assert!(overflow.definition.contains("2147478648"));
    let bare = &classes[4];
    assert_eq!(bare.bucket, "bnd_bare_reference_allele_row");
    assert!(bare
        .perl_citation
        .contains("StructuralVariationOverlap.pm:64-77"));
    assert!(bare
        .definition
        .contains("reference allele has no variant effect"));
    assert!(bare.definition.contains("set aside"));
    for class in classes {
        assert!(
            class.definition.ends_with('.') && class.definition.split_whitespace().count() >= 12,
            "{}",
            class.bucket
        );
        assert!(class.perl_citation.contains(".pm"), "{}", class.bucket);
    }
    let distinct_buckets: HashSet<&str> = classes.iter().map(|class| class.bucket).collect();
    assert_eq!(distinct_buckets.len(), classes.len());
    let distinct_taxonomy_classes: HashSet<&str> =
        classes.iter().map(|class| class.taxonomy_class).collect();
    assert_eq!(distinct_taxonomy_classes.len(), classes.len());
}

#[test]
fn every_class_states_a_definition_and_cites_its_source_in_report_order() {
    let report_order = [
        &SET_PAIR_CLASSES[..],
        &ROW_PAIR_CLASSES[..],
        &ONE_SIDED_CLASSES[..],
    ]
    .concat();
    assert_eq!(buckets(DIVERGENCE_CLASSES_116_2), report_order);
    assert_eq!(
        PAIR_RULE_COUNT,
        SET_PAIR_CLASSES.len() + ROW_PAIR_CLASSES.len()
    );
    let (pair_classes, _) = DIVERGENCE_CLASSES_116_2.split_at(PAIR_RULE_COUNT);
    let registry: Vec<(&str, bool)> = EXCLUSION_REGISTRY
        .iter()
        .map(|rule| (rule.bucket, rule.excludes))
        .collect();
    let stated: Vec<(&str, bool)> = pair_classes
        .iter()
        .map(|class| (class.bucket, class.excludes))
        .collect();
    assert_eq!(
        stated, registry,
        "a pair class is set aside in the report where its rule sets it aside"
    );
    for class in DIVERGENCE_CLASSES_116_2 {
        assert!(class.definition.ends_with('.'), "{}", class.bucket);
        assert!(class.perl_citation.contains(".pm"), "{}", class.bucket);
    }
}

#[test]
fn the_taxonomy_page_names_every_class_of_both_references() {
    let page =
        std::fs::read_to_string(TAXONOMY_PAGE).unwrap_or_else(|e| panic!("{TAXONOMY_PAGE}: {e}"));
    let every_class = DIVERGENCE_CLASSES_116_2
        .iter()
        .map(|class| class.bucket)
        .chain([SV_TRANSCRIPT_SELECTION_BUCKET]);
    for bucket in every_class {
        assert!(page.contains(&format!("`{bucket}`")), "{bucket}");
    }
    for release in ["115.2", "116.2"] {
        assert!(page.contains(&format!("`{release}`")), "{release}");
    }
    let heading = "The pair classes that read the two rows of a pair";
    let (_, after) = page
        .split_once(&format!(". {heading}\n"))
        .unwrap_or_else(|| panic!("the page has no section {heading:?}"));
    let section = after
        .split_once("\n#")
        .map_or(after, |(section, _)| section);
    let listed: Vec<&str> = section
        .lines()
        .filter_map(|line| line.strip_prefix("- [`"))
        .map(|rest| rest.split_once('`').map_or(rest, |(name, _)| name))
        .collect();
    let row_aware: Vec<&str> = row_pair_rules().iter().map(|rule| rule.bucket).collect();
    assert_eq!(row_aware, ROW_PAIR_CLASSES);
    assert_eq!(
        listed, row_aware,
        "the section lists the row-aware pair rules, in their order, and no other class"
    );
}
