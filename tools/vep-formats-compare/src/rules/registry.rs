// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The registry classes.

use crate::shared::EXCLUDING_RULES;

/// The registry class (`REGISTRY_CLASSES`, restricted to `classes`) whose excluding
/// rule accepts the pair of consequence sets, or `None`: the SNP/indel comparator's
/// `EXCLUDING_RULES` decide, so the pairs set aside here are the pairs its adjusted figure
/// sets aside.
pub fn registry_class(
    expected_consequence: &str,
    actual_consequence: &str,
    allele: &str,
    classes: &[&str],
) -> Option<String> {
    EXCLUDING_RULES
        .iter()
        .find(|rule| {
            classes.contains(&rule.taxonomy_class)
                && rule.matches(expected_consequence, actual_consequence, allele)
        })
        .map(|rule| rule.taxonomy_class.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    const BOTH: &[&str] = &["start_cooccurrence_swap", "splice_family_swap"];
    const BOTH_STARTS: &str = "5_prime_UTR_variant,start_lost,start_retained_variant";
    const ONE_START: &str = "5_prime_UTR_variant,start_retained_variant";

    #[test]
    fn a_pair_is_the_named_class_whose_excluding_rule_accepts_it() {
        let class = |expected: &str, actual: &str, allele: &str, classes: &[&str]| {
            registry_class(expected, actual, allele, classes)
        };
        let swap = Some("start_cooccurrence_swap".to_string());
        assert_eq!(class(BOTH_STARTS, ONE_START, "G", BOTH), swap);
        let lost = "5_prime_UTR_variant,start_lost";
        assert_eq!(class(BOTH_STARTS, lost, "deletion", BOTH), swap);
        assert_eq!(class(BOTH_STARTS, lost, "G", BOTH), None);
        // The rules are asked in their own order, whatever the order of the names.
        assert_eq!(
            class(
                BOTH_STARTS,
                ONE_START,
                "G",
                &["splice_family_swap", "start_cooccurrence_swap"]
            ),
            swap
        );
        let family = Some("splice_family_swap".to_string());
        let shifted = "frameshift_variant";
        for region in [
            "frameshift_variant,splice_region_variant",
            "frameshift_variant,intron_variant,splice_region_variant",
        ] {
            assert_eq!(class(shifted, region, "ACAC", BOTH), family, "{region}");
        }
        assert_eq!(
            class(
                "intron_variant",
                "splice_donor_variant",
                "A",
                &["splice_family_swap"]
            ),
            None
        );
    }

    #[test]
    fn a_class_that_is_not_named_takes_no_pair() {
        let (shifted, region) = (
            "frameshift_variant",
            "frameshift_variant,splice_region_variant",
        );
        assert_eq!(registry_class(shifted, region, "A", &[]), None);
        assert_eq!(
            registry_class(shifted, region, "A", &["start_cooccurrence_swap"]),
            None
        );
        assert_eq!(
            registry_class(BOTH_STARTS, ONE_START, "G", &["splice_family_swap"]),
            None
        );
    }
}
