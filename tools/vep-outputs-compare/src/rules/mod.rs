// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The intended-divergence rules of the adjusted F1: the pair rules over consequence-set
//! swaps, the one-sided rules over rows one engine wrote and the other did not, the
//! collectors that feed them from the sorted merge, and the pass that applies them to
//! `discordant.tsv`.

pub mod breakend;
pub mod filter;
pub mod gvcf;
pub mod one_sided;
pub mod pair;
pub mod row;
pub mod scope;
pub mod table;

use row::PairRow;

/// One pair rule: its reporting bucket, whether it classifies a pair into that bucket and
/// whether the adjusted F1 removes the pairs it recognises. The two roles are independent:
/// of the rules that read the two consequence sets, only the start co-emission rule does both.
pub struct ExclusionRule {
    pub bucket: &'static str,
    pub classifies: bool,
    pub excludes: bool,
    predicate: Predicate,
}

enum Predicate {
    Sets(fn(&str, &str) -> bool),
    SetsAndAllele(fn(&str, &str, &str) -> bool),
    /// A row-aware rule: the reference's and vep-rs's row of the pair.
    Rows(fn(&PairRow, &PairRow) -> bool),
}

impl ExclusionRule {
    /// A rule over the two consequence sets and the allele on the pair's tuples; a row-aware
    /// rule matches nothing here.
    pub fn matches(&self, perl_csq: &str, rust_csq: &str, allele: &str) -> bool {
        match self.predicate {
            Predicate::Sets(f) => f(perl_csq, rust_csq),
            Predicate::SetsAndAllele(f) => f(perl_csq, rust_csq, allele),
            Predicate::Rows(_) => false,
        }
    }

    /// Whether the rule reads the two rows of the pair in the two outputs.
    pub fn row_aware(&self) -> bool {
        matches!(self.predicate, Predicate::Rows(_))
    }

    /// A row-aware rule on the rows of the two outputs that carry the pair's tuples: every
    /// row of the reference's tuple against every row of vep-rs's, so the rows of one key
    /// stand or fall together. False for any other rule and for a side without a row.
    pub fn matches_rows(&self, perl: &[PairRow], rust: &[PairRow]) -> bool {
        match self.predicate {
            Predicate::Rows(f) => {
                !perl.is_empty()
                    && !rust.is_empty()
                    && perl.iter().all(|e| rust.iter().all(|a| f(e, a)))
            }
            _ => false,
        }
    }
}

/// The pair rules in classification order, most specific first, since classification is
/// first-match; `covered_splice_region_swap` is the exclusion policy of the splice family,
/// a property of the pair that cuts across its four buckets, so it excludes without
/// classifying.
/// The row-aware rules of the 116.2 reference come last, in their own classification
/// order: each differs from its reference's set in a stop codon term, or in a splice term the
/// reference alone holds, which no rule before them admits, and no pair satisfies two of them.
pub const EXCLUSION_REGISTRY: [ExclusionRule; 12] = [
    ExclusionRule {
        bucket: "intronic_overcall_swap",
        classifies: true,
        excludes: false,
        predicate: Predicate::Sets(pair::is_intronic_overcall_swap),
    },
    ExclusionRule {
        bucket: "ppt_overcall_swap",
        classifies: true,
        excludes: false,
        predicate: Predicate::Sets(pair::is_ppt_overcall_swap),
    },
    ExclusionRule {
        bucket: "frameshift_splice_cascade_swap",
        classifies: true,
        excludes: false,
        predicate: Predicate::Sets(pair::is_frameshift_splice_cascade_swap),
    },
    ExclusionRule {
        bucket: "splice_lastwrite_swap",
        classifies: true,
        excludes: false,
        predicate: Predicate::Sets(pair::is_splice_overcall_swap),
    },
    ExclusionRule {
        bucket: "start_cooccurrence_swap",
        classifies: true,
        excludes: true,
        predicate: Predicate::SetsAndAllele(pair::is_start_cooccurrence_swap),
    },
    ExclusionRule {
        bucket: "covered_splice_region_swap",
        classifies: false,
        excludes: true,
        predicate: Predicate::Sets(pair::is_covered_splice_region_swap),
    },
    ExclusionRule {
        bucket: "stop_codon_overlap_genomic_span",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::stop_codon_overlap_genomic_span_pair),
    },
    ExclusionRule {
        bucket: "stop_retained_without_stop_codon_test",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::stop_retained_without_stop_codon_test_pair),
    },
    ExclusionRule {
        bucket: "stop_retained_three_base_stop_codon_test",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::stop_retained_three_base_stop_codon_test_pair),
    },
    ExclusionRule {
        bucket: "stop_codon_removed_without_stop_codon_test",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::stop_codon_removed_without_stop_codon_test_pair),
    },
    ExclusionRule {
        bucket: "stop_codon_overlap_insertion_range",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::stop_codon_overlap_insertion_range_pair),
    },
    ExclusionRule {
        bucket: "splice_term_past_the_record_span",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::splice_term_past_the_record_span_pair),
    },
];

/// The row-aware rules of the registry in classification order. They belong to the 116.2
/// reference and are asked, after every other rule, for a pair none of those classifies or
/// excludes.
pub fn row_pair_rules() -> Vec<&'static ExclusionRule> {
    let registry: &'static [ExclusionRule] = &EXCLUSION_REGISTRY;
    registry.iter().filter(|rule| rule.row_aware()).collect()
}

/// The bucket of the first classifying rule a pair matches.
pub fn classify(perl_csq: &str, rust_csq: &str, allele: &str) -> Option<&'static str> {
    EXCLUSION_REGISTRY
        .iter()
        .filter(|rule| rule.classifies)
        .find(|rule| rule.matches(perl_csq, rust_csq, allele))
        .map(|rule| rule.bucket)
}

/// Whether any excluding rule matches the pair.
pub fn excluded(perl_csq: &str, rust_csq: &str, allele: &str) -> bool {
    EXCLUSION_REGISTRY
        .iter()
        .filter(|rule| rule.excludes)
        .any(|rule| rule.matches(perl_csq, rust_csq, allele))
}

/// The blanks `int()` drops around a number: those of `str.isspace()` but U+001C to U+001F,
/// which `str.strip()` drops and `int()` refuses.
fn is_int_space(c: char) -> bool {
    crate::pyfmt::is_py_space(c) && !('\u{1c}'..='\u{1f}').contains(&c)
}

/// `int()` of a Python `str` restricted to ASCII: the blanks `int()` drops around the number
/// dropped, an optional sign, decimal digits with single underscores between them; None
/// where Python raises `ValueError` or the value is outside `i128`.
pub fn py_int(text: &str) -> Option<i128> {
    let body = text.trim_matches(is_int_space);
    let (negative, digits) = match body.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, body.strip_prefix('+').unwrap_or(body)),
    };
    if digits.is_empty()
        || digits.starts_with('_')
        || digits.ends_with('_')
        || digits.contains("__")
    {
        return None;
    }
    let mut value: i128 = 0;
    for b in digits.bytes() {
        match b {
            b'_' => {}
            b'0'..=b'9' => {
                value = value.checked_mul(10)?.checked_add(i128::from(b - b'0'))?;
            }
            _ => return None,
        }
    }
    Some(if negative { -value } else { value })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_integer_is_read_between_the_blanks_python_drops_and_no_others() {
        assert_eq!(py_int(" \t1_000\n"), Some(1000));
        assert_eq!(py_int("\u{a0}+12\u{3000}"), Some(12));
        assert_eq!(py_int("\u{2003}-12\u{85}"), Some(-12));
        for refused in ['\u{1c}', '\u{1d}', '\u{1e}', '\u{1f}'] {
            assert_eq!(py_int(&format!("{refused}12")), None, "{refused:?}");
            assert_eq!(py_int(&format!("12{refused}")), None, "{refused:?}");
            assert_eq!(py_int(&format!("\u{a0}{refused}12")), None, "{refused:?}");
        }
    }
}
