// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The intended-divergence rules of the adjusted F1: the pair rules over consequence-set
//! swaps, the one-sided rules over rows one engine wrote and the other did not, the
//! collectors that feed them from the sorted merge, and the pass that applies them to
//! `discordant.tsv`.

pub mod breakend;
pub mod cache_transcripts;
pub mod coding_end;
/// What stands at the two ends of a cache transcript's coding sequence, and the coding terms
/// vep-rs owes on a row there: the source files of `tools/vep-formats-compare`, compiled
/// here so that the comparison tools read a cache record by one definition. A root that
/// exports nothing reads the names of that tool's public items by the lints of private ones.
#[allow(dead_code, clippy::enum_variant_names)]
#[path = "../../../vep-formats-compare/src/coding_ends.rs"]
pub mod coding_ends;
pub mod filter;
pub mod gvcf;
pub mod one_sided;
pub mod pair;
pub mod reference_bases;
pub mod row;
pub mod scope;
pub mod shared_dbid;
pub mod shared_stable_id;
pub mod shared_synonym;
pub mod table;

use coding_end::AtACodingEnd;
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
    /// A row-aware rule that asks the coding ends of the pair's transcript: the two rows and
    /// what the comparison read of the transcript's cache record.
    RowsAtACodingEnd(fn(&PairRow, &PairRow, &AtACodingEnd<'_>) -> bool),
}

impl ExclusionRule {
    /// A rule over the two consequence sets and the allele on the pair's tuples; a row-aware
    /// rule matches nothing here.
    pub fn matches(&self, perl_csq: &str, rust_csq: &str, allele: &str) -> bool {
        match self.predicate {
            Predicate::Sets(f) => f(perl_csq, rust_csq),
            Predicate::SetsAndAllele(f) => f(perl_csq, rust_csq, allele),
            Predicate::Rows(_) | Predicate::RowsAtACodingEnd(_) => false,
        }
    }

    /// Whether the rule reads the two rows of the pair in the two outputs.
    pub fn row_aware(&self) -> bool {
        matches!(
            self.predicate,
            Predicate::Rows(_) | Predicate::RowsAtACodingEnd(_)
        )
    }

    /// A row-aware rule on the rows of the two outputs that carry the pair's tuples, in a
    /// comparison that was given no cache (`matches_rows_at`).
    #[allow(dead_code)]
    pub fn matches_rows(&self, perl: &[PairRow], rust: &[PairRow]) -> bool {
        self.matches_rows_at(perl, rust, &AtACodingEnd::NO_CACHE)
    }

    /// A row-aware rule on the rows of the two outputs that carry the pair's tuples: every
    /// row of the reference's tuple against every row of vep-rs's, so the rows of one key
    /// stand or fall together, a rule that asks the coding ends of the pair's transcript
    /// with `at`. False for any other rule and for a side without a row.
    pub fn matches_rows_at(
        &self,
        perl: &[PairRow],
        rust: &[PairRow],
        at: &AtACodingEnd<'_>,
    ) -> bool {
        let every = |f: &dyn Fn(&PairRow, &PairRow) -> bool| {
            !perl.is_empty()
                && !rust.is_empty()
                && perl.iter().all(|e| rust.iter().all(|a| f(e, a)))
        };
        match self.predicate {
            Predicate::Rows(f) => every(&f),
            Predicate::RowsAtACodingEnd(f) => every(&|e, a| f(e, a, at)),
            _ => false,
        }
    }
}

/// The pair rules in classification order, most specific first, since classification is
/// first-match; `covered_splice_region_swap` is the exclusion policy of the splice family,
/// a property of the pair that cuts across its four buckets, so it excludes without
/// classifying.
/// The row-aware rules of the 116.2 reference come last, in their own classification
/// order: each differs from its reference's set in a stop codon term or a start codon term, in
/// `frameshift_variant` on a row whose reference residues open with a stop, in a splice term
/// the reference alone holds, in `synonymous_variant` against `missense_variant`, in a
/// `start_lost` vep-rs alone holds, or in a `start_retained_variant` vep-rs alone holds
/// beside no `start_lost`, which no rule before them admits. One pair satisfies two
/// of them: a stop codon the inserted bases keep before unknown bases, in the place of the
/// reference's `frameshift_variant`, is a pair of `stop_codon_overlap_insertion_range` by its
/// terms and of `stop_codon_settled_before_unknown_bases` by its codons and residues, and is
/// counted under the first. The five
/// classes of a coding end stand first among them, the stop end's before the start end's, so
/// a pair that carries a charged term of both ends is counted under the stop end's class; each
/// asks the reading of the transcript's cache record and takes nothing in a comparison that
/// was given no cache. The five stop codon classes after them, and the splice class for the
/// rest of a difference it shares with one of the ten, ask the same reading beside a cache
/// and stand without one as they stood. `stop_codon_settled_before_unknown_bases` and
/// `frameshift_inside_a_held_stop_codon`, between those five and the splice class, and the
/// last three, the genetic code 26 class, `start_codon_removed_without_five_prime_utr` and
/// `start_codon_kept_without_five_prime_utr`, read the two rows alone.
/// A pair whose allele opens with its whole reference allele is put to the
/// row-aware rules also where a set rule classifies it without setting it aside, and is counted
/// under the row-aware rule that takes it.
pub const EXCLUSION_REGISTRY: [ExclusionRule; 22] = [
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
        bucket: "structural_deletion_stop_lost_without_stop_codon",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(
            row::structural_deletion_stop_lost_without_stop_codon_pair,
        ),
    },
    ExclusionRule {
        bucket: "stop_codon_completed_by_polyadenylation",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(row::stop_codon_completed_by_polyadenylation_pair),
    },
    ExclusionRule {
        bucket: "stop_codon_overlap_without_stop_codon",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(row::stop_codon_overlap_without_stop_codon_pair),
    },
    ExclusionRule {
        bucket: "stop_codon_behind_last_residue_without_stop_codon",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(
            row::stop_codon_behind_last_residue_without_stop_codon_pair,
        ),
    },
    ExclusionRule {
        bucket: "start_codon_overlap_without_start_codon",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(row::start_codon_overlap_without_start_codon_pair),
    },
    ExclusionRule {
        bucket: "stop_codon_overlap_genomic_span",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(row::stop_codon_overlap_genomic_span_pair_at),
    },
    ExclusionRule {
        bucket: "stop_retained_without_stop_codon_test",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(row::stop_retained_without_stop_codon_test_pair_at),
    },
    ExclusionRule {
        bucket: "stop_retained_three_base_stop_codon_test",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(
            row::stop_retained_three_base_stop_codon_test_pair_at,
        ),
    },
    ExclusionRule {
        bucket: "stop_codon_removed_without_stop_codon_test",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(
            row::stop_codon_removed_without_stop_codon_test_pair_at,
        ),
    },
    ExclusionRule {
        bucket: "stop_codon_overlap_insertion_range",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(row::stop_codon_overlap_insertion_range_pair_at),
    },
    ExclusionRule {
        bucket: "stop_codon_settled_before_unknown_bases",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::stop_codon_settled_before_unknown_bases_pair),
    },
    ExclusionRule {
        bucket: "frameshift_inside_a_held_stop_codon",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::frameshift_inside_a_held_stop_codon_pair),
    },
    ExclusionRule {
        bucket: "splice_term_past_the_record_span",
        classifies: true,
        excludes: true,
        predicate: Predicate::RowsAtACodingEnd(row::splice_term_past_the_record_span_pair_at),
    },
    ExclusionRule {
        bucket: "pachysolen_nuclear_code_read_as_standard",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::pachysolen_nuclear_code_read_as_standard_pair),
    },
    ExclusionRule {
        bucket: "start_codon_removed_without_five_prime_utr",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::start_codon_removed_without_five_prime_utr_pair),
    },
    ExclusionRule {
        bucket: "start_codon_kept_without_five_prime_utr",
        classifies: true,
        excludes: true,
        predicate: Predicate::Rows(row::start_codon_kept_without_five_prime_utr_pair),
    },
];

/// The row-aware rules of the registry in classification order. They belong to the 116.2
/// reference and are asked, after every other rule, for a pair none of those classifies or
/// excludes, and for the pair of an insertion written after its span that a set rule only
/// counts (`EXCLUSION_REGISTRY`).
pub fn row_pair_rules() -> Vec<&'static ExclusionRule> {
    let registry: &'static [ExclusionRule] = &EXCLUSION_REGISTRY;
    registry.iter().filter(|rule| rule.row_aware()).collect()
}

/// The bucket of the set rule that stands on a start codon: beside a cache it takes a pair
/// only on a transcript whose record shows one (`row::start_codon_stands`).
pub const START_COOCCURRENCE_SWAP: &str = "start_cooccurrence_swap";

/// The bucket of the first classifying rule a pair matches.
#[allow(dead_code)]
pub fn classify(perl_csq: &str, rust_csq: &str, allele: &str) -> Option<&'static str> {
    classify_where(perl_csq, rust_csq, allele, &|_| true)
}

/// Whether any excluding rule matches the pair.
#[allow(dead_code)]
pub fn excluded(perl_csq: &str, rust_csq: &str, allele: &str) -> bool {
    excluded_where(perl_csq, rust_csq, allele, &|_| true)
}

/// The bucket of the first classifying rule that matches a pair and that `stands` admits;
/// `stands` is asked of the rules that match the pair's sets.
pub fn classify_where(
    perl_csq: &str,
    rust_csq: &str,
    allele: &str,
    stands: &dyn Fn(&ExclusionRule) -> bool,
) -> Option<&'static str> {
    EXCLUSION_REGISTRY
        .iter()
        .filter(|rule| rule.classifies)
        .find(|rule| rule.matches(perl_csq, rust_csq, allele) && stands(rule))
        .map(|rule| rule.bucket)
}

/// Whether any excluding rule that `stands` admits matches the pair.
pub fn excluded_where(
    perl_csq: &str,
    rust_csq: &str,
    allele: &str,
    stands: &dyn Fn(&ExclusionRule) -> bool,
) -> bool {
    EXCLUSION_REGISTRY
        .iter()
        .filter(|rule| rule.excludes)
        .any(|rule| rule.matches(perl_csq, rust_csq, allele) && stands(rule))
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

#[cfg(test)]
#[path = "mod_registry_tests.rs"]
mod registry_tests;
