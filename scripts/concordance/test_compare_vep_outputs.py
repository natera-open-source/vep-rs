#!/usr/bin/env python3
"""Unit tests for intended-divergence filter helpers in compare_vep_outputs.

Run with: ``python3 -m pytest -q scripts/concordance/test_compare_vep_outputs.py``
Or: ``python3 scripts/concordance/test_compare_vep_outputs.py``
"""

from __future__ import annotations

import csv
import os
import tempfile
import unittest
from collections import Counter
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

from compare_vep_outputs import (
    CLASSIFICATION_CHAIN,
    EXCLUDING_RULES,
    EXCLUSION_REGISTRY,
    MODE_REQUIRED,
    PER_CLASS_MIN_TUPLES,
    PerClassCounts,
    _is_covered_splice_region_swap,
    _is_frameshift_splice_cascade_swap,
    _is_intronic_overcall_swap,
    _is_ppt_overcall_swap,
    _is_splice_overcall_swap,
    _is_start_cooccurrence_swap,
    is_structural_allele,
    _round_half_up,
    _sort_unique_stream_to_file,
    adjusted_f1_str,
    count_terms_in_discordance,
    count_terms_in_sorted_keys,
    derive_per_class_counts,
    filter_snp_indel_intended_divergences,
    normalize_consequence_set,
    parse_args,
    per_class_floor,
    pool_per_class_counts,
    read_per_class_by_dataset_csv,
    read_per_term_exclusions,
    require_mode_args,
    write_per_class_by_dataset_csv,
    write_per_class_pooled_csv,
)


def _norm(csq: str) -> str:
    """Normalize a consequence set exactly as the comparator does.

    This delegates rather than reimplementing: a local copy of the sort-and-dedupe
    rule can drift from production and let a predicate test pass against ordering
    the comparator would never produce.
    """
    return normalize_consequence_set(csq)


class PptOvercallSwapTests(unittest.TestCase):
    """_is_ppt_overcall_swap covers the four pure-intronic overcall shapes."""

    # --- positive cases ---------------------------------------------------

    def test_intron_ppt_to_intron_ppt_region(self) -> None:
        """Rust adds splice_region to intron + PPT."""
        self.assertTrue(
            _is_ppt_overcall_swap(
                _norm("intron_variant,splice_polypyrimidine_tract_variant"),
                _norm(
                    "intron_variant,splice_polypyrimidine_tract_variant,"
                    "splice_region_variant"
                ),
            )
        )

    def test_intron_region_to_intron_ppt_region(self) -> None:
        """Rust adds PPT to intron + splice_region."""
        self.assertTrue(
            _is_ppt_overcall_swap(
                _norm("intron_variant,splice_region_variant"),
                _norm(
                    "intron_variant,splice_polypyrimidine_tract_variant,"
                    "splice_region_variant"
                ),
            )
        )

    def test_intron_only_to_intron_ppt(self) -> None:
        """Rust adds PPT to bare intron_variant."""
        self.assertTrue(
            _is_ppt_overcall_swap(
                _norm("intron_variant"),
                _norm("intron_variant,splice_polypyrimidine_tract_variant"),
            )
        )

    def test_intron_only_to_intron_region(self) -> None:
        """Rust adds splice_region to bare intron_variant."""
        self.assertTrue(
            _is_ppt_overcall_swap(
                _norm("intron_variant"),
                _norm("intron_variant,splice_region_variant"),
            )
        )

    def test_add_both_ppt_and_region(self) -> None:
        """Rust adds BOTH PPT and splice_region at once: a match."""
        self.assertTrue(
            _is_ppt_overcall_swap(
                _norm("intron_variant"),
                _norm(
                    "intron_variant,splice_polypyrimidine_tract_variant,"
                    "splice_region_variant"
                ),
            )
        )

    # --- negative cases ---------------------------------------------------

    def test_no_intron_on_perl_side_is_rejected(self) -> None:
        """Missing intron_variant on the Perl side is NOT a PPT overcall."""
        self.assertFalse(
            _is_ppt_overcall_swap(
                _norm("splice_region_variant"),
                _norm("intron_variant,splice_region_variant"),
            )
        )

    def test_rust_adds_non_ppt_splice_subterm_is_rejected(self) -> None:
        """Rust adds splice_donor_variant, which is not in {PPT, splice_region}."""
        self.assertFalse(
            _is_ppt_overcall_swap(
                _norm("intron_variant"),
                _norm("intron_variant,splice_donor_variant"),
            )
        )

    def test_perl_has_extra_term_is_rejected(self) -> None:
        """Perl has a term vep-rs lacks: the wrong direction for this filter."""
        self.assertFalse(
            _is_ppt_overcall_swap(
                _norm("intron_variant,splice_region_variant"),
                _norm("intron_variant"),
            )
        )

    def test_identical_sets_rejected(self) -> None:
        """Identical sets are not a swap at all."""
        self.assertFalse(
            _is_ppt_overcall_swap(
                _norm("intron_variant,splice_region_variant"),
                _norm("intron_variant,splice_region_variant"),
            )
        )

    def test_unrelated_coding_swap_rejected(self) -> None:
        """Unrelated missense/synonymous swap is not a PPT overcall."""
        self.assertFalse(
            _is_ppt_overcall_swap(
                _norm("missense_variant"),
                _norm("synonymous_variant"),
            )
        )


class SpliceOvercallSwapTests(unittest.TestCase):
    """_is_splice_overcall_swap catches the frameshift and donor cascade shapes."""

    def test_frameshift_stop_plus_intron_acceptor(self) -> None:
        self.assertTrue(
            _is_splice_overcall_swap(
                _norm("frameshift_variant,stop_gained"),
                _norm(
                    "frameshift_variant,intron_variant,"
                    "splice_acceptor_variant,stop_gained"
                ),
            )
        )

    def test_frameshift_plus_donor_group(self) -> None:
        self.assertTrue(
            _is_splice_overcall_swap(
                _norm("frameshift_variant"),
                _norm(
                    "frameshift_variant,intron_variant,"
                    "splice_donor_5th_base_variant,splice_donor_variant"
                ),
            )
        )

    def test_frameshift_stop_plus_donor_group(self) -> None:
        self.assertTrue(
            _is_splice_overcall_swap(
                _norm("frameshift_variant,stop_gained"),
                _norm(
                    "frameshift_variant,intron_variant,"
                    "splice_donor_5th_base_variant,splice_donor_variant,"
                    "stop_gained"
                ),
            )
        )

    def test_donor_plus_intron_and_5th_base(self) -> None:
        self.assertTrue(
            _is_splice_overcall_swap(
                _norm("splice_donor_variant"),
                _norm(
                    "intron_variant,splice_donor_5th_base_variant,splice_donor_variant"
                ),
            )
        )

    def test_coding_downgrade_not_matched(self) -> None:
        """Coding downgrade: Perl=frameshift, Rust=coding_sequence_variant.
        Not caught: a term downgrade is a real divergence, not an intended
        one.
        """
        self.assertFalse(
            _is_splice_overcall_swap(
                _norm("frameshift_variant"),
                _norm("coding_sequence_variant,intron_variant,splice_acceptor_variant"),
            )
        )

    def test_unrelated_missense_to_synonymous_not_matched(self) -> None:
        self.assertFalse(
            _is_splice_overcall_swap(
                _norm("missense_variant"),
                _norm("synonymous_variant"),
            )
        )


class StartCooccurrenceSwapTests(unittest.TestCase):
    """The start co-emission mask, which decides on the allele kind as well as the sets."""

    def test_sequence_variant_drops_start_lost(self) -> None:
        """On a sequence allele Perl's retention check has certified the ATG intact, so
        the erroneous member is start_lost and vep-rs emits the set without it."""
        self.assertTrue(
            _is_start_cooccurrence_swap(
                _norm("frameshift_variant,start_lost,start_retained_variant"),
                _norm("frameshift_variant,start_retained_variant"),
                "-",
            )
        )
        self.assertTrue(
            _is_start_cooccurrence_swap(
                _norm("5_prime_UTR_variant,start_lost,start_retained_variant"),
                _norm("5_prime_UTR_variant,start_retained_variant"),
                "TGGCGGCGCG",
            )
        )

    def test_sequence_variant_keeping_start_lost_is_not_the_shape(self) -> None:
        """Dropping start_retained_variant on a sequence allele keeps the erroneous member."""
        self.assertFalse(
            _is_start_cooccurrence_swap(
                _norm("start_lost,start_retained_variant"),
                _norm("start_lost"),
                "-",
            )
        )

    def test_structural_allele_drops_start_retained(self) -> None:
        """On a structural allele the retention check is skipped, so the erroneous member
        is start_retained_variant and vep-rs keeps start_lost."""
        self.assertTrue(
            _is_start_cooccurrence_swap(
                _norm("feature_truncation,frameshift_variant,start_lost,start_retained_variant"),
                _norm("feature_truncation,frameshift_variant,start_lost"),
                "deletion",
            )
        )
        self.assertFalse(
            _is_start_cooccurrence_swap(
                _norm("start_lost,start_retained_variant"),
                _norm("start_retained_variant"),
                "deletion",
            )
        )

    def test_start_retained_only_not_matched(self) -> None:
        self.assertFalse(
            _is_start_cooccurrence_swap(
                _norm("start_retained_variant"),
                _norm("start_lost"),
                "-",
            )
        )

    def test_is_structural_allele(self) -> None:
        for a in ("-", "A", "TGGCGGCGCG", "acgt", "N"):
            self.assertFalse(is_structural_allele(a), a)
        for a in ("deletion", "duplication", "tandem_repeat", "insertion", "inversion"):
            self.assertTrue(is_structural_allele(a), a)


class DispatchOrderingTests(unittest.TestCase):
    """Confirm the classification chain's nesting, which its ORDER depends on.

    Each narrower predicate's positives must also be positives of the broader one
    that follows it, so the chain's first-match dispatch routes a pair into its
    narrowest bucket. That nesting is also why an arm may be emptied of its
    exclusion but must NOT be deleted: a deleted arm falls through to the generic
    `splice_lastwrite_swap` arm, which re-absorbs its pairs and mislabels them.
    """

    def test_ppt_is_subset_of_splice_overcall(self) -> None:
        """Every PPT positive must also be a splice-overcall positive."""
        ppt_positives = [
            (
                "intron_variant,splice_polypyrimidine_tract_variant",
                "intron_variant,splice_polypyrimidine_tract_variant,"
                "splice_region_variant",
            ),
            ("intron_variant", "intron_variant,splice_region_variant"),
            (
                "intron_variant",
                "intron_variant,splice_polypyrimidine_tract_variant",
            ),
            (
                "intron_variant,splice_region_variant",
                "intron_variant,splice_polypyrimidine_tract_variant,"
                "splice_region_variant",
            ),
        ]
        for perl, rust in ppt_positives:
            pn, rn = _norm(perl), _norm(rust)
            self.assertTrue(
                _is_ppt_overcall_swap(pn, rn),
                msg=f"PPT should match: {perl} -> {rust}",
            )
            self.assertTrue(
                _is_splice_overcall_swap(pn, rn),
                msg=(
                    f"Splice overcall must also match (ordering must prefer "
                    f"PPT): {perl} -> {rust}"
                ),
            )

    def test_cascade_is_subset_of_splice_overcall(self) -> None:
        """Every cascade positive must also be a splice-overcall positive.

        The cascade arm has no fixture anywhere else in this file, so without this
        the nesting that makes its arm safe to empty rather than delete is
        unpinned: were the cascade predicate to admit a shape the generic arm
        rejects, deleting the arm would silently drop those pairs from the report
        instead of relabelling them.
        """
        cascade_positives = [
            (
                "frameshift_variant,stop_gained",
                "frameshift_variant,intron_variant,splice_acceptor_variant,stop_gained",
            ),
            (
                "frameshift_variant",
                "frameshift_variant,intron_variant,"
                "splice_donor_5th_base_variant,splice_donor_variant",
            ),
            (
                "frameshift_variant,stop_gained",
                "frameshift_variant,intron_variant,splice_donor_5th_base_variant,"
                "splice_donor_variant,stop_gained",
            ),
        ]
        for perl, rust in cascade_positives:
            pn, rn = _norm(perl), _norm(rust)
            self.assertTrue(
                _is_frameshift_splice_cascade_swap(pn, rn),
                msg=f"cascade should match: {perl} -> {rust}",
            )
            self.assertTrue(
                _is_splice_overcall_swap(pn, rn),
                msg=(
                    f"splice overcall must also match, or the cascade arm cannot "
                    f"safely be emptied: {perl} -> {rust}"
                ),
            )

    def test_intronic_is_subset_of_splice_overcall(self) -> None:
        """The intronic arm nests inside the generic arm too, for the same reason."""
        for perl, rust in (
            ("intron_variant", "intron_variant,splice_region_variant"),
            (
                "intron_variant,splice_polypyrimidine_tract_variant",
                "intron_variant,splice_polypyrimidine_tract_variant,"
                "splice_region_variant",
            ),
        ):
            pn, rn = _norm(perl), _norm(rust)
            self.assertTrue(_is_intronic_overcall_swap(pn, rn))
            self.assertTrue(_is_splice_overcall_swap(pn, rn))


class CoveredSpliceRegionSwapTests(unittest.TestCase):
    """`_is_covered_splice_region_swap`: the ONE shape `_intron_effects` produces.

    Only line 215 of `_intron_effects` assigns a return value rather than a guarded
    `= 1`, so the only consequence term Perl's last-write-wins can drop is
    `splice_region_variant`. Every other splice term's flag is set-once. These tests
    are the boundary of that argument: additions of any other splice term are
    NEGATIVES no matter which reporting bucket they land in.
    """

    # --- positive cases ---------------------------------------------------

    def test_bare_intron_gains_splice_region(self) -> None:
        """A bare intron gaining splice_region, inside the `intronic_overcall_swap` bucket, is covered."""
        self.assertTrue(
            _is_covered_splice_region_swap(
                _norm("intron_variant"),
                _norm("intron_variant,splice_region_variant"),
            )
        )

    def test_ppt_set_gains_only_splice_region(self) -> None:
        """Covered: Perl already carries PPT, Rust adds splice_region alone."""
        self.assertTrue(
            _is_covered_splice_region_swap(
                _norm("intron_variant,splice_polypyrimidine_tract_variant"),
                _norm(
                    "intron_variant,splice_polypyrimidine_tract_variant,"
                    "splice_region_variant"
                ),
            )
        )

    def test_coding_backbone_gains_splice_region(self) -> None:
        """The shape is not restricted to intronic sets; a coding pair qualifies."""
        self.assertTrue(
            _is_covered_splice_region_swap(
                _norm("missense_variant"),
                _norm("missense_variant,splice_region_variant"),
            )
        )

    def test_intron_variant_admitted_alongside_splice_region(self) -> None:
        """A skipped region withholds the intronic flag and splice_region together."""
        self.assertTrue(
            _is_covered_splice_region_swap(
                _norm("missense_variant"),
                _norm("intron_variant,missense_variant,splice_region_variant"),
            )
        )

    # --- negative cases ---------------------------------------------------

    def test_ppt_addition_is_not_covered(self) -> None:
        """The PPT flag is set-once (lines 156/160/193/205), so its absence
        from Perl is not last-write-wins."""
        self.assertFalse(
            _is_covered_splice_region_swap(
                _norm("intron_variant"),
                _norm("intron_variant,splice_polypyrimidine_tract_variant"),
            )
        )

    def test_ppt_and_region_addition_is_not_covered(self) -> None:
        """Adding PPT alongside splice_region fails too: PPT is not droppable."""
        self.assertFalse(
            _is_covered_splice_region_swap(
                _norm("intron_variant"),
                _norm(
                    "intron_variant,splice_polypyrimidine_tract_variant,"
                    "splice_region_variant"
                ),
            )
        )

    def test_cascade_addition_is_not_covered(self) -> None:
        """The frameshift cascade adds acceptor/donor terms behind set-once flags
        (lines 177/181/185/197), which line 215 cannot clear."""
        self.assertFalse(
            _is_covered_splice_region_swap(
                _norm("frameshift_variant,stop_gained"),
                _norm(
                    "frameshift_variant,intron_variant,splice_acceptor_variant,"
                    "stop_gained"
                ),
            )
        )

    def test_donor_subterm_addition_is_not_covered(self) -> None:
        self.assertFalse(
            _is_covered_splice_region_swap(
                _norm("intron_variant"),
                _norm("intron_variant,splice_donor_5th_base_variant"),
            )
        )

    def test_intron_variant_alone_is_not_covered(self) -> None:
        """`splice_region_variant` is required; the intronic flag alone is not the
        mechanism's signature."""
        self.assertFalse(
            _is_covered_splice_region_swap(
                _norm("missense_variant"),
                _norm("intron_variant,missense_variant"),
            )
        )

    def test_perl_side_loss_is_not_covered(self) -> None:
        """A term Perl carries and Rust lacks is a vep-rs downgrade, whatever else
        Rust adds."""
        self.assertFalse(
            _is_covered_splice_region_swap(
                _norm("missense_variant,splice_donor_variant"),
                _norm("missense_variant,splice_region_variant"),
            )
        )

    def test_identical_sets_are_not_covered(self) -> None:
        self.assertFalse(
            _is_covered_splice_region_swap(
                _norm("intron_variant,splice_region_variant"),
                _norm("intron_variant,splice_region_variant"),
            )
        )


class ExclusionRegistryTests(unittest.TestCase):
    """The registry is the single source of truth for the exclusion rules.

    Its two structural invariants are what let classification and exclusion be
    independent without the mask becoming unattributable.
    """

    def test_every_excluding_rule_nests_inside_a_classifying_rule(self) -> None:
        """Otherwise a pair could be excluded while landing in no bucket, and the
        per-bucket `excluded` counts would stop summing to `excluded_perl`. The
        comparator raises on that; this pins the property that makes it
        unreachable."""
        corpus = [
            ("intron_variant", "intron_variant,splice_region_variant"),
            (
                "intron_variant,splice_polypyrimidine_tract_variant",
                "intron_variant,splice_polypyrimidine_tract_variant,"
                "splice_region_variant",
            ),
            ("missense_variant", "missense_variant,splice_region_variant"),
            (
                "missense_variant",
                "intron_variant,missense_variant,splice_region_variant",
            ),
            ("start_lost,start_retained_variant", "start_retained_variant"),
            (
                "5_prime_UTR_variant,start_lost,start_retained_variant",
                "5_prime_UTR_variant,start_retained_variant",
            ),
        ]
        for perl, rust in corpus:
            pn, rn = _norm(perl), _norm(rust)
            for rule in EXCLUDING_RULES:
                if not rule.matches(pn, rn, "-"):
                    continue
                self.assertTrue(
                    any(c.matches(pn, rn, "-") for c in CLASSIFICATION_CHAIN),
                    msg=(
                        f"{rule.bucket} matched {perl} -> {rust} but no "
                        f"classification rule does; the pair would be excluded "
                        f"with no bucket to attribute it to"
                    ),
                )

    def test_covered_shape_nests_inside_the_generic_splice_arm(self) -> None:
        """The specific instance of the above that the mask depends on."""
        for perl, rust in (
            ("intron_variant", "intron_variant,splice_region_variant"),
            ("missense_variant", "missense_variant,splice_region_variant"),
        ):
            pn, rn = _norm(perl), _norm(rust)
            self.assertTrue(_is_covered_splice_region_swap(pn, rn))
            self.assertTrue(_is_splice_overcall_swap(pn, rn))

    def test_the_two_excluding_rules_are_disjoint(self) -> None:
        """`covered` requires a non-empty Rust-side addition; `start` requires an
        empty one, so no pair is double-counted by the `any()` in the dispatch."""
        for perl, rust in (
            ("start_lost,start_retained_variant", "start_retained_variant"),
            ("intron_variant", "intron_variant,splice_region_variant"),
        ):
            pn, rn = _norm(perl), _norm(rust)
            matches = [r.bucket for r in EXCLUDING_RULES if r.matches(pn, rn, "-")]
            self.assertLessEqual(len(matches), 1, msg=f"{perl} -> {rust}: {matches}")

    def test_registry_buckets_are_unique(self) -> None:
        buckets = [r.bucket for r in EXCLUSION_REGISTRY]
        self.assertEqual(len(buckets), len(set(buckets)))

    def test_only_the_start_rule_both_classifies_and_excludes(self) -> None:
        both = [r.bucket for r in EXCLUSION_REGISTRY if r.classifies and r.excludes]
        self.assertEqual(both, ["start_cooccurrence_swap"])

    def test_every_registry_rule_declares_a_taxonomy_class(self) -> None:
        """The taxonomy class is what maps a bucket onto a row of
        `manuscript/data/discordance_taxonomy.csv`; an empty one would make that
        mapping vacuous rather than failing."""
        for rule in EXCLUSION_REGISTRY:
            self.assertTrue(rule.taxonomy_class, msg=rule.bucket)


class FrameshiftSpliceCascadeSwapTests(unittest.TestCase):
    """_is_frameshift_splice_cascade_swap covers the coding-preserving cascade shapes."""

    # --- positive cases ---------------------------------------------------

    def test_frameshift_stop_gained_plus_intron_acceptor(self) -> None:
        """frameshift+stop_gained preserved, Rust adds intron + acceptor."""
        self.assertTrue(
            _is_frameshift_splice_cascade_swap(
                _norm("frameshift_variant,stop_gained"),
                _norm(
                    "frameshift_variant,intron_variant,"
                    "splice_acceptor_variant,stop_gained"
                ),
            )
        )

    def test_frameshift_plus_intron_and_donor_subterms(self) -> None:
        """frameshift preserved, Rust adds intron + donor sub-terms."""
        self.assertTrue(
            _is_frameshift_splice_cascade_swap(
                _norm("frameshift_variant"),
                _norm(
                    "frameshift_variant,intron_variant,"
                    "splice_donor_5th_base_variant,splice_donor_variant"
                ),
            )
        )

    def test_frameshift_stop_gained_plus_donor_cascade(self) -> None:
        """frameshift+stop_gained preserved, Rust adds full donor cascade."""
        self.assertTrue(
            _is_frameshift_splice_cascade_swap(
                _norm("frameshift_variant,stop_gained"),
                _norm(
                    "frameshift_variant,intron_variant,"
                    "splice_donor_5th_base_variant,"
                    "splice_donor_variant,stop_gained"
                ),
            )
        )

    def test_inframe_deletion_preserved_with_splice_cascade(self) -> None:
        """The same cascade shape with inframe_deletion instead of frameshift."""
        self.assertTrue(
            _is_frameshift_splice_cascade_swap(
                _norm("inframe_deletion"),
                _norm(
                    "inframe_deletion,intron_variant,"
                    "splice_donor_5th_base_variant,splice_donor_variant"
                ),
            )
        )

    # --- negative cases ---------------------------------------------------

    def test_downgrade_frameshift_to_coding_seq_is_bug_not_divergence(self) -> None:
        """vep-rs downgrades frameshift -> coding_sequence_variant: a real
        divergence, which must NOT be filtered."""
        self.assertFalse(
            _is_frameshift_splice_cascade_swap(
                _norm(
                    "frameshift_variant,intron_variant,"
                    "splice_acceptor_variant,splice_donor_5th_base_variant,"
                    "splice_donor_variant"
                ),
                _norm(
                    "coding_sequence_variant,intron_variant,"
                    "splice_acceptor_variant,splice_donor_5th_base_variant,"
                    "splice_donor_variant"
                ),
            )
        )

    def test_inframe_downgrade_to_coding_seq_is_bug(self) -> None:
        """vep-rs downgrades inframe_deletion -> coding_sequence_variant."""
        self.assertFalse(
            _is_frameshift_splice_cascade_swap(
                _norm("inframe_deletion,intron_variant,splice_donor_variant"),
                _norm("coding_sequence_variant,intron_variant,splice_donor_variant"),
            )
        )

    def test_missense_with_splice_region_not_cascade(self) -> None:
        """Simple missense+splice_region - no frameshift-like backbone, owned
        by _is_splice_overcall_swap. Must NOT match cascade helper."""
        self.assertFalse(
            _is_frameshift_splice_cascade_swap(
                _norm("missense_variant"),
                _norm("missense_variant,splice_region_variant"),
            )
        )

    def test_no_coding_backbone_does_not_match(self) -> None:
        """Pure splice swap with no coding backbone - not a cascade swap."""
        self.assertFalse(
            _is_frameshift_splice_cascade_swap(
                _norm("splice_donor_variant"),
                _norm("intron_variant,splice_donor_variant"),
            )
        )


class IntronicOvercallSwapTests(unittest.TestCase):
    """_is_intronic_overcall_swap covers the four pure-intronic shapes."""

    # --- positive cases ---------------------------------------------------

    def test_intron_ppt_adds_splice_region(self) -> None:
        """Rust adds splice_region to intron+PPT."""
        self.assertTrue(
            _is_intronic_overcall_swap(
                _norm("intron_variant,splice_polypyrimidine_tract_variant"),
                _norm(
                    "intron_variant,splice_polypyrimidine_tract_variant,"
                    "splice_region_variant"
                ),
            )
        )

    def test_intron_region_adds_ppt(self) -> None:
        """Rust adds PPT to intron+splice_region."""
        self.assertTrue(
            _is_intronic_overcall_swap(
                _norm("intron_variant,splice_region_variant"),
                _norm(
                    "intron_variant,splice_polypyrimidine_tract_variant,"
                    "splice_region_variant"
                ),
            )
        )

    def test_intron_only_adds_ppt(self) -> None:
        """Rust adds PPT to bare intron_variant."""
        self.assertTrue(
            _is_intronic_overcall_swap(
                _norm("intron_variant"),
                _norm("intron_variant,splice_polypyrimidine_tract_variant"),
            )
        )

    def test_intron_only_adds_splice_region(self) -> None:
        """Rust adds splice_region to bare intron_variant."""
        self.assertTrue(
            _is_intronic_overcall_swap(
                _norm("intron_variant"),
                _norm("intron_variant,splice_region_variant"),
            )
        )

    # --- negative cases ---------------------------------------------------

    def test_coding_term_disqualifies(self) -> None:
        """Any coding term on either side disqualifies."""
        self.assertFalse(
            _is_intronic_overcall_swap(
                _norm("intron_variant"),
                _norm("intron_variant,frameshift_variant"),
            )
        )

    def test_utr_term_disqualifies(self) -> None:
        """UTR term on either side disqualifies (not purely intronic)."""
        self.assertFalse(
            _is_intronic_overcall_swap(
                _norm("5_prime_UTR_variant,intron_variant"),
                _norm("5_prime_UTR_variant,intron_variant,splice_region_variant"),
            )
        )

    def test_missing_intron_anchor(self) -> None:
        """Without intron_variant on Perl side, pattern doesn't anchor."""
        self.assertFalse(
            _is_intronic_overcall_swap(
                _norm("splice_region_variant"),
                _norm(
                    "intron_variant,splice_region_variant,"
                    "splice_polypyrimidine_tract_variant"
                ),
            )
        )

    def test_rust_subset_of_perl_not_match(self) -> None:
        """Rust must be a strict superset of Perl (not a subset)."""
        self.assertFalse(
            _is_intronic_overcall_swap(
                _norm("intron_variant,splice_region_variant"),
                _norm("intron_variant"),
            )
        )


# No test class covers a stop_lost/3'UTR trade-off mask because the comparator
# has none (see the note in compare_vep_outputs.py). A test importing a
# predicate the comparator lacks fails the MODULE at import, and a test file
# that cannot be imported reports no failures.


class SortCollationTests(unittest.TestCase):
    """`_sort_unique_stream_to_file` must collate by byte, whatever the caller's locale.

    The sort-merge that computes the intersection compares keys by codepoint, so both
    sides must agree on ordering. Under a UTF-8 locale, `sort` orders "1:117329" before
    "10:110755"; under C collation, and in Python's `sorted`, the order is reversed. A
    disagreement makes the merge walk past matching keys, which can only LOSE
    intersections and therefore depresses F1 with no error raised.

    This fails if the implementation uses `env.setdefault("LC_ALL", "C")`, which a
    caller-exported locale overrides.
    """

    # Chosen so byte order and numeric order disagree: "10:" < "1:" bytewise,
    # because "0" (0x30) sorts before ":" (0x3a).
    KEYS = ("1:117329", "10:110755")

    def _sorted_under(self, locale_value: str | None) -> list[str]:
        prior = os.environ.get("LC_ALL")
        if locale_value is None:
            os.environ.pop("LC_ALL", None)
        else:
            os.environ["LC_ALL"] = locale_value
        try:
            with tempfile.TemporaryDirectory() as d:
                out = Path(d) / "keys.sorted"
                _sort_unique_stream_to_file(iter(self.KEYS), out)
                return out.read_text(encoding="utf-8").split()
        finally:
            if prior is None:
                os.environ.pop("LC_ALL", None)
            else:
                os.environ["LC_ALL"] = prior

    def test_byte_order_under_utf8_locale(self) -> None:
        """An exported UTF-8 locale must not change the collation."""
        self.assertEqual(self._sorted_under("en_US.UTF-8"), ["10:110755", "1:117329"])

    def test_byte_order_with_no_locale_exported(self) -> None:
        self.assertEqual(self._sorted_under(None), ["10:110755", "1:117329"])

    def test_matches_python_sorted(self) -> None:
        """The external sort and Python's fallback must agree, since either may run."""
        self.assertEqual(self._sorted_under("en_US.UTF-8"), sorted(self.KEYS))


class SubstitutionRejectionTests(unittest.TestCase):
    """The splice masks must reject SUBSTITUTIONS, not just accept term-loss.

    Both predicates cite Perl's `_intron_effects` last-write-wins as the mechanism.
    That mechanism reaches exactly one hash key, `splice_region`
    (BaseTranscriptVariationAllele.pm:215 is the only plain assignment among the
    subroutine's fourteen writes). `start_splice_site` and `end_splice_site`, behind
    splice_donor_variant and splice_acceptor_variant, are guarded write-1 flags that
    no later differing region can clear. So Perl CANNOT lose those terms this way,
    and a pair where Perl carries one and vep-rs does not is a vep-rs downgrade.

    """

    def test_donor_downgraded_to_region_is_rejected(self) -> None:
        """Perl splice_donor -> vep-rs splice_region is a downgrade, not term-loss."""
        self.assertFalse(
            _is_splice_overcall_swap(
                _norm("missense_variant,splice_donor_variant"),
                _norm("missense_variant,splice_region_variant"),
            )
        )

    def test_acceptor_downgraded_to_region_is_rejected(self) -> None:
        self.assertFalse(
            _is_splice_overcall_swap(
                _norm("intron_variant,splice_acceptor_variant"),
                _norm("intron_variant,splice_region_variant"),
            )
        )

    def test_cascade_substitution_is_rejected(self) -> None:
        """The frameshift-cascade sibling must reject the same shape."""
        self.assertFalse(
            _is_frameshift_splice_cascade_swap(
                _norm("frameshift_variant,splice_acceptor_variant"),
                _norm("frameshift_variant,splice_region_variant"),
            )
        )

    def test_cascade_intron_only_addition_is_rejected(self) -> None:
        """Adding only intron_variant is no splice mechanism at all.

        Both predicates require a splice term; a pair whose only difference is
        `intron_variant` matches neither.
        """
        self.assertFalse(
            _is_frameshift_splice_cascade_swap(
                _norm("frameshift_variant"),
                _norm("frameshift_variant,intron_variant"),
            )
        )

    def test_pure_addition_still_masked(self) -> None:
        """The legitimate shape is unaffected: vep-rs adds, Perl loses nothing."""
        self.assertTrue(
            _is_splice_overcall_swap(
                _norm("frameshift_variant,stop_gained"),
                _norm(
                    "frameshift_variant,intron_variant,"
                    "splice_acceptor_variant,stop_gained"
                ),
            )
        )

    def test_start_cooccurrence_rejects_extra_rust_terms(self) -> None:
        """vep-rs adding arbitrary terms is not the co-emission shape.

        The predicate promises "the same consequence set minus the erroneous start
        term"; `pset - rset` alone would leave the vep-rs side unconstrained.
        """
        self.assertFalse(
            _is_start_cooccurrence_swap(
                _norm("start_lost,start_retained_variant"),
                _norm("start_retained_variant,frameshift_variant"),
                "-",
            )
        )

    def test_start_cooccurrence_positive_still_matched(self) -> None:
        self.assertTrue(
            _is_start_cooccurrence_swap(
                _norm("start_lost,start_retained_variant"),
                _norm("start_retained_variant"),
                "-",
            )
        )


class ExclusionBoundTests(unittest.TestCase):
    """Exclusions per side cannot exceed that side's discordant rows.

    Adjusted F1 subtracts these counts from the raw denominators, so an over-count
    inflates F1 and nothing else in the pipeline would notice. A key omitting
    `feature_type` while the scored tuple includes it, so two rows differing only by
    feature type collapsed onto one key and the nested loop counted their cross
    product: 4 exclusions from 2 real tuples.
    """

    HEADER = (
        "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set"
    )

    def _run(
        self, rows: list[tuple[str, ...]]
    ) -> tuple[int, int, dict[str, dict[str, int]]]:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "discordant.tsv"
            with p.open("w", encoding="utf-8") as f:
                f.write(self.HEADER + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")
            return filter_snp_indel_intended_divergences(p)

    def test_two_feature_types_on_one_key_do_not_cross_multiply(self) -> None:
        loc, allele, feat = "21:100", "A", "ENST1"
        perl_csq = _norm("intron_variant")
        rows = [
            ("f", "missing_in_rust", loc, allele, feat, "Transcript", perl_csq),
            ("f", "missing_in_rust", loc, allele, feat, "RegulatoryFeature", perl_csq),
            (
                "f",
                "extra_in_rust",
                loc,
                allele,
                feat,
                "Transcript",
                _norm("intron_variant,splice_region_variant"),
            ),
            (
                "f",
                "extra_in_rust",
                loc,
                allele,
                feat,
                "RegulatoryFeature",
                _norm("intron_variant,splice_polypyrimidine_tract_variant"),
            ),
        ]
        ep, er, buckets = self._run(rows)
        # Both rows are CLASSIFIED, one into each intronic-family bucket, and only
        # the splice_region one is EXCLUDED. That divergence between two rows a
        # single policy would treat identically IS the policy.
        self.assertEqual(buckets["intronic_overcall_swap"], {"count": 2, "excluded": 1})
        self.assertEqual(
            ep, 1, "only the splice_region row is covered; the PPT row is not"
        )
        self.assertEqual(er, 1)
        # A 3-field key would collapse both feature types onto one key and evaluate
        # 4 pairs, two of which are covered, so the count would read 2 rather than 1.
        self.assertLessEqual(ep, 2, "one exclusion per real Perl row at most")

    def test_exclusions_never_exceed_available_rows(self) -> None:
        """The backstop assertion holds on a normal one-per-key input."""
        rows = [
            (
                "f",
                "missing_in_rust",
                "21:200",
                "A",
                "ENST2",
                "Transcript",
                _norm("intron_variant"),
            ),
            (
                "f",
                "extra_in_rust",
                "21:200",
                "A",
                "ENST2",
                "Transcript",
                _norm("intron_variant,splice_region_variant"),
            ),
        ]
        ep, er, _ = self._run(rows)
        self.assertEqual((ep, er), (1, 1))


class ClassifyWithoutExcludingTests(unittest.TestCase):
    """A bucket may count a pair without the adjusted F1 removing it.

    The classify-versus-exclude policy, expressed end-to-end through the file
    reader rather than at the predicate level: naming a divergence shape and
    crediting vep-rs for it are separate acts. Each case below is a real shape
    from the discordance lists.
    """

    HEADER = (
        "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set"
    )

    def _run(self, perl_csq: str, rust_csq: str):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "discordant.tsv"
            with p.open("w", encoding="utf-8") as f:
                f.write(self.HEADER + "\n")
                for source, csq in (
                    ("missing_in_rust", perl_csq),
                    ("extra_in_rust", rust_csq),
                ):
                    f.write(
                        "\t".join(
                            (
                                "f",
                                source,
                                "21:100",
                                "A",
                                "ENST1",
                                "Transcript",
                                _norm(csq),
                            )
                        )
                        + "\n"
                    )
            return filter_snp_indel_intended_divergences(p)

    def _assert(self, perl_csq, rust_csq, bucket, excluded):
        ep, er, buckets = self._run(perl_csq, rust_csq)
        self.assertEqual(
            buckets[bucket],
            {"count": 1, "excluded": 1 if excluded else 0},
            msg=f"{perl_csq} -> {rust_csq} in {bucket}",
        )
        self.assertEqual((ep, er), (1, 1) if excluded else (0, 0))
        self.assertEqual(
            sum(v["excluded"] for v in buckets.values()),
            ep,
            "per-bucket excluded counts must sum to the reported exclusions",
        )

    def test_ppt_addition_is_counted_and_not_excluded(self) -> None:
        """Pure-intronic, so it classifies into `intronic_overcall_swap`."""
        self._assert(
            "intron_variant",
            "intron_variant,splice_polypyrimidine_tract_variant",
            "intronic_overcall_swap",
            excluded=False,
        )

    def test_ppt_addition_on_a_coding_set_is_counted_and_not_excluded(self) -> None:
        """The same addition on a set carrying a UTR term leaves the pure-intronic
        arm and lands in `ppt_overcall_swap`, also unexcluded."""
        self._assert(
            "5_prime_UTR_variant,intron_variant",
            "5_prime_UTR_variant,intron_variant,splice_polypyrimidine_tract_variant",
            "ppt_overcall_swap",
            excluded=False,
        )

    def test_cascade_addition_is_counted_and_not_excluded(self) -> None:
        self._assert(
            "frameshift_variant,stop_gained",
            "frameshift_variant,intron_variant,splice_acceptor_variant,stop_gained",
            "frameshift_splice_cascade_swap",
            excluded=False,
        )

    def test_splice_subterm_addition_is_counted_and_not_excluded(self) -> None:
        self._assert(
            "missense_variant",
            "missense_variant,splice_donor_5th_base_variant",
            "splice_lastwrite_swap",
            excluded=False,
        )

    def test_covered_splice_region_addition_is_counted_and_excluded(self) -> None:
        """The one shape that excludes, and it does so from inside a bucket
        whose own rule declares `excludes=False`."""
        self._assert(
            "missense_variant",
            "missense_variant,splice_region_variant",
            "splice_lastwrite_swap",
            excluded=True,
        )

    def test_start_cooccurrence_is_counted_and_excluded(self) -> None:
        self._assert(
            "start_lost,start_retained_variant",
            "start_retained_variant",
            "start_cooccurrence_swap",
            excluded=True,
        )


# Per-consequence-class F1


def _round_half_up_consumer_side(x: float, places: int) -> str:
    """The consumer-side ROUND_HALF_UP rule, reimplemented on purpose.

    Every other helper in this file delegates to production rather than copying it.
    This one is the deliberate exception: it is the rounding rule a reader re-deriving
    a published value uses, and the property under test is that the comparator agrees
    with it. Delegating to the comparator's own
    `_round_half_up` would make the test tautological.
    """
    return str(Decimal(str(x)).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))


class RoundHalfUpTests(unittest.TestCase):
    """`_round_half_up` matches the consumer-side rule, including on exact halves."""

    def test_matches_consumer_side_rule_on_exact_half(self) -> None:
        # 0.1234565 at 6 places is the case the two rounding rules disagree on.
        self.assertEqual(_round_half_up(0.1234565, 6), _round_half_up_consumer_side(0.1234565, 6))

    def test_half_rounds_up_not_to_even(self) -> None:
        """Python's own formatting rounds 2.675 DOWN at 2 places; this must not."""
        self.assertEqual(_round_half_up(0.125, 2), "0.13")
        self.assertNotEqual(_round_half_up(0.125, 2), f"{0.125:.2f}")

    def test_pads_to_fixed_width(self) -> None:
        self.assertEqual(_round_half_up(1.0, 6), "1.000000")
        self.assertEqual(_round_half_up(0.0, 6), "0.000000")

    def test_near_one_keeps_six_decimals(self) -> None:
        """The published values sit at 0.99999x, so the width must not collapse."""
        self.assertEqual(_round_half_up(0.9999284, 6), "0.999928")


class PerClassCountsTests(unittest.TestCase):
    """The two derivations, and the guards that stop a wrong pairing publishing."""

    def test_intersection_is_vep_total_less_vep_only(self) -> None:
        c = PerClassCounts(term="missense_variant", vep_tuples=100, vep_only=7, engine_only=3)
        self.assertEqual(c.intersection, 93)

    def test_engine_tuples_is_intersection_plus_engine_only(self) -> None:
        c = PerClassCounts(term="missense_variant", vep_tuples=100, vep_only=7, engine_only=3)
        self.assertEqual(c.engine_tuples, 96)

    def test_precision_and_recall_use_their_own_denominators(self) -> None:
        c = PerClassCounts(term="t", vep_tuples=100, vep_only=7, engine_only=3)
        self.assertAlmostEqual(c.precision, 93 / 96)
        self.assertAlmostEqual(c.recall, 93 / 100)

    def test_f1_is_two_i_over_v_plus_e_exactly(self) -> None:
        """The direct form, bit-for-bit, because a consumer re-derives it that way.

        `2PR/(P+R)` is algebraically the same and can differ in the last float bits,
        which would put the published table and the CSV on opposite sides of a
        sixth-decimal boundary although every number in both is correct.
        """
        c = PerClassCounts(term="t", vep_tuples=100, vep_only=7, engine_only=3)
        self.assertEqual(c.f1(), 2 * 93 / (100 + 96))

    def test_perfect_concordance_is_one(self) -> None:
        c = PerClassCounts(term="t", vep_tuples=50, vep_only=0, engine_only=0)
        self.assertEqual(c.f1_str(), "1.000000")

    def test_term_vep_never_emits_is_zero_not_an_error(self) -> None:
        """An engine emitting a term VEP never emits is a finding, not a crash.

        Its recall denominator is zero, so recall is reported as 0.0 rather than
        raising; the tuple threshold is what keeps it out of the reported floor.
        """
        c = PerClassCounts(term="invented_term", vep_tuples=0, vep_only=0, engine_only=4)
        self.assertEqual(c.intersection, 0)
        self.assertEqual(c.engine_tuples, 4)
        self.assertEqual(c.recall, 0.0)
        self.assertEqual(c.f1_str(), "0.000000")

    def test_both_sides_empty_gives_zero_rather_than_dividing_by_zero(self) -> None:
        c = PerClassCounts(term="t", vep_tuples=0, vep_only=0, engine_only=0)
        self.assertEqual(c.f1(), 0.0)
        self.assertEqual(c.precision, 0.0)

    def test_vep_only_exceeding_vep_total_is_rejected(self) -> None:
        """The mismatched-archive guard: the report belongs to a different run.

        Without it the intersection goes NEGATIVE and F1 exceeds 1, or is produced
        from two unrelated runs, and nothing in the emitted CSV shows it.
        """
        with self.assertRaises(AssertionError) as ctx:
            PerClassCounts(term="synonymous_variant", vep_tuples=1, vep_only=5, engine_only=0)
        self.assertIn("synonymous_variant", str(ctx.exception))
        self.assertIn("5 VEP-only", str(ctx.exception))

    def test_negative_count_is_rejected(self) -> None:
        with self.assertRaises(AssertionError):
            PerClassCounts(term="t", vep_tuples=10, vep_only=-1, engine_only=0)
        with self.assertRaises(AssertionError):
            PerClassCounts(term="t", vep_tuples=10, vep_only=0, engine_only=-2)

    def test_vep_only_equal_to_total_is_allowed(self) -> None:
        """Every tuple of a term discordant is legal, and gives F1 0."""
        c = PerClassCounts(term="t", vep_tuples=9, vep_only=9, engine_only=0)
        self.assertEqual(c.intersection, 0)
        self.assertEqual(c.f1_str(), "0.000000")


class CountTermsInSortedKeysTests(unittest.TestCase):
    """Per-term totals come from the deduplicated key file, not the raw output."""

    def _write(self, tmp: Path, lines: list[str]) -> Path:
        path = tmp / "keys.sorted.txt"
        path.write_text("".join(f"{ln}\n" for ln in lines), encoding="utf-8")
        return path

    def test_multi_term_tuple_counts_toward_every_term(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = self._write(
                Path(td),
                ["21:100\tA\tENST1\tTranscript\tintron_variant,splice_region_variant"],
            )
            self.assertEqual(
                count_terms_in_sorted_keys(path),
                Counter({"intron_variant": 1, "splice_region_variant": 1}),
            )

    def test_repeated_term_across_tuples_accumulates(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = self._write(
                Path(td),
                [
                    "21:100\tA\tENST1\tTranscript\tmissense_variant",
                    "21:200\tC\tENST1\tTranscript\tmissense_variant",
                ],
            )
            self.assertEqual(count_terms_in_sorted_keys(path)["missense_variant"], 2)

    def test_blank_lines_are_skipped(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "keys.txt"
            path.write_text(
                "21:100\tA\tENST1\tTranscript\tmissense_variant\n\n", encoding="utf-8"
            )
            self.assertEqual(count_terms_in_sorted_keys(path)["missense_variant"], 1)

    def test_wrong_field_count_raises_rather_than_misparsing(self) -> None:
        """A four-field line would otherwise put a FEATURE name in the term column.

        Splitting on the last tab and trusting it is how a schema change becomes a
        term inventory full of transcript ids that still looks like a valid CSV.
        """
        with tempfile.TemporaryDirectory() as td:
            path = self._write(Path(td), ["21:100\tA\tENST1\tmissense_variant"])
            with self.assertRaises(ValueError):
                count_terms_in_sorted_keys(path)

    def test_empty_file_gives_empty_counter(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "empty.txt"
            path.write_text("", encoding="utf-8")
            self.assertEqual(count_terms_in_sorted_keys(path), Counter())

    def test_counts_the_deduplicated_stream_the_comparator_scores(self) -> None:
        """End-to-end against `_sort_unique_stream_to_file`, which is the contract.

        A duplicated tuple must contribute once, because the aggregate F1 this
        measurement decomposes is scored on the sort-unique'd key set.
        """
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            keys = [
                "21:100\tA\tENST1\tTranscript\tmissense_variant",
                "21:100\tA\tENST1\tTranscript\tmissense_variant",
                "21:200\tC\tENST1\tTranscript\tmissense_variant",
            ]
            out = tmp / "sorted.txt"
            _sort_unique_stream_to_file(iter(keys), out)
            self.assertEqual(count_terms_in_sorted_keys(out)["missense_variant"], 2)


class CountTermsInDiscordanceTests(unittest.TestCase):
    """One-sided per-term counts, keyed on the report's `source` column."""

    COLUMNS = (
        "file_name",
        "source",
        "location",
        "allele",
        "feature",
        "feature_type",
        "consequence_set",
    )

    def _write(self, tmp: Path, rows: list[tuple[str, str]], columns=None) -> Path:
        path = tmp / "discordant.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(columns or self.COLUMNS)
            for i, (source, csq) in enumerate(rows):
                writer.writerow(["f.txt", source, f"21:{i}", "A", "ENST1", "Transcript", csq])
        return path

    def test_splits_the_two_sides(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = self._write(
                Path(td),
                [
                    ("missing_in_rust", "intron_variant"),
                    ("extra_in_rust", "intron_variant,splice_region_variant"),
                ],
            )
            vep_only, engine_only = count_terms_in_discordance(path)
            self.assertEqual(vep_only, Counter({"intron_variant": 1}))
            self.assertEqual(
                engine_only,
                Counter({"intron_variant": 1, "splice_region_variant": 1}),
            )

    def test_unknown_source_is_rejected(self) -> None:
        """A dropped side reads as perfect concordance on that side."""
        with tempfile.TemporaryDirectory() as td:
            path = self._write(Path(td), [("both_ways", "intron_variant")])
            with self.assertRaises(ValueError) as ctx:
                count_terms_in_discordance(path)
            self.assertIn("both_ways", str(ctx.exception))

    def test_missing_consequence_column_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "d.tsv"
            path.write_text("file_name\tsource\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                count_terms_in_discordance(path)
            self.assertIn("consequence_set", str(ctx.exception))

    def test_term_order_in_the_report_does_not_matter(self) -> None:
        """The report is normalised on read, so an unsorted set counts identically."""
        with tempfile.TemporaryDirectory() as td:
            path = self._write(
                Path(td), [("extra_in_rust", "splice_region_variant,intron_variant")]
            )
            _, engine_only = count_terms_in_discordance(path)
            self.assertEqual(engine_only["intron_variant"], 1)
            self.assertEqual(engine_only["splice_region_variant"], 1)

    def test_header_only_report_gives_two_empty_counters(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = self._write(Path(td), [])
            self.assertEqual(count_terms_in_discordance(path), (Counter(), Counter()))


class DerivePerClassCountsTests(unittest.TestCase):
    """Assembly over the union of terms, ordered as the published table is."""

    def test_union_of_all_three_inputs(self) -> None:
        rows = derive_per_class_counts(
            Counter({"a": 10, "b": 5}),
            Counter({"a": 1}),
            Counter({"c": 2}),
        )
        self.assertEqual({r.term for r in rows}, {"a", "b", "c"})

    def test_engine_only_term_is_kept_with_zero_vep_tuples(self) -> None:
        rows = {r.term: r for r in derive_per_class_counts(
            Counter({"a": 10}), Counter(), Counter({"c": 2})
        )}
        self.assertEqual(rows["c"].vep_tuples, 0)
        self.assertEqual(rows["c"].engine_only, 2)

    def test_ordered_by_descending_vep_tuples_then_term(self) -> None:
        rows = derive_per_class_counts(
            Counter({"low": 1, "high": 99, "tie_b": 5, "tie_a": 5}), Counter(), Counter()
        )
        self.assertEqual([r.term for r in rows], ["high", "tie_a", "tie_b", "low"])

    def test_all_empty_gives_no_rows(self) -> None:
        self.assertEqual(derive_per_class_counts(Counter(), Counter(), Counter()), [])


class PoolPerClassCountsTests(unittest.TestCase):
    """Pooling sums the three independent counts, then re-derives."""

    def test_counts_sum_across_datasets(self) -> None:
        pooled = pool_per_class_counts(
            [
                PerClassCounts(term="a", vep_tuples=100, vep_only=2, engine_only=1),
                PerClassCounts(term="a", vep_tuples=50, vep_only=1, engine_only=4),
            ]
        )
        self.assertEqual(len(pooled), 1)
        self.assertEqual(pooled[0].vep_tuples, 150)
        self.assertEqual(pooled[0].vep_only, 3)
        self.assertEqual(pooled[0].engine_only, 5)

    def test_pooling_the_independents_equals_pooling_the_derived(self) -> None:
        """The property that makes storing only three counts safe.

        Both derivations are linear, so summing the independents and re-deriving must
        equal summing the derived values directly. If it ever did not, the CSV would
        carry two disagreeing truths for the same term.
        """
        parts = [
            PerClassCounts(term="a", vep_tuples=100, vep_only=2, engine_only=1),
            PerClassCounts(term="a", vep_tuples=50, vep_only=1, engine_only=4),
            PerClassCounts(term="a", vep_tuples=7, vep_only=0, engine_only=0),
        ]
        pooled = pool_per_class_counts(parts)[0]
        self.assertEqual(pooled.intersection, sum(p.intersection for p in parts))
        self.assertEqual(pooled.engine_tuples, sum(p.engine_tuples for p in parts))

    def test_pooled_f1_is_not_the_mean_of_the_parts(self) -> None:
        """A micro-average, deliberately: a big dataset must dominate a tiny one."""
        pooled = pool_per_class_counts(
            [
                PerClassCounts(term="a", vep_tuples=1_000_000, vep_only=0, engine_only=0),
                PerClassCounts(term="a", vep_tuples=2, vep_only=2, engine_only=0),
            ]
        )[0]
        self.assertEqual(pooled.f1_str(), "0.999999")

    def test_pooled_rows_are_ordered_by_descending_vep_tuples(self) -> None:
        pooled = pool_per_class_counts(
            [
                PerClassCounts(term="small", vep_tuples=1, vep_only=0, engine_only=0),
                PerClassCounts(term="big", vep_tuples=9, vep_only=0, engine_only=0),
            ]
        )
        self.assertEqual([r.term for r in pooled], ["big", "small"])


class PerClassFloorTests(unittest.TestCase):
    """The floor, the threshold that qualifies it, and the count below it."""

    def test_floor_is_the_minimum_among_eligible_terms(self) -> None:
        rows = [
            PerClassCounts(term="good", vep_tuples=10_000, vep_only=1, engine_only=0),
            PerClassCounts(term="worst", vep_tuples=10_000, vep_only=100, engine_only=0),
        ]
        floor = per_class_floor(rows, min_tuples=1000)
        self.assertEqual(floor.term, "worst")
        self.assertEqual(floor.n_eligible, 2)
        self.assertEqual(floor.n_below_threshold, 0)

    def test_a_term_below_the_threshold_cannot_become_the_floor(self) -> None:
        """The whole point of the threshold: 1 discordant tuple in 3 is noise.

        Without it the reported floor would be 0.666667 on a term appearing three
        times, which says nothing about whether a class is systematically mis-called.
        """
        rows = [
            PerClassCounts(term="solid", vep_tuples=10_000, vep_only=1, engine_only=0),
            PerClassCounts(term="rare", vep_tuples=3, vep_only=1, engine_only=0),
        ]
        floor = per_class_floor(rows, min_tuples=1000)
        self.assertEqual(floor.term, "solid")
        self.assertEqual(floor.n_below_threshold, 1)

    def test_no_eligible_term_yields_no_floor_rather_than_a_wrong_one(self) -> None:
        rows = [PerClassCounts(term="rare", vep_tuples=3, vep_only=1, engine_only=0)]
        floor = per_class_floor(rows, min_tuples=1000)
        self.assertIsNone(floor.f1)
        self.assertIsNone(floor.term)
        self.assertEqual(floor.n_eligible, 0)
        self.assertEqual(floor.n_below_threshold, 1)

    def test_ties_break_on_the_term_name_for_reproducibility(self) -> None:
        """Several terms at the same six-decimal value is the expected case."""
        rows = [
            PerClassCounts(term="zeta", vep_tuples=10_000, vep_only=5, engine_only=0),
            PerClassCounts(term="alpha", vep_tuples=10_000, vep_only=5, engine_only=0),
        ]
        self.assertEqual(per_class_floor(rows, min_tuples=1000).term, "alpha")

    def test_term_at_exactly_the_threshold_is_eligible(self) -> None:
        rows = [PerClassCounts(term="edge", vep_tuples=1000, vep_only=1, engine_only=0)]
        self.assertEqual(per_class_floor(rows, min_tuples=1000).n_eligible, 1)

    def test_threshold_default_is_the_pre_committed_value(self) -> None:
        """Pinned so an edit that loosens it has to change a test, not a default.

        Raising the threshold after seeing the data selects the value that makes the
        floor look best.
        """
        self.assertEqual(PER_CLASS_MIN_TUPLES, 1000)
        rows = [PerClassCounts(term="t", vep_tuples=999, vep_only=1, engine_only=0)]
        self.assertEqual(per_class_floor(rows).n_eligible, 0)

    def test_empty_input_yields_no_floor(self) -> None:
        floor = per_class_floor([])
        self.assertIsNone(floor.f1)
        self.assertEqual(floor.n_below_threshold, 0)


class PerClassCsvTests(unittest.TestCase):
    """The two CSVs: their schemas, the round trip, and the append guard."""

    ROWS = [
        PerClassCounts(term="missense_variant", vep_tuples=100, vep_only=3, engine_only=1),
        PerClassCounts(term="intron_variant", vep_tuples=40, vep_only=0, engine_only=2),
    ]
    PINS = {"vep-rs": ("sweep-20260101T000000Z", "deadbeef"), "fastvep": ("sweep-20260101T000000Z", "v0.0.0")}

    def test_by_dataset_round_trip_preserves_the_three_counts(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "by_dataset.csv"
            write_per_class_by_dataset_csv(
                path, self.ROWS, dataset="s01", assembly="GRCh37",
                engine="vep-rs", append=False,
                sweep_id="sweep-20260101T000000Z", binary_md5="deadbeef",
            )
            back, pins = read_per_class_by_dataset_csv(path)
            self.assertEqual(set(back), {"vep-rs"})
            self.assertEqual(pins, {"vep-rs": ("sweep-20260101T000000Z", "deadbeef")})
            self.assertEqual(sorted(back["vep-rs"], key=lambda r: r.term),
                             sorted(self.ROWS, key=lambda r: r.term))

    def test_append_adds_rows_without_a_second_header(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "by_dataset.csv"
            write_per_class_by_dataset_csv(
                path, self.ROWS, dataset="s01", assembly="GRCh37",
                engine="vep-rs", append=False,
                sweep_id="sweep-20260101T000000Z", binary_md5="deadbeef",
            )
            write_per_class_by_dataset_csv(
                path, self.ROWS, dataset="s02", assembly="GRCh38",
                engine="fastvep", append=True,
                sweep_id="sweep-20260101T000000Z", binary_md5="v0.0.0",
            )
            text = path.read_text(encoding="utf-8")
            self.assertEqual(text.count("dataset,assembly,engine"), 1)
            back, pins = read_per_class_by_dataset_csv(path)
            self.assertEqual(set(back), {"vep-rs", "fastvep"})
            self.assertEqual(pins["fastvep"], ("sweep-20260101T000000Z", "v0.0.0"))

    def test_append_into_a_foreign_header_is_rejected(self) -> None:
        """An append under a mismatched header puts values in the wrong columns.

        The result is a CSV that parses, carries plausible integers and means
        something else entirely.
        """
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "by_dataset.csv"
            path.write_text("term,f1\nmissense_variant,1.0\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                write_per_class_by_dataset_csv(
                    path, self.ROWS, dataset="s01", assembly="GRCh37",
                    engine="vep-rs", append=True,
                sweep_id="sweep-20260101T000000Z", binary_md5="deadbeef",
            )
            self.assertIn("cannot append", str(ctx.exception))

    def test_append_to_a_missing_file_creates_it_with_a_header(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "nested" / "by_dataset.csv"
            write_per_class_by_dataset_csv(
                path, self.ROWS, dataset="s01", assembly="GRCh37",
                engine="vep-rs", append=True,
                sweep_id="sweep-20260101T000000Z", binary_md5="deadbeef",
            )
            self.assertTrue(path.read_text(encoding="utf-8").startswith("dataset,"))

    def test_pooled_csv_carries_every_column_a_consumer_reads_by_name(self) -> None:
        """A contract test on the pooled CSV's column names.

        A reader of the pooled CSV selects these eight columns by name, so renaming one
        here breaks that reader rather than only this file.
        """
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "pooled.csv"
            write_per_class_pooled_csv(path, {"vep-rs": list(self.ROWS)}, self.PINS)
            header = next(csv.reader(path.open(encoding="utf-8")))
            for column in (
                "engine", "term", "vep_tuples", "engine_tuples", "intersection", "f1",
                "sweep_id", "binary_md5",
            ):
                self.assertIn(column, header)

    def test_pooled_csv_f1_matches_a_consumer_re_derivation(self) -> None:
        """The emitted `f1` string must equal ROUND_HALF_UP(2I/(V+E), 6) exactly.

        Tested here so a rounding mismatch fails on the comparator's own suite
        rather than only once the CSV reaches a consumer.
        """
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "pooled.csv"
            write_per_class_pooled_csv(path, {"vep-rs": list(self.ROWS)}, self.PINS)
            for row in csv.DictReader(path.open(encoding="utf-8")):
                v, e, i = (
                    int(row["vep_tuples"]),
                    int(row["engine_tuples"]),
                    int(row["intersection"]),
                )
                self.assertEqual(row["f1"], _round_half_up_consumer_side(2 * i / (v + e), 6))

    def test_pooled_csv_puts_vep_rs_first(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "pooled.csv"
            write_per_class_pooled_csv(
                path, {"fastvep": list(self.ROWS), "vep-rs": list(self.ROWS)}, self.PINS
            )
            engines = [r["engine"] for r in csv.DictReader(path.open(encoding="utf-8"))]
            self.assertEqual(engines[0], "vep-rs")
            self.assertEqual(engines[-1], "fastvep")

    def test_read_rejects_a_csv_missing_a_count_column(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "bad.csv"
            path.write_text("engine,term,vep_tuples,sweep_id,binary_md5\nvep-rs,t,1,s,b\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                read_per_class_by_dataset_csv(path)
            self.assertIn("engine_only", str(ctx.exception))

    def test_read_propagates_the_mismatched_archive_guard(self) -> None:
        """A hand-edited CSV cannot smuggle past the guard the counter enforces."""
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "bad.csv"
            path.write_text(
                "dataset,assembly,engine,sweep_id,binary_md5,term,vep_tuples,vep_only,engine_only\n"
                "s01,GRCh37,vep-rs,s,b,t,1,5,0\n",
                encoding="utf-8",
            )
            with self.assertRaises(AssertionError):
                read_per_class_by_dataset_csv(path)

    def test_read_rejects_one_engine_pinned_to_two_generations(self) -> None:
        """Two binaries of one engine pooled together would publish a floor of neither."""
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "mixed.csv"
            path.write_text(
                "dataset,assembly,engine,sweep_id,binary_md5,term,vep_tuples,vep_only,engine_only\n"
                "s01,GRCh37,vep-rs,sweep-A,aaaa,t,10,1,0\n"
                "s02,GRCh38,vep-rs,sweep-B,bbbb,t,10,1,0\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError) as ctx:
                read_per_class_by_dataset_csv(path)
            self.assertIn("more than one (sweep_id, binary_md5)", str(ctx.exception))

    def test_pooled_rows_carry_their_engines_pin(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "pooled.csv"
            write_per_class_pooled_csv(
                path, {"fastvep": list(self.ROWS), "vep-rs": list(self.ROWS)}, self.PINS
            )
            rows = list(csv.DictReader(path.open(encoding="utf-8")))
            self.assertEqual({(r["engine"], r["sweep_id"], r["binary_md5"]) for r in rows},
                             {("vep-rs", "sweep-20260101T000000Z", "deadbeef"), ("fastvep", "sweep-20260101T000000Z", "v0.0.0")})


class PerClassAdjustedF1Tests(unittest.TestCase):
    """The pooled CSV's adj_f1 column: the per-term adjusted-F1 rule, filled for one engine only."""

    ROWS = PerClassCsvTests.ROWS
    PINS = PerClassCsvTests.PINS
    PER_TERM = (
        "term\tvep_only_raw\tengine_only_raw\tvep_only_excluded\tengine_only_excluded"
        "\tvep_only_open\tengine_only_open\n"
        "missense_variant\t3\t1\t2\t1\t1\t0\n"
    )

    def test_adjusted_f1_removes_the_excluded_pairs_from_both_sides(self) -> None:
        """adj = 2I / (2I + (vep_only - ex_v) + (engine_only - ex_e)), the per-term adjusted-F1 rule."""
        row = self.ROWS[0]  # I = 97, vep_only 3, engine_only 1
        self.assertEqual(adjusted_f1_str(row, (0, 0)), row.f1_str())
        self.assertEqual(adjusted_f1_str(row, (2, 1)), _round_half_up_consumer_side(2 * 97 / (2 * 97 + 1 + 0), 6))
        self.assertEqual(adjusted_f1_str(row, (3, 1)), "1.000000")

    def test_an_exclusion_larger_than_the_one_sided_count_is_refused(self) -> None:
        """A per_term.tsv from another run cannot be paired with these counts."""
        with self.assertRaises(AssertionError):
            adjusted_f1_str(self.ROWS[0], (4, 0))

    def test_read_per_term_exclusions_needs_the_two_excluded_columns(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            good = Path(td) / "per_term.tsv"
            good.write_text(self.PER_TERM, encoding="utf-8")
            self.assertEqual(read_per_term_exclusions(good), {"missense_variant": (2, 1)})
            bad = Path(td) / "bad.tsv"
            bad.write_text("term\tvep_only_raw\nmissense_variant\t3\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                read_per_term_exclusions(bad)
            self.assertIn("engine_only_excluded", str(ctx.exception))

    def test_pooled_csv_fills_adj_f1_for_the_named_engine_only(self) -> None:
        """vep-rs rows carry adj_f1; a term absent from per_term.tsv reads its raw F1; fastVEP is empty."""
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "pooled.csv"
            write_per_class_pooled_csv(
                path,
                {"vep-rs": list(self.ROWS), "fastvep": list(self.ROWS)},
                self.PINS,
                {"vep-rs": {"missense_variant": (2, 1)}},
            )
            rows = list(csv.DictReader(path.open(encoding="utf-8")))
            self.assertIn("adj_f1", rows[0])
            by = {(r["engine"], r["term"]): r for r in rows}
            self.assertEqual(by[("vep-rs", "missense_variant")]["adj_f1"], adjusted_f1_str(self.ROWS[0], (2, 1)))
            self.assertEqual(
                (by[("vep-rs", "missense_variant")]["vep_only_excluded"],
                 by[("vep-rs", "missense_variant")]["engine_only_excluded"]),
                ("2", "1"),
            )
            self.assertEqual(by[("vep-rs", "intron_variant")]["adj_f1"], self.ROWS[1].f1_str())
            self.assertEqual(by[("vep-rs", "intron_variant")]["vep_only_excluded"], "0")
            for col in ("adj_f1", "vep_only_excluded", "engine_only_excluded"):
                self.assertEqual(by[("fastvep", "missense_variant")][col], "")
                self.assertEqual(by[("fastvep", "intron_variant")][col], "")
            # The adjusted value re-derives from the row's own integers, like f1 does.
            r = by[("vep-rs", "missense_variant")]
            inter = int(r["intersection"])
            denom = 2 * inter + (int(r["vep_only"]) - int(r["vep_only_excluded"])) + (
                int(r["engine_only"]) - int(r["engine_only_excluded"])
            )
            self.assertEqual(r["adj_f1"], _round_half_up_consumer_side(2 * inter / denom, 6))

    def test_pooled_csv_without_a_per_term_file_leaves_adj_f1_empty(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "pooled.csv"
            write_per_class_pooled_csv(path, {"vep-rs": list(self.ROWS)}, self.PINS)
            self.assertEqual({r["adj_f1"] for r in csv.DictReader(path.open(encoding="utf-8"))}, {""})


class PerClassModeArgTests(unittest.TestCase):
    """Mode selection and the per-mode required-argument check."""

    def test_default_mode_is_compare(self) -> None:
        args = parse_args(
            ["--perl-dir", "p", "--rust-dir", "r", "--report-dir", "o"]
        )
        self.assertEqual(args.mode, "compare")
        require_mode_args(args)  # must not raise

    def test_by_consequence_class_flag_selects_the_mode(self) -> None:
        args = parse_args(
            [
                "--by-consequence-class",
                "--perl-dir", "p",
                "--discordant-in", "d.tsv",
                "--dataset", "s01",
                "--class-csv", "c.csv",
                "--sweep-id", "sweep-20260101T000000Z",
                "--binary-md5", "deadbeef",
            ]
        )
        self.assertEqual(args.mode, "by-consequence-class")
        require_mode_args(args)

    def test_per_class_mode_requires_the_pins(self) -> None:
        """A per-class CSV without a sweep and binary pin cannot be traced to a run."""
        args = parse_args(
            ["--by-consequence-class", "--perl-dir", "p", "--discordant-in", "d.tsv",
             "--dataset", "s01", "--class-csv", "c.csv"]
        )
        with self.assertRaises(SystemExit) as ctx:
            require_mode_args(args)
        self.assertIn("--sweep-id", str(ctx.exception))
        self.assertIn("--binary-md5", str(ctx.exception))

    def test_per_class_mode_does_not_require_rust_dir(self) -> None:
        """It derives from archived artefacts, so no engine output is needed."""
        self.assertNotIn("rust_dir", MODE_REQUIRED["by-consequence-class"])

    def test_missing_required_argument_names_the_mode_and_the_flag(self) -> None:
        args = parse_args(["--by-consequence-class", "--perl-dir", "p"])
        with self.assertRaises(SystemExit) as ctx:
            require_mode_args(args)
        message = str(ctx.exception)
        self.assertIn("by-consequence-class", message)
        self.assertIn("--discordant-in", message)

    def test_compare_mode_still_requires_all_three_directories(self) -> None:
        args = parse_args(["--perl-dir", "p"])
        with self.assertRaises(SystemExit) as ctx:
            require_mode_args(args)
        self.assertIn("--rust-dir", str(ctx.exception))

    def test_pool_mode_requires_only_the_two_csv_paths(self) -> None:
        args = parse_args(
            ["--pool-consequence-class", "--class-csv", "c.csv", "--pooled-csv", "p.csv",
             "--per-term", "shapes/per_term.tsv"]
        )
        self.assertEqual(args.mode, "pool-consequence-class")
        self.assertEqual(args.per_term, "shapes/per_term.tsv")
        self.assertEqual(args.per_term_engine, "vep-rs")
        require_mode_args(args)

    def test_the_two_mode_flags_are_mutually_exclusive(self) -> None:
        with self.assertRaises(SystemExit):
            parse_args(["--by-consequence-class", "--pool-consequence-class"])

    def test_engine_label_is_constrained_to_the_published_spellings(self) -> None:
        """`vep-rs` with a hyphen is the spelling consumers filter on.

        A free-text engine label would let a typo emit a CSV with zero rows under
        the expected name, which reads as a missing measurement rather than a
        mislabelled one.
        """
        with self.assertRaises(SystemExit):
            parse_args(["--engine", "vep_rs"])
        self.assertEqual(parse_args(["--engine", "vep-rs"]).engine, "vep-rs")

    def test_min_tuples_defaults_to_the_pre_committed_threshold(self) -> None:
        self.assertEqual(parse_args([]).min_tuples, PER_CLASS_MIN_TUPLES)


if __name__ == "__main__":
    unittest.main(verbosity=2)
