#!/usr/bin/env python3
"""Unit tests for intended-divergence filter helpers in compare_vep_outputs.

Run with: ``python3 -m pytest -q scripts/concordance/test_compare_vep_outputs.py``
Or: ``python3 scripts/concordance/test_compare_vep_outputs.py``
"""

from __future__ import annotations

import csv
import dataclasses
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path
from unittest import mock

import compare_vep_outputs as cvo
from compare_vep_outputs import (
    CLASSIFICATION_CHAIN,
    EXCLUDING_RULES,
    EXCLUSION_REGISTRY,
    MODE_REQUIRED,
    ExclusionRule,
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


class KeyFileOrderTests(unittest.TestCase):
    """`compare_sorted_key_files` walks both key files forward once and takes each to
    ascend with every key once.

    Where a key sorts below the key before it, a walk that went on would step past a key
    both files hold and count it as missing and as extra at once, which lowers F1 with no
    error raised. The walk stops there instead, naming the file and the line.
    """

    def _compare(self, perl: list[str], rust: list[str]) -> tuple:
        from compare_vep_outputs import compare_sorted_key_files

        with tempfile.TemporaryDirectory() as d:
            self.perl_path, self.rust_path = Path(d) / "perl.keys", Path(d) / "rust.keys"
            self.perl_path.write_text("".join(k + "\n" for k in perl), encoding="utf-8")
            self.rust_path.write_text("".join(k + "\n" for k in rust), encoding="utf-8")
            return compare_sorted_key_files(
                self.perl_path,
                self.rust_path,
                file_name="f",
                discordant_out=None,
                track_consequence_buckets=False,
            )

    def test_a_key_below_the_key_before_it_in_the_first_file_names_that_file_and_line(self) -> None:
        with self.assertRaises(ValueError) as raised:
            self._compare(["a1", "c3", "b2", "d4"], ["a1", "b2", "c3", "d4"])
        message = str(raised.exception)
        self.assertIn(f"{self.perl_path}: line 3:", message)
        self.assertNotIn(str(self.rust_path), message)

    def test_a_key_below_the_key_before_it_in_the_second_file_names_that_file_and_line(self) -> None:
        with self.assertRaises(ValueError) as raised:
            self._compare(["a1", "b2", "c3", "d4", "e5"], ["a1", "b2", "c3", "e5", "d4"])
        message = str(raised.exception)
        self.assertIn(f"{self.rust_path}: line 5:", message)
        self.assertNotIn(str(self.perl_path), message)

    def test_a_key_out_of_order_after_a_key_only_the_first_file_holds_names_that_file_and_line(self) -> None:
        with self.assertRaises(ValueError) as raised:
            self._compare(["a1", "b2", "d4", "c3"], ["a1", "x9"])
        self.assertIn(f"{self.perl_path}: line 4:", str(raised.exception))

    def test_a_key_out_of_order_after_a_key_only_the_second_file_holds_names_that_file_and_line(self) -> None:
        with self.assertRaises(ValueError) as raised:
            self._compare(["a1", "x9"], ["a1", "b2", "d4", "c3"])
        self.assertIn(f"{self.rust_path}: line 4:", str(raised.exception))

    def test_a_key_written_twice_names_the_line_of_its_second_copy(self) -> None:
        with self.assertRaises(ValueError) as raised:
            self._compare(["a1", "b2", "b2", "c3"], ["a1", "b2", "c3"])
        self.assertIn(f"{self.perl_path}: line 3:", str(raised.exception))

    def test_files_that_ascend_are_counted_to_their_last_keys(self) -> None:
        counts = self._compare(["a1", "b2", "d4", "e5", "f6"], ["a1", "c3", "d4"])
        self.assertEqual(counts[:5], (5, 3, 2, 3, 1))


class RunStoppedOnKeyOrderTests(unittest.TestCase):
    """A compare run that stops on a key file out of order keeps the files it stopped on.

    The program's own sort writes every key file in order, so each run here that must stop
    puts a stand-in `sort` first on the path, which writes the keys descending from a given
    call on. The stopped run leaves the temporary directory of its key files in place and
    names it in the error, and the discordance report written so far is renamed with a
    `.partial` suffix, so no file of the report directory carries a name read as a result.
    Every other ending removes the temporary directory.
    """

    COLS = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
    RESULT_NAMES = ("discordant.tsv", "discordant_open.tsv", "summary.json", "summary.md")

    STAND_IN_SORT = """#!/bin/sh
n=$(( $(cat "$STAND_IN_SORT_CALLS" 2>/dev/null || echo 0) + 1 ))
echo "$n" > "$STAND_IN_SORT_CALLS"
if [ "$n" -ge "$STAND_IN_SORT_FROM" ]; then
    if [ "$STAND_IN_SORT_MODE" = fail ]; then
        cat > /dev/null
        echo "stand-in sort fails" >&2
        exit 3
    fi
    exec "{sort}" "$@" -r
fi
exec "{sort}" "$@"
"""

    @staticmethod
    def _line(loc: str, feature: str) -> str:
        return f"v1\t{loc}\t-\tg1\t{feature}\tTranscript\tintron_variant\n"

    def setUp(self) -> None:
        sort = shutil.which("sort")
        if sort is None:
            self.skipTest("no sort program to stand in for")
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        for sub in ("perl", "rust", "report", "bin"):
            (self.tmp / sub).mkdir()
        shared = self._line("21:100", "tx1") + self._line("21:200", "tx2") + self._line("21:300", "tx3")
        (self.tmp / "perl" / "output.txt").write_text(self.COLS + shared, encoding="utf-8")
        (self.tmp / "rust" / "output.txt").write_text(
            self.COLS + shared + self._line("21:900", "tx9"), encoding="utf-8"
        )
        stand_in = self.tmp / "bin" / "sort"
        stand_in.write_text(self.STAND_IN_SORT.replace("{sort}", sort), encoding="utf-8")
        stand_in.chmod(0o755)
        self.report = self.tmp / "report"

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _run(self, *extra: str, mode: str | None = None, start: int = 1) -> subprocess.CompletedProcess:
        import compare_vep_outputs

        env = dict(os.environ)
        if mode is not None:
            env["PATH"] = f"{self.tmp / 'bin'}{os.pathsep}{env['PATH']}"
            env["STAND_IN_SORT_MODE"] = mode
            env["STAND_IN_SORT_FROM"] = str(start)
            env["STAND_IN_SORT_CALLS"] = str(self.tmp / "sort.calls")
        return subprocess.run(
            [
                sys.executable,
                compare_vep_outputs.__file__,
                "--perl-dir", str(self.tmp / "perl"),
                "--rust-dir", str(self.tmp / "rust"),
                "--report-dir", str(self.report),
                "--assembly", "GRCh37",
                *extra,
            ],
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )

    def _kept_directory(self) -> Path:
        kept = sorted(self.report.glob("compare_vep_tmp_*"))
        self.assertEqual(len(kept), 1, f"temporary directories left in the report directory: {kept}")
        return kept[0]

    def _descending_keys(self, side: str, transform=None) -> list[str]:
        from compare_vep_outputs import iter_vep_key_lines

        keys = set(iter_vep_key_lines(self.tmp / side / "output.txt", consequence_transform=transform))
        return sorted(keys, reverse=True)

    def test_a_run_stopped_on_a_key_file_out_of_order_keeps_its_key_files_and_names_their_directory(self) -> None:
        run = self._run(mode="descending")
        self.assertNotEqual(run.returncode, 0)
        error = run.stderr.strip().splitlines()[-1]
        kept = self._kept_directory()
        self.assertIn(f"{kept / 'output.txt.perl.keys.sorted.txt'}: line 2:", error)
        self.assertIn(f"the key files are kept in {kept}", error)
        for side in ("perl", "rust"):
            keys = (kept / f"output.txt.{side}.keys.sorted.txt").read_text(encoding="utf-8").splitlines()
            self.assertEqual(keys, self._descending_keys(side))

    def test_a_run_stopped_on_a_key_file_out_of_order_leaves_no_file_under_a_result_name(self) -> None:
        run = self._run(mode="descending")
        self.assertNotEqual(run.returncode, 0)
        for name in self.RESULT_NAMES:
            self.assertFalse((self.report / name).exists(), name)
        partial = self.report / "discordant.tsv.partial"
        self.assertTrue(partial.read_text(encoding="utf-8").startswith("file_name\tsource\tlocation\t"))
        self.assertIn(f"the rows written before the stop in {partial}", run.stderr.strip().splitlines()[-1])

    def test_a_stop_in_the_context_collapsed_comparison_keeps_its_key_files_too(self) -> None:
        from compare_vep_outputs import collapse_transcript_context_terms

        run = self._run("--discordant-tsv", "custom.tsv", mode="descending", start=3)
        self.assertNotEqual(run.returncode, 0)
        error = run.stderr.strip().splitlines()[-1]
        kept = self._kept_directory()
        self.assertIn(f"{kept / 'output.txt.perl.keys.collapsed.sorted.txt'}: line 2:", error)
        self.assertIn(f"the key files are kept in {kept}", error)
        for side in ("perl", "rust"):
            keys = (kept / f"output.txt.{side}.keys.collapsed.sorted.txt").read_text(encoding="utf-8").splitlines()
            self.assertEqual(keys, self._descending_keys(side, collapse_transcript_context_terms))
        for name in ("custom.tsv", "custom_open.tsv", "summary.json", "summary.md"):
            self.assertFalse((self.report / name).exists(), name)
        self.assertTrue((self.report / "custom.tsv.partial").exists())

    def test_a_run_that_ends_well_removes_its_temporary_directory(self) -> None:
        run = self._run()
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(sorted(p.name for p in self.report.iterdir()), sorted(self.RESULT_NAMES))

    def test_a_run_that_ends_on_another_error_removes_its_temporary_directory(self) -> None:
        run = self._run(mode="fail", start=2)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("sort -u failed", run.stderr.strip().splitlines()[-1])
        self.assertEqual([p.name for p in self.report.iterdir()], ["discordant.tsv"])

    def test_a_value_error_from_elsewhere_in_the_comparison_removes_its_temporary_directory(self) -> None:
        import compare_vep_outputs

        # main() lifts the interpreter's bound on an integer's digits for the whole process, so
        # a test that calls it in this process puts the bound back.
        self.addCleanup(sys.set_int_max_str_digits, sys.get_int_max_str_digits())
        argv = [
            "compare_vep_outputs.py",
            "--perl-dir", str(self.tmp / "perl"),
            "--rust-dir", str(self.tmp / "rust"),
            "--report-dir", str(self.report),
            "--assembly", "GRCh37",
        ]
        with (
            mock.patch.object(sys, "argv", argv),
            mock.patch.object(
                compare_vep_outputs, "compare_sorted_key_files", side_effect=ValueError("another error")
            ),
            self.assertRaises(ValueError) as raised,
        ):
            compare_vep_outputs.main()
        self.assertIs(type(raised.exception), ValueError)
        self.assertEqual([p.name for p in self.report.iterdir()], ["discordant.tsv"])


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
    inflates F1 and nothing else in the pipeline would notice. A key that omits
    `feature_type` while the scored tuple includes it collapses two rows differing
    only by feature type onto one key, and a nested loop over that key counts their
    cross product: 4 exclusions from 2 real tuples.
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
            return filter_snp_indel_intended_divergences(p, reference_release="115.2")

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
            return filter_snp_indel_intended_divergences(p, reference_release="115.2")

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


class OpenDiscordantTsvTests(unittest.TestCase):
    """The open subset of a discordance report: every row no excluding rule removed.

    Both rows of an open pair and every unpaired row survive; both rows of an excluded
    pair go. The file is the adjusted residual enumerated, so its row count must equal
    the one-sided total the adjusted F1 is computed on.
    """

    HEADER = (
        "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set"
    )

    @staticmethod
    def _row(source: str, loc: str, allele: str, feature: str, csq: str) -> str:
        return "\t".join(("f", source, loc, allele, feature, "Transcript", _norm(csq)))

    def _rows(self) -> list[str]:
        return [
            # Covered splice_region pair: excluded.
            self._row("missing_in_rust", "21:100", "A", "ENST1", "missense_variant"),
            self._row("extra_in_rust", "21:100", "A", "ENST1", "missense_variant,splice_region_variant"),
            # Start co-emission pair: excluded.
            self._row("missing_in_rust", "21:200", "-", "ENST2", "start_lost,start_retained_variant"),
            self._row("extra_in_rust", "21:200", "-", "ENST2", "start_retained_variant"),
            # PPT addition: classified, not excluded, so open.
            self._row("missing_in_rust", "21:300", "T", "ENST3", "intron_variant"),
            self._row("extra_in_rust", "21:300", "T", "ENST3", "intron_variant,splice_polypyrimidine_tract_variant"),
            # Unpaired rows, one per side.
            self._row("missing_in_rust", "21:400", "G", "ENST4", "downstream_gene_variant"),
            self._row("extra_in_rust", "21:500", "C", "ENST5", "upstream_gene_variant"),
        ]

    def _run(self, rows: list[str], open_out: bool = True):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "discordant.tsv"
            p.write_text(self.HEADER + "\n" + "".join(r + "\n" for r in rows), encoding="utf-8")
            o = Path(d) / "discordant_open.tsv"
            result = filter_snp_indel_intended_divergences(
                p, open_out=o if open_out else None, reference_release="115.2"
            )
            return result, o.read_text(encoding="utf-8") if o.exists() else None

    def test_open_rows_are_the_unexcluded_pairs_plus_the_unpaired_rows(self) -> None:
        rows = self._rows()
        (ep, er, _), text = self._run(rows)
        self.assertEqual((ep, er), (2, 2))
        assert text is not None
        self.assertEqual(text.splitlines(), [self.HEADER] + rows[4:])

    def test_row_count_is_the_adjusted_one_sided_total(self) -> None:
        rows = self._rows()
        (ep, er, _), text = self._run(rows)
        assert text is not None
        n_perl = sum("\tmissing_in_rust\t" in r for r in rows)
        n_rust = sum("\textra_in_rust\t" in r for r in rows)
        self.assertEqual(len(text.splitlines()) - 1, (n_perl - ep) + (n_rust - er))

    def test_file_order_and_layout_follow_the_source_report(self) -> None:
        """A row is written back byte for byte, in the position it had, so the open
        file is a line-subset of the report and greps the same way."""
        rows = self._rows()
        rows.reverse()
        _, text = self._run(rows)
        assert text is not None
        source_lines = [self.HEADER] + rows
        positions = [source_lines.index(line) for line in text.splitlines()]
        self.assertEqual(positions, sorted(positions))
        self.assertEqual(len(positions), 5)

    def test_nothing_is_written_without_a_path(self) -> None:
        (ep, er, buckets), text = self._run(self._rows(), open_out=False)
        self.assertIsNone(text)
        self.assertEqual((ep, er), (2, 2))
        self.assertEqual(sum(v["excluded"] for v in buckets.values()), 2)

    def test_a_fully_masked_report_leaves_only_the_header(self) -> None:
        _, text = self._run(self._rows()[:4])
        self.assertEqual(text, self.HEADER + "\n")

    def test_a_header_only_report_gives_a_header_only_open_file(self) -> None:
        (ep, er, _), text = self._run([])
        self.assertEqual((ep, er), (0, 0))
        self.assertEqual(text, self.HEADER + "\n")


class CompareModeOpenTsvEndToEndTests(unittest.TestCase):
    """Compare mode writes ``discordant_open.tsv`` beside ``discordant.tsv`` and the
    two agree with ``summary.json``: raw rows with the raw one-sided total, open
    rows with the adjusted one."""

    COLS = (
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
    )

    @staticmethod
    def _line(loc: str, allele: str, feature: str, csq: str) -> str:
        return f"v\t{loc}\t{allele}\tG\t{feature}\tTranscript\t{csq}\n"

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        for sub in ("perl", "rust", "report"):
            (self.tmp / sub).mkdir()
        shared = self._line("21:1", "A", "ENST0", "intron_variant")
        (self.tmp / "perl" / "output.txt").write_text(
            self.COLS
            + shared
            + self._line("21:100", "A", "ENST1", "missense_variant")
            + self._line("21:300", "T", "ENST3", "intron_variant")
            + self._line("21:400", "G", "ENST4", "downstream_gene_variant"),
            encoding="utf-8",
        )
        (self.tmp / "rust" / "output.txt").write_text(
            self.COLS
            + shared
            + self._line("21:100", "A", "ENST1", "missense_variant,splice_region_variant")
            + self._line("21:300", "T", "ENST3", "intron_variant,splice_polypyrimidine_tract_variant"),
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _run(self, *extra: str) -> dict:
        import json
        import sys
        from compare_vep_outputs import main

        self.addCleanup(sys.set_int_max_str_digits, sys.get_int_max_str_digits())
        argv = sys.argv
        sys.argv = [
            "compare_vep_outputs.py",
            "--perl-dir", str(self.tmp / "perl"),
            "--rust-dir", str(self.tmp / "rust"),
            "--report-dir", str(self.tmp / "report"),
            "--assembly", "GRCh37",
            *extra,
        ]
        try:
            self.assertEqual(main(), 0)
        finally:
            sys.argv = argv
        return json.loads((self.tmp / "report" / "summary.json").read_text(encoding="utf-8"))

    def test_open_file_rows_equal_the_adjusted_one_sided_total(self) -> None:
        adj = self._run()["aggregate"]["adjusted"]
        raw_lines = (self.tmp / "report" / "discordant.tsv").read_text().splitlines()
        open_lines = (self.tmp / "report" / "discordant_open.tsv").read_text().splitlines()
        self.assertEqual(adj["excluded_perl_tuples"], 1)
        self.assertEqual(len(raw_lines) - 1, 5)
        self.assertEqual(
            len(open_lines) - 1,
            (adj["perl_tuple_count"] - adj["intersection_count"])
            + (adj["rust_tuple_count"] - adj["intersection_count"]),
        )
        self.assertEqual(open_lines[0], raw_lines[0])
        self.assertTrue(set(open_lines) <= set(raw_lines))
        self.assertFalse(any("\t21:100\t" in line for line in open_lines))
        self.assertEqual(sum("\t21:300\t" in line for line in open_lines), 2)
        self.assertEqual(sum("\t21:400\t" in line for line in open_lines), 1)

    def test_open_file_is_named_after_the_discordant_tsv_flag(self) -> None:
        self._run("--discordant-tsv", "custom.tsv")
        self.assertTrue((self.tmp / "report" / "custom.tsv").exists())
        self.assertTrue((self.tmp / "report" / "custom_open.tsv").exists())
        self.assertFalse((self.tmp / "report" / "discordant_open.tsv").exists())


class CanonicalContigsFromCacheTests(unittest.TestCase):
    """Under an assembly that is not GRCh37 or GRCh38 the canonical set is the one
    ``--canonical-contigs`` names, the cache's ``contig_order`` as run_concordance.sh
    passes it; bare, the flag keeps the human set for the human assemblies and is
    off, with a warning, for any other."""

    COLS = (
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
    )
    CONTIG_ORDER = ("2L", "2R", "3L", "3R", "4", "X", "Y", "mitochondrion_genome")

    @staticmethod
    def _line(loc: str, feature: str, csq: str) -> str:
        return f"v\t{loc}\tA\tG\t{feature}\tTranscript\t{csq}\n"

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        for sub in ("perl", "rust", "report", "cache"):
            (self.tmp / sub).mkdir()
        shared = (
            self._line("2L:100", "FBtr0", "missense_variant")
            + self._line("X:200", "FBtr1", "intron_variant")
            + self._line("mitochondrion_genome:300", "FBtr2", "synonymous_variant")
        )
        (self.tmp / "perl" / "output.txt").write_text(
            self.COLS + shared + self._line("rDNA:400", "FBtr3", "intron_variant"),
            encoding="utf-8",
        )
        (self.tmp / "rust" / "output.txt").write_text(
            self.COLS
            + shared
            + self._line("Unmapped_Scaffold_8_D1580_D1567:500", "FBtr4", "intron_variant"),
            encoding="utf-8",
        )
        (self.tmp / "cache" / "info.json").write_text(
            json.dumps({"species": "drosophila_melanogaster", "assembly": "BDGP6.54",
                        "contig_order": self.CONTIG_ORDER}),
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _run(self, *extra: str) -> tuple[dict, str]:
        import contextlib
        import io
        import sys
        from compare_vep_outputs import main

        self.addCleanup(sys.set_int_max_str_digits, sys.get_int_max_str_digits())
        argv = sys.argv
        sys.argv = [
            "compare_vep_outputs.py",
            "--perl-dir", str(self.tmp / "perl"),
            "--rust-dir", str(self.tmp / "rust"),
            "--report-dir", str(self.tmp / "report"),
            *extra,
        ]
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                self.assertEqual(main(), 0)
        finally:
            sys.argv = argv
        summary = json.loads((self.tmp / "report" / "summary.json").read_text(encoding="utf-8"))
        return summary, out.getvalue()

    def test_canonical_contigs_come_from_the_caches_contig_order_for_a_non_human_assembly(
        self,
    ) -> None:
        # The list run_concordance.sh derives from the cache's info.json.
        info = json.loads((self.tmp / "cache" / "info.json").read_text(encoding="utf-8"))
        contig_list = ",".join(info["contig_order"])
        summary, out = self._run("--assembly", "BDGP6.54", "--canonical-contigs", contig_list)
        self.assertTrue(summary["canonical_contigs_filter_active"])
        self.assertEqual(summary["canonical_contigs_set"], sorted(self.CONTIG_ORDER))
        raw = summary["aggregate"]
        # The rDNA and scaffold rows are outside the set: three tuples each side, all shared.
        self.assertEqual(
            (raw["perl_tuple_count"], raw["rust_tuple_count"], raw["intersection_count"]),
            (3, 3, 3),
        )
        self.assertNotIn("WARNING: --canonical-contigs", out)

        # The same set from a file, one name per line, a comment line ignored.
        listing = self.tmp / "contigs.txt"
        listing.write_text("# the cache's karyotype\n" + "\n".join(self.CONTIG_ORDER) + "\n")
        from_file, _ = self._run("--assembly", "BDGP6.54", "--canonical-contigs", str(listing))
        self.assertEqual(from_file["canonical_contigs_set"], sorted(self.CONTIG_ORDER))
        self.assertEqual(from_file["aggregate"]["intersection_count"], 3)

        # Bare under a non-human assembly: no set is derivable, the filter is off and
        # the one-sided rows count, with the warning naming the assembly and the remedy.
        bare, out = self._run("--assembly", "BDGP6.54")
        self.assertFalse(bare["canonical_contigs_filter_active"])
        self.assertIsNone(bare["canonical_contigs_set"])
        self.assertEqual(
            (bare["aggregate"]["perl_tuple_count"], bare["aggregate"]["rust_tuple_count"]),
            (4, 4),
        )
        self.assertIn("--assembly is BDGP6.54, which names no canonical set", out)
        self.assertIn("--canonical-contigs <file or list>", out)

        # Bare under a human assembly keeps the human set.
        from compare_vep_outputs import CANONICAL_CONTIGS

        human, _ = self._run("--assembly", "GRCh38")
        self.assertEqual(human["canonical_contigs_set"], sorted(CANONICAL_CONTIGS))

    def test_a_canonical_contigs_value_naming_no_contig_is_a_usage_error(self) -> None:
        from compare_vep_outputs import parse_canonical_contigs

        self.assertEqual(parse_canonical_contigs("2L, X mitochondrion_genome"),
                         frozenset({"2L", "X", "mitochondrion_genome"}))
        empty = self.tmp / "empty.txt"
        empty.write_text("# nothing\n\n")
        with self.assertRaises(SystemExit):
            parse_canonical_contigs(str(empty))
        with self.assertRaises(SystemExit):
            parse_canonical_contigs(",")


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


class VepTotalsCacheTests(unittest.TestCase):
    """The per-term VEP totals cache: keyed on the reference files and the contig
    filter, served only on an exact match, recomputed and rewritten otherwise."""

    COLS = (
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
    )

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        (self.tmp / "perl").mkdir()
        self.gt = self.tmp / "perl" / "output.txt"
        self.gt.write_text(
            self.COLS
            + "v\t21:100\tA\tG\tENST1\tTranscript\tmissense_variant\n"
            + "v\t21:200\tT\tG\tENST2\tTranscript\tintron_variant,splice_region_variant\n"
            + "v\tKI270728.1:5\tC\tG\tENST3\tTranscript\tintron_variant\n",
            encoding="utf-8",
        )
        self.disc = self.tmp / "discordant.tsv"
        self.disc.write_text(
            "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n"
            "output.txt\tmissing_in_rust\t21:100\tA\tENST1\tTranscript\tmissense_variant\n",
            encoding="utf-8",
        )
        self.cache = self.tmp / "cache" / "s01.csv"

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _run(self, *extra: str) -> tuple[str, str]:
        import contextlib
        import io
        import sys
        from compare_vep_outputs import main

        class_csv = self.tmp / "class.csv"
        self.addCleanup(sys.set_int_max_str_digits, sys.get_int_max_str_digits())
        argv = sys.argv
        sys.argv = [
            "compare_vep_outputs.py", "--by-consequence-class",
            "--perl-dir", str(self.tmp / "perl"), "--discordant-in", str(self.disc),
            "--dataset", "s01", "--class-csv", str(class_csv),
            "--sweep-id", "sweep-x", "--binary-md5", "abcd1234",
            "--assembly", "GRCh37", "--vep-totals-cache", str(self.cache), *extra,
        ]
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                self.assertEqual(main(), 0)
        finally:
            sys.argv = argv
        return out.getvalue(), class_csv.read_text(encoding="utf-8")

    def test_fingerprint_reads_name_size_mtime_and_contig_filter(self) -> None:
        from compare_vep_outputs import CANONICAL_CONTIGS, vep_totals_fingerprint

        fp = vep_totals_fingerprint([self.gt], CANONICAL_CONTIGS)
        st = self.gt.stat()
        self.assertIn(f"output.txt:{st.st_size}:{st.st_mtime_ns}", fp)
        self.assertNotEqual(fp, vep_totals_fingerprint([self.gt], None))
        os.utime(self.gt, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
        self.assertNotEqual(fp, vep_totals_fingerprint([self.gt], CANONICAL_CONTIGS))

    def test_read_returns_none_on_missing_stale_or_corrupt_cache(self) -> None:
        from compare_vep_outputs import read_vep_totals_cache, write_vep_totals_cache

        self.assertIsNone(read_vep_totals_cache(self.cache, "fp"))
        write_vep_totals_cache(self.cache, "fp", Counter({"missense_variant": 3}))
        self.assertEqual(read_vep_totals_cache(self.cache, "fp"), Counter({"missense_variant": 3}))
        self.assertIsNone(read_vep_totals_cache(self.cache, "other"))
        self.cache.write_text("# vep_totals_fingerprint=fp\nterm,vep_tuples\nx,notanumber\n")
        self.assertIsNone(read_vep_totals_cache(self.cache, "fp"))
        self.assertFalse(self.cache.with_name("s01.csv.tmp").exists())

    def test_second_run_is_served_from_the_cache_with_identical_output(self) -> None:
        out1, csv1 = self._run()
        self.assertIn("cache miss, wrote", out1)
        self.assertTrue(self.cache.is_file())
        cached = self.cache.read_text(encoding="utf-8")
        # The alt-contig tuple is outside the canonical filter, so intron_variant is 1.
        self.assertIn("intron_variant,1\n", cached)
        self.assertIn("missense_variant,1\n", cached)
        self.assertIn("splice_region_variant,1\n", cached)
        out2, csv2 = self._run()
        self.assertIn("cache hit", out2)
        self.assertEqual(csv1, csv2)

    def test_a_touched_reference_invalidates_the_cache(self) -> None:
        self._run()
        st = self.gt.stat()
        os.utime(self.gt, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
        out, _ = self._run()
        self.assertIn("cache miss", out)

    def test_a_different_contig_filter_invalidates_the_cache(self) -> None:
        """Totals counted without the filter carry the alt-contig tuple; serving
        them to a filtered run would inflate every term's VEP total."""
        self._run()
        out, _ = self._run("--no-canonical-contigs")
        self.assertIn("cache miss", out)
        self.assertIn("intron_variant,2\n", self.cache.read_text(encoding="utf-8"))

    def test_no_cache_flag_neither_reads_nor_writes(self) -> None:
        import contextlib
        import io
        import sys
        from compare_vep_outputs import main

        self.addCleanup(sys.set_int_max_str_digits, sys.get_int_max_str_digits())
        argv = sys.argv
        sys.argv = [
            "compare_vep_outputs.py", "--by-consequence-class",
            "--perl-dir", str(self.tmp / "perl"), "--discordant-in", str(self.disc),
            "--dataset", "s01", "--class-csv", str(self.tmp / "class.csv"),
            "--sweep-id", "sweep-x", "--binary-md5", "abcd1234", "--assembly", "GRCh37",
        ]
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                self.assertEqual(main(), 0)
        finally:
            sys.argv = argv
        self.assertNotIn("VEP per-term totals: cache", out.getvalue())
        self.assertFalse(self.cache.exists())


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

    def test_a_missing_input_vcf_or_input_dir_is_refused_before_any_pairing(self) -> None:
        """A path given to `--input-vcf` that is not a file, or to `--input-dir` that is not a
        directory, is a usage error (exit code 2) before any output is read, so a mistyped path
        cannot run the comparison without the input it names."""
        import contextlib
        import io

        with tempfile.TemporaryDirectory() as d:
            for flag, path in (("--input-vcf", Path(d) / "absent.vcf"), ("--input-dir", Path(d) / "absent"),
                               ("--input-vcf", Path(d)), ("--input-dir", Path(d) / "absent.vcf")):
                if flag == "--input-dir" and path.name == "absent.vcf":
                    path.write_text("", encoding="utf-8")
                with self.subTest(flag=flag, path=path.name), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as raised:
                    parse_args(["--perl-dir", "p", "--rust-dir", "r", "--report-dir", "o", flag, str(path)])
                self.assertEqual(raised.exception.code, 2)

    def test_min_tuples_defaults_to_the_pre_committed_threshold(self) -> None:
        self.assertEqual(parse_args([]).min_tuples, PER_CLASS_MIN_TUPLES)


class SvTranscriptSelectionTests(unittest.TestCase):
    """The structural-variant comparator's transcript-selection class on a whole-genome
    input: a record above --max_sv_size that Perl marks skipped is annotated against
    its batch's loaded regions alone, so the transcripts vep-rs names beyond Perl's are
    set aside, under the two guards (Perl named at least one transcript; vep-rs named
    more), and a Perl-only transcript is set aside only when vep-rs named at least as
    many. The scope map is what the sorted merge collects: every transcript each engine
    named on the record, matched rows included.
    """

    HEADER = (
        "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set"
    )
    GIANT = "18:62909529-80261528"  # 17.35 Mb, the span of a 1000 Genomes <DEL>
    SMALL = "18:62909529-62919528"  # 10 kb

    def _run(
        self,
        rows: list[tuple[str, ...]],
        scope: dict[tuple[str, str], tuple[set[str], set[str]]] | None,
    ):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "discordant.tsv"
            o = Path(d) / "discordant_open.tsv"
            with p.open("w", encoding="utf-8") as f:
                f.write(self.HEADER + "\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")
            result = filter_snp_indel_intended_divergences(
                p, open_out=o, scope_features=scope, reference_release="115.2"
            )
            open_rows = o.read_text(encoding="utf-8").splitlines()[1:]
            return result, open_rows

    def _row(self, source: str, loc: str, feature: str, csq: str = "transcript_ablation") -> tuple[str, ...]:
        return ("f", source, loc, "deletion", feature, "Transcript", _norm(csq))

    def test_vep_rs_surplus_transcripts_on_a_giant_record_are_set_aside(self) -> None:
        # Perl named 4 transcripts (all matched, so none is in discordant.tsv); vep-rs
        # named those 4 plus 3 more.
        scope = {(self.GIANT, "deletion"): ({"T1", "T2", "T3", "T4"}, {"T1", "T2", "T3", "T4", "T5", "T6", "T7"})}
        rows = [self._row("extra_in_rust", self.GIANT, t) for t in ("T5", "T6", "T7")]
        (ep, er, buckets), open_rows = self._run(rows, scope)
        self.assertEqual((ep, er), (0, 3))
        self.assertEqual(buckets["sv_transcript_selection"], {"count": 3, "excluded": 3, "perl": 0, "rust": 3})
        self.assertEqual(open_rows, [], "every surplus row leaves the adjusted residual")

    def test_a_pair_on_a_transcript_both_engines_named_stays_charged(self) -> None:
        scope = {(self.GIANT, "deletion"): ({"T1", "T2"}, {"T1", "T2", "T3"})}
        rows = [
            self._row("missing_in_rust", self.GIANT, "T1", "intron_variant"),
            self._row("extra_in_rust", self.GIANT, "T1", "intron_variant,stop_lost"),
            self._row("extra_in_rust", self.GIANT, "T3"),
        ]
        (ep, er, buckets), open_rows = self._run(rows, scope)
        self.assertEqual((ep, er), (0, 1), "T3 is set aside; the T1 pair is a consequence difference, not a selection one")
        self.assertEqual(len(open_rows), 2)
        self.assertEqual(buckets["sv_transcript_selection"], {"count": 3, "excluded": 1, "perl": 0, "rust": 1},
                         "all three rows sit on the record in scope; one is the class")

    def test_direction_guard_keeps_a_vep_rs_under_annotation_charged(self) -> None:
        # Perl named 5 transcripts, vep-rs 2 of them: Perl's extra rows are not the class.
        scope = {(self.GIANT, "deletion"): ({"T1", "T2", "T3", "T4", "T5"}, {"T1", "T2"})}
        rows = [self._row("missing_in_rust", self.GIANT, t) for t in ("T3", "T4", "T5")]
        (ep, er, buckets), open_rows = self._run(rows, scope)
        self.assertEqual((ep, er), (0, 0))
        self.assertEqual(len(open_rows), 3)
        self.assertEqual(buckets["sv_transcript_selection"], {"count": 3, "excluded": 0, "perl": 0, "rust": 0})

    def test_perl_only_transcripts_set_aside_when_vep_rs_named_at_least_as_many(self) -> None:
        # Perl's batch loaded a region vep-rs's set does not reach, and vep-rs named as many.
        scope = {(self.GIANT, "deletion"): ({"T1", "T2", "T9"}, {"T1", "T2", "T3"})}
        rows = [self._row("missing_in_rust", self.GIANT, "T9"), self._row("extra_in_rust", self.GIANT, "T3")]
        (ep, er, buckets), open_rows = self._run(rows, scope)
        self.assertEqual((ep, er), (1, 0), "vep-rs named 3 against Perl's 3, so the Perl-only transcript is the class; vep-rs did not name MORE, so its own extra stays")
        self.assertEqual(open_rows, ["\t".join(self._row("extra_in_rust", self.GIANT, "T3"))])

    def test_a_record_perl_never_annotated_is_outside_the_class(self) -> None:
        scope = {(self.GIANT, "deletion"): (set(), {"T1", "T2"})}
        rows = [self._row("extra_in_rust", self.GIANT, t) for t in ("T1", "T2")]
        (ep, er, buckets), open_rows = self._run(rows, scope)
        self.assertEqual((ep, er), (0, 0))
        self.assertEqual(len(open_rows), 2)

    def test_a_record_below_max_sv_size_never_enters_the_scope(self) -> None:
        # The merge records nothing for it, so the exclusion pass never sees it.
        scope: dict = {}
        rows = [self._row("extra_in_rust", self.SMALL, t) for t in ("T5", "T6")]
        (ep, er, buckets), open_rows = self._run(rows, scope)
        self.assertEqual((ep, er), (0, 0))
        self.assertNotIn("sv_transcript_selection", buckets)
        self.assertEqual(len(open_rows), 2)

    def test_merge_collects_the_scope_from_both_streams(self) -> None:
        from compare_vep_outputs import compare_sorted_key_files

        def key(loc: str, allele: str, feat: str, csq: str) -> str:
            return f"{loc}\t{allele}\t{feat}\tTranscript\t{csq}"

        perl = sorted([key(self.GIANT, "deletion", "T1", "intron_variant"), key(self.GIANT, "deletion", "T2", "intron_variant"),
                       key(self.SMALL, "deletion", "S1", "intron_variant"), key("18:100", "A", "TA", "missense_variant")])
        rust = sorted([key(self.GIANT, "deletion", "T1", "intron_variant"), key(self.GIANT, "deletion", "T3", "transcript_ablation"),
                       key(self.SMALL, "deletion", "S1", "intron_variant"), key("18:200-200", "CPX", "C1", "intron_variant")])
        with tempfile.TemporaryDirectory() as d:
            pp, rp = Path(d) / "p.txt", Path(d) / "r.txt"
            pp.write_text("\n".join(perl) + "\n", encoding="utf-8")
            rp.write_text("\n".join(rust) + "\n", encoding="utf-8")
            scope: dict = {}
            compare_sorted_key_files(pp, rp, file_name="f", discordant_out=None, track_consequence_buckets=False, scope_features=scope)
        self.assertEqual(scope, {
            (self.GIANT, "deletion"): ({"T1", "T2"}, {"T1", "T3"}),
            ("18:200-200", "CPX"): (set(), {"C1"}),
        }, "the giant record and the unsupported type are in scope, the 10 kb deletion and the SNV are not; the matched T1 counts on both sides")

    def test_a_typeless_symbolic_alt_is_in_scope_under_either_spelling(self) -> None:
        """A bare `<CPX>` or `<CTX>` prints its ALT as written and the same ALT under
        `SVTYPE=CPX` prints `CPX` (Parser/VCF.pm:469-481); every spelling is the
        `vep_skip` route, a `<DEL>` below the size limit is not."""
        from compare_vep_outputs import _batch_dependent_annotation_scope

        for allele in ("CPX", "CTX", "<CPX>", "<CTX>"):
            self.assertTrue(_batch_dependent_annotation_scope("21:800-900", allele), allele)
        self.assertFalse(_batch_dependent_annotation_scope("21:800-900", "deletion"))


# The reference release selector (Ensembl VEP 115.2 or 116.2)


class ReferenceRegistryTests(unittest.TestCase):
    """`REFERENCE_REGISTRIES` carries one registry per Ensembl VEP release the comparator
    scores against. Both share the pair rules (`EXCLUSION_REGISTRY`, the paper's, by
    identity); they differ in the one-sided classes."""

    def test_both_releases_share_the_paper_pair_rules_by_identity(self) -> None:
        from compare_vep_outputs import REFERENCE_REGISTRIES, REFERENCE_RELEASES

        self.assertEqual(tuple(REFERENCE_REGISTRIES), REFERENCE_RELEASES)
        for release, registry in REFERENCE_REGISTRIES.items():
            self.assertIs(registry.pair_rules, EXCLUSION_REGISTRY, release)
            self.assertEqual(registry.release, release)

    def test_the_116_2_one_sided_rules_are_the_skipped_record_and_the_thirteen_record_shapes(self) -> None:
        from compare_vep_outputs import (
            BND_BARE_REFERENCE_ALLELE_ROW_BUCKET,
            BND_MATE_32_BIT_BUFFER_BUCKET,
            BND_OWN_ROWS_LOST_BUCKET,
            BND_SYNONYM_MATE_ROWS_LOST_BUCKET,
            CNV_TR_EMPTY_EXPANSION_BUCKET,
            CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET,
            DOT_BREAKEND_PAIR_BUCKET,
            FASTA_NAMED_SLICE_INTERGENIC_BUCKET,
            GVCF_VARIANT_SITE_DROPPED_BUCKET,
            MISSING_ALT_ANNOTATED_BUCKET,
            NON_INTEGER_END_BUCKET,
            NON_NUMERIC_POS_NUMIFIED_BUCKET,
            REFERENCE_116_2_ONE_SIDED_RULES,
            REFERENCE_REGISTRIES,
            REFERENCE_SKIPPED_RECORD_BUCKET,
            STAR_ALONE_SVTYPE_ANNOTATED_BUCKET,
        )

        rules = REFERENCE_REGISTRIES["116.2"].one_sided_rules
        self.assertIs(rules, REFERENCE_116_2_ONE_SIDED_RULES)
        self.assertEqual(
            [r.bucket for r in rules],
            [REFERENCE_SKIPPED_RECORD_BUCKET, BND_OWN_ROWS_LOST_BUCKET,
             BND_SYNONYM_MATE_ROWS_LOST_BUCKET, FASTA_NAMED_SLICE_INTERGENIC_BUCKET,
             BND_BARE_REFERENCE_ALLELE_ROW_BUCKET, BND_MATE_32_BIT_BUFFER_BUCKET,
             GVCF_VARIANT_SITE_DROPPED_BUCKET, DOT_BREAKEND_PAIR_BUCKET, MISSING_ALT_ANNOTATED_BUCKET,
             NON_INTEGER_END_BUCKET, CNV_TR_EMPTY_EXPANSION_BUCKET, STAR_ALONE_SVTYPE_ANNOTATED_BUCKET,
             NON_NUMERIC_POS_NUMIFIED_BUCKET, CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET],
        )
        self.assertEqual([r.excludes for r in rules],
                         [True, False, False, False, True, False, True, True, False, True, False, False, True, False],
                         "this comparator reads no cache and opens the input only to recognise gVCF variant "
                         "sites: the lost-row breakend shapes, the missing ALT, the empty expansion, the "
                         "`*`-only ALT and the dropped `RN=0` allele are judged only by the structural-variant "
                         "comparator")
        self.assertEqual([r.excludes_with_cache for r in rules],
                         [True, True, True, False, True, False, True, True, True, True, True, True, True, True],
                         "the structural-variant comparator sets every class aside but the FASTA shape and the "
                         "buffer a coordinate overflows")
        for bucket in (NON_INTEGER_END_BUCKET, BND_MATE_32_BIT_BUFFER_BUCKET):
            self.assertIn("SVLEN" if bucket == NON_INTEGER_END_BUCKET else "span's end", rules[[r.bucket for r in rules].index(bucket)].definition)
        overflow = rules[5]
        self.assertEqual(overflow.bucket, "bnd_mate_32bit_buffer_intergenic")
        self.assertIn("2147478648", overflow.definition)
        bare = rules[4]
        self.assertEqual(bare.bucket, "bnd_bare_reference_allele_row")
        self.assertIn("StructuralVariationOverlap.pm:64-77", bare.perl_citation)
        self.assertIn("reference allele has no variant effect", bare.definition)
        self.assertIn("set aside", bare.definition)
        for r in rules:
            self.assertTrue(r.definition.endswith(".") and len(r.definition.split()) >= 12, r.bucket)
            self.assertRegex(r.perl_citation, r"\.pm", r.bucket)
        self.assertEqual(len({r.bucket for r in rules}), len(rules))
        self.assertEqual(len({r.taxonomy_class for r in rules}), len(rules))

    def test_the_115_2_registry_describes_the_scope_class_it_implements(self) -> None:
        from compare_vep_outputs import REFERENCE_REGISTRIES, SV_TRANSCRIPT_SELECTION_BUCKET

        (rule,) = REFERENCE_REGISTRIES["115.2"].one_sided_rules
        self.assertEqual(rule.bucket, SV_TRANSCRIPT_SELECTION_BUCKET)
        self.assertTrue(rule.excludes)

    def test_every_class_of_a_release_has_a_definition_and_a_citation_in_the_report(self) -> None:
        from compare_vep_outputs import REFERENCE_REGISTRIES, divergence_class_definitions

        for release, registry in REFERENCE_REGISTRIES.items():
            classes = divergence_class_definitions(release)
            expected = [r.bucket for r in registry.pair_rules] + [r.bucket for r in registry.one_sided_rules]
            self.assertEqual(list(classes), expected, release)
            for bucket, c in classes.items():
                self.assertTrue(c["definition"].endswith("."), bucket)
                self.assertIn(".pm", c["perl_citation"], bucket)
                self.assertIn(c["excludes"], (True, False))

    def test_the_filter_refuses_the_other_releases_collector(self) -> None:
        from compare_vep_outputs import ReferenceRecordScope

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "discordant.tsv"
            p.write_text("file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n")
            with self.assertRaises(ValueError):
                filter_snp_indel_intended_divergences(p, reference_release="115.2", record_scope=ReferenceRecordScope())
            with self.assertRaises(ValueError):
                filter_snp_indel_intended_divergences(p, reference_release="116.2", scope_features={})
            with self.assertRaises(ValueError):
                filter_snp_indel_intended_divergences(p, reference_release="117.0")

    def test_the_taxonomy_page_names_every_class_of_both_registries(self) -> None:
        """`docs/concordance-provenance/116.2-divergence-classes.md` is the 116.2 registry's
        taxonomy: every bucket the comparator reports under either release is named on it,
        so a class cannot be added here without its definition reaching the page. The
        row-aware pair rules have a section of their own, whose list names them in the
        registry's order and names no class the registry lacks."""
        from compare_vep_outputs import REFERENCE_REGISTRIES

        page = (Path(__file__).resolve().parents[2] / "docs" / "concordance-provenance"
                / "116.2-divergence-classes.md").read_text(encoding="utf-8")
        for registry in REFERENCE_REGISTRIES.values():
            for rule in (*registry.pair_rules, *registry.row_pair_rules, *registry.one_sided_rules):
                self.assertIn(f"`{rule.bucket}`", page, rule.bucket)
        for release in REFERENCE_REGISTRIES:
            self.assertIn(f"`{release}`", page)
        heading = "The pair classes that read the two rows of a pair"
        _, found, after = page.partition(f". {heading}\n")
        self.assertTrue(found, f"the page has no section {heading!r}")
        listed = [line[4:].partition("`")[0] for line in after.partition("\n#")[0].splitlines()
                  if line.startswith("- [`")]
        self.assertEqual(listed, [rule.bucket for registry in REFERENCE_REGISTRIES.values()
                                  for rule in registry.row_pair_rules],
                         "the section lists the registry's row-aware pair rules, in their order")

    def test_the_flag_defaults_to_116_2_and_names_its_choices(self) -> None:
        from compare_vep_outputs import DEFAULT_REFERENCE_RELEASE

        self.assertEqual(DEFAULT_REFERENCE_RELEASE, "116.2")
        self.assertEqual(parse_args([]).reference_release, "116.2")
        self.assertEqual(parse_args(["--reference-release", "115.2"]).reference_release, "115.2")
        with self.assertRaises(SystemExit):
            parse_args(["--reference-release", "116.1"])


class ReferenceRecordScopeTests(unittest.TestCase):
    """The 116.2 collector reads the sorted merge once: a Location whose rows are all
    vep-rs's and structural is a candidate for a record the reference dropped; a Location with any
    reference row, under any allele, is not."""

    @staticmethod
    def key(loc: str, allele: str, feat: str, csq: str, ftype: str = "Transcript") -> str:
        return f"{loc}\t{allele}\t{feat}\t{ftype}\t{csq}"

    def merge(self, perl: list[str], rust: list[str]):
        from compare_vep_outputs import ReferenceRecordScope, compare_sorted_key_files

        scope = ReferenceRecordScope()
        with tempfile.TemporaryDirectory() as d:
            pp, rp = Path(d) / "p.txt", Path(d) / "r.txt"
            pp.write_text("".join(k + "\n" for k in sorted(perl)), encoding="utf-8")
            rp.write_text("".join(k + "\n" for k in sorted(rust)), encoding="utf-8")
            compare_sorted_key_files(pp, rp, file_name="f", discordant_out=None,
                                     track_consequence_buckets=False, record_scope=scope)
        return scope

    def test_a_structural_record_with_no_reference_row_is_absent(self) -> None:
        rust = [self.key("21:800-900", "CPX", "C3", "intron_variant"),
                self.key("21:800-900", "CPX", "C4", "coding_sequence_variant"),
                self.key("21:30000001", "deletion", "TX_DEL", "feature_truncation"),
                self.key("21:1", "A", "TX_SNV", "intron_variant")]
        perl = [self.key("21:1", "A", "TX_SNV", "intron_variant")]
        self.assertEqual(self.merge(perl, rust).absent, {("f", "21:800-900"), ("f", "21:30000001")})

    def test_a_reference_row_under_another_allele_keeps_the_record_present(self) -> None:
        """The same Location carries the reference's `copy_number_variation` rows and vep-rs's
        `deletion` and `duplication` rows: a naming divergence, not a skipped record."""
        perl = [self.key("21:30165801-30166000", "copy_number_variation", "TX_FIRST", "3_prime_UTR_variant")]
        rust = [self.key("21:30165801-30166000", "deletion", "TX_FIRST", "feature_truncation,3_prime_UTR_variant"),
                self.key("21:30165801-30166000", "duplication", "TX_FIRST", "feature_elongation,3_prime_UTR_variant")]
        self.assertEqual(self.merge(perl, rust).absent, set())

    def test_sequence_alleles_and_tandem_repeats_never_make_a_record_absent(self) -> None:
        rust = [self.key("21:100", "A", "TX_FIRST", "missense_variant"),
                self.key("21:900-950", "tandem_repeat", "T7", "feature_elongation")]
        self.assertEqual(self.merge([], rust).absent, set())

    def test_a_mate_position_that_overflows_the_buffer_is_recorded_from_the_reference_side(self) -> None:
        """Ensembl VEP 116.2 writes `N[21:2147478648[` beside an intergenic `N.` and turns the
        buffer's other records intergenic; the position is at `BND_MATE_32_BIT_LIMIT`
        (2^31 - 5000), while 2147478647 and a vep-rs-only bracket at the limit are not the
        shape."""
        from compare_vep_outputs import BND_MATE_32_BIT_LIMIT

        self.assertEqual(BND_MATE_32_BIT_LIMIT, 2147478648)
        perl = [self.key("21:30425651", "N.", "-", "intergenic_variant", "-"),
                self.key("21:30425651", "N[21:2147478648[", "-", "intergenic_variant", "-"),
                self.key("21:40000001", "N[21:2147478647[", "T1", "intron_variant")]
        rust = [self.key("21:30425651", "N.", "T2", "upstream_gene_variant"),
                self.key("21:50000001", "N[21:2147478648[", "T3", "intron_variant")]
        scope = self.merge(perl, rust)
        self.assertEqual(scope.mate_overflow, {("f", "21:30425651")})

    def test_the_last_location_of_the_file_is_closed(self) -> None:
        rust = [self.key("21:999999-1999999", "CPX", "C1", "intron_variant")]
        perl = [self.key("21:1", "A", "TX_SNV", "intron_variant")]
        self.assertEqual(self.merge(perl, rust).absent, {("f", "21:999999-1999999")})

    def test_breakend_records_carry_each_sides_alleles_and_the_mate_chromosomes(self) -> None:
        perl = [self.key("21:26037668", "N[22:26037668[", "M1", "intron_variant"),
                self.key("21:11400001", "N.", "-", "intergenic_variant", "-"),
                self.key("21:11400001", "N[22:14000000[", "-", "intergenic_variant", "-")]
        rust = [self.key("21:26037668", "N[22:26037668[", "M1", "intron_variant"),
                self.key("21:26037668", "chromosome_breakpoint", "O1", "upstream_gene_variant")]
        scope = self.merge(perl, rust)
        self.assertEqual(scope.breakends[("f", "21:26037668")],
                         ({"N[22:26037668["}, {"N[22:26037668[", "chromosome_breakpoint"}, {"22"}))
        self.assertEqual(scope.breakends[("f", "21:11400001")], ({"N.", "N[22:14000000["}, set(), {"22"}))
        self.assertEqual(scope.absent, set(), "every breakend record here has a reference row")


class Apply116OneSidedRulesTests(unittest.TestCase):
    """The 116.2 one-sided pass over discordant rows: the skipped record is set aside, the
    two lost-row breakend shapes and the FASTA shape are counted and left charged, and the
    reference-allele rows of a bare breakpoint, set aside by the key stream before the pass,
    are reported from the scope's count."""

    HEADER = "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set"

    def row(self, source: str, loc: str, allele: str, feature: str, csq: str) -> tuple[str, ...]:
        ftype = "-" if feature == "-" else "Transcript"
        return ("f", source, loc, allele, feature, ftype, _norm(csq))

    def run_pass(self, rows, scope):
        with tempfile.TemporaryDirectory() as d:
            p, o = Path(d) / "discordant.tsv", Path(d) / "open.tsv"
            p.write_text(self.HEADER + "\n" + "".join("\t".join(r) + "\n" for r in rows), encoding="utf-8")
            result = filter_snp_indel_intended_divergences(p, open_out=o, reference_release="116.2", record_scope=scope)
            return result, o.read_text(encoding="utf-8").splitlines()[1:]

    def scope(self, absent=(), breakends=None):
        from compare_vep_outputs import ReferenceRecordScope

        s = ReferenceRecordScope()
        s.absent = {("f", loc) for loc in absent}
        s.breakends = breakends or {}
        return s

    def test_every_vep_rs_row_of_a_skipped_record_is_set_aside(self) -> None:
        rows = [self.row("extra_in_rust", "21:800-900", "CPX", "C3", "intron_variant"),
                self.row("extra_in_rust", "21:800-900", "CPX", "C4", "coding_sequence_variant"),
                self.row("extra_in_rust", "21:30000001", "deletion", "TX_DEL", "feature_truncation"),
                self.row("extra_in_rust", "21:600-700", "<NON_REF>", "-", "intergenic_variant"),
                self.row("missing_in_rust", "21:400", "G", "TX_DOWN", "downstream_gene_variant")]
        (ep, er, buckets), open_rows = self.run_pass(rows, self.scope(absent=["21:800-900", "21:30000001", "21:600-700"]))
        self.assertEqual((ep, er), (0, 4))
        self.assertEqual(buckets["reference_skipped_record"], {"count": 4, "excluded": 4, "records": 3, "perl": 0, "rust": 4})
        self.assertEqual(len(open_rows), 1, "the reference-only row stays")

    def test_a_record_the_reference_kept_is_not_set_aside_by_the_116_2_pass(self) -> None:
        """The 115.2 batch class does not apply: surplus vep-rs transcripts on a record the
        reference annotated stay charged, whatever the span."""
        rows = [self.row("extra_in_rust", "21:18457513-32484691", "N.", "T3", "transcript_ablation")]
        (ep, er, buckets), open_rows = self.run_pass(rows, self.scope())
        self.assertEqual((ep, er, len(open_rows)), (0, 0, 1))
        self.assertEqual(buckets["reference_skipped_record"], {"count": 0, "excluded": 0, "records": 0, "perl": 0, "rust": 0})
        self.assertNotIn("sv_transcript_selection", buckets)

    def test_own_end_rows_the_reference_lost_are_counted_not_set_aside(self) -> None:
        loc = "21:26037668"
        breakends = {("f", loc): ({"N[22:26037668["}, {"N[22:26037668[", "chromosome_breakpoint"}, {"22"})}
        rows = [self.row("extra_in_rust", loc, "chromosome_breakpoint", "O1", "upstream_gene_variant"),
                self.row("extra_in_rust", loc, "chromosome_breakpoint", "O2", "intron_variant")]
        (ep, er, buckets), open_rows = self.run_pass(rows, self.scope(breakends=breakends))
        self.assertEqual((ep, er), (0, 0))
        self.assertEqual(buckets["bnd_own_chromosome_rows_lost"], {"count": 2, "excluded": 0, "records": 1, "perl": 0, "rust": 0})
        self.assertEqual(len(open_rows), 2)

    def test_own_end_rows_beside_a_reference_own_end_row_are_not_the_shape(self) -> None:
        loc = "21:26037668"
        breakends = {("f", loc): ({"N[22:26037668[", "N."}, {"N[22:26037668[", "N."}, {"22"})}
        rows = [self.row("extra_in_rust", loc, "N.", "O2", "intron_variant")]
        (_, _, buckets), _ = self.run_pass(rows, self.scope(breakends=breakends))
        self.assertEqual(buckets["bnd_own_chromosome_rows_lost"]["count"], 0)

    def test_a_same_chromosome_mate_is_not_the_own_end_shape(self) -> None:
        loc = "21:24808769-35740566"
        breakends = {("f", loc): ({"N[21:35740567["}, {"N[21:35740567[", "N."}, {"21"})}
        rows = [self.row("extra_in_rust", loc, "N.", "O2", "intron_variant")]
        (_, _, buckets), _ = self.run_pass(rows, self.scope(breakends=breakends))
        self.assertEqual(buckets["bnd_own_chromosome_rows_lost"]["count"], 0)

    def test_mate_rows_under_a_synonym_named_mate_are_counted(self) -> None:
        loc = "21:26037668"
        syn = "N[NC_000022.11:26037668["
        breakends = {("f", loc): ({"N."}, {"N.", syn}, {"NC_000022.11"})}
        rows = [self.row("extra_in_rust", loc, syn, "M1", "intron_variant")]
        (_, _, buckets), open_rows = self.run_pass(rows, self.scope(breakends=breakends))
        self.assertEqual(buckets["bnd_synonym_mate_rows_lost"], {"count": 1, "excluded": 0, "records": 1, "perl": 0, "rust": 0})
        self.assertEqual(len(open_rows), 1)

    def test_a_canonical_mate_name_is_not_the_synonym_shape(self) -> None:
        loc = "21:26037668"
        breakends = {("f", loc): ({"N."}, {"N.", "N[chr22:26037668["}, {"22"})}
        rows = [self.row("extra_in_rust", loc, "N[chr22:26037668[", "M1", "intron_variant")]
        (_, _, buckets), _ = self.run_pass(rows, self.scope(breakends=breakends))
        self.assertEqual(buckets["bnd_synonym_mate_rows_lost"]["count"], 0)

    def test_a_mitochondrial_mate_written_m_is_a_synonym_of_the_caches_mt(self) -> None:
        """The slice test strips `chr` and nothing else, so `M` is not the cache's `MT`."""
        loc = "21:26037668"
        breakends = {("f", loc): ({"N."}, {"N.", "N[chrM:100["}, {"M"})}
        rows = [self.row("extra_in_rust", loc, "N[chrM:100[", "MT_TX", "upstream_gene_variant")]
        (_, _, buckets), _ = self.run_pass(rows, self.scope(breakends=breakends))
        self.assertEqual(buckets["bnd_synonym_mate_rows_lost"]["count"], 1)
        breakends = {("f", loc): ({"N."}, {"N.", "N[MT:100["}, {"MT"})}
        rows = [self.row("extra_in_rust", loc, "N[MT:100[", "MT_TX", "upstream_gene_variant")]
        (_, _, buckets), _ = self.run_pass(rows, self.scope(breakends=breakends))
        self.assertEqual(buckets["bnd_synonym_mate_rows_lost"]["count"], 0)

    def test_the_reference_allele_rows_of_a_bare_breakpoint_are_set_aside_and_counted(self) -> None:
        """Ensembl VEP 116.2 writes a bare `<BND>`'s rows twice, ALLELE_NUM 1 for the reference
        allele N of `N/<BND>` and 2 for the ALT (StructuralVariationOverlap.pm:64-77); vep-rs
        writes each once. The key stream carries each tuple once per side, so the merge counts
        3 reference tuples, 3 vep-rs tuples and 3 in the intersection (the reference-allele rows
        enter neither figure), the scope counts the 2 rows the stream dropped, and the one-sided
        pass reports them under the class, all excluded, with nothing subtracted from the
        reference's tuple count."""
        from compare_vep_outputs import (
            ReferenceRecordScope,
            _sort_unique_stream_to_file,
            compare_sorted_key_files,
            count_bare_breakpoint_reference_rows,
            iter_vep_key_lines,
        )

        header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
        own = "bare_bnd\t21:30425651\tchromosome_breakpoint\tG\tENST00000390690\tTranscript\tcoding_sequence_variant\n"
        up = "bare_bnd\t21:30425651\tchromosome_breakpoint\tG\tENST00000334068\tTranscript\tupstream_gene_variant\n"
        mate = "paired\t21:30425651\tN[21:30425700[\tG\tENST00000390690\tTranscript\tcoding_sequence_variant\n"
        scope = ReferenceRecordScope()
        with tempfile.TemporaryDirectory() as d:
            perl, rust = Path(d) / "perl.txt", Path(d) / "rust.txt"
            perl.write_text(header + own + own + up + up + mate + mate, encoding="utf-8")
            rust.write_text(header + own + up + mate, encoding="utf-8")
            ps, rs = Path(d) / "p.keys", Path(d) / "r.keys"
            _sort_unique_stream_to_file(iter_vep_key_lines(perl), ps)
            _sort_unique_stream_to_file(iter_vep_key_lines(rust), rs)
            perl_count, rust_count, intersection, missing, extra, _, _ = compare_sorted_key_files(
                ps, rs, file_name="f", discordant_out=None, track_consequence_buckets=False, record_scope=scope
            )
            scope.record_bare_reference_rows("f", count_bare_breakpoint_reference_rows(perl))
        self.assertEqual((perl_count, rust_count, intersection, missing, extra), (3, 3, 3, 0, 0))
        self.assertEqual(scope.bare_reference_rows, {("f", "21:30425651"): 2},
                         "the bracket record's duplicate mate rows are not the shape")
        (ep, er, buckets), open_rows = self.run_pass([], scope)
        self.assertEqual((ep, er, open_rows), (0, 0, []))
        self.assertEqual(buckets["bnd_bare_reference_allele_row"], {"count": 2, "excluded": 2, "records": 1, "perl": 0, "rust": 0})

    def test_the_reference_rows_of_a_record_whose_mate_overflows_the_buffer_are_counted(self) -> None:
        """The reference's rows at the breakend record's Location are counted and stay
        charged; vep-rs's rows there and the reference's rows at the buffer's other records
        are not the shape, since the outputs do not tell which records shared the buffer."""
        rows = [self.row("missing_in_rust", "21:30425651", "N.", "-", "intergenic_variant"),
                self.row("missing_in_rust", "21:30425651", "N[21:2147478648[", "-", "intergenic_variant"),
                self.row("extra_in_rust", "21:30425651", "N.", "T2", "feature_truncation,coding_sequence_variant"),
                self.row("missing_in_rust", "21:30425651-30425899", "duplication", "-", "intergenic_variant")]
        scope = self.scope()
        scope.mate_overflow = {("f", "21:30425651")}
        (ep, er, buckets), open_rows = self.run_pass(rows, scope)
        self.assertEqual((ep, er, len(open_rows)), (0, 0, 4))
        self.assertEqual(buckets["bnd_mate_32bit_buffer_intergenic"], {"count": 2, "excluded": 0, "records": 1, "perl": 0, "rust": 0})

    def test_a_gvcf_variant_site_the_reference_dropped_is_set_aside_from_its_own_name(self) -> None:
        """Ensembl VEP 116.2 drops `21 30425650 . N G,<NON_REF>` whole (`G/NON_REF is not a
        supported structural variant type`); vep-rs annotates the `G` allele, so its rows
        carry a sequence allele at a Location the reference wrote nothing for. The record's
        generated name, `21_30425650_N/G/<NON_REF>`, read from vep-rs's output, admits the
        Location, and every row there is set aside under the class; a record named by an ID
        gives no such signal here and is left to the structural-variant comparator."""
        from compare_vep_outputs import (
            ReferenceRecordScope,
            _sort_unique_stream_to_file,
            compare_sorted_key_files,
            gvcf_site_locations,
            iter_vep_key_lines,
        )

        header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
        site = ("21_30425650_N/G/<NON_REF>\t21:30425650\tG\tG\tENST00000390690\tTranscript\tcoding_sequence_variant\n"
                "21_30425650_N/G/<NON_REF>\t21:30425650\tG\tG\tENST00000334068\tTranscript\tupstream_gene_variant\n")
        named = "rs1\t21:40000000\tT\tG\tENST_X\tTranscript\tmissense_variant\n"
        block = "21_50000000_N/<NON_REF>\t21:50000001-50000100\t<NON_REF>\t-\t-\t-\tintergenic_variant\n"
        scope = ReferenceRecordScope()
        with tempfile.TemporaryDirectory() as d:
            perl, rust = Path(d) / "perl.txt", Path(d) / "rust.txt"
            perl.write_text(header, encoding="utf-8")
            rust.write_text(header + site + named + block, encoding="utf-8")
            self.assertEqual(gvcf_site_locations(rust), {"21:30425650"})
            self.assertEqual(gvcf_site_locations(Path(d) / "absent.txt"), set())
            ps, rs = Path(d) / "p.keys", Path(d) / "r.keys"
            _sort_unique_stream_to_file(iter_vep_key_lines(perl), ps)
            _sort_unique_stream_to_file(iter_vep_key_lines(rust), rs)
            scope.record_gvcf_sites("f", gvcf_site_locations(rust))
            compare_sorted_key_files(ps, rs, file_name="f", discordant_out=None,
                                     track_consequence_buckets=False, record_scope=scope)
        self.assertEqual(scope.gvcf_absent, {("f", "21:30425650")})
        self.assertEqual(scope.absent, {("f", "21:50000001-50000100")}, "the lone block is the skipped record")
        rows = [self.row("extra_in_rust", "21:30425650", "G", "ENST00000390690", "coding_sequence_variant"),
                self.row("extra_in_rust", "21:30425650", "G", "ENST00000334068", "upstream_gene_variant"),
                self.row("extra_in_rust", "21:40000000", "G", "ENST_X", "missense_variant"),
                self.row("extra_in_rust", "21:50000001-50000100", "<NON_REF>", "-", "intergenic_variant")]
        (ep, er, buckets), open_rows = self.run_pass(rows, scope)
        self.assertEqual((ep, er), (0, 3))
        self.assertEqual(buckets["gvcf_variant_site_dropped"], {"count": 2, "excluded": 2, "records": 1, "perl": 0, "rust": 2})
        self.assertEqual(buckets["reference_skipped_record"]["excluded"], 1)
        self.assertEqual([r.split("\t")[2] for r in open_rows], ["21:40000000"])

    def test_a_gvcf_variant_site_is_recognised_by_its_shape_from_the_input(self) -> None:
        """A record named by an ID shows no `<NON_REF>` in vep-rs's output, so the input VCF
        names the site by its shape: `<NON_REF>` beside a real allele, whatever the ID. The
        Locations vep-rs writes for the site follow its parser's sequence path
        (`gvcf_site_rows`): `ATG AG` loses its anchor base and is minimised from both ends to
        the deletion of `T` at POS + 1, so the ID'd row `rs1 21:30425651 -` is set aside under
        the class once the input is read, and stays charged without it. A list of blocks
        alone, or `*` beside a block, is no site."""
        from compare_vep_outputs import (
            ReferenceRecordScope,
            _sort_unique_stream_to_file,
            compare_sorted_key_files,
            gvcf_site_locations,
            gvcf_site_locations_from_vcf,
            input_vcf_for,
            is_gvcf_variant_site,
            iter_vep_key_lines,
        )

        self.assertTrue(is_gvcf_variant_site(["G", "<NON_REF>"]))
        self.assertTrue(is_gvcf_variant_site(["<non_ref>", "AT"]))
        self.assertFalse(is_gvcf_variant_site(["<NON_REF>"]))
        self.assertFalse(is_gvcf_variant_site(["<NON_REF>", "<NON_REF>"]))
        self.assertFalse(is_gvcf_variant_site(["*", "<NON_REF>"]))
        self.assertFalse(is_gvcf_variant_site(["G", "T"]))
        # A block beside a symbolic allele is an unsupported structural list (the
        # skipped-record class); under a true SVTYPE the block is an allele of a structural
        # record in both engines (Ensembl VEP 116.2 writes two `deletion` rows per transcript
        # for `G,<NON_REF>` with `SVTYPE=DEL`, as vep-rs does).
        self.assertFalse(is_gvcf_variant_site(["<DEL>", "<NON_REF>"]))
        self.assertFalse(is_gvcf_variant_site(["G", "<DEL>", "<NON_REF>"]))
        self.assertFalse(is_gvcf_variant_site(["G", "<NON_REF>"], "DEL"))
        self.assertTrue(is_gvcf_variant_site(["G", "<NON_REF>"], None))
        self.assertTrue(is_gvcf_variant_site(["G", "<NON_REF>"], "0"))
        header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
        named = "rs1\t21:30425651\t-\tG\tENST00000390690\tTranscript\tframeshift_variant\n"
        blocks = "21_50000000_N/<NON_REF>/<NON_REF>\t21:50000001-50000100\t<NON_REF>\t-\t-\t-\tintergenic_variant\n"
        vcf = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
               "21\t30425650\trs1\tATG\tAG,<NON_REF>\t.\t.\t.\n"
               "21\t50000000\t.\tN\t<NON_REF>,<NON_REF>\t.\t.\tEND=50000100\n"
               "21\t60000000\t.\tN\t<NON_REF>\t.\t.\tEND=60000100\n"
               "chr22\t100\t.\tA\tG,<NON_REF>\t.\t.\t.\n")
        with tempfile.TemporaryDirectory() as d:
            perl, rust = Path(d) / "s.txt", Path(d) / "rust" / "s.txt"
            rust.parent.mkdir()
            perl.write_text(header, encoding="utf-8")
            rust.write_text(header + named + blocks, encoding="utf-8")
            inputs = Path(d) / "inputs"
            inputs.mkdir()
            (inputs / "s.vcf").write_text(vcf, encoding="utf-8")
            import gzip

            with gzip.open(inputs / "z.vcf.gz", "wt", encoding="utf-8") as handle:
                handle.write(vcf)
            self.assertEqual(gvcf_site_locations(rust), set(), "a block list and an ID'd row name no site")
            for path in (inputs / "s.vcf", inputs / "z.vcf.gz"):
                self.assertEqual(gvcf_site_locations_from_vcf(path), {"21:30425651", "chr22:100"}, path.name)
            self.assertEqual(gvcf_site_locations_from_vcf(inputs / "s.vcf", frozenset({"21"})), {"21:30425651"})
            self.assertEqual(input_vcf_for(perl, inputs, None), inputs / "s.vcf")
            self.assertEqual(input_vcf_for(Path(d) / "z.txt", inputs, None), inputs / "z.vcf.gz")
            self.assertIsNone(input_vcf_for(Path(d) / "none.txt", inputs, None))
            self.assertEqual(input_vcf_for(perl, inputs, inputs / "z.vcf.gz"), inputs / "z.vcf.gz")
            self.assertIsNone(input_vcf_for(perl, None, None))
            ps, rs = Path(d) / "p.keys", Path(d) / "r.keys"
            _sort_unique_stream_to_file(iter_vep_key_lines(perl), ps)
            _sort_unique_stream_to_file(iter_vep_key_lines(rust), rs)
            for with_input in (False, True):
                scope = ReferenceRecordScope()
                sites = gvcf_site_locations(rust)
                if with_input:
                    sites |= gvcf_site_locations_from_vcf(inputs / "s.vcf")
                scope.record_gvcf_sites("f", sites)
                compare_sorted_key_files(ps, rs, file_name="f", discordant_out=None,
                                         track_consequence_buckets=False, record_scope=scope)
                rows = [self.row("extra_in_rust", "21:30425651", "-", "ENST00000390690", "frameshift_variant"),
                        self.row("extra_in_rust", "21:50000001-50000100", "<NON_REF>", "-", "intergenic_variant")]
                (ep, er, buckets), open_rows = self.run_pass(rows, scope)
                if with_input:
                    self.assertEqual(scope.gvcf_absent, {("f", "21:30425651")})
                    self.assertEqual((ep, er, open_rows), (0, 2, []))
                    self.assertEqual(buckets["gvcf_variant_site_dropped"]["rust"], 1)
                else:
                    self.assertEqual(scope.gvcf_absent, set())
                    self.assertEqual([r.split("\t")[2] for r in open_rows], ["21:30425651"])
                self.assertEqual(buckets["reference_skipped_record"]["excluded"], 1, "the block list is a skipped record")
            args = parse_args(["--perl-dir", "p", "--rust-dir", "r", "--report-dir", "o", "--input-dir", str(inputs)])
            self.assertEqual((args.input_dir, args.input_vcf), (str(inputs), None))
            args = parse_args(["--perl-dir", "p", "--rust-dir", "r", "--report-dir", "o", "--input-vcf", str(inputs / "s.vcf")])
            self.assertEqual((args.input_dir, args.input_vcf), (None, str(inputs / "s.vcf")))

    def test_a_gvcf_site_is_read_at_the_pos_vep_rs_numifies(self) -> None:
        """vep-rs reads a POS in any spelling of an exact integer that Rust's `f64` parser
        accepts (`3.042565e7`, `30425650.0`, `+30425650`, `30425650.`, ` 30425650 `), so the gVCF
        site `ATG AG,<NON_REF>` at each is the deletion at `21:30425651`; it fails the line for
        `30_425_650`, `-30425650`, `abc`, fullwidth digits, `0x1D04A1A`, `inf` and `1e400`
        (`Invalid POS`), so those records name no site. Only the unpadded digit string is read
        exactly (`parse_coordinate`); a padded one goes through the double, so ` 72057594037927937 `
        is 72057594037927936 and the ceiling with blanks around it is refused."""
        from compare_vep_outputs import VEP_RS_COORDINATE_CEILING, gvcf_site_locations_from_vcf, vep_rs_pos

        for text in ("30425650", "3.042565e7", "30425650.0", "+30425650", "30425650.", " 30425650 ", "030425650"):
            self.assertEqual(vep_rs_pos(text), 30425650, repr(text))
        self.assertEqual(vep_rs_pos("3e7"), 30000000)
        self.assertEqual(VEP_RS_COORDINATE_CEILING, 2**62 - 1)
        self.assertEqual(vep_rs_pos(str(VEP_RS_COORDINATE_CEILING)), VEP_RS_COORDINATE_CEILING)
        self.assertEqual(vep_rs_pos(str(2**56 + 1)), 2**56 + 1)
        self.assertEqual(vep_rs_pos(f" {2**56 + 1} "), 2**56)
        for text in (".", "", "abc", "-30425650", "30_425_650", "３０４２５６５０", "0x1D04A1A", "inf", "nan", "1e400",
                     "30425650.5", "++30425650", str(VEP_RS_COORDINATE_CEILING + 1), f" {VEP_RS_COORDINATE_CEILING} "):
            self.assertIsNone(vep_rs_pos(text), repr(text))
        read = ("3.042565e7", "30425650.0", "+30425650", "30425650.", " 30425650 ")
        refused = ("30_425_650", "-30425650", "abc", "３０４２５６５０", "0x1D04A1A", "inf", "1e400")
        header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

        def site(pos: str, vid: str) -> str:
            return f"21\t{pos}\t{vid}\tATG\tAG,<NON_REF>\t.\t.\t.\n"

        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "pos.vcf"
            path.write_text(header + "".join(site(pos, f"p{i}") for i, pos in enumerate(read + refused)), encoding="utf-8")
            self.assertEqual(gvcf_site_locations_from_vcf(path), {"21:30425651"})
            path.write_text(header + "".join(site(pos, f"r{i}") for i, pos in enumerate(refused)), encoding="utf-8")
            self.assertEqual(gvcf_site_locations_from_vcf(path), set())

    def test_a_gvcf_site_allele_is_written_in_upper_case(self) -> None:
        """vep-rs writes `21 30425650 x A g,<NON_REF> . . .` and `21 30425650 x a g,<non_ref>
        . . .` as `21:30425650 G`; the derived allele is upper-cased as the parser writes it."""
        from compare_vep_outputs import gvcf_site_rows

        self.assertEqual(gvcf_site_rows(30425650, "A", ["g", "<NON_REF>"]), [("G", 30425650, 30425650)])
        self.assertEqual(gvcf_site_rows(30425650, "a", ["g", "<non_ref>"]), [("G", 30425650, 30425650)])
        self.assertEqual(gvcf_site_rows(30425650, "atg", ["ag", "<NON_REF>"]), [("-", 30425651, 30425651)])

    def test_a_star_beside_a_gvcf_sites_indel_keeps_the_anchor_trimmed_alleles(self) -> None:
        """vep-rs minimises the one allele beside the block only when it is the one ALT there: for
        `21 30425650 x ATT AT,*,<NON_REF> . . .` it writes `21:30425651-30425652 T`, for
        `AT ATT,*,<NON_REF>` `21:30425651 TT` and for `ATTG ATC,*,<NON_REF>` `21:30425651-30425653
        TC`, the anchor-trimmed alleles Ensembl VEP 116.2 writes for the same records without the
        block (`minimise_alleles` leaves the three-part string as it is, Parser.pm:978-984); a `.`
        beside the indel is compared and never shares REF's first base, so `ATT AT,.,<NON_REF>`
        keeps its anchor (`21:30425650-30425652 AT`, as Ensembl VEP 116.2 writes `ATT AT,.`); and
        without either the allele is minimised from both ends (`21:30425652 -`,
        `21:30425651-30425652 T`, `21:30425652-30425653 C`)."""
        from compare_vep_outputs import gvcf_site_rows

        self.assertEqual(gvcf_site_rows(30425650, "ATT", ["AT", "*", "<NON_REF>"]), [("T", 30425651, 30425652)])
        self.assertEqual(gvcf_site_rows(30425650, "AT", ["ATT", "*", "<NON_REF>"]), [("TT", 30425651, 30425651)])
        self.assertEqual(gvcf_site_rows(30425650, "ATTG", ["ATC", "*", "<NON_REF>"]), [("TC", 30425651, 30425653)])
        self.assertEqual(gvcf_site_rows(30425650, "ATT", ["AT", ".", "<NON_REF>"]), [("AT", 30425650, 30425652)])
        self.assertEqual(gvcf_site_rows(30425650, "ATT", ["AT", "<NON_REF>"]), [("-", 30425652, 30425652)])
        self.assertEqual(gvcf_site_rows(30425650, "AT", ["ATT", "<NON_REF>"]), [("T", 30425652, 30425651)])
        self.assertEqual(gvcf_site_rows(30425650, "ATTG", ["ATC", "<NON_REF>"]), [("C", 30425652, 30425653)])

    def test_a_star_counts_as_an_allele_of_one_base_in_a_gvcf_sites_anchor_trim(self) -> None:
        """The anchor trim counts a `*` as an allele of one base and leaves it out of the first-base
        comparison (Parser/VCF.pm:296-306), so vep-rs writes `21 30425650 x AG AC,*,<NON_REF> . . .`
        at `21:30425651 C`, `AGC ATT,*,<NON_REF>` at `21:30425651-30425652 TT`, `AG *,AC,<NON_REF>` at
        `21:30425651 C` and `AG AC,AT,*,<NON_REF>` at `21:30425651 C` and `T`, the rows Ensembl VEP
        116.2 writes for the records without the block; the trim needs the shared base (`AG TC,*`
        stays at `21:30425650-30425651 TC`), a one-base REF has nothing to trim (`A G,*`, `A *,G`,
        `A G,T,*` stay at 30425650), without the `*` the same-length record is no indel
        (`AG AC,AT` stays `21:30425650-30425651`), a `.` is compared and never shares the base
        (`AG AC,.` and `AGC ATT,.` keep their anchors), and a `<*>` beside an indel allele is left
        out of the comparison and kept whole (`A AG,<*>` is `G` and `<*>` at `21:30425650-30425651`,
        as the binary writes them)."""
        from compare_vep_outputs import gvcf_site_rows, sequence_location

        def rows(ref: str, alts: str) -> list[tuple[str, str]]:
            return sorted((sequence_location("21", start, end), allele)
                          for allele, start, end in gvcf_site_rows(30425650, ref, alts.split(",")))

        self.assertEqual(rows("A", "G,*,<NON_REF>"), [("21:30425650", "G")])
        self.assertEqual(rows("AG", "AC,*,<NON_REF>"), [("21:30425651", "C")])
        self.assertEqual(rows("AG", "TC,*,<NON_REF>"), [("21:30425650-30425651", "TC")])
        self.assertEqual(rows("AGC", "ATT,*,<NON_REF>"), [("21:30425651-30425652", "TT")])
        self.assertEqual(rows("AG", "A,*,<NON_REF>"), [("21:30425651", "-")])
        self.assertEqual(rows("A", "AG,*,<NON_REF>"), [("21:30425650-30425651", "G")])
        self.assertEqual(rows("AG", "*,AC,<NON_REF>"), [("21:30425651", "C")])
        self.assertEqual(rows("AG", "AC,AT,*,<NON_REF>"), [("21:30425651", "C"), ("21:30425651", "T")])
        self.assertEqual(rows("A", "*,G,<NON_REF>"), [("21:30425650", "G")])
        self.assertEqual(rows("AG", "*,A,<NON_REF>"), [("21:30425651", "-")])
        self.assertEqual(rows("A", "G,T,*,<NON_REF>"), [("21:30425650", "G"), ("21:30425650", "T")])
        self.assertEqual(rows("AG", "AC,AT,<NON_REF>"), [("21:30425650-30425651", "AC"), ("21:30425650-30425651", "AT")])
        self.assertEqual(rows("AG", "AC,.,<NON_REF>"), [("21:30425650-30425651", "AC")])
        self.assertEqual(rows("AGC", "ATT,.,<NON_REF>"), [("21:30425650-30425652", "ATT")])
        self.assertEqual(rows("A", "AG,<*>,<NON_REF>"), [("21:30425650-30425651", "<*>"), ("21:30425650-30425651", "G")])

    def test_gvcf_site_rows_follow_the_parsers_sequence_path(self) -> None:
        """vep-rs trims an indel record's anchor base across every real allele, then minimises
        a lone allele from both ends (`-` for an emptied allele); an insertion ends one base
        before its start; same-length alleles keep POS and their bases; a multi-allelic site
        keeps its anchor-trimmed alleles."""
        from compare_vep_outputs import gvcf_site_rows, sequence_location

        self.assertEqual(gvcf_site_rows(100, "A", ["G", "<NON_REF>"]), [("G", 100, 100)])
        self.assertEqual(gvcf_site_rows(100, "AT", ["GC", "<NON_REF>"]), [("GC", 100, 101)])
        self.assertEqual(gvcf_site_rows(100, "A", ["AT", "<NON_REF>"]), [("T", 101, 100)])
        self.assertEqual(gvcf_site_rows(100, "N", ["NA", "<NON_REF>"]), [("A", 101, 100)])
        self.assertEqual(gvcf_site_rows(100, "ATG", ["A", "<NON_REF>"]), [("-", 101, 102)])
        self.assertEqual(gvcf_site_rows(100, "ATG", ["AG", "<NON_REF>"]), [("-", 101, 101)])
        self.assertEqual(gvcf_site_rows(100, "A", ["AT", "ATT", "<NON_REF>"]), [("T", 101, 100), ("TT", 101, 100)])
        self.assertEqual(gvcf_site_rows(100, "A", ["G", "T", "<NON_REF>"]), [("G", 100, 100), ("T", 100, 100)])
        self.assertEqual(gvcf_site_rows(100, "A", ["<NON_REF>", "<NON_REF>"]), [])
        self.assertEqual(gvcf_site_rows(100, "A", ["*", "G", "<NON_REF>"]), [("G", 100, 100)])
        self.assertEqual(sequence_location("21", 101, 100), "21:100-101")
        self.assertEqual(sequence_location("21", 101, 102), "21:101-102")
        self.assertEqual(sequence_location("chr21", 100, 100), "chr21:100")

    def test_a_trailing_empty_alt_field_is_no_field(self) -> None:
        """Ensembl VEP's parser splits the ALT column in list context (BaseVCF4.pm:489
        `[split(',',$alt_allele)]`), which drops trailing empty fields and keeps the middle
        ones: `split(',', "<DEL>,")` has one field and `split(',', "G,,T,,")` is `G`, the empty
        field, `T`. So `A G,` is a plain SNV, no site; `A G,<NON_REF>,` is the site `A G,<NON_REF>`
        with the one row `G` at POS; `<NON_REF>,` is a block alone. The generated name of an
        ID-less record carries the ALT fields after the REF, so a trailing `/` is the same empty
        field."""
        from compare_vep_outputs import (
            gvcf_site_locations,
            gvcf_site_locations_from_vcf,
            gvcf_site_rows,
            is_gvcf_variant_site,
            perl_fields,
            perl_split,
        )

        self.assertEqual(perl_split("<DEL>,", ","), ["<DEL>"])
        self.assertEqual(perl_split("G,,T,,", ","), ["G", "", "T"])
        self.assertEqual(perl_split("", ","), [])
        self.assertEqual(perl_split(",", ","), [])
        self.assertEqual(perl_split("N/G/<NON_REF>/", "/"), ["N", "G", "<NON_REF>"])
        self.assertEqual(perl_fields(["G", "<NON_REF>", "", ""]), ["G", "<NON_REF>"])
        self.assertFalse(is_gvcf_variant_site(["G", ""]))
        self.assertTrue(is_gvcf_variant_site(["G", "<NON_REF>", ""]))
        self.assertFalse(is_gvcf_variant_site(["<NON_REF>", ""]), "a block alone is no site")
        self.assertEqual(gvcf_site_rows(100, "A", ["G", ""]), [("G", 100, 100)])
        self.assertEqual(gvcf_site_rows(100, "A", ["G", "<NON_REF>", ""]), [("G", 100, 100)])
        header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
        rows = ("21_30425650_A/G/<NON_REF>/\t21:30425650\tG\tG\tENST_1\tTranscript\tmissense_variant\n"
                "21_40000000_A/G/\t21:40000000\tG\tG\tENST_2\tTranscript\tmissense_variant\n")
        vcf = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
               "21\t30425650\trs1\tA\tG,<NON_REF>,\t.\t.\t.\n"
               "21\t40000000\trs2\tA\tG,\t.\t.\t.\n"
               "21\t50000000\trs3\tA\t<NON_REF>,\t.\t.\tEND=50000100\n")
        with tempfile.TemporaryDirectory() as d:
            rust, inputs = Path(d) / "s.txt", Path(d) / "s.vcf"
            rust.write_text(header + rows, encoding="utf-8")
            inputs.write_text(vcf, encoding="utf-8")
            self.assertEqual(gvcf_site_locations(rust), {"21:30425650"})
            self.assertEqual(gvcf_site_locations_from_vcf(inputs), {"21:30425650"})

    def test_an_info_value_is_read_to_a_second_equals_sign(self) -> None:
        """`my ($key,$value) = split('=',$info)` (BaseVCF4.pm:582) is a two-element list
        assignment, so a value ends at a second `=`: both engines write the SNV rows of
        `21 30425650 x N G . . SVTYPE==DEL;END=30425949` at `21:30425650` (the SVTYPE is the empty
        value, false, and the record takes the sequence path) and the mate of
        `21 30425650 x N <BND> . . SVTYPE=BND;CHR2=21;END2=30425700=5` as `N[21:30425700[`. A
        `G,<NON_REF>` site under `SVTYPE==DEL` is a gVCF variant site."""
        from compare_vep_outputs import (
            _info_written,
            gvcf_site_locations_from_vcf,
            is_gvcf_variant_site,
        )

        self.assertEqual(_info_written("SVTYPE==DEL;END=30425949", "SVTYPE"), "")
        self.assertEqual(_info_written("SVTYPE==DEL;END=30425949", "END"), "30425949")
        self.assertEqual(_info_written("SVTYPE=BND;CHR2=21;END2=30425700=5", "END2"), "30425700")
        self.assertEqual(_info_written("END2=30425700=5=7", "END2"), "30425700")
        self.assertIsNone(_info_written("SVTYPE==DEL;END", "END"))
        self.assertTrue(is_gvcf_variant_site(["G", "<NON_REF>"], _info_written("SVTYPE==DEL", "SVTYPE")))
        with tempfile.TemporaryDirectory() as d:
            vcf = Path(d) / "s.vcf"
            vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                           "21\t30425650\trs1\tA\tG,<NON_REF>\t.\t.\tSVTYPE==DEL;END=30425949\n", encoding="utf-8")
            self.assertEqual(gvcf_site_locations_from_vcf(vcf), {"21:30425650"})

    def test_a_dot_breakend_pair_and_a_non_integer_end_are_set_aside(self) -> None:
        """Ensembl VEP 116.2 writes `21:30425650-30425651 chromosome_breakpoint` for `N N.`
        without SVTYPE where vep-rs writes the point `21:30425651`: the two rows on one
        feature are a pair, both set aside; the reference's row on a feature vep-rs did not
        name stays charged. A reference row whose Location carries `END=30425949.5` as
        written is set aside on its own."""
        rows = [self.row("missing_in_rust", "21:30425650-30425651", "chromosome_breakpoint", "T2", "coding_sequence_variant"),
                self.row("extra_in_rust", "21:30425651", "chromosome_breakpoint", "T2", "coding_sequence_variant"),
                self.row("missing_in_rust", "21:30425650-30425651", "chromosome_breakpoint", "T1", "upstream_gene_variant"),
                self.row("missing_in_rust", "21:30425651-30425949.5", "deletion", "T2", "feature_truncation,coding_sequence_variant"),
                self.row("missing_in_rust", "21:30425651-3.0425949e7", "deletion", "T1", "upstream_gene_variant"),
                self.row("extra_in_rust", "21:30425651", "chromosome_breakpoint", "T3", "intron_variant")]
        (ep, er, buckets), open_rows = self.run_pass(rows, self.scope())
        self.assertEqual((ep, er), (3, 1))
        self.assertEqual(buckets["dot_breakend_without_svtype_pair"], {"count": 2, "excluded": 2, "records": 2, "perl": 1, "rust": 1})
        self.assertEqual(buckets["non_integer_end_annotated"], {"count": 2, "excluded": 2, "records": 2, "perl": 2, "rust": 0})
        self.assertEqual(sorted(r.split("\t")[4] for r in open_rows), ["T1", "T3"])
        # Two vep-rs point rows under one feature (two consequence sets) count the
        # reference's pair row once.
        rows = [self.row("missing_in_rust", "21:30425650-30425651", "chromosome_breakpoint", "T2", "coding_sequence_variant"),
                self.row("extra_in_rust", "21:30425651", "chromosome_breakpoint", "T2", "coding_sequence_variant"),
                self.row("extra_in_rust", "21:30425651", "chromosome_breakpoint", "T2", "coding_sequence_variant,feature_truncation")]
        (ep, er, buckets), open_rows = self.run_pass(rows, self.scope())
        self.assertEqual((ep, er, open_rows), (1, 2, []))
        self.assertEqual(buckets["dot_breakend_without_svtype_pair"], {"count": 3, "excluded": 3, "records": 2, "perl": 1, "rust": 2})

    def test_a_location_is_keyed_by_its_numbers(self) -> None:
        """Ensembl VEP prints an INFO/END as written (`21:30425651-030425949`,
        `21:30425651-+30425949`); the key stream carries the number, so the rows agree with
        vep-rs's `21:30425651-30425949`. A fractional coordinate stays as written."""
        from compare_vep_outputs import canonical_location, iter_vep_key_lines

        for text, canonical in [("21:30425651-030425949", "21:30425651-30425949"),
                                ("21:30425651-+30425949", "21:30425651-30425949"),
                                ("21:+30425650", "21:30425650"),
                                ("21:30425651-30425949.5", "21:30425651-30425949.5"),
                                ("21:30425651", "21:30425651"), ("MT", "MT")]:
            self.assertEqual(canonical_location(text), canonical, text)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "perl.txt"
            path.write_text("x\t21:30425651-030425949\tdeletion\tG\tT1\tTranscript\tupstream_gene_variant\n", encoding="utf-8")
            self.assertEqual(list(iter_vep_key_lines(path)), ["21:30425651-30425949\tdeletion\tT1\tTranscript\tupstream_gene_variant"])

    def test_a_reference_intergenic_row_against_vep_rs_transcripts_is_counted(self) -> None:
        rows = [self.row("missing_in_rust", "21:100", "A", "-", "intergenic_variant"),
                self.row("extra_in_rust", "21:100", "A", "TX_FIRST", "missense_variant"),
                self.row("extra_in_rust", "21:100", "A", "TX_SECOND", "intron_variant"),
                self.row("missing_in_rust", "21:200", "T", "-", "intergenic_variant"),
                self.row("extra_in_rust", "21:200", "T", "-", "intergenic_variant")]
        (ep, er, buckets), open_rows = self.run_pass(rows, self.scope())
        self.assertEqual((ep, er), (0, 0))
        self.assertEqual(buckets["fasta_named_slice_intergenic"], {"count": 1, "excluded": 0, "records": 1, "perl": 0, "rust": 0})
        self.assertEqual(len(open_rows), 5)


    def reference_only_scope(self, *locations: str):
        s = self.scope()
        s.reference_only = {("f", loc) for loc in locations}
        return s

    def test_a_negative_first_coordinate_is_the_numified_pos_class(self) -> None:
        """Ensembl VEP 116.2 writes `21:-4-295 deletion - intergenic_variant` for `21 -5 x N <DEL>
        . . SVTYPE=DEL;SVLEN=-300`, `21:-2147483647--2147483348 deletion` for POS -2147483648
        and `21:-5 G - intergenic_variant` for `21 -5 x A G`; vep-rs refuses each line, so the
        rows are set aside. A reference-only row at a start of 1 (`21:1-30425949`, POS `.`
        under --max_sv_size -1) is not told from a record at POS 0 and stays charged here."""
        rows = [self.row("missing_in_rust", "21:-4-295", "deletion", "-", "intergenic_variant"),
                self.row("missing_in_rust", "21:-2147483647--2147483348", "deletion", "-", "intergenic_variant"),
                self.row("missing_in_rust", "21:-5", "G", "-", "intergenic_variant"),
                self.row("missing_in_rust", "21:1-30425949", "deletion", "ENST00000262354", "transcript_ablation")]
        scope = self.reference_only_scope("21:-4-295", "21:-2147483647--2147483348", "21:-5", "21:1-30425949")
        (ep, er, buckets), open_rows = self.run_pass(rows, scope)
        self.assertEqual((ep, er), (3, 0))
        self.assertEqual(buckets["non_numeric_pos_numified"], {"count": 3, "excluded": 3, "records": 3, "perl": 3, "rust": 0})
        self.assertEqual(buckets["non_integer_end_annotated"]["count"], 0, "a negative coordinate is an integer")
        self.assertEqual([r.split("\t")[2] for r in open_rows], ["21:1-30425949"])

    def test_the_dot_breakend_pair_covers_the_own_and_mate_alleles_of_an_info_mate(self) -> None:
        """`21 30425650 x N N. . . CHR2=21;END2=30425700`: Ensembl VEP 116.2 writes the pair
        `21:30425650-30425651` under `N.` and `N[21:30425700[` on ENST00000334068 and
        ENST00000390690; vep-rs the point `21:30425651` under the same alleles. The eight rows
        are four pairs; a point under a sequence allele pairs with nothing."""
        pair, point = "21:30425650-30425651", "21:30425651"
        rows = []
        for allele in ("N.", "N[21:30425700["):
            for feature, csq in (("ENST00000334068", "upstream_gene_variant"), ("ENST00000390690", "feature_truncation,coding_sequence_variant")):
                rows.append(self.row("missing_in_rust", pair, allele, feature, csq))
                rows.append(self.row("extra_in_rust", point, allele, feature, csq))
        rows.append(self.row("missing_in_rust", pair, "ACG", "ENST00000390690", "inframe_insertion"))
        rows.append(self.row("extra_in_rust", point, "ACG", "ENST00000390690", "inframe_insertion"))
        (ep, er, buckets), open_rows = self.run_pass(rows, self.scope())
        self.assertEqual((ep, er), (4, 4))
        self.assertEqual(buckets["dot_breakend_without_svtype_pair"], {"count": 8, "excluded": 8, "records": 2, "perl": 4, "rust": 4})
        self.assertEqual(sorted(r.split("\t")[3] for r in open_rows), ["ACG", "ACG"])

    def test_the_dot_breakend_pair_takes_the_spans_start_or_the_base_after_it(self) -> None:
        """The reference's start is POS + 1 when the ALT's length differs from the REF's and POS
        when it does not, its end POS + len(REF) - 1 (BaseVCF4.pm:264-275, 297-317), and vep-rs
        writes the point POS + 1 for every shape. Ensembl VEP 116.2 writes `21 30425650 x NAC .N
        . . .` at `21:30425651-30425652` and `21 30425650 x NACG .N . . .` at `21:30425651-30425653`
        (the point is the span's start), `21 30425650 x NAC NA. . . .` and `21 30425650 x NAC .AC
        . . .` at `21:30425650-30425652` and `21 30425650 x NACG NAC. . . .` at
        `21:30425650-30425653` (the point is the base after the start), each as
        `chromosome_breakpoint` on ENST00000334068 (`upstream_gene_variant`) and ENST00000390690
        (`coding_sequence_variant`); vep-rs writes `21:30425651 chromosome_breakpoint` on both for
        each, and each pair is set aside. A point two or more bases past the span's start is no
        pair: no REF length produces that shape."""
        from compare_vep_outputs import dot_breakend_points

        self.assertEqual(dot_breakend_points(("21", 30425650, 30425651)), ["21:30425650", "21:30425651"])
        self.assertEqual(dot_breakend_points(("21", 30425651, 30425653)), ["21:30425651", "21:30425652"])
        self.assertEqual(dot_breakend_points(("21", 30425650, 30425653)), ["21:30425650", "21:30425651"])
        for span in ("21:30425651-30425652", "21:30425651-30425653", "21:30425650-30425652", "21:30425650-30425653"):
            with self.subTest(span=span):
                rows = [self.row("missing_in_rust", span, "chromosome_breakpoint", "ENST00000334068", "upstream_gene_variant"),
                        self.row("missing_in_rust", span, "chromosome_breakpoint", "ENST00000390690", "coding_sequence_variant"),
                        self.row("extra_in_rust", "21:30425651", "chromosome_breakpoint", "ENST00000334068", "upstream_gene_variant"),
                        self.row("extra_in_rust", "21:30425651", "chromosome_breakpoint", "ENST00000390690", "coding_sequence_variant")]
                (ep, er, buckets), open_rows = self.run_pass(rows, self.scope())
                self.assertEqual((ep, er, open_rows), (2, 2, []))
                self.assertEqual(buckets["dot_breakend_without_svtype_pair"], {"count": 4, "excluded": 4, "records": 2, "perl": 2, "rust": 2})
        rows = [self.row("missing_in_rust", "21:30425649-30425651", "chromosome_breakpoint", "ENST00000390690", "coding_sequence_variant"),
                self.row("extra_in_rust", "21:30425651", "chromosome_breakpoint", "ENST00000390690", "coding_sequence_variant")]
        (ep, er, buckets), open_rows = self.run_pass(rows, self.scope())
        self.assertEqual((ep, er, len(open_rows)), (0, 0, 2))
        self.assertEqual(buckets["dot_breakend_without_svtype_pair"]["count"], 0)

    def test_a_mate_prefix_the_reference_does_not_resolve_is_the_synonym_shape(self) -> None:
        """`N[Chr21:30425700[` on CHROM 21: Ensembl VEP 116.2 writes the own end only, vep-rs the
        mate rows too (`21:30425651 N[Chr21:30425700[ ENST00000334068 upstream_gene_variant`),
        counted under the synonym class; `N[chr21:30425700[` is resolved and both engines agree."""
        loc = "21:30425651"
        breakends = {("f", loc): ({"N."}, {"N.", "N[Chr21:30425700["}, {"21"})}
        rows = [self.row("extra_in_rust", loc, "N[Chr21:30425700[", "ENST00000334068", "upstream_gene_variant"),
                self.row("extra_in_rust", loc, "N[Chr21:30425700[", "ENST00000390690", "feature_truncation,coding_sequence_variant")]
        (_, _, buckets), open_rows = self.run_pass(rows, self.scope(breakends=breakends))
        self.assertEqual(buckets["bnd_synonym_mate_rows_lost"], {"count": 2, "excluded": 0, "records": 1, "perl": 0, "rust": 0})
        self.assertEqual(len(open_rows), 2)
        breakends = {("f", loc): ({"N.", "N[chr21:30425700["}, {"N.", "N[chr21:30425700["}, {"21"})}
        rows = [self.row("extra_in_rust", loc, "N[chr21:30425700[", "ENST00000418755", "downstream_gene_variant")]
        (_, _, buckets), _ = self.run_pass(rows, self.scope(breakends=breakends))
        self.assertEqual(buckets["bnd_synonym_mate_rows_lost"]["count"], 0)


class DroppedRecordNameTests(unittest.TestCase):
    """The skipped-record class sets a row aside because the reference dropped its input
    record, not because the reference has no row at one Location.

    Both engines print the input record's name in the first column. A Location at which the
    engine under test wrote a structural row and the reference nothing is set aside only
    when no reference row of the same file carries the name of a row there. The rule is one
    for every engine: the rows below are written as vep-rs writes a structural record (the
    reference's own start and class word) and as an engine that prints the record from POS
    under its symbolic allele writes it. Every row is synthetic.
    """

    COLS = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
    SHARED = "snv_1\t21:100\t-\tg1\tT_SHARED\tTranscript\tintron_variant\n"

    @staticmethod
    def line(name: str, loc: str, allele: str, feature: str, csq: str = "feature_truncation") -> str:
        return f"{name}\t{loc}\t{allele}\tg1\t{feature}\tTranscript\t{csq}\n"

    def compare(self, perl: str, rust: str, files: dict[str, tuple[str, str]] | None = None):
        """Run compare mode over one output pair (or ``files``: name to both sides' rows) and
        return the class's statistics with the open discordance rows."""
        from compare_vep_outputs import main

        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, True)
        for sub in ("perl", "rust", "report"):
            (tmp / sub).mkdir()
        for name, (perl_rows, rust_rows) in (files or {"output.txt": (perl, rust)}).items():
            (tmp / "perl" / name).write_text(self.COLS + self.SHARED + perl_rows, encoding="utf-8")
            (tmp / "rust" / name).write_text(self.COLS + self.SHARED + rust_rows, encoding="utf-8")
        self.addCleanup(sys.set_int_max_str_digits, sys.get_int_max_str_digits())
        argv = sys.argv
        sys.argv = ["compare_vep_outputs.py", "--perl-dir", str(tmp / "perl"),
                    "--rust-dir", str(tmp / "rust"), "--report-dir", str(tmp / "report"),
                    "--assembly", "GRCh38", "--glob", "*.txt"]
        try:
            with mock.patch("builtins.print"):
                self.assertEqual(main(), 0)
        finally:
            sys.argv = argv
        adjusted = json.loads((tmp / "report" / "summary.json").read_text())["aggregate"]["adjusted"]
        open_rows = (tmp / "report" / "discordant_open.tsv").read_text().splitlines()[1:]
        return adjusted["excluded_categories"]["reference_skipped_record"], open_rows

    def test_a_symbolic_record_the_reference_wrote_one_base_on_is_not_set_aside(self) -> None:
        """The engine writes the record at POS under `<DEL>`, the reference at POS + 1 under
        `deletion`: the reference has no row at the engine's Location and wrote the record."""
        perl = self.line("sv_kept", "21:5001-5300", "deletion", "T1")
        rust = self.line("sv_kept", "21:5000-5300", "<DEL>", "T1") + self.line("sv_kept", "21:5000-5300", "<DEL>", "T2")
        bucket, open_rows = self.compare(perl, rust)
        self.assertEqual(bucket, {"count": 0, "excluded": 0, "records": 0, "perl": 0, "rust": 0})
        self.assertEqual(sum("\t21:5000-5300\t" in row for row in open_rows), 2, "the engine's rows stay charged")
        self.assertEqual(sum("\t21:5001-5300\t" in row for row in open_rows), 1, "and so does the reference's")

    def test_a_symbolic_record_the_reference_wrote_no_row_for_is_set_aside(self) -> None:
        rust = self.line("sv_dropped", "21:7000-27000100", "<DEL>", "T1") + self.line("sv_dropped", "21:7000-27000100", "<DEL>", "T2")
        bucket, open_rows = self.compare("", rust)
        self.assertEqual(bucket, {"count": 2, "excluded": 2, "records": 1, "perl": 0, "rust": 2})
        self.assertEqual(open_rows, [])

    def test_a_record_written_as_vep_rs_writes_it_follows_the_same_rule(self) -> None:
        """Under the reference's class word: a record the reference wrote with another end
        stays charged, and a record no reference row names is set aside."""
        perl = self.line("sv_kept", "21:5001-5300", "deletion", "T1")
        rust = self.line("sv_kept", "21:5001-5400", "deletion", "T1") + self.line("sv_dropped", "21:800-900", "CPX", "T3")
        bucket, open_rows = self.compare(perl, rust)
        self.assertEqual(bucket, {"count": 1, "excluded": 1, "records": 1, "perl": 0, "rust": 1})
        self.assertEqual(sum("\t21:5001-5400\t" in row for row in open_rows), 1)
        self.assertFalse(any("\t21:800-900\t" in row for row in open_rows))

    def test_a_record_without_an_id_is_told_by_its_generated_name(self) -> None:
        """Ensembl VEP names a record without an ID from its own line, `CHROM_POS_REF/ALT`
        with the POS as the input gives it. An engine that generates the same name and prints
        the record from POS is charged for a record the reference wrote one base on, and a
        record under a generated name no reference row carries is set aside."""
        perl = self.line("21_5000_k/<DEL>", "21:5001-5300", "deletion", "T1")
        rust = self.line("21_5000_k/<DEL>", "21:5000-5300", "<DEL>", "T1") + self.line("21_9000_k/<DEL>", "21:9000-29001000", "<DEL>", "T4")
        bucket, open_rows = self.compare(perl, rust)
        self.assertEqual(bucket, {"count": 1, "excluded": 1, "records": 1, "perl": 0, "rust": 1})
        self.assertEqual(sum("\t21:5000-5300\t" in row for row in open_rows), 1)
        self.assertFalse(any("\t21:9000-29001000\t" in row for row in open_rows))

    def test_a_sequence_allele_of_a_dropped_record_goes_with_it(self) -> None:
        rust = self.line("sv_mixed", "21:800-900", "CPX", "T3") + self.line("sv_mixed", "21:800-900", "-", "T3", "intron_variant")
        bucket, open_rows = self.compare("", rust)
        self.assertEqual(bucket, {"count": 2, "excluded": 2, "records": 1, "perl": 0, "rust": 2})
        self.assertEqual(open_rows, [])

    def test_a_location_shared_with_a_record_the_reference_wrote_stays_charged(self) -> None:
        """Two records at one Location of the engine, one the reference wrote elsewhere: the
        Location stands for both, and no row of a record the reference wrote is set aside."""
        perl = self.line("sv_kept", "21:5001-5300", "deletion", "T1")
        rust = self.line("sv_kept", "21:5000-5300", "<DEL>", "T1") + self.line("sv_dropped", "21:5000-5300", "<DUP>", "T2")
        bucket, open_rows = self.compare(perl, rust)
        self.assertEqual(bucket["excluded"], 0)
        self.assertEqual(sum("\t21:5000-5300\t" in row for row in open_rows), 2)

    def test_a_row_that_names_no_record_is_not_set_aside(self) -> None:
        """An empty first column, `.` and `-` name no input record, so nothing shows that the
        reference dropped it."""
        for name in ("", ".", "-"):
            with self.subTest(name=name):
                bucket, open_rows = self.compare("", self.line(name, "21:7000-27000100", "<DEL>", "T1"))
                self.assertEqual(bucket["excluded"], 0)
                self.assertEqual(len(open_rows), 1)

    def test_a_name_is_read_without_the_blanks_around_it(self) -> None:
        """The first column is read stripped, as the key's columns are: a name with blanks
        around it in one output is the name without them in the other."""
        perl = self.line("sv_kept", "21:5001-5300", "deletion", "T1")
        bucket, open_rows = self.compare(perl, self.line(" sv_kept ", "21:5000-5300", "<DEL>", "T1"))
        self.assertEqual(bucket["excluded"], 0)
        self.assertEqual(sum("\t21:5000-5300\t" in row for row in open_rows), 1)

    def test_a_name_is_looked_for_in_the_reference_output_of_the_same_file(self) -> None:
        """A reference row of another file under the same name is another input's record."""
        files = {
            "a.txt": ("", self.line("sv_1", "21:7000-27000100", "<DEL>", "T1")),
            "b.txt": (self.line("sv_1", "21:5001-5300", "deletion", "T1"), self.line("sv_1", "21:5001-5300", "deletion", "T1")),
        }
        bucket, open_rows = self.compare("", "", files=files)
        self.assertEqual(bucket, {"count": 1, "excluded": 1, "records": 1, "perl": 0, "rust": 1})
        self.assertEqual(open_rows, [])

    def test_a_reference_row_the_comparison_does_not_read_names_no_record(self) -> None:
        """The name is looked for on the rows the comparison keys: a reference row on a contig
        outside the canonical set is not one."""
        perl = self.line("sv_1", "GL000195.1:5001-5300", "deletion", "T1")
        bucket, _ = self.compare(perl, self.line("sv_1", "21:7000-27000100", "<DEL>", "T1"))
        self.assertEqual(bucket, {"count": 1, "excluded": 1, "records": 1, "perl": 0, "rust": 1})

    def test_outputs_with_no_candidate_location_are_not_read_again(self) -> None:
        from compare_vep_outputs import reference_named_locations

        with mock.patch("compare_vep_outputs.iter_record_names") as reader:
            self.assertEqual(reference_named_locations(Path("perl.txt"), Path("rust.txt"), set()), set())
        reader.assert_not_called()

    def test_the_named_locations_are_those_holding_a_name_the_reference_wrote(self) -> None:
        from compare_vep_outputs import ReferenceRecordScope, reference_named_locations

        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, True)
        perl, rust = tmp / "perl.txt", tmp / "rust.txt"
        perl.write_text(self.COLS + self.line("sv_kept", "21:5001-5300", "deletion", "T1"), encoding="utf-8")
        rust.write_text(
            self.COLS
            + self.line("sv_kept", "21:5000-5300", "<DEL>", "T1")
            + self.line("sv_dropped", "21:7000-27000100", "<DEL>", "T1")
            + self.line("sv_kept", "21:1-5", "<DEL>", "T9"),
            encoding="utf-8",
        )
        candidates = {"21:5000-5300", "21:7000-27000100"}
        self.assertEqual(reference_named_locations(perl, rust, candidates), {"21:5000-5300"})
        scope = ReferenceRecordScope()
        scope.absent = {("f", loc) for loc in candidates} | {("g", "21:5000-5300")}
        self.assertEqual(scope.absent_locations("f"), candidates)
        scope.record_reference_named("f", {"21:5000-5300"})
        self.assertEqual(scope.absent, {("f", "21:7000-27000100"), ("g", "21:5000-5300")})


class NumifiedCoordinateTests(unittest.TestCase):
    """The reference's rows at a coordinate it numified from a POS vep-rs refused, recognised
    without the input from the Location alone: a negative first coordinate
    (`non_numeric_pos_numified`) at a Location vep-rs wrote nothing for."""

    @staticmethod
    def key(loc: str, allele: str, feat: str, csq: str, ftype: str = "Transcript") -> str:
        return f"{loc}\t{allele}\t{feat}\t{ftype}\t{csq}"

    def merge(self, perl: list[str], rust: list[str]):
        from compare_vep_outputs import ReferenceRecordScope, compare_sorted_key_files

        scope = ReferenceRecordScope()
        with tempfile.TemporaryDirectory() as d:
            pp, rp = Path(d) / "p.txt", Path(d) / "r.txt"
            pp.write_text("".join(k + "\n" for k in sorted(perl)), encoding="utf-8")
            rp.write_text("".join(k + "\n" for k in sorted(rust)), encoding="utf-8")
            compare_sorted_key_files(pp, rp, file_name="f", discordant_out=None,
                                     track_consequence_buckets=False, record_scope=scope)
        return scope

    def test_a_location_with_a_negative_coordinate_is_parsed_and_is_an_integer(self) -> None:
        """Ensembl VEP 116.2 writes `21:-4-295` for `21 -5 x N <DEL> . . SVTYPE=DEL;SVLEN=-300`
        and `21:-2147483647--2147483348` for POS -2147483648: integer coordinates, outside the
        non-integer class, which keeps `21:30425651-NaN`, `21:30425949.0` and `21: 30425949`."""
        from compare_vep_outputs import _location_span, _non_integer_end

        self.assertEqual(_location_span("21:-4-295"), ("21", -4, 295))
        self.assertEqual(_location_span("21:-2147483647--2147483348"), ("21", -2147483647, -2147483348))
        self.assertEqual(_location_span("21:-5"), ("21", -5, -5))
        self.assertEqual(_location_span("21:18446744073709551617"), ("21", 2**64 + 1, 2**64 + 1))
        self.assertIsNone(_location_span("21:30425651-NaN"))
        for loc in ("21:-4-295", "21:-2147483647--2147483348", "21:-5", "21:30425651-30425949", "21:30425651"):
            self.assertFalse(_non_integer_end(loc), loc)
        for loc in ("21:30425651-NaN", "21:30425949.0", "21: 30425949", "21:30425949 ", "21:3e7", "21:30425651-30425950.7"):
            self.assertTrue(_non_integer_end(loc), loc)

    def test_the_scope_records_the_locations_only_the_reference_wrote(self) -> None:
        perl = [self.key("21:-4-295", "deletion", "-", "intergenic_variant", "-"),
                self.key("21:100", "A", "T1", "missense_variant"),
                self.key("21:200", "G", "T2", "intron_variant")]
        rust = [self.key("21:100", "A", "T1", "missense_variant"),
                self.key("21:200", "G", "T2", "upstream_gene_variant")]
        scope = self.merge(perl, rust)
        self.assertEqual(scope.reference_only, {("f", "21:-4-295")},
                         "a Location vep-rs wrote under another consequence is not reference-only")

    def test_a_reference_coordinate_past_the_buffers_32_bit_bound_is_an_overflow(self) -> None:
        """`21:30425651-2147483647 deletion - intergenic_variant` (Ensembl VEP 116.2 for
        `21 30425650 x N <DEL> . . SVTYPE=DEL;END=2147483647` under --max_sv_size -1) and
        `21:2147483649 N.` (POS 2147483648) enter the buffer's bounds as a mate at 2147478648
        does; a vep-rs-only Location there and an end at 2147478647 do not."""
        perl = [self.key("21:30425651-2147483647", "deletion", "-", "intergenic_variant", "-"),
                self.key("21:2147483649", "N.", "-", "intergenic_variant", "-"),
                self.key("21:30425651-2147478647", "deletion", "-", "intergenic_variant", "-")]
        rust = [self.key("21:2147483649", "N[21:30425700[", "T1", "upstream_gene_variant"),
                self.key("21:40000001-2147483647", "deletion", "T2", "transcript_ablation")]
        scope = self.merge(perl, rust)
        self.assertEqual(scope.mate_overflow, {("f", "21:30425651-2147483647"), ("f", "21:2147483649")})

    def test_a_reference_coordinate_at_the_buffers_bound_is_counted_through_the_digit_guard(self) -> None:
        """A coordinate reaches `BND_MATE_32_BIT_LIMIT` (2147478648, ten digits) only with ten or more
        digits, so a reference Location is parsed for the bound only when one of its coordinates has
        that many. The bound itself is counted as a point (`21:2147478648`) and as a span's end
        (`21:30425651-2147478648`, `2:5-2147478648`); one below it (`21:2147478647`), a nine-digit
        coordinate (`21:30425651-999999999`) and a ten-digit negative one
        (`21:-2147483647--2147483348`) are not."""
        from compare_vep_outputs import BND_MATE_32_BIT_LIMIT

        self.assertEqual(len(str(BND_MATE_32_BIT_LIMIT)), 10)
        perl = [self.key("21:2147478648", "N.", "-", "intergenic_variant", "-"),
                self.key("21:30425651-2147478648", "deletion", "-", "intergenic_variant", "-"),
                self.key("2:5-2147478648", "deletion", "-", "intergenic_variant", "-"),
                self.key("21:2147478647", "N.", "-", "intergenic_variant", "-"),
                self.key("21:30425651-999999999", "deletion", "-", "intergenic_variant", "-"),
                self.key("21:-2147483647--2147483348", "deletion", "-", "intergenic_variant", "-")]
        scope = self.merge(perl, [])
        self.assertEqual(scope.mate_overflow,
                         {("f", "21:2147478648"), ("f", "21:30425651-2147478648"), ("f", "2:5-2147478648")})

    def test_the_mate_chromosome_as_the_reference_reads_it(self) -> None:
        from compare_vep_outputs import reference_mate_chromosome

        for allele, name in (("N[chr21:30425700[", "21"), ("N[Chr21:30425700[", "Chr21"), ("N[CHR21:30425700[", "CHR21"),
                             ("N[21:30425700[", "21"), ("]NC_000022.11:26037668]N", "NC_000022.11"), ("N[chrM:100[", "M")):
            self.assertEqual(reference_mate_chromosome(allele), name, allele)
        self.assertIsNone(reference_mate_chromosome("N."))
        self.assertIsNone(reference_mate_chromosome("deletion"))


class ReferenceReleaseFixtureTests(unittest.TestCase):
    """Compare mode on the committed fixture under each release.

    `testdata/reference_release/expected_115_2/` is the report the 115.2 registry writes on
    this fixture, the published comparator's own output; under `--reference-release 115.2`
    every report file must equal it byte for byte (summary.json after its run-specific
    timestamp and paths), so the published figures reproduce. Under 116.2 the record the
    reference wrote nothing for is set aside whole and the 115.2 scope class is not
    applied.
    """

    FIXTURE = Path(__file__).resolve().parent / "testdata" / "reference_release"
    VOLATILE = ("generated_at_utc", "perl_dir", "rust_dir", "report_dir")

    def _run(self, release: str) -> Path:
        import shutil
        import sys
        from compare_vep_outputs import main

        out = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, out, True)
        self.addCleanup(sys.set_int_max_str_digits, sys.get_int_max_str_digits())
        argv = sys.argv
        sys.argv = ["compare_vep_outputs.py", "--perl-dir", str(self.FIXTURE / "perl"),
                    "--rust-dir", str(self.FIXTURE / "rust"), "--report-dir", str(out),
                    "--assembly", "GRCh37", "--reference-release", release]
        try:
            self.assertEqual(main(), 0)
        finally:
            sys.argv = argv
        return out

    def test_an_output_without_an_input_vcf_under_input_dir_warns_once(self) -> None:
        """`--input-dir` holding no VCF named after an output's stem writes one warning per such
        stem to stderr and the comparison goes on, recognising a gVCF variant site from an
        ID-less record's generated name alone."""
        import contextlib
        import io
        import shutil
        import sys

        from compare_vep_outputs import main, no_input_vcf_warning

        out, inputs = Path(tempfile.mkdtemp()), Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, out, True)
        self.addCleanup(shutil.rmtree, inputs, True)
        self.addCleanup(sys.set_int_max_str_digits, sys.get_int_max_str_digits())
        argv = sys.argv
        sys.argv = ["compare_vep_outputs.py", "--perl-dir", str(self.FIXTURE / "perl"),
                    "--rust-dir", str(self.FIXTURE / "rust"), "--report-dir", str(out),
                    "--assembly", "GRCh37", "--input-dir", str(inputs)]
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err):
                self.assertEqual(main(), 0)
        finally:
            sys.argv = argv
        self.assertEqual(err.getvalue().splitlines(), [no_input_vcf_warning("output", inputs.resolve())])
        self.assertEqual(no_input_vcf_warning("s07", Path("in")),
                         "WARNING: [compare_vep_outputs] no input VCF for s07 under in; "
                         "gVCF variant sites are recognised from an ID-less record's generated name only")

    def test_115_2_reproduces_the_published_report_byte_for_byte(self) -> None:
        import json

        out = self._run("115.2")
        expected = self.FIXTURE / "expected_115_2"
        for name in ("summary.md", "discordant.tsv", "discordant_open.tsv"):
            with self.subTest(file=name):
                self.assertEqual((out / name).read_bytes(), (expected / name).read_bytes())
        got = json.loads((out / "summary.json").read_text())
        want = json.loads((expected / "summary.json").read_text())
        for k in self.VOLATILE:
            got[k] = want[k]
        self.assertEqual(got, want)
        self.assertNotIn("reference_release", got["aggregate"]["adjusted"])

    def test_116_2_sets_aside_the_skipped_record_and_not_the_batch_class(self) -> None:
        import json

        out = self._run("116.2")
        adj = json.loads((out / "summary.json").read_text())["aggregate"]["adjusted"]
        self.assertEqual(adj["reference_release"], "116.2")
        self.assertEqual(adj["excluded_categories"]["reference_skipped_record"],
                         {"count": 2, "excluded": 2, "records": 1, "perl": 0, "rust": 2})
        self.assertNotIn("sv_transcript_selection", adj["excluded_categories"])
        self.assertEqual((adj["excluded_perl_tuples"], adj["excluded_rust_tuples"]), (2, 4))
        self.assertEqual(set(adj["divergence_classes"]) >= {"reference_skipped_record", "splice_lastwrite_swap"}, True)
        open_rows = (out / "discordant_open.tsv").read_text().splitlines()[1:]
        self.assertFalse(any("\t21:800-900\t" in r for r in open_rows), "the skipped <CPX> is set aside")
        self.assertEqual(sum("\t21:18457513-32484691\t" in r for r in open_rows), 2,
                         "the surplus transcripts on a record the reference kept stay charged")
        self.assertTrue(any("\ttandem_repeat\t" in r for r in open_rows))
        md = (out / "summary.md").read_text()
        self.assertIn("- Reference release: Ensembl VEP 116.2", md)
        self.assertIn("## Divergence classes", md)
        self.assertIn("| `reference_skipped_record` | sets aside |", md)


class GvcfSiteRowsFixtureTests(unittest.TestCase):
    """`scripts/parity/gvcf_site_rows_fixture.tsv` holds the (Location, Allele) rows vep-rs wrote
    for 300 generated gVCF variant sites (`scripts/parity/make_gvcf_site_fixture.py`: a REF of 1
    to 4 bases, 1 to 3 real ALTs of 1 to 5 bases, `<NON_REF>` anywhere in the list, the anchor
    base shared in every second site). The comparator's prediction of those rows is held to the
    binary's output site by site."""

    FIXTURE = Path(__file__).resolve().parent.parent / "parity" / "gvcf_site_rows_fixture.tsv"

    def test_gvcf_site_rows_predict_every_row_the_binary_wrote(self) -> None:
        """For every site, `gvcf_site_rows` over the ALT column, placed by `sequence_location` and
        read through `canonical_location`, is exactly the fixture's (Location, Allele) set: no
        allele the binary wrote is missed and none is predicted that it did not write."""
        from compare_vep_outputs import (
            canonical_location,
            gvcf_site_rows,
            sequence_location,
        )

        written: dict[tuple[int, str, str], set[tuple[str, str]]] = {}
        with self.FIXTURE.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, ["pos", "ref", "alt_column", "allele", "location"])
            for row in reader:
                site = (int(row["pos"]), row["ref"], row["alt_column"])
                written.setdefault(site, set()).add((row["location"], row["allele"]))
        self.assertEqual(len(written), 300)
        self.assertTrue(all(25_000_000 <= pos <= 46_000_000 for pos, _, _ in written))
        self.assertTrue(all("<NON_REF>" in alt_column.split(",") for _, _, alt_column in written))
        mismatches = []
        for (pos, ref, alt_column), rows in sorted(written.items()):
            predicted = {(canonical_location(sequence_location("21", start, end)), allele)
                         for allele, start, end in gvcf_site_rows(pos, ref, alt_column.split(","))}
            if predicted != rows:
                mismatches.append((pos, ref, alt_column, sorted(predicted), sorted(rows)))
        self.assertEqual(mismatches, [])


ROW_COLS = (
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\t"
    "CDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
)


def stop_codon_row(allele: str, consequence: str, cdna: str, cds: str, location: str = "19:13208046-13208050",
                   feature_type: str = "Transcript", name: str = "1993671", feature: str = "ENST00000635895") -> str:
    """One fabricated row of the default format, as the format comparator's stop codon tests
    write it."""
    return (f"{name}\t{location}\t{allele}\tENSG1\t{feature}\t{feature_type}\t{consequence}\t{cdna}\t{cds}"
            "\t-\t-\t-\t-\tIMPACT=LOW;STRAND=-1\n")


def compare_row_outputs(test: unittest.TestCase, perl: str, rust: str, *options: str) -> tuple[dict, list[str]]:
    """Run compare mode over one fabricated output pair, with the options a concordance run
    passes, and return the summary's adjusted block with the open discordance rows."""
    from compare_vep_outputs import main

    tmp = Path(tempfile.mkdtemp())
    test.addCleanup(shutil.rmtree, tmp, True)
    for sub, rows in (("perl", perl), ("rust", rust)):
        (tmp / sub).mkdir()
        (tmp / sub / "output.txt").write_text(ROW_COLS + rows, encoding="utf-8")
    test.addCleanup(sys.set_int_max_str_digits, sys.get_int_max_str_digits())
    argv = ["compare_vep_outputs.py", "--perl-dir", str(tmp / "perl"), "--rust-dir", str(tmp / "rust"),
            "--report-dir", str(tmp / "report"), "--assembly", "GRCh38", *options]
    with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print"):
        test.assertEqual(main(), 0)
    adjusted = json.loads((tmp / "report" / "summary.json").read_text())["aggregate"]["adjusted"]
    return adjusted, (tmp / "report" / "discordant_open.tsv").read_text().splitlines()[1:]


@contextmanager
def row_pair_rules(*rules: ExclusionRule) -> Iterator[None]:
    """The 116.2 registry with ``rules`` as its row-aware pair rules."""
    registry = dataclasses.replace(cvo.REFERENCE_REGISTRIES["116.2"], row_pair_rules=rules)
    with mock.patch.dict(cvo.REFERENCE_REGISTRIES, {"116.2": registry}), \
            mock.patch.dict(cvo.PAIR_RULE_DEFINITIONS, {r.bucket: "A fabricated shape." for r in rules}), \
            mock.patch.dict(cvo.PAIR_RULE_CITATIONS, {r.taxonomy_class: "Fabricated.pm:1" for r in rules}):
        yield


def fabricated_row_rule(predicate: Callable[..., bool]) -> ExclusionRule:
    return ExclusionRule(bucket="fabricated_row_shape", predicate=predicate, classifies=True, excludes=True,
                         taxonomy_class="fabricated_row_shape", row_aware=True)


class RowAwarePairRuleTests(unittest.TestCase):
    """A pair rule may read, beside the two consequence sets and the allele, the first
    column, the Location, the cDNA_position and the CDS_position of the two rows of an
    unmatched pair, fetched from the two outputs the comparator is given. Every row is
    fabricated."""

    def test_a_row_aware_rule_is_handed_the_named_columns_of_the_two_rows_of_an_unmatched_pair(self) -> None:
        seen: list[tuple] = []

        def shape(e, a) -> bool:
            seen.append((e, a))
            return True

        name = "19_13208045_GACCTG/G"
        matched = stop_codon_row("-", "missense_variant", "10-12", "4-6", feature="ENST00000000001")
        perl = matched + stop_codon_row("-", "stop_retained_variant", "?-7057", "-", name=name)
        rust = matched + stop_codon_row("-", "stop_lost", "?-7057/7140", "-/7023", name=name)
        with row_pair_rules(fabricated_row_rule(shape)):
            adjusted, open_rows = compare_row_outputs(self, perl, rust)
        ((e, a),) = seen
        self.assertEqual(
            (e.key, e.location, e.allele, e.consequence, e.fields, e.indel),
            ("19:13208046-13208050|-|ENST00000635895|Transcript", "19:13208046-13208050", "-",
             "stop_retained_variant", {"cDNA_position": "?-7057", "CDS_position": "-"}, True),
        )
        self.assertEqual((a.consequence, a.fields, a.indel),
                         ("stop_lost", {"cDNA_position": "?-7057/7140", "CDS_position": "-/7023"}, True))
        self.assertEqual(adjusted["excluded_categories"]["fabricated_row_shape"], {"count": 1, "excluded": 1})
        self.assertEqual(list(adjusted["excluded_categories"])[len(CLASSIFICATION_CHAIN)], "fabricated_row_shape")
        self.assertEqual(list(adjusted["divergence_classes"])[len(EXCLUSION_REGISTRY)], "fabricated_row_shape")
        self.assertEqual((adjusted["f1"], adjusted["excluded_perl_tuples"], adjusted["excluded_rust_tuples"]), (1.0, 1, 1))
        self.assertEqual(open_rows, [])

    def test_a_record_under_an_identifier_leaves_the_shape_to_the_location(self) -> None:
        seen: list[tuple] = []
        perl = stop_codon_row("-", "stop_retained_variant", "7055-7059", "6978-6982")
        rust = stop_codon_row("-", "coding_sequence_variant", "7055-7059", "6978-6982")
        with row_pair_rules(fabricated_row_rule(lambda e, a: seen.append((e.indel, a.indel, e.location)) or False)):
            adjusted, open_rows = compare_row_outputs(self, perl, rust)
        self.assertEqual(seen, [(None, None, "19:13208046-13208050")])
        self.assertEqual(adjusted["excluded_categories"]["fabricated_row_shape"], {"count": 0, "excluded": 0})
        self.assertEqual(len(open_rows), 2)

    def test_only_an_unmatched_pair_no_earlier_rule_takes_reaches_a_row_aware_rule(self) -> None:
        seen: list[tuple] = []
        matched = stop_codon_row("-", "missense_variant", "10-12", "4-6", feature="ENST00000000001")
        perl = (matched
                + stop_codon_row("-", "stop_retained_variant", "20-22", "14-16", feature="ENST00000000002")
                + stop_codon_row("-", "intron_variant", "-", "-", feature="ENST00000000004"))
        rust = (matched
                + stop_codon_row("-", "stop_lost", "20-22", "14-16", feature="ENST00000000003")
                + stop_codon_row("-", "intron_variant,splice_region_variant", "-", "-", feature="ENST00000000004"))
        with row_pair_rules(fabricated_row_rule(lambda e, a: seen.append((e, a)) or True)):
            adjusted, open_rows = compare_row_outputs(self, perl, rust)
        self.assertEqual(seen, [], "a matched tuple, a tuple one engine alone writes and a pair an earlier rule takes")
        self.assertEqual(adjusted["excluded_categories"]["fabricated_row_shape"], {"count": 0, "excluded": 0})
        self.assertEqual(adjusted["excluded_categories"]["intronic_overcall_swap"], {"count": 1, "excluded": 1})
        self.assertEqual(sorted(row.split("\t")[4] for row in open_rows), ["ENST00000000002", "ENST00000000003"])

    def test_the_rows_of_one_key_stand_or_fall_together(self) -> None:
        perl = (stop_codon_row("-", "stop_retained_variant", "10-12", "4-6", name="first")
                + stop_codon_row("-", "stop_retained_variant", "20-22", "14-16", name="second")
                + stop_codon_row("-", "stop_retained_variant", "10-12", "4-6", name="third"))
        rust = stop_codon_row("-", "stop_lost", "10-12", "4-6", name="first")
        for refused, expected in ((None, ({"count": 1, "excluded": 1}, 0)), ("20-22", ({"count": 0, "excluded": 0}, 2))):
            seen: list[str] = []

            def shape(e, a, refused=refused, seen=seen) -> bool:
                seen.append(e.fields["cDNA_position"])
                return e.fields["cDNA_position"] != refused

            with row_pair_rules(fabricated_row_rule(shape)):
                adjusted, open_rows = compare_row_outputs(self, perl, rust)
            self.assertEqual(seen, ["10-12", "20-22"], "each distinct row of the reference's tuple, once")
            self.assertEqual((adjusted["excluded_categories"]["fabricated_row_shape"], len(open_rows)), expected, refused)

    def test_a_row_of_fewer_than_nine_columns_is_put_to_no_row_aware_rule(self) -> None:
        seen: list[tuple] = []
        short = "1993671\t19:13208046-13208050\t-\tENSG1\tENST00000635895\tTranscript\t{}\t?-7057\n"
        cols = "\t".join(ROW_COLS.split("\t")[:8]) + "\n"
        with row_pair_rules(fabricated_row_rule(lambda e, a: seen.append((e, a)) or True)), \
                mock.patch(f"{__name__}.ROW_COLS", cols):
            adjusted, open_rows = compare_row_outputs(self, short.format("stop_retained_variant"), short.format("stop_lost"))
        self.assertEqual(seen, [])
        self.assertEqual((adjusted["excluded_categories"]["fabricated_row_shape"], len(open_rows)), ({"count": 0, "excluded": 0}, 2))

    def test_without_a_row_aware_rule_no_output_is_read_again_and_the_report_is_the_registry_s(self) -> None:
        perl = stop_codon_row("-", "stop_retained_variant", "?-7057", "-")
        rust = stop_codon_row("-", "stop_lost", "?-7057", "-")
        with row_pair_rules(), mock.patch("compare_vep_outputs.fetch_pair_rows", side_effect=AssertionError("read")):
            adjusted, open_rows = compare_row_outputs(self, perl, rust)
        self.assertEqual(len(open_rows), 2)
        chain = [r.bucket for r in CLASSIFICATION_CHAIN] + [r.bucket for r in cvo.REFERENCE_116_2_ONE_SIDED_RULES]
        self.assertEqual(list(adjusted["excluded_categories"]), chain)
        self.assertEqual(adjusted["divergence_classes"], json.loads(json.dumps(cvo.divergence_class_definitions("116.2"))))

    def test_the_115_2_registry_asks_no_row_aware_rule(self) -> None:
        seen: list[tuple] = []
        self.assertEqual(cvo.REFERENCE_REGISTRIES["115.2"].row_pair_rules, ())
        perl = stop_codon_row("-", "stop_retained_variant", "?-7057", "-")
        rust = stop_codon_row("-", "stop_lost", "?-7057", "-")
        with row_pair_rules(fabricated_row_rule(lambda e, a: seen.append((e, a)) or True)):
            adjusted, open_rows = compare_row_outputs(self, perl, rust, "--reference-release", "115.2")
        self.assertEqual((seen, len(open_rows)), ([], 2))
        self.assertNotIn("fabricated_row_shape", adjusted["excluded_categories"])

    def test_a_row_aware_rule_matches_nothing_without_the_rows(self) -> None:
        rule = fabricated_row_rule(lambda e, a: True)
        self.assertFalse(rule.matches("stop_retained_variant", "stop_lost", "-"))
        self.assertFalse(rule.matches("stop_retained_variant", "stop_lost", "-", ([], [])))
        row = cvo.pair_row(stop_codon_row("-", "stop_lost", "?-7057", "-").rstrip("\n").split("\t")[:9])
        self.assertTrue(rule.matches("stop_retained_variant", "stop_lost", "-", ([row], [row])))
        self.assertEqual(cvo.report_class_definitions("115.2"), cvo.divergence_class_definitions("115.2"))
        registry = cvo.REFERENCE_REGISTRIES["116.2"]
        self.assertEqual(
            list(cvo.report_class_definitions("116.2")),
            [r.bucket for r in (*registry.pair_rules, *registry.row_pair_rules, *registry.one_sided_rules)],
        )

    def test_the_fetch_reads_the_two_outputs_alone_and_only_for_the_pairs_it_is_given(self) -> None:
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, True)
        rows = []
        for sub, source, term in (("perl", "missing_in_rust", "stop_retained_variant"), ("rust", "extra_in_rust", "stop_lost")):
            (tmp / sub).mkdir()
            (tmp / sub / "output.txt").write_text(
                ROW_COLS + stop_codon_row("-", term, "?-7057", "-", location="19:013208046-+13208050"), encoding="utf-8")
            rows.append({"file_name": "output.txt", "source": source, "location": "19:13208046-13208050", "allele": "-",
                         "feature": "ENST00000635895", "feature_type": "Transcript", "consequence_set": term})
        self.assertEqual(cvo.fetch_pair_rows(rows, [], tmp / "absent", tmp / "absent"), {})
        fetched = cvo.fetch_pair_rows(rows, [(0, 1)], tmp / "perl", tmp / "rust")
        self.assertEqual(sorted(fetched), [0, 1])
        self.assertEqual([r.location for r in fetched[0]], ["19:013208046-+13208050"], "the Location as the row writes it")
        self.assertEqual(fetched[1][0].consequence, "stop_lost")


class RowPairRuleChecks:
    """The four checks of one row-aware pair rule of the 116.2 registry, on the fabricated
    rows the format comparator's test of the same class builds: the rule's pair is set
    aside; the pair with one further differing term, and the pair with a position column
    outside the rule's condition, stay differences; and beside the pair a tuple vep-rs lacks
    and a tuple of other terms stay differences."""

    bucket: str
    reference: str
    candidate: str
    further_term: str
    other_column: tuple[str, str]

    MATCHED = stop_codon_row("-", "missense_variant", "10-12", "4-6", location="21:5000-5002", feature="ENST00000000001")

    def scored(self, perl: str, rust: str) -> tuple[dict, float, list[str]]:
        adjusted, open_rows = compare_row_outputs(self, self.MATCHED + perl, self.MATCHED + rust)
        return adjusted["excluded_categories"][self.bucket], adjusted["f1"], open_rows

    def test_the_pair_of_rows_is_set_aside(self) -> None:
        self.assertEqual(self.scored(self.reference, self.candidate), ({"count": 1, "excluded": 1}, 1.0, []))

    def test_one_further_differing_term_stays_a_difference(self) -> None:
        stats, f1, open_rows = self.scored(self.reference, self.further_term)
        self.assertEqual((stats, len(open_rows)), ({"count": 0, "excluded": 0}, 2))
        self.assertLess(f1, 1.0)

    def test_a_column_outside_the_condition_stays_a_difference(self) -> None:
        stats, f1, open_rows = self.scored(*self.other_column)
        self.assertEqual((stats, len(open_rows)), ({"count": 0, "excluded": 0}, 2))
        self.assertLess(f1, 1.0)

    def test_a_missing_tuple_and_a_wrong_tuple_stay_differences(self) -> None:
        missing = stop_codon_row("-", "intron_variant", "-", "-", location="21:5000-5002", feature="ENST00000000002")
        wrong = "21:7000-7002"
        stats, f1, open_rows = self.scored(
            self.reference + missing + stop_codon_row("-", "missense_variant", "30-32", "24-26", location=wrong),
            self.candidate + stop_codon_row("-", "synonymous_variant", "30-32", "24-26", location=wrong),
        )
        self.assertEqual(stats, {"count": 1, "excluded": 1})
        self.assertEqual(sorted(row.split("\t")[2] for row in open_rows), ["21:5000-5002", wrong, wrong])
        self.assertLess(f1, 1.0)


class StopCodonOverlapGenomicSpanRuleTests(RowPairRuleChecks, unittest.TestCase):
    """`stop_codon_overlap_genomic_span`: the reference's coding_sequence_variant against
    vep-rs's stop_retained_variant beside splice_region_variant, on a row whose stop codon
    test is made."""

    bucket = "stop_codon_overlap_genomic_span"
    terms = "splice_region_variant,{},NMD_transcript_variant"
    reference = stop_codon_row("-", terms.format("coding_sequence_variant"), "1705", "1630", "3:36993663")
    candidate = stop_codon_row("-", terms.format("stop_retained_variant"), "1705", "1630", "3:36993663")
    further_term = stop_codon_row("-", terms.format("stop_retained_variant,stop_lost"), "1705", "1630", "3:36993663")
    other_column = (
        stop_codon_row("-", terms.format("coding_sequence_variant"), "1705-?", "1630-?", "3:36993663"),
        stop_codon_row("-", terms.format("stop_retained_variant"), "1705-?", "1630-?", "3:36993663"),
    )


class StopRetainedWithoutStopCodonTestRuleTests(RowPairRuleChecks, unittest.TestCase):
    """`stop_retained_without_stop_codon_test`: the reference's stop_retained_variant against
    nothing, coding_sequence_variant or frameshift_variant, on rows whose stop codon test is
    not made."""

    bucket = "stop_retained_without_stop_codon_test"
    terms = "splice_acceptor_variant,{}"
    reference = stop_codon_row("-", terms.format("stop_retained_variant"), "?-7057", "-")
    candidate = stop_codon_row("-", terms.format("coding_sequence_variant"), "?-7057", "-")
    further_term = stop_codon_row("-", terms.format("coding_sequence_variant,frameshift_variant"), "?-7057", "-")
    other_column = (reference, stop_codon_row("-", terms.format("coding_sequence_variant"), "7055-7059", "6978-6982"))


class StopRetainedThreeBaseStopCodonTestRuleTests(RowPairRuleChecks, unittest.TestCase):
    """`stop_retained_three_base_stop_codon_test`: the reference's stop_retained_variant
    against coding_sequence_variant, or against frameshift_variant with or without
    stop_gained, on rows whose stop codon test is made."""

    bucket = "stop_retained_three_base_stop_codon_test"
    terms = "{},3_prime_UTR_variant"
    location = "2:190970701-190970707"
    reference = stop_codon_row("-", terms.format("stop_retained_variant"), "2558-2564", "2249-?", location)
    candidate = stop_codon_row("-", terms.format("coding_sequence_variant"), "2558-2564", "2249-?", location)
    further_term = stop_codon_row("-", terms.format("coding_sequence_variant,frameshift_variant"), "2558-2564", "2249-?", location)
    other_column = (reference, stop_codon_row("-", terms.format("coding_sequence_variant"), "2558-?", "2249-?", location))


class StopCodonRemovedWithoutStopCodonTestRuleTests(RowPairRuleChecks, unittest.TestCase):
    """`stop_codon_removed_without_stop_codon_test`: the reference's stop_retained_variant
    against vep-rs's stop_lost, on rows without a cDNA position at one end or a CDS start."""

    bucket = "stop_codon_removed_without_stop_codon_test"
    terms = "splice_acceptor_variant,{}"
    reference = stop_codon_row("-", terms.format("stop_retained_variant"), "?-7057", "-")
    candidate = stop_codon_row("-", terms.format("stop_lost"), "?-7057", "-")
    further_term = stop_codon_row("-", terms.format("stop_lost,frameshift_variant"), "?-7057", "-")
    other_column = (
        stop_codon_row("-", terms.format("stop_retained_variant"), "7055-7059", "6978-6982"),
        stop_codon_row("-", terms.format("stop_lost"), "7055-7059", "6978-6982"),
    )


class StopCodonOverlapInsertionRangeRuleTests(RowPairRuleChecks, unittest.TestCase):
    """`stop_codon_overlap_insertion_range`: the reference's coding_sequence_variant, or its
    frameshift_variant with or without stop_gained, against vep-rs's stop_retained_variant on
    an insertion's row with both cDNA positions and a CDS start."""

    bucket = "stop_codon_overlap_insertion_range"
    location = "19:13208046-13208047"
    reference = stop_codon_row("TTGA", "inframe_insertion,coding_sequence_variant", "7054-7055", "6977-6978", location)
    candidate = stop_codon_row("TTGA", "inframe_insertion,stop_retained_variant", "7054-7055", "6977-6978", location)
    further_term = stop_codon_row("TTGA", "inframe_insertion,stop_retained_variant,stop_lost", "7054-7055", "6977-6978", location)
    other_column = (
        stop_codon_row("TTGA", "inframe_insertion,coding_sequence_variant", "7054-?", "6977-?", location),
        stop_codon_row("TTGA", "inframe_insertion,stop_retained_variant", "7054-?", "6977-?", location),
    )


class SpliceTermPastTheRecordSpanRuleTests(RowPairRuleChecks, unittest.TestCase):
    """`splice_term_past_the_record_span`: the reference's splice_donor_variant and
    intron_variant against vep-rs's splice_region_variant, on a row under an identifier, which
    carries no record: its three positions lie in the transcript's sequence under an allele of
    five."""

    bucket = "splice_term_past_the_record_span"
    location = "7:100998-101000"
    reference = stop_codon_row("GTTCA", "splice_donor_variant,frameshift_variant,intron_variant", "998-1000", "898-900", location)
    candidate = stop_codon_row("GTTCA", "frameshift_variant,splice_region_variant", "998-1000", "898-900", location)
    further_term = stop_codon_row("GTTCA", "stop_gained,frameshift_variant,splice_region_variant", "998-1000", "898-900", location)
    other_column = (
        stop_codon_row("GTTCA", "splice_donor_variant,frameshift_variant,intron_variant", "999-?", "899-?", location),
        stop_codon_row("GTTCA", "frameshift_variant,splice_region_variant", "999-?", "899-?", location),
    )


def stop_codon_pair_row(row: str) -> cvo.PairRow:
    """The row a row-aware pair rule reads of one fabricated row of the default format."""
    return cvo.pair_row(row.rstrip("\n").split("\t")[:9])


class InsertionRangeRuleRowTests(unittest.TestCase):
    """`stop_codon_overlap_insertion_range` takes an insertion between two neighbouring
    positions and no wider row: a Location of two neighbouring positions on every row; where
    the first column spells the record's alleles, a reference allele that is empty once the
    shared first base is dropped; under an identifier, an allele of another length than two.
    Under an identifier two bases replaced by a sequence of another length have an
    insertion's Location and allele, and are set aside with the insertions."""

    TWO = "19:13208046-13208047"
    # (first column, Location, allele, the shape)
    WIDER = (
        ("19_13208046_AC/GTT", TWO, "GTT", "two bases replaced by three"),
        ("19_13208046_AC/AGT", TWO, "AGT", "two bases replaced by three, the first base shared"),
        ("19_13208046_AC/A/T", TWO, "A", "a deletion written as a sequence"),
        ("19_13208046_A/AC/T", "19:13208046", "AC", "an insertion kept with its shared first base"),
        ("19_13208046_-/TGA", "19:13208046", "TGA", "an insertion's alleles at one position"),
        ("19_13208046_AC/A/T", TWO, "GTT", "an allele no alternate names beside a reference allele of two"),
        ("1993671", TWO, "CG", "a two-base allele under an identifier"),
    )
    INSERTED = (
        ("19_13208046_-/TGA", TWO, "TGA", "an insertion, the reference allele empty"),
        ("19_13208046_A/ATGA", TWO, "TGA", "an insertion spelt with its shared first base"),
        ("19_13208046_-/CG/T", TWO, "CG", "an insertion of two bases, the alleles spelt"),
        ("19_13208046_-/CG/T", TWO, "TGA", "an allele no alternate names beside an empty reference allele"),
        ("1993671", TWO, "TGA", "an insertion under an identifier"),
        ("1993671", TWO, "G", "under an identifier, two bases replaced by one"),
        ("1993671", TWO, "CGT", "under an identifier, two bases replaced by three"),
    )

    @staticmethod
    def rows(name: str, location: str, allele: str) -> tuple[str, str]:
        return tuple(  # type: ignore[return-value]
            stop_codon_row(allele, terms, "7054-7055", "6977-6978", location, name=name)
            for terms in ("coding_sequence_variant", "stop_retained_variant")
        )

    def taken(self, name: str, location: str, allele: str) -> bool:
        reference, candidate = self.rows(name, location, allele)
        return cvo.stop_codon_overlap_insertion_range_pair(stop_codon_pair_row(reference), stop_codon_pair_row(candidate))

    def test_an_insertion_between_two_positions_is_taken_and_a_wider_row_stays_a_difference(self) -> None:
        self.assertEqual([shape for *row, shape in self.INSERTED if not self.taken(*row)], [])
        self.assertEqual([shape for *row, shape in self.WIDER if self.taken(*row)], [])

    def test_the_comparison_sets_the_insertions_aside_and_leaves_the_wider_rows_open(self) -> None:
        perl = rust = ""
        for i, (name, location, allele, _) in enumerate((*self.WIDER, *self.INSERTED)):
            chromosome, _, span = location.partition(":")
            start, _, end = span.partition("-")
            moved = f"{chromosome}:{int(start) + 100 * i}" + (f"-{int(end) + 100 * i}" if end else "")
            reference, candidate = self.rows(name, moved, allele)
            perl, rust = perl + reference, rust + candidate
        adjusted, open_rows = compare_row_outputs(self, perl, rust)
        self.assertEqual(adjusted["excluded_categories"]["stop_codon_overlap_insertion_range"],
                         {"count": len(self.INSERTED), "excluded": len(self.INSERTED)})
        self.assertEqual(len(open_rows), 2 * len(self.WIDER))

    def test_the_row_carries_whether_its_alleles_make_an_insertion(self) -> None:
        self.assertIn("insertion", [f.name for f in dataclasses.fields(cvo.PairRow)])
        for (name, location, allele, shape), insertion in zip(
            (*self.WIDER, *self.INSERTED),
            (False, False, False, False, True, False, None, True, True, True, True, None, None, None),
        ):
            self.assertIs(stop_codon_pair_row(self.rows(name, location, allele)[0]).insertion, insertion, shape)
        self.assertIs(cvo.reference_allele_is_empty("A", ["AT", "*"]), True)
        self.assertIs(cvo.reference_allele_is_empty("A", ["AT", "AG"]), True)
        self.assertIs(cvo.reference_allele_is_empty("TCT", ["TCTC"]), False, "a base stays past the shared first one")
        self.assertIs(cvo.reference_allele_is_empty("-", ["T"], chopped=True), True)
        self.assertIs(cvo.reference_allele_is_empty("A", ["AT"], chopped=True), False)
        self.assertIs(cvo.reference_allele_is_empty("GACCTG", ["G"]), False)
        self.assertIs(cvo.reference_allele_is_empty("A", ["T", "G"]), False)


class ThreeBaseRuleSubsetTests(unittest.TestCase):
    """`stop_retained_three_base_stop_codon_test` takes the reference's stop_retained_variant
    alone, or beside one of inframe_insertion, inframe_deletion and protein_altering_variant,
    against frameshift_variant with or without stop_gained. Two or three of them stay a
    difference."""

    SEVERAL = (
        "inframe_insertion,inframe_deletion",
        "inframe_insertion,protein_altering_variant",
        "inframe_deletion,protein_altering_variant",
        "inframe_insertion,inframe_deletion,protein_altering_variant",
    )
    CANDIDATES = ("frameshift_variant", "stop_gained,frameshift_variant")

    @staticmethod
    def row(terms: str, location: str = "2:190970701-190970707") -> str:
        return stop_codon_row("-", terms, "961-967", "264-270", location)

    def taken(self, reference: str, candidate: str) -> bool:
        return cvo.stop_retained_three_base_stop_codon_test_pair(
            stop_codon_pair_row(self.row(reference)), stop_codon_pair_row(self.row(candidate)))

    def test_the_retained_term_is_taken_alone_or_beside_one_read_away_term_and_no_more(self) -> None:
        for candidate in self.CANDIDATES:
            self.assertTrue(self.taken("stop_retained_variant", candidate), candidate)
            for one in sorted(cvo.FRAMESHIFT_READ_TERMS):
                self.assertTrue(self.taken(f"{one},stop_retained_variant", candidate), (one, candidate))
        taken = [(several, candidate) for several in self.SEVERAL for candidate in self.CANDIDATES
                 if self.taken(f"{several},stop_retained_variant", candidate)]
        self.assertEqual(taken, [])

    def test_the_comparison_sets_none_or_one_aside_and_leaves_several_open(self) -> None:
        references = ["stop_retained_variant", *(f"{one},stop_retained_variant" for one in sorted(cvo.FRAMESHIFT_READ_TERMS)),
                      *(f"{several},stop_retained_variant" for several in self.SEVERAL)]
        perl = "".join(self.row(terms, f"2:{1000 * (i + 1)}-{1000 * (i + 1) + 6}") for i, terms in enumerate(references))
        rust = "".join(self.row("frameshift_variant", f"2:{1000 * (i + 1)}-{1000 * (i + 1) + 6}") for i in range(len(references)))
        adjusted, open_rows = compare_row_outputs(self, perl, rust)
        self.assertEqual(adjusted["excluded_categories"]["stop_retained_three_base_stop_codon_test"], {"count": 4, "excluded": 4})
        self.assertEqual(len(open_rows), 2 * len(self.SEVERAL))

    def test_the_definitions_name_the_term_alone_or_beside_one_and_the_row_under_an_identifier(self) -> None:
        definition = cvo.PAIR_RULE_DEFINITIONS["stop_retained_three_base_stop_codon_test"]
        self.assertIn(
            "where the reference writes it alone or beside one of inframe_insertion, inframe_deletion and "
            "protein_altering_variant, frameshift_variant with or without stop_gained in the place of those terms.",
            definition,
        )
        self.assertIn(
            "Under an identifier a row is read as an insertion's by a Location of two neighbouring positions and an "
            "Allele of another length than two",
            cvo.PAIR_RULE_DEFINITIONS["stop_codon_overlap_insertion_range"],
        )


def row_pair_predicate(test: unittest.TestCase, bucket: str) -> Callable[..., bool]:
    """The predicate of the row-aware pair rule the 116.2 registry holds under ``bucket``."""
    predicates = {rule.bucket: rule.predicate for rule in cvo.REFERENCE_REGISTRIES["116.2"].row_pair_rules}
    test.assertIn(bucket, predicates)
    return predicates[bucket]


def splice_term_rows(reference: str, candidate: str, *, allele: str = "GTTCA", location: str = "7:100998-101000",
                     cdna: str = "998-1000", cds: str = "898-900", candidate_cdna: str | None = None,
                     candidate_cds: str | None = None, **columns: str) -> tuple[str, str]:
    """The reference's and vep-rs's fabricated row of one key: the two sets of terms on three
    positions that lie in the transcript's sequence under an allele of five, unless an
    argument says otherwise. ``columns`` are `stop_codon_row`'s."""
    return (
        stop_codon_row(allele, reference, cdna, cds, location, **columns),
        stop_codon_row(allele, candidate, candidate_cdna or cdna, candidate_cds or cds, location, **columns),
    )


# What differs from `splice_term_rows`'s row on two neighbouring positions of one exon, and on
# the first three positions of an exon.
TWO_NEIGHBOURS = {"location": "7:100999-101000", "cdna": "999-1000", "cds": "899-900"}
EXON_START = {"cdna": "901-903", "cds": "801-803"}


class SpliceTermPastTheRecordSpanRowTests(unittest.TestCase):
    """Without its record, `splice_term_past_the_record_span` takes a Transcript row whose
    Location holds three or more positions under an allele of more nucleotide letters, both
    rows carrying the same cDNA and CDS positions and two cDNA positions as far apart as the
    Location's ends, where every term of the six the reference holds (the five intron-window
    splice terms and intron_variant) is gone from vep-rs's set and none of the six is in it,
    with nothing in their place or splice_region_variant where a donor, acceptor, donor
    region or fifth base term is gone. Every other row stays a difference, a Location of two
    neighbouring positions included. Every row stands under an identifier, so none carries
    its record's alleles."""

    BUCKET = "splice_term_past_the_record_span"
    DONOR = "splice_donor_variant,frameshift_variant"
    FIFTH = "splice_donor_variant,frameshift_variant,splice_donor_5th_base_variant"
    TRACT = "frameshift_variant,splice_polypyrimidine_tract_variant"
    INTRON = "frameshift_variant,intron_variant"
    BACK = "frameshift_variant,splice_region_variant"
    NOTHING = "frameshift_variant"
    # (the shape, the reference's terms, vep-rs's terms, what differs from `splice_term_rows`'s row)
    TAKEN = (
        ("the donor term gone and the region term back", DONOR, BACK, {}),
        ("the donor term gone and nothing back", DONOR, NOTHING, {}),
        ("the donor and the donor region term gone", DONOR + ",splice_donor_region_variant", BACK, {}),
        ("the donor and the fifth base term gone", FIFTH, NOTHING, {"allele": "GTTCAGGAT"}),
        ("the acceptor term gone at an exon's first positions", "splice_acceptor_variant,frameshift_variant", BACK, EXON_START),
        ("the acceptor and the tract term gone", "splice_acceptor_variant," + TRACT, BACK, EXON_START),
        ("the tract term alone gone and nothing back", TRACT, NOTHING, EXON_START),
        ("the intron term gone with the donor and the donor region term", DONOR + ",splice_donor_region_variant,intron_variant", BACK, {}),
        ("the intron term gone with the donor and the fifth base term", FIFTH + ",intron_variant", NOTHING, {"allele": "GTTCAGGAT"}),
        ("the intron term gone with the acceptor term", "splice_acceptor_variant," + INTRON, BACK, EXON_START),
        ("the intron term alone gone and nothing back", INTRON, NOTHING, {}),
        ("four positions under an allele of five", DONOR, BACK, {"location": "7:100997-101000", "cdna": "997-1000", "cds": "897-900"}),
        ("positions with the transcript's lengths appended", DONOR, BACK, {"cdna": "998-1000/2400", "cds": "898-900/1800"}),
    )
    REFUSED = (
        ("a span across an exon's end", DONOR, BACK, {"cdna": "999-?", "cds": "899-?"}),
        ("an intron span", "splice_donor_variant,intron_variant", "intron_variant", {"cdna": "-", "cds": "-"}),
        ("an allele as long as the span", DONOR, NOTHING, {"allele": "GTT"}),
        ("an allele shorter than the span", DONOR, NOTHING, {"allele": "GT"}),
        ("a Location of one position", DONOR, NOTHING, {"location": "7:101000", "cdna": "1000", "cds": "900", "allele": "GTT"}),
        ("a deletion", DONOR, NOTHING, {"allele": "-"}),
        ("an insertion between two neighbouring positions", DONOR, NOTHING, {**TWO_NEIGHBOURS, "allele": "GTT"}),
        ("two replaced bases", DONOR, BACK, TWO_NEIGHBOURS),
        ("one of the five splice terms kept", FIFTH, "frameshift_variant,splice_donor_5th_base_variant", {}),
        ("the intron term kept where the donor term is gone", DONOR + ",intron_variant", BACK + ",intron_variant", {}),
        ("the intron term kept and nothing back", DONOR + ",intron_variant", INTRON, {}),
        ("the donor region term come", FIFTH, "frameshift_variant,splice_donor_region_variant", {}),
        ("the region term come where the tract term alone is gone", TRACT, BACK, EXON_START),
        ("the region term come where the intron term alone is gone", INTRON, BACK, {}),
        ("other cDNA positions on vep-rs's row", DONOR, BACK, {"candidate_cdna": "498-500"}),
        ("other CDS positions on vep-rs's row", DONOR, BACK, {"candidate_cds": "398-400"}),
        ("cDNA positions further apart than the Location's ends", DONOR, BACK, {"cdna": "997-1000", "cds": "897-900"}),
        ("the terms the other way round", BACK, DONOR, {}),
        ("a further term gone", DONOR, "splice_region_variant", {}),
        ("a further term come", DONOR, "stop_gained," + BACK, {}),
        ("a row of another feature type", DONOR, BACK, {"feature_type": "RegulatoryFeature"}),
    )

    def taken(self, reference: str, candidate: str, changed: dict[str, str]) -> bool:
        e, a = (stop_codon_pair_row(row) for row in splice_term_rows(reference, candidate, **changed))
        return row_pair_predicate(self, self.BUCKET)(e, a)

    def test_a_longer_allele_over_three_positions_of_the_transcript_is_taken_and_no_other_row(self) -> None:
        self.assertEqual([shape for shape, *row in self.TAKEN if not self.taken(*row)], [])
        self.assertEqual([shape for shape, *row in self.REFUSED if self.taken(*row)], [])

    def test_the_comparison_sets_the_taken_rows_aside_and_leaves_the_others_open(self) -> None:
        perl = rust = ""
        for i, (_, reference, candidate, changed) in enumerate((*self.TAKEN, *self.REFUSED)):
            e, a = splice_term_rows(reference, candidate, feature=f"ENST{i:011d}", **changed)
            perl, rust = perl + e, rust + a
        adjusted, open_rows = compare_row_outputs(self, perl, rust)
        self.assertEqual(adjusted["excluded_categories"][self.BUCKET], {"count": len(self.TAKEN), "excluded": len(self.TAKEN)})
        self.assertEqual(sorted(row.split("\t")[4] for row in open_rows),
                         sorted(2 * [f"ENST{i:011d}" for i in range(len(self.TAKEN), len(self.TAKEN) + len(self.REFUSED))]))


class RecordAlleleReferenceTests(unittest.TestCase):
    """`replaced_reference` reads, by position from a record's own alleles, the reference bases
    an output row's allele stands in for; `uploaded_position` reads the record's POS from a
    first column that spells the alleles; and `pair_row` carries the bases of such a row.
    Every record is fabricated."""

    # (the shape, POS, REF, the alternates, Location, Allele, the reference bases)
    READ = (
        ("three bases replaced by five", 100998, "ACG", ("GTTCA",), "7:100998-101000", "GTTCA", "ACG"),
        ("three bases replaced by five behind a shared first base", 100997, "TACG", ("TGTTCA",), "7:100998-101000", "GTTCA", "ACG"),
        ("one base replaced by one", 101000, "A", ("G",), "7:101000", "G", "A"),
        ("the second of two alternates", 100998, "ACG", ("TTT", "GTTCA"), "7:100998-101000", "GTTCA", "ACG"),
        ("a pure insertion", 100999, "A", ("AGTT",), "7:100999-101000", "GTT", ""),
        ("a pure insertion behind two shared bases", 100998, "AC", ("ACGTT",), "7:100999-101000", "GTT", ""),
        ("a deletion", 100997, "TACG", ("T",), "7:100998-101000", "-", "ACG"),
        ("lower-case alleles", 100998, "acg", ("gttca",), "7:100998-101000", "GTTCA", "ACG"),
    )
    # (the shape, POS, REF, the alternates, Location, Allele)
    NOT_READ = (
        ("an alternate that writes no such row", 100998, "ACG", ("GTTCAT",), "7:100998-101000", "GTTCA"),
        ("a Location past the record's bases", 100998, "ACG", ("GTTCA",), "7:100999-101001", "GTTCA"),
        ("two alternates that read the row differently", 100999, "AC", ("GTT", "AGTTC"), "7:100999-101000", "GTT"),
        ("an allele that is not a sequence", 100998, "ACG", ("<DEL>",), "7:100998-101000", "deletion"),
        ("a Location that names no whole positions", 100998, "ACG", ("GTTCA",), "7:100998-101000.5", "GTTCA"),
        ("a row without a Location", 100998, "ACG", ("GTTCA",), None, "GTTCA"),
    )

    def test_the_bases_an_allele_stands_in_for_are_read_by_position_from_the_record_s_alleles(self) -> None:
        read = {shape: cvo.replaced_reference(*record) for shape, *record, _ in self.READ}
        self.assertEqual(read, {shape: bases for shape, *_, bases in self.READ})
        self.assertEqual([shape for shape, *record in self.NOT_READ if cvo.replaced_reference(*record) is not None], [])

    def test_a_first_column_that_spells_the_alleles_spells_the_record_s_position_before_them(self) -> None:
        self.assertEqual(cvo.uploaded_position("7_100998_ACG/GTTCA"), 100998)
        self.assertEqual(cvo.uploaded_position("7_100997_tacg/t/-"), 100997)
        self.assertEqual(cvo.uploaded_position("chrUn_made_up_100998_A/T"), 100998)
        for identifier in ("1993671", "rs1993671", "7_100998", "made_up_7"):
            self.assertIsNone(cvo.uploaded_position(identifier), identifier)

    def test_a_row_carries_the_bases_its_first_column_spells_and_none_under_an_identifier(self) -> None:
        self.assertIn("reference", [f.name for f in dataclasses.fields(cvo.PairRow)])

        def row(name: str, location: str, allele: str) -> cvo.PairRow:
            return stop_codon_pair_row(stop_codon_row(allele, "frameshift_variant", "998-1000", "898-900", location, name=name))

        for shape, position, reference, alternates, location, allele, bases in self.READ:
            spelt = f"7_{position}_{reference}/{'/'.join(alternates)}"
            self.assertEqual(row(spelt, location, allele).reference, bases, shape)
            self.assertIsNone(row("1993671", location, allele).reference, shape)
        self.assertIsNone(row("7_100998_ACG/GTTCAT", "7:100998-101000", "GTTCA").reference, "an alternate that writes no such row")


class RecordReferencesFromTheInputTests(unittest.TestCase):
    """`read_record_references` gives each row that carries no reference bases the bases its
    allele stands in for, read by position from the records of the input VCF, plain or
    gzipped. Every record is fabricated."""

    HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    # (CHROM, POS, REF, ALT)
    RECORDS = (
        ("7", 100998, "ACG", "GTTCA"),
        ("chr7", 200997, "TACG", "TGTTCA"),
        ("7", 300999, "A", "AGTT"),
        ("7", 400997, "TACG", "T"),
        ("7", 500999, "AC", "GTT"),
        ("7", 500999, "A", "AGTT"),
        ("7", 700998, "ACG", "TTT,GTTCA"),
    )
    # (the shape, Location, Allele, the reference bases the records give)
    ROWS = (
        ("a replacement", "7:100998-101000", "GTTCA", "ACG"),
        ("a replacement behind a shared first base, the record's chromosome under a prefix", "7:200998-201000", "GTTCA", "ACG"),
        ("an insertion, the row's chromosome under a prefix", "chr7:300999-301000", "GTT", ""),
        ("a deletion", "7:400998-401000", "-", "ACG"),
        ("the second alternate of a record", "7:700998-701000", "GTTCA", "ACG"),
        ("a row two records read differently", "7:500999-501000", "GTT", None),
        ("a row whose record lies elsewhere", "7:600998-601000", "GTTCA", None),
        ("a row at a record's place that the record does not write", "7:100998-101000", "GTTC", None),
        ("a row of another chromosome than its record's", "8:100998-101000", "GTTCA", None),
    )

    def rows(self) -> list[cvo.PairRow]:
        return [stop_codon_pair_row(stop_codon_row(allele, "frameshift_variant", "998-1000", "898-900", location))
                for _, location, allele, _ in self.ROWS]

    def written(self, name: str) -> Path:
        """The records as an input VCF of the given name, gzipped under a ``.gz`` suffix."""
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, True)
        opener = gzip.open if name.endswith(".gz") else open
        with opener(tmp / name, "wt", encoding="utf-8") as handle:
            handle.write(self.HEADER)
            for i, (chrom, position, reference, alternates) in enumerate(self.RECORDS):
                handle.write(f"{chrom}\t{position}\tmade_up_{i}\t{reference}\t{alternates}\t.\t.\t.\n")
        return tmp / name

    def test_rows_under_an_identifier_take_their_record_s_bases_from_a_plain_and_a_gzipped_vcf(self) -> None:
        for name in ("input.vcf", "input.vcf.gz"):
            rows = self.rows()
            self.assertEqual([row.reference for row in rows], [None] * len(self.ROWS), name)
            cvo.read_record_references(self.written(name), rows)
            self.assertEqual({shape: row.reference for (shape, *_), row in zip(self.ROWS, rows)},
                             {shape: bases for shape, *_, bases in self.ROWS}, name)

    def test_a_row_that_carries_its_bases_keeps_them_and_the_file_is_not_opened_when_no_row_asks(self) -> None:
        records = self.written("input.vcf")
        absent = records.parent / "absent.vcf"
        spelt = stop_codon_pair_row(stop_codon_row("GTT", "frameshift_variant", "999-1000", "899-900", "7:500999-501000",
                                                   name="7_500999_AC/GTT"))
        self.assertEqual(spelt.reference, "AC")
        for path in (absent, records):
            cvo.read_record_references(path, [spelt])
            self.assertEqual(spelt.reference, "AC", path.name)
        cvo.read_record_references(absent, [])
        with self.assertRaises(FileNotFoundError):
            cvo.read_record_references(absent, self.rows())


class SpliceTermRecordRegionTests(unittest.TestCase):
    """The three readings `splice_term_past_the_record_span` makes of a record: the regions in
    which its two alleles differ (`differing_regions`), the terms those regions and the span
    give beside one intron (`intron_rows`), and the places of that intron a row's cDNA
    positions allow (`intron_placements`). Every offset counts from the span's first base."""

    DONOR, ACCEPTOR, FIFTH = "splice_donor_variant", "splice_acceptor_variant", "splice_donor_5th_base_variant"
    DONOR_REGION, TRACT = "splice_donor_region_variant", "splice_polypyrimidine_tract_variant"
    REGION, INTRON = "splice_region_variant", "intron_variant"
    # (the shape, the reference bases, the allele, the regions)
    DIFFERING = (
        ("every index differs and the longer allele runs on", "ACG", "GTTCA", [(0, 4)]),
        ("two indices alike before the longer allele runs on", "ACGT", "ACTTAG", [(2, 2), (4, 5)]),
        ("three runs of indices", "ACGTAC", "ATGTCCGG", [(1, 1), (4, 4), (6, 7)]),
        ("a shorter allele", "ACGT", "AT", [(1, 3)]),
        ("one reference base", "A", "GTT", [(0, 0)]),
        ("an allele of one base", "ACG", "T", [(0, 2)]),
        ("two alleles alike", "ACG", "ACG", [(0, 2)]),
    )
    # (the shape, the regions, the span's bases, the intron's first and last offset, whether the strand is the
    # forward one, whether another intron of the transcript is a frameshift intron, the reference's terms, the
    # candidate's terms)
    INTRON_ROWS = (
        ("an intron right after a span of three",
         ((0, 5),), 3, (3, 102), True, False, (DONOR, DONOR_REGION, INTRON), (REGION,)),
        ("the same intron on the reverse strand",
         ((0, 5),), 3, (3, 102), False, False, (ACCEPTOR, INTRON), (REGION,)),
        ("regions that end on the intron's second base, on the reverse strand",
         ((0, 4),), 3, (3, 102), False, False, (ACCEPTOR,), (REGION,)),
        ("an intron three bases after a span of two",
         ((0, 5),), 2, (4, 103), True, False, (DONOR,), (REGION,)),
        ("an intron four bases after a span of two",
         ((0, 5),), 2, (5, 104), True, False, (), ()),
        ("the donor region term withheld beside the fifth base term",
         ((0, 7),), 2, (2, 101), True, False, (DONOR, FIFTH, INTRON), (REGION,)),
        ("a span on an intron's third and fourth base",
         ((0, 3),), 2, (-2, 97), True, False, (FIFTH, INTRON), (DONOR_REGION, INTRON)),
        ("a span across an exon's end",
         ((0, 6),), 2, (1, 100), True, False, (DONOR, FIFTH, INTRON), (DONOR, INTRON)),
        ("a span of two before an intron's 5' end on the reverse strand",
         ((0, 3),), 2, (-96, 3), False, False, (DONOR, DONOR_REGION, INTRON), (DONOR_REGION, INTRON)),
        ("the tract off the exons",
         ((0, 14),), 3, (-69, 30), True, False, (TRACT, INTRON), (INTRON,)),
        ("the same span under the wide exons, twelve bases and more from each",
         ((0, 14),), 3, (-69, 30), True, True, (TRACT, INTRON), (INTRON,)),
        ("a span that ends four bases before an intron's end",
         ((0, 6),), 3, (-93, 6), True, False, (ACCEPTOR, TRACT, INTRON), (REGION, TRACT, INTRON)),
        ("the same span under the wide exons, which hold it",
         ((0, 6),), 3, (-93, 6), True, True, (ACCEPTOR, INTRON), (REGION, INTRON)),
        ("a short intron whose two ends the span reads",
         ((0, 23),), 3, (-2, 21), True, False, (ACCEPTOR, FIFTH, TRACT, INTRON), (FIFTH, INTRON)),
        ("a frameshift intron, whose regions set nothing",
         ((0, 5),), 3, (3, 10), True, False, (), ()),
        ("an intron of fourteen bases",
         ((0, 5),), 3, (3, 16), True, False, (DONOR, DONOR_REGION, INTRON), (REGION,)),
        ("a region before a frameshift intron and one on it",
         ((0, 2), (5, 6)), 3, (5, 12), True, False, (REGION,), (REGION,)),
        ("the same two regions beside a long intron",
         ((0, 2), (5, 6)), 3, (5, 104), True, False, (DONOR,), (REGION,)),
    )
    # Cases of the shape of `INTRON_ROWS` on which splice_region_variant follows the order of the regions: the
    # reference's row holds it where the last region read before one meets the intron's first two or last two
    # bases lies in a window of the term, the candidate's where any region does.
    REGION_TERM = (
        ("two runs, the first in a window of the region term and the last in none",
         ((5, 5), (8, 10)), 8, (-1, 98), True, False, (INTRON,), (REGION, INTRON)),
        ("two runs, the last in a window of the region term and the first in none",
         ((2, 2), (8, 10)), 8, (-85, 14), True, False, (REGION, TRACT, INTRON), (REGION, TRACT, INTRON)),
        ("five runs in a twenty-base intron, the last in a window",
         ((0, 0), (2, 2), (4, 4), (6, 6), (8, 10)), 8, (-6, 13), True, False, (REGION, TRACT, INTRON), (REGION, TRACT, INTRON)),
        ("the same five runs before the tract of a longer intron, the last in no window",
         ((0, 0), (2, 2), (4, 4), (6, 6), (8, 10)), 8, (-6, 24), True, False, (TRACT, INTRON), (REGION, INTRON)),
        ("a run on the intron's first two bases between two runs in a window, the donor term withholding the region term",
         ((0, 0), (2, 2), (6, 8)), 6, (2, 101), True, False, (DONOR, FIFTH, INTRON), (DONOR, DONOR_REGION, INTRON)),
    )

    def test_two_alleles_differ_in_runs_of_indices_and_one_base_gives_the_span(self) -> None:
        regions = {shape: cvo.differing_regions(reference, allele) for shape, reference, allele, _ in self.DIFFERING}
        self.assertEqual(regions, {shape: expected for shape, *_, expected in self.DIFFERING})

    def test_the_regions_give_the_reference_s_terms_and_the_span_the_candidate_s_beside_one_intron(self) -> None:
        rows = {shape: cvo.intron_rows(regions, span, *intron, forward, wide_exons)
                for shape, regions, span, intron, forward, wide_exons, _, _ in self.INTRON_ROWS}
        self.assertEqual(rows, {shape: (frozenset(expected), frozenset(actual)) for shape, *_, expected, actual in self.INTRON_ROWS})

    def test_the_last_region_read_decides_the_reference_s_region_term_and_any_region_the_candidate_s(self) -> None:
        rows = {shape: cvo.intron_rows(regions, span, *intron, forward, wide_exons)
                for shape, regions, span, intron, forward, wide_exons, _, _ in self.REGION_TERM}
        self.assertEqual(rows, {shape: (frozenset(expected), frozenset(actual)) for shape, *_, expected, actual in self.REGION_TERM})

    def test_two_numbered_cdna_positions_as_far_apart_as_the_span_put_the_intron_after_the_span_or_before_it(self) -> None:
        places = list(cvo.intron_placements(3, 2, ("998", "1000", "898", "900")))
        after = [place for place in places if place[1] > 0]
        self.assertEqual({(forward, start) for forward, start, _ in after},
                         {(forward, start) for forward in (True, False) for start in range(3, 10)}, "one to seven bases on")
        self.assertTrue(all(end >= start for _, start, end in after))
        self.assertIn((True, 3, 3), after, "an intron of one base")
        self.assertIn((False, 9, 10), after, "a two-base intron seven bases on")
        self.assertTrue(any(end - start > 10**6 for _, start, end in after), "an intron longer than any window")
        before = [place for place in places if place[1] <= 0]
        self.assertEqual(sorted(before), sorted((forward, start, end) for forward in (True, False)
                                                for start in range(-16, 0) for end in range(start, 0)),
                         "an intron of sixteen bases or fewer that ends before the span's first base")
        self.assertIn((True, -2, -1), before, "a two-base intron that ends beside the span")
        self.assertEqual(list(cvo.intron_placements(3, 2, ("997", "1000", "897", "900"))), [], "numbers further apart than the span")

    def test_two_numbered_cdna_positions_closer_than_the_span_put_an_intron_of_the_missing_length_inside_it(self) -> None:
        self.assertEqual(sorted(cvo.intron_placements(5, 2, ("998", "1000", "898", "900"))),
                         [(False, 1, 2), (False, 2, 3), (True, 1, 2), (True, 2, 3)], "two bases of five, the span's ends kept")
        self.assertEqual(sorted(cvo.intron_placements(3, 2, ("999", "1000", "899", "900"))),
                         [(False, 1, 1), (True, 1, 1)], "the middle base of three")
        self.assertEqual(list(cvo.intron_placements(2, 2, ("1000", "1000", "900", "900"))), [], "no base inside a span of two")

    def test_a_half_numbered_cdna_position_puts_one_end_of_the_intron_inside_the_span(self) -> None:
        places = set(cvo.intron_placements(2, 5, ("1000", "?", "900", "?")))
        self.assertEqual({start for forward, start, _ in places if forward}, {1})
        self.assertEqual({end for forward, _, end in places if not forward}, {0})
        self.assertTrue(all(end >= 1 if forward else start <= 0 for forward, start, end in places))
        mirrored = set(cvo.intron_placements(2, 5, ("?", "901", "?", "801")))
        self.assertEqual({end for forward, _, end in mirrored if forward}, {0})
        self.assertEqual({start for forward, start, _ in mirrored if not forward}, {1})
        three = cvo.intron_placements(3, 5, ("999", "?", "899", "?"))
        self.assertEqual({start for forward, start, _ in three if forward}, {1, 2})

    def test_an_unnumbered_cdna_position_puts_the_span_inside_the_intron(self) -> None:
        places = list(cvo.intron_placements(2, 3, ("-", "-", "-", "-")))
        self.assertEqual({forward for forward, _, _ in places}, {True, False})
        self.assertTrue(all(start <= 0 and end >= 1 for _, start, end in places))
        self.assertIn((True, 0, 1), places, "an intron that is the span")
        self.assertIn((False, -17, 42), places, "the furthest start and end a window reaches")
        self.assertNotIn((True, 0, 43), places)
        self.assertEqual(sum(forward and start == 0 and end > 10**6 for forward, start, end in places), 1, "a far end, once")
        reached = {end for _, _, end in cvo.intron_placements(2, 30, ("-", "-", "-", "-"))}
        self.assertIn(1 + 30 + 18, reached, "an end the regions reach")
        self.assertNotIn(1 + 30 + 19, reached)


@dataclasses.dataclass(frozen=True)
class RecordPair:
    """One fabricated pair of rows with its record and the two verdicts of
    `splice_term_past_the_record_span`: with the record's alleles, and without them."""

    kind: str
    position: int
    reference: str
    alternate: str
    location: str
    perl: str
    rust: str
    with_record: bool
    without_record: bool


def record_row(position: int, reference: str, alternate: str) -> tuple[str, str]:
    """The positions of the Location and the Allele of the row a VCF record of one alternate
    writes: the shared first base dropped where the two alleles differ in length, as
    `reference_allele_is_empty` reads a record, an insertion between its two neighbouring
    positions and a deletion under ``-``."""
    if len(reference) != len(alternate) and reference[0] == alternate[0]:
        position, reference, alternate = position + 1, reference[1:], alternate[1:]
    if not reference:
        return f"{position - 1}-{position}", alternate
    end = position + len(reference) - 1
    return (f"{position}-{end}" if end > position else str(position)), alternate or "-"


class SpliceTermPastTheRecordSpanKindTests(unittest.TestCase):
    """One fabricated pair of rows with its record for each kind of record on which the
    reference's set and vep-rs's differ by a term of the six, or by the splice_region_variant
    an insertion beside an exon's end loses, with the two verdicts
    `splice_term_past_the_record_span` gives it: with the record's alleles, and without them.

    Each pair stands in its own block of 10,000 positions of chromosome 7, on one of two
    layouts of the block. On the forward strand's an exon ends at position 1000 (cDNA 1000,
    CDS 900), its intron runs from 1001 to 1100 and the next exon starts at 1101 (cDNA 1001).
    On the reverse strand's an exon starts at 1000 (cDNA 901, CDS 801) and runs down the
    block, the intron before it runs from its 3' end at 1001 to its 5' end at 1100, and the
    exon before that ends at 1101 (cDNA 900). A two-base intron runs from 1001 to 1002 and a
    twenty-base intron from 1001 to 1020, each in the place of the layout's intron.

    With the record a pair is set aside where an allele longer than its two or more reference
    bases gives, at one place of the intron, the reference's terms from its regions and
    vep-rs's from its span, which lies inside an exon, across its end or in the intron.
    Without the record only three or more positions of the transcript's sequence are set
    aside, vep-rs's row holding none of the six. A pure insertion, an allele that is not the
    longer one, a term gone whose window holds a base of the span, a differing term outside
    the class and differing positions stay differences either way."""

    BUCKET = "splice_term_past_the_record_span"
    DONOR, ACCEPTOR, FIFTH = "splice_donor_variant", "splice_acceptor_variant", "splice_donor_5th_base_variant"
    DONOR_REGION, TRACT = "splice_donor_region_variant", "splice_polypyrimidine_tract_variant"
    REGION, SHIFT, INTRON, CODING = "splice_region_variant", "frameshift_variant", "intron_variant", "coding_sequence_variant"
    # The cDNA and CDS positions of two and of three bases: the last of the forward layout's exon, the first of
    # the reverse layout's, and those that end two bases short of the intron on each layout; and the positions
    # of a span that holds no base of the transcript's sequence.
    END2, END3, START2, START3 = ("999-1000", "899-900"), ("998-1000", "898-900"), ("901-902", "801-802"), ("901-903", "801-803")
    SHORT2, SHORT3 = ("997-998", "897-898"), ("996-998", "896-898")
    REVERSE_SHORT2, REVERSE_SHORT3, NONE = ("903-904", "803-804"), ("903-905", "803-805"), ("-", "-")
    # The reference's terms on a span before the forward layout's intron whose regions reach its third base, and
    # those whose regions reach its fifth; the terms vep-rs writes with the region term back, and with nothing.
    THIRD, PAST_FIFTH = (DONOR, SHIFT, DONOR_REGION, INTRON), (DONOR, SHIFT, FIFTH, INTRON)
    BACK, NOTHING = (SHIFT, REGION), (SHIFT,)
    # (the kind, the record's POS in its block, REF, ALT, the cDNA and CDS positions with vep-rs's cDNA position
    # where it is another, the reference's terms, vep-rs's terms, set aside with the record, set aside without it)
    PAIRS = (
        ("two bases at an exon's end replaced by five, the region term back", 999, "AC", "GTTCA", END2, THIRD, BACK, True, False),
        ("two bases at an exon's end replaced by five, nothing back", 999, "AC", "GTTCA", END2, THIRD, NOTHING, False, False),
        ("three bases at an exon's end replaced by six, the region term back", 998, "ACG", "GTTCAG", END3, THIRD, BACK, True, True),
        ("three bases at an exon's end replaced by six, nothing back", 998, "ACG", "GTTCAG", END3, THIRD, NOTHING, False, True),
        ("two bases at an exon's end replaced by six, the region term back", 999, "AC", "GTTCAG", END2, THIRD, BACK, True, False),
        ("two bases at an exon's end replaced by six, nothing back", 999, "AC", "GTTCAG", END2, THIRD, NOTHING, False, False),
        ("three bases at an exon's end replaced by seven, the region term back", 998, "ACG", "GTTCAGG", END3, THIRD, BACK, True, True),
        ("three bases at an exon's end replaced by seven, nothing back", 998, "ACG", "GTTCAGG", END3, THIRD, NOTHING, False, True),
        ("two bases at an exon's end replaced by eight, the fifth base term gone with the donor's",
         999, "AC", "GTTCAGGA", END2, PAST_FIFTH, BACK, True, False),
        ("two bases at an exon's end replaced by eight, nothing back", 999, "AC", "GTTCAGGA", END2, PAST_FIFTH, NOTHING, False, False),
        ("three bases at an exon's end replaced by nine, the fifth base term gone with the donor's",
         998, "ACG", "GTTCAGGAT", END3, PAST_FIFTH, BACK, True, True),
        ("three bases at an exon's end replaced by nine, nothing back", 998, "ACG", "GTTCAGGAT", END3, PAST_FIFTH, NOTHING, False, True),
        ("two bases one short of an exon's end replaced by six", 998, "AC", "GTTCAG", ("998-999", "898-899"), THIRD, BACK, True, False),
        ("two bases two short of an exon's end replaced by seven, the region term back",
         997, "AC", "GTTCAGG", SHORT2, THIRD, BACK, True, False),
        ("two bases two short of an exon's end replaced by seven, nothing back", 997, "AC", "GTTCAGG", SHORT2, THIRD, NOTHING, False, False),
        ("three bases two short of an exon's end replaced by eight, the region term back",
         996, "ACG", "GTTCAGGA", SHORT3, THIRD, BACK, True, True),
        ("three bases two short of an exon's end replaced by eight, nothing back",
         996, "ACG", "GTTCAGGA", SHORT3, THIRD, NOTHING, False, True),
        ("two bases two short of an exon's end replaced by eight, the region term back",
         997, "AC", "GTTCAGGA", SHORT2, THIRD, BACK, True, False),
        ("two bases two short of an exon's end replaced by eight, nothing back", 997, "AC", "GTTCAGGA", SHORT2, THIRD, NOTHING, False, False),
        ("two bases at an exon's start on the reverse strand replaced by four, the region term back",
         999, "AC", "GTTC", START2, (ACCEPTOR, SHIFT), BACK, True, False),
        ("two bases at an exon's start on the reverse strand replaced by four, nothing back",
         999, "AC", "GTTC", START2, (ACCEPTOR, SHIFT), NOTHING, False, False),
        ("three bases at an exon's start on the reverse strand replaced by five, the region term back",
         998, "ACG", "GTTCA", START3, (ACCEPTOR, SHIFT), BACK, True, True),
        ("three bases at an exon's start on the reverse strand replaced by five, nothing back",
         998, "ACG", "GTTCA", START3, (ACCEPTOR, SHIFT), NOTHING, False, True),
        ("two bases two short of an exon's start on the reverse strand replaced by six, the region term back",
         997, "AC", "GTTCAG", REVERSE_SHORT2, (ACCEPTOR, SHIFT), BACK, True, False),
        ("two bases two short of an exon's start on the reverse strand replaced by six, nothing back",
         997, "AC", "GTTCAG", REVERSE_SHORT2, (ACCEPTOR, SHIFT), NOTHING, False, False),
        ("three bases two short of an exon's start on the reverse strand replaced by seven, the region term back",
         996, "ACG", "GTTCAGG", REVERSE_SHORT3, (ACCEPTOR, SHIFT), BACK, True, True),
        ("three bases two short of an exon's start on the reverse strand replaced by seven, nothing back",
         996, "ACG", "GTTCAGG", REVERSE_SHORT3, (ACCEPTOR, SHIFT), NOTHING, False, True),
        ("three bases two short of an exon's start on the reverse strand replaced by nine, the intron term gone with the acceptor's",
         996, "ACG", "GTTCAGGAT", REVERSE_SHORT3, (ACCEPTOR, SHIFT, INTRON), BACK, True, True),
        ("three bases two short of an exon's start on the reverse strand replaced by nine, nothing back",
         996, "ACG", "GTTCAGGAT", REVERSE_SHORT3, (ACCEPTOR, SHIFT, INTRON), NOTHING, False, True),
        ("the two exon bases after a two-base intron replaced by six, the fifth base term gone and the donor region term back",
         1003, "AC", "GTTCAG", ("1001-1002", "901-902"), (FIFTH, SHIFT), (DONOR_REGION, SHIFT), True, False),
        ("three bases whose last lies six before a two-base intron on the reverse strand replaced by six, the fifth base term gone",
         993, "ACG", "GTTCAG", ("906-908", "806-808"), (FIFTH, SHIFT), BACK, True, True),
        ("two bases of an intron on the reverse strand replaced by four, the donor term gone beside a kept donor region term",
         1097, "AC", "GTTC", NONE, (DONOR, DONOR_REGION, INTRON), (DONOR_REGION, INTRON), True, False),
        ("two bases of an intron on the reverse strand replaced by eight, the donor and the fifth base term gone",
         1093, "AC", "GTTCAGGA", NONE, (DONOR, FIFTH, INTRON), (REGION, INTRON), True, False),
        ("two bases of an intron on the reverse strand replaced by eight, nothing back",
         1093, "AC", "GTTCAGGA", NONE, (DONOR, FIFTH, INTRON), (INTRON,), False, False),
        ("three bases of an intron replaced by seven, the acceptor term gone beside a kept tract term",
         1094, "ACG", "GTTCAGG", NONE, (ACCEPTOR, TRACT, INTRON), (REGION, TRACT, INTRON), True, False),
        ("three bases of an intron replaced by seven, nothing back",
         1094, "ACG", "GTTCAGG", NONE, (ACCEPTOR, TRACT, INTRON), (TRACT, INTRON), False, False),
        ("three bases of an intron replaced by seven on a transcript with a frameshift intron, the acceptor term gone and no tract term",
         1094, "ACG", "GTTCAGG", NONE, (ACCEPTOR, INTRON), (REGION, INTRON), True, False),
        ("an intron's third and fourth base replaced by four, the fifth base term gone and the donor region term back",
         1003, "AC", "GTTC", NONE, (FIFTH, INTRON), (DONOR_REGION, INTRON), True, False),
        ("an intron's third and fourth base replaced by four, nothing back",
         1003, "AC", "GTTC", NONE, (FIFTH, INTRON), (INTRON,), False, False),
        ("two bases across an exon's end replaced by seven, the fifth base term gone beside a kept donor term",
         1000, "AC", "GTTCAGG", ("1000-?", "900-?"), (DONOR, FIFTH, CODING, INTRON), (DONOR, CODING, INTRON), True, False),
        ("an intron's second and third base replaced by five, the fifth base term gone and the donor region term back",
         1002, "AC", "GTTCA", NONE, (DONOR, FIFTH, INTRON), (DONOR, DONOR_REGION, INTRON), True, False),
        ("the same two rows without the donor region term, the rows of an intron's first two bases",
         1002, "AC", "GTTCA", NONE, (DONOR, FIFTH, INTRON), (DONOR, INTRON), True, False),
        ("two bases of an intron on the reverse strand replaced by four, the fifth base term gone and the donor region term back",
         1094, "AC", "GTTC", NONE, (FIFTH, INTRON), (DONOR_REGION, INTRON), True, False),
        ("two bases of an intron on the reverse strand replaced by four, the fifth base term gone and nothing back",
         1094, "AC", "GTTC", NONE, (FIFTH, INTRON), (INTRON,), False, False),
        ("two bases of an intron on the reverse strand replaced by five, the fifth base term gone and the region term back",
         1093, "AC", "GTTCA", NONE, (FIFTH, INTRON), (REGION, INTRON), True, False),
        ("two bases of an intron on the reverse strand replaced by five, the fifth base term gone and nothing back",
         1093, "AC", "GTTCA", NONE, (FIFTH, INTRON), (INTRON,), False, False),
        ("an intron's first two bases replaced by four, the donor region term gone beside a kept donor term",
         1001, "AC", "GTTC", NONE, (DONOR, DONOR_REGION, INTRON), (DONOR, INTRON), True, False),
        ("two bases across an exon's end replaced by five, the donor region term gone beside a kept donor term",
         1000, "AC", "GTTCA", ("1000-?", "900-?"), (DONOR, DONOR_REGION, CODING, INTRON), (DONOR, CODING, INTRON), True, False),
        ("two bases of an intron on the reverse strand replaced by four, the donor region term gone and the region term back",
         1092, "AC", "GTTC", NONE, (DONOR_REGION, INTRON), (REGION, INTRON), True, False),
        ("two bases of an intron on the reverse strand replaced by four, the donor region term gone and nothing back",
         1092, "AC", "GTTC", NONE, (DONOR_REGION, INTRON), (INTRON,), False, False),
        ("three bases of an intron replaced by fifteen, the tract term gone",
         1070, "ACG", "GTTCAGGATCCTTAG", NONE, (TRACT, INTRON), (INTRON,), True, False),
        ("an intron's last two bases on the reverse strand replaced by five, the tract term gone beside a kept acceptor term",
         1001, "AC", "GTTCA", NONE, (ACCEPTOR, TRACT, INTRON), (ACCEPTOR, INTRON), True, False),
        ("three bases across an exon's end replaced by eight, the intron term kept on the intron's first base",
         999, "ACG", "GTTCAGGA", ("999-?", "899-?"), (DONOR, FIFTH, CODING, INTRON), (DONOR, CODING, INTRON), True, False),
        ("three bases across an exon's end replaced by eight, the intron term gone from the intron's first base",
         999, "ACG", "GTTCAGGA", ("999-?", "899-?"), (DONOR, FIFTH, CODING, INTRON), (DONOR, CODING), False, False),
        ("two bases across an exon's end replaced by four with one base alike, the donor region term gone",
         1000, "AC", "TCTG", ("1000-?", "900-?"), (DONOR, DONOR_REGION, CODING, INTRON), (DONOR, CODING, INTRON), True, False),
        ("two bases across an exon's end replaced by four with one base alike, the donor term gone from the donor's base",
         1000, "AC", "TCTG", ("1000-?", "900-?"), (DONOR, DONOR_REGION, CODING, INTRON), (CODING, INTRON), False, False),
        ("eight bases of an intron on the reverse strand replaced by ten, the donor term gone from the donor's base",
         1092, "ACGTACGT", "GTGTATCTAG", NONE, (DONOR, DONOR_REGION, INTRON), (DONOR_REGION, INTRON), False, False),
        ("three bases of an intron on the reverse strand replaced by four, the tract term gone from the tract's base",
         1001, "ACG", "TCGA", NONE, (ACCEPTOR, TRACT, INTRON), (ACCEPTOR, INTRON), False, False),
        ("two bases across an exon's start on the reverse strand replaced by eight, the acceptor term gone from the acceptor's base",
         1000, "AC", "TCGTTCAG", ("?-901", "?-801"), (ACCEPTOR, CODING, INTRON), (CODING, INTRON), False, False),
        ("eight bases of a twenty-base intron replaced by eleven alike at four of them, the tract term gone from the tract's bases",
         1007, "ACGTACGT", "GCTTGCTTGAT", NONE, (REGION, TRACT, INTRON), (REGION, INTRON), False, False),
        ("an insertion of three between an intron's second and third base",
         1002, "A", "AGTT", NONE, (REGION, INTRON), (INTRON,), False, False),
        ("an insertion of three between the third and the second base from an intron's 3' end, the region term gone",
         1002, "T", "TGTT", NONE, (REGION, TRACT, INTRON), (TRACT, INTRON), False, False),
        ("an insertion of three between the third and the second base from an intron's 3' end, the tract term gone",
         1002, "T", "TGTT", NONE, (REGION, TRACT, INTRON), (REGION, INTRON), False, False),
        ("an insertion of four between an intron's last base and an exon's first",
         1100, "A", "AGTTC", ("?-1001", "?-901"), (REGION, CODING, INTRON), (CODING, INTRON), False, False),
        ("an insertion of four between an exon's last base and an intron's first on the reverse strand",
         1100, "A", "AGTTC", ("900-?", "800-?"), (REGION, CODING, INTRON), (CODING, INTRON), False, False),
        ("an insertion of three between a two-base intron's last base and an exon's first",
         1002, "A", "AGTT", ("?-1001", "?-901"), (REGION, CODING, INTRON), (CODING, INTRON), False, False),
        ("an insertion of three between an exon's last base and a two-base intron's first on the reverse strand",
         1002, "A", "AGTT", ("900-?", "800-?"), (REGION, CODING, INTRON), (CODING, INTRON), False, False),
        ("an insertion of five between the eighteenth and the seventeenth base from an intron's 3' end",
         1083, "A", "AGTTCA", NONE, (TRACT, INTRON), (INTRON,), False, False),
        ("an insertion of five between the seventeenth and the eighteenth base from an intron's 3' end on the reverse strand",
         1017, "A", "AGTTCA", NONE, (TRACT, INTRON), (INTRON,), False, False),
        ("an insertion of three between two bases of one exon",
         990, "A", "AGTT", ("990-991", "890-891"), (DONOR, SHIFT), NOTHING, False, False),
        ("an insertion of four between an exon's last two bases",
         999, "A", "AGTTC", END2, (DONOR, SHIFT, REGION), BACK, False, False),
        ("three bases replaced by three", 998, "ACG", "GTT", END3, (DONOR, SHIFT), NOTHING, False, False),
        ("four bases replaced by three", 997, "ACGA", "GTT", ("997-1000", "897-900"), (DONOR, SHIFT), NOTHING, False, False),
        ("three bases deleted", 997, "TACG", "T", END3, (DONOR, SHIFT), NOTHING, False, False),
        ("one base replaced by three", 1000, "A", "GTT", ("1000", "900"), (DONOR, SHIFT), NOTHING, False, False),
        ("a further term gone", 998, "ACG", "GTTCAG", END3, THIRD, (REGION,), False, False),
        ("a further term come", 998, "ACG", "GTTCAG", END3, THIRD, ("stop_gained", SHIFT, REGION), False, False),
        ("the intron term kept inside an exon", 998, "ACG", "GTTCAG", END3, THIRD, (SHIFT, REGION, INTRON), False, False),
        ("the terms the other way round", 998, "ACG", "GTTCAG", END3, BACK, THIRD, False, False),
        ("another cDNA position on vep-rs's row", 998, "ACG", "GTTCAG", (*END3, "498-500"), THIRD, BACK, False, False),
    )
    # Pairs whose first column spells the alleles, which every reading of the row has: (the kind, POS in its
    # block, REF, ALT, the cDNA and CDS positions, the reference's terms, vep-rs's terms, set aside)
    SPELT = (
        ("two bases at an exon's end replaced by five, the first column spelling the alleles",
         999, "AC", "GTTCA", END2, THIRD, BACK, True),
        ("two bases of an intron on the reverse strand replaced by eight, the first column spelling the alleles",
         1093, "AC", "GTTCAGGA", NONE, (DONOR, FIFTH, INTRON), (REGION, INTRON), True),
        ("three bases at an exon's start on the reverse strand replaced by five behind a shared first base, the first column spelling the alleles",
         997, "TACG", "TGTTCA", START3, (ACCEPTOR, SHIFT), BACK, True),
        ("an insertion of three between an intron's second and third base, the first column spelling the alleles",
         1002, "A", "AGTT", NONE, (REGION, INTRON), (INTRON,), False),
        ("an insertion of three between two bases of one exon, the first column spelling the alleles",
         990, "A", "AGTT", ("990-991", "890-891"), (DONOR, SHIFT), NOTHING, False),
    )

    @classmethod
    def built(cls) -> list[RecordPair]:
        """Every pair of the two tables, each in its own block, under an identifier or, for a
        pair of `SPELT`, under the first column Ensembl VEP writes for a record without one."""
        pairs = []
        for index, (kind, place, reference, alternate, positions, expected, actual, *verdicts) in enumerate((*cls.PAIRS, *cls.SPELT)):
            position = 100000 + 10000 * index + place
            span, allele = record_row(position, reference, alternate)
            cdna, cds, *other = positions
            name = f"7_{position}_{reference}/{alternate}" if len(verdicts) == 1 else f"made_up_{index}"
            perl, rust = splice_term_rows(",".join(expected), ",".join(actual), allele=allele, location=f"7:{span}", cdna=cdna,
                                          cds=cds, candidate_cdna=other[0] if other else None, name=name)
            pairs.append(RecordPair(kind, position, reference, alternate, f"7:{span}", perl, rust, verdicts[0], verdicts[-1]))
        return pairs

    def test_the_rule_gives_each_kind_its_verdict_with_the_record_and_without_it(self) -> None:
        predicate = row_pair_predicate(self, self.BUCKET)
        verdicts = {}
        for pair in self.built():
            e, a = stop_codon_pair_row(pair.perl), stop_codon_pair_row(pair.rust)
            bases = cvo.replaced_reference(pair.position, pair.reference, [pair.alternate], pair.location, e.allele)
            self.assertIsNotNone(bases, pair.kind)
            if e.reference is None:
                without_record = predicate(e, a)
                e, a = dataclasses.replace(e, reference=bases), dataclasses.replace(a, reference=bases)
                verdicts[pair.kind] = (predicate(e, a), without_record)
            else:
                self.assertEqual((e.reference, a.reference), (bases, bases), pair.kind)
                verdicts[pair.kind] = (predicate(e, a),) * 2
        self.assertEqual(verdicts, {pair.kind: (pair.with_record, pair.without_record) for pair in self.built()})
        self.assertEqual(len(verdicts), len(self.PAIRS) + len(self.SPELT))

    def test_the_comparison_sets_a_kind_aside_by_the_input_records_where_it_is_given_them(self) -> None:
        pairs = self.built()
        # Three positions under an identifier at a place where the input holds a record of another alternate
        # allele: no record writes the row, which is read as a row without one.
        position = 100000 + 10000 * len(pairs) + 998
        location = f"7:{position}-{position + 2}"
        perl, rust = splice_term_rows(",".join(self.THIRD), ",".join(self.BACK), allele="GTTCAG", location=location)
        pairs.append(RecordPair("three positions whose record holds another alternate allele", position, "ACG", "GTTCAGGA",
                                location, perl, rust, True, True))
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, True)
        records = ""
        for pair in pairs:
            name = pair.perl.split("\t", 1)[0]
            identifier = "." if cvo.uploaded_alleles(name) else name
            records += f"7\t{pair.position}\t{identifier}\t{pair.reference}\t{pair.alternate}\t.\t.\t.\n"
        (tmp / "input.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + records, encoding="utf-8")
        perl, rust = "".join(pair.perl for pair in pairs), "".join(pair.rust for pair in pairs)
        for options, verdict in ((("--input-vcf", str(tmp / "input.vcf")), "with_record"), ((), "without_record")):
            adjusted, open_rows = compare_row_outputs(self, perl, rust, *options)
            taken = sum(getattr(pair, verdict) for pair in pairs)
            self.assertEqual(adjusted["excluded_categories"][self.BUCKET], {"count": taken, "excluded": taken}, verdict)
            self.assertEqual(sorted(row.split("\t")[2] for row in open_rows),
                             sorted(2 * [pair.location for pair in pairs if not getattr(pair, verdict)]), verdict)
        # A pair the input alone sets aside has a Location of two positions or a cDNA position that is not two
        # numbers: without the input such a row stays open.
        by_the_input_alone = [stop_codon_pair_row(pair.perl) for pair in pairs if pair.with_record and not pair.without_record]
        self.assertTrue(by_the_input_alone)
        for row in by_the_input_alone:
            numbered = row.fields["cDNA_position"].replace("-", "").isdigit()
            self.assertTrue(cvo.location_span(row) == 2 or not numbered, row.location)


# Two bases across an exon's end replaced by seven, the first column spelling the alleles: the
# span lacks a cDNA position at its end, and the intron's fifth base lies past it.
ACROSS_AN_EXON_END = {"allele": "GTTCAGG", "location": "7:101000-101001", "cdna": "1000-?", "cds": "900-?",
                      "name": "7_101000_AC/GTTCAGG"}


class SpliceTermBesideAStopCodonClassTests(unittest.TestCase):
    """A row can differ by the terms of a stop codon class beside those of
    `splice_term_past_the_record_span`. The class takes such a pair where the stop codon
    class's rule takes the two rows read with vep-rs's splice and intron terms on both, and
    its own reading takes them with the reference's other terms on both. One further
    differing term outside both classes stays a difference, and a pair that differs by the
    stop codon class's terms alone is that class's."""

    BUCKET = "splice_term_past_the_record_span"
    STOP_CODON_BUCKET = "stop_codon_removed_without_stop_codon_test"
    REFERENCE = "splice_donor_variant,splice_donor_5th_base_variant,stop_retained_variant,intron_variant"
    CANDIDATE = "splice_donor_variant,stop_lost,intron_variant"
    FURTHER = "splice_donor_variant,stop_lost,3_prime_UTR_variant,intron_variant"
    STOP_CODON_ALONE = "splice_donor_variant,stop_retained_variant,intron_variant"

    def verdicts(self, reference: str, candidate: str) -> tuple[bool, bool]:
        """Whether this class's rule, and the stop codon class's, takes the pair of rows."""
        e, a = (stop_codon_pair_row(row) for row in splice_term_rows(reference, candidate, **ACROSS_AN_EXON_END))
        self.assertEqual((e.reference, a.reference), ("AC", "AC"))
        return row_pair_predicate(self, self.BUCKET)(e, a), row_pair_predicate(self, self.STOP_CODON_BUCKET)(e, a)

    def counted(self, reference: str, candidate: str) -> tuple[int, int, int]:
        """The pairs the comparison sets aside under this class's bucket and under the stop
        codon class's, and the rows it leaves open."""
        adjusted, open_rows = compare_row_outputs(self, *splice_term_rows(reference, candidate, **ACROSS_AN_EXON_END))
        counts = []
        for bucket in (self.BUCKET, self.STOP_CODON_BUCKET):
            stats = adjusted["excluded_categories"][bucket]
            self.assertEqual(stats["count"], stats["excluded"], bucket)
            counts.append(stats["excluded"])
        return counts[0], counts[1], len(open_rows)

    def test_a_pair_that_differs_by_the_terms_of_both_classes_is_set_aside_under_this_class(self) -> None:
        self.assertEqual(self.verdicts(self.REFERENCE, self.CANDIDATE), (True, False))
        self.assertEqual(self.counted(self.REFERENCE, self.CANDIDATE), (1, 0, 0))

    def test_one_further_term_outside_both_classes_stays_a_difference(self) -> None:
        self.assertEqual(self.verdicts(self.REFERENCE, self.FURTHER), (False, False))
        self.assertEqual(self.counted(self.REFERENCE, self.FURTHER), (0, 0, 2))

    def test_a_pair_that_differs_by_the_stop_codon_class_s_terms_alone_is_that_class_s(self) -> None:
        self.assertEqual(self.verdicts(self.STOP_CODON_ALONE, self.CANDIDATE), (False, True))
        self.assertEqual(self.counted(self.STOP_CODON_ALONE, self.CANDIDATE), (0, 1, 0))


WITH_START_LOST = "frameshift_variant,start_lost,start_retained_variant"
WITHOUT_START_LOST = "frameshift_variant,start_retained_variant"


def fabricated_row(location: str, allele: str, feature: str, consequence: str = "feature_truncation") -> str:
    """One fabricated row of an output: the Location, the allele, the feature and the
    consequence set."""
    return f"v\t{location}\t{allele}\tg1\t{feature}\tTranscript\t{consequence}\n"


def buffer_bound_pairs(digits: int) -> dict[str, tuple[str, str]]:
    """The reference's intergenic row and vep-rs's transcript row of one fabricated record,
    for each coordinate the buffer bound's class reads, a run of ``digits`` nines: the
    record's own position, a span's end and a breakend mate's position."""
    coordinate = "9" * digits
    records = {
        "position": (f"21:{coordinate}", "N[5:30000000["),
        "end": (f"21:5000-{coordinate}", "deletion"),
        "mate": ("21:30425651", f"N[5:{coordinate}["),
    }
    return {
        f"{shape} of {digits} digits": (fabricated_row(location, allele, "-"), fabricated_row(location, allele, "T1"))
        for shape, (location, allele) in records.items()
    }


# A bracket breakend allele whose mate position, the bound 2147478648, is written in the
# decimal digits of another script than `0` to `9`, on a chromosome no cache names.
def other_digits(number: int) -> str:
    """``number`` written in the decimal digits of another script than ``0`` to ``9``."""
    return "".join(chr(0x0660 + int(digit)) for digit in str(number))


OTHER_DIGITS_MATE = f"N[Chr21:{other_digits(2147478648)}["

# Fabricated output pairs, the reference's rows and vep-rs's: nothing left to compare and
# its neighbours, then coordinates around the widths an integer type holds and past the
# length of an integer's text CPython converts by default.
FABRICATED_PAIRS: dict[str, tuple[str, str]] = {
    "no tuple": ("", ""),
    "every tuple shared": (fabricated_row("21:100", "-", "T1", "intron_variant"),) * 2,
    "every tuple set aside as a pair": (
        fabricated_row("21:100", "-", "T1", WITH_START_LOST),
        fabricated_row("21:100", "-", "T1", WITHOUT_START_LOST),
    ),
    "every tuple set aside on one side": ("", fabricated_row("21:1000-2000", "deletion", "T1")),
    "one tuple left": (
        fabricated_row("21:100", "-", "T1", WITH_START_LOST),
        fabricated_row("21:100", "-", "T1", WITHOUT_START_LOST) + fabricated_row("21:300", "-", "T2", "intron_variant"),
    ),
    "one tuple shared, the others set aside": (
        fabricated_row("21:100", "-", "T1", WITH_START_LOST) + fabricated_row("21:500", "-", "T3", "intron_variant"),
        fabricated_row("21:100", "-", "T1", WITHOUT_START_LOST) + fabricated_row("21:500", "-", "T3", "intron_variant"),
    ),
    **{name: pair for digits in (9, 10, 19, 20, 21, 39, 46, 4301) for name, pair in buffer_bound_pairs(digits).items()},
    "a mate position in other digits at the bound": (
        fabricated_row("21:30425651", OTHER_DIGITS_MATE, "-"),
        fabricated_row("21:30425651", OTHER_DIGITS_MATE, "T1"),
    ),
    "a mate position in other digits beside the record's own end": (
        fabricated_row("21:30425651", "N.", "-"),
        fabricated_row("21:30425651", "N.", "-") + fabricated_row("21:30425651", OTHER_DIGITS_MATE, "T1"),
    ),
    "a Location in other digits at the bound": (
        fabricated_row(f"21:{other_digits(2147478648)}", "N[5:30000000[", "-"),
        fabricated_row(f"21:{other_digits(2147478648)}", "N[5:30000000[", "T1"),
    ),
    "a position in other digits beside the same in 0 to 9": (
        fabricated_row(f"21:{other_digits(100)}", "-", "T1", "intron_variant"),
        fabricated_row("21:100", "-", "T1", "intron_variant"),
    ),
    "a superscript digit in a Location": (
        fabricated_row("21:10\u00b2", "-", "T1", "intron_variant"),
        fabricated_row("21:102", "-", "T1", "intron_variant"),
    ),
    "a dot breakend pair at a start of 40 digits": (
        fabricated_row(f"21:{'9' * 40}-{'9' * 46}", ".N", "-"),
        fabricated_row(f"21:{'9' * 40}", ".N", "-"),
    ),
    "a dot breakend pair after a start of 40 digits": (
        fabricated_row(f"21:{'9' * 40}-{'9' * 46}", ".N", "-"),
        fabricated_row(f"21:1{'0' * 40}", ".N", "-"),
    ),
}


def span_rows(location: str) -> tuple[str, str]:
    """One fabricated structural record at ``location``: the reference's rows and vep-rs's,
    which name one transcript alike and one each alone."""

    def rows(*features: str) -> str:
        return "".join(fabricated_row(location, "deletion", feature) for feature in features)

    return rows("T1", "T2"), rows("T1", "T3")


# Fabricated spans past the size the 115.2 scope class reads, the end after a blank: one an
# integer's reading drops (a no-break space) and one it refuses (U+001C).
FABRICATED_SPANS_115_2: dict[str, tuple[str, str]] = {
    "a span past the scope size": span_rows("21:5000-30000000"),
    "a span end after a no-break space": span_rows("21:5000-\u00a030000000"),
    "a span end after a blank an integer's reading refuses": span_rows("21:5000-\x1c30000000"),
}


class NothingLeftToCompareTests(unittest.TestCase):
    """The adjusted figures of an output pair whose every differing tuple a class sets aside:
    precision, recall and F1 read 0.0 where a side holds no tuple, the value `prf` gives a
    ratio over nothing. Every row is fabricated."""

    def adjusted(self, name: str) -> dict:
        return compare_row_outputs(self, *FABRICATED_PAIRS[name])[0]

    def assert_figures(self, adjusted: dict, tuples: int, excluded: tuple[int, int], figure: float) -> None:
        counts = (adjusted["perl_tuple_count"], adjusted["rust_tuple_count"], adjusted["intersection_count"])
        self.assertEqual(counts, (tuples, tuples, tuples))
        self.assertEqual((adjusted["excluded_perl_tuples"], adjusted["excluded_rust_tuples"]), excluded)
        self.assertEqual((adjusted["precision"], adjusted["recall"], adjusted["f1"]), (figure, figure, figure))

    def test_a_pair_whose_every_tuple_is_set_aside_reads_an_adjusted_f1_of_zero(self) -> None:
        self.assert_figures(self.adjusted("every tuple set aside as a pair"), 0, (1, 1), 0.0)
        self.assert_figures(self.adjusted("every tuple set aside on one side"), 0, (0, 1), 0.0)

    def test_a_pair_without_a_tuple_reads_an_adjusted_f1_of_zero(self) -> None:
        self.assert_figures(self.adjusted("no tuple"), 0, (0, 0), 0.0)

    def test_one_shared_tuple_beside_tuples_all_set_aside_reads_an_adjusted_f1_of_one(self) -> None:
        self.assert_figures(self.adjusted("one tuple shared, the others set aside"), 1, (1, 1), 1.0)

    def test_a_ratio_over_nothing_is_zero(self) -> None:
        self.assertEqual(cvo.prf(0, 0, 0), (0.0, 0.0, 0.0))
        self.assertEqual(cvo.prf(0, 0, 3), (0.0, 0.0, 0.0))
        self.assertEqual(cvo.prf(0, 3, 0), (0.0, 0.0, 0.0))


class BufferBoundCoordinateLengthTests(unittest.TestCase):
    """The buffer bound's class reads a coordinate by its value, whatever its length. Every
    row is fabricated."""

    def counted(self, digits: int) -> dict[str, tuple[int, int]]:
        counts = {}
        for name, pair in buffer_bound_pairs(digits).items():
            stats = compare_row_outputs(self, *pair)[0]["excluded_categories"][cvo.BND_MATE_32_BIT_BUFFER_BUCKET]
            counts[name.split()[0]] = (stats["count"], stats["excluded"])
        return counts

    def test_a_coordinate_at_or_past_the_bound_is_counted_whatever_its_length(self) -> None:
        for digits in (10, 19, 20, 21, 39, 46, 4301):
            with self.subTest(digits=digits):
                self.assertEqual(self.counted(digits), {"position": (1, 0), "end": (1, 0), "mate": (1, 0)})

    def test_a_coordinate_below_the_bound_is_not_counted(self) -> None:
        self.assertEqual(self.counted(9), {"position": (0, 0), "end": (0, 0), "mate": (0, 0)})

    def test_a_mate_position_in_other_digits_than_0_to_9_names_no_breakend_mate(self) -> None:
        self.assertIsNone(cvo.bracket_mate(OTHER_DIGITS_MATE))
        self.assertIsNone(cvo.reference_mate_chromosome(OTHER_DIGITS_MATE))
        for name, bucket in (
            ("a mate position in other digits at the bound", cvo.BND_MATE_32_BIT_BUFFER_BUCKET),
            ("a mate position in other digits beside the record's own end", cvo.BND_SYNONYM_MATE_ROWS_LOST_BUCKET),
        ):
            with self.subTest(pair=name):
                stats = compare_row_outputs(self, *FABRICATED_PAIRS[name])[0]["excluded_categories"][bucket]
                self.assertEqual((stats["count"], stats["excluded"]), (0, 0))


class TenDigitsAloneTests(unittest.TestCase):
    """Every reader of a coordinate takes the digits ``0`` to ``9`` alone: a position, a span's
    end, a mate's position or a Location written in the decimal digits of another script is
    no number to it. One test a reader; every row is fabricated."""

    @staticmethod
    def row(location: str, allele: str, cdna: str = "-") -> cvo.PairRow:
        return cvo.pair_row(["rs1", location, allele, "g1", "T1", "Transcript", "intron_variant", cdna, "-"])

    def test_the_location_span_reader(self) -> None:
        self.assertEqual(cvo._location_span("21:5-9"), ("21", 5, 9))
        self.assertIsNone(cvo._location_span(f"21:{other_digits(5)}-{other_digits(9)}"))

    def test_the_reader_of_a_coordinate_that_is_no_integer(self) -> None:
        self.assertFalse(cvo._non_integer_end("21:5-9"))
        self.assertTrue(cvo._non_integer_end(f"21:{other_digits(5)}-{other_digits(9)}"))

    def test_the_rewriting_of_a_location(self) -> None:
        self.assertEqual(cvo.canonical_location("21:005-+9"), "21:5-9")
        written = f"21:{other_digits(5)}-{other_digits(9)}"
        self.assertEqual(cvo.canonical_location(written), written)
        self.assertEqual(cvo.canonical_location("21:10\u00b2"), "21:10\u00b2")

    def test_the_span_size_reader(self) -> None:
        self.assertTrue(cvo._span_exceeds_max_sv_size("21:5-2147478648"))
        self.assertFalse(cvo._span_exceeds_max_sv_size(f"21:5-{other_digits(2147478648)}"))
        self.assertTrue(cvo._span_exceeds_max_sv_size("21:5-\u00a02147478648\u2003"))
        self.assertFalse(cvo._span_exceeds_max_sv_size("21:5-\x1c2147478648"))

    def test_the_mate_position_reader(self) -> None:
        self.assertEqual(cvo.bracket_mate("N[5:30000000["), ("5", 30000000))
        self.assertIsNone(cvo.bracket_mate(f"N[5:{other_digits(30000000)}["))

    def test_the_reader_of_vep_rs_s_position(self) -> None:
        self.assertEqual(cvo.vep_rs_pos("30425650"), 30425650)
        self.assertIsNone(cvo.vep_rs_pos(other_digits(30425650)))

    def test_the_readers_of_a_generated_record_name(self) -> None:
        self.assertEqual(cvo.uploaded_alleles("21_100_A/G"), ("A", ["G"]))
        self.assertEqual(cvo.uploaded_position("21_100_A/G"), 100)
        name = f"21_{other_digits(100)}_A/G"
        self.assertIsNone(cvo.uploaded_alleles(name))
        self.assertIsNone(cvo.uploaded_position(name))
        self.assertIsNotNone(cvo._GENERATED_NAME_RE.match("21_100_A/G,<NON_REF>"))
        self.assertIsNone(cvo._GENERATED_NAME_RE.match(f"21_{other_digits(100)}_A/G,<NON_REF>"))

    def test_the_row_readers_of_a_location_s_ends(self) -> None:
        spans_five = f"21:{other_digits(1)}-{other_digits(5)}"
        self.assertTrue(cvo.indel_row(self.row("21:1-5", "AC")))
        self.assertFalse(cvo.indel_row(self.row(spans_five, "AC")))
        self.assertTrue(cvo.length_changing_row(self.row("21:1-5", "RY")))
        self.assertFalse(cvo.length_changing_row(self.row(spans_five, "RY")))
        self.assertTrue(cvo.insertion_row(self.row("21:1-2", "AC")))
        self.assertFalse(cvo.insertion_row(self.row(f"21:{other_digits(1)}-{other_digits(2)}", "AC")))

    def test_the_readers_of_a_row_s_transcript_positions(self) -> None:
        numbered = list(cvo.intron_placements(3, 20, ("10", "12", "4", "6")))
        unnumbered = list(cvo.intron_placements(3, 20, ("?", "?", "4", "6")))
        self.assertNotEqual(numbered, unnumbered)
        self.assertEqual(list(cvo.intron_placements(3, 20, (other_digits(10), other_digits(12), "4", "6"))), unnumbered)
        self.assertTrue(cvo.longer_allele_inside_one_exon(self.row("21:100-102", "ACGTA", "10-12")))
        other = f"{other_digits(10)}-{other_digits(12)}"
        self.assertFalse(cvo.longer_allele_inside_one_exon(self.row("21:100-102", "ACGTA", other)))

    def test_the_reader_of_the_input_records(self) -> None:
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, True)
        records = tmp / "input.vcf"
        records.write_text("#CHROM\tPOS\tID\tREF\tALT\n21\t100\trs1\tACG\tA\n", encoding="utf-8")
        asked = self.row("21:101-102", "-")
        cvo.read_record_references(records, [asked])
        self.assertEqual(asked.reference, "CG")
        other = self.row(f"21:{other_digits(101)}-{other_digits(102)}", "-")
        cvo.read_record_references(records, [other])
        self.assertIsNone(other.reference)


OUTPUTS_TOOL_NAME = "vep-outputs-compare"


def built_outputs_tool() -> Path | None:
    """The built Rust comparator: the file ``VEP_OUTPUTS_COMPARE`` names, else the binary under
    the cargo target directory or under the tool's own."""
    named = os.environ.get("VEP_OUTPUTS_COMPARE")
    if named:
        return Path(named)
    targets = [Path(__file__).resolve().parents[2] / "tools" / OUTPUTS_TOOL_NAME / "target"]
    if os.environ.get("CARGO_TARGET_DIR"):
        targets.insert(0, Path(os.environ["CARGO_TARGET_DIR"]))
    for target in targets:
        for profile in ("release", "debug"):
            if (target / profile / OUTPUTS_TOOL_NAME).is_file():
                return target / profile / OUTPUTS_TOOL_NAME
    return None


class OutputsToolWritesTheSameReportsTests(unittest.TestCase):
    """The Rust comparator beside this program over the fabricated pairs: the same exit
    status, standard output, standard error and report files, byte for byte. The test runs
    where the tool is built and skips elsewhere."""

    def reports(self, tmp: Path, program: list[str], *options: str) -> dict[str, object]:
        """What one program writes over the pair under ``tmp``, the run's time masked."""
        report = tmp / "report"
        run = subprocess.run(
            [*program, "--perl-dir", str(tmp / "perl"), "--rust-dir", str(tmp / "rust"),
             "--report-dir", str(report), "--assembly", "GRCh38", *options],
            capture_output=True, check=False,
        )
        written: dict[str, object] = {"exit status": run.returncode, "stdout": run.stdout, "stderr": run.stderr}
        for file in sorted(report.iterdir()):
            content = file.read_bytes()
            if file.name == "summary.json":
                stamp = json.loads(content)["generated_at_utc"].encode()
                self.assertEqual(content.count(stamp), 1)
                content = content.replace(stamp, b"")
            written[file.name] = content
        shutil.rmtree(report)
        return written

    def test_the_tool_writes_this_program_s_reports_over_the_fabricated_pairs(self) -> None:
        tool = built_outputs_tool()
        if tool is None:
            self.skipTest(f"no {OUTPUTS_TOOL_NAME} binary: set VEP_OUTPUTS_COMPARE or build tools/{OUTPUTS_TOOL_NAME}")
        cases = [(name, pair, ()) for name, pair in FABRICATED_PAIRS.items()]
        cases += [(name, pair, ("--reference-release", "115.2")) for name, pair in FABRICATED_SPANS_115_2.items()]
        for name, (perl, rust), options in cases:
            with self.subTest(pair=name):
                tmp = Path(tempfile.mkdtemp())
                self.addCleanup(shutil.rmtree, tmp, True)
                for side, rows in (("perl", perl), ("rust", rust)):
                    (tmp / side).mkdir()
                    (tmp / side / "output.txt").write_text(ROW_COLS + rows, encoding="utf-8")
                own = self.reports(tmp, [sys.executable, cvo.__file__], *options)
                self.assertEqual(own["exit status"], 0, own["stderr"])
                self.assertEqual(sorted(own)[:4], ["discordant.tsv", "discordant_open.tsv", "exit status", "stderr"])
                self.assertEqual(self.reports(tmp, [str(tool), *options]), own)


if __name__ == "__main__":
    unittest.main(verbosity=2)
