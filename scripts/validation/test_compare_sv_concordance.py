#!/usr/bin/env python3
"""Tests for the SV concordance comparator.

`compare_sv_concordance.py` produces every reported SV concordance figure and the
per-class exclusion inventory behind them. Every test here pins a measurable
number rather than an abstract style property.

Run:
    python3 -m pytest scripts/validation/test_compare_sv_concordance.py -q
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import compare_sv_concordance as C


def T(loc: str, allele: str, feature: str, csq: str, ftype: str = "Transcript"):
    """Build one annotation tuple with the comparator's own normalization."""
    return (
        C.normalize_location(loc),
        C.normalize_allele(allele),
        feature,
        ftype,
        C.normalize_consequence_set(csq),
    )


def write_vep(path: Path, rows: list[tuple[str, ...]], header: bool = True) -> None:
    """Write a minimal 17-column VEP tab output."""
    with open(path, "w") as fh:
        if header:
            fh.write("## fastVEP output\n")
            fh.write(
                "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t"
                "Consequence\tcDNA_position\tCDS_position\tProtein_position\t"
                "Amino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\n"
            )
        fh.writelines(
            "\t".join([uploaded, loc, allele, gene, feature, ftype, csq] + ["-"] * 10)
            + "\n"
            for uploaded, loc, allele, gene, feature, ftype, csq in rows
        )


# F1 arithmetic, hand-computed


class MetricsArithmeticTests(unittest.TestCase):
    def test_perfect_agreement(self):
        t = {T("21:100", "<DEL>", "ENST1", "feature_truncation")}
        m = C.compute_metrics([], t, set(), set(t), set())
        self.assertEqual((m.precision, m.recall, m.f1), (1.0, 1.0, 1.0))

    def test_hand_computed_f1(self):
        # Perl 4 tuples, Rust 3, intersection 2.
        # P = 2/3, R = 2/4 = 0.5, F1 = 2*(2/3)*(0.5)/((2/3)+0.5) = 0.571428...
        perl = {T("21:1", "A", f"E{i}", "intron_variant") for i in range(4)}
        rust = {T("21:1", "A", f"E{i}", "intron_variant") for i in range(2)} | {
            T("21:1", "A", "E9", "intron_variant")
        }
        m = C.compute_metrics([], perl, set(), rust, set())
        self.assertEqual(m.perl_tuples, 4)
        self.assertEqual(m.rust_tuples, 3)
        self.assertEqual(m.intersection, 2)
        self.assertAlmostEqual(m.f1, 0.5714285714285714, places=12)

    def test_both_empty_scores_zero_not_one(self):
        """A comparison with no tuples on either side must NOT read as perfect.

        Two missing engine outputs otherwise produce F1 = 1.0 on zero tuples,
        which is indistinguishable in a report from genuine agreement.
        """
        m = C.compute_metrics([], set(), set(), set(), set())
        self.assertEqual(m.f1, 0.0)

    def test_one_sided_output_scores_zero(self):
        perl = {T("21:1", "A", "E1", "intron_variant")}
        m = C.compute_metrics([], perl, set(), set(), set())
        self.assertEqual(m.f1, 0.0)

    def test_recompute_prf_matches_compute_metrics(self):
        perl = {T("21:1", "A", f"E{i}", "intron_variant") for i in range(7)}
        rust = {T("21:1", "A", f"E{i}", "intron_variant") for i in range(5)}
        m = C.compute_metrics([], perl, set(), rust, set())
        expected = m.f1
        m.recompute_prf()
        self.assertAlmostEqual(m.f1, expected, places=15)


# Malformed or truncated output must fail, not silently shrink


class ParseRobustnessTests(unittest.TestCase):
    def setUp(self):
        import tempfile

        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_wellformed_file_parses(self):
        p = self.tmp / "ok.txt"
        write_vep(
            p,
            [
                (
                    "v1",
                    "21:100",
                    "<DEL>",
                    "G1",
                    "ENST1",
                    "Transcript",
                    "feature_truncation",
                ),
                (
                    "v2",
                    "21:200",
                    "<DUP>",
                    "G2",
                    "ENST2",
                    "Transcript",
                    "feature_elongation",
                ),
            ],
        )
        tuples, locs = C.parse_vep_output(p)
        self.assertEqual(len(tuples), 2)
        self.assertEqual(locs, {"21:100", "21:200"})

    def test_truncated_midline_raises(self):
        """A file cut mid-line must raise, not return a smaller tuple set.

        Skipping the short row yields a plausible-but-smaller set, and a smaller
        set with an unchanged intersection reads as BETTER concordance.
        """
        p = self.tmp / "trunc.txt"
        write_vep(
            p,
            [
                (
                    "v1",
                    "21:100",
                    "<DEL>",
                    "G1",
                    "ENST1",
                    "Transcript",
                    "feature_truncation",
                )
            ],
        )
        with open(p, "a") as fh:
            fh.write("v2\t21:200\t<DUP>\tG2\n")  # only 4 columns
        with self.assertRaises(ValueError) as ctx:
            C.parse_vep_output(p)
        self.assertIn("expected at least 7", str(ctx.exception))

    def test_missing_file_returns_empty_without_raising(self):
        tuples, locs = C.parse_vep_output(self.tmp / "does_not_exist.txt")
        self.assertEqual((tuples, locs), (set(), set()))

    def test_header_and_blank_lines_skipped(self):
        p = self.tmp / "hdr.txt"
        write_vep(
            p,
            [
                (
                    "v1",
                    "21:100",
                    "<DEL>",
                    "G1",
                    "ENST1",
                    "Transcript",
                    "feature_truncation",
                )
            ],
        )
        with open(p, "a") as fh:
            fh.write("\n")
        tuples, _ = C.parse_vep_output(p)
        self.assertEqual(len(tuples), 1)


# The transcript-selection mask must be direction-aware


class IntendedDivergenceDirectionTests(unittest.TestCase):
    # A span above Perl's --max_sv_size default (10 Mb), which is the only regime
    # where the masked defect can occur: below it Perl loads every region the
    # variant overlaps, so no batch dependence exists.
    GIANT = "21:1000000-21000000"
    # A span below that threshold. Perl requests every region it needs here.
    SMALL = "21:1000-2000"

    def test_masks_perl_under_annotation_the_documented_defect(self):
        """Perl names FEWER transcripts on a giant SV (skipped region loading).

        This is the shape the filter exists for: vep-rs annotates the complete
        transcript set, Perl sees only what a neighbouring variant loaded.
        """
        csq = "feature_truncation"
        perl = {T(self.GIANT, "<DEL>", "ENST_A", csq)}
        rust = {T(self.GIANT, "<DEL>", f"ENST_{c}", csq) for c in "BCDE"} | perl
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        self.assertEqual(len(ex_rust), 4, "vep-rs extra transcripts should be masked")
        self.assertEqual(len(ex_perl), 0)

    def test_does_not_mask_vep_rs_false_positive_on_small_sv(self):
        """A spurious vep-rs transcript below --max_sv_size must NOT be masked.

        Perl annotated this variant correctly and completely; a transcript vep-rs
        named and Perl did not is a vep-rs false positive, and the batch-dependent
        mechanism cannot reach a variant this small, so without the scope gate the
        false positive would cost vep-rs nothing.
        """
        csq = "feature_truncation"
        perl = {T(self.SMALL, "<DEL>", f"ENST_{c}", csq) for c in "AB"}
        rust = perl | {T(self.SMALL, "<DEL>", "ENST_SPURIOUS", csq)}
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        self.assertEqual(
            len(ex_rust),
            0,
            "a vep-rs-only transcript on a sub-threshold SV is a false positive, "
            "not Perl under-annotation, and must stay in both scales",
        )
        self.assertEqual(len(ex_perl), 0)

    def test_does_not_mask_vep_rs_false_positive_at_point_locus(self):
        """Same, at a point location, where the span is 0 by construction."""
        csq = "feature_truncation"
        perl = {T("21:5000", "A", f"ENST_{c}", csq) for c in "AB"}
        rust = perl | {T("21:5000", "A", "ENST_SPURIOUS", csq)}
        ex_rust, _ = C.filter_intended_divergences(perl, rust)
        self.assertEqual(len(ex_rust), 0)

    # Spans that straddle Perl's --max_sv_size default. The existing GIANT/SMALL
    # cases sit 20 Mb above and 10 Mb below it, so neither exercises the boundary.
    # The end coordinate is start + span: a span of 10,000,000 from start 1000 ends
    # at 10,001,000, NOT 11,000,000.
    AT_THRESHOLD = "21:1000-10001000"  # span exactly 10,000,000
    OVER_THRESHOLD = "21:1000-10001001"  # span 10,000,001

    def test_span_exactly_at_max_sv_size_is_not_over_threshold(self):
        """Perl skips at `$len > $max_sv_size`, so equality is NOT a skip.

        `Parser.pm:493` compares `end - start` strictly greater than the 10,000,000
        default from `Config.pm:310`. A gate using `>=` would mask a class of
        variants Perl annotated completely.
        """
        self.assertFalse(C._span_exceeds_max_sv_size(self.AT_THRESHOLD))

    def test_span_one_over_max_sv_size_is_over_threshold(self):
        """One base past the threshold is where Perl's skip begins."""
        self.assertTrue(C._span_exceeds_max_sv_size(self.OVER_THRESHOLD))

    def test_equal_transcript_counts_are_asymmetric_by_design(self):
        """Equal-sized disjoint sets: the two arms deliberately differ.

        The vep-rs arm requires vep-rs to have named STRICTLY more, so it
        excludes nothing here: with equal counts there is no evidence Perl loaded
        fewer, and a vep-rs-only transcript could equally be a false positive.

        The Perl arm keeps a `>=` contract (it skips only when vep-rs named
        strictly FEWER), so equal counts still mask, treating disjoint
        equal-size sets as a buffer artifact. That asymmetry is intentional and
        conservative in both directions: the tighter test governs the side where a
        wrong answer would flatter vep-rs.
        """
        csq = "feature_truncation"
        perl = {T(self.GIANT, "<DEL>", f"ENST_P{i}", csq) for i in range(3)}
        rust = {T(self.GIANT, "<DEL>", f"ENST_R{i}", csq) for i in range(3)}
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        self.assertEqual(len(ex_rust), 0, "vep-rs arm requires strictly more")
        self.assertEqual(len(ex_perl), 3, "Perl arm masks on equal counts")

    def test_full_location_key_separates_records_sharing_a_start(self):
        """Two records at one start with different ends are distinct variants.

        Keying on the start alone merged their transcript sets, so a transcript
        named on one could suppress a divergence on the other.
        """
        csq = "feature_truncation"
        a = "21:1000000-21000000"  # giant, in scope
        b = "21:1000000-1000500"  # same start, small, out of scope
        perl = {T(a, "<DEL>", "ENST_A", csq), T(b, "<DEL>", "ENST_B", csq)}
        rust = perl | {T(b, "<DEL>", "ENST_SPURIOUS_ON_B", csq)}
        ex_rust, _ = C.filter_intended_divergences(perl, rust)
        self.assertEqual(
            len(ex_rust),
            0,
            "the spurious transcript belongs to the SMALL record and must not "
            "inherit the giant record's scope",
        )

    def test_does_not_mask_vep_rs_under_annotation(self):
        """vep-rs names FEWER transcripts: the OPPOSITE shape, an open vep-rs gap.

        Reproduces the shape of the 25 kb <DEL> at 21:23699593-23724683 (Perl 122
        transcripts vs vep-rs 2). Masking this direction converts a coverage gap
        into an excluded transcript-selection divergence and inflates adjusted F1.
        """
        csq = "feature_truncation,coding_sequence_variant"
        loc = "21:23699593-23724683"
        perl = {T(loc, "<DEL>", f"ENST_P{i}", csq) for i in range(122)}
        rust = {T(loc, "<DEL>", f"ENST_R{i}", csq) for i in range(2)}
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        self.assertEqual(
            len(ex_perl),
            0,
            "Perl-only tuples must NOT be masked when vep-rs named fewer "
            "transcripts; that direction is a vep-rs gap, not a Perl defect",
        )

    def test_masks_perl_only_tuples_when_vep_rs_names_more(self):
        """The Perl-side arm still fires in the defect regime.

        Perl named 2 transcripts on a giant SV, vep-rs named 4 disjoint ones. The
        Perl-only tuples are the ones Perl's incidental cache load produced.
        """
        csq = "feature_truncation"
        perl = {T(self.GIANT, "<DEL>", f"ENST_P{i}", csq) for i in range(2)}
        rust = {T(self.GIANT, "<DEL>", f"ENST_R{i}", csq) for i in range(4)}
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        self.assertEqual(len(ex_perl), 2)
        self.assertEqual(len(ex_rust), 4)

    def test_no_mask_when_other_engine_never_saw_the_variant(self):
        """A variant only one engine annotated at all is a real miss."""
        perl = {T("21:1", "<DEL>", "ENST_A", "feature_truncation")}
        rust = {T("21:999", "<DEL>", "ENST_B", "feature_truncation")}
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        self.assertEqual(len(ex_rust), 0)
        self.assertEqual(len(ex_perl), 0)

    def test_same_feature_different_consequence_is_not_a_selection_divergence(self):
        """Same transcript, different terms: a consequence swap, not selection."""
        loc = "21:1000-2000"
        perl = {T(loc, "<DEL>", "ENST_A", "feature_truncation")}
        rust = {T(loc, "<DEL>", "ENST_A", "transcript_ablation")}
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        self.assertEqual(len(ex_rust), 0)
        self.assertEqual(len(ex_perl), 0)


# Consequence-normalization mask


class CrossChromosomeMaskTests(unittest.TestCase):
    """Perl naming a transcript that is not on the variant's chromosome.

    Given a chr21 structural variant whose `INFO/CHR2` reads `chr21` while `CHROM` reads
    `21`, Ensembl VEP r115.2 annotates it against chromosome 22 transcripts at numerically
    similar coordinates. On `gnomAD-SV_v3_DEL_chr21_55f8f183` it names 122 transcripts of
    which 120 resolve to chr22. Unlike every other mask here the criterion is objectively
    checkable rather than mechanistic, and the authority is the vep-rs JSON cache both
    engines were run against.
    """

    TX_BY_CHR = {"21": {"ENST_CHR21_A", "ENST_CHR21_B"}, "22": {"ENST_CHR22_X"}}

    def test_perl_naming_an_off_chromosome_transcript_is_excluded(self) -> None:
        loc = "21:23699593-23724683"
        perl = {
            T(loc, "deletion", "ENST_CHR21_A", "transcript_ablation"),
            T(loc, "deletion", "ENST_CHR22_X", "transcript_ablation"),
        }
        rust = {T(loc, "deletion", "ENST_CHR21_A", "transcript_ablation")}
        ex_rust, ex_perl, _ = C.filter_cross_chromosome_divergences(perl, rust, self.TX_BY_CHR)
        self.assertEqual(ex_rust, set(), "the mask is Perl-side only")
        self.assertEqual(
            ex_perl, {T(loc, "deletion", "ENST_CHR22_X", "transcript_ablation")}
        )

    def test_a_transcript_on_the_variants_chromosome_is_never_excluded(self) -> None:
        """The mask must not reach a genuine divergence. `ENST_CHR21_B` is a real chr21
        transcript vep-rs missed, which is a vep-rs gap and belongs in its denominator."""
        loc = "21:1000-2000"
        perl = {T(loc, "deletion", "ENST_CHR21_B", "feature_truncation")}
        _, ex_perl, _ = C.filter_cross_chromosome_divergences(perl, set(), self.TX_BY_CHR)
        self.assertEqual(ex_perl, set())

    def test_bracket_notation_breakend_alleles_are_skipped(self) -> None:
        """A mate coordinate legitimately names another chromosome, so a chr22 transcript
        on a `N[chr22:...[` allele is correct output, not the defect."""
        loc = "21:31916312"
        perl = {T(loc, "N[chr22:44739371[", "ENST_CHR22_X", "intergenic_variant")}
        _, ex_perl, _ = C.filter_cross_chromosome_divergences(perl, set(), self.TX_BY_CHR)
        self.assertEqual(ex_perl, set())

    def test_a_chr_prefixed_location_resolves_to_the_same_chromosome(self) -> None:
        """The two corpora disagree on the prefix, and the mask must not depend on it."""
        for loc in ("21:5000-6000", "chr21:5000-6000"):
            perl = {T(loc, "deletion", "ENST_CHR22_X", "transcript_ablation")}
            _, ex_perl, _ = C.filter_cross_chromosome_divergences(perl, set(), self.TX_BY_CHR)
            self.assertEqual(len(ex_perl), 1, msg=loc)

    def test_an_unknown_chromosome_is_not_adjudicated(self) -> None:
        """No index for the chromosome means the cache cannot judge, so nothing is
        excluded. Excluding here would mask on absence of evidence."""
        perl = {T("7:5000-6000", "deletion", "ENST_CHR22_X", "transcript_ablation")}
        _, ex_perl, _ = C.filter_cross_chromosome_divergences(perl, set(), self.TX_BY_CHR)
        self.assertEqual(ex_perl, set())

    def test_an_empty_authority_disables_the_mask(self) -> None:
        """With no cache the mask excludes nothing rather than guessing. The report
        records that it did not run, so a zero cannot read as a clean result."""
        perl = {T("21:5000-6000", "deletion", "ENST_CHR22_X", "transcript_ablation")}
        _, ex_perl, _ = C.filter_cross_chromosome_divergences(perl, set(), {})
        self.assertEqual(ex_perl, set())

    def test_vep_rs_producing_the_shape_raises_instead_of_being_masked(self) -> None:
        """vep-rs queries a per-chromosome index, so it cannot name an off-chromosome
        transcript. If it ever does that is a vep-rs defect, and masking it as a VEP
        defect would hide it: the whole point of this mask is that VEP is the one at
        fault.
        """
        loc = "21:5000-6000"
        rust = {T(loc, "deletion", "ENST_CHR22_X", "feature_truncation")}
        perl = {T(loc, "deletion", "ENST_CHR21_A", "feature_truncation")}
        with self.assertRaises(AssertionError) as ctx:
            C.filter_cross_chromosome_divergences(perl, rust, self.TX_BY_CHR)
        self.assertIn("vep-rs bug", str(ctx.exception))

    def test_a_non_vep_rs_scored_engine_is_counted_not_raised(self) -> None:
        """A non-vep-rs engine builds its own transcript set from its own annotation
        release, so it legitimately names transcripts this cache does not hold for the
        chromosome. Raising there scored NOTHING at all, which is what kept this filter
        off the fastVEP arm and left that arm's adjusted F1 on a different exclusion set.
        The tuple is counted and stays in raw AND adjusted; only the Perl side is masked.
        """
        loc = "21:5000-6000"
        rust = {T(loc, "deletion", "ENST_NEWER_RELEASE", "feature_truncation")}
        perl = {T(loc, "deletion", "ENST_CHR22_X", "feature_truncation")}
        ex_rust, ex_perl, absent = C.filter_cross_chromosome_divergences(
            perl, rust, self.TX_BY_CHR, "fastvep"
        )
        self.assertEqual(ex_rust, set(), "no scored-engine tuple is ever excluded here")
        self.assertEqual(
            ex_perl,
            {T(loc, "deletion", "ENST_CHR22_X", "feature_truncation")},
            "the Perl-side mask must still fire for a non-vep-rs scored engine",
        )
        self.assertEqual(absent, 1)

    def test_the_scored_engine_default_keeps_the_vep_rs_assertion(self) -> None:
        """Omitting the argument must not silently relax the assertion: the default is the
        strict vep-rs behaviour, so a caller that passes no engine keeps the detector."""
        loc = "21:5000-6000"
        rust = {T(loc, "deletion", "ENST_CHR22_X", "feature_truncation")}
        with self.assertRaises(AssertionError):
            C.filter_cross_chromosome_divergences(set(), rust, self.TX_BY_CHR)

    def test_non_transcript_features_are_ignored(self) -> None:
        """Regulatory features and motif ids are not `ENST` and carry no chromosome
        claim this mask can check."""
        perl = {T("21:5000-6000", "deletion", "ENSR00000123456", "regulatory_region_variant")}
        _, ex_perl, _ = C.filter_cross_chromosome_divergences(perl, set(), self.TX_BY_CHR)
        self.assertEqual(ex_perl, set())

    def test_cache_loader_reads_both_shard_layouts(self) -> None:
        """Storable-derived caches nest transcripts under a chromosome key; native ones
        are a bare list. The loader must read both, since the SV suites run on the
        Storable-derived caches."""
        import json
        import tempfile

        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "transcripts"
            (root / "21").mkdir(parents=True)
            (root / "22").mkdir(parents=True)
            (root / "21" / "0-1000.json").write_text(
                json.dumps([{"stable_id": "ENST_A"}, {"stable_id": "ENST_B"}])
            )
            (root / "22" / "0-1000.json").write_text(
                json.dumps({"22": [{"stable_id": "ENST_C"}, None]})
            )
            got = C.load_transcripts_by_chromosome(d)
        self.assertEqual(got, {"21": {"ENST_A", "ENST_B"}, "22": {"ENST_C"}})

    def test_cache_loader_returns_empty_for_a_missing_directory(self) -> None:
        self.assertEqual(C.load_transcripts_by_chromosome("/nonexistent/cache"), {})


class SvFilterSetTests(unittest.TestCase):
    """The comparator's masks are exactly the `filter_*_divergences` definitions.

    Every mask is a function of that name shape and every function of that name
    shape is a mask, so the set of names is the set of masks: two full masks
    (giant-SV transcript selection, cross-chromosome annotation) and one partial
    (<CNV:TR>). A mask over pairs on which the two engines AGREE is unreachable by
    construction and would hide a divergence rather than mask one, so no
    consequence-normalisation mask exists.
    """

    def test_the_sv_filter_set_is_exactly_the_three_expected(self) -> None:
        """Pinned as a set rather than a count, so adding or removing a filter fails
        here with its name, beside the taxonomy row it would need."""
        self.assertEqual(
            {n for n in dir(C) if n.startswith("filter_") and n.endswith("_divergences")},
            {
                "filter_intended_divergences",
                "filter_cnv_tr_expansion_divergences",
                "filter_cross_chromosome_divergences",
            },
        )


# Cross-file tuple collisions


class CrossFileAggregationTests(unittest.TestCase):
    def test_identical_tuple_in_two_files_collapses_under_union(self):
        """Demonstrates the aggregation hazard the file-summed counts address.

        Two different input VCFs of one assembly can emit a byte-identical tuple.
        A set union keeps one, while `total_input` sums both, so an aggregate row
        built from both mixes two aggregations, and on the real corpus the
        collapse moves raw F1 at the third decimal on GRCh38.
        """
        shared = T("21:33036920", "A.", "ENST00000609934", "intron_variant")
        file_a = {shared, T("21:1", "A", "E1", "intron_variant")}
        file_b = {shared, T("21:2", "A", "E2", "intron_variant")}
        union = file_a | file_b
        filesum = len(file_a) + len(file_b)
        self.assertEqual(filesum, 4)
        self.assertEqual(len(union), 3)
        self.assertEqual(filesum - len(union), 1)

    def test_f1_differs_between_the_two_aggregations(self):
        shared = T("21:100", "<DEL>", "ENST_S", "feature_truncation")
        perl_a = {shared}
        perl_b = {shared}
        rust_a = {shared}
        rust_b = {shared, T("21:200", "<DEL>", "ENST_X", "feature_truncation")}
        union = C.compute_metrics([], perl_a | perl_b, set(), rust_a | rust_b, set())
        # File-summed: Perl 2, Rust 3, intersection 2.
        p, r = 2 / 3, 2 / 2
        filesum_f1 = 2 * p * r / (p + r)
        self.assertNotAlmostEqual(union.f1, filesum_f1, places=6)


# End-to-end: BOTH aggregate rows must use the file-summed basis.
#
# The two tests above assert the arithmetic PREMISE of the file-summed counts by
# hand and never call main(), so they stay green if main() drops its file-sum
# override for either row. These drive the real main() over a two-file corpus and
# assert the aggregation of the emitted report, so removing either override fails
# here.


class AggregateBasisEndToEndTests(unittest.TestCase):
    """main() must report BOTH overall and adjusted on the file-summed basis."""

    def setUp(self):
        import tempfile

        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        for sub in ("in", "perl", "rust", "out"):
            (self.tmp / sub).mkdir()

        # One tuple deliberately emitted by BOTH files: the cross-file duplicate
        # whose double-counting is the whole point of the two aggregations.
        dup = (
            "dupvar",
            "21:500",
            "<DEL>",
            "G0",
            "ENST_DUP",
            "Transcript",
            "feature_truncation",
        )
        # A tuple only Perl emits, and one only vep-rs emits, so the intersection
        # is a strict subset of both sides and F1 is neither 0 nor 1.
        perl_only = (
            "pv",
            "21:600",
            "<DEL>",
            "G1",
            "ENST_P",
            "Transcript",
            "feature_truncation",
        )
        rust_only = (
            "rv",
            "21:700",
            "<DUP>",
            "G2",
            "ENST_R",
            "Transcript",
            "feature_elongation",
        )

        for name, extra_perl, extra_rust in (
            ("fileA", [perl_only], [rust_only]),
            ("fileB", [], []),
        ):
            # A minimal VCF so the file is discovered and total_input is non-zero.
            with open(self.tmp / "in" / f"{name}.vcf", "w") as fh:
                fh.write("##fileformat=VCFv4.2\n")
                fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                fh.write("21\t500\t.\tG\t<DEL>\t.\tPASS\tEND=600\n")
            write_vep(self.tmp / "perl" / f"{name}.txt", [dup] + extra_perl)
            write_vep(self.tmp / "rust" / f"{name}.txt", [dup] + extra_rust)

        self.report = self._run_main()

    def tearDown(self):
        self._tmp.cleanup()

    def _run_main(self) -> dict:
        import json

        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py",
            "--assembly",
            "grch37",
            "--input-dir",
            str(self.tmp / "in"),
            "--perl-dir",
            str(self.tmp / "perl"),
            "--rust-dir",
            str(self.tmp / "rust"),
            "--output-dir",
            str(self.tmp / "out"),
        ]
        try:
            C.main()
        except SystemExit as exc:  # main() exits 0 on success
            if exc.code not in (0, None):
                raise
        finally:
            sys.argv = argv
        with open(self.tmp / "out" / "concordance_report.json") as fh:
            return json.load(fh)

    # The corpus: Perl emits dup twice + perl_only once = 3 file-summed, 2 unique.
    # vep-rs emits dup twice + rust_only once = 3 file-summed, 2 unique.
    # The intersection is dup, appearing in both files = 2 file-summed, 1 unique.

    def test_overall_uses_file_summed_counts(self):
        ov = self.report["overall"]
        self.assertEqual(ov["perl_tuples"], 3, "raw Perl count is not the file-sum")
        self.assertEqual(ov["rust_tuples"], 3, "raw vep-rs count is not the file-sum")
        self.assertEqual(ov["intersection"], 2, "raw intersection is not the file-sum")

    def test_overall_reports_the_union_alongside(self):
        ov = self.report["overall"]
        self.assertEqual(ov["perl_tuples_union"], 2)
        self.assertEqual(ov["rust_tuples_union"], 2)
        self.assertEqual(ov["intersection_union"], 1)
        self.assertEqual(ov["cross_file_duplicate_perl"], 1)
        self.assertEqual(ov["cross_file_duplicate_rust"], 1)

    def test_adjusted_uses_the_same_basis_as_raw(self):
        """The adjusted row must sit on the same basis as the raw one.

        No exclusion predicate fires on this fixture, so with both rows on one
        basis the adjusted counts equal the raw counts exactly. Under the union
        basis they would instead be the smaller deduplicated figures.
        """
        ov, adj = self.report["overall"], self.report["adjusted"]
        self.assertEqual(
            (adj["perl_tuples"], adj["rust_tuples"], adj["intersection"]),
            (ov["perl_tuples"], ov["rust_tuples"], ov["intersection"]),
            "adjusted aggregate is on a different basis from raw",
        )

    def test_adjusted_f1_equals_raw_f1_when_nothing_is_excluded(self):
        """Two bases in one row let adjusted F1 beat raw F1 with an empty mask.

        That is the reader-visible symptom: the Adj column looks better than Raw
        purely from the aggregation mismatch, not from any excluded divergence.
        """
        adj = self.report["adjusted"]
        self.assertEqual(adj["excluded_divergences"], 0)
        self.assertEqual(adj["excluded_cnvtr_rust"], 0)
        # The report carries no `excluded_csq_rust` or Perl twin, because there is no
        # BND non-coding normalisation layer. Their absence is asserted so an added
        # field cannot quietly start reporting a mask the taxonomy lacks.
        self.assertNotIn("excluded_csq_rust", adj)
        self.assertNotIn("excluded_csq_perl", adj)
        self.assertAlmostEqual(adj["f1"], self.report["overall"]["f1"], places=6)

    def test_adjusted_also_reports_its_union_counts(self):
        adj = self.report["adjusted"]
        for key in (
            "perl_tuples_union",
            "rust_tuples_union",
            "intersection_union",
            "cross_file_duplicate_perl",
            "cross_file_duplicate_rust",
        ):
            self.assertIn(key, adj, f"adjusted row is missing {key}")
        self.assertEqual(adj["cross_file_duplicate_perl"], 1)

    def test_per_file_rows_are_unchanged_by_the_aggregate_override(self):
        """Per-file rows must stay per-file: no union fields, no summed totals."""
        for section in ("per_file", "per_file_adjusted"):
            for name, m in self.report[section].items():
                self.assertEqual(
                    m["perl_tuples"],
                    2 if name == "fileA" else 1,
                    f"{section}/{name} per-file Perl count changed",
                )
                self.assertNotIn(
                    "perl_tuples_union",
                    m,
                    f"{section}/{name} leaked an aggregate-only union field",
                )


# Normalization helpers


class NormalizationTests(unittest.TestCase):
    def test_chr_prefix_stripped(self):
        self.assertEqual(C.normalize_location("chr21:100"), "21:100")

    def test_chrM_becomes_MT(self):
        self.assertEqual(C.normalize_location("chrM:100"), "MT:100")

    def test_bnd_bracket_chr_stripped(self):
        self.assertEqual(C.normalize_allele("N[chr21:1234["), "N[21:1234[")

    def test_consequence_set_sorted_and_deduped(self):
        self.assertEqual(
            C.normalize_consequence_set("intron_variant,stop_gained,intron_variant"),
            "intron_variant,stop_gained",
        )

    def test_consequence_order_does_not_affect_identity(self):
        a = T("21:1", "A", "E1", "b_term,a_term")
        b = T("21:1", "A", "E1", "a_term,b_term")
        self.assertEqual(a, b)

    def test_variant_classification(self):
        self.assertEqual(C.classify_variant("A", "T"), "SNV")
        self.assertEqual(C.classify_variant("N", "<DEL>"), "Symbolic_DEL")
        self.assertEqual(C.classify_variant("N", "<DUP:TANDEM>"), "Symbolic_DUP")
        self.assertEqual(C.classify_variant("N", "<CN=3>"), "Symbolic_CNV")
        self.assertEqual(C.classify_variant("N", "<INS:ME:ALU>"), "Mobile_Element")
        self.assertEqual(C.classify_variant("N", "N[21:100["), "BND")
        self.assertEqual(C.classify_variant("N", "<NON_REF>"), "NON_REF")
        self.assertEqual(C.classify_variant("A", "*"), "Spanning")


# Symmetry / invariance properties


class PropertyTests(unittest.TestCase):
    def test_precision_recall_swap_under_engine_swap(self):
        perl = {T("21:1", "A", f"E{i}", "intron_variant") for i in range(5)}
        rust = {T("21:1", "A", f"E{i}", "intron_variant") for i in range(3)}
        a = C.compute_metrics([], perl, set(), rust, set())
        b = C.compute_metrics([], rust, set(), perl, set())
        self.assertAlmostEqual(a.precision, b.recall, places=15)
        self.assertAlmostEqual(a.recall, b.precision, places=15)
        self.assertAlmostEqual(a.f1, b.f1, places=15)

    def test_mask_never_changes_the_intersection(self):
        """A swap-mode mask removes only discordants, so agreement is fixed."""
        csq = "feature_truncation"
        loc = "21:1000-2000"
        shared = T(loc, "<DEL>", "ENST_SHARED", csq)
        perl = {shared, T(loc, "<DEL>", "ENST_P", csq)}
        rust = {shared, T(loc, "<DEL>", "ENST_R", csq)}
        base = C.compute_metrics([], perl, set(), rust, set())
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        adj = C.compute_metrics([], perl - ex_perl, set(), rust - ex_rust, set())
        self.assertEqual(base.intersection, adj.intersection)
        self.assertIn(shared, perl - ex_perl)
        self.assertIn(shared, rust - ex_rust)

    def test_mask_cannot_lower_f1(self):
        csq = "feature_truncation"
        loc = "21:1000-2000"
        shared = T(loc, "<DEL>", "ENST_SHARED", csq)
        perl = {shared} | {T(loc, "<DEL>", f"ENST_P{i}", csq) for i in range(3)}
        rust = {shared} | {T(loc, "<DEL>", f"ENST_R{i}", csq) for i in range(3)}
        base = C.compute_metrics([], perl, set(), rust, set())
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        adj = C.compute_metrics([], perl - ex_perl, set(), rust - ex_rust, set())
        self.assertGreaterEqual(adj.f1, base.f1)


# fastVEP symbolic-allele normalisation (--scored-engine fastvep)
#
# Perl VEP and vep-rs write an SV's Allele column as its SO class string
# (`deletion`, `tandem_duplication`, `Alu_insertion`, ...); fastVEP writes the
# raw VCF token (`<DEL>`, `<DUP:TANDEM>`, `<INS:ME:ALU>`, ...). Without the mapping,
# every symbolic-allele tuple on the fastVEP arm is a guaranteed miss whatever its
# consequence set, so that arm's SV F1 measures a representation difference rather
# than annotation agreement. The mapping is derived from vep-rs's own
# `classify_symbolic_alt` + `VariantClass::so_term` + `display_allele`, and is gated
# on the scored engine so it can never touch the Perl reference or the vep-rs arm.


class FastvepSymbolicAlleleTests(unittest.TestCase):
    # The exact table the comparator implements, token -> Perl/vep-rs Allele string.
    # Pinned as a whole so a change to any row fails here by name. Keys are the
    # canonical upper-case spellings; case-insensitivity is tested separately.
    MAPPING = {
        "<DEL>": "deletion",
        "<DEL:ME>": "mobile_element_deletion",
        "<DEL:ME:ALU>": "Alu_deletion",
        "<DEL:ME:L1>": "LINE1_deletion",
        "<DEL:ME:LINE1>": "LINE1_deletion",
        "<DEL:ME:SVA>": "SVA_deletion",
        "<DEL:ME:HERV>": "HERV_deletion",
        "<INS>": "insertion",
        "<INS:ME>": "mobile_element_insertion",
        "<INS:ME:ALU>": "Alu_insertion",
        "<INS:ME:L1>": "LINE1_insertion",
        "<INS:ME:LINE1>": "LINE1_insertion",
        "<INS:ME:SVA>": "SVA_insertion",
        "<INS:ME:HERV>": "HERV_insertion",
        "<INS:ME:LINE>": "mobile_element_insertion",
        "<DUP>": "duplication",
        "<DUP:TANDEM>": "tandem_duplication",
        "<DUP:ISP>": "duplication",
        "<INV>": "inversion",
        "<CNV>": "copy_number_variation",
        "<CNV:TR>": "tandem_repeat",
        "<CN0>": "deletion",
        "<CN=0>": "deletion",
        "<CN1>": "copy_number_variation",
        "<CN2>": "duplication",
        "<CN=2>": "duplication",
        "<CN3>": "copy_number_variation",
        "<CN9>": "copy_number_variation",
        "<CPX>": "CPX",
        "<CPX:DEL>": "CPX",
    }
    # Tokens the mapping must leave exactly as written, and why (see the mapping
    # function's docstring): both engines write the literal, or the Perl allele is
    # synthesised from INFO fields a per-row token does not carry.
    IDENTITY = ["<NON_REF>", "<*>", "<BND>", "<CTX>", "<STR>", "A", "ATAT", "-", "deletion"]

    def setUp(self):
        import tempfile

        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _row(self, allele: str, csq: str = "feature_truncation"):
        return ("v1", "21:1000-2000", allele, "G1", "ENST_A", "Transcript", csq)

    # (a) the headline defect: <DEL> vs deletion

    def test_fastvep_del_matches_perl_deletion_under_the_fastvep_engine(self):
        """A fastVEP `<DEL>` row and a Perl `deletion` row are the SAME annotation
        and must intersect once the scored engine is fastVEP."""
        perl_p, fast_p = self.tmp / "perl.txt", self.tmp / "fastvep.txt"
        write_vep(perl_p, [self._row("deletion")])
        write_vep(fast_p, [self._row("<DEL>")])
        perl, _ = C.parse_vep_output(perl_p)
        fast, _ = C.parse_vep_output(fast_p, engine="fastvep")
        self.assertEqual(perl, fast, "the two rows must normalise to one tuple")
        self.assertEqual(C.compute_metrics([], perl, set(), fast, set()).f1, 1.0)

    def test_fastvep_del_does_not_match_without_the_scored_engine_flag(self):
        """Without `--scored-engine fastvep` the raw token is kept, so the pair is a
        miss. That is the behaviour for every engine other than fastVEP: mapping a
        vep-rs or Perl `<DEL>` would hide a real allele-column bug in whichever
        engine wrote it."""
        perl_p, fast_p = self.tmp / "perl.txt", self.tmp / "fastvep.txt"
        write_vep(perl_p, [self._row("deletion")])
        write_vep(fast_p, [self._row("<DEL>")])
        perl, _ = C.parse_vep_output(perl_p)
        default, _ = C.parse_vep_output(fast_p)
        as_vep_rs, _ = C.parse_vep_output(fast_p, engine="vep-rs")
        self.assertEqual(perl & default, set())
        self.assertEqual(perl & as_vep_rs, set())
        self.assertEqual(C.compute_metrics([], perl, set(), default, set()).f1, 0.0)

    def test_the_perl_reference_is_never_mapped(self):
        """`normalize_allele` with the default or `perl` engine returns the token as
        written. The mapping is a claim about fastVEP's OUTPUT FORMAT, and applying
        it to the reference would make a Perl row that wrongly carried a raw token
        match a correct row."""
        for engine in ("vep-rs", "perl", "VEP-RS"):
            for tok in self.MAPPING:
                self.assertEqual(C.normalize_allele(tok, engine), tok, msg=(engine, tok))

    # (b) each mapped token, one test per row family so a failure names its token

    def _assert_maps(self, token: str, expected: str) -> None:
        self.assertEqual(C.fastvep_symbolic_allele_to_so_class(token), expected)
        self.assertEqual(C.normalize_allele(token, "fastvep"), expected)

    def test_maps_del(self):
        self._assert_maps("<DEL>", "deletion")

    def test_maps_ins(self):
        self._assert_maps("<INS>", "insertion")

    def test_maps_dup(self):
        self._assert_maps("<DUP>", "duplication")

    def test_maps_dup_tandem_before_generic_dup(self):
        """`<DUP:TANDEM>` starts with `<DUP`; the tandem rule must win."""
        self._assert_maps("<DUP:TANDEM>", "tandem_duplication")

    def test_maps_dup_subtype_other_than_tandem_to_duplication(self):
        """`<DUP:ISP>` (in the corpus, 50 records) has no SO term of its own: Perl's
        `get_SO_term` strips the `:ISP` and lands on DUP."""
        self._assert_maps("<DUP:ISP>", "duplication")

    def test_maps_inv(self):
        self._assert_maps("<INV>", "inversion")

    def test_maps_cnv(self):
        self._assert_maps("<CNV>", "copy_number_variation")

    def test_maps_cnv_tr_before_generic_cnv(self):
        """`<CNV:TR>` starts with `<CNV`; the tandem-repeat rule must win, and it must
        produce the exact marker `filter_cnv_tr_expansion_divergences` keys on."""
        self._assert_maps("<CNV:TR>", C._CNV_TR_RUST_ALLELE)
        self.assertEqual(C._CNV_TR_RUST_ALLELE, "tandem_repeat")

    def test_maps_cn0_to_deletion(self):
        """Perl: `$abbrev = "DEL" if $type =~ /^<?CN=?0>?$/`."""
        self._assert_maps("<CN0>", "deletion")
        self._assert_maps("<CN=0>", "deletion")

    def test_maps_cn2_to_duplication(self):
        """Perl: `$abbrev = "DUP" if $type =~ /^<?CN=?2>?$/`."""
        self._assert_maps("<CN2>", "duplication")
        self._assert_maps("<CN=2>", "duplication")

    def test_maps_other_copy_numbers_to_the_generic_class(self):
        for n in (1, 3, 4, 5, 6, 7, 8, 9, 12):
            self._assert_maps(f"<CN{n}>", "copy_number_variation")
            self._assert_maps(f"<CN={n}>", "copy_number_variation")

    def test_unparseable_copy_number_falls_through_to_the_generic_class(self):
        """vep-rs: a `parse::<u32>()` failure on `<CNx>` falls through to
        CopyNumberVariation rather than erroring; mirror it."""
        self._assert_maps("<CNX>", "copy_number_variation")
        self._assert_maps("<CN>", "copy_number_variation")
        self._assert_maps("<CN=X>", "copy_number_variation")

    def test_maps_generic_mobile_element_insertion(self):
        self._assert_maps("<INS:ME>", "mobile_element_insertion")

    def test_maps_alu_insertion(self):
        self._assert_maps("<INS:ME:ALU>", "Alu_insertion")

    def test_maps_line1_insertion_from_both_spellings(self):
        """Perl aliases `L1` to `LINE1` before the subtype lookup."""
        self._assert_maps("<INS:ME:L1>", "LINE1_insertion")
        self._assert_maps("<INS:ME:LINE1>", "LINE1_insertion")

    def test_maps_sva_insertion(self):
        self._assert_maps("<INS:ME:SVA>", "SVA_insertion")

    def test_maps_herv_insertion(self):
        self._assert_maps("<INS:ME:HERV>", "HERV_insertion")

    def test_unknown_mobile_element_subtype_is_the_generic_class(self):
        """`<INS:ME:LINE>` (in the corpus, 200 records) is NOT `LINE1`: Perl's
        `@mobile_elements` grep is exact, so the subtype stays `ME`."""
        self._assert_maps("<INS:ME:LINE>", "mobile_element_insertion")
        self._assert_maps("<INS:ME:UNKNOWN>", "mobile_element_insertion")

    def test_maps_mobile_element_deletion_before_generic_del(self):
        """`<DEL:ME>` starts with `<DEL`; the mobile-element rule must win."""
        self._assert_maps("<DEL:ME>", "mobile_element_deletion")
        self._assert_maps("<DEL:ME:ALU>", "Alu_deletion")
        self._assert_maps("<DEL:ME:SVA>", "SVA_deletion")

    def test_maps_cpx_to_the_bracket_stripped_abbreviation(self):
        """Perl has no SO term for CPX and writes the stripped type; vep-rs's
        ComplexStructural display rule does the same, dropping any `:subtype`."""
        self._assert_maps("<CPX>", "CPX")
        self._assert_maps("<CPX:DEL>", "CPX")

    def test_the_whole_mapping_table_is_pinned(self):
        """Every row at once, so the table in the function docstring, this test and the
        implementation cannot drift apart silently."""
        got = {tok: C.fastvep_symbolic_allele_to_so_class(tok) for tok in self.MAPPING}
        self.assertEqual(got, self.MAPPING)

    def test_mapped_values_are_exactly_the_perl_so_terms(self):
        """The value set must be Perl's `%SO_TERMS` values (plus the ME subtype names
        and CPX), never a raw token: a mapped string that Perl never writes is a
        guaranteed miss dressed up as normalisation."""
        perl_so_terms = {
            "insertion",
            "mobile_element_insertion",
            "Alu_insertion",
            "HERV_insertion",
            "LINE1_insertion",
            "SVA_insertion",
            "deletion",
            "mobile_element_deletion",
            "Alu_deletion",
            "HERV_deletion",
            "LINE1_deletion",
            "SVA_deletion",
            "tandem_repeat",
            "tandem_duplication",
            "duplication",
            "copy_number_variation",
            "inversion",
            "CPX",
        }
        self.assertEqual(set(self.MAPPING.values()), perl_so_terms)

    # tokens that must pass through unchanged

    def test_identity_tokens_are_unchanged(self):
        for tok in self.IDENTITY:
            self.assertEqual(C.fastvep_symbolic_allele_to_so_class(tok), tok, msg=tok)

    def test_non_ref_and_star_are_literal_in_both_engines(self):
        """Perl has no SO term for either and vep-rs `display_allele` returns the raw
        allele, so identity IS the match; mapping them would break it."""
        self._assert_maps("<NON_REF>", "<NON_REF>")
        self._assert_maps("<*>", "<*>")

    def test_symbolic_bnd_is_not_mapped(self):
        """Perl builds `N[chr:pos[` and `N.` alleles for `<BND>` from CHR2/END2/SVLEN.
        No per-row token carries those, so a mapping here could only invent a string
        Perl never writes."""
        self._assert_maps("<BND>", "<BND>")

    def test_an_unknown_symbolic_token_is_left_alone_not_guessed(self):
        """vep-rs classifies an unrecognised token by INFO/SVTYPE, which the comparator
        cannot see. A miss must stay a miss."""
        self._assert_maps("<CTX>", "<CTX>")
        self._assert_maps("<STR>", "<STR>")
        self._assert_maps("<FOO:BAR>", "<FOO:BAR>")

    def test_mapping_is_idempotent_on_already_normalised_values(self):
        """Feeding an SO class string back through the mapping must be a no-op, so the
        order of normalisation steps can never matter."""
        for value in set(self.MAPPING.values()):
            self.assertEqual(C.fastvep_symbolic_allele_to_so_class(value), value)

    def test_token_matching_is_case_insensitive(self):
        """vep-rs upper-cases before classifying and Perl's regexes carry `/i`."""
        self._assert_maps("<del>", "deletion")
        self._assert_maps("<Dup:Tandem>", "tandem_duplication")
        self._assert_maps("<ins:me:alu>", "Alu_insertion")
        self._assert_maps("<cnv:tr>", "tandem_repeat")

    def test_engine_name_matching_is_case_insensitive(self):
        for engine in ("fastvep", "fastVEP", "FASTVEP", " fastvep "):
            self.assertEqual(C.normalize_allele("<DEL>", engine), "deletion", msg=engine)

    def test_breakend_bracket_normalisation_still_applies_under_fastvep(self):
        """The symbolic-token mapping is IN ADDITION to `_strip_chr_in_bracket`, not
        instead of it."""
        self.assertEqual(C.normalize_allele("N[chr21:1234[", "fastvep"), "N[21:1234[")
        self.assertEqual(C.normalize_allele("N[chr21:1234[", "vep-rs"), "N[21:1234[")

    def test_literal_bases_are_untouched_under_fastvep(self):
        for allele in ("A", "ACGT", "-", "ATATATATATATATAT"):
            self.assertEqual(C.normalize_allele(allele, "fastvep"), allele)


class FastvepEngineEndToEndTests(unittest.TestCase):
    """`main()` must hand `--scored-engine` to the scored side's parser and ONLY that
    side. The unit tests above cannot see whether main() threads the flag through, and
    a comparator that mapped nothing would still exit 0 with a plausible report."""

    def setUp(self):
        import tempfile

        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        for sub in ("in", "perl", "rust", "out"):
            (self.tmp / sub).mkdir()
        with open(self.tmp / "in" / "05_symbolic_del_ins.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("21\t1000\t.\tG\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000\n")
            fh.write("21\t3000\t.\tG\t<INS:ME:ALU>\t.\tPASS\tSVTYPE=INS\n")
        def row(allele: str, loc: str, feature: str, csq: str):
            return ("v", loc, allele, "G", feature, "Transcript", csq)

        rows_perl = [
            row("deletion", "21:1001-2000", "ENST_A", "feature_truncation"),
            row("Alu_insertion", "21:3001", "ENST_B", "intron_variant"),
        ]
        rows_fast = [
            row("<DEL>", "21:1001-2000", "ENST_A", "feature_truncation"),
            row("<INS:ME:ALU>", "21:3001", "ENST_B", "intron_variant"),
        ]
        write_vep(self.tmp / "perl" / "05_symbolic_del_ins.txt", rows_perl)
        write_vep(self.tmp / "rust" / "05_symbolic_del_ins.txt", rows_fast)

    def tearDown(self):
        self._tmp.cleanup()

    def _run_main(self, *extra: str) -> dict:
        import json

        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py",
            "--assembly",
            "grch38",
            "--input-dir",
            str(self.tmp / "in"),
            "--perl-dir",
            str(self.tmp / "perl"),
            "--rust-dir",
            str(self.tmp / "rust"),
            "--output-dir",
            str(self.tmp / "out"),
            *extra,
        ]
        try:
            C.main()
        except SystemExit as exc:
            if exc.code not in (0, None):
                raise
        finally:
            sys.argv = argv
        with open(self.tmp / "out" / "concordance_report.json") as fh:
            return json.load(fh)

    def test_scored_engine_fastvep_makes_raw_tokens_match(self):
        report = self._run_main("--scored-engine", "fastvep")
        self.assertEqual(report["overall"]["intersection"], 2)
        self.assertEqual(report["overall"]["f1"], 1.0)
        self.assertEqual(report["adjusted"]["scored_engine"], "fastvep")

    def test_default_engine_leaves_raw_tokens_as_misses(self):
        """The same two files without the flag score 0, the correct outcome for a
        vep-rs output that wrongly carried raw tokens."""
        report = self._run_main()
        self.assertEqual(report["overall"]["intersection"], 0)
        self.assertEqual(report["overall"]["f1"], 0.0)

    def test_the_perl_side_is_parsed_without_the_mapping(self):
        """A Perl reference that itself carried a raw token must NOT be normalised into
        a match: swap the two directories and the mapped side is now the reference."""
        report = self._run_main(
            "--perl-dir",
            str(self.tmp / "rust"),
            "--rust-dir",
            str(self.tmp / "perl"),
            "--scored-engine",
            "fastvep",
        )
        self.assertEqual(report["overall"]["intersection"], 0)


# <CNV:TR> class total: the denominator of a PARTIAL mask must be measured
#
# `excluded_cnvtr_*` alone reads as a full exclusion, so the comparator emits the
# class total on the same tuple sets the mask ran on.


class CnvTrSwapPopulationTests(unittest.TestCase):
    """The <CNV:TR> literal-versus-symbolic pairing and the representation mask.

    vep-rs reads a record symbolically over the whole reference run (Allele
    `tandem_repeat`); VEP expands RUS x RUC and trims the common prefix, so its tuple
    sits at the run's 3' end: a gain is an insertion of whole units between the run's
    last base and the next (Location `end-(end+1)`), a loss a deletion of whole units
    ending at the run's last base (Allele `-`).
    """

    TR = "tandem_repeat"

    @staticmethod
    def _run(end: int, feature: str, rust_csq: str):
        return T(f"21:{end - 17}-{end}", "tandem_repeat", feature, rust_csq)

    @staticmethod
    def _gain(end: int, feature: str, perl_csq: str, units: str = "CACA"):
        return T(f"21:{end}-{end + 1}", units, feature, perl_csq)

    @staticmethod
    def _loss(end: int, feature: str, perl_csq: str, deleted: int = 4):
        return T(f"21:{end - deleted + 1}-{end}", "-", feature, perl_csq)

    def test_pair_kinds(self):
        run = "21:1001-1018"
        self.assertEqual(C._cnv_tr_pair_kind("21:1018-1019", "CACA", run), "gain")
        self.assertEqual(C._cnv_tr_pair_kind("21:1015-1018", "-", run), "loss")
        # A deletion ending elsewhere, an insertion not at the run's end, or a
        # symbolic allele is not this record's literal reading.
        self.assertIsNone(C._cnv_tr_pair_kind("21:1015-1017", "-", run))
        self.assertIsNone(C._cnv_tr_pair_kind("21:1010-1011", "CACA", run))
        self.assertIsNone(C._cnv_tr_pair_kind("21:1001-1018", "deletion", run))
        self.assertIsNone(C._cnv_tr_pair_kind("22:1018-1019", "CACA", run))

    def test_representation_pairs_are_masked_and_consequence_disagreements_are_not(self):
        # Identical sets: pure location and allele difference.
        self.assertTrue(C._is_cnv_tr_representation_pair("intron_variant", "intron_variant", "gain"))
        # The symbolic reading's direction term alone.
        self.assertTrue(C._is_cnv_tr_representation_pair(
            "3_prime_UTR_variant", "3_prime_UTR_variant,feature_elongation", "gain"))
        self.assertTrue(C._is_cnv_tr_representation_pair(
            "3_prime_UTR_variant", "3_prime_UTR_variant,feature_truncation", "loss"))
        # The literal coding effect against the symbolic region term.
        self.assertTrue(C._is_cnv_tr_representation_pair(
            "inframe_insertion", "coding_sequence_variant,feature_elongation", "gain"))
        self.assertTrue(C._is_cnv_tr_representation_pair(
            "frameshift_variant", "coding_sequence_variant,feature_truncation", "loss"))
        # The wrong direction term, a start term, or a splice term is a disagreement.
        self.assertFalse(C._is_cnv_tr_representation_pair(
            "intron_variant", "intron_variant,feature_truncation", "gain"))
        self.assertFalse(C._is_cnv_tr_representation_pair(
            "coding_sequence_variant,start_lost,start_retained_variant",
            "coding_sequence_variant,feature_truncation", "loss"))
        self.assertFalse(C._is_cnv_tr_representation_pair(
            "splice_region_variant,intron_variant", "intron_variant,feature_elongation", "gain"))
        self.assertFalse(C._is_cnv_tr_representation_pair(
            "inframe_insertion", "feature_elongation", "gain"))

    def test_total_counts_every_class_tuple_not_only_the_masked_pairs(self):
        """Three pairs in the class; two differ only in representation and are masked,
        one carries a start term on the literal side and stays charged."""
        masked_gain = (self._gain(2000, "ENST_A", "intron_variant"), self._run(2000, "ENST_A", "intron_variant"))
        masked_loss = (self._loss(3000, "ENST_B", "inframe_deletion"),
                       self._run(3000, "ENST_B", "coding_sequence_variant,feature_truncation"))
        charged = (self._loss(4000, "ENST_C", "coding_sequence_variant,start_lost,start_retained_variant"),
                   self._run(4000, "ENST_C", "coding_sequence_variant,feature_truncation"))
        perl = {masked_gain[0], masked_loss[0], charged[0]}
        rust = {masked_gain[1], masked_loss[1], charged[1]}
        self.assertEqual(C.count_cnv_tr_swap_population(perl, rust), (3, 3))
        ex_rust, ex_perl = C.filter_cnv_tr_expansion_divergences(perl, rust)
        self.assertEqual((len(ex_rust), len(ex_perl)), (2, 2))
        self.assertNotIn(charged[0], ex_perl)
        self.assertNotIn(charged[1], ex_rust)

    def test_total_is_zero_when_the_corpus_has_no_tandem_repeats(self):
        perl = {T("21:1000-2000", "deletion", "ENST_A", "feature_truncation")}
        rust = {T("21:1000-2000", "deletion", "ENST_B", "feature_truncation")}
        self.assertEqual(C.count_cnv_tr_swap_population(perl, rust), (0, 0))

    def test_a_concordant_tandem_repeat_tuple_is_not_a_swap(self):
        """A `tandem_repeat` tuple both engines emitted is agreement, not a class
        member: the population is one-sided tuples only."""
        shared = T("21:1000-2000", self.TR, "ENST_A", "feature_elongation")
        self.assertEqual(C.count_cnv_tr_swap_population({shared}, {shared}), (0, 0))

    def test_a_symbolic_tuple_with_no_literal_partner_counts_on_its_side_only(self):
        """The vep-rs side is defined by the allele; the Perl side has no allele
        signature of its own (a loss is `-`), so it is counted through the pairing. A
        Perl tuple away from the run's end is not this record's literal reading."""
        perl = {T("21:1001-1004", "-", "ENST_A", "intron_variant")}
        rust = {self._run(2000, "ENST_A", "intron_variant")}
        self.assertEqual(C.count_cnv_tr_swap_population(perl, rust), (1, 0))
        self.assertEqual(C.filter_cnv_tr_expansion_divergences(perl, rust), (set(), set()))

    def test_excluded_is_a_subset_of_total_on_each_side(self):
        pairs = [
            (self._gain(1000 * i, f"ENST_{i}", "intron_variant"), self._run(1000 * i, f"ENST_{i}", "intron_variant"))
            for i in range(2, 7)
        ] + [(self._loss(9000, "ENST_X", "start_lost,intron_variant"), self._run(9000, "ENST_X", "intron_variant"))]
        perl = {p for p, _ in pairs}
        rust = {r for _, r in pairs}
        tot_rust, tot_perl = C.count_cnv_tr_swap_population(perl, rust)
        ex_rust, ex_perl = C.filter_cnv_tr_expansion_divergences(perl, rust)
        self.assertLessEqual(len(ex_rust), tot_rust)
        self.assertLessEqual(len(ex_perl), tot_perl)
        self.assertEqual((tot_rust, tot_perl, len(ex_rust), len(ex_perl)), (6, 6, 5, 5))

    def test_the_counter_is_not_a_filter(self):
        """It excludes nothing, so it must not carry the `filter_*_divergences` name
        shape that marks the masks."""
        self.assertFalse(
            [n for n in dir(C) if n.startswith("filter_") and "swap_population" in n]
        )


class TranscriptSelectionOwnCountEndToEndTests(unittest.TestCase):
    """The report states the transcript-selection layer's own size beside the union.

    ``excluded_rust_extra`` is the union of the two vep-rs-side layers, so on its own
    it cannot say how large the transcript-selection class is. The fixture puts one
    file under each layer: a giant deletion on which vep-rs names three transcripts
    Perl lacks, and a ``<CNV:TR>`` gain pair differing only in representation.
    """

    SHARED = (
        "v0", "21:100-200", "deletion", "G", "ENST_S", "Transcript", "feature_truncation"
    )

    @staticmethod
    def _row(uploaded: str, loc: str, allele: str, feature: str, csq: str):
        return (uploaded, loc, allele, "G", feature, "Transcript", csq)

    def setUp(self):
        import json
        import tempfile

        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        for sub in ("in", "perl", "rust", "out"):
            (self.tmp / sub).mkdir()
        with open(self.tmp / "in" / "07_cnv_repeat.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("21\t1984\t.\tA\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=16\n")
        with open(self.tmp / "in" / "09_breakends.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("21\t1000000\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=21000000\n")
        giant = "21:1000000-21000000"
        write_vep(
            self.tmp / "perl" / "07_cnv_repeat.txt",
            [self.SHARED, self._row("v1", "21:2000-2001", "ATAT", "ENST_A", "intron_variant")],
        )
        write_vep(
            self.tmp / "rust" / "07_cnv_repeat.txt",
            [
                self.SHARED,
                self._row("v1", "21:1985-2000", "tandem_repeat", "ENST_A", "intron_variant,feature_elongation"),
            ],
        )
        write_vep(
            self.tmp / "perl" / "09_breakends.txt",
            [self._row("g1", giant, "deletion", "ENST_P", "feature_truncation")],
        )
        write_vep(
            self.tmp / "rust" / "09_breakends.txt",
            [self._row("g1", giant, "deletion", f"ENST_{c}", "feature_truncation") for c in "PQRS"],
        )
        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py",
            "--assembly",
            "grch37",
            "--input-dir",
            str(self.tmp / "in"),
            "--perl-dir",
            str(self.tmp / "perl"),
            "--rust-dir",
            str(self.tmp / "rust"),
            "--output-dir",
            str(self.tmp / "out"),
        ]
        try:
            C.main()
        except SystemExit as exc:
            if exc.code not in (0, None):
                raise
        finally:
            sys.argv = argv
        with open(self.tmp / "out" / "concordance_report.json") as fh:
            self.report = json.load(fh)

    def tearDown(self):
        self._tmp.cleanup()

    def test_own_count_sits_beside_the_union(self):
        adj = self.report["adjusted"]
        self.assertEqual(adj["excluded_transcript_selection_rust"], 3)
        self.assertEqual(adj["excluded_transcript_selection_perl"], 0)
        self.assertEqual(adj["excluded_cnvtr_rust"], 1)
        self.assertEqual(adj["transcript_selection_cnvtr_overlap_rust"], 0)
        self.assertEqual(
            adj["excluded_rust_extra"],
            adj["excluded_transcript_selection_rust"]
            + adj["excluded_cnvtr_rust"]
            - adj["transcript_selection_cnvtr_overlap_rust"],
        )

    def test_per_file_split_names_only_the_files_the_layer_touched(self):
        self.assertEqual(
            self.report["adjusted"]["excluded_transcript_selection_by_file"],
            {"09_breakends": {"rust": 3, "perl": 0}},
        )

    def test_own_count_keys_are_present_when_the_layer_excludes_nothing(self):
        import json

        write_vep(self.tmp / "perl" / "09_breakends.txt", [self.SHARED])
        write_vep(self.tmp / "rust" / "09_breakends.txt", [self.SHARED])
        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py",
            "--assembly",
            "grch37",
            "--input-dir",
            str(self.tmp / "in"),
            "--perl-dir",
            str(self.tmp / "perl"),
            "--rust-dir",
            str(self.tmp / "rust"),
            "--output-dir",
            str(self.tmp / "out"),
        ]
        try:
            C.main()
        except SystemExit as exc:
            if exc.code not in (0, None):
                raise
        finally:
            sys.argv = argv
        with open(self.tmp / "out" / "concordance_report.json") as fh:
            adj = json.load(fh)["adjusted"]
        self.assertEqual(adj["excluded_transcript_selection_rust"], 0)
        self.assertEqual(adj["excluded_transcript_selection_perl"], 0)
        self.assertEqual(adj["excluded_transcript_selection_by_file"], {})


class CnvTrTotalReportEndToEndTests(unittest.TestCase):
    """The JSON report must carry the class total beside the excluded count."""

    SHARED = (
        "v0", "21:100-200", "deletion", "G", "ENST_S", "Transcript", "feature_truncation"
    )

    @staticmethod
    def _row(uploaded: str, loc: str, allele: str, feature: str, csq: str):
        return (uploaded, loc, allele, "G", feature, "Transcript", csq)

    def setUp(self):
        import tempfile

        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        for sub in ("in", "perl", "rust", "out"):
            (self.tmp / sub).mkdir()
        with open(self.tmp / "in" / "07_cnv_repeat.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("21\t1984\t.\tA\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=16\n")
            fh.write("21\t2984\t.\tA\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=16\n")
        # Pair 1: a gain differing only in representation, masked. Pair 2: a loss
        # whose literal reading carries a start term vep-rs lacks, in the class, not
        # masked. Plus one unrelated concordant tuple so F1 is neither 0 nor 1.
        tr = "tandem_repeat"
        write_vep(
            self.tmp / "perl" / "07_cnv_repeat.txt",
            [
                self.SHARED,
                self._row("v1", "21:2000-2001", "ATAT", "ENST_A", "intron_variant"),
                self._row("v2", "21:2997-3000", "-", "ENST_B", "intron_variant,start_lost"),
            ],
        )
        write_vep(
            self.tmp / "rust" / "07_cnv_repeat.txt",
            [
                self.SHARED,
                self._row("v1", "21:1985-2000", tr, "ENST_A", "intron_variant,feature_elongation"),
                self._row("v2", "21:2985-3000", tr, "ENST_B", "intron_variant,feature_truncation"),
            ],
        )
        import json

        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py",
            "--assembly",
            "grch37",
            "--input-dir",
            str(self.tmp / "in"),
            "--perl-dir",
            str(self.tmp / "perl"),
            "--rust-dir",
            str(self.tmp / "rust"),
            "--output-dir",
            str(self.tmp / "out"),
        ]
        try:
            C.main()
        except SystemExit as exc:
            if exc.code not in (0, None):
                raise
        finally:
            sys.argv = argv
        with open(self.tmp / "out" / "concordance_report.json") as fh:
            self.report = json.load(fh)
        self.markdown = (self.tmp / "out" / "concordance_report.md").read_text()

    def tearDown(self):
        self._tmp.cleanup()

    def test_json_carries_the_class_total_beside_the_excluded_count(self):
        adj = self.report["adjusted"]
        self.assertEqual(adj["excluded_cnvtr_rust"], 1)
        self.assertEqual(adj["excluded_cnvtr_perl"], 1)
        self.assertEqual(adj["cnvtr_swap_pairs_total_rust"], 2)
        self.assertEqual(adj["cnvtr_swap_pairs_total_perl"], 2)

    def test_json_carries_the_per_file_split(self):
        self.assertEqual(
            self.report["adjusted"]["cnvtr_swap_pairs_total_by_file"],
            {"07_cnv_repeat": {"rust": 2, "perl": 2}},
        )

    def test_markdown_caption_states_n_of_m(self):
        """A caption carrying only the numerator reads as a full exclusion."""
        self.assertIn("(1 of 2 pairs GRCH37, 07_cnv_repeat only)", self.markdown)

    def test_total_keys_are_present_even_when_zero(self):
        """A corpus with no <CNV:TR> records must state a zero total, not omit the key,
        so a reader cannot mistake absence for a class that was never measured."""
        write_vep(self.tmp / "perl" / "07_cnv_repeat.txt", [self.SHARED])
        write_vep(self.tmp / "rust" / "07_cnv_repeat.txt", [self.SHARED])
        import json

        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py",
            "--assembly",
            "grch37",
            "--input-dir",
            str(self.tmp / "in"),
            "--perl-dir",
            str(self.tmp / "perl"),
            "--rust-dir",
            str(self.tmp / "rust"),
            "--output-dir",
            str(self.tmp / "out"),
        ]
        try:
            C.main()
        except SystemExit as exc:
            if exc.code not in (0, None):
                raise
        finally:
            sys.argv = argv
        with open(self.tmp / "out" / "concordance_report.json") as fh:
            adj = json.load(fh)["adjusted"]
        self.assertEqual(adj["cnvtr_swap_pairs_total_rust"], 0)
        self.assertEqual(adj["cnvtr_swap_pairs_total_perl"], 0)
        self.assertEqual(adj["cnvtr_swap_pairs_total_by_file"], {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
