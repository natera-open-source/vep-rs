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
import tempfile
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

    # Spans that straddle Perl's --max_sv_size default. The GIANT span is twice
    # the default and the SMALL span is 1 kb, so neither exercises the boundary.
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

        A key on the start alone merges their transcript sets, so a transcript
        named on one can suppress a divergence on the other.
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
    similar coordinates, so nearly every transcript it names resolves to chr22. Unlike
    every other mask here the criterion is objectively checkable rather than mechanistic,
    and the authority is the vep-rs JSON cache both engines were run against.
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
        chromosome. Raising there would leave the engine unscored, so the tuple is
        counted and stays in raw AND adjusted; only the Perl side is masked.
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

    def test_there_is_no_consequence_divergence_filter(self) -> None:
        self.assertFalse(
            hasattr(C, "filter_consequence_divergences"),
            "a BND non-coding normalisation mask would exclude pairs on which the two "
            "engines agree; the `filter_*_divergences` name shape marks a real mask",
        )

    def test_there_is_no_normalization_predicate(self) -> None:
        self.assertFalse(
            hasattr(C, "_is_nc_normalization_swap"),
            "the predicate has no caller without the filter",
        )

    def test_the_sv_filter_set_is_exactly_the_expected_one(self) -> None:
        """The 115.2 masks: two full masks (transcript selection, cross-chromosome
        annotation), each with a second arm on the same class (the mate allele of a giant
        breakend; the vep-rs intergenic row of a cross-chromosome record), one partial
        (<CNV:TR>), the pair filter that applies the SNP/indel registry's excluding rules,
        the gVCF reference-block filter, and the mate-side local-read pair filter. The
        116.2 masks: the record the reference skipped, the two breakend shapes of the
        chromosome-aware overlap, the reference-allele rows of a bare breakpoint, the gVCF
        variant site, the dot-breakend pair, the annotated missing ALT, the non-integer
        coordinate, the empty tandem-repeat expansion, the `*`-only ALT under a true SVTYPE,
        the numified POS and the dropped `RN=0` allele.

        Pinned as a set rather than a count, so adding or removing a filter fails
        here with its name, beside the taxonomy row it would need."""
        self.assertEqual(
            {n for n in dir(C) if n.startswith("filter_")},
            {
                "filter_intended_divergences",
                "filter_giant_breakend_mate_divergences",
                "filter_cnv_tr_expansion_divergences",
                "filter_cross_chromosome_divergences",
                "filter_cross_chromosome_orphan_intergenic",
                "filter_registry_swap_pairs",
                "filter_non_ref_batch_divergences",
                "filter_breakend_mate_local_read_pairs",
                "filter_reference_skipped_records",
                "filter_bnd_own_rows_lost",
                "filter_bnd_synonym_mate_rows",
                "filter_bnd_bare_reference_allele_rows",
                "filter_gvcf_variant_sites",
                "filter_dot_breakend_pairs",
                "filter_missing_alt_annotated",
                "filter_non_integer_end_rows",
                "filter_cnv_tr_empty_expansion",
                "filter_star_alone_svtype_annotated",
                "filter_non_numeric_pos_numified",
                "filter_cnv_tr_rn_zero_allele_dropped",
            },
        )

    def test_every_filter_belongs_to_a_release(self) -> None:
        """Every `filter_*` mask is applied against at least one reference release, and
        every name a release lists is a mask; the 115.2 list is the published
        comparator's eight in its order, the 116.2 list drops the class-C arms and the
        two class-D masks and adds the 116.2 masks."""
        filters = {n for n in dir(C) if n.startswith("filter_")}
        listed = {m for masks in C.MASKS_BY_RELEASE.values() for m in masks}
        self.assertEqual(listed, filters)
        self.assertEqual(
            C.MASKS_BY_RELEASE["115.2"],
            (
                "filter_intended_divergences",
                "filter_cnv_tr_expansion_divergences",
                "filter_cross_chromosome_divergences",
                "filter_cross_chromosome_orphan_intergenic",
                "filter_giant_breakend_mate_divergences",
                "filter_registry_swap_pairs",
                "filter_non_ref_batch_divergences",
                "filter_breakend_mate_local_read_pairs",
            ),
        )
        retired = set(C.MASKS_BY_RELEASE["115.2"]) - set(C.MASKS_BY_RELEASE["116.2"])
        self.assertEqual(retired, set(C.MASKS_RETIRED_IN_116_2), "every retired mask carries its reason")
        self.assertEqual(
            set(C.MASKS_BY_RELEASE["116.2"]) - set(C.MASKS_BY_RELEASE["115.2"]),
            {"filter_reference_skipped_records", "filter_bnd_own_rows_lost", "filter_bnd_synonym_mate_rows",
             "filter_bnd_bare_reference_allele_rows", "filter_gvcf_variant_sites", "filter_dot_breakend_pairs",
             "filter_missing_alt_annotated", "filter_non_integer_end_rows", "filter_cnv_tr_empty_expansion",
             "filter_star_alone_svtype_annotated", "filter_non_numeric_pos_numified",
             "filter_cnv_tr_rn_zero_allele_dropped"},
        )

    def test_the_release_block_describes_every_116_2_rule_with_this_comparators_verdict(self) -> None:
        """Every one-sided rule of the 116.2 registry has a count in this comparator's report
        and carries its own cache-authority verdict; a rule added to the registry without a
        mask here is refused at report time, by name."""
        from unittest import mock

        block = C.release_report_block(C.MaskResult(release="116.2"), {}, True)
        classes = block["divergence_classes"]
        self.assertEqual(list(classes), [r.bucket for r in C.REFERENCE_116_2_ONE_SIDED_RULES])
        for rule in C.REFERENCE_116_2_ONE_SIDED_RULES:
            self.assertEqual(classes[rule.bucket]["excludes"], rule.excludes_with_cache, rule.bucket)
            self.assertIn("excluded", classes[rule.bucket])
        self.assertEqual([c["excludes"] for c in classes.values()],
                         [True, True, True, False, True, False, True, True, True, True, True, True, True, True])
        self.assertEqual(block["bnd_bare_reference_allele_rows"], 0)
        self.assertEqual(block["bnd_mate_32bit_buffer_records"], 0)
        self.assertEqual(block["excluded_gvcf_variant_site_rust"], 0)
        with_rows = C.release_report_block(C.MaskResult(release="116.2", bare_reference_rows=3), {}, True)
        self.assertEqual(with_rows["bnd_bare_reference_allele_rows"], 3)
        self.assertEqual(with_rows["divergence_classes"]["bnd_bare_reference_allele_row"]["excluded"], 3)
        sixth = C.REFERENCE_116_2_ONE_SIDED_RULES[0].__class__(
            bucket="a_sixth_shape", definition="A shape.", excludes=False, taxonomy_class="x", perl_citation="X.pm:1"
        )
        with mock.patch.object(C, "REFERENCE_116_2_ONE_SIDED_RULES", C.REFERENCE_116_2_ONE_SIDED_RULES + (sixth,)):
            with self.assertRaisesRegex(AssertionError, "a_sixth_shape"):
                C.release_report_block(C.MaskResult(release="116.2"), {}, True)

    def test_the_bare_breakpoint_reference_allele_rows_are_set_aside_and_counted(self) -> None:
        """Ensembl VEP 116.2 writes a bare `<BND>`'s rows twice, ALLELE_NUM 1 for the reference
        allele N of `N/<BND>` and 2 for the ALT (StructuralVariationOverlap.pm:64-77). The
        mask reads the raw rows and returns the reference's second row per feature, the rows
        `parse_vep_output` dropped from the tuple set, so the set-aside count is 2 while the
        tuple set carries 3 `chromosome_breakpoint` tuples; a bracket record's own and mate
        rows, which differ in Allele, count nothing, and a missing file counts zero."""
        import tempfile
        from pathlib import Path

        rows = [
            "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence",
            "bare\t21:30425651\tchromosome_breakpoint\tG1\tENST1\tTranscript\tcoding_sequence_variant",
            "bare\t21:30425651\tchromosome_breakpoint\tG1\tENST1\tTranscript\tcoding_sequence_variant",
            "bare\t21:30425651\tchromosome_breakpoint\tG2\tENST2\tTranscript\tupstream_gene_variant",
            "bare\t21:30425651\tchromosome_breakpoint\tG2\tENST2\tTranscript\tupstream_gene_variant",
            "paired\t21:30425651\tN.\tG1\tENST1\tTranscript\tfeature_truncation,coding_sequence_variant",
            "paired\t21:30425651\tN[21:30425700[\tG1\tENST1\tTranscript\tfeature_truncation,coding_sequence_variant",
            "single\t21:30425651\tchromosome_breakpoint\tG3\tENST3\tTranscript\tintron_variant",
        ]
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "perl.txt"
            path.write_text("\n".join(rows) + "\n")
            self.assertEqual(C.filter_bnd_bare_reference_allele_rows([path]), 2)
            tuples, _ = C.parse_vep_output(path)
            self.assertEqual(len([t for t in tuples if t[1] == "chromosome_breakpoint"]), 3)
            self.assertEqual(C.filter_bnd_bare_reference_allele_rows([Path(d) / "absent.txt"]), 0)
            self.assertEqual(C.filter_bnd_bare_reference_allele_rows([path, path]), 4)
            masks = C.apply_masks(tuples, tuples, "116.2", {}, {}, {}, "vep-rs", [], perl_paths=[path])
            self.assertEqual(masks.bare_reference_rows, 2)
            self.assertEqual(masks.all_perl(), set(), "rows, not tuples: nothing leaves the tuple arithmetic")


# Cross-file tuple collisions


class CrossFileAggregationTests(unittest.TestCase):
    def test_identical_tuple_in_two_files_collapses_under_union(self):
        """Demonstrates the aggregation hazard the file-summed counts address.

        Two different input VCFs of one assembly can emit a byte-identical tuple.
        A set union keeps one, while `total_input` sums both, so an aggregate row
        built from both mixes two aggregations.
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
            "--reference-release",
            "115.2",
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
        self.assertEqual(C.classify_variant("N", "<BND>"), "BND")
        self.assertEqual(C.classify_variant("N", "<NON_REF>"), "NON_REF")
        self.assertEqual(C.classify_variant("A", "*"), "Spanning")

    def test_per_type_table_keys_a_ranged_record_by_its_printed_location(self):
        """A `<NON_REF>` record at POS 38589780 with END 38590008 prints as
        `21:38589781-38590008`; keyed by `chrom:POS` alone its tuples reach no type row.
        The join through the record ID puts them in the NON_REF row and types the
        discordant rows."""
        variants = [C.VcfVariant("21", 38589780, "synth_nonref_0188", "A", "<NON_REF>", "NON_REF", "21:38589780")]
        loc = "21:38589781-38590008"
        perl = {T(loc, "<NON_REF>", "ENST_A", "intron_variant")}
        rust = {T(loc, "<NON_REF>", "-", "intergenic_variant", ftype="-")}
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "perl.txt"
            write_vep(out, [("synth_nonref_0188", loc, "<NON_REF>", "-", "ENST_A", "Transcript", "intron_variant")])
            types = C.location_types_from_outputs(variants, out)
        self.assertEqual(types, {"21:38589780": "NON_REF", loc: "NON_REF"})
        without = C.compute_per_type_metrics(variants, perl, {loc}, rust, {loc})["NON_REF"]
        with_map = C.compute_per_type_metrics(variants, perl, {loc}, rust, {loc}, types)["NON_REF"]
        self.assertEqual((without.perl_tuples, without.rust_tuples), (0, 0))
        self.assertEqual((with_map.perl_tuples, with_map.rust_tuples, with_map.intersection), (1, 1, 0))
        rows = C.collect_discordants("10_special_alleles", variants, perl, rust, types)
        self.assertEqual({r.variant_type for r in rows}, {"NON_REF"})
        self.assertEqual({r.variant_type for r in C.collect_discordants("10_special_alleles", variants, perl, rust)}, {"Unknown"})

    def test_a_record_without_an_id_keeps_its_position_key(self):
        variants = [C.VcfVariant("21", 100, ".", "A", "T", "SNV", "21:100")]
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "perl.txt"
            write_vep(out, [(".", "21:100", "T", "-", "ENST_A", "Transcript", "missense_variant")])
            self.assertEqual(C.location_types_from_outputs(variants, out), {"21:100": "SNV"})


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
        """`<DUP:ISP>`, one of the synthetic corpus's categories, has no SO term of its
        own: Perl's `get_SO_term` strips the `:ISP` and lands on DUP."""
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
        """`<INS:ME:LINE>`, one of the synthetic corpus's categories, is NOT `LINE1`:
        Perl's `@mobile_elements` grep is exact, so the subtype stays `ME`."""
        self._assert_maps("<INS:ME:LINE>", "mobile_element_insertion")
        self._assert_maps("<INS:ME:UNKNOWN>", "mobile_element_insertion")

    def test_maps_mobile_element_deletion_before_generic_del(self):
        """`<DEL:ME>` starts with `<DEL`; the mobile-element rule must win."""
        self._assert_maps("<DEL:ME>", "mobile_element_deletion")
        self._assert_maps("<DEL:ME:ALU>", "Alu_deletion")
        self._assert_maps("<DEL:ME:SVA>", "SVA_deletion")

    def test_maps_cpx_to_the_bracket_stripped_abbreviation(self):
        """Perl has no SO term for CPX and keeps the type string, `<CPX>` as written or
        `CPX` from INFO/SVTYPE; `record_alleles` lists both spellings, so the fastVEP token
        maps to the bare abbreviation whatever its `:subtype`."""
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
            "--reference-release",
            "115.2",
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
        # A sequence allele elsewhere on the run, or right after it, is the literal
        # reading of a record the expansion does not trim (a multi-allelic record, a
        # unit that is not the reference's repeat, an `N` unit, a disagreeing SVLEN):
        # the whole-run substitution, a partial one, and the units after the run.
        self.assertEqual(C._cnv_tr_pair_kind("21:1010-1011", "CACA", run), "literal")
        self.assertEqual(C._cnv_tr_pair_kind("21:1001-1018", "CACACACACACACACACACA", run), "literal")
        self.assertEqual(C._cnv_tr_pair_kind("21:1019-1024", "CAC", run), "literal")
        self.assertEqual(C._cnv_tr_pair_kind("21:1001-1018", "NNNN", run), "literal")
        # The reference's own symbolic row at the run, written without a FASTA.
        self.assertEqual(C._cnv_tr_pair_kind("21:1001-1018", "tandem_repeat", run), "unexpanded")
        self.assertIsNone(C._cnv_tr_pair_kind("21:1001-1017", "tandem_repeat", run))
        # A deletion ending elsewhere, a symbolic allele of another class, a sequence
        # allele beyond the base after the run, or another chromosome is not a reading.
        self.assertIsNone(C._cnv_tr_pair_kind("21:1015-1017", "-", run))
        self.assertIsNone(C._cnv_tr_pair_kind("21:1001-1018", "deletion", run))
        self.assertIsNone(C._cnv_tr_pair_kind("21:1020-1021", "CACA", run))
        self.assertIsNone(C._cnv_tr_pair_kind("22:1018-1019", "CACA", run))

    def test_a_reading_never_starts_before_the_run(self):
        """The reference's literal allele starts at the run (`$start` of the record) or after
        it: trimming a shared prefix moves the start forward, never back. An upper-case allele
        whose span begins before the run is another record's, whatever its letters: the raw
        abbreviation Ensembl VEP 115.2 writes for a `<CPX>` record (`CPX`) over a 22 kb span that
        covers an 18-base run on the same transcript, or a sequence allele of an earlier run.
        Neither pairs with the run's symbolic row, and the mask leaves both sides charged."""
        run = "21:27268731-27268748"
        self.assertIsNone(C._cnv_tr_pair_kind("21:27254789-27276862", "CPX", run))
        self.assertIsNone(C._cnv_tr_pair_kind("21:27268730-27268748", "CACACACACACACACACACA", run))
        self.assertEqual(C._cnv_tr_pair_kind("21:27268731-27268748", "CACACACACACACACACACA", run), "literal")
        perl = {T("21:27254789-27276862", "CPX", "ENST_A", "intron_variant")}
        rust = {T(run, "tandem_repeat", "ENST_A", "feature_elongation,intron_variant")}
        self.assertEqual(C.filter_cnv_tr_expansion_divergences(perl, rust), (set(), set()))
        self.assertEqual(C.count_cnv_tr_swap_population(perl, rust), (1, 0))

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
        # A literal substitution reads any coding term where the symbolic reading names
        # the region with either direction term (Ensembl VEP 116.2 with a FASTA writes
        # `AGCAGCAGCAGC protein_altering_variant` and `AGCAGC inframe_deletion` for a
        # two-ALT <CNV:TR>, `NNNN frameshift_variant` for `RUS=.`, and
        # `5_prime_UTR_variant,coding_sequence_variant` for a 300-base unit substitution).
        self.assertTrue(C._is_cnv_tr_representation_pair(
            "protein_altering_variant", "coding_sequence_variant,feature_elongation", "literal"))
        self.assertTrue(C._is_cnv_tr_representation_pair(
            "missense_variant", "coding_sequence_variant,feature_truncation", "literal"))
        self.assertTrue(C._is_cnv_tr_representation_pair(
            "5_prime_UTR_variant,coding_sequence_variant", "5_prime_UTR_variant,coding_sequence_variant", "literal"))
        self.assertFalse(C._is_cnv_tr_representation_pair(
            "splice_region_variant,missense_variant", "coding_sequence_variant,feature_elongation", "literal"))
        # Without a FASTA the reference's symbolic row labels every record a gain: a loss
        # vep-rs reads as the deletion of the run's tail pairs with it on the direction and
        # coding terms alone (Ensembl VEP 116.2 writes `feature_elongation,coding_sequence_variant`
        # for `SVLEN=9;RUS=AGC;RUC=2` without `--fasta`).
        self.assertTrue(C._is_cnv_tr_representation_pair(
            "feature_elongation,coding_sequence_variant", "feature_truncation,inframe_deletion", "unexpanded"))
        self.assertTrue(C._is_cnv_tr_representation_pair(
            "feature_elongation,coding_sequence_variant", "feature_elongation,coding_sequence_variant", "unexpanded"))
        self.assertFalse(C._is_cnv_tr_representation_pair(
            "feature_elongation,coding_sequence_variant", "feature_truncation,intron_variant", "unexpanded"))

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
            "--reference-release",
            "115.2",
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
            "--reference-release",
            "115.2",
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
            "--reference-release",
            "115.2",
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
            "--reference-release",
            "115.2",
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


class AdjustedDiscordantTsvEndToEndTests(unittest.TestCase):
    """``discordant_adjusted.tsv`` lists exactly the one-sided tuples the adjusted
    row counts: every mask layer removes its pairs from it, and nothing else does.

    The corpus puts one tuple under each layer and two outside every layer:

    * ``09_breakends``: a giant deletion on which vep-rs names three transcripts Perl
      lacks (transcript-selection layer, vep-rs side) and Perl names one chr22
      transcript (cross-chromosome layer, Perl side); on a 1 kb deletion below the
      transcript-selection scope Perl names a chr21 transcript vep-rs lacks, which
      stays: a genuine vep-rs gap.
    * ``07_cnv_repeat``: a ``<CNV:TR>`` gain pair differing only in representation
      (the literal-expansion layer, both sides) and a loss pair whose literal reading
      carries ``start_lost``, which stays on both sides.
    """

    TX_BY_CHR = {"21": ["ENST_P", "ENST_Q", "ENST_R", "ENST_S", "ENST_GAP", "ENST_A", "ENST_B"],
                 "22": ["ENST_CHR22"]}
    GIANT = "21:1000000-21000000"
    SMALL = "21:5000-6000"

    @staticmethod
    def _row(uploaded: str, loc: str, allele: str, feature: str, csq: str):
        return (uploaded, loc, allele, "G", feature, "Transcript", csq)

    def setUp(self):
        import json
        import tempfile

        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        for sub in ("in", "perl", "rust", "out", "cache"):
            (self.tmp / sub).mkdir()
        for chrom, ids in self.TX_BY_CHR.items():
            shard = self.tmp / "cache" / "transcripts" / chrom
            shard.mkdir(parents=True)
            (shard / "0-1.json").write_text(json.dumps([{"stable_id": i} for i in ids]))
        with open(self.tmp / "in" / "07_cnv_repeat.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("21\t1984\t.\tA\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=16\n")
            fh.write("21\t2984\t.\tA\t<CNV:TR>\t.\tPASS\tSVTYPE=CNV;SVLEN=16\n")
        with open(self.tmp / "in" / "09_breakends.vcf", "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("21\t1000000\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=21000000\n")
            fh.write("21\t5000\t.\tA\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=6000\n")
        tr = "tandem_repeat"
        write_vep(
            self.tmp / "perl" / "07_cnv_repeat.txt",
            [
                self._row("v1", "21:2000-2001", "ATAT", "ENST_A", "intron_variant"),
                self._row("v2", "21:2997-3000", "-", "ENST_B", "intron_variant,start_lost"),
            ],
        )
        write_vep(
            self.tmp / "rust" / "07_cnv_repeat.txt",
            [
                self._row("v1", "21:1985-2000", tr, "ENST_A", "intron_variant,feature_elongation"),
                self._row("v2", "21:2985-3000", tr, "ENST_B", "intron_variant,feature_truncation"),
            ],
        )
        write_vep(
            self.tmp / "perl" / "09_breakends.txt",
            [
                self._row("g1", self.GIANT, "deletion", "ENST_P", "feature_truncation"),
                self._row("g1", self.GIANT, "deletion", "ENST_CHR22", "feature_truncation"),
                self._row("g2", self.SMALL, "deletion", "ENST_GAP", "feature_truncation"),
            ],
        )
        write_vep(
            self.tmp / "rust" / "09_breakends.txt",
            [self._row("g1", self.GIANT, "deletion", f"ENST_{c}", "feature_truncation") for c in "PQRS"],
        )
        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py",
            "--reference-release",
            "115.2",
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
            "--vep-rs-cache",
            str(self.tmp / "cache"),
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
        self.raw_lines = (self.tmp / "out" / "discordant.tsv").read_text().splitlines()
        self.adj_lines = (self.tmp / "out" / "discordant_adjusted.tsv").read_text().splitlines()

    def tearDown(self):
        self._tmp.cleanup()

    def test_every_layer_fired_on_the_fixture(self):
        """The corpus must exercise all three layers, or the subtraction below is
        vacuous."""
        adj = self.report["adjusted"]
        self.assertEqual(adj["excluded_transcript_selection_rust"], 3)
        self.assertEqual(adj["excluded_cross_chromosome_perl"], 1)
        self.assertEqual((adj["excluded_cnvtr_rust"], adj["excluded_cnvtr_perl"]), (1, 1))

    def test_same_header_as_the_raw_file(self):
        self.assertEqual(self.adj_lines[0], self.raw_lines[0])

    def test_row_count_is_the_adjusted_one_sided_total(self):
        adj = self.report["adjusted"]
        self.assertEqual(len(self.adj_lines) - 1, adj["only_perl"] + adj["only_rust"])
        self.assertEqual(
            len(self.adj_lines) - 1,
            sum(m["only_perl"] + m["only_rust"] for m in self.report["per_file_adjusted"].values()),
        )

    def test_raw_file_still_carries_every_one_sided_tuple(self):
        ov = self.report["overall"]
        self.assertEqual(len(self.raw_lines) - 1, ov["only_perl"] + ov["only_rust"])
        self.assertGreater(len(self.raw_lines), len(self.adj_lines))

    def test_rows_are_a_subset_of_the_raw_rows_in_raw_order(self):
        raw_index = {line: i for i, line in enumerate(self.raw_lines)}
        positions = [raw_index[line] for line in self.adj_lines]
        self.assertEqual(positions, sorted(positions))

    def test_exactly_the_unmasked_tuples_survive(self):
        survivors = {tuple(line.split("\t")[2:6]) for line in self.adj_lines[1:]}
        self.assertEqual(
            survivors,
            {
                ("only_perl", self.SMALL, "deletion", "ENST_GAP"),
                ("only_perl", "21:2997-3000", "-", "ENST_B"),
                ("only_rust", "21:2985-3000", "tandem_repeat", "ENST_B"),
            },
        )
        masked = {("only_rust", self.GIANT, "deletion", f"ENST_{c}") for c in "QRS"}
        masked |= {
            ("only_perl", self.GIANT, "deletion", "ENST_CHR22"),
            ("only_perl", "21:2000-2001", "ATAT", "ENST_A"),
            ("only_rust", "21:1985-2000", "tandem_repeat", "ENST_A"),
        }
        raw_keys = {tuple(line.split("\t")[2:6]) for line in self.raw_lines[1:]}
        self.assertEqual(raw_keys, survivors | masked)

    def test_file_is_written_empty_with_a_header_when_nothing_survives(self):
        write_vep(self.tmp / "perl" / "07_cnv_repeat.txt", [])
        write_vep(self.tmp / "rust" / "07_cnv_repeat.txt", [])
        rows = [self._row("g1", self.GIANT, "deletion", "ENST_P", "feature_truncation")]
        write_vep(self.tmp / "perl" / "09_breakends.txt", rows)
        write_vep(self.tmp / "rust" / "09_breakends.txt", rows)
        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py",
            "--reference-release",
            "115.2",
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
        lines = (self.tmp / "out" / "discordant_adjusted.tsv").read_text().splitlines()
        self.assertEqual(lines, [self.raw_lines[0]])


# The vep-rs side of the cross-chromosome mask


class CrossChromosomeOrphanIntergenicTests(unittest.TestCase):
    """The intergenic row vep-rs writes for a record whose every Perl tuple the
    cross-chromosome mask removed.

    The exemplar is the gnomAD GRCh37 record `21:16155574-16155648 deletion`: Perl named
    ENST00000380289, a chrX transcript numerically overlapping the record, because a
    neighbouring breakend's mate loaded the chrX region; alone it writes
    `intergenic_variant`, which is what vep-rs wrote.
    """

    TX_BY_CHR = {"21": {"ENST_21_A", "ENST_21_B"}, "X": {"ENST00000380289"}, "22": {"ENST_22_X"}}
    LOC = "21:16155574-16155648"

    def _xchr(self, perl, rust):
        _, ex_perl, _ = C.filter_cross_chromosome_divergences(perl, rust, self.TX_BY_CHR)
        return ex_perl

    def test_the_orphan_intergenic_row_is_excluded(self) -> None:
        perl = {T(self.LOC, "deletion", "ENST00000380289", "intron_variant")}
        rust = {T(self.LOC, "deletion", "-", "intergenic_variant", ftype="-")}
        ex = C.filter_cross_chromosome_orphan_intergenic(perl, rust, self._xchr(perl, rust))
        self.assertEqual(ex, rust)

    def test_every_perl_tuple_must_have_been_masked(self) -> None:
        """One surviving on-chromosome Perl transcript makes the vep-rs intergenic row a
        genuine divergence (vep-rs missed a chr21 transcript), so nothing is excluded."""
        perl = {
            T(self.LOC, "deletion", "ENST00000380289", "intron_variant"),
            T(self.LOC, "deletion", "ENST_21_A", "intron_variant"),
        }
        rust = {T(self.LOC, "deletion", "-", "intergenic_variant", ftype="-")}
        ex_perl = self._xchr(perl, rust)
        self.assertEqual(len(ex_perl), 1, "only the chrX tuple is cross-chromosome")
        self.assertEqual(C.filter_cross_chromosome_orphan_intergenic(perl, rust, ex_perl), set())

    def test_a_record_perl_wrote_nothing_for_is_not_touched(self) -> None:
        rust = {T(self.LOC, "deletion", "-", "intergenic_variant", ftype="-")}
        self.assertEqual(C.filter_cross_chromosome_orphan_intergenic(set(), rust, set()), set())

    def test_a_vep_rs_transcript_row_beside_the_intergenic_row_blocks_it(self) -> None:
        """vep-rs's set must be exactly the intergenic row; anything else at the key is a
        shape this arm does not describe and stays charged."""
        perl = {T(self.LOC, "deletion", "ENST00000380289", "intron_variant")}
        rust = {
            T(self.LOC, "deletion", "-", "intergenic_variant", ftype="-"),
            T(self.LOC, "deletion", "ENST_21_B", "intron_variant"),
        }
        self.assertEqual(C.filter_cross_chromosome_orphan_intergenic(perl, rust, self._xchr(perl, rust)), set())

    def test_a_non_intergenic_vep_rs_only_tuple_is_never_excluded_here(self) -> None:
        perl = {T(self.LOC, "deletion", "ENST00000380289", "intron_variant")}
        rust = {T(self.LOC, "deletion", "ENST_21_B", "intron_variant")}
        self.assertEqual(C.filter_cross_chromosome_orphan_intergenic(perl, rust, self._xchr(perl, rust)), set())

    def test_no_cross_chromosome_exclusions_means_no_orphans(self) -> None:
        perl = {T(self.LOC, "deletion", "ENST_21_A", "intron_variant")}
        rust = {T(self.LOC, "deletion", "-", "intergenic_variant", ftype="-")}
        self.assertEqual(C.filter_cross_chromosome_orphan_intergenic(perl, rust, set()), set())


# The mate-allele arm of the transcript-selection mask


class GiantBreakendMateAlleleTests(unittest.TestCase):
    """The exemplar is `gnomAD-SV_v3_BND_chr21_3f3f16e7` (GRCh38): Location
    `21:14604132-41770788`, mate `N[chr21:41770788[`. vep-rs writes 24 mate-allele rows;
    Perl writes none under that allele, only rows under the local allele `N.`. With a
    record loading the mate's region Perl writes exactly the 20 rows whose transcript lies
    within 5 kb of the mate; the other 4 are transcripts near the LOCAL breakend, which
    Perl's distance gate never admits on the mate allele.
    """

    LOC = "21:14604132-41770788"
    MATE = "N[chr21:41770788["
    NEAR = {  # transcript: (start, end, consequence set) within 5 kb of 41770788
        "ENST00000332512": (41739373, 41767052, "coding_transcript_variant"),
        "ENST00000352483": (41739369, 41767089, "coding_transcript_variant"),
        "ENST00000686319": (41767163, 41770302, "non_coding_transcript_variant"),
        "ENST00000700893": (41767163, 41774864, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000701256": (41767176, 41776073, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000701512": (41767167, 41770305, "non_coding_transcript_variant"),
        "ENST00000827190": (41764826, 41776429, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827191": (41765005, 41776073, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827192": (41767164, 41776073, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827193": (41768372, 41776433, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827194": (41767173, 41775164, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827195": (41767815, 41776234, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827196": (41769001, 41776073, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827197": (41767926, 41774864, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827198": (41768819, 41774858, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827199": (41768990, 41774864, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827200": (41769079, 41774864, "feature_truncation,non_coding_transcript_exon_variant,intron_variant"),
        "ENST00000827201": (41767163, 41770305, "non_coding_transcript_variant"),
        "ENST00000827202": (41767816, 41770307, "non_coding_transcript_variant"),
        "ENST00000827413": (41768671, 41769359, "non_coding_transcript_variant"),
    }
    FAR = {  # transcripts near the local breakend at 14.6 Mb: vep-rs rows Perl never writes
        "ENST00000400562": (14593720, 14658821, "intron_variant,non_coding_transcript_variant"),
        "ENST00000435732": (14591930, 14643133, "intron_variant,non_coding_transcript_variant"),
        "ENST00000647101": (14485228, 14658821, "intron_variant"),
        "ENST00000749056": (14585012, 14612700, "intron_variant,non_coding_transcript_variant"),
    }

    def spans(self):
        return {("21", tx): ("21", s, e) for tx, (s, e, _c) in {**self.NEAR, **self.FAR}.items()}

    def rust(self):
        return {T(self.LOC, self.MATE, tx, c) for tx, (_s, _e, c) in {**self.NEAR, **self.FAR}.items()}

    def perl_local(self):
        return {T(self.LOC, "N.", tx, "feature_truncation,non_coding_transcript_variant") for tx in ("ENST00000435315", "ENST00000448152")}

    def test_exactly_the_transcripts_within_5kb_of_the_mate_are_excluded(self) -> None:
        ex = C.filter_giant_breakend_mate_divergences(self.perl_local(), self.rust(), self.spans())
        self.assertEqual({t[2] for t in ex}, set(self.NEAR))
        self.assertEqual(len(ex), 20)

    def test_a_perl_row_under_the_mate_allele_hands_the_record_to_the_transcript_selection_arm(self) -> None:
        """Once Perl names any transcript under the mate allele the transcript-selection arms apply
        and this one must stay out, or the two would double-claim."""
        perl = self.perl_local() | {T(self.LOC, self.MATE, "ENST00000332512", "coding_transcript_variant")}
        self.assertEqual(C.filter_giant_breakend_mate_divergences(perl, self.rust(), self.spans()), set())

    def test_a_record_perl_never_saw_is_not_touched(self) -> None:
        self.assertEqual(C.filter_giant_breakend_mate_divergences(set(), self.rust(), self.spans()), set())

    def test_a_span_within_max_sv_size_is_out_of_scope(self) -> None:
        loc = "21:25729803-25735486"
        rust = {T(loc, "N[chr21:25729802[", "ENST00000400090", "feature_truncation,intron_variant")}
        perl = {T(loc, "N.", "ENST00000400090", "feature_truncation,intron_variant")}
        spans = {("21", "ENST00000400090"): ("21", 25725000, 25740000)}
        self.assertEqual(C.filter_giant_breakend_mate_divergences(perl, rust, spans), set())

    def test_a_mate_on_another_chromosome_than_the_feature_is_out_of_scope(self) -> None:
        rust = {T(self.LOC, "N[chr22:41770788[", "ENST00000332512", "coding_transcript_variant")}
        self.assertEqual(C.filter_giant_breakend_mate_divergences(self.perl_local(), rust, self.spans()), set())

    def test_a_non_bracket_allele_is_out_of_scope(self) -> None:
        rust = {T(self.LOC, "N.", "ENST00000332512", "coding_transcript_variant")}
        self.assertEqual(C.filter_giant_breakend_mate_divergences(self.perl_local(), rust, self.spans()), set())

    def test_no_cache_excludes_nothing(self) -> None:
        self.assertEqual(C.filter_giant_breakend_mate_divergences(self.perl_local(), self.rust(), {}), set())

    def test_bracket_mate_parsing(self) -> None:
        self.assertEqual(C._bracket_mate("N[chr21:41770788["), ("21", 41770788))
        self.assertEqual(C._bracket_mate("[5:37717217[A"), ("5", 37717217))
        self.assertEqual(C._bracket_mate("]13:123456]T"), ("13", 123456))
        self.assertIsNone(C._bracket_mate("N."))
        self.assertIsNone(C._bracket_mate("deletion"))


# The vep_skip ALT-type scope of the transcript-selection mask


class VepSkipAlleleScopeTests(unittest.TestCase):
    """`CPX` and `CTX` reach the transcript-selection arms at any span.

    Perl flags both `vep_skip` for want of a Sequence Ontology term and loads no region
    for them, so their annotation is whatever the batch loaded: alone, every such record
    writes `intergenic_variant`, and with its regions loaded Perl writes exactly the set
    vep-rs writes. The exemplar is the GRCh38 gnomAD record `21:14923959-16588080 CPX`,
    on which vep-rs names more transcripts than Perl.
    """

    def test_a_vep_skip_allele_is_in_scope_at_any_span(self) -> None:
        self.assertTrue(C._batch_dependent_annotation_scope("21:14923959-16588080", "CPX"))
        self.assertTrue(C._batch_dependent_annotation_scope("21:39955356", "CTX"))
        self.assertTrue(C._batch_dependent_annotation_scope("21:1-20000002", "deletion"))
        self.assertFalse(C._batch_dependent_annotation_scope("21:38530770-38531198", "deletion"))
        self.assertFalse(C._batch_dependent_annotation_scope("21:38530770-38531198", "<*>"))

    def test_a_non_ref_allele_is_out_of_scope(self) -> None:
        """A gVCF reference block gets no transcript tuple from vep-rs, so Perl's
        batch-derived transcript tuples on it are not a transcript-set difference this
        arm describes; they belong to a filter of their own."""
        self.assertFalse(C._batch_dependent_annotation_scope("21:38530770-38531198", "<NON_REF>"))
        loc = "21:38530770-38531198"
        perl = {T(loc, "<NON_REF>", "ENST_A", "intron_variant")}
        rust = {T(loc, "<NON_REF>", "-", "intergenic_variant", ftype="-")}
        self.assertEqual(C.filter_intended_divergences(perl, rust), (set(), set()))

    def test_perl_naming_fewer_transcripts_on_a_cpx_record_is_excluded(self) -> None:
        loc = "21:14923959-16588080"
        perl = {T(loc, "CPX", "ENST_A", "intron_variant")}
        rust = {T(loc, "CPX", "ENST_A", "intron_variant"), T(loc, "CPX", "ENST_B", "non_coding_transcript_variant")}
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        self.assertEqual(ex_rust, {T(loc, "CPX", "ENST_B", "non_coding_transcript_variant")})
        self.assertEqual(ex_perl, set())

    def test_perl_intergenic_beside_vep_rs_transcripts_on_a_cpx_record(self) -> None:
        """Perl's batch loaded nothing for the record and it wrote the intergenic row;
        vep-rs named the overlapping transcripts. Both sides are the masked shape."""
        loc = "21:22678893-22680544"
        perl = {T(loc, "CPX", "-", "intergenic_variant", ftype="-")}
        rust = {T(loc, "CPX", "ENST_A", "intron_variant"), T(loc, "CPX", "ENST_B", "upstream_gene_variant")}
        ex_rust, ex_perl = C.filter_intended_divergences(perl, rust)
        self.assertEqual(ex_rust, rust)
        self.assertEqual(ex_perl, perl)

    def test_the_direction_guard_holds_for_a_vep_skip_allele(self) -> None:
        """vep-rs naming FEWER transcripts than Perl on a CPX record is not the defect."""
        loc = "21:14923959-16588080"
        perl = {T(loc, "CPX", "ENST_A", "intron_variant"), T(loc, "CPX", "ENST_B", "intron_variant")}
        rust = {T(loc, "CPX", "ENST_A", "intron_variant")}
        self.assertEqual(C.filter_intended_divergences(perl, rust), (set(), set()))

    def test_a_pair_on_a_vep_skip_record_is_left_alone(self) -> None:
        """Same transcript, different terms: a consequence divergence, not a
        transcript-selection one."""
        loc = "21:37386786-37420291"
        perl = {T(loc, "CPX", "ENST_A", "intron_variant,splice_polypyrimidine_tract_variant")}
        rust = {T(loc, "CPX", "ENST_A", "intron_variant")}
        self.assertEqual(C.filter_intended_divergences(perl, rust), (set(), set()))

    def test_a_supported_small_allele_stays_out_of_scope(self) -> None:
        loc = "21:23699593-23724683"
        perl = {T(loc, "deletion", "ENST_A", "transcript_ablation")}
        rust = {T(loc, "deletion", "ENST_A", "transcript_ablation"), T(loc, "deletion", "ENST_B", "transcript_ablation")}
        self.assertEqual(C.filter_intended_divergences(perl, rust), (set(), set()))


# The SNP/indel registry's excluding rules on this comparator's pairs


class RegistrySwapPairTests(unittest.TestCase):
    """Both members of a pair one of `compare_vep_outputs.EXCLUDING_RULES` matches.

    Exemplars: `21:26965975-26966647 deletion ENST00000284987` (GRCh38 gnomAD), Perl
    `5_prime_UTR_variant,feature_truncation,frameshift_variant,start_lost,start_retained_variant`
    against vep-rs without `start_retained_variant`; `21:36265246 - ENST00000486278`
    (GRCh37 ClinVar SV), Perl `frameshift_variant,start_lost,start_retained_variant`
    against vep-rs without `start_lost`; and `21:31663791-31663794 CA ENST00000270142`
    (GRCh38 `02_mnp_complex`), vep-rs adding `splice_region_variant` to Perl's
    `frameshift_variant`.
    """

    def test_structural_start_co_emission_pair_is_excluded_both_sides(self) -> None:
        loc = "21:26965975-26966647"
        p = T(loc, "deletion", "ENST00000284987", "5_prime_UTR_variant,feature_truncation,frameshift_variant,start_lost,start_retained_variant")
        r = T(loc, "deletion", "ENST00000284987", "5_prime_UTR_variant,feature_truncation,frameshift_variant,start_lost")
        ex_rust, ex_perl, buckets = C.filter_registry_swap_pairs({p}, {r})
        self.assertEqual((ex_rust, ex_perl), ({r}, {p}))
        self.assertEqual(buckets["start_cooccurrence_swap"], 1)
        self.assertEqual(buckets["covered_splice_region_swap"], 0)

    def test_sequence_start_co_emission_pair_drops_the_other_member(self) -> None:
        p = T("21:36265246", "-", "ENST00000486278", "frameshift_variant,start_lost,start_retained_variant")
        r = T("21:36265246", "-", "ENST00000486278", "frameshift_variant,start_retained_variant")
        ex_rust, ex_perl, _ = C.filter_registry_swap_pairs({p}, {r})
        self.assertEqual((ex_rust, ex_perl), ({r}, {p}))

    def test_the_wrong_member_for_the_allele_kind_is_not_excluded(self) -> None:
        """A structural allele on which vep-rs dropped `start_lost` instead is a vep-rs
        divergence, not the co-emission shape (`_is_start_cooccurrence_swap`)."""
        loc = "21:26965975-26966647"
        p = T(loc, "deletion", "ENST00000284987", "5_prime_UTR_variant,feature_truncation,frameshift_variant,start_lost,start_retained_variant")
        r = T(loc, "deletion", "ENST00000284987", "5_prime_UTR_variant,feature_truncation,frameshift_variant,start_retained_variant")
        self.assertEqual(C.filter_registry_swap_pairs({p}, {r})[:2], (set(), set()))

    def test_both_start_terms_missing_on_vep_rs_is_not_the_shape(self) -> None:
        loc = "21:39660415-39660694"
        p = T(loc, "Alu_insertion", "ENST00000852645", "5_prime_UTR_variant,coding_sequence_variant,feature_elongation,start_lost,start_retained_variant")
        r = T(loc, "Alu_insertion", "ENST00000852645", "5_prime_UTR_variant,coding_sequence_variant,feature_elongation")
        self.assertEqual(C.filter_registry_swap_pairs({p}, {r})[:2], (set(), set()))

    def test_covered_splice_region_pair_is_excluded_both_sides(self) -> None:
        loc = "21:31663791-31663794"
        p = T(loc, "CA", "ENST00000270142", "frameshift_variant")
        r = T(loc, "CA", "ENST00000270142", "frameshift_variant,splice_region_variant")
        ex_rust, ex_perl, buckets = C.filter_registry_swap_pairs({p}, {r})
        self.assertEqual((ex_rust, ex_perl), ({r}, {p}))
        self.assertEqual(buckets["covered_splice_region_swap"], 1)

    def test_a_polypyrimidine_tract_pair_is_not_the_covered_shape(self) -> None:
        loc = "21:29079809-29085808"
        p = T(loc, "deletion", "ENST00000341618", "intron_variant,splice_polypyrimidine_tract_variant")
        r = T(loc, "deletion", "ENST00000341618", "intron_variant")
        self.assertEqual(C.filter_registry_swap_pairs({p}, {r})[:2], (set(), set()))

    def test_an_extra_vep_rs_term_beside_the_co_emission_rejects_the_pair(self) -> None:
        loc = "21:26965975-26966647"
        p = T(loc, "deletion", "ENST00000284987", "feature_truncation,start_lost,start_retained_variant")
        r = T(loc, "deletion", "ENST00000284987", "feature_truncation,start_lost,intron_variant")
        self.assertEqual(C.filter_registry_swap_pairs({p}, {r})[:2], (set(), set()))

    def test_one_sided_tuples_are_never_pairs(self) -> None:
        p = T("21:1-2", "deletion", "ENST_A", "start_lost,start_retained_variant")
        r = T("21:1-2", "deletion", "ENST_B", "start_lost")
        self.assertEqual(C.filter_registry_swap_pairs({p}, {r})[:2], (set(), set()))


# Perl's batch-derived transcript tuples on a gVCF reference block


class NonRefBatchTests(unittest.TestCase):
    """`<NON_REF>` records: vep-rs writes one intergenic row, Perl whatever its batch loaded.

    The exemplar is `synth_nonref_0188` (`21:38589781-38590008 <NON_REF>`, GRCh37): alone,
    Perl writes `intergenic_variant`; in its batch it writes transcript rows for the regions
    its neighbours loaded.
    """

    LOC = "21:38589781-38590008"

    def intergenic(self, loc=None):
        return T(loc or self.LOC, "<NON_REF>", "-", "intergenic_variant", ftype="-")

    def test_perl_transcript_tuples_and_the_displaced_intergenic_row_are_excluded(self) -> None:
        perl = {T(self.LOC, "<NON_REF>", "ENST00000440629", "upstream_gene_variant"),
                T(self.LOC, "<NON_REF>", "ENST00000454482", "intron_variant,non_coding_transcript_variant")}
        rust = {self.intergenic()}
        ex_rust, ex_perl = C.filter_non_ref_batch_divergences(perl, rust)
        self.assertEqual(ex_perl, perl)
        self.assertEqual(ex_rust, rust)

    def test_a_matching_intergenic_row_stays_on_both_sides(self) -> None:
        """Perl's batch loaded nothing: both engines wrote the intergenic row and it is
        concordant; nothing to exclude."""
        self.assertEqual(C.filter_non_ref_batch_divergences({self.intergenic()}, {self.intergenic()}), (set(), set()))

    def test_a_vep_rs_transcript_tuple_on_a_non_ref_record_is_left_charged(self) -> None:
        """vep-rs writing anything but the intergenic row is outside the shape whatever
        the Perl side holds, so a synthetic vep-rs transcript tuple stays in both
        denominators."""
        perl = {T(self.LOC, "<NON_REF>", "ENST00000440629", "upstream_gene_variant")}
        rust = {self.intergenic(), T(self.LOC, "<NON_REF>", "ENST00000440629", "intron_variant")}
        self.assertEqual(C.filter_non_ref_batch_divergences(perl, rust), (set(), set()))
        rust2 = {T(self.LOC, "<NON_REF>", "ENST00000440629", "upstream_gene_variant")}
        self.assertEqual(C.filter_non_ref_batch_divergences(perl, rust2), (set(), set()))

    def test_another_allele_is_never_touched(self) -> None:
        perl = {T(self.LOC, "deletion", "ENST00000440629", "upstream_gene_variant")}
        rust = {T(self.LOC, "deletion", "-", "intergenic_variant", ftype="-")}
        self.assertEqual(C.filter_non_ref_batch_divergences(perl, rust), (set(), set()))

    def test_a_record_perl_wrote_nothing_for_is_not_touched(self) -> None:
        self.assertEqual(C.filter_non_ref_batch_divergences(set(), {self.intergenic()}), (set(), set()))


# A point breakend's mate-side row, read by Perl at the local coordinate


class BreakendMateLocalReadTests(unittest.TestCase):
    """The exemplar is `1:14374334 ]21:33034355]A` on ENST00000270142 (GRCh37, chromosome
    21 forward strand, CDS 33032083-33040891): the mate lies in intron 1 and the local
    coordinate on chromosome 1 selects nothing, so Perl writes `feature_truncation` and
    the mate breakend supports `feature_truncation,intron_variant`.
    """

    MODEL = C.TranscriptModel(
        chrom="21", start=33031935, end=33041244, strand=1, biotype="protein_coding",
        exons=((33031935, 33032154), (33036103, 33036199), (33038762, 33038831), (33039571, 33039688), (33040784, 33041244)),
        cds=(33032083, 33040891), mirna_mature=(),
    )
    MODELS = {("21", "ENST00000270142"): MODEL}
    LOC, ALLELE, TX = "1:14374334", "]21:33034355]A", "ENST00000270142"

    def pair(self, perl_csq, rust_csq):
        return {T(self.LOC, self.ALLELE, self.TX, perl_csq)}, {T(self.LOC, self.ALLELE, self.TX, rust_csq)}

    def test_the_exemplar_pair_is_excluded_on_both_sides(self) -> None:
        perl, rust = self.pair("feature_truncation", "feature_truncation,intron_variant")
        ex_rust, ex_perl, kinds = C.filter_breakend_mate_local_read_pairs(perl, rust, self.MODELS)
        self.assertEqual((ex_rust, ex_perl), (rust, perl))
        self.assertEqual(kinds, {"cross_chromosome": 1, "same_chromosome": 0})

    def test_a_vep_rs_over_annotation_is_left_charged(self) -> None:
        perl, rust = self.pair("feature_truncation", "feature_truncation,intron_variant,splice_polypyrimidine_tract_variant")
        self.assertEqual(C.filter_breakend_mate_local_read_pairs(perl, rust, self.MODELS)[:2], (set(), set()))

    def test_a_vep_rs_under_annotation_is_left_charged(self) -> None:
        perl, rust = self.pair("feature_truncation", "intron_variant")
        self.assertEqual(C.filter_breakend_mate_local_read_pairs(perl, rust, self.MODELS)[:2], (set(), set()))

    def test_a_perl_set_the_local_read_does_not_explain_is_left_charged(self) -> None:
        perl, rust = self.pair("feature_truncation,coding_sequence_variant", "feature_truncation,intron_variant")
        self.assertEqual(C.filter_breakend_mate_local_read_pairs(perl, rust, self.MODELS)[:2], (set(), set()))

    def test_a_same_chromosome_mate_beyond_five_kb_of_the_local_breakend_counts_as_same_chromosome(self) -> None:
        perl = {T("21:27348903", "A[21:33034355[", self.TX, "feature_truncation")}
        rust = {T("21:27348903", "A[21:33034355[", self.TX, "feature_truncation,intron_variant")}
        ex_rust, ex_perl, kinds = C.filter_breakend_mate_local_read_pairs(perl, rust, self.MODELS)
        self.assertEqual((ex_rust, ex_perl), (rust, perl))
        self.assertEqual(kinds["same_chromosome"], 1)

    def test_a_row_the_local_coordinate_also_selects_is_not_a_mate_side_row(self) -> None:
        perl = {T("21:33034000", "N[21:33034355[", self.TX, "feature_truncation")}
        rust = {T("21:33034000", "N[21:33034355[", self.TX, "feature_truncation,intron_variant")}
        self.assertEqual(C.filter_breakend_mate_local_read_pairs(perl, rust, self.MODELS)[:2], (set(), set()))

    def test_coinciding_reads_are_not_a_pair_to_exclude(self) -> None:
        """When the mate and the local coordinate derive the same set the two engines
        agree; a divergence there is something else and stays charged."""
        gstt1 = C.TranscriptModel(chrom="22", start=24376133, end=24384284, strand=-1, biotype="protein_coding",
                                  exons=((24376133, 24376617), (24376822, 24376998), (24379361, 24379511), (24381700, 24381787), (24384120, 24384284)),
                                  cds=(24376423, 24384231), mirna_mature=())
        perl = {T("21:24372025", "N[22:24372025[", "ENST00000248935", "downstream_gene_variant")}
        rust = {T("21:24372025", "N[22:24372025[", "ENST00000248935", "intergenic_variant")}
        self.assertEqual(C.filter_breakend_mate_local_read_pairs(perl, rust, {("22", "ENST00000248935"): gstt1})[:2], (set(), set()))

    def test_one_sided_rows_are_never_pairs(self) -> None:
        perl = {T(self.LOC, self.ALLELE, self.TX, "feature_truncation")}
        self.assertEqual(C.filter_breakend_mate_local_read_pairs(perl, set(), self.MODELS)[:2], (set(), set()))

    def test_no_cache_excludes_nothing(self) -> None:
        perl, rust = self.pair("feature_truncation", "feature_truncation,intron_variant")
        self.assertEqual(C.filter_breakend_mate_local_read_pairs(perl, rust, {})[:2], (set(), set()))

    def test_a_transcript_absent_from_the_cache_is_not_adjudicated(self) -> None:
        perl, rust = self.pair("feature_truncation", "feature_truncation,intron_variant")
        self.assertEqual(C.filter_breakend_mate_local_read_pairs(perl, rust, {("22", "ENST_OTHER"): self.MODEL})[:2], (set(), set()))


# The reference release selector (Ensembl VEP 115.2 or 116.2)


def INTERGENIC(loc: str, allele: str):
    return T(loc, allele, "-", "intergenic_variant", ftype="-")


def REC(pos: int, ref: str, alt: str, info: str, vid: str = "r", chrom: str = "21") -> "C.VcfVariant":
    """An input record as parse_vcf builds it."""
    return C.VcfVariant(chrom, pos, vid, ref, alt, C.classify_variant(ref, alt), f"{chrom}:{pos}", info)


class ReferenceSkippedRecordTests(unittest.TestCase):
    """116 drops a `vep_skip` record before the input buffer (Parser.pm 116.2:491-501,
    821-827; Parser/VCF.pm 116.2:477-481, 555-558, 575) and writes no row for it; vep-rs
    annotates it. The record is named by the input (`perl_skip_route`), confirmed by the
    reference's absence under its own Location start and alleles (`record_alleles`), and its
    vep-rs rows under those alleles are the class, whatever else shares the start."""

    def test_every_vep_rs_row_of_a_skipped_record_is_excluded(self) -> None:
        records = [REC(50000, "N", "<CPX>", "END=60000;SVTYPE=CPX", "cpx"),
                   REC(30000000, "N", "<DEL>", "SVTYPE=DEL", "del_no_end"),
                   REC(70000, "N", "<NON_REF>", "END=70100", "gvcf"),
                   REC(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL", "kept")]
        rust = {
            T("21:50001-60000", "CPX", "ENST_A", "intron_variant"),
            T("21:50001-60000", "CPX", "ENST_B", "coding_sequence_variant"),
            T("21:30000001", "deletion", "ENST_C", "feature_truncation"),
            T("21:70001-70100", "<NON_REF>", "-", "intergenic_variant", ftype="-"),
            T("21:1001-2000", "deletion", "ENST_D", "feature_truncation"),
        }
        perl = {T("21:1001-2000", "deletion", "ENST_D", "feature_truncation")}
        excluded, records_n = C.filter_reference_skipped_records(perl, rust, records)
        self.assertEqual(excluded, {t for t in rust if t[0] != "21:1001-2000"})
        self.assertEqual(records_n, 3)
        self.assertEqual(C.check_skip_routes(records, perl, rust), [])

    def test_a_kept_record_at_the_same_start_neither_hides_nor_joins_a_skipped_one(self) -> None:
        """A `<NON_REF>` block beside a kept multi-span record at one start (the shape of the
        synthetic special-allele set), and a `<CPX>` sharing its whole Location string with a
        kept deletion (gnomAD-SV): the kept record's reference rows are under its own alleles,
        so the skipped record is confirmed, and only the skipped record's rows are set aside."""
        records = [REC(27330442, "N", "<NON_REF>", "END=27330631", "nonref"),
                   REC(27330442, "N", "<DEL>", "END=27330700;SVTYPE=DEL", "multi_span"),
                   REC(16450356, "N", "<CPX>", "END=16456471;SVTYPE=CPX", "cpx"),
                   REC(16450356, "N", "<DEL>", "END=16456471;SVTYPE=DEL", "del")]
        perl = {T("21:27330443-27330700", "deletion", "ENST_M", "intron_variant"),
                T("21:16450357-16456471", "deletion", "ENST_D", "feature_truncation")}
        rust = perl | {T("21:27330443-27330631", "<NON_REF>", "-", "intergenic_variant", ftype="-"),
                       T("21:16450357-16456471", "CPX", "ENST_D", "coding_sequence_variant"),
                       T("21:16450357-16456471", "CPX", "ENST_E", "intron_variant")}
        excluded, n = C.filter_reference_skipped_records(perl, rust, records)
        self.assertEqual(excluded, rust - perl)
        self.assertEqual(n, 2)
        self.assertEqual(C.check_skip_routes(records, perl, rust), [])

    def test_a_skipped_breakend_record_is_matched_by_its_own_end_and_mate_strings(self) -> None:
        """A gnomAD-SV breakend over the cap: vep-rs writes the own end `N.` and the mate
        string built from CHR2 and END2 (or END when a caller writes no END2); a kept
        breakend at the same start with another mate stays."""
        records = [REC(18457512, "N", "<BND>", "END=18457513;SVTYPE=BND;SVLEN=14027179;CHR2=21;POS2=32484691;END2=32484692", "giant"),
                   REC(16316682, "N", "<BND>", "CHR2=chr21;END=16316682;SVLEN=15787275;SVTYPE=BND", "giant_no_end2"),
                   REC(18457512, "N", "N[22:5000[", "SVTYPE=BND", "kept_bnd")]
        perl = {T("21:18457513", "N[22:5000[", "ENST_K", "intron_variant")}
        rust = perl | {T("21:18457513-32484691", "N.", "ENST_1", "intron_variant"),
                       T("21:18457513-32484691", "N[21:32484692[", "ENST_2", "downstream_gene_variant"),
                       T("21:16316683-32103957", "N.", "ENST_3", "intron_variant"),
                       T("21:16316683-32103957", "N[21:16316682[", "ENST_4", "intron_variant")}
        excluded, n = C.filter_reference_skipped_records(perl, rust, records)
        self.assertEqual(excluded, rust - perl)
        self.assertEqual(n, 2)
        self.assertEqual(C.check_skip_routes(records, perl, rust), [])

    def test_a_mixed_alt_list_is_set_aside_whole_including_its_sequence_allele_rows(self) -> None:
        """`<DEL>,T` has no Sequence Ontology term, so the reference skips the line; vep-rs
        writes a `deletion` row and a `T` row at the Location, and both are the class."""
        records = [REC(30165800, "N", "<DEL>,T", "END=30166000", "mixed")]
        rust = {T("21:30165801-30166000", "deletion", "ENST_1", "feature_truncation"),
                T("21:30165801-30166000", "T", "ENST_1", "missense_variant")}
        self.assertEqual(C.filter_reference_skipped_records(set(), rust, records), (rust, 1))

    def test_a_record_the_reference_kept_is_never_touched(self) -> None:
        """A `<DEL>,<DUP>` record is kept as `copy_number_variation`; vep-rs's `deletion` and
        `duplication` rows at its Location are a naming divergence, not a skipped record, and a
        record on no route with reference rows is likewise untouched."""
        records = [REC(30165700, "N", "<DEL>,<DUP>", "END=30165900", "cnv"),
                   REC(18457512, "N", "<BND>", "SVTYPE=BND;CHR2=21;END2=32484692", "bnd_in_cap")]
        perl = {T("21:30165701-30165900", "copy_number_variation", "ENST_1", "3_prime_UTR_variant"),
                T("21:18457513", "N.", "ENST_2", "intron_variant")}
        rust = {T("21:30165701-30165900", "deletion", "ENST_1", "feature_truncation,3_prime_UTR_variant"),
                T("21:30165701-30165900", "duplication", "ENST_1", "feature_elongation,3_prime_UTR_variant"),
                T("21:18457513", "N.", "ENST_2", "intron_variant"),
                T("21:18457513", "N.", "ENST_3", "transcript_ablation")}
        self.assertEqual(C.filter_reference_skipped_records(perl, rust, records), (set(), 0))
        self.assertEqual(C.check_skip_routes(records, perl, rust), [])

    def test_a_tandem_repeat_record_is_outside_the_class(self) -> None:
        """The reference keeps a <CNV:TR> and writes it at the run's 3' end, so its vep-rs
        rows never share a Location with the reference's and belong to the <CNV:TR> class."""
        records = [REC(899, "N", "<CNV:TR>", "END=950;SVLEN=51;RUS=AC;RUC=10", "tr")]
        rust = {T("21:900-950", "tandem_repeat", "ENST_T", "feature_elongation")}
        perl = {T("21:950-951", "ACG", "ENST_T", "inframe_insertion")}
        self.assertEqual(C.filter_reference_skipped_records(perl, rust, records), (set(), 0))
        self.assertEqual(C.check_skip_routes(records, perl, rust), [])

    def test_the_mask_without_input_records_is_refused(self) -> None:
        with self.assertRaisesRegex(ValueError, "input_variants"):
            C.apply_masks(set(), set(), "116.2", {}, {}, {}, "vep-rs")


class RecordAllelesTests(unittest.TestCase):
    """`record_alleles` names every Allele string the engines write for a record."""

    def test_symbolic_types_with_and_without_a_term(self) -> None:
        self.assertEqual(C.record_alleles(REC(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL")), {"deletion"})
        self.assertEqual(C.record_alleles(REC(1000, "N", "<CPX>", "END=2000;SVTYPE=CPX")), {"CPX", "<CPX>"})
        self.assertEqual(C.record_alleles(REC(1000, "N", "<NON_REF>", "END=1100")), {"NON_REF", "<NON_REF>"})
        self.assertEqual(C.record_alleles(REC(1000, "N", "<DEL>,<DUP>", "END=2000")), {"deletion", "duplication", "copy_number_variation"})
        self.assertEqual(C.record_alleles(REC(1000, "N", "<INS>,<DEL>", "END=2000")), {"insertion", "deletion"})
        self.assertEqual(C.record_alleles(REC(1000, "N", "<DEL>,T", "END=2000")), {"deletion", "T"})
        self.assertEqual(C.record_alleles(REC(1000, "N", "<INS:ME:ALU>", "SVLEN=300")), {"Alu_insertion"})

    def test_breakend_records(self) -> None:
        self.assertEqual(
            C.record_alleles(REC(18457512, "N", "<BND>", "SVTYPE=BND;CHR2=21;POS2=32484691;END2=32484692;END=18457513")),
            {"chromosome_breakpoint", "N.", ".N", "N[21:32484692[", "N[21:32484691[", "N[21:18457513["},
        )
        self.assertEqual(
            C.record_alleles(REC(100, "A", "A[chr22:5000[", "SVTYPE=BND")),
            {"A[22:5000[", "N.", ".N", "A.", ".A", "chromosome_breakpoint"},
        )
        self.assertEqual(C.record_alleles(REC(100, "T", ".T", "SVTYPE=BND")), {".T", "chromosome_breakpoint", "N.", ".N", "T."})


class InfoValueTests(unittest.TestCase):
    def test_an_info_value_is_read_to_a_second_equals_sign(self) -> None:
        """`my ($key,$value) = split('=',$info)` (BaseVCF4.pm:582) is a two-element list
        assignment, so a value ends at a second `=`. `21 30425650 x N G . . SVTYPE==DEL;END=30425949`
        has the empty SVTYPE, false, and takes the sequence path: both engines write
        `21:30425650 G ENST00000390690 coding_sequence_variant` and `21:30425650 G ENST00000334068
        upstream_gene_variant`, the rows of `21 30425650 x A G . . SVTYPE=DEL;END=30425949`
        (`missense_variant` on ENST00000390690) beside them, and neither record is on a skip route.
        The mate both engines write for `21 30425650 x N <BND> . . SVTYPE=BND;CHR2=21;END2=30425700=5`
        is `N[21:30425700[` (`21:30425651 N.` and `21:30425651 N[21:30425700[` on both transcripts),
        the string `record_alleles` names."""
        self.assertEqual(C._info_fields("SVTYPE==DEL;END=30425949"), {"SVTYPE": "", "END": "30425949"})
        self.assertEqual(C._info_fields("SVTYPE=BND;CHR2=21;END2=30425700=5")["END2"], "30425700")
        self.assertEqual(C._info_fields("END2=30425700=5=7;CHR2"), {"END2": "30425700", "CHR2": ""})
        eq_eq = REC(30425650, "N", "G", "SVTYPE==DEL;END=30425949", "eq_eq")
        seq = REC(30425650, "A", "G", "SVTYPE=DEL;END=30425949", "seq")
        self.assertTrue(C.perl_sequence_path(eq_eq))
        self.assertIsNone(C.perl_skip_route(eq_eq))
        self.assertEqual(C.vep_rs_record_rows(eq_eq), {("21:30425650", "G")})
        rows = {T("21:30425650", "G", "ENST00000390690", "coding_sequence_variant"),
                T("21:30425650", "G", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425650", "G", "ENST00000390690", "missense_variant")}
        self.assertEqual(C.check_skip_routes([eq_eq, seq], rows, rows), [])
        self.assertEqual(C.filter_reference_skipped_records(rows, rows, [eq_eq, seq]), (set(), 0))
        bnd = REC(30425650, "N", "<BND>", "SVTYPE=BND;CHR2=21;END2=30425700=5", "bnd")
        self.assertEqual(C.record_alleles(bnd), {"chromosome_breakpoint", "N.", ".N", "N[21:30425700["})
        mate_rows = {T("21:30425651", allele, feature, csq)
                     for allele in ("N.", "N[21:30425700[")
                     for feature, csq in (("ENST00000334068", "upstream_gene_variant"),
                                          ("ENST00000390690", "feature_truncation,coding_sequence_variant"))}
        self.assertEqual(C.check_skip_routes([bnd], mate_rows, mate_rows), [])
        self.assertEqual(ADJUSTED(mate_rows, mate_rows, [bnd])[:2], (set(), set()))


class RecordLocationTests(unittest.TestCase):
    def test_start_is_the_base_after_pos_and_end_follows_get_end(self) -> None:
        self.assertEqual(C.record_location(REC(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL")), ("21", 1001, 2000))
        self.assertEqual(C.record_location(REC(1000, "N", "<DEL>", "END=2000;SVLEN=-500;SVTYPE=DEL")), ("21", 1001, 1500), "SVLEN before END")
        self.assertEqual(C.record_location(REC(30165850, "N", "<DEL>", "SVTYPE=DEL")), ("21", 30165851, 30165851))
        self.assertEqual(C.record_location(REC(18457512, "N", "<BND>", "END=18457513;SVLEN=14027179;CHR2=21;END2=32484692")), ("21", 18457513, 32484691))
        self.assertEqual(C.record_location(REC(18457512, "N", "<BND>", "END=18457513;CHR2=21;END2=32484692")), ("21", 18457513, 18457513), "a breakend's END is discarded")

    def test_an_infinite_svlen_or_end_is_a_saturated_integer_end(self) -> None:
        """Perl reads `1e400` as `Inf` (`abs("1e400")` is `Inf`) and Ensembl VEP 116.2 skips the
        record as oversize (`variant size (Inf) is bigger than --max_sv_size`); the end derived
        for the record is saturated at the machine word and stays an integer, so the record's
        Location is computed, never an overflow that ends the comparison."""
        saturated = 2**63 - 1
        self.assertEqual(C.record_location(REC(1000, "N", "<DUP>", "SVTYPE=DUP;SVLEN=1e400")), ("21", 1001, 1000 + saturated))
        self.assertEqual(C.record_location(REC(1000, "N", "<DUP>", "SVTYPE=DUP;END=1e400")), ("21", 1001, saturated))
        self.assertEqual(C.record_location(REC(1000, "N", "<DUP>", "SVTYPE=DUP;END=-1e400")), ("21", 1001, -saturated))
        for chrom, _start, end in (C.record_location(REC(1000, "N", "<DUP>", info)) for info in ("SVTYPE=DUP;SVLEN=1e400", "SVTYPE=DUP;END=1e400")):
            self.assertIsInstance(end, int, chrom)
        oversize = REC(1000, "N", "<DUP>", "SVTYPE=DUP;SVLEN=1e400")
        rust = {T("21:1001-9223372036854775807", "duplication", "ENST_1", "intron_variant")}
        self.assertEqual(C.filter_reference_skipped_records(set(), rust, [oversize]), (rust, 1))
        self.assertEqual(C.check_skip_routes([oversize], set(), rust), [])


class BndOwnRowsLostTests(unittest.TestCase):
    """A bracket breakend whose mate is on another chromosome loses its own-end rows
    when a mate-chromosome transcript, loaded for an earlier record of the same buffer, names
    its slice first (AnnotationType/Transcript.pm 116.2:122-143 with InputBuffer.pm
    116.2:191-207, 223-234 and AnnotationSource.pm 116.2:235-258)."""

    LOC = "21:26037668"
    MATE = "N[22:26037668["
    TX_BY_CHR = {"21": {"ENST21_NEAR", "ENST21_FAR"}, "22": {"ENST22_M"}}
    SPANS = {
        ("21", "ENST21_NEAR"): ("21", 26030000, 26035000),   # 2,668 bases from the breakpoint
        ("21", "ENST21_FAR"): ("21", 26000000, 26010000),    # 27,668 bases away
        ("22", "ENST22_M"): ("22", 26036000, 26040000),
    }

    def rust(self, *own_features: str) -> set:
        return {T(self.LOC, self.MATE, "ENST22_M", "intron_variant")} | {
            T(self.LOC, "N.", f, "upstream_gene_variant") for f in own_features
        }

    def test_own_end_rows_within_the_admission_window_are_excluded(self) -> None:
        perl = {T(self.LOC, self.MATE, "ENST22_M", "intron_variant")}
        got = C.filter_bnd_own_rows_lost(perl, self.rust("ENST21_NEAR"), self.TX_BY_CHR, self.SPANS)
        self.assertEqual(got, {T(self.LOC, "N.", "ENST21_NEAR", "upstream_gene_variant")})

    def test_vep_rs_own_allele_string_may_be_the_class_term(self) -> None:
        """vep-rs writes the local end of a bracket breakend as `chromosome_breakpoint`."""
        perl = {T(self.LOC, self.MATE, "ENST22_M", "intron_variant")}
        rust = {T(self.LOC, self.MATE, "ENST22_M", "intron_variant"),
                T(self.LOC, "chromosome_breakpoint", "ENST21_NEAR", "upstream_gene_variant")}
        self.assertEqual(len(C.filter_bnd_own_rows_lost(perl, rust, self.TX_BY_CHR, self.SPANS)), 1)

    def test_a_transcript_beyond_the_window_stays_charged(self) -> None:
        """The reference's own rule admits a breakend within 5,000 bases of the transcript
        (StructuralVariationOverlap.pm 116:144-160); a farther row is vep-rs's alone."""
        perl = {T(self.LOC, self.MATE, "ENST22_M", "intron_variant")}
        got = C.filter_bnd_own_rows_lost(perl, self.rust("ENST21_NEAR", "ENST21_FAR"), self.TX_BY_CHR, self.SPANS)
        self.assertEqual({t[2] for t in got}, {"ENST21_NEAR"})

    def test_the_reference_keeping_any_own_end_row_is_not_the_shape(self) -> None:
        perl = {T(self.LOC, self.MATE, "ENST22_M", "intron_variant"), INTERGENIC(self.LOC, "N.")}
        self.assertEqual(C.filter_bnd_own_rows_lost(perl, self.rust("ENST21_NEAR"), self.TX_BY_CHR, self.SPANS), set())

    def test_a_same_chromosome_mate_is_not_the_shape(self) -> None:
        loc, mate = "21:24808769-35740566", "N[21:35740567["
        perl = {T(loc, mate, "ENST21_NEAR", "intron_variant")}
        rust = perl | {T(loc, "N.", "ENST21_NEAR", "intron_variant")}
        self.assertEqual(C.filter_bnd_own_rows_lost(perl, rust, self.TX_BY_CHR, self.SPANS), set())

    def test_the_reference_dropping_the_mate_rows_too_is_not_the_shape(self) -> None:
        """No reference Transcript row under the bracket allele means the record was not
        reached through its mate; that is the skipped-record class or a real divergence."""
        perl = {INTERGENIC(self.LOC, self.MATE)}
        self.assertEqual(C.filter_bnd_own_rows_lost(perl, self.rust("ENST21_NEAR"), self.TX_BY_CHR, self.SPANS), set())

    def test_no_cache_excludes_nothing(self) -> None:
        perl = {T(self.LOC, self.MATE, "ENST22_M", "intron_variant")}
        self.assertEqual(C.filter_bnd_own_rows_lost(perl, self.rust("ENST21_NEAR"), {}, {}), set())


class BndSynonymMateRowsTests(unittest.TestCase):
    """A mate whose chromosome carries a prefix the mate-allele test does not resolve (`Chr21`
    for `21`) loses its mate rows in the reference alone (StructuralVariationOverlap.pm
    116:86-102, 144-160; Utils/VariationEffect.pm 116:111-124); a mate named by an accession
    on another chromosome passes the overlap test (InputBuffer.pm 116.2:348-357) and fails the
    slice test (AnnotationType/Transcript.pm 116.2:134), while vep-rs resolves it through the
    synonym table and writes the mate rows, which are set aside."""

    LOC = "21:26037668"
    SYN = "N[NC_000022.11:26037668["
    TX_BY_CHR = {"21": {"ENST21_OWN"}, "22": {"ENST22_M"}}

    def test_mate_side_rows_under_a_synonym_named_mate_are_excluded(self) -> None:
        perl = {T(self.LOC, "N.", "ENST21_OWN", "upstream_gene_variant")}
        rust = perl | {T(self.LOC, self.SYN, "ENST22_M", "intron_variant")}
        self.assertEqual(C.filter_bnd_synonym_mate_rows(perl, rust, self.TX_BY_CHR), {T(self.LOC, self.SYN, "ENST22_M", "intron_variant")})

    def test_a_cache_named_mate_is_not_the_shape(self) -> None:
        perl = {T(self.LOC, "N.", "ENST21_OWN", "upstream_gene_variant")}
        rust = perl | {T(self.LOC, "N[22:26037668[", "ENST22_M", "intron_variant")}
        self.assertEqual(C.filter_bnd_synonym_mate_rows(perl, rust, self.TX_BY_CHR), set())

    def test_an_own_chromosome_transcript_under_the_mate_allele_stays_charged(self) -> None:
        perl = {T(self.LOC, "N.", "ENST21_OWN", "upstream_gene_variant")}
        rust = perl | {T(self.LOC, self.SYN, "ENST21_OWN", "intron_variant")}
        self.assertEqual(C.filter_bnd_synonym_mate_rows(perl, rust, self.TX_BY_CHR), set())

    def test_the_reference_writing_the_mate_allele_is_not_the_shape(self) -> None:
        perl = {T(self.LOC, self.SYN, "ENST22_M", "intron_variant")}
        rust = perl | {T(self.LOC, self.SYN, "ENST22_X", "downstream_gene_variant")}
        self.assertEqual(C.filter_bnd_synonym_mate_rows(perl, rust, self.TX_BY_CHR), set())

    def test_a_record_the_reference_wrote_nothing_for_is_the_skipped_class_not_this_one(self) -> None:
        rust = {T(self.LOC, self.SYN, "ENST22_M", "intron_variant")}
        self.assertEqual(C.filter_bnd_synonym_mate_rows(set(), rust, self.TX_BY_CHR), set())

    def test_no_cache_excludes_nothing(self) -> None:
        perl = {T(self.LOC, "N.", "ENST21_OWN", "upstream_gene_variant")}
        rust = perl | {T(self.LOC, self.SYN, "ENST22_M", "intron_variant")}
        self.assertEqual(C.filter_bnd_synonym_mate_rows(perl, rust, {}), set())


class FastaNamedSliceIntergenicTests(unittest.TestCase):
    def test_counts_records_whose_only_reference_row_is_intergenic_against_vep_rs_transcripts(self) -> None:
        perl = {INTERGENIC("21:100", "A"), INTERGENIC("21:200", "T"), T("21:300", "G", "ENST_3", "intron_variant")}
        rust = {T("21:100", "A", "ENST_1", "missense_variant"), T("21:100", "A", "ENST_1b", "intron_variant"),
                INTERGENIC("21:200", "T"), T("21:300", "G", "ENST_3", "intron_variant")}
        self.assertEqual(C.count_fasta_named_slice_intergenic(perl, rust), 1)

    def test_a_mate_past_the_buffers_32_bit_bound_is_counted_from_the_input_record(self) -> None:
        """Ensembl VEP 116.2 turns the buffer intergenic for a mate at 2147478648 named by
        a bracket ALT or by CHR2 with END2 (else END); 2147478647 and a mate on another
        route are not the shape. The reference's tuples at the record's Location are
        counted with it."""
        V = C.VcfVariant
        records = [
            V("21", 30425650, "bracket", "N", "N[21:2147478648[", "BND", "21:30425650", "SVTYPE=BND"),
            V("21", 30425650, "chr2_end2", "N", "<BND>", "BND", "21:30425650", "SVTYPE=BND;CHR2=21;END2=2147478648"),
            V("21", 30425650, "chr2_end", "N", "<BND>", "BND", "21:30425650", "SVTYPE=BND;CHR2=21;END=2147478648"),
            V("21", 30425650, "below", "N", "<BND>", "BND", "21:30425650", "SVTYPE=BND;CHR2=21;END2=2147478647"),
            V("21", 30425650, "pos2", "N", "<BND>", "BND", "21:30425650", "SVTYPE=BND;CHR2=21;POS2=2147478648"),
        ]
        perl = {INTERGENIC("21:30425651", "N."), INTERGENIC("21:30425651", "N[21:2147478648["),
                T("21:40000001", "N.", "ENST_1", "intron_variant")}
        self.assertEqual(C.count_bnd_mate_over_32_bit(records, perl), (3, 2))
        self.assertEqual(C.count_bnd_mate_over_32_bit(records[3:], perl), (0, 0))

    def test_a_gvcf_variant_site_is_its_own_class_apart_from_the_skipped_records(self) -> None:
        """Ensembl VEP 116.2 drops `A G,<NON_REF>` as `G/NON_REF is not a supported structural
        variant type`; vep-rs annotates the `G` allele at POS. The site's vep-rs tuples are set
        aside under `gvcf_variant_site_dropped`, and the skipped-record mask, which names
        every other route, leaves them alone."""
        records = [REC(30425650, "A", "G,<NON_REF>", "."), REC(50000, "N", "<CPX>", "END=60000;SVTYPE=CPX")]
        perl: set = set()
        rust = {T("21:30425650", "G", "ENST_1", "missense_variant"), T("21:30425650", "G", "ENST_2", "upstream_gene_variant"),
                T("21:50001-60000", "CPX", "ENST_3", "intron_variant")}
        excluded, sites = C.filter_gvcf_variant_sites(perl, rust, records)
        self.assertEqual((sites, excluded), (1, {t for t in rust if t[1] == "G"}))
        skipped, records_n = C.filter_reference_skipped_records(perl, rust, records)
        self.assertEqual((records_n, skipped), (1, {T("21:50001-60000", "CPX", "ENST_3", "intron_variant")}))
        self.assertEqual(C.check_skip_routes(records, perl, rust), [])
        masks = C.apply_masks(perl, rust, "116.2", {}, {}, {}, "vep-rs", records)
        self.assertEqual((masks.skipped_records, masks.gvcf_records), (1, 1), "the mask run names the site once")

    def test_a_dot_breakend_pair_is_set_aside_on_both_sides(self) -> None:
        """`N N.` without SVTYPE: the reference's `21:30425650-30425651 chromosome_breakpoint`
        row and vep-rs's `21:30425651` row on one feature are a pair; a reference row on a
        feature vep-rs did not name, or a point vep-rs wrote under a bracket allele, is not."""
        perl = {T("21:30425650-30425651", "chromosome_breakpoint", "ENST_1", "coding_sequence_variant"),
                T("21:30425650-30425651", "chromosome_breakpoint", "ENST_2", "upstream_gene_variant")}
        rust = {T("21:30425651", "chromosome_breakpoint", "ENST_1", "coding_sequence_variant"),
                T("21:30425651", "N[21:30425700[", "ENST_2", "upstream_gene_variant")}
        excluded_rust, excluded_perl = C.filter_dot_breakend_pairs(perl, rust)
        self.assertEqual(excluded_rust, {T("21:30425651", "chromosome_breakpoint", "ENST_1", "coding_sequence_variant")})
        self.assertEqual(excluded_perl, {T("21:30425650-30425651", "chromosome_breakpoint", "ENST_1", "coding_sequence_variant")})

    def test_a_missing_alt_the_reference_annotated_is_set_aside_from_the_input(self) -> None:
        """ALT `.` with `SVTYPE=BND` is a `chromosome_breakpoint` at POS + 1 and with
        `SVTYPE=DEL;END` a deletion over END in Ensembl VEP 116.2; vep-rs writes nothing. The
        reference's tuples at the record's Location are set aside; a `.` without SVTYPE
        produced no reference row and names nothing."""
        records = [REC(30425650, "N", ".", "SVTYPE=BND"), REC(40000000, "N", ".", "SVTYPE=DEL;END=40000100"),
                   REC(50000000, "N", ".", ".")]
        perl = {T("21:30425651", "chromosome_breakpoint", "ENST_1", "coding_sequence_variant"),
                T("21:40000001-40000100", "deletion", "ENST_2", "intron_variant"),
                T("21:50000001", "deletion", "ENST_3", "intron_variant")}
        self.assertEqual(C.filter_missing_alt_annotated(perl, records), {
            T("21:30425651", "chromosome_breakpoint", "ENST_1", "coding_sequence_variant"),
            T("21:40000001-40000100", "deletion", "ENST_2", "intron_variant")})

    def test_a_non_integer_end_row_is_set_aside_and_a_text_end_is_read_as_its_number(self) -> None:
        """`END=30425949.5` and `END=3.0425949e7` are annotated and printed as written by
        Ensembl VEP 116.2 (`21:30425651-30425949.5`), where vep-rs skips the record;
        `END=030425949` and `END=+30425949` name the same span as `30425949` and are
        normalised to it, so no row is set aside for them."""
        perl = {T("21:30425651-30425949.5", "deletion", "ENST_1", "intron_variant"),
                T("21:30425651-3.0425949e7", "deletion", "ENST_2", "intron_variant"),
                T("21:30425651-30425949", "deletion", "ENST_3", "intron_variant")}
        self.assertEqual(C.filter_non_integer_end_rows(perl, set(), []), ({t for t in perl if t[2] != "ENST_3"}, set()))
        for text in ("21:30425651-030425949", "21:30425651-+30425949", "chr21:30425651-30425949"):
            self.assertEqual(C.normalize_location(text), "21:30425651-30425949", text)
        self.assertEqual(C.normalize_location("21:30425651-30425949.5"), "21:30425651-30425949.5")

    def test_an_empty_tandem_repeat_expansion_is_set_aside_only_when_the_reference_wrote_nothing(self) -> None:
        """With a FASTA, Ensembl VEP 116.2 writes nothing for a `<CNV:TR>` whose expansion is
        empty or unchanged (`RN=.`, `RUC=abc`, `RUC=-2`, `RUC=3` on a three-unit run) and no
        warning, and vep-rs with a FASTA writes nothing for it either; without one both
        engines write the symbolic rows over the run. vep-rs's symbolic rows are set aside
        only when the reference wrote nothing at or beside the run."""
        run = "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;"
        empties = [REC(30425587, "N", "<CNV:TR>", run + "RUC=3"), REC(30425587, "N", "<CNV:TR>", run + "RN=.;RUC=4"),
                   REC(30425587, "N", "<CNV:TR>", run + "RUC=abc"), REC(30425587, "N", "<CNV:TR>", run + "RUC=-2")]
        gain = REC(30425587, "N", "<CNV:TR>", run + "RUC=4")
        for v in empties:
            self.assertTrue(C.cnv_tr_expansion_is_empty(v), v.info)
        self.assertFalse(C.cnv_tr_expansion_is_empty(gain))
        self.assertFalse(C.cnv_tr_expansion_is_empty(REC(30425587, "N", "<CNV:TR>,<CNV:TR>", run + "RUC=3,3")))
        # `$info->{RUC} || $info->{RB}`: a RUC Perl reads as false falls through to RB.
        self.assertFalse(C.cnv_tr_expansion_is_empty(REC(30425587, "N", "<CNV:TR>", run + "RUC=0;RB=12")))
        self.assertTrue(C.cnv_tr_expansion_is_empty(REC(30425587, "N", "<CNV:TR>", run + "RUC=0;RB=9")))
        rust = {T("21:30425588-30425596", "tandem_repeat", "ENST_1", "feature_elongation,coding_sequence_variant")}
        self.assertEqual(C.filter_cnv_tr_empty_expansion(set(), rust, empties[:1]), rust)
        symbolic_reference = {T("21:30425588-30425596", "tandem_repeat", "ENST_1", "feature_elongation,coding_sequence_variant")}
        self.assertEqual(C.filter_cnv_tr_empty_expansion(symbolic_reference, rust, empties[:1]), set())
        literal_reference = {T("21:30425596-30425597", "AGC", "ENST_1", "inframe_insertion")}
        self.assertEqual(C.filter_cnv_tr_empty_expansion(literal_reference, rust, [gain]), set())

    def test_an_expansion_count_is_truncated_as_perls_repetition_truncates_it(self) -> None:
        """`$seq x $num` (Parser/VCF.pm:419) repeats the unit an integer number of times, so
        `RUC=3.5` on a three-unit `AGC` run rebuilds the run (empty) and `RUC=0.5` repeats
        nothing (empty), while `RUC=4.5` is a gain; a `.` unit is one `N` base; with RN the
        bases of every repeat sequence are summed and a two-sequence allele of eleven bases
        is not the nine-base run, and a unit of no bases (`RUS=,AC` with RB) divides by zero
        in Perl, so the record is not judged."""
        run = "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;"
        for info, empty in (
            (run + "RUC=3.5", True),
            (run + "RUC=3.9", True),
            (run + "RUC=0.5", True),
            (run + "RUC=4.5", False),
            (run + "RUC=2.5", False),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RUS=.;RUC=9", True),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RUS=.;RUC=10", False),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RN=2;RUS=AGC,T;RUC=3,2", False),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RN=2;RUS=AGC,T;RUC=0,0.9", True),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RN=2;RUS=AGC;RUC=3,2", False),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RUS=,AC;RB=12", False),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RUS=,AC;RUC=0;RB=12", False),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=0;RB=10.5", True),
            ("SVTYPE=CNV;END=30425596;SVLEN=9,12;RUS=AGC;RUC=4", True),
            # `$info->{RN} ?` and `$info->{RUS} &&` (`:385`, `:394`): a `0` or empty RN is one
            # sequence, a `0` RUS no expansion.
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RN=0;RUS=AGC;RUC=4", False),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RN=;RUS=AGC;RUC=3", True),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RUS=0;RUC=3", False),
            ("SVTYPE=CNV;END=30425596;SVLEN=9;RN=-1;RUS=AGC,T;RUC=3,2", True),
        ):
            with self.subTest(info=info):
                self.assertEqual(C.cnv_tr_expansion_is_empty(REC(30425587, "N", "<CNV:TR>", info)), empty)

    def test_a_repeat_unit_of_no_bases_is_empty_under_ruc_and_unjudged_under_rb(self) -> None:
        """`$repeat .= $seq x $num` (Parser/VCF.pm:419) with the empty unit of `RUS=,AGC` and
        `RUC=5` is `"" x 5`, nothing, so with a FASTA Ensembl VEP 116.2 builds the allele string
        `N/` and writes nothing for the record: an empty expansion. `RUS=,AC;RB=12` reaches
        `shift(@RB) / length($seq)` (`:414`) with a length of zero and the reference run ends
        (`Illegal division by zero`), so no reference output can hold the record and it is not
        judged. With `RN=2;RUC=5,2` the second unit repeats `AGC` twice, so the allele is not
        empty."""
        run = "SVTYPE=CNV;END=30425596;SVLEN=9;"
        for info, empty in (
            (run + "RUS=,AGC;RUC=5", True),
            (run + "RUS=,AC;RB=12", False),
            (run + "RUS=,AGC;RN=2;RUC=5,2", False),
            (run + "RUS=,AGC;RN=2;RUC=5,0", True),
            (run + "RUS=,AC;RUC=0;RB=12", False),
        ):
            with self.subTest(info=info):
                self.assertEqual(C.cnv_tr_expansion_is_empty(REC(30425587, "N", "<CNV:TR>", info)), empty)

    def test_a_gvcf_site_allele_is_upper_cased_as_the_parser_writes_it(self) -> None:
        """Ensembl VEP 116.2 drops `21 30425650 x A g,<NON_REF> . . .` and `21 30425650 x a
        g,<non_ref> . . .` (`g/NON_REF is not a supported structural variant type`); vep-rs
        writes the allele in upper case, `21:30425650 G ENST00000390690 missense_variant` and
        `21:30425650 G ENST00000334068 upstream_gene_variant`, and the site's rows are set aside
        under those spellings."""
        rust = {T("21:30425650", "G", "ENST00000390690", "missense_variant"),
                T("21:30425650", "G", "ENST00000334068", "upstream_gene_variant")}
        for ref, alt in (("A", "g,<NON_REF>"), ("a", "g,<non_ref>")):
            with self.subTest(ref=ref, alt=alt):
                record = REC(30425650, ref, alt, ".")
                self.assertEqual(C.perl_skip_route(record), "gvcf_variant_site")
                self.assertEqual(C.filter_gvcf_variant_sites(set(), rust, [record]), (rust, 1))
                self.assertEqual(ADJUSTED(set(), rust, [record])[:2], (set(), set()))

    def test_a_gvcf_site_is_set_aside_at_the_alleles_the_parser_writes(self) -> None:
        """vep-rs places a gVCF site's alleles as its parser does (`gvcf_site_rows`): `ATG AG`
        is the deletion of `T` at POS + 1 (anchor trimmed, then minimised from both ends), an
        insertion list keeps its anchor-trimmed alleles at the pair, and a block list is a
        skipped record, not a site."""
        records = [REC(30425650, "ATG", "AG,<NON_REF>", "."), REC(40000000, "A", "AT,ATT,<NON_REF>", "."),
                   REC(50000000, "N", "<NON_REF>,<NON_REF>", "END=50000100")]
        rust = {T("21:30425651", "-", "ENST_1", "frameshift_variant"), T("21:30425651", "G", "ENST_1", "missense_variant"),
                T("21:40000000-40000001", "T", "ENST_2", "frameshift_variant"),
                T("21:40000000-40000001", "TT", "ENST_2", "inframe_insertion"),
                T("21:50000001-50000100", "<NON_REF>", "-", "intergenic_variant", "-")}
        excluded, sites = C.filter_gvcf_variant_sites(set(), rust, records)
        self.assertEqual(sites, 2)
        self.assertEqual(excluded, {T("21:30425651", "-", "ENST_1", "frameshift_variant"),
                                    T("21:40000000-40000001", "T", "ENST_2", "frameshift_variant"),
                                    T("21:40000000-40000001", "TT", "ENST_2", "inframe_insertion")})
        self.assertEqual(C.perl_skip_route(records[2]), "unsupported_type")
        self.assertEqual(C.perl_skip_route(REC(1000, "N", "*,<NON_REF>", ".")), "unsupported_type")
        self.assertEqual(C.perl_skip_route(REC(1000, "N", "<DEL>,<NON_REF>", "END=1100")), "unsupported_type")
        # Under a true SVTYPE the block is a deletion allele in both engines: the record is kept.
        self.assertIsNone(C.perl_skip_route(REC(1000, "N", "G,<NON_REF>", "SVTYPE=DEL;END=1100")))
        self.assertEqual(C.perl_skip_route(REC(1000, "N", "G,<NON_REF>", "SVTYPE=FOO;END=1100")), "unsupported_type")
        self.assertEqual(C.perl_skip_route(REC(1000, "N", "G,<NON_REF>", "SVTYPE=;END=1100")), "gvcf_variant_site")
        skipped, n = C.filter_reference_skipped_records(set(), rust, records)
        self.assertEqual((n, skipped), (1, {T("21:50000001-50000100", "<NON_REF>", "-", "intergenic_variant", "-")}))

    def test_perl_numification_and_the_structural_end(self) -> None:
        for text, value in (("300", 300), ("-300", -300), ("300abc", 300), (" 3e2", 300), ("abc", 0), ("", 0), (".", 0), ("300.7", 300.7)):
            self.assertEqual(C.perl_numeric_value(text), value, text)
        # `abs("1e400")` is `Inf` to Perl; the reading saturates at the machine word with its
        # sign, as vep-rs's does, so every later `int()` of a coordinate holds.
        self.assertEqual(C.perl_numeric_value("1e400"), 2**63 - 1)
        self.assertEqual(C.perl_numeric_value("-1e400"), -(2**63 - 1))
        self.assertEqual(C.perl_numeric_value("1.7e308"), 1.7e308, "the largest finite reading is not saturated")
        self.assertEqual(C.perl_structural_end(1000, "1e400", None, False), 1000 + 2**63 - 1)
        self.assertEqual(C.perl_structural_end(1000, None, "1e400", False), 2**63 - 1)
        self.assertEqual(C.perl_structural_end(1000, "-300", "2000", False), 1300)
        self.assertEqual(C.perl_structural_end(1000, "abc", "2000", False), 1000)
        self.assertEqual(C.perl_structural_end(1000, None, "2000", False), 2000)
        self.assertEqual(C.perl_structural_end(1000, None, "2000.5", False), 2000.5)
        self.assertIsNone(C.perl_structural_end(1000, None, "abc", False))
        self.assertIsNone(C.perl_structural_end(1000, None, "", False))
        self.assertEqual(C.perl_structural_end(1000, None, "2000", True), 1001)
        self.assertEqual(C.perl_structural_end(1000, None, "0", True), 0)
        self.assertIsNone(C.perl_structural_end(1000, None, "", True))
        self.assertEqual(C.perl_structural_end(1000, "100", "0", True), 1100)

    def test_a_vep_rs_intergenic_row_beside_its_transcripts_is_not_the_shape(self) -> None:
        perl = {INTERGENIC("21:100", "A")}
        rust = {INTERGENIC("21:100", "A"), T("21:100", "A", "ENST_1", "missense_variant")}
        self.assertEqual(C.count_fasta_named_slice_intergenic(perl, rust), 0)


class ReferenceReleaseFixtureTests(unittest.TestCase):
    """The comparator on the committed fixture under each release.

    `testdata/reference_release/perl/` is a 115.2-shaped reference (it carries batch rows
    for a `<CPX>` and a `<NON_REF>` record); `expected_115_2/` is the report the 115.2
    registry writes on it, and every output file must equal it byte for byte, so the
    published figures reproduce. `perl_116/` is the same reference as 116.2 writes it, the
    two skipped records dropped: under 116.2 the class-C arms and the class-D masks are not
    applied and the records the reference skipped are set aside whole.
    """

    FIXTURE = Path(__file__).resolve().parent / "testdata" / "reference_release"
    FILES = ("concordance_report.json", "concordance_report.md", "discordant.tsv", "discordant_adjusted.tsv")

    def _run(self, release: str, perl: str = "perl") -> tuple[Path, int, str]:
        import io
        import shutil
        import tempfile
        from contextlib import redirect_stderr

        out = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, out, True)
        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py", "--reference-release", release, "--assembly", "grch37",
            "--input-dir", str(self.FIXTURE / "in"), "--perl-dir", str(self.FIXTURE / perl),
            "--rust-dir", str(self.FIXTURE / "rust"), "--vep-rs-cache", str(self.FIXTURE / "cache"),
            "--output-dir", str(out),
        ]
        rc, err = 0, io.StringIO()
        try:
            with redirect_stderr(err):
                C.main()
        except SystemExit as exc:
            rc = exc.code or 0
        finally:
            sys.argv = argv
        return out, rc, err.getvalue()

    def test_115_2_reproduces_the_published_report_byte_for_byte(self) -> None:
        out, rc, _ = self._run("115.2")
        self.assertEqual(rc, 0)
        for name in self.FILES:
            with self.subTest(file=name):
                self.assertEqual((out / name).read_bytes(), (self.FIXTURE / "expected_115_2" / name).read_bytes())

    def test_116_2_sets_aside_the_skipped_records_whole_and_leaves_the_retired_classes_charged(self) -> None:
        import json

        out, rc, _ = self._run("116.2", perl="perl_116")
        self.assertEqual(rc, 0)
        adj = json.loads((out / "concordance_report.json").read_text())["adjusted"]
        self.assertEqual(adj["reference_release"], "116.2")
        self.assertEqual(adj["masks_applied"], list(C.MASKS_BY_RELEASE["116.2"]))
        self.assertEqual(set(adj["masks_not_applied"]), set(C.MASKS_RETIRED_IN_116_2))
        # The <CPX> the reference dropped (2 rows), the gVCF block (1), the <CPX> with no
        # reference row (2) and the incomplete deletion (1): four records, six rows.
        self.assertEqual((adj["excluded_reference_skipped_rust"], adj["reference_skipped_records"]), (6, 4))
        self.assertEqual(adj["excluded_reference_skipped_by_file"], {"sv": {"rust": 6, "records": 4}})
        for key in ("excluded_cross_chromosome_perl", "excluded_cross_chromosome_rust",
                    "excluded_transcript_selection_rust", "excluded_non_ref_batch_perl"):
            self.assertEqual(adj[key], 0, key)
        self.assertTrue(adj["breakend_masks_check_ran"])
        self.assertEqual(set(adj["divergence_classes"]), {
            "reference_skipped_record", "bnd_own_chromosome_rows_lost",
            "bnd_synonym_mate_rows_lost", "fasta_named_slice_intergenic",
            "bnd_bare_reference_allele_row", "bnd_mate_32bit_buffer_intergenic",
            "gvcf_variant_site_dropped", "dot_breakend_without_svtype_pair", "missing_alt_annotated",
            "non_integer_end_annotated", "cnv_tr_empty_expansion_dropped", "star_alone_svtype_annotated",
            "non_numeric_pos_numified", "cnv_tr_rn_zero_allele_dropped",
        })
        for key in ("excluded_non_integer_end_rust", "excluded_star_alone_svtype_perl", "excluded_non_numeric_pos_perl",
                    "excluded_cnv_tr_rn_zero_rust"):
            self.assertEqual(adj[key], 0, key)
        self.assertEqual((adj["bnd_mate_32bit_buffer_records"], adj["bnd_mate_32bit_buffer_reference_tuples"]), (0, 0))
        self.assertEqual(adj["bnd_bare_reference_allele_rows"], 0)
        surviving = (out / "discordant_adjusted.tsv").read_text().splitlines()[1:]
        self.assertEqual(len(surviving), 3)
        self.assertFalse(any("\t21:50001-60000\t" in l or "\t21:30000001\t" in l or "\t21:70001-70100\t" in l for l in surviving))
        self.assertTrue(
            any("ENST22_X" in l for l in surviving),
            "an other-chromosome transcript row (filter_cross_chromosome_divergences, retired) is charged under 116.2",
        )
        self.assertTrue(any("ENST21_G3" in l for l in surviving), "the batch-class surplus on a kept record is charged")
        md = (out / "concordance_report.md").read_text()
        self.assertIn("- **Reference**: Ensembl VEP 116.2; masks applied: `filter_cnv_tr_expansion_divergences`", md)
        self.assertIn("- **116.2 classes**: `reference_skipped_record` 6 vep-rs rows on 4 records set aside;", md)

    def test_116_2_refuses_a_reference_that_kept_a_record_on_a_skip_route(self) -> None:
        """The 115.2-shaped reference wrote batch rows for a `<CPX>` and a `<NON_REF>`: under
        116.2 those records are on the unsupported-type route, so the cross-check names
        them and the comparison stops instead of scoring over a reference that did not
        skip what 116.2 skips (a run under --dont_skip, or an earlier release)."""
        out, rc, err = self._run("116.2", perl="perl")
        self.assertEqual(rc, 1)
        self.assertIn("ERROR: [compare_sv_concordance] skipped-record cross-check: cpx_kept N><CPX> at 21:45716816 is on the unsupported_type route", err)
        self.assertIn("gvcf_block N><NON_REF> at 21:70001 is on the unsupported_type route", err)
        self.assertFalse((out / "concordance_report.json").exists(), "nothing is scored over a failed cross-check")

    def test_116_2_sets_aside_the_reference_allele_rows_of_a_bare_breakpoint(self) -> None:
        """A bare `<BND>` the reference wrote twice per feature (ALLELE_NUM 1, the reference
        allele N of `N/<BND>`, and 2) and vep-rs once: the two reference-allele rows are set
        aside under `bnd_bare_reference_allele_row`, the tuple sets agree (raw F1 1.0, adjusted
        1.0), and the report names the rows set aside, per file and overall."""
        import json
        import shutil
        import tempfile

        fixture = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, fixture, True)
        for sub in ("in", "perl", "rust"):
            (fixture / sub).mkdir()
        shutil.copytree(self.FIXTURE / "cache", fixture / "cache")
        (fixture / "in" / "bare.vcf").write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "21\t26032000\tbare_bnd\tN\t<BND>\t.\tPASS\tSVTYPE=BND\n"
        )
        header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
        own = "bare_bnd\t21:26032001\tchromosome_breakpoint\tG\tENST21_OWN\tTranscript\tintron_variant\n"
        near = "bare_bnd\t21:26032001\tchromosome_breakpoint\tG\tENST21_G3\tTranscript\tupstream_gene_variant\n"
        (fixture / "perl" / "bare.txt").write_text(header + own + own + near + near)
        (fixture / "rust" / "bare.txt").write_text(header + own + near)
        out = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, out, True)
        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py", "--reference-release", "116.2", "--assembly", "grch37",
            "--input-dir", str(fixture / "in"), "--perl-dir", str(fixture / "perl"),
            "--rust-dir", str(fixture / "rust"), "--vep-rs-cache", str(fixture / "cache"),
            "--output-dir", str(out),
        ]
        try:
            C.main()
        finally:
            sys.argv = argv
        report = json.loads((out / "concordance_report.json").read_text())
        self.assertEqual((report["overall"]["perl_tuples"], report["overall"]["rust_tuples"]), (2, 2))
        self.assertEqual(report["overall"]["f1"], 1.0)
        adj = report["adjusted"]
        self.assertEqual(adj["f1"], 1.0)
        self.assertEqual(adj["bnd_bare_reference_allele_rows"], 2)
        self.assertEqual(adj["bnd_bare_reference_allele_rows_by_file"], {"bare": 2})
        self.assertIn("filter_bnd_bare_reference_allele_rows", adj["masks_applied"])
        cls = adj["divergence_classes"]["bnd_bare_reference_allele_row"]
        self.assertEqual((cls["excludes"], cls["excluded"]), (True, 2))
        self.assertEqual(adj["excluded_perl_extra"], 0, "the rows left the tuple sets before either figure")
        self.assertEqual((adj["perl_tuples"], adj["intersection"]), (2, 2))
        md = (out / "concordance_report.md").read_text()
        self.assertIn("`bnd_bare_reference_allele_row` 2 reference-allele rows set aside", md)

    def test_116_2_exits_0_on_a_dot_breakend_whose_point_lies_inside_the_reference_span(self) -> None:
        """`21 30425650 d3eq_trail NAC NA. . . .`, `21 30425650 d3eq_lead NAC .AC . . .` and
        `21 30425650 d4eq_trail NACG NAC. . . .`: Ensembl VEP 116.2 writes the first two at
        `21:30425650-30425652` and the third at `21:30425650-30425653`, vep-rs each at `21:30425651`,
        as `chromosome_breakpoint` on ENST00000334068 (`upstream_gene_variant`) and ENST00000390690
        (`coding_sequence_variant`). Each record in its own file, the two rows a side of each pair
        under `dot_breakend_without_svtype_pair`, the skipped-record cross-check passes, the
        comparator exits 0 and no one-sided row survives the masks."""
        import json
        import shutil
        import tempfile

        fixture = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, fixture, True)
        for sub in ("in", "perl", "rust"):
            (fixture / sub).mkdir()
        shutil.copytree(self.FIXTURE / "cache", fixture / "cache")
        header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n"
        features = (("ENSG00000186971", "ENST00000334068", "upstream_gene_variant"),
                    ("ENSG00000240432", "ENST00000390690", "coding_sequence_variant"))
        for vid, ref, alt, span in (("d3eq_trail", "NAC", "NA.", "21:30425650-30425652"),
                                    ("d3eq_lead", "NAC", ".AC", "21:30425650-30425652"),
                                    ("d4eq_trail", "NACG", "NAC.", "21:30425650-30425653")):
            (fixture / "in" / f"{vid}.vcf").write_text(
                f"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t30425650\t{vid}\t{ref}\t{alt}\t.\t.\t.\n"
            )
            for side, location in (("perl", span), ("rust", "21:30425651")):
                (fixture / side / f"{vid}.txt").write_text(header + "".join(
                    f"{vid}\t{location}\tchromosome_breakpoint\t{gene}\t{feature}\tTranscript\t{csq}\n"
                    for gene, feature, csq in features))
        out = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, out, True)
        argv = sys.argv
        sys.argv = [
            "compare_sv_concordance.py", "--reference-release", "116.2", "--assembly", "grch38",
            "--input-dir", str(fixture / "in"), "--perl-dir", str(fixture / "perl"),
            "--rust-dir", str(fixture / "rust"), "--vep-rs-cache", str(fixture / "cache"),
            "--output-dir", str(out),
        ]
        rc = 0
        try:
            C.main()
        except SystemExit as exc:
            rc = exc.code or 0
        finally:
            sys.argv = argv
        self.assertEqual(rc, 0)
        report = json.loads((out / "concordance_report.json").read_text())
        self.assertEqual((report["overall"]["perl_tuples"], report["overall"]["rust_tuples"], report["overall"]["intersection"]), (6, 6, 0))
        for name, per_file in report["per_file"].items():
            with self.subTest(file=name):
                self.assertEqual((per_file["perl_tuples"], per_file["rust_tuples"], per_file["intersection"]), (2, 2, 0))
        self.assertIn("filter_dot_breakend_pairs", report["adjusted"]["masks_applied"])
        self.assertTrue(report["adjusted"]["divergence_classes"]["dot_breakend_without_svtype_pair"]["excludes"])
        self.assertEqual((out / "discordant.tsv").read_text().count("\n"), 13, "the header and the twelve one-sided rows")
        self.assertEqual((out / "discordant_adjusted.tsv").read_text().count("\n"), 1, "the header alone: every row paired")

    def test_the_115_2_report_carries_no_release_block(self) -> None:
        import json

        out, _, _ = self._run("115.2")
        adj = json.loads((out / "concordance_report.json").read_text())["adjusted"]
        self.assertNotIn("reference_release", adj)
        self.assertNotIn("masks_applied", adj)
        self.assertNotIn("Reference", (out / "concordance_report.md").read_text())

    def test_the_default_release_is_116_2_and_an_unknown_one_is_refused(self) -> None:
        self.assertEqual(C.DEFAULT_REFERENCE_RELEASE, "116.2")
        self.assertEqual(C.REFERENCE_RELEASES, ("115.2", "116.2"))
        with self.assertRaises(ValueError):
            C.apply_masks(set(), set(), "117.0", {}, {}, {}, "vep-rs")


class SkipRouteCrossCheckTests(unittest.TestCase):
    """`perl_skip_route` is the reference's own drop logic on the input record (Parser/VCF.pm
    116.2:237-245, 468-481, 555-558; Parser.pm:491-501, 663-728; ensembl-io BaseVCF4
    get_start/get_end), and `check_skip_routes` holds the records and the outputs to each
    other."""

    @staticmethod
    def rec(pos: int, ref: str, alt: str, info: str, vid: str = "r") -> "C.VcfVariant":
        return C.VcfVariant("21", pos, vid, ref, alt, C.classify_variant(ref, alt), f"21:{pos}", info)

    def test_the_four_skip_routes_and_the_kept_shapes(self) -> None:
        R = self.rec
        self.assertEqual(C.perl_skip_route(R(18457512, "N", "<BND>", "END=18457513;SVTYPE=BND;SVLEN=14027179;CHR2=21;END2=32484692")), "oversize")
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<DEL>", "END=12000000;SVTYPE=DEL")), "oversize")
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<DEL>", "END=10001001;SVTYPE=DEL")), None, "a span of exactly 10,000,000 is not over the cap")
        # `variant size (Inf) is bigger than --max_sv_size`: Perl reads `1e400` as `Inf`.
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<DUP>", "SVTYPE=DUP;SVLEN=1e400")), "oversize")
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<DUP>", "SVTYPE=DUP;END=1e400")), "oversize")
        self.assertEqual(C.perl_skip_route(R(45716815, "N", "<CPX>", "END=45730272;SVTYPE=CPX")), "unsupported_type")
        self.assertEqual(C.perl_skip_route(R(70000, "N", "<NON_REF>", "END=70100")), "unsupported_type")
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<CTX>", "SVTYPE=CTX")), "unsupported_type")
        self.assertEqual(C.perl_skip_route(R(30165800, "N", "<INS>,<DEL>", "END=30166000")), "unsupported_type")
        self.assertEqual(C.perl_skip_route(R(30165820, "N", "<INV>,<DUP>", "END=30166020")), "unsupported_type")
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<DEL>,T", "END=2000")), "unsupported_type", "a mixed ALT list has no term")
        self.assertEqual(C.perl_skip_route(R(30165850, "N", "<DEL>", "SVTYPE=DEL")), "incomplete_deletion")
        self.assertEqual(C.perl_skip_route(R(2000, "N", "<DEL>", "END=1500;SVTYPE=DEL")), "incomplete_deletion")
        # Perl numifies SVLEN (`abs('abc')` is 0, `SVLEN=` and `SVLEN=.` are defined and 0), so
        # a deletion with any of them ends before its start; a `<DUP>` with them is the pair
        # and kept; an END that is not a number, or one before POS, is refused by
        # `validate_vf` for every class (`invalid_span`), a deletion naming its own route
        # first; a `G,<NON_REF>` site is its own route.
        for info in ("SVTYPE=DEL;SVLEN=abc", "SVTYPE=DEL;SVLEN=", "SVTYPE=DEL;SVLEN=.;END=2000", "SVTYPE=DEL;SVLEN=0;END=2000"):
            with self.subTest(info=info):
                self.assertEqual(C.perl_skip_route(R(1000, "N", "<DEL>", info)), "incomplete_deletion")
        for info in ("SVTYPE=DUP;SVLEN=abc", "SVTYPE=DUP;SVLEN=", "SVTYPE=DUP;SVLEN=.", "SVTYPE=DUP;END=1000"):
            with self.subTest(info=info):
                self.assertIsNone(C.perl_skip_route(R(1000, "N", "<DUP>", info)))
        for alt, info in (("<DUP>", "SVTYPE=DUP;END=abc"), ("<DUP>", "SVTYPE=DUP;END=0"), ("<DUP>", "SVTYPE=DUP;END=900"),
                          ("<BND>", "SVTYPE=BND;END=0"), ("<BND>", "SVTYPE=BND;END="), ("<DUP>", "SVTYPE=DUP;END=")):
            with self.subTest(alt=alt, info=info):
                self.assertEqual(C.perl_skip_route(R(1000, "N", alt, info)), "invalid_span")
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<DEL>", "SVTYPE=DEL;END=abc")), "incomplete_deletion")
        self.assertIsNone(C.perl_skip_route(R(1000, "N", "<DUP>", "SVTYPE=DUP;END=1500.5")), "a fractional END is annotated")
        self.assertEqual(C.perl_skip_route(R(1000, "N", "G,<NON_REF>", ".")), "gvcf_variant_site")
        self.assertEqual(C.perl_skip_route(R(1000, "A", "G,<NON_REF>", ".")), "gvcf_variant_site")
        self.assertIsNone(C.perl_skip_route(R(1000, "N", "G", "SVTYPE=DEL;END=2000")), "an N reference base routes the record to the structural path")
        self.assertIsNone(C.perl_skip_route(R(1000, "N", "<del>", "SVTYPE=DUP;END=2000")), "a lowercase ALT is outside the VCF 4.4 forms, so SVTYPE names the type")
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<del>", ".")), "unsupported_type")
        for alt, info in (("<DEL>,<DUP>", "END=2000"), ("<DEL>", "END=2000;SVTYPE=DEL"), ("<CN0>", "END=2000"),
                          ("<CN0>,<CN2>", "END=2000"), ("<INS:ME:ALU>", "SVLEN=300"), ("<DUP:TANDEM>", "END=5000"),
                          ("<CNV:TR>", "END=1050;SVLEN=50;RUS=AC;RUC=10"), ("N[22:50000000[", "SVTYPE=BND"),
                          ("N[21:50000000[", "SVTYPE=BND;END=50000000"), ("<INS>,<INS>", "SVLEN=300")):
            with self.subTest(alt=alt):
                self.assertIsNone(C.perl_skip_route(R(1000, "N", alt, info)))
        self.assertIsNone(C.perl_skip_route(R(1000, "A", "T", "AF=0.1")), "a sequence variant takes the sequence path")
        self.assertIsNone(C.perl_skip_route(R(1000, "AT", "A", ".")))

    def test_so_abbreviations_follow_get_so_term(self) -> None:
        for sv_type, abbrev in (("<INS:ME:L1>", "INS_LINE1"), ("<DEL:ME:SVA>", "DEL_SVA"), ("<INS:ME>", "INS_ME"),
                                ("<DUP:TANDEM>", "TDUP"), ("<CNV:TR>", "TREP"), ("<CN3>", "CNV"), ("<CN0>", "DEL"),
                                ("<CN2>", "DUP"), ("<CNV>", "CNV"), ("]13:123456]T", "BND"), (".N", "BND"),
                                ("<INV>", "INV"), ("DUP", "DUP"), ("<INS>/<INS>", "INS"), ("<DEL>/<DUP:TANDEM>", "CNV"),
                                ("<DUP:INT>", "DUP")):
            with self.subTest(sv_type=sv_type):
                self.assertEqual(C.perl_so_abbrev(sv_type), abbrev)
        for sv_type in ("<CPX>", "<NON_REF>", "CTX", "<INS>/<DEL>", "<DEL>/T"):
            with self.subTest(sv_type=sv_type):
                self.assertIsNone(C.perl_so_abbrev(sv_type))

    def test_a_trailing_empty_alt_field_is_no_field(self) -> None:
        """Ensembl VEP's parser splits the ALT column in list context (BaseVCF4.pm:489
        `[split(',',$alt_allele)]`), which drops a trailing empty field: `N <DEL>,` with
        `SVTYPE=DEL;END=1100` is the one allele `<DEL>`, a deletion the reference annotates.
        `get_SO_term` normalises a repeated type over the same list-context split, on `/` or
        `,` (Parser.pm:670), so `DEL/` is the one part `DEL` and `SVTYPE=DEL/` names a deletion.
        Every reader of the ALT column sees the kept fields: `<DEL>,` is one symbolic deletion
        to `classify_variant`, `record_alleles` and `record_location`, `N.,` a breakend whose
        END is discarded, `<CNV:TR>,` a single-ALT tandem repeat, `AG,<NON_REF>,` the gVCF site
        `AG,<NON_REF>` and a bracket mate list with a trailing field still names its mate."""
        R = self.rec
        self.assertIsNone(C.perl_skip_route(R(1000, "N", "<DEL>,", "SVTYPE=DEL;END=1100")))
        self.assertIsNone(C.perl_skip_route(R(1000, "N", "N.", "SVTYPE=DEL/;END=1100")))
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<DEL>,,", "END=1100")), None, "two trailing empty fields are none")
        self.assertEqual(C.perl_skip_route(R(1000, "N", "<DEL>,,T", "END=1100")), "unsupported_type", "a middle empty field is kept")
        self.assertEqual(C.perl_so_abbrev("<DEL>/"), "DEL")
        self.assertEqual(C.perl_so_abbrev("DEL/"), "DEL")
        self.assertEqual(C.perl_so_abbrev("<DEL>//<DUP>"), "CNV")
        self.assertEqual(C.classify_variant("N", "<DEL>,"), "Symbolic_DEL")
        self.assertEqual(C.classify_variant("N", "<DEL>,<DUP>"), "Multi_allelic")
        self.assertEqual(C.record_alleles(R(1000, "N", "<DEL>,", "END=2000;SVTYPE=DEL")), {"deletion"})
        self.assertEqual(C.record_location(R(1000, "N", "N.,", "END=2000")), ("21", 1001, 1001))
        self.assertEqual(C.record_location(R(1000, "N", "<DEL>,", "END=2000")), ("21", 1001, 2000))
        run = "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=3"
        self.assertTrue(C.cnv_tr_expansion_is_empty(R(30425587, "N", "<CNV:TR>,", run)))
        rust = {T("21:30425651", "-", "ENST_1", "frameshift_variant")}
        self.assertEqual(C.filter_gvcf_variant_sites(set(), rust, [R(30425650, "ATG", "AG,<NON_REF>,", ".")]), (rust, 1))
        self.assertEqual(C.filter_gvcf_variant_sites(set(), rust, [R(30425650, "ATG", "AG,", ".")]), (set(), 0))
        records, rows = C.count_bnd_mate_over_32_bit([R(1000, "N", "N[21:2147483647[,", "SVTYPE=BND")], {T("21:1001", "N[21:2147483647[", "ENST_1", "intron_variant")})
        self.assertEqual((records, rows), (1, 1))

    def test_an_absent_location_no_skip_route_record_explains_is_a_problem(self) -> None:
        """vep-rs structural rows at a Location the reference wrote nothing for, with no
        skip-route record starting there: a record dropped on a route the registry does not
        know (`start > end+1`), or a record missing from the input."""
        kept = self.rec(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL", vid="kept")
        rust = {T("21:1001-2000", "deletion", "ENST_D", "feature_truncation"), T("21:9999", "deletion", "ENST_X", "intron_variant")}
        problems = C.check_skip_routes([kept], set(), rust)
        self.assertEqual(len(problems), 2)
        self.assertTrue(all("no input record on a skip route starts there" in p for p in problems))

    def test_a_skip_route_record_with_reference_rows_under_its_own_allele_is_a_problem(self) -> None:
        cpx = self.rec(45716815, "N", "<CPX>", "END=45730272;SVTYPE=CPX", vid="cpx")
        perl = {T("21:45716816-45730272", "CPX", "ENST_C1", "intron_variant")}
        problems = C.check_skip_routes([cpx], perl, set())
        self.assertEqual(len(problems), 1)
        self.assertIn("cpx N><CPX> at 21:45716816 is on the unsupported_type route, yet the reference wrote 1 row(s) under its own allele(s) (CPX)", problems[0])

    def test_reference_rows_of_a_kept_record_at_the_same_start_are_not_the_skipped_records(self) -> None:
        cpx = self.rec(16450356, "N", "<CPX>", "END=16456471;SVTYPE=CPX", vid="cpx")
        dele = self.rec(16450356, "N", "<DEL>", "END=16456471;SVTYPE=DEL", vid="del")
        perl = {T("21:16450357-16456471", "deletion", "ENST_D", "feature_truncation")}
        self.assertEqual(C.check_skip_routes([cpx, dele], perl, perl), [])

    def test_twins_sharing_start_end_and_allele_are_a_problem(self) -> None:
        """A `<DEL>,T` record (no Sequence Ontology term, skipped) and a kept `<DEL>` at the
        same start both produce a `deletion` row; the derived end tells them apart, and
        when it cannot the comparison stops rather than guess whose rows are whose."""
        mixed = self.rec(1000, "N", "<DEL>,T", "END=2000", vid="mixed")
        kept_shorter = self.rec(1000, "N", "<DEL>", "END=2000;SVLEN=-500;SVTYPE=DEL", vid="kept_shorter")
        perl = {T("21:1001-1500", "deletion", "ENST_D", "feature_truncation")}
        rust = perl | {T("21:1001-2000", "deletion", "ENST_E", "transcript_ablation"), T("21:1001-2000", "T", "ENST_E", "missense_variant")}
        self.assertEqual(C.check_skip_routes([mixed, kept_shorter], perl, rust), [])
        self.assertEqual(C.filter_reference_skipped_records(perl, rust, [mixed, kept_shorter]), (rust - perl, 1))
        twin = self.rec(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL", vid="twin")
        twin_perl = {T("21:1001-2000", "deletion", "ENST_D", "feature_truncation")}
        problems = C.check_skip_routes([mixed, twin], twin_perl, twin_perl | {T("21:1001-2000", "T", "ENST_D", "missense_variant")})
        self.assertEqual(problems, ["mixed N><DEL>,T and twin N><DEL> share 21:1001-2000 and an allele; the skipped one cannot be told from the kept one"])

    def test_a_kept_sequence_record_carrying_svtype_and_end_is_not_the_twin_of_a_skipped_record(self) -> None:
        """`21 30425650 x A G . . SVTYPE=DEL;END=30425949` is a sequence record to both engines (REF
        and ALT a run of `ACGT`), written by both at `21:30425650 G ENST00000390690 missense_variant`
        and `21:30425650 G ENST00000334068 upstream_gene_variant`; `21 30425650 x N G . .
        SVTYPE=DEL ;END=30425949` is skipped by both (`DEL  is not a supported structural variant
        type`). The sequence record is keyed at POS under its own allele, not at the span
        `21:30425651-30425949` its INFO would give a structural record, so the skipped record is
        confirmed and the cross-check holds, beside the kept `21 30425650 x N <DEL> . .
        SVTYPE=DEL;END=30425949` (`21:30425651-30425949 deletion` on both transcripts) as well.
        A sequence record still contests a skipped record at its own start with its own end: a
        kept `A T` at 1001 beside a skipped `<DEL>,T` over 1001-2000 is told apart by the end,
        and beside one whose END is 1001 the two cannot be."""
        kept = REC(30425650, "A", "G", "SVTYPE=DEL;END=30425949", "seq")
        skipped = REC(30425650, "N", "G", "SVTYPE=DEL ;END=30425949", "space")
        control = REC(30425650, "N", "<DEL>", "SVTYPE=DEL;END=30425949", "del")
        rows = {T("21:30425650", "G", "ENST00000390690", "missense_variant"),
                T("21:30425650", "G", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425651-30425949", "deletion", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425651-30425949", "deletion", "ENST00000390690",
                  "feature_truncation,coding_sequence_variant,5_prime_UTR_variant")}
        self.assertTrue(C.perl_sequence_path(kept))
        self.assertFalse(C.perl_sequence_path(skipped))
        self.assertEqual(C.perl_skip_route(skipped), "unsupported_type")
        self.assertEqual(C.check_skip_routes([kept, skipped, control], rows, rows), [])
        self.assertEqual(C.filter_reference_skipped_records(rows, rows, [kept, skipped, control]), (set(), 1))
        perl_open, rust_open, masks = ADJUSTED(rows, rows, [kept, skipped, control])
        self.assertEqual((perl_open, rust_open, masks.skipped_records, masks.skipped_problems), (set(), set(), 1, []))
        snv = REC(1001, "A", "T", ".", "snv")
        mixed = REC(1000, "N", "<DEL>,T", "END=2000", "mixed")
        shared = {T("21:1001", "T", "ENST_D", "missense_variant")}
        mixed_rows = {T("21:1001-2000", "deletion", "ENST_E", "transcript_ablation"), T("21:1001-2000", "T", "ENST_E", "missense_variant")}
        self.assertEqual(C.check_skip_routes([snv, mixed], shared, shared | mixed_rows), [])
        self.assertEqual(C.filter_reference_skipped_records(shared, shared | mixed_rows, [snv, mixed]), (mixed_rows, 1))
        twin = REC(1000, "N", "<DEL>,T", "END=1001", "twin")
        self.assertEqual(C.check_skip_routes([snv, twin], shared, shared | {T("21:1001", "deletion", "ENST_D", "feature_truncation")}),
                         ["twin N><DEL>,T and snv A>T share 21:1001-1001 and an allele; the skipped one cannot be told from the kept one"])

    def test_an_infinite_span_is_oversize_and_its_pair_rows_are_set_aside(self) -> None:
        """Ensembl VEP 116.2 drops `21 30425650 x N <INS> . . SVTYPE=INS;SVLEN=Inf` as oversize
        (`variant size (Inf) is bigger than --max_sv_size (10000000)`), as it drops
        `21 30425650 x N <DUP:TANDEM> . . SVTYPE=DUP;SVLEN=Inf` and `21 30425587 x N <CNV:TR> . .
        SVTYPE=CNV;SVLEN=Inf;RUS=AGC;RUC=4` in every run; vep-rs reads the length as 0 and writes
        the pair `21:30425650-30425651 insertion ENST00000390690 coding_sequence_variant,feature_elongation`
        with `21:30425650-30425651 insertion ENST00000334068 upstream_gene_variant`
        (`tandem_duplication` for the second, the pair `21:30425587-30425588 tandem_repeat` for
        the third), set aside under the class. `SVLEN=NaN` and `END=NaN` keep a record (every
        comparison on NaN is false; the reference writes `21:30425651-NaN`), `END=Inf` is
        oversize, and a breakend's `END=Inf` is discarded with its END."""
        cases = [
            (REC(30425650, "N", "<INS>", "SVTYPE=INS;SVLEN=Inf", "ins"),
             {T("21:30425650-30425651", "insertion", "ENST00000390690", "coding_sequence_variant,feature_elongation"),
              T("21:30425650-30425651", "insertion", "ENST00000334068", "upstream_gene_variant")}),
            (REC(30425650, "N", "<DUP:TANDEM>", "SVTYPE=DUP;SVLEN=Inf", "tdup"),
             {T("21:30425650-30425651", "tandem_duplication", "ENST00000390690", "coding_sequence_variant,feature_elongation"),
              T("21:30425650-30425651", "tandem_duplication", "ENST00000334068", "upstream_gene_variant")}),
            (REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;SVLEN=Inf;RUS=AGC;RUC=4", "tr"),
             {T("21:30425587-30425588", "tandem_repeat", "ENST00000390690", "coding_sequence_variant,feature_elongation"),
              T("21:30425587-30425588", "tandem_repeat", "ENST00000334068", "upstream_gene_variant")}),
        ]
        for record, rust in cases:
            with self.subTest(alt=record.alt):
                self.assertEqual(C.perl_skip_route(record), "oversize")
                self.assertIsNone(C.perl_skip_route(record, -1), "without a cap the reference run dies on the span")
                self.assertEqual(C.filter_reference_skipped_records(set(), rust, [record]), (rust, 1))
                self.assertEqual(C.check_skip_routes([record], set(), rust), [])
        self.assertIsNone(C.perl_skip_route(REC(30425650, "N", "<DEL>", "SVTYPE=DEL;SVLEN=NaN")))
        self.assertIsNone(C.perl_skip_route(REC(30425650, "N", "<DEL>", "SVTYPE=DEL;END=NaN")))
        self.assertEqual(C.perl_skip_route(REC(30425650, "N", "<DEL>", "SVTYPE=DEL;END=Inf")), "oversize")
        self.assertEqual(C.perl_skip_route(REC(30425650, "N", "<DEL>", "SVTYPE=DEL;SVLEN=-Infinity")), "oversize")
        self.assertIsNone(C.perl_skip_route(REC(30425650, "N", "<BND>", "SVTYPE=BND;END=Inf")))

    def test_agreement_yields_no_problem(self) -> None:
        cpx = self.rec(45716815, "N", "<CPX>", "END=45730272;SVTYPE=CPX", vid="cpx")
        kept = self.rec(1000, "N", "<DEL>", "END=2000;SVTYPE=DEL", vid="kept")
        perl = {T("21:1001-2000", "deletion", "ENST_D", "feature_truncation")}
        rust = perl | {T("21:45716816-45730272", "CPX", "ENST_C", "intron_variant")}
        self.assertEqual(C.check_skip_routes([cpx, kept], perl, rust), [])


def VCF(*records: tuple[str, ...]) -> list[C.VcfVariant]:
    """The records of a VCF written from its eight-column data lines, read as `parse_vcf` reads them."""
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "in.vcf"
        path.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "".join("\t".join(record) + "\n" for record in records),
            encoding="utf-8",
        )
        return C.parse_vcf(path)


def ADJUSTED(perl: set, rust: set, records: list, max_sv_size: int = C._MAX_SV_SIZE) -> tuple[set, set, C.MaskResult]:
    """The one-sided tuples every 116.2 mask leaves in place, each side, with the mask result."""
    masks = C.apply_masks(perl, rust, "116.2", {}, {}, {}, "vep-rs", records, max_sv_size)
    return (perl - rust) - masks.all_perl(), (rust - perl) - masks.all_rust(), masks


class NumifiedPosTests(unittest.TestCase):
    """`parse_vcf` reads POS as Ensembl VEP 116.2 numifies it (BaseVCF4.pm:263 `$start++` on
    the structural path), and the class over a POS the reference annotated where vep-rs
    refused the line."""

    def test_a_pos_in_another_spelling_is_read_as_its_number(self) -> None:
        """`21 3e7 x N <DEL> . . SVTYPE=DEL;END=30425949` is the record at 30000000 (Ensembl VEP
        116.2 writes `21:30000001-30425949 deletion ENST00000286808 transcript_ablation`);
        `30425949.0`, ` 30425949`, `30425949 `, `+30425949` and `030425949` are 30425949;
        `abc`, `.` and fullwidth digits are 0; a digit string past 2^63 is exact."""
        records = VCF(
            ("21", "3e7", "exp", "N", "<DEL>", ".", ".", "SVTYPE=DEL;END=30425949"),
            ("21", "30425949.0", "float0", "A", "G", ".", ".", "."),
            ("21", " 30425949", "lead", "A", "G", ".", ".", "."),
            ("21", "30425949 ", "trail", "A", "G", ".", ".", "."),
            ("21", "+30425949", "plus", "N", "<DEL>", ".", ".", "SVTYPE=DEL;END=30425949"),
            ("21", "030425949", "zeropad", "N", "<DEL>", ".", ".", "SVTYPE=DEL;END=30425949"),
            ("21", "abc", "abc", "A", "G", ".", ".", "."),
            ("21", ".", "dot", "N", "<DEL>", ".", ".", "SVTYPE=DEL;END=30425949"),
            ("21", "３０４２５９４９", "fullwidth", "N", "<DEL>", ".", ".", "SVTYPE=DEL;END=30425949"),
            ("21", "18446744073709551616", "huge", "N", "<DEL>", ".", ".", "SVTYPE=DEL;SVLEN=-300"),
        )
        self.assertEqual([v.pos for v in records], [30000000, 30425949, 30425949, 30425949, 30425949, 30425949, 0, 0, 0, 2**64])
        self.assertEqual([v.pos_text for v in records[:3]], ["3e7", "30425949.0", " 30425949"])
        self.assertEqual((records[0].location, records[0].variant_type), ("21:30000000", "Symbolic_DEL"))
        for text, value in (("3e7", 30000000), ("30425949.", 30425949), ("3.0425949E7", 30425949), ("1e3", 1000),
                            ("-5", -5), ("30_425_949", 30), ("0x1D04A1D", 0), ("NaN", 0), ("1e400", 2**63 - 1),
                            ("4611686018427387904", 2**62)):
            self.assertEqual(C.perl_pos(text), value, text)
        for text in ("30425949", "+30425949", "030425949", "3e7", "30425949.0", "30425949.", " 30425949", "30425949 ",
                     "3.0425949E7", "1e3", None):
            self.assertTrue(C.perl_pos_is_a_number(text), text)
        for text in (".", "", "0x1D04A1D", "３０４２５９４９", "30_425_949", "-5", "-2147483648", "NaN", "Inf", "e5",
                     "30425949,30425950"):
            self.assertFalse(C.perl_pos_is_a_number(text), text)

    def test_the_reference_rows_of_a_numified_pos_are_set_aside_from_the_input(self) -> None:
        """Ensembl VEP 116.2 numifies a POS that is not a number: under --max_sv_size -1
        `21 . x N <DEL> . . SVTYPE=DEL;END=30425949` is annotated from base 1 (`21:1-30425949
        deletion ENST00000262354 transcript_ablation` and a row for every other transcript the
        span reaches, every reference row of the record set aside) and `21 30_425_949 x N <DEL>
        . . SVTYPE=DEL;END=30425949` from 31 (`21:31-30425949 deletion ...`); `21 -5 x N <DEL>
        . . SVTYPE=DEL;SVLEN=-300` is `21:-4-295 deletion - intergenic_variant`, `21 -2147483648
        x N <DEL> . . SVTYPE=DEL;SVLEN=-300` `21:-2147483647--2147483348 deletion -
        intergenic_variant` and the sequence record `21 -5 x A G` `21:-5 G - intergenic_variant`;
        vep-rs refuses every line. `21 0 x N <DEL> . . SVTYPE=DEL;SVLEN=-300` is `21:1-300` in
        both engines and not the class, and a negative coordinate is an integer, outside
        `non_integer_end_annotated`."""
        records = VCF(
            ("21", ".", "dot", "N", "<DEL>", ".", ".", "SVTYPE=DEL;END=30425949"),
            ("21", "30_425_949", "underscore", "N", "<DEL>", ".", ".", "SVTYPE=DEL;END=30425949"),
            ("21", "-5", "neg", "N", "<DEL>", ".", ".", "SVTYPE=DEL;SVLEN=-300"),
            ("21", "-2147483648", "negbig", "N", "<DEL>", ".", ".", "SVTYPE=DEL;SVLEN=-300"),
            ("21", "-5", "neg_snv", "A", "G", ".", ".", "."),
            ("21", "0", "zero", "N", "<DEL>", ".", ".", "SVTYPE=DEL;SVLEN=-300"),
        )
        numified = {
            T("21:1-30425949", "deletion", "ENST00000262354", "transcript_ablation"),
            T("21:1-30425949", "deletion", "ENST00000282964", "transcript_ablation"),
            T("21:31-30425949", "deletion", "ENST00000262354", "transcript_ablation"),
            INTERGENIC("21:-4-295", "deletion"),
            INTERGENIC("21:-2147483647--2147483348", "deletion"),
            INTERGENIC("21:-5", "G"),
        }
        shared = {INTERGENIC("21:1-300", "deletion")}
        self.assertEqual(C.filter_non_numeric_pos_numified(numified | shared, shared, records), numified)
        self.assertEqual(C.filter_non_integer_end_rows(numified | shared, shared, records), (set(), set()))
        perl_open, rust_open, masks = ADJUSTED(numified | shared, shared, records, -1)
        self.assertEqual((perl_open, rust_open, masks.numified_pos_perl, masks.skipped_problems), (set(), set(), numified, []))


class StarAloneSvtypeTests(unittest.TestCase):
    def test_the_reference_rows_of_a_star_only_alt_under_a_true_svtype_are_set_aside(self) -> None:
        """Ensembl VEP 116.2 annotates `21 30425650 x N * . . SVTYPE=DEL;END=30425949` as
        `21:30425651-30425949 deletion ENST00000334068 upstream_gene_variant` and
        `21:30425651-30425949 deletion ENST00000390690 feature_truncation,coding_sequence_variant,
        5_prime_UTR_variant`, `A *` the same, and `21 30425650 x N * . . SVTYPE=BND` as
        `21:30425651 chromosome_breakpoint` on both transcripts (ALLELE_NUM 1 and 2, one tuple
        each); vep-rs writes no row (`No alternate allele`). `21 30425650 x N * . . .` has no
        reference row, and a record with a real ALT at the Location is not touched."""
        records = [REC(30425650, "N", "*", "SVTYPE=DEL;END=30425949", "star"),
                   REC(30425650, "A", "*", "SVTYPE=DEL;END=30425949", "star_a"),
                   REC(30425650, "N", "*", "SVTYPE=BND", "star_bnd"),
                   REC(30425650, "N", "*,*", "SVTYPE=DEL;END=30425949", "two_stars"),
                   REC(30425650, "N", "*", ".", "bare"),
                   REC(11000000, "N", "<DEL>", "SVTYPE=DEL;END=11000299", "kept")]
        deletion = {T("21:30425651-30425949", "deletion", "ENST00000334068", "upstream_gene_variant"),
                    T("21:30425651-30425949", "deletion", "ENST00000390690",
                      "feature_truncation,coding_sequence_variant,5_prime_UTR_variant")}
        breakpoint = {T("21:30425651", "chromosome_breakpoint", "ENST00000334068", "upstream_gene_variant"),
                      T("21:30425651", "chromosome_breakpoint", "ENST00000390690", "coding_sequence_variant")}
        kept = {INTERGENIC("21:11000001-11000299", "deletion")}
        self.assertEqual(C.filter_star_alone_svtype_annotated(deletion | breakpoint | kept, records), deletion | breakpoint)
        self.assertEqual(C.filter_star_alone_svtype_annotated(deletion, records[4:]), set())
        perl_open, rust_open, masks = ADJUSTED(deletion | breakpoint | kept, kept, records)
        self.assertEqual((perl_open, rust_open, masks.star_alone_perl), (set(), set(), deletion | breakpoint))


class UnsupportedAltAlleleTests(unittest.TestCase):
    """vep-rs annotates a symbolic ALT whose type Ensembl VEP 116.2 has no term for under the
    ALT's own class (`classify_symbolic_alt`, `display_allele`); `record_alleles` names that
    spelling so the record's rows are set aside under `reference_skipped_record`."""

    def test_the_class_term_vep_rs_prints_for_an_unsupported_alt(self) -> None:
        for alt, term in (("<del>", "deletion"), ("<Del>", "deletion"), ("<dEL>", "deletion"), ("<DEL:>", "deletion"),
                          ("<DEL >", "deletion"), ("<DEL\u00a0>", "deletion"), ("<DEL\u200b>", "deletion"),
                          ("<Ins>", "insertion"), ("<INS:>", "insertion"), ("<Bnd>", "chromosome_breakpoint"),
                          ("<CN=X>", "copy_number_variation"), ("<CNX>", "copy_number_variation"),
                          ("<CN>", "copy_number_variation"), ("<Dup>", "duplication"), ("<Inv>", "inversion"),
                          ("<<DEL>>", "<<DEL>>"), ("<:DEL>", "<:DEL>"), ("<CPX>", "<CPX>"), ("<CPX:X>", "<CPX:X>"),
                          ("<DÉL>", "<DÉL>"), ("<ＤＥＬ>", "<ＤＥＬ>"), ("< DEL >", "< DEL >")):
            self.assertEqual(C.vep_rs_unsupported_alt_allele(alt, None), term, alt)
        self.assertEqual(C.vep_rs_unsupported_alt_allele("<FOO>", "DEL"), "deletion", "INFO/SVTYPE names the class of a type the ALT does not")
        self.assertEqual(C.vep_rs_unsupported_alt_allele("<FOO>", "foo"), "foo", "INFO/SVTYPE is the type string when it names no class")
        self.assertEqual(C.vep_rs_unsupported_alt_allele("<CPX>", "CPX"), "CPX")
        self.assertEqual(C.vep_rs_unsupported_alt_allele("<FOO>", "0"), "<FOO>", "a false SVTYPE names nothing, as Perl reads `0`")
        self.assertEqual(C.vep_rs_unsupported_alt_allele("<DEL:X>", "CPX"), "deletion", "a VCF 4.4 form keeps its own class")
        self.assertEqual(C.record_alleles(REC(30425650, "N", "<Del>", "END=30425949")), {"Del", "<Del>", "deletion"})
        self.assertEqual(C.record_alleles(REC(30425650, "N", "<CN=X>", "END=30425949")), {"CN=X", "<CN=X>", "copy_number_variation"})
        self.assertEqual(C.record_alleles(REC(30425650, "N", "<:DEL>", "END=30425949")), {"", "<>", "<:DEL>", "-"},
                         "the tab format writes the empty Allele as `-`")
        self.assertEqual(C.record_alleles(REC(30425650, "N", "<DEL>,<FOO>", "END=30425949")), {"deletion", "FOO", "<FOO>", "<DEL>/<FOO>"},
                         "the typeless ALT's rows print the list as written, the DEL's its class")
        self.assertEqual(C.record_alleles(REC(30425650, "N", "<FOO>,<BAR>", "SVTYPE=baz;END=30425949")), {"FOO", "<FOO>", "BAR", "<BAR>", "baz"},
                         "INFO/SVTYPE names a list that opens outside the VCF 4.4 forms")
        self.assertEqual(C.record_alleles(REC(30425650, "N", "<<DEL>>,<DEL>", "END=30425949")),
                         {"deletion", "<<DEL>>", "<<DEL>>/<DEL>"},
                         "`<<DEL>>` alone is a DEL to the reference, the list `DEL/DEL` has no term, and vep-rs prints the list as written on the doubled ALT's rows")
        self.assertEqual(C.record_alleles(REC(30425650, "N", "<DEL>,T", "SVTYPE=INS;END=30425949")), {"deletion", "insertion", "T"},
                         "a mixed list under a type the table lacks takes INFO/SVTYPE's class in vep-rs")
        self.assertEqual(C.record_alleles(REC(1000, "N", "<CPX>", "END=2000;SVTYPE=CPX")), {"CPX", "<CPX>"}, "the raw spelling stays")

    def test_the_rows_of_an_unsupported_alt_are_set_aside_under_the_class_term(self) -> None:
        """`21 30425650 x N <Del> . . END=30425949` is skipped by Ensembl VEP 116.2 (`Del is not a
        supported structural variant type`); vep-rs writes `21:30425651-30425949 deletion
        ENST00000334068 upstream_gene_variant` and `21:30425651-30425949 deletion ENST00000390690
        feature_truncation,coding_sequence_variant,5_prime_UTR_variant`, `<CN=X>` the same rows
        under `copy_number_variation` (`coding_sequence_variant,5_prime_UTR_variant` on the
        coding transcript), `<Bnd>` `21:30425651 chromosome_breakpoint` rows, `<<DEL>>,<DEL>`
        a `<<DEL>>/<DEL>` row (the list as written) beside a `deletion` row per transcript, and
        `21 30425650 x N <:DEL> . . END=30425949` the ALT as written in the default and tab
        formats (`21:30425651-30425949 <:DEL> ENST00000390690
        5_prime_UTR_variant,coding_sequence_variant`); every row is set aside and the
        cross-check holds."""
        span = "21:30425651-30425949"
        cases = [
            (REC(30425650, "N", "<:DEL>", "END=30425949", "colon"),
             {T(span, "<:DEL>", "ENST00000334068", "upstream_gene_variant"),
              T(span, "<:DEL>", "ENST00000390690", "coding_sequence_variant,5_prime_UTR_variant")}),
            (REC(30425650, "N", "<Del>", "END=30425949", "lower"),
             {T(span, "deletion", "ENST00000334068", "upstream_gene_variant"),
              T(span, "deletion", "ENST00000390690", "feature_truncation,coding_sequence_variant,5_prime_UTR_variant")}),
            (REC(30425650, "N", "<DEL >", "SVTYPE=DEL;END=30425949", "blank"),
             {T(span, "deletion", "ENST00000334068", "upstream_gene_variant"),
              T(span, "deletion", "ENST00000390690", "feature_truncation,coding_sequence_variant,5_prime_UTR_variant")}),
            (REC(30425650, "N", "<CN=X>", "END=30425949", "cn"),
             {T(span, "copy_number_variation", "ENST00000334068", "upstream_gene_variant"),
              T(span, "copy_number_variation", "ENST00000390690", "coding_sequence_variant,5_prime_UTR_variant")}),
            (REC(30425650, "N", "<Bnd>", ".", "bnd"),
             {T("21:30425651", "chromosome_breakpoint", "ENST00000334068", "upstream_gene_variant"),
              T("21:30425651", "chromosome_breakpoint", "ENST00000390690", "coding_sequence_variant")}),
            (REC(30425650, "N", "<<DEL>>,<DEL>", "END=30425949", "doubled"),
             {T(span, "<<DEL>>/<DEL>", "ENST00000334068", "upstream_gene_variant"),
              T(span, "deletion", "ENST00000334068", "upstream_gene_variant"),
              T(span, "<<DEL>>/<DEL>", "ENST00000390690", "coding_sequence_variant,5_prime_UTR_variant"),
              T(span, "deletion", "ENST00000390690", "feature_truncation,coding_sequence_variant,5_prime_UTR_variant")}),
        ]
        for record, rust in cases:
            with self.subTest(alt=record.alt):
                self.assertEqual(C.perl_skip_route(record), "unsupported_type")
                perl_open, rust_open, masks = ADJUSTED(set(), rust, [record])
                self.assertEqual((perl_open, rust_open, masks.skipped_rust, masks.skipped_records, masks.skipped_problems),
                                 (set(), set(), rust, 1, []))


class IncompleteDeletionPairTests(unittest.TestCase):
    def test_an_incomplete_deletions_pair_rows_are_set_aside(self) -> None:
        """Ensembl VEP 116.2 drops `21 30425650 x N <DEL> . . SVTYPE=DEL;END=30425650` (`deletion
        looks incomplete`, Parser/VCF.pm:555-558), as it drops `SVLEN=0`, `SVLEN=abc`, `SVLEN=.`,
        `SVLEN=`, `SVLEN=e5`, `SVLEN=0x1D04A1D` and each with `END=30425949` beside it; vep-rs
        writes the pair `21:30425650-30425651 deletion ENST00000334068 upstream_gene_variant` and
        `21:30425650-30425651 deletion ENST00000390690 coding_sequence_variant`, printed low end
        first. `21 +30425949 x N <DEL> . . SVTYPE=DEL;END=30425949` is the same shape at
        `21:30425949-30425950`."""
        pair = {T("21:30425650-30425651", "deletion", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425650-30425651", "deletion", "ENST00000390690", "coding_sequence_variant")}
        for info in ("SVTYPE=DEL;END=30425650", "SVTYPE=DEL;SVLEN=0", "SVTYPE=DEL;SVLEN=abc", "SVTYPE=DEL;SVLEN=.",
                     "SVTYPE=DEL;SVLEN=", "SVTYPE=DEL;SVLEN=e5", "SVTYPE=DEL;SVLEN=0x1D04A1D", "SVTYPE=DEL;SVLEN=.;END=30425949",
                     "SVTYPE=DEL;SVLEN=;END=30425949", "SVTYPE=DEL;SVLEN=abc;END=30425949"):
            with self.subTest(info=info):
                record = REC(30425650, "N", "<DEL>", info)
                self.assertEqual(C.perl_skip_route(record), "incomplete_deletion")
                self.assertEqual(C.filter_reference_skipped_records(set(), pair, [record]), (pair, 1))
                self.assertEqual(C.check_skip_routes([record], set(), pair), [])
        (plus,) = VCF(("21", "+30425949", "plus", "N", "<DEL>", ".", ".", "SVTYPE=DEL;END=30425949"))
        plus_pair = {T("21:30425949-30425950", "deletion", "ENST00000334068", "upstream_gene_variant"),
                     T("21:30425949-30425950", "deletion", "ENST00000390690", "upstream_gene_variant")}
        self.assertEqual(C.filter_reference_skipped_records(set(), plus_pair, [plus]), (plus_pair, 1))
        self.assertEqual(C.check_skip_routes([plus], set(), plus_pair), [])


class DotBreakendInfoMateTests(unittest.TestCase):
    def test_a_dot_breakend_with_an_info_mate_pairs_under_its_own_and_mate_alleles(self) -> None:
        """`21 30425650 x N N. . . CHR2=21;END2=30425700`: Ensembl VEP 116.2 writes the pair
        `21:30425650-30425651` under `N.` and under `N[21:30425700[` on ENST00000334068
        (`upstream_gene_variant`) and ENST00000390690 (`feature_truncation,coding_sequence_variant`);
        vep-rs writes the point `21:30425651` under the same alleles. The eight rows are set
        aside as pairs and the skipped-record cross-check holds, as it does for the bare
        `21 30425650 x N N. . . .` (`chromosome_breakpoint` on both sides)."""
        perl = {T("21:30425650-30425651", allele, feature, csq)
                for allele in ("N.", "N[21:30425700[")
                for feature, csq in (("ENST00000334068", "upstream_gene_variant"),
                                     ("ENST00000390690", "feature_truncation,coding_sequence_variant"))}
        rust = {T("21:30425651", t[1], t[2], t[4]) for t in perl}
        self.assertEqual(C.filter_dot_breakend_pairs(perl, rust), (rust, perl))
        record = REC(30425650, "N", "N.", "CHR2=21;END2=30425700")
        self.assertEqual(C.check_skip_routes([record], perl, rust), [])
        perl_open, rust_open, masks = ADJUSTED(perl, rust, [record])
        self.assertEqual((perl_open, rust_open, len(masks.dot_pair_perl), len(masks.dot_pair_rust)), (set(), set(), 4, 4))
        bare_perl = {T("21:30425650-30425651", "chromosome_breakpoint", "ENST00000334068", "upstream_gene_variant"),
                     T("21:30425650-30425651", "chromosome_breakpoint", "ENST00000390690", "coding_sequence_variant")}
        bare_rust = {T("21:30425651", t[1], t[2], t[4]) for t in bare_perl}
        self.assertEqual(C.check_skip_routes([REC(30425650, "N", "N.", ".")], bare_perl, bare_rust), [])
        self.assertEqual(C.filter_dot_breakend_pairs(bare_perl, bare_rust), (bare_rust, bare_perl))
        sequence = {T("21:30425650-30425651", "ACG", "ENST00000390690", "inframe_insertion")}
        self.assertEqual(C.filter_dot_breakend_pairs(sequence, {T("21:30425651", "ACG", "ENST00000390690", "inframe_insertion")}), (set(), set()),
                         "a sequence allele is not a breakend")

    def test_a_dot_breakend_pairs_at_the_spans_start_or_the_base_after_it(self) -> None:
        """The reference's start is POS + 1 when the ALT's length differs from the REF's and POS
        when it does not, its end POS + len(REF) - 1 (BaseVCF4.pm:264-275, 297-316); vep-rs writes
        the point `21:30425651 chromosome_breakpoint` on both transcripts for every shape. Ensembl
        VEP 116.2 writes `21 30425650 x NAC .N . . .` as `21:30425651-30425652 chromosome_breakpoint
        ENST00000334068 upstream_gene_variant` and `21:30425651-30425652 chromosome_breakpoint
        ENST00000390690 coding_sequence_variant`, `21 30425650 x NACG .N . . .` the same rows at
        `21:30425651-30425653` and `21 30425650 x NACGT N. . . .` at `21:30425651-30425654`, the
        point each span's start; `21 30425650 x NAC NA. . . .` and `21 30425650 x NAC .AC . . .` at
        `21:30425650-30425652` and `21 30425650 x NACG NAC. . . .` at `21:30425650-30425653`, the
        point the base after each span's start; `21 30425650 x N N. . . .` and `21 30425650 x NA N.
        . . .` at `21:30425650-30425651`, the point its end. Each pairs and the skipped-record
        cross-check holds. A point two or more bases past the span's start pairs with nothing: no
        REF length produces that shape."""
        rust = {T("21:30425651", "chromosome_breakpoint", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425651", "chromosome_breakpoint", "ENST00000390690", "coding_sequence_variant")}
        for ref, alt, span in (("NAC", ".N", "21:30425651-30425652"), ("NACG", ".N", "21:30425651-30425653"),
                               ("NACGT", "N.", "21:30425651-30425654"), ("NAC", "NA.", "21:30425650-30425652"),
                               ("NAC", ".AC", "21:30425650-30425652"), ("NACG", "NAC.", "21:30425650-30425653"),
                               ("N", "N.", "21:30425650-30425651"), ("NA", "N.", "21:30425650-30425651")):
            with self.subTest(ref=ref, alt=alt):
                record = REC(30425650, ref, alt, ".")
                perl = {T(span, t[1], t[2], t[4]) for t in rust}
                self.assertEqual(C.filter_dot_breakend_pairs(perl, rust), (rust, perl))
                self.assertEqual(C.check_skip_routes([record], perl, rust), [])
                perl_open, rust_open, masks = ADJUSTED(perl, rust, [record])
                self.assertEqual((perl_open, rust_open, len(masks.dot_pair_perl), len(masks.dot_pair_rust)), (set(), set(), 2, 2))
        wide_end = {T("21:30425649-30425651", t[1], t[2], t[4]) for t in rust}
        self.assertEqual(C.filter_dot_breakend_pairs(wide_end, rust), (set(), set()))


class SynonymMatePrefixTests(unittest.TestCase):
    def test_a_mate_prefix_the_reference_does_not_resolve_loses_the_mate_rows_on_the_own_chromosome(self) -> None:
        """`21 30425650 x N N[Chr21:30425700[ . . SVTYPE=BND`: Ensembl VEP 116.2 writes the own end
        only, `21:30425651 N. ENST00000334068 upstream_gene_variant` and `21:30425651 N.
        ENST00000390690 feature_truncation,coding_sequence_variant`; vep-rs adds the mate rows
        `21:30425651 N[Chr21:30425700[` on both transcripts, which are set aside. With
        `N[chr21:30425700[` both engines write the four rows and nothing is set aside, and a
        `N[NC_000022.11:...[` mate keeps its own-chromosome rule."""
        tx = {"21": {"ENST00000334068", "ENST00000390690"}, "22": {"ENST22"}}
        own = {T("21:30425651", "N.", "ENST00000334068", "upstream_gene_variant"),
               T("21:30425651", "N.", "ENST00000390690", "feature_truncation,coding_sequence_variant")}
        mate = {T("21:30425651", "N[Chr21:30425700[", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425651", "N[Chr21:30425700[", "ENST00000390690", "feature_truncation,coding_sequence_variant")}
        self.assertEqual(C.reference_mate_chromosome("N[Chr21:30425700["), "Chr21")
        self.assertEqual(C.reference_mate_chromosome("N[chr21:30425700["), "21")
        self.assertEqual(C.filter_bnd_synonym_mate_rows(own, own | mate, tx), mate)
        lower = {T(t[0], "N[chr21:30425700[", t[2], t[4]) for t in mate}
        self.assertEqual(C.filter_bnd_synonym_mate_rows(own | lower, own | lower, tx), set())
        synonym_own = {T("21:30425651", "N[NC_000022.11:26037668[", "ENST00000390690", "intron_variant")}
        synonym_other = {T("21:30425651", "N[NC_000022.11:26037668[", "ENST22", "intron_variant")}
        self.assertEqual(C.filter_bnd_synonym_mate_rows(own, own | synonym_own | synonym_other, tx), synonym_other)


class BufferCoordinateOverflowTests(unittest.TestCase):
    def test_a_span_end_or_own_position_past_the_buffers_32_bit_bound_is_counted(self) -> None:
        """Under --max_sv_size -1, Ensembl VEP 116.2 writes one row, `21:30425651-2147483647
        deletion - intergenic_variant`, for `21 30425650 x N <DEL> . . SVTYPE=DEL;END=2147483647`,
        as for `SVLEN=2147483647` (`21:30425651-2177909297`) and `SVLEN=-2147483648`
        (`21:30425651-2177909298`), where vep-rs writes every transcript from the record to the
        chromosome's end; `21 2147483648 x N N[21:30425700[ . . SVTYPE=BND` is `21:2147483649 N.`
        and `21:2147483649 N[21:30425700[`, both intergenic, as POS 2147483647 is at 2147483648.
        The counter counts the records and the reference's rows at their Locations; under the
        default cap the span records are skipped and not counted, and an end at 2147478647 is
        below the bound."""
        records = [REC(30425650, "N", "<DEL>", "SVTYPE=DEL;END=2147483647", "end"),
                   REC(30425650, "N", "<DEL>", "SVTYPE=DEL;SVLEN=2147483647", "svlen"),
                   REC(30425650, "N", "<DEL>", "SVTYPE=DEL;SVLEN=-2147483648", "negative"),
                   REC(2147483648, "N", "N[21:30425700[", "SVTYPE=BND", "pos"),
                   REC(2147483647, "N", "N[21:30425700[", "SVTYPE=BND", "pos_minus_one"),
                   REC(30425650, "N", "<DEL>", "SVTYPE=DEL;END=2147478647", "below")]
        perl = {INTERGENIC("21:30425651-2147483647", "deletion"), INTERGENIC("21:30425651-2177909297", "deletion"),
                INTERGENIC("21:30425651-2177909298", "deletion"),
                INTERGENIC("21:2147483649", "N."), INTERGENIC("21:2147483649", "N[21:30425700["),
                INTERGENIC("21:2147483648", "N."), INTERGENIC("21:2147483648", "N[21:30425700["),
                INTERGENIC("21:30425651-2147478647", "deletion")}
        self.assertEqual(C.count_bnd_mate_over_32_bit(records, perl, -1), (5, 7))
        self.assertEqual(C.count_bnd_mate_over_32_bit(records, perl), (2, 4))
        masks = C.apply_masks(perl, set(), "116.2", {}, {}, {}, "vep-rs", records, -1)
        self.assertEqual((masks.mate_overflow_records, masks.mate_overflow_perl_tuples), (5, 7))


class TandemRepeatReadingTests(unittest.TestCase):
    """The <CNV:TR> readings of Ensembl VEP 116.2 with a FASTA against vep-rs's symbolic row:
    each test is one record, its reference rows and its vep-rs rows, and every one-sided
    tuple is set aside under `cnv_tr_expansion_divergences` (the literal-versus-symbolic
    pairing), `cnv_tr_empty_expansion_dropped` (the empty or unchanged expansion) or, for a
    fractional SVLEN, `non_integer_end_annotated`."""

    RUN_A = "21:30425588-30425596"
    RUN_B = "21:30425756-30425763"
    RUN_C = "21:30430410-30430417"

    def rust_a(self, csq: str) -> set:
        return {T(self.RUN_A, "tandem_repeat", "ENST00000334068", "upstream_gene_variant"),
                T(self.RUN_A, "tandem_repeat", "ENST00000390690", csq)}

    def test_an_unchanged_expansion_over_three_repeat_sequences_is_empty(self) -> None:
        """`21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RN=3;RUS=AGC,AGC,AGC;RUC=1,1,1`
        rebuilds the nine-base run: with a FASTA Ensembl VEP 116.2 writes nothing; vep-rs writes
        `21:30425588-30425596 tandem_repeat ENST00000390690 feature_elongation,coding_sequence_variant`
        and the upstream row, both set aside."""
        record = REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RN=3;RUS=AGC,AGC,AGC;RUC=1,1,1")
        rust = self.rust_a("feature_elongation,coding_sequence_variant")
        self.assertTrue(C.cnv_tr_expansion_is_empty(record))
        self.assertEqual(ADJUSTED(set(), rust, [record])[:2], (set(), set()))

    def test_an_unchanged_expansion_whose_second_sequence_repeats_nothing_is_empty(self) -> None:
        """`21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RN=2;RUS=AGC,AG;RUC=3,0`:
        three `AGC` and no `AG` rebuild the run; Ensembl VEP 116.2 writes nothing with a FASTA,
        vep-rs the symbolic gain rows, set aside."""
        record = REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RN=2;RUS=AGC,AG;RUC=3,0")
        self.assertTrue(C.cnv_tr_expansion_is_empty(record))
        self.assertEqual(ADJUSTED(set(), self.rust_a("feature_elongation,coding_sequence_variant"), [record])[:2], (set(), set()))

    def test_a_unit_of_no_bases_under_ruc_is_empty(self) -> None:
        """`21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RUS=,AGC;RUC=4,4`: the one
        repeat sequence (no RN) takes the empty unit, `"" x 4`; Ensembl VEP 116.2 writes nothing
        with a FASTA, vep-rs `feature_truncation,inframe_deletion` with the upstream row, set
        aside."""
        record = REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=,AGC;RUC=4,4")
        self.assertTrue(C.cnv_tr_expansion_is_empty(record))
        self.assertEqual(ADJUSTED(set(), self.rust_a("feature_truncation,inframe_deletion"), [record])[:2], (set(), set()))

    def test_a_fullwidth_count_repeats_nothing(self) -> None:
        """`21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=４`: Perl
        numifies the fullwidth digit to 0, so Ensembl VEP 116.2 writes nothing with a FASTA;
        vep-rs reads the empty expansion as the loss of the run (`feature_truncation,inframe_deletion`),
        set aside."""
        record = REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC;RUC=４")
        self.assertEqual(C.perl_numeric_value("４"), 0)
        self.assertTrue(C.cnv_tr_expansion_is_empty(record))
        self.assertEqual(ADJUSTED(set(), self.rust_a("feature_truncation,inframe_deletion"), [record])[:2], (set(), set()))

    def _literal_pairs(self, record: C.VcfVariant, perl_allele: str, perl_csq: str, rust_csq: str, run: str = RUN_A) -> None:
        perl = {T(run, perl_allele, "ENST00000334068", "upstream_gene_variant"), T(run, perl_allele, "ENST00000390690", perl_csq)}
        rust = {T(run, "tandem_repeat", "ENST00000334068", "upstream_gene_variant"), T(run, "tandem_repeat", "ENST00000390690", rust_csq)}
        self.assertEqual(C._cnv_tr_pair_kind(run, perl_allele, run), "literal")
        perl_open, rust_open, masks = ADJUSTED(perl, rust, [record])
        self.assertEqual((perl_open, rust_open, masks.cnvtr_perl, masks.cnvtr_rust), (set(), set(), perl, rust))

    def test_a_digit_unit_is_a_literal_reading(self) -> None:
        """`21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RUS=123;RUC=4`: Ensembl VEP
        116.2 writes `21:30425588-30425596 123123123123 ENST00000390690 coding_sequence_variant`
        (and the upstream row) with a FASTA; vep-rs `tandem_repeat ...
        feature_elongation,coding_sequence_variant`; the four rows pair."""
        self._literal_pairs(REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=123;RUC=4"),
                            "123123123123", "coding_sequence_variant", "feature_elongation,coding_sequence_variant")

    def test_an_iupac_unit_is_a_literal_reading(self) -> None:
        """`21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RUS=R;RUC=10`: Ensembl VEP
        116.2 writes `21:30425588-30425596 RRRRRRRRRR ENST00000390690 frameshift_variant`; vep-rs
        `feature_elongation,coding_sequence_variant`; the rows pair."""
        self._literal_pairs(REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=R;RUC=10"),
                            "RRRRRRRRRR", "frameshift_variant", "feature_elongation,coding_sequence_variant")

    def test_a_star_unit_is_a_literal_reading(self) -> None:
        """`21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RUS=*;RUC=9`: Ensembl VEP
        116.2 writes `21:30425588-30425596 ********* ENST00000390690 coding_sequence_variant`;
        vep-rs `feature_elongation,coding_sequence_variant`; the rows pair, and the reference's
        rows at the run keep `cnv_tr_empty_expansion_dropped` from the record."""
        self._literal_pairs(REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=*;RUC=9"),
                            "*********", "coding_sequence_variant", "feature_elongation,coding_sequence_variant")

    def test_an_x_unit_is_a_literal_reading(self) -> None:
        """`21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RUS=X;RUC=10`: Ensembl VEP
        116.2 writes `21:30425588-30425596 XXXXXXXXXX ENST00000390690 frameshift_variant`; vep-rs
        `feature_elongation,coding_sequence_variant`; the rows pair."""
        self._literal_pairs(REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=X;RUC=10"),
                            "XXXXXXXXXX", "frameshift_variant", "feature_elongation,coding_sequence_variant")

    def test_a_loss_read_as_a_protein_altering_substitution_pairs(self) -> None:
        """`21 30425755 x N <CNV:TR> . . SVTYPE=CNV;END=30425763;SVLEN=8;RN=3;RUS=AG,AG,C;RUC=1,1,1`:
        `AGAGC` against the eight-base run trims to `21:30425760-30425763 C`, which Ensembl VEP
        116.2 writes as `protein_altering_variant` on ENST00000390690 (and `upstream_gene_variant`
        on ENST00000334068); vep-rs reads the loss, `21:30425756-30425763 tandem_repeat
        feature_truncation,inframe_deletion`; the rows pair."""
        record = REC(30425755, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425763;SVLEN=8;RN=3;RUS=AG,AG,C;RUC=1,1,1")
        perl = {T("21:30425760-30425763", "C", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425760-30425763", "C", "ENST00000390690", "protein_altering_variant")}
        rust = {T(self.RUN_B, "tandem_repeat", "ENST00000334068", "upstream_gene_variant"),
                T(self.RUN_B, "tandem_repeat", "ENST00000390690", "feature_truncation,inframe_deletion")}
        self.assertEqual(C._cnv_tr_pair_kind("21:30425760-30425763", "C", self.RUN_B), "literal")
        self.assertEqual(ADJUSTED(perl, rust, [record])[:2], (set(), set()))

    def test_a_loss_read_as_a_stop_gained_substitution_pairs(self) -> None:
        """`21 30430409 x N <CNV:TR> . . SVTYPE=CNV;END=30430417;SVLEN=8;RN=3;RUS=CT,AG,C;RUC=1,1,1`:
        Ensembl VEP 116.2 writes `21:30430412-30430417 AGC ENST00000334068 stop_gained` (and
        `upstream_gene_variant` on ENST00000390690); vep-rs `21:30430410-30430417 tandem_repeat
        feature_truncation,inframe_deletion`; the rows pair."""
        record = REC(30430409, "N", "<CNV:TR>", "SVTYPE=CNV;END=30430417;SVLEN=8;RN=3;RUS=CT,AG,C;RUC=1,1,1")
        perl = {T("21:30430412-30430417", "AGC", "ENST00000334068", "stop_gained"),
                T("21:30430412-30430417", "AGC", "ENST00000390690", "upstream_gene_variant")}
        rust = {T(self.RUN_C, "tandem_repeat", "ENST00000334068", "feature_truncation,inframe_deletion"),
                T(self.RUN_C, "tandem_repeat", "ENST00000390690", "upstream_gene_variant")}
        self.assertEqual(ADJUSTED(perl, rust, [record])[:2], (set(), set()))

    def test_a_dot_unit_under_rb_pairs_as_n_bases(self) -> None:
        """`21 30425755 x N <CNV:TR> . . SVTYPE=CNV;END=30425763;SVLEN=8;RUS=.;RB=5`: Ensembl VEP
        116.2 writes `21:30425756-30425763 NNNNN ENST00000390690 coding_sequence_variant`; vep-rs
        the loss `feature_truncation,inframe_deletion`; the rows pair."""
        self._literal_pairs(REC(30425755, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425763;SVLEN=8;RUS=.;RB=5"),
                            "NNNNN", "coding_sequence_variant", "feature_truncation,inframe_deletion", run=self.RUN_B)

    def test_an_svlen_numified_to_thirty_pairs_as_a_literal_reading(self) -> None:
        """`21 30425587 x N <CNV:TR> . . SVTYPE=CNV;SVLEN=30_425_949;RUS=AGC;RUC=4`: Perl reads the
        SVLEN as 30, a thirty-base run to 30425617 whose prefix the twelve-base allele shares,
        and Ensembl VEP 116.2 writes `21:30425597-30425617 AGC ENST00000390690
        protein_altering_variant`; vep-rs `21:30425588-30425617 tandem_repeat
        feature_truncation,inframe_deletion`; the rows pair."""
        record = REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;SVLEN=30_425_949;RUS=AGC;RUC=4")
        self.assertEqual(C.record_location(record), ("21", 30425588, 30425617))
        perl = {T("21:30425597-30425617", "AGC", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425597-30425617", "AGC", "ENST00000390690", "protein_altering_variant")}
        rust = {T("21:30425588-30425617", "tandem_repeat", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425588-30425617", "tandem_repeat", "ENST00000390690", "feature_truncation,inframe_deletion")}
        self.assertEqual(ADJUSTED(perl, rust, [record])[:2], (set(), set()))

    def test_a_loss_over_the_stop_codon_pairs_with_and_without_a_fasta(self) -> None:
        """`21 30430409 x N <CNV:TR> . . SVTYPE=CNV;SVLEN=4000;RUS=CT;RUC=5`: with a FASTA Ensembl
        VEP 116.2 writes `21:30430418-30434408 C ENST00000334068 stop_retained_variant,3_prime_UTR_variant`
        (with `upstream_gene_variant` on ENST00000390690 and `downstream_gene_variant` on
        ENST00000418755); without one `21:30430410-30434409 tandem_repeat ENST00000334068
        coding_sequence_variant,3_prime_UTR_variant`; vep-rs writes `21:30430410-30434409
        tandem_repeat ENST00000334068 stop_lost,feature_truncation,coding_sequence_variant,3_prime_UTR_variant`
        in both runs. Every row pairs."""
        record = REC(30430409, "N", "<CNV:TR>", "SVTYPE=CNV;SVLEN=4000;RUS=CT;RUC=5")
        run = "21:30430410-30434409"
        rust = {T(run, "tandem_repeat", "ENST00000334068", "stop_lost,feature_truncation,coding_sequence_variant,3_prime_UTR_variant"),
                T(run, "tandem_repeat", "ENST00000390690", "upstream_gene_variant"),
                T(run, "tandem_repeat", "ENST00000418755", "downstream_gene_variant")}
        fasta = {T("21:30430418-30434408", "C", "ENST00000334068", "stop_retained_variant,3_prime_UTR_variant"),
                 T("21:30430418-30434408", "C", "ENST00000390690", "upstream_gene_variant"),
                 T("21:30430418-30434408", "C", "ENST00000418755", "downstream_gene_variant")}
        self.assertEqual(ADJUSTED(fasta, rust, [record])[:2], (set(), set()))
        symbolic = {T(run, "tandem_repeat", "ENST00000334068", "coding_sequence_variant,3_prime_UTR_variant"),
                    T(run, "tandem_repeat", "ENST00000390690", "upstream_gene_variant"),
                    T(run, "tandem_repeat", "ENST00000418755", "downstream_gene_variant")}
        self.assertEqual(ADJUSTED(symbolic, rust, [record])[:2], (set(), set()))

    def test_a_two_alt_record_whose_parts_are_the_run_and_nothing_is_empty(self) -> None:
        """`21 30425756 x N <CNV:TR>,<CNV:TR> . . SVTYPE=CNV;END=30425763;SVLEN=8;RN=2,0;RUS=AG,AG;RUC=2,2`:
        the first allele's two sequences rebuild the eight-base run and the second expands to
        nothing, so the allele string is `AGAGAGAG/AGAGAGAG/`, its trailing part dropped and
        the rest the reference run; with a FASTA Ensembl VEP 116.2 writes nothing. vep-rs
        writes the first allele's `21:30425756-30425763 tandem_repeat ENST00000390690
        coding_sequence_variant,feature_elongation` (with `ENST00000334068
        upstream_gene_variant`) and the second's loss rows, all set aside under
        `cnv_tr_empty_expansion_dropped`; the CT run
        (`21 30430409 x N <CNV:TR>,<CNV:TR> . . SVTYPE=CNV;END=30430417;SVLEN=8;RN=2,0;RUS=CT,CT;RUC=2,2`)
        is the same shape. An empty part before a non-empty one (`RN=0,2`) is the loss the
        reference keeps, and a second allele whose count the RUC list has run out of beside a
        twelve-base first allele (`RUS=AGC,AGC;RUC=4`) is a gain the reference keeps."""
        record = REC(30425755, "N", "<CNV:TR>,<CNV:TR>", "SVTYPE=CNV;END=30425763;SVLEN=8;RN=2,0;RUS=AG,AG;RUC=2,2")
        rust = {T(self.RUN_B, "tandem_repeat", "ENST00000334068", "upstream_gene_variant"),
                T(self.RUN_B, "tandem_repeat", "ENST00000390690", "coding_sequence_variant,feature_elongation"),
                T(self.RUN_B, "tandem_repeat", "ENST00000390690", "frameshift_variant,feature_truncation")}
        self.assertTrue(C.cnv_tr_expansion_is_empty(record))
        self.assertTrue(C.cnv_tr_expansion_is_empty(
            REC(30430409, "N", "<CNV:TR>,<CNV:TR>", "SVTYPE=CNV;END=30430417;SVLEN=8;RN=2,0;RUS=CT,CT;RUC=2,2")))
        self.assertTrue(C.cnv_tr_expansion_is_empty(REC(30425587, "N", "<CNV:TR>,<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RN=.,1;RUS=AGC;RUC=3")))
        self.assertFalse(C.cnv_tr_expansion_is_empty(
            REC(30425755, "N", "<CNV:TR>,<CNV:TR>", "SVTYPE=CNV;END=30425763;SVLEN=8;RN=0,2;RUS=AG,AG;RUC=2,2")))
        self.assertFalse(C.cnv_tr_expansion_is_empty(REC(30425587, "N", "<CNV:TR>,<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC,AGC;RUC=4")))
        perl_open, rust_open, masks = ADJUSTED(set(), rust, [record])
        self.assertEqual((perl_open, rust_open, masks.cnvtr_empty_rust), (set(), set(), rust))

    def test_a_fractional_svlen_expanded_under_a_fasta_leaves_no_reference_row(self) -> None:
        """With a FASTA, Ensembl VEP 116.2 skips `21 30425587 x N <CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9.5;RUS=AGC;RUC=4`
        (`Alleles look like an insertion (-/AGC) but coordinates are not start = end + 1
        (START=30425597, END=30425596.5)`, Parser.pm:564-571), as it skips the AG run
        (`21 30425756 x N <CNV:TR> . . SVTYPE=CNV;END=30425763;SVLEN=8.5;RUS=AG;RUC=5`), the CT run
        (`SVLEN=8.5;RUS=CT;RUC=5` at 30430409) and `21 30425587 x N <CNV:TR> . . SVTYPE=CNV;SVLEN=NaN;RUS=AGC;RUC=4`;
        vep-rs reads 9.5 as 9 and NaN as 0 and writes `21:30425588-30425596 tandem_repeat
        ENST00000390690 coding_sequence_variant,feature_elongation` with `ENST00000334068
        upstream_gene_variant` (the pair `21:30425587-30425588 tandem_repeat` for NaN), set aside
        under `non_integer_end_annotated` when the reference wrote nothing beside the run. Without a FASTA the
        reference writes `21:30425588-30425596.5 tandem_repeat` and the text pairing holds; a
        `<DEL>` with `SVLEN=300.7` and no reference row is not the shape."""
        for text in ("9.5", "300.7", "-8.5", "NaN", "-Inf", "Infinity", "1e400"):
            self.assertTrue(C.perl_number_is_not_an_integer(text), text)
        for text in ("9", "9.0", "3e2", "+30425949", "abc", ".", "", "４"):
            self.assertFalse(C.perl_number_is_not_an_integer(text), text)
        record = REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9.5;RUS=AGC;RUC=4")
        rust = self.rust_a("coding_sequence_variant,feature_elongation")
        self.assertEqual(C.filter_non_integer_end_rows(set(), rust, [record]), (set(), rust))
        perl_open, rust_open, masks = ADJUSTED(set(), rust, [record])
        self.assertEqual((perl_open, rust_open, masks.non_integer_end_rust), (set(), set(), rust))
        for pos, info, run in ((30425755, "SVTYPE=CNV;END=30425763;SVLEN=8.5;RUS=AG;RUC=5", self.RUN_B),
                               (30430409, "SVTYPE=CNV;END=30430417;SVLEN=8.5;RUS=CT;RUC=5", self.RUN_C)):
            with self.subTest(info=info):
                rows = {T(run, "tandem_repeat", "ENST00000390690", "coding_sequence_variant,feature_elongation")}
                self.assertEqual(C.filter_non_integer_end_rows(set(), rows, [REC(pos, "N", "<CNV:TR>", info)]), (set(), rows))
        symbolic = {T("21:30425588-30425596.5", "tandem_repeat", "ENST00000390690", "feature_elongation,coding_sequence_variant")}
        self.assertEqual(C.filter_non_integer_end_rows(symbolic, rust, [record]), (symbolic, rust))
        nan = REC(30425587, "N", "<CNV:TR>", "SVTYPE=CNV;SVLEN=NaN;RUS=AGC;RUC=4")
        pair = {T("21:30425587-30425588", "tandem_repeat", "ENST00000390690", "coding_sequence_variant,feature_elongation"),
                T("21:30425587-30425588", "tandem_repeat", "ENST00000334068", "upstream_gene_variant")}
        self.assertIsNone(C.perl_skip_route(nan))
        self.assertEqual(C.filter_non_integer_end_rows(set(), pair, [nan]), (set(), pair))
        self.assertEqual(ADJUSTED(set(), pair, [nan])[:2], (set(), set()))
        deletion = REC(30425650, "N", "<DEL>", "SVTYPE=DEL;SVLEN=300.7")
        deletion_rows = {T("21:30425651-30425950", "deletion", "ENST00000390690", "feature_truncation,coding_sequence_variant,5_prime_UTR_variant")}
        self.assertEqual(C.filter_non_integer_end_rows(set(), deletion_rows, [deletion]), (set(), set()))

    def test_the_collapsed_symbolic_row_of_a_multi_alt_record_pairs_every_vep_rs_allele(self) -> None:
        """`21 30425587 x N <CNV:TR>,<CNV:TR>,<CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RN=1,1,1;RUS=AGC,AGC,AGC;RUC=4,2,3`
        without a FASTA: Ensembl VEP 116.2 writes the symbolic `21:30425588-30425596 tandem_repeat
        ENST00000390690 feature_elongation,coding_sequence_variant` for each of the three alleles,
        one tuple; vep-rs writes the gain for the first and third and the loss
        `feature_truncation,inframe_deletion` for the second, so the loss pairs with the tuple
        both engines wrote and only the vep-rs side is set aside."""
        record = REC(30425587, "N", "<CNV:TR>,<CNV:TR>,<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RN=1,1,1;RUS=AGC,AGC,AGC;RUC=4,2,3")
        perl = self.rust_a("feature_elongation,coding_sequence_variant")
        loss = T(self.RUN_A, "tandem_repeat", "ENST00000390690", "feature_truncation,inframe_deletion")
        rust = perl | {loss}
        self.assertEqual(C.filter_cnv_tr_expansion_divergences(perl, rust), ({loss}, set()))
        self.assertEqual(C.count_cnv_tr_swap_population(perl, rust), (1, 0))
        perl_open, rust_open, masks = ADJUSTED(perl, rust, [record])
        self.assertEqual((perl_open, rust_open, masks.cnvtr_rust, masks.cnvtr_perl), (set(), set(), {loss}, set()))
        two = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", "SVTYPE=CNV;END=30425596;SVLEN=9;RUS=AGC,AGC;RUC=4,2")
        self.assertEqual(ADJUSTED(perl, rust, [two])[:2], (set(), set()))


class RnZeroAlleleDroppedTests(unittest.TestCase):
    """With a FASTA, Ensembl VEP 116.2 drops one allele of a two-ALT <CNV:TR> whose other allele
    has `RN=0` (Parser/VCF.pm:411, 428; Parser.pm:881, 978); vep-rs writes both."""

    RUN = "21:30425588-30425596"
    INFO = "SVTYPE=CNV;END=30425596;SVLEN=9;"

    def rust(self) -> set:
        return {T(self.RUN, "tandem_repeat", "ENST00000334068", "upstream_gene_variant"),
                T(self.RUN, "tandem_repeat", "ENST00000390690", "feature_truncation,inframe_deletion"),
                T(self.RUN, "tandem_repeat", "ENST00000390690", "feature_elongation,coding_sequence_variant")}

    def test_a_leading_rn_zero_drops_the_real_allele(self) -> None:
        """`21 30425587 x N <CNV:TR>,<CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=4`:
        Ensembl VEP 116.2 keeps the empty allele as `21:30425588-30425596 - ENST00000390690
        inframe_deletion` (and `- ENST00000334068 upstream_gene_variant`) and drops the
        twelve-base gain; vep-rs writes the loss of the run for the first allele
        (`feature_truncation,inframe_deletion`) and the gain for the second
        (`feature_elongation,coding_sequence_variant`). The loss pairs under
        `cnv_tr_expansion_divergences`; the gain is the dropped allele's row, set aside."""
        record = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + "RN=0,1;RUS=AGC;RUC=4")
        perl = {T(self.RUN, "-", "ENST00000334068", "upstream_gene_variant"), T(self.RUN, "-", "ENST00000390690", "inframe_deletion")}
        gain = T(self.RUN, "tandem_repeat", "ENST00000390690", "feature_elongation,coding_sequence_variant")
        self.assertEqual(C.cnv_tr_rn_zero_dropped_direction(record), "gain")
        self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(perl, self.rust(), [record]), {gain})
        perl_open, rust_open, masks = ADJUSTED(perl, self.rust(), [record])
        self.assertEqual((perl_open, rust_open, masks.rn_zero_rust), (set(), set(), {gain}))

    def test_a_leading_rn_zero_beside_an_unchanged_real_allele_drops_the_allele_vep_rs_reads_as_a_gain(self) -> None:
        """`21 30425587 x N <CNV:TR>,<CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=3`:
        the real allele rebuilds the nine-base run (`AGCAGCAGC/-` is the allele string), so
        Ensembl VEP 116.2 keeps the empty allele, `21:30425588-30425596 - ENST00000390690
        inframe_deletion` and `- ENST00000334068 upstream_gene_variant`, and drops the unchanged
        one; vep-rs writes the loss of the run for the first allele (`tandem_repeat ENST00000390690
        feature_truncation,inframe_deletion`, pairing with the `-` rows under
        `cnv_tr_expansion_divergences`) and reads the unchanged second allele as a gain
        (`tandem_repeat ENST00000390690 feature_elongation,coding_sequence_variant`, the upstream
        row collapsing with the first allele's), whose row is the dropped allele's, set aside."""
        record = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + "RN=0,1;RUS=AGC;RUC=3")
        perl = {T(self.RUN, "-", "ENST00000334068", "upstream_gene_variant"), T(self.RUN, "-", "ENST00000390690", "inframe_deletion")}
        gain = T(self.RUN, "tandem_repeat", "ENST00000390690", "feature_elongation,coding_sequence_variant")
        self.assertFalse(C.cnv_tr_expansion_is_empty(record))
        self.assertEqual(C.cnv_tr_rn_zero_dropped_direction(record), "gain")
        self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(perl, self.rust(), [record]), {gain})
        perl_open, rust_open, masks = ADJUSTED(perl, self.rust(), [record])
        self.assertEqual((perl_open, rust_open, masks.rn_zero_rust), (set(), set(), {gain}))

    def test_a_trailing_rn_zero_drops_the_rn_zero_allele(self) -> None:
        """`21 30425587 x N <CNV:TR>,<CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RN=1,0;RUS=AGC;RUC=4`:
        Ensembl VEP 116.2 writes the gain only, `21:30425596-30425597 AGC ENST00000390690
        inframe_insertion` (and `AGC ENST00000334068 upstream_gene_variant`); vep-rs writes the
        gain for the first allele and the loss of the run for the second. The gain pairs under
        `cnv_tr_expansion_divergences`; the loss is the dropped allele's row, set aside.
        `RN=2,0;RUS=AGC,AGC;RUC=2,2` is the same shape."""
        loss = T(self.RUN, "tandem_repeat", "ENST00000390690", "feature_truncation,inframe_deletion")
        perl = {T("21:30425596-30425597", "AGC", "ENST00000334068", "upstream_gene_variant"),
                T("21:30425596-30425597", "AGC", "ENST00000390690", "inframe_insertion")}
        for info in ("RN=1,0;RUS=AGC;RUC=4", "RN=2,0;RUS=AGC,AGC;RUC=2,2"):
            with self.subTest(info=info):
                record = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + info)
                self.assertEqual(C.cnv_tr_rn_zero_dropped_direction(record), "loss")
                perl_open, rust_open, masks = ADJUSTED(perl, self.rust(), [record])
                self.assertEqual((perl_open, rust_open, masks.rn_zero_rust), (set(), set(), {loss}))

    def test_a_vep_rs_run_with_a_fasta_writes_the_dropped_allele_literally(self) -> None:
        """`21 30425587 x N <CNV:TR>,<CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RN=0,1;RUS=AGC;RUC=4`
        with a FASTA on both sides: vep-rs writes Ensembl VEP 116.2's two `-` rows for the
        `RN=0` allele and, for the dropped real allele, `21:30425588-30425596 AGCAGCAGCAGC
        ENST00000390690 protein_altering_variant` and `AGCAGCAGCAGC ENST00000334068
        upstream_gene_variant`, the whole allele over the run (the multi-allelic string is not
        trimmed). The two literal rows are the dropped allele's, set aside; a literal row under
        another allele, or at another location, is not."""
        record = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + "RN=0,1;RUS=AGC;RUC=4")
        perl = {T(self.RUN, "-", "ENST00000334068", "upstream_gene_variant"), T(self.RUN, "-", "ENST00000390690", "inframe_deletion")}
        literal = {T(self.RUN, "AGCAGCAGCAGC", "ENST00000334068", "upstream_gene_variant"),
                   T(self.RUN, "AGCAGCAGCAGC", "ENST00000390690", "protein_altering_variant")}
        rust = perl | literal
        self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(perl, rust, [record]), literal)
        perl_open, rust_open, masks = ADJUSTED(perl, rust, [record])
        self.assertEqual((perl_open, rust_open, masks.rn_zero_rust), (set(), set(), literal))
        other = T(self.RUN, "AGCAGCAGCAGCAGC", "ENST00000390690", "protein_altering_variant")
        elsewhere = T("21:30425596-30425597", "AGCAGCAGCAGC", "ENST00000390690", "inframe_insertion")
        self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(perl, rust | {other, elsewhere}, [record]), literal)
        # The unchanged real allele behind a leading `RN=0` (`RUC=3`) gets no row from a run
        # with a FASTA; a row under the run's own sequence would be the allele's literal.
        unchanged = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + "RN=0,1;RUS=AGC;RUC=3")
        self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(perl, perl, [unchanged]), set())
        run_itself = T(self.RUN, "AGCAGCAGC", "ENST00000390690", "protein_altering_variant")
        self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(perl, perl | {run_itself}, [unchanged]), {run_itself})

    def test_the_dropped_allele_is_expanded_as_the_parser_expands_it(self) -> None:
        """The literal is the real allele's repeat sequences, each unit repeated by its count:
        RB over the unit length (`RB=12` is four `AGC`), a count truncated as Perl's `x`
        truncates it (`RUC=4.9` is four), two repeat sequences joined in order (`RN=0,2`),
        `.` read as `N`, the string upper-cased, and a loss (`RUC=2`, `AGCAGC`: with a FASTA
        vep-rs writes `21:30425588-30425596 AGCAGC ENST00000390690 inframe_deletion` beside
        the `-` rows, Ensembl VEP 116.2 the `-` rows alone). Behind a trailing `RN=0` the
        dropped allele is the `RN=0` one and only its symbolic `feature_truncation` row is set
        aside, and a record outside the class sets nothing aside."""
        perl = {T(self.RUN, "-", "ENST00000334068", "upstream_gene_variant"), T(self.RUN, "-", "ENST00000390690", "inframe_deletion")}
        cases = (
            ("RN=0,1;RUS=AGC;RB=12", "AGCAGCAGCAGC"),
            ("RN=0,1;RUS=AGC;RUC=4.9", "AGCAGCAGCAGC"),
            ("RN=0,2;RUS=AGC,AG;RUC=3,2", "AGCAGCAGCAGAG"),
            ("RN=0,1;RUS=.;RUC=12", "NNNNNNNNNNNN"),
            ("RN=0,1;RUS=agc;RUC=4", "AGCAGCAGCAGC"),
            ("RN=0,1;RUS=AGC;RUC=2", "AGCAGC"),
        )
        for info, literal in cases:
            with self.subTest(info=info):
                record = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + info)
                row = T(self.RUN, literal, "ENST00000390690", "protein_altering_variant")
                decoy = T(self.RUN, literal + "AGC", "ENST00000390690", "protein_altering_variant")
                self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(perl, perl | {row, decoy}, [record]), {row})
        gain_rows = {T("21:30425596-30425597", "AGC", "ENST00000334068", "upstream_gene_variant"),
                     T("21:30425596-30425597", "AGC", "ENST00000390690", "inframe_insertion")}
        loss = T(self.RUN, "tandem_repeat", "ENST00000390690", "feature_truncation,inframe_deletion")
        row = T(self.RUN, "AGCAGCAGCAGC", "ENST00000390690", "protein_altering_variant")
        for info in ("RN=1,0;RUS=AGC;RUC=4", "RN=2,0;RUS=AGC,AGC;RUC=2,2"):
            with self.subTest(info=info):
                record = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + info)
                self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(gain_rows, gain_rows | {loss, row}, [record]), {loss})
        for info in ("RN=1,1;RUS=AGC,AGC;RUC=4,2", "RN=0,0;RUS=AGC;RUC=4", "RN=.,1;RUS=AGC;RUC=4"):
            with self.subTest(info=info):
                record = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + info)
                self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(perl, perl | {row}, [record]), set())

    def test_shapes_outside_the_class_are_not_judged(self) -> None:
        """A real allele that is a loss (`RN=0,1;RUS=AGC;RUC=2`) shares the direction term with
        the `RN=0` allele; an unchanged real allele beside a trailing empty part
        (`21 30425587 x N <CNV:TR>,<CNV:TR> . . SVTYPE=CNV;END=30425596;SVLEN=9;RN=1,0;RUS=AGC;RUC=3`,
        the allele string `AGCAGCAGC/AGCAGCAGC/`) has no reference row at all with a FASTA and both
        vep-rs rows (`feature_elongation,coding_sequence_variant` for the first allele,
        `feature_truncation,inframe_deletion` for the second) are `cnv_tr_empty_expansion_dropped`;
        `RN=.,1` leaves no empty part (`:410`); a single-ALT record, an RN with two empty parts, and
        a run without a FASTA (symbolic reference rows) are not the shape."""
        for info in ("RN=0,1;RUS=AGC;RUC=2", "RN=1,0;RUS=AGC;RUC=3", "RN=.,1;RUS=AGC;RUC=4", "RN=0,0;RUS=AGC;RUC=4",
                     "RN=1,1;RUS=AGC,AGC;RUC=4,2", "RUS=AGC;RUC=4"):
            with self.subTest(info=info):
                self.assertIsNone(C.cnv_tr_rn_zero_dropped_direction(REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + info)))
        trailing = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + "RN=1,0;RUS=AGC;RUC=3")
        self.assertTrue(C.cnv_tr_expansion_is_empty(trailing))
        perl_open, rust_open, masks = ADJUSTED(set(), self.rust(), [trailing])
        self.assertEqual((perl_open, rust_open, masks.cnvtr_empty_rust, masks.rn_zero_rust), (set(), set(), self.rust(), set()))
        self.assertIsNone(C.cnv_tr_rn_zero_dropped_direction(REC(30425587, "N", "<CNV:TR>", self.INFO + "RN=0;RUS=AGC;RUC=4")))
        record = REC(30425587, "N", "<CNV:TR>,<CNV:TR>", self.INFO + "RN=0,1;RUS=AGC;RUC=4")
        symbolic = {T(self.RUN, "tandem_repeat", "ENST00000390690", "feature_elongation,coding_sequence_variant")}
        self.assertEqual(C.filter_cnv_tr_rn_zero_allele_dropped(symbolic, self.rust(), [record]), set())


if __name__ == "__main__":
    unittest.main(verbosity=2)

