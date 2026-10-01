#!/usr/bin/env python3
"""Tests for the breakend mate-context derivation.

Every expected set is derived by hand from the transcript structure quoted in the test,
and every exemplar is a record of the structural-variant corpora whose Perl and vep-rs
rows are known.

Run:
    python3 -m pytest scripts/validation/test_breakend_mate_context.py -q
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import breakend_mate_context as B


def exons(spec: str) -> tuple[tuple[int, int], ...]:
    return tuple(tuple(int(x) for x in e.split("-")) for e in spec.split(";"))


# ENST00000270142 (chromosome 21, forward strand, protein coding): five exons, CDS
# 33032083-33040891.
APP_LIKE = B.TranscriptModel(
    chrom="21", start=33031935, end=33041244, strand=1, biotype="protein_coding",
    exons=exons("33031935-33032154;33036103-33036199;33038762-33038831;33039571-33039688;33040784-33041244"),
    cds=(33032083, 33040891), mirna_mature=(),
)
# ENST00000248935 GSTT1 (chromosome 22, reverse strand, protein coding), CDS 24376423-24384231.
GSTT1 = B.TranscriptModel(
    chrom="22", start=24376133, end=24384284, strand=-1, biotype="protein_coding",
    exons=exons("24376133-24376617;24376822-24376998;24379361-24379511;24381700-24381787;24384120-24384284"),
    cds=(24376423, 24384231), mirna_mature=(),
)
# ENST00000265107 WDR70 (chromosome 5, forward strand, protein coding), CDS 37379368-37752573.
WDR70 = B.TranscriptModel(
    chrom="5", start=37379318, end=37753435, strand=1, biotype="protein_coding",
    exons=exons("37379318-37379392;37379489-37379554;37381602-37381685;37392000-37392120;37396375-37396570;"
                "37437922-37437981;37443239-37443372;37479834-37479987;37516514-37516590;37605064-37605238;"
                "37697655-37697754;37701058-37701142;37702949-37703087;37721115-37721215;37722855-37722934;"
                "37724934-37725050;37726883-37727045;37752486-37753435"),
    cds=(37379368, 37752573), mirna_mature=(),
)
# ENST00000385128 MIR125B2 (chromosome 21, forward strand, miRNA), one exon, mature 17-38 and 54-75.
MIR = B.TranscriptModel(
    chrom="21", start=16590237, end=16590325, strand=1, biotype="miRNA",
    exons=exons("16590237-16590325"), cds=None, mirna_mature=((17, 38), (54, 75)),
)
# ENST00000687248 PTCHD1-AS (chromosome X, reverse strand, lncRNA), nine exons, no CDS.
PTCHD1_AS = B.TranscriptModel(
    chrom="X", start=22227405, end=23293174, strand=-1, biotype="lncRNA",
    exons=exons("22227405-22227667;22253939-22254065;22258013-22258136;22259755-22259842;22389663-22389720;"
                "22743229-22743338;23064038-23064121;23270051-23270176;23293014-23293174"),
    cds=None, mirna_mature=(),
)


class ModelTests(unittest.TestCase):
    def test_introns_are_the_gaps_between_sorted_exons(self) -> None:
        self.assertEqual(APP_LIKE.introns[0], (33032155, 33036102))
        self.assertEqual(len(APP_LIKE.introns), 4)
        self.assertEqual(MIR.introns, ())

    def test_genomic_cds_is_mapped_through_the_pairs_in_both_orientations(self) -> None:
        """Forward: cDNA 149 in a first exon of 220 bases starting at 33031935 is
        33032083. Reverse: cDNA 1 of a reverse exon 24384120-24384284 is 24384284, so a
        coding start at cDNA 54 is 24384231."""
        fwd = {"variation_effect_feature_cache": {"mapper": {
            "cdna_coding_start": 149, "cdna_coding_end": 613,
            "pairs": [{"from_start": 1, "from_end": 220, "to_start": 33031935, "to_end": 33032154, "ori": 1},
                      {"from_start": 221, "from_end": 317, "to_start": 33036103, "to_end": 33036199, "ori": 1},
                      {"from_start": 318, "from_end": 387, "to_start": 33038762, "to_end": 33038831, "ori": 1},
                      {"from_start": 388, "from_end": 505, "to_start": 33039571, "to_end": 33039688, "ori": 1},
                      {"from_start": 506, "from_end": 966, "to_start": 33040784, "to_end": 33041244, "ori": 1}]}}}
        self.assertEqual(B._genomic_cds(fwd), (33032083, 33040891))
        rev = {"variation_effect_feature_cache": {"mapper": {
            "cdna_coding_start": 54, "cdna_coding_end": 195,
            "pairs": [{"from_start": 1, "from_end": 165, "to_start": 24384120, "to_end": 24384284, "ori": -1},
                      {"from_start": 166, "from_end": 253, "to_start": 24381700, "to_end": 24381787, "ori": -1}]}}}
        self.assertEqual(B._genomic_cds(rev), (24381758, 24384231))

    def test_a_transcript_without_coding_bounds_has_no_cds(self) -> None:
        self.assertIsNone(B._genomic_cds({"variation_effect_feature_cache": {"mapper": {"pairs": []}}}))
        self.assertIsNone(B._genomic_cds({}))

    def test_models_are_keyed_by_chromosome_and_id(self) -> None:
        """A pseudoautosomal transcript is filed under X and Y with different coordinates;
        one key per chromosome keeps both."""
        with tempfile.TemporaryDirectory() as d:
            for chrom, start in (("X", 200000), ("Y", 150000)):
                p = Path(d) / "transcripts" / chrom
                p.mkdir(parents=True)
                (p / "1-1000000.json").write_text(json.dumps([{
                    "stable_id": "ENST_PAR", "start": start, "end": start + 999, "strand": 1, "biotype": "protein_coding",
                    "exons": [{"start": start, "end": start + 999}],
                    "variation_effect_feature_cache": {"mapper": {"pairs": [], "pair_count": 0}},
                }]))
            models = B.load_transcript_models(d)
        self.assertEqual(models[("X", "ENST_PAR")].start, 200000)
        self.assertEqual(models[("Y", "ENST_PAR")].start, 150000)

    def test_a_storable_shaped_shard_is_read(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "transcripts" / "21"
            p.mkdir(parents=True)
            (p / "1-1000000.json").write_text(json.dumps({"21": [None, {
                "stable_id": "ENST_A", "start": 10, "end": 20, "strand": -1, "biotype": "lncRNA",
                "exons": [{"start": 10, "end": 20}], "attributes": [{"code": "miRNA", "value": "3-7"}],
            }]}))
            models = B.load_transcript_models(d)
        self.assertEqual(models[("21", "ENST_A")].mirna_mature, ((3, 7),))
        self.assertIsNone(models[("21", "ENST_A")].cds)


class RegionTermTests(unittest.TestCase):
    def test_cds_utr_by_strand(self) -> None:
        self.assertEqual(B.region_terms(33032100, APP_LIKE), {"coding_sequence_variant"})
        self.assertEqual(B.region_terms(33032000, APP_LIKE), {"5_prime_UTR_variant"})
        self.assertEqual(B.region_terms(33041000, APP_LIKE), {"3_prime_UTR_variant"})
        # reverse strand: bases beyond the CDS end are the 5' UTR
        self.assertEqual(B.region_terms(24384250, GSTT1), {"5_prime_UTR_variant"})
        self.assertEqual(B.region_terms(24376300, GSTT1), {"3_prime_UTR_variant"})

    def test_intron_terms(self) -> None:
        """Inside intron 1 of APP_LIKE (33032155-33036102): the first two and last two
        bases are the splice sites and carry no intron term; the forward-strand
        polypyrimidine window is the acceptor side, 16 to 2 bases before the intron end."""
        self.assertEqual(B.region_terms(33034355, APP_LIKE), {"intron_variant"})
        self.assertEqual(B.region_terms(33032155, APP_LIKE), set())
        self.assertEqual(B.region_terms(33036102, APP_LIKE), set())
        self.assertEqual(B.region_terms(33036102 - 5, APP_LIKE), {"intron_variant", "splice_polypyrimidine_tract_variant"})
        self.assertEqual(B.region_terms(33036102 - 17, APP_LIKE), {"intron_variant"})
        # reverse strand: the acceptor is the intron START side
        self.assertEqual(B.region_terms(24376618 + 5, GSTT1), {"intron_variant", "splice_polypyrimidine_tract_variant"})

    def test_non_coding_intron_and_exon(self) -> None:
        self.assertEqual(B.region_terms(22500000, PTCHD1_AS), {"intron_variant", "non_coding_transcript_variant"})
        self.assertEqual(B.region_terms(22227500, PTCHD1_AS), {"non_coding_transcript_exon_variant"})

    def test_nmd_transcript_adds_its_term(self) -> None:
        nmd = B.TranscriptModel(chrom="21", start=100, end=1000, strand=1, biotype="nonsense_mediated_decay",
                                exons=((100, 200), (800, 1000)), cds=(150, 900), mirna_mature=())
        self.assertEqual(B.region_terms(500, nmd), {"intron_variant", "NMD_transcript_variant"})
        self.assertEqual(B.region_terms(160, nmd), {"coding_sequence_variant", "NMD_transcript_variant"})

    def test_mature_mirna_hit_stands_alone(self) -> None:
        """cDNA 17-38 of a forward single-exon miRNA starting at 16590237 is 16590253-16590274."""
        self.assertEqual(B.region_terms(16590260, MIR), {"mature_miRNA_variant"})
        self.assertEqual(B.region_terms(16590240, MIR), {"non_coding_transcript_exon_variant"})


class BreakendContextTests(unittest.TestCase):
    def test_inside_with_the_mate_inside_carries_feature_truncation(self) -> None:
        self.assertEqual(B.breakend_context(33034355, APP_LIKE, mate_inside=True), {"feature_truncation", "intron_variant"})

    def test_upstream_downstream_are_strand_aware(self) -> None:
        # 4,108 bases before the reverse-strand GSTT1 start: downstream
        self.assertEqual(B.breakend_context(24372025, GSTT1, mate_inside=False), {"downstream_gene_variant"})
        self.assertEqual(B.breakend_context(24384284 + 100, GSTT1, mate_inside=False), {"upstream_gene_variant"})
        self.assertEqual(B.breakend_context(33031935 - 100, APP_LIKE, mate_inside=False), {"upstream_gene_variant"})

    def test_beyond_five_kb_carries_only_the_truncation_term_when_the_mate_is_inside(self) -> None:
        self.assertEqual(B.breakend_context(14374334, APP_LIKE, mate_inside=True), {"feature_truncation"})
        self.assertEqual(B.breakend_context(14374334, APP_LIKE, mate_inside=False), set())

    def test_a_mature_mirna_hit_takes_no_truncation_term(self) -> None:
        self.assertEqual(B.breakend_context(16590260, MIR, mate_inside=True), {"mature_miRNA_variant"})


class MateSideDerivationTests(unittest.TestCase):
    def test_cross_chromosome_row_perl_read_at_the_local_coordinate(self) -> None:
        """`1:14374334 ]21:33034355]A` on ENST00000270142: the mate sits in intron 1, the
        local coordinate on chromosome 1 is nowhere near, so Perl writes only
        `feature_truncation` and the mate supports `feature_truncation,intron_variant`."""
        got = B.mate_side_derivation("1:14374334", "]21:33034355]A", APP_LIKE)
        self.assertEqual(got, (frozenset({"feature_truncation", "intron_variant"}), frozenset({"feature_truncation"})))

    def test_the_local_coordinate_landing_in_a_cds_exon_of_the_mate_transcript(self) -> None:
        """`21:37392051 [5:37717217[A` on WDR70: the mate is in intron 13
        (37703088-37721114), the local coordinate lands numerically in exon 4
        (37392000-37392120), inside the CDS: Perl's read is
        `feature_truncation,coding_sequence_variant`."""
        got = B.mate_side_derivation("21:37392051", "[5:37717217[A", WDR70)
        self.assertEqual(got, (frozenset({"feature_truncation", "intron_variant"}), frozenset({"feature_truncation", "coding_sequence_variant"})))

    def test_a_far_local_coordinate_leaves_only_the_truncation_term(self) -> None:
        got = B.mate_side_derivation("21:5100001", "[5:37717217[T", WDR70)
        self.assertEqual(got[1], frozenset({"feature_truncation"}))

    def test_a_far_local_coordinate_beside_an_outside_mate_defaults_to_intergenic(self) -> None:
        """`21:9500001 N[22:24372025[` on GSTT1: the mate is 4,108 bases downstream, the
        local coordinate fires nothing, and Perl writes `intergenic_variant` on the
        Transcript row."""
        got = B.mate_side_derivation("21:9500001", "N[22:24372025[", GSTT1)
        self.assertEqual(got, (frozenset({"downstream_gene_variant"}), frozenset({"intergenic_variant"})))

    def test_equal_numeric_coordinates_give_equal_sets(self) -> None:
        """`21:24372025 N[22:24372025[`: the two reads coincide, which is why Perl's row is
        right on this record."""
        mate_set, local_set = B.mate_side_derivation("21:24372025", "N[22:24372025[", GSTT1)
        self.assertEqual(mate_set, local_set)
        self.assertEqual(mate_set, frozenset({"downstream_gene_variant"}))

    def test_a_row_the_local_coordinate_selects_is_not_a_mate_side_row(self) -> None:
        self.assertIsNone(B.mate_side_derivation("21:33034000", "N[21:33034355[", APP_LIKE))

    def test_a_ranged_location_or_a_non_bracket_allele_is_not_a_mate_side_row(self) -> None:
        self.assertIsNone(B.mate_side_derivation("21:14604132-41770788", "N[21:41770788[", APP_LIKE))
        self.assertIsNone(B.mate_side_derivation("1:14374334", "deletion", APP_LIKE))

    def test_a_mate_beyond_five_kb_of_the_transcript_is_not_a_mate_side_row(self) -> None:
        self.assertIsNone(B.mate_side_derivation("1:14374334", "]21:33050000]A", APP_LIKE))

    def test_a_transcript_on_another_chromosome_than_the_mate_is_not_a_mate_side_row(self) -> None:
        self.assertIsNone(B.mate_side_derivation("1:14374334", "]22:33034355]A", APP_LIKE))

    def test_helpers(self) -> None:
        self.assertEqual(B.bracket_mate("N[chr21:41770788["), ("21", 41770788))
        self.assertEqual(B.bracket_mate("]13:123456]T"), ("13", 123456))
        self.assertIsNone(B.bracket_mate("N."))
        self.assertEqual(B.point_location("chr21:5"), ("21", 5))
        self.assertIsNone(B.point_location("21:1-2"))


if __name__ == "__main__":
    unittest.main()
