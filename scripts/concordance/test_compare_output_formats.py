#!/usr/bin/env python3
"""Unit tests for the parsers + key projection in compare_output_formats.

Run with: ``python3 -m pytest -q scripts/concordance/test_compare_output_formats.py``
Or: ``python3 scripts/concordance/test_compare_output_formats.py``

These lock down the behaviors the output-format-parity verdict depends on:
  - the three parsers (vcf / tab / json) each emit the expected 5-tuple set,
  - ``norm_consequence`` normalizes the ``&`` (VCF) and ``,`` (tab/json)
    within-entry term separators to the same canonical string,
  - ``project_key`` collapses the format-specific coordinate in ``calls`` mode
    (the load-bearing claim behind "vcf == tab == json IDENTICAL") while
    ``strict`` mode keeps it, and ``calls`` still discriminates by contig.

The fixtures follow the engine output conventions of the vep-rs formatters in
``crates/vep-io/src/output/``: VCF joins multi-value CSQ fields with ``&`` and
joins CSQ entries with ``,`` (``vcf_output.rs``); VEP tab joins consequence
terms with ``,`` and puts ``IMPACT=`` first in the ``;``-delimited Extra column
(``tab.rs``); JSON uses ``seq_region_name`` / ``variant_allele`` /
``consequence_terms`` (``json.rs``).
"""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from compare_output_formats import (
    _impact_from_extra,
    norm_consequence,
    norm_contig,
    parse_json,
    parse_tab,
    parse_vcf,
    project_key,
)


class _FixtureMixin:
    """Provides a per-test temp dir and a writer that returns the path."""

    def setUp(self) -> None:  # noqa: N802 (unittest API)
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self) -> None:  # noqa: N802 (unittest API)
        self._tmp.cleanup()

    def write(self, name: str, content: str) -> Path:
        p = self.tmp / name
        p.write_text(content, encoding="utf-8")
        return p

    def write_gz(self, name: str, content: str) -> Path:
        p = self.tmp / name
        with gzip.open(p, "wt", encoding="utf-8") as fh:
            fh.write(content)
        return p


class NormContigTests(unittest.TestCase):
    def test_strip_chr_prefix(self) -> None:
        self.assertEqual(norm_contig("chr21"), "21")
        self.assertEqual(norm_contig("21"), "21")

    def test_chrM_maps_to_MT(self) -> None:
        self.assertEqual(norm_contig("chrM"), "MT")
        self.assertEqual(norm_contig("M"), "MT")

    def test_MT_passthrough(self) -> None:
        self.assertEqual(norm_contig("MT"), "MT")


class NormConsequenceTests(unittest.TestCase):
    def test_ampersand_separator(self) -> None:
        """VCF joins terms with '&'."""
        self.assertEqual(
            norm_consequence("missense_variant&splice_region_variant"),
            "missense_variant,splice_region_variant",
        )

    def test_comma_separator(self) -> None:
        """VEP tab/json join terms with ','."""
        self.assertEqual(
            norm_consequence("splice_region_variant,missense_variant"),
            "missense_variant,splice_region_variant",
        )

    def test_order_canonicalized(self) -> None:
        """Same set, different order -> same canonical string."""
        self.assertEqual(
            norm_consequence("b_variant&a_variant"),
            norm_consequence("a_variant,b_variant"),
        )

    def test_different_sets_do_not_collapse(self) -> None:
        """A strict subset must NOT canonicalize to the same string."""
        self.assertNotEqual(
            norm_consequence("missense_variant"),
            norm_consequence("missense_variant&splice_region_variant"),
        )

    def test_empty_terms_dropped(self) -> None:
        self.assertEqual(norm_consequence("a&&b"), "a,b")


class ImpactFromExtraTests(unittest.TestCase):
    def test_first_key(self) -> None:
        self.assertEqual(_impact_from_extra("IMPACT=MODIFIER;STRAND=1"), "MODIFIER")

    def test_not_first_key(self) -> None:
        self.assertEqual(_impact_from_extra("SOMEKEY=x;IMPACT=HIGH"), "HIGH")

    def test_absent_returns_empty(self) -> None:
        self.assertEqual(_impact_from_extra("STRAND=1;DISTANCE=2164"), "")

    def test_substring_key_not_matched(self) -> None:
        """The (?:^|;) anchor must reject a substring like XIMPACT=."""
        self.assertEqual(_impact_from_extra("XIMPACT=LOW"), "")


class ParseVcfTests(_FixtureMixin, unittest.TestCase):
    def test_multi_allelic_with_ampersand_consequence(self) -> None:
        body = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t.\tA\tT,AGT\t.\t.\t"
            "CSQ=T|missense_variant&splice_region_variant|MODERATE|S|G|"
            "Transcript|ENST1|pc,AGT|frameshift_variant|HIGH|S|G|"
            "Transcript|ENST1|pc\n"
        )
        got = parse_vcf(self.write("ma.vcf", body))
        self.assertEqual(
            got,
            {
                (
                    "1:100",
                    "T",
                    "ENST1",
                    "missense_variant,splice_region_variant",
                    "MODERATE",
                ),
                ("1:100", "AGT", "ENST1", "frameshift_variant", "HIGH"),
            },
        )

    def test_header_overrides_default_field_order(self) -> None:
        fmt = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature"
        body = (
            f'##INFO=<ID=CSQ,Number=.,Type=String,Description="x Format: {fmt}">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t.\tA\tT\t.\t.\t"
            "CSQ=T|stop_gained|HIGH|S|G|Transcript|ENSTH\n"
        )
        got = parse_vcf(self.write("hdr.vcf", body))
        self.assertEqual(got, {("1:100", "T", "ENSTH", "stop_gained", "HIGH")})

    def test_default_field_order_and_chr_strip(self) -> None:
        body = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr21\t500\t.\tA\tG\t.\t.\t"
            "CSQ=G|synonymous_variant|LOW|S|Gn|Transcript|ENSTd|pc\n"
        )
        got = parse_vcf(self.write("def.vcf", body))
        self.assertEqual(got, {("21:500", "G", "ENSTd", "synonymous_variant", "LOW")})

    def test_record_without_csq_skipped(self) -> None:
        body = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t.\tA\tT\t.\t.\tAC=1;AN=2\n"
        )
        self.assertEqual(parse_vcf(self.write("nocsq.vcf", body)), set())

    def test_short_record_skipped(self) -> None:
        body = "#CHROM\tPOS\tID\tREF\tALT\n1\t100\t.\tA\tT\n"
        self.assertEqual(parse_vcf(self.write("short.vcf", body)), set())

    def test_empty_file(self) -> None:
        self.assertEqual(parse_vcf(self.write("empty.vcf", "")), set())

    def test_header_only(self) -> None:
        body = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        self.assertEqual(parse_vcf(self.write("hdronly.vcf", body)), set())

    def test_gzip_input(self) -> None:
        body = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t.\tA\tT\t.\t.\t"
            "CSQ=T|missense_variant|MODERATE|S|G|Transcript|ENSTz|pc\n"
        )
        got = parse_vcf(self.write_gz("z.vcf.gz", body))
        self.assertEqual(got, {("1:100", "T", "ENSTz", "missense_variant", "MODERATE")})


class ParseTabTests(_FixtureMixin, unittest.TestCase):
    HEADER = "#Uploaded_variation\tLocation\tAllele\tFeature\tConsequence\tExtra\n"

    def test_basic_row(self) -> None:
        body = self.HEADER + (
            "rs1\tchr1:100\tT\tENST1\tmissense_variant\tIMPACT=MODERATE;STRAND=1\n"
        )
        got = parse_tab(self.write("b.tab", body))
        self.assertEqual(got, {("1:100", "T", "ENST1", "missense_variant", "MODERATE")})

    def test_impact_absent_in_extra(self) -> None:
        body = self.HEADER + ("rs1\t1:100\tT\tENST1\tsynonymous_variant\tSTRAND=1\n")
        got = parse_tab(self.write("noimp.tab", body))
        self.assertEqual(got, {("1:100", "T", "ENST1", "synonymous_variant", "")})

    def test_no_header_skips_all_rows(self) -> None:
        body = "rs1\tchr1:100\tT\tENST1\tmissense_variant\tIMPACT=HIGH\n"
        self.assertEqual(parse_tab(self.write("nohdr.tab", body)), set())

    def test_empty_file(self) -> None:
        self.assertEqual(parse_tab(self.write("e.tab", "")), set())


class ParseJsonTests(_FixtureMixin, unittest.TestCase):
    def test_line_form_with_variant_allele(self) -> None:
        body = (
            '{"seq_region_name":"chr1","start":100,"allele_string":"A/T",'
            '"transcript_consequences":[{"transcript_id":"ENST1",'
            '"variant_allele":"T","consequence_terms":["missense_variant"],'
            '"impact":"MODERATE"}]}\n'
        )
        got = parse_json(self.write("line.json", body))
        self.assertEqual(got, {("1:100", "T", "ENST1", "missense_variant", "MODERATE")})

    def test_array_form(self) -> None:
        body = (
            '[{"seq_region_name":"1","start":200,"allele_string":"C/G",'
            '"transcript_consequences":[{"transcript_id":"ENSTa",'
            '"variant_allele":"G","consequence_terms":["stop_gained"],'
            '"impact":"HIGH"}]}]'
        )
        got = parse_json(self.write("arr.json", body))
        self.assertEqual(got, {("1:200", "G", "ENSTa", "stop_gained", "HIGH")})

    def test_variant_allele_fallback_indel(self) -> None:
        """variant_allele absent -> last segment of REF/ALT allele_string."""
        body = (
            '{"seq_region_name":"1","start":100,"allele_string":"G/GT",'
            '"transcript_consequences":[{"transcript_id":"ENSTi",'
            '"consequence_terms":["frameshift_variant"],"impact":"HIGH"}]}\n'
        )
        got = parse_json(self.write("indel.json", body))
        self.assertEqual(got, {("1:100", "GT", "ENSTi", "frameshift_variant", "HIGH")})

    def test_variant_allele_fallback_snv(self) -> None:
        body = (
            '{"seq_region_name":"1","start":100,"allele_string":"A/T",'
            '"transcript_consequences":[{"transcript_id":"ENSTs",'
            '"consequence_terms":["missense_variant"],"impact":"MODERATE"}]}\n'
        )
        got = parse_json(self.write("snv.json", body))
        self.assertEqual(got, {("1:100", "T", "ENSTs", "missense_variant", "MODERATE")})

    def test_empty_file(self) -> None:
        self.assertEqual(parse_json(self.write("e.json", "")), set())


class ProjectKeyTests(unittest.TestCase):
    """The load-bearing parity-verdict semantics."""

    # Two tuples that differ ONLY in the coordinate (the VCF-POS-vs-VEP-span
    # representation difference): same contig, allele, feature, consequence,
    # impact.
    A = ("1:100", "T", "ENST1", "missense_variant", "MODERATE")
    B = ("1:101-102", "T", "ENST1", "missense_variant", "MODERATE")

    def test_calls_mode_collapses_coordinate(self) -> None:
        self.assertEqual(len(project_key({self.A, self.B}, "calls")), 1)

    def test_calls_mode_key_form_drops_coordinate_keeps_contig(self) -> None:
        self.assertEqual(
            project_key({self.A}, "calls"),
            {("1", "T", "ENST1", "missense_variant", "MODERATE")},
        )

    def test_strict_mode_keeps_distinct(self) -> None:
        self.assertEqual(len(project_key({self.A, self.B}, "strict")), 2)

    def test_calls_mode_still_discriminates_contig(self) -> None:
        """Dropping the coordinate must NOT drop the contig."""
        c1 = ("1:100", "T", "ENST1", "missense_variant", "MODERATE")
        c2 = ("2:100", "T", "ENST1", "missense_variant", "MODERATE")
        self.assertEqual(len(project_key({c1, c2}, "calls")), 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
