"""Tests for the golden-corpus scripts on synthetic inputs, plus checks of the committed corpora."""
from __future__ import annotations

import json
import hashlib
import gzip
import os
import random
import re
import shutil
import struct
import subprocess
import sys
import zlib
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_golden_corpus as bgc  # noqa: E402
import extract_reference as xr  # noqa: E402
import make_custom_fixtures as mcf  # noqa: E402
import prune_json_cache as pjc  # noqa: E402


def test_auto_names_follow_vep_rules():
    # The first name is release 116's, the raw line as CHR_POS_REF/ALTs; the rest are
    # release 115's minimised names.
    assert bgc.vep_auto_names("21", 100, "A", ["G"]) == ["21_100_A/G"]
    assert bgc.vep_auto_names("21", 100, "AT", ["A"]) == ["21_100_AT/A", "21_101_AT/A"]  # bi-allelic indel: raw alleles, trimmed start
    assert bgc.vep_auto_names("21", 100, "A", ["AC"]) == ["21_100_A/AC", "21_101_A/AC"]
    assert bgc.vep_auto_names("21", 100, "A", ["G", "T"]) == ["21_100_A/G/T"]  # multi-allelic SNV: no chop
    assert bgc.vep_auto_names("21", 100, "TAG", ["TA", "T"]) == ["21_100_TAG/TA/T", "21_101_TAG/TA/T"]  # multi-allelic indel: raw alleles, chopped start
    assert bgc.vep_auto_names("21", 100, "ATT", ["AT"]) == ["21_100_ATT/AT", "21_102_ATT/AT"]  # repeat: chop, then trim the shared T
    assert bgc.vep_auto_names("21", 100, "C", ["CT", "G"]) == ["21_100_C/CT/G"]  # mixed first bases: no chop
    assert bgc.vep_auto_names("21", 100, "T", ["<CN0>"]) == ["21_100_T/<CN0>", "21_100_<CN0>", "21_101_<CN0>"]
    assert bgc.vep_auto_names("21", 100, "AT", ["TG"]) == ["21_100_AT/TG"]  # equal length: raw at POS


def test_record_names_take_the_first_id_element_perl_reads_as_true():
    """The GRCh38-vcf-edge-cases corpus proves the rule: `rs1;rs2` names its record `rs1`,
    while `.`, `;rs2`, `0;rs2`, an empty ID column and a multi-allelic `0` leave Ensembl VEP
    to name the record from its line. An empty or `0` element also offers the name built from
    the parsed alleles in the case written, the one Ensembl VEP prints for the
    GRCh38-vcf-check-ref corpus's `21 30319403 0 tca t` (`21_30319404_ca/-`)."""
    assert bgc.record_names("21\t30289183\trs1;rs2\tC\tG\t.\t.\t.") == ["rs1"]
    assert bgc.record_names("21\t30289183\t.\tC\tG,T\t.\t.\t.") == ["21_30289183_C/G/T"]
    assert bgc.record_names("21\t30289183\t;rs2\tC\tG\t.\t.\t.") == ["21_30289183_C/G"]
    assert bgc.record_names("21\t30289183\t0;rs2\tC\tT\t.\t.\t.") == ["21_30289183_C/T"]
    assert bgc.record_names("21\t30396602\t\tTTGA\tT\t.\t.\t.") == [
        "21_30396602_TTGA/T",
        "21_30396603_TTGA/T",
        "21_30396603_TGA/-",
    ]
    assert bgc.record_names("21\t30396602\t0\tTTGA\tT,TTGATC\t.\t.\t.") == [
        "21_30396602_TTGA/T/TTGATC",
        "21_30396603_TTGA/T/TTGATC",
        "21_30396603_TGA/-/TGATC",
    ]
    assert bgc.record_names("21\t30319403\t0\ttca\tt\t.\t.\t.") == [
        "21_30319403_tca/t",
        "21_30319404_tca/t",
        "21_30319404_ca/-",
    ]
    # A shared trailing base is trimmed from the alleles but moves no start; a `.` element
    # never takes the parsed name, and a symbolic allele has none.
    assert bgc.record_names("21\t100\t\tattg\tatg\t.\t.\t.") == ["21_100_attg/atg", "21_102_attg/atg", "21_102_t/-"]
    assert bgc.record_names("21\t30319403\t.\ttca\tt\t.\t.\t.") == ["21_30319403_tca/t", "21_30319404_tca/t"]
    assert bgc.record_names("21\t100\t0\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=200") == bgc.vep_auto_names("21", 100, "N", ["<DEL>"])


def test_resolve_records_finds_a_record_under_every_name_it_can_carry(tmp_path: Path):
    """`select` resolves a wanted Uploaded_variation to its input record by the same names
    `record_names` lists, an ID-bearing record by its auto-names as well."""
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t30319403\t0\ttca\tt\t.\t.\t.\n"
        "21\t30289183\trs1;rs2\tC\tG\t.\t.\t.\n"
        "21\t30289184\trs3\tC\tG\t.\t.\t.\n",
        encoding="utf-8",
    )
    wanted = {"21_30319404_ca/-", "rs1", "21_30289184_C/G"}
    resolved = bgc.resolve_records([("s", vcf)], {"s": wanted})
    assert {name for name, _chrom in resolved.found["s"]} == wanted
    assert resolved.found["s"][("21_30319404_ca/-", "21")] == ["21\t30319403\t0\ttca\tt\t.\t.\t."]


def test_variant_class_inference():
    assert bgc.variant_class("21:100", "G") == "snv"
    assert bgc.variant_class("21:101", "-") == "deletion"
    assert bgc.variant_class("21:100-101", "C") == "insertion"
    assert bgc.variant_class("21:100-105", "-") == "deletion"
    assert bgc.variant_class("21:100-105", "deletion") == "symbolic_deletion"
    assert bgc.variant_class("21:100-105", "Alu_insertion") == "symbolic_Alu_insertion"
    assert bgc.variant_class("21:100-101", "TG") == "insertion"
    assert bgc.variant_class("21:100", "TG") == "mnv_or_complex"
    assert bgc.variant_class("21:100", "A[21:200[") == "breakend"
    assert bgc.variant_class("21:100", "]MT:200]G") == "breakend"  # the mate's contig name is not inserted sequence
    assert bgc.variant_class("21:100", "[21:200[TTGA") == "breakend_insertion"
    assert bgc.variant_class("21:100", ".G") == "single_breakend"


def test_a_two_base_substitution_is_classed_from_its_input_record(tmp_path: Path, capsys):
    """A two-base substitution and a two-base insertion print the same Location and Allele, so
    such a row is classed with its input record: the substitution is `mnv_or_complex`, a
    stratum and a class exemplar of its own, counted so in the manifest's `variant_classes`
    and `family_strata`, in the rows the record-class fill adds and by `check`; the insertion
    stays `insertion`, as does a row read with no input record."""
    ref = tmp_path / "ref.txt"
    ref.write_text(
        _row("mnv1", "21:100-101", "TG", "ENST1", "coding_sequence_variant")
        + _row("ins1", "21:200-201", "TG", "ENST1", "coding_sequence_variant")
        + _row("21_300_AC/TG", "21:300-301", "TG", "ENST1", "coding_sequence_variant"),
        encoding="utf-8",
    )
    vcf = tmp_path / "in.vcf"
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    insertion = "21\t200\tins1\tA\tATG\t.\tPASS\t.\n"
    vcf.write_text(header + "21\t100\tmnv1\tAC\tTG\t.\tPASS\t.\n" + insertion + "21\t300\t.\tAC\tTG\t.\tPASS\t.\n", encoding="utf-8")

    # One exemplar per combination: the second class enters as a stratum of its own.
    out = tmp_path / "corpus"
    assert _select(ref, vcf, out, "--covers", "symbol", "--family-columns", "IMPACT=") == 0
    m = json.loads((out / "manifest.json").read_text())
    assert m["variant_classes"] == {"insertion": 1, "mnv_or_complex": 1}
    assert m["family_strata"] == {"IMPACT=MODIFIER;insertion": 1, "IMPACT=MODIFIER;mnv_or_complex": 1}
    assert {r["vcf"].split("\t")[3] for r in m["records"]} == {"A", "AC"}
    check = ["check", "--family-columns", "IMPACT=", "--corpus"]
    assert bgc.main([*check, str(out), "--suite", f"s01={ref}={vcf}"]) == 0

    # A corpus of the insertion alone lacks the substitution's stratum.
    ins_ref, ins_vcf, lacking = tmp_path / "ins.txt", tmp_path / "ins.vcf", tmp_path / "lacking"
    ins_ref.write_text(_row("ins1", "21:200-201", "TG", "ENST1", "coding_sequence_variant"), encoding="utf-8")
    ins_vcf.write_text(header + insertion, encoding="utf-8")
    assert _select(ins_ref, ins_vcf, lacking, "--covers", "symbol", "--family-columns", "IMPACT=") == 0
    assert json.loads((lacking / "manifest.json").read_text())["variant_classes"] == {"insertion": 1}
    assert bgc.main([*check, str(lacking), "--suite", f"s01={ref}={vcf}"]) == 1

    # With family columns `check` reads the input, and refuses an unreadable one after the
    # errors found before it.
    capsys.readouterr()
    absent, matrix = tmp_path / "absent.vcf", _matrix(tmp_path, "other")
    with pytest.raises(SystemExit) as refusal:
        bgc.main([*check, str(out), "--suite", f"s01={ref}={absent}", "--matrix", str(matrix)])
    assert str(refusal.value).startswith(f"ERROR: [build_golden_corpus] {absent} is not a readable text file (")
    assert capsys.readouterr().err == f"ERROR: [build_golden_corpus] covers names symbol, which is not a row of {matrix}\n"

    # The record-class fill adds every row of the records it takes, read with the input too.
    def motif(name: str, location: str, allele: str, feature: str) -> str:
        return _feature_row(name, location, allele, feature, "MotifFeature", "TF_binding_site_variant", "IMPACT=MODIFIER;STRAND=1")

    fill_ref, fill_vcf, fill = tmp_path / "fill.txt", tmp_path / "fill.vcf", tmp_path / "fill"
    fill_ref.write_text(
        motif("snv1", "21:300", "G", "ENSM1")
        + "".join(motif(name, location, "TG", feature) for name, location in (("subA", "21:100-101"), ("subB", "21:500-501")) for feature in ("ENSM1", "ENSM2")),
        encoding="utf-8",
    )
    fill_vcf.write_text(header + "21\t100\tsubA\tAC\tTG\t.\tPASS\t.\n21\t300\tsnv1\tA\tG\t.\tPASS\t.\n21\t500\tsubB\tAC\tTG\t.\tPASS\t.\n", encoding="utf-8")
    assert _select(fill_ref, fill_vcf, fill, "--covers", "symbol", "--record-classes", "regulatory_class", "--record-class-exemplars", "4") == 0
    filled = json.loads((fill / "manifest.json").read_text())
    assert filled["variant_classes"] == {"mnv_or_complex": 3, "snv": 1}
    assert filled["record_classes"]["motif_mnv"] == 2

    alts = bgc.substitution_alts(vcf)
    assert alts == {("21", "mnv1"): [(100, "TG")], ("21", "21_300_AC/TG"): [(300, "TG")]}
    assert [r.klass for r in bgc.iter_rows(ref, "s01", substitutions=alts)] == ["mnv_or_complex", "insertion", "mnv_or_complex"]
    assert [r.klass for r in bgc.iter_rows(ref, "s01")] == ["insertion"] * 3
    assert bgc.variant_class("21:100-101", "TG", True) == "mnv_or_complex"
    assert bgc.holds_allele("100-101", "tg", [(100, "TG")])
    assert not bgc.holds_allele("100-101", "CA", [(100, "TG")])
    assert not bgc.holds_allele("100-101", "TGA", [(100, "TG")]), "an allele of another length"
    assert not bgc.holds_allele("101-102", "G", [(100, "TG")]), "an allele shorter than its span"
    assert not bgc.holds_allele("x-y", "TG", [(100, "TG")])
    # Bases the alleles share are dropped: the anchor base of a multi-allelic record with an
    # indel (ATG>AG,ACG prints 301-302 CG), either end under --minimal (AAAG>AGGG prints
    # 401-402 GG). A record that only shares the row's name, here an insertion, is not read.
    shared_ref, shared_vcf = tmp_path / "shared.txt", tmp_path / "shared.vcf"
    shared_ref.write_text(
        "".join(_row("rs1", location, allele, "ENST1", "coding_sequence_variant")
                for location, allele in (("21:100-101", "TG"), ("21:200-201", "TG"), ("21:301-302", "CG"), ("21:401-402", "GG"))),
        encoding="utf-8",
    )
    shared_vcf.write_text(
        header + "".join(f"21\t{pos}\trs1\t{ref_allele}\t{alt}\t.\tPASS\t.\n"
                         for pos, ref_allele, alt in ((100, "AC", "TG"), (200, "A", "ATG"), (300, "ATG", "AG,ACG,<DEL>,*"), (400, "AAAG", "AGGG"))),
        encoding="utf-8",
    )
    shared = bgc.substitution_alts(shared_vcf)
    assert shared == {("21", "rs1"): [(100, "TG"), (300, "ACG"), (400, "AGGG")]}
    assert [r.klass for r in bgc.iter_rows(shared_ref, "s01", substitutions=shared)] == ["mnv_or_complex", "insertion", "mnv_or_complex", "mnv_or_complex"]
    # The motif allele shape reads a row its record's pairs do not match from the row alone.
    trimmed = bgc.Row("s01", "rs1", "21:401-402", "GG", "ENSM1", "MotifFeature", "TF_binding_site_variant", "1", "", ("-", "-", "-"), substitution=True)
    assert trimmed.klass == "mnv_or_complex" and bgc.allele_shape(trimmed, [("AAAG", "AGGG")]) == "indel"


def _row(uv, loc, allele, feat, cons, extra="IMPACT=MODIFIER;STRAND=1"):
    return "\t".join([uv, loc, allele, "ENSG1", feat, "Transcript", cons, "-", "-", "-", "-", "-", "-", extra]) + "\n"


def test_select_builds_corpus_and_manifest(tmp_path: Path):
    ref = tmp_path / "ref.txt"
    ref.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + _row("rs1", "21:100", "G", "ENST2", "intron_variant")
        + _row("21_201_AT/A", "21:201", "-", "ENST1", "frameshift_variant", "IMPACT=HIGH;STRAND=1;FLAGS=cds_end_NF")
        + _row("21_300_A/G/T", "21:300", "G", "ENST3", "synonymous_variant")
        + _row("21_300_A/G/T", "21:300", "T", "ENST3", "missense_variant"),
        encoding="utf-8",
    )
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\trs1\tA\tG\t.\tPASS\t.\n"
        "21\t150\trs9\tC\tT\t.\tPASS\t.\n"
        "21\t200\t.\tAT\tA\t.\tPASS\t.\n"
        "21\t300\t.\tA\tG,T\t.\tPASS\t.\n"
        "21\t400\trs4\tCA\tC,CGG\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    out = tmp_path / "corpus"
    rc = bgc.main(["select", "--release", "115", "--assembly", "GRCh37", "--suite", f"s01={ref}={vcf}", "--out", str(out), "--k", "2", "--seed", "3", "--multi-allelic-per-suite", "0", "--covers", "symbol"])
    assert rc == 0
    m = json.loads((out / "manifest.json").read_text())
    assert set(m["combinations"]) == {"missense_variant", "intron_variant", "frameshift_variant", "synonymous_variant"}
    assert m["uncovered_combinations"] == [] and m["unresolved_exemplars"] == []
    records = [r["vcf"] for r in m["records"]]
    assert records == ["21\t100\trs1\tA\tG\t.\tPASS\t.", "21\t200\t.\tAT\tA\t.\tPASS\t.", "21\t300\t.\tA\tG,T\t.\tPASS\t."]
    body = (out / "variants.vcf").read_text().splitlines()
    assert body[0] == "##fileformat=VCFv4.2" and body[1].startswith("#CHROM") and len(body) == 5
    assert "rs9" not in "".join(body)  # the unselected record is not carried

    # With multi-allelic sampling on, the SNV-only and the indel multi-allelic
    # records both join the corpus, tagged by stratum, without claiming combinations.
    out2 = tmp_path / "corpus2"
    rc = bgc.main(["select", "--release", "115", "--assembly", "GRCh37", "--suite", f"s01={ref}={vcf}", "--out", str(out2), "--k", "2", "--seed", "3", "--multi-allelic-per-suite", "2", "--covers", "symbol"])
    assert rc == 0
    m2 = json.loads((out2 / "manifest.json").read_text())
    by_vcf = {r["vcf"]: r for r in m2["records"]}
    assert by_vcf["21\t300\t.\tA\tG,T\t.\tPASS\t."]["strata"] == ["multi_allelic_snv"]
    assert by_vcf["21\t400\trs4\tCA\tC,CGG\t.\tPASS\t."]["strata"] == ["multi_allelic_indel"]
    assert by_vcf["21\t400\trs4\tCA\tC,CGG\t.\tPASS\t."]["combinations"] == []


def test_select_writes_the_species_and_the_cache_dir_follows_it(tmp_path: Path):
    """`--species` lands in the manifest and names the Ensembl-layout cache directory; a
    manifest without the key is homo_sapiens, and a species that is not a production name
    is refused with the manifest's other shape errors."""
    ref = tmp_path / "ref.txt"
    ref.write_text(_row("rs1", "IV:100", "G", "YDL001", "missense_variant"), encoding="utf-8")
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nIV\t100\trs1\tA\tG\t.\tPASS\t.\n", encoding="utf-8")
    common = ["--suite", f"s01={ref}={vcf}", "--k", "1", "--multi-allelic-per-suite", "0", "--covers", "symbol"]
    yeast = tmp_path / "yeast"
    assert bgc.main(["select", "--release", "116", "--assembly", "R64-1-1", "--species", "saccharomyces_cerevisiae", "--out", str(yeast)] + common) == 0
    m = json.loads((yeast / "manifest.json").read_text())
    assert m["species"] == "saccharomyces_cerevisiae" and m["assembly"] == "R64-1-1"
    assert bgc.corpus_cache_dir(yeast, m) == yeast / "json_cache"
    m["cache_layout"] = "ensembl"
    assert bgc.corpus_cache_dir(yeast, m) == yeast / "cache_root" / "saccharomyces_cerevisiae" / "116_R64-1-1"
    human = tmp_path / "human"
    assert bgc.main(["select", "--release", "116", "--assembly", "GRCh38", "--out", str(human)] + common) == 0
    m = json.loads((human / "manifest.json").read_text())
    assert m["species"] == "homo_sapiens"
    del m["species"]
    m["cache_layout"] = "ensembl"
    assert bgc.corpus_cache_dir(human, m) == human / "cache_root" / "homo_sapiens" / "116_GRCh38"
    assert bgc.manifest_shape_error(m) is None
    m["species"] = "Saccharomyces cerevisiae"
    assert bgc.manifest_shape_error(m) == 'species is "Saccharomyces cerevisiae", not an Ensembl species name'
    m["species"] = 7
    assert bgc.manifest_shape_error(m) == "species is 7, not an Ensembl species name"
    assert bgc.main(["select", "--release", "116", "--assembly", "GRCh38", "--species", "Homo sapiens", "--out", str(tmp_path / "bad")] + common) == 1
    assert not (tmp_path / "bad").exists()
    # A regeneration over the yeast manifest keeps its species without --species.
    assert bgc.main(["select", "--release", "116", "--assembly", "R64-1-1", "--out", str(yeast)] + common) == 0
    assert json.loads((yeast / "manifest.json").read_text())["species"] == "saccharomyces_cerevisiae"


def _feature_row(uv, loc, allele, feat, ftype, cons, extra="IMPACT=MODIFIER"):
    return "\t".join([uv, loc, allele, "-", feat, ftype, cons, "-", "-", "-", "-", "-", "-", extra]) + "\n"


def test_regulatory_class_key_covers_every_class(tmp_path: Path):
    """`--record-classes regulatory_class` classes each record by its regulatory rows, adds
    records until every class the reference shows has the asked exemplars, records the classes
    per record and per class, and `check` recomputes the coverage from the reference."""
    rf = "regulatory_region_variant"
    tf = "TF_binding_site_variant"
    ref = tmp_path / "ref.txt"
    ref.write_text(
        # Two records of one shape: the combination sampling takes one, and `rf_only` reaches its
        # third exemplar through the other.
        _row("rs_rf", "21:1000", "G", "ENST1", "intron_variant")
        + _feature_row("rs_rf", "21:1000", "G", "ENSR1", "RegulatoryFeature", rf)
        + _row("rs_rf2", "21:1100", "G", "ENST1", "intron_variant")
        + _feature_row("rs_rf2", "21:1100", "G", "ENSR1", "RegulatoryFeature", rf)
        + _row("rs_mf", "21:2000", "G", "ENST1", "intron_variant")
        + _feature_row("rs_mf", "21:2000", "G", "ENSM1", "MotifFeature", tf, "IMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM1")
        + _row("rs_both", "21:3000", "G", "ENST1", "intron_variant")
        + _feature_row("rs_both", "21:3000", "G", "ENSR2", "RegulatoryFeature", rf)
        + _feature_row("rs_both", "21:3000", "G", "ENSM2", "MotifFeature", tf, "IMPACT=MODIFIER;STRAND=-1;MOTIF_NAME=ENSPFM2")
        + _feature_row("rs_regonly", "21:4000", "G", "ENSR3", "RegulatoryFeature", rf)
        + _feature_row("rs_regonly", "21:4000", "G", "-", "-", "intergenic_variant")
        + _row("21_5000_AT/GC", "21:5000-5001", "GC", "ENST1", "intron_variant")
        + _feature_row("21_5000_AT/GC", "21:5000-5001", "GC", "ENSM3", "MotifFeature", tf, "IMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM3")
        + _row("21_6000_AT/A", "21:6001", "-", "ENST1", "intron_variant")
        + _feature_row("21_6000_AT/A", "21:6001", "-", "ENSR4", "RegulatoryFeature", rf)
        + _feature_row("21_6000_AT/A", "21:6001", "-", "ENSM4", "MotifFeature", tf, "IMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM4")
        + _row("rs_plain", "21:7000", "G", "ENST1", "intron_variant")
        # A multi-allelic record the stratum sampling carries: its classes are read by its names.
        + _feature_row("21_9000_A/G/T", "21:9000", "G", "ENSM5", "MotifFeature", tf, "IMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM5")
        + _feature_row("21_9000_A/G/T", "21:9000", "G", "-", "-", "intergenic_variant")
        + _feature_row("21_9000_A/G/T", "21:9000", "T", "ENSM5", "MotifFeature", tf, "IMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM5")
        + _feature_row("21_9000_A/G/T", "21:9000", "T", "-", "-", "intergenic_variant"),
        encoding="utf-8",
    )
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t1000\trs_rf\tA\tG\t.\tPASS\t.\n"
        "21\t1100\trs_rf2\tA\tG\t.\tPASS\t.\n"
        "21\t2000\trs_mf\tA\tG\t.\tPASS\t.\n"
        "21\t3000\trs_both\tA\tG\t.\tPASS\t.\n"
        "21\t4000\trs_regonly\tA\tG\t.\tPASS\t.\n"
        "21\t5000\t.\tAT\tGC\t.\tPASS\t.\n"
        "21\t6000\t.\tAT\tA\t.\tPASS\t.\n"
        "21\t7000\trs_plain\tA\tG\t.\tPASS\t.\n"
        "21\t9000\t.\tA\tG,T\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    out = tmp_path / "corpus"
    common = ["--release", "116", "--assembly", "GRCh38", "--suite", f"reg={ref}={vcf}", "--out", str(out),
              "--k", "1", "--seed", "1", "--multi-allelic-per-suite", "1", "--covers", "regulatory"]
    assert bgc.main(["select", *common, "--record-classes", "regulatory_class", "--record-class-exemplars", "3"]) == 0
    m = json.loads((out / "manifest.json").read_text())
    assert m["record_class_key"] == "regulatory_class" and m["record_class_exemplars"] == 3
    classes = {f[2] if f[2] != "." else f[1]: r["record_classes"] for r in m["records"] for f in [r["vcf"].split("\t")]}
    assert classes["rs_rf"] == ["rf_only"] and classes["rs_rf2"] == ["rf_only"]
    assert classes["rs_mf"] == ["mf_only", "motif_snv"]
    assert classes["rs_both"] == ["motif_reverse", "motif_snv", "rf_and_mf"]
    assert classes["rs_regonly"] == ["regulatory_only", "rf_only"]
    assert classes["5000"] == ["mf_only", "motif_mnv"]
    assert classes["6000"] == ["motif_indel", "rf_and_mf"]
    assert classes["9000"] == ["mf_only", "motif_snv", "regulatory_only"]
    assert "rs_plain" not in classes or classes["rs_plain"] == []
    assert m["record_classes"] == {"mf_only": 3, "motif_indel": 1, "motif_mnv": 1, "motif_reverse": 1, "motif_snv": 3,
                                   "regulatory_only": 2, "rf_and_mf": 2, "rf_only": 3}
    assert m["uncovered_record_classes"] == []
    # With one exemplar per class the multi-allelic record enters through its stratum alone,
    # with no exemplar row, and its classes are still read from the reference by its names.
    one = tmp_path / "one"
    assert bgc.main(["select", *(str(one) if a == str(out) else a for a in common), "--record-classes", "regulatory_class",
                     "--record-class-exemplars", "1"]) == 0
    sampled = next(r for r in json.loads((one / "manifest.json").read_text())["records"] if r["vcf"].startswith("21\t9000\t"))
    assert sampled["combinations"] == [] and sampled["strata"] == ["multi_allelic_snv"]
    assert sampled["record_classes"] == ["mf_only", "motif_snv", "regulatory_only"]
    assert bgc.regulatory_classes([], []) == set()
    # A two-base substitution and a two-base insertion print one row shape: the input's alleles
    # tell them apart, and a row matching no pair reads as the row alone does.
    mnv_row = bgc.Row("reg", "x", "21:5000-5001", "GC", "ENSM3", "MotifFeature", "TF_binding_site_variant", "1", "", ("-", "-", "-"))
    assert bgc.allele_shape(mnv_row, [("AT", "GC")]) == "mnv"
    assert bgc.allele_shape(mnv_row, [("A", "AGC")]) == "indel"
    assert bgc.allele_shape(mnv_row, []) == "indel"
    assert bgc.allele_shape(mnv_row, [("A", "<DEL>"), ("AT", "GC")]) == "mnv"
    del_row = bgc.Row("reg", "x", "21:6001", "-", "ENSM4", "MotifFeature", "TF_binding_site_variant", "1", "", ("-", "-", "-"))
    assert bgc.allele_shape(del_row, [("AT", "A")]) == "indel"
    # A multi-allelic record with an indel loses the first base of every allele and nothing more.
    multi = [("AT", "A"), ("AT", "ATT")]
    tt_row = bgc.Row("reg", "x", "21:1001", "TT", "ENSM4", "MotifFeature", "TF_binding_site_variant", "1", "", ("-", "-", "-"))
    assert bgc.allele_shape(tt_row, multi) == "indel"
    assert bgc.allele_shape(bgc.Row("reg", "x", "21:1001", "-", "ENSM4", "MotifFeature", "TF_binding_site_variant", "1", "", ("-", "-", "-")), multi) == "indel"
    cg_row = bgc.Row("reg", "x", "21:1001-1002", "CG", "ENSM4", "MotifFeature", "TF_binding_site_variant", "1", "", ("-", "-", "-"))
    assert bgc.allele_shape(cg_row, [("ATG", "AG"), ("ATG", "ACG")]) == "mnv"
    g_row = bgc.Row("reg", "x", "21:1000", "G", "ENSM4", "MotifFeature", "TF_binding_site_variant", "1", "", ("-", "-", "-"))
    assert bgc.allele_shape(g_row, [("A", "AT"), ("A", "G")]) == "snv"
    star_row = bgc.Row("reg", "x", "21:1001", "*", "ENSM4", "MotifFeature", "TF_binding_site_variant", "1", "", ("-", "-", "-"))
    assert bgc.allele_shape(star_row, [("AT", "A"), ("AT", "*")]) == "indel"
    assert bgc.main(["check", "--corpus", str(out), "--suite", f"reg={ref}={vcf}"]) == 0

    # A class the reference shows and no corpus record carries fails the check.
    wider = tmp_path / "wider.txt"
    wider.write_text(ref.read_text() + _feature_row("rs_extra", "21:8000", "G", "ENSM9", "MotifFeature", tf, "IMPACT=MODIFIER;STRAND=-1")
                     + _feature_row("rs_extra", "21:8000", "G", "-", "-", "intergenic_variant"), encoding="utf-8")
    m["records"] = [r for r in m["records"] if "rs_both" not in r["vcf"]]
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    assert bgc.main(["check", "--corpus", str(out), "--suite", f"reg={wider}={vcf}"]) == 1

    with pytest.raises(SystemExit, match="--record-classes wants one of regulatory_class"):
        bgc.main(["select", *common, "--record-classes", "splice_class"])
    m["record_class_key"] = "splice_class"
    assert bgc.manifest_shape_error(m) == "record_class_key is not one of regulatory_class"


def test_multi_allelic_strata():
    assert bgc.multi_allelic_stratum("A", ["G", "T"]) == "multi_allelic_snv"
    assert bgc.multi_allelic_stratum("CA", ["C", "CGG"]) == "multi_allelic_indel"  # chopped: A/-/GG, nothing left to trim
    assert bgc.multi_allelic_stratum("CA", ["C", "CAA"]) == "multi_allelic_indel_trimmable"  # chopped: A/-/AA, A shared
    # chopped TCACACA/TACACACACACACACACA -> CACACA/ACACACACACACACACA share their last base
    assert bgc.multi_allelic_stratum("TCACACA", ["TACACACACACACACACA", "T"]) == "multi_allelic_indel_trimmable"
    assert bgc.multi_allelic_stratum("CTTT", ["C", "CTT"]) == "multi_allelic_indel_trimmable"  # TTT/-/TT share T


def _record(vcf: str, **extra) -> dict:
    return {"suite": "s01", "vcf": vcf, "combinations": [], **extra}


def _manifest(corpus: Path, records: list[dict]) -> None:
    corpus.mkdir(exist_ok=True)
    (corpus / "manifest.json").write_text(json.dumps({"records": records}), encoding="utf-8")


def test_classify_marks_documented_and_unexplained_divergences(tmp_path: Path):
    """Each divergence names the records its key resolves to; the ID-less indel resolves
    only through release 116's raw-line name (`21_400_AT/A`, where release 115 wrote
    `21_401_AT/A`); a key only vep-rs emits for a record VEP annotated is an
    `unexplained_residual`."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\trs1\tA\tG\t.\tPASS\t."),
        _record("21\t200\trs2\tA\tT\t.\tPASS\t."),
        _record("21\t300\trs3\tA\tC\t.\tPASS\t."),
        _record("21\t400\t.\tAT\tA\t.\tPASS\t."),
    ])
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + _row("rs2", "21:200", "T", "ENST2", "start_lost,start_retained_variant")
        + _row("rs3", "21:300", "C", "ENST3", "intron_variant")
        + _row("21_400_AT/A", "21:401", "-", "ENST4", "frameshift_variant"),
        encoding="utf-8",
    )
    rs.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + _row("rs2", "21:200", "T", "ENST2", "start_retained_variant")
        + _row("rs3", "21:300", "C", "ENST3", "intron_variant,splice_region_variant")
        + _row("21_400_AT/A", "21:401", "-", "ENST4", "frameshift_variant")
        + _row("21_400_AT/A", "21:401", "-", "ENST5", "intron_variant"),
        encoding="utf-8",
    )
    rc = bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])
    assert rc == 0
    m = json.loads((corpus / "manifest.json").read_text())
    by_loc = {d["location"]: (d["expected_divergence"], d["record_indices"]) for d in m["divergences"]}
    assert by_loc["21:200"] == ("start_cooccurrence_swap", [1])
    assert by_loc["21:300"] == ("splice_family_swap", [2])
    assert "21:100" not in by_loc and "21:401" not in by_loc
    assert m["divergence_summary"] == {"start_cooccurrence_swap": 1, "splice_family_swap": 1}
    assert m["vep_rs_only_tuples"] == [{
        "location": "21:401", "allele": "-", "feature": "ENST5", "feature_type": "Transcript",
        "record_indices": [3], "vep_rs_consequence_sets": ["intron_variant"],
        "expected_divergence": "unexplained_residual",
    }]
    assert m["vep_rs_only_summary"] == {"unexplained_residual": 1}
    assert not any("reference_rows" in r for r in m["records"])


WITH_REFERENCE = "vep_rs_consequence_sets_with_reference_sequence"


def _reference_key_corpus(tmp_path: Path) -> tuple[Path, Path, Path]:
    """A two-record corpus with the two outputs classify reads: both engines' rows differ on
    each record, by a pair the registry classes name."""
    corpus = tmp_path / "c"
    _manifest(corpus, [_record("21\t200\trs2\tA\tT\t.\tPASS\t."), _record("21\t300\trs3\tA\tC\t.\tPASS\t.")])
    vep, rs = tmp_path / "vep.txt", tmp_path / "rs.txt"
    vep.write_text(
        _row("rs2", "21:200", "T", "ENST2", "start_lost,start_retained_variant") + _row("rs3", "21:300", "C", "ENST3", "intron_variant"),
        encoding="utf-8",
    )
    rs.write_text(
        _row("rs2", "21:200", "T", "ENST2", "start_retained_variant") + _row("rs3", "21:300", "C", "ENST3", "intron_variant,splice_region_variant"),
        encoding="utf-8",
    )
    return corpus, vep, rs


def test_classify_carries_a_with_reference_key_through_a_regeneration(tmp_path: Path, capsys):
    """classify writes no `vep_rs_consequence_sets_with_reference_sequence` of its own, and a
    key written by hand survives it: the entry classify writes for the same Location, Allele,
    Feature, records and Ensembl VEP set carries the key unchanged, every other entry and both
    summaries as a run without the key writes them. A flag set's own list carries its key the
    same way. A run whose outputs give no such entry, or whose regenerated sets the key
    repeats, stops with the reason and leaves the manifest as it was."""
    corpus, vep, rs = _reference_key_corpus(tmp_path)
    classify = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    path = corpus / "manifest.json"
    assert bgc.main(classify) == 0
    first = json.loads(path.read_text())
    assert [d["location"] for d in first["divergences"]] == ["21:200", "21:300"]
    assert not any(WITH_REFERENCE in d for d in first["divergences"])

    stated = json.loads(path.read_text())
    stated["divergences"][0][WITH_REFERENCE] = ["start_lost,start_retained_variant"]
    path.write_text(json.dumps(stated), encoding="utf-8")
    assert bgc.main(classify) == 0
    second = json.loads(path.read_text())
    assert second["divergences"][0] == {**first["divergences"][0], WITH_REFERENCE: ["start_lost,start_retained_variant"]}
    assert list(second["divergences"][0])[-1] == WITH_REFERENCE
    assert second["divergences"][1] == first["divergences"][1]
    assert {k: v for k, v in second.items() if k != "divergences"} == {k: v for k, v in first.items() if k != "divergences"}
    carried = path.read_text()

    # A flag set's own list: its key is carried by `--flag-set`, and the manifest's own stays.
    with_set = json.loads(carried)
    with_set["flag_sets"] = [{"name": "again", "flags": ["--again"]}]
    path.write_text(json.dumps(with_set), encoding="utf-8")
    assert bgc.main([*classify, "--flag-set", "again"]) == 0
    in_set = json.loads(path.read_text())
    assert not any(WITH_REFERENCE in d for d in in_set["flag_sets"][0]["divergences"])
    in_set["flag_sets"][0]["divergences"][1][WITH_REFERENCE] = ["intron_variant"]
    path.write_text(json.dumps(in_set), encoding="utf-8")
    assert bgc.main([*classify, "--flag-set", "again"]) == 0
    in_set = json.loads(path.read_text())
    assert in_set["flag_sets"][0]["divergences"][1][WITH_REFERENCE] == ["intron_variant"]
    assert WITH_REFERENCE not in in_set["flag_sets"][0]["divergences"][0]
    assert in_set["divergences"] == second["divergences"]
    path.write_text(carried, encoding="utf-8")
    capsys.readouterr()

    # The run's outputs agree on the key: no entry is written to carry the statement on.
    agreeing = tmp_path / "agreeing.txt"
    agreeing.write_text(rs.read_text().replace("\tstart_retained_variant\t", "\tstart_lost,start_retained_variant\t"), encoding="utf-8")
    assert bgc.main([*classify[:-1], str(agreeing)]) == 1
    assert capsys.readouterr().err == (
        f"ERROR: [build_golden_corpus] the entry on 21:200 T ENST2 carries {WITH_REFERENCE}, and this run writes no entry of its "
        "Location, Allele, Feature, records and Ensembl VEP set to carry it on: remove the key, or classify the "
        "outputs of a run without a reference sequence\n"
    )
    assert path.read_text() == carried

    # The regenerated sets are the stated ones: the key states no second outcome.
    repeating = json.loads(carried)
    repeating["divergences"][0][WITH_REFERENCE] = ["start_lost"]
    path.write_text(json.dumps(repeating), encoding="utf-8")
    kept = path.read_text()
    lost = tmp_path / "lost.txt"
    lost.write_text(rs.read_text().replace("\tstart_retained_variant\t", "\tstart_lost\t"), encoding="utf-8")
    assert bgc.main([*classify[:-1], str(lost)]) == 1
    assert capsys.readouterr().err == (
        f"ERROR: [build_golden_corpus] the entry on 21:200 T ENST2: {WITH_REFERENCE} "
        "repeats vep_rs_consequence_sets and states no second outcome\n"
    )
    assert path.read_text() == kept


def test_a_malformed_with_reference_key_is_refused_by_every_reader(tmp_path: Path):
    """The key is a non-empty sorted list of distinct consequence sets, other than the entry's
    `vep_rs_consequence_sets`, that names Ensembl VEP's set alone or not at all, on a
    `divergences` entry: the builder refuses any other shape wherever it reads a manifest,
    its own run's list and a flag set's alike, before anything is written, and the format
    comparator refuses the same shapes for the same reasons."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    assert cvf.WITH_REFERENCE_SETS == bgc.WITH_REFERENCE_SETS == WITH_REFERENCE
    corpus, vep, rs = _reference_key_corpus(tmp_path)
    classify = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    path = corpus / "manifest.json"
    assert bgc.main(classify) == 0
    good = json.loads(path.read_text())
    entry = good["divergences"][0]
    vep_set, own = entry["vep_consequence_set"], entry["vep_rs_consequence_sets"]
    assert (vep_set, own) == ("start_lost,start_retained_variant", ["start_retained_variant"])
    for value in ([vep_set], ["intron_variant"], ["intron_variant", "start_lost"]):
        stated = {**entry, WITH_REFERENCE: value}
        assert bgc.with_reference_sets_error(stated) is None and cvf.with_reference_sets(stated) == value
    assert bgc.with_reference_sets_error(entry) is None and cvf.with_reference_sets(entry) is None
    shapes = [
        (vep_set, "is not a non-empty list"),
        ([], "is not a non-empty list"),
        ({"vep_rs_consequence_sets": [vep_set]}, "is not a non-empty list"),
        ([vep_set, None], "holds an item that is not a sorted, comma-joined set of terms"),
        (["start_retained_variant,start_lost"], "holds an item that is not a sorted, comma-joined set of terms"),
        ([""], "holds an item that is not a sorted, comma-joined set of terms"),
        (["start_lost", "intron_variant"], "is not a sorted list of distinct sets"),
        (["intron_variant", "intron_variant"], "is not a sorted list of distinct sets"),
        (own, "repeats vep_rs_consequence_sets and states no second outcome"),
        (["intron_variant", vep_set], "names Ensembl VEP's set beside another"),
    ]
    commands = [classify, ["covers", str(corpus), "--list", "--matrix", str(_matrix(tmp_path, "symbol"))]]
    for value, reason in shapes:
        stated = {**entry, WITH_REFERENCE: value}
        assert bgc.with_reference_sets_error(stated) == reason
        with pytest.raises(cvf.ComparatorError, match=re.escape(f"the manifest's item on 21:200 T ENST2: {WITH_REFERENCE} {reason}")):
            cvf.with_reference_sets(stated)
        for owner, manifest in (
            ("the manifest's", {**good, "divergences": [stated, *good["divergences"][1:]]}),
            ("flag set other", {**good, "flag_sets": [{"name": "other", "divergences": [stated]}]}),
        ):
            path.write_text(json.dumps(manifest), encoding="utf-8")
            written = path.read_text()
            for command in commands:
                with pytest.raises(SystemExit, match=re.escape(f"{owner} divergences entry on 21:200 T ENST2: {WITH_REFERENCE} {reason}")):
                    bgc.main(command)
                assert path.read_text() == written
    # The key belongs on a `divergences` entry: the other two lists state no Ensembl VEP set.
    for section in ("vep_rs_only_tuples", "field_divergences"):
        path.write_text(json.dumps({**good, section: [{"location": "21:200", "allele": "T", "feature": "ENST9", WITH_REFERENCE: [vep_set]}]}), encoding="utf-8")
        with pytest.raises(SystemExit, match=re.escape(f"the manifest's {section} entry on 21:200 T ENST9: {WITH_REFERENCE} belongs on a divergences entry")):
            bgc.main(classify)


def test_classify_leaves_a_reference_row_without_a_vep_rs_row_unexplained(tmp_path: Path):
    """A row Ensembl VEP writes on a key vep-rs writes no row for has no documented class,
    whatever the shape of its allele: a symbolic allele's row, a breakend's and a sequence
    allele's are each an `unexplained_residual` that lists no vep-rs set."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=500"),
        _record("21\t600\tbnd1\tA\tA[21:900[\t.\tPASS\tSVTYPE=BND"),
        _record("21\t700\trs1\tA\tG\t.\tPASS\t."),
    ])
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    shared = (
        _row("sv1", "21:101-500", "deletion", "ENST1", "transcript_ablation")
        + _row("bnd1", "21:600", "A[21:900[", "ENST1", "feature_truncation")
        + _row("rs1", "21:700", "G", "ENST1", "intron_variant")
    )
    vep.write_text(
        shared
        + _row("sv1", "21:101-500", "deletion", "ENST2", "intron_variant")
        + _row("bnd1", "21:600", "A[21:900[", "ENST2", "feature_truncation")
        + _row("rs1", "21:700", "G", "ENST2", "intron_variant"),
        encoding="utf-8",
    )
    rs.write_text(shared, encoding="utf-8")
    rc = bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])
    assert rc == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert [
        (d["location"], d["allele"], d["feature"], d["record_indices"], d["vep_rs_consequence_sets"],
         d["expected_divergence"])
        for d in m["divergences"]
    ] == [
        ("21:101-500", "deletion", "ENST2", [0], [], "unexplained_residual"),
        ("21:600", "A[21:900[", "ENST2", [1], [], "unexplained_residual"),
        ("21:700", "G", "ENST2", [2], [], "unexplained_residual"),
    ]
    assert m["divergence_summary"] == {"unexplained_residual": 3}
    assert m["vep_rs_only_tuples"] == []


def test_classify_names_the_motif_score_and_cell_type_field_classes(tmp_path: Path):
    """A `MOTIF_SCORE_CHANGE` difference on a reverse-strand motif row whose allele spans two
    or more bases is `motif_score_reverse_strand_mnv`, an empty Ensembl VEP `CELL_TYPE` beside
    a vep-rs activity under `--cell_type` is `cell_type_comma_names`, and every other field
    difference stays an `unexplained_residual`."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\trs_rev\tCT\tAG\t.\tPASS\t."),
        _record("21\t200\trs_fwd\tCT\tAG\t.\tPASS\t."),
        _record("21\t300\trs_snv\tA\tG\t.\tPASS\t."),
        _record("21\t400\trs_rf\tA\tG\t.\tPASS\t."),
    ])
    m = json.loads((corpus / "manifest.json").read_text())
    m["flag_sets"] = [{"name": "comma", "flags": ["--regulatory", "--cell_type", "A549_(m,_58_y)"]}]
    (corpus / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    motif = "IMPACT=MODIFIER;STRAND={};MOTIF_NAME=ENSPFM1;MOTIF_POS=9;HIGH_INF_POS=N;MOTIF_SCORE_CHANGE={}"
    tf = "TF_binding_site_variant"
    rf = "regulatory_region_variant"
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(
        _feature_row("rs_rev", "21:100-101", "AG", "ENSM1", "MotifFeature", tf, motif.format(-1, "-0.092"))
        + _feature_row("rs_fwd", "21:200-201", "AG", "ENSM2", "MotifFeature", tf, motif.format(1, "-0.016"))
        + _feature_row("rs_snv", "21:300", "G", "ENSM3", "MotifFeature", tf, motif.format(-1, "-0.089"))
        + _feature_row("rs_rf", "21:400", "G", "ENSR1", "RegulatoryFeature", rf, "IMPACT=MODIFIER;CELL_TYPE="),
        encoding="utf-8",
    )
    rs.write_text(
        _feature_row("rs_rev", "21:100-101", "AG", "ENSM1", "MotifFeature", tf, motif.format(-1, "-0.136"))
        + _feature_row("rs_fwd", "21:200-201", "AG", "ENSM2", "MotifFeature", tf, motif.format(1, "-0.017"))
        + _feature_row("rs_snv", "21:300", "G", "ENSM3", "MotifFeature", tf, motif.format(-1, "-0.090"))
        + _feature_row("rs_rf", "21:400", "G", "ENSR1", "RegulatoryFeature", rf, "IMPACT=MODIFIER;CELL_TYPE=A549_(m,_58_y):ACTIVE"),
        encoding="utf-8",
    )
    args = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
            "--field", "MOTIF_SCORE_CHANGE", "--field", "CELL_TYPE"]
    assert bgc.main([*args, "--flag-set", "comma"]) == 0
    fs = json.loads((corpus / "manifest.json").read_text())["flag_sets"][0]
    classes = {d["location"]: d["expected_divergence"] for d in fs["field_divergences"]}
    assert classes == {"21:100-101": "motif_score_reverse_strand_mnv", "21:200-201": "unexplained_residual",
                       "21:300": "unexplained_residual", "21:400": "cell_type_comma_names"}
    assert fs["field_divergence_summary"] == {"cell_type_comma_names": 1, "motif_score_reverse_strand_mnv": 1, "unexplained_residual": 2}
    # Without --cell_type among the run's flags the empty reference value is not the comma class.
    assert bgc.main([*args, "--flag=--regulatory"]) == 0
    top = json.loads((corpus / "manifest.json").read_text())
    assert {d["location"]: d["expected_divergence"] for d in top["field_divergences"]}["21:400"] == "unexplained_residual"


def test_classify_names_the_rejoined_record_plugin_field_class(tmp_path: Path):
    """A covered plugin's field (`rejoined_record_plugin_field`) that differs on a `MINIMISED` row
    of a multi-allelic record under `--minimal` is `minimal_plugin_sees_the_rejoined_record`; the
    same difference on a bi-allelic record, on a multi-allelic record the flag left unsplit (no
    `MINIMISED` on its rows), under a run without `--minimal`, or in a field no covered plugin
    writes stays an `unexplained_residual`."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\trs_multi\tC\tCCA,T\t.\tPASS\t."),
        _record("21\t200\trs_bi\tC\tT\t.\tPASS\t."),
        _record("21\t300\trs_unsplit\tG\tT,C\t.\tPASS\t."),
    ])
    m = json.loads((corpus / "manifest.json").read_text())
    m["flag_sets"] = [
        {"name": "minimal", "flags": ["--plugin", "gnomADc,cov.tsv.gz", "--minimal", "--allele_number"]},
        {"name": "plain", "flags": ["--plugin", "gnomADc,cov.tsv.gz", "--allele_number"]},
    ]
    (corpus / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    intron = "intron_variant"
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(
        _feature_row("rs_multi", "21:100", "CA", "ENST1", "Transcript", intron, "ALLELE_NUM=1;MINIMISED=1;gnomAD_mean_cov=1.7938;NOTE=x")
        + _feature_row("rs_multi", "21:100", "T", "ENST1", "Transcript", intron, "ALLELE_NUM=2;MINIMISED=1;gnomAD_mean_cov=1.7938;NOTE=x")
        + _feature_row("rs_bi", "21:200", "T", "ENST1", "Transcript", intron, "ALLELE_NUM=1;MINIMISED=1;gnomAD_mean_cov=2.0000;NOTE=x")
        + _feature_row("rs_unsplit", "21:300", "T", "ENST1", "Transcript", intron, "ALLELE_NUM=1;gnomAD_mean_cov=3.0000;NOTE=x"),
        encoding="utf-8",
    )
    rs.write_text(
        _feature_row("rs_multi", "21:100", "CA", "ENST1", "Transcript", intron, "ALLELE_NUM=1;MINIMISED=1;gnomAD_mean_cov=1.7937;NOTE=y")
        + _feature_row("rs_multi", "21:100", "T", "ENST1", "Transcript", intron, "ALLELE_NUM=2;MINIMISED=1;gnomAD_mean_cov=1.7938;NOTE=x")
        + _feature_row("rs_bi", "21:200", "T", "ENST1", "Transcript", intron, "ALLELE_NUM=1;MINIMISED=1;gnomAD_mean_cov=2.5000;NOTE=x")
        + _feature_row("rs_unsplit", "21:300", "T", "ENST1", "Transcript", intron, "ALLELE_NUM=1;gnomAD_mean_cov=3.5000;NOTE=x"),
        encoding="utf-8",
    )
    args = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
            "--field", "gnomAD_mean_cov", "--field", "NOTE"]
    assert bgc.main([*args, "--flag-set", "minimal"]) == 0
    fs = json.loads((corpus / "manifest.json").read_text())["flag_sets"][0]
    classes = {(d["location"], d["allele"], d["field"]): d["expected_divergence"] for d in fs["field_divergences"]}
    assert classes == {
        ("21:100", "CA", "gnomAD_mean_cov"): "minimal_plugin_sees_the_rejoined_record",
        ("21:100", "CA", "NOTE"): "unexplained_residual",
        ("21:200", "T", "gnomAD_mean_cov"): "unexplained_residual",
        ("21:300", "T", "gnomAD_mean_cov"): "unexplained_residual",
    }
    assert fs["field_divergence_summary"] == {"minimal_plugin_sees_the_rejoined_record": 1, "unexplained_residual": 3}
    assert bgc.main([*args, "--flag-set", "plain"]) == 0
    plain = json.loads((corpus / "manifest.json").read_text())["flag_sets"][1]
    assert {d["expected_divergence"] for d in plain["field_divergences"]} == {"unexplained_residual"}


# Made-up scores of two symbols at the deletion 21 5000 GCA G, whose rows stand at
# 21:5001-5002 with the allele written -.
SPLICEAI_SCORES = {"GENEA": "0.15|0.08|0.06|0.17|-6|-41|-17|-9", "GENEB": "0.00|1.00|0.10|0.05|-25|-23|37|17"}
SPLICEAI_WHOLE = f"GENEA|{SPLICEAI_SCORES['GENEA']},G|GENEB|{SPLICEAI_SCORES['GENEB']}"


def _spliceai_score_file(directory: Path, join: bool = True, variant: tuple[str, str, str] = ("5000", "GCA", "G"),
                         records: list[tuple[str, ...]] | None = None) -> Path:
    """A made-up indel score file in the fixture generator's layout, tabix-indexed and each
    line a BGZF block of its own: one record per symbol at 21 5000 GCA G, or, under `join`,
    the one record the generator's `--join-entries` form writes there. `variant` is the
    records' POS, REF and ALT: under a REF and an ALT of one base each the file is the
    single-base one. `records` names the file's record lines instead, each the entries of one
    line at the variant, a symbol and its eight values."""
    sys.path.insert(0, str(HERE.parent / "fixtures"))
    import plugin_score_fixtures as generator

    def block(data: bytes) -> bytes:
        deflater = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = deflater.compress(data) + deflater.flush()
        return (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", len(body) + 25)
                + body + struct.pack("<II", zlib.crc32(data), len(data)))

    position, reference, alternate = variant
    lines = ["\t".join(["21", position, ".", reference, alternate, ".", ".", f"SpliceAI={alternate}|{symbol}|{scores}"]) for symbol, scores in SPLICEAI_SCORES.items()]
    if records is not None:
        lines = ["\t".join(["21", position, ".", reference, alternate, ".", ".", "SpliceAI=" + ",".join(f"{alternate}|{entry}" for entry in entries)]) for entries in records]
    elif join:
        lines, _ = generator.in_record_forms([(f"21:{position}", line.split("\t")) for line in lines], generator.RecordForms(join_entries={f"21:{position}": None}))
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / ("snv.synthetic.vcf.gz" if len(reference) == len(alternate) == 1 else "indel.synthetic.vcf.gz")
    header = block(("\n".join(generator.SPLICEAI_HEADER) + "\n").encode())
    body = b"".join(block(f"{line}\n".encode()) for line in lines)
    path.write_bytes(header + body + block(b""))
    index = b"TBI\x01" + struct.pack("<8i", 1, 0, 1, 2, 2, ord("#"), 0, 3) + b"21\0"
    index += struct.pack("<iIi2Qi", 1, 0, 1, len(header) << 16, (len(header) + len(body)) << 16, 0)
    Path(f"{path}.tbi").write_bytes(gzip.compress(index))
    return path


def test_classify_names_the_spliceai_multi_entry_cell_class_under_the_score_files(tmp_path: Path):
    """A made-up score record of two entries, GENEA's then GENEB's: Ensembl VEP writes the
    record's whole value on the GENEA transcript's row and nothing on the GENEB transcript's,
    vep-rs each gene's own entry and a flag on both. Under `--spliceai-table` the three
    differing fields are `spliceai_multi_entry_cell`, on rows that state their gene as `SYMBOL`
    and on rows whose gene the corpus's cache names; without the files, under a file that
    scores each gene in a record of its own, and on rows whose gene nothing names, each is an
    `unexplained_residual`. So is a vep-rs value that is the other gene's entry, with the flag
    beside it, and a first-entry row on which Ensembl VEP writes no value."""
    corpus = tmp_path / "c"
    _manifest(corpus, [_record("21\t5000\tdel_two_genes\tGCA\tG\t.\t.\t.")])
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    intron = "intron_variant"

    def rows(path: Path, *values: str, symbols: bool = True) -> Path:
        path.write_text("".join(
            _feature_row("del_two_genes", "21:5001-5002", "-", feature, "Transcript", intron,
                         "IMPACT=MODIFIER" + (f";SYMBOL={gene}" if symbols else "") + (f";{value}" if value else ""))
            for (feature, gene), value in zip((("ENST1", "GENEA"), ("ENST2", "GENEB")), values)
        ), encoding="utf-8")
        return path

    vep = rows(tmp_path / "vep.txt", f"SpliceAI_pred={SPLICEAI_WHOLE};SpliceAI_cutoff=FAIL", "")
    rs = rows(tmp_path / "rs.txt", f"SpliceAI_pred={first};SpliceAI_cutoff=FAIL", f"SpliceAI_pred={second};SpliceAI_cutoff=PASS")
    joined = _spliceai_score_file(tmp_path / "joined")
    single = _spliceai_score_file(tmp_path / "single", join=False)

    def classes(rs_output: Path, *tables: Path, vep_output: Path = vep) -> dict:
        args = ["classify", "--corpus", str(corpus), "--vep-default", str(vep_output), "--vep-rs-output", str(rs_output),
                "--field", "SpliceAI_pred", "--field", "SpliceAI_cutoff", *(word for table in tables for word in ("--spliceai-table", str(table)))]
        assert bgc.main(args) == 0
        manifest = json.loads((corpus / "manifest.json").read_text())
        return {(d["feature"], d["field"]): d["expected_divergence"] for d in manifest["field_divergences"] if d["field"] in manifest["fields_compared"]}

    keys = (("ENST1", "SpliceAI_pred"), ("ENST2", "SpliceAI_pred"), ("ENST2", "SpliceAI_cutoff"))
    assert classes(rs) == dict.fromkeys(keys, "unexplained_residual")
    assert json.loads((corpus / "manifest.json").read_text())["field_divergence_summary"] == {"unexplained_residual": 3}
    assert classes(rs, joined) == dict.fromkeys(keys, "spliceai_multi_entry_cell")
    assert classes(rs, single) == dict.fromkeys(keys, "unexplained_residual")
    other = rows(tmp_path / "other.txt", f"SpliceAI_pred={first};SpliceAI_cutoff=FAIL", f"SpliceAI_pred={first};SpliceAI_cutoff=PASS")
    assert classes(other, joined) == {keys[0]: "spliceai_multi_entry_cell", keys[1]: "unexplained_residual", keys[2]: "unexplained_residual"}
    silent = rows(tmp_path / "silent.txt", "", "")
    assert classes(rs, joined, vep_output=silent) == {
        ("ENST1", "SpliceAI_pred"): "unexplained_residual", ("ENST1", "SpliceAI_cutoff"): "unexplained_residual",
        keys[1]: "spliceai_multi_entry_cell", keys[2]: "spliceai_multi_entry_cell",
    }
    # Rows that state no SYMBOL: nothing names their gene until the corpus holds a cache that does.
    bare_vep = rows(tmp_path / "bare-vep.txt", f"SpliceAI_pred={SPLICEAI_WHOLE};SpliceAI_cutoff=FAIL", "", symbols=False)
    bare_rs = rows(tmp_path / "bare-rs.txt", f"SpliceAI_pred={first};SpliceAI_cutoff=FAIL", f"SpliceAI_pred={second};SpliceAI_cutoff=PASS", symbols=False)
    assert classes(bare_rs, joined, vep_output=bare_vep) == dict.fromkeys(keys, "unexplained_residual")
    shards = corpus / "json_cache" / "transcripts" / "21"
    shards.mkdir(parents=True)
    (shards / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST1", "gene_symbol": "GENEA"}, {"stable_id": "ENST2", "gene_symbol": "GENEB"}]), encoding="utf-8")
    assert classes(bare_rs, joined, vep_output=bare_vep) == dict.fromkeys(keys, "spliceai_multi_entry_cell")


def test_classify_names_the_spliceai_multi_entry_cell_class_on_the_rows_of_a_single_base_record(tmp_path: Path):
    """A made-up score record of a substitution, one base for one base, of two entries, GENEA's
    then GENEB's, whose rows stand at the one position with the alternate base as their allele:
    Ensembl VEP writes the record's whole value on the GENEA transcript's row and nothing on the
    GENEB transcript's, vep-rs each gene's own entry. Under `--spliceai-table` both differing
    fields are `spliceai_multi_entry_cell`; without the files, and under a file that scores each
    gene in a record of its own, each is an `unexplained_residual`. So is a vep-rs value that is
    the other gene's entry, and an Ensembl VEP value that is the first entry joined with itself."""
    corpus = tmp_path / "c"
    _manifest(corpus, [_record("21\t7000\tsnv_two_genes\tA\tC\t.\t.\t.")])
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())

    def rows(path: Path, *values: str) -> Path:
        path.write_text("".join(
            _feature_row("snv_two_genes", "21:7000", "C", feature, "Transcript", "intron_variant",
                         f"IMPACT=MODIFIER;SYMBOL={gene}" + (f";SpliceAI_pred={value}" if value else ""))
            for (feature, gene), value in zip((("ENST1", "GENEA"), ("ENST2", "GENEB")), values)
        ), encoding="utf-8")
        return path

    vep = rows(tmp_path / "vep.txt", f"{first},C|{second}", "")
    rs = rows(tmp_path / "rs.txt", first, second)
    joined = _spliceai_score_file(tmp_path / "joined", variant=("7000", "A", "C"))
    single = _spliceai_score_file(tmp_path / "single", join=False, variant=("7000", "A", "C"))

    def classes(rs_output: Path, *tables: Path, vep_output: Path = vep) -> dict:
        args = ["classify", "--corpus", str(corpus), "--vep-default", str(vep_output), "--vep-rs-output", str(rs_output),
                "--field", "SpliceAI_pred", *(word for table in tables for word in ("--spliceai-table", str(table)))]
        assert bgc.main(args) == 0
        manifest = json.loads((corpus / "manifest.json").read_text())
        return {d["feature"]: d["expected_divergence"] for d in manifest["field_divergences"] if d["field"] in manifest["fields_compared"]}

    assert classes(rs) == {"ENST1": "unexplained_residual", "ENST2": "unexplained_residual"}
    assert classes(rs, joined) == {"ENST1": "spliceai_multi_entry_cell", "ENST2": "spliceai_multi_entry_cell"}
    assert classes(rs, single) == {"ENST1": "unexplained_residual", "ENST2": "unexplained_residual"}
    other = rows(tmp_path / "other.txt", first, first)
    assert classes(other, joined) == {"ENST1": "spliceai_multi_entry_cell", "ENST2": "unexplained_residual"}
    twice = rows(tmp_path / "twice.txt", f"{first},C|{first}", "")
    assert classes(rs, joined, vep_output=twice) == {"ENST1": "unexplained_residual", "ENST2": "spliceai_multi_entry_cell"}


def test_classify_names_the_spliceai_multi_entry_cell_class_under_every_record_line_at_the_allele(tmp_path: Path):
    """A made-up score file whose record line of two entries, GENEA's then GENEB's, stands more
    than once at the substitution's allele, or beside a line of the same entries in the other
    order or of GENEA's entry alone: Ensembl VEP writes a row the whole value of the last line
    its gene leads, or nothing where it leads none, vep-rs the gene's entry. Under
    `--spliceai-table` every `SpliceAI_pred` the two fill differently is
    `spliceai_multi_entry_cell`; where GENEA's line of one entry stands last both write that
    entry on its row, which differs in nothing. Where the lines give a gene entries that
    differ in text, vep-rs writes the entry of the last line that names the gene, and its
    rows are `spliceai_multi_entry_cell` under two and three lines of the same genes with
    other scores and under two lines of other scores that each gene leads one of. The row of
    a gene no entry names, and a row on which vep-rs's value is an earlier line's entry, are
    each an `unexplained_residual`, and so is every differing row without the files."""
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    changed = "GENEA|0.31|0.00|0.44|0.09|3|-2|11|-30"
    other, third, last = f"GENEB|{SPLICEAI_SCORES['GENEA']}", f"GENEA|{SPLICEAI_SCORES['GENEB']}", "GENEB|0.31|0.00|0.44|0.09|3|-2|11|-30"
    both = (first, second)
    taken, left = "spliceai_multi_entry_cell", "unexplained_residual"
    genes = (("ENST1", "GENEA"), ("ENST2", "GENEB"), ("ENST3", "GENEX"))

    def whole(*values: str) -> str:
        return ",C|".join(values)

    for name, lines, vep_values, rs_values, classes in (
        ("one-line", [both], (whole(*both), ""), both, {"ENST1": taken, "ENST2": taken}),
        ("twice", [both, both], (whole(*both), ""), both, {"ENST1": taken, "ENST2": taken}),
        ("three-times", [both, both, both], (whole(*both), ""), both, {"ENST1": taken, "ENST2": taken}),
        ("other-order", [both, both[::-1]], (whole(*both), whole(*both[::-1])), both, {"ENST1": taken, "ENST2": taken}),
        ("one-entry-last", [both, (first,)], (first, ""), both, {"ENST2": taken}),
        ("one-entry-first", [(first,), both], (whole(*both), ""), both, {"ENST1": taken, "ENST2": taken}),
        ("gene-unnamed", [both, both], (whole(*both), "", ""), (*both, second), {"ENST1": taken, "ENST2": taken, "ENST3": left}),
        ("other-scores", [both, (changed, other)], (whole(changed, other), ""), (changed, other), {"ENST1": taken, "ENST2": taken}),
        ("other-scores-three-times", [both, (changed, other), (third, last)], (whole(third, last), ""), (third, last), {"ENST1": taken, "ENST2": taken}),
        ("other-scores-other-order", [both, (other, changed)], (whole(*both), whole(other, changed)), (changed, other), {"ENST1": taken, "ENST2": taken}),
        ("earlier-line", [both, (changed, other)], (whole(changed, other), ""), (first, other), {"ENST1": left, "ENST2": taken}),
    ):
        directory = tmp_path / name
        directory.mkdir()
        corpus = directory / "c"
        _manifest(corpus, [_record("21\t7000\tsnv_two_genes\tA\tC\t.\t.\t.")])
        outputs = []
        for side, values in (("vep", vep_values), ("rs", rs_values)):
            outputs.append(directory / f"{side}.txt")
            outputs[-1].write_text("".join(
                _feature_row("snv_two_genes", "21:7000", "C", feature, "Transcript", "intron_variant",
                             f"IMPACT=MODIFIER;SYMBOL={gene}" + (f";SpliceAI_pred={value}" if value else ""))
                for (feature, gene), value in zip(genes, values)
            ), encoding="utf-8")
        table = _spliceai_score_file(directory / "scores", variant=("7000", "A", "C"), records=lines)
        for files, expected in ((["--spliceai-table", str(table)], classes), ([], dict.fromkeys(classes, left))):
            args = ["classify", "--corpus", str(corpus), "--vep-default", str(outputs[0]), "--vep-rs-output", str(outputs[1]), "--field", "SpliceAI_pred", *files]
            assert bgc.main(args) == 0
            manifest = json.loads((corpus / "manifest.json").read_text())
            found = {d["feature"]: d["expected_divergence"] for d in manifest["field_divergences"] if d["field"] in manifest["fields_compared"]}
            assert found == expected, (name, files)


def _spliceai_derived_entries(tmp_path: Path, vep_values: tuple[str, str], rs_values: tuple[str, str], fields: tuple[str, ...]) -> dict:
    """The entries `classify` writes beside the compared fields, by (Feature, field), for the
    deletion's two rows under the score file of the record of two entries: each (vep, vep_rs)."""
    corpus = tmp_path / "c"
    _manifest(corpus, [_record("21\t5000\tdel_two_genes\tGCA\tG\t.\t.\t.")])

    def rows(path: Path, values: tuple[str, str]) -> Path:
        path.write_text("".join(
            _feature_row("del_two_genes", "21:5001-5002", "-", feature, "Transcript", "intron_variant", f"IMPACT=MODIFIER;SYMBOL={gene}" + (f";{value}" if value else ""))
            for (feature, gene), value in zip((("ENST1", "GENEA"), ("ENST2", "GENEB")), values)
        ), encoding="utf-8")
        return path

    args = ["classify", "--corpus", str(corpus), "--vep-default", str(rows(tmp_path / "vep.txt", vep_values)),
            "--vep-rs-output", str(rows(tmp_path / "rs.txt", rs_values)), "--spliceai-table", str(_spliceai_score_file(tmp_path / "joined")),
            *(word for field in fields for word in ("--field", field))]
    assert bgc.main(args) == 0
    manifest = json.loads((corpus / "manifest.json").read_text())
    assert {d["expected_divergence"] for d in manifest["field_divergences"]} == {"spliceai_multi_entry_cell"}
    return {(d["feature"], d["field"]): (d["vep"], d["vep_rs"]) for d in manifest["field_divergences"] if d["field"] not in fields}


def test_classify_writes_the_member_fields_of_a_row_of_the_spliceai_class_for_the_vcf_format(tmp_path: Path):
    """The VCF format writes a SpliceAI value as nine fields whatever the run's options, and the
    golden harness reads a documented field by its format's name. For a row written under
    `spliceai_multi_entry_cell` in a run whose default format writes the one field, `classify`
    also writes each member field whose two parts differ: the last position alone on the first
    entry's gene's row, where Ensembl VEP's ninth part runs on into the second entry's allele,
    and all nine on the second entry's gene's row. A run that compares the nine fields itself
    gains none."""
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    derived = _spliceai_derived_entries(tmp_path, (f"SpliceAI_pred={SPLICEAI_WHOLE}", ""), (f"SpliceAI_pred={first}", f"SpliceAI_pred={second}"), ("SpliceAI_pred",))
    members = {key: value for key, value in derived.items() if key[1] != "spliceai"}
    assert members == {
        ("ENST1", "SpliceAI_pred_DP_DL"): ("-9,G", "-9"),
        **{("ENST2", field): ("-", part) for field, part in zip(bgc.SPLICEAI_MEMBER_FIELDS, second.split("|"))},
    }

    def pairs(value: str) -> str:
        return ";".join(f"{field}={part}" for field, part in zip(bgc.SPLICEAI_MEMBER_FIELDS, value.split("|")))

    split = tmp_path / "split"
    split.mkdir()
    derived = _spliceai_derived_entries(split, (pairs(SPLICEAI_WHOLE), ""), (pairs(first), pairs(second)), bgc.SPLICEAI_MEMBER_FIELDS)
    assert set(derived) == {("ENST1", "spliceai"), ("ENST2", "spliceai")}


def test_classify_writes_the_json_object_of_a_row_of_the_spliceai_class_as_the_harness_prints_one(tmp_path: Path):
    """The JSON format writes a SpliceAI value as one object, which the golden harness prints
    with its members in the byte order of their names and no spaces. For a row written under
    `spliceai_multi_entry_cell`, `classify` also writes that object for both engines: a number
    for each member whose text is one, written by its value as both engines write it (`0.00`
    as `0`, `1.00` as `1`, `0.10` as `0.1`), Ensembl VEP's last position a string where it runs
    on into the second entry's allele, the cutoff flag among the members of a row that carries
    one, and no object for a row without the value."""
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    derived = _spliceai_derived_entries(
        tmp_path, (f"SpliceAI_pred={SPLICEAI_WHOLE};SpliceAI_cutoff=FAIL", ""),
        (f"SpliceAI_pred={first};SpliceAI_cutoff=FAIL", f"SpliceAI_pred={second};SpliceAI_cutoff=PASS"), ("SpliceAI_pred", "SpliceAI_cutoff"),
    )
    own = '{"DP_AG":-6,"DP_AL":-41,"DP_DG":-17,"DP_DL":-9,"DS_AG":0.15,"DS_AL":0.08,"DS_DG":0.06,"DS_DL":0.17,"SYMBOL":"GENEA","SpliceAI_cutoff":"FAIL"}'
    assert derived[("ENST1", "spliceai")] == (own.replace('"DP_DL":-9', '"DP_DL":"-9,G"'), own)
    assert derived[("ENST2", "spliceai")] == (
        "-", '{"DP_AG":-25,"DP_AL":-23,"DP_DG":37,"DP_DL":17,"DS_AG":0,"DS_AL":1,"DS_DG":0.1,"DS_DL":0.05,"SYMBOL":"GENEB","SpliceAI_cutoff":"PASS"}',
    )


def test_classify_names_the_edited_model_hgvsc_class(tmp_path: Path):
    """An `HGVSc` whose two notations differ in their cDNA coordinates alone, by one of the edit
    distances of a transcript the corpus cache files with an `_rna_edit` that inserts or deletes
    bases, is `refseq_edited_model_hgvsc_offset`: the exon coordinate moved by the 3-base
    deletion, the exon end moved by the two insertions' running sum, and the exon end a 9-base
    deletion moves by nine while the start slides into the exon. The same shape on a model whose edit substitutes
    one base for one, on a model without an edit, with another reference name or allele text,
    or a coordinate off by one on a model edited by three stays an `unexplained_residual`. The
    cache's indel-edited models are read as Ensembl VEP reads the attribute value: a two-field
    value deletes, a three-field value whose bases span its coordinates substitutes, any other
    inserts; a model's distances are the running sums of its edits in start order."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\trs_del\tAT\tA\t.\tPASS\t."),
        _record("21\t200\trs_sub\tA\tG\t.\tPASS\t."),
        _record("21\t300\trs_plain\tA\tG\t.\tPASS\t."),
        _record("21\t400\trs_name\tA\tG\t.\tPASS\t."),
        _record("21\t500\trs_allele\tA\tG\t.\tPASS\t."),
        _record("21\t600\trs_two\tA\tG\t.\tPASS\t."),
        _record("21\t700\trs_off\tA\tG\t.\tPASS\t."),
        _record("21\t800\trs_slide\tAT\tA\t.\tPASS\t."),
    ])
    (corpus / "json_cache" / "transcripts" / "21").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([
        {"stable_id": "NM_1.1", "start": 50, "end": 900, "strand": 1, "attributes": [{"code": "_rna_edit", "value": "50 52 "}]},
        {"stable_id": "NM_2.1", "start": 50, "end": 900, "strand": 1, "attributes": [{"code": "_rna_edit", "value": "70 70 T"}]},
        {"stable_id": "NM_3.1", "start": 50, "end": 900, "strand": 1, "attributes": []},
        {"stable_id": "NM_4.1", "start": 50, "end": 900, "strand": 1, "attributes": [{"code": "_rna_edit", "value": "10 9 AAA"}]},
        {"stable_id": "NM_5.1", "start": 50, "end": 900, "strand": 1, "attributes": [{"code": "_rna_edit", "value": "10 9 AAA"}]},
        {"stable_id": "NM_6.1", "start": 50, "end": 900, "strand": 1, "attributes": [
            {"code": "_rna_edit", "value": "900 899 AAAAAAAAAAAAAAAAAA"}, {"code": "_rna_edit", "value": "20 19 G"}]},
        {"stable_id": "NM_7.1", "start": 50, "end": 900, "strand": 1, "attributes": [{"code": "_rna_edit", "value": "10 9 AAA"}]},
        {"stable_id": "NM_8.1", "start": 50, "end": 900, "strand": 1, "attributes": [{"code": "_rna_edit", "value": "3387 3395 "}]},
    ]), encoding="utf-8")
    assert bgc.indel_edited_models(corpus / "json_cache") == {
        "NM_1.1": frozenset({3}), "NM_4.1": frozenset({3}), "NM_5.1": frozenset({3}),
        "NM_6.1": frozenset({1, 19}), "NM_7.1": frozenset({3}), "NM_8.1": frozenset({9}),
    }
    assert bgc.indel_edited_models(tmp_path / "absent") == {}
    hgvsc = "IMPACT=MODIFIER;STRAND=1;HGVSc={}"
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(
        _row("rs_del", "21:101", "-", "NM_1.1", "intron_variant", hgvsc.format("NM_1.1:c.247+601del"))
        + _row("rs_sub", "21:200", "G", "NM_2.1", "intron_variant", hgvsc.format("NM_2.1:c.247+601T>A"))
        + _row("rs_plain", "21:300", "G", "NM_3.1", "intron_variant", hgvsc.format("NM_3.1:c.247+601T>A"))
        + _row("rs_name", "21:400", "G", "NM_4.1", "intron_variant", hgvsc.format("NM_4.1:c.247+601T>A"))
        + _row("rs_allele", "21:500", "G", "NM_5.1", "intron_variant", hgvsc.format("NM_5.1:c.247+601T>A"))
        + _row("rs_two", "21:600", "G", "NM_6.1", "intron_variant", hgvsc.format("NM_6.1:c.56+4A>G"))
        + _row("rs_off", "21:700", "G", "NM_7.1", "intron_variant", hgvsc.format("NM_7.1:c.5A>G"))
        + _row("rs_slide", "21:801", "-", "NM_8.1", "intron_variant", hgvsc.format("NM_8.1:c.3546_3553+23del")),
        encoding="utf-8",
    )
    rs.write_text(
        _row("rs_del", "21:101", "-", "NM_1.1", "intron_variant", hgvsc.format("NM_1.1:c.250+601del"))
        + _row("rs_sub", "21:200", "G", "NM_2.1", "intron_variant", hgvsc.format("NM_2.1:c.250+601T>A"))
        + _row("rs_plain", "21:300", "G", "NM_3.1", "intron_variant", hgvsc.format("NM_3.1:c.250+601T>A"))
        + _row("rs_name", "21:400", "G", "NM_4.1", "intron_variant", hgvsc.format("NM_4.2:c.250+601T>A"))
        + _row("rs_allele", "21:500", "G", "NM_5.1", "intron_variant", hgvsc.format("NM_5.1:c.250+601T>C"))
        + _row("rs_two", "21:600", "G", "NM_6.1", "intron_variant", hgvsc.format("NM_6.1:c.57+4A>G"))
        + _row("rs_off", "21:700", "G", "NM_7.1", "intron_variant", hgvsc.format("NM_7.1:c.6A>G"))
        + _row("rs_slide", "21:801", "-", "NM_8.1", "intron_variant", hgvsc.format("NM_8.1:c.3544+2_3544+32del")),
        encoding="utf-8",
    )
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "HGVSc"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert {(d["location"], d["feature"]): d["expected_divergence"] for d in m["field_divergences"]} == {
        ("21:101", "NM_1.1"): "refseq_edited_model_hgvsc_offset", ("21:200", "NM_2.1"): "unexplained_residual",
        ("21:300", "NM_3.1"): "unexplained_residual", ("21:400", "NM_4.1"): "unexplained_residual",
        ("21:500", "NM_5.1"): "unexplained_residual", ("21:600", "NM_6.1"): "refseq_edited_model_hgvsc_offset",
        ("21:700", "NM_7.1"): "unexplained_residual", ("21:801", "NM_8.1"): "refseq_edited_model_hgvsc_offset",
    }
    assert m["field_divergence_summary"] == {"refseq_edited_model_hgvsc_offset": 3, "unexplained_residual": 5}


def test_classify_names_the_edited_model_hgvs_offset_class(tmp_path: Path):
    """An `HGVS_OFFSET` Ensembl VEP leaves empty and vep-rs fills with a non-zero integer on a
    transcript the corpus cache files with an `_rna_edit` that inserts or deletes bases is
    `refseq_edited_model_hgvs_offset`, a reverse-strand model's negative offset included; the same
    offset on a model without such an edit, a zero, a non-integer, and two filled offsets that
    differ stay an `unexplained_residual`."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\trs_edited\tT\tTGCG\t.\tPASS\t."),
        _record("21\t200\trs_plain\tT\tTGCG\t.\tPASS\t."),
        _record("21\t300\trs_zero\tT\tTGCG\t.\tPASS\t."),
        _record("21\t400\trs_text\tT\tTGCG\t.\tPASS\t."),
        _record("21\t500\trs_both\tT\tTGCG\t.\tPASS\t."),
        _record("21\t600\trs_minus\tT\tTGCG\t.\tPASS\t."),
    ])
    (corpus / "json_cache" / "transcripts" / "21").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([
        {"stable_id": "NM_1.1", "start": 50, "end": 900, "strand": 1, "attributes": [{"code": "_rna_edit", "value": "160 174 "}]},
        {"stable_id": "NM_2.1", "start": 50, "end": 900, "strand": 1, "attributes": [{"code": "_rna_edit", "value": "70 70 T"}]},
    ]), encoding="utf-8")
    extra = "IMPACT=MODERATE;STRAND=1;BAM_EDIT=OK;HGVSc=NM_1.1:c.986_988dup"
    rows = [("rs_edited", "21:100-101", "NM_1.1"), ("rs_plain", "21:200-201", "NM_2.1"), ("rs_zero", "21:300-301", "NM_1.1"),
            ("rs_text", "21:400-401", "NM_1.1"), ("rs_both", "21:500-501", "NM_1.1"), ("rs_minus", "21:600-601", "NM_1.1")]
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(
        _row(name, loc, "GCG", feature, "inframe_insertion", extra + (";HGVS_OFFSET=3" if name == "rs_both" else ""))
        for name, loc, feature in rows), encoding="utf-8")
    offsets = {"rs_edited": "21", "rs_plain": "21", "rs_zero": "0", "rs_text": "dup", "rs_both": "21", "rs_minus": "-7"}
    rs.write_text("".join(
        _row(name, loc, "GCG", feature, "inframe_insertion", f"{extra};HGVS_OFFSET={offsets[name]}")
        for name, loc, feature in rows), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert {(d["location"], d["feature"]): d["expected_divergence"] for d in m["field_divergences"]} == {
        ("21:100-101", "NM_1.1"): "refseq_edited_model_hgvs_offset", ("21:200-201", "NM_2.1"): "unexplained_residual",
        ("21:300-301", "NM_1.1"): "unexplained_residual", ("21:400-401", "NM_1.1"): "unexplained_residual",
        ("21:500-501", "NM_1.1"): "unexplained_residual", ("21:600-601", "NM_1.1"): "refseq_edited_model_hgvs_offset",
    }
    assert m["field_divergence_summary"] == {"refseq_edited_model_hgvs_offset": 2, "unexplained_residual": 4}


def test_cdna_coordinate_order_walks_the_transcript():
    """The 5' UTR counts down to c.-1, the coding sequence up, the 3' UTR up from c.*1; an
    intron offset follows its exon base, `-M` before it and `+M` after; `*-M`, the stop codon's
    last base plus a negative intron distance printed without that base's coordinate, sits
    after every coding coordinate and before `*1`."""
    coordinates = ["-5-2", "-5", "-5+2", "-1", "1", "3-1", "3", "3+1", "300", "*-3", "*1", "*12-3", "*12"]
    assert sorted(coordinates, key=bgc.cdna_coordinate_order) == coordinates
    assert bgc.cdna_coordinate_order("*-3") == (1, 0, -3)
    assert bgc.cdna_coordinate_order("-228+19") == (-1, -228, 19)


def test_classify_names_the_shifted_indel_hgvsc_classes(tmp_path: Path):
    """An `HGVSc`, `HGVS_OFFSET` or HGVSIntronOffset field Ensembl VEP leaves empty where vep-rs
    writes an insertion's or deletion's notation beside a non-zero offset is
    `hgvs_shifted_indel_written_span_gate`; two deletions of one transcript ending together, the
    reference's starting later, on a row both engines give a non-zero offset, are
    `hgvs_shifted_deletion_clamped_span`, across the 5' UTR, coding, intronic and 3' UTR
    numberings and down to a one-base reference span. A substitution's notation the reference
    lacks, a notation without an offset, another transcript's notation, a row annotated upstream,
    two deletions ending apart, a candidate starting later, a clamped pair without an offset and
    a clamped pair on an insertion's row stay an `unexplained_residual`."""
    corpus = tmp_path / "c"
    names = ["rs_ins", "rs_del", "rs_last", "rs_sub", "rs_nooff", "rs_clamp", "rs_utr", "rs_intron", "rs_one", "rs_end", "rs_later", "rs_cold",
             "rs_other", "rs_outside", "rs_insdel"]
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\t{name}\tT\tTA\t.\tPASS\t.") for i, name in enumerate(names)])
    base = "IMPACT=MODIFIER;STRAND=1"
    vep_extra = {
        "rs_ins": base, "rs_del": base, "rs_last": base, "rs_sub": base, "rs_nooff": base,
        "rs_clamp": f"{base};HGVSc=ENST1.1:n.34_43del;HGVS_OFFSET=33",
        "rs_utr": f"{base};HGVSc=ENST1.1:c.-355_-327del;HGVS_OFFSET=6",
        "rs_intron": f"{base};HGVSc=ENST1.1:c.-301_-228+19del;HGVS_OFFSET=6",
        "rs_one": f"{base};HGVSc=ENST1.1:n.14del;HGVS_OFFSET=13",
        "rs_end": f"{base};HGVSc=ENST1.1:n.34_43del;HGVS_OFFSET=33",
        "rs_later": f"{base};HGVSc=ENST1.1:n.34_43del;HGVS_OFFSET=33",
        "rs_cold": f"{base};HGVSc=ENST1.1:n.34_43del",
        "rs_other": base, "rs_outside": base,
        "rs_insdel": f"{base};HGVSc=ENST1.1:n.34_43del;HGVS_OFFSET=33",
    }
    rs_extra = {
        "rs_ins": f"{base};HGVSc=ENST1.1:n.15dup;HGVS_OFFSET=17",
        "rs_del": f"{base};HGVSc=ENST1.1:n.5_25del;HGVS_OFFSET=33",
        "rs_last": f"{base};HGVSc=ENST1.1:c.*2333_*2540del;HGVS_OFFSET=74;HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=0",
        "rs_sub": f"{base};HGVSc=ENST1.1:n.16A>T;HGVS_OFFSET=17",
        "rs_nooff": f"{base};HGVSc=ENST1.1:n.15dup",
        "rs_clamp": f"{base};HGVSc=ENST1.1:n.23_43del;HGVS_OFFSET=33",
        "rs_utr": f"{base};HGVSc=ENST1.1:c.-361_-327del;HGVS_OFFSET=6",
        "rs_intron": f"{base};HGVSc=ENST1.1:c.-307_-228+19del;HGVS_OFFSET=6",
        "rs_one": f"{base};HGVSc=ENST1.1:n.13_14del;HGVS_OFFSET=13",
        "rs_end": f"{base};HGVSc=ENST1.1:n.23_42del;HGVS_OFFSET=33",
        "rs_later": f"{base};HGVSc=ENST1.1:n.35_43del;HGVS_OFFSET=33",
        "rs_cold": f"{base};HGVSc=ENST1.1:n.23_43del",
        "rs_other": f"{base};HGVSc=ENST2.1:n.15dup;HGVS_OFFSET=17",
        "rs_outside": f"{base};HGVSc=ENST1.1:n.15dup;HGVS_OFFSET=17",
        "rs_insdel": f"{base};HGVSc=ENST1.1:n.23_43del;HGVS_OFFSET=33",
    }
    allele = {"rs_ins": "A", "rs_sub": "A", "rs_nooff": "A", "rs_other": "A", "rs_outside": "A", "rs_insdel": "A"}
    consequence = {"rs_outside": "upstream_gene_variant"}
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    for path, extras in ((vep, vep_extra), (rs, rs_extra)):
        path.write_text("".join(
            _row(name, f"21:{100 * (i + 1)}-{100 * (i + 1) + 1}", allele.get(name, "-"), "ENST1",
                 consequence.get(name, "non_coding_transcript_exon_variant"), extras[name])
            for i, name in enumerate(names)), encoding="utf-8")
    fields = ["--field", "HGVSc", "--field", "HGVS_OFFSET", "--field", "HGVS_IntronStartOffset", "--field", "HGVS_IntronEndOffset"]
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs), *fields]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {(d["location"], d["field"]): d["expected_divergence"] for d in m["field_divergences"]}
    gate, clamped = "hgvs_shifted_indel_written_span_gate", "hgvs_shifted_deletion_clamped_span"
    assert classes == {
        ("21:100-101", "HGVSc"): gate, ("21:100-101", "HGVS_OFFSET"): gate,
        ("21:200-201", "HGVSc"): gate, ("21:200-201", "HGVS_OFFSET"): gate,
        ("21:300-301", "HGVSc"): gate, ("21:300-301", "HGVS_OFFSET"): gate,
        ("21:300-301", "HGVS_IntronStartOffset"): gate, ("21:300-301", "HGVS_IntronEndOffset"): gate,
        ("21:400-401", "HGVSc"): "unexplained_residual", ("21:400-401", "HGVS_OFFSET"): "unexplained_residual",
        ("21:500-501", "HGVSc"): "unexplained_residual",
        ("21:600-601", "HGVSc"): clamped, ("21:700-701", "HGVSc"): clamped, ("21:800-801", "HGVSc"): clamped,
        ("21:900-901", "HGVSc"): clamped,
        ("21:1000-1001", "HGVSc"): "unexplained_residual", ("21:1100-1101", "HGVSc"): "unexplained_residual",
        ("21:1200-1201", "HGVSc"): "unexplained_residual",
        ("21:1300-1301", "HGVSc"): "unexplained_residual", ("21:1300-1301", "HGVS_OFFSET"): "unexplained_residual",
        ("21:1400-1401", "HGVSc"): "unexplained_residual", ("21:1400-1401", "HGVS_OFFSET"): "unexplained_residual",
        ("21:1500-1501", "HGVSc"): "unexplained_residual",
    }
    assert m["field_divergence_summary"] == {gate: 8, clamped: 4, "unexplained_residual": 11}


def test_classify_names_the_clipped_reference_against_whole_alternate_class(tmp_path: Path):
    """An `HGVSc` on a sequence allele's row whose reference cDNA position has an open side,
    Ensembl VEP's a notation of the row's transcript and vep-rs's empty or a notation of the same
    transcript and numbering whose every coordinate lies within the row's span plus one of every
    coordinate of Ensembl VEP's, is `hgvs_clipped_reference_against_whole_alternate`: the pairs
    the committed corpora, the examples and the class page's multi-allelic records print, of
    every kind the whole alternate reads as against the bases inside, and a `del` candidate
    against the empty notation Ensembl VEP writes where the whole alternate equals the bases
    inside. An allele shorter than its span that deletes bases over the edge has no vep-rs
    notation, and Ensembl VEP's `del`, `delins` or substitution against it is the class. A
    closed cDNA position, a deletion's row, another transcript or numbering, a coordinate
    beyond the reach, an empty reference against a non-deletion, and an empty candidate
    against a reference inversion or, on the row of an allele as long as its span, against
    anything but a reference `ins` or `dup` stay an `unexplained_residual`."""
    corpus = tmp_path / "c"
    cases = {
        # name: (allele, span, cDNA position, VEP HGVSc, vep-rs HGVSc, class)
        "rs_delins": ("TGCTGGGCGCC", 11, "?-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST1.6:c.-47A>C", True),
        "rs_flank": ("CGGCATGCTGGGAGCC", 16, "?-6", "ENST1.6:c.-49-1_-49insCGGCATGCTG", "-", True),
        "rs_last": ("TGGCTTATAAAGTGAA", 16, "8485-?", "ENST1.6:c.*1461_*1465delinsCTTATAAAGTGAA", "ENST1.6:c.*1461G>C", True),
        "rs_inv": ("ACCC", 4, "?-2", "ENST1.6:c.-49_-48delinsACCC", "ENST1.6:c.-49_-48inv", True),
        "rs_sub_ins": ("GGA", 2, "?-1", "ENST1.6:c.-49_-48insGA", "ENST1.6:c.-49_-48insA", True),
        "rs_sub_del": ("GGG", 4, "?-3", "ENST1.6:c.-47A>G", "ENST1.6:c.-47del", True),
        "rs_one_ins": ("CA", 2, "147-?", "ENST1.6:c.146_147insC", "ENST1.6:c.147A>C", True),
        "rs_ins_past": ("GGGAC", 4, "?-3", "ENST1.6:c.-47delinsGAC", "ENST1.6:c.-47_-46insC", True),
        "rs_anchor_del": ("T", 4, "?-2", "ENST1.6:c.-49_-48delinsT", "-", True),
        "rs_sub_two_del": ("TG", 4, "?-2", "ENST1.6:c.-49G>T", "ENST1.6:c.-49_-48del", True),
        "rs_sub_delins": ("GA", 3, "?-2", "ENST1.6:c.-48G>A", "ENST1.6:c.-49_-48delinsA", True),
        "rs_delins_wider": ("GAAA", 4, "?-2", "ENST1.6:c.-48delinsAAA", "ENST1.6:c.-49_-48delinsAA", True),
        "rs_ins_delins": ("GGAA", 4, "?-2", "ENST1.6:c.-48_-47insAA", "ENST1.6:c.-49_-48delinsAA", True),
        "rs_wrap_del": ("TGG", 4, "?-2", "ENST1.6:c.-49-1_-49insT", "ENST1.6:c.-48del", True),
        "rs_wrap_ins": ("TGGGG", 4, "?-2", "ENST1.6:c.-49-1_-49insTGG", "ENST1.6:c.-48_-47insG", True),
        "rs_dup_ins": ("GGGG", 3, "?-2", "ENST1.6:c.-49_-48dup", "ENST1.6:c.-48_-47insG", True),
        "rs_ins_dup": ("GGGGG", 3, "?-2", "ENST1.6:c.-48_-47insGGG", "ENST1.6:c.-49_-48dup", True),
        "rs_silent_del": ("GG", 3, "?-2", "-", "ENST1.6:c.-48del", True),
        "rs_silent_span": ("GG", 3, "?-2", "-", "ENST1.6:c.-49_-48del", True),
        "rs_silent_sub": ("GG", 3, "?-2", "-", "ENST1.6:c.-48A>C", False),
        "rs_silent_ins": ("GG", 3, "?-2", "-", "ENST1.6:c.-49_-48insA", False),
        "rs_shorter_del": ("G", 3, "?-2", "ENST1.6:c.-48del", "ENST1.6:c.-49_-48del", True),
        "rs_del_far": ("G", 3, "?-2", "ENST1.6:c.-48del", "ENST1.6:c.-43del", False),
        "rs_dup_outside": ("GG", 4, "?-1", "ENST1.6:c.-49dup", "-", True),
        "rs_refdel_over_edge": ("GG", 4, "?-1", "ENST1.6:c.-49del", "-", True),
        "rs_refinv_over_edge": ("CC", 4, "?-2", "ENST1.6:c.-49_-48inv", "-", False),
        "rs_refsub_over_edge": ("GA", 4, "?-2", "ENST1.6:c.-49G>A", "-", True),
        "rs_refdel_whole_span": ("GG", 2, "?-1", "ENST1.6:c.-49del", "-", False),
        "rs_refsub_whole_span": ("GA", 2, "?-2", "ENST1.6:c.-49G>A", "-", False),
        "rs_refdelins_closed": ("T", 4, "1-2", "ENST1.6:c.-49_-48delinsT", "-", False),
        "rs_closed": ("TGCTGGGCGCC", 11, "1-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST1.6:c.-47A>C", False),
        "rs_deletion": ("-", 11, "?-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST1.6:c.-47A>C", False),
        "rs_other": ("TGCTGGGCGCC", 11, "?-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST2.6:c.-47A>C", False),
        "rs_numbering": ("TGCTGGGCGCC", 11, "?-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST1.6:n.-47A>C", False),
        "rs_reach": ("TGCTGGGCGCC", 11, "?-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST1.6:c.-37A>C", True),
        "rs_far": ("TGCTGGGCGCC", 11, "?-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST1.6:c.-36A>C", False),
        "rs_utr3": ("TGCTGGGCGCC", 11, "?-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST1.6:c.*5A>C", False),
        "rs_empty": ("TGCTGGGCGCC", 11, "?-6", "-", "ENST1.6:c.-47A>C", False),
    }
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\t{name}\tT\tTA\t.\tPASS\t.") for i, name in enumerate(cases)])

    def location(i: int, span: int) -> str:
        return f"21:{100 * (i + 1)}-{100 * (i + 1) + span - 1}"

    def row(name: str, i: int, hgvsc: str) -> str:
        allele, span, cdna, _, _, _ = cases[name]
        extra = "IMPACT=MODIFIER;STRAND=1" + ("" if hgvsc == "-" else f";HGVSc={hgvsc}")
        return "\t".join([name, location(i, span), allele, "ENSG1", "ENST1", "Transcript", "5_prime_UTR_variant", cdna, "-", "-", "-", "-", "-", extra]) + "\n"

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, i, case[3]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    rs.write_text("".join(row(name, i, case[4]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs), "--field", "HGVSc"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {d["location"]: d["expected_divergence"] for d in m["field_divergences"]}
    klass = "hgvs_clipped_reference_against_whole_alternate"
    assert classes == {
        location(i, case[1]): (klass if case[5] else "unexplained_residual") for i, case in enumerate(cases.values())
    }
    assert m["field_divergence_summary"] == {klass: 24, "unexplained_residual": 14}
    assert {d["vep_rs"] for d in m["field_divergences"] if d["expected_divergence"] == klass} == {
        "ENST1.6:c.-47A>C", "-", "ENST1.6:c.*1461G>C", "ENST1.6:c.-49_-48inv", "ENST1.6:c.-49_-48insA", "ENST1.6:c.-47del", "ENST1.6:c.147A>C",
        "ENST1.6:c.-47_-46insC", "ENST1.6:c.-48del", "ENST1.6:c.-49_-48del", "ENST1.6:c.-49_-48delinsA", "ENST1.6:c.-49_-48delinsAA",
        "ENST1.6:c.-48_-47insG", "ENST1.6:c.-49_-48dup", "ENST1.6:c.-37A>C"}
    assert sum(1 for d in m["field_divergences"] if d["expected_divergence"] == klass and d["vep_rs"] == "-") == 5


def test_classify_clipped_reference_against_whole_alternate_accepts_a_masked_reference_base(tmp_path: Path):
    """A reference FASTA masked (`N`) over the row's window gives either engine's substitution an
    `N` reference base: vep-rs's `c.-47N>C` against Ensembl VEP's `delins`, and Ensembl VEP's
    `c.-47N>G` against vep-rs's `del`, are `hgvs_clipped_reference_against_whole_alternate` as the
    `A` and `G` reference bases are. An `N` alternate base is not a notation of the class, and the
    pair stays an `unexplained_residual`."""
    corpus = tmp_path / "c"
    cases = {
        # name: (allele, span, cDNA position, VEP HGVSc, vep-rs HGVSc, class)
        "rs_masked": ("TGCTGGGCGCC", 11, "?-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST1.6:c.-47N>C", True),
        "rs_masked_ref": ("GGG", 4, "?-3", "ENST1.6:c.-47N>G", "ENST1.6:c.-47del", True),
        "rs_masked_alt": ("TGCTGGGCGCC", 11, "?-6", "ENST1.6:c.-49_-47delinsTGCTGGGC", "ENST1.6:c.-47A>N", False),
    }
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\t{name}\tT\tTA\t.\tPASS\t.") for i, name in enumerate(cases)])

    def location(i: int, span: int) -> str:
        return f"21:{100 * (i + 1)}-{100 * (i + 1) + span - 1}"

    def row(name: str, i: int, hgvsc: str) -> str:
        allele, span, cdna, _, _, _ = cases[name]
        extra = f"IMPACT=MODIFIER;STRAND=1;HGVSc={hgvsc}"
        return "\t".join([name, location(i, span), allele, "ENSG1", "ENST1", "Transcript", "5_prime_UTR_variant", cdna, "-", "-", "-", "-", "-", extra]) + "\n"

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, i, case[3]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    rs.write_text("".join(row(name, i, case[4]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs), "--field", "HGVSc"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    klass = "hgvs_clipped_reference_against_whole_alternate"
    assert {d["location"]: d["expected_divergence"] for d in m["field_divergences"]} == {
        location(i, case[1]): (klass if case[5] else "unexplained_residual") for i, case in enumerate(cases.values())
    }
    assert m["field_divergence_summary"] == {klass: 2, "unexplained_residual": 1}


def test_classify_clipped_reference_against_whole_alternate_names_the_intron_offset_fields_beside_its_notation(tmp_path: Path):
    """The HGVSIntronOffset plugin's two fields that Ensembl VEP alone writes beside an `ins` or
    `dup` the class takes against a vep-rs row without an `HGVSc` are
    `hgvs_clipped_reference_against_whole_alternate` where they are the intron parts of the
    notation's coordinates as the plugin reads them: the `ins` the clip leaves from its 3'
    flank first, which prints second unless the 5' flank lies past the stop codon, a `dup` in
    its printed order, as is the `del`, `delins` or substitution of an allele shorter than its
    span.
    The two values the other way round, a row whose cDNA position is closed, a vep-rs row that
    writes a field and the `delins` of an allele as long as its span stay
    `unexplained_residual`, and the builder reads each fabricated row as the comparator's rule
    does."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    corpus = tmp_path / "c"
    feature = "ENST1"
    cases = {
        # name: (cDNA position, VEP HGVSc, VEP start and end offsets, vep-rs extra, class)
        "rs_star": ("553-?", "c.*234+1_*234insT", ("1", "0"), "", True),
        "rs_turned": ("553-?", "c.76-1_76insT", ("0", "-1"), "", True),
        "rs_dup": ("553-?", "c.76-2_76-1dup", ("-2", "-1"), "", True),
        "rs_other_way": ("553-?", "c.76-1_76insT", ("-1", "0"), "", False),
        "rs_closed": ("553-557", "c.*234+1_*234insT", ("1", "0"), "", False),
        "rs_written": ("553-?", "c.*234+1_*234insT", ("1", "0"), ";HGVS_IntronEndOffset=0", False),
        "rs_delins": ("553-?", "c.*234+1_*234delinsTG", ("1", "0"), "", False),
        "rs_shorter_delins": ("553-?", "c.*233+2_*234delinsT", ("2", "0"), "", True),
        "rs_shorter_del": ("553-?", "c.*232_*234del", ("0", "0"), "", True),
        "rs_shorter_del_other_way": ("553-?", "c.*233+2_*234del", ("0", "2"), "", False),
    }
    shorter = {name for name in cases if name.startswith("rs_shorter")}
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\t{name}\tT\tTA\t.\tPASS\t.") for i, name in enumerate(cases)])

    def row(name: str, i: int, extra: str) -> str:
        location = f"21:{100 * (i + 1)}-{100 * (i + 1) + 4}"
        columns = [name, location, "TG" if name in shorter else "TGCAT", "ENSG1", feature, "Transcript", "3_prime_UTR_variant", cases[name][0], "-", "-", "-", "-", "-"]
        return "\t".join([*columns, f"IMPACT=MODIFIER;STRAND=-1{extra}"]) + "\n"

    def vep_extra(case: tuple) -> str:
        return f";HGVSc={feature}.1:{case[1]};HGVS_IntronStartOffset={case[2][0]};HGVS_IntronEndOffset={case[2][1]}"

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, i, vep_extra(case)) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    rs.write_text("".join(row(name, i, case[3]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    fields = ["--field", "HGVSc", "--field", "HGVS_IntronStartOffset", "--field", "HGVS_IntronEndOffset"]
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs), *fields]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    klass = "hgvs_clipped_reference_against_whole_alternate"
    named = {(d["location"], d["field"]): d["expected_divergence"] == klass for d in m["field_divergences"]}
    header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
    for i, (name, case) in enumerate(cases.items()):
        location = f"21:{100 * (i + 1)}-{100 * (i + 1) + 4}"
        e, a = (cvf.parse_default(header + row(name, i, extra))[0] for extra in (vep_extra(case), case[3]))
        for field in cvf.HGVS_INTRON_OFFSET_FIELDS:
            if cvf.field_of(a, (field,)) == cvf.field_of(e, (field,)):
                continue
            assert named[location, field] == case[4], (name, field)
            assert cvf.field_class_accepts(klass, e, field, cvf.field_of(e, (field,)), cvf.field_of(a, (field,)), a) == case[4], (name, field)


def test_classify_clipped_reference_against_whole_alternate_names_the_intron_offset_fields_on_vep_rs_side(tmp_path: Path):
    """The HGVSIntronOffset plugin's two fields that vep-rs alone writes, beside its `ins` of the
    bases a longer allele adds inside the span on a reverse-strand transcript's row whose cDNA
    position is open at its start and whose Ensembl VEP `HGVSc` is empty, are
    `hgvs_clipped_reference_against_whole_alternate`, while the `HGVSc` of such a row is not
    the class's; the differing field of a row whose two `ins` notations the class takes as a
    pair is the class's too, each side's fields the intron parts of its own notation. A `del`
    on vep-rs's side, an `ins` of another length or away from the bases inside, a
    forward-strand transcript's row, a reference `dup` beside vep-rs's fields, an `HGVS_OFFSET`
    on the row
    and Ensembl VEP's two values the other way round stay `unexplained_residual`, and the
    builder reads each fabricated row, its `HGVSc` included, as the comparator's rule does."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    corpus = tmp_path / "c"
    feature = "ENST1"
    long_allele, added = "ACGTACGTAGACGTACGTAG", "ACGTACGTAG"
    zeros = ";HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=0"
    inside = f";HGVSc={feature}.1:n.7_8ins{added}{zeros}"
    outside_ins = f";HGVSc={feature}.1:c.-529-1_-529insACGTACGTACGT"
    paired = f";HGVSc={feature}.1:c.-528_-527insACGTACG{zeros}"
    cases = {
        # name: (cDNA position, STRAND, allele, VEP extra, vep-rs extra, the fields' class, the HGVSc's class)
        "rs_inside": ("?-7", "-1", long_allele, "", inside, True, False),
        "rs_forward_end": ("7-?", "1", long_allele, "", inside, False, False),
        "rs_genomic_start": ("?-7", "1", long_allele, "", inside, False, False),
        "rs_del": ("?-7", "-1", long_allele, "", f";HGVSc={feature}.1:n.1del{zeros}", False, True),
        "rs_fewer": ("?-7", "-1", long_allele, "", inside.replace(f"ins{added}", f"ins{added[:-1]}"), False, False),
        "rs_past_the_end": ("?-7", "-1", long_allele, "", inside.replace("n.7_8", "n.8_9"), False, False),
        # The row's four fields are those of `hgvs_shifted_indel_written_span_gate`, asked first.
        "rs_offset": ("?-7", "-1", long_allele, "", inside + ";HGVS_OFFSET=2", False, None),
        "rs_both": ("?-2", "1", "ACGTACGACGTACG", outside_ins + ";HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=-1", paired, True, True),
        "rs_both_turned": ("?-2", "1", "ACGTACGACGTACG", outside_ins + ";HGVS_IntronStartOffset=-1;HGVS_IntronEndOffset=1", paired, False, True),
        "rs_both_dup": ("?-2", "1", "ACGTACGACGTACG", f";HGVSc={feature}.1:c.-529-1dup;HGVS_IntronStartOffset=-1;HGVS_IntronEndOffset=-1", paired, False, True),
        "rs_both_del": ("?-2", "1", "ACGTACGACGTACG", outside_ins + ";HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=-1",
                        f";HGVSc={feature}.1:c.-528del{zeros}", False, True),
    }
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\t{name}\tT\tTA\t.\tPASS\t.") for i, name in enumerate(cases)])

    def row(name: str, i: int, extra: str) -> str:
        position, strand, allele = cases[name][:3]
        location = f"21:{100 * (i + 1)}-{100 * (i + 1) + len(allele) // 2 - 1}"
        columns = [name, location, allele, "ENSG1", feature, "Transcript", "5_prime_UTR_variant", position, "-", "-", "-", "-", "-"]
        return "\t".join([*columns, f"IMPACT=MODIFIER;STRAND={strand}{extra}"]) + "\n"

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, i, case[3]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    rs.write_text("".join(row(name, i, case[4]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    fields = ["--field", "HGVSc", "--field", "HGVS_OFFSET", "--field", "HGVS_IntronStartOffset", "--field", "HGVS_IntronEndOffset"]
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs), *fields]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    klass = "hgvs_clipped_reference_against_whole_alternate"
    named = {(d["location"], d["field"]): d["expected_divergence"] == klass for d in m["field_divergences"]}
    header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
    for i, (name, case) in enumerate(cases.items()):
        e, a = (cvf.parse_default(header + row(name, i, extra))[0] for extra in (case[3], case[4]))
        for field, expected in (("HGVSc", case[6]), *((f, case[5]) for f in cvf.HGVS_INTRON_OFFSET_FIELDS)):
            if cvf.field_of(a, (field,)) == cvf.field_of(e, (field,)) or expected is None:
                continue
            assert named[e.location, field] == expected, (name, field)
            assert cvf.field_class_accepts(klass, e, field, cvf.field_of(e, (field,)), cvf.field_of(a, (field,)), a) == expected, (name, field)


def test_classify_clipped_reference_against_whole_alternate_names_the_end_offset_before_a_forward_strand_first_base(tmp_path: Path):
    """The HGVSIntronOffset plugin's end offset that Ensembl VEP writes as -1 beside its `ins`
    from one base before a forward-strand transcript's first base, where vep-rs writes 0
    beside its duplication of the bases the allele adds, is
    `hgvs_clipped_reference_against_whole_alternate` where the two rows show the clamp at the
    first base: the `ins` writes the allele less as many of its last bases as the cDNA
    position's closed end counts, fewer than the span holds, and the first base's coordinate
    lies that count and vep-rs's `HGVS_OFFSET` before the duplication's last base. The
    `HGVS_OFFSET` of such a row, and its `HGVSc` where the duplication lies past the span's
    reach, stay `unexplained_residual`. So does the end offset of a reverse-strand
    transcript's row, of a row whose closed end is as large as its span, beside the `ins` the
    clip leaves at a true intron, which writes the added bases alone, beside an `ins` of
    fewer bases, a reference `HGVS_OFFSET`, a duplication of another count, across the start
    codon or at another count from the first base, and where either row carries `BAM_EDIT`; the
    start offset of a plain insertion at an exon's first base, whose pair is -1 and 0, stays
    one too. The builder reads each fabricated row as the comparator's rule does."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    corpus = tmp_path / "c"
    feature = "ENST1"
    allele = "GATTCACGATTCAC"
    first = f";HGVSc={feature}.1:c.-529-1_-529ins"
    turned = ";HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=-1"
    zeros = ";HGVS_IntronStartOffset=0;HGVS_IntronEndOffset=0"
    vep = first + allele[:12] + turned
    carried = f";HGVSc={feature}.1:c.-517_-511dup;HGVS_OFFSET=17{zeros}"
    plain = f";HGVSc={feature}.1:c.784-1_784ins{allele[:7]};HGVS_IntronStartOffset=-1;HGVS_IntronEndOffset=0"
    cases = {
        # name: (span, cDNA position, STRAND and shared extra, allele, VEP extra, vep-rs extra, the end offset's class)
        "rs_first_base": (7, "?-2", "1", allele, vep, carried, True),
        "rs_not_carried": (7, "?-2", "1", allele[:8], first + allele[:6] + turned, f";HGVSc={feature}.1:c.-528dup{zeros}", True),
        "rs_reverse": (7, "?-2", "-1", allele, vep, carried, False),
        "rs_whole_span_inside": (7, "?-7", "1", allele, first + allele[:7] + turned, carried.replace("=17", "=12"), False),
        "rs_true_intron": (7, "?-2", "1", allele, f";HGVSc={feature}.1:c.88-1_88ins{allele[5:12]}{turned}",
                           f";HGVSc={feature}.1:c.100_106dup;HGVS_OFFSET=17{zeros}", False),
        "rs_fewer": (7, "?-2", "1", allele, first + allele[:11] + turned, carried, False),
        "rs_reference_offset": (7, "?-2", "1", allele, vep.replace(turned, ";HGVS_OFFSET=2" + turned), carried, False),
        "rs_six": (7, "?-2", "1", allele, vep, carried.replace("-517_", "-516_"), False),
        "rs_one_short": (7, "?-2", "1", allele, vep, carried.replace("=17", "=16"), False),
        "rs_across_the_start": (7, "?-2", "1", allele, vep, carried.replace("c.-517_-511dup", "c.5_11dup").replace("=17", "=538"), False),
        "rs_edited": (7, "?-2", "1;BAM_EDIT=OK", allele, vep, carried, False),
        "rs_edited_reference_alone": (7, "?-2", "1", allele, ";BAM_EDIT=OK" + vep, carried, False),
        "rs_edited_vep_rs_alone": (7, "?-2", "1", allele, vep, ";BAM_EDIT=OK" + carried, False),
        "rs_plain_insertion": (2, "900-901", "1", allele[:7], plain, f";HGVSc={feature}.1:c.790_796dup;HGVS_OFFSET=13{zeros}", False),
    }
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\t{name}\tT\tTA\t.\tPASS\t.") for i, name in enumerate(cases)])

    def row(name: str, i: int, extra: str) -> str:
        span, position, strand, bases = cases[name][:4]
        columns = [name, f"21:{100 * (i + 1)}-{100 * (i + 1) + span - 1}", bases, "ENSG1", feature, "Transcript", "5_prime_UTR_variant",
                   position, "-", "-", "-", "-", "-"]
        return "\t".join([*columns, f"IMPACT=MODIFIER;STRAND={strand}{extra}"]) + "\n"

    vep_output = tmp_path / "vep.txt"
    rs_output = tmp_path / "rs.txt"
    vep_output.write_text("".join(row(name, i, case[4]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    rs_output.write_text("".join(row(name, i, case[5]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    fields = ["--field", "HGVSc", "--field", "HGVS_OFFSET", "--field", "HGVS_IntronStartOffset", "--field", "HGVS_IntronEndOffset"]
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep_output), "--vep-rs-output", str(rs_output), *fields]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    klass = "hgvs_clipped_reference_against_whole_alternate"
    named = {(d["location"], d["field"]): d["expected_divergence"] for d in m["field_divergences"]}
    header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
    asked = 0
    for i, (name, case) in enumerate(cases.items()):
        e, a = (cvf.parse_default(header + row(name, i, extra))[0] for extra in (case[4], case[5]))
        for field in cvf.HGVS_INTRON_OFFSET_FIELDS:
            if cvf.field_of(a, (field,)) == cvf.field_of(e, (field,)):
                continue
            expected = case[6] and field == "HGVS_IntronEndOffset"
            asked += 1
            assert (named[e.location, field] == klass) == expected, (name, field)
            assert cvf.field_class_accepts(klass, e, field, cvf.field_of(e, (field,)), cvf.field_of(a, (field,)), a) == expected, (name, field)
    assert asked == len(cases)
    location = "21:100-106"
    assert named[location, "HGVSc"] == named[location, "HGVS_OFFSET"] == "unexplained_residual"


def test_clipped_body_is_the_comparators_clipped_notation_pattern():
    """The builder's `HGVSC_CLIPPED_BODY` and the comparator's `CLIPPED_NOTATION_RE` read the
    same notations, so a row the comparator masks under the clipped class is one the builder
    classifies; a change to either alphabet is made in both."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    assert bgc.HGVSC_CLIPPED_BODY.pattern == cvf.CLIPPED_NOTATION_RE.pattern
    assert bgc.HGVSC_CLIPPED_BODY.flags == cvf.CLIPPED_NOTATION_RE.flags


def test_classify_names_the_added_bases_written_outside_the_transcript_class(tmp_path: Path):
    """An `HGVSc` and an `HGVS_OFFSET` that vep-rs writes and Ensembl VEP leaves empty, on the
    row of inserted bases both engines annotate upstream or downstream of the transcript, are
    `hgvs_added_bases_written_outside_the_transcript`: a `dup` or an `ins` of the row's
    transcript, of as many bases as the allele holds or adds to its span, beside an offset larger
    than the row's DISTANCE on an upstream row, or of its size where the allele adds its bases
    to a span it keeps, negative on a reverse strand, and a `dup` beside
    no offset on a downstream row. A row with a term inside the transcript is the written-span
    gate's. A reference HGVSc or HGVS_OFFSET, an upstream row with no offset, with one within
    the DISTANCE or with one of the other strand's sign, a downstream row with an offset or an
    `ins`, a notation of another kind, count or transcript, an allele of no more bases than its
    span, a deletion's row and a row without a DISTANCE stay an `unexplained_residual`, and the
    builder takes each fabricated row as the comparator's class does."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    klass, gate, residual = "hgvs_added_bases_written_outside_the_transcript", "hgvs_shifted_indel_written_span_gate", "unexplained_residual"
    upstream, downstream = "upstream_gene_variant", "downstream_gene_variant"
    corpus = tmp_path / "c"
    cases = {
        # name: (term, STRAND, DISTANCE, allele, positions, VEP's own fields, vep-rs HGVSc, vep-rs HGVS_OFFSET, class)
        "rs_upstream": (upstream, "1", "9", "TA", 2, "", "ENST1.1:c.-103_-102dup", "18", klass),
        "rs_beside": (upstream, "1", "0", "T", 2, "", "ENST1.1:n.14dup", "14", klass),
        "rs_reverse": (upstream, "-1", "5", "T", 2, "", "ENST1.1:n.14dup", "-14", klass),
        "rs_intron": (upstream, "1", "5", "T", 2, "", "ENST1.1:c.20+6dup", "26", klass),
        "rs_ins": (upstream, "1", "5", "T", 2, "", "ENST1.1:c.20+8_20+9insT", "28", klass),
        "rs_added": (upstream, "1", "5", "GTT", 2, "", "ENST1.1:n.14dup", "14", klass),
        "rs_added_wider": (upstream, "1", "5", "GTTA", 3, "", "ENST1.1:n.14dup", "14", klass),
        "rs_downstream": (downstream, "1", "0", "T", 2, "", "ENST1.1:n.2601dup", "", klass),
        "rs_inside": ("intron_variant", "1", "", "T", 2, "", "ENST1.1:n.14dup", "14", gate),
        "rs_reference_hgvsc": (upstream, "1", "5", "T", 2, ";HGVSc=ENST1.1:n.1_2insT", "ENST1.1:n.14dup", "14", residual),
        "rs_reference_offset": (upstream, "1", "5", "T", 2, ";HGVS_OFFSET=2", "ENST1.1:n.14dup", "14", residual),
        "rs_no_offset": (upstream, "1", "5", "T", 2, "", "ENST1.1:n.14dup", "", residual),
        "rs_within_distance": (upstream, "1", "5", "T", 2, "", "ENST1.1:n.14dup", "5", residual),
        "rs_kept_at_the_distance": (upstream, "1", "5", "GT", 1, "", "ENST1.1:n.1dup", "5", klass),
        "rs_kept_within_distance": (upstream, "1", "5", "GT", 1, "", "ENST1.1:n.1dup", "4", residual),
        "rs_added_at_the_distance": (upstream, "1", "5", "GTT", 2, "", "ENST1.1:n.1dup", "5", klass),
        "rs_three_at_the_distance": (upstream, "1", "5", "GTT", 2, "", "ENST1.1:n.1_3dup", "5", residual),
        "rs_other_sign": (upstream, "1", "5", "T", 2, "", "ENST1.1:n.14dup", "-14", residual),
        "rs_reverse_unsigned": (upstream, "-1", "5", "T", 2, "", "ENST1.1:n.14dup", "14", residual),
        "rs_downstream_offset": (downstream, "1", "0", "T", 2, "", "ENST1.1:n.2601dup", "1", residual),
        "rs_downstream_ins": (downstream, "1", "0", "T", 2, "", "ENST1.1:n.2601_2602insT", "", residual),
        "rs_del": (upstream, "1", "5", "T", 2, "", "ENST1.1:n.14del", "14", residual),
        "rs_dup_count": (upstream, "1", "5", "T", 2, "", "ENST1.1:n.13_14dup", "14", residual),
        "rs_ins_count": (upstream, "1", "5", "T", 2, "", "ENST1.1:c.20+8_20+9insTT", "28", residual),
        "rs_other": (upstream, "1", "5", "T", 2, "", "ENST2.1:n.14dup", "14", residual),
        "rs_kept_two": (upstream, "1", "5", "GT", 2, "", "ENST1.1:n.14dup", "14", residual),
        "rs_substitution": (upstream, "1", "5", "T", 1, "", "ENST1.1:n.14dup", "14", residual),
        "rs_deletion": (upstream, "1", "5", "-", 2, "", "ENST1.1:n.14dup", "14", residual),
        "rs_no_distance": (upstream, "1", "", "T", 2, "", "ENST1.1:n.14dup", "14", residual),
    }
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\t{name}\tT\tTA\t.\tPASS\t.") for i, name in enumerate(cases)])

    def location(i: int, positions: int) -> str:
        start = 100 * (i + 1)
        return f"21:{start}" if positions == 1 else f"21:{start}-{start + positions - 1}"

    def row(name: str, i: int, own: str) -> str:
        term, strand, distance, allele, positions = cases[name][:5]
        extra = "IMPACT=MODIFIER" + (f";DISTANCE={distance}" if distance else "") + f";STRAND={strand}{own}"
        return _row(name, location(i, positions), allele, "ENST1", term, extra)

    def vep_rs_own(case: tuple) -> str:
        return f";HGVSc={case[6]}" + (f";HGVS_OFFSET={case[7]}" if case[7] else "")

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, i, case[5]) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    rs.write_text("".join(row(name, i, vep_rs_own(case)) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "HGVSc", "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {(d["location"], d["field"]): d["expected_divergence"] for d in m["field_divergences"]}
    expected = {}
    for i, case in enumerate(cases.values()):
        expected[location(i, case[4]), "HGVSc"] = case[8]
        if case[7]:
            expected[location(i, case[4]), "HGVS_OFFSET"] = case[8]
    assert classes == expected
    assert m["field_divergence_summary"] == {klass: 19, gate: 2, residual: 34}
    assert {d["vep_rs"] for d in m["field_divergences"] if d["expected_divergence"] == klass and d["field"] == "HGVS_OFFSET"} == {
        "18", "14", "-14", "26", "28", "5"}
    header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
    for i, (name, case) in enumerate(cases.items()):
        e, a = (cvf.parse_default(header + row(name, i, own))[0] for own in (case[5], vep_rs_own(case)))
        for field in ("HGVSc", "HGVS_OFFSET"):
            ev, av = cvf.field_of(e, (field,)), cvf.field_of(a, (field,))
            if ev != av:
                assert cvf.field_class_accepts(klass, e, field, ev, av, a) == (case[8] == klass), (name, field)


def test_classify_added_bases_written_outside_the_transcript_holds_the_intron_offset_fields_to_its_notation(tmp_path: Path):
    """In a run that compares the HGVSIntronOffset plugin's two fields, a row of
    `hgvs_added_bases_written_outside_the_transcript` carries on vep-rs's side the intron offsets
    of its notation's coordinates in their printed order, and each of the two that differs from
    Ensembl VEP's, which writes a pair of its own or neither (HGVSIntronOffset.pm:62), is the
    class beside the `HGVSc` and the `HGVS_OFFSET`; a row on which neither engine writes a
    field, as in a run without the plugin, is the class in its `HGVSc` and `HGVS_OFFSET`.
    Fields that are not the notation's, one field alone on either side and a vep-rs row that
    writes neither beside an Ensembl VEP row that writes both stay `unexplained_residual` in
    every field, and the builder takes each fabricated row as the comparator's class does."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    klass = "hgvs_added_bases_written_outside_the_transcript"
    upstream, downstream = "upstream_gene_variant", "downstream_gene_variant"
    corpus = tmp_path / "c"
    feature = "ENST1"
    cases = {
        # name: (term, STRAND, vep-rs HGVSc, vep-rs start and end offsets, vep-rs HGVS_OFFSET, VEP start and end offsets, class)
        "rs_dup": (upstream, "1", "n.14dup", ("0", "0"), "14", None, True),
        "rs_same_pair": (upstream, "1", "n.14dup", ("0", "0"), "14", ("0", "0"), True),
        "rs_other_pair": (upstream, "1", "n.14dup", ("0", "0"), "14", ("1", "0"), True),
        "rs_reverse": (upstream, "-1", "n.14dup", ("0", "0"), "-14", ("0", "-1"), True),
        "rs_intron_dup": (upstream, "1", "c.20+6dup", ("6", "6"), "26", None, True),
        "rs_intron_ins": (upstream, "1", "c.20+8_20+9insT", ("8", "9"), "28", None, True),
        "rs_before_exon": (upstream, "1", "c.21-9_21-8insT", ("-9", "-8"), "28", ("0", "0"), True),
        "rs_downstream": (downstream, "1", "n.2601dup", ("0", "0"), "", ("1", "0"), True),
        "rs_neither_writes": (upstream, "1", "n.14dup", None, "14", None, True),
        "rs_not_its_coordinates": (upstream, "1", "c.20+8_20+9insT", ("9", "9"), "28", None, False),
        "rs_fields_not_written": (upstream, "1", "n.14dup", None, "14", ("0", "0"), False),
        "rs_one_field_written": (upstream, "1", "n.14dup", ("0", None), "14", None, False),
        "rs_one_reference_field": (upstream, "1", "n.14dup", ("0", "0"), "14", ("1", None), False),
    }
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\t{name}\tT\tTA\t.\tPASS\t.") for i, name in enumerate(cases)])

    def row(name: str, i: int, own: str) -> str:
        term, strand = cases[name][:2]
        return _row(name, f"21:{100 * (i + 1)}-{100 * (i + 1) + 1}", "T", feature, term, f"IMPACT=MODIFIER;DISTANCE=5;STRAND={strand}{own}")

    def pair(offsets: tuple | None) -> str:
        return "".join(f";{field}={value}" for field, value in zip(bgc.HGVS_INTRON_OFFSET_FIELDS, offsets or ()) if value is not None)

    def vep_rs_own(case: tuple) -> str:
        return f";HGVSc={feature}.1:{case[2]}" + (f";HGVS_OFFSET={case[4]}" if case[4] else "") + pair(case[3])

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, i, pair(case[5])) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    rs.write_text("".join(row(name, i, vep_rs_own(case)) for i, (name, case) in enumerate(cases.items())), encoding="utf-8")
    fields = ("HGVSc", "HGVS_OFFSET", *bgc.HGVS_INTRON_OFFSET_FIELDS)
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     *(word for field in fields for word in ("--field", field))]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    named = {(d["location"], d["field"]): d["expected_divergence"] == klass for d in m["field_divergences"]}
    header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"
    asked = 0
    for i, (name, case) in enumerate(cases.items()):
        location = f"21:{100 * (i + 1)}-{100 * (i + 1) + 1}"
        e, a = (cvf.parse_default(header + row(name, i, own))[0] for own in (pair(case[5]), vep_rs_own(case)))
        for field in fields:
            ev, av = cvf.field_of(e, (field,)), cvf.field_of(a, (field,))
            if ev == av:
                continue
            asked += 1
            assert named[location, field] == case[6], (name, field)
            assert cvf.field_class_accepts(klass, e, field, ev, av, a) == case[6], (name, field)
    assert asked == len(named) == 43
    assert m["field_divergence_summary"] == {klass: 28, "unexplained_residual": 15}


def test_exonic_run_counts_a_notations_bases_as_the_comparator_does():
    """The builder's `exonic_run` counts the bases from one cDNA coordinate to another as the
    comparator's does, none where a coordinate carries an intron offset or the two lie in
    different regions of the transcript, and the builder takes a coordinate's intron offset as
    the comparator's `coordinate_intron_offset` does, so a notation the comparator counts for
    `hgvs_added_bases_written_outside_the_transcript` is one the builder counts."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    assert bgc.exonic_run("-103", "-95") == 9
    assert bgc.exonic_run("14", "14") == 1
    assert bgc.exonic_run("*5", "*9") == 5
    assert bgc.exonic_run("20+6", "20+6") is None
    assert bgc.exonic_run("-3", "2") is None
    assert bgc.exonic_run("20", "*1") is None
    coordinates = ["-103", "-95", "-49-1", "-5+2", "-1", "1", "14", "20+6", "21-9", "300", "*1", "*5", "*9", "*12-3"]
    for first in coordinates:
        assert str(bgc.cdna_coordinate_order(first)[2]) == cvf.coordinate_intron_offset(first), first
        for last in coordinates:
            assert bgc.exonic_run(first, last) == cvf.exonic_run(first, last), (first, last)


def _intron_end_transcript(number: int, strand: int, coding: tuple[int, int], *exons: tuple[int, int]) -> dict:
    """A cache transcript `ENST<number>` whose exons, in the transcript's order, lie `exons` past
    `number` times 10,000, coding from the first to the second cDNA position of `coding`."""
    base = number * 10_000
    return {"stable_id": f"ENST{number:011d}", "strand": strand, "cdna_coding_start": coding[0], "cdna_coding_end": coding[1],
            "exons": [{"start": base + start, "end": base + end} for start, end in exons]}


# Transcripts 1 to 14 (forward) and 21 (reverse) have three exons of 100 bases and code from
# their 51st base to their 250th, so the second exon opens on coding base 51 behind an intron of
# 100 bases (forward) or 200 (reverse). Transcript 22 codes from its 151st base to its 250th: its
# second exon opens in the 5' UTR and its fourth in the 3' UTR. Transcripts 23 and 24 code over
# their whole length behind introns of 13, 2 and 14 bases, the fourth exon opening on coding
# base 286.
_INTRON_END_TRANSCRIPTS = [
    *(_intron_end_transcript(number, 1, (51, 250), (1000, 1099), (1200, 1299), (1500, 1599)) for number in range(1, 15)),
    _intron_end_transcript(21, -1, (51, 250), (1500, 1599), (1200, 1299), (1000, 1099)),
    _intron_end_transcript(22, 1, (151, 250), (1000, 1099), (1200, 1299), (1400, 1499), (1600, 1699)),
    *(_intron_end_transcript(number, 1, (1, 371), (1000, 1099), (1113, 1199), (1202, 1299), (1314, 1399)) for number in (23, 24)),
]


def test_classify_names_the_intron_end_insertion_class(tmp_path: Path):
    """An `HGVSp` Ensembl VEP writes and vep-rs leaves empty, on an insertion written between an
    intron's last base and the first coding base of the next exon, is
    `hgvsp_intron_end_insertion` for each form of the row's `HGVSc`: a copy of the intron's last
    base, a duplication of two or more intron positions ending on it, numbered from either side
    of the intron, an insertion of one base or of bases that do not end in a pair that ends an
    intron, and an insertion of two or more bases that end in either such pair, the pair alone
    included, on either strand. A row both engines give an `HGVS_OFFSET`, an exonic `HGVSc`, an
    intron's end before a 5' UTR or a 3' UTR position, two `HGVSc` values that differ, the end
    of an intron of 13 bases and of 2 bases, a transcript the cache does not hold, a vep-rs
    notation that is filled and differs, a reference without a notation and a Location whose
    span lies beside no gap stay an `unexplained_residual`, the insertion ending in a pair
    among them where both engines give an `HGVS_OFFSET` and at the end of either short intron,
    as does every row of a run that does not compare `HGVSc` and `HGVS_OFFSET` beside `HGVSp`."""
    corpus = tmp_path / "c"
    protein = "ENSP00000000001.1:p.Gly17ArgfsTer5"
    cases = {
        # name: (transcript, the lower position beside the gap, allele, HGVSc body, class)
        "rec_copy": (1, 11199, "G", "c.51-1dup", True),
        "rec_duplication": (2, 21199, "CAG", "c.51-3_51-1dup", True),
        "rec_first_half": (23, 231313, "TTCTTCAG", "c.285+7_286-1dup", True),
        "rec_one_base": (3, 31199, "T", "c.51-1_51insT", True),
        "rec_insertion": (4, 41199, "TAT", "c.51-1_51insTAT", True),
        "rec_reverse": (21, 211299, "C", "c.51-1dup", True),
        "rec_major_pair": (5, 51199, "TAG", "c.51-1_51insTAG", True),
        "rec_minor_pair": (6, 61199, "TAC", "c.51-1_51insTAC", True),
        "rec_pair_alone": (14, 141199, "AG", "c.51-1_51insAG", True),
        "rec_pair_reverse": (21, 211099, "GTA", "c.151-1_151insTAC", True),
        "rec_pair_moved": (13, 131199, "TAG", "c.51-1_51insTAG", False),
        "rec_pair_13_bases": (24, 241112, "TAG", "c.101-1_101insTAG", False),
        "rec_pair_2_bases": (24, 241201, "TAC", "c.188-1_188insTAC", False),
        "rec_moved": (7, 71199, "G", "c.51-1dup", False),
        "rec_exonic": (8, 81249, "T", "c.100_101insT", False),
        "rec_utr5": (22, 221199, "G", "c.-50-1dup", False),
        "rec_utr3": (22, 221599, "G", "c.*51-1dup", False),
        "rec_two_notations": (9, 91199, "G", "c.51-1dup", False),
        "rec_13_bases": (23, 231112, "G", "c.101-1dup", False),
        "rec_2_bases": (23, 231201, "G", "c.188-1dup", False),
        "rec_absent": (99, 991199, "G", "c.51-1dup", False),
        "rec_other_protein": (10, 101199, "G", "c.51-1dup", False),
        "rec_no_reference": (11, 111199, "G", "c.51-1dup", False),
        "rec_wide": (12, 121199, "G", "c.51-1dup", False),
    }
    # The fields a case's Ensembl VEP row and vep-rs row write otherwise, and the Location of the
    # case written over three positions.
    changed = {
        "rec_moved": ({"HGVS_OFFSET": "1"}, {"HGVS_OFFSET": "1"}),
        "rec_pair_moved": ({"HGVS_OFFSET": "1"}, {"HGVS_OFFSET": "1"}),
        "rec_two_notations": ({}, {"HGVSc": "ENST00000000009.1:c.51-2_51-1dup"}),
        "rec_other_protein": ({}, {"HGVSp": protein.replace("fsTer5", "fsTer9")}),
        "rec_no_reference": ({"HGVSp": ""}, {"HGVSp": protein}),
    }
    locations = {name: f"1:{position}-{position + 1}" for name, (_, position, *_) in cases.items()}
    locations["rec_wide"] = "1:121201-121203"
    _manifest(corpus, [_record(f"1\t{position}\t{name}\tG\tG{allele}\t.\tPASS\t.") for name, (_, position, allele, *_) in cases.items()])
    (corpus / "json_cache" / "transcripts" / "1").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "1" / "1-1000000.json").write_text(json.dumps(_INTRON_END_TRANSCRIPTS), encoding="utf-8")

    def row(name: str, side: int) -> str:
        number, _, allele, body, _ = cases[name]
        feature = f"ENST{number:011d}"
        fields = {"HGVSc": f"{feature}.1:{body}", "HGVSp": "" if side else protein, "HGVS_OFFSET": "", **changed.get(name, ({}, {}))[side]}
        extra = f"IMPACT=HIGH;STRAND={-1 if number == 21 else 1}" + "".join(f";{field}={value}" for field, value in fields.items() if value)
        return _row(name, locations[name], allele, feature, "frameshift_variant,splice_region_variant", extra)

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, 0) for name in cases), encoding="utf-8")
    rs.write_text("".join(row(name, 1) for name in cases), encoding="utf-8")
    command = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main([*command, "--field", "HGVSc", "--field", "HGVSp", "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    klass = "hgvsp_intron_end_insertion"
    assert {(d["location"], d["field"]): d["expected_divergence"] for d in m["field_divergences"]} == {
        **{(locations[name], "HGVSp"): (klass if case[4] else "unexplained_residual") for name, case in cases.items()},
        (locations["rec_two_notations"], "HGVSc"): "unexplained_residual",
    }
    assert m["field_divergence_summary"] == {klass: 10, "unexplained_residual": 15}
    assert {d["vep_rs"] for d in m["field_divergences"] if d["expected_divergence"] == klass} == {"-"}
    # The class reads the row's HGVSc and HGVS_OFFSET from the compared fields: a run that
    # compares HGVSp alone names nothing.
    assert bgc.main([*command, "--field", "HGVSp"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert m["field_divergence_summary"] == {"unexplained_residual": 24}
    key = ("1:11199-11200", "G", "ENST00000000001", "Transcript", "frameshift_variant,splice_region_variant")
    fields = {"HGVSc": "ENST00000000001.1:c.51-1dup", "HGVSp": protein, "HGVS_OFFSET": "-"}
    gaps = bgc.intron_end_gaps(corpus / "json_cache")
    assert gaps[("1", "ENST00000000001")] == {11199, 11499} and gaps[("1", "ENST00000000021")] == {211099, 211299}
    assert gaps[("1", "ENST00000000022")] == {221399} and gaps[("1", "ENST00000000023")] == {231313}

    def named(perl_row: dict, rust_row: dict, intron_ends: dict | None = gaps, name: str = "HGVSp", key: tuple = key) -> str | None:
        found = bgc.classify_field_divergence(name, key, perl_row[name], rust_row[name], "1", [], perl_row=perl_row, rust_row=rust_row, intron_ends=intron_ends)
        return found[0] if found else None

    assert named(fields, {**fields, "HGVSp": "-"}) == klass
    assert named(fields, {**fields, "HGVSp": "-"}, key=("chr" + key[0], *key[1:])) == klass
    assert named(fields, {**fields, "HGVSp": "-"}, None) is None and named(fields, {**fields, "HGVSp": "-"}, {}) is None
    assert named(fields, {**fields, "HGVSp": "-", "HGVS_OFFSET": "1"}) is None
    assert named({**fields, "HGVS_OFFSET": "1"}, {**fields, "HGVSp": "-"}) is None
    assert named(fields, {**fields, "HGVSp": "-"}, name="HGVSc") is None
    without_offset = {field: value for field, value in fields.items() if field != "HGVS_OFFSET"}
    assert named(without_offset, {**without_offset, "HGVSp": "-"}) is None


# Chromosome 1 transcripts of one exon layout: an exon of three bases, 201 to 203, between an
# intron from 101 to 200 and one from 204 to 300, on the forward strand (ENST1) and on the
# reverse strand (ENST2); ENST3 has one intron from 101 to 300.
_MAXENTSCAN_EXONS = [{"start": 1, "end": 100}, {"start": 201, "end": 203}, {"start": 301, "end": 400}]
_MAXENTSCAN_TRANSCRIPTS = [
    {"stable_id": "ENST1", "strand": 1, "exons": _MAXENTSCAN_EXONS},
    {"stable_id": "ENST2", "strand": -1, "exons": _MAXENTSCAN_EXONS},
    {"stable_id": "ENST3", "strand": 1, "exons": [_MAXENTSCAN_EXONS[0], _MAXENTSCAN_EXONS[2]]},
]
# A donor window's nine bases and an acceptor window's twenty-three.
_MAXENTSCAN_DONOR_TEXT, _MAXENTSCAN_ACCEPTOR_TEXT = "AC" * 4 + "G", "AC" * 11 + "G"


def test_classify_names_maxentscan_overlapping_windows_where_the_cache_places_the_record_in_two_windows(tmp_path: Path):
    """`classify --field` names `maxentscan_overlapping_windows` on a `MaxEntScan_` field of a row
    whose one-base record lies in the acceptor window of one intron of its transcript and the
    donor window of the next, read from the corpus cache: the three scores where both engines
    write a number, and under `verbose` the two window texts where Ensembl VEP writes nine bases
    and vep-rs twenty-three. A row of the same shape whose record lies in one window or in none,
    on a transcript without such an exon or over several bases, stays an `unexplained_residual`,
    as does every such row where the corpus holds no cache."""
    klass = "maxentscan_overlapping_windows"
    corpus = tmp_path / "c"
    places = {"rec_two": ("1:202", "ENST1"), "rec_two_reverse": ("1:202", "ENST2"), "rec_acceptor_alone": ("1:199", "ENST1"),
              "rec_donor_alone": ("1:205", "ENST1"), "rec_no_window": ("1:150", "ENST1"), "rec_one_intron": ("1:202", "ENST3")}
    _manifest(corpus, [_record(f"1\t{location.split(':')[1]}\t{name}\tA\tG\t.\tPASS\t.") for name, (location, _) in places.items()])
    scores = (("8.100", "2.300", "5.800"), ("6.400", "1.100", "5.300"))
    texts = (_MAXENTSCAN_DONOR_TEXT, _MAXENTSCAN_ACCEPTOR_TEXT)
    names = ["MaxEntScan_ref", "MaxEntScan_alt", "MaxEntScan_diff", "MaxEntScan_ref_seq", "MaxEntScan_alt_seq"]

    def output(side: int) -> str:
        values = [*scores[side], texts[side], texts[side][:-1] + "T"]
        extra = "IMPACT=LOW;STRAND=1" + "".join(f";{name}={value}" for name, value in zip(names, values))
        return "".join(_row(name, location, "G", feature, "splice_region_variant", extra) for name, (location, feature) in places.items())

    vep, rs = tmp_path / "vep.txt", tmp_path / "rs.txt"
    vep.write_text(output(0), encoding="utf-8")
    rs.write_text(output(1), encoding="utf-8")
    command = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
               *(word for name in names for word in ("--field", name))]
    assert bgc.main(command) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert m["field_divergence_summary"] == {"unexplained_residual": 30}
    (corpus / "json_cache" / "transcripts" / "1").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "1" / "1-1000000.json").write_text(json.dumps(_MAXENTSCAN_TRANSCRIPTS), encoding="utf-8")
    assert bgc.main(command) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert m["field_divergence_summary"] == {klass: 10, "unexplained_residual": 20}
    assert {(d["feature"], d["location"]) for d in m["field_divergences"] if d["expected_divergence"] == klass} == {("ENST1", "1:202"), ("ENST2", "1:202")}
    assert {d["field"] for d in m["field_divergences"] if d["expected_divergence"] == klass} == set(names)
    # Under the scores alone the class names the three numbers of each of the two rows.
    assert bgc.main(command[: command.index("--field") + 6]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert m["field_divergence_summary"] == {klass: 6, "unexplained_residual": 12}

    positions = bgc.two_window_positions(corpus / "json_cache")
    assert positions == {("1", "ENST1"): {201, 202, 203}, ("1", "ENST2"): {201, 202, 203}}
    key = ("1:202", "G", "ENST1", "Transcript", "splice_region_variant")

    def named(name: str, vep_value: str, rs_value: str, key: tuple = key, two_windows: dict | None = positions) -> str | None:
        found = bgc.classify_field_divergence(name, key, vep_value, rs_value, "1", [], two_windows=two_windows)
        return found[0] if found else None

    assert named("MaxEntScan_ref", "8.100", "6.400") == klass
    assert named("MaxEntScan_ref", "8.100", "6.400", key=("chr1:202", "G", "ENST1.3", *key[3:])) == klass
    assert named("MaxEntScan_ref", "8.100", "6.400", two_windows=None) is None and named("MaxEntScan_ref", "8.100", "6.400", two_windows={}) is None
    assert named("MaxEntScan_ref", "8.100", "-") is None and named("MaxEntScan_ref", "8.100", "abc") is None
    assert named("MaxEntScan_ref", "8.100", "Inf") is None
    assert named("MaxEntScan_ref", "8.100", "6.400", key=("1:201-203", *key[1:])) is None
    assert named("MaxEntScan_ref", "8.100", "6.400", key=("2:202", *key[1:])) is None
    assert named("MES-SWA_donor_ref", "8.100", "6.400") is None
    assert named("MaxEntScan_ref_seq", *texts) == klass
    assert named("MaxEntScan_ref_seq", texts[1], texts[0]) is None and named("MaxEntScan_alt_seq", texts[0], texts[0][::-1]) is None
    assert named("MaxEntScan_ref_seq", texts[0], "-") is None


def test_two_window_positions_are_the_comparators():
    """The positions the builder reads in two canonical windows of a cache transcript are those
    the format comparator's `SpliceWindows` holds, on either strand, beside an exon of one to six
    bases and beside an intron of two bases."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    for strand in (1, -1):
        for width in range(1, 7):
            exons = [{"start": 1, "end": 100}, {"start": 201, "end": 200 + width}, {"start": 301 + width, "end": 400}]
            tx = {"stable_id": "ENST1", "strand": strand, "exons": exons}
            found = bgc.transcript_two_window_positions(tx)
            assert (len(found) > 0) == (width <= 5), (strand, width)
            windows = cvf.SpliceWindows.__new__(cvf.SpliceWindows)
            windows.input_lines, windows._chromosomes = [], {"1"}
            gaps = [(101, 200), (201 + width, 300 + width)]
            windows._introns = {"1": {"ENST1": (strand, tuple(gaps if strand == 1 else reversed(gaps)))}}
            held = {position for position in range(150, 350)
                    if windows.holds(cvf.Entry(f"1:{position}|G|ENST1|Transcript", None, f"1:{position}", "G", "ENST1", "x"))}
            assert held == found, (strand, width)
    short_intron = {"strand": 1, "exons": [{"start": 1, "end": 100}, {"start": 103, "end": 200}, {"start": 301, "end": 400}]}
    assert bgc.transcript_two_window_positions(short_intron) == set()
    assert bgc.transcript_two_window_positions({"strand": 1, "exons": [{"start": "x", "end": 100}]}) == set()


def test_intron_end_forms_and_gaps_are_the_comparators():
    """The builder's four `HGVSc` forms of `hgvsp_intron_end_insertion`, the pairs that end an
    intron and the bound of a short intron are the comparator's, and the two read one
    transcript's gaps alike, so a row the comparator masks under the class is one the builder
    classifies; a change to either is made in both."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    assert bgc.INTRON_END_COPY_BODY.pattern == cvf.INTRON_END_COPY_RE.pattern
    assert bgc.INTRON_END_DUPLICATION_BODY.pattern == cvf.INTRON_END_DUPLICATION_RE.pattern
    assert bgc.INTRON_END_INSERTION_BODY.pattern == cvf.INTRON_END_INSERTION_RE.pattern
    assert bgc.INTRON_END_PAIRS == cvf.INTRON_END_PAIRS and bgc.SHORT_INTRON_BASES == cvf.SHORT_INTRON_BASES
    for transcript in _INTRON_END_TRANSCRIPTS:
        assert bgc.transcript_intron_end_gaps(transcript) == cvf.transcript_intron_end_gaps(transcript), transcript["stable_id"]
    bodies = ("c.51-1dup", "c.51-3_51-1dup", "c.50+7_51-1dup", "c.51-1_51insT", "c.51-1_51insTAT", "c.51-1_51insTAG",
              "c.51-1_51insTAC", "c.51-1_51-1dup", "c.49+7_51-1dup", "c.-50-1dup", "c.*51-1dup", "c.100dup", "n.51-1dup")
    assert [bgc.intron_end_insertion_form(body) for body in bodies] == [cvf.intron_end_insertion_form(body) for body in bodies]
    assert [bgc.intron_end_insertion_form(body) for body in bodies[:7]] == ["copy", "duplication", "duplication", "insertion", "insertion", "pair", "pair"]
    helpers = ("copy_of_the_intron_last_base", "duplication_of_the_intron_last_bases", "insertion_that_cannot_end_the_intron",
               "insertion_that_ends_in_an_intron_end_pair")
    for helper in helpers:
        assert [getattr(bgc, helper)(body) for body in bodies] == [getattr(cvf, helper)(body) for body in bodies], helper
    assert [body for body in bodies if bgc.insertion_that_ends_in_an_intron_end_pair(body)] == list(bodies[5:7])
    assert [body for body in bodies if bgc.insertion_that_cannot_end_the_intron(body)] == list(bodies[3:5])


def test_classify_names_the_intron_start_insertion_class(tmp_path: Path):
    """An `HGVSp` Ensembl VEP writes and vep-rs leaves empty, on an insertion the 3' rule leaves
    between a coding exon's last base and the first base of the intron after it, is
    `hgvsp_intron_start_insertion`: written in that gap or inside the exon and carried to it by
    one `HGVS_OFFSET` on both rows, on either strand, before an intron of 13 bases and of 2
    bases, and with that offset on a row whose Location names the gap itself, as a run under
    `--shift_genomic 1` prints it. A Location with an offset that is at the gap neither as named
    nor moved, a notation on both rows, a reference without one, an `HGVSc` of another form (a
    duplication, an insertion between two coding positions, two different numbers, a 5' UTR or a
    3' UTR position), two `HGVSc` values that differ, offsets that differ or are zero, a Location
    inside the exon without an offset, at the intron's end, at an exon end outside the coding
    region or over a span beside no gap, and a transcript the cache does not hold stay an
    `unexplained_residual`, as does every row of a run that does not compare `HGVSc` and
    `HGVS_OFFSET` beside `HGVSp`."""
    corpus = tmp_path / "c"
    protein = "ENSP00000000001.1:p.Gly17ArgfsTer5"
    cases = {
        # name: (transcript, the lower position the record is written beside, allele, HGVSc body, class)
        "rec_gap": (1, 11099, "T", "c.50_50+1insT", True),
        "rec_carried": (2, 21097, "TAT", "c.50_50+1insTAT", True),
        "rec_reverse": (21, 211499, "ATA", "c.50_50+1insTAT", True),
        "rec_reverse_carried": (21, 211201, "A", "c.150_150+1insT", True),
        "rec_13_bases": (23, 231099, "T", "c.100_100+1insT", True),
        "rec_2_bases": (23, 231199, "TAT", "c.187_187+1insTAT", True),
        "rec_both_notations": (3, 31099, "TAT", "c.50_50+1insTAT", False),
        "rec_no_reference": (4, 41099, "TAT", "c.50_50+1insTAT", False),
        "rec_duplication": (5, 51099, "G", "c.50dup", False),
        "rec_coding_positions": (6, 61099, "TAT", "c.49_50insTAT", False),
        "rec_two_numbers": (7, 71099, "TAT", "c.50_51+1insTAT", False),
        "rec_utr5": (22, 221099, "TAT", "c.-51_-51+1insTAT", False),
        "rec_utr3": (22, 221499, "TAT", "c.*50_*50+1insTAT", False),
        "rec_two_notations": (8, 81099, "TAT", "c.50_50+1insTAT", False),
        "rec_offsets_differ": (9, 91097, "TAT", "c.50_50+1insTAT", False),
        "rec_zero_offset": (10, 101099, "TAT", "c.50_50+1insTAT", False),
        "rec_inside_exon": (11, 111097, "TAT", "c.50_50+1insTAT", False),
        "rec_intron_end": (12, 121199, "TAT", "c.50_50+1insTAT", False),
        "rec_utr_exon_end": (22, 221099, "T", "c.50_50+1insT", False),
        "rec_absent": (99, 991099, "TAT", "c.50_50+1insTAT", False),
        "rec_wide": (13, 131099, "TAT", "c.50_50+1insTAT", False),
        "rec_located_at_gap": (14, 141099, "TAT", "c.50_50+1insTAT", True),
        "rec_reverse_located_at_gap": (21, 211199, "ATA", "c.150_150+1insTAT", True),
        "rec_fits_neither": (14, 141095, "T", "c.50_50+1insT", False),
    }
    # The fields a case's Ensembl VEP row and vep-rs row write otherwise, and the Location of the
    # case written over three positions.
    changed = {
        "rec_carried": ({"HGVS_OFFSET": "2"}, {"HGVS_OFFSET": "2"}),
        "rec_reverse_carried": ({"HGVS_OFFSET": "-2"}, {"HGVS_OFFSET": "-2"}),
        "rec_both_notations": ({}, {"HGVSp": protein.replace("fsTer5", "fsTer9")}),
        "rec_no_reference": ({"HGVSp": ""}, {"HGVSp": protein}),
        "rec_two_notations": ({}, {"HGVSc": "ENST00000000008.1:c.50_50+1insT"}),
        "rec_offsets_differ": ({"HGVS_OFFSET": "2"}, {"HGVS_OFFSET": "3"}),
        "rec_zero_offset": ({"HGVS_OFFSET": "0"}, {"HGVS_OFFSET": "0"}),
        "rec_located_at_gap": ({"HGVS_OFFSET": "2"}, {"HGVS_OFFSET": "2"}),
        "rec_reverse_located_at_gap": ({"HGVS_OFFSET": "-2"}, {"HGVS_OFFSET": "-2"}),
        "rec_fits_neither": ({"HGVS_OFFSET": "2"}, {"HGVS_OFFSET": "2"}),
    }
    names = list(cases)
    locations = {name: f"1:{position}-{position + 1}" for name, (_, position, *_) in cases.items()}
    locations["rec_wide"] = "1:131101-131103"
    _manifest(corpus, [_record(f"1\t{position}\t{name}\tG\tG{allele}\t.\tPASS\t.") for name, (_, position, allele, *_) in cases.items()])
    (corpus / "json_cache" / "transcripts" / "1").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "1" / "1-1000000.json").write_text(json.dumps(_INTRON_END_TRANSCRIPTS), encoding="utf-8")

    def row(name: str, side: int) -> str:
        number, _, allele, body, _ = cases[name]
        feature = f"ENST{number:011d}"
        fields = {"HGVSc": f"{feature}.1:{body}", "HGVSp": "" if side else protein, "HGVS_OFFSET": "", **changed.get(name, ({}, {}))[side]}
        extra = f"IMPACT=HIGH;STRAND={-1 if number == 21 else 1}" + "".join(f";{field}={value}" for field, value in fields.items() if value)
        return _row(name, locations[name], allele, feature, "frameshift_variant,splice_region_variant", extra)

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, 0) for name in cases), encoding="utf-8")
    rs.write_text("".join(row(name, 1) for name in cases), encoding="utf-8")
    command = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main([*command, "--field", "HGVSc", "--field", "HGVSp", "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    klass = "hgvsp_intron_start_insertion"
    assert all(len(d["record_indices"]) == 1 for d in m["field_divergences"])
    assert {(names[d["record_indices"][0]], d["field"]): d["expected_divergence"] for d in m["field_divergences"]} == {
        **{(name, "HGVSp"): (klass if case[4] else "unexplained_residual") for name, case in cases.items()},
        ("rec_two_notations", "HGVSc"): "unexplained_residual",
        ("rec_offsets_differ", "HGVS_OFFSET"): "unexplained_residual",
    }
    assert m["field_divergence_summary"] == {klass: 8, "unexplained_residual": 18}
    assert {d["vep_rs"] for d in m["field_divergences"] if d["expected_divergence"] == klass} == {"-"}
    # The class reads the row's HGVSc and HGVS_OFFSET from the compared fields: a run that
    # compares HGVSp alone names nothing.
    assert bgc.main([*command, "--field", "HGVSp"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert m["field_divergence_summary"] == {"unexplained_residual": 24}
    gaps = bgc.intron_start_gaps(corpus / "json_cache")
    assert gaps[("1", "ENST00000000001")] == {11099, 11299} and gaps[("1", "ENST00000000021")] == {211199, 211499}
    assert gaps[("1", "ENST00000000022")] == {221299} and gaps[("1", "ENST00000000023")] == {231099, 231199, 231299}
    key = ("1:11097-11098", "TAT", "ENST00000000001", "Transcript", "frameshift_variant,splice_region_variant")
    fields = {"HGVSc": "ENST00000000001.1:c.50_50+1insTAT", "HGVSp": protein, "HGVS_OFFSET": "2"}

    def named(perl_row: dict, rust_row: dict, intron_starts: dict | None = gaps, name: str = "HGVSp", key: tuple = key) -> str | None:
        found = bgc.classify_field_divergence(name, key, perl_row[name], rust_row[name], "1", [], perl_row=perl_row, rust_row=rust_row, intron_starts=intron_starts)
        return found[0] if found else None

    assert named(fields, {**fields, "HGVSp": "-"}) == klass
    assert named(fields, {**fields, "HGVSp": "-"}, key=("chr" + key[0], key[1], key[2] + ".4", *key[3:])) == klass
    assert named(fields, {**fields, "HGVSp": "-"}, key=("1:11099-11100", *key[1:])) == klass
    assert named(fields, {**fields, "HGVSp": "-"}, key=("1:11095-11096", *key[1:])) is None
    assert named(fields, {**fields, "HGVSp": "-"}, None) is None and named(fields, {**fields, "HGVSp": "-"}, {}) is None
    assert named(fields, {**fields, "HGVSp": "-", "HGVS_OFFSET": "-"}) is None
    assert named({**fields, "HGVS_OFFSET": "-"}, {**fields, "HGVSp": "-", "HGVS_OFFSET": "-"}) is None
    assert named({**fields, "HGVS_OFFSET": "-2"}, {**fields, "HGVSp": "-", "HGVS_OFFSET": "-2"}) is None
    assert named(fields, {**fields, "HGVSp": "-"}, name="HGVSc") is None
    without_offset = {field: value for field, value in fields.items() if field != "HGVS_OFFSET"}
    assert named(without_offset, {**without_offset, "HGVSp": "-"}) is None


def test_intron_start_form_and_gaps_are_the_comparators():
    """The builder's `HGVSc` form of `hgvsp_intron_start_insertion` and its reading of a
    transcript's coding exons and of the gaps in which a coding exon ends before an intron are the
    comparator's, so a row the comparator masks under the class is one the builder classifies; a
    change to either is made in both."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    assert bgc.INTRON_START_INSERTION_BODY.pattern == cvf.INTRON_START_INSERTION_RE.pattern
    for transcript in _INTRON_END_TRANSCRIPTS:
        assert bgc.coding_transcript_exons(transcript) == cvf.coding_transcript_exons(transcript), transcript["stable_id"]
        assert bgc.transcript_intron_start_gaps(transcript) == cvf.transcript_intron_start_gaps(transcript), transcript["stable_id"]
    bodies = ("c.50_50+1insT", "c.50_50+1insTAT", "c.50_51+1insTAT", "c.49_50insTAT", "c.50dup", "c.-51_-51+1insTAT",
              "c.*50_*50+1insTAT", "c.50_50+2insTAT", "c.50-1_50insTAT", "n.50_50+1insTAT")
    assert [bgc.insertion_after_the_exon_last_base(body) for body in bodies] == [cvf.insertion_after_the_exon_last_base(body) for body in bodies]
    assert [body for body in bodies if bgc.insertion_after_the_exon_last_base(body)] == list(bodies[:2])


def test_classify_names_the_intron_classes_for_an_insertion_held_with_a_shared_base(tmp_path: Path):
    """A record that carries another allele beside an insertion holds the insertion with the
    reference base it shares, so its row's Location is that base's one position and its Allele the
    base followed by the inserted bases. An `HGVSp` Ensembl VEP writes and vep-rs leaves empty on
    such a row is `hgvsp_intron_start_insertion` where the clipped allele lies after a coding
    exon's last base and `hgvsp_intron_end_insertion` where it lies before an exon's first coding
    base, on a reverse-strand transcript included. A span beside no gap of the class and an Allele
    whose base count does not fit the bases the `HGVSc` inserts stay an `unexplained_residual`."""
    corpus = tmp_path / "c"
    protein = "ENSP00000000001.1:p.Gly17ArgfsTer5"
    start, end = "hgvsp_intron_start_insertion", "hgvsp_intron_end_insertion"
    cases = {
        # name: (transcript, the record's position, allele, HGVSc body, class)
        "rec_start_held": (1, 11099, "GTAT", "c.50_50+1insTAT", start),
        "rec_start_beside_no_gap": (2, 21097, "GTAT", "c.50_50+1insTAT", None),
        "rec_start_count": (3, 31099, "GTAT", "c.50_50+1insT", None),
        "rec_end_held": (21, 211299, "GATA", "c.51-1_51insTAT", end),
        "rec_end_beside_no_gap": (21, 211298, "GATA", "c.51-1_51insTAT", None),
        "rec_end_count": (21, 211099, "GATA", "c.151-1_151insT", None),
    }
    names = list(cases)
    _manifest(corpus, [_record(f"1\t{position}\t{name}\tG\t{allele},C\t.\tPASS\t.") for name, (_, position, allele, *_) in cases.items()])
    (corpus / "json_cache" / "transcripts" / "1").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "1" / "1-1000000.json").write_text(json.dumps(_INTRON_END_TRANSCRIPTS), encoding="utf-8")

    def row(name: str, side: int) -> str:
        number, position, allele, body, _ = cases[name]
        feature = f"ENST{number:011d}"
        extra = f"IMPACT=HIGH;STRAND={-1 if number == 21 else 1};HGVSc={feature}.1:{body}" + ("" if side else f";HGVSp={protein}")
        return _row(name, f"1:{position}", allele, feature, "frameshift_variant,splice_region_variant", extra)

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, 0) for name in cases), encoding="utf-8")
    rs.write_text("".join(row(name, 1) for name in cases), encoding="utf-8")
    command = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main([*command, "--field", "HGVSc", "--field", "HGVSp", "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert all(len(d["record_indices"]) == 1 for d in m["field_divergences"])
    assert {(names[d["record_indices"][0]], d["field"]): d["expected_divergence"] for d in m["field_divergences"]} == {
        (name, "HGVSp"): (case[4] or "unexplained_residual") for name, case in cases.items()
    }
    assert m["field_divergence_summary"] == {start: 1, end: 1, "unexplained_residual": 4}


def test_held_insertion_gaps_are_the_comparators():
    """The builder reads the gaps of an insertion a record holds with reference bases it shares,
    the bases an `ins` notation inserts and a base string on the other strand as the comparator
    does, so a held row the comparator masks under an intron class is one the builder classifies;
    a change to either is made in both."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    held = [
        (100, 1, "GTAT", "TAT", "1"), (100, 1, "TATG", "TAT", "1"), (100, 1, "GATA", "TAT", "-1"),
        (100, 1, "GATA", "TAT", "-"), (100, 1, "GATA", "TAT", "1"), (100, 1, "GGGG", "GGG", "1"),
        (100, 2, "CGTAT", "TAT", "1"), (100, 1, "GTAT", "T", "1"), (100, 1, "TAT", "TAT", "1"), (100, 1, "GTAT", "", "1"),
    ]
    assert [bgc.held_insertion_gaps(*case) for case in held] == [cvf.held_insertion_gaps(*case) for case in held]
    assert [bgc.held_insertion_gaps(*case) for case in held[:6]] == [[100], [99], [100], [100], [], [99, 100]]
    for bases in ("TAT", "GTAT", "GGGG", "CGTAT"):
        assert bgc.reverse_complement(bases) == cvf.reverse_complement(bases)
    bodies = ("c.51-1_51insTAT", "c.50_50+1insT", "c.51-1dup", "c.50dup")
    patterns = ((bgc.INTRON_END_INSERTION_BODY, cvf.INTRON_END_INSERTION_RE), (bgc.INTRON_START_INSERTION_BODY, cvf.INTRON_START_INSERTION_RE))
    for builder_pattern, comparator_pattern in patterns:
        assert [bgc.notation_inserted_bases(builder_pattern, body) for body in bodies] == [cvf.notation_inserted_bases(comparator_pattern, body) for body in bodies]


def test_classify_places_a_row_by_its_location_span_for_the_intron_classes(tmp_path: Path):
    """A record that carries another allele beside the insertion keeps its own Location, one
    position for a one-base REF, and either intron class is named by a gap of that span. Under
    `--minimal` the row's Allele is the insertion alone: the bases an `ins` notation inserts,
    turned by the offset where the 3' rule moved the record, or as many bases as a `dup` notation
    duplicates. A held Allele spelled on the reverse strand is read reverse-complemented. A
    minimised Allele that is not the notation's bases, a `dup` row whose Allele's base count does
    not fit, a duplication that starts in the intron's first half, a moved Allele that holds the
    shared base and a span beside no gap stay an `unexplained_residual`."""
    corpus = tmp_path / "c"
    protein = "ENSP00000000001.1:p.Gly17ArgfsTer5"
    start, end = "hgvsp_intron_start_insertion", "hgvsp_intron_end_insertion"
    cases = {
        # name: (transcript, the Location's position, allele, HGVSc body, class)
        "rec_end_minimised": (23, 231313, "TAT", "c.286-1_286insTAT", end),
        "rec_end_minimised_reverse": (21, 211099, "ATA", "c.151-1_151insTAT", end),
        "rec_end_minimised_copy": (1, 11199, "G", "c.51-1dup", end),
        "rec_end_minimised_duplication": (2, 21199, "CAG", "c.51-3_51-1dup", end),
        "rec_end_reverse_strand_record": (21, 211299, "TATC", "c.51-1_51insTAT", end),
        "rec_start_minimised": (3, 31099, "TAT", "c.50_50+1insTAT", start),
        "rec_start_minimised_past_the_gap": (10, 101100, "TAT", "c.50_50+1insTAT", start),
        "rec_start_minimised_moved": (4, 41097, "TAT", "c.50_50+1insTTA", start),
        "rec_start_reverse_strand_record": (5, 51099, "ATAC", "c.50_50+1insTAT", start),
        "rec_other_bases": (6, 61099, "CCG", "c.50_50+1insTAT", None),
        "rec_count": (7, 71199, "CAG", "c.51-1dup", None),
        "rec_held_moved": (11, 111099, "GTAT", "c.50_50+1insTAT", None),
        "rec_first_half": (24, 241313, "TTCTTCAG", "c.285+7_286-1dup", None),
        "rec_below_the_gap": (8, 81098, "TAT", "c.50_50+1insTAT", None),
        "rec_above_the_gap": (24, 241315, "TAT", "c.286-1_286insTAT", None),
    }
    names = list(cases)
    _manifest(corpus, [_record(f"1\t{position}\t{name}\tG\tG{allele},C\t.\tPASS\t.") for name, (_, position, allele, *_) in cases.items()])
    (corpus / "json_cache" / "transcripts" / "1").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "1" / "1-1000000.json").write_text(json.dumps(_INTRON_END_TRANSCRIPTS), encoding="utf-8")

    def row(name: str, side: int) -> str:
        number, position, allele, body, _ = cases[name]
        feature = f"ENST{number:011d}"
        extra = f"IMPACT=HIGH;STRAND={-1 if number == 21 else 1};HGVSc={feature}.1:{body}" + ("" if side else f";HGVSp={protein}")
        extra += ";HGVS_OFFSET=2" if name in ("rec_start_minimised_moved", "rec_held_moved") else ""
        return _row(name, f"1:{position}", allele, feature, "frameshift_variant,splice_region_variant", extra)

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, 0) for name in cases), encoding="utf-8")
    rs.write_text("".join(row(name, 1) for name in cases), encoding="utf-8")
    command = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main([*command, "--field", "HGVSc", "--field", "HGVSp", "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert all(len(d["record_indices"]) == 1 for d in m["field_divergences"])
    assert {(names[d["record_indices"][0]], d["field"]): d["expected_divergence"] for d in m["field_divergences"]} == {
        (name, "HGVSp"): (case[4] or "unexplained_residual") for name, case in cases.items()
    }
    assert m["field_divergence_summary"] == {end: 5, start: 4, "unexplained_residual": 6}


def test_classify_places_a_row_of_two_adjacent_positions_by_each_gap_of_its_span(tmp_path: Path):
    """A record whose reference allele spans two bases and that carries another allele beside the
    insertion keeps both bases under `--minimal`: the row's Location is their two adjacent
    positions and its Allele the insertion alone. Either intron class is named by a gap of that
    span, its own, as with an insertion written alone, the one before its first position or the one
    after its last, moved by the offset where the 3' rule moved the record. A Location none of
    whose three gaps is a gap of the class, and an Allele that is not the insertion at a Location
    whose own gap is none, stay an `unexplained_residual`."""
    corpus = tmp_path / "c"
    protein = "ENSP00000000001.1:p.Gly17ArgfsTer5"
    start, end = "hgvsp_intron_start_insertion", "hgvsp_intron_end_insertion"
    cases = {
        # name: (transcript, the Location's first position, allele, HGVSc body, HGVS_OFFSET, class)
        "rec_end_after_the_span": (23, 231312, "TAT", "c.286-1_286insTAT", "", end),
        "rec_end_before_the_span": (24, 241314, "TAT", "c.286-1_286insTAT", "", end),
        "rec_start_after_the_span": (1, 11098, "TAT", "c.50_50+1insTAT", "", start),
        "rec_start_before_the_span": (2, 21100, "TAT", "c.50_50+1insTAT", "", start),
        "rec_start_moved": (3, 31096, "TTA", "c.50_50+1insTTA", "2", start),
        "rec_end_no_gap": (23, 231310, "TAT", "c.286-1_286insTAT", "", None),
        "rec_start_no_gap": (4, 41101, "TAT", "c.50_50+1insTAT", "", None),
        "rec_start_other_bases": (5, 51098, "CCG", "c.50_50+1insTAT", "", None),
    }
    names = list(cases)
    _manifest(corpus, [_record(f"1\t{first}\t{name}\tCG\tCG{allele},C\t.\tPASS\t.") for name, (_, first, allele, *_) in cases.items()])
    (corpus / "json_cache" / "transcripts" / "1").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "1" / "1-1000000.json").write_text(json.dumps(_INTRON_END_TRANSCRIPTS), encoding="utf-8")

    def row(name: str, side: int) -> str:
        number, first, allele, body, offset, _ = cases[name]
        feature = f"ENST{number:011d}"
        extra = f"IMPACT=HIGH;STRAND=1;HGVSc={feature}.1:{body}" + ("" if side else f";HGVSp={protein}") + (f";HGVS_OFFSET={offset}" if offset else "")
        return _row(name, f"1:{first}-{first + 1}", allele, feature, "frameshift_variant,splice_region_variant", extra)

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, 0) for name in cases), encoding="utf-8")
    rs.write_text("".join(row(name, 1) for name in cases), encoding="utf-8")
    command = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main([*command, "--field", "HGVSc", "--field", "HGVSp", "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert all(len(d["record_indices"]) == 1 for d in m["field_divergences"])
    assert {(names[d["record_indices"][0]], d["field"]): d["expected_divergence"] for d in m["field_divergences"]} == {
        (name, "HGVSp"): (case[5] or "unexplained_residual") for name, case in cases.items()
    }
    assert m["field_divergence_summary"] == {end: 2, start: 3, "unexplained_residual": 3}


def test_lone_insertion_gaps_are_the_comparators():
    """The builder reads the gaps of a row whose Allele is the insertion alone, and the count of
    positions a `dup` notation duplicates, as the comparator does, so a minimised row the
    comparator masks under an intron class is one the builder classifies; a change to either is
    made in both."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    lone = [
        (100, 1, "TAT", "TAT"), (100, 1, "ATA", "TAT"), (100, 3, "TAT", "TAT"), (100, 1, "CCG", "TAT"),
        (100, 1, "TAT", "TTA", 0, 2), (100, 1, "TAT", "TTA", 0, -5), (100, 1, "TAT", "TAT", 0, 2),
        (100, 1, "ATA", "TAT", 0, -2), (100, 1, "ATT", "TTA", 0, 2),
        (100, 1, "CAG", "", 3), (100, 1, "CAG", "", 1), (100, 1, "CAG", "", 3, 1), (100, 0, "TAT", "TAT"),
        (100, 1, "DELETION", "", 8), (100, 1, "-", "", 1),
    ]
    assert [bgc.lone_insertion_gaps(*case) for case in lone] == [cvf.lone_insertion_gaps(*case) for case in lone]
    assert [bgc.lone_insertion_gaps(*case) for case in lone[:4]] == [[99, 100], [99, 100], [99, 100, 101, 102], []]
    assert [bgc.lone_insertion_gaps(*case) for case in lone[4:9]] == [[99, 100], [99, 100], [99, 100], [99, 100], []]
    bodies = ("c.51-1dup", "c.51-3_51-1dup", "c.50+7_51-1dup", "c.51+3_51-1dup", "c.52-3_51-1dup", "c.51-1_51insTAT", "c.50_50+1insTAT")
    assert [bgc.notation_duplicated_positions(body) for body in bodies] == [cvf.notation_duplicated_positions(body) for body in bodies] == [1, 3, 0, 0, 0, 0, 0]


def test_classify_names_the_intron_start_class_for_an_allele_printed_as_the_rule_moved_it(tmp_path: Path):
    """A run under `--shift_3prime 1` prints a moved record's Allele as the 3' rule moved it, the
    bases the row's `HGVSc` inserts as they stand, beside the offset: a minimised row of a record
    that carries another allele is `hgvsp_intron_start_insertion` with the record's own Location,
    which the offset carries to the gap, and with the moved Location a run under
    `--shift_genomic 1` prints, on a forward-strand transcript with a positive offset and on a
    reverse-strand one with a negative offset. An Allele that is the notation's bases under
    neither reading and a span whose gaps are a gap of the class neither as named nor moved stay
    an `unexplained_residual`."""
    corpus = tmp_path / "c"
    protein = "ENSP00000000001.1:p.Gly17ArgfsTer5"
    klass = "hgvsp_intron_start_insertion"
    cases = {
        # name: (transcript, the Location's position, allele, HGVS_OFFSET, class)
        "rec_forward": (1, 11097, "TTA", "2", True),
        "rec_forward_located_moved": (2, 21099, "TTA", "2", True),
        "rec_reverse": (21, 211501, "TAA", "-2", True),
        "rec_reverse_located_moved": (21, 211499, "TAA", "-2", True),
        "rec_neither_reading": (3, 31097, "ATT", "2", False),
        "rec_neither_place": (4, 41095, "TTA", "2", False),
    }
    names = list(cases)
    _manifest(corpus, [_record(f"1\t{position}\t{name}\tG\tG{allele},C\t.\tPASS\t.") for name, (_, position, allele, *_) in cases.items()])
    (corpus / "json_cache" / "transcripts" / "1").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "1" / "1-1000000.json").write_text(json.dumps(_INTRON_END_TRANSCRIPTS), encoding="utf-8")

    def row(name: str, side: int) -> str:
        number, position, allele, offset, _ = cases[name]
        feature = f"ENST{number:011d}"
        extra = f"IMPACT=HIGH;STRAND={-1 if number == 21 else 1};HGVSc={feature}.1:c.50_50+1insTTA" + ("" if side else f";HGVSp={protein}")
        return _row(name, f"1:{position}", allele, feature, "frameshift_variant,splice_region_variant", f"{extra};HGVS_OFFSET={offset}")

    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(row(name, 0) for name in cases), encoding="utf-8")
    rs.write_text("".join(row(name, 1) for name in cases), encoding="utf-8")
    command = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main([*command, "--field", "HGVSc", "--field", "HGVSp", "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert all(len(d["record_indices"]) == 1 for d in m["field_divergences"])
    assert {(names[d["record_indices"][0]], d["field"]): d["expected_divergence"] for d in m["field_divergences"]} == {
        (name, "HGVSp"): (klass if case[4] else "unexplained_residual") for name, case in cases.items()
    }
    assert m["field_divergence_summary"] == {klass: 4, "unexplained_residual": 2}


def test_classify_names_the_frequency_filter_zero_value_class(tmp_path: Path):
    """Under `--check_frequency` Ensembl VEP drops records by frequency with no warning: a
    record both engines drop is rowless by the run's nature, the keys of one vep-rs keeps
    whose FREQS list a frequency of exactly 0 are `frequency_filter_zero_value`, one vep-rs
    keeps with no zero in its FREQS is refused as a zero-row record without a warning, and
    without the flag the dropped records are refused the same way."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\trs1\tA\tG\t.\tPASS\t."),
        _record("21\t200\trs2\tA\tT\t.\tPASS\t."),
        _record("21\t300\trs3\tA\tG\t.\tPASS\t."),
        _record("21\t400\trs4\tA\tC\t.\tPASS\t."),
    ])
    m = json.loads((corpus / "manifest.json").read_text())
    m["flag_sets"] = [{"name": "include", "flags": ["--check_frequency", "--freq_pop", "gnomADe_AFR", "--freq_freq", "0.05",
                                                    "--freq_gt_lt", "lt", "--freq_filter", "include"]}]
    (corpus / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    common = "IMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE_AFR:{}:{}"
    vep = tmp_path / "vep.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "intron_variant", common.format("G", "0.01")), encoding="utf-8")
    kept = (_row("rs1", "21:100", "G", "ENST1", "intron_variant", common.format("G", "0.01"))
            + _row("rs3", "21:300", "G", "ENST3", "intron_variant", common.format("G", "0"))
            + _row("rs3", "21:300", "G", "ENST4", "upstream_gene_variant", "IMPACT=MODIFIER;DISTANCE=12;STRAND=1;FREQS=GNOMADE_AFR:G:0"))
    rs = tmp_path / "rs.txt"
    rs.write_text(kept, encoding="utf-8")
    args = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main([*args, "--flag-set", "include"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    fs = m["flag_sets"][0]
    assert [(t["feature"], t["record_indices"], t["expected_divergence"]) for t in fs["vep_rs_only_tuples"]] == [
        ("ENST3", [2], "frequency_filter_zero_value"),
        ("ENST4", [2], "frequency_filter_zero_value"),
    ]
    assert fs["vep_rs_only_summary"] == {"frequency_filter_zero_value": 2}
    assert fs["divergences"] == [] and not any("reference_rows" in r for r in m["records"])
    assert "vep_rs_only_tuples" not in m

    rs.write_text(kept + _row("rs4", "21:400", "C", "ENST5", "intron_variant", common.format("C", "0.3")), encoding="utf-8")
    assert bgc.main([*args, "--flag-set", "include"]) == 1
    rs.write_text(kept, encoding="utf-8")
    assert bgc.main(args) == 1


def test_lists_a_zero_frequency_reads_the_last_field_of_each_item():
    assert bgc.lists_a_zero_frequency("GNOMADE_AFR:G:0")
    assert bgc.lists_a_zero_frequency("AF:A:0.5625,AF:C:0.0")
    assert not bgc.lists_a_zero_frequency("AF:A:NA,AF:G:1")
    assert not bgc.lists_a_zero_frequency("-")
    assert not bgc.lists_a_zero_frequency("AF:-:6.365e-06")


def test_classify_names_the_organelle_peptide_path_class(tmp_path: Path):
    """An `HGVSp` difference on a transcript the corpus cache files under a codon table other
    than the standard one is `organelle_peptide_paths_translation_table` when the row's
    consequence set and the descriptions place it in one of the four bare-translation paths:
    any frameshift row (a stop distance, a `Ter` at an immediate stop, a lost-start shape),
    an extension count on a `stop_lost` row (the deletion of an empty reference window
    among them), a duplication read as an insertion. The same pairs on a standard-code
    transcript, a substitution, a `delins` and a plain deletion on the organelle transcript
    (their residues come through `peptide`, table-aware on both engines) and a pair naming
    two proteins stay `unexplained_residual`."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("MT\t3313\tfs\tAT\tA\t.\tPASS\t."),
        _record("MT\t3399\tdup\tA\tAATA\t.\tPASS\t."),
        _record("MT\t7443\text\tA\tG\t.\tPASS\t."),
        _record("MT\t3308\tsub\tT\tC\t.\tPASS\t."),
        _record("MT\t3400\tname\tC\tCG\t.\tPASS\t."),
        _record("MT\t4260\tstop\tCG\tC\t.\tPASS\t."),
        _record("MT\t4470\tlost\tGT\tG\t.\tPASS\t."),
        _record("MT\t5000\tdelins\tTGCTA\tTA\t.\tPASS\t."),
        _record("21\t100\tnuclear\tAT\tA\t.\tPASS\t."),
    ])
    shard = corpus / "json_cache" / "transcripts" / "MT"
    shard.mkdir(parents=True)
    (shard / "1-1000000.json").write_text(json.dumps([
        {"stable_id": "ENSTMT1", "variation_effect_feature_cache": {"codon_table": 2}},
        {"stable_id": "ENSTMT2", "variation_effect_feature_cache": {"codon_table": "2"}},
    ]), encoding="utf-8")
    (corpus / "json_cache" / "transcripts" / "21").mkdir()
    (corpus / "json_cache" / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([
        {"stable_id": "ENST1", "variation_effect_feature_cache": {"codon_table": 1}},
    ]), encoding="utf-8")
    hgvsp = "IMPACT=HIGH;STRAND=1;HGVSp={}"
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(
        _row("fs", "MT:3314", "-", "ENSTMT1", "frameshift_variant", hgvsp.format("ENSPMT1.1:p.Met3ArgfsTer12"))
        + _row("dup", "MT:3399-3400", "ATA", "ENSTMT1", "inframe_insertion", hgvsp.format("ENSPMT1.1:p.Met31_Gln32insMet"))
        + _row("ext", "MT:7443", "G", "ENSTMT2", "stop_lost", hgvsp.format("ENSPMT2.1:p.Ter514GlyextTer1"))
        + _row("sub", "MT:3308", "C", "ENSTMT1", "missense_variant", hgvsp.format("ENSPMT1.1:p.Met1Thr"))
        + _row("name", "MT:3400-3401", "G", "ENSTMT1", "frameshift_variant", hgvsp.format("ENSPMT1.1:p.Gln32ArgfsTer20"))
        + _row("stop", "MT:4261", "-", "ENSTMT1", "frameshift_variant", hgvsp.format("ENSPMT1.1:p.Asp3Ter"))
        + _row("lost", "MT:4471", "-", "ENSTMT2", "frameshift_variant,start_lost", hgvsp.format("ENSPMT2.1:p.Thr2_?1"))
        + _row("delins", "MT:5001-5003", "-", "ENSTMT1", "inframe_deletion", hgvsp.format("ENSPMT1.1:p.Ala4delinsGly"))
        + _row("nuclear", "21:101", "-", "ENST1", "frameshift_variant", hgvsp.format("ENSP1.1:p.Leu5ArgfsTer12")),
        encoding="utf-8",
    )
    rs.write_text(
        _row("fs", "MT:3314", "-", "ENSTMT1", "frameshift_variant", hgvsp.format("ENSPMT1.1:p.Met3Ter"))
        + _row("dup", "MT:3399-3400", "ATA", "ENSTMT1", "inframe_insertion", hgvsp.format("ENSPMT1.1:p.Met31dup"))
        + _row("ext", "MT:7443", "G", "ENSTMT2", "stop_lost", hgvsp.format("ENSPMT2.1:p.Ter514GlyextTer3"))
        + _row("sub", "MT:3308", "C", "ENSTMT1", "missense_variant", hgvsp.format("ENSPMT1.1:p.Met1Ile"))
        + _row("name", "MT:3400-3401", "G", "ENSTMT1", "frameshift_variant", hgvsp.format("ENSPMT1.2:p.Gln32ArgfsTer5"))
        + _row("stop", "MT:4261", "-", "ENSTMT1", "frameshift_variant", hgvsp.format("ENSPMT1.1:p.Arg2Ter"))
        + _row("lost", "MT:4471", "-", "ENSTMT2", "frameshift_variant,start_lost", hgvsp.format("ENSPMT2.1:p.Lys3_?1"))
        + _row("delins", "MT:5001-5003", "-", "ENSTMT1", "inframe_deletion", hgvsp.format("ENSPMT1.1:p.Ala4delinsSer"))
        + _row("nuclear", "21:101", "-", "ENST1", "frameshift_variant", hgvsp.format("ENSP1.1:p.Leu5Ter")),
        encoding="utf-8",
    )
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "HGVSp"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {d["location"]: d["expected_divergence"] for d in m["field_divergences"]}
    organelle = "organelle_peptide_paths_translation_table"
    assert classes == {"MT:3314": organelle, "MT:3399-3400": organelle, "MT:7443": organelle,
                       "MT:4261": organelle, "MT:4471": organelle,
                       "MT:3308": "unexplained_residual", "MT:3400-3401": "unexplained_residual",
                       "MT:5001-5003": "unexplained_residual", "21:101": "unexplained_residual"}
    assert m["field_divergence_summary"] == {organelle: 5, "unexplained_residual": 4}
    assert bgc.organelle_codon_tables(corpus) == {"ENSTMT1": 2, "ENSTMT2": 2}
    shape = bgc.peptide_description_shape
    assert [shape(d) for d in ("Met31dup", "Met31_Gln32insMet", "Met31_Gln32dup", "Ser513del",
                               "Ala4delinsGly", "Met3ArgfsTer12", "Met1?")] == [
        "dup", "ins", "dup", "other", "other", "other", "other"]
    assert not bgc.organelle_peptide_path_divergence("stop_lost", "Ter514Gly", "Ter514Trp")
    assert bgc.organelle_peptide_path_divergence("stop_lost", "Ter4delextTer?", "Ter4delextTer2")
    assert bgc.organelle_peptide_path_divergence("protein_altering_variant", "Met31dup", "Met31_Gln32insMet")
    assert not bgc.organelle_peptide_path_divergence("inframe_deletion", "Ser5_Lys6del", "Ser5del")


# The transcript the Downstream plugin is tested on in crates/vep-builtins/src/plugins/downstream.rs:
# the codons ATG AAA ACT CTT GCT GTC CGT GGT TAA (`MKTLAVRG*`) and a 3' UTR that reads GGG CCC AAA
# TAG (`GPK*`) from its first base.
DOWNSTREAM_CDS = "ATGAAAACTCTTGCTGTCCGTGGTTAA"
DOWNSTREAM_UTR = "GGGCCCAAATAGCCC"
# A second transcript, ATG TAC GGA AAA CTG TAA (`MYGKL*`), whose second codon becomes a stop when
# its third base goes.
DOWNSTREAM_STOP_CDS = "ATGTACGGAAAACTGTAA"
# A third transcript, ATG AAA CTG ACC TAG (`MKLT*`), filed under codon table 2 with a 3' UTR that
# opens with A: when the third base of its second codon goes, the read crosses TGA, tryptophan
# under that table and a stop under the standard one, and ends at AGA, a stop under table 2 alone.
DOWNSTREAM_TABLE_TWO_CDS = "ATGAAACTGACCTAG"
# The Extra column of a frameshift row the Downstream plugin annotated: the row's STRAND and the
# plugin's two fields.
DOWNSTREAM_EXTRA = "IMPACT=HIGH;STRAND={};DownstreamProtein={};ProteinLengthChange={}"
START_PAST_OWN_CODON = "downstream_start_past_own_codon"
UTR_FRAME_PAST_CODING_END = "downstream_utr_frame_past_coding_end"


def _coding_row(uv, loc, allele, feat, cons, cds, extra):
    """`_row` with the row's `CDS_position` column set to `cds`."""
    columns = _row(uv, loc, allele, feat, cons, extra).rstrip("\n").split("\t")
    columns[8] = cds
    return "\t".join(columns) + "\n"


def _downstream_field_classes(tmp_path: Path, rows: dict[str, tuple[str, str, str, str]],
                              alleles: dict[str, str] | None = None,
                              terms: dict[str, str] | None = None) -> tuple[dict[tuple[str, str], str], dict[str, int]]:
    """`classify` on `DownstreamProtein`, `ProteinLengthChange` and `HGVS_OFFSET` over one record
    per row, each a (transcript, `CDS_position`, Ensembl VEP Extra, vep-rs Extra) under its
    record's name, its allele `-` and its term `frameshift_variant` unless `alleles` or `terms`
    names another. The corpus cache holds the first transcript as `ENST1`, as `ENST4` and, under
    codon table 7, as `ENST5`, the second transcript as `ENST2` and the third, under codon table
    2, as `ENST3`. Returns the class of each field divergence by (record name, field) and the
    run's `field_divergence_summary`."""
    corpus = tmp_path / "c"
    names = list(rows)
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\t{name}\tAT\tA\t.\tPASS\t.") for i, name in enumerate(names)])
    shard = corpus / "json_cache" / "transcripts" / "21"
    shard.mkdir(parents=True)
    first = {"translateable_seq": DOWNSTREAM_CDS, "three_prime_utr": DOWNSTREAM_UTR}
    (shard / "1-1000000.json").write_text(json.dumps([
        {"stable_id": "ENST1", "variation_effect_feature_cache": {**first, "codon_table": 1}},
        {"stable_id": "ENST2", "variation_effect_feature_cache": {**first, "translateable_seq": DOWNSTREAM_STOP_CDS}},
        {"stable_id": "ENST3", "variation_effect_feature_cache": {
            "translateable_seq": DOWNSTREAM_TABLE_TWO_CDS, "three_prime_utr": "ACC", "codon_table": "2"}},
        {"stable_id": "ENST4", "variation_effect_feature_cache": first},
        {"stable_id": "ENST5", "variation_effect_feature_cache": {**first, "codon_table": 7}},
    ]), encoding="utf-8")
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    for path, engine in ((vep, 2), (rs, 3)):
        path.write_text("".join(
            _coding_row(name, f"21:{100 * (i + 1) + 1}", (alleles or {}).get(name, "-"), rows[name][0],
                        (terms or {}).get(name, "frameshift_variant"), rows[name][1], rows[name][engine])
            for i, name in enumerate(names)), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "DownstreamProtein", "--field", "ProteinLengthChange", "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {(names[d["record_indices"][0]], d["field"]): d["expected_divergence"] for d in m["field_divergences"]}
    assert len(classes) == len(m["field_divergences"])
    return classes, m["field_divergence_summary"]


def test_codon_table_amino_acids_are_the_codon_modules_genetic_codes():
    """`CODON_TABLE_AMINO_ACIDS` carries, id for id, the amino-acid string of every table `CODES`
    of crates/vep-core/src/codon.rs holds, and none for an id that list leaves empty (0, 7, 8 and
    17 to 20)."""
    source = (HERE.parent.parent / "crates" / "vep-core" / "src" / "codon.rs").read_text(encoding="utf-8")
    start = source.index("const CODES:")
    block = source[start:source.index("];", start)]
    entries = re.findall(r'^    (?:None|code\(\s*b"([^"]{64})")', block, flags=re.MULTILINE)
    assert len(entries) == 25
    assert {table: amino_acids for table, amino_acids in enumerate(entries) if amino_acids} == bgc.CODON_TABLE_AMINO_ACIDS
    assert sorted(bgc.GENETIC_CODES) == [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24]
    assert all(len(code) == 64 for code in bgc.GENETIC_CODES.values())


def test_translate_codons_reads_codons_as_the_codon_table_module_does():
    """Whole codons are read in order under the named table and a trailing one or two bases are
    dropped; a stop stays `*`; `U` reads as `T` and lower case as upper; a codon with an ambiguous
    letter is the one amino acid all its readings code under that table, `B` for aspartate or
    asparagine, `Z` for glutamate or glutamine and `X` otherwise; `---` is a gap, and a codon with
    any other character that is no nucleotide letter is `X`. Table 2 reads TGA as tryptophan, AGA
    and AGG as stops and ATA as methionine, where the standard table reads a stop, arginine twice
    and isoleucine."""
    translate = bgc.translate_codons
    assert translate(DOWNSTREAM_CDS, 1) == "MKTLAVRG*"
    assert translate(DOWNSTREAM_STOP_CDS, 1) == "MYGKL*"
    assert translate(DOWNSTREAM_TABLE_TWO_CDS, 1) == translate(DOWNSTREAM_TABLE_TWO_CDS, 2) == "MKLT*"
    assert translate("ATGA", 1) == translate("ATGAA", 1) == "M" and translate("AT", 1) == ""
    assert translate("augCTGtaa", 1) == "ML*"
    assert [translate(codon, 1) for codon in ("CTN", "TAR", "TRA", "RAY", "SAA", "NAT", "NNN")] == ["L", "*", "*", "B", "Z", "X", "X"]
    assert [translate(codon, 1) for codon in ("---", "A-G", "A.G")] == ["-", "X", "X"]
    codons = ("TGA", "AGA", "AGG", "ATA", "TGR", "AGR")
    assert [translate(codon, 1) for codon in codons] == ["*", "R", "R", "I", "X", "R"]
    assert [translate(codon, 2) for codon in codons] == ["W", "*", "*", "M", "W", "*"]


def test_downstream_sequences_reads_the_coding_transcripts_of_the_corpus_cache(tmp_path: Path):
    """Each transcript carrying a `translateable_seq` comes with that sequence, its
    `three_prime_utr` (the empty sequence when the cache spells none) and its `codon_table`, a
    number or the converted cache's string, 1 when the cache names none or 0; a transcript with
    no coding sequence, or whose table is no number, is left out."""
    corpus = tmp_path / "c"
    _manifest(corpus, [])
    shard = corpus / "json_cache" / "transcripts" / "21"
    shard.mkdir(parents=True)
    coding = {"translateable_seq": DOWNSTREAM_CDS}
    (shard / "1-1000000.json").write_text(json.dumps([
        {"stable_id": "ENST1", "variation_effect_feature_cache": {**coding, "three_prime_utr": DOWNSTREAM_UTR, "codon_table": 1}},
        {"stable_id": "ENST2", "variation_effect_feature_cache": {**coding, "three_prime_utr": None}},
        {"stable_id": "ENST3", "variation_effect_feature_cache": {**coding, "codon_table": "2"}},
        {"stable_id": "ENST4", "variation_effect_feature_cache": {**coding, "codon_table": 0}},
        {"stable_id": "ENST5", "variation_effect_feature_cache": {**coding, "codon_table": "0"}},
        {"stable_id": "ENST6", "variation_effect_feature_cache": {**coding, "codon_table": "vertebrate"}},
        {"stable_id": "ENST7", "variation_effect_feature_cache": {"three_prime_utr": DOWNSTREAM_UTR}},
        {"stable_id": "ENST8"},
    ]), encoding="utf-8")
    assert bgc.downstream_sequences(corpus) == {
        "ENST1": (DOWNSTREAM_CDS, DOWNSTREAM_UTR, 1), "ENST2": (DOWNSTREAM_CDS, "", 1), "ENST3": (DOWNSTREAM_CDS, "", 2),
        "ENST4": (DOWNSTREAM_CDS, "", 1), "ENST5": (DOWNSTREAM_CDS, "", 1),
    }


def test_classify_names_the_downstream_class_of_a_deletion_cut_to_the_last_coding_base(tmp_path: Path):
    """A deletion with no inserted base from a codon's third base to the last coding base is
    `downstream_utr_frame_past_coding_end`. Deleting bases 24 to 27 leaves 23 bases ending GG:
    Ensembl VEP's boundary, base 24, lies past them, so it reads the UTR alone, GGG CCC AAA TAG
    (`GPK`), and vep-rs reads GG with the UTR, GGG GGC CCA AAT AGC (`GGPNS`). Deleting base 27
    leaves the stop codon's TA, which takes the UTR's first base: TAG, a stop, so vep-rs writes
    the empty value where Ensembl VEP reads the UTR alone."""
    extra = DOWNSTREAM_EXTRA.format
    classes, summary = _downstream_field_classes(tmp_path, {
        "to_the_end": ("ENST1", "24-27", extra(1, "GPK", 4), extra(1, "GGPNS", 4)),
        "last_base": ("ENST1", "27", extra(1, "GPK", 0), extra(1, "", 0)),
    })
    assert classes == {("to_the_end", "DownstreamProtein"): UTR_FRAME_PAST_CODING_END,
                       ("last_base", "DownstreamProtein"): UTR_FRAME_PAST_CODING_END}
    assert summary == {UTR_FRAME_PAST_CODING_END: 2}


def test_classify_names_the_downstream_class_of_a_deletion_from_a_codons_third_base(tmp_path: Path):
    """A deletion whose first base is a codon's third and that leaves coding bases after it is
    `downstream_start_past_own_codon`: vep-rs starts at the changed codon and Ensembl VEP at the
    codon after it. Deleting base 6 keeps the residue, the codon AAA taking the next A: vep-rs
    reads AAA CTC TTG CTG TCC GTG GTT AAG GGC CCA AAT AGC (`KLLLSVVKGPNS`) and Ensembl VEP starts
    at CTC. Deleting bases 6 and 7 changes it, AAC then TCT TGC TGT CCG TGG TTA AGG GCC CAA ATA
    GCC (`NSCCPWLRAQIA`), and Ensembl VEP starts at TCT. Deleting base 6 of the second transcript
    turns TAC into TAG, a stop, so vep-rs writes the empty value and Ensembl VEP reads past the
    stop, GAA AAC TGT AAG GGC CCA AAT AGC (`ENCKGPNS`)."""
    extra = DOWNSTREAM_EXTRA.format
    classes, summary = _downstream_field_classes(tmp_path, {
        "kept_residue": ("ENST1", "6", extra(1, "LLLSVVKGPNS", 5), extra(1, "KLLLSVVKGPNS", 5)),
        "changed_residue": ("ENST1", "6-7", extra(1, "SCCPWLRAQIA", 5), extra(1, "NSCCPWLRAQIA", 5)),
        "stop": ("ENST2", "6", extra(1, "ENCKGPNS", -4), extra(1, "", -4)),
    })
    assert classes == {("kept_residue", "DownstreamProtein"): START_PAST_OWN_CODON,
                       ("changed_residue", "DownstreamProtein"): START_PAST_OWN_CODON,
                       ("stop", "DownstreamProtein"): START_PAST_OWN_CODON}
    assert summary == {START_PAST_OWN_CODON: 3}


def test_classify_reads_a_downstream_class_on_the_printed_span_alone(tmp_path: Path):
    """The two classes of a `DownstreamProtein` are read on the `CDS_position` the row prints,
    whatever its `HGVS_OFFSET`, since vep-rs computes on that span. Bases 11 and 12 are both T:
    a row printing base 11, a codon's second, carries the pair of a deletion of base 12, a
    codon's third, where vep-rs would read CTG CTG TCC GTG GTT AAG GGC CCA AAT AGC
    (`LLSVVKGPNS`) and Ensembl VEP start one codon later. No rule gives that pair on base 11,
    with an offset of 1 or 2 or none, on either strand, so on a transcript whose cache has no
    mapper pairs to move the span by the row stays `unexplained_residual`."""
    extra = DOWNSTREAM_EXTRA.format
    pair = ("LSVVKGPNS", "LLSVVKGPNS")
    classes, summary = _downstream_field_classes(tmp_path, {
        "shifted": ("ENST1", "11", *(extra(1, value, 5) + ";HGVS_OFFSET=1" for value in pair)),
        "shifted_reverse": ("ENST4", "11", *(extra(-1, value, 5) + ";HGVS_OFFSET=-1" for value in pair)),
        "unshifted": ("ENST1", "11", *(extra(1, value, 5) for value in pair)),
        "shifted_further": ("ENST1", "11", *(extra(1, value, 5) + ";HGVS_OFFSET=2" for value in pair)),
    })
    assert classes == {(name, "DownstreamProtein"): "unexplained_residual"
                       for name in ("shifted", "shifted_reverse", "unshifted", "shifted_further")}
    assert summary == {"unexplained_residual": 4}


def test_classify_reads_a_downstream_allele_on_the_transcripts_strand(tmp_path: Path):
    """A deletion with inserted bases on a reverse-strand transcript puts the reverse complement
    of the row's allele into the coding sequence: the allele AC replaces bases 6 to 8 with GT,
    so vep-rs reads AAG TTC TTG CTG TCC GTG GTT AAG GGC CCA AAT AGC (`KFLLSVVKGPNS`) and Ensembl
    VEP starts at TTC. The same row on a forward-strand transcript, where the allele goes in as
    written, gives another pair and stays `unexplained_residual`."""
    extra = DOWNSTREAM_EXTRA.format
    classes, summary = _downstream_field_classes(tmp_path, {
        "reverse": ("ENST4", "6-8", extra(-1, "FLLSVVKGPNS", 5), extra(-1, "KFLLSVVKGPNS", 5)),
        "forward": ("ENST1", "6-8", extra(1, "FLLSVVKGPNS", 5), extra(1, "KFLLSVVKGPNS", 5)),
    }, alleles={"reverse": "AC", "forward": "AC"})
    assert classes == {("reverse", "DownstreamProtein"): START_PAST_OWN_CODON,
                       ("forward", "DownstreamProtein"): "unexplained_residual"}
    assert summary == {START_PAST_OWN_CODON: 1, "unexplained_residual": 1}


def test_classify_leaves_a_downstream_pair_the_sequences_do_not_give_unexplained(tmp_path: Path):
    """The pair of a deletion of base 6 stays `unexplained_residual` on a row whose span opens on
    a codon's first or second base, has a gap for an end or runs past the coding sequence; with
    either value changed or absent; on a row without `frameshift_variant` or with a splice term;
    and on a transcript the cache lacks. A `ProteinLengthChange` that differs is in neither class,
    on a row whose `DownstreamProtein` is."""
    extra = DOWNSTREAM_EXTRA.format
    pair = (extra(1, "LLLSVVKGPNS", 5), extra(1, "KLLLSVVKGPNS", 5))
    absent = "IMPACT=HIGH;STRAND=1"
    classes, summary = _downstream_field_classes(tmp_path, {
        "first_base": ("ENST1", "4", *pair),
        "second_base": ("ENST1", "5", *pair),
        "gap": ("ENST1", "6-?", *pair),
        "past_the_end": ("ENST1", "24-28", extra(1, "GPK", 4), extra(1, "GGPNS", 4)),
        "other_vep_rs_value": ("ENST1", "6", pair[0], extra(1, "KLLLSVVKGPN", 5)),
        "other_vep_value": ("ENST1", "6", extra(1, "LLSVVKGPNS", 5), pair[1]),
        "no_vep_rs_value": ("ENST1", "6", pair[0], absent),
        "no_vep_value": ("ENST1", "6", absent, pair[1]),
        "inframe": ("ENST1", "6", *pair),
        "splice": ("ENST1", "6", *pair),
        "uncached": ("ENST9", "6", *pair),
        "length": ("ENST1", "6", pair[0], extra(1, "KLLLSVVKGPNS", 6)),
    }, terms={"inframe": "inframe_deletion", "splice": "frameshift_variant,splice_region_variant"})
    unexplained = ["first_base", "second_base", "gap", "past_the_end", "other_vep_rs_value", "other_vep_value",
                   "no_vep_rs_value", "no_vep_value", "inframe", "splice", "uncached"]
    assert classes == {
        **{(name, "DownstreamProtein"): "unexplained_residual" for name in unexplained},
        ("no_vep_rs_value", "ProteinLengthChange"): "unexplained_residual",
        ("no_vep_value", "ProteinLengthChange"): "unexplained_residual",
        ("length", "DownstreamProtein"): START_PAST_OWN_CODON,
        ("length", "ProteinLengthChange"): "unexplained_residual",
    }
    assert summary == {"unexplained_residual": 14, START_PAST_OWN_CODON: 1}


def test_classify_reads_a_downstream_pair_under_the_transcripts_codon_table(tmp_path: Path):
    """Both values are read under the transcript's codon table. Deleting base 6 of the third
    transcript, filed under table 2, leaves AAC TGA CCT AG and the UTR's A: under table 2 TGA is
    tryptophan and AGA a stop, so vep-rs reads `NWP` and Ensembl VEP, one codon later, `WP`, and
    the row is `downstream_start_past_own_codon`. Under the standard table TGA is a stop, which
    gives vep-rs `N` and Ensembl VEP the empty value: the row carrying that pair stays
    `unexplained_residual`. So does the first transcript's pair on a transcript filed under table
    7, an id Bio::Tools::CodonTable holds no code for."""
    extra = DOWNSTREAM_EXTRA.format
    classes, summary = _downstream_field_classes(tmp_path, {
        "table_two": ("ENST3", "6", extra(1, "WP", 0), extra(1, "NWP", 0)),
        "standard_table_values": ("ENST3", "6", extra(1, "", -2), extra(1, "N", -2)),
        "table_seven": ("ENST5", "6", extra(1, "LLLSVVKGPNS", 5), extra(1, "KLLLSVVKGPNS", 5)),
    })
    assert classes == {("table_two", "DownstreamProtein"): START_PAST_OWN_CODON,
                       ("standard_table_values", "DownstreamProtein"): "unexplained_residual",
                       ("table_seven", "DownstreamProtein"): "unexplained_residual"}
    assert summary == {START_PAST_OWN_CODON: 1, "unexplained_residual": 2}


HGVS_SHIFTED_CDS_COORDINATES = "plugin_hgvs_shifted_cds_coordinates"
# Mapper pairs (first and last cDNA base, first and last genomic base, orientation) that lay a
# transcript of three 5' UTR bases, a coding sequence and DOWNSTREAM_UTR on the genome from base
# 1001. On one exon coding base N is genomic base 1003 + N; on the reverse strand the first
# transcript's coding base N is genomic base 1043 - N.
ONE_EXON = [(1, 45, 1001, 1045, 1)]
REVERSE_EXON = [(1, 45, 1001, 1045, -1)]
# Two exons around a two-base intron after coding base 12: coding base 13 is genomic base 1018.
TWO_EXONS = [(1, 15, 1001, 1015, 1), (16, 45, 1018, 1047, 1)]
# Two pairs that touch on the genome between bases 1015 and 1016 while the model carries coding
# base 13 between them, and two that touch between bases 1014 and 1015 around coding bases 12
# and 13.
ONE_BASE_SEAM = [(1, 15, 1001, 1015, 1), (17, 45, 1016, 1044, 1)]
TWO_BASE_SEAM = [(1, 14, 1001, 1014, 1), (17, 45, 1015, 1043, 1)]
# The second transcript on one exon: its 18 coding bases end at cDNA base 21.
STOP_EXON = [(1, 36, 1001, 1036, 1)]


def _mapped_downstream_classes(tmp_path: Path, rows: dict[str, tuple],
                               coding: dict[str, str] | None = None,
                               without: dict[str, tuple[str, ...]] | None = None) -> tuple[dict[tuple[str, str], str], dict[str, int]]:
    """`classify` on `DownstreamProtein`, `ProteinLengthChange` and `HGVS_OFFSET` over one record
    per row, each a (mapper pairs, POS, REF, ALT, `CDS_position`, Ensembl VEP Extra, vep-rs Extra)
    under its record's name. REF and ALT share their first base, and the row's Location and
    Allele are the ones Ensembl VEP prints for such a record. Every row has a transcript of its
    own in the corpus cache, `ENST1` for the first row, laid on the genome by the row's pairs
    (`None` leaves it out of the cache): three 5' UTR bases, the first transcript's coding
    sequence, or the one `coding` names for the row, with its peptide and the 3' UTR, less the
    cache keys `without` names for the row. Returns as `_downstream_field_classes` does."""
    corpus = tmp_path / "c"
    names = list(rows)
    _manifest(corpus, [_record(f"21\t{rows[name][1]}\t{name}\t{rows[name][2]}\t{rows[name][3]}\t.\tPASS\t.") for name in names])
    transcripts = []
    outputs: tuple[list[str], list[str]] = ([], [])
    for i, name in enumerate(names):
        pairs, pos, ref, alt, cds, *extras = rows[name]
        sequence = (coding or {}).get(name, DOWNSTREAM_CDS)
        if pairs is not None:
            cache = {"translateable_seq": sequence, "three_prime_utr": DOWNSTREAM_UTR, "codon_table": "1",
                     "peptide": bgc.translate_codons(sequence, 1).rstrip("*"),
                     "mapper": {"pairs": [dict(zip(bgc.MAPPER_PAIR_KEYS, map(str, pair))) for pair in pairs]}}
            for key in (without or {}).get(name, ()):
                del cache[key]
            transcripts.append({"stable_id": f"ENST{i + 1}", "cdna_coding_start": 4, "cdna_coding_end": 3 + len(sequence),
                                "variation_effect_feature_cache": cache})
        low, high = sorted((pos + 1, pos + len(ref) - 1))
        location = f"21:{low}" if low == high else f"21:{low}-{high}"
        for output, extra in zip(outputs, extras):
            output.append(_coding_row(name, location, alt[1:] or "-", f"ENST{i + 1}", "frameshift_variant", cds, extra))
    shard = corpus / "json_cache" / "transcripts" / "21"
    shard.mkdir(parents=True)
    (shard / "1-1000000.json").write_text(json.dumps(transcripts), encoding="utf-8")
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    for path, output in zip((vep, rs), outputs):
        path.write_text("".join(output), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "DownstreamProtein", "--field", "ProteinLengthChange", "--field", "HGVS_OFFSET"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {(names[d["record_indices"][0]], d["field"]): d["expected_divergence"] for d in m["field_divergences"]}
    assert len(classes) == len(m["field_divergences"])
    return classes, m["field_divergence_summary"]


def _shifted(strand: int, offset: int, vep: tuple[str, int], vep_rs: tuple[str, int]) -> tuple[str, str]:
    """The two Extra columns of a frameshift row whose `HGVS_OFFSET` both engines print: Ensembl
    VEP's with its (`DownstreamProtein`, `ProteinLengthChange`), then vep-rs's."""
    return tuple(DOWNSTREAM_EXTRA.format(strand, *values) + f";HGVS_OFFSET={offset}" for values in (vep, vep_rs))


def test_record_span_is_the_span_the_vcf_parser_gives_a_record():
    """A record's start is POS and its end the last REF base; where an ALT differs from REF in
    length and REF and the ALTs open with one base, the start moves one base on
    (Parser/VCF.pm:295-336), so an insertion's start lies one past its end. An ALT naming `*`
    is left out of the comparison among several alleles; one ALT opening with another base, or
    a record whose alleles all have REF's length, keeps POS. `row_record_span` picks the one
    record span a Location prints, an insertion's two flanks lower first."""
    span = bgc.record_span
    assert span("21\t100\t.\tAT\tA\t.\t.\t.") == (101, 101)
    assert span("21\t100\t.\tATTT\tA\t.\t.\t.") == (101, 103)
    assert span("21\t100\t.\tA\tATT\t.\t.\t.") == (101, 100)
    assert span("21\t100\t.\tA\tG\t.\t.\t.") == (100, 100)
    assert span("21\t100\t.\tAT\tGC\t.\t.\t.") == (100, 101)
    assert span("21\t100\t.\tAT\tG\t.\t.\t.") == (100, 101)
    assert span("21\t100\t.\tAT\tA,ATT\t.\t.\t.") == (101, 101)
    assert span("21\t100\t.\tAT\tA,*\t.\t.\t.") == (101, 101)
    assert span("21\t100\t.\tAT\tA,G\t.\t.\t.") == (100, 101)
    assert span("21\t100\t.\tA\tG,T\t.\t.\t.") == (100, 100)
    insertion, deletion, two_bases = ("21\t100\t.\tA\tATT\t.\t.\t.", "21\t100\t.\tAT\tA\t.\t.\t.", "21\t100\t.\tAT\tGC\t.\t.\t.")
    assert bgc.row_record_span("21:100-101", [insertion]) == (101, 100)
    assert bgc.row_record_span("21:100-101", [two_bases]) == (100, 101)
    assert bgc.row_record_span("21:101", [insertion, deletion]) == (101, 101)
    assert bgc.row_record_span("21:101", [deletion, deletion]) == (101, 101)
    assert bgc.row_record_span("21:100-101", [insertion, two_bases]) is None
    assert bgc.row_record_span("21:102", [insertion, deletion]) is None
    assert bgc.row_record_span("21:101", []) is None


def test_hgvs_shifted_allele_turns_the_written_allele_as_the_shift_does():
    """`shift_feature_seqs` moves the first base of the allele as written to its end once per
    base of the shift on a forward-strand transcript, and the allele's length less the shift
    times on a reverse-strand one, not at all where that is not positive
    (TranscriptVariationAllele.pm:1330-1345); a deletion has no bases to turn."""
    turned = bgc.hgvs_shifted_allele
    assert [turned("GTGA", offset, False) for offset in (1, 2, 4, 5, -1)] == ["TGAG", "GAGT", "GTGA", "TGAG", "TGAG"]
    assert [turned("GTGA", offset, True) for offset in (-1, -3, -4, -5, 1)] == ["AGTG", "TGAG", "GTGA", "GTGA", "AGTG"]
    assert turned("", 3, False) == turned("", -3, True) == ""


def test_coding_model_places_a_genomic_span_on_the_coding_sequence():
    """`CodingModel.cds_pair` gives a span of bases the CDS coordinates of its two ends, in the
    transcript's order on either strand, and an insertion the bases between its flanks'
    coordinates (Mapper.pm:504-558): a start one past its end between neighbours, the bases the
    model carries where two pairs touch on the genome, and the place against its one mapped
    flank at an exon's edge. A span with an end in an intron, a UTR or outside the transcript
    has no pair, nor has an insertion with a flank in a UTR where the model's bases between
    its flanks open or close the coding sequence (TranscriptMapper.pm:441-443), and the start
    exon's phase moves every coordinate by the bases it pads."""
    def model(pairs, coding_start=4, coding_end=30):
        return bgc.CodingModel(tuple(pairs), coding_start, coding_end, 8)

    forward, reverse, exons = model(ONE_EXON), model(REVERSE_EXON), model(TWO_EXONS)
    assert forward.cdna(1001) == 1 and forward.cdna(1045) == 45 and forward.cdna(1046) is None
    assert reverse.cdna(1045) == 1 and reverse.cdna(1001) == 45
    assert forward.cds_pair(1007, 1007, False, 0) == (4, 4)
    assert forward.cds_pair(1007, 1009, False, 0) == (4, 6)
    assert forward.cds_pair(1007, 1009, False, 2) == (6, 8)
    assert reverse.cds_pair(1037, 1039, True, 0) == (4, 6)
    assert forward.cds_pair(1007, 1006, False, 0) == (4, 3)
    assert reverse.cds_pair(1040, 1039, True, 0) == (4, 3)
    assert forward.cds_pair(1004, 1003, False, 0) is None and forward.cds_pair(1004, 1004, False, 0) == (1, 1)
    assert forward.cds_pair(1030, 1030, False, 0) == (27, 27) and forward.cds_pair(1031, 1031, False, 0) is None
    assert forward.cds_pair(1031, 1030, False, 0) is None and forward.cds_pair(1030, 1029, False, 0) == (27, 26)
    assert forward.cds_pair(1045, 1046, False, 0) is None
    assert exons.cds_pair(1015, 1018, False, 0) == (12, 13) and exons.cds_pair(1015, 1016, False, 0) is None
    assert exons.cds_pair(1016, 1015, False, 0) == (13, 12) == exons.cds_pair(1018, 1017, False, 0)
    assert exons.cds_pair(1017, 1016, False, 0) is None
    assert model(ONE_BASE_SEAM).cds_pair(1016, 1015, False, 0) == (13, 13)
    assert model(TWO_BASE_SEAM).cds_pair(1015, 1014, False, 0) == (12, 13)
    assert model(TWO_BASE_SEAM, 14, 30).cds_pair(1015, 1014, False, 0) == (2, 3)
    assert model(TWO_BASE_SEAM, 15, 30).cds_pair(1015, 1014, False, 0) is None
    assert model(TWO_BASE_SEAM, 4, 17).cds_pair(1015, 1014, False, 0) == (12, 13)
    assert model(TWO_BASE_SEAM, 4, 16).cds_pair(1015, 1014, False, 0) is None


def test_downstream_coding_models_reads_the_mapped_coding_transcripts_of_the_corpus_cache(tmp_path: Path):
    """Each transcript carrying a `cdna_coding_start`, a `cdna_coding_end` and mapper pairs comes,
    under its chromosome and stable id, with them as numbers, the converted cache's strings
    included, and the length of its `peptide`, `None` when the cache holds none; a transcript
    lacking one of the three, or with a pair lacking a key or holding a value that is no number,
    is left out."""
    corpus = tmp_path / "c"
    _manifest(corpus, [])
    shard = corpus / "json_cache" / "transcripts" / "21"
    shard.mkdir(parents=True)
    pair = dict(zip(bgc.MAPPER_PAIR_KEYS, ONE_EXON[0]))
    text = {key: str(value) for key, value in pair.items()}
    span = {"cdna_coding_start": 4, "cdna_coding_end": 30}
    (shard / "1-1000000.json").write_text(json.dumps([
        {"stable_id": "ENST1", **span, "variation_effect_feature_cache": {"mapper": {"pairs": [pair]}, "peptide": "MKTLAVRG"}},
        {"stable_id": "ENST2", "cdna_coding_start": "4", "cdna_coding_end": "30",
         "variation_effect_feature_cache": {"mapper": {"pairs": [text, text]}}},
        {"stable_id": "ENST3", "cdna_coding_start": 4, "variation_effect_feature_cache": {"mapper": {"pairs": [pair]}}},
        {"stable_id": "ENST4", **span, "variation_effect_feature_cache": {"mapper": {"pair_count": 1}}},
        {"stable_id": "ENST5", **span, "variation_effect_feature_cache": {"mapper": {"pairs": [{**pair, "ori": "+"}]}}},
        {"stable_id": "ENST6", **span, "variation_effect_feature_cache": {"mapper": {"pairs": [{"from_start": 1}]}}},
        {"stable_id": "ENST7", "cdna_coding_start": None, "cdna_coding_end": None,
         "variation_effect_feature_cache": {"mapper": {"pairs": [pair]}}},
        {"stable_id": "ENST8", **span},
    ]), encoding="utf-8")
    assert bgc.downstream_coding_models(corpus) == {
        ("21", "ENST1"): bgc.CodingModel((ONE_EXON[0],), 4, 30, 8),
        ("21", "ENST2"): bgc.CodingModel((ONE_EXON[0], ONE_EXON[0]), 4, 30, None),
    }


def test_classify_names_the_class_of_a_downstream_value_read_at_the_hgvs_shifted_span(tmp_path: Path):
    """Under `--hgvs` Ensembl VEP's plugin reads the span the notation's 3' shift leaves, the
    record's own moved by the row's `HGVS_OFFSET`, and vep-rs the record's own: where that
    changes Ensembl VEP's value the row is `plugin_hgvs_shifted_cds_coordinates`. Bases 4 to 7
    are all A. A deletion of base 4 is written at base 7, past the codon's end: both sequences
    are ATG AAA CTC TTG CTG TCC GTG GTT AAG GGC CCA AAT AGC, which vep-rs reads from AAA
    (`KLLLSVVKGPNS`) and Ensembl VEP from CTC. An A inserted before base 4 is written after
    base 7: ATG AAA AAC TCT TGC TGT CCG TGG TTA AGG GCC CAA ATA GCC, read from AAA by vep-rs
    and from AAC by Ensembl VEP. TC inserted after base 8 is written as CT after base 11, the
    allele turned once per base of the shift: ATG AAA ACT CTC TTG CTG ..., read from ACT by
    vep-rs and from CTC by Ensembl VEP. A deletion of base 11, a codon's second base, is
    written at base 12, the same codon's third: both rules read CTG CTG TCC ... on base 11
    (`LLSVVKGPNS`) and Ensembl VEP's starts one codon later on base 12."""
    classes, summary = _mapped_downstream_classes(tmp_path, {
        "deletion": (ONE_EXON, 1006, "GA", "G", "4", *_shifted(1, 3, ("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5))),
        "insertion": (ONE_EXON, 1006, "G", "GA", "3-4", *_shifted(1, 4, ("NSCCPWLRAQIA", 6), ("KNSCCPWLRAQIA", 6))),
        "turned_allele": (ONE_EXON, 1011, "C", "CTC", "8-9", *_shifted(1, 3, ("LLLSVVKGPNS", 6), ("TLLLSVVKGPNS", 6))),
        "inside_one_codon": (ONE_EXON, 1013, "CT", "C", "11", *_shifted(1, 1, ("LSVVKGPNS", 5), ("LLSVVKGPNS", 5))),
    })
    assert classes == {(name, "DownstreamProtein"): HGVS_SHIFTED_CDS_COORDINATES
                       for name in ("deletion", "insertion", "turned_allele", "inside_one_codon")}
    assert summary == {HGVS_SHIFTED_CDS_COORDINATES: 4}


def test_classify_reads_a_shifted_downstream_span_on_a_reverse_strand_transcript(tmp_path: Path):
    """A reverse-strand row prints a negative `HGVS_OFFSET`, the shift along the genome, and its
    allele turns its length less the shift times before it is read on the transcript's strand.
    The deletion of coding base 4, genomic base 1039, is written at coding base 7, three genomic
    bases lower. GTGA, TCAC on the transcript, inserted between coding bases 11 and 12 is
    written one base on as CACT, which is GTGA turned three times, AGTG, on the transcript's
    strand: ATG AAA ACT CTT CAC TGC TGT CCG TGG TTA AGG GCC CAA ATA GCC, read from CTT by vep-rs
    (`LHCCPWLRAQIA`) and from CAC by Ensembl VEP. The same row with the offset's sign turned
    moves the other way along the transcript, gives another sequence and stays
    `unexplained_residual`."""
    insertion = (("HCCPWLRAQIA", 7), ("LHCCPWLRAQIA", 7))
    classes, summary = _mapped_downstream_classes(tmp_path, {
        "deletion": (REVERSE_EXON, 1038, "TT", "T", "4", *_shifted(-1, -3, ("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5))),
        "insertion": (REVERSE_EXON, 1031, "A", "AGTGA", "11-12", *_shifted(-1, -1, *insertion)),
        "other_direction": (REVERSE_EXON, 1031, "A", "AGTGA", "11-12", *_shifted(-1, 1, *insertion)),
    })
    assert classes == {("deletion", "DownstreamProtein"): HGVS_SHIFTED_CDS_COORDINATES,
                       ("insertion", "DownstreamProtein"): HGVS_SHIFTED_CDS_COORDINATES,
                       ("other_direction", "DownstreamProtein"): "unexplained_residual"}
    assert summary == {HGVS_SHIFTED_CDS_COORDINATES: 2, "unexplained_residual": 1}


def test_classify_moves_a_shifted_downstream_span_through_the_transcripts_mapper_pairs(tmp_path: Path):
    """The shifted span is the record's genomic span moved by the offset and placed on the
    coding sequence by the transcript's mapper pairs, so a shift across an intron moves the CDS
    pair by fewer bases than the offset. A T inserted between coding bases 11 and 12, both T,
    is written three genomic bases on, past the last exon base and a two-base intron, before
    coding base 13: ATG AAA ACT CTT TGC TGT CCG TGG TTA AGG GCC CAA ATA GCC, read from CTT by
    vep-rs (`LCCPWLRAQIA`) and from TGC by Ensembl VEP. An offset of 2 leaves both flanks in
    the intron, where the plugin reads the record's own span, and on one exon the offset of 3
    moves the pair three coding bases: no rule gives the pair on either row."""
    values = (("CCPWLRAQIA", 6), ("LCCPWLRAQIA", 6))
    classes, summary = _mapped_downstream_classes(tmp_path, {
        "across_the_intron": (TWO_EXONS, 1014, "T", "TT", "11-12", *_shifted(1, 3, *values)),
        "into_the_intron": (TWO_EXONS, 1014, "T", "TT", "11-12", *_shifted(1, 2, *values)),
        "one_exon": (ONE_EXON, 1014, "T", "TT", "11-12", *_shifted(1, 3, *values)),
    })
    assert classes == {("across_the_intron", "DownstreamProtein"): HGVS_SHIFTED_CDS_COORDINATES,
                       ("into_the_intron", "DownstreamProtein"): "unexplained_residual",
                       ("one_exon", "DownstreamProtein"): "unexplained_residual"}
    assert summary == {HGVS_SHIFTED_CDS_COORDINATES: 1, "unexplained_residual": 2}


def test_classify_reads_a_downstream_insertion_on_a_seam_of_the_model(tmp_path: Path):
    """Where two mapper pairs touch on the genome and the model carries bases between them, an
    insertion at that point replaces those bases (Mapper.pm:504-558, Downstream.pm:108-110). A T
    inserted between coding bases 11 and 12 is written one base on, at a seam that carries
    coding base 13: there Ensembl VEP puts the T in that base's place, ATG AAA ACT CTT TCT GTC
    CGT GGT TAA, a sequence of the reference's length whose protein from TCT is `SVRG`, where
    vep-rs inserts it at the record's own point and reads 14 residues, 6 more than the
    peptide's 8: both fields are `plugin_hgvs_shifted_cds_coordinates`. An A inserted at a seam
    that carries coding bases 12 and 13 prints that span in order, a codon's third base first:
    ATG AAA ACT CTA CTG TCC GTG GTT AAG GGC CCA AAT AGC, read from CTA by vep-rs and from CTG
    by Ensembl VEP, `downstream_start_past_own_codon` on the row's own span. The same column on
    a model without the seam is an insertion's two flanks, which no rule gives that pair for."""
    on_the_seam = (DOWNSTREAM_EXTRA.format(1, "LSVVKGPNS", 5), DOWNSTREAM_EXTRA.format(1, "LLSVVKGPNS", 5))
    classes, summary = _mapped_downstream_classes(tmp_path, {
        "onto_the_seam": (ONE_BASE_SEAM, 1014, "T", "TT", "11-12", *_shifted(1, 1, ("SVRG", 0), ("LCCPWLRAQIA", 6))),
        "on_the_seam": (TWO_BASE_SEAM, 1014, "T", "TA", "12-13", *on_the_seam),
        "no_seam": (ONE_EXON, 1015, "T", "TA", "12-13", *on_the_seam),
    })
    assert classes == {("onto_the_seam", "DownstreamProtein"): HGVS_SHIFTED_CDS_COORDINATES,
                       ("onto_the_seam", "ProteinLengthChange"): HGVS_SHIFTED_CDS_COORDINATES,
                       ("on_the_seam", "DownstreamProtein"): START_PAST_OWN_CODON,
                       ("no_seam", "DownstreamProtein"): "unexplained_residual"}
    assert summary == {HGVS_SHIFTED_CDS_COORDINATES: 2, START_PAST_OWN_CODON: 1, "unexplained_residual": 1}


def test_classify_keeps_the_downstream_class_of_a_span_whose_shifted_reading_is_its_own(tmp_path: Path):
    """A row whose own span opens on a codon's third base and whose shifted reading gives
    Ensembl VEP's own-span value stays `downstream_start_past_own_codon`: the deletion of base
    6 is written at base 7 and Ensembl VEP reads from CTC on either, vep-rs from AAA. Where the
    shift and the start rule both change the value the row is
    `plugin_hgvs_shifted_cds_coordinates`: on the second transcript bases 9 to 12 are all A, and
    the deletion of base 9, a codon's third, is written at base 12; the sequence is ATG TAC GGA
    AAC TGT AAG GGC CCA AAT AGC, which vep-rs reads from GGA (`GNCKGPNS`), Ensembl VEP's rule
    from AAC on base 9 and from TGT on base 12. The value that rule gives on base 9 is not what
    Ensembl VEP writes for the row, so a row carrying it stays `unexplained_residual`."""
    classes, summary = _mapped_downstream_classes(tmp_path, {
        "own_value": (ONE_EXON, 1008, "AA", "A", "6", *_shifted(1, 1, ("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5))),
        "both_rules": (STOP_EXON, 1011, "GA", "G", "9", *_shifted(1, 3, ("CKGPNS", 5), ("GNCKGPNS", 5))),
        "start_rule_alone": (STOP_EXON, 1011, "GA", "G", "9", *_shifted(1, 3, ("NCKGPNS", 5), ("GNCKGPNS", 5))),
    }, coding={"both_rules": DOWNSTREAM_STOP_CDS, "start_rule_alone": DOWNSTREAM_STOP_CDS})
    assert classes == {("own_value", "DownstreamProtein"): START_PAST_OWN_CODON,
                       ("both_rules", "DownstreamProtein"): HGVS_SHIFTED_CDS_COORDINATES,
                       ("start_rule_alone", "DownstreamProtein"): "unexplained_residual"}
    assert summary == {START_PAST_OWN_CODON: 1, HGVS_SHIFTED_CDS_COORDINATES: 1, "unexplained_residual": 1}


def test_classify_records_the_second_class_of_a_downstream_value_two_causes_change(tmp_path: Path):
    """Where the shift and the start rule both change a `DownstreamProtein`, the entry is
    `plugin_hgvs_shifted_cds_coordinates` and names `downstream_start_past_own_codon` under
    `second_expected_divergence`, the class the row's own span has. On the second transcript
    the deletion of base 9, a codon's third, is written at base 12: vep-rs reads from GGA,
    Ensembl VEP's rule from TGT on base 12 and from AAC on base 9, three values. An entry one
    cause changes carries no second class: the deletion of base 4 written at base 7, where
    vep-rs's value is Ensembl VEP's on base 4, and the deletion of base 6 written at base 7,
    where Ensembl VEP's rule gives one value on both. The summary counts each entry once,
    under its own class."""
    names = ["both_causes", "flag_alone", "start_rule_alone"]
    one_codon_late = (("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5))
    _, summary = _mapped_downstream_classes(tmp_path, {
        "both_causes": (STOP_EXON, 1011, "GA", "G", "9", *_shifted(1, 3, ("CKGPNS", 5), ("GNCKGPNS", 5))),
        "flag_alone": (ONE_EXON, 1006, "GA", "G", "4", *_shifted(1, 3, *one_codon_late)),
        "start_rule_alone": (ONE_EXON, 1008, "AA", "A", "6", *_shifted(1, 1, *one_codon_late)),
    }, coding={"both_causes": DOWNSTREAM_STOP_CDS})
    entries = json.loads((tmp_path / "c" / "manifest.json").read_text())["field_divergences"]
    assert {names[d["record_indices"][0]]: (d["expected_divergence"], d.get("second_expected_divergence")) for d in entries} == {
        "both_causes": (HGVS_SHIFTED_CDS_COORDINATES, START_PAST_OWN_CODON),
        "flag_alone": (HGVS_SHIFTED_CDS_COORDINATES, None),
        "start_rule_alone": (START_PAST_OWN_CODON, None),
    }
    assert {names[d["record_indices"][0]] for d in entries if "second_expected_divergence" in d} == {"both_causes"}
    assert summary == {HGVS_SHIFTED_CDS_COORDINATES: 2, START_PAST_OWN_CODON: 1}


def test_classify_leaves_a_downstream_pair_the_hgvs_shift_does_not_give_unexplained(tmp_path: Path):
    """The pair of a deletion of base 4 written at base 7 stays `unexplained_residual` on a row
    without an `HGVS_OFFSET`; with an offset that ends in the 3' UTR, which does not reach the
    plugin; with one that leaves the span in its own codon, where the shifted reading is the
    own one; with either value changed or absent; on a transcript the cache lacks or gives no
    mapper pairs for. A `ProteinLengthChange` has no class with another value on either side,
    nor on a transcript whose cache holds no peptide, beside a `DownstreamProtein` that has."""
    extra = DOWNSTREAM_EXTRA.format
    pair = (("LLLSVVKGPNS", 5), ("KLLLSVVKGPNS", 5))
    deletion = (ONE_EXON, 1006, "GA", "G", "4")
    absent = "IMPACT=HIGH;STRAND=1;HGVS_OFFSET=3"
    seam = (ONE_BASE_SEAM, 1014, "T", "TT", "11-12")
    classes, summary = _mapped_downstream_classes(tmp_path, {
        "no_offset": (*deletion, *(extra(1, *values) for values in pair)),
        "into_the_utr": (*deletion, *_shifted(1, 24, *pair)),
        "own_codon": (*deletion, *_shifted(1, 1, *pair)),
        "other_vep_value": (*deletion, *_shifted(1, 3, ("LLLSVVKGPN", 5), pair[1])),
        "other_vep_rs_value": (*deletion, *_shifted(1, 3, pair[0], ("KLLLSVVKGPN", 5))),
        "no_vep_value": (*deletion, absent, _shifted(1, 3, *pair)[1]),
        "no_vep_rs_value": (*deletion, _shifted(1, 3, *pair)[0], absent),
        "uncached": (None, *deletion[1:], *_shifted(1, 3, *pair)),
        "no_mapper": (*deletion, *_shifted(1, 3, *pair)),
        "other_vep_length": (*seam, *_shifted(1, 1, ("SVRG", 1), ("LCCPWLRAQIA", 6))),
        "other_vep_rs_length": (*seam, *_shifted(1, 1, ("SVRG", 0), ("LCCPWLRAQIA", 7))),
        "no_peptide": (*seam, *_shifted(1, 1, ("SVRG", 0), ("LCCPWLRAQIA", 6))),
    }, without={"no_mapper": ("mapper",), "no_peptide": ("peptide",)})
    unexplained = ["no_offset", "into_the_utr", "own_codon", "other_vep_value", "other_vep_rs_value", "no_vep_value",
                   "no_vep_rs_value", "uncached", "no_mapper"]
    lengths = ["no_vep_value", "no_vep_rs_value", "other_vep_length", "other_vep_rs_length", "no_peptide"]
    assert classes == {
        **{(name, "DownstreamProtein"): "unexplained_residual" for name in unexplained},
        **{(name, "ProteinLengthChange"): "unexplained_residual" for name in lengths},
        **{(name, "DownstreamProtein"): HGVS_SHIFTED_CDS_COORDINATES
           for name in ("other_vep_length", "other_vep_rs_length", "no_peptide")},
    }
    assert summary == {"unexplained_residual": 14, HGVS_SHIFTED_CDS_COORDINATES: 3}


def test_classify_leaves_out_the_masked_regulatory_rows(tmp_path: Path):
    """Under `masks: ["regulatory_feature_rows"]` the reference's RegulatoryFeature and
    MotifFeature rows are not VEP-only keys; without the mask they are unexplained
    residuals, and an unknown class or a `masks` that is not a list of names is refused."""
    regulatory = "rs1\t21:100\tG\t-\tENSR1\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
    motif = "rs1\t21:100\tG\t-\tENSM1\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
    vep = tmp_path / "vep.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant") + regulatory + motif, encoding="utf-8")
    rs = tmp_path / "rs.txt"
    rs.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    for masks, expected in (([], {"unexplained_residual": 2}), (["regulatory_feature_rows"], {})):
        corpus = tmp_path / f"c{len(masks)}"
        corpus.mkdir()
        (corpus / "manifest.json").write_text(json.dumps({"masks": masks, "records": [_record("21\t100\trs1\tA\tG\t.\tPASS\t.")]}), encoding="utf-8")
        assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]) == 0
        m = json.loads((corpus / "manifest.json").read_text())
        assert m["divergence_summary"] == expected and m["vep_rs_only_summary"] == {}
    corpus = tmp_path / "unknown"
    corpus.mkdir()
    (corpus / "manifest.json").write_text(json.dumps({"masks": ["splice_rows"], "records": []}), encoding="utf-8")
    with pytest.raises(SystemExit, match="unknown classes \\['splice_rows'\\]"):
        bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])
    (corpus / "manifest.json").write_text(json.dumps({"masks": "regulatory_feature_rows", "records": []}), encoding="utf-8")
    with pytest.raises(SystemExit, match="masks is not a list of names"):
        bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])


def test_classify_marks_reference_skipped_records_and_their_tuples(tmp_path: Path):
    """A record no VEP row names carries VEP's reason as `reference_warning`; classify writes
    `reference_rows` 0 on it, classes every key only vep-rs emits for it
    `reference_skipped_record`, and removes a `reference_rows` another run left on a record
    that has rows."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\trs1\tA\tG\t.\tPASS\t.", reference_rows=1),
        _record("21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900",
                reference_warning="CPX is not a supported structural variant type"),
    ])
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rs.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant")
        + _row("sv1", "21:501-900", "CPX", "ENST2", "feature_truncation")
        + _row("sv1", "21:501-900", "CPX", "ENST3", "feature_truncation,intron_variant"),
        encoding="utf-8",
    )
    rc = bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])
    assert rc == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert "reference_rows" not in m["records"][0]
    assert list(m["records"][1].items())[-2:] == [
        ("reference_rows", 0), ("reference_warning", "CPX is not a supported structural variant type"),
    ]
    assert [(t["feature"], t["record_indices"], t["expected_divergence"]) for t in m["vep_rs_only_tuples"]] == [
        ("ENST2", [1], "reference_skipped_record"),
        ("ENST3", [1], "reference_skipped_record"),
    ]
    assert m["vep_rs_only_summary"] == {"reference_skipped_record": 2}
    assert m["divergences"] == [] and m["divergence_summary"] == {}


def test_classify_refuses_a_zero_row_record_without_a_warning(tmp_path: Path, capsys):
    """A record no VEP row names and no `reference_warning` explains is an error: its name
    may be unresolved, or the corpus builder did not record VEP's reason."""
    corpus = tmp_path / "c"
    _manifest(corpus, [_record("21\t100\trs1\tA\tG\t.\tPASS\t."), _record("21\t150\trs9\tC\tT\t.\tPASS\t.")])
    before = (corpus / "manifest.json").read_text()
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rs.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rc = bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])
    assert rc == 1
    assert capsys.readouterr().err == "ERROR: [build_golden_corpus] record 1 (21 150 rs9 C T) has no VEP row and no reference_warning\n"
    assert (corpus / "manifest.json").read_text() == before, "a refused run writes nothing"


def test_classify_keeps_an_unwarned_non_variant_record(tmp_path: Path, capsys):
    """A record whose ALT is `.` carries no row by its nature (under `--allow_non_variant` VEP
    writes its VCF line and JSON object and no default row), so one without a warning is
    neither refused nor marked skipped; one VEP warned about is skipped as any other."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\trs1\tA\tG\t.\tPASS\t."),
        _record("21\t200\tnonvar\tTG\t.\t.\tPASS\t."),
        _record("21\t300\tnonvar_bad\tG\t.\t.\tPASS\t.",
                reference_warning="Specified reference allele G does not match Ensembl reference allele"),
    ])
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rs.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rc = bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])
    assert rc == 0
    assert "reference-skipped records [2]" in capsys.readouterr().out
    records = json.loads((corpus / "manifest.json").read_text())["records"]
    assert "reference_rows" not in records[1]
    assert records[2]["reference_rows"] == 0


def test_classify_keeps_a_record_the_reference_carried_without_a_row(tmp_path: Path, capsys):
    """A record whose JSON object the reference wrote with `most_severe_consequence` `?` has no
    allele to annotate (`C C`, `C <DEL:*>`, `A *`), so one without a warning is kept, not refused
    or marked skipped; the object is matched by the record's first five columns, so a twin
    input carrying sample columns names it too; a record the JSON lacks is refused, and a run
    whose expected directory has no JSON file falls back to the ALT `.` rule."""
    corpus = tmp_path / "c"
    _manifest(corpus, [
        _record("21\t100\trs1\tA\tG\t.\tPASS\t."),
        _record("21\t200\talt_eq_ref\tC\tC\t.\tPASS\t."),
        _record("21\t300\tlone_delstar\tC\t<DEL:*>\t.\tPASS\t."),
        _record("21\t400\tstar_only\tA\t*\t.\tPASS\t."),
    ])

    def write_expected_json(objects: list[dict]) -> None:
        (corpus / "expected").mkdir(exist_ok=True)
        with gzip.open(corpus / "expected" / "json.jsonl.gz", "wt", encoding="utf-8") as fh:
            fh.write("".join(json.dumps(o) + "\n" for o in objects))

    objects = [
        {"id": "rs1", "input": "21\t100\trs1\tA\tG\t.\tPASS\t.", "most_severe_consequence": "missense_variant"},
        {"id": "alt_eq_ref", "input": "21\t200\talt_eq_ref\tC\tC\t.\tPASS\t.\tGT\t0/1", "most_severe_consequence": "?"},
        {"id": "lone_delstar", "input": "21\t300\tlone_delstar\tC\t<DEL:*>\t.\tPASS\t.", "most_severe_consequence": "?"},
    ]
    write_expected_json(objects)
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rs.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    classify = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main(classify) == 1
    assert capsys.readouterr().err == (
        "ERROR: [build_golden_corpus] record 3 (21 400 star_only A *) has no VEP row and no reference_warning\n"
    )
    objects.append({"id": "star_only", "input": "21\t400\tstar_only\tA\t*\t.\tPASS\t.", "most_severe_consequence": "?"})
    write_expected_json(objects)
    assert bgc.main(classify) == 0
    assert "reference-skipped records []" in capsys.readouterr().out
    records = json.loads((corpus / "manifest.json").read_text())["records"]
    assert all("reference_rows" not in r and "reference_warning" not in r for r in records)
    (corpus / "expected" / "json.jsonl.gz").unlink()
    assert bgc.main(classify) == 1


def test_classify_reads_a_flag_sets_own_json_and_refuses_an_unreadable_one(tmp_path: Path, capsys):
    """A flag set's carried records come from its own expected JSON when it compares that
    format, else from the default run's; a JSON file that is not gzip, or holds a line that is
    not JSON or not a JSON object, is refused with the tool's own error instead of a traceback."""
    corpus = tmp_path / "c"
    corpus.mkdir()
    (corpus / "manifest.json").write_text(json.dumps({
        "records": [_record("21\t100\trs1\tA\tG\t.\tPASS\t."), _record("21\t200\talt_eq_ref\tC\tC\t.\tPASS\t.")],
        "flag_sets": [
            {"name": "own", "flags": ["--pick"], "expected_dir": "expected_own", "covers": []},
            {"name": "shared", "flags": ["--pick"], "expected_dir": "expected", "covers": []},
        ],
    }), encoding="utf-8")
    (corpus / "expected").mkdir()
    (corpus / "expected_own").mkdir()
    carried = {"id": "alt_eq_ref", "input": "21\t200\talt_eq_ref\tC\tC\t.\tPASS\t.", "most_severe_consequence": "?"}
    with gzip.open(corpus / "expected_own" / "json.jsonl.gz", "wt", encoding="utf-8") as fh:
        fh.write(json.dumps(carried) + "\n")
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rs.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    classify = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main([*classify, "--flag-set", "own"]) == 0
    assert bgc.main([*classify, "--flag-set", "shared"]) == 1
    assert "record 1 (21 200 alt_eq_ref C C) has no VEP row" in capsys.readouterr().err
    (corpus / "expected" / "json.jsonl.gz").write_bytes(b"not gzip")
    with pytest.raises(SystemExit, match=r"ERROR: \[build_golden_corpus\] .*json\.jsonl\.gz is not a readable JSON output"):
        bgc.main([*classify, "--flag-set", "shared"])
    with gzip.open(corpus / "expected" / "json.jsonl.gz", "wt", encoding="utf-8") as fh:
        fh.write("{not json\n")
    with pytest.raises(SystemExit, match=r"is not a readable JSON output"):
        bgc.main(classify)
    with gzip.open(corpus / "expected" / "json.jsonl.gz", "wt", encoding="utf-8") as fh:
        fh.write(json.dumps(carried) + "\n[]\n")
    with pytest.raises(SystemExit, match=r"is not a readable JSON output \(line 2 is not a JSON object\)"):
        bgc.main(classify)


def test_classify_keeps_a_record_a_row_filter_left_rowless_on_both_sides(tmp_path: Path, capsys):
    """Under `--coding_only` or `--no_intergenic` a record with no coding (or no intergenic)
    row has no row on either side, so a flag set or a default run under one keeps such a
    record without a warning, not refused and not marked skipped, with no JSON to carry it; the
    same record under no row filter, or with vep-rs rows and no reference row, is refused."""
    corpus = tmp_path / "c"
    corpus.mkdir()
    records = [_record("21\t100\trs1\tA\tG\t.\tPASS\t."), _record("21\t150\trs9\tC\tT\t.\tPASS\t.")]
    (corpus / "manifest.json").write_text(json.dumps({
        "records": records,
        "flag_sets": [
            {"name": "coding-only", "flags": ["--coding_only"], "expected_dir": "expected_coding-only", "covers": []},
            {"name": "no-intergenic", "flags": ["--no_intergenic"], "expected_dir": "expected_no-intergenic", "covers": []},
            {"name": "pick", "flags": ["--pick"], "expected_dir": "expected_pick", "covers": []},
        ],
    }), encoding="utf-8")
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rs.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    classify = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    for flag_set in ("coding-only", "no-intergenic"):
        assert bgc.main([*classify, "--flag-set", flag_set]) == 0
        assert "reference-skipped records []" in capsys.readouterr().out
    assert bgc.main([*classify, "--flag=--coding_only"]) == 0
    capsys.readouterr()
    m = json.loads((corpus / "manifest.json").read_text())
    assert all("reference_rows" not in r and "reference_warning" not in r for r in m["records"])
    assert m["flags"] == ["--coding_only"] and m["divergences"] == [] and m["vep_rs_only_tuples"] == []
    assert all(s["divergences"] == [] for s in m["flag_sets"] if s["name"] != "pick")
    refused = "ERROR: [build_golden_corpus] record 1 (21 150 rs9 C T) has no VEP row and no reference_warning\n"
    assert bgc.main([*classify, "--flag-set", "pick"]) == 1
    assert capsys.readouterr().err == refused
    rs.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant") + _row("rs9", "21:150", "T", "ENST2", "missense_variant"),
        encoding="utf-8",
    )
    assert bgc.main([*classify, "--flag-set", "coding-only"]) == 1
    assert capsys.readouterr().err == refused


def test_classify_numbers_a_flag_sets_own_input_records(tmp_path: Path, capsys):
    """A flag set over another `input` runs that file's records, numbered as the golden tests
    number them (every non-header line, in order), so its divergences and skipped records
    resolve by those ordinals; a line the manifest records, a twin of `variants.vcf` in
    another encoding (CRLF line ends, gzip), keeps its entry, so its `reference_warning` holds
    and `reference_rows` lands in the manifest; a line that is no VCF record (a space-delimited
    or an Ensembl-format re-encoding of `variants.vcf`) is the manifest record at its ordinal,
    and refused when the file's record count is not the manifest's; a line of a VCF input that
    is neither a twin nor a VCF record is refused with its ordinal; a record the reference
    skipped in an input the manifest does not record is refused; an unreadable input is
    refused with its path."""
    corpus = tmp_path / "c"
    corpus.mkdir()
    records = [
        _record("21\t100\trs1\tA\tG\t.\tPASS\t."),
        _record("21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900", reference_warning="CPX is not a supported structural variant type"),
    ]
    (corpus / "manifest.json").write_text(json.dumps({
        "records": records,
        "flag_sets": [
            {"name": "structural", "flags": ["--coding_only"], "input": "structural.vcf", "expected_dir": "expected_structural", "covers": []},
            {"name": "crlf", "flags": [], "input": "variants.crlf.vcf", "expected_dir": "expected_crlf", "covers": []},
            {"name": "missing", "flags": [], "input": "absent.vcf", "expected_dir": "expected_missing", "covers": []},
            {"name": "space", "flags": ["--delimiter", " "], "input": "variants.space.vcf", "expected_dir": "expected_space", "covers": []},
            {"name": "ensembl", "flags": ["--format", "ensembl"], "input": "variants.ensembl.txt", "expected_dir": "expected_ensembl", "covers": []},
            {"name": "gz", "flags": [], "input": "variants.vcf.gz", "expected_dir": "expected_gz", "covers": []},
            {"name": "short", "flags": ["--format", "ensembl"], "input": "one.ensembl.txt", "expected_dir": "expected_short", "covers": []},
        ],
    }), encoding="utf-8")
    (corpus / "structural.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t1000\tsv_a\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=1900\n"
        "21\t3000\tsv_b\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=3900\n"
        "21\t5000\tsv_c\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=5900\n",
        encoding="utf-8",
    )
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(_row("sv_b", "21:3001-3900", "deletion", "ENST2", "coding_sequence_variant,feature_truncation"), encoding="utf-8")
    rs.write_text(
        _row("sv_b", "21:3001-3900", "deletion", "ENST2", "coding_sequence_variant,feature_truncation")
        + _row("sv_b", "21:3001-3900", "deletion", "ENST3", "feature_truncation"),
        encoding="utf-8",
    )
    classify = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main([*classify, "--flag-set", "structural"]) == 0
    assert "reference-skipped records []" in capsys.readouterr().out
    m = json.loads((corpus / "manifest.json").read_text())
    structural = m["flag_sets"][0]
    assert structural["divergences"] == []
    assert [(t["feature"], t["record_indices"], t["expected_divergence"]) for t in structural["vep_rs_only_tuples"]] == [
        ("ENST3", [1], "unexplained_residual"),
    ]
    assert all("reference_rows" not in r for r in m["records"])
    # The twin input carries the manifest's records: the skipped one keeps its warning.
    (corpus / "variants.crlf.vcf").write_bytes(
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\n"
        b"21\t100\trs1\tA\tG\t.\tPASS\t.\r\n21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900\r\n"
    )
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rs.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant") + _row("sv1", "21:501-900", "CPX", "ENST2", "feature_truncation"),
        encoding="utf-8",
    )
    assert bgc.main([*classify, "--flag-set", "crlf"]) == 0
    assert "reference-skipped records [1]" in capsys.readouterr().out
    m = json.loads((corpus / "manifest.json").read_text())
    assert m["records"][1]["reference_rows"] == 0
    assert m["flag_sets"][1]["vep_rs_only_summary"] == {"reference_skipped_record": 1}
    # A space-delimited or Ensembl-format re-encoding of variants.vcf holds no VCF record and
    # takes the manifest's records in order; the gzipped copy is a twin read through gzip.
    (corpus / "variants.space.vcf").write_text(
        "#CHROM POS ID REF ALT QUAL FILTER INFO\n21 100 rs1 A G . PASS .\n21 500 sv1 N <CPX> . PASS SVTYPE=CPX;END=900\n",
        encoding="utf-8",
    )
    (corpus / "variants.ensembl.txt").write_text("21\t100\t100\tA/G\t1\trs1\n21\t500\t900\tCPX\t1\tsv1\n", encoding="utf-8")
    with gzip.open(corpus / "variants.vcf.gz", "wt", encoding="utf-8") as fh:
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t500\tsv1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;END=900\n")
    for name in ("space", "ensembl", "gz"):
        assert bgc.main([*classify, "--flag-set", name]) == 0
        assert "reference-skipped records [1]" in capsys.readouterr().out
    m = json.loads((corpus / "manifest.json").read_text())
    assert all(s["vep_rs_only_summary"] == {"reference_skipped_record": 1} for s in m["flag_sets"] if s["name"] in ("space", "ensembl", "gz"))
    (corpus / "one.ensembl.txt").write_text("21\t100\t100\tA/G\t1\trs1\n", encoding="utf-8")
    with pytest.raises(SystemExit, match=r"one\.ensembl\.txt holds no VCF record and its 1 records are not the manifest's 2"):
        bgc.main([*classify, "--flag-set", "short"])
    # A standalone VCF of short records (five or seven columns) numbers its own records.
    (corpus / "own.short.vcf").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n21\t7000\tsnv_a\tC\tT\n21\t8000\tsnv_b\tG\tA\t.\tPASS\n",
        encoding="utf-8",
    )
    m = json.loads((corpus / "manifest.json").read_text())
    m["flag_sets"].append({"name": "own-short", "flags": [], "input": "own.short.vcf", "expected_dir": "expected_own-short", "covers": []})
    (corpus / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    vep.write_text(_row("snv_a", "21:7000", "T", "ENST7", "missense_variant") + _row("snv_b", "21:8000", "A", "ENST8", "synonymous_variant"), encoding="utf-8")
    rs.write_text(_row("snv_a", "21:7000", "T", "ENST7", "missense_variant") + _row("snv_b", "21:8000", "A", "ENST8", "missense_variant"), encoding="utf-8")
    assert bgc.main([*classify, "--flag-set", "own-short"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    own_short = next(s for s in m["flag_sets"] if s["name"] == "own-short")
    assert [(d["feature"], d["record_indices"], d["expected_divergence"]) for d in own_short["divergences"]] == [
        ("ENST8", [1], "unexplained_residual"),
    ]
    assert "reference_rows" not in m["records"][0]
    # A line of a VCF input that is neither a twin nor a VCF record is refused with its ordinal.
    (corpus / "own.short.vcf").write_text("21\t7000\tsnv_a\tC\tT\n21\t.\tsnv_x\tG\tA\n", encoding="utf-8")
    with pytest.raises(SystemExit, match=r"own\.short\.vcf record 1 \(.*\) is no VCF record in a VCF input"):
        bgc.main([*classify, "--flag-set", "own-short"])
    # A skipped record of an input the manifest does not record has no warning to carry.
    vep.write_text(_row("sv_b", "21:3001-3900", "deletion", "ENST2", "coding_sequence_variant,feature_truncation"), encoding="utf-8")
    rs.write_text(
        _row("sv_b", "21:3001-3900", "deletion", "ENST2", "coding_sequence_variant,feature_truncation")
        + _row("sv_c", "21:5001-5900", "deletion", "ENST4", "feature_truncation"),
        encoding="utf-8",
    )
    assert bgc.main([*classify, "--flag-set", "structural"]) == 1
    assert capsys.readouterr().err == "ERROR: [build_golden_corpus] record 2 (21 5000 sv_c N <DEL>) has no VEP row and no reference_warning\n"
    with pytest.raises(SystemExit, match=r"absent\.vcf is not a readable text file"):
        bgc.main([*classify, "--flag-set", "missing"])


def test_classify_refuses_a_warned_record_that_has_rows(tmp_path: Path, capsys):
    corpus = tmp_path / "c"
    _manifest(corpus, [_record("21\t100\trs1\tA\tG\t.\tPASS\t.", reference_warning="deletion looks incomplete")])
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rs.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    rc = bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])
    assert rc == 1
    assert capsys.readouterr().err == (
        "ERROR: [build_golden_corpus] record 0 (21 100 rs1 A G) has VEP rows and reference_warning 'deletion looks incomplete'\n"
    )


def test_classify_names_the_mate_side_local_read(tmp_path: Path):
    """A point breakend's mate-side Transcript row: VEP reads it at the local coordinate
    (`feature_truncation` alone when the local coordinate lands nowhere near the mate
    transcript), vep-rs at the mate (`feature_truncation,intron_variant` for a mate in
    intron 1 of a transcript filed in the corpus cache), under either cache layout the
    manifest's `cache_layout` names."""
    shard = json.dumps([{
        "stable_id": "ENST1", "start": 33031935, "end": 33041244, "strand": 1, "biotype": "protein_coding",
        "exons": [{"start": 33031935, "end": 33032154}, {"start": 33036103, "end": 33036199}],
        "variation_effect_feature_cache": {"mapper": {"pairs": [
            {"from_start": 1, "from_end": 220, "to_start": 33031935, "to_end": 33032154, "ori": 1},
            {"from_start": 221, "from_end": 317, "to_start": 33036103, "to_end": 33036199, "ori": 1},
        ], "pair_count": 2, "cdna_coding_start": 50, "cdna_coding_end": 300}},
    }])
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(
        _row("b1", "1:14374334", "]21:33034355]A", "ENST1", "feature_truncation")
        + _row("b2", "1:14374334", "]21:33034355]A", "ENST2", "feature_truncation"),
        encoding="utf-8",
    )
    rs.write_text(
        _row("b1", "1:14374334", "]21:33034355]A", "ENST1", "feature_truncation,intron_variant")
        + _row("b2", "1:14374334", "]21:33034355]A", "ENST2", "feature_truncation,intron_variant"),
        encoding="utf-8",
    )
    layouts = [
        (tmp_path / "c", Path("json_cache"), {"records": []}),
        (tmp_path / "e", Path("cache_root", "homo_sapiens", "116_GRCh38"),
         {"release": "116", "assembly": "GRCh38", "cache_layout": "ensembl", "records": []}),
    ]
    for corpus, cache, manifest in layouts:
        (corpus / cache / "transcripts" / "21").mkdir(parents=True)
        (corpus / cache / "transcripts" / "21" / "33000001-34000000.json").write_text(shard, encoding="utf-8")
        (corpus / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        rc = bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])
        assert rc == 0
        m = json.loads((corpus / "manifest.json").read_text())
        by_feature = {d["feature"]: d["expected_divergence"] for d in m["divergences"]}
        assert by_feature["ENST1"] == "breakend_mate_local_read", cache
        # a transcript the corpus cache does not carry is not adjudicated
        assert by_feature["ENST2"] == "unexplained_residual", cache


def _positioned_row(uv, loc, allele, feat, cons, cdna, cds, ftype="Transcript"):
    columns = [uv, loc, allele, "ENSG1", feat, ftype, cons, cdna, cds, "-", "-", "-", "-", "IMPACT=LOW;STRAND=1"]
    return "\t".join(columns) + "\n"


def test_classify_names_the_stop_codon_test_classes(tmp_path: Path):
    """`stop_retained_without_stop_codon_test`: the reference's `stop_retained_variant` against
    nothing, `coding_sequence_variant` or `frameshift_variant` in vep-rs's set, on an insertion's
    or deletion's Transcript row whose cDNA position is missing whole or at one end, which has no
    CDS start, or whose allele holds a letter that names no base. `stop_codon_overlap_genomic_span`:
    the reference's `coding_sequence_variant` against vep-rs's `stop_retained_variant` beside
    `splice_region_variant`, on such a row with its positions.
    `stop_retained_three_base_stop_codon_test`: on a row whose test is made, the reference's
    `stop_retained_variant` against `coding_sequence_variant`, or it, alone or beside one in-frame
    or protein-altering term, against `frameshift_variant` with or without `stop_gained`.
    `stop_codon_removed_without_stop_codon_test`: the reference's `stop_retained_variant`
    against vep-rs's `stop_lost` on a row without a cDNA position at one end or a CDS start.
    Another term in the place, a substitution's row, a row without the splice term and a
    regulatory row stay an `unexplained_residual`."""
    acceptor = "splice_acceptor_variant,"
    region = "splice_region_variant,"
    cases = [
        # (location, allele, cDNA, CDS, reference terms, vep-rs terms, class)
        ("21:101-105", "-", "?-7057", "-", acceptor + "stop_retained_variant", acceptor + "coding_sequence_variant", bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST),
        ("21:201-205", "-", "-", "-", "stop_retained_variant,3_prime_UTR_variant", "coding_sequence_variant,3_prime_UTR_variant", bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST),
        ("21:301-305", "-", "1705-?", "1630-?", acceptor + "stop_retained_variant", acceptor + "frameshift_variant", bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST),
        ("21:401-9405", "-", "120-7100", "-", "start_lost,stop_retained_variant", "start_lost", bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST),
        ("21:501-502", "ACNT", "7038-7039", "6961-6962", "stop_retained_variant", "frameshift_variant", bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST),
        ("21:601", "-", "1705", "1630", region + "coding_sequence_variant", region + "stop_retained_variant", bgc.STOP_CODON_OVERLAP_GENOMIC_SPAN),
        ("21:701-705", "-", "7055-7059", "6978-6982", acceptor + "stop_retained_variant", acceptor + "coding_sequence_variant", bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST),
        ("21:801-805", "-", "7055-7059", "6978-?", "stop_retained_variant,3_prime_UTR_variant", "coding_sequence_variant,3_prime_UTR_variant", bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST),
        ("21:1601-1607", "-", "961-967", "264-270", "inframe_deletion,stop_retained_variant", "frameshift_variant", bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST),
        ("21:1701-1702", "ACGTACGT", "1213-1214", "1100-1101", "inframe_insertion,stop_retained_variant", "stop_gained,frameshift_variant", bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST),
        ("21:1801-1802", "ACGTACGT", "1213-1214", "1100-1101", "protein_altering_variant,stop_retained_variant", "stop_gained,frameshift_variant", bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST),
        ("21:1901-1907", "-", "961-967", "264-270", "inframe_deletion,stop_retained_variant", "coding_sequence_variant", bgc.UNEXPLAINED),
        ("21:2001-2007", "-", "961-967", "264-270", "missense_variant,stop_retained_variant", "frameshift_variant", bgc.UNEXPLAINED),
        ("21:2101-2107", "-", "961-?", "264-?", "inframe_deletion,stop_retained_variant", "frameshift_variant", bgc.UNEXPLAINED),
        ("21:2201-2205", "-", "?-7057", "-", acceptor + "stop_retained_variant", acceptor + "stop_lost", bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST),
        ("21:2301-9405", "-", "120-7100", "-", "start_lost,stop_retained_variant", "stop_lost,start_lost", bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST),
        ("21:2401-2413", "A", "241-?", "242-?", "stop_retained_variant,3_prime_UTR_variant", "stop_lost,3_prime_UTR_variant", bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST),
        ("21:2501-2505", "-", "7055-7059", "6978-?", "stop_retained_variant,3_prime_UTR_variant", "stop_lost,3_prime_UTR_variant", bgc.UNEXPLAINED),
        ("21:2601-2605", "-", "?-7057", "-", acceptor + "stop_retained_variant", acceptor + "stop_lost,frameshift_variant", bgc.UNEXPLAINED),
        ("21:901-905", "-", "?-7057", "-", acceptor + "stop_retained_variant", acceptor + "inframe_deletion", bgc.UNEXPLAINED),
        ("21:1001", "A", "-", "-", "stop_retained_variant", "coding_sequence_variant", bgc.UNEXPLAINED),
        ("21:1101", "-", "1705", "1630", "coding_sequence_variant", "stop_retained_variant", bgc.UNEXPLAINED),
        ("21:1201", "-", "1705-?", "1630-?", region + "coding_sequence_variant", region + "stop_retained_variant", bgc.UNEXPLAINED),
        ("21:1301", "-", "1705", "1630", region + "stop_retained_variant", region + "coding_sequence_variant", bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST),
        # Under --total_length a missing position prints its `-` before the `/total` suffix.
        ("21:1501-1505", "-", "-/7140", "-/7023", acceptor + "stop_retained_variant", acceptor + "coding_sequence_variant", bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST),
    ]
    corpus = tmp_path / "c"
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\trs{i}\tAC\tA\t.\tPASS\t.") for i in range(len(cases) + 1)])
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    regulatory = ("21:1401-1405", "-", "-", "-", "stop_retained_variant", "coding_sequence_variant")
    vep.write_text(
        "".join(_positioned_row(f"rs{i}", loc, allele, f"ENST{i}", perl, cdna, cds) for i, (loc, allele, cdna, cds, perl, _, _) in enumerate(cases))
        + _positioned_row(f"rs{len(cases)}", *regulatory[:2], "ENSR1", regulatory[4], *regulatory[2:4], ftype="RegulatoryFeature"),
        encoding="utf-8",
    )
    rs.write_text(
        "".join(_positioned_row(f"rs{i}", loc, allele, f"ENST{i}", rust, cdna, cds) for i, (loc, allele, cdna, cds, _, rust, _) in enumerate(cases))
        + _positioned_row(f"rs{len(cases)}", *regulatory[:2], "ENSR1", regulatory[5], *regulatory[2:4], ftype="RegulatoryFeature"),
        encoding="utf-8",
    )
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert {d["location"]: d["expected_divergence"] for d in m["divergences"]} == {
        **{loc: klass for loc, _, _, _, _, _, klass in cases}, regulatory[0]: bgc.UNEXPLAINED,
    }
    assert m["divergence_summary"] == {
        bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST: 6, bgc.STOP_CODON_OVERLAP_GENOMIC_SPAN: 1,
        bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST: 6, bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST: 3,
        bgc.UNEXPLAINED: 10,
    }

    # A flag set's record-level row (`--summary`): the one row of a record names no feature and
    # carries every feature's terms. It takes the first class or the removed-codon class where
    # the manifest's own run documents a Transcript row of the record under one of the two and
    # its sets differ in that class's way, and stays unexplained for any other record.
    m["flag_sets"] = [{"name": "summary", "flags": ["--summary"]}]
    (corpus / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    summary_vep = tmp_path / "summary_vep.txt"
    summary_rs = tmp_path / "summary_rs.txt"
    union = "splice_acceptor_variant,{}intron_variant"
    removed = next(i for i, case in enumerate(cases) if case[0] == "21:2201-2205")
    others = "".join(
        _positioned_row(f"rs{i}", loc, allele, "-", union.format(""), "-", "-", ftype="-")
        for i, (loc, allele, *_) in enumerate([*cases, regulatory]) if i not in (0, 6, removed)
    )
    summary_vep.write_text(
        _positioned_row("rs0", "21:101-105", "-", "-", union.format("stop_retained_variant,coding_sequence_variant,"), "-", "-", ftype="-")
        + _positioned_row("rs6", "21:701-705", "-", "-", union.format("stop_retained_variant,"), "-", "-", ftype="-")
        + _positioned_row(f"rs{removed}", "21:2201-2205", "-", "-", union.format("stop_retained_variant,"), "-", "-", ftype="-") + others,
        encoding="utf-8",
    )
    summary_rs.write_text(
        _positioned_row("rs0", "21:101-105", "-", "-", union.format("coding_sequence_variant,"), "-", "-", ftype="-")
        + _positioned_row("rs6", "21:701-705", "-", "-", union.format("coding_sequence_variant,"), "-", "-", ftype="-")
        + _positioned_row(f"rs{removed}", "21:2201-2205", "-", "-", union.format("stop_lost,"), "-", "-", ftype="-") + others,
        encoding="utf-8",
    )
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(summary_vep), "--vep-rs-output", str(summary_rs),
                     "--flag-set", "summary"]) == 0
    fs = json.loads((corpus / "manifest.json").read_text())["flag_sets"][0]
    assert {d["location"]: d["expected_divergence"] for d in fs["divergences"]} == {
        "21:101-105": bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST, "21:701-705": bgc.UNEXPLAINED,
        "21:2201-2205": bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST,
    }


def _classified(tmp_path: Path, cases: list[tuple[str, str, str, str, str, str, str]]) -> dict:
    """The manifest `classify` writes into `tmp_path / "c"` for rows of (location, allele, cDNA, CDS,
    reference terms, vep-rs terms, class)."""
    corpus = tmp_path / "c"
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\trs{i}\tAC\tA\t.\tPASS\t.") for i in range(len(cases))])
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text("".join(_positioned_row(f"rs{i}", loc, allele, f"ENST{i}", perl, cdna, cds) for i, (loc, allele, cdna, cds, perl, _, _) in enumerate(cases)), encoding="utf-8")
    rs.write_text("".join(_positioned_row(f"rs{i}", loc, allele, f"ENST{i}", rust, cdna, cds) for i, (loc, allele, cdna, cds, _, rust, _) in enumerate(cases)), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert {d["location"]: d["expected_divergence"] for d in m["divergences"]} == {loc: klass for loc, _, _, _, _, _, klass in cases}
    return m


def test_classify_names_the_insertion_range_class(tmp_path: Path):
    """`stop_codon_overlap_insertion_range`: on an insertion's Transcript row with both cDNA
    positions and a CDS start, the reference's `coding_sequence_variant` away from
    `splice_region_variant`, or its `frameshift_variant` with or without `stop_gained`, against
    vep-rs's `stop_retained_variant`, the second with or without `inframe_insertion`; an allele
    holding `N` is one. Beside `splice_region_variant` the first pair is the split codon's. A
    row without the CDS start or a cDNA end, a deletion's row, another term in the place and an
    allele of another letter stay an `unexplained_residual`."""
    region = "inframe_insertion,splice_region_variant,"
    kept = "inframe_insertion,stop_retained_variant"
    cases = [
        ("21:101-102", "TTGA", "7054-7055", "6977-6978", "inframe_insertion,coding_sequence_variant", kept, bgc.STOP_CODON_OVERLAP_INSERTION_RANGE),
        ("21:201-202", "A", "7054-7055", "6977-6978", "stop_gained,frameshift_variant", kept, bgc.STOP_CODON_OVERLAP_INSERTION_RANGE),
        ("21:301-302", "TNGA", "7054-7055/7140", "6977-6978/7023", "frameshift_variant", "stop_retained_variant", bgc.STOP_CODON_OVERLAP_INSERTION_RANGE),
        ("21:401-402", "TTGAT", "7054-7055", "6978-?", "frameshift_variant", kept, bgc.STOP_CODON_OVERLAP_INSERTION_RANGE),
        ("21:501-502", "TTGA", "7054-7055", "6977-6978", region + "coding_sequence_variant", region + "stop_retained_variant", bgc.STOP_CODON_OVERLAP_GENOMIC_SPAN),
        ("21:601-602", "TTGA", "7054-7055", "?-6978", "inframe_insertion,coding_sequence_variant", kept, bgc.UNEXPLAINED),
        ("21:701-702", "TTGA", "7054-?", "6977-?", "frameshift_variant", kept, bgc.UNEXPLAINED),
        ("21:801-804", "-", "7055-7058", "6978-6981", "coding_sequence_variant", "stop_retained_variant", bgc.UNEXPLAINED),
        ("21:901-902", "TTGA", "7054-7055", "6977-6978", "frameshift_variant", "protein_altering_variant,stop_retained_variant", bgc.UNEXPLAINED),
        ("21:1001-1002", "TTGA", "7054-7055", "6977-6978", "coding_sequence_variant", kept, bgc.UNEXPLAINED),
        ("21:1101-1102", "TRGA", "7054-7055", "6977-6978", "frameshift_variant", kept, bgc.UNEXPLAINED),
    ]
    m = _classified(tmp_path, cases)
    assert m["divergence_summary"] == {bgc.STOP_CODON_OVERLAP_INSERTION_RANGE: 4, bgc.STOP_CODON_OVERLAP_GENOMIC_SPAN: 1, bgc.UNEXPLAINED: 6}

    # A flag set's record-level row (`--summary`) names no feature and carries every feature's
    # terms. It takes the class where the manifest's own run documents a Transcript row of the
    # record under it and its sets differ in the class's way, and stays unexplained otherwise.
    corpus = tmp_path / "c"
    m["flag_sets"] = [{"name": "summary", "flags": ["--summary"]}]
    (corpus / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    union = "{}downstream_gene_variant"
    differing = {
        0: ("inframe_insertion,coding_sequence_variant,", kept + ","),
        1: ("stop_gained,frameshift_variant,", kept + ","),
        5: ("inframe_insertion,coding_sequence_variant,", kept + ","),
        8: ("frameshift_variant,", "protein_altering_variant,stop_retained_variant,"),
    }
    summary_vep = tmp_path / "summary_vep.txt"
    summary_rs = tmp_path / "summary_rs.txt"
    for path, side in ((summary_vep, 0), (summary_rs, 1)):
        path.write_text("".join(
            _positioned_row(f"rs{i}", loc, "-", "-", union.format(differing.get(i, ("", ""))[side]), "-", "-", ftype="-")
            for i, (loc, *_) in enumerate(cases)
        ), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(summary_vep), "--vep-rs-output", str(summary_rs),
                     "--flag-set", "summary"]) == 0
    fs = json.loads((corpus / "manifest.json").read_text())["flag_sets"][0]
    assert {d["location"]: d["expected_divergence"] for d in fs["divergences"]} == {
        "21:101-102": bgc.STOP_CODON_OVERLAP_INSERTION_RANGE, "21:201-202": bgc.STOP_CODON_OVERLAP_INSERTION_RANGE,
        "21:601-602": bgc.UNEXPLAINED, "21:901-902": bgc.UNEXPLAINED,
    }


def test_classify_names_the_splice_term_past_the_record_span_class(tmp_path: Path):
    """`splice_term_past_the_record_span` on rows no record of the run writes, which are read
    without their record's alleles: on a Transcript row whose Location holds three positions
    under a longer allele, with two cDNA positions as far apart as the Location's ends, the
    reference's donor, acceptor, fifth base, donor region, tract or intron term, one or more,
    against a set of vep-rs's without them, with or without `splice_region_variant` in their
    place; positions carrying the transcript's lengths are read the same.
    `splice_region_variant` where the tract term or the intron term alone is gone,
    `splice_donor_region_variant` in the place of the fifth base term, two neighbouring
    positions, a span across an exon's end, an intron span, an allele as long as the span, a
    deletion, a kept term of the six, a further term and the terms the other way round stay an
    `unexplained_residual`."""
    klass = "splice_term_past_the_record_span"
    donor, back, shift = "splice_donor_variant,frameshift_variant", "frameshift_variant,splice_region_variant", "frameshift_variant"
    tract = "frameshift_variant,splice_polypyrimidine_tract_variant"
    cases = [
        ("21:101-103", "GTTCA", "998-1000", "898-900", donor, back, klass),
        ("21:201-203", "GTTCA", "998-1000", "898-900", donor, shift, klass),
        ("21:301-303", "GTTCAGGAT", "998-1000/2400", "898-900/1800", donor + ",splice_donor_5th_base_variant", back, klass),
        ("21:401-403", "GTTCA", "901-903", "801-803", "splice_acceptor_variant," + tract, back, klass),
        ("21:501-503", "GTTCA", "901-903", "801-803", tract, shift, klass),
        ("21:601-603", "GTTCA", "901-903", "801-803", tract, back, bgc.UNEXPLAINED),
        ("21:701-702", "GTTCA", "999-1000", "899-900", donor, back, bgc.UNEXPLAINED),
        ("21:801-803", "GTTCA", "999-?", "899-?", donor, back, bgc.UNEXPLAINED),
        ("21:901-903", "GTTCA", "-", "-", "splice_donor_variant,intron_variant", "intron_variant", bgc.UNEXPLAINED),
        ("21:1001-1003", "GTT", "998-1000", "898-900", donor, shift, bgc.UNEXPLAINED),
        ("21:1101-1103", "-", "998-1000", "898-900", donor, shift, bgc.UNEXPLAINED),
        ("21:1201-1203", "GTTCAGGAT", "998-1000", "898-900", donor + ",splice_donor_5th_base_variant", shift + ",splice_donor_5th_base_variant", bgc.UNEXPLAINED),
        ("21:1301-1303", "GTTCA", "998-1000", "898-900", donor, "stop_gained," + back, bgc.UNEXPLAINED),
        ("21:1401-1403", "GTTCA", "998-1000", "898-900", back, donor, bgc.UNEXPLAINED),
        ("21:1501-1503", "GTTCAG", "998-1000", "898-900", "splice_donor_region_variant,intron_variant," + donor, back, klass),
        ("21:1601-1603", "GTTCAG", "998-1000", "898-900", "intron_variant," + shift, shift, klass),
        ("21:1701-1703", "GTTCAG", "998-1000", "898-900", "intron_variant," + shift, back, bgc.UNEXPLAINED),
        ("21:1801-1803", "GTTCAGGAT", "998-1000", "898-900", "splice_donor_5th_base_variant," + shift, "splice_donor_region_variant," + shift, bgc.UNEXPLAINED),
        ("21:1901-1903", "GTTCAG", "998-1000", "898-900", "intron_variant," + donor, "intron_variant," + back, bgc.UNEXPLAINED),
    ]
    m = _classified(tmp_path, cases)
    assert m["divergence_summary"] == {bgc.SPLICE_TERM_PAST_THE_RECORD_SPAN: 7, bgc.UNEXPLAINED: 12}
    assert bgc.SPLICE_TERM_PAST_THE_RECORD_SPAN == klass

    # A flag set's record-level row (`--summary`) names no feature and carries no position of
    # one transcript, which the class reads: the class names no such row, whatever class the
    # manifest's own run documents the record's Transcript row under.
    corpus = tmp_path / "c"
    m["flag_sets"] = [{"name": "summary", "flags": ["--summary"]}]
    (corpus / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    union = "{}downstream_gene_variant"
    differing = {0: (donor + ",", back + ","), 1: (donor + ",", shift + ","), 6: (donor + ",", back + ",")}
    summary_vep = tmp_path / "summary_vep.txt"
    summary_rs = tmp_path / "summary_rs.txt"
    for path, side in ((summary_vep, 0), (summary_rs, 1)):
        path.write_text("".join(
            _positioned_row(f"rs{i}", loc, allele, "-", union.format(differing.get(i, ("", ""))[side]), "-", "-", ftype="-")
            for i, (loc, allele, *_) in enumerate(cases)
        ), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(summary_vep), "--vep-rs-output", str(summary_rs),
                     "--flag-set", "summary"]) == 0
    fs = json.loads((corpus / "manifest.json").read_text())["flag_sets"][0]
    assert {d["location"]: d["expected_divergence"] for d in fs["divergences"]} == {
        "21:101-103": bgc.UNEXPLAINED, "21:201-203": bgc.UNEXPLAINED, "21:701-702": bgc.UNEXPLAINED,
    }


# Rows of one made-up transcript each, by the record that names them: (VCF POS, ID, REF, ALT,
# Location, allele). An ID of `.` leaves the row under the name Ensembl VEP 116 writes from
# the record's own columns.
_INSERTION_ROWS = [
    (100, ".", "A", "ATGA", "21:100-101", "TGA"),
    (200, ".", "A", "ACG", "21:200-201", "CG"),
    (300, "rs3", "A", "ATGA", "21:300-301", "TGA"),
    (400, "rs4", "AC", "G", "21:400-401", "G"),
    (500, "rs5", "AC", "CGT", "21:500-501", "CGT"),
    (600, "rs6", "AC", "CG", "21:600-601", "CG"),
    (700, ".", "AC", "GTT", "21:700-701", "GTT"),
    (800, ".", "AC", "AGT", "21:800-801", "AGT"),
    (900, ".", "AC", "A,T", "21:900-901", "A"),
    (1000, ".", "A", "AC,T", "21:1000", "AC"),
    (1100, ".", "AC", "A,T", "21:1100-1101", "GTT"),
]
_READ_AWAY_SETS = [
    "", "inframe_insertion,", "inframe_deletion,", "protein_altering_variant,",
    "inframe_insertion,inframe_deletion,", "inframe_insertion,protein_altering_variant,",
    "inframe_deletion,protein_altering_variant,", "inframe_insertion,inframe_deletion,protein_altering_variant,",
]
# A stop codon row: (VCF POS, ID, REF, ALT, Location, allele, cDNA, CDS, reference terms, vep-rs
# terms, vep-rs's cDNA, vep-rs's CDS).
StopCodonRow = tuple[int, str, str, str, str, str, str, str, str, str, str, str]


def _insertion_rows() -> list[StopCodonRow]:
    """The rows of `_INSERTION_ROWS` with the insertion class's first pair of terms and both
    cDNA positions and a CDS start on both sides."""
    rows = []
    for row in _INSERTION_ROWS:
        cdna, cds = ("7054-7055", "6977-6978") if "-" in row[4] else ("7054", "6977")
        rows.append((*row, cdna, cds, "coding_sequence_variant", "stop_retained_variant", cdna, cds))
    return rows


def _read_away_rows() -> list[StopCodonRow]:
    """A deletion's row per set of `_READ_AWAY_SETS` beside the reference's
    `stop_retained_variant` against `frameshift_variant`, with and without `stop_gained`."""
    rows = []
    for i, beside in enumerate(_READ_AWAY_SETS):
        for j, shift in enumerate(("frameshift_variant", "stop_gained,frameshift_variant")):
            pos = 2000 + 100 * (2 * i + j)
            rows.append((pos, f"rs{pos}", "ACGTACGT", "A", f"21:{pos + 1}-{pos + 7}", "-", "961-967", "264-270",
                         beside + "stop_retained_variant", shift, "961-967", "264-270"))
    return rows


def _further_rows() -> list[StopCodonRow]:
    """Rows on which the two rows of a pair, or the record's alleles and the Location, say
    different things: an insertion of two bases whose alleles the first column spells and a
    one-position row under an identifier, each on the three-base class's terms, and for each
    class a pair whose vep-rs row lacks the positions the reference's row has."""
    retained, shift = "stop_retained_variant", "frameshift_variant"
    return [
        (5000, ".", "A", "ACG", "21:5000-5001", "CG", "961-962", "264-265", retained, shift, "961-962", "264-265"),
        (5100, "rs5100", "A", "ACGT", "21:5100", "CGT", "961", "264", retained, shift, "961", "264"),
        (5200, "rs5200", "ACGTACGT", "A", "21:5201-5207", "-", "961-967", "264-270", retained, shift, "-", "-"),
        (5300, "rs5300", "A", "ATGA", "21:5300-5301", "TGA", "7054-7055", "6977-6978", "coding_sequence_variant", retained, "7054-?", "6977-?"),
    ]


def _classified_rows(tmp_path: Path, rows: list[StopCodonRow], elsewhere: bool = False) -> dict[str, str]:
    """Location -> the class `classify` writes for each row. The corpus's records are the rows'
    own; with `elsewhere` each lies 100000 positions past its row, so that no record of the
    run writes a row."""
    corpus = tmp_path / "c"
    records = [f"21\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t." for pos, vid, ref, alt, *_ in rows]
    names = [bgc.record_names(vcf)[0] for vcf in records]
    if elsewhere:
        records = [f"21\t{pos + 100000}\t{vid}\t{ref}\t{alt}\t.\tPASS\t." for pos, vid, ref, alt, *_ in rows]
    _manifest(corpus, [_record(vcf) for vcf in records])
    for path, terms, cdna, cds in ((tmp_path / "vep.txt", 8, 6, 7), (tmp_path / "rs.txt", 9, 10, 11)):
        path.write_text("".join(
            _positioned_row(name, row[4], row[5], f"ENST{i}", row[terms], row[cdna], row[cds]) for i, (name, row) in enumerate(zip(names, rows))
        ), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(tmp_path / "vep.txt"), "--vep-rs-output", str(tmp_path / "rs.txt")]) == 0
    return {d["location"]: d["expected_divergence"] for d in json.loads((corpus / "manifest.json").read_text())["divergences"]}


def test_classify_takes_an_insertion_between_two_positions_alone_for_the_insertion_range_class(tmp_path: Path):
    """`stop_codon_overlap_insertion_range` is the class of an insertion between two neighbouring
    positions and of no wider row: where the first column spells the record's alleles the
    reference allele is empty once the shared first base is dropped, an insertion of two bases
    included; under an identifier the allele has another length than two, so two bases replaced
    by one or by three are read as an insertion. Two bases replaced by three with or without a
    shared first base, a deletion written as a sequence, an allele no alternate names beside a
    reference allele of two and an insertion kept with its shared first base at one position,
    each with its alleles spelt, and a two-base allele under an identifier stay an
    `unexplained_residual`."""
    assert _classified_rows(tmp_path, _insertion_rows()) == {
        **{location: bgc.STOP_CODON_OVERLAP_INSERTION_RANGE for location in ("21:100-101", "21:200-201", "21:300-301", "21:400-401", "21:500-501")},
        **{location: bgc.UNEXPLAINED for location in ("21:600-601", "21:700-701", "21:800-801", "21:900-901", "21:1000", "21:1100-1101")},
    }


def test_classify_takes_the_retained_term_alone_or_beside_one_read_away_term(tmp_path: Path):
    """`stop_retained_three_base_stop_codon_test` is the class of the reference's
    `stop_retained_variant`, alone or beside one of `inframe_insertion`, `inframe_deletion` and
    `protein_altering_variant`, against `frameshift_variant` with or without `stop_gained`.
    Two or three of the three terms beside it stay an `unexplained_residual`."""
    rows = _read_away_rows()
    several = {row[4] for row in rows if row[8].count(",") > 1}
    assert len(rows) == 16 and len(several) == 8
    assert _classified_rows(tmp_path, rows) == {
        row[4]: bgc.UNEXPLAINED if row[4] in several else bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST for row in rows
    }


def _other_class_rows() -> list[StopCodonRow]:
    """Rows of the three stop codon classes that read both rows of a pair and the record's own
    alleles: for each class a pair both of whose rows meet it, a pair whose vep-rs row does
    not, and for the genomic span's class an insertion of two bases whose alleles the first
    column spells and a one-position row under an identifier."""
    retained, coding, lost = "stop_retained_variant", "coding_sequence_variant", "stop_lost"
    region = "splice_region_variant,"
    deleted = ("ACGTACGT", "A")
    return [
        (6000, "rs6000", *deleted, "21:6001-6007", "-", "1705-1711", "1630-1636", region + coding, region + retained, "1705-1711", "1630-1636"),
        (6100, "rs6100", *deleted, "21:6101-6107", "-", "1705-1711", "1630-1636", region + coding, region + retained, "1705-?", "1630-?"),
        (6200, ".", "A", "ACG", "21:6200-6201", "CG", "1705-1706", "1630-1631", region + coding, region + retained, "1705-1706", "1630-1631"),
        (6300, "rs6300", "A", "ACGT", "21:6300", "CGT", "1705", "1630", region + coding, region + retained, "1705", "1630"),
        (7000, "rs7000", *deleted, "21:7001-7007", "-", "?-7057", "-", retained, coding, "?-7057", "-"),
        (7100, "rs7100", *deleted, "21:7101-7107", "-", "?-7057", "-", retained, coding, "7051-7057", "6974-6980"),
        (8000, "rs8000", *deleted, "21:8001-8007", "-", "?-7057", "-", retained, lost, "?-7057", "-"),
        (8100, "rs8100", *deleted, "21:8101-8107", "-", "?-7057", "-", retained, lost, "7051-7057", "6974-6980"),
    ]


def test_classify_reads_both_rows_and_the_records_alleles_for_the_genomic_span_class(tmp_path: Path):
    """`stop_codon_overlap_genomic_span` is the class of a pair both of whose rows show the stop
    codon test made, an insertion of two bases whose alleles the first column spells included.
    A pair whose vep-rs row lacks a cDNA end and a one-position row under an identifier stay an
    `unexplained_residual`."""
    rows = [row for row in _other_class_rows() if row[0] < 7000]
    assert _classified_rows(tmp_path, rows) == {
        "21:6001-6007": bgc.STOP_CODON_OVERLAP_GENOMIC_SPAN, "21:6101-6107": bgc.UNEXPLAINED,
        "21:6200-6201": bgc.STOP_CODON_OVERLAP_GENOMIC_SPAN, "21:6300": bgc.UNEXPLAINED,
    }


def test_classify_reads_both_rows_for_the_untested_class(tmp_path: Path):
    """`stop_retained_without_stop_codon_test` is the class of a pair both of whose rows show the
    stop codon test not made; a pair whose vep-rs row carries both cDNA positions and a CDS
    start stays an `unexplained_residual`."""
    rows = [row for row in _other_class_rows() if 7000 <= row[0] < 8000]
    assert _classified_rows(tmp_path, rows) == {
        "21:7001-7007": bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST, "21:7101-7107": bgc.UNEXPLAINED,
    }


def test_classify_reads_both_rows_for_the_removed_codon_class(tmp_path: Path):
    """`stop_codon_removed_without_stop_codon_test` is the class of a pair neither of whose rows
    carries both cDNA positions and a CDS start; a pair whose vep-rs row carries them stays an
    `unexplained_residual`."""
    rows = [row for row in _other_class_rows() if row[0] >= 8000]
    assert _classified_rows(tmp_path, rows) == {
        "21:8001-8007": bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST, "21:8101-8107": bgc.UNEXPLAINED,
    }


def test_classify_asks_every_row_of_a_tuple_two_rows_carry(tmp_path: Path):
    """A tuple that two rows of an output carry is named under a stop codon class only when
    every reference row against every vep-rs row meets it, as the key comparator's rule takes
    the pair: a row whose first column spells two bases replaced by three beside a row under
    an identifier leaves the tuple an `unexplained_residual` in either order, and so does a
    row without its positions beside one with them, on the reference's side or on vep-rs's.
    Two rows of insertions are named under the insertion class."""
    coding, retained, shift = "coding_sequence_variant", "stop_retained_variant", "frameshift_variant"
    tested, cut = ("961-967", "264-270"), ("961-?", "264-?")
    replaced, identified, inserted = (".", "AC", "GTT"), ("rs1", "AC", "GTT"), (".", "A", "ATGA")
    # (VCF POS, Location, allele, reference terms, vep-rs terms, the reference's rows, vep-rs's rows, class), a row
    # being its record's (ID, REF, ALT) and its two positions.
    cases = [
        (100, "21:100-101", "GTT", coding, retained, [(replaced, tested), (identified, tested)], None, bgc.UNEXPLAINED),
        (200, "21:200-201", "GTT", coding, retained, [(identified, tested), (replaced, tested)], None, bgc.UNEXPLAINED),
        (300, "21:300-301", "TGA", coding, retained, [(inserted, tested), (("rs3", "A", "ATGA"), tested)], None, bgc.STOP_CODON_OVERLAP_INSERTION_RANGE),
        (400, "21:401-407", "-", retained, shift, [(("rs4", "ACGTACGT", "A"), cut), (("rs5", "ACGTACGT", "A"), tested)],
         [(("rs4", "ACGTACGT", "A"), tested), (("rs5", "ACGTACGT", "A"), tested)], bgc.UNEXPLAINED),
        (500, "21:501-507", "-", retained, shift, [(("rs6", "ACGTACGT", "A"), tested), (("rs7", "ACGTACGT", "A"), tested)],
         [(("rs6", "ACGTACGT", "A"), ("-", "-")), (("rs7", "ACGTACGT", "A"), tested)], bgc.UNEXPLAINED),
    ]
    corpus = tmp_path / "c"
    records: dict[str, None] = {}
    sides: dict[str, list[str]] = {"vep": [], "rs": []}
    for i, (pos, location, allele, reference, candidate, reference_rows, candidate_rows, _) in enumerate(cases):
        for side, terms, rows in (("vep", reference, reference_rows), ("rs", candidate, candidate_rows or reference_rows)):
            for (vid, ref, alt), (cdna, cds) in rows:
                vcf = f"21\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t."
                records[vcf] = None
                sides[side].append(_positioned_row(bgc.record_names(vcf)[0], location, allele, f"ENST{i}", terms, cdna, cds))
    _manifest(corpus, [_record(vcf) for vcf in records])
    for side, rows in sides.items():
        (tmp_path / f"{side}.txt").write_text("".join(rows), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(tmp_path / "vep.txt"), "--vep-rs-output", str(tmp_path / "rs.txt")]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert {d["location"]: d["expected_divergence"] for d in m["divergences"]} == {case[1]: case[7] for case in cases}


def test_the_builder_and_the_key_comparator_admit_the_same_rows_to_the_stop_codon_classes(tmp_path: Path):
    """One list of rows through the key comparator's predicates of the five stop codon classes
    and through `classify`: the builder names a row under a class exactly where the
    comparator's rule of that name takes the row's pair, the two rows of a pair read as the
    comparator reads them."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_outputs as cvo

    rows = [*_insertion_rows(), *_read_away_rows(), *_further_rows(), *_other_class_rows()]
    rules = {r.bucket: r.predicate for r in cvo.REFERENCE_REGISTRIES["116.2"].row_pair_rules
             if r.bucket != bgc.SPLICE_TERM_PAST_THE_RECORD_SPAN}
    assert set(rules) == {
        bgc.STOP_CODON_OVERLAP_GENOMIC_SPAN, bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST, bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST,
        bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST, bgc.STOP_CODON_OVERLAP_INSERTION_RANGE,
    }
    comparator = {}
    for i, (pos, vid, ref, alt, location, allele, cdna, cds, reference, candidate, candidate_cdna, candidate_cds) in enumerate(rows):
        name = bgc.record_names(f"21\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.")[0]
        e = cvo.pair_row([name, location, allele, "ENSG1", f"ENST{i}", "Transcript", reference, cdna, cds])
        a = cvo.pair_row([name, location, allele, "ENSG1", f"ENST{i}", "Transcript", candidate, candidate_cdna, candidate_cds])
        comparator[location] = next((bucket for bucket, predicate in rules.items() if predicate(e, a)), bgc.UNEXPLAINED)
    assert _classified_rows(tmp_path, rows) == comparator
    assert set(comparator.values()) == {*rules, bgc.UNEXPLAINED}
    assert {location: comparator[location] for location in ("21:5000-5001", "21:5100", "21:5201-5207", "21:5300-5301")} == {
        "21:5000-5001": bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST, "21:5100": bgc.UNEXPLAINED,
        "21:5201-5207": bgc.UNEXPLAINED, "21:5300-5301": bgc.UNEXPLAINED,
    }


# A pair of rows of the splice term class or of its nearest shapes, with its record: (VCF POS,
# REF, ALT, Location, allele, cDNA, CDS, reference terms, vep-rs terms, vep-rs's cDNA, vep-rs's
# CDS, whether `splice_term_past_the_record_span` takes the pair with the record's alleles,
# whether it takes it without them).
SpliceTermPair = tuple[int, str, str, str, str, str, str, str, str, str, str, bool, bool]


def _splice_term_pairs() -> list[SpliceTermPair]:
    """Pairs of the splice term class and its nearest shapes. With the record's alleles the
    class takes two replaced bases inside an exon, three, three short of the intron on the
    reverse strand's layout with the intron term gone, a span across an exon's end that keeps
    the donor and intron terms, an intron span, a donor region term gone beside a kept donor
    term, three bases whose rows differ by `stop_retained_variant` against `stop_lost`
    without a CDS start as well, and three bases behind a first base the record's two
    alleles share; it refuses three replaced bases whose vep-rs row lacks the region term
    the span gives, a pure insertion, a record of equal lengths, a donor term gone where the
    span holds a base of its window, a further term and other cDNA positions on vep-rs's
    row. Without them it takes the pairs of three positions in the transcript's sequence,
    the one without the region term among them, and refuses two positions, the span across
    the exon's end and the intron spans."""
    donor, region, intron = "splice_donor_variant", "splice_donor_region_variant", "intron_variant"
    fifth, back, shift, coding = "splice_donor_5th_base_variant", "splice_region_variant", "frameshift_variant", "coding_sequence_variant"
    past = f"{donor},{region},{intron}"
    two, three, none = ("999-1000", "899-900"), ("998-1000", "898-900"), ("-", "-")
    return [
        (9000, "AC", "GTTCA", "21:9000-9001", "GTTCA", *two, f"{past},{shift}", f"{shift},{back}", *two, True, False),
        (9100, "ACG", "GTTCAG", "21:9100-9102", "GTTCAG", *three, f"{past},{shift}", f"{shift},{back}", *three, True, True),
        (9200, "ACG", "GTTCAG", "21:9200-9202", "GTTCAG", *three, f"{past},{shift}", shift, *three, False, True),
        (9300, "ACG", "GTTCAGGAT", "21:9300-9302", "GTTCAGGAT", "903-905", "803-805",
         f"splice_acceptor_variant,{intron},{shift}", f"{shift},{back}", "903-905", "803-805", True, True),
        (9400, "ACG", "GTTCAGGA", "21:9400-9402", "GTTCAGGA", "999-?", "899-?",
         f"{donor},{fifth},{coding},{intron}", f"{donor},{coding},{intron}", "999-?", "899-?", True, False),
        (9500, "ACG", "GTTCA", "21:9500-9502", "GTTCA", *none, f"{donor},{fifth},{intron}", f"{donor},{region},{intron}", *none, True, False),
        (9600, "AC", "GTTC", "21:9600-9601", "GTTC", *none, past, f"{donor},{intron}", *none, True, False),
        (9700, "A", "AGTT", "21:9700-9701", "GTT", *two, f"{donor},inframe_insertion", f"inframe_insertion,{back}", *two, False, False),
        (9800, "ACG", "GTT", "21:9800-9802", "GTT", *three, f"{donor},missense_variant", f"missense_variant,{back}", *three, False, False),
        (9900, "AC", "GTTC", "21:9900-9901", "GTTC", "1000-?", "900-?", f"{past},{coding}", f"{coding},{intron},{back}", "1000-?", "900-?", False, False),
        (10000, "ACG", "GTTCAG", "21:10000-10002", "GTTCAG", *three, f"{past},{shift}", f"stop_gained,{shift},{back}", *three, False, False),
        (10100, "ACG", "GTTCAG", "21:10100-10102", "GTTCAG", *three, f"{past},{shift}", f"{shift},{back}", "498-500", "898-900", False, False),
        (10200, "ACG", "GTTCAG", "21:10200-10202", "GTTCAG", "2398-2400", "-",
         f"{past},stop_retained_variant,3_prime_UTR_variant", f"stop_lost,3_prime_UTR_variant,{back}", "2398-2400", "-", True, True),
        (10300, "AACG", "AGTTCAG", "21:10301-10303", "GTTCAG", *three, f"{past},{shift}", f"{shift},{back}", *three, True, True),
    ]


# The readings of a pair's record: (name, whether the first column spells the record's alleles,
# whether the run or the comparison holds the record of a row under an identifier).
_RECORD_READINGS = (("spelt", True, False), ("records", False, True), ("without", False, False))


def _splice_term_rows(pairs: list[SpliceTermPair], spelt: bool) -> list[StopCodonRow]:
    """The pairs as stop codon rows, each under the name Ensembl VEP 116 writes from its
    record's own columns (`spelt`) or under an identifier."""
    return [(pair[0], "." if spelt else f"rs{pair[0]}", *pair[1:11]) for pair in pairs]


def _records_vcf(path: Path, rows: list[StopCodonRow]) -> Path:
    """An input VCF of the rows' records."""
    path.write_text("".join(f"21\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\n" for pos, vid, ref, alt, *_ in rows), encoding="utf-8")
    return path


def test_the_builder_and_the_key_comparator_admit_the_same_rows_to_the_splice_term_class(tmp_path: Path):
    """One list of pairs (`_splice_term_pairs`) through the key comparator's row-aware
    predicates, the five stop codon classes' and `splice_term_past_the_record_span`'s, and
    through `classify`, in each reading of the record: its alleles spelt in the first column,
    a row under an identifier whose record the run holds (the comparator's rows take the same
    records from an input VCF), and a row under an identifier that no record of the run
    writes. In each the builder names a pair under the splice term class exactly where the
    comparator's rule of that name takes it, and under no stop codon class."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_outputs as cvo

    pairs = _splice_term_pairs()
    rules = {r.bucket: r.predicate for r in cvo.REFERENCE_REGISTRIES["116.2"].row_pair_rules}
    assert set(rules) == {
        bgc.STOP_CODON_OVERLAP_GENOMIC_SPAN, bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST, bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST,
        bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST, bgc.STOP_CODON_OVERLAP_INSERTION_RANGE, bgc.SPLICE_TERM_PAST_THE_RECORD_SPAN,
    }
    for reading, spelt, held in _RECORD_READINGS:
        directory = tmp_path / reading
        directory.mkdir()
        rows = _splice_term_rows(pairs, spelt)
        builder = _classified_rows(directory, rows, elsewhere=not (spelt or held))
        sides = [[cvo.pair_row(line.split("\t")[:9]) for line in (directory / name).read_text(encoding="utf-8").splitlines()]
                 for name in ("vep.txt", "rs.txt")]
        if held:
            records = _records_vcf(directory / "records.vcf", rows)
            for side in sides:
                cvo.read_record_references(records, side)
        comparator = {e.location: next((bucket for bucket, predicate in rules.items() if predicate(e, a)), bgc.UNEXPLAINED)
                      for e, a in zip(*sides)}
        assert builder == comparator, reading
        taken = 12 if reading == "without" else 11
        assert comparator == {pair[3]: bgc.SPLICE_TERM_PAST_THE_RECORD_SPAN if pair[taken] else bgc.UNEXPLAINED for pair in pairs}, reading


def test_every_python_reading_of_the_splice_term_class_admits_the_same_pairs(tmp_path: Path):
    """One list of pairs (`_splice_term_pairs`) through each Python reading of
    `splice_term_past_the_record_span`: the key comparator's registry rule on its `pair_row`
    rows, the format comparator's predicate of the class on its default-format entries, the
    structural-variant comparator's registry mask over the two outputs and this builder's
    `classify`. With the record's alleles spelt in the first column each admits the pairs
    the class takes with the record. Under an identifier the rows that are given the pairs'
    records by position (`read_record_references`), the entries the format comparator gives
    them (`read_entry_references`) and the builder whose run holds them admit the same
    pairs; the structural-variant comparator, which is given no input records, admits those
    the class takes without the record, as every reading does where no record writes the
    rows."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    sys.path.insert(0, str(HERE.parent / "validation"))
    import compare_sv_concordance as csc
    import compare_vep_formats as cvf
    import compare_vep_outputs as cvo

    klass = "splice_term_past_the_record_span"
    rule = next((r for r in cvo.REFERENCE_REGISTRIES["116.2"].row_pair_rules if r.bucket == klass), None)
    assert rule is not None, f"the key comparator's 116.2 registry holds no row-aware rule {klass}"
    pairs = _splice_term_pairs()
    with_the_record = {pair[3]: pair[11] for pair in pairs}
    without_the_record = {pair[3]: pair[12] for pair in pairs}
    assert {(pair[11], pair[12]) for pair in pairs} == {(True, True), (True, False), (False, True), (False, False)}
    header = "#" + "\t".join(bgc.DEFAULT_COLUMNS) + "\n"
    for reading, spelt, held in _RECORD_READINGS:
        directory = tmp_path / reading
        directory.mkdir()
        rows = _splice_term_rows(pairs, spelt)
        expected = without_the_record if reading == "without" else with_the_record
        builder = {location: named == klass for location, named in _classified_rows(directory, rows, elsewhere=not (spelt or held)).items()}
        outputs = [directory / "vep.txt", directory / "rs.txt"]
        key_rows = [[cvo.pair_row(line.split("\t")[:9]) for line in path.read_text(encoding="utf-8").splitlines()] for path in outputs]
        entries = [cvf.parse_default(header + path.read_text(encoding="utf-8")) for path in outputs]
        if held:
            records = _records_vcf(directory / "records.vcf", rows)
            for side in key_rows:
                cvo.read_record_references(records, side)
            cvf.read_entry_references(*entries, records)
        key_comparator = {e.location: rule.predicate(e, a) for e, a in zip(*key_rows)}
        format_comparator = {e.location: cvf.TERM_PAIR_CLASSES[klass](e, a) for e, a in zip(*entries)}
        perl, rust = (csc.parse_vep_output(path)[0] for path in outputs)
        _, set_aside, counts = csc.filter_registry_swap_pairs(perl, rust, (rule,), outputs[:1], outputs[1:])
        structural_variant_comparator = {location: location in {t[0] for t in set_aside} for location in expected}
        read_from_the_outputs = with_the_record if spelt else without_the_record
        assert key_comparator == expected, reading
        assert format_comparator == expected, reading
        assert builder == expected, reading
        assert structural_variant_comparator == read_from_the_outputs, reading
        assert counts[klass] == sum(read_from_the_outputs.values()), reading


def test_classify_reads_a_splice_term_row_with_the_alleles_of_the_record_that_writes_it(tmp_path: Path):
    """`classify` reads a row with the alleles of the record of its run that writes it, found
    by position (`record_reference`), where the row's first column spells none. Two replaced
    bases inside an exon take `splice_term_past_the_record_span` under an identifier, behind
    a `chr` prefix on the record, and under the name Ensembl VEP writes from the record's own
    columns. An insertion between the same two positions, whose row has the same columns,
    three replaced bases whose vep-rs row lacks the region term the span gives, and two
    replaced bases beside an insertion of the same allele at their Location stay an
    `unexplained_residual`. A row whose record lies elsewhere is read without it: three
    positions take the class and two do not. A flag set over `variants.vcf` in another
    format reads the manifest's records, and one over an input of its own that file's."""
    klass = bgc.SPLICE_TERM_PAST_THE_RECORD_SPAN
    gone, back = "splice_donor_variant,splice_donor_region_variant,intron_variant,frameshift_variant", "frameshift_variant,splice_region_variant"
    two, three = ("999-1000", "899-900"), ("998-1000", "898-900")
    # (the record's first five columns, the row's name, Location, allele, cDNA and CDS, vep-rs's terms, class)
    cases = [
        ("21\t9000\trs1\tAC\tGTTCA", "rs1", "21:9000-9001", "GTTCA", two, back, klass),
        ("chr21\t9100\trs2\tAC\tGTTCA", "rs2", "21:9100-9101", "GTTCA", two, back, klass),
        ("21\t9200\t.\tAC\tGTTCA", "21_9200_AC/GTTCA", "21:9200-9201", "GTTCA", two, back, klass),
        ("21\t9300\trs4\tA\tAGTTCA", "rs4", "21:9300-9301", "GTTCA", two, back, bgc.UNEXPLAINED),
        ("21\t9400\trs5\tACG\tGTTCAG", "rs5", "21:9400-9402", "GTTCAG", three, "frameshift_variant", bgc.UNEXPLAINED),
        ("21\t9500\trs6\tAC\tGTTCA", "rs6", "21:9500-9501", "GTTCA", two, back, bgc.UNEXPLAINED),
        ("21\t9500\trs7\tA\tAGTTCA", "rs7", "21:9500-9501", "GTTCA", two, back, bgc.UNEXPLAINED),
        ("21\t600\trs8\tACG\tGTTCAG", "rs8", "21:9600-9602", "GTTCAG", three, "frameshift_variant", klass),
        ("21\t700\trs9\tAC\tGTTCA", "rs9", "21:9700-9701", "GTTCA", two, back, bgc.UNEXPLAINED),
    ]
    corpus = tmp_path / "c"
    _manifest(corpus, [_record(f"{vcf}\t.\tPASS\t.") for vcf, *_ in cases])
    vep, rs = tmp_path / "vep.txt", tmp_path / "rs.txt"
    vep.write_text("".join(_positioned_row(name, location, allele, f"ENST{i}", gone, *positions)
                           for i, (_, name, location, allele, positions, _, _) in enumerate(cases)), encoding="utf-8")
    rs.write_text("".join(_positioned_row(name, location, allele, f"ENST{i}", terms, *positions)
                          for i, (_, name, location, allele, positions, terms, _) in enumerate(cases)), encoding="utf-8")
    classify = ["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)]
    assert bgc.main(classify) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert [(d["feature"], d["record_indices"], d["expected_divergence"]) for d in m["divergences"]] == [
        (f"ENST{i}", [i], case[6]) for i, case in enumerate(cases)
    ]
    assert m["divergence_summary"] == {klass: 4, bgc.UNEXPLAINED: 5}

    replaced, inserted = _record("chr21\t9000\trs1\tAC\tGTTCA\t.\tPASS\t."), _record("21\t9000\trs2\tA\tAGTTCA\t.\tPASS\t.")
    assert bgc.record_reference([replaced], "21:9000-9001", "GTTCA") == "AC"
    assert bgc.record_reference([inserted], "chr21:9000-9001", "GTTCA") == ""
    assert bgc.record_reference([replaced, inserted], "21:9000-9001", "GTTCA") is None
    assert bgc.record_reference([replaced], "21:9100-9101", "GTTCA") is None
    assert bgc.record_reference([replaced], "22:9000-9001", "GTTCA") is None

    # A flag set over `variants.vcf` in Ensembl format holds no VCF record: its rows are read
    # with the manifest's records. One over its own VCF is read with that file's lines.
    m["flag_sets"] = [
        {"name": "ensembl", "flags": ["--format", "ensembl"], "input": "variants.ensembl.txt"},
        {"name": "own", "flags": ["--pick"], "input": "own.vcf"},
    ]
    (corpus / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    columns = [vcf.split("\t") for vcf, *_ in cases]
    (corpus / "variants.ensembl.txt").write_text(
        "".join(f"{chrom}\t{pos}\t{int(pos) + len(ref) - 1}\t{ref}/{alt}\t1\t{vid}\n" for chrom, pos, vid, ref, alt in columns), encoding="utf-8")
    assert bgc.main([*classify, "--flag-set", "ensembl"]) == 0
    ensembl = json.loads((corpus / "manifest.json").read_text())["flag_sets"][0]
    assert [(d["record_indices"], d["expected_divergence"]) for d in ensembl["divergences"]] == [([i], case[6]) for i, case in enumerate(cases)]
    (corpus / "own.vcf").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t8000\trs20\tAC\tGTTCA\t.\tPASS\t.\n21\t8100\trs21\tA\tAGTTCA\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    own = [("rs20", "21:8000-8001", klass), ("rs21", "21:8100-8101", bgc.UNEXPLAINED)]
    vep.write_text("".join(_positioned_row(name, location, "GTTCA", f"ENST{i}", gone, *two) for i, (name, location, _) in enumerate(own)), encoding="utf-8")
    rs.write_text("".join(_positioned_row(name, location, "GTTCA", f"ENST{i}", back, *two) for i, (name, location, _) in enumerate(own)), encoding="utf-8")
    assert bgc.main([*classify, "--flag-set", "own"]) == 0
    fs = json.loads((corpus / "manifest.json").read_text())["flag_sets"][1]
    assert [(d["location"], d["record_indices"], d["expected_divergence"]) for d in fs["divergences"]] == [
        (location, [i], named) for i, (_, location, named) in enumerate(own)
    ]


def test_classify_leaves_the_stop_codon_test_out_for_want_of_the_cds_start(tmp_path: Path):
    """Ensembl VEP leaves its stop codon test out for a span without a cDNA position at one end
    or a CDS start (Utils/VariationEffect.pm:1508: `return 0 unless $cdna_start && $cdna_end &&
    $cds_start`). A deletion from the 5' UTR to the stop codon's last base, whose CDS position
    reads `?-N`, takes the two classes of a test that is not made; a row with the CDS start and
    no CDS end is one whose test is made."""
    m = _classified(tmp_path, [
        ("21:1201-9405", "-", "21-99", "?-69", "start_lost,stop_retained_variant", "stop_lost,start_lost", bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST),
        ("21:1301-9405", "-", "21-99/105", "?-69/69", "start_lost,stop_retained_variant", "stop_lost,start_lost", bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST),
        ("21:1401-9405", "-", "21-99", "?-69", "start_lost,stop_retained_variant", "start_lost", bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST),
        ("21:1501-1505", "-", "21-99", "60-?", "stop_retained_variant,3_prime_UTR_variant", "stop_lost,3_prime_UTR_variant", bgc.UNEXPLAINED),
        ("21:1601-1605", "-", "21-99", "60-?", "stop_retained_variant,3_prime_UTR_variant", "coding_sequence_variant,3_prime_UTR_variant", bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST),
    ])
    assert m["divergence_summary"] == {
        bgc.STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST: 2, bgc.STOP_RETAINED_WITHOUT_STOP_CODON_TEST: 1,
        bgc.STOP_RETAINED_THREE_BASE_STOP_CODON_TEST: 1, bgc.UNEXPLAINED: 1,
    }


def test_select_focus_terms_take_more_exemplars(tmp_path: Path):
    """A combination containing a focus term takes --focus-k exemplars; every other
    combination takes --k."""
    ref = tmp_path / "ref.txt"
    rows = ""
    vcf_lines = ""
    for i in range(1, 6):
        rows += _row(f"rs{i}", f"21:{i * 100}", "G", "ENST1", "inframe_insertion")
        rows += _row(f"rs{i}", f"21:{i * 100}", "G", "ENST2", "intron_variant")
        vcf_lines += f"21\t{i * 100}\trs{i}\tA\tG\t.\tPASS\t.\n"
    ref.write_text(rows, encoding="utf-8")
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + vcf_lines, encoding="utf-8")
    out = tmp_path / "corpus"
    rc = bgc.main(["select", "--release", "115", "--assembly", "GRCh37", "--suite", f"s01={ref}={vcf}", "--out", str(out),
                   "--k", "1", "--focus-term", "inframe_insertion", "--focus-k", "4", "--seed", "1", "--multi-allelic-per-suite", "0", "--covers", "symbol"])
    assert rc == 0
    m = json.loads((out / "manifest.json").read_text())
    assert m["focus_terms"] == ["inframe_insertion"] and m["focus_exemplars_per_combination"] == 4
    assert m["combinations"]["inframe_insertion"]["exemplar_records"] == 4
    # A record's combinations are the ones it was chosen for, so intron_variant keeps --k.
    assert m["combinations"]["intron_variant"]["exemplar_records"] == 1
    assert 4 <= len(m["records"]) <= 5
    # Without --flag or --fasta-name the manifest carries neither key.
    assert "flags" not in m and "fasta" not in m


def test_select_with_a_cache_takes_the_exemplar_that_costs_the_fewest_bytes(tmp_path: Path):
    """Two records show the same combination in the same stratum; with --cache the one whose
    reachable transcripts deflate smaller is chosen, and the manifest records the flank."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    heavy = {"stable_id": "ENST_HEAVY", "start": "100", "end": "900", "variation_effect_feature_cache": {"peptide": "MKV" * 4000}}
    light = {"stable_id": "ENST_LIGHT", "start": "500100", "end": "500900"}
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([heavy, light]))
    ref = tmp_path / "ref.txt"
    ref.write_text(_row("rs_heavy", "21:500", "G", "ENST_HEAVY", "missense_variant") + _row("rs_light", "21:500500", "G", "ENST_LIGHT", "missense_variant"), encoding="utf-8")
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t500\trs_heavy\tA\tG\t.\tPASS\t.\n21\t500500\trs_light\tA\tG\t.\tPASS\t.\n", encoding="utf-8")
    common = ["select", "--release", "116", "--assembly", "GRCh38", "--suite", f"sel={ref}={vcf}", "--k", "1", "--multi-allelic-per-suite", "0", "--covers", "symbol"]
    picked = set()
    for seed in range(1, 7):
        out = tmp_path / f"plain{seed}"
        assert bgc.main(common + ["--out", str(out), "--seed", str(seed)]) == 0
        picked.add(json.loads((out / "manifest.json").read_text())["records"][0]["vcf"].split("\t")[2])
        out = tmp_path / f"costed{seed}"
        assert bgc.main(common + ["--out", str(out), "--seed", str(seed), "--cache", str(cache)]) == 0
        m = json.loads((out / "manifest.json").read_text())
        assert [r["vcf"].split("\t")[2] for r in m["records"]] == ["rs_light"]
        assert m["cache_cost_flank"] == 5000
    assert picked == {"rs_heavy", "rs_light"}, "without a cache the tie is broken by the seed"
    cost = bgc.CacheCost(cache, 5000)
    far = bgc.Row("sel", "x", "21:600000", "G", "ENST_LIGHT", "Transcript", "intron_variant", "1", "", ("-", "-", "-"))
    assert cost.row_cost(far) == 0, "nothing within the flank"
    mate = bgc.Row("sel", "x", "21:600000", "N[21:400[", "ENST_HEAVY", "Transcript", "feature_truncation", "1", "", ("-", "-", "-"))
    assert cost.row_cost(mate) == cost.span_cost("21", 400, 600000) > cost.span_cost("21", 500500, 500500) > 0


def test_select_keeps_an_exemplar_of_every_variant_class(tmp_path: Path):
    """A class no chosen exemplar shows (here a breakend with inserted sequence, outranked by
    a cheaper SNV for the same combination) still gets one record, and the manifest lists
    the classes with their exemplar rows."""
    ref = tmp_path / "ref.txt"
    ref.write_text(
        _row("rs1", "21:100", "G", "ENST1", "feature_truncation")
        + _row("bnd1", "21:300000", "[21:300900[TTGC", "ENST2", "feature_truncation"),
        encoding="utf-8",
    )
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\trs1\tA\tG\t.\tPASS\t.\n"
        "21\t300000\tbnd1\tC\t[21:300900[TTGC\t.\tPASS\tSVTYPE=BND\n",
        encoding="utf-8",
    )
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([
        {"stable_id": "ENST1", "start": "50", "end": "150"},
        {"stable_id": "ENST2", "start": "299000", "end": "301000", "variation_effect_feature_cache": {"peptide": "MKV" * 400}},
    ]))
    out = tmp_path / "corpus"
    rc = bgc.main(["select", "--release", "116", "--assembly", "GRCh38", "--suite", f"sel={ref}={vcf}", "--out", str(out),
                   "--k", "1", "--seed", "1", "--multi-allelic-per-suite", "0", "--cache", str(cache), "--covers", "symbol"])
    assert rc == 0
    m = json.loads((out / "manifest.json").read_text())
    by_id = {r["vcf"].split("\t")[2]: r for r in m["records"]}
    assert by_id["rs1"]["combinations"] == ["feature_truncation"]
    assert by_id["bnd1"]["combinations"] == ["feature_truncation"]
    assert m["variant_classes"] == {"breakend_insertion": 1, "snv": 1}
    assert m["combinations"]["feature_truncation"]["exemplar_records"] == 2
    # Without the class guarantee the cheaper SNV alone would carry the combination.
    rows = list(bgc.iter_rows(ref, "sel"))
    chosen = bgc.choose_exemplars({"feature_truncation": rows}, 1, random.Random(1), 0, (), None, bgc.CacheCost(cache, 5000).row_cost)
    assert [r.uploaded_variation for r in chosen["feature_truncation"]] == ["rs1"]


def _matrix(tmp_path: Path, *items: str) -> Path:
    """A parity matrix with the given items (the seven columns of tests/parity/matrix.tsv)."""
    path = tmp_path / "matrix.tsv"
    rows = ["\t".join([item, "Config.pm", "implemented", "open", "", "", ""]) for item in items]
    path.write_text("item\tperl_reference\ttarget\tstate\ttest1\ttest2\tinterim\n" + "".join(r + "\n" for r in rows), encoding="utf-8")
    return path


def _two_record_suite(tmp_path: Path) -> tuple[Path, Path]:
    ref = tmp_path / "ref.txt"
    ref.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant") + _row("rs2", "21:200", "T", "ENST2", "intron_variant"), encoding="utf-8")
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n", encoding="utf-8")
    return ref, vcf


def _select(ref: Path, vcf: Path, out: Path, *extra: str) -> int:
    return bgc.main(["select", "--release", "116", "--assembly", "GRCh38", "--suite", f"s01={ref}={vcf}", "--out", str(out),
                     "--k", "1", "--seed", "1", "--multi-allelic-per-suite", "0", *extra])


def test_select_writes_covers(tmp_path: Path):
    """`--covers` (repeatable) is the manifest's `covers`, in order, placed where the hand-written
    lists sit: after the run options, before `suites`."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    assert _select(ref, vcf, out, "--covers", "symbol", "--covers", "max_af", "--covers", "symbol", "--flag=--symbol") == 0
    m = json.loads((out / "manifest.json").read_text())
    assert m["covers"] == ["symbol", "max_af"]
    keys = list(m)
    assert keys.index("flags") < keys.index("covers") < keys.index("suites")
    assert "flag_sets" not in m


def test_select_refuses_a_corpus_without_covers(tmp_path: Path, capsys):
    """A corpus that names no matrix row certifies nothing, so `select` writes neither file."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    assert _select(ref, vcf, out) == 1
    assert capsys.readouterr().err == (
        "ERROR: [build_golden_corpus] select writes no manifest without --covers: "
        "name the parity-matrix rows the corpus exercises\n"
    )
    assert not out.exists()


def test_select_refuses_a_cache_without_transcripts(tmp_path: Path):
    """A `--cache` with no `transcripts/` directory is refused by name before any reference is
    read, so `select` writes nothing."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    cache = tmp_path / "cache"
    cache.mkdir()
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, out, "--covers", "symbol", "--cache", str(cache))
    assert str(refusal.value) == f"ERROR: [build_golden_corpus] {cache} has no transcripts/ directory"
    assert not out.exists()
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, out, "--covers", "symbol", "--cache", str(tmp_path / "absent"))
    assert str(refusal.value) == f"ERROR: [build_golden_corpus] {tmp_path / 'absent'} has no transcripts/ directory"
    # A reference that cannot be read is refused the same way, not with a traceback.
    with pytest.raises(SystemExit) as refusal:
        _select(tmp_path / "absent.txt", vcf, out, "--covers", "symbol")
    assert str(refusal.value).startswith(f"ERROR: [build_golden_corpus] {tmp_path / 'absent.txt'} is not a readable text file (")


def test_select_preserves_covers_on_regeneration(tmp_path: Path):
    """Regenerating over an existing manifest keeps its `covers` (union with any new `--covers`),
    its `masks`, its `reference_warnings`, its `plugin_files` and its `flag_sets`, in that order
    before `suites`, so a regenerated corpus keeps certifying the rows it did, expecting the
    run-level warnings it did and naming the plugin data files its reference was generated over;
    a set's divergence keys name the old records by ordinal, so they go the way the manifest's own
    do; a mask class `classify` does not know, or a `reference_warnings` that is not a list of
    texts, is refused before anything is written."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    matrix = _matrix(tmp_path, "symbol", "max_af", "pick")
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    assert bgc.main(["flag-set", str(out), "--name", "pick", "--flag=--pick", "--covers", "pick", "--matrix", str(matrix)]) == 0
    m = json.loads((out / "manifest.json").read_text())
    entry = dict(m["flag_sets"][0])
    m["masks"] = ["regulatory_feature_rows"]
    notice = "could not fetch sequence for tandem repeats; consequence calculation for tandem repeats is less precise when using --offline without --fasta"
    m["reference_warnings"] = [notice]
    plugin_files = {"LoFTEE": ["plugin_data/human_ancestor.fa.gz", "plugin_data/gerp.bw"]}
    m["plugin_files"] = plugin_files
    m["divergences"] = [{"location": "21:100", "record_indices": [0]}]
    m["flag_sets"][0].update({"divergences": [{"location": "21:100", "record_indices": [0]}], "divergence_summary": {}})
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    assert _select(ref, vcf, out) == 0, "the existing list satisfies the covers requirement"
    m = json.loads((out / "manifest.json").read_text())
    assert m["covers"] == ["symbol"]
    assert m["masks"] == ["regulatory_feature_rows"]
    assert m["reference_warnings"] == [notice]
    assert m["plugin_files"] == plugin_files
    assert m["flag_sets"] == [entry] and "divergences" not in m
    assert _select(ref, vcf, out, "--covers", "max_af", "--covers", "symbol") == 0
    m = json.loads((out / "manifest.json").read_text())
    assert m["covers"] == ["symbol", "max_af"]
    assert m["masks"] == ["regulatory_feature_rows"] and m["flag_sets"] == [entry]
    assert m["reference_warnings"] == [notice]
    assert m["plugin_files"] == plugin_files
    keys = list(m)
    assert (keys.index("covers") < keys.index("masks") < keys.index("reference_warnings")
            < keys.index("plugin_files") < keys.index("flag_sets") < keys.index("suites"))
    m["masks"] = ["splice_rows"]
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    (out / "variants.vcf").unlink()
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, out)
    assert str(refusal.value) == (f"ERROR: [build_golden_corpus] {out / 'manifest.json'} masks names unknown classes "
                                  f"['splice_rows'] (the classes are {', '.join(bgc.MASK_CLASSES)})")
    assert json.loads((out / "manifest.json").read_text()) == m and not (out / "variants.vcf").exists(), "nothing is written"
    m["masks"] = ["regulatory_feature_rows"]
    m["reference_warnings"] = notice
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    with pytest.raises(SystemExit, match="reference_warnings is not a list of texts"):
        _select(ref, vcf, out)
    assert json.loads((out / "manifest.json").read_text()) == m and not (out / "variants.vcf").exists(), "nothing is written"


def test_select_matrix_and_empty_names(tmp_path: Path, capsys):
    """`select --matrix` refuses a `--covers` name that is not a row before any file is written,
    and an empty name is refused with or without a matrix."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    matrix = _matrix(tmp_path, "symbol")
    assert _select(ref, vcf, out, "--covers", "symbol", "--covers", "ghost", "--matrix", str(matrix)) == 1
    assert capsys.readouterr().err == f"ERROR: [build_golden_corpus] covers names ghost, which is not a row of {matrix}\n"
    assert not out.exists()
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, out, "--covers", "")
    assert str(refusal.value) == "ERROR: [build_golden_corpus] --covers wants a non-empty row name"
    assert not out.exists()
    assert _select(ref, vcf, out, "--covers", "symbol", "--matrix", str(matrix)) == 0


def _one_stratum_suite(tmp_path: Path) -> tuple[Path, Path]:
    """Two SNVs of one combination on one strand: one stratum, so `--k 1` samples one of them."""
    ref = tmp_path / "ref.txt"
    ref.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant") + _row("rs2", "21:200", "T", "ENST1", "missense_variant"), encoding="utf-8")
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n", encoding="utf-8")
    return ref, vcf


def _record_ids(out: Path) -> list[str]:
    return [r["vcf"].split("\t")[2] for r in json.loads((out / "manifest.json").read_text())["records"]]


def test_select_include_carries_the_named_record_outside_the_sampling(tmp_path: Path, capsys):
    """`--include sid=name` carries the record whatever the sampling chose, with every combination
    its rows show, and records the names under `included_records`; a name with no reference row is
    refused before anything is written, as are a suite not among `--suite` and a malformed argument."""
    ref, vcf = _one_stratum_suite(tmp_path)
    out = tmp_path / "corpus"
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    sampled = _record_ids(out)
    assert len(sampled) == 1
    other = ({"rs1", "rs2"} - set(sampled)).pop()
    assert _select(ref, vcf, out, "--covers", "symbol", "--include", f"s01={other}") == 0
    m = json.loads((out / "manifest.json").read_text())
    assert sorted(_record_ids(out)) == ["rs1", "rs2"]
    assert m["included_records"] == {"s01": [other]}
    combinations_by_id = {r["vcf"].split("\t")[2]: r["combinations"] for r in m["records"]}
    assert combinations_by_id[other] == ["missense_variant"]
    assert "kept_records" not in m
    keys = list(m)
    assert keys.index("exemplars_per_combination") < keys.index("included_records") < keys.index("max_span")
    (out / "variants.vcf").unlink()
    assert _select(ref, vcf, out, "--include", "s01=ghost", "--include", "s01=rs1") == 1
    assert capsys.readouterr().err == "ERROR: [build_golden_corpus] --include s01=ghost has no row in the s01 reference\n"
    assert not (out / "variants.vcf").exists() and json.loads((out / "manifest.json").read_text()) == m, "nothing is written"
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, out, "--include", "s02=rs1")
    assert str(refusal.value) == "ERROR: [build_golden_corpus] --include s02=rs1 names suite 's02', which is not among --suite"
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, out, "--include", "rs1")
    assert str(refusal.value) == "ERROR: [build_golden_corpus] --include wants sid=name, got 'rs1'"


def test_select_keep_records_keeps_every_record_of_the_manifest(tmp_path: Path, capsys):
    """A regeneration under `--keep-records` keeps each record of the manifest being regenerated,
    rowed again from its suite's reference, records their count as `kept_records` and carries the
    manifest's `included_records`; the summary counts the distinct named records; a kept record with
    no reference row is refused before anything is written, a kept record of a suite not among
    `--suite`, a run with no manifest to regenerate and a malformed `included_records` are refused."""
    ref, vcf = _one_stratum_suite(tmp_path)
    out = tmp_path / "corpus"
    assert _select(ref, vcf, out, "--covers", "symbol", "--k", "2") == 0
    assert sorted(_record_ids(out)) == ["rs1", "rs2"]
    assert _select(ref, vcf, out) == 0
    assert len(_record_ids(out)) == 1, "the sampling alone drops one of the two"
    assert _select(ref, vcf, out, "--k", "2") == 0
    capsys.readouterr()
    assert _select(ref, vcf, out, "--keep-records") == 0
    assert capsys.readouterr().out.rstrip().endswith(", named records 2")
    m = json.loads((out / "manifest.json").read_text())
    assert sorted(_record_ids(out)) == ["rs1", "rs2"]
    assert m["kept_records"] == 2 and "included_records" not in m
    assert _select(ref, vcf, out, "--keep-records", "--max-records", "1") == 0
    assert sorted(_record_ids(out)) == ["rs1", "rs2"], "a kept record survives the record cap"
    capsys.readouterr()
    assert _select(ref, vcf, out, "--keep-records", "--include", "s01=rs1") == 0
    assert capsys.readouterr().out.rstrip().endswith(", named records 2"), "a record both kept and included counts once"
    assert json.loads((out / "manifest.json").read_text())["included_records"] == {"s01": ["rs1"]}
    assert _select(ref, vcf, out, "--keep-records") == 0
    m = json.loads((out / "manifest.json").read_text())
    assert m["included_records"] == {"s01": ["rs1"]}, "a regeneration that keeps the records keeps their included names"
    assert _select(ref, vcf, out, "--keep-records", "--include", "s01=rs2") == 0
    assert json.loads((out / "manifest.json").read_text())["included_records"] == {"s01": ["rs1", "rs2"]}, "carried and new names merge"
    assert _select(ref, vcf, out) == 0
    assert "included_records" not in json.loads((out / "manifest.json").read_text()), "a sampling run carries no included names"
    m["included_records"] = {"s01": ["rs1"], "s02": ["rs9"]}
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, out, "--keep-records")
    assert str(refusal.value) == "ERROR: [build_golden_corpus] carried include s02=rs9 names suite 's02', which is not among --suite"
    m["included_records"] = {"s01": ["rs1", "rs9"]}
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    (out / "variants.vcf").unlink()
    assert _select(ref, vcf, out, "--keep-records") == 1
    assert capsys.readouterr().err == "ERROR: [build_golden_corpus] carried include s01=rs9 has no row in the s01 reference\n"
    assert not (out / "variants.vcf").exists() and json.loads((out / "manifest.json").read_text()) == m, "nothing is written"
    for malformed in (["rs1"], {"s01": []}, {"": ["rs1"]}, {"s01": [""]}):
        m["included_records"] = malformed
        (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
        with pytest.raises(SystemExit) as refusal:
            _select(ref, vcf, out, "--keep-records")
        assert str(refusal.value) == f"ERROR: [build_golden_corpus] {out / 'manifest.json'} included_records is not a map of suite to record names"
    m["included_records"] = {"s01": ["rs1"]}
    m["records"].append({"suite": "s01", "vcf": "21\t300\trs3\tA\tG\t.\tPASS\t.", "combinations": []})
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    assert _select(ref, vcf, out, "--keep-records") == 1
    assert capsys.readouterr().err == "ERROR: [build_golden_corpus] kept record 21 300 rs3 A G has no row in the s01 reference\n"
    assert not (out / "variants.vcf").exists() and json.loads((out / "manifest.json").read_text()) == m, "nothing is written"
    m["records"][-1] = {"suite": "s02", "vcf": "21\t300\trs3\tA\tG\t.\tPASS\t.", "combinations": []}
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, out, "--keep-records")
    assert str(refusal.value) == "ERROR: [build_golden_corpus] kept record 21 300 rs3 A G names suite 's02', which is not among --suite"
    fresh = tmp_path / "fresh"
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, fresh, "--covers", "symbol", "--keep-records")
    assert str(refusal.value) == f"ERROR: [build_golden_corpus] --keep-records wants a manifest to regenerate in {fresh}"
    assert not fresh.exists()


def test_plugin_files_shape_is_held_at_every_read():
    """A corpus that compares a built-in plugin with its Ensembl VEP module names the data
    files the plugin reads under `plugin_files`, as paths inside the corpus, so a reference
    is tied to them; `manifest_shape_error` refuses every shape `plugin_files_shape_error`
    does, with its reason."""
    good = {"LoFTEE": ["plugin_data/human_ancestor.fa.gz", "plugin_data/gerp.bw"]}
    assert bgc.plugin_files_shape_error(good) is None
    assert bgc.manifest_shape_error({"flags": [], "plugin_files": good}) is None
    for value, reason in (
        (["plugin_data/gerp.bw"], "plugin_files is not a map of plugin names"),
        ({"": ["plugin_data/gerp.bw"]}, "plugin_files is not a map of plugin names"),
        ({"LoFTEE": "plugin_data/gerp.bw"}, "plugin_files for LoFTEE is not a list of paths"),
        ({"LoFTEE": ["plugin_data/gerp.bw", ""]}, "plugin_files for LoFTEE is not a list of paths"),
        ({"LoFTEE": ["/data/gerp.bw"]}, "plugin_files for LoFTEE names paths outside the corpus ['/data/gerp.bw']"),
        ({"LoFTEE": ["../gerp.bw"]}, "plugin_files for LoFTEE names paths outside the corpus ['../gerp.bw']"),
        ({"LoFTEE": ["plugin data/gerp.bw"]},
         "plugin_files for LoFTEE names paths outside the corpus ['plugin data/gerp.bw']"),
    ):
        assert bgc.plugin_files_shape_error(value) == reason
        assert bgc.manifest_shape_error({"flags": [], "plugin_files": value}) == reason


def test_manifest_shapes_are_refused_by_every_subcommand(tmp_path: Path):
    """A `covers` that is not a list of names, or a `flag_sets` that is not a list of named
    objects, is refused wherever a manifest is read, with the reason, before anything is
    written; `--flag-set ''` names no set rather than the manifest's own list."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    matrix = _matrix(tmp_path, "symbol")
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    manifest = out / "manifest.json"
    good = manifest.read_text()
    commands = [
        ["select", "--release", "116", "--assembly", "GRCh38", "--suite", f"s01={ref}={vcf}", "--out", str(out),
         "--k", "1", "--seed", "1", "--multi-allelic-per-suite", "0"],
        ["check", "--corpus", str(out), "--suite", f"s01={ref}={vcf}"],
        ["covers", str(out), "--list", "--matrix", str(matrix)],
        ["flag-set", str(out), "--list", "--matrix", str(matrix)],
        ["classify", "--corpus", str(out), "--vep-default", str(ref), "--vep-rs-output", str(ref)],
    ]
    shapes = [
        ('"covers": "symbol"', "covers is not a list of names"),
        ('"covers": [["symbol"]]', "covers is not a list of names"),
        ('"covers": ["symbol", ""]', "covers is not a list of names"),
        ('"covers": ["symbol"], "masks": "regulatory_feature_rows"', "masks is not a list of names"),
        ('"covers": ["symbol"], "masks": [["regulatory_feature_rows"]]', "masks is not a list of names"),
        ('"covers": ["symbol"], "masks": ["order_only", "splice_rows"]',
         f"masks names unknown classes ['splice_rows'] (the classes are {', '.join(bgc.MASK_CLASSES)})"),
        ('"covers": ["symbol"], "reference_warnings": "no FASTA"', "reference_warnings is not a list of texts"),
        ('"covers": ["symbol"], "reference_warnings": ["no FASTA", ""]', "reference_warnings is not a list of texts"),
        ('"covers": ["symbol"], "plugin_files": ["plugin_data/gerp.bw"]', "plugin_files is not a map of plugin names"),
        ('"covers": ["symbol"], "plugin_files": {"LoFTEE": ["../gerp.bw"]}',
         "plugin_files for LoFTEE names paths outside the corpus ['../gerp.bw']"),
        ('"covers": ["symbol"], "flag_sets": "pick"', "flag_sets is not a list"),
        ('"covers": ["symbol"], "flag_sets": ["pick"]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"covers": []}]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "default"}]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "-x"}]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a"}, {"name": "a"}]', "names flag set a twice"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "covers": "x"}]', "flag set a covers is not a list of names"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "masks": "order_only"}]', "flag set a masks is not a list of names"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "masks": ["order_only", "splice_rows"]}]',
         f"flag set a masks names unknown classes ['splice_rows'] (the classes are {', '.join(bgc.MASK_CLASSES)})"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "reference_warnings": "no FASTA"}]', "flag set a reference_warnings is not a list of texts"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "reference_warnings": [["no FASTA"]]}]', "flag set a reference_warnings is not a list of texts"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "reference_warnings": ["no FASTA", " "]}]', "flag set a reference_warnings is not a list of texts"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "flags": "--pick"}]', "flag set a flags is not a list of strings"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "input": 5}]', "flag set a input is not a path inside the corpus directory"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "input": "../x.vcf"}]', "flag set a input is not a path inside the corpus directory"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "expected_dir": ""}]', "flag set a expected_dir is not a path inside the corpus directory"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "stdin": "true"}]', "flag set a stdin is not true or false"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "unordered": 1}]', "flag set a unordered is not true or false"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "formats": []}]', "flag set a formats is not a non-empty subset of default,tab,vcf,json"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "formats": "vcf"}]', "flag set a formats is not a non-empty subset of default,tab,vcf,json"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "stdin": true, "flags": ["--pick"]}]', "flag set a runs on stdin without --format among its flags"),
        ('"covers": ["symbol"], "cache_layout": "perl"', 'cache_layout is "perl", not "ensembl"'),
        ('"covers": ["symbol"], "cache_layout": null', 'cache_layout is null, not "ensembl"'),
        ('"covers": ["symbol"], "cache_layout": "ensembl", "release": "116"', "cache_layout ensembl wants release and assembly"),
        ('"covers": ["symbol"], "fasta_in_cache": "yes"', 'fasta_in_cache is "yes", not a boolean'),
        ('"covers": ["symbol"], "fasta_in_cache": true', "fasta_in_cache wants fasta"),
        ('"covers": ["symbol"], "keep_samples": "8"', "keep_samples is not a count"),
        ('"covers": ["symbol"], "keep_samples": -1', "keep_samples is not a count"),
        ('"covers": ["symbol"], "keep_samples": true', "keep_samples is not a count"),
        ('"covers": ["symbol"], "keep_samples": 2', "samples does not name 2 kept samples"),
        ('"covers": ["symbol"], "keep_samples": 2, "samples": ["a"]', "samples does not name 2 kept samples"),
        ('"covers": ["symbol"], "keep_samples": 1, "samples": [""]', "samples does not name 1 kept samples"),
        ('"covers": ["symbol"], "required_genotype_classes": "het"', "required_genotype_classes is not a list of genotype classes"),
        ('"covers": ["symbol"], "required_genotype_classes": ["het", "ghost"]', "required_genotype_classes is not a list of genotype classes"),
    ]
    for body, reason in shapes:
        manifest.write_text('{"records": [], ' + body + "}", encoding="utf-8")
        for argv in commands:
            with pytest.raises(SystemExit) as refusal:
                bgc.main(argv)
            assert str(refusal.value) == f"ERROR: [build_golden_corpus] {manifest} {reason}", argv[0]
    # Every mask class passes at both levels, the harness-only classes included, and so does
    # an empty list on a set (its run masks nothing).
    accepted = json.loads(good)
    accepted["masks"] = ["order_only", "regulatory_feature_rows"]
    accepted["flag_sets"] = [
        {"name": "a", "masks": ["order_only"]},
        {"name": "b", "masks": []},
        {"name": "c", "masks": list(bgc.HARNESS_ONLY_MASK_CLASSES)},
    ]
    assert set(bgc.HARNESS_ONLY_MASK_CLASSES) == {
        "unloaded_module_plugin", "conservation_short_name_spelling", "pc_description_text",
        "file_cache_writes_no_source_field",
    }
    manifest.write_text(json.dumps(accepted), encoding="utf-8")
    assert bgc.main(commands[1]) == 0
    manifest.write_text("[]", encoding="utf-8")
    with pytest.raises(SystemExit) as refusal:
        bgc.main(commands[1])
    assert str(refusal.value) == f"ERROR: [build_golden_corpus] {manifest} is not a JSON object"
    manifest.write_text("{", encoding="utf-8")
    with pytest.raises(SystemExit) as refusal:
        bgc.main(commands[0])
    assert str(refusal.value).startswith(f"ERROR: [build_golden_corpus] {manifest} is not readable JSON (")
    manifest.write_text(good, encoding="utf-8")
    with pytest.raises(SystemExit) as refusal:
        bgc.main(["covers", str(out), "--flag-set", "", "--list", "--matrix", str(matrix)])
    assert str(refusal.value) == "ERROR: [build_golden_corpus] the manifest has no flag set ''"
    with pytest.raises(SystemExit) as refusal:
        bgc.main(["covers", str(out), "--add", "symbol", "--matrix", str(tmp_path / "none.tsv")])
    assert str(refusal.value) == f"ERROR: [build_golden_corpus] {tmp_path / 'none.tsv'} is not a file: name the parity matrix with --matrix"
    with pytest.raises(SystemExit) as refusal:
        bgc.main(["check", "--corpus", str(out), "--suite", f"s01={ref}={vcf}", "--matrix", ""])
    assert str(refusal.value) == "ERROR: [build_golden_corpus] . is not a file: name the parity matrix with --matrix"
    (tmp_path / "binary.tsv").write_bytes(b"item\t\xff\xfe\n")
    with pytest.raises(SystemExit) as refusal:
        bgc.main(["covers", str(out), "--add", "symbol", "--matrix", str(tmp_path / "binary.tsv")])
    assert str(refusal.value).startswith(f"ERROR: [build_golden_corpus] {tmp_path / 'binary.tsv'} is not a readable text file (")
    assert manifest.read_text() == good


def test_covers_subcommand_adds_a_row_and_refuses_a_non_row(tmp_path: Path, capsys):
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    matrix = _matrix(tmp_path, "symbol", "max_af", "flag_matrix_concordance")
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    before = (out / "manifest.json").read_text()
    assert bgc.main(["covers", str(out), "--add", "no_such_row", "--add", "max_af", "--matrix", str(matrix)]) == 1
    assert capsys.readouterr().err == f"ERROR: [build_golden_corpus] covers names no_such_row, which is not a row of {matrix}\n"
    assert (out / "manifest.json").read_text() == before, "a refused run writes nothing"
    assert bgc.main(["covers", str(out), "--add", "flag_matrix_concordance", "--add", "symbol", "--matrix", str(matrix)]) == 0
    assert bgc.main(["covers", str(out), "--list", "--matrix", str(matrix)]) == 0
    assert capsys.readouterr().out == "symbol\nflag_matrix_concordance\n"
    after = (out / "manifest.json").read_text()
    # Only the list changed: the manifest round-trips byte for byte around it.
    assert after.replace('"covers": [\n  "symbol",\n  "flag_matrix_concordance"\n ]', '"covers": [\n  "symbol"\n ]') == before
    assert bgc.main(["flag-set", str(out), "--name", "af", "--flag=--max_af", "--matrix", str(matrix)]) == 0
    before = (out / "manifest.json").read_text()
    assert bgc.main(["covers", str(out), "--flag-set", "af", "--add", "ghost", "--add", "max_af", "--matrix", str(matrix)]) == 1
    assert capsys.readouterr().err == f"ERROR: [build_golden_corpus] flag set af covers names ghost, which is not a row of {matrix}\n"
    assert (out / "manifest.json").read_text() == before
    assert bgc.main(["covers", str(out), "--flag-set", "af", "--add", "max_af", "--matrix", str(matrix)]) == 0
    m = json.loads((out / "manifest.json").read_text())
    assert m["flag_sets"][0]["covers"] == ["max_af"] and m["covers"] == ["symbol", "flag_matrix_concordance"]
    with pytest.raises(SystemExit):
        bgc.main(["covers", str(out), "--flag-set", "none", "--list", "--matrix", str(matrix)])


def test_check_refuses_covers_outside_the_matrix(tmp_path: Path, capsys):
    """`check --matrix` reads every `covers` list, the manifest's and each flag set's."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    matrix = _matrix(tmp_path, "symbol", "max_af")
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    check = ["check", "--corpus", str(out), "--suite", f"s01={ref}={vcf}", "--matrix", str(matrix)]
    assert bgc.main(check) == 0
    m = json.loads((out / "manifest.json").read_text())
    m["covers"].append("retired_row")
    m["flag_sets"] = [{"name": "af", "flags": ["--max_af"], "covers": ["max_af", "ghost"]}]
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    assert bgc.main(check) == 1
    captured = capsys.readouterr()
    assert captured.err == (
        f"ERROR: [build_golden_corpus] covers names retired_row, which is not a row of {matrix}\n"
        f"ERROR: [build_golden_corpus] flag set af covers names ghost, which is not a row of {matrix}\n"
    )
    assert captured.out.splitlines()[-1] == "observed 2, covered by corpus 2, missing 0"
    assert bgc.main(check[:-2]) == 0, "without --matrix the names are not read"


def test_family_columns_stratify_exemplars(tmp_path: Path, capsys):
    """Two rows of one combination and one (class, strand) stratum differ only in the family
    key (HGVSc present or absent; IMPACT value): with `--k 1` alone one survives, with
    `--family-columns` every stratum keeps an exemplar, the manifest records the key, and
    `check --family-columns` passes on that corpus and names the stratum a corpus lacks."""
    ref = tmp_path / "ref.txt"
    ref.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant", "IMPACT=MODERATE;STRAND=1;HGVSc=ENST1.1:c.5C>A")
        + _row("rs2", "21:200", "T", "ENST1", "missense_variant", "IMPACT=MODERATE;STRAND=1")
        + _row("rs3", "21:300", "A", "ENST1", "missense_variant", "IMPACT=HIGH;STRAND=1;HGVSc=ENST1.1:c.9C>A")
        + _row("rs4", "21:400", "-", "ENST1", "frameshift_variant", "IMPACT=HIGH;STRAND=1;HGVSc=ENST1.1:c.12del"),
        encoding="utf-8",
    )
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n21\t300\trs3\tC\tA\t.\tPASS\t.\n21\t400\trs4\tCA\tC\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    plain = tmp_path / "plain"
    assert _select(ref, vcf, plain, "--covers", "symbol") == 0
    m = json.loads((plain / "manifest.json").read_text())
    assert m["combinations"]["missense_variant"]["exemplar_records"] == 1 and "family_columns" not in m
    family = tmp_path / "family"
    assert _select(ref, vcf, family, "--covers", "symbol", "--family-columns", "HGVSc,IMPACT=") == 0
    m = json.loads((family / "manifest.json").read_text())
    assert m["family_columns"] == ["HGVSc", "IMPACT="]
    assert m["family_strata"] == {
        "HGVSc+;IMPACT=HIGH;deletion": 1, "HGVSc+;IMPACT=HIGH;snv": 1, "HGVSc+;IMPACT=MODERATE;snv": 1, "HGVSc-;IMPACT=MODERATE;snv": 1,
    }
    assert m["combinations"]["missense_variant"]["exemplar_records"] == 3
    check = ["check", "--suite", f"s01={ref}={vcf}", "--family-columns", "HGVSc,IMPACT="]
    assert bgc.main([*check, "--corpus", str(family)]) == 0
    assert capsys.readouterr().out.splitlines()[-1] == "family strata observed 4, covered by corpus 4, missing 0"
    assert bgc.main([*check, "--corpus", str(plain)]) == 1
    out = capsys.readouterr().out.splitlines()
    assert out[-3] == "family strata observed 4, covered by corpus 2, missing 2"
    # The three missense rows tie on every ordering key, so `--seed 1` alone decides which one the
    # plain corpus keeps: rs2, whose stratum is the one snv stratum the corpus covers.
    assert [r["vcf"].split("\t")[2] for r in json.loads((plain / "manifest.json").read_text())["records"]] == ["rs2", "rs4"]
    assert out[-2:] == [
        "  MISSING family stratum HGVSc+;IMPACT=HIGH;snv",
        "  MISSING family stratum HGVSc+;IMPACT=MODERATE;snv",
    ]
    with pytest.raises(SystemExit):
        bgc.main([*check[:-1], "HGVSc;IMPACT", "--corpus", str(plain)])
    # The record cap keeps a stratum's sole exemplar as it keeps a combination's: four strata
    # survive a cap of three, where the same cap over four records chosen without the key
    # (`--k 3` takes the three missense rows) leaves three.
    capped = tmp_path / "capped"
    assert _select(ref, vcf, capped, "--covers", "symbol", "--family-columns", "HGVSc,IMPACT=", "--max-records", "3") == 0
    assert len(json.loads((capped / "manifest.json").read_text())["records"]) == 4
    plain_capped = tmp_path / "plain-capped"
    assert _select(ref, vcf, plain_capped, "--covers", "symbol", "--k", "3", "--max-records", "3") == 0
    assert len(json.loads((plain_capped / "manifest.json").read_text())["records"]) == 3


def test_custom_match_key_covers_every_class(tmp_path: Path, capsys):
    """`custom_match=CE/CO` reads the columns a `--custom` file's `type=exact` source (CE) and
    `type=overlap` source (CO) write (`<name>` the matched record names, `<name>_<field>` one
    value per record, comma-joined: OutputFactory.pm:1297-1310 at release/116.2) and classes a
    row as one of the five classes; `select --family-columns` keeps an exemplar per class and
    `check --family-columns` names the classes a corpus lacks. The exact source matches by
    allele (File/VCF.pm:405-420) and keeps the matching allele's value of an allele-specific
    field (File/VCF.pm:227-245), where the overlap source lists every allele's (:247-260): a
    field with more values than records whose exact-source value is shorter is the multi-allelic
    record's mark; two matched records with one value each, a FILTER list, or a value with a
    comma the exact source lists whole are not."""
    assert bgc.parse_family_columns("HGVSc,custom_match=CE/CO") == ["HGVSc", "custom_match=CE/CO"]
    for bad in ("custom_match=CE", "custom_match=CE/CO/X", "custom_match=/CO"):
        with pytest.raises(SystemExit):
            bgc.parse_family_columns(bad)
    klass = lambda extra, loc="21:100", allele="G": bgc.custom_match_class(bgc.parse_extra(extra), "CE", "CO", loc, allele)
    seen = {
        klass("CE=exact_1;CO=exact_1;CE_AF=0.10;CO_AF=0.10"),
        klass("CO=overlap_1;CO_AF=0.20"),
        klass("IMPACT=MODIFIER"),
        klass("CE=multi_1;CO=multi_1;CE_AF=0.30;CO_AF=0.30,0.40"),
        klass("CO=exact_1,overlap_1;CO_AF=0.10,0.20", "21:500-900", "deletion"),
    }
    assert seen == set(bgc.CUSTOM_MATCH_CLASSES) == {"exact", "overlap", "miss", "multi_allelic", "sv_span"}
    assert klass("CE=exact_1;CO=exact_1;CE_AF=0.10;CO_AF=0.10") == "exact"
    assert klass("CO=overlap_1;CO_AF=0.20") == "overlap"
    assert klass("IMPACT=MODIFIER") == "miss"
    assert klass("CE=multi_1;CO=multi_1;CE_AF=0.30;CO_AF=0.30,0.40") == "multi_allelic"
    assert klass("CO=exact_1,overlap_1;CO_AF=0.10,0.20", "21:500-900", "deletion") == "sv_span"
    assert klass("CO=sv_1", "21:500", "A[21:900[") == "sv_span"
    assert klass("CE=sv_1;CO=sv_1", "21:500-900", "deletion") == "exact", "an exact structural match is exact"
    assert klass("CE=exact_2;CO=exact_2,overlap_2;CE_AF=0.5;CO_AF=0.5,0.6") == "exact", "two records, one value each"
    assert klass("CE=exact_3;CO=exact_3;CE_FILTER=q10,s50;CO_FILTER=q10,s50") == "exact", "FILTER lists filters, not alleles"
    assert klass("CE=exact_4;CO=exact_4;CE_CLNSIG=a,b;CO_CLNSIG=a,b") == "exact", "a comma the exact source lists whole"
    assert klass("CE=-;CO=overlap_5") == "overlap", "a dash is an absent value"
    assert bgc.family_classes(bgc.parse_extra("CO=o;IMPACT=HIGH"), ["custom_match=CE/CO", "IMPACT="], "21:100", "G") == (
        "custom_match=overlap", "IMPACT=HIGH",
    )

    ref = tmp_path / "ref.txt"
    ref.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant", "IMPACT=MODERATE;STRAND=1;CE=exact_1;CO=exact_1;CE_AF=0.10;CO_AF=0.10")
        + _row("rs2", "21:200", "T", "ENST1", "missense_variant", "IMPACT=MODERATE;STRAND=1;CO=overlap_1;CO_AF=0.20")
        + _row("rs3", "21:300", "A", "ENST1", "missense_variant", "IMPACT=MODERATE;STRAND=1")
        + _row("rs4", "21:400", "C", "ENST1", "missense_variant", "IMPACT=MODERATE;STRAND=1;CE=multi_1;CO=multi_1;CE_AF=0.30;CO_AF=0.30,0.40")
        + _row("rs5", "21:500-900", "deletion", "ENST1", "transcript_ablation", "IMPACT=HIGH;STRAND=1;CO=exact_1,overlap_1;CO_AF=0.10,0.20"),
        encoding="utf-8",
    )
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n21\t300\trs3\tC\tA\t.\tPASS\t.\n"
        "21\t400\trs4\tA\tC\t.\tPASS\t.\n21\t499\trs5\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=900\n",
        encoding="utf-8",
    )
    every = tmp_path / "every"
    assert _select(ref, vcf, every, "--covers", "custom", "--family-columns", "custom_match=CE/CO") == 0
    m = json.loads((every / "manifest.json").read_text())
    assert m["family_columns"] == ["custom_match=CE/CO"]
    assert m["family_strata"] == {
        "custom_match=exact;snv": 1, "custom_match=miss;snv": 1, "custom_match=multi_allelic;snv": 1,
        "custom_match=overlap;snv": 1, "custom_match=sv_span;symbolic_deletion": 1,
    }
    check = ["check", "--suite", f"s01={ref}={vcf}", "--family-columns", "custom_match=CE/CO"]
    assert bgc.main([*check, "--corpus", str(every)]) == 0
    assert capsys.readouterr().out.splitlines()[-1] == "family strata observed 5, covered by corpus 5, missing 0"
    plain = tmp_path / "plain"
    assert _select(ref, vcf, plain, "--covers", "custom") == 0
    assert bgc.main([*check, "--corpus", str(plain)]) == 1
    out = capsys.readouterr().out.splitlines()
    assert out[-4] == "family strata observed 5, covered by corpus 2, missing 3"
    # The four missense rows tie on every ordering key, so the seed alone decides which one the
    # plain corpus keeps; the three classes it does not keep are the missing strata.
    snv_class = {"rs1": "exact", "rs2": "overlap", "rs3": "miss", "rs4": "multi_allelic"}
    kept = {r["vcf"].split("\t")[2] for r in json.loads((plain / "manifest.json").read_text())["records"]}
    kept_snv = [snv_class[name] for name in kept if name in snv_class]
    assert len(kept_snv) == 1 and "rs5" in kept
    missing = sorted(f"custom_match={c};snv" for c in snv_class.values() if c != kept_snv[0])
    assert out[-3:] == [f"  MISSING family stratum {stratum}" for stratum in missing]



def test_refseq_family_columns_stratify_exemplars(tmp_path: Path, capsys):
    """The RefSeq and merged cache columns stratify a corpus: `REFSEQ_MATCH=` by the code list
    (one stratum per distinct list, so every `rseq_*` code keeps an exemplar), `SOURCE=` and
    `BAM_EDIT=` by value, `REFSEQ_OFFSET` by presence (Ensembl VEP prints it only when the
    misalignment offset is non-zero, OutputFactory.pm:1769 at release/116.2) and
    `GIVEN_REF/USED_REF` by whether the two agree (`agreement_class`: `same`, `differ`, or `-`
    when the run wrote neither, as under `--use_given_ref`). Six missense rows of one (class,
    strand) stratum differ in those columns alone: `--k 1` keeps one, the family key keeps every
    stratum, `check --family-columns` passes on that corpus and names the strata the plain one
    lacks."""
    assert bgc.parse_family_columns("REFSEQ_MATCH=,SOURCE=,BAM_EDIT=,REFSEQ_OFFSET,GIVEN_REF/USED_REF") == [
        "REFSEQ_MATCH=", "SOURCE=", "BAM_EDIT=", "REFSEQ_OFFSET", "GIVEN_REF/USED_REF",
    ]
    for bad in ("GIVEN_REF/", "/USED_REF", "A/B/C", "GIVEN_REF/USED_REF="):
        with pytest.raises(SystemExit):
            bgc.parse_family_columns(bad)
    agree = lambda extra: bgc.agreement_class(bgc.parse_extra(extra), "GIVEN_REF", "USED_REF")
    assert agree("GIVEN_REF=A;USED_REF=A") == "same"
    assert agree("GIVEN_REF=A;USED_REF=G") == "differ"
    assert agree("IMPACT=MODIFIER") == "-" and agree("GIVEN_REF=A") == "-" and agree("GIVEN_REF=A;USED_REF=-") == "-"
    assert set(bgc.AGREEMENT_CLASSES) == {"same", "differ", "-"}
    assert bgc.family_classes(
        bgc.parse_extra("REFSEQ_MATCH=rseq_mrna_nonmatch,rseq_5p_mismatch;BAM_EDIT=OK;GIVEN_REF=C;USED_REF=T"),
        ["REFSEQ_MATCH=", "SOURCE=", "BAM_EDIT=", "REFSEQ_OFFSET", "GIVEN_REF/USED_REF"],
    ) == ("REFSEQ_MATCH=rseq_mrna_nonmatch,rseq_5p_mismatch", "SOURCE=-", "BAM_EDIT=OK", "REFSEQ_OFFSET-", "GIVEN_REF/USED_REF=differ")

    ref = tmp_path / "ref.txt"
    base = "IMPACT=MODERATE;STRAND=1;"
    ref.write_text(
        _row("rs1", "21:100", "G", "NM_1.1", "missense_variant", base + "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;BAM_EDIT=OK;GIVEN_REF=A;USED_REF=A")
        + _row("rs2", "21:200", "T", "NM_1.1", "missense_variant", base + "REFSEQ_MATCH=rseq_mrna_nonmatch,rseq_5p_mismatch;BAM_EDIT=OK;GIVEN_REF=C;USED_REF=T")
        + _row("rs3", "21:300", "A", "NM_1.1", "missense_variant", base + "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;BAM_EDIT=FAILED;GIVEN_REF=C;USED_REF=C")
        + _row("rs4", "21:400", "C", "NM_1.1", "missense_variant", base + "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;BAM_EDIT=OK;GIVEN_REF=A;USED_REF=A;REFSEQ_OFFSET=-3")
        + _row("rs5", "21:500", "G", "NM_1.1", "missense_variant", base + "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;SOURCE=RefSeq;BAM_EDIT=OK;GIVEN_REF=A;USED_REF=A")
        + _row("rs6", "21:600", "T", "ENST1", "missense_variant", base + "SOURCE=Ensembl;GIVEN_REF=A;USED_REF=A"),
        encoding="utf-8",
    )
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\trs1\tA\tG\t.\tPASS\t.\n21\t200\trs2\tC\tT\t.\tPASS\t.\n21\t300\trs3\tC\tA\t.\tPASS\t.\n"
        "21\t400\trs4\tA\tC\t.\tPASS\t.\n21\t500\trs5\tA\tG\t.\tPASS\t.\n21\t600\trs6\tA\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    key = "REFSEQ_MATCH=,SOURCE=,BAM_EDIT=,REFSEQ_OFFSET,GIVEN_REF/USED_REF"
    every = tmp_path / "every"
    assert _select(ref, vcf, every, "--covers", "refseq", "--family-columns", key) == 0
    m = json.loads((every / "manifest.json").read_text())
    assert m["family_columns"] == key.split(",")
    assert m["family_strata"] == {
        "REFSEQ_MATCH=-;SOURCE=Ensembl;BAM_EDIT=-;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=same;snv": 1,
        "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;SOURCE=-;BAM_EDIT=FAILED;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=same;snv": 1,
        "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;SOURCE=-;BAM_EDIT=OK;REFSEQ_OFFSET+;GIVEN_REF/USED_REF=same;snv": 1,
        "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;SOURCE=-;BAM_EDIT=OK;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=same;snv": 1,
        "REFSEQ_MATCH=rseq_mrna_match,rseq_ens_match_wt;SOURCE=RefSeq;BAM_EDIT=OK;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=same;snv": 1,
        "REFSEQ_MATCH=rseq_mrna_nonmatch,rseq_5p_mismatch;SOURCE=-;BAM_EDIT=OK;REFSEQ_OFFSET-;GIVEN_REF/USED_REF=differ;snv": 1,
    }
    assert m["combinations"]["missense_variant"]["exemplar_records"] == 6
    check = ["check", "--suite", f"s01={ref}={vcf}", "--family-columns", key]
    assert bgc.main([*check, "--corpus", str(every)]) == 0
    assert capsys.readouterr().out.splitlines()[-1] == "family strata observed 6, covered by corpus 6, missing 0"
    plain = tmp_path / "plain"
    assert _select(ref, vcf, plain, "--covers", "refseq") == 0
    assert len(json.loads((plain / "manifest.json").read_text())["records"]) == 1
    assert bgc.main([*check, "--corpus", str(plain)]) == 1
    out = capsys.readouterr().out.splitlines()
    assert out[-6] == "family strata observed 6, covered by corpus 1, missing 5"
    assert all(line.startswith("  MISSING family stratum ") for line in out[-5:])


CUSTOM_FIXTURE_VCF = (
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "21\t100\trs1\tA\tG\t.\tPASS\t.\n"  # exact
    "21\t200\trs2\tTCA\tT\t.\tPASS\t.\n"  # overlap
    "21\t300\trs3\tC\tA,T\t.\tPASS\t.\n"  # miss
    "21\t300\trs3b\tC\tG\t.\tPASS\t.\n"  # the locus is written once
    "21\t400\trs4\tA\tAT\t.\tPASS\t.\n"  # multi_allelic
    "21\t500\trs5\tG\tC\t.\tPASS\t.\n"  # exact again
    "21\t599\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=900\n"  # sv_span, inside 600..899
    "21\t1000\tbnd1\tA\tA[22:5000[\t.\tPASS\tSVTYPE=BND\n"  # no span
    "22\t50\trs6\tT\tC\t.\tPASS\t.\n"  # overlap, the next chromosome
)


def test_make_custom_fixtures_records_follow_the_class_cycle():
    """A corpus's distinct sequence loci take the classes exact, overlap, miss, multi_allelic in
    turn (a second variant at a locus gets no record), an overlap record keeps POS and REF and
    changes the first ALT's last base (`TCA>T` to `TCA>A`, so `type=exact` sees another allele
    and `type=overlap` the same span: File/VCF.pm:385-420 at release/116.2), a multi-allelic
    record lists the first ALT then that allele with one AF each, a structural variant with END
    gets one single-base record inside its span at a position no sequence variant holds, a
    breakend none; records sort by chromosome then position for tabix."""
    records = [line.split("\t") for line in CUSTOM_FIXTURE_VCF.splitlines() if not line.startswith("#")]
    sites, counts = mcf.derive_sites(records)
    assert counts == {"exact": 2, "overlap": 2, "miss": 1, "multi_allelic": 1, "sv_span": 1}
    assert [s.vcf_line() for s in sites] == [
        "21\t100\texact_1\tA\tG\t.\tPASS\tAF=0.001;CLASS=exact",
        "21\t200\toverlap_1\tTCA\tA\t.\tPASS\tAF=0.002;CLASS=overlap",
        "21\t400\tmulti_1\tA\tAT,AA\t.\tPASS\tAF=0.003,0.503;CLASS=multi_allelic",
        "21\t500\texact_2\tG\tC\t.\tPASS\tAF=0.004;CLASS=exact",
        "21\t600\tsv_1\tN\tA\t.\tPASS\tAF=0.006;CLASS=sv_span",
        "22\t50\toverlap_2\tT\tA\t.\tPASS\tAF=0.005;CLASS=overlap",
    ]
    # The BED score is the record's ordinal, the AF in thousandths.
    assert [s.bed_line() for s in sites][:2] == ["21\t99\t100\texact_1\t1\t+", "21\t199\t202\toverlap_1\t2\t+"]
    assert [s.bedgraph_line() for s in sites][1] == "21\t199\t200\t0.002"
    # The UCSC copies spell the chromosome as UCSC does and score the second record with its
    # ordinal in hundredths minus one (phyloP) and one minus its AF (phastCons).
    assert sites[1].bed_line(mcf.ucsc_chrom(sites[1].chrom)) == "chr21\t199\t202\toverlap_1\t2\t+"
    assert sites[1].bedgraph_line("chr21", sites[1].phylop()) == "chr21\t199\t200\t-0.98"
    assert sites[1].bedgraph_line("chr21", sites[1].phastcons()) == "chr21\t199\t200\t0.998"
    assert mcf.ucsc_chrom("21") == "chr21" and mcf.ucsc_chrom("MT") == "chrM" and mcf.ucsc_chrom("chr1") == "chr1"
    # A promoter window five bases either side of every record, every third without a gene id,
    # the fifth with a second window twenty bases either side, and the enhancer over the first
    # record, sorted by chromosome and start; a window starts at base 1 at the earliest.
    assert mcf.promoter_lines(sites) == [
        "21\t.\tpromoter_window\t95\t105\t.\t+\t.\tID=PW1;gene_id=ENSG90000000001",
        "21\t.\tenhancer\t95\t105\t.\t+\t.\tID=EH1",
        "21\t.\tpromoter_window\t195\t207\t.\t+\t.\tID=PW2;gene_id=ENSG90000000002",
        "21\t.\tpromoter_window\t395\t405\t.\t+\t.\tID=PW3",
        "21\t.\tpromoter_window\t495\t505\t.\t+\t.\tID=PW4;gene_id=ENSG90000000004",
        "21\t.\tpromoter_window\t580\t620\t.\t+\t.\tID=PW5w;gene_id=ENSG90000000005",
        "21\t.\tpromoter_window\t595\t605\t.\t+\t.\tID=PW5;gene_id=ENSG90000000005",
        "22\t.\tpromoter_window\t45\t55\t.\t+\t.\tID=PW6",
    ]
    early = [mcf.Site("1", 3, "exact_1", "A", ("G",), ("0.001",), "exact")]
    assert mcf.promoter_lines(early)[0] == "1\t.\tpromoter_window\t1\t8\t.\t+\t.\tID=PW1;gene_id=ENSG90000000001"
    # A single-base SNV's other allele avoids REF as well as the original ALT.
    assert mcf.other_allele("G", "A") == "C" and mcf.other_allele("A", "G") == "C" and mcf.other_allele("C", "A") == "G"
    assert mcf.other_allele("T", "TCA") == "A" and mcf.other_allele("AT", "A") == "AA"
    assert mcf.structural_end(10, "N", "<DEL>", "SVTYPE=DEL;END=20") == 20
    assert mcf.structural_end(10, "N", "<DEL>", "SVTYPE=DEL;SVLEN=-5") == 15
    assert mcf.structural_end(10, "ACGTACGT", "<CN0>", ".") == 17
    assert mcf.structural_end(10, "N", "N[22:5[", "SVTYPE=BND") is None
    assert mcf.structural_end(10, "N", ".N", "") is None
    assert mcf.contig_length("GRCh38", "chr21", 5) == 46709983 and mcf.contig_length("GRCh37", "21", 5) == 48129895
    assert mcf.contig_length("GRCh38", "NC_000021.9", 5000) == 5000


def test_make_custom_fixtures_writes_indexed_files(tmp_path: Path, capsys):
    """The tool writes the VCF, BED and GFF sources as BGZF (the gzip magic with the extra-field
    flag, `1f 8b 08 04`) with their tabix indexes beside them, the bigWig with its magic,
    records the tool versions and every file's sha256 under provenance.json's custom_fixtures
    key while keeping the record's other keys, and refuses to overwrite without --force."""
    missing = [tool for tool in mcf.TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"{', '.join(missing)} not on PATH")
    corpus = tmp_path / "corpus"
    corpus.mkdir()
    (corpus / "variants.vcf").write_text(CUSTOM_FIXTURE_VCF, encoding="utf-8")
    (corpus / "manifest.json").write_text(json.dumps({"assembly": "GRCh38"}), encoding="utf-8")
    (corpus / "provenance.json").write_text(json.dumps({"vep_image": "x", "input": "hand-built records"}), encoding="utf-8")
    assert mcf.main(["--corpus", str(corpus)]) == 0
    assert capsys.readouterr().out.strip().endswith("loci per class: exact 2, overlap 2, miss 1, multi_allelic 1, sv_span 1")
    out = corpus / "custom"
    for name in ("sites.vcf.gz", "regions.bed.gz", "promoters.gff.gz", "clinvar_source.vcf.gz", "regions_chr.bed.gz"):
        assert (out / name).read_bytes()[:4] == b"\x1f\x8b\x08\x04", name
        assert (out / f"{name}.tbi").is_file(), name
        assert (out / f"{name}.tbi").read_bytes()[:2] == b"\x1f\x8b", name
    with gzip.open(out / "clinvar_source.vcf.gz", "rt", encoding="utf-8") as fh:
        clinvar_lines = fh.read().splitlines()
    assert clinvar_lines[:2] == ["##fileformat=VCFv4.2", "##source=ClinVar"]
    with gzip.open(out / "regions_chr.bed.gz", "rt", encoding="utf-8") as fh:
        assert fh.read().splitlines()[:2] == ["chr21\t99\t100\texact_1\t1\t+", "chr21\t199\t202\toverlap_1\t2\t+"]
    for name in (mcf.PHYLOP_BW, mcf.PHASTCONS_BW):
        assert (out / name).read_bytes()[:4] == b"\x26\xfc\x8f\x88", name
    with gzip.open(out / "promoters.gff.gz", "rt", encoding="utf-8") as fh:
        gff_lines = fh.read().splitlines()
    assert gff_lines[:2] == ["##gff-version 3", "21\t.\tpromoter_window\t95\t105\t.\t+\t.\tID=PW1;gene_id=ENSG90000000001"]
    assert len(gff_lines) == 9
    assert (out / "scores.bw").read_bytes()[:4] == b"\x26\xfc\x8f\x88"
    with gzip.open(out / "sites.vcf.gz", "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2" and "##contig=<ID=21>" in lines and "##contig=<ID=22>" in lines
    assert [line.split("\t")[2] for line in lines if not line.startswith("#")] == [
        "exact_1", "overlap_1", "multi_1", "exact_2", "sv_1", "overlap_2",
    ]
    with gzip.open(out / "regions.bed.gz", "rt", encoding="utf-8") as fh:
        assert fh.read().splitlines()[0] == "21\t99\t100\texact_1\t1\t+"
    assert sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_file()) == sorted(mcf.FILES), "no scratch file survives"
    record = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))
    assert record["vep_image"] == "x" and record["input"] == "hand-built records"
    fixtures = record["custom_fixtures"]
    assert fixtures["generator"] == "scripts/golden/make_custom_fixtures.py"
    assert fixtures["variants_vcf_sha256"] == hashlib.sha256(CUSTOM_FIXTURE_VCF.encode()).hexdigest()
    assert fixtures["assembly"] == "GRCh38"
    assert fixtures["loci"] == {"exact": 2, "overlap": 2, "miss": 1, "multi_allelic": 1, "sv_span": 1}
    assert sorted(fixtures["files"]) == sorted(f"custom/{name}" for name in mcf.FILES)
    for name, digest in fixtures["files"].items():
        assert hashlib.sha256((corpus / name).read_bytes()).hexdigest() == digest, name
    assert sorted(fixtures["tools"]) == ["bgzip", "bigtools", "tabix"]
    assert all(version.split()[-1][0].isdigit() for version in fixtures["tools"].values()), fixtures["tools"]
    with pytest.raises(SystemExit, match="pass --force"):
        mcf.main(["--corpus", str(corpus)])
    assert mcf.main(["--corpus", str(corpus), "--force"]) == 0
    assert json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["custom_fixtures"]["files"] == fixtures["files"], \
        "a regeneration from the same input writes the same bytes"


def test_genotype_labels_follow_vep_genotype_reading():
    """Each GT's classes follow how Ensembl VEP reads a genotype for `--individual`: the bits
    split on `|` when the GT has one, else on `/` (BaseVCF4.pm 944-949; the `\\` of
    Parser/VCF.pm 627 never meets one, so `1\\2` is a single bit and no allele index), `.` bits
    dropped before translation (BaseVCF4.pm 946-949), so `./1` is a one-allele genotype (HOM,
    not HET) and `./.` yields nothing, an allele containing `*` outside the non-reference set
    (Parser/VCF.pm 635) so an individual whose only non-reference allele is `*` or `<*>` is
    reference-homozygous, a haploid `1` HOM, and an index read as a number (`00` is `0`)."""
    labels = bgc.genotype_labels
    assert labels("0|1", ["G"], "r") == {"het"}
    assert labels("1/0", ["G"], "r") == {"het"}
    assert labels("1|1", ["G"], "r") == {"hom_alt"}
    assert labels("0/0", ["G"], "r") == {"hom_ref"}
    assert labels("./.", ["G"], "r") == {"missing"}
    assert labels(".", ["G"], "r") == {"missing"}
    assert labels("./1", ["G"], "r") == {"partially_missing", "hom_alt"}
    assert labels(".|0", ["G"], "r") == {"partially_missing", "hom_ref"}
    assert labels("1", ["G"], "r") == {"haploid", "hom_alt"}
    assert labels("0", ["G"], "r") == {"haploid", "hom_ref"}
    assert labels("1/2", ["G", "T"], "r") == {"multi_allelic"}
    assert labels("0/1/2", ["G", "T"], "r") == {"multi_allelic"}
    assert labels("1/1/1", ["G", "T"], "r") == {"hom_alt"}
    assert labels("0/00", ["G"], "r") == {"hom_ref"}
    assert labels("01/1", ["G", "T"], "r") == {"hom_alt"}
    assert labels("0/2", ["G", "*"], "r") == {"star", "hom_ref"}
    assert labels("2/2", ["G", "*"], "r") == {"star", "hom_ref"}
    assert labels("1/2", ["G", "*"], "r") == {"star", "het"}
    assert labels("2", ["G", "*"], "r") == {"haploid", "star", "hom_ref"}
    assert labels("./2", ["G", "*"], "r") == {"partially_missing", "star", "hom_ref"}
    assert labels("0/2", ["G", "<*>"], "r") == {"star", "hom_ref"}
    for gt in ("0/3", "a/b", "0/-1", "1/", "", "1\\2", "0/1|2"):
        with pytest.raises(SystemExit) as refusal:
            labels(gt, ["G", "T"], "21 100 rs1 A G,T")
        assert str(refusal.value) == f"ERROR: [build_golden_corpus] record 21 100 rs1 A G,T has a GT {gt!r} outside its alleles"


def test_genotype_labels_refuse_a_non_ascii_digit():
    """A GT bit is an ASCII allele index. `str.isdigit` alone accepts a superscript `²` or a
    circled `①`, for which `int` raises, and an Arabic-Indic `٣`, which `int` reads as 3 (inside
    three alleles); each is refused as a bit outside the record's alleles, never a traceback,
    through `genotype_labels` and through a record's kept sample."""
    record = "21 100 rs1 A G,T,C"
    for gt in ("0/²", "①", "٣/0"):
        with pytest.raises(SystemExit) as refusal:
            bgc.genotype_labels(gt, ["G", "T", "C"], record)
        assert str(refusal.value) == f"ERROR: [build_golden_corpus] record {record} has a GT {gt!r} outside its alleles"
    with pytest.raises(SystemExit) as refusal:
        bgc.record_genotype_classes(["21", "100", "rs1", "A", "G,T,C", ".", "PASS", ".", "GT", "0/0", "0/²"], 2)
    assert str(refusal.value) == f"ERROR: [build_golden_corpus] record {record} has a GT '0/²' outside its alleles"


def test_record_genotype_classes_add_the_record_shapes():
    """A record shows every class any kept sample shows, plus the shapes of the whole kept set:
    every genotype missing (Ensembl VEP's per-individual list is then empty and the record
    yields nothing, Parser/VCF.pm 370-376), no non-reference carrier (every individual copy
    is reference-homozygous, 640-650), and one reference-homozygous sample with the rest
    missing (`individual_zyg`'s `hom_ref_samples` has that one entry, 723-725). GT is read at
    its FORMAT index, a bare `.` sample is missing, and columns past the kept ones are not read.
    A record Ensembl VEP reads as a structural variant (alleles not all plain bases, with a
    valued INFO `SVTYPE` or a symbolic, breakend or single-breakend allele, Parser/VCF.pm
    236-244) is annotated before the per-individual step and shows no class; `<*>` is not such
    an allele, and a bare `SVTYPE` flag or `SVTYPE=` is not such a value
    (`test_structural_record_reads_svtype_as_ensembl_vep_does`)."""
    classes = bgc.record_genotype_classes
    fixed = ["21", "100", "rs1", "A", "G", ".", "PASS", "."]
    structural = [
        ["21", "100", "sv1", "A", "<CN0>", ".", "PASS", ".", "GT", "0|1", "1|1"],
        ["21", "100", "sv2", "N", "<INS:ME:ALU>", ".", "PASS", "SVTYPE=INS", "GT", "0/1", "0/0"],
        ["21", "100", "sv3", "A", "N", ".", "PASS", "SVTYPE=DEL;END=5000", "GT", "0/1", "0/0"],
        ["21", "100", "b1", "A", "A[21:200[", ".", "PASS", ".", "GT", "0/1", "0/0"],
        ["21", "100", "b2", "A", ".G", ".", "PASS", ".", "GT", "0/1", "0/0"],
    ]
    assert [classes(f, 2) for f in structural] == [set()] * 5
    assert classes(["21", "100", "sv5", "A", "G", ".", "PASS", "SVTYPE=DEL", "GT", "0/1", "0/0"], 2) == {"het", "hom_ref"}, "plain bases: no structural variant"
    assert classes(["21", "100", "sv6", "N", "NC", ".", "PASS", "SVTYPE=", "GT", "0/1", "0/0"], 2) == {"het", "hom_ref"}, "an empty SVTYPE is none"
    assert classes(["21", "100", "sv4", "N", "NC", ".", "PASS", "SVTYPE", "GT", "0/1", "0/0"], 2) == {"het", "hom_ref"}, "a bare SVTYPE flag is none"
    assert classes(["21", "100", "g1", "A", "G,<*>", ".", "PASS", ".", "GT", "0/2", "1/2"], 2) == {"star", "hom_ref", "het"}
    assert classes([*fixed, "GT", "./.", "./.", "./."], 3) == {"missing", "all_missing"}
    assert classes([*fixed, "GT", "0/0", "./.", "./."], 3) == {"hom_ref", "missing", "no_carrier", "one_hom_ref_rest_missing"}
    assert classes([*fixed, "GT", "0/0", "0/0", "./."], 3) == {"hom_ref", "missing", "no_carrier"}
    assert classes([*fixed, "GT", "0/0", "0/1", "./."], 3) == {"hom_ref", "het", "missing"}
    assert classes([*fixed, "GT:DP", "0/0:9", "0/0:9", "0/1:9"], 2) == {"hom_ref", "no_carrier"}, "the third sample is not kept"
    assert classes([*fixed, "DP:GT", "9:1", ".:.", ".", "9:0"], 4) == {"haploid", "hom_alt", "hom_ref", "missing"}
    assert classes([*fixed, "GT", "0/0", "0/0"], 2) == {"hom_ref", "no_carrier"}
    with pytest.raises(SystemExit) as refusal:
        classes([*fixed, "DP", "9", "9"], 2)
    assert str(refusal.value) == "ERROR: [build_golden_corpus] record 21 100 rs1 A G has no GT in its FORMAT"
    with pytest.raises(SystemExit) as refusal:
        classes([*fixed, "GT", "0/0"], 2)
    assert str(refusal.value) == "ERROR: [build_golden_corpus] record 21 100 rs1 A G has 1 sample columns, 2 are kept"
    with pytest.raises(SystemExit) as refusal:
        classes(fixed, 1)
    assert str(refusal.value) == "ERROR: [build_golden_corpus] record 21 100 rs1 A G has 0 sample columns, 1 are kept"


def test_structural_record_reads_svtype_as_ensembl_vep_does():
    """`SVTYPE` makes a record structural only with a value Perl reads as true. ensembl-io's
    `get_info` (BaseVCF4.pm 576-589) splits each INFO item on `=` and stores the second piece,
    so a bare `SVTYPE` flag and `SVTYPE=` store nothing, `SVTYPE=0` stores Perl's false `0`,
    and of a repeated key the last one counts; `SVTYPE=0.0` is true. Plain `ACGT` alleles are
    never structural, and the allele test runs on the comma-joined ALT string, so a single
    breakend counts at that string's ends only."""
    structural = bgc.is_structural_record
    for info in ("SVTYPE", "END=5;SVTYPE", "SVTYPE=", "SVTYPE==", "SVTYPE=0", "SVTYPE=0=X", "SVTYPE=DEL;SVTYPE=0", "."):
        assert not structural("N", ["NC"], info), info
    for info in ("SVTYPE=DEL", "SVTYPE=DEL=X", "SVTYPE=0.0", "SVTYPE=0;SVTYPE=DEL", "END=5;SVTYPE=INS"):
        assert structural("N", ["NC"], info), info
    assert not structural("A", ["G"], "SVTYPE=DEL")
    assert structural("A", ["<CN0>"], ".") and not structural("A", ["<*>"], ".")
    assert structural("A", [".T", "G"], ".") and not structural("A", ["G", ".T"], ".")
    assert structural("A", ["G", "T."], ".") and not structural("A", ["T.", "G"], ".")


SAMPLES = [f"S{i}" for i in range(1, 13)]


def _sample_vcf(path: Path, samples: list[str] = SAMPLES) -> Path:
    """Twelve-sample records whose first eight samples show every genotype class, with a
    het carried only by the dropped samples (rs8) and GT second in one FORMAT (rs7)."""

    def rec(pos: int, vid: str, ref: str, alt: str, fmt: str, kept: str, dropped: str = "0/0") -> str:
        """One record: eight kept genotypes, space-separated, then four of `dropped`."""
        return "\t".join(["21", str(pos), vid, ref, alt, ".", "PASS", ".", fmt, *kept.split(" "), *[dropped] * 4]) + "\n"

    body = (
        rec(100, "rs1", "A", "G", "GT:DP", "0|1:9 1|1:9 0|0:9 ./.:. ./1:9 0|0:9 0|0:9 0|0:9")
        + rec(200, "rs2", "C", "T", "GT", "0/0 0/0 0/0 0/0 0/0 0/0 0/0 0/0")
        + rec(300, "rs3", "A", "G,T", "GT", "1/2 0/0 0/0 ./. ./. 0/0 0/0 0/0")
        + rec(400, "rs4", "G", "A,*", "GT", "0/2 1/2 2/2 0/0 0/0 0/0 0/0 0/0")
        + rec(500, "rs5", "T", "C", "GT", "./. ./. ./. ./. ./. ./. ./. ./.")
        + rec(600, "rs6", "G", "C", "GT", "./. 0/0 ./. ./. ./. ./. ./. ./.")
        + rec(700, "rs7", "A", "T", "DP:GT", "9:1 9:0 9:1 .:. 9:1 9:0 9:0 9:0", "9:0")
        + rec(800, "rs8", "C", "A", "GT", "0/0 0/0 0/0 0/0 0/0 0/0 0/0 0/0", "0/1")
    )
    path.write_text("##fileformat=VCFv4.2\n" + "\t".join([*bgc.FIXED_COLUMNS, "FORMAT", *samples]) + "\n" + body, encoding="utf-8")
    return path


def test_select_keeps_samples_and_every_genotype_class(tmp_path: Path, capsys):
    """With `--keep-samples 8` every record and the header keep `FORMAT` and the first eight
    samples; the reference's two exemplars show six classes and the other five each take a
    record no consequence set selected; every record lists its classes; without the option
    the corpus carries the eight fixed columns, the two exemplars and no genotype key."""
    ref = tmp_path / "ref.txt"
    ref.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant") + _row("rs2", "21:200", "T", "ENST2", "intron_variant"), encoding="utf-8")
    vcf = _sample_vcf(tmp_path / "in.vcf")
    out = tmp_path / "corpus"
    assert _select(ref, vcf, out, "--covers", "individual", "--keep-samples", "8", "--genotype-exemplars", "1") == 0
    assert capsys.readouterr().out.splitlines()[-1].endswith(", genotype classes observed 11, covered 11")
    body = (out / "variants.vcf").read_text(encoding="utf-8").splitlines()
    assert body[1] == "\t".join([*bgc.FIXED_COLUMNS, "FORMAT", *SAMPLES[:8]])
    assert [len(line.split("\t")) for line in body[2:]] == [17] * 7
    m = json.loads((out / "manifest.json").read_text())
    assert (m["keep_samples"], m["samples"], m["genotype_exemplars_per_class"], m["required_genotype_classes"]) == (8, SAMPLES[:8], 1, [])
    keys = list(m)
    assert keys.index("cache_cost_flank" if "cache_cost_flank" in m else "max_span") < keys.index("keep_samples") < keys.index("covers")
    assert keys.index("variant_classes") < keys.index("genotype_classes") < keys.index("records")
    by_id = {r["vcf"].split("\t")[2]: r for r in m["records"]}
    assert list(by_id) == ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7"], "rs8 shows nothing the others lack"
    assert by_id["rs1"]["vcf"] == "21\t100\trs1\tA\tG\t.\tPASS\t.\tGT:DP\t0|1:9\t1|1:9\t0|0:9\t./.:.\t./1:9\t0|0:9\t0|0:9\t0|0:9"
    assert by_id["rs1"]["combinations"] == ["missense_variant"]
    assert by_id["rs1"]["genotype_classes"] == ["het", "hom_alt", "hom_ref", "missing", "partially_missing"]
    assert by_id["rs2"]["genotype_classes"] == ["hom_ref", "no_carrier"]
    assert by_id["rs3"]["genotype_classes"] == ["hom_ref", "missing", "multi_allelic"]
    assert by_id["rs4"]["genotype_classes"] == ["het", "hom_ref", "star"]
    assert by_id["rs5"]["genotype_classes"] == ["all_missing", "missing"]
    assert by_id["rs6"]["genotype_classes"] == ["hom_ref", "missing", "no_carrier", "one_hom_ref_rest_missing"]
    assert by_id["rs7"]["genotype_classes"] == ["haploid", "hom_alt", "hom_ref", "missing"]
    assert all(by_id[i]["combinations"] == [] for i in ("rs3", "rs4", "rs5", "rs6", "rs7"))
    assert m["genotype_classes"] == {
        "all_missing": 1, "haploid": 1, "het": 2, "hom_alt": 2, "hom_ref": 6, "missing": 5, "multi_allelic": 1,
        "no_carrier": 2, "one_hom_ref_rest_missing": 1, "partially_missing": 1, "star": 1,
    }
    assert set(m["genotype_classes"]) == set(bgc.GENOTYPE_CLASSES)
    assert m["uncovered_genotype_classes"] == [] and keys.index("uncovered_combinations") < keys.index("uncovered_genotype_classes")
    check = ["check", "--corpus", str(out), "--suite", f"s01={ref}={vcf}"]
    assert bgc.main(check) == 0
    assert capsys.readouterr().out.splitlines()[-1] == "genotype classes observed 11, covered by corpus 11, missing 0"
    # Two records per class where the inputs have them; a class's sole record survives the cap.
    two = tmp_path / "two"
    assert _select(ref, vcf, two, "--covers", "individual", "--keep-samples", "8", "--genotype-exemplars", "2") == 0
    m2 = json.loads((two / "manifest.json").read_text())
    assert m2["genotype_exemplars_per_class"] == 2 and 7 <= len(m2["records"]) <= 8
    source_records = {"het": 2, "hom_alt": 2, "hom_ref": 7, "missing": 5, "no_carrier": 3}
    for c, n in source_records.items():
        assert m2["genotype_classes"][c] >= min(n, 2), c
    assert all(m2["genotype_classes"][c] == 1 for c in ("all_missing", "haploid", "multi_allelic", "one_hom_ref_rest_missing", "partially_missing", "star"))
    capped = tmp_path / "capped"
    assert _select(ref, vcf, capped, "--covers", "individual", "--keep-samples", "8", "--genotype-exemplars", "1", "--max-records", "3") == 0
    assert len(json.loads((capped / "manifest.json").read_text())["records"]) == 7
    plain = tmp_path / "plain"
    assert _select(ref, vcf, plain, "--covers", "individual") == 0
    assert capsys.readouterr().out.splitlines()[-1].endswith("variant classes 1 (snv 2)")
    body = (plain / "variants.vcf").read_text(encoding="utf-8").splitlines()
    assert body[1] == "\t".join(bgc.FIXED_COLUMNS) and [len(line.split("\t")) for line in body[2:]] == [8, 8]
    mp = json.loads((plain / "manifest.json").read_text())
    assert [r["vcf"].split("\t")[2] for r in mp["records"]] == ["rs1", "rs2"]
    assert not any(k in mp for k in ("keep_samples", "samples", "genotype_classes", "uncovered_genotype_classes"))
    assert "genotype_classes" not in mp["records"][0]


def test_select_with_kept_samples_is_the_same_under_every_hash_seed(tmp_path: Path):
    """The genotype-class reservoirs draw from the one seeded generator in a fixed class order,
    so interpreters with different string-hash seeds (`PYTHONHASHSEED`) write the same corpus
    from the same inputs and `--seed`; here five classes have several candidate records and no
    exemplar, so a draw order that followed set iteration would pick different records."""
    rng = random.Random(7)
    genotypes = ["0/0", "0/0", "0/0", "0/1", "1/1", "./.", "1/2", "0/2"]
    records = []
    for i in range(300):
        alt = "G,T" if i % 5 == 0 else "G"
        gts = [rng.choice(genotypes if alt == "G,T" else genotypes[:6]) for _ in range(8)]
        records.append("\t".join(["21", str(1000 + i * 10), f"v{i}", "A", alt, ".", "PASS", ".", "GT", *gts]))
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n" + "\t".join([*bgc.FIXED_COLUMNS, "FORMAT", *SAMPLES[:8]]) + "\n" + "\n".join(records) + "\n", encoding="utf-8")
    ref = tmp_path / "ref.txt"
    ref.write_text(_row("v3", "21:1030", "G", "ENST1", "missense_variant") + _row("v8", "21:1080", "G", "ENST2", "intron_variant"), encoding="utf-8")
    outputs = set()
    for hash_seed in ("0", "1", "2"):
        out = tmp_path / f"corpus-{hash_seed}"
        argv = [sys.executable, str(HERE / "build_golden_corpus.py"), "select", "--release", "116", "--assembly", "GRCh38",
                "--suite", f"s01={ref}={vcf}", "--out", str(out), "--k", "1", "--seed", "1", "--multi-allelic-per-suite", "0",
                "--covers", "individual", "--keep-samples", "8", "--genotype-exemplars", "2"]
        run = subprocess.run(argv, env={**os.environ, "PYTHONHASHSEED": hash_seed}, capture_output=True, text=True, check=False)
        assert run.returncode == 0, run.stderr
        outputs.add(((out / "variants.vcf").read_bytes(), (out / "manifest.json").read_bytes()))
    assert len(outputs) == 1
    m = json.loads((tmp_path / "corpus-0" / "manifest.json").read_text())
    assert m["uncovered_genotype_classes"] == [] and len(m["records"]) > 2, "the class fill chose records beyond the two exemplars"


def test_check_refuses_a_corpus_lacking_a_genotype_class_or_its_sample_columns(tmp_path: Path, capsys):
    """`check` recomputes the classes the inputs show among the manifest's samples and names
    each one no record shows; a required class is missing only when no record shows it; and
    `variants.vcf` must carry `FORMAT` and the manifest's samples on its header and on every
    record."""
    ref = tmp_path / "ref.txt"
    ref.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    vcf = _sample_vcf(tmp_path / "in.vcf")
    out = tmp_path / "corpus"
    assert _select(ref, vcf, out, "--covers", "individual", "--keep-samples", "8", "--genotype-exemplars", "1") == 0
    check = ["check", "--corpus", str(out), "--suite", f"s01={ref}={vcf}"]
    assert bgc.main(check) == 0
    capsys.readouterr()
    manifest = out / "manifest.json"
    m = json.loads(manifest.read_text())
    m["records"] = [r for r in m["records"] if r["vcf"].split("\t")[2] not in ("rs5", "rs7")]
    manifest.write_text(json.dumps(m), encoding="utf-8")
    missing_report = [
        "genotype classes observed 11, covered by corpus 9, missing 2",
        "  MISSING genotype class all_missing",
        "  MISSING genotype class haploid",
    ]
    assert bgc.main(check) == 1
    assert capsys.readouterr().out.splitlines()[-3:] == missing_report
    assert bgc.main([*check, "--require-genotype-class", "haploid", "--require-genotype-class", "star"]) == 1
    assert capsys.readouterr().out.splitlines()[-3:] == missing_report, "a required class the corpus shows is not missing"
    with pytest.raises(SystemExit) as refusal:
        bgc.main([*check, "--require-genotype-class", "hom-ref"])
    assert str(refusal.value) == f"ERROR: [build_golden_corpus] --require-genotype-class wants one of {', '.join(bgc.GENOTYPE_CLASSES)}, got 'hom-ref'"
    variants = out / "variants.vcf"
    lines = variants.read_text(encoding="utf-8").splitlines()
    lines[1] = "\t".join(lines[1].split("\t")[:16])
    lines[2] = "\t".join(lines[2].split("\t")[:9])
    variants.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert bgc.main(check) == 1
    captured = capsys.readouterr()
    assert captured.err == (
        "ERROR: [build_golden_corpus] variants.vcf header names FORMAT S1 S2 S3 S4 S5 S6 S7, the manifest keeps FORMAT S1 S2 S3 S4 S5 S6 S7 S8\n"
        "ERROR: [build_golden_corpus] variants.vcf record 21 100 rs1 A G has 9 columns, the header 17\n"
    )
    variants.write_text("##fileformat=VCFv4.2\n" + "\n".join(lines[2:]) + "\n", encoding="utf-8")
    assert bgc.main(check) == 1
    assert capsys.readouterr().err == (
        "ERROR: [build_golden_corpus] variants.vcf record 21 100 rs1 A G has 9 columns, the header 17\n"
        "ERROR: [build_golden_corpus] variants.vcf has no #CHROM header line\n"
    )
    # An input that keeps other samples than the manifest's cannot be the corpus's source.
    other = _sample_vcf(tmp_path / "other.vcf", [*SAMPLES[:7], "X8", *SAMPLES[8:]])
    with pytest.raises(SystemExit) as refusal:
        bgc.main(["check", "--corpus", str(out), "--suite", f"s01={ref}={other}"])
    assert str(refusal.value) == (
        f"ERROR: [build_golden_corpus] {other} keeps samples S1,S2,S3,S4,S5,S6,S7,X8 where the manifest keeps "
        "S1,S2,S3,S4,S5,S6,S7,S8: every input must keep the same samples"
    )
    # A plain corpus has no sample columns to check and refuses a required class.
    plain = tmp_path / "plain"
    assert _select(ref, vcf, plain, "--covers", "individual") == 0
    with pytest.raises(SystemExit) as refusal:
        bgc.main(["check", "--corpus", str(plain), "--suite", f"s01={ref}={vcf}", "--require-genotype-class", "het"])
    assert str(refusal.value) == "ERROR: [build_golden_corpus] --require-genotype-class wants a corpus built with --keep-samples"


def test_check_names_an_unreadable_variants_vcf(tmp_path: Path, capsys):
    """A VCF the tool cannot read, the corpus's `variants.vcf` or a suite input, is one error
    naming the file, never a traceback, under `check` and under `select`; the errors `check`
    found before the refusal are printed ahead of it. A `.gz` input that is not gzip, or is
    cut short, fails on its first line rather than at the open and is refused the same way."""
    ref = tmp_path / "ref.txt"
    ref.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    vcf = _sample_vcf(tmp_path / "in.vcf")
    out = tmp_path / "corpus"
    assert _select(ref, vcf, out, "--covers", "individual", "--keep-samples", "8", "--genotype-exemplars", "1") == 0
    capsys.readouterr()
    variants = out / "variants.vcf"
    variants.unlink()
    check = ["check", "--corpus", str(out), "--suite", f"s01={ref}={vcf}"]
    for shape in ("missing", "directory"):
        if shape == "directory":
            variants.mkdir()
        with pytest.raises(SystemExit) as refusal:
            bgc.main(check)
        assert str(refusal.value).startswith(f"ERROR: [build_golden_corpus] {variants} is not a readable text file ("), shape
    variants.rmdir()
    variants.write_text(vcf.read_text(encoding="utf-8").replace("\tS8\t", "\tX8\t"), encoding="utf-8")
    with pytest.raises(SystemExit) as refusal:
        bgc.main(["check", "--corpus", str(out), "--suite", f"s01={ref}={tmp_path / 'absent.vcf'}"])
    assert str(refusal.value).startswith(f"ERROR: [build_golden_corpus] {tmp_path / 'absent.vcf'} is not a readable text file (")
    assert capsys.readouterr().err == (
        "ERROR: [build_golden_corpus] variants.vcf header names FORMAT S1 S2 S3 S4 S5 S6 S7 X8 S9 S10 S11 S12, "
        "the manifest keeps FORMAT S1 S2 S3 S4 S5 S6 S7 S8\n"
        + "".join(f"ERROR: [build_golden_corpus] variants.vcf record 21 {pos} rs{i} {ref_alt} has 21 columns, the header 17\n"
                  for i, (pos, ref_alt) in enumerate([(100, "A G"), (200, "C T"), (300, "A G,T"), (400, "G A,*"), (500, "T C"), (600, "G C"), (700, "A T"), (800, "C A")], start=1))
    ), "the column errors come before the refusal"
    with pytest.raises(SystemExit) as refusal:
        _select(ref, tmp_path / "absent.vcf", tmp_path / "none", "--covers", "individual")
    assert str(refusal.value).startswith(f"ERROR: [build_golden_corpus] {tmp_path / 'absent.vcf'} is not a readable text file (")
    gzipped = gzip.compress(vcf.read_bytes())
    for name, blob in (("plain.vcf.gz", vcf.read_bytes()), ("cut.vcf.gz", gzipped[: len(gzipped) // 2])):
        bad = tmp_path / name
        bad.write_bytes(blob)
        with pytest.raises(SystemExit) as refusal:
            _select(ref, bad, tmp_path / "none", "--covers", "individual", "--keep-samples", "8")
        assert str(refusal.value).startswith(f"ERROR: [build_golden_corpus] {bad} is not a readable text file ("), name
        with pytest.raises(SystemExit) as refusal:
            bgc.main(["check", "--corpus", str(out), "--suite", f"s01={ref}={bad}"])
        assert str(refusal.value).startswith(f"ERROR: [build_golden_corpus] {bad} is not a readable text file ("), name
        capsys.readouterr()
    assert not (tmp_path / "none").exists()


def test_select_requires_named_genotype_classes_and_one_sample_set(tmp_path: Path, capsys):
    """A required class the inputs never show leaves `select` red with the class listed, as a
    required class does on `check`; a second input keeping other samples, an input with too few,
    an unknown class, a negative count and a requirement without kept samples are refused
    before anything is written."""
    ref = tmp_path / "ref.txt"
    ref.write_text(_row("rs1", "21:100", "G", "ENST1", "missense_variant"), encoding="utf-8")
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        "21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n21\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t1/1\t./.\n",
        encoding="utf-8",
    )
    out = tmp_path / "corpus"
    assert _select(ref, vcf, out, "--covers", "individual", "--keep-samples", "2", "--require-genotype-class", "haploid",
                   "--require-genotype-class", "het", "--require-genotype-class", "haploid") == 1
    assert capsys.readouterr().out.splitlines()[-1].endswith(", genotype classes observed 4, covered 4")
    m = json.loads((out / "manifest.json").read_text())
    assert m["required_genotype_classes"] == ["haploid", "het"] and m["uncovered_genotype_classes"] == ["haploid"]
    assert m["genotype_classes"] == {"het": 1, "hom_alt": 1, "hom_ref": 1, "missing": 1}
    check = ["check", "--corpus", str(out), "--suite", f"s01={ref}={vcf}"]
    assert bgc.main(check) == 1
    assert capsys.readouterr().out.splitlines()[-2:] == ["genotype classes observed 4, covered by corpus 4, missing 1", "  MISSING genotype class haploid"]
    m["required_genotype_classes"] = []
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    assert bgc.main(check) == 0
    before = {p.name: p.read_text() for p in out.iterdir()}
    other = tmp_path / "other.vcf"
    other.write_text(vcf.read_text().replace("\tA\tB\n", "\tA\tC\n"), encoding="utf-8")
    refusals = [
        (["--suite", f"s02={ref}={other}", "--keep-samples", "2"],
         f"{other} keeps samples A,C where an earlier input keeps A,B: every input must keep the same samples"),
        (["--keep-samples", "3"], f"{vcf} has 2 sample columns, 3 are kept"),
        (["--keep-samples", "-1"], "--keep-samples wants a count of sample columns, got -1"),
        (["--keep-samples", "2", "--genotype-exemplars", "-1"], "--genotype-exemplars wants a count of records per class, got -1"),
        (["--keep-samples", "2", "--require-genotype-class", "hom-ref"],
         f"--require-genotype-class wants one of {', '.join(bgc.GENOTYPE_CLASSES)}, got 'hom-ref'"),
        (["--require-genotype-class", "het"], "--require-genotype-class wants --keep-samples: a corpus without sample columns shows no genotype"),
    ]
    headerless = tmp_path / "headerless.vcf"
    headerless.write_text("21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n", encoding="utf-8")
    no_format = tmp_path / "no_format.vcf"
    no_format.write_text(vcf.read_text(encoding="utf-8").replace("\tFORMAT\tA\tB\n", "\tA\tB\tC\n"), encoding="utf-8")
    # A header-less input beside a headed one, and a header without FORMAT, bind no sample to
    # a name: each is refused by name, in either position.
    refusals += [
        (["--suite", f"s02={ref}={headerless}", "--keep-samples", "2"], f"{headerless} has no #CHROM header line naming its samples"),
        (["--suite", f"s02={ref}={no_format}", "--keep-samples", "2"], f"{no_format} has no FORMAT column, so no sample columns"),
    ]
    for extra, message in refusals:
        with pytest.raises(SystemExit) as refusal:
            _select(ref, vcf, out, "--covers", "individual", *extra)
        assert str(refusal.value) == f"ERROR: [build_golden_corpus] {message}"
    assert {p.name: p.read_text() for p in out.iterdir()} == before
    with pytest.raises(SystemExit) as refusal:
        _select(ref, headerless, tmp_path / "none", "--covers", "individual", "--keep-samples", "2", "--suite", f"s02={ref}={vcf}")
    assert str(refusal.value) == f"ERROR: [build_golden_corpus] {headerless} has no #CHROM header line naming its samples"
    empty = tmp_path / "empty.vcf"
    empty.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    with pytest.raises(SystemExit) as refusal:
        _select(ref, empty, tmp_path / "none", "--covers", "individual", "--keep-samples", "2")
    assert str(refusal.value) == "ERROR: [build_golden_corpus] no input has a #CHROM header line naming its samples"
    assert not (tmp_path / "none").exists()


def test_flag_set_subcommand_writes_an_entry(tmp_path: Path, capsys):
    """An entry carries the whole flag list of its run, its input, stdin and unordered modes,
    the formats compared, its expected directory, its covers and, for a reference generated
    against the Ensembl database, its `perl_mode`; a set is rewritten by name; the inputs the
    harness cannot run are refused before anything is written."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    matrix = _matrix(tmp_path, "symbol", "pick", "per_gene")
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    capsys.readouterr()
    assert bgc.main(["flag-set", str(out), "--name", "pick", "--flag=--pick", "--flag=--symbol", "--covers", "pick", "--matrix", str(matrix)]) == 0
    assert capsys.readouterr().out == "flag set pick: 2 flags, input variants.vcf, formats default,tab,vcf,json, expected expected_pick, covers 1\n"
    m = json.loads((out / "manifest.json").read_text())
    assert m["flag_sets"] == [{
        "name": "pick", "flags": ["--pick", "--symbol"], "input": "variants.vcf", "stdin": False, "unordered": False,
        "formats": ["default", "tab", "vcf", "json"], "expected_dir": "expected_pick", "covers": ["pick"],
    }]
    (out / "variants.ensembl").write_text("21\t100\t100\tA/G\t+\n", encoding="utf-8")
    assert bgc.main(["flag-set", str(out), "--name", "ensembl-stdin", "--flag=--format", "--flag=ensembl", "--input", "variants.ensembl",
                     "--stdin", "--format", "default,json", "--expected-dir", "expected_ens", "--matrix", str(matrix)]) == 0
    # A run whose reference orders a record's rows by Ensembl VEP's hash order is marked unordered.
    assert bgc.main(["flag-set", str(out), "--name", "pick", "--flag=--per_gene", "--unordered", "--covers", "per_gene",
                     "--matrix", str(matrix)]) == 0
    assert capsys.readouterr().out.endswith(
        "flag set pick: 1 flags, input variants.vcf, unordered, formats default,tab,vcf,json, expected expected_pick, covers 1\n")
    m = json.loads((out / "manifest.json").read_text())
    assert [fs["name"] for fs in m["flag_sets"]] == ["pick", "ensembl-stdin"]
    assert m["flag_sets"][0]["flags"] == ["--per_gene"] and m["flag_sets"][0]["covers"] == ["per_gene"]
    assert m["flag_sets"][0]["unordered"] is True
    assert m["flag_sets"][1] == {
        "name": "ensembl-stdin", "flags": ["--format", "ensembl"], "input": "variants.ensembl", "stdin": True,
        "unordered": False, "formats": ["default", "json"], "expected_dir": "expected_ens", "covers": [],
    }
    # `--format=<format>` is the other spelling of the format a stdin run must name.
    assert bgc.main(["flag-set", str(out), "--name", "eq", "--flag=--format=vcf", "--stdin", "--matrix", str(matrix)]) == 0
    capsys.readouterr()
    # A reference generated against the Ensembl database is marked, and only then.
    assert bgc.main(["flag-set", str(out), "--name", "db", "--flag=--ga4gh_vrs", "--format", "default,tab,json",
                     "--perl-mode", "database", "--matrix", str(matrix)]) == 0
    assert capsys.readouterr().out == ("flag set db: 1 flags, input variants.vcf, reference database, "
                                       "formats default,tab,json, expected expected_db, covers 0\n")
    m = json.loads((out / "manifest.json").read_text())
    assert m["flag_sets"][3]["perl_mode"] == "database"
    assert "perl_mode" not in m["flag_sets"][2]
    assert bgc.manifest_shape_error(m) is None
    m["flag_sets"][3]["perl_mode"] = "cache"
    assert bgc.manifest_shape_error(m) == "flag set db perl_mode is not one of offline, database, gtf, gff"
    # A reference Ensembl VEP wrote reading its transcripts from a GTF or GFF is marked the same
    # way, and the manifest's own run may carry the mode, which `select` keeps on a regeneration
    # and refuses when it is not one it knows.
    assert bgc.main(["flag-set", str(out), "--name", "gtf", "--flag=--symbol", "--perl-mode", "gtf", "--matrix", str(matrix)]) == 0
    assert capsys.readouterr().out == ("flag set gtf: 1 flags, input variants.vcf, reference gtf, "
                                       "formats default,tab,vcf,json, expected expected_gtf, covers 0\n")
    m = json.loads((out / "manifest.json").read_text())
    assert m["flag_sets"][4]["perl_mode"] == "gtf"
    m["perl_mode"] = "gff"
    assert bgc.manifest_shape_error(m) is None
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    capsys.readouterr()
    m = json.loads((out / "manifest.json").read_text())
    assert m["perl_mode"] == "gff"
    m["perl_mode"] = "cache"
    assert bgc.manifest_shape_error(m) == "perl_mode is not one of offline, database, gtf, gff"
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    with pytest.raises(SystemExit, match="perl_mode is not one of"):
        bgc.main(["flag-set", str(out), "--list"])
    del m["perl_mode"]
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    assert bgc.main(["flag-set", str(out), "--list"]) == 0
    assert capsys.readouterr().out == (
        "pick\t--per_gene\tvariants.vcf\tfalse\ttrue\tdefault,tab,vcf,json\texpected_pick\tper_gene\n"
        "ensembl-stdin\t--format ensembl\tvariants.ensembl\ttrue\tfalse\tdefault,json\texpected_ens\t\n"
        "eq\t--format=vcf\tvariants.vcf\ttrue\tfalse\tdefault,tab,vcf,json\texpected_eq\t\n"
        "db\t--ga4gh_vrs\tvariants.vcf\tfalse\tfalse\tdefault,tab,json\texpected_db\t\n"
        "gtf\t--symbol\tvariants.vcf\tfalse\tfalse\tdefault,tab,vcf,json\texpected_gtf\t\n"
    )
    before = (out / "manifest.json").read_text()
    refused = [
        (["--name", "default"], "flag set name 'default': letters, digits, _ and -, and not the manifest's own run, default"),
        (["--name", "a/b"], "flag set name 'a/b': letters, digits, _ and -, and not the manifest's own run, default"),
        (["--name", "x", "--format", "default,parquet"], "--format wants a non-empty subset of default,tab,vcf,json, got 'default,parquet'"),
        # A set comparing no format would certify its covers without a comparison.
        (["--name", "x", "--format", ","], "--format wants a non-empty subset of default,tab,vcf,json, got ','"),
        (["--name", "x", "--input", "missing.vcf"], f"input missing.vcf is not a file in {out}"),
        (["--name", "x", "--covers", "ghost"], f"flag set x covers names ghost, which is not a row of {matrix}"),
        # Ensembl VEP detects no format on STDIN, so a stdin run without --format has no reference.
        (["--name", "x", "--stdin"], "a stdin run names its input format: add --flag=--format --flag=<format>"),
        ([], "flag-set wants --name (or --list)"),
    ]
    for extra, message in refused:
        assert bgc.main(["flag-set", str(out), *extra, "--matrix", str(matrix)]) == 1
        assert capsys.readouterr().err == f"ERROR: [build_golden_corpus] {message}\n"
    # A path that leaves the corpus directory, or an empty covers name, is refused as the tool exits.
    outside = str(tmp_path / "outside.vcf")
    Path(outside).write_text("", encoding="utf-8")
    for extra, message in [
        (["--input", outside], f"--input wants a path inside the corpus directory, got {outside!r}"),
        (["--input", "../outside.vcf"], "--input wants a path inside the corpus directory, got '../outside.vcf'"),
        (["--expected-dir", "../expected_x"], "--expected-dir wants a path inside the corpus directory, got '../expected_x'"),
        (["--covers", ""], "--covers wants a non-empty row name"),
    ]:
        with pytest.raises(SystemExit) as refusal:
            bgc.main(["flag-set", str(out), "--name", "x", *extra, "--matrix", str(matrix)])
        assert str(refusal.value) == f"ERROR: [build_golden_corpus] {message}"
    assert (out / "manifest.json").read_text() == before


def test_flag_set_rewrite_keeps_the_sets_hand_written_keys(tmp_path: Path, capsys):
    """Rewriting a set by name sets the options' keys where they sit and keeps every other key of
    the entry (`masks` and `reference_warnings` are written by hand and have no option); the
    divergence keys describe the run before the rewrite and go, and `perl_mode` follows the
    option like the rest."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    matrix = _matrix(tmp_path, "symbol", "domains", "per_gene")
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    assert bgc.main(["flag-set", str(out), "--name", "domains", "--flag=--domains", "--covers", "domains",
                     "--perl-mode", "database", "--matrix", str(matrix)]) == 0
    m = json.loads((out / "manifest.json").read_text())
    entry = m["flag_sets"][0]
    notice = "found tandem repeats with different references per alternative allele: SVLEN=9,6; only using reference with largest size"
    # The hand-written keys sit before `covers`, as a hand-edited manifest places them.
    m["flag_sets"][0] = {
        **{k: v for k, v in entry.items() if k != "covers"},
        "masks": ["order_only"], "reference_warnings": [notice], "covers": entry["covers"],
        "divergences": [{"location": "21:100", "record_indices": [0]}], "divergence_summary": {},
        "note": "kept as written",
    }
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    capsys.readouterr()
    assert bgc.main(["flag-set", str(out), "--name", "domains", "--flag=--domains", "--flag=--per_gene", "--unordered",
                     "--covers", "domains", "--covers", "per_gene", "--matrix", str(matrix)]) == 0
    assert capsys.readouterr().out == ("flag set domains: 2 flags, input variants.vcf, unordered, "
                                       "formats default,tab,vcf,json, expected expected_domains, covers 2\n")
    m = json.loads((out / "manifest.json").read_text())
    assert m["flag_sets"] == [{
        "name": "domains", "flags": ["--domains", "--per_gene"], "input": "variants.vcf", "stdin": False,
        "unordered": True, "formats": ["default", "tab", "vcf", "json"], "expected_dir": "expected_domains",
        "masks": ["order_only"], "reference_warnings": [notice], "covers": ["domains", "per_gene"],
        "note": "kept as written",
    }]
    assert bgc.manifest_shape_error(m) is None
    # A kept key is still held to its shape at the next read, so a rewrite never writes
    # past a `masks` or `reference_warnings` the harness would refuse.
    for body, reason in (
        ({"masks": ["splice_rows"]},
         f"flag set domains masks names unknown classes ['splice_rows'] (the classes are {', '.join(bgc.MASK_CLASSES)})"),
        ({"masks": "order_only"}, "flag set domains masks is not a list of names"),
        ({"reference_warnings": [notice, ""]}, "flag set domains reference_warnings is not a list of texts"),
        ({"reference_warnings": notice}, "flag set domains reference_warnings is not a list of texts"),
    ):
        broken = json.loads((out / "manifest.json").read_text())
        broken["flag_sets"][0].update(body)
        (out / "manifest.json").write_text(json.dumps(broken), encoding="utf-8")
        with pytest.raises(SystemExit) as refusal:
            bgc.main(["flag-set", str(out), "--name", "domains", "--flag=--domains", "--matrix", str(matrix)])
        assert str(refusal.value) == f"ERROR: [build_golden_corpus] {out / 'manifest.json'} {reason}"
        assert json.loads((out / "manifest.json").read_text()) == broken, "nothing is written"


def test_classify_applies_the_runs_own_masks(tmp_path: Path):
    """`classify` masks a run as the golden harness does: under `order_only` the `DOMAINS` and
    `FREQS` values of `--field` agree as multisets of their items while an item that differs
    is still a field divergence, and a flag set's own `masks` (an empty list included) replace
    the manifest's for `classify --flag-set`, the regulatory class among them."""

    def row(domains: str) -> str:
        return _row("rs1", "21:100", "G", "ENST1", "missense_variant", f"IMPACT=MODERATE;DOMAINS={domains}")

    regulatory = "rs1\t21:100\tG\t-\tENSR1\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
    vep = tmp_path / "vep.txt"
    vep.write_text(row("Pfam:PF00059,PANTHER:PTHR14789") + regulatory, encoding="utf-8")
    reordered = tmp_path / "reordered.txt"
    reordered.write_text(row("PANTHER:PTHR14789,Pfam:PF00059"), encoding="utf-8")
    changed = tmp_path / "changed.txt"
    changed.write_text(row("PANTHER:PTHR14789,Pfam:PF00060"), encoding="utf-8")
    corpus = tmp_path / "c"
    corpus.mkdir()
    manifest = corpus / "manifest.json"

    def classify(rs: Path, flag_set: str | None = None) -> dict:
        flag_set_option = [] if flag_set is None else ["--flag-set", flag_set]
        assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                         "--field", "DOMAINS", *flag_set_option]) == 0
        m = json.loads(manifest.read_text())
        target = m if flag_set is None else bgc.flag_set_named(m, flag_set)
        return {"fields": target["field_divergence_summary"], "vep_only": target["divergence_summary"]}

    records = [_record("21\t100\trs1\tA\tG\t.\tPASS\t.")]
    sets = [{"name": "plain", "masks": []}, {"name": "regulatory", "masks": ["regulatory_feature_rows"]}, {"name": "inherits"}]
    manifest.write_text(json.dumps({"masks": ["order_only"], "records": records, "flag_sets": sets}), encoding="utf-8")
    assert classify(reordered) == {"fields": {}, "vep_only": {"unexplained_residual": 1}}
    assert classify(changed) == {"fields": {"unexplained_residual": 1}, "vep_only": {"unexplained_residual": 1}}
    assert classify(reordered, "inherits") == {"fields": {}, "vep_only": {"unexplained_residual": 1}}
    assert classify(reordered, "plain") == {"fields": {"unexplained_residual": 1}, "vep_only": {"unexplained_residual": 1}}
    assert classify(reordered, "regulatory") == {"fields": {"unexplained_residual": 1}, "vep_only": {}}
    manifest.write_text(json.dumps({"records": records}), encoding="utf-8")
    assert classify(reordered) == {"fields": {"unexplained_residual": 1}, "vep_only": {"unexplained_residual": 1}}


def test_order_only_compares_the_zyg_list_as_a_multiset():
    """Under `order_only` the `ZYG` list of `--individual_zyg` agrees in any order of its
    `sample:ZYG` items (Ensembl VEP prints it in hash order), while an item that differs,
    or the same list without the mask, is a field divergence."""
    assert bgc.fields_agree("ZYG", "S7:HOMREF,S4:HET,S1:HOM", "S1:HOM,S4:HET,S7:HOMREF", True)
    assert not bgc.fields_agree("ZYG", "S7:HOMREF,S4:HET,S1:HOM", "S1:HOM,S4:HETREF,S7:HOMREF", True)
    assert not bgc.fields_agree("ZYG", "S7:HOMREF,S4:HET,S1:HOM", "S1:HOM,S4:HET,S7:HOMREF", False)
    assert "ZYG" in bgc.ORDER_ONLY_COLUMNS


def test_classify_writes_a_flag_sets_divergences_into_its_entry(tmp_path: Path, capsys):
    """`classify --flag-set` documents that set's run: its divergence keys land in the entry, the
    manifest's own stay as they were, and the default run's `--flag` option is refused."""
    corpus = tmp_path / "c"
    corpus.mkdir()
    (corpus / "manifest.json").write_text(json.dumps({
        "records": [{"suite": "s01", "vcf": "21\t100\trs1\tA\tG\t.\tPASS\t.", "combinations": []}],
        "covers": ["symbol"],
        "divergences": [{"location": "kept"}],
        "flag_sets": [{"name": "pick", "flags": ["--pick"], "covers": ["pick"]}],
    }), encoding="utf-8")
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(_row("rs1", "21:100", "G", "ENST1", "start_lost,start_retained_variant"), encoding="utf-8")
    rs.write_text(_row("rs1", "21:100", "G", "ENST1", "start_retained_variant"), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs), "--flag-set", "pick", "--flag=--pick"]) == 1
    assert capsys.readouterr().err == (
        "ERROR: [build_golden_corpus] --flag and --fasta-name describe the default run; a flag set's flags are its own entry's\n"
    )
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs), "--flag-set", "pick"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert m["divergences"] == [{"location": "kept"}] and m["covers"] == ["symbol"]
    fs = m["flag_sets"][0]
    assert fs["covers"] == ["pick"] and fs["divergence_summary"] == {"start_cooccurrence_swap": 1}
    assert [d["location"] for d in fs["divergences"]] == ["21:100"] and fs["vep_rs_only_tuples"] == []
    with pytest.raises(SystemExit):
        bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs), "--flag-set", "none"])


def test_classify_records_field_divergences_on_agreeing_keys(tmp_path: Path):
    """A named field that differs on a key whose consequence sets agree is a
    field divergence (never documented, so `unexplained_residual`); a key whose sets
    differ is a consequence divergence and its fields are not compared."""
    corpus = tmp_path / "c"
    corpus.mkdir()
    (corpus / "manifest.json").write_text(json.dumps({"records": []}), encoding="utf-8")
    vep = tmp_path / "vep.txt"
    rs = tmp_path / "rs.txt"
    vep.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant", "IMPACT=MODERATE;STRAND=1;HGVSc=ENST1.1:c.5C>A;HGVSp=ENSP1.1:p.Ala2Asp")
        + _row("rs2", "21:200", "-", "ENST2", "inframe_insertion", "IMPACT=MODERATE;STRAND=1;HGVSc=ENST2.1:c.9_10insAGA;HGVSp=ENSP2.1:p.Gly3_Lys4insArg")
        + _row("rs3", "21:300", "-", "ENST3", "start_lost,start_retained_variant", "IMPACT=HIGH;STRAND=1;HGVSp=ENSP3.1:p.Met1?"),
        encoding="utf-8",
    )
    rs.write_text(
        _row("rs1", "21:100", "G", "ENST1", "missense_variant", "IMPACT=MODERATE;STRAND=1;HGVSc=ENST1.1:c.5C>A;HGVSp=ENSP1.1:p.Ala2Asp")
        + _row("rs2", "21:200", "-", "ENST2", "inframe_insertion", "IMPACT=MODERATE;STRAND=1;HGVSc=ENST2.1:c.9_10insAGA;HGVSp=-")
        + _row("rs3", "21:300", "-", "ENST3", "start_retained_variant", "IMPACT=LOW;STRAND=1;HGVSp=ENSP3.1:p.Met1del"),
        encoding="utf-8",
    )
    rc = bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs), "--field", "HGVSc", "--field", "HGVSp",
                   "--flag=--hgvs", "--fasta-name", "reference.fa.gz"])
    assert rc == 0
    m = json.loads((corpus / "manifest.json").read_text())
    assert m["fields_compared"] == ["HGVSc", "HGVSp"]
    # The run options the harness replays on both engines travel with the manifest.
    assert m["flags"] == ["--hgvs"] and m["fasta"] == "reference.fa.gz"
    assert [(d["location"], d["field"], d["vep"], d["vep_rs"], d["expected_divergence"]) for d in m["field_divergences"]] == [
        ("21:200", "HGVSp", "ENSP2.1:p.Gly3_Lys4insArg", "-", "unexplained_residual"),
    ]
    assert m["field_divergence_summary"] == {"unexplained_residual": 1}
    assert {d["location"]: d["expected_divergence"] for d in m["divergences"]} == {"21:300": "start_cooccurrence_swap"}


def test_extract_reference_writes_the_contig_and_its_index(tmp_path: Path):
    src = tmp_path / "genome.fa"
    src.write_text(">1 dna:chromosome\nACGTACGTAC\nGGGG\n>21 dna:chromosome chromosome:GRCh37:21:1:14:1 REF\nAACCGGTTAA\nCCGN\n>22\nTTTT\n", encoding="utf-8")
    out = tmp_path / "reference.fa"
    rc = xr.main([str(src), "--contig", "21", "--out", str(out)])
    assert rc == 0
    assert out.read_text(encoding="utf-8") == ">21\nAACCGGTTAA\nCCGN\n"
    assert Path(str(out) + ".fai").read_text(encoding="utf-8") == "21\t14\t4\t10\t11\n"
    with pytest.raises(SystemExit):
        xr.main([str(src), "--contig", "MT", "--out", str(tmp_path / "none.fa")])


def test_prune_keeps_reachable_transcripts_in_every_shard(tmp_path: Path):
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    t_a = {"stable_id": "ENST_A", "start": "10000", "end": "20000", "dbID": "1"}
    t_b = {"stable_id": "ENST_B", "start": "990000", "end": "1010000", "dbID": "2"}  # straddles the 1M boundary
    t_c = {"stable_id": "ENST_C", "start": "500000", "end": "500500", "dbID": "3"}
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([t_a, t_b, t_c]))
    (cache / "transcripts" / "21" / "1000001-2000000.json").write_text(json.dumps([t_b]))
    (cache / "transcripts" / "22").mkdir()
    (cache / "transcripts" / "22" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_D", "start": "5", "end": "9", "dbID": "4"}]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t24000\t.\tA\tG\t.\t.\t.\n21\t995000\t.\tT\t<DEL>\t.\t.\tEND=996000\n")
    out = tmp_path / "out"
    rc = pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--assembly", "GRCh37", "--cache-version", "115"])
    assert rc == 0
    s1 = json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())
    s2 = json.loads((out / "transcripts" / "21" / "1000001-2000000.json").read_text())
    assert [t["stable_id"] for t in s1] == ["ENST_A", "ENST_B"]  # A within 5 kb flank, C unreachable
    assert [t["stable_id"] for t in s2] == ["ENST_B"]  # the straddling copy survives in both shards
    assert not (out / "transcripts" / "22").exists()
    info = json.loads((out / "info.json").read_text())
    assert info["assembly"] == "GRCh37" and info["cache_version"] == 115
    for name in ("1-1000000.json.gz", "1-1000000.json"):
        assert name.removesuffix(".gz") == (name[: -len(".gz")] if name.endswith(".gz") else name) == "1-1000000.json"  # noqa: FURB188


def test_prune_keeps_the_same_symbol_hgnc_donors_of_an_unnamed_model(tmp_path: Path, capsys):
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    reached = {"stable_id": "NR_A", "start": "10000", "end": "20000", "dbID": "1", "source_cache": "RefSeq",
               "gene_symbol": "ANKRD30BP2", "gene_symbol_source": "EntrezGene"}
    donor = {"stable_id": "ENST_DONOR", "start": "500000", "end": "500500", "dbID": "2", "source_cache": "Ensembl",
             "gene_symbol": "ANKRD30BP2", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:16620"}
    same_symbol_unnamed = {"stable_id": "ENST_UNNAMED", "start": "600000", "end": "600500", "dbID": "3",
                           "source_cache": "Ensembl", "gene_symbol": "ANKRD30BP2"}
    other_symbol = {"stable_id": "ENST_OTHER", "start": "700000", "end": "700500", "dbID": "4", "source_cache": "Ensembl",
                    "gene_symbol": "MRPS6", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:14051"}
    named_reached = {"stable_id": "ENST_B", "start": "1200000", "end": "1210000", "dbID": "5", "source_cache": "Ensembl",
                     "gene_symbol": "MRPS6", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:14051"}
    unnamed_reached = {"stable_id": "NM_B", "start": "1200000", "end": "1210000", "dbID": "6", "source_cache": "RefSeq",
                       "gene_symbol": "MRPS6", "gene_symbol_source": "EntrezGene"}
    spare_donor = {"stable_id": "ENST_SPARE", "start": "1500000", "end": "1500500", "dbID": "7", "source_cache": "Ensembl",
                   "gene_symbol": "MRPS6", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:14051"}
    far_donor = {"stable_id": "ENST_FAR", "start": "1600000", "end": "1600500", "dbID": "8", "source_cache": "Ensembl",
                 "gene_symbol": "ANKRD30BP2", "gene_symbol_source": "HGNC", "gene_hgnc_id": "HGNC:16620"}
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([reached, donor, same_symbol_unnamed, other_symbol]))
    (cache / "transcripts" / "21" / "1000001-2000000.json").write_text(json.dumps([named_reached, unnamed_reached, spare_donor, far_donor]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t24000\t.\tA\tG\t.\t.\t.\n21\t1205000\t.\tT\tC\t.\t.\t.\n")
    out = tmp_path / "out"
    rc = pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--assembly", "GRCh38", "--cache-version", "116"])
    assert rc == 0
    s1 = json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())
    s2 = json.loads((out / "transcripts" / "21" / "1000001-2000000.json").read_text())
    # The donor stays in the shard's order; a same-symbol model without an id and a model of
    # another symbol are not donors; a symbol a reached model already carries with its id pulls
    # nothing in, and a donor is read from the unnamed model's own shard alone.
    assert [t["stable_id"] for t in s1] == ["NR_A", "ENST_DONOR"]
    assert [t["stable_id"] for t in s2] == ["ENST_B", "NM_B"]
    assert json.loads(capsys.readouterr().out)["hgnc_donors"] == 1


def test_prune_regulatory_keeps_the_reachable_features_and_declares_them(tmp_path: Path, capsys):
    """`--regulatory` writes the regulatory and motif features within the flank of a record
    span, in shards named as in the source, the binding matrices the kept motifs name, and
    `info.json` declares `regulatory` with the schema and cell-type list only when a shard was
    written; without the flag the tree is left out and the flag reads false."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh38", "cache_version": 116,
                                                 "regulatory": True, "regulatory_schema": 1, "cell_types": ["A549", "HUVEC"]}))
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "10000", "end": "20000", "dbID": "1"}]))
    (cache / "regulatory" / "21").mkdir(parents=True)
    near = {"stable_id": "ENSR_NEAR", "start": 26000, "end": 27000, "feature_type": "RegulatoryFeature", "biotype": "enhancer", "activity": "AI"}
    far = {"stable_id": "ENSR_FAR", "start": 500000, "end": 500500, "feature_type": "RegulatoryFeature", "biotype": "promoter", "activity": "IA"}
    motif = {"stable_id": "ENSM_IN", "start": 23990, "end": 24010, "strand": 1, "feature_type": "MotifFeature", "matrix": "ENSPFM1"}
    far_motif = {"stable_id": "ENSM_FAR", "start": 600000, "end": 600010, "strand": 1, "feature_type": "MotifFeature", "matrix": "ENSPFM2"}
    (cache / "regulatory" / "21" / "1-1000000.json").write_text(json.dumps([near, far, motif, far_motif]))
    (cache / "regulatory" / "binding_matrices.json").write_text(json.dumps(
        {"ENSPFM1": {"transcription_factors": ["TF"], "frequencies": []}, "ENSPFM2": {"transcription_factors": ["TG"], "frequencies": []}}))
    with gzip.open(cache / "regulatory" / "21" / "1000001-2000000.json.gz", "wt") as fh:
        fh.write(json.dumps([{"stable_id": "ENSR_OTHER", "start": 1500000, "end": 1500100, "strand": 0, "feature_type": "RegulatoryFeature"}]))
    (cache / "regulatory" / "22").mkdir()
    (cache / "regulatory" / "22" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENSR_22", "start": 5, "end": 9, "strand": 0, "feature_type": "RegulatoryFeature"}]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t24000\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    rc = pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--regulatory", "--gzip"])
    assert rc == 0
    summary = json.loads(capsys.readouterr().out.strip().splitlines()[-1])
    assert summary["regulatory_features_per_chromosome"] == {"21": 2} and summary["regulatory_shards"] == 1
    assert summary["binding_matrices"] == 1
    with gzip.open(out / "regulatory" / "21" / "1-1000000.json.gz", "rt") as fh:
        kept = json.load(fh)
    assert [f["stable_id"] for f in kept] == ["ENSR_NEAR", "ENSM_IN"]  # the far feature is beyond the 5 kb flank
    assert kept[1]["matrix"] == "ENSPFM1"
    assert json.loads((out / "regulatory" / "binding_matrices.json").read_text()) == {"ENSPFM1": {"transcription_factors": ["TF"], "frequencies": []}}
    assert not (out / "regulatory" / "21" / "1000001-2000000.json.gz").exists()
    assert not (out / "regulatory" / "22").exists()
    info = json.loads((out / "info.json").read_text())
    assert info["regulatory"] is True and info["regulatory_schema"] == 1 and info["cell_types"] == ["A549", "HUVEC"]
    manifest = (out / "CHECKSUMS.sha256").read_text()
    assert "regulatory/21/1-1000000.json.gz" in manifest and "regulatory/binding_matrices.json" in manifest

    # A record reaching the far feature alone keeps that shard and an empty table: the runtime
    # opens the table of every regulatory tree.
    far_vcf = tmp_path / "far.vcf"
    far_vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t500100\t.\tA\tG\t.\t.\t.\n")
    far_out = tmp_path / "far_out"
    assert pjc.main(["--cache", str(cache), "--out", str(far_out), "--vcf", str(far_vcf), "--regulatory"]) == 0
    assert [f["stable_id"] for f in json.loads((far_out / "regulatory" / "21" / "1-1000000.json").read_text())] == ["ENSR_FAR"]
    assert json.loads((far_out / "regulatory" / "binding_matrices.json").read_text()) == {}

    # A --perl-info run carries the schema and cell-type list from the source's own info.json,
    # since Ensembl VEP's info.txt has neither; a source lacking them is refused.
    info_txt = tmp_path / "info.txt"
    info_txt.write_text("species\thomo_sapiens\nassembly\tGRCh38\nregulatory\t1\n")
    perl_out = tmp_path / "perl_out"
    assert pjc.main(["--cache", str(cache), "--out", str(perl_out), "--vcf", str(vcf), "--regulatory",
                     "--perl-info", str(info_txt), "--cache-version", "116"]) == 0
    perl_info = json.loads((perl_out / "info.json").read_text())
    assert perl_info["regulatory_schema"] == 1 and perl_info["cell_types"] == ["A549", "HUVEC"]
    bare = tmp_path / "bare_cache"
    shutil.copytree(cache, bare)
    (bare / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh38", "cache_version": 116, "regulatory": True}))
    assert pjc.main(["--cache", str(bare), "--out", str(tmp_path / "bare_out"), "--vcf", str(vcf), "--regulatory"]) == 2
    assert "lacks regulatory_schema, cell_types" in capsys.readouterr().err

    plain = tmp_path / "plain"
    rc = pjc.main(["--cache", str(cache), "--out", str(plain), "--vcf", str(vcf)])
    assert rc == 0
    assert not (plain / "regulatory").exists()
    info = json.loads((plain / "info.json").read_text())
    assert info["regulatory"] is False and "regulatory_schema" not in info and "cell_types" not in info

    rc = pjc.main(["--cache", str(tmp_path / "none"), "--out", str(tmp_path / "x"), "--vcf", str(vcf), "--regulatory",
                   "--assembly", "GRCh38", "--cache-version", "116"])
    assert rc == 2


def test_prune_writes_the_checksum_manifest_last(tmp_path: Path, capsys):
    """The output ends with `CHECKSUMS.sha256`: one `sha256sum` text-mode line per file of the
    output, the manifest itself excepted, paths relative to the output in byte order, so
    `sha256sum --check` accepts it; `--info-only` refreshes it after rewriting `info.json`, and
    `--manifest-only` writes it alone over an assembled cache directory."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(
        json.dumps([{"stable_id": "ENST_A", "start": "10000", "end": "20000", "dbID": "1"}]))
    (cache / "chr_synonyms.txt").write_text("21\tchr21\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t12000\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--assembly", "GRCh38",
                     "--cache-version", "116", "--gzip"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["manifest_files"] == 3
    checksums = out / "CHECKSUMS.sha256"
    lines = checksums.read_text().splitlines()
    assert [line[66:] for line in lines] == ["chr_synonyms.txt", "info.json", "transcripts/21/1-1000000.json.gz"]
    for line in lines:
        digest, path = line[:64], line[66:]
        assert line[64:66] == "  "
        assert hashlib.sha256((out / path).read_bytes()).hexdigest() == digest
    if shutil.which("sha256sum"):
        check = subprocess.run(["sha256sum", "--check", "--strict", "CHECKSUMS.sha256"], cwd=out, capture_output=True, text=True, check=False)
        assert check.returncode == 0, check.stdout + check.stderr
    before = checksums.read_text()
    assert pjc.main(["--info-only", "--out", str(out), "--assembly", "GRCh38", "--cache-version", "117"]) == 0
    after = checksums.read_text()
    assert after != before
    assert hashlib.sha256((out / "info.json").read_bytes()).hexdigest() in after
    (out / "transcripts" / "21" / "1000001-2000000.json").write_text("[]")
    assert pjc.main(["--manifest-only", "--out", str(out)]) == 0
    assert json.loads(capsys.readouterr().out) == {"manifest_files": 4}
    assert "transcripts/21/1000001-2000000.json" in checksums.read_text()
    assert pjc.main(["--manifest-only", "--out", str(tmp_path / "not-a-cache")]) == 2
    assert "holds no info.json" in capsys.readouterr().err


def test_manifest_lists_files_through_directory_links_in_path_byte_order(tmp_path: Path):
    """The writer lists a file reached through a symbolic link to a directory, as the Rust writer
    does, and orders the lines by the bytes of the paths alone: a name that is a proper prefix of
    another whose next byte sorts below a newline comes first, as it does in the Rust writer's
    byte order. A tree nested deeper than `MANIFEST_MAX_DEPTH`, a link to an ancestor, is refused
    with the tool's error instead of walked without end."""
    out = tmp_path / "out"
    (out / "transcripts" / "21").mkdir(parents=True)
    (out / "info.json").write_text("{}\n", encoding="utf-8")
    (out / "info.json\tcopy").write_text("{}\n", encoding="utf-8")
    (out / "transcripts" / "21" / "1-1000000.json").write_text("[]", encoding="utf-8")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (elsewhere / "1-1000000_var.gz").write_bytes(b"gz")
    (out / "variations").mkdir()
    os.symlink(elsewhere, out / "variations" / "21", target_is_directory=True)
    assert pjc.write_manifest(out) == 4
    lines = (out / "CHECKSUMS.sha256").read_text(encoding="utf-8").splitlines()
    assert [line[66:] for line in lines] == [
        "info.json", "info.json\tcopy", "transcripts/21/1-1000000.json", "variations/21/1-1000000_var.gz",
    ]
    assert lines[3][:64] == hashlib.sha256(b"gz").hexdigest()
    os.symlink(out, out / "transcripts" / "loop", target_is_directory=True)
    with pytest.raises(SystemExit) as refusal:
        pjc.write_manifest(out)
    assert str(refusal.value).startswith("ERROR: [prune_json_cache] ")
    assert str(refusal.value).endswith(f": directories nest more than {pjc.MANIFEST_MAX_DEPTH} levels deep")
    (out / "transcripts" / "loop").unlink()
    os.symlink(out / "absent", out / "dangling")
    with pytest.raises(SystemExit) as refusal:
        pjc.write_manifest(out)
    assert str(refusal.value) == f"ERROR: [prune_json_cache] {out / 'dangling'}: a symbolic link to nothing"


def test_prune_writes_contig_order_and_the_kept_chromosomes_synonyms(tmp_path: Path, capsys):
    """`contig_order` lists the output's own chromosome directories in the converter's order
    (autosomes by number, X, Y, MT, then the rest by name), whatever the source advertised, and
    `chr_synonyms.txt` keeps the lines naming a chromosome the output holds."""
    cache = tmp_path / "cache"
    for chrom in ["21", "X", "3", "KI270728.1", "MT", "22"]:
        (cache / "transcripts" / chrom).mkdir(parents=True)
        (cache / "transcripts" / chrom / "1-1000000.json").write_text(
            json.dumps([{"stable_id": f"ENST_{chrom}", "start": "100", "end": "150", "dbID": "1"}]))
    (cache / "info.json").write_text(json.dumps({
        "species": "homo_sapiens", "assembly": "GRCh38", "cache_version": 116, "source_versions": {},
        "contig_order": ["1", "2", "3", "21", "22", "X", "Y", "MT", "KI270728.1"]}))
    (cache / "chr_synonyms.txt").write_text(
        "21\tNC_000021.9\n21\tchr21\nX\tchrX\n22\tNC_000022.11\nMT\tchrM\nKI270728.1\tchr21_KI270728v1_random\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "X\t120\t.\tA\tG\t.\t.\t.\n21\t120\t.\tA\tG\t.\t.\t.\n3\t120\t.\tA\tG\t.\t.\t.\n"
                   "KI270728.1\t120\t.\tA\tG\t.\t.\t.\nMT\t120\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf)]) == 0
    summary = json.loads(capsys.readouterr().out)
    info = json.loads((out / "info.json").read_text())
    assert info["contig_order"] == ["3", "21", "X", "MT", "KI270728.1"]
    assert (out / "chr_synonyms.txt").read_text().splitlines() == [
        "21\tNC_000021.9", "21\tchr21", "X\tchrX", "MT\tchrM", "KI270728.1\tchr21_KI270728v1_random"]
    assert summary["synonym_lines"] == 5
    # A source without the file writes none, and an info-only output names no contig.
    (cache / "chr_synonyms.txt").unlink()
    bare = tmp_path / "bare"
    assert pjc.main(["--cache", str(cache), "--out", str(bare), "--vcf", str(vcf)]) == 0
    capsys.readouterr()
    assert not (bare / "chr_synonyms.txt").exists()
    info_only = tmp_path / "info_only"
    assert pjc.main(["--info-only", "--out", str(info_only), "--assembly", "GRCh38", "--cache-version", "116"]) == 0
    assert "contig_order" not in json.loads((info_only / "info.json").read_text())


def test_prune_resolves_input_chromosome_names_as_the_runtime_does(tmp_path: Path, capsys):
    """Records on `NC_000021.9` (a synonym the cache's file links), `chr21` and `CHR21` (the
    `chr` prefix removed, in any case, with no synonym line naming them), `Y` (the prefix added,
    for a directory the cache names `chrY`) and `M` (read as `MT`) reach the cache's `21`,
    `chrY` and `MT` directories; a name that resolves to nothing keeps nothing."""
    cache = tmp_path / "cache"
    for chrom in ["21", "chrY", "MT"]:
        (cache / "transcripts" / chrom).mkdir(parents=True)
        (cache / "transcripts" / chrom / "1-1000000.json").write_text(
            json.dumps([{"stable_id": f"ENST_{chrom}", "start": "100", "end": "150", "dbID": "1"},
                        {"stable_id": f"ENST_{chrom}_far", "start": "900000", "end": "900100", "dbID": "2"}]))
    (cache / "chr_synonyms.txt").write_text("21\tNC_000021.9\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "NC_000021.9\t120\t.\tA\tG\t.\t.\t.\nchr21\t130\t.\tA\tG\t.\t.\t.\nCHR21\t140\t.\tA\tG\t.\t.\t.\n"
                   "Y\t120\t.\tA\tG\t.\t.\t.\nM\t120\t.\tA\tG\t.\t.\t.\nNC_000022.11\t120\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--assembly", "GRCh38", "--cache-version", "116"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["spans_by_chromosome"] == {"21": 3, "MT": 1, "NC_000022.11": 1, "chrY": 1}
    assert summary["transcripts_per_chromosome"] == {"21": 1, "MT": 1, "chrY": 1}
    assert [t["stable_id"] for t in json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())] == ["ENST_21"]


def test_prune_resolves_names_through_the_runs_synonyms_file(tmp_path: Path, capsys):
    """Under `--synonyms` the records resolve through that file alone, as Ensembl VEP reads a
    `--synonyms` file in place of the cache's `chr_synonyms.txt`: `MT` reaches the cache's
    `mitochondrion_genome` through the run's line and the cache's own `NC_000021.9` line does
    not resolve; the copied `chr_synonyms.txt` stays the cache's; a missing file is refused."""
    cache = tmp_path / "cache"
    for chrom in ["21", "mitochondrion_genome"]:
        (cache / "transcripts" / chrom).mkdir(parents=True)
        (cache / "transcripts" / chrom / "1-1000000.json").write_text(
            json.dumps([{"stable_id": f"T_{chrom}", "start": "100", "end": "150", "dbID": "1"}]))
    (cache / "chr_synonyms.txt").write_text("21\tNC_000021.9\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "MT\t120\t.\tA\tG\t.\t.\t.\nNC_000021.9\t120\t.\tA\tG\t.\t.\t.\n21\t130\t.\tA\tG\t.\t.\t.\n")
    synonyms = tmp_path / "synonyms.txt"
    synonyms.write_text("mitochondrion_genome\tMT\n")
    common = ["--cache", str(cache), "--vcf", str(vcf), "--assembly", "BDGP6.54", "--cache-version", "116"]
    assert pjc.main([*common, "--out", str(tmp_path / "plain")]) == 0
    plain = json.loads(capsys.readouterr().out)
    assert plain["spans_by_chromosome"] == {"21": 2, "MT": 1}
    assert plain["transcripts_per_chromosome"] == {"21": 1}
    out = tmp_path / "out"
    assert pjc.main([*common, "--out", str(out), "--synonyms", str(synonyms)]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["spans_by_chromosome"] == {"21": 1, "NC_000021.9": 1, "mitochondrion_genome": 1}
    assert summary["transcripts_per_chromosome"] == {"21": 1, "mitochondrion_genome": 1}
    assert (out / "chr_synonyms.txt").read_text() == "21\tNC_000021.9\n"
    assert pjc.main([*common, "--out", str(tmp_path / "out2"), "--synonyms", str(tmp_path / "absent.txt")]) == 2
    assert not (tmp_path / "out2" / "transcripts").exists()


def test_breakend_mate_positions_count_as_spans(tmp_path: Path, capsys):
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\t.\tN\tN]13:5000]\t.\t.\tSVTYPE=BND\n"
        "21\t200\t.\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=22;POS2=700;END2=701\n"
        "21\t300\t.\tN\t<BND>\t.\t.\tSVTYPE=BND;CHR2=21;POS2=900;END2=901\n"
    )
    spans = pjc.record_spans(vcf)
    assert spans["13"] == [(5000, 5000)]
    assert spans["22"] == [(700, 701)]
    assert (300, 901) in spans["21"], "an intra-chromosomal pair keeps the span between the breakends"
    # The summary counts those spans, not records: three records, one mate span each, plus
    # the pair span of the same-chromosome breakend.
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "150", "end": "250"}]))
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out"), "--vcf", str(vcf), "--assembly", "GRCh37", "--cache-version", "116"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["spans_by_chromosome"] == {"13": 1, "21": 5, "22": 1}
    assert summary["transcripts_per_chromosome"] == {"21": 1} and summary["shards"] == 1


def test_prune_reads_a_bgz_vcf_and_refuses_a_missing_one(tmp_path: Path, capsys):
    """A `.bgz` input (BGZF is gzip) yields the spans of the same records stored plain or as
    `.gz`, and a `--vcf` path that is not a file is refused by name with status 2."""
    records = (
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t100\trs1\tA\tG\t.\tPASS\t.\n"
        "21\t24000\tdel1\tT\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=25000\n"
        "21\t300\tbnd1\tN\tN]22:5000]\t.\tPASS\tSVTYPE=BND\n"
    ).encode("utf-8")
    plain = tmp_path / "variants.vcf"
    plain.write_bytes(records)
    gz = tmp_path / "variants.vcf.gz"
    with gzip.open(gz, "wb") as handle:
        handle.write(records)
    bgz = tmp_path / "variants.vcf.bgz"
    with bgz.open("wb") as handle:
        writer = pjc.BgzfWriter(handle)
        writer.write(records)
        writer.close()
    spans = pjc.record_spans(plain)
    assert spans == {"21": [(100, 100), (24000, 25000), (300, 300)], "22": [(5000, 5000)]}
    assert pjc.record_spans(gz) == spans and pjc.record_spans(bgz) == spans
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "150", "end": "250"}]))
    common = ["--cache", str(cache), "--out", str(tmp_path / "out"), "--assembly", "GRCh38", "--cache-version", "116"]
    assert pjc.main([*common, "--vcf", str(bgz)]) == 0
    assert json.loads(capsys.readouterr().out)["spans_by_chromosome"] == {c: len(s) for c, s in sorted(spans.items())}
    missing = tmp_path / "absent.vcf.bgz"
    assert pjc.main([*common, "--vcf", str(bgz), "--vcf", str(missing)]) == 2
    assert capsys.readouterr().err == f"ERROR: [prune_json_cache] --vcf {missing} is not a file\n"
    assert pjc.main([*common, "--vcf", str(tmp_path)]) == 2
    assert capsys.readouterr().err == f"ERROR: [prune_json_cache] --vcf {tmp_path} is not a file\n"


def test_info_json_from_perl_info_txt(tmp_path: Path, capsys):
    info_txt = tmp_path / "info.txt"
    info_txt.write_text(
        "species\thomo_sapiens\nassembly\tGRCh37\nsift\tb\npolyphen\tb\n"
        "source_gencode\tGENCODE 19\nsource_assembly\tGRCh37.p13\nsource_1000genomes\tphase3\n"
        "variation_cols\tchr,variation_name,failed\nvar_type\ttabix\nregulatory\t1\nbuild\t-\n"
        "cell_types\tA549_(m,_58_y),B_cell_(f)\n"
    )
    out = tmp_path / "out"
    rc = pjc.main(["--info-only", "--out", str(out), "--perl-info", str(info_txt), "--cache-version", "115"])
    assert rc == 0
    info = json.loads((out / "info.json").read_text())
    assert info["species"] == "homo_sapiens" and info["assembly"] == "GRCh37" and info["cache_version"] == 115
    assert info["source_versions"] == {"gencode": "GENCODE 19", "assembly": "GRCh37.p13", "1000genomes": "phase3"}
    # No variation slice was written, so the keys that advertise one are dropped.
    assert "variation_cols" not in info and "var_type" not in info
    assert "cell_types" not in info
    # The flag follows the output's own tree, not the source's line: no shard, no regulatory data.
    assert info["regulatory"] is False and "build" not in info
    assert not (out / "transcripts").exists()
    (out / "regulatory" / "21").mkdir(parents=True)
    (out / "regulatory" / "21" / "1-1000000.json").write_text("[]")
    # A tree needs its schema and cell-type list, which no info.txt carries: a refresh keeps
    # the ones the metadata being replaced holds and refuses when neither side has them.
    assert pjc.main(["--info-only", "--out", str(out), "--perl-info", str(info_txt), "--cache-version", "115"]) == 2
    assert "carries regulatory_schema, cell_types" in capsys.readouterr().err
    (out / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh37", "cache_version": 115,
                                               "regulatory": True, "regulatory_schema": 1, "cell_types": ["A549"]}))
    assert pjc.main(["--info-only", "--out", str(out), "--perl-info", str(info_txt), "--cache-version", "115"]) == 0
    refreshed = json.loads((out / "info.json").read_text())
    assert refreshed["regulatory"] is True and refreshed["regulatory_schema"] == 1 and refreshed["cell_types"] == ["A549"]
    # Each missing key is filled alone: an info.txt carrying the schema keeps it as the whole
    # number the runtime reads, and the cell-type list comes from the metadata being replaced.
    schema_txt = tmp_path / "schema_info.txt"
    schema_txt.write_text(info_txt.read_text() + "regulatory_schema\t1\n")
    (out / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh37", "cache_version": 115,
                                               "regulatory": True, "cell_types": ["A549"]}))
    assert pjc.main(["--info-only", "--out", str(out), "--perl-info", str(schema_txt), "--cache-version", "115"]) == 0
    mixed = json.loads((out / "info.json").read_text())
    assert mixed["regulatory"] is True and mixed["regulatory_schema"] == 1 and mixed["cell_types"] == ["A549"]
    assert isinstance(mixed["regulatory_schema"], int)
    bad_schema_txt = tmp_path / "bad_schema_info.txt"
    bad_schema_txt.write_text(info_txt.read_text() + "regulatory_schema\tone\n")
    assert pjc.main(["--info-only", "--out", str(out), "--perl-info", str(bad_schema_txt), "--cache-version", "115"]) == 2
    assert capsys.readouterr().err == f"ERROR: [prune_json_cache] --perl-info {bad_schema_txt}: regulatory_schema 'one' is not a whole number\n"
    # --regulatory names a slice to cut, which an --info-only run does not.
    assert pjc.main(["--info-only", "--regulatory", "--out", str(out), "--perl-info", str(info_txt), "--cache-version", "115"]) == 2


VAR_COLS = ["chr", "variation_name", "failed", "somatic", "start", "end", "allele_string", "strand", "clin_sig",
            "phenotype_or_disease", "clin_sig_allele", "clin_sig_ref_allele", "clinical_impact", "pubmed",
            "var_synonyms", "AF", "gnomADg"]


def _all_vars_line(chrom: str, name: str, start: int, alleles: str = "A/G") -> bytes:
    fields = [chrom, name, ".", ".", str(start), ".", alleles] + ["."] * (len(VAR_COLS) - 7)
    return ("\t".join(fields) + "\n").encode()


def _write_all_vars(path: Path, lines: list[bytes]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        writer = pjc.BgzfWriter(raw)
        for line in lines:
            writer.write(line)
        writer.close()


def _bgzf_blocks(path: Path) -> dict[int, bytes]:
    """Decompressed bgzf blocks by their file offset."""
    data = path.read_bytes()
    blocks, offset = {}, 0
    while offset < len(data):
        assert data[offset:offset + 4] == b"\x1f\x8b\x08\x04", "a bgzf member with the extra field"
        block_size = int.from_bytes(data[offset + 16:offset + 18], "little") + 1
        payload = data[offset + 18:offset + block_size - 8]
        blocks[offset] = zlib.decompress(payload, -15)
        offset += block_size
    return blocks


def _read_chunk(blocks: dict[int, bytes], begin: int, end: int) -> bytes:
    offsets = sorted(blocks)
    out = b""
    for block_offset in offsets:
        if block_offset < begin >> 16 or block_offset > end >> 16:
            continue
        block = blocks[block_offset]
        lo = (begin & 0xFFFF) if block_offset == begin >> 16 else 0
        hi = (end & 0xFFFF) if block_offset == end >> 16 else len(block)
        out += block[lo:hi]
    return out


def _parse_csi(path: Path) -> tuple[dict, dict[int, list[tuple[int, int, int]]]]:
    """The header fields and, for the single reference, bin -> [(loffset, begin, end), ...]."""
    data = gzip.open(path, "rb").read()
    assert data[:4] == b"CSI\x01"
    min_shift, depth, l_aux = struct.unpack("<iii", data[4:16])
    aux = data[16:16 + l_aux]
    fmt, seq_col, beg_col, end_col, meta, skip, l_nm = struct.unpack("<7i", aux[:28])
    names = aux[28:28 + l_nm].split(b"\0")[:-1]
    offset = 16 + l_aux
    n_ref, = struct.unpack("<i", data[offset:offset + 4])
    offset += 4
    header = {"min_shift": min_shift, "depth": depth, "format": fmt, "seq_col": seq_col, "beg_col": beg_col,
              "end_col": end_col, "meta": meta, "skip": skip, "names": names, "n_ref": n_ref}
    bins: dict[int, list[tuple[int, int, int]]] = {}
    for _ in range(n_ref):
        n_bin, = struct.unpack("<i", data[offset:offset + 4])
        offset += 4
        for _ in range(n_bin):
            b, loffset, n_chunk = struct.unpack("<IQi", data[offset:offset + 16])
            offset += 16
            for _ in range(n_chunk):
                begin, end = struct.unpack("<QQ", data[offset:offset + 16])
                offset += 16
                bins.setdefault(b, []).append((loffset, begin, end))
    return header, bins


def _reg2bins(beg: int, end: int, min_shift: int, depth: int) -> list[int]:
    found = []
    for level in range(depth + 1):
        shift = min_shift + 3 * (depth - level)
        first = ((1 << (3 * level)) - 1) // 7
        found.extend(range(first + (beg >> shift), first + ((end - 1) >> shift) + 1))
    return found


def _csi_query(gz: Path, start: int) -> list[str]:
    """The record names an index-driven reader returns for a one-base query at `start`."""
    header, bins = _parse_csi(Path(str(gz) + ".csi"))
    blocks = _bgzf_blocks(gz)
    names = []
    for b in _reg2bins(start - 1, start, header["min_shift"], header["depth"]):
        for _loffset, begin, end in bins.get(b, []):
            for line in _read_chunk(blocks, begin, end).split(b"\n"):
                fields = line.split(b"\t")
                if len(fields) > 4 and int(fields[header["beg_col"] - 1]) == start:
                    names.append(fields[1].decode())
    return names


def test_prune_variation_slice_reindexes(tmp_path: Path):
    src = tmp_path / "perl-cache"
    lines = [_all_vars_line("21", "rs1", 100), _all_vars_line("21", "rs2", 150), _all_vars_line("21", "COSV9900001", 200, "COSMIC_MUTATION"),
             _all_vars_line("21", "rs3", 200), _all_vars_line("21", "rs3b", 201), _all_vars_line("21", "CM990001", 201, "HGMD_MUTATION"),
             _all_vars_line("21", "rs4", 202),
             _all_vars_line("21", "rs5", 5_000_000), _all_vars_line("21", "rs6", 5_000_001), _all_vars_line("21", "rs7", 5_000_003),
             _all_vars_line("21", "rs8", 9_000_000)]
    # A long record spans two bgzf blocks, as a record with many synonyms does.
    long_record = ("\t".join(["21", "rs9", ".", ".", "9000001", ".", "A/T"] + ["."] * 7 + ["x" * 70_000] + ["."] * 2) + "\n").encode()
    lines.append(long_record)
    lines.append(_all_vars_line("21", "rs10", 9_000_002))
    _write_all_vars(src / "21" / "all_vars.gz", lines)
    _write_all_vars(src / "22" / "all_vars.gz", [_all_vars_line("22", "rs22", 100)])
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "150", "end": "250"}]))
    info_txt = tmp_path / "info.txt"
    info_txt.write_text("species\thomo_sapiens\nassembly\tGRCh38\nvariation_cols\t" + ",".join(VAR_COLS) + "\nvar_type\ttabix\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "21\t200\t.\tA\tG\t.\t.\t.\n"
        "21\t5000002\t.\tA\tG\t.\t.\t.\n"
        "21\t9000001\t.\tA\tT\t.\t.\t.\n"
    )
    out = tmp_path / "out"
    rc = pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116",
                   "--variations", str(src), "--exclude-variation", r"^COSV\d", "--exclude-variation", r"^C[MID]\d{6}"])
    assert rc == 0
    gz = out / "variations" / "21" / "all_vars.gz"
    kept = [line.split(b"\t")[1].decode() for line in gzip.open(gz, "rb").read().splitlines()]
    # Starts within one base of a record span survive (199..201, 5000001..5000003, 9000000..9000002),
    # the COSMIC and HGMD records are dropped, chromosome 22 has no record span.
    assert kept == ["rs3", "rs3b", "rs6", "rs7", "rs8", "rs9", "rs10"]
    assert not (out / "variations" / "22").exists()
    header, bins = _parse_csi(Path(str(gz) + ".csi"))
    assert header["names"] == [b"21"] and header["n_ref"] == 1
    assert (header["seq_col"], header["beg_col"], header["end_col"]) == (1, 5, 5), "the index names the start column as both bounds"
    assert (header["format"], header["meta"], header["skip"]) == (0, ord("#"), 0)
    assert (header["min_shift"], header["depth"]) == (pjc.CSI_MIN_SHIFT, pjc.CSI_DEPTH)
    meta_bin = pjc.bin_first(pjc.CSI_DEPTH + 1) + 1
    assert bins[meta_bin][1][1:] == (7, 0), "the pseudo-bin counts the records written"
    # The first bin's `loffset` is the first record's virtual offset, 0: htslib's `hts_itr_query`
    # starts each chunk at the `loffset` of the bin covering the query start, so a later offset
    # hides the first record from a query at its position.
    first_bin = pjc.bin_first(pjc.CSI_DEPTH) + (199 >> pjc.CSI_MIN_SHIFT)
    assert bins[first_bin][0][0] == 0
    for start, expected in [(200, ["rs3"]), (201, ["rs3b"]), (202, []), (5_000_001, ["rs6"]), (5_000_003, ["rs7"]), (9_000_001, ["rs9"]),
                            (9_000_002, ["rs10"]), (100, []), (150, [])]:
        assert _csi_query(gz, start) == expected, start
    info = json.loads((out / "info.json").read_text())
    assert info["variation_cols"] == VAR_COLS and info["var_type"] == "tabix"
    # Without --variations the same run writes no slice and no variation keys.
    out2 = tmp_path / "out2"
    assert pjc.main(["--cache", str(cache), "--out", str(out2), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116"]) == 0
    info2 = json.loads((out2 / "info.json").read_text())
    assert "variation_cols" not in info2 and "var_type" not in info2 and not (out2 / "variations").exists()
    # A --variations directory without a file for any record chromosome, and a negative flank,
    # are refused rather than written as an empty slice.
    empty = tmp_path / "empty-variations"
    (empty / "22").mkdir(parents=True)
    out3 = tmp_path / "out3"
    assert pjc.main(["--cache", str(cache), "--out", str(out3), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116",
                     "--variations", str(empty)]) == 2
    assert not (out3 / "info.json").exists() and not (out3 / "transcripts").exists()
    assert pjc.main(["--cache", str(cache), "--out", str(out3), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116",
                     "--variations", str(src), "--variation-flank", "-1"]) == 2


def test_prune_variation_slice_takes_the_columns_of_the_variation_cache(tmp_path: Path):
    """A slice cut from a native cache beside an Ensembl cache's metadata declares the native
    cache's columns: `--variation-info` replaces `variation_cols` and `var_type` of the info being
    written, and the index is built on the start column of those columns. The option without
    `--variations`, either option under `--info-only`, an info.json without a `variation_cols`
    list, or one whose top level is not an object, is refused; so is a slice whose records carry
    a field count other than the declared columns, or whose file is not gzip or does not inflate,
    before anything is written, while a chromosome file with no record is written as the empty
    indexed slice the Ensembl VEP dumper writes, not refused."""
    native = tmp_path / "native"
    native_cols = VAR_COLS[:15] + ["AF"]
    line = ("\t".join(["21", "rs1", ".", ".", "200", ".", "A/G"] + ["."] * (len(native_cols) - 7)) + "\n").encode()
    _write_all_vars(native / "variations" / "21" / "all_vars.gz", [line])
    (native / "info.json").write_text(json.dumps({"species": "homo_sapiens", "assembly": "GRCh38", "cache_version": 116,
                                                   "variation_cols": native_cols, "var_type": "tabix", "source_versions": {}}))
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "150", "end": "250"}]))
    info_txt = tmp_path / "info.txt"
    info_txt.write_text("species\thomo_sapiens\nassembly\tGRCh38\nsource_dbSNP\t156\nvariation_cols\t" + ",".join(VAR_COLS) + "\nvar_type\ttabix\n")
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t200\t.\tA\tG\t.\t.\t.\n")
    common = ["--cache", str(cache), "--vcf", str(vcf), "--perl-info", str(info_txt), "--cache-version", "116"]
    native_info = ["--variation-info", str(native / "info.json")]
    out = tmp_path / "out"
    assert pjc.main([*common, "--out", str(out), "--variations", str(native / "variations"), *native_info]) == 0
    info = json.loads((out / "info.json").read_text())
    assert info["variation_cols"] == native_cols and info["var_type"] == "tabix"
    assert info["source_versions"] == {"dbSNP": "156"}, "the rest of the metadata is still the Ensembl cache's"
    gz = out / "variations" / "21" / "all_vars.gz"
    assert _csi_query(gz, 200) == ["rs1"]
    assert pjc.main([*common, "--out", str(tmp_path / "out2"), *native_info]) == 2
    (tmp_path / "bare.json").write_text(json.dumps({"species": "homo_sapiens"}))
    assert pjc.main([*common, "--out", str(tmp_path / "out3"), "--variations", str(native / "variations"), "--variation-info", str(tmp_path / "bare.json")]) == 2
    assert not (tmp_path / "out3" / "info.json").exists()
    (tmp_path / "list.json").write_text(json.dumps([]))
    assert pjc.main([*common, "--out", str(tmp_path / "out4"), "--variations", str(native / "variations"), "--variation-info", str(tmp_path / "list.json")]) == 2
    assert not (tmp_path / "out4" / "info.json").exists()
    # The Ensembl cache's columns over the native files: the record width disagrees, so the
    # declaration is refused rather than written over records it does not describe.
    assert pjc.main([*common, "--out", str(tmp_path / "out5"), "--variations", str(native / "variations")]) == 2
    assert not (tmp_path / "out5").exists() or not any((tmp_path / "out5").iterdir())
    assert pjc.main(["--info-only", "--out", str(tmp_path / "out6"), "--perl-info", str(info_txt), "--cache-version", "116", *native_info]) == 2
    assert not (tmp_path / "out6" / "info.json").exists()
    assert pjc.main(["--info-only", "--out", str(tmp_path / "out7"), "--perl-info", str(info_txt), "--cache-version", "116",
                     "--variations", str(native / "variations"), *native_info]) == 2
    assert not (tmp_path / "out7" / "info.json").exists()
    # A chromosome file with no record yields the empty bgzf file and an index naming the sequence
    # with no bin, and the metadata keeps the variation keys; a file that is not gzip is refused.
    empty = tmp_path / "empty"
    _write_all_vars(empty / "variations" / "21" / "all_vars.gz", [])
    assert pjc.main([*common, "--out", str(tmp_path / "out8"), "--variations", str(empty / "variations"), *native_info]) == 0
    empty_slice = tmp_path / "out8" / "variations" / "21" / "all_vars.gz"
    assert empty_slice.read_bytes() == pjc.BGZF_EOF
    header, bins = _parse_csi(Path(str(empty_slice) + ".csi"))
    assert header["names"] == [b"21"] and header["n_ref"] == 1 and bins == {}
    assert _csi_query(empty_slice, 200) == []
    info8 = json.loads((tmp_path / "out8" / "info.json").read_text())
    assert info8["variation_cols"] == native_cols and info8["var_type"] == "tabix"
    plain = tmp_path / "plain"
    (plain / "variations" / "21").mkdir(parents=True)
    (plain / "variations" / "21" / "all_vars.gz").write_bytes(b"21\trs1\t.\t.\t200\n")
    assert pjc.main([*common, "--out", str(tmp_path / "out9"), "--variations", str(plain / "variations"), *native_info]) == 2
    assert not (tmp_path / "out9" / "transcripts").exists()
    # A gzip header over a corrupt deflate stream raises zlib.error, not OSError or EOFError.
    corrupt = tmp_path / "corrupt"
    (corrupt / "variations" / "21").mkdir(parents=True)
    (corrupt / "variations" / "21" / "all_vars.gz").write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03" + b"\xff\xff\xff\xff\x00garbage" * 4)
    assert pjc.main([*common, "--out", str(tmp_path / "out10"), "--variations", str(corrupt / "variations"), *native_info]) == 2
    assert not (tmp_path / "out10" / "transcripts").exists()


def test_bgzf_writer_blocks_and_virtual_offsets(tmp_path: Path):
    path = tmp_path / "x.gz"
    with path.open("wb") as raw:
        writer = pjc.BgzfWriter(raw)
        assert writer.tell() == 0
        writer.write(b"a" * 1000)
        assert writer.tell() == 1000
        writer.write(b"b" * pjc.BGZF_BLOCK_BYTES)
        first_block_size = writer.block_offset
        assert first_block_size > 0 and writer.tell() == (first_block_size << 16) | 1000
        writer.close()
    blocks = _bgzf_blocks(path)
    assert [len(b) for b in blocks.values()] == [pjc.BGZF_BLOCK_BYTES, 1000, 0]
    assert gzip.open(path, "rb").read() == b"a" * 1000 + b"b" * pjc.BGZF_BLOCK_BYTES
    assert path.read_bytes().endswith(pjc.BGZF_EOF)


def test_reg2bin_follows_the_binning_scheme():
    # Bins at the finest level start at (8^8 - 1) / 7 and cover 16 kb windows; a span across a
    # window boundary moves up one level.
    finest = pjc.bin_first(pjc.CSI_DEPTH)
    assert pjc.reg2bin(0, 1) == finest
    assert pjc.reg2bin(16_384, 16_385) == finest + 1
    assert pjc.reg2bin(16_383, 16_385) == pjc.bin_first(pjc.CSI_DEPTH - 1)
    assert pjc.bin_level(0) == 0 and pjc.bin_level(finest) == pjc.CSI_DEPTH
    assert pjc.bin_bottom(finest + 3) == 3 and pjc.bin_bottom(pjc.bin_first(pjc.CSI_DEPTH - 1) + 1) == 8
    assert pjc.variation_windows([(100, 100), (103, 110), (300, 300)], 1) == [(99, 111), (299, 301)]
    assert pjc.variation_windows([(1, 5)], 3) == [(1, 8)]


def test_gzip_option_writes_gz_shards_readable_as_source(tmp_path: Path):
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([{"stable_id": "ENST_A", "start": "10000", "end": "20000"}]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t12000\t.\tA\tG\t.\t.\t.\n")
    out1 = tmp_path / "out1"
    assert pjc.main(["--cache", str(cache), "--out", str(out1), "--vcf", str(vcf), "--gzip", "--assembly", "GRCh37", "--cache-version", "115"]) == 0
    gz = out1 / "transcripts" / "21" / "1-1000000.json.gz"
    assert gz.is_file()
    assert pjc.read_shard(gz)[0]["stable_id"] == "ENST_A"
    # A gzipped cache is itself a valid pruning source.
    out2 = tmp_path / "out2"
    assert pjc.main(["--cache", str(out1), "--out", str(out2), "--vcf", str(vcf), "--assembly", "GRCh37", "--cache-version", "115"]) == 0
    assert json.loads((out2 / "transcripts" / "21" / "1-1000000.json").read_text())[0]["stable_id"] == "ENST_A"


def test_peptide_is_kept_and_cannot_be_dropped(tmp_path: Path):
    """The reference peptide's sequence edits are read from the cached peptide, so a pruned
    cache keeps it and `--drop peptide` is refused."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    tx = {"stable_id": "ENST_A", "start": "10000", "end": "20000",
          "variation_effect_feature_cache": {"peptide": "MLW", "sorted_exons": [], "translateable_seq": "GTGTTATGG"}}
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([tx]))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t12000\t.\tA\tG\t.\t.\t.\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--drop", "sorted_exons", "--assembly", "GRCh37", "--cache-version", "115"]) == 0
    kept = json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())[0]["variation_effect_feature_cache"]
    assert kept["peptide"] == "MLW" and "sorted_exons" not in kept
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out2"), "--vcf", str(vcf), "--drop", "peptide", "--assembly", "GRCh37", "--cache-version", "115"]) == 2


def test_drop_except_keeps_the_dropped_keys_on_the_named_transcripts(tmp_path: Path, capsys):
    """`--drop-except <file>` names the transcripts that keep every `--drop` key, so a corpus
    carries protein features and prediction matrices for the transcripts its records test and
    for no other; the option needs a `--drop` key to apply to, and a file that cannot be read
    or decoded is refused with the path and the reason, not a traceback."""
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    features = [{"start": 1, "end": 9, "hseqname": "PF00080", "analysis": "Pfam"}]
    matrices = {"sift": {"analysis": "sift", "predictions_data": "H4sI"}}
    vefc = {"peptide": "MLW", "sorted_exons": [], "protein_features": features, "protein_function_predictions": matrices}
    transcripts = [
        {"stable_id": "ENST_KEPT", "start": "10000", "end": "20000", "variation_effect_feature_cache": vefc},
        {"stable_id": "ENST_DROPPED", "start": "12000", "end": "22000", "variation_effect_feature_cache": vefc},
    ]
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps(transcripts))
    vcf = tmp_path / "v.vcf"
    vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t15000\t.\tA\tG\t.\t.\t.\n")
    keep = tmp_path / "keep.txt"
    keep.write_text("ENST_KEPT\n\n")
    out = tmp_path / "out"
    assert pjc.main(["--cache", str(cache), "--out", str(out), "--vcf", str(vcf), "--drop", "sorted_exons",
                     "--drop", "protein_features", "--drop", "protein_function_predictions",
                     "--drop-except", str(keep), "--assembly", "GRCh38", "--cache-version", "116"]) == 0
    assert json.loads(capsys.readouterr().out)["drop_except"] == 1
    by_id = {t["stable_id"]: t["variation_effect_feature_cache"]
             for t in json.loads((out / "transcripts" / "21" / "1-1000000.json").read_text())}
    assert by_id["ENST_KEPT"]["protein_features"] == features
    assert by_id["ENST_KEPT"]["protein_function_predictions"] == matrices
    assert by_id["ENST_KEPT"]["sorted_exons"] == []
    assert set(by_id["ENST_DROPPED"]) == {"peptide"}
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out2"), "--vcf", str(vcf),
                     "--drop-except", str(keep), "--assembly", "GRCh38", "--cache-version", "116"]) == 2
    assert "--drop-except needs at least one --drop key" in capsys.readouterr().err
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out3"), "--vcf", str(vcf), "--drop", "sorted_exons",
                     "--drop-except", str(tmp_path / "missing.txt"), "--assembly", "GRCh38", "--cache-version", "116"]) == 2
    assert f"--drop-except {tmp_path / 'missing.txt'}: " in capsys.readouterr().err
    not_utf8 = tmp_path / "latin1.txt"
    not_utf8.write_bytes(b"\xffENST_KEPT\n")
    assert pjc.main(["--cache", str(cache), "--out", str(tmp_path / "out4"), "--vcf", str(vcf), "--drop", "sorted_exons",
                     "--drop-except", str(not_utf8), "--assembly", "GRCh38", "--cache-version", "116"]) == 2
    assert f"--drop-except {not_utf8}: 'utf-8' codec can't decode" in capsys.readouterr().err


def _corpora() -> list[Path]:
    root = Path(__file__).resolve().parents[2] / "tests" / "golden"
    return sorted(p for p in root.glob("*/*") if (p / "manifest.json").is_file())


@pytest.mark.parametrize("corpus", _corpora(), ids=lambda p: f"{p.parent.name}/{p.name}")
def test_committed_expected_files_match_their_provenance_digests(corpus: Path):
    """The expected outputs are stored gzipped; their uncompressed bytes must hash to
    what provenance.json recorded when VEP produced them, for the manifest's own run
    and for every flag set it names, each under its own record in provenance.json's
    `flag_sets`. A tool corpus (a manifest with `tool_cases`) records a companion
    tool's outputs instead: its `sha256` map is keyed by each file's path below
    `expected/`, every case's record repeats its digest, and it carries a cache only
    when its manifest names one (a tool that annotates, like `vep recode`); a tool
    that reads an annotated output alone holds none."""
    prov = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))
    manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
    runs = [("expected", prov["sha256"])]
    for flag_set in manifest.get("flag_sets") or []:
        set_name = flag_set["name"]
        record = (prov.get("flag_sets") or {}).get(set_name)
        assert record is not None, f"{corpus.name}: flag set {set_name} has no provenance record"
        runs.append((flag_set.get("expected_dir", f"expected_{set_name}"), record["sha256"]))
    for expected_dir, digests in runs:
        for name, digest in digests.items():
            gz = corpus / expected_dir / (name + ".gz")
            plain = corpus / expected_dir / name
            data = gzip.open(gz, "rb").read() if gz.is_file() else plain.read_bytes()
            assert hashlib.sha256(data).hexdigest() == digest, f"{corpus.name}/{expected_dir}/{name} differs from its recorded digest"
    # The directory is the assembly, optionally suffixed (`GRCh37-hgvs`, `R64-1-1-yeast`);
    # an assembly may itself carry hyphens, so the name is matched against the manifest's.
    assembly = manifest["assembly"]
    assert corpus.name == assembly or corpus.name.startswith(assembly + "-"), (
        f"{corpus.name}: the directory does not start with the manifest's assembly {assembly}"
    )
    # The set directory names the Ensembl VEP release the reference was generated with; the
    # manifest's release is the cache version both engines passed, the same number for an
    # Ensembl species and the division's own for an Ensembl Genomes one (Ensembl Plants 63
    # runs on the Ensembl 116 API), which the record carries as cache_version; a tool corpus
    # runs no annotation and names the set's release.
    if "tool_cases" in manifest:
        assert manifest["release"] == corpus.parent.name
    else:
        vep_release = prov["vep_version_reported"].rsplit(":", 1)[-1].strip().split(".")[0]
        assert vep_release == corpus.parent.name, (
            f"{corpus.name}: the reference was generated with Ensembl VEP {vep_release}, "
            f"not the set's {corpus.parent.name}"
        )
        assert manifest["release"] == str(prov["cache_version"]), (
            f"{corpus.name}: the manifest's release {manifest['release']} is not the record's "
            f"cache_version {prov['cache_version']}"
        )
    if "tool_cases" in manifest:
        for case in manifest["tool_cases"]:
            record = prov["tool_cases"][case["name"]]
            assert record["expected"] == case["expected"], f"{corpus.name}: tool case {case['name']}"
            digest = hashlib.sha256((corpus / case["expected"]).read_bytes()).hexdigest()
            assert digest == record["sha256"] == prov["sha256"][case["expected"].removeprefix("expected/")], (
                f"{corpus.name}/{case['expected']} differs from its recorded digest"
            )
        if "cache" not in manifest:
            assert not (corpus / "json_cache").exists(), f"{corpus.name}: a tool corpus without a cache in its manifest holds none"
            return
    cache = bgc.corpus_cache_dir(corpus, manifest)
    info = json.loads((cache / "info.json").read_text(encoding="utf-8"))
    assert info["assembly"] == assembly and info["source_versions"], "info.json must carry the cache source versions"
    if "fasta" in manifest:
        assert (corpus / manifest["fasta"]).is_file() and (corpus / (manifest["fasta"].removesuffix(".gz") + ".fai")).is_file()


def _corpora_with_peptides() -> list[Path]:
    """The corpora whose cache directory, located by the manifest's `cache_layout`, holds
    `transcripts/`."""
    corpora = []
    for corpus in _corpora():
        manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
        cache = bgc.corpus_cache_dir(corpus, manifest)
        if (cache / "transcripts").is_dir():
            corpora.append(corpus)
    return corpora


@pytest.mark.parametrize("corpus", _corpora_with_peptides(), ids=lambda p: f"{p.parent.name}/{p.name}")
def test_committed_caches_carry_a_peptide_for_every_coding_transcript(corpus: Path):
    """Every committed golden cache carries the peptide beside each CDS, so the format-parity
    tests exercise the sequence-edit overlay the measurement caches exercise."""
    manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
    cache = bgc.corpus_cache_dir(corpus, manifest)
    missing = 0
    for shard in (cache / "transcripts").glob("*/*"):
        for t in pjc.read_shard(shard):
            v = (t or {}).get("variation_effect_feature_cache") or {}
            if v.get("translateable_seq") and not v.get("peptide"):
                missing += 1
    assert missing == 0, f"{corpus}: {missing} coding transcripts without a peptide"


def test_corrected_behaviour_classes_agree_across_the_builder_the_harness_and_the_comparator():
    """The classes of the corrected behaviours are one list in this builder's
    `CORRECTED_BEHAVIOUR_CLASSES` and the format comparator's registry, field, warnings and
    manifest-applied classes, and the golden harness's constant of the same name
    (crates/vep-cli/tests/common/mod.rs) names every class of that list a committed corpus
    manifest's `masks` carries: a class the harness lacks is one no golden corpus names in its
    `masks`, the one place the harness refuses an unknown name, so the manifests, the harness
    and the comparator accept the same names wherever one is used."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    comparator = {*cvf.REGISTRY_CLASSES, *cvf.FIELD_CLASSES, *cvf.PLUGIN_WARNING_CLASSES, *cvf.MANIFEST_CLASSES}
    assert comparator == set(bgc.CORRECTED_BEHAVIOUR_CLASSES)
    harness = (HERE.parent.parent / "crates" / "vep-cli" / "tests" / "common" / "mod.rs").read_text(encoding="utf-8")
    start = harness.index("const CORRECTED_BEHAVIOUR_CLASSES: &[&str] = &[")
    block = harness[start:harness.index("];", start)]
    names = {line.strip().strip(",").strip('"') for line in block.splitlines()[1:] if line.strip().startswith('"')}
    assert names <= set(bgc.CORRECTED_BEHAVIOUR_CLASSES)
    named_in_manifests: set[str] = set()
    for manifest in sorted((HERE.parent.parent / "tests" / "golden").rglob("manifest.json")):
        documented = json.loads(manifest.read_text(encoding="utf-8"))
        for run in (documented, *documented.get("flag_sets", [])):
            named_in_manifests.update(run.get("masks", []))
    assert not (set(bgc.CORRECTED_BEHAVIOUR_CLASSES) - names) & named_in_manifests
    assert len(bgc.CORRECTED_BEHAVIOUR_CLASSES) == len(set(bgc.CORRECTED_BEHAVIOUR_CLASSES))
    assert bgc.masks_shape_error(list(bgc.CORRECTED_BEHAVIOUR_CLASSES), "masks") is None


def _harness_mask_classes() -> list[str]:
    """The classes the golden harness's `MASK_CLASSES` block (crates/vep-cli/tests/common/mod.rs)
    names: the block lists constants, each resolved to the string its `pub const` declares."""
    harness = (HERE.parent.parent / "crates" / "vep-cli" / "tests" / "common" / "mod.rs").read_text(encoding="utf-8")
    constants = dict(re.findall(r'pub const ([A-Z_]+): &str = "([a-z0-9_]+)";', harness))
    start = harness.index("const MASK_CLASSES: &[&str] = &[")
    block = harness[start:harness.index("];", start)]
    names = [line.strip().rstrip(",") for line in block.splitlines()[1:] if line.strip()]
    unresolved = [name for name in names if name not in constants]
    assert not unresolved, f"the harness's MASK_CLASSES names constants without a string: {unresolved}"
    return [constants[name] for name in names]


def test_harness_mask_classes_are_the_builder_mask_classes():
    """The golden harness's `MASK_CLASSES`, the classes that change how it compares, are the
    builder's `REGULATORY_FEATURE_ROWS`, `ORDER_ONLY` and `HARNESS_ONLY_MASK_CLASSES`, so a class
    one side names and the other lacks fails here: a harness class the builder lacks would
    otherwise refuse a committed manifest at the next `classify`, `check` or `select`, and a
    builder class the harness lacks would be accepted into a manifest the harness then refuses."""
    harness = _harness_mask_classes()
    builder = [bgc.REGULATORY_FEATURE_ROWS, bgc.ORDER_ONLY, *bgc.HARNESS_ONLY_MASK_CLASSES]
    assert set(harness) - set(bgc.MASK_CLASSES) == set(), "harness mask classes the builder lacks"
    assert set(harness) == set(builder), "harness and builder harness-only mask classes differ"
    assert len(harness) == len(set(harness)) and len(builder) == len(set(builder))


def test_every_committed_manifest_passes_the_builder_shape_check():
    """Every committed corpus manifest under `tests/golden/`, with each of its flag sets, is a
    shape `manifest_shape_error` accepts, so `classify`, `check` and `select` read every corpus
    the harness reads. A tool corpus (a manifest with `tool_cases`, which the harness and the
    builder never read since it has no `variants.vcf`) names a companion tool's mask classes,
    the format comparator's recoder classes, and is held to the same shape with those set aside."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    recoder_classes = {*cvf.RECODER_MASKS, cvf.ARCHIVE_DBSNP_SYNONYMS, cvf.IDENTIFIER_INDEL_NOTATIONS}
    manifests = sorted((HERE.parent.parent / "tests" / "golden").rglob("manifest.json"))
    assert manifests, "no committed manifests found"
    refused: list[str] = []
    for path in manifests:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        corpus = f"{path.parent.parent.name}/{path.parent.name}"
        if "tool_cases" in manifest:
            assert not (path.parent / "variants.vcf").exists(), f"{corpus}: a tool corpus holds no variants.vcf"
            unknown = sorted(set(manifest.get("masks", [])) - recoder_classes)
            if unknown:
                refused.append(f"{corpus}: masks names classes the recoder comparison lacks {unknown}")
            manifest = {key: value for key, value in manifest.items() if key != "masks"}
        else:
            assert (path.parent / "variants.vcf").is_file(), f"{corpus}: a corpus the harness reads holds variants.vcf"
        reason = bgc.manifest_shape_error(manifest)
        if reason:
            refused.append(f"{corpus}: {reason}")
    assert not refused, "committed manifests the builder refuses:\n" + "\n".join(refused)


def _unexplained_entries(manifest: dict) -> list[str]:
    """The places of a manifest, its own run and each flag set's, that record a row or a field
    under `unexplained_residual`: an entry of a `classify` list or a key of a summary."""
    runs = [("the manifest's run", manifest)]
    runs += [(f"flag set {fs.get('name')}", fs) for fs in manifest.get("flag_sets") or []]
    found: list[str] = []
    for owner, run in runs:
        for key in ("divergences", "vep_rs_only_tuples", "field_divergences"):
            count = sum(1 for entry in run.get(key) or [] if entry.get("expected_divergence") == bgc.UNEXPLAINED)
            if count:
                found.append(f"{owner}: {count} in {key}")
        for key in ("divergence_summary", "vep_rs_only_summary", "field_divergence_summary"):
            if bgc.UNEXPLAINED in (run.get(key) or {}):
                found.append(f"{owner}: {key}")
    return found


def test_no_committed_manifest_carries_an_unexplained_residual():
    """No committed corpus manifest records a row or a field under `unexplained_residual`, in
    its own run or in a flag set. The golden harness applies a `divergences`,
    `vep_rs_only_tuples` or `field_divergences` entry whatever its class, so this check is the
    one that refuses a difference `classify` could not explain."""
    entry = {"location": "21:100", "expected_divergence": bgc.UNEXPLAINED}
    assert _unexplained_entries({"divergences": [{"expected_divergence": "splice_family_swap"}]}) == []
    assert _unexplained_entries({
        "divergences": [entry], "divergence_summary": {bgc.UNEXPLAINED: 1},
        "flag_sets": [{
            "name": "hgvs", "vep_rs_only_tuples": [entry], "field_divergences": [entry, entry],
            "vep_rs_only_summary": {bgc.UNEXPLAINED: 1}, "field_divergence_summary": {bgc.UNEXPLAINED: 2},
        }],
    }) == [
        "the manifest's run: 1 in divergences", "the manifest's run: divergence_summary",
        "flag set hgvs: 1 in vep_rs_only_tuples", "flag set hgvs: 2 in field_divergences",
        "flag set hgvs: vep_rs_only_summary", "flag set hgvs: field_divergence_summary",
    ]
    manifests = sorted((HERE.parent.parent / "tests" / "golden").rglob("manifest.json"))
    assert manifests, "no committed manifests found"
    unexplained = [
        f"{path.parent.parent.name}/{path.parent.name}: {place}"
        for path in manifests
        for place in _unexplained_entries(json.loads(path.read_text(encoding="utf-8")))
    ]
    assert not unexplained, "committed manifests that carry unexplained_residual:\n" + "\n".join(unexplained)


def test_classify_names_the_stop_codon_protein_comparison_class(tmp_path: Path):
    """An `HGVSp` at the stop codon of one protein, Ensembl VEP's description read from the
    codon window and vep-rs's the one the two proteins give, is
    `hgvsp_stop_codon_protein_comparison` in each of its nine shapes; the nearest shapes stay
    an `unexplained_residual`, and the format comparator's class of the same name reads every
    pair the same way."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    klass = "hgvsp_stop_codon_protein_comparison"
    pairs = [
        ("Ter284delinsTer", "Ter284=", "inframe_insertion,stop_retained_variant", "", True),  # shape A
        ("Ter291delinsTer", "Ter291=", "stop_retained_variant", "", True),  # shape A
        ("Ter284delinsIleTer", "Ter284IleextTer22", "stop_lost", "", True),  # shape B
        ("Ter284delinsIleTer", "Ter284IleextTer?", "stop_lost", "", True),  # shape B
        ("Ter284delinsIleLeuValProTer", "Ter284IleextTer?", "stop_lost", "", True),  # shape B
        ("Ter284delinsIleLeuValPro", "Ter284IleextTer9", "stop_lost,inframe_insertion", "", True),  # shape B
        ("Ter284delinsIleLeuValPro", "Ter284IleextTer4", "stop_lost", "", True),  # shape B
        ("Ter284delinsAsxSerGlxTrpArgPheGlyAspValHisTer", "Ter284AsxextTer10", "stop_lost,inframe_insertion", "", True),  # shape C
        ("Ter284delinsAsxTer", "Ter284AsxextTer1", "stop_lost", "", True),  # shape C
        ("Asp283_Ter284insHis", "Ter284HisextTer1", "inframe_insertion", "", True),  # shape D
        ("Asp283_Ter284insHisHisAsxGly", "Ter284HisextTer4", "inframe_insertion", "", True),  # shape D
        ("Ter284=", "Ter284LysextTer47", "stop_lost", "", True),  # shape E
        ("Ter284=", "Ter284LysextTer?", "stop_lost,NMD_transcript_variant", "", True),  # shape E
        ("Trp283_Ter284delinsAsx", "Trp283AsxfsTer4", "stop_lost,protein_altering_variant,NMD_transcript_variant", "", True),  # shape F
        ("Trp283_Ter284delinsAsx", "Trp283AsxfsTer?", "stop_lost,protein_altering_variant", "", True),  # shape F
        ("Trp281_Ter284delinsAsxGlyHis", "Trp281AsxfsTer24", "stop_lost", "", True),  # shape F
        ("Trp281_Ter284delinsAsxGlyHis", "Trp281AsxfsTer9", "stop_lost,inframe_deletion", "", True),  # shape F
        ("Ter557AlaextTer?", "Ter557AlaextTer35", "stop_lost", "", True),  # shape G
        ("Ter591AlaextTer?", "Ter591AlaextTer1", "stop_lost,NMD_transcript_variant", "", True),  # shape G
        ("Ter284delextTer12", "Ter284LysextTer12", "stop_lost,inframe_deletion", "", True),  # shape H
        ("Ter230delextTer11", "Ter230ArgextTer11", "stop_lost,NMD_transcript_variant", "", True),  # shape H
        ("Ter220delextTer?", "Ter220ArgextTer?", "stop_lost,inframe_deletion", "", True),  # shape H
        ("Asp282_Trp283insTrpTer", "Ter284=", "inframe_insertion,stop_retained_variant", "", True),  # shape I
        ("Ter284=", "Ter284delinsTer", "inframe_insertion,stop_retained_variant", "", False),
        ("Ter284IleextTer22", "Ter284delinsIleTer", "stop_lost", "", False),
        ("Ter284AsxextTer10", "Ter284delinsAsxSerGlxTrpArgPheGlyAspValHisTer", "stop_lost,inframe_insertion", "", False),
        ("Ter284HisextTer1", "Asp283_Ter284insHis", "inframe_insertion", "", False),
        ("Ter284LysextTer47", "Ter284=", "stop_lost", "", False),
        ("Ter284delinsTer", "Ter285=", "stop_retained_variant", "", False),
        ("Ter284delinsIleTer", "Ter284=", "stop_lost", "", False),
        ("Ter284delinsIleTer", "Ter284LeuextTer22", "stop_lost", "", False),
        ("Ter284delinsIleTer", "Ter285IleextTer22", "stop_lost", "", False),
        ("Ter284delinsIleLeuValProTer", "Ter284IleextTer3", "stop_lost", "", False),
        ("Ter284delinsAsxSerGlxTrpArgPheGlyAspValHisTer", "Ter284AsxextTer9", "stop_lost", "", False),
        ("Asp283_Ter284insHisHisAsxGly", "Ter284HisextTer5", "inframe_insertion", "", False),
        ("Asp283_Ter284insHisHisAsxGly", "Ter284GlnextTer4", "inframe_insertion", "", False),
        ("Asp282_Ter284insHisHisAsxGly", "Ter284HisextTer4", "inframe_insertion", "", False),
        ("Ter284=", "Ter284LysextTer47", "splice_acceptor_variant,stop_retained_variant", "", False),
        ("Ter284=", "Ter284LysextTer?", "stop_lost", "cds_end_NF", False),
        ("Ter284=", "Ter284LysextTer?", "stop_lost", "cds_start_NF,cds_end_NF", False),
        ("Ter1=", "Ter1LysextTer5", "stop_lost", "cds_start_NF", False),
        ("", "Ter284LysextTer12", "stop_lost", "", False),
        ("Asp283_Ter284delinsPhePro", "Asp283ProfsTer9", "stop_lost", "", False),
        ("Ter284LysextTer12", "Ter284LysextTer13", "stop_lost", "", False),
        ("Ter284delinsIleTerTer", "Ter284IleextTer2", "stop_lost", "", False),
        ("Trp283AsxfsTer4", "Trp283_Ter284delinsAsx", "stop_lost", "", False),
        ("Trp283_Ter284delinsAsx", "Trp282AsxfsTer4", "stop_lost", "", False),
        ("Trp283_Ter284delinsAsx", "Pro283AsxfsTer4", "stop_lost", "", False),
        ("Trp281_Ter284delinsAsxGlyHis", "Trp281GlyfsTer24", "stop_lost", "", False),
        ("Trp281_Ter284delinsAsxGlyHis", "Trp281AsxfsTer3", "stop_lost", "", False),
        ("Trp283_Ter284delinsAsx", "Trp283AsxfsTer4", "inframe_deletion,stop_retained_variant", "", False),
        ("Trp283_Ter284delinsAsxTer", "Trp283AsxfsTer2", "stop_lost", "", False),
        ("Trp283_Ter284delinsAsx", "Trp283AsxfsTer4", "stop_lost", "cds_end_NF", False),
        ("Trp284_Ter284delinsAsx", "Trp284AsxfsTer4", "stop_lost", "", False),
        ("Trp281_Ter284delinsAsxGlyHis", "Trp281AsxfsTer4", "stop_lost,inframe_deletion", "", False),
        ("Trp283_Ter284delinsAsx", "Trp283AsxfsTer2", "stop_lost", "", False),
        ("Ter557AlaextTer35", "Ter557AlaextTer?", "stop_lost", "", False),
        ("Ter557AlaextTer?", "Ter557GlxextTer35", "stop_lost", "", False),
        ("Ter557AlaextTer?", "Ter558AlaextTer35", "stop_lost", "", False),
        ("Ter557AlaextTer?", "Ter557AlaextTer35", "frameshift_variant", "", False),
        ("Ter557AlaextTer?", "Ter557AlaextTer35", "stop_lost", "cds_end_NF", False),
        ("Ter230ArgextTer11", "Ter230delextTer11", "stop_lost", "", False),
        ("Ter230delextTer11", "Ter230ArgextTer12", "stop_lost", "", False),
        ("Ter230delextTer?", "Ter230ArgextTer11", "stop_lost", "", False),
        ("Ter230delextTer11", "Ter231ArgextTer11", "stop_lost", "", False),
        ("Ter230delextTer11", "Ter230ArgextTer11", "inframe_deletion", "", False),
        ("Ter230del", "Ter230ArgextTer11", "stop_lost", "", False),
    ]
    corpus = tmp_path / "c"
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\trs{i}\tAC\tA\t.\tPASS\t.") for i in range(len(pairs))])

    def rows(side: int) -> str:
        text = []
        for i, pair in enumerate(pairs):
            extra = "IMPACT=HIGH;STRAND=1" + (f";FLAGS={pair[3]}" if pair[3] else "")
            if pair[side]:
                extra += f";HGVSp=ENSP1.1:p.{pair[side].replace('=', '%3D')}"
            text.append(_row(f"rs{i}", f"21:{100 * (i + 1) + 1}", "-", "ENST1", pair[2], extra))
        return "".join(text)

    vep, rs = tmp_path / "vep.txt", tmp_path / "rs.txt"
    vep.write_text(rows(0), encoding="utf-8")
    rs.write_text(rows(1), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "HGVSp"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {d["location"]: d["expected_divergence"] for d in m["field_divergences"]}
    assert classes == {
        f"21:{100 * (i + 1) + 1}": klass if pair[4] else "unexplained_residual" for i, pair in enumerate(pairs)
    }
    taken = sum(1 for pair in pairs if pair[4])
    assert m["field_divergence_summary"] == {klass: taken, "unexplained_residual": len(pairs) - taken}
    for reference, candidate, terms, flags, expected in pairs:
        ev = f"ENSP1.1:p.{reference}" if reference else "-"
        assert cvf.stop_codon_protein_comparison(terms, flags, ev, f"ENSP1.1:p.{candidate}") is expected, (reference, candidate)


def test_classify_names_a_residue_and_a_stop_inserted_before_the_last_residue(tmp_path: Path):
    """An `HGVSp` in which Ensembl VEP writes one residue and a stop inserted between two
    neighbouring residues, the residue the second neighbour's, on a row whose terms hold
    `stop_retained_variant`, and vep-rs the stop unchanged at the residue after that neighbour,
    is `hgvsp_stop_codon_protein_comparison`. The sides exchanged, another inserted residue,
    residues that are not neighbours, a flank that is a stop, a stop at another residue, more
    than one residue or none before the inserted stop, a row whose terms hold no
    `stop_retained_variant`, a row whose `FLAGS` name `cds_end_NF`, one value alone inside
    prediction parentheses, and either description against no notation on vep-rs's side stay
    an `unexplained_residual`, and the format comparator's class reads every pair the same
    way. The function reads the terms from the row's key and the `FLAGS` from the reference
    row."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    klass = "hgvsp_stop_codon_protein_comparison"
    inserted, kept = "Asp282_Trp283insTrpTer", "inframe_insertion,stop_retained_variant"
    pairs = [
        (inserted, "Ter284=", kept, "", True),
        (inserted, "Ter284=", "stop_retained_variant,NMD_transcript_variant", "", True),
        ("Met1_Trp2insTrpTer", "Ter3=", kept, "", True),
        ("Ter284=", inserted, kept, "", False),  # the sides exchanged
        ("Asp282_Trp283insAspTer", "Ter284=", kept, "", False),  # another inserted residue
        ("Asp281_Trp283insTrpTer", "Ter284=", kept, "", False),  # residues that are not neighbours
        ("Asp283_Trp283insTrpTer", "Ter284=", kept, "", False),  # one residue number twice
        ("Ter282_Trp283insTrpTer", "Ter284=", kept, "", False),  # a flank that is a stop
        (inserted, "Ter283=", kept, "", False),  # the stop at the second neighbour's residue
        (inserted, "Ter285=", kept, "", False),  # the stop two residues past the second neighbour
        ("Asp282_Trp283insTrpGlyTer", "Ter285=", kept, "", False),  # two residues before the inserted stop
        ("Asp282_Trp283insTer", "Ter283=", kept, "", False),  # no residue before the inserted stop
        (inserted, "Ter284=", "inframe_insertion", "", False),  # no stop_retained_variant
        (inserted, "Ter284=", "inframe_insertion,stop_gained", "", False),  # no stop_retained_variant
        (inserted, "Ter284=", kept, "cds_end_NF", False),  # a coding sequence that no stop codon closes
        (f"({inserted})", "Ter284=", kept, "", False),  # Ensembl VEP's value alone inside prediction parentheses
        (inserted, "(Ter284=)", kept, "", False),  # vep-rs's value alone inside prediction parentheses
        (inserted, "", kept, "", False),  # no notation on vep-rs's side
        ("Ter284=", "", kept, "", False),  # the stop unchanged against no notation
    ]
    corpus = tmp_path / "c"
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\trs{i}\tAC\tA\t.\tPASS\t.") for i in range(len(pairs))])

    def rows(side: int) -> str:
        text = []
        for i, pair in enumerate(pairs):
            extra = "IMPACT=MODERATE;STRAND=-1" + (f";FLAGS={pair[3]}" if pair[3] else "")
            extra += f";HGVSp={value(pair[side]).replace('=', '%3D')}" if pair[side] else ""
            text.append(_row(f"rs{i}", f"21:{100 * (i + 1) + 1}", "-", "ENST1", pair[2], extra))
        return "".join(text)

    def value(description: str) -> str:
        return f"ENSP1.1:p.{description}" if description else ""

    vep, rs = tmp_path / "vep.txt", tmp_path / "rs.txt"
    vep.write_text(rows(0), encoding="utf-8")
    rs.write_text(rows(1), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "HGVSp"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {d["location"]: d["expected_divergence"] for d in m["field_divergences"]}
    assert classes == {
        f"21:{100 * (i + 1) + 1}": klass if pair[4] else "unexplained_residual" for i, pair in enumerate(pairs)
    }
    for reference, candidate, terms, flags, expected in pairs:
        assert cvf.stop_codon_protein_comparison(terms, flags, value(reference), value(candidate)) is expected, (reference, candidate)

    key = ("21:101", "-", "ENST1", "Transcript", kept)
    reference, candidate = f"ENSP1.1:p.{inserted}", "ENSP1.1:p.Ter284%3D"
    assert bgc.stop_codon_protein_comparison(key, reference, candidate, {"FLAGS": ""})
    assert bgc.stop_codon_protein_comparison(key, reference, candidate, None)
    assert not bgc.stop_codon_protein_comparison(key, reference, candidate, {"FLAGS": "cds_start_NF,cds_end_NF"})
    assert not bgc.stop_codon_protein_comparison((*key[:4], "inframe_insertion"), reference, candidate, {"FLAGS": ""})
    assert not bgc.stop_codon_protein_comparison(key, reference, "ENSP2.1:p.Ter284%3D", {"FLAGS": ""})
    for none in ("", "-"):
        assert not bgc.stop_codon_protein_comparison(key, reference, none, {"FLAGS": ""}), none
        assert not bgc.stop_codon_protein_comparison(key, candidate, none, {"FLAGS": ""}), none


def test_classify_names_the_stop_named_on_incomplete_coding_sequence_class(tmp_path: Path):
    """An `HGVSp` in which Ensembl VEP names the stop alone on a row that shows an incomplete
    coding sequence, where vep-rs writes none, is `hgvsp_stop_named_on_incomplete_coding_sequence`; the
    nearest shapes stay an `unexplained_residual`, and the format comparator's class of the
    same name reads every pair the same way."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    klass = "hgvsp_stop_named_on_incomplete_coding_sequence"
    pairs = [
        ("Ter305=", "", "incomplete_terminal_codon_variant,coding_sequence_variant", "cds_end_NF", True),  # the stop unchanged, the flag and the term
        ("Ter305=", "", "incomplete_terminal_codon_variant,coding_sequence_variant", "", True),  # the stop unchanged, the term alone
        ("Ter185del", "", "protein_altering_variant,incomplete_terminal_codon_variant", "cds_end_NF", True),  # the stop deleted
        ("Ter185del", "", "frameshift_variant", "cds_start_NF,cds_end_NF", True),  # the stop deleted, the flag alone
        ("Ter224Val", "", "incomplete_terminal_codon_variant,coding_sequence_variant", "", True),  # the stop replaced by one residue
        ("Ter176HisextTer?", "", "frameshift_variant", "cds_end_NF", True),  # an extension of unknown count
        ("Ter176HisextTer12", "", "frameshift_variant", "cds_end_NF", True),  # an extension with a count
        ("Ter1=", "", "coding_sequence_variant", "cds_start_NF", True),  # the stop unchanged at residue 1
        ("Ter1=", "", "coding_sequence_variant", "cds_start_NF,cds_end_NF", True),  # the stop unchanged at residue 1, both flags
        ("Ter305=", "", "coding_sequence_variant", "", False),  # no flag and no term
        ("Ter1=", "", "coding_sequence_variant", "", False),  # the stop unchanged at residue 1 on an unflagged transcript
        ("Ter5=", "", "coding_sequence_variant", "cds_start_NF", False),  # the start flag at another residue than 1
        ("Ter1Val", "", "coding_sequence_variant", "cds_start_NF", False),  # the start flag on another description at residue 1
        ("Asp223_Ter224del", "", "frameshift_variant", "cds_end_NF", False),  # a range deleted up to the stop
        ("Ter224delinsValHis", "", "protein_altering_variant,incomplete_terminal_codon_variant", "cds_end_NF", False),  # a deletion-insertion of the stop alone
        ("Ter1_Val2insHis", "", "inframe_insertion", "cds_start_NF", False),  # an insertion that names a stop at residue 1
        ("Asp223ValfsTer4", "", "frameshift_variant", "cds_end_NF", False),  # a notation that names no stop alone
        ("Ter305=", "Ter305HisextTer?", "stop_lost", "cds_end_NF", False),  # a candidate that writes a notation
        ("", "Ter305=", "incomplete_terminal_codon_variant,coding_sequence_variant", "cds_end_NF", False),  # the sides exchanged
    ]
    corpus = tmp_path / "c"
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\trs{i}\tAC\tA\t.\tPASS\t.") for i in range(len(pairs))])

    def rows(side: int) -> str:
        text = []
        for i, pair in enumerate(pairs):
            extra = "IMPACT=HIGH;STRAND=1" + (f";FLAGS={pair[3]}" if pair[3] else "")
            if pair[side]:
                extra += f";HGVSp=ENSP1.1:p.{pair[side].replace('=', '%3D')}"
            text.append(_row(f"rs{i}", f"21:{100 * (i + 1) + 1}", "-", "ENST1", pair[2], extra))
        return "".join(text)

    vep, rs = tmp_path / "vep.txt", tmp_path / "rs.txt"
    vep.write_text(rows(0), encoding="utf-8")
    rs.write_text(rows(1), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "HGVSp"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {d["location"]: d["expected_divergence"] for d in m["field_divergences"]}
    assert classes == {
        f"21:{100 * (i + 1) + 1}": klass if pair[4] else "unexplained_residual" for i, pair in enumerate(pairs)
    }
    taken = sum(1 for pair in pairs if pair[4])
    assert m["field_divergence_summary"] == {klass: taken, "unexplained_residual": len(pairs) - taken}
    for reference, candidate, terms, flags, expected in pairs:
        ev = f"ENSP1.1:p.{reference}" if reference else "-"
        av = f"ENSP1.1:p.{candidate}" if candidate else "-"
        assert cvf.stop_named_on_incomplete_coding_sequence(terms, flags, ev, av) is expected, (reference, candidate)
    # Beyond the row: a transcript whose coding sequence the corpus cache holds as ending
    # inside a codon.
    key = ("21:101", "-", "ENST1", "Transcript", "frameshift_variant,stop_lost")
    described = f"ENSP1.1:p.{pairs[5][0]}"
    assert not bgc.stop_named_on_incomplete_coding_sequence(key, described, "-", {"FLAGS": "-"})
    assert bgc.stop_named_on_incomplete_coding_sequence(key, described, "-", {"FLAGS": "-"}, frozenset({"ENST1"}))
    assert not bgc.stop_named_on_incomplete_coding_sequence(key, described, "-", {"FLAGS": "-"}, frozenset({"ENST2"}))


def test_classify_reads_the_stop_named_alone_inside_prediction_parentheses(tmp_path: Path):
    """Under `--hgvsp_use_prediction` Ensembl VEP writes a description inside a pair of
    parentheses. The stop named alone inside one pair, on a row that shows an incomplete coding
    sequence where vep-rs writes none, is `hgvsp_stop_named_on_incomplete_coding_sequence` as it
    is without them. A candidate that writes a notation, more than the stop alone inside the
    pair, a row with neither the flag nor the term and a parenthesis without its pair stay an
    `unexplained_residual`, and the format comparator's class reads every pair the same way."""
    sys.path.insert(0, str(HERE.parent / "concordance"))
    import compare_vep_formats as cvf

    klass = "hgvsp_stop_named_on_incomplete_coding_sequence"
    incomplete = "incomplete_terminal_codon_variant,coding_sequence_variant"
    pairs = [
        ("(Ter305=)", "", incomplete, "cds_end_NF", True),  # the stop unchanged, the flag and the term
        ("(Ter185del)", "", "frameshift_variant", "cds_start_NF,cds_end_NF", True),  # the stop deleted, the flag alone
        ("(Ter224Val)", "", incomplete, "", True),  # the stop replaced by one residue, the term alone
        ("(Ter176HisextTer?)", "", "frameshift_variant", "cds_end_NF", True),  # an extension of unknown count
        ("(Ter1=)", "", "coding_sequence_variant", "cds_start_NF", True),  # the stop unchanged at residue 1
        ("(Ter305=)", "(Ter305HisextTer?)", incomplete, "cds_end_NF", False),  # a candidate that writes a notation
        ("(Asp223_Ter224del)", "", incomplete, "cds_end_NF", False),  # a range deleted up to the stop
        ("(Ter224delinsValHis)", "", incomplete, "cds_end_NF", False),  # a deletion-insertion of the stop alone
        ("(Ter305=)", "", "coding_sequence_variant", "", False),  # no flag and no term
        ("(Ter5=)", "", "coding_sequence_variant", "cds_start_NF", False),  # the start flag at another residue than 1
        ("(Ter305=", "", incomplete, "cds_end_NF", False),  # an opening parenthesis alone
        ("Ter305=)", "", incomplete, "cds_end_NF", False),  # a closing parenthesis alone
        ("((Ter305=))", "", incomplete, "cds_end_NF", False),  # a second pair around the first
    ]
    corpus = tmp_path / "c"
    _manifest(corpus, [_record(f"21\t{100 * (i + 1)}\trs{i}\tAC\tA\t.\tPASS\t.") for i in range(len(pairs))])

    def rows(side: int) -> str:
        text = []
        for i, pair in enumerate(pairs):
            extra = "IMPACT=HIGH;STRAND=1" + (f";FLAGS={pair[3]}" if pair[3] else "")
            if pair[side]:
                extra += f";HGVSp=ENSP1.1:p.{pair[side].replace('=', '%3D')}"
            text.append(_row(f"rs{i}", f"21:{100 * (i + 1) + 1}", "-", "ENST1", pair[2], extra))
        return "".join(text)

    vep, rs = tmp_path / "vep.txt", tmp_path / "rs.txt"
    vep.write_text(rows(0), encoding="utf-8")
    rs.write_text(rows(1), encoding="utf-8")
    assert bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs),
                     "--field", "HGVSp"]) == 0
    m = json.loads((corpus / "manifest.json").read_text())
    classes = {d["location"]: d["expected_divergence"] for d in m["field_divergences"]}
    assert classes == {
        f"21:{100 * (i + 1) + 1}": klass if pair[4] else "unexplained_residual" for i, pair in enumerate(pairs)
    }
    taken = sum(1 for pair in pairs if pair[4])
    assert m["field_divergence_summary"] == {klass: taken, "unexplained_residual": len(pairs) - taken}
    for reference, candidate, terms, flags, expected in pairs:
        ev = f"ENSP1.1:p.{reference}"
        av = f"ENSP1.1:p.{candidate}" if candidate else "-"
        assert cvf.stop_named_on_incomplete_coding_sequence(terms, flags, ev, av) is expected, (reference, candidate)
        key = ("21:101", "-", "ENST1", "Transcript", terms)
        assert bgc.stop_named_on_incomplete_coding_sequence(key, ev, av, {"FLAGS": flags or "-"}) is expected, (reference, candidate)
    # Beyond the row: a transcript whose coding sequence the corpus cache holds as ending
    # inside a codon.
    key = ("21:101", "-", "ENST1", "Transcript", "frameshift_variant,stop_lost")
    described = f"ENSP1.1:p.{pairs[3][0]}"
    assert not bgc.stop_named_on_incomplete_coding_sequence(key, described, "-", {"FLAGS": "-"})
    assert bgc.stop_named_on_incomplete_coding_sequence(key, described, "-", {"FLAGS": "-"}, frozenset({"ENST1"}))
