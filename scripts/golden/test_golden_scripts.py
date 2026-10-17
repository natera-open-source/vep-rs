"""Tests for the golden-corpus generator and the cache pruner on synthetic inputs."""
from __future__ import annotations

import json
import hashlib
import gzip
import random
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_golden_corpus as bgc  # noqa: E402
import extract_reference as xr  # noqa: E402
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
    intron 1 of a transcript filed in the corpus cache)."""
    corpus = tmp_path / "c"
    (corpus / "json_cache" / "transcripts" / "21").mkdir(parents=True)
    (corpus / "json_cache" / "transcripts" / "21" / "33000001-34000000.json").write_text(json.dumps([{
        "stable_id": "ENST1", "start": 33031935, "end": 33041244, "strand": 1, "biotype": "protein_coding",
        "exons": [{"start": 33031935, "end": 33032154}, {"start": 33036103, "end": 33036199}],
        "variation_effect_feature_cache": {"mapper": {"pairs": [
            {"from_start": 1, "from_end": 220, "to_start": 33031935, "to_end": 33032154, "ori": 1},
            {"from_start": 221, "from_end": 317, "to_start": 33036103, "to_end": 33036199, "ori": 1},
        ], "pair_count": 2, "cdna_coding_start": 50, "cdna_coding_end": 300}},
    }]), encoding="utf-8")
    (corpus / "manifest.json").write_text(json.dumps({"records": []}), encoding="utf-8")
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
    rc = bgc.main(["classify", "--corpus", str(corpus), "--vep-default", str(vep), "--vep-rs-output", str(rs)])
    assert rc == 0
    m = json.loads((corpus / "manifest.json").read_text())
    by_feature = {d["feature"]: d["expected_divergence"] for d in m["divergences"]}
    assert by_feature["ENST1"] == "breakend_mate_local_read"
    # a transcript the corpus cache does not carry is not adjudicated
    assert by_feature["ENST2"] == "unexplained_residual"


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


def test_select_preserves_covers_on_regeneration(tmp_path: Path):
    """Regenerating over an existing manifest keeps its `covers` (union with any new `--covers`)
    and its `flag_sets`, so a regenerated corpus keeps certifying the rows it did; a set's
    divergence keys name the old records by ordinal, so they go the way the manifest's own do."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    matrix = _matrix(tmp_path, "symbol", "max_af", "pick")
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    assert bgc.main(["flag-set", str(out), "--name", "pick", "--flag=--pick", "--covers", "pick", "--matrix", str(matrix)]) == 0
    m = json.loads((out / "manifest.json").read_text())
    entry = dict(m["flag_sets"][0])
    m["divergences"] = [{"location": "21:100", "record_indices": [0]}]
    m["flag_sets"][0].update({"divergences": [{"location": "21:100", "record_indices": [0]}], "divergence_summary": {}})
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    assert _select(ref, vcf, out) == 0, "the existing list satisfies the covers requirement"
    m = json.loads((out / "manifest.json").read_text())
    assert m["covers"] == ["symbol"]
    assert m["flag_sets"] == [entry] and "divergences" not in m
    assert _select(ref, vcf, out, "--covers", "max_af", "--covers", "symbol") == 0
    m = json.loads((out / "manifest.json").read_text())
    assert m["covers"] == ["symbol", "max_af"]
    assert list(m).index("covers") < list(m).index("flag_sets") < list(m).index("suites")


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
        ('"covers": ["symbol"], "flag_sets": "pick"', "flag_sets is not a list"),
        ('"covers": ["symbol"], "flag_sets": ["pick"]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"covers": []}]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "default"}]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "-x"}]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a"}, {"name": "a"}]', "names flag set a twice"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "covers": "x"}]', "flag set a covers is not a list of names"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "flags": "--pick"}]', "flag set a flags is not a list of strings"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "input": 5}]', "flag set a input is not a path inside the corpus directory"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "input": "../x.vcf"}]', "flag set a input is not a path inside the corpus directory"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "expected_dir": ""}]', "flag set a expected_dir is not a path inside the corpus directory"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "stdin": "true"}]', "flag set a stdin is not true or false"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "formats": []}]', "flag set a formats is not a non-empty subset of default,tab,vcf,json"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "formats": "vcf"}]', "flag set a formats is not a non-empty subset of default,tab,vcf,json"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "stdin": true, "flags": ["--pick"]}]', "flag set a runs on stdin without --format among its flags"),
    ]
    for body, reason in shapes:
        manifest.write_text('{"records": [], ' + body + "}", encoding="utf-8")
        for argv in commands:
            with pytest.raises(SystemExit) as refusal:
                bgc.main(argv)
            assert str(refusal.value) == f"ERROR: [build_golden_corpus] {manifest} {reason}", argv[0]
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


def test_flag_set_subcommand_writes_an_entry(tmp_path: Path, capsys):
    """An entry carries the whole flag list of its run, its input and stdin mode, the formats
    compared, its expected directory and its covers; a set is rewritten by name; the inputs
    the harness cannot run are refused before anything is written."""
    ref, vcf = _two_record_suite(tmp_path)
    out = tmp_path / "corpus"
    matrix = _matrix(tmp_path, "symbol", "pick", "per_gene")
    assert _select(ref, vcf, out, "--covers", "symbol") == 0
    capsys.readouterr()
    assert bgc.main(["flag-set", str(out), "--name", "pick", "--flag=--pick", "--flag=--symbol", "--covers", "pick", "--matrix", str(matrix)]) == 0
    assert capsys.readouterr().out == "flag set pick: 2 flags, input variants.vcf, formats default,tab,vcf,json, expected expected_pick, covers 1\n"
    m = json.loads((out / "manifest.json").read_text())
    assert m["flag_sets"] == [{
        "name": "pick", "flags": ["--pick", "--symbol"], "input": "variants.vcf", "stdin": False,
        "formats": ["default", "tab", "vcf", "json"], "expected_dir": "expected_pick", "covers": ["pick"],
    }]
    (out / "variants.ensembl").write_text("21\t100\t100\tA/G\t+\n", encoding="utf-8")
    assert bgc.main(["flag-set", str(out), "--name", "ensembl-stdin", "--flag=--format", "--flag=ensembl", "--input", "variants.ensembl",
                     "--stdin", "--format", "default,json", "--expected-dir", "expected_ens", "--matrix", str(matrix)]) == 0
    assert bgc.main(["flag-set", str(out), "--name", "pick", "--flag=--per_gene", "--covers", "per_gene", "--matrix", str(matrix)]) == 0
    m = json.loads((out / "manifest.json").read_text())
    assert [fs["name"] for fs in m["flag_sets"]] == ["pick", "ensembl-stdin"]
    assert m["flag_sets"][0]["flags"] == ["--per_gene"] and m["flag_sets"][0]["covers"] == ["per_gene"]
    assert m["flag_sets"][1] == {
        "name": "ensembl-stdin", "flags": ["--format", "ensembl"], "input": "variants.ensembl", "stdin": True,
        "formats": ["default", "json"], "expected_dir": "expected_ens", "covers": [],
    }
    # `--format=<format>` is the other spelling of the format a stdin run must name.
    assert bgc.main(["flag-set", str(out), "--name", "eq", "--flag=--format=vcf", "--stdin", "--matrix", str(matrix)]) == 0
    capsys.readouterr()
    assert bgc.main(["flag-set", str(out), "--list"]) == 0
    assert capsys.readouterr().out == (
        "pick\t--per_gene\tvariants.vcf\tfalse\tdefault,tab,vcf,json\texpected_pick\tper_gene\n"
        "ensembl-stdin\t--format ensembl\tvariants.ensembl\ttrue\tdefault,json\texpected_ens\t\n"
        "eq\t--format=vcf\tvariants.vcf\ttrue\tdefault,tab,vcf,json\texpected_eq\t\n"
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


def test_info_json_from_perl_info_txt(tmp_path: Path):
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
    assert info["variation_cols"] == ["chr", "variation_name", "failed"]
    assert "cell_types" not in info
    assert info["regulatory"] is True and info["var_type"] == "tabix" and "build" not in info
    assert not (out / "transcripts").exists()


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


def _corpora() -> list[Path]:
    root = Path(__file__).resolve().parents[2] / "tests" / "golden"
    return sorted(p for p in root.glob("*/*") if (p / "manifest.json").is_file())


@pytest.mark.parametrize("corpus", _corpora(), ids=lambda p: f"{p.parent.name}/{p.name}")
def test_committed_expected_files_match_their_provenance_digests(corpus: Path):
    """The expected outputs are stored gzipped; their uncompressed bytes must hash to
    what provenance.json recorded when VEP produced them."""
    prov = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))
    expected_dir = corpus / "expected"
    for name, digest in prov["sha256"].items():
        gz = expected_dir / (name + ".gz")
        plain = expected_dir / name
        data = gzip.open(gz, "rb").read() if gz.is_file() else plain.read_bytes()
        assert hashlib.sha256(data).hexdigest() == digest, f"{corpus.name}/{name} differs from its recorded digest"
    manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
    # The directory is the assembly, optionally suffixed (`GRCh37-hgvs`).
    assembly = corpus.name.split("-", 1)[0]
    assert manifest["release"] == corpus.parent.name and manifest["assembly"] == assembly
    info = json.loads((corpus / "json_cache" / "info.json").read_text(encoding="utf-8"))
    assert info["assembly"] == assembly and info["source_versions"], "info.json must carry the cache source versions"
    if "fasta" in manifest:
        assert (corpus / manifest["fasta"]).is_file() and (corpus / (manifest["fasta"].removesuffix(".gz") + ".fai")).is_file()


def _corpora_with_peptides() -> list[Path]:
    return [p for p in _corpora() if (p / "json_cache" / "transcripts").is_dir()]


@pytest.mark.parametrize("corpus", _corpora_with_peptides(), ids=lambda p: f"{p.parent.name}/{p.name}")
def test_committed_caches_carry_a_peptide_for_every_coding_transcript(corpus: Path):
    """Every committed golden cache carries the peptide beside each CDS, so the format-parity
    tests exercise the sequence-edit overlay the measurement caches exercise."""
    missing = 0
    for shard in (corpus / "json_cache" / "transcripts").glob("*/*"):
        for t in pjc.read_shard(shard):
            v = (t or {}).get("variation_effect_feature_cache") or {}
            if v.get("translateable_seq") and not v.get("peptide"):
                missing += 1
    assert missing == 0, f"{corpus}: {missing} coding transcripts without a peptide"
