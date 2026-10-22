"""Tests for the golden-corpus generator and the cache pruner on synthetic inputs."""
from __future__ import annotations

import json
import hashlib
import gzip
import os
import random
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
    """Regenerating over an existing manifest keeps its `covers` (union with any new `--covers`),
    its `masks`, its `reference_warnings` and its `flag_sets`, in that order before `suites`, so a
    regenerated corpus keeps certifying the rows it did and expecting the run-level warnings it did;
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
    m["divergences"] = [{"location": "21:100", "record_indices": [0]}]
    m["flag_sets"][0].update({"divergences": [{"location": "21:100", "record_indices": [0]}], "divergence_summary": {}})
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    assert _select(ref, vcf, out) == 0, "the existing list satisfies the covers requirement"
    m = json.loads((out / "manifest.json").read_text())
    assert m["covers"] == ["symbol"]
    assert m["masks"] == ["regulatory_feature_rows"]
    assert m["reference_warnings"] == [notice]
    assert m["flag_sets"] == [entry] and "divergences" not in m
    assert _select(ref, vcf, out, "--covers", "max_af", "--covers", "symbol") == 0
    m = json.loads((out / "manifest.json").read_text())
    assert m["covers"] == ["symbol", "max_af"]
    assert m["masks"] == ["regulatory_feature_rows"] and m["flag_sets"] == [entry]
    assert m["reference_warnings"] == [notice]
    keys = list(m)
    assert keys.index("covers") < keys.index("masks") < keys.index("reference_warnings") < keys.index("flag_sets") < keys.index("suites")
    m["masks"] = ["splice_rows"]
    (out / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    (out / "variants.vcf").unlink()
    with pytest.raises(SystemExit) as refusal:
        _select(ref, vcf, out)
    assert str(refusal.value) == (f"ERROR: [build_golden_corpus] {out / 'manifest.json'} masks names unknown classes "
                                  "['splice_rows'] (the classes are regulatory_feature_rows and order_only)")
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
         "masks names unknown classes ['splice_rows'] (the classes are regulatory_feature_rows and order_only)"),
        ('"covers": ["symbol"], "reference_warnings": "no FASTA"', "reference_warnings is not a list of texts"),
        ('"covers": ["symbol"], "reference_warnings": ["no FASTA", ""]', "reference_warnings is not a list of texts"),
        ('"covers": ["symbol"], "flag_sets": "pick"', "flag_sets is not a list"),
        ('"covers": ["symbol"], "flag_sets": ["pick"]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"covers": []}]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "default"}]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "-x"}]', "has a flag_sets entry without a name (letters, digits, _ and -, not default)"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a"}, {"name": "a"}]', "names flag set a twice"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "covers": "x"}]', "flag set a covers is not a list of names"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "masks": "order_only"}]', "flag set a masks is not a list of names"),
        ('"covers": ["symbol"], "flag_sets": [{"name": "a", "masks": ["order_only", "splice_rows"]}]',
         "flag set a masks names unknown classes ['splice_rows'] (the classes are regulatory_feature_rows and order_only)"),
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
    # Both mask classes pass at both levels, an empty list on a set (its run masks nothing) too.
    accepted = json.loads(good)
    accepted["masks"] = ["order_only", "regulatory_feature_rows"]
    accepted["flag_sets"] = [{"name": "a", "masks": ["order_only"]}, {"name": "b", "masks": []}]
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
    assert capsys.readouterr().out == ("flag set db: 1 flags, input variants.vcf, reference from the database, "
                                       "formats default,tab,json, expected expected_db, covers 0\n")
    m = json.loads((out / "manifest.json").read_text())
    assert m["flag_sets"][3]["perl_mode"] == "database"
    assert "perl_mode" not in m["flag_sets"][2]
    assert bgc.manifest_shape_error(m) is None
    m["flag_sets"][3]["perl_mode"] = "cache"
    assert bgc.manifest_shape_error(m) == "flag set db perl_mode is not one of offline, database"
    assert bgc.main(["flag-set", str(out), "--list"]) == 0
    assert capsys.readouterr().out == (
        "pick\t--per_gene\tvariants.vcf\tfalse\ttrue\tdefault,tab,vcf,json\texpected_pick\tper_gene\n"
        "ensembl-stdin\t--format ensembl\tvariants.ensembl\ttrue\tfalse\tdefault,json\texpected_ens\t\n"
        "eq\t--format=vcf\tvariants.vcf\ttrue\tfalse\tdefault,tab,vcf,json\texpected_eq\t\n"
        "db\t--ga4gh_vrs\tvariants.vcf\tfalse\tfalse\tdefault,tab,json\texpected_db\t\n"
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
         "flag set domains masks names unknown classes ['splice_rows'] (the classes are regulatory_feature_rows and order_only)"),
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
    # No variation slice was written, so the keys that advertise one are dropped.
    assert "variation_cols" not in info and "var_type" not in info
    assert "cell_types" not in info
    assert info["regulatory"] is True and "build" not in info
    assert not (out / "transcripts").exists()


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
    before anything is written, while a chromosome file with no record is left out, not refused."""
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
    # A chromosome file with no record is skipped, not refused; a file that is not gzip is refused.
    empty = tmp_path / "empty"
    _write_all_vars(empty / "variations" / "21" / "all_vars.gz", [])
    assert pjc.main([*common, "--out", str(tmp_path / "out8"), "--variations", str(empty / "variations"), *native_info]) == 0
    assert not (tmp_path / "out8" / "variations").exists()
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
    `flag_sets`."""
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
    # The directory is the assembly, optionally suffixed (`GRCh37-hgvs`).
    assembly = corpus.name.split("-", 1)[0]
    assert manifest["release"] == corpus.parent.name and manifest["assembly"] == assembly
    cache = corpus / "json_cache"
    if manifest.get("cache_layout") == "ensembl":
        cache = corpus / "cache_root" / "homo_sapiens" / f"{manifest['release']}_{manifest['assembly']}"
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
        cache = corpus / "json_cache"
        if manifest.get("cache_layout") == "ensembl":
            cache = corpus / "cache_root" / "homo_sapiens" / f"{manifest['release']}_{manifest['assembly']}"
        if (cache / "transcripts").is_dir():
            corpora.append(corpus)
    return corpora


@pytest.mark.parametrize("corpus", _corpora_with_peptides(), ids=lambda p: f"{p.parent.name}/{p.name}")
def test_committed_caches_carry_a_peptide_for_every_coding_transcript(corpus: Path):
    """Every committed golden cache carries the peptide beside each CDS, so the format-parity
    tests exercise the sequence-edit overlay the measurement caches exercise."""
    manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
    cache = corpus / "json_cache"
    if manifest.get("cache_layout") == "ensembl":
        cache = corpus / "cache_root" / "homo_sapiens" / f"{manifest['release']}_{manifest['assembly']}"
    missing = 0
    for shard in (cache / "transcripts").glob("*/*"):
        for t in pjc.read_shard(shard):
            v = (t or {}).get("variation_effect_feature_cache") or {}
            if v.get("translateable_seq") and not v.get("peptide"):
                missing += 1
    assert missing == 0, f"{corpus}: {missing} coding transcripts without a peptide"
