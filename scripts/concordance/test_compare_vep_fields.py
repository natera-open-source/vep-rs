"""Tests for compare_vep_fields.py: the field-level comparison on the key intersection."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_vep_fields as cvf  # noqa: E402

HEADER = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n"


def row(uv, loc, allele, gene, feat, ftype, cons, cdna="-", cds="-", prot="-", aa="-", codons="-", existing="-", extra="IMPACT=MODIFIER;STRAND=1"):
    return "\t".join([uv, loc, allele, gene, feat, ftype, cons, cdna, cds, prot, aa, codons, existing, extra]) + "\n"


@pytest.fixture
def pair(tmp_path: Path):
    perl = tmp_path / "perl.txt"
    rust = tmp_path / "rust.txt"
    perl.write_text(
        HEADER
        + row("rs1", "1:100", "T", "ENSG1", "ENST1", "Transcript", "missense_variant", "10", "5", "2", "A/V", "gCc/gTc", "-", "IMPACT=MODERATE;STRAND=1")
        + row("rs2", "1:200", "G", "ENSG1", "ENST1", "Transcript", "intron_variant", extra="IMPACT=MODIFIER;STRAND=1")
        + row("rs3", "1:300", "A", "ENSG2", "ENST2", "Transcript", "upstream_gene_variant", extra="IMPACT=MODIFIER;DISTANCE=1500;STRAND=-1")
        + row("rs4", "1:400", "C", "ENSG2", "ENST2", "Transcript", "synonymous_variant", "40", "35", "12", "L", "ctG/ctA", "-", "IMPACT=LOW;STRAND=-1")
        + row("rs5", "1:500", "T", "ENSG3", "ENST3", "Transcript", "intron_variant,splice_region_variant")
        + row("rs5", "1:500", "T", "ENSG3", "ENST3", "Transcript", "splice_region_variant,intron_variant"),  # duplicate key, different term order
        encoding="utf-8",
    )
    rust.write_text(
        HEADER
        + row("1_100_C/T", "1:100", "T", "ENSG1", "ENST1", "Transcript", "missense_variant", "10", "5", "2", "A/V", "gCc/gTc", "-", "IMPACT=MODERATE;STRAND=1")  # Uploaded_variation differs
        + row("rs2", "1:200", "G", "ENSG1", "ENST1", "Transcript", "intron_variant", extra="IMPACT=MODIFIER;STRAND=1")
        + row("rs3", "1:300", "A", "ENSG2", "ENST2", "Transcript", "upstream_gene_variant", extra="IMPACT=MODIFIER;DISTANCE=1499;STRAND=-1")  # DISTANCE differs
        + row("rs4", "1:400", "C", "ENSG2", "ENST2", "Transcript", "synonymous_variant", "40", "35", "12", "L", "ctG/ctA", "-", "IMPACT=LOW;STRAND=-1;FLAGS=cds_end_NF")  # FLAGS only on rust side
        + row("rs6", "1:600", "G", "ENSG4", "ENST4", "Transcript", "intron_variant"),  # scored-only key
        encoding="utf-8",
    )
    return perl, rust


def test_key_totals_and_duplicates(pair):
    perl, rust = pair
    report = cvf.compare(perl, rust, examples=5, tmp_dir=None)
    t = report["totals"]
    assert t == {
        "reference_keys": 5,
        "scored_keys": 5,
        "intersection": 4,
        "reference_only": 1,  # rs5
        "scored_only": 1,  # rs6
        "reference_duplicate_keys": 1,  # rs5's second row collapses onto the first
        "scored_duplicate_keys": 0,
    }


def test_plain_fields_compare_only_on_intersection(pair):
    perl, rust = pair
    report = cvf.compare(perl, rust, examples=5, tmp_dir=None)
    f = report["fields"]
    assert f["Uploaded_variation"]["compared"] == 4
    assert f["Uploaded_variation"]["differ"] == 1
    assert f["Uploaded_variation"]["top_mismatches"] == [{"reference": "rs1", "scored": "1_100_C/T", "count": 1}]
    for name in ("Gene", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation"):
        assert f[name] == {"compared": 4, "equal": 4, "differ": 0, "agreement": 1.0, "distinct_mismatch_shapes": 0, "top_mismatches": []}, name


def test_extra_keys_compared_individually_with_absence_counted(pair):
    perl, rust = pair
    report = cvf.compare(perl, rust, examples=5, tmp_dir=None)
    e = report["extra"]
    assert e["IMPACT"]["equal"] == 4 and e["IMPACT"]["compared"] == 4
    assert e["STRAND"]["equal"] == 4
    # DISTANCE is present on one intersecting tuple on both sides and differs there.
    assert e["DISTANCE"] == {
        "compared": 1, "equal": 0, "differ": 1, "agreement": 0.0, "distinct_mismatch_shapes": 1,
        "top_mismatches": [{"reference": "1500", "scored": "1499", "count": 1}],
    }
    # FLAGS is emitted by the scored side only, so the reference value is recorded as absent.
    assert e["FLAGS"]["top_mismatches"] == [{"reference": cvf.ABSENT, "scored": "cds_end_NF", "count": 1}]


def test_consequence_order_is_normalised_in_the_key():
    assert cvf.normalize_consequence_set("splice_region_variant,intron_variant") == "intron_variant,splice_region_variant"
    assert cvf.normalize_consequence_set("a,a,b") == "a,b"


def test_cli_writes_report_and_exits_zero(pair, tmp_path: Path):
    perl, rust = pair
    out = tmp_path / "r.json"
    proc = subprocess.run(
        [sys.executable, str(HERE / "compare_vep_fields.py"), "--perl", str(perl), "--rust", str(rust), "--out", str(out)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    report = json.loads(out.read_text())
    assert report["totals"]["intersection"] == 4
    assert "Uploaded_variation" in proc.stdout and "Extra:DISTANCE" in proc.stdout


def test_missing_input_is_a_single_line_error(tmp_path: Path):
    proc = subprocess.run(
        [sys.executable, str(HERE / "compare_vep_fields.py"), "--perl", str(tmp_path / "nope"), "--rust", str(tmp_path / "nope"), "--out", str(tmp_path / "r.json")],
        capture_output=True, text=True,
    )
    assert proc.returncode == 2
    assert proc.stderr.startswith("ERROR: [compare_vep_fields]")
