"""Tests for compare_plugin_outputs.py's end-to-end paths.

The comparator's `--json-output` branch aggregates per-plugin counters that the
per-file results carry; a counter the aggregate does not initialise raises KeyError on
the first file and the run dies after printing a passing report. These tests drive the
script as a subprocess, the way the harness does, on two-record fixtures.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).with_name("compare_plugin_outputs.py")
HEADER = "\t".join(
    [
        "#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
        "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
        "Codons", "Existing_variation", "Extra",
    ]
)


def _row(loc: str, allele: str, feature: str, extra: str) -> str:
    return "\t".join(
        ["v", loc, allele, "ENSG1", feature, "Transcript", "missense_variant", "10", "10",
         "4", "A/T", "gcc/Acc", "-", extra]
    )


def _write(path: Path, rows: list[str]) -> Path:
    path.write_text(HEADER + "\n" + "\n".join(rows) + "\n")
    return path


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True, check=False
    )


@pytest.fixture
def pair(tmp_path: Path) -> tuple[Path, Path]:
    perl = _write(
        tmp_path / "perl.txt",
        [_row("21:1000", "T", "ENST1", "IMPACT=MODERATE;REVEL_score=0.512"),
         _row("21:2000", "G", "ENST2", "IMPACT=MODERATE;REVEL_score=0.100")],
    )
    rust = _write(
        tmp_path / "rust.txt",
        [_row("21:1000", "T", "ENST1", "IMPACT=MODERATE;REVEL_score=0.512"),
         _row("21:2000", "G", "ENST2", "IMPACT=MODERATE;REVEL_score=0.100")],
    )
    return perl, rust


def test_json_output_aggregates_carry_annotated_counts(pair: tuple[Path, Path], tmp_path: Path) -> None:
    perl, rust = pair
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVEL",
                "--json-output", str(out))
    assert proc.returncode == 0, proc.stderr
    data = json.loads(out.read_text())
    agg = data["aggregates"]["REVEL"]
    assert agg["perl_annotated"] == 2 and agg["rust_annotated"] == 2
    assert agg["match"] == 2 and agg["mismatch"] == 0
    assert agg["pass"] is True


def test_json_output_flags_a_mismatch_and_a_rust_only_record(tmp_path: Path) -> None:
    perl = _write(tmp_path / "perl.txt",
                  [_row("21:1000", "T", "ENST1", "REVEL_score=0.512")])
    rust = _write(tmp_path / "rust.txt",
                  [_row("21:1000", "T", "ENST1", "REVEL_score=0.900"),
                   _row("21:3000", "C", "ENST3", "REVEL_score=0.300")])
    out = tmp_path / "result.json"
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVEL",
                "--json-output", str(out))
    assert proc.returncode != 0
    agg = json.loads(out.read_text())["aggregates"]["REVEL"]
    assert agg["mismatch"] == 1 and agg["rust_only"] == 1 and agg["pass"] is False


def test_unknown_plugin_name_is_refused(pair: tuple[Path, Path]) -> None:
    perl, rust = pair
    proc = _run("--perl-output", str(perl), "--rust-output", str(rust), "--plugins", "REVELL")
    assert proc.returncode == 2
    assert "unknown plugin" in proc.stderr
