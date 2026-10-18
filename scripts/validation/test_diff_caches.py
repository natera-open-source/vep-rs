"""Tests for the variation comparison of diff_caches.py.

Run from the repository root:

    python3 -m pytest scripts/validation/test_diff_caches.py
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import diff_caches  # noqa: E402

COLUMNS = [
    "chr",
    "variation_name",
    "failed",
    "somatic",
    "start",
    "end",
    "allele_string",
    "strand",
    "clin_sig",
    "phenotype_or_disease",
    "clin_sig_allele",
    "clin_sig_ref_allele",
    "clinical_impact",
    "pubmed",
    "var_synonyms",
]

PERL_COLUMNS = COLUMNS + ["AF", "AFR", "AMR", "EAS", "EUR", "SAS", "gnomADe", "gnomADg"]
NATIVE_COLUMNS = COLUMNS + ["AF", "AFR", "AMR", "EAS", "EUR", "SAS"]


def write_cache(root: Path, columns: list[str], lines: dict[str, list[str]]) -> Path:
    """Write ``info.json`` with ``columns`` and one ``all_vars.gz`` per chromosome."""
    root.mkdir()
    (root / "info.json").write_text(
        json.dumps({"variation_cols": columns, "var_type": "tabix"})
    )
    for chromosome, records in lines.items():
        directory = root / "variations" / chromosome
        directory.mkdir(parents=True)
        with gzip.open(directory / "all_vars.gz", "wt") as fh:
            fh.write("".join(record + "\n" for record in records))
    return root


def test_parse_all_vars_line_drops_empty_cells():
    record = diff_caches.parse_all_vars_line(
        "21\trs1\t.\t1\t100\t.\tA/G\t.\tpathogenic\t1\t.\t.\t.\t.\t.\n", COLUMNS
    )
    assert record == {
        "chr": "21",
        "variation_name": "rs1",
        "somatic": "1",
        "start": "100",
        "allele_string": "A/G",
        "clin_sig": "pathogenic",
        "phenotype_or_disease": "1",
    }


def test_compare_variations_reports_fixed_fields_and_shared_frequencies():
    baseline = {"start": "100", "allele_string": "A/G", "AF": "G:0.5", "gnomADe": "G:0.1"}
    builder = {"start": "100", "allele_string": "A/G", "AF": "G:0.5"}
    assert diff_caches.compare_variations(baseline, builder, ["AF"]) == {}
    builder["AF"] = "G:0.4"
    assert diff_caches.compare_variations(baseline, builder, ["AF"]) == {
        "frequencies": "content_differs"
    }
    builder["AF"] = "G:0.5"
    builder["allele_string"] = "A/T"
    assert diff_caches.compare_variations(baseline, builder, ["AF"]) == {
        "allele_string": {"baseline": "A/G", "builder": "A/T"}
    }


def test_diff_variations_matches_records_by_chromosome_and_name(tmp_path: Path):
    perl = write_cache(
        tmp_path / "perl",
        PERL_COLUMNS,
        {
            "21": [
                "21\trs1\t.\t.\t100\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.1\t.",
                "21\trs2\t.\t.\t201\t200\t-/T\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.",
                "21\tCOSV1\t.\t1\t300\t.\tCOSMIC_MUTATION\t.\t.\t1\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.",
            ],
            "22": ["22\trs3\t.\t.\t100\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t."],
        },
    )
    native = write_cache(
        tmp_path / "native",
        NATIVE_COLUMNS,
        {
            "21": [
                "21\trs1\t.\t.\t100\t.\tA/G\t.\t.\t.\t.\t.\t.\t.\t.\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.5\tG:0.5",
                "21\trs2\t.\t.\t200\t199\t-/T\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.",
                "21\trs4\t.\t.\t400\t.\tC/T\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.",
            ],
        },
    )
    report = diff_caches.diff_variations(perl, native)
    assert report["total_baseline"] == 4
    assert report["total_builder"] == 3
    assert report["shared"] == 2
    assert report["only_baseline"] == 2
    assert report["only_builder"] == 1
    assert report["matched"] == 1
    assert report["differed"] == 1
    assert report["field_mismatches"] == {"start": 1, "end": 1}
    assert report["frequency_columns_compared"] == ["AF", "AFR", "AMR", "EAS", "EUR", "SAS"]
    assert report["match_rate"] == 0.5


def test_diff_variations_without_variation_files_is_empty(tmp_path: Path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    report = diff_caches.diff_variations(a, b)
    assert report["shared"] == 0
    assert report["match_rate"] == 0.0
    assert report["frequency_columns_compared"] == []
