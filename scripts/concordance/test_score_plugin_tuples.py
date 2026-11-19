"""Tests for score_plugin_tuples.py: the plugin tuple counts per source, engine and projection.

Every test plants its own registry rows, class table and predicates, so none depends on
which plugins the registry of compare_plugin_outputs.py lists or on what a predicate of
compare_vep_formats.py accepts; the tests of ``plugin_value_follows_stop_codon_term_pair``
and of ``plugin_value_follows_splice_term_past_the_record_span`` score under the comparator's
own classes and predicate, and the tests of the committed class table read it against the
comparator's own classes and the registry's own plugins.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import random
import re
import shutil
import signal
import struct
import subprocess
import sys
import zlib
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_plugin_outputs
import compare_vep_formats as formats
import score_plugin_tuples as spt

SCRIPT = HERE / "score_plugin_tuples.py"
COMMITTED_CLASS_TABLE = HERE / "plugin_class_table.tsv"
# The comparator's own classes and predicate, which ``planted`` replaces in every test.
COMPARATOR_FIELD_CLASSES = formats.FIELD_CLASSES
COMPARATOR_FIELD_CLASS_ACCEPTS = formats.field_class_accepts
PAIRING_TOOL_NAME = "vep-plugin-pairs"
HEADER = spt.DEFAULT_HEADER + "\n"
COUNT_KEYS = ("perl", "rust", "intersection", "adj_perl", "adj_rust", "adj_intersection")
RAW_KEYS = COUNT_KEYS[:3]
PLANTED = {
    "Planted": ["PLANT_score", "PLANT_label"],
    "PlantedDynamic": ["PlantedDynamic_a", "PlantedDynamic_b"],
    "LOEUF": ["LOEUF"],
    "MaxEntScan": ["MaxEntScan_diff"],
    "GWAS": ["GWAS_pmid"],
    "NearestExonJB": ["NearestExonJB"],
    "UTRAnnotator": ["5UTR_consequence"],
    "REVEL": ["REVEL"],
}
PLANTED_CLASSES = (
    "planted_engine_zero",
    "planted_two_numbers",
    "planted_reference_value",
    "planted_row_cdna_one",
    "planted_sibling_value",
)
CLASS_TABLE = (
    "#class\tplugin\tfield names\tmet or not met\tcommit\n"
    "order_only\tNearestExonJB\tNearestExonJB\tmet\t-\n"
    "planted_engine_zero\tLOEUF\tLOEUF\tmet\t-\n"
    "planted_two_numbers\tMaxEntScan\tMaxEntScan_*\tmet\t-\n"
    "planted_reference_value\tGWAS\tGWAS_*\tmet\t-\n"
    "planted_row_cdna_one\tUTRAnnotator\t5UTR_consequence\tmet\t-\n"
    "\n"
    "planted_sibling_value\tREVEL\tREVEL\tnot met\t-\n"
)


def plant_registry(monkeypatch, registry: dict[str, list[str]]) -> None:
    names = [name for name in registry if name not in compare_plugin_outputs.NON_PLUGIN_FIELD_SETS]
    monkeypatch.setattr(compare_plugin_outputs, "PLUGIN_FIELDS", registry)
    monkeypatch.setattr(compare_plugin_outputs, "ALL_PLUGIN_NAMES", names)


def is_number(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


@pytest.fixture(autouse=True)
def planted(monkeypatch, tmp_path_factory):
    """Plants the registry rows, the predicates and a default class table with no row, and
    returns the calls the scorer makes:
    (class, reference entry, field, reference value, engine value, engine entry, record rows)."""
    calls = []

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None):
        calls.append((klass, e, name, expected_value, actual_value, a, record_rows))
        expected = None if formats.blank(expected_value) else expected_value
        actual = None if formats.blank(actual_value) else actual_value
        if klass == "planted_engine_zero":
            return expected is None and actual == "0"
        if klass == "planted_two_numbers":
            return expected is not None and actual is not None and is_number(expected) and is_number(actual)
        if klass == "planted_reference_value":
            return expected is not None
        if klass == "planted_row_cdna_one":
            return e.fields.get("cDNA_position") == "1"
        assert klass == "planted_sibling_value", klass
        return (
            expected is not None
            and record_rows is not None
            and any(
                other.allele != e.allele and other.feature == e.feature and other.fields.get(name) == expected
                for other in record_rows
            )
        )

    def plugin_of_field(field, plugins):
        writers = [name for name, fields in compare_plugin_outputs.PLUGIN_FIELDS.items() if field in fields and name in plugins]
        return writers[0] if len(writers) == 1 else None

    plant_registry(monkeypatch, dict(PLANTED))
    monkeypatch.setattr(compare_plugin_outputs, "plugin_of_field", plugin_of_field, raising=False)
    monkeypatch.setattr(compare_plugin_outputs, "PLUGIN_ALIASES", {}, raising=False)
    monkeypatch.setattr(formats, "CORE_COLUMNS", {"impact"})
    monkeypatch.setattr(formats, "FIELD_CLASSES", dict.fromkeys(PLANTED_CLASSES, ()))
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    # The committed table names the comparator's own classes and the registry's own plugins.
    default_table = tmp_path_factory.mktemp("default") / "classes.tsv"
    default_table.write_text("", encoding="utf-8")
    monkeypatch.setattr(spt, "DEFAULT_CLASS_TABLE", default_table, raising=False)
    return calls


def built_pairing_tool() -> Path | None:
    """The built pairing tool: the file ``VEP_PLUGIN_PAIRS`` names, else the binary under the
    cargo target directory or under the tool's own."""
    named = os.environ.get("VEP_PLUGIN_PAIRS")
    if named:
        return Path(named)
    targets = [HERE.parents[1] / "tools" / PAIRING_TOOL_NAME / "target"]
    if os.environ.get("CARGO_TARGET_DIR"):
        targets.insert(0, Path(os.environ["CARGO_TARGET_DIR"]))
    for target in targets:
        for profile in ("release", "debug"):
            if (target / profile / PAIRING_TOOL_NAME).is_file():
                return target / profile / PAIRING_TOOL_NAME
    return None


# With ``SCORE_PLUGIN_TUPLES_PATH=pairing-tool`` every test that scores through ``score``
# below scores through the pairing tool, so the module holds both paths to one set of
# small outputs; a test of the script's own extraction or sorts runs on its own path alone.
THROUGH_PAIRING_TOOL = os.environ.get("SCORE_PLUGIN_TUPLES_PATH") == "pairing-tool"
own_path_alone = pytest.mark.skipif(THROUGH_PAIRING_TOOL, reason="reads the script's own extraction and sorts")


@pytest.fixture
def pairing_tool() -> Path:
    tool = built_pairing_tool()
    if tool is None:
        pytest.skip(f"no {PAIRING_TOOL_NAME} binary: set VEP_PLUGIN_PAIRS or build tools/{PAIRING_TOOL_NAME}")
    return tool


def registry_plugins() -> list[str]:
    return list(compare_plugin_outputs.ALL_PLUGIN_NAMES)


def row(location: str, allele: str, feature: str, extra: str = "-", cdna: str = "-", uploaded: str = "v") -> str:
    columns = [uploaded, location, allele, "ENSG1", feature, "Transcript", "missense_variant"]
    return "\t".join([*columns, cdna, "-", "-", "-", "-", "-", extra]) + "\n"


def output(path: Path, *rows: str) -> Path:
    path.write_text(HEADER + "".join(rows), encoding="utf-8")
    return path


def tuples(path: Path, *lines: tuple[str, str, str, str, str]) -> Path:
    path.write_text("".join("\t".join(line) + "\n" for line in lines), encoding="utf-8")
    return path


def class_table(tmp_path: Path, content: str | None) -> Path | None:
    if content is None:
        return None
    path = tmp_path / "classes.tsv"
    path.write_text(content, encoding="utf-8")
    return path


def score(
    tmp_path: Path,
    reference: Path,
    engines: dict[str, Path],
    custom: dict[str, str] | None = None,
    residuals: Path | None = None,
    table: str | None = CLASS_TABLE,
    applied: list[str] | None = None,
    populations: dict[str, int] | None = None,
    plugins: list[str] | None = None,
    catalog: Path | None = None,
    run_input: Path | None = None,
    sweep: list[Path] | None = None,
    cache: Path | None = None,
    tool: Path | None = None,
    threads: int | None = None,
    revel: Path | None = None,
    revel_no_match: bool = False,
    spliceai: list[Path] | None = None,
) -> dict:
    if tool is None and THROUGH_PAIRING_TOOL:
        tool = built_pairing_tool()
        assert tool is not None, f"SCORE_PLUGIN_TUPLES_PATH=pairing-tool needs a built {PAIRING_TOOL_NAME}"
    return spt.score(
        reference,
        list(engines.items()),
        registry_plugins() if plugins is None else plugins,
        custom or {},
        tmp_path / "work",
        sort_memory="16M",
        residuals_dir=residuals,
        class_table=class_table(tmp_path, table),
        applied=applied,
        populations=populations,
        gwas_catalog=catalog,
        revel_table=revel,
        revel_no_match=revel_no_match,
        utrannotator_sweep=sweep,
        json_cache=cache,
        run_input=run_input,
        pairing_tool=tool,
        pairing_threads=threads,
        spliceai_table=spliceai,
    )


def residual_lines(directory: Path, source: str, engine: str = "engine", projection: str = "feature") -> list[str]:
    return (directory / f"{source}.{engine}.{projection}.tsv").read_text(encoding="utf-8").splitlines()


def entry(report: dict, source: str, engine: str = "engine", projection: str = "feature") -> dict:
    found = [
        item
        for item in report["scores"]
        if (item["source"], item["engine"], item["projection"]) == (source, engine, projection)
    ]
    assert len(found) == 1, (source, engine, projection, report["scores"])
    return found[0]


def counts(item: dict) -> tuple[int, ...]:
    return tuple(item[key] for key in COUNT_KEYS)


REFERENCE_ROWS = (
    row("1:100", "T", "ENST1", "IMPACT=MODERATE;PLANT_score=0.5;PLANT_label=benign;LOEUF=0.31"),
    row("1:100", "T", "ENST2", "IMPACT=MODERATE;PLANT_score=0.5;MaxEntScan_diff=1.5"),
    row("1:100", "T", "ENST3", "IMPACT=MODIFIER"),
    row("1:200-201", "-", "ENST1", "PLANT_score=0.25;PlantedDynamic_a=x;PlantedDynamic_b=y;Custom=id1;Custom_SIG=p"),
    row("10:5", "G", "-", "PLANT_label=-;GWAS_pmid=17"),
)


def test_reference_against_itself_pairs_every_tuple_of_every_source(tmp_path, planted):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    report = score(tmp_path, reference, {"engine": reference}, {"Custom": "Custom"})
    assert {item["source"] for item in report["scores"]} == {
        "Planted",
        "PlantedDynamic",
        "LOEUF",
        "MaxEntScan",
        "GWAS",
        "Custom",
    }
    for item in report["scores"]:
        assert len(set(counts(item))) == 1, item
        assert not any(item["masked"].values()) and item["not_applied"] == {}
        assert item["text_only_differences"] == 0
        assert all(type(item[key]) is int for key in COUNT_KEYS)
    assert counts(entry(report, "Planted")) == (4, 4, 4, 4, 4, 4)
    assert counts(entry(report, "Planted", projection="allele")) == (3, 3, 3, 3, 3, 3)
    assert counts(entry(report, "PlantedDynamic")) == (2, 2, 2, 2, 2, 2)
    assert counts(entry(report, "Custom")) == (2, 2, 2, 2, 2, 2)
    assert report["sources"]["Custom"] == ["Custom", "Custom_SIG"]
    assert report["not_populated"] == ["NearestExonJB", "REVEL", "UTRAnnotator"]
    assert report["tuple"]["allele"] == ["Location", "Allele", "field", "value"]
    assert planted == []


def test_documented_class_difference_moves_the_raw_counts_alone(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    changed = REFERENCE_ROWS[1].replace("MaxEntScan_diff=1.5", "MaxEntScan_diff=-2.25")
    engine = output(tmp_path / "engine.txt", REFERENCE_ROWS[0], changed, *REFERENCE_ROWS[2:])
    report = score(tmp_path, reference, {"engine": engine})
    for projection in ("feature", "allele"):
        item = entry(report, "MaxEntScan", projection=projection)
        assert counts(item) == (1, 1, 0, 0, 0, 0)
        assert item["masked"] == {"planted_two_numbers": 1}
    assert counts(entry(report, "Planted")) == (4, 4, 4, 4, 4, 4)


def test_documented_class_sets_aside_a_tuple_one_side_alone_writes(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    filled = REFERENCE_ROWS[2].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;LOEUF=0")
    engine = output(tmp_path / "engine.txt", *REFERENCE_ROWS[:2], filled, *REFERENCE_ROWS[3:])
    item = entry(score(tmp_path, reference, {"engine": engine}), "LOEUF")
    assert counts(item) == (1, 2, 1, 1, 1, 1)
    assert item["masked"] == {"planted_engine_zero": 1}


def test_undocumented_difference_moves_the_raw_and_the_adjusted_counts(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    changed = REFERENCE_ROWS[0].replace("PLANT_label=benign", "PLANT_label=pathogenic")
    engine = output(tmp_path / "engine.txt", changed, *REFERENCE_ROWS[1:])
    report = score(tmp_path, reference, {"engine": engine})
    for projection, expected in (("feature", (4, 4, 3, 4, 4, 3)), ("allele", (3, 3, 2, 3, 3, 2))):
        item = entry(report, "Planted", projection=projection)
        assert counts(item) == expected
        assert item["masked"] == {} and item["not_applied"] == {}


def test_class_is_not_applied_to_another_plugins_field(tmp_path, monkeypatch, planted):
    registry = {"LOEUF": ["LOEUF", "LOEUF_b"], "Other": ["LOEUF_twin"]}
    plant_registry(monkeypatch, registry)
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "IMPACT=LOW"))
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "LOEUF=0;LOEUF_b=0;LOEUF_twin=0"))
    residuals = tmp_path / "residuals"
    table = "planted_engine_zero\tLOEUF\tLOEUF,LOEUF_t*\n"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=table)
    for projection, columns in (("feature", "1:100\tT\tENST1"), ("allele", "1:100\tT")):
        own = entry(report, "LOEUF", projection=projection)
        assert counts(own) == (0, 2, 0, 0, 1, 0)
        assert own["masked"] == {"planted_engine_zero": 1}
        assert residual_lines(residuals, "LOEUF", projection=projection) == [f"engine\t{columns}\tLOEUF_b\t0"]
        other = entry(report, "Other", projection=projection)
        assert counts(other) == (0, 1, 0, 0, 1, 0)
        assert other["masked"] == {} and other["not_applied"] == {}
        assert residual_lines(residuals, "Other", projection=projection) == [f"engine\t{columns}\tLOEUF_twin\t0"]
    assert [(call[0], call[2]) for call in planted] == [("planted_engine_zero", "LOEUF")]


def test_field_two_plugins_write_counts_once_for_the_runs_own_plugin(tmp_path, monkeypatch):
    registry = {"dbNSFP": ["SIFT_score", "SIFT_pred"], "TranscriptAnnotator": ["SIFT_score", "Comment"]}
    plant_registry(monkeypatch, registry)
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "SIFT_score=0.1;Comment=kept;IMPACT=LOW"))
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "SIFT_score=0.9;Comment=kept"))
    runs = {
        ("TranscriptAnnotator",): ({"TranscriptAnnotator": (2, 2, 1)}, {}, []),
        ("dbNSFP",): ({"dbNSFP": (1, 1, 0)}, {"Comment": 1}, []),
        ("TranscriptAnnotator", "dbNSFP"): ({"TranscriptAnnotator": (1, 1, 1)}, {"SIFT_score": 1}, ["dbNSFP"]),
    }
    for plugins, (expected, unclaimed, not_populated) in runs.items():
        residuals = tmp_path / "-".join(plugins)
        report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=None, plugins=list(plugins))
        assert report["plugins"] == list(plugins)
        found = {item["source"]: tuple(item[key] for key in RAW_KEYS) for item in report["scores"] if item["projection"] == "feature"}
        assert found == expected
        assert report["not_populated"] == not_populated
        assert report["reference"]["unclaimed_keys"] == unclaimed
        assert report["engines"][0]["unclaimed_keys"] == unclaimed
        assert sorted(path.name.split(".")[0] for path in residuals.iterdir()) == sorted([*expected] * 2)
        if len(plugins) == 1:
            assert residual_lines(residuals, plugins[0]) == [
                "reference\t1:100\tT\tENST1\tSIFT_score\t0.1",
                "engine\t1:100\tT\tENST1\tSIFT_score\t0.9",
            ]


def test_field_is_asked_of_the_registry_among_the_runs_plugins_alone(tmp_path, monkeypatch):
    asked = []

    def plugin_of_field(field, plugins=None):
        asked.append((field, plugins))
        return "Planted" if field.startswith("depth_") and field[6:].isdigit() else None

    monkeypatch.setattr(compare_plugin_outputs, "plugin_of_field", plugin_of_field)
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "depth_10=1;depth_x=2;PLANT_score=3"))
    report = score(tmp_path, reference, {"engine": reference}, table=None, plugins=["Planted", "LOEUF", "Planted"])
    assert report["plugins"] == ["Planted", "LOEUF"]
    assert report["sources"] == {"Planted": ["depth_10"]}
    assert report["reference"]["unclaimed_keys"] == {"PLANT_score": 1, "depth_x": 1}
    assert report["not_populated"] == ["LOEUF"]
    assert sorted(asked) == [("PLANT_score", ("Planted", "LOEUF")), ("depth_10", ("Planted", "LOEUF")), ("depth_x", ("Planted", "LOEUF"))]


def test_run_names_its_sources(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    with pytest.raises(spt.ScorerError, match="--plugins names NoSuchPlugin, no plugin of the registry"):
        score(tmp_path, reference, {"engine": reference}, plugins=["Planted", "NoSuchPlugin"])
    with pytest.raises(spt.ScorerError, match="no source to score"):
        score(tmp_path, reference, {"engine": reference}, plugins=[], table=None)
    report = score(tmp_path, reference, {"engine": reference}, {"Custom": "Custom"}, plugins=[], table=None)
    assert [item["source"] for item in report["scores"]] == ["Custom", "Custom"]
    assert report["plugins"] == [] and report["not_populated"] == []
    assert report["reference"]["unclaimed_keys"] == {
        "GWAS_pmid": 1,
        "LOEUF": 1,
        "MaxEntScan_diff": 1,
        "PLANT_label": 1,
        "PLANT_score": 3,
        "PlantedDynamic_a": 1,
        "PlantedDynamic_b": 1,
    }


def test_custom_source_claims_only_what_no_plugin_of_the_run_claims(tmp_path):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "PLANT_score=1;PLANT_other=2;PLANT_label=x"))
    custom = {"Wide": "PLANT_", "Wider": "PLANT"}
    report = score(tmp_path, reference, {"engine": reference}, custom, plugins=["Planted"])
    assert report["sources"] == {"Planted": ["PLANT_label", "PLANT_score"], "Wide": ["PLANT_other"]}
    assert counts(entry(report, "Planted")) == (2, 2, 2, 2, 2, 2)
    assert counts(entry(report, "Wide")) == (1, 1, 1, 1, 1, 1)
    assert report["not_populated"] == ["Wider"]


def test_keys_no_source_claims_are_listed_per_input(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "IMPACT=LOW;Mystery=1;PLANT_score=1;Bare;Blank=-"),
        row("1:200", "T", "ENST1", "impact=LOW;Mystery=2;Other=x"),
    )
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "PLANT_score=1;Sidecar=7"))
    listed = tuples(
        tmp_path / "engine.tsv",
        ("1:100", "T", "ENST1", "PLANT_score", "1"),
        ("1:100", "T", "ENST1", "IMPACT", "LOW"),
        ("1:100", "T", "ENST1", "Sidecar", "7"),
        ("1:100", "T", "ENST1", "Sidecar", "-"),
    )
    report = score(tmp_path, reference, {"engine": engine, "listed": listed})
    assert report["reference"]["unclaimed_keys"] == {"Mystery": 2, "Other": 1}
    assert report["engines"][0]["unclaimed_keys"] == {"Sidecar": 1}
    assert report["engines"][1]["unclaimed_keys"] == {"Sidecar": 1}
    assert report["engines"][1]["tuples_outside_every_source"] == 2


def test_unclaimed_keys_that_read_alike_keep_every_count(tmp_path):
    reference = tmp_path / "reference.txt"
    written = (HEADER + row("1:100", "T", "ENST1", "KA=1;KB=1;KB=2;PLANT_score=1")).encode()
    reference.write_bytes(written.replace(b"KA", b"K\xff").replace(b"KB", b"K\xfe"))
    report = score(tmp_path, reference, {"engine": reference})
    assert report["reference"]["unclaimed_keys"] == {"K�": 3}
    assert counts(entry(report, "Planted")) == (1, 1, 1, 1, 1, 1)


def test_allele_tuples_find_their_verdicts_whatever_bytes_their_values_carry(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "GWAS_pmid=x"),
        row("1:100", "T", "ENST2", "GWAS_pmid=x\x01"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "IMPACT=LOW"),
        row("1:100", "T", "ENST2", "IMPACT=LOW"),
    )
    report = score(tmp_path, reference, {"engine": engine})
    for projection in ("feature", "allele"):
        item = entry(report, "GWAS", projection=projection)
        assert counts(item) == (2, 0, 0, 0, 0, 0)
        assert item["masked"] == {"planted_reference_value": 2 if projection == "feature" else 1}


def test_missing_and_wrong_values_on_a_one_allele_row_are_residual_under_a_declared_class(tmp_path, planted):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "REVEL=0.51", uploaded="rec1"),
        row("1:200", "T", "ENST1", "REVEL=0.51", uploaded="rec2"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "IMPACT=LOW", uploaded="rec1"),
        row("1:200", "T", "ENST1", "REVEL=0.99", uploaded="rec2"),
    )
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals)
    item = entry(report, "REVEL")
    assert counts(item) == (2, 1, 0, 2, 1, 0)
    assert item["masked"] == {"planted_sibling_value": 0} and item["not_applied"] == {}
    assert item["residual_rows"] == 3
    assert residual_lines(residuals, "REVEL") == [
        "reference\t1:100\tT\tENST1\tREVEL\t0.51",
        "reference\t1:200\tT\tENST1\tREVEL\t0.51",
        "engine\t1:200\tT\tENST1\tREVEL\t0.99",
    ]
    assert [(call[0], call[3], call[4]) for call in planted] == [
        ("planted_sibling_value", "0.51", formats.ABSENT),
        ("planted_sibling_value", "0.51", "0.99"),
    ]
    assert [call[6] for call in planted] == [None, None]


def test_predicate_is_handed_the_reference_rows_of_the_same_record(tmp_path, planted):
    reference = output(
        tmp_path / "reference.txt",
        row("1:300", "A", "ENST1", "REVEL=0.7", uploaded="1_300_C/A/G"),
        row("1:300", "A", "ENST2", "REVEL=0.7", uploaded="1_300_C/A/G"),
        row("1:300", "G", "ENST1", "REVEL=0.7", uploaded="1_300_C/A/G"),
        row("1:300", "C", "ENST1", "REVEL=0.2", uploaded="1_300_G/C"),
        row("1:300", "T", "ENST1", "REVEL=0.7", uploaded="1_300_G/T"),
        row("1:400", "G", "ENST1", "REVEL=0.7", uploaded="1_400_A/G"),
        row("1:500", "A", "ENST1", "REVEL=0.7", uploaded="1_500_C/A/G"),
        row("1:500", "A", "ENST1", "IMPACT=LOW", uploaded="1_500_C/A/G"),
        row("1:500", "G", "ENST1", "REVEL=0.7", uploaded="1_500_C/A/G"),
        row("1:600", "A", "ENST1", "REVEL=0.7", uploaded="rs6"),
        row("1:600", "G", "ENST1", "REVEL=0.7", uploaded="rs6"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:300", "A", "ENST1", "REVEL=0.7", uploaded="other-name"),
        row("1:300", "A", "ENST2", "REVEL=0.7", uploaded="other-name"),
        row("1:300", "G", "ENST1", "IMPACT=LOW", uploaded="other-name"),
        row("1:300", "C", "ENST1", "REVEL=0.2", uploaded="other-name"),
        row("1:300", "T", "ENST1", "IMPACT=LOW", uploaded="other-name"),
        row("1:400", "G", "ENST1", "IMPACT=LOW", uploaded="other-name"),
        row("1:500", "A", "ENST1", "REVEL=0.7", uploaded="other-name"),
        row("1:500", "G", "ENST1", "IMPACT=LOW", uploaded="other-name"),
        row("1:600", "A", "ENST1", "REVEL=0.7", uploaded="other-name"),
        row("1:600", "G", "ENST1", "IMPACT=LOW", uploaded="other-name"),
    )
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals)
    item = entry(report, "REVEL")
    assert counts(item) == (10, 5, 5, 9, 5, 5)
    assert item["masked"] == {"planted_sibling_value": 1}
    assert residual_lines(residuals, "REVEL") == [
        "reference\t1:300\tT\tENST1\tREVEL\t0.7",
        "reference\t1:400\tG\tENST1\tREVEL\t0.7",
        "reference\t1:500\tG\tENST1\tREVEL\t0.7",
        "reference\t1:600\tG\tENST1\tREVEL\t0.7",
    ]
    handed = {(call[1].location, call[1].allele): call[6] for call in planted}
    assert [(other.allele, other.feature) for other in handed["1:300", "G"]] == [("A", "ENST1"), ("G", "ENST1")]
    assert [other.fields["Uploaded_variation"] for other in handed["1:300", "T"]] == ["1_300_G/T"]
    assert [other.allele for other in handed["1:400", "G"]] == ["G"]
    assert handed["1:500", "G"] is None
    assert handed["1:600", "G"] is None
    assert all(call[5].fields["Uploaded_variation"] == "other-name" for call in planted)
    allele = entry(report, "REVEL", projection="allele")
    assert counts(allele) == (9, 4, 4, 8, 4, 4)
    assert allele["masked"] == {"planted_sibling_value": 1}


def test_class_is_not_consulted_without_a_row_that_wrote_the_value_or_none(tmp_path, planted):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "LOEUF=0.31", uploaded="v1"),
        row("1:200", "T", "ENST1", "GWAS_pmid=17", uploaded="v3"),
        row("1:200", "T", "ENST1", "GWAS_pmid=18", uploaded="v4"),
        row("1:300", "T", "ENST1", "LOEUF=", uploaded="v5"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "LOEUF=0.31", uploaded="v1"),
        row("1:100", "T", "ENST1", "LOEUF=0", uploaded="v2"),
        row("1:200", "T", "ENST1", "GWAS_pmid=17", uploaded="v3"),
        row("1:300", "T", "ENST1", "LOEUF=0", uploaded="v5"),
    )
    listed = tuples(
        tmp_path / "engine.tsv",
        ("1:100", "T", "ENST1", "LOEUF", "0.31"),
        ("1:100", "T", "ENST1", "LOEUF", "0"),
        ("1:200", "T", "ENST1", "GWAS_pmid", "17"),
        ("1:300", "T", "ENST1", "LOEUF", "0"),
    )
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine, "listed": listed}, residuals=residuals)
    for ident in ("engine", "listed"):
        zero = entry(report, "LOEUF", engine=ident)
        assert counts(zero) == (1, 3, 1, 1, 2, 1)
        assert zero["masked"] == {"planted_engine_zero": 1}
        assert residual_lines(residuals, "LOEUF", ident) == ["engine\t1:100\tT\tENST1\tLOEUF\t0"]
        left = entry(report, "GWAS", engine=ident)
        assert counts(left) == (2, 1, 1, 2, 1, 1)
        assert left["masked"] == {"planted_reference_value": 0}
        assert residual_lines(residuals, "GWAS", ident) == ["reference\t1:200\tT\tENST1\tGWAS_pmid\t18"]
    asked = [(call[0], call[1].location, call[1].fields["LOEUF"], call[3], call[4]) for call in planted]
    assert asked == [("planted_engine_zero", "1:300", "", formats.ABSENT, "0")] * 2


def test_classes_named_on_the_command_line_are_the_ones_applied(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    changed = REFERENCE_ROWS[1].replace("MaxEntScan_diff=1.5", "MaxEntScan_diff=-2.25")
    filled = REFERENCE_ROWS[2].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;LOEUF=0")
    engine = output(tmp_path / "engine.txt", REFERENCE_ROWS[0], changed, filled, *REFERENCE_ROWS[3:])
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, applied=["planted_two_numbers"])
    assert [item["class"] for item in report["classes"] if item["applied"]] == ["planted_two_numbers"]
    assert report["classes"][1] == {
        "class": "planted_engine_zero",
        "plugins": {"LOEUF": ["LOEUF"]},
        "applied": False,
        "further_columns": ["met", "-"],
    }
    table = tmp_path / "classes.tsv"
    assert report["class_table"] == {"path": str(table), "sha256": hashlib.sha256(table.read_bytes()).hexdigest()}
    for projection, columns in (("feature", "1:100\tT\tENST3"), ("allele", "1:100\tT")):
        applied = entry(report, "MaxEntScan", projection=projection)
        assert counts(applied) == (1, 1, 0, 0, 0, 0)
        assert applied["masked"] == {"planted_two_numbers": 1} and applied["not_applied"] == {}
        held = entry(report, "LOEUF", projection=projection)
        assert tuple(held[key] for key in RAW_KEYS) == (1, 2, 1)
        assert not any(key.startswith("adj_") for key in held)
        assert held["masked"] == {} and held["not_applied"] == {"planted_engine_zero": 1}
        assert held["residual_rows"] == 1
        assert residual_lines(residuals, "LOEUF", projection=projection) == [f"engine\t{columns}\tLOEUF\t0"]
    nothing = score(tmp_path, reference, {"engine": engine}, applied=[])
    assert not any(item["applied"] for item in nothing["classes"])
    assert entry(nothing, "MaxEntScan")["not_applied"] == {"planted_two_numbers": 1}
    with pytest.raises(spt.ScorerError, match="--classes names planted_missing"):
        score(tmp_path, reference, {"engine": engine}, applied=["planted_missing"])


def test_applied_class_is_asked_before_one_that_is_not(tmp_path, monkeypatch):
    plant_registry(monkeypatch, {"GWAS": ["GWAS_pmid"]})
    table = "planted_reference_value\tGWAS\tGWAS_*\nplanted_two_numbers\tGWAS\tGWAS_pmid\n"
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "GWAS_pmid=17"))
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "GWAS_pmid=18"))
    first = entry(score(tmp_path, reference, {"engine": engine}, table=table), "GWAS")
    assert first["masked"] == {"planted_reference_value": 1, "planted_two_numbers": 0}
    second = entry(score(tmp_path, reference, {"engine": engine}, table=table, applied=["planted_two_numbers"]), "GWAS")
    assert second["masked"] == {"planted_two_numbers": 1} and second["not_applied"] == {}
    assert counts(second) == (1, 1, 0, 0, 0, 0)


def test_no_class_table_declares_no_class(tmp_path, planted):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    filled = REFERENCE_ROWS[2].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;LOEUF=0")
    engine = output(tmp_path / "engine.txt", *REFERENCE_ROWS[:2], filled, *REFERENCE_ROWS[3:])
    report = score(tmp_path, reference, {"engine": engine}, table=None)
    assert report["classes"] == [] and report["class_table"] is None
    item = entry(report, "LOEUF")
    assert counts(item) == (1, 2, 1, 1, 2, 1)
    assert item["masked"] == {} and item["not_applied"] == {}
    assert planted == []


@pytest.mark.parametrize(
    ("table", "message"),
    [
        ("planted_engine_zero\tLOEUF\n", "line 1 lacks a class, its plugins and their field names"),
        ("#header\nno_such_class\tLOEUF\tLOEUF\n", "line 2 names no_such_class, no class of"),
        ("planted_engine_zero\tLOEUF\tLOEUF\nplanted_engine_zero\tGWAS\tGWAS_*\n", "line 2 declares planted_engine_zero"),
        ("planted_engine_zero\tLOEUF;GWAS\tLOEUF\n", "line 1 gives 2 plugins and 1 field lists"),
        ("planted_engine_zero\tNoSuchPlugin\tLOEUF\n", "line 1 names 'NoSuchPlugin', no source"),
        ("planted_engine_zero\tLOEUF\t ,\n", "line 1 names 'LOEUF', no source, or no field of it"),
        ("planted_engine_zero\tLOEUF; LOEUF\tLOEUF;LOEUF_b\n", "line 1 names 'LOEUF' twice"),
    ],
)
def test_malformed_class_table_is_refused(tmp_path, table, message):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    with pytest.raises(spt.ScorerError, match=message):
        score(tmp_path, reference, {"engine": reference}, table=table)
    assert not (tmp_path / "work").exists()


def test_class_table_that_is_not_text_is_refused(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    table = tmp_path / "classes.tsv"
    table.write_bytes(b"planted_engine_zero\tLOEUF\tLOEUF\xff\n")
    with pytest.raises(spt.ScorerError, match="classes.tsv: the class table is not UTF-8 text"):
        spt.score(reference, [("engine", reference)], ["LOEUF"], {}, tmp_path / "work", class_table=table)
    assert not (tmp_path / "work").exists()


def test_star_is_the_class_tables_only_wildcard(tmp_path, monkeypatch):
    plant_registry(monkeypatch, {"LOEUF": ["LOEUF", "LOEUF_b", "LOEU?", "L.EUF"]})
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "IMPACT=LOW"))
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "LOEUF=0;LOEUF_b=0;LOEU?=0;L.EUF=0"))
    residuals = tmp_path / "residuals"
    table = "planted_engine_zero\tLOEUF\tLOEU?,[A-Z]OEUF,LO*F_b,L.EUF\n"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=table)
    item = entry(report, "LOEUF")
    assert counts(item) == (0, 4, 0, 0, 1, 0)
    assert item["masked"] == {"planted_engine_zero": 3}
    assert residual_lines(residuals, "LOEUF") == ["engine\t1:100\tT\tENST1\tLOEUF\t0"]
    assert report["classes"][0]["plugins"] == {"LOEUF": ["LOEU?", "[A-Z]OEUF", "LO*F_b", "L.EUF"]}


def test_class_of_two_plugins_reads_each_plugins_own_fields(tmp_path):
    table = "planted_engine_zero\tLOEUF;GWAS\tLOEUF;GWAS_pmid\n"
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "IMPACT=LOW"))
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "LOEUF=0;GWAS_pmid=0"))
    report = score(tmp_path, reference, {"engine": engine}, table=table)
    assert report["classes"] == [
        {
            "class": "planted_engine_zero",
            "plugins": {"LOEUF": ["LOEUF"], "GWAS": ["GWAS_pmid"]},
            "applied": True,
            "further_columns": [],
        }
    ]
    for source in ("LOEUF", "GWAS"):
        assert counts(entry(report, source)) == (0, 1, 0, 0, 0, 0)
        assert entry(report, source)["masked"] == {"planted_engine_zero": 1}


def test_residual_files_are_empty_when_a_reference_scores_against_itself(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    residuals = tmp_path / "kept" / "residuals"
    report = score(tmp_path, reference, {"engine": reference}, {"Custom": "Custom"}, residuals)
    written = sorted(path.name for path in residuals.iterdir())
    assert written == sorted(f"{item['source']}.engine.{item['projection']}.tsv" for item in report["scores"])
    assert len(written) == 12
    assert all(path.read_bytes() == b"" for path in residuals.iterdir())
    assert all(item["residual_rows"] == 0 for item in report["scores"])


def test_residual_file_has_one_line_per_side_for_an_undocumented_difference(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    changed = REFERENCE_ROWS[0].replace("PLANT_label=benign", "PLANT_label=pathogenic")
    lower = REFERENCE_ROWS[3].replace("PLANT_score=0.25", "PLANT_score=0.125")
    engine = output(tmp_path / "engine.txt", changed, *REFERENCE_ROWS[1:3], lower, REFERENCE_ROWS[4])
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals)
    assert residual_lines(residuals, "Planted") == [
        "reference\t1:100\tT\tENST1\tPLANT_label\tbenign",
        "engine\t1:100\tT\tENST1\tPLANT_label\tpathogenic",
        "engine\t1:200-201\t-\tENST1\tPLANT_score\t0.125",
        "reference\t1:200-201\t-\tENST1\tPLANT_score\t0.25",
    ]
    assert residual_lines(residuals, "Planted", projection="allele") == [
        "reference\t1:100\tT\tPLANT_label\tbenign",
        "engine\t1:100\tT\tPLANT_label\tpathogenic",
        "engine\t1:200-201\t-\tPLANT_score\t0.125",
        "reference\t1:200-201\t-\tPLANT_score\t0.25",
    ]
    assert entry(report, "Planted")["residual_rows"] == 4
    assert entry(report, "Planted", projection="allele")["residual_rows"] == 4
    assert residual_lines(residuals, "LOEUF") == []
    assert entry(report, "LOEUF")["residual_rows"] == 0


def test_residual_file_has_no_line_for_a_documented_difference(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    changed = REFERENCE_ROWS[1].replace("MaxEntScan_diff=1.5", "MaxEntScan_diff=-2.25")
    filled = REFERENCE_ROWS[2].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;LOEUF=0")
    engine = output(tmp_path / "engine.txt", REFERENCE_ROWS[0], changed, filled, *REFERENCE_ROWS[3:])
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals)
    for source in ("MaxEntScan", "LOEUF"):
        for projection in ("feature", "allele"):
            item = entry(report, source, projection=projection)
            assert any(item["masked"].values()) and item["intersection"] < item["rust"]
            assert item["residual_rows"] == 0
            assert residual_lines(residuals, source, projection=projection) == []


def test_residuals_need_names_a_file_name_can_carry(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    with pytest.raises(spt.ScorerError, match="cannot be part of a residual file name"):
        score(tmp_path, reference, {"../engine": reference}, residuals=tmp_path / "residuals")
    assert not (tmp_path / "residuals").exists()
    assert counts(entry(score(tmp_path, reference, {"../engine": reference}), "Planted", engine="../engine"))


def test_residual_file_names_may_not_collide(tmp_path):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "X_a=1;Y_a=1"))
    engines = [("B.C", reference), ("C", reference)]
    custom = {"A": "X_", "A.B": "Y_"}
    work = tmp_path / "work"
    with pytest.raises(spt.ScorerError, match="share a residual file name"):
        spt.score(reference, engines, [], custom, work, sort_memory="16M", residuals_dir=tmp_path / "kept")
    with pytest.raises(spt.ScorerError, match="share a residual file name"):
        spt.score(reference, [engines[1]] * 2, [], custom, work, sort_memory="16M", residuals_dir=tmp_path / "kept")
    assert not (tmp_path / "kept").exists()
    assert len(spt.score(reference, engines, [], custom, work, sort_memory="16M")["scores"]) == 8


def test_residuals_directory_of_an_earlier_run_is_refused(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    residuals = tmp_path / "residuals"
    score(tmp_path, reference, {"engine": reference}, residuals=residuals)
    earlier = sorted(path.name for path in residuals.iterdir())
    with pytest.raises(spt.ScorerError, match="holds .tsv files of an earlier run"):
        score(tmp_path, reference, {"later": reference}, residuals=residuals)
    assert sorted(path.name for path in residuals.iterdir()) == earlier


@pytest.mark.skipif(os.geteuid() == 0, reason="the superuser writes in any directory")
def test_residuals_directory_that_cannot_be_written_is_refused_before_scoring(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    residuals = tmp_path / "residuals"
    residuals.mkdir(mode=0o555)
    with pytest.raises(spt.ScorerError, match="cannot write in --residuals"):
        score(tmp_path, reference, {"engine": reference}, residuals=residuals)
    assert not any((tmp_path / "work").iterdir())


def test_items_in_another_order_are_set_aside_as_order_only(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "NearestExonJB=ENSE1|5|start|90,ENSE2|5|end|70"),
        row("1:300", "C", "ENST1", "NearestExonJB=ENSE3|1|start|10"),
        row("1:400", "C", "ENST1", "NearestExonJB=ENSE4|1|start|10"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "NearestExonJB=ENSE2|5|end|70,ENSE1|5|start|90"),
        row("1:300", "C", "ENST1", "NearestExonJB=ENSE3|2|start|10"),
        row("1:400", "C", "ENST1", "IMPACT=LOW"),
    )
    item = entry(score(tmp_path, reference, {"engine": engine}), "NearestExonJB")
    assert counts(item) == (3, 2, 0, 2, 1, 0)
    assert item["masked"] == {"order_only": 1}


def test_field_class_is_consulted_only_where_both_sides_wrote_the_row(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "GWAS_pmid=17"),
        row("1:100", "T", "ENST2", "GWAS_pmid=18"),
    )
    without_the_field = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "IMPACT=LOW"))
    item = entry(score(tmp_path, reference, {"engine": without_the_field}), "GWAS")
    assert counts(item) == (2, 0, 0, 1, 0, 0)
    assert item["masked"] == {"planted_reference_value": 1}


def test_field_class_reads_the_row_that_wrote_the_value(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "5UTR_consequence=first"),
        row("1:100", "T", "ENST1", "5UTR_consequence=second", cdna="1"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "5UTR_consequence=third"),
        row("1:100", "T", "ENST1", "IMPACT=LOW"),
    )
    item = entry(score(tmp_path, reference, {"engine": engine}), "UTRAnnotator")
    assert counts(item) == (2, 1, 0, 1, 1, 0)
    assert item["masked"] == {"planted_row_cdna_one": 1}


def test_row_that_wrote_no_value_is_read_for_one_unpaired_tuple_alone(tmp_path):
    key = ("1:100", "T", "ENST1")
    once = output(tmp_path / "once.txt", row(*key, "IMPACT=LOW"))
    twice = output(tmp_path / "twice.txt", row(*key, "IMPACT=LOW"), row(*key, "IMPACT=LOW"))
    zeros = output(tmp_path / "zeros.txt", row(*key, "LOEUF=0"), row(*key, "LOEUF=0"))
    listed = tuples(tmp_path / "zeros.tsv", (*key, "LOEUF", "0"), (*key, "LOEUF", "0"))
    residuals = tmp_path / "residuals"
    report = score(tmp_path, once, {"engine": zeros, "listed": listed}, residuals=residuals)
    for ident in ("engine", "listed"):
        item = entry(report, "LOEUF", engine=ident)
        assert counts(item) == (0, 2, 0, 0, 1, 0)
        assert item["masked"] == {"planted_engine_zero": 1}
        assert residual_lines(residuals, "LOEUF", ident) == ["engine\t1:100\tT\tENST1\tLOEUF\t0"]
        allele = entry(report, "LOEUF", engine=ident, projection="allele")
        assert counts(allele) == (0, 1, 0, 0, 1, 0)
        assert allele["masked"] == {"planted_engine_zero": 0}
    both = score(tmp_path, twice, {"engine": zeros, "listed": listed})
    for ident in ("engine", "listed"):
        assert counts(entry(both, "LOEUF", engine=ident)) == (0, 2, 0, 0, 0, 0)
    pmids = output(tmp_path / "pmids.txt", row(*key, "GWAS_pmid=17"), row(*key, "GWAS_pmid=17"))
    written = score(tmp_path, pmids, {"once": once, "twice": twice}, residuals=tmp_path / "written")
    left = entry(written, "GWAS", engine="once")
    assert counts(left) == (2, 0, 0, 1, 0, 0)
    assert left["masked"] == {"planted_reference_value": 1}
    assert residual_lines(tmp_path / "written", "GWAS", "once") == ["reference\t1:100\tT\tENST1\tGWAS_pmid\t17"]
    assert counts(entry(written, "GWAS", engine="twice")) == (2, 0, 0, 0, 0, 0)
    distinct = output(tmp_path / "distinct.txt", row(*key, "GWAS_pmid=17"), row(*key, "GWAS_pmid=18"))
    apart = score(tmp_path, distinct, {"once": once, "twice": twice}, residuals=tmp_path / "apart")
    assert residual_lines(tmp_path / "apart", "GWAS", "once") == ["reference\t1:100\tT\tENST1\tGWAS_pmid\t18"]
    assert counts(entry(apart, "GWAS", engine="twice")) == (2, 0, 0, 0, 0, 0)


@pytest.mark.parametrize("positions", [("1", "-"), ("-", "1")])
def test_each_unpaired_copy_of_a_value_is_judged_with_a_row_that_wrote_it(tmp_path, planted, positions):
    key = ("1:100", "T", "ENST1")
    reference = output(
        tmp_path / "reference.txt", *(row(*key, "5UTR_consequence=x", cdna=position) for position in positions)
    )
    engine = output(tmp_path / "engine.txt", row(*key, "IMPACT=LOW"), row(*key, "IMPACT=LOW"))
    residuals = tmp_path / "residuals"
    item = entry(score(tmp_path, reference, {"engine": engine}, residuals=residuals), "UTRAnnotator")
    assert counts(item) == (2, 0, 0, 1, 0, 0)
    assert item["masked"] == {"planted_row_cdna_one": 1}
    assert residual_lines(residuals, "UTRAnnotator") == ["reference\t1:100\tT\tENST1\t5UTR_consequence\tx"]
    assert [call[1].fields["cDNA_position"] for call in planted] == list(positions)


def test_copies_of_a_value_from_two_records_are_each_handed_their_own_records_rows(tmp_path, planted):
    reference = output(
        tmp_path / "reference.txt",
        row("1:300", "A", "ENST1", "REVEL=0.7", uploaded="1_300_C/A/G"),
        row("1:300", "G", "ENST1", "REVEL=0.7", uploaded="1_300_C/A/G"),
        row("1:300", "G", "ENST1", "REVEL=0.7", uploaded="1_300_C/G"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:300", "A", "ENST1", "REVEL=0.7"),
        row("1:300", "G", "ENST1", "IMPACT=LOW"),
        row("1:300", "G", "ENST1", "IMPACT=LOW"),
    )
    residuals = tmp_path / "residuals"
    item = entry(score(tmp_path, reference, {"engine": engine}, residuals=residuals), "REVEL")
    assert counts(item) == (3, 1, 1, 2, 1, 1)
    assert item["masked"] == {"planted_sibling_value": 1}
    assert residual_lines(residuals, "REVEL") == ["reference\t1:300\tG\tENST1\tREVEL\t0.7"]
    handed = [[other.fields["Uploaded_variation"] for other in call[6]] for call in planted]
    assert handed == [["1_300_C/A/G", "1_300_C/A/G"], ["1_300_C/G"]]


def test_values_left_under_a_key_two_records_share_pair_by_record(tmp_path, monkeypatch):
    """Two input records that print one Location and Allele share a key on a Feature. Their
    values left pair by record, the rows' ``Uploaded_variation``, whatever order the bytes
    of the values give: a value with the other side's value of its record, else against its
    record's row that wrote none, never with another record's value. A tuple file states no
    rows, so its values pair in their order."""
    judged = []

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None):
        judged.append((expected_value, actual_value, e.fields["Uploaded_variation"], a.fields.get("Uploaded_variation")))
        return expected_value == actual_value[::-1]

    monkeypatch.setattr(formats, "FIELD_CLASSES", {"planted_reversed_value": ()})
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    table = "planted_reversed_value\tUTRAnnotator\t5UTR_consequence\tmet\t-\norder_only\tNearestExonJB\tNearestExonJB\tmet\t-\n"
    both, one, items = ("1:100", "T", "ENST1"), ("1:200", "T", "ENST1"), ("1:300", "C", "ENST1")
    reference = output(
        tmp_path / "reference.txt",
        row(*both, "5UTR_consequence=ay", uploaded="rs1"),
        row(*both, "5UTR_consequence=bx", uploaded="rs2"),
        row(*one, "IMPACT=LOW", uploaded="rs3"),
        row(*one, "5UTR_consequence=cz", uploaded="rs4"),
        row(*items, "NearestExonJB=ENSE1|1|start|10,ENSE4|4|end|40", uploaded="rs5"),
        row(*items, "NearestExonJB=ENSE2|2|start|20,ENSE3|3|end|30", uploaded="rs6"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row(*both, "5UTR_consequence=ya", uploaded="rs1"),
        row(*both, "5UTR_consequence=xb", uploaded="rs2"),
        row(*one, "5UTR_consequence=zc", uploaded="rs3"),
        row(*one, "IMPACT=LOW", uploaded="rs4"),
        row(*items, "NearestExonJB=ENSE4|4|end|40,ENSE1|1|start|10", uploaded="rs5"),
        row(*items, "NearestExonJB=ENSE3|3|end|30,ENSE2|2|start|20", uploaded="rs6"),
    )
    listed = tuples(tmp_path / "engine.tsv", (*both, "5UTR_consequence", "ya"), (*both, "5UTR_consequence", "xb"))
    residuals = tmp_path / "residuals"
    plugins = ["UTRAnnotator", "NearestExonJB"]
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=table, plugins=plugins)
    item = entry(report, "UTRAnnotator")
    assert counts(item) == (3, 3, 0, 1, 1, 0)
    assert item["masked"] == {"planted_reversed_value": 1}
    assert residual_lines(residuals, "UTRAnnotator") == [
        "reference\t1:200\tT\tENST1\t5UTR_consequence\tcz",
        "engine\t1:200\tT\tENST1\t5UTR_consequence\tzc",
    ]
    assert sorted(judged) == [
        (formats.ABSENT, "zc", "rs3", "rs3"),
        ("ay", "ya", "rs1", "rs1"),
        ("bx", "xb", "rs2", "rs2"),
        ("cz", formats.ABSENT, "rs4", "rs4"),
    ]
    in_order = entry(report, "NearestExonJB")
    assert (counts(in_order), in_order["masked"]) == ((2, 2, 0, 0, 0, 0), {"order_only": 1})
    judged.clear()
    without_rows = entry(score(tmp_path, reference, {"listed": listed}, table=table, plugins=plugins), "UTRAnnotator", "listed")
    assert (counts(without_rows), without_rows["masked"]) == ((3, 2, 0, 3, 2, 0), {"planted_reversed_value": 0})
    assert [pair[:2] for pair in judged] == [("ay", "xb"), ("bx", "ya"), ("cz", formats.ABSENT)]


def test_value_left_stands_against_its_records_row_among_the_rows_that_wrote_none(tmp_path, monkeypatch):
    """Where two rows of a side under a shared key write no value of the field, a value the
    other side has left stands against the row of its own record among them, the second as
    readily as the first. A key that the rows of one side alone share pairs by record as
    well: the value whose record the other side wrote no row for is put to no class."""
    judged = []

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None):
        judged.append((expected_value, actual_value, e.fields["Uploaded_variation"], a.fields.get("Uploaded_variation")))
        return expected_value == actual_value[::-1]

    monkeypatch.setattr(formats, "FIELD_CLASSES", {"planted_reversed_value": ()})
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    table = "planted_reversed_value\tUTRAnnotator\t5UTR_consequence\tmet\t-\n"
    three, one_side = ("1:100", "T", "ENST1"), ("1:200", "T", "ENST1")
    reference = output(
        tmp_path / "reference.txt",
        row(*three, "IMPACT=LOW", uploaded="rs1"),
        row(*three, "IMPACT=LOW", uploaded="rs2"),
        row(*three, "5UTR_consequence=cz", uploaded="rs3"),
        row(*one_side, "5UTR_consequence=ay", uploaded="rs4"),
        row(*one_side, "5UTR_consequence=bx", uploaded="rs5"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row(*three, "IMPACT=LOW", uploaded="rs1"),
        row(*three, "5UTR_consequence=yb", uploaded="rs2"),
        row(*three, "IMPACT=LOW", uploaded="rs3"),
        row(*one_side, "5UTR_consequence=xb", uploaded="rs5"),
    )
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=table, plugins=["UTRAnnotator"])
    item = entry(report, "UTRAnnotator")
    assert counts(item) == (3, 2, 0, 2, 1, 0)
    assert item["masked"] == {"planted_reversed_value": 1}
    assert residual_lines(residuals, "UTRAnnotator") == [
        "reference\t1:100\tT\tENST1\t5UTR_consequence\tcz",
        "engine\t1:100\tT\tENST1\t5UTR_consequence\tyb",
        "reference\t1:200\tT\tENST1\t5UTR_consequence\tay",
    ]
    assert sorted(judged) == [
        (formats.ABSENT, "yb", "rs2", "rs2"),
        ("bx", "xb", "rs5", "rs5"),
        ("cz", formats.ABSENT, "rs3", "rs3"),
    ]


def test_tuple_file_states_no_rows_so_every_reference_row_is_read(tmp_path, planted):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "GWAS_pmid=17"),
        row("1:100", "T", "ENST2", "GWAS_pmid=18"),
    )
    engine = tuples(
        tmp_path / "engine.tsv",
        ("1:100", "T", "ENST2", "PLANT_label", "kept"),
        ("1:900", "A", "ENST9", "GWAS_pmid", "3"),
    )
    item = entry(score(tmp_path, reference, {"engine": engine}), "GWAS")
    assert counts(item) == (2, 1, 0, 0, 1, 0)
    assert item["masked"] == {"planted_reference_value": 2}
    assert [call[5].fields for call in planted] == [
        {"Location": "1:100", "Allele": "T", "Feature": "ENST1"},
        {"Location": "1:100", "Allele": "T", "Feature": "ENST2", "PLANT_label": "kept"},
    ]


@pytest.mark.parametrize(
    ("one", "other", "agree"),
    [
        (b"1.5", b"1.500001", True),
        (b"1.5", b"1.5000011", False),
        (b"1.5", b"1.499999", True),
        (b"1.5", b"1.4999989", False),
        (b"0.5", b"0.5000005", True),
        (b"0.5", b"0.500001", False),
        (b"0.5", b"0.4999995", True),
        (b"0.5", b"0.49999949", False),
        (b"-0.5", b"-0.5000005", True),
        (b"-0.5", b"0.5", False),
        (b"4E-24", b"4.000001E-24", True),
        (b"4E-24", b"4.00001E-24", False),
        (b"-0.000001", b"0", False),
        (b"0", b"1e-30", False),
        (b"0", b"0.0", True),
        (b"-0", b"0", True),
        (b"1000000", b"1000001", False),
        (b"1000000", b"1000000.000001", True),
        (b"1", b"2", False),
        (b"100", b"1e2", True),
        (b"1e-3", b"0.001", True),
        (b"5E-8", b"0.00000005", True),
        (b"1e999", b"1e999", True),
        (b"1e999", b"-1e999", False),
        (b"12345678901234567890.0000001", b"12345678901234567890.0000012", False),
        (b"12345678901234567890.0000001", b"12345678901234567890.0000011", True),
        (b"0.0000010000000000000000000000000001", b"0", False),
        (b"1.", b"1", True),
        (b".5", b"0.5", True),
        (b"1e1000", b"1e1000", False),
        (b"nan", b"nan", False),
        (b"inf", b"inf", False),
        (b"1_0", b"10", False),
        (b" 1", b"1", False),
        (b"0x10", b"16", False),
        (b"1,2", b"1,2.0", False),
        (b"abc", b"abc", False),
    ],
)
def test_numbers_agree_within_the_tolerance_in_decimal_arithmetic(one, other, agree):
    assert spt.numbers_agree(one, other) is agree
    assert spt.numbers_agree(other, one) is agree


def test_tolerance_is_the_published_one():
    assert Decimal("0.000001") == spt.TOLERANCE


SMALL_NUMBERS_APART = (
    ("3E-24", "6E-9"),
    ("7E-9", "8E-10"),
    ("3E-11", "4E-18"),
    ("3E-18", "3E-8"),
    ("4E-7", "7E-7"),
    ("2E-16", "4E-8"),
    ("3E-9", "9E-7"),
)
WHOLE_NUMBERS_ONE_APART = (("1", "2"), ("0", "1"), ("1000000", "1000001"), ("250", "251"), ("-7", "-8"))
SPELLINGS_OF_ONE_NUMBER = (("5E-9", "5e-09"), ("5E-9", "0.000000005"), ("4E-24", "4.0E-24"), ("1.0", "1"), ("100", "1e2"), ("0.10", "0.1"))


def scored_pairs(tmp_path: Path, pairs: tuple[tuple[str, str], ...]) -> tuple[dict, list[str]]:
    """The Planted entry and the residual lines of a run whose rows each write one pair of
    values, the reference's and the engine's, under one key."""
    reference = output(tmp_path / "reference.txt", *(row(f"1:{at}", "T", "ENST1", f"PLANT_score={one}") for at, (one, _) in enumerate(pairs, 100)))
    engine = output(tmp_path / "engine.txt", *(row(f"1:{at}", "T", "ENST1", f"PLANT_score={other}") for at, (_, other) in enumerate(pairs, 100)))
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=None)
    assert report["value_rule"] == {
        "numeric_tolerance": "0.000001",
        "numeric_tolerance_applies": ["to the difference as a share of the larger number", "to the difference"],
    }
    return entry(report, "Planted"), residual_lines(residuals, "Planted")


@pytest.mark.parametrize("pairs", [SMALL_NUMBERS_APART, WHOLE_NUMBERS_ONE_APART], ids=["two small numbers", "whole numbers one apart"])
def test_numbers_close_in_size_alone_or_one_apart_are_not_equal(tmp_path, pairs):
    for one, other in pairs:
        assert not spt.numbers_agree(one.encode(), other.encode()), (one, other)
    item, lines = scored_pairs(tmp_path, pairs)
    assert counts(item) == (len(pairs), len(pairs), 0, len(pairs), len(pairs), 0)
    assert (item["text_only_differences"], item["residual_rows"]) == (0, 2 * len(pairs))
    assert len(lines) == 2 * len(pairs)


def test_two_spellings_of_one_number_are_a_text_only_pair(tmp_path):
    for one, other in SPELLINGS_OF_ONE_NUMBER:
        assert spt.numbers_agree(one.encode(), other.encode()), (one, other)
    item, lines = scored_pairs(tmp_path, SPELLINGS_OF_ONE_NUMBER)
    pairs = len(SPELLINGS_OF_ONE_NUMBER)
    assert counts(item) == (pairs, pairs, pairs, pairs, pairs, pairs)
    assert (item["text_only_differences"], item["residual_rows"], lines) == (pairs, 0, [])


def test_value_rule_at_its_boundary_and_text_only_differences(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "PLANT_score=1.5"),
        row("1:200", "T", "ENST1", "PLANT_score=1.5"),
        row("1:300", "T", "ENST1", "PLANT_score=0.50"),
        row("1:400", "T", "ENST1", "PLANT_score=0.5;PLANT_label=A"),
        row("1:500", "T", "ENST1", "PLANT_score=0.5"),
        row("1:600", "T", "ENST1", "PLANT_score=0.5"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "PLANT_score=1.500001"),
        row("1:200", "T", "ENST1", "PLANT_score=1.5000011"),
        row("1:300", "T", "ENST1", "PLANT_score=0.5"),
        row("1:400", "T", "ENST1", "PLANT_score=0.5;PLANT_label=a"),
        row("1:500", "T", "ENST1", "PLANT_score=0.5000005"),
        row("1:600", "T", "ENST1", "PLANT_score=0.500001"),
    )
    report = score(tmp_path, reference, {"engine": engine, "same": reference})
    item = entry(report, "Planted")
    assert counts(item) == (7, 7, 4, 7, 7, 4)
    assert item["text_only_differences"] == 3
    assert entry(report, "Planted", engine="same")["text_only_differences"] == 0
    allele = entry(report, "Planted", projection="allele")
    assert counts(allele) == (7, 7, 4, 7, 7, 4)
    assert allele["text_only_differences"] == 3


def test_values_under_one_key_pair_one_to_one(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "PLANT_score=0.1234561"),
        row("1:100", "T", "ENST1", "PLANT_score=0.1234562"),
    )
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "PLANT_score=0.12345615"))
    item = entry(score(tmp_path, reference, {"engine": engine}), "Planted")
    assert counts(item) == (2, 1, 1, 2, 1, 1)
    assert item["text_only_differences"] == 1


def test_numbers_pair_as_many_as_the_tolerance_allows():
    assert spt.pair_values([b"1.5", b"1.500001"], [b"1.5000005", b"14.999992e-1"]) == (0, 2, [], [])
    assert spt.pair_values([b"0.5", b"0.500001"], [b"0.5000005", b"4.999992e-1"]) == (0, 1, [b"0.500001"], [b"4.999992e-1"])
    assert spt.pair_values([b"a", b"1", b"3"], [b"3.0000001", b"b", b"0.9999995", b"7"]) == (0, 2, [b"a"], [b"b", b"7"])
    assert spt.pair_values([b"1", b"x"], [b"1", b"1.0", b"x"]) == (2, 0, [], [b"1.0"])
    assert spt.pair_values([b"1", b"3"], [b"2.9999999"]) == (0, 1, [b"1"], [])
    assert spt.pair_values([b"2.9999999"], [b"1", b"3"]) == (0, 1, [], [b"1"])
    assert spt.pair_values([b"1.0", b"1.0"], [b"1", b"1.00"]) == (0, 2, [], [])


def test_value_held_more_than_once_pairs_as_often_as_both_sides_hold_it():
    assert spt.pair_values([b"a", b"a", b"b"], [b"a", b"b", b"b"]) == (2, 0, [b"a"], [b"b"])
    assert spt.pair_values([b"1.0", b"1.0"], [b"1"]) == (0, 1, [b"1.0"], [])
    assert spt.pair_values([b"1", b"1"], [b"1", b"1.0"]) == (1, 1, [], [])
    assert spt.pair_values([b"x", b"1", b"x", b"1"], [b"1.00", b"x"]) == (1, 1, [b"x", b"1"], [])
    assert spt.pair_values([b"a", b"a"], []) == (0, 0, [b"a", b"a"], [])


def test_tuple_file_engine(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    engine = tuples(
        tmp_path / "engine.tsv",
        ("1:100", "T", "ENST1", "PLANT_score", "0.5"),
        ("1:100", "T", "ENST1", "PLANT_label", "benign"),
        ("1:100", "T", "ENST2", "PLANT_score", "0.50"),
        ("1:200-201", "-", "ENST1", "PLANT_score", "0.25"),
        ("1:200-201", "-", "ENST1", "PLANT_label", "-"),
        ("1:200-201", "-", "ENST1", "IMPACT", "MODIFIER"),
        ("1:200-201", "-", "ENST1", "PlantedDynamic_a", "x"),
    )
    report = score(tmp_path, reference, {"engine": engine})
    item = entry(report, "Planted")
    assert counts(item) == (4, 4, 4, 4, 4, 4)
    assert item["text_only_differences"] == 1
    assert counts(entry(report, "PlantedDynamic")) == (2, 1, 1, 2, 1, 1)
    described = report["engines"][0]
    assert described["format"] == "tuples"
    assert described["tuples"] == 7
    assert described["tuples_outside_every_source"] == 1
    assert described["sha256"] == hashlib.sha256(engine.read_bytes()).hexdigest()


def test_allele_projection_drops_the_feature_and_deduplicates(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "PLANT_score=0.9"),
        row("1:100", "T", "ENST2", "PLANT_score=0.9"),
        row("1:100", "T", "ENST3", "IMPACT=LOW"),
        row("1:100", "G", "ENST1", "PLANT_score=0.2"),
    )
    engine = tuples(
        tmp_path / "engine.tsv",
        ("1:100", "T", "ENST1", "PLANT_score", "0.9"),
        ("1:100", "T", "ENST3", "PLANT_score", "0.9"),
        ("1:100", "T", "OWN_MODEL.1", "PLANT_score", "0.9"),
        ("1:100", "G", "OWN_MODEL.1", "PLANT_score", "0.2"),
        ("1:100", "C", "OWN_MODEL.1", "PLANT_score", "0.7"),
    )
    report = score(tmp_path, reference, {"engine": engine})
    assert counts(entry(report, "Planted")) == (3, 5, 1, 3, 5, 1)
    assert counts(entry(report, "Planted", projection="allele")) == (2, 3, 2, 2, 3, 2)


def test_allele_tuple_is_set_aside_when_every_unpaired_feature_tuple_was(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "IMPACT=LOW"),
        row("1:100", "T", "ENST2", "IMPACT=LOW"),
        row("1:200", "C", "ENST1", "IMPACT=LOW"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "LOEUF=0"),
        row("1:100", "T", "ENST2", "LOEUF=0"),
        row("1:200", "C", "ENST1", "LOEUF=0"),
        row("1:200", "C", "ENST7", "LOEUF=0"),
    )
    report = score(tmp_path, reference, {"engine": engine})
    feature = entry(report, "LOEUF")
    assert counts(feature) == (0, 4, 0, 0, 1, 0)
    assert feature["masked"] == {"planted_engine_zero": 3}
    allele = entry(report, "LOEUF", projection="allele")
    assert counts(allele) == (0, 2, 0, 0, 1, 0)
    assert allele["masked"] == {"planted_engine_zero": 1}
    assert report["not_populated"] == sorted(PLANTED)


def test_allele_tuple_under_an_applied_and_a_held_class_rests_on_the_held_one(tmp_path, monkeypatch):
    plant_registry(monkeypatch, {"GWAS": ["GWAS_pmid"]})
    table = "planted_row_cdna_one\tGWAS\tGWAS_pmid\nplanted_reference_value\tGWAS\tGWAS_pmid\n"
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "GWAS_pmid=17", cdna="1"),
        row("1:100", "T", "ENST2", "GWAS_pmid=17"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "IMPACT=LOW", cdna="1"),
        row("1:100", "T", "ENST2", "IMPACT=LOW"),
    )
    report = score(tmp_path, reference, {"engine": engine}, table=table, applied=["planted_row_cdna_one"])
    feature = entry(report, "GWAS")
    assert tuple(feature[key] for key in RAW_KEYS) == (2, 0, 0)
    assert feature["masked"] == {"planted_row_cdna_one": 1}
    assert feature["not_applied"] == {"planted_reference_value": 1}
    assert feature["residual_rows"] == 1
    allele = entry(report, "GWAS", projection="allele")
    assert tuple(allele[key] for key in RAW_KEYS) == (1, 0, 0)
    assert allele["masked"] == {"planted_row_cdna_one": 0}
    assert allele["not_applied"] == {"planted_reference_value": 1}
    assert allele["residual_rows"] == 1
    assert not any(key.startswith("adj_") for key in (*feature, *allele))


def test_several_unpaired_values_of_one_allele_each_find_their_class(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "MaxEntScan_diff=1.5"),
        row("1:100", "T", "ENST2", "MaxEntScan_diff=2.5"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "MaxEntScan_diff=-1.25"),
        row("1:100", "T", "ENST2", "MaxEntScan_diff=-2.25"),
    )
    report = score(tmp_path, reference, {"engine": engine})
    for projection in ("feature", "allele"):
        item = entry(report, "MaxEntScan", projection=projection)
        assert counts(item) == (2, 2, 0, 0, 0, 0)
        assert item["residual_rows"] == 0


def test_allele_tuple_whose_feature_tuples_all_paired_stays_in_both_counts(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "LOEUF=0"),
        row("1:100", "T", "ENST2", "IMPACT=LOW"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "LOEUF=0.0"),
        row("1:100", "T", "ENST2", "LOEUF=0"),
    )
    report = score(tmp_path, reference, {"engine": engine})
    feature = entry(report, "LOEUF")
    assert counts(feature) == (1, 2, 1, 1, 1, 1)
    assert feature["masked"] == {"planted_engine_zero": 1}
    assert feature["text_only_differences"] == 1
    allele = entry(report, "LOEUF", projection="allele")
    assert counts(allele) == (1, 2, 1, 1, 2, 1)
    assert allele["masked"] == {"planted_engine_zero": 0}
    assert allele["text_only_differences"] == 0


def test_custom_source_groups_the_fields_under_its_prefix(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "ClinVar=7105;ClinVar_CLNSIG=Benign;ClinVarious=1;gnomAD_AF=0.01"),
    )
    engine = tuples(
        tmp_path / "engine.tsv",
        ("1:100", "T", "ENST1", "ClinVar", "7105"),
        ("1:100", "T", "ENST1", "ClinVar_CLNSIG", "Pathogenic"),
        ("1:100", "T", "ENST1", "gnomAD_AF", "0.0100"),
    )
    report = score(tmp_path, reference, {"engine": engine}, {"ClinVar": "ClinVar", "gnomAD": "gnomAD_"})
    assert counts(entry(report, "ClinVar")) == (3, 2, 1, 3, 2, 1)
    assert counts(entry(report, "gnomAD")) == (1, 1, 1, 1, 1, 1)
    assert entry(report, "gnomAD")["text_only_differences"] == 1
    assert report["sources"] == {"ClinVar": ["ClinVar", "ClinVar_CLNSIG", "ClinVarious"], "gnomAD": ["gnomAD_AF"]}


def test_custom_source_may_not_take_a_registry_plugins_name():
    with pytest.raises(spt.ScorerError, match="--custom-source Planted names a plugin of the registry"):
        spt.Sources([], {"Planted": "PLANT_"})


def test_non_plugin_field_set_of_the_registry_is_no_source(tmp_path, monkeypatch):
    non_plugin = next(iter(compare_plugin_outputs.NON_PLUGIN_FIELD_SETS))
    plant_registry(monkeypatch, {non_plugin: ["HGVSc"], "Planted": ["PLANT_score"]})
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "HGVSc=ENST1:c.1A>T;PLANT_score=1"))
    report = score(tmp_path, reference, {"engine": reference}, table=None)
    assert report["plugins"] == ["Planted"]
    assert {item["source"] for item in report["scores"]} == {"Planted"}
    assert report["not_populated"] == []
    with pytest.raises(spt.ScorerError, match=f"--plugins names {non_plugin}, no plugin of the registry"):
        score(tmp_path, reference, {"engine": reference}, table=None, plugins=[non_plugin])


def test_source_only_the_engine_writes_is_scored_and_listed_as_not_populated(tmp_path):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "PLANT_score=1"))
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "PLANT_score=1;PlantedDynamic_a=x"))
    report = score(tmp_path, reference, {"engine": engine}, {"Unused": "Unused_"})
    assert counts(entry(report, "PlantedDynamic")) == (0, 1, 0, 0, 1, 0)
    assert report["not_populated"] == sorted({*PLANTED, "Unused"} - {"Planted"})


def test_gzip_input_and_the_digest_of_the_file_as_stored(tmp_path):
    plain = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    packed = tmp_path / "reference.txt.gz"
    packed.write_bytes(gzip.compress(plain.read_bytes()))
    report = score(tmp_path, packed, {"engine": plain})
    assert report["reference"]["sha256"] == hashlib.sha256(packed.read_bytes()).hexdigest()
    assert report["engines"][0]["sha256"] == hashlib.sha256(plain.read_bytes()).hexdigest()
    assert report["reference"]["rows"] == report["engines"][0]["rows"] == len(REFERENCE_ROWS)
    assert counts(entry(report, "Planted")) == (4, 4, 4, 4, 4, 4)


@pytest.mark.parametrize("packed", [False, True])
def test_reader_digest_covers_the_file_whenever_it_is_taken(tmp_path, packed):
    rng = random.Random(116)
    data = "".join(f"1:{n}\tT\tENST1\tPLANT_label\t{rng.randbytes(50).hex()}\n" for n in range(40_000)).encode()
    path = tmp_path / "input"
    path.write_bytes(gzip.compress(data) if packed else data)
    assert path.stat().st_size > 1 << 20
    reader = spt.Reader(path)
    try:
        assert next(iter(reader.lines)) == data[: data.index(b"\n") + 1]
        assert reader.sha256() == hashlib.sha256(path.read_bytes()).hexdigest()
    finally:
        reader.close()


def test_rows_that_share_a_key_are_never_collapsed(tmp_path):
    twice = output(tmp_path / "twice.txt", *REFERENCE_ROWS, *REFERENCE_ROWS)
    once = output(tmp_path / "once.txt", *REFERENCE_ROWS)
    residuals = tmp_path / "residuals"
    report = score(tmp_path, twice, {"once": once, "twice": twice}, residuals=residuals)
    assert report["reference"]["rows"] == 2 * len(REFERENCE_ROWS)
    single = entry(report, "Planted", engine="once")
    assert counts(single) == (8, 4, 4, 8, 4, 4)
    assert (single["residual_rows"], single["shared_keys"]) == (4, 3)
    assert residual_lines(residuals, "Planted", "once") == [
        "reference\t1:100\tT\tENST1\tPLANT_label\tbenign",
        "reference\t1:100\tT\tENST1\tPLANT_score\t0.5",
        "reference\t1:100\tT\tENST2\tPLANT_score\t0.5",
        "reference\t1:200-201\t-\tENST1\tPLANT_score\t0.25",
    ]
    double = entry(report, "Planted", engine="twice")
    assert counts(double) == (8, 8, 8, 8, 8, 8)
    assert (double["residual_rows"], double["shared_keys"]) == (0, 3)
    for ident in ("once", "twice"):
        allele = entry(report, "Planted", engine=ident, projection="allele")
        assert counts(allele) == (3, 3, 3, 3, 3, 3)
        assert allele["shared_keys"] == 0
    alone = entry(score(tmp_path, once, {"engine": once}), "Planted")
    assert (counts(alone), alone["shared_keys"]) == ((4, 4, 4, 4, 4, 4), 0)


def test_values_under_a_shared_key_pair_as_often_as_each_side_writes_them(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=a", uploaded="v1"),
        row("1:100", "T", "ENST1", "PLANT_score=0.7;PLANT_label=a", uploaded="v2"),
        row("1:100", "T", "ENST1", "PLANT_score=0.7", uploaded="v3"),
        row("1:100", "T", "ENST2", "PLANT_score=1.0", uploaded="v1"),
        row("1:100", "T", "ENST2", "PLANT_score=1.0", uploaded="v2"),
        row("1:900", "G", "ENST1", "PLANT_score=1"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "PLANT_score=0.7;PLANT_label=a", uploaded="v2"),
        row("1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=a", uploaded="v1"),
        row("1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=a", uploaded="v3"),
        row("1:100", "T", "ENST2", "PLANT_score=1", uploaded="v1"),
        row("1:100", "T", "ENST2", "PLANT_score=1", uploaded="v2"),
        row("1:900", "G", "ENST1", "PLANT_score=1"),
    )
    listed = tuples(
        tmp_path / "engine.tsv",
        ("1:100", "T", "ENST1", "PLANT_score", "0.7"),
        ("1:100", "T", "ENST1", "PLANT_label", "a"),
        ("1:100", "T", "ENST1", "PLANT_score", "0.5"),
        ("1:100", "T", "ENST1", "PLANT_label", "a"),
        ("1:100", "T", "ENST1", "PLANT_score", "0.5"),
        ("1:100", "T", "ENST1", "PLANT_label", "a"),
        ("1:100", "T", "ENST2", "PLANT_score", "1"),
        ("1:100", "T", "ENST2", "PLANT_score", "1"),
        ("1:900", "G", "ENST1", "PLANT_score", "1"),
    )
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine, "listed": listed}, residuals=residuals, table=None)
    for ident in ("engine", "listed"):
        item = entry(report, "Planted", engine=ident)
        assert counts(item) == (8, 9, 7, 8, 9, 7)
        assert (item["text_only_differences"], item["residual_rows"], item["shared_keys"]) == (2, 3, 2)
        assert residual_lines(residuals, "Planted", ident) == [
            "engine\t1:100\tT\tENST1\tPLANT_label\ta",
            "engine\t1:100\tT\tENST1\tPLANT_score\t0.5",
            "reference\t1:100\tT\tENST1\tPLANT_score\t0.7",
        ]
        allele = entry(report, "Planted", engine=ident, projection="allele")
        assert counts(allele) == (5, 5, 5, 5, 5, 5)
        assert (allele["text_only_differences"], allele["shared_keys"]) == (1, 1)
        assert residual_lines(residuals, "Planted", ident, "allele") == []


def test_percent_encoded_value_equals_its_decoded_form_in_a_tuple_file(tmp_path):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "PLANT_label=a%3Bb%3Dc"))
    engine = tuples(tmp_path / "engine.tsv", ("1:100", "T", "ENST1", "PLANT_label", "a;b=c"))
    assert counts(entry(score(tmp_path, reference, {"engine": engine}), "Planted")) == (1, 1, 1, 1, 1, 1)


def test_location_coordinates_compare_as_numbers(tmp_path):
    reference = output(
        tmp_path / "reference.txt",
        row("21:100-0200", "T", "ENST1", "PLANT_score=1"),
        row("21:²", "T", "ENST1", "PLANT_score=1"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("21:100-200", "T", "ENST1", "PLANT_score=1"),
        row("21:²", "T", "ENST1", "PLANT_score=1"),
    )
    assert counts(entry(score(tmp_path, reference, {"engine": engine}), "Planted")) == (2, 2, 2, 2, 2, 2)


def test_headerless_default_output_and_format_detection(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    bare = tmp_path / "engine.txt"
    bare.write_text("".join(REFERENCE_ROWS), encoding="utf-8")
    report = score(tmp_path, reference, {"engine": bare})
    assert report["engines"][0]["format"] == "default"
    assert counts(entry(report, "Planted")) == (4, 4, 4, 4, 4, 4)
    headed = tmp_path / "headed.tsv"
    headed.write_text(HEADER + "1:100\tT\tENST1\tPLANT_score\t0.5\n", encoding="utf-8")
    with pytest.raises(spt.ScorerError, match="line 2 has 5 columns"):
        score(tmp_path, reference, {"engine": headed})


def test_malformed_inputs_are_refused(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    short = tmp_path / "short.tsv"
    short.write_text("1:100\tT\tENST1\tPLANT_score\t0.5\n1:100\tT\tENST1\tPLANT_score\n", encoding="utf-8")
    with pytest.raises(spt.ScorerError, match="line 2 has 4 columns"):
        score(tmp_path, reference, {"engine": short})
    with pytest.raises(spt.ScorerError, match="line 1 has 5 columns; a default-format row has 14$"):
        score(tmp_path, tuples(tmp_path / "t.tsv", ("1:100", "T", "ENST1", "PLANT_score", "1")), {"engine": reference})
    with pytest.raises(spt.ScorerError, match="no data row"):
        score(tmp_path, output(tmp_path / "empty.txt"), {"engine": reference})
    with pytest.raises(spt.ScorerError, match="no engine"):
        score(tmp_path, reference, {})
    odd = tmp_path / "odd.txt"
    odd.write_text(HEADER + REFERENCE_ROWS[0] + "#Uploaded_variation\tLocation\n", encoding="utf-8")
    with pytest.raises(spt.ScorerError, match="line 3, the column header, names 2 columns, not 14"):
        score(tmp_path, reference, {"engine": odd})
    cut = tmp_path / "cut.txt.gz"
    cut.write_bytes(gzip.compress(reference.read_bytes())[:-12])
    with pytest.raises(spt.ScorerError, match="cut.txt.gz: cannot read the gzip stream"):
        score(tmp_path, cut, {"engine": reference})
    assert not any((tmp_path / "work").iterdir())


@pytest.mark.parametrize("rows", [1, 40_000])
def test_failed_sort_is_reported_with_its_own_message(tmp_path, rows):
    reference = output(
        tmp_path / "reference.txt", *(row(f"1:{n}", "T", "ENST1", "PLANT_score=1") for n in range(rows))
    )
    with pytest.raises(spt.ScorerError, match="sort exited [0-9]+ writing input0"):
        spt.score(reference, [("engine", reference)], ["Planted"], {}, tmp_path / "work", sort_memory="no-such-size")
    assert not any((tmp_path / "work").iterdir())


def test_merge_orders_bytes_whatever_the_callers_locale(tmp_path, monkeypatch):
    monkeypatch.setenv("LC_ALL", "en_US.UTF-8")
    rows = [
        row(location, allele, feature, "PLANT_score=1;PLANT_label=x")
        for location in ("1:117329", "10:110755", "1:1173290", "1:117329-117330", "X:5", "1_KI270706v1_random:9")
        for allele in ("-", "A", "a", "AT")
        for feature in ("-", "ENST1", "ENST10", "ENST1.2")
    ]
    reference = output(tmp_path / "reference.txt", *rows)
    engine = output(tmp_path / "engine.txt", *reversed(rows))
    item = entry(score(tmp_path, reference, {"engine": engine}), "Planted")
    assert counts(item) == (192, 192, 192, 192, 192, 192)


def test_sorted_file_that_is_not_ascending_ends_the_run_in_the_error_line(tmp_path, capsys, monkeypatch):
    """The merge and the cursors read a sorted file forward alone, so past a line that sorts
    below the one before it they would miss the lines it belongs with and leave their values
    without a pair or a class: such a line is a refusal that names the file and the line.
    ``sort`` orders a line that opens the next one first whatever byte follows it there, and
    a line may repeat."""
    head = b"1:100\tT\tPLANT_label\t"
    sorted_lines = [head + b"a\n", head + b"a\x01b\n", head + b"a\x01b\n", head + b"b\n"]
    in_order = tmp_path / "in-order"
    in_order.write_bytes(b"".join(sorted_lines))
    assert list(spt.ascending(in_order)) == sorted_lines
    for name, lines, number in (
        ("descending", [b"1:200\tT\tPLANT_label\ta\n", head + b"a\n"], 2),
        ("opened-line-last", [head + b"a\n", head + b"a\x01b\n", head + b"a\n"], 3),
    ):
        path = tmp_path / name
        path.write_bytes(b"".join(lines))
        message = re.escape(f"{name}: line {number} sorts below line {number - 1}: the file is not in ascending order")
        cursor = spt.Cursor(path)
        with pytest.raises(spt.ScorerError, match=message):
            cursor.starting(b"9:1\t")
        cursor.close()
        for reference, engine in ((path, in_order), (in_order, path)):
            with pytest.raises(spt.ScorerError, match=message):
                spt.merge(reference, engine, defaultdict(int), lambda *values: None, lambda head: None)
    # A whole run whose ``sort`` writes its files in descending order.
    backwards = tmp_path / "sort"
    backwards.write_text('#!/bin/sh\nexec sort -r "$@"\n', encoding="utf-8")
    backwards.chmod(0o755)
    monkeypatch.setattr(spt.shutil, "which", lambda name: str(backwards))
    reference = output(tmp_path / "reference.txt", *(row(f"1:{n}00", "T", "ENST1", f"PLANT_score={n}") for n in (1, 2)))
    engine = output(tmp_path / "engine.txt", *(row(f"1:{n}00", "T", "ENST1", f"PLANT_score={n + 2}") for n in (1, 2)))
    tool = ["--pairing-tool", str(built_pairing_tool())] if THROUGH_PAIRING_TOOL else []
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={engine}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
        "--residuals", str(tmp_path / "residuals"), *tool,
    )  # fmt: skip
    assert (code, lines) == (2, [])
    walked = r"input[01]\.(feature|verdicts)"
    assert re.fullmatch(
        rf"ERROR: \[score_plugin_tuples\] {walked}: line 2 sorts below line 1: the file is not in ascending order\n", error
    )
    assert not (tmp_path / "report.json").exists() and not any((tmp_path / "work").iterdir())
    assert not any((tmp_path / "residuals").glob("*.tsv"))


@pytest.mark.parametrize("side", ["reference", "engine"])
def test_value_with_no_row_under_its_key_among_the_sorted_rows_ends_the_run_in_the_error_line(tmp_path, monkeypatch, side):
    """A value of a default-format input comes from a row under its key, so a class that
    reads the rows finds one there. Where the rows of the key are not read from the sorted
    rows, as when they lie out of a forward-only reader's way, the run ends in a named error
    in place of leaving the value without a class."""
    rows_of = {"reference": ("input0.rows", "input1.rows.reference"), "engine": ("input1.rows", "input1.rows.engine")}

    class RowsOutOfTheWay(spt.Cursor):
        def __init__(self, path):
            super().__init__(path)
            self.hides = path.name in rows_of[side]

        def starting(self, prefix):
            found = super().starting(prefix)
            return [] if self.hides and prefix.startswith(b"1:200\t") else found

    reference = output(
        tmp_path / "reference.txt", row("1:100", "T", "ENST1", "GWAS_pmid=17"), row("1:200", "T", "ENST1", "GWAS_pmid=17")
    )
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "GWAS_pmid=17"), row("1:200", "T", "ENST1", "GWAS_pmid=18"))
    residuals = tmp_path / "residuals"
    item = entry(score(tmp_path, reference, {"engine": engine}, residuals=residuals), "GWAS")
    assert (counts(item), item["masked"]) == ((2, 2, 1, 1, 1, 1), {"planted_reference_value": 1})
    for path in residuals.iterdir():
        path.unlink()
    monkeypatch.setattr(spt, "Cursor", RowsOutOfTheWay)
    written = {"reference": reference, "engine": engine}[side]
    message = (
        f"{written}: no row under Location 1:200, Allele T, Feature ENST1 was read from the sorted rows,"
        " where the input writes a value of GWAS_pmid"
    )
    with pytest.raises(spt.ScorerError, match=re.escape(message)):
        score(tmp_path, reference, {"engine": engine}, residuals=residuals)
    assert not any((tmp_path / "work").iterdir()) and not any(residuals.glob("*.tsv"))


def expected_counts(reference: list[tuple], engine: list[tuple], drop_feature: bool) -> tuple[int, int, int, int, list[str]]:
    """Multiset arithmetic over (Location, Allele, Feature, field, value) tuples whose only
    numerically equal texts are ``0.1`` and ``0.10``, each tuple once where the Feature is
    dropped: the two sides' counts, the pairs, the keys under which a side holds a field more
    than once, and the unpaired tuples as residual lines in the tuples' order."""
    sides = []
    for side in (reference, engine):
        grouped: dict[tuple, Counter] = defaultdict(Counter)
        for location, allele, feature, name, value in side:
            held = grouped[(location, allele, name) if drop_feature else (location, allele, feature, name)]
            held[value] = 1 if drop_feature else held[value] + 1
        sides.append(grouped)
    paired = 0
    left = []
    for key in set(sides[0]) | set(sides[1]):
        ref_values, eng_values = sides[0].get(key, Counter()), sides[1].get(key, Counter())
        paired += (ref_values & eng_values).total()
        ref_left, eng_left = ref_values - eng_values, eng_values - ref_values
        for one, other in (("0.1", "0.10"), ("0.10", "0.1")):
            as_numbers = min(ref_left[one], eng_left[other])
            paired += as_numbers
            ref_left[one] -= as_numbers
            eng_left[other] -= as_numbers
        left += [((*key, value), "reference") for value in ref_left.elements()]
        left += [((*key, value), "engine") for value in eng_left.elements()]
    lines = ["\t".join((side, *columns)) for columns, side in sorted(left, key=lambda item: "\t".join(item[0]).encode())]
    shared = {key[:-1] for side in sides for key, values in side.items() if values.total() > 1}
    totals = [sum(values.total() for values in side.values()) for side in sides]
    return totals[0], totals[1], paired, len(shared), lines


def test_random_outputs_agree_with_multiset_arithmetic(tmp_path):
    rng = random.Random(116)
    locations = ["1:5", "1:50", "1:5-6", "10:5", "2:5", "MT:5"]
    alleles = ["A", "-", "AT", "T"]
    features = ["ENST1", "ENST10", "ENST2", "-"]
    values = ["a", "b", "0.1", "0.10", "7"]

    def side() -> tuple[list[str], list[tuple]]:
        rows, written = [], []
        for _ in range(400):
            key = (rng.choice(locations), rng.choice(alleles), rng.choice(features))
            pairs = [(name, rng.choice(values)) for name in ("PLANT_score", "PLANT_label") if rng.random() < 0.7]
            rows.append(row(*key, ";".join(f"{name}={value}" for name, value in pairs) or "-"))
            written += [(*key, name, value) for name, value in pairs]
        return rows, written

    reference_rows, reference_tuples = side()
    engine_rows, engine_tuples = side()
    reference = output(tmp_path / "reference.txt", *reference_rows)
    as_output = output(tmp_path / "engine.txt", *engine_rows)
    as_tuples = tuples(tmp_path / "engine.tsv", *engine_tuples)
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"output": as_output, "tuples": as_tuples}, residuals=residuals)
    for projection, drop_feature in (("feature", False), ("allele", True)):
        perl, rust, paired, shared, lines = expected_counts(reference_tuples, engine_tuples, drop_feature)
        assert paired < min(perl, rust) and shared > 0
        assert len(lines) == perl + rust - 2 * paired
        for engine in ("output", "tuples"):
            item = entry(report, "Planted", engine=engine, projection=projection)
            assert counts(item) == (perl, rust, paired, perl, rust, paired), (engine, projection)
            assert (item["residual_rows"], item["shared_keys"]) == (len(lines), shared)
            assert residual_lines(residuals, "Planted", engine, projection) == lines


def run_main(capsys, *arguments: str) -> tuple[int, list[str], str]:
    code = spt.main(list(arguments))
    captured = capsys.readouterr()
    return code, captured.out.splitlines(), captured.err


def test_summary_prints_each_class_count_and_marks_none_and_above_population(tmp_path, capsys):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    changed = REFERENCE_ROWS[1].replace("MaxEntScan_diff=1.5", "MaxEntScan_diff=-2.25")
    filled = REFERENCE_ROWS[2].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;LOEUF=0")
    engine = output(tmp_path / "engine.txt", REFERENCE_ROWS[0], changed, filled, *REFERENCE_ROWS[3:])
    table = class_table(tmp_path, CLASS_TABLE.replace("GWAS_*\tmet\t-\n", "GWAS_*\n"))
    report_path = tmp_path / "report.json"
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={engine}", "--report", str(report_path),
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--class-table", str(table),
        "--classes", "planted_two_numbers, planted_reference_value,order_only",
        "--class-population", "planted_two_numbers=0", "--class-population", "planted_reference_value=5",
        "--plugins", ", ".join(PLANTED), "--engine-binary", f"engine={'ab' * 32}", "--engine-commit", "engine=0123abc",
    )  # fmt: skip
    assert (code, error) == (0, "")
    assert "keys of the reference that no source claims: Custom (1), Custom_SIG (1)" in lines
    assert "keys of the engine engine that no source claims: Custom (1), Custom_SIG (1)" in lines
    exceeded = "class planted_two_numbers set aside 1 (above its reachable population of 0) [met, -]"
    assert f"MaxEntScan\tengine\tfeature\t{exceeded}" in lines
    assert "GWAS\tengine\tfeature\tclass planted_reference_value set aside 0 (none)" in lines
    assert "GWAS\tengine\tallele\tclass planted_reference_value set aside 0 (none)" in lines
    held = "class planted_engine_zero is not applied and would set aside 1: raw counts only [met, -]"
    assert f"LOEUF\tengine\tfeature\t{held}" in lines
    assert "LOEUF\tengine\tfeature\tperl=1 rust=2 intersection=1 text_only_differences=0 residual_rows=1 shared_keys=0" in lines
    assert not any("Planted\t" in line and "class" in line for line in lines)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert [item for item in report["classes"] if "population" in item] == [
        {
            "class": "planted_two_numbers",
            "plugins": {"MaxEntScan": ["MaxEntScan_*"]},
            "applied": True,
            "further_columns": ["met", "-"],
            "population": 0,
        },
        {
            "class": "planted_reference_value",
            "plugins": {"GWAS": ["GWAS_*"]},
            "applied": True,
            "further_columns": [],
            "population": 5,
        },
    ]
    assert entry(report, "MaxEntScan")["masked"] == {"planted_two_numbers": 1}
    assert report["plugins"] == list(PLANTED)
    described = report["engines"][0]
    assert (described["id"], described["binary_sha256"], described["commit"]) == ("engine", "ab" * 32, "0123abc")
    assert "binary_sha256" not in report["reference"]


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["--classes", "planted_engine_zero"], "--classes names planted_engine_zero, a class the class table does not"),
        (["--class-table", "{missing}"], "not a file: {missing}"),
        (["--class-table", "{table}", "--class-population", "planted_engine_zero=many"], "takes CLASS=N"),
        (["--class-table", "{table}", "--class-population", "no_such_class=1"], "--class-population names no_such_class"),
        (["--plugins", "Planted,Nope"], "--plugins names Nope, no plugin of the registry"),
        (["--plugins", ""], "no source to score"),
        (["--engine-binary", "engine=abc"], "--engine-binary takes ID=SHA256, not engine=abc"),
        (["--engine-commit", "engine=not-hex"], "--engine-commit takes ID=COMMIT, not engine=not-hex"),
        (["--engine-commit", "other=0123abc"], "other is no engine of this run"),
    ],
)
def test_class_options_are_refused_with_one_error_line(tmp_path, capsys, arguments, message):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    names = {"missing": str(tmp_path / "missing.tsv"), "table": str(class_table(tmp_path, CLASS_TABLE))}
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={reference}", "--report", str(tmp_path / "report.json"),
        "--work-dir", str(tmp_path / "work"), "--plugins", "Planted",
        *(argument.format(**names) for argument in arguments),
    )  # fmt: skip
    assert (code, lines) == (2, [])
    assert error.startswith("ERROR: [score_plugin_tuples] ") and error.count("\n") == 1
    assert message.format(**names) in error
    assert not (tmp_path / "report.json").exists() and not (tmp_path / "work").exists()


def test_plugin_is_named_by_its_registry_name_or_an_alias_the_registry_lists(tmp_path, monkeypatch):
    monkeypatch.setattr(compare_plugin_outputs, "PLUGIN_ALIASES", {"Loaded": "Planted"})
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "PLANT_score=1"))
    report = score(tmp_path, reference, {"engine": reference}, plugins=["Loaded", "Planted"], table=None)
    assert report["plugins"] == ["Planted"]
    assert counts(entry(report, "Planted")) == (1, 1, 1, 1, 1, 1)


def test_field_the_registry_gives_a_plugin_is_scored_though_named_like_a_core_key(tmp_path, monkeypatch):
    plant_registry(monkeypatch, {**PLANTED, "Planted": [*PLANTED["Planted"], "IMPACT"]})
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "IMPACT=LOW;PLANT_score=1;LOEUF=0.1"))
    with_the_plugin = score(tmp_path, reference, {"engine": reference}, plugins=["Planted", "LOEUF"])
    assert with_the_plugin["sources"]["Planted"] == ["IMPACT", "PLANT_score"]
    assert counts(entry(with_the_plugin, "Planted")) == (2, 2, 2, 2, 2, 2)
    assert with_the_plugin["reference"]["unclaimed_keys"] == {}
    without_it = score(tmp_path, reference, {"engine": reference}, plugins=["LOEUF"])
    assert without_it["sources"] == {"LOEUF": ["LOEUF"]}
    assert without_it["reference"]["unclaimed_keys"] == {"PLANT_score": 1}


def test_run_command_is_handed_to_the_registry_which_can_refuse_the_run(tmp_path, capsys, monkeypatch):
    asked = []
    refusal = "the run loads Planted and sets --planted: PLANT_label is both the plugin's column and Ensembl VEP's own key"

    def plugin_of_field(field, plugins, flags=None):
        asked.append(flags)
        if flags is not None and field == "PLANT_label" and "--planted" in flags:
            raise ValueError(refusal)
        return "Planted" if field.startswith("PLANT_") else None

    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "PLANT_score=1;PLANT_label=x"))
    report_path = tmp_path / "report.json"

    def run(name, *more):
        asked.clear()
        return run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
            "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(report_path),
            "--residuals", str(tmp_path / name), *more,
        )  # fmt: skip

    code, _, error = run("unasked")
    assert (code, error) == (0, "")
    assert json.loads(report_path.read_text(encoding="utf-8"))["run_command"] is None
    monkeypatch.setattr(compare_plugin_outputs, "plugin_of_field", plugin_of_field)
    code, _, error = run("handed", "--run-command", "vep --plugin 'Planted,a file' --offline")
    assert (code, error) == (0, "")
    words = ["vep", "--plugin", "Planted,a file", "--offline"]
    assert set(asked) == {tuple(words)}
    assert json.loads(report_path.read_text(encoding="utf-8"))["run_command"] == words
    code, lines, error = run("refused", "--run-command", "vep --planted")
    assert (code, lines) == (2, [])
    assert error == f"ERROR: [score_plugin_tuples] {refusal}\n"
    assert not any((tmp_path / "work").iterdir()) and not any((tmp_path / "refused").iterdir())
    code, lines, error = run("unbalanced", "--run-command", "vep 'a file")
    assert (code, lines) == (2, [])
    assert error == """ERROR: [score_plugin_tuples] --run-command takes the run's command line, not "vep 'a file"\n"""


def test_allele_tuple_is_credited_to_the_first_class_in_table_order(tmp_path, monkeypatch):
    plant_registry(monkeypatch, {"GWAS": ["GWAS_pmid"]})
    table = "planted_row_cdna_one\tGWAS\tGWAS_pmid\nplanted_reference_value\tGWAS\tGWAS_pmid\n"
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "GWAS_pmid=17"),
        row("1:100", "T", "ENST2", "GWAS_pmid=17", cdna="1"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "IMPACT=LOW"),
        row("1:100", "T", "ENST2", "IMPACT=LOW", cdna="1"),
    )
    report = score(tmp_path, reference, {"engine": engine}, table=table)
    assert entry(report, "GWAS")["masked"] == {"planted_row_cdna_one": 1, "planted_reference_value": 1}
    allele = entry(report, "GWAS", projection="allele")
    assert counts(allele) == (1, 0, 0, 0, 0, 0)
    assert allele["masked"] == {"planted_row_cdna_one": 1, "planted_reference_value": 0}
    held = entry(score(tmp_path, reference, {"engine": engine}, table=table, applied=[]), "GWAS", projection="allele")
    assert held["not_applied"] == {"planted_row_cdna_one": 1}


@own_path_alone
def test_no_rows_are_kept_when_no_class_is_declared(tmp_path, monkeypatch):
    sorted_files = []
    sinks = spt.Sorter.sinks

    def recording(self, *paths, **options):
        sorted_files.extend(path.name for path in paths)
        return sinks(self, *paths, **options)

    monkeypatch.setattr(spt.Sorter, "sinks", recording)
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    score(tmp_path, reference, {"engine": reference}, table=None)
    assert sorted_files == ["input0.feature", "input0.allele", "input1.feature", "input1.allele", "input1.verdicts"]
    sorted_files.clear()
    score(tmp_path, reference, {"engine": reference})
    assert "input0.rows" in sorted_files and "input1.rows" in sorted_files


def test_every_input_is_checked_before_any_is_scored(tmp_path, monkeypatch):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    narrow = tmp_path / "narrow.tsv"
    narrow.write_text("1:100\tT\tENST1\n", encoding="utf-8")
    cut = tmp_path / "cut.txt.gz"
    cut.write_bytes(gzip.compress(reference.read_bytes())[:20])
    headed = tmp_path / "headed.txt"
    headed.write_text("##a comment\n#Uploaded_variation\tLocation\n" + REFERENCE_ROWS[0], encoding="utf-8")
    monkeypatch.setattr(
        spt, "extract", lambda *arguments, **options: pytest.fail("an input was read before all were checked")
    )
    with pytest.raises(spt.ScorerError, match="headed.txt: line 2, the column header, names 2 columns, not 14"):
        score(tmp_path, reference, {"engine": reference, "later": headed}, residuals=tmp_path / "residuals")
    with pytest.raises(spt.ScorerError, match="narrow.tsv: line 1 has 3 columns; a default-format row has 14, a tuple"):
        score(tmp_path, reference, {"engine": reference, "later": narrow}, residuals=tmp_path / "residuals")
    with pytest.raises(spt.ScorerError, match="cut.txt.gz: cannot read the gzip stream"):
        score(tmp_path, reference, {"engine": reference, "later": cut}, residuals=tmp_path / "residuals")
    with pytest.raises(spt.ScorerError, match="narrow.tsv: line 1 has 3 columns; a default-format row has 14$"):
        score(tmp_path, narrow, {"engine": reference})
    assert not (tmp_path / "work").exists() and not (tmp_path / "residuals").exists()


def test_distinct_unclaimed_keys_are_bounded(tmp_path, monkeypatch):
    monkeypatch.setattr(spt, "RESOLVED_LIMIT", 2)
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "One=1;Two=2;PLANT_score=1"))
    assert score(tmp_path, reference, {"engine": reference})["reference"]["unclaimed_keys"] == {"One": 1, "Two": 1}
    wide = output(tmp_path / "wide.txt", row("1:100", "T", "ENST1", "One=1;Two=2;Three=3;PLANT_score=1"))
    listed = tuples(tmp_path / "wide.tsv", *(("1:100", "T", "ENST1", name, "1") for name in ("One", "Two", "Three")))
    for path in (wide, listed):
        with pytest.raises(spt.ScorerError, match=f"{path.name}: more than 2 distinct keys that no source claims"):
            score(tmp_path, reference, {"engine": path})
        assert not any((tmp_path / "work").iterdir())


def test_interrupt_ends_in_one_error_line_and_leaves_nothing_behind(tmp_path, capsys, monkeypatch):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    changed = REFERENCE_ROWS[0].replace("PLANT_label=benign", "PLANT_label=pathogenic")
    engine = output(tmp_path / "engine.txt", changed, *REFERENCE_ROWS[1:])
    merge = spt.merge
    walks = []

    def interrupted(reference_path, *others):
        walks.append(reference_path.suffix)
        if len(walks) == 4:
            os.kill(os.getpid(), signal.SIGTERM)
        return merge(reference_path, *others)

    monkeypatch.setattr(spt, "merge", interrupted)
    report_path = tmp_path / "report.json"

    def callers(number, frame):
        pytest.fail("the caller's handler ran during the run")

    before = signal.signal(signal.SIGTERM, callers)
    try:
        code, lines, error = run_main(
            capsys,
            "--reference", str(reference), "--engine", f"first={engine}", "--engine", f"second={engine}",
            "--plugins", "Planted", "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M",
            "--report", str(report_path), "--residuals", str(tmp_path / "residuals"),
        )  # fmt: skip
        restored = signal.getsignal(signal.SIGTERM)
    finally:
        signal.signal(signal.SIGTERM, before)
    assert (code, lines) == (130, [])
    assert error == "ERROR: [score_plugin_tuples] interrupted before the report was written\n"
    assert walks == [".feature", ".allele", ".feature", ".allele"]
    assert not report_path.exists()
    assert not any((tmp_path / "work").iterdir())
    assert not any((tmp_path / "residuals").iterdir())
    assert restored is callers


INTERRUPT, TERMINATION = signal.SIGINT, signal.SIGTERM


@pytest.mark.parametrize(
    "points",
    [
        pytest.param([(spt, "canonical", TERMINATION)], id="a termination while an input is read"),
        pytest.param([(spt.Sink, "close", INTERRUPT)], id="an interrupt while a sort is waited for"),
        pytest.param(
            [(spt, "canonical", TERMINATION), (spt.Sink, "reap", TERMINATION)],
            id="a second termination while the sorts are stopped",
        ),
        pytest.param(
            [(spt, "canonical", INTERRUPT), (spt.Sink, "reap", TERMINATION)],
            id="a termination after an interrupt while the sorts are stopped",
        ),
        pytest.param(
            [(spt, "canonical", TERMINATION), (shutil, "rmtree", INTERRUPT)],
            id="an interrupt after a termination while the scratch is removed",
        ),
        pytest.param(
            [(spt, "canonical", INTERRUPT), (shutil, "rmtree", INTERRUPT)],
            id="a second interrupt while the scratch is removed",
        ),
    ],
)
def test_stop_signal_stops_every_sort_and_a_repeat_changes_nothing(tmp_path, capsys, monkeypatch, points):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    started = []
    popen = subprocess.Popen

    def recording(*arguments, **options):
        started.append(popen(*arguments, **options))
        return started[-1]

    def signalling(real, number):
        calls = []

        def signalled(*arguments, **options):
            calls.append(arguments)
            if len(calls) == 1:
                os.kill(os.getpid(), number)
            return real(*arguments, **options)

        return signalled

    monkeypatch.setattr(subprocess, "Popen", recording)
    for owner, name, number in points:
        monkeypatch.setattr(owner, name, signalling(getattr(owner, name), number))
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
    )  # fmt: skip
    assert (code, lines) == (130, [])
    assert error == "ERROR: [score_plugin_tuples] interrupted before the report was written\n"
    assert [process.returncode for process in started] == [-signal.SIGKILL, -signal.SIGKILL]
    assert not any((tmp_path / "work").iterdir())


def test_signal_the_process_ignores_stays_ignored(tmp_path, capsys, monkeypatch):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    canonical = spt.canonical

    def signalled(location):
        os.kill(os.getpid(), signal.SIGINT)
        return canonical(location)

    monkeypatch.setattr(spt, "canonical", signalled)
    before = signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        code, _, error = run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
            "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
        )  # fmt: skip
        left = signal.getsignal(signal.SIGINT)
    finally:
        signal.signal(signal.SIGINT, before)
    assert (code, error) == (0, "")
    assert left is signal.SIG_IGN


def test_termination_while_the_report_is_put_together_leaves_no_residual_file(tmp_path, capsys, monkeypatch):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    input_entry = spt.input_entry

    def signalled(extracted):
        os.kill(os.getpid(), signal.SIGTERM)
        return input_entry(extracted)

    monkeypatch.setattr(spt, "input_entry", signalled)
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
        "--residuals", str(tmp_path / "residuals"),
    )  # fmt: skip
    assert (code, lines) == (130, [])
    assert error == "ERROR: [score_plugin_tuples] interrupted before the report was written\n"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["reference.txt", "residuals", "work"]
    assert not any((tmp_path / "residuals").iterdir())


def test_stop_signal_while_the_scratch_is_removed_after_scoring_leaves_nothing_behind(tmp_path, capsys, monkeypatch):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    rmtree = shutil.rmtree
    removals = []

    def signalled(*arguments, **options):
        removals.append(arguments)
        if len(removals) == 1:
            os.kill(os.getpid(), signal.SIGTERM)
        return rmtree(*arguments, **options)

    monkeypatch.setattr(shutil, "rmtree", signalled)
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
        "--residuals", str(tmp_path / "residuals"),
    )  # fmt: skip
    assert (code, lines) == (130, [])
    assert error == "ERROR: [score_plugin_tuples] interrupted before the report was written\n"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["reference.txt", "residuals", "work"]
    assert not any((tmp_path / "work").iterdir())
    assert not any((tmp_path / "residuals").iterdir())


def test_scratch_that_cannot_be_removed_is_named_and_leaves_no_residual_file(tmp_path, capsys, monkeypatch):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)

    def busy(*arguments, **options):
        raise OSError(16, "Device or resource busy")

    monkeypatch.setattr(shutil, "rmtree", busy)
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
        "--residuals", str(tmp_path / "residuals"),
    )  # fmt: skip
    assert (code, lines) == (2, [])
    assert error == "ERROR: [score_plugin_tuples] [Errno 16] Device or resource busy\n"
    assert not (tmp_path / "report.json").exists()
    assert not any((tmp_path / "residuals").iterdir())


@pytest.mark.parametrize("earlier", [None, b"an earlier report\n"], ids=["at a free path", "over an earlier report"])
def test_report_that_cannot_be_written_is_named_and_leaves_no_residual_file(tmp_path, capsys, monkeypatch, earlier):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    report_path = tmp_path / "report.json"
    if earlier is not None:
        report_path.write_bytes(earlier)

    def full(self, data, *arguments, **options):
        assert sorted(path.name for path in (tmp_path / "residuals").iterdir()) == [
            "Planted.engine.allele.tsv",
            "Planted.engine.feature.tsv",
        ]
        self.write_bytes(data.encode()[:20])
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(Path, "write_text", full)
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(report_path),
        "--residuals", str(tmp_path / "residuals"),
    )  # fmt: skip
    assert (code, lines) == (2, [])
    assert error == f"ERROR: [score_plugin_tuples] cannot write the report {report_path}: No space left on device\n"
    assert not any((tmp_path / "residuals").iterdir())
    kept = [] if earlier is None else ["report.json"]
    assert sorted(path.name for path in tmp_path.iterdir()) == ["reference.txt", *kept, "residuals", "work"]
    assert earlier is None or report_path.read_bytes() == earlier


def test_report_replaces_an_earlier_one_and_takes_its_mode(tmp_path, capsys):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    report_path = tmp_path / "report.json"
    report_path.write_bytes(b"an earlier report\n")
    report_path.chmod(0o600)
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(report_path),
    )  # fmt: skip
    assert (code, error) == (0, "")
    assert lines[0].startswith("Planted\tengine\tfeature\tperl=4 rust=4 intersection=4 ")
    assert json.loads(report_path.read_bytes())["plugins"] == ["Planted"]
    assert report_path.stat().st_mode & 0o777 == 0o600
    assert sorted(path.name for path in tmp_path.iterdir()) == ["reference.txt", "report.json", "work"]


@pytest.mark.skipif(os.geteuid() == 0, reason="the superuser writes in any directory")
@pytest.mark.parametrize("protected", ["the report", "its directory"])
def test_earlier_report_that_cannot_be_replaced_is_refused_before_scoring(tmp_path, capsys, protected):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    report_path = tmp_path / "reports" / "report.json"
    report_path.parent.mkdir()
    report_path.write_bytes(b"an earlier report\n")
    locked = report_path if protected == "the report" else report_path.parent
    mode = locked.stat().st_mode & 0o777
    locked.chmod(mode & 0o555)
    try:
        code, lines, error = run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
            "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(report_path),
        )  # fmt: skip
    finally:
        locked.chmod(mode)
    assert (code, lines) == (2, [])
    assert error == f"ERROR: [score_plugin_tuples] cannot write the report {report_path}\n"
    assert report_path.read_bytes() == b"an earlier report\n"
    assert not (tmp_path / "work").exists()


@pytest.mark.parametrize("shape", ["a link", "a pipe"])
def test_report_path_that_is_no_regular_file_is_written_as_it_is(tmp_path, capsys, shape):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    report_path = tmp_path / "report.json"
    if shape == "a link":
        (tmp_path / "linked.json").write_bytes(b"an earlier report\n")
        report_path.symlink_to("linked.json")
    else:
        os.mkfifo(report_path)
    # A pipe takes a writer only while a reader holds it open, and this reader never blocks.
    reader = os.open(report_path, os.O_RDONLY | os.O_NONBLOCK)
    try:
        code, lines, error = run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
            "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(report_path),
        )  # fmt: skip
        written = os.read(reader, 1 << 20)
    finally:
        os.close(reader)
    assert (code, error) == (0, "")
    assert lines[0].startswith("Planted\tengine\tfeature\tperl=4 rust=4 intersection=4 ")
    assert json.loads(written)["plugins"] == ["Planted"]
    assert report_path.is_symlink() if shape == "a link" else report_path.is_fifo()


@pytest.mark.parametrize("stream", [None, io.StringIO])
def test_summary_needs_no_standard_output_it_can_reconfigure(tmp_path, monkeypatch, stream):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    monkeypatch.setattr(sys, "stdout", stream and stream())
    arguments = [
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
    ]  # fmt: skip
    assert spt.main(arguments) == 0
    assert (tmp_path / "report.json").exists()
    if stream:
        assert sys.stdout.getvalue().startswith("Planted\tengine\tfeature\tperl=4 rust=4 intersection=4 ")


def test_main_runs_off_the_main_thread(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    arguments = [
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
    ]  # fmt: skip
    with ThreadPoolExecutor(max_workers=1) as pool:
        assert pool.submit(spt.main, arguments).result() == 0
    assert (tmp_path / "report.json").exists()


def test_main_runs_under_a_termination_handler_that_was_not_set_from_python(tmp_path, monkeypatch):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    handlers = []

    # Such a handler reads as None, which ``signal.signal`` hands back and refuses to set.
    def refusing_none(number, handler):
        handlers.append(handler)
        if handler is None:
            raise TypeError("signal handler must be signal.SIG_IGN, signal.SIG_DFL, or a callable object")

    monkeypatch.setattr(signal, "getsignal", lambda number: None)
    monkeypatch.setattr(signal, "signal", refusing_none)
    arguments = [
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
    ]  # fmt: skip
    assert spt.main(arguments) == 0
    assert handlers == []
    assert (tmp_path / "report.json").exists()


def test_plugins_of_a_run_are_resolved_only_by_a_registry_with_its_alias_table(tmp_path, monkeypatch):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "PLANT_score=1;Custom=id1"))
    monkeypatch.delattr(compare_plugin_outputs, "PLUGIN_ALIASES")
    with pytest.raises(AttributeError, match="PLUGIN_ALIASES"):
        score(tmp_path, reference, {"engine": reference}, plugins=["Planted"], table=None)
    custom = score(tmp_path, reference, {"engine": reference}, plugins=[], custom={"Custom": "Custom"}, table=None)
    assert counts(entry(custom, "Custom")) == (1, 1, 1, 1, 1, 1)


CATALOG_TABLE = "planted_catalog_place\tGWAS\tGWAS_*\tmet\t-\nplanted_engine_zero\tLOEUF\tLOEUF\tmet\t-\n"


def plant_catalog(monkeypatch) -> tuple[list, list]:
    """Plants a reader of a catalog of ``location<TAB>allele`` lines and a class that accepts
    a row at one of its places, and returns the paths read and, per predicate call, the
    class, the row's Location and the keywords handed beyond the record's rows."""
    reads, handed = [], []

    class Catalog:
        def __init__(self, places):
            self.places = places

        @classmethod
        def from_catalog(cls, path):
            reads.append(path)
            lines = path.read_text(encoding="utf-8").splitlines()
            if lines[:1] != ["location\tallele"]:
                raise formats.ComparatorError(f"{path}: not a GWAS catalog associations file")
            return cls({tuple(line.split("\t")) for line in lines[1:]})

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None, **more):
        handed.append((klass, e.location, more))
        if klass == "planted_catalog_place":
            return "gwas_pairs" in more and (e.location, e.allele) in more["gwas_pairs"].places
        return klass == "planted_engine_zero" and formats.blank(expected_value) and actual_value == "0"

    monkeypatch.setattr(formats, "GwasCatalogPairs", Catalog, raising=False)
    monkeypatch.setattr(formats, "FIELD_CLASSES", dict.fromkeys((*PLANTED_CLASSES, "planted_catalog_place"), ()))
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    return reads, handed


def catalog_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    catalog = tmp_path / "associations.tsv"
    catalog.write_text("location\tallele\n1:100\tT\n", encoding="utf-8")
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "GWAS_pmid=17"),
        row("1:200", "T", "ENST1", "GWAS_pmid=18"),
        row("1:300", "T", "ENST1", "IMPACT=LOW"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "IMPACT=LOW"),
        row("1:200", "T", "ENST1", "IMPACT=LOW"),
        row("1:300", "T", "ENST1", "LOEUF=0"),
    )
    return catalog, reference, engine


def test_gwas_catalog_is_read_once_and_handed_to_the_classes_of_gwas_alone(tmp_path, monkeypatch):
    reads, handed = plant_catalog(monkeypatch)
    catalog, reference, engine = catalog_inputs(tmp_path)
    engines = {"engine": engine, "listed": tuples(tmp_path / "engine.tsv", ("1:300", "T", "ENST1", "LOEUF", "0"))}
    plugins = ["GWAS", "LOEUF"]
    without = score(tmp_path, reference, engines, residuals=tmp_path / "without", table=CATALOG_TABLE, plugins=plugins)
    assert without["gwas_catalog"] is None and reads == []
    assert all(more == {} for _, _, more in handed)
    for ident in engines:
        assert counts(entry(without, "GWAS", engine=ident)) == (2, 0, 0, 2, 0, 0)
        assert entry(without, "GWAS", engine=ident)["masked"] == {"planted_catalog_place": 0}
        assert len(residual_lines(tmp_path / "without", "GWAS", ident)) == 2
    handed.clear()
    report = score(
        tmp_path, reference, engines, residuals=tmp_path / "with", table=CATALOG_TABLE, plugins=plugins, catalog=catalog
    )
    assert report["gwas_catalog"] == {"path": str(catalog), "sha256": hashlib.sha256(catalog.read_bytes()).hexdigest()}
    assert reads == [catalog]
    for ident in engines:
        item = entry(report, "GWAS", engine=ident)
        assert counts(item) == (2, 0, 0, 1, 0, 0)
        assert item["masked"] == {"planted_catalog_place": 1}
        assert residual_lines(tmp_path / "with", "GWAS", ident) == ["reference\t1:200\tT\tENST1\tGWAS_pmid\t18"]
        allele = entry(report, "GWAS", engine=ident, projection="allele")
        assert (counts(allele), allele["masked"]) == ((2, 0, 0, 1, 0, 0), {"planted_catalog_place": 1})
        zero = entry(report, "LOEUF", engine=ident)
        assert (counts(zero), zero["masked"]) == ((0, 1, 0, 0, 0, 0), {"planted_engine_zero": 1})
    of_gwas = [more for klass, _, more in handed if klass == "planted_catalog_place"]
    assert len(of_gwas) == 4 and all(list(more) == ["gwas_pairs"] for more in of_gwas)
    assert len({id(more["gwas_pairs"]) for more in of_gwas}) == 1
    assert [more for klass, _, more in handed if klass == "planted_engine_zero"] == [{}, {}]


def test_class_of_gwas_that_is_not_applied_is_handed_the_catalog_and_listed_as_not_applied(tmp_path, monkeypatch):
    _, handed = plant_catalog(monkeypatch)
    catalog, reference, engine = catalog_inputs(tmp_path)
    engines = {"engine": engine, "listed": tuples(tmp_path / "engine.tsv", ("1:300", "T", "ENST1", "LOEUF", "0"))}
    residuals, held = tmp_path / "residuals", ["planted_engine_zero"]
    report = score(
        tmp_path, reference, engines, residuals=residuals, table=CATALOG_TABLE, applied=held, catalog=catalog
    )
    for ident in engines:
        for projection in ("feature", "allele"):
            item = entry(report, "GWAS", engine=ident, projection=projection)
            assert tuple(item[key] for key in RAW_KEYS) == (2, 0, 0)
            assert not any(key.startswith("adj_") for key in item)
            assert item["masked"] == {} and item["not_applied"] == {"planted_catalog_place": 1}
            assert item["residual_rows"] == 2
        assert residual_lines(residuals, "GWAS", ident) == [
            "reference\t1:100\tT\tENST1\tGWAS_pmid\t17",
            "reference\t1:200\tT\tENST1\tGWAS_pmid\t18",
        ]
    of_gwas = [more for klass, _, more in handed if klass == "planted_catalog_place"]
    assert len(of_gwas) == 4 and all(list(more) == ["gwas_pairs"] for more in of_gwas)


def test_summary_says_that_classes_of_gwas_set_nothing_aside_without_the_catalog(tmp_path, capsys, monkeypatch):
    plant_catalog(monkeypatch)
    catalog, reference, engine = catalog_inputs(tmp_path)
    table = class_table(tmp_path, CATALOG_TABLE)
    notice = "no --gwas-catalog: a class of GWAS that reads the catalog sets nothing aside"

    def run(name, *more):
        return run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={engine}", "--plugins", "GWAS,LOEUF",
            "--class-table", str(table), "--work-dir", str(tmp_path / f"work-{name}"), "--sort-memory", "16M",
            "--report", str(tmp_path / "report.json"), "--residuals", str(tmp_path / name), *more,
        )  # fmt: skip

    code, lines, error = run("without")
    assert (code, error) == (0, "")
    assert lines[-1] == notice
    assert "GWAS\tengine\tfeature\tclass planted_catalog_place set aside 0 (none) [met, -]" in lines
    code, lines, error = run("with", "--gwas-catalog", str(catalog))
    assert (code, error) == (0, "")
    assert notice not in lines
    assert "GWAS\tengine\tfeature\tclass planted_catalog_place set aside 1 [met, -]" in lines
    code, lines, error = run("unapplied", "--classes", "planted_engine_zero")
    assert (code, error) == (0, "")
    assert notice not in lines
    code, lines, error = run("unloaded", "--plugins", "LOEUF")
    assert (code, error) == (0, "")
    assert notice not in lines
    code, lines, error = run("missing", "--gwas-catalog", str(tmp_path / "nowhere.tsv"))
    assert (code, lines) == (2, [])
    assert error == f"ERROR: [score_plugin_tuples] not a file: {tmp_path / 'nowhere.tsv'}\n"
    other = tmp_path / "other.tsv"
    other.write_text("a\tb\n", encoding="utf-8")
    code, lines, error = run("refused", "--gwas-catalog", str(other))
    assert (code, lines) == (2, [])
    assert error == f"ERROR: [score_plugin_tuples] {other}: not a GWAS catalog associations file\n"
    assert not any((tmp_path / name).exists() for name in ("missing", "work-missing", "refused", "work-refused"))


SWEEP_TABLE = "planted_swept_terms\tUTRAnnotator\t5UTR_consequence\tmet\t-\nplanted_engine_zero\tLOEUF\tLOEUF\tmet\t-\n"


def plant_sweep(monkeypatch) -> tuple[list, list]:
    """Plants a reader of seeded runs' outputs and a class that accepts a pair of terms those
    runs both write at the row's Location, and returns the lists of paths read and, per
    predicate call, the class, the row's Location and the keywords handed beyond the rows."""
    reads, handed = [], []

    class Sweep:
        def __init__(self, terms):
            self.terms = terms

        @classmethod
        def from_outputs(cls, paths):
            reads.append(list(paths))
            terms = set()
            for path in paths:
                written = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
                if HEADER.encode() not in written:
                    raise formats.ComparatorError(f"{path}: not a default-format output")
                terms |= {(e.location, e.fields.get("5UTR_consequence")) for e in formats.parse_default(written.decode())}
            return cls(terms)

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None, **more):
        handed.append((klass, e.location, more))
        if klass == "planted_swept_terms":
            terms = more["utrannotator_sweep"].terms if "utrannotator_sweep" in more else set()
            return {(e.location, expected_value), (e.location, actual_value)} <= terms
        return klass == "planted_engine_zero" and formats.blank(expected_value) and actual_value == "0"

    monkeypatch.setattr(formats, "UtrannotatorSeedSweep", Sweep, raising=False)
    monkeypatch.setattr(formats, "FIELD_CLASSES", dict.fromkeys((*PLANTED_CLASSES, "planted_swept_terms"), ()))
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    return reads, handed


def sweep_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """A reference and an engine that write different terms at two Locations, and a directory
    of seeded runs that write both terms at the first Location alone, beside files that are
    no output."""
    def rows(*terms: str, extra: str = "") -> list[str]:
        return [row(f"1:{at}", "T", "ENST1", f"5UTR_consequence={term}{extra}") for at, term in zip((100, 200), terms)]

    reference = output(tmp_path / "reference.txt", *rows("lost", "lost"), row("1:300", "T", "ENST1", "IMPACT=LOW"))
    engine = output(tmp_path / "engine.txt", *rows("gained", "gained"), row("1:300", "T", "ENST1", "LOEUF=0"))
    runs = tmp_path / "sweep"
    (runs / "tools").mkdir(parents=True)
    output(runs / "seed-a.txt", *rows("lost", "lost"))
    commented = "## ENSEMBL VARIANT EFFECT PREDICTOR\n## Column descriptions:\n" + HEADER + "".join(rows("gained", "other"))
    (runs / "seed-b.txt.gz").write_bytes(gzip.compress(commented.encode()))
    output(runs / "tools" / "seed-c.txt", *rows("gained", "gained"))
    (runs / "seed-a.command.txt").write_text("vep --plugin UTRAnnotator\n", encoding="utf-8")
    (runs / "input.vcf").write_text(VCF_HEADER, encoding="utf-8")
    (runs / "seed-a.log").write_bytes(b"")
    (runs / "seeds.bin").write_bytes(bytes(range(256)) * 400)
    return reference, engine, runs


@pytest.mark.parametrize("as_directory", [True, False], ids=["a directory", "each output"])
def test_seeded_runs_are_read_once_and_handed_to_the_classes_of_utrannotator_alone(tmp_path, monkeypatch, as_directory):
    reads, handed = plant_sweep(monkeypatch)
    reference, engine, runs = sweep_inputs(tmp_path)
    outputs = [runs / "seed-a.txt", runs / "seed-b.txt.gz"]
    engines = {"engine": engine, "again": engine}
    plugins = ["UTRAnnotator", "LOEUF"]
    without = score(tmp_path, reference, engines, residuals=tmp_path / "without", table=SWEEP_TABLE, plugins=plugins)
    assert without["utrannotator_sweep"] is None and reads == []
    assert all(more == {} for _, _, more in handed)
    assert counts(entry(without, "UTRAnnotator")) == (2, 2, 0, 2, 2, 0)
    assert entry(without, "UTRAnnotator")["masked"] == {"planted_swept_terms": 0}
    handed.clear()
    report = score(
        tmp_path, reference, engines, residuals=tmp_path / "with", table=SWEEP_TABLE, plugins=plugins,
        sweep=[runs] if as_directory else outputs,
    )  # fmt: skip
    assert report["utrannotator_sweep"] == [
        {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in outputs
    ]
    assert reads == [outputs]
    for ident in engines:
        item = entry(report, "UTRAnnotator", engine=ident)
        assert (counts(item), item["masked"]) == ((2, 2, 0, 1, 1, 0), {"planted_swept_terms": 1})
        assert residual_lines(tmp_path / "with", "UTRAnnotator", ident) == [
            "engine\t1:200\tT\tENST1\t5UTR_consequence\tgained",
            "reference\t1:200\tT\tENST1\t5UTR_consequence\tlost",
        ]
        zero = entry(report, "LOEUF", engine=ident)
        assert (counts(zero), zero["masked"]) == ((0, 1, 0, 0, 0, 0), {"planted_engine_zero": 1})
    swept = [more for klass, _, more in handed if klass == "planted_swept_terms"]
    assert len(swept) == 4 and all(list(more) == ["utrannotator_sweep"] for more in swept)
    assert len({id(more["utrannotator_sweep"]) for more in swept}) == 1
    assert [more for klass, _, more in handed if klass == "planted_engine_zero"] == [{}, {}]


def test_summary_says_that_classes_of_utrannotator_set_nothing_aside_without_the_seeded_runs(tmp_path, capsys, monkeypatch):
    plant_sweep(monkeypatch)
    reference, engine, runs = sweep_inputs(tmp_path)
    table = class_table(tmp_path, SWEEP_TABLE)
    notice = "no --utrannotator-sweep: a class of UTRAnnotator sets aside no row it reads the seeded runs for"

    def run(name, *more):
        return run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={engine}", "--plugins", "UTRAnnotator,LOEUF",
            "--class-table", str(table), "--work-dir", str(tmp_path / f"work-{name}"), "--sort-memory", "16M",
            "--report", str(tmp_path / "report.json"), "--residuals", str(tmp_path / name), *more,
        )  # fmt: skip

    code, lines, error = run("without")
    assert (code, error) == (0, "")
    assert lines[-1] == notice
    code, lines, error = run("with", "--utrannotator-sweep", str(runs / "seed-a.txt"), "--utrannotator-sweep", str(runs / "seed-b.txt.gz"))
    assert (code, error) == (0, "")
    assert notice not in lines
    assert "UTRAnnotator\tengine\tfeature\tclass planted_swept_terms set aside 1 [met, -]" in lines
    code, lines, error = run("unloaded", "--plugins", "LOEUF")
    assert (code, error) == (0, "") and notice not in lines
    for name, path, message in (
        ("missing", tmp_path / "nowhere.txt", f"not a file: {tmp_path / 'nowhere.txt'}"),
        ("bare", runs / "tools" / "..", None),
        ("empty", runs / "tools" / "none", f"{runs / 'tools' / 'none'}: no default-format output in the directory"),
        ("refused", runs / "input.vcf", f"{runs / 'input.vcf'}: not a default-format output"),
    ):
        if name == "empty":
            path.mkdir()
            (path / "notes.txt").write_text("no output here\n", encoding="utf-8")
        code, lines, error = run(name, "--utrannotator-sweep", str(path))
        if message is None:
            assert (code, error) == (0, "") and notice not in lines
            continue
        assert (code, lines) == (2, [])
        assert error == f"ERROR: [score_plugin_tuples] {message}\n"
        assert not any((tmp_path / left).exists() for left in (name, f"work-{name}"))


UTR_STRUCTURAL = "utrannotator_structural_allele_text"
UTR_SPLICED = "5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;5UTR_annotation=alt_type=uORF:ref_type=uORF"
UTR_SPLICED_REORDERED = "5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;5UTR_annotation=ref_type=uORF:alt_type=uORF"
UTR_LOST_START = "5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;5UTR_annotation=type=uORF"


def utr_structural_row(location: str, allele: str, *fields: str) -> str:
    """A 5' UTR row under ``allele`` whose Extra carries ``fields`` after its IMPACT."""
    return row(location, allele, "ENST1", ";".join(("IMPACT=MODIFIER", *fields)), uploaded="sv")


@pytest.mark.parametrize("seeded", [False, True])
def test_utrannotator_fields_on_a_structural_allele_row_are_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch, seeded):
    """``utrannotator_structural_allele_text`` under the committed table's own row, decided by
    the comparator's predicate, unplanted: the two fields the reference alone fills on a row
    whose Allele is a structural class word are set aside, on each row of a key two rows
    share, with the seeded runs given or not. A sequence allele's row, a breakend's own
    string, a row the engine fills, a value with no engine row left to stand against, the
    plugin's ``Existing_uORFs`` and another plugin's field stay residuals, and both paths
    write the same report and residual files."""
    declared = [line for line in COMMITTED_CLASS_TABLE.read_text(encoding="utf-8").splitlines() if line.startswith(f"{UTR_STRUCTURAL}\t")]
    assert len(declared) == 1 and declared[0].split("\t")[1:4] == ["UTRAnnotator", "5UTR_consequence,5UTR_annotation", "met"]
    plugins = {"UTRAnnotator": ["5UTR_annotation", "5UTR_consequence", "Existing_uORFs"], "LOEUF": ["LOEUF"]}
    plant_registry(monkeypatch, plugins)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    existing, two_rows, one_row_less = "Existing_uORFs=6", "21:7000-7100", "21:8000-8100"
    reference = output(
        tmp_path / "reference.txt",
        utr_structural_row("21:1000-1040", "deletion", UTR_SPLICED, existing),
        utr_structural_row("21:2000-2100", "duplication", UTR_SPLICED),
        utr_structural_row("21:3000-3004", "-", UTR_SPLICED),
        utr_structural_row("21:4000-4001", "insertion"),
        utr_structural_row("21:5000-5100", "inversion", UTR_SPLICED, existing),
        utr_structural_row("21:6000-6001", "Alu_insertion", UTR_SPLICED, "LOEUF=0.31"),
        utr_structural_row(two_rows, "copy_number_variation", UTR_SPLICED),
        utr_structural_row(two_rows, "copy_number_variation", UTR_SPLICED_REORDERED),
        utr_structural_row(one_row_less, "copy_number_variation", UTR_SPLICED),
        utr_structural_row(one_row_less, "copy_number_variation", UTR_SPLICED),
        utr_structural_row("21:9000", "N[21:9500[", UTR_SPLICED),
    )
    engine = output(
        tmp_path / "engine.txt",
        utr_structural_row("21:1000-1040", "deletion", existing),
        utr_structural_row("21:2000-2100", "duplication", UTR_LOST_START),
        utr_structural_row("21:3000-3004", "-"),
        utr_structural_row("21:4000-4001", "insertion", UTR_SPLICED),
        utr_structural_row("21:5000-5100", "inversion"),
        utr_structural_row("21:6000-6001", "Alu_insertion"),
        utr_structural_row(two_rows, "copy_number_variation"),
        utr_structural_row(two_rows, "copy_number_variation"),
        utr_structural_row(one_row_less, "copy_number_variation"),
        utr_structural_row("21:9000", "N[21:9500[")
    )  # fmt: skip
    scored = {"table": declared[0] + "\n", "plugins": list(plugins), "sweep": [reference] if seeded else None}
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, **scored)
    spliced = dict(pair.split("=", 1) for pair in UTR_SPLICED.split(";"))
    lost = dict(pair.split("=", 1) for pair in UTR_LOST_START.split(";"))

    def left(side: str, location: str, allele: str, fields: dict[str, str], feature: str = "ENST1\t") -> list[str]:
        return [f"{side}\t{location}\t{allele}\t{feature}{name}\t{value}" for name, value in sorted(fields.items())]

    stays = [
        ("engine", "21:2000-2100", "duplication", lost), ("reference", "21:2000-2100", "duplication", spliced),
        ("reference", "21:3000-3004", "-", spliced), ("engine", "21:4000-4001", "insertion", spliced),
        ("reference", "21:5000-5100", "inversion", {"Existing_uORFs": "6"}),
        ("reference", one_row_less, "copy_number_variation", spliced), ("reference", "21:9000", "N[21:9500[", spliced),
    ]  # fmt: skip
    item = entry(report, "UTRAnnotator")
    assert counts(item) == (22, 5, 1, 10, 5, 1)
    assert item["masked"] == {UTR_STRUCTURAL: 10}
    assert sorted(residual_lines(residuals, "UTRAnnotator")) == sorted(line for stay in stays for line in left(*stay))
    by_allele = entry(report, "UTRAnnotator", projection="allele")
    assert by_allele["masked"] == {UTR_STRUCTURAL: 8}
    assert sorted(residual_lines(residuals, "UTRAnnotator", projection="allele")) == sorted(line for stay in stays for line in left(*stay, feature=""))
    assert entry(report, "LOEUF")["masked"] == {}
    assert residual_lines(residuals, "LOEUF") == ["reference\t21:6000-6001\tAlu_insertion\tENST1\tLOEUF\t0.31"]
    tool = built_pairing_tool()
    if tool is not None and not THROUGH_PAIRING_TOOL:
        (tmp_path / "paired").mkdir()
        paired = score(tmp_path / "paired", reference, {"engine": engine}, residuals=tmp_path / "paired-residuals", tool=tool, **scored)
        assert paired["scores"] == report["scores"]
        for written in sorted(residuals.iterdir()):
            assert (tmp_path / "paired-residuals" / written.name).read_bytes() == written.read_bytes(), written.name


UTR_SPAN_EFFECT = "5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant;5UTR_annotation=AltStop=True:FrameWithCDS=inFrame"
UTR_SPAN_EFFECT_OTHER = "5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant;5UTR_annotation=AltStop=False:FrameWithCDS=inFrame"


def utr_span_row(location: str, allele: str, cdna: str, *fields: str, terms: str = "feature_truncation,5_prime_UTR_variant") -> str:
    """A structural record's row under ``allele`` with its cDNA positions and the Consequence
    ``terms``, whose Extra carries ``fields`` after its IMPACT."""
    text = row(location, allele, "ENST1", ";".join(("IMPACT=MODIFIER", *fields)), cdna=cdna, uploaded="sv")
    assert text.count("\tmissense_variant\t") == 1
    return text.replace("\tmissense_variant\t", f"\t{terms}\t")


@pytest.mark.parametrize("seeded", [False, True])
def test_utrannotator_fields_the_engine_fills_on_a_deletion_inside_one_utr_exon_are_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch, seeded):
    """``utrannotator_structural_allele_text`` under the committed table's own row, decided by
    the comparator's predicate, unplanted: the two fields the engine fills on the row of a
    symbolic deletion inside one exon of the 5' UTR are set aside, beside an empty reference,
    beside the reference's own term and, where the term is equal, for the annotation alone.
    A row the engine fills under another class word, across two exons, into the intron, with
    cDNA positions open at an end, with the gained stop or with a frameshift on a span that
    is a multiple of three stays a residual, as do the plugin's ``Existing_uORFs`` and
    another plugin's field, and both paths write the same report and residual files."""
    declared = [line for line in COMMITTED_CLASS_TABLE.read_text(encoding="utf-8").splitlines() if line.startswith(f"{UTR_STRUCTURAL}\t")]
    assert len(declared) == 1 and declared[0].split("\t")[1:4] == ["UTRAnnotator", "5UTR_consequence,5UTR_annotation", "met"]
    plugins = {"UTRAnnotator": ["5UTR_annotation", "5UTR_consequence", "Existing_uORFs"], "LOEUF": ["LOEUF"]}
    plant_registry(monkeypatch, plugins)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    stop_gain = "5UTR_consequence=5_prime_UTR_uORF_stop_codon_gain_variant;5UTR_annotation=ref_type=uORF"
    frameshift = "5UTR_consequence=5_prime_UTR_uORF_frameshift_variant;5UTR_annotation=ref_type=uORF"
    shapes = [
        # location, allele, cDNA positions, the reference's fields, the engine's fields, Consequence terms or None
        ("21:1001-1040", "deletion", "101-140", (), (UTR_SPAN_EFFECT,), None),
        ("21:2001-2040", "deletion", "101-140", (UTR_SPLICED,), (UTR_SPAN_EFFECT,), None),
        ("21:3001-3040", "mobile_element_deletion", "101-140", (UTR_SPAN_EFFECT_OTHER,), (UTR_SPAN_EFFECT,), None),
        ("21:4001-4040", "deletion", "101-140", ("Existing_uORFs=6",), (UTR_SPAN_EFFECT, "LOEUF=0.31"), None),
        ("21:5001-5040", "duplication", "101-140", (), (UTR_SPAN_EFFECT,), None),
        ("21:6001-6040", "deletion", "101-120", (), (UTR_SPAN_EFFECT,), None),
        ("21:7001-7040", "deletion", "?-140", (UTR_SPLICED,), (UTR_SPAN_EFFECT,), None),
        ("21:8001-8040", "deletion", "101-140", (), (UTR_SPAN_EFFECT,), "feature_truncation,5_prime_UTR_variant,intron_variant"),
        ("21:9001-9040", "deletion", "101-140", (), (stop_gain,), None),
        ("21:9101-9139", "deletion", "101-139", (), (frameshift,), None),
        ("21:9201-9240", "deletion", "101-140", (), (frameshift,), None),
        ("21:9301-9340", "deletion", "101-140", (UTR_SPLICED,), (), None),
    ]
    sides = []
    for index in (3, 4):
        sides.append([
            utr_span_row(shape[0], shape[1], shape[2], *shape[index], **({} if shape[5] is None else {"terms": shape[5]}))
            for shape in shapes
        ])
    reference, engine = output(tmp_path / "reference.txt", *sides[0]), output(tmp_path / "engine.txt", *sides[1])
    scored = {"table": declared[0] + "\n", "plugins": list(plugins), "sweep": [reference] if seeded else None}
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, **scored)

    def left(side: str, location: str, allele: str, fields: tuple[str, ...]) -> list[str]:
        pairs = dict(pair.split("=", 1) for text in fields for pair in text.split(";"))
        return [f"{side}\t{location}\t{allele}\tENST1\t{name}\t{value}" for name, value in sorted(pairs.items()) if name != "LOEUF"]

    stays = [
        ("reference", "21:4001-4040", "deletion", ("Existing_uORFs=6",)),
        ("engine", "21:5001-5040", "duplication", (UTR_SPAN_EFFECT,)),
        ("engine", "21:6001-6040", "deletion", (UTR_SPAN_EFFECT,)),
        ("engine", "21:7001-7040", "deletion", (UTR_SPAN_EFFECT,)), ("reference", "21:7001-7040", "deletion", (UTR_SPLICED,)),
        ("engine", "21:8001-8040", "deletion", (UTR_SPAN_EFFECT,)),
        ("engine", "21:9001-9040", "deletion", (stop_gain,)),
        ("engine", "21:9101-9139", "deletion", (frameshift,)),
    ]  # fmt: skip
    item = entry(report, "UTRAnnotator")
    # Set aside: two fields each on the rows at 1001, 2001, 4001, 9201 and 9301, the annotation at 3001.
    assert item["masked"] == {UTR_STRUCTURAL: 11}
    assert sorted(residual_lines(residuals, "UTRAnnotator")) == sorted(line for stay in stays for line in left(*stay))
    assert residual_lines(residuals, "LOEUF") == ["engine\t21:4001-4040\tdeletion\tENST1\tLOEUF\t0.31"]
    tool = built_pairing_tool()
    if tool is not None and not THROUGH_PAIRING_TOOL:
        (tmp_path / "paired").mkdir()
        paired = score(tmp_path / "paired", reference, {"engine": engine}, residuals=tmp_path / "paired-residuals", tool=tool, **scored)
        assert paired["scores"] == report["scores"]
        for written in sorted(residuals.iterdir()):
            assert (tmp_path / "paired-residuals" / written.name).read_bytes() == written.read_bytes(), written.name


UTR_FRAMESHIFT = "5UTR_consequence=5_prime_UTR_uORF_frameshift_variant"
UTR_REFERENCE_KEYS = "5UTR_annotation=Evidence=False:KozakContext=GGGATGC:KozakStrength=Moderate:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51"
UTR_REFERENCE_KEYS_AND_TYPE = "5UTR_annotation=Evidence=False:KozakContext=GGGATGC:KozakStrength=Moderate:alt_type=uORF:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51"
UTR_TEXT_LENGTH = "5UTR_annotation=alt_type_length=39:KozakContext=GGGATGC:KozakStrength=Moderate:ref_type_length=51:ref_type=uORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=207"


@pytest.mark.parametrize("seeded", [False, True])
def test_utrannotator_fields_the_engine_fills_on_an_insertion_of_a_stated_length_are_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch, seeded):
    """``utrannotator_structural_allele_text`` under the committed table's own row, decided by
    the comparator's predicate, unplanted: on the row of a plain symbolic insertion whose
    stated length is no multiple of three, the annotation the engine writes beside the
    frameshift, of the six keys read from the reference UTR with or without ``alt_type=uORF``
    in each effect, is set aside beside the reference's own frameshift, and both fields
    beside another of the plugin's terms and beside an empty reference, whatever other
    Consequence terms the row names. An effect with ``alt_type_length``, another
    ``alt_type``, five of the six keys or an eighth key, another term, the frameshift beside
    no annotation, a span that is a multiple of three, a mobile element's word, a
    duplication, an inversion and a row without ``5_prime_UTR_variant`` stay residuals, and
    both paths write the same report and residual files."""
    declared = [line for line in COMMITTED_CLASS_TABLE.read_text(encoding="utf-8").splitlines() if line.startswith(f"{UTR_STRUCTURAL}\t")]
    assert len(declared) == 1 and declared[0].split("\t")[1:4] == ["UTRAnnotator", "5UTR_consequence,5UTR_annotation", "met"]
    plugins = {"UTRAnnotator": ["5UTR_annotation", "5UTR_consequence", "Existing_uORFs"]}
    plant_registry(monkeypatch, plugins)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    text_length, six, six_and_type = (UTR_FRAMESHIFT, UTR_TEXT_LENGTH), (UTR_FRAMESHIFT, UTR_REFERENCE_KEYS), (UTR_FRAMESHIFT, UTR_REFERENCE_KEYS_AND_TYPE)
    two_effects = f"{UTR_REFERENCE_KEYS},{UTR_REFERENCE_KEYS_AND_TYPE.partition('=')[2]}"
    stop_loss = "5UTR_consequence=5_prime_UTR_uORF_stop_codon_loss_variant"
    other_effects = {
        "21:6001-6113": f"{UTR_REFERENCE_KEYS_AND_TYPE}:alt_type_length=66",
        "21:6201-6313": UTR_REFERENCE_KEYS.replace(":ref_StartDistanceToCDS", ":alt_type=OutOfFrame_oORF:ref_StartDistanceToCDS"),
        "21:6401-6513": UTR_REFERENCE_KEYS.replace("Evidence=False:", ""),
        "21:6601-6713": f"{UTR_REFERENCE_KEYS_AND_TYPE}:type=uORF",
    }
    assert not {UTR_REFERENCE_KEYS, UTR_REFERENCE_KEYS_AND_TYPE} & set(other_effects.values())
    other_words = {"21:9001-9113": "Alu_insertion", "21:9201-9313": "mobile_element_insertion", "21:9401-9513": "duplication", "21:9601-9713": "inversion"}
    nmd, intron, no_utr = (f"feature_elongation,{terms}" for terms in ("5_prime_UTR_variant,NMD_transcript_variant", "5_prime_UTR_variant,intron_variant", "intron_variant"))
    shapes = [
        # location, allele, the reference's fields, the engine's fields, Consequence terms or None
        ("21:1001-1113", "insertion", text_length, six, None),
        ("21:2001-2113", "insertion", text_length, (UTR_FRAMESHIFT, two_effects), None),
        ("21:3001-3113", "insertion", (UTR_LOST_START,), six, None),
        ("21:4001-4113", "insertion", (), six_and_type, None),
        ("21:5001-5113", "insertion", text_length, six, nmd),
        ("21:5201-5313", "insertion", text_length, six, intron),
        *((location, "insertion", text_length, (UTR_FRAMESHIFT, effect), None) for location, effect in other_effects.items()),
        ("21:7001-7113", "insertion", text_length, (stop_loss, UTR_REFERENCE_KEYS), None),
        ("21:7201-7313", "insertion", text_length, (UTR_FRAMESHIFT,), None),
        ("21:8001-8114", "insertion", text_length, six, None),
        *((location, word, text_length, six, None) for location, word in other_words.items()),
        ("21:9801-9913", "insertion", text_length, six, no_utr),
    ]
    sides = []
    for index in (2, 3):
        sides.append([
            utr_span_row(shape[0], shape[1], "101-102", *shape[index], terms="feature_elongation,5_prime_UTR_variant" if shape[4] is None else shape[4])
            for shape in shapes
        ])
    reference, engine = output(tmp_path / "reference.txt", *sides[0]), output(tmp_path / "engine.txt", *sides[1])
    scored = {"table": declared[0] + "\n", "plugins": list(plugins), "sweep": [reference] if seeded else None}
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, **scored)

    def left(side: str, location: str, allele: str, *fields: str) -> list[str]:
        pairs = dict(pair.split("=", 1) for text in fields for pair in text.split(";"))
        return [f"{side}\t{location}\t{allele}\tENST1\t{name}\t{value}" for name, value in sorted(pairs.items())]

    # Beside the reference's frameshift the consequence is shared, and the annotations stay.
    stays = [
        *(line for location, effect in other_effects.items() for line in (("reference", location, "insertion", UTR_TEXT_LENGTH), ("engine", location, "insertion", effect))),
        ("reference", "21:7001-7113", "insertion", *text_length), ("engine", "21:7001-7113", "insertion", stop_loss, UTR_REFERENCE_KEYS),
        ("reference", "21:7201-7313", "insertion", UTR_TEXT_LENGTH),
        ("reference", "21:8001-8114", "insertion", UTR_TEXT_LENGTH), ("engine", "21:8001-8114", "insertion", UTR_REFERENCE_KEYS),
        *(line for location, word in other_words.items() for line in (("reference", location, word, UTR_TEXT_LENGTH), ("engine", location, word, UTR_REFERENCE_KEYS))),
        ("reference", "21:9801-9913", "insertion", UTR_TEXT_LENGTH), ("engine", "21:9801-9913", "insertion", UTR_REFERENCE_KEYS),
    ]  # fmt: skip
    item = entry(report, "UTRAnnotator")
    # Set aside: the annotation on the rows at 1001, 2001, 5001 and 5201, two fields each at 3001 and 4001.
    assert item["masked"] == {UTR_STRUCTURAL: 8}
    assert sorted(residual_lines(residuals, "UTRAnnotator")) == sorted(line for stay in stays for line in left(*stay))
    tool = built_pairing_tool()
    if tool is not None and not THROUGH_PAIRING_TOOL:
        (tmp_path / "paired").mkdir()
        paired = score(tmp_path / "paired", reference, {"engine": engine}, residuals=tmp_path / "paired-residuals", tool=tool, **scored)
        assert paired["scores"] == report["scores"]
        for written in sorted(residuals.iterdir()):
            assert (tmp_path / "paired-residuals" / written.name).read_bytes() == written.read_bytes(), written.name


UTR_START_LOST = "5UTR_consequence=5_prime_UTR_premature_start_codon_loss_variant;5UTR_annotation=CapDistanceToStart=63:type=uORF"
UTR_START_GAINED = "5UTR_consequence=5_prime_UTR_premature_start_codon_gain_variant;5UTR_annotation=CapDistanceToStart=64:type=uORF"


def test_row_under_an_identifier_several_input_lines_carry_is_handed_the_record_that_writes_it(tmp_path, monkeypatch):
    """``utrannotator_gained_start_window`` under the committed table's own row, decided by the
    comparator's predicate, unplanted: an insertion's row under an identifier is handed the
    one input line that writes it at its Location with its allele, where several lines carry
    the identifier, and its two fields are set aside as where one line alone carries it.
    Two lines that both write the row, lines of which none does and a run without ``--input``
    leave both fields residuals, and two records that share the identifier and the Location
    are each set aside by their own line."""
    klass = "utrannotator_gained_start_window"
    declared = [line for line in COMMITTED_CLASS_TABLE.read_text(encoding="utf-8").splitlines() if line.startswith(f"{klass}\t")]
    assert len(declared) == 1
    plant_registry(monkeypatch, {"UTRAnnotator": ["5UTR_annotation", "5UTR_consequence"]})
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    scored = {"table": declared[0] + "\n", "plugins": ["UTRAnnotator"]}
    own, other_allele, elsewhere = ("1", "100", "rs1", "C", "CCAT"), ("1", "100", "rs1", "C", "CCAG"), ("1", "150", "rs1", "C", "CCAT")
    grid = {
        "one line": ((own,), 2),
        "two lines at one place with different alleles": ((other_allele, own), 2),
        "two lines at different places": ((elsewhere, own), 2),
        "two lines that both write the row": ((own, own), 0),
        "a line whose two alternates write the row in two ways, beside the insertion's line": ((("1", "100", "rs1", "CA", "CAT,CCATA"), own), 0),
        "two lines of which none writes the row": ((other_allele, elsewhere), 0),
        "no input": (None, 0),
    }
    reference = output(tmp_path / "reference.txt", row("1:100-101", "CAT", "ENST1", UTR_START_LOST, cdna="64-65", uploaded="rs1"))
    engine = output(tmp_path / "engine.txt", row("1:100-101", "CAT", "ENST1", UTR_START_GAINED, cdna="64-65", uploaded="rs1"))
    for number, (what, (lines, keys)) in enumerate(grid.items()):
        work = tmp_path / f"grid{number}"
        work.mkdir()
        run_input = None if lines is None else vcf(work / "input.vcf", *lines)
        item = entry(score(work, reference, {"engine": engine}, residuals=work / "residuals", run_input=run_input, **scored), "UTRAnnotator")
        assert (item["masked"], item["residual_rows"]) == ({klass: keys}, 4 - 2 * keys), what
    # Two records under one identifier at one Location: each row reads its own line.
    reference = output(
        tmp_path / "reference.txt",
        row("1:100-101", "CAT", "ENST1", UTR_START_LOST, cdna="64-65", uploaded="rs1"),
        row("1:100-101", "CAG", "ENST1", UTR_START_LOST, cdna="64-65", uploaded="rs1"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100-101", "CAT", "ENST1", UTR_START_GAINED, cdna="64-65", uploaded="rs1"),
        row("1:100-101", "CAG", "ENST1", UTR_START_GAINED, cdna="64-65", uploaded="rs1"),
    )
    work = tmp_path / "two-records"
    work.mkdir()
    report = score(work, reference, {"engine": engine}, residuals=work / "residuals", run_input=vcf(work / "input.vcf", own, other_allele), **scored)
    item = entry(report, "UTRAnnotator")
    assert (counts(item), item["masked"]) == ((4, 4, 0, 0, 0, 0), {klass: 4})
    assert residual_lines(work / "residuals", "UTRAnnotator") == []
    tool = built_pairing_tool()
    if tool is not None and not THROUGH_PAIRING_TOOL:
        paired = tmp_path / "paired"
        paired.mkdir()
        again = score(paired, reference, {"engine": engine}, residuals=paired / "residuals", run_input=vcf(paired / "input.vcf", own, other_allele), tool=tool, **scored)
        assert again["scores"] == report["scores"]


def test_row_at_a_two_base_location_under_an_identifier_is_handed_its_input_line(tmp_path, monkeypatch):
    plant_record_reader(monkeypatch)
    handed = []

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None):
        handed.append(record_rows)
        return any((r.fields.get("REF"), r.fields.get("ALT")) == ("C", "CT") for r in record_rows or ())

    monkeypatch.setattr(formats, "FIELD_CLASSES", {"planted_inserted_length": ()})
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    table = "planted_inserted_length\tUTRAnnotator\t5UTR_consequence\n"
    reference = output(tmp_path / "reference.txt", row("1:100-101", "T", "ENST1", "5UTR_consequence=lost", cdna="64-65", uploaded="ins1"))
    engine = output(tmp_path / "engine.txt", row("1:100-101", "T", "ENST1", "5UTR_consequence=gained", cdna="64-65", uploaded="ins1"))
    run_input = vcf(tmp_path / "input.vcf", ("1", "100", "ins1", "C", "CT"), ("1", "100", "other", "C", "G"))
    plugins = ["UTRAnnotator"]
    without = entry(score(tmp_path, reference, {"engine": engine}, table=table, plugins=plugins), "UTRAnnotator")
    assert (counts(without), without["masked"], handed) == ((1, 1, 0, 1, 1, 0), {"planted_inserted_length": 0}, [None])
    item = entry(score(tmp_path, reference, {"engine": engine}, table=table, plugins=plugins, run_input=run_input), "UTRAnnotator")
    assert (counts(item), item["masked"]) == ((1, 1, 0, 0, 0, 0), {"planted_inserted_length": 1})
    row_entry, line_entry = handed[1]
    assert (row_entry.location, row_entry.allele, row_entry.fields["cDNA_position"]) == ("1:100-101", "T", "64-65")
    assert line_entry.fields == {"CHROM": "1", "POS": "100", "Uploaded_variation": "ins1", "REF": "C", "ALT": "CT"}


DOWNSTREAM_TABLE = "planted_recomputed_value\tDownstream\tDownstreamProtein\tmet\t-\nplanted_engine_zero\tLOEUF\tLOEUF\tmet\t-\n"


def plant_downstream(monkeypatch) -> tuple[list, list]:
    """Plants a reader of a candidate cache and a class that accepts a pair when it is handed
    that reader, and returns the caches read and, per predicate call, the class and the
    keywords handed beyond the rows."""
    reads, handed = [], []

    class Sequences:
        def __init__(self, cache):
            if not (cache / "transcripts").is_dir():
                raise formats.ComparatorError(f"--json-cache: {cache} has no transcripts/ directory")
            reads.append(cache)

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None, **more):
        handed.append((klass, more))
        if klass == "planted_recomputed_value":
            return "downstream_sequences" in more
        return klass == "planted_engine_zero" and formats.blank(expected_value) and actual_value == "0"

    plant_registry(monkeypatch, {"Downstream": ["DownstreamProtein"], "LOEUF": ["LOEUF"]})
    monkeypatch.setattr(formats, "DownstreamSequences", Sequences, raising=False)
    monkeypatch.setattr(formats, "FIELD_CLASSES", dict.fromkeys((*PLANTED_CLASSES, "planted_recomputed_value"), ()))
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    return reads, handed


def downstream_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """A reference and an engine that write different values of a Downstream field at one
    Location, and a candidate cache directory."""
    reference = output(tmp_path / "reference.txt", row("1:100", "-", "ENST1", "DownstreamProtein=KL"), row("1:300", "T", "ENST1", "IMPACT=LOW"))
    engine = output(tmp_path / "engine.txt", row("1:100", "-", "ENST1", "DownstreamProtein=RKL"), row("1:300", "T", "ENST1", "LOEUF=0"))
    cache = tmp_path / "cache"
    (cache / "transcripts").mkdir(parents=True)
    return reference, engine, cache


def test_candidate_cache_is_read_once_and_handed_to_the_classes_of_downstream_alone(tmp_path, monkeypatch):
    reads, handed = plant_downstream(monkeypatch)
    reference, engine, cache = downstream_inputs(tmp_path)
    engines = {"engine": engine, "again": engine}
    plugins = ["Downstream", "LOEUF"]
    without = score(tmp_path, reference, engines, table=DOWNSTREAM_TABLE, plugins=plugins)
    assert without["json_cache"] is None and reads == []
    assert handed and all(more == {} for _, more in handed)
    assert counts(entry(without, "Downstream")) == (1, 1, 0, 1, 1, 0)
    assert entry(without, "Downstream")["masked"] == {"planted_recomputed_value": 0}
    handed.clear()
    report = score(tmp_path, reference, engines, table=DOWNSTREAM_TABLE, plugins=plugins, cache=cache)
    assert report["json_cache"] == {"path": str(cache)}
    assert reads == [cache]
    for ident in engines:
        item = entry(report, "Downstream", engine=ident)
        assert (counts(item), item["masked"]) == ((1, 1, 0, 0, 0, 0), {"planted_recomputed_value": 1})
        zero = entry(report, "LOEUF", engine=ident)
        assert (counts(zero), zero["masked"]) == ((0, 1, 0, 0, 0, 0), {"planted_engine_zero": 1})
    of_downstream = [more for klass, more in handed if klass == "planted_recomputed_value"]
    assert of_downstream and all(list(more) == ["downstream_sequences", "declared"] for more in of_downstream)
    assert len({id(more["downstream_sequences"]) for more in of_downstream}) == 1
    assert all(more == {} for klass, more in handed if klass == "planted_engine_zero")


def test_a_class_of_downstream_is_handed_the_classes_the_run_applies(tmp_path, monkeypatch):
    """A Downstream value two causes change is set aside only beside the second cause's class
    (``downstream_class`` of ``compare_vep_formats.py``), so a class of the plugin is handed
    the applied classes as ``declared``: every declared one by default, the named ones under
    ``--classes``."""
    _, handed = plant_downstream(monkeypatch)
    reference, engine, cache = downstream_inputs(tmp_path)
    plugins = ["Downstream", "LOEUF"]

    def declared_sets() -> set[frozenset[str]]:
        return {frozenset(more["declared"]) for klass, more in handed if klass == "planted_recomputed_value"}

    score(tmp_path, reference, {"engine": engine}, table=DOWNSTREAM_TABLE, plugins=plugins, cache=cache)
    assert declared_sets() == {frozenset({"planted_recomputed_value", "planted_engine_zero"})}
    handed.clear()
    score(tmp_path, reference, {"engine": engine}, table=DOWNSTREAM_TABLE, plugins=plugins, cache=cache, applied=["planted_recomputed_value"])
    assert declared_sets() == {frozenset({"planted_recomputed_value"})}


def test_summary_says_that_classes_of_downstream_set_nothing_aside_without_the_candidate_cache(tmp_path, capsys, monkeypatch):
    plant_downstream(monkeypatch)
    reference, engine, cache = downstream_inputs(tmp_path)
    table = class_table(tmp_path, DOWNSTREAM_TABLE)
    notice = "no --json-cache: a class of Downstream that reads the candidate cache sets nothing aside"

    def run(name, *more):
        return run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={engine}", "--plugins", "Downstream,LOEUF",
            "--class-table", str(table), "--work-dir", str(tmp_path / f"work-{name}"), "--sort-memory", "16M",
            "--report", str(tmp_path / "report.json"), "--residuals", str(tmp_path / name), *more,
        )  # fmt: skip

    code, lines, error = run("without")
    assert (code, error) == (0, "")
    assert lines[-1] == notice
    code, lines, error = run("with", "--json-cache", str(cache))
    assert (code, error) == (0, "")
    assert notice not in lines
    assert "Downstream\tengine\tfeature\tclass planted_recomputed_value set aside 1 [met, -]" in lines
    for name, path, message in (
        ("missing", tmp_path / "nowhere", f"not a directory: {tmp_path / 'nowhere'}"),
        ("refused", tmp_path, f"--json-cache: {tmp_path} has no transcripts/ directory"),
    ):
        code, lines, error = run(name, "--json-cache", str(path))
        assert (code, lines) == (2, [])
        assert error == f"ERROR: [score_plugin_tuples] {message}\n"
        assert not any((tmp_path / left).exists() for left in (name, f"work-{name}"))


EMPTY_VALUE_TABLE = "".join(
    f"{klass}\tDownstream\tDownstreamProtein\tmet\t-\n" for klass in ("planted_other_class", *formats.DOWNSTREAM_CLASSES)
)


def test_classes_of_downstream_are_handed_an_empty_value_as_a_value(tmp_path, monkeypatch):
    """A row that writes ``DownstreamProtein=`` writes no tuple, so the other side's value is
    left with none to pair. The two classes of Downstream read the empty value as the
    plugin's own, and are handed it as the empty text; a ``-``, a row without the key and a
    tuple file's missing tuple stay absent, and any other class is handed absent for all."""
    handed = []
    first, second = formats.DOWNSTREAM_CLASSES

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None, **more):
        handed.append((klass, expected_value, actual_value))
        if {expected_value, actual_value} & {"-", formats.ABSENT}:
            return False
        return (klass, "") in ((first, expected_value), (second, actual_value))

    plant_registry(monkeypatch, {"Downstream": ["DownstreamProtein"]})
    monkeypatch.setattr(formats, "FIELD_CLASSES", dict.fromkeys(("planted_other_class", first, second), ()))
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    key = ("-", "ENST1")
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", *key, "DownstreamProtein="),
        row("1:200", *key, "DownstreamProtein=KL"),
        row("1:300", *key, "DownstreamProtein=KL"),
        row("1:400", *key, "DownstreamProtein=-"),
        row("1:500", *key, "DownstreamProtein="),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", *key, "DownstreamProtein=K"),
        row("1:200", *key, "DownstreamProtein="),
        row("1:300", *key, "IMPACT=LOW"),
        row("1:400", *key, "DownstreamProtein=R"),
        row("1:500", *key, "DownstreamProtein="),
    )
    listed = tuples(tmp_path / "engine.tsv", ("1:100", *key, "DownstreamProtein", "K"), ("1:400", *key, "DownstreamProtein", "R"))
    residuals = tmp_path / "residuals"
    report = score(
        tmp_path, reference, {"engine": engine, "listed": listed}, residuals=residuals, table=EMPTY_VALUE_TABLE, plugins=["Downstream"]
    )
    left = ["reference\t1:300\t-\tENST1\tDownstreamProtein\tKL", "engine\t1:400\t-\tENST1\tDownstreamProtein\tR"]
    item = entry(report, "Downstream")
    assert counts(item) == (2, 2, 0, 1, 1, 0)
    assert item["masked"] == {"planted_other_class": 0, first: 1, second: 1}
    assert residual_lines(residuals, "Downstream") == left
    # A tuple file states no empty value: its side of the second row stays absent.
    item = entry(report, "Downstream", engine="listed")
    assert counts(item) == (2, 2, 0, 2, 1, 0)
    assert item["masked"] == {"planted_other_class": 0, first: 1, second: 0}
    assert residual_lines(residuals, "Downstream", "listed") == ["reference\t1:200\t-\tENST1\tDownstreamProtein\tKL", *left]
    of_downstream = {(expected, actual) for klass, expected, actual in handed if klass in (first, second)}
    assert of_downstream == {("", "K"), ("KL", ""), ("KL", formats.ABSENT), (formats.ABSENT, "R")}
    of_other = {(expected, actual) for klass, expected, actual in handed if klass == "planted_other_class"}
    assert of_other == {(formats.ABSENT, "K"), ("KL", formats.ABSENT), (formats.ABSENT, "R")}


VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def test_the_class_of_the_shifted_cds_coordinates_is_handed_an_empty_value_as_a_value(tmp_path, monkeypatch):
    """``plugin_hgvs_shifted_cds_coordinates`` reads a ``DownstreamProtein`` as the two classes of
    Downstream do, so a side whose row writes the field with an empty value is handed to it as
    the empty text, on either side; a class the comparator does not count among
    ``TRANSCRIPT_READING_CLASSES`` is handed absent for the same rows."""
    handed = []
    shifted = formats.PLUGIN_HGVS_SHIFTED_CDS_COORDINATES
    assert shifted in formats.TRANSCRIPT_READING_CLASSES and shifted not in formats.DOWNSTREAM_CLASSES

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None, **more):
        handed.append((klass, expected_value, actual_value))
        return klass == shifted and "" in (expected_value, actual_value) and formats.ABSENT not in (expected_value, actual_value)

    plant_registry(monkeypatch, {"Downstream": ["DownstreamProtein"]})
    monkeypatch.setattr(formats, "FIELD_CLASSES", dict.fromkeys(("planted_other_class", shifted), ()))
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    key = ("-", "ENST1")
    reference = output(tmp_path / "reference.txt", row("1:100", *key, "DownstreamProtein="), row("1:200", *key, "DownstreamProtein=KL"))
    engine = output(tmp_path / "engine.txt", row("1:100", *key, "DownstreamProtein=K"), row("1:200", *key, "DownstreamProtein="))
    table = f"planted_other_class\tDownstream\tDownstreamProtein\tmet\t-\n{shifted}\tDownstream\tDownstreamProtein\tmet\t-\n"
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=table, plugins=["Downstream"])
    item = entry(report, "Downstream")
    assert handed == [
        ("planted_other_class", formats.ABSENT, "K"),
        (shifted, "", "K"),
        ("planted_other_class", "KL", formats.ABSENT),
        (shifted, "KL", ""),
    ]
    assert item["masked"] == {"planted_other_class": 0, shifted: 2}
    assert counts(item) == (1, 1, 0, 0, 0, 0)
    assert residual_lines(residuals, "Downstream") == []


def test_a_later_alleles_downstream_value_is_set_aside_where_the_input_states_the_records_alleles(tmp_path, monkeypatch):
    """The format comparison's own rule on a made-up record of two alternates, one base or two
    inserted after coding base 4 of a transcript whose next three bases repeat it: the later
    allele's reference row prints the CDS pair three bases on, the candidate's the record's
    own, and ``plugin_hgvs_shifted_cds_coordinates`` takes the row's value only as a later
    allele's. The rule is handed the record's rows with the input's line, which states the
    alleles in order, so both alleles' values are set aside; without ``--input`` the later
    allele's stays."""
    shifted = formats.PLUGIN_HGVS_SHIFTED_CDS_COORDINATES
    plant_registry(monkeypatch, {"Downstream": ["DownstreamProtein"]})
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    coding = "ATGAAAACTCTTGCTGTCCGTGGTTAA"
    model = {"pairs": [{"from_start": 1, "from_end": 42, "to_start": 1001, "to_end": 1042, "ori": 1}]}
    transcript = {
        "stable_id": "ENST1", "strand": 1, "cdna_coding_start": 1, "cdna_coding_end": len(coding), "exons": [{"start": 1001, "end": 1042}],
        "variation_effect_feature_cache": {"translateable_seq": coding, "three_prime_utr": "GGGCCCAAATAGCCC", "codon_table": 1, "mapper": model},
    }  # fmt: skip
    cache = tmp_path / "cache"
    (cache / "transcripts" / "21").mkdir(parents=True)
    (cache / "transcripts" / "21" / "1-1000000.json").write_text(json.dumps([transcript]), encoding="utf-8")

    def frameshift(allele: str, cds: str, protein: str) -> str:
        columns = ["v", "21:1004-1005", allele, "ENSG1", "ENST1", "Transcript", "frameshift_variant", "-", cds, "-", "-", "-", "-"]
        return "\t".join([*columns, f"IMPACT=HIGH;STRAND=1;HGVS_OFFSET=3;DownstreamProtein={protein};ProteinLengthChange=6"]) + "\n"

    reference = output(tmp_path / "reference.txt", frameshift("A", "4-5", "NSCCPWLRAQIA"), frameshift("AA", "7-8", "KLLLSVVKGPNS"))
    engine = output(tmp_path / "engine.txt", frameshift("A", "4-5", "KNSCCPWLRAQIA"), frameshift("AA", "4-5", "KKLLLSVVKGPNS"))
    table = f"{shifted}\tDownstream\tDownstreamProtein\tmet\t-\n"
    scored = {"table": table, "plugins": ["Downstream"], "cache": cache}
    without = entry(score(tmp_path, reference, {"engine": engine}, **scored), "Downstream")
    assert without["masked"] == {shifted: 1}
    residuals = tmp_path / "residuals"
    run_input = vcf(tmp_path / "input.vcf", ("21", "1004", "v", "A", "AA,AAA"))
    item = entry(score(tmp_path, reference, {"engine": engine}, residuals=residuals, run_input=run_input, **scored), "Downstream")
    assert item["masked"] == {shifted: 2}
    assert residual_lines(residuals, "Downstream") == []


def coding_row(location: str, allele: str, terms: str, cdna: str, cds: str, extra: str, uploaded: str = "v") -> str:
    columns = [uploaded, location, allele, "ENSG1", "ENST1", "Transcript", terms]
    return "\t".join([*columns, cdna, cds, "-", "-", "-", "-", extra]) + "\n"


def test_plugin_value_that_follows_a_stop_codon_term_pair_is_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch):
    """``plugin_value_follows_stop_codon_term_pair`` is declared like any class of
    ``FIELD_CLASSES`` and decided by the comparator's own predicate, unplanted: a Downstream,
    NMD or LOFTEE value one side alone writes is set aside on a row whose two term sets are a
    stop codon pair the class reads (the untested stop codon, the insertion that keeps the
    stop codon, the three-base test), on whichever side's terms let the plugin run. A
    ``LoF_flags`` the table's row leaves out and the same values on a row of no stop codon
    pair stay residuals."""
    klass = formats.PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR
    plugins = {"Downstream": ["DownstreamProtein", "ProteinLengthChange"], "NMD": ["NMD"],
               "LoFTEE": ["LoF", "LoF_filter", "LoF_flags", "LoF_info"]}
    plant_registry(monkeypatch, plugins)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    untested = ("19:13208046-13208050", "-"), ("1705-?", "1630-?")
    kept_stop = ("19:13208146-13208147", "TTGA"), ("7054-7055", "6977-6978")
    tested = ("2:190970701-190970707", "-"), ("961-967", "264-270")
    no_pair = ("7:500-502", "-"), ("100-102", "40-42")
    escaping = "NMD=NMD_escaping_variant"
    reference = output(
        tmp_path / "reference.txt",
        coding_row(*untested[0], "stop_retained_variant", *untested[1], "IMPACT=LOW"),
        coding_row(*kept_stop[0], "stop_gained,frameshift_variant", *kept_stop[1],
                   f"IMPACT=HIGH;DownstreamProtein=MK;ProteinLengthChange=-3;{escaping};LoF=LC;LoF_info=PERCENTILE:0.5"),
        coding_row(*tested[0], "stop_retained_variant", *tested[1], "IMPACT=LOW"),
        coding_row(*no_pair[0], "inframe_deletion", *no_pair[1], "IMPACT=MODERATE"),
    )
    engine = output(
        tmp_path / "engine.txt",
        coding_row(*untested[0], "frameshift_variant", *untested[1],
                   f"IMPACT=HIGH;DownstreamProtein=MKV;ProteinLengthChange=-12;{escaping};LoF=HC;LoF_filter=END_TRUNC;"
                   "LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.97"),
        coding_row(*kept_stop[0], "inframe_insertion,stop_retained_variant", *kept_stop[1], "IMPACT=MODERATE"),
        coding_row(*tested[0], "frameshift_variant", *tested[1], f"IMPACT=HIGH;{escaping}"),
        coding_row(*no_pair[0], "frameshift_variant", *no_pair[1], f"IMPACT=HIGH;DownstreamProtein=MR;{escaping}"),
    )
    table = f"{klass}\tDownstream;NMD;LoFTEE\tDownstreamProtein,ProteinLengthChange;NMD;LoF,LoF_info,LoF_filter\tmet\t-\n"
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=table, plugins=list(plugins))
    expected = {
        "Downstream": ((2, 3, 0, 0, 1, 0), 4, ["engine\t7:500-502\t-\tENST1\tDownstreamProtein\tMR"]),
        "NMD": ((1, 3, 0, 0, 1, 0), 3, ["engine\t7:500-502\t-\tENST1\tNMD\tNMD_escaping_variant"]),
        "LoFTEE": ((2, 4, 0, 0, 1, 0), 5, ["engine\t19:13208046-13208050\t-\tENST1\tLoF_flags\tSINGLE_EXON"]),
    }
    for source, (counted, keys, left) in expected.items():
        item = entry(report, source)
        assert counts(item) == counted, source
        assert item["masked"] == {klass: keys}, source
        assert residual_lines(residuals, source) == left, source


def test_loftee_flags_that_follow_a_stop_codon_term_pair_are_set_aside_where_the_class_table_names_the_field(tmp_path, monkeypatch):
    """``LoF_flags`` is a field ``plugin_value_follows_stop_codon_term_pair`` reads, decided
    by the comparator's own predicate, unplanted: under a class table whose row names the
    field, a flag list one side alone writes is set aside on a row whose two term sets are a
    stop codon pair the class reads, on whichever side's terms let LOFTEE run, and the same
    list on a row of no stop codon pair stays a residual. Under a row that leaves the field
    out every flag list stays a residual."""
    klass = formats.PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR
    plugins = {"LoFTEE": ["LoF", "LoF_filter", "LoF_flags", "LoF_info"]}
    plant_registry(monkeypatch, plugins)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    untested = ("19:13208046-13208050", "-"), ("1705-?", "1630-?")
    kept_stop = ("19:13208146-13208147", "TTGA"), ("7054-7055", "6977-6978")
    no_pair = ("7:500-502", "-"), ("100-102", "40-42")
    reference = output(
        tmp_path / "reference.txt",
        coding_row(*untested[0], "stop_retained_variant", *untested[1], "IMPACT=LOW"),
        coding_row(*kept_stop[0], "stop_gained,frameshift_variant", *kept_stop[1], "IMPACT=HIGH;LoF=HC;LoF_flags=SINGLE_EXON"),
        coding_row(*no_pair[0], "inframe_deletion", *no_pair[1], "IMPACT=MODERATE"),
    )
    engine = output(
        tmp_path / "engine.txt",
        coding_row(*untested[0], "frameshift_variant", *untested[1], "IMPACT=HIGH;LoF=HC;LoF_flags=PHYLOCSF_WEAK"),
        coding_row(*kept_stop[0], "inframe_insertion,stop_retained_variant", *kept_stop[1], "IMPACT=MODERATE"),
        coding_row(*no_pair[0], "frameshift_variant", *no_pair[1], "IMPACT=HIGH;LoF=HC;LoF_flags=SINGLE_EXON"),
    )
    off_pair = ["engine\t7:500-502\t-\tENST1\tLoF\tHC", "engine\t7:500-502\t-\tENST1\tLoF_flags\tSINGLE_EXON"]
    on_pair = [
        "engine\t19:13208046-13208050\t-\tENST1\tLoF_flags\tPHYLOCSF_WEAK",
        "reference\t19:13208146-13208147\tTTGA\tENST1\tLoF_flags\tSINGLE_EXON",
    ]
    expected = {
        "LoF,LoF_info,LoF_filter,LoF_flags": ((2, 4, 0, 0, 2, 0), 4, off_pair),
        "LoF,LoF_info,LoF_filter": ((2, 4, 0, 1, 3, 0), 2, on_pair + off_pair),
    }
    for number, (fields, (counted, keys, left)) in enumerate(expected.items()):
        run = tmp_path / str(number)
        run.mkdir()
        report = score(run, reference, {"engine": engine}, residuals=run / "residuals",
                       table=f"{klass}\tLoFTEE\t{fields}\tmet\t-\n", plugins=list(plugins))
        item = entry(report, "LoFTEE")
        assert counts(item) == counted, fields
        assert item["masked"] == {klass: keys}, fields
        assert residual_lines(run / "residuals", "LoFTEE") == left, fields


SPLICE_SPAN_PLUGINS = {"Downstream": ["DownstreamProtein", "ProteinLengthChange"], "NMD": ["NMD"],
                       "LoFTEE": ["LoF", "LoF_filter", "LoF_flags", "LoF_info"]}
SPLICE_SPAN_TABLE = (
    f"{formats.PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN}\tDownstream;NMD;LoFTEE\t"
    "DownstreamProtein,ProteinLengthChange;NMD;LoF,LoF_info,LoF_filter,LoF_flags\tmet\t-\n"
)
# Per plugin, what the class leaves of the outputs of ``splice_span_outputs`` scored with their
# input: the counts, the keys set aside and the residual lines.
SPLICE_SPAN_SCORES = {
    "Downstream": ((0, 4, 0, 0, 2, 0), 2, [
        "engine\t7:300999-301000\tACGTA\tENST1\tDownstreamProtein\tMKL",
        "engine\t7:300999-301000\tACGTA\tENST1\tProteinLengthChange\t-3",
    ]),
    "NMD": ((6, 2, 2, 3, 2, 2), 3, ["reference\t7:400998-401000\tACGTAC\tENST1\tNMD\tNMD_escaping_variant"]),
    "LoFTEE": ((11, 3, 2, 4, 3, 2), 7, [
        "reference\t7:400998-401000\tACGTAC\tENST1\tLoF\tHC",
        "engine\t7:500998-501000\tACGTAC\tENST1\tLoF_info\tPERCENTILE:0.5",
        "reference\t7:500998-501000\tACGTAC\tENST1\tLoF_info\tPERCENTILE:0.9",
    ]),
}  # fmt: skip
# The same scored without the input, where the row under the identifier ``rs90`` is read
# without its record: its NMD and LOFTEE values stay residuals.
SPLICE_SPAN_SCORES_WITHOUT_INPUT = {
    "Downstream": SPLICE_SPAN_SCORES["Downstream"],
    "NMD": ((6, 2, 2, 4, 2, 2), 2, [*SPLICE_SPAN_SCORES["NMD"][2], "reference\t7:600999-601000\tCGTAC\tENST1\tNMD\tNMD_escaping_variant"]),
    "LoFTEE": ((11, 3, 2, 5, 3, 2), 6, [*SPLICE_SPAN_SCORES["LoFTEE"][2], "reference\t7:600999-601000\tCGTAC\tENST1\tLoF\tHC"]),
}  # fmt: skip


def splice_span_outputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    """A made-up reference and engine output of seven rows, and the input VCF of the run.
    Three are rows of ``splice_term_past_the_record_span`` read without their record, three
    positions in the transcript's sequence under a longer allele and the identifier ``v``,
    which no line of the input carries: one on which the reference alone holds the donor
    term, its row's only truncating term, and writes NMD and four LOFTEE fields; one on which
    it holds the acceptor term beside frameshift_variant, and the engine, whose row holds no
    splice term, alone writes the two Downstream fields; and one on which both sides write a
    ``LoF_info``, with two values. The values of the first two stand as well, under ``v``, on
    a row of two neighbouring positions and on a row whose span crosses an exon's end. Two
    more rows of two positions inside an exon are the class's with their record alone, the
    reference writing NMD and ``LoF`` beside its donor term on each: one under the identifier
    ``rs90``, whose record the input holds, and one whose first column spells the alleles."""
    in_frame = "splice_donor_variant,protein_altering_variant", "protein_altering_variant,splice_region_variant"
    shifted = "splice_acceptor_variant,frameshift_variant", "frameshift_variant"
    taken = ("7:100998-101000", "ACGTAC"), ("998-1000", "898-900")
    taken_shifted = ("7:200998-201000", "ACGTA"), ("498-500", "398-400")
    two_neighbours = ("7:300999-301000", "ACGTA"), ("999-1000", "899-900")
    across = ("7:400998-401000", "ACGTAC"), ("999-?", "899-?")
    both_write = ("7:500998-501000", "ACGTAC"), ("298-300", "198-200")
    named = ("7:600999-601000", "CGTAC"), ("999-1000", "899-900")
    spelt = ("7:700999-701000", "CGTAC"), ("999-1000", "899-900")
    spelt_alleles = "7_700998_ACG/ACGTAC"
    escaping = "NMD=NMD_escaping_variant"
    reference = output(
        tmp_path / "reference.txt",
        coding_row(*taken[0], in_frame[0], *taken[1],
                   f"IMPACT=HIGH;{escaping};LoF=HC;LoF_filter=END_TRUNC;LoF_flags=SINGLE_EXON;LoF_info=PERCENTILE:0.97"),
        coding_row(*taken_shifted[0], shifted[0], *taken_shifted[1], f"IMPACT=HIGH;{escaping};LoF=HC"),
        coding_row(*two_neighbours[0], shifted[0], *two_neighbours[1], f"IMPACT=HIGH;{escaping};LoF=HC"),
        coding_row(*across[0], in_frame[0], *across[1], f"IMPACT=HIGH;{escaping};LoF=HC"),
        coding_row(*both_write[0], in_frame[0], *both_write[1], "IMPACT=HIGH;LoF=HC;LoF_info=PERCENTILE:0.9"),
        coding_row(*named[0], in_frame[0], *named[1], f"IMPACT=HIGH;{escaping};LoF=HC", uploaded="rs90"),
        coding_row(*spelt[0], in_frame[0], *spelt[1], f"IMPACT=HIGH;{escaping};LoF=HC", uploaded=spelt_alleles),
    )
    engine = output(
        tmp_path / "engine.txt",
        coding_row(*taken[0], in_frame[1], *taken[1], "IMPACT=MODERATE"),
        coding_row(*taken_shifted[0], shifted[1], *taken_shifted[1],
                   f"IMPACT=HIGH;DownstreamProtein=MKV;ProteinLengthChange=-12;{escaping};LoF=HC"),
        coding_row(*two_neighbours[0], shifted[1], *two_neighbours[1],
                   f"IMPACT=HIGH;DownstreamProtein=MKL;ProteinLengthChange=-3;{escaping};LoF=HC"),
        coding_row(*across[0], in_frame[1], *across[1], "IMPACT=MODERATE"),
        coding_row(*both_write[0], in_frame[1], *both_write[1], "IMPACT=MODERATE;LoF_info=PERCENTILE:0.5"),
        coding_row(*named[0], in_frame[1], *named[1], "IMPACT=MODERATE", uploaded="rs90"),
        coding_row(*spelt[0], in_frame[1], *spelt[1], "IMPACT=MODERATE", uploaded=spelt_alleles),
    )
    run_input = vcf(tmp_path / "input.vcf", ("7", "600998", "rs90", "ACG", "ACGTAC"), ("7", "700998", ".", "ACG", "ACGTAC"))
    return reference, engine, run_input


def test_plugin_value_that_follows_a_splice_term_past_the_record_span_is_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch):
    """``plugin_value_follows_splice_term_past_the_record_span`` is declared like any class of
    ``FIELD_CLASSES`` and decided by the comparator's own predicate, unplanted: on a row of
    ``splice_term_past_the_record_span`` an NMD or LOFTEE value the reference alone writes
    beside its donor term and a Downstream value the engine alone writes beside
    frameshift_variant and no splice term are set aside. The class reads a row with its
    record where the row's first column spells the alleles and, under an identifier, where
    the lines of ``--input`` that carry the identifier give the row by position; scored
    without the input, the row under the identifier is read without its record and its
    values stay residuals. The same values on a row of two neighbouring positions and on a
    row whose span crosses an exon's end, under an identifier no input line carries, and a
    ``LoF_info`` both sides write with two values stay residuals either way."""
    klass = formats.PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN
    assert klass == "plugin_value_follows_splice_term_past_the_record_span"
    assert klass in formats.RECORD_REFERENCE_CLASSES
    plant_registry(monkeypatch, SPLICE_SPAN_PLUGINS)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    reference, engine, run_input = splice_span_outputs(tmp_path)
    for name, given, scores in (("given", run_input, SPLICE_SPAN_SCORES), ("without", None, SPLICE_SPAN_SCORES_WITHOUT_INPUT)):
        residuals = tmp_path / f"residuals-{name}"
        report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=SPLICE_SPAN_TABLE,
                       plugins=list(SPLICE_SPAN_PLUGINS), run_input=given)
        for source, (counted, keys, left) in scores.items():
            item = entry(report, source)
            assert counts(item) == counted, (name, source)
            assert item["masked"] == {klass: keys}, (name, source)
            assert residual_lines(residuals, source) == left, (name, source)


def vcf(path: Path, *records: tuple[str, str, str, str, str]) -> Path:
    path.write_text(VCF_HEADER + "".join("\t".join((*record, ".", ".", ".")) + "\n" for record in records), encoding="utf-8")
    return path


LOST_START_PLUGINS = {"UTRAnnotator": ["5UTR_consequence", "5UTR_annotation", "Existing_uORFs"]}
LOST_START_TERM = "5_prime_UTR_premature_start_codon_loss_variant"
LOST_START_FRAMESHIFT = "5_prime_UTR_uORF_frameshift_variant"
LOST_START_GAIN = "5_prime_UTR_premature_start_codon_gain_variant"


def lost_start_table() -> str:
    """Three of the plugin's classes that read the two fields, in the committed table's order."""
    fields = "UTRAnnotator\t5UTR_consequence,5UTR_annotation\tmet\t-\n"
    return "".join(f"{klass}\t{fields}" for klass in (
        formats.UTRANNOTATOR_GAINED_START_WINDOW, formats.UTRANNOTATOR_CONSEQUENCE_PRECEDENCE,
        formats.UTRANNOTATOR_LOST_START_AFTER_INSERTION,
    ))


def utr_fields(consequence: str = "-", annotation: str = "-") -> str:
    return f"IMPACT=MODIFIER;5UTR_annotation={annotation};5UTR_consequence={consequence};Existing_uORFs=1"


def lost_start_outputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    """A made-up reference and engine output of six rows, each an insertion's on a transcript
    of its own, the input VCF of the run and Ensembl VEP's outputs under two hash seeds. On
    five the reference reports lost the start codon at the row's lower cDNA position, the
    codon directly after the insertion: the engine reports nothing (``whole``), a frameshift
    the seeded runs write on the row too (``swept``), a frameshift on a row no seeded run
    holds (``unswept``), and a start gained inside its allele under an identifier one input
    line carries (``gaining``) and under one that two carry (``shared``). On the sixth
    (``parted``) the reference reports the codon one base before that position, the
    insertion standing inside it on a row that states no strand, and the engine reports
    nothing."""
    lost = "CapDistanceToStart={}:type=uORF"
    shift = "alt_type=uORF:ref_type=uORF"
    gain = "CapDistanceToStart={}:type=uORF"
    rows = (
        ("whole", "1:1000-1001", "CCA", "ENST1", "146-147", (LOST_START_TERM, lost.format(146)), ("-", "-")),
        ("swept", "1:2000-2001", "C", "ENST2", "40-41", (LOST_START_TERM, lost.format(40)), (LOST_START_FRAMESHIFT, shift)),
        ("unswept", "1:3000-3001", "C", "ENST3", "50-51", (LOST_START_TERM, lost.format(50)), (LOST_START_FRAMESHIFT, shift)),
        ("gaining", "1:4000-4001", "CCATGC", "ENST4", "300-301", (LOST_START_TERM, lost.format(300)), (LOST_START_GAIN, gain.format(302))),
        ("shared", "1:5000-5001", "CCATGC", "ENST5", "300-301", (LOST_START_TERM, lost.format(300)), (LOST_START_GAIN, gain.format(302))),
        ("parted", "1:6000-6001", "CCA", "ENST6", "146-147", (LOST_START_TERM, lost.format(145)), ("-", "-")),
    )
    reference = output(tmp_path / "reference.txt", *(row(loc, allele, feature, utr_fields(*ref), cdna=cdna, uploaded=name)
                                                     for name, loc, allele, feature, cdna, ref, _ in rows))
    engine = output(tmp_path / "engine.txt", *(row(loc, allele, feature, utr_fields(*eng), cdna=cdna, uploaded=name)
                                               for name, loc, allele, feature, cdna, _, eng in rows))
    run_input = vcf(
        tmp_path / "input.vcf",
        ("1", "1000", "whole", "G", "GCCA"), ("1", "2000", "swept", "G", "GC"), ("1", "3000", "unswept", "G", "GC"),
        ("1", "4000", "gaining", "G", "GCCATGC"), ("1", "5000", "shared", "G", "GCCATGC"), ("1", "5000", "shared", "G", "GT"),
        ("1", "6000", "parted", "G", "GCCA"),
    )
    swept = rows[1]
    seeds = []
    for seed, fields in (("a", swept[5]), ("b", swept[6])):
        name, loc, allele, feature, cdna, _, _ = swept
        seeds.append(output(tmp_path / f"seed-{seed}.txt", row(loc, allele, feature, utr_fields(*fields), cdna=cdna, uploaded=name)))
    return reference, engine, run_input, *seeds


def test_lost_start_after_an_insertion_is_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch):
    """``utrannotator_lost_start_after_insertion`` is declared like any class of
    ``FIELD_CLASSES`` and decided by the comparator's own predicate, unplanted, after the
    plugin's two classes that stand before it in the table. With the run's input and the
    seeded runs, the two rows whose gained start the window class explains from the record's
    alleles, one of them under an identifier two input lines carry, and the row whose two
    terms the seeded runs write are credited to those classes, and the two rows neither
    explains to this one: the engine without a term and the frameshift no seeded run holds.
    Without the input and the seeded runs all five are this class's. The row whose
    insertion stands inside the codon, which states no strand to read its allele on, stays
    a residual either way."""
    klass = formats.UTRANNOTATOR_LOST_START_AFTER_INSERTION
    assert klass == "utrannotator_lost_start_after_insertion"
    assert klass not in formats.RECORD_READING_CLASSES
    plant_registry(monkeypatch, LOST_START_PLUGINS)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    reference, engine, run_input, seed_a, seed_b = lost_start_outputs(tmp_path)
    window, precedence = formats.UTRANNOTATOR_GAINED_START_WINDOW, formats.UTRANNOTATOR_CONSEQUENCE_PRECEDENCE
    parted = [
        "reference\t1:6000-6001\tCCA\tENST6\t5UTR_annotation\tCapDistanceToStart=145:type=uORF",
        f"reference\t1:6000-6001\tCCA\tENST6\t5UTR_consequence\t{LOST_START_TERM}",
    ]
    for name, given, sweep, masked in (
        ("given", run_input, [seed_a, seed_b], {window: 4, precedence: 2, klass: 4}),
        ("without", None, None, {window: 0, precedence: 0, klass: 10}),
    ):
        residuals = tmp_path / f"residuals-{name}"
        report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=lost_start_table(),
                       plugins=list(LOST_START_PLUGINS), run_input=given, sweep=sweep)
        item = entry(report, "UTRAnnotator")
        assert item["masked"] == masked, name
        assert counts(item) == (18, 14, 6, 8, 6, 6), name
        assert residual_lines(residuals, "UTRAnnotator") == parted, name


def moved_start_outputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    """A made-up reference and engine output of seven rows, each an insertion's on a
    transcript of its own, the input VCF of the run and Ensembl VEP's outputs under two hash
    seeds. On every row the reference reports lost a start codon one or two bases before the
    row's lower cDNA position, the insertion standing inside it. On five the allele, read on
    the transcript's strand, ends with the codon's first base or its first two, so the codon
    stands whole further on: the engine reports nothing (``moved``, and ``moved_two`` at two
    bases), a frameshift the seeded runs write on the row too, on a reverse-strand
    transcript (``swept``), a start gained inside its allele, which holds a start codon
    (``making``), and a start gained at the lost codon's index plus the allele's length, the
    moved codon itself (``named``). On ``parted`` the allele ends with neither, and
    ``strandless`` is ``moved`` on a row that states no strand."""
    lost = "CapDistanceToStart={}:type=uORF"
    shift = "alt_type=uORF:ref_type=uORF"
    rows = (
        ("moved", "1:1000-1001", "CCA", "1", "ENST1", "146-147", (LOST_START_TERM, lost.format(145)), ("-", "-")),
        ("swept", "1:2000-2001", "TGG", "-1", "ENST2", "146-147", (LOST_START_TERM, lost.format(145)), (LOST_START_FRAMESHIFT, shift)),
        ("moved_two", "1:3000-3001", "CAT", "1", "ENST3", "40-41", (LOST_START_TERM, lost.format(38)), ("-", "-")),
        ("making", "1:4000-4001", "CATGA", "1", "ENST4", "300-301", (LOST_START_TERM, lost.format(299)), (LOST_START_GAIN, lost.format(301))),
        ("named", "1:5000-5001", "CCA", "1", "ENST5", "300-301", (LOST_START_TERM, lost.format(299)), (LOST_START_GAIN, lost.format(302))),
        ("parted", "1:6000-6001", "CCC", "1", "ENST6", "146-147", (LOST_START_TERM, lost.format(145)), ("-", "-")),
        ("strandless", "1:7000-7001", "CCA", None, "ENST7", "146-147", (LOST_START_TERM, lost.format(145)), ("-", "-")),
    )

    def fields(strand: str | None, values: tuple[str, str]) -> str:
        return utr_fields(*values).replace("IMPACT=MODIFIER;", f"IMPACT=MODIFIER;STRAND={strand};" if strand else "IMPACT=MODIFIER;", 1)

    reference = output(tmp_path / "reference.txt", *(row(loc, allele, feature, fields(strand, ref), cdna=cdna, uploaded=name)
                                                     for name, loc, allele, strand, feature, cdna, ref, _ in rows))
    engine = output(tmp_path / "engine.txt", *(row(loc, allele, feature, fields(strand, eng), cdna=cdna, uploaded=name)
                                               for name, loc, allele, strand, feature, cdna, _, eng in rows))
    run_input = vcf(
        tmp_path / "input.vcf",
        ("1", "1000", "moved", "G", "GCCA"), ("1", "2000", "swept", "G", "GTGG"), ("1", "3000", "moved_two", "G", "GCAT"),
        ("1", "4000", "making", "G", "GCATGA"), ("1", "5000", "named", "G", "GCCA"), ("1", "6000", "parted", "G", "GCCC"),
        ("1", "7000", "strandless", "G", "GCCA"),
    )
    name, loc, allele, strand, feature, cdna, ref, eng = rows[1]
    seeds = [output(tmp_path / f"seed-{seed}.txt", row(loc, allele, feature, fields(strand, values), cdna=cdna, uploaded=name))
             for seed, values in (("a", ref), ("b", eng))]
    return reference, engine, run_input, *seeds


def test_lost_start_at_a_codon_an_insertion_only_moves_on_is_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch):
    """``utrannotator_lost_start_after_insertion`` takes a row on which the reference reports
    lost a start codon the insertion stands inside and only moves on, decided by the
    comparator's own predicate, unplanted, after the plugin's two classes that stand before
    it in the table. With the run's input and the seeded runs, the gained start of a codon
    the inserted bases make is the window class's and the frameshift the seeded runs write
    the precedence class's; the engine without a term, at one base and at two, and the
    gained start that names the moved codon, which the window class refuses, are this
    class's. Without the input and the seeded runs all five are. The row whose allele does
    not end with the codon's first bases and the row that states no strand stay residuals
    either way."""
    klass = formats.UTRANNOTATOR_LOST_START_AFTER_INSERTION
    plant_registry(monkeypatch, LOST_START_PLUGINS)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    reference, engine, run_input, seed_a, seed_b = moved_start_outputs(tmp_path)
    window, precedence = formats.UTRANNOTATOR_GAINED_START_WINDOW, formats.UTRANNOTATOR_CONSEQUENCE_PRECEDENCE
    left = [
        "reference\t1:6000-6001\tCCC\tENST6\t5UTR_annotation\tCapDistanceToStart=145:type=uORF",
        f"reference\t1:6000-6001\tCCC\tENST6\t5UTR_consequence\t{LOST_START_TERM}",
        "reference\t1:7000-7001\tCCA\tENST7\t5UTR_annotation\tCapDistanceToStart=145:type=uORF",
        f"reference\t1:7000-7001\tCCA\tENST7\t5UTR_consequence\t{LOST_START_TERM}",
    ]
    for name, given, sweep, masked in (
        ("given", run_input, [seed_a, seed_b], {window: 2, precedence: 2, klass: 6}),
        ("without", None, None, {window: 0, precedence: 0, klass: 10}),
    ):
        residuals = tmp_path / f"residuals-{name}"
        report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=lost_start_table(),
                       plugins=list(LOST_START_PLUGINS), run_input=given, sweep=sweep)
        item = entry(report, "UTRAnnotator")
        assert item["masked"] == masked, name
        assert counts(item) == (21, 13, 7, 11, 7, 7), name
        assert residual_lines(residuals, "UTRAnnotator") == left, name


def remade_start_table() -> str:
    """Four of the plugin's classes that read the two fields, in the committed table's order."""
    fields = "UTRAnnotator\t5UTR_consequence,5UTR_annotation\tmet\t-\n"
    return lost_start_table() + f"{formats.UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE}\t{fields}"


def remade_start_outputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    """A made-up reference and engine output of seven rows, each an insertion's on a
    transcript of its own, the input VCF of the run and Ensembl VEP's outputs under two hash
    seeds. On five the inserted bases, read on the transcript's strand, open with the last
    bases of a start codon one or two bases before the row's lower cDNA position and leave
    it at its own index. The reference reports that codon gained on ``kept`` (the engine
    reports nothing), on ``kept_rev``, a reverse-strand transcript (the engine a frameshift
    the seeded runs write on the row too), on ``shifted`` at two bases (the engine a
    frameshift no seeded run holds) and on ``making`` (the engine a start gained further
    on, inside its allele). On ``swept`` both report the frameshift, the engine's listing
    one ORF more, and a seeded run reports the gained start. On ``parted`` the allele opens
    with neither, and ``strandless`` is ``kept`` on a row that states no strand."""
    gain = "CapDistanceToStart={}:type=uORF"
    shift = "alt_type=uORF:ref_type=uORF"
    own = f"{shift},alt_type=inFrame_oORF:ref_type=uORF"
    rows = (
        ("kept", "1:1000-1001", "TGC", "1", "ENST1", "146-147", (LOST_START_GAIN, gain.format(145)), ("-", "-")),
        ("kept_rev", "1:2000-2001", "GCA", "-1", "ENST2", "146-147", (LOST_START_GAIN, gain.format(145)), (LOST_START_FRAMESHIFT, shift)),
        ("shifted", "1:3000-3001", "G", "1", "ENST3", "40-41", (LOST_START_GAIN, gain.format(38)), (LOST_START_FRAMESHIFT, own)),
        ("swept", "1:4000-4001", "TG", "1", "ENST4", "300-301", (LOST_START_FRAMESHIFT, shift), (LOST_START_FRAMESHIFT, own)),
        ("making", "1:5000-5001", "TGATGC", "1", "ENST5", "300-301", (LOST_START_GAIN, gain.format(299)), (LOST_START_GAIN, gain.format(302))),
        ("parted", "1:6000-6001", "CCC", "1", "ENST6", "146-147", (LOST_START_GAIN, gain.format(145)), ("-", "-")),
        ("strandless", "1:7000-7001", "TGC", None, "ENST7", "146-147", (LOST_START_GAIN, gain.format(145)), ("-", "-")),
    )

    def fields(strand: str | None, values: tuple[str, str]) -> str:
        return utr_fields(*values).replace("IMPACT=MODIFIER;", f"IMPACT=MODIFIER;STRAND={strand};" if strand else "IMPACT=MODIFIER;", 1)

    def line(index: int, values: tuple[str, str]) -> str:
        name, loc, allele, strand, feature, cdna, _, _ = rows[index]
        return row(loc, allele, feature, fields(strand, values), cdna=cdna, uploaded=name)

    reference = output(tmp_path / "reference.txt", *(line(i, entry[6]) for i, entry in enumerate(rows)))
    engine = output(tmp_path / "engine.txt", *(line(i, entry[7]) for i, entry in enumerate(rows)))
    run_input = vcf(
        tmp_path / "input.vcf",
        ("1", "1000", "kept", "G", "GTGC"), ("1", "2000", "kept_rev", "G", "GGCA"), ("1", "3000", "shifted", "G", "GG"),
        ("1", "4000", "swept", "G", "GTG"), ("1", "5000", "making", "G", "GTGATGC"), ("1", "6000", "parted", "G", "GCCC"),
        ("1", "7000", "strandless", "G", "GTGC"),
    )
    seed_a = output(tmp_path / "seed-a.txt", line(1, rows[1][6]), line(3, rows[3][6]))
    seed_b = output(tmp_path / "seed-b.txt", line(1, rows[1][7]), line(3, (LOST_START_GAIN, gain.format(299))))
    return reference, engine, run_input, seed_a, seed_b


def test_gained_start_of_a_codon_left_at_its_own_index_is_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch):
    """``utrannotator_gained_start_remade_in_place`` is declared like any class of
    ``FIELD_CLASSES`` and decided by the comparator's own predicate, unplanted, after the
    plugin's three classes that stand before it in the table. With the run's input and the
    seeded runs: the frameshift the seeded runs write beside the gained start is the
    precedence class's; the engine without a term, the frameshift no seeded run holds, the
    annotation of the start the engine gains further on, and the annotation of the row on
    which both report the frameshift and a seeded run the gained start are this class's.
    Without the seeded runs the first frameshift is this class's too and the last annotation
    stays a residual on both sides. The row whose allele does not open with the codon's last
    bases and the row that states no strand stay residuals either way."""
    klass = formats.UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE
    assert klass == "utrannotator_gained_start_remade_in_place"
    assert klass not in formats.RECORD_READING_CLASSES
    plant_registry(monkeypatch, LOST_START_PLUGINS)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    reference, engine, run_input, seed_a, seed_b = remade_start_outputs(tmp_path)
    window, precedence, lost = (formats.UTRANNOTATOR_GAINED_START_WINDOW, formats.UTRANNOTATOR_CONSEQUENCE_PRECEDENCE,
                                formats.UTRANNOTATOR_LOST_START_AFTER_INSERTION)
    left = [
        "reference\t1:6000-6001\tCCC\tENST6\t5UTR_annotation\tCapDistanceToStart=145:type=uORF",
        f"reference\t1:6000-6001\tCCC\tENST6\t5UTR_consequence\t{LOST_START_GAIN}",
        "reference\t1:7000-7001\tTGC\tENST7\t5UTR_annotation\tCapDistanceToStart=145:type=uORF",
        f"reference\t1:7000-7001\tTGC\tENST7\t5UTR_consequence\t{LOST_START_GAIN}",
    ]
    swept = [
        "reference\t1:4000-4001\tTG\tENST4\t5UTR_annotation\talt_type=uORF:ref_type=uORF",
        "engine\t1:4000-4001\tTG\tENST4\t5UTR_annotation\talt_type=uORF:ref_type=uORF,alt_type=inFrame_oORF:ref_type=uORF",
    ]
    for name, given, sweep, masked, counted, residual in (
        ("given", run_input, [seed_a, seed_b], {window: 0, precedence: 2, lost: 0, klass: 6}, (21, 15, 9, 13, 9, 9), left),
        ("without", None, None, {window: 0, precedence: 0, lost: 0, klass: 7}, (21, 15, 9, 14, 10, 9), [*swept, *left]),
    ):
        residuals = tmp_path / f"residuals-{name}"
        report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=remade_start_table(),
                       plugins=list(LOST_START_PLUGINS), run_input=given, sweep=sweep)
        item = entry(report, "UTRAnnotator")
        assert item["masked"] == masked, name
        assert counts(item) == counted, name
        assert sorted(residual_lines(residuals, "UTRAnnotator")) == sorted(residual), name


def stop_boundary_table() -> str:
    """Five of the plugin's classes that read its two fields or the annotation alone, in the
    committed table's order."""
    return remade_start_table() + f"{formats.UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY}\tUTRAnnotator\t5UTR_annotation\tmet\t-\n"


STOP_BOUNDARY_THROUGH = "alt_type=uORF:alt_type_length=150:ref_StartDistanceToCDS=120:ref_type=OutOfFrame_oORF:ref_type_length=NA"
STOP_BOUNDARY_ENDS = "alt_type=uORF:alt_type_length=45:ref_StartDistanceToCDS=87:ref_type=uORF:ref_type_length=36"


def stop_boundary_outputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    """A made-up reference and engine output of five rows, each an insertion's on a
    transcript of its own, the input VCF of the run and Ensembl VEP's outputs under two hash
    seeds. On every row both report the frameshift and the engine lists one uORF more, whose
    stop codon begins 54 bases before the coding start. The row's HGVSc places the base
    directly after the insertion there on ``boundary``, on ``shifted_back`` once its
    ``HGVS_OFFSET`` of two bases is taken back, and on ``swept``, where a seeded run reports
    a start codon gained at the index the inserted base leaves it at; one base nearer the
    coding start on ``elsewhere``; and ``unstated`` carries no HGVSc."""
    rows = (
        ("boundary", "1:1000-1001", "A", "ENST1", "HGVSc=ENST1.1:c.-55_-54insA"),
        ("shifted_back", "1:2000-2001", "A", "ENST2", "HGVSc=ENST2.1:c.-53_-52insA;HGVS_OFFSET=2"),
        ("swept", "1:3000-3001", "G", "ENST3", "HGVSc=ENST3.1:c.-55_-54insG"),
        ("elsewhere", "1:4000-4001", "A", "ENST4", "HGVSc=ENST4.1:c.-54_-53insA"),
        ("unstated", "1:5000-5001", "A", "ENST5", ""),
    )

    def line(index: int, annotation: str, consequence: str = LOST_START_FRAMESHIFT) -> str:
        name, loc, allele, feature, notation = rows[index]
        fields = utr_fields(consequence, annotation).replace("IMPACT=MODIFIER;", "IMPACT=MODIFIER;STRAND=1;", 1)
        return row(loc, allele, feature, fields + (f";{notation}" if notation else ""), cdna="146-147", uploaded=name)

    longer = f"{STOP_BOUNDARY_THROUGH},{STOP_BOUNDARY_ENDS}"
    reference = output(tmp_path / "reference.txt", *(line(i, STOP_BOUNDARY_THROUGH) for i in range(len(rows))))
    engine = output(tmp_path / "engine.txt", *(line(i, longer) for i in range(len(rows))))
    run_input = vcf(tmp_path / "input.vcf", *(
        ("1", loc.partition(":")[2].partition("-")[0], name, "G", f"G{allele}") for name, loc, allele, _, _ in rows
    ))
    seed_a = output(tmp_path / "seed-a.txt", line(2, STOP_BOUNDARY_THROUGH))
    seed_b = output(tmp_path / "seed-b.txt", line(2, "CapDistanceToStart=144:type=uORF", LOST_START_GAIN))
    return reference, engine, run_input, seed_a, seed_b


def test_frameshift_list_at_a_stop_codon_boundary_is_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch):
    """``utrannotator_frameshift_at_stop_boundary`` is declared like any class of
    ``FIELD_CLASSES`` and decided by the comparator's own predicate, unplanted, after the
    plugin's four classes that stand before it in the table; it reads no record. With the
    seeded runs, the annotation of the row on which one of them reports the start codon the
    inserted base leaves in place is the remade-in-place class's; without them it is this
    class's, as the two rows whose HGVSc places the insertion at the added uORF's stop codon
    are either way. The row whose HGVSc places it elsewhere and the row without an HGVSc
    stay residuals on both sides."""
    klass = formats.UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY
    assert klass == "utrannotator_frameshift_at_stop_boundary"
    assert klass not in formats.RECORD_READING_CLASSES
    assert COMPARATOR_FIELD_CLASSES[klass] == formats.UTRANNOTATOR_ANNOTATION_FIELDS
    plant_registry(monkeypatch, LOST_START_PLUGINS)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    reference, engine, run_input, seed_a, seed_b = stop_boundary_outputs(tmp_path)
    window, precedence, lost, remade = (
        formats.UTRANNOTATOR_GAINED_START_WINDOW, formats.UTRANNOTATOR_CONSEQUENCE_PRECEDENCE,
        formats.UTRANNOTATOR_LOST_START_AFTER_INSERTION, formats.UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE,
    )
    left = [
        line
        for loc, feature in (("1:4000-4001", "ENST4"), ("1:5000-5001", "ENST5"))
        for line in (
            f"reference\t{loc}\tA\t{feature}\t5UTR_annotation\t{STOP_BOUNDARY_THROUGH}",
            f"engine\t{loc}\tA\t{feature}\t5UTR_annotation\t{STOP_BOUNDARY_THROUGH},{STOP_BOUNDARY_ENDS}",
        )
    ]
    reports = {}
    for name, given, sweep, masked in (
        ("given", run_input, [seed_a, seed_b], {window: 0, precedence: 0, lost: 0, remade: 1, klass: 2}),
        ("without", None, None, {window: 0, precedence: 0, lost: 0, remade: 0, klass: 3}),
    ):
        residuals = tmp_path / f"residuals-{name}"
        report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=stop_boundary_table(),
                       plugins=list(LOST_START_PLUGINS), run_input=given, sweep=sweep)
        item = entry(report, "UTRAnnotator")
        assert item["masked"] == masked, name
        assert sorted(residual_lines(residuals, "UTRAnnotator")) == sorted(left), name
        reports[name] = counts(item)
    # The same keys either way: of the fifteen each side writes, ten are equal, three are set
    # aside on both sides and two are left a side.
    assert reports["given"] == reports["without"] == (15, 15, 10, 12, 12, 10)


def before_stop_table() -> str:
    """Six of the plugin's classes that read its two fields or the annotation alone, in the
    committed table's order."""
    fields = "UTRAnnotator\t5UTR_consequence,5UTR_annotation\tmet\t-\n"
    return stop_boundary_table() + f"{formats.UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON}\t{fields}"


BEFORE_STOP_SHIFTED = f"{STOP_BOUNDARY_THROUGH},KozakContext=GCCATGG:KozakStrength=Strong:{STOP_BOUNDARY_ENDS}"
BEFORE_STOP_LOST = "AltStop=True:AltStopDistanceToCDS={}:FrameWithCDS=outOfFrame:KozakContext=GCCATGG:KozakStrength=Strong"
LOST_START_STOP_LOSS = "5_prime_UTR_uORF_stop_codon_loss_variant"


def before_stop_outputs(tmp_path: Path) -> tuple[Path, Path, Path, Path, Path]:
    """A made-up reference and engine output of four rows, each of two bases on a transcript
    of its own between cDNA positions 146 and 147, the input VCF of the run and Ensembl VEP's
    outputs under two hash seeds. On every row the reference reports the frameshift and
    lists a uORF that starts 87 bases before the coding start and is 45 bases long once
    shifted, and the engine reports the lost stop. The engine's group names that uORF's
    next stop codon, 47 bases before the coding start, on ``before``, on ``swept``, where
    the seeded runs write both terms on the row, and on ``replaced``, whose record replaces
    one base; it names a stop codon one base off on ``elsewhere``."""
    rows = (
        ("before", "1:1000-1001", "ENST1", 47, ("G", "GAT")),
        ("swept", "1:2000-2001", "ENST2", 47, ("G", "GAT")),
        ("elsewhere", "1:3000-3001", "ENST3", 48, ("G", "GAT")),
        ("replaced", "1:4000-4001", "ENST4", 47, ("GC", "GAT")),
    )

    def line(index: int, engine: bool) -> str:
        name, loc, feature, distance, _ = rows[index]
        values = (LOST_START_STOP_LOSS, BEFORE_STOP_LOST.format(distance)) if engine else (LOST_START_FRAMESHIFT, BEFORE_STOP_SHIFTED)
        fields = utr_fields(*values).replace("IMPACT=MODIFIER;", "IMPACT=MODIFIER;STRAND=1;", 1)
        return row(loc, "AT", feature, fields, cdna="146-147", uploaded=name)

    reference = output(tmp_path / "reference.txt", *(line(i, False) for i in range(len(rows))))
    engine = output(tmp_path / "engine.txt", *(line(i, True) for i in range(len(rows))))
    run_input = vcf(tmp_path / "input.vcf", *(
        ("1", loc.partition(":")[2].partition("-")[0], name, *alleles) for name, loc, _, _, alleles in rows
    ))
    seed_a = output(tmp_path / "seed-a.txt", line(1, False))
    seed_b = output(tmp_path / "seed-b.txt", line(1, True))
    return reference, engine, run_input, seed_a, seed_b


def test_lost_stop_written_before_the_stop_codon_is_set_aside_by_the_comparators_own_rule(tmp_path, monkeypatch):
    """``utrannotator_lost_stop_written_before_stop_codon`` is declared like any class of
    ``FIELD_CLASSES`` and decided by the comparator's own predicate, unplanted, after the
    plugin's five classes that stand before it in the table. With the run's input and the
    seeded runs: the row on which those runs write both terms is the precedence class's,
    the row whose record states a reference allele that is not empty stays a residual, and
    the row whose group names the listed uORF's next stop codon is this class's with both
    fields. Without them the three are this class's. The row whose group names another
    stop codon stays a residual on both sides either way."""
    klass = formats.UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON
    assert klass == "utrannotator_lost_stop_written_before_stop_codon"
    assert klass not in formats.RECORD_READING_CLASSES
    assert COMPARATOR_FIELD_CLASSES[klass] == formats.UTRANNOTATOR_CONSEQUENCE_FIELDS + formats.UTRANNOTATOR_ANNOTATION_FIELDS
    plant_registry(monkeypatch, LOST_START_PLUGINS)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    reference, engine, run_input, seed_a, seed_b = before_stop_outputs(tmp_path)
    window, precedence, lost, remade, boundary = (
        formats.UTRANNOTATOR_GAINED_START_WINDOW, formats.UTRANNOTATOR_CONSEQUENCE_PRECEDENCE,
        formats.UTRANNOTATOR_LOST_START_AFTER_INSERTION, formats.UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE,
        formats.UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY,
    )

    def left(loc: str, feature: str, distance: int) -> list[str]:
        return [
            f"reference\t{loc}\tAT\t{feature}\t5UTR_annotation\t{BEFORE_STOP_SHIFTED}",
            f"reference\t{loc}\tAT\t{feature}\t5UTR_consequence\t{LOST_START_FRAMESHIFT}",
            f"engine\t{loc}\tAT\t{feature}\t5UTR_annotation\t{BEFORE_STOP_LOST.format(distance)}",
            f"engine\t{loc}\tAT\t{feature}\t5UTR_consequence\t{LOST_START_STOP_LOSS}",
        ]

    elsewhere, replaced = left("1:3000-3001", "ENST3", 48), left("1:4000-4001", "ENST4", 47)
    for name, given, sweep, masked, counted, residual in (
        ("given", run_input, [seed_a, seed_b], {window: 0, precedence: 2, lost: 0, remade: 0, boundary: 0, klass: 2},
         (12, 12, 4, 8, 8, 4), [*elsewhere, *replaced]),
        ("without", None, None, {window: 0, precedence: 0, lost: 0, remade: 0, boundary: 0, klass: 6}, (12, 12, 4, 6, 6, 4), elsewhere),
    ):
        residuals = tmp_path / f"residuals-{name}"
        report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=before_stop_table(),
                       plugins=list(LOST_START_PLUGINS), run_input=given, sweep=sweep)
        item = entry(report, "UTRAnnotator")
        assert item["masked"] == masked, name
        assert counts(item) == counted, name
        assert sorted(residual_lines(residuals, "UTRAnnotator")) == sorted(residual), name


def plant_record_reader(monkeypatch) -> list:
    """Plants a reader of a record's rows with the interface of ``RecordRows`` and returns,
    per reader built, the (Location, Allele) of the entries and the input lines it was given."""
    built = []

    class Records:
        def __init__(self, entries, input_lines=None):
            built.append(([(e.location, e.allele) for e in entries], list(input_lines or ())))
            self.entries = entries
            columns = [line.split("\t") for line in input_lines or ()]
            names = [line[2] for line in columns]
            self.lines = {line[2]: line for line in columns if names.count(line[2]) == 1}

        def of(self, e):
            name = e.fields["Uploaded_variation"]
            if name not in self.lines:
                return None
            rows = [r for r in self.entries if (r.location, r.fields["Uploaded_variation"]) == (e.location, name)]
            fields = dict(zip(("CHROM", "POS", "Uploaded_variation", "REF", "ALT"), self.lines[name]))
            return [*rows, formats.Entry(f"{e.location}|record", None, e.location, "", "", "", fields)]

    monkeypatch.setattr(formats, "RecordRows", Records, raising=False)
    return built


def named_rows(tmp_path: Path) -> tuple[Path, Path, Path]:
    other = row("1:600", "A", "ENST2", "REVEL=0.7", uploaded="rs66")
    reference = output(
        tmp_path / "reference.txt",
        *(row(location, allele, "ENST1", "REVEL=0.7", uploaded=name) for location, name in NAMED for allele in "AG"),
        other,
    )
    engine = output(
        tmp_path / "engine.txt",
        *(
            row(location, allele, "ENST1", "REVEL=0.7" if allele == "A" else "IMPACT=LOW", uploaded=name)
            for location, name in NAMED
            for allele in "AG"
        ),
        other,
    )
    paired = [*((location, "A", "ENST1", "REVEL", "0.7") for location, _ in NAMED), ("1:600", "A", "ENST2", "REVEL", "0.7")]
    return reference, engine, tuples(tmp_path / "engine.tsv", *paired)


NAMED = (("1:600", "rs6"), ("1:800", "rs8"), ("1:900", "rs9"))
NAMED_RECORDS = (
    ("1", "600", "rs6", "C", "A,G"),
    ("1", "600", "rs66", "C", "A"),
    ("1", "800", "rs8", "C", "A"),
    ("1", "800", "rs8", "C", "G"),
    ("1", "50", ".", "C", "T"),
)


@pytest.mark.parametrize("zipped", [False, True], ids=["plain", "gzip"])
def test_input_names_the_record_of_rows_under_an_identifier_that_one_line_carries(tmp_path, monkeypatch, planted, zipped):
    built = plant_record_reader(monkeypatch)
    reference, engine, listed = named_rows(tmp_path)
    run_input = vcf(tmp_path / "input.vcf", *NAMED_RECORDS)
    if zipped:
        run_input = tmp_path / "input.vcf.gz"
        run_input.write_bytes(gzip.compress((tmp_path / "input.vcf").read_bytes()))
    engines = {"engine": engine, "listed": listed}
    left = [f"reference\t{location}\tG\tENST1\tREVEL\t0.7" for location, _ in NAMED]
    without = score(tmp_path, reference, engines, residuals=tmp_path / "without")
    assert without["run_input"] is None and built == []
    for ident in engines:
        assert counts(entry(without, "REVEL", engine=ident)) == (7, 4, 4, 7, 4, 4)
        assert residual_lines(tmp_path / "without", "REVEL", ident) == left
    assert all(call[6] is None for call in planted)
    planted.clear()
    report = score(tmp_path, reference, engines, residuals=tmp_path / "with", run_input=run_input)
    assert report["run_input"] == {
        "path": str(run_input),
        "sha256": hashlib.sha256(run_input.read_bytes()).hexdigest(),
        "records": 5,
    }
    for ident in engines:
        item = entry(report, "REVEL", engine=ident)
        assert counts(item) == (7, 4, 4, 6, 4, 4)
        assert item["masked"] == {"planted_sibling_value": 1}
        assert residual_lines(tmp_path / "with", "REVEL", ident) == left[1:]
    handed = {call[1].location: call[6] for call in planted}
    assert [(r.allele, r.fields.get("ALT")) for r in handed["1:600"]] == [("A", None), ("G", None), ("", "A,G")]
    assert handed["1:800"] is None and handed["1:900"] is None
    assert built == [
        ([("1:600", "A"), ("1:600", "G"), ("1:600", "A")], ["1\t600\trs6\tC\tA,G", "1\t600\trs66\tC\tA"]),
        ([("1:800", "A"), ("1:800", "G")], ["1\t800\trs8\tC\tA", "1\t800\trs8\tC\tG"]),
        ([("1:900", "A"), ("1:900", "G")], []),
    ] * 2


def test_allele_order_of_a_record_under_an_identifier_comes_from_its_input_line(tmp_path, monkeypatch):
    plant_record_reader(monkeypatch)

    def accepts(klass, e, name, expected_value, actual_value, a=None, record_rows=None):
        return any("ALT" in r.fields and e.allele != r.fields["ALT"].split(",")[0] for r in record_rows or ())

    monkeypatch.setattr(formats, "FIELD_CLASSES", {"planted_later_alternate": ()})
    monkeypatch.setattr(formats, "field_class_accepts", accepts)
    table = "planted_later_alternate\tREVEL\tREVEL\n"
    reference = output(tmp_path / "reference.txt", *(row("1:100", allele, "ENST1", "REVEL=0.7", uploaded="rs1") for allele in "GT"))
    engine = output(tmp_path / "engine.txt", *(row("1:100", allele, "ENST1", "REVEL=0.9", uploaded="rs1") for allele in "GT"))
    run_input = vcf(tmp_path / "input.vcf", ("1", "100", "rs1", "A", "G,T"))
    without = entry(score(tmp_path, reference, {"engine": engine}, table=table), "REVEL")
    assert (counts(without), without["masked"]) == ((2, 2, 0, 2, 2, 0), {"planted_later_alternate": 0})
    residuals = tmp_path / "residuals"
    item = entry(score(tmp_path, reference, {"engine": engine}, table=table, residuals=residuals, run_input=run_input), "REVEL")
    assert (counts(item), item["masked"]) == ((2, 2, 0, 1, 1, 0), {"planted_later_alternate": 1})
    assert residual_lines(residuals, "REVEL") == [
        "reference\t1:100\tG\tENST1\tREVEL\t0.7",
        "engine\t1:100\tG\tENST1\tREVEL\t0.9",
    ]


def test_input_index_finds_every_line_that_carries_an_identifier(tmp_path):
    """A line carries the first ``;``-joined element of its ID column: the lines of ``rs1;rs2``
    are found under ``rs1`` with the lines of ``rs1`` alone, and under no other name."""
    rng = random.Random(7)
    names = ["rs1", "rs10", "rs100", "rs1;rs2", "rs2", "a", "a b", "é", "zz", "."]
    records = [(str(rng.randint(1, 3)), str(rng.randint(1, 500)), rng.choice(names), "C", rng.choice("AGT")) for _ in range(600)]
    carried = defaultdict(list)
    for chromosome, position, name, reference, alternate in records:
        first = name.split(";")[0]
        carried[first].append(f"{chromosome}\t{position}\t{first}\t{reference}\t{alternate}")
    work = tmp_path / "work"
    work.mkdir()
    index = spt.InputIndex(vcf(tmp_path / "input.vcf", *records), spt.Sorter(work, "16M", None), work)
    try:
        assert index.records == 600
        for name in names[:-1]:
            assert sorted(index.named(name.encode())) == sorted(carried[name]), name
        for absent in (".", "rs", "rs11", "rs1;", "", "zzz", "0"):
            assert index.named(absent.encode()) == [], absent
    finally:
        index.close()
    unnamed = work / "unnamed"
    unnamed.mkdir()
    index = spt.InputIndex(vcf(tmp_path / "unnamed.vcf", ("1", "5", ".", "C", "T")), spt.Sorter(unnamed, "16M", None), unnamed)
    try:
        assert (index.records, index.named(b"rs1")) == (1, [])
    finally:
        index.close()


def test_input_index_files_a_line_under_the_first_identifier_of_its_id_column(tmp_path):
    """A record whose ID column joins several identifiers is found under the first alone, the
    name both engines print for its rows, and is handed on with that identifier as its third
    column: under no later identifier, not under the whole column and not under the first
    identifier with a bare semicolon."""
    records = [
        ("1", "100", "rs1;rs2", "C", "T"),
        ("1", "200", "rs3;rs4;rs5;rs6;rs7", "C", "G"),
        ("2", "300", "rs2", "C", "A"),
        ("2", "400", "rs9;rs1", "C", "A"),
    ]
    work = tmp_path / "work"
    work.mkdir()
    index = spt.InputIndex(vcf(tmp_path / "input.vcf", *records), spt.Sorter(work, "16M", None), work)
    try:
        assert index.records == 4
        assert index.named(b"rs1") == ["1\t100\trs1\tC\tT"]
        assert index.named(b"rs3") == ["1\t200\trs3\tC\tG"]
        assert index.named(b"rs2") == ["2\t300\trs2\tC\tA"]
        assert index.named(b"rs9") == ["2\t400\trs9\tC\tA"]
        for absent in ("rs1;rs2", "rs1;", ";rs2", "rs4", "rs7", "rs3;rs4;rs5;rs6;rs7", "rs9;rs1"):
            assert index.named(absent.encode()) == [], absent
    finally:
        index.close()


def test_input_index_files_no_line_under_a_first_element_that_names_nothing(tmp_path):
    """A first element that is a dot, empty or ``0`` is no identifier in either engine's
    output, so its line is found under no name, the elements behind it included."""
    records = [
        ("3", "500", ".;rs10", "C", "T"),
        ("3", "600", ";rs11", "C", "T"),
        ("3", "700", "0;rs12", "C", "T"),
        ("3", "800", "0", "C", "T"),
        ("3", "900", "", "C", "T"),
        ("3", "950", "rs13", "C", "T"),
    ]
    work = tmp_path / "work"
    work.mkdir()
    index = spt.InputIndex(vcf(tmp_path / "input.vcf", *records), spt.Sorter(work, "16M", None), work)
    try:
        assert (index.records, index.named(b"rs13")) == (6, ["3\t950\trs13\tC\tT"])
        for absent in (".", "", "0", "rs10", "rs11", "rs12", ".;rs10", ";rs11", "0;rs12"):
            assert index.named(absent.encode()) == [], absent
    finally:
        index.close()


def first_identifier_rows(tmp_path, monkeypatch) -> tuple[str, dict, Path, Path]:
    """The window class under the committed table's own row with the comparator's predicate,
    unplanted, and an insertion's reference and engine rows under the identifier ``rs1``."""
    klass = "utrannotator_gained_start_window"
    declared = [line for line in COMMITTED_CLASS_TABLE.read_text(encoding="utf-8").splitlines() if line.startswith(f"{klass}\t")]
    assert len(declared) == 1
    plant_registry(monkeypatch, {"UTRAnnotator": ["5UTR_annotation", "5UTR_consequence"]})
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    reference = output(tmp_path / "reference.txt", row("1:100-101", "CAT", "ENST1", UTR_START_LOST, cdna="64-65", uploaded="rs1"))
    engine = output(tmp_path / "engine.txt", row("1:100-101", "CAT", "ENST1", UTR_START_GAINED, cdna="64-65", uploaded="rs1"))
    return klass, {"table": declared[0] + "\n", "plugins": ["UTRAnnotator"]}, reference, engine


def test_row_under_the_first_identifier_of_a_several_identifier_line_is_handed_its_record(tmp_path, monkeypatch):
    """An insertion's row prints the first identifier of its record's ID column, so the input
    line that joins two identifiers, or five, is its record, and its two fields are set aside
    as under a line of that identifier alone. A line that carries the row's identifier as a
    later element, or behind a first element that is a dot, empty or ``0``, is no record of
    the row, and both fields stay residuals."""
    klass, scored, reference, engine = first_identifier_rows(tmp_path, monkeypatch)
    grid = {
        "one identifier": ("rs1", 2),
        "two identifiers": ("rs1;rs2", 2),
        "five identifiers": ("rs1;rs2;rs3;rs4;rs5", 2),
        "the identifier as a later element": ("rs2;rs1", 0),
        "a first element that is a dot": (".;rs1", 0),
        "a first element that is empty": (";rs1", 0),
        "a first element that is zero": ("0;rs1", 0),
    }
    for number, (what, (name, keys)) in enumerate(grid.items()):
        work = tmp_path / f"grid{number}"
        work.mkdir()
        run_input = vcf(work / "input.vcf", ("1", "100", name, "C", "CCAT"))
        item = entry(score(work, reference, {"engine": engine}, residuals=work / "residuals", run_input=run_input, **scored), "UTRAnnotator")
        assert (item["masked"], item["residual_rows"]) == ({klass: keys}, 4 - 2 * keys), what
    tool = built_pairing_tool()
    if tool is not None and not THROUGH_PAIRING_TOOL:
        paired = tmp_path / "paired"
        paired.mkdir()
        run_input = vcf(paired / "input.vcf", ("1", "100", "rs1;rs2;rs3;rs4;rs5", "C", "CCAT"))
        again = entry(score(paired, reference, {"engine": engine}, run_input=run_input, tool=tool, **scored), "UTRAnnotator")
        assert (again["masked"], again["residual_rows"]) == ({klass: 2}, 0)


def test_lines_that_share_a_first_identifier_are_told_apart_by_place_and_allele(tmp_path, monkeypatch):
    """A line of one identifier and a line that joins that identifier with another stand under
    one name, as two joined lines with one first identifier do. A row under the name is handed
    the one of them that writes it at its Location with its allele, whichever of the two it
    is and whether the other stands elsewhere or spells another allele; where both write the
    row nothing is handed and both fields stay residuals."""
    klass, scored, reference, engine = first_identifier_rows(tmp_path, monkeypatch)
    own, joined = ("1", "100", "rs1", "C", "CCAT"), ("1", "100", "rs1;rs2", "C", "CCAT")
    grid = {
        "the joined line writes the row, the plain line stands elsewhere": ((("1", "150", "rs1", "CA", "C"), joined), 2),
        "the joined line writes the row, the plain line spells another allele": ((("1", "100", "rs1", "C", "CCAG"), joined), 2),
        "the plain line writes the row, the joined line stands elsewhere": ((own, ("1", "150", "rs1;rs2", "CA", "C")), 2),
        "the plain line writes the row, the joined line spells another allele": ((own, ("1", "100", "rs1;rs2", "C", "CCAG")), 2),
        "both lines write the row": ((own, joined), 0),
        "two joined lines of which one writes the row": ((("1", "150", "rs1;rs3", "CA", "C"), joined), 2),
        "two joined lines that both write the row": ((("1", "100", "rs1;rs3", "C", "CCAT"), joined), 0),
    }
    for number, (what, (lines, keys)) in enumerate(grid.items()):
        work = tmp_path / f"grid{number}"
        work.mkdir()
        item = entry(score(work, reference, {"engine": engine}, residuals=work / "residuals", run_input=vcf(work / "input.vcf", *lines), **scored), "UTRAnnotator")
        assert (item["masked"], item["residual_rows"]) == ({klass: keys}, 4 - 2 * keys), what
    work = tmp_path / "index"
    work.mkdir()
    index = spt.InputIndex(vcf(tmp_path / "shared.vcf", own, joined, ("1", "150", "rs1;rs3", "CA", "C")), spt.Sorter(work, "16M", None), work)
    try:
        assert sorted(index.named(b"rs1")) == ["1\t100\trs1\tC\tCCAT", "1\t100\trs1\tC\tCCAT", "1\t150\trs1\tCA\tC"]
        assert index.named(b"rs2") == index.named(b"rs3") == []
    finally:
        index.close()


def test_identifier_that_is_not_utf8_names_the_lines_that_carry_its_bytes(tmp_path, monkeypatch):
    built = plant_record_reader(monkeypatch)
    name, shown = b"rs\xfe", "rs\ufffd"

    def written(path: Path, *rows: str) -> Path:
        path.write_bytes((HEADER + "".join(rows)).encode().replace(b"NAME", name))
        return path

    reference = written(tmp_path / "reference.txt", *(row("1:600", allele, "ENST1", "REVEL=0.7", uploaded="NAME") for allele in "AG"))
    engine = written(
        tmp_path / "engine.txt",
        row("1:600", "A", "ENST1", "REVEL=0.7", uploaded="NAME"),
        row("1:600", "G", "ENST1", "IMPACT=LOW", uploaded="NAME"),
    )
    # A record elsewhere whose identifier is the text the bytes of ``name`` decode to.
    other = f"2\t500\t{shown}\tC\tA,G".encode()
    for label, carriers, aside in (("one", ["C\tA,G"], 1), ("two", ["C\tA", "C\tG"], 0), ("none", [], 0)):
        lines = [b"1\t600\t" + name + b"\t" + alleles.encode() for alleles in carriers]
        run_input = tmp_path / f"{label}.vcf"
        run_input.write_bytes(VCF_HEADER.encode() + b"".join(line + b"\t.\t.\t.\n" for line in (*lines, other)))
        item = entry(score(tmp_path, reference, {"engine": engine}, run_input=run_input), "REVEL")
        assert item["masked"] == {"planted_sibling_value": aside}, label
        assert built[-1][1] == [f"1\t600\t{shown}\t{alleles}" for alleles in carriers], label


def test_input_that_is_no_vcf_or_cannot_be_read_is_refused_with_nothing_left_behind(tmp_path, capsys, monkeypatch, planted):
    plant_record_reader(monkeypatch)
    reference, engine, _ = named_rows(tmp_path)
    table = class_table(tmp_path, CLASS_TABLE)
    empty = tmp_path / "empty.vcf"
    empty.write_text(VCF_HEADER + "1\t5\trs1\tC\n", encoding="utf-8")
    cut = tmp_path / "cut.vcf.gz"
    cut.write_bytes(gzip.compress(vcf(tmp_path / "whole.vcf", *NAMED_RECORDS).read_bytes())[:-12])
    for name, path, message in (
        ("missing", tmp_path / "nowhere.vcf", f"not a file: {tmp_path / 'nowhere.vcf'}"),
        ("empty", empty, f"{empty}: no VCF record, a data line of five or more columns"),
        ("cut", cut, f"{cut}: cannot read the gzip stream: "),
    ):
        code, lines, error = run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={engine}", "--plugins", "REVEL",
            "--class-table", str(table), "--work-dir", str(tmp_path / f"work-{name}"), "--sort-memory", "16M",
            "--report", str(tmp_path / "report.json"), "--residuals", str(tmp_path / name), "--input", str(path),
        )  # fmt: skip
        assert (code, lines) == (2, [])
        assert error.startswith(f"ERROR: [score_plugin_tuples] {message}") and error.count("\n") == 1
        assert not (tmp_path / "report.json").exists()
        for left in (tmp_path / name, tmp_path / f"work-{name}"):
            assert not left.exists() or not any(left.iterdir())


def command(*arguments: str, **environment: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *arguments],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, **environment},
    )


def test_summary_escapes_what_the_stream_cannot_encode(tmp_path):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "Custom=id1;Kéy=1"))
    done = command(
        "--reference", str(reference), "--engine", f"own={reference}", "--custom-source", "Custom=Custom",
        "--work-dir", str(tmp_path / "work"), "--report", str(tmp_path / "report.json"), "--sort-memory", "16M",
        PYTHONIOENCODING="ascii",
    )  # fmt: skip
    assert done.returncode == 0, done.stderr
    assert "keys of the reference that no source claims: K\\xe9y (1)" in done.stdout.splitlines()
    report = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert report["reference"]["unclaimed_keys"] == {"Kéy": 1}


def test_help_prints_the_usage_line_by_line():
    done = command("-h")
    assert done.returncode == 0
    assert "    python3 score_plugin_tuples.py --reference <vep_output.txt[.gz]> \\\n" in done.stdout
    assert "        [--work-dir <dir>] --report <report.json> [--residuals <dir>]\n" in done.stdout


def test_command_line_writes_the_report(tmp_path):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "Custom=id1;Custom_SIG=p"))
    engine = tuples(tmp_path / "engine.tsv", ("1:100", "T", "ENST1", "Custom", "id1"))
    report_path = tmp_path / "reports" / "report.json"
    done = command(
        "--reference", str(reference), "--engine", f"own={reference}", "--engine", f"other={engine}",
        "--custom-source", "Custom=Custom", "--work-dir", str(tmp_path / "work"), "--report", str(report_path),
        "--sort-memory", "16M", "--sort-threads", "1", "--residuals", str(tmp_path / "kept" / "residuals"),
    )  # fmt: skip
    assert done.returncode == 0, done.stderr
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert [engine["id"] for engine in report["engines"]] == ["own", "other"]
    assert report["class_table"]["path"] == str(COMMITTED_CLASS_TABLE)
    assert counts(entry(report, "Custom", engine="own")) == (2, 2, 2, 2, 2, 2)
    assert counts(entry(report, "Custom", engine="other", projection="allele")) == (2, 1, 1, 2, 1, 1)
    assert "Custom\tother\tallele\tperl=2 rust=1 intersection=1" in done.stdout
    assert "text_only_differences=0 residual_rows=1 shared_keys=0\n" in done.stdout
    assert not any((tmp_path / "work").iterdir())
    residuals = tmp_path / "kept" / "residuals"
    assert sorted(path.name for path in residuals.iterdir()) == [
        "Custom.other.allele.tsv",
        "Custom.other.feature.tsv",
        "Custom.own.allele.tsv",
        "Custom.own.feature.tsv",
    ]
    assert residual_lines(residuals, "Custom", "other") == ["reference\t1:100\tT\tENST1\tCustom_SIG\tp"]
    assert residual_lines(residuals, "Custom", "other", "allele") == ["reference\t1:100\tT\tCustom_SIG\tp"]
    assert residual_lines(residuals, "Custom", "own") == []


@pytest.mark.parametrize(
    ("unwritable", "shape"),
    [("--report", "below a file"), ("--report", "a directory"), ("--residuals", "below a file")],
)
def test_unwritable_output_path_is_refused_before_scoring(tmp_path, unwritable, shape):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "Custom=id1"))
    paths = {"--report": tmp_path / "report.json", "--residuals": tmp_path / "residuals"}
    paths[unwritable] = reference / "below-a-file" if shape == "below a file" else tmp_path
    done = command(
        "--reference", str(reference), "--engine", f"own={reference}", "--custom-source", "Custom=Custom",
        "--work-dir", str(tmp_path / "work"), "--report", str(paths["--report"]),
        "--residuals", str(paths["--residuals"]),
    )  # fmt: skip
    assert done.returncode == 2
    assert done.stdout == ""
    lines = done.stderr.splitlines()
    assert len(lines) == 1 and lines[0].startswith("ERROR: [score_plugin_tuples] ")
    assert not (tmp_path / "report.json").exists()
    assert not (tmp_path / "work").exists() or unwritable == "--residuals"
    assert not any((tmp_path / "work").glob("*"))


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["--engine", "own={reference}"], "not a file: {missing}"),
        (["--engine", "{reference}"], "--engine takes ID=PATH"),
        (["--engine", "a={reference}", "--engine", "a={reference}"], "--engine names a twice"),
        (["--engine", "a={reference}", "--custom-source", "Custom"], "--custom-source takes NAME=PREFIX"),
    ],
)
def test_command_line_refusals_are_one_error_line(tmp_path, arguments, message):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "Custom=id1"))
    missing = tmp_path / "missing.txt"
    names = {"reference": str(reference), "missing": str(missing)}
    target = missing if "missing" in message else reference
    done = command(
        "--reference", str(target), "--report", str(tmp_path / "report.json"),
        *(argument.format(**names) for argument in arguments),
    )  # fmt: skip
    assert done.returncode == 2
    assert done.stdout == ""
    lines = done.stderr.splitlines()
    assert len(lines) == 1 and lines[0].startswith("ERROR: [score_plugin_tuples] ")
    assert message.format(**names) in lines[0]
    assert not (tmp_path / "report.json").exists()


def committed_rows() -> list[list[str]]:
    """The committed class table's rows as their cells, its header and blank lines left out."""
    assert COMMITTED_CLASS_TABLE.is_file(), f"{COMMITTED_CLASS_TABLE} is no file"
    lines = COMMITTED_CLASS_TABLE.read_text(encoding="utf-8").splitlines()
    rows = [line.split("\t") for line in lines if line.strip() and not line.startswith("#")]
    assert rows, f"{COMMITTED_CLASS_TABLE} declares no class"
    return rows


def test_committed_class_table_names_classes_of_the_comparator_alone():
    names = [cells[0] for cells in committed_rows()]
    known = {formats.ORDER_ONLY, *COMPARATOR_FIELD_CLASSES}
    assert [name for name in names if name not in known] == []


def test_committed_class_table_is_taken_whole_by_the_scorers_reader(monkeypatch):
    rows = committed_rows()
    # The comparator's own classes and the registry's own plugins, in place of the planted ones.
    monkeypatch.undo()
    sources = spt.Sources(list(compare_plugin_outputs.ALL_PLUGIN_NAMES), {})
    classes = spt.Classes(COMMITTED_CLASS_TABLE, None, {}, sources)
    assert list(classes.declared) == [cells[0] for cells in rows]
    for name, plugins, groups, *further in rows:
        fields = [tuple(group.split(",")) for group in groups.split(";")]
        assert classes.declared[name] == dict(zip(plugins.split(";"), fields, strict=True)), name
        assert classes.further[name] == further, name


def test_committed_class_table_declares_no_class_twice():
    times = Counter(cells[0] for cells in committed_rows())
    assert [name for name, count in times.items() if count > 1] == []


def test_committed_class_table_marks_every_row_met_or_not_met():
    rows = committed_rows()
    header = COMMITTED_CLASS_TABLE.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert header[:4] == ["#class", "plugin", "field names", "met or not met"]
    for cells in rows:
        assert len(cells) == len(header), cells[0]
        assert cells[3] in ("met", "not met"), cells[0]


def default_table_run(tmp_path: Path, *arguments: str) -> dict:
    """The report of a command-line run that scores one LOEUF row against itself under the
    comparator's own classes and the registry's own plugins."""
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "LOEUF=0.31"))
    tool = ["--pairing-tool", str(built_pairing_tool())] if THROUGH_PAIRING_TOOL else []
    done = command(
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "LOEUF",
        "--work-dir", str(tmp_path / "work"), "--sort-memory", "16M", "--report", str(tmp_path / "report.json"),
        *tool, *arguments,
    )  # fmt: skip
    assert (done.returncode, done.stderr) == (0, ""), done.stderr
    return json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))


def applied_classes(report: dict) -> set[str]:
    return {item["class"] for item in report["classes"] if item["applied"]}


OTHER_CLASS_TABLE = "loeuf_zero_minimum\tLOEUF\tLOEUF\norder_only\tNearestExonJB\tNearestExonJB\tnot met\t-\n"


def test_command_line_reads_the_committed_class_table_by_default(tmp_path):
    report = default_table_run(tmp_path)
    assert report["class_table"] is not None, "the run read no class table"
    assert report["class_table"]["path"] == str(COMMITTED_CLASS_TABLE)
    assert report["class_table"]["sha256"] == hashlib.sha256(COMMITTED_CLASS_TABLE.read_bytes()).hexdigest()
    assert [item["class"] for item in report["classes"]] == [cells[0] for cells in committed_rows()]


def test_command_line_class_table_and_classes_decide_over_the_default(tmp_path):
    other = class_table(tmp_path, OTHER_CLASS_TABLE)
    report = default_table_run(tmp_path, "--class-table", str(other), "--classes", "order_only")
    assert report["class_table"]["path"] == str(other)
    assert [item["class"] for item in report["classes"]] == ["loeuf_zero_minimum", "order_only"]
    assert applied_classes(report) == {"order_only"}
    assert applied_classes(default_table_run(tmp_path, "--class-table", str(other), "--classes", "")) == set()
    rows = committed_rows()
    named = {rows[0][0], next((cells[0] for cells in rows if cells[3] == "not met"), rows[-1][0])}
    report = default_table_run(tmp_path, "--classes", ",".join(sorted(named)))
    assert report["class_table"]["path"] == str(COMMITTED_CLASS_TABLE)
    assert applied_classes(report) == named
    assert applied_classes(default_table_run(tmp_path, "--classes", "")) == set()


def test_command_line_applies_by_default_every_row_but_those_that_read_not_met(tmp_path):
    met = {cells[0] for cells in committed_rows() if cells[3] != "not met"}
    assert applied_classes(default_table_run(tmp_path)) == met
    other = class_table(tmp_path, OTHER_CLASS_TABLE)
    assert applied_classes(default_table_run(tmp_path, "--class-table", str(other))) == {"loeuf_zero_minimum"}


def on_both_paths(tmp_path: Path, tool: Path, reference: Path, engines: dict[str, Path], **options) -> tuple[dict, dict]:
    """Scores on the script's own path and through the pairing tool, on one thread and on
    several, and returns the report and the residual files, by name, once every run wrote
    the same bytes of both."""
    written = []
    for name, binary, threads in (("own", None, None), ("one", tool, 1), ("several", tool, 5)):
        residuals = tmp_path / f"residuals-{name}"
        report = score(tmp_path, reference, engines, residuals=residuals, tool=binary, threads=threads, **options)
        files = {path.name: path.read_bytes() for path in sorted(residuals.iterdir())}
        written.append((json.dumps(report, indent=2), files))
        shutil.rmtree(residuals)
        assert not any((tmp_path / "work").iterdir())
    assert written[1] == written[0] and written[2] == written[0]
    return json.loads(written[0][0]), written[0][1]


def test_pairing_tool_writes_the_report_and_residual_files_of_the_scripts_own_path(tmp_path, pairing_tool, planted):
    rng = random.Random(1162)
    locations = ["1:5", "1:50", "1:5-06", "1:5-6", "10:5", "MT:+5", "21:²", " 2:5 "]
    alleles = ["A", "-", "AT", " T"]
    features = ["ENST1", "ENST10", "ENST2 ", "-"]
    written = {
        "PLANT_score": ["0.1", "0.10", "7", "7.0000001", "-", ""],
        "PLANT_label": ["a", "b", "a%3Bb%3Dc", "x\x01y", "%3d"],
        "LOEUF": ["0", "0.31"],
        "MaxEntScan_diff": ["1.5", "-2", "n"],
        "GWAS_pmid": ["17", "18"],
        "NearestExonJB": ["a,b", "b,a", "a"],
        "5UTR_consequence": ["gained", "lost"],
        "REVEL": ["0.2", "0.4"],
        "IMPACT": ["LOW"],
        "Custom_SIG": ["p", "q"],
        "Stray": ["1"],
        "Bare": ["yes"],
    }

    def side(count: int) -> tuple[list[str], list[tuple]]:
        rows, listed = [], []
        for _ in range(count):
            key = (rng.choice(locations), rng.choice(alleles), rng.choice(features))
            pairs = [(name, rng.choice(values)) for name, values in written.items() if rng.random() < 0.4]
            extra = ";".join("Bare" if name == "Bare" else f"{name}={value}" for name, value in pairs)
            rows.append(row(*key, extra or "-", cdna=rng.choice("12"), uploaded=rng.choice(["v", "1_5_A/T"])))
            listed += [(*key, name, formats.decode_extra(value)) for name, value in pairs if name != "Bare"]
        return rows, listed

    reference = output(tmp_path / "reference.txt", *side(300)[0])
    engine_rows, _ = side(300)
    plain = tmp_path / "engine.txt"
    plain.write_text(HEADER + "".join(engine_rows[:150]) + "\r\n##a comment\r\n" + "".join(engine_rows[150:]))
    packed = tmp_path / "engine.txt.gz"
    packed.write_bytes(gzip.compress(plain.read_bytes()[:900]) + gzip.compress(plain.read_bytes()[900:]) + bytes(40))
    listed = tuples(tmp_path / "engine.tsv", *side(300)[1])
    engines = {"output": plain, "packed": packed, "listed": listed}
    calls = []
    for table in (CLASS_TABLE, None):
        planted.clear()
        report, files = on_both_paths(tmp_path, pairing_tool, reference, engines, custom={"Custom": "Custom"}, table=table)
        per_run = len(planted) // 3
        assert planted[:per_run] == planted[per_run : 2 * per_run] == planted[2 * per_run :]
        calls.append(per_run)
        assert len(files) == 2 * len(engines) * len(report["sources"])
        assert sum(map(len, files.values())) > 0
        assert report["engines"][2]["format"] == "tuples" and report["engines"][2]["tuples_outside_every_source"] > 0
        assert report["engines"][1]["sha256"] == hashlib.sha256(packed.read_bytes()).hexdigest()
        assert report["engines"][0]["unclaimed_keys"] == report["engines"][1]["unclaimed_keys"]
        assert list(report["engines"][0]["unclaimed_keys"]) == ["Stray"]
        assert all(item["shared_keys"] > 0 for item in report["scores"] if item["source"] == "Planted")
    assert calls[0] > 0 and calls[1] == 0
    masked = [item["masked"] for item in json.loads(json.dumps(report))["scores"]]
    assert masked and not any(masked)


def test_pairing_tool_leaves_a_missing_an_extra_and_a_changed_value_to_the_rules(tmp_path, pairing_tool):
    reference = output(
        tmp_path / "reference.txt",
        row("1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=kept"),
        row("1:200", "T", "ENST1", "PLANT_score=0.25;PLANT_label=missing"),
        row("1:300", "T", "ENST1", "PLANT_score=1;PLANT_label=before"),
        row("1:400", "T", "ENST1", "PLANT_score=0.1"),
    )
    engine = output(
        tmp_path / "engine.txt",
        row("1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=kept;PlantedDynamic_a=extra"),
        row("1:200", "T", "ENST1", "PLANT_score=0.25"),
        row("1:300", "T", "ENST1", "PLANT_score=1;PLANT_label=after"),
        row("1:400", "T", "ENST1", "PLANT_score=0.10"),
    )
    report, files = on_both_paths(tmp_path, pairing_tool, reference, {"engine": engine}, table=None)
    planted_entry = entry(report, "Planted")
    assert counts(planted_entry) == (7, 6, 5, 7, 6, 5)
    assert (planted_entry["text_only_differences"], planted_entry["residual_rows"]) == (1, 3)
    assert files["Planted.engine.feature.tsv"].decode().splitlines() == [
        "reference\t1:200\tT\tENST1\tPLANT_label\tmissing",
        "engine\t1:300\tT\tENST1\tPLANT_label\tafter",
        "reference\t1:300\tT\tENST1\tPLANT_label\tbefore",
    ]
    assert counts(entry(report, "PlantedDynamic")) == (0, 1, 0, 0, 1, 0)
    assert files["PlantedDynamic.engine.feature.tsv"] == b"engine\t1:100\tT\tENST1\tPlantedDynamic_a\textra\n"
    assert files["PlantedDynamic.engine.allele.tsv"] == b"engine\t1:100\tT\tPlantedDynamic_a\textra\n"


def test_pairing_tool_refusal_is_the_scripts_own_error(tmp_path, pairing_tool, monkeypatch):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    narrow = tmp_path / "narrow.txt"
    narrow.write_text(HEADER + REFERENCE_ROWS[0] + "1:100\tT\tENST1\n", encoding="utf-8")
    odd = tmp_path / "odd.txt"
    odd.write_text(HEADER + REFERENCE_ROWS[0] + "#Uploaded_variation\tLocation\n", encoding="utf-8")
    wide = output(tmp_path / "wide.txt", row("1:100", "T", "ENST1", "One=1;Two=2;Three=3;PLANT_score=1"))
    monkeypatch.setattr(spt, "RESOLVED_LIMIT", 2)

    def refusal(engine: Path, tool: Path | None) -> str:
        with pytest.raises(spt.ScorerError) as raised:
            score(tmp_path, reference, {"engine": reference, "later": engine}, tool=tool)
        assert not any((tmp_path / "work").iterdir())
        return str(raised.value)

    for engine, message in (
        (narrow, "narrow.txt: line 3 has 3 columns; a default-format row has 14, a tuple file line 5"),
        (odd, "odd.txt: line 3, the column header, names 2 columns, not 14"),
        (wide, "wide.txt: more than 2 distinct keys that no source claims"),
    ):
        assert refusal(engine, pairing_tool) == refusal(engine, None) == f"{engine.parent}/{message}"
    # A path can carry a line boundary, can read as the tool's help flag and, where the file
    # system takes them, can carry bytes that are no text; the tool can be named without a directory.
    tool = pairing_tool.absolute()
    with monkeypatch.context() as within:
        within.setenv("VEP_PLUGIN_PAIRS", str(tool))
        within.chdir(tmp_path)
        Path("tool-here").symlink_to(tool)
        for name in ("nar\x0crow\N{LINE SEPARATOR}.txt", "-h", "nar\udcffrow.txt"):
            try:
                Path(name).write_bytes(narrow.read_bytes())
            except OSError:
                continue
            message = f"{name}: line 3 has 3 columns; a default-format row has 14, a tuple file line 5"
            assert refusal(Path(name), Path("tool-here")) == refusal(Path(name), None) == message
    registry = compare_plugin_outputs.plugin_of_field

    def refusing(field, plugins):
        if field == "PLANT_label":
            raise ValueError("the run's flags make PLANT_label a key of both engines")
        return registry(field, plugins)

    monkeypatch.setattr(compare_plugin_outputs, "plugin_of_field", refusing)
    assert refusal(reference, pairing_tool) == refusal(reference, None)
    assert refusal(reference, None) == "the run's flags make PLANT_label a key of both engines"


def test_pairing_tool_that_does_not_run_or_ends_without_its_error_line_is_named(tmp_path):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    with pytest.raises(spt.ScorerError, match="cannot run the pairing tool .*no-such-tool: No such file"):
        score(tmp_path, reference, {"engine": reference}, tool=tmp_path / "no-such-tool")
    failing = tmp_path / "failing-tool"
    failing.write_text("#!/bin/sh\necho out of luck >&2\nexit 7\n", encoding="utf-8")
    failing.chmod(0o755)
    with pytest.raises(spt.ScorerError, match="the pairing tool exited 7: out of luck$"):
        score(tmp_path, reference, {"engine": reference}, tool=failing)
    assert not any((tmp_path / "work").iterdir())


def test_pairing_tool_counts_that_its_remainder_does_not_bear_out_are_refused(tmp_path, pairing_tool, monkeypatch):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    handed = spt.Pairing.engine

    def miscounted(self, index):
        against, engine, counted = handed(self, index)
        counted["feature"]["values"][b"PLANT_score".hex()][0] += 1
        return against, engine, counted

    monkeypatch.setattr(spt.Pairing, "engine", miscounted)
    with pytest.raises(spt.ScorerError, match="the pairing tool's counts of PLANT_score and its remainder disagree"):
        score(tmp_path, reference, {"engine": reference}, tool=pairing_tool)
    assert not any((tmp_path / "work").iterdir())


def test_interrupt_while_the_pairing_tool_runs_stops_it_and_leaves_nothing_behind(
    tmp_path, capsys, monkeypatch, pairing_tool
):
    reference = output(tmp_path / "reference.txt", *REFERENCE_ROWS)
    started = []
    popen = subprocess.Popen

    def recorded(*arguments, **options):
        started.append(popen(*arguments, **options))
        return started[-1]

    def interrupted(self, field):
        os.kill(os.getpid(), signal.SIGINT)

    monkeypatch.setattr(spt.subprocess, "Popen", recorded)
    monkeypatch.setattr(spt.Sources, "of", interrupted)
    code, lines, error = run_main(
        capsys,
        "--reference", str(reference), "--engine", f"engine={reference}", "--plugins", "Planted",
        "--work-dir", str(tmp_path / "work"), "--report", str(tmp_path / "report.json"),
        "--residuals", str(tmp_path / "residuals"), "--pairing-tool", str(pairing_tool),
    )  # fmt: skip
    assert (code, lines) == (130, [])
    assert error == "ERROR: [score_plugin_tuples] interrupted before the report was written\n"
    assert [run.args[0] for run in started] == [str(pairing_tool)] and started[0].returncode == -signal.SIGKILL
    assert not (tmp_path / "report.json").exists()
    assert not any((tmp_path / "work").iterdir()) and not any((tmp_path / "residuals").iterdir())


def test_command_line_names_the_pairing_tool_and_its_threads(tmp_path, pairing_tool):
    reference = output(tmp_path / "reference.txt", row("1:100", "T", "ENST1", "Custom=id1;Custom_SIG=p;Other=1"))
    engine = output(tmp_path / "engine.txt", row("1:100", "T", "ENST1", "Custom=id1;Custom_SIG=q;Other=1"))
    written = []
    for name, extra in (("own", []), ("tool", ["--pairing-tool", str(pairing_tool), "--pairing-threads", "2"])):
        done = command(
            "--reference", str(reference), "--engine", f"engine={engine}", "--custom-source", "Custom=Custom",
            "--work-dir", str(tmp_path / "work"), "--report", str(tmp_path / f"{name}.json"),
            "--residuals", str(tmp_path / name), *extra,
        )  # fmt: skip
        assert (done.returncode, done.stderr) == (0, "")
        files = {path.name: path.read_bytes() for path in sorted((tmp_path / name).iterdir())}
        written.append((done.stdout, (tmp_path / f"{name}.json").read_bytes(), files))
        assert not any((tmp_path / "work").iterdir())
    assert written[0] == written[1]
    assert "keys of the reference that no source claims: Other (1)\n" in written[1][0]
    assert written[1][2]["Custom.engine.feature.tsv"] == (
        b"reference\t1:100\tT\tENST1\tCustom_SIG\tp\nengine\t1:100\tT\tENST1\tCustom_SIG\tq\n"
    )
    missing = command(
        "--reference", str(reference), "--engine", f"engine={engine}", "--custom-source", "Custom=Custom",
        "--report", str(tmp_path / "report.json"), "--pairing-tool", str(tmp_path / "no-such-tool"),
    )  # fmt: skip
    assert (missing.returncode, missing.stdout) == (2, "")
    assert missing.stderr == f"ERROR: [score_plugin_tuples] not a file: {tmp_path / 'no-such-tool'}\n"


SIBLING_ALLELE_CASES = (
    (formats.ALPHAMISSENSE_SIBLING_ALLELE_ROWS, "AlphaMissense", ("am_class", "am_pathogenicity"),
     ("likely_benign", "0.1079"), ("likely_pathogenic", "0.9")),
    (formats.REVEL_SIBLING_ALLELE_ROWS, "REVEL", ("REVEL",), ("0.025",), ("0.027",)),
    (formats.SPLICEVAULT_SIBLING_ALLELE_ROWS, "SpliceVault", ("SpliceVault_site_type", "SpliceVault_site_pos"),
     ("donor", "5"), ("acceptor", "7")),
)  # fmt: skip
SPELT_RECORD = "3_7000200_A/C/T"
NAMED_RECORD = ("3", "7000200", "rs77", "A", "C,T")


def two_alternate_rows(path: Path, first: str, second: str, uploaded: str = SPELT_RECORD, location: str = "3:7000200") -> Path:
    """A made-up record of two alternates, C and T, on one transcript, both giving one
    amino-acid change: its two rows, ``first`` and ``second`` the Extra of each allele's."""
    columns = "ENSG1\tENST1\tTranscript\tmissense_variant\t-\t-\t-\tM/L\t-\t-"
    rows = (f"{uploaded}\t{location}\t{allele}\t{columns}\t{extra}\n" for allele, extra in (("C", first), ("T", second)))
    return output(path, *rows)


def plugin_extra(names: tuple[str, ...], values: tuple[str, ...]) -> str:
    return ";".join(["IMPACT=MODERATE", *(f"{name}={value}" for name, value in zip(names, values))])


@pytest.mark.parametrize("klass, plugin, names, copied, own", SIBLING_ALLELE_CASES, ids=[case[1] for case in SIBLING_ALLELE_CASES])
def test_sibling_allele_class_is_handed_the_engines_rows_of_the_same_record(tmp_path, monkeypatch, klass, plugin, names, copied, own):
    """A made-up record of two alternates on one transcript, both giving one amino-acid
    change: the reference writes the first allele's values on both alleles' rows, the engine
    each allele's own. The comparator's own predicate, unplanted, takes the second allele's
    values as its sibling's only where the engine's row of the first allele carries them too,
    which it reads from the engine's rows of the record: so the scorer, handing them, sets the
    second allele's values aside, under an ``Uploaded_variation`` that spells the alleles and
    under an identifier one line of the input carries. Where the engine's first row carries
    other values, and under an identifier without the input, nothing is set aside. Each
    verdict is the format comparison's own on the same two outputs."""
    plant_registry(monkeypatch, {plugin: list(names)})
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    table = f"{klass}\t{plugin}\t{','.join(names)}\tmet\t-\n"
    run_input = vcf(tmp_path / "input.vcf", NAMED_RECORD)
    line = "\t".join((*NAMED_RECORD, ".", ".", "."))
    held = len(names)
    every, half = 2 * held, held
    cases = (
        ("spelt", SPELT_RECORD, copied, False, (every, every, half, half, half, half), 1),
        ("first row differs", SPELT_RECORD, own, False, (every, every, 0, every, every, 0), 0),
        ("named", NAMED_RECORD[2], copied, True, (every, every, half, half, half, half), 1),
        ("named, no input", NAMED_RECORD[2], copied, False, (every, every, half, every, every, half), 0),
    )
    for label, uploaded, first, with_input, counted, aside in cases:
        reference = two_alternate_rows(tmp_path / "reference.txt", plugin_extra(names, copied), plugin_extra(names, copied), uploaded)
        engine = two_alternate_rows(tmp_path / "engine.txt", plugin_extra(names, first), plugin_extra(names, own), uploaded)
        report = score(
            tmp_path, reference, {"engine": engine}, table=table, plugins=[plugin], run_input=run_input if with_input else None
        )
        for projection in spt.PROJECTIONS:
            item = entry(report, plugin, projection=projection)
            assert (counts(item), item["masked"]) == (counted, {klass: aside * held}), (label, projection)
            assert item["residual_rows"] == counted[3] + counted[4] - 2 * counted[5], (label, projection)
        documented = formats.Documented.from_manifest(None, {klass})
        compared = formats.compare_entries(
            formats.parse_default(reference.read_text(encoding="utf-8")),
            formats.parse_default(engine.read_text(encoding="utf-8")),
            documented,
            input_lines=[line] if with_input else None,
        )
        item = entry(report, plugin)
        accepted = sum(compared["columns"].get(name, {}).get("masked", {}).get(klass, 0) for name in names)
        assert item["masked"] == {klass: accepted} and documented.applied.get(klass, 0) == aside, label
        assert (item["residual_rows"] == 0) == (compared["verdict"] == "green"), label


def test_pairing_tool_hands_the_engines_rows_of_a_record_to_the_sibling_allele_classes(tmp_path, monkeypatch, pairing_tool):
    """The engine's row of a record's other allele pairs as equal text, so it is no part of
    the remainder: the pairing tool hands it over with the engine's rows of the Location, and
    both paths set the same values aside under each sibling-allele class, on a record whose
    ``Uploaded_variation`` spells its alleles and on one under an identifier, and leave the
    same residual files: the rows of a record whose first allele's row differs as well."""
    plant_registry(monkeypatch, {plugin: list(names) for _, plugin, names, _, _ in SIBLING_ALLELE_CASES})
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    table = "".join(f"{klass}\t{plugin}\t{','.join(names)}\tmet\t-\n" for klass, plugin, names, _, _ in SIBLING_ALLELE_CASES)
    copied = ";".join(plugin_extra(names, values).partition(";")[2] for _, _, names, values, _ in SIBLING_ALLELE_CASES)
    own = ";".join(plugin_extra(names, values).partition(";")[2] for _, _, names, _, values in SIBLING_ALLELE_CASES)
    named = ("3", "7000500", "rs78", "A", "C,T")
    records = ((SPELT_RECORD, "3:7000200"), (named[2], "3:7000500"), ("3_7000900_A/C/T", "3:7000900"))

    def rows(path: Path, first_of_the_last: str, second: str) -> Path:
        written = [
            two_alternate_rows(path, first_of_the_last if location == "3:7000900" else copied, second, uploaded, location)
            .read_text(encoding="utf-8")
            .removeprefix(HEADER)
            for uploaded, location in records
        ]
        return output(path, *written)

    reference = rows(tmp_path / "reference.txt", copied, copied)
    engine = rows(tmp_path / "engine.txt", own, own)
    run_input = vcf(tmp_path / "input.vcf", named)
    report, files = on_both_paths(tmp_path, pairing_tool, reference, {"engine": engine}, table=table, run_input=run_input)
    for klass, plugin, names, _, _ in SIBLING_ALLELE_CASES:
        held = len(names)
        item = entry(report, plugin)
        assert (counts(item), item["masked"]) == ((6 * held, 6 * held, 2 * held, 4 * held, 4 * held, 2 * held), {klass: 2 * held}), plugin
        left = files[f"{plugin}.engine.feature.tsv"].splitlines()
        assert len(left) == item["residual_rows"] == 4 * held
        assert all(line.split(b"\t")[1] == b"3:7000900" for line in left)


REVEL_OTHER = formats.REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS
REVEL_CLASSES_TABLE = f"{formats.REVEL_SIBLING_ALLELE_ROWS}\tREVEL\tREVEL\tmet\t-\n{REVEL_OTHER}\tREVEL\tREVEL\tmet\t-\n"
REVEL_NOTICE = "no --revel-table: a class of REVEL that reads the table sets nothing aside"


def revel_table_file(directory: Path, *rows: tuple[str, str, str, str, str, str]) -> Path:
    """A made-up REVEL table on chromosome 3, tabix-indexed on its GRCh38 column and each line
    a BGZF block of its own: (position, alternate allele, aaref, aaalt, score, transcript ids)
    each, the reference allele A."""
    def block(data: bytes) -> bytes:
        deflater = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = deflater.compress(data) + deflater.flush()
        return (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", len(body) + 25)
                + body + struct.pack("<II", zlib.crc32(data), len(data)))

    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "revel.tsv.gz"
    header = block(b"#chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid\n")
    body = b"".join(block(f"3\t999\t{position}\tA\t{alternate}\t{aaref}\t{aaalt}\t{score}\t{ids}\n".encode())
                    for position, alternate, aaref, aaalt, score, ids in rows)
    path.write_bytes(header + body + block(b""))
    index = b"TBI\x01" + struct.pack("<8i", 1, 0, 1, 3, 3, ord("#"), 0, 2) + b"3\0"
    index += struct.pack("<iIi2Qi", 1, 0, 1, len(header) << 16, (len(header) + len(body)) << 16, 0)
    Path(f"{path}.tbi").write_bytes(gzip.compress(index))
    return path


def revel_row(uploaded: str, location: str, allele: str, feature: str, score: str) -> str:
    """A missense row giving M/L, with ``score`` as its REVEL or without the field."""
    columns = f"ENSG1\t{feature}\tTranscript\tmissense_variant\t-\t-\t-\tM/L\t-\t-"
    return f"{uploaded}\t{location}\t{allele}\t{columns}\tIMPACT=MODERATE" + (f";REVEL={score}" if score else "") + "\n"


def plant_revel(monkeypatch) -> None:
    plant_registry(monkeypatch, {"REVEL": ["REVEL"]})
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)


def test_revel_table_is_handed_to_the_classes_of_revel_which_set_nothing_aside_without_it(tmp_path, monkeypatch):
    """Three made-up records on one transcript, each row giving M/L, against a made-up table
    whose row at each position lists the transcript and names another reference amino acid.
    The reference writes that row's score and the engine nothing. With ``--revel-table`` the
    comparator's own predicate, unplanted, sets each value aside: on a record whose
    ``Uploaded_variation`` spells its alleles, on one under an identifier one input line
    carries, and on a one-base row under an identifier with no input, which the class reads
    without its record. Without the table nothing is set aside. The report names the table
    with its digest, and each verdict is the format comparison's own on the same outputs."""
    plant_revel(monkeypatch)
    table = revel_table_file(tmp_path / "table", ("7000200", "T", "V", "L", "0.300", "ENST1"), ("7000500", "T", "V", "L", "0.310", "ENST1"),
                             ("7000900", "T", "V", "L", "0.320", "ENST1"), ("7001200", "T", "M", "L", "0.330", "ENST1"))
    records = (("3_7000200_A/T", "3:7000200", "0.300"), ("rs78", "3:7000500", "0.310"), ("rs79", "3:7000900", "0.320"),
               ("3_7001200_A/T", "3:7001200", "0.330"))
    reference = output(tmp_path / "reference.txt", *(revel_row(uploaded, location, "T", "ENST1", score) for uploaded, location, score in records))
    engine = output(tmp_path / "engine.txt", *(revel_row(uploaded, location, "T", "ENST1", "") for uploaded, location, _ in records))
    run_input = vcf(tmp_path / "input.vcf", ("3", "7000500", "rs78", "A", "T"))
    without = score(tmp_path, reference, {"engine": engine}, table=REVEL_CLASSES_TABLE, plugins=["REVEL"], run_input=run_input)
    assert without["revel_table"] is None
    for projection in spt.PROJECTIONS:
        item = entry(without, "REVEL", projection=projection)
        assert (counts(item), item["residual_rows"]) == ((4, 0, 0, 4, 0, 0), 4), projection
        assert item["masked"] == {formats.REVEL_SIBLING_ALLELE_ROWS: 0, REVEL_OTHER: 0}
    report = score(tmp_path, reference, {"engine": engine}, table=REVEL_CLASSES_TABLE, plugins=["REVEL"], run_input=run_input, revel=table)
    assert report["revel_table"] == {"path": str(table), "sha256": hashlib.sha256(table.read_bytes()).hexdigest(), "no_match": False}
    for projection in spt.PROJECTIONS:
        item = entry(report, "REVEL", projection=projection)
        assert (counts(item), item["residual_rows"]) == ((4, 0, 0, 1, 0, 0), 1), projection
        assert item["masked"] == {formats.REVEL_SIBLING_ALLELE_ROWS: 0, REVEL_OTHER: 3}
    documented = formats.Documented.from_manifest(None, {REVEL_OTHER})
    compared = formats.compare_entries(
        formats.parse_default(reference.read_text(encoding="utf-8")), formats.parse_default(engine.read_text(encoding="utf-8")),
        documented, input_lines=["3\t7000500\trs78\tA\tT\t.\t.\t."], revel_table=formats.RevelTable(table),
    )
    assert compared["columns"]["REVEL"]["masked"] == {REVEL_OTHER: 3} and compared["columns"]["REVEL"]["differing"] == 1
    # A table whose rows list another transcript explains the values only for a run that set
    # the plugin's no_match.
    unlisted = revel_table_file(tmp_path / "unlisted", *((position, "T", "V", "L", value, "ENST9") for _, position, value in
                                                         ((u, loc.partition(":")[2], s) for u, loc, s in records[:3])))
    matched = score(tmp_path, reference, {"engine": engine}, table=REVEL_CLASSES_TABLE, plugins=["REVEL"], revel=unlisted)
    assert entry(matched, "REVEL")["masked"][REVEL_OTHER] == 0
    unmatched = score(tmp_path, reference, {"engine": engine}, table=REVEL_CLASSES_TABLE, plugins=["REVEL"], revel=unlisted, revel_no_match=True)
    assert unmatched["revel_table"]["no_match"] is True and entry(unmatched, "REVEL")["masked"][REVEL_OTHER] == 3


def test_both_paths_set_the_score_of_another_alternates_row_aside_under_the_table(tmp_path, monkeypatch, pairing_tool):
    """A made-up record of two alternates, C and T, on one transcript: the table's row of C
    names the T row's alternate amino acid beside another reference amino acid and comes
    first, so the reference writes its score on the T row; the engine writes nothing there,
    and neither writes a value on the C row. The class reads the record's alternates, which
    the name spells or the one input line of the identifier states, and the script's own
    path and the pairing tool leave the same report and residual files."""
    plant_revel(monkeypatch)
    table = revel_table_file(tmp_path / "table", ("7000200", "C", "V", "L", "0.016", "ENST1"), ("7000200", "T", "I", "L", "0.029", "ENST1"),
                             ("7000500", "C", "V", "L", "0.017", "ENST1"), ("7000500", "T", "I", "L", "0.028", "ENST1"),
                             ("7000900", "C", "M", "L", "0.018", "ENST1"))
    named = ("3", "7000500", "rs78", "A", "C,T")
    records = ((SPELT_RECORD, "3:7000200", "0.016"), (named[2], "3:7000500", "0.017"), ("3_7000900_A/C/T", "3:7000900", "0.018"))

    def rows(path: Path, with_scores: bool) -> Path:
        return output(path, *(
            revel_row(uploaded, location, allele, "ENST1", score if with_scores and allele == "T" else "")
            for uploaded, location, score in records for allele in "CT"
        ))

    reference, engine = rows(tmp_path / "reference.txt", True), rows(tmp_path / "engine.txt", False)
    run_input = vcf(tmp_path / "input.vcf", named)
    report, files = on_both_paths(tmp_path, pairing_tool, reference, {"engine": engine}, table=REVEL_CLASSES_TABLE,
                                  plugins=["REVEL"], run_input=run_input, revel=table)
    item = entry(report, "REVEL")
    assert (counts(item), item["masked"]) == ((3, 0, 0, 1, 0, 0), {formats.REVEL_SIBLING_ALLELE_ROWS: 0, REVEL_OTHER: 2})
    assert files["REVEL.engine.feature.tsv"] == b"reference\t3:7000900\tT\tENST1\tREVEL\t0.018\n"


def test_summary_says_that_a_class_of_revel_sets_nothing_aside_without_the_table(tmp_path, capsys, monkeypatch):
    plant_revel(monkeypatch)
    table = revel_table_file(tmp_path / "table", ("7000200", "T", "V", "L", "0.300", "ENST1"))
    reference = output(tmp_path / "reference.txt", revel_row("3_7000200_A/T", "3:7000200", "T", "ENST1", "0.300"))
    engine = output(tmp_path / "engine.txt", revel_row("3_7000200_A/T", "3:7000200", "T", "ENST1", ""))
    classes = class_table(tmp_path, REVEL_CLASSES_TABLE)

    def run(name, *more):
        return run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={engine}", "--plugins", "REVEL",
            "--class-table", str(classes), "--work-dir", str(tmp_path / f"work-{name}"), "--sort-memory", "16M",
            "--report", str(tmp_path / "report.json"), "--residuals", str(tmp_path / name), *more,
        )  # fmt: skip

    code, lines, error = run("without")
    assert (code, error) == (0, "")
    assert lines[-1] == REVEL_NOTICE
    assert f"REVEL\tengine\tfeature\tclass {REVEL_OTHER} set aside 0 (none) [met, -]" in lines
    code, lines, error = run("with", "--revel-table", str(table))
    assert (code, error) == (0, "")
    assert REVEL_NOTICE not in lines
    assert f"REVEL\tengine\tfeature\tclass {REVEL_OTHER} set aside 1 [met, -]" in lines
    code, lines, error = run("no-match", "--revel-table", str(table), "--revel-no-match")
    assert (code, error) == (0, "")
    assert json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))["revel_table"]["no_match"] is True
    code, lines, error = run("missing", "--revel-table", str(tmp_path / "nowhere.tsv.gz"))
    assert (code, lines) == (2, [])
    assert error == f"ERROR: [score_plugin_tuples] not a file: {tmp_path / 'nowhere.tsv.gz'}\n"
    Path(f"{table}.tbi").unlink()
    code, lines, error = run("unindexed", "--revel-table", str(table))
    assert (code, lines) == (2, [])
    assert error.startswith(f"ERROR: [score_plugin_tuples] --revel-table: cannot read the tabix index {table}.tbi")


SPLICEAI_MULTI = "spliceai_multi_entry_cell"
SPLICEAI_CLASSES_TABLE = f"{SPLICEAI_MULTI}\tSpliceAI\tSpliceAI_pred,SpliceAI_pred_*,SpliceAI_cutoff\tmet\t-\n"
SPLICEAI_NOTICE = f"no --spliceai-table: {SPLICEAI_MULTI} sets nothing aside"
SPLICEAI_CACHE_NOTICE = f"no --json-cache: {SPLICEAI_MULTI} sets aside no row that writes no SYMBOL"
# Made-up scores of two symbols at the deletion 3 7000199 GCA G, whose rows stand at
# 3:7000200-7000201 with the allele written -.
SPLICEAI_SCORES = {"GENEA": "0.15|0.08|0.06|0.17|-6|-41|-17|-9", "GENEB": "0.02|0.12|0.10|0.05|-25|-23|37|17"}
SPLICEAI_WHOLE = f"GENEA|{SPLICEAI_SCORES['GENEA']},G|GENEB|{SPLICEAI_SCORES['GENEB']}"
# The POS, REF and ALT of a made-up substitution, whose rows stand at 3:7000300 with the allele C.
SPLICEAI_BASE = ("7000300", "A", "C")


def spliceai_score_file(directory: Path, join: bool = True, variant: tuple[str, str, str] = ("7000199", "GCA", "G"),
                        records: list[tuple[str, ...]] | None = None, lines: list[str] | None = None) -> Path:
    """A made-up indel score file in the fixture generator's layout, tabix-indexed and each
    line a BGZF block of its own: one record per symbol at 3 7000199 GCA G, or, under
    ``join``, the one record the generator's ``--join-entries`` form writes there.
    ``variant`` is the records' POS, REF and ALT: under a REF and an ALT of one base each
    the file is the single-base one. ``records`` names the file's record lines instead, each
    the entries of one line at the variant, a symbol and its eight values, and ``lines`` gives
    the file's record lines as they stand, in the order of their positions."""
    sys.path.insert(0, str(HERE.parent / "fixtures"))
    import plugin_score_fixtures as generator

    def block(data: bytes) -> bytes:
        deflater = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = deflater.compress(data) + deflater.flush()
        return (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", len(body) + 25)
                + body + struct.pack("<II", zlib.crc32(data), len(data)))

    position, reference, alternate = variant
    if lines is not None:
        records, join = None, False
    else:
        lines = ["\t".join(["3", position, ".", reference, alternate, ".", ".", f"SpliceAI={alternate}|{symbol}|{scores}"]) for symbol, scores in SPLICEAI_SCORES.items()]
    if records is not None:
        lines = ["\t".join(["3", position, ".", reference, alternate, ".", ".", "SpliceAI=" + ",".join(f"{alternate}|{entry}" for entry in entries)]) for entries in records]
    elif join:
        lines, _ = generator.in_record_forms([(f"3:{position}", line.split("\t")) for line in lines], generator.RecordForms(join_entries={f"3:{position}": None}))
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / ("snv.synthetic.vcf.gz" if len(reference) == len(alternate) == 1 else "indel.synthetic.vcf.gz")
    header = block(("\n".join(generator.SPLICEAI_HEADER) + "\n").encode())
    body = b"".join(block(f"{line}\n".encode()) for line in lines)
    path.write_bytes(header + body + block(b""))
    index = b"TBI\x01" + struct.pack("<8i", 1, 0, 1, 2, 2, ord("#"), 0, 2) + b"3\0"
    index += struct.pack("<iIi2Qi", 1, 0, 1, len(header) << 16, (len(header) + len(body)) << 16, 0)
    Path(f"{path}.tbi").write_bytes(gzip.compress(index))
    return path


def spliceai_cache(directory: Path) -> Path:
    """A made-up candidate cache on chromosome 3: ENST1 a transcript of GENEA, ENST2 of GENEB
    and ENST3 of a gene no score entry names."""
    shards = directory / "transcripts" / "3"
    shards.mkdir(parents=True)
    genes = (("ENST1", "GENEA"), ("ENST2", "GENEB"), ("ENST3", "GENEX"))
    (shards / "7000001-8000000.json").write_text(json.dumps([{"stable_id": feature, "gene_symbol": gene} for feature, gene in genes]), encoding="utf-8")
    return directory


def spliceai_row(feature: str, value: str, flag: str = "", base: bool = False) -> str:
    """A row of the deletion on ``feature`` with ``value`` as its SpliceAI_pred and ``flag`` as
    its cutoff flag, or without either; under ``base`` a row of the substitution
    ``SPLICEAI_BASE``, which stands at the one position with the alternate base as its allele."""
    extra = "IMPACT=MODIFIER" + (f";SpliceAI_pred={value}" if value else "") + (f";SpliceAI_cutoff={flag}" if flag else "")
    if base:
        return row(f"3:{SPLICEAI_BASE[0]}", SPLICEAI_BASE[2], feature, extra, uploaded="snv_two_genes")
    return row("3:7000200-7000201", "-", feature, extra, uploaded="del_two_genes")


def plant_spliceai(monkeypatch) -> None:
    plant_registry(monkeypatch, {"SpliceAI": ["SpliceAI_pred", "SpliceAI_cutoff"]})
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)


def test_spliceai_table_is_handed_to_the_classes_of_spliceai_which_set_nothing_aside_without_it(tmp_path, monkeypatch):
    """A made-up score record of two entries, GENEA's then GENEB's, and the rows of a
    transcript of each gene: the reference writes the record's whole value on the first
    gene's row and nothing on the second's, the engine each gene's own entry and a flag on
    both. With ``--spliceai-table`` and the candidate cache, which names each row's gene, the
    comparator's own predicate, unplanted, sets the three tuples aside that the record
    explains; without the files, without the cache, and under a file that scores each gene in
    a record of its own, nothing is set aside. The report names the files by path alone."""
    plant_spliceai(monkeypatch)
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    reference = output(tmp_path / "reference.txt", spliceai_row("ENST1", SPLICEAI_WHOLE, "FAIL"), spliceai_row("ENST2", ""))
    engine = output(tmp_path / "engine.txt", spliceai_row("ENST1", first, "FAIL"), spliceai_row("ENST2", second, "PASS"))
    without = score(tmp_path, reference, {"engine": engine}, table=SPLICEAI_CLASSES_TABLE, plugins=["SpliceAI"])
    assert without["spliceai_table"] is None
    item = entry(without, "SpliceAI")
    assert (counts(item), item["residual_rows"], item["masked"]) == ((2, 4, 1, 2, 4, 1), 4, {SPLICEAI_MULTI: 0})
    joined = spliceai_score_file(tmp_path / "joined")
    cache = spliceai_cache(tmp_path / "cache")
    unnamed = score(tmp_path, reference, {"engine": engine}, table=SPLICEAI_CLASSES_TABLE, plugins=["SpliceAI"], spliceai=[joined])
    assert entry(unnamed, "SpliceAI")["masked"] == {SPLICEAI_MULTI: 0} and entry(unnamed, "SpliceAI")["residual_rows"] == 4
    report = score(tmp_path, reference, {"engine": engine}, table=SPLICEAI_CLASSES_TABLE, plugins=["SpliceAI"], spliceai=[joined], cache=cache)
    assert report["spliceai_table"] == [{"path": str(joined)}]
    item = entry(report, "SpliceAI")
    assert (counts(item), item["residual_rows"], item["masked"]) == ((2, 4, 1, 1, 1, 1), 0, {SPLICEAI_MULTI: 3})
    item = entry(report, "SpliceAI", projection="allele")
    assert (item["adj_perl"], item["adj_rust"], item["adj_intersection"], item["residual_rows"]) == (1, 1, 1, 0)
    single = spliceai_score_file(tmp_path / "single", join=False)
    apart = score(tmp_path, reference, {"engine": engine}, table=SPLICEAI_CLASSES_TABLE, plugins=["SpliceAI"], spliceai=[single], cache=cache)
    assert entry(apart, "SpliceAI")["masked"] == {SPLICEAI_MULTI: 0} and entry(apart, "SpliceAI")["residual_rows"] == 4


def test_spliceai_rows_off_the_shape_of_a_record_of_several_entries_stay_residual(tmp_path, monkeypatch):
    """Under the score file of the record of two entries and the cache that names each row's
    gene: an engine value that is the other gene's entry, or the row's own entry with one
    score changed, a reference that writes no value on the first entry's gene's row, and an
    entry on the row of a gene no entry names, each leave their tuples residual."""
    plant_spliceai(monkeypatch)
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    joined = spliceai_score_file(tmp_path / "joined")
    cache = spliceai_cache(tmp_path / "cache")
    whole_row, silent_row = spliceai_row("ENST1", SPLICEAI_WHOLE), spliceai_row("ENST2", "")
    for name, reference_rows, engine_rows, masked, residual in (
        ("other-entry", (whole_row, silent_row), (spliceai_row("ENST1", first), spliceai_row("ENST2", first)), 1, 1),
        ("changed", (whole_row, silent_row), (spliceai_row("ENST1", first), spliceai_row("ENST2", second.replace("0.02", "0.03"))), 1, 1),
        ("silent-first", (spliceai_row("ENST1", ""), silent_row), (spliceai_row("ENST1", first), spliceai_row("ENST2", second)), 1, 1),
        ("unnamed-gene", (whole_row, silent_row, spliceai_row("ENST3", "")),
         (spliceai_row("ENST1", first), spliceai_row("ENST2", second), spliceai_row("ENST3", second)), 2, 1),
    ):
        directory = tmp_path / name
        directory.mkdir()
        report = score(directory, output(directory / "reference.txt", *reference_rows), {"engine": output(directory / "engine.txt", *engine_rows)},
                       table=SPLICEAI_CLASSES_TABLE, plugins=["SpliceAI"], spliceai=[joined], cache=cache)
        item = entry(report, "SpliceAI")
        assert (item["masked"], item["residual_rows"]) == ({SPLICEAI_MULTI: masked}, residual), name


def test_spliceai_rows_of_a_single_base_record_of_several_entries_are_set_aside_and_stay_residual_off_its_shape(tmp_path, monkeypatch):
    """A made-up score record of a substitution, one base for one base, of two entries,
    GENEA's then GENEB's, and the rows of a transcript of each gene at the one position with
    the alternate base as their allele: the reference writes the record's whole value on the
    first gene's row and nothing on the second's, the engine each gene's own entry. Under the
    score file and the cache that names each row's gene the tuples of both rows are set
    aside; without the file, and under a file that scores each gene in a record of its own,
    none is. An engine value that is the other gene's entry, and a reference value that is
    the first entry joined with itself, leave their tuples residual."""
    plant_spliceai(monkeypatch)
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    whole = f"{first},{SPLICEAI_BASE[2]}|{second}"
    joined = spliceai_score_file(tmp_path / "joined", variant=SPLICEAI_BASE)
    single = spliceai_score_file(tmp_path / "single", join=False, variant=SPLICEAI_BASE)
    cache = spliceai_cache(tmp_path / "cache")
    for name, reference_values, engine_values, files, masked, residual in (
        ("several-entries", (whole, ""), (first, second), [joined], 2, 0),
        ("no-file", (whole, ""), (first, second), None, 0, 3),
        ("one-entry-records", (whole, ""), (first, second), [single], 0, 3),
        ("other-entry", (whole, ""), (first, first), [joined], 1, 1),
        ("first-twice", (f"{first},{SPLICEAI_BASE[2]}|{first}", ""), (first, second), [joined], 1, 2),
    ):
        directory = tmp_path / name
        directory.mkdir()
        reference, engine = (
            output(directory / f"{side}.txt", *(spliceai_row(feature, value, base=True) for feature, value in zip(("ENST1", "ENST2"), values)))
            for side, values in (("reference", reference_values), ("engine", engine_values))
        )
        report = score(directory, reference, {"engine": engine}, table=SPLICEAI_CLASSES_TABLE, plugins=["SpliceAI"], spliceai=files, cache=cache)
        item = entry(report, "SpliceAI")
        assert (item["masked"], item["residual_rows"]) == ({SPLICEAI_MULTI: masked}, residual), name


def test_spliceai_rows_are_set_aside_under_every_record_line_at_the_allele(tmp_path, monkeypatch):
    """A made-up score file whose record line of two entries, GENEA's then GENEB's, stands
    more than once at the substitution's allele, or beside a line of the same entries in the
    other order or of GENEA's entry alone: the reference writes a row the whole value of the
    last line its gene leads, or nothing where it leads none, the engine the gene's entry.
    Under the score file and the cache that names each row's gene every tuple the two fill
    differently is set aside; where GENEA's line of one entry stands last both write that
    entry on its row and the other row's tuple alone is set aside. Where the lines give a gene
    entries that differ in text, the engine writes the entry of the last line that names the
    gene, and its tuples are set aside under two and three lines of the same genes with other
    scores and under two lines of other scores that each gene leads one of. The tuple of a
    row whose gene no entry names, and the tuples of a row whose engine value is an earlier
    line's entry, stay residual, and without the file nothing is set aside."""
    plant_spliceai(monkeypatch)
    cache = spliceai_cache(tmp_path / "cache")
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    changed = "GENEA|0.31|0.00|0.44|0.09|3|-2|11|-30"
    other, third, last = f"GENEB|{SPLICEAI_SCORES['GENEA']}", f"GENEA|{SPLICEAI_SCORES['GENEB']}", "GENEB|0.31|0.00|0.44|0.09|3|-2|11|-30"
    both = (first, second)

    def whole(*values: str) -> str:
        return f",{SPLICEAI_BASE[2]}|".join(values)

    for name, lines, reference_values, engine_values, masked, residual in (
        ("one-line", [both], (whole(*both), ""), both, 2, 0),
        ("twice", [both, both], (whole(*both), ""), both, 2, 0),
        ("three-times", [both, both, both], (whole(*both), ""), both, 2, 0),
        ("other-order", [both, both[::-1]], (whole(*both), whole(*both[::-1])), both, 2, 0),
        ("one-entry-last", [both, (first,)], (first, ""), both, 1, 0),
        ("one-entry-first", [(first,), both], (whole(*both), ""), both, 2, 0),
        ("gene-unnamed", [both, both], (whole(*both), "", ""), (*both, second), 2, 1),
        ("other-scores", [both, (changed, other)], (whole(changed, other), ""), (changed, other), 2, 0),
        ("other-scores-three-times", [both, (changed, other), (third, last)], (whole(third, last), ""), (third, last), 2, 0),
        ("other-scores-other-order", [both, (other, changed)], (whole(*both), whole(other, changed)), (changed, other), 2, 0),
        ("earlier-line", [both, (changed, other)], (whole(changed, other), ""), (first, other), 1, 2),
    ):
        directory = tmp_path / name
        directory.mkdir()
        file = spliceai_score_file(directory / "scores", variant=SPLICEAI_BASE, records=lines)
        reference, engine = (
            output(directory / f"{side}.txt", *(spliceai_row(feature, value, base=True) for feature, value in zip(("ENST1", "ENST2", "ENST3"), values)))
            for side, values in (("reference", reference_values), ("engine", engine_values))
        )
        for files, expected in (([file], (masked, residual)), (None, (0, None))):
            report = score(directory, reference, {"engine": engine}, table=SPLICEAI_CLASSES_TABLE, plugins=["SpliceAI"], spliceai=files, cache=cache)
            item = entry(report, "SpliceAI")
            assert item["masked"] == {SPLICEAI_MULTI: expected[0]}, (name, files)
            assert item["residual_rows"] == expected[1] if files else item["residual_rows"] > 0, (name, files)


# The rows each ground of ``formats.SPLICEAI_GROUNDS`` refuses in ``spliceai_refused_rows``.
SPLICEAI_REFUSED = {
    "no_candidate_value": 1, "no_span": 1, "other_info": 1, "several_references": 1, "no_record": 2, "no_gene": 1, "gene_not_named": 1,
    "lone_entry": 1, "candidate_value": 1, "reference_not_silent": 1, "reference_value": 1, "reference_flag": 1,
}


def spliceai_refused_rows() -> tuple[list[str], list[str], list[str]]:
    """Made-up record lines on chromosome 3 and the reference's and the engine's rows over
    them: for each ground of the class one row that meets it first, two for the ground of no
    record line, then the two rows of a line of GENEA's entry and GENEB's, which the class
    takes. The rows under the two last grounds and under ``reference_not_silent`` differ in
    their value and in their flag."""
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    changed = "GENEA|0.31|0.00|0.44|0.09|3|-2|11|-30"
    whole = f"{first},C|{second}"

    def record(position: int, *entries: str, reference: str = "A", more: str = "") -> str:
        return "\t".join(["3", str(position), ".", reference, "C", ".", ".", "SpliceAI=" + ",".join(f"C|{entry}" for entry in entries) + more])

    def pairs(value: str, flag: str = "") -> str:
        return (f";SpliceAI_pred={value}" if value else "") + (f";SpliceAI_cutoff={flag}" if flag else "")

    cases = (
        ("3:7000410", "ENST1", [record(7000410, first, second)], pairs(whole), ""),
        ("3", "ENST1", [], pairs(whole), pairs(first)),
        ("3:7000430", "ENST1", [record(7000430, first, second, more=";DB")], pairs(whole), pairs(first)),
        ("3:7000440", "ENST1", [record(7000440, first, second), record(7000440, first, second, reference="G")], pairs(whole), pairs(first)),
        ("3:7000450", "ENST1", [], pairs(whole), pairs(first)),
        ("3:7000455", "ENST1", [], pairs(whole), pairs(first)),
        ("3:7000460", "ENST9", [record(7000460, first, second)], pairs(whole), pairs(first)),
        ("3:7000470", "ENST3", [record(7000470, first, second)], "", pairs(second)),
        ("3:7000490", "ENST1", [record(7000490, first)], pairs(first), pairs(changed)),
        ("3:7000500", "ENST1", [record(7000500, first, second), record(7000500, changed, second)], pairs(f"{changed},C|{second}"), pairs(first)),
        ("3:7000510", "ENST1", [record(7000510, second, first)], pairs(second, "FAIL"), pairs(first, "PASS")),
        ("3:7000520", "ENST1", [record(7000520, first, second)], pairs(second), pairs(first)),
        ("3:7000530", "ENST1", [record(7000530, first, second)], pairs(whole, "FAIL"), pairs(first, "PASS")),
        ("3:7000550", "ENST1", [record(7000550, first, second)], pairs(whole), pairs(first)),
        ("3:7000550", "ENST2", [], "", pairs(second)),
    )
    rows = [
        [row(location, "C", feature, "IMPACT=MODIFIER" + written[side], uploaded=f"snv_{at}") for at, (location, feature, _, *written) in enumerate(cases)]
        for side in (0, 1)
    ]
    return [line for _, _, lines, _, _ in cases for line in lines], rows[0], rows[1]


def test_summary_counts_the_rows_the_spliceai_class_refuses_by_ground_under_the_score_files(tmp_path, capsys, monkeypatch):
    """Under ``--spliceai-table`` the summary's last lines count the rows
    ``spliceai_multi_entry_cell`` was asked about and refused and, a line a ground, the refused
    rows under the first ground met for each: counts alone, with no row, gene or value among
    them, a row asked about for two fields counted once. Every ground refuses a made-up row
    here; the row under the ground of the candidate's value writes the entry of an earlier
    line than the last that names its gene. The report carries the same counts under ``spliceai_refusals``;
    without the score files nothing is counted or printed, and rows the class takes whole
    print the one line."""
    plant_spliceai(monkeypatch)
    records, reference_rows, engine_rows = spliceai_refused_rows()
    scores = spliceai_score_file(tmp_path / "scores", variant=SPLICEAI_BASE, lines=records)
    classes = class_table(tmp_path, SPLICEAI_CLASSES_TABLE)
    cache = spliceai_cache(tmp_path / "cache")

    def run(name: str, reference: list[str], engine: list[str], *more: str) -> tuple[list[str], dict]:
        directory = tmp_path / name
        directory.mkdir()
        code, lines, error = run_main(
            capsys,
            "--reference", str(output(directory / "reference.txt", *reference)), "--engine", f"engine={output(directory / 'engine.txt', *engine)}",
            "--plugins", "SpliceAI", "--json-cache", str(cache), "--class-table", str(classes), "--work-dir", str(directory / "work"),
            "--sort-memory", "16M", "--report", str(directory / "report.json"), *more,
        )  # fmt: skip
        assert (code, error) == (0, "")
        return lines, json.loads((directory / "report.json").read_text(encoding="utf-8"))

    assert list(SPLICEAI_REFUSED) == list(formats.SPLICEAI_GROUNDS)
    lines, report = run("with", reference_rows, engine_rows, "--spliceai-table", str(scores))
    counted = [
        f"{SPLICEAI_MULTI} under --spliceai-table, engine engine: rows asked about 15, rows refused 13",
        *(f"{SPLICEAI_MULTI}, engine engine: rows refused {rows}: {formats.SPLICEAI_GROUNDS[ground]}" for ground, rows in SPLICEAI_REFUSED.items()),
    ]
    assert lines[-len(counted):] == counted
    assert counted[1].endswith(": rows refused 1: the candidate's row writes no SpliceAI value")
    assert counted[5].endswith(": rows refused 2: no record line of the score files stands at the row's variant")
    assert counted[7].endswith(": rows refused 1: no entry of a record line at the row's variant names the row's gene")
    assert counted[9].endswith(": rows refused 1: the candidate's value is not the gene's entry of the last record line that names it")
    assert not any(word in line for line in counted for word in ("GENE", "ENST", "7000", "0.1", "snv_"))
    assert f"SpliceAI\tengine\tfeature\tclass {SPLICEAI_MULTI} set aside 2 [met, -]" in lines[:-len(counted)]
    assert not any("rows refused" in line or "rows asked about" in line for line in lines[:-len(counted)])
    assert report["spliceai_refusals"] == [{"engine": "engine", "asked": 15, "refused": 13, "grounds": SPLICEAI_REFUSED}]
    lines, report = run("without", reference_rows, engine_rows)
    assert lines[-1] == SPLICEAI_NOTICE and report["spliceai_refusals"] is None
    assert not any("rows refused" in line or "rows asked about" in line for line in lines)
    lines, report = run("taken", reference_rows[-2:], engine_rows[-2:], "--spliceai-table", str(scores))
    assert lines[-1] == f"{SPLICEAI_MULTI} under --spliceai-table, engine engine: rows asked about 2, rows refused 0"
    assert report["spliceai_refusals"] == [{"engine": "engine", "asked": 2, "refused": 0, "grounds": {}}]
    lines, report = run("equal", reference_rows, reference_rows, "--spliceai-table", str(scores))
    assert report["spliceai_refusals"] == [] and not any("rows asked about" in line for line in lines)


def test_summary_says_that_the_spliceai_class_sets_nothing_aside_without_the_score_files(tmp_path, capsys, monkeypatch):
    plant_spliceai(monkeypatch)
    joined = spliceai_score_file(tmp_path / "joined")
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    reference = output(tmp_path / "reference.txt", spliceai_row("ENST1", SPLICEAI_WHOLE), spliceai_row("ENST2", ""))
    engine = output(tmp_path / "engine.txt", spliceai_row("ENST1", first), spliceai_row("ENST2", second))
    classes = class_table(tmp_path, SPLICEAI_CLASSES_TABLE)

    cache = spliceai_cache(tmp_path / "cache")

    def run(name, *more):
        return run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={engine}", "--plugins", "SpliceAI", "--json-cache", str(cache),
            "--class-table", str(classes), "--work-dir", str(tmp_path / f"work-{name}"), "--sort-memory", "16M",
            "--report", str(tmp_path / "report.json"), "--residuals", str(tmp_path / name), *more,
        )  # fmt: skip

    code, lines, error = run("without")
    assert (code, error) == (0, "")
    assert lines[-1] == SPLICEAI_NOTICE
    assert f"SpliceAI\tengine\tfeature\tclass {SPLICEAI_MULTI} set aside 0 (none) [met, -]" in lines
    code, lines, error = run("with", "--spliceai-table", str(joined), "--spliceai-table", str(spliceai_score_file(tmp_path / "single", join=False)))
    assert (code, error) == (0, "")
    assert SPLICEAI_NOTICE not in lines
    assert f"SpliceAI\tengine\tfeature\tclass {SPLICEAI_MULTI} set aside 2 [met, -]" not in lines
    code, lines, error = run("joined", "--spliceai-table", str(joined))
    assert (code, error) == (0, "")
    assert f"SpliceAI\tengine\tfeature\tclass {SPLICEAI_MULTI} set aside 2 [met, -]" in lines
    assert json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))["spliceai_table"] == [{"path": str(joined)}]
    code, lines, error = run("missing", "--spliceai-table", str(tmp_path / "nowhere.synthetic.vcf.gz"))
    assert (code, lines) == (2, [])
    assert error == f"ERROR: [score_plugin_tuples] not a file: {tmp_path / 'nowhere.synthetic.vcf.gz'}\n"
    Path(f"{joined}.tbi").unlink()
    code, lines, error = run("unindexed", "--spliceai-table", str(joined))
    assert (code, lines) == (2, [])
    assert error.startswith(f"ERROR: [score_plugin_tuples] --spliceai-table: cannot read the tabix index {joined}.tbi")


GOLDEN_CORPORA = HERE.parent.parent / "tests/golden/116"
# Per SpliceAI corpus and set (``base`` the manifest's own run): the keys under which the
# class sets a tuple aside in the feature projection and in the allele projection. Each
# corpus's deletion and substitution two genes share are each scored by one record of two
# entries, a further substitution in their span by a record line of two entries that stands
# twice, another by a line of three entries that stands three times, and another by two
# lines of two entries whose scores differ; the rows of the first entry's gene differ in
# one field, those of the second entry's in every field the run writes, and the allele
# projection counts a field once.
SPLICEAI_CORPUS_COUNTS = {
    ("GRCh37-plugin-spliceai", "base"): (20, 5),
    ("GRCh37-plugin-spliceai", "cutoff"): (30, 7),
    ("GRCh37-plugin-spliceai", "split-output"): (100, 45),
    ("GRCh37-plugin-spliceai", "split-output-cutoff"): (110, 47),
    ("GRCh38-plugin-spliceai", "base"): (20, 5),
    ("GRCh38-plugin-spliceai", "cutoff"): (25, 8),
    ("GRCh38-plugin-spliceai", "split-output"): (60, 45),
    ("GRCh38-plugin-spliceai", "split-output-cutoff"): (65, 48),
}


def corpus_engine_output(target: dict, expected: Path, path: Path) -> Path:
    """vep-rs's default-format output of a corpus run as the corpus's committed files state
    it: Ensembl VEP's expected rows with each field difference the run's manifest entry
    documents under a compared field written as vep-rs's value."""
    documented: dict[tuple[str, str, str], dict[str, str]] = defaultdict(dict)
    for item in target["field_divergences"]:
        if item["field"] in target["fields_compared"]:
            documented[item["location"], item["allele"], item["feature"]][item["field"]] = item["vep_rs"]
    lines = []
    with gzip.open(expected, "rt", encoding="utf-8") as handle:
        for line in handle:
            cells = line.rstrip("\n").split("\t")
            fields = None if line.startswith("#") else documented.get((cells[1], cells[2], cells[4]))
            if fields:
                pairs = dict(pair.partition("=")[::2] for pair in cells[13].split(";")) if cells[13] != "-" else {}
                cells[13] = ";".join(f"{name}={value}" for name, value in {**pairs, **fields}.items())
                line = "\t".join(cells) + "\n"
            lines.append(line)
    path.write_text("".join(lines), encoding="utf-8")
    return path


def test_summary_says_that_the_spliceai_class_names_no_gene_from_a_cache_without_one(tmp_path, capsys, monkeypatch):
    plant_spliceai(monkeypatch)
    joined = spliceai_score_file(tmp_path / "joined")
    first, second = (f"{symbol}|{scores}" for symbol, scores in SPLICEAI_SCORES.items())
    reference = output(tmp_path / "reference.txt", spliceai_row("ENST1", SPLICEAI_WHOLE), spliceai_row("ENST2", ""))
    engine = output(tmp_path / "engine.txt", spliceai_row("ENST1", first), spliceai_row("ENST2", second))
    classes = class_table(tmp_path, SPLICEAI_CLASSES_TABLE)
    cache = spliceai_cache(tmp_path / "cache")

    def run(name, *more):
        return run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={engine}", "--plugins", "SpliceAI",
            "--class-table", str(classes), "--work-dir", str(tmp_path / f"work-{name}"), "--sort-memory", "16M",
            "--report", str(tmp_path / "report.json"), "--residuals", str(tmp_path / name), *more,
        )  # fmt: skip

    code, lines, error = run("table-alone", "--spliceai-table", str(joined))
    assert (code, error) == (0, "")
    # The notice stands before the counts of the rows the class refused, both for want of a named gene.
    assert lines[-3:] == [
        SPLICEAI_CACHE_NOTICE,
        f"{SPLICEAI_MULTI} under --spliceai-table, engine engine: rows asked about 2, rows refused 2",
        f"{SPLICEAI_MULTI}, engine engine: rows refused 2: {formats.SPLICEAI_GROUNDS['no_gene']}",
    ]
    assert f"SpliceAI\tengine\tfeature\tclass {SPLICEAI_MULTI} set aside 0 (none) [met, -]" in lines
    code, lines, error = run("both", "--spliceai-table", str(joined), "--json-cache", str(cache))
    assert (code, error) == (0, "")
    assert SPLICEAI_CACHE_NOTICE not in lines
    assert f"SpliceAI\tengine\tfeature\tclass {SPLICEAI_MULTI} set aside 2 [met, -]" in lines
    code, lines, error = run("neither")
    assert (code, error) == (0, "")
    assert lines[-1] == SPLICEAI_NOTICE
    assert SPLICEAI_CACHE_NOTICE not in lines


def test_the_spliceai_corpora_apply_the_class_of_a_record_of_several_entries_at_their_count(tmp_path, capsys, monkeypatch):
    """Each run of the two SpliceAI corpora, scored as the plugin gate scores it: Ensembl
    VEP's committed output against vep-rs's as the manifest documents it, under the SpliceAI
    rows of the committed class table and the corpus's own score files. The class sets aside
    exactly the keys of ``SPLICEAI_CORPUS_COUNTS`` and nothing is left, and under a table whose
    row reads ``not met`` the same tuples are residual."""
    members = [f"SpliceAI_pred_{member}" for member in formats.SPLICEAI_MEMBERS]
    plant_registry(monkeypatch, {"SpliceAI": ["SpliceAI_pred", "SpliceAI_cutoff", *members]})
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    rows = [line for line in COMMITTED_CLASS_TABLE.read_text(encoding="utf-8").splitlines(keepends=True) if line.split("\t")[1:2] == ["SpliceAI"]]
    met = class_table(tmp_path, "".join(rows))
    (tmp_path / "unmet").mkdir()
    unmet = class_table(tmp_path / "unmet", "".join(line.replace("\tmet\t", "\tnot met\t") if line.startswith(SPLICEAI_MULTI) else line for line in rows))

    def scored(reference: Path, engine: Path, directory: Path, tables: list[str], work: Path, table: Path, label: str) -> dict:
        code, _, error = run_main(
            capsys, "--reference", str(reference), "--engine", f"vep-rs={engine}", "--plugins", "SpliceAI", "--class-table", str(table),
            "--json-cache", str(directory / "json_cache"), *tables, "--work-dir", str(work / f"work-{label}"), "--sort-memory", "16M",
            "--report", str(work / f"{label}.json"),
        )  # fmt: skip
        assert (code, error) == (0, ""), (work, error)
        return json.loads((work / f"{label}.json").read_text(encoding="utf-8"))

    for (corpus, name), (feature, allele) in SPLICEAI_CORPUS_COUNTS.items():
        directory = GOLDEN_CORPORA / corpus
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        target = manifest if name == "base" else next(entry for entry in manifest["flag_sets"] if entry["name"] == name)
        reference = directory / ("expected" if name == "base" else target["expected_dir"]) / "default.txt.gz"
        work = tmp_path / corpus / name
        work.mkdir(parents=True)
        engine = corpus_engine_output(target, reference, work / "vep-rs.txt")
        tables = [word for kind in ("snv", "indel") for file in sorted((directory / "plugin_data").glob(f"spliceai_scores.raw.{kind}.*.vcf.gz"))
                  for word in ("--spliceai-table", str(file))]
        applied = scored(reference, engine, directory, tables, work, met, "met")
        left = scored(reference, engine, directory, tables, work, unmet, "unmet")
        for projection, count in (("feature", feature), ("allele", allele)):
            item = entry(applied, "SpliceAI", engine="vep-rs", projection=projection)
            assert item["masked"][SPLICEAI_MULTI] == count, (corpus, name, projection)
            assert item["adj_perl"] == item["adj_rust"] == item["adj_intersection"] and item["residual_rows"] == 0, (corpus, name, projection)
            item = entry(left, "SpliceAI", engine="vep-rs", projection=projection)
            assert item["not_applied"] == {SPLICEAI_MULTI: count} and item["residual_rows"] > 0, (corpus, name, projection)


def test_pairing_tool_hands_the_rows_of_a_splice_term_past_the_record_span_to_its_plugin_value_class(tmp_path, monkeypatch, pairing_tool):
    """``plugin_value_follows_splice_term_past_the_record_span`` reads both rows of a key: their
    first column, Location, Allele, Consequence and cDNA and CDS positions. The pairing tool
    hands them over with the remainder, and the record of a row under an identifier is read
    from the input on either path, so both paths set the same values aside on the rows of
    ``splice_term_past_the_record_span``, with the record from the first column or from the
    input and without it, leave the same values on the rows it refuses and write the same
    report and residual files, with the input and without it."""
    klass = formats.PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN
    plant_registry(monkeypatch, SPLICE_SPAN_PLUGINS)
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)
    reference, engine, run_input = splice_span_outputs(tmp_path)
    for given, scores in ((run_input, SPLICE_SPAN_SCORES), (None, SPLICE_SPAN_SCORES_WITHOUT_INPUT)):
        report, files = on_both_paths(
            tmp_path, pairing_tool, reference, {"engine": engine}, table=SPLICE_SPAN_TABLE, plugins=list(SPLICE_SPAN_PLUGINS),
            run_input=given,
        )
        for source, (counted, keys, left) in scores.items():
            item = entry(report, source)
            assert (counts(item), item["masked"]) == (counted, {klass: keys}), (given, source)
            assert files[f"{source}.engine.feature.tsv"].decode().splitlines() == left, (given, source)


SKIPPED_RECORD_CLASS = "plugin_value_on_reference_skipped_record"
SKIPPED_RECORD_PLUGINS = {"DosageSensitivity": ["pHaplo", "pTriplo"], "UTRAnnotator": ["Existing_uORFs", "5UTR_consequence"]}
SKIPPED_RECORD_TABLE = f"{SKIPPED_RECORD_CLASS}\tDosageSensitivity;UTRAnnotator\tpHaplo,pTriplo;Existing_uORFs\tmet\t-\n"
# The fabricated grid of the class: per cell its name, the reference's rows, the engine's rows,
# and whether the engine's values there are set aside with the input records and without them.
SKIPPED_RECORD_GRID = (
    ("a structural record over the size limit", [],
     [row("21:7000001-27000000", "deletion", "ENST1", "pHaplo=0.9;pTriplo=0.2", uploaded="sv_over"),
      row("21:7000001-27000000", "deletion", "ENST2", "pHaplo=0.5;Existing_uORFs=2", uploaded="sv_over")], True, True),
    ("a structural record under the limit", [], [row("21:30000001-30000900", "duplication", "ENST1", "pHaplo=0.7", uploaded="sv_under")], True, True),
    ("a sequence allele beside the symbolic one", [],
     [row("21:37000001-37000002", "deletion", "ENST1", "pHaplo=0.11", uploaded="sv_mixed"),
      row("21:37000001-37000002", "T", "ENST1", "pHaplo=0.12", uploaded="sv_mixed")], True, True),
    ("every row a sequence allele", [], [row("21:31000000", "T", "ENST1", "pHaplo=0.6", uploaded="snv_only")], False, False),
    ("the allele tandem_repeat", [], [row("21:32000001-32000050", "tandem_repeat", "ENST1", "pHaplo=0.4", uploaded="tr_1")], False, False),
    ("the reference wrote a row at the Location", [row("21:33000001-33000500", "deletion", "ENST1", "IMPACT=MODIFIER", uploaded="sv_kept")],
     [row("21:33000001-33000500", "deletion", "ENST1", "pHaplo=0.3", uploaded="sv_kept")], False, False),
    ("a reference row carries the record's name elsewhere", [row("21:34000001-34000500", "deletion", "ENST1", "IMPACT=MODIFIER", uploaded="sv_moved")],
     [row("21:34000000-34000500", "deletion", "ENST1", "pHaplo=0.8", uploaded="sv_moved")], False, False),
    ("an engine row at the Location names no record", [],
     [row("21:35000001-35900000", "deletion", "ENST1", "pHaplo=0.2", uploaded="sv_named"),
      row("21:35000001-35900000", "deletion", "ENST2", "pHaplo=0.25", uploaded=".")], False, False),
    ("a field the class does not read", [], [row("21:36000001-36900000", "deletion", "ENST1", "5UTR_consequence=gained", uploaded="sv_field")], False, False),
    ("a gVCF variant site under an identifier", [],
     [row("21:38000000", "T", "ENST1", "pHaplo=0.13", uploaded="gv_1"), row("21:38000000", "<NON_REF>", "ENST1", "pHaplo=0.14", uploaded="gv_1")], False, True),
    ("a gVCF variant site under a generated name", [],
     [row("21:39000000", "T", "ENST1", "pHaplo=0.15", uploaded="21_39000000_A/T/<NON_REF>"),
      row("21:39000000", "<NON_REF>", "ENST1", "pHaplo=0.16", uploaded="21_39000000_A/T/<NON_REF>")], False, False),
)  # fmt: skip
SKIPPED_RECORD_INPUT = (
    ("21", "7000000", "sv_over", "N", "<DEL>"), ("21", "30000000", "sv_under", "N", "<DUP>"), ("21", "37000000", "sv_mixed", "N", "<DEL>,NT"),
    ("21", "38000000", "gv_1", "A", "T,<NON_REF>"), ("21", "39000000", ".", "A", "T,<NON_REF>"),
)  # fmt: skip


def plant_skipped_record_class(monkeypatch) -> None:
    """Plants the two plugins' fields and the class's name, with a predicate that takes no
    pair: the scorer decides the class without one."""
    plant_registry(monkeypatch, SKIPPED_RECORD_PLUGINS)
    monkeypatch.setattr(formats, "FIELD_CLASSES", {SKIPPED_RECORD_CLASS: ()})
    monkeypatch.setattr(formats, "field_class_accepts", lambda *rows, **more: False)


def skipped_record_outputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    kept = row("21:100", "T", "ENST9", "pHaplo=0.01", uploaded="snv_1")
    reference = output(tmp_path / "reference.txt", kept, row("21:40000001-40000500", "deletion", "ENST1", "pHaplo=0.1", uploaded="sv_ref"),
                       *(line for _, rows, _, _, _ in SKIPPED_RECORD_GRID for line in rows))
    engine = output(tmp_path / "engine.txt", kept, *(line for _, _, rows, _, _ in SKIPPED_RECORD_GRID for line in rows))
    return reference, engine, vcf(tmp_path / "input.vcf", *SKIPPED_RECORD_INPUT)


def skipped_record_cells_set_aside(residuals: Path) -> list[bool]:
    """Per cell of the grid, whether no value of its engine rows is left in the residual files."""
    left = {tuple(line.split("\t")[1:4]) for source in SKIPPED_RECORD_PLUGINS for line in residual_lines(residuals, source)}
    return [not any(tuple(line.split("\t")[i] for i in (1, 2, 4)) in left for line in rows) for _, _, rows, _, _ in SKIPPED_RECORD_GRID]


@pytest.mark.parametrize("with_input", [True, False])
def test_plugin_value_on_a_record_the_reference_dropped_is_set_aside_on_the_grid_and_nothing_wider(tmp_path, monkeypatch, with_input):
    """A value a default-format engine alone writes at a Location of a record the reference
    dropped is the class's: the engine wrote a structural row there and the reference none,
    and no engine row there names no record or one a reference row carries. Every row at such
    a Location goes with it, a sequence allele beside the symbolic one included, whatever the
    record's span. A Location of sequence alleles alone, ``tandem_repeat``, a Location the
    reference wrote a row at, a field the class does not read and a gVCF variant site, which
    the input's lines show under an identifier and a generated name shows without them, stay
    residuals, and so does a value the reference alone writes."""
    plant_skipped_record_class(monkeypatch)
    reference, engine, run_input = skipped_record_outputs(tmp_path)
    residuals = tmp_path / "residuals"
    report = score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=SKIPPED_RECORD_TABLE,
                   plugins=list(SKIPPED_RECORD_PLUGINS), run_input=run_input if with_input else None)
    expected = [with_records if with_input else without for _, _, _, with_records, without in SKIPPED_RECORD_GRID]
    assert skipped_record_cells_set_aside(residuals) == expected
    aside = 7 if with_input else 9
    assert entry(report, "DosageSensitivity")["masked"] == {SKIPPED_RECORD_CLASS: aside - 1}
    assert entry(report, "UTRAnnotator")["masked"] == {SKIPPED_RECORD_CLASS: 1}
    assert "reference\t21:40000001-40000500\tdeletion\tENST1\tpHaplo\t0.1" in residual_lines(residuals, "DosageSensitivity")


@pytest.mark.parametrize("with_input", [True, False])
def test_plugin_value_class_takes_the_rows_the_key_comparison_sets_aside_as_a_skipped_record(tmp_path, monkeypatch, with_input):
    """The class reads a dropped record as ``compare_vep_outputs.py`` does: over the grid's
    two outputs that program sets aside under ``reference_skipped_record`` exactly the engine
    rows whose values the scorer sets aside, with the input records and without them."""
    import compare_vep_outputs

    plant_skipped_record_class(monkeypatch)
    reference, engine, run_input = skipped_record_outputs(tmp_path)
    residuals = tmp_path / "residuals"
    score(tmp_path, reference, {"engine": engine}, residuals=residuals, table=SKIPPED_RECORD_TABLE,
          plugins=list(SKIPPED_RECORD_PLUGINS), run_input=run_input if with_input else None)
    scored = skipped_record_cells_set_aside(residuals)
    for side, written in (("perl", reference), ("rust", engine)):
        (tmp_path / side).mkdir()
        shutil.copy(written, tmp_path / side / "output.txt")
    words = ["compare_vep_outputs.py", "--perl-dir", str(tmp_path / "perl"), "--rust-dir", str(tmp_path / "rust"),
             "--report-dir", str(tmp_path / "report"), "--assembly", "GRCh38", "--glob", "*.txt"]
    monkeypatch.setattr(sys, "argv", words + (["--input-vcf", str(run_input)] if with_input else []))
    # main() lifts the interpreter's bound on an integer's digits for the whole process, so
    # the test puts the bound back.
    digits = sys.get_int_max_str_digits()
    try:
        assert compare_vep_outputs.main() == 0
    finally:
        sys.set_int_max_str_digits(digits)
    adjusted = json.loads((tmp_path / "report" / "summary.json").read_text(encoding="utf-8"))["aggregate"]["adjusted"]
    one_sided, left_open = (
        {tuple(line.split("\t")[2:5]) for line in (tmp_path / "report" / name).read_text(encoding="utf-8").splitlines()[1:]}
        for name in ("discordant.tsv", "discordant_open.tsv")
    )
    rows_of = [{tuple(line.split("\t")[i] for i in (1, 2, 4)) for line in rows} for _, _, rows, _, _ in SKIPPED_RECORD_GRID]
    # A gVCF variant site its generated name or the input shows is set aside too, under its own class.
    other_class = rows_of[10] | (rows_of[9] if with_input else set())
    compared = [bool(rows) and rows <= one_sided - left_open - other_class for rows in rows_of]
    compared[8] = False  # the row of the field the class does not read is a dropped record's row all the same
    assert scored == compared
    assert adjusted["excluded_categories"]["reference_skipped_record"]["excluded"] == sum(len(rows) for rows, aside in zip(rows_of, compared) if aside) + 1
    assert adjusted["excluded_categories"]["gvcf_variant_site_dropped"]["excluded"] == len(other_class)


def test_plugin_value_on_a_dropped_record_is_left_for_a_tuple_file_engine_and_read_alike_on_both_paths(tmp_path, monkeypatch, pairing_tool):
    """A tuple file states no rows, so none of its values is at a dropped record's Location;
    and the script's own path and the pairing tool write the same report and residual files
    over the grid, the reference stored plain or gzip."""
    plant_skipped_record_class(monkeypatch)
    reference, engine, run_input = skipped_record_outputs(tmp_path)
    listed = tuples(tmp_path / "engine.tsv", ("21:7000001-27000000", "deletion", "ENST1", "pHaplo", "0.9"))
    item = entry(score(tmp_path, reference, {"engine": listed}, table=SKIPPED_RECORD_TABLE, plugins=list(SKIPPED_RECORD_PLUGINS)), "DosageSensitivity")
    assert item["masked"] == {SKIPPED_RECORD_CLASS: 0}
    stored = tmp_path / "reference.txt.gz"
    stored.write_bytes(gzip.compress(reference.read_bytes(), mtime=0))
    reports = []
    for given in (reference, stored):
        report, files = on_both_paths(tmp_path, pairing_tool, given, {"engine": engine}, table=SKIPPED_RECORD_TABLE,
                                      plugins=list(SKIPPED_RECORD_PLUGINS), run_input=run_input)
        reports.append((entry(report, "DosageSensitivity")["masked"], files))
    assert reports[0] == reports[1] and reports[0][0] == {SKIPPED_RECORD_CLASS: 6}


def test_dropped_record_reads_the_reference_for_the_names_at_a_location_the_scope_took(tmp_path, monkeypatch):
    """The reference's first column is read only for a Location at which the engine wrote a
    structural row and the reference none, and a name it carries on a row of a gzip output
    gives the Location back."""
    dropped = [row("21:7000001-27000000", "deletion", "ENST1", uploaded="sv_over")]
    stored = tmp_path / "reference.txt.gz"
    stored.write_bytes(gzip.compress((HEADER + row("21:5", "T", "ENST1", uploaded="sv_else")).encode(), mtime=0))
    read = []
    named = spt.reference_named_locations
    monkeypatch.setattr(spt, "reference_named_locations", lambda *given: read.append(given[2]) or named(*given))
    assert spt.dropped_by_the_reference([], dropped, stored, [])
    assert not spt.dropped_by_the_reference([], [row("21:5", "T", "ENST1", uploaded="snv")], stored, [])
    assert not spt.dropped_by_the_reference(dropped, dropped, stored, [])
    assert read == [{"21:7000001-27000000"}, set(), set()]
    stored.write_bytes(gzip.compress((HEADER + row("21:9000-9500", "deletion", "ENST1", uploaded="sv_over")).encode(), mtime=0))
    assert not spt.dropped_by_the_reference([], dropped, stored, [])


MAXENTSCAN_CLASS = "maxentscan_overlapping_windows"
MAXENTSCAN_TABLE = f"{MAXENTSCAN_CLASS}\tMaxEntScan\tMaxEntScan_*\tmet\t-\t-\tsets nothing aside without --json-cache\n"
# A donor window's nine bases and an acceptor window's twenty-three.
MAXENTSCAN_DONOR_TEXT, MAXENTSCAN_ACCEPTOR_TEXT = "AC" * 4 + "G", "AC" * 11 + "G"


def unplant_maxentscan(monkeypatch) -> None:
    """The registry's MaxEntScan fields under the comparator's own classes and predicate."""
    plant_registry(monkeypatch, {"MaxEntScan": ["MaxEntScan_ref", "MaxEntScan_alt", "MaxEntScan_diff", "MaxEntScan_ref_seq", "MaxEntScan_alt_seq"]})
    monkeypatch.setattr(formats, "FIELD_CLASSES", COMPARATOR_FIELD_CLASSES)
    monkeypatch.setattr(formats, "field_class_accepts", COMPARATOR_FIELD_CLASS_ACCEPTS)


def maxentscan_inputs(tmp_path: Path, verbose: bool = False) -> tuple[Path, Path, Path]:
    """A reference, an engine and a candidate cache of chromosome 1 whose ENST1 (forward) and
    ENST2 (reverse) have an exon of three bases, 201 to 203, between an intron from 101 to 200
    and one from 204 to 300, so that 201 to 203 lie in the acceptor window of one intron and
    the donor window of the next on both, 199 in ENST1's acceptor window alone and 150 in no
    window. The two sides write the scores of two different windows, and under ``verbose``
    the window texts, on both transcripts at 202 and on ENST1 at 199, and agree at 150."""
    exons = [{"start": 1, "end": 100}, {"start": 201, "end": 203}, {"start": 301, "end": 400}]
    transcripts = [{"stable_id": "ENST1", "strand": 1, "exons": exons}, {"stable_id": "ENST2", "strand": -1, "exons": exons}]
    cache = tmp_path / "cache"
    (cache / "transcripts" / "1").mkdir(parents=True, exist_ok=True)
    (cache / "transcripts" / "1" / "1-1000000.json").write_text(json.dumps(transcripts), encoding="utf-8")

    def extra(scores: tuple[str, str, str], window: str) -> str:
        fields = dict(zip(("MaxEntScan_ref", "MaxEntScan_alt", "MaxEntScan_diff"), scores))
        if verbose:
            fields.update(MaxEntScan_ref_seq=window, MaxEntScan_alt_seq=window[:-1] + "T")
        return ";".join(f"{name}={value}" for name, value in fields.items())

    donor = extra(("8.100", "2.300", "5.800"), MAXENTSCAN_DONOR_TEXT)
    acceptor = extra(("6.400", "1.100", "5.300"), MAXENTSCAN_ACCEPTOR_TEXT)
    places = (("1:202", "ENST1"), ("1:202", "ENST2"), ("1:199", "ENST1"))
    agreed = row("1:150", "G", "ENST1", donor)
    reference = output(tmp_path / "reference.txt", *(row(location, "G", feature, donor) for location, feature in places), agreed)
    engine = output(tmp_path / "engine.txt", *(row(location, "G", feature, acceptor) for location, feature in places), agreed)
    return reference, engine, cache


@pytest.mark.parametrize("verbose", [False, True])
def test_maxentscan_overlapping_windows_is_applied_where_the_cache_places_the_record_in_two_windows(tmp_path, monkeypatch, verbose):
    """``maxentscan_overlapping_windows`` is decided by the comparator's own predicate,
    unplanted. With the candidate cache the class sets aside the fields of the rows whose
    record lies in two windows of their transcript, the three scores and, under verbose, the
    two window texts, counted once a key and field; the row of the same shape in one window
    stays a difference. Without the cache nothing is set aside."""
    unplant_maxentscan(monkeypatch)
    reference, engine, cache = maxentscan_inputs(tmp_path, verbose)
    fields = 5 if verbose else 3
    report = score(tmp_path, reference, {"engine": engine}, table=MAXENTSCAN_TABLE, plugins=["MaxEntScan"], cache=cache, residuals=tmp_path / "with")
    assert report["json_cache"] == {"path": str(cache)}
    feature, allele = entry(report, "MaxEntScan"), entry(report, "MaxEntScan", projection="allele")
    assert feature["masked"] == {MAXENTSCAN_CLASS: 2 * fields} and allele["masked"] == {MAXENTSCAN_CLASS: fields}
    assert counts(feature) == (4 * fields, 4 * fields, fields, 2 * fields, 2 * fields, fields)
    left = residual_lines(tmp_path / "with", "MaxEntScan")
    assert len(left) == 2 * fields and {line.split("\t")[1] for line in left} == {"1:199"}
    without = score(tmp_path, reference, {"engine": engine}, table=MAXENTSCAN_TABLE, plugins=["MaxEntScan"])
    assert without["json_cache"] is None
    assert entry(without, "MaxEntScan")["masked"] == {MAXENTSCAN_CLASS: 0} and entry(without, "MaxEntScan", projection="allele")["masked"] == {MAXENTSCAN_CLASS: 0}
    assert counts(entry(without, "MaxEntScan")) == (4 * fields, 4 * fields, fields, 4 * fields, 4 * fields, fields)


def test_summary_says_that_the_class_of_maxentscan_sets_nothing_aside_without_the_candidate_cache(tmp_path, capsys, monkeypatch):
    unplant_maxentscan(monkeypatch)
    reference, engine, cache = maxentscan_inputs(tmp_path, verbose=True)
    table = class_table(tmp_path, MAXENTSCAN_TABLE)
    notice = "no --json-cache: a class of MaxEntScan that reads the candidate cache sets nothing aside"

    def run(name, *more):
        return run_main(
            capsys,
            "--reference", str(reference), "--engine", f"engine={engine}", "--plugins", "MaxEntScan",
            "--class-table", str(table), "--work-dir", str(tmp_path / f"work-{name}"), "--sort-memory", "16M",
            "--report", str(tmp_path / "report.json"), *more,
        )  # fmt: skip

    code, lines, error = run("without")
    assert (code, error) == (0, "")
    assert lines[-1] == notice
    assert f"MaxEntScan\tengine\tfeature\tclass {MAXENTSCAN_CLASS} set aside 0 (none) [met, -, -, sets nothing aside without --json-cache]" in lines
    code, lines, error = run("with", "--json-cache", str(cache))
    assert (code, error) == (0, "")
    assert notice not in lines
    assert f"MaxEntScan\tengine\tfeature\tclass {MAXENTSCAN_CLASS} set aside 10 [met, -, -, sets nothing aside without --json-cache]" in lines
    assert f"MaxEntScan\tengine\tallele\tclass {MAXENTSCAN_CLASS} set aside 5 [met, -, -, sets nothing aside without --json-cache]" in lines


GOLDEN_CORPORA = HERE.parents[1] / "tests" / "golden" / "116"
# The keys under which the class sets a field aside on each run of the two MaxEntScan fixture
# corpora, in the feature projection and in the allele projection: GRCh38 carries one row in
# two windows and GRCh37 two rows of one record, each with three scores and, under verbose,
# two window texts.
MAXENTSCAN_CORPUS_COUNTS = {
    ("GRCh37-plugin-maxentscan-fixture", "base"): (6, 3),
    ("GRCh37-plugin-maxentscan-fixture", "swa-ncss"): (6, 3),
    ("GRCh37-plugin-maxentscan-fixture", "swa-ncss-verbose"): (10, 5),
    ("GRCh38-plugin-maxentscan-fixture", "base"): (3, 3),
    ("GRCh38-plugin-maxentscan-fixture", "swa-ncss"): (3, 3),
    ("GRCh38-plugin-maxentscan-fixture", "swa-ncss-verbose"): (5, 5),
}


def documented_candidate(corpus: Path, run: str, directory: Path, changed_score: bool = False) -> tuple[Path, Path, dict]:
    """A corpus run's committed Ensembl VEP default-format output unpacked, the candidate output
    the manifest documents (the same rows with vep-rs's value in place of Ensembl VEP's in every
    field the run's ``field_divergences`` name) and the run's manifest entry. With
    ``changed_score`` the candidate also writes another ``MaxEntScan_ref`` on the first row the
    manifest documents no field of."""
    manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
    target = manifest if run == "base" else next(item for item in manifest["flag_sets"] if item["name"] == run)
    expected = corpus / ("expected" if run == "base" else target["expected_dir"]) / "default.txt.gz"
    documented: dict[tuple[str, str, str], dict[str, tuple[str, str]]] = {}
    for item in target["field_divergences"]:
        documented.setdefault((item["location"], item["allele"], item["feature"]), {})[item["field"]] = (item["vep"], item["vep_rs"])
    lines = gzip.decompress(expected.read_bytes()).decode("utf-8").splitlines(keepends=True)
    written = []
    for line in lines:
        columns = line.rstrip("\n").split("\t")
        if not line.startswith("#"):
            fields = documented.get((columns[1], columns[2], columns[4]), {})
            pairs = [item.partition("=") for item in columns[13].split(";")]
            if not fields and changed_score and any(name == "MaxEntScan_ref" for name, _, _ in pairs):
                score = next(value for name, _, value in pairs if name == "MaxEntScan_ref")
                fields, changed_score = {"MaxEntScan_ref": (score, "99.999")}, False
            for name, (reference_value, candidate_value) in fields.items():
                assert (name, "=", reference_value) in pairs, (columns[1], name)
            columns[13] = ";".join(f"{name}={fields[name][1]}" if name in fields else f"{name}{sign}{value}" for name, sign, value in pairs)
        written.append("\t".join(columns) + "\n")
    directory.mkdir(parents=True, exist_ok=True)
    reference, candidate = directory / "reference.txt", directory / "candidate.txt"
    reference.write_text("".join(lines), encoding="utf-8")
    candidate.write_text("".join(written), encoding="utf-8")
    return reference, candidate, target


@pytest.mark.parametrize(("corpus", "run"), sorted(MAXENTSCAN_CORPUS_COUNTS))
def test_committed_class_table_applies_maxentscan_overlapping_windows_on_the_fixture_corpora_at_their_count(tmp_path, monkeypatch, corpus, run):
    """Under the committed class table and the corpus's own cache, the scorer sets aside every
    field a MaxEntScan fixture corpus documents under ``maxentscan_overlapping_windows``, at
    the count the corpus carries, and the adjusted counts agree; without the cache the same
    fields are differences, and a score changed on a row off the two windows stays one."""
    # The comparator's own classes and the registry's own plugins, in place of the planted ones.
    monkeypatch.undo()
    by_feature, by_allele = MAXENTSCAN_CORPUS_COUNTS[(corpus, run)]
    tool = built_pairing_tool() if THROUGH_PAIRING_TOOL else None

    def scored(name: str, cache: Path | None, changed_score: bool = False) -> tuple[dict, dict, dict]:
        reference, candidate, target = documented_candidate(GOLDEN_CORPORA / corpus, run, tmp_path / name, changed_score)
        report = spt.score(
            reference, [("vep-rs", candidate)], ["MaxEntScan"], {}, tmp_path / name / "work", sort_memory="16M",
            class_table=COMMITTED_CLASS_TABLE, met_alone=True, json_cache=cache, pairing_tool=tool,
        )  # fmt: skip
        return entry(report, "MaxEntScan", "vep-rs"), entry(report, "MaxEntScan", "vep-rs", "allele"), target

    cache = GOLDEN_CORPORA / corpus / "json_cache"
    feature, allele, target = scored("with", cache)
    assert target["field_divergence_summary"] == {MAXENTSCAN_CLASS: by_feature} and len(target["field_divergences"]) == by_feature
    assert feature["masked"] == {MAXENTSCAN_CLASS: by_feature} and allele["masked"] == {MAXENTSCAN_CLASS: by_allele}
    for item in (feature, allele):
        assert item["adj_perl"] == item["adj_rust"] == item["adj_intersection"] and item["residual_rows"] == 0
        assert item["not_applied"] == {} and item["perl"] == item["rust"] > item["intersection"]
    assert feature["perl"] - feature["intersection"] == by_feature
    without, _, _ = scored("without", None)
    assert without["masked"] == {MAXENTSCAN_CLASS: 0} and without["adj_perl"] - without["adj_intersection"] == by_feature
    changed, _, _ = scored("changed", cache, changed_score=True)
    assert changed["masked"] == {MAXENTSCAN_CLASS: by_feature} and changed["adj_perl"] - changed["adj_intersection"] == 1
    assert changed["residual_rows"] > 0
