"""Tests for the release figure's values layer, on records each test writes into a temporary directory and
on the fixture plugin record beside this file: the values file's bytes and keys, the four panels' values,
the texts, the registry, every refusal by its message, and the pairing of the committed figure with the
README. No test draws: the image is checked through its values file and its stamp."""

from __future__ import annotations

import ast
import copy
import hashlib
import itertools
import json
import math
import re
import shutil
import struct
import subprocess
import sys
import zlib
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
import render_release_figure as figure

SCRIPT = HERE / "render_release_figure.py"
REGISTRY = HERE / "comparators.tsv"
FIXTURE = HERE / "testdata" / "plugin_record_fixture.json"
PLUGINS = json.loads(FIXTURE.read_text(encoding="utf-8"))
# What the fixture's support blocks read, as the texts state it: the released tool's count, then the comparator's.
READS = " and ".join(f"{block['read']} of {block['of']}" for block in PLUGINS["support"])
# The synthetic records state the fixture's release and tools, so the set draws without --draft.
RELEASE = PLUGINS["release"]
COMPARATOR_VERSION = next(e["engine_version"] for e in PLUGINS["engines"] if e["id"] == "fastvep")
REFERENCE_VERSION = PLUGINS["perl_reference"][0]["ensembl_vep"]
IMAGE_TAG = f"ensemblorg/ensembl-vep:release_{REFERENCE_VERSION}"
HEADLINE = "--exclude-readthrough --exclude-artifact"
PROVENANCE = "docs/concordance-provenance"
ARCHES = ("arm64", "x86_64")
POOLED = "The four whole-genome datasets, pooled"
# The plugin panel's group of rows scored on real data, named for what the fixture's rows of the released tool sum:
# the ClinVar dataset and the population dataset of the row's assembly, whose labels open with one word.
REAL_GROUP = "Genome-wide ClinVar and Population"
# The runs a wall-time cell of the synthetic records holds: the reference's, and the two other tools'.
REFERENCE_RUNS, TOOL_RUNS = 1, 3
RUNS_CLAUSE = f"each wall-time cell of Ensembl VEP holds {REFERENCE_RUNS} run, against {TOOL_RUNS} for vep-rs and fastVEP"
IMAGE_LINE = re.compile(r"^!\[[^\]]*\]\(docs/assets/comparison\.png\)\s*$", re.MULTILINE)


def nine_places(matched: int, reference: int, engine: int) -> str:
    """2I/(A+B) rounded half up to nine places, in integer arithmetic."""
    scaled = (4 * matched * 10**9 + reference + engine) // (2 * (reference + engine))
    return f"{scaled // 10**9}.{scaled % 10**9:09d}"


def cell(matched: int, reference: int, engine: int) -> dict:
    """A release record's concordance cell whose adjusted counts are all the matched count."""
    return {"intersection": matched, "perl": reference, "rust": engine,
            "adj_intersection": matched, "adj_perl": matched, "adj_rust": matched}


def file_counts(matched: int, reference: int, engine: int) -> dict:
    """One file's counts of a structural-variant cell."""
    return {"intersection": matched, "perl_tuples": reference, "engine_tuples": engine}


def release_record(release: str = RELEASE) -> dict:
    scored = [
        {"dataset": "pop_a", "label": "Population A, whole genome", "assembly": "GRCh38", "intersection": 99,
         "perl_tuples": 100, "vep_rs_tuples": 100, "excluded_perl": 1, "excluded_rust": 1},
        {"dataset": "pop_b", "label": "Population B, whole genome", "assembly": "GRCh37", "intersection": 7,
         "perl_tuples": 8, "vep_rs_tuples": 8, "excluded_perl": 1, "excluded_rust": 1},
    ]
    medians = {("s01", "arm64"): 1.5, ("s01", "x86_64"): 2.345, ("s02", "arm64"): 3.0, ("s02", "x86_64"): 4.125}
    aggregates = [{"arch": arch, "suite_id": cell_id, "engine": "vep_rs", "median_sec": median, "n_clones": TOOL_RUNS}
                  for (cell_id, arch), median in medians.items()]
    aggregates.append({"arch": "arm64", "suite_id": "s01", "engine": "vep_rs", "median_sec": 99.0, "n_clones": TOOL_RUNS, "target_cpu": "generic"})
    return {
        "record": "release-addendum", "release": release, "source_commit": "1" * 40,
        "perl_reference": {"snv_indel": {"image": IMAGE_TAG}, "sv": {"image": IMAGE_TAG}},
        "concordance": {
            "s01": cell(3, 4, 4), "s02": cell(9, 10, 10),
            "s07": {**cell(1, 2, 2), "per_file": {"one": file_counts(1, 1, 1), "two": file_counts(0, 1, 1)}},
            "s08": {**cell(4, 4, 4), "per_file": {"one": file_counts(3, 3, 3), "two": file_counts(1, 1, 1)}},
        },
        "population_concordance": {
            "per_dataset": scored,
            "pooled": {"label": "All two", "intersection": 106, "perl_tuples": 108, "vep_rs_tuples": 108,
                       "excluded_perl": 2, "excluded_rust": 2},
        },
        "per_cell_aggregates": aggregates,
    }


def population_record(release: str = RELEASE, comparator: str = COMPARATOR_VERSION, args: str = HEADLINE) -> dict:
    def sums(base: float, runs: int | None = None) -> list[dict]:
        """A tool's sums; each entry states the runs of its cells where `runs` is given, as the tools' entries do."""
        return [{"dataset": dataset, "arch": arch, "sum_of_cell_medians_sec": base * factor + step} | ({} if runs is None else {"n_clones_per_cell": runs})
                for dataset, factor in (("pop_a", 10), ("pop_b", 1)) for arch, step in zip(ARCHES, (0.4, 0.5))]

    return {
        "record": "population-measurements", "release": release,
        "datasets": [{"dataset": "pop_a", "label": "Population A, whole genome", "assembly": "GRCh38"},
                     {"dataset": "pop_b", "label": "Population B, whole genome", "assembly": "GRCh37"}],
        "wall_time": {
            "vep_rs": {"engine_version": release, "n_clones_per_cell": TOOL_RUNS, "per_dataset": sums(10, TOOL_RUNS)},
            "fastvep": {"engine_version": comparator, "engine_extra_args": args, "n_clones_per_cell": TOOL_RUNS, "per_dataset": sums(100, TOOL_RUNS)},
            "ensembl_vep": {"engine_version": f"release_{REFERENCE_VERSION}", "n_clones_per_cell": REFERENCE_RUNS, "per_dataset": sums(1000)},
        },
        "comparator_concordance": {
            "engine": "fastvep", "engine_version": comparator, "engine_extra_args": args,
            "per_dataset": [
                {"dataset": "pop_a", "intersection": 50, "perl_tuples": 100, "engine_tuples": 100, "excluded_perl": 0, "excluded_rust": 0},
                {"dataset": "pop_b", "intersection": 6, "perl_tuples": 8, "engine_tuples": 8, "excluded_perl": 1, "excluded_rust": 1},
            ],
            "pooled": {"label": "All two", "intersection": 56, "perl_tuples": 108, "engine_tuples": 108},
        },
    }


def dataset_record(release: str = RELEASE, args: str | None = HEADLINE, matched: int = 1) -> dict:
    record = {
        "record": "comparator-suites", "engine": "fastvep", "engine_version": COMPARATOR_VERSION,
        "scored_with_release": {"version": release},
        "concordance": {name: {"intersection": matched, "perl": 4, "engine_tuples": 4,
                               "adjusted": {"intersection": matched, "perl": 4, "engine_tuples": 4}} for name in ("s01", "s02", "s07", "s08")},
    }
    if args is not None:
        record["engine_extra_args"] = args
    return record


def walltime_record(release: str = RELEASE, args: str | None = HEADLINE) -> dict:
    engine = {"version": COMPARATOR_VERSION}
    if args is not None:
        engine["engine_extra_args"] = args
    return {
        "record": "comparator-walltime", "release": release,
        "engines": {"fastvep": engine, "perl": {"version": f"release_{REFERENCE_VERSION}"}},
        "per_cell_aggregates": [
            {"engine": name, "suite_id": cell_id, "arch": arch, "median_sec": base + step, "n_clones": runs}
            for name, base, runs in (("fastvep", 20.0, TOOL_RUNS), ("perl", 700.0, REFERENCE_RUNS)) for cell_id in ("s01", "s02")
            for step, arch in enumerate(ARCHES)
        ],
    }


def record_set(**replaced) -> dict:
    """The records of one release by file name: the four synthetic ones and the fixture plugin record."""
    records = {
        f"2001-01-01-release-{RELEASE}.json": release_record(),
        f"2001-01-01-population-{RELEASE}.json": population_record(),
        f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json": dataset_record(),
        f"2001-01-01-comparator-walltime-{RELEASE}.json": walltime_record(),
        f"2001-01-01-plugins-{RELEASE}.json": copy.deepcopy(PLUGINS),
    }
    for kind, record in replaced.items():
        name = next(name for name in records if f"-{kind}-" in name)
        if record is None:
            del records[name]
        else:
            records[name] = record
    return records


def write_records(directory: Path, records: dict) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    for name, record in records.items():
        (directory / name).write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    return directory


def run(tmp_path: Path, records: dict, *flags: str, registry: Path = REGISTRY, out: str = "out/values.json"):
    """The script on `records` with --values-only; returns the finished process and the values path."""
    directory = write_records(tmp_path / "records", records)
    values = tmp_path / out
    done = subprocess.run(
        [sys.executable, str(SCRIPT), "--repo", str(tmp_path / "repo"), "--records", str(directory), "--registry", str(registry),
         "--values", str(values), "--values-only", *flags],
        capture_output=True, text=True, check=False,
    )
    return done, values


def values_of(tmp_path: Path, records: dict | None = None, *flags: str, registry: Path = REGISTRY) -> dict:
    done, values = run(tmp_path, record_set() if records is None else records, *flags, registry=registry)
    assert done.returncode == 0, done.stderr
    return json.loads(values.read_text(encoding="utf-8"))


def refusal(tmp_path: Path, records: dict, *flags: str, registry: Path = REGISTRY) -> str:
    """The one line the script refuses `records` with; no values file is written."""
    done, values = run(tmp_path, records, *flags, registry=registry)
    assert done.returncode == 1 and not values.exists(), (done.returncode, done.stdout)
    assert done.stderr.count("\n") == 1 and done.stderr.startswith("ERROR: [render_release_figure] ")
    return done.stderr.strip().removeprefix("ERROR: [render_release_figure] ")


def panel(document: dict, panel_id: str) -> dict:
    return next(p for p in document["panels"] if p["id"] == panel_id)


def rows_of(document: dict, panel_id: str) -> dict:
    """A panel's rows by (group label, row label)."""
    return {(g["label"], r["label"]): r for g in panel(document, panel_id)["groups"] for r in g["rows"]}


def marks(row: dict) -> dict:
    return {(v["tool"], v["kind"]): v["value"] for v in row["values"]}


def registry_with(tmp_path: Path, edit) -> Path:
    """A copy of the registry whose rows `edit` rewrote: it receives {id: {column: cell}} and may add rows."""
    header, *lines = REGISTRY.read_text(encoding="utf-8").splitlines()
    columns = header.split("\t")
    rows = {line.split("\t")[0]: dict(zip(columns, line.split("\t"))) for line in lines}
    edit(rows)
    path = tmp_path / "registry.tsv"
    path.write_text("\n".join([header, *("\t".join(row[c] for c in columns) for row in rows.values())]) + "\n", encoding="utf-8")
    return path


# Numbers and keys


def test_f1_is_rounded_half_up_to_nine_places_from_the_counts():
    assert figure.f1(3, 4, 4, "x") == "0.750000000"
    assert figure.f1(5, 5, 5, "x") == "1.000000000"
    assert figure.f1(1, 2 * 10**9, 2 * 10**9, "x") == "0.000000001"
    assert figure.f1(269568, 269598, 269598, "x") == nine_places(269568, 269598, 269598) == "0.999888723"


def test_wall_times_are_rounded_half_up_from_the_decimal_text():
    assert figure.seconds(2378.345, 2, "x") == "2378.35"
    assert figure.seconds(121.125, 2, "x") == "121.13"
    assert figure.seconds(507.5, 0, "x") == "508"
    assert figure.seconds(7, 2, "x") == "7.00"


def test_two_runs_on_the_same_records_write_the_same_bytes(tmp_path):
    first, first_values = run(tmp_path / "a", record_set())
    second, second_values = run(tmp_path / "b", record_set())
    assert first.returncode == second.returncode == 0, first.stderr + second.stderr
    assert first_values.read_bytes() == second_values.read_bytes()
    assert first_values.read_bytes().endswith(b"\n")


def test_values_file_carries_the_fixed_keys_and_every_record_with_its_digest(tmp_path):
    records = record_set()
    document = values_of(tmp_path, records)
    assert list(document) == ["image", "canvas", "alt_text", "caption", "notes", "sources", "tools", "records", "panels"]
    assert document["image"] == "docs/assets/comparison.png"
    assert [p["id"] for p in document["panels"]] == ["a", "b", "c", "d"]
    written = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in (tmp_path / "records").iterdir()}
    assert {r["file"]: r["sha256"] for r in document["records"]} == written
    files = {r["file"] for r in document["records"]}
    for p in document["panels"]:
        for group in p["groups"]:
            for row in group["rows"]:
                for value in row["values"]:
                    # An F1 carries the counts it is computed from, and prints what they give at nine places.
                    counted = p["id"] in ("a", "c") and value["kind"] != "text"
                    assert set(value) == {"tool", "kind", "value", "record", "key"} | ({"counts"} if counted else set())
                    assert value["kind"] in ("bar", "marker", "text") and value["key"]
                    # A value's source is a record of the set, or two of them where it draws on two.
                    assert all(name in files for name in value["record"].split(" and "))
                    if counted:
                        assert all(type(number) is int for number in value["counts"]) and value["value"] == nine_places(*value["counts"])


def test_canvas_is_four_by_five(tmp_path):
    canvas = values_of(tmp_path)["canvas"]
    assert canvas == {"width_px": 1080, "height_px": 1350}
    assert canvas["width_px"] * 5 == canvas["height_px"] * 4


# The four panels


def test_panel_a_states_each_dataset_from_its_counts(tmp_path):
    rows = rows_of(values_of(tmp_path), "a")
    assert [group for group, _ in rows] == (["Genome-wide ClinVar"] * 2 + ["Whole-genome population datasets"] * 2
                                            + [POOLED] + ["Structural variants"] * 2)
    assert marks(rows[("Genome-wide ClinVar", "ClinVar GRCh37")]) == {
        ("vep-rs", "bar"): "0.750000000", ("vep-rs", "marker"): "1.000000000", ("fastvep", "bar"): "0.250000000",
        ("fastvep", "marker"): "0.250000000"}
    assert marks(rows[("Whole-genome population datasets", "Population A GRCh38")]) == {
        ("vep-rs", "bar"): "0.990000000", ("vep-rs", "marker"): "1.000000000", ("fastvep", "bar"): "0.500000000",
        ("fastvep", "marker"): "0.500000000"}
    # The pool sums the four whole-genome rows' integer counts, as the README pools them: for vep-rs
    # 3 + 9 + 99 + 7 matched of 4 + 10 + 100 + 8 on each side, for fastVEP 1 + 1 + 50 + 6 of 4 + 4 + 100 + 8.
    pooled = rows[(POOLED, "All four")]
    assert marks(pooled) == {("vep-rs", "bar"): nine_places(118, 122, 122), ("vep-rs", "marker"): "1.000000000", ("fastvep", "bar"): "0.500000000",
                             ("fastvep", "marker"): nine_places(58, 115, 115)}
    assert nine_places(118, 122, 122) == "0.967213115"
    cells = "concordance.s01, concordance.s02"
    assert [v["key"] for v in pooled["values"]] == [
        f"summed over population_concordance.per_dataset (intersection, perl_tuples, vep_rs_tuples) and {cells} (intersection, perl, rust)",
        f"summed over population_concordance.per_dataset (intersection, perl_tuples - excluded_perl, vep_rs_tuples - excluded_rust) and {cells} "
        + "(adj_intersection, adj_perl, adj_rust)",
        "summed over comparator_concordance.per_dataset (intersection, perl_tuples, engine_tuples) and, in "
        + f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json, {cells} (intersection, perl, engine_tuples)",
        "summed over comparator_concordance.per_dataset (intersection, perl_tuples - excluded_perl, engine_tuples - excluded_rust) and, in "
        + f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json, the adjusted blocks of {cells} (intersection, perl, engine_tuples)",
    ]
    assert [v["record"] for v in pooled["values"]] == [f"2001-01-01-release-{RELEASE}.json"] * 2 + [
        f"2001-01-01-population-{RELEASE}.json and 2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json"] * 2
    assert marks(rows[("Structural variants", "SV per-VCF, 2 files pooled GRCh37")]) == {
        ("vep-rs", "bar"): "0.500000000", ("vep-rs", "marker"): "1.000000000", ("fastvep", "bar"): "0.250000000",
        ("fastvep", "marker"): "0.250000000"}
    bar = next(v for v in rows[("Genome-wide ClinVar", "ClinVar GRCh37")]["values"] if v["tool"] == "vep-rs" and v["kind"] == "bar")
    assert bar["key"] == "concordance.s01: intersection, perl, rust" and bar["record"] == f"2001-01-01-release-{RELEASE}.json"


def test_panel_b_states_wall_time_per_dataset_and_architecture(tmp_path):
    rows = rows_of(values_of(tmp_path), "b")
    assert list(rows) == [(group, arch) for group in ("ClinVar GRCh37", "ClinVar GRCh38", "Population A GRCh38", "Population B GRCh37")
                          for arch in ("ARM", "x86")]
    # A median is stated to two places and a sum to the second, both rounded half up; the release record's
    # aggregate of a build for another target is never read.
    assert marks(rows[("ClinVar GRCh37", "ARM")]) == {("vep-rs", "marker"): "1.50", ("fastvep", "marker"): "20.00", ("ensembl-vep", "marker"): "700.00"}
    assert marks(rows[("ClinVar GRCh37", "x86")])[("vep-rs", "marker")] == "2.35"
    assert marks(rows[("ClinVar GRCh38", "x86")])[("vep-rs", "marker")] == "4.13"
    assert marks(rows[("Population A GRCh38", "x86")]) == {("vep-rs", "marker"): "101", ("fastvep", "marker"): "1001", ("ensembl-vep", "marker"): "10001"}
    assert marks(rows[("Population B GRCh37", "ARM")]) == {("vep-rs", "marker"): "10", ("fastvep", "marker"): "100", ("ensembl-vep", "marker"): "1000"}


def test_the_key_of_the_released_tools_clinvar_wall_time_names_the_one_entry_without_a_cpu_target(tmp_path):
    release, document = release_record(), values_of(tmp_path)
    rows = rows_of(document, "b")
    selector = re.compile(r"^per_cell_aggregates\[suite_id=(s0[12]),arch=(\w+),no target_cpu\]\.median_sec$")
    checked = []
    for (group, arch_label), row in rows.items():
        if not group.startswith("ClinVar"):
            continue
        ours = next(value for value in row["values"] if (value["tool"], value["kind"]) == ("vep-rs", "marker"))
        # The key names the cell, the architecture and the condition, and exactly one entry of the release record
        # fits all three: its median is the value.
        cell_id, arch = selector.match(ours["key"]).groups()
        fits = [a for a in release["per_cell_aggregates"] if (a["suite_id"], a["arch"]) == (cell_id, arch)]
        named = [a for a in fits if "target_cpu" not in a]
        assert len(named) == 1 and Decimal(ours["value"]) == Decimal(str(named[0]["median_sec"])).quantize(Decimal("0.01"), "ROUND_HALF_UP")
        checked.append((cell_id, arch, len(fits)))
        # A comparator's key names its engine, cell and architecture, which fit one entry of its record.
        theirs = next(value for value in row["values"] if (value["tool"], value["kind"]) == ("fastvep", "marker"))
        assert theirs["key"] == f"per_cell_aggregates[engine=fastvep,suite_id={cell_id},arch={arch}].median_sec"
    # The record holds a second entry of one cell and architecture, built for a named CPU target: without the
    # condition the key would fit both, and the figure draws the entry that names none.
    assert checked == [("s01", "arm64", 2), ("s01", "x86_64", 1), ("s02", "arm64", 1), ("s02", "x86_64", 1)]
    other = next(a for a in release["per_cell_aggregates"] if "target_cpu" in a)
    assert (other["suite_id"], other["arch"], other["median_sec"]) == ("s01", "arm64", 99.0)
    assert marks(rows[("ClinVar GRCh37", "ARM")])[("vep-rs", "marker")] == "1.50"
    # A sum of a population dataset is one entry of its list by dataset and architecture, and its key names the dataset and the architecture alone.
    summed = next(value for value in rows[("Population A GRCh38", "ARM")]["values"] if value["tool"] == "vep-rs")
    assert summed["key"] == "wall_time.vep_rs.per_dataset[dataset=pop_a,arch=arm64].sum_of_cell_medians_sec"


def plugin_sums(list_name: str, plugin: str, assembly: str, engine: str, projection: str, keys=("intersection", "perl", "rust"),
                population: bool = True):
    """The fixture's counts under `keys` summed over the entries of one row, engine and projection, or None
    where the fixture holds no such entry or one of them lacks a key. A row of `real_data` sums its entries on
    the population datasets too, as the figure does, unless `population` is false: the ClinVar entries alone."""
    lists = (list_name, "population_data") if population and list_name == "real_data" else (list_name,)
    entries = [e for name in lists for e in PLUGINS[name] if (e["plugin"], e["assembly"], e["engine"], e["projection"]) == (plugin, assembly, engine, projection)]
    return [sum(e[k] for e in entries) for k in keys] if entries and all(k in e for e in entries for k in keys) else None


def test_fixture_record_drives_panel_c(tmp_path):
    document = values_of(tmp_path)
    groups = panel(document, "c")["groups"]
    assert [g["label"] for g in groups] == [
        REAL_GROUP, "Committed table", "Synthetic fixtures (data that cannot be redistributed)", "Custom annotation (--custom sources)"]
    lists = dict(zip((g["label"] for g in groups), ("real_data", "committed_table", "synthetic_fixtures", "custom_annotation")))
    gate_rows = [(r["plugin"], r["assembly"], r["scored_on"]) for r in PLUGINS["rows"]]
    drawn = [r["label"] for g in groups[:3] for r in g["rows"]]
    assert sorted(drawn) == sorted(f"{p} {a}" + (" synthetic" if kind == "synthetic_fixture" else "") for p, a, kind in gate_rows)
    absent = {(a["plugin"], a["engine"]) for a in PLUGINS["absent"]}
    checked = {"bar": 0, "marker": 0, "text": 0}
    for group in groups:
        for row in group["rows"]:
            plugin, assembly = row["label"].removesuffix(" synthetic").rsplit(" ", 1)
            found = marks(row)
            for engine, kind, projection, keys in (
                ("vep-rs", "bar", "feature", ("intersection", "perl", "rust")),
                ("vep-rs", "marker", "feature", ("adj_intersection", "adj_perl", "adj_rust")),
                ("fastvep", "bar", "feature", ("intersection", "perl", "rust")),
                ("fastvep", "marker", "feature", ("adj_intersection", "adj_perl", "adj_rust")),
            ):
                sums = plugin_sums(lists[group["label"]], plugin, assembly, engine, projection, keys)
                if sums is None:
                    assert (engine, kind) not in found, row["label"]
                else:
                    assert found[(engine, kind)] == nine_places(*sums), row["label"]
                    checked[kind] += 1
            if ("fastvep", "text") in found:
                assert (plugin.split(" (")[0], "fastvep") in absent and found[("fastvep", "text")] == "no such source"
                assert ("fastvep", "bar") not in found
                checked["text"] += 1
    assert checked["bar"] > 57 and checked["marker"] >= 57 and checked["text"] > 0

    rows = rows_of(document, "c")
    # An absent source, a synthetic-fixture row, a committed-table row and a comparator's per-feature bar.
    assert marks(rows[(REAL_GROUP, "AncestralAllele GRCh37")])[("fastvep", "text")] == "no such source"
    synthetic = rows[("Synthetic fixtures (data that cannot be redistributed)", "dbNSFP GRCh37 synthetic")]
    assert marks(synthetic)[("vep-rs", "bar")] == "1.000000000"
    assert next(v["key"] for v in synthetic["values"] if v["kind"] == "bar" and v["tool"] == "vep-rs").endswith("summed over 2 entries")
    assert all(label.endswith(" synthetic") for group, label in rows if group.startswith("Synthetic fixtures"))
    assert not any(label.endswith("synthetic") for group, label in rows if not group.startswith("Synthetic fixtures"))
    assert ("vep-rs", "bar") in marks(rows[("Committed table", "TranscriptAnnotator GRCh37")])
    alpha = marks(rows[(REAL_GROUP, "AlphaMissense GRCh37")])
    assert alpha[("vep-rs", "bar")] == nine_places(*plugin_sums("real_data", "AlphaMissense", "GRCh37", "vep-rs", "feature")) == "0.999940193"
    assert alpha[("vep-rs", "marker")] == "1.000000000"
    assert alpha[("fastvep", "bar")] == nine_places(263768, 269598, 263791) and ("fastvep", "marker") not in alpha


# Rows scored on real data: ClinVar and the population datasets of the row's assembly

ADJUSTED = ("adj_intersection", "adj_perl", "adj_rust")
PLUGINS_NAME = f"2001-01-01-plugins-{RELEASE}.json"
POPULATION_NAME = f"2001-01-01-population-{RELEASE}.json"
# The sentence of both texts on the datasets each tool's rows cover, over the fixture; `{grch38}` are the
# population datasets of GRCh38 that the released tool's rows sum.
COVERAGE = (
    "On a row of a plugin verified on real data, vep-rs's counts are summed over the genome-wide ClinVar dataset of the row's genome build "
    "and whole-genome population datasets of that build (GRCh37: Population B; GRCh38: {grch38}), except on gnomADMt GRCh38 (ClinVar alone); "
    "fastVEP's counts are of the genome-wide ClinVar dataset of the row's genome build alone."
)
UNEVEN_OPENING = "A row's figure is summed over the sets"


def real_entries(plugins: dict, plugin: str, assembly: str, engine: str = "vep-rs", projection: str = "feature") -> list[dict]:
    """The entries a row scored on real data sums for an engine: its ClinVar entries, then its entries on the
    population datasets."""
    return [e for name in ("real_data", "population_data") for e in plugins.get(name, [])
            if (e["plugin"], e["assembly"], e["engine"], e["projection"]) == (plugin, assembly, engine, projection)]


def summed(entries: list[dict], keys=("intersection", "perl", "rust")) -> list[int]:
    return [sum(entry[key] for entry in entries) for key in keys]


def population_entry(plugin: str, assembly: str, dataset: str, counts: tuple[int, int, int], engine: str = "vep-rs", projection: str = "feature") -> dict:
    """A fabricated entry of the plugin record's `population_data`: `counts` are matched, reference and engine,
    and each adjusted count is the matched count."""
    matched, reference, tool = counts
    return {"plugin": plugin, "assembly": assembly, "dataset": dataset, "family": "fabricated", "engine": engine, "projection": projection,
            "scored_on": "population", "perl": reference, "rust": tool, "intersection": matched, "adj_perl": matched, "adj_rust": matched,
            "adj_intersection": matched, "masked": {}, "residual_rows": 0, "text_only_differences": 0}


def second_dataset() -> tuple[dict, dict]:
    """The record set with a second population dataset on GRCh38, `pop_c`: the population record lists it with
    its wall-time sums, and the plugin record scores on it every plugin it scores on `pop_a`, with other counts,
    and states it for the plugin it states `pop_a` for. Returns the records and the plugin record among them."""
    records = record_set()
    population, plugins = records[POPULATION_NAME], records[PLUGINS_NAME]
    population["datasets"].append({"dataset": "pop_c", "label": "Population C, whole genome", "assembly": "GRCh38"})
    for block in population["wall_time"].values():
        block["per_dataset"] += [{**entry, "dataset": "pop_c"} for entry in block["per_dataset"] if entry["dataset"] == "pop_a"]
    for entry in [e for e in plugins["population_data"] if e["dataset"] == "pop_a"]:
        more = {key: entry[key] + 1000 for key in ("perl", "rust")} | {key: entry[key] + 900 for key in ("intersection", *ADJUSTED)}
        plugins["population_data"].append({**entry, "dataset": "pop_c", **more})
    plugins["population_absent"] += [{**item, "dataset": "pop_c"} for item in plugins["population_absent"] if item["dataset"] == "pop_a"]
    return records, plugins


def plugin_row_of(document: dict, label: str) -> dict:
    """The plugin panel's row of one label, whatever its group is named."""
    (row,) = [row for (_, name), row in rows_of(document, "c").items() if name == label]
    return row


def value_of(row: dict, tool: str, kind: str) -> dict:
    return next(value for value in row["values"] if (value["tool"], value["kind"]) == (tool, kind))


def test_a_row_scored_on_real_data_sums_the_released_tools_counts_over_clinvar_and_every_population_dataset_of_its_assembly(tmp_path):
    records, plugins = second_dataset()
    document = values_of(tmp_path, records)
    row = plugin_row_of(document, "AlphaMissense GRCh38")
    entries = real_entries(plugins, "AlphaMissense", "GRCh38")
    assert [entry["dataset"] for entry in entries] == ["clinvar_grch38", "pop_a", "pop_c"]
    bar, marker = value_of(row, "vep-rs", "bar"), value_of(row, "vep-rs", "marker")
    # The bar is the F1 of the three datasets' raw counts summed, the marker that of their adjusted counts summed.
    assert bar["counts"] == summed(entries) and bar["value"] == nine_places(*summed(entries))
    assert marker["counts"] == summed(entries, ADJUSTED) and marker["value"] == "1.000000000"
    assert bar["counts"] != summed(entries[:1]) and marker["counts"] != summed(entries[:1], ADJUSTED)
    # Each key names every entry its figure sums, and how many they are.
    selectors = " and ".join(
        f"{name}[plugin=AlphaMissense,assembly=GRCh38,dataset={dataset},engine=vep-rs,projection=feature]"
        for name, dataset in (("real_data", "clinvar_grch38"), ("population_data", "pop_a"), ("population_data", "pop_c")))
    assert bar["key"] == f"{selectors}: intersection, perl, rust, summed over 3 entries"
    assert marker["key"] == f"{selectors}: adj_intersection, adj_perl, adj_rust, summed over 3 entries"
    # A row of the other assembly sums its own two datasets, and no per-allele entry is summed anywhere.
    other = real_entries(plugins, "AlphaMissense", "GRCh37")
    assert [entry["dataset"] for entry in other] == ["clinvar_grch37", "pop_b"]
    assert value_of(plugin_row_of(document, "AlphaMissense GRCh37"), "vep-rs", "bar")["counts"] == summed(other)
    assert len(real_entries(plugins, "AlphaMissense", "GRCh38", projection="allele")) == 3


def test_a_summed_row_has_no_marker_where_one_of_its_entries_holds_no_adjusted_counts(tmp_path):
    plugins = copy.deepcopy(PLUGINS)
    entry = next(e for e in plugins["population_data"] if (e["plugin"], e["assembly"], e["projection"]) == ("Blosum62", "GRCh37", "feature"))
    for key in ADJUSTED:
        del entry[key]
    done, values = run(tmp_path, record_set(plugins=plugins))
    assert done.returncode == 0, done.stderr
    found = marks(plugin_row_of(json.loads(values.read_text(encoding="utf-8")), "Blosum62 GRCh37"))
    assert found[("vep-rs", "bar")] == nine_places(*summed(real_entries(plugins, "Blosum62", "GRCh37"))) and ("vep-rs", "marker") not in found
    assert done.stdout.count(f"note: {PLUGINS_NAME} holds no adjusted counts of vep-rs on Blosum62 GRCh37; the figure omits that marker\n") == 1


def test_a_comparators_row_reads_the_datasets_the_record_holds_for_it_and_no_other(tmp_path):
    document = values_of(tmp_path)
    plain = plugin_row_of(document, "LOEUF GRCh38")
    clinvar = real_entries(PLUGINS, "LOEUF", "GRCh38", "fastvep")
    # The fixture scores the comparator on ClinVar alone: its bar and marker are that entry's, beside the released tool's sum.
    assert [entry["dataset"] for entry in clinvar] == ["clinvar_grch38"]
    for kind, keys in (("bar", ("intersection", "perl", "rust")), ("marker", ADJUSTED)):
        value = value_of(plain, "fastvep", kind)
        assert value["counts"] == summed(clinvar, keys)
        assert value["key"] == f"real_data[plugin=LOEUF,assembly=GRCh38,dataset=clinvar_grch38,engine=fastvep,projection=feature]: {', '.join(keys)}"
    assert value_of(plain, "vep-rs", "bar")["key"].count("population_data[") == 1
    # An entry the record holds for the comparator on a population dataset is summed into its row, and into no other.
    plugins = copy.deepcopy(PLUGINS)
    plugins["population_data"].append(population_entry("LOEUF", "GRCh38", "pop_a", (700, 1000, 800), engine="fastvep"))
    more = values_of(tmp_path / "more", record_set(plugins=plugins))
    bar = value_of(plugin_row_of(more, "LOEUF GRCh38"), "fastvep", "bar")
    assert bar["counts"] == [held + added for held, added in zip(summed(clinvar), (700, 1000, 800))] and bar["key"].endswith("summed over 2 entries")
    assert value_of(plugin_row_of(more, "LOEUF GRCh37"), "fastvep", "bar") == value_of(plugin_row_of(document, "LOEUF GRCh37"), "fastvep", "bar")
    assert ("fastVEP's counts are summed over the genome-wide ClinVar dataset of the row's genome build and whole-genome population datasets of "
            "that build (GRCh38: Population A), except on ") in more["caption"]


def test_one_sentence_of_each_text_says_which_datasets_the_rows_of_each_tool_cover(tmp_path):
    document = values_of(tmp_path)
    sentence = COVERAGE.format(grch38="Population A")
    for text in (document["caption"], document["alt_text"]):
        assert text.count(sentence) == 1 and text.count("On a row of a plugin verified on real data") == 1
    # In the caption it follows the groups' clauses, which leave the datasets to it.
    caption = document["caption"]
    assert f"marked synthetic; custom-annotation sources on the genome-wide ClinVar datasets. {sentence} The bar is each tool's" in caption
    assert "plugins verified on real data on" not in caption
    # It is drawn from the entries the rows sum and names each dataset as the population record labels it.
    (source,) = [source for source in document["sources"] if source["text"] == sentence]
    assert source == {"text": sentence, "record": f"{PLUGINS_NAME} and {POPULATION_NAME}", "key": figure.COVERAGE_KEY}
    records, _ = second_dataset()
    both = values_of(tmp_path / "second", records)
    assert both["caption"].count(COVERAGE.format(grch38="Population A and Population C")) == 1
    # The sentence on rows whose tools are summed over different numbers of sets counts the sets of a row's own
    # key: a population dataset is none of them, so no row scored on real data is named there.
    uneven = [source["text"] for source in document["sources"] if source["text"].startswith(UNEVEN_OPENING)]
    assert len(uneven) == 1 and not any(f"{row['plugin']} {row['assembly']} (" in uneven[0] for row in PLUGINS["rows"] if row["scored_on"] == "clinvar")
    assert [source["text"] for source in both["sources"] if source["text"].startswith(UNEVEN_OPENING)] == uneven


def test_a_row_of_the_released_tool_that_covers_fewer_datasets_needs_the_records_statement_outside_a_draft(tmp_path):
    # The fixture states that gnomADMt has no count on the population dataset of its assembly: the row is drawn
    # from ClinVar alone, and the sentence names it with what it covers.
    document = values_of(tmp_path)
    bar = value_of(plugin_row_of(document, "gnomADMt GRCh38"), "vep-rs", "bar")
    assert bar["counts"] == summed(real_entries(PLUGINS, "gnomADMt", "GRCh38")) and "population_data[" not in bar["key"]
    assert "except on gnomADMt GRCh38 (ClinVar alone)" in document["caption"]
    # Without the statement the row is refused, and so is a row whose entry on the dataset the record lacks.
    unstated, lacking = copy.deepcopy(PLUGINS), copy.deepcopy(PLUGINS)
    unstated["population_absent"] = []
    lacking["population_data"] = [e for e in lacking["population_data"] if (e["plugin"], e["assembly"]) != ("NMD", "GRCh37")]
    for name, plugins, label, dataset in (("unstated", unstated, "gnomADMt GRCh38", "Population A"), ("lacking", lacking, "NMD GRCh37", "Population B")):
        missing = (f"{PLUGINS_NAME} holds neither a per-feature entry of vep-rs on {label} for {dataset} in population_data nor a statement in "
                   "population_absent that the plugin has no count there")
        assert refusal(tmp_path / name, record_set(plugins=plugins)) == f"{missing}; a row scored on real data covers every population dataset of its assembly"
        # A draft draws the row from the datasets the record holds, with a printed note, and the sentence names it.
        done, values = run(tmp_path / name / "draft", record_set(plugins=plugins), "--draft")
        assert done.returncode == 0, done.stderr
        assert done.stdout.count(f"note: {missing}; the figure draws that row from the datasets the record holds\n") == 1
        drawn = json.loads(values.read_text(encoding="utf-8"))
        assert value_of(plugin_row_of(drawn, label), "vep-rs", "bar")["counts"] == summed(real_entries(plugins, *label.split(" ")))
        assert f"{label} (ClinVar alone)" in drawn["caption"]


@pytest.mark.parametrize("flags", [(), ("--draft",)], ids=["final", "draft"])
@pytest.mark.parametrize(("edit", "message"), [
    (lambda plugins: plugins["population_data"].append(population_entry("AlphaMissense", "GRCh38", "pop_x", (1, 1, 1))),
     f"{PLUGINS_NAME} population_data[36] is on pop_x, which {POPULATION_NAME} does not list among its datasets"),
    (lambda plugins: plugins["population_data"].append(population_entry("AlphaMissense", "GRCh37", "pop_a", (1, 1, 1))),
     f"{PLUGINS_NAME} population_data[36] states GRCh37 for pop_a, which {POPULATION_NAME} lists on GRCh38"),
    (lambda plugins: plugins["population_absent"].append({"plugin": "NMD", "assembly": "GRCh38", "dataset": "pop_x", "reason": "none"}),
     f"{PLUGINS_NAME} population_absent[1] is on pop_x, which {POPULATION_NAME} does not list among its datasets"),
    (lambda plugins: plugins["population_absent"].append({"plugin": "NMD", "assembly": "GRCh37", "dataset": "pop_a", "reason": "none"}),
     f"{PLUGINS_NAME} population_absent[1] states GRCh37 for pop_a, which {POPULATION_NAME} lists on GRCh38"),
    (lambda plugins: plugins["population_absent"].append({"plugin": "NMD", "assembly": "GRCh38", "dataset": "pop_a", "reason": "none"}),
     f"{PLUGINS_NAME} population_absent[1] states that NMD has no count on pop_a, and population_data holds an entry of it"),
    (lambda plugins: plugins["population_absent"][0].pop("reason"), f"{PLUGINS_NAME} population_absent[0] states no reason"),
    (lambda plugins: plugins["population_absent"][0].pop("dataset"), f"{PLUGINS_NAME} population_absent[0] names no plugin, assembly and dataset"),
], ids=["entry-unlisted-dataset", "entry-other-assembly", "statement-unlisted-dataset", "statement-other-assembly", "statement-beside-entry",
        "statement-no-reason", "statement-unnamed"])
def test_an_entry_or_a_statement_on_a_population_dataset_that_the_records_do_not_bear_out_is_refused_by_name(tmp_path, edit, message, flags):
    plugins = copy.deepcopy(PLUGINS)
    edit(plugins)
    assert refusal(tmp_path, record_set(plugins=plugins), *flags) == message


def test_the_group_of_rows_scored_on_real_data_names_in_plain_words_the_datasets_the_released_tools_rows_sum(tmp_path):
    # A dataset the record labels in several releases goes by the leading words its labels share.
    assert figure.dataset_families(["gnomAD v4.1", "gnomAD v2.1.1", "1000 Genomes Phase 3", "1000 Genomes high-coverage"]) == ["gnomAD", "1000 Genomes"]
    assert figure.dataset_families(["Population A"]) == ["Population A"] and figure.dataset_families([]) == []
    document = values_of(tmp_path)
    assert panel(document, "c")["groups"][0]["label"] == REAL_GROUP == "Genome-wide ClinVar and Population"
    # Datasets whose labels share no opening word are named each by its own label, in the population record's order.
    records = record_set()
    next(d for d in records[POPULATION_NAME]["datasets"] if d["dataset"] == "pop_b")["label"] = "Cohort B, whole genome"
    apart = values_of(tmp_path / "apart", records)
    assert panel(apart, "c")["groups"][0]["label"] == "Genome-wide ClinVar, Population A and Cohort B"
    assert "(GRCh37: Cohort B; GRCh38: Population A)" in apart["caption"]


def test_a_record_without_an_entry_on_a_population_dataset_is_refused_and_a_draft_draws_it_from_clinvar_alone_with_its_notes(tmp_path):
    plugins = copy.deepcopy(PLUGINS)
    del plugins["population_data"], plugins["population_absent"]
    real = [row for row in plugins["rows"] if row["scored_on"] == "clinvar"]
    labels = {"GRCh38": "Population A", "GRCh37": "Population B"}
    assert refusal(tmp_path, record_set(plugins=plugins)).startswith(
        f"{PLUGINS_NAME} holds neither a per-feature entry of vep-rs on {real[0]['plugin']} {real[0]['assembly']} for {labels[real[0]['assembly']]} in population_data")
    done, values = run(tmp_path / "draft", record_set(plugins=plugins), "--draft")
    assert done.returncode == 0, done.stderr
    document = json.loads(values.read_text(encoding="utf-8"))
    # One note per row; every row is drawn from its ClinVar entries alone, under the key it has without a sum.
    notes = [line for line in done.stdout.splitlines() if line.endswith("; the figure draws that row from the datasets the record holds")]
    assert len(notes) == len(real) == 35
    for row in real:
        bar = value_of(plugin_row_of(document, f"{row['plugin']} {row['assembly']}"), "vep-rs", "bar")
        assert bar["counts"] == plugin_sums("real_data", row["plugin"], row["assembly"], "vep-rs", "feature", population=False)
        assert bar["key"].endswith("projection=feature]: intersection, perl, rust") and "population_data[" not in bar["key"]
    # The group is named for ClinVar alone, and the sentence says of each tool that its rows cover ClinVar alone.
    assert panel(document, "c")["groups"][0]["label"] == "Genome-wide ClinVar"
    alone = "counts are of the genome-wide ClinVar dataset of the row's genome build alone"
    assert f"On a row of a plugin verified on real data, vep-rs's {alone}; fastVEP's {alone}." in document["caption"]


def test_no_per_allele_tick_is_drawn_on_any_panel_or_named_in_the_legend_or_the_caption(tmp_path):
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    # The record holds the comparator's per-allele entries, and the figure reads none of them: no value of any
    # panel is a tick, and no key names the allele projection.
    assert any(entry["engine"] == "fastvep" and entry["projection"] == "allele" for entry in PLUGINS["real_data"])
    values = [value for p in document["panels"] for group in p["groups"] for row in group["rows"] for value in row["values"]]
    assert {value["kind"] for value in values} == {"bar", "marker", "text"}
    assert not any("projection=allele" in value["key"] for value in values)
    # A comparator's bar on a plugin row is its per-feature figure, as the released tool's is.
    alpha = rows_of(document, "c")[(REAL_GROUP, "AlphaMissense GRCh37")]
    bar = next(value for value in alpha["values"] if (value["tool"], value["kind"]) == ("fastvep", "bar"))
    assert "engine=fastvep,projection=feature" in bar["key"] and bar["counts"] == plugin_sums("real_data", "AlphaMissense", "GRCh37", "fastvep", "feature")
    # The legend, the strings drawn, the caption and the alternative text name no tick and no per-allele figure.
    assert "tick" not in {entry["mark"] for entry in figure.legend_entries(document)}
    assert not any("tick" in kinds for kinds in figure.LEGEND_KINDS)
    words = [*figure.drawn_texts(document, registry), document["caption"], document["alt_text"], *document["notes"]]
    assert not any(word in text for text in words for word in (" tick", "per-allele", "per allele"))
    # The drawing has no branch for a tick: a row's marks are bars, markers and their dashed lines alone.
    row = f1_row(("vep-rs", "bar", "0.990000000"), ("vep-rs", "marker", "1.000000000"), ("fastvep", "bar", "0.500000000"))
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, 35.0)
    assert [name for name, _ in axes.calls] == ["barh", "plot", "plot", "barh"]
    assert not hasattr(figure, "TICK_ORDER") and '"tick"' not in SCRIPT.read_text(encoding="utf-8")


def test_the_plugin_panel_shows_each_tool_of_the_record_on_every_gate_row(tmp_path):
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    gate_rows = [row for group in panel(document, "c")["groups"][:3] for row in group["rows"]]
    # Every gate row holds the released tool's bar and, for the comparator, its bar or the record's reason.
    assert len(gate_rows) == len(PLUGINS["rows"]) == 57
    for row in gate_rows:
        found = marks(row)
        assert ("vep-rs", "bar") in found and (("fastvep", "bar") in found) != (("fastvep", "text") in found), row["label"]
    absent = [row for row in gate_rows if ("fastvep", "text") in marks(row)]
    assert len(absent) == 43 and {marks(row)[("fastvep", "text")] for row in absent} == {"no such source"}
    # An absent plugin is drawn as a plain cross in the comparator's colour on the comparator's lane, and not as words.
    axes, across = Recorder(), object()
    figure.draw_absent(axes, across, [(0, absent[0]), (1, gate_rows[0])], LANES, FILLS, 261.0)
    ((at, cross),) = axes.drawn("plot")
    assert at[1] == [0.22] and not axes.drawn("text") and not axes.drawn("annotate")
    assert (cross["marker"], cross["color"], cross["transform"], cross["linestyle"]) == ("x", FILLS["fastvep"], across, "none")
    # The legend names the cross once, beside the other marks of an F1, and the words stand nowhere else.
    entries = figure.legend_entries(document)
    assert [(entry["tool"], entry["label"]) for entry in entries if entry["mark"] == "absent"] == [("fastvep", "fastVEP: no such source")]
    assert [entry["mark"] for entry in figure.legend_rows(entries)[1]] == ["raw", "exact", "adjusted", "absent"]
    drawn = figure.drawn_texts(document, registry)
    assert drawn.count("fastVEP: no such source") == 1 and "no such source" not in drawn
    assert 'a cross marks a plugin for which the comparator has no counterpart ("no such source")' in document["caption"]

    # A gate row that neither scores the comparator nor lists the plugin as absent for it is refused by name; a
    # draft notes it and draws no mark of the comparator on that row.
    plugins = copy.deepcopy(PLUGINS)
    plugins["absent"] = [item for item in plugins["absent"] if item["plugin"] != "AncestralAllele"]
    name = f"2001-01-01-plugins-{RELEASE}.json"
    missing = f"{name} neither scores fastVEP per feature on AncestralAllele GRCh37 nor lists the plugin as absent for it"
    assert refusal(tmp_path / "neither", record_set(plugins=plugins)) == f"{missing}; every gate row shows each tool the record scores"
    done, values = run(tmp_path / "draft", record_set(plugins=plugins), "--draft")
    assert done.returncode == 0, done.stderr
    assert f"note: {missing}; the figure draws no mark of fastVEP on that row\n" in done.stdout
    row = rows_of(json.loads(values.read_text(encoding="utf-8")), "c")[(REAL_GROUP, "AncestralAllele GRCh37")]
    assert {value["tool"] for value in row["values"]} == {"vep-rs"}

    # One plain mark has one meaning: two reasons for one comparator's absent plugins are refused by name.
    plugins = copy.deepcopy(PLUGINS)
    plugins["absent"][1]["reason"] = "not ported"
    assert refusal(tmp_path / "reasons", record_set(plugins=plugins)) == (
        f"{name} absent states 2 reasons for fastvep (no such source and not ported); the figure draws one plain mark for a tool's "
        "absent plugins and names one reason for it")


def test_the_cross_of_an_absent_plugin_stands_in_the_room_left_of_the_tick_that_reads_zero(tmp_path):
    document = values_of(tmp_path)
    rows = [row for group in panel(document, "c")["groups"] for row in group["rows"]]
    absent = next(row for row in rows if ("fastvep", "text") in marks(row))
    for width in (float(figure.MIN_PLOT), 261.0, 330.0):
        axes, across = Recorder(), object()
        figure.draw_absent(axes, across, [(0, absent)], LANES, FILLS, width)
        ((at, cross),) = axes.drawn("plot")
        # The cross is placed by its share of the plot's width, in the middle of the room, whatever that width.
        centre = at[0][0] * width
        reach = (cross["markersize"] / 2 + cross["markeredgewidth"] * math.sqrt(0.125)) * figure.DPI / 72
        assert cross["transform"] is across and reach == pytest.approx(figure.absent_reach())
        # Its arms end a pixel and a half or more right of the axis line and short of the tick that reads 0: the
        # whole cross stands where no value has a place, and cannot be taken for a low score.
        assert centre - reach >= 1.5 and centre + reach < figure.f1_room()
        scale = figure.f1_scale(width, 1e-7)
        assert pixels_from_the_edge(1.0, scale) - (centre + reach) > 0
        # It reads as a cross at the README's width, 830 pixels of 1080: six pixels or more across there, with a
        # stroke of a pixel and a half.
        assert 2 * reach * 830 / 1080 >= 6.0 and cross["markeredgewidth"] * figure.DPI / 72 * 830 / 1080 >= 1.5
        # It stands short of the bar on the lane above at the least height a line takes.
        assert reach < (LANES["fastvep"] - LANES["vep-rs"] - 0.3 / 2) * figure.MIN_PITCH
    # A bar of the comparator's, at any figure, ends at or right of that tick: no mark of a figure stands in the room.
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, f1_row(("fastvep", "bar", "0.000000000")))], LANES, FILLS, 0.3, 1e-7, scale)
    ((bar_at, bar),) = axes.drawn("barh")
    assert pixels_from_the_edge(bar["left"] + bar_at[1], scale) == pytest.approx(figure.f1_room())
    assert not hasattr(figure, "ABSENT_INSET")


def test_a_cross_takes_eight_pixel_columns_and_leaves_one_of_ground_before_the_gridline_of_the_tick_that_reads_zero(monkeypatch):
    room, place, columns = figure.f1_room(), figure.absent_place(), figure.absent_columns()
    # The drawing library sets its own cross between two pixels, the ends of its strokes on whole pixels where the
    # stroke rounds to an even number of them, and a stroke's corners reach past its end: four columns to either
    # side of the centre, eight in all, which are six pixels or more at the README's width.
    stroke = figure.ABSENT_STROKE * figure.DPI / 72
    ends = math.floor(figure.ABSENT_SIZE / 2 * figure.DPI / 72 + 0.5)
    assert (figure.ABSENT_MARK, figure.ABSENT_SIZE, figure.ABSENT_STROKE) == ("x", 3.3, 1.0) and math.floor(stroke + 0.5) % 2 == 0
    assert ends == 3 and columns == math.ceil(ends + stroke * math.sqrt(0.125)) == 4 and 2 * columns * 830 / 1080 >= 6.0
    # Its centre stands a whole number of pixels left of the tick that reads 0: its four columns and one of ground.
    assert figure.ABSENT_CLEAR == 1 and room - place == pytest.approx(columns + figure.ABSENT_CLEAR, abs=1e-9)

    def column(at: float) -> int:
        """The pixel the drawing library sets a line or a marker on: its place rounded to a whole pixel."""
        return math.floor(at + 0.5)

    apart = {10: 0, 11: 0}
    for step in range(1000):
        # The plot's left edge on each thousandth of a pixel: the columns of the axis line, the gridline and the cross.
        edge = 192 + step / 1000
        axis, gridline, centre = column(edge), column(edge + room), column(edge + place)
        inked = range(centre - columns, centre + columns)
        # Eight columns of ink, one of ground, then the gridline's, wherever the edge falls.
        assert len(inked) == 8 and inked[-1] + figure.ABSENT_CLEAR + 1 == gridline, edge
        # The room sets the two lines eleven columns apart, or ten: the cross never stands on the axis line's column,
        # and at eleven one column lies between the two, the one the axis line's soft edge falls on.
        assert inked[0] - axis == gridline - axis - 9 >= 1, edge
        apart[gridline - axis] += 1
    assert apart[11] > 900 and apart[10] == 1000 - apart[11] > 0
    # The room is the one that holds them: a marker's half width and four pixels, a pixel more than a cross whose
    # arms end at the tick takes.
    assert figure.ROOM_CLEAR == 4.0 and room == pytest.approx(figure.marker_reach("adjusted") + 4.0)

    # The cross is drawn at that place as the library's own cross, on the plot and in the legend.
    row = {"label": "a row", "values": [{"tool": "fastvep", "kind": "text", "value": "no such source"}]}
    for width in (float(figure.MIN_PLOT), 261.0, 330.0):
        axes, across = Recorder(), object()
        figure.draw_absent(axes, across, [(0, row)], LANES, FILLS, width)
        ((at, cross),) = axes.drawn("plot")
        assert at[0][0] * width == pytest.approx(place) and "fillstyle" not in cross
        assert (cross["marker"], cross["markersize"], cross["markeredgewidth"]) == ("x", figure.ABSENT_SIZE, figure.ABSENT_STROKE)
    lines, patches = type(sys)("matplotlib.lines"), type(sys)("matplotlib.patches")
    lines.Line2D, patches.Patch = type("Line2D", (Handle,), {}), type("Patch", (Handle,), {})
    for module in (type(sys)("matplotlib"), lines, patches):
        monkeypatch.setitem(sys.modules, module.__name__, module)
    shown = figure.legend_handle({"mark": "absent", "tool": "fastvep", "label": "fastVEP: no such source"}, FILLS["fastvep"]).style
    assert (shown["marker"], shown["markersize"], shown["markeredgewidth"], shown["color"]) == ("x", figure.ABSENT_SIZE, figure.ABSENT_STROKE, FILLS["fastvep"])
    assert not hasattr(figure, "absent_mark")


def test_panel_d_equals_the_records_support_block(tmp_path):
    document = values_of(tmp_path)
    rows = [row for group in panel(document, "d")["groups"] for row in group["rows"]]
    support = {block["engine"]: block for block in PLUGINS["support"]}
    assert [row["label"] for row in rows] == ["vep-rs", "fastVEP"]
    for row, engine in zip(rows, ("vep-rs", "fastvep")):
        block = support[engine]
        assert marks(row) == {(engine, "bar"): str(block["read"]), (engine, "text"): f"{block['read']} of {block['of']}"}
        assert [v["key"] for v in row["values"]] == [f"support[engine={engine}].read", f"support[engine={engine}].read, of"]


def test_panel_d_names_the_tools_it_draws_in_its_title_and_carries_the_note_of_its_counts_under_it(tmp_path):
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    counts = panel(document, "d")
    # The title names the tools the panel draws, each a row of it, and not the reference, which has no row.
    assert [row["label"] for group in counts["groups"] for row in group["rows"]] == ["vep-rs", "fastVEP"]
    assert counts["title"] == "Plugins read by vep-rs and fastVEP" and "Ensembl VEP" not in counts["title"]
    # The note that says what the counts are of, with the number of plugins the reference lists, is the panel's
    # own: it is no line of the canvas's notes, and it is drawn once, with the panel.
    supported = next(block["of"] for block in PLUGINS["support"] if block["engine"] == "vep-rs")
    assert counts["note"] == (f"Counts are of the {supported} plugins vep-rs supports, of Ensembl VEP's {PLUGINS['not_ported']['of']}. "
                              "Plugins vep-rs does not port are not counted; the README lists them by reason.")
    assert document["notes"] == [] and [p["id"] for p in document["panels"] if "note" in p] == ["d"]
    drawn = figure.drawn_texts(document, registry)
    assert drawn.count(counts["note"]) == 1 and drawn.index(counts["note"]) > drawn.index(counts["title"])
    # The drawing sets the note under the panel of counts, at that panel's column, before any note of the canvas.
    source = SCRIPT.read_text(encoding="utf-8")
    under = source.index('draw_count(sheet, document, tools, panels["d"], lefts[1], tops[1], columns[1])')
    assert under < source.index("sheet.text(lefts[1], tops[1] + NOTE_GAP + index * NOTE_PITCH, line") < source.index("note_top = max(tops) + NOTE_GAP")
    # A draft keeps its line among the canvas's notes, and the panel keeps its own.
    draft = values_of(tmp_path / "draft", record_set(plugins=ungated()), "--draft")
    assert draft["notes"] == ["Draft: not drawn from one release's final records."] and panel(draft, "d")["note"] == counts["note"]


def test_the_note_under_the_panel_of_counts_reads_one_plugin_where_the_released_tool_supports_one(tmp_path):
    listed = PLUGINS["not_ported"]["of"]
    tail = f" vep-rs supports, of Ensembl VEP's {listed}. Plugins vep-rs does not port are not counted; the README lists them by reason."
    # The note counts the plugins the released tool supports, from the record's support block: "plugin" for one.
    for supported, noun in ((1, "plugin"), (2, "plugins")):
        plugins = copy.deepcopy(PLUGINS)
        for block, read in zip(plugins["support"], (supported, 0)):
            block.update(read=read, of=supported)
        document = values_of(tmp_path / noun, record_set(plugins=plugins))
        note = panel(document, "d")["note"]
        assert note == f"Counts are of the {supported} {noun}{tail}" and "1 plugins" not in note, supported
        assert figure.drawn_texts(document, figure.read_registry(REGISTRY)).count(note) == 1
        assert f"(d) The plugins each tool reads, of the {supported} that vep-rs supports: vep-rs {supported} and fastVEP 0; " in document["caption"]


def test_panel_d_is_drawn_only_when_the_released_tool_and_every_comparator_have_a_count(tmp_path):
    # With a count for each, the released tool and the comparator have a row; the reference has none and its note names it.
    document = values_of(tmp_path)
    rows = [row for group in panel(document, "d")["groups"] for row in group["rows"]]
    assert [row["label"] for row in rows] == ["vep-rs", "fastVEP"]
    assert {value["tool"] for row in rows for value in row["values"]} == {"vep-rs", "fastvep"}
    assert panel(document, "d")["note"].startswith("Counts are of the ") and "Ensembl VEP's" in panel(document, "d")["note"]
    # Without a count for one of them the panel is left out, with one line in the caption and no note of counts.
    line = "The panel of the plugins each tool reads is left out: the records hold no count for {}."
    for kept, uncounted in (("vep-rs", "fastVEP"), ("fastvep", "vep-rs"), (None, "vep-rs and fastVEP")):
        plugins = copy.deepcopy(PLUGINS)
        plugins["support"] = [block for block in plugins["support"] if block["engine"] == kept]
        document = values_of(tmp_path / str(kept), record_set(plugins=plugins))
        assert [p["id"] for p in document["panels"]] == ["a", "b", "c"]
        assert document["caption"].count(line.format(uncounted)) == 1 and "(d)" not in document["caption"]
        assert document["caption"].count("The panel of the plugins each tool reads is left out") == 1
        assert document["alt_text"].startswith("Three panels comparing ") and "(d)" not in document["alt_text"]
        assert not any(note.startswith("Counts are of") for note in document["notes"]) and not any("note" in p for p in document["panels"])
        assert "(d)" not in figure.drawn_texts(document, figure.read_registry(REGISTRY))


def test_panels_c_and_d_are_drawn_only_with_a_plugin_record(tmp_path):
    document = values_of(tmp_path, record_set(plugins=None))
    assert [p["id"] for p in document["panels"]] == ["a", "b"]
    assert document["notes"] == []
    assert "(c)" not in document["caption"] and document["alt_text"].startswith("Two panels comparing ")


def test_a_missing_figure_of_a_comparator_loses_its_bar_with_a_printed_note(tmp_path):
    scored = dataset_record()
    del scored["concordance"]["s07"]
    done, values = run(tmp_path, record_set(fastvep=scored))
    assert done.returncode == 0, done.stderr
    assert "note: no record holds fastVEP's concordance on SV per-VCF, 2 files pooled GRCh37; the figure omits that bar\n" in done.stdout
    rows = rows_of(json.loads(values.read_text(encoding="utf-8")), "a")
    assert marks(rows[("Structural variants", "SV per-VCF, 2 files pooled GRCh37")]) == {("vep-rs", "bar"): "0.500000000", ("vep-rs", "marker"): "1.000000000"}
    assert ("fastvep", "bar") in marks(rows[("Structural variants", "SV per-VCF, 2 files pooled GRCh38")])


def test_a_comparator_without_every_whole_genome_row_loses_its_pooled_bar_with_a_printed_note(tmp_path):
    scored = dataset_record()
    del scored["concordance"]["s02"]
    done, values = run(tmp_path, record_set(fastvep=scored))
    assert done.returncode == 0, done.stderr
    assert "note: no record holds fastVEP's concordance on ClinVar GRCh38; the figure omits that bar\n" in done.stdout
    assert "note: no record holds fastVEP's concordance on every dataset of All four; the figure omits that bar\n" in done.stdout
    rows = rows_of(json.loads(values.read_text(encoding="utf-8")), "a")
    assert marks(rows[(POOLED, "All four")]) == {("vep-rs", "bar"): "0.967213115", ("vep-rs", "marker"): "1.000000000"}


def test_a_dataset_only_the_comparator_is_scored_on_is_noted_and_stays_out_of_its_pool(tmp_path):
    population = population_record()
    population["comparator_concordance"]["per_dataset"].append({"dataset": "pop_c", "intersection": 1, "perl_tuples": 9, "engine_tuples": 9})
    done, values = run(tmp_path, record_set(population=population))
    assert done.returncode == 0, done.stderr
    assert (f"note: 2001-01-01-population-{RELEASE}.json scores fastVEP on pop_c, which the release record's population_concordance does not "
            "list; the figure omits it\n") in done.stdout
    document = json.loads(values.read_text(encoding="utf-8"))
    assert panel(document, "a") == panel(values_of(tmp_path / "plain"), "a")
    assert marks(rows_of(document, "a")[(POOLED, "All four")])[("fastvep", "bar")] == "0.500000000"


# The registry


def test_registry_lists_the_three_tools():
    tools = figure.read_registry(REGISTRY)
    assert [(t["id"], t["name"], t["role"]) for t in tools] == [
        ("vep-rs", "vep-rs", "released"), ("fastvep", "fastVEP", "comparator"), ("ensembl-vep", "Ensembl VEP", "reference")]
    assert [t["fill"] for t in tools] == ["#0072B2", "#E69F00", "#009E73"]
    assert [(t["concordance"], t["speed"]) for t in tools] == [(True, True), (True, True), (False, True)]
    assert [t["headline"] for t in tools] == ["", HEADLINE, ""]


def test_a_registry_row_that_forbids_a_concordance_panel_keeps_the_tool_out_of_panels_a_and_c(tmp_path):
    registry = registry_with(tmp_path, lambda rows: rows["fastvep"].update(concordance_panel="no"))
    document = values_of(tmp_path, registry=registry)
    for panel_id in ("a", "c"):
        tools = {v["tool"] for g in panel(document, panel_id)["groups"] for r in g["rows"] for v in r["values"]}
        assert tools == {"vep-rs"}, panel_id
    assert "fastvep" in {v["tool"] for g in panel(document, "b")["groups"] for r in g["rows"] for v in r["values"]}
    assert "fastvep" in {v["tool"] for g in panel(document, "d")["groups"] for r in g["rows"] for v in r["values"]}


def test_a_registry_row_that_forbids_a_speed_panel_keeps_the_tool_out_of_panel_b(tmp_path):
    registry = registry_with(tmp_path, lambda rows: rows["ensembl-vep"].update(speed_panel="no"))
    document = values_of(tmp_path, registry=registry)
    assert {v["tool"] for g in panel(document, "b")["groups"] for r in g["rows"] for v in r["values"]} == {"vep-rs", "fastvep"}


def test_no_tool_outside_the_registry_is_named(tmp_path):
    records = record_set()
    population = next(r for name, r in records.items() if "-population-" in name)
    population["wall_time"]["othertool"] = copy.deepcopy(population["wall_time"]["fastvep"])
    plugins = next(r for name, r in records.items() if "-plugins-" in name)
    plugins["engines"].append({**plugins["engines"][1], "id": "othertool", "name": "OtherTool"})
    plugins["support"].append({**plugins["support"][1], "engine": "othertool"})
    done, values = run(tmp_path, records)
    assert done.returncode == 0, done.stderr
    assert "times othertool, which the registry does not list; the figure omits it" in done.stdout
    assert "scores othertool, which the registry does not list; the figure omits it" in done.stdout
    text = values.read_text(encoding="utf-8")
    document = json.loads(text)
    registry = {tool["id"] for tool in figure.read_registry(REGISTRY)}
    assert {tool["id"] for tool in document["tools"]} == registry
    assert {v["tool"] for p in document["panels"] for g in p["groups"] for r in g["rows"] for v in r["values"]} == registry
    assert "othertool" not in text.lower()


def test_an_engine_of_a_wall_time_record_that_the_registry_does_not_list_is_noted_and_left_out(tmp_path):
    timed = walltime_record()
    timed["engines"]["othertool"] = {"version": "v9.9.9"}
    timed["per_cell_aggregates"] += [{"engine": "othertool", "suite_id": cell_id, "arch": arch, "median_sec": 5.0}
                                     for cell_id in ("s01", "s02") for arch in ARCHES]
    done, values = run(tmp_path, record_set(walltime=timed))
    assert done.returncode == 0, done.stderr
    assert (f"note: 2001-01-01-comparator-walltime-{RELEASE}.json times othertool, which the registry does not list; the figure omits it\n"
            in done.stdout)
    document, plain = json.loads(values.read_text(encoding="utf-8")), values_of(tmp_path / "plain")
    assert (document["panels"], document["tools"], document["caption"]) == (plain["panels"], plain["tools"], plain["caption"])


def test_the_headline_run_is_read_when_the_records_hold_two_runs_of_a_tool(tmp_path):
    records = record_set()
    records[f"2001-01-02-comparator-fastvep-{COMPARATOR_VERSION}-no-extra-args.json"] = dataset_record(args=None, matched=4)
    records[f"2001-01-02-comparator-walltime-{RELEASE}-no-extra-args.json"] = walltime_record(args=None)
    document = values_of(tmp_path, records)
    assert [r["file"] for r in document["records"]] == [
        f"2001-01-01-release-{RELEASE}.json", f"2001-01-01-population-{RELEASE}.json",
        f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json", f"2001-01-01-comparator-walltime-{RELEASE}.json",
        f"2001-01-01-plugins-{RELEASE}.json"]
    assert marks(rows_of(document, "a")[("Genome-wide ClinVar", "ClinVar GRCh37")])[("fastvep", "bar")] == "0.250000000"


# The texts


def test_every_mark_names_its_tool_and_the_version_its_record_states(tmp_path):
    document = values_of(tmp_path)
    assert document["tools"] == [
        {"id": "vep-rs", "name": "vep-rs", "version": RELEASE, "colour": "#0072B2"},
        {"id": "fastvep", "name": "fastVEP", "version": COMPARATOR_VERSION, "colour": "#E69F00"},
        {"id": "ensembl-vep", "name": "Ensembl VEP", "version": REFERENCE_VERSION, "colour": "#009E73"},
    ]
    tools = figure.read_registry(REGISTRY)
    legend = {entry["mark"]: entry["label"] for entry in figure.legend_entries(document)}
    labels = [entry["label"] for entry in figure.legend_entries(document)]
    assert labels[:3] == [f"vep-rs {RELEASE}", f"fastVEP {COMPARATOR_VERSION}", f"Ensembl VEP {REFERENCE_VERSION}"]
    assert legend["exact"] == "Adjusted F1, exactly 1"
    assert (legend["filled"], legend["open"]) == ("ARM", "x86") and legend["adjusted"] == "Adjusted F1, below 1"
    for panel_id in ("a", "c"):
        assert figure.axis_lines(document, tools, panel(document, panel_id))[0] == "F1 against Ensembl VEP"
    script = SCRIPT.read_text(encoding="utf-8")
    assert not re.search(r"\d+\.\d+\.\d+", script), "a version is typed in the script"
    for version in (RELEASE, COMPARATOR_VERSION, REFERENCE_VERSION):
        assert version not in script and version.removeprefix("v") not in script, f"{version} is typed in the script"


def test_a_legend_entry_names_its_tool_and_version_and_not_the_place_of_its_bar(tmp_path):
    # Two tools with a bar in every F1 row.
    document = values_of(tmp_path)
    entries = figure.legend_entries(document)
    assert [entry["label"] for entry in entries if entry["mark"] in figure.LEGEND_KINDS[0]] == [
        f"vep-rs {RELEASE}", f"fastVEP {COMPARATOR_VERSION}", f"Ensembl VEP {REFERENCE_VERSION}"]
    # No entry says where a bar stands in its row: the colour says which bar is which.
    assert not any(word in entry["label"] for entry in entries if entry["mark"] in figure.LEGEND_KINDS[0] for word in ("bar", "circle", "upper", "lower", "top"))
    assert not any(word in entry["label"] for entry in entries for word in ("upper", "lower", "top"))


def test_an_adjusted_f1_below_one_has_its_own_legend_entry(tmp_path):
    release = release_record()
    release["concordance"]["s01"]["adj_perl"] = 4
    document = values_of(tmp_path, record_set(release=release))
    assert marks(rows_of(document, "a")[("Genome-wide ClinVar", "ClinVar GRCh37")])[("vep-rs", "marker")] == nine_places(3, 4, 3)
    assert marks(rows_of(document, "a")[(POOLED, "All four")])[("vep-rs", "marker")] == nine_places(118, 119, 118)
    # The markers of one tool are named as its own and wear its colour; with a comparator's beside them they are
    # named for what they state, in no tool's colour.
    alone = values_of(tmp_path / "alone", record_set(release=release), registry=registry_with(
        tmp_path, lambda rows: rows["fastvep"].update(concordance_panel="no")))
    legend = {entry["mark"]: entry for entry in figure.legend_entries(alone)}
    assert (legend["adjusted"]["label"], legend["exact"]["label"]) == ("vep-rs adjusted F1, below 1", "vep-rs adjusted F1, exactly 1")
    assert legend["adjusted"]["tool"] == legend["exact"]["tool"] == "vep-rs"
    legend = {entry["mark"]: entry for entry in figure.legend_entries(document)}
    assert (legend["adjusted"]["label"], legend["exact"]["label"]) == ("Adjusted F1, below 1", "Adjusted F1, exactly 1")
    assert legend["adjusted"]["tool"] is legend["exact"]["tool"] is None


def test_every_f1_prints_what_its_counts_give_at_nine_places_and_is_placed_from_the_counts(tmp_path):
    # One tuple of ten thousand million unmatched: the raw F1 prints as 1.000000000 at nine places, as the README
    # prints it, and is not 1. The adjusted counts are equal, so the adjusted F1 is exactly 1.
    release = release_record()
    many = 10**10
    release["concordance"]["s01"] = {"intersection": many - 1, "perl": many, "rust": many,
                                     "adj_intersection": many - 1, "adj_perl": many - 1, "adj_rust": many - 1}
    scored = dataset_record()
    scored["concordance"]["s01"] = {"intersection": many - 1, "perl": many, "engine_tuples": many}
    document = values_of(tmp_path, record_set(release=release, fastvep=scored))
    row = rows_of(document, "a")[("Genome-wide ClinVar", "ClinVar GRCh37")]
    bar, marker, theirs = row["values"]
    assert (bar["tool"], bar["kind"], bar["value"], bar["counts"]) == ("vep-rs", "bar", "1.000000000", [many - 1, many, many])
    assert (marker["kind"], marker["value"], marker["counts"]) == ("marker", "1.000000000", [many - 1] * 3)
    # A comparator's row prints what the rule gives.
    assert (theirs["tool"], theirs["value"], theirs["counts"]) == ("fastvep", "1.000000000", [many - 1, many, many])
    # A mark's place is its counts', never its printed figure's: the bars stand 1e-10 short of 1 and are not
    # full, and the marker, whose counts are equal, stands at the floor.
    assert figure.discordance(bar) == figure.discordance(theirs) == Fraction(1, many) and figure.discordance(marker) == 0
    assert figure.marker_shape(marker) == "exact"
    floor = figure.discord_floor(row["values"])
    assert figure.discord(bar, floor) == pytest.approx(1e-10) and figure.discord(bar, floor) > floor == figure.discord(marker, floor)


@pytest.mark.parametrize("flags", [(), ("--draft",)], ids=["final", "draft"])
def test_adjusted_counts_of_the_released_tool_that_are_unequal_and_round_to_one_are_refused_by_name(tmp_path, flags):
    # The README's renderer refuses such a row of the released version: the figure would print an exact 1 for it.
    release = release_record()
    many = 10**10
    release["concordance"]["s01"].update(adj_intersection=many - 1, adj_perl=many, adj_rust=many)
    assert figure.f1(many - 1, many, many, "x") == "1.000000000"
    assert refusal(tmp_path, record_set(release=release), *flags) == (
        f"2001-01-01-release-{RELEASE}.json concordance.s01 adjusted counts are unequal and their F1 rounds to 1 at nine places "
        f"(reference {many}, vep-rs {many}, matched {many - 1}); the released version's adjusted F1 is exactly 1 or the record is "
        "not a shipped state")
    # Equal counts are exactly 1, and unequal counts that do not round to 1 print what they give.
    assert figure.released_adjusted_f1(5, 5, 5, "x", "vep-rs") == "1.000000000"
    assert figure.released_adjusted_f1(3, 4, 4, "x", "vep-rs") == "0.750000000"
    assert figure.released_adjusted_f1(many - 20, many, many, "x", "vep-rs") == "0.999999998"


def test_a_comparators_adjusted_marker_is_drawn_from_its_records_and_left_out_where_they_hold_none(tmp_path):
    scored, population = f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json", f"2001-01-01-population-{RELEASE}.json"
    document = values_of(tmp_path)
    rows = rows_of(document, "a")
    theirs = {label: next(v for v in row["values"] if (v["tool"], v["kind"]) == ("fastvep", "marker")) for (_, label), row in rows.items()}
    # Every row of panel (a): a ClinVar or structural-variant cell from the `adjusted` block of the comparator's
    # dataset record, a population dataset from its counts less the excluded ones, the pool from their sums.
    assert list(theirs) == [label for _, label in rows]
    cell = theirs["ClinVar GRCh37"]
    assert (cell["value"], cell["counts"], cell["record"]) == ("0.250000000", [1, 4, 4], scored)
    assert cell["key"] == "concordance.s01.adjusted: intersection, perl, engine_tuples"
    assert theirs["SV per-VCF, 2 files pooled GRCh38"]["key"] == "concordance.s08.adjusted: intersection, perl, engine_tuples"
    dataset = theirs["Population B GRCh37"]
    assert (dataset["value"], dataset["counts"], dataset["record"]) == (nine_places(6, 7, 7), [6, 7, 7], population)
    assert dataset["key"] == "comparator_concordance.per_dataset[dataset=pop_b]: intersection, perl_tuples - excluded_perl, engine_tuples - excluded_rust"
    assert (theirs["All four"]["value"], theirs["All four"]["counts"]) == (nine_places(58, 115, 115), [58, 115, 115])
    assert nine_places(6, 7, 7) == "0.857142857" and nine_places(58, 115, 115) == "0.504347826"
    # On the plugin rows a comparator's marker is read from the adjusted counts of its entries, summed as its bar's are.
    loeuf = next(v for v in rows_of(document, "c")[(REAL_GROUP, "LOEUF GRCh37")]["values"] if (v["tool"], v["kind"]) == ("fastvep", "marker"))
    counts = plugin_sums("real_data", "LOEUF", "GRCh37", "fastvep", "feature", ("adj_intersection", "adj_perl", "adj_rust"))
    assert (loeuf["value"], loeuf["counts"], loeuf["record"]) == (nine_places(*counts), counts, f"2001-01-01-plugins-{RELEASE}.json")
    assert loeuf["key"] == ("real_data[plugin=LOEUF,assembly=GRCh37,dataset=clinvar_grch37,engine=fastvep,projection=feature]: "
                            "adj_intersection, adj_perl, adj_rust")
    # The marker and the dashed line from its bar's end are drawn as the released tool's are, in the comparator's colour.
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, rows[("Whole-genome population datasets", "Population B GRCh37")])], LANES, FILLS, 0.3, 1e-7, 35.0)
    joins = joins_of(axes)
    markers = [(args, style) for args, style in axes.drawn("plot") if style.get("marker")]
    assert [(args[1], style["color"]) for args, style in joins] == [([-0.22, -0.22], FILLS["vep-rs"]), ([0.22, 0.22], FILLS["fastvep"])]
    assert joins[1][0][0] == drawn_join(0.25, 1 / 7, "adjusted", 35.0) and joins[1][0][0][0] == pytest.approx(0.25)
    # Its body stands three pixels past its bar's end: a solid stub, as the released tool's line of that length is.
    assert joins[1][1]["linestyle"] == figure.join_dashes(35.0 * math.log10(0.25 / figure.join_drawn(0.25, 1 / 7, "adjusted", 35.0)[1])) == "solid"
    assert [(args[1], style["marker"], style["markerfacecolor"]) for args, style in markers] == [
        ([-0.22], ">", FILLS["vep-rs"]), ([0.22], "D", FILLS["fastvep"])]
    # The caption says whose a marker is. One entry of the fixture, the comparator's per-feature AlphaMissense on
    # GRCh37, holds no adjusted counts: that marker is left out with a printed note and the caption's one line,
    # which a record that holds them all does not carry.
    left_out = "fastVEP's adjusted F1 is left out where no record holds its adjusted counts."
    caption = document["caption"]
    assert "(a) Raw F1 per dataset on a log(1 - F1) axis; a marker is its tool's adjusted F1, joined to the end of its bar by a line where the two differ." in caption
    assert "The bar is each tool's per-feature raw F1, the marker its tool's adjusted F1, joined to the end of its bar" in caption
    assert caption.count(left_out) == 1
    done, _ = run(tmp_path / "noted", record_set())
    assert done.stdout.count("note: no record holds fastVEP's adjusted counts on AlphaMissense GRCh37; the figure omits that marker\n") == 1
    whole = copy.deepcopy(PLUGINS)
    for entry in whole["real_data"]:
        if entry["engine"] == "fastvep" and "adj_perl" not in entry:
            entry.update(adj_intersection=entry["intersection"], adj_perl=entry["perl"], adj_rust=entry["rust"])
    document = values_of(tmp_path / "whole", record_set(plugins=whole))
    assert left_out not in document["caption"]
    scored = [row for row in rows_of(document, "c").values() if ("fastvep", "bar") in marks(row)]
    assert len(scored) == 16 and all(("fastvep", "marker") in marks(row) for row in scored)

    # A cell without its `adjusted` block and a dataset without its excluded counts: the marker is left out on
    # those rows and on the pool, each with a printed note, and the caption says so in one line.
    partial, counted = dataset_record(), population_record()
    del partial["concordance"]["s07"]["adjusted"]
    del counted["comparator_concordance"]["per_dataset"][0]["excluded_rust"]
    done, values = run(tmp_path / "partial", record_set(fastvep=partial, population=counted))
    assert done.returncode == 0, done.stderr
    for label in ("SV per-VCF, 2 files pooled GRCh37", "Population A GRCh38", "All four"):
        assert f"note: no record holds fastVEP's adjusted counts on {label}; the figure omits that marker\n" in done.stdout
    document = json.loads(values.read_text(encoding="utf-8"))
    held = [label for (_, label), row in rows_of(document, "a").items() if ("fastvep", "marker") in marks(row)]
    assert held == ["ClinVar GRCh37", "ClinVar GRCh38", "Population B GRCh37", "SV per-VCF, 2 files pooled GRCh38"]
    assert all(("fastvep", "bar") in marks(row) for row in rows_of(document, "a").values())
    assert document["caption"].count(left_out) == 1

    # Records that hold no adjusted counts of the comparator at all: its bars alone, the released tool's markers
    # named as its own, and the one line.
    bare, uncounted, unadjusted = dataset_record(), population_record(), copy.deepcopy(PLUGINS)
    for held_cell in bare["concordance"].values():
        del held_cell["adjusted"]
    for row in uncounted["comparator_concordance"]["per_dataset"]:
        del row["excluded_perl"], row["excluded_rust"]
    for name in ("real_data", "synthetic_fixtures", "custom_annotation"):
        for entry in unadjusted[name]:
            if entry["engine"] == "fastvep":
                for key in ("adj_intersection", "adj_perl", "adj_rust"):
                    entry.pop(key, None)
    document = values_of(tmp_path / "bare", record_set(fastvep=bare, population=uncounted, plugins=unadjusted))
    assert not any(("fastvep", "marker") in marks(row) for panel_id in ("a", "c") for row in rows_of(document, panel_id).values())
    assert all(("fastvep", "bar") in marks(row) or ("fastvep", "text") in marks(row) for row in rows_of(document, "c").values())
    assert "the marker is vep-rs's adjusted F1, joined to the end of its bar by a line where the two differ." in document["caption"]
    assert "The bar is each tool's per-feature raw F1, the marker vep-rs's adjusted F1, joined to the end of its bar" in document["caption"]
    assert document["caption"].count(left_out) == 1
    legend = {entry["mark"]: entry["label"] for entry in figure.legend_entries(document)}
    assert legend["exact"] == "vep-rs adjusted F1, exactly 1"


def test_the_comparator_carries_its_adjusted_marker_joined_to_its_bar_by_the_dashed_line_on_both_f1_panels(tmp_path):
    document = values_of(tmp_path)
    # A dataset row of panel (a) and a plugin row of panel (c): on each the comparator has its bar, its marker
    # and the dashed line between them, on its own lane, in its own colour.
    for panel_id, key in (("a", ("Genome-wide ClinVar", "ClinVar GRCh37")), ("c", (REAL_GROUP, "LOEUF GRCh37")),
                          ("c", ("Synthetic fixtures (data that cannot be redistributed)", "dbNSFP GRCh37 synthetic"))):
        row = rows_of(document, panel_id)[key]
        found = {(value["tool"], value["kind"]): value for value in row["values"]}
        bar, marker = found[("fastvep", "bar")], found[("fastvep", "marker")]
        axes = Recorder()
        figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, 35.0)
        joins = [args for args, style in joins_of(axes) if style["color"] == FILLS["fastvep"]]
        (drawn,) = [(args, style) for args, style in axes.drawn("plot") if style.get("marker") and style["markerfacecolor"] == FILLS["fastvep"]]
        # The marker is one of the released tool's two symbols, by the comparator's own adjusted counts.
        shape = "exact" if len(set(marker["counts"])) == 1 else "adjusted"
        assert drawn[0] == ([pytest.approx(figure.discord(marker, 1e-7))], [0.22]) and drawn[1]["marker"] == figure.MARKERS[shape][0], key
        assert drawn[1]["markeredgecolor"] == figure.GROUND and drawn[1]["zorder"] > figure.BAR_ORDER
        # The dashed line runs from the end of the comparator's bar to its marker's body where the two do not touch.
        tip = figure.discord(bar, 1e-7)
        end = figure.join_end(tip, figure.discord(marker, 1e-7), shape, 35.0)
        assert joins == ([] if end is None else [(drawn_join(tip, figure.discord(marker, 1e-7), shape, 35.0), [0.22, 0.22])]), key
    # The arrow at exactly 1 and the diamond below 1, as the released tool's: by the counts, never by the tool.
    rows = [(0, f1_row(("vep-rs", "bar", "0.990000000"), ("vep-rs", "marker", "1.000000000"), ("fastvep", "bar", "0.900000000"), ("fastvep", "marker", "1.000000000"))),
            (1, f1_row(("vep-rs", "bar", "0.990000000"), ("vep-rs", "marker", "0.999000000"), ("fastvep", "bar", "0.100000000"), ("fastvep", "marker", "0.900000000")))]
    axes = Recorder()
    figure.draw_f1_marks(axes, rows, LANES, FILLS, 0.3, 1e-7, 35.0)
    assert [(args[1], style["marker"], style["markerfacecolor"]) for args, style in axes.drawn("plot") if style.get("marker")] == [
        ([-0.22], ">", FILLS["vep-rs"]), ([0.22], ">", FILLS["fastvep"]), ([1 - 0.22], "D", FILLS["vep-rs"]), ([1 + 0.22], "D", FILLS["fastvep"])]
    assert len(joins_of(axes)) == 4
    # Every plugin row on which the comparator has a figure carries its marker where the row's entries hold the
    # adjusted counts, and no row carries a marker of the comparator without its bar.
    lists = {REAL_GROUP: "real_data", "Committed table": "committed_table",
             "Synthetic fixtures (data that cannot be redistributed)": "synthetic_fixtures", "Custom annotation (--custom sources)": "custom_annotation"}
    scored = 0
    for (group, label), row in rows_of(document, "c").items():
        plugin, assembly = label.removesuffix(" synthetic").rsplit(" ", 1)
        held = plugin_sums(lists[group], plugin, assembly, "fastvep", "feature", ("adj_intersection", "adj_perl", "adj_rust"))
        assert (("fastvep", "marker") in marks(row)) == (held is not None), label
        assert ("fastvep", "bar") in marks(row) or ("fastvep", "marker") not in marks(row), label
        scored += ("fastvep", "bar") in marks(row)
    assert scored == 16 and sum(("fastvep", "marker") in marks(row) for row in rows_of(document, "c").values()) == 15


def test_the_caption_says_what_a_comparators_adjusted_figure_applies_only_in_the_words_of_a_record(tmp_path):
    name = f"2001-01-01-population-{RELEASE}.json"
    # No key of these records states what a comparator's adjusted figure applies: the caption says nothing on it,
    # in no words of the script's own.
    document = values_of(tmp_path)
    assert any(("fastvep", "marker") in marks(row) for row in rows_of(document, "a").values())
    assert not any(word in document["caption"] for word in ("adjusted F1 applies", "not modelled", "record's note reads"))
    assert not any(source["key"] == "comparator_concordance.note" for source in document["sources"])
    assert not hasattr(figure, "ADJUSTED_RULES")
    # The population record's note on its comparator block states it: the caption quotes the sentence of the note
    # that speaks of the adjusted figure, once, and the values file names the key it is drawn from.
    population = population_record()
    population["comparator_concordance"]["note"] = ("The comparator's output on every dataset, scored against the reference's archived output. "
                                                    "The adjusted figure sets aside the same documented classes for every tool; where a tool "
                                                    "writes none of their shapes it equals the raw one. A residual is reported, not gated.")
    sentence = "The adjusted figure sets aside the same documented classes for every tool; where a tool writes none of their shapes it equals the raw one."
    document = values_of(tmp_path / "noted", record_set(population=population))
    quoted = f'On a comparator\'s adjusted F1 the record\'s note reads: "{sentence}"'
    assert document["caption"].count(quoted) == 1 and document["caption"].count(sentence) == 1
    assert {"text": sentence, "record": name, "key": "comparator_concordance.note"} in document["sources"]
    assert sentence in population["comparator_concordance"]["note"]
    # A note that says nothing of the adjusted figure, and one that is no text, are quoted nowhere.
    for index, note in enumerate(("The comparator's output on every dataset. A residual is reported, not gated.", {"adjusted": "by class"})):
        population["comparator_concordance"]["note"] = note
        document = values_of(tmp_path / f"other{index}", record_set(population=population))
        assert "record's note reads" not in document["caption"] and not any(s["key"] == "comparator_concordance.note" for s in document["sources"])
    # A comparator without an adjusted marker anywhere has no sentence on its adjusted figure, whatever the note.
    population["comparator_concordance"]["note"] = "The adjusted figure sets aside the same documented classes for every tool."
    registry = registry_with(tmp_path, lambda rows: rows["fastvep"].update(concordance_panel="no"))
    document = values_of(tmp_path / "unmarked", record_set(population=population), registry=registry)
    assert "record's note reads" not in document["caption"]


def test_a_draft_quotes_no_note_of_a_record_on_a_comparators_adjusted_figure(tmp_path):
    sentence = "The adjusted figure sets aside the same documented classes for every tool; where a tool writes none of their shapes it equals the raw one."
    population = population_record()
    population["comparator_concordance"]["note"] = f"The comparator's output on every dataset. {sentence} A residual is reported, not gated."

    def noted(document: dict) -> tuple[int, int, int]:
        """How often the caption quotes a note, how often it holds the sentence, and the sources that name the note's key."""
        return (document["caption"].count("record's note reads"), document["caption"].count(sentence),
                sum(source["key"] == "comparator_concordance.note" for source in document["sources"]))

    # Records of one release that passed their gate: the caption quotes the note's sentence once, with its key.
    final = values_of(tmp_path / "final", record_set(population=population))
    assert noted(final) == (1, 1, 1) and not final["caption"].startswith(figure.DRAFT_LINE)
    # A draft quotes none, whatever makes it a draft: the same records under --draft, an ungated plugin record and
    # records of two releases. The sentence is in neither text and no source names the note's key.
    mixed = mixed_releases()
    next(record for name, record in mixed.items() if "-population-" in name)["comparator_concordance"]["note"] = population["comparator_concordance"]["note"]
    for name, records in (("same", record_set(population=population)), ("ungated", record_set(population=population, plugins=ungated())), ("mixed", mixed)):
        draft = values_of(tmp_path / name, records, "--draft")
        assert draft["caption"].startswith(figure.DRAFT_LINE) and noted(draft) == (0, 0, 0), name
        assert sentence not in draft["alt_text"] and "sets aside" not in draft["caption"], name
        # The comparator's adjusted markers are drawn all the same, and the other sentences on them stand.
        assert any(("fastvep", "marker") in marks(row) for row in rows_of(draft, "a").values()), name
        assert "a marker is its tool's adjusted F1" in draft["caption"], name


def test_a_comparators_adjusted_marker_is_left_out_where_the_record_states_that_it_is_not_determined(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"
    left_out = "fastVEP's adjusted F1 is left out where no record holds its adjusted counts."

    def stated(*counts: int) -> dict:
        """The fixture with the comparator's two per-feature AlphaMissense entries stating, GRCh37 then GRCh38, how
        many of its values meet the condition of a declared class that is not applied."""
        plugins = copy.deepcopy(PLUGINS)
        entries = [e for e in plugins["real_data"] if (e["plugin"], e["engine"], e["projection"]) == ("AlphaMissense", "fastvep", "feature")]
        assert [e["assembly"] for e in entries] == ["GRCh37", "GRCh38"] and ["adj_perl" in e for e in entries] == [False, True]
        for entry, met in zip(entries, counts):
            entry["adjusted_undetermined"] = {"class": "a_declared_class", "values_meeting_condition": met}
        return plugins

    key = "real_data[plugin=AlphaMissense,assembly={},dataset=clinvar_{},engine=fastvep,projection=feature]: adjusted_undetermined.values_meeting_condition"
    # Three values meet the condition on GRCh37, none on GRCh38: the first row's marker is left out and the caption
    # says in one line that the adjusted F1 is not determined there; the second row's marker is drawn from the
    # entry's adjusted counts, since a class that sets nothing aside leaves the figure as it is.
    done, values = run(tmp_path / "one", record_set(plugins=stated(3, 0)))
    assert done.returncode == 0, done.stderr
    document = json.loads(values.read_text(encoding="utf-8"))
    rows = rows_of(document, "c")
    first, second = (marks(rows[(REAL_GROUP, f"AlphaMissense {assembly}")]) for assembly in ("GRCh37", "GRCh38"))
    assert ("fastvep", "bar") in first and ("fastvep", "marker") not in first
    counts = plugin_sums("real_data", "AlphaMissense", "GRCh38", "fastvep", "feature", ("adj_intersection", "adj_perl", "adj_rust"))
    assert second[("fastvep", "marker")] == nine_places(*counts)
    assert done.stdout.count(f"note: {name} states that fastVEP's adjusted F1 is not determined on AlphaMissense GRCh37; the figure omits that marker\n") == 1
    line = "fastVEP's adjusted F1 is not determined on AlphaMissense GRCh37."
    assert document["caption"].count(line) == 1 and document["caption"].count("is not determined") == 1 and left_out not in document["caption"]
    assert {"text": line, "record": name, "key": key.format("GRCh37", "grch37")} in document["sources"]
    # One or more on both rows: both markers are left out, the entry's adjusted counts notwithstanding, and the
    # caption still carries one line, which names both rows and both record keys.
    document = values_of(tmp_path / "both", record_set(plugins=stated(1, 2)))
    rows = rows_of(document, "c")
    assert not any(("fastvep", "marker") in marks(rows[(REAL_GROUP, f"AlphaMissense {assembly}")]) for assembly in ("GRCh37", "GRCh38"))
    line = "fastVEP's adjusted F1 is not determined on AlphaMissense GRCh37 and AlphaMissense GRCh38."
    assert document["caption"].count(line) == 1 and document["caption"].count("is not determined") == 1
    assert {"text": line, "record": name, "key": f"{key.format('GRCh37', 'grch37')}; {key.format('GRCh38', 'grch38')}"} in document["sources"]
    assert "its marker on 14 of them" in document["alt_text"]
    # None on either row: nothing is undetermined. The row whose entry holds adjusted counts carries its marker,
    # and the one whose entry holds none is left out as any row without adjusted counts is.
    document = values_of(tmp_path / "none", record_set(plugins=stated(0, 0)))
    rows = rows_of(document, "c")
    assert ("fastvep", "marker") not in marks(rows[(REAL_GROUP, "AlphaMissense GRCh37")])
    assert ("fastvep", "marker") in marks(rows[(REAL_GROUP, "AlphaMissense GRCh38")])
    assert "is not determined" not in document["caption"] and document["caption"].count(left_out) == 1
    # A statement without its count is refused by name.
    plugins = stated(3, 0)
    next(e for e in plugins["real_data"] if "adjusted_undetermined" in e)["adjusted_undetermined"] = {"class": "a_declared_class"}
    assert refusal(tmp_path / "uncounted", record_set(plugins=plugins)).endswith(".adjusted_undetermined carries no count under values_meeting_condition: None")


def test_the_alternative_text_ends_a_comparators_clause_cleanly_where_no_row_carries_a_cross(tmp_path):
    # Gate rows on which the comparator is scored, and no other: no row carries a cross, and one of its entries
    # holds no adjusted counts, so it has fewer markers than bars.
    plugins = copy.deepcopy(PLUGINS)
    scored = {(e["plugin"], e["assembly"]) for name in ("real_data", "synthetic_fixtures") for e in plugins[name] if e["engine"] == "fastvep"}
    plugins["rows"] = [row for row in plugins["rows"] if (row["plugin"], row["assembly"]) in scored]
    document = values_of(tmp_path, record_set(plugins=plugins))
    rows = list(rows_of(document, "c").values())
    bars, markers = (sum(("fastvep", kind) in marks(row) for row in rows) for kind in ("bar", "marker"))
    assert bars == len(rows) and markers == bars - 1 and not any(("fastvep", "text") in marks(row) for row in rows)
    alt = document["alt_text"]
    assert f"fastVEP has a figure on {bars} of the {len(rows)} rows, its marker on {markers} of them; vep-rs's adjusted F1 is exactly 1 on" in alt
    assert ",;" not in alt and "a cross" not in alt
    # With crosses the clause on markers keeps its comma before the clause on them.
    assert "its marker on 15 of them, and a cross (no such source) on 43; " in values_of(tmp_path / "crossed")["alt_text"]


def test_a_row_whose_adjusted_counts_are_all_zero_has_no_marker_and_the_caption_says_why(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"
    left_out = "fastVEP's adjusted F1 is left out where no record holds its adjusted counts."
    # Every tuple of the comparator's LOEUF row and of the released tool's REVEL row on GRCh37 is set aside.
    plugins = copy.deepcopy(PLUGINS)
    # The released tool's row sums its entry on the population dataset too, and that one is set aside whole as well.
    for engine, plugin, summed in (("fastvep", "LOEUF", 1), ("vep-rs", "REVEL", 2)):
        entries = [e for name in ("real_data", "population_data") for e in plugins[name]
                   if (e["engine"], e["plugin"], e["assembly"], e["projection"]) == (engine, plugin, "GRCh37", "feature")]
        assert len(entries) == summed
        for entry in entries:
            entry.update(adj_intersection=0, adj_perl=0, adj_rust=0)
    done, values = run(tmp_path, record_set(plugins=plugins))
    assert done.returncode == 0, done.stderr
    document = json.loads(values.read_text(encoding="utf-8"))
    rows = rows_of(document, "c")
    # Each row keeps its bar and has no adjusted marker: no F1 is left to place.
    loeuf, revel = marks(rows[(REAL_GROUP, "LOEUF GRCh37")]), marks(rows[(REAL_GROUP, "REVEL GRCh37")])
    assert ("fastvep", "bar") in loeuf and ("fastvep", "marker") not in loeuf and ("vep-rs", "bar") in revel and ("vep-rs", "marker") not in revel
    for tool, label in (("fastVEP", "LOEUF GRCh37"), ("vep-rs", "REVEL GRCh37")):
        assert done.stdout.count(f"note: {name} sets aside every tuple of {tool} on {label}; the figure omits that marker\n") == 1
        line = f"{tool}'s adjusted F1 is left out on {label}, where every tuple of the row is set aside and none is left to score."
        assert document["caption"].count(line) == 1
        (source,) = [source for source in document["sources"] if source["text"] == line]
        ending = "projection=feature]: adj_intersection, adj_perl, adj_rust" + (", summed over 2 entries" if tool == "vep-rs" else "")
        assert source["record"] == name and source["key"].endswith(ending)
    # The line on adjusted counts that no record holds speaks of the one entry without them and of no other.
    assert document["caption"].count(left_out) == 1
    assert "no record holds fastVEP's adjusted counts on LOEUF GRCh37" not in done.stdout


def test_note_line_is_in_the_values_file_and_among_the_drawn_strings(tmp_path):
    document = values_of(tmp_path)
    supported = next(block["of"] for block in PLUGINS["support"] if block["engine"] == "vep-rs")
    note = (f"Counts are of the {supported} plugins vep-rs supports, of Ensembl VEP's {PLUGINS['not_ported']['of']}. "
            "Plugins vep-rs does not port are not counted; the README lists them by reason.")
    assert document["notes"] == [] and panel(document, "d")["note"] == note
    drawn = figure.drawn_texts(document, figure.read_registry(REGISTRY))
    assert drawn.count(note) == 1
    assert {"(a)", "(b)", "(c)", "(d)", "fastVEP: no such source", "Synthetic fixtures (data that cannot be redistributed)"} <= set(drawn)


def test_the_alternative_text_alone_says_that_the_plugin_panel_continues_in_the_second_column(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"
    document = values_of(tmp_path)
    by_id = {p["id"]: p for p in document["panels"]}
    groups = by_id["c"]["groups"]
    # The split is a function of the panels alone, and the drawing splits the panel by it.
    first = figure.plugin_split(by_id)
    assert 0 < first < len(groups) and "first = plugin_split(panels)" in SCRIPT.read_text(encoding="utf-8")
    assert figure.plugin_split({key: by_id[key] for key in ("a", "b")}) == 0
    assert figure.plugin_split({**by_id, "c": {**by_id["c"], "groups": groups[:1]}}) == 1
    # The alternative text names the group from which the panel continues, once, with its key among the sources;
    # the caption says nothing of columns: the image's own title says that the second part continues the first.
    clause = f"The plugin panel runs over both columns: its groups from {groups[first]['label']} on continue in the second column, under the same title."
    assert document["alt_text"].count(clause) == 1 and document["alt_text"].count("second column") == 1
    assert "column" not in document["caption"]
    (source,) = [source for source in document["sources"] if source["text"] == clause]
    assert source["record"] == name and source["key"].startswith("rows: scored_on")
    # A figure without a plugin panel has no such sentence.
    assert "second column" not in values_of(tmp_path / "two", record_set(plugins=None))["alt_text"]


def test_the_key_of_the_sentence_on_the_plugin_panels_columns_names_the_order_the_figure_gives_its_groups(tmp_path):
    document = values_of(tmp_path)
    (source,) = [source for source in document["sources"] if source["text"].startswith("The plugin panel runs over both columns:")]
    assert source["key"] == figure.CONTINUES_KEY == (
        "rows: scored_on (the plugin panel's groups in the order the figure gives them, which the layout splits over the two columns)")
    assert "record's order" not in source["key"]
    # The order is the figure's own, one group per kind of data in the order of PLUGIN_GROUPS and the custom sources
    # last, whatever order the record's rows stand in.
    kinds = {kind: label for kind, _, label, _ in figure.PLUGIN_GROUPS} | {figure.REAL_KIND: REAL_GROUP}
    ordered = [kinds[kind] for kind in kinds if any(row["scored_on"] == kind for row in PLUGINS["rows"])] + [figure.CUSTOM_GROUP[1]]
    assert [group["label"] for group in panel(document, "c")["groups"]] == ordered
    assert ordered[:3] == [REAL_GROUP, "Committed table", "Synthetic fixtures (data that cannot be redistributed)"]
    turned = copy.deepcopy(PLUGINS)
    turned["rows"].reverse()
    assert list(dict.fromkeys(row["scored_on"] for row in turned["rows"])) != list(dict.fromkeys(row["scored_on"] for row in PLUGINS["rows"]))
    again = values_of(tmp_path / "turned", record_set(plugins=turned))
    assert [group["label"] for group in panel(again, "c")["groups"]] == ordered
    assert [source["key"] for source in again["sources"] if source["text"].startswith("The plugin panel runs over both columns:")] == [figure.CONTINUES_KEY]


def record_sources(document: dict) -> list[dict]:
    """The sources of a values file that name one key of the plugin record's references or of its list of plugins not ported."""
    return [source for source in document["sources"] if source["key"].startswith(("perl_reference", "not_ported"))]


def test_the_texts_say_from_the_records_own_counts_what_the_image_leaves_unsaid(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"
    document = values_of(tmp_path)
    texts = (document["caption"], document["alt_text"])
    files = {record["file"] for record in document["records"]}
    lists = {"clinvar": "real_data", "synthetic_fixture": "synthetic_fixtures", "committed_table": "committed_table"}
    # The rows on which the two tools' figures are summed over different numbers of the row's sets, counted here
    # from the fixture's entries: both texts name them with each tool's number, once, and the values file names the key.
    uneven = []
    for row in PLUGINS["rows"]:
        ours, theirs = (sum((e["plugin"], e["assembly"], e["engine"], e["projection"]) == (row["plugin"], row["assembly"], engine, "feature")
                            for e in PLUGINS[lists[row["scored_on"]]]) for engine in ("vep-rs", "fastvep"))
        if theirs and ours != theirs:
            uneven.append(f"{row['plugin']} {row['assembly']} (vep-rs {ours}, fastVEP {theirs})")
    assert uneven
    sets = (f"A row's figure is summed over the sets of the row that the record holds for the tool, and the tools differ in them on "
            f"{len(uneven)} rows: {figure.and_list(uneven)}.")
    # The row names of the plugin panel against its plugins: a plugin scored under two names counts as one.
    names = list(dict.fromkeys(row["plugin"] for row in PLUGINS["rows"]))
    held = {}
    for row_name in names:
        held.setdefault(row_name.split(" (")[0], []).append(row_name)
    doubled = [f"{plugin} is scored under {figure.and_list(under)}" for plugin, under in held.items() if len(under) > 1]
    assert doubled and len(names) == len(held) + len(doubled)
    named = f"The plugin panel's {len(PLUGINS['rows'])} gate rows carry {len(names)} row names for {len(held)} plugins: {figure.and_list(doubled)}."
    # Where an F1 of exactly 1 stands, with how many of the released tool's adjusted markers stand there.
    ours = [value for panel_id in ("a", "c") for row in rows_of(document, panel_id).values() for value in row["values"]
            if (value["tool"], value["kind"]) == ("vep-rs", "marker")]
    exact = ("An F1 of exactly 1 has no place on a scale of 1 - F1: it is drawn where the axis reads 1, past a break in the axis line, and "
             f"{sum(len(set(value['counts'])) == 1 for value in ours)} of the {len(ours)} adjusted markers of vep-rs stand there.")
    for sentence in (sets, named, exact):
        assert [text.count(sentence) for text in texts] == [1, 1], sentence
        (source,) = [source for source in document["sources"] if source["text"] == sentence]
        assert source["key"] and all(file in files for file in source["record"].split(" and "))
    assert next(source for source in document["sources"] if source["text"] == named) == {"text": named, "record": name, "key": "rows: plugin"}
    # The panel of counts states each tool's count in both texts, as its values hold it.
    assert f"a bar each: {READS}." in document["alt_text"]
    assert f"supports: vep-rs {PLUGINS['support'][0]['read']} and fastVEP {PLUGINS['support'][1]['read']}; plugins" in document["caption"]
    # A sentence whose record holds nothing to say is left out: without a plugin record the texts speak of no row
    # names and no sets, and with the comparator kept out of the concordance panels no row differs in its sets.
    two = values_of(tmp_path / "two", record_set(plugins=None))
    assert not any(word in text for text in (two["caption"], two["alt_text"]) for word in ("row names", "differ in them"))
    assert two["caption"].count("it is drawn where the axis reads 1, past a break in the axis line") == 1
    alone = values_of(tmp_path / "alone", registry=registry_with(tmp_path, lambda rows: rows["fastvep"].update(concordance_panel="no")))
    assert "differ in them" not in alone["caption"] and alone["caption"].count(named) == 1


def test_the_sentence_on_differing_sets_reads_one_row_where_exactly_one_row_differs(tmp_path):
    lead = "A row's figure is summed over the sets of the row that the record holds for the tool, and the tools differ in them on "

    def feature_sets(plugins: dict, plugin: str, assembly: str, engine: str) -> list[dict]:
        return [e for e in plugins["synthetic_fixtures"] if (e["plugin"], e["assembly"], e["engine"], e["projection"]) == (plugin, assembly, engine, "feature")]

    # The fixture's rows on which the released tool is summed over two sets and the comparator over one.
    uneven = [(row["plugin"], row["assembly"]) for row in PLUGINS["rows"] if row["scored_on"] == "synthetic_fixture"
              and len(feature_sets(PLUGINS, row["plugin"], row["assembly"], "fastvep")) == 1
              and len(feature_sets(PLUGINS, row["plugin"], row["assembly"], "vep-rs")) == 2]
    assert len(uneven) == 3
    several = values_of(tmp_path / "several")
    listed = [f"{plugin} {assembly} (vep-rs 2, fastVEP 1)" for plugin, assembly in uneven]
    assert f"{lead}3 rows: {figure.and_list(listed)}." in several["caption"]
    # All but the first brought to one set of the released tool: one row differs, and both texts say "1 row".
    plugins = copy.deepcopy(PLUGINS)
    for plugin, assembly in uneven[1:]:
        plugins["synthetic_fixtures"].remove(feature_sets(plugins, plugin, assembly, "vep-rs")[1])
    one = values_of(tmp_path / "one", record_set(plugins=plugins))
    sentence = f"{lead}1 row: {listed[0]}."
    for text in (one["caption"], one["alt_text"]):
        assert text.count(sentence) == 1 and "1 rows" not in text and text.count("differ in them") == 1
    assert [source["key"] for source in one["sources"] if source["text"] == sentence] == [figure.UNEVEN_KEY]


def test_the_sentence_on_row_names_reads_one_plugin_where_the_rows_are_of_one_plugin(tmp_path):
    # The fixture's plugin that is scored under two row names, and a record whose gate rows are that plugin's alone.
    held: dict[str, list[str]] = {}
    for row in PLUGINS["rows"]:
        under = held.setdefault(row["plugin"].split(" (")[0], [])
        if row["plugin"] not in under:
            under.append(row["plugin"])
    ((plugin, names),) = [(plugin, names) for plugin, names in held.items() if len(names) > 1]
    several = values_of(tmp_path / "several")
    assert f"row names for {len(held)} plugins: {plugin} is scored under {figure.and_list(names)}." in several["caption"] and len(held) > 1
    plugins = copy.deepcopy(PLUGINS)
    plugins["rows"] = [row for row in plugins["rows"] if row["plugin"] in names]
    one = values_of(tmp_path / "one", record_set(plugins=plugins))
    sentence = (f"The plugin panel's {len(plugins['rows'])} gate rows carry {len(names)} row names for 1 plugin: "
                f"{plugin} is scored under {figure.and_list(names)}.")
    assert len(plugins["rows"]) >= len(names) == 2
    for text in (one["caption"], one["alt_text"]):
        assert text.count(sentence) == 1 and "1 plugins" not in text and text.count("row names for") == 1
    assert [source["key"] for source in one["sources"] if source["text"] == sentence] == ["rows: plugin"]


def test_the_sentence_on_an_f1_of_exactly_one_reads_stands_where_one_marker_stands_there(tmp_path):
    lead = "it is drawn where the axis reads 1, past a break in the axis line, and "

    def inexact(*cells: str, datasets: bool) -> dict:
        """The release record with one more adjusted reference tuple on each of `cells`, and on the population
        datasets and their pool when `datasets`: those adjusted F1 fall below 1."""
        release = release_record()
        for name in cells:
            release["concordance"][name]["adj_perl"] += 1
        if datasets:
            for block in (*release["population_concordance"]["per_dataset"], release["population_concordance"]["pooled"]):
                block["excluded_perl"] = 0
        return release

    def texts(directory: Path, release: dict) -> tuple[dict, list[dict]]:
        document = values_of(directory, record_set(plugins=None, release=release))
        ours = [value for row in rows_of(document, "a").values() for value in row["values"] if (value["tool"], value["kind"]) == ("vep-rs", "marker")]
        return document, ours

    # One of the released tool's seven adjusted markers stands at exactly 1: both texts read "stands".
    document, ours = texts(tmp_path / "one", inexact("s01", "s02", "s07", datasets=True))
    assert len(ours) == 7 and sum(len(set(value["counts"])) == 1 for value in ours) == 1
    sentence = f"{lead}1 of the 7 adjusted markers of vep-rs stands there."
    for text in (document["caption"], document["alt_text"]):
        assert text.count(sentence) == 1 and "stand there" not in text
    assert [source["key"] for source in document["sources"] if source["text"].endswith(sentence)] == [figure.EXACT_KEY.format(released="vep-rs")]
    # None, two and all seven read "stand".
    for name, release, exact in (("none", inexact("s01", "s02", "s07", "s08", datasets=True), 0), ("two", inexact("s01", "s02", datasets=True), 2),
                                 ("all", release_record(), 7)):
        document, ours = texts(tmp_path / name, release)
        assert sum(len(set(value["counts"])) == 1 for value in ours) == exact, name
        for text in (document["caption"], document["alt_text"]):
            assert text.count(f"{lead}{exact} of the 7 adjusted markers of vep-rs stand there.") == 1 and "stands there" not in text, name


def test_the_sentence_on_an_f1_of_exactly_one_reads_one_adjusted_marker_where_the_released_tool_has_one(tmp_path):
    lead = "An F1 of exactly 1 has no place on a scale of 1 - F1: it is drawn where the axis reads 1, past a break in the axis line, and "
    # The noun agrees with the number of the released tool's adjusted markers, the verb with how many stand there.
    assert figure.exact_sentence(1, 1, "vep-rs") == f"{lead}1 of the 1 adjusted marker of vep-rs stands there."
    assert figure.exact_sentence(0, 1, "vep-rs") == f"{lead}0 of the 1 adjusted marker of vep-rs stand there."
    assert figure.exact_sentence(1, 7, "vep-rs") == f"{lead}1 of the 7 adjusted markers of vep-rs stands there."
    assert figure.exact_sentence(7, 7, "vep-rs") == f"{lead}7 of the 7 adjusted markers of vep-rs stand there."
    assert not any("1 adjusted markers" in figure.exact_sentence(exact, 1, "vep-rs") for exact in (0, 1))
    # Both texts carry that sentence, with the counts of their own values, and no other wording of it.
    document = values_of(tmp_path, record_set(plugins=None))
    ours = [value for row in rows_of(document, "a").values() for value in row["values"] if (value["tool"], value["kind"]) == ("vep-rs", "marker")]
    sentence = figure.exact_sentence(sum(len(set(value["counts"])) == 1 for value in ours), len(ours), "vep-rs")
    assert len(ours) == 7 and sentence.endswith("of the 7 adjusted markers of vep-rs stand there.")
    for text in (document["caption"], document["alt_text"]):
        assert text.count(sentence) == 1 and text.count("adjusted marker") == 1
    assert SCRIPT.read_text(encoding="utf-8").count("EXACT_SENTENCE.format(") == 1


def test_every_figure_and_sentence_of_the_texts_that_no_value_holds_names_the_record_key_it_is_drawn_from(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"
    document = values_of(tmp_path)
    stand_ins = "LoFTEE's reference ran with stand-ins for the MaxEntScan files its module bundles, on the Ensembl VEP image with samtools added."
    listed = str(PLUGINS["not_ported"]["of"])
    # The caption's sentence on the reference's stand-ins and the note's count of the plugins the reference lists are
    # no panel's value: the values file names the record and the key each is drawn from.
    assert record_sources(document) == [
        {"text": stand_ins, "record": name, "key": "perl_reference[3].stand_ins.plugin, perl_reference[3].derived_image.package"},
        {"text": listed, "record": name, "key": "not_ported.of"},
    ]
    reference = PLUGINS["perl_reference"][3]
    assert (reference["stand_ins"]["plugin"], reference["derived_image"]["package"]) == ("LoFTEE", "samtools")
    files = {record["file"] for record in document["records"]}
    note = panel(document, "d")["note"]
    assert all(source["key"] and all(file in files for file in source["record"].split(" and ")) for source in document["sources"])
    assert document["caption"].count(stand_ins) == 1 and f"of Ensembl VEP's {listed}." in note
    # Every number of the note is a value of its panel or one of those sources.
    held = {number for row in panel(document, "d")["groups"][0]["rows"] for value in row["values"] for number in re.findall(r"\d+", value["value"])}
    assert set(re.findall(r"\d+", note)) <= held | {source["text"] for source in document["sources"]}
    # A record that holds no such key gives no such text: without stand-ins the sentence is not written, and
    # without the panel of counts the count of the reference's plugins is not stated.
    plain = copy.deepcopy(PLUGINS)
    del plain["perl_reference"][3]["stand_ins"]
    document = values_of(tmp_path / "plain", record_set(plugins=plain))
    assert "stand-ins" not in document["caption"] and record_sources(document) == [{"text": listed, "record": name, "key": "not_ported.of"}]
    uncounted = copy.deepcopy(PLUGINS)
    uncounted["support"] = []
    document = values_of(tmp_path / "uncounted", record_set(plugins=uncounted))
    assert [source["key"] for source in record_sources(document)] == ["perl_reference[3].stand_ins.plugin, perl_reference[3].derived_image.package"]
    assert listed not in re.findall(r"\d+", document["caption"])
    assert record_sources(values_of(tmp_path / "two", record_set(plugins=None))) == []


def test_every_number_of_sets_the_texts_print_stands_in_the_values_file_with_the_entries_it_counts(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"
    document = values_of(tmp_path)
    lists = {"clinvar": "real_data", "synthetic_fixture": "synthetic_fixtures", "committed_table": "committed_table"}
    engines = {"vep-rs": "vep-rs", "fastVEP": "fastvep"}
    selector = re.compile(r"^(\w+)\[plugin=([^,]+),assembly=([^,]+),(dataset|corpus)=([^,]+),engine=([^,]+),projection=feature\]: the number of entries$")
    counted = [source for source in document["sources"] if source["key"].endswith(": the number of entries")]
    # The rows whose tools are summed over different numbers of sets, and each tool's number, counted here from
    # the fixture's entries: one source a number, in the order the sentence prints them.
    expected = []
    for row in PLUGINS["rows"]:
        held = {tool: [e for e in PLUGINS[lists[row["scored_on"]]]
                       if (e["plugin"], e["assembly"], e["engine"], e["projection"]) == (row["plugin"], row["assembly"], engine, "feature")]
                for tool, engine in engines.items()}
        if held["fastVEP"] and len(held["vep-rs"]) != len(held["fastVEP"]):
            expected += [(row["plugin"], row["assembly"], tool, len(entries)) for tool, entries in held.items()]
    assert len(expected) == 6 and [source["text"] for source in counted] == [str(number) for *_, number in expected]
    (sentence,) = [source["text"] for source in document["sources"] if source["key"] == figure.UNEVEN_KEY]
    assert re.findall(r"(?:vep-rs|fastVEP) (\d+)", sentence) == [source["text"] for source in counted]
    for source, (plugin, assembly, tool, number) in zip(counted, expected):
        # The key names the entries the number counts: the record's list, the row and the tool's engine.
        list_name, key_plugin, key_assembly, target, target_value, engine = selector.match(source["key"]).groups()
        assert (key_plugin, key_assembly, engine, source["record"]) == (plugin, assembly, engines[tool], name)
        entries = [e for e in PLUGINS[list_name] if (e["plugin"], e["assembly"], e.get(target), e["engine"], e["projection"])
                   == (plugin, assembly, target_value, engine, "feature")]
        assert len(entries) == number == int(source["text"])
        assert f"{plugin} {assembly} (" in sentence and f"{tool} {number}" in sentence
    # The number of plugins the reference lists, which the note under the panel of counts prints, is a source too.
    listed = {"text": str(PLUGINS["not_ported"]["of"]), "record": name, "key": "not_ported.of"}
    assert document["sources"].count(listed) == 1 and f"of Ensembl VEP's {listed['text']}." in panel(document, "d")["note"]
    # Where no row's tools differ in their sets no number of sets is printed, and none is kept.
    alone = values_of(tmp_path / "alone", registry=registry_with(tmp_path, lambda rows: rows["fastvep"].update(concordance_panel="no")))
    assert "differ in them" not in alone["caption"] and not any(source["key"].endswith(": the number of entries") for source in alone["sources"])


def test_the_wall_time_panel_names_no_tool_and_ties_no_name_to_a_circle(tmp_path):
    document = values_of(tmp_path)
    registry = figure.read_registry(REGISTRY)
    circles = [value["tool"] for group in panel(document, "b")["groups"] for row in group["rows"] for value in row["values"]]
    assert {tool: circles.count(tool) for tool in ("vep-rs", "fastvep", "ensembl-vep")} == {"vep-rs": 8, "fastvep": 8, "ensembl-vep": 8}
    # Eight circles a tool, each in its tool's colour, and no tool's name among the strings the panel draws: a
    # name stands alone only as a row of the panel of counts, and nowhere without that panel.
    drawn = figure.drawn_texts(document, registry)
    assert [drawn.count(name) for name in ("vep-rs", "fastVEP", "Ensembl VEP")] == [1, 1, 0]
    assert [row["label"] for group in panel(document, "d")["groups"] for row in group["rows"]] == ["vep-rs", "fastVEP"]
    assert not {"vep-rs", "fastVEP", "Ensembl VEP"} & set(figure.drawn_texts(values_of(tmp_path / "two", record_set(plugins=None)), registry))
    # The legend names each tool once, with its version and its colour.
    named = [entry for entry in figure.legend_entries(document) if entry["tool"] and entry["mark"] in figure.LEGEND_KINDS[0]]
    assert [entry["tool"] for entry in named] == ["vep-rs", "fastvep", "ensembl-vep"]
    assert all(tool["name"] in entry["label"] and tool["version"] in entry["label"] for entry, tool in zip(named, document["tools"]))
    # Nothing ties a name to a circle: the panel's drawing has no hairline to draw.
    assert not hasattr(figure, "draw_leader") and "hairline" not in SCRIPT.read_text(encoding="utf-8")


def test_caption_names_every_group_present_and_every_tool_with_its_version(tmp_path):
    caption = values_of(tmp_path)["caption"]
    assert caption.startswith(f"**vep-rs {RELEASE} against Ensembl VEP {REFERENCE_VERSION}, with fastVEP {COMPARATOR_VERSION} as comparator.** (a) ")
    clauses = {
        REAL_GROUP: "plugins verified on real data",
        "Committed table": "TranscriptAnnotator on its committed table",
        "Synthetic fixtures (data that cannot be redistributed)": "plugins whose data cannot be redistributed on synthetic fixtures, marked synthetic",
        "Custom annotation (--custom sources)": "custom-annotation sources on the genome-wide ClinVar datasets",
    }
    # The clauses name the groups; the sentence after them says which datasets the rows scored on real data cover.
    assert f"(c) Raw F1 per plugin: {'; '.join(clauses.values())}. On a row of a plugin verified on real data, " in caption
    assert caption.count(". The bar is each tool's per-feature raw F1") == 1
    assert "LoFTEE's reference ran with stand-ins for the MaxEntScan files its module bundles, on the Ensembl VEP image with samtools added." in caption
    supported = next(block["of"] for block in PLUGINS["support"] if block["engine"] == "vep-rs")
    assert (f"(d) The plugins each tool reads, of the {supported} that vep-rs supports: vep-rs {PLUGINS['support'][0]['read']} and "
            f"fastVEP {PLUGINS['support'][1]['read']}; plugins vep-rs does not port are not counted and are listed by reason in the README.") in caption
    assert caption.count("vep-rs matches Ensembl VEP exactly except for the documented divergences.") == 1
    assert caption.endswith("vep-rs matches Ensembl VEP exactly except for the documented divergences. "
                            "Every value is in the records listed under Provenance records.")
    assert caption.isprintable() and "\n" not in caption

    # A record without synthetic-fixture rows and custom sources draws neither group and its caption names neither.
    plugins = copy.deepcopy(PLUGINS)
    plugins["rows"] = [row for row in plugins["rows"] if row["scored_on"] != "synthetic_fixture"]
    plugins["synthetic_fixtures"], plugins["custom_annotation"] = [], []
    document = values_of(tmp_path / "fewer", record_set(plugins=plugins))
    present = [group["label"] for group in panel(document, "c")["groups"]]
    assert present == [REAL_GROUP, "Committed table"]
    for label, clause in clauses.items():
        assert (clause in document["caption"]) == (label in present), label


def test_the_caption_of_the_plugin_panel_calls_its_bar_raw_f1_in_the_legends_own_words(tmp_path):
    document = values_of(tmp_path)
    caption = document["caption"]
    c = caption[caption.index("(c) "):caption.index("(d) ")]
    # The legend's entry for a bar reads "Bar: raw F1", and the caption of the plugin panel names the bar so: one
    # figure under one name, counted per feature, and no other name for it.
    legend = [entry["label"] for entry in figure.legend_entries(document) if entry["mark"] == "raw"]
    assert legend == [figure.RAW_LEGEND] == ["Bar: raw F1"]
    assert c.startswith("(c) Raw F1 per plugin: ") and "The bar is each tool's per-feature raw F1, the marker its tool's adjusted F1, " in c
    assert "figure" not in c and c.count("raw F1") == 1 and c.count("Raw F1") == 1
    # The bar it names is the row's raw counts at the feature projection.
    row = rows_of(document, "c")[(REAL_GROUP, "LOEUF GRCh37")]
    bar = next(value for value in row["values"] if (value["tool"], value["kind"]) == ("vep-rs", "bar"))
    counts = plugin_sums("real_data", "LOEUF", "GRCh37", "vep-rs", "feature", ("intersection", "perl", "rust"))
    assert bar["counts"] == counts and bar["key"].endswith("projection=feature]: intersection, perl, rust, summed over 2 entries")


def test_the_caption_says_once_that_on_an_f1_axis_a_value_is_a_place_and_not_a_length(tmp_path):
    sentence = ("On an F1 axis a value is a place, not a length: every bar starts at the left edge of its plot, left of where an F1 "
                "of 0 stands, and ends at the place of its F1.")
    # Once, after the panels' own sentences, whatever panels the records give and in a draft too.
    for name, records, flags in (("four", record_set(), ()), ("two", record_set(plugins=None), ()), ("draft", record_set(plugins=ungated()), ("--draft",))):
        caption = values_of(tmp_path / name, records, *flags)["caption"]
        assert caption.count(sentence) == 1 and caption.count("a place, not a length") == 1, name
        last = max(caption.index(f"({panel_id}) ") for panel_id in "abcd" if f"({panel_id}) " in caption)
        assert caption.index(sentence) > last, name
    # What it says is what is drawn: a bar of an F1 of 0 has the room's length, and one of 0.5 ends at its place.
    scale = figure.f1_scale(261.0, 1e-7)
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, f1_row(("vep-rs", "bar", "0.000000000"), ("fastvep", "bar", "0.500000000")))], LANES, FILLS, 0.3, 1e-7, scale)
    (none_at, none), (half_at, half) = axes.drawn("barh")
    assert none["left"] == half["left"] == figure.f1_edge(scale) > 1.0
    assert (none["left"] + none_at[1], half["left"] + half_at[1]) == (pytest.approx(1.0), pytest.approx(0.5))


def test_alternative_text_is_one_line_the_image_syntax_carries_and_says_what_each_panel_holds_as_it_is_drawn(tmp_path):
    lead = f"comparing vep-rs {RELEASE} and fastVEP {COMPARATOR_VERSION} against Ensembl VEP {REFERENCE_VERSION}."
    document = values_of(tmp_path)
    alt = document["alt_text"]
    # The plugin panel's counts are of its own rows: those with the comparator's bar, with its marker, with its
    # cross, and those on which the released tool's adjusted counts are equal.
    rows = list(rows_of(document, "c").values())
    scored, marked, crossed = (sum(("fastvep", kind) in marks(row) for row in rows) for kind in ("bar", "marker", "text"))
    exact = sum(len(set(value["counts"])) == 1 for row in rows for value in row["values"] if (value["tool"], value["kind"]) == ("vep-rs", "marker"))
    assert (len(rows), scored, marked, crossed, exact) == (59, 16, 15, 43, 59)
    assert alt.startswith(
        f"Four panels {lead} (a) Raw F1 per dataset as a bar for each of vep-rs and fastVEP, with each tool's adjusted F1 as a marker "
        "joined to the end of the bar by a line where the two differ. "
        "(b) Wall time per dataset and architecture of vep-rs, fastVEP and Ensembl VEP, a circle each, filled for ARM and open for x86; "
        f"from left to right, the circles of every dataset and architecture are those of vep-rs, fastVEP and Ensembl VEP; {RUNS_CLAUSE}. "
        "(c) Raw F1 per plugin as a bar for each of vep-rs and fastVEP, with each tool's adjusted F1 as a marker joined to the end of the "
        f"bar by a line where the two differ: fastVEP has a figure on {scored} of the {len(rows)} rows, its marker on {marked} of them, and a cross "
        f"(no such source) on {crossed}; vep-rs's adjusted F1 is exactly 1 on {exact} of the {len(rows)} rows. "
        f"(d) The plugins read by vep-rs and fastVEP, a bar each: {READS}.")
    assert alt.isprintable() and not set(alt) & set("[]") and "tick" not in alt
    # The result is computed from the record: one adjusted F1 of the released tool below 1 on a plugin row is one row fewer.
    plugins = copy.deepcopy(PLUGINS)
    next(e for e in plugins["real_data"] if e["engine"] == "vep-rs" and e["projection"] == "feature")["adj_perl"] += 1
    assert f"vep-rs's adjusted F1 is exactly 1 on {exact - 1} of the {len(rows)} rows." in values_of(tmp_path / "below", record_set(plugins=plugins))["alt_text"]
    # A comparator whose every figure carries its marker has no clause on its markers.
    whole = copy.deepcopy(PLUGINS)
    for entry in whole["real_data"]:
        if entry["engine"] == "fastvep" and "adj_perl" not in entry:
            entry.update(adj_intersection=entry["intersection"], adj_perl=entry["perl"], adj_rust=entry["rust"])
    assert (f"fastVEP has a figure on {scored} of the {len(rows)} rows and a cross (no such source) on {crossed}; "
            in values_of(tmp_path / "whole", record_set(plugins=whole))["alt_text"])
    # Each panel's sentence names the tools that panel holds a value of, and no other: a comparator kept out of the
    # concordance panels is named for the wall times and the counts alone, and the markers are the released tool's.
    registry = registry_with(tmp_path, lambda rows: rows["fastvep"].update(concordance_panel="no"))
    alt = values_of(tmp_path / "kept-out", registry=registry)["alt_text"]
    assert alt.startswith(f"Four panels {lead} (a) Raw F1 per dataset as a bar for vep-rs, with its adjusted F1 as a marker joined to the end of "
                   "the bar by a line where the two differ. "
                   "(b) Wall time per dataset and architecture of vep-rs, fastVEP and Ensembl VEP, a circle each, filled for ARM and open for x86; "
                   f"from left to right, the circles of every dataset and architecture are those of vep-rs, fastVEP and Ensembl VEP; {RUNS_CLAUSE}. "
                   "(c) Raw F1 per plugin as a bar for vep-rs, with its adjusted F1 as a marker joined to the end of the bar by a "
                   f"line where the two differ: vep-rs's adjusted F1 is exactly 1 on {exact} of the {len(rows)} rows. "
                   f"(d) The plugins read by vep-rs and fastVEP, a bar each: {READS}.")
    # A comparator without adjusted counts on the datasets: the markers of that panel are named as the released tool's.
    bare, uncounted = dataset_record(), population_record()
    for held_cell in bare["concordance"].values():
        del held_cell["adjusted"]
    for row in uncounted["comparator_concordance"]["per_dataset"]:
        del row["excluded_perl"], row["excluded_rust"]
    alt = values_of(tmp_path / "bare", record_set(fastvep=bare, population=uncounted))["alt_text"]
    assert "(a) Raw F1 per dataset as a bar for each of vep-rs and fastVEP, with vep-rs's adjusted F1 as a marker joined" in alt


def test_each_count_the_alternative_text_gives_of_the_plugin_panel_is_counted_from_its_values_and_kept_among_the_sources(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"

    def panel_counts(document: dict) -> list[tuple[str, str]]:
        """The entries of `sources` that count rows of the plugin panel by what they hold, each as its text and its
        key; each names the record the panel is drawn from."""
        found = [source for source in document["sources"] if source["key"].startswith("panels[id=c]: the rows ")]
        assert all(source["record"] == name for source in found)
        return [(source["text"], source["key"]) for source in found]

    def recount(document: dict) -> tuple[int, int, int, int, int]:
        """Counted again from the panel's values: its rows, those with a figure, a marker and a cross of the
        comparator, and those on which the released tool's adjusted counts are equal."""
        rows = list(rows_of(document, "c").values())
        held = [sum(("fastvep", kind) in marks(row) for row in rows) for kind in ("bar", "marker", "text")]
        exact = sum(len(set(value["counts"])) == 1 for row in rows for value in row["values"] if (value["tool"], value["kind"]) == ("vep-rs", "marker"))
        return len(rows), *held, exact

    def keys(*kinds: str) -> list[str]:
        return [figure.PANEL_EXACT_KEY.format(panel="c", tool="vep-rs") if kind == "exact" else
                figure.PANEL_ROWS_KEY.format(panel="c", tool="fastvep", kind=kind) for kind in kinds]

    assert keys("bar", "exact") == ["panels[id=c]: the rows that hold a value of tool=fastvep, kind=bar: their number",
                                    "panels[id=c]: the rows whose value of tool=vep-rs, kind=marker has equal counts: their number"]
    # On these records the sentence prints four counts of the panel's rows: the comparator's figures, its markers,
    # which are fewer, its crosses, and the rows on which the released tool's adjusted F1 is exactly 1. Each is
    # the count of the panel's own values and stands among the sources in the sentence's order, with what it counts.
    document = values_of(tmp_path)
    rows, scored, marked, crossed, exact = recount(document)
    assert 0 < marked < scored and crossed and exact
    assert (f": fastVEP has a figure on {scored} of the {rows} rows, its marker on {marked} of them, and a cross (no such source) on {crossed}; "
            f"vep-rs's adjusted F1 is exactly 1 on {exact} of the {rows} rows. (d) ") in document["alt_text"]
    assert panel_counts(document) == list(zip(map(str, (scored, marked, crossed, exact)), keys("bar", "marker", "text", "exact")))
    # The counts stand in the alternative text alone: their entries, with the one of the panel's number of rows among
    # them, follow the caption's, whose last is the sentence on an F1 of exactly 1.
    at = [index for index, source in enumerate(document["sources"]) if source["key"].startswith("panels[id=c]: ")]
    after = next(index for index, source in enumerate(document["sources"]) if source["key"] == figure.EXACT_KEY.format(released="vep-rs"))
    assert at == list(range(after + 1, after + 6))

    # A panel that changes changes its counts and their entries with it: in a draft whose record lists two plugins
    # fewer as absent for the comparator, their rows carry no cross, and the count of crosses and its entry say so.
    unlisted = copy.deepcopy(PLUGINS)
    unlisted["absent"] = [item for item in unlisted["absent"] if item["plugin"] not in ("GWAS", "LoFTEE")]
    draft = values_of(tmp_path / "draft", record_set(plugins=unlisted), "--draft")
    fewer = recount(draft)
    assert fewer == (rows, scored, marked, crossed - 2, exact) and f"and a cross (no such source) on {crossed - 2}; " in draft["alt_text"]
    assert panel_counts(draft) == list(zip(map(str, fewer[1:]), keys("bar", "marker", "text", "exact")))

    # A count the sentence does not print has no entry: a comparator kept out of the concordance panels has no
    # clause, and the one count left is the released tool's.
    registry = registry_with(tmp_path, lambda rows: rows["fastvep"].update(concordance_panel="no"))
    alone = values_of(tmp_path / "kept-out", registry=registry)
    assert "fastVEP has a figure" not in alone["alt_text"] and panel_counts(alone) == [(str(exact), keys("exact")[0])]
    # Without the plugin record there is no plugin panel, and no count of one.
    assert panel_counts(values_of(tmp_path / "none", record_set(plugins=None))) == []


def test_each_number_of_rows_the_alternative_text_prints_is_counted_from_its_panels_rows_and_kept_among_the_sources(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"

    def kept(document: dict, panel_id: str) -> list[tuple[str, str, str]]:
        """The entries of `sources` that count rows of a panel, each as its text, its record and its key."""
        return [(source["text"], source["record"], source["key"]) for source in document["sources"]
                if source["key"].startswith(f"panels[id={panel_id}]: ")]

    # The plugin panel. The sentence closes two of its counts with the number of the panel's rows: that number is
    # the count of the panel's rows, and it is kept once, where the sentence first prints it, after the first count.
    document = values_of(tmp_path)
    rows = len(rows_of(document, "c"))
    total = "panels[id=c]: its rows: their number"
    assert total == figure.PANEL_TOTAL_KEY.format(panel="c") and rows > len(PLUGINS["rows"])
    assert document["alt_text"].count(f" of the {rows} rows") == 2
    assert [entry for entry in kept(document, "c") if entry[2] == total] == [kept(document, "c")[1]] == [(str(rows), name, total)]
    # A panel of other rows prints and keeps their number: the gate rows alone, and one row, which reads "row".
    plain = copy.deepcopy(PLUGINS)
    plain["custom_annotation"] = []
    gate = values_of(tmp_path / "gate", record_set(plugins=plain))
    assert len(rows_of(gate, "c")) == len(PLUGINS["rows"]) and kept(gate, "c")[1] == (str(len(PLUGINS["rows"])), name, total)
    assert gate["alt_text"].count(f" of the {len(PLUGINS['rows'])} rows") == 2
    plain["rows"] = [row for row in plain["rows"] if (row["plugin"], row["assembly"]) == ("LOEUF", "GRCh37")]
    one = values_of(tmp_path / "one", record_set(plugins=plain))
    assert one["alt_text"].count(" of the 1 row") == 2 and kept(one, "c")[1] == ("1", name, total) and len(kept(one, "c")) == 3
    # Where the released tool's count is the first to print it, the number follows that count; no plugin panel, none.
    registry = registry_with(tmp_path, lambda rows: rows["fastvep"].update(concordance_panel="no"))
    alone = kept(values_of(tmp_path / "kept-out", registry=registry), "c")
    assert [key for _, _, key in alone] == [figure.PANEL_EXACT_KEY.format(panel="c", tool="vep-rs"), total] and alone[1] == (str(rows), name, total)
    assert kept(values_of(tmp_path / "none", record_set(plugins=None)), "c") == []

    # The wall-time panel. One order of the tools' circles on every row prints no number of rows and keeps none.
    assert "the circles of every dataset and architecture are those of" in document["alt_text"] and kept(document, "b") == []

    def orders(document: dict) -> dict[tuple[str, ...], int]:
        """Counted again from the panel's values: the rows by the order of their tools from the shortest time."""
        found: dict[tuple[str, ...], int] = {}
        for row in rows_of(document, "b").values():
            order = tuple(value["tool"] for value in sorted(row["values"], key=lambda value: Decimal(value["value"])))
            found[order] = found.get(order, 0) + 1
        return found

    def entries(document: dict) -> list[tuple[str, str, str]]:
        """What the clause keeps: each order's number of rows as the clause prints the orders, the number of rows
        that hold a circle after the first, each with the records of the panel's values."""
        record = " and ".join(dict.fromkeys(value["record"] for row in rows_of(document, "b").values() for value in row["values"]))
        found = orders(document)
        each = [(str(number), record, figure.PANEL_ORDER_KEY.format(panel="b", tools=", ".join(f"tool={tool}" for tool in order)))
                for order, number in found.items()]
        return [each[0], (str(sum(found.values())), record, figure.PANEL_TIMED_KEY.format(panel="b")), *each[1:]]

    # A comparator slower than the reference on one dataset stands right of it on that dataset's two rows: the
    # clause names both orders, each with its number of rows, and each number is kept with the order it counts.
    population = population_record()
    for block in population["wall_time"]["fastvep"]["per_dataset"]:
        if block["dataset"] == "pop_b":
            block["sum_of_cell_medians_sec"] += 5000
    two = values_of(tmp_path / "two-orders", record_set(population=population))
    assert orders(two) == {("vep-rs", "fastvep", "ensembl-vep"): 6, ("vep-rs", "ensembl-vep", "fastvep"): 2}
    assert "are those of vep-rs, fastVEP and Ensembl VEP on 6 of the 8 and vep-rs, Ensembl VEP and fastVEP on 2 of the 8; " in two["alt_text"]
    assert kept(two, "b") == entries(two) and [text for text, _, _ in kept(two, "b")] == ["6", "8", "2"]
    assert [key for _, _, key in kept(two, "b")] == [
        "panels[id=b]: the rows whose values of kind=marker, from the smallest, are of tool=vep-rs, tool=fastvep, tool=ensembl-vep: their number",
        "panels[id=b]: the rows that hold a value of kind=marker: their number",
        "panels[id=b]: the rows whose values of kind=marker, from the smallest, are of tool=vep-rs, tool=ensembl-vep, tool=fastvep: their number"]
    assert kept(two, "b")[0][1] == " and ".join(f"2001-01-01-{kind}-{RELEASE}.json" for kind in ("release", "comparator-walltime", "population"))
    # The entries stand where the alternative text prints them: after the caption's and before the plugin panel's.
    at = [index for index, source in enumerate(two["sources"]) if source["key"].startswith("panels[id=")]
    assert at == list(range(at[0], at[0] + 8)) and [two["sources"][index]["key"][:13] for index in at] == ["panels[id=b]:"] * 3 + ["panels[id=c]:"] * 5
    # Rows that hold fewer tools are an order of their own, with its number of rows.
    fewer = values_of(tmp_path / "fewer", record_set(walltime=None))
    assert orders(fewer) == {("vep-rs",): 4, ("vep-rs", "fastvep", "ensembl-vep"): 4}
    assert kept(fewer, "b") == entries(fewer) and [text for text, _, _ in kept(fewer, "b")] == ["4", "8", "4"]
    # One tool on every row has no order to name, and no number of rows.
    registry = registry_with(tmp_path, lambda rows: [rows[tool].update(speed_panel="no") for tool in ("fastvep", "ensembl-vep")])
    assert kept(values_of(tmp_path / "single", registry=registry), "b") == []


def test_the_alternative_text_names_the_order_of_the_tools_circles_on_a_row_of_the_wall_time_panel(tmp_path):
    lead = "(b) Wall time per dataset and architecture of vep-rs, fastVEP and Ensembl VEP, a circle each, filled for ARM and open for x86"
    names = {"vep-rs": "vep-rs", "fastvep": "fastVEP", "ensembl-vep": "Ensembl VEP"}

    def orders(document: dict) -> list[list[str]]:
        """Per row of the wall-time panel, its tools from the shortest time; equal times keep the tools' own order."""
        return [[names[value["tool"]] for value in sorted(row["values"], key=lambda value: Decimal(value["value"]))]
                for row in rows_of(document, "b").values()]

    # The panel tells its tools apart by colour alone. One order holds on every row of these records: the
    # alternative text names it once, after the panel's marks, and the caption says nothing of it.
    document = values_of(tmp_path)
    assert orders(document) == [["vep-rs", "fastVEP", "Ensembl VEP"]] * 8
    assert f"{lead}; from left to right, the circles of every dataset and architecture are those of vep-rs, fastVEP and Ensembl VEP; {RUNS_CLAUSE}. (c) " in document["alt_text"]
    assert document["alt_text"].count("from left to right") == 1 and "from left to right" not in document["caption"]
    # The order is each row's own, from its values: a comparator slower than the reference on one dataset stands
    # right of it on that dataset's two rows, and the text names both orders, each with its number of rows.
    population = population_record()
    for block in population["wall_time"]["fastvep"]["per_dataset"]:
        if block["dataset"] == "pop_b":
            block["sum_of_cell_medians_sec"] += 5000
    document = values_of(tmp_path / "two-orders", record_set(population=population))
    assert orders(document).count(["vep-rs", "Ensembl VEP", "fastVEP"]) == 2 and orders(document).count(["vep-rs", "fastVEP", "Ensembl VEP"]) == 6
    assert (f"{lead}; from left to right, the circles of a dataset and architecture are those of vep-rs, fastVEP and Ensembl VEP on 6 of the 8 "
            f"and vep-rs, Ensembl VEP and fastVEP on 2 of the 8; {RUNS_CLAUSE}. (c) ") in document["alt_text"]
    # Two tools with one time on a row keep the tools' own order there, and the one order is still named for every row.
    population = population_record()
    slowest = {(block["dataset"], block["arch"]): block["sum_of_cell_medians_sec"] for block in population["wall_time"]["ensembl_vep"]["per_dataset"]}
    for block in population["wall_time"]["fastvep"]["per_dataset"]:
        if block["dataset"] == "pop_b":
            block["sum_of_cell_medians_sec"] = slowest[(block["dataset"], block["arch"])]
    document = values_of(tmp_path / "equal", record_set(population=population))
    level = [row for row in rows_of(document, "b").values() if marks(row)[("fastvep", "marker")] == marks(row)[("ensembl-vep", "marker")]]
    assert len(level) == 2 and orders(document) == [["vep-rs", "fastVEP", "Ensembl VEP"]] * 8
    assert f"{lead}; from left to right, the circles of every dataset and architecture are those of vep-rs, fastVEP and Ensembl VEP; {RUNS_CLAUSE}. (c) " in document["alt_text"]
    # Rows that hold fewer tools are named with the tools they hold: without the comparators' wall-time record the
    # ClinVar rows carry the released tool's circle alone.
    document = values_of(tmp_path / "fewer", record_set(walltime=None))
    assert orders(document).count(["vep-rs"]) == 4 and orders(document).count(["vep-rs", "fastVEP", "Ensembl VEP"]) == 4
    assert ("; from left to right, the circles of a dataset and architecture are those of vep-rs on 4 of the 8 and vep-rs, fastVEP and Ensembl VEP "
            f"on 4 of the 8; {RUNS_CLAUSE}. (c) ") in document["alt_text"]
    # One tool on every row has no order to name.
    registry = registry_with(tmp_path, lambda rows: [rows[name].update(speed_panel="no") for name in ("fastvep", "ensembl-vep")])
    alone = values_of(tmp_path / "alone", registry=registry)["alt_text"]
    assert (f"(b) Wall time per dataset and architecture of vep-rs, a circle each, filled for ARM and open for x86; each wall-time cell holds {TOOL_RUNS} "
            "runs. (c) ") in alone and "from left to right" not in alone


def test_the_alternative_text_reads_of_the_one_row_where_the_plugin_panel_holds_one_row(tmp_path):
    lead = "with each tool's adjusted F1 as a marker joined to the end of the bar by a line where the two differ: "
    # The plugin panel's counts are of its rows: with two rows and more the text reads "rows".
    plugins = copy.deepcopy(PLUGINS)
    plugins["custom_annotation"] = []
    scored = [row for row in plugins["rows"] if (row["plugin"], row["assembly"]) == ("LOEUF", "GRCh37")]
    crossed = [row for row in plugins["rows"] if any((item["plugin"], item["engine"]) == (row["plugin"], "fastvep") for item in plugins["absent"])][:1]
    assert len(scored) == 1 and len(crossed) == 1
    plugins["rows"] = scored + crossed
    two = values_of(tmp_path / "two", record_set(plugins=plugins))["alt_text"]
    assert (f"{lead}fastVEP has a figure on 1 of the 2 rows and a cross (no such source) on 1; vep-rs's adjusted F1 is exactly 1 on 2 of the 2 rows. (d) ") in two
    # One row, scored for both tools: "of the 1 row", for the comparator's clause and for the released tool's.
    plugins["rows"] = scored
    one = values_of(tmp_path / "one", record_set(plugins=plugins))
    assert sum(len(group["rows"]) for group in panel(one, "c")["groups"]) == 1
    assert f"{lead}fastVEP has a figure on 1 of the 1 row; vep-rs's adjusted F1 is exactly 1 on 1 of the 1 row. (d) " in one["alt_text"]
    assert "1 rows" not in one["alt_text"]
    # One row with the comparator's cross alone.
    plugins["rows"] = crossed
    alone = values_of(tmp_path / "crossed", record_set(plugins=plugins))["alt_text"]
    assert "fastVEP has a figure on 0 of the 1 row and a cross (no such source) on 1; vep-rs's adjusted F1 is exactly 1 on 1 of the 1 row. (d) " in alone
    assert "1 rows" not in alone


def test_the_wall_time_panel_says_which_lines_are_medians_and_which_are_sums(tmp_path):
    document = values_of(tmp_path)
    wall = panel(document, "b")
    # The ClinVar lines are medians and the population datasets' lines sums of medians, by the keys of their values.
    for group in wall["groups"]:
        ending = ".median_sec" if group["label"].startswith("ClinVar ") else ".sum_of_cell_medians_sec"
        assert all(value["key"].endswith(ending) for row in group["rows"] for value in row["values"]), group["label"]
    assert figure.median_lines(wall) == (["ClinVar"], True)
    # The axis says so under its first line, among the drawn strings, and the caption says it in full.
    lines = figure.axis_lines(document, figure.read_registry(REGISTRY), wall)
    assert lines[1] == "ClinVar: medians; other datasets: sums of medians" and len(lines) == 2
    assert lines[1] in figure.drawn_texts(document, figure.read_registry(REGISTRY))
    assert ("(b) Wall time per dataset and architecture: a ClinVar line is the median of its runs, any other the sum of its dataset's "
            f"cell medians; {RUNS_CLAUSE}.") in document["caption"]


def test_the_caption_says_a_sum_of_cell_medians_and_what_a_cell_is_in_the_records_own_words(tmp_path):
    name = f"2001-01-01-population-{RELEASE}.json"
    # The values' key is the record's: a sum of cell medians. The caption uses that word, and no other of its own.
    document = values_of(tmp_path)
    sums = [value["key"] for group in panel(document, "b")["groups"] for row in group["rows"] for value in row["values"] if "sum_of" in value["key"]]
    assert sums and all(key.endswith(".sum_of_cell_medians_sec") for key in sums)
    assert f"any other the sum of its dataset's cell medians; {RUNS_CLAUSE}." in document["caption"] and "per-chromosome" not in document["caption"]
    # A record without a methodology note gives no gloss, and the script says so on its output.
    done, _ = run(tmp_path / "noted", record_set())
    assert f"note: {name} carries no methodology note that says what a dataset's sum is of; the caption omits it\n" in done.stdout
    assert not any(source["key"] == "methodology_note" for source in document["sources"])
    # What a cell is stands in the record's methodology note: the caption quotes the clause that says what a
    # dataset's sum is of, and the values file names the key it is drawn from.
    population = population_record()
    population["methodology_note"] = ("Five machines per chromosome per architecture, each running one timed run; cell quantiles by one "
                                      "estimator; a dataset's figure is the sum of its per-chromosome medians. The inputs are the ones scored.")
    document = values_of(tmp_path / "glossed", record_set(population=population))
    clause = "a dataset's figure is the sum of its per-chromosome medians"
    assert ("(b) Wall time per dataset and architecture: a ClinVar line is the median of its runs, any other the sum of its dataset's "
            f'cell medians (the record\'s methodology note: "{clause}"); {RUNS_CLAUSE}. (c) ') in document["caption"]
    assert document["sources"][0] == {"text": clause, "record": name, "key": "methodology_note"}
    assert clause in population["methodology_note"] and document["caption"].count(clause) == 1
    # A note that speaks of no sum of medians, or one that is no text, is quoted nowhere.
    for index, note in enumerate(("Five machines per chromosome; quantiles by one estimator.", ["a dataset's figure is the sum of its medians"])):
        population["methodology_note"] = note
        document = values_of(tmp_path / f"other{index}", record_set(population=population))
        assert "methodology note" not in document["caption"] and not any(source["key"] == "methodology_note" for source in document["sources"])
    # A panel of medians alone, or of sums alone, says neither word.
    assert figure.median_lines({"groups": []}) == ([], False)


def test_the_texts_say_how_many_runs_a_wall_time_cell_holds_from_the_numbers_the_records_state(tmp_path):
    release, population, walltime = (f"2001-01-01-{kind}-{RELEASE}.json" for kind in ("release", "population", "comparator-walltime"))

    def wall_texts(document: dict) -> list[str]:
        """The wall-time panel's sentence in the caption and in the alternative text."""
        return [text[text.index("(b) "):text.index(" (c) ")] for text in (document["caption"], document["alt_text"])]

    def run_sources(document: dict) -> list[tuple[str, str, list[str]]]:
        """The entries of `sources` that hold a number of runs: its text, its record and its keys."""
        return [(source["text"], source["record"], source["key"].split("; ")) for source in document["sources"]
                if source["key"].endswith((".n_clones", ".n_clones_per_cell"))]

    # A cell of the reference holds one run and a cell of either other tool as many as the other's: one clause, the
    # same in both texts, closes the panel's sentence, and the script prints no note on it.
    done, values = run(tmp_path, record_set())
    document = json.loads(values.read_text(encoding="utf-8"))
    caption, alt = wall_texts(document)
    assert caption.endswith(f"cell medians; {RUNS_CLAUSE}.") and alt.endswith(f"those of vep-rs, fastVEP and Ensembl VEP; {RUNS_CLAUSE}.")
    assert document["caption"].count(RUNS_CLAUSE) == document["alt_text"].count(RUNS_CLAUSE) == 1 and "number of runs" not in done.stdout
    # Each number is read beside a value: `n_clones` of a ClinVar entry, `n_clones_per_cell` of a sum's entry and,
    # where the entry states none, as the reference's do, of its engine's block. The values file keeps each number
    # with its keys, in the clause's order, tool by tool and record by record.
    cells = [f"suite_id={cell_id},arch={arch}" for cell_id in ("s01", "s02") for arch in ARCHES]
    sums = [f"per_dataset[dataset={dataset},arch={arch}]" for dataset in ("pop_a", "pop_b") for arch in ARCHES]
    assert run_sources(document) == [
        (str(REFERENCE_RUNS), walltime, [f"per_cell_aggregates[engine=perl,{cell}].n_clones" for cell in cells]),
        (str(REFERENCE_RUNS), population, ["wall_time.ensembl_vep.n_clones_per_cell"]),
        (str(TOOL_RUNS), release, [f"per_cell_aggregates[{cell},no target_cpu].n_clones" for cell in cells]),
        (str(TOOL_RUNS), population, [f"wall_time.vep_rs.{entry}.n_clones_per_cell" for entry in sums]),
        (str(TOOL_RUNS), walltime, [f"per_cell_aggregates[engine=fastvep,{cell}].n_clones" for cell in cells]),
        (str(TOOL_RUNS), population, [f"wall_time.fastvep.{entry}.n_clones_per_cell" for entry in sums]),
    ]

    # Where every tool's cells hold as many runs, the clause names no tool.
    level_wall, level_population = walltime_record(), population_record()
    for entry in level_wall["per_cell_aggregates"]:
        entry["n_clones"] = TOOL_RUNS
    level_population["wall_time"]["ensembl_vep"]["n_clones_per_cell"] = TOOL_RUNS
    document = values_of(tmp_path / "level", record_set(walltime=level_wall, population=level_population))
    assert all(text.endswith(f"; each wall-time cell holds {TOOL_RUNS} runs.") for text in wall_texts(document))
    assert [text for text, _, _ in run_sources(document)] == [str(TOOL_RUNS)] * 6

    # Another number in a record is the number the clause gives: each number with its tools, the fewest first, and
    # a span for a tool whose cells differ. A sum's own entry is read before its engine's block, which here still
    # states the first number.
    other_wall, other_population = walltime_record(), population_record()
    for entry in other_wall["per_cell_aggregates"]:
        if entry["engine"] == "fastvep":
            entry["n_clones"] = TOOL_RUNS - 1
    for entry in other_population["wall_time"]["fastvep"]["per_dataset"]:
        entry["n_clones_per_cell"] = TOOL_RUNS - 1
    for entry in other_population["wall_time"]["vep_rs"]["per_dataset"]:
        entry["n_clones_per_cell"] = TOOL_RUNS + 4
    document = values_of(tmp_path / "other", record_set(walltime=other_wall, population=other_population))
    clause = (f"each wall-time cell of Ensembl VEP holds {REFERENCE_RUNS} run, against {TOOL_RUNS - 1} for fastVEP and "
              f"{TOOL_RUNS} to {TOOL_RUNS + 4} for vep-rs")
    assert all(text.endswith(f"; {clause}.") for text in wall_texts(document))
    assert [(text, record) for text, record, _ in run_sources(document)] == [
        (str(REFERENCE_RUNS), walltime), (str(REFERENCE_RUNS), population), (str(TOOL_RUNS - 1), walltime), (str(TOOL_RUNS - 1), population),
        (f"{TOOL_RUNS} to {TOOL_RUNS + 4}", release), (f"{TOOL_RUNS} to {TOOL_RUNS + 4}", population)]

    # A record that states no number for a cell, or none that is a whole number of one run or more, gives no
    # clause: both texts leave the runs of a cell out, the values file keeps no number of runs, and a note is printed.
    for index, runs in enumerate((None, 0, True, "5", 2.5)):
        bare = walltime_record()
        if runs is None:
            del bare["per_cell_aggregates"][-1]["n_clones"]
        else:
            bare["per_cell_aggregates"][-1]["n_clones"] = runs
        done, values = run(tmp_path / f"bare{index}", record_set(walltime=bare))
        assert done.returncode == 0, done.stderr
        assert f"note: {walltime} states no number of runs for a wall-time cell of Ensembl VEP; the texts leave the runs of a cell out\n" in done.stdout
        document = json.loads(values.read_text(encoding="utf-8"))
        caption, alt = wall_texts(document)
        assert caption.endswith("cell medians.") and alt.endswith("those of vep-rs, fastVEP and Ensembl VEP.") and run_sources(document) == []
        assert "wall-time cell" not in document["caption"] + document["alt_text"]
    blockless = population_record()
    del blockless["wall_time"]["ensembl_vep"]["n_clones_per_cell"]
    done, values = run(tmp_path / "blockless", record_set(population=blockless))
    assert f"note: {population} states no number of runs for a wall-time cell of Ensembl VEP; the texts leave the runs of a cell out\n" in done.stdout
    assert run_sources(json.loads(values.read_text(encoding="utf-8"))) == []


def test_the_wall_time_axis_says_which_way_is_better(tmp_path):
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    # The F1 axes say that further right is better; the wall-time axis, where it is the other way, says so too.
    assert figure.axis_lines(document, registry, panel(document, "a"))[1].endswith("further right is better")
    assert figure.axis_lines(document, registry, panel(document, "b"))[0] == "Wall time (s), log scale; further left is better"


def test_the_f1_axis_is_labelled_f1_against_the_reference_and_the_legend_says_that_a_bar_is_the_raw_f1(tmp_path):
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    # The axis carries raw bars and adjusted markers, so its label names neither kind, and no version.
    for panel_id in ("a", "c"):
        lines = figure.axis_lines(document, registry, panel(document, panel_id))
        assert lines == ["F1 against Ensembl VEP", "log(1 − F1) scale; further right is better"]
        assert not any(word in line for line in lines for word in ("Raw", "raw", "djusted", REFERENCE_VERSION))
    drawn = figure.drawn_texts(document, registry)
    assert drawn.count("F1 against Ensembl VEP") == 2 and not any(text.startswith("Raw F1") for text in drawn)
    # The legend says what a bar is, once, for every tool, first among the marks of an F1 and in no tool's colour.
    entries = figure.legend_entries(document)
    assert [(entry["tool"], entry["label"]) for entry in entries if entry["mark"] == "raw"] == [(None, "Bar: raw F1")]
    assert figure.legend_rows(entries)[1][0]["mark"] == "raw" and drawn.count("Bar: raw F1") == 1
    # A panel of counts and wall times alone has no bar of an F1 to explain.
    bare = copy.deepcopy(document)
    bare["panels"] = [p for p in bare["panels"] if p["id"] in ("b", "d")]
    assert "raw" not in {entry["mark"] for entry in figure.legend_entries(bare)}


def test_no_dagger_stands_in_the_image_and_the_caption_says_nothing_of_a_tools_gene_models(tmp_path):
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    # The registry has no column on a tool's gene models, and the values file states nothing of them.
    assert "gene_models" not in figure.REGISTRY_COLUMNS and "gene_models" not in REGISTRY.read_text(encoding="utf-8")
    assert not any("own_gene_models" in tool for tool in document["tools"]) and not hasattr(figure, "OWN_GENE_MODELS")
    assert all(set(tool) == {"id", "name", "version", "colour"} for tool in document["tools"])
    # No mark and no word of the image says so: the legend's entry of the tool is its name and version alone, the
    # legend has no row of words under its tools, and no string drawn carries a dagger.
    rows = figure.legend_rows(figure.legend_entries(document))
    assert [entry["label"] for entry in rows[0]] == [f"vep-rs {RELEASE}", f"fastVEP {COMPARATOR_VERSION}", f"Ensembl VEP {REFERENCE_VERSION}"]
    assert "gloss" not in {entry["mark"] for row in rows for entry in row} and not any("gloss" in kinds for kinds in figure.LEGEND_KINDS)
    texts = [*figure.drawn_texts(document, registry), document["caption"], document["alt_text"], *document["notes"]]
    assert not any(word in text for text in texts for word in ("†", "dagger", "gene model"))
    assert not hasattr(figure, "DAGGER_GLOSS") and not hasattr(figure, "DAGGER_LEGEND")
    # A registry that carries such a column is refused by name, as any other header is.
    old = tmp_path / "with-gene-models.tsv"
    old.write_text(REGISTRY.read_text(encoding="utf-8").replace("speed_panel\theadline_extra_args", "speed_panel\tgene_models\theadline_extra_args"), encoding="utf-8")
    with pytest.raises(SystemExit) as refused:
        figure.read_registry(old)
    assert "with-gene-models.tsv opens with the columns" in str(refused.value) and "speed_panel, gene_models, headline_extra_args" in str(refused.value)


def test_the_concordance_sentence_stands_only_over_records_whose_every_adjusted_f1_of_the_released_tool_is_one(tmp_path):
    sentence = "vep-rs matches Ensembl VEP exactly except for the documented divergences."
    unshipped = "This is a draft of a state that is not shipped: not every adjusted F1 of vep-rs is exactly 1."
    following = " Every value is in the records"
    # Every adjusted F1 of the released tool exactly 1: the sentence stands, once, whatever a comparator's are.
    document = values_of(tmp_path)
    ours = [value for panel_id in ("a", "c") for group in panel(document, panel_id)["groups"] for row in group["rows"]
            for value in row["values"] if value["kind"] == "marker"]
    assert all(len(set(value["counts"])) == 1 for value in ours if value["tool"] == "vep-rs")
    assert any(len(set(value["counts"])) > 1 for value in ours if value["tool"] == "fastvep")
    assert document["caption"].count(sentence) == 1 and unshipped not in document["caption"]
    assert f"{sentence}{following}" in document["caption"]
    # One adjusted F1 of the released tool below 1, on a dataset or on a plugin: the sentence is not printed, and
    # in its place the caption says that the render is a draft of a state that is not shipped.
    release = release_record()
    release["concordance"]["s08"]["adj_perl"] += 1
    plugins = copy.deepcopy(PLUGINS)
    entry = next(e for e in plugins["real_data"] if e["engine"] == "vep-rs" and e["projection"] == "feature")
    entry["adj_perl"] += 1
    for name, records in (("dataset", record_set(release=release)), ("plugin", record_set(plugins=plugins))):
        caption = values_of(tmp_path / name, records)["caption"]
        assert sentence not in caption and caption.count(unshipped) == 1, name
        assert f"{unshipped}{following}" in caption, name


def test_a_structural_variant_row_says_that_its_figure_pools_its_files(tmp_path):
    # The row's F1 is the cell's counts, which are its files' counts pooled: the label says so with their number.
    rows = rows_of(values_of(tmp_path), "a")
    assert [label for group, label in rows if group == "Structural variants"] == [
        "SV per-VCF, 2 files pooled GRCh37", "SV per-VCF, 2 files pooled GRCh38"]
    bar = rows[("Structural variants", "SV per-VCF, 2 files pooled GRCh37")]["values"][0]
    assert (bar["value"], bar["counts"], bar["key"]) == ("0.500000000", [1, 2, 2], "concordance.s07: intersection, perl, rust")
    # It is no mean of the files' figures: files that read 1 and 0.25 pool to 0.5, and their mean is 0.625.
    release = release_record()
    release["concordance"]["s07"] = {**cell(3, 6, 6), "per_file": {"one": file_counts(2, 2, 2), "two": file_counts(1, 4, 4)}}
    assert (nine_places(2, 2, 2), nine_places(1, 4, 4)) == ("1.000000000", "0.250000000")
    pooled = rows_of(values_of(tmp_path / "pooled", record_set(release=release)), "a")[("Structural variants", "SV per-VCF, 2 files pooled GRCh37")]
    assert marks(pooled)[("vep-rs", "bar")] == nine_places(3, 6, 6) == "0.500000000"
    # A cell whose files do not sum to its counts is refused by name: the label would say what the figure is not.
    release["concordance"]["s07"]["per_file"]["two"] = file_counts(1, 4, 5)
    assert refusal(tmp_path / "unsummed", record_set(release=release)) == (
        f"2001-01-01-release-{RELEASE}.json concordance.s07 counts 3 matched of 6 and 6, and its 2 files sum to 3 of 6 and 7; the row's "
        "figure is the files' counts pooled")
    # A cell that lists no files is named without their number, and a file without its counts is refused.
    del release["concordance"]["s07"]["per_file"]
    assert ("Structural variants", "SV per-VCF GRCh37") in rows_of(values_of(tmp_path / "unlisted", record_set(release=release)), "a")
    # The name is short: with the sixteen files of a release it takes the room of a dataset's name with its assembly.
    assert len("SV per-VCF, 16 files pooled GRCh38") <= len("1000 Genomes high-coverage GRCh38") + 1
    release["concordance"]["s07"]["per_file"] = {"one": {}}
    assert refusal(tmp_path / "uncounted", record_set(release=release)) == (
        f"2001-01-01-release-{RELEASE}.json concordance.s07.per_file carries no count under intersection: None")


def test_a_pooled_figure_that_draws_on_two_records_names_both_as_its_source(tmp_path):
    release, population = f"2001-01-01-release-{RELEASE}.json", f"2001-01-01-population-{RELEASE}.json"
    scored = f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json"
    document = values_of(tmp_path)
    pooled = rows_of(document, "a")[(POOLED, "All four")]["values"]
    # The released tool's pool is summed inside the release record. The comparator's sums the population record's
    # datasets and the ClinVar cells of its dataset record: its bar and its marker name both files, as their keys do.
    assert [(value["tool"], value["kind"], value["record"]) for value in pooled] == [
        ("vep-rs", "bar", release), ("vep-rs", "marker", release),
        ("fastvep", "bar", f"{population} and {scored}"), ("fastvep", "marker", f"{population} and {scored}")]
    assert all(scored in value["key"] and "comparator_concordance.per_dataset" in value["key"] for value in pooled[2:])
    # Every other value of the figure names the one record it is read from.
    listed = {record["file"] for record in document["records"]}
    others = [value["record"] for p in document["panels"] for group in p["groups"] for row in group["rows"] for value in row["values"]
              if value not in pooled[2:]]
    assert set(others) <= listed and {population, scored} <= listed


def test_the_pooled_row_says_once_in_full_what_it_pools(tmp_path):
    document = values_of(tmp_path)
    # The group over the pooled row says in full that it pools the four whole-genome datasets; the row keeps the
    # README's name for it.
    pooled = panel(document, "a")["groups"][2]
    assert (pooled["label"], [row["label"] for row in pooled["rows"]]) == ("The four whole-genome datasets, pooled", ["All four"])
    drawn = figure.drawn_texts(document, figure.read_registry(REGISTRY))
    assert drawn.count("The four whole-genome datasets, pooled") == 1 and drawn.count("All four") == 1
    assert document["caption"].count("All four is the pool of the four whole-genome datasets, the ClinVar and the population datasets.") == 1


# Refusals


def test_registry_of_more_than_eight_tools_is_refused_naming_small_multiples(tmp_path):
    def nine(rows):
        for index in range(6):
            rows[f"tool{index}"] = {**rows["fastvep"], "id": f"tool{index}", "record_keys": f"tool{index}"}

    message = refusal(tmp_path, record_set(), registry=registry_with(tmp_path, nine))
    assert message == ("registry.tsv lists 9 tools; past 8 a comparison is drawn as small multiples, one panel per tool, "
                       "which this figure is not")


@pytest.mark.parametrize("count", [4, 8])
def test_registry_of_more_tools_than_the_figure_has_fills_is_refused_by_name(tmp_path, count):
    def more(rows):
        for index in range(count - 3):
            rows[f"tool{index}"] = {**rows["fastvep"], "id": f"tool{index}", "record_keys": f"tool{index}"}

    message = refusal(tmp_path, record_set(), registry=registry_with(tmp_path, more))
    assert message == f"registry.tsv lists {count} tools and the figure holds 3 fills, one per tool; a tool past them has no colour of its own"


@pytest.mark.parametrize(("edit", "message"), [
    (lambda rows: rows["fastvep"].update(role="rival"), "registry.tsv line 3 carries role 'rival'; it is one of released, comparator, reference"),
    (lambda rows: rows["fastvep"].update(colour_slot="1"), "registry.tsv repeats a colour slot"),
    (lambda rows: rows["fastvep"].update(colour_slot="9"), "registry.tsv line 3 carries colour_slot '9'; it is one of 1, 2, 3"),
    (lambda rows: rows["fastvep"].update(role="released"), "registry.tsv lists 2 tools of the role released; the figure reads one"),
    (lambda rows: rows.pop("ensembl-vep"), "registry.tsv lists 0 tools of the role reference; the figure reads one"),
])
def test_malformed_registry_is_refused_by_name(tmp_path, edit, message):
    assert refusal(tmp_path, record_set(), registry=registry_with(tmp_path, edit)) == message


def test_missing_record_is_refused_by_name(tmp_path):
    assert refusal(tmp_path, record_set(release=None)) == "no release record is among the records"
    assert refusal(tmp_path / "b", record_set(population=None)) == "no population record is among the records"


def test_ambiguous_record_is_refused_by_name(tmp_path):
    records = record_set()
    records[f"2001-01-02-population-{RELEASE}.json"] = population_record()
    assert refusal(tmp_path, records) == (
        f"2 population records are given (2001-01-01-population-{RELEASE}.json, 2001-01-02-population-{RELEASE}.json); the figure reads one")
    records = record_set()
    records[f"2001-01-02-comparator-fastvep-{COMPARATOR_VERSION}.json"] = dataset_record()
    assert refusal(tmp_path / "b", records) == (
        f"2 comparator dataset records hold a run of fastVEP (2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json, "
        f"2001-01-02-comparator-fastvep-{COMPARATOR_VERSION}.json) and 2 of them ran it with '{HEADLINE}'; a panel reads one run")


def test_a_named_file_that_is_no_record_is_refused_by_name(tmp_path):
    stray = tmp_path / "stray.json"
    stray.write_text('{"record": "something-else"}\n', encoding="utf-8")
    done = subprocess.run([sys.executable, str(SCRIPT), "--repo", str(tmp_path), "--records", str(stray), "--values", str(tmp_path / "v.json"),
                           "--values-only"], capture_output=True, text=True, check=False)
    assert done.returncode == 1
    assert done.stderr == "ERROR: [render_release_figure] stray.json is no record the figure reads: a release, population, comparator or plugin record\n"


def test_a_file_of_a_records_directory_that_is_no_object_is_passed_over(tmp_path):
    plain, plain_values = run(tmp_path / "plain", record_set())
    strays = tmp_path / "strays" / "records"
    strays.mkdir(parents=True)
    (strays / "a-list.json").write_text("[1, 2]\n", encoding="utf-8")
    (strays / "a-number.json").write_text("7\n", encoding="utf-8")
    # So is an object whose `record` is no string.
    (strays / "a-list-for-a-record.json").write_text('{"record": [1, 2]}\n', encoding="utf-8")
    done, values = run(tmp_path / "strays", record_set())
    assert plain.returncode == done.returncode == 0, plain.stderr + done.stderr
    assert values.read_bytes() == plain_values.read_bytes()
    # Named, the same file is refused as any file that is no record.
    named = subprocess.run([sys.executable, str(SCRIPT), "--repo", str(tmp_path / "repo"), "--records", str(strays / "a-list.json"),
                            "--values", str(tmp_path / "v.json"), "--values-only"], capture_output=True, text=True, check=False)
    assert named.returncode == 1
    assert named.stderr == "ERROR: [render_release_figure] a-list.json is no record the figure reads: a release, population, comparator or plugin record\n"


def test_a_record_lacking_a_field_the_released_tool_reads_is_refused_by_name(tmp_path):
    release = release_record()
    del release["concordance"]["s02"]
    assert refusal(tmp_path, record_set(release=release)) == f"2001-01-01-release-{RELEASE}.json concordance carries no s02"
    release = release_record()
    release["concordance"]["s01"]["perl"] = "4"
    assert refusal(tmp_path / "b", record_set(release=release)) == f"2001-01-01-release-{RELEASE}.json concordance.s01 carries no count under perl: '4'"
    release = release_record()
    release["per_cell_aggregates"] = [a for a in release["per_cell_aggregates"] if "target_cpu" in a or a["suite_id"] != "s01"]
    assert refusal(tmp_path / "c", record_set(release=release)) == f"2001-01-01-release-{RELEASE}.json has no per_cell_aggregates entry for s01 on arm64"


def test_two_wall_time_aggregates_of_one_cell_and_architecture_are_refused_by_name(tmp_path):
    release = release_record()
    release["per_cell_aggregates"].append({"arch": "x86_64", "suite_id": "s02", "engine": "vep_rs", "median_sec": 5.0})
    assert refusal(tmp_path, record_set(release=release)) == (
        f"2001-01-01-release-{RELEASE}.json has 2 per_cell_aggregates entries for s02 on x86_64; the figure reads one")
    timed = walltime_record()
    timed["per_cell_aggregates"].append({"engine": "perl", "suite_id": "s01", "arch": "arm64", "median_sec": 650.0})
    assert refusal(tmp_path / "b", record_set(walltime=timed)) == (
        f"2001-01-01-comparator-walltime-{RELEASE}.json has 2 per_cell_aggregates entries of perl for s01 on arm64; the figure reads one")


def test_two_wall_time_sums_of_one_dataset_and_architecture_are_refused_for_every_tool(tmp_path):
    for engine, name in (("fastvep", "fastVEP"), ("ensembl_vep", "Ensembl VEP"), ("vep_rs", "vep-rs")):
        population = population_record()
        population["wall_time"][engine]["per_dataset"].append({"dataset": "pop_a", "arch": "arm64", "sum_of_cell_medians_sec": 77.0})
        assert refusal(tmp_path / engine, record_set(population=population)) == (
            f"2001-01-01-population-{RELEASE}.json has 2 wall_time sums of {name} for pop_a on arm64; the figure reads one")
    population = population_record()
    population["wall_time"]["vep_rs"]["per_dataset"] = [s for s in population["wall_time"]["vep_rs"]["per_dataset"] if s["dataset"] != "pop_b"]
    assert refusal(tmp_path / "none", record_set(population=population)) == (
        f"2001-01-01-population-{RELEASE}.json has no wall_time sum of vep-rs for pop_b on arm64")


@pytest.mark.parametrize("entry", [{"intersection": 1, "perl_tuples": 9, "engine_tuples": 9}, {"dataset": ["pop_d"]}, {"dataset": ""}, "pop_d"],
                         ids=["no-dataset", "a-list", "empty", "no-object"])
def test_a_comparator_entry_that_names_no_dataset_is_refused_by_name(tmp_path, entry):
    population = population_record()
    population["comparator_concordance"]["per_dataset"].append(entry)
    assert refusal(tmp_path, record_set(population=population)) == (
        f"2001-01-01-population-{RELEASE}.json comparator_concordance.per_dataset holds an entry that names no dataset")


def test_a_record_whose_release_plugin_or_kind_of_data_is_no_text_is_refused_by_name(tmp_path):
    release = release_record()
    release["release"] = [RELEASE]
    assert refusal(tmp_path / "release", record_set(release=release)) == (
        f"2001-01-01-release-{RELEASE}.json states a release that is no text: ['{RELEASE}']")
    scored = dataset_record()
    scored["scored_with_release"]["version"] = 6
    assert refusal(tmp_path / "scored", record_set(fastvep=scored)) == (
        f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json states a release that is no text: 6")
    plugins = copy.deepcopy(PLUGINS)
    plugins["absent"][0]["plugin"] = ["AncestralAllele"]
    assert refusal(tmp_path / "absent", record_set(plugins=plugins)) == f"2001-01-01-plugins-{RELEASE}.json absent[0] names no plugin and engine"
    plugins = copy.deepcopy(PLUGINS)
    plugins["rows"][3]["scored_on"] = ["clinvar"]
    assert refusal(tmp_path / "rows", record_set(plugins=plugins)) == (
        f"2001-01-01-plugins-{RELEASE}.json rows[3] is scored on ['clinvar']; the figure groups clinvar, committed_table, synthetic_fixture")


def test_a_comparator_run_with_other_extra_arguments_is_refused(tmp_path):
    message = refusal(tmp_path, record_set(fastvep=dataset_record(args=None)))
    assert message.startswith("the records ran fastVEP with different extra arguments (") and message.endswith("); its marks come from one configuration")
    assert f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json: none" in message


def test_one_tool_at_two_versions_is_refused(tmp_path):
    message = refusal(tmp_path, record_set(population=population_record(comparator="v0.0.9")))
    assert message == f"the records name fastVEP at {COMPARATOR_VERSION} and v0.0.9; a tool is drawn at one version, and --draft draws a marked draft"


def ungated(engine="vep-rs", passed=False) -> dict:
    plugins = copy.deepcopy(PLUGINS)
    plugins["gate"].update(engine=engine, passed=passed)
    return plugins


@pytest.mark.parametrize(("engine", "passed"), [("vep-rs", False), (None, False), ("fastvep", True)])
def test_a_plugin_record_that_did_not_pass_its_gate_for_the_released_tool_is_refused(tmp_path, engine, passed):
    message = refusal(tmp_path, record_set(plugins=ungated(engine, passed)))
    assert message == (f"2001-01-01-plugins-{RELEASE}.json did not pass its gate for vep-rs (gate engine {engine!r}, passed {passed!r}); "
                       "--draft draws a marked draft from it")


def mixed_releases() -> dict:
    """The fixture plugin record beside four records of another release."""
    other = "v0.0.1"
    return {
        f"2001-01-01-release-{other}.json": release_record(other),
        f"2001-01-01-population-{other}.json": population_record(other),
        f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json": dataset_record(other),
        f"2001-01-01-comparator-walltime-{other}.json": walltime_record(other),
        "2001-01-01-plugins.json": copy.deepcopy(PLUGINS),
    }


def test_records_that_do_not_all_state_one_release_are_refused(tmp_path):
    message = refusal(tmp_path, mixed_releases())
    assert message.startswith("the records do not all state one release (2001-01-01-release-v0.0.1.json states v0.0.1; ")
    assert message.endswith(f"2001-01-01-plugins.json states {RELEASE}); --draft draws a marked draft from them")
    unstated = dataset_record()
    del unstated["scored_with_release"]
    message = refusal(tmp_path / "b", record_set(fastvep=unstated))
    assert f"2001-01-01-comparator-fastvep-{COMPARATOR_VERSION}.json states none" in message


@pytest.mark.parametrize("records", [record_set(plugins=ungated()), mixed_releases()], ids=["ungated", "mixed"])
def test_draft_draws_an_ungated_record_and_mixed_releases_marked_as_a_draft(tmp_path, records):
    document = values_of(tmp_path, records, "--draft")
    line = "Draft: not drawn from one release's final records."
    assert document["caption"].startswith(f"{line} **vep-rs ") and document["alt_text"].startswith(f"{line} Four panels comparing ")
    assert document["notes"][0] == line and line in figure.drawn_texts(document, figure.read_registry(REGISTRY))
    # No panel of a draft names the release of its records, and no row a version: the legend alone names each
    # tool with the version, or versions, its records state.
    registry = figure.read_registry(REGISTRY)
    assert [p["title"] for p in document["panels"]] == [
        "Concordance per dataset", "Wall time per dataset and architecture", "Concordance per plugin", "Plugins read by vep-rs and fastVEP"]
    assert [row["label"] for group in panel(document, "d")["groups"] for row in group["rows"]] == ["vep-rs", "fastVEP"]
    named = [entry["label"] for entry in figure.legend_entries(document) if entry["mark"] in figure.LEGEND_KINDS[0]]
    assert named == [f"{tool['name']} {tool['version']}" for tool in document["tools"]]
    versions = {version for tool in document["tools"] for version in tool["version"].split(" and ")}
    assert versions >= {RELEASE, COMPARATOR_VERSION, REFERENCE_VERSION}
    drawn = figure.drawn_texts(document, registry)
    assert sorted(text for text in drawn if any(version in text for version in versions) or "records)" in text) == sorted(named)


def test_a_gated_set_of_one_release_carries_no_draft_mark(tmp_path):
    document = values_of(tmp_path)
    assert "Draft" not in json.dumps(document)
    assert [p["title"] for p in document["panels"]] == [
        "Concordance per dataset", "Wall time per dataset and architecture", "Concordance per plugin", "Plugins read by vep-rs and fastVEP"]


@pytest.mark.parametrize("flag", ["--values", "--png"])
def test_a_draft_is_never_written_under_docs_assets(tmp_path, flag):
    directory = write_records(tmp_path / "records", record_set(plugins=ungated()))
    outputs = {"--values": tmp_path / "out" / "v.json", "--png": tmp_path / "out" / "v.png"}
    outputs[flag] = tmp_path / "docs" / "assets" / outputs[flag].name
    done = subprocess.run(
        [sys.executable, str(SCRIPT), "--repo", str(tmp_path), "--records", str(directory), "--draft",
         "--values", str(outputs["--values"]), "--png", str(outputs["--png"])], capture_output=True, text=True, check=False)
    assert done.returncode == 1
    assert done.stderr == f"ERROR: [render_release_figure] --draft writes no file under docs/assets/: {outputs[flag]}\n"
    assert not (tmp_path / "docs").exists() and not (tmp_path / "out").exists()


def test_draft_without_named_outputs_is_refused_since_the_defaults_lie_under_docs_assets(tmp_path):
    directory = write_records(tmp_path / "records", record_set(plugins=ungated()))
    done = subprocess.run([sys.executable, str(SCRIPT), "--repo", str(tmp_path), "--records", str(directory), "--draft", "--values-only"],
                          capture_output=True, text=True, check=False)
    assert done.returncode == 1 and "--draft writes no file under docs/assets/" in done.stderr


# The drawing's checks


class Label:
    """A text artist's surface, as the drawing's checks read it."""

    def __init__(self, content: str, visible: bool = True):
        self.content, self.visible = content, visible

    def get_text(self) -> str:
        return self.content

    def get_visible(self) -> bool:
        return self.visible

    def get_window_extent(self, renderer):
        return (renderer, self.content)


class Tick:
    def __init__(self, location: float, content: str):
        self.location, self.label1, self.label2 = location, Label(content), Label("")

    def get_loc(self) -> float:
        return self.location


class Axis:
    def __init__(self, view: tuple[float, float], major: list[Tick], minor: list[Tick] | None = None):
        self.view, self.major, self.minor = view, major, minor or []

    def get_view_interval(self) -> tuple[float, float]:
        return self.view

    def get_major_ticks(self) -> list[Tick]:
        return self.major

    def get_minor_ticks(self) -> list[Tick]:
        return self.minor


class Axes:
    def __init__(self, xaxis: Axis, yaxis: Axis):
        self.xaxis, self.yaxis = xaxis, yaxis


def test_the_drawing_checks_pass_over_tick_labels_outside_their_axis_view():
    # A reversed logarithmic axis from 1 to 4.5e-8 with decade ticks, one past each end of the view, and a
    # linear axis from 0 to 31 whose locator also made a tick at 35.
    reversed_axis = Axis((1.0, 4.5e-8), [Tick(10.0, "past the start"), Tick(1.0, "0"), Tick(1e-6, "0.999999"), Tick(1e-8, "past the end")])
    linear_axis = Axis((0.0, 31.0), [Tick(0.0, "0"), Tick(30.0, "30"), Tick(35.0, "35")], [Tick(-1.0, "minor past the start")])
    axes = Axes(reversed_axis, linear_axis)
    free = [Label("a title"), Label("   "), Label("hidden", visible=False)]
    texts = [label for axis in (reversed_axis, linear_axis) for tick in (*axis.major, *axis.minor) for label in (tick.label1, tick.label2)] + free
    placed = figure.placed_texts(texts, [axes], "renderer")
    assert [artist.get_text() for artist, _ in placed] == ["0", "0.999999", "0", "30", "a title"]
    assert all(box == ("renderer", artist.get_text()) for artist, box in placed)
    assert figure.hidden_tick_labels([axes]) == {
        id(label) for tick in (reversed_axis.major[0], reversed_axis.major[3], linear_axis.major[2], linear_axis.minor[0])
        for label in (tick.label1, tick.label2)}


class Ruler:
    """A sheet that measures a text at six pixels a character."""

    def measure(self, content: str, **style) -> float:
        return 6.0 * len(content)


def test_a_group_label_stands_in_the_label_column_and_the_plugin_panels_two_plots_share_one_width():
    ruler = Ruler()
    short = [{"label": "Group", "rows": [{"label": "a row label", "values": []}]}]
    long_group = [{"label": "A group label longer than any of its rows", "rows": [{"label": "row", "values": []}]}]
    gap = figure.label_column(ruler, short) - 6.0 * len("a row label")
    assert gap > 0
    # A group label wider than every row label widens the column, so it ends before the plot begins.
    assert figure.label_column(ruler, long_group) == 6.0 * len("A group label longer than any of its rows") + gap
    # The two columns take the whole width, and what each leaves its plot is the same, so the plots share a scale.
    first, second = figure.column_widths(ruler, 1000.0, (short, long_group))
    assert first + second + figure.GUTTER == pytest.approx(1000.0)
    assert first - figure.label_column(ruler, short) == pytest.approx(second - figure.label_column(ruler, long_group))
    assert first < second
    # With the plugin panel in one column, or absent, the columns are equal.
    assert figure.column_widths(ruler, 1000.0, (short, [])) == ((1000.0 - figure.GUTTER) / 2,) * 2
    assert figure.column_widths(ruler, 1000.0, ([], [])) == ((1000.0 - figure.GUTTER) / 2,) * 2


def test_a_plot_narrower_than_the_minimum_is_refused_by_name():
    assert figure.plot_box(16.0, 400.0, 150.0, "A panel") == (166.0, 400.0 - 150.0 - figure.PLOT_PAD)
    assert figure.plot_box(16.0, 150.0 + figure.PLOT_PAD + figure.MIN_PLOT, 150.0, "A panel")[1] == figure.MIN_PLOT
    with pytest.raises(SystemExit) as narrow:
        figure.plot_box(16.0, 250.0, 150.0, "A panel")
    assert str(narrow.value) == (
        f"ERROR: [render_release_figure] the drawn figure leaves the plot of 'A panel' {250 - 150 - figure.PLOT_PAD} pixels beside its "
        f"labels, under the minimum of {figure.MIN_PLOT}")
    # Labels as wide as the column, or wider, leave no room: the refusal states no negative width.
    for width in (150.0, 100.0):
        with pytest.raises(SystemExit) as none:
            figure.plot_box(16.0, width, 150.0, "A panel")
        assert str(none.value) == (
            "ERROR: [render_release_figure] the drawn figure leaves the plot of 'A panel' no room beside its labels, under the minimum of "
            f"{figure.MIN_PLOT}")


def test_wall_time_axis_starts_at_the_decade_under_its_smallest_value_and_names_its_decades_in_short_form():
    # The axis starts at the decade under its smallest value and ends a fifth of a decade past its largest.
    low, high = figure.seconds_limits([1.79, 12.84, 265212.0])
    assert low == 1.0 and high == pytest.approx(265212.0 * 10**0.2) and 265212.0 < high < 10**6
    # Its lowest tick stands at its start, so no value lies before the lowest tick.
    ticks, labels = figure.seconds_ticks(low, high)
    assert ticks == [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0] and labels == ["1", "10", "100", "1K", "10K", "100K"]
    assert ticks[0] == low <= 1.79
    # A value that lies on a decade starts the axis; values under one second start it under one.
    assert figure.seconds_limits([10.0, 40.0])[0] == 10.0 and figure.seconds_limits([99.9, 10000.0])[0] == 10.0
    assert figure.seconds_limits([0.5, 3.0])[0] == 0.1 and figure.seconds_limits([0.00015, 3.0])[0] == 0.0001
    assert [figure.time_label(10.0**exponent) for exponent in range(-2, 8)] == ["0.01", "0.1", "1", "10", "100", "1K", "10K", "100K", "1M", "10M"]


def test_wall_time_axis_has_ticks_however_narrow_its_values():
    # Decades where two or more fall inside the limits.
    assert figure.seconds_ticks(1.13, 420000.0) == ([10.0, 100.0, 1000.0, 10000.0, 100000.0], ["10", "100", "1K", "10K", "100K"])
    # Values within a factor of four: one decade inside, so the 1, 2 and 5 of each decade.
    assert figure.seconds_ticks(*figure.seconds_limits([12.0, 40.0])) == ([10.0, 20.0, 50.0], ["10", "20", "50"])
    assert figure.seconds_ticks(*figure.seconds_limits([1200.0, 4000.0])) == ([1000.0, 2000.0, 5000.0], ["1K", "2K", "5K"])
    assert figure.seconds_ticks(*figure.seconds_limits([21.0, 24.0])) == ([10.0, 20.0], ["10", "20"])
    # Narrower than that, and under the 2 of their decade: the axis's start and the largest value at two figures.
    assert figure.seconds_ticks(*figure.seconds_limits([10.5, 12.0])) == ([10.0, 12.0], ["10", "12"])


def f1_value(tool: str, kind: str, value: str) -> dict:
    """An F1 value whose counts give `value` exactly: so many of a thousand million matched."""
    return {"tool": tool, "kind": kind, "value": value, "counts": [int(Decimal(value) * 10**9), 10**9, 10**9]}


def floor_of(*values: str) -> float:
    """The floor of an F1 axis that holds `values`."""
    return figure.discord_floor([f1_value("vep-rs", "bar", value) for value in values])


def test_the_floor_is_the_first_decade_under_the_smallest_value():
    # No value between 0.9999999 and 1: the floor is the seventh decade, and the ticks before it are six.
    assert floor_of("0.750000000", "0.999998400", "1.000000000") == 1e-7
    assert floor_of("1.000000000") == floor_of() == 1e-7
    assert figure.f1_tick_labels(figure.f1_ticks(1e-7), 1e-7) == ["0", "0.9", "0.99", "0.999", "0.9999", "0.99999", "0.999999", "1"]
    # The smallest 1 - F1 is 3e-7: the first decade under it is the seventh, not a whole decade further.
    assert floor_of("0.999999700", "0.990000000") == 1e-7
    # A value past the seventh decade takes the first decade under it, and one that lies exactly on a decade the
    # decade below, as far as the counts go: nine places hold nothing nearer 1 than the ninth, the counts may.
    assert floor_of("0.999999950", "0.750000000") == floor_of("0.999999985") == 1e-8
    assert floor_of("0.999999900", "1.000000000") == 1e-8
    assert floor_of("0.999999990") == 1e-9
    assert floor_of("0.999999999") == 1e-10
    assert figure.decade_of(Fraction(286, 10**9)) == -7 and figure.decade_of(Fraction(1, 10**7)) == -7 and figure.decade_of(Fraction(999, 10**10)) == -8
    unequal = {"tool": "vep-rs", "kind": "bar", "value": "1.000000000", "counts": [10**10 - 1, 10**10, 10**10]}
    assert figure.discord_floor([unequal]) == 1e-11
    assert figure.f1_ticks(1e-10) == [1.0, 0.1, 0.01, 0.001, 0.0001, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
    assert figure.f1_tick_labels(figure.f1_ticks(1e-10), 1e-10)[-3:] == ["0.99999999", "0.999999999", "1"]
    # Every F1 below 1 is drawn short of the floor and no further from it than the last tick before it; an F1 of 1 at the floor.
    for value in ("0.999999900", "0.999999950", "0.999999999", "0.999998400", "0.999999700", "0.750000000"):
        floor = floor_of(value)
        place = figure.discord(f1_value("vep-rs", "bar", value), floor)
        assert place > floor and place == pytest.approx(1 - float(value), rel=1e-6)
        assert floor == 1e-7 or figure.f1_ticks(floor)[-2] >= place
        assert figure.discord(f1_value("vep-rs", "bar", "1.000000000"), floor) == floor


def test_wall_time_circles_that_would_meet_are_set_apart_along_their_line():
    clear, pitch = 10.6, 44.0
    centre = 1.5 * pitch
    upper, lower = centre - 0.18 * pitch, centre + 0.18 * pitch
    # Circles a decade apart keep their lanes, the two lanes of one time too.
    spread = [(100.0, upper, centre, "A"), (170.0, upper, centre, "B"), (100.0, lower, centre, "A"), (170.0, lower, centre, "B")]
    assert figure.circle_heights(spread, clear, pitch / 2) == [upper, upper, lower, lower]
    # Two tools 6 percent apart at 70 pixels a decade lie under two pixels apart: the second is set clear of
    # the first, away from the line's centre, each still at its own time.
    step = 70.0 * math.log10(1.06)
    assert step < 2.0
    first, second = figure.circle_heights([(100.0, upper, centre, "A"), (100.0 + step, upper, centre, "B")], clear, pitch / 2)
    assert first == upper and second < upper
    assert math.hypot(step, second - first) == pytest.approx(clear)
    # The same on both lanes: four circles in one column, each clear of every other and inside its line.
    column = [(100.0, upper, centre, "A"), (100.0 + step, upper, centre, "B"), (100.0, lower, centre, "A"), (100.0 + step, lower, centre, "B")]
    heights = figure.circle_heights(column, clear, pitch / 2)
    assert heights[0] == upper and heights[2] == lower and heights[1] < upper and heights[3] > lower
    assert all(abs(y - centre) <= pitch / 2 for y in heights)
    for one in range(4):
        for other in range(one):
            assert math.hypot(column[one][0] - column[other][0], heights[one] - heights[other]) >= clear - 1e-6
    # A line too low for the column has no height left for its third circle, which is refused by name.
    low = 30.0
    tight = [(100.0, 1.5 * low - 0.18 * low, 1.5 * low, "A (ARM)"), (100.0, 1.5 * low - 0.18 * low, 1.5 * low, "B (ARM)"),
             (100.0, 1.5 * low + 0.18 * low, 1.5 * low, "A (x86)")]
    with pytest.raises(SystemExit) as refused:
        figure.circle_heights(tight, clear, low / 2)
    assert str(refused.value) == (
        "ERROR: [render_release_figure] the drawn figure has no room on its line to set the circle of A (x86) clear of the circles it meets")


def test_f1_tick_labels_name_the_floor():
    # The floor is where an F1 of exactly 1 is drawn, and its tick says so.
    assert figure.f1_tick_labels([1.0, 0.1, 1e-6, 1e-7]) == ["0", "0.9", "0.999999", "1"]


def test_tick_labels_closer_than_a_labels_height_are_thinned_keeping_the_one_nearer_one():
    # Eight ticks, the decades 36 pixels apart and the floor 34 past the last; labels 7 pixels a character and 10 high.
    labels = ["0", "0.9", "0.99", "0.999", "0.9999", "0.99999", "0.999999", "1"]
    places = [0.0, 36.0, 72.0, 108.0, 144.0, 180.0, 216.0, 250.0]
    widths = [7.0 * len(label) for label in labels]
    shown = figure.fitted_labels(places, widths, 10.0)
    # From the end that reads 1, a label is set when it stands its own height clear of every label set before it.
    assert [label for label, keep in zip(labels, shown) if keep] == ["0", "0.9", "0.999", "0.99999", "1"]
    kept = [index for index, keep in enumerate(shown) if keep]
    for position in range(1, len(kept)):
        left, right = kept[position - 1], kept[position]
        assert (places[right] - widths[right] / 2) - (places[left] + widths[left] / 2) >= 10.0
    # Of two labels 5 pixels apart the one nearer 1 is kept; at a label's height apart both are.
    assert figure.fitted_labels([0.0, 15.0], [10.0, 10.0], 10.0) == [False, True]
    assert figure.fitted_labels([0.0, 20.0], [10.0, 10.0], 10.0) == [True, True]
    # Three in a row, each too close to the next: the one nearest 1 and the one two before it.
    assert figure.fitted_labels([0.0, 15.0, 30.0], [10.0, 10.0, 10.0], 10.0) == [True, False, True]


class Rise(Ruler):
    """A sheet whose glyphs are ten pixels high."""

    def rise(self, content: str, **style) -> float:
        return 10.0


def test_the_tick_that_reads_zero_keeps_its_label_at_every_floor_and_width():
    # A label next to the one that reads 0 and too close to it gives way: the axis says where an F1 of 0 stands.
    assert figure.fitted_labels([0.0, 14.0, 60.0], [7.0, 21.0, 7.0], 10.0) == [True, False, True]
    assert figure.fitted_labels([0.0, 14.0, 60.0, 120.0], [7.0, 21.0, 28.0, 7.0], 10.0) == [True, False, True, True]
    # Two labels alone that would meet: the one that reads 1 is kept. One label alone is set, and none is none.
    assert figure.fitted_labels([0.0, 15.0], [10.0, 10.0], 10.0) == [False, True]
    assert figure.fitted_labels([40.0], [7.0], 10.0) == [True] and figure.fitted_labels([], [], 10.0) == []
    # From the seventh decade to the tenth, on the narrowest plot the canvas admits and on wide ones, both ends
    # of every F1 axis are labelled, and every label set stands its glyphs' height clear of the others.
    sheet = Rise()
    for floor in (1e-7, 1e-8, 1e-9, 1e-10):
        for width in (float(figure.MIN_PLOT), 240.0, 261.0, 330.0):
            labels = figure.f1_axis_labels(sheet, [width], floor)
            assert (labels[0], labels[-1]) == ("0", "1"), (floor, width, labels)
            scale = figure.f1_scale(width, floor)
            spans = [(scale * -math.log10(tick) - 3.0 * len(label), scale * -math.log10(tick) + 3.0 * len(label))
                     for tick, label in zip(figure.f1_ticks(floor), labels) if label]
            assert all(right[0] - left[1] >= 10.0 - 1e-9 for left, right in itertools.pairwise(spans)), (floor, width, labels)


def test_the_f1_panels_label_one_set_of_decades_whatever_their_widths():
    sheet, floor = Rise(), 1e-7
    narrow, wide = 150.0, 330.0
    # A decade takes the plot's width less the room left of F1 = 0, over the decades from there to a little past the floor.
    assert figure.f1_scale(wide, floor) == pytest.approx((wide - figure.f1_room()) / (7 + math.log10(1 / figure.FLOOR_MARGIN)))
    # Each plot alone would label the decades its own width has room for.
    assert figure.f1_axis_labels(sheet, [narrow], floor) == ["0", "", "0.99", "", "", "0.99999", "", "1"]
    assert figure.f1_axis_labels(sheet, [wide], floor) == ["0", "0.9", "0.99", "0.999", "0.9999", "", "0.999999", "1"]
    # Together they carry one set, the one that fits the narrowest, in whatever order the plots are given.
    assert figure.f1_axis_labels(sheet, [narrow, wide], floor) == figure.f1_axis_labels(sheet, [wide, narrow], floor) == [
        "0", "", "0.99", "", "", "0.99999", "", "1"]
    # A label a tick, blank where it is thinned, so every plot draws every tick.
    assert len(figure.f1_axis_labels(sheet, [narrow, wide], 1e-10)) == len(figure.f1_ticks(1e-10))


def test_a_tool_without_a_bar_has_a_circle_in_the_legend(tmp_path):
    entries = figure.legend_entries(values_of(tmp_path))
    assert [(entry["mark"], entry["tool"]) for entry in entries[:3]] == [("fill_dot", "vep-rs"), ("fill_dot", "fastvep"), ("dot", "ensembl-vep")]


class Handle:
    """A legend handle of the drawing library, as the legend's marks are made: its arguments and its style."""

    def __init__(self, *args, **style):
        self.args, self.style = args, style


def test_a_circle_in_the_tools_colour_stands_beside_each_tools_swatch_in_the_legend(tmp_path, monkeypatch):
    document = values_of(tmp_path)
    # A tool with bars in an F1 or count panel and circles in the wall-time panel shows both in the legend; the
    # reference, which has circles alone, shows its circle.
    entries = figure.legend_entries(document)
    assert [(entry["mark"], entry["tool"], entry["label"]) for entry in entries[:3]] == [
        ("fill_dot", "vep-rs", f"vep-rs {RELEASE}"), ("fill_dot", "fastvep", f"fastVEP {COMPARATOR_VERSION}"),
        ("dot", "ensembl-vep", f"Ensembl VEP {REFERENCE_VERSION}")]
    assert [entry["mark"] for entry in figure.legend_rows(entries)[0]] == ["fill_dot", "fill_dot", "dot"]
    assert [entry["mark"] for entry in entries[3:]] == ["raw", "exact", "adjusted", "absent", "filled", "open"]
    # A tool kept out of the wall-time panel has no circle there, and none in the legend.
    unclocked = values_of(tmp_path / "unclocked", registry=registry_with(tmp_path, lambda rows: rows["fastvep"].update(speed_panel="no")))
    assert [(entry["mark"], entry["tool"]) for entry in figure.legend_entries(unclocked)[:3]] == [
        ("fill_dot", "vep-rs"), ("fill", "fastvep"), ("dot", "ensembl-vep")]
    # The panel of wall times holds eight circles a tool, no name, no other kind of value.
    wall = [value for group in panel(document, "b")["groups"] for row in group["rows"] for value in row["values"]]
    assert {value["kind"] for value in wall} == {"marker"} and len(wall) == 24
    # The pair is a swatch and a circle, both in the tool's colour, the circle as the wall-time panel draws a
    # filled one; on such a row the reference's circle has a blank where a swatch would stand.
    lines, patches = type(sys)("matplotlib.lines"), type(sys)("matplotlib.patches")
    lines.Line2D, patches.Patch = type("Line2D", (Handle,), {}), type("Patch", (Handle,), {})
    for name, module in (("matplotlib", type(sys)("matplotlib")), ("matplotlib.lines", lines), ("matplotlib.patches", patches)):
        monkeypatch.setitem(sys.modules, name, module)
    colour = FILLS["fastvep"]
    swatch, circle = figure.legend_handle(entries[1], colour)
    assert type(swatch).__name__ == "Patch" and swatch.style == {"facecolor": colour, "edgecolor": "none"}
    assert type(circle).__name__ == "Line2D" and circle.style["marker"] == "o" and circle.style["markersize"] == figure.CIRCLE
    assert circle.style["markerfacecolor"] == circle.style["markeredgecolor"] == colour and circle.style["markeredgewidth"] == figure.CIRCLE_EDGE
    alone = figure.legend_handle(entries[2], "#009E73")
    assert type(alone).__name__ == "Line2D" and alone.style["markerfacecolor"] == "#009E73"
    assert type(figure.legend_handle({"mark": "fill", "tool": "fastvep", "label": "a tool"}, colour)).__name__ == "Patch"
    assert figure.LEGEND_PAIR > 2 * figure.LEGEND_HANDLE
    # The row is set as runs of entries whose marks take one width, side by side: the two pairs, then the
    # reference's circle with a single mark's width, so no blank stands before its circle.
    runs = figure.legend_runs(figure.legend_rows(entries)[0])
    assert [[entry["mark"] for entry in run] for run in runs] == [["fill_dot", "fill_dot"], ["dot"]]
    assert [[entry["tool"] for entry in run] for run in runs] == [["vep-rs", "fastvep"], ["ensembl-vep"]]
    assert figure.legend_runs(figure.legend_rows(entries)[1]) == [figure.legend_rows(entries)[1]]
    mixed = [{"mark": mark} for mark in ("dot", "fill_dot", "fill", "fill", "fill_dot")]
    assert [[entry["mark"] for entry in run] for run in figure.legend_runs(mixed)] == [["dot"], ["fill_dot"], ["fill", "fill"], ["fill_dot"]]
    assert figure.legend_runs([]) == []
    # Each run is moved to its place in the figure's own coordinates, as the legend's frame is set, so the saved
    # window carries the runs and their frame together.
    assert ("legend.set_bbox_to_anchor((left / sheet.width, 1 - below / sheet.height), transform=sheet.figure.transFigure)"
            in SCRIPT.read_text(encoding="utf-8"))


class Inked(Ruler):
    """A sheet that records each text it sets with its colour; a text without one is set in the first ink."""

    def __init__(self):
        self.placed = {}

    def text(self, x: float, y: float, content: str, **style):
        self.placed[content] = style.get("color", figure.INK)


class Box:
    def __init__(self, x0: float, y0: float, x1: float, y1: float):
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1


# A grey under the contrast of text against the white ground.
FAINT = "#898781"


def test_the_word_synthetic_and_every_text_meet_the_contrast_of_text():
    # The ratio of two relative luminances, each plus 0.05: 21 for black on white, 1 for a colour on itself.
    assert figure.ink_contrast("#000000") == pytest.approx(21.0) and figure.ink_contrast("#ffffff") == pytest.approx(1.0)
    assert figure.ink_contrast("#767676") == pytest.approx(4.54, abs=0.005) and figure.ink_contrast("#777777") < 4.5
    assert figure.ink_contrast("#ffffff", "#0072B2") == figure.ink_contrast("#0072B2") == pytest.approx(5.19, abs=0.01)
    # The two inks text is set in are over the contrast of text; a paler grey is under it.
    assert figure.MIN_CONTRAST == 4.5 and figure.ink_contrast(FAINT) == pytest.approx(3.6, abs=0.05)
    assert figure.ink_contrast(figure.INK) > figure.ink_contrast(figure.INK_SECONDARY) > 7.0
    # A synthetic row's label is set as its name and the word synthetic, each in an ink over the contrast of text.
    sheet = Inked()
    groups = [{"label": "Synthetic fixtures", "rows": [{"label": "CADD GRCh38 synthetic", "values": []}, {"label": "REVEL GRCh37", "values": []}]}]
    figure.draw_rows(sheet, Recorder(), 16.0, 200.0, 100.0, 20.0, groups)
    assert set(sheet.placed) == {"Synthetic fixtures", "CADD GRCh38", "synthetic", "REVEL GRCh37"}
    assert all(figure.ink_contrast(colour) >= 4.5 for colour in sheet.placed.values())
    assert sheet.placed["synthetic"] != sheet.placed["CADD GRCh38"] == figure.INK
    # A drawing with a text under that contrast is refused by name, as one with small type or a text over another is.
    clear = [("a title", 8.0, figure.INK, Box(0, 0, 40, 10)), ("synthetic", 7.0, figure.INK_SECONDARY, Box(50, 0, 90, 10))]
    figure.check_texts(clear)
    for texts, message in (
        ([*clear, ("synthetic", 7.0, FAINT, Box(100, 0, 140, 10))],
         "the drawn figure sets 'synthetic' at a contrast of 3.6 to 1 against its ground, under the minimum of 4.5"),
        ([*clear, ("a note", 6.5, figure.INK, Box(100, 0, 140, 10))], "the drawn figure sets 'a note' at 6.5 points, under the minimum of 7"),
        ([*clear, ("a label", 7.0, figure.INK, Box(85, 5, 120, 15))], "the drawn figure sets 'synthetic' over 'a label'"),
    ):
        with pytest.raises(SystemExit) as refused:
            figure.check_texts(texts)
        assert str(refused.value) == f"ERROR: [render_release_figure] {message}"
    # Texts that touch by a pixel are not one over the other.
    figure.check_texts([*clear, ("a label", 7.0, figure.INK, Box(89.5, 0, 120, 10))])


def test_every_colour_a_text_is_set_in_meets_the_contrast_of_text_against_its_ground():
    # The bound is measured on the colours the script sets, a text's and its ground's: each ink against the ground.
    assert figure.TEXT_INKS == (figure.INK, figure.INK_SECONDARY) and figure.MIN_CONTRAST == 4.5
    for ink in figure.TEXT_INKS:
        assert figure.ink_contrast(ink, figure.GROUND) >= figure.MIN_CONTRAST, ink
    # Every call that sets a text names its colour as one of those inks, never a literal or a tool's colour.
    source = SCRIPT.read_text(encoding="utf-8")
    colours = [
        keyword.value for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr in ("text", "annotate", "legend", "wrap", "measure", "rise")
        for keyword in node.keywords if keyword.arg in ("color", "labelcolor")
    ]
    assert len(colours) >= 6 and all(isinstance(colour, ast.Name) and colour.id in ("INK", "INK_SECONDARY") for colour in colours)
    # A text set without a colour takes the first ink: the sheet's default, and the defaults of the drawing library.
    assert 'style.setdefault("color", INK)' in source
    assert all(figure.STYLE[key] == figure.INK for key in ("text.color", "axes.labelcolor", "xtick.labelcolor", "ytick.labelcolor"))
    # Every text stands on the one ground, which fills the canvas and every plot; no fill is named apart from it.
    assert figure.STYLE["figure.facecolor"] == figure.STYLE["axes.facecolor"] == figure.GROUND and '"white"' not in source
    # The bound has teeth: a tool's colour would fail it as an ink.
    assert figure.ink_contrast(figure.FILLS["2"], figure.GROUND) < figure.MIN_CONTRAST


def test_the_legend_is_grouped_by_kind(tmp_path):
    # A release with an adjusted F1 below 1 on one dataset draws both adjusted entries.
    release = release_record()
    release["concordance"]["s01"]["adj_perl"] = 4
    entries = figure.legend_entries(values_of(tmp_path, record_set(release=release)))
    rows = figure.legend_rows(entries)
    # The tools on one row, the marks of an F1 on one row with the two adjusted entries side by side, the
    # architectures on one.
    assert [[entry["mark"] for entry in row] for row in rows] == [
        ["fill_dot", "fill_dot", "dot"], ["raw", "exact", "adjusted", "absent"], ["filled", "open"]]
    assert [entry["tool"] for entry in rows[0]] == ["vep-rs", "fastvep", "ensembl-vep"]
    assert [entry["label"] for entry in rows[1][:3]] == ["Bar: raw F1", "Adjusted F1, exactly 1", "Adjusted F1, below 1"]
    assert [entry["label"] for entry in rows[2]] == ["ARM", "x86"]
    # Every entry stands in one row, and a kind without an entry has no row.
    assert sorted(entry["label"] for row in rows for entry in row) == sorted(entry["label"] for entry in entries)
    assert figure.legend_rows([entry for entry in entries if entry["mark"] not in ("raw", "exact", "adjusted", "absent")]) == [
        rows[0], rows[2]]
    assert figure.legend_rows([]) == []


class Page(Ruler):
    """A sheet that records where each text is set."""

    width = 1080

    def __init__(self):
        self.placed = []

    def text(self, x: float, y: float, content: str, **style):
        self.placed.append((content, x))


def test_an_axis_label_is_centred_on_its_plot_and_kept_on_the_canvas():
    page = Page()
    line = "a label wider than the plot it stands under"
    figure.draw_axis_label(page, 500.0, 0.0, [line])
    figure.draw_axis_label(page, 20.0, 0.0, [line])
    figure.draw_axis_label(page, 1070.0, 0.0, [line])
    half = 6.0 * len(line) / 2
    assert [x for _, x in page.placed] == [500.0, figure.MARGIN / 2 + half, 1080 - figure.MARGIN / 2 - half]


def test_the_plugin_panels_second_part_carries_no_second_letter_and_is_titled_as_the_first_parts_continued(tmp_path):
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    plugins = panel(document, "c")
    # The first part, or the only one, carries the panel's letter and title; a part that continues it carries
    # the title as continued and no letter.
    assert figure.part_heading(plugins, False) == ("c", "Concordance per plugin")
    assert figure.part_heading(plugins, True) == (None, "Concordance per plugin, continued")
    assert figure.part_heading(panel(document, "a"), False) == ("a", "Concordance per dataset")
    # Drawn, the first is a letter and the title beside it, the second the title alone, from the column's left edge.
    page = Page()
    figure.draw_title(page, 16.0, 100.0, *figure.part_heading(plugins, False))
    figure.draw_title(page, 544.0, 300.0, *figure.part_heading(plugins, True))
    assert [content for content, _ in page.placed] == ["(c)", "Concordance per plugin", "Concordance per plugin, continued"]
    assert [x for _, x in page.placed] == [16.0, 16.0 + 6.0 * len("(c)") + 6, 544.0]
    # The drawing calls for the continued heading over the second part alone, and refuses a figure that sets a
    # panel's letter twice or not at all.
    source = SCRIPT.read_text(encoding="utf-8")
    assert 'heading = part_heading(panels["c"], bool(index and parts[0]))' in source and '", continued"' in source
    assert source.count("CONTINUED") == 2 and source.count('", continued"') == 1
    letters = [f"({p['id']})" for p in document["panels"]]
    drawn = figure.drawn_texts(document, registry)
    figure.check_letters([*drawn, "Concordance per plugin, continued"], letters)
    for strings, times in (([*drawn, "(c)"], 2), ([text for text in drawn if text != "(c)"], 0)):
        with pytest.raises(SystemExit) as refused:
            figure.check_letters(strings, letters)
        assert str(refused.value) == f"ERROR: [render_release_figure] the drawn figure sets the letter (c) {times} times; a panel carries its letter once"
    # The caption and the alternative text follow it: each names the panel once, by its one letter, the
    # alternative text counts four panels, and neither speaks of a second panel of plugins.
    assert [drawn.count(letter) for letter in letters] == [1, 1, 1, 1] and "Concordance per plugin, continued" not in drawn
    for text in (document["caption"], document["alt_text"]):
        assert [text.count(f"{letter} ") for letter in letters] == [1, 1, 1, 1] and "continued" not in text
    assert document["alt_text"].startswith("Four panels comparing ") and len(document["panels"]) == 4


class Recorder:
    """An axes that records what is asked of it: each call's name and arguments, and its keywords beside them."""

    def __init__(self):
        self.calls, self.styles = [], []

    def __getattr__(self, name):
        def record(*args, **kwargs):
            self.calls.append((name, args))
            self.styles.append(kwargs)

        return record

    def drawn(self, name: str) -> list[tuple[tuple, dict]]:
        """The arguments and keywords of every call of `name`, in order."""
        return [(args, style) for (called, args), style in zip(self.calls, self.styles) if called == name]


def f1_row(*values: tuple[str, str, str]) -> dict:
    """An F1 row of (tool, kind, value) values."""
    return {"label": "a row", "values": [f1_value(tool, kind, value) for tool, kind, value in values]}


LANES = {"vep-rs": -0.22, "fastvep": 0.22}
FILLS = {"vep-rs": "#0072B2", "fastvep": "#E69F00"}


def joins_of(axes) -> list:
    """The joins an axes was asked to draw, as (arguments, keywords): the lines that set a cap for a dash."""
    return [(args, style) for args, style in axes.drawn("plot") if not style.get("marker") and "dash_capstyle" in style]


def feet_of(axes) -> list:
    """The solid stretches an axes was asked to draw on into a marker's body, past a join's end."""
    return [(args, style) for args, style in axes.drawn("plot") if not style.get("marker") and "dash_capstyle" not in style]


def drawn_join(tip: float, place: float, shape: str, scale: float) -> list:
    """The two ends a dashed line is drawn between, to compare with a recorded call."""
    return [pytest.approx(value) for value in figure.join_drawn(tip, place, shape, scale)]


def placed(tool: str, kind: str, distance: float) -> dict:
    """An F1 value whose counts give a 1 - F1 of `distance` to fifteen places, and an F1 of exactly 1 for 0."""
    whole = 10**15
    return {"tool": tool, "kind": kind, "value": "", "counts": [whole - round(distance * whole), whole, whole]}


def test_a_marker_over_a_bars_tip_is_drawn_whole_and_on_top_with_an_outline_in_the_grounds_colour():
    # Raw a little over two pixels short of adjusted on the marker's own lane, and a comparator's tip under the marker one lane down.
    row = f1_row(("vep-rs", "bar", "0.999993000"), ("vep-rs", "marker", "0.999994000"), ("fastvep", "bar", "0.999994000"))
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, 35.0)
    # Two bars and the marker with its line, and nothing more: no gap is cut into the marker and no bar is drawn
    # again over it.
    assert [name for name, _ in axes.calls] == ["barh", "plot", "barh"]
    bars, ((marker_at, marker),) = axes.drawn("barh"), [(args, style) for args, style in axes.drawn("plot") if style.get("marker")]
    assert [style["color"] for _, style in bars] == [FILLS["vep-rs"], FILLS["fastvep"]]
    assert marker_at == ([pytest.approx(1 - 0.999994, rel=1e-6)], [-0.22])
    # The marker lies over every bar, whole, in its tool's colour, with a thin outline in the ground's colour.
    assert marker["zorder"] > max(style["zorder"] for _, style in bars)
    assert marker["markerfacecolor"] == FILLS["vep-rs"] and marker["markeredgecolor"] == figure.GROUND == "#ffffff"
    assert 0 < marker["markeredgewidth"] <= 0.6 and marker["clip_on"] is False
    # A marker of exactly 1 over the full bar of another tool is drawn the same way.
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, f1_row(("vep-rs", "bar", "0.900000000"), ("vep-rs", "marker", "1.000000000"),
                                           ("fastvep", "bar", "1.000000000")))], LANES, FILLS, 0.3, 1e-7, 35.0)
    assert [name for name, _ in axes.calls] == ["barh", "plot", "plot", "barh"]
    assert axes.styles[2]["markeredgecolor"] == figure.GROUND and axes.styles[2]["zorder"] > axes.styles[3]["zorder"]


def pixels_from_the_edge(place: float, scale: float) -> float:
    """The pixels a place on an F1 axis lies right of the plot's left edge, at `scale` pixels a decade."""
    return scale * math.log10(figure.f1_edge(scale) / place)


def test_every_f1_plots_left_edge_stands_a_fixed_room_left_of_the_tick_that_reads_zero_and_every_bar_starts_there():
    # The room is the widest marker's half width, its outline included, and four pixels more: the diamond's
    # half diagonal, which is wider than the triangle from its centre to its tip.
    reaches = {shape: figure.marker_reach(shape) for shape in figure.MARKERS}
    diamond = figure.MARKERS["adjusted"][1]
    assert reaches["adjusted"] == pytest.approx((diamond + figure.MARKER_EDGE) * math.sqrt(0.5) * figure.DPI / 72)
    assert reaches["adjusted"] > reaches["exact"] > 0
    assert figure.ROOM_CLEAR == 4.0 and figure.f1_room() == pytest.approx(reaches["adjusted"] + 4.0)
    # The room is the same number of pixels whatever the plot's width and the floor: the tick that reads 0 stands
    # that far right of the left edge, and the plot's width runs from the edge to a little past the floor.
    for width, floor in ((float(figure.MIN_PLOT), 1e-7), (261.0, 1e-7), (330.0, 1e-10)):
        scale = figure.f1_scale(width, floor)
        assert figure.f1_ticks(floor)[0] == 1.0 < figure.f1_edge(scale)
        assert pixels_from_the_edge(1.0, scale) == pytest.approx(figure.f1_room())
        assert pixels_from_the_edge(floor * figure.FLOOR_MARGIN, scale) == pytest.approx(width)
    # The axis starts at that edge, and every bar starts there too, on every row: at a figure of 0, at a full bar.
    assert "axes.set_xlim(f1_edge(scale), floor * FLOOR_MARGIN)" in SCRIPT.read_text(encoding="utf-8")
    scale = figure.f1_scale(261.0, 1e-7)
    rows = [(0, f1_row(("vep-rs", "bar", "0.990000000"), ("fastvep", "bar", "0.000000000"))),
            (1, f1_row(("vep-rs", "bar", "1.000000000"), ("fastvep", "bar", "0.500000000")))]
    axes = Recorder()
    figure.draw_f1_marks(axes, rows, LANES, FILLS, 0.3, 1e-7, scale)
    bars = axes.drawn("barh")
    assert len(bars) == 4 and all(style["left"] == figure.f1_edge(scale) for _, style in bars)
    # A bar of an F1 of 0 is the room long: on these axes a value is a place, not a length.
    assert pixels_from_the_edge(bars[1][1]["left"] + bars[1][0][1], scale) == pytest.approx(figure.f1_room())


def test_a_marker_stands_whole_inside_the_plot_and_a_bar_is_seen_left_of_its_own_marker_on_every_row():
    # From an F1 of 0 to one near a quarter, as a comparator's on the structural-variant rows, on the narrowest plot
    # the canvas admits and on a wide one: the marker's place is its counts', and the bar ends there.
    for scale in (figure.f1_scale(float(figure.MIN_PLOT), 1e-7), figure.f1_scale(261.0, 1e-7)):
        for value in ("0.000000000", "0.010000000", "0.230000000", "0.250000000", "0.500000000"):
            row = f1_row(("fastvep", "bar", value), ("fastvep", "marker", value))
            axes = Recorder()
            figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
            ((bar_at, bar),) = axes.drawn("barh")
            ((marker_at, marker),) = [(args, style) for args, style in axes.drawn("plot") if style.get("marker")]
            shape = figure.marker_shape(row["values"][1])
            assert (marker["marker"], marker["markersize"], marker["markeredgewidth"]) == (*figure.MARKERS[shape], figure.MARKER_EDGE)
            start, tip = pixels_from_the_edge(bar["left"], scale), pixels_from_the_edge(bar["left"] + bar_at[1], scale)
            centre, reach = pixels_from_the_edge(marker_at[0][0], scale), figure.marker_reach(shape)
            assert centre == pytest.approx(scale * -math.log10(1 - float(value)) + figure.f1_room()) == pytest.approx(tip)
            # The marker stands whole inside the plot, four pixels or more clear of the axis line at the left edge.
            assert centre - reach >= 4.0 - 1e-9
            # The bar starts at the edge, so four pixels or more of it are seen left of its own marker.
            assert start == pytest.approx(0.0, abs=1e-9) and (centre - reach) - start >= 4.0 - 1e-9


def test_a_bar_ends_at_its_counts_place_and_no_least_length_draws_it_longer():
    scale = 35.0
    # Figures near 0 lie close together on this axis, and each bar still ends at its own place: none is drawn out
    # to a least length, which would draw different figures alike.
    figures = ("0.010000000", "0.060000000", "0.120000000", "0.150000000", "0.500000000")
    rows = [(index, f1_row(("vep-rs", "bar", "0.990000000"), ("fastvep", "bar", value))) for index, value in enumerate(figures)]
    axes = Recorder()
    figure.draw_f1_marks(axes, rows, LANES, FILLS, 0.3, 1e-7, scale)
    theirs = [(args, style) for args, style in axes.drawn("barh") if style["color"] == FILLS["fastvep"]]
    ends = [style["left"] + args[1] for args, style in theirs]
    assert ends == [pytest.approx(1 - float(value)) for value in figures]
    assert ends == sorted(ends, reverse=True) and len({round(end, 9) for end in ends}) == len(figures)
    assert not hasattr(figure, "MIN_BAR") and "least" not in figure.draw_f1_marks.__doc__
    # The dashed line to a marker starts at the end of the bar, its counts' place.
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, f1_row(("fastvep", "bar", figures[0]), ("fastvep", "marker", "0.500000000")))], LANES, FILLS, 0.3, 1e-7, scale)
    ((join_at, _),) = joins_of(axes)
    tip = 1 - float(figures[0])
    assert join_at[0] == drawn_join(tip, 0.5, "adjusted", scale) and join_at[0][0] == pytest.approx(tip)


def test_a_dashed_line_is_drawn_wherever_a_marker_and_its_bars_end_do_not_touch(tmp_path):
    scale = 35.0
    apart = f1_row(("vep-rs", "bar", "0.990000000"), ("vep-rs", "marker", "1.000000000"), ("fastvep", "bar", "0.900000000"))
    below = f1_row(("vep-rs", "bar", "0.990000000"), ("vep-rs", "marker", "0.999000000"))
    full = f1_row(("vep-rs", "bar", "1.000000000"), ("vep-rs", "marker", "1.000000000"))
    lone = f1_row(("vep-rs", "marker", "1.000000000"), ("fastvep", "bar", "0.900000000"))
    # A marker 8.5 pixels from its bar's end, three pixels past its body's reach, and one 4 pixels from it,
    # whose body stands a pixel and more over the bar's end.
    near = f1_row(("fastvep", "bar", "0.998250000"), ("fastvep", "marker", "0.999000000"))
    touching = f1_row(("fastvep", "bar", "0.998700000"), ("fastvep", "marker", "0.999000000"))
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, apart), (1, below), (2, full), (3, lone), (4, near), (5, touching)], LANES, FILLS, 0.3, 1e-7, scale)
    joins = joins_of(axes)
    reach = {shape: figure.marker_body(shape, 0) for shape in figure.MARKERS}
    # One line for each marker that does not touch its own tool's bar's end, along the tool's lane, from the bar's
    # end to the marker's body and a pixel into it: to a marker of exactly 1 as to one below 1, for a comparator as
    # for the released tool. A marker at the end of a full bar, one whose body stands over its bar's end and one
    # without a bar have none.
    assert [(args, style["color"]) for args, style in joins] == [
        (([pytest.approx(0.01), pytest.approx(1e-7 * 10 ** ((reach["exact"] - figure.JOIN_LAP) / scale))], [-0.22, -0.22]), FILLS["vep-rs"]),
        (([pytest.approx(0.01), pytest.approx(0.001 * 10 ** ((reach["adjusted"] - figure.JOIN_LAP) / scale))], [1 - 0.22, 1 - 0.22]), FILLS["vep-rs"]),
        (([pytest.approx(0.00175), pytest.approx(0.001 * 10 ** ((reach["adjusted"] - figure.JOIN_LAP) / scale))], [4 + 0.22, 4 + 0.22]), FILLS["fastvep"])]
    assert figure.join_end(0.0013, 0.001, "adjusted", scale) is None and figure.join_end(1e-7, 1e-7, "exact", scale) is None
    # The stretch three pixels long carries ink from end to end: a solid stub, which starts at the bar's end.
    assert 3.0 < scale * math.log10(0.00175 / figure.join_end(0.00175, 0.001, "adjusted", scale)) < 3.1
    assert joins[2][1]["linestyle"] == "solid" and joins[2][1]["solid_capstyle"] == "butt"
    for _, style in joins[:2]:
        offset, (dash, gap) = style["linestyle"]
        assert offset == 0 and dash >= 1.5 and gap > 0 and style["dash_capstyle"] == "butt" and 0 < style["linewidth"] <= 1.0
        # A dash with its gap takes under five pixels, so a line of a few pixels reads as dashed; the line keeps a
        # pixel and a half of ink at the README's width, where the image is drawn at 830 pixels of 1080.
        assert (dash + gap) * style["linewidth"] * figure.DPI / 72 < 5.0 and gap * style["linewidth"] * figure.DPI / 72 >= 1.5
        assert style["linewidth"] * figure.DPI / 72 * 830 / 1080 >= 1.4
        # It is drawn over the marker's outline, which it crosses to reach the marker's body.
        assert style["zorder"] > figure.MARKER_ORDER > figure.BAR_ORDER
    # A marker that stands left of its bar's end has no line.
    assert figure.join_end(0.001, 0.01, "adjusted", scale) is None
    assert len([style for _, style in axes.drawn("plot") if style.get("marker")]) == 6
    assert "a marker is its tool's adjusted F1, joined to the end of its bar by a line where the two differ." in values_of(tmp_path)["caption"]


def test_a_marker_touches_its_bar_where_its_body_inside_its_outline_stands_over_the_bars_end_and_the_line_runs_to_the_body():
    scale = 35.0
    edge = figure.MARKER_EDGE * figure.DPI / 72
    shape_reach = {shape: figure.MARKER_SPANS[figure.MARKERS[shape][0]] * figure.MARKERS[shape][1] * figure.DPI / 72 for shape in figure.MARKERS}
    # A marker's body is what its outline leaves of it: on a flat side the outline covers half its width, at a
    # right corner the mitre's reach inward, and at the triangle's tip the reach of a sharper mitre.
    assert figure.marker_body("exact", 0) == pytest.approx(shape_reach["exact"] - edge / 2)
    assert figure.marker_body("exact", 1) == pytest.approx(shape_reach["exact"] - edge / 2 / math.sin(math.atan(0.5)))
    for side in (0, 1):
        assert figure.marker_body("adjusted", side) == pytest.approx(shape_reach["adjusted"] - edge / 2 / math.sin(math.radians(45)))
    for shape in figure.MARKERS:
        assert 0 < figure.marker_body(shape, 0) < shape_reach[shape] < figure.marker_reach(shape)
    # A diamond 6.2 pixels from its bar's end: its outline reaches the bar's end and its body does not. The two do
    # not touch, and the line is drawn: from the bar's end through the outline to the body.
    body = figure.marker_body("adjusted", 0)
    apart = scale * math.log10(0.0015 / 0.001)
    assert body < apart < shape_reach["adjusted"]
    row = f1_row(("fastvep", "bar", "0.998500000"), ("fastvep", "marker", "0.999000000"))
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
    ((join_at, join),) = joins_of(axes)
    ((_, marker),) = [(args, style) for args, style in axes.drawn("plot") if style.get("marker")]
    end = figure.join_end(0.0015, 0.001, "adjusted", scale)
    assert scale * math.log10(end / 0.001) == pytest.approx(body)
    # The outline lies over the bar's end, so the line is drawn from a pixel before the outline, over the bar and
    # in the bar's colour: it crosses the whole outline, and no ground is left between the bar and the body.
    ((_, bar),) = axes.drawn("barh")
    assert scale * math.log10(join_at[0][0] / 0.001) == pytest.approx(figure.marker_reach("adjusted") + figure.JOIN_LAP)
    assert figure.marker_reach("adjusted") + figure.JOIN_LAP > apart
    assert join["color"] == bar["color"] and join["linestyle"] == "solid"
    # It is drawn over the marker and its outline, in the colour of the body, and on into the body: a marker and a
    # line's ends are each set on whole pixels, up to half a pixel from their places, so a pixel's lap meets the body.
    assert join["zorder"] > marker["zorder"] and join["color"] == marker["markerfacecolor"] != marker["markeredgecolor"]
    assert scale * math.log10(end / join_at[0][1]) == pytest.approx(figure.JOIN_LAP) and 1.0 <= figure.JOIN_LAP < min(
        figure.marker_body(shape, side) for shape in figure.MARKERS for side in (0, 1))
    assert join_at[0] == drawn_join(0.0015, 0.001, "adjusted", scale) and 0.001 < join_at[0][1] < end < 0.0015 < join_at[0][0]
    # The two touch where the body stands half a pixel or more over the bar's end, and only there: a hundredth of
    # a pixel short of that the line is drawn, however short, for the diamond and for the triangle, and it ends
    # at the body's edge.
    for shape, place in (("adjusted", 0.001), ("exact", 1e-7)):
        reach = figure.marker_body(shape, 0)
        assert figure.join_end(place * 10 ** ((reach - figure.JOIN_TOUCH - 0.01) / scale), place, shape, scale) is None
        for stretch in (-figure.JOIN_TOUCH + 0.01, 0.0, 0.01):
            short = figure.join_end(place * 10 ** ((reach + stretch) / scale), place, shape, scale)
            assert short == pytest.approx(place * 10 ** (reach / scale)), (shape, stretch)
    # A marker left of its bar's end is not joined, and one that touches its bar's end is not.
    assert figure.join_end(0.001, 0.01, "exact", scale) is None and figure.join_drawn(0.001, 0.01, "exact", scale) is None
    assert figure.join_drawn(0.0013, 0.001, "adjusted", scale) is None


def test_the_line_is_drawn_unless_the_markers_body_stands_half_a_pixel_or_more_over_its_bars_end():
    scale = 35.0
    assert figure.JOIN_TOUCH == 0.5 and figure.JOIN_LAP - figure.JOIN_TOUCH >= 0.5

    def joins(tool: str, shape: str, place: float, over: float) -> list:
        """The lines drawn for a marker whose body stands `over` pixels over its bar's end: each line's two ends
        in pixels from the bar's end, towards the marker, and its style."""
        tip = place * 10 ** ((figure.marker_body(shape, 0) - over) / scale)
        row = {"label": "a row", "values": [placed(tool, "bar", tip), placed(tool, "marker", 0 if shape == "exact" else place)]}
        axes = Recorder()
        figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
        return [([scale * math.log10(tip / end) for end in args[0]], style) for args, style in joins_of(axes)]

    for tool, (shape, place) in itertools.product(("vep-rs", "fastvep"), (("adjusted", 0.001), ("exact", 1e-7))):
        body = figure.marker_body(shape, 0)
        outline = figure.marker_reach(shape) - body
        # Half a pixel or more over the bar's end, up to a marker that stands on it: the two touch, and no line.
        for over in (0.51, 1.0, 3.0, body):
            assert joins(tool, shape, place, over) == [], (tool, shape, over)
        # Under half a pixel over it the line is drawn, a solid stub in the tool's colour. A marker is set on a
        # whole pixel, up to half a pixel from its place: the stub starts over the bar, half a pixel and more
        # before the outline, and ends half a pixel and more inside the body, so no ground shows before the body.
        for over in (0.49, 0.25, 0.01):
            (((start, stop), style),) = joins(tool, shape, place, over)
            assert style["linestyle"] == "solid" and style["color"] == FILLS[tool], (tool, shape, over)
            assert start <= -(over + outline) - 0.5 and stop >= 0.5 - over, (tool, shape, over)
            assert start == pytest.approx(-(over + outline + figure.JOIN_LAP), abs=1e-6) and stop == pytest.approx(figure.JOIN_LAP - over, abs=1e-6)


def test_a_marker_left_of_its_bars_end_stands_alone_without_a_line():
    scale = 35.0
    place = 0.01
    body = figure.marker_body("adjusted", 0)
    for tool in ("vep-rs", "fastvep"):
        # An adjusted F1 under the raw one: the diamond's place lies left of its bar's end, from a fraction of a
        # pixel to the plot's width. No line is drawn at any distance: a line there would lie over the bar.
        for left in (0.3, 3.0, body, body + 0.01, 6.5, 35.0, 200.0):
            tip = place * 10 ** (-left / scale)
            assert figure.join_end(tip, place, "adjusted", scale) is None and figure.join_drawn(tip, place, "adjusted", scale) is None, left
            row = {"label": "a row", "values": [placed(tool, "bar", tip), placed(tool, "marker", place)]}
            axes = Recorder()
            figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
            assert [name for name, _ in axes.calls] == ["barh", "plot"] and joins_of(axes) == feet_of(axes) == [], (tool, left)
            ((bar_at, bar),), ((marker_at, marker),) = axes.drawn("barh"), axes.drawn("plot")
            # The bar ends at its counts' place and the marker stands at its own, whole and over the bar.
            assert bar["left"] + bar_at[1] == pytest.approx(tip) and marker_at == ([pytest.approx(place)], [LANES[tool]]), (tool, left)
            assert marker["marker"] == "D" and marker["markerfacecolor"] == FILLS[tool] and marker["zorder"] > bar["zorder"], (tool, left)
            assert scale * math.log10(place / tip) == pytest.approx(left)
        # As far right of its bar's end as those stand left of it, a marker past the touch of the two has its line.
        for right in (body, body + 0.01, 6.5, 35.0, 200.0):
            tip = place * 10 ** (right / scale)
            row = {"label": "a row", "values": [placed(tool, "bar", tip), placed(tool, "marker", place)]}
            axes = Recorder()
            figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
            assert len(joins_of(axes)) == 1 and len(axes.drawn("plot")) == 3 and figure.join_end(tip, place, "adjusted", scale) is not None, (tool, right)


def test_a_join_ends_on_its_markers_body_on_every_pixel_row_it_inks():
    scale = 35.0
    half_line = figure.JOIN_WIDTH * figure.DPI / 72 / 2
    assert figure.MARKER_TAPERS == {">": 0.0, "D": 1.0} and set(figure.MARKER_TAPERS) == {marker for marker, _ in figure.MARKERS.values()}

    def drawn(tool: str, shape: str, place: float, stretch: float) -> tuple[list, list]:
        """The joins and the feet drawn for a marker whose body stands `stretch` pixels past its bar's end, each as
        its two ends in pixels right of the body's edge on the lane's line, with its style."""
        body = figure.marker_body(shape, 0)
        tip = place * 10 ** ((body + stretch) / scale)
        row = {"label": "a row", "values": [placed(tool, "bar", tip), placed(tool, "marker", 0 if shape == "exact" else place)]}
        axes = Recorder()
        figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
        into = [[([body - scale * math.log10(end / place) for end in args[0]], args[1], style) for args, style in found]
                for found in (joins_of(axes), feet_of(axes))]
        return into[0], into[1]

    for tool, stretch in itertools.product(("vep-rs", "fastvep"), (-0.3, 0.4, 3.0, 12.0, 60.0, 200.0)):
        # A diamond's body narrows towards its bar: one pixel row off the lane's line its edge stands a pixel
        # further from the bar's end. From half a pixel before the body's edge, inside the join's last dash, a solid
        # stretch in the join's colour, width and order is drawn on into the body, on the lane's line.
        ((join_at, lane, join),), ((foot_at, foot_lane, foot),) = drawn(tool, "adjusted", 0.001, stretch)
        assert foot_lane == lane == [LANES[tool]] * 2 and foot["linestyle"] == "solid" and foot["solid_capstyle"] == "butt", (tool, stretch)
        assert (foot["color"], foot["linewidth"], foot["zorder"]) == (join["color"], join["linewidth"], join["zorder"]) == (FILLS[tool], figure.JOIN_WIDTH, figure.JOIN_ORDER)
        assert foot_at == [pytest.approx(-figure.JOIN_LAP / 2, abs=1e-6), pytest.approx(figure.JOIN_LAP + half_line + 0.5, abs=1e-6)], (tool, stretch)
        # It starts a pixel and a half before the join's end. A join's dashes are set from its start, which is set
        # on a whole pixel, up to half a pixel from its place: the last dash, two pixels long or longer, still
        # holds the foot's start, which is set on a whole pixel too.
        assert join_at[1] - foot_at[0] == pytest.approx(1.5 * figure.JOIN_LAP, abs=1e-6)
        if join["linestyle"] != "solid":
            assert join["linestyle"][1][0] * figure.JOIN_WIDTH * figure.DPI / 72 >= 1.5 * figure.JOIN_LAP + 0.5, (tool, stretch)
        # A line's ink reaches half its width to either side of its pixel row, which lies up to half a pixel from
        # the lane's line: on the outermost row it inks the body begins that much further on. A marker and a
        # line's end are each set on a whole pixel, up to half a pixel from their places, and the foot still ends
        # on the body there.
        outermost = half_line + 0.5
        assert foot_at[1] - 0.5 - 0.5 >= figure.MARKER_TAPERS["D"] * outermost - 1e-9
        # It ends inside the body, short of the marker's place, where the body is wider than the line's reach.
        assert foot_at[1] + 0.5 + 0.5 < figure.marker_body("adjusted", 0) and figure.MARKER_TAPERS["D"] * (foot_at[1] - 1.0) >= outermost - 1e-9
        # The flat side of a triangle's body is as near the bar on every row: its join has no foot.
        (_,), feet = drawn(tool, "exact", 1e-7, stretch)
        assert feet == [], (tool, stretch)
    # A marker that touches its bar, and one left of its bar's end, have no join and so no foot.
    for stretch in (-0.6, -3.0, -20.0):
        assert drawn("fastvep", "adjusted", 0.001, stretch) == ([], [])
    assert figure.join_foot(1e-7, "exact", scale) is None
    start, stop = figure.join_foot(0.001, "adjusted", scale)
    assert 0.001 < stop < figure.join_drawn(0.01, 0.001, "adjusted", scale)[1] < figure.join_end(0.01, 0.001, "adjusted", scale) < start


def test_the_line_is_drawn_at_any_gap_and_is_a_solid_stub_where_the_stretch_is_shorter_than_one_dash():
    scale = 35.0
    width = figure.JOIN_WIDTH * figure.DPI / 72
    dash, gap = (length * width for length in figure.JOIN_DASHES)
    lap = figure.JOIN_LAP

    def line_of(stretch: float, shape: str, place: float):
        """The dashed line of a marker whose body stands `stretch` pixels past its bar's end, as it is asked of
        the axes: its two ends in pixels from the bar's end, towards the marker, and its style; None without one."""
        tip = place * 10 ** ((figure.marker_body(shape, 0) + stretch) / scale)
        row = {"label": "a row", "values": [placed("fastvep", "bar", tip), placed("fastvep", "marker", 0 if shape == "exact" else place)]}
        axes = Recorder()
        figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
        joins = joins_of(axes)
        ((_, bar),) = axes.drawn("barh")
        if not joins:
            return None
        ((ends, lane), style), = joins
        assert lane == [0.22, 0.22] and style["color"] == bar["color"] == FILLS["fastvep"]
        return [scale * math.log10(tip / end) for end in ends], style

    for shape, place in (("adjusted", 0.001), ("exact", 1e-7)):
        outline = figure.marker_reach(shape) - figure.marker_body(shape, 0)
        # No line where the body stands half a pixel or more over the bar's end.
        assert line_of(-figure.JOIN_TOUCH - 0.01, shape, place) is None and line_of(-2.0, shape, place) is None
        # The line is drawn at any gap, from a body that stands under half a pixel over the bar's end on, and under
        # one dash it is a solid stub with square ends: it is drawn on from the bar's end over the stretch and one
        # pixel into the body, no further.
        for stretch in (-figure.JOIN_TOUCH + 0.01, -0.01, 0.01, 0.3, 1.0, 2.0, dash - 0.01):
            (start, stop), style = line_of(stretch, shape, place)
            assert style["linestyle"] == "solid" and style["solid_capstyle"] == "butt", (shape, stretch)
            assert stop == pytest.approx(stretch + lap, abs=1e-6) and lap < figure.marker_body(shape, 0)
            # Where the marker's outline lies over the bar's end, or within a pixel of it, the stub is drawn from a
            # pixel before the outline, over the bar; from the bar's end otherwise.
            assert start == pytest.approx(min(0.0, stretch - outline - lap), abs=1e-6), (shape, stretch)
        # A stretch too short for a first dash that shows before a gap and a last dash of their least lengths is a
        # stub too, whatever its length.
        least = figure.JOIN_LEAST / figure.JOIN_DASHES[0] * (2 * dash + gap) - lap
        finest = least - (figure.JOIN_LEAST * width - figure.JOIN_SHOWN)
        assert dash < 4.0 < finest < 5.0 < least
        for stretch in (dash, 4.0, finest - 0.01):
            (start, stop), style = line_of(stretch, shape, place)
            assert style["linestyle"] == "solid" and start == pytest.approx(0.0, abs=1e-6) and stop == pytest.approx(stretch + lap, abs=1e-6)
        # From there to the stretch that holds two dashes of the least length the pattern is finer: a shorter
        # first dash before a gap and a last dash of their least lengths.
        for stretch in (finest + 0.01, 5.0, 5.5, least - 0.01):
            (start, stop), style = line_of(stretch, shape, place)
            offset, (on, off) = style["linestyle"]
            first = on - offset
            assert on == pytest.approx(figure.JOIN_LEAST) and 0 < first < on and style["dash_capstyle"] == "butt", (shape, stretch)
            assert start == pytest.approx(0.0, abs=1e-6) and stop == pytest.approx(stretch + lap, abs=1e-6)
            assert (first + off + on) * width == pytest.approx(stop, abs=1e-6) and on / off == pytest.approx(figure.JOIN_DASHES[0] / figure.JOIN_DASHES[1])
            # Two dashes show: the first a pixel and a half or more, the gap as in the usual pattern, and of the
            # last what stands before the body.
            assert first * width >= figure.JOIN_SHOWN - 1e-9 >= 1.5 - 1e-9 and off * width >= 1.5 and on * width - lap >= 1.8, (shape, stretch)
        # A longer one is dashed: it starts with a dash at the bar's end and ends with one where it is drawn to, a
        # pixel inside the body, so no ground stands before the body. Its dashes and gaps are fitted alike, no dash
        # under half again the line's width and no gap under a pixel and a half, and a dash's length past the lap
        # is seen before the body.
        fits = []
        for stretch in (least + 0.01, 8.0, 12.0, 20.0, 33.3, 100.0, 220.0):
            (start, stop), style = line_of(stretch, shape, place)
            offset, (on, off) = style["linestyle"]
            assert offset == 0 and style["dash_capstyle"] == "butt" and start == pytest.approx(0.0, abs=1e-6), (shape, stretch)
            assert stop == pytest.approx(stretch + lap, abs=1e-6)
            gaps = (stop / width - on) / (on + off)
            assert gaps == pytest.approx(round(gaps), abs=1e-6) and round(gaps) >= 1, (shape, stretch)
            assert on / off == pytest.approx(figure.JOIN_DASHES[0] / figure.JOIN_DASHES[1])
            assert on >= figure.JOIN_LEAST - 1e-9 >= 1.5 - 1e-9 and off * width >= 1.5 and on * width - lap >= 1.8
            # One dash and gap more would press the dashes under their least length.
            assert stop / width / ((round(gaps) + 2) * figure.JOIN_DASHES[0] + (round(gaps) + 1) * figure.JOIN_DASHES[1]) < figure.JOIN_LEAST / figure.JOIN_DASHES[0]
            fits.append(on / figure.JOIN_DASHES[0])
        # The shortest dashed line holds two dashes of the least length, and a long line's dashes are near the pair's own.
        assert fits[0] == pytest.approx(figure.JOIN_LEAST / figure.JOIN_DASHES[0], abs=2e-3) and max(fits) < 1.52 and all(0.93 < fit < 1.03 for fit in fits[-2:])
    # The style is one repeating pair, which a drawing library keeps whole at any length.
    assert figure.join_dashes(500.0)[1] == pytest.approx(tuple(figure.join_dashes(500.0)[1][:2])) and len(figure.join_dashes(500.0)[1]) == 2


def test_a_join_too_short_for_the_usual_dashes_shows_two_finer_ones_or_is_a_stub_and_the_texts_call_it_a_line(tmp_path):
    width = figure.JOIN_WIDTH * figure.DPI / 72
    least_dash, least_gap = figure.JOIN_LEAST * width, figure.JOIN_LEAST * figure.JOIN_DASHES[1] / figure.JOIN_DASHES[0] * width
    usual = 2 * least_dash + least_gap
    finest = figure.JOIN_SHOWN + least_gap + least_dash
    assert figure.JOIN_SHOWN == 1.5 and 5.8 < finest < 6.0 and 7.2 < usual < 7.3

    def dashes(length: float) -> list[tuple[float, float]]:
        """The stretches of ink of a join drawn `length` pixels on from its bar's end, in pixels from that end."""
        style = figure.join_dashes(length)
        if style == "solid":
            return [(0.0, length)]
        offset, (on, off) = style[0] * width, (style[1][0] * width, style[1][1] * width)
        inked, at = [], -offset
        while at < length - 1e-9:
            if at + on > 1e-9:
                inked.append((max(at, 0.0), min(at + on, length)))
            at += on + off
        return inked

    # The usual pattern from the length that holds two dashes of the least length: it starts with a dash and ends
    # with one, and so does the finer pattern of a shorter line, which holds exactly two.
    for length in (usual + 0.01, 9.0, 30.0, 150.0):
        inked = dashes(length)
        assert len(inked) >= 2 and inked[0][0] == 0.0 and inked[-1][1] == pytest.approx(length), length
    for length in (finest + 0.01, 6.0, 6.5, 7.0, usual - 0.01):
        (first, last) = dashes(length)
        assert first[0] == 0.0 and last[1] == pytest.approx(length) and last[1] - last[0] == pytest.approx(least_dash), length
        # The first dash is a pixel and a half or longer, the gap keeps its least length, over a pixel and a half,
        # and so does the last dash, of which what stands before the marker's body is seen.
        assert first[1] >= figure.JOIN_SHOWN - 1e-9 and last[0] - first[1] == pytest.approx(least_gap) and least_gap >= 1.5, length
        assert least_dash - figure.JOIN_LAP >= 1.8
    # The finer pattern runs into the usual one where the line holds two dashes of the least length.
    below, above = figure.join_dashes(usual - 1e-6), figure.join_dashes(usual + 1e-6)
    assert below[0] == pytest.approx(0.0, abs=1e-5) and below[1] == pytest.approx(above[1], abs=1e-5) and above[0] == 0
    # A shorter line cannot show two dashes at the image's size, its first dash under a pixel and a half: it is
    # a solid stub.
    for length in (0.5, 1.5, 3.0, 5.0, finest - 0.01):
        assert figure.join_dashes(length) == "solid" and dashes(length) == [(0.0, length)], length
    # Such a join is drawn: a diamond whose body stands under half a pixel past its bar's end is joined by a stub.
    scale = 35.0
    tip = 0.001 * 10 ** ((figure.marker_body("adjusted", 0) + 0.44) / scale)
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, {"label": "a row", "values": [placed("fastvep", "bar", tip), placed("fastvep", "marker", 0.001)]})], LANES, FILLS, 0.3, 1e-7, scale)
    ((_, stub),) = joins_of(axes)
    assert stub["linestyle"] == "solid"
    # So both texts call the join a line, for each F1 panel, and neither calls it dashed.
    document = values_of(tmp_path)
    for text, phrase in ((document["caption"], "joined to the end of its bar by a line where the two differ"),
                         (document["alt_text"], "as a marker joined to the end of the bar by a line where the two differ")):
        assert text.count(phrase) == 2 and "dashed" not in text
    assert figure.JOIN_CLAUSE == "joined to the end of its bar by a line where the two differ"


def test_the_dashed_line_joins_every_tools_bar_to_its_marker_and_the_caption_explains_it_while_the_legend_names_no_line(tmp_path):
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    # Every row of both F1 panels, drawn: each marker whose tool has a bar on the row has one dashed line, in its
    # tool's colour and on its tool's lane, from the bar's end to the marker. The comparator's as the released tool's.
    joined, touching = {"vep-rs": 0, "fastvep": 0}, {"vep-rs": 0, "fastvep": 0}
    for panel_id in ("a", "c"):
        for row in rows_of(document, panel_id).values():
            axes = Recorder()
            figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, 35.0)
            joins = joins_of(axes)
            bars = {value["tool"]: value for value in row["values"] if value["kind"] == "bar"}
            markers = [value for value in row["values"] if value["kind"] == "marker" and value["tool"] in bars]
            for marker in markers:
                tool = marker["tool"]
                end = figure.join_end(figure.discord(bars[tool], 1e-7), figure.discord(marker, 1e-7), figure.marker_shape(marker), 35.0)
                lines = [args for args, style in joins if style["color"] == FILLS[tool]]
                # A marker that touches its bar's end has no line; any other has one, from the bar's end to its body.
                tip = figure.discord(bars[tool], 1e-7)
                assert lines == ([] if end is None else [(drawn_join(tip, figure.discord(marker, 1e-7), figure.marker_shape(marker), 35.0), [LANES[tool]] * 2)]), row["label"]
                joined[tool] += end is not None
                touching[tool] += end is None
            assert len(joins) == sum(figure.join_end(figure.discord(bars[m["tool"]], 1e-7), figure.discord(m, 1e-7), figure.marker_shape(m), 35.0) is not None
                                     for m in markers), row["label"]
    assert joined["vep-rs"] + touching["vep-rs"] == 7 + 59 and joined["fastvep"] + touching["fastvep"] == 7 + 15
    assert joined["vep-rs"] > 0 and joined["fastvep"] > 0 and touching["fastvep"] > 0
    # The legend names no line: it has no entry for the dashed line, no entry's words speak of a line, and no
    # string of the image does. The marks of an F1 it names are the bar, the two markers and the cross.
    entries = figure.legend_entries(document)
    assert "join" not in {entry["mark"] for entry in entries} and not any("join" in kinds for kinds in figure.LEGEND_KINDS)
    assert [entry["mark"] for entry in figure.legend_rows(entries)[1]] == ["raw", "exact", "adjusted", "absent"]
    assert not any(word in text for text in figure.drawn_texts(document, registry) for word in ("ashed", " line"))
    assert not hasattr(figure, "JOIN_LEGEND")
    # The caption keeps its words on the line, for panel (a) and for panel (c).
    caption = document["caption"]
    a, c = caption[caption.index("(a) "):caption.index("(b) ")], caption[caption.index("(c) "):caption.index("(d) ")]
    assert "a marker is its tool's adjusted F1, joined to the end of its bar by a line where the two differ." in a
    assert "the marker its tool's adjusted F1, joined to the end of its bar by a line where the two differ, as in (a);" in c
    assert caption.count("by a line where the two differ") == 2 and "dashed" not in caption


def test_a_comparators_dashed_line_is_drawn_by_the_released_tools_rule_and_nothing_lengthens_it():
    scale = figure.f1_scale(261.0, 1e-7)
    # The same two figures on both tools' lanes, from far apart to equal, near 1 and near 0.
    pairs = (("0.900000000", "1.000000000"), ("0.990000000", "0.999900000"), ("0.998350000", "0.999000000"),
             ("0.999000000", "0.999100000"), ("0.250000000", "0.250200000"), ("0.500000000", "0.500000000"))
    lined = []
    for raw, adjusted in pairs:
        row = f1_row(("vep-rs", "bar", raw), ("vep-rs", "marker", adjusted), ("fastvep", "bar", raw), ("fastvep", "marker", adjusted))
        axes = Recorder()
        figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
        joins = joins_of(axes)
        markers = [(args, style) for args, style in axes.drawn("plot") if style.get("marker")]
        place = figure.discord(row["values"][2], 1e-7), figure.discord(row["values"][3], 1e-7)
        end = figure.join_end(*place, figure.marker_shape(row["values"][3]), scale)
        # The bar ends and the marker stands at the places of the row's counts, with or without a line.
        assert [style["left"] + args[1] for args, style in axes.drawn("barh")] == [pytest.approx(place[0])] * 2
        assert [args[0] for args, _ in markers] == [[pytest.approx(place[1])]] * 2
        lined.append(end is not None)
        if end is None:
            # The marker reaches its bar's end: the two touch, and neither tool has a line.
            assert joins == []
            continue
        # One rule for both tools: the same ends along the axis, each on its tool's lane, in its tool's colour, and
        # the same width, dashes and order.
        (ours, ours_style), (theirs, theirs_style) = joins
        assert ours[0] == theirs[0] == drawn_join(*place, figure.marker_shape(row["values"][3]), scale)
        assert (ours[1], theirs[1]) == ([-0.22, -0.22], [0.22, 0.22])
        assert (ours_style["color"], theirs_style["color"]) == (FILLS["vep-rs"], FILLS["fastvep"])
        assert {key: value for key, value in ours_style.items() if key != "color"} == {key: value for key, value in theirs_style.items() if key != "color"}
        # Nothing lengthens a line: it ends at the marker's body, and what is drawn of it past that edge lies
        # inside the body, short of the marker's centre; what is drawn of it before the bar's end lies over the bar.
        assert place[1] < ours[0][1] < end <= place[0] <= ours[0][0]
        assert scale * math.log10(ours[0][0] / place[1]) <= max(scale * math.log10(place[0] / place[1]), figure.marker_reach("adjusted") + figure.JOIN_LAP) + 1e-9
    # Far apart, and two pixels past the body's reach, the line is drawn; where the body reaches the bar's end,
    # or the marker stands on it, none is, for either tool.
    assert lined == [True, True, True, False, False, False]


def test_at_the_least_row_pitch_a_dashed_line_stands_clear_of_the_vertex_of_a_marker_on_the_next_lane():
    half_line = figure.JOIN_WIDTH * figure.DPI / 72 / 2
    # A marker's half height: the diamond's half diagonal and half the triangle's base.
    half_heights = {"adjusted": math.sqrt(0.5) * figure.MARKERS["adjusted"][1] * figure.DPI / 72, "exact": 0.5 * figure.MARKERS["exact"][1] * figure.DPI / 72}
    assert figure.MARKERS["adjusted"][0] == "D" and figure.MARKERS["exact"][0] == ">"
    tallest = max(half_heights.values())
    # Every pitch a figure is drawn at, in half pixels, and the pitches under the least down to 16 pixels.
    tightest = 16.0
    pitches = [tightest + half / 2 for half in range(round((figure.MAX_PITCH - tightest) * 2) + 1)]
    assert pitches[0] == tightest < figure.MIN_PITCH and figure.MIN_PITCH in pitches and pitches[-1] == figure.MAX_PITCH
    widened = []
    for pitch in pitches:
        bar, offsets = figure.lane_offsets(["vep-rs", "fastvep"], pitch)
        apart = (offsets["fastvep"] - offsets["vep-rs"]) * pitch
        # The edge of a dashed line on one lane stands at or past the vertex of the tallest marker on the other.
        assert apart - half_line >= tallest - 1e-9, pitch
        # The lanes lie about their row's line with their bars inside the row, a bar's height and the gap apart
        # wherever that keeps the line clear.
        assert bar == 0.3 and offsets["vep-rs"] == -offsets["fastvep"] and apart / pitch + bar <= 1.0, pitch
        if 0.44 * pitch - half_line >= tallest:
            assert apart == pytest.approx(0.44 * pitch), pitch
        else:
            widened.append(pitch)
            assert apart == pytest.approx(tallest + half_line), pitch
    # At 16 pixels alone the plain step would bring the line onto the diamond's vertex, by under a tenth of a
    # pixel, and the lanes stand further apart there. No pitch a figure is drawn at needs it: from the least
    # pitch on the lanes keep a bar's height and the gap.
    assert widened == [tightest] and 0 < tallest + half_line - 0.44 * tightest < 0.1
    # One lane stands on its row's line. More lanes than two leave a row no room to move them apart: they keep
    # a bar's height and the gap between them.
    assert figure.lane_offsets(["vep-rs"], figure.MIN_PITCH) == (0.3, {"vep-rs": 0.0})
    for count in range(3, figure.MAX_TOOLS):
        bar, offsets = figure.lane_offsets([str(lane) for lane in range(count)], figure.MIN_PITCH)
        assert bar == pytest.approx(0.72 / count) and offsets["1"] - offsets["0"] == pytest.approx(bar + 0.14), count
    # The panels' lanes are these.
    assert "bar, offsets = lane_offsets(lanes, pitch)" in SCRIPT.read_text(encoding="utf-8")


def test_rows_are_never_set_tighter_than_eighteen_pixels_and_a_figure_with_more_rows_is_refused_for_a_taller_canvas():
    assert (figure.MIN_PITCH, figure.MAX_PITCH) == (18.0, 24.0)
    # The pitch is the widest, in half pixels, at which every column's lines fit its room, and no wider than the most.
    assert figure.row_pitch([(1000.0, 30), (1000.0, 40)]) == 24.0
    assert figure.row_pitch([(1010.0, 49), (1010.0, 37)]) == 20.5 and figure.row_pitch([(1010.0, 37), (1010.0, 49)]) == 20.5
    assert figure.row_pitch([(900.0, 30), (1010.0, 49)]) == 20.5 and figure.row_pitch([(18.0 * 55, 55)]) == 18.0
    # More rows in the same room are set closer only down to the least pitch: no pitch under it is ever given.
    room = 1010.0
    fitting = [lines for lines in range(30, 80) if room / lines >= figure.MIN_PITCH]
    assert fitting[-1] == 56 and all(figure.MIN_PITCH <= figure.row_pitch([(room, lines)]) <= figure.MAX_PITCH for lines in fitting)
    assert all(figure.row_pitch([(room, lines)]) * lines <= room for lines in fitting)
    # One line more is refused by name, in either column, and nothing is drawn tighter: the refusal names the
    # taller canvas such a figure takes.
    for columns in ([(room, 57)], [(room, 40), (room, 57)], [(room, 57), (room, 40)]):
        with pytest.raises(SystemExit) as refused:
            figure.row_pitch(columns)
        assert str(refused.value) == (
            "ERROR: [render_release_figure] the figure's rows do not fit the canvas: 57 lines in a column at 17.5 pixels each, under the least "
            "pitch of 18; a figure with more rows takes a taller canvas, never a tighter pitch")
    # The drawing takes its pitch from this rule alone.
    assert SCRIPT.read_text(encoding="utf-8").count("pitch = row_pitch(") == 1
    # What the least pitch keeps, for the two lanes of a row: a dashed line on one lane stands clear of the tip of
    # the outline of a diamond on the other, and that tip stands clear of the next row's bar. Half a pixel
    # tighter the line would lie over the tip.
    half_line = figure.JOIN_WIDTH * figure.DPI / 72 / 2
    tip = figure.marker_reach("adjusted")

    def clear(pitch: float) -> tuple[float, float]:
        """The pixels between a line and the outline's tip on the other lane, and between that tip and the next row's bar."""
        bar, offsets = figure.lane_offsets(["vep-rs", "fastvep"], pitch)
        apart = (offsets["fastvep"] - offsets["vep-rs"]) * pitch
        return apart - half_line - tip, (pitch - apart) - bar * pitch / 2 - tip

    for half in range(round((figure.MAX_PITCH - figure.MIN_PITCH) * 2) + 1):
        assert min(clear(figure.MIN_PITCH + half / 2)) >= 0, half
    assert clear(figure.MIN_PITCH - 0.5)[0] < 0


def test_the_texts_say_that_a_marker_is_joined_to_its_bar_where_the_two_differ_and_at_its_bars_end_means_they_are_equal_as_drawn(tmp_path):
    sentence = "A marker standing at the end of its bar means that the tool's raw and adjusted F1 are equal as drawn."
    place = "On an F1 axis a value is a place, not a length:"
    joined = "by a line where the two differ"
    # The caption says it once, right after the sentence on places, and the alternative text once, right after its
    # panels' sentences, whatever panels the records give and in a draft too.
    for name, records, flags in (("four", record_set(), ()), ("two", record_set(plugins=None), ()), ("draft", record_set(plugins=ungated()), ("--draft",))):
        document = values_of(tmp_path / name, records, *flags)
        caption, alt = document["caption"], document["alt_text"]
        assert caption.count(sentence) == 1 and caption.count("equal as drawn") == 1 and "agree at this scale" not in caption, name
        assert caption.index(place) < caption.index(sentence) and f"and ends at the place of its F1. {sentence} " in caption, name
        assert alt.count(sentence) == 1 and alt.count("equal as drawn") == 1 and "agree at this scale" not in alt, name
        last = max(alt.index(f"({panel_id}) ") for panel_id in "abcd" if f"({panel_id}) " in alt)
        assert alt.index(sentence) > last and f". {sentence}" in alt, name
        # Both texts say of every F1 panel that a marker is joined to its bar by a line where the two differ,
        # and neither says that every marker is joined.
        panels = sum(f"({panel_id}) " in caption for panel_id in "ac")
        assert caption.count(f"joined to the end of its bar {joined}") == panels == caption.count("by a line") and "dashed" not in caption, name
        assert alt.count(f"as a marker joined to the end of the bar {joined}") == panels == alt.count("by a line") and "dashed" not in alt, name
    # The sentence names no tool and no mark but the marker and the bar, and it is no string of the image.
    document = values_of(tmp_path / "drawn")
    assert sentence not in figure.drawn_texts(document, figure.read_registry(REGISTRY)) and "line" not in sentence
    # What the texts say is what is drawn. A marker within its own half width of its bar's end stands on it without
    # a line: the two figures are equal as drawn.
    scale = figure.f1_scale(261.0, 1e-7)
    row = f1_row(("fastvep", "bar", "0.999000000"), ("fastvep", "marker", "0.999100000"))
    apart = abs(pixels_from_the_edge(figure.discord(row["values"][1], 1e-7), scale) - pixels_from_the_edge(figure.discord(row["values"][0], 1e-7), scale))
    assert 0 < apart < figure.marker_reach("adjusted")
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
    assert joins_of(axes) == [] and feet_of(axes) == []
    # One whose figure differs from its bar's by more is joined to the bar's end, and the fixture's rows hold both.
    row = f1_row(("fastvep", "bar", "0.990000000"), ("fastvep", "marker", "0.999100000"))
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, 1e-7, scale)
    assert len(joins_of(axes)) == 1
    lined = set()
    for panel_id in ("a", "c"):
        for row in rows_of(document, panel_id).values():
            bars = {value["tool"]: value for value in row["values"] if value["kind"] == "bar"}
            lined |= {figure.join_end(figure.discord(bars[m["tool"]], 1e-7), figure.discord(m, 1e-7), figure.marker_shape(m), scale) is not None
                      for m in row["values"] if m["kind"] == "marker" and m["tool"] in bars}
    assert lined == {True, False}


def test_every_row_carries_its_adjusted_marker_at_exactly_one_too(tmp_path):
    full = f1_row(("vep-rs", "bar", "1.000000000"), ("vep-rs", "marker", "1.000000000"), ("fastvep", "bar", "0.900000000"))
    short = f1_row(("vep-rs", "bar", "0.750000000"), ("vep-rs", "marker", "1.000000000"), ("fastvep", "bar", "1.000000000"))
    below = f1_row(("vep-rs", "bar", "0.990000000"), ("vep-rs", "marker", "0.999000000"))
    assert [figure.marker_shape(row["values"][1]) for row in (full, short, below)] == ["exact", "exact", "adjusted"]
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, full), (1, short), (2, below)], LANES, FILLS, 0.3, 1e-7, 35.0)
    assert len(axes.drawn("barh")) == 5
    # One marker a row, on its tool's lane: a triangle at the floor, at the end of a full bar as beside a shorter
    # one, and a diamond below 1. No row is left to read as unscored, and no value is written at a bar's end.
    markers = [(args, style) for args, style in axes.drawn("plot") if style.get("marker")]
    assert [(args[0], args[1], style["marker"]) for args, style in markers] == [
        ([1e-7], [-0.22], ">"), ([1e-7], [1 - 0.22], ">"), ([pytest.approx(0.001)], [2 - 0.22], "D")]
    assert axes.drawn("annotate") == []
    # Every adjusted F1 a record set holds is drawn, most of them at exactly 1, and the legend names the triangle.
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    held = [value for panel_id in ("a", "c") for group in panel(document, panel_id)["groups"] for row in group["rows"]
            for value in row["values"] if value["kind"] == "marker"]
    shapes = [figure.marker_shape(value) for value in held]
    assert len(held) > 57 and shapes.count("exact") > 57 and set(shapes) <= {"exact", "adjusted"}
    assert "1" not in figure.drawn_texts(document, registry)
    assert "exact" in {entry["mark"] for entry in figure.legend_entries(document)}


def test_a_row_at_the_maximum_carries_no_number_and_only_equal_counts_are_the_maximum(tmp_path):
    many = 10**10
    rounded = {"tool": "fastvep", "kind": "bar", "value": "1.000000000", "counts": [many - 1, many, many]}
    row = {"label": "a row", "values": [f1_value("vep-rs", "bar", "1.000000000"), f1_value("vep-rs", "marker", "1.000000000"), rounded]}
    floor = figure.discord_floor(row["values"])
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, floor, 35.0)
    # The row at the maximum: its bar runs to the floor with the triangle at its end, and no number is written.
    (full_at, full), (short_at, short) = axes.drawn("barh")
    ((marker_at, marker),) = [(args, style) for args, style in axes.drawn("plot") if style.get("marker")]
    assert full["left"] + full_at[1] == pytest.approx(floor) and (marker_at[0], marker["marker"]) == ([floor], ">")
    assert not axes.drawn("annotate") and not axes.drawn("text")
    # Only equal counts are the maximum. A bar whose unequal counts round to 1 ends short of the floor, at its
    # counts' place, and its value stays the figure nine places give.
    assert short["left"] + short_at[1] == pytest.approx(1e-10) and floor == 1e-11
    assert rounded["value"] == figure.f1(*rounded["counts"], "x") == "1.000000000"
    # No value of a record set's F1 panels is among the strings drawn, the 1 of a full bar least of all.
    document = values_of(tmp_path)
    drawn = figure.drawn_texts(document, figure.read_registry(REGISTRY))
    values = {value["value"] for panel_id in ("a", "c") for group in panel(document, panel_id)["groups"] for row in group["rows"]
              for value in row["values"] if value["kind"] != "text"}
    assert "1.000000000" in values and not values & set(drawn) and "1" not in drawn


def test_a_figure_that_prints_one_from_unequal_counts_has_no_label_stands_short_of_the_floor_and_wears_the_diamond(tmp_path):
    # A comparator's cell with one tuple of ten thousand million unmatched, raw and adjusted alike.
    many = 10**10
    scored = dataset_record()
    scored["concordance"]["s01"] = {"intersection": many - 1, "perl": many, "engine_tuples": many,
                                    "adjusted": {"intersection": many - 1, "perl": many, "engine_tuples": many}}
    document = values_of(tmp_path, record_set(fastvep=scored))
    row = rows_of(document, "a")[("Genome-wide ClinVar", "ClinVar GRCh37")]
    found = {(value["tool"], value["kind"]): value for value in row["values"]}
    bar, marker = found[("fastvep", "bar")], found[("fastvep", "marker")]
    # Its figure is in the values file as nine places give it, with the counts that say it is not 1.
    assert bar["value"] == marker["value"] == "1.000000000" and bar["counts"] == marker["counts"] == [many - 1, many, many]
    # Drawn, the bar ends short of the floor, at its counts' place, and its marker is the diamond of a figure below 1.
    floor = figure.discord_floor(row["values"])
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, row)], LANES, FILLS, 0.3, floor, 35.0)
    ((bar_at, drawn_bar),) = [(args, style) for args, style in axes.drawn("barh") if style["color"] == FILLS["fastvep"]]
    ((marker_at, drawn_marker),) = [(args, style) for args, style in axes.drawn("plot") if style.get("markerfacecolor") == FILLS["fastvep"]]
    assert drawn_bar["left"] + bar_at[1] == pytest.approx(1e-10) and marker_at[0] == [pytest.approx(1e-10)] and floor == 1e-11
    assert figure.marker_shape(marker) == "adjusted" and drawn_marker["marker"] == figure.MARKERS["adjusted"][0] == "D"
    # No label is drawn for it: no text beside the row, and its figure is no string of the image.
    assert not axes.drawn("annotate") and not axes.drawn("text")
    drawn = figure.drawn_texts(document, figure.read_registry(REGISTRY))
    assert "1.000000000" not in drawn and "1" not in drawn and not any("1.000000000" in text for text in drawn)


def test_a_cross_is_drawn_only_where_the_record_lists_the_plugin_as_absent_for_the_comparator(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"
    document = values_of(tmp_path)
    listed = {item["plugin"] for item in PLUGINS["absent"] if item["engine"] == "fastvep"}
    assert {"GWAS", "LoFTEE"} <= listed
    # Every cross is a plugin the record lists as absent for the comparator, and its key names that entry.
    crossed = set()
    for (_, label), row in rows_of(document, "c").items():
        plugin = label.removesuffix(" synthetic").rsplit(" ", 1)[0].split(" (")[0]
        for value in row["values"]:
            if value["kind"] == "text":
                assert plugin in listed and (value["tool"], value["value"]) == ("fastvep", "no such source"), label
                assert value["key"] == f"absent[plugin={plugin},engine=fastvep].reason" and value["record"] == name
                assert ("fastvep", "bar") not in marks(row) and ("fastvep", "marker") not in marks(row)
                crossed.add(plugin)
    assert crossed == listed & {row["plugin"].split(" (")[0] for row in PLUGINS["rows"]} and {"GWAS", "LoFTEE"} <= crossed
    # A record that does not list a plugin draws no cross for it: the row is refused, and in a draft it carries
    # no mark of the comparator, with a printed note.
    unlisted = copy.deepcopy(PLUGINS)
    unlisted["absent"] = [item for item in unlisted["absent"] if item["plugin"] not in ("GWAS", "LoFTEE")]
    assert refusal(tmp_path / "unlisted", record_set(plugins=unlisted)) == (
        f"{name} neither scores fastVEP per feature on GWAS GRCh38 nor lists the plugin as absent for it; every gate row shows each tool "
        "the record scores")
    done, values = run(tmp_path / "draft", record_set(plugins=unlisted), "--draft")
    assert done.returncode == 0, done.stderr
    rows = rows_of(json.loads(values.read_text(encoding="utf-8")), "c")
    for label in ("GWAS GRCh38", "LoFTEE GRCh38"):
        assert {value["tool"] for value in rows[(REAL_GROUP, label)]["values"]} == {"vep-rs"}
        assert f"note: {name} neither scores fastVEP per feature on {label} nor lists the plugin as absent for it" in done.stdout
    # A plugin the record lists as absent and also scores for the comparator draws the figure, and no cross.
    both = copy.deepcopy(PLUGINS)
    both["absent"].append({"plugin": "LOEUF", "engine": "fastvep", "reason": "no such source"})
    loeuf = marks(rows_of(values_of(tmp_path / "both", record_set(plugins=both)), "c")[(REAL_GROUP, "LOEUF GRCh37")])
    assert ("fastvep", "bar") in loeuf and ("fastvep", "text") not in loeuf


def test_the_reason_the_crosses_stand_for_is_one_sources_entry_that_names_where_it_is_read(tmp_path):
    name = f"2001-01-01-plugins-{RELEASE}.json"
    document = values_of(tmp_path)
    crossed = [value for row in rows_of(document, "c").values() for value in row["values"] if value["kind"] == "text"]
    assert len(crossed) == 43 and {value["value"] for value in crossed} == {"no such source"}
    # One entry for all the crosses of the comparator: its text is the reason, its key the record's list of
    # absent plugins by the comparator's engine.
    entry = {"text": "no such source", "record": name, "key": "absent[engine=fastvep].reason (one reason for every plugin the record lists as absent for the engine)"}
    assert [source for source in document["sources"] if source["key"].startswith("absent[")] == [entry]
    assert entry["key"] == figure.ABSENT_KEY.format(engine="fastvep")
    # The record bears it out: every plugin it lists as absent for that engine states this reason, and each
    # cross's own key names one of those entries.
    listed = [item for item in PLUGINS["absent"] if item["engine"] == "fastvep"]
    assert listed and {item["reason"] for item in listed} == {entry["text"]}
    assert {value["key"] for value in crossed} <= {f"absent[plugin={item['plugin']},engine=fastvep].reason" for item in listed}
    # Both texts give the words the entry holds, the caption in double quotes.
    assert document["caption"].count('a cross marks a plugin for which the comparator has no counterpart ("no such source")') == 1
    assert document["alt_text"].count("a cross (no such source) on 43") == 1
    # A panel without a cross quotes no reason and keeps no such entry.
    scored = copy.deepcopy(PLUGINS)
    scored["rows"] = [row for row in scored["rows"] if (row["plugin"], row["assembly"]) == ("LOEUF", "GRCh37")]
    scored["custom_annotation"] = []
    plain = values_of(tmp_path / "scored", record_set(plugins=scored))
    assert "a cross marks" not in plain["caption"] and not any(source["key"].startswith("absent[") for source in plain["sources"])


def test_no_row_of_an_f1_panel_carries_a_digit_of_its_figure_and_no_row_or_axis_label_a_version(tmp_path):
    registry = figure.read_registry(REGISTRY)
    document = values_of(tmp_path)
    drawn = figure.drawn_texts(document, registry)
    # No figure of a dataset row or of a plugin row is a string of the image: a row's figure is read from its marks.
    for panel_id in ("a", "c"):
        figures = {value["value"] for row in rows_of(document, panel_id).values() for value in row["values"] if value["kind"] != "text"}
        assert figures and not figures & set(drawn) and not any(figure_text in text for figure_text in figures for text in drawn), panel_id
    # No row of the panel of counts and no axis label names a version: the legend's tool entries alone do.
    versions = {version for tool in document["tools"] for version in tool["version"].split(" and ")}
    assert [row["label"] for group in panel(document, "d")["groups"] for row in group["rows"]] == ["vep-rs", "fastVEP"]
    labels = [line for p in document["panels"] for line in figure.axis_lines(document, registry, p)]
    assert not any(version in line for version in versions for line in labels)
    named = [entry["label"] for entry in figure.legend_entries(document) if entry["mark"] in figure.LEGEND_KINDS[0]]
    assert sorted(text for text in drawn if any(version in text for version in versions)) == sorted(named)
    # The triangle at the end of a full bar stays: a row at the maximum is its full bar and that triangle.
    full = f1_row(("vep-rs", "bar", "1.000000000"), ("vep-rs", "marker", "1.000000000"))
    axes = Recorder()
    figure.draw_f1_marks(axes, [(0, full)], LANES, FILLS, 0.3, 1e-7, 35.0)
    ((bar_at, bar),) = axes.drawn("barh")
    ((marker_at, marker),) = [(args, style) for args, style in axes.drawn("plot") if style.get("marker")]
    assert bar["left"] + bar_at[1] == pytest.approx(1e-7) and (marker_at[0], marker["marker"]) == ([1e-7], ">")


# From the widest F1 plot to the narrowest the canvas admits, at the seventh decade and at the tenth.
@pytest.mark.parametrize(("floor", "scale", "height"), [
    (1e-7, 35.0, 400.0), (1e-7, figure.MIN_PLOT / 7.35, 120.0), (1e-10, figure.MIN_PLOT / 10.35, 900.0)], ids=["wide", "narrow", "tenth"])
def test_the_f1_axis_carries_a_break_mark_before_the_tick_that_reads_one(floor, scale, height):
    def before(place: float) -> float:
        """The pixels a place on the axis lies before the tick that reads 1."""
        return scale * math.log10(place / floor)

    axes, along = Recorder(), object()
    figure.draw_axis_break(axes, along, floor, scale, height)
    (cut_at, cut), *strokes = axes.drawn("plot")
    # The axis line is cut, in the ground's colour and over the line, between the last decade's tick and the tick
    # that reads 1: 1 is no decade of 1 - F1, and the cut says the scale ends before it.
    assert cut["color"] == figure.GROUND and cut_at[1] == [0, 0] and cut["zorder"] > 2.5 and cut["linewidth"] > 0.6
    ends = [before(place) for place in cut_at[0]]
    assert 0 < ends[1] < ends[0] < scale and ends[0] - ends[1] == pytest.approx(figure.BREAK_WIDTH)
    # A leaning stroke in ink marks each end of the cut, crossing the axis line from under it to over it.
    assert len(strokes) == 2
    for end, ((along_axis, heights), stroke) in zip(ends, strokes):
        low, high = (before(place) for place in along_axis)
        assert stroke["color"] == figure.INK and stroke["zorder"] > cut["zorder"]
        assert (low + high) / 2 == pytest.approx(end) and low > high and 0 < high and low < scale
        assert heights == [pytest.approx(-figure.BREAK_REACH / height), pytest.approx(figure.BREAK_REACH / height)]
    # Every part is placed by the axis's values and the plot's height, and is drawn outside the plot too.
    assert all(style["transform"] is along and style["clip_on"] is False for _, style in axes.drawn("plot"))
    # The mark is large enough to read as a break at the README's width: the cut leaves three pixels or more of
    # the axis line open between the strokes, and a stroke rises four or more over the line and as far under it.
    stroke = strokes[0][1]["linewidth"] * figure.DPI / 72
    assert figure.BREAK_WIDTH - stroke >= 3.0 and figure.BREAK_REACH >= 4.0
    # A stroke leans as far along the axis as it rises over it, so under the axis line it is not taken for a tick.
    assert 2 * figure.BREAK_LEAN >= figure.BREAK_REACH


def test_wall_time_axis_draws_no_minor_ticks():
    axes = Recorder()
    ticks = figure.set_time_axis(axes, [1.79, 265212.0])
    assert ticks == [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]
    assert [name for name, _ in axes.calls] == ["set_xscale", "set_xlim", "set_xticks", "set_xticklabels", "minorticks_off"]
    assert axes.calls[0][1] == ("log",) and axes.calls[1][1] == figure.seconds_limits([1.79, 265212.0])
    assert axes.calls[2][1] == (ticks,) and axes.calls[3][1] == (["1", "10", "100", "1K", "10K", "100K"],)


# The committed figure and the README


def png(width: int, height: int, text: dict[str, str]) -> bytes:
    """A blank grayscale PNG of the given size with `text` as its text chunks."""
    def chunk(kind: bytes, body: bytes) -> bytes:
        return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body))

    out = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 0, 0, 0, 0))
    for keyword, value in text.items():
        out += chunk(b"tEXt", keyword.encode("latin-1") + b"\x00" + value.encode("latin-1"))
    return out + chunk(b"IDAT", zlib.compress(b"\x00" * ((width + 1) * height))) + chunk(b"IEND", b"")


def test_png_helpers_read_the_size_and_the_stamp():
    data = png(1080, 1350, {"values_sha256": "ab" * 32, "Title": "x"})
    assert figure.png_size(data) == (1080, 1350)
    assert figure.png_text(data) == {"values_sha256": "ab" * 32, "Title": "x"}
    with pytest.raises(SystemExit, match=r"ERROR: \[render_release_figure\] the image is not a PNG"):
        figure.png_size(b"GIF89a")


def tables(readme: str) -> dict[str, list[dict[str, str]]]:
    """The tables of the README's current-release section by the heading above them, each row by its header cell."""
    start = re.search(r"^## Current release: .*$", readme, re.MULTILINE)
    assert start, "README.md has no current-release section"
    section = readme[start.end():]
    following = re.search(r"^## ", section, re.MULTILINE)
    section = section[: following.start()] if following else section
    out: dict[str, list[dict[str, str]]] = {}
    heading, header = "", None
    for line in section.splitlines():
        if line.startswith("#"):
            heading, header = line.lstrip("# ").strip(), None
        elif line.startswith("|"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if header is None:
                header = cells
            elif not set(line) <= set("|-: "):
                out.setdefault(heading, []).append(dict(zip(header, cells)))
        else:
            header = None
    return out


def one(rows: list[dict[str, str]], **cells: str) -> dict[str, str]:
    found = [row for row in rows if all(row.get(column) == text for column, text in cells.items())]
    assert len(found) == 1, f"{len(found)} README rows match {cells}"
    return found[0]


def readme_mismatches(readme: str, document: dict) -> list[str]:
    """Every value of the values file that differs from the README table cell it corresponds to. The pooled row
    of panel (a) is the README's row of the same name, which states both assemblies."""
    found = tables(readme)
    labels = {tool["id"]: f"{tool['name']} {tool['version']}" for tool in document["tools"]}
    comparators = {tool["id"] for tool in figure.read_registry(REGISTRY) if tool["role"] == "comparator"}
    out = []

    def compare(where: str, value: str, shown: str) -> None:
        if shown.replace(",", "") != value:
            out.append(f"{where}: the figure states {value}, the README {shown}")

    for group in panel(document, "a")["groups"]:
        for row in group["rows"]:
            name, assembly = row["label"].rsplit(" ", 1)
            for value in row["values"]:
                column = {"bar": "Raw F1", "marker": "Adjusted F1"}[value["kind"]]
                if group["label"] == "Structural variants":
                    shown = one(found["Structural variants"], Engine=labels[value["tool"]], Assembly=assembly)
                elif group["label"].endswith(" whole-genome datasets, pooled"):
                    shown = one(found["Whole-genome concordance"], Engine=labels[value["tool"]], Dataset=row["label"])
                else:
                    shown = one(found["Whole-genome concordance"], Engine=labels[value["tool"]], Dataset=name, Assembly=assembly)
                compare(f"(a) {row['label']} {value['tool']} {value['kind']}", value["value"], shown[column])
    for group in panel(document, "b")["groups"]:
        name, assembly = group["label"].rsplit(" ", 1)
        for row in group["rows"]:
            shown = one(found["Whole-genome wall time"], Dataset=name, Assembly=assembly, Arch=row["label"])
            for value in row["values"]:
                column = next(c for c in shown if c.startswith(f"{labels[value['tool']]} (s)"))
                compare(f"(b) {group['label']} {row['label']} {value['tool']}", value["value"], shown[column].split(" ")[0])
    if any(p["id"] == "c" for p in document["panels"]):
        data = {REAL_GROUP: "ClinVar", "Committed table": "committed table",
                "Synthetic fixtures (data that cannot be redistributed)": "synthetic fixture", "Custom annotation (--custom sources)": "ClinVar"}
        rows = [row for heading in ("Plugins", "Synthetic fixtures (data that cannot be redistributed)", "Custom annotation (`--custom` sources)")
                for row in found.get(heading, [])]
        for group in panel(document, "c")["groups"]:
            for row in group["rows"]:
                plugin, assembly = row["label"].removesuffix(" synthetic").rsplit(" ", 1)
                for value in row["values"]:
                    engine = labels[value["tool"]]
                    if value["tool"] in comparators and value["kind"] != "text":
                        engine += " (per feature)"
                    matching = [r for r in rows if r["Engine"] == engine and r["Plugin or source"] == plugin and r["Assembly"] == assembly
                                and r["Data"].split(" (")[0] == data[group["label"]]]
                    assert len(matching) == 1, f"{len(matching)} README rows state {engine} on {row['label']}"
                    column = "Adjusted F1" if value["kind"] == "marker" else "Raw F1"
                    # The plugin tables state a row on its ClinVar dataset: a figure that sums entries of the
                    # population datasets with it has its row there and no cell that states it.
                    if f"{figure.POPULATION_LIST}[" in value["key"]:
                        continue
                    compare(f"(c) {row['label']} {value['tool']} {value['kind']}", value["value"], matching[0][column])
        support = found["Plugins each tool reads"]
        for row in (r for g in panel(document, "d")["groups"] for r in g["rows"]):
            # A row names its tool; the README's row names the tool with its version, as the legend does.
            shown = one(support, Engine=labels[row["values"][0]["tool"]])
            column = next(c for c in shown if c.startswith("Plugins read, of "))
            read = shown[column]
            for value in row["values"]:
                compare(f"(d) {row['label']} {value['kind']}", value["value"], read if value["kind"] == "bar" else f"{read} of {column.rsplit(' ', 1)[1]}")
    return out


def pair_problems(repo: Path) -> list[str]:
    """What is wrong with a tree's figure: without the README's image line neither file is in the tree; with it
    both are, the values file recomputes byte for byte from the committed records, the PNG is the canvas's
    size and carries the values file's digest, and every value equals its README cell."""
    readme = (repo / "README.md").read_text(encoding="utf-8")
    image, values = repo / figure.IMAGE, repo / figure.VALUES
    if not IMAGE_LINE.search(readme):
        return [f"{path.relative_to(repo)} is in the tree and README.md shows no image line for the figure" for path in (image, values) if path.exists()]
    missing = [f"README.md shows the figure and {path.relative_to(repo)} is not in the tree" for path in (image, values) if not path.exists()]
    if missing:
        return missing
    problems = []
    tools = figure.read_registry(REGISTRY)
    document, _ = figure.build_values(tools, figure.select(figure.discover(repo, None), tools, False), False)
    committed = values.read_bytes()
    if figure.values_bytes(document) != committed:
        problems.append("the values file does not recompute from the committed records")
    data = image.read_bytes()
    if figure.png_text(data).get(figure.STAMP_KEY) != hashlib.sha256(committed).hexdigest():
        problems.append("the PNG's stamp is not the values file's digest")
    if figure.png_size(data) != (document["canvas"]["width_px"], document["canvas"]["height_px"]):
        problems.append("the PNG is not the canvas's size")
    return problems + readme_mismatches(readme, json.loads(committed))


def test_the_tree_holds_the_figure_exactly_when_the_readme_shows_it():
    assert pair_problems(ROOT) == []


def committed_tree(tmp_path: Path) -> Path:
    """A copy of the repository's README and of the records its release names, with the figure's image line
    added when the README has none, and the pair written from those records."""
    repo = tmp_path / "repo"
    (repo / PROVENANCE).mkdir(parents=True)
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for record in figure.discover(ROOT, None):
        shutil.copy(record["path"], repo / PROVENANCE / record["path"].name)
    if not IMAGE_LINE.search(readme):
        readme = re.sub(r"^(## Current release: .*\n)", r"\1\n![The release's figure](docs/assets/comparison.png)\n", readme, count=1, flags=re.MULTILINE)
    (repo / "README.md").write_text(readme, encoding="utf-8")
    done = subprocess.run([sys.executable, str(SCRIPT), "--repo", str(repo), "--values-only"], capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr
    body = (repo / figure.VALUES).read_bytes()
    (repo / figure.IMAGE).write_bytes(png(1080, 1350, {"values_sha256": hashlib.sha256(body).hexdigest()}))
    return repo


def test_figure_of_the_committed_records_equals_the_readme_tables_cell_for_cell(tmp_path):
    repo = committed_tree(tmp_path)
    assert pair_problems(repo) == []
    document = json.loads((repo / figure.VALUES).read_text(encoding="utf-8"))
    assert all(record["file"].startswith(f"{PROVENANCE}/") for record in document["records"])
    assert sum(len(row["values"]) for p in document["panels"] for g in p["groups"] for row in g["rows"]) >= 50


def test_pair_check_names_a_stale_values_file_a_wrong_stamp_and_a_value_off_its_cell(tmp_path):
    repo = committed_tree(tmp_path)
    values, image = repo / figure.VALUES, repo / figure.IMAGE
    body = values.read_bytes()
    document = json.loads(body)
    first = document["panels"][0]["groups"][0]["rows"][0]
    first["values"][0]["value"] = "0.123456789"
    values.write_bytes(figure.values_bytes(document))
    problems = pair_problems(repo)
    assert "the values file does not recompute from the committed records" in problems
    assert "the PNG's stamp is not the values file's digest" in problems
    assert any(p.startswith(f"(a) {first['label']} vep-rs bar: the figure states 0.123456789, the README ") for p in problems)
    values.write_bytes(body)
    image.write_bytes(png(1080, 1080, {"values_sha256": hashlib.sha256(body).hexdigest()}))
    assert pair_problems(repo) == ["the PNG is not the canvas's size"]


def test_pair_check_refuses_a_figure_the_readme_does_not_show_and_a_shown_figure_that_is_missing(tmp_path):
    repo = committed_tree(tmp_path)
    readme = repo / "README.md"
    shown = readme.read_text(encoding="utf-8")
    (repo / figure.IMAGE).unlink()
    assert pair_problems(repo) == ["README.md shows the figure and docs/assets/comparison.png is not in the tree"]
    readme.write_text(IMAGE_LINE.sub("", shown), encoding="utf-8")
    assert pair_problems(repo) == ["docs/assets/comparison.values.json is in the tree and README.md shows no image line for the figure"]
    (repo / figure.VALUES).unlink()
    assert pair_problems(repo) == []


def plugin_tables() -> str:
    """The README's plugin tables for the fixture record, every figure from the record's counts in integer
    arithmetic: the rows of each engine per gate row, a row scored on real data from its ClinVar entries
    alone, and the support table."""
    header = "| Engine | Plugin or source | Data | Assembly | Raw F1 | Adjusted F1 |\n| --- | --- | --- | --- | --- | --- |\n"
    kinds = {"clinvar": ("real_data", "ClinVar"), "committed_table": ("committed_table", "committed table"),
             "synthetic_fixture": ("synthetic_fixtures", "synthetic fixture")}
    absent = {(a["plugin"], a["engine"]) for a in PLUGINS["absent"]}

    def lines(list_name: str, data: str, plugin: str, assembly: str) -> str:
        out = ""
        for engine, label in (("vep-rs", f"vep-rs {RELEASE}"), ("fastvep", f"fastVEP {COMPARATOR_VERSION}")):
            scored = False
            for projection in ("feature", "allele") if engine == "fastvep" else ("feature",):
                raw = plugin_sums(list_name, plugin, assembly, engine, projection, population=False)
                if raw is None:
                    continue
                scored = True
                adjusted = (plugin_sums(list_name, plugin, assembly, engine, projection, ("adj_intersection", "adj_perl", "adj_rust"), population=False)
                            if projection == "feature" else None)
                name = label if engine == "vep-rs" else f"{label} (per {projection})"
                out += f"| {name} | {plugin} | {data} | {assembly} | {nine_places(*raw)} | {nine_places(*adjusted) if adjusted else ''} |\n"
            if not scored and (plugin.split(" (")[0], engine) in absent and list_name != "custom_annotation":
                out += f"| {label} | {plugin} | {data} | {assembly} | no such source | |\n"
        return out

    body = {kind: "" for kind in kinds}
    for row in PLUGINS["rows"]:
        list_name, data = kinds[row["scored_on"]]
        body[row["scored_on"]] += lines(list_name, data, row["plugin"], row["assembly"])
    custom = "".join(lines("custom_annotation", "ClinVar", plugin, assembly)
                     for plugin, assembly in dict.fromkeys((e["plugin"], e["assembly"]) for e in PLUGINS["custom_annotation"]))
    support = {block["engine"]: block for block in PLUGINS["support"]}
    reads = f"| Engine | Plugins read, of {support['vep-rs']['of']} |\n| --- | --- |\n"
    reads += f"| vep-rs {RELEASE} | {support['vep-rs']['read']} |\n| fastVEP {COMPARATOR_VERSION} | {support['fastvep']['read']} |\n"
    return (f"### Plugins\n\n{header}{body['clinvar']}{body['committed_table']}\n"
            f"#### Synthetic fixtures (data that cannot be redistributed)\n\n{header}{body['synthetic_fixture']}\n"
            f"#### Custom annotation (`--custom` sources)\n\n{header}{custom}\n#### Plugins each tool reads\n\n{reads}\n")


SYNTHETIC_README = f"""# Project

## Current release: {RELEASE}

![The release's figure](docs/assets/comparison.png)

### Whole-genome concordance

| Engine | Dataset | Assembly | Raw F1 | Adjusted F1 |
| --- | --- | --- | --- | --- |
| vep-rs {RELEASE} | Population A | GRCh38 | 0.990000000 | 1.000000000 |
| vep-rs {RELEASE} | Population B | GRCh37 | 0.875000000 | 1.000000000 |
| vep-rs {RELEASE} | ClinVar | GRCh37 | 0.750000000 | 1.000000000 |
| vep-rs {RELEASE} | ClinVar | GRCh38 | 0.900000000 | 1.000000000 |
| vep-rs {RELEASE} | All four | GRCh37 and GRCh38 | 0.967213115 | 1.000000000 |
| fastVEP {COMPARATOR_VERSION} | Population A | GRCh38 | 0.500000000 | 0.500000000 |
| fastVEP {COMPARATOR_VERSION} | Population B | GRCh37 | 0.750000000 | 0.857142857 |
| fastVEP {COMPARATOR_VERSION} | ClinVar | GRCh37 | 0.250000000 | 0.250000000 |
| fastVEP {COMPARATOR_VERSION} | ClinVar | GRCh38 | 0.250000000 | 0.250000000 |
| fastVEP {COMPARATOR_VERSION} | All four | GRCh37 and GRCh38 | 0.500000000 | 0.504347826 |

### Whole-genome wall time

| Dataset | Assembly | Arch | vep-rs {RELEASE} (s) | fastVEP {COMPARATOR_VERSION} (s) | Ensembl VEP {REFERENCE_VERSION} (s) |
| --- | --- | --- | --- | --- | --- |
| Population A | GRCh38 | ARM | 100 | 1,000 (10.0x) | 10,000 (100x) |
| Population A | GRCh38 | x86 | 101 | 1,001 (9.9x) | 10,001 (99.0x) |
| Population B | GRCh37 | ARM | 10 | 100 (9.6x) | 1,000 (96.2x) |
| Population B | GRCh37 | x86 | 11 | 101 (9.6x) | 1,001 (95.3x) |
| ClinVar | GRCh37 | ARM | 1.50 | 20.00 (13.3x) | 700.00 (467x) |
| ClinVar | GRCh37 | x86 | 2.35 | 21.00 (9.0x) | 701.00 (299x) |
| ClinVar | GRCh38 | ARM | 3.00 | 20.00 (6.7x) | 700.00 (233x) |
| ClinVar | GRCh38 | x86 | 4.13 | 21.00 (5.1x) | 701.00 (170x) |

### Structural variants

| Engine | Assembly | Raw F1 | Adjusted F1 |
| --- | --- | --- | --- |
| vep-rs {RELEASE} | GRCh37 | 0.500000000 | 1.000000000 |
| vep-rs {RELEASE} | GRCh38 | 1.000000000 | 1.000000000 |
| fastVEP {COMPARATOR_VERSION} | GRCh37 | 0.250000000 | 0.250000000 |
| fastVEP {COMPARATOR_VERSION} | GRCh38 | 0.250000000 | 0.250000000 |

"""


def synthetic_tree(tmp_path: Path) -> Path:
    """A tree whose README states the synthetic records' figures in hand-written cells and the fixture record's
    in tables computed here, with the pair written from those records."""
    repo = tmp_path / "repo"
    write_records(repo / PROVENANCE, record_set())
    (repo / "README.md").write_text("placeholder\n\n## Current release: " + RELEASE + "\n", encoding="utf-8")
    done = subprocess.run([sys.executable, str(SCRIPT), "--repo", str(repo), "--values-only"], capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr
    body = (repo / figure.VALUES).read_bytes()
    (repo / "README.md").write_text(SYNTHETIC_README + plugin_tables() + "## Next\n", encoding="utf-8")
    (repo / figure.IMAGE).write_bytes(png(1080, 1350, {"values_sha256": hashlib.sha256(body).hexdigest()}))
    return repo


def test_every_value_of_all_four_panels_equals_its_readme_cell(tmp_path):
    repo = synthetic_tree(tmp_path)
    assert pair_problems(repo) == []
    readme = repo / "README.md"
    text = readme.read_text(encoding="utf-8")
    support = next(block for block in PLUGINS["support"] if block["engine"] == "fastvep")
    for old, new, problem in (
        (f"| vep-rs {RELEASE} | All four | GRCh37 and GRCh38 | 0.967213115 |", f"| vep-rs {RELEASE} | All four | GRCh37 and GRCh38 | 0.967213114 |",
         "(a) All four vep-rs bar: the figure states 0.967213115, the README 0.967213114"),
        ("| ClinVar | GRCh37 | x86 | 2.35 |", "| ClinVar | GRCh37 | x86 | 2.34 |", "(b) ClinVar GRCh37 x86 vep-rs: the figure states 2.35, the README 2.34"),
        (f"| fastVEP {COMPARATOR_VERSION} (per feature) | AlphaMissense | ClinVar | GRCh37 | 0.989026770 |",
         f"| fastVEP {COMPARATOR_VERSION} (per feature) | AlphaMissense | ClinVar | GRCh37 | 0.989026769 |",
         "(c) AlphaMissense GRCh37 fastvep bar: the figure states 0.989026770, the README 0.989026769"),
        (f"| fastVEP {COMPARATOR_VERSION} | AncestralAllele | ClinVar | GRCh37 | no such source |", f"| fastVEP {COMPARATOR_VERSION} | AncestralAllele | ClinVar | GRCh37 | 0.5 |",
         "(c) AncestralAllele GRCh37 fastvep text: the figure states no such source, the README 0.5"),
        (f"| fastVEP {COMPARATOR_VERSION} | {support['read']} |", f"| fastVEP {COMPARATOR_VERSION} | {support['read'] + 1} |",
         f"(d) fastVEP bar: the figure states {support['read']}, the README {support['read'] + 1}"),
    ):
        assert text.count(old) == 1, old
        readme.write_text(text.replace(old, new), encoding="utf-8")
        assert problem in pair_problems(repo)
    # The released tool's row sums the population dataset of its assembly with ClinVar, and the README's row states
    # ClinVar alone: the row is there, and no cell of it is held to the figure's value.
    alone = nine_places(*plugin_sums("real_data", "AlphaMissense", "GRCh37", "vep-rs", "feature", population=False))
    summed = nine_places(*plugin_sums("real_data", "AlphaMissense", "GRCh37", "vep-rs", "feature"))
    old = f"| vep-rs {RELEASE} | AlphaMissense | ClinVar | GRCh37 | {alone} |"
    assert alone != summed and text.count(old) == 1
    readme.write_text(text.replace(old, old.replace(alone, "0.500000000")), encoding="utf-8")
    assert pair_problems(repo) == []
    readme.write_text(text.replace(next(line for line in text.splitlines(keepends=True) if line.startswith(old)), ""), encoding="utf-8")
    with pytest.raises(AssertionError, match="0 README rows state vep-rs"):
        pair_problems(repo)
