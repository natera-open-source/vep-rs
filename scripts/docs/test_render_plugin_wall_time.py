"""Tests for `render_plugin_wall_time.py`: on trees the tests write, the section it gives a fabricated record
and each input it refuses; on this tree, the committed Benchmarks page against its record and the record
against the release record and the plugin support table."""

from __future__ import annotations

import copy
import json
import re
import sys
from decimal import Decimal
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
import render_plugin_wall_time as render

VERSION = "1.2.3"
RECORD_NAME = f"2001-01-01-plugin-walltime-v{VERSION}.json"
BEFORE = "# Benchmarks\n\nA paragraph the script leaves alone.\n\n"
AFTER = "\nA closing line the script leaves alone.\n"


def runs(ensembl, vep_rs):
    """A cell's runs, the same on both architectures but for one second more on x86_64."""
    return {
        "arm64": {"ensembl-vep": {"wall_seconds": ensembl}, "vep-rs": {"wall_seconds": vep_rs}},
        "x86_64": {"ensembl-vep": {"wall_seconds": [w + 1 for w in ensembl]}, "vep-rs": {"wall_seconds": [w + 1 for w in vep_rs]}},
    }


def record():
    """A record of the release: two rows without a plugin and two plugins, one of them under `--hgvs`."""
    return {
        "record": "plugin-walltime",
        "release": f"v{VERSION}",
        "arches": [
            {"arch": "arm64", "instance_type": "a1.large", "cpu": "AWS Graviton9", "host_vcpu_count": 8, "local_nvme": True},
            {"arch": "x86_64", "instance_type": "b1.large", "cpu": "Intel x86_64", "host_vcpu_count": 8, "local_nvme": True},
        ],
        "datasets": {"GRCh38": {"label": "ClinVar", "variants": 1234567}},
        "engines": {
            "ensembl-vep": {"name": "Ensembl VEP", "version": "999.1", "discarded_full_runs": 0,
                            "arguments": "vep --buffer_size 500 --fork 4 --stats_text --quiet <the row's arguments>"},
            "vep-rs": {"name": "vep-rs", "version": f"v{VERSION}", "discarded_full_runs": 1,
                       "arguments": "vep --fork 4 --buffer_size 500 --no_stats <the row's arguments>"},
        },
        "loftee_note": "A note the record carries.",
        "cells": [
            {"plugin": None, "hgvs": False, "assembly": "GRCh38", "runs": runs([1000.005], [1.0, 3.0, 2.004, 9.0, 0.5])},
            {"plugin": None, "hgvs": True, "assembly": "GRCh38", "runs": runs([2000.0], [4.0, 4.0, 5.0, 6.0, 7.0])},
            {"plugin": "First", "hgvs": False, "assembly": "GRCh38", "runs": runs([12345.675], [2.675])},
            {"plugin": "Second", "hgvs": True, "assembly": "GRCh38", "runs": runs([3000.0], [7.0])},
        ],
    }


def tree(tmp_path, records=None, page=None):
    """A repository of the release: the README's heading, the page with its markers, and the records given."""
    (tmp_path / "README.md").write_text(f"# vep-rs\n\n## Current release: v{VERSION}\n\nFigures.\n", encoding="utf-8")
    directory = tmp_path / render.RECORDS
    directory.mkdir(parents=True)
    for name, data in ({RECORD_NAME: record()} if records is None else records).items():
        (directory / name).write_text(json.dumps(data), encoding="utf-8")
    path = tmp_path / render.PAGE
    path.parent.mkdir(parents=True)
    path.write_text(BEFORE + f"{render.BEGIN}\n{render.END}\n" + AFTER if page is None else page, encoding="utf-8")
    return path


SECTION = f"""
## Wall time with each plugin

Seconds for one annotation run of the genome-wide ClinVar dataset of each assembly: 1,234,567 variants on GRCh38. \
A plugin's row is the run of the first row, `No plugin`, with that one plugin added. \
Second, which reads `--hgvs`, is added to the run under `--hgvs`. \
The 2 plugins are those the [Plugins](#plugins) tables score on ClinVar: a plugin whose data cannot be redistributed \
is not timed, and neither is TranscriptAnnotator, which reads a table of the user's own.

| Plugin | Assembly | Ensembl VEP 999.1, ARM (s) | vep-rs v{VERSION}, ARM (s) | Ensembl VEP 999.1, x86 (s) | vep-rs v{VERSION}, x86 (s) |
| --- | --- | --- | --- | --- | --- |
| No plugin | GRCh38 | 1,000.01 | 2.00 | 1,001.01 | 3.00 |
| No plugin, `--hgvs` | GRCh38 | 2,000.00 | 5.00 | 2,001.00 | 6.00 |
| First | GRCh38 | 12,345.68 | 2.68 | 12,346.68 | 3.68 |
| Second | GRCh38 | 3,000.00 | 7.00 | 3,001.00 | 8.00 |

**Method.** For a row and architecture both engines ran on one machine, `a1.large` (ARM Graviton9) or `b1.large` (x86 Intel), \
8 vCPU with local NVMe, each with `--fork 4` and `--buffer_size 500`; the timed span is the engine's process alone. \
Ensembl VEP's figure is one timed run. \
vep-rs's is the median of 5 timed runs on the rows without a plugin, and one timed run on every other row, after one discarded full run.

Ensembl VEP ran with `--stats_text`, as its concordance reference runs do, where the ClinVar rows of \
[Whole-genome wall time](#whole-genome-wall-time) time it with `--no_stats`. A note the record carries. \
Every timed run is in [`docs/concordance-provenance/{RECORD_NAME}`](../concordance-provenance/{RECORD_NAME}).

"""


def test_the_section_is_written_between_the_markers_and_the_rest_of_the_page_is_left(tmp_path, capsys):
    path = tree(tmp_path)
    render.main(["--repo", str(tmp_path)])
    assert capsys.readouterr().out == f"{render.PAGE}: plugin wall-time section of v{VERSION}\n"
    assert path.read_text(encoding="utf-8") == BEFORE + render.BEGIN + "\n" + SECTION + render.END + "\n" + AFTER


def test_a_second_run_changes_no_byte_and_the_check_passes(tmp_path):
    path = tree(tmp_path)
    render.main(["--repo", str(tmp_path)])
    written = path.read_bytes()
    render.main(["--repo", str(tmp_path)])
    assert path.read_bytes() == written
    render.main(["--repo", str(tmp_path), "--check"])


def test_the_check_writes_nothing_and_names_the_first_line_that_differs(tmp_path):
    path = tree(tmp_path)
    before = path.read_bytes()
    with pytest.raises(SystemExit) as excinfo:
        render.main(["--repo", str(tmp_path), "--check"])
    assert excinfo.value.code == (
        f"ERROR: [render_plugin_wall_time] {render.PAGE} line 6 is not the line the plugin wall-time record of v{VERSION} gives; "
        "scripts/docs/render_plugin_wall_time.py writes the section"
    )
    assert path.read_bytes() == before


def test_a_figure_is_the_median_of_its_runs_rounded_half_up_from_the_decimal_text(tmp_path):
    data = record()
    data["cells"] = [{"plugin": None, "hgvs": False, "assembly": "GRCh38", "runs": runs([0.005, 0.015], [1.005, 1.005, 7.0])}]
    for engine in data["engines"].values():
        engine["discarded_full_runs"] = 0
    tree(tmp_path, {RECORD_NAME: data})
    page, _ = render.build(tmp_path)
    assert "| No plugin | GRCh38 | 0.01 | 1.01 | 1.01 | 2.01 |\n" in page
    assert "Ensembl VEP's figure is the median of 2 timed runs. vep-rs's is the median of 3 timed runs.\n" in page


def test_a_tree_without_a_record_of_the_release_gets_an_empty_section(tmp_path):
    other = record()
    other["release"] = "v9.9.9"
    filled = BEFORE + render.BEGIN + "\n\nA section of another release.\n\n" + render.END + "\n" + AFTER
    path = tree(tmp_path, {RECORD_NAME: other, "2001-01-01-release.json": {"record": "release-addendum", "release": f"v{VERSION}"}}, filled)
    render.main(["--repo", str(tmp_path)])
    assert path.read_text(encoding="utf-8") == BEFORE + f"{render.BEGIN}\n{render.END}\n" + AFTER


def test_a_record_without_the_statistics_flag_or_a_note_states_neither(tmp_path):
    data = record()
    del data["loftee_note"]
    data["engines"]["ensembl-vep"]["arguments"] = "vep --buffer_size 500 --fork 4 --no_stats"
    tree(tmp_path, {RECORD_NAME: data})
    page, _ = render.build(tmp_path)
    assert f"\n\nEvery timed run is in [`docs/concordance-provenance/{RECORD_NAME}`]" in page
    assert "--stats_text" not in page


def without(key):
    def change(data):
        del data[key]

    return change


def set_runs(cell, engine, walls):
    def change(data):
        for arch in ("arm64", "x86_64"):
            data["cells"][cell]["runs"][arch][engine]["wall_seconds"] = walls

    return change


def other_fork(data):
    data["engines"]["vep-rs"]["arguments"] = "vep --fork 8 --buffer_size 500"


def no_run_on_an_architecture(data):
    del data["cells"][2]["runs"]["x86_64"]["vep-rs"]


def counts_differ_by_architecture(data):
    data["cells"][2]["runs"]["x86_64"]["vep-rs"]["wall_seconds"] = [1.0, 2.0]


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (without("cells"), f"{RECORD_NAME} has no `cells`"),
        (without("arches"), f"{RECORD_NAME} has no `arches`"),
        (no_run_on_an_architecture, f"{RECORD_NAME}: the cell of First on GRCh38, x86_64 has no `vep-rs`"),
        (set_runs(2, "vep-rs", []), f"{RECORD_NAME}: the cell of First on GRCh38, arm64, vep-rs: `wall_seconds` is not a list of numbers"),
        (set_runs(2, "vep-rs", ["2.5"]), f"{RECORD_NAME}: the cell of First on GRCh38, arm64, vep-rs: `wall_seconds` is not a list of numbers"),
        (other_fork, f"{RECORD_NAME}: the engines do not share one `--fork`"),
        (counts_differ_by_architecture, f"{RECORD_NAME}: vep-rs has different run counts on the architectures of First"),
        (set_runs(2, "vep-rs", [1.0, 2.0, 3.0]), f"{RECORD_NAME}: vep-rs's run counts [1, 3, 5] are not one count, or one count beside single runs"),
        (set_runs(1, "vep-rs", [4.0]), f"{RECORD_NAME}: vep-rs's rows without a plugin do not all hold 5 timed runs"),
    ],
    ids=["no-cells", "no-arches", "no-run-of-an-engine", "no-timed-run", "a-run-that-is-text", "another-fork", "counts-by-architecture",
         "three-run-counts", "a-bare-row-with-one-run"],
)
def test_a_record_the_section_cannot_be_written_from_is_refused(tmp_path, change, message):
    data = copy.deepcopy(record())
    change(data)
    tree(tmp_path, {RECORD_NAME: data})
    with pytest.raises(SystemExit) as excinfo:
        render.build(tmp_path)
    assert excinfo.value.code == f"ERROR: [render_plugin_wall_time] {message}"


def test_two_records_of_the_release_are_refused_by_name(tmp_path):
    tree(tmp_path, {RECORD_NAME: record(), "2001-01-02-plugin-walltime.json": record()})
    with pytest.raises(SystemExit) as excinfo:
        render.build(tmp_path)
    assert excinfo.value.code == (
        f"ERROR: [render_plugin_wall_time] 2 files of {render.RECORDS} are the plugin wall-time record of v{VERSION}: "
        f"{RECORD_NAME}, 2001-01-02-plugin-walltime.json"
    )


@pytest.mark.parametrize(
    "page",
    [BEFORE + AFTER, BEFORE + render.BEGIN + "\n" + AFTER, BEFORE + f"{render.END}\n{render.BEGIN}\n", BEFORE + f"{render.BEGIN}\n{render.BEGIN}\n{render.END}\n"],
    ids=["no-marker", "no-end", "end-before-begin", "two-begins"],
)
def test_a_page_without_its_two_markers_in_order_is_refused(tmp_path, page):
    tree(tmp_path, page=page)
    with pytest.raises(SystemExit) as excinfo:
        render.build(tmp_path)
    assert excinfo.value.code == (
        f"ERROR: [render_plugin_wall_time] {render.PAGE} does not hold one `{render.BEGIN}` line followed by one `{render.END}` line"
    )


def test_a_readme_without_one_current_release_heading_is_refused(tmp_path):
    tree(tmp_path)
    (tmp_path / "README.md").write_text("# vep-rs\n", encoding="utf-8")
    with pytest.raises(SystemExit) as excinfo:
        render.build(tmp_path)
    assert excinfo.value.code == "ERROR: [render_plugin_wall_time] README.md has 0 '## Current release: vX.Y.Z' headings; one is read"


# This tree


def committed_record():
    version = render.current_release(ROOT / "README.md")
    found = render.release_record(ROOT, version)
    assert found is not None, f"no plugin wall-time record of v{version}"
    return found


def test_the_committed_benchmarks_page_is_the_page_its_record_gives():
    page, _ = render.build(ROOT)
    assert (ROOT / render.PAGE).read_text(encoding="utf-8") == page


def test_the_committed_record_times_the_binaries_and_classes_the_release_record_measured():
    name, data = committed_record()
    releases = [json.loads(path.read_text(encoding="utf-8")) for path in sorted((ROOT / render.RECORDS).glob("*.json"))]
    (release,) = [r for r in releases if isinstance(r, dict) and r.get("record") == "release-addendum" and r.get("release") == data["release"]]
    for arch in data["arches"]:
        (stated,) = [a for a in release["arches"] if a["arch"] == arch["arch"]]
        assert {key: stated[key] for key in ("instance_type", "cpu", "host_vcpu_count", "local_nvme")} == {
            key: arch[key] for key in ("instance_type", "cpu", "host_vcpu_count", "local_nvme")
        }, name
        assert data["engines"]["vep-rs"]["binary_sha256"][arch["arch"]] == release["binary_proof"][arch["arch"]]["measured"]["sha256"], name
    assert data["engines"]["vep-rs"]["version"] == data["release"]


def support_table():
    """The plugin support table's rows as plugin -> (assemblies, status)."""
    rows = {}
    for line in (ROOT / "docs" / "plugins.md").read_text(encoding="utf-8").splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split(" | ")]
        if line.startswith("| ") and len(cells) == 6 and cells[4].startswith(("supported", "not ported")):
            name = re.sub(r"^\[([^\]]+)\]\([^)]*\)$", r"\1", cells[0])
            rows[name] = ([assembly.strip() for assembly in cells[3].split(",")], cells[4])
    return rows


def test_the_committed_record_times_the_plugins_verified_on_real_data_and_no_other():
    _, data = committed_record()
    table = support_table()
    real = {(plugin, assembly) for plugin, (assemblies, status) in table.items() if status == "supported, verified on real data"
            for assembly in assemblies}
    timed = [(cell["plugin"], cell["assembly"]) for cell in data["cells"] if cell["plugin"] is not None]
    assert len(real) > 30
    assert len(timed) == len(set(timed))
    assert set(timed) == {(plugin, assembly) for plugin, assembly in real if plugin != "TranscriptAnnotator"}
    assert not {plugin for plugin, _ in timed} & {plugin for plugin, (_, status) in table.items() if "synthetic" in status}


def test_the_committed_record_lists_its_rows_in_the_table_s_order_and_every_run_of_both_engines():
    name, data = committed_record()
    cells = data["cells"]
    keys = [(cell["plugin"] is not None, (cell["plugin"] or "").lower(), cell["hgvs"], cell["assembly"]) for cell in cells]
    assert keys == sorted(keys), name
    assert {(cell["hgvs"], cell["assembly"]) for cell in cells if cell["plugin"] is None} == {
        (hgvs, assembly) for hgvs in (False, True) for assembly in data["datasets"]
    }
    for cell in cells:
        assert cell["output_rows"] > data["datasets"][cell["assembly"]]["variants"], cell
        for arch in data["arches"]:
            for engine in render.ENGINES:
                walls = render.cell_runs(cell, arch["arch"], engine, name)
                assert all(wall > Decimal(0) for wall in walls), cell
                assert cell["hgvs"] == ("--hgvs" in cell["arguments"][engine].split()), cell


def test_the_release_wall_time_runs_time_ensembl_vep_without_its_statistics_report():
    """The section's sentence on `--stats_text` names the flag the release's own wall-time runs pass."""
    script = (ROOT / "scripts" / "concordance" / "run_clone_measurement.sh").read_text(encoding="utf-8")
    assert "--no_headers --no_stats --quiet --force_overwrite" in script
    assert render.STATISTICS_FLAG not in script
