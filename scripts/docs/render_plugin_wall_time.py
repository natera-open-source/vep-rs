#!/usr/bin/env python3
"""Write the plugin wall-time section of the Benchmarks page from the release's plugin wall-time record.

The page, `docs/benchmarks/current-release.md`, holds the section between the lines
`<!-- plugin-wall-time: begin -->` and `<!-- plugin-wall-time: end -->`. This script writes everything
between them and leaves the rest of the page as it is.

Inputs. The README's `## Current release: vX.Y.Z` heading names the release. The record is the file of
`docs/concordance-provenance/` that states `record: plugin-walltime` and that release or, where the tree
holds none of that release, the newest release before it that has one: the table's headers name the
versions measured. A tree with no such record gets an empty section.

The section. A heading, a paragraph that says what a row is, a table, a Method paragraph and a paragraph
of notes that ends by naming the record. The table has one row per cell of the record, in the record's
order: the plugin (`No plugin` for a cell without one, with `--hgvs` named where the cell passes it), the
assembly, and per architecture each engine's figure in seconds, the median of the cell's `wall_seconds`
of that engine and architecture, rounded half up to two decimals from the record's decimal text, with
thousands separators. Every count, version, machine class, option value and file name in the prose is
read from the record, and so are the number of timed runs behind an engine's figures and the full runs
it discarded before them. The notes say that Ensembl VEP wrote its statistics report where its command
in the record passes `--stats_text`, and carry the record's `loftee_note` where it has one.

Refused, each as one `ERROR: [render_plugin_wall_time] ...` line: no heading or several, two records of
one release, a record whose `release` is not `vX.Y.Z`, a page without the two marker lines in order, a record without a field the section reads, a
cell without a run of an engine on an architecture, run counts the Method paragraph has no words for,
and engines that differ in `--fork` or `--buffer_size`. The same tree gives the same bytes. Stdlib only.

usage: render_plugin_wall_time.py [--repo DIR] [--check]

--check writes nothing and exits 1, naming the first line that differs, when the page is not what the
script writes.
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMPONENT = "render_plugin_wall_time"
RECORDS = "docs/concordance-provenance"
PAGE = "docs/benchmarks/current-release.md"
KIND = "plugin-walltime"
BEGIN, END = "<!-- plugin-wall-time: begin -->", "<!-- plugin-wall-time: end -->"
HEADING = re.compile(r"^## Current release: v(\d+\.\d+\.\d+)\s*$", re.MULTILINE)
RELEASE = re.compile(r"v\d+\.\d+\.\d+")
# The engines in column order, and the architectures as the page's other tables name them.
ENGINES = ("ensembl-vep", "vep-rs")
ARCH_LABELS = {"arm64": "ARM", "x86_64": "x86"}
# Ensembl VEP's flag for its statistics report, which the release's own wall-time runs switch off.
STATISTICS_FLAG = "--stats_text"


def fail(message: str) -> None:
    sys.exit(f"ERROR: [{COMPONENT}] {message}")


def field(block, key: str, where: str):
    if not isinstance(block, dict) or key not in block:
        fail(f"{where} has no `{key}`")
    return block[key]


def current_release(readme: Path) -> str:
    """The version the README's current-release heading names."""
    found = HEADING.findall(readme.read_text(encoding="utf-8")) if readme.is_file() else []
    if len(found) != 1:
        fail(f"README.md has {len(found)} '## Current release: vX.Y.Z' headings; one is read")
    return found[0]


def release_number(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def release_record(repo: Path, version: str) -> tuple[str, dict] | None:
    """The plugin wall-time record the page shows as (file name, data): the release's own, else that of the
    newest release before it, or None when the tree holds neither."""
    found: dict[tuple[int, ...], list[tuple[str, dict]]] = {}
    for path in sorted((repo / RECORDS).glob("*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"), parse_float=Decimal)
        except (OSError, ValueError):
            fail(f"{RECORDS}/{path.name} is not readable JSON")
        if not isinstance(data, dict) or data.get("record") != KIND:
            continue
        stated = data.get("release")
        if not isinstance(stated, str) or not RELEASE.fullmatch(stated):
            fail(f"{path.name}: `release` is not vX.Y.Z")
        found.setdefault(release_number(stated[1:]), []).append((path.name, data))
    for number, records in found.items():
        if len(records) > 1:
            fail(f"{len(records)} files of {RECORDS} are the plugin wall-time record of v{'.'.join(map(str, number))}: "
                 + ", ".join(name for name, _ in records))
    shown = [number for number in found if number <= release_number(version)]
    return found[max(shown)][0] if shown else None


def seconds(value: Decimal) -> str:
    """A figure in seconds: two decimals, rounded half up, with thousands separators."""
    return f"{value.quantize(Decimal('0.01'), rounding=ROUND_HALF_UP):,}"


def count(n: int) -> str:
    return f"{n:,}"


def and_list(items: list[str]) -> str:
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def row_label(cell: dict) -> str:
    plugin = field(cell, "plugin", "a cell")
    if plugin is not None:
        return plugin
    return "No plugin, `--hgvs`" if field(cell, "hgvs", "a cell") else "No plugin"


def cell_runs(cell: dict, arch: str, engine: str, name: str) -> list[Decimal]:
    where = f"{name}: the cell of {row_label(cell)} on {field(cell, 'assembly', 'a cell')}"
    runs = field(field(field(cell, "runs", where), arch, f"{where}, `runs`"), engine, f"{where}, {arch}")
    walls = field(runs, "wall_seconds", f"{where}, {arch}, {engine}")
    if not isinstance(walls, list) or not walls or not all(isinstance(w, (Decimal, int)) and not isinstance(w, bool) for w in walls):
        fail(f"{where}, {arch}, {engine}: `wall_seconds` is not a list of numbers")
    return [Decimal(w) for w in walls]


def shared_option(engines: dict, option: str, name: str) -> str:
    """The value both engines' commands give `option`."""
    values = set()
    for engine in ENGINES:
        m = re.search(rf"{re.escape(option)} (\S+)", field(field(engines, engine, f"{name}, `engines`"), "arguments", f"{name}, {engine}"))
        values.add(m[1] if m else None)
    if len(values) != 1 or None in values:
        fail(f"{name}: the engines do not share one `{option}`")
    return values.pop()


def run_words(cells: list[dict], arches: list[str], engine: str, name: str) -> str:
    """How many timed runs stand behind the engine's figures, in the Method paragraph's words."""
    counts = {}
    for cell in cells:
        per_arch = {len(cell_runs(cell, arch, engine, name)) for arch in arches}
        if len(per_arch) != 1:
            fail(f"{name}: {engine} has different run counts on the architectures of {row_label(cell)}")
        counts.setdefault(per_arch.pop(), []).append(cell)

    def of(n: int) -> str:
        return "one timed run" if n == 1 else f"the median of {n} timed runs"

    if len(counts) == 1:
        return of(next(iter(counts)))
    if len(counts) != 2 or 1 not in counts:
        fail(f"{name}: {engine}'s run counts {sorted(counts)} are not one count, or one count beside single runs")
    many = max(counts)
    bare = [cell for cell in cells if cell["plugin"] is None]
    if any(cell not in counts[many] for cell in bare):
        fail(f"{name}: {engine}'s rows without a plugin do not all hold {many} timed runs")
    plugins = list(dict.fromkeys(cell["plugin"] for cell in counts[many] if cell["plugin"] is not None))
    if any(cell["plugin"] in plugins for cell in counts[1]):
        fail(f"{name}: {engine} holds {many} timed runs on one assembly of a plugin and one on another")
    where = and_list(["on the rows without a plugin"] + [f"on {plugin}" for plugin in plugins])
    return f"{of(many)} {where}, and one timed run on every other row"


def discarded_words(engine_block: dict, engine: str, name: str) -> str:
    """The clause for the full runs an engine made and discarded before its timed runs, empty for none."""
    discarded = field(engine_block, "discarded_full_runs", f"{name}, {engine}")
    if discarded == 0:
        return ""
    return ", after one discarded full run" if discarded == 1 else f", after {discarded} discarded full runs"


def class_label(arch: dict, name: str) -> str:
    """A machine class as the page's other tables write it: `c9gd.8xlarge` (ARM Graviton5)."""
    label = ARCH_LABELS.get(field(arch, "arch", f"{name}, `arches`"))
    if label is None:
        fail(f"{name}: no label for the architecture {arch['arch']}")
    cpu = field(arch, "cpu", f"{name}, `arches`").replace("AWS ", "").replace(" x86_64", "")
    return f"`{field(arch, 'instance_type', name)}` ({label} {cpu})"


def section(name: str, record: dict) -> list[str]:
    """The lines between the markers for the record."""
    engines = field(record, "engines", name)
    titles = {engine: f"{field(engines[engine], 'name', name)} {field(engines[engine], 'version', name)}" for engine in ENGINES if engine in engines}
    if len(titles) != len(ENGINES):
        fail(f"{name}: `engines` does not hold {and_list(list(ENGINES))}")
    arches = field(record, "arches", name)
    arch_ids = [field(arch, "arch", f"{name}, `arches`") for arch in arches]
    cells = field(record, "cells", name)
    datasets = field(record, "datasets", name)
    assemblies = list(dict.fromkeys(field(cell, "assembly", "a cell") for cell in cells))
    labels = {field(field(datasets, a, f"{name}, `datasets`"), "label", name) for a in assemblies}
    if len(labels) != 1:
        fail(f"{name}: the assemblies' datasets differ in `label`")
    plugins = list(dict.fromkeys(cell["plugin"] for cell in cells if field(cell, "plugin", "a cell") is not None))
    hgvs_plugins = list(dict.fromkeys(cell["plugin"] for cell in cells if cell["plugin"] is not None and field(cell, "hgvs", "a cell")))
    variants = and_list([f"{count(field(datasets[a], 'variants', name))} variants on {a}" for a in assemblies])
    vcpus = {field(arch, "host_vcpu_count", name) for arch in arches}
    if len(vcpus) != 1 or not all(field(arch, "local_nvme", name) for arch in arches):
        fail(f"{name}: the machine classes differ in vCPUs or lack local NVMe")

    dataset = labels.pop()
    intro = [
        f"Seconds for one annotation run of the genome-wide {dataset} dataset of each assembly: {variants}.",
        "A plugin's row is the run of the first row, `No plugin`, with that one plugin added.",
    ]
    if hgvs_plugins:
        verb = "reads" if len(hgvs_plugins) == 1 else "read"
        intro.append(f"{and_list(hgvs_plugins)}, which {verb} `--hgvs`, is added to the run under `--hgvs`.")
    intro.append(
        f"The {len(plugins)} plugins timed are verified on real data, and the [Plugins](#plugins) tables score them on "
        f"{dataset}: a plugin whose data cannot be redistributed is not timed, and neither is TranscriptAnnotator, which "
        "reads a table of the user's own."
    )

    classes = " or ".join(class_label(arch, name) for arch in arches)
    header = ["Plugin", "Assembly"] + [f"{titles[engine]}, {ARCH_LABELS[arch]} (s)" for arch in arch_ids for engine in ENGINES]
    table = ["| " + " | ".join(header) + " |", "| " + " | ".join("---" for _ in header) + " |"]
    for cell in cells:
        figures = [seconds(statistics.median(cell_runs(cell, arch, engine, name))) for arch in arch_ids for engine in ENGINES]
        table.append("| " + " | ".join([row_label(cell), cell["assembly"]] + figures) + " |")

    fork, buffer_size = shared_option(engines, "--fork", name), shared_option(engines, "--buffer_size", name)
    method = [
        (
            f"**Method.** For a row and architecture both engines ran on one machine, {classes}, {vcpus.pop()} vCPU with "
            f"local NVMe, each with `--fork {fork}` and `--buffer_size {buffer_size}`; the timed span is the engine's process alone."
        )
    ]
    for n, engine in enumerate(ENGINES):
        subject = f"{engines[engine]['name']}'s figure is" if n == 0 else f"{engines[engine]['name']}'s is"
        method.append(f"{subject} {run_words(cells, arch_ids, engine, name)}{discarded_words(engines[engine], engine, name)}.")
    notes = []
    if STATISTICS_FLAG in field(engines["ensembl-vep"], "arguments", name).split():
        notes.append(
            f"{engines['ensembl-vep']['name']} ran with `{STATISTICS_FLAG}`, as its concordance reference runs do, where the "
            "ClinVar rows of [Whole-genome wall time](#whole-genome-wall-time) time it with `--no_stats`."
        )
    if "loftee_note" in record:
        notes.append(record["loftee_note"])
    notes.append(f"Every timed run is in [`{RECORDS}/{name}`](../concordance-provenance/{name}).")
    return ["## Wall time with each plugin", "", " ".join(intro), "", *table, "", " ".join(method), "", " ".join(notes)]


def build(repo: Path) -> tuple[str, str]:
    """The page as the script writes it, and the release its record states, the README's where it has none."""
    version = current_release(repo / "README.md")
    path = repo / PAGE
    if not path.is_file():
        fail(f"no {PAGE}")
    lines = path.read_text(encoding="utf-8").split("\n")
    marks = [i for i, line in enumerate(lines) if line in (BEGIN, END)]
    if len(marks) != 2 or lines[marks[0]] != BEGIN or lines[marks[1]] != END:
        fail(f"{PAGE} does not hold one `{BEGIN}` line followed by one `{END}` line")
    found = release_record(repo, version)
    body = section(*found) if found else []
    inner = [""] + body + [""] if body else []
    return "\n".join(lines[: marks[0] + 1] + inner + lines[marks[1] :]), found[1]["release"][1:] if found else version


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Write the plugin wall-time section of the Benchmarks page from its record.")
    ap.add_argument("--repo", type=Path, default=HERE.parents[1], help="repository root")
    ap.add_argument("--check", action="store_true", help="write nothing; exit 1 when the page differs")
    args = ap.parse_args(argv)
    page, version = build(args.repo)
    path = args.repo / PAGE
    if args.check:
        have, want = path.read_text(encoding="utf-8").split("\n"), page.split("\n")
        for n, (a, b) in enumerate(zip(have + [None] * len(want), want + [None] * len(have)), 1):
            if a != b:
                fail(f"{PAGE} line {n} is not the line the plugin wall-time record of v{version} gives; "
                     "scripts/docs/render_plugin_wall_time.py writes the section")
        return
    path.write_text(page, encoding="utf-8")
    print(f"{PAGE}: plugin wall-time section of v{version}")


if __name__ == "__main__":
    main()
