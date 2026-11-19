#!/usr/bin/env python3
"""Write the documentation home page, `docs/index.md`, from its template and the release's records.

Template. `home_page.template.md` beside this script holds the page's text. It opens with one line per
panel of the release figure, `<<panel ID>>CLAUSE`: the words the figure's alternative text gives that
panel. The line `<<page>>` ends them and the page follows. In the page `<<NAME>>` is a place this script
fills, and `<<if NAME>>` to `<<end NAME>>` is a part written only when it holds; an `if` or `end` mark
alone on its line leaves no line. The site's own tokens, in double braces, pass through as they are.

Inputs. The README's `## Current release: vX.Y.Z` heading names the release. Its records are the files of
`docs/concordance-provenance/` that state it under `release`: one population record
(`record: population-measurements`), one release record (`release-addendum`) and, where the tree holds
one, the plugin record (`record_type: plugins`). The release figure's pair, `docs/assets/comparison.png`
and `docs/assets/comparison.values.json`, is read where the tree holds it. `comparators.tsv` beside this
script gives each tool's role and the names the records file it under.

Places, each from one record:
- `speed-up`: over the population record's `datasets` and the architectures arm64 and x86_64, the largest
  quotient of Ensembl VEP's `sum_of_cell_medians_sec` by vep-rs's (the `per_dataset` rows under
  `wall_time`), from the record's decimal text, rounded half up to a whole number.
- `gnomad-label`, `gnomad-variants`, `gnomad-arm-seconds`: of the dataset `gnomad_v4.1_genomes`, its
  `label` without `, whole genome`, its `variants`, and vep-rs's arm64 sum rounded half up to whole
  seconds.
- `reference-release`: the one Ensembl VEP release the image tags under the release record's
  `perl_reference` name.
- `whole-genome-datasets`: as a word, the population record's datasets plus the genome-wide ClinVar cells
  of the release record's `concordance`.
- `plugins-read`, inside the part `plugins`: the plugins the released tool reads in the plugin record's
  `support`. The part holds when the tree has a plugin record of the release.
- inside the part `figure`, from the values file: `figure-notes` (the lines drawn inside the image: its
  `notes`, then the `note` of each panel that carries one, each followed by a space), `released`,
  `reference` and `comparators` (each tool's name and version, by the registry's roles), `compared` (the
  released tool and the comparators), `panel-count` and `panels` (the clause of each panel the file holds,
  in its order). The parts `one-comparator` and `several-comparators` hold by the number of comparators
  among the file's tools. The part `figure` holds when the tree has the pair.

Every count and every figure in seconds carries thousands separators.

Refused, each as one `ERROR: [render_home_page] ...` line: no heading or several, a record of a kind
missing or doubled, a record without a field a place reads, one file of the figure's pair without the
other, a values file that names another image or whose released tool is not at the README's release, a
panel the template has no clause for, a mark the script does not know, and a place outside the part that
gives it a value. The same tree gives the same bytes. Stdlib only.

usage: render_home_page.py [--repo DIR] [--check | --values]

--check writes nothing and exits 1, naming the first line that differs, when `docs/index.md` is not what
the script writes. --values writes nothing and prints each place's value with the file and key it is
read from, tab-separated.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from decimal import ROUND_HALF_UP, Decimal, InvalidOperation, localcontext
from pathlib import Path

COMPONENT = "render_home_page"
HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "home_page.template.md"
REGISTRY = HERE / "comparators.tsv"
PAGE = "docs/index.md"
RECORDS = "docs/concordance-provenance"
IMAGE = "docs/assets/comparison.png"
VALUES = "docs/assets/comparison.values.json"
HEADING = re.compile(r"^## Current release: v(\d+\.\d+\.\d+)\s*$", re.MULTILINE)

KINDS = {"population-measurements": "population", "release-addendum": "release"}
KIND_NAMES = {
    "population": "population record (record: population-measurements)",
    "release": "release record (record: release-addendum)",
    "plugins": "plugin record (record_type: plugins)",
}
ARCHES = ("arm64", "x86_64")
ARM = "arm64"
ENGINE = "vep_rs"
REFERENCE_ENGINE = "ensembl_vep"
SUM = "sum_of_cell_medians_sec"
GNOMAD = "gnomad_v4.1_genomes"
WHOLE_GENOME_SUFFIX = ", whole genome"
CLINVAR_CELLS = ("s01", "s02")
ROLES = ("released", "comparator", "reference")
NUMBER_WORDS = (
    "zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve",
)

PANEL_LINE = re.compile(r"^<<panel ([^<>\s]+)>>(.+)$")
PAGE_MARK = "<<page>>"
BLOCK_PART = re.compile(r"^<<if ([a-z-]+)>>\n(.*?)^<<end \1>>\n", re.MULTILINE | re.DOTALL)
INLINE_PART = re.compile(r"<<if ([a-z-]+)>>(.*?)<<end \1>>")
MARK = re.compile(r"<<([^<>\n]*)>>")
# The characters Markdown reads as markup inside a sentence or an image's alternative text.
MARKUP = re.compile(r"([\\`*_\[\]<>])")
CONTROL = re.compile(r"[\x00-\x1f]")


def fail(message: str) -> None:
    sys.exit(f"ERROR: [{COMPONENT}] {message}")


def field(block, key: str, where: str):
    """The value under `key` of a record block, refused by name when it is absent or the block is no
    object."""
    value = block.get(key) if isinstance(block, dict) else None
    if value is None:
        fail(f"{where} carries no {key}")
    return value


def integer(value, key: str, where: str) -> int:
    """A record's integer, refused by name when it is not one: `bool` is a subclass of `int`."""
    if isinstance(value, bool) or not isinstance(value, int):
        fail(f"{where} carries {key} that is not an integer: {value!r}")
    return value


def seconds(value, where: str) -> Decimal:
    """A record's wall time as a Decimal of its own text, refused by name when it is not a finite positive
    number."""
    try:
        x = Decimal(str(value))
    except InvalidOperation:
        fail(f"{where} carries {SUM} that is not a number: {value!r}")
    if isinstance(value, bool) or not x.is_finite() or x <= 0:
        fail(f"{where} carries {SUM} that is not a positive number: {value!r}")
    return x


def half_up(x: Decimal) -> int:
    return int(x.quantize(Decimal(1), rounding=ROUND_HALF_UP))


def count_word(n: int, what: str) -> str:
    if not 0 < n < len(NUMBER_WORDS):
        fail(f"no word for a count of {n} {what}")
    return NUMBER_WORDS[n]


def and_list(items: list[str]) -> str:
    if len(items) <= 1:
        return "".join(items)
    return f"{', '.join(items[:-1])} and {items[-1]}"


def markdown_text(text: str) -> str:
    """Text from a record as a Markdown sentence shows it: every markup character behind a backslash."""
    return MARKUP.sub(r"\\\1", text)


def current_release(readme: Path) -> str:
    """The version the README's current-release heading names."""
    if not readme.is_file():
        fail(f"no {readme}")
    found = HEADING.findall(readme.read_text(encoding="utf-8"))
    if len(found) != 1:
        fail(f"{readme.name} has {len(found)} '## Current release: vX.Y.Z' headings; one is read")
    return found[0]


def load_json(path: Path, name: str):
    try:
        return json.loads(path.read_text(encoding="utf-8"), parse_float=Decimal)
    except (OSError, ValueError):
        fail(f"{name} is not readable JSON")


def release_records(repo: Path, version: str) -> dict[str, tuple[str, dict] | None]:
    """The records of the provenance directory that state the release, by kind, each as (path under the
    repository, data). One population record and one release record are required; the plugin record is
    None when the tree holds none; two of a kind are refused by name."""
    found: dict[str, list[tuple[str, dict]]] = {kind: [] for kind in KIND_NAMES}
    for path in sorted((repo / RECORDS).glob("*.json")):
        name = f"{RECORDS}/{path.name}"
        data = load_json(path, name)
        if not isinstance(data, dict) or data.get("release") != f"v{version}":
            continue
        named = data.get("record")
        kind = "plugins" if data.get("record_type") == "plugins" else KINDS.get(named) if isinstance(named, str) else None
        if kind is not None:
            found[kind].append((name, data))
    chosen: dict[str, tuple[str, dict] | None] = {}
    for kind, records in found.items():
        if len(records) > 1:
            fail(f"{len(records)} files of {RECORDS} are the {KIND_NAMES[kind]} of v{version}: "
                 + ", ".join(Path(name).name for name, _ in records))
        if not records and kind != "plugins":
            fail(f"no file of {RECORDS} is the {KIND_NAMES[kind]} of v{version}")
        chosen[kind] = records[0] if records else None
    return chosen


def read_registry(path: Path) -> dict[str, dict]:
    """The registry's tools by id, each {role, keys}; refused by name without its columns or with a role
    outside the registry's three."""
    if not path.is_file():
        fail(f"no registry {path}")
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    header = lines[0].split("\t") if lines else []
    missing = [column for column in ("id", "role", "record_keys") if column not in header]
    if missing:
        fail(f"{path.name} has no column {', '.join(missing)}")
    tools = {}
    for number, line in enumerate(lines[1:], start=2):
        row = dict(zip(header, (cell.strip() for cell in line.split("\t"))))
        if len(row) != len(header) or row["role"] not in ROLES:
            fail(f"{path.name} line {number} is no registry row with a role of {', '.join(ROLES)}")
        tools[row["id"]] = {"role": row["role"], "keys": tuple(row["record_keys"].split())}
    return tools


def dataset_sums(record: dict, name: str, engine: str, datasets: list[str]) -> dict[tuple[str, str], Decimal]:
    """Per (dataset, architecture), the engine's `sum_of_cell_medians_sec`; a pair the record times twice,
    or not at all, is refused by name."""
    where = f"{name} wall_time.{engine}.per_dataset"
    rows = field(field(field(record, "wall_time", name), engine, f"{name} wall_time"), "per_dataset", f"{name} wall_time.{engine}")
    if not isinstance(rows, list):
        fail(f"{where} is not a list")
    sums: dict[tuple[str, str], Decimal] = {}
    for index, row in enumerate(rows):
        row_where = f"{where}[{index}]"
        pair = (str(field(row, "dataset", row_where)), str(field(row, "arch", row_where)))
        if pair in sums:
            fail(f"{where} times {pair[0]} on {pair[1]} twice")
        sums[pair] = seconds(field(row, SUM, row_where), row_where)
    for dataset in datasets:
        for arch in ARCHES:
            if (dataset, arch) not in sums:
                fail(f"{where} has no row of {dataset} on {arch}")
    return sums


def row_key(engine: str, dataset: str, arch: str) -> str:
    return f"wall_time.{engine}.per_dataset[dataset={dataset}, arch={arch}].{SUM}"


def population_places(name: str, record: dict) -> tuple[dict[str, tuple], int]:
    """The places the population record fills, each (text, file, key, rule), and its count of datasets."""
    entries = field(record, "datasets", name)
    if not isinstance(entries, list) or not entries:
        fail(f"{name} lists no datasets")
    datasets = [str(field(entry, "dataset", f"{name} datasets[{index}]")) for index, entry in enumerate(entries)]
    if len(set(datasets)) != len(datasets):
        fail(f"{name} names a dataset twice")
    engine = dataset_sums(record, name, ENGINE, datasets)
    reference = dataset_sums(record, name, REFERENCE_ENGINE, datasets)
    with localcontext() as context:
        context.prec = 50
        quotient, dataset, arch = max((reference[pair] / engine[pair], *pair) for pair in ((d, a) for d in datasets for a in ARCHES))
    places = {"speed-up": (
        f"{half_up(quotient):,}", name, f"{row_key(REFERENCE_ENGINE, dataset, arch)} / {row_key(ENGINE, dataset, arch)}",
        f"the largest quotient over the datasets and {' and '.join(ARCHES)}, {quotient:.3f}, rounded half up, thousands separators",
    )}
    if GNOMAD not in datasets:
        fail(f"{name} datasets holds no {GNOMAD}")
    index = datasets.index(GNOMAD)
    where = f"{name} datasets[{index}]"
    selector = f"datasets[dataset={GNOMAD}]"
    label = str(field(entries[index], "label", where)).replace(WHOLE_GENOME_SUFFIX, "")
    variants = integer(field(entries[index], "variants", where), "variants", where)
    arm = engine[(GNOMAD, ARM)]
    places["gnomad-label"] = (markdown_text(label), name, f"{selector}.label", f"without '{WHOLE_GENOME_SUFFIX}'")
    places["gnomad-variants"] = (f"{variants:,}", name, f"{selector}.variants", "thousands separators")
    places["gnomad-arm-seconds"] = (
        f"{half_up(arm):,}", name, row_key(ENGINE, GNOMAD, ARM), f"{arm} rounded half up to whole seconds, thousands separators",
    )
    return places, len(datasets)


def release_places(name: str, record: dict, population_datasets: int) -> dict[str, tuple]:
    """The places the release record fills: the reference release, and the count of whole-genome datasets
    with the population record's."""
    families = field(record, "perl_reference", name)
    if not isinstance(families, dict) or not families:
        fail(f"{name} perl_reference names no reference")
    releases = []
    for family, block in families.items():
        image = block.get("image") if isinstance(block, dict) else None
        if not isinstance(image, str) or ":" not in image:
            fail(f"{name} perl_reference.{family} names no image tag")
        release = image.rsplit(":", 1)[1].removeprefix("release_")
        if release not in releases:
            releases.append(release)
    if len(releases) != 1:
        fail(f"{name} perl_reference names the releases {and_list(releases)}; the page names one")
    cells = field(record, "concordance", name)
    missing = [cell for cell in CLINVAR_CELLS if not isinstance(cells, dict) or cell not in cells]
    if missing:
        fail(f"{name} concordance holds no cell {', '.join(missing)}")
    count = population_datasets + len(CLINVAR_CELLS)
    return {
        "reference-release": (markdown_text(releases[0]), name, "perl_reference.*.image", "the tag without 'release_'"),
        "whole-genome-datasets": (
            count_word(count, "whole-genome datasets"), name, f"concordance.{', '.join(CLINVAR_CELLS)}",
            f"the population record's {population_datasets} datasets and these cells, as a word",
        ),
    }


def plugin_places(name: str, record: dict, tools: dict[str, dict]) -> dict[str, tuple]:
    """The place the plugin record fills: the plugins the released tool reads."""
    keys = next((tool["keys"] for tool in tools.values() if tool["role"] == "released"), ())
    blocks = field(record, "support", name)
    found = [block for block in blocks if isinstance(block, dict) and block.get("engine") in keys] if isinstance(blocks, list) else []
    if len(found) != 1:
        fail(f"{name} support carries {len(found)} blocks of the released tool; one is read")
    where = f"{name} support[engine={found[0]['engine']}]"
    read = integer(field(found[0], "read", where), "read", where)
    return {"plugins-read": (f"{read:,}", name, f"support[engine={found[0]['engine']}].read", "thousands separators")}


def figure_places(
    repo: Path, version: str, tools: dict[str, dict], clauses: dict[str, str], registry: str,
) -> tuple[dict[str, tuple], int] | None:
    """The places the figure's values file fills and its count of comparators, or None when the tree holds
    neither file of the pair. `registry` is the registry's file name, for the messages."""
    present = [name for name in (IMAGE, VALUES) if (repo / name).is_file()]
    if not present:
        return None
    if len(present) == 1:
        fail(f"the tree holds {present[0]} without {VALUES if present[0] == IMAGE else IMAGE}; the page shows the figure from both")
    document = load_json(repo / VALUES, VALUES)
    if field(document, "image", VALUES) != IMAGE:
        fail(f"{VALUES} names the image {document['image']!r}; the page shows {IMAGE}")
    labels: dict[str, list[str]] = {role: [] for role in ROLES}
    released_version = None
    entries = field(document, "tools", VALUES)
    for index, entry in enumerate(entries if isinstance(entries, list) else []):
        where = f"{VALUES} tools[{index}]"
        tool_id, name, tool_version = (str(field(entry, key, where)) for key in ("id", "name", "version"))
        if tool_id not in tools:
            fail(f"{where} names {tool_id}, which {registry} does not list")
        role = tools[tool_id]["role"]
        labels[role].append(markdown_text(f"{name} {tool_version}"))
        if role == "released":
            released_version = tool_version
    for role in ("released", "reference"):
        if len(labels[role]) != 1:
            fail(f"{VALUES} tools names {len(labels[role])} tools of the role {role}; the page names one")
    if released_version != f"v{version}":
        fail(f"{VALUES} draws {labels['released'][0]}; the README names v{version}")
    panels = field(document, "panels", VALUES)
    ids = [str(field(panel, "id", f"{VALUES} panels[{index}]")) for index, panel in enumerate(panels if isinstance(panels, list) else [])]
    if not ids:
        fail(f"{VALUES} holds no panel")
    for panel_id in ids:
        if panel_id not in clauses:
            fail(f"{VALUES} holds the panel {panel_id}, and {TEMPLATE.name} has no '<<panel {panel_id}>>' clause")
    notes = document.get("notes") or []
    if not isinstance(notes, list) or not all(isinstance(note, str) and not CONTROL.search(note) for note in notes):
        fail(f"{VALUES} notes is not a list of lines")
    lines = list(notes)
    for index, panel in enumerate(panels):
        note = panel.get("note")
        if note is None:
            continue
        if not isinstance(note, str) or CONTROL.search(note):
            fail(f"{VALUES} panels[{index}] note is not a line")
        lines.append(note)
    word = count_word(len(ids), "panels").capitalize()
    by_role = f"tools, by the role of each id in {registry}"
    places = {
        "figure-notes": (
            "".join(f"{markdown_text(line)} " for line in lines), VALUES, "notes, panels[].note",
            "each line followed by a space, the file's notes before the panels'",
        ),
        "released": (labels["released"][0], VALUES, by_role, "name and version"),
        "reference": (labels["reference"][0], VALUES, by_role, "name and version"),
        "comparators": (and_list(labels["comparator"]), VALUES, by_role, "names and versions"),
        "compared": (and_list(labels["released"] + labels["comparator"]), VALUES, by_role, "the released tool, then the comparators"),
        "panel-count": (f"{word} panel{'' if len(ids) == 1 else 's'}", VALUES, "panels", "the count as a word"),
        "panels": (and_list([clauses[panel_id] for panel_id in ids]), VALUES, "panels[].id", f"the clause {TEMPLATE.name} gives each id"),
    }
    return places, len(labels["comparator"])


def read_template(text: str) -> tuple[dict[str, str], str]:
    """The template's panel clauses by panel id, and its page."""
    head, mark, page = text.partition(f"{PAGE_MARK}\n")
    if not mark:
        fail(f"{TEMPLATE.name} has no '{PAGE_MARK}' line")
    clauses: dict[str, str] = {}
    for line in head.splitlines():
        found = PANEL_LINE.match(line)
        if line.strip() and not found:
            fail(f"{TEMPLATE.name} holds a line above '{PAGE_MARK}' that is no '<<panel ID>>CLAUSE': {line!r}")
        if found:
            clauses[found.group(1)] = found.group(2).strip()
    return clauses, page


def fill(page: str, parts: dict[str, bool], places: dict[str, tuple]) -> str:
    """The page with each part kept or dropped and each place filled. Refused by name: a mark that is
    neither a part nor a place, wherever it stands, a part without its pair, and a place without a value in
    this tree."""
    for mark in MARK.finditer(page):
        kind, _, name = mark.group(1).partition(" ")
        if kind in ("if", "end") and name:
            if name not in parts:
                fail(f"{TEMPLATE.name} names the part '{mark.group(0)}', which the script does not know")
        elif mark.group(1) not in places:
            fail(f"{TEMPLATE.name} holds the mark '{mark.group(0)}', which the script does not fill")

    def part(match: re.Match) -> str:
        return match.group(2) if parts[match.group(1)] else ""

    def place(match: re.Match) -> str:
        name = match.group(1)
        if name not in places:
            fail(f"{TEMPLATE.name} holds '{match.group(0)}' without its pair")
        if places[name] is None:
            fail(f"{TEMPLATE.name} holds the place '<<{name}>>' outside its part, and this tree gives it no value")
        return places[name][0]

    for pattern in (BLOCK_PART, INLINE_PART):
        previous = None
        while previous != page:
            previous, page = page, pattern.sub(part, page)
    return MARK.sub(place, page)


def build(repo: Path, template: str, registry: Path) -> tuple[str, dict[str, bool], dict[str, tuple], str]:
    """The page for the tree under `repo`, the parts that hold, every place with its source, and the
    release the README names."""
    version = current_release(repo / "README.md")
    clauses, page = read_template(template)
    tools = read_registry(registry)
    records = release_records(repo, version)
    places: dict[str, tuple | None] = dict.fromkeys(
        ("plugins-read", "figure-notes", "released", "reference", "comparators", "compared", "panel-count", "panels"),
    )
    population, datasets = population_places(*records["population"])
    places.update(population)
    places.update(release_places(*records["release"], datasets))
    if records["plugins"] is not None:
        places.update(plugin_places(*records["plugins"], tools))
    figure = figure_places(repo, version, tools, clauses, registry.name)
    comparators = 0
    if figure is not None:
        places.update(figure[0])
        comparators = figure[1]
    parts = {
        "plugins": records["plugins"] is not None, "figure": figure is not None,
        "one-comparator": comparators == 1, "several-comparators": comparators > 1,
    }
    return fill(page, parts, places), parts, places, version


def values_table(parts: dict[str, bool], places: dict[str, tuple | None]) -> str:
    rows = ["place\tvalue\tfile\tkey\trule"]
    rows += [f"if {name}\t{'holds' if holds else 'left out'}\t\t\t" for name, holds in parts.items()]
    rows += ["\t".join((name, *value)) for name, value in places.items() if value is not None]
    return "\n".join(rows) + "\n"


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description="Write docs/index.md from its template and the release's records.")
    ap.add_argument("--repo", type=Path, default=HERE.parents[1], help="repository root")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="write nothing; exit 1 when the page differs")
    mode.add_argument("--values", action="store_true", help="write nothing; print each place's value and source")
    args = ap.parse_args(argv)
    if not TEMPLATE.is_file():
        fail(f"no template {TEMPLATE}")
    page, parts, places, version = build(args.repo, TEMPLATE.read_text(encoding="utf-8"), REGISTRY)
    out = args.repo / PAGE
    if args.values:
        sys.stdout.write(values_table(parts, places))
        return
    if args.check:
        held = out.read_bytes().splitlines(keepends=True) if out.is_file() else []
        written = page.encode("utf-8").splitlines(keepends=True)
        for number in range(max(len(held), len(written))):
            if held[number:number + 1] != written[number:number + 1]:
                fail(f"{PAGE} line {number + 1} is not the line the records of v{version} give; "
                     "scripts/docs/render_home_page.py writes the page")
        print(f"{PAGE} is the page the records of v{version} give")
        return
    out.write_bytes(page.encode("utf-8"))
    print(f"wrote {out} ({len(page.encode('utf-8'))} bytes) for v{version}")


if __name__ == "__main__":
    main()
