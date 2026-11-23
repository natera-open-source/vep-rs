#!/usr/bin/env python3
"""Check the entry tables of docs/intended-divergences.md against a table of trigger frequencies.

usage: check_divergence_order.py [--reference VERSION] [--reproductions DIR] PAGE TABLE

The page's entries are the rows of its entry tables, the Markdown tables whose
header is ``| Entry | In one sentence | Tier and count | Verified against | Page |``;
each sits under a category heading, a ``##`` heading or a ``###`` heading inside
one, and the first cell of a row links the entry's key in code font. The
``Verified against`` cell names the Ensembl VEP release the entry's reproduction
was last run against, as a release number such as ``116.3``; with
``--reference`` every cell must name that release. The table is
tab-separated, ``#`` lines being comments, with at least the columns ``entry``,
``tier``, ``records_total``, ``records_with_shape``, ``per_million``,
``corpora_with_shape``, ``flag_sets_reaching`` and ``basis``. A row whose
``basis`` starts ``output rows over the population cuts`` counts the records
of the cuts with a triggering output row of the entry's stated command, and
its cell states the rate as ``per million by output rows``.

The section headed ``Differences only malformed input reaches`` holds the
entries whose input class is ``malformed``, under headings that repeat those of
sections above it one level down: a category is a ``###`` heading there and a
subsection of it a ``####`` heading; a ``####`` line anywhere else on the page
is body text. An entry's input class is the ``input_class`` of its record
``<entry>.toml`` under ``--reproductions`` (default: ``reproductions/`` beside
this script).

The check passes when the keys of the entry tables equal the table's ``entry``
column one to one, when each row's ``Tier and count`` cell is the rendering
``render_count`` gives its table row, when the rows under each heading stand in
the table's order (tier first; within tier 1 the rate per million falling,
within tier 2 the flag sets reaching the shape then the corpora carrying it
falling, then the key in byte order; within tier 3 the key in byte order), when
the subsections of a category that list entries stand in alphabetical order,
letter case ignored, after a first subsection named ``Several plugins``, each
headed once, and when no category lists entries both directly and in
subsections. The same
rules hold inside the malformed-input section, which lists its entries after
every other entry table and under its ``###`` headings alone. Each of those
headings, with entries or without, repeats the heading of a ``##`` section
above it, once and in the order they have there, and lists its entries as that
section does, directly or in subsections; every ``####`` heading under it,
with entries or without, repeats the heading of a ``###`` subsection of that
section. Every ``malformed`` entry
stands in that section with tier 3, no ``observed`` or ``unobserved`` entry
stands there, and every tier 1 or tier 2 entry is ``observed``; a tier 3 entry
may be any of the three. The page heads that section once at most; its opening
prose, the lines above its first ``###`` heading, writes no entry of the
section in code font, and every key it links there is an entry of the tables
outside it.
Every violation is printed on a line of its own; the exit status is 1 when
there is any and 0 otherwise. Stdlib only; Python 3.11 or later for
``tomllib``.
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

REPRODUCTIONS = Path(__file__).resolve().parent / "reproductions"
INPUT_CLASSES = ("observed", "unobserved", "malformed")
# The heading of the section that holds the `malformed` entries, without any `{#id}` attribute.
MALFORMED_SECTION = "Differences only malformed input reaches"
HEADING = re.compile(r"^(?P<marks>#{2,4}) (?P<title>.+)$")
HEADING_ID = re.compile(r"\s*\{#[^}]*\}$")
ENTRY_TABLE_CELLS = ["Entry", "In one sentence", "Tier and count", "Verified against", "Page"]
ENTRY_TABLE_HEADER = "| " + " | ".join(ENTRY_TABLE_CELLS) + " |"
ENTRY_CELL = re.compile(r"^\[`(?P<key>[a-z0-9_]+)`\]\([^)]+\)$")
# A key in code font in prose, linked as the first cell of its row writes it, or bare.
NAMED_KEY = re.compile(r"(?P<link>\[)?`(?P<key>[a-z0-9_]+)`(?(link)\]\([^)]+\))")
COUNT_CELL = ENTRY_TABLE_CELLS.index("Tier and count")
VERIFIED_CELL = ENTRY_TABLE_CELLS.index("Verified against")
# The `Verified against` cell: an Ensembl VEP release number.
RELEASE = re.compile(r"^[0-9]+(?:\.[0-9]+)*$")
TIERS = (1, 2, 3)
# A rendered `Tier and count` cell: the tier, the figure, and the flag sets (or tool runs) and corpora that reach it.
CELL = re.compile(r"^(?P<tier>[1-3]) · (?P<figure>[^·]+?) · (?P<tail>.+)$")
RATE = re.compile(r"^(?P<rate>[0-9][0-9,]*(?:\.[0-9]+)?) per million(?: by output rows)?$")
TAIL = re.compile(r"^(?P<flag_sets>[0-9][0-9,]*) (?:flag sets?|tool runs?), (?P<corpora>[0-9][0-9,]*) corp(?:us|ora)$")
TABLE_COLUMNS = (
    "entry",
    "tier",
    "records_total",
    "records_with_shape",
    "per_million",
    "corpora_with_shape",
    "flag_sets_reaching",
    "basis",
)
POPULATION_BASIS = "population cut records"
# A shape that is a property of the annotation is counted from the output rows of the entry's stated command over the same cuts, and its rate says so.
OUTPUT_ROW_BASIS = "output rows over the population cuts"
# The committed runs of `vep filter` and `vep recode` are not flag sets, so a tool entry's tail counts tool runs.
TOOL_RUN_BASIS = "committed runs of the"
# The unit a row's figure counts, by the start of the table's `basis` text, for a shape counted outside the population cuts.
UNITS = (
    ("corpora records on MT", "MT records"),
    ("corpora records on Y", "Y records"),
    ("corpora records with sample columns", "records with sample columns"),
    ("VCF files", "VCF files"),
    ("corpora ensembl-format lines", "Ensembl-format lines"),
    ("corpora hgvs-format lines", "HGVS lines"),
    ("corpora spdi-format lines", "SPDI lines"),
    ("corpora region-format lines", "region lines"),
    ("corpora id-format lines", "identifier lines"),
    ("rows of the committed tables", "table rows"),
    ("transcripts of the committed caches", "cache transcripts"),
    ("regulatory features of the committed caches", "regulatory features"),
    ("epigenome names of the committed regulatory caches", "epigenome names"),
    ("committed flag sets", "flag sets"),
    (TOOL_RUN_BASIS, "tool runs"),
    ("BGZF blocks", "BGZF blocks"),
    ("compressed files", "compressed files"),
)
FIRST_SUBSECTION = "Several plugins"


@dataclass
class Row:
    key: str
    cells: list[str]
    line: int


@dataclass
class Section:
    title: str
    line: int
    rows: list[Row] = field(default_factory=list)
    subsections: list[Section] = field(default_factory=list)
    # The keys in code font in the opening prose of the malformed-input section, each with its line and whether it is linked; empty for every other section.
    named: list[tuple[str, int, bool]] = field(default_factory=list)


def split_cells(line: str) -> list[str]:
    """The cells of a Markdown table row, a `\\|` inside a cell kept as a pipe."""
    inner = re.sub(r"(?<!\\)\|$", "", line.strip().removeprefix("|"))
    return [c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", inner)]


def heading_text(title: str) -> str:
    """A heading's words, without a trailing `{#id}` attribute."""
    return HEADING_ID.sub("", title)


def is_malformed_section(section: Section) -> bool:
    return heading_text(section.title) == MALFORMED_SECTION


def walk(section: Section) -> Iterator[Section]:
    """A section and every subsection below it, in page order."""
    yield section
    for sub in section.subsections:
        yield from walk(sub)


def leaves(sections: list[Section]) -> list[tuple[str, Section]]:
    """Every section that lists entries, in page order, with the name a message gives it: its heading, led by the malformed-input section's when it stands inside that section."""
    out = []
    for top in sections:
        lead = f"{MALFORMED_SECTION}, " if is_malformed_section(top) else ""
        for section in walk(top):
            if section.rows:
                name = section.title if section is top else lead + section.title
                out.append((name, section))
    return out


def class_counts(classes: Iterable[str]) -> str:
    """How many entries name each input class, as a summary line states it."""
    named = list(classes)
    return ", ".join(f"{named.count(input_class)} {input_class}" for input_class in INPUT_CLASSES)


def read_classes(directory: Path = REPRODUCTIONS) -> dict[str, str]:
    """The input class of every record of the reproduction registry, by entry key; a record that is not TOML or names none of the three classes is an error."""
    classes: dict[str, str] = {}
    for path in sorted(directory.glob("*.toml")):
        try:
            with open(path, "rb") as fh:
                input_class = tomllib.load(fh).get("input_class")
        except (tomllib.TOMLDecodeError, UnicodeDecodeError) as exc:
            raise ValueError(f"{path}: {exc}") from exc
        if input_class not in INPUT_CLASSES:
            raise ValueError(f"{path}: input_class is {input_class!r}, not one of {', '.join(INPUT_CLASSES)}")
        classes[path.stem] = input_class
    return classes


def parse_page(text: str) -> list[Section]:
    """The `##` sections of the page with their entry rows and `###` subsections and, inside the malformed-input section, the keys its own opening prose writes in code font and the `####` subsections of each `###` category; a table without a heading above it is an error.

    An entry table opens on its header row and runs to the next blank line or section heading. A `####` line outside the malformed-input section is body text: it opens no section, so the rows under it stay with the section open above it, and directly under a row it is refused as a row.
    """
    sections: list[Section] = []
    # The open `###` section, which a `####` heading of the malformed-input section goes under.
    category: Section | None = None
    current: Section | None = None
    in_table = False
    for n, line in enumerate(text.splitlines(), 1):
        heading = HEADING.match(line)
        level = len(heading.group("marks")) if heading else 0
        if level == 4 and not (sections and is_malformed_section(sections[-1])):
            level = 0
        if level:
            current = Section(heading.group("title").strip(), n)
            if level == 2:
                sections.append(current)
                category = None
            elif level == 3:
                if not sections:
                    raise ValueError(f"line {n}: a subsection before any section")
                sections[-1].subsections.append(current)
                category = current
            elif category is None:
                raise ValueError(f"line {n}: a `####` heading outside a `###` subsection")
            else:
                category.subsections.append(current)
            in_table = False
        elif not in_table and line.lstrip(" ").startswith("|") and len(line) - len(line.lstrip(" ")) <= 3 and split_cells(line) == ENTRY_TABLE_CELLS:
            if current is None:
                raise ValueError(f"line {n}: an entry table before any heading")
            in_table = True
        elif in_table and line.strip():
            # A table runs to the next blank line; every line until then is a row, however it is indented and with or without its outer pipes.
            cells = split_cells(line)
            if all(set(c) <= {"-", ":"} and "-" in c for c in cells):
                continue
            if len(cells) != len(ENTRY_TABLE_CELLS):
                raise ValueError(f"line {n}: an entry row with {len(cells)} cells")
            m = ENTRY_CELL.match(cells[0])
            if not m:
                raise ValueError(f"line {n}: the first cell is not a linked key in code font: {cells[0]!r}")
            current.rows.append(Row(m.group("key"), cells, n))
        else:
            in_table = False
            # `current` is the `##` section until its first `###` heading opens, so a subsection that repeats the title names nothing.
            if sections and current is sections[-1] and is_malformed_section(current):
                current.named += [(m.group("key"), n, bool(m.group("link"))) for m in NAMED_KEY.finditer(line)]
    return sections


def read_table(path: str) -> dict[str, dict[str, str]]:
    """The table's rows by entry key; the header is the first line that is not a comment."""
    rows: dict[str, dict[str, str]] = {}
    header: list[str] | None = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                header = line.split("\t")
                missing = [c for c in TABLE_COLUMNS if c not in header]
                if missing:
                    raise ValueError(f"{path}: the header lacks {', '.join(missing)}")
                continue
            cells = line.split("\t")
            if len(cells) != len(header):
                raise ValueError(f"{path}: a row with {len(cells)} cells under a header of {len(header)}")
            row = dict(zip(header, cells))
            if row["entry"] in rows:
                raise ValueError(f"{path}: {row['entry']} appears twice")
            rows[row["entry"]] = row
    if header is None:
        raise ValueError(f"{path}: no header line")
    return rows


def format_rate(per_million: float) -> str:
    if per_million >= 1000:
        return f"{per_million:,.0f}"
    if per_million >= 1:
        return f"{per_million:.1f}"
    return f"{per_million:.2g}"


def plural(n: int, one: str, many: str) -> str:
    return f"{n:,} {one if n == 1 else many}"


def render_count(row: dict[str, str]) -> str:
    """The `Tier and count` cell of a table row: the tier, the figure that placed it, and the flag sets (tool runs for a tool entry) and corpora that reach it."""
    tier = tier_of(row)
    total = int(row["records_total"])
    with_shape = int(row["records_with_shape"])
    corpora = int(row["corpora_with_shape"])
    flag_sets = int(row["flag_sets_reaching"])
    basis = row["basis"]
    by_output_rows = basis.startswith(OUTPUT_ROW_BASIS)
    if by_output_rows or basis.startswith(POPULATION_BASIS):
        figure = f"{format_rate(float(row['per_million']))} per million" if with_shape else "none on the cuts"
        if by_output_rows and with_shape:
            figure += " by output rows"
    else:
        unit = next((u for prefix, u in UNITS if basis.startswith(prefix)), None)
        if unit is None:
            raise ValueError(f"{row['entry']}: no unit for the basis {basis!r}")
        figure = f"{with_shape:,} of {total:,} {unit}" if total else f"no {unit} to count"
    if basis.startswith(TOOL_RUN_BASIS):
        reaching = plural(flag_sets, "tool run", "tool runs")
    else:
        reaching = plural(flag_sets, "flag set", "flag sets")
    tail = "no flag needed" if tier == 1 else f"{reaching}, {plural(corpora, 'corpus', 'corpora')}"
    return f"{tier} · {figure} · {tail}"


def tier_of(row: dict[str, str]) -> int:
    tier = int(row["tier"])
    if tier not in TIERS:
        raise ValueError(f"{row['entry']}: tier {row['tier']} is not 1, 2 or 3")
    return tier


def order_key(row: dict[str, str]) -> tuple:
    """The rule's sort key: tier 1 by rate, tier 2 by flag sets then corpora then key, tier 3 by key; keys compare in byte order."""
    tier = tier_of(row)
    if tier == 1:
        return (1, -float(row["per_million"]), row["entry"])
    if tier == 2:
        return (2, -int(row["flag_sets_reaching"]), -int(row["corpora_with_shape"]), row["entry"])
    return (3, 0, row["entry"])


def cell_order_key(key: str, cell: str) -> tuple:
    """The rule's sort key read from a rendered `Tier and count` cell, for a check of the page alone."""
    m = CELL.match(cell)
    if not m:
        raise ValueError(f"{key}: the cell {cell!r} is not a rendered tier and count")
    tier = int(m.group("tier"))
    if tier == 1:
        rate = RATE.match(m.group("figure"))
        if not rate:
            raise ValueError(f"{key}: a tier 1 cell states no rate per million: {cell!r}")
        # The rendered rate is rounded, so entries it cannot tell apart keep the page's order.
        return (1, -float(rate.group("rate").replace(",", "")))
    if tier == 2:
        tail = TAIL.match(m.group("tail"))
        if not tail:
            raise ValueError(f"{key}: a tier 2 cell names no flag sets and corpora: {cell!r}")
        return (2, -int(tail.group("flag_sets").replace(",", "")), -int(tail.group("corpora").replace(",", "")), key)
    return (3, 0, key)


def category_errors(category: Section) -> list[str]:
    """A category's subsections stand in alphabetical order after `Several plugins`, each headed once, and it lists its entries directly or in subsections, never both."""
    errors = []
    titles = [heading_text(section.title) for section in category.subsections if section.rows]
    # A subsection headed twice splits its table, and the order rule reads each part alone.
    if len(set(titles)) != len(titles):
        errors.append(f"line {category.line}: {category.title!r} lists a subsection twice")
    if FIRST_SUBSECTION in titles[1:]:
        errors.append(f"line {category.line}: {FIRST_SUBSECTION!r} is not the first subsection of {category.title!r}")
    elif titles and titles[0] == FIRST_SUBSECTION:
        titles = titles[1:]
    if titles != sorted(titles, key=str.casefold):
        errors.append(f"line {category.line}: the subsections of {category.title!r} are not in alphabetical order")
    if category.rows and any(section.rows for section in category.subsections):
        errors.append(f"line {category.line}: {category.title!r} lists entries both directly and in subsections")
    return errors


def check_layout(sections: list[Section]) -> list[str]:
    """Every violation of the page's layout: the subsections of each category, and the malformed-input section, which the page heads once at most and which lists its entries after every other entry table and under its `###` headings alone; each of those headings, with entries or without, repeats the heading of a `##` section above it, once and in the order they have there, and lists its entries as that section does, directly or in subsections; every `####` heading under it, with entries or without, repeats the heading of a `###` subsection of that section."""
    headed = [top for top in sections if is_malformed_section(top)]
    malformed = headed[0] if headed else None
    errors = [error for top in sections if top is not malformed for error in category_errors(top)]
    if malformed is None:
        return errors
    errors += [f"line {top.line}: {MALFORMED_SECTION!r} is headed twice" for top in headed[1:]]
    place = sections.index(malformed)
    if malformed.rows:
        errors.append(f"line {malformed.line}: {MALFORMED_SECTION!r} lists entries outside a category")
    above = {heading_text(top.title): top for top in sections[:place]}
    category_titles = [heading_text(category.title) for category in malformed.subsections]
    if len(set(category_titles)) != len(category_titles):
        errors.append(f"line {malformed.line}: {MALFORMED_SECTION!r} lists a category twice")
    for category, title in zip(malformed.subsections, category_titles):
        errors += category_errors(category)
        if title not in above:
            errors.append(f"line {category.line}: {category.title!r} is no category of the tables above {MALFORMED_SECTION!r}")
            continue
        grouped = any(section.rows for section in category.subsections)
        # A category that lists its entries directly above lists them directly here: a `####` heading would split its table and hide its order.
        if above[title].rows and grouped:
            errors.append(f"line {category.line}: {category.title!r} lists entries in subsections, which the category above does not have")
            continue
        if category.rows and not grouped and any(section.rows for section in above[title].subsections):
            errors.append(f"line {category.line}: {category.title!r} lists entries directly, which the category above lists in subsections")
        groups = {heading_text(section.title) for section in above[title].subsections}
        errors += [
            f"line {section.line}: {section.title!r} is no subsection of {category.title!r} in the tables above {MALFORMED_SECTION!r}"
            for section in category.subsections
            if heading_text(section.title) not in groups
        ]
    order = list(above)
    places = [order.index(title) for title in category_titles if title in above]
    if places != sorted(places):
        errors.append(f"line {malformed.line}: the categories of {MALFORMED_SECTION!r} do not stand in the order they have above it")
    for top in sections[place + 1 :]:
        if any(section.rows for section in walk(top)):
            errors.append(f"line {top.line}: {top.title!r} lists entries after {MALFORMED_SECTION!r}, the last section that does")
    return errors


def check_named(sections: list[Section]) -> list[str]:
    """Every key in code font in the opening prose of the malformed-input section that is an entry of that section, and every key linked there that is no entry of the tables outside it, as one line each."""
    malformed = next((top for top in sections if is_malformed_section(top)), None)
    if malformed is None:
        return []
    inside = {row.key for section in walk(malformed) for row in section.rows}
    outside = {row.key for top in sections if top is not malformed for section in walk(top) for row in section.rows}
    errors = []
    for key, line, linked in malformed.named:
        if key in inside:
            errors.append(f"line {line}: {key}, which the opening of {MALFORMED_SECTION!r} names as an entry of the tables above, stands in that section")
        elif linked and key not in outside:
            errors.append(f"line {line}: {key}, which the opening of {MALFORMED_SECTION!r} names, is no entry of the page")
    return errors


def check_classes(sections: list[Section], classes: dict[str, str]) -> list[str]:
    """Every violation of the input classes against the page's sections and tiers, as one line each: a tier 1 or tier 2 entry is `observed`, a `malformed` entry stands in the malformed-input section and no other entry does, and a tier 3 entry may be any of the three. A row whose key `classes` lacks is not judged, and one whose cell states no tier is held to its section alone."""
    errors = []
    for top in sections:
        inside = is_malformed_section(top)
        for section in walk(top):
            for row in section.rows:
                input_class = classes.get(row.key)
                if input_class is None:
                    continue
                cell = CELL.match(row.cells[COUNT_CELL])
                if cell and cell.group("tier") != "3" and input_class != "observed":
                    errors.append(f"line {row.line}: {row.key} is tier {cell.group('tier')} with input_class {input_class!r}: a tier 1 or tier 2 entry is observed")
                if input_class == "malformed" and not inside:
                    errors.append(f"line {row.line}: {row.key} is malformed and stands under {section.title!r}, outside {MALFORMED_SECTION!r}")
                if input_class != "malformed" and inside:
                    errors.append(f"line {row.line}: {row.key} stands under {MALFORMED_SECTION!r} with input_class {input_class!r}")
    return errors


def check(page_text: str, table: dict[str, dict[str, str]], reference: str | None = None, classes: dict[str, str] | None = None) -> list[str]:
    """Every violation of the page against the table, as one line each; with `reference`, every `Verified against` cell must name it, and with `classes`, the registry's input class of each entry, every entry must have one and stand where its class places it."""
    errors: list[str] = []
    try:
        sections = parse_page(page_text)
    except ValueError as exc:
        return [str(exc)]
    listing = leaves(sections)
    rows = [row for _, section in listing for row in section.rows]
    seen: set[str] = set()
    for row in rows:
        if row.key in seen:
            errors.append(f"line {row.line}: {row.key} is listed twice")
        seen.add(row.key)
        if row.key not in table:
            errors.append(f"line {row.line}: {row.key} has no row in the table")
    errors += [f"{key} is in the table and not on the page" for key in table if key not in seen]
    for row in rows:
        verified = row.cells[VERIFIED_CELL]
        if not RELEASE.match(verified):
            errors.append(f"line {row.line}: {row.key} is verified against {verified!r}, not a release number")
        elif reference is not None and verified != reference:
            errors.append(f"line {row.line}: {row.key} is verified against {verified}, the reference is {reference}")
        if row.key not in table:
            continue
        rendered = render_count(table[row.key])
        if row.cells[COUNT_CELL] != rendered:
            errors.append(f"line {row.line}: {row.key} states {row.cells[COUNT_CELL]!r}, the table gives {rendered!r}")
    for name, section in listing:
        known = [row for row in section.rows if row.key in table]
        expected = sorted(known, key=lambda row: order_key(table[row.key]))
        for actual, wanted in zip(known, expected):
            if actual.key != wanted.key:
                errors.append(f"line {actual.line}: under {name!r} {actual.key} precedes {wanted.key}, which the rule places first")
                break
    errors += check_layout(sections)
    errors += check_named(sections)
    if classes is not None:
        errors += [f"line {row.line}: {row.key} has no input class in the registry" for row in rows if row.key not in classes]
        errors += check_classes(sections, classes)
    return errors


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("page", help="docs/intended-divergences.md")
    ap.add_argument("table", help="the tab-separated trigger-frequency table")
    ap.add_argument("--reference", help="the Ensembl VEP release every `Verified against` cell must name")
    ap.add_argument("--reproductions", type=Path, default=REPRODUCTIONS, help="the reproduction registry, whose records name each entry's input class")
    args = ap.parse_args(argv)
    try:
        with open(args.page, encoding="utf-8") as fh:
            page_text = fh.read()
        table = read_table(args.table)
        classes = read_classes(args.reproductions)
        errors = check(page_text, table, args.reference, classes)
    except (OSError, ValueError) as exc:
        print(f"ERROR: [divergence-order] {exc}", file=sys.stderr)
        return 1
    for err in errors:
        print(f"ERROR: [divergence-order] {err}", file=sys.stderr)
    if errors:
        return 1
    keys = [row.key for _, section in leaves(parse_page(page_text)) for row in section.rows]
    print(f"{len(keys)} entries match the table's {len(table)} rows and stand in its order: {class_counts(classes[key] for key in keys)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
