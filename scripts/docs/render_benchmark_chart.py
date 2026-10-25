#!/usr/bin/env python3
"""Render the current release's whole-genome wall-time chart as an SVG.

Input: the README's `## Current release: vX.Y.Z` heading names the release, and the
one `docs/concordance-provenance/<date>-population-vX.Y.Z.json` of that version is
its population record. The chart reads that record's `datasets` (the four whole-genome
population datasets) and, per engine under `wall_time` (`vep_rs`, `ensembl_vep`,
`fastvep`), the ARM `per_cell` wall-time medians, summed per dataset as the README's
tables sum them, plus the protocol fields `n_clones_per_cell` and `fork` and each
engine's `engine_version`. No figure is typed here; a missing or ambiguous record,
dataset or field is refused by name, and a dataset the record did not time with
fastVEP loses that bar with a printed note.

Output: `docs/assets/benchmark.svg`, horizontal bars on a logarithmic time axis, one
group per dataset, every bar labelled with its wall time, a legend and a title naming
the protocol. The same record gives byte-identical output.

Layouts: the root holds two, each a nested `<svg>` with its own viewBox. The narrow one
puts each dataset's name above its bars and shows below BREAKPOINT px of viewport width;
the wide one puts the names beside the bars and shows from BREAKPOINT px. An embedded
`<style>` shows one of them and gives the root that layout's aspect ratio; a style
element in an inlined SVG applies to the whole page, so its selectors name only the
chart's own classes. Every label is set at FONT units, which renders at 12 px or more
wherever the page shows its layout.

Page and file: the file is inlined into a Markdown page, so it has no blank line and no
colour of its own. Bars carry the classes `bar-vep-rs`, `bar-ensembl-vep` and
`bar-fastvep` for the page's stylesheet, and text and gridlines use the page's `--fg`.
The root's width and height are the wide layout's, the size the file takes where no page
sizes it, as when it is opened on its own or shown through an `<img>` element.

Staleness stamp: the root carries `data-release="vX.Y.Z"`, repeated in `<desc>`, with
the record's release, so the site build can refuse a chart whose release differs from
the README heading of the same tree. Stdlib only.

usage: render_benchmark_chart.py [--repo DIR] [--out FILE]
"""

from __future__ import annotations

import argparse
import html
import json
import math
import re
import sys
from decimal import ROUND_HALF_UP, Decimal, InvalidOperation
from pathlib import Path

ARCH = "arm64"
ARCH_LABEL = "ARM"
# Record key, display name and bar class, in bar and legend order.
SERIES = (
    ("vep_rs", "vep-rs", "bar-vep-rs"),
    ("fastvep", "fastVEP", "bar-fastvep"),
    ("ensembl_vep", "Ensembl VEP", "bar-ensembl-vep"),
)
OPTIONAL_ENGINES = {"fastvep"}
BASE_ENGINE = "vep_rs"
HEADING = re.compile(r"^## Current release: v(\d+\.\d+\.\d+)\s*$", re.MULTILINE)
WHOLE_GENOME_SUFFIX = ", whole genome"
MINUTES_FLOOR = Decimal(600)

CLASS = "benchmark-chart"
WIDE_WIDTH = 680
NARROW_WIDTH = 340
# Viewport width from which the page's text column is wide enough for the wide layout's
# labels to render at 12 px or more. Below it the root is capped, so a column wider than
# a phone does not scale the narrow layout's labels past the page's body text.
BREAKPOINT = 640
NARROW_MAX_WIDTH = 400
# Inset of the top, right and bottom edges. Text starts one unit from the left edge, in line with the
# page's text, so its anti-aliased first glyph is not clipped by the chart's edge.
PAD = 4
INSET = 1
BAR = 10
# The gap keeps a label's comma, which descends below the bar it is centred on, off the next bar.
GAP = 3
GROUP_GAP = 8
TIP_RADIUS = 4
FONT = 15
LINE = 17
TITLE_LINE = 19
SWATCH = 12
LEGEND_GAP = 20
VALUE_GAP = 5
LABEL_CHARS = 16
# Height of digits and capitals, and depth of descenders, as fractions of the font size.
CAP = 0.72
DESCENT = 0.24
# Width per character as a fraction of the font size, generous for the sans-serif faces a
# page may set, so a reserved margin fits a label the chart cannot measure.
EM = 0.62
SECONDARY = "0.75"
GRID_OPACITY = "0.25"


def fail(message: str) -> None:
    sys.exit(f"ERROR: [render_benchmark_chart] {message}")


def required(block, keys: tuple[str, ...], where: str) -> tuple:
    """The values under `keys` of a record block, refused by name when any is absent or the block is
    not an object."""
    values = block if isinstance(block, dict) else {}
    missing = [k for k in keys if values.get(k) is None]
    if missing:
        fail(f"{where} carries no {', '.join(missing)}")
    return tuple(values[k] for k in keys)


def integer(value, key: str, where: str) -> int:
    """A record's integer, refused by name when it is not one: `bool` is a subclass of `int`."""
    if isinstance(value, bool) or not isinstance(value, int):
        fail(f"{where} carries {key} that is not an integer: {value!r}")
    return value


def seconds(value, where: str) -> Decimal:
    """A record's `median_sec` as a Decimal, refused by name when it is not a finite positive number;
    a logarithmic axis has no place for zero."""
    try:
        x = Decimal(str(value))
    except InvalidOperation:
        fail(f"{where} carries a median_sec that is not a number: {value!r}")
    if not x.is_finite() or x <= 0:
        fail(f"{where} carries a median_sec that is not a positive number: {value!r}")
    return x


def current_release(readme: Path) -> str:
    """The version the README's current-release heading names; the heading is read once."""
    if not readme.exists():
        fail(f"no {readme}")
    found = HEADING.findall(readme.read_text(encoding="utf-8"))
    if len(found) != 1:
        fail(f"{readme} has {len(found)} '## Current release: vX.Y.Z' headings; one is read")
    return found[0]


def population_record(repo: Path, version: str) -> Path:
    """The one docs/concordance-provenance/<date>-population-v<version>.json; none or several fail."""
    pattern = f"*-population-v{version}.json"
    records = sorted((repo / "docs" / "concordance-provenance").glob(pattern))
    if not records:
        fail(f"no docs/concordance-provenance/{pattern}")
    if len(records) > 1:
        fail(f"{len(records)} records match docs/concordance-provenance/{pattern}: "
             + ", ".join(r.name for r in records))
    return records[0]


def load_record(path: Path, version: str) -> dict:
    record = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(record, dict):
        fail(f"{path.name} is not an object")
    (release,) = required(record, ("release",), path.name)
    if release != f"v{version}":
        fail(f"{path.name} records release {release!r}; the README names v{version}")
    return record


def read_datasets(record: dict, name: str) -> list[dict]:
    """The record's whole-genome datasets as {key, label, assembly, chromosomes}, in record order."""
    (entries,) = required(record, ("datasets",), name)
    if not isinstance(entries, list) or len(entries) != 4:
        fail(f"{name} lists {len(entries) if isinstance(entries, list) else 'no'} datasets; the chart draws four")
    out = []
    for index, entry in enumerate(entries):
        where = f"{name} datasets[{index}]"
        key, label, assembly, chromosomes = required(entry, ("dataset", "label", "assembly", "chromosomes"), where)
        if not isinstance(chromosomes, list) or not chromosomes:
            fail(f"{where} carries no chromosome list")
        if len(set(chromosomes)) != len(chromosomes):
            fail(f"{where} repeats a chromosome")
        out.append({"key": str(key), "label": str(label), "assembly": str(assembly),
                    "chromosomes": {str(c) for c in chromosomes}})
    if len({d["key"] for d in out}) != 4:
        fail(f"{name} names a dataset twice")
    return out


def engine_times(block, engine: str, datasets: list[dict], name: str) -> dict[str, Decimal]:
    """Per dataset key, the sum of the engine's ARM per-chromosome medians. A dataset is summed only when
    its cells cover its chromosome list exactly once each; a repeated chromosome is refused, a partial
    cover is refused, and no cover at all leaves the dataset out for the caller to judge."""
    where = f"{name} wall_time.{engine}"
    (cells,) = required(block, ("per_cell",), where)
    if not isinstance(cells, list):
        fail(f"{where}.per_cell is not a list")
    by_dataset: dict[str, dict[str, Decimal]] = {}
    for index, cell in enumerate(cells):
        cell_where = f"{where}.per_cell[{index}]"
        key, arch, chromosome, median = required(cell, ("dataset", "arch", "chromosome", "median_sec"), cell_where)
        if arch != ARCH:
            continue
        medians = by_dataset.setdefault(str(key), {})
        if str(chromosome) in medians:
            fail(f"{where} times chromosome {chromosome} of {key} on {ARCH} twice")
        medians[str(chromosome)] = seconds(median, cell_where)
    sums: dict[str, Decimal] = {}
    for dataset in datasets:
        medians = by_dataset.get(dataset["key"], {})
        if not medians:
            continue
        if set(medians) != dataset["chromosomes"]:
            expected = len(dataset["chromosomes"])
            fail(f"{where} covers {len(medians)} of {expected} chromosomes of {dataset['key']} on {ARCH}")
        sums[dataset["key"]] = sum(medians.values(), Decimal(0))
    return sums


def engine_label(block, engine: str, display: str, name: str) -> str:
    """The legend text: the display name and the record's engine version, a `release_` prefix dropped."""
    (version,) = required(block, ("engine_version",), f"{name} wall_time.{engine}")
    version = str(version)
    return f"{display} {version.removeprefix('release_')}"


def read_chart(record: dict, name: str, version: str) -> dict:
    """Everything the chart draws, read from the record: datasets, per-series labels and times, the base
    series' machine count, every series measured on another count, the thread count, and notes on
    omitted bars."""
    datasets = read_datasets(record, name)
    (wall_time, fork) = required(record, ("wall_time", "fork"), name)
    fork = integer(fork, "fork", name)
    series = []
    clones: dict[str, int] = {}
    notes: list[str] = []
    for engine, display, css_class in SERIES:
        where = f"{name} wall_time.{engine}"
        block = wall_time.get(engine) if isinstance(wall_time, dict) else None
        if block is None:
            if engine in OPTIONAL_ENGINES:
                notes.append(f"note: {name} has no wall_time.{engine}; the chart omits {display}")
                continue
            fail(f"{name} carries no wall_time.{engine}")
        times = engine_times(block, engine, datasets, name)
        missing = [d["key"] for d in datasets if d["key"] not in times]
        if missing and engine not in OPTIONAL_ENGINES:
            fail(f"{where} has no {ARCH} wall time for {', '.join(missing)}")
        for key in missing:
            notes.append(f"note: {name} has no {display} {ARCH} wall time for {key}; the chart omits that bar")
        if not times:
            notes.append(f"note: {name} times {display} on no dataset; the chart omits the series")
            continue
        (n_clones,) = required(block, ("n_clones_per_cell",), where)
        clones[engine] = integer(n_clones, "n_clones_per_cell", where)
        series.append({"engine": engine, "display": display, "css_class": css_class,
                       "label": engine_label(block, engine, display, name), "times": times})
    base_clones = clones[BASE_ENGINE]
    other_clones = [(s["display"], clones[s["engine"]]) for s in series if clones[s["engine"]] != base_clones]
    return {"version": version, "datasets": datasets, "series": series, "fork": fork,
            "n_clones": base_clones, "other_clones": other_clones, "notes": notes}


def significant_3(x: Decimal) -> str:
    """A ratio at three significant figures, rounded half up: 422, 98.0, 5.30."""
    if x >= 100:
        return str(x.quantize(Decimal(1), rounding=ROUND_HALF_UP))
    if x >= 10:
        return str(x.quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))
    return str(x.quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))


def ratio_clause(chart: dict, base: dict, series: dict) -> str:
    """`<low> to <high> times faster than <engine label>`: the range over the datasets the series was timed
    on, with their count when fewer than all. The base series covers every dataset, so the range is never
    empty."""
    ratios = sorted(series["times"][k] / base["times"][k] for k in series["times"])
    low, high = significant_3(ratios[0]), significant_3(ratios[-1])
    span = low if low == high else f"{low} to {high}"
    scope = "" if len(ratios) == len(chart["datasets"]) else f" on the {len(ratios)} datasets it was timed on"
    return f"{span} times faster than {series['label']}{scope}"


def headline(chart: dict) -> str:
    """The `aria-label`: the base series' speed-up range against every other series."""
    base = next(s for s in chart["series"] if s["engine"] == BASE_ENGINE)
    clauses = " and ".join(ratio_clause(chart, base, s) for s in chart["series"] if s is not base)
    return f"{base['display']} annotates the four whole-genome datasets {clauses} on {ARCH_LABEL}."


def display_unit(values: list[Decimal]) -> tuple[str, Decimal]:
    """Seconds, the record's unit, unless every value exceeds ten minutes, then minutes."""
    if min(values) > MINUTES_FLOOR:
        return "min", Decimal(60)
    return "s", Decimal(1)


def value_label(x: Decimal, unit: str) -> str:
    """Seconds, whole from 10 s, to one decimal from 1 s and to two below, or minutes to one decimal; rounded
    half up with thousands separators."""
    if unit == "s":
        places = Decimal(1) if x >= 10 else Decimal("0.1") if x >= 1 else Decimal("0.01")
    else:
        places = Decimal("0.1")
    return f"{x.quantize(places, rounding=ROUND_HALF_UP):,} {unit}"


def tick_label(decade: int) -> str:
    """A decade tick in short form: 10, 100, 1K, 10K, 100K, 1M."""
    for power, suffix in ((6, "M"), (3, "K")):
        if decade >= power:
            return f"{10 ** (decade - power)}{suffix}"
    return f"{10 ** decade}"


def words(content: str) -> list[str]:
    """The words of `content`, a number kept with the word after it."""
    units: list[str] = []
    for word in content.split():
        if units and units[-1].isdigit():
            units[-1] += " " + word
        else:
            units.append(word)
    return units


def wrap(content: str, max_chars: int) -> list[str]:
    """Greedy word wrap; a word longer than the limit stands on its own line."""
    lines: list[str] = []
    for word in words(content):
        if lines and len(lines[-1]) + 1 + len(word) <= max_chars:
            lines[-1] += " " + word
        else:
            lines.append(word)
    return lines


def balanced_wrap(content: str, max_chars: int) -> list[str]:
    """The wrap with as few lines as the greedy wrap at `max_chars` whose longest line is shortest."""
    count = len(wrap(content, max_chars))
    limit = max(len(word) for word in words(content))
    while len(wrap(content, limit)) > count:
        limit += 1
    return wrap(content, limit)


def est(content: str) -> float:
    return len(content) * EM * FONT


def chars_within(width: float) -> int:
    return int(width / (EM * FONT))


def num(x: float) -> str:
    return f"{x:.1f}".rstrip("0").rstrip(".")


def text(
    x: float, y: float, content: str, anchor: str = "start", weight: str = "", opacity: str = "", suffix: str = ""
) -> str:
    """A label in the page's text colour; `suffix`, when given, follows it in secondary ink."""
    attrs = [f'x="{num(x)}"', f'y="{num(y)}"', f'font-size="{FONT}"']
    if anchor != "start":
        attrs.append(f'text-anchor="{anchor}"')
    if weight:
        attrs.append(f'font-weight="{weight}"')
    if opacity:
        attrs.append(f'fill-opacity="{opacity}"')
    attrs.append('fill="var(--fg)"')
    body = html.escape(content, quote=False)
    if suffix:
        body += f' <tspan fill-opacity="{SECONDARY}">{html.escape(suffix, quote=False)}</tspan>'
    return f"<text {' '.join(attrs)}>{body}</text>"


def bar(x: float, y: float, width: float, css_class: str) -> str:
    """A bar square at the axis and rounded at its tip; one shorter than two radii is a plain rectangle."""
    if width <= 2 * TIP_RADIUS:
        return f'<rect class="{css_class}" x="{num(x)}" y="{num(y)}" width="{num(width)}" height="{BAR}"/>'
    r = TIP_RADIUS
    straight = num(width - r)
    return (
        f'<path class="{css_class}" d="M{num(x)} {num(y)}h{straight}a{r} {r} 0 0 1 {r} {r}'
        f"v{BAR - 2 * r}a{r} {r} 0 0 1 -{r} {r}h-{straight}z\"/>"
    )


def gridline(x: float, top: float, bottom: float) -> str:
    return (
        f'<line x1="{num(x)}" y1="{num(top)}" x2="{num(x)}" y2="{num(bottom)}" '
        f'stroke="var(--fg)" stroke-opacity="{GRID_OPACITY}" stroke-width="1"/>'
    )


def axis(chart: dict) -> dict:
    """The values in the unit shown, keyed by (engine, dataset key), and the decades the logarithmic axis
    spans, the first one strictly below the smallest value so every bar has a length."""
    values = [v for s in chart["series"] for v in s["times"].values()]
    unit, divisor = display_unit(values)
    shown = {(s["engine"], k): v / divisor for s in chart["series"] for k, v in s["times"].items()}
    smallest, largest = min(shown.values()), max(shown.values())
    low = math.floor(math.log10(smallest))
    if Decimal(10) ** low >= smallest:
        low -= 1
    high = max(math.ceil(math.log10(largest)), low + 1)
    return {"unit": unit, "shown": shown, "low": low, "high": high}


def fraction(ax: dict, value: Decimal) -> float:
    """Where a value sits along the axis: 0 at its first decade, 1 at its last."""
    return (math.log10(value) - ax["low"]) / (ax["high"] - ax["low"])


def decade_x(ax: dict, decade: int, left: float, plot: float) -> float:
    return left + (decade - ax["low"]) / (ax["high"] - ax["low"]) * plot


def plot_width(ax: dict, left: float, right: float) -> float:
    """The widest plot from `left` whose last tick label, and every value label set VALUE_GAP past its bar's
    tip, end by `right`."""
    limits = [right - left - est(tick_label(ax["high"])) / 2]
    for value in ax["shown"].values():
        limits.append((right - left - VALUE_GAP - est(value_label(value, ax["unit"]))) / fraction(ax, value))
    return min(limits)


def short_name(dataset: dict) -> str:
    return dataset["label"].replace(WHOLE_GENOME_SUFFIX, "")


def title_text(chart: dict) -> str:
    """The protocol line: the base series' machines per cell, with every series measured on another
    count named beside it."""
    machines = f"{chart['n_clones']} machines"
    if chart["other_clones"]:
        others = ", ".join(f"{display} over {count}" for display, count in chart["other_clones"])
        machines += f" ({others})"
    return f"Wall time, sum of per-chromosome medians over {machines} at {chart['fork']} threads, {ARCH_LABEL}"


def heading(chart: dict, width: float) -> tuple[list[str], float]:
    """The title in balanced lines, then the legend in as many rows as `width` needs; returns them and the
    baseline of the last legend row."""
    out = []
    lines = balanced_wrap(title_text(chart), chars_within(width - INSET - PAD))
    for index, line in enumerate(lines):
        out.append(text(INSET, PAD + FONT + index * TITLE_LINE, line, weight="600"))
    baseline = PAD + FONT + (len(lines) - 1) * TITLE_LINE + 22
    x = float(INSET)
    for s in chart["series"]:
        entry = SWATCH + 6 + est(s["label"])
        if x > INSET and x + entry > width - PAD:
            x = float(INSET)
            baseline += LINE + 5
        out.append(
            f'<rect class="{s["css_class"]}" x="{num(x)}" y="{num(baseline - CAP * FONT / 2 - SWATCH / 2)}" '
            f'width="{SWATCH}" height="{SWATCH}" rx="2"/>'
        )
        out.append(text(x + SWATCH + 6, baseline, s["label"]))
        x += entry + LEGEND_GAP
    return out, baseline


def bars(chart: dict, ax: dict, dataset: dict, left: float, plot: float, top: float) -> list[str]:
    """One dataset's bars from `top`, one row per series, each labelled at its tip; a series the record did
    not time on the dataset leaves its row empty."""
    out = []
    for row, s in enumerate(chart["series"]):
        value = ax["shown"].get((s["engine"], dataset["key"]))
        if value is None:
            continue
        y = top + row * (BAR + GAP)
        tip = left + fraction(ax, value) * plot
        out.append(bar(left, y, tip - left, s["css_class"]))
        out.append(text(tip + VALUE_GAP, y + BAR / 2 + CAP * FONT / 2, value_label(value, ax["unit"])))
    return out


def footer(ax: dict, left: float, plot: float, top: float) -> tuple[list[str], int]:
    """The decade tick labels below `top`, at every decade, or at every nth from the first where adjacent
    labels would touch; then the axis title, right-aligned on the tick row when it fits left of the first
    tick label, else centred on the plot one row down. Returns them and the layout's height."""
    decades = range(ax["low"], ax["high"] + 1)
    spacing = plot / (ax["high"] - ax["low"])
    widest = max(est(tick_label(d)) for d in decades)
    step = 1
    while spacing * step < widest + 6:
        step += 1
    tick_baseline = top + LINE + 1
    out = [
        text(decade_x(ax, d, left, plot), tick_baseline, tick_label(d), anchor="middle", opacity=SECONDARY)
        for d in decades[::step]
    ]
    axis_title = f"{'Seconds' if ax['unit'] == 's' else 'Minutes'}, log scale"
    title_end = left - est(tick_label(ax["low"])) / 2 - 8
    if est(axis_title) <= title_end - INSET:
        out.append(text(title_end, tick_baseline, axis_title, anchor="end", opacity=SECONDARY))
        last_baseline = tick_baseline
    else:
        last_baseline = tick_baseline + LINE + 4
        out.append(text(left + plot / 2, last_baseline, axis_title, anchor="middle", opacity=SECONDARY))
    return out, math.ceil(last_baseline + DESCENT * FONT + PAD)


def group_height(chart: dict) -> int:
    return len(chart["series"]) * BAR + (len(chart["series"]) - 1) * GAP


def wide_labels(dataset: dict) -> list[tuple[str, str]]:
    """A dataset's name beside its bars as (line, suffix) rows: the name wrapped, the assembly as a
    secondary-ink suffix on the last line when the name wrapped, else on a line of its own."""
    lines = wrap(short_name(dataset), LABEL_CHARS)
    if len(lines) == 1:
        return [(lines[0], ""), ("", dataset["assembly"])]
    return [(line, "") for line in lines[:-1]] + [(lines[-1], dataset["assembly"])]


def wide_layout(chart: dict, ax: dict) -> tuple[list[str], int]:
    """Each dataset's name right-aligned beside its bars, gridlines across the plot; returns the elements
    and the layout's height. A group takes the height of its bars or of its name block, whichever is
    taller, so blocks never overlap however few bars a group has."""
    out, legend_baseline = heading(chart, WIDE_WIDTH)
    labels = [wide_labels(d) for d in chart["datasets"]]
    left = max(est(f"{line} {suffix}".strip()) for rows in labels for line, suffix in rows) + 10
    plot = plot_width(ax, left, WIDE_WIDTH - PAD)
    group = group_height(chart)
    blocks = [(len(rows) - 1) * LINE + CAP * FONT + DESCENT * FONT for rows in labels]
    slots = [max(group, block) for block in blocks]
    top = legend_baseline + math.ceil(DESCENT * FONT) + 8
    bottom = top + sum(slots) + (len(slots) - 1) * GROUP_GAP
    out += [gridline(decade_x(ax, d, left, plot), top - 4, bottom + 4) for d in range(ax["low"], ax["high"] + 1)]
    slot_top = top
    for dataset, rows, slot in zip(chart["datasets"], labels, slots):
        centre = slot_top + slot / 2
        first_baseline = centre - ((len(rows) - 1) * LINE - CAP * FONT) / 2
        for index, (line, suffix) in enumerate(rows):
            y = first_baseline + index * LINE
            if line:
                out.append(text(left - 10, y, line, anchor="end", suffix=suffix))
            else:
                out.append(text(left - 10, y, suffix, anchor="end", opacity=SECONDARY))
        out += bars(chart, ax, dataset, left, plot, centre - group / 2)
        slot_top += slot + GROUP_GAP
    rest, height = footer(ax, left, plot, bottom + 2)
    return out + rest, height


def narrow_layout(chart: dict, ax: dict) -> tuple[list[str], int]:
    """Each dataset's name, then its assembly in secondary ink, on a line above its bars, so the bars take the
    whole width; gridlines span each group's bars only, so no line crosses a name. Returns the elements and
    the layout's height."""
    out, legend_baseline = heading(chart, NARROW_WIDTH)
    left = est(tick_label(ax["low"])) / 2
    plot = plot_width(ax, left, NARROW_WIDTH - PAD)
    group = group_height(chart)
    max_chars = chars_within(NARROW_WIDTH - INSET - PAD)
    baseline = legend_baseline + 30
    bottom = baseline
    for dataset in chart["datasets"]:
        name = short_name(dataset)
        if len(name) + 1 + len(dataset["assembly"]) <= max_chars:
            lines = [(name, dataset["assembly"])]
        else:
            lines = [(line, "") for line in wrap(name, max_chars)] + [("", dataset["assembly"])]
        for index, (line, assembly) in enumerate(lines):
            y = baseline + index * LINE
            if line:
                out.append(text(INSET, y, line, suffix=assembly))
            else:
                out.append(text(INSET, y, assembly, opacity=SECONDARY))
        group_top = baseline + (len(lines) - 1) * LINE + math.ceil(DESCENT * FONT) + 2
        bottom = group_top + group
        decades = range(ax["low"], ax["high"] + 1)
        out += [gridline(decade_x(ax, d, left, plot), group_top - 2, bottom + 2) for d in decades]
        out += bars(chart, ax, dataset, left, plot, group_top)
        baseline = bottom + math.ceil(CAP * FONT) + 8
    rest, height = footer(ax, left, plot, bottom + 2)
    return out + rest, height


def render(chart: dict) -> str:
    """The SVG as text: the complete root tag alone on the first line, no blank line, a trailing newline."""
    ax = axis(chart)
    narrow, narrow_height = narrow_layout(chart, ax)
    wide, wide_height = wide_layout(chart, ax)
    root, narrow_class, wide_class = f".{CLASS}", f".{CLASS} .{CLASS}-narrow", f".{CLASS} .{CLASS}-wide"
    style = (
        f"{root}{{display:block;max-width:min(100%,{NARROW_MAX_WIDTH}px);height:auto;"
        f"aspect-ratio:{NARROW_WIDTH}/{narrow_height}}}"
        f"{narrow_class},{wide_class}{{width:100%;height:100%}}"
        f"{wide_class}{{display:none}}"
        f"@media (min-width:{BREAKPOINT}px){{{root}{{max-width:100%;aspect-ratio:{WIDE_WIDTH}/{wide_height}}}"
        f"{narrow_class}{{display:none}}{wide_class}{{display:inline}}}}"
    )
    legend_text = ", ".join(s["label"] for s in chart["series"])
    desc = (
        f"Release v{chart['version']}: whole-genome wall time per dataset on {ARCH_LABEL} for {legend_text}, "
        "from the release's population record."
    )
    out = [
        (f'<svg xmlns="http://www.w3.org/2000/svg" class="{CLASS}" width="{WIDE_WIDTH}" height="{wide_height}" '
         f'role="img" aria-label="{html.escape(headline(chart))}" data-release="v{chart["version"]}">'),
        f"<desc>{html.escape(desc, quote=False)}</desc>",
        f"<style>{style}</style>",
        f'<svg class="{CLASS}-narrow" viewBox="0 0 {NARROW_WIDTH} {narrow_height}">',
        *narrow,
        "</svg>",
        f'<svg class="{CLASS}-wide" viewBox="0 0 {WIDE_WIDTH} {wide_height}">',
        *wide,
        "</svg>",
        "</svg>",
    ]
    return "\n".join(out) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description="Render the current release's wall-time chart as an SVG.")
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2], help="repository root")
    ap.add_argument("--out", type=Path, help="output path; default docs/assets/benchmark.svg under --repo")
    args = ap.parse_args()
    version = current_release(args.repo / "README.md")
    path = population_record(args.repo, version)
    chart = read_chart(load_record(path, version), path.name, version)
    svg = render(chart)
    if "\n\n" in svg:
        fail("the rendered SVG contains a blank line")
    out = args.out or args.repo / "docs" / "assets" / "benchmark.svg"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(svg, encoding="utf-8")
    for note in chart["notes"]:
        print(note)
    print(f"wrote {out} ({len(svg.encode('utf-8'))} bytes) for v{version} from {path.name}")


if __name__ == "__main__":
    main()
