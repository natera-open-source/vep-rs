#!/usr/bin/env python3
"""Section-by-section comparison of two VEP statistics reports.

``--expected F --actual F [--tolerance N] [--mask-row LABEL ...] [--report F]``
    One Ensembl VEP statistics report (the ``<output>_summary.txt`` that ``--stats_text``
    writes) against a candidate's, a JSON report with a ``verdict``, exit 0 on ``green``
    and 1 on ``red``. Either file may be gzipped (a ``.gz`` name).

The grammar read is ``Stats.pm:806-831`` at release/116.2: a ``[Section]`` line opens each
section, a row is ``key<TAB>value``, sections are separated by one blank line, and the
HTML tags of the three tabular sections are stripped before they are written. The three
tabular sections come first and in this order: ``[VEP run statistics]`` (``Stats.pm:624-681``),
``[Data version]`` (``:695-700``; the one line ``No data version information available.``
when the cache names no data versions) and ``[General statistics]`` (``:716-740``, a row whose
value is undefined is written as its key alone). Every later section is a chart
(``:505-609``): ``key<TAB>count`` rows in the chart's sort order.

The rules:

* A tabular section compares row by row: the same keys in the same order, every value
  verbatim, except the masked rows, whose values no second run reproduces and which
  compare by key alone: ``Command line options``, ``Start time``, ``End time``, ``Run time``,
  ``Input file``, ``Input data`` and ``Output file``; in ``Annotation sources`` every path
  (a token carrying ``/``) is replaced by ``[PATH]`` before the values compare, so
  ``Cache: /a/cache`` and ``Cache: [PATH]/116_GRCh38`` agree. ``--mask-row LABEL``
  (repeatable) masks the value of one more row wherever its label occurs.
* A chart section compares as a set of ``(key, count)`` pairs, since ``sort_keys``
  (``Stats.pm:760-790``) breaks ties in hash order: two keys with one count, or two
  non-numeric chromosome names sharing a first character, come out in an order that
  differs from run to run. Under ``--tolerance N`` (N at least 1, for a comparison over
  inputs on which the documented divergences apply) the two counts of a key
  may differ by up to N, a key absent on one side counting as 0 there; at the default
  tolerance 0 the two sets must be equal.
* The two files must carry the same sections in the same order: a chart Ensembl VEP writes
  only under a flag (``[SIFT summary]``), only with a coding row (``[Position in protein]``)
  or only with a FASTA (the per-chromosome distributions) is expected of the candidate
  under the same run.

The report: ``mode`` and ``format`` (both ``stats``), ``verdict``, ``reasons`` (one line per
section that differs), ``expected`` and ``actual`` (the paths), ``tolerance``,
``masked_rows``, ``sections`` (each file's section titles in order and whether the shared
ones agree in order), ``compared`` (sections and rows compared) and ``differences``, a map
from section title to the list of its differences (``kind`` one of ``value``, ``count``,
``missing``, ``unexpected``, ``order``, ``duplicate``, ``section_missing``,
``section_unexpected``; ``row``, ``expected`` and ``actual`` where they apply).
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import zlib
from pathlib import Path

COMPONENT = "compare_vep_stats"
MODE = "stats"

RUN_STATISTICS = "VEP run statistics"
DATA_VERSION = "Data version"
GENERAL_STATISTICS = "General statistics"
TABULAR_SECTIONS = (RUN_STATISTICS, DATA_VERSION, GENERAL_STATISTICS)

# The rows of [VEP run statistics] whose values differ between any two runs
# (Stats.pm:645-663): the command line names the engine's own options and paths, the
# times and the run time are the run's own, the input and output rows are paths.
MASKED_ROWS = (
    "Command line options",
    "Start time",
    "End time",
    "Run time",
    "Input file",
    "Input data",
    "Output file",
)
# Stats.pm:666-677: `Cache: <dir>`, `<db> on <host>` and `Custom: <file> (<type>)` joined
# by `; `; the directory and the files are paths.
PATH_ROW = "Annotation sources"
MASKED = "<masked>"
PATH = "[PATH]"
# A token carrying a slash, up to a blank or the `; ` that joins two sources.
PATH_TOKEN_RE = re.compile(r"[^\s;]*/[^\s;]*")
SECTION_RE = re.compile(r"^\[(.*)\]$")


class ComparatorError(Exception):
    """A file this comparator cannot read or an option it cannot apply."""


def read_text(path: Path) -> str:
    """The text of ``path``, gunzipped when its name ends in ``.gz``."""
    try:
        if str(path).endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8") as fh:
                return fh.read()
        return Path(path).read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError, EOFError, zlib.error) as exc:
        raise ComparatorError(f"cannot read {path}: {exc}") from exc


class Section:
    """One ``[title]`` block: its rows as ``(key, value)`` pairs in file order, a row
    without a tab carrying ``None`` as its value."""

    def __init__(self, title: str) -> None:
        self.title = title
        self.rows: list[tuple[str, str | None]] = []

    @property
    def is_chart(self) -> bool:
        return self.title not in TABULAR_SECTIONS

    def ordered_keys(self) -> list[str]:
        return [key for key, _ in self.rows]


def parse_report(text: str) -> list[Section]:
    """The sections of a statistics report in file order. A non-blank line before the
    first section header is refused: the file is not a statistics report."""
    sections: list[Section] = []
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r")
        if not line.strip():
            continue
        header = SECTION_RE.match(line)
        if header:
            sections.append(Section(header.group(1)))
            continue
        if not sections:
            raise ComparatorError(f"line {number} precedes the first [section] header: {line!r}")
        key, sep, value = line.partition("\t")
        sections[-1].rows.append((key, value if sep else None))
    return sections


def mask_value(key: str, value: str | None, masked_rows: frozenset[str]) -> str | None:
    """The value as it compares: ``<masked>`` for a masked row, paths replaced in
    ``Annotation sources``, verbatim otherwise."""
    if key in masked_rows:
        return MASKED
    if key == PATH_ROW and value is not None:
        return PATH_TOKEN_RE.sub(PATH, value)
    return value


def rows_by_key(section: Section, side: str, differences: list[dict]) -> dict[str, str | None]:
    """The section's rows as a map, every key printed twice recorded as a difference."""
    out: dict[str, str | None] = {}
    for key, value in section.rows:
        if key in out:
            differences.append({"kind": "duplicate", "row": key, "side": side})
            continue
        out[key] = value
    return out


def compare_tabular(expected: Section, actual: Section, masked_rows: frozenset[str]) -> list[dict]:
    """Row by row: the same keys in the same order, the masked values set aside."""
    differences: list[dict] = []
    e_rows = rows_by_key(expected, "expected", differences)
    a_rows = rows_by_key(actual, "actual", differences)
    for key in e_rows:
        if key not in a_rows:
            differences.append({"kind": "missing", "row": key, "expected": e_rows[key]})
    for key in a_rows:
        if key not in e_rows:
            differences.append({"kind": "unexpected", "row": key, "actual": a_rows[key]})
    shared_e = [k for k in e_rows if k in a_rows]
    shared_a = [k for k in a_rows if k in e_rows]
    if shared_e != shared_a:
        differences.append({"kind": "order", "expected": shared_e, "actual": shared_a})
    for key in shared_e:
        e_value = mask_value(key, e_rows[key], masked_rows)
        a_value = mask_value(key, a_rows[key], masked_rows)
        if e_value != a_value:
            differences.append({"kind": "value", "row": key, "expected": e_value, "actual": a_value})
    return differences


def count_of(value: str | None) -> int | str | None:
    """A chart row's count as an integer, or the text as written when it is not one."""
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return value


def compare_chart(
    expected: Section, actual: Section, tolerance: int, masked_rows: frozenset[str]
) -> list[dict]:
    """As a set of ``(key, count)`` pairs; under a tolerance of at least 1 a key's two
    counts may differ by up to the tolerance, an absent key counting as 0."""
    differences: list[dict] = []
    e_rows = rows_by_key(expected, "expected", differences)
    a_rows = rows_by_key(actual, "actual", differences)
    keys = list(e_rows) + [k for k in a_rows if k not in e_rows]
    for key in keys:
        if key in masked_rows:
            continue
        e_count = count_of(e_rows.get(key))
        a_count = count_of(a_rows.get(key))
        if key in e_rows and key in a_rows:
            kind = "count"
            within = e_count == a_count or (
                isinstance(e_count, int) and isinstance(a_count, int) and abs(e_count - a_count) <= tolerance
            )
        else:
            kind = "missing" if key in e_rows else "unexpected"
            present = e_count if key in e_rows else a_count
            within = tolerance >= 1 and isinstance(present, int) and present <= tolerance
        if not within:
            differences.append({"kind": kind, "row": key, "expected": e_count, "actual": a_count})
    return differences


def compare_reports(
    expected_text: str, actual_text: str, tolerance: int = 0, mask_rows: tuple[str, ...] = ()
) -> dict:
    """The report of the two statistics texts (see the module docstring)."""
    if tolerance < 0:
        raise ComparatorError(f"--tolerance {tolerance} is negative")
    masked_rows = frozenset(MASKED_ROWS) | frozenset(mask_rows)
    expected = parse_report(expected_text)
    actual = parse_report(actual_text)
    e_titles = [s.title for s in expected]
    a_titles = [s.title for s in actual]
    differences: dict[str, list[dict]] = {}
    reasons: list[str] = []
    for side, titles in (("expected", e_titles), ("actual", a_titles)):
        if not titles:
            reasons.append(f"the {side} file carries no [section] header")
        for title in titles:
            repeated = f"section [{title}] appears twice in the {side} file"
            if titles.count(title) > 1 and repeated not in reasons:
                reasons.append(repeated)
    for title in e_titles:
        if title not in a_titles:
            differences[title] = [{"kind": "section_missing"}]
            reasons.append(f"section [{title}] only in the expected file")
    for title in a_titles:
        if title not in e_titles:
            differences[title] = [{"kind": "section_unexpected"}]
            reasons.append(f"section [{title}] only in the actual file")
    shared_e = [t for t in dict.fromkeys(e_titles) if t in a_titles]
    shared_a = [t for t in dict.fromkeys(a_titles) if t in e_titles]
    order_agrees = shared_e == shared_a
    if not order_agrees:
        reasons.append("the shared sections come in another order in the actual file")
    # A repeated title is reported above; on either side its first copy is the one compared.
    actual_by_title: dict[str, Section] = {}
    for section in actual:
        actual_by_title.setdefault(section.title, section)
    rows_compared = 0
    compared_titles: set[str] = set()
    for section in expected:
        other = actual_by_title.get(section.title)
        if other is None or section.title in compared_titles:
            continue
        compared_titles.add(section.title)
        if section.is_chart:
            found = compare_chart(section, other, tolerance, masked_rows)
        else:
            found = compare_tabular(section, other, masked_rows)
        rows_compared += len(section.rows)
        if found:
            differences[section.title] = found
            noun = "difference" if len(found) == 1 else "differences"
            reasons.append(f"[{section.title}]: {len(found)} {noun}")
    return {
        "mode": MODE,
        "format": MODE,
        "verdict": "red" if reasons else "green",
        "reasons": reasons,
        "tolerance": tolerance,
        "masked_rows": sorted(masked_rows),
        "sections": {"expected": e_titles, "actual": a_titles, "order_agrees": order_agrees},
        "compared": {"sections": len(shared_e), "rows": rows_compared},
        "differences": differences,
    }


def compare_files(expected: Path, actual: Path, tolerance: int = 0, mask_rows: tuple[str, ...] = ()) -> dict:
    report = compare_reports(read_text(expected), read_text(actual), tolerance, mask_rows)
    report["expected"] = str(expected)
    report["actual"] = str(actual)
    return report


def write_report(report: dict, path: Path | None) -> None:
    if path is None:
        return
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    except OSError as exc:
        raise ComparatorError(f"cannot write {path}: {exc.strerror or exc}") from exc


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="compare_vep_stats.py",
        description="Section-by-section comparison of two VEP statistics reports; see the module docstring for the rules.",
    )
    ap.add_argument("--expected", type=Path, required=True, help="the Ensembl VEP statistics report")
    ap.add_argument("--actual", type=Path, required=True, help="the candidate's statistics report")
    ap.add_argument(
        "--tolerance",
        type=int,
        default=0,
        help="the count difference a chart row may show (default 0: the sets must be equal)",
    )
    ap.add_argument(
        "--mask-row",
        action="append",
        default=[],
        metavar="LABEL",
        help="a further row label whose value is set aside (repeatable)",
    )
    ap.add_argument("--report", type=Path, help="where the JSON report is written")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = compare_files(args.expected, args.actual, args.tolerance, tuple(args.mask_row))
        write_report(report, args.report)
    except ComparatorError as exc:
        print(f"ERROR: [{COMPONENT}] {exc}", file=sys.stderr)
        return 2
    print(f"{MODE}: {report['verdict']}")
    for reason in report["reasons"]:
        print(f"  {reason}")
    return 0 if report["verdict"] == "green" else 1


if __name__ == "__main__":
    sys.exit(main())
