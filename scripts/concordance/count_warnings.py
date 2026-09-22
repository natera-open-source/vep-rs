#!/usr/bin/env python3
"""Aggregate per-suite warning counts from engine stderr logs.

Reads ``<work_dir>/{perl,rust}/<base>.stderr.log`` (one file per benchmark
suite per engine), classifies each non-blank line as ERROR / WARN / INFO /
OTHER via simple regex, and emits a single JSON summary at
``<work_dir>/reports/warning_counts.json`` plus a Markdown section written to
``--out-md``.

Captures verbose Perl + Rust output, classifies it, and surfaces counts so
a reader can scan for run health without trawling raw logs.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Iterable

# Order matters: ERROR > WARN > INFO. First match wins.
PATTERNS = [
    (
        "ERROR",
        re.compile(r"\b(?:ERROR|FATAL|panic(?:ked)?|Traceback)\b", re.IGNORECASE),
    ),
    ("WARN", re.compile(r"\b(?:WARN(?:ING)?|deprecated|skipping)\b", re.IGNORECASE)),
    # INFO is anchored on log-prefix patterns only; mid-sentence keywords like
    # "loaded"/"wrote"/"reading" fall to OTHER, which is not surveilled.
    ("INFO", re.compile(r"^(?:\[INFO\]|INFO[:\s])", re.IGNORECASE)),
]


def classify(line: str) -> str:
    for label, pat in PATTERNS:
        if pat.search(line):
            return label
    return "OTHER"


def count_file(path: Path) -> dict[str, int]:
    counts = {"ERROR": 0, "WARN": 0, "INFO": 0, "OTHER": 0, "lines": 0}
    if not path.exists():
        return counts
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        counts["lines"] += 1
        counts[classify(line)] += 1
    return counts


def discover_logs(work_dir: Path) -> dict[str, dict[str, Path]]:
    """Return ``{base_name: {engine: stderr_path}}`` for every (engine, base) pair found.

    Supports two layouts:

    1. ``run_concordance.sh`` ad-hoc layout: ``<work_dir>/{perl,rust}/<base>.stderr.log``.
    2. Per-suite layout, as ``run_concordance.sh --work-dir <work_dir>/suites/<sid>``
       writes it: ``<work_dir>/suites/<sid>/{perl,rust}/<base>.stderr*`` (covers both
       the single-file SNP case and the multi-VCF SV case).
    """
    found: dict[str, dict[str, Path]] = {}

    # Layout 1: flat <work_dir>/{perl,rust}/<base>.stderr.log
    for engine in ("perl", "rust"):
        engine_dir = work_dir / engine
        if not engine_dir.is_dir():
            continue
        for path in sorted(engine_dir.glob("*.stderr.log")):
            base = path.name[: -len(".stderr.log")]
            found.setdefault(base, {})[engine] = path

    # Layout 2: <work_dir>/suites/<sid>/{perl,rust}/*.stderr{,.log}
    suites_dir = work_dir / "suites"
    if suites_dir.is_dir():
        for sid_dir in sorted(p for p in suites_dir.iterdir() if p.is_dir()):
            sid = sid_dir.name
            for engine in ("perl", "rust"):
                engine_dir = sid_dir / engine
                if not engine_dir.is_dir():
                    continue
                # Both .stderr.log (SNP and indel suites) and .stderr (SV-loop) variants.
                for pattern in ("*.stderr.log", "*.stderr"):
                    for path in sorted(engine_dir.glob(pattern)):
                        base_local = path.name
                        for suffix in (".stderr.log", ".stderr"):
                            if base_local.endswith(suffix):
                                base_local = base_local[: -len(suffix)]
                                break
                        # Disambiguate per suite when multiple SV files share a base.
                        full_key = f"{sid}/{base_local}"
                        found.setdefault(full_key, {})[engine] = path
    return found


def render_markdown(report: dict) -> str:
    lines = ["## Warning counts (from engine stderr)\n"]
    suites = report.get("suites", {})
    if not suites:
        lines.append("No engine stderr logs found.\n")
        return "\n".join(lines)
    lines.append("| Suite | Engine | ERROR | WARN | INFO | Other | Lines |")
    lines.append("| ----- | ------ | -----:| ----:| ----:| -----:| -----:|")
    for base in sorted(suites):
        for engine in ("perl", "rust"):
            c = suites[base].get(engine)
            if not c:
                continue
            lines.append(
                f"| {base} | {engine} | {c['ERROR']} | {c['WARN']} | {c['INFO']} | "
                f"{c['OTHER']} | {c['lines']} |"
            )
    totals = report.get("totals", {})
    lines.append("")
    lines.append(
        f"**Totals**: ERROR={totals.get('ERROR', 0)}, WARN={totals.get('WARN', 0)}, "
        f"INFO={totals.get('INFO', 0)}, Other={totals.get('OTHER', 0)}, "
        f"Lines={totals.get('lines', 0)}"
    )
    return "\n".join(lines) + "\n"


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="Concordance work directory; expects perl/ and rust/ subdirs with *.stderr.log",
    )
    parser.add_argument(
        "--out-json",
        type=Path,
        default=None,
        help="Output JSON path. Defaults to <work_dir>/reports/warning_counts.json.",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=None,
        help="Optional Markdown output path. When set, also writes a renderable section.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    work_dir = args.work_dir.resolve()
    out_json = args.out_json or (work_dir / "reports" / "warning_counts.json")
    out_json.parent.mkdir(parents=True, exist_ok=True)

    pairs = discover_logs(work_dir)
    suites: dict[str, dict[str, dict[str, int]]] = {}
    totals = {"ERROR": 0, "WARN": 0, "INFO": 0, "OTHER": 0, "lines": 0}

    for base, engines in pairs.items():
        suites[base] = {}
        for engine, path in engines.items():
            counts = count_file(path)
            suites[base][engine] = counts
            for key, val in counts.items():
                totals[key] = totals.get(key, 0) + val

    report = {
        "work_dir": str(work_dir),
        "suite_count": len(suites),
        "totals": totals,
        "suites": suites,
    }
    out_json.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    if args.out_md is not None:
        args.out_md.parent.mkdir(parents=True, exist_ok=True)
        args.out_md.write_text(render_markdown(report), encoding="utf-8")

    print(
        f"warning_counts: {len(suites)} suites, "
        f"ERROR={totals['ERROR']} WARN={totals['WARN']} "
        f"INFO={totals['INFO']} OTHER={totals['OTHER']} "
        f"-> {out_json}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
