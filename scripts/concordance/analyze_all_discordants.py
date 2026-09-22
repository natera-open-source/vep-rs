#!/usr/bin/env python3
"""Aggregate classified discordants across several concordance suites into a pattern map.

Finds every discordant.tsv produced by compare_vep_outputs.py (or
compare_sv_concordance.py), runs classify_discordants.py on each, then builds
a cross-suite pattern map showing which discordance categories appear in which
suites and at what volume.

Input layout, under --output-dir:
  suites/<suite-id>/report/discordant.tsv
      One per suite. run_concordance.sh writes its report as
      <work-dir>/reports/discordant.tsv, so a run whose --work-dir is
      suites/<suite-id> lands in the right place; both the `report/` and
      `reports/` spellings are read.
  dashboard.json (optional)
      {"suites": [{"id": <suite-id>, "name": <display name>,
                   "severity": "blocking" | "major" | "minor"}, ...]}
      Supplies display names and the per-suite severity that sets the
      sample limit. Absent, names are the suite ids and every suite is
      minor.

Outputs, under --output-dir/analysis/ unless --output names another path:
  - pattern_map.json: machine-readable category x suite matrix with samples
  - pattern_map.md: human-readable markdown table

Usage:
  scripts/concordance/analyze_all_discordants.py \
    --output-dir "${VEP_WORK_DIR}/concordance_run" \
    [--repo-dir .] \
    [--output pattern_map.json]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


# Per-severity sample limits for representative discordant rows. A suite's
# severity comes from dashboard.json; a suite without one is minor.
SAMPLE_LIMITS = {"blocking": 10, "major": 5, "minor": 5}
DEFAULT_SAMPLE_LIMIT = 5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate classified discordants across all concordance suites"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Base output directory containing per-suite report directories",
    )
    parser.add_argument(
        "--repo-dir",
        default=None,
        help="Repository root (default: two levels up from this script)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output path for pattern_map.json (default: <output-dir>/analysis/pattern_map.json)",
    )
    return parser.parse_args()


def find_discordant_files(output_dir: Path) -> list[tuple[str, Path]]:
    """Find every discordant.tsv under suites/*/report/ or suites/*/reports/."""
    results = []
    suites_dir = output_dir / "suites"
    if not suites_dir.is_dir():
        return results
    for suite_dir in sorted(suites_dir.iterdir()):
        if not suite_dir.is_dir():
            continue
        for report_dir in ("report", "reports"):
            tsv = suite_dir / report_dir / "discordant.tsv"
            if tsv.exists():
                results.append((suite_dir.name, tsv))
                break
    return results


def load_dashboard(output_dir: Path) -> dict[str, dict]:
    """Load dashboard.json and return a mapping of suite_id -> suite info."""
    dashboard_path = output_dir / "dashboard.json"
    if not dashboard_path.exists():
        return {}
    with open(dashboard_path) as f:
        data = json.load(f)
    return {s["id"]: s for s in data.get("suites", [])}


def run_classify(classify_script: Path, discordant_tsv: Path, output_tsv: Path) -> bool:
    """Run classify_discordants.py on a single discordant.tsv."""
    import subprocess

    cmd = [
        sys.executable,
        str(classify_script),
        str(discordant_tsv),
        "-o",
        str(output_tsv),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(
            f"  WARN: classify_discordants.py failed for {discordant_tsv}: "
            f"{result.stderr.strip()[:200]}",
            file=sys.stderr,
        )
        return False
    return True


def read_classified_rows(path: Path) -> list[dict[str, str]]:
    """Read a classified TSV and return all rows as dicts."""
    if not path.exists():
        return []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def build_pattern_map(
    suite_rows: dict[str, list[dict[str, str]]],
    suite_info: dict[str, dict],
) -> dict:
    """Build the cross-suite pattern map structure."""
    # Collect all categories and per-suite counts
    category_suite_counts: dict[str, dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    category_samples: dict[str, list[dict[str, str]]] = defaultdict(list)

    all_suite_ids = sorted(suite_rows.keys())

    for suite_id in all_suite_ids:
        rows = suite_rows[suite_id]
        severity = suite_info.get(suite_id, {}).get("severity", "minor")
        sample_limit = SAMPLE_LIMITS.get(severity, DEFAULT_SAMPLE_LIMIT)

        per_category_count: dict[str, int] = defaultdict(int)
        per_category_sampled: dict[str, int] = defaultdict(int)

        for row in rows:
            category = row.get("mismatch_type", "unclassified")
            per_category_count[category] += 1
            category_suite_counts[category][suite_id] += 1

            if per_category_sampled[category] < sample_limit:
                sample = {
                    "suite": suite_id,
                    "location": row.get("location", ""),
                    "allele": row.get("allele", ""),
                    "feature": row.get("feature", ""),
                    "consequence_set": row.get(
                        "consequence_set", row.get("consequence", "")
                    ),
                    "source": row.get("source", row.get("direction", "")),
                }
                category_samples[category].append(sample)
                per_category_sampled[category] += 1

    # Build sorted category list
    categories = []
    for cat_name in sorted(
        category_suite_counts.keys(),
        key=lambda c: sum(category_suite_counts[c].values()),
        reverse=True,
    ):
        suites_map = dict(category_suite_counts[cat_name])
        total = sum(suites_map.values())
        categories.append(
            {
                "name": cat_name,
                "suites": suites_map,
                "total": total,
                "samples": category_samples[cat_name],
            }
        )

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "suite_ids": all_suite_ids,
        "suite_names": {
            sid: suite_info.get(sid, {}).get("name", sid) for sid in all_suite_ids
        },
        "categories": categories,
    }


def write_pattern_map_markdown(
    path: Path,
    pattern_map: dict,
) -> None:
    """Write a markdown table summarizing the pattern map."""
    suite_ids = pattern_map["suite_ids"]
    suite_names = pattern_map["suite_names"]
    categories = pattern_map["categories"]

    # Column headers abbreviate the assembly in each suite's display name.
    short_names = {}
    for sid in suite_ids:
        name = suite_names.get(sid, sid)
        short = name.replace("GRCh37", "37").replace("GRCh38", "38")
        short_names[sid] = short.strip()

    lines = [
        "# Discordant Pattern Map",
        "",
        f"Generated: {pattern_map['generated_at']}",
        "",
    ]

    # Header row
    header = "| Category |"
    sep = "| --- |"
    for sid in suite_ids:
        header += f" {short_names[sid]} |"
        sep += " ---: |"
    header += " Total |"
    sep += " ---: |"
    lines.append(header)
    lines.append(sep)

    # Data rows
    for cat in categories:
        row = f"| {cat['name']} |"
        for sid in suite_ids:
            count = cat["suites"].get(sid, 0)
            row += f" {count:,} |" if count else " |"
        row += f" {cat['total']:,} |"
        lines.append(row)

    # Totals row
    totals_row = "| **Total** |"
    grand_total = 0
    for sid in suite_ids:
        col_total = sum(cat["suites"].get(sid, 0) for cat in categories)
        totals_row += f" **{col_total:,}** |" if col_total else " |"
        grand_total += col_total
    totals_row += f" **{grand_total:,}** |"
    lines.append(totals_row)

    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()

    if args.repo_dir:
        repo_dir = Path(args.repo_dir).resolve()
    else:
        repo_dir = Path(__file__).resolve().parent.parent.parent

    classify_script = repo_dir / "scripts" / "validation" / "classify_discordants.py"
    if not classify_script.exists():
        print(
            f"ERROR: classify_discordants.py not found at {classify_script}",
            file=sys.stderr,
        )
        return 1

    # Find all discordant files
    disc_files = find_discordant_files(output_dir)
    if not disc_files:
        print(
            f"WARN: No discordant.tsv files found under {output_dir}/suites/*/report(s)/",
            file=sys.stderr,
        )
        return 0

    # Load dashboard for suite metadata
    suite_info = load_dashboard(output_dir)

    print(f"Found {len(disc_files)} discordant files")

    # Classify each suite's discordants
    suite_rows: dict[str, list[dict[str, str]]] = {}
    for suite_id, disc_tsv in disc_files:
        classified_tsv = disc_tsv.parent / "discordant_classified.tsv"
        print(f"  Classifying {suite_id}...")
        if run_classify(classify_script, disc_tsv, classified_tsv):
            rows = read_classified_rows(classified_tsv)
            suite_rows[suite_id] = rows
            print(f"    {len(rows)} discordants classified")
        else:
            print(f"    WARN: skipping {suite_id} (classification failed)")

    if not suite_rows:
        print("WARN: No suites were successfully classified", file=sys.stderr)
        return 0

    # Build pattern map
    pattern_map = build_pattern_map(suite_rows, suite_info)

    # Write outputs
    analysis_dir = output_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    if args.output:
        json_path = Path(args.output).resolve()
    else:
        json_path = analysis_dir / "pattern_map.json"

    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(pattern_map, indent=2), encoding="utf-8")

    md_path = json_path.with_suffix(".md")
    write_pattern_map_markdown(md_path, pattern_map)

    # Summary
    total_discordants = sum(cat["total"] for cat in pattern_map["categories"])
    print(
        f"\nPattern map: {len(pattern_map['categories'])} categories, "
        f"{total_discordants:,} total discordants across {len(suite_rows)} suites"
    )
    print(f"JSON: {json_path}")
    print(f"Markdown: {md_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
