#!/usr/bin/env python3
"""Structurally compare two JSON cache directories (baseline vs builder).

Compares transcripts, variations, and info.json field-by-field, producing
a JSON report with match rates and detailed mismatch breakdowns.

A region shard is ``<cache>/<transcripts|variations>/<chr>/<start>-<end>.json``
or the same name with a ``.gz`` suffix; the two sides are matched on the name
with any ``.gz`` removed, so a gzipped cache (``tests/golden/*/*/json_cache``)
compares against a plain one. A shard is a JSON array of records. A JSON object
is accepted too: its values are either records keyed by id or arrays of records
keyed by chromosome.

Usage:
    python3 scripts/validation/diff_caches.py \
      --baseline-dir /path/to/perl_cache \
      --builder-dir /path/to/builder_cache \
      --output-report /path/to/diff_report.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path


# Transcript comparison

TRANSCRIPT_COMPARE_FIELDS = [
    "start",
    "end",
    "strand",
    "biotype",
    "source",
    "gene_stable_id",
    "gene_symbol",
    "hgnc_id",
    "canonical",
]

TRANSCRIPT_VEFC_FIELDS = [
    "introns",
    "sorted_exons",
    "five_prime_utr",
    "three_prime_utr",
    "codon_table",
    "protein_features",
    "seq_edits",
]


def _count_list_field(transcript: dict, field: str) -> int | None:
    """Return the length of a list/dict field, or None if absent."""
    val = transcript.get(field)
    if val is None:
        # Check inside variation_effect_feature_cache / vefc
        vefc = transcript.get("variation_effect_feature_cache") or transcript.get(
            "vefc", {}
        )
        val = vefc.get(field)
    if isinstance(val, (list, dict)):
        return len(val)
    return None


def compare_transcripts(baseline: dict, builder: dict) -> dict:
    """Compare two transcript dicts. Returns dict of mismatched fields."""
    mismatches = {}
    for field in TRANSCRIPT_COMPARE_FIELDS:
        bval = baseline.get(field)
        cval = builder.get(field)
        if bval != cval:
            mismatches[field] = {"baseline": bval, "builder": cval}

    # Compare exon/intron/mapper counts
    b_vefc = baseline.get("variation_effect_feature_cache") or baseline.get("vefc", {})
    c_vefc = builder.get("variation_effect_feature_cache") or builder.get("vefc", {})

    b_exons = b_vefc.get("sorted_exons") or baseline.get("sorted_exons")
    c_exons = c_vefc.get("sorted_exons") or builder.get("sorted_exons")
    b_exon_count = len(b_exons) if isinstance(b_exons, list) else 0
    c_exon_count = len(c_exons) if isinstance(c_exons, list) else 0
    if b_exon_count != c_exon_count:
        mismatches["exon_count"] = {"baseline": b_exon_count, "builder": c_exon_count}

    b_introns = b_vefc.get("introns") or baseline.get("introns")
    c_introns = c_vefc.get("introns") or builder.get("introns")
    b_intron_count = len(b_introns) if isinstance(b_introns, list) else 0
    c_intron_count = len(c_introns) if isinstance(c_introns, list) else 0
    if b_intron_count != c_intron_count:
        mismatches["intron_count"] = {
            "baseline": b_intron_count,
            "builder": c_intron_count,
        }

    b_mapper = b_vefc.get("mapper") or baseline.get("mapper")
    c_mapper = c_vefc.get("mapper") or builder.get("mapper")
    b_mapper_count = len(b_mapper) if isinstance(b_mapper, (list, dict)) else 0
    c_mapper_count = len(c_mapper) if isinstance(c_mapper, (list, dict)) else 0
    if b_mapper_count != c_mapper_count:
        mismatches["mapper_pair_count"] = {
            "baseline": b_mapper_count,
            "builder": c_mapper_count,
        }

    # Compare VEFC-level fields
    for field in TRANSCRIPT_VEFC_FIELDS:
        bval = b_vefc.get(field)
        cval = c_vefc.get(field)
        if bval != cval:
            # For large fields, just note length difference
            if isinstance(bval, (list, dict)) and isinstance(cval, (list, dict)):
                if len(bval) != len(cval):
                    mismatches[f"vefc_{field}_count"] = {
                        "baseline": len(bval),
                        "builder": len(cval),
                    }
                else:
                    mismatches[f"vefc_{field}"] = "content_differs"
            else:
                mismatches[f"vefc_{field}"] = {"baseline": bval, "builder": cval}

    return mismatches


def _read_json(path: Path):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def _shard_records(data, wrapper_key: str) -> list[dict]:
    """Flatten a shard into its record dicts, whatever container it uses.

    Accepts a bare array, an object wrapping the array under ``wrapper_key``, an
    object of arrays keyed by chromosome, or an object of records keyed by id.
    Null entries and non-dict values are dropped.
    """
    if isinstance(data, dict):
        data = data.get(wrapper_key, data)
    if isinstance(data, list):
        return [r for r in data if isinstance(r, dict)]
    if isinstance(data, dict):
        records: list[dict] = []
        for value in data.values():
            if isinstance(value, list):
                records.extend(r for r in value if isinstance(r, dict))
            elif isinstance(value, dict):
                records.append(value)
        return records
    return []


def load_region_transcripts(path: Path) -> dict[str, dict]:
    """Load a region JSON file and return {stable_id: transcript_dict}."""
    if not path.exists():
        return {}
    result = {}
    for tr in _shard_records(_read_json(path), "transcripts"):
        sid = tr.get("stable_id") or tr.get("id", "")
        if sid:
            result[sid] = tr
    return result


# Variation comparison

VARIATION_COMPARE_FIELDS = [
    "start",
    "end",
    "allele_string",
    "strand",
    "clin_sig",
    "minor_allele",
    "minor_allele_freq",
]


def compare_variations(baseline: dict, builder: dict) -> dict:
    """Compare two variation dicts. Returns dict of mismatched fields."""
    mismatches = {}
    for field in VARIATION_COMPARE_FIELDS:
        bval = baseline.get(field)
        cval = builder.get(field)
        if bval != cval:
            mismatches[field] = {"baseline": bval, "builder": cval}

    # Compare frequencies (nested dict)
    b_freqs = baseline.get("frequencies", {})
    c_freqs = builder.get("frequencies", {})
    if b_freqs != c_freqs:
        mismatches["frequencies"] = "content_differs"

    return mismatches


def load_region_variations(path: Path) -> dict[str, dict]:
    """Load a region variation JSON file and return {variation_name: variation_dict}."""
    if not path.exists():
        return {}
    result = {}
    for var in _shard_records(_read_json(path), "variations"):
        name = var.get("variation_name") or var.get("name", "")
        if name:
            result[name] = var
    return result


# Directory walking


def find_region_files(cache_dir: Path, subdir: str) -> list[Path]:
    """Find all region JSON files (plain or gzipped) under cache_dir/{subdir}/{chr}/."""
    base = cache_dir / subdir
    if not base.is_dir():
        return []
    return sorted(list(base.rglob("*.json")) + list(base.rglob("*.json.gz")))


def relative_key(path: Path, cache_dir: Path, subdir: str) -> str:
    """Return the path relative to {cache_dir}/{subdir}/, without any .gz suffix."""
    base = cache_dir / subdir
    key = str(path.relative_to(base))
    return key[: -len(".gz")] if key.endswith(".gz") else key


# Main comparison


def diff_transcripts(baseline_dir: Path, builder_dir: Path) -> dict:
    """Compare all transcript region files between two cache dirs."""
    baseline_files = find_region_files(baseline_dir, "transcripts")
    builder_files = find_region_files(builder_dir, "transcripts")

    baseline_keys = {relative_key(f, baseline_dir, "transcripts"): f for f in baseline_files}
    builder_keys = {relative_key(f, builder_dir, "transcripts"): f for f in builder_files}

    all_keys = sorted(set(baseline_keys) | set(builder_keys))

    total_baseline = 0
    total_builder = 0
    total_shared = 0
    total_only_baseline = 0
    total_only_builder = 0
    total_matched = 0
    total_differed = 0
    field_mismatch_counts: Counter[str] = Counter()
    region_summaries: list[dict] = []

    for key in all_keys:
        b_path = baseline_keys.get(key)
        c_path = builder_keys.get(key)

        b_transcripts = load_region_transcripts(b_path) if b_path else {}
        c_transcripts = load_region_transcripts(c_path) if c_path else {}

        b_ids = set(b_transcripts)
        c_ids = set(c_transcripts)
        shared = b_ids & c_ids
        only_b = b_ids - c_ids
        only_c = c_ids - b_ids

        total_baseline += len(b_ids)
        total_builder += len(c_ids)
        total_shared += len(shared)
        total_only_baseline += len(only_b)
        total_only_builder += len(only_c)

        region_matched = 0
        region_differed = 0
        for sid in sorted(shared):
            diffs = compare_transcripts(b_transcripts[sid], c_transcripts[sid])
            if diffs:
                region_differed += 1
                for field in diffs:
                    field_mismatch_counts[field] += 1
            else:
                region_matched += 1

        total_matched += region_matched
        total_differed += region_differed

        if only_b or only_c or region_differed > 0:
            region_summaries.append(
                {
                    "region": key,
                    "baseline_count": len(b_ids),
                    "builder_count": len(c_ids),
                    "shared": len(shared),
                    "only_baseline": sorted(only_b)[:10],
                    "only_builder": sorted(only_c)[:10],
                    "matched": region_matched,
                    "differed": region_differed,
                }
            )

    denominator = total_shared if total_shared > 0 else 1
    match_rate = total_matched / denominator

    return {
        "total_baseline": total_baseline,
        "total_builder": total_builder,
        "shared": total_shared,
        "only_baseline": total_only_baseline,
        "only_builder": total_only_builder,
        "matched": total_matched,
        "differed": total_differed,
        "field_mismatches": dict(field_mismatch_counts.most_common()),
        "match_rate": round(match_rate, 6),
        "regions_with_diffs": region_summaries[:50],
    }


def diff_variations(baseline_dir: Path, builder_dir: Path) -> dict:
    """Compare all variation region files between two cache dirs."""
    baseline_files = find_region_files(baseline_dir, "variations")
    builder_files = find_region_files(builder_dir, "variations")

    baseline_keys = {relative_key(f, baseline_dir, "variations"): f for f in baseline_files}
    builder_keys = {relative_key(f, builder_dir, "variations"): f for f in builder_files}

    all_keys = sorted(set(baseline_keys) | set(builder_keys))

    total_baseline = 0
    total_builder = 0
    total_shared = 0
    total_only_baseline = 0
    total_only_builder = 0
    total_matched = 0
    total_differed = 0
    field_mismatch_counts: Counter[str] = Counter()

    for key in all_keys:
        b_path = baseline_keys.get(key)
        c_path = builder_keys.get(key)

        b_vars = load_region_variations(b_path) if b_path else {}
        c_vars = load_region_variations(c_path) if c_path else {}

        b_names = set(b_vars)
        c_names = set(c_vars)
        shared = b_names & c_names
        only_b = b_names - c_names
        only_c = c_names - b_names

        total_baseline += len(b_names)
        total_builder += len(c_names)
        total_shared += len(shared)
        total_only_baseline += len(only_b)
        total_only_builder += len(only_c)

        for name in shared:
            diffs = compare_variations(b_vars[name], c_vars[name])
            if diffs:
                total_differed += 1
                for field in diffs:
                    field_mismatch_counts[field] += 1
            else:
                total_matched += 1

    denominator = total_shared if total_shared > 0 else 1
    match_rate = total_matched / denominator

    return {
        "total_baseline": total_baseline,
        "total_builder": total_builder,
        "shared": total_shared,
        "only_baseline": total_only_baseline,
        "only_builder": total_only_builder,
        "matched": total_matched,
        "differed": total_differed,
        "field_mismatches": dict(field_mismatch_counts.most_common()),
        "match_rate": round(match_rate, 6),
    }


def diff_info_json(baseline_dir: Path, builder_dir: Path) -> dict:
    """Compare info.json files between two cache dirs."""
    b_path = baseline_dir / "info.json"
    c_path = builder_dir / "info.json"

    if not b_path.exists() and not c_path.exists():
        return {"matches": True, "diffs": [], "note": "neither cache has info.json"}

    if not b_path.exists():
        return {"matches": False, "diffs": ["info.json missing in baseline"]}
    if not c_path.exists():
        return {"matches": False, "diffs": ["info.json missing in builder"]}

    with open(b_path) as fh:
        b_info = json.load(fh)
    with open(c_path) as fh:
        c_info = json.load(fh)

    if b_info == c_info:
        return {"matches": True, "diffs": []}

    diffs = []
    all_keys = sorted(set(b_info) | set(c_info))
    for key in all_keys:
        bval = b_info.get(key)
        cval = c_info.get(key)
        if bval != cval:
            diffs.append(
                {"field": key, "baseline": str(bval)[:200], "builder": str(cval)[:200]}
            )

    return {"matches": False, "diffs": diffs}


# CLI


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Structurally compare two JSON cache directories."
    )
    parser.add_argument(
        "--baseline-dir",
        required=True,
        help="Path to the baseline (Perl-derived) JSON cache directory",
    )
    parser.add_argument(
        "--builder-dir",
        required=True,
        help="Path to the builder-generated JSON cache directory",
    )
    parser.add_argument(
        "--output-report",
        required=True,
        help="Path for the output JSON diff report",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    baseline_dir = Path(args.baseline_dir).resolve()
    builder_dir = Path(args.builder_dir).resolve()
    output_report = Path(args.output_report).resolve()

    if not baseline_dir.is_dir():
        print(
            f"ERROR: [diff_caches] baseline directory not found: {baseline_dir}",
            file=sys.stderr,
        )
        return 1
    if not builder_dir.is_dir():
        print(
            f"ERROR: [diff_caches] builder directory not found: {builder_dir}",
            file=sys.stderr,
        )
        return 1

    print(f"Comparing caches:")
    print(f"  Baseline: {baseline_dir}")
    print(f"  Builder:  {builder_dir}")

    print("  Comparing transcripts...")
    transcripts = diff_transcripts(baseline_dir, builder_dir)
    print(
        f"    {transcripts['shared']} shared, "
        f"{transcripts['matched']} matched, "
        f"{transcripts['differed']} differed, "
        f"rate={transcripts['match_rate']:.4f}"
    )

    print("  Comparing variations...")
    variations = diff_variations(baseline_dir, builder_dir)
    print(
        f"    {variations['shared']} shared, "
        f"{variations['matched']} matched, "
        f"{variations['differed']} differed, "
        f"rate={variations['match_rate']:.4f}"
    )

    print("  Comparing info.json...")
    info_json = diff_info_json(baseline_dir, builder_dir)
    print(f"    matches={info_json['matches']}")

    report = {
        "baseline_dir": str(baseline_dir),
        "builder_dir": str(builder_dir),
        "transcripts": transcripts,
        "variations": variations,
        "info_json": info_json,
    }

    output_report.parent.mkdir(parents=True, exist_ok=True)
    with open(output_report, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")

    print(f"\nReport written to: {output_report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
