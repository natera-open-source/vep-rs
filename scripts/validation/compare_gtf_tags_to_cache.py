#!/usr/bin/env python3
"""Compare the transcript attributes a GTF yields with those a JSON cache carries.

`vep-cache-builder --gtf` turns the GTF `tag` values on each `transcript` line
into the four cache attributes VEP keeps (`basic` or `gencode_basic` ->
`gencode_basic`, `gencode_primary`, `cds_start_NF`, `cds_end_NF`); VEP's own
cache dumper stores the same four from the Ensembl core database. This script reads a (gzipped) Ensembl GTF the way the
builder does and a JSON transcript cache (plain or gzipped shards) the way the
runtime does, and reports, per attribute code, how many transcripts agree on its
presence and lists every disagreement.

Usage::

    python3 scripts/validation/compare_gtf_tags_to_cache.py \
        --gtf Homo_sapiens.GRCh38.115.gtf.gz --cache <json_cache_dir> \
        [--out report.json] [--max-mismatches 50]

Exit status is 0 when every code agrees on every transcript present on both
sides, 1 otherwise, 2 on a usage error. Transcripts present on only one side
are counted, not treated as disagreements.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from collections import Counter, defaultdict

KEPT_TAGS = (
    ("basic", "gencode_basic"),
    ("gencode_basic", "gencode_basic"),
    ("gencode_primary", "gencode_primary"),
    ("cds_start_NF", "cds_start_NF"),
    ("cds_end_NF", "cds_end_NF"),
)
CODES = tuple(dict.fromkeys(code for _, code in KEPT_TAGS))


def _open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt", encoding="utf-8")


def read_gtf(path: str) -> dict[str, set[str]]:
    """Unversioned transcript id to the kept attribute codes its `transcript` line implies."""
    out: dict[str, set[str]] = {}
    tag_to_code = dict(KEPT_TAGS)
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            transcript_id = None
            codes: set[str] = set()
            for attr in fields[8].split(";"):
                attr = attr.strip()
                if " " not in attr:
                    continue
                key, value = attr.split(" ", 1)
                value = value.strip().strip('"')
                if key == "transcript_id":
                    transcript_id = value.split(".")[0]
                elif key == "tag" and value in tag_to_code:
                    codes.add(tag_to_code[value])
            if transcript_id is not None:
                out.setdefault(transcript_id, set()).update(codes)
    return out


def _iter_shards(cache_dir: str):
    root = os.path.join(cache_dir, "transcripts")
    if not os.path.isdir(root):
        root = cache_dir
    for dirpath, _dirs, files in os.walk(root):
        for name in sorted(files):
            if name.endswith(".json") or name.endswith(".json.gz"):
                yield os.path.join(dirpath, name)


def _transcripts(payload):
    if isinstance(payload, list):
        yield from (t for t in payload if t)
    elif isinstance(payload, dict):
        for value in payload.values():
            if isinstance(value, list):
                yield from (t for t in value if t)


def read_cache(cache_dir: str) -> dict[str, set[str]]:
    """Unversioned stable id to the kept attribute codes present on the transcript.

    A transcript stored in several shards is read once; the copies are identical.
    """
    out: dict[str, set[str]] = {}
    for shard in _iter_shards(cache_dir):
        with _open_text(shard) as fh:
            payload = json.load(fh)
        for tx in _transcripts(payload):
            stable_id = str(tx.get("stable_id", "")).split(".")[0]
            if not stable_id or stable_id in out:
                continue
            codes = {
                str(a.get("code"))
                for a in (tx.get("attributes") or [])
                if isinstance(a, dict) and a.get("code") in CODES
            }
            out[stable_id] = codes
    return out


def compare(gtf: dict[str, set[str]], cache: dict[str, set[str]], max_mismatches: int) -> dict:
    shared = sorted(set(gtf) & set(cache))
    agree = Counter()
    mismatches: dict[str, list] = defaultdict(list)
    for tid in shared:
        for code in CODES:
            in_gtf = code in gtf[tid]
            in_cache = code in cache[tid]
            if in_gtf == in_cache:
                agree[code] += 1
            elif len(mismatches[code]) < max_mismatches:
                mismatches[code].append({"transcript": tid, "gtf": in_gtf, "cache": in_cache})
    disagree = {code: len(shared) - agree[code] for code in CODES}
    return {
        "transcripts": {
            "gtf": len(gtf),
            "cache": len(cache),
            "shared": len(shared),
            "gtf_only": len(set(gtf) - set(cache)),
            "cache_only": len(set(cache) - set(gtf)),
        },
        "per_code": {
            code: {
                "agree": agree[code],
                "disagree": disagree[code],
                "present_in_gtf": sum(1 for t in shared if code in gtf[t]),
                "present_in_cache": sum(1 for t in shared if code in cache[t]),
                "examples": mismatches[code],
            }
            for code in CODES
        },
        "exact_agreement": all(v == 0 for v in disagree.values()),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--cache", required=True, help="JSON cache directory (holds transcripts/<chr>/)")
    ap.add_argument("--out", help="write the JSON report here as well as printing a summary")
    ap.add_argument("--max-mismatches", type=int, default=50, help="examples kept per code")
    args = ap.parse_args(argv)
    if not os.path.exists(args.gtf) or not os.path.isdir(args.cache):
        print("ERROR: [compare_gtf_tags_to_cache] --gtf must be a file and --cache a directory", file=sys.stderr)
        return 2
    report = compare(read_gtf(args.gtf), read_cache(args.cache), args.max_mismatches)
    t = report["transcripts"]
    print(f"transcripts: gtf {t['gtf']:,} cache {t['cache']:,} shared {t['shared']:,} "
          f"gtf_only {t['gtf_only']:,} cache_only {t['cache_only']:,}")
    for code, row in report["per_code"].items():
        print(f"  {code:<16} agree {row['agree']:>9,}  disagree {row['disagree']:>7,}  "
              f"present gtf {row['present_in_gtf']:>8,} cache {row['present_in_cache']:>8,}")
        for ex in row["examples"][:5]:
            print(f"    {ex['transcript']}: gtf={ex['gtf']} cache={ex['cache']}")
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=1)
            fh.write("\n")
    return 0 if report["exact_agreement"] else 1


if __name__ == "__main__":
    sys.exit(main())
