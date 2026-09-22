#!/usr/bin/env python3
"""Compare two VEP Storable cache dumps transcript by transcript.

`compare_storable_cache_versions.pl` writes one line per distinct transcript of a cache
(`stable_id contig start end md5_full md5_content`). This script joins two such dumps on
`stable_id` and reports the ids in one cache only, the ids in both whose full dump is
identical, those identical in content once database bookkeeping and number quoting are
normalised (the Perl script's second digest), and those that differ in content. It writes
the id lists beside a summary JSON and, with `--parity-csv`, extends the matching assembly
row of the transcript-set parity CSV with the counts (columns are added on first use and
overwritten on a re-run).

Usage:
  compare_storable_cache_versions.py --version-a 113 --dump-a 113.tsv --version-b 115 --dump-b 115.tsv
      --out-dir <dir> [--parity-csv manuscript/data/transcript_set_parity.csv --assembly GRCh37]
"""

from __future__ import annotations

import argparse
import csv
import re
import json
import sys
from pathlib import Path


def load(path: Path) -> dict[str, tuple[str, str, str, str, str]]:
    out: dict[str, tuple[str, str, str, str, str]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 6:
                raise SystemExit(f"ERROR: [compare_storable_cache_versions] malformed dump line in {path}: {line!r}")
            out[cols[0]] = (cols[1], cols[2], cols[3], cols[4], cols[5])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version-a", required=True)
    ap.add_argument("--dump-a", required=True, type=Path)
    ap.add_argument("--version-b", required=True)
    ap.add_argument("--dump-b", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--parity-csv", type=Path)
    ap.add_argument("--sweep-id", help="id of the run that produced the two dumps (sweep-<UTC>), written to the parity row")
    ap.add_argument("--assembly")
    a = ap.parse_args()
    if a.sweep_id and not re.fullmatch(r"sweep-\d{8}T\d{6}Z", a.sweep_id):
        print(f"ERROR: [compare_storable_cache_versions] sweep id must be sweep-<UTC>: {a.sweep_id}", file=sys.stderr)
        return 2
    if bool(a.parity_csv) != bool(a.assembly):
        print("ERROR: [compare_storable_cache_versions] --parity-csv and --assembly go together", file=sys.stderr)
        return 2

    da, db = load(a.dump_a), load(a.dump_b)
    only_a = sorted(set(da) - set(db))
    only_b = sorted(set(db) - set(da))
    both = sorted(set(da) & set(db))
    identical_full = [i for i in both if da[i][3] == db[i][3]]
    identical_stripped_only = [i for i in both if da[i][3] != db[i][3] and da[i][4] == db[i][4]]
    differing = [i for i in both if da[i][4] != db[i][4]]
    va, vb = a.version_a, a.version_b
    a.out_dir.mkdir(parents=True, exist_ok=True)
    for name, ids in ((f"only_in_{va}", only_a), (f"only_in_{vb}", only_b),
                      (f"identical_after_strip_{va}_{vb}", identical_stripped_only), (f"differing_{va}_{vb}", differing)):
        with open(a.out_dir / f"{name}.tsv", "w", encoding="utf-8") as fh:
            fh.write("transcript_id\tcontig\tstart\tend\n")
            for i in ids:
                src = da.get(i) or db[i]
                fh.write("\t".join((i, src[0], src[1], src[2])) + "\n")
    summary = {
        "version_a": va, "version_b": vb,
        f"storable_{va}_transcripts": len(da), f"storable_{vb}_transcripts": len(db),
        "in_both": len(both), f"only_in_{va}": len(only_a), f"only_in_{vb}": len(only_b),
        "identical_full_dump": len(identical_full), "identical_content_only": len(identical_stripped_only),
        "differing_content": len(differing),
        "content_digest_normalises": ["dbID", "adaptor", "created_date", "modified_date", "_gene_version", "integer scalars stored as strings"],
        "content_identical": not only_a and not only_b and not differing,
    }
    (a.out_dir / f"summary_{va}_vs_{vb}.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[compare_storable_cache_versions] {va}: {len(da):,} transcripts; {vb}: {len(db):,}; both {len(both):,}; "
          f"only {va} {len(only_a):,}; only {vb} {len(only_b):,}; identical {len(identical_full):,}; "
          f"identical after strip {len(identical_stripped_only):,}; differing {len(differing):,}")

    if a.parity_csv:
        with open(a.parity_csv, newline="") as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
            header = list(reader.fieldnames or [])
        extra = {
            f"storable_{va}_transcripts": len(da), f"storable_{vb}_transcripts": len(db),
            f"storable_{va}_{vb}_identical_transcripts": len(identical_full) + len(identical_stripped_only),
            f"storable_{va}_{vb}_differing_transcripts": len(differing),
            f"storable_{va}_only": len(only_a), f"storable_{vb}_only": len(only_b),
        }
        if a.sweep_id:
            extra[f"storable_{va}_{vb}_sweep_id"] = a.sweep_id
        for col in extra:
            if col not in header:
                header.append(col)
        hit = [r for r in rows if r["assembly"] == a.assembly]
        if len(hit) != 1:
            print(f"ERROR: [compare_storable_cache_versions] {a.parity_csv} has {len(hit)} rows for {a.assembly}", file=sys.stderr)
            return 1
        for k, v in extra.items():
            hit[0][k] = str(v)
        tmp = a.parity_csv.with_suffix(".csv.tmp")
        with open(tmp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=header, lineterminator="\n")
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in header})
        tmp.replace(a.parity_csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
