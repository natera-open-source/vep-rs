#!/usr/bin/env python3
"""Field-level agreement between two VEP default-format outputs.

`compare_vep_outputs.py` scores the consequence key alone: (Location, Allele,
Feature, Feature_type, consequence set). This script takes the tuples both
engines emit on that key and compares every other column of the default output
on each of them: Uploaded_variation, Gene, cDNA_position, CDS_position,
Protein_position, Amino_acids, Codons, Existing_variation, and each KEY=VALUE
pair of the Extra column (IMPACT, DISTANCE, STRAND, FLAGS by default; every key
either side emits is compared). A pair whose key is on one side only is counted
and not compared, so a field disagreement is never conflated with a
consequence disagreement.

Usage:
    python3 compare_vep_fields.py --perl <vep_output.txt> --rust <vep_rs_output.txt> \
        --out <report.json> [--examples 20] [--tmp-dir DIR]

Both inputs stream through `sort` under C collation so the merge holds one line
per side in memory; a key repeated within one side keeps its first row and is
counted under `*_duplicate_keys`. Values are compared verbatim: `-` and an
empty string are different values, and no normalisation is applied to any
column, so a systematic formatting difference surfaces as its own example.

Report: per field, the compared count, the equal count, the differing count
and the most frequent (reference value, scored value) pairs.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import IO, Iterator

PLAIN_FIELDS = (
    "Uploaded_variation",
    "Gene",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
)
# Column index of each plain field in the 14-column default output.
PLAIN_INDEX = {
    "Uploaded_variation": 0,
    "Gene": 3,
    "cDNA_position": 7,
    "CDS_position": 8,
    "Protein_position": 9,
    "Amino_acids": 10,
    "Codons": 11,
    "Existing_variation": 12,
}
EXTRA_INDEX = 13
ABSENT = "<absent>"


def normalize_consequence_set(raw: str) -> str:
    terms = [t for t in raw.strip().split(",") if t]
    return ",".join(sorted(set(terms)))


def iter_rows(path: Path) -> Iterator[str]:
    """Yield one tab-joined line per output row: the 5 key columns, then the 9 value columns.

    Rows with fewer than 7 columns are skipped as the key comparator skips them;
    a row with 7 to 13 columns pads the missing value columns with an empty string.
    """
    opener = _opener(path)
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            while len(fields) < 14:
                fields.append("")
            key = (
                fields[1].strip(),
                fields[2].strip(),
                fields[4].strip(),
                fields[5].strip(),
                normalize_consequence_set(fields[6]),
            )
            values = [fields[PLAIN_INDEX[f]].strip() for f in PLAIN_FIELDS]
            values.append(fields[EXTRA_INDEX].strip())
            yield "\t".join(key) + "\t" + "\t".join(values)


def _opener(path: Path):
    if path.suffix == ".gz":
        import gzip

        return gzip.open
    return open


def sort_to_file(rows: Iterator[str], out_path: Path) -> None:
    """Sort rows bytewise so equal keys are adjacent and both sides order alike.

    The key occupies the first five tab-separated fields and a tab sorts below
    every printable byte, so a bytewise sort of the whole line orders the key
    prefix exactly as a field-wise sort would; the merge then compares key
    strings in Python, whose codepoint order equals C collation for this ASCII.
    """
    sort_exe = shutil.which("sort")
    if sort_exe is None:
        out_path.write_text("".join(f"{r}\n" for r in sorted(rows)), encoding="utf-8")
        return
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    env.pop("LC_COLLATE", None)
    with out_path.open("w", encoding="utf-8", newline="") as out_handle:
        proc = subprocess.Popen(
            [sort_exe, "-S", "25%"],
            stdin=subprocess.PIPE,
            stdout=out_handle,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
        assert proc.stdin is not None
        try:
            for row in rows:
                proc.stdin.write(row)
                proc.stdin.write("\n")
        finally:
            proc.stdin.close()
        proc.wait()
        err = proc.stderr.read() if proc.stderr else ""
        if proc.returncode != 0:
            raise RuntimeError(f"sort failed (rc={proc.returncode}): {err.strip()}")


def parse_extra(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if not raw or raw == "-":
        return out
    for part in raw.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = ""
    return out


class FieldStats:
    def __init__(self, examples: int) -> None:
        self.compared = 0
        self.equal = 0
        self.examples_cap = examples
        self.mismatch_shapes: Counter[tuple[str, str]] = Counter()

    def observe(self, perl_value: str, rust_value: str) -> None:
        self.compared += 1
        if perl_value == rust_value:
            self.equal += 1
        else:
            self.mismatch_shapes[(perl_value, rust_value)] += 1

    def report(self) -> dict:
        differ = self.compared - self.equal
        return {
            "compared": self.compared,
            "equal": self.equal,
            "differ": differ,
            "agreement": round(self.equal / self.compared, 6) if self.compared else None,
            "distinct_mismatch_shapes": len(self.mismatch_shapes),
            "top_mismatches": [
                {"reference": p, "scored": r, "count": n}
                for (p, r), n in self.mismatch_shapes.most_common(self.examples_cap)
            ],
        }


def _next_unique(handle: IO[str], pending: list[str | None]) -> tuple[str | None, list[str] | None, int]:
    """Return (key, values, duplicates_skipped) for the next distinct key, or (None, None, n) at EOF.

    `pending` holds a one-line lookahead so a run of equal keys can be collapsed
    to its first row without re-reading.
    """
    line = pending[0]
    pending[0] = None
    if line is None:
        line = handle.readline()
    if not line:
        return None, None, 0
    parts = line.rstrip("\n").split("\t")
    key = "\t".join(parts[:5])
    values = parts[5:]
    dups = 0
    while True:
        nxt = handle.readline()
        if not nxt:
            break
        if "\t".join(nxt.rstrip("\n").split("\t")[:5]) == key:
            dups += 1
            continue
        pending[0] = nxt
        break
    return key, values, dups


def merge(perl_sorted: Path, rust_sorted: Path, examples: int) -> dict:
    stats = {f: FieldStats(examples) for f in PLAIN_FIELDS}
    extra_stats: dict[str, FieldStats] = {}
    totals = {
        "reference_keys": 0,
        "scored_keys": 0,
        "intersection": 0,
        "reference_only": 0,
        "scored_only": 0,
        "reference_duplicate_keys": 0,
        "scored_duplicate_keys": 0,
    }
    with perl_sorted.open("r", encoding="utf-8") as ph, rust_sorted.open("r", encoding="utf-8") as rh:
        p_pending: list[str | None] = [None]
        r_pending: list[str | None] = [None]
        pk, pv, pd = _next_unique(ph, p_pending)
        rk, rv, rd = _next_unique(rh, r_pending)
        totals["reference_duplicate_keys"] += pd
        totals["scored_duplicate_keys"] += rd
        while pk is not None or rk is not None:
            if rk is None or (pk is not None and pk < rk):
                totals["reference_keys"] += 1
                totals["reference_only"] += 1
                pk, pv, pd = _next_unique(ph, p_pending)
                totals["reference_duplicate_keys"] += pd
            elif pk is None or rk < pk:
                totals["scored_keys"] += 1
                totals["scored_only"] += 1
                rk, rv, rd = _next_unique(rh, r_pending)
                totals["scored_duplicate_keys"] += rd
            else:
                totals["reference_keys"] += 1
                totals["scored_keys"] += 1
                totals["intersection"] += 1
                assert pv is not None and rv is not None
                for i, f in enumerate(PLAIN_FIELDS):
                    stats[f].observe(pv[i], rv[i])
                p_extra = parse_extra(pv[len(PLAIN_FIELDS)] if len(pv) > len(PLAIN_FIELDS) else "")
                r_extra = parse_extra(rv[len(PLAIN_FIELDS)] if len(rv) > len(PLAIN_FIELDS) else "")
                for k in sorted(set(p_extra) | set(r_extra)):
                    extra_stats.setdefault(k, FieldStats(examples)).observe(
                        p_extra.get(k, ABSENT), r_extra.get(k, ABSENT)
                    )
                pk, pv, pd = _next_unique(ph, p_pending)
                rk, rv, rd = _next_unique(rh, r_pending)
                totals["reference_duplicate_keys"] += pd
                totals["scored_duplicate_keys"] += rd
    return {
        "totals": totals,
        "fields": {f: stats[f].report() for f in PLAIN_FIELDS},
        "extra": {k: extra_stats[k].report() for k in sorted(extra_stats)},
    }


def compare(perl: Path, rust: Path, examples: int, tmp_dir: Path | None) -> dict:
    with tempfile.TemporaryDirectory(prefix="vep_fields_", dir=tmp_dir) as tmp:
        ps = Path(tmp) / "reference.sorted"
        rs = Path(tmp) / "scored.sorted"
        sort_to_file(iter_rows(perl), ps)
        sort_to_file(iter_rows(rust), rs)
        report = merge(ps, rs, examples)
    report["inputs"] = {"reference": str(perl), "scored": str(rust)}
    report["key"] = ["Location", "Allele", "Feature", "Feature_type", "consequence_set"]
    report["comparison"] = "verbatim per field on tuples both sides emit under the same key"
    return report


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perl", required=True, type=Path, help="reference engine default-format output")
    ap.add_argument("--rust", required=True, type=Path, help="scored engine default-format output")
    ap.add_argument("--out", required=True, type=Path, help="JSON report path")
    ap.add_argument("--examples", type=int, default=20, help="mismatch shapes kept per field")
    ap.add_argument("--tmp-dir", type=Path, default=None, help="directory for the sorted intermediates")
    args = ap.parse_args(argv)
    for p in (args.perl, args.rust):
        if not p.is_file():
            print(f"ERROR: [compare_vep_fields] not a file: {p}", file=sys.stderr)
            return 2
    report = compare(args.perl, args.rust, args.examples, args.tmp_dir)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    t = report["totals"]
    print(
        f"keys: reference {t['reference_keys']:,} scored {t['scored_keys']:,} "
        f"intersection {t['intersection']:,} (reference-only {t['reference_only']:,}, "
        f"scored-only {t['scored_only']:,})"
    )
    for name, block in list(report["fields"].items()) + [(f"Extra:{k}", v) for k, v in report["extra"].items()]:
        agr = block["agreement"]
        print(f"  {name:<24} equal {block['equal']:>12,} / {block['compared']:>12,}  agreement {agr}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
