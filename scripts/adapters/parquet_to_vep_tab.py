#!/usr/bin/env python3
"""Convert vep-rs Parquet output back to the `--tab` format.

The Parquet schema is the five key columns (`chrom`, `pos`, `end`, `ref`, `alt`)
followed by exactly the `--tab` columns, with `-` stored as NULL and `pos`, `end`,
`DISTANCE` and `STRAND` stored as integers. This adapter projects the tab columns
back out: NULL becomes `-`, integers become text, and a nested-shape file is
unnested so every (variant x consequence) row is one line again. The header is
the `#`-prefixed column line; the `##` meta lines of a live run are not stored
in Parquet and are not reproduced.

Rows come out in Parquet order (chrom, pos, ref, alt), not input order; compare
against a `--tab` file as a multiset of lines.

Shape detection: nested (LIST columns) vs flat (scalar) is probed from the file's
schema. A Hive-partitioned directory (`<dir>/chrom=<c>/*.parquet`) or a single
file are both accepted.

Usage::

    python3 scripts/adapters/parquet_to_vep_tab.py --input out.parquet [--output out.tab]
    python3 scripts/adapters/parquet_to_vep_tab.py out.parquet    # to stdout

Requires ``duckdb`` on PATH. No other Python deps (pure stdlib + subprocess).
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import shutil
import subprocess
import sys

KEY_COLUMNS = ["chrom", "pos", "end", "ref", "alt"]


def _duckdb_available() -> bool:
    return shutil.which("duckdb") is not None


def _escape_sql_literal(s: str) -> str:
    return s.replace("'", "''")


def _source(parquet_path: str) -> str:
    """A `read_parquet` expression over a file or a partitioned directory.

    Hive partition columns are not re-derived from the path (they are stored in
    the files as VARCHAR), so `chrom` keeps its stored type.
    """
    if os.path.isdir(parquet_path):
        glob = os.path.join(parquet_path, "**", "*.parquet")
        return f"read_parquet('{_escape_sql_literal(glob)}', hive_partitioning=false)"
    return f"read_parquet('{_escape_sql_literal(parquet_path)}', hive_partitioning=false)"


def _run_duckdb(sql: str) -> str:
    proc = subprocess.run(
        ["duckdb", "-csv", "-noheader", "-c", sql],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        print(f"ERROR: [parquet_to_vep_tab] duckdb failed: {proc.stderr.strip()}", file=sys.stderr)
        raise SystemExit(2)
    return proc.stdout


def _schema(parquet_path: str) -> list[tuple[str, str]]:
    """(column name, DuckDB type) pairs in stored order."""
    out = _run_duckdb(f"DESCRIBE SELECT * FROM {_source(parquet_path)}")
    rows = list(csv.reader(io.StringIO(out)))
    return [(r[0], r[1]) for r in rows if r]


def _build_sql(parquet_path: str, schema: list[tuple[str, str]]) -> str:
    """SQL producing the tab columns as text, one line per consequence row."""
    tab_cols = [c for c, _ in schema if c not in KEY_COLUMNS]
    types = dict(schema)
    nested = any(types[c].endswith("[]") for c in tab_cols)
    if nested:
        # Explode every LIST column in lock-step; scalar per-variant columns repeat.
        list_cols = [c for c in tab_cols if types[c].endswith("[]")]
        first = list_cols[0]
        exprs = []
        for c in tab_cols:
            if c in list_cols:
                exprs.append(f'COALESCE(CAST("{c}"[i] AS VARCHAR), \'-\') AS "{c}"')
            else:
                exprs.append(f'COALESCE(CAST("{c}" AS VARCHAR), \'-\') AS "{c}"')
        return (
            f"SELECT {', '.join(exprs)} FROM {_source(parquet_path)}, "
            f'generate_series(1, greatest(len("{first}"), 1)) AS g(i) '
            "ORDER BY chrom, pos, ref, alt, i"
        )
    exprs = [f'COALESCE(CAST("{c}" AS VARCHAR), \'-\') AS "{c}"' for c in tab_cols]
    return f"SELECT {', '.join(exprs)} FROM {_source(parquet_path)} ORDER BY chrom, pos, ref, alt"


def convert(parquet_path: str, out) -> int:
    schema = _schema(parquet_path)
    tab_cols = [c for c, _ in schema if c not in KEY_COLUMNS]
    if not tab_cols:
        print("ERROR: [parquet_to_vep_tab] no tab columns found in the Parquet schema", file=sys.stderr)
        return 2
    sql = _build_sql(parquet_path, schema)
    proc = subprocess.run(
        ["duckdb", "-c", f"COPY ({sql}) TO '/dev/stdout' (FORMAT CSV, DELIMITER '\\t', HEADER false, QUOTE '', ESCAPE '')"],
        check=False,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        print(f"ERROR: [parquet_to_vep_tab] duckdb failed: {proc.stderr.strip()}", file=sys.stderr)
        return 2
    out.write("#" + "\t".join(tab_cols) + "\n")
    out.write(proc.stdout)
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", nargs="?", help="Parquet file or partitioned directory")
    ap.add_argument("--input", "-i", dest="input_path", help="Parquet file or partitioned directory")
    ap.add_argument("--output", "-o", default="-", help="output path (default stdout)")
    args = ap.parse_args(argv)
    path = args.input_path or args.path
    if not path:
        ap.error("a Parquet path is required")
    if not _duckdb_available():
        print("ERROR: [parquet_to_vep_tab] duckdb CLI not found on PATH", file=sys.stderr)
        return 2
    if args.output == "-":
        return convert(path, sys.stdout)
    with open(args.output, "w", encoding="utf-8") as fh:
        return convert(path, fh)


if __name__ == "__main__":
    sys.exit(main())
