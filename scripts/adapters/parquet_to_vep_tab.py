#!/usr/bin/env python3
"""Project vep-rs Parquet output back to the `--tab` format.

A Parquet output is one file, a directory of `part-*.parquet` files, or a directory of
`chrom=<c>/` directories holding them. Every file carries the key columns (`contig_index`,
`chrom`, `pos`, `end`, `ref`, `alt`) and then exactly the `--tab` columns, each stored under
the type the writer's cell contract declares (the module documentation of
`crates/vep-io/src/output/parquet`). This adapter inverts that contract cell by cell:

* a NULL cell prints `-`; a string prints as stored (so the allele columns keep a stored
  `-`, an empty allele there, and an empty cell stays empty);
* an integer prints in decimal; a double prints as the shortest decimal that reads back to
  the same value (`0.095`, never `0.0950`; `0.000098`, never `9.8e-05`);
* a boolean prints its column's spelling: `YES` for `CANONICAL`, `1` for `PICK` and
  `GENCODE_PRIMARY` (false prints `-` for these three), `Y` and `N` for `HIGH_INF_POS`,
  `1` and `0` inside a list (`SOMATIC`, `PHENO`);
* a list prints its elements `,`-joined, a NULL element as `-`;
* a positions struct (`EXON`, `INTRON`, `cDNA_position`, `CDS_position`,
  `Protein_position`) prints `start` when `start` and `end` agree, `start-end` otherwise
  with `?` for a NULL member, `-` when both are NULL, then `/total` when `total` is set.

The output is the `#`-prefixed header line of the tab columns, then one line per row in
file order, the files of a directory in sorted path order; the `##` meta lines of a live
run are not stored in Parquet and are not reproduced. Rows are in the file's order, not
the input's: compare against a `--tab` file as a multiset of lines.

`--canonical` prints the canonical form instead: the header line, then the data lines
sorted, every DOUBLE cell as its shortest round-trip decimal and every INT32 cell as its
integer. The input may then also be a `--tab` output (plain or gzip), whose `##` lines are
dropped and whose DOUBLE-typed cells (the population frequency columns of
`DOUBLE_LIST_COLUMNS`, `MAX_AF` and `MOTIF_SCORE_CHANGE`) and INT32-typed cells
(`INT32_COLUMNS`) are reprinted the same way; two outputs with equal canonical forms hold
the same rows.

Usage::

    python3 scripts/adapters/parquet_to_vep_tab.py out.parquet [--output out.tab]
    python3 scripts/adapters/parquet_to_vep_tab.py --canonical out.parquet
    python3 scripts/adapters/parquet_to_vep_tab.py --canonical out.tab

Requires ``pyarrow`` for a Parquet input.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re
import sys
from collections.abc import Iterable, Iterator
from decimal import Decimal
from pathlib import Path
from typing import IO

COMPONENT = "parquet_to_vep_tab"

KEY_COLUMNS = ("contig_index", "chrom", "pos", "end", "ref", "alt")

# The `--tab` spelling of a BOOLEAN column's values: true, then false.
BOOLEAN_SPELLINGS = {
    "CANONICAL": ("YES", "-"),
    "PICK": ("1", "-"),
    "GENCODE_PRIMARY": ("1", "-"),
    "HIGH_INF_POS": ("Y", "N"),
}
# The spelling of a boolean inside a list and of a BOOLEAN column of no listed spelling.
FLAG_SPELLING = ("1", "0")

# The cache columns behind the population frequency fields, in output-field order; the
# frequency column of each is `<column>_AF`, `AF` for `AF`.
FREQUENCY_POPULATIONS = (
    "AF",
    "AFR",
    "AMR",
    "ASN",
    "EAS",
    "EUR",
    "SAS",
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_MID",
    "gnomADe_NFE",
    "gnomADe_REMAINING",
    "gnomADe_SAS",
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_REMAINING",
    "gnomADg_SAS",
)
DOUBLE_LIST_COLUMNS = frozenset(
    "AF" if population == "AF" else f"{population}_AF" for population in FREQUENCY_POPULATIONS
)
DOUBLE_COLUMNS = frozenset({"MAX_AF", "MOTIF_SCORE_CHANGE"})
INT32_COLUMNS = frozenset(
    {"STRAND", "DISTANCE", "ALLELE_NUM", "TSL", "HGVS_OFFSET", "OverlapBP", "MOTIF_POS"}
)
INT32_RANGE = range(-(2**31), 2**31)

NULL_TEXT = "-"
PARQUET_MAGIC = b"PAR1"
GZIP_MAGIC = b"\x1f\x8b"

# A decimal number as the writer's cell contract parses it: an optional sign, digits with
# an optional fraction or a fraction alone, an optional exponent; or an infinity or NaN.
_DOUBLE = re.compile(
    r"^[+-]?(?:(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|inf|infinity|nan)$", re.IGNORECASE
)
# An integer as the writer's INT32 parse reads it: an optional sign and ASCII digits, no blank.
_INT32 = re.compile(r"^([+-]?)([0-9]+)\Z")


class AdapterError(Exception):
    """A failure the command reports on stderr and exits 2 on."""


def shortest_decimal(value: float) -> str:
    """The shortest decimal that reads back to ``value``, in plain notation: ``0.095`` for
    0.095, ``0.000098`` for 9.8e-05, ``1`` for 1.0, ``-0`` for -0.0.

    Of two shortest decimals equally close to the value (a 17-digit case), the writer's
    formatter prints the one of larger magnitude where Python's ``repr`` prints the even one,
    so that tie is re-settled here."""
    if math.isnan(value):
        return "NaN"
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    shortest = Decimal(repr(value))
    exact = Decimal(value)
    unit = Decimal(1).scaleb(shortest.as_tuple().exponent)
    away = shortest - unit if value < 0 else shortest + unit
    if abs(away - exact) == abs(shortest - exact):
        shortest = away
    return format(shortest.normalize(), "f")


def parse_double(token: str) -> float | None:
    """``token`` as a double under the contract's number grammar, ``None`` outside it."""
    if _DOUBLE.match(token) is None:
        return None
    return float(token)


def canonical_number(token: str) -> str:
    """A DOUBLE cell's token reprinted as its shortest decimal; a token outside the grammar
    stays as written."""
    value = parse_double(token)
    return token if value is None else shortest_decimal(value)


def canonical_int32(token: str) -> str:
    """An INT32 cell's token reprinted as its integer (``007`` is ``7``), as the writer stores
    it; a token that is no INT32 (``2147483648``, ``1.5``) is stored as NULL and prints so."""
    match = _INT32.match(token)
    if match is None:
        return NULL_TEXT
    sign, digits = match.groups()
    digits = digits.lstrip("0") or "0"
    # Ten digits hold every INT32, so a longer run is outside the type before it is
    # converted, which keeps the conversion within the interpreter's digit limit.
    if len(digits) > 10:
        return NULL_TEXT
    value = int(sign + digits)
    if value not in INT32_RANGE:
        return NULL_TEXT
    return str(value)


def canonical_tab_cell(column: str, text: str) -> str:
    """A ``--tab`` cell in canonical form: a DOUBLE or ``LIST<DOUBLE>`` cell with each number
    as its shortest decimal, an INT32 cell as its integer, every other cell unchanged."""
    if text == NULL_TEXT:
        return text
    if column in DOUBLE_COLUMNS:
        return canonical_number(text)
    if column in DOUBLE_LIST_COLUMNS:
        return ",".join(canonical_number(token) for token in text.split(","))
    if column in INT32_COLUMNS:
        return canonical_int32(text)
    return text


def _flag_text(value: bool, spelling: tuple[str, str]) -> str:
    return spelling[0] if value else spelling[1]


def _positions_text(members: dict) -> str:
    start, end, total = members.get("start"), members.get("end"), members.get("total")
    if start is None and end is None:
        text = NULL_TEXT
    elif start == end:
        text = str(start)
    else:
        text = f"{'?' if start is None else start}-{'?' if end is None else end}"
    if total is not None:
        text += f"/{total}"
    return text


def _element_text(value) -> str:
    if value is None:
        return NULL_TEXT
    if isinstance(value, bool):
        return _flag_text(value, FLAG_SPELLING)
    if isinstance(value, float):
        return shortest_decimal(value)
    return str(value)


def cell_text(column: str, value) -> str:
    """The ``--tab`` cell of one stored value, as ``pyarrow`` hands it back: ``None`` for
    NULL, ``str``, ``int``, ``float``, ``bool``, a ``list`` of those or a ``dict`` of the
    positions members."""
    if value is None:
        return NULL_TEXT
    if isinstance(value, bool):
        return _flag_text(value, BOOLEAN_SPELLINGS.get(column, FLAG_SPELLING))
    if isinstance(value, float):
        return shortest_decimal(value)
    if isinstance(value, list):
        return ",".join(_element_text(element) for element in value)
    if isinstance(value, dict):
        return _positions_text(value)
    return str(value)


def parquet_files(path: Path) -> list[Path]:
    """The Parquet files of an output: the file itself, or every ``.parquet`` file under the
    directory (a ``chrom=<c>/`` directory included) in sorted path order."""
    if path.is_file():
        return [path]
    if not path.is_dir():
        reason = "no such file or directory" if not path.exists() else "neither a regular file nor a directory"
        raise AdapterError(f"{path}: {reason}")
    files = sorted(p for p in path.rglob("*.parquet") if p.is_file())
    if not files:
        raise AdapterError(f"{path}: no .parquet file under the directory")
    return files


def tab_columns(schema) -> list[str]:
    """The tab columns of a file's schema: every field after the leading key columns."""
    names = schema.names
    skip = 0
    while skip < len(names) and names[skip] in KEY_COLUMNS:
        skip += 1
    columns = names[skip:]
    if not columns:
        raise AdapterError("no tab columns in the Parquet schema")
    return columns


def _pyarrow_parquet():
    try:
        import pyarrow.parquet as pq
    except ImportError as error:
        raise AdapterError(f"pyarrow is not installed ({error}); install it to read Parquet") from error
    return pq


def parquet_rows(path: Path) -> tuple[list[str], Iterator[str]]:
    """The tab columns of a Parquet output and its rows as ``--tab`` lines, in file order."""
    pq = _pyarrow_parquet()
    files = parquet_files(path)
    readers = [pq.ParquetFile(file) for file in files]
    columns = tab_columns(readers[0].schema_arrow)
    for file, reader in zip(files[1:], readers[1:]):
        if tab_columns(reader.schema_arrow) != columns:
            raise AdapterError(f"{file}: its tab columns differ from {files[0]}'s")

    def lines() -> Iterator[str]:
        for reader in readers:
            for batch in reader.iter_batches(columns=columns):
                cells = [batch.column(i).to_pylist() for i in range(len(columns))]
                for row in zip(*cells):
                    yield "\t".join(cell_text(column, value) for column, value in zip(columns, row))

    return columns, lines()


def _read_text(path: Path) -> str:
    with open(path, "rb") as fh:
        head = fh.read(len(GZIP_MAGIC))
    opener = gzip.open if head == GZIP_MAGIC else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        return fh.read()


def tab_rows(path: Path) -> tuple[list[str], Iterator[str]]:
    """The columns of a ``--tab`` output and its data lines; the ``##`` lines are dropped. A
    line ends at ``\\n`` alone (``\\r`` before it dropped), so a cell holding another
    separator character stays in its row."""
    header = None
    data: list[str] = []
    for line in _read_text(path).split("\n"):
        line = line.rstrip("\r")
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            if header is not None:
                raise AdapterError(f"{path}: two header lines")
            header = line[1:].split("\t")
            continue
        if line:
            data.append(line)
    if header is None:
        raise AdapterError(f"{path}: no # header line; a tab output needs its header")
    return header, iter(data)


def is_parquet(path: Path) -> bool:
    """Whether ``path`` is a directory or a file opening with the Parquet magic bytes."""
    if path.is_dir():
        return True
    with open(path, "rb") as fh:
        return fh.read(len(PARQUET_MAGIC)) == PARQUET_MAGIC


def canonical_lines(columns: list[str], lines: Iterable[str]) -> list[str]:
    """The canonical data lines: every DOUBLE cell as its shortest decimal and every INT32
    cell as its integer, sorted; a cell past the last header column stays as written."""
    canonical = []
    for line in lines:
        cells = line.split("\t")
        for index, column in enumerate(columns[: len(cells)]):
            cells[index] = canonical_tab_cell(column, cells[index])
        canonical.append("\t".join(cells))
    canonical.sort()
    return canonical


def write_tab(path: Path, out: IO[str], canonical: bool) -> None:
    """Writes the ``--tab`` projection of the output at ``path`` to ``out``."""
    if not path.exists():
        raise AdapterError(f"{path}: no such file or directory")
    if not (path.is_file() or path.is_dir()):
        raise AdapterError(f"{path}: neither a regular file nor a directory")
    parquet = is_parquet(path)
    if canonical and not parquet:
        columns, lines = tab_rows(path)
    elif not parquet:
        raise AdapterError(f"{path}: not a Parquet file")
    else:
        columns, lines = parquet_rows(path)
    out.write("#" + "\t".join(columns) + "\n")
    if canonical:
        lines = canonical_lines(columns, lines)
    for line in lines:
        out.write(line + "\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", nargs="?", help="Parquet file or directory (with --canonical, a --tab output too)")
    ap.add_argument("--input", "-i", dest="input_path", help="the same as the positional path")
    ap.add_argument("--output", "-o", default="-", help="output path (default stdout)")
    ap.add_argument(
        "--canonical",
        action="store_true",
        help="the canonical form: the header, then the data lines sorted with every DOUBLE cell "
        "as its shortest round-trip decimal and every INT32 cell as its integer; the input may "
        "be a --tab output",
    )
    args = ap.parse_args(argv)
    path = args.input_path or args.path
    if not path:
        ap.error("a Parquet path is required")
    try:
        if args.output == "-":
            write_tab(Path(path), sys.stdout, args.canonical)
            sys.stdout.flush()
        else:
            with open(args.output, "w", encoding="utf-8") as fh:
                write_tab(Path(path), fh, args.canonical)
    except BrokenPipeError:
        # The reader closed the pipe (`| head`): nothing more is owed, and stdout is pointed
        # at /dev/null so the interpreter's flush at exit does not raise again.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 0
    except (AdapterError, OSError, ValueError, EOFError) as error:
        print(f"ERROR: [{COMPONENT}] {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
