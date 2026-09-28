#!/usr/bin/env python3
"""Extract one contig from a reference FASTA into its own indexed FASTA.

A golden corpus whose records were annotated with `--hgvs` needs the reference
both engines shifted their indels against: HGVS moves an insertion or deletion to
its most 3' position while the next reference base matches, so the corpus ships
the complete contig its records sit on (chromosome 21 for the GRCh37 HGVS corpus)
rather than a windowed extract, and the golden test passes the decompressed file to
`--fasta`.

The contig is written with the source's line width, byte for byte the sequence the
source carries, and a samtools-style `.fai` (name, length, offset, bases per line,
bytes per line) is written beside it so an indexed reader opens it without
re-indexing. Compress the FASTA afterwards (`gzip -9 reference.fa`); the `.fai`
describes the plain file and is committed as is.

Usage:
    extract_reference.py <reference.fa> --contig 21 --out reference.fa
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def extract(source: Path, contig: str, out: Path) -> tuple[int, int, int]:
    """Write `contig` from `source` to `out` and `out.fai`; return (length, line_bases, line_bytes)."""
    length = 0
    line_bases = 0
    line_bytes = 0
    found = False
    with source.open("rb") as src, out.open("wb") as dst:
        header = f">{contig}".encode()
        in_contig = False
        for line in src:
            if line.startswith(b">"):
                if in_contig:
                    break
                name = line[1:].split()[0] if len(line) > 1 else b""
                in_contig = name == contig.encode()
                if in_contig:
                    found = True
                    dst.write(header + b"\n")
                continue
            if not in_contig:
                continue
            bases = line.rstrip(b"\r\n")
            if not bases:
                continue
            if line_bases == 0:
                line_bases = len(bases)
                line_bytes = len(line) if line.endswith(b"\n") else len(line) + 1
            dst.write(bases + b"\n")
            length += len(bases)
    if not found:
        out.unlink(missing_ok=True)
        raise SystemExit(
            f"ERROR: [extract_reference] contig {contig!r} not found in {source}"
        )
    offset = len(header) + 1
    Path(str(out) + ".fai").write_text(
        f"{contig}\t{length}\t{offset}\t{line_bases}\t{line_bases + 1}\n",
        encoding="utf-8",
    )
    return length, line_bases, line_bases + 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("source", type=Path, help="reference FASTA holding the contig")
    ap.add_argument(
        "--contig", required=True, help="sequence name as the FASTA header spells it"
    )
    ap.add_argument(
        "--out",
        required=True,
        type=Path,
        help="output FASTA; its .fai is written beside it",
    )
    args = ap.parse_args(argv)
    length, line_bases, line_bytes = extract(args.source, args.contig, args.out)
    print(
        f"{args.contig}: {length} bases, {line_bases} bases per line ({line_bytes} bytes), written to {args.out}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
