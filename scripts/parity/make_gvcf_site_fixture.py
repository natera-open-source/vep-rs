#!/usr/bin/env python3
"""Write the gVCF variant-site fixture `gvcf_site_rows_fixture.tsv` beside this script.

usage: make_gvcf_site_fixture.py --vep BIN --json-cache DIR --out TSV [--seed N] [--sites N]

A gVCF variant site lists `<NON_REF>` beside its real alleles (`A G,<NON_REF>`). vep-rs reads
the site on its sequence path, and the comparators predict the (Location, Allele) it writes for
each real allele (`gvcf_site_rows` in scripts/concordance/compare_vep_outputs.py). This fixture
holds what the binary wrote for generated sites, so the test of that prediction holds it to the
binary without running one.

The sites are drawn from `--seed`: on chromosome 21 between 25,000,000 and 46,000,000, each at
its own position, a REF of 1 to 4 bases over ACGT, 1 to 3 real ALTs of 1 to 5 bases over ACGT
that differ from REF and from each other, every ALT sharing REF's first base in every second
site, and `<NON_REF>` at a random position of the ALT list; no `*`, no `.`, no empty field. The
binary annotates them from a JSON cache holding chromosome 21 of GRCh38 (`--offline
--json_cache DIR --assembly GRCh38 --force_overwrite --no_stats --tab --allele_number`), and
the fixture gets one row per (site, Location, Allele) the binary wrote, sorted, under the
columns `pos`, `ref`, `alt_column`, `allele`, `location`. A site the binary wrote no row for is
an error, never a silent gap.

Regenerate with the defaults after a change to the parser's sequence path:

    scripts/parity/make_gvcf_site_fixture.py --vep target/debug/vep --json-cache DIR \\
        --out scripts/parity/gvcf_site_rows_fixture.tsv
"""

from __future__ import annotations

import argparse
import csv
import random
import subprocess
import sys
import tempfile
from pathlib import Path

CHROM = "21"
POSITION_RANGE = (25_000_000, 46_000_000)
BASES = "ACGT"
NON_REF = "<NON_REF>"
COLUMNS = ("pos", "ref", "alt_column", "allele", "location")
VEP_FLAGS = ("--offline", "--assembly", "GRCh38", "--force_overwrite", "--no_stats", "--tab", "--allele_number")


def draw_sites(rng: random.Random, count: int) -> list[tuple[int, str, str]]:
    """``count`` sites as (pos, ref, alt_column), sorted by position, each position drawn once."""
    sites: dict[int, tuple[str, str]] = {}
    while len(sites) < count:
        pos = rng.randint(*POSITION_RANGE)
        if pos in sites:
            continue
        ref = "".join(rng.choice(BASES) for _ in range(rng.randint(1, 4)))
        anchored = len(sites) % 2 == 0
        alts: list[str] = []
        wanted = rng.randint(1, 3)
        while len(alts) < wanted:
            length = rng.randint(1, 5)
            alt = "".join(rng.choice(BASES) for _ in range(length))
            if anchored:
                alt = ref[0] + alt[1:]
            if alt != ref and alt not in alts:
                alts.append(alt)
        alts.insert(rng.randint(0, len(alts)), NON_REF)
        sites[pos] = (ref, ",".join(alts))
    return [(pos, ref, alt_column) for pos, (ref, alt_column) in sorted(sites.items())]


def write_vcf(path: Path, sites: list[tuple[int, str, str]]) -> dict[str, tuple[int, str, str]]:
    """Write the sites as VCF records named `s<n>`; returns the name of each site."""
    names: dict[str, tuple[int, str, str]] = {}
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for n, (pos, ref, alt_column) in enumerate(sites, 1):
            name = f"s{n:04d}"
            names[name] = (pos, ref, alt_column)
            handle.write(f"{CHROM}\t{pos}\t{name}\t{ref}\t{alt_column}\t.\t.\t.\n")
    return names


def run_vep(vep: Path, json_cache: Path, vcf: Path, out: Path) -> None:
    command = [str(vep), *VEP_FLAGS, "--json_cache", str(json_cache), "-i", str(vcf), "-o", str(out)]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        sys.exit(f"ERROR: [make_gvcf_site_fixture] {vep} exited {result.returncode}:\n{result.stderr}")


def written_rows(
    output: Path, names: dict[str, tuple[int, str, str]]
) -> list[tuple[int, str, str, str, str]]:
    """The distinct (pos, ref, alt_column, allele, location) rows of the binary's tab output."""
    rows: set[tuple[int, str, str, str, str]] = set()
    seen: set[str] = set()
    with output.open(encoding="utf-8") as handle:
        header: list[str] | None = None
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = fields
                continue
            if header is None:
                sys.exit("ERROR: [make_gvcf_site_fixture] the output has no column header")
            record = dict(zip(header, fields))
            name = record["#Uploaded_variation"]
            if name not in names:
                sys.exit(f"ERROR: [make_gvcf_site_fixture] the output names a record the input lacks: {name}")
            seen.add(name)
            pos, ref, alt_column = names[name]
            rows.add((pos, ref, alt_column, record["Allele"], record["Location"]))
    missing = sorted(set(names) - seen)
    if missing:
        sys.exit(f"ERROR: [make_gvcf_site_fixture] no row written for {len(missing)} site(s): {' '.join(missing)}")
    return sorted(rows)


def write_fixture(path: Path, rows: list[tuple[int, str, str, str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--vep", required=True, type=Path, help="the vep-rs binary")
    parser.add_argument("--json-cache", required=True, type=Path, help="a JSON cache holding chromosome 21 of GRCh38")
    parser.add_argument("--out", required=True, type=Path, help="the fixture TSV to write")
    parser.add_argument("--seed", type=int, default=7919, help="the seed the sites are drawn from")
    parser.add_argument("--sites", type=int, default=300, help="the number of sites")
    args = parser.parse_args(argv)
    sites = draw_sites(random.Random(args.seed), args.sites)
    with tempfile.TemporaryDirectory() as tmp:
        vcf, output = Path(tmp) / "sites.vcf", Path(tmp) / "sites.txt"
        names = write_vcf(vcf, sites)
        run_vep(args.vep, args.json_cache, vcf, output)
        rows = written_rows(output, names)
    write_fixture(args.out, rows)
    print(f"wrote {len(rows)} rows for {len(sites)} sites to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
