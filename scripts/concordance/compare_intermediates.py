#!/usr/bin/env python3
"""Compare Perl VEP intermediates with Rust VEP output to pinpoint divergences.

Joins Perl tracer TSV (from trace_perl_intermediates.pl) with Rust VEP output
on (location, allele, transcript_id) and reports which fields diverge.

Usage:
    python3 compare_intermediates.py \
        --perl-trace perl_intermediates.tsv \
        --rust-output rust_output.txt \
        [--discordant-only]

Output: per-row comparison + aggregate field-divergence stats.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass
class PerlRow:
    location: str
    allele: str
    transcript_id: str
    cds_start: str
    cds_end: str
    translation_start: str
    translation_end: str
    ref_codon: str
    alt_codon: str
    ref_peptide: str
    alt_peptide: str
    consequence_set: str


@dataclass
class RustRow:
    location: str
    allele: str
    transcript_id: str
    consequence_set: str
    amino_acids: str  # "X/Y" format
    codons: str       # "aGt/aTt" format
    cds_position: str
    protein_position: str


def parse_perl_trace(path: Path) -> dict[tuple[str, str, str], PerlRow]:
    """Parse Perl tracer TSV into dict keyed by (location, allele, transcript_id)."""
    rows: dict[tuple[str, str, str], PerlRow] = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for rec in reader:
            key = (rec["location"], rec["allele"], rec["transcript_id"])
            rows[key] = PerlRow(**{k: rec.get(k, "") for k in PerlRow.__dataclass_fields__})
    return rows


def parse_rust_output(path: Path) -> dict[tuple[str, str, str], RustRow]:
    """Parse Rust VEP default output into dict keyed by (location, allele, transcript_id)."""
    rows: dict[tuple[str, str, str], RustRow] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 14:
                continue

            # VEP default output columns:
            # 0: Uploaded_variation, 1: Location, 2: Allele, 3: Gene,
            # 4: Feature, 5: Feature_type, 6: Consequence,
            # 7: cDNA_position, 8: CDS_position, 9: Protein_position,
            # 10: Amino_acids, 11: Codons, 12: Existing_variation, 13: Extra
            location = parts[1]
            allele = parts[2]
            transcript_id = parts[4]
            consequence = parts[6]
            cds_position = parts[8]
            protein_position = parts[9]
            amino_acids = parts[10]
            codons = parts[11]

            key = (location, allele, transcript_id)
            rows[key] = RustRow(
                location=location,
                allele=allele,
                transcript_id=transcript_id,
                consequence_set=consequence,
                amino_acids=amino_acids,
                codons=codons,
                cds_position=cds_position,
                protein_position=protein_position,
            )
    return rows


def normalize_consequences(cons: str) -> set[str]:
    """Normalize consequence string to a set of terms."""
    return set(cons.replace(",", "&").split("&")) - {""}


def compare(
    perl_rows: dict[tuple[str, str, str], PerlRow],
    rust_rows: dict[tuple[str, str, str], RustRow],
    discordant_only: bool,
    out: TextIO,
) -> None:
    """Compare Perl and Rust rows, reporting field-level divergences."""

    field_diverge_counts: Counter[str] = Counter()
    consequence_diverge_patterns: Counter[str] = Counter()
    total_compared = 0
    total_discordant = 0
    missing_in_rust = 0
    missing_in_perl = 0

    all_keys = set(perl_rows.keys()) | set(rust_rows.keys())

    print(f"Perl rows: {len(perl_rows)}", file=out)
    print(f"Rust rows: {len(rust_rows)}", file=out)
    print(f"Unique keys: {len(all_keys)}", file=out)
    print(file=out)

    for key in sorted(all_keys):
        perl = perl_rows.get(key)
        rust = rust_rows.get(key)

        if perl is None:
            missing_in_perl += 1
            continue
        if rust is None:
            missing_in_rust += 1
            continue

        total_compared += 1

        # Compare consequences
        perl_cons = normalize_consequences(perl.consequence_set)
        rust_cons = normalize_consequences(rust.consequence_set)
        cons_match = perl_cons == rust_cons

        # Compare fields
        divergent_fields: list[str] = []

        if not cons_match:
            divergent_fields.append("consequence")

        # Compare CDS position
        rust_cds = rust.cds_position.split("-")[0] if rust.cds_position != "-" else ""
        if perl.cds_start and rust_cds and perl.cds_start != rust_cds:
            divergent_fields.append("cds_start")

        # Compare codons
        if perl.ref_codon and rust.codons != "-":
            rust_codons = rust.codons.split("/")
            if len(rust_codons) == 2:
                if perl.ref_codon.upper() != rust_codons[0].upper():
                    divergent_fields.append("ref_codon")
                if perl.alt_codon.upper() != rust_codons[1].upper():
                    divergent_fields.append("alt_codon")

        # Compare peptides
        if perl.ref_peptide and rust.amino_acids != "-":
            rust_aas = rust.amino_acids.split("/")
            if len(rust_aas) == 2:
                if perl.ref_peptide != rust_aas[0]:
                    divergent_fields.append("ref_peptide")
                if perl.alt_peptide != rust_aas[1]:
                    divergent_fields.append("alt_peptide")

        if divergent_fields:
            total_discordant += 1
            for f in divergent_fields:
                field_diverge_counts[f] += 1

            if not cons_match:
                pattern = f"perl={','.join(sorted(perl_cons))} rust={','.join(sorted(rust_cons))}"
                consequence_diverge_patterns[pattern] += 1

        if not discordant_only or divergent_fields:
            loc, allele, tr = key
            status = "DISCORD" if divergent_fields else "match"
            fields_str = ",".join(divergent_fields) if divergent_fields else ""
            print(
                f"{status}\t{loc}\t{allele}\t{tr}\t{fields_str}\t"
                f"perl_cons={perl.consequence_set}\trust_cons={rust.consequence_set}\t"
                f"perl_cds={perl.cds_start}-{perl.cds_end}\t"
                f"perl_pep={perl.ref_peptide}/{perl.alt_peptide}\t"
                f"rust_aa={rust.amino_acids}\t"
                f"perl_codon={perl.ref_codon}/{perl.alt_codon}\t"
                f"rust_codon={rust.codons}",
                file=out,
            )

    # Summary
    print(file=out)
    print("=" * 60, file=out)
    print("SUMMARY", file=out)
    print("=" * 60, file=out)
    print(f"Total compared:    {total_compared}", file=out)
    print(f"Total discordant:  {total_discordant}", file=out)
    print(f"Missing in Rust:   {missing_in_rust}", file=out)
    print(f"Missing in Perl:   {missing_in_perl}", file=out)
    print(file=out)

    print("Field divergence counts:", file=out)
    for field, count in field_diverge_counts.most_common():
        pct = 100.0 * count / total_discordant if total_discordant else 0
        print(f"  {field:20s}  {count:6d}  ({pct:.1f}%)", file=out)

    print(file=out)
    print("Top consequence divergence patterns:", file=out)
    for pattern, count in consequence_diverge_patterns.most_common(20):
        print(f"  {count:6d}  {pattern}", file=out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--perl-trace", required=True, help="Perl tracer TSV output")
    parser.add_argument("--rust-output", required=True, help="Rust VEP default output")
    parser.add_argument(
        "--discordant-only",
        action="store_true",
        help="Only print discordant rows",
    )
    parser.add_argument("--output", default="-", help="Output file (default: stdout)")
    args = parser.parse_args()

    perl_rows = parse_perl_trace(Path(args.perl_trace))
    rust_rows = parse_rust_output(Path(args.rust_output))

    if args.output == "-":
        compare(perl_rows, rust_rows, args.discordant_only, sys.stdout)
    else:
        with open(args.output, "w") as f:
            compare(perl_rows, rust_rows, args.discordant_only, f)


if __name__ == "__main__":
    main()
