#!/usr/bin/env python3
"""Prepare benchmark VCFs for concordance runs.

This script:
1) Decompresses .vcf.gz (or copies .vcf) files.
2) Optionally normalizes chromosome naming (chr1 <-> 1).
3) Filters to canonical contigs: one universal Ensembl-style
   (1-22, X, Y, MT) for both GRCh37 and GRCh38. UCSC-style chr-prefix
   inputs are auto-stripped at filter time so all engines see the same
   Ensembl-style canonical VCF.
4) Optionally truncates to N variants for smoke tests.
5) Writes a manifest with provenance, counts, per-file sha256, and per-file
   chr_stripped flag for audit.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


# Canonical contig set:
# universal Ensembl-style for both GRCh37 and GRCh38. The assembly difference
# is in the genome-sequence content, NOT the contig naming. UCSC-style chr-prefix
# inputs are auto-stripped at filter time so the membership check always operates
# on Ensembl-style names. M and MT are equivalent mitochondrial names; both
# included for compatibility.
#
# This single constant must match scripts/concordance/run_concordance.sh
# CANONICAL_CONTIGS_LIST and scripts/concordance/compare_vep_outputs.py
# CANONICAL_CONTIGS frozenset.
CANONICAL_CONTIGS = frozenset([str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"])


@dataclass
class FileSummary:
    input_file: str
    output_file: str
    detected_input_chrom_style: str
    detected_cache_chrom_style: str
    normalization_mode: str
    headers_written: int
    variants_written: int
    variants_dropped_non_canonical: int = 0
    output_sha256: str = ""
    chr_stripped: bool = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare benchmark VCF inputs")
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing benchmark .vcf.gz or .vcf files",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write prepared .vcf files",
    )
    parser.add_argument(
        "--chrom-mode",
        choices=["auto", "keep", "strip", "add"],
        default="auto",
        help="Chromosome normalization mode (default: auto). Independent of "
        "the canonical-contigs filter; if --canonical-contigs is active, "
        "UCSC-style inputs are auto-stripped regardless of this mode.",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Converted JSON cache directory (used only by --chrom-mode auto)",
    )
    parser.add_argument(
        "--max-variants",
        type=int,
        default=None,
        help="If set, keep only first N variant lines per input file",
    )
    parser.add_argument(
        "--assembly",
        choices=["GRCh37", "GRCh38"],
        default=None,
        help="Reference assembly metadata (recorded in manifest). Optional under "
        "one set, since the canonical contigs are the same for both assemblies. "
        "Pass it for audit-trail clarity.",
    )
    parser.add_argument(
        "--canonical-contigs",
        action="store_true",
        default=True,
        help="Filter input VCFs to canonical contigs (DEFAULT). "
        "Universal Ensembl-style (1-22,X,Y,MT) for both assemblies. UCSC-style "
        "chr-prefix inputs are auto-stripped at filter time. Drops variants on "
        "alt contigs (NT_*, NW_*, decoys), so that asymmetric alt-contig handling "
        "between engines cannot inflate F1 divergences.",
    )
    parser.add_argument(
        "--no-canonical-contigs",
        dest="canonical_contigs",
        action="store_false",
        help="Opt out of canonical-contigs filtering and annotate every contig "
        "the input carries, alt contigs included.",
    )
    return parser.parse_args()


def list_input_vcfs(input_dir: Path) -> list[Path]:
    vcfs = sorted(input_dir.glob("*.vcf")) + sorted(input_dir.glob("*.vcf.gz"))
    if not vcfs:
        raise FileNotFoundError(f"No .vcf or .vcf.gz files found in {input_dir}")
    return vcfs


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def first_variant_chrom(vcf_path: Path) -> Optional[str]:
    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0]:
                return fields[0]
    return None


def detect_style(contig: Optional[str]) -> str:
    if not contig:
        return "unknown"
    return "chr" if contig.lower().startswith("chr") else "plain"


def detect_cache_style(cache_dir: Optional[Path]) -> str:
    if cache_dir is None:
        return "unknown"

    for sub in ("transcripts", "variations"):
        root = cache_dir / sub
        if not root.exists() or not root.is_dir():
            continue
        dirs = sorted(p for p in root.iterdir() if p.is_dir())
        if not dirs:
            continue
        return detect_style(dirs[0].name)
    return "unknown"


def resolve_mode(requested: str, input_style: str, cache_style: str) -> str:
    if requested != "auto":
        return requested

    if input_style == "unknown" or cache_style == "unknown":
        return "keep"

    if input_style == cache_style:
        return "keep"

    if input_style == "chr" and cache_style == "plain":
        return "strip"
    if input_style == "plain" and cache_style == "chr":
        return "add"
    return "keep"


def normalize_chrom(chrom: str, mode: str) -> str:
    if mode == "strip":
        if chrom.lower().startswith("chr"):
            return chrom[3:]
        return chrom
    if mode == "add":
        if chrom.lower().startswith("chr"):
            return chrom
        return f"chr{chrom}"
    return chrom


def strip_chr(chrom: str) -> str:
    """Strip UCSC chr- prefix to get Ensembl-style contig name.

    chr1 -> 1, chrX -> X, chrM -> M (M and MT both retained as canonical).
    Plain Ensembl-style inputs are returned unchanged.
    """
    if chrom.lower().startswith("chr"):
        return chrom[3:]
    return chrom


def prepared_name(path: Path) -> str:
    if path.name.endswith(".vcf.gz"):
        return path.name[: -len(".gz")]
    if path.name.endswith(".vcf"):
        return path.name
    return f"{path.name}.vcf"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def prepare_single_file(
    src: Path,
    dst: Path,
    mode: str,
    max_variants: Optional[int],
    canonical_contigs: Optional[frozenset[str]],
) -> tuple[int, int, int, bool]:
    """Prepare one VCF: decompress, normalize chrom naming, optionally filter to
    canonical contigs (universal Ensembl-style; auto-strip UCSC chr-).

    Returns (headers_written, variants_written, variants_dropped, chr_stripped).
    """
    header_count = 0
    variant_count = 0
    dropped_count = 0
    chr_stripped = False
    dst.parent.mkdir(parents=True, exist_ok=True)

    with open_text(src) as reader, dst.open("w", encoding="utf-8") as writer:
        for line in reader:
            if line.startswith("#"):
                writer.write(line)
                header_count += 1
                continue

            if max_variants is not None and variant_count >= max_variants:
                break

            fields = line.rstrip("\n").split("\t")
            if fields:
                fields[0] = normalize_chrom(fields[0], mode)
                # Auto-strip the UCSC chr- prefix BEFORE the canonical
                # membership check so the filter operates on Ensembl-style names.
                # Track per-file chr_stripped flag for the manifest audit trail.
                if canonical_contigs is not None and fields[0].lower().startswith(
                    "chr"
                ):
                    fields[0] = strip_chr(fields[0])
                    chr_stripped = True
            if (
                canonical_contigs is not None
                and fields
                and fields[0] not in canonical_contigs
            ):
                dropped_count += 1
                continue
            writer.write("\t".join(fields) + "\n")
            variant_count += 1

    return header_count, variant_count, dropped_count, chr_stripped


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    cache_dir = Path(args.cache_dir).resolve() if args.cache_dir else None

    # One universal Ensembl-style canonical set for both assemblies. The
    # frozenset does not depend on --assembly; --assembly is recorded in the
    # manifest for audit only.
    canonical_contigs: Optional[frozenset[str]] = (
        CANONICAL_CONTIGS if args.canonical_contigs else None
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    vcfs = list_input_vcfs(input_dir)
    cache_style = detect_cache_style(cache_dir)
    summaries: list[FileSummary] = []

    for src in vcfs:
        input_style = detect_style(first_variant_chrom(src))
        mode = resolve_mode(args.chrom_mode, input_style, cache_style)
        dst = output_dir / prepared_name(src)
        headers, variants, dropped, chr_stripped = prepare_single_file(
            src, dst, mode, args.max_variants, canonical_contigs
        )
        out_sha = sha256_file(dst)

        summaries.append(
            FileSummary(
                input_file=str(src),
                output_file=str(dst),
                detected_input_chrom_style=input_style,
                detected_cache_chrom_style=cache_style,
                normalization_mode=mode,
                headers_written=headers,
                variants_written=variants,
                variants_dropped_non_canonical=dropped,
                output_sha256=out_sha,
                chr_stripped=chr_stripped,
            )
        )
        canonical_note = ""
        if canonical_contigs is not None:
            stripped_note = " chr_stripped=true" if chr_stripped else ""
            canonical_note = f", canonical_dropped={dropped}{stripped_note}"
        print(
            f"prepared {src.name} -> {dst.name} "
            f"(mode={mode}, variants={variants}, headers={headers}{canonical_note})"
        )

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "requested_chrom_mode": args.chrom_mode,
        "cache_dir": str(cache_dir) if cache_dir else None,
        "max_variants": args.max_variants,
        "assembly": args.assembly,
        "canonical_contigs_policy": "universal_ensembl_style"
        if args.canonical_contigs
        else "disabled",
        "canonical_contigs_filter_active": args.canonical_contigs,
        "canonical_contigs_set": (
            sorted(canonical_contigs) if canonical_contigs is not None else None
        ),
        "files": [asdict(summary) for summary in summaries],
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"wrote manifest: {manifest_path}")


if __name__ == "__main__":
    main()
