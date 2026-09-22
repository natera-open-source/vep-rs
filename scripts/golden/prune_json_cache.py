#!/usr/bin/env python3
"""Prune a vep-rs JSON transcript cache to the transcripts a set of VCF records can reach.

vep-rs annotates a variant against every transcript whose span, extended by the
up/downstream flank (5,000 bp on both sides), overlaps the variant. A fixture cache therefore
needs exactly those transcripts and nothing else. This script reads one or more VCFs, keeps
every transcript record of the source cache whose extended span overlaps any record span
(symbolic alleles use END or SVLEN; breakend mates count as one-base spans on their own
chromosome), and writes the survivors into shards named exactly as in the source so a
transcript the source stores in two shards is kept in both. Chromosome directories with no
surviving transcript are omitted. `--gzip` writes each shard as `<start>-<end>.json.gz`, and
`--drop <key>` removes a `variation_effect_feature_cache` entry vep-rs can do without;
`peptide` is not one of them (the reference peptide's sequence edits are read from it, so a
cache without it prints a non-ATG start codon's residue as the codon translates rather than
as VEP's edited `M`) and is refused.

`info.json` is the cache metadata vep-rs prints in its output headers. It is built from the
Ensembl VEP cache's `info.txt` when `--perl-info` names one (every `source_<key>` line
becomes a `source_versions` entry, `variation_cols` becomes a list, `cell_types` is dropped), else
copied from the source cache, else written from `--species`, `--assembly` and
`--cache-version`. `--info-only` writes just that file.

Usage:
    prune_json_cache.py --cache <src-dir> --out <dst-dir> --vcf <a.vcf> [--vcf <b.vcf.gz>]
        [--flank 5000] [--gzip] [--drop sorted_exons ...] [--perl-info <cache>/info.txt]
        [--species homo_sapiens --assembly GRCh37 --cache-version 115]
    prune_json_cache.py --info-only --out <dst-dir> --perl-info <cache>/info.txt --cache-version 115

The source layout is `<dir>/transcripts/<chr>/<start>-<end>.json[.gz]`, each shard a JSON
array of transcript objects with integer-like `start`/`end` strings.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

BND_MATE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def opener(path: Path):
    return gzip.open if str(path).endswith(".gz") else open


def record_spans(vcf: Path) -> dict[str, list[tuple[int, int]]]:
    """Per chromosome, the [start, end] genomic spans a record and its breakend mates cover."""
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with opener(vcf)(vcf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, info = f[0], int(f[1]), f[3], f[4], f[7] if len(f) > 7 else ""
            end = pos + max(len(ref) - 1, 0)
            for kv in info.split(";"):
                if kv.startswith("END="):
                    try:
                        end = max(end, int(kv[4:]))
                    except ValueError:
                        pass
                elif kv.startswith("SVLEN="):
                    try:
                        end = max(end, pos + abs(int(kv[6:].split(",")[0])))
                    except ValueError:
                        pass
            spans[chrom].append((pos, end))
            # Breakend mates: bracket notation in ALT, or the CHR2/POS2/END2 INFO
            # keys of a symbolic <BND>. A mate on the same chromosome is annotated
            # over the whole breakend-pair span, so that span is kept too.
            mates: list[tuple[str, int, int]] = []
            for m in BND_MATE.finditer(alt):
                mates.append((m.group(1).removeprefix("chr"), int(m.group(2)), int(m.group(2))))
            kv = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
            if "CHR2" in kv and ("POS2" in kv or "END2" in kv):
                try:
                    p2 = int(kv.get("POS2", kv.get("END2")))
                    e2 = int(kv.get("END2", kv.get("POS2")))
                    mates.append((kv["CHR2"].removeprefix("chr"), min(p2, e2), max(p2, e2)))
                except ValueError:
                    pass
            for mate_chr, m_start, m_end in mates:
                spans[mate_chr].append((m_start, m_end))
                if mate_chr == chrom:
                    spans[chrom].append((min(pos, m_start), max(end, m_end)))
    return spans


def overlaps(t_start: int, t_end: int, spans: list[tuple[int, int]], flank: int) -> bool:
    lo, hi = t_start - flank, t_end + flank
    return any(s <= hi and e >= lo for s, e in spans)


def read_shard(shard: Path) -> list:
    with opener(shard)(shard, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def prune(cache: Path, out: Path, spans: dict[str, list[tuple[int, int]]], flank: int, compress: bool = False, drop: tuple[str, ...] = ()) -> dict:
    kept: dict[str, int] = {}
    shards_written = 0
    for chrom, chrom_spans in sorted(spans.items()):
        src_dir = cache / "transcripts" / chrom
        if not src_dir.is_dir():
            continue
        chrom_spans = sorted(chrom_spans)
        shards = sorted(list(src_dir.glob("*.json")) + list(src_dir.glob("*.json.gz")))
        for shard in shards:
            records = read_shard(shard)
            survivors = [t for t in records if overlaps(int(t["start"]), int(t["end"]), chrom_spans, flank)]
            if not survivors:
                continue
            for t in survivors:
                vefc = t.get("variation_effect_feature_cache")
                if isinstance(vefc, dict):
                    for key in drop:
                        vefc.pop(key, None)
            name = shard.name[: -len(".gz")] if shard.name.endswith(".gz") else shard.name
            dst = out / "transcripts" / chrom / (name + ".gz" if compress else name)
            dst.parent.mkdir(parents=True, exist_ok=True)
            payload = json.dumps(survivors, separators=(",", ":"), ensure_ascii=False) + "\n"
            if compress:
                # mtime 0 keeps the archive byte-stable across regenerations.
                with gzip.GzipFile(dst, "wb", compresslevel=9, mtime=0) as fh:
                    fh.write(payload.encode("utf-8"))
            else:
                dst.write_text(payload, encoding="utf-8")
            kept[chrom] = kept.get(chrom, 0) + len(survivors)
            shards_written += 1
    return {"transcripts_per_chromosome": kept, "shards": shards_written}


LIST_KEYS = ("variation_cols",)
FLAG_KEYS = ("regulatory",)
# Cell-type names contain commas, so the line cannot be split into a list; VEP stores it as
# one string and vep-rs does not read it.
SKIP_KEYS = ("cell_types",)


def info_from_perl_info_txt(text: str, cache_version: int | None) -> dict:
    """Translate an Ensembl VEP cache `info.txt` into vep-rs's `info.json` object.

    Follows VEP's own reader: `source_<key>` lines are version data, `variation_cols` is a
    comma-separated list, a value of `-` means unset, and every other line is a scalar.
    """
    info: dict = {"species": None, "assembly": None, "cache_version": cache_version, "source_versions": {}}
    for line in text.splitlines():
        if not line or line.startswith("#") or "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        if key in SKIP_KEYS:
            continue
        if key.startswith("source_"):
            info["source_versions"][key[len("source_"):]] = value
        elif key in LIST_KEYS:
            info[key] = [v for v in value.split(",") if v]
        elif key in FLAG_KEYS:
            info[key] = value not in ("", "0", "-")
        elif value != "-":
            info[key] = value
    if info.get("cache_version") is None and "cache_version" in info:
        info["cache_version"] = int(info["cache_version"])
    return info


def write_info(out: Path, args) -> int:
    dest = out / "info.json"
    if args.perl_info is not None:
        info = info_from_perl_info_txt(Path(args.perl_info).read_text(encoding="utf-8"), args.cache_version)
        if args.assembly is not None:
            info["assembly"] = args.assembly
        if info.get("species") is None:
            info["species"] = args.species
        if info.get("assembly") is None or info.get("cache_version") is None:
            print("ERROR: [prune_json_cache] --perl-info lacks assembly or cache version; pass --assembly/--cache-version", file=sys.stderr)
            return 2
        dest.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return 0
    info_src = args.cache / "info.json" if args.cache is not None else None
    if info_src is not None and info_src.is_file():
        dest.write_text(info_src.read_text(encoding="utf-8"), encoding="utf-8")
        return 0
    if args.assembly is None or args.cache_version is None:
        print("ERROR: [prune_json_cache] source has no info.json; pass --perl-info or --assembly and --cache-version", file=sys.stderr)
        return 2
    info = {"species": args.species, "assembly": args.assembly, "cache_version": args.cache_version, "source_versions": {}}
    dest.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--vcf", action="append", type=Path, default=[])
    ap.add_argument("--flank", type=int, default=5000)
    ap.add_argument("--gzip", action="store_true", help="write shards as .json.gz")
    ap.add_argument("--drop", action="append", default=[], metavar="KEY",
                    help="omit this variation_effect_feature_cache key from every transcript; "
                         "vep-rs rebuilds sorted_exons from exons, and protein_features and "
                         "protein_function_predictions serve only --domains, --sift and --polyphen; "
                         "peptide is refused")
    ap.add_argument("--perl-info", type=Path, default=None, help="Ensembl VEP cache info.txt to derive info.json from")
    ap.add_argument("--info-only", action="store_true", help="write info.json and nothing else")
    ap.add_argument("--species", default="homo_sapiens")
    ap.add_argument("--assembly", default=None)
    ap.add_argument("--cache-version", type=int, default=None)
    args = ap.parse_args(argv)
    if "peptide" in args.drop:
        print("ERROR: [prune_json_cache] --drop peptide: vep-rs reads the reference peptide's sequence edits from it", file=sys.stderr)
        return 2
    args.out.mkdir(parents=True, exist_ok=True)
    if args.info_only:
        return write_info(args.out, args)
    if args.cache is None or not (args.cache / "transcripts").is_dir():
        print(f"ERROR: [prune_json_cache] {args.cache} has no transcripts/ directory", file=sys.stderr)
        return 2
    if not args.vcf:
        print("ERROR: [prune_json_cache] at least one --vcf is required", file=sys.stderr)
        return 2
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for vcf in args.vcf:
        for chrom, s in record_spans(vcf).items():
            spans[chrom].extend(s)
    summary = prune(args.cache, args.out, spans, args.flank, args.gzip, tuple(args.drop))
    rc = write_info(args.out, args)
    if rc:
        return rc
    summary["records_by_chromosome"] = {c: len(v) for c, v in sorted(spans.items())}
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
