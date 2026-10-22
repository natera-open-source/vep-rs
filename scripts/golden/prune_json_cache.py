#!/usr/bin/env python3
"""Prune a vep-rs JSON transcript cache to the transcripts a set of VCF records can reach.

vep-rs annotates a variant against every transcript whose span, extended by the
up/downstream flank (5,000 bp on both sides), overlaps the variant. A fixture cache therefore
needs exactly those transcripts and nothing else. This script reads one or more VCFs, keeps
every transcript record of the source cache whose extended span overlaps any record span
(symbolic alleles use END or SVLEN; breakend mates count as one-base spans on their own
chromosome; a record's chromosome name is resolved to the cache's as the runtime resolves it,
through the cache's `chr_synonyms.txt` and the `chr` prefix), and writes the survivors into
shards named exactly as in the source so a transcript the source stores in two shards is kept
in both. Chromosome directories with no
surviving transcript are omitted. The summary printed at the end counts, per chromosome,
the transcripts kept and the spans searched: one span per record, one per breakend mate and
one per same-chromosome breakend pair covering both ends, so a chromosome's spans outnumber
its records whenever it carries breakends. `--gzip` writes each shard as
`<start>-<end>.json.gz`, and `--drop <key>` removes a `variation_effect_feature_cache` entry
vep-rs can do without; `peptide` is not one of them (the reference peptide's sequence edits
are read from it, so a cache without it prints a non-ATG start codon's residue as the codon
translates rather than as VEP's edited `M`) and is refused. `--drop-except <file>` names the
transcripts (one stable id per line) that keep every `--drop` key, so a corpus can carry the
protein features and prediction matrices of the transcripts its records test and no others.

`--variations <dir>` prunes the variation files beside the transcripts: `<dir>/<chr>/all_vars.gz`
(an Ensembl VEP cache directory, or a converted cache's `variations/`) is read for every
chromosome a record touches, and the records whose start lies within `--variation-flank`
(1 bp by default, the window vep-rs and Ensembl VEP query) of a record span are written as
`<dst>/variations/<chr>/all_vars.gz`, bgzip-compressed and indexed by a `.csi` on the start
column, the layout and index the Ensembl VEP dumper writes. The start column is read from
`variation_cols` of the `info.json` being written. `--exclude-variation <regex>` (repeatable)
drops every record whose name matches, so a slice can leave out identifiers that are not
redistributable; the summary counts the records kept and excluded per chromosome.
`--variation-info <info.json>` names the metadata of the cache the variation files belong to
when it is not the cache the rest of `info.json` describes (a native cache's `variations/`
beside an Ensembl cache's `info.txt`): its `variation_cols` and `var_type` are written in
place of the others', so the columns `info.json` declares are the columns the slice has.

`info.json` is the cache metadata vep-rs prints in its output headers. It is built from the
Ensembl VEP cache's `info.txt` when `--perl-info` names one (every `source_<key>` line
becomes a `source_versions` entry, `variation_cols` becomes a list, `cell_types` is dropped), else
copied from the source cache, else written from `--species`, `--assembly` and
`--cache-version`. `variation_cols` and `var_type` are kept only when the output carries a
variation slice, so a cache's metadata never advertises variation data its directory lacks, and
`contig_order` names the output's own chromosome directories in the converter's order. A source
cache's `chr_synonyms.txt` is copied beside `info.json`, kept to the lines that name a
chromosome the output holds. `--info-only` writes just that file, then refreshes the manifest.

The output ends with the checksum manifest `CHECKSUMS.sha256` the cache tools write: one
`sha256sum` line per file of the output directory, the manifest itself excepted, paths relative
to the directory in byte order, so `sha256sum --check CHECKSUMS.sha256` run in the directory
verifies the slice and vep-rs checks each file it reads against it. `--manifest-only` writes
that file alone for a cache directory assembled by other means.

Usage:
    prune_json_cache.py --cache <src-dir> --out <dst-dir> --vcf <a.vcf> [--vcf <b.vcf.gz>]
        [--flank 5000] [--gzip] [--drop sorted_exons ...] [--drop-except <ids.txt>]
        [--perl-info <cache>/info.txt]
        [--variations <dir> [--variation-flank 1] [--exclude-variation <regex> ...]
         [--variation-info <cache>/info.json]]
        [--species homo_sapiens --assembly GRCh37 --cache-version 116]
    prune_json_cache.py --info-only --out <dst-dir> --perl-info <cache>/info.txt --cache-version 116
    prune_json_cache.py --manifest-only --out <cache-dir>

The source layout is `<dir>/transcripts/<chr>/<start>-<end>.json[.gz]`, each shard a JSON
array of transcript objects with integer-like `start`/`end` strings.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import struct
import sys
import zlib
from bisect import bisect_right
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


def prune(cache: Path, out: Path, spans: dict[str, list[tuple[int, int]]], flank: int, compress: bool = False, drop: tuple[str, ...] = (),
          drop_except: frozenset[str] = frozenset()) -> dict:
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
                if isinstance(vefc, dict) and t.get("stable_id") not in drop_except:
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


# bgzip writes blocks of at most this many uncompressed bytes.
BGZF_BLOCK_BYTES = 0xFF00
# The 28-byte empty block bgzip ends every file with.
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
# The binning parameters of the `.csi` files in the Ensembl VEP variation caches.
CSI_MIN_SHIFT = 14
CSI_DEPTH = 8
# The distance, in compressed-file bytes between the blocks holding a bin's first and last
# chunk, below which the bin's chunks are folded into its parent bin.
MIN_MARKER_DIST = 0x10000
# The tabix header's `format` field for a generic 1-based file (`tabix` without `-0`).
TABIX_GENERIC_ONE_BASED = 0
VARIATION_CHR_COLUMN = "chr"
VARIATION_NAME_COLUMN = "variation_name"
VARIATION_START_COLUMN = "start"
VARIATION_COLUMNS_NEEDED = (VARIATION_CHR_COLUMN, VARIATION_NAME_COLUMN, VARIATION_START_COLUMN)


def bgzf_block(data: bytes) -> bytes:
    """One bgzf block: a gzip member with the `BC` extra field holding its own size."""
    deflater = zlib.compressobj(6, zlib.DEFLATED, -15)
    payload = deflater.compress(data) + deflater.flush()
    block_size = 18 + len(payload) + 8
    header = struct.pack("<4BIBBH2BHH", 0x1F, 0x8B, 8, 4, 0, 0, 0xFF, 6, ord("B"), ord("C"), 2, block_size - 1)
    return header + payload + struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))


class BgzfWriter:
    """Writes a bgzf file block by block and reports virtual offsets, as `bgzf_tell` does."""

    def __init__(self, handle):
        self.handle = handle
        self.buffer = bytearray()
        self.block_offset = 0

    def tell(self) -> int:
        """The virtual offset of the next byte: the block's file offset shifted by 16, plus the
        offset inside the block."""
        return (self.block_offset << 16) | len(self.buffer)

    def write(self, data: bytes) -> None:
        self.buffer += data
        while len(self.buffer) >= BGZF_BLOCK_BYTES:
            self._flush(BGZF_BLOCK_BYTES)

    def _flush(self, size: int) -> None:
        block = bgzf_block(bytes(self.buffer[:size]))
        self.handle.write(block)
        self.block_offset += len(block)
        del self.buffer[:size]

    def close(self) -> None:
        if self.buffer:
            self._flush(len(self.buffer))
        self.handle.write(BGZF_EOF)


def reg2bin(beg: int, end: int, min_shift: int = CSI_MIN_SHIFT, depth: int = CSI_DEPTH) -> int:
    """The bin of a 0-based half-open interval (`hts_reg2bin`)."""
    end -= 1
    level, shift, offset = depth, min_shift, ((1 << (depth * 3)) - 1) // 7
    while level > 0:
        if beg >> shift == end >> shift:
            return offset + (beg >> shift)
        level -= 1
        shift += 3
        offset -= 1 << (level * 3)
    return 0


def bin_first(level: int) -> int:
    return ((1 << (level * 3)) - 1) // 7


def bin_level(b: int) -> int:
    level = 0
    while b:
        b = (b - 1) >> 3
        level += 1
    return level


def bin_bottom(b: int, depth: int = CSI_DEPTH) -> int:
    """The first window, at the finest level, a bin covers (`hts_bin_bot`)."""
    level = bin_level(b)
    return (b - bin_first(level)) << ((depth - level) * 3)


class CsiBuilder:
    """Builds the CSI index `tabix -C -s <seq> -b <start> -e <start>` writes for one sequence
    per file: records keyed by their start alone, binned with `CSI_MIN_SHIFT` and `CSI_DEPTH`,
    the chunks of consecutive records of one bin merged, small bins folded into their parents,
    and each bin's `loffset` taken from the linear index the way htslib fills it."""

    def __init__(self, names: list[str], seq_col: int, start_col: int):
        self.names = names
        self.seq_col = seq_col
        self.start_col = start_col
        # Per reference: bin -> [[chunk_begin, chunk_end], ...], linear index, record count.
        self.bins: list[dict[int, list[list[int]]]] = [{} for _ in names]
        # None marks a window no record has reached yet; the first record's virtual offset
        # is 0, so 0 cannot.
        self.linear: list[list[int | None]] = [[] for _ in names]
        self.counts = [0] * len(names)
        self.first_offset = [None] * len(names)
        self.open_run: tuple[int, int, int] | None = None  # (tid, bin, run begin offset)

    def _insert_chunk(self, tid: int, b: int, begin: int, end: int) -> None:
        chunks = self.bins[tid].setdefault(b, [])
        if chunks and chunks[-1][1] >> 16 == begin >> 16:
            chunks[-1][1] = end
        else:
            chunks.append([begin, end])

    def add(self, tid: int, start: int, offset: int) -> None:
        """Adds a record at 1-based `start` whose bytes begin at virtual offset `offset`."""
        beg, end = start - 1, start
        b = reg2bin(beg, end)
        linear = self.linear[tid]
        window = beg >> CSI_MIN_SHIFT
        if window >= len(linear):
            linear.extend([None] * (window + 1 - len(linear)))
        if linear[window] is None:
            linear[window] = offset
        if self.open_run is None or self.open_run[:2] != (tid, b):
            if self.open_run is not None:
                run_tid, run_bin, run_begin = self.open_run
                self._insert_chunk(run_tid, run_bin, run_begin, offset)
            self.open_run = (tid, b, offset)
        if self.first_offset[tid] is None:
            self.first_offset[tid] = offset
        self.counts[tid] += 1

    def finish(self, final_offset: int) -> bytes:
        if self.open_run is not None:
            run_tid, run_bin, run_begin = self.open_run
            self._insert_chunk(run_tid, run_bin, run_begin, final_offset)
            self.open_run = None
        n_bins = bin_first(CSI_DEPTH + 1)
        meta_bin = n_bins + 1
        out = bytearray(b"CSI\x01")
        names_blob = b"".join(name.encode() + b"\0" for name in self.names)
        aux = struct.pack("<7i", TABIX_GENERIC_ONE_BASED, self.seq_col, self.start_col, self.start_col, ord("#"), 0, len(names_blob)) + names_blob
        out += struct.pack("<iii", CSI_MIN_SHIFT, CSI_DEPTH, len(aux)) + aux
        out += struct.pack("<i", len(self.names))
        for tid, bins in enumerate(self.bins):
            linear = self.linear[tid]
            for i, entry in enumerate(linear):
                if entry is None:
                    linear[i] = linear[i - 1] if i else self.first_offset[tid]
            loffset = {}
            for b in bins:
                bottom = bin_bottom(b)
                loffset[b] = linear[bottom] if bottom < len(linear) else 0
            # `compress_binning`: fold a bin whose chunks span less than MIN_MARKER_DIST of
            # compressed file into its parent, deepest level first, when the parent holds records;
            # then sort bin 0 and merge adjacent chunks that start in the same block.
            for level in range(CSI_DEPTH, 0, -1):
                first, next_first = bin_first(level), bin_first(level + 1)
                for b in sorted(k for k in bins if first <= k < next_first):
                    chunks = bins[b]
                    if level < CSI_DEPTH and len(chunks) > 1:
                        chunks.sort()
                    if (chunks[-1][1] >> 16) - (chunks[0][0] >> 16) < MIN_MARKER_DIST:
                        parent = (b - 1) >> 3
                        if parent in bins:
                            bins[parent].extend(chunks)
                            del bins[b]
            if 0 in bins:
                bins[0].sort()
            for chunks in bins.values():
                merged = [chunks[0]]
                for chunk in chunks[1:]:
                    if merged[-1][1] >> 16 >= chunk[0] >> 16:
                        merged[-1][1] = max(merged[-1][1], chunk[1])
                    else:
                        merged.append(chunk)
                chunks[:] = merged
            if self.counts[tid]:
                bins[meta_bin] = [[self.first_offset[tid], final_offset], [self.counts[tid], 0]]
            out += struct.pack("<i", len(bins))
            for b in sorted(bins):
                chunks = bins[b]
                out += struct.pack("<IQi", b, loffset.get(b, 0), len(chunks))
                for begin, end in chunks:
                    out += struct.pack("<QQ", begin, end)
        out += struct.pack("<Q", 0)
        return bytes(out)


def variation_windows(spans: list[tuple[int, int]], flank: int) -> list[tuple[int, int]]:
    """The sorted, merged start windows the spans of one chromosome query."""
    windows = sorted((max(s - flank, 1), e + flank) for s, e in spans)
    merged: list[tuple[int, int]] = []
    for low, high in windows:
        if merged and low <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], high))
        else:
            merged.append((low, high))
    return merged


def start_in_windows(start: int, windows: list[tuple[int, int]], lows: list[int]) -> bool:
    """Whether `start` lies in one of the sorted, disjoint `windows`, whose lower bounds are `lows`."""
    i = bisect_right(lows, start)
    return i > 0 and start <= windows[i - 1][1]


def prune_variations(src: Path, out: Path, spans: dict[str, list[tuple[int, int]]], flank: int,
                     columns: list[str], exclude: list[re.Pattern]) -> dict:
    """Writes `<out>/variations/<chr>/all_vars.gz` and its `.csi` for every chromosome of
    `spans` whose `<src>/<chr>/all_vars.gz` exists: the records whose start lies within `flank`
    of a span, in file order, minus those whose name matches a pattern of `exclude`. `columns`
    is the cache's `variation_cols`, which the caller has checked names every column of
    `VARIATION_COLUMNS_NEEDED`."""
    chr_col = columns.index(VARIATION_CHR_COLUMN)
    name_col = columns.index(VARIATION_NAME_COLUMN)
    start_col = columns.index(VARIATION_START_COLUMN)
    last_field = max(name_col, start_col)
    summary: dict[str, dict] = {}
    for chrom, chrom_spans in sorted(spans.items()):
        source = src / chrom / "all_vars.gz"
        if not source.is_file():
            continue
        windows = variation_windows(chrom_spans, flank)
        lows = [w[0] for w in windows]
        kept = excluded = 0
        dst_dir = out / "variations" / chrom
        dst_dir.mkdir(parents=True, exist_ok=True)
        dst = dst_dir / "all_vars.gz"
        index_dst = dst_dir / "all_vars.gz.csi"
        index = CsiBuilder([chrom], chr_col + 1, start_col + 1)
        with gzip.open(source, "rb") as reader, dst.open("wb") as raw:
            writer = BgzfWriter(raw)
            for line in reader:
                fields = line.split(b"\t", last_field + 1)
                start = int(fields[start_col])
                if start > windows[-1][1]:
                    break
                if not start_in_windows(start, windows, lows):
                    continue
                if any(pattern.search(fields[name_col].decode()) for pattern in exclude):
                    excluded += 1
                    continue
                index.add(0, start, writer.tell())
                writer.write(line)
                kept += 1
            writer.close()
            final_offset = writer.tell()
        if kept == 0:
            dst.unlink()
            index_dst.unlink(missing_ok=True)
            if not any(dst_dir.iterdir()):
                dst_dir.rmdir()
        else:
            with index_dst.open("wb") as raw:
                index_writer = BgzfWriter(raw)
                index_writer.write(index.finish(final_offset))
                index_writer.close()
        summary[chrom] = {"kept": kept, "excluded": excluded, "windows": len(windows)}
    variations_dir = out / "variations"
    if variations_dir.is_dir() and not any(variations_dir.iterdir()):
        variations_dir.rmdir()
    return summary


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
    if cache_version is not None:
        info["cache_version"] = cache_version
    elif isinstance(info.get("cache_version"), str):
        info["cache_version"] = int(info["cache_version"])
    return info


VARIATION_INFO_KEYS = ("variation_cols", "var_type")


def cache_info(args) -> dict | int:
    """The `info.json` object for the output: from `--perl-info`, else the source cache's own,
    else the command line; an exit status when the options do not name one."""
    if args.perl_info is not None:
        info = info_from_perl_info_txt(Path(args.perl_info).read_text(encoding="utf-8"), args.cache_version)
        if args.assembly is not None:
            info["assembly"] = args.assembly
        if info.get("species") is None:
            info["species"] = args.species
        if info.get("assembly") is None or info.get("cache_version") is None:
            print("ERROR: [prune_json_cache] --perl-info lacks assembly or cache version; pass --assembly/--cache-version", file=sys.stderr)
            return 2
        return info
    info_src = args.cache / "info.json" if args.cache is not None else None
    if info_src is not None and info_src.is_file():
        return json.loads(info_src.read_text(encoding="utf-8"))
    if args.assembly is None or args.cache_version is None:
        print("ERROR: [prune_json_cache] source has no info.json; pass --perl-info or --assembly and --cache-version", file=sys.stderr)
        return 2
    return {"species": args.species, "assembly": args.assembly, "cache_version": args.cache_version, "source_versions": {}}


def has_variation_files(out: Path) -> bool:
    return any((out / "variations").glob("*/all_vars.gz")) or any((out / "variations").glob("*/*_var.gz"))


def write_info(out: Path, info: dict) -> None:
    """Writes `info.json`, without the variation keys when the output holds no variation file,
    and with `contig_order` naming the output's own chromosome directories (autosomes by
    number, then X, Y and MT, then every other contig by name, the converter's order), or
    without the key when the output holds no transcript directory."""
    info = dict(info)
    if not has_variation_files(out):
        for key in VARIATION_INFO_KEYS:
            info.pop(key, None)

    def contig_key(name: str) -> tuple[int, int, str]:
        if name.isdigit():
            return (0, int(name), "")
        tier = {"X": 1, "Y": 2, "MT": 3}.get(name, 4)
        return (tier, 0, name)

    transcripts = out / "transcripts"
    chromosomes: list[str] = []
    if transcripts.is_dir():
        chromosomes = [p.name for p in transcripts.iterdir() if p.is_dir()]
    if chromosomes:
        info["contig_order"] = sorted(chromosomes, key=contig_key)
    else:
        info.pop("contig_order", None)
    (out / "info.json").write_text(json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")


MANIFEST_FILE = "CHECKSUMS.sha256"
# Directory levels the walk descends below the output directory before it treats the tree as a
# symbolic-link loop, the limit the Rust writer applies; a cache is three levels deep.
MANIFEST_MAX_DEPTH = 16


def write_manifest(out: Path) -> int:
    """Writes `<out>/CHECKSUMS.sha256`: one `<sha256 hex>  <path relative to out>` line per regular
    file under `out`, through symbolic links to files and to directories, the manifest itself
    excepted, in byte order of the paths with `/` between components, the shape `sha256sum --check`
    reads and the cache tools write. Returns the number of files listed; exits with the tool's
    error when directories nest more than `MANIFEST_MAX_DEPTH` levels deep or a symbolic link
    points to nothing, as the cache tools' writer fails on both."""
    entries: list[tuple[str, str]] = []
    for root, _dirnames, filenames in os.walk(out, followlinks=True):
        depth = len(Path(root).relative_to(out).parts)
        if depth > MANIFEST_MAX_DEPTH:
            raise SystemExit(f"ERROR: [prune_json_cache] {root}: directories nest more than {MANIFEST_MAX_DEPTH} levels deep")
        for name in filenames:
            path = Path(root) / name
            if not path.is_file():
                if path.is_symlink():
                    raise SystemExit(f"ERROR: [prune_json_cache] {path}: a symbolic link to nothing")
                continue
            relative = path.relative_to(out).as_posix()
            if relative == MANIFEST_FILE:
                continue
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1 << 20), b""):
                    digest.update(chunk)
            entries.append((relative, digest.hexdigest()))
    entries.sort(key=lambda entry: entry[0].encode("utf-8"))
    (out / MANIFEST_FILE).write_text("".join(f"{digest}  {relative}\n" for relative, digest in entries), encoding="utf-8")
    return len(entries)


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
    ap.add_argument("--drop-except", type=Path, default=None, metavar="FILE",
                    help="a file of transcript stable ids, one per line, that keep every --drop key")
    ap.add_argument("--perl-info", type=Path, default=None, help="Ensembl VEP cache info.txt to derive info.json from")
    ap.add_argument("--variations", type=Path, default=None,
                    help="a directory holding <chr>/all_vars.gz: its records within --variation-flank of a record span are "
                         "written as a bgzip and csi-indexed slice under variations/")
    ap.add_argument("--variation-flank", type=int, default=1,
                    help="how far from a record span a variation record's start may lie and still be kept (default 1)")
    ap.add_argument("--exclude-variation", action="append", default=[], metavar="REGEX",
                    help="drop every variation record whose name matches (repeatable)")
    ap.add_argument("--variation-info", type=Path, default=None, metavar="INFO_JSON",
                    help="the info.json of the cache the --variations files belong to, whose variation_cols and "
                         "var_type replace those of the info being written")
    ap.add_argument("--info-only", action="store_true",
                    help="write info.json, refresh CHECKSUMS.sha256 and nothing else")
    ap.add_argument("--manifest-only", action="store_true",
                    help="write CHECKSUMS.sha256 over the files already under --out and nothing else")
    ap.add_argument("--species", default="homo_sapiens")
    ap.add_argument("--assembly", default=None)
    ap.add_argument("--cache-version", type=int, default=None)
    args = ap.parse_args(argv)
    if "peptide" in args.drop:
        print("ERROR: [prune_json_cache] --drop peptide: vep-rs reads the reference peptide's sequence edits from it", file=sys.stderr)
        return 2
    if args.manifest_only:
        if not (args.out / "info.json").is_file():
            print(f"ERROR: [prune_json_cache] --manifest-only: {args.out} holds no info.json, so it is not a cache directory", file=sys.stderr)
            return 2
        print(json.dumps({"manifest_files": write_manifest(args.out)}))
        return 0
    args.out.mkdir(parents=True, exist_ok=True)
    info = cache_info(args)
    if isinstance(info, int):
        return info
    if args.variation_info is not None and args.variations is None:
        print("ERROR: [prune_json_cache] --variation-info describes the files of --variations; pass both", file=sys.stderr)
        return 2
    if args.info_only:
        if args.variations is not None:
            print("ERROR: [prune_json_cache] --info-only writes info.json and refreshes CHECKSUMS.sha256, nothing else; --variations and --variation-info describe a slice this run does not cut", file=sys.stderr)
            return 2
        write_info(args.out, info)
        write_manifest(args.out)
        return 0
    if args.cache is None or not (args.cache / "transcripts").is_dir():
        print(f"ERROR: [prune_json_cache] {args.cache} has no transcripts/ directory", file=sys.stderr)
        return 2
    if not args.vcf:
        print("ERROR: [prune_json_cache] at least one --vcf is required", file=sys.stderr)
        return 2
    columns: list[str] = []
    exclude: list[re.Pattern] = []
    if args.variations is not None:
        if args.variation_flank < 0:
            print(f"ERROR: [prune_json_cache] --variation-flank {args.variation_flank}: a distance is not negative", file=sys.stderr)
            return 2
        if args.variation_info is not None:
            try:
                variation_info = json.loads(args.variation_info.read_text(encoding="utf-8"))
            except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
                print(f"ERROR: [prune_json_cache] --variation-info {args.variation_info}: {getattr(error, 'strerror', None) or error}", file=sys.stderr)
                return 2
            if not isinstance(variation_info, dict) or not isinstance(variation_info.get("variation_cols"), list):
                print(f"ERROR: [prune_json_cache] --variation-info {args.variation_info} declares no variation_cols list", file=sys.stderr)
                return 2
            for key in VARIATION_INFO_KEYS:
                info.pop(key, None)
                if key in variation_info:
                    info[key] = variation_info[key]
        columns = info.get("variation_cols") or []
        if not isinstance(columns, list):
            print(f"ERROR: [prune_json_cache] --variations needs variation_cols as a list of column names in the info being written, not {columns!r}", file=sys.stderr)
            return 2
        missing = [name for name in VARIATION_COLUMNS_NEEDED if name not in columns]
        if missing:
            print(f"ERROR: [prune_json_cache] --variations needs variation_cols naming {', '.join(missing)} in the info being written (--perl-info or the source info.json)", file=sys.stderr)
            return 2
        try:
            exclude = [re.compile(pattern) for pattern in args.exclude_variation]
        except re.error as error:
            print(f"ERROR: [prune_json_cache] --exclude-variation {error.pattern!r}: {error}", file=sys.stderr)
            return 2
    cache_chromosomes = {p.name for p in (args.cache / "transcripts").iterdir() if p.is_dir()}
    synonyms_src = args.cache / "chr_synonyms.txt"
    synonym_lines: list[str] = []
    if synonyms_src.is_file():
        synonym_lines = synonyms_src.read_text(encoding="utf-8").splitlines()
    linked: dict[str, list[str]] = defaultdict(list)
    for line in synonym_lines:
        names = line.split()
        for synonym in names[1:]:
            linked[names[0]].append(synonym)
            linked[synonym].append(names[0])

    def cache_name(chrom: str) -> str:
        """The cache's own name for a record's chromosome, resolved as the runtime resolves an
        input name (Ensembl VEP's `get_source_chr_name`): the name itself, else a linked synonym
        the cache holds, else the name with `chr` removed or added, else a bare `M` as `MT`."""
        if chrom in cache_chromosomes:
            return chrom
        for synonym in linked.get(chrom, []):
            if synonym in cache_chromosomes:
                return synonym
        if chrom[:3].lower() == "chr":
            stripped = chrom[3:]
            if stripped in cache_chromosomes:
                return stripped
            return chrom
        if f"chr{chrom}" in cache_chromosomes:
            return f"chr{chrom}"
        if chrom == "M" and "MT" in cache_chromosomes:
            return "MT"
        return chrom

    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for vcf in args.vcf:
        for chrom, s in record_spans(vcf).items():
            spans[cache_name(chrom)].extend(s)
    # Checked before `prune()` writes the transcript shards, so a refused run leaves no partial output.
    if args.variations is not None and not any((args.variations / chrom / "all_vars.gz").is_file() for chrom in spans):
        print(f"ERROR: [prune_json_cache] --variations {args.variations}: no <chr>/all_vars.gz for any chromosome of the records ({', '.join(sorted(spans))})", file=sys.stderr)
        return 2
    if args.variations is not None:
        # The runtime zips variation_cols and a record's fields by position and reads a shorter
        # record silently, so a wrong declaration mislabels every value.
        for chrom in sorted(spans):
            source = args.variations / chrom / "all_vars.gz"
            if not source.is_file():
                continue
            try:
                with gzip.open(source, "rb") as reader:
                    first_record = reader.readline()
            except (OSError, EOFError, zlib.error) as error:
                print(f"ERROR: [prune_json_cache] --variations {source}: {getattr(error, 'strerror', None) or error}", file=sys.stderr)
                return 2
            if not first_record:
                continue
            field_count = first_record.count(b"\t") + 1
            if field_count != len(columns):
                print(f"ERROR: [prune_json_cache] --variations {source} carries {field_count} fields where variation_cols names {len(columns)} columns", file=sys.stderr)
                return 2
    drop_except: frozenset[str] = frozenset()
    if args.drop_except is not None:
        if not args.drop:
            print("ERROR: [prune_json_cache] --drop-except needs at least one --drop key to keep", file=sys.stderr)
            return 2
        try:
            id_lines = args.drop_except.read_text(encoding="utf-8").splitlines()
        except (OSError, UnicodeDecodeError) as error:
            print(f"ERROR: [prune_json_cache] --drop-except {args.drop_except}: {getattr(error, 'strerror', None) or error}", file=sys.stderr)
            return 2
        drop_except = frozenset(line.strip() for line in id_lines if line.strip())
    summary = prune(args.cache, args.out, spans, args.flank, args.gzip, tuple(args.drop), drop_except)
    if drop_except:
        summary["drop_except"] = len(drop_except)
    if args.variations is not None:
        summary["variations_per_chromosome"] = prune_variations(
            args.variations, args.out, spans, args.variation_flank, columns, exclude)
    # The cache's chr_synonyms.txt, kept to the lines naming a chromosome the output holds,
    # so the slice resolves the same input names for its chromosomes as the whole cache does.
    if synonyms_src.is_file():
        kept_chromosomes = set(summary["transcripts_per_chromosome"])
        kept_lines = [line for line in synonym_lines if any(name in kept_chromosomes for name in line.split())]
        (args.out / "chr_synonyms.txt").write_text("".join(f"{line}\n" for line in kept_lines), encoding="utf-8")
        summary["synonym_lines"] = len(kept_lines)
    write_info(args.out, info)
    summary["manifest_files"] = write_manifest(args.out)
    summary["spans_by_chromosome"] = {c: len(v) for c, v in sorted(spans.items())}
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
