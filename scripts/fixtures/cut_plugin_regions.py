#!/usr/bin/env python3
"""Cut a plugin's open data files to the regions a golden corpus reads.

A corpus that compares a built-in plugin with its Ensembl VEP module on real,
redistributable data (LoFTEE's ancestral-allele FASTA and GERP bigWig) cannot
commit the whole files, so it commits a cut: the same contigs at their full
length, so every coordinate stays genomic, with the bases a run never reads
masked to `N` or left out of the bigWig. The regions come from the corpus itself:
the transcripts of its `json_cache` (spans and exons) and the records of its
`variants.vcf`. A record on a sequence the cache holds no transcript of keeps
its own region of that sequence's contig: the contig is looked up under the
record's name in the reference's index, else in the ancestral FASTA's index,
else among the contigs the bigWig declares, and is written under the name it
has there.

    reference.fa.gz (+ .fa.gz.fai, .fa.gz.gzi, .fa.fai)
        the run's reference, `--reference` with every base outside each transcript
        span plus `--flank` and each record plus `--record-flank` masked to N; the
        manifest names it under `fasta`
    plugin_data/human_ancestor.fa.gz (+ .fa.gz.fai, .fa.gz.gzi, .fa.fai)
        `--ancestor` masked the same way around the records alone: the plugin
        reads the ancestral base at the variant and nowhere else
    plugin_data/gerp.bw (or `--gerp-name`)
        `--gerp` restricted to each transcript's exons plus `--exon-flank`, the
        bases LoFTEE's END_TRUNC walk sums (`--gerp-regions exons`, the default),
        or to each record plus `--record-flank`, the span a per-variant plugin
        such as Conservation reads (`--gerp-regions records`); written without
        zoom levels so a summary read of the file totals the per-base values
        exactly, as an interval read does. `--gerp-name` sets the file name a
        plugin derives from its directory, LoFTEE's
        `gerp_conservation_scores.homo_sapiens.GRCh38.bw`

Each of the three source files is optional and only the ones given are cut; a
corpus whose plugin reads the bigWig alone passes `--gerp` alone. Each masked
FASTA keeps the contig lengths of its own index. The cut bigWig declares the
contig lengths of the reference's index, else of the ancestral FASTA's index,
else the lengths the source bigWig declares. A cut of the reference alone, the
cut of a corpus that reads no plugin data, writes no `plugin_data` directory and
runs without `bigtools`.

The FASTA indexes are written in samtools' format (name, length, offset, bases
per line, bytes per line) for the plain and the bgzipped file, and `bgzip -i`
writes the `.gzi`. The bigWig is cut with `bigtools bigwigtobedgraph` over the
region list and rebuilt with `bigtools bedgraphtobigwig --nzooms 0`. The versions
the tools of a cut report (`bgzip` where a FASTA is given, `bigtools` where the
bigWig is) and the sha256 of every file written go under `plugin_fixtures`
in the corpus's `provenance.json`, beside the entries other generators record
there under a file's name; those entries and every other key of the record are
kept. An existing fixture file is overwritten only under `--force`. A run prints
one JSON line: the transcripts and records it read, the bases the masked
reference keeps (`reference_bases_kept`) and the bases of the regions the bigWig
is cut to (`gerp_bases`), each 0 for a file the run does not cut, and the files
written.

Usage:
    cut_plugin_regions.py --corpus tests/golden/116/<corpus> [--reference <genome.fa>]
        [--ancestor <human_ancestor.fa>] [--gerp <gerp.bw>] [--gerp-regions exons|records]
        [--gerp-name gerp.bw] [--flank 5000] [--record-flank 2000] [--exon-flank 50] [--force]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

LINE_BASES = 60


def fail(message: str) -> SystemExit:
    return SystemExit(f"ERROR: [cut_plugin_regions] {message}")


def read_fai(path: Path) -> dict[str, tuple[int, int, int, int]]:
    """`name -> (length, offset, bases per line, bytes per line)` of a samtools index."""
    index: dict[str, tuple[int, int, int, int]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        name, length, offset, line_bases, line_bytes = line.split("\t")[:5]
        index[name] = (int(length), int(offset), int(line_bases), int(line_bytes))
    return index


def fetch(fasta: Path, entry: tuple[int, int, int, int], start: int, end: int) -> bytes:
    """The bases of 1-based inclusive `[start, end]` from an indexed plain FASTA."""
    length, offset, line_bases, line_bytes = entry
    start = max(start, 1)
    end = min(end, length)
    if end < start:
        return b""
    first = start - 1
    with fasta.open("rb") as handle:
        handle.seek(offset + first // line_bases * line_bytes + first % line_bases)
        wanted = end - first
        line_end_bytes = line_bytes - line_bases
        # Every line end the read may cross, at the file's own width (CRLF counts two).
        raw = handle.read(wanted + (wanted // line_bases + 2) * line_end_bytes)
    return raw.replace(b"\n", b"").replace(b"\r", b"")[:wanted]


def merge(regions: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Sorted, non-overlapping 1-based inclusive regions."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(regions):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def read_transcripts(cache: Path) -> dict[str, list[dict]]:
    """The cached transcripts per chromosome directory of a JSON cache."""
    transcripts: dict[str, list[dict]] = {}
    for chrom_dir in sorted((cache / "transcripts").iterdir()):
        for shard in sorted(chrom_dir.iterdir()):
            opener = gzip.open if shard.suffix == ".gz" else open
            with opener(shard, "rt", encoding="utf-8") as handle:
                transcripts.setdefault(chrom_dir.name, []).extend(json.load(handle))
    return transcripts


def read_records(vcf: Path) -> dict[str, list[tuple[int, int]]]:
    """Each record's span per chromosome: POS to the farther of POS plus the REF length
    minus one and INFO `END` (a symbolic allele's end)."""
    spans: dict[str, list[tuple[int, int]]] = {}
    for line in vcf.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        chrom, pos, ref = fields[0], int(fields[1]), fields[3]
        info = fields[7] if len(fields) > 7 else ""
        end = pos + len(ref) - 1
        for item in info.split(";"):
            if item.startswith("END="):
                end = max(end, int(item[4:]))
        spans.setdefault(chrom, []).append((pos, end))
    return spans


def read_synonyms(cache: Path) -> dict[str, list[str]]:
    """The names a JSON cache's `chr_synonyms.txt` links to each name, in file order; none
    for a cache without the file."""
    linked: dict[str, list[str]] = {}
    path = cache / "chr_synonyms.txt"
    if path.is_file():
        for line in path.read_text(encoding="utf-8").splitlines():
            names = line.split()
            for synonym in names[1:]:
                linked.setdefault(names[0], []).append(synonym)
                linked.setdefault(synonym, []).append(names[0])
    return linked


def spellings(chrom: str, linked: dict[str, list[str]]) -> list[str]:
    """The names a record's chromosome may stand under, in the order they are tried: as
    written, each name the cache's synonyms file links to it, less a leading `chr` in any
    case, with `chr` in front, the spellings both engines resolve."""
    candidates = [chrom, *linked.get(chrom, [])]
    if chrom[:3].lower() == "chr":
        candidates.append(chrom[3:])
    candidates.append(f"chr{chrom}")
    return candidates


def cache_chromosome(chrom: str, transcripts: dict[str, list[dict]], linked: dict[str, list[str]]) -> str | None:
    """The cache's name for a record's chromosome: the first of its spellings the cache
    holds a transcript under."""
    return next((name for name in spellings(chrom, linked) if name in transcripts), None)


def contig_name(chrom: str, index: dict[str, tuple[int, int, int, int]], label: str) -> str:
    """The FASTA's name for a cache chromosome: as written, or with a `chr` prefix."""
    for name in (chrom, f"chr{chrom}"):
        if name in index:
            return name
    raise fail(f"{label} has no contig named {chrom} or chr{chrom}")


def write_masked_fasta(
    source: Path, out: Path, regions: dict[str, list[tuple[int, int]]], label: str
) -> None:
    """`out` (plain) holds every contig of `regions` at its full length, the bases
    inside the regions copied from `source` and every other base N, plus the
    plain index beside it, the bgzipped copy with its `.fai` and `.gzi`."""
    index = read_fai(Path(f"{source}.fai"))
    fai_lines = []
    offset = 0
    with out.open("wb") as handle:
        for chrom in sorted(regions, key=lambda c: (len(c), c)):
            name = contig_name(chrom, index, label)
            entry = index[name]
            length = entry[0]
            header = f">{chrom}\n".encode()
            handle.write(header)
            offset += len(header)
            fai_lines.append(f"{chrom}\t{length}\t{offset}\t{LINE_BASES}\t{LINE_BASES + 1}")
            sequence = bytearray(b"N" * length)
            for start, end in merge(regions[chrom]):
                start, end = max(start, 1), min(end, length)
                if end >= start:
                    sequence[start - 1 : end] = fetch(source, entry, start, end)
            for i in range(0, length, LINE_BASES):
                line = bytes(sequence[i : i + LINE_BASES]) + b"\n"
                handle.write(line)
                offset += len(line)
    Path(f"{out}.fai").write_text("".join(f"{line}\n" for line in fai_lines), encoding="utf-8")
    compressed = Path(f"{out}.gz")
    run(["bgzip", "-f", "-i", "-I", f"{compressed}.gzi", "-c", str(out)], compressed)
    shutil.copyfile(f"{out}.fai", f"{compressed}.fai")


def run(cmd: list[str], stdout: Path | None = None) -> None:
    try:
        if stdout is None:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        else:
            with stdout.open("wb") as handle:
                subprocess.run(cmd, check=True, stdout=handle)
    except FileNotFoundError as e:
        raise fail(f"{cmd[0]} is not on PATH ({e})") from e
    except subprocess.CalledProcessError as e:
        raise fail(f"{' '.join(cmd)} exited {e.returncode}") from e


def tool_version(tool: str) -> str:
    try:
        out = subprocess.run([tool, "--version"], check=True, capture_output=True, text=True)
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        raise fail(f"{tool} --version failed ({e})") from e
    return (out.stdout or out.stderr).strip().splitlines()[0]


def tools_of(args: argparse.Namespace) -> list[str]:
    """The tools the given files are cut with: `bgzip` for either FASTA, `bigtools` for the bigWig."""
    tools = []
    if args.reference is not None or args.ancestor is not None:
        tools.append("bgzip")
    if args.gerp is not None:
        tools.append("bigtools")
    return tools


def write_gerp(source: Path, out: Path, regions: dict[str, list[tuple[int, int]]], sizes: dict[str, int]) -> None:
    """`out` holds the values of `source` inside `regions` alone, no zoom levels."""
    chroms = sorted(regions, key=lambda c: (len(c), c))
    with tempfile.TemporaryDirectory(prefix="cut_plugin_regions.") as scratch_name:
        scratch = Path(scratch_name)
        bed = scratch / "regions.bed"
        bed.write_text(
            "".join(f"{chrom}\t{max(start, 1) - 1}\t{end}\n" for chrom in chroms for start, end in merge(regions[chrom])),
            encoding="utf-8",
        )
        bedgraph = scratch / "gerp.bedgraph"
        run(["bigtools", "bigwigtobedgraph", str(source), str(bedgraph), str(bed)])
        lines = bedgraph.read_text(encoding="utf-8").splitlines()
        lines.sort(key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
        sorted_bedgraph = scratch / "gerp.sorted.bedgraph"
        sorted_bedgraph.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")
        chrom_sizes = scratch / "chrom.sizes"
        chrom_sizes.write_text("".join(f"{chrom}\t{sizes[chrom]}\n" for chrom in chroms), encoding="utf-8")
        run(["bigtools", "bedgraphtobigwig", "--nzooms", "0", str(sorted_bedgraph), str(chrom_sizes), str(out)])


def bigwig_sizes(bigwig: Path) -> dict[str, int]:
    """The contig lengths a bigWig declares (`bigtools bigwiginfo --chroms`: one
    indented `<name> <id> <length>` line per contig after the header lines)."""
    try:
        out = subprocess.run(["bigtools", "bigwiginfo", "--chroms", str(bigwig)], check=True, capture_output=True, text=True)
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        raise fail(f"bigtools bigwiginfo --chroms {bigwig} failed ({e})") from e
    sizes: dict[str, int] = {}
    for line in out.stdout.splitlines():
        fields = line.split()
        if line.startswith(("\t", " ")) and len(fields) == 3 and fields[2].isdigit():
            sizes[fields[0]] = int(fields[2])
    if not sizes:
        raise fail(f"{bigwig} declares no contigs")
    return sizes


def count_bases(regions: dict[str, list[tuple[int, int]]], sizes: dict[str, int]) -> int:
    """The bases the merged regions cover inside their contigs, flanks past an end left out."""
    bases = 0
    for chrom, items in regions.items():
        for start, end in merge(items):
            start, end = max(start, 1), min(end, sizes[chrom])
            if end >= start:
                bases += end - start + 1
    return bases


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def update_provenance(corpus: Path, files: list[Path], versions: dict[str, str], args: argparse.Namespace) -> None:
    """Rewrites the cut's own keys of `plugin_fixtures` in the corpus's provenance.json and
    keeps every other key there, each another generator's entry for one file, and every
    other key of the record."""
    path = corpus / "provenance.json"
    record = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    earlier = record.get("plugin_fixtures")
    fixtures = earlier if isinstance(earlier, dict) else {}
    fixtures.update(
        {
            "generator": "scripts/fixtures/cut_plugin_regions.py",
            "flank": args.flank,
            "record_flank": args.record_flank,
            "exon_flank": args.exon_flank,
            "gerp_regions": args.gerp_regions,
            "tools": versions,
            "files": {str(f.relative_to(corpus)): sha256(f) for f in sorted(files)},
        }
    )
    record["plugin_fixtures"] = fixtures
    path.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--corpus", required=True, type=Path)
    parser.add_argument("--reference", type=Path, help="the genome FASTA (plain, with its .fai)")
    parser.add_argument("--ancestor", type=Path, help="the ancestral-allele FASTA (plain, with its .fai)")
    parser.add_argument("--gerp", type=Path, help="the GERP bigWig")
    parser.add_argument(
        "--gerp-regions", choices=("exons", "records"), default="exons",
        help="the bigWig regions: each cached exon plus --exon-flank, or each record plus --record-flank",
    )
    parser.add_argument(
        "--gerp-name", default="gerp.bw",
        help="the cut bigWig's file name under plugin_data, the name a plugin derives from its directory",
    )
    parser.add_argument("--flank", type=int, default=5000, help="bases kept around each transcript span")
    parser.add_argument("--record-flank", type=int, default=2000, help="bases kept around each record")
    parser.add_argument("--exon-flank", type=int, default=50, help="bases kept around each exon in the bigWig")
    parser.add_argument("--force", action="store_true", help="overwrite existing fixture files")
    args = parser.parse_args(argv)

    corpus = args.corpus
    cache = corpus / "json_cache"
    vcf = corpus / "variants.vcf"
    if args.reference is None and args.ancestor is None and args.gerp is None:
        raise fail("nothing to cut; give --reference, --ancestor or --gerp")
    inputs = [cache / "transcripts", vcf]
    for fasta in (args.reference, args.ancestor):
        if fasta is not None:
            inputs += [fasta, Path(f"{fasta}.fai")]
    if args.gerp is not None:
        inputs.append(args.gerp)
    for path in inputs:
        if not path.exists():
            raise fail(f"{path} not found")

    plugin_data = corpus / "plugin_data"
    reference_out = corpus / "reference.fa"
    ancestor_out = plugin_data / "human_ancestor.fa"
    if "/" in args.gerp_name or args.gerp_name in ("", ".", ".."):
        raise fail(f"--gerp-name {args.gerp_name!r} is not a file name")
    gerp_out = plugin_data / args.gerp_name
    outputs: list[Path] = []
    if args.reference is not None:
        outputs += [Path(f"{reference_out}.gz"), Path(f"{reference_out}.gz.fai"), Path(f"{reference_out}.gz.gzi"), Path(f"{reference_out}.fai")]
    if args.ancestor is not None:
        outputs += [Path(f"{ancestor_out}.gz"), Path(f"{ancestor_out}.gz.fai"), Path(f"{ancestor_out}.gz.gzi"), Path(f"{ancestor_out}.fai")]
    if args.gerp is not None:
        outputs.append(gerp_out)
    existing = [p for p in outputs if p.exists()]
    if existing and not args.force:
        raise fail(f"{existing[0]} exists; pass --force to overwrite the fixtures")
    if args.ancestor is not None or args.gerp is not None:
        plugin_data.mkdir(exist_ok=True)

    transcripts = read_transcripts(cache)
    records = read_records(vcf)
    # Contig lengths: the reference's index, else the ancestral FASTA's, else the bigWig's own.
    if args.reference is not None or args.ancestor is not None:
        if args.reference is not None:
            fasta, label = args.reference, "the reference"
        else:
            fasta, label = args.ancestor, "the ancestral FASTA"
        fasta_index = read_fai(Path(f"{fasta}.fai"))
        sizes = {chrom: fasta_index[contig_name(chrom, fasta_index, label)][0] for chrom in transcripts}
        lengths = {name: entry[0] for name, entry in fasta_index.items()}
        cut_from = {contig_name(chrom, fasta_index, label): chrom for chrom in transcripts}
    else:
        lengths, label = bigwig_sizes(args.gerp), "the bigWig"
        for chrom in transcripts:
            if chrom not in lengths:
                raise fail(f"the bigWig has no contig named {chrom}")
        sizes = {chrom: lengths[chrom] for chrom in transcripts}
        cut_from = {chrom: chrom for chrom in transcripts}
    span_regions: dict[str, list[tuple[int, int]]] = {}
    exon_regions: dict[str, list[tuple[int, int]]] = {}
    for chrom, items in transcripts.items():
        for t in items:
            span_regions.setdefault(chrom, []).append((int(t["start"]) - args.flank, int(t["end"]) + args.flank))
            for exon in t.get("exons") or []:
                exon_regions.setdefault(chrom, []).append((int(exon["start"]) - args.exon_flank, int(exon["end"]) + args.exon_flank))
    record_regions: dict[str, list[tuple[int, int]]] = {}
    linked = read_synonyms(cache)
    for record_chrom, spans in records.items():
        chrom = cache_chromosome(record_chrom, transcripts, linked)
        if chrom is None:
            # The cache holds no transcript under any spelling of the record's sequence. The
            # record's region is cut on the contig that the first spelling found in the file of
            # the contig lengths names there: a cache chromosome where that contig is the one
            # the chromosome is cut from, else a contig of its own under the file's name.
            contig = next((name for name in spellings(record_chrom, linked) if name in lengths), None)
            if contig is None:
                raise fail(f"variants.vcf names chromosome {record_chrom}, which the cache does not hold and {label} has no contig for")
            chrom = cut_from.get(contig, contig)
            sizes.setdefault(chrom, lengths[contig])
        for start, end in spans:
            region = (start - args.record_flank, end + args.record_flank)
            span_regions.setdefault(chrom, []).append(region)
            record_regions.setdefault(chrom, []).append(region)

    if args.reference is not None:
        write_masked_fasta(args.reference, reference_out, span_regions, "the reference")
        reference_out.unlink()
    if args.ancestor is not None:
        write_masked_fasta(args.ancestor, ancestor_out, record_regions, "the ancestral FASTA")
        ancestor_out.unlink()
    gerp_regions = exon_regions if args.gerp_regions == "exons" else record_regions
    if args.gerp is not None:
        write_gerp(args.gerp, gerp_out, gerp_regions, sizes)

    versions = {tool: tool_version(tool) for tool in tools_of(args)}
    update_provenance(corpus, outputs, versions, args)
    kept = count_bases(span_regions, sizes) if args.reference is not None else 0
    gerp_bases = count_bases(gerp_regions, sizes) if args.gerp is not None else 0
    print(json.dumps({"transcripts": sum(len(v) for v in transcripts.values()), "records": sum(len(v) for v in records.values()),
                      "reference_bases_kept": kept, "gerp_bases": gerp_bases, "files": [str(p.relative_to(corpus)) for p in outputs]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
