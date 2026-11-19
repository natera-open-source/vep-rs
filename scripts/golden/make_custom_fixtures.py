#!/usr/bin/env python3
"""Write a golden corpus's custom annotation sources from its own variants.

A corpus whose flag sets exercise `--custom` needs source files whose records sit on
its variants, named in the manifest as `{corpus}/custom/<file>` so both engines read
one copy. This tool derives them from `variants.vcf` alone, deterministically, so a
regeneration from the same input writes the same records:

    custom/sites.vcf.gz (+ .tbi)    a VCF source: one record per sequence variant locus
                                    in a cycle of classes, plus one inside every
                                    structural variant's span
    custom/regions.bed.gz (+ .tbi)  a BED source: the span of every sites.vcf.gz record
    custom/scores.bw                a bigWig source: the first base of every record,
                                    scored with its AF
    custom/promoters.gff.gz (+ .tbi) a GFF source for gff_type=gencode_promoter: a
                                    promoter_window around every record, a second wider
                                    one around every fifth, and one enhancer line
    custom/clinvar_source.vcf.gz (+ .tbi)  sites.vcf.gz with a `##source=ClinVar` header
                                    line, the header under which Ensembl VEP keeps a
                                    multi-allelic INFO value whole (File/VCF.pm:197-199)
    custom/regions_chr.bed.gz (+ .tbi)  regions.bed.gz with UCSC chromosome names (chr21,
                                    chrM), so a source's names reach the input's through
                                    the synonym and prefix rules
    custom/ucsc/hg38/phyloP100way/hg38.phyloP100way.bw and
    custom/ucsc/hg38/phastCons100way/hg38.phastCons100way.bw
                                    the tree --phyloP 100 and --phastCons 100 read under
                                    --ucsc_data_root custom/ucsc/ --ucsc_assembly hg38: one
                                    base per record under UCSC chromosome names, scored
                                    with the record's ordinal in hundredths minus one
                                    (phyloP, so a score may be negative) and with one
                                    minus its AF (phastCons)

A sequence variant (REF and ALT alleles all bases) takes its class from its ordinal among
the corpus's distinct sequence loci, modulo four, in this order; a later variant at a
locus already written gets no record of its own and shares the first one's hits:

    exact          the variant's REF and first ALT, so a `type=exact` source and a
                   `type=overlap` source both match (ID exact_<n>)
    overlap        the same POS and REF with the first ALT's last base changed, so only a
                   `type=overlap` source matches (ID overlap_<n>)
    miss           no record, so neither matches
    multi_allelic  the first ALT and that changed allele as two ALT alleles with one AF
                   per allele (Number=A), so a `type=exact` source keeps the matching
                   allele's AF and a `type=overlap` source lists both (ID multi_<n>)

A structural variant (a symbolic ALT, with INFO END or SVLEN naming its end) gets one
single-base record at the first position inside its span that no sequence variant
occupies (ID sv_<n>, class sv_span): a `type=overlap` source matches the structural
variant, a `type=exact` source needs equal coordinates and does not. A breakend has no
span and gets none, nor does a record at POS 0, before the first base, which no BED or
bigWig interval can start at. Every record carries INFO AF (Number=A) and CLASS
(Number=1), and FILTER PASS.

`build_golden_corpus.py check --family-columns custom_match=<exact>/<overlap>` reads the
classes a reference realised over these sources; two variants at one locus, or a
sequence variant inside a structural variant's span, share each other's hits, so the
realised class of a row can differ from the class its record was written for.

The promoter windows reach five bases either side of each record and are named PW<n> in
file order, each tied to a fabricated gene id except every third, which carries none, so
the source reports `.`; every fifth record also lies inside a second window twenty bases
either side, so one variant lists two records; the enhancer line over the first record is
of a type the mode does not read, so a run raises its one warning.

bgzip and tabix write the compressed sources and their indexes, `bigtools bedgraphtobigwig`
the bigWig, each found on PATH; the versions they report, the sha256 of every file written
and the record count per class go under `custom_fixtures` in the corpus's provenance.json,
every other key of which is kept. An existing fixture file is overwritten only under
--force.

Usage:
    make_custom_fixtures.py --corpus tests/golden/116/<corpus> [--out <dir>] [--force]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

GENERATOR = "scripts/golden/make_custom_fixtures.py"
CLASS_CYCLE = ("exact", "overlap", "miss", "multi_allelic")
SV_CLASS = "sv_span"
CLASSES = (*CLASS_CYCLE, SV_CLASS)
PHYLOP_BW = "ucsc/hg38/phyloP100way/hg38.phyloP100way.bw"
PHASTCONS_BW = "ucsc/hg38/phastCons100way/hg38.phastCons100way.bw"
FILES = (
    "sites.vcf.gz", "sites.vcf.gz.tbi", "regions.bed.gz", "regions.bed.gz.tbi", "scores.bw",
    "promoters.gff.gz", "promoters.gff.gz.tbi",
    "clinvar_source.vcf.gz", "clinvar_source.vcf.gz.tbi", "regions_chr.bed.gz", "regions_chr.bed.gz.tbi",
    PHYLOP_BW, PHASTCONS_BW,
)
TOOLS = ("bgzip", "tabix", "bigtools")
PROMOTER_FLANK = 5
WIDE_PROMOTER_FLANK = 20
BASES = "ACGT"
SEQUENCE = re.compile(r"^[ACGTNacgtn]+$")
# Primary assembly contig lengths, the chromosome sizes the bigWig's header carries; a
# contig outside the table takes the end of its last interval.
CONTIG_LENGTHS = {
    "GRCh38": {
        "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259,
        "6": 170805979, "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422,
        "11": 135086622, "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189,
        "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
        "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415, "MT": 16569,
    },
    "GRCh37": {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260,
        "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747,
        "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392,
        "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566, "MT": 16569,
    },
}


def fail(message: str) -> SystemExit:
    return SystemExit(f"ERROR: [make_custom_fixtures] {message}")


@dataclass(frozen=True)
class Site:
    """One record of sites.vcf.gz: `pos` is 1-based, `afs` one value per ALT allele."""

    chrom: str
    pos: int
    name: str
    ref: str
    alts: tuple[str, ...]
    afs: tuple[str, ...]
    klass: str

    @property
    def end(self) -> int:
        """The last reference base the record covers, 1-based."""
        return self.pos + len(self.ref) - 1

    @property
    def ordinal(self) -> int:
        """The record's ordinal among the records written, the AF in thousandths."""
        return round(float(self.afs[0]) * 1000)

    def vcf_line(self) -> str:
        info = f"AF={','.join(self.afs)};CLASS={self.klass}"
        return "\t".join([self.chrom, str(self.pos), self.name, self.ref, ",".join(self.alts), ".", "PASS", info])

    def bed_line(self, chrom: str | None = None) -> str:
        return "\t".join([chrom or self.chrom, str(self.pos - 1), str(self.end), self.name, str(self.ordinal), "+"])

    def bedgraph_line(self, chrom: str | None = None, value: str | None = None) -> str:
        return "\t".join([chrom or self.chrom, str(self.pos - 1), str(self.pos), value or self.afs[0]])

    def phylop(self) -> str:
        return f"{self.ordinal / 100 - 1:.2f}"

    def phastcons(self) -> str:
        return f"{1 - float(self.afs[0]):.3f}"


def ucsc_chrom(chrom: str) -> str:
    """A chromosome name as the UCSC files spell it: `chr` prefixed, the mitochondrion `chrM`."""
    if chrom.lower().startswith("chr"):
        return chrom
    return "chrM" if chrom == "MT" else f"chr{chrom}"


def other_allele(alt: str, ref: str) -> str:
    """`alt` with its last base changed to the first base of ACGT that differs from it, and
    from a single-base REF when the allele is one base, so the allele stays a variant."""
    last = alt[-1].upper()
    avoid = {last}
    if len(alt) == 1 and len(ref) == 1:
        avoid.add(ref.upper())
    replacement = next(b for b in BASES if b not in avoid)
    return alt[:-1] + replacement


def af(ordinal: int) -> str:
    """The AF of the `ordinal`th record written: 0.001, 0.002, ...; a multi-allelic record's
    second allele takes the value five hundred records on, so in a file of up to 499 records
    every value is a frequency below one and no two coincide."""
    return f"{0.001 * ordinal:.3f}"


def structural_end(pos: int, ref: str, alt: str, info: str) -> int | None:
    """The last base of a structural variant's span: INFO END, else POS plus |SVLEN|, else the
    end of a REF longer than one base; None for a breakend or an ALT with no span."""
    if "[" in alt or "]" in alt or alt.startswith(".") or alt.endswith("."):
        return None
    fields = dict(part.split("=", 1) for part in info.split(";") if "=" in part)
    if fields.get("END", "").lstrip("-").isdigit():
        return int(fields["END"])
    svlen = fields.get("SVLEN", "").split(",")[0]
    if svlen.lstrip("-").isdigit() and int(svlen) != 0:
        return pos + abs(int(svlen))
    if len(ref) > 1:
        return pos + len(ref) - 1
    return None


def derive_sites(records: list[list[str]]) -> tuple[list[Site], dict[str, int]]:
    """The sites of a corpus's data lines (VCF columns split on tabs) in the order the module
    docstring states, and the count of loci per class, the `miss` loci among them."""
    sites: list[Site] = []
    loci: set[tuple[str, int]] = set()
    counts = {klass: 0 for klass in CLASSES}
    ordinal = 0
    structural: list[tuple[str, int, int]] = []
    for fields in records:
        if len(fields) < 8:
            raise fail(f"a data line has {len(fields)} columns, not the 8 of a VCF record: {fields!r}")
        chrom, pos, ref, alt_column, info = fields[0], int(fields[1]), fields[3], fields[4], fields[7]
        if pos < 1:
            continue
        alts = alt_column.split(",")
        if SEQUENCE.match(ref) and all(SEQUENCE.match(a) for a in alts):
            if (chrom, pos) in loci:
                continue
            loci.add((chrom, pos))
            klass = CLASS_CYCLE[ordinal % len(CLASS_CYCLE)]
            ordinal += 1
            counts[klass] += 1
            if klass == "miss":
                continue
            class_ordinal, site_ordinal = counts[klass], len(sites) + 1
            first = alts[0]
            if klass == "exact":
                sites.append(Site(chrom, pos, f"exact_{class_ordinal}", ref, (first,), (af(site_ordinal),), klass))
            elif klass == "overlap":
                sites.append(Site(chrom, pos, f"overlap_{class_ordinal}", ref, (other_allele(first, ref),), (af(site_ordinal),), klass))
            else:
                second = other_allele(first, ref)
                sites.append(Site(chrom, pos, f"multi_{class_ordinal}", ref, (first, second), (af(site_ordinal), af(500 + site_ordinal)), klass))
        else:
            end = structural_end(pos, ref, alts[0], info)
            if end is not None and end > pos:
                structural.append((chrom, pos, end))
    for chrom, pos, end in structural:
        inside = next((p for p in range(pos + 1, end) if (chrom, p) not in loci), None)
        if inside is None:
            continue
        loci.add((chrom, inside))
        counts[SV_CLASS] += 1
        sites.append(Site(chrom, inside, f"sv_{counts[SV_CLASS]}", "N", ("A",), (af(len(sites) + 1),), SV_CLASS))
    # Chromosome names in byte order, not the corpus's contig order: `bigtools bedgraphtobigwig`
    # refuses a bedGraph whose chromosomes are not in `sort -k1,1` order (so `10` before `2`);
    # tabix takes any chromosome order.
    sites.sort(key=lambda s: (s.chrom, s.pos))
    return sites, counts


def gff_line(chrom: str, feature_type: str, start: int, end: int, attributes: str) -> str:
    return "\t".join([chrom, ".", feature_type, str(start), str(end), ".", "+", ".", attributes])


def promoter_lines(sites: list[Site]) -> list[str]:
    """The records of promoters.gff.gz in the order the module docstring states, sorted by
    chromosome and start for tabix; a window never starts before base 1."""
    def window(site: Site, flank: int, feature_type: str, attributes: str) -> tuple[str, str, int, int, str]:
        return (site.chrom, feature_type, max(1, site.pos - flank), site.end + flank, attributes)

    windows: list[tuple[str, str, int, int, str]] = []
    for ordinal, site in enumerate(sites, start=1):
        gene = f";gene_id=ENSG{90_000_000_000 + ordinal:011d}" if ordinal % 3 else ""
        windows.append(window(site, PROMOTER_FLANK, "promoter_window", f"ID=PW{ordinal}{gene}"))
        if ordinal % 5 == 0:
            windows.append(window(site, WIDE_PROMOTER_FLANK, "promoter_window", f"ID=PW{ordinal}w{gene}"))
    windows.append(window(sites[0], PROMOTER_FLANK, "enhancer", "ID=EH1"))
    windows.sort(key=lambda w: (w[0], w[2], w[3]))
    return [gff_line(*columns) for columns in windows]


def read_records(vcf: Path) -> list[list[str]]:
    lines = vcf.read_text(encoding="utf-8").splitlines()
    return [line.split("\t") for line in lines if line and not line.startswith("#")]


def contig_length(assembly: str, chrom: str, last_end: int) -> int:
    table = CONTIG_LENGTHS.get(assembly, {})
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return max(table.get(name, 0), last_end)


def tool_version(tool: str) -> str:
    """The first line `<tool> --version` prints."""
    if shutil.which(tool) is None:
        raise fail(f"{tool} is not on PATH; bgzip, tabix and bigtools write the sources and their indexes")
    out = subprocess.run([tool, "--version"], capture_output=True, text=True, check=False)
    lines = (out.stdout or out.stderr).strip().splitlines()
    if out.returncode != 0 or not lines:
        raise fail(f"{tool} --version exited {out.returncode} and reported no version")
    return lines[0]


def run(cmd: list[str]) -> None:
    out = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if out.returncode != 0:
        raise fail(f"{' '.join(cmd)} exited {out.returncode}: {out.stderr.strip()}")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fixtures(sites: list[Site], out: Path, assembly: str) -> None:
    """Writes every file of `FILES` under `out` through bgzip, tabix and bigtools, from plain
    files in a scratch directory."""
    out.mkdir(parents=True, exist_ok=True)
    chroms = list(dict.fromkeys(s.chrom for s in sites))
    with tempfile.TemporaryDirectory(dir=out) as scratch_name:
        scratch = Path(scratch_name)
        vcf = scratch / "sites.vcf"
        class_description = "The class of the corpus variant the record was written for: exact, overlap, multi_allelic or sv_span"
        header = [
            "##fileformat=VCFv4.2",
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency, one value per alternate allele">',
            f'##INFO=<ID=CLASS,Number=1,Type=String,Description="{class_description}">',
            *(f"##contig=<ID={c}>" for c in chroms),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ]
        records = [s.vcf_line() for s in sites]
        vcf.write_text("".join(f"{line}\n" for line in [*header, *records]), encoding="utf-8")
        clinvar = scratch / "clinvar_source.vcf"
        clinvar.write_text("".join(f"{line}\n" for line in [header[0], "##source=ClinVar", *header[1:], *records]), encoding="utf-8")
        bed = scratch / "regions.bed"
        bed.write_text("".join(f"{s.bed_line()}\n" for s in sites), encoding="utf-8")
        bed_chr = scratch / "regions_chr.bed"
        bed_chr.write_text("".join(f"{s.bed_line(ucsc_chrom(s.chrom))}\n" for s in sites), encoding="utf-8")
        bedgraph = scratch / "scores.bedgraph"
        bedgraph.write_text("".join(f"{s.bedgraph_line()}\n" for s in sites), encoding="utf-8")
        phylop = scratch / "phylop.bedgraph"
        phylop.write_text("".join(f"{s.bedgraph_line(ucsc_chrom(s.chrom), s.phylop())}\n" for s in sites), encoding="utf-8")
        phastcons = scratch / "phastcons.bedgraph"
        phastcons.write_text("".join(f"{s.bedgraph_line(ucsc_chrom(s.chrom), s.phastcons())}\n" for s in sites), encoding="utf-8")
        gff = scratch / "promoters.gff"
        gff.write_text("".join(f"{line}\n" for line in ["##gff-version 3", *promoter_lines(sites)]), encoding="utf-8")
        sizes = scratch / "chrom.sizes"
        last_end = {c: max(s.pos for s in sites if s.chrom == c) for c in chroms}
        sizes.write_text("".join(f"{c}\t{contig_length(assembly, c, last_end[c])}\n" for c in chroms), encoding="utf-8")
        ucsc_sizes = scratch / "ucsc.chrom.sizes"
        ucsc_sizes.write_text("".join(f"{ucsc_chrom(c)}\t{contig_length(assembly, c, last_end[c])}\n" for c in chroms), encoding="utf-8")
        for plain, preset in ((vcf, "vcf"), (clinvar, "vcf"), (bed, "bed"), (bed_chr, "bed"), (gff, "gff")):
            compressed = scratch / f"{plain.name}.gz"
            with compressed.open("wb") as handle:
                bgzip_run = subprocess.run(["bgzip", "-c", str(plain)], stdout=handle, stderr=subprocess.PIPE, text=True, check=False)
            if bgzip_run.returncode != 0:
                raise fail(f"bgzip -c {plain.name} exited {bgzip_run.returncode}: {bgzip_run.stderr.strip()}")
            run(["tabix", "-p", preset, str(compressed)])
        for graph, chrom_sizes, name in ((bedgraph, sizes, "scores.bw"), (phylop, ucsc_sizes, PHYLOP_BW), (phastcons, ucsc_sizes, PHASTCONS_BW)):
            (scratch / name).parent.mkdir(parents=True, exist_ok=True)
            run(["bigtools", "bedgraphtobigwig", str(graph), str(chrom_sizes), str(scratch / name)])
        for name in FILES:
            (out / name).parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(scratch / name), str(out / name))


def update_provenance(corpus: Path, out: Path, counts: dict[str, int], versions: dict[str, str], assembly: str) -> None:
    """Rewrites the `custom_fixtures` key of the corpus's provenance.json and keeps every other."""
    path = corpus / "provenance.json"
    record = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    if not isinstance(record, dict):
        raise fail(f"{path}: the record is not a JSON object")
    rel = out.relative_to(corpus) if out.is_relative_to(corpus) else out
    record["custom_fixtures"] = {
        "generator": GENERATOR,
        "variants_vcf_sha256": sha256(corpus / "variants.vcf"),
        "assembly": assembly or None,
        "loci": counts,
        "files": {str(rel / name): sha256(out / name) for name in FILES},
        "tools": versions,
    }
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    tmp.replace(path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", required=True, help="the corpus directory (variants.vcf, manifest.json)")
    ap.add_argument("--out", default=None, help="where the sources are written (default <corpus>/custom)")
    ap.add_argument("--force", action="store_true", help="overwrite existing source files")
    args = ap.parse_args(argv)
    corpus = Path(args.corpus)
    vcf = corpus / "variants.vcf"
    if not vcf.is_file():
        raise fail(f"corpus has no variants.vcf: {corpus}")
    out = Path(args.out) if args.out else corpus / "custom"
    existing = [name for name in FILES if (out / name).exists()]
    if existing and not args.force:
        raise fail(f"{out} already holds {', '.join(existing)}; pass --force to overwrite")
    assembly = ""
    manifest = corpus / "manifest.json"
    if manifest.is_file():
        assembly = str(json.loads(manifest.read_text(encoding="utf-8")).get("assembly") or "")
    sites, counts = derive_sites(read_records(vcf))
    if not sites:
        raise fail(f"{vcf} has no sequence variant and no structural variant with a span; nothing to write")
    versions = {tool: tool_version(tool) for tool in TOOLS}
    write_fixtures(sites, out, assembly)
    update_provenance(corpus, out, counts, versions, assembly)
    print(f"wrote {len(FILES)} files under {out}; loci per class: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
