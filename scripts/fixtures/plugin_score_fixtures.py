#!/usr/bin/env python3
"""Write a golden corpus's synthetic plugin data file from its own variants.

The score files of PrimateAI, ClinPred, BayesDel and SpliceVault are not redistributed here. Their
plugins are verified on made-up fixtures instead: a file with the preamble, header and
column layout its shape lists below, bgzip-compressed and tabix-indexed, whose values are
fabricated. This tool writes one such file into `<corpus>/data/` from `variants.vcf`,
deterministically from `--seed`, so a regeneration from the same input and seed writes the same
bytes and the expected outputs committed beside it stay valid. The loftee_gerp shape writes the
conservation table LoFTEE's `master` branch reads, with made-up scores, from the corpus's
`json_cache` instead.

Shapes:

    primateai   eleven `#` preamble lines, then the `#chr ... pos ... alt ... primateDL_score`
                header on line twelve (the line the plugin's example command prefixes with `#`,
                PrimateAI.pm:65-69); one row per alternate base at every single-nucleotide site;
                index `tabix -s 1 -b 2 -e 2`
    clinpred    the `#chr Start Ref Alt ClinPred_Score` header; one row per alternate base
                at every single-nucleotide site; index `tabix -s 1 -b 2 -e 2`
    bayesdel    a `#` header line of five names; one row per alternate base at every
                single-nucleotide site, plus one row per allele that is not a single base
                change (an insertion, a deletion, a longer substitution) in the
                representation the plugin compares (the alleles as Ensembl VEP's parser
                trims them, `-` for the empty allele, at the position the plugin queries);
                index `tabix -s 1 -b 2 -e 2`
    splicevault a twelve-column `#chrom pos ref alt transcript_id ...` header;
                rows at the sites `--transcripts` names, one per
                transcript id listed for the site: at a single-nucleotide site one row per
                alternate base, at an insertion or deletion one row with the alleles as the
                VCF writes them (the anchor base kept);
                the events column carries the `Top<rank>:<type>;<impact>;<percent>;<frame>`
                groups joined by `|`, with a percent above one thousand written with a space
                as its thousands separator; index `tabix -s 1 -b 2 -e 2`
    cadd_snv    the synthetic CADD SNV file: a first line `##CADD <assembly>-v1.7 layout: ...`
                that names the build (`--assembly`) and says the scores are fabricated,
                the `#Chrom Pos Ref Alt RawScore
                PHRED` header, one row per alternate base at every single-nucleotide site,
                a site being a record whose alleles are one base after the common bases at
                both ends are trimmed (the plugin trims both sides before it compares, so a
                record such as `ACG/ATG` is scored by the `C/T` row one base along);
                index `tabix -s 1 -b 2 -e 2`
    cadd_indel  the synthetic CADD indel file: the same first line and header; one row per
                insertion or deletion as the VCF writes it (the anchor base kept),
                beside a decoy row at the same anchor one base
                longer than the allele; a substitution of several bases gets no row;
                index `tabix -s 1 -b 2 -e 2`
    cadd_sv     the synthetic CADD-SV file: a first line `## CADD-SV GRCh38 layout: ...` that says the scores are fabricated, and
                a seven-column header (`#Chrom Start End Type
                Name CADD-SV_PHRED-score CADD-SV_Raw-score`); one row per `<DEL>`, `<INS>`
                or `<DUP>` record with `Start` the VCF position (the fixture's
                zero-based start, the base before the event) and `End` from `INFO/END`
                (`INFO/SVLEN` when there is no `END`; the base after the start for an
                insertion without either), beside two decoys the plugin must not pick:
                the same span under another type, and the same type one base along;
                a record of any other symbolic class gets no row; index
                `tabix -0 -s 1 -b 2 -e 3`
    spliceai_snv
                a SpliceAI score VCF: `##fileformat`, the `##INFO=<ID=SpliceAI,...>`
                line whose description names the tool version (`SpliceAIv1.3`), the eight
                VCF columns; at every single-nucleotide site one row per alternate base and
                per gene symbol `--symbols` lists for the site, its INFO
                `SpliceAI=ALT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL` (delta
                scores with two decimals, delta positions as integers); index `tabix -p vcf`
    spliceai_indel
                the same layout for the alleles that are not a single base change, one row
                per allele and per listed symbol, the alleles as the VCF writes them (the
                anchor base kept); index `tabix -p vcf`
    loftee_gerp the conservation table of LoFTEE's `master` branch, which runs
                `tabix <gerp_file> chr:start-end` and adds up the third tab-separated field
                of every row returned (gerp_dist.pl:92-96): a first line that says the
                scores are fabricated, a `#chrom pos score` header, then one row a base
                (chromosome, position, score) over every exon of the corpus's `json_cache`
                transcripts with `--exon-flank` bases either side;
                index `tabix -s 1 -b 2 -e 2`

A site is every record of `variants.vcf` whose REF and ALT are bases (a symbolic, `*` or
breakend allele writes nothing, and a record carrying a symbolic or breakend allele is a
structural variant that writes nothing, except in the cadd_sv shape, which reads those
records alone), one row per ALT allele of a multi-allelic record.
The alternate bases a site gets are the three bases other than its reference, so the row
the plugin must pick sits beside two it must not. A site named by `--miss chr:pos` (the
position the plugin queries, the VCF position for a single base change) gets no row, so the
corpus carries a variant the file does not score. Scores derive from the seed, the shape
and the row's chromosome, position, reference and alternate allele (and, for splicevault,
its transcript id; for the spliceai shapes, its gene symbol; for cadd_sv, its start, end
and type in place of the position and alleles) alone, so a row's score does not move when
another row is added.

The loftee_gerp shape reads no record: its rows are the bases of the cached transcripts'
exons, the bases the plugin's END_TRUNC walk totals, under the cache's chromosome names. A
base's score, two decimals between -3 and 5, derives from the seed, the chromosome as written
and the position alone, so tables written from one seed over two caches carry the same score
at every base both hold. `--exon-flank N` (default 50) keeps N scored bases either side of
every exon, so a read one base past an exon's end moves a total. A base named by `--miss
chr:pos` gets no row.

The splicevault shape reads `--transcripts <file>`, a two-column tab-separated file of
`chr:pos` (the input's spelling and VCF position) and one transcript id per line: the
plugin annotates a transcript row from the file rows naming its stable id, so the file
names which transcripts of the corpus get a row at each site, and an id no transcript of
the corpus carries is a row the plugin must skip under its transcript match. The spliceai
shapes read `--symbols <file>` of the same form with a gene symbol per line: the plugin
annotates a transcript row from the file row naming its gene symbol, so a symbol no
transcript at the site carries is a row the plugin must skip.

Chromosome names are written as the input spells them (`--chr-style bare`), or with the
`chr` prefix added (`--chr-style chr`); both engines resolve either spelling against the
input's.
`--mitochondrial-name NAME` (`M`, `MT`, `chrM` or `chrMT`) writes the rows of the
mitochondrial contig under NAME, whatever the input spells it, for a table whose files name
that contig otherwise than the input does. Every other contig keeps its spelling, the rows
stay in contig and position order, and a row's score derives from the input's spelling
whichever name is written, so the option moves no score.

`--header-names NAME,NAME,...` (the bayesdel and clinpred shapes, whose header is one line
of column names the plugins read by place) writes the given names in the header line in
place of the layout's, one per column, behind the layout's `#`, for a table whose header
spells or names its columns otherwise than the shape's own header.

`--second-stream-at CHR:POS` writes the file as two block-gzip streams joined, the second
opening with the first row at CHR:POS (the contig and position as the file writes them), so
the first stream's end marker, a block that holds no bytes, stands inside the file before
that row, as in a file made by joining two compressed files. The rows and their scores are
those of the file written as one stream, and tabix indexes the joined file.

The spliceai shapes take four options that write a site's records in another form of the
layout. Each names a site in the input's spelling by its VCF position and is repeatable, a
site stands under one of them alone, and an entry keeps the scores of its own record:

    --join-entries CHR:POS[=N]    each allele's records, one per listed symbol, joined into
                                  records whose `SpliceAI` key holds N entries each (every
                                  entry without N), comma-separated in the records' order,
                                  the form the key's `Number=.` declaration admits for a
                                  site inside several genes
    --no-information CHR:POS      one record per allele whose INFO column is the
                                  missing-value `.`
    --join-alternates CHR:POS     spliceai_snv: one record whose ALT column lists the site's
                                  alternate bases, comma-separated, its key holding every
                                  entry of every alternate
    --any-base-reference CHR:POS  spliceai_snv: the site's records with `N` as their
                                  reference allele

`--repeat-lines CHR:POS=N` (the spliceai shapes; repeatable) writes every record line of the
site N times, the copies next to each other, whichever form its lines take: the plain one,
or the form of one of the four options above, beside which a site may stand under this
option.

`--rescore-lines CHR:POS=N` (the spliceai shapes; repeatable) writes every record line of
the site N times in the same way, the first as its form writes it and each further copy
with the scores of every entry drawn again from the seed, the entry's record and the copy's
place, so the copies of a line name the same alleles and symbols with other scores. A site
stands under this option or under `--repeat-lines`, not both.

`--entry-order CHR:POS=I[:D],I[:D],...` (the spliceai shapes; repeatable) writes one record
per allele of the site, its `SpliceAI` key holding the listed entries in that order. Entry
`I` is the allele's record of place `I`, the records counted from 1 in the file's order (by
symbol), with the scores that record carries; entry `I:D` names the same allele and symbol
with scores drawn again from the seed, the record and the draw `D` (draw 1 is the record's
own), so one symbol can stand in several entries of a record with equal or with other
scores. A site stands under this option or under one of the four above, not both, and may
stand under `--repeat-lines` or `--rescore-lines` beside it.

The file is written under `<corpus>/data/`, or under the directory `--dir` names (a corpus
that passes `--dir_plugins {plugins}` keeps its fixtures under `plugin_data/`, the directory
that placeholder names for vep-rs).

bgzip and tabix write the file and its index, each found on PATH; the versions they report,
the sha256 of both files (and of the transcripts file), the seed, the shape and the row
counts go under `plugin_fixtures.<name>` in the corpus's provenance.json, every other key
of which is kept; a mitochondrial name, header names, the second stream's site, the sites
of the other record forms and the sites whose lines stand more than once, each with its
count, are recorded there when given, and so is the exon flank of a loftee_gerp table. An
existing fixture file is overwritten only under --force.

Usage:
    plugin_score_fixtures.py --corpus tests/golden/116/<corpus> --shape primateai \\
        --name PrimateAI_scores_v0.2_GRCh38_sorted_synthetic.tsv.bgz --seed 38 [--chr-style chr] \\
        [--miss 21:19268555]... [--force]
    plugin_score_fixtures.py --corpus tests/golden/116/<corpus> --shape cadd_snv \\
        --name whole_genome_SNVs_synthetic.tsv.gz --seed 17 --assembly GRCh38 [--miss 21:5030000]...
    plugin_score_fixtures.py --corpus tests/golden/116/<corpus> --shape bayesdel \\
        --name BayesDel_170824_addAF_all_scores_synthetic.txt.gz --seed 5 \\
        [--mitochondrial-name M] [--header-names Chr,Start,Ref,Alt,Score]
    plugin_score_fixtures.py --corpus tests/golden/116/<corpus> --shape cadd_sv \\
        --name 1000G_phase3_SVs_synthetic.tsv.gz --seed 17 --assembly GRCh38
    plugin_score_fixtures.py --corpus tests/golden/116/<corpus> --shape splicevault \\
        --name SpliceVault_data_GRCh38_synthetic.tsv.gz --seed 6 --chr-style chr \\
        --transcripts tests/golden/116/<corpus>/data/splicevault_transcripts.tsv
    plugin_score_fixtures.py --corpus tests/golden/116/<corpus> --shape spliceai_snv \\
        --name spliceai_scores.raw.snv.hg38.synthetic.vcf.gz --seed 9 --dir plugin_data \\
        --symbols tests/golden/116/<corpus>/plugin_data/spliceai_symbols.tsv \\
        [--join-entries 21:17619500[=2]]... [--no-information 21:18279000]... \\
        [--join-alternates 21:21324500]... [--any-base-reference 21:26387000]... \\
        [--entry-order 21:17620200=1,1:2]... \\
        [--repeat-lines 21:17619500=2]... [--rescore-lines 21:17619900=2]... \\
        [--second-stream-at 21:18279000]
    plugin_score_fixtures.py --corpus tests/golden/116/<corpus> --shape loftee_gerp \\
        --name GERP_scores.synthetic.txt.gz --seed 37 --dir plugin_data [--exon-flank 50]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from itertools import groupby
from pathlib import Path

GENERATOR = "scripts/fixtures/plugin_score_fixtures.py"
SPLICEAI_SHAPES = ("spliceai_snv", "spliceai_indel")
LOFTEE_GERP_SHAPE = "loftee_gerp"
SHAPES = ("primateai", "clinpred", "bayesdel", "splicevault", "cadd_snv", "cadd_indel", "cadd_sv", *SPLICEAI_SHAPES, LOFTEE_GERP_SHAPE)
CADD_SHAPES = ("cadd_snv", "cadd_indel", "cadd_sv")
CADD_ASSEMBLIES = ("GRCh37", "GRCh38")
# The shapes whose plugin splits a row into its columns by place and reads no header
# name (BayesDel.pm:151, ClinPred.pm:127).
HEADER_NAME_SHAPES = ("bayesdel", "clinpred")
MITOCHONDRIAL_NAMES = ("M", "MT", "chrM", "chrMT")
TOOLS = ("bgzip", "tabix")
BASES = "ACGT"
# The reference allele `--any-base-reference` writes: the letter for any base.
ANY_BASE = "N"
SITE = re.compile(r"[^:\s]+:[0-9]+")
SEQUENCE = re.compile(r"^[ACGTacgt]+$")
# The alternate alleles that make a record a structural variant (Parser/VCF.pm:241),
# tested on the ALT column as written.
STRUCTURAL = re.compile(r"[<\[\]][^*]+[>\]\[]|^\.\w+|\w+\.$")
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

PRIMATEAI_PREAMBLE = [
    "# Synthetic PrimateAI score file for a vep-rs golden corpus.",
    "# The scores are fabricated from a seed and carry no meaning.",
    "# The layout follows the preparation steps of PrimateAI.pm:57-69:",
    "# eleven lines, comments in this made-up file, the column header",
    "# behind a '#' on line twelve, rows sorted by chromosome and position.",
    "# Columns: chr pos ref alt refAA altAA strand_1pos_0neg",
    "# trinucleotide_context UCSC_gene ExAC_coverage primateDL_score.",
    "# The plugin reads pos, alt and primateDL_score by header name.",
    "# chr may carry the UCSC prefix or the Ensembl name.",
    "# The file is bgzip-compressed and indexed with tabix -s 1 -b 2 -e 2.",
    "# Generated by scripts/fixtures/plugin_score_fixtures.py.",
]
PRIMATEAI_HEADER = (
    "#chr\tpos\tref\talt\trefAA\taltAA\tstrand_1pos_0neg\ttrinucleotide_context"
    "\tUCSC_gene\tExAC_coverage\tprimateDL_score"
)
CLINPRED_HEADER = "#chr\tStart\tRef\tAlt\tClinPred_Score"
BAYESDEL_HEADER = "#Chr\tPos\tRef\tAlt\tBayesDel_addAF"
BAYESDEL_RANGE = (-1.29334, 0.75731)  # the range of the score as the plugin's documentation states it (BayesDel.pm:40)
SPLICEVAULT_HEADER = (
    "#chrom\tpos\tref\talt\ttranscript_id\tsplicevault_site_type\tsplicevault_site_loc"
    "\tspliceai_delta\tsplicevault_out_of_frame_events\tsplicevault_top4_events"
    "\tsplicevault_site_sample_count\tsplicevault_site_max_depth"
)
# The first line of the CADD SNV and indel fixtures names the build each engine looks for
# in the header (CADD.pm:151) and says the file is synthetic. The CADD-SV fixture's header
# holds seven columns.
CADD_MARKER = (
    "##CADD {assembly}-v1.7 layout: a synthetic score file for a vep-rs golden corpus; "
    "the scores are fabricated from a seed and carry no meaning"
)
CADD_HEADER = "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED"
CADD_SV_MARKER = "## CADD-SV GRCh38 layout: a synthetic score file for a vep-rs golden corpus; the scores are fabricated from a seed and carry no meaning"
CADD_SV_HEADER = "#Chrom\tStart\tEnd\tType\tName\tCADD-SV_PHRED-score\tCADD-SV_Raw-score"
CADD_RAW_RANGE = (-3.0, 10.0)
CADD_PHRED_RANGE = (0.0, 50.0)
CADD_SV_RAW_RANGE = (-1.5, 2.5)
CADD_SV_PHRED_RANGE = (0.0, 30.0)
# The symbolic alleles the CADD-SV fixture scores, each under the class word of its rows.
CADD_SV_TYPES = {"<DEL>": "DEL", "<INS>": "INS", "<DUP>": "DUP"}
# How tabix indexes each shape: the start column alone, the zero-based span of the
# CADD-SV file, or the VCF preset of the SpliceAI files, under which a record spans its
# reference allele, so a deletion is found from its deleted bases as well as from its
# anchor base.
INDEX_ARGS = {
    "cadd_sv": ["-0", "-s", "1", "-b", "2", "-e", "3"],
    "spliceai_snv": ["-p", "vcf"],
    "spliceai_indel": ["-p", "vcf"],
}
DEFAULT_INDEX_ARGS = ["-s", "1", "-b", "2", "-e", "2"]
LOFTEE_GERP_MARKER = (
    "## LoFTEE GERP layout: a synthetic conservation table for a vep-rs golden corpus; "
    "the scores are fabricated from a seed and carry no meaning"
)
LOFTEE_GERP_HEADER = "#chrom\tpos\tscore"
# The fabricated score's bounds in hundredths. Their mean is one a base, so walks of a few
# hundred bases total on both sides of the 180 the `master` branch compares with (LoF.pm:268).
LOFTEE_GERP_HUNDREDTHS = (-300, 500)
DEFAULT_EXON_FLANK = 50
SPLICEVAULT_SITE_TYPES = ("Donor_loss", "Acceptor_loss")
SPLICEVAULT_EVENT_TYPES = ("ES", "CD", "CA")
SPLICEVAULT_FRAMES = ("Frameshift", "inFrame")
SPLICEAI_HEADER = [
    "##fileformat=VCFv4.2",
    '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 synthetic fixture, '
    "fabricated values. Format: "
    'ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
]


def fail(message: str) -> SystemExit:
    return SystemExit(f"ERROR: [plugin_score_fixtures] {message}")


@dataclass(frozen=True)
class Allele:
    """One allele of the corpus in the representation the plugins compare: `pos` is the
    position the plugin queries, `ref` and `alt` the alleles as Ensembl VEP's parser trims
    them, `-` for an empty one."""

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1 and self.ref != "-" and self.alt != "-"


def natural_key(chrom: str) -> tuple[int, int | str]:
    """Chromosomes in the order `sort -V` gives them: numbers first, then names."""
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return (0, int(name)) if name.isdigit() else (1, name)


def trimmed(chrom: str, pos: int, ref: str, alts: list[str]) -> list[Allele]:
    """The sequence alleles of one VCF record as Ensembl VEP stores them, each at the
    position the plugins query: the smaller of the record's start and end, which is the
    anchor base for an insertion. A record with an alternate allele of another length than
    the reference loses one shared leading base from every allele and its start moves past
    it (Parser/VCF.pm:295-336; with several alternate alleles only when all of them share
    it, a `*` allele not voting); a record with one alternate allele then loses every
    shared leading and trailing base (Parser.pm:860-881, Sequence.pm:1004-1035). Bases are
    compared as written and upper-cased afterwards, an empty allele is `-`, and alleles of
    one length are never trimmed. Nothing for a structural variant (a symbolic or breakend
    allele), for an allele that is not a sequence, or for one that equals the reference."""
    if STRUCTURAL.search(",".join(alts)) or not SEQUENCE.match(ref):
        return []
    raw_ref, raw_alts = ref, alts
    start, end = pos, pos + len(ref) - 1
    is_indel = any(len(alt) != len(ref) or alt[:1] in ("D", "I") for alt in alts)
    voters = [ref, *(alt for alt in alts if len(alts) == 1 or "*" not in alt)]
    if is_indel and len({allele[:1] for allele in voters}) == 1:
        ref, alts = ref[1:] or "-", [alt[1:] or "-" for alt in alts]
        start += 1
    if len(alts) == 1 and len(raw_ref) != len(raw_alts[0]) and "-" not in raw_ref + raw_alts[0]:
        alt = alts[0]
        while ref and alt and ref[0] == alt[0]:
            ref, alt, start = ref[1:], alt[1:], start + 1
        while ref and alt and ref[-1] == alt[-1]:
            ref, alt, end = ref[:-1], alt[:-1], end - 1
        alts = [alt or "-"]
        ref = ref or "-"
    alleles: list[Allele] = []
    for raw_alt, alt in zip(raw_alts, alts):
        if SEQUENCE.match(raw_alt) and alt.upper() != ref.upper():
            alleles.append(Allele(chrom, min(start, end), ref.upper(), alt.upper()))
    return alleles


def read_alleles(vcf: Path) -> list[Allele]:
    alleles: list[Allele] = []
    for line in vcf.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise fail(f"{vcf}: a record has fewer than five columns: {line!r}")
        chrom, pos, ref, alts = fields[0], int(fields[1]), fields[3], fields[4].split(",")
        alleles.extend(trimmed(chrom, pos, ref, alts))
    return alleles


def spell(chrom: str, style: str) -> str:
    if style == "chr" and not chrom.lower().startswith("chr"):
        return f"chr{chrom}"
    return chrom


def score_rng(seed: int, shape: str, chrom: str, pos: int, ref: str, alt: str) -> random.Random:
    """One generator per row, seeded from the row's identity, so a score depends on
    nothing but the seed and the row."""
    return random.Random(f"{seed}:{shape}:{chrom}:{pos}:{ref}:{alt}")


def primateai_row(rng: random.Random, chrom: str, pos: int, ref: str, alt: str) -> str:
    ref_aa = rng.choice(AMINO_ACIDS)
    alt_aa = rng.choice(AMINO_ACIDS.replace(ref_aa, ""))
    strand = rng.choice("01")
    context = rng.choice(BASES) + ref + rng.choice(BASES)
    gene = f"NM_{rng.randrange(100000, 999999)}.{rng.randrange(1, 5)}"
    coverage = f"{rng.uniform(0, 100):.2f}"
    score = f"{rng.uniform(0, 1):.4f}"
    return "\t".join([chrom, str(pos), ref, alt, ref_aa, alt_aa, strand, context, gene, coverage, score])


def clinpred_row(rng: random.Random, chrom: str, pos: int, ref: str, alt: str) -> str:
    return "\t".join([chrom, str(pos), ref, alt, f"{rng.uniform(0, 1):.6f}"])


def bayesdel_row(rng: random.Random, chrom: str, pos: int, ref: str, alt: str) -> str:
    return "\t".join([chrom, str(pos), ref, alt, f"{rng.uniform(*BAYESDEL_RANGE):.5f}"])


ROW_WRITERS = {"primateai": primateai_row, "clinpred": clinpred_row, "bayesdel": bayesdel_row}
HEADERS = {
    "primateai": [*PRIMATEAI_PREAMBLE, PRIMATEAI_HEADER],
    "clinpred": [CLINPRED_HEADER],
    "bayesdel": [BAYESDEL_HEADER],
    "splicevault": [SPLICEVAULT_HEADER],
    "spliceai_snv": SPLICEAI_HEADER,
    "spliceai_indel": SPLICEAI_HEADER,
    LOFTEE_GERP_SHAPE: [LOFTEE_GERP_MARKER, LOFTEE_GERP_HEADER],
}


def rows_for(alleles: list[Allele], shape: str, seed: int, chr_style: str, misses: set[str]) -> tuple[list[str], dict[str, int]]:
    """The data rows of a shape, sorted as tabix needs them, with the count of rows per
    class: `snv` (three rows per site, one per alternate base), `non_snv` (one row per
    allele that is not a single base change, bayesdel alone) and `miss` (sites left
    without a row)."""
    counts = {"snv": 0, "non_snv": 0, "miss": 0}
    keyed: dict[tuple, str] = {}
    write = ROW_WRITERS[shape]
    snv_sites: set[tuple[str, int, str]] = set()
    missed: set[tuple[str, int]] = set()
    for allele in alleles:
        if f"{allele.chrom}:{allele.pos}" in misses:
            missed.add((allele.chrom, allele.pos))
            continue
        chrom = spell(allele.chrom, chr_style)
        if allele.is_snv:
            snv_sites.add((chrom, allele.pos, allele.ref))
        elif shape == "bayesdel":
            key = (natural_key(chrom), allele.pos, allele.ref, allele.alt)
            if key not in keyed:
                rng = score_rng(seed, shape, chrom, allele.pos, allele.ref, allele.alt)
                keyed[key] = write(rng, chrom, allele.pos, allele.ref, allele.alt)
                counts["non_snv"] += 1
    for chrom, pos, ref in snv_sites:
        for alt in BASES:
            if alt == ref:
                continue
            key = (natural_key(chrom), pos, ref, alt)
            if key not in keyed:
                rng = score_rng(seed, shape, chrom, pos, ref, alt)
                keyed[key] = write(rng, chrom, pos, ref, alt)
                counts["snv"] += 1
    counts["miss"] = len(missed)
    return [keyed[key] for key in sorted(keyed)], counts


@dataclass(frozen=True)
class Record:
    """One sequence allele of a VCF record as the record writes it: `pos` the VCF position,
    `ref` and `alt` the alleles as written, upper-cased, the anchor base of an insertion or
    deletion kept."""

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1


def read_records(vcf: Path) -> list[Record]:
    """The sequence alleles of every record, one per ALT, as the record writes them; a
    structural variant, a `*` allele, a non-sequence allele and an allele equal to the
    reference write nothing, as in `trimmed`."""
    records: list[Record] = []
    for line in vcf.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise fail(f"{vcf}: a record has fewer than five columns: {line!r}")
        chrom, pos, ref, alts = fields[0], int(fields[1]), fields[3], fields[4].split(",")
        if STRUCTURAL.search(",".join(alts)) or not SEQUENCE.match(ref):
            continue
        for alt in alts:
            if SEQUENCE.match(alt) and alt.upper() != ref.upper():
                records.append(Record(chrom, pos, ref.upper(), alt.upper()))
    return records


def read_site_labels(path: Path, label: str) -> dict[str, list[str]]:
    """The labels (`label` names their kind: a transcript id, a gene symbol) listed per
    `chr:pos` site, in file order, each once; a blank line or a `#` line is skipped."""
    sites: dict[str, list[str]] = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 2 or not re.fullmatch(r"[^:\s]+:[0-9]+", fields[0]) or not fields[1].strip():
            raise fail(f"{path}:{number}: a line is chr:pos<TAB>{label}, got {line!r}")
        labels = sites.setdefault(fields[0], [])
        if fields[1] not in labels:
            labels.append(fields[1])
    if not sites:
        raise fail(f"{path}: no site listed")
    return sites


def read_transcripts(path: Path) -> dict[str, list[str]]:
    """The transcript ids listed per `chr:pos` site ([`read_site_labels`])."""
    return read_site_labels(path, "transcript_id")


def read_symbols(path: Path) -> dict[str, list[str]]:
    """The gene symbols listed per `chr:pos` site ([`read_site_labels`])."""
    return read_site_labels(path, "SYMBOL")


def splicevault_percent(rng: random.Random) -> str:
    """A percent of supporting samples in one of the spellings the generator writes: a
    fraction below one with three decimals, an integer up to one hundred, or an integer in
    the thousands with a space as its thousands separator (above 100 where an event is seen
    in more samples than the annotated splicing, SpliceVault.pm:52-54)."""
    shape = rng.randrange(12)
    if shape == 0:
        return f"{rng.uniform(0.001, 0.999):.3f}%"
    if shape == 1:
        return f"{rng.randrange(1000, 99999):,}%".replace(",", " ")
    return f"{rng.randrange(1, 100)}%"


def splicevault_row(rng: random.Random, chrom: str, pos: int, ref: str, alt: str, transcript: str) -> str:
    site_type = rng.choice(SPLICEVAULT_SITE_TYPES)
    site = f"{chrom}:{pos + rng.randrange(-40, 41)}"
    delta = f"{rng.uniform(0.2, 0.99):.2f}"
    events = []
    frameshifts = 0
    for rank in range(1, rng.choice((2, 3, 4, 4, 4)) + 1):
        kind = rng.choice(SPLICEVAULT_EVENT_TYPES)
        if kind == "ES":
            first = rng.randrange(2, 12)
            impact = f"{first}-{first + 1}" if rng.randrange(3) else str(first)
        else:
            impact = f"{rng.choice('+-')}{rng.randrange(1, 500)}"
        frame = rng.choice(SPLICEVAULT_FRAMES)
        frameshifts += frame == "Frameshift"
        events.append(f"Top{rank}:{kind};{impact};{splicevault_percent(rng)};{frame}")
    out_of_frame = f"Frameshift:{frameshifts}/{len(events)}"
    count = str(rng.randrange(100, 200000))
    depth = str(rng.randrange(5, 400))
    return "\t".join([chrom, str(pos), ref, alt, transcript, site_type, site, delta, out_of_frame, "|".join(events), count, depth])


def splicevault_rows(records: list[Record], transcripts: dict[str, list[str]], seed: int, chr_style: str, misses: set[str]) -> tuple[list[str], dict[str, int]]:
    """The rows of the splicevault shape, sorted as tabix needs them (by chromosome,
    position, then reference, alternate allele and transcript id, so a site's rows come in
    one order on every run), with the count of rows per class: `snv` (one per alternate
    base per transcript at a single-nucleotide site), `non_snv` (one per transcript for an
    allele of another length, written as the VCF writes it) and `miss` (sites named by
    `--miss`, left without a row)."""
    counts = {"snv": 0, "non_snv": 0, "miss": 0}
    keyed: dict[tuple, str] = {}
    missed: set[str] = set()
    for record in records:
        site = f"{record.chrom}:{record.pos}"
        if site in misses:
            missed.add(site)
            continue
        chrom = spell(record.chrom, chr_style)
        alts = [alt for alt in BASES if alt != record.ref] if record.is_snv else [record.alt]
        for alt in alts:
            for transcript in transcripts.get(site, []):
                key = (natural_key(chrom), record.pos, record.ref, alt, transcript)
                if key in keyed:
                    continue
                rng = score_rng(seed, "splicevault", chrom, record.pos, record.ref, f"{alt}:{transcript}")
                keyed[key] = splicevault_row(rng, chrom, record.pos, record.ref, alt, transcript)
                counts["snv" if record.is_snv else "non_snv"] += 1
    counts["miss"] = len(missed)
    return [keyed[key] for key in sorted(keyed)], counts


def spliceai_row(rng: random.Random, chrom: str, pos: int, ref: str, alt: str, symbol: str) -> str:
    """One VCF record of a SpliceAI fixture: four delta scores with two decimals and four
    delta positions. Half the fixture's rows keep every score below 0.2, so a corpus run
    under a cutoff shows both verdicts of the flag."""
    ceiling = 0.2 if rng.random() < 0.5 else 1.0
    scores = [f"{rng.uniform(0, ceiling):.2f}" for _ in range(4)]
    positions = [str(rng.randrange(-50, 51)) for _ in range(4)]
    info = "SpliceAI=" + "|".join([alt, symbol, *scores, *positions])
    return "\t".join([chrom, str(pos), ".", ref, alt, ".", ".", info])


@dataclass(frozen=True)
class RecordForms:
    """The sites, in the input's spelling, whose records a spliceai shape writes in another
    form of the layout: `join_entries` (a site and the entries per record, `None` for every
    entry in one), `no_information`, `join_alternates`, `any_base_reference` and
    `entry_order` (a site and the entries of its one record per allele, each the place of a
    symbol among the allele's records and the draw of its scores); and
    `repeat_lines`, the sites whose record lines stand more than once, each with the times
    a line stands, in whichever form the site's lines take, and `rescore_lines`, the sites
    whose further copies of a line carry scores drawn again."""

    join_entries: dict[str, int | None] = field(default_factory=dict)
    no_information: frozenset[str] = frozenset()
    join_alternates: frozenset[str] = frozenset()
    any_base_reference: frozenset[str] = frozenset()
    entry_order: dict[str, tuple[tuple[int, int], ...]] = field(default_factory=dict)
    repeat_lines: dict[str, int] = field(default_factory=dict)
    rescore_lines: dict[str, int] = field(default_factory=dict)

    def of(self, site: str) -> str | None:
        """The form a site's records take, by its option's name."""
        for name in ("join_entries", "no_information", "join_alternates", "any_base_reference", "entry_order"):
            if site in getattr(self, name):
                return name
        return None

    def sites(self) -> dict[str, str]:
        """Every named site with the option that names it, the option of its form where the
        site stands under one beside `--repeat-lines` or `--rescore-lines`."""
        named = {site: f"--{name.replace('_', '-')}" for name in ("join_entries", "no_information", "join_alternates", "any_base_reference", "entry_order") for site in getattr(self, name)}
        return {**{site: "--repeat-lines" for site in self.repeat_lines}, **{site: "--rescore-lines" for site in self.rescore_lines}, **named}


def record_forms_for(
    shape: str,
    join_entries: list[str],
    no_information: list[str],
    join_alternates: list[str],
    any_base_reference: list[str],
    repeat_lines: list[str] | None = None,
    rescore_lines: list[str] | None = None,
    entry_order: list[str] | None = None,
) -> RecordForms | None:
    """The record forms the five options name and the sites `--repeat-lines` names, `None`
    when none is given; refused for a shape that is not a spliceai one, for the two
    single-base forms under the indel shape, for a site that is not `chr:pos`, for an entry
    count below two and for a site under two of the five options. A site may stand under
    `--repeat-lines` or `--rescore-lines` beside one of them; a word of either that is not
    `chr:pos=N`, a count below two, a site named twice by one of the two and a site named
    by both are refused. A word of `--entry-order` is `chr:pos=I[:D],I[:D],...`: two or
    more entries, each the place `I` of a symbol and the draw `D` of its scores (the
    record's own without one), both counted from 1; any other word and a site named twice
    are refused."""
    given = {"--join-entries": join_entries, "--no-information": no_information, "--join-alternates": join_alternates, "--any-base-reference": any_base_reference}
    counted = {"--repeat-lines": repeat_lines or [], "--rescore-lines": rescore_lines or []}
    if not any(given.values()) and not any(counted.values()) and not entry_order:
        return None
    if shape not in SPLICEAI_SHAPES:
        raise fail(f"{', '.join([*given, '--entry-order', *counted])} are given with --shape spliceai_snv or spliceai_indel and with no other shape")
    for option in ("--join-alternates", "--any-base-reference"):
        if given[option] and shape != "spliceai_snv":
            raise fail(f"{option} is given with --shape spliceai_snv: it writes the records of a single-base site")
    joined: dict[str, int | None] = {}
    seen: dict[str, str] = {}
    for option, sites in given.items():
        for named in sites:
            site, _, count = named.partition("=") if option == "--join-entries" else (named, "", "")
            if not SITE.fullmatch(site) or (count and not count.isdigit()):
                raise fail(f"{option} names a site as chr:pos{'[=N]' if option == '--join-entries' else ''}, got {named!r}")
            if count and int(count) < 2:
                raise fail(f"{option} {named}: a record of several entries holds two or more")
            if seen.setdefault(site, option) != option:
                raise fail(f"{site} stands under {seen[site]} and {option}; a site takes one form")
            if option == "--join-entries":
                joined[site] = int(count) if count else None
    ordered: dict[str, tuple[tuple[int, int], ...]] = {}
    for named in entry_order or []:
        site, _, word = named.partition("=")
        slots = [slot.partition(":") for slot in word.split(",")] if word else []
        if not SITE.fullmatch(site) or not slots or any(not place.isdecimal() or (colon and not draw.isdecimal()) for place, colon, draw in slots):
            raise fail(f"--entry-order names a site and its entries as chr:pos=I[:D],I[:D],..., got {named!r}")
        entries = tuple((int(place), int(draw) if colon else 1) for place, colon, draw in slots)
        if len(entries) < 2:
            raise fail(f"--entry-order {named}: a record of several entries holds two or more")
        if any(place < 1 or draw < 1 for place, draw in entries):
            raise fail(f"--entry-order {named}: a place and a draw count from 1")
        if seen.setdefault(site, "--entry-order") != "--entry-order":
            raise fail(f"{site} stands under {seen[site]} and --entry-order; a site takes one form")
        if site in ordered:
            raise fail(f"--entry-order names {site} twice; a site's entries stand in one order")
        ordered[site] = entries
    times: dict[str, dict[str, int]] = {option: {} for option in counted}
    for option, words in counted.items():
        for named in words:
            site, _, count = named.partition("=")
            if not SITE.fullmatch(site) or not count.isdigit():
                raise fail(f"{option} names a site and the times its lines stand as chr:pos=N, got {named!r}")
            if int(count) < 2:
                raise fail(f"{option} {named}: a line that stands more than once stands twice or more")
            if site in times[option]:
                raise fail(f"{option} names {site} twice; a site's lines stand one number of times")
            times[option][site] = int(count)
    for site in times["--rescore-lines"]:
        if site in times["--repeat-lines"]:
            raise fail(f"{site} stands under --repeat-lines and --rescore-lines; the copies of a site's lines are written one way")
    return RecordForms(joined, frozenset(no_information), frozenset(join_alternates), frozenset(any_base_reference), ordered, times["--repeat-lines"], times["--rescore-lines"])


def in_record_forms(
    ordered: list[tuple[str, list[str]]],
    forms: RecordForms,
    rescored: Callable[[list[str], int], str] | None = None,
    drawn: Callable[[list[str], int], str] | None = None,
) -> tuple[list[str], dict[str, int]]:
    """The rows of `ordered` (each a site in the input's spelling and a record's eight
    cells, in file order) with the sites `forms` names written in their form, and the
    records written per form. An entry is a record's value behind `SpliceAI=`, so its
    allele, symbol and scores are those of the record it comes from. Every line of a site
    under `repeat_lines` then stands that many times, the copies next to each other, and a
    form counts the lines the file holds in it, `repeat_lines` the lines of its sites. A
    site under `rescore_lines` holds its lines as many times in the same places, each copy
    behind the first being `rescored(cells, place)`, the line of those cells with its
    entries' scores drawn again for that place (two for the second copy). A site under `entry_order` holds one record per
    allele whose entries are the ones its word lists, in that order: entry `I:D` is the
    value of the allele's record of place `I`, the records standing in the order of
    `ordered`, as that record writes it for draw 1 and else `drawn(cells, D)`, the entry of
    those cells with its scores drawn again for that draw, so one symbol can stand in
    several entries of a record with equal or with other scores. A named site
    without a record, a site joined into records of several entries whose alleles
    have one record each, and an entry whose place is past the allele's records are
    refused."""
    rows: list[str] = []
    written = {name: 0 for name in ("join_entries", "no_information", "join_alternates", "any_base_reference", "entry_order", "repeat_lines", "rescore_lines") if getattr(forms, name)}
    with_rows = {site for site, _ in ordered}
    for site, option in forms.sites().items():
        if site not in with_rows:
            raise fail(f"{option} names {site}, a site with no record in this file")

    def entry(cells: list[str]) -> str:
        return cells[7].removeprefix("SpliceAI=")

    for site, of_site in groupby(ordered, key=lambda item: item[0]):
        records = [cells for _, cells in of_site]
        form = forms.of(site)
        first, in_form = len(rows), written.get(form, 0)
        if form is None:
            rows += ["\t".join(cells) for cells in records]
        elif form == "any_base_reference":
            rows += ["\t".join([*cells[:3], ANY_BASE, *cells[4:]]) for cells in records]
            written[form] += len(records)
        elif form == "join_alternates":
            alts = list(dict.fromkeys(cells[4] for cells in records))
            rows.append("\t".join([*records[0][:4], ",".join(alts), *records[0][5:7], "SpliceAI=" + ",".join(entry(cells) for cells in records)]))
            written[form] += 1
        else:
            for _, of_allele in groupby(records, key=lambda cells: (cells[3], cells[4])):
                same_allele = list(of_allele)
                if form == "no_information":
                    rows.append("\t".join([*same_allele[0][:7], "."]))
                    written[form] += 1
                    continue
                if form == "entry_order":
                    if drawn is None:
                        raise fail(f"--entry-order names {site}; an entry of a further draw takes its scores from the shape's seed")
                    listed = []
                    for place, draw in forms.entry_order[site]:
                        if place > len(same_allele):
                            raise fail(f"--entry-order names entry {place} of {site}, whose alleles have {len(same_allele)} record{'' if len(same_allele) == 1 else 's'} each")
                        cells = same_allele[place - 1]
                        listed.append(entry(cells) if draw == 1 else drawn(cells, draw))
                    rows.append("\t".join([*same_allele[0][:7], "SpliceAI=" + ",".join(listed)]))
                    written[form] += 1
                    continue
                if len(same_allele) < 2:
                    raise fail(f"--join-entries names {site}, whose alleles have one record each; a record of several entries needs two symbols listed for the site")
                size = forms.join_entries[site] or len(same_allele)
                for start in range(0, len(same_allele), size):
                    group = same_allele[start : start + size]
                    rows.append("\t".join([*group[0][:7], "SpliceAI=" + ",".join(entry(cells) for cells in group)]))
                    written[form] += 1
        times = forms.repeat_lines.get(site)
        if times:
            rows[first:] = [line for line in rows[first:] for _ in range(times)]
            written["repeat_lines"] += len(rows) - first
        copies = forms.rescore_lines.get(site)
        if copies:
            if rescored is None:
                raise fail(f"--rescore-lines names {site}; a further copy of a line takes its scores from the shape's seed")
            rows[first:] = [line if place == 1 else rescored(line.split("\t"), place) for line in rows[first:] for place in range(1, copies + 1)]
            written["rescore_lines"] += len(rows) - first
        if (times or copies) and form is not None:
            written[form] = in_form + (written[form] - in_form) * (times or copies)
    return rows, written


def spliceai_rows(records: list[Record], symbols: dict[str, list[str]], shape: str, seed: int, chr_style: str, misses: set[str], forms: RecordForms | None = None) -> tuple[list[str], dict[str, int]]:
    """The rows of a spliceai shape, sorted as tabix needs them (by chromosome, position,
    then reference, alternate allele and symbol), with the count of rows per class: `snv`
    (the `spliceai_snv` shape: one per alternate base per symbol at a single-nucleotide
    site), `non_snv` (the `spliceai_indel` shape: one per symbol for an allele of another
    length, written as the VCF writes it) and `miss` (sites named by `--miss`, left without
    a row). A record of the other shape's class writes nothing.
    Under `forms` the named sites' records are written in their
    form ([`in_record_forms`]): the class counts the rows then written, and each form given
    counts its records under its own name. A further copy of a line under `rescore_lines`
    takes each entry from a row seeded as the entry's own record is, with the copy's place
    behind the symbol, and an entry of a further draw under `entry_order` from a row seeded
    with the word `entry` and the draw behind the symbol."""
    counts = {"snv": 0, "non_snv": 0, "miss": 0}
    keyed: dict[tuple, str] = {}
    site_of: dict[tuple, str] = {}
    missed: set[str] = set()
    wants_snv = shape == "spliceai_snv"
    for record in records:
        if record.is_snv != wants_snv:
            continue
        site = f"{record.chrom}:{record.pos}"
        if site in misses:
            missed.add(site)
            continue
        chrom = spell(record.chrom, chr_style)
        alts = [alt for alt in BASES if alt != record.ref] if record.is_snv else [record.alt]
        for alt in alts:
            for symbol in symbols.get(site, []):
                key = (natural_key(chrom), record.pos, record.ref, alt, symbol)
                if key in keyed:
                    continue
                rng = score_rng(seed, shape, chrom, record.pos, record.ref, f"{alt}:{symbol}")
                keyed[key] = spliceai_row(rng, chrom, record.pos, record.ref, alt, symbol)
                site_of[key] = site
                counts["snv" if record.is_snv else "non_snv"] += 1
    counts["miss"] = len(missed)
    if forms is None:
        return [keyed[key] for key in sorted(keyed)], counts

    def rescored(cells: list[str], place: int) -> str:
        if not cells[7].startswith("SpliceAI="):
            return "\t".join(cells)
        chrom, pos, ref = cells[0], int(cells[1]), cells[3]
        entries = []
        for entry in cells[7].removeprefix("SpliceAI=").split(","):
            alt, symbol = entry.split("|")[:2]
            rng = score_rng(seed, shape, chrom, pos, ref, f"{alt}:{symbol}:{place}")
            entries.append(spliceai_row(rng, chrom, pos, ref, alt, symbol).split("\t")[7].removeprefix("SpliceAI="))
        return "\t".join([*cells[:7], "SpliceAI=" + ",".join(entries)])

    def drawn(cells: list[str], draw: int) -> str:
        chrom, pos, ref = cells[0], int(cells[1]), cells[3]
        alt, symbol = cells[7].removeprefix("SpliceAI=").split("|")[:2]
        rng = score_rng(seed, shape, chrom, pos, ref, f"{alt}:{symbol}:entry{draw}")
        return spliceai_row(rng, chrom, pos, ref, alt, symbol).split("\t")[7].removeprefix("SpliceAI=")

    rows, written = in_record_forms([(site_of[key], keyed[key].split("\t")) for key in sorted(keyed)], forms, rescored, drawn)
    counts["snv" if wants_snv else "non_snv"] = len(rows)
    return rows, {**counts, **written}


def cadd_row(rng: random.Random, chrom: str, pos: int, ref: str, alt: str) -> str:
    raw = f"{rng.uniform(*CADD_RAW_RANGE):.6f}"
    phred = f"{rng.uniform(*CADD_PHRED_RANGE):.3f}"
    return "\t".join([chrom, str(pos), ref, alt, raw, phred])


def trimmed_both_ends(ref: str, alt: str) -> tuple[str, str, int]:
    """The alleles without the bases they share at the start and then the end, with the
    offset of the first base kept, as `get_matched_variant_alleles` trims both sides of a
    match (ensembl-variation Utils/Sequence.pm:1004-1035); `-` for an allele trimmed away."""
    offset = 0
    while ref and alt and ref[0] == alt[0]:
        ref, alt, offset = ref[1:], alt[1:], offset + 1
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return ref or "-", alt or "-", offset


def cadd_snv_rows(alleles: list[Allele], seed: int, chr_style: str, misses: set[str]) -> tuple[list[str], dict[str, int]]:
    """The rows of the cadd_snv shape: three per single-nucleotide site, the site of an
    allele being its one changed base after the common bases at both ends are trimmed
    (`-` alleles never trim, so an insertion or deletion is no site); a `--miss` site
    (the allele's position before that trim) gets none."""
    counts = {"snv": 0, "non_snv": 0, "miss": 0}
    keyed: dict[tuple, str] = {}
    sites: set[tuple[str, int, str]] = set()
    missed: set[tuple[str, int]] = set()
    for allele in alleles:
        if f"{allele.chrom}:{allele.pos}" in misses:
            missed.add((allele.chrom, allele.pos))
            continue
        ref, alt, offset = trimmed_both_ends(allele.ref, allele.alt)
        if len(ref) == 1 and len(alt) == 1 and ref != "-" and alt != "-":
            sites.add((spell(allele.chrom, chr_style), allele.pos + offset, ref))
    for chrom, pos, ref in sites:
        for alt in BASES:
            if alt == ref:
                continue
            key = (natural_key(chrom), pos, ref, alt)
            if key not in keyed:
                keyed[key] = cadd_row(score_rng(seed, "cadd_snv", chrom, pos, ref, alt), chrom, pos, ref, alt)
                counts["snv"] += 1
    counts["miss"] = len(missed)
    return [keyed[key] for key in sorted(keyed)], counts


def cadd_indel_rows(records: list[Record], seed: int, chr_style: str, misses: set[str]) -> tuple[list[str], dict[str, int]]:
    """The rows of the cadd_indel shape: one per insertion or deletion as the VCF writes it,
    and a decoy at the same anchor whose allele is one base longer (the inserted bases
    with their last base repeated; the deleted bases with one more); a substitution of
    any length gets no row; a `--miss` site (the VCF position) gets none."""
    counts = {"snv": 0, "non_snv": 0, "miss": 0}
    keyed: dict[tuple, str] = {}
    missed: set[str] = set()
    for record in records:
        if len(record.ref) == len(record.alt):
            continue
        site = f"{record.chrom}:{record.pos}"
        if site in misses:
            missed.add(site)
            continue
        chrom = spell(record.chrom, chr_style)
        if len(record.alt) > len(record.ref):
            decoy = (record.ref, record.alt + record.alt[-1])
        else:
            decoy = (record.ref + record.ref[-1], record.alt)
        for ref, alt in ((record.ref, record.alt), decoy):
            key = (natural_key(chrom), record.pos, ref, alt)
            if key not in keyed:
                keyed[key] = cadd_row(score_rng(seed, "cadd_indel", chrom, record.pos, ref, alt), chrom, record.pos, ref, alt)
                counts["non_snv"] += 1
    counts["miss"] = len(missed)
    return [keyed[key] for key in sorted(keyed)], counts


@dataclass(frozen=True)
class Structural:
    """One structural record the CADD-SV fixture scores: `start` the VCF position, the
    fixture's zero-based start, `end` the last base of the event, `kind` the class
    word of its rows."""

    chrom: str
    start: int
    end: int
    kind: str


def read_structural(vcf: Path) -> list[Structural]:
    """The `<DEL>`, `<INS>` and `<DUP>` records of the VCF, one per record, with the end
    `INFO/END` names, else the start plus `INFO/SVLEN`'s size, else, for an insertion, the
    base after the start; a deletion with neither writes nothing, as Ensembl VEP skips it
    (Parser/VCF.pm:555-557), and so does a duplication with neither, which Ensembl VEP keeps
    as a one-base event at the base after its position (BaseVCF4.pm:310-311). A record of
    any other allele writes nothing."""
    records: list[Structural] = []
    for line in vcf.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise fail(f"{vcf}: a record has fewer than five columns: {line!r}")
        chrom, pos, alt = fields[0], int(fields[1]), fields[4]
        kind = CADD_SV_TYPES.get(alt)
        if kind is None:
            continue
        info_items = fields[7].split(";") if len(fields) > 7 else []
        info = dict(item.split("=", 1) for item in info_items if "=" in item)
        if "END" in info:
            end = int(info["END"])
        elif "SVLEN" in info:
            end = pos + abs(int(info["SVLEN"].split(",")[0]))
        elif kind == "INS":
            end = pos + 1
        else:
            continue
        records.append(Structural(chrom, pos, end, kind))
    return records


def cadd_sv_row(rng: random.Random, chrom: str, start: int, end: int, kind: str) -> str:
    name = f"fixture_{kind}_{chrom}:{start}-{end}"
    phred = f"{rng.uniform(*CADD_SV_PHRED_RANGE):.3f}"
    raw = f"{rng.uniform(*CADD_SV_RAW_RANGE):.15f}"
    return "\t".join([chrom, str(start), str(end), kind, name, phred, raw])


def cadd_sv_rows(records: list[Structural], seed: int, chr_style: str, misses: set[str]) -> tuple[list[str], dict[str, int]]:
    """The rows of the cadd_sv shape, sorted by chromosome, start, end and type: each
    record's own row and two decoys, the same span under another type (a deletion's
    duplication and the reverse, an insertion's deletion) and the same type one base along;
    a `--miss` site (the VCF position) gets none. Rows are counted as `sv` (the records'
    own) and `decoy`."""
    counts = {"sv": 0, "decoy": 0, "miss": 0}
    keyed: dict[tuple, str] = {}
    missed: set[str] = set()
    other = {"DEL": "DUP", "DUP": "DEL", "INS": "DEL"}
    for record in records:
        site = f"{record.chrom}:{record.start}"
        if site in misses:
            missed.add(site)
            continue
        chrom = spell(record.chrom, chr_style)
        rows = [
            (record.start, record.end, record.kind, True),
            (record.start, record.end, other[record.kind], False),
            (record.start + 1, record.end, record.kind, False),
        ]
        for start, end, kind, own in rows:
            key = (natural_key(chrom), start, end, kind)
            if key not in keyed:
                rng = score_rng(seed, "cadd_sv", chrom, start, f"{end}", kind)
                keyed[key] = cadd_sv_row(rng, chrom, start, end, kind)
                counts["sv" if own else "decoy"] += 1
    counts["miss"] = len(missed)
    return [keyed[key] for key in sorted(keyed)], counts


def read_exon_regions(cache: Path, flank: int) -> dict[str, list[tuple[int, int]]]:
    """Per chromosome directory of a JSON cache, the exons of every cached transcript with
    `flank` bases either side, merged: sorted 1-based inclusive spans that neither overlap
    nor touch, none starting below 1. A chromosome without an exon is left out."""
    regions: dict[str, list[tuple[int, int]]] = {}
    for chrom_dir in sorted(path for path in (cache / "transcripts").iterdir() if path.is_dir()):
        spans: list[tuple[int, int]] = []
        for shard in sorted(chrom_dir.iterdir()):
            opener = gzip.open if shard.suffix == ".gz" else open
            with opener(shard, "rt", encoding="utf-8") as handle:
                transcripts = json.load(handle)
            for transcript in transcripts:
                for exon in (transcript or {}).get("exons") or []:
                    spans.append((max(int(exon["start"]) - flank, 1), int(exon["end"]) + flank))
        merged: list[tuple[int, int]] = []
        for start, end in sorted(spans):
            if merged and start <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        if merged:
            regions[chrom_dir.name] = merged
    return regions


def loftee_gerp_score(seed: int, chrom: str, pos: int) -> str:
    """The fabricated score of one base, two decimals within `LOFTEE_GERP_HUNDREDTHS`: a
    function of the seed, the chromosome as the file writes it and the position."""
    hundredths = score_rng(seed, LOFTEE_GERP_SHAPE, chrom, pos, "", "").randint(*LOFTEE_GERP_HUNDREDTHS)
    return f"{hundredths / 100:.2f}"


def loftee_gerp_rows(regions: dict[str, list[tuple[int, int]]], seed: int, chr_style: str, misses: set[str]) -> tuple[list[str], dict[str, int]]:
    """The rows of the loftee_gerp shape, one per base of `regions`, each chromosome's rows
    together and in position order as tabix needs them; a base `--miss` names in the cache's
    spelling gets none. Counted as `bases` (the rows), `regions` and `miss`."""
    counts = {"bases": 0, "regions": 0, "miss": 0}
    rows: list[str] = []
    for chrom in sorted(regions, key=natural_key):
        name = spell(chrom, chr_style)
        for start, end in regions[chrom]:
            counts["regions"] += 1
            for pos in range(start, end + 1):
                if f"{chrom}:{pos}" in misses:
                    counts["miss"] += 1
                    continue
                rows.append(f"{name}\t{pos}\t{loftee_gerp_score(seed, name, pos)}")
    counts["bases"] = len(rows)
    return rows, counts


def with_mitochondrial_name(rows: list[str], name: str) -> list[str]:
    """`rows` with the mitochondrial contig, in any of its spellings, written as `name`:
    the first column alone changes, so a row keeps the score its own spelling seeds, and
    the rows are put back in contig and position order, each contig's rows together, as
    tabix needs them. One position the input names under two of the spellings is refused,
    since the file holds the contig once."""
    spellings = {known.upper() for known in MITOCHONDRIAL_NAMES}
    seen: dict[int, str] = {}
    keyed: list[tuple[tuple, str]] = []
    for row in rows:
        chrom, pos, rest = row.split("\t", 2)
        if chrom.upper() in spellings:
            if seen.setdefault(int(pos), chrom) != chrom:
                raise fail(f"--mitochondrial-name {name}: position {pos} of the mitochondrial contig has rows as {seen[int(pos)]} and as {chrom}; the file holds the contig once")
            chrom = name
        keyed.append(((natural_key(chrom), int(pos)), f"{chrom}\t{pos}\t{rest}"))
    return [row for _, row in sorted(keyed, key=lambda item: item[0])]


def header_names_for(shape: str, given: str) -> list[str]:
    """The names `--header-names` gives a shape's header line, one per column of the
    layout; refused for a shape whose header the plugin reads by name or that has several
    header lines, for another count than the layout's, and for a name that is empty,
    holds white space or opens with the `#` the layout writes itself."""
    if shape not in HEADER_NAME_SHAPES:
        raise fail(f"--header-names is given with --shape {' or '.join(HEADER_NAME_SHAPES)} and with no other shape")
    names = given.split(",")
    layout = HEADERS[shape][0].split("\t")
    if len(names) != len(layout):
        raise fail(f"--header-names: {len(names)} given, the {shape} header has {len(layout)} columns")
    for name in names:
        if not name or name.startswith("#") or re.search(r"\s", name):
            raise fail(f"--header-names gives a column the name {name!r}; a name is not empty, holds no white space and does not open with #")
    return names


def headers_for(shape: str, assembly: str | None, names: list[str] | None = None) -> list[str]:
    """The header lines of a shape: the CADD shapes name the build in their first line,
    and `names` stand in the one header line in place of the layout's, behind its `#`."""
    if shape in ("cadd_snv", "cadd_indel"):
        return [CADD_MARKER.format(assembly=assembly), CADD_HEADER]
    if shape == "cadd_sv":
        return [CADD_SV_MARKER, CADD_SV_HEADER]
    if names:
        prefix = "#" if HEADERS[shape][0].startswith("#") else ""
        return [prefix + "\t".join(names)]
    return HEADERS[shape]


def tool_version(tool: str) -> str:
    """The first line `<tool> --version` prints."""
    if shutil.which(tool) is None:
        raise fail(f"{tool} is not on PATH; bgzip and tabix write the file and its index")
    out = subprocess.run([tool, "--version"], capture_output=True, text=True, check=False)
    lines = (out.stdout or out.stderr).strip().splitlines()
    if out.returncode != 0 or not lines:
        raise fail(f"{tool} --version exited {out.returncode} and reported no version")
    return lines[0]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def second_stream_line(lines: list[str], site: str) -> int:
    """The index in `lines` of the first row at `site`, the contig and position as the file
    writes them, where `--second-stream-at` opens the second stream; refused for a site
    that is not `chr:pos` and for one no row stands at."""
    if not SITE.fullmatch(site):
        raise fail(f"--second-stream-at names a site as chr:pos, got {site!r}")
    chrom, pos = site.split(":")
    for index, line in enumerate(lines):
        if not line.startswith("#") and line.split("\t")[:2] == [chrom, pos]:
            return index
    raise fail(f"--second-stream-at names {site}, where the file has no row")


def write_fixture(lines: list[str], out: Path, index_args: list[str] | None = None, second_stream: int | None = None) -> None:
    """Writes `lines` through bgzip to `out` and tabix to `out.tbi` (`-s 1 -b 2 -e 2` unless
    `index_args` names the columns), from plain files in a scratch directory beside the
    target. Under `second_stream` the lines from that index on are a second bgzip stream
    written behind the first, whose end marker then stands inside the file."""
    out.parent.mkdir(parents=True, exist_ok=True)
    streams = [lines] if second_stream is None else [lines[:second_stream], lines[second_stream:]]
    with tempfile.TemporaryDirectory(dir=out.parent) as scratch_name:
        scratch = Path(scratch_name)
        compressed = scratch / out.name
        with compressed.open("wb") as handle:
            for number, stream in enumerate(streams):
                plain = scratch / f"fixture.{number}.tsv"
                plain.write_text("".join(f"{line}\n" for line in stream), encoding="utf-8")
                bgzip_run = subprocess.run(["bgzip", "-c", str(plain)], stdout=handle, stderr=subprocess.PIPE, text=True, check=False)
                if bgzip_run.returncode != 0:
                    raise fail(f"bgzip -c exited {bgzip_run.returncode}: {bgzip_run.stderr.strip()}")
        tabix_run = subprocess.run(["tabix", "-f", *(index_args or DEFAULT_INDEX_ARGS), str(compressed)], capture_output=True, text=True, check=False)
        if tabix_run.returncode != 0:
            raise fail(f"tabix exited {tabix_run.returncode}: {tabix_run.stderr.strip()}")
        shutil.move(str(compressed), str(out))
        shutil.move(f"{compressed}.tbi", f"{out}.tbi")


def update_provenance(corpus: Path, out: Path, args: argparse.Namespace, counts: dict[str, int], versions: dict[str, str]) -> None:
    """Rewrites `plugin_fixtures.<file name>` in the corpus's provenance.json and keeps every
    other key, other fixtures of the corpus included. A splicevault entry also names its
    transcripts file, a spliceai entry its symbols file (its resolved path, relative to the
    corpus when it lies inside it), each with its sha256."""
    path = corpus / "provenance.json"
    record = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    if not isinstance(record, dict):
        raise fail(f"{path}: the record is not a JSON object")
    fixtures = record.setdefault("plugin_fixtures", {})
    if not isinstance(fixtures, dict):
        raise fail(f"{path}: plugin_fixtures is not a JSON object")
    rel = out.relative_to(corpus)
    entry = {
        "generator": GENERATOR,
        "shape": args.shape,
        "seed": args.seed,
        "chr_style": args.chr_style,
        "misses": sorted(args.miss or []),
        "variants_vcf_sha256": sha256(corpus / "variants.vcf"),
        "rows": counts,
        "files": {str(rel): sha256(out), f"{rel}.tbi": sha256(Path(f"{out}.tbi"))},
        "tools": versions,
    }
    if args.assembly:
        entry["assembly"] = args.assembly
    if args.exon_flank is not None:
        entry["exon_flank"] = args.exon_flank
    if args.mitochondrial_name:
        entry["mitochondrial_name"] = args.mitochondrial_name
    if args.header_names:
        entry["header_names"] = args.header_names.split(",")
    if args.second_stream_at:
        entry["second_stream_at"] = args.second_stream_at
    forms = {name: sorted(given) for name in ("join_entries", "no_information", "join_alternates", "any_base_reference", "entry_order", "repeat_lines", "rescore_lines") if (given := getattr(args, name))}
    if forms:
        entry["record_forms"] = forms
    for key, given in (("transcripts", args.transcripts), ("symbols", args.symbols)):
        if given:
            listed = Path(given)
            shown = listed.resolve()
            if shown.is_relative_to(corpus.resolve()):
                shown = shown.relative_to(corpus.resolve())
            entry[key] = {str(shown): sha256(listed)}
    fixtures[out.name] = entry
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    tmp.replace(path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", required=True, help="the corpus directory holding variants.vcf")
    ap.add_argument("--shape", required=True, choices=SHAPES, help="the file layout written")
    ap.add_argument("--name", required=True, help="the file name under <corpus>/<dir>/ (its .tbi sits beside it)")
    ap.add_argument("--dir", "--data-dir", dest="dir", default="data", help="the directory inside the corpus the file is written to (default: data)")
    ap.add_argument("--seed", required=True, type=int, help="the seed every score derives from")
    ap.add_argument("--chr-style", default="bare", choices=("bare", "chr"), help="chromosome spelling: the input's, or with the chr prefix")
    ap.add_argument("--mitochondrial-name", choices=MITOCHONDRIAL_NAMES, help="the name the file gives the mitochondrial contig, whatever the input spells it")
    ap.add_argument("--header-names", metavar="NAME,NAME,...", help="bayesdel, clinpred: the column names of the header line, in place of the layout's")
    ap.add_argument("--miss", action="append", metavar="CHR:POS", help="a site (the input's spelling) left without a row (repeatable)")
    ap.add_argument("--transcripts", metavar="FILE", help="splicevault: the chr:pos<TAB>transcript_id lines naming the rows written")
    ap.add_argument("--symbols", metavar="FILE", help="spliceai_snv, spliceai_indel: the chr:pos<TAB>SYMBOL lines naming the rows written")
    ap.add_argument("--join-entries", action="append", metavar="CHR:POS[=N]", help="spliceai shapes: the site's records joined into records of N entries (every entry without N) (repeatable)")
    ap.add_argument("--no-information", action="append", metavar="CHR:POS", help="spliceai shapes: the site's records written with the missing-value INFO (repeatable)")
    ap.add_argument("--join-alternates", action="append", metavar="CHR:POS", help="spliceai_snv: the site's records joined into one whose ALT lists the alternate bases (repeatable)")
    ap.add_argument("--any-base-reference", action="append", metavar="CHR:POS", help="spliceai_snv: the site's records written with N as their reference allele (repeatable)")
    ap.add_argument("--entry-order", action="append", metavar="CHR:POS=I[:D],I[:D],...", help="spliceai shapes: the site's records joined into one record per allele whose entries are the listed ones in that order, I the place of a symbol among the allele's records and D the draw of its scores (the record's own without D) (repeatable)")
    ap.add_argument("--repeat-lines", action="append", metavar="CHR:POS=N", help="spliceai shapes: every record line of the site, in its plain or its other form, written N times next to each other (repeatable)")
    ap.add_argument("--rescore-lines", action="append", metavar="CHR:POS=N", help="spliceai shapes: every record line of the site written N times next to each other, each copy behind the first with its entries' scores drawn again (repeatable)")
    ap.add_argument("--second-stream-at", metavar="CHR:POS", help="the file's first row at this site opens a second bgzip stream joined behind the first")
    ap.add_argument("--assembly", choices=CADD_ASSEMBLIES, help="cadd shapes: the build the first line names (cadd_sv is GRCh38 alone)")
    ap.add_argument("--exon-flank", type=int, metavar="N", help=f"loftee_gerp: the scored bases kept either side of every exon (default: {DEFAULT_EXON_FLANK})")
    ap.add_argument("--force", action="store_true", help="overwrite an existing fixture file")
    args = ap.parse_args(argv)

    corpus = Path(args.corpus)
    vcf = corpus / "variants.vcf"
    if not vcf.is_file():
        raise fail(f"corpus has no variants.vcf: {corpus}")
    if "/" in args.name or not args.name:
        raise fail(f"--name is a file name, not a path: {args.name!r}")
    if "/" in args.dir or args.dir in ("", ".."):
        raise fail(f"--dir is a directory name inside the corpus, not a path: {args.dir!r}")
    for miss in args.miss or []:
        if not re.fullmatch(r"[^:\s]+:[0-9]+", miss):
            raise fail(f"--miss names a site as chr:pos, got {miss!r}")
    if (args.shape == "splicevault") != bool(args.transcripts):
        raise fail("--transcripts is given with --shape splicevault and with no other shape")
    if (args.shape in SPLICEAI_SHAPES) != bool(args.symbols):
        raise fail("--symbols is given with --shape spliceai_snv or spliceai_indel and with no other shape")
    for option, given in (("--transcripts", args.transcripts), ("--symbols", args.symbols)):
        if given and not Path(given).is_file():
            raise fail(f"{option} file not found: {given}")
    if (args.shape in CADD_SHAPES) != bool(args.assembly):
        raise fail("--assembly is given with a cadd shape and with no other shape")
    if args.shape == "cadd_sv" and args.assembly != "GRCh38":
        raise fail("cadd_sv is written for GRCh38 alone")
    cache = corpus / "json_cache"
    if args.shape == LOFTEE_GERP_SHAPE:
        if args.exon_flank is None:
            args.exon_flank = DEFAULT_EXON_FLANK
        if args.exon_flank < 0:
            raise fail(f"--exon-flank is a count of bases, got {args.exon_flank}")
        if not (cache / "transcripts").is_dir():
            raise fail(f"corpus has no json_cache/transcripts, whose exons the loftee_gerp shape scores: {corpus}")
    elif args.exon_flank is not None:
        raise fail("--exon-flank is given with --shape loftee_gerp and with no other shape")
    header_names = header_names_for(args.shape, args.header_names) if args.header_names is not None else None
    forms = record_forms_for(args.shape, args.join_entries or [], args.no_information or [], args.join_alternates or [], args.any_base_reference or [], args.repeat_lines or [], args.rescore_lines or [], args.entry_order or [])
    out = corpus / args.dir / args.name
    if out.exists() and not args.force:
        raise fail(f"{out} exists; pass --force to overwrite it")
    versions = {tool: tool_version(tool) for tool in TOOLS}

    misses = set(args.miss or [])
    if args.shape == "splicevault":
        rows, counts = splicevault_rows(read_records(vcf), read_transcripts(Path(args.transcripts)), args.seed, args.chr_style, misses)
    elif args.shape in SPLICEAI_SHAPES:
        rows, counts = spliceai_rows(read_records(vcf), read_symbols(Path(args.symbols)), args.shape, args.seed, args.chr_style, misses, forms)
    elif args.shape == "cadd_snv":
        rows, counts = cadd_snv_rows(read_alleles(vcf), args.seed, args.chr_style, misses)
    elif args.shape == "cadd_indel":
        rows, counts = cadd_indel_rows(read_records(vcf), args.seed, args.chr_style, misses)
    elif args.shape == "cadd_sv":
        rows, counts = cadd_sv_rows(read_structural(vcf), args.seed, args.chr_style, misses)
    elif args.shape == LOFTEE_GERP_SHAPE:
        rows, counts = loftee_gerp_rows(read_exon_regions(cache, args.exon_flank), args.seed, args.chr_style, misses)
        if not rows:
            raise fail(f"{cache}: no exon base to score (no cached exon, or every base missed)")
    else:
        rows, counts = rows_for(read_alleles(vcf), args.shape, args.seed, args.chr_style, misses)
    if not rows:
        raise fail(f"{vcf}: no site to score (no sequence allele, or every site missed)")
    if args.mitochondrial_name:
        rows = with_mitochondrial_name(rows, args.mitochondrial_name)
    lines = [*headers_for(args.shape, args.assembly, header_names), *rows]
    second_stream = second_stream_line(lines, args.second_stream_at) if args.second_stream_at else None
    write_fixture(lines, out, INDEX_ARGS.get(args.shape), second_stream)
    update_provenance(corpus, out, args, counts, versions)
    summary = ", ".join(f"{'missed sites' if key == 'miss' else key} {value}" for key, value in counts.items())
    print(f"wrote {out.relative_to(corpus)} and its .tbi: {len(rows)} rows ({summary}), seed {args.seed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
