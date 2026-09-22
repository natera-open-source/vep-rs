#!/usr/bin/env python3
"""Generate 13 synthetic VCF files covering 41 VCF variant categories.

Produces 15,400 variants per assembly across chromosome 21 gene regions, with
REF alleles read from the reference FASTA, for vep-rs structural variant
classification testing. Generation is seeded (``RANDOM_SEED``), so the same
FASTA and assembly reproduce the same files.

Only chromosome 21 is loaded from the FASTA. The mate record of an
inter-chromosomal breakend (``09_breakends.vcf``, chromosomes 1 to 5) therefore
carries a placeholder reference base, ``A``, rather than the base at that
position; every other REF, including the anchor base of a symbolic allele, is
the reference sequence.

Usage:
    python3 scripts/validation/generate_sv_test_vcfs.py \
        --fasta /path/to/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
        --output-dir /path/to/output/
"""

import argparse
import bisect
import os
import random
import re
import subprocess
import sys
from dataclasses import dataclass

# Constants

CONTIG = "21"
RANDOM_SEED = 42

# Contig lengths by assembly (for VCF header)
CONTIG_LENGTHS_BY_ASSEMBLY = {
    "GRCh37": {
        "21": 48129895,
        "1": 249250621,
        "2": 243199373,
        "3": 198022430,
        "4": 191154276,
        "5": 180915260,
    },
    "GRCh38": {
        "21": 46709983,
        "1": 248956422,
        "2": 242193529,
        "3": 198295559,
        "4": 190214555,
        "5": 181538259,
    },
}

# Gene regions on chromosome 21 by assembly
GENE_REGIONS_BY_ASSEMBLY = {
    "GRCh37": {
        "SOD1": (33031935, 33041243),
        "APP": (27252861, 27543138),
        "DYRK1A": (38739398, 38887511),
        "KCNJ6": (38348088, 38655076),
        "DSCAM": (41382995, 42217961),
    },
    "GRCh38": {
        "SOD1": (31659633, 31668931),
        "APP": (25880535, 26171128),
        "DYRK1A": (37365573, 37526358),
        "KCNJ6": (37607373, 38121345),
        "DSCAM": (40010999, 40847158),
    },
}

# Module-level default; main() rebinds it to the selected assembly.
GENE_REGIONS = GENE_REGIONS_BY_ASSEMBLY["GRCh37"]

# Bases for generating random sequences
BASES = "ACGT"

# Tandem repeats placed by file 07: a primitive unit of 2 to 6 bases repeated perfectly at
# least five times in the reference. The alternate allele is kept below 5,000 bases, the
# fallback expansion cap in `Parser/VCF.pm` (`param('max_sv_size') || 5000`), so every
# record's allele is written out on the reference rather than kept symbolic.
TR_UNIT_LENGTHS = range(2, 7)
TR_MIN_COPIES = 5
TR_MAX_ALT_BASES = 4000

# Current assembly (set by main())
ASSEMBLY = "GRCh37"


def vcf_meta(fileformat: str = "VCFv4.3") -> str:
    """Build assembly-aware VCF meta-header."""
    contigs = CONTIG_LENGTHS_BY_ASSEMBLY[ASSEMBLY]
    lines = [
        f"##fileformat={fileformat}",
        "##source=generate_sv_test_vcfs.py",
        f"##reference={ASSEMBLY}",
    ]
    for cid, clen in contigs.items():
        lines.append(f"##contig=<ID={cid},length={clen}>")
    return "\n".join(lines)


# INFO field definitions used across the generated files
INFO_DEFS = {
    "END": '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
    "SVTYPE": '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    "SVLEN": '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    "CIPOS": '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">',
    "CIEND": '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">',
    "IMPRECISE": '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
    "MATEID": '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">',
    "EVENT": '##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">',
    "CN": '##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number of segment">',
    "RN": '##INFO=<ID=RN,Number=A,Type=Integer,Description="Total number of repeat sequences in this allele">',
    "RUS": '##INFO=<ID=RUS,Number=.,Type=String,Description="Repeat unit sequence of the corresponding repeat sequence">',
    "RUL": '##INFO=<ID=RUL,Number=.,Type=Integer,Description="Repeat unit length of the corresponding repeat sequence">',
    "RUC": '##INFO=<ID=RUC,Number=.,Type=Float,Description="Repeat unit count of corresponding repeat sequence">',
    "RB": '##INFO=<ID=RB,Number=.,Type=Integer,Description="Total number of bases in the corresponding repeat sequence">',
    "MEINFO": '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info: name,start,end,polarity">',
    "LEN": '##INFO=<ID=LEN,Number=1,Type=Integer,Description="Reference block length">',
    "SVCLAIM": '##INFO=<ID=SVCLAIM,Number=1,Type=String,Description="SV claim type: D=abundance, J=adjacency, DJ=both">',
    "PARID": '##INFO=<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">',
}

# ALT definitions
ALT_DEFS = {
    "DEL": '##ALT=<ID=DEL,Description="Deletion">',
    "INS": '##ALT=<ID=INS,Description="Insertion">',
    "DUP": '##ALT=<ID=DUP,Description="Duplication">',
    "DUP:TANDEM": '##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">',
    "DUP:ISP": '##ALT=<ID=DUP:ISP,Description="Interspersed duplication">',
    "INV": '##ALT=<ID=INV,Description="Inversion">',
    "CNV": '##ALT=<ID=CNV,Description="Copy number variable region">',
    "CNV:TR": '##ALT=<ID=CNV:TR,Description="Tandem repeat copy number variant">',
    "INS:ME": '##ALT=<ID=INS:ME,Description="Mobile element insertion">',
    "INS:ME:LINE": '##ALT=<ID=INS:ME:LINE,Description="LINE mobile element insertion">',
    "INS:ME:ALU": '##ALT=<ID=INS:ME:ALU,Description="Alu mobile element insertion">',
    "INS:ME:SVA": '##ALT=<ID=INS:ME:SVA,Description="SVA mobile element insertion">',
    "DEL:ME": '##ALT=<ID=DEL:ME,Description="Mobile element deletion">',
    "CPX": '##ALT=<ID=CPX,Description="Complex structural variant">',
    "NON_REF": '##ALT=<ID=NON_REF,Description="Non-reference allele (gVCF)">',
}

# CN0-CN9 ALT defs
for i in range(10):
    ALT_DEFS[f"CN{i}"] = f'##ALT=<ID=CN{i},Description="Copy number {i}">'


# Helpers


@dataclass
class VcfRecord:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: str = "."
    filt: str = "PASS"
    info: str = "."

    def to_line(self) -> str:
        return f"{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}\t{self.qual}\t{self.filt}\t{self.info}"


class RefFasta:
    """Fetch reference bases with full-chromosome pre-loading.

    Loads the entire primary contig (chr21, ~48MB) into memory once so all
    lookups are O(1) string slicing with zero subprocess overhead.
    """

    def __init__(self, fasta_path: str, contig: str = CONTIG):
        self.fasta_path = fasta_path
        self._seq: str = ""
        self._contig = contig
        self._load_contig(contig)

    def _load_contig(self, contig: str) -> None:
        """Load entire contig into memory via single samtools call."""
        result = subprocess.run(
            ["samtools", "faidx", self.fasta_path, contig],
            capture_output=True,
            text=True,
            check=True,
        )
        self._seq = "".join(result.stdout.strip().split("\n")[1:]).upper()
        print(f"  Loaded {contig}: {len(self._seq):,} bp")

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch bases from chrom:start-end (1-based inclusive)."""
        if chrom == self._contig and 1 <= start and end <= len(self._seq):
            return self._seq[start - 1 : end]
        # Fallback for other contigs
        region = f"{chrom}:{start}-{end}"
        result = subprocess.run(
            ["samtools", "faidx", self.fasta_path, region],
            capture_output=True,
            text=True,
            check=True,
        )
        return "".join(result.stdout.strip().split("\n")[1:]).upper()

    def fetch_base(self, chrom: str, pos: int) -> str:
        return self.fetch(chrom, pos, pos)


def random_seq(length: int) -> str:
    return "".join(random.choice(BASES) for _ in range(length))


def is_primitive_unit(unit: str) -> bool:
    """False when the unit is itself a repetition of a shorter unit (ATAT of AT, AAA of A)."""
    return not any(
        len(unit) % k == 0 and unit == unit[:k] * (len(unit) // k)
        for k in range(1, len(unit))
    )


def find_tandem_repeats(
    ref: RefFasta, regions: dict[str, tuple[int, int]]
) -> list[tuple[int, str, int]]:
    """Every maximal perfect tandem repeat inside `regions`, as (first base, unit, copies).

    A run is a primitive unit of `TR_UNIT_LENGTHS` bases repeated at least `TR_MIN_COPIES`
    times without interruption; the first base is 1-based, a partial trailing copy is not
    counted, and a run found under two unit lengths is kept under the shorter one.
    """
    runs: dict[int, tuple[str, int]] = {}
    for start, end in sorted(regions.values()):
        seq = ref.fetch(CONTIG, start, end)
        for k in TR_UNIT_LENGTHS:
            pattern = re.compile(r"([ACGT]{%d})\1{%d,}" % (k, TR_MIN_COPIES - 1))
            for m in pattern.finditer(seq):
                unit = m.group(1)
                if not is_primitive_unit(unit):
                    continue
                runs.setdefault(start + m.start(), (unit, len(m.group(0)) // k))
    return sorted((pos, unit, copies) for pos, (unit, copies) in runs.items())


def mutate_base(base: str) -> str:
    """Return a different base."""
    alts = [b for b in BASES if b != base]
    return random.choice(alts)


def spread_positions(regions: dict[str, tuple[int, int]], count: int) -> list[int]:
    """Generate sorted, unique positions spread across gene regions."""
    positions = set()
    region_list = list(regions.values())
    while len(positions) < count:
        start, end = random.choice(region_list)
        # Stay within region, leave margin
        pos = random.randint(start + 10, end - 10)
        positions.add(pos)
    return sorted(positions)


def spread_positions_with_margin(
    regions: dict[str, tuple[int, int]], count: int, min_gap: int = 1
) -> list[int]:
    """Generate sorted unique positions with minimum gap between them."""
    positions = spread_positions(regions, count * 3)
    # Filter to maintain minimum gap
    result = [positions[0]]
    for p in positions[1:]:
        if p - result[-1] >= min_gap and len(result) < count:
            result.append(p)
    # Relax the gap constraint, with a retry limit, if the count falls short
    attempts = 0
    while len(result) < count and attempts < 50000:
        start, end = random.choice(list(regions.values()))
        pos = random.randint(start + 100, end - 100)
        if all(abs(pos - r) >= min_gap for r in result):
            result.append(pos)
        attempts += 1
    # If gap constraint is too tight, just fill with spread positions (no gap)
    if len(result) < count:
        extras = spread_positions(regions, count - len(result) + 100)
        for p in extras:
            if p not in result:
                result.append(p)
            if len(result) >= count:
                break
    return sorted(result[:count])


def build_vcf_header(
    info_keys: list[str],
    alt_keys: list[str],
    extra_meta: list[str] | None = None,
    fileformat: str = "VCFv4.3",
) -> str:
    """Build VCF header string with specified INFO and ALT definitions."""
    lines = [vcf_meta(fileformat)]
    for k in sorted(info_keys):
        if k in INFO_DEFS:
            lines.append(INFO_DEFS[k])
    for k in sorted(alt_keys):
        if k in ALT_DEFS:
            lines.append(ALT_DEFS[k])
    if extra_meta:
        lines.extend(extra_meta)
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return "\n".join(lines)


def write_vcf(path: str, header: str, records: list[VcfRecord]) -> int:
    """Write VCF file. Returns record count."""
    # Sort by position
    records.sort(key=lambda r: (r.chrom, r.pos))
    with open(path, "w") as f:
        f.write(header + "\n")
        for rec in records:
            f.write(rec.to_line() + "\n")
    return len(records)


def validate_vcf(path: str) -> bool:
    """Validate VCF with bcftools. Returns True if valid."""
    result = subprocess.run(
        ["bcftools", "view", "--no-version", "-H", path],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"  WARNING: bcftools validation failed for {path}")
        print(f"  stderr: {result.stderr[:500]}")
        return False
    return True


# File generators


def gen_01_snv(ref: RefFasta, output_dir: str) -> tuple[str, int, list[str]]:
    """File 1: 1,000 SNVs across all gene regions."""
    fname = "01_snv.vcf"
    path = os.path.join(output_dir, fname)
    count = 1000

    positions = spread_positions(GENE_REGIONS, count)
    records = []
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        alt_base = mutate_base(ref_base)
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_snv_{i+1:04d}",
                ref=ref_base,
                alt=alt_base,
            )
        )

    header = build_vcf_header([], [])
    n = write_vcf(path, header, records)
    return fname, n, ["Cat1:SNV"]


def gen_02_mnp_complex(ref: RefFasta, output_dir: str) -> tuple[str, int, list[str]]:
    """File 2: 1,000 MNPs + 1,000 complex substitutions."""
    fname = "02_mnp_complex.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    # Category 2: MNPs (equal ref/alt length, 2-5bp)
    positions = spread_positions_with_margin(GENE_REGIONS, 1000, min_gap=10)
    for i, pos in enumerate(positions):
        length = random.randint(2, 5)
        ref_seq = ref.fetch(CONTIG, pos, pos + length - 1)
        if "N" in ref_seq:
            ref_seq = random_seq(length)
        # Mutate at least 2 positions
        alt_list = list(ref_seq)
        mut_positions = random.sample(range(length), min(2, length))
        for mp in mut_positions:
            alt_list[mp] = mutate_base(alt_list[mp])
        alt_seq = "".join(alt_list)
        if alt_seq == ref_seq:
            alt_list[0] = mutate_base(alt_list[0])
            alt_seq = "".join(alt_list)
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_mnp_{i+1:04d}",
                ref=ref_seq,
                alt=alt_seq,
            )
        )

    # Category 5: Complex substitutions (different ref/alt length, both > 1bp)
    positions = spread_positions_with_margin(GENE_REGIONS, 1000, min_gap=10)
    for i, pos in enumerate(positions):
        ref_len = random.randint(2, 6)
        alt_len = random.randint(2, 6)
        while alt_len == ref_len:
            alt_len = random.randint(2, 6)
        ref_seq = ref.fetch(CONTIG, pos, pos + ref_len - 1)
        if "N" in ref_seq:
            ref_seq = random_seq(ref_len)
        alt_seq = random_seq(alt_len)
        # Make sure alt differs
        while alt_seq == ref_seq:
            alt_seq = random_seq(alt_len)
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_cpx_sub_{i+1:04d}",
                ref=ref_seq,
                alt=alt_seq,
            )
        )

    header = build_vcf_header([], [])
    n = write_vcf(path, header, records)
    return fname, n, ["Cat2:MNP", "Cat5:ComplexSub"]


def gen_03_small_indels(ref: RefFasta, output_dir: str) -> tuple[str, int, list[str]]:
    """File 3: 1,000 small deletions + 1,000 small insertions (1-50bp)."""
    fname = "03_small_indels.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    # Category 3: Small deletions
    positions = spread_positions_with_margin(GENE_REGIONS, 1000, min_gap=60)
    for i, pos in enumerate(positions):
        del_len = random.randint(1, 50)
        # VCF convention: include padding base before deletion
        ref_seq = ref.fetch(CONTIG, pos, pos + del_len)
        if "N" in ref_seq:
            ref_seq = random_seq(del_len + 1)
        alt_seq = ref_seq[0]  # padding base only
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_del_sm_{i+1:04d}",
                ref=ref_seq,
                alt=alt_seq,
            )
        )

    # Category 4: Small insertions
    positions = spread_positions_with_margin(GENE_REGIONS, 1000, min_gap=10)
    for i, pos in enumerate(positions):
        ins_len = random.randint(1, 50)
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        alt_seq = ref_base + random_seq(ins_len)
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_ins_sm_{i+1:04d}",
                ref=ref_base,
                alt=alt_seq,
            )
        )

    header = build_vcf_header([], [])
    n = write_vcf(path, header, records)
    return fname, n, ["Cat3:SmallDel", "Cat4:SmallIns"]


def gen_04_large_explicit(ref: RefFasta, output_dir: str) -> tuple[str, int, list[str]]:
    """File 4: 200 large explicit deletions + 200 large explicit insertions."""
    fname = "04_large_explicit.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    # Category 6: Large explicit deletions (50bp - 10kb)
    # Use wider spacing since these are large
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=500)
    for i, pos in enumerate(positions):
        del_len = random.randint(50, 10000)
        # Fetch real reference sequence (including padding base)
        ref_seq = ref.fetch(CONTIG, pos, pos + del_len)
        if len(ref_seq) < del_len:
            continue  # skip if near contig boundary
        alt_seq = ref_seq[0]
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_del_lg_{i+1:04d}",
                ref=ref_seq,
                alt=alt_seq,
            )
        )

    # Category 7: Large explicit insertions (50bp - 1kb)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=100)
    for i, pos in enumerate(positions):
        ins_len = random.randint(50, 1000)
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        alt_seq = ref_base + random_seq(ins_len)
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_ins_lg_{i+1:04d}",
                ref=ref_base,
                alt=alt_seq,
            )
        )

    header = build_vcf_header([], [])
    n = write_vcf(path, header, records)
    return fname, n, ["Cat6:LargeExplicitDel", "Cat7:LargeExplicitIns"]


def gen_05_symbolic_del_ins(
    ref: RefFasta, output_dir: str
) -> tuple[str, int, list[str]]:
    """File 5: 1,000 symbolic DEL + 1,000 symbolic INS."""
    fname = "05_symbolic_del_ins.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    # Category 8: Symbolic <DEL> (1kb - 1Mb)
    positions = spread_positions_with_margin(GENE_REGIONS, 1000, min_gap=100)
    for i, pos in enumerate(positions):
        sv_len = random.randint(1000, 1000000)
        end = pos + sv_len
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_del_sym_{i+1:04d}",
                ref=ref_base,
                alt="<DEL>",
                info=f"SVTYPE=DEL;END={end};SVLEN=-{sv_len}",
            )
        )

    # Category 9: Symbolic <INS> (varying sizes)
    positions = spread_positions_with_margin(GENE_REGIONS, 1000, min_gap=100)
    for i, pos in enumerate(positions):
        sv_len = random.randint(50, 500000)
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_ins_sym_{i+1:04d}",
                ref=ref_base,
                alt="<INS>",
                info=f"SVTYPE=INS;SVLEN={sv_len}",
            )
        )

    header = build_vcf_header(["END", "SVTYPE", "SVLEN"], ["DEL", "INS"])
    n = write_vcf(path, header, records)
    return fname, n, ["Cat8:SymbolicDEL", "Cat9:SymbolicINS"]


def gen_06_symbolic_dup_inv(
    ref: RefFasta, output_dir: str
) -> tuple[str, int, list[str]]:
    """File 6: DUP, DUP:TANDEM, DUP:ISP, INV."""
    fname = "06_symbolic_dup_inv.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    # Category 10: <DUP> (1,000)
    positions = spread_positions_with_margin(GENE_REGIONS, 1000, min_gap=100)
    for i, pos in enumerate(positions):
        sv_len = random.randint(500, 500000)
        end = pos + sv_len
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_dup_{i+1:04d}",
                ref=ref_base,
                alt="<DUP>",
                info=f"SVTYPE=DUP;END={end};SVLEN={sv_len}",
            )
        )

    # Category 11: <DUP:TANDEM> (200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=500)
    for i, pos in enumerate(positions):
        sv_len = random.randint(100, 100000)
        end = pos + sv_len
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_dup_tan_{i+1:04d}",
                ref=ref_base,
                alt="<DUP:TANDEM>",
                info=f"SVTYPE=DUP;END={end};SVLEN={sv_len}",
            )
        )

    # Category 12: <DUP:ISP> (50)
    positions = spread_positions_with_margin(GENE_REGIONS, 50, min_gap=1000)
    for i, pos in enumerate(positions):
        sv_len = random.randint(500, 50000)
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_dup_isp_{i+1:04d}",
                ref=ref_base,
                alt="<DUP:ISP>",
                info=f"SVTYPE=DUP;SVLEN={sv_len}",
            )
        )

    # Category 13: <INV> (1,000)
    positions = spread_positions_with_margin(GENE_REGIONS, 1000, min_gap=100)
    for i, pos in enumerate(positions):
        sv_len = random.randint(500, 500000)
        end = pos + sv_len
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_inv_{i+1:04d}",
                ref=ref_base,
                alt="<INV>",
                info=f"SVTYPE=INV;END={end};SVLEN={sv_len}",
            )
        )

    header = build_vcf_header(
        ["END", "SVTYPE", "SVLEN"],
        ["DUP", "DUP:TANDEM", "DUP:ISP", "INV"],
    )
    n = write_vcf(path, header, records)
    return (
        fname,
        n,
        ["Cat10:DUP", "Cat11:DUP:TANDEM", "Cat12:DUP:ISP", "Cat13:INV"],
    )


def gen_07_cnv_repeat(ref: RefFasta, output_dir: str) -> tuple[str, int, list[str]]:
    """File 7: CNV, CN0-CN9, CNV:TR."""
    fname = "07_cnv_repeat.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    # Category 14: <CNV> (500)
    positions = spread_positions_with_margin(GENE_REGIONS, 500, min_gap=200)
    for i, pos in enumerate(positions):
        sv_len = random.randint(1000, 500000)
        end = pos + sv_len
        cn = random.randint(0, 10)
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_cnv_{i+1:04d}",
                ref=ref_base,
                alt="<CNV>",
                info=f"SVTYPE=CNV;END={end};CN={cn}",
            )
        )

    # Category 15: <CN0> through <CN9> (20 each = 200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=200)
    for i, pos in enumerate(positions):
        cn_val = i % 10  # cycle through 0-9
        sv_len = random.randint(5000, 200000)
        end = pos + sv_len
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_cn{cn_val}_{(i // 10)+1:04d}",
                ref=ref_base,
                alt=f"<CN{cn_val}>",
                info=f"SVTYPE=CNV;END={end};CN={cn_val}",
            )
        )

    # Category 16: <CNV:TR> (200), each on a tandem repeat the reference carries. POS is
    # the base before the run, END its last base and SVLEN its length, so the reference
    # allele VEP reads over SVLEN is the run itself; RUS is the run's unit and RUC the
    # alternate copy count, above the reference count on even records (an expansion) and
    # below it on odd ones (a contraction), with RB = RUC x RUL.
    #
    # Files 08 to 13 draw from the same generator after this category, so three
    # per-record draws (a choice among eight, a uniform and a randint) are consumed and
    # discarded here, and the categories after this one see the same random stream
    # whatever this category draws; the repeat-specific choices come from a generator
    # of their own.
    runs = find_tandem_repeats(ref, GENE_REGIONS)
    run_positions = [r[0] for r in runs]
    tr_rng = random.Random(RANDOM_SEED * 1000 + 16)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=200)
    chosen: list[tuple[int, str, int]] = []
    for pos in positions:
        random.randrange(8)
        random.uniform(5.0, 200.0)
        random.randint(100, 5000)
        # The nearest run to the drawn position at least 200 bases from every run chosen,
        # whose reference count leaves room for an expansion under the allele cap.
        k = bisect.bisect_left(run_positions, pos)
        lo, hi = k - 1, k
        while lo >= 0 or hi < len(runs):
            if hi >= len(runs) or (lo >= 0 and pos - run_positions[lo] <= run_positions[hi] - pos):
                cand, lo = runs[lo], lo - 1
            else:
                cand, hi = runs[hi], hi + 1
            run_pos, unit, copies = cand
            if copies + 1 > TR_MAX_ALT_BASES // len(unit):
                continue
            if all(abs(run_pos - c[0]) >= 200 for c in chosen):
                chosen.append(cand)
                break
        else:
            raise SystemExit(
                "ERROR: [generate_sv_test_vcfs] no tandem repeat left near position "
                f"{pos} for a <CNV:TR> record"
            )
    for i, (run_pos, unit, copies) in enumerate(sorted(chosen)):
        rul = len(unit)
        sv_len = copies * rul
        pos = run_pos - 1
        end = pos + sv_len
        if i % 2 == 0:
            ruc = copies + tr_rng.randint(1, min(200, TR_MAX_ALT_BASES // rul - copies))
        else:
            ruc = copies - tr_rng.randint(1, copies - 1)
        rb = ruc * rul
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_cnv_tr_{i+1:04d}",
                ref=ref_base,
                alt="<CNV:TR>",
                info=f"SVTYPE=CNV;END={end};SVLEN={sv_len};RN=1;RUS={unit};RUL={rul};RUC={ruc};RB={rb}",
            )
        )

    cn_alts = [f"CN{i}" for i in range(10)]
    header = build_vcf_header(
        ["END", "SVTYPE", "SVLEN", "CN", "RN", "RUS", "RUL", "RUC", "RB"],
        ["CNV", "CNV:TR"] + cn_alts,
    )
    n = write_vcf(path, header, records)
    return fname, n, ["Cat14:CNV", "Cat15:CN0-CN9", "Cat16:CNV:TR"]


def gen_08_mobile_elements(
    ref: RefFasta, output_dir: str
) -> tuple[str, int, list[str]]:
    """File 8: Mobile element insertions and deletions."""
    fname = "08_mobile_elements.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    def make_meinfo(name: str) -> str:
        start = random.randint(0, 500)
        end = start + random.randint(100, 6000)
        polarity = random.choice(["+", "-"])
        return f"MEINFO={name},{start},{end},{polarity}"

    # Category 17: <INS:ME> generic (200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=100)
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_ins_me_{i+1:04d}",
                ref=ref_base,
                alt="<INS:ME>",
                info=f"SVTYPE=INS;{make_meinfo('UNKNOWN')}",
            )
        )

    # Category 18: <INS:ME:LINE> (200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=100)
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        line_type = random.choice(["L1", "L1HS", "L1PA2", "L1PA3"])
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_ins_line_{i+1:04d}",
                ref=ref_base,
                alt="<INS:ME:LINE>",
                info=f"SVTYPE=INS;{make_meinfo(line_type)}",
            )
        )

    # Category 19: <INS:ME:ALU> (200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=100)
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        alu_type = random.choice(["AluYa5", "AluYb8", "AluY", "AluSx"])
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_ins_alu_{i+1:04d}",
                ref=ref_base,
                alt="<INS:ME:ALU>",
                info=f"SVTYPE=INS;{make_meinfo(alu_type)}",
            )
        )

    # Category 20: <INS:ME:SVA> (100)
    positions = spread_positions_with_margin(GENE_REGIONS, 100, min_gap=200)
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_ins_sva_{i+1:04d}",
                ref=ref_base,
                alt="<INS:ME:SVA>",
                info=f"SVTYPE=INS;{make_meinfo('SVA_E')}",
            )
        )

    # Category 21: <DEL:ME> (100)
    positions = spread_positions_with_margin(GENE_REGIONS, 100, min_gap=200)
    for i, pos in enumerate(positions):
        sv_len = random.randint(100, 6000)
        end = pos + sv_len
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        me_name = random.choice(["L1", "AluY", "SVA_D", "L1HS"])
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_del_me_{i+1:04d}",
                ref=ref_base,
                alt="<DEL:ME>",
                info=f"SVTYPE=DEL;END={end};SVLEN=-{sv_len};{make_meinfo(me_name)}",
            )
        )

    header = build_vcf_header(
        ["END", "SVTYPE", "SVLEN", "MEINFO"],
        ["INS:ME", "INS:ME:LINE", "INS:ME:ALU", "INS:ME:SVA", "DEL:ME"],
    )
    n = write_vcf(path, header, records)
    return (
        fname,
        n,
        [
            "Cat17:INS:ME",
            "Cat18:INS:ME:LINE",
            "Cat19:INS:ME:ALU",
            "Cat20:INS:ME:SVA",
            "Cat21:DEL:ME",
        ],
    )


def gen_09_breakends(ref: RefFasta, output_dir: str) -> tuple[str, int, list[str]]:
    """File 9: All BND orientations, unmatched, with inserted seq, inter-chromosomal."""
    fname = "09_breakends.vcf"
    path = os.path.join(output_dir, fname)

    records = []
    other_chroms = ["1", "2", "3", "4", "5"]

    def get_ref_base(pos: int) -> str:
        b = ref.fetch_base(CONTIG, pos)
        return b if b != "N" else "A"

    # Helper: create a BND pair with proper MATEID linkage
    def make_bnd_pair(
        idx: int,
        prefix: str,
        pos_a: int,
        pos_b: int,
        chrom_a: str,
        chrom_b: str,
        orientation: str,
        inserted: str = "",
    ) -> list[VcfRecord]:
        """Create a BND pair.
        orientation: 'ff' (t[p[), 'fr' (t]p]), 'rf' (]p]t), 'rr' ([p[t)
        """
        id_a = f"{prefix}_{idx:04d}_A"
        id_b = f"{prefix}_{idx:04d}_B"
        # A mate off the loaded contig gets the placeholder base described in the
        # module docstring; an N at a loaded position gets the same base.
        ref_a = ref.fetch_base(chrom_a, pos_a) if chrom_a == CONTIG else "N"
        ref_b = ref.fetch_base(chrom_b, pos_b) if chrom_b == CONTIG else "N"
        if ref_a == "N":
            ref_a = "A"
        if ref_b == "N":
            ref_b = "A"

        t = inserted  # inserted sequence between breakpoint and join

        if orientation == "ff":
            # Forward-after: t[chr:pos[
            alt_a = f"{ref_a}{t}[{chrom_b}:{pos_b}["
            alt_b = f"{ref_b}{t}[{chrom_a}:{pos_a}["
        elif orientation == "fr":
            # Reverse-after: t]chr:pos]
            alt_a = f"{ref_a}{t}]{chrom_b}:{pos_b}]"
            alt_b = f"[{chrom_a}:{pos_a}[{t}{ref_b}"
        elif orientation == "rf":
            # Reverse-before: ]chr:pos]t
            alt_a = f"]{chrom_b}:{pos_b}]{t}{ref_a}"
            alt_b = f"]{chrom_a}:{pos_a}]{t}{ref_b}"
        elif orientation == "rr":
            # Forward-before: [chr:pos[t
            alt_a = f"[{chrom_b}:{pos_b}[{t}{ref_a}"
            alt_b = f"{ref_b}{t}]{chrom_a}:{pos_a}]"
        else:
            raise ValueError(f"Unknown orientation: {orientation}")

        return [
            VcfRecord(
                chrom=chrom_a,
                pos=pos_a,
                id=id_a,
                ref=ref_a,
                alt=alt_a,
                info=f"SVTYPE=BND;MATEID={id_b}",
            ),
            VcfRecord(
                chrom=chrom_b,
                pos=pos_b,
                id=id_b,
                ref=ref_b,
                alt=alt_b,
                info=f"SVTYPE=BND;MATEID={id_a}",
            ),
        ]

    # Category 22: Forward-after t[chr:pos[ (250 variants = 125 pairs)
    positions = spread_positions_with_margin(GENE_REGIONS, 250, min_gap=50)
    for i in range(0, 250, 2):
        if i + 1 >= len(positions):
            break
        pair = make_bnd_pair(
            i // 2 + 1,
            "synth_bnd_ff",
            positions[i],
            positions[i + 1],
            CONTIG,
            CONTIG,
            "ff",
        )
        records.extend(pair)

    # Category 23: Reverse-after t]chr:pos] (250 = 125 pairs)
    positions = spread_positions_with_margin(GENE_REGIONS, 250, min_gap=50)
    for i in range(0, 250, 2):
        if i + 1 >= len(positions):
            break
        pair = make_bnd_pair(
            i // 2 + 1,
            "synth_bnd_fr",
            positions[i],
            positions[i + 1],
            CONTIG,
            CONTIG,
            "fr",
        )
        records.extend(pair)

    # Category 24: Reverse-before ]chr:pos]t (250 = 125 pairs)
    positions = spread_positions_with_margin(GENE_REGIONS, 250, min_gap=50)
    for i in range(0, 250, 2):
        if i + 1 >= len(positions):
            break
        pair = make_bnd_pair(
            i // 2 + 1,
            "synth_bnd_rf",
            positions[i],
            positions[i + 1],
            CONTIG,
            CONTIG,
            "rf",
        )
        records.extend(pair)

    # Category 25: Forward-before [chr:pos[t (250 = 125 pairs)
    positions = spread_positions_with_margin(GENE_REGIONS, 250, min_gap=50)
    for i in range(0, 250, 2):
        if i + 1 >= len(positions):
            break
        pair = make_bnd_pair(
            i // 2 + 1,
            "synth_bnd_rr",
            positions[i],
            positions[i + 1],
            CONTIG,
            CONTIG,
            "rr",
        )
        records.extend(pair)

    # Category 26: Unmatched/single BND (100)
    positions = spread_positions_with_margin(GENE_REGIONS, 100, min_gap=50)
    for i, pos in enumerate(positions):
        ref_base = get_ref_base(pos)
        if i % 2 == 0:
            alt = f"{ref_base}."
        else:
            alt = f".{ref_base}"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_bnd_single_{i+1:04d}",
                ref=ref_base,
                alt=alt,
                info="SVTYPE=BND",
            )
        )

    # Category 27: BND with inserted sequence (100 = 50 pairs)
    positions = spread_positions_with_margin(GENE_REGIONS, 100, min_gap=50)
    for i in range(0, 100, 2):
        if i + 1 >= len(positions):
            break
        ins_seq = random_seq(random.randint(1, 20))
        orientation = random.choice(["ff", "fr", "rf", "rr"])
        pair = make_bnd_pair(
            i // 2 + 1,
            "synth_bnd_ins",
            positions[i],
            positions[i + 1],
            CONTIG,
            CONTIG,
            orientation,
            inserted=ins_seq,
        )
        records.extend(pair)

    # Category 28: Inter-chromosomal BND (200 = 100 pairs)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=50)
    for i in range(0, 200, 2):
        if i + 1 >= len(positions):
            break
        other_chrom = random.choice(other_chroms)
        other_pos = random.randint(10000000, 100000000)
        orientation = random.choice(["ff", "fr", "rf", "rr"])
        pair = make_bnd_pair(
            i // 2 + 1,
            "synth_bnd_inter",
            positions[i],
            other_pos,
            CONTIG,
            other_chrom,
            orientation,
        )
        records.extend(pair)

    header = build_vcf_header(
        ["SVTYPE", "MATEID", "EVENT"],
        [],
    )
    n = write_vcf(path, header, records)
    return (
        fname,
        n,
        [
            "Cat22:BND_ff",
            "Cat23:BND_fr",
            "Cat24:BND_rf",
            "Cat25:BND_rr",
            "Cat26:BND_single",
            "Cat27:BND_ins",
            "Cat28:BND_inter",
        ],
    )


def gen_10_special_alleles(
    ref: RefFasta, output_dir: str
) -> tuple[str, int, list[str]]:
    """File 10: Spanning deletions, no-variant, gVCF NON_REF, unspecified."""
    fname = "10_special_alleles.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    # Category 29: Spanning deletions ALT=* (500)
    positions = spread_positions_with_margin(GENE_REGIONS, 500, min_gap=10)
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_span_{i+1:04d}",
                ref=ref_base,
                alt="*",
            )
        )

    # Category 30: No-variant / ref-only ALT=. (200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=10)
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_refonly_{i+1:04d}",
                ref=ref_base,
                alt=".",
            )
        )

    # Category 31: gVCF <NON_REF> with END (200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=500)
    for i, pos in enumerate(positions):
        end = pos + random.randint(10, 500)
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_nonref_{i+1:04d}",
                ref=ref_base,
                alt="<NON_REF>",
                info=f"END={end}",
            )
        )

    # Category 32: Unspecified <*> with END (100)
    positions = spread_positions_with_margin(GENE_REGIONS, 100, min_gap=500)
    for i, pos in enumerate(positions):
        end = pos + random.randint(10, 1000)
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_unspec_{i+1:04d}",
                ref=ref_base,
                alt="<*>",
                info=f"END={end}",
            )
        )

    header = build_vcf_header(
        ["END"],
        ["NON_REF"],
        ['##ALT=<ID=*,Description="Unspecified alternate allele">'],
    )
    n = write_vcf(path, header, records)
    return (
        fname,
        n,
        ["Cat29:SpanDel", "Cat30:RefOnly", "Cat31:NON_REF", "Cat32:Unspecified"],
    )


def gen_11_multi_allelic(ref: RefFasta, output_dir: str) -> tuple[str, int, list[str]]:
    """File 11: Multi-allelic with mixed symbolic + explicit, and spanning."""
    fname = "11_multi_allelic.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    # Category 33: Mixed multi-allelic ALT=G,<DEL> (500)
    positions = spread_positions_with_margin(GENE_REGIONS, 500, min_gap=50)
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        snv_alt = mutate_base(ref_base)
        sv_len = random.randint(100, 50000)
        end = pos + sv_len
        # Mix different symbolic types
        sym_choices = ["<DEL>", "<DUP>", "<INV>", "<CNV>"]
        symbolic = random.choice(sym_choices)
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_multi_mix_{i+1:04d}",
                ref=ref_base,
                alt=f"{snv_alt},{symbolic}",
                info=f"END={end};SVTYPE=DEL;SVLEN=-{sv_len}",
            )
        )

    # Category 34: Multi-allelic with spanning ALT=G,* (200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=50)
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        snv_alt = mutate_base(ref_base)
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_multi_span_{i+1:04d}",
                ref=ref_base,
                alt=f"{snv_alt},*",
            )
        )

    header = build_vcf_header(
        ["END", "SVTYPE", "SVLEN"],
        ["DEL", "DUP", "INV", "CNV"],
        ['##ALT=<ID=*,Description="Unspecified alternate allele">'],
    )
    n = write_vcf(path, header, records)
    return fname, n, ["Cat33:MultiMixed", "Cat34:MultiSpanning"]


def gen_12_complex_imprecise(
    ref: RefFasta, output_dir: str
) -> tuple[str, int, list[str]]:
    """File 12: Complex SVs and imprecise breakpoints."""
    fname = "12_complex_imprecise.vcf"
    path = os.path.join(output_dir, fname)

    records = []

    # Category 35: Complex SVs <CPX> with EVENT field (100 = ~33 events of 3 records)
    positions = spread_positions_with_margin(GENE_REGIONS, 100, min_gap=200)
    event_idx = 0
    i = 0
    while i < len(positions) - 2 and event_idx < 34:
        event_id = f"cpx_event_{event_idx+1:04d}"
        # Each complex SV event has 3 records linked by EVENT
        for j in range(3):
            if i + j >= len(positions):
                break
            pos = positions[i + j]
            ref_base = ref.fetch_base(CONTIG, pos)
            if ref_base == "N":
                ref_base = "A"
            sv_len = random.randint(100, 50000)
            end = pos + sv_len
            records.append(
                VcfRecord(
                    chrom=CONTIG,
                    pos=pos,
                    id=f"synth_cpx_{event_idx+1:04d}_{j+1}",
                    ref=ref_base,
                    alt="<CPX>",
                    info=f"SVTYPE=CPX;END={end};SVLEN={sv_len};EVENT={event_id}",
                )
            )
        i += 3
        event_idx += 1

    # Fill up remaining to reach ~100 total CPX records
    extra_positions = spread_positions_with_margin(
        GENE_REGIONS,
        max(0, 100 - len([r for r in records if "cpx" in r.id])),
        min_gap=200,
    )
    for k, pos in enumerate(extra_positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        sv_len = random.randint(100, 50000)
        end = pos + sv_len
        event_id = f"cpx_event_extra_{k+1:04d}"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_cpx_extra_{k+1:04d}",
                ref=ref_base,
                alt="<CPX>",
                info=f"SVTYPE=CPX;END={end};SVLEN={sv_len};EVENT={event_id}",
            )
        )

    # Category 36: Imprecise breakpoints (200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=200)
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        sv_len = random.randint(1000, 500000)
        end = pos + sv_len
        cipos_left = -random.randint(10, 500)
        cipos_right = random.randint(10, 500)
        ciend_left = -random.randint(10, 500)
        ciend_right = random.randint(10, 500)
        sv_type = random.choice(["DEL", "DUP", "INV"])
        alt = f"<{sv_type}>"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_imprecise_{i+1:04d}",
                ref=ref_base,
                alt=alt,
                info=f"IMPRECISE;SVTYPE={sv_type};END={end};SVLEN=-{sv_len};CIPOS={cipos_left},{cipos_right};CIEND={ciend_left},{ciend_right}",
            )
        )

    header = build_vcf_header(
        ["END", "SVTYPE", "SVLEN", "CIPOS", "CIEND", "IMPRECISE", "EVENT"],
        ["CPX", "DEL", "DUP", "INV"],
    )
    n = write_vcf(path, header, records)
    return fname, n, ["Cat35:CPX", "Cat36:Imprecise"]


def gen_13_vcf45_features(ref: RefFasta, output_dir: str) -> tuple[str, int, list[str]]:
    """File 13: VCF 4.5 features -- reference blocks, SVCLAIM, PARID, high-arity, boundary."""
    fname = "13_vcf45_features.vcf"
    path = os.path.join(output_dir, fname)

    records = []
    contig_len = CONTIG_LENGTHS_BY_ASSEMBLY[ASSEMBLY][CONTIG]

    # Category 37: <*> reference blocks with LEN INFO (200)
    positions = spread_positions_with_margin(GENE_REGIONS, 200, min_gap=500)
    for i, pos in enumerate(positions):
        block_len = random.randint(50, 5000)
        end = pos + block_len
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_refblock_{i+1:04d}",
                ref=ref_base,
                alt="<*>",
                info=f"END={end};LEN={block_len}",
            )
        )

    # Category 38: SVCLAIM field (200)
    # SVCLAIM=D (100), SVCLAIM=J (50), SVCLAIM=DJ (50)
    svclaim_specs = [("D", 100), ("J", 50), ("DJ", 50)]
    claim_idx = 0
    for claim_type, count in svclaim_specs:
        positions = spread_positions_with_margin(GENE_REGIONS, count, min_gap=100)
        for i, pos in enumerate(positions):
            sv_kind = random.choice(["DEL", "DUP", "INS"])
            sv_len = random.randint(500, 200000)
            end = pos + sv_len
            ref_base = ref.fetch_base(CONTIG, pos)
            if ref_base == "N":
                ref_base = "A"
            info = f"SVTYPE={sv_kind};END={end};SVLEN={'-' if sv_kind == 'DEL' else ''}{sv_len};SVCLAIM={claim_type}"
            records.append(
                VcfRecord(
                    chrom=CONTIG,
                    pos=pos,
                    id=f"synth_svclaim_{claim_type.lower()}_{claim_idx+1:04d}",
                    ref=ref_base,
                    alt=f"<{sv_kind}>",
                    info=info,
                )
            )
            claim_idx += 1

    # Category 39: PARID field -- BND pairs using PARID instead of MATEID (100 = 50 pairs)
    positions = spread_positions_with_margin(GENE_REGIONS, 100, min_gap=50)
    orientations = ["ff", "fr", "rf", "rr"]
    for i in range(0, 100, 2):
        if i + 1 >= len(positions):
            break
        pos_a = positions[i]
        pos_b = positions[i + 1]
        pair_idx = i // 2 + 1
        id_a = f"synth_parid_{pair_idx:04d}_A"
        id_b = f"synth_parid_{pair_idx:04d}_B"
        ref_a = ref.fetch_base(CONTIG, pos_a)
        ref_b = ref.fetch_base(CONTIG, pos_b)
        if ref_a == "N":
            ref_a = "A"
        if ref_b == "N":
            ref_b = "A"

        orientation = random.choice(orientations)
        if orientation == "ff":
            alt_a = f"{ref_a}[{CONTIG}:{pos_b}["
            alt_b = f"{ref_b}[{CONTIG}:{pos_a}["
        elif orientation == "fr":
            alt_a = f"{ref_a}]{CONTIG}:{pos_b}]"
            alt_b = f"[{CONTIG}:{pos_a}[{ref_b}"
        elif orientation == "rf":
            alt_a = f"]{CONTIG}:{pos_b}]{ref_a}"
            alt_b = f"]{CONTIG}:{pos_a}]{ref_b}"
        else:  # rr
            alt_a = f"[{CONTIG}:{pos_b}[{ref_a}"
            alt_b = f"{ref_b}]{CONTIG}:{pos_a}]"

        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos_a,
                id=id_a,
                ref=ref_a,
                alt=alt_a,
                info=f"SVTYPE=BND;PARID={id_b}",
            )
        )
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos_b,
                id=id_b,
                ref=ref_b,
                alt=alt_b,
                info=f"SVTYPE=BND;PARID={id_a}",
            )
        )

    # Category 40: High-arity multi-allelic (100 records with 5-8 ALT alleles)
    positions = spread_positions_with_margin(GENE_REGIONS, 100, min_gap=50)
    symbolic_pool = ["<DEL>", "<DUP>", "<INV>", "<INS>", "<CNV>"]
    for i, pos in enumerate(positions):
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        num_alts = random.randint(5, 8)
        alts = []
        # Always include at least one SNV, one indel, and one symbolic
        other_bases = [b for b in BASES if b != ref_base]
        alts.append(random.choice(other_bases))  # SNV
        ins_len = random.randint(1, 20)
        alts.append(ref_base + random_seq(ins_len))  # insertion
        alts.append(random.choice(symbolic_pool))  # symbolic
        # Fill remaining with mix
        while len(alts) < num_alts:
            choice = random.randint(0, 2)
            if choice == 0:
                # Another SNV
                b = random.choice(other_bases)
                if b not in alts:
                    alts.append(b)
                else:
                    alts.append(ref_base + random_seq(random.randint(1, 10)))
            elif choice == 1:
                # Indel
                alts.append(ref_base + random_seq(random.randint(1, 30)))
            else:
                # Symbolic
                sym = random.choice(symbolic_pool)
                if sym not in alts:
                    alts.append(sym)
                else:
                    alts.append(ref_base + random_seq(random.randint(2, 15)))
        alts = alts[:num_alts]
        # Build INFO -- if any symbolic SV, add END
        sv_len = random.randint(500, 50000)
        end = pos + sv_len
        has_symbolic = any(a.startswith("<") for a in alts)
        info = f"END={end};SVTYPE=DEL;SVLEN=-{sv_len}" if has_symbolic else "."
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_higharity_{i+1:04d}",
                ref=ref_base,
                alt=",".join(alts),
                info=info,
            )
        )

    # Category 41: Contig-boundary variants (50)
    # 25 at position 1, 25 near contig end
    for i in range(25):
        pos = 1
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        alt_base = mutate_base(ref_base)
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_boundary_start_{i+1:04d}",
                ref=ref_base,
                alt=alt_base,
            )
        )

    for i in range(25):
        margin = random.randint(1, 100)
        pos = contig_len - margin
        ref_base = ref.fetch_base(CONTIG, pos)
        if ref_base == "N":
            ref_base = "A"
        alt_base = mutate_base(ref_base)
        records.append(
            VcfRecord(
                chrom=CONTIG,
                pos=pos,
                id=f"synth_boundary_end_{i+1:04d}",
                ref=ref_base,
                alt=alt_base,
            )
        )

    header = build_vcf_header(
        ["END", "SVTYPE", "SVLEN", "LEN", "SVCLAIM", "PARID"],
        ["DEL", "DUP", "INS", "INV", "CNV"],
        ['##ALT=<ID=*,Description="Unspecified alternate allele">'],
        fileformat="VCFv4.5",
    )
    n = write_vcf(path, header, records)
    return (
        fname,
        n,
        [
            "Cat37:RefBlock",
            "Cat38:SVCLAIM",
            "Cat39:PARID",
            "Cat40:HighArity",
            "Cat41:ContigBoundary",
        ],
    )


# Main

GENERATORS = {
    "01_snv": gen_01_snv,
    "02_mnp_complex": gen_02_mnp_complex,
    "03_small_indels": gen_03_small_indels,
    "04_large_explicit": gen_04_large_explicit,
    "05_symbolic_del_ins": gen_05_symbolic_del_ins,
    "06_symbolic_dup_inv": gen_06_symbolic_dup_inv,
    "07_cnv_repeat": gen_07_cnv_repeat,
    "08_mobile_elements": gen_08_mobile_elements,
    "09_breakends": gen_09_breakends,
    "10_special_alleles": gen_10_special_alleles,
    "11_multi_allelic": gen_11_multi_allelic,
    "12_complex_imprecise": gen_12_complex_imprecise,
    "13_vcf45_features": gen_13_vcf45_features,
}


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic VCF files for SV variant classification testing.",
    )
    parser.add_argument(
        "--fasta",
        default=None,
        help=(
            "Path to an indexed reference FASTA (.fa + .fai) for the chosen "
            "assembly. Required: REF alleles are read from it, so there is no "
            "usable default. Set $VEP_FASTA_GRCH37 or $VEP_FASTA_GRCH38 to avoid "
            "passing it every time."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory for VCF files. Default: tests/sv_validation/grch37 or tests/sv_validation/grch38",
    )
    parser.add_argument(
        "--assembly",
        default="GRCh37",
        choices=["GRCh37", "GRCh38"],
        help="Reference genome assembly (default: GRCh37)",
    )
    parser.add_argument(
        "--file",
        default=None,
        help="Generate only this file (e.g., '01_snv'). Default: all files.",
    )
    args = parser.parse_args()

    # Set assembly-dependent defaults
    if args.fasta is None:
        env_var = f"VEP_FASTA_{args.assembly.upper()}"
        args.fasta = os.environ.get(env_var)
        if not args.fasta:
            parser.error(
                f"--fasta is required (or set ${env_var}): REF alleles are read "
                f"from the reference FASTA for {args.assembly}."
            )
    if args.output_dir is None:
        args.output_dir = f"tests/sv_validation/{args.assembly.lower()}"

    # Set gene regions and assembly for selected assembly
    global GENE_REGIONS, ASSEMBLY
    GENE_REGIONS = GENE_REGIONS_BY_ASSEMBLY[args.assembly]
    ASSEMBLY = args.assembly

    # Validate inputs
    if not os.path.isfile(args.fasta):
        print(f"ERROR: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.fasta + ".fai"):
        print(f"ERROR: FASTA index not found: {args.fasta}.fai", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    random.seed(RANDOM_SEED)

    print("Loading reference chromosome...")
    ref_fasta = RefFasta(args.fasta, contig=CONTIG)
    print()

    # Determine which generators to run
    if args.file:
        key = args.file.replace(".vcf", "")
        if key not in GENERATORS:
            print(
                f"ERROR: Unknown file '{key}'. Available: {', '.join(GENERATORS.keys())}",
                file=sys.stderr,
            )
            sys.exit(1)
        gens = {key: GENERATORS[key]}
    else:
        gens = GENERATORS

    total_variants = 0
    all_categories = []
    results = []

    print(f"Generating {len(gens)} VCF file(s) in {args.output_dir}\n")
    print(f"{'File':<35} {'Variants':>8}  Categories")
    print("-" * 80)

    for name, gen_func in gens.items():
        fname, count, categories = gen_func(ref_fasta, args.output_dir)
        total_variants += count
        all_categories.extend(categories)

        # bcftools validation
        vcf_path = os.path.join(args.output_dir, fname)
        valid = validate_vcf(vcf_path)
        status = "OK" if valid else "WARN"

        cat_str = ", ".join(categories)
        print(f"  {fname:<33} {count:>8}  {cat_str}  [{status}]")
        results.append((fname, count, categories, valid))

    print("-" * 80)
    print(f"  {'TOTAL':<33} {total_variants:>8}  {len(all_categories)} categories")
    print()

    # Summary
    failed = [r for r in results if not r[3]]
    if failed:
        print(f"WARNING: {len(failed)} file(s) had validation warnings:")
        for fname, _, _, _ in failed:
            print(f"  - {fname}")
    else:
        print("All files passed bcftools validation.")


if __name__ == "__main__":
    main()
