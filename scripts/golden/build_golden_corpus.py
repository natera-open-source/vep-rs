#!/usr/bin/env python3
"""Build a golden corpus: exemplar input records for every consequence-set combination a
VEP reference output contains, keyed by Ensembl release and assembly.

The corpus is the fixture behind `crates/vep-cli/tests/golden.rs` and
`crates/vep-cli/tests/format_parity.rs`: a small VCF whose records, annotated by VEP and by
vep-rs against the same transcript cache, must yield the same rows in every output format.

Subcommands:

  select    Stream one or more VEP default-format reference outputs, enumerate the distinct
            consequence sets they contain, choose up to K exemplar records per set
            (stratified by variant class and transcript strand, and by the family key of
            `--family-columns` when given; when `--cache` names the JSON cache the fixture
            is pruned from, the records whose transcripts cost the fixture cache the fewest
            bytes come first, then the records that exercise transcript flags, coordinate
            ranges and unmappable ends; every variant class and every family stratum seen
            keeps at least one exemplar), resolve each exemplar back to its input VCF
            record, and write `variants.vcf` plus `manifest.json`. The manifest's `covers`
            list names the parity-matrix rows the corpus exercises: `--covers <item>`
            (repeatable) writes it, each name a row of `--matrix` when one is given, a corpus
            regenerated over an existing manifest keeps that manifest's list, its `masks`
            (each a class `classify` knows, else the run is refused before anything is
            written), its `reference_warnings` (the texts Ensembl VEP wrote about the
            run rather than about one record) and its `flag_sets` (each set's run
            definition and `covers`; the divergence keys `classify --flag-set` wrote
            describe the old records and are dropped, as the manifest's own are), and a
            run with neither
            writes nothing. `--include <sid>=<name>` (repeatable) names a record of suite
            `sid`'s input, by its ID or by VEP's name for an ID-less record, that the
            corpus carries whatever the sampling chooses: its rows in the suite's
            reference join the chosen exemplars after the sampling, so a record chosen
            for what it co-locates with rather than for its consequence set stays in
            even when another record of the same stratum is sampled; the names are
            recorded as `included_records`. `--keep-records` keeps every record of the
            manifest being regenerated the same way, each resolved again from its
            suite's input and rowed again from its suite's reference, so a regeneration
            loses no record (the count is recorded as `kept_records`, and the manifest's
            `included_records` are carried into the new map as named records). A named or
            kept record is an exemplar of every combination its reference rows show, so
            `variant_classes` and `family_strata` count each of those rows, and it survives
            `--max-records`; one with no row in its suite's reference is an error, and
            nothing is written, so a corpus whose manifest carries a record the reference
            skipped (`reference_warning`) is regenerated without `--keep-records`, which
            keeps no record by name: that record returns only when the sampling picks it.
            `--keep-samples N` keeps `FORMAT` and the first N sample
            columns of every record and of the `#CHROM` header, and adds records until
            every genotype class the inputs show among those samples has
            `--genotype-exemplars` records (`--require-genotype-class` names a class the
            corpus must carry whether or not the inputs show it).
  classify  Compare a VEP default-format output for `variants.vcf` with a vep-rs output for
            the same file and record, per exemplar row VEP emits that vep-rs does not, the
            documented divergence class that explains it (or `unexplained_residual`), so the
            golden test asserts the documented difference instead of equality. A record no
            VEP row names is one VEP dropped before annotation: it must carry
            `reference_warning` (VEP's reason, from its warnings file, recorded when the
            corpus is built), it gets `reference_rows` 0, and every key only vep-rs emits
            for it is classed `reference_skipped_record`; any other key only vep-rs emits is
            `unexplained_residual`. A zero-row record without the warning, or a warned record
            with rows, is an error, unless the record carries no row by its nature: a
            non-variant record (ALT `.`), one the run's expected JSON names as a feature
            with no allele to annotate (`most_severe_consequence` `?`: every ALT a `*` or
            `<DEL:*>`, equal to the reference or a repeat) or, under a row filter of the run
            (`--coding_only`, `--no_intergenic`), one vep-rs writes no row for either. The
            records are the run's input's, numbered as the golden tests number them: the
            manifest's for `variants.vcf`, and for a flag set over another `input` that
            file's lines, a line the manifest records (a twin of `variants.vcf`) keeping
            its entry and a file in another format or delimiter, `variants.vcf` re-encoded
            record for record, taking the manifest's records in order. A run whose `masks` list
            `regulatory_feature_rows` has its VEP rows of `Feature_type` RegulatoryFeature
            and MotifFeature left out first, the rows the golden test sets aside for a run
            under `--regulatory`; one whose `masks` list `order_only` has its `DOMAINS`,
            `FREQS` and `ZYG` values compared as multisets of their items under `--field`,
            as the golden test compares them. `covers` and `flag_sets` are kept; with `--flag-set
            <name>` the divergence keys are written into that flag set instead of the
            manifest's top level, so each flag set documents its own run, under its own
            `masks` when its entry carries the key and the manifest's otherwise.
  check     Recompute the observed combination space from the reference outputs and fail
            when a combination has no exemplar in the manifest; with `--matrix` fail when a
            `covers` name (the manifest's or a flag set's) is not a matrix row; with
            `--family-columns` fail when a family stratum the references show has no
            exemplar among the corpus's records; on a corpus built with `--keep-samples`
            fail when `variants.vcf` does not carry `FORMAT` and the manifest's samples on
            its header and every record, or when a genotype class the inputs show among
            those samples (or one the manifest or `--require-genotype-class` requires) has
            no record in the corpus.
  covers    List (`--list`) or extend (`--add <item>`, repeatable) the `covers` of a corpus,
            or of one of its flag sets (`--flag-set <name>`); a name that is not a row of
            `--matrix` is refused and nothing is written.
  flag-set  Write one `flag_sets` entry: a run of the corpus under its own `flags` beside
            the manifest's default run, with its own expected outputs. An entry carries
            `name`, `flags` (the whole flag list of that run), `input` (a file in the corpus,
            `variants.vcf` by default), `stdin` (the input piped on standard input, with
            `-i STDIN`; its `flags` must name `--format`, which Ensembl VEP needs on
            STDIN), `unordered` (the reference's rows within a record come out in
            Ensembl VEP's hash order, so the golden tests compare the rows and not
            their order), `formats` (the output formats compared, all four by default),
            `expected_dir` (`expected_<name>` by default), `perl_mode` (written as
            `database` when the reference was generated against the Ensembl database
            rather than offline) and `covers`; `input` and `expected_dir` are paths inside
            the corpus directory. A set that exists has those keys rewritten from the
            options and keeps every other key of its entry, its hand-written `masks` (the
            classes of its run alone, each one of `regulatory_feature_rows` and
            `order_only`, as the manifest's) and `reference_warnings` (the texts Ensembl
            VEP wrote about its run alone) among them; its divergence keys describe the
            run before the rewrite and go, so `classify --flag-set` runs after it.

Row keys follow `compare_vep_outputs.py`: (Location, Allele, Feature, Feature_type,
consequence set), consequence terms sorted and comma-joined. A record is identified by
VEP's (Uploaded_variation, Location) pair, which every row of one input record shares.

The family key of `--family-columns` is one class per named Extra column plus the variant
class: a bare name (`HGVSc`) classes a row by the column's presence (`HGVSc+`, `HGVSc-`), a
name ending in `=` (`IMPACT=`) by its value (`IMPACT=HIGH`); `check --family-columns`
asserts every stratum the references show is covered, which is the per-family coverage
claim behind a flag family's corpus.

The genotype classes of `--keep-samples` are read from the kept samples' GT values as
Ensembl VEP reads a genotype for `--individual` (ensembl-io BaseVCF4.pm 941-951 and
Parser/VCF.pm 627-651 at release 116): per sample `het`, `hom_alt`, `hom_ref`, `missing`,
`partially_missing`, `haploid`, `star` and `multi_allelic`; per record `all_missing`,
`no_carrier` and `one_hom_ref_rest_missing` (`genotype_labels` and
`record_genotype_classes` define each; a structural record shows none, since Ensembl VEP
annotates it without genotypes). A record shows every class any kept sample shows, the
manifest lists each record's classes and the count of records per class, and `check`
asserts every class the inputs show is covered.

Selection is deterministic under `--seed`: reservoir sampling bounds the candidates kept per
combination and the same seed reproduces the same corpus from the same references.
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import re
import sys
import zlib
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Iterator

GENERATOR_VERSION = "2"
RESERVOIR = 64
FAMILY_RESERVOIR = 8
SEQUENCE_CLASSES = ("snv", "deletion", "insertion", "mnv_or_complex")
OUTPUT_FORMATS = ("default", "tab", "vcf", "json")
PERL_MODES = ("offline", "database")
DEFAULT_INPUT = "variants.vcf"
# The manifest `masks` class under which the reference's regulatory rows are left out, and
# the `Feature_type` values of those rows.
REGULATORY_FEATURE_ROWS = "regulatory_feature_rows"
REGULATORY_FEATURE_TYPES = frozenset({"RegulatoryFeature", "MotifFeature"})
# The `masks` class under which the set-valued columns agree as a multiset of their items
# in any order (Ensembl VEP prints DOMAINS in the cache's protein_features order,
# BaseTranscriptVariation.pm:634-643, FREQS in hash order, OutputFactory.pm:1152-1153, and
# the ZYG list of --individual_zyg in hash order, OutputFactory.pm:983), and those columns;
# the default format joins the items with `,` (Utils.pm:185-198).
ORDER_ONLY = "order_only"
ORDER_ONLY_COLUMNS = ("DOMAINS", "FREQS", "ZYG")
# Every class a `masks` list, the manifest's or a flag set's, may name.
MASK_CLASSES = (REGULATORY_FEATURE_ROWS, ORDER_ONLY)
CLASSIFY_KEYS = ("divergences", "vep_rs_only_tuples", "divergence_summary", "vep_rs_only_summary",
                 "fields_compared", "field_divergences", "field_divergence_summary")
DEFAULT_COLUMNS = (
    "Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
    "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "Extra",
)
FIXED_COLUMNS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")


def opener(path: Path):
    """`gzip.open` for a `.gz` or `.bgz` (BGZF is gzip), else `open`."""
    return gzip.open if str(path).endswith((".gz", ".bgz")) else open


def normalize_consequence_set(raw: str) -> str:
    return ",".join(sorted({t for t in raw.strip().split(",") if t}))


def parse_extra(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if not raw or raw == "-":
        return out
    for part in raw.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        elif part:
            out[part] = ""
    return out


def variant_class(location: str, allele: str) -> str:
    """Infer VEP's variant shape from the Location and Allele columns of one row: one of
    `SEQUENCE_CLASSES` for a sequence allele; `breakend`, or `breakend_insertion` when the
    bracket allele carries inserted sequence (more than the anchor base outside the brackets);
    `single_breakend`; and `symbolic_<word>` for the class word VEP prints for a symbolic
    allele (`symbolic_deletion`, `symbolic_Alu_insertion`). Each is annotated by a path of its
    own, so each is a stratum of its own."""
    if "[" in allele or "]" in allele:
        bases = re.sub(r"[\[\]][^\[\]]*[\[\]]", "", allele)
        return "breakend_insertion" if len(bases) > 1 else "breakend"
    if allele.startswith(".") or allele.endswith("."):
        return "single_breakend"
    if allele and not re.fullmatch(r"[ACGTNacgtn*-]+", allele):
        return f"symbolic_{allele}"
    span = location.split(":", 1)[1] if ":" in location else location
    if allele == "-":
        return "deletion"
    if "-" in span:
        a, b = span.split("-", 1)
        if a.isdigit() and b.isdigit() and int(b) == int(a) + 1:
            return "insertion"
        return "mnv_or_complex"
    return "snv" if len(allele) == 1 else "mnv_or_complex"


def parse_family_columns(arg: str | None) -> list[str]:
    """`--family-columns HGVSc,IMPACT=`: a column name (presence) or a name ending in `=` (value)."""
    specs = [s.strip() for s in (arg or "").split(",") if s.strip()]
    for spec in specs:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+=?", spec):
            raise SystemExit(f"ERROR: [build_golden_corpus] --family-columns wants NAME or NAME=, got {spec!r}")
    return specs


def family_classes(extra: dict[str, str], specs: list[str]) -> tuple[str, ...]:
    """One class per spec: `NAME+` or `NAME-` for presence, `NAME=<value>` for a value (`-` when absent)."""
    classes = []
    for spec in specs:
        if spec.endswith("="):
            name = spec[:-1]
            classes.append(f"{name}={extra.get(name) or '-'}")
        else:
            present = spec in extra and extra[spec] not in ("", "-")
            classes.append(f"{spec}{'+' if present else '-'}")
    return tuple(classes)


SAMPLE_CLASSES = ("het", "hom_alt", "hom_ref", "missing", "partially_missing", "haploid", "star", "multi_allelic")
RECORD_CLASSES = ("all_missing", "no_carrier", "one_hom_ref_rest_missing")
GENOTYPE_CLASSES = SAMPLE_CLASSES + RECORD_CLASSES
STRUCTURAL_ALLELE = re.compile(r"[<\[\]][^*]+[>\]\[]|^\.\w+|\w+\.$")


def genotype_labels(gt: str, alts: list[str], record: str) -> set[str]:
    """The classes one sample's GT shows at a record, read as Ensembl VEP reads a genotype
    for `--individual` (ensembl-io BaseVCF4.pm 941-951 and Parser/VCF.pm 627-651, release
    116): the bits split on `|` when the GT has one, else on `/`; `.` bits dropped before
    anything else is read; an allele containing `*` outside the non-reference set. `missing`
    (no called bit) stands alone; otherwise `partially_missing` when a bit was dropped,
    `haploid` for a one-bit genotype, `star` for a `*` allele, and one of `multi_allelic`
    (two or more distinct non-reference alleles), `hom_ref` (none), `hom_alt` (one, every
    called bit the same) or `het`. A bit that is not an index of the record's alleles is
    refused, naming `record`."""
    bits = gt.split("|") if "|" in gt else gt.split("/")
    called_bits = [b for b in bits if b != "."]
    if not called_bits:
        return {"missing"}
    for b in called_bits:
        # str.isdigit alone accepts a superscript such as "²", which int() rejects.
        if not (b.isascii() and b.isdigit()) or int(b) > len(alts):
            raise failure(f"record {record} has a GT {gt!r} outside its alleles")
    labels: set[str] = set()
    if len(called_bits) < len(bits):
        labels.add("partially_missing")
    if len(bits) == 1:
        labels.add("haploid")
    called = {int(b) for b in called_bits}
    star = {i + 1 for i, a in enumerate(alts) if "*" in a}
    if star & called:
        labels.add("star")
    non_ref = called - {0} - star
    if len(non_ref) > 1:
        labels.add("multi_allelic")
    elif not non_ref:
        labels.add("hom_ref")
    elif len(called) == 1:
        labels.add("hom_alt")
    else:
        labels.add("het")
    return labels


def is_structural_record(ref: str, alts: list[str], info: str) -> bool:
    """Whether Ensembl VEP reads a record as a structural variant, which it annotates before
    the per-individual step and so without any individual's genotype (Parser/VCF.pm 236-244):
    its alleles are not all plain bases and it carries an INFO `SVTYPE` whose value Perl reads
    as true, or a symbolic, breakend or single-breakend allele (`<*>` is none of these). The
    value is read as ensembl-io's `get_info` stores it (BaseVCF4.pm 576-589): the piece after
    the first `=`, so a bare `SVTYPE` flag and `SVTYPE=` hold none, `SVTYPE=0` holds Perl's
    false `0`, and of a repeated key the last one counts. The allele test runs on the
    comma-joined ALT string as Ensembl VEP's does, so a single-breakend allele counts only at
    the string's start or end."""
    if re.fullmatch(r"[ACGT]+", "".join([ref, *alts])):
        return False
    svtype = ""
    for part in info.split(";"):
        pieces = part.split("=")
        if pieces[0] == "SVTYPE":
            svtype = pieces[1] if len(pieces) > 1 else ""
    if svtype not in ("", "0"):
        return True
    return STRUCTURAL_ALLELE.search(",".join(alts)) is not None


def record_genotype_classes(fields: list[str], keep: int) -> set[str]:
    """The genotype classes one record shows among its first `keep` samples: every kept
    sample's `genotype_labels`, plus `all_missing` (no kept sample has a called bit),
    `no_carrier` (a called bit somewhere but no non-reference allele in any kept sample)
    and `one_hom_ref_rest_missing` (one kept sample `hom_ref`, every other one `missing`).
    A structural record shows no class (`is_structural_record`)."""
    record = " ".join(fields[:5])
    if len(fields) < 9 + keep:
        raise failure(f"record {record} has {max(len(fields) - 9, 0)} sample columns, {keep} are kept")
    alts = fields[4].split(",")
    if is_structural_record(fields[3], alts, fields[7]):
        return set()
    keys = fields[8].split(":")
    if "GT" not in keys:
        raise failure(f"record {record} has no GT in its FORMAT")
    gt_index = keys.index("GT")
    classes: set[str] = set()
    missing = hom_ref = 0
    for sample in fields[9:9 + keep]:
        values = sample.split(":")
        labels = genotype_labels(values[gt_index] if gt_index < len(values) else ".", alts, record)
        classes |= labels
        missing += "missing" in labels
        hom_ref += "hom_ref" in labels
    if missing == keep:
        classes.add("all_missing")
    elif missing + hom_ref == keep:
        classes.add("no_carrier")
        if hom_ref == 1:
            classes.add("one_hom_ref_rest_missing")
    return classes


def genotype_class_names(names: list[str] | None) -> list[str]:
    """The `--require-genotype-class` names of a command line, each a genotype class, in order
    and without repeats."""
    names = list(dict.fromkeys(names or []))
    for n in names:
        if n not in GENOTYPE_CLASSES:
            raise failure(f"--require-genotype-class wants one of {', '.join(GENOTYPE_CLASSES)}, got {n!r}")
    return names


def vcf_lines(path: Path, columns: int | None = None) -> Iterator[tuple[bool, list[str]]]:
    """Every `#CHROM` header line (`True`, every column) and every record of at least eight
    fields (`False`) of a VCF, split on tabs; `##` lines are skipped. With `columns` a record
    is split into that many fields and one more holding the rest of the line, so a source of
    thousands of samples is not split beyond the columns read. A file that cannot be opened
    or read (a `.gz` that is not gzip fails on its first line, not at the open) is refused
    with its path."""
    try:
        with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    yield True, line.rstrip("\n").split("\t")
                    continue
                fields = line.rstrip("\n").split("\t", -1 if columns is None else columns)
                if len(fields) >= 8:
                    yield False, fields
    except (OSError, EOFError, zlib.error) as exc:
        raise failure(f"{path} is not a readable text file ({exc})") from None


def kept_samples(header: list[str], keep: int, path: Path, samples: list[str], source: str = "an earlier input") -> list[str]:
    """The first `keep` sample names of a `#CHROM` header line, refused when the line lacks
    `FORMAT` or that many samples, or when `samples` (the names `source` keeps) names others."""
    if len(header) < 9 or header[8] != "FORMAT":
        raise failure(f"{path} has no FORMAT column, so no sample columns")
    if len(header) < 9 + keep:
        raise failure(f"{path} has {len(header) - 9} sample columns, {keep} are kept")
    kept = header[9:9 + keep]
    if samples and kept != samples:
        raise failure(f"{path} keeps samples {','.join(kept)} where {source} keeps {','.join(samples)}: "
                      "every input must keep the same samples")
    return kept


@dataclass
class Row:
    suite: str
    uploaded_variation: str
    location: str
    allele: str
    feature: str
    feature_type: str
    consequence_set: str
    strand: str
    flags: str
    positions: tuple[str, str, str]
    family: tuple[str, ...] = ()

    @property
    def record_key(self) -> tuple[str, str]:
        return (self.uploaded_variation, self.location)

    @property
    def klass(self) -> str:
        return variant_class(self.location, self.allele)

    @property
    def family_stratum(self) -> str:
        """The family classes and the variant class, `;`-joined: the key `--family-columns` covers."""
        return ";".join((*self.family, self.klass))

    def preference(self) -> int:
        score = 0
        if self.flags:
            score += 4
        if any("?" in p for p in self.positions):
            score += 3
        if any("-" in p and p != "-" for p in self.positions):
            score += 2
        return score

    @property
    def span(self) -> int:
        """Genomic width of the Location, in bases (0 for a single position)."""
        coords = self.location.split(":", 1)[-1]
        if "-" not in coords:
            return 0
        a, b = coords.split("-", 1)
        try:
            return abs(int(b) - int(a))
        except ValueError:
            return 0


def iter_rows(path: Path, suite: str, family_specs: list[str] | None = None) -> Iterator[Row]:
    with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            while len(f) < 14:
                f.append("-")
            extra = parse_extra(f[13])
            yield Row(
                suite=suite,
                uploaded_variation=f[0],
                location=f[1],
                allele=f[2],
                feature=f[4],
                feature_type=f[5],
                consequence_set=normalize_consequence_set(f[6]),
                strand=extra.get("STRAND", ""),
                flags=extra.get("FLAGS", ""),
                positions=(f[7], f[8], f[9]),
                family=family_classes(extra, family_specs) if family_specs else (),
            )


def parse_suite_arg(arg: str) -> tuple[str, Path, Path]:
    """`sid=reference.txt=input.vcf[.gz]`."""
    parts = arg.split("=")
    if len(parts) != 3:
        raise SystemExit(f"ERROR: [build_golden_corpus] --suite wants sid=reference=input, got {arg!r}")
    return parts[0], Path(parts[1]), Path(parts[2])


BND_MATE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
SHARD_NAME = re.compile(r"^(\d+)-(\d+)\.json(\.gz)?$")


class CacheCost:
    """Compressed bytes a row adds to the fixture cache `prune_json_cache.py` cuts for it:
    every transcript of the source cache whose span, extended by `flank`, overlaps the row's
    location (and a bracket breakend's mate position), costed at its own deflate size. A
    shard is read on first touch; the source stores a transcript in every shard it reaches,
    so only the shards overlapping the flanked location are read."""

    def __init__(self, cache: Path, flank: int) -> None:
        self.flank = flank
        self.shards: dict[str, list[tuple[int, int, Path]]] = {}
        for chrom_dir in (cache / "transcripts").iterdir():
            self.shards[chrom_dir.name] = sorted(
                (int(m.group(1)), int(m.group(2)), shard)
                for shard in chrom_dir.iterdir()
                if (m := SHARD_NAME.match(shard.name))
            )
        self.loaded: dict[Path, list[tuple[int, int, int]]] = {}
        self.memo: dict[tuple[str, int, int], int] = {}

    def transcripts(self, shard: Path) -> list[tuple[int, int, int]]:
        if shard not in self.loaded:
            with opener(shard)(shard, "rt", encoding="utf-8") as handle:
                records = json.load(handle)
            self.loaded[shard] = [
                (int(t["start"]), int(t["end"]),
                 len(zlib.compress(json.dumps(t, separators=(",", ":"), ensure_ascii=False).encode("utf-8"))))
                for t in records
            ]
        return self.loaded[shard]

    def span_cost(self, chrom: str, start: int, end: int) -> int:
        key = (chrom, start, end)
        if key not in self.memo:
            lo, hi = start - self.flank, end + self.flank
            total = 0
            for shard_start, shard_end, shard in self.shards.get(chrom, ()):
                if shard_start <= hi and shard_end >= lo:
                    total += sum(b for t_start, t_end, b in self.transcripts(shard) if t_start <= hi and t_end >= lo)
            self.memo[key] = total
        return self.memo[key]

    def row_cost(self, row: Row) -> int:
        chrom, _, coords = row.location.partition(":")
        a, _, b = coords.partition("-")
        try:
            start, end = int(a), int(b or a)
        except ValueError:
            return 0
        start, end = min(start, end), max(start, end)
        total = self.span_cost(chrom, start, end)
        for m in BND_MATE.finditer(row.allele):
            mate_chrom, mate_pos = m.group(1).removeprefix("chr"), int(m.group(2))
            if mate_chrom == chrom:
                # A same-chromosome mate is annotated over the whole breakend-pair span.
                total = self.span_cost(chrom, min(start, mate_pos), max(end, mate_pos))
            else:
                total += self.span_cost(mate_chrom, mate_pos, mate_pos)
        return total


# ----------------------------------------------------------------------------- select
def reservoir_step(rng: random.Random, bucket: list, item, offered: int, capacity: int) -> None:
    """One reservoir step: `bucket` stays a uniform sample of at most `capacity` of the
    `offered` items seen so far for its key."""
    if len(bucket) < capacity:
        bucket.append(item)
    else:
        j = rng.randrange(offered)
        if j < capacity:
            bucket[j] = item


def reservoir_pass(
    suites: list[tuple[str, Path, Path]],
    seed: int,
    family_specs: list[str] | None = None,
    named: dict[str, set[str]] | None = None,
) -> tuple[dict[str, list[Row]], Counter[str], dict[str, Counter[str]], dict[str, list[Row]], dict[tuple[str, str], list[Row]]]:
    """One streaming pass: per consequence set, a bounded uniform sample of rows plus counts,
    per family stratum (when `family_specs` is given) a smaller sample of its own, so a
    stratum too rare to reach any consequence set's sample still has candidates, and every
    row of the records `named` (per suite, their Uploaded_variation names), keyed by suite
    and name, which the samples' bounds never touch."""
    rng = random.Random(seed)
    reservoir: dict[str, list[Row]] = defaultdict(list)
    seen: Counter[str] = Counter()
    per_suite: dict[str, Counter[str]] = defaultdict(Counter)
    families: dict[str, list[Row]] = defaultdict(list)
    family_seen: Counter[str] = Counter()
    named_rows: dict[tuple[str, str], list[Row]] = defaultdict(list)
    for sid, ref, _inp in suites:
        wanted_names = (named or {}).get(sid, set())
        for row in iter_rows(ref, sid, family_specs):
            c = row.consequence_set
            seen[c] += 1
            per_suite[sid][c] += 1
            reservoir_step(rng, reservoir[c], row, seen[c], RESERVOIR)
            if family_specs:
                stratum = row.family_stratum
                family_seen[stratum] += 1
                reservoir_step(rng, families[stratum], row, family_seen[stratum], FAMILY_RESERVOIR)
            if row.uploaded_variation in wanted_names:
                named_rows[(sid, row.uploaded_variation)].append(row)
    return reservoir, seen, per_suite, families, named_rows


def exemplar_order(row: Row, max_span: int, cost: Callable[[Row], int] | None) -> tuple[int, int, int]:
    """Sort key of a candidate row: no wider than `max_span` (when > 0) first, then the fewest
    cache bytes under `cost`, then the highest preference."""
    return (1 if max_span and row.span > max_span else 0, cost(row) if cost else 0, -row.preference())


def ensure_class_exemplars(
    chosen: dict[str, list[Row]],
    reservoir: dict[str, list[Row]],
    max_span: int,
    cost: Callable[[Row], int] | None,
) -> dict[str, list[Row]]:
    """Every variant class the references show keeps an exemplar: a class no chosen row has
    takes, from the whole reservoir, the row that sorts first in `exemplar_order`, added to the
    combination it was sampled for."""
    covered = {r.klass for rows in chosen.values() for r in rows}
    best: dict[str, tuple[tuple[int, int, int], str, Row]] = {}
    for combo, rows in reservoir.items():
        for r in rows:
            if r.klass in covered:
                continue
            key = exemplar_order(r, max_span, cost)
            if r.klass not in best or key < best[r.klass][0]:
                best[r.klass] = (key, combo, r)
    for _key, combo, r in best.values():
        chosen.setdefault(combo, []).append(r)
    return chosen


def ensure_family_exemplars(
    chosen: dict[str, list[Row]],
    families: dict[str, list[Row]],
    max_span: int,
    cost: Callable[[Row], int] | None,
) -> dict[str, list[Row]]:
    """Every family stratum the references show keeps an exemplar: a stratum no chosen row
    has takes, from its own sample, the row that sorts first in `exemplar_order`, added to
    the combination it was sampled for."""
    covered = {r.family_stratum for rows in chosen.values() for r in rows}
    for stratum, rows in sorted(families.items()):
        if stratum in covered:
            continue
        r = min(rows, key=lambda row: exemplar_order(row, max_span, cost))
        chosen.setdefault(r.consequence_set, []).append(r)
    return chosen


def choose_exemplars(
    reservoir: dict[str, list[Row]],
    k: int,
    rng: random.Random,
    max_span: int = 0,
    focus_terms: Iterable[str] = (),
    focus_k: int | None = None,
    cost: Callable[[Row], int] | None = None,
) -> dict[str, list[Row]]:
    """Up to K rows per combination, distinct (class, strand, family classes) first, in
    `exemplar_order`.

    Rows wider than `max_span` bases (when > 0) come after every narrower row, so a
    combination takes a wide structural variant only when nothing narrower shows it:
    each wide record pulls every transcript of its span into the fixture cache. With
    `cost` (bytes the row adds to the fixture cache, `CacheCost.row_cost`) the cheapest
    rows come next, then the higher preference; without it the preference alone orders
    the rows and ties fall in random order.

    A combination containing any of `focus_terms` takes up to `focus_k` rows instead of
    `k`, so a corpus can carry many exemplars of the consequence shapes it exists to test
    beside one exemplar of every other shape.
    """
    focus = set(focus_terms)
    chosen: dict[str, list[Row]] = {}
    for combo in sorted(reservoir):
        want = focus_k if focus_k is not None and focus & set(combo.split(",")) else k
        cands = list(reservoir[combo])
        rng.shuffle(cands)
        cands.sort(key=lambda r: exemplar_order(r, max_span, cost))
        picked: list[Row] = []
        strata_seen: set[tuple[str, str, tuple[str, ...]]] = set()
        for r in cands:
            stratum = (r.klass, r.strand, r.family)
            if stratum in strata_seen:
                continue
            picked.append(r)
            strata_seen.add(stratum)
            if len(picked) == want:
                break
        for r in cands:
            if len(picked) == want:
                break
            if r not in picked:
                picked.append(r)
        chosen[combo] = picked
    return chosen


def vep_auto_names(chrom: str, pos: int, ref: str, alts: list[str]) -> list[str]:
    """Reproduce VEP's auto-generated Uploaded_variation for an ID-less VCF record, in the
    two forms releases 116 and 115 write.

    Release 116 (ensembl-vep PR 1941) names the record from the raw line,
    `CHR_POS_REF/ALT1/ALT2`, and that name comes first. Release 115's names
    follow: the alleles are the RAW REF and ALTs joined by `/` (VEP's
    `nontrimmed_allele_string`) and the start is the variant's after parsing: unchanged
    when no allele differs in length from REF; for a bi-allelic indel the first shared
    base is chopped and the start then advances over the shared leading bases
    (`trim_sequences`, whose trailing trim moves only the end); for a multi-allelic
    indel the first base is chopped when every allele (ignoring `*`) shares it. Symbolic
    and breakend alleles: `chr_start_<ALT...>` with start POS or POS+1 (both are
    returned; which one VEP used depends on whether the REF base was padding).
    """
    alleles = [ref] + alts
    raw = "/".join(alleles)
    names = [f"{chrom}_{pos}_{raw}"]
    if any(a.startswith("<") or "[" in a or "]" in a for a in alts):
        joined = "/".join(alts)
        names += [f"{chrom}_{pos}_{joined}", f"{chrom}_{pos + 1}_{joined}"]
    elif len(alts) == 1:
        r, a, start = ref, alts[0], pos
        if len(r) != len(a):
            if r[:1] == a[:1]:
                r, a, start = r[1:] or "-", a[1:] or "-", start + 1
            while r and a and r[0] == a[0]:
                r, a, start = r[1:], a[1:], start + 1
            names.append(f"{chrom}_{start}_{raw}")
    elif any(len(a) != len(ref) for a in alts):
        firsts = {x[:1] for x in alleles if "*" not in x}
        if len(firsts) == 1:
            names.append(f"{chrom}_{pos + 1}_{raw}")
    return list(dict.fromkeys(names))


MULTI_ALLELIC_STRATA = ("multi_allelic_snv", "multi_allelic_indel", "multi_allelic_indel_trimmable")


def multi_allelic_stratum(ref: str, alts: list[str]) -> str:
    """Which multi-allelic stratum a record belongs to: SNV-only alleles; alleles of
    differing length; or such alleles that VEP's bi-allelic minimisation WOULD trim
    further after the shared first base is chopped (a shared prefix or suffix left
    between REF and some ALT), the shape on which a per-allele minimiser and VEP's
    untouched multi-allelic record disagree on positions and distances."""
    if not any(len(a) != len(ref) for a in alts):
        return "multi_allelic_snv"
    alleles = [ref] + alts
    if len({x[:1] for x in alleles}) == 1:
        alleles = [x[1:] or "-" for x in alleles]
    r = alleles[0]
    for a in alleles[1:]:
        if r == "-" or a == "-" or len(r) == len(a):
            continue
        if r[0] == a[0] or r[-1] == a[-1]:
            return "multi_allelic_indel_trimmable"
    return "multi_allelic_indel"


@dataclass
class ResolvedInputs:
    """What one scan of the input VCFs yields: per suite the record lines whose
    Uploaded_variation (ID, or VEP's auto-name) was wanted; the `#CHROM` header line the
    corpus writes (eight columns, or `FORMAT` and the kept samples after them); per suite
    and multi-allelic stratum a uniform sample of records; and, with kept samples, the names
    kept, how many records show each genotype class and a uniform sample of records per
    class."""
    found: dict[str, dict[tuple[str, str], list[str]]]
    header: str
    multi: dict[str, dict[str, list[str]]]
    samples: list[str]
    genotype_seen: Counter[str]
    genotype_pool: dict[str, list[tuple[str, str]]]


def resolve_records(
    inputs: list[tuple[str, Path]],
    wanted: dict[str, set[str]],
    multi_per_suite: int = 0,
    seed: int = 1,
    keep_samples: int = 0,
    genotype_per_class: int = 0,
) -> ResolvedInputs:
    """Scan each input VCF once. A record line is its first eight columns, or with
    `keep_samples` its first eight, `FORMAT` and the first `keep_samples` sample columns.
    The multi-allelic sample (`multi_per_suite` sequence records per suite and stratum:
    SNV-only alleles; alleles of differing length) and the genotype-class sample
    (`genotype_per_class` records per class the kept samples show) are records no
    consequence set can select for: the first exercises VEP's treatment of a multi-allelic
    record as one variant, the second the genotype shapes of `--individual`."""
    rng = random.Random(seed)
    found: dict[str, dict[tuple[str, str], list[str]]] = defaultdict(lambda: defaultdict(list))
    multi: dict[str, dict[str, list[str]]] = defaultdict(lambda: {k: [] for k in MULTI_ALLELIC_STRATA})
    multi_seen: Counter[tuple[str, str]] = Counter()
    header: str | None = None
    samples: list[str] = []
    genotype_seen: Counter[str] = Counter()
    genotype_pool: dict[str, list[tuple[str, str]]] = defaultdict(list)
    kept_columns = 9 + keep_samples if keep_samples else 8
    for sid, path in inputs:
        want = wanted.get(sid, set())
        if not want and not multi_per_suite and not keep_samples:
            continue
        named = False
        for is_header, f in vcf_lines(path, kept_columns):
            if is_header:
                if keep_samples:
                    samples = kept_samples(f, keep_samples, path, samples)
                    named = True
                if header is None:
                    header = "\t".join(f[:kept_columns])
                continue
            if keep_samples and not named:
                raise failure(f"{path} has no #CHROM header line naming its samples")
            chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
            alts = alt.split(",")
            record = "\t".join(f[:kept_columns])
            if keep_samples:
                # Sorted so the draws fall on the classes in one order in every interpreter.
                for c in sorted(record_genotype_classes(f, keep_samples)):
                    genotype_seen[c] += 1
                    if genotype_per_class:
                        reservoir_step(rng, genotype_pool[c], (sid, record), genotype_seen[c], genotype_per_class)
            if multi_per_suite and len(alts) > 1 and not any(a.startswith("<") or "[" in a or "]" in a or a == "*" for a in alts):
                stratum = multi_allelic_stratum(ref, alts)
                multi_seen[(sid, stratum)] += 1
                reservoir_step(rng, multi[sid][stratum], record, multi_seen[(sid, stratum)], multi_per_suite)
            if not want:
                continue
            names: list[str] = []
            if vid != ".":
                names.append(vid.split(";")[0])
            names.extend(vep_auto_names(chrom, pos, ref, alts))
            for n in names:
                if n in want:
                    found[sid][(n, chrom)].append(record)
                    break
    if keep_samples and not samples:
        raise failure("no input has a #CHROM header line naming its samples")
    return ResolvedInputs(found, header or "\t".join(FIXED_COLUMNS), multi, samples, genotype_seen, genotype_pool)


def run_options(args: argparse.Namespace) -> dict:
    """The manifest keys the golden harness passes to both engines: `flags` (extra
    command-line flags such as `--hgvs`) and `fasta` (a gzipped reference beside the
    manifest, decompressed and passed to `--fasta`)."""
    options: dict = {}
    if getattr(args, "flag", None):
        options["flags"] = list(args.flag)
    if getattr(args, "fasta_name", None):
        options["fasta"] = args.fasta_name
    return options


def failure(message: str) -> SystemExit:
    return SystemExit(f"ERROR: [build_golden_corpus] {message}")


def is_name_list(value: object) -> bool:
    return isinstance(value, list) and all(isinstance(n, str) and n.strip() for n in value)


def masks_shape_error(masks: object, owner: str) -> str | None:
    """The reason a `masks` list (`owner` names it: the manifest's own or a flag set's)
    cannot be applied: not a list of names, or a class other than `MASK_CLASSES`."""
    if not is_name_list(masks):
        return f"{owner} is not a list of names"
    unknown = sorted(set(masks) - set(MASK_CLASSES))
    if unknown:
        return f"{owner} names unknown classes {unknown} (the classes are {' and '.join(MASK_CLASSES)})"
    return None


FLAG_SET_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*")


def is_flag_set_name(value: object) -> bool:
    return isinstance(value, str) and value != "default" and FLAG_SET_NAME.fullmatch(value) is not None


def is_corpus_path(value: object) -> bool:
    """A non-empty path relative to the corpus directory with no `..` segment."""
    if not isinstance(value, str) or not value:
        return False
    path = Path(value)
    return not path.is_absolute() and ".." not in path.parts


def names_format(flags: list[str]) -> bool:
    """Whether a flag list carries `--format` (as its own token or as `--format=<format>`)."""
    return any(f == "--format" or f.startswith("--format=") for f in flags)


def flag_set_shape_error(fs: object) -> str | None:
    """The reason one `flag_sets` entry cannot be run: the shapes `flag-set` writes, so a
    hand-written entry is held to the same rules; its hand-written `masks` (classes of
    `MASK_CLASSES`) and `reference_warnings` (non-empty texts) are held to the rules the
    manifest's own are, as the golden harness holds them."""
    name = fs.get("name") if isinstance(fs, dict) else None
    if not is_flag_set_name(name):
        return "has a flag_sets entry without a name (letters, digits, _ and -, not default)"
    if "flags" in fs and not (isinstance(fs["flags"], list) and all(isinstance(f, str) for f in fs["flags"])):
        return f"flag set {name} flags is not a list of strings"
    for key in ("input", "expected_dir"):
        if key in fs and not is_corpus_path(fs[key]):
            return f"flag set {name} {key} is not a path inside the corpus directory"
    for key in ("stdin", "unordered"):
        if key in fs and not isinstance(fs[key], bool):
            return f"flag set {name} {key} is not true or false"
    if "formats" in fs and not (isinstance(fs["formats"], list) and fs["formats"]
                                and all(f in OUTPUT_FORMATS for f in fs["formats"])):
        return f"flag set {name} formats is not a non-empty subset of {','.join(OUTPUT_FORMATS)}"
    if fs.get("stdin") and not names_format(fs.get("flags") or []):
        return f"flag set {name} runs on stdin without --format among its flags"
    if "perl_mode" in fs and fs["perl_mode"] not in PERL_MODES:
        return f"flag set {name} perl_mode is not one of {', '.join(PERL_MODES)}"
    for key, kind in (("covers", "names"), ("reference_warnings", "texts")):
        if key in fs and not is_name_list(fs[key]):
            return f"flag set {name} {key} is not a list of {kind}"
    if "masks" in fs:
        return masks_shape_error(fs["masks"], f"flag set {name} masks")
    return None


def manifest_shape_error(manifest: dict) -> str | None:
    """The reason a manifest's `covers`, `masks`, `reference_warnings`, `flag_sets`, `cache_layout`,
    `fasta_in_cache`, `included_records` or sample keys cannot be read: `covers` and
    `reference_warnings` must be lists of non-empty strings, `masks` a list of `MASK_CLASSES`
    (`masks_shape_error`), `flag_sets` a
    list of entries of the shape `flag_set_shape_error` holds them to, with distinct names,
    `cache_layout`, when present, `ensembl` (the cache at `cache_root/homo_sapiens/<release>_<assembly>`,
    so `release` and `assembly` must be named) rather than the default `json_cache/`, `fasta_in_cache`,
    when present, a boolean that is true only beside a `fasta` (the golden run then finds the reference
    inside its cache directory instead of receiving `--fasta`), `keep_samples` a count,
    `samples` the names of that many kept samples when the count is above zero, and
    `included_records`, when present, a map of non-empty suite names to non-empty lists of
    record names."""
    if "covers" in manifest and not is_name_list(manifest["covers"]):
        return "covers is not a list of names"
    if "masks" in manifest:
        reason = masks_shape_error(manifest["masks"], "masks")
        if reason:
            return reason
    if "reference_warnings" in manifest and not is_name_list(manifest["reference_warnings"]):
        return "reference_warnings is not a list of texts"
    keep = manifest.get("keep_samples", 0)
    if not isinstance(keep, int) or isinstance(keep, bool) or keep < 0:
        return "keep_samples is not a count"
    if keep and not (is_name_list(manifest.get("samples")) and len(manifest["samples"]) == keep):
        return f"samples does not name {keep} kept samples"
    required = manifest.get("required_genotype_classes", [])
    if not (isinstance(required, list) and all(c in GENOTYPE_CLASSES for c in required)):
        return "required_genotype_classes is not a list of genotype classes"
    flag_sets = manifest.get("flag_sets", [])
    if not isinstance(flag_sets, list):
        return "flag_sets is not a list"
    names: set[str] = set()
    for fs in flag_sets:
        reason = flag_set_shape_error(fs)
        if reason:
            return reason
        if fs["name"] in names:
            return f"names flag set {fs['name']} twice"
        names.add(fs["name"])
    if "cache_layout" in manifest and manifest["cache_layout"] != "ensembl":
        return f'cache_layout is {json.dumps(manifest["cache_layout"])}, not "ensembl"'
    if manifest.get("cache_layout") == "ensembl" and not all(isinstance(manifest.get(k), str) and manifest[k] for k in ("release", "assembly")):
        return "cache_layout ensembl wants release and assembly"
    included = manifest.get("included_records", {})
    if not (isinstance(included, dict)
            and all(isinstance(sid, str) and sid.strip() and entry and is_name_list(entry) for sid, entry in included.items())):
        return "included_records is not a map of suite to record names"
    if "fasta_in_cache" in manifest and not isinstance(manifest["fasta_in_cache"], bool):
        return f'fasta_in_cache is {json.dumps(manifest["fasta_in_cache"])}, not a boolean'
    if manifest.get("fasta_in_cache") and not (isinstance(manifest.get("fasta"), str) and manifest["fasta"]):
        return "fasta_in_cache wants fasta"
    return None


def read_manifest(corpus: Path) -> dict:
    """The corpus manifest, refused unless it is a JSON object whose `covers`, `flag_sets` and
    `cache_layout` have the shapes `manifest_shape_error` describes."""
    path = corpus / "manifest.json"
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise failure(f"{path} is not readable JSON ({exc})") from None
    if not isinstance(manifest, dict):
        raise failure(f"{path} is not a JSON object")
    reason = manifest_shape_error(manifest)
    if reason:
        raise failure(f"{path} {reason}")
    return manifest


def write_manifest(corpus: Path, manifest: dict) -> None:
    (corpus / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")


def matrix_items(path: Path) -> set[str]:
    """The `item` column of a parity matrix (tab-separated, `item` first)."""
    if not path.is_file():
        raise failure(f"{path} is not a file: name the parity matrix with --matrix")
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        raise failure(f"{path} is not a readable text file ({exc})") from None
    if not lines or lines[0].split("\t", 1)[0] != "item":
        raise failure(f"{path} is not a parity matrix (first column must be item)")
    return {line.split("\t", 1)[0] for line in lines[1:] if line.strip()}


def covers_owner(flag_set: str | None) -> str:
    """How an error names a `covers` list: the manifest's own, or that of the flag set `flag_set`."""
    return f"flag set {flag_set} covers" if flag_set is not None else "covers"


def covers_lists(manifest: dict) -> list[tuple[str, list]]:
    """Every `covers` list of a manifest with its owner: the manifest's own, then each flag set's."""
    lists = [(covers_owner(None), manifest.get("covers") or [])]
    for fs in manifest.get("flag_sets") or []:
        lists.append((covers_owner(fs["name"]), fs.get("covers") or []))
    return lists


def names_outside_matrix(names: Iterable[str], items: set[str], matrix: Path, owner: str) -> list[str]:
    """One error line per name that is not a row of the matrix, naming the list it came from."""
    return [f"{owner} names {n}, which is not a row of {matrix}" for n in names if n not in items]


def covers_names(args_names: list[str] | None) -> list[str]:
    """The `--covers` names of a command line, each non-empty, in order and without repeats."""
    names = list(dict.fromkeys(args_names or []))
    if not is_name_list(names):
        raise failure("--covers wants a non-empty row name")
    return names


def flag_set_named(manifest: dict, name: str) -> dict:
    for fs in manifest.get("flag_sets") or []:
        if fs["name"] == name:
            return fs
    raise failure(f"the manifest has no flag set {name!r}")


def path_inside_corpus(value: str, option: str) -> str:
    """A file or directory name relative to the corpus directory: no absolute path, no `..` segment."""
    if not is_corpus_path(value):
        raise failure(f"{option} wants a path inside the corpus directory, got {value!r}")
    return value


@dataclass
class NamedRecord:
    """A record the corpus carries outside the sampling: its suite, the label an error names it
    by (the `--include` argument, a kept record's first five columns, or a carried
    `included_records` name) and the Uploaded_variation names its reference rows can carry."""
    suite: str
    label: str
    names: tuple[str, ...]


def named_records(args: argparse.Namespace, existing: dict, suite_ids: set[str]) -> list[NamedRecord]:
    """The records `--include` names and, under `--keep-records`, every record of the manifest
    being regenerated and every name of its `included_records`, so a carried name is held to the
    checks a `--include` gets; a record of a suite that is not among `--suite` is refused."""
    named: list[NamedRecord] = []
    for arg in args.include or []:
        sid, sep, name = arg.partition("=")
        if not (sep and sid and name):
            raise failure(f"--include wants sid=name, got {arg!r}")
        named.append(NamedRecord(sid, f"--include {arg}", (name,)))
    if args.keep_records:
        if not existing:
            raise failure(f"--keep-records wants a manifest to regenerate in {args.out}")
        for rec in existing.get("records") or []:
            label = "kept record " + " ".join(rec["vcf"].split("\t")[:5])
            named.append(NamedRecord(rec["suite"], label, tuple(record_names(rec["vcf"]))))
        for sid, names in (existing.get("included_records") or {}).items():
            named.extend(NamedRecord(sid, f"carried include {sid}={name}", (name,)) for name in names)
    for n in named:
        if n.suite not in suite_ids:
            raise failure(f"{n.label} names suite {n.suite!r}, which is not among --suite")
    return named


def included_by_suite(includes: list[str], carried: dict[str, list[str]] | None = None) -> dict[str, list[str]]:
    """The `--include` arguments as suite -> sorted names, merged with the `included_records` a
    regeneration under `--keep-records` carries from the manifest it keeps: the new manifest's
    `included_records`."""
    by_suite: dict[str, set[str]] = defaultdict(set)
    for sid, names in (carried or {}).items():
        by_suite[sid].update(names)
    for arg in includes:
        sid, _sep, name = arg.partition("=")
        by_suite[sid].add(name)
    return {sid: sorted(names) for sid, names in sorted(by_suite.items())}


def cmd_select(args: argparse.Namespace) -> int:
    suites = [parse_suite_arg(s) for s in args.suite]
    out = Path(args.out)
    existing = read_manifest(out) if (out / "manifest.json").is_file() else {}
    named = named_records(args, existing, {sid for sid, _ref, _inp in suites})
    covers = list(dict.fromkeys([*(existing.get("covers") or []), *covers_names(args.covers)]))
    masks = existing.get("masks") or []
    reference_warnings = existing.get("reference_warnings") or []
    # A set's divergence keys name the old records by ordinal; `classify --flag-set` rebuilds them.
    flag_sets = [{k: v for k, v in fs.items() if k not in CLASSIFY_KEYS} for fs in existing.get("flag_sets") or []]
    if not covers:
        print("ERROR: [build_golden_corpus] select writes no manifest without --covers: "
              "name the parity-matrix rows the corpus exercises", file=sys.stderr)
        return 1
    if args.matrix is not None:
        matrix = Path(args.matrix)
        errors = names_outside_matrix(covers, matrix_items(matrix), matrix, covers_owner(None))
        for e in errors:
            print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
        if errors:
            return 1
    family_specs = parse_family_columns(args.family_columns)
    keep_samples = args.keep_samples
    required = genotype_class_names(args.require_genotype_class)
    if keep_samples < 0:
        raise failure(f"--keep-samples wants a count of sample columns, got {keep_samples}")
    if args.genotype_exemplars < 0:
        raise failure(f"--genotype-exemplars wants a count of records per class, got {args.genotype_exemplars}")
    if required and not keep_samples:
        raise failure("--require-genotype-class wants --keep-samples: a corpus without sample columns shows no genotype")
    names_by_suite: dict[str, set[str]] = defaultdict(set)
    for n in named:
        names_by_suite[n.suite].update(n.names)
    rng = random.Random(args.seed)
    reservoir, seen, per_suite, families, named_rows = reservoir_pass(suites, args.seed, family_specs, names_by_suite)
    unrowed = [n for n in named if not any((n.suite, name) in named_rows for name in n.names)]
    for n in unrowed:
        print(f"ERROR: [build_golden_corpus] {n.label} has no row in the {n.suite} reference", file=sys.stderr)
    if unrowed:
        return 1
    cost = CacheCost(Path(args.cache), args.flank).row_cost if args.cache else None
    chosen = choose_exemplars(reservoir, args.k, rng, args.max_span, args.focus_term or (), args.focus_k, cost)
    chosen = ensure_class_exemplars(chosen, reservoir, args.max_span, cost)
    if family_specs:
        chosen = ensure_family_exemplars(chosen, families, args.max_span, cost)
    for rows in named_rows.values():
        for r in rows:
            bucket = chosen.setdefault(r.consequence_set, [])
            if r not in bucket:
                bucket.append(r)
    wanted: dict[str, set[str]] = defaultdict(set)
    for combo, rows in chosen.items():
        for r in rows:
            wanted[r.suite].add(r.uploaded_variation)
    resolved = resolve_records(
        [(sid, inp) for sid, _ref, inp in suites], wanted, args.multi_allelic_per_suite, args.seed,
        keep_samples, args.genotype_exemplars,
    )
    records: dict[str, dict] = {}
    missing: list[str] = []
    for combo, rows in chosen.items():
        for r in rows:
            chrom = r.location.split(":", 1)[0]
            lines = resolved.found.get(r.suite, {}).get((r.uploaded_variation, chrom), [])
            if not lines:
                missing.append(f"{r.suite}\t{r.uploaded_variation}\t{r.location}")
                continue
            for line in lines:
                key = f"{r.suite}\t{line}"
                rec = records.setdefault(key, {"suite": r.suite, "vcf": line, "combinations": [], "exemplar_rows": []})
                if combo not in rec["combinations"]:
                    rec["combinations"].append(combo)
                rec["exemplar_rows"].append({
                    "uploaded_variation": r.uploaded_variation, "location": r.location, "allele": r.allele,
                    "feature": r.feature, "variant_class": r.klass, "strand": r.strand,
                    **({"family_stratum": r.family_stratum} if family_specs else {}),
                })
    for sid in sorted(resolved.multi):
        for stratum, lines in sorted(resolved.multi[sid].items()):
            for line in lines:
                key = f"{sid}\t{line}"
                rec = records.setdefault(key, {"suite": sid, "vcf": line, "combinations": [], "exemplar_rows": []})
                rec.setdefault("strata", []).append(stratum)
    if keep_samples:
        for rec in records.values():
            rec["genotype_classes"] = sorted(record_genotype_classes(rec["vcf"].split("\t"), keep_samples))
        showing = Counter(c for rec in records.values() for c in rec["genotype_classes"])
        for c in sorted(resolved.genotype_seen):
            for sid, line in resolved.genotype_pool[c]:
                if showing[c] >= args.genotype_exemplars:
                    break
                key = f"{sid}\t{line}"
                if key not in records:
                    classes = sorted(record_genotype_classes(line.split("\t"), keep_samples))
                    records[key] = {"suite": sid, "vcf": line, "combinations": [], "exemplar_rows": [], "genotype_classes": classes}
                    showing.update(classes)
    if len(records) > args.max_records:
        # The first record of every combination, family stratum and genotype class and every
        # named record survive the cap, which those alone may exceed; the rest fill in sorted
        # record order, so the cut is deterministic.
        by_group: dict[str, list[str]] = defaultdict(list)
        for key, rec in records.items():
            for c in rec["combinations"]:
                by_group[f"combination {c}"].append(key)
            for row in rec["exemplar_rows"]:
                if "family_stratum" in row:
                    by_group[f"stratum {row['family_stratum']}"].append(key)
            for c in rec.get("genotype_classes", ()):
                by_group[f"genotype class {c}"].append(key)
        keep = {keys[0] for keys in by_group.values()}
        for key, rec in records.items():
            suite_names = names_by_suite.get(rec["suite"], ())
            if any(row["uploaded_variation"] in suite_names for row in rec["exemplar_rows"]):
                keep.add(key)
        for key in sorted(records):
            if len(keep) >= args.max_records:
                break
            keep.add(key)
        records = {k: v for k, v in records.items() if k in keep}

    def sort_key(item):
        f = item[1]["vcf"].split("\t")
        chrom = f[0]
        chrom_rank = (0, int(chrom)) if chrom.isdigit() else (1, {"X": 1, "Y": 2, "MT": 3}.get(chrom, 9))
        return (chrom_rank, int(f[1]), f[3], f[4])

    ordered = sorted(records.items(), key=sort_key)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "variants.vcf").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(resolved.header + "\n")
        for _key, rec in ordered:
            fh.write(rec["vcf"] + "\n")
    covered = {c for _k, rec in ordered for c in rec["combinations"]}
    genotype_classes = Counter(c for _k, rec in ordered for c in rec.get("genotype_classes", ()))
    uncovered_genotype_classes = sorted((set(resolved.genotype_seen) | set(required)) - set(genotype_classes))
    carried = existing.get("included_records") if args.keep_records else None
    included = included_by_suite(args.include or [], carried)
    manifest = {
        "release": args.release,
        "assembly": args.assembly,
        "generator_version": GENERATOR_VERSION,
        "seed": args.seed,
        "exemplars_per_combination": args.k,
        **({"focus_terms": sorted(args.focus_term), "focus_exemplars_per_combination": args.focus_k} if args.focus_term else {}),
        **({"included_records": included} if included else {}),
        **({"kept_records": len(existing.get("records") or [])} if args.keep_records else {}),
        "max_span": args.max_span,
        **({"cache_cost_flank": args.flank} if args.cache else {}),
        **({"keep_samples": keep_samples, "samples": resolved.samples, "genotype_exemplars_per_class": args.genotype_exemplars,
            "required_genotype_classes": required} if keep_samples else {}),
        **run_options(args),
        "covers": covers,
        **({"masks": masks} if masks else {}),
        **({"reference_warnings": reference_warnings} if reference_warnings else {}),
        **({"flag_sets": flag_sets} if flag_sets else {}),
        # File names only: the directory each suite was staged in belongs to the build
        # host, not to the corpus, and the tests read neither field.
        "suites": {sid: {"reference": Path(ref).name, "input": Path(inp).name, "rows": sum(per_suite[sid].values()), "combinations": len(per_suite[sid])} for sid, ref, inp in suites},
        "combinations": {c: {"observed_rows": seen[c], "exemplar_records": sum(1 for _k, rec in ordered if c in rec["combinations"])} for c in sorted(seen)},
        "variant_classes": dict(sorted(Counter(r.klass for rows in chosen.values() for r in rows).items())),
        **({"family_columns": family_specs,
            "family_strata": dict(sorted(Counter(r.family_stratum for rows in chosen.values() for r in rows).items()))}
           if family_specs else {}),
        **({"genotype_classes": dict(sorted(genotype_classes.items()))} if keep_samples else {}),
        "records": [
            {"suite": rec["suite"], "vcf": rec["vcf"], "combinations": sorted(rec["combinations"]),
             **({"strata": rec["strata"]} if rec.get("strata") else {}),
             **({"genotype_classes": rec["genotype_classes"]} if keep_samples else {})}
            for _k, rec in ordered
        ],
        "multi_allelic_per_suite": args.multi_allelic_per_suite,
        "uncovered_combinations": sorted(set(seen) - covered),
        **({"uncovered_genotype_classes": uncovered_genotype_classes} if keep_samples else {}),
        "unresolved_exemplars": missing,
    }
    write_manifest(out, manifest)
    classes = manifest["variant_classes"]
    class_counts = ", ".join(f"{k} {v}" for k, v in classes.items())
    genotype_note = (f", genotype classes observed {len(resolved.genotype_seen)}, covered {len(genotype_classes)}"
                     if keep_samples else "")
    # A record both included and kept is two entries of `named` and one key of `named_rows`.
    named_note = f", named records {len(named_rows)}" if named else ""
    print(f"combinations observed {len(seen)}, covered {len(covered)}, records {len(ordered)}, unresolved exemplars {len(missing)}, "
          f"variant classes {len(classes)} ({class_counts}){genotype_note}{named_note}")
    return 0 if not missing and not manifest["uncovered_combinations"] and not uncovered_genotype_classes else 1


# --------------------------------------------------------------------------- classify
def run_masks(manifest: dict, target: dict) -> list[str]:
    """The mask classes of a run, as the golden harness reads them: the flag set's own `masks`
    when its entry (`target`) carries the key, else the manifest's; `target` is the manifest
    itself for the default run."""
    if "masks" in target:
        return list(target["masks"])
    return list(manifest.get("masks") or [])


def masked_feature_types(masks: Iterable[str]) -> frozenset[str]:
    """The `Feature_type` values a run's `masks` leave out of the reference: the regulatory
    rows' under `regulatory_feature_rows`; `order_only` leaves out no row (`fields_agree`
    gives it its meaning)."""
    return REGULATORY_FEATURE_TYPES if REGULATORY_FEATURE_ROWS in masks else frozenset()


def fields_agree(name: str, vep: str, vep_rs: str, order_only: bool) -> bool:
    """Whether two values of the Extra field `name` agree: exactly, or, for a column of
    `ORDER_ONLY_COLUMNS` under `order_only`, as a multiset of their comma-separated items;
    an absent value (`-`) agrees with itself alone."""
    if order_only and name in ORDER_ONLY_COLUMNS and "-" not in (vep, vep_rs):
        return Counter(vep.split(",")) == Counter(vep_rs.split(","))
    return vep == vep_rs


def load_tuples(path: Path, masked: frozenset[str] = frozenset()) -> tuple[dict[tuple[str, str, str, str], set[str]], dict[tuple[str, str, str, str], set[str]]]:
    """(Location, Allele, Feature, Feature_type) -> consequence-set strings, and -> Uploaded_variation
    names, leaving out the rows whose Feature_type is in `masked`."""
    out: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    names: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    for r in iter_rows(path, "x"):
        if r.feature_type in masked:
            continue
        key = (r.location, r.allele, r.feature, r.feature_type)
        out[key].add(r.consequence_set)
        names[key].add(r.uploaded_variation)
    return out, names


def load_fields(path: Path, fields: list[str], masked: frozenset[str] = frozenset()) -> dict[tuple[str, str, str, str, str], dict[str, str]]:
    """(Location, Allele, Feature, Feature_type, consequence set) -> the named Extra fields ('-' when absent),
    leaving out the rows whose Feature_type is in `masked`."""
    out: dict[tuple[str, str, str, str, str], dict[str, str]] = {}
    with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            while len(f) < 14:
                f.append("-")
            if f[5] in masked:
                continue
            extra = parse_extra(f[13])
            key = (f[1], f[2], f[4], f[5], normalize_consequence_set(f[6]))
            out.setdefault(key, {name: extra.get(name, "-") for name in fields})
    return out



def record_names(vcf_line: str) -> list[str]:
    """The Uploaded_variation names a VCF record line can carry: its ID, else VEP's auto-names."""
    f = vcf_line.split("\t")
    chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
    return [vid.split(";")[0]] if vid != "." else vep_auto_names(chrom, pos, ref, alt.split(","))


def record_name_index(records: list[dict]) -> dict[str, list[int]]:
    """Uploaded_variation (the ID, else VEP's auto-name) -> ordinals of the corpus records it names."""
    index: dict[str, list[int]] = defaultdict(list)
    for i, rec in enumerate(records):
        for n in record_names(rec["vcf"]):
            index[n].append(i)
    return index


REFERENCE_SKIPPED = "reference_skipped_record"
UNEXPLAINED = "unexplained_residual"
# The flags that drop a record's rows by what they overlap (OutputFactory.pm:482-500 and
# :540-558), so a record of a run under one may have no row on either side.
ROW_FILTERS = frozenset({"--coding_only", "--no_intergenic"})
# A VCF REF column: bases or `.`, never the `/`-joined allele string or the class word of an
# Ensembl-format line, so a line of five or more columns with one is a VCF record.
VCF_REF = re.compile(r"[ACGTNacgtn.]+")


def first_five_columns(line: str) -> str:
    """A VCF line's CHROM, POS, ID, REF and ALT, the columns a twin of the line (another
    encoding, sample columns added or dropped) shares with it."""
    return "\t".join(line.split("\t")[:5])


def carried_records(records: list[dict], json_path: Path) -> set[int]:
    """Ordinals of the records the reference carried with no row: the JSON output at
    `json_path`, when the file exists, has an object for the record whose
    `most_severe_consequence` is `?`, the object VEP writes for a feature with no overlap
    allele (every ALT a `*` or `<DEL:*>`, equal to the reference or a repeat of an earlier
    one; VariationFeatureOverlap.pm:459-491) and for a record a row filter (`--coding_only`,
    `--no_intergenic`) left no consequence of (JSON.pm:387-388). A record is the object whose
    `input` opens with its first five columns, so a twin input carrying sample columns names
    it too. A file that is not gzip, or holds a line that is not a JSON object, is refused
    with its path."""
    if not json_path.is_file():
        return set()
    by_head: dict[str, list[int]] = defaultdict(list)
    for i, rec in enumerate(records):
        by_head[first_five_columns(rec["vcf"])].append(i)
    carried: set[int] = set()
    try:
        with gzip.open(json_path, "rt", encoding="utf-8", errors="replace") as fh:
            for number, line in enumerate(fh, 1):
                if not line.strip():
                    continue
                obj = json.loads(line)
                if not isinstance(obj, dict):
                    raise failure(f"{json_path} is not a readable JSON output (line {number} is not a JSON object)")
                if obj.get("most_severe_consequence") != "?":
                    continue
                head = first_five_columns(str(obj.get("input", "")))
                carried.update(by_head.get(head, []))
    except (OSError, EOFError, zlib.error, json.JSONDecodeError) as exc:
        raise failure(f"{json_path} is not a readable JSON output ({exc})") from None
    return carried


def run_records(corpus: Path, manifest: dict, target: dict) -> list[dict]:
    """The records of a run, numbered as the golden harness numbers them (every line of the
    run's input that is not a header, in order): the manifest's `records` for a run over
    `variants.vcf`; for a flag set over another `input`, that file's lines. A line sharing
    its first five columns with a manifest record is that record (a twin of `variants.vcf`
    in another encoding keeps its `reference_warning` and takes `reference_rows`). A file
    holding a VCF record (`vcf_record`: five tab-separated columns or more, a POS of digits,
    a REF of bases) is a VCF whose other lines are records of their own, which carry no
    `reference_warning`, so a record the reference skipped in such an input is refused as
    unexplained. A file holding none is `variants.vcf` re-encoded in another format or
    delimiter (`--format hgvs`, an Ensembl-format line, a space-delimited VCF line), record
    for record, and each line is the manifest record at its ordinal; such a file whose
    record count is not the manifest's is refused with its path, as is a line of a VCF input
    that is neither a twin nor a VCF record, and a file that cannot be read."""
    input_name = target.get("input") or DEFAULT_INPUT
    if input_name == DEFAULT_INPUT:
        return manifest["records"]
    by_head = {first_five_columns(rec["vcf"]): rec for rec in manifest["records"]}
    path = corpus / input_name
    try:
        with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
            lines = [line.rstrip("\n") for line in handle if not line.startswith("#")]
    except (OSError, EOFError, zlib.error) as exc:
        raise failure(f"{path} is not a readable text file ({exc})") from None

    def vcf_record(line: str) -> bool:
        columns = line.split("\t")
        return len(columns) >= 5 and columns[1].isdecimal() and VCF_REF.fullmatch(columns[3]) is not None

    if any(vcf_record(line) for line in lines):
        records: list[dict] = []
        for i, line in enumerate(lines):
            manifest_record = by_head.get(first_five_columns(line))
            if manifest_record is not None:
                records.append(manifest_record)
            elif vcf_record(line):
                records.append({"vcf": line})
            else:
                raise failure(f"{path} record {i} ({line[:40]!r}) is no VCF record in a VCF input")
        return records
    if len(lines) != len(manifest["records"]):
        raise failure(f"{path} holds no VCF record and its {len(lines)} records are not the manifest's "
                      f"{len(manifest['records'])} re-encoded")
    return manifest["records"]


def reference_skipped_records(
    records: list[dict],
    perl_names: dict[tuple[str, str, str, str], set[str]],
    by_name: dict[str, list[int]],
    carried: set[int] = frozenset(),
) -> tuple[set[int], list[str]]:
    """Ordinals of the records no VEP row names (VEP dropped them before annotation), and the
    records that contradict their `reference_warning`: a zero-row record without one (its name
    unresolved, or the reason unrecorded) and a warned record with rows. A record that carries
    no row by its nature is kept, not skipped, when it has no warning: a non-variant record (ALT
    `.`, written under `--allow_non_variant` as its VCF line and JSON object and no default row,
    Parser/VCF.pm:259-267) and a record of `carried`, one the reference's JSON names as a
    feature with no allele to annotate (`carried_records`) or one a row filter of the run
    left without a row on either side (`cmd_classify`)."""
    with_rows = {i for names in perl_names.values() for n in names for i in by_name.get(n, [])}
    skipped = set(range(len(records))) - with_rows
    problems: list[str] = []
    for i, rec in enumerate(records):
        warned = isinstance(rec.get("reference_warning"), str)
        columns = rec["vcf"].split("\t")
        label = " ".join(columns[:5])
        non_variant = len(columns) > 4 and columns[4] == "."
        if i in skipped and not warned:
            if non_variant or i in carried:
                skipped.discard(i)
            else:
                problems.append(f"record {i} ({label}) has no VEP row and no reference_warning")
        elif warned and i not in skipped:
            problems.append(f"record {i} ({label}) has VEP rows and reference_warning {rec['reference_warning']!r}")
    return skipped, problems


def classify_pair(perl_csq: str, rust_csq: str, allele: str) -> str | None:
    """First matching rule of the key comparator's registry, by taxonomy class."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "concordance"))
    import compare_vep_outputs as cvo  # noqa: E402

    for rule in cvo.EXCLUSION_REGISTRY:
        if rule.matches(perl_csq, rust_csq, allele):
            return rule.taxonomy_class
    return None


def classify_mate_side_pair(
    location: str, allele: str, feature: str, perl_csq: str, rust_csq: str, models: dict
) -> str | None:
    """`breakend_mate_local_read` when VEP's set is the local-coordinate read of a point
    breakend's mate-side Transcript row and vep-rs's set is the mate-coordinate read, the
    two differing; the derivation is the SV comparator's `breakend_mate_context`, over the
    corpus's own cache."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc  # noqa: E402

    mate = bmc.bracket_mate(allele)
    if mate is None:
        return None
    model = models.get((mate[0], feature))
    if model is None:
        return None
    derived = bmc.mate_side_derivation(location, allele, model)
    if derived is None:
        return None
    mate_set, local_set = derived
    if mate_set == local_set:
        return None
    if set(perl_csq.split(",")) == local_set and set(rust_csq.split(",")) == mate_set:
        return "breakend_mate_local_read"
    return None


def load_corpus_models(corpus: Path) -> dict:
    """The corpus's pruned cache, located by the manifest's `cache_layout`, as
    `breakend_mate_context` transcript models."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc  # noqa: E402

    manifest = read_manifest(corpus)
    cache = corpus / "json_cache"
    if manifest.get("cache_layout") == "ensembl":
        cache = corpus / "cache_root" / "homo_sapiens" / f"{manifest['release']}_{manifest['assembly']}"
    return bmc.load_transcript_models(cache)


def cmd_classify(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = read_manifest(corpus)
    target = flag_set_named(manifest, args.flag_set) if args.flag_set is not None else manifest
    if args.flag_set is not None and run_options(args):
        print("ERROR: [build_golden_corpus] --flag and --fasta-name describe the default run; "
              "a flag set's flags are its own entry's", file=sys.stderr)
        return 1
    masks = run_masks(manifest, target)
    masked = masked_feature_types(masks)
    perl, perl_names = load_tuples(Path(args.vep_default), masked)
    rust, rust_names = load_tuples(Path(args.vep_rs_output))
    records = run_records(corpus, manifest, target)
    by_name = record_name_index(records)
    # The default run's JSON stands in for a flag set whose expected directory has none:
    # `carried_records` names a record by its first five columns, which a twin input shares.
    reference_json = corpus / "expected" / "json.jsonl.gz"
    if args.flag_set is not None:
        expected_dir = target.get("expected_dir", f"expected_{args.flag_set}")
        own_json = corpus / expected_dir / "json.jsonl.gz"
        if own_json.is_file():
            reference_json = own_json
    rowless = carried_records(records, reference_json)
    if args.flag_set is not None:
        flags = target.get("flags") or []
    else:
        flags = run_options(args).get("flags") or manifest.get("flags") or []
    if ROW_FILTERS.intersection(flags):
        # A row filter leaves a record with no coding (or no intergenic) row rowless on both
        # sides; the vep-rs side is read here, the reference side by the skip rule.
        with_rust_rows = {i for names in rust_names.values() for n in names for i in by_name.get(n, [])}
        rowless |= set(range(len(records))) - with_rust_rows
    skipped, problems = reference_skipped_records(records, perl_names, by_name, rowless)
    if problems:
        for p in problems:
            print(f"ERROR: [build_golden_corpus] {p}", file=sys.stderr)
        return 1
    for i, rec in enumerate(records):
        rec.pop("reference_rows", None)
        if i in skipped:
            # The two markers are the record's last two keys, in this order.
            rec["reference_rows"] = 0
            rec["reference_warning"] = rec.pop("reference_warning")
    models = load_corpus_models(corpus)

    def records_for(key: tuple[str, str, str, str]) -> list[int]:
        found: set[int] = set()
        for n in perl_names.get(key, set()) | rust_names.get(key, set()):
            found.update(by_name.get(n, []))
        return sorted(found)

    divergences: list[dict] = []
    for key, perl_sets in sorted(perl.items()):
        rust_sets = rust.get(key, set())
        for cs in sorted(perl_sets - rust_sets):
            klass = None
            for rs in sorted(rust_sets):
                klass = classify_pair(cs, rs, key[1]) or classify_mate_side_pair(
                    key[0], key[1], key[2], cs, rs, models
                )
                if klass:
                    break
            if klass is None and not rust_sets:
                klass = "vep_only_transcript" if variant_class(key[0], key[1]) not in SEQUENCE_CLASSES else None
            divergences.append({
                "location": key[0], "allele": key[1], "feature": key[2], "feature_type": key[3],
                "record_indices": records_for(key),
                "vep_consequence_set": cs, "vep_rs_consequence_sets": sorted(rust_sets),
                "expected_divergence": klass or UNEXPLAINED,
            })
    extra_rust: list[dict] = []
    for k, v in sorted(rust.items()):
        if k in perl:
            continue
        indices = records_for(k)
        extra_rust.append({
            "location": k[0], "allele": k[1], "feature": k[2], "feature_type": k[3],
            "record_indices": indices, "vep_rs_consequence_sets": sorted(v),
            "expected_divergence": REFERENCE_SKIPPED if indices and skipped.issuperset(indices) else UNEXPLAINED,
        })
    # Classify's own keys are rewritten from scratch, in one order, so a key a previous run
    # wrote with other inputs (another `--field` list) does not survive this one.
    for key in CLASSIFY_KEYS:
        target.pop(key, None)
    target["divergences"] = divergences
    target["vep_rs_only_tuples"] = extra_rust
    target["divergence_summary"] = dict(Counter(d["expected_divergence"] for d in divergences))
    target["vep_rs_only_summary"] = dict(Counter(t["expected_divergence"] for t in extra_rust))
    field_note = ""
    if args.field:
        # Field values are compared only on keys whose consequence sets agree; a key whose
        # terms differ is already a documented (or unexplained) divergence above and the
        # golden test skips its fields, which is where the start co-emission keys land (VEP
        # prints `Met1?` wherever its start_lost predicate holds, vep-rs the peptide change
        # wherever it keeps start_retained_variant). No field difference on an agreeing key
        # is documented, so each one is an `unexplained_residual` the golden test fails on.
        perl_fields = load_fields(Path(args.vep_default), args.field, masked)
        rust_fields = load_fields(Path(args.vep_rs_output), args.field)
        field_divergences: list[dict] = []
        for key, pvals in sorted(perl_fields.items()):
            rvals = rust_fields.get(key)
            if rvals is None:
                continue
            for name in args.field:
                if not fields_agree(name, pvals[name], rvals[name], ORDER_ONLY in masks):
                    field_divergences.append({
                        "location": key[0], "allele": key[1], "feature": key[2], "feature_type": key[3],
                        "record_indices": records_for(key[:4]),
                        "consequence_set": key[4], "field": name, "vep": pvals[name], "vep_rs": rvals[name],
                        "expected_divergence": UNEXPLAINED,
                    })
        target["fields_compared"] = list(args.field)
        target["field_divergences"] = field_divergences
        target["field_divergence_summary"] = dict(Counter(d["expected_divergence"] for d in field_divergences))
        field_note = f"; field divergences {len(field_divergences)} ({target['field_divergence_summary']})"
    manifest.update(run_options(args))
    write_manifest(corpus, manifest)
    print(f"divergent VEP rows {len(divergences)} ({target['divergence_summary']}); "
          f"vep-rs-only tuples {len(extra_rust)} ({target['vep_rs_only_summary']}); "
          f"reference-skipped records {sorted(skipped)}{field_note}")
    return 0


# ------------------------------------------------------------------------------ check
def sample_column_errors(corpus: Path, manifest: dict) -> list[str]:
    """One line per way `variants.vcf` fails to carry the manifest's kept samples: a `#CHROM`
    line whose columns are not the eight fixed ones, `FORMAT` and the manifest's `samples`, or
    a record whose column count differs from that header's; a file that cannot be read is
    refused by `vcf_lines`."""
    keep = manifest["keep_samples"]
    expected = [*FIXED_COLUMNS, "FORMAT", *manifest["samples"]]
    errors: list[str] = []
    has_header = False
    for is_header, f in vcf_lines(corpus / "variants.vcf"):
        if is_header:
            has_header = True
            if f != expected:
                errors.append(f"variants.vcf header names {' '.join(f[8:]) or 'no sample columns'}, the manifest keeps FORMAT {' '.join(manifest['samples'])}")
        elif len(f) != 9 + keep:
            errors.append(f"variants.vcf record {' '.join(f[:5])} has {len(f)} columns, the header {9 + keep}")
    if not has_header:
        errors.append("variants.vcf has no #CHROM header line")
    return errors


def genotype_coverage(
    suites: list[tuple[str, Path, Path]], manifest: dict, required: list[str]
) -> tuple[set[str], set[str], list[str]]:
    """The genotype classes the suite inputs show among the manifest's kept samples, the ones
    the corpus records show, and those of the first set or of `required` the corpus lacks."""
    keep = manifest["keep_samples"]
    observed: set[str] = set()
    for _sid, _ref, inp in suites:
        samples: list[str] = []
        for is_header, f in vcf_lines(inp, 9 + keep):
            if is_header:
                samples = kept_samples(f, keep, inp, manifest["samples"], "the manifest")
            else:
                observed |= record_genotype_classes(f, keep)
        if not samples:
            raise failure(f"{inp} has no #CHROM header line naming its samples")
    covered = {c for rec in manifest["records"] for c in record_genotype_classes(rec["vcf"].split("\t"), keep)}
    return observed, covered, sorted((observed | set(required)) - covered)


def cmd_check(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = read_manifest(corpus)
    errors: list[str] = []
    if args.matrix is not None:
        matrix = Path(args.matrix)
        items = matrix_items(matrix)
        for owner, names in covers_lists(manifest):
            errors += names_outside_matrix(names, items, matrix, owner)
    suites = [parse_suite_arg(s) for s in args.suite]
    family_specs = parse_family_columns(args.family_columns)
    required = genotype_class_names(args.require_genotype_class)
    if required and not manifest.get("keep_samples"):
        raise failure("--require-genotype-class wants a corpus built with --keep-samples")
    by_name = record_name_index(manifest["records"])
    records = manifest["records"]

    def is_corpus_row(r: Row) -> bool:
        chrom = r.location.split(":", 1)[0]
        return any(records[i]["suite"] == r.suite and records[i]["vcf"].split("\t", 1)[0] == chrom
                   for i in by_name.get(r.uploaded_variation, []))

    observed: set[str] = set()
    strata_observed: set[str] = set()
    strata_covered: set[str] = set()
    for sid, ref, _inp in suites:
        for r in iter_rows(ref, sid, family_specs):
            observed.add(r.consequence_set)
            if family_specs:
                strata_observed.add(r.family_stratum)
                if is_corpus_row(r):
                    strata_covered.add(r.family_stratum)
    covered = {c for rec in records for c in rec["combinations"]}
    missing = sorted(observed - covered)
    print(f"observed {len(observed)}, covered by corpus {len(covered & observed)}, missing {len(missing)}")
    for m in missing[:50]:
        print("  MISSING", m)
    missing_strata = sorted(strata_observed - strata_covered)
    if family_specs:
        print(f"family strata observed {len(strata_observed)}, covered by corpus {len(strata_covered)}, "
              f"missing {len(missing_strata)}")
        for m in missing_strata[:50]:
            print("  MISSING family stratum", m)
    missing_classes: list[str] = []
    if manifest.get("keep_samples"):
        try:
            errors += sample_column_errors(corpus, manifest)
            observed_classes, covered_classes, missing_classes = genotype_coverage(
                suites, manifest, [*manifest.get("required_genotype_classes", []), *required]
            )
        except SystemExit:
            # The refusal that ends the run comes after the errors found before it.
            for e in errors:
                print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
            raise
        print(f"genotype classes observed {len(observed_classes)}, covered by corpus {len(covered_classes)}, "
              f"missing {len(missing_classes)}")
        for m in missing_classes:
            print("  MISSING genotype class", m)
    for e in errors:
        print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
    return 1 if missing or missing_strata or missing_classes or errors else 0


def cmd_covers(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = read_manifest(corpus)
    target = flag_set_named(manifest, args.flag_set) if args.flag_set is not None else manifest
    covers = list(target.get("covers") or [])
    if args.add:
        added = covers_names(args.add)
        matrix = Path(args.matrix)
        errors = names_outside_matrix(added, matrix_items(matrix), matrix, covers_owner(args.flag_set))
        for e in errors:
            print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
        if errors:
            return 1
        covers = list(dict.fromkeys([*covers, *added]))
        target["covers"] = covers
        write_manifest(corpus, manifest)
    if args.list or not args.add:
        for name in covers:
            print(name)
    return 0


def cmd_flag_set(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = read_manifest(corpus)
    sets = manifest.setdefault("flag_sets", [])
    if args.list:
        for fs in sets:
            columns = [
                fs["name"],
                " ".join(fs.get("flags") or []),
                fs.get("input", DEFAULT_INPUT),
                str(bool(fs.get("stdin"))).lower(),
                str(bool(fs.get("unordered"))).lower(),
                ",".join(fs.get("formats") or OUTPUT_FORMATS),
                fs.get("expected_dir", f"expected_{fs['name']}"),
                ",".join(fs.get("covers") or []),
            ]
            print("\t".join(columns))
        return 0
    if not args.name:
        print("ERROR: [build_golden_corpus] flag-set wants --name (or --list)", file=sys.stderr)
        return 1
    if not is_flag_set_name(args.name):
        print(f"ERROR: [build_golden_corpus] flag set name {args.name!r}: letters, digits, _ and -, "
              "and not the manifest's own run, default", file=sys.stderr)
        return 1
    formats = [f.strip() for f in (args.format or ",".join(OUTPUT_FORMATS)).split(",") if f.strip()]
    bad = [f for f in formats if f not in OUTPUT_FORMATS]
    if bad or not formats:
        print(f"ERROR: [build_golden_corpus] --format wants a non-empty subset of {','.join(OUTPUT_FORMATS)}, "
              f"got {args.format!r}", file=sys.stderr)
        return 1
    flags = list(args.flag or [])
    if args.stdin and not names_format(flags):
        # Ensembl VEP detects no format on STDIN and dies without --format, so a stdin run's
        # reference exists only for a flag list that names one.
        print("ERROR: [build_golden_corpus] a stdin run names its input format: add --flag=--format --flag=<format>",
              file=sys.stderr)
        return 1
    path_inside_corpus(args.input, "--input")
    if not (corpus / args.input).is_file():
        print(f"ERROR: [build_golden_corpus] input {args.input} is not a file in {corpus}", file=sys.stderr)
        return 1
    expected_dir = path_inside_corpus(args.expected_dir or f"expected_{args.name}", "--expected-dir")
    covers = covers_names(args.covers)
    if covers:
        matrix = Path(args.matrix)
        errors = names_outside_matrix(covers, matrix_items(matrix), matrix, covers_owner(args.name))
        for e in errors:
            print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
        if errors:
            return 1
    written = {
        "name": args.name,
        "flags": flags,
        "input": args.input,
        "stdin": bool(args.stdin),
        "unordered": bool(args.unordered),
        "formats": formats,
        "expected_dir": expected_dir,
        "covers": covers,
    }
    if args.perl_mode != "offline":
        written["perl_mode"] = args.perl_mode
    for i, fs in enumerate(sets):
        if fs["name"] == args.name:
            # The options set their keys where they sit and the entry's other keys (`masks`,
            # `reference_warnings`, any hand-written key) stay. The divergence keys go, since
            # they describe the run before the rewrite (`classify --flag-set` rewrites them), and
            # so does `perl_mode`, which `written` carries only for a database reference.
            entry = {k: v for k, v in fs.items() if k not in CLASSIFY_KEYS and k != "perl_mode"}
            entry.update(written)
            sets[i] = entry
            break
    else:
        entry = written
        sets.append(entry)
    write_manifest(corpus, manifest)
    print(f"flag set {args.name}: {len(entry['flags'])} flags, input {entry['input']}"
          f"{' on stdin' if entry['stdin'] else ''}{', unordered' if entry['unordered'] else ''}"
          f"{', reference from the database' if 'perl_mode' in entry else ''}, "
          f"formats {','.join(formats)}, expected {entry['expected_dir']}, "
          f"covers {len(covers)}")
    return 0


def add_run_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--flag", action="append",
                        help="a command-line flag both engines ran with, written --flag=--hgvs (repeatable), recorded as `flags`")
    parser.add_argument("--fasta-name", default=None,
                        help="file name of the gzipped reference beside the manifest, recorded as `fasta`")


MATRIX_HELP = "the parity matrix whose `item` column every `covers` name must be in"
DEFAULT_MATRIX = "tests/parity/matrix.tsv"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--release", required=True)
    s.add_argument("--assembly", required=True)
    s.add_argument("--suite", action="append", required=True, help="sid=reference.txt=input.vcf[.gz]")
    s.add_argument("--out", required=True)
    s.add_argument("--k", type=int, default=3)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--max-records", type=int, default=1500)
    s.add_argument("--multi-allelic-per-suite", type=int, default=4,
                   help="multi-allelic sequence records sampled per suite and stratum (SNV-only, indel)")
    s.add_argument("--max-span", type=int, default=1_000_000,
                   help="prefer exemplars no wider than this many bases; 0 disables the preference")
    s.add_argument("--focus-term", action="append",
                   help="a consequence term; combinations containing it take --focus-k exemplars")
    s.add_argument("--focus-k", type=int, default=None, help="exemplars per focus combination")
    s.add_argument("--include", action="append", metavar="SID=NAME",
                   help="a record of suite SID's input, by its ID or VEP's name for an ID-less record, carried whatever "
                        "the sampling chooses (repeatable), recorded as `included_records`")
    s.add_argument("--keep-records", action="store_true",
                   help="keep every record of the manifest being regenerated, resolved and rowed again from its suite; "
                        "the count is recorded as `kept_records`")
    s.add_argument("--cache", default=None,
                   help="the JSON cache the fixture cache is pruned from; exemplars whose transcripts add fewer "
                        "bytes to it come first, ahead of the flag and position preferences")
    s.add_argument("--flank", type=int, default=5000,
                   help="the transcript flank the pruner keeps around each record, in bases (its --flank)")
    s.add_argument("--covers", action="append",
                   help="a parity-matrix row the corpus exercises (repeatable), recorded as `covers`; "
                        "required unless the manifest being regenerated already carries a list")
    s.add_argument("--family-columns", default=None,
                   help="Extra columns of the family key, comma-separated: NAME stratifies by presence, NAME= by value; "
                        "every stratum seen keeps an exemplar")
    s.add_argument("--keep-samples", type=int, default=0,
                   help="keep FORMAT and the first N sample columns of every record and of the #CHROM header "
                        "(default: none), and sample records by the genotype classes those samples show")
    s.add_argument("--genotype-exemplars", type=int, default=8,
                   help="records showing each genotype class the inputs show among the kept samples, when they have that many")
    s.add_argument("--require-genotype-class", action="append", metavar="CLASS",
                   help=f"a genotype class the corpus must carry, shown by the inputs or not (repeatable; one of {', '.join(GENOTYPE_CLASSES)})")
    s.add_argument("--matrix", default=None, help=MATRIX_HELP)
    add_run_options(s)
    s.set_defaults(fn=cmd_select)
    c = sub.add_parser("classify")
    c.add_argument("--corpus", required=True)
    c.add_argument("--vep-default", required=True, help="VEP default-format output for variants.vcf")
    c.add_argument("--vep-rs-output", required=True, help="vep-rs default-format output for variants.vcf")
    c.add_argument("--field", action="append",
                   help="an Extra field (HGVSc, HGVSp) compared on every key whose consequence sets agree")
    c.add_argument("--flag-set", default=None,
                   help="write the divergence keys into this flag set's entry instead of the manifest's top level")
    add_run_options(c)
    c.set_defaults(fn=cmd_classify)
    k = sub.add_parser("check")
    k.add_argument("--corpus", required=True)
    k.add_argument("--suite", action="append", required=True)
    k.add_argument("--matrix", default=None, help=MATRIX_HELP)
    k.add_argument("--family-columns", default=None,
                   help="the family key (see select); every stratum the suites show must have an exemplar")
    k.add_argument("--require-genotype-class", action="append", metavar="CLASS",
                   help="a genotype class the corpus must carry, beside the manifest's required_genotype_classes (repeatable)")
    k.set_defaults(fn=cmd_check)
    v = sub.add_parser("covers")
    v.add_argument("corpus")
    v.add_argument("--add", action="append", help="a parity-matrix row to add (repeatable)")
    v.add_argument("--list", action="store_true", help="print the list, one name per line")
    v.add_argument("--flag-set", default=None, help="the flag set whose list is read or extended")
    v.add_argument("--matrix", default=DEFAULT_MATRIX, help=f"{MATRIX_HELP} (default: {DEFAULT_MATRIX})")
    v.set_defaults(fn=cmd_covers)
    f = sub.add_parser(
        "flag-set",
        description="Write one flag_sets entry of the corpus manifest: a run of the corpus under its own flags. "
                    "A set that exists has the keys the options below set rewritten and keeps every other key "
                    "of its entry (its hand-written `masks` and `reference_warnings` among them); its divergence "
                    "keys go, so `classify --flag-set` runs after it.",
    )
    f.add_argument("corpus")
    f.add_argument("--name", default=None, help="the set's name; `expected_<name>` holds its outputs unless --expected-dir")
    f.add_argument("--flag", action="append", help="a command-line flag of the run, written --flag=--pick (repeatable)")
    f.add_argument("--input", default=DEFAULT_INPUT, help=f"the input file inside the corpus (default: {DEFAULT_INPUT})")
    f.add_argument("--stdin", action="store_true", help="pipe the input on standard input and pass -i STDIN")
    f.add_argument("--unordered", action="store_true",
                   help="the reference's row order within a record is Ensembl VEP's hash order (--per_gene gene "
                        "buckets, --pick_allele allele buckets), so the runner compares the entry set and not the order")
    f.add_argument("--perl-mode", default="offline", choices=PERL_MODES,
                   help="how the set's reference was generated: offline (Ensembl VEP on its cache, the default) "
                        "or database (connected to the Ensembl database, for a column Ensembl VEP writes only with "
                        "one; the default and tab headers then compare by their column descriptions alone)")
    f.add_argument("--format", default=None, help=f"the output formats compared, comma-separated (default: all of {','.join(OUTPUT_FORMATS)})")
    f.add_argument("--expected-dir", default=None, help="the directory of the set's expected outputs")
    f.add_argument("--covers", action="append", help="a parity-matrix row the set exercises (repeatable)")
    f.add_argument("--matrix", default=DEFAULT_MATRIX, help=f"{MATRIX_HELP} (default: {DEFAULT_MATRIX})")
    f.add_argument("--list", action="store_true", help="print the sets, one tab-separated line each")
    f.set_defaults(fn=cmd_flag_set)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
