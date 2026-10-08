#!/usr/bin/env python3
"""Compare Perl VEP and vep-rs outputs across a directory of SV VCFs.

The default input is the synthetic corpus under ``tests/sv_validation/``; a
real-world SV set is scored the same way through ``--input-dir``. Produces
per-variant-type concordance metrics with JSON, Markdown, and TSV reports of
discordant annotation tuples. ``discordant.tsv`` lists every one-sided tuple;
``discordant_adjusted.tsv``, same columns, lists the ones every mask leaves in
place, so its row count is the adjusted row's ``only_perl + only_rust``.

Usage:
    python3 scripts/validation/compare_sv_concordance.py \
      --input-dir tests/sv_validation/ \
      --perl-dir tmp/sv_perl/ \
      --rust-dir tmp/sv_rust/ \
      --output-dir tmp/sv_concordance/

Tuple key and normalisation
---------------------------
A tuple is ``(Location, Allele, Feature, Feature_type, Consequence)``. Before two
rows are compared, each column is normalised exactly as far as the two engines are
known to differ in REPRESENTATION while agreeing in MEANING, and no further:

* Location: the ``chr`` prefix is stripped and ``M`` becomes ``MT``
  (`normalize_location`). Perl preserves the input's chromosome spelling; vep-rs
  writes the Ensembl form.
* Allele: the ``chr`` prefix inside breakend brackets is stripped
  (`_strip_chr_in_bracket`) for every engine. With ``--scored-engine fastvep`` the
  scored side's raw symbolic tokens (``<DEL>``, ``<DUP:TANDEM>``, ``<INS:ME:ALU>``,
  ``<CN0>``, ...) are additionally mapped to the Sequence Ontology class strings
  Perl VEP and vep-rs write in that column (``deletion``, ``tandem_duplication``,
  ``Alu_insertion``, ``deletion``, ...); see `fastvep_symbolic_allele_to_so_class`.
  The Perl side and the vep-rs side are never mapped.
* Consequence: the comma-separated term set is sorted and de-duplicated
  (`normalize_consequence_set`).

Feature and Feature_type are compared verbatim.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

# The SNP/indel comparator's exclusion registry lives beside this script's directory;
# the breakend mate-context derivation shares this one.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "concordance"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from breakend_mate_context import (  # noqa: E402
    TranscriptModel,
    iter_cache_transcripts,
    mate_side_derivation,
    model_from_cache_record,
)
from compare_vep_outputs import (  # noqa: E402  (sibling script, one definition for both comparators)
    BND_OWN_ROWS_LOST_BUCKET,
    BND_SYNONYM_MATE_ROWS_LOST_BUCKET,
    DEFAULT_REFERENCE_RELEASE,
    FASTA_NAMED_SLICE_INTERGENIC_BUCKET,
    REFERENCE_116_2_ONE_SIDED_RULES,
    REFERENCE_RELEASES,
    REFERENCE_SKIPPED_RECORD_BUCKET,
    _MAX_SV_SIZE,
    _batch_dependent_annotation_scope,
    _location_span,
    _span_exceeds_max_sv_size,
    is_own_breakend_allele,
    is_structural_allele,
)


# Variant type classification

_ACGT = set("ACGT")


def vcf_basename(vcf_path: Path) -> str:
    """Extract the bare name from a VCF path, stripping .vcf and .gz suffixes."""
    name = vcf_path.stem  # "01_snv.vcf" for .vcf.gz, "01_snv" for .vcf
    name = name.removesuffix(".vcf")
    return name


def normalize_location(loc: str) -> str:
    """Strip chr prefix from Location for comparison (Rust normalizes, Perl preserves)."""
    parts = loc.split(":", 1)
    if len(parts) < 2:
        return loc
    chrom = parts[0]
    rest = parts[1]
    # Strip chr prefix (case-insensitive)
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    # Ensembl convention: M -> MT
    if chrom.upper() == "M":
        chrom = "MT"
    return f"{chrom}:{rest}"


def _strip_chr_in_bracket(allele: str) -> str:
    """Strip chr prefix inside BND bracket notation for allele comparison.

    Perl preserves the original VCF chromosome in BND alleles (e.g., N[chr21:1234[)
    while Rust normalizes to Ensembl style (e.g., N[21:1234[). Normalize both to
    the chr-stripped form so these match in concordance tuples.
    """

    return re.sub(r"(\[|\])chr", r"\g<1>", allele)


# fastVEP symbolic-allele normalisation (--scored-engine fastvep only)

# The `--scored-engine` value that turns the symbolic-allele mapping on. Compared
# case-insensitively; every other engine name (the default `vep-rs`, `perl`) leaves
# the Allele column as written.
FASTVEP_ENGINE = "fastvep"

# Mobile-element subtypes Perl VEP names in the Allele column, keyed by the
# upper-cased third segment of `<INS:ME:...>` / `<DEL:ME:...>`. Mirrors
# `me_display_allele` in crates/vep-core/src/variant.rs, which mirrors the
# `@mobile_elements` list and the `L1 -> LINE1` alias in Perl's
# `Parser.pm::get_SO_term`. Any other subtype falls back to the generic
# `mobile_element_<operation>`, as both engines do.
_ME_SUBTYPE_NAMES = {
    "ALU": "Alu",
    "L1": "LINE1",
    "LINE1": "LINE1",
    "SVA": "SVA",
    "HERV": "HERV",
}


def _is_fastvep(engine: str) -> bool:
    return engine.strip().lower() == FASTVEP_ENGINE


def fastvep_symbolic_allele_to_so_class(allele: str) -> str:
    """Map a raw symbolic VCF ALT token to the Allele string Perl VEP writes for it.

    Perl VEP writes an SV's Allele column as its Sequence Ontology class term
    (`OutputFactory.pm`: ``$vfoa->{breakend}->{string} || $svf->class_SO_term``), and
    vep-rs reproduces that through `VariantClass::so_term` plus `display_allele` in
    crates/vep-core/src/variant.rs. fastVEP instead writes the raw VCF
    token, so without this mapping every symbolic-allele tuple on the fastVEP arm is
    a guaranteed miss whatever its consequence set: ``<DEL>`` can never equal
    ``deletion``.

    The rules below are vep-rs's own classification (`classify_symbolic_alt` in
    crates/vep-cli/src/vcf_parser.rs, checked in this order) followed by its display
    rule, so the mapping is identical to what the scored vep-rs arm already writes:

        token (case-insensitive)   vep-rs VariantClass       Allele written
        -------------------------  ------------------------  --------------------------
        <INS:ME>                   MobileElementInsertion    mobile_element_insertion
        <INS:ME:ALU>               MobileElementInsertion    Alu_insertion
        <INS:ME:L1>, <INS:ME:LINE1>  MobileElementInsertion  LINE1_insertion
        <INS:ME:SVA>               MobileElementInsertion    SVA_insertion
        <INS:ME:HERV>              MobileElementInsertion    HERV_insertion
        <INS:ME:other>             MobileElementInsertion    mobile_element_insertion
        <DEL:ME[:x]>               MobileElementDeletion     as above, with _deletion
        <DEL...>                   StructuralDeletion        deletion
        <INS...>                   StructuralInsertion       insertion
        <DUP:TANDEM...>            TandemDuplication         tandem_duplication
        <DUP...>                   Duplication               duplication
        <INV...>                   Inversion                 inversion
        <CNV:TR...>                TandemRepeat              tandem_repeat
        <CN0> / <CN=0>             StructuralDeletion        deletion
        <CN2> / <CN=2>             Duplication               duplication
        <CNn> (any other n)        CopyNumberVariation       copy_number_variation
        <CNV...>, <CN...>          CopyNumberVariation       copy_number_variation
        <CPX[:x]>                  ComplexStructural         CPX

    Left UNCHANGED, deliberately:

    * ``<NON_REF>`` and ``<*>``: both engines write the literal token (Perl has no SO
      term for either; vep-rs `display_allele` returns the raw allele), so identity
      is already the match.
    * ``<BND>``: Perl synthesises bracket-notation and single-breakend alleles from
      ``INFO/CHR2``, ``END2``/``END`` and ``SVLEN``, none of which a per-row output
      token carries, so no token-level mapping can reproduce them.
    * Any other ``<...>`` token: vep-rs classifies it by ``INFO/SVTYPE`` when present
      and only falls back to the stripped token, and the comparator cannot see
      INFO. A miss stays a miss rather than inventing a Perl string.
    * Everything that is not a ``<...>`` token: literal bases, breakend notation,
      ``-``, and the SO class strings themselves are returned as-is, so applying the
      mapping to an already-normalised value is a no-op.

    Two Perl behaviours are out of reach of a per-row mapping and score as
    divergences on the fastVEP arm: a multi-allelic ``<CN0>,<CN2>`` record, which
    Perl coalesces into ONE ``copy_number_variation`` feature, and ``<CNV:TR>`` run
    with ``--fasta``, where Perl materialises the literal repeat sequence (the
    class `filter_cnv_tr_expansion_divergences` masks for vep-rs).
    """
    if len(allele) < 3 or not (allele.startswith("<") and allele.endswith(">")):
        return allele
    upper = allele.upper()
    inner = upper[1:-1]

    if upper in ("<*>", "<NON_REF>", "<BND>"):
        return allele

    # Mobile elements before generic DEL/INS: `<DEL:ME>` starts with `<DEL`.
    if upper.startswith("<INS:ME") or upper.startswith("<DEL:ME"):
        parts = inner.split(":")
        operation = "deletion" if parts[0] == "DEL" else "insertion"
        if len(parts) >= 3 and parts[2] in _ME_SUBTYPE_NAMES:
            return f"{_ME_SUBTYPE_NAMES[parts[2]]}_{operation}"
        return f"mobile_element_{operation}"
    if upper.startswith("<DEL"):
        return "deletion"
    if upper.startswith("<INS"):
        return "insertion"
    if upper.startswith("<DUP:TANDEM"):
        return "tandem_duplication"
    if upper.startswith("<DUP"):
        return "duplication"
    if upper.startswith("<INV"):
        return "inversion"
    if upper.startswith("<CNV:TR"):
        return "tandem_repeat"
    # `<CNn>` / `<CN=n>`: Perl maps CN0 -> deletion and CN2 -> duplication, every
    # other copy number to the generic class. An unparseable number falls through
    # to the generic rule below, as vep-rs's `parse::<u32>()` failure does.
    if upper.startswith("<CN") and not upper.startswith("<CNV"):
        num_str = inner[2:].removeprefix("=")
        if num_str.isascii() and num_str.isdigit():
            cn = int(num_str)
            if cn == 0:
                return "deletion"
            if cn == 2:
                return "duplication"
            return "copy_number_variation"
    if upper.startswith("<CNV") or upper.startswith("<CN"):
        return "copy_number_variation"
    if upper.startswith("<CPX"):
        # Perl strips the angle brackets and any `:subtype`; vep-rs's
        # ComplexStructural display rule does the same.
        return "CPX"
    return allele


def normalize_allele(allele: str, engine: str = "vep-rs") -> str:
    """Normalize allele string for concordance comparison.

    `engine` names the engine that wrote the row. Breakend bracket normalisation
    applies to every engine; the symbolic-token mapping applies ONLY when `engine`
    is fastVEP (`FASTVEP_ENGINE`), because vep-rs and Perl already write the SO
    class strings and mapping them again would be a no-op at best and a masking
    hazard at worst.
    """
    allele = _strip_chr_in_bracket(allele)
    if _is_fastvep(engine):
        allele = fastvep_symbolic_allele_to_so_class(allele)
    return allele


def classify_variant(ref: str, alt: str) -> str:
    """Classify a VCF variant by REF/ALT into one of the defined types."""
    # Multi-allelic (comma in ALT) -- check first since other rules assume
    # a single ALT allele.
    if "," in alt:
        return "Multi_allelic"

    # Spanning deletion
    if alt == "*":
        return "Spanning"

    # Reference-only
    if alt == ".":
        return "RefOnly"

    # NON_REF / <*>
    if alt.startswith("<NON_REF") or alt == "<*>":
        return "NON_REF"

    # Symbolic alleles (angle-bracket notation)
    if alt.startswith("<"):
        upper = alt.upper()
        if upper.startswith("<DEL"):
            return "Symbolic_DEL"
        if upper.startswith("<INS"):
            # <INS:ME:*> handled below, but plain <INS> first
            if "ME" in upper:
                return "Mobile_Element"
            return "Symbolic_INS"
        if upper.startswith("<DUP"):
            return "Symbolic_DUP"
        if upper.startswith("<INV"):
            return "Symbolic_INV"
        if upper.startswith("<CNV") or upper.startswith("<CN"):
            return "Symbolic_CNV"
        if upper.startswith("<CPX"):
            return "Complex_SV"
        # A symbolic breakend (`<BND>` with INFO/CHR2 and END2) is a breakend record like
        # its bracket-notation form.
        if upper.startswith("<BND"):
            return "BND"
        return "Other"

    # Mobile element (ALT contains ME but not in angle brackets -- rare)
    if "ME" in alt.upper() and alt.startswith("<"):
        return "Mobile_Element"

    # Breakend notation
    if "[" in alt or "]" in alt:
        return "BND"
    # Single breakend: ALT is e.g. "A." or ".A"
    if alt.endswith(".") and len(alt) > 1 and alt[:-1].isalpha():
        return "BND"
    if alt.startswith(".") and len(alt) > 1 and alt[1:].isalpha():
        return "BND"

    ref_upper = ref.upper()
    alt_upper = alt.upper()
    ref_is_acgt = all(c in _ACGT for c in ref_upper)
    alt_is_acgt = all(c in _ACGT for c in alt_upper)

    if ref_is_acgt and alt_is_acgt:
        lr, la = len(ref_upper), len(alt_upper)
        if lr == 1 and la == 1:
            return "SNV"
        if lr == la and lr > 1:
            return "MNP"
        if lr > la:
            return "Small_DEL"
        if lr < la:
            return "Small_INS"
        # Same length > 1 already covered by MNP
        return "Complex"

    return "Other"


# VCF parsing


@dataclass
class VcfVariant:
    chrom: str
    pos: int
    vid: str
    ref: str
    alt: str
    variant_type: str
    location: str  # "chrom:pos" for matching to VEP output
    info: str = ""  # the INFO column as written, read by the skip-route cross-check


def parse_vcf(path: Path) -> list[VcfVariant]:
    """Parse a VCF file (plain or bgzipped) and return classified variants."""
    variants: list[VcfVariant] = []
    if path.suffix == ".gz":
        # Use subprocess zcat for bgzip compatibility (Python gzip can't read bgzip multi-block)
        import subprocess as _sp

        _proc = _sp.Popen(
            ["zcat", str(path)], stdout=_sp.PIPE, stderr=_sp.DEVNULL, text=True
        )
        fh = _proc.stdout
    else:
        fh = open(path, "rt")
    with fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos_str, vid, ref, alt = (
                parts[0],
                parts[1],
                parts[2],
                parts[3],
                parts[4],
            )
            pos = int(pos_str)
            vtype = classify_variant(ref, alt)
            location = normalize_location(f"{chrom}:{pos}")
            info = parts[7] if len(parts) > 7 else ""
            variants.append(VcfVariant(chrom, pos, vid, ref, alt, vtype, location, info))
    return variants


# VEP output parsing

# Comparison tuple: (Location, Allele, Feature, Feature_type, Consequence)
AnnotTuple = tuple[str, str, str, str, str]


def normalize_consequence_set(raw: str) -> str:
    """Sort + dedupe a comma-separated consequence set.

    Mirrors ``compare_vep_outputs.py::normalize_consequence_set`` so the SV arm and
    the SNP/indel arm define a tuple identically. Without it the two F1 columns are
    not the same metric: two engines emitting the same terms in a different order
    would count as concordant on one arm and discordant on the other.

    Both engines order terms identically by SO rank, so no discordant pair differs
    only by ordering and this moves no number; it exists because relying on that
    agreement is a latent dependency on undocumented output ordering.
    """
    terms = [t.strip() for t in raw.split(",") if t.strip()]
    return ",".join(sorted(set(terms)))


def parse_vep_output(
    path: Path, engine: str = "vep-rs"
) -> tuple[set[AnnotTuple], set[str]]:
    """Parse VEP tab-delimited output. Returns (set of tuples, set of Locations).

    `engine` names the engine that produced `path` and is handed to
    `normalize_allele`; only fastVEP (`FASTVEP_ENGINE`) changes anything, mapping its
    raw symbolic ALT tokens to the SO class strings the other two engines write. The
    Perl ground truth is always parsed with the default, so the mapping can never
    touch the reference side.

    A missing file yields empty sets, which callers treat as "engine produced no
    output for this VCF"; `compute_metrics` scores a both-sides-empty comparison
    as F1 = 0.0 rather than as perfect agreement on zero tuples.

    Every other defect raises. A malformed row (fewer than 7 columns) is a hard
    error, because the alternative -- skipping it -- turns a truncated file into a
    smaller-but-plausible tuple set, and a smaller set with the same intersection
    reads as BETTER concordance.
    """
    tuples: set[AnnotTuple] = set()
    locations: set[str] = set()

    if not path.exists():
        return tuples, locations

    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line or line.startswith("##") or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 7:
                # A short row is malformed input, not a row to skip: skipped
                # silently, a file truncated mid-line scores as high concordance
                # on fewer tuples.
                raise ValueError(
                    f"{path}: line {lineno} has {len(cols)} columns, expected at "
                    f"least 7. A truncated or malformed engine output must fail "
                    f"the comparison rather than silently reduce the tuple set."
                )
            # col0=Uploaded_variation, col1=Location, col2=Allele,
            # col3=Gene, col4=Feature, col5=Feature_type, col6=Consequence
            location = normalize_location(cols[1])
            allele = normalize_allele(cols[2], engine)
            feature = cols[4]
            feature_type = cols[5]
            consequence = normalize_consequence_set(cols[6])
            locations.add(location)
            tuples.add((location, allele, feature, feature_type, consequence))

    return tuples, locations


def location_types_from_outputs(
    input_variants: list[VcfVariant], *output_paths: Path
) -> dict[str, str]:
    """Map every Location the engines print to the input record's variant type.

    A record's ``chrom:POS`` matches only a point Location: VEP prints a symbolic or
    breakend record at ``chrom:(POS+1)-END`` (``chrom:POS+1`` for a point breakend), so
    a ranged Location keyed by position alone reaches no type row and its discordant
    rows read ``Unknown``. The join is through the ``Uploaded_variation`` column, the
    record's ID, which both engines print beside the Location; a record without an ID
    (``.``) keeps the ``chrom:POS`` key. A multi-allelic record's alleles share its ID
    and its type.
    """
    type_by_id: dict[str, str] = {}
    types: dict[str, str] = {}
    for v in input_variants:
        types[v.location] = v.variant_type
        if v.vid and v.vid != ".":
            type_by_id.setdefault(v.vid, v.variant_type)
    for path in output_paths:
        if not path.exists():
            continue
        with open(path) as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t", 3)
                if len(cols) < 3:
                    continue
                vtype = type_by_id.get(cols[0])
                if vtype is not None:
                    types.setdefault(normalize_location(cols[1]), vtype)
    return types


# Metrics


@dataclass
class Metrics:
    total_input: int = 0
    perl_variants: int = 0
    rust_variants: int = 0
    perl_tuples: int = 0
    rust_tuples: int = 0
    intersection: int = 0
    only_perl: int = 0
    only_rust: int = 0
    precision: float = 0.0
    recall: float = 0.0
    f1: float = 0.0
    # Aggregate rows only: the union-based counts, kept beside the file-summed
    # ones so a reader can see how many tuples appear in more than one input VCF
    # instead of having the difference silently absorbed by a set union.
    perl_tuples_union: int | None = None
    rust_tuples_union: int | None = None
    intersection_union: int | None = None
    cross_file_duplicate_perl: int | None = None
    cross_file_duplicate_rust: int | None = None

    def recompute_prf(self) -> None:
        """Recompute precision/recall/F1 from the current count fields."""
        self.precision = (
            self.intersection / self.rust_tuples if self.rust_tuples else 0.0
        )
        self.recall = self.intersection / self.perl_tuples if self.perl_tuples else 0.0
        self.f1 = (
            2 * self.precision * self.recall / (self.precision + self.recall)
            if (self.precision + self.recall)
            else 0.0
        )


def compute_metrics(
    input_variants: list[VcfVariant],
    perl_tuples: set[AnnotTuple],
    perl_locations: set[str],
    rust_tuples: set[AnnotTuple],
    rust_locations: set[str],
) -> Metrics:
    m = Metrics()
    m.total_input = len(input_variants)

    # Count input variants that have at least one annotation line
    input_locs = {v.location for v in input_variants}
    m.perl_variants = len(input_locs & perl_locations)
    m.rust_variants = len(input_locs & rust_locations)

    m.perl_tuples = len(perl_tuples)
    m.rust_tuples = len(rust_tuples)
    m.intersection = len(perl_tuples & rust_tuples)
    m.only_perl = len(perl_tuples - rust_tuples)
    m.only_rust = len(rust_tuples - perl_tuples)

    if m.rust_tuples > 0:
        m.precision = m.intersection / m.rust_tuples
    if m.perl_tuples > 0:
        m.recall = m.intersection / m.perl_tuples
    if m.precision + m.recall > 0:
        m.f1 = 2 * m.precision * m.recall / (m.precision + m.recall)

    return m


# Discordant tuple collection


@dataclass
class DiscordantRecord:
    source_file: str
    variant_type: str
    direction: str  # "only_perl" or "only_rust"
    location: str
    allele: str
    feature: str
    feature_type: str
    consequence: str


def collect_discordants(
    source_file: str,
    input_variants: list[VcfVariant],
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    location_types: dict[str, str] | None = None,
) -> list[DiscordantRecord]:
    """One row per one-sided tuple, typed through `location_types` (the printed
    Locations, see `location_types_from_outputs`) when given, else by ``chrom:POS``."""
    records: list[DiscordantRecord] = []
    loc_to_type: dict[str, str] = dict(location_types or {})
    for v in input_variants:
        loc_to_type.setdefault(v.location, v.variant_type)

    for t in sorted(perl_tuples - rust_tuples):
        loc, allele, feature, ft, csq = t
        vtype = loc_to_type.get(loc, "Unknown")
        records.append(
            DiscordantRecord(
                source_file, vtype, "only_perl", loc, allele, feature, ft, csq
            )
        )

    for t in sorted(rust_tuples - perl_tuples):
        loc, allele, feature, ft, csq = t
        vtype = loc_to_type.get(loc, "Unknown")
        records.append(
            DiscordantRecord(
                source_file, vtype, "only_rust", loc, allele, feature, ft, csq
            )
        )

    return records


# Per-variant-type metrics


def compute_per_type_metrics(
    input_variants: list[VcfVariant],
    perl_tuples: set[AnnotTuple],
    perl_locations: set[str],
    rust_tuples: set[AnnotTuple],
    rust_locations: set[str],
    location_types: dict[str, str] | None = None,
) -> dict[str, Metrics]:
    """Group input variants by type and compute metrics for each group.

    A tuple belongs to a type through its Location: through `location_types` (the
    printed Locations mapped by record ID, see `location_types_from_outputs`) when
    given, else through the record's ``chrom:POS`` alone, which only a point Location
    matches.
    """
    # Group variants by type
    type_variants: dict[str, list[VcfVariant]] = defaultdict(list)
    for v in input_variants:
        type_variants[v.variant_type].append(v)
    printed_locs: dict[str, set[str]] = defaultdict(set)
    for loc, vtype in (location_types or {}).items():
        printed_locs[vtype].add(loc)

    type_metrics: dict[str, Metrics] = {}
    for vtype, variants in sorted(type_variants.items()):
        # Filter tuples to only those whose Location matches a variant of this type
        type_locs = {v.location for v in variants} | printed_locs.get(vtype, set())

        perl_filtered = {t for t in perl_tuples if t[0] in type_locs}
        rust_filtered = {t for t in rust_tuples if t[0] in type_locs}
        perl_locs_filtered = perl_locations & type_locs
        rust_locs_filtered = rust_locations & type_locs

        type_metrics[vtype] = compute_metrics(
            variants,
            perl_filtered,
            perl_locs_filtered,
            rust_filtered,
            rust_locs_filtered,
        )

    return type_metrics


# Divergence masks


# Transcript span type: (chromosome, start, end) as the cache files it.
TranscriptSpan = tuple[str, int, int]


def load_cache_authority(
    cache_dir: str | Path,
) -> tuple[
    dict[str, set[str]], dict[tuple[str, str], TranscriptSpan], dict[tuple[str, str], TranscriptModel]
]:
    """Read the vep-rs JSON cache once and return three views of it.

    The first maps each chromosome to the transcript ids the cache holds for it; the
    second maps ``(chromosome, transcript id)`` to the transcript's ``(chromosome, start,
    end)``; the third maps the same key to the transcript's structure (`TranscriptModel`:
    exons, strand, biotype, genomic CDS bounds, mature miRNA ranges). The chromosome is
    part of the key because the pseudoautosomal transcripts are filed under both X and Y
    with their own coordinates on each. Both engines are
    run against this one cache, so a transcript it files under chr21 is a chr21
    transcript and its structure is the structure both engines annotated against,
    whatever either engine's output says.

    Returns three empty mappings when the directory is absent or unreadable, and the
    callers treat that as "cannot judge" rather than "nothing to exclude".
    """
    by_chr: dict[str, set[str]] = {}
    spans: dict[tuple[str, str], TranscriptSpan] = {}
    models: dict[tuple[str, str], TranscriptModel] = {}
    for chrom, tx in iter_cache_transcripts(cache_dir):
        by_chr.setdefault(chrom, set()).add(tx["stable_id"])
        try:
            spans[(chrom, tx["stable_id"])] = (chrom, int(tx["start"]), int(tx["end"]))
        except (KeyError, TypeError, ValueError):
            pass
        model = model_from_cache_record(chrom, tx)
        if model is not None:
            models[(chrom, tx["stable_id"])] = model
    return by_chr, spans, models


def load_transcript_index(
    cache_dir: str | Path,
) -> tuple[dict[str, set[str]], dict[tuple[str, str], TranscriptSpan]]:
    """The chromosome and span views of the cache; see `load_cache_authority`."""
    by_chr, spans, _models = load_cache_authority(cache_dir)
    return by_chr, spans


def load_transcripts_by_chromosome(cache_dir: str | Path) -> dict[str, set[str]]:
    """Map each chromosome to the transcript ids the vep-rs JSON cache holds for it.

    The authority for "which chromosome is this transcript on"; see
    `load_transcript_index`, whose first element this is.
    """
    return load_transcript_index(cache_dir)[0]


def _chromosome_of_location(loc: str) -> str | None:
    """The chromosome a VEP Location names, normalised to the cache's convention."""
    if ":" not in loc:
        return None
    chrom = loc.split(":", 1)[0]
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def _allele_names_another_chromosome(allele: str) -> bool:
    """A bracket-notation breakend allele carries its mate's coordinate.

    `N[chr22:44739371[` legitimately reaches chr22, so a transcript there is not an
    error. Only non-breakend alleles are candidates for the cross-chromosome filter.
    """
    return "[" in allele or "]" in allele


def filter_cross_chromosome_divergences(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    transcripts_by_chr: dict[str, set[str]],
    scored_engine: str = "vep-rs",
) -> tuple[set[AnnotTuple], set[AnnotTuple], int]:
    """Exclude Perl tuples naming a transcript that is not on the variant's chromosome.

    THE DEFECT. Given a structural variant whose `CHROM` is `21` and whose `INFO/CHR2`
    reads `chr21`, Ensembl VEP r115.2 annotates it against transcripts on chromosome 22
    at numerically similar coordinates: for a 25 kb `<DEL>` at `21:23699593-23724683`,
    VEP names 122 transcripts of which 120 resolve to chr22 by Ensembl's own id lookup
    and only 2 to chr21. A deletion on one chromosome cannot truncate a transcript on
    another, so those 120 are wrong output.

    IT IS NOT ASSEMBLY-SPECIFIC. The GRCh38 gnomAD SV VCF carries `CHR2=chr21` against
    an unprefixed `CHROM` while the GRCh37 file carries `CHR2=21`, and the defect fires
    on both when Perl runs against a full-genome cache; a chr21-only Perl cache holds no
    off-chromosome transcript for the defect to resolve to, so it cannot fire there
    however bad the lookup is. Resolving the GRCh37 Perl-only transcript ids against
    Ensembl's own GRCh37 GFF3 places them on chrX, chrY and chr22, the chromosomes
    following 21 in Ensembl's ordering.

    THE MECHANISM. A neighbouring inter-chromosomal breakend in the same batch loads the
    cache regions around its mate on the other chromosome (ensembl-vep
    AnnotationSource.pm:249-257); the batch's interval tree holds every record by numeric
    coordinate with no chromosome test (InputBuffer.pm:345-366), so
    `get_overlapping_vfs` (InputBuffer.pm:284-330) hands the other chromosome's
    transcripts every record whose numbers overlap them; and a non-breakend record gets
    its overlap allele without the breakend distance gate that compares chromosome
    names (ensembl-variation StructuralVariationOverlap.pm:64-72 against :73-87). Alone
    in a batch, the same record writes `intergenic_variant`. Contig naming plays no
    part, which is why the GRCh37 file's `CHR2=21` fires it as readily as GRCh38's.

    Why this is masked rather than charged to vep-rs: vep-rs names exactly the transcripts
    its per-chromosome index holds for the variant's chromosome, which is correct, and
    raw F1 otherwise charges it for VEP's error. Unlike every other mask in this file the
    criterion is objectively checkable rather than mechanistic: a transcript either is or
    is not on the variant's chromosome.

    ONE-SIDED BY NATURE, and asserted to be so WHEN THE SCORED ENGINE IS vep-rs. vep-rs
    cannot produce this shape, because it queries a per-chromosome index built from the
    very cache passed here; if it ever does, that is a vep-rs bug and masking it would
    hide the bug. So a vep-rs-side match raises rather than excludes.

    `scored_engine` exists because that assertion is about vep-rs specifically, not about
    the mask. A third-party engine builds its own transcript set from its own annotation
    release, so it legitimately names transcripts this cache does not hold for the
    chromosome, and raising there would score nothing at all. For any engine other than
    vep-rs the shape is COUNTED and returned as the third element rather than raised, and
    those tuples stay in both raw and adjusted F1: they belong to the scored engine, and no
    Perl defect excuses them.

    Bracket-notation breakend alleles are skipped: their mate coordinate legitimately names
    another chromosome.
    """
    excluded_perl: set[AnnotTuple] = set()
    unresolvable_scored = 0
    if not transcripts_by_chr:
        return set(), excluded_perl, 0

    for source, tuples in (("rust", rust_tuples - perl_tuples), ("perl", perl_tuples - rust_tuples)):
        for tup in tuples:
            loc, allele, feature, _ftype, _csq = tup
            if not feature.startswith("ENST") or _allele_names_another_chromosome(allele):
                continue
            chrom = _chromosome_of_location(loc)
            local = transcripts_by_chr.get(chrom or "")
            # No index for that chromosome means the cache cannot adjudicate it.
            if not local or feature in local:
                continue
            if source == "perl":
                excluded_perl.add(tup)
            elif scored_engine == "vep-rs":
                raise AssertionError(
                    f"vep-rs named {feature}, which the cache does not hold for "
                    f"chromosome {chrom} ({loc}). vep-rs queries a per-chromosome index, "
                    f"so this shape should be unreachable; it is a vep-rs bug and must not "
                    f"be masked as a VEP defect."
                )
            else:
                unresolvable_scored += 1

    return set(), excluded_perl, unresolvable_scored


_INTERGENIC_TUPLE_TAIL = ("-", "-", "intergenic_variant")


def filter_cross_chromosome_orphan_intergenic(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    excluded_perl_cross_chromosome: set[AnnotTuple],
) -> set[AnnotTuple]:
    """The vep-rs side of the cross-chromosome mask: the ``intergenic_variant`` row
    that stands where Perl's off-chromosome transcripts stood.

    THE SHAPE. For a non-breakend structural variant Perl annotated ONLY against
    transcripts on other chromosomes (every one of its tuples for the record is in
    `excluded_perl_cross_chromosome`), and vep-rs, finding no transcript within 5 kb
    on the record's own chromosome, wrote the single row
    ``(location, allele, -, -, intergenic_variant)``. After the Perl-side mask the
    record has no Perl tuple left, so the vep-rs row is charged as a one-sided
    divergence although it is the row Perl itself writes for the record when no
    other record loads an off-chromosome region: ensembl-variation
    StructuralVariationFeature.pm:674-689 creates the intergenic overlap exactly when
    the transcript overlap list is empty, and ensembl-vep AnnotationSource.pm:249-257
    is where a neighbouring breakend's mate loads the other chromosome's regions that
    InputBuffer.pm:284-330 then matches against this record by numeric position alone.

    THE PREDICATE, on a vep-rs-only tuple: its feature is ``-`` and its consequence
    ``intergenic_variant``; Perl has at least one tuple at the same Location and Allele;
    every Perl tuple at that key is in `excluded_perl_cross_chromosome`; and vep-rs's
    whole set at that key is this one tuple. A record with any surviving Perl tuple, a
    record Perl wrote nothing for, or a vep-rs set with anything beside the intergenic
    row is never touched: those stay charged where they fall.

    Returns the vep-rs tuples to exclude. Perl-side only by construction of the inputs,
    since `excluded_perl_cross_chromosome` is what `filter_cross_chromosome_divergences`
    returned; with no cache that set is empty and this excludes nothing.
    """
    if not excluded_perl_cross_chromosome:
        return set()
    perl_by_variant: dict[tuple[str, str], set[AnnotTuple]] = defaultdict(set)
    for t in perl_tuples:
        perl_by_variant[(t[0], t[1])].add(t)
    rust_by_variant: dict[tuple[str, str], set[AnnotTuple]] = defaultdict(set)
    for t in rust_tuples:
        rust_by_variant[(t[0], t[1])].add(t)
    excluded: set[AnnotTuple] = set()
    for t in rust_tuples - perl_tuples:
        if t[2:] != _INTERGENIC_TUPLE_TAIL:
            continue
        key = (t[0], t[1])
        perl_here = perl_by_variant.get(key)
        if not perl_here or not perl_here <= excluded_perl_cross_chromosome:
            continue
        if rust_by_variant[key] != {t}:
            continue
        excluded.add(t)
    return excluded




def filter_intended_divergences(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
) -> tuple[set[AnnotTuple], set[AnnotTuple]]:
    """Identify transcript-selection divergences in BOTH directions (bidirectional filter).

    The class is ``sv_nondeterministic_transcript_selection``: for a variant above
    Perl's --max_sv_size, the transcripts Perl annotates depend on which cache regions
    its input buffer happened to load, while vep-rs annotates every overlapping
    transcript deterministically. The divergence runs in both directions:
      - Rust-only tuples: transcripts vep-rs annotated that Perl's buffer never loaded
      - Perl-only tuples: transcripts Perl annotated that vep-rs did not name
        (because Perl's buffer loaded a different transcript set for the variant)

    Returns (excluded_rust_only, excluded_perl_only), sets to exclude from adjusted
    metrics. Only excludes tuples where:
      - The tuple is in one output but not the other
      - The OTHER engine DID annotate the same variant (same Location + Allele)
      - But the other engine used a DIFFERENT set of transcripts (Feature not in
        the other engine's set for that variant)
      - The direction and span guards below both admit it
    """
    rust_only = rust_tuples - perl_tuples
    perl_only = perl_tuples - rust_tuples

    # Build: (full_location, allele) -> set of features seen in each engine.
    #
    # Keyed on the FULL location. The start position alone is a position rather
    # than a variant identity, and two records sharing a start with different ends
    # would collapse together. `filter_cnv_tr_expansion_divergences` keys on the end
    # coordinate for the mirror-image reason.
    perl_features_by_variant: dict[tuple[str, str], set[str]] = defaultdict(set)
    for loc, allele, feature, _ft, _csq in perl_tuples:
        perl_features_by_variant[(loc, allele)].add(feature)

    rust_features_by_variant: dict[tuple[str, str], set[str]] = defaultdict(set)
    for loc, allele, feature, _ft, _csq in rust_tuples:
        rust_features_by_variant[(loc, allele)].add(feature)

    # vep-rs-only: exclude only where the defect can reach.
    #
    # Two guards, both load-bearing:
    #
    # (a) DIRECTION. The masked defect is Perl seeing FEWER transcripts than
    #     vep-rs, because a variant above --max_sv_size is skipped during
    #     cache-region loading and annotated against whatever its batch happened
    #     to load. Without a direction test this arm would exclude every vep-rs-only
    #     tuple on any variant Perl also touched, which is exactly the shape of a
    #     vep-rs FALSE POSITIVE: a transcript vep-rs named and Perl did not. The
    #     Perl-side arm below carries its mirror of this guard.
    #
    # (b) SCOPE. The mechanism requires `vep_skip`: a span above Perl's
    #     --max_sv_size default (10 Mb, ensembl-vep Config.pm:310) or an ALT type
    #     Perl has no Sequence Ontology term for (`_VEP_SKIP_UNSUPPORTED_ALLELES`).
    #     Otherwise Perl requests every region a variant overlaps, so no batch
    #     dependence exists and nothing is masked: a spurious transcript on a
    #     100 bp deletion or a 1 bp point locus stays charged to vep-rs.
    #
    # Known blind spot, deliberately left: a giant INTER-CHROMOSOMAL breakend
    # carries its mate coordinate in the Allele bracket, not the Location, so its
    # location span reads 0 and this gate does not admit it. That is the
    # conservative direction (it masks less), and the residual stays counted
    # against vep-rs.
    excluded_rust: set[AnnotTuple] = set()
    for t in rust_only:
        loc, allele, feature, _ft, _csq = t
        perl_feats = perl_features_by_variant.get((loc, allele))
        if perl_feats is None or feature in perl_feats:
            continue
        rust_feats = rust_features_by_variant.get((loc, allele), set())
        if len(rust_feats) <= len(perl_feats):
            # vep-rs did not name MORE transcripts here, so Perl under-annotation
            # is not the explanation. Leave it in both raw and adjusted F1.
            continue
        if not _batch_dependent_annotation_scope(loc, allele):
            # Below --max_sv_size and with a supported ALT type Perl loads every
            # region it needs, so a transcript-set difference here is not the
            # masked defect.
            continue
        excluded_rust.add(t)

    # Perl-only: exclude if vep-rs annotated the same variant but with different
    # transcripts AND vep-rs named at least as many transcripts as Perl did for
    # that variant.
    #
    # The direction guard is load-bearing. The defect this filter masks is Perl's
    # buffer-dependent cache loading, in which Perl sees only the transcripts a
    # neighbouring variant happened to load and therefore names FEWER than vep-rs,
    # which always annotates the complete set for the chromosome. Without a
    # direction test the same predicate also absorbs the opposite shape, vep-rs
    # naming FEWER transcripts than Perl, which is not a Perl defect (on the 25 kb
    # `<DEL>` at 21:23699593-23724683 in gnomAD SV, Perl names 122 transcripts and
    # vep-rs names 2, an under-annotation this filter must not mask); masking it
    # would convert a vep-rs coverage gap into an excluded transcript-selection
    # divergence and inflate adjusted F1. This arm carries the same SCOPE gate as its sibling
    # above and keys on the full location for the same identity reason.
    excluded_perl: set[AnnotTuple] = set()
    for t in perl_only:
        loc, allele, feature, _ft, _csq = t
        rust_feats = rust_features_by_variant.get((loc, allele))
        if rust_feats is None or feature in rust_feats:
            continue
        perl_feats = perl_features_by_variant.get((loc, allele), set())
        if len(rust_feats) < len(perl_feats):
            # vep-rs under-annotated this variant relative to Perl: the opposite
            # of the masked defect. Leave it in both raw and adjusted F1.
            continue
        if not _batch_dependent_annotation_scope(loc, allele):
            continue
        excluded_perl.add(t)

    return excluded_rust, excluded_perl


# Perl's breakend-to-feature distance gate (ensembl-variation StructuralVariationOverlap.pm
# :130-146, MAX_DISTANCE_FROM_TRANSCRIPT in Utils/VariationEffect.pm:60).
_BREAKEND_FEATURE_DISTANCE = 5000

_BRACKET_MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def _bracket_mate(allele: str) -> tuple[str, int] | None:
    """The mate ``(chromosome, position)`` a bracket-notation breakend allele names."""
    m = _BRACKET_MATE_RE.search(allele)
    if not m:
        return None
    chrom = m.group(1)
    return (chrom[3:] if chrom.lower().startswith("chr") else chrom), int(m.group(2))


def filter_giant_breakend_mate_divergences(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    transcript_spans: dict[tuple[str, str], TranscriptSpan],
) -> set[AnnotTuple]:
    """The mate-allele arm of the transcript-selection mask.

    THE SHAPE. A breakend record whose Location span exceeds --max_sv_size carries two
    alleles in the output: the local breakend (``N.``) and the mate (``N[chr:pos[``).
    Perl flags the record `vep_skip` (Parser.pm:493-497) and loads no region for it,
    neither for the span nor for the mate, because AnnotationSource.pm:238 skips the
    record before the breakend loop at :249-257. Both alleles are then annotated only
    against transcripts other records of the batch loaded. `filter_intended_divergences`
    reaches the local allele, where Perl names some transcripts and vep-rs more; it
    cannot reach the mate allele when Perl names none under it, because its arms key
    on the (Location, Allele) pair and require the other engine to have annotated it.
    Perl's own mate allele for such a record appears the moment another record of the
    batch loads the mate's region, and it then carries exactly the transcripts within
    `_BREAKEND_FEATURE_DISTANCE` of the mate (StructuralVariationOverlap.pm:73-87).

    THE PREDICATE, on a vep-rs-only tuple: its Allele is bracket notation whose mate
    the cache places on the same chromosome as its feature; the Location span exceeds
    --max_sv_size; Perl has at least one tuple at the Location (it saw the record) and
    none at the (Location, Allele); and the mate lies within
    `_BREAKEND_FEATURE_DISTANCE` of the feature's cached span, so Perl's own gate would
    have created the allele. A mate farther from the feature than that is a vep-rs
    tuple Perl never writes and stays charged to vep-rs.

    Returns the vep-rs tuples to exclude. With no `transcript_spans` (no cache) it
    excludes nothing.
    """
    if not transcript_spans:
        return set()
    perl_alleles_by_location: dict[str, set[str]] = defaultdict(set)
    for loc, allele, _f, _ft, _c in perl_tuples:
        perl_alleles_by_location[loc].add(allele)
    excluded: set[AnnotTuple] = set()
    for t in rust_tuples - perl_tuples:
        loc, allele, feature, _ft, _c = t
        mate = _bracket_mate(allele)
        if mate is None or not _span_exceeds_max_sv_size(loc):
            continue
        perl_alleles = perl_alleles_by_location.get(loc)
        if not perl_alleles or allele in perl_alleles:
            continue
        span = transcript_spans.get((mate[0], feature))
        if span is None:
            continue
        if not (span[1] - _BREAKEND_FEATURE_DISTANCE <= mate[1] <= span[2] + _BREAKEND_FEATURE_DISTANCE):
            continue
        excluded.add(t)
    return excluded


def filter_registry_swap_pairs(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
) -> tuple[set[AnnotTuple], set[AnnotTuple], dict[str, int]]:
    """Apply the SNP/indel comparator's excluding rules to this comparator's pairs.

    A pair is a (Location, Allele, Feature, Feature_type) present on both sides with
    different consequence sets; each side contributes exactly one tuple per key. The
    rules are `EXCLUDING_RULES` of ``compare_vep_outputs.py``, imported so the two
    comparators cannot drift: the start co-emission (VariationEffect.pm:884-893 fires
    ``start_lost`` on any structural span over the start codon while :1037 makes
    ``start_retained_variant`` fire on every structural allele without reading
    sequence) and the covered splice-region shape (BaseTranscriptVariationAllele.pm:215).
    Structural alleles reach the sequence-variant defects through the same predicates,
    and a sequence allele in an SV corpus reaches them as it does in the SNP/indel one.

    Returns (excluded_rust, excluded_perl, matches_by_bucket): both members of every
    matching pair, and how many pairs each rule matched.
    """
    from compare_vep_outputs import EXCLUDING_RULES  # noqa: E402  (sibling script)

    perl_by_key: dict[tuple[str, str, str, str], list[AnnotTuple]] = defaultdict(list)
    rust_by_key: dict[tuple[str, str, str, str], list[AnnotTuple]] = defaultdict(list)
    for t in perl_tuples - rust_tuples:
        perl_by_key[t[:4]].append(t)
    for t in rust_tuples - perl_tuples:
        rust_by_key[t[:4]].append(t)
    excluded_rust: set[AnnotTuple] = set()
    excluded_perl: set[AnnotTuple] = set()
    by_bucket: dict[str, int] = {rule.bucket: 0 for rule in EXCLUDING_RULES}
    for key, perl_side in perl_by_key.items():
        rust_side = rust_by_key.get(key)
        if not rust_side or len(perl_side) != 1 or len(rust_side) != 1:
            continue
        (pt,), (rt,) = perl_side, rust_side
        for rule in EXCLUDING_RULES:
            if rule.matches(pt[4], rt[4], key[1]):
                excluded_perl.add(pt)
                excluded_rust.add(rt)
                by_bucket[rule.bucket] += 1
                break
    return excluded_rust, excluded_perl, by_bucket


_NON_REF_ALLELE = "<NON_REF>"


def filter_non_ref_batch_divergences(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
) -> tuple[set[AnnotTuple], set[AnnotTuple]]:
    """Perl's batch-derived transcript tuples on a gVCF reference block.

    THE DEFECT, Perl's own contradiction. Perl warns `NON_REF is not a supported
    structural variant type` and keeps the record with `vep_skip` (Parser/VCF.pm:477-481,
    :575), loads no cache region for it (AnnotationSource.pm:238) and bounds nothing by
    it (:143), yet annotates it against whatever regions its batch loaded
    (AnnotationType/Transcript.pm:108 over the batch's numeric interval tree,
    InputBuffer.pm:284-330). Alone in a batch, a `<NON_REF>` record writes exactly one
    row, `intergenic_variant`, on both assemblies (StructuralVariationFeature.pm:674-689:
    the intergenic overlap exists when the transcript overlap list is empty); in a batch
    with other records it carries transcript rows for whatever regions those records
    loaded, and the key set moves with `--fork`. vep-rs writes that one row for every
    `<NON_REF>` record, Perl's isolated
    output. Supporting, not load-bearing: a reference block states that the sample
    matches the reference, so no consequence is called for.

    `<*>`, the same gVCF construct written the other way, is outside this filter: it is
    not `vep_skip`, so both engines annotate it as a span and their rows agree.

    THE PREDICATE, on a (Location, Allele) key whose Allele is `<NON_REF>` and whose
    vep-rs set is exactly the intergenic row: every Perl Transcript tuple at the key is
    excluded, and the vep-rs intergenic row is excluded when Perl wrote Transcript tuples
    and no intergenic row there (Perl's own isolated output displaced by its batch). A
    record on which vep-rs wrote anything else is never touched, and a key where Perl
    also wrote the intergenic row keeps it on both sides, matched.

    Returns (excluded_rust, excluded_perl).
    """
    perl_by_variant: dict[tuple[str, str], set[AnnotTuple]] = defaultdict(set)
    for t in perl_tuples:
        if t[1] == _NON_REF_ALLELE:
            perl_by_variant[(t[0], t[1])].add(t)
    rust_by_variant: dict[tuple[str, str], set[AnnotTuple]] = defaultdict(set)
    for t in rust_tuples:
        if t[1] == _NON_REF_ALLELE:
            rust_by_variant[(t[0], t[1])].add(t)
    excluded_rust: set[AnnotTuple] = set()
    excluded_perl: set[AnnotTuple] = set()
    for key, rust_here in rust_by_variant.items():
        intergenic = key + _INTERGENIC_TUPLE_TAIL
        if rust_here != {intergenic}:
            continue
        perl_here = perl_by_variant.get(key, set())
        perl_transcripts = {t for t in perl_here if t[2] != "-"}
        if not perl_transcripts:
            continue
        excluded_perl |= perl_transcripts
        if intergenic not in perl_here:
            excluded_rust.add(intergenic)
    return excluded_rust, excluded_perl


def filter_breakend_mate_local_read_pairs(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    transcript_models: dict[tuple[str, str], TranscriptModel],
) -> tuple[set[AnnotTuple], set[AnnotTuple], dict[str, int]]:
    """A point breakend's mate-side Transcript row, read by Perl at the local coordinate.

    THE SHAPE. Perl builds a Transcript row for every transcript within 5 kb of a
    breakend's mate on the mate chromosome (StructuralVariationOverlap.pm:73-87) and then
    evaluates every positional predicate on it with the LOCAL variation feature
    (BaseVariationFeatureOverlapAllele.pm:257,273); only `feature_truncation` reads the
    mate (Utils/VariationEffect.pm:358). The row therefore carries the region term the
    local coordinate lands on when applied numerically to the mate transcript, or
    `intergenic_variant` on a Transcript row when nothing fires. Moving the local
    position of one record while its mate stays put moves the row's terms; on a
    chromosome-21 coordinate against a chromosome-5 transcript Perl computes a CDS
    position and a protein position. vep-rs evaluates the mate row at the mate
    coordinate, the position the row describes.

    THE PREDICATE, on a pair (one tuple per side at a Location, Allele, Feature,
    Feature_type): the row is a mate-side row of a point breakend
    (`mate_side_derivation`), Perl's set equals the local-coordinate derivation, vep-rs's
    set equals the mate-coordinate derivation, and the two derivations differ. Both
    members are excluded. A vep-rs set that differs from the mate-coordinate derivation
    in any term stays charged, as does a Perl set the local read does not explain.

    Returns (excluded_rust, excluded_perl, pairs_by_kind) with the pairs counted under
    ``cross_chromosome`` and ``same_chromosome``. With no `transcript_models` (no cache)
    it excludes nothing.
    """
    by_kind = {"cross_chromosome": 0, "same_chromosome": 0}
    if not transcript_models:
        return set(), set(), by_kind
    perl_by_key: dict[tuple[str, str, str, str], list[AnnotTuple]] = defaultdict(list)
    rust_by_key: dict[tuple[str, str, str, str], list[AnnotTuple]] = defaultdict(list)
    for t in perl_tuples - rust_tuples:
        perl_by_key[t[:4]].append(t)
    for t in rust_tuples - perl_tuples:
        rust_by_key[t[:4]].append(t)
    excluded_rust: set[AnnotTuple] = set()
    excluded_perl: set[AnnotTuple] = set()
    for key, perl_side in perl_by_key.items():
        rust_side = rust_by_key.get(key)
        if not rust_side or len(perl_side) != 1 or len(rust_side) != 1:
            continue
        loc, allele, feature, _ftype = key
        mate = _bracket_mate(allele)
        if mate is None:
            continue
        model = transcript_models.get((mate[0], feature))
        if model is None:
            continue
        derived = mate_side_derivation(loc, allele, model)
        if derived is None:
            continue
        mate_set, local_set = derived
        if mate_set == local_set:
            continue
        (pt,), (rt,) = perl_side, rust_side
        if set(pt[4].split(",")) != local_set or set(rt[4].split(",")) != mate_set:
            continue
        excluded_perl.add(pt)
        excluded_rust.add(rt)
        local_chrom = _chromosome_of_location(loc)
        by_kind["same_chromosome" if local_chrom == mate[0] else "cross_chromosome"] += 1
    return excluded_rust, excluded_perl, by_kind


# <CNV:TR> literal-versus-symbolic mask

# vep-rs writes this literal string in the Allele column for <CNV:TR> records.
_CNV_TR_RUST_ALLELE = "tandem_repeat"

# The direction term the symbolic reading adds, by the direction of the record's copy change.
_CNV_TR_DIRECTION_TERM = {"gain": "feature_elongation", "loss": "feature_truncation"}

# The coding effect the literal reading names where the symbolic reading names the region:
# an in-frame or frameshifting indel of whole repeat units inside the CDS reads as
# `coding_sequence_variant` plus the direction term on the symbolic side.
_CNV_TR_LITERAL_CODING_TERMS = {
    "gain": frozenset({"inframe_insertion", "frameshift_variant"}),
    "loss": frozenset({"inframe_deletion", "frameshift_variant"}),
}
_CNV_TR_SYMBOLIC_CODING_TERM = "coding_sequence_variant"


def _cnv_tr_pair_kind(perl_loc: str, perl_allele: str, rust_loc: str) -> str | None:
    """Whether a Perl tuple is the literal reading of the tandem repeat a vep-rs tuple
    reads symbolically, and which direction: ``"gain"``, ``"loss"`` or None.

    VEP expands ``RUS x RUC`` against the reference run ``SVLEN`` spans and trims the
    common prefix, so its tuple sits at the run's 3' end: a gain is an insertion of
    whole units between the run's last base and the next (Location ``end-(end+1)``,
    Allele the inserted bases), a loss a deletion of whole units ending at the run's
    last base (Location inside the run, Allele ``-``). vep-rs's tuple spans the run.
    """
    ps, rs = _location_span(perl_loc), _location_span(rust_loc)
    if ps is None or rs is None or ps[0] != rs[0]:
        return None
    _, run_start, run_end = rs
    _, p_start, p_end = ps
    if perl_allele == "-" and p_end == run_end and p_start >= run_start and p_start <= p_end:
        return "loss"
    if perl_allele and all(c in _ACGT for c in perl_allele) and p_start == run_end and p_end == run_end + 1:
        return "gain"
    return None


def _is_cnv_tr_representation_pair(perl_csq: str, rust_csq: str, kind: str) -> bool:
    """Whether a paired literal and symbolic reading of one tandem repeat differ only
    in representation.

    The two term sets must agree once the symbolic reading's direction term
    (`feature_elongation` on a gain, `feature_truncation` on a loss) is set aside, with
    one further correspondence inside the CDS: the literal reading names the coding
    effect of the unit indel (`inframe_insertion`, `inframe_deletion`,
    `frameshift_variant`) where the symbolic reading names the region
    (`coding_sequence_variant`). Any other difference, a start or splice term one side
    carries and the other lacks, is a consequence disagreement and stays charged.
    """
    pset = {t for t in perl_csq.split(",") if t}
    rset = {t for t in rust_csq.split(",") if t}
    perl_extra = pset - rset
    rust_extra = rset - pset - {_CNV_TR_DIRECTION_TERM[kind]}
    if rust_extra == {_CNV_TR_SYMBOLIC_CODING_TERM}:
        return bool(perl_extra) and perl_extra <= _CNV_TR_LITERAL_CODING_TERMS[kind]
    return not rust_extra and not perl_extra


def _cnv_tr_pairs(
    perl_only: set[AnnotTuple], rust_only: set[AnnotTuple]
) -> list[tuple[AnnotTuple, AnnotTuple, str]]:
    """Every (Perl tuple, vep-rs tuple, kind) pair of one tandem repeat on one feature:
    the vep-rs tuple carries the symbolic allele and the Perl tuple is its literal
    reading by :func:`_cnv_tr_pair_kind`."""
    rust_by_key: dict[tuple[str, str, str], list[AnnotTuple]] = defaultdict(list)
    for t in rust_only:
        loc, allele, feature, ftype, _csq = t
        if allele != _CNV_TR_RUST_ALLELE:
            continue
        span = _location_span(loc)
        if span is None:
            continue
        rust_by_key[(span[0], feature, ftype)].append(t)
    pairs: list[tuple[AnnotTuple, AnnotTuple, str]] = []
    if not rust_by_key:
        return pairs
    for pt in perl_only:
        ploc, pallele, pfeature, pftype, _pcsq = pt
        span = _location_span(ploc)
        if span is None:
            continue
        for rt in rust_by_key.get((span[0], pfeature, pftype), ()):
            kind = _cnv_tr_pair_kind(ploc, pallele, rt[0])
            if kind is not None:
                pairs.append((pt, rt, kind))
    return pairs


def filter_cnv_tr_expansion_divergences(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
) -> tuple[set[AnnotTuple], set[AnnotTuple]]:
    """Identify the <CNV:TR> literal-versus-symbolic pairs whose readings differ only
    in representation, on both sides.

    A pair is a Perl tuple and a vep-rs tuple of one tandem repeat on one feature
    (:func:`_cnv_tr_pairs`); it is excluded when
    :func:`_is_cnv_tr_representation_pair` holds for its consequence sets.

    Returns (excluded_rust_only, excluded_perl_only).
    """
    rust_only = rust_tuples - perl_tuples
    perl_only = perl_tuples - rust_tuples
    excluded_rust: set[AnnotTuple] = set()
    excluded_perl: set[AnnotTuple] = set()
    for pt, rt, kind in _cnv_tr_pairs(perl_only, rust_only):
        if _is_cnv_tr_representation_pair(pt[4], rt[4], kind):
            excluded_perl.add(pt)
            excluded_rust.add(rt)
    return excluded_rust, excluded_perl


def count_cnv_tr_swap_population(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
) -> tuple[int, int]:
    """Size of the <CNV:TR> literal-versus-symbolic class on each side, BEFORE the mask.

    Returns ``(rust_total, perl_total)``: the vep-rs-only tuples whose Allele is the
    symbolic ``tandem_repeat`` marker, and the Perl-only tuples that are the literal
    reading of one of them (:func:`_cnv_tr_pairs`). Those are the two populations
    :func:`filter_cnv_tr_expansion_divergences` draws its pairs from, so
    ``excluded_cnvtr_rust <= rust_total`` and ``excluded_cnvtr_perl <= perl_total`` by
    construction, and each ratio is the mask's coverage of the class on that side.

    This is the denominator behind the "N of M pairs" statement of the mask and behind
    ``manuscript/data/discordance_taxonomy.csv``'s ``n_tuples`` column. The report
    carries it as ``cnvtr_swap_pairs_total_rust`` / ``cnvtr_swap_pairs_total_perl``
    (plus a per-file split) beside the excluded counts, so the fraction is measured by
    the same run that measures its numerator rather than read off ``discordant.tsv``.

    A Perl tuple has no allele signature of its own on a loss (its Allele is ``-``), so
    the Perl side is counted through the pairing; a vep-rs ``tandem_repeat`` tuple with
    no literal partner still counts on its side.

    Not a ``filter_*_divergences`` function: it excludes nothing, and that prefix is
    reserved for the functions that remove tuples from a side.
    """
    rust_only = rust_tuples - perl_tuples
    perl_only = perl_tuples - rust_tuples
    rust_total = sum(1 for t in rust_only if t[1] == _CNV_TR_RUST_ALLELE)
    perl_total = len({pt for pt, _rt, _kind in _cnv_tr_pairs(perl_only, rust_only)})
    return rust_total, perl_total


# The reference release selector (Ensembl VEP 115.2 or 116.2)

# The masks `apply_masks` runs against each reference release.
MASKS_BY_RELEASE: dict[str, tuple[str, ...]] = {
    "115.2": (
        "filter_intended_divergences",
        "filter_cnv_tr_expansion_divergences",
        "filter_cross_chromosome_divergences",
        "filter_cross_chromosome_orphan_intergenic",
        "filter_giant_breakend_mate_divergences",
        "filter_registry_swap_pairs",
        "filter_non_ref_batch_divergences",
        "filter_breakend_mate_local_read_pairs",
    ),
    "116.2": (
        "filter_cnv_tr_expansion_divergences",
        "filter_registry_swap_pairs",
        "filter_breakend_mate_local_read_pairs",
        "filter_reference_skipped_records",
        "filter_bnd_own_rows_lost",
        "filter_bnd_synonym_mate_rows",
    ),
}
# Why each 115.2 mask is not applied against 116.2, for the report.
MASKS_RETIRED_IN_116_2: dict[str, str] = {
    "filter_intended_divergences": "116 drops a vep_skip record before the input buffer; the reference writes no row to compare transcript sets against",
    "filter_giant_breakend_mate_divergences": "116 drops a breakend pair over --max_sv_size whole; its mate allele has no reference row",
    "filter_non_ref_batch_divergences": "116 drops a <NON_REF> record before the input buffer",
    "filter_cross_chromosome_divergences": "116.1's chromosome-aware overlap never pairs a record with a transcript on another chromosome",
    "filter_cross_chromosome_orphan_intergenic": "with no off-chromosome rows to displace it, the reference writes the intergenic row itself",
}


def _tuples_by_location(tuples: set[AnnotTuple]) -> dict[str, set[AnnotTuple]]:
    by_loc: dict[str, set[AnnotTuple]] = defaultdict(set)
    for t in tuples:
        by_loc[t[0]].add(t)
    return by_loc


# The reference's own skip routes, ported for the cross-check below. `_SO_TERMS` is
# ensembl-variation Utils/Config.pm:1417-1438 at release 116; the abbreviation rules are
# Parser.pm 116.2:663-726 (`get_SO_term`), the type derivation Parser/VCF.pm:468-481, the
# routing of a record to the structural path Parser/VCF.pm:236-246, the coordinates
# ensembl-io BaseVCF4.pm `get_start` / `get_end` (the base before the event is the
# start; SVLEN sets the end before END does; a breakend's END is discarded).
_SO_TERMS = {
    "INS": "insertion", "INS_ME": "mobile_element_insertion", "INS_ALU": "Alu_insertion",
    "INS_HERV": "HERV_insertion", "INS_LINE1": "LINE1_insertion", "INS_SVA": "SVA_insertion",
    "DEL": "deletion", "DEL_ME": "mobile_element_deletion", "DEL_ALU": "Alu_deletion",
    "DEL_HERV": "HERV_deletion", "DEL_LINE1": "LINE1_deletion", "DEL_SVA": "SVA_deletion",
    "TREP": "tandem_repeat", "TDUP": "tandem_duplication", "DUP": "duplication",
    "CNV": "copy_number_variation", "INV": "inversion", "BND": "chromosome_breakpoint",
}
_SV_ROUTED = re.compile(r"[<\[\]][^*]+[>\]\[]|^\.\w+|\w+\.$")
_VCF44_SYMBOLIC = re.compile(r"^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)")


def perl_so_abbrev(sv_type: str) -> str | None:
    """The Sequence Ontology abbreviation `get_SO_term` derives for a type string, or None
    when the reference has no term for it (the type is then `vep_skip`)."""
    if re.search(r"[/,]", sv_type):
        parts = re.split(r"[/,]", sv_type)
        normalized = [re.sub(r"^<|>$", "", part.strip()).upper() for part in parts]
        if normalized and all(n == normalized[0] for n in normalized):
            sv_type = parts[0]
    me = re.search(r"(INS|DEL):(ME):?([A-Z0-9]+)?", sv_type, re.I)
    if me:
        abbrev, subtype = me.group(1).upper(), me.group(2).upper()
        element = me.group(3).upper() if me.group(3) else None
        if element is not None:
            element = "LINE1" if element == "L1" else element
            if element in ("ALU", "HERV", "LINE1", "SVA"):
                subtype = element
        abbrev = f"{abbrev}_{subtype}"
    elif re.search("DEL", sv_type, re.I) and re.search("DUP", sv_type, re.I):
        abbrev = "CNV"
    elif re.search("DUP:TANDEM", sv_type, re.I):
        abbrev = "TDUP"
    elif re.search("CNV:TR", sv_type, re.I):
        abbrev = "TREP"
    elif re.search(r"CN=?[0-9]", sv_type, re.I):
        abbrev = "CNV"
        if re.match(r"^<?CN=?0>?$", sv_type):
            abbrev = "DEL"
        if re.match(r"^<?CN=?2>?$", sv_type):
            abbrev = "DUP"
    elif re.search("CNV", sv_type, re.I):
        abbrev = "CNV"
    elif re.search(r"[\[\]]|^\.|\.$", sv_type):
        abbrev = "BND"
    elif re.search(r"^<|>$", sv_type):
        abbrev = re.sub(r":.+", "", re.sub(r"<|>", "", sv_type))
    else:
        abbrev = sv_type
    return abbrev if abbrev in _SO_TERMS else None


def _info_fields(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in info.split(";"):
        if item and item != ".":
            key, _, value = item.partition("=")
            out[key] = value
    return out


def perl_skip_route(v: VcfVariant, max_sv_size: int = _MAX_SV_SIZE) -> str | None:
    """The route by which Ensembl VEP 116.2 drops a VCF record before its input buffer, or
    None when it keeps the record: ``unsupported_type`` (Parser/VCF.pm:477-481),
    ``incomplete_deletion`` (:555-557), ``oversize`` (Parser.pm:491-500 over the span
    `get_end` derives: SVLEN first, else END, else the start; a breakend's END is discarded).
    ``max_sv_size`` is the cap the reference ran with; -1 is no cap, as for the flag. A
    record the parser routes to the sequence path is kept."""
    alts = v.alt.split(",")
    info = _info_fields(v.info)
    if re.fullmatch(r"[ACGT]+", v.ref + "".join(alts)) or not (info.get("SVTYPE") or _SV_ROUTED.search(",".join(alts))):
        return None
    sv_type = "/".join(alts)
    if info.get("SVTYPE") and not _VCF44_SYMBOLIC.match(sv_type):
        sv_type = info["SVTYPE"]
    abbrev = perl_so_abbrev(sv_type)
    if abbrev is None:
        return "unsupported_type"
    start = v.pos + 1
    svlen, end_field = info.get("SVLEN"), info.get("END")
    if svlen not in (None, ""):
        end = start + abs(int(svlen.split(",")[0])) - 1
    elif end_field not in (None, "") and abbrev != "BND":
        end = int(end_field)
    else:
        end = start
    if "del" in _SO_TERMS[abbrev] and (start > end or (svlen in (None, "") and end_field in (None, ""))):
        return "incomplete_deletion"
    if max_sv_size != -1 and end - start > max_sv_size:
        return "oversize"
    return None


_OWN_END_ALLELES = frozenset({"N.", ".N", "chromosome_breakpoint"})


def record_alleles(v: VcfVariant) -> set[str]:
    """Every Allele string either engine writes for a VCF record, normalised as the
    comparator normalises them: the Sequence Ontology class of each symbolic ALT and of the
    whole ALT list (`copy_number_variation` for a DEL beside a DUP), the raw abbreviation of
    a type with no term as both engines spell it (`CPX`, `<NON_REF>`), each bracket breakend
    string with `chr` stripped inside the brackets, the own-end strings of a breakend record
    (`N.`, `.N`, `chromosome_breakpoint`, the REF base with the dot), the mate string a
    `<BND>` record's CHR2/END2 give, and a sequence ALT as written."""
    alleles: set[str] = set()
    info = _info_fields(v.info)
    alts = v.alt.split(",")
    for alt in alts:
        if alt.startswith("<") and alt.endswith(">"):
            abbrev = perl_so_abbrev(alt)
            if abbrev:
                alleles.add(_SO_TERMS[abbrev])
            else:
                raw = re.sub(r":.+", "", alt[1:-1])
                alleles.update({raw, f"<{raw}>", alt})
        elif _BRACKET_MATE_RE.search(alt):
            alleles.add(normalize_allele(alt))
            alleles |= _OWN_END_ALLELES | {f"{v.ref}.", f".{v.ref}"}
        elif re.fullmatch(r"\.?[A-Za-z]+\.?", alt) and (alt.startswith(".") or alt.endswith(".")):
            alleles.update({alt, "chromosome_breakpoint"})
        else:
            alleles.add(alt)
    whole = perl_so_abbrev(info["SVTYPE"] if info.get("SVTYPE") and not _VCF44_SYMBOLIC.match("/".join(alts)) else "/".join(alts))
    if whole:
        alleles.add(_SO_TERMS[whole])
    if whole == "BND" or info.get("SVTYPE") == "BND":
        alleles |= _OWN_END_ALLELES | {f"{v.ref}.", f".{v.ref}"}
        chr2 = info.get("CHR2")
        if chr2:
            # The mate string the parser builds from INFO (Parser/VCF.pm:488-497): the mate
            # position is END2 when present, else the record's END, which a breakend record
            # from a caller that writes no END2 carries as its mate position.
            chr2 = chr2[3:] if chr2.lower().startswith("chr") else chr2
            for pos2 in (info.get("END2"), info.get("POS2"), info.get("END")):
                if pos2:
                    alleles.add(f"N[{chr2}:{pos2}[")
    return alleles


def record_location(v: VcfVariant) -> tuple[str, int, int]:
    """The (chromosome, start, end) both engines write for a structural record: the base
    after POS, and the end `get_end` derives (SVLEN first, else END, else the start; a
    breakend's END discarded)."""
    info = _info_fields(v.info)
    chrom = v.location.split(":", 1)[0]
    start = v.pos + 1
    svlen, end_field = info.get("SVLEN"), info.get("END")
    is_bnd = "[" in v.alt or "]" in v.alt or v.alt.startswith(".") or v.alt.endswith(".") or info.get("SVTYPE") == "BND"
    if svlen not in (None, ""):
        end = start + abs(int(svlen.split(",")[0])) - 1
    elif end_field not in (None, "") and not is_bnd:
        end = int(end_field)
    else:
        end = start
    return chrom, start, end


def skipped_record_rows(
    input_variants: list[VcfVariant],
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    max_sv_size: int = _MAX_SV_SIZE,
) -> tuple[set[AnnotTuple], int, list[str]]:
    """The records the 116 reference dropped as `vep_skip`, their vep-rs rows, and every
    disagreement between the records and the outputs.

    THE SHAPE. 116.2 writes at least one row for every record it keeps and no row in any
    format for a record it marks `vep_skip`; vep-rs annotates all of them. A record is
    named by the input: `perl_skip_route` over its ALT and INFO. It is CONFIRMED skipped
    when the reference wrote no row under the record's own Location start and an allele the
    record produces (`record_alleles`); a kept record sharing the start, or the whole
    Location string, with a skipped one therefore never hides it and is never mistaken for
    it. When two records at one start produce the same allele, the end of the Location
    tells them apart.

    SET ASIDE: every vep-rs tuple at a confirmed record's start under one of its alleles
    (the same end rule applies). `<CNV:TR>` records are outside the class on both sides:
    the reference keeps them and writes them at another Location
    (`filter_cnv_tr_expansion_divergences`).

    PROBLEMS, one line each, make the rule's blind spots loud: a skip-route record whose
    own rows the reference wrote (a reference run under --dont_skip or another
    --max_sv_size, declared with --reference-max-sv-size); a vep-rs structural row at a
    Location the reference wrote nothing for that no skip-route record explains (a record
    dropped on a route 116.2 shares with 115.2, `start > end+1`); two records at one start
    whose alleles and Locations coincide, one skipped and one kept.

    The SNP/indel comparator has no input VCF and applies the output-only rule
    (`ReferenceRecordScope`): a Location with vep-rs structural rows and no reference row.

    Returns (excluded_rust, records, problems).
    """
    def located(tuples: set[AnnotTuple]) -> dict[tuple[str, int], list[tuple[AnnotTuple, int]]]:
        out: dict[tuple[str, int], list[tuple[AnnotTuple, int]]] = defaultdict(list)
        for t in tuples:
            span = _location_span(t[0])
            if span is not None:
                out[(span[0], span[1])].append((t, span[2]))
        return out

    perl_at = located(perl_tuples)
    rust_at = located(rust_tuples)
    records_at: dict[tuple[str, int], list[tuple[VcfVariant, str | None, set[str], int]]] = defaultdict(list)
    for v in input_variants:
        chrom, start, end = record_location(v)
        records_at[(chrom, start)].append((v, perl_skip_route(v, max_sv_size), record_alleles(v), end))

    def rows_of(bucket: list[tuple[AnnotTuple, int]], alleles: set[str], end: int, contested: bool) -> list[AnnotTuple]:
        return [t for t, t_end in bucket if t[1] in alleles and (not contested or t_end == end)]

    excluded: set[AnnotTuple] = set()
    problems: list[str] = []
    confirmed_starts: set[tuple[str, int]] = set()
    confirmed = 0
    for key, records in records_at.items():
        for v, route, alleles, end in records:
            if route is None:
                continue
            # A kept record at this start producing one of this record's alleles: the end
            # decides whose rows are whose, and equal ends cannot be told apart.
            rivals = [(rv, ra, re_) for rv, rr, ra, re_ in records if rv is not v and rr is None and ra & alleles]
            contested = bool(rivals)
            if any(re_ == end for _, _, re_ in rivals):
                twin = next(rv for rv, _, re_ in rivals if re_ == end)
                problems.append(
                    f"{v.vid} {v.ref}>{v.alt} and {twin.vid} {twin.ref}>{twin.alt} share {key[0]}:{key[1]}-{end} and an allele; "
                    "the skipped one cannot be told from the kept one"
                )
                continue
            kept = rows_of(perl_at.get(key, []), alleles, end, contested)
            if kept:
                problems.append(
                    f"{v.vid} {v.ref}>{v.alt} at {key[0]}:{key[1]} is on the {route} route, yet the reference wrote "
                    f"{len(kept)} row(s) under its own allele(s) ({', '.join(sorted({t[1] for t in kept}))}): a reference run "
                    "under --dont_skip or another --max_sv_size (pass --reference-max-sv-size)"
                )
                continue
            confirmed_starts.add(key)
            confirmed += 1
            excluded.update(rows_of(rust_at.get(key, []), alleles, end, contested))
    # The output-only rule's residue: a Location with vep-rs structural rows and no
    # reference row that no confirmed record explains.
    perl_locations = {t[0] for t in perl_tuples}
    for t in rust_tuples:
        span = _location_span(t[0])
        if (t[0] not in perl_locations and is_structural_allele(t[1]) and t[1] != _CNV_TR_RUST_ALLELE
                and span is not None and (span[0], span[1]) not in confirmed_starts):
            problems.append(f"{t[0]} {t[1]}: the reference wrote nothing at this Location and no input record on a skip route starts there")
    return excluded, confirmed, sorted(set(problems))


def filter_reference_skipped_records(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    input_variants: list[VcfVariant],
    max_sv_size: int = _MAX_SV_SIZE,
) -> tuple[set[AnnotTuple], int]:
    """The mask over `skipped_record_rows`: the vep-rs tuples of every record the reference
    dropped as `vep_skip`, and the number of such records."""
    excluded, records, _problems = skipped_record_rows(input_variants, perl_tuples, rust_tuples, max_sv_size)
    return excluded, records


def check_skip_routes(
    input_variants: list[VcfVariant],
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    max_sv_size: int = _MAX_SV_SIZE,
) -> list[str]:
    """The problems of `skipped_record_rows`; empty when the records and the outputs agree."""
    return skipped_record_rows(input_variants, perl_tuples, rust_tuples, max_sv_size)[2]


def _breakpoint_position(loc: str) -> int | None:
    span = _location_span(loc)
    return None if span is None else span[1]


def filter_bnd_own_rows_lost(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    transcripts_by_chr: dict[str, set[str]],
    transcript_spans: dict[tuple[str, str], TranscriptSpan],
) -> set[AnnotTuple]:
    """A bracket breakend whose own-end rows the reference lost to its mate chromosome.

    THE SHAPE. A record's slice is named by the first transcript that reaches it
    (AnnotationType/Transcript.pm 116.2:122, `$vf->{slice} ||= $tr->{slice}`), and the
    transcripts of a region already in memory come first (AnnotationSource.pm
    116.2:119-134). When a mate-chromosome region is in memory from the previous buffer
    and the record's own region is not, a mate-chromosome transcript names the slice,
    every own-chromosome transcript then fails the slice test (:128-139: neither the
    record's slice nor any breakend's is on the transcript's chromosome) and the
    reference writes nothing under the own-end allele while keeping the mate rows. The
    result depends on the order of the input, which no annotation should.

    THE PREDICATE, on a Location whose bracket allele names a mate on another
    chromosome: the reference wrote a Transcript tuple under the bracket allele and no
    tuple under any own-end allele; then every vep-rs-only Transcript tuple under an
    own-end allele whose transcript the cache files on the Location's chromosome within
    `_BREAKEND_FEATURE_DISTANCE` of the breakpoint, the reference's own admission rule
    (StructuralVariationOverlap.pm 116:144-160), is excluded. A transcript farther away
    than that, or on another chromosome, stays charged. With no cache it excludes
    nothing.
    """
    if not transcripts_by_chr or not transcript_spans:
        return set()
    perl_by_loc = _tuples_by_location(perl_tuples)
    excluded: set[AnnotTuple] = set()
    for loc, rust_here in _tuples_by_location(rust_tuples - perl_tuples).items():
        own = _chromosome_of_location(loc)
        pos = _breakpoint_position(loc)
        if own is None or pos is None:
            continue
        perl_here = perl_by_loc.get(loc, set())
        mates = {m for t in perl_here | rust_here if (m := _bracket_mate(t[1])) is not None}
        if not any(m[0] != own for m in mates):
            continue
        if any(is_own_breakend_allele(t[1]) for t in perl_here):
            continue
        if not any(_bracket_mate(t[1]) is not None and t[2] != "-" for t in perl_here):
            continue
        local = transcripts_by_chr.get(own, set())
        for t in rust_here:
            if not is_own_breakend_allele(t[1]) or t[2] == "-" or t[2] not in local:
                continue
            span = transcript_spans.get((own, t[2]))
            if span is None:
                continue
            if span[1] - _BREAKEND_FEATURE_DISTANCE <= pos <= span[2] + _BREAKEND_FEATURE_DISTANCE:
                excluded.add(t)
    return excluded


def filter_bnd_synonym_mate_rows(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    transcripts_by_chr: dict[str, set[str]],
) -> set[AnnotTuple]:
    """A bracket breakend whose mate is written with a chromosome synonym.

    THE SHAPE. The overlap test resolves the mate's name through the synonym table
    (InputBuffer.pm 116.2:348-358, `get_source_chr_name`), so a mate-chromosome
    transcript receives the record; the slice test then compares the breakend's raw name,
    `chr` stripped, with the cache's (AnnotationType/Transcript.pm 116.2:134), and
    `NC_000002.12` is not `2`, so the reference writes nothing under the mate allele
    while keeping the record's own rows.

    THE PREDICATE, on a Location with a bracket allele whose mate name, `chr` stripped,
    is not a chromosome the cache files: the reference wrote at least one tuple at the
    Location and none under that bracket allele; then every vep-rs-only Transcript tuple
    under that allele whose transcript the cache does not file on the Location's own
    chromosome is excluded (an own-chromosome transcript under the mate allele is not
    the shape). With no cache it excludes nothing.
    """
    if not transcripts_by_chr:
        return set()
    perl_by_loc = _tuples_by_location(perl_tuples)
    excluded: set[AnnotTuple] = set()
    for loc, rust_here in _tuples_by_location(rust_tuples - perl_tuples).items():
        own = _chromosome_of_location(loc)
        perl_here = perl_by_loc.get(loc)
        if own is None or not perl_here:
            continue
        perl_alleles = {t[1] for t in perl_here}
        local = transcripts_by_chr.get(own, set())
        for t in rust_here:
            mate = _bracket_mate(t[1])
            if mate is None or mate[0] in transcripts_by_chr or t[1] in perl_alleles:
                continue
            if t[2] != "-" and t[2] not in local:
                excluded.add(t)
    return excluded


def count_fasta_named_slice_intergenic(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
) -> int:
    """Records whose only reference row is `intergenic_variant` while vep-rs wrote
    Transcript rows and no intergenic row: the shape of a slice named after a FASTA
    whose sequence names differ from the cache's (Parser.pm 116.2:749 through
    BaseVEP.pm `get_slice`, every transcript failing AnnotationType/Transcript.pm:141).
    The condition is the run's configuration, so the count is reported and nothing is
    set aside; a run that shows it is re-run with matching sequence names."""
    perl_by_variant: dict[tuple[str, str], set[AnnotTuple]] = defaultdict(set)
    for t in perl_tuples:
        perl_by_variant[(t[0], t[1])].add(t)
    rust_by_variant: dict[tuple[str, str], set[AnnotTuple]] = defaultdict(set)
    for t in rust_tuples:
        rust_by_variant[(t[0], t[1])].add(t)
    n = 0
    for key, perl_here in perl_by_variant.items():
        rust_here = rust_by_variant.get(key)
        if not rust_here:
            continue
        if perl_here == {key + _INTERGENIC_TUPLE_TAIL} and all(t[2] != "-" for t in rust_here):
            n += 1
    return n


@dataclass
class MaskResult:
    """Every mask's exclusions on one pair of tuple sets, under one reference release.

    A layer a release does not apply is an empty set here, so the arithmetic below
    (the union of the vep-rs-side layers, the union of the reference-side layers) is
    the same under both releases and the report states which layers were applied.
    """

    release: str
    tsel_rust: set[AnnotTuple] = field(default_factory=set)
    tsel_perl: set[AnnotTuple] = field(default_factory=set)
    cnvtr_rust: set[AnnotTuple] = field(default_factory=set)
    cnvtr_perl: set[AnnotTuple] = field(default_factory=set)
    xchr_perl: set[AnnotTuple] = field(default_factory=set)
    xchr_unresolvable_scored: int = 0
    xchr_orphan_rust: set[AnnotTuple] = field(default_factory=set)
    mate_rust: set[AnnotTuple] = field(default_factory=set)
    swap_rust: set[AnnotTuple] = field(default_factory=set)
    swap_perl: set[AnnotTuple] = field(default_factory=set)
    swap_by_bucket: dict[str, int] = field(default_factory=dict)
    non_ref_rust: set[AnnotTuple] = field(default_factory=set)
    non_ref_perl: set[AnnotTuple] = field(default_factory=set)
    mate_read_rust: set[AnnotTuple] = field(default_factory=set)
    mate_read_perl: set[AnnotTuple] = field(default_factory=set)
    mate_read_by_kind: dict[str, int] = field(default_factory=dict)
    skipped_rust: set[AnnotTuple] = field(default_factory=set)
    skipped_records: int = 0
    skipped_problems: list[str] = field(default_factory=list)
    bnd_own_rust: set[AnnotTuple] = field(default_factory=set)
    bnd_synonym_rust: set[AnnotTuple] = field(default_factory=set)
    fasta_intergenic_records: int = 0

    def all_rust(self) -> set[AnnotTuple]:
        return (
            self.tsel_rust | self.cnvtr_rust | self.xchr_orphan_rust | self.mate_rust
            | self.swap_rust | self.non_ref_rust | self.mate_read_rust
            | self.skipped_rust | self.bnd_own_rust | self.bnd_synonym_rust
        )

    def all_perl(self) -> set[AnnotTuple]:
        return (
            self.tsel_perl | self.cnvtr_perl | self.xchr_perl | self.swap_perl
            | self.non_ref_perl | self.mate_read_perl
        )


def apply_masks(
    perl_tuples: set[AnnotTuple],
    rust_tuples: set[AnnotTuple],
    reference_release: str,
    transcripts_by_chr: dict[str, set[str]],
    transcript_spans: dict[tuple[str, str], TranscriptSpan],
    transcript_models: dict[tuple[str, str], TranscriptModel],
    scored_engine: str,
    input_variants: list[VcfVariant] | None = None,
    max_sv_size: int = _MAX_SV_SIZE,
) -> MaskResult:
    """Run the masks of `MASKS_BY_RELEASE[reference_release]` on one pair of tuple sets.

    The layers of the 115.2 release are the published comparator's, in its order;
    `filter_cross_chromosome_orphan_intergenic` reads `filter_cross_chromosome_divergences`'
    result. Under 116.2 the class-C arms and the class-D masks are not called at all, so a
    reference that kept a `vep_skip` record (a run under --dont_skip) is scored as what it
    is, and the 116.2 classes are applied instead; the skipped-record class reads the
    input records (``input_variants``) and records its disagreements with the outputs in
    ``skipped_problems`` for the caller to stop on.
    """
    if reference_release not in MASKS_BY_RELEASE:
        raise ValueError(f"unknown reference release {reference_release!r}; one of {REFERENCE_RELEASES}")
    r = MaskResult(release=reference_release)
    masks = MASKS_BY_RELEASE[reference_release]
    if "filter_intended_divergences" in masks:
        r.tsel_rust, r.tsel_perl = filter_intended_divergences(perl_tuples, rust_tuples)
    if "filter_cnv_tr_expansion_divergences" in masks:
        r.cnvtr_rust, r.cnvtr_perl = filter_cnv_tr_expansion_divergences(perl_tuples, rust_tuples)
    if "filter_cross_chromosome_divergences" in masks:
        _, r.xchr_perl, r.xchr_unresolvable_scored = filter_cross_chromosome_divergences(
            perl_tuples, rust_tuples, transcripts_by_chr, scored_engine
        )
    if "filter_cross_chromosome_orphan_intergenic" in masks:
        r.xchr_orphan_rust = filter_cross_chromosome_orphan_intergenic(
            perl_tuples, rust_tuples, r.xchr_perl
        )
    if "filter_giant_breakend_mate_divergences" in masks:
        r.mate_rust = filter_giant_breakend_mate_divergences(perl_tuples, rust_tuples, transcript_spans)
    if "filter_registry_swap_pairs" in masks:
        r.swap_rust, r.swap_perl, r.swap_by_bucket = filter_registry_swap_pairs(perl_tuples, rust_tuples)
    if "filter_non_ref_batch_divergences" in masks:
        r.non_ref_rust, r.non_ref_perl = filter_non_ref_batch_divergences(perl_tuples, rust_tuples)
    if "filter_breakend_mate_local_read_pairs" in masks:
        r.mate_read_rust, r.mate_read_perl, r.mate_read_by_kind = filter_breakend_mate_local_read_pairs(
            perl_tuples, rust_tuples, transcript_models
        )
    if "filter_reference_skipped_records" in masks:
        if input_variants is None:
            raise ValueError("the skipped-record mask reads the input records; pass input_variants")
        r.skipped_rust, r.skipped_records, r.skipped_problems = skipped_record_rows(
            input_variants, perl_tuples, rust_tuples, max_sv_size
        )
    if "filter_bnd_own_rows_lost" in masks:
        r.bnd_own_rust = filter_bnd_own_rows_lost(
            perl_tuples, rust_tuples, transcripts_by_chr, transcript_spans
        )
    if "filter_bnd_synonym_mate_rows" in masks:
        r.bnd_synonym_rust = filter_bnd_synonym_mate_rows(perl_tuples, rust_tuples, transcripts_by_chr)
    if reference_release == "116.2":
        r.fasta_intergenic_records = count_fasta_named_slice_intergenic(perl_tuples, rust_tuples)
    return r


def release_report_block(
    overall: MaskResult,
    per_file: dict[str, MaskResult],
    cache_supplied: bool,
) -> dict:
    """The 116.2 additions to the report's ``adjusted`` block: the release, the masks
    applied and the 115.2 masks retired with the reason, each 116.2 class's exclusions
    (overall and per file) and its definition, and whether the cache authority the two
    breakend masks need was supplied."""
    # The counts each class contributes, keyed by bucket; a rule without one is a rule this
    # comparator has no mask for, which the report refuses to describe as applied.
    excluded_by_bucket = {
        REFERENCE_SKIPPED_RECORD_BUCKET: len(overall.skipped_rust),
        BND_OWN_ROWS_LOST_BUCKET: len(overall.bnd_own_rust),
        BND_SYNONYM_MATE_ROWS_LOST_BUCKET: len(overall.bnd_synonym_rust),
        FASTA_NAMED_SLICE_INTERGENIC_BUCKET: 0,
    }
    missing = [r.bucket for r in REFERENCE_116_2_ONE_SIDED_RULES if r.bucket not in excluded_by_bucket]
    if missing:
        raise AssertionError(f"one-sided rule(s) without a mask in this comparator: {missing}")
    return {
        "reference_release": overall.release,
        "masks_applied": list(MASKS_BY_RELEASE[overall.release]),
        "masks_not_applied": dict(MASKS_RETIRED_IN_116_2),
        "excluded_reference_skipped_rust": len(overall.skipped_rust),
        "reference_skipped_records": overall.skipped_records,
        "excluded_reference_skipped_by_file": {
            k: {"rust": len(v.skipped_rust), "records": v.skipped_records}
            for k, v in sorted(per_file.items())
            if v.skipped_rust
        },
        "excluded_bnd_own_rows_lost_rust": len(overall.bnd_own_rust),
        "excluded_bnd_own_rows_lost_by_file": {
            k: len(v.bnd_own_rust) for k, v in sorted(per_file.items()) if v.bnd_own_rust
        },
        "excluded_bnd_synonym_mate_rows_rust": len(overall.bnd_synonym_rust),
        "excluded_bnd_synonym_mate_rows_by_file": {
            k: len(v.bnd_synonym_rust) for k, v in sorted(per_file.items()) if v.bnd_synonym_rust
        },
        "fasta_named_slice_intergenic_records": overall.fasta_intergenic_records,
        "breakend_masks_check_ran": cache_supplied,
        "divergence_classes": {
            rule.bucket: {
                "definition": rule.definition,
                # This comparator's verdict is the rule's under a cache authority.
                "excludes": rule.excludes_with_cache,
                "excluded": excluded_by_bucket[rule.bucket],
                "taxonomy_class": rule.taxonomy_class,
                "perl_citation": rule.perl_citation,
            }
            for rule in REFERENCE_116_2_ONE_SIDED_RULES
        },
    }


# Report generation


def metrics_to_dict(m: Metrics) -> dict:
    d = {
        "total_input": m.total_input,
        "perl_variants": m.perl_variants,
        "rust_variants": m.rust_variants,
        "perl_tuples": m.perl_tuples,
        "rust_tuples": m.rust_tuples,
        "intersection": m.intersection,
        "only_perl": m.only_perl,
        "only_rust": m.only_rust,
        "precision": round(m.precision, 6),
        "recall": round(m.recall, 6),
        "f1": round(m.f1, 6),
    }
    # Aggregate rows carry the union-based counts alongside the file-summed ones
    # so the cross-file duplicate count is auditable from the report rather than
    # only visible to whoever re-runs the comparator.
    if m.perl_tuples_union is not None:
        d["perl_tuples_union"] = m.perl_tuples_union
        d["rust_tuples_union"] = m.rust_tuples_union
        d["intersection_union"] = m.intersection_union
        d["cross_file_duplicate_perl"] = m.cross_file_duplicate_perl
        d["cross_file_duplicate_rust"] = m.cross_file_duplicate_rust
    return d


def write_json_report(
    output_path: Path,
    file_metrics: dict[str, Metrics],
    per_type_metrics: dict[str, Metrics],
    overall: Metrics,
    adjusted_overall: Metrics | None = None,
    adjusted_file_metrics: dict[str, Metrics] | None = None,
    excluded_count: int = 0,
    excluded_perl_count: int = 0,
    excluded_cnvtr_rust_count: int = 0,
    excluded_cnvtr_perl_count: int = 0,
    excluded_cnvtr_by_file: dict[str, int] | None = None,
    cnvtr_swap_pairs_total_rust: int = 0,
    cnvtr_swap_pairs_total_perl: int = 0,
    cnvtr_swap_pairs_total_by_file: dict[str, dict[str, int]] | None = None,
    excluded_xchr_perl_count: int = 0,
    excluded_xchr_by_file: dict[str, int] | None = None,
    cross_chromosome_check_ran: bool = False,
    scored_engine: str = "vep-rs",
    scored_absent_from_cache: int = 0,
    excluded_transcript_selection_rust_count: int = 0,
    excluded_transcript_selection_perl_count: int = 0,
    excluded_transcript_selection_by_file: dict[str, dict[str, int]] | None = None,
    transcript_selection_cnvtr_overlap_rust: int = 0,
    excluded_xchr_orphan_rust_count: int = 0,
    excluded_xchr_orphan_rust_by_file: dict[str, int] | None = None,
    excluded_mate_rust_count: int = 0,
    excluded_mate_rust_by_file: dict[str, int] | None = None,
    excluded_swap_rust_count: int = 0,
    excluded_swap_perl_count: int = 0,
    excluded_swap_by_bucket: dict[str, int] | None = None,
    excluded_swap_by_file: dict[str, dict[str, int]] | None = None,
    transcript_selection_xchr_overlap_perl: int = 0,
    excluded_non_ref_rust_count: int = 0,
    excluded_non_ref_perl_count: int = 0,
    excluded_non_ref_by_file: dict[str, dict[str, int]] | None = None,
    excluded_mate_read_rust_count: int = 0,
    excluded_mate_read_perl_count: int = 0,
    excluded_mate_read_by_kind: dict[str, int] | None = None,
    excluded_mate_read_by_file: dict[str, dict[str, int]] | None = None,
    release_block: dict | None = None,
) -> None:
    """Write the JSON report.

    The ``adjusted`` block carries, beside the F1 arithmetic, every count the mask
    accounting needs to be reconstructed: the excluded counts per layer
    (``excluded_*``) and, for the one PARTIALLY masked SV class, the class total it
    was excluded from (``cnvtr_swap_pairs_total_rust`` / ``_perl``, with
    ``cnvtr_swap_pairs_total_by_file`` giving ``{"rust": n, "perl": n}`` per input
    file). Read without the total, ``excluded_cnvtr_*`` presents a partial mask as a
    full exclusion. ``excluded_rust_extra`` is the UNION of the vep-rs-side layers, so
    the transcript-selection layer's own size is written beside it as
    ``excluded_transcript_selection_rust`` (with ``_perl`` and ``_by_file``), and the
    number of tuples both vep-rs-side layers claim as
    ``transcript_selection_cnvtr_overlap_rust``; the union equals the two layer sizes
    less that overlap.

    Five further layers write their own sizes: ``excluded_cross_chromosome_rust`` (the
    vep-rs intergenic row of a record whose Perl tuples the cross-chromosome mask
    removed), ``excluded_transcript_selection_mate_rust`` (the mate allele of a
    breakend above --max_sv_size that Perl never wrote),
    ``excluded_registry_swap_rust`` / ``_perl`` with ``excluded_registry_swap_by_bucket``
    (both members of every pair an SNP/indel excluding rule matched),
    ``excluded_non_ref_batch_perl`` / ``_rust`` (Perl's batch-derived transcript tuples on
    a `<NON_REF>` record and the vep-rs intergenic row they displaced) and
    ``excluded_breakend_mate_local_read_rust`` / ``_perl`` with ``_by_kind`` (both members
    of every mate-side pair Perl read at the local coordinate, split by whether the mate
    is on the local chromosome). The transcript-selection Perl arm and the
    cross-chromosome mask can claim one tuple together; that count is
    ``transcript_selection_cross_chromosome_overlap_perl``. ``cross_chromosome_check_ran``
    also gates the mate-allele arm and the mate-side local-read pairs, since all three
    read the cache.

    ``release_block`` (116.2 and later): ``reference_release``, ``masks_applied``,
    ``masks_not_applied`` with the reason each 115.2 mask is retired, the 116.2 classes'
    exclusions overall and per file, ``breakend_masks_check_ran`` and every class's
    definition under ``divergence_classes``. Absent from a 115.2 report.
    """
    report = {
        "overall": metrics_to_dict(overall),
        "per_file": {k: metrics_to_dict(v) for k, v in file_metrics.items()},
        "per_variant_type": {
            k: metrics_to_dict(v) for k, v in per_type_metrics.items()
        },
    }
    if adjusted_overall is not None:
        adj_dict = metrics_to_dict(adjusted_overall)
        adj_dict["excluded_divergences"] = excluded_count
        adj_dict["excluded_rust_extra"] = excluded_count
        adj_dict["excluded_perl_extra"] = excluded_perl_count
        adj_dict["excluded_cnvtr_rust"] = excluded_cnvtr_rust_count
        adj_dict["excluded_cnvtr_perl"] = excluded_cnvtr_perl_count
        adj_dict["excluded_transcript_selection_rust"] = excluded_transcript_selection_rust_count
        adj_dict["excluded_transcript_selection_perl"] = excluded_transcript_selection_perl_count
        adj_dict["transcript_selection_cnvtr_overlap_rust"] = transcript_selection_cnvtr_overlap_rust
        # The class TOTAL the two excluded counts are a fraction of. Always written,
        # zero included, so a report from a corpus with no <CNV:TR> records states
        # that rather than omitting the key.
        adj_dict["cnvtr_swap_pairs_total_rust"] = cnvtr_swap_pairs_total_rust
        adj_dict["cnvtr_swap_pairs_total_perl"] = cnvtr_swap_pairs_total_perl
        adj_dict["excluded_cross_chromosome_perl"] = excluded_xchr_perl_count
        adj_dict["excluded_cross_chromosome_rust"] = excluded_xchr_orphan_rust_count
        adj_dict["excluded_transcript_selection_mate_rust"] = excluded_mate_rust_count
        adj_dict["excluded_registry_swap_rust"] = excluded_swap_rust_count
        adj_dict["excluded_registry_swap_perl"] = excluded_swap_perl_count
        adj_dict["excluded_registry_swap_by_bucket"] = dict(excluded_swap_by_bucket or {})
        adj_dict["excluded_non_ref_batch_rust"] = excluded_non_ref_rust_count
        adj_dict["excluded_non_ref_batch_perl"] = excluded_non_ref_perl_count
        adj_dict["excluded_breakend_mate_local_read_rust"] = excluded_mate_read_rust_count
        adj_dict["excluded_breakend_mate_local_read_perl"] = excluded_mate_read_perl_count
        adj_dict["excluded_breakend_mate_local_read_by_kind"] = dict(excluded_mate_read_by_kind or {})
        adj_dict["transcript_selection_cross_chromosome_overlap_perl"] = (
            transcript_selection_xchr_overlap_perl
        )
        # Recorded so a report cannot be read as "no cross-chromosome tuples" when it
        # actually means "no cache was supplied to adjudicate them".
        adj_dict["cross_chromosome_check_ran"] = cross_chromosome_check_ran
        # Which engine's output was scored, and how many of its tuples named a transcript
        # the cache does not hold for the variant's chromosome. That count is zero by
        # construction for vep-rs (the filter raises instead), and non-zero for a
        # third-party engine on its own transcript set. Recorded rather than dropped so a
        # non-vep-rs run states what it tolerated; those tuples stay in raw AND adjusted.
        adj_dict["scored_engine"] = scored_engine
        adj_dict["scored_engine_transcripts_absent_from_cache"] = scored_absent_from_cache
        if excluded_xchr_by_file:
            adj_dict["excluded_cross_chromosome_by_file"] = {
                k: v for k, v in sorted(excluded_xchr_by_file.items()) if v
            }
        if excluded_cnvtr_by_file is not None:
            adj_dict["excluded_cnvtr_by_file"] = {
                k: v for k, v in sorted(excluded_cnvtr_by_file.items()) if v
            }
        if excluded_transcript_selection_by_file is not None:
            adj_dict["excluded_transcript_selection_by_file"] = {
                k: v
                for k, v in sorted(excluded_transcript_selection_by_file.items())
                if v["rust"] or v["perl"]
            }
        if cnvtr_swap_pairs_total_by_file is not None:
            adj_dict["cnvtr_swap_pairs_total_by_file"] = {
                k: v
                for k, v in sorted(cnvtr_swap_pairs_total_by_file.items())
                if v["rust"] or v["perl"]
            }
        if excluded_xchr_orphan_rust_by_file:
            adj_dict["excluded_cross_chromosome_rust_by_file"] = {
                k: v for k, v in sorted(excluded_xchr_orphan_rust_by_file.items()) if v
            }
        if excluded_mate_rust_by_file:
            adj_dict["excluded_transcript_selection_mate_by_file"] = {
                k: v for k, v in sorted(excluded_mate_rust_by_file.items()) if v
            }
        if excluded_swap_by_file:
            adj_dict["excluded_registry_swap_by_file"] = {
                k: v
                for k, v in sorted(excluded_swap_by_file.items())
                if v["rust"] or v["perl"]
            }
        if excluded_non_ref_by_file:
            adj_dict["excluded_non_ref_batch_by_file"] = {
                k: v for k, v in sorted(excluded_non_ref_by_file.items()) if v["rust"] or v["perl"]
            }
        if excluded_mate_read_by_file:
            adj_dict["excluded_breakend_mate_local_read_by_file"] = {
                k: v
                for k, v in sorted(excluded_mate_read_by_file.items())
                if v["rust"] or v["perl"]
            }
        if release_block:
            adj_dict.update(release_block)
        report["adjusted"] = adj_dict
    if adjusted_file_metrics is not None:
        report["per_file_adjusted"] = {
            k: metrics_to_dict(v) for k, v in adjusted_file_metrics.items()
        }
    with open(output_path, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")


def write_markdown_report(
    output_path: Path,
    file_metrics: dict[str, Metrics],
    per_type_metrics: dict[str, Metrics],
    overall: Metrics,
    adjusted_overall: Metrics | None = None,
    adjusted_file_metrics: dict[str, Metrics] | None = None,
    excluded_count: int = 0,
    excluded_perl_count: int = 0,
    excluded_cnvtr_rust_count: int = 0,
    excluded_cnvtr_by_file: dict[str, int] | None = None,
    assembly: str = "",
    cnvtr_swap_pairs_total_rust: int = 0,
    reference_release: str | None = None,
    masks_applied: tuple[str, ...] | None = None,
    release_masks: MaskResult | None = None,
) -> None:
    lines: list[str] = []
    lines.append("# SV Concordance Report")
    lines.append("")

    # Overall summary
    lines.append("## Overall")
    lines.append("")
    if reference_release is not None:
        lines.append(
            f"- **Reference**: Ensembl VEP {reference_release}; masks applied: "
            + ", ".join(f"`{m}`" for m in (masks_applied or ()))
        )
    if release_masks is not None:
        # The release's own classes, so the headline's "Rust-extra" total below reads
        # beside the part of it that is skipped records.
        lines.append(
            f"- **{reference_release} classes**: `{REFERENCE_SKIPPED_RECORD_BUCKET}` "
            f"{len(release_masks.skipped_rust):,} vep-rs rows on {release_masks.skipped_records:,} records set aside; "
            f"`{BND_OWN_ROWS_LOST_BUCKET}` {len(release_masks.bnd_own_rust):,} set aside; "
            f"`{BND_SYNONYM_MATE_ROWS_LOST_BUCKET}` {len(release_masks.bnd_synonym_rust):,} set aside; "
            f"`{FASTA_NAMED_SLICE_INTERGENIC_BUCKET}` {release_masks.fasta_intergenic_records:,} records counted, none set aside"
        )
    total_excluded = excluded_count + excluded_perl_count
    if adjusted_overall is not None and total_excluded > 0:
        lines.append(
            f"- **F1 (adjusted)**: {adjusted_overall.f1:.4f}  "
            f"*(excluded {excluded_count:,} Rust-extra + {excluded_perl_count:,} Perl-extra transcript divergences)*"
        )
    lines.append(f"- **F1 (raw)**: {overall.f1:.4f}")
    lines.append(f"- **Total input variants**: {overall.total_input}")
    lines.append(f"- **Perl tuples**: {overall.perl_tuples}")
    lines.append(f"- **Rust tuples**: {overall.rust_tuples}")
    lines.append(f"- **Intersection**: {overall.intersection}")
    lines.append(f"- **Only Perl**: {overall.only_perl}")
    lines.append(f"- **Only Rust**: {overall.only_rust}")
    if excluded_cnvtr_rust_count > 0:
        files = ", ".join(
            sorted(k for k, v in (excluded_cnvtr_by_file or {}).items() if v)
        )
        asm_label = f" {assembly.upper()}" if assembly else ""
        # "N of M": the mask is PARTIAL, and a caption carrying only N reads as a
        # full exclusion. M is the class total from `count_cnv_tr_swap_population`.
        of_total = (
            f" of {cnvtr_swap_pairs_total_rust:,}" if cnvtr_swap_pairs_total_rust else ""
        )
        lines.append("")
        lines.append(
            f"SV adjusted F1 additionally excludes the `<CNV:TR>` literal-expansion "
            f"class ({excluded_cnvtr_rust_count:,}{of_total} pairs{asm_label}, "
            f"{files} only). Raw F1 unaffected."
        )
    lines.append("")

    # Per-file table
    lines.append("## Per-File Metrics")
    lines.append("")
    if adjusted_file_metrics is not None:
        lines.append(
            "| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |"
        )
        lines.append(
            "|------|-------|-------------|-------------|-----------|--------|--------|"
        )
        for fname, m in sorted(file_metrics.items()):
            adj_m = adjusted_file_metrics.get(fname, m)
            adj_f1 = f"{adj_m.f1:.4f}" if adj_m.f1 != m.f1 else ""
            lines.append(
                f"| {fname} | {m.total_input} "
                f"| {m.perl_tuples} | {m.rust_tuples} | {m.intersection} "
                f"| {m.f1:.4f} | {adj_f1} |"
            )
    else:
        lines.append(
            "| File | Input | Perl Tuples | Rust Tuples | Intersect | Precision | Recall | F1 |"
        )
        lines.append(
            "|------|-------|-------------|-------------|-----------|-----------|--------|------|"
        )
        for fname, m in sorted(file_metrics.items()):
            lines.append(
                f"| {fname} | {m.total_input} "
                f"| {m.perl_tuples} | {m.rust_tuples} | {m.intersection} "
                f"| {m.precision:.4f} | {m.recall:.4f} | {m.f1:.4f} |"
            )
    lines.append("")

    # Per-variant-type table
    lines.append("## Per-Variant-Type Metrics")
    lines.append("")
    lines.append(
        "| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |"
    )
    lines.append(
        "|--------------|-------|-------------|-------------|-----------|-----------|------|"
    )
    for vtype, m in sorted(per_type_metrics.items()):
        flag = " **" if m.f1 < 0.99 and m.f1 > 0 else ""
        lines.append(
            f"| {vtype} | {m.total_input} "
            f"| {m.perl_tuples} | {m.rust_tuples} "
            f"| {m.only_perl} | {m.only_rust} "
            f"| {m.f1:.4f}{flag} |"
        )
    lines.append("")

    with open(output_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_discordant_tsv(output_path: Path, records: list[DiscordantRecord]) -> None:
    header = "source_file\tvariant_type\tdirection\tlocation\tallele\tfeature\tfeature_type\tconsequence"
    with open(output_path, "w") as fh:
        fh.write(header + "\n")
        fh.writelines(
            f"{r.source_file}\t{r.variant_type}\t{r.direction}\t"
            f"{r.location}\t{r.allele}\t{r.feature}\t{r.feature_type}\t{r.consequence}\n"
            for r in records
        )


# Main


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare Perl VEP and vep-rs outputs across a directory of SV VCFs."
    )
    parser.add_argument(
        "--assembly",
        default="grch37",
        choices=["grch37", "grch38"],
        help="Assembly to compare (default: grch37). Sets default paths for input/perl/rust/output dirs.",
    )
    parser.add_argument(
        "--input-dir",
        default=None,
        help="Directory with input VCF files (default: tests/sv_validation/{assembly})",
    )
    parser.add_argument(
        "--perl-dir",
        default=None,
        help="Directory with Perl VEP output .txt files (default: tmp/sv_perl/{assembly})",
    )
    parser.add_argument(
        "--rust-dir",
        default=None,
        help="Directory with vep-rs output .txt files (default: tmp/sv_rust/{assembly})",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for concordance reports (default: tmp/sv_concordance/{assembly})",
    )
    parser.add_argument(
        "--vep-rs-cache",
        default=None,
        help=(
            "vep-rs JSON cache directory, the authority for which chromosome a transcript "
            "is on. Enables the cross-chromosome mask, which removes Perl tuples naming a "
            "transcript the variant's chromosome does not contain. Omit it and that mask "
            "is INACTIVE; the report records which, so an unadjudicated run cannot read "
            "as a clean one."
        ),
    )
    parser.add_argument(
        "--scored-engine",
        default="vep-rs",
        help=(
            "Which engine produced --rust-dir. Two things read it. (1) The cross-chromosome "
            "filter's scored-engine sanity assertion: vep-rs naming a transcript the cache "
            "does not hold for the variant's chromosome is a vep-rs bug and raises, while a "
            "third-party engine building its own transcript set from a different annotation "
            "release does it legitimately and is counted instead (reported as "
            "scored_engine_transcripts_absent_from_cache). (2) Allele normalisation: for "
            "'fastvep' the scored side's raw symbolic ALT tokens (<DEL>, <DUP:TANDEM>, "
            "<INS:ME:ALU>, <CN0>, ...) are mapped to the SO class strings Perl VEP and "
            "vep-rs write in the Allele column (deletion, tandem_duplication, "
            "Alu_insertion, deletion, ...), without which every symbolic-allele tuple on "
            "that arm is a guaranteed miss. The Perl side is never mapped. Pass the "
            "engine's name, e.g. fastvep, to score a non-vep-rs output."
        ),
    )
    parser.add_argument(
        "--reference-release",
        choices=REFERENCE_RELEASES,
        default=DEFAULT_REFERENCE_RELEASE,
        help=(
            "The Ensembl VEP release that produced --perl-dir, which selects the masks "
            "(MASKS_BY_RELEASE). Under 115.2 they are the published comparator's; under "
            "116.2 (default) the class-C arms and the two cross-chromosome masks are not "
            "applied, a record the reference dropped as vep_skip is set aside whole, and the "
            "two breakend shapes of 116.1's chromosome-aware overlap are set aside under "
            "--vep-rs-cache. The report names the release and the masks it applied. "
            "run_clone_measurement.sh derives the value from the reference set's provenance."
        ),
    )
    parser.add_argument(
        "--reference-max-sv-size",
        type=int,
        default=_MAX_SV_SIZE,
        help=(
            "The --max_sv_size the reference ran with (default 10000000, Ensembl VEP's own; -1 "
            "when it disabled the cap). Read by the 116.2 skipped-record cross-check, which "
            "refuses a reference that kept a record on a skip route or dropped one that is not."
        ),
    )
    args = parser.parse_args()

    # Set assembly-dependent defaults
    asm = args.assembly
    if args.input_dir is None:
        args.input_dir = f"tests/sv_validation/{asm}"
    if args.perl_dir is None:
        args.perl_dir = f"tmp/sv_perl/{asm}"
    if args.rust_dir is None:
        args.rust_dir = f"tmp/sv_rust/{asm}"
    if args.output_dir is None:
        args.output_dir = f"tmp/sv_concordance/{asm}"

    input_dir = Path(args.input_dir)
    perl_dir = Path(args.perl_dir)
    rust_dir = Path(args.rust_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.is_dir():
        print(
            f"ERROR: [compare_sv_concordance] input-dir not found: {input_dir}",
            file=sys.stderr,
        )
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    # The transcript-to-chromosome authority for the cross-chromosome mask. A supplied
    # path that yields no chromosomes is a hard error rather than a silent skip: the mask
    # would then report zero exclusions, which is indistinguishable from a clean run.
    transcripts_by_chr: dict[str, set[str]] = {}
    transcript_spans: dict[tuple[str, str], TranscriptSpan] = {}
    transcript_models: dict[tuple[str, str], TranscriptModel] = {}
    if args.vep_rs_cache:
        transcripts_by_chr, transcript_spans, transcript_models = load_cache_authority(
            args.vep_rs_cache
        )
        if not transcripts_by_chr:
            print(
                f"ERROR: [compare_sv_concordance] --vep-rs-cache "
                f"{args.vep_rs_cache} yielded no transcripts; the cross-chromosome mask "
                f"cannot adjudicate and a zero count would read as a clean run",
                file=sys.stderr,
            )
            sys.exit(1)
        print(
            f"cross-chromosome mask active over "
            f"{len(transcripts_by_chr)} chromosomes "
            f"({sum(len(v) for v in transcripts_by_chr.values()):,} transcripts)"
        )

    # Find input VCF files (all .vcf.gz, then fall back to .vcf)
    vcf_files = sorted(input_dir.glob("*.vcf.gz"))
    if not vcf_files:
        vcf_files = sorted(input_dir.glob("*.vcf"))
    if not vcf_files:
        print(
            f"ERROR: [compare_sv_concordance] no VCF files found in {input_dir}",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Found {len(vcf_files)} input VCF files")

    # Accumulators for overall + per-type
    all_perl_tuples: set[AnnotTuple] = set()
    all_rust_tuples: set[AnnotTuple] = set()
    all_perl_locations: set[str] = set()
    all_rust_locations: set[str] = set()
    all_input_variants: list[VcfVariant] = []
    all_discordants: list[DiscordantRecord] = []
    file_metrics: dict[str, Metrics] = {}
    # Per-file-summed tuple counts, the aggregation that matches `total_input`.
    # See the accumulation comment below for why these exist alongside the unions.
    perl_tuples_filesum = 0
    rust_tuples_filesum = 0
    intersection_filesum = 0

    # Per-file parsed data, kept for the adjusted metrics computation
    file_data: dict[
        str,
        tuple[
            list[VcfVariant], set[AnnotTuple], set[str], set[AnnotTuple], set[str], dict[str, str]
        ],
    ] = {}
    all_location_types: dict[str, str] = {}


    for vcf_path in vcf_files:
        basename = vcf_basename(vcf_path)
        print(f"  Processing {basename}...")

        # Parse input VCF
        input_variants = parse_vcf(vcf_path)

        # Parse VEP outputs
        perl_path = perl_dir / f"{basename}.txt"
        rust_path = rust_dir / f"{basename}.txt"

        # The Perl reference is parsed with the default engine so the fastVEP
        # symbolic-allele mapping can only ever touch the scored side.
        perl_tuples, perl_locs = parse_vep_output(perl_path)
        rust_tuples, rust_locs = parse_vep_output(rust_path, engine=args.scored_engine)

        if not perl_path.exists():
            print(f"    WARNING: Perl output missing: {perl_path}")
        if not rust_path.exists():
            print(f"    WARNING: Rust output missing: {rust_path}")

        # Per-file metrics
        m = compute_metrics(
            input_variants, perl_tuples, perl_locs, rust_tuples, rust_locs
        )
        file_metrics[basename] = m
        print(
            f"    Input={m.total_input} Perl={m.perl_tuples} Rust={m.rust_tuples} "
            f"Intersect={m.intersection} P={m.precision:.4f} R={m.recall:.4f} F1={m.f1:.4f}"
        )

        # The printed Locations of this file's records, by type, for the per-type table
        # and the discordant lists.
        location_types = location_types_from_outputs(input_variants, perl_path, rust_path)
        all_location_types.update(location_types)

        # Collect discordants
        disc = collect_discordants(
            basename, input_variants, perl_tuples, rust_tuples, location_types
        )
        all_discordants.extend(disc)

        # Store per-file data for adjusted metrics
        file_data[basename] = (
            input_variants,
            perl_tuples,
            perl_locs,
            rust_tuples,
            rust_locs,
            location_types,
        )

        # Accumulate for overall.
        #
        # The tuple key carries no source filename, so a tuple emitted from two
        # different VCFs of the same assembly collapses to one member of these
        # unions, while `total_input` (a list extend) keeps summing; an aggregate
        # row built from both mixes two incompatible aggregations, and the
        # cross-file duplicates move raw F1 at the sixth decimal on GRCh37 and the
        # third on GRCh38.
        #
        # The unions stay as-is because every mask predicate consumes bare 5-tuples;
        # the file-scoped counts below are what BOTH aggregate rows use, so the two
        # bases are reported side by side rather than silently mixed.
        all_perl_tuples |= perl_tuples
        all_rust_tuples |= rust_tuples
        all_perl_locations |= perl_locs
        all_rust_locations |= rust_locs
        all_input_variants.extend(input_variants)
        perl_tuples_filesum += len(perl_tuples)
        rust_tuples_filesum += len(rust_tuples)
        intersection_filesum += len(perl_tuples & rust_tuples)

    # Overall metrics (raw). Computed on the file-scoped sums so the tuple counts
    # aggregate the same way `total_input` does; the union-based figures are
    # reported alongside as `*_union` so the difference is visible rather than
    # implicit.
    overall = compute_metrics(
        all_input_variants,
        all_perl_tuples,
        all_perl_locations,
        all_rust_tuples,
        all_rust_locations,
    )
    overall.perl_tuples_union = overall.perl_tuples
    overall.rust_tuples_union = overall.rust_tuples
    overall.intersection_union = overall.intersection
    overall.cross_file_duplicate_perl = perl_tuples_filesum - overall.perl_tuples
    overall.cross_file_duplicate_rust = rust_tuples_filesum - overall.rust_tuples
    overall.perl_tuples = perl_tuples_filesum
    overall.rust_tuples = rust_tuples_filesum
    overall.intersection = intersection_filesum
    overall.only_perl = perl_tuples_filesum - intersection_filesum
    overall.only_rust = rust_tuples_filesum - intersection_filesum
    overall.recompute_prf()

    # Per-variant-type metrics (across all files)
    per_type = compute_per_type_metrics(
        all_input_variants,
        all_perl_tuples,
        all_perl_locations,
        all_rust_tuples,
        all_rust_locations,
        all_location_types,
    )

    # Adjusted metrics: the masks of the selected reference release (`MASKS_BY_RELEASE`),
    # run once on the union sets and once per file.
    masks = apply_masks(
        all_perl_tuples,
        all_rust_tuples,
        args.reference_release,
        transcripts_by_chr,
        transcript_spans,
        transcript_models,
        args.scored_engine,
        all_input_variants,
        args.reference_max_sv_size,
    )
    if masks.skipped_problems:
        # The skipped-record class against the input records: a disagreement is a
        # reference run under other flags, a record dropped on a route the registry does
        # not know, or two records the outputs cannot tell apart, and the comparison
        # stops rather than scoring over it.
        for line in masks.skipped_problems:
            print(f"ERROR: [compare_sv_concordance] skipped-record cross-check: {line}", file=sys.stderr)
        sys.exit(1)
    excluded_rust_overall = masks.tsel_rust
    excluded_perl_overall = masks.tsel_perl
    cnvtr_excluded_rust, cnvtr_excluded_perl = masks.cnvtr_rust, masks.cnvtr_perl
    # The class TOTAL those two exclusions are a fraction of, measured on the same
    # union sets the mask ran on, so the reported "N of M" shares one basis.
    cnvtr_total_rust, cnvtr_total_perl = count_cnv_tr_swap_population(
        all_perl_tuples, all_rust_tuples
    )
    xchr_excluded_perl, xchr_unresolvable_scored = masks.xchr_perl, masks.xchr_unresolvable_scored
    xchr_orphan_rust = masks.xchr_orphan_rust
    mate_excluded_rust = masks.mate_rust
    swap_excluded_rust, swap_excluded_perl, swap_by_bucket = (
        masks.swap_rust, masks.swap_perl, masks.swap_by_bucket
    )
    non_ref_excluded_rust, non_ref_excluded_perl = masks.non_ref_rust, masks.non_ref_perl
    mate_read_rust, mate_read_perl, mate_read_by_kind = (
        masks.mate_read_rust, masks.mate_read_perl, masks.mate_read_by_kind
    )
    # The transcript-selection layer's own size, kept before the union so the report
    # states the class exactly beside the union the adjusted denominators subtract.
    tsel_rust_count = len(excluded_rust_overall)
    tsel_perl_count = len(excluded_perl_overall)
    tsel_overlap_rust = len(excluded_rust_overall & cnvtr_excluded_rust)
    tsel_xchr_overlap_perl = len(excluded_perl_overall & xchr_excluded_perl)
    excluded_rust_overall = masks.all_rust()
    excluded_perl_overall = masks.all_perl()
    excluded_count = len(excluded_rust_overall)
    excluded_perl_count = len(excluded_perl_overall)
    adj_rust_tuples = all_rust_tuples - excluded_rust_overall
    adj_perl_tuples = all_perl_tuples - excluded_perl_overall
    adjusted_overall = compute_metrics(
        all_input_variants,
        adj_perl_tuples,
        all_perl_locations,
        adj_rust_tuples,
        all_rust_locations,
    )

    # Per-file adjusted metrics (from stored per-file data, no re-parsing)
    adjusted_file_metrics: dict[str, Metrics] = {}
    # The one-sided tuples the adjusted row counts, collected from the same per-file
    # sets it is computed from so the file and the row cannot disagree.
    adjusted_discordants: list[DiscordantRecord] = []
    cnvtr_rust_by_file: dict[str, int] = {}
    cnvtr_perl_by_file: dict[str, int] = {}
    cnvtr_total_by_file: dict[str, dict[str, int]] = {}
    xchr_perl_by_file: dict[str, int] = {}
    xchr_orphan_rust_by_file: dict[str, int] = {}
    mate_rust_by_file: dict[str, int] = {}
    swap_by_file: dict[str, dict[str, int]] = {}
    non_ref_by_file: dict[str, dict[str, int]] = {}
    mate_read_by_file: dict[str, dict[str, int]] = {}
    tsel_by_file: dict[str, dict[str, int]] = {}
    masks_by_file: dict[str, MaskResult] = {}
    adj_perl_filesum = 0
    adj_rust_filesum = 0
    adj_intersection_filesum = 0
    for basename, (variants_f, perl_t, perl_l, rust_t, rust_l, loc_types_f) in file_data.items():
        masks_f = apply_masks(
            perl_t,
            rust_t,
            args.reference_release,
            transcripts_by_chr,
            transcript_spans,
            transcript_models,
            args.scored_engine,
            variants_f,
            args.reference_max_sv_size,
        )
        masks_by_file[basename] = masks_f
        cnvtr_rust_by_file[basename] = len(masks_f.cnvtr_rust)
        cnvtr_perl_by_file[basename] = len(masks_f.cnvtr_perl)
        tsel_by_file[basename] = {"rust": len(masks_f.tsel_rust), "perl": len(masks_f.tsel_perl)}
        tot_rust_f, tot_perl_f = count_cnv_tr_swap_population(perl_t, rust_t)
        cnvtr_total_by_file[basename] = {"rust": tot_rust_f, "perl": tot_perl_f}
        xchr_perl_by_file[basename] = len(masks_f.xchr_perl)
        xchr_orphan_rust_by_file[basename] = len(masks_f.xchr_orphan_rust)
        mate_rust_by_file[basename] = len(masks_f.mate_rust)
        swap_by_file[basename] = {
            "rust": len(masks_f.swap_rust), "perl": len(masks_f.swap_perl), **masks_f.swap_by_bucket
        }
        non_ref_by_file[basename] = {"rust": len(masks_f.non_ref_rust), "perl": len(masks_f.non_ref_perl)}
        mate_read_by_file[basename] = {
            "rust": len(masks_f.mate_read_rust), "perl": len(masks_f.mate_read_perl), **masks_f.mate_read_by_kind
        }
        adj_rust_f = rust_t - masks_f.all_rust()
        adj_perl_f = perl_t - masks_f.all_perl()
        adjusted_file_metrics[basename] = compute_metrics(
            variants_f, adj_perl_f, perl_l, adj_rust_f, rust_l
        )
        adjusted_discordants.extend(
            collect_discordants(basename, variants_f, adj_perl_f, adj_rust_f, loc_types_f)
        )
        adj_perl_filesum += len(adj_perl_f)
        adj_rust_filesum += len(adj_rust_f)
        adj_intersection_filesum += len(adj_perl_f & adj_rust_f)

    # Put the ADJUSTED aggregate on the same file-summed basis as the raw one, so a
    # suite's two F1 columns share one denominator. File-summing is the aggregation
    # rule throughout: counts are pooled and F1 computed once from the pooled
    # totals, the same rule that makes the SNP/indel micro-F1 the sum of the six
    # per-suite rows.
    #
    # The masks operate on unions (every predicate consumes bare 5-tuples),
    # so the union figures are preserved as `*_union` for audit, exactly as the
    # raw row does.
    adjusted_overall.perl_tuples_union = adjusted_overall.perl_tuples
    adjusted_overall.rust_tuples_union = adjusted_overall.rust_tuples
    adjusted_overall.intersection_union = adjusted_overall.intersection
    adjusted_overall.cross_file_duplicate_perl = (
        adj_perl_filesum - adjusted_overall.perl_tuples
    )
    adjusted_overall.cross_file_duplicate_rust = (
        adj_rust_filesum - adjusted_overall.rust_tuples
    )
    adjusted_overall.perl_tuples = adj_perl_filesum
    adjusted_overall.rust_tuples = adj_rust_filesum
    adjusted_overall.intersection = adj_intersection_filesum
    adjusted_overall.only_perl = adj_perl_filesum - adj_intersection_filesum
    adjusted_overall.only_rust = adj_rust_filesum - adj_intersection_filesum
    adjusted_overall.recompute_prf()

    # Write reports
    json_path = output_dir / "concordance_report.json"
    md_path = output_dir / "concordance_report.md"
    tsv_path = output_dir / "discordant.tsv"
    adjusted_tsv_path = output_dir / "discordant_adjusted.tsv"

    # The 115.2 report is the published comparator's and keeps its exact shape; a later
    # release's report names itself, the masks it applied and the 116.2 classes.
    release_reported = args.reference_release != "115.2"
    release_block = (
        release_report_block(masks, masks_by_file, bool(transcripts_by_chr)) if release_reported else None
    )
    write_json_report(
        json_path,
        file_metrics,
        per_type,
        overall,
        adjusted_overall=adjusted_overall,
        adjusted_file_metrics=adjusted_file_metrics,
        excluded_count=excluded_count,
        excluded_perl_count=excluded_perl_count,
        excluded_cnvtr_rust_count=len(cnvtr_excluded_rust),
        excluded_cnvtr_perl_count=len(cnvtr_excluded_perl),
        excluded_cnvtr_by_file=cnvtr_rust_by_file,
        excluded_transcript_selection_rust_count=tsel_rust_count,
        excluded_transcript_selection_perl_count=tsel_perl_count,
        excluded_transcript_selection_by_file=tsel_by_file,
        transcript_selection_cnvtr_overlap_rust=tsel_overlap_rust,
        cnvtr_swap_pairs_total_rust=cnvtr_total_rust,
        cnvtr_swap_pairs_total_perl=cnvtr_total_perl,
        cnvtr_swap_pairs_total_by_file=cnvtr_total_by_file,
        excluded_xchr_perl_count=len(xchr_excluded_perl),
        excluded_xchr_by_file=xchr_perl_by_file,
        cross_chromosome_check_ran=bool(transcripts_by_chr),
        scored_engine=args.scored_engine,
        scored_absent_from_cache=xchr_unresolvable_scored,
        excluded_xchr_orphan_rust_count=len(xchr_orphan_rust),
        excluded_xchr_orphan_rust_by_file=xchr_orphan_rust_by_file,
        excluded_mate_rust_count=len(mate_excluded_rust),
        excluded_mate_rust_by_file=mate_rust_by_file,
        excluded_swap_rust_count=len(swap_excluded_rust),
        excluded_swap_perl_count=len(swap_excluded_perl),
        excluded_swap_by_bucket=swap_by_bucket,
        excluded_swap_by_file=swap_by_file,
        transcript_selection_xchr_overlap_perl=tsel_xchr_overlap_perl,
        excluded_non_ref_rust_count=len(non_ref_excluded_rust),
        excluded_non_ref_perl_count=len(non_ref_excluded_perl),
        excluded_non_ref_by_file=non_ref_by_file,
        excluded_mate_read_rust_count=len(mate_read_rust),
        excluded_mate_read_perl_count=len(mate_read_perl),
        excluded_mate_read_by_kind=mate_read_by_kind,
        excluded_mate_read_by_file=mate_read_by_file,
        release_block=release_block,
    )
    write_markdown_report(
        md_path,
        file_metrics,
        per_type,
        overall,
        adjusted_overall=adjusted_overall,
        adjusted_file_metrics=adjusted_file_metrics,
        excluded_count=excluded_count,
        excluded_perl_count=excluded_perl_count,
        excluded_cnvtr_rust_count=len(cnvtr_excluded_rust),
        excluded_cnvtr_by_file=cnvtr_rust_by_file,
        assembly=asm,
        cnvtr_swap_pairs_total_rust=cnvtr_total_rust,
        reference_release=args.reference_release if release_reported else None,
        masks_applied=MASKS_BY_RELEASE[args.reference_release] if release_reported else None,
        release_masks=masks if release_reported else None,
    )
    write_discordant_tsv(tsv_path, all_discordants)
    write_discordant_tsv(adjusted_tsv_path, adjusted_discordants)

    print()
    total_excluded = excluded_count + excluded_perl_count
    if total_excluded > 0:
        print(
            f"Overall (adjusted): P={adjusted_overall.precision:.4f} "
            f"R={adjusted_overall.recall:.4f} F1={adjusted_overall.f1:.4f}  "
            f"(excluded {excluded_count:,} Rust-extra + {excluded_perl_count:,} Perl-extra transcript divergences)"
        )
    print(
        f"Overall (raw):      P={overall.precision:.4f} R={overall.recall:.4f} F1={overall.f1:.4f}"
    )
    print(f"  Perl tuples: {overall.perl_tuples}")
    print(f"  Rust tuples: {overall.rust_tuples}")
    print(f"  Intersection: {overall.intersection}")
    print(f"  Only Perl: {overall.only_perl}")
    print(f"  Only Rust: {overall.only_rust}")
    if total_excluded > 0:
        print(f"  Excluded (Rust-extra transcripts): {excluded_count}")
        print(f"  Excluded (Perl-extra transcripts): {excluded_perl_count}")
    if cnvtr_excluded_rust or cnvtr_excluded_perl or cnvtr_total_rust or cnvtr_total_perl:
        print(
            f"  Excluded (<CNV:TR> literal expansion): "
            f"{len(cnvtr_excluded_rust)} of {cnvtr_total_rust} Rust / "
            f"{len(cnvtr_excluded_perl)} of {cnvtr_total_perl} Perl class tuples"
        )
        for fname, n in sorted(cnvtr_rust_by_file.items()):
            tot = cnvtr_total_by_file.get(fname, {"rust": 0, "perl": 0})
            if n or tot["rust"] or tot["perl"]:
                print(f"    {fname}: {n} of {tot['rust']} Rust / {tot['perl']} Perl")
    print()
    print("Per-variant-type summary:")
    for vtype, m in sorted(per_type.items()):
        flag = " <--" if m.f1 < 0.99 else ""
        print(f"  {vtype:20s}  n={m.total_input:5d}  F1={m.f1:.4f}{flag}")
    print()
    # Per-file summary with adjusted column
    if total_excluded > 0:
        print("Per-file summary:")
        for fname in sorted(file_metrics.keys()):
            raw_f1 = file_metrics[fname].f1
            adj_f1 = adjusted_file_metrics[fname].f1
            delta = adj_f1 - raw_f1
            if delta > 0.0001:
                print(
                    f"  {fname:35s}  Raw F1={raw_f1:.4f}  Adj F1={adj_f1:.4f}  "
                    f"(+{delta:.4f})"
                )
        print()
    if release_reported:
        print(
            f"Reference release {args.reference_release}: "
            f"{len(masks.skipped_rust)} vep-rs tuple(s) on {masks.skipped_records} record(s) the reference skipped; "
            f"{len(masks.bnd_own_rust)} own-end breakend row(s) and {len(masks.bnd_synonym_rust)} synonym-mate row(s) set aside; "
            f"{masks.fasta_intergenic_records} FASTA-named-slice record(s) counted"
        )
    print(f"Reports written to: {output_dir}/")
    print(f"  {json_path.name}")
    print(f"  {md_path.name}")
    print(f"  {tsv_path.name} ({len(all_discordants)} discordant records)")
    print(
        f"  {adjusted_tsv_path.name} ({len(adjusted_discordants)} one-sided records "
        f"surviving every mask)"
    )


if __name__ == "__main__":
    main()
