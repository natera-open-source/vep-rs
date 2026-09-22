#!/usr/bin/env python3
"""Compare Perl VEP and vep-rs outputs across a directory of SV VCFs.

The default input is the synthetic corpus under ``tests/sv_validation/``; a
real-world SV set is scored the same way through ``--input-dir``. Produces
per-variant-type concordance metrics with JSON, Markdown, and TSV reports of
discordant annotation tuples.

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
from dataclasses import dataclass
from pathlib import Path

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
            variants.append(VcfVariant(chrom, pos, vid, ref, alt, vtype, location))
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
) -> list[DiscordantRecord]:
    records: list[DiscordantRecord] = []
    # Build location -> variant_type lookup
    loc_to_type: dict[str, str] = {}
    for v in input_variants:
        loc_to_type[v.location] = v.variant_type

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
) -> dict[str, Metrics]:
    """Group input variants by type and compute metrics for each group."""
    # Group variants by type
    type_variants: dict[str, list[VcfVariant]] = defaultdict(list)
    for v in input_variants:
        type_variants[v.variant_type].append(v)

    type_metrics: dict[str, Metrics] = {}
    for vtype, variants in sorted(type_variants.items()):
        # Filter tuples to only those whose Location matches a variant of this type
        type_locs = {v.location for v in variants}

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

# Perl VEP's --max_sv_size default (ensembl-vep Config.pm:310). A variant longer
# than this is flagged `vep_skip` at parse time (Parser.pm:493-497), excluded from
# cache-region loading (AnnotationSource.pm:238) but NOT from the annotation pass,
# which is the whole mechanism the transcript-selection mask exists for. Below the
# threshold Perl requests every region a variant overlaps, uncapped
# (AnnotationSource.pm:194-200), so no batch dependence arises and nothing should
# be masked.
_MAX_SV_SIZE = 10_000_000


def _span_exceeds_max_sv_size(loc: str) -> bool:
    """True if ``loc``'s span exceeds Perl's --max_sv_size default.

    ``loc`` is ``chrom:start-end`` for a ranged record or ``chrom:pos`` for a point
    one. A point location yields a span of 0 and is therefore never in scope, which
    is correct for everything except a giant inter-chromosomal breakend, whose mate
    coordinate lives in the Allele bracket rather than the Location. That case is
    documented at the call sites and left unmasked, the conservative direction.

    VEP's insertion convention puts ``end`` before ``start``, so the magnitude is
    taken rather than the signed difference.
    """
    parts = loc.split(":", 1)
    if len(parts) < 2:
        return False
    coords = parts[1]
    if "-" not in coords:
        return False
    start_s, _, end_s = coords.partition("-")
    try:
        return abs(int(end_s) - int(start_s)) > _MAX_SV_SIZE
    except ValueError:
        return False


def load_transcripts_by_chromosome(cache_dir: str | Path) -> dict[str, set[str]]:
    """Map each chromosome to the transcript ids the vep-rs JSON cache holds for it.

    The authority for "which chromosome is this transcript on". Both engines are run
    against this one cache, so a transcript it files under chr21 is a chr21 transcript
    and a transcript absent from chr21 is not one, whatever either engine's output says.

    Returns an empty mapping when the directory is absent or unreadable, and the caller
    treats that as "cannot judge" rather than "nothing to exclude".
    """
    root = Path(cache_dir) / "transcripts"
    out: dict[str, set[str]] = {}
    if not root.is_dir():
        return out
    for chr_dir in sorted(root.iterdir()):
        if not chr_dir.is_dir():
            continue
        ids: set[str] = set()
        for shard in sorted(chr_dir.glob("*.json")):
            try:
                data = json.loads(shard.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            # Storable-derived caches use {"chr": [tx, null, ...]}; native ones a list.
            if isinstance(data, dict):
                data = [x for v in data.values() if v for x in v if x]
            if not isinstance(data, list):
                continue
            for tx in data:
                if isinstance(tx, dict) and tx.get("stable_id"):
                    ids.add(tx["stable_id"])
        if ids:
            out[chr_dir.name] = ids
    return out


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

    THE MECHANISM IS NOT ESTABLISHED. The destinations suggest a lookup running past
    chr21's cache slice, but that is a hypothesis; what is established is that the
    contig-naming trigger cannot be the whole story, because the GRCh37 file's `CHR2=21`
    matches its `CHROM` and the defect fires there anyway.

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
    # (b) SCOPE. The mechanism requires a span above Perl's --max_sv_size default
    #     (10 Mb, ensembl-vep Config.pm:310). Below it Perl requests every region a
    #     variant overlaps, so no batch dependence exists and nothing is masked: a
    #     spurious transcript on a 100 bp deletion or a 1 bp point locus stays
    #     charged to vep-rs.
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
        if not _span_exceeds_max_sv_size(loc):
            # Below --max_sv_size Perl loads every region it needs, so a
            # transcript-set difference here is not the masked defect.
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
        if not _span_exceeds_max_sv_size(loc):
            continue
        excluded_perl.add(t)

    return excluded_rust, excluded_perl


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


def _location_span(loc: str) -> tuple[str, int, int] | None:
    """``chrom:start-end`` or ``chrom:pos`` as (chrom, start, end); None when unparseable."""
    m = re.fullmatch(r"([^:]+):(\d+)(?:-(\d+))?", loc)
    if not m:
        return None
    start = int(m.group(2))
    return m.group(1), start, int(m.group(3)) if m.group(3) else start


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
) -> None:
    lines: list[str] = []
    lines.append("# SV Concordance Report")
    lines.append("")

    # Overall summary
    lines.append("## Overall")
    lines.append("")
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
    if args.vep_rs_cache:
        transcripts_by_chr = load_transcripts_by_chromosome(args.vep_rs_cache)
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
        tuple[list[VcfVariant], set[AnnotTuple], set[str], set[AnnotTuple], set[str]],
    ] = {}

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

        # Collect discordants
        disc = collect_discordants(basename, input_variants, perl_tuples, rust_tuples)
        all_discordants.extend(disc)

        # Store per-file data for adjusted metrics
        file_data[basename] = (
            input_variants,
            perl_tuples,
            perl_locs,
            rust_tuples,
            rust_locs,
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
    )

    # Adjusted metrics: three layers of divergence exclusion.
    # Layer 1: Transcript-selection divergences (different transcripts for same variant)
    excluded_rust_overall, excluded_perl_overall = filter_intended_divergences(
        all_perl_tuples, all_rust_tuples
    )
    # Layer 2: <CNV:TR> literal-expansion vs symbolic annotation. Perl materializes
    # the symbolic tandem-repeat allele and re-runs point-variant splice predicates
    # on the expanded bases; vep-rs keeps it symbolic as feature_elongation.
    # 07_cnv_repeat only.
    cnvtr_excluded_rust, cnvtr_excluded_perl = filter_cnv_tr_expansion_divergences(
        all_perl_tuples, all_rust_tuples
    )
    # The class TOTAL those two exclusions are a fraction of, measured on the same
    # union sets the mask ran on, so the reported "N of M" shares one basis.
    cnvtr_total_rust, cnvtr_total_perl = count_cnv_tr_swap_population(
        all_perl_tuples, all_rust_tuples
    )
    # Layer 3: Perl naming transcripts that are not on the variant's chromosome.
    # Perl-side only, and `filter_cross_chromosome_divergences` raises rather than
    # excluding if vep-rs ever produces the shape.
    _, xchr_excluded_perl, xchr_unresolvable_scored = filter_cross_chromosome_divergences(
        all_perl_tuples, all_rust_tuples, transcripts_by_chr, args.scored_engine
    )
    # The transcript-selection layer's own size, kept before the union so the report
    # states the class exactly beside the union the adjusted denominators subtract.
    tsel_rust_count = len(excluded_rust_overall)
    tsel_perl_count = len(excluded_perl_overall)
    tsel_overlap_rust = len(excluded_rust_overall & cnvtr_excluded_rust)
    excluded_rust_overall |= cnvtr_excluded_rust
    excluded_perl_overall |= cnvtr_excluded_perl | xchr_excluded_perl
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
    cnvtr_rust_by_file: dict[str, int] = {}
    cnvtr_perl_by_file: dict[str, int] = {}
    cnvtr_total_by_file: dict[str, dict[str, int]] = {}
    xchr_perl_by_file: dict[str, int] = {}
    tsel_by_file: dict[str, dict[str, int]] = {}
    adj_perl_filesum = 0
    adj_rust_filesum = 0
    adj_intersection_filesum = 0
    for basename, (variants_f, perl_t, perl_l, rust_t, rust_l) in file_data.items():
        excluded_rust_f, excluded_perl_f = filter_intended_divergences(perl_t, rust_t)
        cnvtr_excl_rust_f, cnvtr_excl_perl_f = filter_cnv_tr_expansion_divergences(
            perl_t, rust_t
        )
        _, xchr_excl_perl_f, _ = filter_cross_chromosome_divergences(
            perl_t, rust_t, transcripts_by_chr, args.scored_engine
        )
        cnvtr_rust_by_file[basename] = len(cnvtr_excl_rust_f)
        cnvtr_perl_by_file[basename] = len(cnvtr_excl_perl_f)
        tsel_by_file[basename] = {"rust": len(excluded_rust_f), "perl": len(excluded_perl_f)}
        tot_rust_f, tot_perl_f = count_cnv_tr_swap_population(perl_t, rust_t)
        cnvtr_total_by_file[basename] = {"rust": tot_rust_f, "perl": tot_perl_f}
        xchr_perl_by_file[basename] = len(xchr_excl_perl_f)
        adj_rust_f = rust_t - excluded_rust_f - cnvtr_excl_rust_f
        adj_perl_f = perl_t - excluded_perl_f - cnvtr_excl_perl_f - xchr_excl_perl_f
        adjusted_file_metrics[basename] = compute_metrics(
            variants_f, adj_perl_f, perl_l, adj_rust_f, rust_l
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
    )
    write_discordant_tsv(tsv_path, all_discordants)

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
    print(f"Reports written to: {output_dir}/")
    print(f"  {json_path.name}")
    print(f"  {md_path.name}")
    print(f"  {tsv_path.name} ({len(all_discordants)} discordant records)")


if __name__ == "__main__":
    main()
