#!/usr/bin/env python3
"""Derive what a point breakend's mate-side Transcript row supports, from the vep-rs JSON cache.

A mate-side row is a Transcript row of a bracket-notation breakend whose Feature lies
within `MAX_DISTANCE_FROM_TRANSCRIPT` of the MATE coordinate on the mate chromosome
while the LOCAL coordinate does not itself select that transcript (another chromosome,
or farther than that distance on the same one). Ensembl VEP builds such a row for
every transcript near the mate (ensembl-variation StructuralVariationOverlap.pm:73-87,
`_close_to_feature`) and then evaluates every positional predicate on it with the
LOCAL variation feature (BaseVariationFeatureOverlapAllele.pm:257,273 hand the local
feature to each predicate; `upstream`/`downstream` at Utils/VariationEffect.pm:444-457
and `_bvfo_preds` at BaseVariationFeatureOverlapAllele.pm:454 read its start and end);
only `feature_truncation` (Utils/VariationEffect.pm:358) reads the mate breakend.

Two derivations, one function: `breakend_context(pos, model, mate_inside)` gives the
term set a point breakend at `pos` supports on the transcript, with `feature_truncation`
decided by whether the MATE lies inside the transcript (the one breakend-aware
predicate). Called with the mate position it is the set the mate breakend supports;
called with the local position it is Ensembl's read, whose empty result becomes the
default `intergenic_variant` on the row. `mate_side_derivation` applies both to a row
and returns them, or None when the row is not a mate-side row of a point breakend.

The transcript structure comes from the same JSON cache both engines annotated against:
exons, strand, biotype, the genomic CDS bounds mapped from the mapper's
`cdna_coding_start`/`cdna_coding_end` through its cDNA-to-genome pairs, and the mature
miRNA ranges of the `miRNA` attribute (cDNA coordinates). Nothing here reads the
Rust crates; a test and the comparator share this one module.

What the derivation does not model, so that a corpus reaching one of these edges is
adjudicated rather than absorbed (the predicate excludes a pair only when both engines'
sets equal the derivations, so an unmodelled edge leaves the pair charged):

- a breakend inside a start or stop codon: Ensembl fires `start_lost` and `stop_lost`
  on span overlap (Utils/VariationEffect.pm:884-893, the stop analogue) and
  `start_retained_variant` through `_ins_del_start_altered` returning 0 on a structural
  allele (:1037); no rule here;
- the frameshift-intron exon lookup: a transcript with an intron of 13 bases or fewer
  stretches its exon lookup by 12 bases (BaseTranscriptVariation.pm:911-914); the exon
  and intron tests here are exact;
- a multi-exon miRNA: the mature range is a cDNA interval that is not read off a
  multi-exon span, so the non-coding exon term stands;
- the exact `MAX_DISTANCE_FROM_TRANSCRIPT` bound: `_close_to_feature` expands the
  feature slice by 5,000 on each side and tests overlap inclusively, matched here by
  `TranscriptModel.near`; a breakend exactly 5,000 bases outside is inside the gate on
  both readings;
- a same-chromosome breakend whose local coordinate is also within the gate of the
  transcript: `mate_side_derivation` returns None (the row is a local-side row) and a
  divergence there stays charged.
"""

from __future__ import annotations

import gzip
import json
import re
from dataclasses import dataclass
from pathlib import Path

# Ensembl's breakend-to-feature distance (Utils/VariationEffect.pm:60, used by
# StructuralVariationOverlap.pm:143) and its up/downstream window, which is the same value.
MAX_DISTANCE_FROM_TRANSCRIPT = 5000

# The intron positions Ensembl's SV path reads as `intron_variant`: strictly inside the
# two-base splice sites. The polypyrimidine tract is the acceptor-side window 3 to 17
# bases into the intron (Utils/VariationEffect.pm `_intron_overlap`, the tract branch).
_SPLICE_SITE = 2
_PPT_INNER = 2
_PPT_OUTER = 16

_BRACKET_MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


@dataclass(frozen=True)
class TranscriptModel:
    """A transcript as the cache files it, with the genomic CDS bounds resolved."""

    chrom: str
    start: int
    end: int
    strand: int
    biotype: str
    exons: tuple[tuple[int, int], ...]
    cds: tuple[int, int] | None
    mirna_mature: tuple[tuple[int, int], ...]

    @property
    def introns(self) -> tuple[tuple[int, int], ...]:
        out = []
        for (s1, e1), (s2, _e2) in zip(self.exons, self.exons[1:]):
            if s2 - 1 >= e1 + 1:
                out.append((e1 + 1, s2 - 1))
        return tuple(out)

    def near(self, pos: int) -> bool:
        return self.start - MAX_DISTANCE_FROM_TRANSCRIPT <= pos <= self.end + MAX_DISTANCE_FROM_TRANSCRIPT

    def inside(self, pos: int) -> bool:
        return self.start <= pos <= self.end


def _int(v) -> int | None:
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _genomic_cds(tx: dict) -> tuple[int, int] | None:
    """Map the mapper's cDNA coding bounds to genomic coordinates through its pairs.

    Each pair maps a cDNA interval (`from_start`..`from_end`) onto a genomic interval
    (`to_start`..`to_end`) in the pair's orientation; a cDNA offset inside the pair lands
    at `to_start + offset` on the forward strand and `to_end - offset` on the reverse.
    The two mapped bounds are returned in ascending genomic order.
    """
    mapper = (tx.get("variation_effect_feature_cache") or {}).get("mapper") or {}
    cs, ce = _int(mapper.get("cdna_coding_start")), _int(mapper.get("cdna_coding_end"))
    if not cs or not ce or cs <= 0 or ce <= 0:
        return None
    gs = ge = None
    for pair in mapper.get("pairs") or []:
        fs, fe = _int(pair.get("from_start")), _int(pair.get("from_end"))
        ts, te, ori = _int(pair.get("to_start")), _int(pair.get("to_end")), _int(pair.get("ori"))
        if None in (fs, fe, ts, te, ori):
            continue
        if fs <= cs <= fe:
            gs = ts + (cs - fs) if ori >= 0 else te - (cs - fs)
        if fs <= ce <= fe:
            ge = ts + (ce - fs) if ori >= 0 else te - (ce - fs)
    if gs is None or ge is None:
        return None
    return (gs, ge) if gs <= ge else (ge, gs)


def _mirna_mature(tx: dict) -> tuple[tuple[int, int], ...]:
    out = []
    for attr in tx.get("attributes") or []:
        if attr.get("code") != "miRNA":
            continue
        m = re.fullmatch(r"\s*(\d+)-(\d+)\s*", str(attr.get("value", "")))
        if m:
            out.append((int(m.group(1)), int(m.group(2))))
    return tuple(out)


def model_from_cache_record(chrom: str, tx: dict) -> TranscriptModel | None:
    """Build a `TranscriptModel` from one transcript record of a cache shard."""
    start, end, strand = _int(tx.get("start")), _int(tx.get("end")), _int(tx.get("strand"))
    if start is None or end is None or strand is None:
        return None
    exons = []
    for e in tx.get("exons") or []:
        s, t = _int(e.get("start")), _int(e.get("end"))
        if s is not None and t is not None:
            exons.append((s, t))
    return TranscriptModel(
        chrom=chrom,
        start=start,
        end=end,
        strand=strand,
        biotype=str(tx.get("biotype") or ""),
        exons=tuple(sorted(exons)),
        cds=_genomic_cds(tx),
        mirna_mature=_mirna_mature(tx),
    )


def iter_cache_transcripts(cache_dir: str | Path):
    """Yield `(chromosome, transcript record)` for every transcript of a vep-rs JSON cache.

    Storable-derived caches file a shard as `{"chr": [tx, null, ...]}`, native ones as a
    list; a shard may be plain or gzipped (`.json.gz`, the golden corpora's form); an
    unreadable shard is skipped.
    """
    root = Path(cache_dir) / "transcripts"
    if not root.is_dir():
        return
    for chr_dir in sorted(root.iterdir()):
        if not chr_dir.is_dir():
            continue
        for shard in sorted(list(chr_dir.glob("*.json")) + list(chr_dir.glob("*.json.gz"))):
            try:
                if shard.suffix == ".gz":
                    with gzip.open(shard, "rt", encoding="utf-8") as fh:
                        data = json.load(fh)
                else:
                    data = json.loads(shard.read_text(encoding="utf-8"))
            except (OSError, ValueError, EOFError):
                continue
            if isinstance(data, dict):
                data = [x for v in data.values() if v for x in v if x]
            if not isinstance(data, list):
                continue
            for tx in data:
                if isinstance(tx, dict) and tx.get("stable_id"):
                    yield chr_dir.name, tx


def load_transcript_models(cache_dir: str | Path) -> dict[tuple[str, str], TranscriptModel]:
    """Every transcript of the cache by `(chromosome, stable id)`.

    The chromosome is part of the key because the pseudoautosomal transcripts are filed
    under both X and Y with their own coordinates on each. Empty when the
    cache is absent.
    """
    models: dict[tuple[str, str], TranscriptModel] = {}
    for chrom, tx in iter_cache_transcripts(cache_dir):
        model = model_from_cache_record(chrom, tx)
        if model is not None:
            models[(chrom, tx["stable_id"])] = model
    return models


def bracket_mate(allele: str) -> tuple[str, int] | None:
    """The mate `(chromosome, position)` a bracket-notation breakend allele names."""
    m = _BRACKET_MATE_RE.search(allele)
    if not m:
        return None
    chrom = m.group(1)
    return (chrom[3:] if chrom.lower().startswith("chr") else chrom), int(m.group(2))


def point_location(loc: str) -> tuple[str, int] | None:
    """`(chromosome, position)` of a point Location; None for a ranged one."""
    if ":" not in loc:
        return None
    chrom, _, coords = loc.partition(":")
    if "-" in coords or not coords.isdigit():
        return None
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return chrom, int(coords)


def _mature_mirna_hit(pos: int, model: TranscriptModel) -> str | None:
    """'hit' or 'miss' for a single-exon miRNA with mature ranges; 'multi-exon' when the
    cDNA position cannot be read off the span; None for any other transcript."""
    if model.biotype != "miRNA" or not model.mirna_mature:
        return None
    if len(model.exons) != 1:
        return "multi-exon"
    cdna = pos - model.start + 1 if model.strand == 1 else model.end - pos + 1
    return "hit" if any(lo <= cdna <= hi for lo, hi in model.mirna_mature) else "miss"


def region_terms(pos: int, model: TranscriptModel) -> set[str]:
    """Terms a point breakend at `pos`, inside the transcript, carries besides
    `feature_truncation`: the exon terms (CDS, 5' or 3' UTR by strand, the mature-miRNA
    term alone, or the non-coding exon term), or the intron terms (`intron_variant`
    strictly inside the splice sites, the polypyrimidine tract window, and the
    non-coding transcript term for a transcript without a CDS that is not NMD), plus
    `NMD_transcript_variant` on a nonsense-mediated-decay transcript."""
    terms: set[str] = set()
    for s, e in model.exons:
        if s <= pos <= e:
            if model.cds:
                cs, ce = model.cds
                if cs <= pos <= ce:
                    terms.add("coding_sequence_variant")
                elif (model.strand == 1 and pos < cs) or (model.strand == -1 and pos > ce):
                    terms.add("5_prime_UTR_variant")
                else:
                    terms.add("3_prime_UTR_variant")
            else:
                if _mature_mirna_hit(pos, model) == "hit":
                    return {"mature_miRNA_variant"}
                terms.add("non_coding_transcript_exon_variant")
            break
    else:
        for s, e in model.introns:
            if s <= pos <= e:
                if s + _SPLICE_SITE <= pos <= e - _SPLICE_SITE:
                    terms.add("intron_variant")
                if (model.strand == 1 and e - _PPT_OUTER <= pos <= e - _PPT_INNER) or (
                    model.strand == -1 and s + _PPT_INNER <= pos <= s + _PPT_OUTER
                ):
                    terms.add("splice_polypyrimidine_tract_variant")
                if not model.cds and model.biotype != "nonsense_mediated_decay":
                    terms.add("non_coding_transcript_variant")
                break
    if model.biotype == "nonsense_mediated_decay":
        terms.add("NMD_transcript_variant")
    return terms


def breakend_context(pos: int, model: TranscriptModel, mate_inside: bool) -> set[str]:
    """The term set a point breakend evaluated at `pos` supports on the transcript.

    `mate_inside` says whether the MATE breakend lies inside the transcript, which is
    what decides `feature_truncation` (Utils/VariationEffect.pm:358 reads the breakend,
    not the evaluated position). Inside the transcript the region terms apply; within
    `MAX_DISTANCE_FROM_TRANSCRIPT` outside it the strand-aware upstream or downstream
    term; beyond that nothing but the truncation term. A mature-miRNA hit stands alone.
    """
    if model.inside(pos):
        terms = region_terms(pos, model)
        if mate_inside and terms != {"mature_miRNA_variant"}:
            terms.add("feature_truncation")
        return terms
    if model.start - MAX_DISTANCE_FROM_TRANSCRIPT <= pos < model.start:
        terms = {"upstream_gene_variant" if model.strand == 1 else "downstream_gene_variant"}
    elif model.end < pos <= model.end + MAX_DISTANCE_FROM_TRANSCRIPT:
        terms = {"downstream_gene_variant" if model.strand == 1 else "upstream_gene_variant"}
    else:
        terms = set()
    if mate_inside:
        terms.add("feature_truncation")
    return terms


def mate_side_derivation(
    loc: str, allele: str, model: TranscriptModel
) -> tuple[frozenset[str], frozenset[str]] | None:
    """`(mate_set, local_set)` for a mate-side row, or None when the row is not one.

    The row is a mate-side row when the Location is a point, the Allele is bracket
    notation, the transcript is on the mate chromosome within the distance gate of the
    mate, and the local coordinate does not itself select the transcript. `local_set`
    is Ensembl's read with its empty result replaced by `intergenic_variant`, the
    default a Transcript row carries when no predicate fires.
    """
    local = point_location(loc)
    mate = bracket_mate(allele)
    if local is None or mate is None:
        return None
    lchr, lpos = local
    mchr, mpos = mate
    if model.chrom != mchr or not model.near(mpos):
        return None
    if lchr == model.chrom and model.near(lpos):
        return None
    mate_inside = model.inside(mpos)
    mate_set = breakend_context(mpos, model, mate_inside)
    local_set = breakend_context(lpos, model, mate_inside) or {"intergenic_variant"}
    return frozenset(mate_set), frozenset(local_set)
