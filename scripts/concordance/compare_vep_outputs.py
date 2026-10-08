#!/usr/bin/env python3
"""Semantic concordance comparison for Perl VEP vs vep-rs outputs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path
from typing import Callable, Iterable, Iterator, TextIO

TupleKey = tuple[str, str, str, str, str]
ContextTermSet = frozenset[str]

TRANSCRIPT_CONTEXT_TERMS: ContextTermSet = frozenset(
    {
        "non_coding_transcript_variant",
        "NMD_transcript_variant",
        "coding_transcript_variant",
    }
)

# Splice terms that Rust may add but Perl drops via _get_differing_regions
# last-write-wins architecture. This is an intended divergence where vep-rs is
# scientifically correct and Perl's behavior is architecturally flawed.
SPLICE_TERMS: frozenset[str] = frozenset(
    {
        "splice_region_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
    }
)


# Adjusted F1: intended-divergence filters for SNP/indel suites


def _is_splice_overcall_swap(perl_csq: str, rust_csq: str) -> bool:
    """Detect consequence swaps caused by Perl's splice last-write-wins.

    Perl's ``_intron_effects`` iterates differing regions produced by
    ``_get_differing_regions`` and *assigns* (not ORs) the splice result for
    each region.  The last region's result overwrites all prior regions.
    When the final differing segment of a complex indel does not overlap
    the exonic splice window, Perl drops splice_region / splice sub-terms
    that earlier segments correctly detected.

    Rust evaluates ALL differing regions and retains the union of splice
    terms, which is scientifically more complete.

    This function returns True when the swap is ONLY in splice-related
    terms (with optional ``intron_variant``), indicating Rust adds splice
    information that Perl's architecture drops.

    Three sub-patterns:
      1. Rust adds splice_region_variant.
      2. Rust adds specific sub-terms like splice_donor, while Perl has
         the weaker splice_region or nothing.
      3. Rust adds intron_variant plus splice sub-terms for large indels
         spanning into introns.
    """
    pset = set(perl_csq.split(","))
    rset = set(rust_csq.split(","))
    diff_p = pset - rset  # in Perl but not Rust
    diff_r = rset - pset  # in Rust but not Perl

    if not diff_r:
        return False  # Rust must be adding something

    splice_and_intron = SPLICE_TERMS | {"intron_variant"}

    # All terms Rust adds must be splice-related (or intron_variant)
    if not diff_r <= splice_and_intron:
        return False

    # STRICT SUBSET REQUIRED. Perl must have NO term Rust lacks.
    #
    # A non-empty diff_p would admit SUBSTITUTIONS, not just term-loss, and the
    # cited Perl mechanism cannot produce them: of the fourteen writes in
    # `_intron_effects`, only `splice_region` is a plain assignment. Every other
    # key, including those behind splice_donor_variant / splice_acceptor_variant,
    # is a guarded write-1 that no later differing region can clear. So Perl
    # cannot lose those terms this way, and a pair where Perl carries one and
    # Rust does not is a Rust downgrade being excluded as a Perl defect.
    if diff_p:
        return False

    # Rust must add at least one splice term (not just intron_variant alone)
    if not (diff_r & SPLICE_TERMS):
        return False

    return True


def _is_ppt_overcall_swap(perl_csq: str, rust_csq: str) -> bool:
    """Detect intronic PPT/splice_region overcall swaps (the four pure-intronic shapes).

    Pattern: both Perl and Rust share ``intron_variant`` as the base term, and
    Rust adds one or both of ``splice_polypyrimidine_tract_variant`` /
    ``splice_region_variant`` to an intronic variant where Perl does not.
    Perl's last-write-wins architecture or stricter sequence-based PPT window
    drops these sub-terms for short intronic indels.

    Rust emits the union of all applicable splice sub-terms.  The resulting
    annotation is scientifically more complete: the variant genuinely sits
    inside the polypyrimidine tract or the 7-bp splice_region flanking
    window.

    This is a *more specific* sibling of :func:`_is_splice_overcall_swap`.
    It must be checked BEFORE the broader filter so that the pure-intronic
    pairs are tagged in their own category bucket for reporting.

    Returns True iff all of:
      * ``intron_variant`` is in both Perl and Rust sets,
      * Perl has no terms Rust lacks (``pset - rset`` is empty),
      * Rust's extra terms are a non-empty subset of
        ``{splice_polypyrimidine_tract_variant, splice_region_variant}``.
    """
    pset = set(perl_csq.split(","))
    rset = set(rust_csq.split(","))

    if "intron_variant" not in pset or "intron_variant" not in rset:
        return False

    diff_p = pset - rset  # Perl-only terms
    diff_r = rset - pset  # Rust-only terms

    if diff_p:
        return False  # Perl must not have any terms Rust dropped

    if not diff_r:
        return False  # Rust must add at least one term

    ppt_and_region = {
        "splice_polypyrimidine_tract_variant",
        "splice_region_variant",
    }
    return diff_r <= ppt_and_region


_START_TERMS: frozenset[str] = frozenset({"start_lost", "start_retained_variant"})


def is_structural_allele(allele: str) -> bool:
    """Whether a tuple's allele column names a structural allele.

    VEP writes a sequence allele as its bases or ``-`` and a symbolic structural
    allele as its Sequence Ontology class (``deletion``, ``duplication``,
    ``tandem_repeat`` and so on), so anything outside ``-`` and the nucleotide
    alphabet is structural.
    """
    return allele != "-" and not set(allele.upper()) <= set("ACGTN")


def _is_start_cooccurrence_swap(perl_csq: str, rust_csq: str, allele: str) -> bool:
    """Detect the contradictory co-emission of start_lost + start_retained_variant.

    Perl emits BOTH ``start_lost`` AND ``start_retained_variant`` for one
    allele and transcript on two routes, and which member is erroneous depends
    on the allele kind:

      * On a sequence variant, ``start_retained_variant`` fires when
        ``_ins_del_start_altered`` (VariationEffect.pm:1028-1075) returns FALSE,
        which it does only when the edited 5'UTR and CDS still begin ``ATG`` at
        the CDS start or the CDS is unchanged (:1073); ``start_lost`` then fires
        anyway, on most pairs through ``_inv_start_altered`` (:906-944), which
        reads the three bases at the OLD 5'UTR length after the edit (:940) and
        so misses a codon the edit re-forms one base over, and on the rest
        through the peptide route (:879-882), which compares the codon-window
        peptides the edit touches. On either route the edited CDS begins ``ATG``,
        re-formed byte-identically one base over (``_inv_start_altered``) or
        standing at the old offset with the downstream frame shifted (the
        peptide route). The start codon is intact, so
        ``start_lost`` is the erroneous member and vep-rs emits the set without it.
      * On a structural allele, ``start_lost`` fires from the overlap of the
        span with the start codon (:884-893) and ``_ins_del_start_altered``
        returns 0 without reading the sequence (:1037), so
        ``start_retained_variant`` (:961) fires on every structural allele over
        a start codon, a deletion that removes the codon included, and is the
        erroneous member; vep-rs emits the set without it.

    Returns True when Perl has both start terms and vep-rs has the same
    consequence set minus the erroneous member for the allele kind.
    """
    pset = set(perl_csq.split(","))
    rset = set(rust_csq.split(","))

    if not _START_TERMS <= pset:
        return False
    # The ONLY difference is the one start term, in BOTH directions: `pset - rset`
    # alone would leave the Rust side unconstrained.
    if rset - pset:
        return False
    dropped = "start_retained_variant" if is_structural_allele(allele) else "start_lost"
    return pset - rset == {dropped}


# No looser start co-emission mask exists. A predicate that admits any pair where
# Perl carries a start term vep-rs lacks absorbs, besides the strict shape, pairs
# where vep-rs drops the OTHER start term for its allele kind (a vep-rs defect, not
# the co-emission), pairs where vep-rs misses start_lost outright
# (P=5UTR,start_lost,start_retained / R=5UTR; P=frameshift,start_lost /
# R=frameshift), unrelated protein-altering-vs-inframe shapes
# (P=5UTR,coding_sequence_variant / R=5UTR,start_lost), and the inverted shape
# (P=frameshift,start_retained / R=frameshift,start_lost) where Perl emits ONLY
# start_retained, a separate predicate-overlap shape. Leaving those unmasked keeps
# the vep-rs divergences visible and countable. The strict
# `_is_start_cooccurrence_swap` is the only start filter, cited to
# VariationEffect.pm:851-881, :947-963 and :1028-1075.


# Coding terms that must be PRESERVED on both sides for the frameshift-splice
# cascade filter. One of these on Perl but not vep-rs (or vice versa) is a term
# DOWNGRADE (frameshift -> coding_sequence), a real divergence that must NOT be
# filtered.
_FRAMESHIFT_CASCADE_CODING_TERMS: frozenset[str] = frozenset(
    {
        "frameshift_variant",
        "inframe_insertion",
        "inframe_deletion",
        "stop_gained",
        "stop_lost",
        "start_lost",
    }
)


def _is_frameshift_splice_cascade_swap(perl_csq: str, rust_csq: str) -> bool:
    """Detect splice-cascade overcalls where coding terms are preserved.

    Splice-cascade overcall shapes:
      * Perl ``frameshift_variant,stop_gained`` ->
            vep-rs ``frameshift_variant,intron_variant,splice_acceptor_variant,stop_gained``
      * Perl ``frameshift_variant`` ->
            vep-rs ``frameshift_variant,intron_variant,splice_donor_5th_base_variant,splice_donor_variant``
      * Perl ``frameshift_variant,stop_gained`` ->
            vep-rs ``frameshift_variant,intron_variant,splice_donor_5th_base_variant,splice_donor_variant,stop_gained``

    Distinct from ``_is_splice_overcall_swap`` because the variant retains one or
    more strong coding terms (frameshift, inframe, stop_gained, stop_lost,
    start_lost) on BOTH sides; vep-rs only adds intronic/splice sub-terms, which
    Perl drops via the ``_intron_effects`` last-write-wins architecture.

    Negative case: do NOT match when vep-rs has DOWNGRADED a coding term (Perl
    ``frameshift_variant`` / vep-rs ``coding_sequence_variant,...``), a real
    divergence; the full set of strong coding terms must be identical on both
    sides.

    Returns True when:
      1. Both sides carry at least one term from
         ``_FRAMESHIFT_CASCADE_CODING_TERMS``.
      2. The strong coding terms are identical on both sides (no downgrade).
      3. All Rust-only extra terms are in ``SPLICE_TERMS | {intron_variant}``.
      4. Rust adds at least one such term (so this really is a swap).
      5. Perl-only extras (if any) are also splice/intron terms (e.g. Perl
         ``splice_donor_variant`` becomes Rust
         ``intron + splice_donor_5th_base + splice_donor``).
    """
    pset = set(perl_csq.split(","))
    rset = set(rust_csq.split(","))

    # (1) + (2): coding backbone must be present and identical on both sides.
    pcoding = pset & _FRAMESHIFT_CASCADE_CODING_TERMS
    rcoding = rset & _FRAMESHIFT_CASCADE_CODING_TERMS
    if not pcoding or not rcoding:
        return False
    if pcoding != rcoding:
        return False

    diff_r = rset - pset
    diff_p = pset - rset
    if not diff_r:
        return False  # Rust must be adding something
    splice_and_intron = SPLICE_TERMS | {"intron_variant"}
    if not diff_r <= splice_and_intron:
        return False

    # STRICT SUBSET REQUIRED, for the reason the sibling predicate documents: the
    # cited mechanism reaches only `splice_region`, so a pair where Perl carries an
    # acceptor/donor/fifth-base term vep-rs lacks is a vep-rs downgrade rather than
    # a Perl term-loss. The cascade shapes above are all pure additions and match.
    if diff_p:
        return False

    # vep-rs must add at least one SPLICE term, not `intron_variant` alone: no
    # splice mechanism explains a pair differing only by `intron_variant`.
    if not (diff_r & SPLICE_TERMS):
        return False

    return True


# Set of terms allowed in an "all-intronic" swap.  Any term outside this set
# (coding, UTR, upstream/downstream) disqualifies the pattern.
_INTRONIC_ALLOWED_TERMS: frozenset[str] = frozenset(
    {
        "intron_variant",
        "splice_region_variant",
        "splice_polypyrimidine_tract_variant",
    }
)


def _is_intronic_overcall_swap(perl_csq: str, rust_csq: str) -> bool:
    """Detect pure-intronic PPT / splice_region overcall swaps.

    Pure-intronic overcall swap patterns:
      * (a) Perl ``intron_variant,splice_polypyrimidine_tract_variant`` ->
            Rust ``intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant``
      * (b) Perl ``intron_variant,splice_region_variant`` ->
            Rust ``intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant``
      * (c) Perl ``intron_variant`` ->
            Rust ``intron_variant,splice_polypyrimidine_tract_variant``
      * (d) Perl ``intron_variant`` ->
            Rust ``intron_variant,splice_region_variant``

    Rust's intronic insertion/non-insertion splice logic adds PPT and/or
    splice_region terms for positions inside the PPT window (2-16 bp from
    acceptor) or the 7 bp splice_region flanking window.  Perl's predicates
    are stricter (possibly sequence-based), so Rust overcalls.

    This helper is MORE SPECIFIC than ``_is_splice_overcall_swap`` and
    ``_is_ppt_overcall_swap``: both sides must contain ONLY terms from
    ``_INTRONIC_ALLOWED_TERMS`` (no coding, no UTR, no upstream/downstream).
    Because it is a strict subset of the splice-overcall surface, it must be
    evaluated BEFORE the generic splice helpers so these pairs are routed into
    the dedicated ``intronic_overcall_swap`` bucket.

    Returns True when:
      1. Perl's set is a subset of ``_INTRONIC_ALLOWED_TERMS``.
      2. Rust's set is a subset of ``_INTRONIC_ALLOWED_TERMS``.
      3. Perl has ``intron_variant`` (anchors the pattern as intronic).
      4. Rust adds at least one term that Perl lacks.
      5. Perl has no terms that Rust lacks (Rust is a superset).
    """
    pset = set(perl_csq.split(","))
    rset = set(rust_csq.split(","))

    if not pset <= _INTRONIC_ALLOWED_TERMS:
        return False
    if not rset <= _INTRONIC_ALLOWED_TERMS:
        return False
    if "intron_variant" not in pset:
        return False

    diff_r = rset - pset
    diff_p = pset - rset
    if not diff_r:
        return False
    if diff_p:
        return False

    return True


def _is_covered_splice_region_swap(perl_csq: str, rust_csq: str) -> bool:
    """The ONE divergence shape Ensembl's ``_intron_effects`` can actually produce.

    ``sub _intron_effects`` (BaseTranscriptVariationAllele.pm:99-224 in
    ensembl-variation release/115 (23c76f60)) performs fourteen writes into
    ``$intron_effects``. Thirteen of them, at lines 139, 149, 156, 160, 172, 177,
    181, 185, 189, 193, 197, 201 and 205, are guarded ``= 1`` assignments that
    nothing clears. Only line 215 assigns a RETURN VALUE, ``_intron_overlap(...)``,
    to ``splice_region``, and that function returns a defined ``0``
    (Utils/VariationEffect.pm:107, same release), so a later boundary intron can
    clobber an earlier ``1``. The ``unless`` at line 216 then freezes
    ``splice_region`` once any region has set ``start_splice_site`` or
    ``end_splice_site``.

    The overwrite can therefore remove exactly one consequence term,
    ``splice_region_variant``. The flags behind
    ``splice_polypyrimidine_tract_variant`` (156, 160, 193, 205),
    ``splice_acceptor_variant`` and ``splice_donor_variant`` (177, 181),
    ``splice_donor_5th_base_variant`` (185, 197),
    ``splice_donor_region_variant`` (189, 201) and ``intron_variant`` (149) are
    all set-once, so Ensembl VEP cannot lose any of them to last-write-wins, and
    a pair where it lacks one is not explained by this mechanism.

    Hence the covered shape, the only shape the adjusted F1 removes from the
    splice family:

      * the Rust-side addition contains ``splice_region_variant``,
      * the Rust-side addition is a subset of
        ``{splice_region_variant, intron_variant}``, and
      * the Perl-side difference is empty.

    ``intron_variant`` is admitted alongside it because line 149 sets the
    intronic flag inside the same per-region loop that the ``next`` at lines 140
    and 173 can skip, so a region Ensembl VEP skipped withholds both terms
    together.

    This is a property of the PAIR and not of any reporting bucket. It cuts
    across all four splice buckets in both directions: of the patterns
    ``_is_intronic_overcall_swap`` documents, (a) and (d) add only
    ``splice_region_variant`` and are covered, while (b) and (c) add the
    polypyrimidine term and are not. That is why exclusion is evaluated
    independently of the classification chain rather than per bucket.
    """
    pset = set(perl_csq.split(","))
    rset = set(rust_csq.split(","))

    # Perl must have no term Rust lacks: a Perl-side loss is a vep-rs downgrade,
    # for the reason `_is_splice_overcall_swap` documents at its own guard.
    if pset - rset:
        return False

    diff_r = rset - pset
    if "splice_region_variant" not in diff_r:
        return False
    return diff_r <= {"splice_region_variant", "intron_variant"}


# Transcript selection on a record Perl marks `vep_skip`: the structural-variant
# comparator's class, applied here because a whole-genome input carries such records.
# Perl's --max_sv_size default is 10,000,000 (ensembl-vep Config.pm:310); a span above it
# is marked skipped (Parser.pm:493-497), as is a symbolic ALT type without a Sequence
# Ontology term (Parser/VCF.pm:477-481, `_VEP_SKIP_UNSUPPORTED_ALLELES`). A skipped
# record loads no cache region of its own (AnnotationSource.pm:238) and is annotated
# against whatever regions its batch loaded, so the transcripts Perl names for it are a
# property of the neighbouring records, while vep-rs names every overlapping transcript.
# `<NON_REF>` takes the same route but is not listed: vep-rs writes no transcript tuple
# for a gVCF reference block, so the arms below have nothing to reach there.
_MAX_SV_SIZE = 10_000_000
_VEP_SKIP_UNSUPPORTED_ALLELES = frozenset({"CPX", "CTX"})
SV_TRANSCRIPT_SELECTION_BUCKET = "sv_transcript_selection"


def _location_span(loc: str) -> tuple[str, int, int] | None:
    """``chrom:start-end`` or ``chrom:pos`` as (chrom, start, end); None when unparseable."""
    m = re.fullmatch(r"([^:]+):(\d+)(?:-(\d+))?", loc)
    if not m:
        return None
    start = int(m.group(2))
    return m.group(1), start, int(m.group(3)) if m.group(3) else start


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


def _batch_dependent_annotation_scope(loc: str, allele: str) -> bool:
    """True where Perl annotates a record only against its batch's loaded regions.

    Two routes reach that state, both through the same `vep_skip` flag: a span above
    --max_sv_size (Parser.pm:493-497) and a symbolic type without a Sequence Ontology
    term (Parser/VCF.pm:477-481, `_VEP_SKIP_UNSUPPORTED_ALLELES`). Below the size limit
    and with a supported type Perl loads every region the record overlaps
    (AnnotationSource.pm:194-200), so a transcript-set difference there is not this
    defect and stays charged.
    """
    return _span_exceeds_max_sv_size(loc) or allele in _VEP_SKIP_UNSUPPORTED_ALLELES


# The features each engine named on every record in `_batch_dependent_annotation_scope`,
# keyed on (location, allele) and collected during the sorted merge, where both engines'
# complete tuple streams are read once.
ScopeFeatures = dict[tuple[str, str], tuple[set[str], set[str]]]


def _record_scope_feature(scope: ScopeFeatures, key_line: str, side: int) -> None:
    """Add the feature of a key line to its record's set for one side (0 Perl, 1 vep-rs)."""
    loc, allele, feature = key_line.split("\t", 3)[:3]
    if _batch_dependent_annotation_scope(loc, allele):
        scope.setdefault((loc, allele), (set(), set()))[side].add(feature)


@dataclass(frozen=True)
class ExclusionRule:
    """One divergence shape: its recogniser, its reporting name, whether the
    adjusted F1 removes it, and the taxonomy class it belongs to.

    ``classifies`` and ``excludes`` are independent because the two roles are
    independent. A rule that classifies assigns a pair to its reporting bucket;
    a rule that excludes removes the pair from both raw denominators. Only
    ``start_cooccurrence_swap`` does both.
    """

    bucket: str
    predicate: Callable[..., bool]
    classifies: bool
    excludes: bool
    taxonomy_class: str
    # A predicate that decides on the allele kind as well as the two consequence
    # sets takes the tuple's allele as its third argument; the others take two.
    allele_aware: bool = False

    def matches(self, perl_csq: str, rust_csq: str, allele: str) -> bool:
        if self.allele_aware:
            return self.predicate(perl_csq, rust_csq, allele)
        return self.predicate(perl_csq, rust_csq)


# The single source of truth for what the adjusted F1 removes and what it merely
# reports, so adding or retiring a rule moves the expected partition here instead
# of in a hardcoded tally kept in step by hand.
#
# ORDER IS THE CLASSIFICATION ORDER, most specific first, because classification
# is first-match:
#   * `intronic_overcall_swap` requires pure-intronic sets on BOTH sides and so is
#     narrower than `ppt_overcall_swap`.
#   * `ppt_overcall_swap` requires the PPT + intron signature and so is narrower
#     than `splice_lastwrite_swap`.
#   * `frameshift_splice_cascade_swap` requires an identical strong-coding backbone
#     on both sides, a signature the generic splice rule does not enforce.
#   * `splice_lastwrite_swap` is the general splice fallback.
#   * `start_cooccurrence_swap` is the strict start co-emission rule.
#
# Every rule except `covered_splice_region_swap` classifies; that one is the mask
# POLICY predicate, which is a property of the pair rather than a bucket, so it
# excludes without owning a bucket of its own.
EXCLUSION_REGISTRY: tuple[ExclusionRule, ...] = (
    ExclusionRule(
        bucket="intronic_overcall_swap",
        predicate=_is_intronic_overcall_swap,
        classifies=True,
        excludes=False,
        taxonomy_class="splice_family_swap",
    ),
    ExclusionRule(
        bucket="ppt_overcall_swap",
        predicate=_is_ppt_overcall_swap,
        classifies=True,
        excludes=False,
        taxonomy_class="splice_family_swap",
    ),
    ExclusionRule(
        bucket="frameshift_splice_cascade_swap",
        predicate=_is_frameshift_splice_cascade_swap,
        classifies=True,
        excludes=False,
        taxonomy_class="splice_family_swap",
    ),
    ExclusionRule(
        bucket="splice_lastwrite_swap",
        predicate=_is_splice_overcall_swap,
        classifies=True,
        excludes=False,
        taxonomy_class="splice_family_swap",
    ),
    ExclusionRule(
        bucket="start_cooccurrence_swap",
        predicate=_is_start_cooccurrence_swap,
        classifies=True,
        excludes=True,
        taxonomy_class="start_cooccurrence_swap",
        allele_aware=True,
    ),
    ExclusionRule(
        bucket="covered_splice_region_swap",
        predicate=_is_covered_splice_region_swap,
        classifies=False,
        excludes=True,
        taxonomy_class="splice_family_swap",
    ),
)

CLASSIFICATION_CHAIN: tuple[ExclusionRule, ...] = tuple(
    r for r in EXCLUSION_REGISTRY if r.classifies
)
EXCLUDING_RULES: tuple[ExclusionRule, ...] = tuple(
    r for r in EXCLUSION_REGISTRY if r.excludes
)


# One sentence per pair bucket, the shape its predicate recognises, for the report.
PAIR_RULE_DEFINITIONS: dict[str, str] = {
    "intronic_overcall_swap": (
        "Both sides are pure-intronic sets and vep-rs adds splice_region_variant or "
        "splice_polypyrimidine_tract_variant, which the reference's last-written region drops."
    ),
    "ppt_overcall_swap": (
        "Both sides carry intron_variant and vep-rs adds splice_polypyrimidine_tract_variant "
        "or splice_region_variant to the same set."
    ),
    "frameshift_splice_cascade_swap": (
        "An identical frameshift or in-frame coding backbone on both sides, vep-rs adding "
        "splice sub-terms the reference's last-written region drops."
    ),
    "splice_lastwrite_swap": (
        "vep-rs adds splice terms, and nothing else, to the reference's set: the reference "
        "stores the splice verdict by assignment and a later differing region overwrote it."
    ),
    "start_cooccurrence_swap": (
        "The reference writes start_lost and start_retained_variant together for one allele; "
        "vep-rs writes the set without the member the allele kind cannot support."
    ),
    "covered_splice_region_swap": (
        "The pairs of the splice family whose only difference is a splice_region_variant "
        "the reference's assignment order dropped; the exclusion policy of that family."
    ),
}
PAIR_RULE_CITATIONS: dict[str, str] = {
    "splice_family_swap": "BaseTranscriptVariationAllele.pm:215; VariationEffect.pm:87-110",
    "start_cooccurrence_swap": "VariationEffect.pm:851-963, 1028-1075",
}


_BRACKET_MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
# The record's own breakend as the two engines write it: `N.` / `.N` (the reference's
# string for the local end of a bracket breakend, and a single breakend's own VCF ALT)
# or `chromosome_breakpoint` (the class term, written for a mateless single breakend by
# the reference and for the local end of a bracket breakend by vep-rs).
_OWN_BREAKEND_RE = re.compile(r"^(?:\.[ACGTN]+|[ACGTN]+\.|chromosome_breakpoint)$")
# The one symbolic record both engines keep but write at different Locations: the
# reference expands a `<CNV:TR>` to literal sequence at the run's 3' end (Parser/VCF.pm:
# 382-431) while vep-rs writes the span under this allele, so its vep-rs rows never have a
# reference row at their Location and are not a skipped record.
_TANDEM_REPEAT_ALLELE = "tandem_repeat"
_INTERGENIC_FEATURE = "-"


def bracket_mate(allele: str) -> tuple[str, int] | None:
    """The mate ``(chromosome, position)`` a bracket-notation breakend allele names, the
    chromosome with any ``chr`` prefix removed as the reference removes it
    (AnnotationType/Transcript.pm:134); None for any other allele."""
    m = _BRACKET_MATE_RE.search(allele)
    if not m:
        return None
    chrom = m.group(1)
    return (chrom[3:] if chrom.lower().startswith("chr") else chrom), int(m.group(2))


def is_own_breakend_allele(allele: str) -> bool:
    """Whether the allele column names a breakend record's own end rather than a mate."""
    return _OWN_BREAKEND_RE.match(allele) is not None


@dataclass(frozen=True)
class OneSidedRule:
    """A divergence shape recognised on rows one engine wrote and the other did not.

    Unlike an `ExclusionRule`, which judges a pair of consequence sets on one
    transcript, a one-sided rule judges a record: what the reference wrote for it
    (nothing, or rows under some alleles) against what vep-rs wrote. ``excludes`` says
    whether this comparator's adjusted F1 removes the rows it recognises;
    ``excludes_with_cache`` whether a comparator with a transcript-to-chromosome
    authority (the structural-variant comparator under its cache) does. A rule that only
    classifies reports a count so the shape is visible without being credited to vep-rs.
    """

    bucket: str
    definition: str
    excludes: bool
    taxonomy_class: str
    perl_citation: str
    excludes_with_cache: bool = False


REFERENCE_SKIPPED_RECORD_BUCKET = "reference_skipped_record"
BND_OWN_ROWS_LOST_BUCKET = "bnd_own_chromosome_rows_lost"
BND_SYNONYM_MATE_ROWS_LOST_BUCKET = "bnd_synonym_mate_rows_lost"
FASTA_NAMED_SLICE_INTERGENIC_BUCKET = "fasta_named_slice_intergenic"

# The three breakend shapes are counted, never set aside: deciding which of the rows the
# reference would have written needs the cache's transcript-to-chromosome authority,
# which this comparator does not read (the structural-variant comparator sets the first
# two aside under it).
REFERENCE_116_2_ONE_SIDED_RULES: tuple[OneSidedRule, ...] = (
    OneSidedRule(
        bucket=REFERENCE_SKIPPED_RECORD_BUCKET,
        definition=(
            "A structural record the reference dropped before its input buffer, marked "
            "vep_skip for a span over --max_sv_size (its SVLEN, else its END), "
            "an ALT type with no Sequence Ontology term, an unsupported multi-allelic list or "
            "a deletion with neither END nor SVLEN, so it has no reference row in any format; "
            "vep-rs's rows for the record are set aside."
        ),
        excludes=True,
        excludes_with_cache=True,
        taxonomy_class="reference_skipped_record",
        perl_citation="Parser.pm:491-500, 821-826; Parser/VCF.pm:477-481, 555-557, 575",
    ),
    OneSidedRule(
        bucket=BND_OWN_ROWS_LOST_BUCKET,
        definition=(
            "A bracket breakend whose mate is on another chromosome loses the rows of its "
            "own end when a mate-chromosome transcript, in memory from the previous buffer, "
            "reaches it first and names its slice: the reference keeps the mate rows and "
            "writes nothing under the own-end allele, an order-dependent result."
        ),
        excludes=False,
        excludes_with_cache=True,
        taxonomy_class="bnd_slice_claimed_by_mate_chromosome",
        perl_citation="AnnotationSource.pm:119-134; AnnotationType/Transcript.pm:122-143",
    ),
    OneSidedRule(
        bucket=BND_SYNONYM_MATE_ROWS_LOST_BUCKET,
        definition=(
            "A bracket breakend whose mate is written with a chromosome synonym passes the "
            "overlap test through the synonym table but fails the slice test, which "
            "compares the raw name with the cache's, so the reference writes nothing under "
            "the mate allele."
        ),
        excludes=False,
        excludes_with_cache=True,
        taxonomy_class="bnd_synonym_mate_slice_mismatch",
        perl_citation="InputBuffer.pm:348-358; AnnotationType/Transcript.pm:134",
    ),
    OneSidedRule(
        bucket=FASTA_NAMED_SLICE_INTERGENIC_BUCKET,
        definition=(
            "Under --check_ref or --lookup_ref with a FASTA whose sequence names differ "
            "from the cache's, the variant's slice is named after the FASTA and every "
            "transcript fails the slice test, so the reference writes intergenic_variant "
            "where vep-rs writes transcript rows; a property of the run, counted and never "
            "set aside."
        ),
        excludes=False,
        taxonomy_class="fasta_named_slice_intergenic",
        perl_citation="Parser.pm:749; BaseVEP.pm get_slice; AnnotationType/Transcript.pm:141",
    ),
)


@dataclass(frozen=True)
class ReferenceRegistry:
    """Every rule the comparator applies against one Ensembl VEP release: the pair rules
    (classification and exclusion of consequence-set swaps) and the one-sided rules
    (record-level classes). The 115.2 registry's one-sided class is the scope pass of
    `filter_snp_indel_intended_divergences`, described here and implemented there, so the
    published figures reproduce from the code that produced them."""

    release: str
    pair_rules: tuple[ExclusionRule, ...]
    one_sided_rules: tuple[OneSidedRule, ...]


# Both registries share the pair rules: the splice assignment and the start co-emission
# they recognise are unchanged in 116 (VariationEffect.pm:851-1075 and
# BaseTranscriptVariationAllele.pm:215).
REFERENCE_REGISTRIES: dict[str, ReferenceRegistry] = {
    "115.2": ReferenceRegistry(
        release="115.2",
        pair_rules=EXCLUSION_REGISTRY,
        one_sided_rules=(
            OneSidedRule(
                bucket=SV_TRANSCRIPT_SELECTION_BUCKET,
                definition=(
                    "A record the reference marks vep_skip stays in its batch and is "
                    "annotated against whatever regions its neighbours loaded; the "
                    "transcripts vep-rs names beyond the reference's are set aside when "
                    "the reference named at least one and vep-rs more, and the reverse "
                    "when vep-rs named at least as many."
                ),
                excludes=True,
                excludes_with_cache=True,
                taxonomy_class="sv_transcript_selection",
                perl_citation="Parser.pm:491-498; Parser/VCF.pm:477-481, 575; AnnotationSource.pm:238; InputBuffer.pm:284-330",
            ),
        ),
    ),
    "116.2": ReferenceRegistry(
        release="116.2",
        pair_rules=EXCLUSION_REGISTRY,
        one_sided_rules=REFERENCE_116_2_ONE_SIDED_RULES,
    ),
}
REFERENCE_RELEASES: tuple[str, ...] = tuple(REFERENCE_REGISTRIES)
DEFAULT_REFERENCE_RELEASE = "116.2"


def divergence_class_definitions(reference_release: str) -> dict[str, dict[str, object]]:
    """Every class of a reference release's registry as the report states it: the
    one-sentence definition, whether the adjusted F1 sets it aside and the Ensembl VEP
    source it cites. Pair rules first, in classification order, then the one-sided rules."""
    registry = REFERENCE_REGISTRIES[reference_release]
    out: dict[str, dict[str, object]] = {}
    for rule in registry.pair_rules:
        out[rule.bucket] = {
            "definition": PAIR_RULE_DEFINITIONS[rule.bucket],
            "excludes": rule.excludes,
            "taxonomy_class": rule.taxonomy_class,
            "perl_citation": PAIR_RULE_CITATIONS[rule.taxonomy_class],
        }
    for one_sided in registry.one_sided_rules:
        out[one_sided.bucket] = {
            "definition": one_sided.definition,
            "excludes": one_sided.excludes,
            "taxonomy_class": one_sided.taxonomy_class,
            "perl_citation": one_sided.perl_citation,
        }
    return out


class ReferenceRecordScope:
    """What the 116.2 one-sided pass needs from the sorted merge, collected in the one
    read of both key streams.

    ``absent``: the (file, Location) pairs at which vep-rs wrote a structural row and
    the reference wrote nothing: a record the reference dropped as `vep_skip`, since
    116.2 writes at least one row (an intergenic one when nothing is in range) for every
    record it keeps and every route to `vep_skip` is structural; `tandem_repeat` is
    outside the test (`_TANDEM_REPEAT_ALLELE`). Every vep-rs row at an absent Location is
    then set aside whatever its allele, a sequence allele beside the symbolic one of a
    mixed ALT list included, because the reference skips the line whole. This is the
    output-only rule: this comparator reads no input VCF. The structural-variant
    comparator reads the input and names the record from it (`skipped_record_rows`:
    the reference's own skip routes over the ALT and INFO, confirmed by the reference's
    absence under the record's own Location start and alleles), so a kept record
    sharing a start or a Location string with a skipped one is told apart there. The
    merge visits keys in sorted order and every key of one Location shares the prefix
    ``Location<TAB>``, so a Location's rows from both streams are contiguous and the
    group closes when the Location changes.

    ``breakends``: for every (file, Location) at which either engine wrote a breakend
    allele, the alleles each side wrote and the mate chromosomes its bracket alleles
    name, for the two breakend classes.
    """

    def __init__(self) -> None:
        self.absent: set[tuple[str, str]] = set()
        self.breakends: dict[tuple[str, str], tuple[set[str], set[str], set[str]]] = {}
        self._group: tuple[str, str] | None = None
        self._reference_wrote = False
        self._rust_wrote_structural = False

    def observe(self, file_name: str, key_line: str, perl: bool, rust: bool) -> None:
        loc, allele = key_line.split("\t", 2)[:2]
        group = (file_name, loc)
        if group != self._group:
            self.close()
            self._group = group
        if perl:
            self._reference_wrote = True
        if rust and is_structural_allele(allele) and allele != _TANDEM_REPEAT_ALLELE:
            self._rust_wrote_structural = True
        mate = bracket_mate(allele)
        if mate is not None or is_own_breakend_allele(allele):
            perl_alleles, rust_alleles, mates = self.breakends.setdefault(
                group, (set(), set(), set())
            )
            if perl:
                perl_alleles.add(allele)
            if rust:
                rust_alleles.add(allele)
            if mate is not None:
                mates.add(mate[0])

    def close(self) -> None:
        """Close the open Location group; called by the merge when its file ends."""
        if self._group is not None and self._rust_wrote_structural and not self._reference_wrote:
            self.absent.add(self._group)
        self._group = None
        self._reference_wrote = False
        self._rust_wrote_structural = False


def _strip_chr(chrom: str) -> str:
    """The reference's own normalisation of a breakend's chromosome name at the slice test
    (AnnotationType/Transcript.pm:134): the ``chr`` prefix removed and nothing else, so
    ``M`` stays ``M`` and is not the cache's ``MT``."""
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def apply_116_2_one_sided_rules(
    rows: list[dict[str, str]],
    masked_rows: set[int],
    scope: ReferenceRecordScope,
) -> dict[str, dict[str, int]]:
    """The 116.2 one-sided pass over the rows no pair rule removed.

    Marks the rows of `REFERENCE_SKIPPED_RECORD_BUCKET`, the one 116.2 class that
    excludes and whose rows are all vep-rs's, in ``masked_rows``; counts the three
    breakend shapes without masking them. Returns the per-bucket statistics: ``count``
    rows recognised, ``excluded`` rows set aside, ``records`` distinct (file, Location)
    pairs the bucket reached.
    """
    stats = {
        rule.bucket: {"count": 0, "excluded": 0, "records": 0}
        for rule in REFERENCE_116_2_ONE_SIDED_RULES
    }
    records: dict[str, set[tuple[str, str]]] = defaultdict(set)

    def hit(bucket: str, idx: int, row: dict[str, str], exclude: bool) -> None:
        stats[bucket]["count"] += 1
        records[bucket].add((row["file_name"], row["location"]))
        if exclude:
            stats[bucket]["excluded"] += 1
            masked_rows.add(idx)

    # The one-sided rows of each side per (file, Location, allele), for the FASTA shape,
    # which reads what the other side wrote under the same allele.
    perl_keys: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    rust_keys: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    for idx, row in enumerate(rows):
        if idx in masked_rows:
            continue
        key = (row["file_name"], row["location"], row["allele"])
        if row["source"] == "missing_in_rust":
            perl_keys[key].append(idx)
        elif row["source"] == "extra_in_rust":
            rust_keys[key].append(idx)

    # The breakend records whose own-end rows the reference lost: it kept rows under a
    # bracket allele naming a mate on another chromosome and wrote nothing under any
    # own-end allele.
    own_rows_lost: set[tuple[str, str]] = set()
    for (file_name, loc), (perl_alleles, _rust_alleles, mates) in scope.breakends.items():
        own = _strip_chr(loc.split(":", 1)[0])
        if (
            any(bracket_mate(a) is not None for a in perl_alleles)
            and not any(is_own_breakend_allele(a) for a in perl_alleles)
            and any(m != own for m in mates)
        ):
            own_rows_lost.add((file_name, loc))

    for idx, row in enumerate(rows):
        if idx in masked_rows or row["source"] != "extra_in_rust":
            continue
        group = (row["file_name"], row["location"])
        allele = row["allele"]
        if group in scope.absent:
            # The reference wrote nothing at the Location; the merge admitted the group
            # only on a structural vep-rs row, and every row of the record goes with it.
            hit(REFERENCE_SKIPPED_RECORD_BUCKET, idx, row, exclude=True)
            continue
        breakend = scope.breakends.get(group)
        if breakend is None or row["feature"] == _INTERGENIC_FEATURE:
            continue
        if is_own_breakend_allele(allele):
            if group in own_rows_lost:
                hit(BND_OWN_ROWS_LOST_BUCKET, idx, row, exclude=False)
            continue
        perl_alleles, _rust_alleles, _mates = breakend
        mate = bracket_mate(allele)
        # A synonym is any mate name that is not a name the cache files a chromosome under:
        # the reference resolved it for the overlap and refused it at the slice test.
        if (
            mate is not None
            and perl_alleles
            and allele not in perl_alleles
            and mate[0] not in CACHE_CONTIG_NAMES
        ):
            hit(BND_SYNONYM_MATE_ROWS_LOST_BUCKET, idx, row, exclude=False)

    # The reference's only row at a key is its intergenic row while vep-rs wrote
    # transcript rows there and no intergenic row: the shape of a slice named after a
    # FASTA whose names differ from the cache's. Counted on the reference's row.
    for key, perl_idx in perl_keys.items():
        rust_idx = rust_keys.get(key)
        if not rust_idx or len(perl_idx) != 1:
            continue
        perl_row = rows[perl_idx[0]]
        if perl_row["feature"] == _INTERGENIC_FEATURE and all(
            rows[i]["feature"] != _INTERGENIC_FEATURE for i in rust_idx
        ):
            hit(FASTA_NAMED_SLICE_INTERGENIC_BUCKET, perl_idx[0], perl_row, exclude=False)

    for bucket in stats:
        stats[bucket]["records"] = len(records[bucket])
    return stats


def filter_snp_indel_intended_divergences(
    discordant_path: Path,
    open_out: Path | None = None,
    scope_features: ScopeFeatures | None = None,
    *,
    reference_release: str,
    record_scope: ReferenceRecordScope | None = None,
) -> tuple[int, int, dict[str, dict[str, int]]]:
    """Read discordant.tsv, classify intended-divergence swap pairs, and exclude
    the ones an Ensembl VEP self-contradiction accounts for.

    ``reference_release`` selects the one-sided pass (`REFERENCE_REGISTRIES`); the pair
    rules are the same under both. Under ``"115.2"`` the pass is the scope class below
    and ``record_scope`` is refused; under ``"116.2"`` it is
    `apply_116_2_one_sided_rules` over ``record_scope`` and ``scope_features`` is
    refused, so a caller cannot run one release's reference against the other's class.

    With ``scope_features`` (from ``compare_sorted_key_files``) the one-sided rows of a
    record Perl marks `vep_skip` take the structural-variant comparator's
    transcript-selection class, under its two guards: a vep-rs-only row on a transcript
    Perl did not name is set aside only when Perl named at least one transcript for the
    record and vep-rs named more; a Perl-only row on a transcript vep-rs did not name is
    set aside only when vep-rs named at least as many. A pair on a transcript both
    engines named stays charged, and a record below --max_sv_size with a supported type
    never enters the scope. The class is reported as ``sv_transcript_selection`` in the
    bucket statistics: ``count`` is every one-sided row on a record in scope, ``excluded``
    the rows set aside, ``perl`` and ``rust`` those rows per side.

    Returns ``(excluded_perl_count, excluded_rust_count, bucket_stats)`` where
    ``bucket_stats`` maps each classification bucket to
    ``{"count": pairs classified here, "excluded": how many of those the
    adjusted F1 removes}``. A swap pair is a
    ``(location, allele, feature, feature_type)`` present in BOTH
    ``missing_in_rust`` and ``extra_in_rust`` with different consequence sets.

    With ``open_out`` set, the rows of ``discordant_path`` that no excluding rule
    removed are written there in their original order and layout: both rows of
    every open pair plus every unpaired row. Its row count is therefore the
    adjusted one-sided total, ``(perl - excluded_perl) + (rust - excluded_rust)``.

    CLASSIFICATION AND EXCLUSION ARE SEPARATE, and that separation is the point.
    Every bucket below reports a divergence shape worth naming; only the shapes
    an Ensembl VEP mechanism produces are removed from the denominators. Naming
    a shape and crediting vep-rs for it are separate acts: a counter fused with
    its exclusion would mask every splice-family pair while the cited mechanism
    explains only a fraction of them.

    The buckets, in first-match classification order, are declared in
    ``EXCLUSION_REGISTRY`` above, which also carries each one's exclusion status
    and taxonomy class. Exclusion is driven by the registry's ``excludes`` rules
    evaluated INDEPENDENTLY of that chain, because
    ``_is_covered_splice_region_swap`` is a property of the pair that cuts across
    four buckets rather than lining up with any one of them.

    Every excluding rule's predicate is a subset of some classifying rule's
    predicate, so an excluded pair always lands in a bucket; the assertion after
    the loop is the backstop for that.

    Two filters are deliberately absent from the registry.

    A 3'UTR stop_lost trade-off ("vep-rs emits stop_lost where Perl emits
    3'UTR-only") has no Perl source-line citation behind it: Perl reads the UTR
    sequence from the cache's Bio::EnsEMBL::Slice, and vep-rs does the same
    (``vefc.{three,five}_prime_utr`` in the JSON cache, read by
    ``compute_{three,five}_prime_utr_sequence`` in
    ``crates/vep-effects/src/coding.rs`` with FASTA as the fallback), so
    cache-only mode reaches the Perl-parity verdict and the shape matches no
    tuple.

    An extended start co-emission predicate, one admitting any pair where Perl
    carries a start term vep-rs lacks, absorbs besides the strict shape: pairs
    where vep-rs drops the OTHER start term; pairs where vep-rs misses
    start_lost outright (P=5UTR,start_lost,start_retained / R=5UTR;
    P=frameshift,start_lost / R=frameshift); unrelated
    protein-altering-vs-inframe shapes (P=5UTR,coding_sequence_variant /
    R=5UTR,start_lost); and the inverted shape (P=frameshift,start_retained /
    R=frameshift,start_lost) where Perl emits ONLY start_retained, a separate
    predicate-overlap shape. Leaving those unmasked keeps the vep-rs
    divergences visible and countable. The strict
    ``_is_start_cooccurrence_swap`` is the only start filter, cited to
    Utils/VariationEffect.pm:851-878, :947-963 and :1028-1075 in
    ensembl-variation release/115 (23c76f60).
    """
    if reference_release not in REFERENCE_REGISTRIES:
        raise ValueError(
            f"unknown reference release {reference_release!r}; one of {REFERENCE_RELEASES}"
        )
    if reference_release == "115.2" and record_scope is not None:
        raise ValueError("record_scope is the 116.2 collector; the 115.2 pass reads scope_features")
    if reference_release == "116.2" and scope_features is not None:
        raise ValueError("scope_features is the 115.2 collector; the 116.2 pass reads record_scope")

    # The key carries feature_type because the SCORED tuple is (location, allele,
    # feature, feature_type, consequence_set): keyed without it, two rows differing
    # only by feature_type collapse onto one key and the nested loop counts their
    # cross product. The bound assertion after the loop is the backstop.
    #
    # Each side holds row indices into `rows` rather than consequence sets, so the
    # open rows can be written back in file order without a second pass.
    perl_by_key: dict[tuple[str, str, str, str], list[int]] = defaultdict(list)
    rust_by_key: dict[tuple[str, str, str, str], list[int]] = defaultdict(list)
    rows: list[dict[str, str]] = []
    n_perl_rows = 0
    n_rust_rows = 0

    with discordant_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = list(reader.fieldnames or ())
        for row in reader:
            key = (
                row["location"],
                row["allele"],
                row["feature"],
                row.get("feature_type", ""),
            )
            idx = len(rows)
            rows.append(row)
            if row["source"] == "missing_in_rust":
                perl_by_key[key].append(idx)
                n_perl_rows += 1
            elif row["source"] == "extra_in_rust":
                rust_by_key[key].append(idx)
                n_rust_rows += 1

    matched_keys = perl_by_key.keys() & rust_by_key.keys()
    excluded_perl = 0
    excluded_rust = 0
    masked_rows: set[int] = set()
    bucket_stats: dict[str, dict[str, int]] = {
        rule.bucket: {"count": 0, "excluded": 0} for rule in CLASSIFICATION_CHAIN
    }

    for key in matched_keys:
        # Each key has exactly 1 csq per side (unique 5-tuples from sorted-merge),
        # so the cross-product is 1x1.
        for pi in perl_by_key[key]:
            pcsq = rows[pi]["consequence_set"]
            for ri in rust_by_key[key]:
                rcsq = rows[ri]["consequence_set"]
                bucket = next(
                    (
                        rule.bucket
                        for rule in CLASSIFICATION_CHAIN
                        if rule.matches(pcsq, rcsq, key[1])
                    ),
                    None,
                )
                excluded = any(rule.matches(pcsq, rcsq, key[1]) for rule in EXCLUDING_RULES)
                if bucket is not None:
                    bucket_stats[bucket]["count"] += 1
                    if excluded:
                        bucket_stats[bucket]["excluded"] += 1
                elif excluded:
                    # Every excluding predicate is a subset of some classifying
                    # predicate, so this is unreachable. It is asserted rather than
                    # ignored because the alternative is a silent divergence between
                    # `excluded_perl` and the per-bucket `excluded` sum, which would
                    # make the reported mask coverage unattributable to any bucket.
                    raise AssertionError(
                        f"pair excluded but classified into no bucket: "
                        f"perl={pcsq!r} rust={rcsq!r}"
                    )
                if excluded:
                    excluded_perl += 1
                    excluded_rust += 1
                    masked_rows.update((pi, ri))

    # Neither side's exclusion count can exceed the discordant rows available on
    # that side. Adjusted F1 subtracts these from the raw denominators, so an
    # over-count inflates F1 silently; nothing else in the pipeline would notice.
    if excluded_perl > n_perl_rows or excluded_rust > n_rust_rows:
        raise AssertionError(
            f"intended-divergence exclusions exceed available discordant rows "
            f"(perl {excluded_perl} > {n_perl_rows} or rust {excluded_rust} > "
            f"{n_rust_rows}); a key is matching more than one consequence set per "
            f"side and the cross product is over-counting"
        )

    # Every excluded pair is attributed to exactly one bucket. A bucket may
    # classify without excluding, so the check compares the per-bucket EXCLUDED
    # sum, not the classified count, against the pairs the mask removed.
    attributed = sum(v["excluded"] for v in bucket_stats.values())
    if attributed != excluded_perl:
        raise AssertionError(
            f"per-bucket excluded counts sum to {attributed} but "
            f"{excluded_perl} pairs were excluded; the mask is not "
            f"attributable to the buckets that report it"
        )

    if scope_features:
        selection = {"count": 0, "excluded": 0, "perl": 0, "rust": 0}
        for idx, row in enumerate(rows):
            if idx in masked_rows:
                continue
            feats = scope_features.get((row["location"], row["allele"]))
            if feats is None:
                continue
            perl_feats, rust_feats = feats
            selection["count"] += 1
            feature = row["feature"]
            if row["source"] == "extra_in_rust":
                # vep-rs named a transcript Perl did not: the class only where Perl
                # annotated the record at all and named fewer transcripts, so a
                # vep-rs false positive on a record Perl annotated in full stays charged.
                if feature in perl_feats or not perl_feats or len(rust_feats) <= len(perl_feats):
                    continue
                excluded_rust += 1
                selection["rust"] += 1
            elif row["source"] == "missing_in_rust":
                # Perl named a transcript vep-rs did not: Perl's batch loaded a region
                # vep-rs's set does not reach only when vep-rs named at least as many,
                # so a vep-rs under-annotation stays charged.
                if feature in rust_feats or not rust_feats or len(rust_feats) < len(perl_feats):
                    continue
                excluded_perl += 1
                selection["perl"] += 1
            else:
                continue
            selection["excluded"] += 1
            masked_rows.add(idx)
        bucket_stats[SV_TRANSCRIPT_SELECTION_BUCKET] = selection

    if record_scope is not None:
        one_sided_stats = apply_116_2_one_sided_rules(rows, masked_rows, record_scope)
        excluded_rust += one_sided_stats[REFERENCE_SKIPPED_RECORD_BUCKET]["excluded"]
        bucket_stats.update(one_sided_stats)

    if open_out is not None:
        with open_out.open("w", encoding="utf-8", newline="") as out:
            out.write("\t".join(columns) + "\n")
            for idx, row in enumerate(rows):
                if idx not in masked_rows:
                    out.write("\t".join(row[c] for c in columns) + "\n")

    return excluded_perl, excluded_rust, bucket_stats


@dataclass
class FileMetrics:
    file_name: str
    perl_tuple_count: int
    rust_tuple_count: int
    intersection_count: int
    precision: float
    recall: float
    f1: float
    missing_in_rust: int
    extra_in_rust: int
    top_missing_consequence_sets: list[tuple[str, int]]
    top_extra_consequence_sets: list[tuple[str, int]]
    context_collapsed_precision: float
    context_collapsed_recall: float
    context_collapsed_f1: float
    context_collapsed_missing_in_rust: int
    context_collapsed_extra_in_rust: int


# What each mode needs. Validated in `main()` rather than with argparse's
# `required=True`, because the three modes need different subsets and a flag marked
# required in one mode has to be optional in the others.
MODE_REQUIRED: dict[str, tuple[str, ...]] = {
    "compare": ("perl_dir", "rust_dir", "report_dir"),
    "by-consequence-class": ("perl_dir", "discordant_in", "dataset", "class_csv", "sweep_id", "binary_md5"),
    "pool-consequence-class": ("class_csv", "pooled_csv"),
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare Perl and Rust VEP outputs")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--by-consequence-class",
        dest="mode",
        action="store_const",
        const="by-consequence-class",
        help="Per-consequence-class mode: count one (dataset, engine)'s per-term VEP "
        "totals from --perl-dir and its per-term one-sided counts from --discordant-in, "
        "then write per-term P/R/F1 and all four counts to --class-csv. Derives from "
        "a run's saved reports; runs no engine.",
    )
    mode.add_argument(
        "--pool-consequence-class",
        dest="mode",
        action="store_const",
        const="pool-consequence-class",
        help="Pool a --class-csv built across datasets into the per-term --pooled-csv, "
        "and report the floor over terms clearing --min-tuples.",
    )
    parser.set_defaults(mode="compare")
    parser.add_argument("--perl-dir", help="Directory with Perl VEP outputs")
    parser.add_argument(
        "--rust-dir", help="Directory with Rust VEP outputs (compare mode only)"
    )
    parser.add_argument("--report-dir", help="Directory for output reports")
    parser.add_argument(
        "--discordant-in",
        help="Per-class mode: path to an EXISTING discordant.tsv to read the per-term "
        "one-sided counts from. Distinct from --discordant-tsv, which names the file "
        "compare mode WRITES inside --report-dir.",
    )
    parser.add_argument(
        "--dataset", help="Per-class mode: dataset label written into each row"
    )
    parser.add_argument(
        "--assembly-label",
        default="",
        help="Per-class mode: assembly label written into each row. Separate from "
        "--assembly, which selects the canonical-contig filter.",
    )
    parser.add_argument(
        "--engine",
        choices=["vep-rs", "fastvep"],
        default="vep-rs",
        help="Per-class mode: which engine the discordance report's extra_in_rust side "
        "belongs to. Constrained to a fixed set because a typo would emit a CSV with no "
        "rows under the expected engine name, which reads as a missing measurement.",
    )
    parser.add_argument(
        "--class-csv",
        help="Per-class mode: per-(dataset, term) CSV to write. Pool mode: the same "
        "file, read as input.",
    )
    parser.add_argument(
        "--sweep-id",
        help="Per-class mode: the sweep whose discordance report is being decomposed, "
        "written into every row so the CSV carries the pin its aggregate rows carry.",
    )
    parser.add_argument(
        "--binary-md5",
        help="Per-class mode: digest (or release tag) of the engine build the discordance "
        "report was produced by, written into every row.",
    )
    parser.add_argument(
        "--vep-totals-cache",
        help="Per-class mode: a CSV caching the per-term VEP totals of --perl-dir. It is "
        "read when it records the same source files (name, size, modification time) "
        "and contig filter, otherwise recomputed and rewritten; those totals are the "
        "one input that does not change between runs against one ground truth.",
    )
    parser.add_argument(
        "--pooled-csv",
        help="Pool mode: pooled per-(engine, term) CSV to write.",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="Per-class mode: append to --class-csv instead of overwriting it, which is "
        "how a multi-dataset loop builds one file.",
    )
    parser.add_argument(
        "--min-tuples",
        type=int,
        default=PER_CLASS_MIN_TUPLES,
        help=f"Pool mode: minimum pooled VEP tuples for a term to be eligible for the "
        f"reported floor (default {PER_CLASS_MIN_TUPLES}, pre-committed). Terms below "
        f"it are still written to the CSV with their counts.",
    )
    parser.add_argument(
        "--per-term",
        help="Pool mode: a per-term TSV of the one-sided counts the mask excludes for "
        "--per-term-engine's run (columns term, vep_only_excluded, engine_only_excluded). "
        "They fill the pooled CSV's adj_f1 column for that engine; every other engine's "
        "adj_f1 is left empty.",
    )
    parser.add_argument(
        "--per-term-engine",
        default="vep-rs",
        help="Pool mode: the engine --per-term describes (default vep-rs).",
    )
    parser.add_argument(
        "--glob",
        default="*.txt",
        help="Filename glob for output pairing (default: *.txt)",
    )
    parser.add_argument(
        "--summary-json",
        default="summary.json",
        help="Summary JSON filename (inside --report-dir)",
    )
    parser.add_argument(
        "--summary-md",
        default="summary.md",
        help="Summary markdown filename (inside --report-dir)",
    )
    parser.add_argument(
        "--discordant-tsv",
        default="discordant.tsv",
        help="Discordant tuple TSV filename (inside --report-dir). The rows no "
        "excluding rule removes are also written beside it as <stem>_open<suffix> "
        "(discordant_open.tsv by default).",
    )
    parser.add_argument(
        "--fail-below-f1",
        type=float,
        default=None,
        help="Exit non-zero when aggregate F1 is below this threshold",
    )
    parser.add_argument(
        "--assembly",
        choices=["GRCh37", "GRCh38"],
        default=None,
        help="Reference assembly. When set with --canonical-contigs (default), tuples "
        "with non-canonical contig in Location are dropped before F1 calc as a "
        "defensive complement to prepare_benchmark_vcfs.py's input filter.",
    )
    parser.add_argument(
        "--canonical-contigs",
        action="store_true",
        default=True,
        help="Defensive contig-allowlist filter (DEFAULT ON). Drops 5-tuples whose "
        "Location contig is not in the canonical set for --assembly. Catches alt-contig "
        "tuples that slip past the input filter (engine cache annotations, etc.).",
    )
    parser.add_argument(
        "--no-canonical-contigs",
        dest="canonical_contigs",
        action="store_false",
        help="Opt out of defensive canonical-contigs filter (rare).",
    )
    parser.add_argument(
        "--reference-release",
        choices=REFERENCE_RELEASES,
        default=DEFAULT_REFERENCE_RELEASE,
        help="The Ensembl VEP release that produced --perl-dir, which selects the "
        "one-sided divergence classes of the adjusted F1 (REFERENCE_REGISTRIES). Under "
        "115.2 the classes are the published ones; under 116.2 (default) a record the "
        "reference dropped as vep_skip is set aside whole and the three breakend shapes "
        "of 116.1's chromosome-aware overlap are counted. The pair rules are the same "
        "under both. run_concordance.sh and run_clone_measurement.sh derive the value "
        "from the reference set's provenance.",
    )
    return parser.parse_args(argv)


def require_mode_args(args: argparse.Namespace) -> None:
    """Fail on a missing argument the selected mode needs, naming the mode.

    argparse cannot express this: the same flag is mandatory in one mode and
    meaningless in another. Raising SystemExit(2) matches argparse's own exit code
    for a usage error.
    """
    missing = [
        "--" + name.replace("_", "-")
        for name in MODE_REQUIRED[args.mode]
        if not getattr(args, name, None)
    ]
    if missing:
        raise SystemExit(
            f"error: mode {args.mode} requires {', '.join(missing)}"
        )


# Canonical contig set: universal
# Ensembl-style for both GRCh37 and GRCh38. Must match prepare_benchmark_vcfs.py
# CANONICAL_CONTIGS frozenset and scripts/concordance/run_concordance.sh
# CANONICAL_CONTIGS_LIST. M and MT are equivalent mitochondrial names.
#
# Rationale: input-prep auto-strips the UCSC chr- prefix so all engine
# outputs emit Ensembl-style Locations. This defensive output-side filter is
# therefore an Ensembl-only membership test on engine output Locations.
CANONICAL_CONTIGS = frozenset([str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"])
# The names the Ensembl caches file chromosomes under, for the breakend synonym test: a
# mate written `M` (or `chrM`) is a synonym of `MT` there, as any other name is.
CACHE_CONTIG_NAMES = CANONICAL_CONTIGS - {"M"}


def canonical_set_for_assembly(assembly: str | None) -> frozenset[str] | None:
    """Return the canonical contig set if the assembly is recognized.

    The same Ensembl-style set is used for both GRCh37 and
    GRCh38 (the assembly difference is in the genome-sequence content, not
    contig naming). Both GRCh37 and GRCh38 return the universal set; any
    other value, including None, returns None and disables the filter.
    """
    if assembly in ("GRCh37", "GRCh38"):
        return CANONICAL_CONTIGS
    return None


def normalize_consequence_set(raw: str) -> str:
    terms = [term.strip() for term in raw.split(",") if term.strip()]
    return ",".join(sorted(set(terms)))


def collapse_transcript_context_terms(raw: str) -> str:
    terms = [term.strip() for term in raw.split(",") if term.strip()]
    kept = sorted({term for term in terms if term not in TRANSCRIPT_CONTEXT_TERMS})
    if kept:
        return ",".join(kept)
    return "__context_only__"


def iter_vep_key_lines(
    path: Path,
    consequence_transform: Callable[[str], str] | None = None,
    canonical_contigs: frozenset[str] | None = None,
) -> Iterator[str]:
    """Yield normalized tuple keys as tab-separated strings.

    Key format (5 columns):
      location, allele, feature, feature_type, consequence_set

    Streaming, one row at a time, so a multi-million-row output never has to
    fit in memory.
    When `canonical_contigs` is provided, tuples whose Location contig is not in
    the set are silently dropped (defensive filter complementing the input-prep
    canonical-contigs filter in prepare_benchmark_vcfs.py).
    """
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 7:
                continue

            location = fields[1].strip()
            allele = fields[2].strip()
            feature = fields[4].strip()
            feature_type = fields[5].strip()
            consequence_set = normalize_consequence_set(fields[6].strip())
            if consequence_transform is not None:
                consequence_set = consequence_transform(consequence_set)

            if not (location and allele and feature):
                continue

            if canonical_contigs is not None:
                # Location format: "chr:pos" or "chr:start-end"; chrom is everything
                # before the first ":".
                contig = location.split(":", 1)[0]
                if contig not in canonical_contigs:
                    continue

            yield "\t".join((location, allele, feature, feature_type, consequence_set))


def _sort_unique_stream_to_file(keys: Iterator[str], out_path: Path) -> None:
    """Write sorted unique keys to `out_path` using external `sort -u` if available."""
    sort_exe = shutil.which("sort")
    if sort_exe is None:
        # Fallback: in-memory for environments without `sort`.
        uniq = sorted(set(keys))
        out_path.write_text("\n".join(uniq) + ("\n" if uniq else ""), encoding="utf-8")
        return

    env = os.environ.copy()
    # Byte-order collation is REQUIRED, not preferred, and must be forced rather
    # than defaulted. The sort-merge downstream compares keys by codepoint, so the
    # two sides have to agree on ordering or the merge walks past matching keys and
    # under-counts the intersection. `setdefault` does not override an LC_ALL the
    # caller exported, which is the whole failure mode: under en_US.UTF-8, `sort`
    # orders "1:117329" before "10:110755" while Python's `sorted` (and C collation)
    # order them the other way. A leaked locale can only LOSE intersections, so it
    # depresses F1 silently rather than erroring.
    env["LC_ALL"] = "C"
    env.pop("LC_COLLATE", None)  # LC_ALL wins over this, but drop it to be explicit
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as out_handle:
        proc = subprocess.Popen(
            [sort_exe, "-u"],
            stdin=subprocess.PIPE,
            stdout=out_handle,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
        assert proc.stdin is not None
        try:
            for key in keys:
                proc.stdin.write(key)
                proc.stdin.write("\n")
        finally:
            proc.stdin.close()
        # Use wait() + manual stderr read instead of communicate() to avoid
        # Python 3.9 ValueError when communicate() tries to flush closed stdin.
        proc.wait()
        err = proc.stderr.read() if proc.stderr else ""
        if proc.returncode != 0:
            raise RuntimeError(f"sort -u failed (rc={proc.returncode}): {err.strip()}")


def prf(
    intersection: int, rust_count: int, perl_count: int
) -> tuple[float, float, float]:
    precision = intersection / rust_count if rust_count else 0.0
    recall = intersection / perl_count if perl_count else 0.0
    if precision + recall == 0.0:
        return precision, recall, 0.0
    return precision, recall, 2 * precision * recall / (precision + recall)


def top_consequence_buckets(
    keys: Iterable[TupleKey], limit: int = 10
) -> list[tuple[str, int]]:
    counter: Counter[str] = Counter()
    for key in keys:
        counter[key[4]] += 1
    return counter.most_common(limit)


# Per-consequence-class F1
#
# Answers "is the aggregate F1 hiding one badly-called consequence class?". The
# quantity that answers it is the FLOOR: the lowest F1 any single Sequence Ontology
# term attains.
#
# It is derived from ARCHIVED artefacts with no engine re-run, which is what makes
# the measurement affordable at all. Three inputs per (dataset, engine):
#
#   * per-term VEP totals, counted over the canonical Perl output's sorted-unique
#     tuple keys -- the same keys the aggregate F1 is scored on;
#   * per-term VEP-only counts, from the run's `discordant.tsv` rows whose `source`
#     is `missing_in_rust`;
#   * per-term engine-only counts, from the same file's `extra_in_rust` rows.
#
# and two derivations:
#
#   intersection_t = VEP_t - VEPonly_t
#   engine_t       = intersection_t + engineonly_t
#
# The first derivation is exact because a tuple is either in the intersection or in
# the VEP-only set, with nothing else to be. The second rests on a property of the
# tuple key: `consequence_set` is one of its five columns, so an intersection tuple
# carries the SAME term set on both sides by construction, and counting term t on
# the VEP side of the intersection counts it on the engine side too. A term-level
# disagreement is never an intersection tuple -- it is a swap PAIR, one row on each
# one-sided list -- so it lands in `VEPonly_t` and `engineonly_t` and never in the
# intersection. That is why the engine's own output is not needed here.
#
# A tuple whose term set holds several terms contributes to each of them, so the
# per-term columns do NOT partition the aggregate and must not be summed back into
# it.

# A per-class F1 over a handful of tuples carries no usable precision (a term
# appearing 3 times with one discordant tuple reads F1 0.67, which is noise, not a
# mis-called class), while at n >= 1,000 a single discordant tuple moves F1 by at
# most 0.001. Terms below it are REPORTED with their counts and excluded from the
# floor only. The threshold is fixed independently of the data: raising it after
# seeing the values would choose the threshold that makes the floor look best.
PER_CLASS_MIN_TUPLES = 1000


def _round_half_up(x: float, places: int) -> str:
    """Render a float at `places` decimals under ROUND_HALF_UP.

    Fixed at ROUND_HALF_UP rather than Python's default banker's rounding so that a
    consumer re-deriving a published per-class F1 from this CSV's own integers gets
    the string the CSV carries, including on an exact-half value.
    """
    return str(Decimal(str(x)).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))


@dataclass(frozen=True)
class PerClassCounts:
    """One Sequence Ontology term's concordance, for one engine.

    Only the three counts are stored, because only they are independent;
    `intersection` and `engine_tuples` derive from them. That matters for pooling:
    summing three independent counts across datasets and re-deriving gives exactly
    the same answer as summing the derived ones, since both derivations are linear,
    and storing one truth instead of two removes any way for them to disagree.
    """

    term: str
    vep_tuples: int
    vep_only: int
    engine_only: int

    def __post_init__(self) -> None:
        # A term cannot have more VEP-only tuples than VEP tuples. When it does, the
        # `discordant.tsv` does not belong to the Perl output it was paired with --
        # a different run, a different canonical-contigs setting, or a different
        # assembly. Raising here is the point: the alternative is a NEGATIVE
        # intersection quietly producing an F1 above 1, or a positive one produced
        # from two unrelated runs, and neither is visible in the emitted CSV.
        if self.vep_only > self.vep_tuples:
            raise AssertionError(
                f"term {self.term!r}: {self.vep_only} VEP-only tuples against "
                f"{self.vep_tuples} VEP tuples. The discordance report does not "
                f"correspond to this Perl output -- check the run id, the assembly "
                f"and the canonical-contigs setting on both."
            )
        for name in ("vep_tuples", "vep_only", "engine_only"):
            if getattr(self, name) < 0:
                raise AssertionError(f"term {self.term!r}: negative {name}")

    @property
    def intersection(self) -> int:
        return self.vep_tuples - self.vep_only

    @property
    def engine_tuples(self) -> int:
        return self.intersection + self.engine_only

    @property
    def precision(self) -> float:
        return self.intersection / self.engine_tuples if self.engine_tuples else 0.0

    @property
    def recall(self) -> float:
        """Intersection over the VEP-side total, or 0.0 when that total is EMPTY.

        `vep_tuples == 0` means VEP never emits this term, so `V_t` is empty and recall is
        undefined rather than zero. 0.0 is reported deliberately: the rows it produces are
        scored-engine-only VCF 4.4 copy-number terms, the per-class CSV carries them with
        `vep_tuples` 0 so a reader sees that VEP has none, and a figure ordering terms by
        the reference engine's discordance does not plot them. Raising here instead would
        break those consumers.
        """
        return self.intersection / self.vep_tuples if self.vep_tuples else 0.0

    def f1(self) -> float:
        """F1 as ``2I / (V + E)``, NOT as ``2PR / (P + R)``.

        The two are algebraically identical and can differ in the last float bits.
        A consumer re-deriving this value does so as
        `2 * intersection / (vep + engine)`, so computing it the same way here
        removes the question of whether a boundary value rounds to a different sixth
        decimal in the published tables and in the CSV.
        """
        denom = self.vep_tuples + self.engine_tuples
        return 2 * self.intersection / denom if denom else 0.0

    def f1_str(self, places: int = 6) -> str:
        return _round_half_up(self.f1(), places)


def count_terms_in_sorted_keys(path: Path) -> Counter[str]:
    """Count, per Sequence Ontology term, the tuples in a sorted-unique key file.

    The file is the output of `_sort_unique_stream_to_file` over
    `iter_vep_key_lines`, so each line is one SCORED tuple and duplicates are
    already gone. Counting the raw engine output instead would double-count any
    line the engine emits twice, and the aggregate F1 does not, so the two
    measurements would disagree on the same data.
    """
    counter: Counter[str] = Counter()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 5:
                raise ValueError(
                    f"{path}: expected 5 tab-separated key fields, got {len(fields)}: "
                    f"{line[:120]!r}"
                )
            for term in fields[4].split(","):
                if term:
                    counter[term] += 1
    return counter


_VEP_TOTALS_CACHE_PREFIX = "# vep_totals_fingerprint="


def vep_totals_fingerprint(
    perl_files: list[Path], canonical_contigs: frozenset[str] | None
) -> str:
    """Identity of a per-term VEP total: the source files' names, sizes and
    modification times plus the contig filter, the only inputs the count depends on."""
    files = ";".join(
        f"{p.name}:{p.stat().st_size}:{p.stat().st_mtime_ns}" for p in perl_files
    )
    contigs = ",".join(sorted(canonical_contigs)) if canonical_contigs is not None else "none"
    return f"{files}|canonical={contigs}"


def read_vep_totals_cache(path: Path, fingerprint: str) -> Counter[str] | None:
    """The cached per-term VEP totals when ``path`` exists and records exactly
    ``fingerprint``; None on a missing, stale or unreadable cache, which the caller
    treats as a miss and recomputes."""
    if not path.is_file():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            first = handle.readline().rstrip("\n")
            if first != _VEP_TOTALS_CACHE_PREFIX + fingerprint:
                return None
            totals: Counter[str] = Counter()
            for row in csv.DictReader(handle):
                totals[row["term"]] = int(row["vep_tuples"])
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return totals


def write_vep_totals_cache(path: Path, fingerprint: str, totals: Counter[str]) -> None:
    """Write the cache for ``read_vep_totals_cache``, atomically, so a concurrent
    reader sees the old file or the new one and never a partial one."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        handle.write(_VEP_TOTALS_CACHE_PREFIX + fingerprint + "\n")
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["term", "vep_tuples"])
        for term in sorted(totals):
            writer.writerow([term, totals[term]])
    os.replace(tmp, path)


def count_terms_in_discordance(path: Path) -> tuple[Counter[str], Counter[str]]:
    """Count per-term one-sided tuples in a `discordant.tsv`.

    Returns ``(vep_only, engine_only)``, keyed on the `source` column:
    `missing_in_rust` rows are tuples VEP emitted and the engine did not,
    `extra_in_rust` rows the reverse. Any other `source` value is a schema change
    and is rejected rather than silently dropped, because a dropped side reads as
    perfect concordance on that side.
    """
    vep_only: Counter[str] = Counter()
    engine_only: Counter[str] = Counter()
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"source", "consequence_set"} - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: discordance report lacks column(s) {sorted(missing)}")
        for row in reader:
            source = row["source"]
            if source == "missing_in_rust":
                target = vep_only
            elif source == "extra_in_rust":
                target = engine_only
            else:
                raise ValueError(
                    f"{path}: unrecognised source {source!r}; expected "
                    f"'missing_in_rust' or 'extra_in_rust'"
                )
            for term in normalize_consequence_set(row["consequence_set"]).split(","):
                if term:
                    target[term] += 1
    return vep_only, engine_only


def derive_per_class_counts(
    vep_totals: Counter[str],
    vep_only: Counter[str],
    engine_only: Counter[str],
) -> list[PerClassCounts]:
    """Assemble one `PerClassCounts` per term seen on any of the three inputs.

    A term present only in `engine_only` is KEPT, with `vep_tuples` 0. That is an
    engine emitting a term VEP never emits anywhere in the dataset, which is a real
    finding and must not vanish; its recall is undefined and it falls below any
    tuple threshold, so it cannot become the published floor by accident.
    """
    terms = set(vep_totals) | set(vep_only) | set(engine_only)
    return sorted(
        (
            PerClassCounts(
                term=term,
                vep_tuples=vep_totals.get(term, 0),
                vep_only=vep_only.get(term, 0),
                engine_only=engine_only.get(term, 0),
            )
            for term in terms
        ),
        key=lambda c: (-c.vep_tuples, c.term),
    )


def pool_per_class_counts(rows: Iterable[PerClassCounts]) -> list[PerClassCounts]:
    """Sum per-term counts across datasets, then re-derive.

    Pooling is a micro-average over the union of the measured tuples, the same
    construction as a tuple-weighted micro-F1. It is the published grain for the
    floor, pre-committed before the values existed: a minimum over up to 41x6
    per-(term, dataset) cells is a worst-case outlier from a tiny denominator rather
    than a property of the engine.
    """
    totals: dict[str, list[int]] = {}
    for row in rows:
        acc = totals.setdefault(row.term, [0, 0, 0])
        acc[0] += row.vep_tuples
        acc[1] += row.vep_only
        acc[2] += row.engine_only
    return sorted(
        (
            PerClassCounts(term=term, vep_tuples=v, vep_only=vo, engine_only=eo)
            for term, (v, vo, eo) in totals.items()
        ),
        key=lambda c: (-c.vep_tuples, c.term),
    )


@dataclass(frozen=True)
class PerClassFloor:
    """The published floor and the two facts that qualify it."""

    f1: str | None
    term: str | None
    n_eligible: int
    n_below_threshold: int


def per_class_floor(
    rows: Iterable[PerClassCounts], min_tuples: int = PER_CLASS_MIN_TUPLES
) -> PerClassFloor:
    """The lowest F1 among terms clearing `min_tuples`, and the term attaining it.

    Ties break on the term name so the named term is reproducible rather than
    dependent on dict ordering; several terms sitting at the same six-decimal value
    is the expected case, not an edge case.
    """
    rows = list(rows)
    eligible = [r for r in rows if r.vep_tuples >= min_tuples]
    below = len(rows) - len(eligible)
    if not eligible:
        return PerClassFloor(None, None, 0, below)
    worst = min(eligible, key=lambda r: (Decimal(r.f1_str()), r.term))
    return PerClassFloor(worst.f1_str(), worst.term, len(eligible), below)


_PER_CLASS_BY_DATASET_COLUMNS = (
    "dataset",
    "assembly",
    "engine",
    "sweep_id",
    "binary_md5",
    "term",
    "vep_tuples",
    "vep_only",
    "engine_only",
    "intersection",
    "engine_tuples",
    "precision",
    "recall",
    "f1",
)

# The pooled CSV's schema. Downstream consumers read `engine`, `term`, `vep_tuples`,
# `engine_tuples`, `intersection` and `f1` from it by name, so renaming any of those
# six breaks a consumer rather than only this file. The last three columns are the
# adjusted view: `vep_only_excluded` and `engine_only_excluded` are the one-sided tuples
# the adjusted-F1 mask excludes per term (from `per_term.tsv`) and `adj_f1` is the F1
# with those removed from both sides, so a reader re-derives it from the row's own
# integers the way `f1` re-derives from the first three counts. All three are filled for
# the engine `--per-term` describes and empty for every other engine.
_PER_CLASS_POOLED_COLUMNS = (
    "engine",
    "sweep_id",
    "binary_md5",
    "term",
    "vep_tuples",
    "vep_only",
    "engine_only",
    "intersection",
    "engine_tuples",
    "precision",
    "recall",
    "f1",
    "vep_only_excluded",
    "engine_only_excluded",
    "adj_f1",
)


def read_per_term_exclusions(path: Path) -> dict[str, tuple[int, int]]:
    """Per term, the (VEP-only, engine-only) one-sided counts the mask excludes.

    Read from a per-term TSV for one run (columns term, vep_only_excluded,
    engine_only_excluded). A term absent from the file has no excluded pair, so its
    adjusted F1 equals its raw F1.
    """
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        needed = {"term", "vep_only_excluded", "engine_only_excluded"}
        missing = needed - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} lacks column(s) {sorted(missing)}")
        return {
            row["term"]: (int(row["vep_only_excluded"]), int(row["engine_only_excluded"]))
            for row in reader
        }


def adjusted_f1_str(counts: PerClassCounts, excluded: tuple[int, int], places: int = 6) -> str:
    """F1 with `excluded` one-sided pairs removed from both sides, rendered like `f1`.

    The intersection is untouched and each side's one-sided count shrinks by the pairs
    the mask excludes, so
    adj = 2I / (2I + (vep_only - ex_v) + (engine_only - ex_e)). An exclusion larger than
    the one-sided count it applies to means the per-term file belongs to another run.
    """
    ex_v, ex_e = excluded
    if ex_v > counts.vep_only or ex_e > counts.engine_only:
        raise AssertionError(
            f"term {counts.term!r}: per_term.tsv excludes {ex_v}/{ex_e} pairs against "
            f"{counts.vep_only}/{counts.engine_only} one-sided tuples; the per-term file "
            "does not belong to this run"
        )
    inter = counts.intersection
    denom = 2 * inter + (counts.vep_only - ex_v) + (counts.engine_only - ex_e)
    return _round_half_up(2 * inter / denom if denom else 0.0, places)


def _per_class_row(counts: PerClassCounts, engine: str, **extra: str) -> dict[str, str]:
    return {
        **extra,
        "engine": engine,
        "term": counts.term,
        "vep_tuples": str(counts.vep_tuples),
        "vep_only": str(counts.vep_only),
        "engine_only": str(counts.engine_only),
        "intersection": str(counts.intersection),
        "engine_tuples": str(counts.engine_tuples),
        "precision": _round_half_up(counts.precision, 6),
        "recall": _round_half_up(counts.recall, 6),
        "f1": counts.f1_str(),
    }


def write_per_class_by_dataset_csv(
    path: Path,
    counts: Iterable[PerClassCounts],
    *,
    dataset: str,
    assembly: str,
    engine: str,
    append: bool,
    sweep_id: str,
    binary_md5: str,
) -> int:
    """Write (or append) one (dataset, engine)'s per-term rows.

    Append mode is how a multi-dataset loop builds one file. The header is written only
    when the file is being created, and an append into a file whose header does not
    match is rejected: a silently mismatched header would put values in the wrong
    columns, which reads as a plausible CSV.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.is_file() and path.stat().st_size > 0
    if append and exists:
        with path.open("r", encoding="utf-8", newline="") as handle:
            header = next(csv.reader(handle), [])
        if tuple(header) != _PER_CLASS_BY_DATASET_COLUMNS:
            raise ValueError(
                f"{path}: cannot append, existing header {header} does not match "
                f"{list(_PER_CLASS_BY_DATASET_COLUMNS)}"
            )
    mode = "a" if (append and exists) else "w"
    written = 0
    with path.open(mode, encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(_PER_CLASS_BY_DATASET_COLUMNS))
        if mode == "w":
            writer.writeheader()
        for row in counts:
            writer.writerow(
                _per_class_row(
                    row, engine, dataset=dataset, assembly=assembly,
                    sweep_id=sweep_id, binary_md5=binary_md5,
                )
            )
            written += 1
    return written


def read_per_class_by_dataset_csv(
    path: Path,
) -> tuple[dict[str, list[PerClassCounts]], dict[str, tuple[str, str]]]:
    """Read the per-(dataset, term) CSV back, grouped by engine, with each engine's pin.

    Every row of one engine must carry the same (sweep_id, binary_md5): a per-class file
    assembled from two generations of one engine would pool two binaries into one
    published floor, so it is refused.
    """
    by_engine: dict[str, list[PerClassCounts]] = defaultdict(list)
    pins: dict[str, set[tuple[str, str]]] = defaultdict(set)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        needed = {"engine", "term", "vep_tuples", "vep_only", "engine_only", "sweep_id", "binary_md5"}
        missing = needed - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: lacks column(s) {sorted(missing)}")
        for row in reader:
            by_engine[row["engine"]].append(
                PerClassCounts(
                    term=row["term"],
                    vep_tuples=int(row["vep_tuples"]),
                    vep_only=int(row["vep_only"]),
                    engine_only=int(row["engine_only"]),
                )
            )
            pins[row["engine"]].add((row["sweep_id"], row["binary_md5"]))
    mixed = {e: sorted(v) for e, v in pins.items() if len(v) != 1}
    if mixed:
        raise ValueError(
            f"{path}: rows of one engine carry more than one (sweep_id, binary_md5) pin, "
            f"so they are not one measurement: {mixed}"
        )
    return dict(by_engine), {e: next(iter(v)) for e, v in pins.items()}


def write_per_class_pooled_csv(
    path: Path,
    pooled: dict[str, list[PerClassCounts]],
    pins: dict[str, tuple[str, str]],
    adjusted: dict[str, dict[str, tuple[int, int]]] | None = None,
) -> int:
    """Write the pooled CSV, one row per (engine, term).

    Engines are emitted in a fixed order with vep-rs first, and terms by descending
    VEP tuple count. Each row carries its engine's (sweep_id, binary_md5) pin.
    `adjusted` maps an engine to its per-term mask exclusions; that engine's rows carry
    `adj_f1` and the others' are empty.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    order = sorted(pooled, key=lambda e: (e != "vep-rs", e))
    adjusted = adjusted or {}
    written = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(_PER_CLASS_POOLED_COLUMNS))
        writer.writeheader()
        for engine in order:
            sweep_id, binary_md5 = pins[engine]
            exclusions = adjusted.get(engine)
            for row in pooled[engine]:
                extra = {
                    "sweep_id": sweep_id, "binary_md5": binary_md5,
                    "vep_only_excluded": "", "engine_only_excluded": "", "adj_f1": "",
                }
                if exclusions is not None:
                    ex_v, ex_e = exclusions.get(row.term, (0, 0))
                    extra["vep_only_excluded"] = str(ex_v)
                    extra["engine_only_excluded"] = str(ex_e)
                    extra["adj_f1"] = adjusted_f1_str(row, (ex_v, ex_e))
                writer.writerow(_per_class_row(row, engine, **extra))
                written += 1
    return written


def write_discordant_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    columns = [
        "file_name",
        "source",
        "location",
        "allele",
        "feature",
        "feature_type",
        "consequence_set",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary_markdown(path: Path, aggregate: dict, files: list[dict]) -> None:
    lines = []
    lines.append("# Concordance Summary")
    lines.append("")
    lines.append("## Aggregate")
    lines.append("")
    lines.append(f"- Perl tuples: {aggregate['perl_tuple_count']}")
    lines.append(f"- Rust tuples: {aggregate['rust_tuple_count']}")
    lines.append(f"- Intersection: {aggregate['intersection_count']}")
    lines.append(f"- Precision: {aggregate['precision']:.6f}")
    lines.append(f"- Recall: {aggregate['recall']:.6f}")
    lines.append(f"- F1: {aggregate['f1']:.6f}")
    adjusted = aggregate.get("adjusted", {})
    if adjusted:
        if "reference_release" in adjusted:
            lines.append(f"- Reference release: Ensembl VEP {adjusted['reference_release']}")
        # The bucket line reports CLASSIFICATION and so sits outside the exclusion
        # gate below: a bucket counts the pairs matching its shape whether or not
        # the mask removes them, and a suite can have buckets with zero exclusions.
        # Building `cats_str` inside the gate and emitting it outside would raise
        # NameError on such a suite rather than print an empty line.
        excl_cats = adjusted.get("excluded_categories", {})
        cats_str = ", ".join(
            f"{k}: {v['count']} ({v['excluded']} excluded)"
            for k, v in sorted(excl_cats.items(), key=lambda x: -x[1]["count"])
            if v["count"]
        )
        if cats_str:
            lines.append(f"- Divergence buckets: {cats_str}")
        if adjusted.get("excluded_perl_tuples", 0) > 0:
            adj_f1 = adjusted.get("f1", 0.0)
            excl_p = adjusted.get("excluded_perl_tuples", 0)
            excl_r = adjusted.get("excluded_rust_tuples", 0)
            lines.append(
                f"- **Adjusted F1: {adj_f1:.6f}** (excluded {excl_p}P + {excl_r}R intended-divergence tuples)"
            )
    collapsed = aggregate.get("context_collapsed", {})
    if collapsed:
        lines.append(
            f"- Context-collapsed F1: {collapsed.get('f1', 0.0):.6f} "
            f"(collapsed terms: {', '.join(collapsed.get('collapsed_terms', []))})"
        )
    classes = adjusted.get("divergence_classes") if adjusted else None
    if classes:
        lines.append("")
        lines.append("## Divergence classes")
        lines.append("")
        lines.append("| Class | Treatment | Definition | Ensembl VEP source |")
        lines.append("|---|---|---|---|")
        for name, c in classes.items():
            verdict = "sets aside" if c["excludes"] else "counts only"
            lines.append(f"| `{name}` | {verdict} | {c['definition']} | {c['perl_citation']} |")
    lines.append("")
    lines.append("## Per File")
    lines.append("")
    lines.append(
        "| File | Perl Tuples | Rust Tuples | Intersection | Precision | Recall | F1 |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for item in files:
        lines.append(
            f"| {item['file_name']} | {item['perl_tuple_count']} | {item['rust_tuple_count']} "
            f"| {item['intersection_count']} | {item['precision']:.6f} "
            f"| {item['recall']:.6f} | {item['f1']:.6f} |"
        )
    lines.append("")
    lines.append("## Per File (Context-Collapsed)")
    lines.append("")
    lines.append("| File | Precision | Recall | F1 | Missing | Extra |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for item in files:
        lines.append(
            f"| {item['file_name']} | {item['context_collapsed_precision']:.6f} "
            f"| {item['context_collapsed_recall']:.6f} "
            f"| {item['context_collapsed_f1']:.6f} "
            f"| {item['context_collapsed_missing_in_rust']} "
            f"| {item['context_collapsed_extra_in_rust']} |"
        )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def _next_line(handle: TextIO) -> str | None:
    line = handle.readline()
    if not line:
        return None
    return line.rstrip("\n")


def compare_sorted_key_files(
    perl_sorted: Path,
    rust_sorted: Path,
    *,
    file_name: str,
    discordant_out: TextIO | None,
    track_consequence_buckets: bool,
    scope_features: ScopeFeatures | None = None,
    record_scope: ReferenceRecordScope | None = None,
) -> tuple[int, int, int, int, int, Counter[str], Counter[str]]:
    """Compare two sorted-unique key files (streaming merge).

    With ``scope_features`` given, every key line on a record in
    ``_batch_dependent_annotation_scope`` adds its feature to that record's set for its
    side, matched and one-sided lines alike, so the exclusion pass can compare the two
    engines' transcript sets for the record.

    With ``record_scope`` given (the 116.2 collector), every key line is observed with the
    side or sides that wrote it, and the collector's open Location group is closed when
    both streams are exhausted.

    Returns:
      perl_count, rust_count, intersection, missing, extra, missing_counter, extra_counter
    """
    perl_count = 0
    rust_count = 0
    intersection = 0
    missing = 0
    extra = 0
    missing_counter: Counter[str] = Counter()
    extra_counter: Counter[str] = Counter()

    extra_tmp: TextIO | None = None
    if discordant_out is not None:
        extra_tmp = tempfile.TemporaryFile(mode="w+", encoding="utf-8", newline="")

    with (
        perl_sorted.open("r", encoding="utf-8", errors="replace") as perl_h,
        rust_sorted.open("r", encoding="utf-8", errors="replace") as rust_h,
    ):
        p = _next_line(perl_h)
        r = _next_line(rust_h)
        while p is not None or r is not None:
            if p is not None and (r is None or p < r):
                perl_count += 1
                missing += 1
                if scope_features is not None:
                    _record_scope_feature(scope_features, p, 0)
                if record_scope is not None:
                    record_scope.observe(file_name, p, perl=True, rust=False)
                if track_consequence_buckets:
                    missing_counter[p.rsplit("\t", 1)[-1]] += 1
                if discordant_out is not None:
                    discordant_out.write(f"{file_name}\tmissing_in_rust\t{p}\n")
                p = _next_line(perl_h)
                continue

            if r is not None and (p is None or r < p):
                rust_count += 1
                extra += 1
                if scope_features is not None:
                    _record_scope_feature(scope_features, r, 1)
                if record_scope is not None:
                    record_scope.observe(file_name, r, perl=False, rust=True)
                if track_consequence_buckets:
                    extra_counter[r.rsplit("\t", 1)[-1]] += 1
                if extra_tmp is not None:
                    extra_tmp.write(f"{file_name}\textra_in_rust\t{r}\n")
                r = _next_line(rust_h)
                continue

            # Equal: intersection
            assert p is not None and r is not None
            perl_count += 1
            rust_count += 1
            intersection += 1
            if scope_features is not None:
                _record_scope_feature(scope_features, p, 0)
                _record_scope_feature(scope_features, p, 1)
            if record_scope is not None:
                record_scope.observe(file_name, p, perl=True, rust=True)
            p = _next_line(perl_h)
            r = _next_line(rust_h)

    if record_scope is not None:
        record_scope.close()

    if extra_tmp is not None and discordant_out is not None:
        extra_tmp.seek(0)
        shutil.copyfileobj(extra_tmp, discordant_out)
        extra_tmp.close()

    return (
        perl_count,
        rust_count,
        intersection,
        missing,
        extra,
        missing_counter,
        extra_counter,
    )


def resolve_canonical_contigs(args: argparse.Namespace) -> frozenset[str] | None:
    """Resolve the defensive contig filter, warning when it is off by omission.

    Factored so the compare and per-class paths cannot drift on it. They must agree:
    per-class VEP totals counted WITHOUT the filter against a discordance report
    produced WITH it would inflate every term's VEP total by the dropped alt-contig
    tuples, which reads as depressed recall on every class at once.
    """
    canonical_contigs = (
        canonical_set_for_assembly(args.assembly) if args.canonical_contigs else None
    )
    if args.canonical_contigs and canonical_contigs is None:
        # Soft-warn rather than hard-fail: the upstream input filter already
        # enforced the canonical-contigs constraint; this script is a defensive
        # complement, and a missing --assembly means the defensive layer is off
        # while the upstream filter applies. Print the warning so the
        # operator sees the asymmetry in the report.
        print(
            "WARNING: --canonical-contigs (default ON) requested but --assembly was "
            "not provided. Defensive output-side filter is OFF for this run. The "
            "upstream input-prep canonical-contigs filter (prepare_benchmark_vcfs.py) "
            "still applies. Pass --assembly GRCh37 or --assembly GRCh38 to enable "
            "the defensive layer, or --no-canonical-contigs to silence this warning.",
        )
    return canonical_contigs


def run_by_consequence_class(args: argparse.Namespace) -> int:
    """Count one (dataset, engine)'s per-term concordance from a run's saved reports."""
    perl_dir = Path(args.perl_dir).resolve()
    discordant_in = Path(args.discordant_in).resolve()
    class_csv = Path(args.class_csv).resolve()
    if not discordant_in.is_file():
        raise FileNotFoundError(f"discordance report not found: {discordant_in}")

    canonical_contigs = resolve_canonical_contigs(args)

    perl_files = sorted(perl_dir.glob(args.glob))
    if not perl_files:
        raise FileNotFoundError(f"No files found in {perl_dir} matching {args.glob}")

    # The VEP-side totals depend on the reference output alone, not on the run, so
    # a cache keyed on the reference files serves every run against the same
    # ground truth.
    cache = Path(args.vep_totals_cache).resolve() if args.vep_totals_cache else None
    fingerprint = vep_totals_fingerprint(perl_files, canonical_contigs)
    vep_totals = read_vep_totals_cache(cache, fingerprint) if cache else None
    if vep_totals is not None:
        print(f"VEP per-term totals: cache hit {cache}")
    else:
        # Sort-unique first, then count. Counting the raw stream would double-count any
        # duplicate line, and the aggregate F1 this measurement decomposes is scored on
        # the deduplicated key set.
        vep_totals = Counter()
        with tempfile.TemporaryDirectory(prefix="per_class_tmp_") as tmp:
            for perl_file in perl_files:
                sorted_keys = Path(tmp) / f"{perl_file.name}.perl.keys.sorted.txt"
                _sort_unique_stream_to_file(
                    iter_vep_key_lines(perl_file, canonical_contigs=canonical_contigs),
                    sorted_keys,
                )
                vep_totals.update(count_terms_in_sorted_keys(sorted_keys))
                sorted_keys.unlink(missing_ok=True)
        if cache:
            write_vep_totals_cache(cache, fingerprint, vep_totals)
            print(f"VEP per-term totals: cache miss, wrote {cache}")

    vep_only, engine_only = count_terms_in_discordance(discordant_in)
    counts = derive_per_class_counts(vep_totals, vep_only, engine_only)

    written = write_per_class_by_dataset_csv(
        class_csv,
        counts,
        dataset=args.dataset,
        assembly=args.assembly_label or (args.assembly or ""),
        engine=args.engine,
        append=args.append,
        sweep_id=args.sweep_id,
        binary_md5=args.binary_md5,
    )

    floor = per_class_floor(counts, args.min_tuples)
    print(f"wrote per-class counts: {class_csv} ({written} rows {args.engine}/{args.dataset})")
    print(
        f"terms: {len(counts)}  VEP tuples: {sum(vep_totals.values())}  "
        f"VEP-only: {sum(vep_only.values())}  engine-only: {sum(engine_only.values())}"
    )
    if floor.f1 is None:
        print(
            f"no term reaches {args.min_tuples} VEP tuples on this dataset, so it "
            f"contributes no per-dataset floor"
        )
    else:
        print(
            f"this dataset's floor: {floor.f1} on {floor.term} "
            f"({floor.n_eligible} terms at or above {args.min_tuples} tuples, "
            f"{floor.n_below_threshold} below). NOT the reported floor, which is "
            f"pooled across datasets -- see --pool-consequence-class."
        )
    return 0


def run_pool_consequence_class(args: argparse.Namespace) -> int:
    """Pool a per-(dataset, term) CSV into the per-term artefact."""
    class_csv = Path(args.class_csv).resolve()
    pooled_csv = Path(args.pooled_csv).resolve()
    if not class_csv.is_file():
        raise FileNotFoundError(f"per-class CSV not found: {class_csv}")

    by_engine, pins = read_per_class_by_dataset_csv(class_csv)
    if not by_engine:
        raise ValueError(f"{class_csv} carries no rows to pool")
    pooled = {engine: pool_per_class_counts(rows) for engine, rows in by_engine.items()}
    adjusted: dict[str, dict[str, tuple[int, int]]] | None = None
    if args.per_term:
        if args.per_term_engine not in pooled:
            raise ValueError(
                f"--per-term-engine {args.per_term_engine!r} has no rows in {class_csv}"
            )
        adjusted = {args.per_term_engine: read_per_term_exclusions(Path(args.per_term).resolve())}
    written = write_per_class_pooled_csv(pooled_csv, pooled, pins, adjusted)
    print(
        f"wrote pooled per-class F1: {pooled_csv} ({written} rows"
        + (f"; adj_f1 filled for {args.per_term_engine})" if adjusted else ")")
    )

    for engine in sorted(pooled, key=lambda e: (e != "vep-rs", e)):
        rows = pooled[engine]
        floor = per_class_floor(rows, args.min_tuples)
        if floor.f1 is None:
            print(
                f"{engine}: no term clears {args.min_tuples} pooled VEP tuples, so no "
                f"floor can be reported under the pre-committed rule"
            )
            continue
        print(
            f"{engine}: floor {floor.f1} on {floor.term}; "
            f"{floor.n_eligible} of {len(rows)} terms at or above {args.min_tuples} "
            f"pooled VEP tuples, {floor.n_below_threshold} below"
        )
    return 0


def main() -> int:
    args = parse_args()
    require_mode_args(args)
    if args.mode == "by-consequence-class":
        return run_by_consequence_class(args)
    if args.mode == "pool-consequence-class":
        return run_pool_consequence_class(args)

    perl_dir = Path(args.perl_dir).resolve()
    rust_dir = Path(args.rust_dir).resolve()
    report_dir = Path(args.report_dir).resolve()
    report_dir.mkdir(parents=True, exist_ok=True)

    canonical_contigs = resolve_canonical_contigs(args)

    perl_files = sorted(perl_dir.glob(args.glob))
    if not perl_files:
        raise FileNotFoundError(f"No files found in {perl_dir} matching {args.glob}")

    per_file_results: list[dict] = []

    perl_total = 0
    rust_total = 0
    intersection_total = 0
    collapsed_perl_total = 0
    collapsed_rust_total = 0
    collapsed_intersection_total = 0
    missing_counter_aggregate: Counter[str] = Counter()
    extra_counter_aggregate: Counter[str] = Counter()
    missing_rust_pairs: list[str] = []
    extra_rust_pairs: list[str] = []

    discordant_path = report_dir / args.discordant_tsv
    discordant_path.parent.mkdir(parents=True, exist_ok=True)
    discordant_out = discordant_path.open("w", encoding="utf-8", newline="")
    # The selected release's one-sided collector. The 115.2 one is keyed without the file:
    # a location names its chromosome, so files cannot collide.
    scope_features: ScopeFeatures | None = {} if args.reference_release == "115.2" else None
    record_scope: ReferenceRecordScope | None = (
        ReferenceRecordScope() if args.reference_release == "116.2" else None
    )
    discordant_out.write(
        "\t".join(
            [
                "file_name",
                "source",
                "location",
                "allele",
                "feature",
                "feature_type",
                "consequence_set",
            ]
        )
        + "\n"
    )

    for perl_file in perl_files:
        rust_file = rust_dir / perl_file.name
        if not rust_file.exists():
            missing_rust_pairs.append(perl_file.name)
            continue

        with tempfile.TemporaryDirectory(
            prefix="compare_vep_tmp_", dir=str(report_dir)
        ) as tmp:
            tmp_dir = Path(tmp)

            # Raw comparison keys
            perl_sorted = tmp_dir / f"{perl_file.name}.perl.keys.sorted.txt"
            rust_sorted = tmp_dir / f"{perl_file.name}.rust.keys.sorted.txt"
            _sort_unique_stream_to_file(
                iter_vep_key_lines(perl_file, canonical_contigs=canonical_contigs),
                perl_sorted,
            )
            _sort_unique_stream_to_file(
                iter_vep_key_lines(rust_file, canonical_contigs=canonical_contigs),
                rust_sorted,
            )

            (
                perl_count,
                rust_count,
                intersection_count,
                missing_count,
                extra_count,
                missing_counter,
                extra_counter,
            ) = compare_sorted_key_files(
                perl_sorted,
                rust_sorted,
                file_name=perl_file.name,
                discordant_out=discordant_out,
                track_consequence_buckets=True,
                scope_features=scope_features,
                record_scope=record_scope,
            )

            # Context-collapsed comparison keys (no discordant TSV, no bucket tracking)
            perl_collapsed_sorted = (
                tmp_dir / f"{perl_file.name}.perl.keys.collapsed.sorted.txt"
            )
            rust_collapsed_sorted = (
                tmp_dir / f"{perl_file.name}.rust.keys.collapsed.sorted.txt"
            )
            _sort_unique_stream_to_file(
                iter_vep_key_lines(
                    perl_file,
                    consequence_transform=collapse_transcript_context_terms,
                    canonical_contigs=canonical_contigs,
                ),
                perl_collapsed_sorted,
            )
            _sort_unique_stream_to_file(
                iter_vep_key_lines(
                    rust_file,
                    consequence_transform=collapse_transcript_context_terms,
                    canonical_contigs=canonical_contigs,
                ),
                rust_collapsed_sorted,
            )
            (
                c_perl_count,
                c_rust_count,
                c_intersection,
                c_missing,
                c_extra,
                _,
                _,
            ) = compare_sorted_key_files(
                perl_collapsed_sorted,
                rust_collapsed_sorted,
                file_name=perl_file.name,
                discordant_out=None,
                track_consequence_buckets=False,
            )

        precision, recall, f1 = prf(intersection_count, rust_count, perl_count)
        collapsed_precision, collapsed_recall, collapsed_f1 = prf(
            c_intersection,
            c_rust_count,
            c_perl_count,
        )

        perl_total += perl_count
        rust_total += rust_count
        intersection_total += intersection_count
        collapsed_perl_total += c_perl_count
        collapsed_rust_total += c_rust_count
        collapsed_intersection_total += c_intersection

        missing_counter_aggregate.update(missing_counter)
        extra_counter_aggregate.update(extra_counter)

        missing_top = missing_counter.most_common(10)
        extra_top = extra_counter.most_common(10)

        metrics = FileMetrics(
            file_name=perl_file.name,
            perl_tuple_count=perl_count,
            rust_tuple_count=rust_count,
            intersection_count=intersection_count,
            precision=precision,
            recall=recall,
            f1=f1,
            missing_in_rust=missing_count,
            extra_in_rust=extra_count,
            top_missing_consequence_sets=missing_top,
            top_extra_consequence_sets=extra_top,
            context_collapsed_precision=collapsed_precision,
            context_collapsed_recall=collapsed_recall,
            context_collapsed_f1=collapsed_f1,
            context_collapsed_missing_in_rust=c_missing,
            context_collapsed_extra_in_rust=c_extra,
        )
        per_file_results.append(metrics.__dict__)

    discordant_out.close()

    for rust_file in sorted(rust_dir.glob(args.glob)):
        if not (perl_dir / rust_file.name).exists():
            extra_rust_pairs.append(rust_file.name)

    aggregate_precision, aggregate_recall, aggregate_f1 = prf(
        intersection_total, rust_total, perl_total
    )
    collapsed_precision, collapsed_recall, collapsed_f1 = prf(
        collapsed_intersection_total, collapsed_rust_total, collapsed_perl_total
    )

    # --- Adjusted F1: exclude intended divergences ---
    open_path = discordant_path.with_name(
        f"{discordant_path.stem}_open{discordant_path.suffix}"
    )
    excl_perl, excl_rust, excl_cats = filter_snp_indel_intended_divergences(
        discordant_path,
        open_out=open_path,
        scope_features=scope_features,
        reference_release=args.reference_release,
        record_scope=record_scope,
    )
    adj_perl = perl_total - excl_perl
    adj_rust = rust_total - excl_rust
    # Intersection unchanged: excluded tuples are in P-only / R-only sets
    adj_precision, adj_recall, adj_f1 = prf(intersection_total, adj_rust, adj_perl)
    total_excluded = excl_perl + excl_rust

    aggregate = {
        "perl_tuple_count": perl_total,
        "rust_tuple_count": rust_total,
        "intersection_count": intersection_total,
        "precision": aggregate_precision,
        "recall": aggregate_recall,
        "f1": aggregate_f1,
        "top_missing_consequence_sets": missing_counter_aggregate.most_common(10),
        "top_extra_consequence_sets": extra_counter_aggregate.most_common(10),
        "context_collapsed": {
            "perl_tuple_count": collapsed_perl_total,
            "rust_tuple_count": collapsed_rust_total,
            "intersection_count": collapsed_intersection_total,
            "precision": collapsed_precision,
            "recall": collapsed_recall,
            "f1": collapsed_f1,
            "collapsed_terms": sorted(TRANSCRIPT_CONTEXT_TERMS),
        },
        "adjusted": {
            "perl_tuple_count": adj_perl,
            "rust_tuple_count": adj_rust,
            "intersection_count": intersection_total,
            "precision": adj_precision,
            "recall": adj_recall,
            "f1": adj_f1,
            "excluded_perl_tuples": excl_perl,
            "excluded_rust_tuples": excl_rust,
            "excluded_categories": excl_cats,
        },
    }
    # The 115.2 report keeps the published shape; a later release's names its reference
    # and defines every class it applied.
    if args.reference_release != "115.2":
        aggregate["adjusted"]["reference_release"] = args.reference_release
        aggregate["adjusted"]["divergence_classes"] = divergence_class_definitions(
            args.reference_release
        )

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "perl_dir": str(perl_dir),
        "rust_dir": str(rust_dir),
        "report_dir": str(report_dir),
        "glob": args.glob,
        "missing_rust_pairs": missing_rust_pairs,
        "extra_rust_pairs": extra_rust_pairs,
        "aggregate": aggregate,
        "files": per_file_results,
        "assembly": args.assembly,
        "canonical_contigs_filter_active": args.canonical_contigs
        and canonical_contigs is not None,
        "canonical_contigs_set": (
            sorted(canonical_contigs) if canonical_contigs is not None else None
        ),
    }

    summary_json_path = report_dir / args.summary_json
    summary_json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_summary_markdown(report_dir / args.summary_md, aggregate, per_file_results)
    if args.reference_release != "115.2":
        print(f"reference release: {args.reference_release}")

    print(f"wrote summary: {summary_json_path}")
    print(f"wrote markdown: {report_dir / args.summary_md}")
    print(f"wrote discordance: {discordant_path}")
    print(f"wrote open discordance: {open_path}")
    print(
        "aggregate: "
        f"precision={aggregate_precision:.6f} "
        f"recall={aggregate_recall:.6f} "
        f"f1={aggregate_f1:.6f}"
    )
    if total_excluded > 0:
        print(
            "adjusted: "
            f"precision={adj_precision:.6f} "
            f"recall={adj_recall:.6f} "
            f"f1={adj_f1:.6f} "
            f"(excluded {excl_perl}P + {excl_rust}R = {total_excluded} "
            f"intended-divergence tuples)"
        )
    print(
        "aggregate (context-collapsed): "
        f"precision={collapsed_precision:.6f} "
        f"recall={collapsed_recall:.6f} "
        f"f1={collapsed_f1:.6f}"
    )

    if args.fail_below_f1 is not None and aggregate_f1 < args.fail_below_f1:
        print(
            f"aggregate f1 {aggregate_f1:.6f} below threshold {args.fail_below_f1:.6f}",
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
