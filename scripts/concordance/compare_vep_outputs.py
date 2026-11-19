#!/usr/bin/env python3
"""Semantic concordance comparison for Perl VEP vs vep-rs outputs."""

from __future__ import annotations

import argparse
import copy
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path
from typing import Callable, Collection, Iterable, Iterator, Protocol, Sequence, TextIO

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


# `_is_start_cooccurrence_swap` is the only start filter. A looser predicate, one
# admitting any pair where Perl carries a start term vep-rs lacks, would also set
# aside pairs where vep-rs misses or drops a start term, which are vep-rs
# divergences.


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
# for a gVCF reference block, so the arms below have nothing to reach there. The Allele
# of a record whose type has no term is Perl's type string (Parser/VCF.pm:469-475):
# the ALT as written (`<CPX>`), or INFO/SVTYPE when that names the type (`CPX`, the
# gnomAD-SV spelling).
_MAX_SV_SIZE = 10_000_000
_VEP_SKIP_UNSUPPORTED_ALLELES = frozenset({"CPX", "CTX", "<CPX>", "<CTX>"})
SV_TRANSCRIPT_SELECTION_BUCKET = "sv_transcript_selection"


def _ten_digits(text: str) -> bool:
    """Whether ``text`` is written in the digits ``0`` to ``9`` alone, the digits a variant
    file writes a coordinate in; ``str.isdigit`` takes the digits of other scripts as well."""
    return text.isascii() and text.isdigit()


def _location_span(loc: str) -> tuple[str, int, int] | None:
    """``chrom:start-end`` or ``chrom:pos`` as (chrom, start, end), a coordinate the reference
    numified below zero included (`21:-4-295`); None when unparseable."""
    m = re.fullmatch(r"([^:]+):(-?\d+)(?:-(-?\d+))?", loc, re.ASCII)
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
    # `int` reads the digits of any script between the blanks it drops; a coordinate is
    # written in `0` to `9`.
    if not (start_s.strip().isascii() and end_s.strip().isascii()):
        return False
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


class RowEntry(Protocol):
    """One output row as a row-aware pair rule reads it: an entry of the format comparator
    (``compare_vep_formats.py``), or the `PairRow` this comparator builds from the row's first
    nine columns. ``key`` ends in the row's Feature_type after its last ``|``, ``fields``
    holds the row's ``cDNA_position`` and ``CDS_position``, and ``indel`` is whether the
    record's own alleles make the allele an insertion or a deletion, ``None`` where the
    first column spells no alleles. ``insertion`` is whether they make every sequence allele
    an insertion, the reference allele empty once the shared first base is dropped
    (``reference_allele_is_empty``), ``None`` where the first column spells no alleles.
    ``reference`` is the reference bases the allele stands in for at the entry's Location, as
    the record's own alleles give them by position (``replaced_reference``): from the first
    column or the format's own allele fields where they spell the alleles, else from the
    input records where the comparison is given them; empty for a pure insertion's entry,
    ``None`` where no record was read."""

    key: str
    location: str | None
    allele: str
    consequence: str
    fields: dict[str, str]
    indel: bool | None
    insertion: bool | None
    reference: str | None


def dash(s: str) -> str:
    """The default format's spelling of absence for an empty VCF or JSON value."""
    return "-" if s == "" else s


def allele_is_indel(reference: str, alternates: Collection[str], allele: str, chopped: bool = False) -> bool | None:
    """Whether ``allele`` is an insertion or a deletion of the record whose reference allele
    is ``reference`` and alternate alleles ``alternates``: the alternate the allele names and
    the reference differ in length, a difference no trimming of shared bases changes. A record
    with one alternate not spelt with a ``*`` has that alternate's shape on every row, since
    Ensembl VEP shortens a lone allele past the shared base (``minimise_alleles``: ``TCT`` to
    ``TCTC`` is written ``C``) and the allele then spells neither the REF nor the ALT. With
    several, Ensembl VEP drops the record's shared leading base only when every alternate not
    spelt with a ``*`` shares it (``Parser/VCF.pm`` lines 306-321 exempt ``<DEL:*>`` as they
    exempt ``*``), so ``A`` to ``AT,G`` keeps ``AT`` as the insertion's allele; a VCF line's
    REF and ALT are read the same way, and the parts of a JSON ``allele_string``, already so
    spelt (``chopped``), as they are, since dropping a base again from ``T/TT/-`` leaves no
    allele spelt ``TT``. ``None`` when the allele is not a sequence (a symbolic or a breakend
    allele), when no alternate names it, or when two alternates of different shapes name it."""
    if not re.fullmatch(r"[ACGTN-]+", allele):
        return None
    trimmed = [a for a in alternates if "*" not in a]
    if len(trimmed) == 1:
        return len(reference.strip("-")) != len(trimmed[0].strip("-"))
    # Parser/VCF.pm:300-323 drop the shared first base only of a record that holds an ALT of
    # another length than its REF, a ``*`` among them.
    holds_indel = any(len(a) != len(reference) for a in alternates)
    trim = not chopped and bool(reference) and holds_indel and all(a[:1] == reference[:1] for a in trimmed)
    ref = reference[1:] if trim else reference
    shapes: set[bool] = set()
    for alternate in alternates:
        alt = alternate[1:] if trim and "*" not in alternate else alternate
        if dash(alt) == allele:
            shapes.add(len(ref.strip("-")) != len(alt.strip("-")))
    return shapes.pop() if len(shapes) == 1 else None


def reference_allele_is_empty(reference: str, alternates: Collection[str], chopped: bool = False) -> bool:
    """Whether the reference allele of a record is empty once the record's shared first base is
    dropped, ``reference`` and ``alternates`` read as ``allele_is_indel`` reads them
    (``Parser/VCF.pm`` lines 300-323 drop the base where every alternate not spelt with a ``*``
    shares it and one ALT has another length than the REF). Every sequence allele of such a
    record is an insertion, which Ensembl VEP stores with a start one above its end. A
    reference allele that keeps a base past the shared first one is not empty here, a lone
    alternate Ensembl VEP shortens further (``trim_sequences``) included."""
    trimmed = [a for a in alternates if "*" not in a]
    holds_indel = any(len(a) != len(reference) for a in alternates)
    trim = not chopped and bool(reference) and holds_indel and all(a[:1] == reference[:1] for a in trimmed)
    return not (reference[1:] if trim else reference).strip("-")


# Ensembl VEP's Uploaded_variation for a VCF record without an ID (OutputFactory.pm:883-931):
# the chromosome, the position and the record's own alleles joined by ``/``.
UPLOADED_ALLELES_RE = re.compile(r"^\S+_[0-9]+_([ACGTNacgtn*-]+(?:/[^/\s]+)+)$")


def uploaded_alleles(uploaded: str) -> tuple[str, list[str]] | None:
    """The reference and alternate alleles an identifier-less ``Uploaded_variation`` spells
    (``21_10542414_A/T/AG``), or ``None`` for an identifier, which names no allele."""
    m = UPLOADED_ALLELES_RE.match(uploaded)
    if m is None:
        return None
    reference, *alternates = m.group(1).split("/")
    return reference, alternates


UPLOADED_POSITION_RE = re.compile(r"^\S+_([0-9]+)_[ACGTNacgtn*-]+(?:/[^/\s]+)+$")


def uploaded_position(uploaded: str) -> int | None:
    """The POS an identifier-less ``Uploaded_variation`` spells before its alleles
    (``21_10542414_A/T/AG``), or ``None`` for an identifier."""
    m = UPLOADED_POSITION_RE.match(uploaded)
    return integer_of(m.group(1), "Uploaded_variation position") if m else None


def replaced_reference(
    position: int, reference: str, alternates: Collection[str], location: str | None, allele: str
) -> str | None:
    """The reference bases an output row's allele stands in for, read by position from a
    record's own alleles. The record whose reference allele ``reference`` starts at
    ``position`` writes ``allele`` at ``location`` where one of its alternates is the
    reference allele with the Location's bases replaced by the allele (``-`` replaces them by
    nothing): those bases are returned. Where an alternate is the reference allele with the
    allele inserted between the Location's two neighbouring positions the row is a pure
    insertion's and the empty string is returned. ``None`` where the record writes no such
    row, where two alternates read it differently, and for an allele or a Location that is
    not a sequence's."""
    if location is None or not re.fullmatch(r"[ACGTN]+|-", allele):
        return None
    start, _, end = location.split(":", 1)[-1].partition("-")
    end = end or start
    if not (re.fullmatch(r"[0-9]+", start) and re.fullmatch(r"[0-9]+", end)):
        return None
    first = integer_of(start, "Location position") - position
    last = integer_of(end, "Location position") - position
    written = "" if allele == "-" else allele
    reference = reference.upper()
    found: set[str] = set()
    for alternate in alternates:
        spelt = alternate.upper()
        if 0 <= first <= last < len(reference) and spelt == reference[:first] + written + reference[last + 1 :]:
            found.add(reference[first : last + 1])
        if (
            written
            and last == first + 1
            and 0 <= first < len(reference)
            and spelt == reference[: first + 1] + written + reference[first + 1 :]
        ):
            found.add("")
    return found.pop() if len(found) == 1 else None


# The field value of an entry that lacks the field.
ABSENT = "<absent>"


def blank(value: str) -> bool:
    """Whether a field value is empty in any format's spelling: ``-``, the empty string or
    absent from the entry."""
    return value in ("", "-", ABSENT)


class ComparatorError(Exception):
    """A file this comparator cannot read or an option it cannot apply."""


def integer_of(digits: str, what: str) -> int:
    """The integer a string of digits writes. ``str.isdigit`` takes digits ``int`` does not
    read (a superscript) and ``int`` refuses a number longer than the interpreter converts:
    either is a malformed input and ends the run with the comparator's error."""
    try:
        return int(digits)
    except ValueError as exc:
        shown = digits if len(digits) <= 24 else digits[:24] + "..."
        raise ComparatorError(f"{what} {shown!r} is not an integer") from exc


def consequence_set(s: str) -> str:
    """Consequence terms as a sorted, de-duplicated, comma-joined set."""
    terms = sorted({t for t in re.split(r"[&,]", s) if t})
    return ",".join(terms)


def feature_kind(feature_type: str) -> str:
    """``JSON.pm:348-349``: the consequence list a feature type's entries go to
    (``Transcript`` to ``transcript``, ``RegulatoryFeature`` to ``regulatory_feature``, an
    undefined type to ``intergenic``); a list name is returned unchanged."""
    kind = (feature_type if feature_type not in ("", "-") else "intergenic").lower()
    if kind == "intergenic" or "_feature" in kind:
        return kind
    return kind.replace("feature", "_feature")


def indel_row(e: RowEntry) -> bool:
    """Whether an entry is an insertion's or a deletion's: an allele of ``-`` is one; else the
    record's own alleles decide where the parser read them (``Entry.indel``: the VCF and JSON
    entries, and a text row whose Uploaded_variation spells the alleles); else the Location
    decides, an indel's spanning more bases than the allele has (an insertion's two flanking
    bases against a one-base allele, a deletion's span against the kept bases; an insertion as
    long as its span, or an insertion at a one-base Location, reads as a substitution here, so
    a row under an identifier may miss the class)."""
    if e.allele == "-":
        return True
    if e.indel is not None:
        return e.indel
    if e.location is None or not re.fullmatch(r"[ACGTN]+", e.allele):
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    if not (_ten_digits(start) and _ten_digits(end)):
        return False
    return integer_of(end, "Location position") - integer_of(start, "Location position") + 1 != len(e.allele)


CODING_SEQUENCE_TERM = "coding_sequence_variant"


STOP_RETAINED_TERM = "stop_retained_variant"


SPLICE_REGION_TERM = "splice_region_variant"


def position_bounds(e: RowEntry, column: str, start_key: str, end_key: str) -> tuple[bool, bool]:
    """Whether the entry carries the start and the end of a transcript position: the two bounds
    of the text formats' ``column`` (``?`` for an end that maps to no transcript base, one number
    standing for both, ``-`` for neither, each read up to the ``/total`` suffix ``--total_length``
    appends, OutputFactory.pm:1726-1739), or the JSON ``start_key`` and ``end_key``, which an
    unmapped end leaves out."""
    position = e.fields.get(column)
    if position is None:
        return not blank(e.fields.get(start_key, ABSENT)), not blank(e.fields.get(end_key, ABSENT))
    position = position.partition("/")[0]
    if blank(position):
        return False, False
    start, sep, end = position.partition("-")
    return start != "?", (end if sep else start) != "?"


def stop_codon_test_coordinates(e: RowEntry) -> bool:
    """Whether the entry carries the coordinates the stop codon test edits at
    (Utils/VariationEffect.pm:1507-1508, ``return 0 unless $cdna_start && $cdna_end &&
    $cds_start``): both cDNA positions and the CDS start. A CDS end alone, the ``?-N`` of a
    span that begins in the 5' UTR, is not one of them."""
    return (
        all(position_bounds(e, "cDNA_position", "cdna_start", "cdna_end"))
        and position_bounds(e, "CDS_position", "cds_start", "cds_end")[0]
    )


def stop_codon_test_not_made(e: RowEntry) -> bool:
    """Whether the entry shows that ``_ins_del_stop_altered_cil`` (Utils/VariationEffect.pm
    :1488-1539) returns before it edits anything: the allele holds a letter other than ``A``,
    ``C``, ``G`` and ``T`` (:1498), the span lacks a cDNA position at one end, or it has no CDS
    start (:1508, ``stop_codon_test_coordinates``)."""
    if e.allele != "-" and not re.fullmatch(r"[ACGT]+", e.allele):
        return True
    return not stop_codon_test_coordinates(e)


# ``_overlaps_stop_codon_cil`` (Utils/VariationEffect.pm:1389-1431) takes the last three
# genomic positions of the coding region as the stop codon (:1418-1422). On a stop codon an
# intron splits, those positions are one or two intron bases beside the codon's last base, and
# the stop bases that end the exon before lie outside them: a deletion of one of those bases
# fails the overlap, ``stop_retained`` is false whatever the edited sequence holds, and
# ``coding_unknown`` (:1616-1649) writes ``coding_sequence_variant``. The candidate reads the
# overlap on the codon's own bases, the genomic bases of the last three coding positions, and
# makes the stop codon test: where the edited sequence holds a stop codon at the old stop
# position the row is ``stop_retained_variant``. The shape: a Transcript row of an insertion or
# deletion with both cDNA positions and a CDS start whose two Consequence sets are equal
# but for the reference's ``coding_sequence_variant`` against the candidate's
# ``stop_retained_variant``, both holding ``splice_region_variant``, which every exon base
# within three of an exon's end carries (``splice_region``, :616-627). A row away from an exon
# end, a row without its positions and the two terms the other way round stay differences.
def stop_codon_overlap_genomic_span_pair(e: RowEntry, a: RowEntry) -> bool:
    """Whether the reference entry ``e`` and the candidate entry ``a`` of one key are a row on
    which Ensembl VEP's genomic stop codon overlap misses a stop base that ends an exon: a
    Transcript entry of an insertion or deletion with both cDNA positions and a CDS start
    on both sides, ``splice_region_variant`` in both sets of terms, and the two sets equal but
    for the reference's ``coding_sequence_variant`` against the candidate's
    ``stop_retained_variant``."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript" or not indel_row(e):
        return False
    expected_terms = set(consequence_set(e.consequence).split(","))
    actual_terms = set(consequence_set(a.consequence).split(","))
    if expected_terms - actual_terms != {CODING_SEQUENCE_TERM} or actual_terms - expected_terms != {STOP_RETAINED_TERM}:
        return False
    if SPLICE_REGION_TERM not in expected_terms:
        return False
    return not any(stop_codon_test_not_made(entry) for entry in (e, a))


def length_changing_row(e: RowEntry) -> bool:
    """Whether an entry is an insertion's or a deletion's as ``indel_row`` reads one, an allele
    of other nucleotide letters than ``A``, ``C``, ``G``, ``T`` and ``N`` included: where the
    record's own alleles say nothing, such an allele against a Location of another length."""
    if indel_row(e):
        return True
    if e.indel is not None or e.location is None or not re.fullmatch(r"[A-Za-z]+", e.allele):
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    return _ten_digits(start) and _ten_digits(end) and int(end) - int(start) + 1 != len(e.allele)


FRAMESHIFT_TERM = "frameshift_variant"


# ``stop_retained`` (Utils/VariationEffect.pm:1290-1324) sends an insertion or deletion whose
# alternate peptide is undefined, empty or holds an ``X`` to its second arm (:1318): the
# genomic overlap with the stop codon (``_overlaps_stop_codon_cil``, :1389-1431) and not
# ``_ins_del_stop_altered_cil`` (:1488-1539). That routine sets its answer to 0 (:1495) and
# returns it before it edits anything for an allele holding a letter that names no base
# (:1498, ``seq_is_unambiguous_dna``) and for a span without a cDNA coordinate at one end or a
# CDS start (:1508), so :1318 reads a test the routine never makes as a stop codon kept and the
# row carries ``stop_retained_variant``: on a deletion that removes the stop codon together
# with a splice site or the transcript's end, among others. The cDNA twin behind ``stop_lost``
# (:1263) reads the same 0 as a stop codon that is not lost, so neither term describes such a
# record there. The candidate makes the test wherever the span has both cDNA coordinates and
# a CDS start, on the three bases the edited sequence holds at the old stop position, and
# reads a span the test cannot be made on as keeping nothing; the row then carries the terms
# its other predicates give: ``coding_sequence_variant`` from ``coding_unknown`` (:1616-1649),
# ``frameshift_variant`` for an inserted allele whose length is no multiple of three, or no
# coding term beside a ``start_lost``. The shape: a Transcript row of an insertion or deletion
# whose two Consequence sets are equal but for the reference's ``stop_retained_variant``,
# against nothing or one of ``coding_sequence_variant`` and ``frameshift_variant`` on the
# candidate's side, both rows showing that the test is not made (``stop_codon_test_not_made``).
# A row with both cDNA positions, a CDS start and an allele of bases (the test is made
# there and both engines read its answer), a candidate that keeps the term, any other term
# beside the two and a substitution's row stay differences. The row's IMPACT follows its
# first term, and under ``--hgvs`` or a plugin that reads ``frameshift_variant`` so do the
# fields those read; the predicate accepts the Consequence alone.
def stop_retained_without_stop_codon_test_pair(e: RowEntry, a: RowEntry) -> bool:
    """Whether the reference entry ``e`` and the candidate entry ``a`` of one key are a row on
    which Ensembl VEP writes ``stop_retained_variant`` for a stop codon test it does not make:
    a Transcript entry of an insertion or deletion, the reference's terms holding
    ``stop_retained_variant`` and the candidate's being the same without it, with nothing,
    ``coding_sequence_variant`` or ``frameshift_variant`` in its place, and both entries
    showing that the test is not made."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript" or not length_changing_row(e):
        return False
    expected_terms = set(consequence_set(e.consequence).split(","))
    actual_terms = set(consequence_set(a.consequence).split(","))
    if expected_terms - actual_terms != {STOP_RETAINED_TERM}:
        return False
    if actual_terms - expected_terms not in (set(), {CODING_SEQUENCE_TERM}, {FRAMESHIFT_TERM}):
        return False
    return stop_codon_test_not_made(e) and stop_codon_test_not_made(a)


STOP_GAINED_TERM = "stop_gained"


# The terms ``frameshift`` (Utils/VariationEffect.pm:1541-1567) reads away where ``stop_retained``
# does not hold: the in-frame pair through ``_get_codon_alleles`` (:1107, :1170) and
# ``protein_altering_variant`` (:390).
FRAMESHIFT_READ_TERMS = frozenset({"inframe_insertion", "inframe_deletion", "protein_altering_variant"})


# ``_ins_del_stop_altered_cil`` (Utils/VariationEffect.pm:1488-1539), where it makes its test,
# writes the allele into the coding sequence followed by the 3' UTR (:1521) and reads three
# bases, the ones at the old stop position (:1528-1535); a stop codon there is its answer "not
# altered", and :1318 writes ``stop_retained_variant``. Nothing asks whether the bases before
# that position are the reference's. A deletion that takes coding bases, the stop codon and
# 3' UTR bases moves later bases into the old stop position, and where those are a stop codon
# the row reads as a stop codon kept although a residue before it changed or an earlier stop
# codon ends the protein; so does an insertion before the stop codon on the reverse strand,
# whose start ``_overlaps_stop_codon_cil`` moves back by the inserted length (:1413). The
# candidate keeps the term only where every codon before the old stop position reads the
# reference's residue as well, and the row then carries the terms its other predicates give:
# ``coding_sequence_variant`` from ``coding_unknown`` (:1616-1649), or ``frameshift_variant``,
# with ``stop_gained`` where the alternate peptide gains a stop, in place of the in-frame or
# protein-altering term ``frameshift`` reads away. The shape: a Transcript row of an insertion
# or deletion, both rows showing that the test is made (both cDNA positions, a CDS start
# and an allele of bases), whose two Consequence sets are equal but for the reference's
# ``stop_retained_variant``, alone against ``coding_sequence_variant``, or alone or beside one
# of those terms against ``frameshift_variant`` with or without ``stop_gained``. A row whose
# test is not made, a candidate that keeps the term, any other term on either side and a
# substitution's row stay differences. The predicate accepts the Consequence alone.
def stop_retained_three_base_stop_codon_test_pair(e: RowEntry, a: RowEntry) -> bool:
    """Whether the reference entry ``e`` and the candidate entry ``a`` of one key are a row on
    which Ensembl VEP writes ``stop_retained_variant`` from the three bases at the old stop
    position alone: a Transcript entry of an insertion or deletion with both cDNA positions, a
    CDS start and an allele of bases on both sides, the reference's terms holding
    ``stop_retained_variant`` and the candidate's ``coding_sequence_variant`` in its place, or,
    where the reference's hold it alone or beside one of ``inframe_insertion``,
    ``inframe_deletion`` and ``protein_altering_variant``, ``frameshift_variant`` with or
    without ``stop_gained`` in the place of those terms."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript" or not length_changing_row(e):
        return False
    expected_terms = set(consequence_set(e.consequence).split(","))
    actual_terms = set(consequence_set(a.consequence).split(","))
    gone, come = expected_terms - actual_terms, actual_terms - expected_terms
    if come == {CODING_SEQUENCE_TERM}:
        if gone != {STOP_RETAINED_TERM}:
            return False
    elif come in ({FRAMESHIFT_TERM}, {STOP_GAINED_TERM, FRAMESHIFT_TERM}):
        beside = gone - {STOP_RETAINED_TERM}
        if STOP_RETAINED_TERM not in gone or len(beside) > 1 or not beside <= FRAMESHIFT_READ_TERMS:
            return False
    else:
        return False
    return not any(stop_codon_test_not_made(entry) for entry in (e, a))


STOP_LOST_TERM = "stop_lost"


# A record whose span holds all three bases of the stop codon and lacks a cDNA coordinate at
# one end or a CDS start takes the codon out with a splice site, an exon or the transcript's
# end. ``stop_lost`` (Utils/VariationEffect.pm:1234-1288) then asks ``_ins_del_stop_altered``,
# which returns 0 without those coordinates (:1376, :1453), and ``stop_retained`` asks its
# twin, which returns the same 0 (:1508) that :1318 reads as a stop codon kept: the row
# carries ``stop_retained_variant`` and no test is made. The candidate reads the sequence the
# record gives: the stop codon is lost unless that sequence shows it kept, so the row carries
# ``stop_lost`` where the given sequence ends before the old stop position or holds no stop
# codon there. The shape: a Transcript row of an insertion or deletion whose two Consequence
# sets are equal but for the reference's ``stop_retained_variant`` against the candidate's
# ``stop_lost``, both rows showing a span without its coordinates: a cDNA position missing
# whole or at one end, or no CDS start. A span that begins in the 5' UTR and ends on the stop
# codon's last base prints the CDS end alone (``?-N``), and :1508 leaves the test out for
# want of the start. A row with those three coordinates, any other term on either side, the
# two terms the other way round and a substitution's row stay differences. The predicate
# accepts the Consequence alone.
def stop_codon_removed_without_stop_codon_test_pair(e: RowEntry, a: RowEntry) -> bool:
    """Whether the reference entry ``e`` and the candidate entry ``a`` of one key are a row on
    which Ensembl VEP writes ``stop_retained_variant`` for a stop codon the record removes with
    no test made: a Transcript entry of an insertion or deletion, both entries lacking a cDNA
    position at one end or a CDS start, and the two sets of terms equal but for the
    reference's ``stop_retained_variant`` against the candidate's ``stop_lost``."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript" or not length_changing_row(e):
        return False
    expected_terms = set(consequence_set(e.consequence).split(","))
    actual_terms = set(consequence_set(a.consequence).split(","))
    if expected_terms - actual_terms != {STOP_RETAINED_TERM} or actual_terms - expected_terms != {STOP_LOST_TERM}:
        return False
    return not any(stop_codon_test_coordinates(entry) for entry in (e, a))


def insertion_row(e: RowEntry) -> bool:
    """Whether an entry's allele is a sequence (``-``, a symbolic or a breakend allele is
    refused) and the record's own alleles mark it an indel where the parser read them
    (``Entry.indel``, true of a deletion spelt as a sequence too, the ``A`` row of ``AT`` to
    ``A,G``, whose ``del`` notation ``clamped_insertion_position`` and ``shifted_insertion_end``
    both refuse), else the Location spans the two bases around an insertion's point."""
    if not re.fullmatch(r"[ACGTN]+", e.allele):
        return False
    if e.indel is not None:
        return e.indel
    if e.location is None:
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    return _ten_digits(start) and _ten_digits(end) and int(end) == int(start) + 1


def inserted_between_neighbours(e: RowEntry) -> bool:
    """Whether an entry ``insertion_row`` takes is an insertion between two neighbouring
    positions, the allele Ensembl VEP stores with a start one above its end and the only one
    whose start ``_overlaps_stop_codon_cil`` moves (Utils/VariationEffect.pm:1408). Its
    Location, where the format writes one, is two neighbouring positions; and, where the parser
    read the record's own alleles, its reference allele is empty once the shared first base is
    dropped (``RowEntry.insertion``), which refuses a deletion spelt as a sequence, a
    replacement by a sequence of another length, an insertion kept with its shared first base
    and an allele none of them names beside a reference allele that keeps a base. Under an
    identifier the allele has another length than two instead, which refuses the two-base substitution of the same Location;
    two bases replaced by a sequence of another length have an insertion's Location and
    allele there and are taken with the insertions."""
    if e.location is not None:
        start, _, end = e.location.split(":", 1)[-1].partition("-")
        if not (re.fullmatch(r"[0-9]+", start) and re.fullmatch(r"[0-9]+", end) and int(end) == int(start) + 1):
            return False
    if e.insertion is not None:
        return e.insertion
    return e.location is not None and len(e.allele) != 2


INFRAME_INSERTION_TERM = "inframe_insertion"


# ``_overlaps_stop_codon_cil`` (Utils/VariationEffect.pm:1389-1431) moves an insertion's start
# by the inserted length before it reads the overlap with the stop codon, forward on the
# forward strand and back on the reverse (:1408-1415). On the forward strand that puts the
# start past the insertion's end, and the overlap holds only for one base inserted between the
# codon's first and second base: an insertion before or inside the stop codon that leaves
# every base before the codon and a stop codon at the old stop position fails it,
# ``stop_retained`` (:1318) is false, and the row carries ``coding_sequence_variant`` from
# ``coding_unknown`` (:1616-1649) or, for an inserted length that is no multiple of three,
# ``frameshift_variant`` with or without ``stop_gained``. The candidate reads the insertion on
# the coding sequence, the same on either strand: it overlaps where it lies before or inside
# the codon and every codon before the old stop position keeps its residue. The stop
# codon test then finds the codon kept and the row is ``stop_retained_variant``; ``frameshift``
# returns 0 for a retained stop codon (:1549), so the ``inframe_insertion`` it reads away
# elsewhere stands beside the term. The shape: a Transcript row of an insertion between two
# neighbouring positions (``inserted_between_neighbours``) whose allele
# holds ``A``, ``C``, ``G``, ``T`` and ``N`` alone, the letters :1408 moves the start for, both
# rows carrying the coordinates the candidate's test edits at (both cDNA positions and a CDS
# start), whose two Consequence sets are equal but for the reference's
# ``coding_sequence_variant`` against the candidate's ``stop_retained_variant``, or the
# reference's ``frameshift_variant`` with or without ``stop_gained`` against the candidate's
# ``stop_retained_variant`` with or without ``inframe_insertion``. A deletion's row, a
# replacement's row where the record's own alleles show it, a two-base allele's row under an
# identifier, a row without those coordinates, the first pair beside ``splice_region_variant`` (the shape
# ``stop_codon_overlap_genomic_span_pair`` reads, on an unsplit codon beside an exon end as
# well), the terms the other way round (the reverse strand's wider range, which
# ``stop_retained_three_base_stop_codon_test_pair`` reads) and any other term on either side
# stay differences. The predicate accepts the Consequence alone.
def stop_codon_overlap_insertion_range_pair(e: RowEntry, a: RowEntry) -> bool:
    """Whether the reference entry ``e`` and the candidate entry ``a`` of one key are a row on
    which Ensembl VEP's genomic stop codon overlap misses an insertion that keeps the stop
    codon: a Transcript entry of an insertion between two neighbouring positions
    (``inserted_between_neighbours``) with both cDNA positions and a CDS start on
    both sides, the two sets of terms equal but for the reference's
    ``coding_sequence_variant`` against the candidate's ``stop_retained_variant`` away from
    ``splice_region_variant``, or the reference's ``frameshift_variant`` with or without
    ``stop_gained`` against the candidate's ``stop_retained_variant`` with or without
    ``inframe_insertion``."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript" or not insertion_row(e):
        return False
    if not inserted_between_neighbours(e):
        return False
    expected_terms = set(consequence_set(e.consequence).split(","))
    actual_terms = set(consequence_set(a.consequence).split(","))
    gone, come = expected_terms - actual_terms, actual_terms - expected_terms
    if gone == {CODING_SEQUENCE_TERM}:
        if come != {STOP_RETAINED_TERM} or SPLICE_REGION_TERM in expected_terms:
            return False
    elif gone in ({FRAMESHIFT_TERM}, {STOP_GAINED_TERM, FRAMESHIFT_TERM}):
        if come not in ({STOP_RETAINED_TERM}, {INFRAME_INSERTION_TERM, STOP_RETAINED_TERM}):
            return False
    else:
        return False
    return all(stop_codon_test_coordinates(entry) for entry in (e, a))


# The five splice terms whose windows lie inside an intron: the two bases at its 5' end, the
# two at its 3' end, its fifth base, its third to sixth base, and the bases 3 to 17 from its
# 3' end (Utils/VariationEffect.pm:558-614).
INTRON_WINDOW_SPLICE_TERMS: frozenset[str] = frozenset(
    {
        "splice_donor_variant",
        "splice_acceptor_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
    }
)
# The four of them beside which ``splice_region`` returns 0 (Utils/VariationEffect.pm:616-627).
SPLICE_REGION_WITHHOLDING_TERMS: frozenset[str] = INTRON_WINDOW_SPLICE_TERMS - {"splice_polypyrimidine_tract_variant"}
# ``within_intron`` (Utils/VariationEffect.pm:629-634) returns the flag ``intronic``; the
# intron itself, its end bases included, is that term's window.
INTRON_TERM = "intron_variant"
# The six terms ``_intron_effects`` sets from the differing regions and the candidate
# writes only where the record's span holds a base of the term's window.
PAST_THE_SPAN_TERMS: frozenset[str] = INTRON_WINDOW_SPLICE_TERMS | {INTRON_TERM}
SPLICE_DONOR_5TH_BASE_TERM = "splice_donor_5th_base_variant"
SPLICE_DONOR_REGION_TERM = "splice_donor_region_variant"
# A gap past which no window of an intron's end reaches a span or its regions.
_FAR = 10**9
# The row-aware rules of the stop codon classes, whose terms a row can differ by beside the
# six: one row can carry both differences.
STOP_CODON_ROW_PAIRS: tuple[Callable[[RowEntry, RowEntry], bool], ...] = (
    stop_codon_overlap_genomic_span_pair,
    stop_retained_without_stop_codon_test_pair,
    stop_retained_three_base_stop_codon_test_pair,
    stop_codon_removed_without_stop_codon_test_pair,
    stop_codon_overlap_insertion_range_pair,
)


def transcript_positions(e: RowEntry) -> tuple[str, str, str, str] | None:
    """The cDNA start and end and the CDS start and end an entry carries, as written (``?``
    for an end that maps to no transcript base, ``-`` for a position the entry lacks), or
    ``None`` where a position is not of that form. A text format's column is read up to the
    ``/total`` suffix ``--total_length`` appends, one number standing for both ends; the JSON
    keys are read as they are."""
    out: list[str] = []
    for column, start_key, end_key in (("cDNA_position", "cdna_start", "cdna_end"), ("CDS_position", "cds_start", "cds_end")):
        position = e.fields.get(column)
        if position is None:
            pair = [dash(str(e.fields.get(key, "-"))) for key in (start_key, end_key)]
        else:
            start, sep, end = position.partition("/")[0].partition("-")
            pair = ["-", "-"] if position in ("", "-") else [start, end if sep else start]
        if not all(re.fullmatch(r"[0-9]+|\?|-", part) for part in pair):
            return None
        out.extend(pair)
    return out[0], out[1], out[2], out[3]


def location_span(e: RowEntry) -> int | None:
    """The number of positions an entry's Location names, or ``None`` where it names no span
    of whole positions."""
    if e.location is None:
        return None
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    end = end or start
    if not (re.fullmatch(r"[0-9]+", start) and re.fullmatch(r"[0-9]+", end)):
        return None
    span = integer_of(end, "Location position") - integer_of(start, "Location position") + 1
    return span if span >= 1 else None


def longer_allele_inside_one_exon(e: RowEntry) -> bool:
    """Whether an entry shows, without its record, a span of three or more bases that lies in
    the transcript's sequence from end to end under an allele of more nucleotide letters than
    the span has bases: its Location names the span, and its two cDNA positions are numbers
    as far apart as the span's two ends, which a span that holds an intron's base does not
    have. A Location of two neighbouring positions is refused: an insertion between them is
    written with the same Location, allele and positions as two replaced bases, and only
    the record tells the two apart."""
    span = location_span(e)
    if span is None or span < 3 or not re.fullmatch(r"[ACGTN]+", e.allele) or len(e.allele) <= span:
        return False
    positions = transcript_positions(e)
    if positions is None or not (_ten_digits(positions[0]) and _ten_digits(positions[1])):
        return False
    return integer_of(positions[1], "cDNA position") - integer_of(positions[0], "cDNA position") + 1 == span


def differing_regions(reference: str, allele: str) -> list[tuple[int, int]]:
    """The regions ``_get_differing_regions`` (VariationFeatureOverlapAllele.pm:405-460)
    cuts a record into, as offsets from the span's first base: where both alleles hold more
    than one base (:418) the runs of indices at which they differ, every index past the
    shorter allele's end differing (:428-450); else one region, the span (:452-454)."""
    if len(reference) <= 1 or len(allele) <= 1:
        return [(0, len(reference) - 1)]
    runs: list[tuple[int, int]] = []
    for index in range(max(len(reference), len(allele))):
        if index < len(reference) and index < len(allele) and reference[index] == allele[index]:
            continue
        if runs and runs[-1][1] == index - 1:
            runs[-1] = (runs[-1][0], index)
        else:
            runs.append((index, index))
    return runs or [(0, len(reference) - 1)]


def intron_rows(
    regions: Sequence[tuple[int, int]], span: int, start: int, end: int, forward: bool, wide_exons: bool
) -> tuple[frozenset[str], frozenset[str]]:
    """The splice and intron terms of the reference's row and of the candidate's for a record
    whose span is the offsets ``0`` to ``span - 1``, on a transcript of the given strand whose
    one intron near the span runs from offset ``start`` to offset ``end``.

    The reference's row is ``_intron_effects`` (BaseTranscriptVariationAllele.pm:99-224) on
    ``regions``, each placed with no bound at the span's end (:117-118). The intron is read
    where the span overlaps it widened by three bases, and its two ends where the span
    overlaps the eleven bases around each (BaseTranscriptVariation.pm:1027-1029); an intron
    of thirteen bases or fewer is a frameshift intron (:1032-1035) and a region that
    overlaps it sets nothing (BaseTranscriptVariationAllele.pm:138-141, :171-174). A region
    in the intron's interior sets the intron term (:143-149) and one in a window of the six
    below its flag (:152-159, :176-204). The flag of ``splice_region`` is assigned for each
    region in turn, whether the region lies in one of its four windows, until a region has
    met the intron's first two or last two bases (:209-216; Utils/VariationEffect.pm:91-94),
    so the last region assigned decides it. The tract term is written only
    off the exons (Constants.pm:971-983), which a transcript with a frameshift intron reads
    twelve bases wider (BaseTranscriptVariation.pm:861; ``wide_exons`` says another intron
    of the transcript is one). The donor region is withheld beside the fifth base and
    ``splice_region_variant`` beside the first four terms (Utils/VariationEffect.pm:597-627).

    The candidate's row keeps a term of the six only where the span itself holds a base of
    a window of that term, the intron being the intron term's window, reads
    ``splice_region`` from any region in one of its four windows, and withholds the donor
    region and ``splice_region_variant`` by its own row's terms."""
    last = span - 1
    in_list = last >= start - 3 and end + 3 >= 0
    at_an_end = (last >= start - 3 and start + 7 >= 0) or (last >= end - 7 and end + 3 >= 0)
    frameshift = end - start <= 12
    read = [r for r in regions if not (frameshift and r[1] >= start and r[0] <= end)]

    def met(first: int, final: int) -> bool:
        return any(r[1] >= first and r[0] <= final for r in read)

    def held(first: int, final: int) -> bool:
        return last >= first and final >= 0

    if forward:
        windows = {
            "splice_donor_variant": (start, start + 1),
            "splice_acceptor_variant": (end - 1, end),
            SPLICE_DONOR_5TH_BASE_TERM: (start + 4, start + 4),
            SPLICE_DONOR_REGION_TERM: (start + 2, start + 5),
            "splice_polypyrimidine_tract_variant": (end - 16, end - 2),
        }
    else:
        windows = {
            "splice_donor_variant": (end - 1, end),
            "splice_acceptor_variant": (start, start + 1),
            SPLICE_DONOR_5TH_BASE_TERM: (end - 4, end - 4),
            SPLICE_DONOR_REGION_TERM: (end - 5, end - 2),
            "splice_polypyrimidine_tract_variant": (start + 2, start + 16),
        }
    off_the_exons = not (
        start > 0 or last > end or ((frameshift or wide_exons) and (start > -12 or last + 12 > end))
    )
    set_by_a_region: set[str] = set()
    held_by_the_span: set[str] = set()
    for term, (first, final) in windows.items():
        read_here = (in_list or at_an_end) and off_the_exons if term.startswith("splice_poly") else at_an_end
        if read_here and met(first, final):
            set_by_a_region.add(term)
        if held(first, final):
            held_by_the_span.add(term)
    if in_list and start + 2 <= end - 2 and met(start + 2, end - 2):
        set_by_a_region.add(INTRON_TERM)
    if held(start, end):
        held_by_the_span.add(INTRON_TERM)
    beside = ((start - 3, start - 1), (start + 2, start + 7), (end - 7, end - 2), (end + 1, end + 3))
    last_assigned = in_any = at_a_site = False
    if at_an_end:
        for region in read:
            in_one = any(region[1] >= first and region[0] <= final for first, final in beside)
            in_any = in_any or in_one
            at_a_site = at_a_site or (region[1] >= start and region[0] <= start + 1) or (region[1] >= end - 1 and region[0] <= end)
            if not at_a_site:
                last_assigned = in_one

    def row(terms: set[str], region: bool) -> frozenset[str]:
        if SPLICE_DONOR_5TH_BASE_TERM in terms:
            terms.discard(SPLICE_DONOR_REGION_TERM)
        if region and not terms & SPLICE_REGION_WITHHOLDING_TERMS:
            terms.add(SPLICE_REGION_TERM)
        return frozenset(terms)

    return row(set(set_by_a_region), last_assigned), row(set_by_a_region & held_by_the_span, in_any)


def intron_placements(span: int, reach: int, positions: tuple[str, str, str, str]) -> Iterator[tuple[bool, int, int]]:
    """The places ``(forward, start, end)`` of one intron beside a span of ``span`` bases that
    an entry's cDNA positions allow, as offsets from the span's first base, for regions that
    run ``reach`` bases past the span's end. A numbered cDNA position is a base of the
    transcript's sequence and any other a base of the intron; the first position is the
    span's first base on the forward strand and its last on the reverse. Two numbers as far
    apart as the span put the intron after the span, one to seven bases on, as far as one of
    its ends is still read, or before it, where an intron shorter than its windows has
    them run on past its last base; two numbers closer than the span put it inside the
    span. A gap wider than any window, and than the regions' reach, is placed once."""
    last = span - 1
    before = [-gap for gap in range(18)] + [-_FAR]
    after = [last + gap for gap in sorted({*range(42), *range(max(42, reach - 20), reach + 19)})] + [_FAR]
    numbered_first, numbered_last = _ten_digits(positions[0]), _ten_digits(positions[1])
    for forward in (True, False):
        first_in_sequence, last_in_sequence = (
            (numbered_first, numbered_last) if forward else (numbered_last, numbered_first)
        )
        if first_in_sequence and last_in_sequence:
            width = integer_of(positions[1], "cDNA position") - integer_of(positions[0], "cDNA position") + 1
            if width == span:
                for gap in range(1, 8):
                    start = last + gap
                    yield from ((forward, start, end) for end in after if end >= start)
                    yield from ((forward, start, start + length) for length in range(42))
                yield from ((forward, start, end) for start in range(-16, 0) for end in range(start, 0))
            elif 0 < width < span:
                length = span - width
                yield from ((forward, start, start + length - 1) for start in range(1, last - length + 1))
        elif first_in_sequence:
            yield from ((forward, start, end) for start in range(1, span) for end in after)
        elif last_in_sequence:
            yield from ((forward, start, end) for end in range(last) for start in before)
        else:
            yield from ((forward, start, end) for start in before for end in after)


def record_span_placement(e: RowEntry, a: RowEntry, reference: str) -> tuple[bool, int, int, bool] | None:
    """A place of the transcript's intron, ``(forward, start, end, wide_exons)`` as
    ``intron_rows`` reads it, at which the record's regions give the reference entry's
    splice and intron terms and the span gives the candidate entry's, or ``None`` where
    the entry's cDNA positions allow no such place. ``reference`` is the reference bases
    the allele stands in for; the places are those of ``intron_placements``."""
    positions = transcript_positions(e)
    if positions is None:
        return None
    read = PAST_THE_SPAN_TERMS | {SPLICE_REGION_TERM}
    expected = frozenset(consequence_set(e.consequence).split(",")) & read
    actual = frozenset(consequence_set(a.consequence).split(",")) & read
    regions = differing_regions(reference, e.allele)
    tract = "splice_polypyrimidine_tract_variant"
    for forward, start, end in intron_placements(len(reference), len(e.allele) - len(reference), positions):
        rows = intron_rows(regions, len(reference), start, end, forward, False)
        if rows == (expected, actual):
            return forward, start, end, False
        if (
            tract in rows[0]
            and tract not in expected
            and intron_rows(regions, len(reference), start, end, forward, True) == (expected, actual)
        ):
            return forward, start, end, True
    return None


# ``_get_differing_regions`` (VariationFeatureOverlapAllele.pm:405-460) compares the two
# alleles index by index to the longer allele's length (:428-429), so every index past the
# reference allele's end differs, and ``_intron_effects`` places each differing region on the
# genome from the record's start with no bound at the span's end
# (BaseTranscriptVariationAllele.pm:117-118). On a record whose alternate allele is the longer
# one, the flags of the five splice terms above and of the intron term (:143-204) are
# therefore set from positions the record does not replace. The candidate writes such a term
# only where a base of the span itself lies in a window of that term, the intron being the
# intron term's window.
#
# With the record's alleles (``RowEntry.reference``: the first column or the format's own
# allele fields, else the input records) the regions are known by position, and a pair is of
# the class where the record is no pure insertion, both alleles hold more than one base and
# the alternate more than the reference, and one intron stands beside the span at a place
# the entry's cDNA positions allow (``intron_placements``) at which the regions give the
# reference's six terms and ``splice_region_variant`` and the span gives the candidate's
# (``intron_rows``): the terms gone are of the six, each without a base of the span in its
# window, and the terms come are the donor region where the fifth base went and
# ``splice_region_variant`` where the terms that withhold it went. No input names the
# transcript's introns, so the place is read off the reference's row: where two places
# give that row and differ in what the span holds, the pair is taken at the one that
# gives the candidate's row.
#
# Without the record a row shows that no base of its span lies in any window where three
# or more positions lie in the transcript's sequence from end to end
# (``longer_allele_inside_one_exon``): there the candidate's set is the reference's without
# every term of the six, with or without ``splice_region_variant`` in the place of a term
# that withholds it. A Location of two neighbouring positions, which an insertion's row
# shares, and a span that holds an intron's base then stay differences. With or without
# the record, two rows whose cDNA or CDS positions differ and the terms the other way
# round stay differences, and so does any other differing term, but for the terms of a
# stop codon class on the same row: there that class's rule reads the rows without this
# difference (``STOP_CODON_ROW_PAIRS``) and this one the rows without that.
def splice_term_past_the_record_span_pair(e: RowEntry, a: RowEntry) -> bool:
    """Whether the reference entry ``e`` and the candidate entry ``a`` of one key are a row on
    which Ensembl VEP writes one or more of the five intron-window splice terms or the intron
    term from positions past the end of the record's span and the candidate does not: a
    Transcript entry with the same cDNA and CDS positions on both sides whose two sets of
    terms differ by terms of the six the reference alone holds and, at most, the donor
    region and ``splice_region_variant`` the candidate alone holds, beside the terms of a
    stop codon class where its rule takes the rest of the difference. With the record's
    alleles the two rows are those of a longer alternate allele that is no pure insertion
    at a place of the intron the cDNA positions allow (``record_span_placement``); without
    them the span lies in the transcript's sequence under a longer allele
    (``longer_allele_inside_one_exon``), the candidate holds none of the six, and
    ``splice_region_variant`` alone comes, where a term that withholds it went."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript":
        return False
    positions = transcript_positions(e)
    if positions is None or positions != transcript_positions(a):
        return False
    expected_terms = set(consequence_set(e.consequence).split(","))
    actual_terms = set(consequence_set(a.consequence).split(","))
    read = PAST_THE_SPAN_TERMS | {SPLICE_REGION_TERM}
    if expected_terms - read != actual_terms - read:
        # The rest of the difference is a stop codon class's: its rule reads the two rows
        # with the candidate's splice and intron terms on both, and this one the two rows
        # with the reference's other terms on both.
        beside = copy.copy(e)
        beside.consequence = ",".join(sorted((expected_terms - read) | (actual_terms & read)))
        if not any(pair(beside, a) for pair in STOP_CODON_ROW_PAIRS):
            return False
        actual_terms = (actual_terms & read) | (expected_terms - read)
        a = copy.copy(a)
        a.consequence = ",".join(sorted(actual_terms))
    gone, come = expected_terms - actual_terms, actual_terms - expected_terms
    if not gone or not gone <= PAST_THE_SPAN_TERMS or not come <= {SPLICE_REGION_TERM, SPLICE_DONOR_REGION_TERM}:
        return False
    reference = e.reference if e.reference is not None else a.reference
    if reference is None:
        if not (longer_allele_inside_one_exon(e) and longer_allele_inside_one_exon(a)):
            return False
        if actual_terms & PAST_THE_SPAN_TERMS:
            return False
        return not come or (come == {SPLICE_REGION_TERM} and bool(gone & SPLICE_REGION_WITHHOLDING_TERMS))
    if location_span(e) != len(reference) or len(reference) < 2:
        return False
    if not re.fullmatch(r"[ACGTN]+", e.allele) or len(e.allele) <= len(reference):
        return False
    return record_span_placement(e, a, reference) is not None


@dataclass
class PairRow:
    """The `RowEntry` of one row of an output."""

    key: str
    location: str | None
    allele: str
    consequence: str
    fields: dict[str, str]
    indel: bool | None = None
    insertion: bool | None = None
    reference: str | None = None


def pair_row(columns: Sequence[str]) -> PairRow:
    """The first nine columns of an output row, Uploaded_variation to CDS_position, in the
    form a row-aware pair rule reads, as the format comparator's ``parse_default`` builds a
    default-format row: the key of Location, Allele, Feature and Feature_type, and the
    allele's shape and the reference bases it stands in for from the record's own alleles
    where the first column spells them."""
    name, location, allele, _, feature, feature_type, consequence, cdna_position, cds_position = columns
    indel = insertion = reference = None
    record_alleles = uploaded_alleles(name)
    if record_alleles is not None:
        reference_allele, alternate_alleles = record_alleles
        indel = allele_is_indel(reference_allele, alternate_alleles, allele)
        insertion = reference_allele_is_empty(reference_allele, alternate_alleles)
        position = uploaded_position(name)
        if position is not None:
            reference = replaced_reference(position, reference_allele, alternate_alleles, location, allele)
    return PairRow(
        key=f"{location}|{allele}|{feature}|{feature_type}",
        location=location,
        allele=allele,
        consequence=consequence,
        fields={"cDNA_position": cdna_position, "CDS_position": cds_position},
        indel=indel,
        insertion=insertion,
        reference=reference,
    )


# The rows of the two outputs that carry the two tuples of one unmatched pair, the
# reference's first.
PairRows = tuple[Sequence[RowEntry], Sequence[RowEntry]]


@dataclass(frozen=True)
class ExclusionRule:
    """One divergence shape: its recogniser, its reporting name, whether the
    adjusted F1 removes it, and the taxonomy class it belongs to.

    ``classifies`` and ``excludes`` are independent because the two roles are
    independent. A rule that classifies assigns a pair to its reporting bucket;
    a rule that excludes removes the pair from both raw denominators. Of the rules
    that read the two consequence sets, only ``start_cooccurrence_swap`` does both.
    """

    bucket: str
    predicate: Callable[..., bool]
    classifies: bool
    excludes: bool
    taxonomy_class: str
    # A predicate that decides on the allele kind as well as the two consequence
    # sets takes the tuple's allele as its third argument; the others take two.
    allele_aware: bool = False
    # A row-aware predicate takes the reference's and vep-rs's row of the pair (`RowEntry`)
    # and matches only where the caller hands the rows over: every row that carries the
    # reference's tuple against every row that carries vep-rs's, so the rows of one key
    # stand or fall together.
    row_aware: bool = False

    def matches(self, perl_csq: str, rust_csq: str, allele: str, rows: PairRows | None = None) -> bool:
        if self.row_aware:
            if rows is None or not rows[0] or not rows[1]:
                return False
            return all(self.predicate(e, a) for e in rows[0] for a in rows[1])
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

# The pair rules of the 116.2 registry that read the two rows of an unmatched pair
# (``row_aware``). ORDER IS THE CLASSIFICATION ORDER here too: they are asked after the
# rules above, for a pair none of those classifies or excludes, and the first match takes it.
# Each differs from its reference's set in a stop codon term, or in a splice term the reference
# alone holds, which no rule above admits, and no pair satisfies two of them, so no order among
# them moves a count. The 115.2 registry holds none, so its report keeps the published shape.
REFERENCE_116_2_ROW_PAIR_RULES: tuple[ExclusionRule, ...] = (
    ExclusionRule(
        bucket="stop_codon_overlap_genomic_span",
        predicate=stop_codon_overlap_genomic_span_pair,
        classifies=True,
        excludes=True,
        taxonomy_class="stop_codon_overlap_genomic_span",
        row_aware=True,
    ),
    ExclusionRule(
        bucket="stop_retained_without_stop_codon_test",
        predicate=stop_retained_without_stop_codon_test_pair,
        classifies=True,
        excludes=True,
        taxonomy_class="stop_retained_without_stop_codon_test",
        row_aware=True,
    ),
    ExclusionRule(
        bucket="stop_retained_three_base_stop_codon_test",
        predicate=stop_retained_three_base_stop_codon_test_pair,
        classifies=True,
        excludes=True,
        taxonomy_class="stop_retained_three_base_stop_codon_test",
        row_aware=True,
    ),
    ExclusionRule(
        bucket="stop_codon_removed_without_stop_codon_test",
        predicate=stop_codon_removed_without_stop_codon_test_pair,
        classifies=True,
        excludes=True,
        taxonomy_class="stop_codon_removed_without_stop_codon_test",
        row_aware=True,
    ),
    ExclusionRule(
        bucket="stop_codon_overlap_insertion_range",
        predicate=stop_codon_overlap_insertion_range_pair,
        classifies=True,
        excludes=True,
        taxonomy_class="stop_codon_overlap_insertion_range",
        row_aware=True,
    ),
    ExclusionRule(
        bucket="splice_term_past_the_record_span",
        predicate=splice_term_past_the_record_span_pair,
        classifies=True,
        excludes=True,
        taxonomy_class="splice_term_past_the_record_span",
        row_aware=True,
    ),
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
    "stop_codon_overlap_genomic_span": (
        "A Transcript row of an insertion or deletion with both cDNA positions and a CDS "
        "start on which the reference writes coding_sequence_variant and vep-rs "
        "stop_retained_variant, both beside splice_region_variant: the reference's genomic "
        "stop codon overlap reads the last three genomic positions of the coding region and "
        "misses a stop base that ends an exon."
    ),
    "stop_retained_without_stop_codon_test": (
        "A Transcript row of an insertion or deletion whose stop codon test the reference "
        "does not make (an allele holding a letter that names no base, or a span without a "
        "cDNA position at one end or a CDS start) and writes stop_retained_variant for; "
        "vep-rs writes the same set without the term, with nothing, coding_sequence_variant "
        "or frameshift_variant in its place."
    ),
    "stop_retained_three_base_stop_codon_test": (
        "A Transcript row of an insertion or deletion whose stop codon test is made, on which "
        "the reference writes stop_retained_variant from the three bases at the old stop "
        "position alone; vep-rs writes coding_sequence_variant in its place, or, where the "
        "reference writes it alone or beside one of inframe_insertion, inframe_deletion and "
        "protein_altering_variant, frameshift_variant with or without stop_gained in the place "
        "of those terms."
    ),
    "stop_codon_removed_without_stop_codon_test": (
        "A Transcript row of an insertion or deletion without a cDNA position at one end or a "
        "CDS start, on which the reference writes stop_retained_variant with no stop codon "
        "test made and vep-rs stop_lost for the stop codon the record removes."
    ),
    "stop_codon_overlap_insertion_range": (
        "A Transcript row of an insertion with both cDNA positions and a CDS start on which "
        "vep-rs writes stop_retained_variant, with or without inframe_insertion, for the stop "
        "codon the insertion keeps, and the reference coding_sequence_variant, or "
        "frameshift_variant with or without stop_gained: the reference's genomic stop codon "
        "overlap moves an insertion's start by the inserted length. Under an identifier a row "
        "is read as an insertion's by a Location of two neighbouring positions and an Allele "
        "of another length than two, so two bases replaced by a sequence of another length "
        "are set aside with the insertions."
    ),
    "splice_term_past_the_record_span": (
        "On a record whose alternate allele is longer than its reference allele, Ensembl VEP "
        "writes a splice donor, splice acceptor, donor fifth base, donor region, "
        "polypyrimidine tract or intron term from positions past the end of the record's span; "
        "where no base of the span lies in a window of that term, vep-rs does not write it. "
        "For the intron term the window is the intron. With the record's alleles, from the "
        "first column or the input records, a row is set aside where the record is no pure "
        "insertion and one intron stands beside the span, at a place the row's cDNA positions "
        "allow, at which the record's regions give the reference's terms and the span gives "
        "vep-rs's, the donor region and splice_region_variant coming where the terms that "
        "withhold them went. Without the record a row is set aside where three or more "
        "positions lie in the transcript's sequence and vep-rs's row holds none of the six "
        "terms; a Location of two neighbouring positions, which an insertion's row shares, "
        "and a span that holds an intron's base then stay differences. A row that differs by "
        "the terms of a stop codon class as well is set aside where that class's rule takes "
        "the rest of the difference."
    ),
}
PAIR_RULE_CITATIONS: dict[str, str] = {
    "splice_family_swap": "BaseTranscriptVariationAllele.pm:215; VariationEffect.pm:87-110",
    "start_cooccurrence_swap": "VariationEffect.pm:851-963, 1028-1075",
    "stop_codon_overlap_genomic_span": "VariationEffect.pm:1389-1431, 1616-1649",
    "stop_retained_without_stop_codon_test": "VariationEffect.pm:1290-1324, 1488-1539",
    "stop_retained_three_base_stop_codon_test": "VariationEffect.pm:1318, 1488-1539",
    "stop_codon_removed_without_stop_codon_test": "VariationEffect.pm:1234-1288, 1290-1324, 1508",
    "stop_codon_overlap_insertion_range": "VariationEffect.pm:1389-1431, 1549",
    "splice_term_past_the_record_span": (
        "VariationFeatureOverlapAllele.pm:405-460; BaseTranscriptVariationAllele.pm:117-118, 143-204"
    ),
}


_BRACKET_MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", re.ASCII)
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


def is_breakend_allele(allele: str) -> bool:
    """Whether the allele column names a breakend's own end or its mate."""
    return is_own_breakend_allele(allele) or bracket_mate(allele) is not None


def dot_breakend_points(span: tuple[str, int, int]) -> list[str]:
    """The Locations at which vep-rs can write the point `c:(P+1)` of a dot-notation breakend
    without SVTYPE whose reference row spans ``span`` (chromosome, low, high): the span's start
    and the base after it. The reference takes the sequence rules for the coordinates
    (BaseVCF4.pm:264-275, 297-317): the start is P + 1 when the ALT's length differs from the
    REF's and P when it does not, the end P + len(REF) - 1, the two written low first. A REF of
    one base (`N N.`) gives `c:P-(P+1)` and the point is the base after the start; a REF of
    three or more bases beside an ALT of another length (`NAC .N`) gives `c:(P+1)-(P+len(REF)-1)`
    and the point is the start; an ALT of the REF's length (`NA N.`, `NAC NA.`, `NACG NAC.`)
    gives `c:P-(P+len(REF)-1)` and the point is the base after the start."""
    chrom, start, _end = span
    return [f"{chrom}:{start}", f"{chrom}:{start + 1}"]


def reference_mate_chromosome(allele: str) -> str | None:
    """The mate chromosome of a bracket-notation breakend allele as the reference resolves it:
    a ``chr`` prefix as written is removed and any other spelling is kept, so `N[chr21:...[`
    names ``21`` and `N[Chr21:...[` names ``Chr21``, a name the cache files no transcript
    under (StructuralVariationOverlap.pm:86-102 compares it with the transcript's); None for
    any other allele."""
    m = _BRACKET_MATE_RE.search(allele)
    return m.group(1).removeprefix("chr") if m else None


@dataclass(frozen=True)
class OneSidedRule:
    """A divergence shape recognised on rows one engine wrote and the other did not.

    Unlike an `ExclusionRule`, which judges a pair of consequence sets on one
    transcript, a one-sided rule judges a record: what the reference wrote for it
    (nothing, or rows under some alleles) against what vep-rs wrote. ``excludes`` says
    whether this comparator's adjusted F1 removes the rows it recognises;
    ``excludes_with_cache`` whether the structural-variant comparator, which reads the input
    records and a cache authority, sets the rows aside. A rule that only classifies reports a
    count so the shape is visible without being credited to vep-rs.
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
BND_BARE_REFERENCE_ALLELE_ROW_BUCKET = "bnd_bare_reference_allele_row"
BND_MATE_32_BIT_BUFFER_BUCKET = "bnd_mate_32bit_buffer_intergenic"
GVCF_VARIANT_SITE_DROPPED_BUCKET = "gvcf_variant_site_dropped"
DOT_BREAKEND_PAIR_BUCKET = "dot_breakend_without_svtype_pair"
MISSING_ALT_ANNOTATED_BUCKET = "missing_alt_annotated"
NON_INTEGER_END_BUCKET = "non_integer_end_annotated"
CNV_TR_EMPTY_EXPANSION_BUCKET = "cnv_tr_empty_expansion_dropped"
STAR_ALONE_SVTYPE_ANNOTATED_BUCKET = "star_alone_svtype_annotated"
NON_NUMERIC_POS_NUMIFIED_BUCKET = "non_numeric_pos_numified"
CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET = "cnv_tr_rn_zero_allele_dropped"
# The smallest coordinate at which Ensembl VEP 116.2 annotates every record of the input
# buffer as intergenic_variant: the coordinate plus the default --distance flank reaches
# 2^31, so the bound moves with the flank.
BND_MATE_32_BIT_LIMIT = 2**31 - 5000
# A coordinate with fewer digits than the bound cannot reach it.
_BND_MATE_32_BIT_LIMIT_DIGITS = len(str(BND_MATE_32_BIT_LIMIT))
# The largest coordinate vep-rs holds; a POS past it fails its line in every input format,
# where the reference numifies the text and annotates the record.
VEP_RS_COORDINATE_CEILING = 2**62 - 1
# The POS spellings vep-rs reads besides a decimal integer: what Rust's `f64` parser accepts,
# ASCII digits with an optional fraction and exponent and an optional `+`, blanks around.
_VEP_RS_POS_NUMBER = re.compile(r"\+?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?", re.ASCII)


def vep_rs_pos(text: str) -> int | None:
    """The POS column as vep-rs's parsers read it: an unpadded decimal integer up to
    `VEP_RS_COORDINATE_CEILING` exactly, or another spelling of an exact non-negative integer
    below the ceiling that Rust's `f64` parser accepts, blanks around it ignored (`3e7`,
    `30425650.0`, `30425650.`, `+30425650`, ` 30425650 `), read as the double it names
    (` 72057594037927937 ` is 72057594037927936); None where vep-rs fails the line (`.`, `abc`,
    `-5`, `0x1D04A1A`, `30_425_650`, fullwidth digits, `inf`, `1e400`, `30425650.5`)."""
    stripped = text.strip(" \t\n\r\x0c")
    if _VEP_RS_POS_NUMBER.fullmatch(stripped) is None:
        return None
    if text.removeprefix("+").isdigit():
        value = int(text)
        return value if value <= VEP_RS_COORDINATE_CEILING else None
    # `VEP_RS_COORDINATE_CEILING` rounds up to 2^62 as a double, so the bound is strict.
    number = float(stripped)
    return int(number) if number.is_integer() and number < 2.0**62 else None


# The two lost-row breakend shapes are counted, never set aside, by this comparator:
# deciding which of the rows the reference would have written needs the cache's
# transcript-to-chromosome authority, which this comparator does not read (the
# structural-variant comparator sets them aside under it).
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
        perl_citation="Parser.pm:491-501, 821-827; Parser/VCF.pm:477-481, 555-558, 575",
    ),
    OneSidedRule(
        bucket=BND_OWN_ROWS_LOST_BUCKET,
        definition=(
            "A bracket breakend whose mate is on another chromosome loses the rows of its "
            "own end when a mate-chromosome transcript, loaded for an earlier record of the "
            "same buffer, reaches it first and names its slice: the reference keeps the mate "
            "rows and writes nothing under the own-end allele, an order-dependent result."
        ),
        excludes=False,
        excludes_with_cache=True,
        taxonomy_class="bnd_own_chromosome_rows_lost",
        perl_citation=(
            "InputBuffer.pm:191-207, 223-234; AnnotationSource.pm:235-258; "
            "AnnotationType/Transcript.pm:122-143"
        ),
    ),
    OneSidedRule(
        bucket=BND_SYNONYM_MATE_ROWS_LOST_BUCKET,
        definition=(
            "A bracket breakend whose mate's chromosome is spelled with a prefix the "
            "mate-allele test does not resolve, Chr21 for 21, loses its mate rows in the "
            "reference and keeps them in vep-rs; a mate named by an accession on another "
            "chromosome loses them in the reference alone too, vep-rs resolving the accession "
            "through the synonym table, and one on the record's own chromosome keeps them in "
            "both engines."
        ),
        excludes=False,
        excludes_with_cache=True,
        taxonomy_class="bnd_synonym_mate_rows_lost",
        perl_citation=(
            "StructuralVariationOverlap.pm:86-102, 144-160; Utils/VariationEffect.pm:111-124; "
            "InputBuffer.pm:348-357; AnnotationType/Transcript.pm:134"
        ),
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
    OneSidedRule(
        bucket=BND_BARE_REFERENCE_ALLELE_ROW_BUCKET,
        definition=(
            "A breakpoint record that names no mate (a bare <BND>, or one with SVLEN, POS2 or "
            "an unreadable END2 and no bracket) is held as N/<BND> and, with no breakend "
            "parsed, one allele is built per slash-separated part, the reference base N "
            "included, so the reference writes every row twice under chromosome_breakpoint, "
            "ALLELE_NUM 1 for the reference allele and 2 for the ALT; a reference allele has "
            "no variant effect, so vep-rs writes the ALT row once and the reference-allele "
            "rows are set aside: each duplicates the ALT row's tuple, so the unique-tuple "
            "key stream drops it before the raw and the adjusted figure, and both "
            "comparators count the rows so set aside per record."
        ),
        excludes=True,
        excludes_with_cache=True,
        taxonomy_class="bnd_bare_reference_allele_row",
        perl_citation="Parser/VCF.pm:499; StructuralVariationOverlap.pm:64-77; OutputFactory.pm:2029",
    ),
    OneSidedRule(
        bucket=BND_MATE_32_BIT_BUFFER_BUCKET,
        definition=(
            "A coordinate of the reference's input buffer that, plus the --distance flank, "
            "reaches 2^31 (a breakend's mate position, a span's end or the record's own "
            "position) turns every record of the buffer intergenic, that record and the "
            "records loaded with it (the first such position is 2147478648 under the "
            "default flank of 5000 and 2147483647 under --distance 1); the coordinates "
            "enter the buffer's bounds with the flank added when the cache regions are "
            "chosen. vep-rs annotates the buffer's records against their transcripts and a "
            "mate at its coordinate. Counted on the reference's rows at the record's "
            "Location and never set aside, since the records it reaches beyond that one "
            "depend on the buffer boundaries."
        ),
        excludes=False,
        taxonomy_class="bnd_mate_32bit_buffer_intergenic",
        perl_citation="AnnotationSource.pm:190, 251-254, 306-309; InputBuffer.pm:386-390",
    ),
    OneSidedRule(
        bucket=GVCF_VARIANT_SITE_DROPPED_BUCKET,
        definition=(
            "A gVCF variant site lists <NON_REF> beside its real allele (`A G,<NON_REF>`); "
            "the <NON_REF> token routes the line to the structural path, where "
            "`G/<NON_REF>` has no Sequence Ontology term, and the reference drops the site "
            "whole, the real allele with it. vep-rs annotates the real allele and reads "
            "<NON_REF> as the absence of a variant. vep-rs's rows at the site are set aside: "
            "recognised here by the record's own name (`CHROM_POS_REF/ALT/<NON_REF>`) at a "
            "Location the reference wrote nothing for, and in the structural-variant "
            "comparator from the input record."
        ),
        excludes=True,
        excludes_with_cache=True,
        taxonomy_class="gvcf_variant_site_dropped",
        perl_citation="Parser/VCF.pm:237-245, 469-481; Parser.pm:663-728",
    ),
    OneSidedRule(
        bucket=DOT_BREAKEND_PAIR_BUCKET,
        definition=(
            "A dot-notation breakend without INFO/SVTYPE (`N.`, `.N`) takes the sequence "
            "rules for its coordinates: the start is POS + 1 when the ALT's length differs "
            "from the REF's and POS when it does not, the end is POS + len(REF) - 1, so the "
            "reference writes the pair `POS-(POS+1)` with the positions of an insertion for "
            "a one-base REF, the span `(POS+1)-(POS+len(REF)-1)` for a longer REF beside an "
            "ALT of another length and `POS-(POS+len(REF)-1)` beside an ALT of its own "
            "length; with `SVTYPE=BND` it writes the single base at POS + 1, as vep-rs "
            "writes it for every REF. The reference's span row and vep-rs's point row at "
            "the span's start or the base after it on the same feature are set aside "
            "together."
        ),
        excludes=True,
        excludes_with_cache=True,
        taxonomy_class="dot_breakend_without_svtype_pair",
        perl_citation="BaseVCF4.pm:255-275, 297-317",
    ),
    OneSidedRule(
        bucket=MISSING_ALT_ANNOTATED_BUCKET,
        definition=(
            "A record whose ALT column is the missing value `.` is routed to the structural "
            "path by INFO/SVTYPE and annotated as an allele: a `chromosome_breakpoint` for "
            "`SVTYPE=BND`, a deletion over END for `SVTYPE=DEL`; without SVTYPE the "
            "reference writes nothing for it, as vep-rs writes nothing in both cases. The "
            "reference's rows are set aside by the structural-variant comparator, which "
            "reads the input; the SNP/indel comparator cannot tell them from a record with "
            "an allele."
        ),
        excludes=False,
        excludes_with_cache=True,
        taxonomy_class="missing_alt_annotated",
        perl_citation="BaseVCF4.pm:486-490; Parser/VCF.pm:237-245; Parser.pm:710-711",
    ),
    OneSidedRule(
        bucket=NON_INTEGER_END_BUCKET,
        definition=(
            "An INFO/END, SVLEN or POS that is not an integer (`END=30425949.5`, "
            "`SVLEN=300.7`, `SVLEN=NaN`, POS `30425949.0` on the sequence path) is taken as "
            "written, passes `looks_like_number`, and the record is annotated with the text "
            "in its Location (`21:30425651-30425950.7`, `21:30425651-NaN`, `21:30425949.0`) "
            "and a fractional DISTANCE; vep-rs refuses an END that is not a coordinate and "
            "reads a SVLEN or POS at its numified value (`21:30425651-30425950`, the pair "
            "`21:30425650-30425651`, `21:30425949`). The reference's rows, recognised by the "
            "Location, are set aside in both comparators, and vep-rs's rows of the same "
            "record by the structural-variant comparator, which reads the input."
        ),
        excludes=True,
        excludes_with_cache=True,
        taxonomy_class="non_integer_end_annotated",
        perl_citation="BaseVCF4.pm:263, 302-308; Parser.pm:475",
    ),
    OneSidedRule(
        bucket=CNV_TR_EMPTY_EXPANSION_BUCKET,
        definition=(
            "With a FASTA, a <CNV:TR> whose expansion yields no alternate allele, or one "
            "equal to the reference run (`RN=.`, a RUC Perl reads as true that gives zero or "
            "fewer repeats, a RUC times unit length equal to SVLEN), is expanded to an allele "
            "string with no alternate and the reference writes nothing for it, without a "
            "warning; vep-rs with a FASTA writes nothing for it either, and without one "
            "writes the symbolic rows over the run, as the reference does without one. "
            "vep-rs's rows are set aside by the structural-variant comparator when the input "
            "record has the shape and the reference wrote nothing at or beside the run; the "
            "SNP/indel comparator does not judge it."
        ),
        excludes=False,
        excludes_with_cache=True,
        taxonomy_class="cnv_tr_empty_expansion_dropped",
        perl_citation="Parser/VCF.pm:382-431, 516-543; Parser.pm:553-592",
    ),
    OneSidedRule(
        bucket=STAR_ALONE_SVTYPE_ANNOTATED_BUCKET,
        definition=(
            "A record whose every ALT is `*` under a true INFO/SVTYPE is routed to the "
            "structural path and annotated from INFO alone, a deletion over END for "
            "`SVTYPE=DEL` or a `chromosome_breakpoint` at POS + 1 for `SVTYPE=BND`, although "
            "`*` names no allele of the record; vep-rs writes no row for it. The reference's "
            "rows are set aside by the structural-variant comparator, which reads the input; "
            "the SNP/indel comparator cannot tell them from a record with an allele."
        ),
        excludes=False,
        excludes_with_cache=True,
        taxonomy_class="star_alone_svtype_annotated",
        perl_citation="Parser/VCF.pm:237-245, 469-475; StructuralVariationOverlap.pm:64-77",
    ),
    OneSidedRule(
        bucket=NON_NUMERIC_POS_NUMIFIED_BUCKET,
        definition=(
            "A POS that is not a non-negative number (`.`, empty, `0x1D04A1D`, fullwidth "
            "digits, `30_425_949`, `-5`) is numified by the structural path's increment, `.` "
            "to 0 and the record annotated from base 1, `30_425_949` from 31, `-5` from -4, "
            "where vep-rs refuses the line; the reference's rows at the record's start are "
            "set aside by the structural-variant comparator from the input, and by the "
            "SNP/indel comparator when the Location's first coordinate is negative, the one "
            "spelling no written POS reaches."
        ),
        excludes=True,
        excludes_with_cache=True,
        taxonomy_class="non_numeric_pos_numified",
        perl_citation="BaseVCF4.pm:255-264; Parser.pm:475-488",
    ),
    OneSidedRule(
        bucket=CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET,
        definition=(
            "With a FASTA, a two-ALT <CNV:TR> of which one allele has `RN=0` loses an "
            "allele: the empty expansion is an empty part of the allele string, which is "
            "not read as multi-allelic, so the first allele alone is kept when the empty "
            "part leads and the trailing empty part is dropped by the split; vep-rs writes "
            "both alleles when the empty part leads, and the structural-variant comparator "
            "sets aside its rows of the "
            "dropped allele: the symbolic rows told by their direction term when the real "
            "allele is a gain, or unchanged beside a leading empty part (vep-rs reads every "
            "symbolic <CNV:TR> allele as a gain), and the literal rows a run with a FASTA "
            "writes for the real allele behind a leading empty part."
        ),
        excludes=False,
        excludes_with_cache=True,
        taxonomy_class="cnv_tr_rn_zero_allele_dropped",
        perl_citation="Parser/VCF.pm:411, 428; Parser.pm:881, 978",
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
    # The pair rules that read the two rows of the pair in the two outputs, asked after
    # ``pair_rules``.
    row_pair_rules: tuple[ExclusionRule, ...] = ()


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
        row_pair_rules=REFERENCE_116_2_ROW_PAIR_RULES,
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


def report_class_definitions(reference_release: str) -> dict[str, dict[str, object]]:
    """Every class of a reference release's registry in the report's order: the pair rules
    of `divergence_class_definitions`, the row-aware pair rules after them in
    classification order, then the one-sided rules."""
    registry = REFERENCE_REGISTRIES[reference_release]
    classes = list(divergence_class_definitions(reference_release).items())
    row_classes = [
        (
            rule.bucket,
            {
                "definition": PAIR_RULE_DEFINITIONS[rule.bucket],
                "excludes": rule.excludes,
                "taxonomy_class": rule.taxonomy_class,
                "perl_citation": PAIR_RULE_CITATIONS[rule.taxonomy_class],
            },
        )
        for rule in registry.row_pair_rules
    ]
    at = len(registry.pair_rules)
    return dict(classes[:at] + row_classes + classes[at:])


class ReferenceRecordScope:
    """What the 116.2 one-sided pass needs from the sorted merge, collected in the one
    read of both key streams.

    ``absent``: the (file, Location) pairs at which vep-rs wrote a structural row and
    the reference wrote nothing, less those `record_reference_named` takes out because a
    row there carries the name of a record the reference wrote a row for
    (`reference_named_locations`): a record the reference dropped as `vep_skip`, since
    116.2 writes at least one row (an intergenic one when nothing is in range) for every
    record it keeps and every route to `vep_skip` is structural; `tandem_repeat` is
    outside the test (`_TANDEM_REPEAT_ALLELE`). Every vep-rs row at an absent Location is
    then set aside whatever its allele, a sequence allele beside the symbolic one of a
    mixed ALT list included, because the reference skips the line whole. This is the
    output-only rule: the ``absent`` test reads no input VCF, and this comparator opens
    one only under `--input-dir` or `--input-vcf`, to recognise gVCF variant sites by
    their shape (``gvcf_sites`` below). The structural-variant comparator reads the input
    and names the record from it (`skipped_record_rows`: the reference's own skip routes
    over the ALT and INFO, confirmed by the reference's absence under the record's own
    Location start and alleles), so a kept record sharing a start or a Location string
    with a skipped one is told apart there. The merge visits keys in sorted order and
    every key of one Location shares the prefix ``Location<TAB>``, so a Location's rows
    from both streams are contiguous and the group closes when the Location changes.

    ``breakends``: for every (file, Location) at which either engine wrote a breakend
    allele, the alleles each side wrote and the mate chromosomes its bracket alleles
    name, for the two breakend classes.

    ``bare_reference_rows``: per (file, Location), the reference-allele rows of a bare
    breakpoint record the unique-tuple key stream dropped, counted on the reference's
    raw output by `count_bare_breakpoint_reference_rows` and handed in through
    `record_bare_reference_rows`, since the merge never sees a duplicate row.

    ``mate_overflow``: the (file, Location) pairs at which the reference wrote a bracket
    allele whose mate position is at least `BND_MATE_32_BIT_LIMIT`, or a Location whose
    start or end is.

    ``reference_only``: the (file, Location) pairs at which the reference wrote and vep-rs
    nothing, under any allele, read by the class that recognises a numified POS from
    the reference's coordinate alone.

    ``gvcf_sites``: the (file, Location) pairs of vep-rs rows at a gVCF variant site, read
    from the input VCF's shape (`gvcf_site_locations_from_vcf`) and from the generated name
    of an ID-less record in vep-rs's output (`gvcf_site_locations`), handed in through
    `record_gvcf_sites`; such a Location the reference wrote nothing for is
    ``gvcf_absent`` rather than ``absent``, whatever its allele.
    """

    def __init__(self) -> None:
        self.absent: set[tuple[str, str]] = set()
        self.gvcf_absent: set[tuple[str, str]] = set()
        self.breakends: dict[tuple[str, str], tuple[set[str], set[str], set[str]]] = {}
        self.bare_reference_rows: dict[tuple[str, str], int] = {}
        self.mate_overflow: set[tuple[str, str]] = set()
        self.gvcf_sites: set[tuple[str, str]] = set()
        self.reference_only: set[tuple[str, str]] = set()
        self._group: tuple[str, str] | None = None
        self._reference_wrote = False
        self._rust_wrote = False
        self._rust_wrote_structural = False

    def observe(self, file_name: str, key_line: str, perl: bool, rust: bool) -> None:
        loc, allele = key_line.split("\t", 2)[:2]
        group = (file_name, loc)
        if group != self._group:
            self.close()
            self._group = group
        if perl:
            self._reference_wrote = True
            coords = loc.partition(":")[2]
            if (
                len(coords) >= _BND_MATE_32_BIT_LIMIT_DIGITS
                and max(map(len, coords.split("-"))) >= _BND_MATE_32_BIT_LIMIT_DIGITS
            ):
                span = _location_span(loc)
                if span is not None and max(span[1], span[2]) >= BND_MATE_32_BIT_LIMIT:
                    self.mate_overflow.add(group)
        if rust:
            self._rust_wrote = True
        if rust and (
            group in self.gvcf_sites
            or (is_structural_allele(allele) and allele != _TANDEM_REPEAT_ALLELE)
        ):
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
                if perl and mate[1] >= BND_MATE_32_BIT_LIMIT:
                    self.mate_overflow.add(group)

    def close(self) -> None:
        """Close the open Location group; called by the merge when its file ends."""
        if self._group is not None and self._rust_wrote_structural and not self._reference_wrote:
            if self._group in self.gvcf_sites:
                self.gvcf_absent.add(self._group)
            else:
                self.absent.add(self._group)
        if self._group is not None and self._reference_wrote and not self._rust_wrote:
            self.reference_only.add(self._group)
        self._group = None
        self._reference_wrote = False
        self._rust_wrote = False
        self._rust_wrote_structural = False

    def record_bare_reference_rows(self, file_name: str, rows_by_location: Counter[str]) -> None:
        """Record one file's reference-allele rows of bare breakpoints, per Location."""
        for loc, rows in rows_by_location.items():
            if rows:
                self.bare_reference_rows[(file_name, loc)] = rows

    def record_gvcf_sites(self, file_name: str, locations: set[str]) -> None:
        """Record one file's gVCF variant-site Locations, read from vep-rs's output."""
        self.gvcf_sites.update((file_name, loc) for loc in locations)

    def absent_locations(self, file_name: str) -> set[str]:
        """One file's Locations in ``absent``."""
        return {loc for name, loc in self.absent if name == file_name}

    def record_reference_named(self, file_name: str, locations: Iterable[str]) -> None:
        """Take out of ``absent`` one file's Locations at which a vep-rs row carries the
        name of a record the reference wrote a row for (`reference_named_locations`)."""
        self.absent.difference_update((file_name, loc) for loc in locations)


# A first column that names no input record.
_NO_RECORD_NAMES = frozenset({"", ".", "-"})


def iter_record_names(
    path: Path,
    canonical_contigs: frozenset[str] | None = None,
    names: frozenset[str] | None = None,
) -> Iterator[tuple[str, str]]:
    """Yield (record name, Location) for every row `iter_vep_key_lines` keys, under the same
    row filter: the name is the first column, the input record's ID or the name the engine
    generates for a record without one, and the Location is read through
    `canonical_location`. With ``names``, only the rows whose name is in it."""
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t", 6)
            if len(fields) < 7:
                continue
            name = fields[0].strip()
            if names is not None and name not in names:
                continue
            location = canonical_location(fields[1].strip())
            if not (location and fields[2].strip() and fields[4].strip()):
                continue
            if canonical_contigs is not None and location.split(":", 1)[0] not in canonical_contigs:
                continue
            yield name, location


def reference_named_locations(
    perl_path: Path,
    rust_path: Path,
    locations: set[str],
    canonical_contigs: frozenset[str] | None = None,
) -> set[str]:
    """Of ``locations``, each a Location of one output pair at which vep-rs wrote and the
    reference nothing, those that are not a record the reference dropped: a vep-rs row there
    carries the name of a record the reference wrote a row for, or names no record.

    Both engines print the input record's name in the first column, so the reference's
    output tells a record it dropped, whose name is on none of its rows, from a record it
    wrote at another Location. An engine that prints a structural record from another start
    than the reference's (POS where the reference prints POS + 1) has no reference row at
    its Location for a record the reference kept, and its rows stay charged. A record
    without an ID carries a generated name, which Ensembl VEP writes from the record's own
    line (`CHROM_POS_REF/ALT`, the POS as the input gives it), so the name tells such a
    record apart for an engine that generates the same one. Two records at one Location
    stand or fall together, as the Location does in ``absent``. Neither output is read
    when ``locations`` is empty.
    """
    if not locations:
        return set()
    names_at: dict[str, set[str]] = defaultdict(set)
    for name, location in iter_record_names(rust_path, canonical_contigs):
        if location in locations:
            names_at[location].add(name)
    wanted = frozenset(set().union(*names_at.values()) - _NO_RECORD_NAMES)
    written = (
        {name for name, _ in iter_record_names(perl_path, canonical_contigs, wanted)}
        if wanted
        else set()
    )
    return {
        location
        for location, names in names_at.items()
        if names & written or names & _NO_RECORD_NAMES
    }


def iter_pair_rows(path: Path, wanted: Collection[str]) -> Iterator[tuple[str, PairRow | None]]:
    """Yield (key line, row) for every row of ``path`` whose key line, as
    `iter_vep_key_lines` writes it under the same row filter, is in ``wanted``: the row as
    `pair_row` reads its first nine columns, or ``None`` for a row of fewer columns, which
    carries no position column."""
    locations = {key_line.split("\t", 1)[0] for key_line in wanted}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 7:
                continue
            location = canonical_location(fields[1].strip())
            if location not in locations:
                continue
            key_line = "\t".join(
                (
                    location,
                    fields[2].strip(),
                    fields[4].strip(),
                    fields[5].strip(),
                    normalize_consequence_set(fields[6].strip()),
                )
            )
            if key_line in wanted:
                yield key_line, pair_row([f.strip() for f in fields[:9]]) if len(fields) >= 9 else None


def fetch_pair_rows(
    discordant_rows: Sequence[dict[str, str]],
    pairs: Iterable[tuple[int, int]],
    perl_dir: Path,
    rust_dir: Path,
) -> dict[int, list[PairRow]]:
    """The rows of the two outputs behind unmatched pairs, for the row-aware pair rules.

    ``pairs`` holds, per pair, the index in ``discordant_rows`` (the rows of
    ``discordant.tsv``) of the reference's tuple and of vep-rs's. Each output file that
    carries one of them is read once more, as `reference_named_locations` reads it, and
    the result maps a row's index to the distinct rows of its output that carry its tuple.
    An index is left out when one of those rows has fewer than nine columns, so no
    row-aware rule reads a pair with such a row. No output is read when ``pairs`` is
    empty, and nothing but the two outputs is read."""
    wanted: dict[Path, dict[str, int]] = defaultdict(dict)
    for perl_index, rust_index in pairs:
        for index, directory in ((perl_index, perl_dir), (rust_index, rust_dir)):
            row = discordant_rows[index]
            key_line = "\t".join(
                (row["location"], row["allele"], row["feature"], row.get("feature_type", ""), row["consequence_set"])
            )
            wanted[directory / row["file_name"]][key_line] = index
    found: dict[int, list[PairRow]] = {}
    short: set[int] = set()
    for path, index_of in wanted.items():
        for key_line, row in iter_pair_rows(path, index_of):
            index = index_of[key_line]
            if row is None:
                short.add(index)
            elif row not in found.setdefault(index, []):
                found[index].append(row)
    return {index: rows for index, rows in found.items() if index not in short}


def read_record_references(path: Path, rows: Iterable[PairRow]) -> None:
    """Give each of ``rows`` that carries no reference bases (``PairRow.reference``) the
    bases its allele stands in for, read by position from the records of the input VCF
    ``path`` (plain or gzipped): a record is asked about a row whose Location starts on one
    of the record's own positions on its chromosome, a ``chr`` prefix set aside
    (`replaced_reference`). A row two records read differently keeps none, and the file is
    not opened when no row asks."""
    asking: dict[tuple[str, int], list[PairRow]] = defaultdict(list)
    for row in rows:
        if row.reference is None and row.location is not None and ":" in row.location:
            chrom, _, coords = canonical_location(row.location).partition(":")
            start = coords.partition("-")[0]
            if _ten_digits(start):
                asking[(_strip_chr(chrom), integer_of(start, "Location position"))].append(row)
    if not asking:
        return
    chroms = {chrom for chrom, _ in asking}
    read: dict[int, set[str]] = defaultdict(set)
    opener = gzip.open if path.suffix in (".gz", ".bgz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t", 5)
            chrom = _strip_chr(fields[0])
            if len(fields) < 5 or chrom not in chroms:
                continue
            position = vep_rs_pos(fields[1])
            if position is None:
                continue
            reference, alternates = fields[3], perl_split(fields[4], ",")
            # An insertion's Location starts on the base before the inserted bases, the
            # record's last where the shared bases are dropped.
            for start in range(position, position + max(len(reference), 1)):
                for row in asking.get((chrom, start), ()):
                    found = replaced_reference(position, reference, alternates, row.location, row.allele)
                    if found is not None:
                        read[id(row)].add(found)
    for group in asking.values():
        for row in group:
            if len(read[id(row)]) == 1:
                row.reference = next(iter(read[id(row)]))


def _strip_chr(chrom: str) -> str:
    """The reference's own normalisation of a breakend's chromosome name at the slice test
    (AnnotationType/Transcript.pm:134): the ``chr`` prefix removed and nothing else, so
    ``M`` stays ``M`` and is not the cache's ``MT``."""
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


_BARE_BREAKPOINT_ALLELE = "chromosome_breakpoint"


def count_bare_breakpoint_reference_rows(
    path: Path, canonical_contigs: frozenset[str] | None = None
) -> Counter[str]:
    """The reference's duplicate rows under `chromosome_breakpoint`, per Location.

    A breakpoint record that names no mate is held as `N/<BND>` and, with no breakend
    parsed, gets one allele per slash-separated part, the reference base `N` included
    (StructuralVariationOverlap.pm:64-77), so every row is written twice, the second
    for the reference allele, under the same key. The second and later rows of a key are
    counted here, read through `iter_vep_key_lines` so that they are the rows the
    unique-tuple key stream drops; a row that function skips, a contig outside
    ``canonical_contigs`` included, counts nothing, as does a missing file.
    """
    rows_by_location: Counter[str] = Counter()
    if not path.exists():
        return rows_by_location
    seen: set[str] = set()
    for key in iter_vep_key_lines(path, canonical_contigs=canonical_contigs):
        location, allele, _ = key.split("\t", 2)
        if allele != _BARE_BREAKPOINT_ALLELE:
            continue
        if key in seen:
            rows_by_location[location] += 1
        else:
            seen.add(key)
    return rows_by_location


# `CHROM_POS_REF/ALT...`, the name VEP gives an ID-less record; the chromosome may carry
# underscores, the position is the first `_digits_` group after it.
_GENERATED_NAME_RE = re.compile(r"^.+?_[0-9]+_(.+)$")


def gvcf_site_locations(path: Path, canonical_contigs: frozenset[str] | None = None) -> set[str]:
    """The Locations of the rows in vep-rs's output whose record is a gVCF variant site:
    an ID-less record is named `CHROM_POS_REF/ALT...` from its own line, so a name whose
    ALT list carries `<NON_REF>` beside another allele is the site (a record with an ID is
    named by the ID and is recognised from the input VCF, `gvcf_site_locations_from_vcf`).
    A missing file yields the empty set."""
    sites: set[str] = set()
    if not path.exists():
        return sites
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 7:
                continue
            name = _GENERATED_NAME_RE.match(fields[0])
            if name is None:
                continue
            if not is_gvcf_variant_site(perl_split(name.group(1), "/")[1:]):
                continue
            location = canonical_location(fields[1].strip())
            if canonical_contigs is not None and location.split(":", 1)[0] not in canonical_contigs:
                continue
            sites.add(location)
    return sites


_NON_REF_BLOCK = "<NON_REF>"


def _is_non_ref_block(alt: str) -> bool:
    return alt.upper() == _NON_REF_BLOCK


def _names_an_allele(alt: str) -> bool:
    """Whether an ALT names a variant: not the `<NON_REF>` block, `*` or `.`."""
    return not _is_non_ref_block(alt) and alt not in ("*", ".")


def perl_true(value: str | None) -> bool:
    """Whether Perl reads an INFO value as true: defined, not empty and not `0`."""
    return value not in (None, "", "0")


def perl_fields(fields: Sequence[str]) -> list[str]:
    """The fields Perl's list-context `split` keeps: the trailing empty fields dropped and the
    middle ones kept, so `G,,T,,` is `G`, the empty field, `T`."""
    kept = list(fields)
    while kept and kept[-1] == "":
        kept.pop()
    return kept


def perl_split(text: str, sep: str) -> list[str]:
    """Perl's list-context `split` of ``text`` on the separator ``sep``, as the reference
    splits the ALT column (BaseVCF4.pm:489): `<DEL>,` is the one field `<DEL>` and the empty
    string has no field."""
    return perl_fields(text.split(sep))


def _info_written(info: str, key: str) -> str | None:
    """The value of `key=value` in an INFO column, Perl's `defined $info->{key}`: `None` for
    an absent key or a bare flag, the last occurrence winning as in the hash Ensembl VEP's
    parser builds (BaseVCF4.pm:581-583), whose two-element list assignment
    `my ($key,$value) = split('=',$info)` keeps the text up to a second `=`
    (`END2=30425700=5` is `30425700`, `SVTYPE==DEL` the empty value)."""
    value = None
    for item in info.split(";"):
        name, sep, text = item.partition("=")
        if name == key:
            value = text.partition("=")[0] if sep else None
    return value


# Parser/VCF.pm:241, the ALT forms that route a line to the structural path.
_STRUCTURAL_ALT_LIST = re.compile(r"[<\[\]][^*]+[>\]\[]|^\.\w+|\w+\.$")


def is_gvcf_variant_site(alts: Sequence[str], svtype: str | None = None) -> bool:
    """Whether a record is a gVCF variant site vep-rs reads on its sequence path: `<NON_REF>`
    beside alleles that name a variant, no true INFO/SVTYPE, and no remaining allele written
    in a structural form (Parser/VCF.pm:237-245: without a true SVTYPE the remaining ALTs
    alone decide the path, since a structural form is never a run of `ACGT`). A list of
    `<NON_REF>` alone is a reference block; `*` or `.` beside the block names no allele; a
    block beside a symbolic allele (`<DEL>,<NON_REF>`) is an unsupported structural list,
    the skipped-record class; under a true SVTYPE every ALT, the block included, is an
    allele of a structural record in both engines. A trailing empty field of the ALT list is
    no field (`perl_fields`)."""
    alts = perl_fields(alts)
    real = [alt for alt in alts if _names_an_allele(alt)]
    if not real or not any(map(_is_non_ref_block, alts)) or perl_true(svtype):
        return False
    return _STRUCTURAL_ALT_LIST.search(",".join(real)) is None


def _minimised(ref: str, alt: str) -> tuple[str, str, int]:
    """A bi-allelic indel minimised as vep-rs minimises it (`trim_alleles`): the shared
    prefix is trimmed while both alleles have bases, then the shared suffix, and an emptied
    allele reads `-`. Returns the two alleles and the bases trimmed in front, by which the
    start moves."""
    shared = 0
    while shared < min(len(ref), len(alt)) and ref[shared].upper() == alt[shared].upper():
        shared += 1
    ref, alt = ref[shared:], alt[shared:]
    while ref and alt and ref[-1].upper() == alt[-1].upper():
        ref, alt = ref[:-1], alt[:-1]
    return ref or "-", alt or "-", shared


def gvcf_site_rows(pos: int, ref: str, alts: Sequence[str]) -> list[tuple[str, int, int]]:
    """The (allele, start, end) vep-rs writes for each real allele of a gVCF variant site,
    the sequence path of its parser: `<NON_REF>`, `*` and `.` name no allele; when any ALT's
    length differs from REF's, a `*` of one base counted (Parser/VCF.pm:296-298), the first
    base REF shares with every ALT is trimmed and the start moves one base, an emptied allele
    reading `-`, an ALT carrying a `*` (`<*>`) kept whole and, beside other ALTs, left out of
    the comparison (:306-316 with several ALTs, :328 with one; a `.` counts and never shares
    the base); a single remaining allele whose length still differs from REF's is minimised,
    its shared prefix then its shared suffix trimmed and the start moved by the prefix; an
    insertion's end is its start minus one; the allele is written in upper case (`a g,<non_ref>`
    is `G`, as Parser.pm:554 writes it). The minimisation needs the one ALT beside the block to
    name the allele: a `*` or `.` beside the indel keeps the anchor-trimmed alleles, as Ensembl
    VEP's `minimise_alleles` leaves a string of three parts as it is (Parser.pm:978-984). A
    trailing empty field of the ALT list is no field (`perl_fields`)."""
    named = [alt for alt in perl_fields(alts) if not _is_non_ref_block(alt)]
    real = [alt for alt in named if _names_an_allele(alt)]
    if not real:
        return []
    is_indel = any(len(alt) != len(ref) for alt in named)
    # With one ALT its own first base decides (:328): a lone `<*>` filtered out would leave
    # nothing to compare, and `all` of nothing anchors. With several, the alleles carrying a
    # `*` are left out of the comparison (:306).
    if len(named) == 1:
        compared = named
    else:
        compared = [alt for alt in named if "*" not in alt]
    anchored = bool(ref) and all(alt[:1].upper() == ref[:1].upper() for alt in compared)
    prefix = 1 if is_indel and anchored else 0
    trimmed_ref = ref[prefix:] or "-"
    rows: list[tuple[str, int, int]] = []
    for alt in real:
        # An allele carrying a `*` (`<*>`) is kept as written (:315).
        allele = alt if "*" in alt else (alt[prefix:] or "-")
        start = pos + prefix
        ref_part = trimmed_ref
        if len(named) == 1 and "-" not in (ref_part, allele) and len(ref_part) != len(allele):
            ref_part, allele, shared = _minimised(ref_part, allele)
            start += shared
        end = start - 1 if ref_part == "-" else start + len(ref_part) - 1
        rows.append((allele.upper(), start, end))
    return rows


def sequence_location(chrom: str, start: int, end: int) -> str:
    """The Location both engines print for a span: `chrom:start`, or `chrom:lo-hi` with an
    insertion's `end = start - 1` written first (`21:30425650-30425651`)."""
    lo, hi = min(start, end), max(start, end)
    return f"{chrom}:{lo}" if lo == hi else f"{chrom}:{lo}-{hi}"


def gvcf_site_locations_from_vcf(
    path: Path, canonical_contigs: frozenset[str] | None = None
) -> set[str]:
    """The Locations vep-rs writes for every gVCF variant site of an input VCF (plain or
    gzipped), recognised by the record's shape (`is_gvcf_variant_site`) whatever its ID, at the
    POS as vep-rs reads it (`vep_rs_pos`); a record whose POS vep-rs refuses has no rows."""
    sites: set[str] = set()
    opener = gzip.open if path.suffix in (".gz", ".bgz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, ref, alt = fields[0], fields[3], fields[4]
            pos = vep_rs_pos(fields[1])
            alts = perl_split(alt, ",")
            svtype = _info_written(fields[7], "SVTYPE") if len(fields) > 7 else None
            if pos is None or not is_gvcf_variant_site(alts, svtype):
                continue
            if canonical_contigs is not None and chrom not in canonical_contigs:
                continue
            for _allele, start, end in gvcf_site_rows(pos, ref, alts):
                sites.add(canonical_location(sequence_location(chrom, start, end)))
    return sites


def input_vcf_for(output: Path, input_dir: Path | None, input_vcf: Path | None) -> Path | None:
    """The input VCF behind an output file: `--input-vcf` when given, else the file in
    `--input-dir` with the output's stem and a `.vcf`, `.vcf.gz` or `.vcf.bgz` suffix."""
    if input_vcf is not None:
        return input_vcf
    if input_dir is None:
        return None
    for suffix in (".vcf", ".vcf.gz", ".vcf.bgz"):
        candidate = input_dir / f"{output.stem}{suffix}"
        if candidate.exists():
            return candidate
    return None


def no_input_vcf_warning(stem: str, input_dir: Path) -> str:
    """The warning written to stderr for an output whose stem `--input-dir` holds no input VCF
    for: the comparison goes on, recognising a gVCF variant site from an ID-less record's
    generated name alone (`gvcf_site_locations`)."""
    return (
        f"WARNING: [compare_vep_outputs] no input VCF for {stem} under {input_dir}; "
        "gVCF variant sites are recognised from an ID-less record's generated name only"
    )


def _non_integer_end(location: str) -> bool:
    """Whether a coordinate of the Location is not an integer (`21:30425651-30425949.5`,
    `21:30425651-NaN`, `21:30425949.0`); a coordinate the reference numified below zero
    (`21:-4-295`) is an integer."""
    _chrom, sep, coords = location.partition(":")
    return bool(sep) and re.fullmatch(r"-?\d+(?:--?\d+)?", coords, re.ASCII) is None


def apply_116_2_one_sided_rules(
    rows: list[dict[str, str]],
    masked_rows: set[int],
    scope: ReferenceRecordScope,
) -> dict[str, dict[str, int]]:
    """The 116.2 one-sided pass over the rows no pair rule removed.

    Marks the set-aside rows in ``masked_rows``: vep-rs's rows of
    `REFERENCE_SKIPPED_RECORD_BUCKET` and `GVCF_VARIANT_SITE_DROPPED_BUCKET`, the
    reference's rows of `NON_INTEGER_END_BUCKET` and of `NON_NUMERIC_POS_NUMIFIED_BUCKET` (a
    negative first coordinate at a Location vep-rs wrote nothing for, the one spelling of a
    numified POS no written POS reaches without the input; a POS numified to 0 prints the
    start 1, as POS 0 does in both engines, and
    `30_425_949` prints 31, so those are the structural-variant comparator's), and both
    rows of a `DOT_BREAKEND_PAIR_BUCKET` pair; counts, without masking them, the two
    lost-row breakend shapes, the FASTA-named slice shape and the reference's rows at a
    record whose coordinate overflows the buffer (`BND_MATE_32_BIT_BUFFER_BUCKET`).
    `STAR_ALONE_SVTYPE_ANNOTATED_BUCKET` and `CNV_TR_RN_ZERO_ALLELE_DROPPED_BUCKET` need
    the input record and count nothing here.
    Each bucket's ``perl`` and ``rust`` say how many of its excluded rows are the
    reference's and vep-rs's. The reference-allele rows of a bare breakpoint
    (`BND_BARE_REFERENCE_ALLELE_ROW_BUCKET`) never reach ``rows``, the unique-tuple key
    stream having dropped them as duplicates, so the class reports the rows the scope
    counted on the reference's raw output, all of them excluded, and nothing is masked
    for it. Returns the per-bucket statistics: ``count`` rows recognised, ``excluded``
    rows set aside, ``records`` distinct (file, Location) pairs the bucket reached.
    """
    stats = {
        rule.bucket: {"count": 0, "excluded": 0, "records": 0, "perl": 0, "rust": 0}
        for rule in REFERENCE_116_2_ONE_SIDED_RULES
    }
    records: dict[str, set[tuple[str, str]]] = defaultdict(set)
    # The bare-breakpoint rows are outside the tuple sets, so they count under neither
    # side's subtraction.
    bare_rows = sum(scope.bare_reference_rows.values())
    stats[BND_BARE_REFERENCE_ALLELE_ROW_BUCKET].update({"count": bare_rows, "excluded": bare_rows})
    records[BND_BARE_REFERENCE_ALLELE_ROW_BUCKET].update(scope.bare_reference_rows)

    def hit(bucket: str, idx: int, row: dict[str, str], exclude: bool) -> None:
        stats[bucket]["count"] += 1
        records[bucket].add((row["file_name"], row["location"]))
        if exclude:
            stats[bucket]["excluded"] += 1
            stats[bucket]["perl" if row["source"] == "missing_in_rust" else "rust"] += 1
            masked_rows.add(idx)

    # The dot-breakend pair: the reference's span row beside vep-rs's point at the span's start
    # or the base after it (`dot_breakend_points`) under one breakend allele
    # (`chromosome_breakpoint` for a bare `N.`, the `N.` and bracket strings of one whose INFO
    # names a mate) on the same feature; the reference's rows whose Location carries a
    # coordinate that is not an integer; and its rows at a Location vep-rs wrote nothing for
    # whose start the reference numified below zero.
    point_rows: dict[tuple[str, str, str, str], list[int]] = defaultdict(list)
    for idx, row in enumerate(rows):
        if idx in masked_rows or row["source"] != "extra_in_rust" or not is_breakend_allele(row["allele"]):
            continue
        span = _location_span(row["location"])
        if span is not None and span[1] == span[2]:
            point_rows[(row["file_name"], row["location"], row["allele"], row["feature"])].append(idx)
    for idx, row in enumerate(rows):
        if idx in masked_rows or row["source"] != "missing_in_rust":
            continue
        if _non_integer_end(row["location"]):
            hit(NON_INTEGER_END_BUCKET, idx, row, exclude=True)
            continue
        span = _location_span(row["location"])
        if (
            span is not None
            and (row["file_name"], row["location"]) in scope.reference_only
            and span[1] < 0
        ):
            hit(NON_NUMERIC_POS_NUMIFIED_BUCKET, idx, row, exclude=True)
            continue
        if not is_breakend_allele(row["allele"]):
            continue
        if span is None or span[1] >= span[2]:
            continue
        for point in dot_breakend_points(span):
            partners = point_rows.pop((row["file_name"], point, row["allele"], row["feature"]), ())
            if partners:
                hit(DOT_BREAKEND_PAIR_BUCKET, idx, row, exclude=True)
                for rust_idx in partners:
                    hit(DOT_BREAKEND_PAIR_BUCKET, rust_idx, rows[rust_idx], exclude=True)
                break

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
        if idx in masked_rows:
            continue
        group = (row["file_name"], row["location"])
        if row["source"] == "missing_in_rust":
            if group in scope.mate_overflow:
                hit(BND_MATE_32_BIT_BUFFER_BUCKET, idx, row, exclude=False)
            continue
        if row["source"] != "extra_in_rust":
            continue
        allele = row["allele"]
        if group in scope.gvcf_absent:
            hit(GVCF_VARIANT_SITE_DROPPED_BUCKET, idx, row, exclude=True)
            continue
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
        mate = reference_mate_chromosome(allele)
        # A synonym is any mate name, read as the reference reads it, that is not a name the
        # cache files a chromosome under: a prefix the mate-allele name test does not
        # resolve, or an accession the overlap test resolves and the slice test refuses.
        if (
            mate is not None
            and perl_alleles
            and allele not in perl_alleles
            and mate not in CACHE_CONTIG_NAMES
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
    outputs: tuple[Path, Path] | None = None,
    input_records: Callable[[str], Path | None] | None = None,
) -> tuple[int, int, dict[str, dict[str, int]]]:
    """Read discordant.tsv, classify intended-divergence swap pairs, and exclude
    the ones an Ensembl VEP self-contradiction accounts for.

    ``outputs`` names the reference's and vep-rs's output directory. With it, a pair no
    rule of ``EXCLUSION_REGISTRY`` classifies or excludes is put to the release's
    row-aware rules (``row_pair_rules``) in their order, the first match taking it, on the
    rows of the two outputs that carry its tuples (`fetch_pair_rows`); their buckets
    follow the pair buckets. Without it, or under a release that registers no such rule,
    no output is read. ``input_records`` names the input VCF behind an output file of that
    name (``--input-vcf`` or ``--input-dir``), or ``None`` for a file without one: the rows of
    such a pair whose first column spells no alleles take their record's from it
    (`read_record_references`).

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

    The registry carries no rule for a 3'UTR stop_lost trade-off (vep-rs
    ``stop_lost`` where Perl emits a 3'UTR term alone) and no start rule looser
    than ``_is_start_cooccurrence_swap``. The first shape matches no tuple,
    because both engines read the UTR sequence from the cache
    (``vefc.{three,five}_prime_utr``, with FASTA as the fallback). The second
    would set aside pairs where vep-rs misses or drops a start term, which are
    vep-rs divergences.
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
    row_rules = REFERENCE_REGISTRIES[reference_release].row_pair_rules
    untaken: list[tuple[int, int]] = []

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
                else:
                    untaken.append((pi, ri))
                if excluded:
                    excluded_perl += 1
                    excluded_rust += 1
                    masked_rows.update((pi, ri))

    # The row-aware rules, on the pairs no rule above took, in the order of the rows in
    # `discordant.tsv`; a row one of them set aside is not paired again.
    row_stats: dict[str, dict[str, int]] = {
        rule.bucket: {"count": 0, "excluded": 0} for rule in row_rules
    }
    if row_rules and untaken and outputs is not None:
        untaken.sort()
        fetched = fetch_pair_rows(rows, untaken, *outputs)
        if input_records is not None:
            by_file: dict[str, list[PairRow]] = defaultdict(list)
            for index, pair_rows_of_index in fetched.items():
                by_file[rows[index]["file_name"]].extend(pair_rows_of_index)
            for file_name, rows_of_file in by_file.items():
                records = input_records(file_name)
                if records is not None:
                    read_record_references(records, rows_of_file)
        for pi, ri in untaken:
            if pi in masked_rows or ri in masked_rows or pi not in fetched or ri not in fetched:
                continue
            pair_rows = (fetched[pi], fetched[ri])
            rule = next(
                (
                    r
                    for r in row_rules
                    if r.matches(rows[pi]["consequence_set"], rows[ri]["consequence_set"], rows[pi]["allele"], pair_rows)
                ),
                None,
            )
            if rule is None:
                continue
            row_stats[rule.bucket]["count"] += 1
            if rule.excludes:
                row_stats[rule.bucket]["excluded"] += 1
                excluded_perl += 1
                excluded_rust += 1
                masked_rows.update((pi, ri))
    bucket_stats.update(row_stats)

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
        excluded_perl += sum(s["perl"] for s in one_sided_stats.values())
        excluded_rust += sum(s["rust"] for s in one_sided_stats.values())
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
        f"reported floor (default {PER_CLASS_MIN_TUPLES}, fixed independently of the "
        f"data). Terms below it are still written to the CSV with their counts.",
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
        "--input-dir",
        help="Directory of the input VCFs (.vcf, .vcf.gz) the outputs were produced from, "
        "paired with each output by stem. A gVCF variant site is then recognised by its "
        "shape (<NON_REF> beside a real allele) whatever its ID; without an input only an "
        "ID-less record's generated name reveals the site.",
    )
    parser.add_argument(
        "--input-vcf",
        help="The one input VCF every output pair was produced from (a single-input run); "
        "takes precedence over --input-dir.",
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
        default=None,
        help="Reference assembly. GRCh37 and GRCh38 select the human canonical contig "
        "set when --canonical-contigs (default) names none; any other assembly's set "
        "comes from the --canonical-contigs value (run_concordance.sh passes the "
        "cache's info.json contig_order).",
    )
    parser.add_argument(
        "--canonical-contigs",
        nargs="?",
        const=True,
        default=True,
        metavar="FILE_OR_LIST",
        help="Defensive contig-allowlist filter (DEFAULT ON). Drops 5-tuples whose "
        "Location contig is not in the canonical set. Bare, the set is the human one "
        "for --assembly GRCh37 or GRCh38; with a value, the set is the comma-separated "
        "names given, or the names in the file named, one per line or comma-separated, "
        "`#` lines ignored (a non-human cache's info.json contig_order). Catches "
        "alt-contig tuples that slip past the input filter (engine cache annotations, etc.).",
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
        "of 116.1's chromosome-aware overlap are counted, as is the reference-allele row "
        "of a bare breakpoint. The pair rules are the same "
        "under both. run_concordance.sh and run_clone_measurement.sh derive the value "
        "from the reference set's provenance.",
    )
    args = parser.parse_args(argv)
    if args.input_vcf is not None and not Path(args.input_vcf).is_file():
        parser.error(f"--input-vcf {args.input_vcf} is not a file")
    if args.input_dir is not None and not Path(args.input_dir).is_dir():
        parser.error(f"--input-dir {args.input_dir} is not a directory")
    return args


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


# The human canonical contig set, Ensembl-style, the same for GRCh37 and GRCh38.
# Must match prepare_benchmark_vcfs.py CANONICAL_CONTIGS frozenset and
# scripts/concordance/run_concordance.sh CANONICAL_CONTIGS_LIST. M and MT are
# equivalent mitochondrial names.
#
# Rationale: input-prep auto-strips the UCSC chr- prefix so all engine
# outputs emit Ensembl-style Locations. This defensive output-side filter is
# therefore an Ensembl-only membership test on engine output Locations.
CANONICAL_CONTIGS = frozenset([str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"])
# The names the Ensembl caches file chromosomes under, for the breakend synonym test: a
# mate written `M` (or `chrM`) is a synonym of `MT` there, as any other name is.
CACHE_CONTIG_NAMES = CANONICAL_CONTIGS - {"M"}


def canonical_set_for_assembly(assembly: str | None) -> frozenset[str] | None:
    """Return the human canonical contig set if the assembly is a human one.

    The same Ensembl-style set is used for both GRCh37 and GRCh38 (the assembly
    difference is in the genome-sequence content, not contig naming). Any other
    value, including None, returns None: another species' set is not derivable
    from its assembly name and comes from `--canonical-contigs <file or list>`.
    """
    if assembly in ("GRCh37", "GRCh38"):
        return CANONICAL_CONTIGS
    return None


def parse_canonical_contigs(value: str) -> frozenset[str]:
    """The contig set a `--canonical-contigs` value names: the lines of the file
    at `value` when one exists (comma- or whitespace-separated names, `#` lines
    ignored), else the comma-separated names in `value` itself. An empty set is
    a usage error: a filter that admits nothing would score every tuple absent."""
    path = Path(value)
    if path.is_file():
        text = " ".join(
            line
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        )
    else:
        text = value
    names = frozenset(name for name in re.split(r"[,\s]+", text) if name)
    if not names:
        raise SystemExit(f"error: --canonical-contigs {value!r} names no contig")
    return names


def normalize_consequence_set(raw: str) -> str:
    terms = [term.strip() for term in raw.split(",") if term.strip()]
    return ",".join(sorted(set(terms)))


def collapse_transcript_context_terms(raw: str) -> str:
    terms = [term.strip() for term in raw.split(",") if term.strip()]
    kept = sorted({term for term in terms if term not in TRANSCRIPT_CONTEXT_TERMS})
    if kept:
        return ",".join(kept)
    return "__context_only__"


def canonical_location(location: str) -> str:
    """The Location with each coordinate written as its number: Ensembl VEP prints an
    INFO/END as given, so `21:30425651-030425949` and `21:30425651-+30425949` name the
    same span as `21:30425651-30425949`. A coordinate that is not an integer (`30425949.5`)
    stays as written."""
    chrom, sep, coords = location.partition(":")
    if not sep:
        return location
    return chrom + ":" + "-".join(_canonical_coordinate(c) for c in coords.split("-"))


def _canonical_coordinate(text: str) -> str:
    digits = text.removeprefix("+")
    return str(int(digits)) if _ten_digits(digits) else text


def iter_vep_key_lines(
    path: Path,
    consequence_transform: Callable[[str], str] | None = None,
    canonical_contigs: frozenset[str] | None = None,
) -> Iterator[str]:
    """Yield normalized tuple keys as tab-separated strings.

    Key format (5 columns):
      location, allele, feature, feature_type, consequence_set
    The Location carries each coordinate as its number (`canonical_location`).

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

            location = canonical_location(fields[1].strip())
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
    floor, fixed independently of the data: a minimum over per-(term, dataset) cells
    is a worst-case outlier from a tiny denominator rather than a property of the
    engine.
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


class KeyFileOrderError(ValueError):
    """A sorted key file holds a key that does not sort after the key before it."""


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

    Each file must ascend in the order the merge compares keys by, with every key once,
    as `_sort_unique_stream_to_file` writes it. At the first key that does not sort after
    the key before it the merge raises `KeyFileOrderError` naming the file and the line: a
    merge that went on would count a key both files hold as missing and as extra at once.

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

    def next_ascending(handle: TextIO, path: Path, previous: str, line: int) -> str | None:
        """Read the key on line `line` of `path`, which must sort after `previous`."""
        key = _next_line(handle)
        if key is not None and key <= previous:
            raise KeyFileOrderError(
                f"{path}: line {line}: the key does not sort after the key before it; "
                "the comparison needs each key file sorted with every key once"
            )
        return key

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
                p = next_ascending(perl_h, perl_sorted, p, perl_count + 1)
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
                r = next_ascending(rust_h, rust_sorted, r, rust_count + 1)
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
            p = next_ascending(perl_h, perl_sorted, p, perl_count + 1)
            r = next_ascending(rust_h, rust_sorted, r, rust_count + 1)

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


@contextmanager
def _key_file_directory(
    report_dir: Path, discordant_out: TextIO, discordant_path: Path
) -> Iterator[str]:
    """Yield a temporary directory under `report_dir` for one comparison's key files.

    The directory is removed when the block ends, by every ending but one. Where the
    comparison stops on a key file out of order (`KeyFileOrderError`), the directory stays
    so the files the error names can be read, and the discordance report written so far is
    closed and renamed with a `.partial` suffix so nothing reads it as a result. The error
    raised in its place names both.
    """
    tmp = tempfile.mkdtemp(prefix="compare_vep_tmp_", dir=str(report_dir))
    keep = False
    try:
        yield tmp
    except KeyFileOrderError as error:
        keep = True
        discordant_out.close()
        partial = discordant_path.with_name(f"{discordant_path.name}.partial")
        discordant_path.replace(partial)
        raise KeyFileOrderError(
            f"{error}; the key files are kept in {tmp} and the rows written before the "
            f"stop in {partial}"
        ) from None
    finally:
        if not keep:
            shutil.rmtree(tmp)


def resolve_canonical_contigs(args: argparse.Namespace) -> frozenset[str] | None:
    """Resolve the defensive contig filter, warning when it is off by omission.

    Factored so the compare and per-class paths cannot drift on it. They must agree:
    per-class VEP totals counted WITHOUT the filter against a discordance report
    produced WITH it would inflate every term's VEP total by the dropped alt-contig
    tuples, which reads as depressed recall on every class at once.

    `--canonical-contigs` is False (off), True (bare: the human set when `--assembly`
    is GRCh37 or GRCh38) or a value naming the set (`parse_canonical_contigs`).
    """
    if args.canonical_contigs is False:
        return None
    if args.canonical_contigs is not True:
        return parse_canonical_contigs(args.canonical_contigs)
    canonical_contigs = canonical_set_for_assembly(args.assembly)
    if canonical_contigs is None:
        # Soft-warn rather than hard-fail: the upstream input filter already
        # enforced the canonical-contigs constraint; this script is a defensive
        # complement, and a missing --assembly means the defensive layer is off
        # while the upstream filter applies. Print the warning so the
        # operator sees the asymmetry in the report.
        print(
            "WARNING: --canonical-contigs (default ON) requested but --assembly is "
            f"{args.assembly or 'not provided'}, which names no canonical set. "
            "Defensive output-side filter is OFF for this run. The upstream "
            "input-prep canonical-contigs filter (prepare_benchmark_vcfs.py) still "
            "applies. Pass --assembly GRCh37 or --assembly GRCh38, or name the set "
            "with --canonical-contigs <file or list> (a cache's info.json "
            "contig_order), to enable the defensive layer, or --no-canonical-contigs "
            "to silence this warning.",
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
                f"floor can be reported under the fixed threshold"
            )
            continue
        print(
            f"{engine}: floor {floor.f1} on {floor.term}; "
            f"{floor.n_eligible} of {len(rows)} terms at or above {args.min_tuples} "
            f"pooled VEP tuples, {floor.n_below_threshold} below"
        )
    return 0


def main() -> int:
    # A coordinate is read whole at any length, where CPython by default refuses to convert
    # an integer's text of more than 4,300 digits.
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
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
    input_dir = Path(args.input_dir).resolve() if args.input_dir else None
    input_vcf_path = Path(args.input_vcf).resolve() if args.input_vcf else None

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

        with _key_file_directory(report_dir, discordant_out, discordant_path) as tmp:
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
            if record_scope is not None:
                record_scope.record_bare_reference_rows(
                    perl_file.name,
                    count_bare_breakpoint_reference_rows(perl_file, canonical_contigs),
                )
                gvcf_sites = gvcf_site_locations(rust_file, canonical_contigs)
                input_vcf = input_vcf_for(perl_file, input_dir, input_vcf_path)
                if input_vcf is not None:
                    gvcf_sites |= gvcf_site_locations_from_vcf(input_vcf, canonical_contigs)
                elif input_dir is not None:
                    print(no_input_vcf_warning(perl_file.stem, input_dir), file=sys.stderr)
                record_scope.record_gvcf_sites(perl_file.name, gvcf_sites)

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
            if record_scope is not None:
                record_scope.record_reference_named(
                    perl_file.name,
                    reference_named_locations(
                        perl_file,
                        rust_file,
                        record_scope.absent_locations(perl_file.name),
                        canonical_contigs,
                    ),
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
        outputs=(perl_dir, rust_dir),
        input_records=lambda name: input_vcf_for(perl_dir / name, input_dir, input_vcf_path),
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
        aggregate["adjusted"]["divergence_classes"] = report_class_definitions(
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
        "canonical_contigs_filter_active": canonical_contigs is not None,
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
