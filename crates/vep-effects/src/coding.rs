// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Coding-region effect analysis: codon extraction, translation, and amino acid comparison.
//!
//! Given a variant position within the CDS, this module extracts the affected codon(s), translates
//! them under the transcript's genetic code (SNV path and `PerlCodingEval`) or the standard one,
//! and compares reference vs alternate amino acids. Handles SNVs, inframe indels, and frameshifts,
//! mirroring the Ensembl VEP `TranscriptVariationAllele` codon/peptide logic for concordance.
//!
//! Ensembl VEP source citations name modules of ensembl-variation release/116
//! (`Bio/EnsEMBL/Variation/...`); `TranscriptMapper`, `Mapper`, `Transcript` and `Slice` are in
//! ensembl core release/116, `AnnotationType/Transcript.pm` and `Parser/...` in ensembl-vep
//! release/116.3, and `Bio::Tools::CodonTable` and `Bio::SeqUtils` in BioPerl.
//!
//! Key functions:
//! - [`get_codon_change`]: codon-level effect for a single variant
//! - `classify_inframe_indel_by_peptide`: inframe insertion/deletion/protein-altering classification
//! - `compute_peptide_alleles`: full CDS translation with trimmed peptide comparison
//! - `compute_codon_window_peptide_alleles`: Ensembl VEP-style local codon-window translation

use vep_core::codon::{
    is_start_codon, is_ter_codon, seq_is_dna, translate_codon, translate_codon_with_table,
};
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

use crate::mapper::CdsSpanBounds;

/// A codon change caused by a variant.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct CodonChange {
    /// Reference codon with the variant base in uppercase (e.g., "gCc").
    pub ref_codon: String,
    /// Alternate codon with the variant base in uppercase (e.g., "gTc").
    pub alt_codon: String,
    /// Reference amino acid single-letter code (e.g., 'A').
    pub ref_amino_acid: u8,
    /// Alternate amino acid single-letter code (e.g., 'V').
    pub alt_amino_acid: u8,
    /// 1-based codon number in the CDS.
    pub codon_number: u64,
    /// Whether this change causes a frameshift.
    pub is_frameshift: bool,
}

/// Compute frameshift status using CDS-mapped variant lengths (Ensembl VEP parity).
///
/// Ensembl VEP's `VariationEffect::frameshift` predicate uses:
///   `var_len = cds_end - cds_start + 1`  (not raw ref allele length)
///   `allele_len = seq_length`             (alt allele length)
///   `return abs(allele_len - var_len) % 3`
///
/// For insertions with CDS-inverted coords (`cds_start > cds_end`), `var_len = 0`.
/// For variants partially spanning introns/UTRs, `var_len` may be smaller than
/// the raw reference allele length, changing the mod-3 result.
pub(crate) fn is_frameshift_cds_aware(variant: &InputVariant, bounds: &CdsSpanBounds) -> bool {
    let vf_nt_len = ((bounds.cds_end as isize) - (bounds.cds_start as isize) + 1).max(0) as usize;
    let alt = variant.alt_allele();
    let allele_len = if alt == b"-" || alt.is_empty() {
        0
    } else {
        alt.len()
    };
    allele_len.abs_diff(vf_nt_len) % 3 != 0
}

/// Sub-classification of inframe indels based on peptide-level analysis.
///
/// Distinguishes between `inframe_insertion`, `inframe_deletion`, and
/// `protein_altering_variant` consequence terms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum InframeIndelKind {
    Insertion,
    Deletion,
    ProteinAltering,
    LeadingStopGained,
}

/// Ensembl VEP's `protein_altering_variant` predicate suppresses protein_altering when
/// the alt peptide contains a gained stop and has prefix/suffix containment with
/// the ref peptide (`$alt_pep =~ /^\Q$ref_pep\E|\Q$ref_pep\E$/`).
fn has_gained_stop_with_containment(ref_pep: &[u8], alt_pep: &[u8]) -> bool {
    !ref_pep.is_empty()
        && alt_pep.contains(&b'*')
        && !ref_pep.contains(&b'*')
        && (alt_pep.starts_with(ref_pep) || alt_pep.ends_with(ref_pep))
}

/// Classify an inframe indel using codon-level heuristics, mirroring
/// Ensembl VEP's `inframe_deletion` / `inframe_insertion` / `protein_altering_variant`
/// predicate logic.
///
/// Ensembl VEP's `inframe_deletion` checks at the **codon string** level:
/// - alt_codon is a prefix or suffix of ref_codon -> inframe_deletion
/// - after trimming shared prefix/suffix, if alt remainder is empty -> inframe_deletion
/// - otherwise -> protein_altering_variant
///
/// Returns `None` if the codon alleles cannot be derived (caller uses fallback).
pub(crate) fn classify_inframe_indel_by_peptide(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
    cds_bounds: Option<&CdsSpanBounds>,
) -> Option<InframeIndelKind> {
    let vefc = transcript.vefc.as_ref()?;
    let translateable_seq = vefc.translateable_seq.as_ref()?;
    let cds = translateable_seq.as_bytes();

    let normalized = normalize_indel_edit_in_transcript_sense(variant, transcript, cds_pos)?;
    let idx = normalized.cds_idx;
    if idx >= cds.len() {
        return None;
    }

    let ref_is_dash = normalized.ref_seq.is_empty();
    let alt_is_dash = normalized.alt_seq.is_empty();

    // A pure deletion, including an allele whose shared flanks reduce it to one
    // (`CGGA/C` of a multi-allelic record, which VEP never minimises), is an
    // in-frame deletion here: the caller has ruled out a frameshift, and VEP's
    // `inframe_deletion` matches an alt codon window that is the ref window minus
    // the deleted bases. The net-deletion branch below is not used for it, because
    // its leading-stop early return would suppress `inframe_deletion` where VEP
    // emits it beside `stop_gained`; the supplemental stop check in
    // consequences.rs adds the stop terms.
    if alt_is_dash {
        return Some(InframeIndelKind::Deletion);
    }

    let alt_seq = normalized.alt_seq;

    // For pure insertions (ref="-"), the insertion point is at cds_pos.
    // The ref codon is the codon containing that position; the alt codon
    // is the ref codon with inserted bases spliced in.
    let ref_len = normalized.ref_seq.len();

    if !ref_is_dash && idx + ref_len > cds.len() {
        return None;
    }

    // Ensembl VEP's _get_codon_alleles returns the codons covering the variant span,
    // aligned to codon boundaries.
    let codon_start = (idx / 3) * 3;
    let end = idx + ref_len.max(1); // for insertions, span at least 1 base for codon calc
    let end_codon = end.div_ceil(3) * 3; // ceil to codon boundary
    let ref_codon_end = end_codon.min(cds.len());

    let ref_codon = &cds[codon_start..ref_codon_end];
    let rel = idx - codon_start;
    let mut alt_codon = Vec::with_capacity(ref_codon.len().saturating_sub(ref_len) + alt_seq.len());
    alt_codon.extend_from_slice(&ref_codon[..rel]);
    alt_codon.extend_from_slice(&alt_seq);
    if rel + ref_len < ref_codon.len() {
        alt_codon.extend_from_slice(&ref_codon[rel + ref_len..]);
    }

    // Uppercase in-place for comparison (avoids allocating new Vecs).
    let mut ref_upper: Vec<u8> = ref_codon.to_vec();
    ref_upper
        .iter_mut()
        .for_each(|b| *b = b.to_ascii_uppercase());
    alt_codon
        .iter_mut()
        .for_each(|b| *b = b.to_ascii_uppercase());
    let alt_upper = alt_codon;

    // Ensembl VEP's `VariationEffect::inframe_insertion` returns 0 when `start_lost`
    // fires, but `inframe_deletion` has no start_lost gate, so the Met-loss
    // check is deferred to the insertion branch and the deletion branch
    // translates its own peptides without gating.
    let normalized_cds_pos = normalized.cds_idx as u64 + 1;
    let start_codon_peps = if normalized_cds_pos <= 3 {
        let ref_pep = translate_to_peptide(&ref_upper);
        let alt_pep = translate_to_peptide(&alt_upper);
        Some((ref_pep, alt_pep))
    } else {
        None
    };

    if alt_upper.len() > ref_upper.len() {
        // Ensembl VEP's inframe_insertion uses peptide comparison (not codon):
        // alt_pep starts/ends with ref_pep -> inframe_insertion.
        let (manual_ref_pep, manual_alt_pep) = match start_codon_peps {
            Some((r, a)) => (r, a),
            None => (
                translate_to_peptide(&ref_upper),
                translate_to_peptide(&alt_upper),
            ),
        };

        // Ensembl VEP's inframe_insertion returns 0 if start_lost(@_); an insertion
        // destroying the Met falls through to protein_altering_variant, itself
        // suppressed by its own `return 0 if start_lost(@_)`. consequences.rs
        // emits the start_lost label, so this routes to ProteinAltering only so
        // that inframe_insertion does not fire.
        if normalized_cds_pos <= 3
            && manual_ref_pep.first() == Some(&b'M')
            && manual_alt_pep.first() != Some(&b'M')
        {
            return Some(InframeIndelKind::ProteinAltering);
        }

        // Gate: Ensembl VEP's inframe_insertion returns 0 if
        // start_retained_variant(@_) && alt_pep ends with ref_pep. This prevents
        // inframe_insertion when the insertion is before/at the start codon and the
        // Met is preserved. Must check before stop truncation (Ensembl VEP order).
        if normalized_cds_pos <= 3
            && manual_ref_pep.first() == Some(&b'M')
            && manual_alt_pep.first() == Some(&b'M')
            && manual_alt_pep.ends_with(&manual_ref_pep)
        {
            return Some(InframeIndelKind::ProteinAltering);
        }

        // The containment check prefers compute_codon_window_peptide_alleles, which
        // builds the full alt CDS like Ensembl VEP's _get_alternate_cds: the manual codon
        // window is too narrow and over-calls protein_altering. Start codon
        // variants (cds_pos <= 3) keep the manual window and its gates.
        //
        // For a pure insertion spanning a codon boundary (translation_start !=
        // translation_end) Ensembl VEP's inverted translation bounds yield a zero-length
        // ref codon window, so containment passes trivially; un-inverted bounds
        // give a wider window that breaks it, hence the empty-ref case.
        let cross_codon_insertion =
            ref_is_dash && cds_bounds.is_some_and(|b| b.translation_start != b.translation_end);
        let (ref_pep, mut alt_pep) = if normalized_cds_pos > 3
            && cross_codon_insertion
            && !normalized.trimmed_transcript_flanks
        {
            // Pure insertion spanning codon boundary: Ensembl VEP has empty ref_pep.
            (Vec::new(), translate_to_peptide(&alt_seq))
        } else if normalized_cds_pos > 3 && !normalized.trimmed_transcript_flanks {
            match cds_bounds
                .and_then(|b| compute_codon_window_peptide_alleles(variant, transcript, b))
            {
                Some(peps) => peps,
                None => (manual_ref_pep, manual_alt_pep),
            }
        } else {
            (manual_ref_pep, manual_alt_pep)
        };

        // Whether the alt peptide starts with a stop decides ProteinAltering versus
        // LeadingStopGained below.
        let alt_starts_with_stop = alt_pep.first() == Some(&b'*') && ref_pep.first() != Some(&b'*');

        // Save the original (pre-truncation) alt_pep. Ensembl VEP's protein_altering_variant
        // predicate checks containment on untrimmed peptides (`$alt_pep =~ /^\Q$ref_pep\E
        // |\Q$ref_pep\E$/`), while inframe_insertion truncates at the first stop before
        // its containment check. When untrimmed containment passes but truncated fails,
        // Ensembl VEP suppresses protein_altering and only stop_gained fires.
        let original_alt_pep = alt_pep.clone();

        // Ensembl VEP: trim everything after first stop in alt_pep for insertion check.
        if let Some(stop_pos) = alt_pep.iter().position(|&b| b == b'*') {
            alt_pep.truncate(stop_pos + 1);
        }

        if !alt_pep.is_empty()
            && (ref_pep.is_empty()
                || alt_pep.starts_with(&ref_pep[..])
                || alt_pep.ends_with(&ref_pep[..]))
        {
            return Some(InframeIndelKind::Insertion);
        }
        // For cross-codon insertions with empty ref_pep, containment is trivially true.
        if ref_pep.is_empty() && !alt_pep.is_empty() {
            return Some(InframeIndelKind::Insertion);
        }

        // Ensembl VEP's protein_altering_variant predicate returns 0 if alt_pep starts
        // with '*'. This suppresses protein_altering and lets stop_gained stand
        // alone. Only applies when containment failed (would have been
        // ProteinAltering); when containment passes, Ensembl VEP's inframe_insertion
        // fires alongside stop_gained.
        if alt_starts_with_stop {
            return Some(InframeIndelKind::LeadingStopGained);
        }

        // When untrimmed containment passes but truncated containment fails
        // (because stop truncation removed the matching suffix), and the alt
        // peptide contains a gained stop, Ensembl VEP suppresses protein_altering.
        if has_gained_stop_with_containment(&ref_pep, &original_alt_pep) {
            return Some(InframeIndelKind::LeadingStopGained);
        }

        return Some(InframeIndelKind::ProteinAltering);
    }

    if alt_upper.len() < ref_upper.len() {
        // Net deletion at codon level.
        // Ensembl VEP's protein_altering_variant predicate returns 0 when alt_pep starts
        // with '*'. For deletions that create a stop at the first alt codon, Ensembl VEP
        // emits only stop_gained without protein_altering or inframe_deletion.
        // Check this before the codon containment test.
        let del_ref_pep = translate_to_peptide(&ref_upper);
        let del_alt_pep = translate_to_peptide(&alt_upper);
        if del_alt_pep.first() == Some(&b'*') && del_ref_pep.first() != Some(&b'*') {
            return Some(InframeIndelKind::LeadingStopGained);
        }

        // Ensembl VEP's inframe_deletion uses codon-level comparison:
        // alt_codon is prefix/suffix of ref_codon -> inframe_deletion.
        if ref_upper.starts_with(&alt_upper) || ref_upper.ends_with(&alt_upper) {
            return Some(InframeIndelKind::Deletion);
        }
        let (_r_trimmed, a_trimmed) = trim_common_flanks(&ref_upper, &alt_upper);
        if a_trimmed.is_empty() {
            return Some(InframeIndelKind::Deletion);
        }

        // When codon-level containment fails but peptide-level containment
        // passes, and the alt peptide contains a gained stop, Ensembl VEP suppresses
        // protein_altering and only stop_gained fires.
        if has_gained_stop_with_containment(&del_ref_pep, &del_alt_pep) {
            return Some(InframeIndelKind::LeadingStopGained);
        }

        return Some(InframeIndelKind::ProteinAltering);
    }

    // Same-length codons with different allele lengths at the DNA level.
    Some(InframeIndelKind::ProteinAltering)
}

pub(crate) fn trim_common_flanks<'a>(a: &'a [u8], b: &'a [u8]) -> (&'a [u8], &'a [u8]) {
    let min_len = a.len().min(b.len());
    let mut prefix = 0usize;
    while prefix < min_len && a[prefix] == b[prefix] {
        prefix += 1;
    }
    let mut suffix = 0usize;
    while suffix < (min_len - prefix) && a[a.len() - 1 - suffix] == b[b.len() - 1 - suffix] {
        suffix += 1;
    }
    let a_end = a.len().saturating_sub(suffix);
    let b_end = b.len().saturating_sub(suffix);
    (&a[prefix..a_end], &b[prefix..b_end])
}

pub(crate) fn is_unambiguous_dna(seq: &[u8]) -> bool {
    seq.iter()
        .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

/// Whether the reference protein holds an `X` among the residues the
/// notation of an inserted peptide of `reach` residues is held against: the
/// `reach` residues before position `start` and after position `end` (both
/// 1-based), and one more on each side. The notation writes the first stop of
/// the inserted residues as `X` (TranscriptVariationAllele.pm:2060-2061) and
/// holds the inserted residues, those past the stop among them, against the
/// residues after the span (:2573-2595) and before it (:2442-2475). That `X`
/// matches an `X` of the reference protein alone, which a coding sequence
/// that opens or ends inside a codon holds; beside one, the notation of an
/// allele of unknown bases would rest on the unknown bases.
pub(crate) fn unknown_residue_within_reach(
    reference: &[u8],
    start: i64,
    end: i64,
    reach: usize,
) -> bool {
    let reach = reach as i64 + 1;
    let from = (start - 1 - reach).clamp(0, reference.len() as i64) as usize;
    let to = (end + reach).clamp(0, reference.len() as i64) as usize;
    from < to && reference[from..to].contains(&b'X')
}

/// Whether a reference allele of at least one base and an alternate allele
/// open with the same base or close with the same base: bases the engines
/// cut from a record of one alternate allele and leave on a record of more.
pub(crate) fn alleles_share_an_end_base(reference: &[u8], alternate: &[u8]) -> bool {
    !reference.is_empty()
        && reference != b"-"
        && (reference.first() == alternate.first() || reference.last() == alternate.last())
}

/// Whether the cut of the bases two alleles share at their ends reads known
/// bases of `alternate` alone: each base it removes and the base it stops
/// at, at each end, unless it runs `reference` out first. Both engines cut
/// the shared bases of a record with one alternate allele, from the left
/// and then from the right (Parser.pm:860-881; `trim_sequences`,
/// Utils/Sequence.pm:1006, 1019), and annotate the span that is left. Where
/// the cut meets an unknown base of `alternate` beside a base of
/// `reference`, one completion of it is cut and another is not, and the
/// two are read at different spans. An allele written `-` is empty.
pub(crate) fn shared_ends_cut_at_known_bases(reference: &[u8], alternate: &[u8]) -> bool {
    fn bases(allele: &[u8]) -> &[u8] {
        if allele == b"-" {
            &[]
        } else {
            allele
        }
    }
    let (mut reference, mut alternate) = (bases(reference), bases(alternate));
    while let (Some(r), Some(a)) = (reference.first(), alternate.first()) {
        if !is_unambiguous_dna(&[*a]) {
            return false;
        }
        if r != a {
            break;
        }
        reference = &reference[1..];
        alternate = &alternate[1..];
    }
    while let (Some(r), Some(a)) = (reference.last(), alternate.last()) {
        if !is_unambiguous_dna(&[*a]) {
            return false;
        }
        if r != a {
            break;
        }
        reference = &reference[..reference.len() - 1];
        alternate = &alternate[..alternate.len() - 1];
    }
    true
}

/// Whether every allele of a record of two or more alternate alleles opens
/// with a known base (an empty allele opens with none). The engines cut the
/// first base of such a record where its reference allele and every
/// alternate allele open with one letter (Parser/VCF.pm:302-321), so an
/// unknown first letter of any of them would decide the span of all.
pub(crate) fn alleles_open_with_known_bases<'a>(
    alleles: impl IntoIterator<Item = &'a [u8]>,
) -> bool {
    alleles.into_iter().all(|allele| {
        allele
            .first()
            .is_none_or(|base| *base == b'-' || is_unambiguous_dna(&[*base]))
    })
}

/// Whether a stop codon stands among the whole codons of `codon`, an
/// allele's codon string, before the first of them that holds a letter
/// outside A, C, G and T. The protein then ends at that stop whatever the
/// unknown bases after it are, so the allele has a peptide: the one every
/// completion of its unknown bases shares through the stop. Ensembl VEP
/// translates no allele that holds such a letter
/// (TranscriptVariationAllele.pm:701).
pub(crate) fn stop_codon_before_unknown_bases(codon: &[u8], table: u8) -> bool {
    for triplet in codon.chunks_exact(3) {
        if !is_unambiguous_dna(triplet) {
            return false;
        }
        if translate_codon_with_table(triplet, table) == b'*' {
            return true;
        }
    }
    false
}

struct NormalizedTranscriptEdit {
    cds_idx: usize,
    ref_seq: Vec<u8>,
    alt_seq: Vec<u8>,
    trimmed_transcript_flanks: bool,
}

fn normalize_indel_edit_in_transcript_sense(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
) -> Option<NormalizedTranscriptEdit> {
    if cds_pos == 0 {
        return None;
    }

    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    // Ensembl VEP's `TranscriptVariationAllele::peptide` / `codon` paths bail out unless
    // the allele sequence is unambiguous DNA. When these helpers return undef,
    // consequence assignment falls back to coding_sequence_variant rather than
    // peptide-based indel terms.
    if (!ref_is_dash && !is_unambiguous_dna(ref_allele))
        || (!alt_is_dash && !is_unambiguous_dna(alt_allele))
    {
        return None;
    }

    let ref_seq = feature_seq(ref_allele, variant.strand, transcript.strand);
    let alt_seq = feature_seq(alt_allele, variant.strand, transcript.strand);

    let raw_ref_len = ref_seq.len();
    let mut cds_idx = (cds_pos - 1) as usize;
    if transcript.strand == vep_core::coordinate::Strand::Reverse && raw_ref_len > 0 {
        cds_idx = cds_idx.checked_sub(raw_ref_len - 1)?;
    }

    let min_len = ref_seq.len().min(alt_seq.len());
    let mut prefix = 0usize;
    while prefix < min_len && ref_seq[prefix] == alt_seq[prefix] {
        prefix += 1;
    }

    let mut suffix = 0usize;
    while suffix < (min_len - prefix)
        && ref_seq[ref_seq.len() - 1 - suffix] == alt_seq[alt_seq.len() - 1 - suffix]
    {
        suffix += 1;
    }

    cds_idx = cds_idx.checked_add(prefix)?;
    let ref_seq = ref_seq[prefix..ref_seq.len().saturating_sub(suffix)].to_vec();
    let alt_seq = alt_seq[prefix..alt_seq.len().saturating_sub(suffix)].to_vec();

    Some(NormalizedTranscriptEdit {
        cds_idx,
        ref_seq,
        alt_seq,
        trimmed_transcript_flanks: prefix > 0 || suffix > 0,
    })
}

pub(crate) fn translate_to_peptide(cds: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(cds.len() / 3);
    for codon in cds.as_chunks::<3>().0 {
        out.push(translate_codon(codon));
    }
    out
}

/// Translate a local codon-window sequence to a peptide, adding a trailing `X` when a partial
/// codon remains and the peptide is not exactly `*` (Ensembl VEP behavior).
fn translate_codon_window_to_peptide(cds_window: &[u8]) -> Vec<u8> {
    let mut pep = Vec::with_capacity(cds_window.len().div_ceil(3));
    let (codons, rem) = cds_window.as_chunks::<3>();
    for codon in codons {
        pep.push(translate_codon(codon));
    }
    if !rem.is_empty() {
        // Ensembl VEP: append X for partial codons unless peptide is exactly "*"
        let is_exact_stop = pep.len() == 1 && pep[0] == b'*';
        if !is_exact_stop {
            pep.push(b'X');
        }
    }
    pep
}

fn trim_common_peptide_flanks<'a>(ref_pep: &'a [u8], alt_pep: &'a [u8]) -> (&'a [u8], &'a [u8]) {
    let mut prefix = 0usize;
    let min_len = ref_pep.len().min(alt_pep.len());
    while prefix < min_len && ref_pep[prefix] == alt_pep[prefix] {
        prefix += 1;
    }

    let mut suffix = 0usize;
    while suffix < (min_len - prefix)
        && ref_pep[ref_pep.len() - 1 - suffix] == alt_pep[alt_pep.len() - 1 - suffix]
    {
        suffix += 1;
    }

    let ref_end = ref_pep.len().saturating_sub(suffix);
    let alt_end = alt_pep.len().saturating_sub(suffix);
    (&ref_pep[prefix..ref_end], &alt_pep[prefix..alt_end])
}

pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| vep_core::codon::complement_base(b))
        .collect()
}

/// The bases of an allele read against the transcript, as `feature_seq` reads
/// them (VariationFeatureOverlapAllele.pm:246-264): the allele as given when
/// the variant's strand is the transcript's, its reverse complement when the two
/// differ and `seq_is_dna` admits it (`:253`, `:315-321`), and no bases for `-`
/// or an empty allele. The variant's strand is the one its alleles are written
/// on (an Ensembl-format or region record names it); its coordinates are
/// genomic whatever the strand, so the CDS index of an edit follows the
/// transcript's strand alone.
pub(crate) fn feature_seq<'a>(
    allele: &'a [u8],
    variant_strand: Strand,
    transcript_strand: Strand,
) -> Cow<'a, [u8]> {
    if allele == b"-" {
        Cow::Borrowed(&[])
    } else if variant_strand != transcript_strand && seq_is_dna(allele) {
        Cow::Owned(reverse_complement(allele))
    } else {
        Cow::Borrowed(allele)
    }
}

/// Rotates `allele`, written on `allele_strand`, as a 3' shift of `shift` bases
/// along a transcript on `transcript_strand` moves it (`shift_feature_seqs`,
/// TranscriptVariationAllele.pm:1330-1345): read on the transcript's strand,
/// the first base moves to the end once per base of the shift, so an allele
/// written on the transcript's strand rotates left by the shift modulo its
/// length and one written on the other strand rotates right by it, the
/// reverse complement of that left rotation. A shift of the allele's length, or
/// of any multiple of it, leaves the allele as written; a longer shift rotates
/// by the remainder, so the bases placed at the shifted position are the ones
/// the variant produces there. Ensembl VEP counts the places from the
/// transcript's strand alone, `shift_length` on a forward-strand transcript and
/// `seq_length - shift_length` on a reverse-strand one (`:1337-1338`), so it
/// rotates an allele written on the reverse strand the wrong way and, along a
/// reverse-strand transcript, rotates nothing once the shift exceeds the
/// allele, two defects that splice bases the variant never produces. Returns
/// the number of left-rotation steps applied, 0 for an allele left as written.
pub(crate) fn rotate_shifted_allele(
    allele: &mut [u8],
    shift: u64,
    allele_strand: Strand,
    transcript_strand: Strand,
) -> usize {
    let length = allele.len() as u64;
    if length == 0 {
        return 0;
    }
    let remainder = (shift % length) as usize;
    if remainder == 0 {
        return 0;
    }
    let steps = if allele_strand == transcript_strand {
        remainder
    } else {
        allele.len() - remainder
    };
    allele.rotate_left(steps);
    steps
}

/// Compute the 5' UTR sequence, preferring cached UTR (`vefc.five_prime_utr`)
/// from the JSON cache and falling back to FASTA only when the cache lacks
/// the field.
///
/// Ensembl VEP's `_ins_del_start_altered` builds `utr_and_translateable = 5'UTR || CDS`
/// and applies the edit at the cDNA anchor before checking whether the ATG
/// triplet at `length($utr->seq)` is still intact. Ensembl VEP reads the UTR from
/// `_variation_effect_feature_cache->{five_prime_utr}` (a `Bio::EnsEMBL::Slice`
/// object whose `->seq()` reads from the cached `primary_seq` region in
/// offline+cache mode). vep-rs mirrors this by extracting the seq string
/// from the Slice during Storable→JSON conversion (see
/// `tools/vep-storable-to-json`, a sequence value) and using
/// `vefc.five_prime_utr` directly, keeping cache-only mode self-contained.
///
/// `fasta` is optional. When `None` and the cache lacks the field, returns
/// `None` and the caller takes its own fallback path.
pub(crate) fn compute_five_prime_utr_sequence(
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<Vec<u8>> {
    let vefc = transcript.vefc.as_ref()?;

    // Ensembl VEP reads the UTR from the cache via a Slice; the extracted seq string is
    // the offline-mode equivalent.
    if let Some(cached) = vefc.five_prime_utr.as_ref() {
        if !cached.is_empty() {
            return Some(cached.as_bytes().to_vec());
        }
    }

    // FASTA fallback, only when the cache lacks the field.
    let fasta = fasta?;
    let mapper = vefc.mapper.as_ref()?;
    let cdna_coding_start = mapper.cdna_coding_start;
    let pairs = &mapper.exon_coord_mapper.pairs;

    if pairs.is_empty() || cdna_coding_start <= 1 {
        return None;
    }

    let utr_cdna_start = 1u64;
    let utr_cdna_end = cdna_coding_start - 1;

    if utr_cdna_start > utr_cdna_end {
        return None;
    }

    let mut utr_seq = Vec::new();
    for pair in pairs {
        let overlap_start = utr_cdna_start.max(pair.from_start);
        let overlap_end = utr_cdna_end.min(pair.from_end);
        if overlap_start > overlap_end {
            continue;
        }

        let offset_start = overlap_start - pair.from_start;
        let offset_end = overlap_end - pair.from_start;

        let (genomic_start, genomic_end) = if pair.ori == 1 {
            (pair.to_start + offset_start, pair.to_start + offset_end)
        } else {
            (pair.to_end - offset_end, pair.to_end - offset_start)
        };

        let exon_seq = fasta.sequence(&transcript.chr, genomic_start, genomic_end)?;

        if pair.ori == -1 {
            utr_seq.extend(reverse_complement(&exon_seq));
        } else {
            utr_seq.extend(exon_seq);
        }
    }

    if utr_seq.is_empty() {
        None
    } else {
        Some(utr_seq)
    }
}

/// Compute the 3' UTR sequence, preferring cached UTR (`vefc.three_prime_utr`)
/// from the JSON cache and falling back to FASTA only when the cache lacks
/// the field.
///
/// Ensembl VEP's `_ins_del_stop_altered` reads the UTR from
/// `_variation_effect_feature_cache->{three_prime_utr}` (a Slice object whose
/// `->seq()` reads from the cached `primary_seq` region in offline+cache mode).
/// vep-rs mirrors this by extracting the seq string from the Slice during
/// Storable→JSON conversion (`tools/vep-storable-to-json`, a sequence
/// value) and using `vefc.three_prime_utr` directly, keeping
/// cache-only mode self-contained.
///
/// `fasta` is optional. When `None` and the cache lacks the field, returns
/// `None` and the caller takes its own fallback path.
#[doc(hidden)]
pub fn compute_three_prime_utr_sequence(
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<Vec<u8>> {
    let vefc = transcript.vefc.as_ref()?;

    // Ensembl VEP reads the UTR from the cache via a Slice; the extracted seq string is
    // the offline-mode equivalent.
    if let Some(cached) = vefc.three_prime_utr.as_ref() {
        if !cached.is_empty() {
            return Some(cached.as_bytes().to_vec());
        }
    }

    // FASTA fallback, only when the cache lacks the field.
    let fasta = fasta?;
    let mapper = vefc.mapper.as_ref()?;
    let cdna_coding_end = mapper.cdna_coding_end;
    let pairs = &mapper.exon_coord_mapper.pairs;

    if pairs.is_empty() || cdna_coding_end == 0 {
        return None;
    }

    let utr_cdna_start = cdna_coding_end + 1;
    let utr_cdna_end = pairs.last()?.from_end;

    if utr_cdna_start > utr_cdna_end {
        return None;
    }

    let mut utr_seq = Vec::new();
    for pair in pairs {
        let overlap_start = utr_cdna_start.max(pair.from_start);
        let overlap_end = utr_cdna_end.min(pair.from_end);
        if overlap_start > overlap_end {
            continue;
        }

        let offset_start = overlap_start - pair.from_start;
        let offset_end = overlap_end - pair.from_start;

        let (genomic_start, genomic_end) = if pair.ori == 1 {
            (pair.to_start + offset_start, pair.to_start + offset_end)
        } else {
            (pair.to_end - offset_end, pair.to_end - offset_start)
        };

        let exon_seq = fasta.sequence(&transcript.chr, genomic_start, genomic_end)?;

        if pair.ori == -1 {
            // Reverse complement for reverse-strand exons
            utr_seq.extend(reverse_complement(&exon_seq));
        } else {
            utr_seq.extend(exon_seq);
        }
    }

    if utr_seq.is_empty() {
        None
    } else {
        Some(utr_seq)
    }
}

/// Replicate Ensembl VEP's `_ins_del_stop_altered` for deletions that span the
/// CDS into the 3' UTR.
///
/// Ensembl VEP's predicate (`VariationEffect::_ins_del_stop_altered`):
/// ```text
/// $utr_and_translateable = $translateable . $utr_seq;
/// substr($utr_and_translateable, $cds_start - 1, ($cdna_end - $cdna_start) + 1) = $vf_feature_seq;
/// return 1 if length($utr_and_translateable) < length($translateable);
/// # Extract codon at the (original) stop position in the edited sequence and translate.
/// my $pep = translate(substr($utr_and_translateable, length($translateable) - 3, 3));
/// return !($pep && $pep eq '*');
/// ```
///
/// Returns:
/// * `Some(true)`: stop is altered (→ Ensembl VEP emits `stop_lost`)
/// * `Some(false)`: stop is not altered (→ Ensembl VEP emits `stop_retained_variant`)
/// * `None`: cannot evaluate (missing FASTA / UTR / CDS data, non-DNA
///   allele, no stop codon standing at the coding end, etc.); the caller
///   takes its own fallback path.
///
/// Required inputs:
/// * `cds_start`: 1-based CDS position where the edit begins (the `cds_start`
///   used by `substr($_, $cds_start - 1, ...)` in Ensembl VEP)
/// * `cds_edit_len`: number of bases the edit removes within the CDS
///   (for a deletion spanning into UTR, this equals `(cdna_end - cdna_start) + 1`
///   clamped by CDS end, i.e. the Perl `substr` length)
/// * `reference_fasta`: optional FASTA fallback when the JSON cache lacks
///   `vefc.three_prime_utr`. The cache-first path is used by default since
///   Ensembl VEP reads UTR from cache in offline mode (see
///   `compute_three_prime_utr_sequence`).
///
/// Caller must have verified that `_overlaps_stop_codon` is true; otherwise
/// Ensembl VEP's `stop_retained` gate fails regardless of this predicate.
pub(crate) fn replicate_ins_del_stop_altered(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_start: u64,
    cds_edit_len: usize,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<bool> {
    let vefc = transcript.vefc.as_ref()?;
    let translateable_seq = vefc.translateable_seq.as_ref()?;
    let cds = translateable_seq.as_bytes();
    if cds.len() < 3 {
        return None;
    }

    // Ensembl VEP's `_overlaps_stop_codon` refuses `cds_end_NF` alone
    // (VariationEffect.pm:1373); here the test is made only where a stop codon
    // stands at the coding end, and `None` leaves the caller on its fallback.
    if stop_codon_at_the_coding_end(transcript) == StopCodonAtTheCodingEnd::None
        || transcript.flags.iter().any(|f| &**f == "cds_start_NF")
    {
        return None;
    }

    if cds_start == 0 {
        return None;
    }
    let idx = (cds_start - 1) as usize;
    if idx >= cds.len() {
        return None;
    }

    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    // Ensembl VEP's `seq_is_unambiguous_dna` gate: `_ins_del_stop_altered` returns 0
    // when either allele contains non-ACGT characters.
    if (!ref_is_dash && !is_unambiguous_dna(ref_allele))
        || (!alt_is_dash && !is_unambiguous_dna(alt_allele))
    {
        return None;
    }

    let alt_seq = feature_seq(alt_allele, variant.strand, transcript.strand);

    // `compute_three_prime_utr_sequence` already reverse-complements for a
    // reverse-strand transcript, so the concatenation below is in reading frame.
    let utr_seq = compute_three_prime_utr_sequence(transcript, reference_fasta)?;
    if utr_seq.is_empty() {
        return None;
    }

    let mut combined: Vec<u8> = Vec::with_capacity(cds.len() + utr_seq.len());
    combined.extend_from_slice(cds);
    combined.extend_from_slice(&utr_seq);

    // Perl's `substr` in Ensembl VEP silently clamps an edit range past the combined sequence;
    // `None` leaves the caller on its fallback instead.
    if idx + cds_edit_len > combined.len() {
        return None;
    }

    // Ensembl VEP: substr($seq, $cds_start - 1, $length) = $vf_feature_seq;
    combined.splice(idx..idx + cds_edit_len, alt_seq.iter().copied());

    // Ensembl VEP: "return 1 if length($utr_and_translateable) < length($translateable);"
    if combined.len() < cds.len() {
        return Some(true);
    }

    // Ensembl VEP: extract the 3-base codon at the original stop position and translate.
    let stop_idx = cds.len() - 3;
    if stop_idx + 3 > combined.len() {
        return None;
    }
    let codon = &combined[stop_idx..stop_idx + 3];
    if codon.len() != 3 {
        return None;
    }
    let aa = translate_codon(codon);
    // `!($pep && $pep eq '*')`: altered iff the codon is not a stop.
    Some(aa != b'*')
}

/// Replicate Ensembl VEP's `_ins_del_start_altered` for deletions/insertions at
/// the 5'UTR/CDS boundary.
///
/// Ensembl VEP's predicate (`VariationEffect::_ins_del_start_altered`):
/// ```text
/// $utr_and_translateable = ($utr ? $utr->seq : '') . $translateable;
/// substr($utr_and_translateable, $cdna_start - 1, ($cdna_end - $cdna_start) + 1) = $vf_feature_seq;
/// if ($utr) {
///     my $atg_start = length($utr->seq);
///     my $new_sc = substr($utr_and_translateable, $atg_start, 3);
///     my $new_utr = substr($utr_and_translateable, 0, length($utr->seq));
///     return 0 if $new_utr eq $utr->seq && $new_sc eq 'ATG';
/// }
/// return 1 if length($utr_and_translateable) < length($translateable);
/// return $translateable ne substr($utr_and_translateable, 0 - length($translateable));
/// ```
///
/// Returns:
/// * `Some(true)`: start is altered (→ Ensembl VEP emits `start_lost`)
/// * `Some(false)`: start is not altered (→ Ensembl VEP keeps only `5_prime_UTR_variant`)
/// * `None`: cannot evaluate (missing FASTA / UTR / CDS data, non-DNA
///   allele, no start codon standing at the coding start, etc.); the caller
///   takes its own fallback path.
///
/// Required inputs:
/// * `cdna_start`: 1-based cDNA position where the edit begins (anchor for
///   the `substr` splice in the combined UTR+CDS string)
/// * `cdna_edit_len`: number of bases the edit removes from the cDNA (equals
///   Ensembl VEP's `$cdna_end - $cdna_start + 1`)
/// * `reference_fasta`: optional FASTA fallback when the JSON cache lacks
///   `vefc.five_prime_utr`. The cache-first path is used by default since
///   Ensembl VEP reads UTR from cache in offline mode (see
///   `compute_five_prime_utr_sequence`).
///
/// Caller must have verified that `_overlaps_start_codon` is true; otherwise
/// Ensembl VEP's `_ins_del_start_altered` returns 0 regardless of this predicate.
pub(crate) fn replicate_ins_del_start_altered(
    variant: &InputVariant,
    transcript: &Transcript,
    cdna_start: u64,
    cdna_edit_len: usize,
    reference_fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<bool> {
    let vefc = transcript.vefc.as_ref()?;
    let translateable_seq = vefc.translateable_seq.as_ref()?;
    let cds = translateable_seq.as_bytes();
    if cds.len() < 3 {
        return None;
    }

    // Ensembl VEP's `_overlaps_start_codon` refuses `cds_start_NF` alone
    // (VariationEffect.pm:975); here the test is made only where a start codon
    // stands at the coding start. The start codon analysis is meaningless for
    // `cds_end_NF`.
    if !start_codon_stands_at_the_coding_start(transcript)
        || transcript.flags.iter().any(|f| &**f == "cds_end_NF")
    {
        return None;
    }

    if cdna_start == 0 {
        return None;
    }

    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    // Ensembl VEP's `seq_is_unambiguous_dna` gate.
    if (!ref_is_dash && !is_unambiguous_dna(ref_allele))
        || (!alt_is_dash && !is_unambiguous_dna(alt_allele))
    {
        return None;
    }

    let alt_seq = feature_seq(alt_allele, variant.strand, transcript.strand);

    // The helper reverse-complements for a reverse-strand transcript, so the
    // concatenation reads 5' to CDS in reading frame.
    let utr_seq = compute_five_prime_utr_sequence(transcript, reference_fasta)?;
    if utr_seq.is_empty() {
        return None;
    }

    let mut combined: Vec<u8> = Vec::with_capacity(utr_seq.len() + cds.len());
    combined.extend_from_slice(&utr_seq);
    combined.extend_from_slice(cds);

    let idx = (cdna_start - 1) as usize;
    // Guard against an edit range that extends past the combined sequence.
    if idx >= combined.len() || idx + cdna_edit_len > combined.len() {
        return None;
    }

    // Ensembl VEP: substr($utr_and_translateable, $cdna_start - 1, $length) = $vf_feature_seq;
    combined.splice(idx..idx + cdna_edit_len, alt_seq.iter().copied());

    // Ensembl VEP: if UTR remained intact and the 3 bases at ATG position are still
    // 'ATG', start is not altered.
    let atg_start = utr_seq.len();
    if atg_start + 3 <= combined.len() {
        let new_utr = &combined[..atg_start];
        let new_sc = &combined[atg_start..atg_start + 3];
        if new_utr == utr_seq.as_slice() && new_sc == b"ATG" {
            return Some(false);
        }
    }

    // Ensembl VEP: sequence shorter → start has been altered.
    if combined.len() < cds.len() {
        return Some(true);
    }

    // Ensembl VEP: `$translateable ne substr($utr_and_translateable, 0 - length($translateable))`,
    // which compares the last `cds.len()` bytes of the edited sequence against the
    // original CDS. If they differ, start is altered.
    let tail_start = combined.len().saturating_sub(cds.len());
    let tail = &combined[tail_start..];
    Some(tail != cds)
}

/// Compute full peptide alleles (ref_pep, alt_pep) for a coding variant.
///
/// Applies the variant to the CDS and translates the full reference and
/// alternate CDS sequences without trimming. This is useful for predicates
/// that need the absolute stop-codon position rather than the minimal local
/// amino-acid change.
pub(crate) fn compute_full_peptide_alleles(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
) -> Option<(Vec<u8>, Vec<u8>)> {
    compute_full_peptide_alleles_impl(variant, transcript, cds_pos, false)
}

/// Like `compute_full_peptide_alleles` but clamps the reference allele length
/// to available CDS bounds instead of rejecting when `idx + ref_len > cds.len()`.
///
/// Used for UTR-spanning deletions where the deletion extends past CDS end.
/// Clamping yields the peptide effect on the CDS portion, allowing
/// stop_lost vs stop_retained classification.
pub(crate) fn compute_full_peptide_alleles_clamped(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
) -> Option<(Vec<u8>, Vec<u8>)> {
    compute_full_peptide_alleles_impl(variant, transcript, cds_pos, true)
}

fn compute_full_peptide_alleles_impl(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
    clamp_ref_len: bool,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let vefc = transcript.vefc.as_ref()?;
    let translateable_seq = vefc.translateable_seq.as_ref()?;
    let cds = translateable_seq.as_bytes();

    if cds_pos == 0 {
        return None;
    }
    let idx = (cds_pos - 1) as usize;
    if idx >= cds.len() {
        return None;
    }

    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    if (!ref_is_dash && !is_unambiguous_dna(ref_allele))
        || (!alt_is_dash && !is_unambiguous_dna(alt_allele))
    {
        return None;
    }

    let alt_seq = feature_seq(alt_allele, variant.strand, transcript.strand);

    let nominal_ref_len = if ref_is_dash { 0 } else { ref_allele.len() };
    let ref_len = if clamp_ref_len {
        nominal_ref_len.min(cds.len() - idx)
    } else {
        if idx + nominal_ref_len > cds.len() {
            return None;
        }
        nominal_ref_len
    };

    let mut alt_cds = Vec::with_capacity(cds.len().saturating_sub(ref_len) + alt_seq.len());
    alt_cds.extend_from_slice(&cds[..idx]);
    alt_cds.extend_from_slice(&alt_seq);
    alt_cds.extend_from_slice(&cds[idx + ref_len..]);

    Some((translate_to_peptide(cds), translate_to_peptide(&alt_cds)))
}

/// Compute trimmed peptide alleles (ref_pep, alt_pep) for a coding variant.
///
/// Applies the variant to the CDS, translates both, and trims common
/// prefix/suffix so that the returned slices represent the minimal
/// amino-acid change, mirroring Ensembl VEP's `_get_peptide_alleles`.
///
/// Returns `None` if the CDS or variant cannot be resolved.
pub(crate) fn compute_peptide_alleles(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let (ref_pep, alt_pep) = compute_full_peptide_alleles(variant, transcript, cds_pos)?;
    let (r, a) = trim_common_peptide_flanks(&ref_pep, &alt_pep);
    Some((r.to_vec(), a.to_vec()))
}

/// Compute peptide alleles over the variant's codon-window, in the style of Ensembl VEP's
/// `TranscriptVariationAllele::codon` / `peptide`.
///
/// This differs from `compute_peptide_alleles` (full-CDS translation) and matches
/// Ensembl VEP's frameshift `stop_gained` behavior, where the gained stop can occur in a
/// later codon inside the local window.
///
/// Ensembl VEP's `codon` method:
///   - Builds the full alternate CDS via `_get_alternate_cds()`
///   - Extracts: `substr($alt_cds, $codon_cds_start-1, $codon_len + ($allele_len - $vf_nt_len))`
///     where `codon_len` = ref codon window size, `allele_len` = alt allele length,
///     `vf_nt_len` = ref allele length in CDS terms.
///
/// Then `peptide` translates full codons and appends `X` for partial remainders.
///
/// For Ensembl VEP parity, callers should derive `bounds` using
/// `mapper::map_genomic_span_to_cds_bounds()` rather than mapping individual endpoints.
pub(crate) fn compute_codon_window_peptide_alleles(
    variant: &InputVariant,
    transcript: &Transcript,
    bounds: &CdsSpanBounds,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let vefc = transcript.vefc.as_ref()?;
    let translateable_seq = vefc.translateable_seq.as_ref()?;
    let cds = translateable_seq.as_bytes();
    if cds.is_empty() {
        return None;
    }

    let alt_allele = variant.alt_allele();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();
    let ref_allele = &variant.ref_allele;
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();

    if (!ref_is_dash && !is_unambiguous_dna(ref_allele))
        || (!alt_is_dash && !is_unambiguous_dna(alt_allele))
    {
        return None;
    }

    let alt_seq = feature_seq(alt_allele, variant.strand, transcript.strand);
    // Ensembl VEP: $allele_len = length of alt variation_feature_seq (in mRNA sense)
    let allele_len = alt_seq.len();

    // Ensembl VEP: $vf_nt_len = cds_end - cds_start + 1.
    // For insertions (cds_start > cds_end in Ensembl VEP), this evaluates to 0.
    let cds_start = bounds.cds_start;
    if cds_start == 0 {
        return None;
    }
    let vf_nt_len_isize = (bounds.cds_end as isize) - (bounds.cds_start as isize) + 1;
    let vf_nt_len: usize = vf_nt_len_isize.max(0) as usize;
    let is_insertion_bounds = bounds.cds_start > bounds.cds_end;
    let is_boundary_insertion =
        is_insertion_bounds && bounds.translation_start != bounds.translation_end;

    // Ensembl VEP: a codon-boundary insertion spans two adjacent peptide coordinates,
    // yielding an empty ref codon window and an alt-only inserted codon/peptide.
    // A same-codon insertion keeps the normal window and can remain
    // protein_altering_variant.
    let (window_start_idx, ref_window_end, alt_window_len) = if is_boundary_insertion {
        let start_idx = (bounds.cds_start - 1) as usize;
        if start_idx > cds.len() {
            return None;
        }
        (start_idx, start_idx, allele_len)
    } else {
        // Ensembl VEP: $codon_cds_start = $tv_tr_start * 3 - 2;  $codon_cds_end = $tv_tr_end * 3;
        let tr_start = bounds.translation_start;
        let tr_end = bounds.translation_end;
        if tr_start == 0 || tr_end == 0 || tr_end < tr_start {
            return None;
        }
        let codon_cds_start = tr_start * 3 - 2; // 1-based CDS position
        let codon_cds_end = tr_end * 3; // 1-based CDS position (inclusive)
        let codon_len = (codon_cds_end - codon_cds_start + 1) as usize;
        let start_idx = (codon_cds_start - 1) as usize; // 0-based
        let end_idx = codon_cds_end as usize; // exclusive
        if start_idx >= cds.len() {
            return None;
        }
        let alt_len =
            ((codon_len as isize) + (allele_len as isize) - (vf_nt_len as isize)).max(0) as usize;
        (start_idx, end_idx.min(cds.len()), alt_len)
    };

    // Ensembl VEP `_get_alternate_cds` splicing:
    //   upstream   = cds[0..cds_start-1]
    //   downstream = cds[cds_end..]   (0-based index, so removes cds_start..cds_end inclusive)
    //
    // For insertions, Ensembl VEP has cds_start > cds_end so downstream starts at cds_start-1 and
    // the inserted sequence is placed at index cds_start-1.
    let splice_start: usize = (cds_start - 1) as usize;
    let splice_end: usize = bounds.cds_end as usize;
    if splice_start > cds.len() || splice_end > cds.len() || splice_end < splice_start {
        return None;
    }

    let mut alt_cds = Vec::with_capacity(cds.len().saturating_sub(vf_nt_len) + allele_len);
    alt_cds.extend_from_slice(&cds[..splice_start]);
    alt_cds.extend_from_slice(&alt_seq);
    alt_cds.extend_from_slice(&cds[splice_end..]);

    // This codon window stops at the CDS end and appends no 3' UTR. The cached UTR
    // (`vefc.three_prime_utr`, see `compute_three_prime_utr_sequence`) is appended
    // instead by `PerlCodingEval::alternate_cds`, the whole port of Ensembl VEP's
    // `_get_alternate_cds`.

    if ref_window_end < window_start_idx {
        return None;
    }
    let ref_window = &cds[window_start_idx..ref_window_end];

    if window_start_idx + alt_window_len > alt_cds.len() {
        return None;
    }
    let alt_window = &alt_cds[window_start_idx..window_start_idx + alt_window_len];

    let ref_pep = translate_codon_window_to_peptide(ref_window);
    let alt_pep = translate_codon_window_to_peptide(alt_window);
    Some((ref_pep, alt_pep))
}

/// Compute the frameshift stop_gained scan limit for the full-CDS fallback path.
///
/// The full-CDS peptide can extend far past Ensembl VEP's local codon window, so the
/// scan is capped to avoid a stop_gained for a stop that Ensembl VEP's narrower
/// codon-window peptide would never reach.
///
/// For net insertions, cap at `1 + floor(net_insertion_nt / 3)` amino acids
/// (approximating the number of fully-formed codons Ensembl VEP's `_get_alternate_cds`
/// window produces from the insertion; partial codons at the boundary
/// translate to 'X' in Ensembl VEP and are not stops).
///
/// For net deletions the window is naturally bounded by the CDS length, so
/// scan the entire alt peptide.
fn frameshift_scan_limit(variant: &InputVariant, alt_pep_len: usize) -> usize {
    let ref_len_nt = if variant.ref_allele == b"-" || variant.ref_allele.is_empty() {
        0
    } else {
        variant.ref_allele.len()
    };
    let alt_allele = variant.alt_allele();
    if alt_allele.len() > ref_len_nt {
        let net_ins_nt = alt_allele.len() - ref_len_nt;
        (1 + net_ins_nt / 3).min(alt_pep_len)
    } else {
        alt_pep_len
    }
}

pub(crate) fn frameshift_stop_gained_in_codon_window(
    variant: &InputVariant,
    transcript: &Transcript,
    bounds: &CdsSpanBounds,
) -> bool {
    let Some((ref_pep, alt_pep)) =
        compute_codon_window_peptide_alleles(variant, transcript, bounds)
    else {
        return false;
    };
    if ref_pep.contains(&b'*') {
        return false;
    }
    // Ensembl VEP `stop_gained` returns true when the alt peptide contains `*` and the
    // ref peptide does not. Ensembl VEP applies no position cap:
    //   $cache->{stop_gained} = (($alt_pep =~ /\*/) and ($ref_pep !~ /\*/));
    //
    // The codon-window peptide is already clipped to Ensembl VEP's window width by
    // `compute_codon_window_peptide_alleles` + `translate_codon_window_to_peptide`
    // (partial boundary codons become 'X' unless the whole-codon prefix is
    // exactly '*'). So scanning the full alt peptide here matches Ensembl VEP exactly.
    alt_pep.contains(&b'*')
}

/// Fallback frameshift stop_gained predicate using full-CDS translation.
///
/// Used when `map_genomic_span_to_cds_bounds()` returns `None` and the
/// codon-window approach is unavailable. Translates the entire ref and alt
/// CDS, then checks the region starting from the variant's amino acid
/// position with the scan-limit logic to avoid finding stops far downstream
/// that Ensembl VEP's narrower codon window would never see.
pub(crate) fn frameshift_stop_gained_full_cds(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
) -> bool {
    let Some((ref_pep, alt_pep)) = compute_full_peptide_alleles(variant, transcript, cds_pos)
    else {
        return false;
    };

    if cds_pos == 0 {
        return false;
    }
    let variant_aa_idx = ((cds_pos - 1) / 3) as usize;
    if variant_aa_idx >= ref_pep.len() || variant_aa_idx >= alt_pep.len() {
        return false;
    }

    let ref_local = &ref_pep[variant_aa_idx..];
    let alt_local = &alt_pep[variant_aa_idx..];

    let scan_limit = frameshift_scan_limit(variant, alt_local.len());

    // If the ref already has a stop in the scan region (e.g., variant is near
    // the terminal stop codon), do not report stop_gained.
    let ref_check_end = scan_limit.min(ref_local.len());
    if ref_local[..ref_check_end].contains(&b'*') {
        return false;
    }

    alt_local[..scan_limit].contains(&b'*')
}

/// Get the codon change for a coding variant.
///
/// For SNVs, extracts the affected codon from the translateable sequence,
/// creates the mutant codon, translates both, and formats the display
/// strings with the variant base in uppercase.
///
/// For indels, determines whether the change is a frameshift (length
/// change not divisible by 3).
pub fn get_codon_change(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
) -> Option<CodonChange> {
    let vefc = transcript.vefc.as_ref()?;
    let translateable_seq = vefc.translateable_seq.as_ref()?;
    let seq_bytes = translateable_seq.as_bytes();

    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();

    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();

    if !ref_is_dash && !alt_is_dash && ref_allele.len() == alt_allele.len() && ref_allele.len() > 1
    {
        return get_mnv_codon_change(variant, seq_bytes, cds_pos, transcript.strand);
    }

    if ref_is_dash || alt_is_dash || ref_allele.len() != 1 || alt_allele.len() != 1 {
        return get_indel_codon_change(variant, seq_bytes, cds_pos, transcript.strand);
    }

    let mut cc = get_snv_codon_change(
        variant,
        seq_bytes,
        cds_pos,
        transcript.strand,
        codon_table_for(transcript),
    )?;
    if !vefc.seq_edits.is_empty() {
        let mut ref_aa = vec![cc.ref_amino_acid];
        let pos = cc.codon_number as i64;
        if apply_seq_edits(&mut ref_aa, &vefc.seq_edits, pos, pos) {
            cc.ref_amino_acid = ref_aa[0];
        }
    } else if let Some(pep_seq) = vefc.peptide.as_deref() {
        let mut ref_aa = [cc.ref_amino_acid];
        overlay_seq_edits(
            &mut ref_aa,
            pep_seq.as_bytes(),
            cc.codon_number as i64,
            seq_bytes.len(),
        );
        cc.ref_amino_acid = ref_aa[0];
    }
    // The reference row of `--process_ref_homs` reads the edited peptide for its
    // allele as well (TranscriptVariationAllele.pm:716, 742-770).
    if variant.reference_allele_row {
        cc.alt_amino_acid = cc.ref_amino_acid;
    }
    Some(cc)
}

/// NCBI translation table for a transcript, defaulting to 1 as Ensembl VEP does.
///
/// Ensembl VEP reads this per transcript in `BaseTranscriptVariation::_codon_table`
/// (`$attrib ? $attrib->value : 1`), the `codon_table` attribute of the
/// transcript's sequence region: 2 on a vertebrate's mitochondrion, 5 on an
/// invertebrate's, 3 on a yeast's, 11 on a plastid, 1 elsewhere.
pub(crate) fn codon_table_for(transcript: &Transcript) -> u8 {
    transcript
        .vefc
        .as_ref()
        .map_or(1, |v| if v.codon_table == 0 { 1 } else { v.codon_table })
}

fn get_snv_codon_change(
    variant: &InputVariant,
    seq_bytes: &[u8],
    cds_pos: u64,
    transcript_strand: vep_core::coordinate::Strand,
    codon_table: u8,
) -> Option<CodonChange> {
    if cds_pos == 0 {
        return None;
    }

    let cds_idx = (cds_pos - 1) as usize; // 0-based index into CDS
    let codon_number = cds_idx / 3 + 1; // 1-based codon number
    let pos_in_codon = cds_idx % 3; // 0, 1, or 2
    let codon_start = cds_idx - pos_in_codon;

    if codon_start + 3 > seq_bytes.len() {
        return None;
    }

    let ref_codon_bytes = &seq_bytes[codon_start..codon_start + 3];
    let mut alt_codon_bytes = ref_codon_bytes.to_vec();

    // The one base as `feature_seq` reads it, complemented in place rather than
    // through a one-byte allocation on this hottest coding path.
    let alt_base = variant.alt_allele()[0];
    let effective_alt = if variant.strand != transcript_strand {
        vep_core::codon::complement_base(alt_base)
    } else {
        alt_base
    };

    alt_codon_bytes[pos_in_codon] = effective_alt.to_ascii_uppercase();

    let ref_aa = translate_codon_with_table(ref_codon_bytes, codon_table);
    let alt_aa = translate_codon_with_table(&alt_codon_bytes, codon_table);

    let ref_display = format_codon_display(ref_codon_bytes, pos_in_codon);
    let alt_display = format_codon_display(&alt_codon_bytes, pos_in_codon);

    Some(CodonChange {
        ref_codon: ref_display,
        alt_codon: alt_display,
        ref_amino_acid: ref_aa,
        alt_amino_acid: alt_aa,
        codon_number: codon_number as u64,
        is_frameshift: false,
    })
}

/// Compute codon change for a multi-nucleotide variant (MNV).
///
/// MNVs are same-length multi-base substitutions (e.g., ref=AC, alt=TG).
/// These may span one or more codons. The full codon window is translated for
/// ref and alt, and the representative AA pair for consequence classification is
/// chosen by priority (stop_gained > stop_lost > missense > synonymous).
fn get_mnv_codon_change(
    variant: &InputVariant,
    seq_bytes: &[u8],
    cds_pos: u64,
    transcript_strand: vep_core::coordinate::Strand,
) -> Option<CodonChange> {
    if cds_pos == 0 {
        return None;
    }

    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();
    let ref_len = ref_allele.len();

    // On the reverse strand cds_pos from the start genomic position is the last
    // CDS base of the variant, so it is moved back to the first affected position.
    let cds_idx = if transcript_strand == vep_core::coordinate::Strand::Reverse {
        let raw = (cds_pos - 1) as usize;
        raw.checked_sub(ref_len - 1)?
    } else {
        (cds_pos - 1) as usize
    };
    let codon_number = cds_idx / 3 + 1;

    let codon_start = (cds_idx / 3) * 3;
    let codon_end = ((cds_idx + ref_len).div_ceil(3) * 3).min(seq_bytes.len());

    if codon_end <= codon_start || cds_idx + ref_len > seq_bytes.len() {
        return None;
    }

    let ref_window = &seq_bytes[codon_start..codon_end];

    let alt_seq = feature_seq(alt_allele, variant.strand, transcript_strand);

    let rel = cds_idx - codon_start;
    let mut alt_window = ref_window.to_vec();
    for (i, &b) in alt_seq.iter().enumerate() {
        if rel + i < alt_window.len() {
            alt_window[rel + i] = b;
        }
    }

    let ref_pep = translate_to_peptide(ref_window);
    let alt_pep = translate_to_peptide(&alt_window);

    if ref_pep.is_empty() || alt_pep.is_empty() {
        return None;
    }

    // Representative AA pair priority: stop change > first non-stop difference >
    // synonymous.
    let min_len = ref_pep.len().min(alt_pep.len());
    let mut rep_ref_aa = ref_pep[0];
    let mut rep_alt_aa = alt_pep[0];
    let mut found = false;

    for i in 0..min_len {
        if alt_pep[i] == b'*' && ref_pep[i] != b'*' {
            rep_ref_aa = ref_pep[i];
            rep_alt_aa = b'*';
            found = true;
            break;
        }
        if ref_pep[i] == b'*' && alt_pep[i] != b'*' {
            rep_ref_aa = b'*';
            rep_alt_aa = alt_pep[i];
            found = true;
            break;
        }
    }

    if !found {
        for i in 0..min_len {
            if ref_pep[i] != alt_pep[i] {
                rep_ref_aa = ref_pep[i];
                rep_alt_aa = alt_pep[i];
                found = true;
                break;
            }
        }
    }

    if !found {
        rep_ref_aa = ref_pep[0];
        rep_alt_aa = alt_pep[0];
    }

    let ref_display = format_mnv_codon_display(ref_window, rel, ref_len);
    let alt_display = format_mnv_codon_display(&alt_window, rel, ref_len);

    Some(CodonChange {
        ref_codon: ref_display,
        alt_codon: alt_display,
        ref_amino_acid: rep_ref_aa,
        alt_amino_acid: rep_alt_aa,
        codon_number: codon_number as u64,
        is_frameshift: false,
    })
}

/// Format a multi-base codon window for display: lowercase with variant positions uppercase.
fn format_mnv_codon_display(window: &[u8], variant_start: usize, variant_len: usize) -> String {
    let mut result = String::with_capacity(window.len());
    for (i, &b) in window.iter().enumerate() {
        if i >= variant_start && i < variant_start + variant_len {
            result.push(b.to_ascii_uppercase() as char);
        } else {
            result.push(b.to_ascii_lowercase() as char);
        }
    }
    result
}

/// Compute codon change for an indel.
///
/// For frameshifts, also computes the first alt amino acid when possible,
/// so that `stop_gained` can be detected (Ensembl VEP parity: `tat/tAat` -> `Y/*`).
fn get_indel_codon_change(
    variant: &InputVariant,
    seq_bytes: &[u8],
    cds_pos: u64,
    transcript_strand: vep_core::coordinate::Strand,
) -> Option<CodonChange> {
    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();

    let ref_len = if ref_allele == b"-" || ref_allele.is_empty() {
        0usize
    } else {
        ref_allele.len()
    };
    let alt_len = if alt_allele == b"-" || alt_allele.is_empty() {
        0usize
    } else {
        alt_allele.len()
    };

    let is_frameshift = ref_len.abs_diff(alt_len) % 3 != 0;

    if cds_pos == 0 {
        return None;
    }

    let cds_idx = (cds_pos - 1) as usize;
    let codon_number = cds_idx / 3 + 1;
    let pos_in_codon = cds_idx % 3;
    let codon_start = cds_idx - pos_in_codon;

    let ref_aa = if codon_start + 3 <= seq_bytes.len() {
        translate_codon(&seq_bytes[codon_start..codon_start + 3])
    } else {
        b'X'
    };

    // The first alt codon of a frameshift detects stop_gained (Ensembl VEP: `tat/tAat`
    // -> `Y/*`).
    let alt_aa = if is_frameshift {
        compute_first_alt_aa_for_indel(seq_bytes, cds_idx, ref_len, variant, transcript_strand)
            .unwrap_or(b'X')
    } else if ref_len == 0 && alt_len > 0 {
        ref_aa
    } else {
        b'X'
    };

    Some(CodonChange {
        ref_codon: if codon_start + 3 <= seq_bytes.len() {
            std::str::from_utf8(&seq_bytes[codon_start..codon_start + 3])
                .unwrap_or("???")
                .to_owned()
        } else {
            String::new()
        },
        alt_codon: String::new(),
        ref_amino_acid: ref_aa,
        alt_amino_acid: alt_aa,
        codon_number: codon_number as u64,
        is_frameshift,
    })
}

/// Compute the first alt amino acid after applying an indel to the CDS.
///
/// Builds the alt CDS locally around the affected codon, then translates
/// the first codon from the variant position. This detects
/// stop_gained/stop_lost in frameshift contexts.
fn compute_first_alt_aa_for_indel(
    cds: &[u8],
    cds_idx: usize,
    ref_len: usize,
    variant: &InputVariant,
    transcript_strand: vep_core::coordinate::Strand,
) -> Option<u8> {
    let pos_in_codon = cds_idx % 3;
    let codon_start = cds_idx - pos_in_codon;
    let codon_end = codon_start.saturating_add(3);

    // A codon needs at least 3 bases after the edit.
    let ref_is_dash = variant.ref_allele == b"-" || variant.ref_allele.is_empty();
    let alt_bases = feature_seq(variant.alt_allele(), variant.strand, transcript_strand);
    let del_len = if ref_is_dash { 0 } else { ref_len };

    if codon_end > cds.len() || cds_idx + del_len > cds.len() {
        return None;
    }

    // Ensembl VEP's frameshift peptide derives from its local `codon()`, which does not
    // pad a short codon with downstream bases after a deletion (`Tta/ta` yields
    // alt peptide `X`), so only bases within the current codon window are
    // stitched: [codon_start..cds_idx] + alt_bases + [post-deletion..codon_end].
    let prefix = &cds[codon_start..cds_idx];
    let del_end_within_codon = (cds_idx + del_len).min(codon_end);
    let suffix_within_codon = &cds[del_end_within_codon..codon_end];

    let mut local_alt =
        Vec::with_capacity(prefix.len() + alt_bases.len() + suffix_within_codon.len());
    local_alt.extend_from_slice(prefix);
    local_alt.extend_from_slice(&alt_bases);
    local_alt.extend_from_slice(suffix_within_codon);

    if local_alt.len() >= 3 {
        Some(translate_codon(&local_alt[0..3]))
    } else {
        None
    }
}

/// Format a codon for display: lowercase bases with the variant position in uppercase.
fn format_codon_display(codon: &[u8], variant_pos: usize) -> String {
    let mut result = String::with_capacity(3);
    for (i, &b) in codon.iter().enumerate() {
        if i == variant_pos {
            result.push(b.to_ascii_uppercase() as char);
        } else {
            result.push(b.to_ascii_lowercase() as char);
        }
    }
    result
}

// Ensembl VEP coding predicates, ported whole.
//
// `perl_coding_terms` reproduces the coding SO terms Ensembl VEP assigns to a
// non-SNV allele on a coding transcript: the `Bio::EnsEMBL::TranscriptMapper`
// projection of the allele (cDNA, CDS and peptide coordinates, including the
// `Bio::EnsEMBL::Mapper::map_insert` treatment of an insertion), the
// `TranscriptVariationAllele::codon` / `peptide` window (alternate CDS with the
// 3' UTR appended, window `codon_len + allele_len - vf_nt_len`, `X` for a
// partial codon), and the predicates of `Utils::VariationEffect` with their
// cache-placeholder recursion (`start_lost` <-> `inframe_insertion`).

use std::borrow::Cow;

use vep_core::consequence::Consequence;
use vep_core::coordinate::Strand;
use vep_core::transcript::{BamEditStatus, ExonCoordMapper, MapperPair, SeqEdit};

/// One element of a `Bio::EnsEMBL::Mapper` result list, in transcript order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PerlMapSeg {
    Coord { start: i64, end: i64 },
    Gap { start: i64, end: i64 },
}

/// `BaseTranscriptVariation` coordinate accessors for one allele span.
///
/// `None` is Perl's `undef`: the first (or last) mapped segment was a gap.
/// `cds_start > cds_end` (and `tl_start > tl_end`) is Ensembl VEP's insertion shape.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PerlSpan {
    pub cdna_start: Option<i64>,
    pub cdna_end: Option<i64>,
    pub cds_start: Option<i64>,
    pub cds_end: Option<i64>,
    pub tl_start: Option<i64>,
    pub tl_end: Option<i64>,
    pub cds_coords: Vec<PerlMapSeg>,
}

fn genomic_order(ecm: &ExonCoordMapper) -> Vec<usize> {
    if ecm.by_genomic_start().len() == ecm.pairs.len() {
        ecm.by_genomic_start().iter().map(|&i| i as usize).collect()
    } else {
        let mut idx: Vec<usize> = (0..ecm.pairs.len()).collect();
        idx.sort_by_key(|&i| ecm.pairs[i].to_start);
        idx
    }
}

/// `Bio::EnsEMBL::Mapper::map_coordinates` (genomic -> cDNA) for a non-insertion span.
fn perl_map_span(ecm: &ExonCoordMapper, start: i64, end: i64, strand: i8) -> Vec<PerlMapSeg> {
    let mut out = Vec::new();
    let mut cur = start;
    let mut last_to_end: Option<i64> = None;
    for i in genomic_order(ecm) {
        let p = &ecm.pairs[i];
        let (ts, te) = (p.to_start as i64, p.to_end as i64);
        if te < start || ts > end {
            continue;
        }
        if cur < ts {
            out.push(PerlMapSeg::Gap {
                start: cur,
                end: ts - 1,
            });
            cur = ts;
        }
        let (fs, fe) = (p.from_start as i64, p.from_end as i64);
        let (mut tstart, mut tend) = (0i64, 0i64);
        if p.ori == 1 {
            tstart = fs + (cur - ts);
        } else {
            tend = fe - (cur - ts);
        }
        if end > te {
            if p.ori == 1 {
                tend = fe;
            } else {
                tstart = fs;
            }
        } else if p.ori == 1 {
            tend = fs + (end - ts);
        } else {
            tstart = fe - (end - ts);
        }
        out.push(PerlMapSeg::Coord {
            start: tstart,
            end: tend,
        });
        last_to_end = Some(te);
        cur = te + 1;
    }
    match last_to_end {
        None => out.push(PerlMapSeg::Gap { start, end }),
        Some(te) if te < end => out.push(PerlMapSeg::Gap { start: te + 1, end }),
        _ => {}
    }
    if strand == -1 {
        out.reverse();
    }
    out
}

/// `Bio::EnsEMBL::Mapper::map_insert`: an insertion (`start == end + 1`) is mapped
/// as the 2 bp region around it, then placed after the exon base 5' of it or
/// before the exon base 3' of it, whichever flank is exonic.
fn perl_map_insert(ecm: &ExonCoordMapper, start: i64, end: i64, strand: i8) -> Vec<PerlMapSeg> {
    let (start, end) = (end, start);
    let mut coords = perl_map_span(ecm, start, end, strand);
    if coords.len() == 1 {
        if let PerlMapSeg::Coord { start: s, end: e } = coords[0] {
            coords[0] = PerlMapSeg::Coord { start: e, end: s };
        }
        return coords;
    }
    if coords.len() != 2 {
        return Vec::new();
    }
    let (c1, c2) = if strand == -1 {
        (coords[1], coords[0])
    } else {
        (coords[0], coords[1])
    };
    let mut out = Vec::new();
    if let PerlMapSeg::Coord { start: s, end: e } = c1 {
        out.push(if strand == -1 {
            PerlMapSeg::Coord {
                start: s,
                end: e - 1,
            }
        } else {
            PerlMapSeg::Coord {
                start: s + 1,
                end: e,
            }
        });
    }
    if let PerlMapSeg::Coord { start: s, end: e } = c2 {
        let c = if strand == -1 {
            PerlMapSeg::Coord {
                start: s + 1,
                end: e,
            }
        } else {
            PerlMapSeg::Coord {
                start: s,
                end: e - 1,
            }
        };
        if strand == -1 {
            out.insert(0, c);
        } else {
            out.push(c);
        }
    }
    out
}

fn perl_genomic2cdna(ecm: &ExonCoordMapper, start: i64, end: i64, strand: i8) -> Vec<PerlMapSeg> {
    if start == end + 1 {
        perl_map_insert(ecm, start, end, strand)
    } else {
        perl_map_span(ecm, start, end, strand)
    }
}

/// Whether the genome is the reference sequence of `transcript`: a model the
/// cache marks `failed`, whose sequence the alignment's correction did not
/// make the aligned one (AnnotationType/Transcript.pm:579, 619-625). A base
/// its mapper holds and the genome lacks is then no base of the transcript.
pub(crate) fn genome_is_reference(transcript: &Transcript) -> bool {
    transcript.bam_edit_status == Some(BamEditStatus::Failed)
}

/// The two halves `map_insert` gives an insertion between two exonic bases
/// (Mapper.pm:522-542), cut down to the first, the point after the base 5'
/// of it, when the mapper numbers the two bases more than one apart. The
/// halves then span the bases the model holds between them, bases a model
/// the genome is the reference of ([`genome_is_reference`]) does not have.
fn insertion_point_on_the_genome(cdna: &mut Vec<PerlMapSeg>) {
    let spans_model_bases = matches!(
        cdna.as_slice(),
        [PerlMapSeg::Coord { start, .. }, PerlMapSeg::Coord { end, .. }] if start <= end
    );
    if spans_model_bases {
        cdna.truncate(1);
    }
}

/// `TranscriptMapper::genomic2cds`: cDNA segments clipped to the coding region, a
/// UTR overhang becoming a gap on that side.
fn perl_genomic2cds(cdna: &[PerlMapSeg], cstart: i64, cend: i64) -> Vec<PerlMapSeg> {
    let mut out = Vec::new();
    for seg in cdna {
        match *seg {
            PerlMapSeg::Gap { .. } => out.push(*seg),
            PerlMapSeg::Coord { start: s, end: e } => {
                if e < cstart || s > cend {
                    out.push(PerlMapSeg::Gap { start: s, end: e });
                    continue;
                }
                let mut cds_s = s - cstart + 1;
                let mut cds_e = e - cstart + 1;
                if s < cstart {
                    out.push(PerlMapSeg::Gap {
                        start: s,
                        end: cstart - 1,
                    });
                    cds_s = 1;
                }
                let mut end_gap = None;
                if e > cend {
                    end_gap = Some(PerlMapSeg::Gap {
                        start: cend + 1,
                        end: e,
                    });
                    cds_e = cend - cstart + 1;
                }
                out.push(PerlMapSeg::Coord {
                    start: cds_s,
                    end: cds_e,
                });
                if let Some(g) = end_gap {
                    out.push(g);
                }
            }
        }
    }
    out
}

/// Project an allele span the way `BaseTranscriptVariation` does: `cdna_start` /
/// `cdna_end` from the first / last cDNA segment, `cds_start` / `cds_end` from the
/// first / last CDS segment plus the start-exon phase, `translation_start` /
/// `translation_end` from `genomic2pep`. On a model the genome is the
/// reference of, an insertion between two bases the mapper numbers more than
/// one apart is projected at the point after the base 5' of it
/// (`insertion_point_on_the_genome`).
///
/// `None` when the transcript has no coding model or no mapper.
pub(crate) fn perl_span(
    transcript: &Transcript,
    span_start: u64,
    span_end: u64,
) -> Option<PerlSpan> {
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let ecm = &mapper.exon_coord_mapper;
    if ecm.pairs.is_empty() || mapper.cdna_coding_start == 0 || mapper.cdna_coding_end == 0 {
        return None;
    }
    transcript.translation.as_ref()?;
    let strand = transcript.strand.as_i8();
    let (start, end) = (span_start as i64, span_end as i64);
    let mut cdna = perl_genomic2cdna(ecm, start, end, strand);
    if start == end + 1 && genome_is_reference(transcript) {
        insertion_point_on_the_genome(&mut cdna);
    }
    let cstart = mapper.cdna_coding_start as i64;
    let cend = mapper.cdna_coding_end as i64;
    let cds = perl_genomic2cds(&cdna, cstart, cend);
    let phase = i64::from(mapper.start_phase.max(0));
    let first_coord = |v: &[PerlMapSeg]| match v.first() {
        Some(PerlMapSeg::Coord { start, .. }) => Some(*start),
        _ => None,
    };
    let last_coord = |v: &[PerlMapSeg]| match v.last() {
        Some(PerlMapSeg::Coord { end, .. }) => Some(*end),
        _ => None,
    };
    let pep = |c: i64| (c + phase + 2).div_euclid(3);
    Some(PerlSpan {
        cdna_start: first_coord(&cdna),
        cdna_end: last_coord(&cdna),
        cds_start: first_coord(&cds).map(|c| c + phase),
        cds_end: last_coord(&cds).map(|c| c + phase),
        tl_start: first_coord(&cds).map(pep),
        tl_end: last_coord(&cds).map(pep),
        cds_coords: cds,
    })
}

/// `Bio::EnsEMBL::Variation::Utils::VariationEffect::overlap`.
fn perl_overlap(f1_start: i64, f1_end: i64, f2_start: i64, f2_end: i64) -> bool {
    f1_end >= f2_start && f1_start <= f2_end
}

/// Whether the span `v_start..=v_end` overlaps a genomic base of the stop
/// codon: a base the mapper pairs `pairs` hold for one of the last
/// `positions` cDNA positions of the coding sequence
/// `cdna_coding_start..=cdna_coding_end`, three of a whole stop codon and the
/// one or two of a stop codon that the coding sequence's last bases begin.
/// Ensembl VEP reads the last three genomic positions of the coding region
/// (VariationEffect.pm:1418-1422). On a stop codon an intron splits those are
/// one or two intron bases beside the codon's last, and miss the stop bases
/// that end the exon before; on a model whose stop codon an `_rna_edit`
/// completes they are the stop bases the genome holds and the one or two
/// sense bases before them. Each run of adjacent stop bases is one range of
/// Ensembl VEP's `overlap` (:1424-1427), so the range an insertion is given,
/// whose start lies past its end, overlaps the run that holds both ends, as
/// it overlaps a whole codon in Ensembl VEP.
fn overlaps_stop_codon_bases(
    pairs: &[MapperPair],
    (cdna_coding_start, cdna_coding_end): (i64, i64),
    positions: i64,
    (v_start, v_end): (i64, i64),
) -> bool {
    let mut bases = [0i64; 3];
    let mut held = 0;
    let first = cdna_coding_end - (positions.clamp(1, 3) - 1);
    for cdna in cdna_coding_start.max(first)..=cdna_coding_end {
        let Some(p) = pairs
            .iter()
            .find(|p| p.from_start as i64 <= cdna && cdna <= p.from_end as i64)
        else {
            continue;
        };
        let offset = cdna - p.from_start as i64;
        bases[held] = if p.ori >= 0 {
            p.to_start as i64 + offset
        } else {
            p.to_end as i64 - offset
        };
        held += 1;
    }
    let bases = &mut bases[..held];
    bases.sort_unstable();
    bases
        .chunk_by(|a, b| b - a == 1)
        .any(|run| perl_overlap(v_start, v_end, run[0], run[run.len() - 1]))
}

/// Whether `edited`, the coding sequence with an allele written in at CDS
/// position `cds_start`, reads the reference's residue at every codon of `cds`
/// before the old stop position, the last three bases of `cds`. A codon whose
/// bases differ keeps its residue where both codons translate to one residue.
/// `edited` is at least as long as `cds`, and `cds` holds a codon.
fn residues_before_stop_kept(cds: &[u8], edited: &[u8], cds_start: i64, table: u8) -> bool {
    let stop_idx = cds.len() - 3;
    let first = ((cds_start - 1).max(0) as usize / 3 * 3).min(stop_idx);
    let whole = first + (stop_idx - first) / 3 * 3;
    let same_residue = |i: usize| {
        let (reference, read) = (&cds[i..i + 3], &edited[i..i + 3]);
        reference.eq_ignore_ascii_case(read) || {
            let residue = perl_translate_codon(read, table);
            residue != b'X' && residue == perl_translate_codon(reference, table)
        }
    };
    (first..whole).step_by(3).all(same_residue)
        && cds[whole..stop_idx].eq_ignore_ascii_case(&edited[whole..stop_idx])
}

/// What the sequence a record gives shows at the old stop position, for a
/// length-changing span that holds the whole stop codon and has no coordinates
/// for the stop codon edit test ([`PerlCodingEval::removed_stop_codon`]).
#[derive(Clone, Copy, PartialEq, Eq)]
enum GivenStopCodon {
    /// No stop codon is shown there, or the given sequence ends before it.
    Lost,
    /// A stop codon is shown there behind the reference's residues, or the
    /// edit can be written where it takes no coding base and no stop base.
    Retained,
    /// A stop codon is shown there and a residue before it differs.
    Changed,
}

/// How many positions a deletion is followed along the reference sequence in
/// each direction for the other places the same edit can be written at.
const EQUIVALENT_SITES_LIMIT: i64 = 10_000;

/// The exons of a transcript on an axis that runs in its own direction: the
/// genomic coordinate on the forward strand, its negative on the reverse.
struct ExonAxis {
    /// `(first position, last position, first cDNA position)` of each exon.
    exons: Vec<(i64, i64, i64)>,
    cdna_coding_start: i64,
    cdna_coding_end: i64,
}

impl ExonAxis {
    /// `None` unless the mapper pairs are the transcript's exons in cDNA order
    /// from position 1, each on `strand` and as long on the genome as in the
    /// cDNA.
    fn new(
        pairs: &[MapperPair],
        strand: i8,
        (cdna_coding_start, cdna_coding_end): (i64, i64),
    ) -> Option<Self> {
        let mut sorted: Vec<&MapperPair> = pairs.iter().collect();
        sorted.sort_by_key(|p| p.from_start);
        let mut exons: Vec<(i64, i64, i64)> = Vec::with_capacity(sorted.len());
        let mut next_cdna = 1;
        for p in sorted {
            let (from_start, from_end) = (p.from_start as i64, p.from_end as i64);
            let (to_start, to_end) = (p.to_start as i64, p.to_end as i64);
            if p.ori != strand
                || from_start != next_cdna
                || from_end - from_start != to_end - to_start
            {
                return None;
            }
            let (first, last) = if strand == 1 {
                (to_start, to_end)
            } else {
                (-to_end, -to_start)
            };
            if exons.last().is_some_and(|&(_, before, _)| before >= first) {
                return None;
            }
            exons.push((first, last, from_start));
            next_cdna = from_end + 1;
        }
        (!exons.is_empty()).then_some(Self {
            exons,
            cdna_coding_start,
            cdna_coding_end,
        })
    }

    /// The cDNA position of the axis position `pos`, `None` off the exons.
    fn cdna(&self, pos: i64) -> Option<i64> {
        self.exons
            .iter()
            .find(|&&(first, last, _)| first <= pos && pos <= last)
            .map(|&(first, _, cdna)| cdna + pos - first)
    }

    /// The axis position of the cDNA position `cdna`.
    fn position(&self, cdna: i64) -> Option<i64> {
        self.exons
            .iter()
            .find(|&&(first, last, start)| start <= cdna && cdna <= start + last - first)
            .map(|&(first, _, start)| first + cdna - start)
    }

    /// The number of exon bases before the axis position `pos`.
    fn exonic_before(&self, pos: i64) -> usize {
        self.exons
            .iter()
            .map(|&(first, last, _)| (last.min(pos - 1) - first + 1).max(0))
            .sum::<i64>() as usize
    }

    /// Of the axis positions `lo..=hi`: the coding bases, the stop codon's
    /// among them, and the exon bases.
    fn taken(&self, lo: i64, hi: i64) -> (i64, i64) {
        let (mut coding, mut exonic) = (0, 0);
        for &(first, last, cdna) in &self.exons {
            let (x, y) = (first.max(lo), last.min(hi));
            if y < x {
                continue;
            }
            exonic += y - x + 1;
            let (c1, c2) = (cdna + x - first, cdna + y - first);
            coding += (c2.min(self.cdna_coding_end) - c1.max(self.cdna_coding_start) + 1).max(0);
        }
        (coding, exonic)
    }
}

/// Whether two stretches of coding sequence of one length, each beginning
/// `filler` bases into a codon, read the same residue at every codon. A codon
/// whose bases differ keeps its residue where both codons translate to one
/// residue; a codon cut short at either end is compared base for base.
fn same_residues(a: &[u8], b: &[u8], filler: usize, table: u8) -> bool {
    let head = ((3 - filler % 3) % 3).min(a.len());
    a[..head].eq_ignore_ascii_case(&b[..head])
        && a[head..].chunks(3).zip(b[head..].chunks(3)).all(|(x, y)| {
            x.eq_ignore_ascii_case(y)
                || (x.len() == 3 && {
                    let residue = perl_translate_codon(x, table);
                    residue != b'X' && residue == perl_translate_codon(y, table)
                })
        })
}

/// Perl-truthy integer: defined and non-zero.
fn truthy(v: Option<i64>) -> Option<i64> {
    v.filter(|&x| x != 0)
}

/// Perl's `substr($s, $off, $len)`: `None` past the end of the string, clipped otherwise.
fn perl_substr(s: &[u8], off: i64, len: i64) -> Option<&[u8]> {
    if off < 0 || off as usize > s.len() {
        return None;
    }
    let off = off as usize;
    let stop = if len < 0 {
        (s.len() as i64 + len).max(off as i64) as usize
    } else {
        (off + len as usize).min(s.len())
    };
    Some(&s[off..stop])
}

/// Perl's 4-arg `substr($s, $off, $len) = $repl`, clipped like Perl does.
fn perl_substr_assign(s: &mut Vec<u8>, off: i64, len: i64, repl: &[u8]) {
    let off = (off.max(0) as usize).min(s.len());
    let stop = (off + len.max(0) as usize).min(s.len());
    s.splice(off..stop, repl.iter().copied());
}

fn is_perl_unambiguous_dna(seq: &[u8]) -> bool {
    !seq.is_empty()
        && seq
            .iter()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'-'))
}

fn iupac_bases(b: u8) -> &'static [u8] {
    match b.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' | b'X' => b"ACGT",
        _ => b"",
    }
}

/// `Bio::Tools::CodonTable::translate` for one codon: an ambiguous codon resolves
/// when every expansion agrees, `B` for Asp/Asn and `Z` for Glu/Gln, else `X`.
fn perl_translate_codon(codon: &[u8], table: u8) -> u8 {
    if codon.len() != 3 {
        return b'X';
    }
    if codon
        .iter()
        .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        return translate_codon_with_table(codon, table);
    }
    let (o0, o1, o2) = (
        iupac_bases(codon[0]),
        iupac_bases(codon[1]),
        iupac_bases(codon[2]),
    );
    if o0.is_empty() || o1.is_empty() || o2.is_empty() {
        return b'X';
    }
    let mut aas: Vec<u8> = Vec::new();
    for &a in o0 {
        for &b in o1 {
            for &c in o2 {
                let aa = translate_codon_with_table(&[a, b, c], table);
                if !aas.contains(&aa) {
                    aas.push(aa);
                }
            }
        }
    }
    match aas.len() {
        1 => aas[0],
        2 if aas.contains(&b'D') && aas.contains(&b'N') => b'B',
        2 if aas.contains(&b'E') && aas.contains(&b'Q') => b'Z',
        _ => b'X',
    }
}

fn perl_translate(seq: &[u8], table: u8) -> Vec<u8> {
    seq.as_chunks::<3>()
        .0
        .iter()
        .map(|c| perl_translate_codon(c, table))
        .collect()
}

/// `TranscriptVariationAllele::peptide` applies the translation's SeqEdits to the
/// reference peptide only, so a selenocysteine `TGA` reads `U`, not `*`. A cache
/// without `seq_edits` still has its `peptide`, which is `Transcript::translate`
/// with the edits applied, so an edited residue is the cached residue wherever
/// it differs from the codon's own translation. Two differences in that string
/// have no SeqEdit behind them and are skipped: position 1, which `translate`
/// rewrites to `M` for any start codon, and a partial terminal codon, which
/// BioPerl resolves when its expansions agree.
///
/// `pep` is the reference codon window's translation for protein positions
/// `tl_start..`; `cds_len` bounds the complete codons.
pub(crate) fn overlay_seq_edits(pep: &mut [u8], pep_seq: &[u8], tl_start: i64, cds_len: usize) {
    for (i, aa) in pep.iter_mut().enumerate() {
        let pos = tl_start + i as i64;
        if pos < 2 || pos as usize * 3 > cds_len {
            continue;
        }
        if let Some(&edited) = pep_seq.get(pos as usize - 1) {
            if edited != *aa {
                *aa = edited;
            }
        }
    }
}

/// The residue an `initial_met` SeqEdit writes over protein position 1 of the
/// reference peptide, read from the cached peptide of a cache without
/// `seq_edits`: `Some(b'M')` when that peptide starts with `M` although the
/// first codon translates otherwise and is not a start codon of the table.
/// `Transcript::translate` writes `M` over any start codon of the table, so
/// there the edit cannot be told from that rewrite and the codon's own residue
/// stands: `None`.
pub(crate) fn initial_met_edit(cds: &[u8], pep_seq: &[u8], table: u8) -> Option<u8> {
    let first = cds.get(..3)?;
    if pep_seq.first() != Some(&b'M')
        || perl_translate_codon(first, table) == b'M'
        || is_start_codon(first, table)
    {
        return None;
    }
    Some(b'M')
}

/// The SeqEdit loop of `TranscriptVariationAllele::peptide` over `pep`, the
/// reference codon window's translation for protein positions
/// `tv_start..=tv_end`: every edit overlapping the span writes its residues over
/// the positions it covers. A position past the end of an edit's replacement
/// (an edit that shortens the protein) makes the peptide undefined, `false`.
pub(crate) fn apply_seq_edits(
    pep: &mut Vec<u8>,
    edits: &[SeqEdit],
    tv_start: i64,
    tv_end: i64,
) -> bool {
    for edit in edits {
        let (se_start, se_end) = (edit.start as i64, edit.end as i64);
        if !perl_overlap(tv_start, tv_end, se_start, se_end) {
            continue;
        }
        let alt = edit.alt_seq.as_bytes();
        for tv_pos in tv_start.max(se_start)..=tv_end.min(se_end) {
            let from = (tv_pos - se_start) as usize;
            if from >= alt.len() {
                return false;
            }
            let at = (tv_pos - tv_start) as usize;
            if at < pep.len() {
                pep[at] = alt[from];
            } else if at == pep.len() {
                pep.push(alt[from]);
            }
        }
    }
    true
}

/// Whether the SeqEdit loop of `TranscriptVariationAllele::peptide` leaves the
/// reference peptide of protein positions `tv_start..=tv_end` undefined, as
/// [`apply_seq_edits`] reports it: an edit over the span covers a position
/// past the end of its replacement (TranscriptVariationAllele.pm:763-765).
pub(crate) fn seq_edits_leave_reference_peptide_undefined(
    edits: &[SeqEdit],
    tv_start: i64,
    tv_end: i64,
) -> bool {
    edits.iter().any(|edit| {
        let (se_start, se_end) = (edit.start as i64, edit.end as i64);
        perl_overlap(tv_start, tv_end, se_start, se_end)
            && (tv_start.max(se_start)..=tv_end.min(se_end))
                .any(|tv_pos| (tv_pos - se_start) as usize >= edit.alt_seq.len())
    })
}

/// The reference-only SeqEdit step of `TranscriptVariationAllele::peptide`, over
/// `pep` as [`apply_seq_edits`] takes it. With the transcript's `seq_edits` that
/// is the loop itself; a cache without them reads the edits from its cached
/// peptide `pep_seq` instead, [`overlay_seq_edits`] past position 1 and
/// [`initial_met_edit`] at it, which can never leave the peptide undefined.
pub(crate) fn edit_reference_peptide(
    pep: &mut Vec<u8>,
    seq_edits: &[SeqEdit],
    pep_seq: &[u8],
    cds: &[u8],
    table: u8,
    tv_start: i64,
    tv_end: i64,
) -> bool {
    if !seq_edits.is_empty() {
        return apply_seq_edits(pep, seq_edits, tv_start, tv_end);
    }
    if !pep.is_empty() {
        overlay_seq_edits(pep, pep_seq, tv_start, cds.len());
        if tv_start == 1 {
            if let Some(m) = initial_met_edit(cds, pep_seq, table) {
                pep[0] = m;
            }
        }
    }
    true
}

/// A cached peptide as `Bio::EnsEMBL::Transcript::translate` writes it from the
/// release 114 core API on: position 1 is `M` whenever the first codon is a start
/// codon of the table, whatever residue that codon translates to on its own. A
/// cache dumped by an earlier API carries the codon's own residue there, and
/// Ensembl VEP's `_peptide` (the string `hgvs_protein` reads for surrounding, frameshift
/// and deletion peptides) is whatever the cache holds; this rewrite reads every
/// cache as the release 114 API and its successors write it. The per-codon translation of
/// [`PerlCodingEval::peptide`] is untouched, as Ensembl VEP's `peptide` is.
fn with_initial_met<'a>(pep: &'a [u8], cds: &[u8], table: u8) -> Cow<'a, [u8]> {
    match (pep.first(), cds.get(..3)) {
        (Some(&first), Some(codon)) if first != b'M' && is_start_codon(codon, table) => {
            let mut owned = pep.to_vec();
            owned[0] = b'M';
            Cow::Owned(owned)
        }
        _ => Cow::Borrowed(pep),
    }
}

/// `Bio::EnsEMBL::Transcript::translate` on the translateable sequence: whole
/// codons only, the terminal stop dropped, a start codon read as `M`. This is
/// the `_peptide` Ensembl VEP compares against in `ref_eq_alt_sequence` when the cache
/// carries no peptide.
fn perl_transcript_peptide(cds: &[u8], table: u8) -> Vec<u8> {
    let whole = &cds[..cds.len() - cds.len() % 3];
    if whole.is_empty() {
        return Vec::new();
    }
    let mut pep = perl_translate(whole, table);
    if is_ter_codon(&whole[whole.len() - 3..], table) {
        pep.pop();
    }
    if pep.first().is_some_and(|&a| a != b'M') && is_start_codon(&whole[..3], table) {
        pep[0] = b'M';
    }
    pep
}

/// The transcript's peptide as `ref_eq_alt_sequence` and the protein notation
/// read it: the one the cache carries, through [`with_initial_met`], or with
/// none cached the translation of the coding sequence
/// ([`perl_transcript_peptide`]).
fn transcript_peptide<'a>(cached: Option<&'a str>, cds: &[u8], table: u8) -> Cow<'a, [u8]> {
    match cached {
        Some(peptide) => with_initial_met(peptide.as_bytes(), cds, table),
        None => Cow::Owned(perl_transcript_peptide(cds, table)),
    }
}

/// What stands at the coding end of a transcript, read from the coding
/// sequence the cache records, in its frame from its first base with the
/// padding the start exon's phase gives it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum StopCodonAtTheCodingEnd {
    /// No stop codon: the transcript carries `cds_end_NF`, its last whole
    /// codon names a residue, or its coding sequence ends inside a codon that
    /// a 3' UTR follows, that adenines do not complete to a stop codon, or
    /// that its record does not show completed.
    None,
    /// The coding sequence is a whole number of codons and its last codon is
    /// a stop codon of the transcript's codon table.
    Whole,
    /// The coding sequence ends `bases` bases, one or two, into a codon on
    /// the transcript's last base, those bases followed by adenines are a
    /// stop codon of the table, and the record shows the model complete
    /// ([`record_shows_the_model_complete`]): polyadenylation completes the
    /// stop codon, and the row is read on the completed coding sequence
    /// ([`completed_by_polyadenylation`]).
    BegunByTheTail { bases: u8 },
}

/// The codon tables under which a reference sequence record annotates a
/// coding sequence that ends inside its stop codon as completed by
/// polyadenylation: the vertebrate, invertebrate, echinoderm and ascidian
/// mitochondrial codes.
const CODON_TABLES_THAT_SHOW_A_COMPLETED_STOP_CODON: [u8; 4] = [2, 5, 9, 13];

/// Whether the record of `transcript` carries a consensus coding sequence
/// id, which names a complete coding sequence.
fn carries_a_consensus_coding_sequence_id(transcript: &Transcript) -> bool {
    transcript
        .ccds
        .as_deref()
        .is_some_and(|id| !id.is_empty() && id != "-")
}

/// Whether the record of `transcript` shows its model complete where the
/// coding sequence ends inside a stop codon on the transcript's last base:
/// the transcript is translated under one of
/// [`CODON_TABLES_THAT_SHOW_A_COMPLETED_STOP_CODON`], or its record carries a
/// consensus coding sequence id.
fn record_shows_the_model_complete(transcript: &Transcript) -> bool {
    CODON_TABLES_THAT_SHOW_A_COMPLETED_STOP_CODON.contains(&codon_table_for(transcript))
        || carries_a_consensus_coding_sequence_id(transcript)
}

/// The stop codon that stands at the coding end of `transcript`
/// ([`StopCodonAtTheCodingEnd`]). A coding sequence with another stop codon
/// inside is closed by its last codon all the same. Ensembl VEP reads the last
/// three coding positions as the stop codon on every transcript without
/// `cds_end_NF`, whatever they hold (VariationEffect.pm:1373-1381, 1400,
/// 1418-1427).
pub(crate) fn stop_codon_at_the_coding_end(transcript: &Transcript) -> StopCodonAtTheCodingEnd {
    if transcript.flags.iter().any(|flag| &**flag == "cds_end_NF") {
        return StopCodonAtTheCodingEnd::None;
    }
    let Some(vefc) = transcript.vefc.as_ref() else {
        return StopCodonAtTheCodingEnd::None;
    };
    let Some(cds) = vefc.translateable_seq.as_deref().map(str::as_bytes) else {
        return StopCodonAtTheCodingEnd::None;
    };
    let table = codon_table_for(transcript);
    let tail = cds.len() % 3;
    if tail == 0 {
        return match cds.len().checked_sub(3) {
            Some(last_codon) if is_ter_codon(&cds[last_codon..], table) => {
                StopCodonAtTheCodingEnd::Whole
            }
            _ => StopCodonAtTheCodingEnd::None,
        };
    }
    // The transcript has no 3' UTR where its last cDNA position is its last
    // coding position.
    let ends_at_the_coding_end = vefc.mapper.as_ref().is_some_and(|mapper| {
        let cdna_length = mapper
            .exon_coord_mapper
            .pairs
            .iter()
            .map(|pair| pair.from_end)
            .max();
        mapper.cdna_coding_end > 0 && cdna_length == Some(mapper.cdna_coding_end)
    });
    let mut completed = [b'A'; 3];
    completed[..tail].copy_from_slice(&cds[cds.len() - tail..]);
    if ends_at_the_coding_end
        && is_ter_codon(&completed, table)
        && record_shows_the_model_complete(transcript)
    {
        StopCodonAtTheCodingEnd::BegunByTheTail { bases: tail as u8 }
    } else {
        StopCodonAtTheCodingEnd::None
    }
}

/// `transcript` with the stop codon that polyadenylation completes carried
/// in its record, for a transcript whose record shows the model complete
/// ([`StopCodonAtTheCodingEnd::BegunByTheTail`]); `None` for every other
/// transcript. The copy's coding sequence is the cached one
/// followed by the one or two adenines that complete its last codon, and its
/// coding end lies that many cDNA positions later. It is the shape of a
/// model whose `_rna_edit` appends the same bases, and every term and cell
/// of a row on `transcript` is read on it
/// (`crate::consequences::calculate_consequences`).
///
/// The poly(A) tail that completes the stop codon goes on past it, so the
/// copy's 3' UTR sequence is adenines: where an allele removes bases at or
/// before the stop codon, the edit tests, the codon windows and the protein
/// notation read on into them, and find the residues the message goes on
/// with and no further stop codon. As many adenines follow the stop codon as
/// the completed coding sequence has bases: an allele inside the transcript
/// removes at most every coding base, so the sequence it leaves still holds
/// adenines up to the end of the stop codon's old place, and no reader meets
/// the end of the sequence where the message goes on.
///
/// No mapper pair holds the completing adenines or the tail, which have no
/// genomic base: no allele lies on them, no row takes a 3' UTR term from
/// them, and the transcript's exons and cDNA length stay the record's.
pub(crate) fn completed_by_polyadenylation(transcript: &Transcript) -> Option<Transcript> {
    let StopCodonAtTheCodingEnd::BegunByTheTail { bases } =
        stop_codon_at_the_coding_end(transcript)
    else {
        return None;
    };
    let adenines = 3 - u64::from(bases);
    let mut completed = transcript.clone();
    let vefc = completed.vefc.as_mut()?;
    let coding = vefc.translateable_seq.as_mut()?;
    for _ in 0..adenines {
        coding.push('A');
    }
    vefc.three_prime_utr = Some("A".repeat(coding.len()));
    vefc.mapper.as_mut()?.cdna_coding_end += adenines;
    completed.cdna_coding_end = completed.cdna_coding_end.map(|end| end + adenines);
    Some(completed)
}

/// Whether a start codon stands at the coding start of `transcript`: it
/// carries no `cds_start_NF`, its first codon is whole, and either its
/// translation opens with Met or its record carries a consensus coding
/// sequence id, which names a complete coding sequence. The translation is
/// the one the cache records, which opens with Met where an edit of residue 1
/// asserts a start on another triplet, or with none cached the coding
/// sequence's own, to which `Transcript::translate` gives Met for a first
/// codon that is a start codon of the table. Ensembl VEP reads the first
/// three coding positions as the start codon on every transcript without the
/// flag, whatever they hold (VariationEffect.pm:975-986).
pub(crate) fn start_codon_stands_at_the_coding_start(transcript: &Transcript) -> bool {
    let Some(vefc) = transcript.vefc.as_ref() else {
        return false;
    };
    let Some(cds) = vefc.translateable_seq.as_deref() else {
        return false;
    };
    let peptide = transcript_peptide(
        vefc.peptide.as_deref(),
        cds.as_bytes(),
        codon_table_for(transcript),
    );
    start_codon_opens(transcript, &peptide)
}

/// The test of [`start_codon_stands_at_the_coding_start`] on `peptide`, the
/// transcript's peptide as [`transcript_peptide`] gives it. A start phase
/// above zero pads the first codon, which is then no codon of the transcript.
fn start_codon_opens(transcript: &Transcript, peptide: &[u8]) -> bool {
    let padded = transcript
        .vefc
        .as_ref()
        .and_then(|vefc| vefc.mapper.as_ref())
        .is_some_and(|mapper| mapper.start_phase > 0);
    !padded
        && (peptide.first() == Some(&b'M') || carries_a_consensus_coding_sequence_id(transcript))
        && !transcript
            .flags
            .iter()
            .any(|flag| &**flag == "cds_start_NF")
}

/// What stands at the two coding ends of `transcript`, each in words, for a
/// reader of cache records that is held to this reading; `None` for a
/// transcript whose record holds no coding sequence.
///
/// The first text is the stop codon at the coding end: `whole` where the
/// coding sequence is a whole number of codons closed by a stop codon of the
/// transcript's codon table; `begun by the tail, 1 base, shown` where it ends
/// one base into a stop codon that adenines complete on the transcript's last
/// base (`2 bases` for two) and the record shows the model complete, by a
/// vertebrate, invertebrate, echinoderm or ascidian mitochondrial codon table
/// or by a consensus coding sequence id; `none shown, flagged` on a transcript
/// that carries `cds_end_NF`; and `none shown` on every other transcript.
///
/// The second text is the start codon at the coding start: `shown` where the
/// transcript carries no `cds_start_NF`, its first codon is whole, and its
/// translation opens with Met or its record carries a consensus coding
/// sequence id; elsewhere `not shown`, followed by `, flagged` on a transcript
/// that carries `cds_start_NF` and by `, padded` on one whose start phase is
/// above zero.
#[doc(hidden)]
pub fn coding_ends_in_words(transcript: &Transcript) -> Option<(String, String)> {
    let coding = transcript.vefc.as_ref()?.translateable_seq.as_deref()?;
    if coding.is_empty() {
        return None;
    }
    let carries = |code: &str| transcript.flags.iter().any(|flag| &**flag == code);
    let stop = match stop_codon_at_the_coding_end(transcript) {
        StopCodonAtTheCodingEnd::Whole => "whole".to_string(),
        StopCodonAtTheCodingEnd::BegunByTheTail { bases } => format!(
            "begun by the tail, {bases} {}, shown",
            if bases == 1 { "base" } else { "bases" }
        ),
        StopCodonAtTheCodingEnd::None if carries("cds_end_NF") => "none shown, flagged".to_string(),
        StopCodonAtTheCodingEnd::None => "none shown".to_string(),
    };
    let mut start = "shown".to_string();
    if !start_codon_stands_at_the_coding_start(transcript) {
        start = "not shown".to_string();
        if carries("cds_start_NF") {
            start.push_str(", flagged");
        }
        let padded = transcript
            .vefc
            .as_ref()
            .and_then(|vefc| vefc.mapper.as_ref())
            .is_some_and(|mapper| mapper.start_phase > 0);
        if padded {
            start.push_str(", padded");
        }
    }
    Some((stop, start))
}

/// Every SO term `perl_coding_terms` can emit; the caller replaces exactly this set.
pub(crate) fn is_perl_coding_term(c: Consequence) -> bool {
    matches!(
        c,
        Consequence::StopGained
            | Consequence::FrameshiftVariant
            | Consequence::StopLost
            | Consequence::StartLost
            | Consequence::InframeInsertion
            | Consequence::InframeDeletion
            | Consequence::MissenseVariant
            | Consequence::ProteinAlteringVariant
            | Consequence::IncompleteTerminalCodonVariant
            | Consequence::StartRetainedVariant
            | Consequence::StopRetainedVariant
            | Consequence::SynonymousVariant
            | Consequence::CodingSequenceVariant
    )
}

/// What the stop codon edit test of [`PerlCodingEval::edited_stop_codon`] finds.
enum EditedStopCodon {
    /// The test is not made: the variant is neither an insertion nor a
    /// deletion, its span lacks a cDNA coordinate at one end or a CDS start, or
    /// the CDS is shorter than a codon.
    NotMade,
    /// The edited sequence is shorter than the CDS (VariationEffect.pm:1524).
    Shorter,
    /// The three bases the edited sequence holds at the old stop position,
    /// and whether every codon before them reads the reference's residue.
    Read { codon: [u8; 3], residues_kept: bool },
}

/// One `TranscriptVariationAllele` evaluation with Ensembl VEP's `_predicate_cache`.
struct PerlCodingEval<'a> {
    transcript: &'a Transcript,
    cds: &'a [u8],
    pep_seq: Cow<'a, [u8]>,
    seq_edits: &'a [SeqEdit],
    utr5: Option<Vec<u8>>,
    utr3: Option<Vec<u8>>,
    table: u8,
    span: PerlSpan,
    /// The transcript's strand, the direction an inserted allele's length
    /// moves the variant start in (`_overlaps_stop_codon_cil`).
    strand: i8,
    /// Whether the alleles were given on the strand opposite the transcript's
    /// and so read reverse-complemented (`feature_seq`).
    opposite_strand: bool,
    /// The reference row of `--process_ref_homs`: Ensembl VEP's reference allele
    /// object itself, whose `is_reference` reads its codon from the translateable
    /// sequence and its peptide with the sequence edits applied
    /// (TranscriptVariationAllele.pm:742-770, 868, 877), so both alleles of the
    /// evaluation take the reference's reading whatever the written bases.
    reference_row: bool,
    /// Whether a start codon stands at the coding start
    /// ([`start_codon_stands_at_the_coding_start`]).
    start_codon_stands: bool,
    /// The stop codon at the coding end ([`stop_codon_at_the_coding_end`]).
    stop_codon: StopCodonAtTheCodingEnd,
    cdna_coding_start: i64,
    cdna_coding_end: i64,
    coding_region: Option<(i64, i64)>,
    frameshift_introns: Vec<(i64, i64)>,
    vf_start: i64,
    vf_end: i64,
    ref_allele: &'a [u8],
    alt_allele: &'a [u8],
    feature_seq: Vec<u8>,
    /// The variant's chromosome and the reference sequence, read for the other
    /// places a deletion can be written at ([`Self::equivalent_sites`]).
    chr: &'a str,
    fasta: Option<&'a vep_fasta::IndexedFasta>,
    seq_length: Option<i64>,
    snp: bool,
    insertion: bool,
    deletion: bool,
    // _predicate_cache
    codon_memo: [Option<Option<Vec<u8>>>; 2],
    peptide_memo: [Option<Option<Vec<u8>>>; 2],
    peps_memo: Option<Option<(Vec<u8>, Vec<u8>)>>,
    overlaps_stop_memo: Option<bool>,
    overlaps_stop_cil_memo: Option<bool>,
    overlaps_start_memo: Option<bool>,
    ins_del_stop_altered_memo: Option<bool>,
    edited_stop_codon_kept_cil_memo: Option<bool>,
    removed_stop_codon_memo: Option<Option<GivenStopCodon>>,
    ins_del_start_altered_memo: Option<bool>,
    inv_start_altered_memo: Option<bool>,
    snp_start_altered_memo: Option<bool>,
    start_lost_memo: Option<bool>,
    stop_lost_memo: Option<bool>,
    stop_retained_memo: Option<bool>,
    stop_gained_memo: Option<bool>,
    /// Whether the alternate peptide of an allele that holds unknown bases is
    /// read for a gained stop codon
    /// ([`Self::gains_a_stop_codon_before_unknown_bases`]).
    reads_unknown_bases: bool,
    /// Whether that peptide is read on a row whose stop codon the stop codon
    /// test finds kept ([`Self::keeps_a_stop_codon_before_unknown_bases`]).
    reads_a_kept_stop_codon: bool,
}

impl<'a> PerlCodingEval<'a> {
    fn allele(&self, is_ref: bool) -> &'a [u8] {
        if is_ref {
            self.ref_allele
        } else {
            self.alt_allele
        }
    }

    /// The bases of the alt allele as the transcript reads them: the first `-`
    /// dropped, the rest reverse-complemented when the alleles are given on the
    /// other strand.
    fn alt_allele_on_transcript(&self) -> Vec<u8> {
        let mut a: Vec<u8> = Vec::with_capacity(self.alt_allele.len());
        let mut dash_removed = false;
        for &b in self.alt_allele {
            if b == b'-' && !dash_removed {
                dash_removed = true;
            } else {
                a.push(b);
            }
        }
        if !a.is_empty() && self.opposite_strand {
            a = reverse_complement(&a);
        }
        a
    }

    /// `_get_alternate_cds` for the alt allele: upstream + allele + downstream,
    /// the incomplete-codon trim (a no-op past 2 bp: `return $seq if
    /// $full_length = $keep_length` assigns), then the 3' UTR.
    fn alternate_cds(&self) -> Option<Vec<u8>> {
        let cds_start = self.span.cds_start?;
        let cds_end = self.span.cds_end?;
        let up_end = ((cds_start - 1).max(0) as usize).min(self.cds.len());
        let down_start = (cds_end.max(0) as usize).min(self.cds.len());
        let a = self.alt_allele_on_transcript();
        let mut alt = Vec::with_capacity(self.cds.len() + a.len() + 8);
        alt.extend_from_slice(&self.cds[..up_end]);
        alt.extend_from_slice(&a);
        alt.extend_from_slice(&self.cds[down_start..]);
        if alt.len() < 3 {
            alt.clear();
        }
        if let Some(utr) = self.utr3.as_ref() {
            alt.extend_from_slice(utr);
        }
        Some(alt)
    }

    /// `TranscriptVariationAllele::codon`. `Some(b"-")` is Ensembl VEP's `'-'`.
    fn codon(&mut self, is_ref: bool) -> Option<Vec<u8>> {
        let slot = usize::from(!is_ref);
        if let Some(m) = &self.codon_memo[slot] {
            return m.clone();
        }
        self.codon_memo[slot] = Some(None);
        let (Some(tl_start), Some(tl_end)) = (truthy(self.span.tl_start), truthy(self.span.tl_end))
        else {
            return None;
        };
        if !seq_is_dna(self.allele(is_ref)) {
            return None;
        }
        let codon_cds_start = tl_start * 3 - 2;
        let codon_cds_end = tl_end * 3;
        let codon_len = codon_cds_end - codon_cds_start + 1;
        let (Some(cds_start), Some(cds_end)) = (self.span.cds_start, self.span.cds_end) else {
            return None;
        };
        let vf_nt_len = cds_end - cds_start + 1;
        // Ensembl VEP builds `_get_alternate_cds` for both alleles but reads the reference
        // codon from the translateable sequence (TranscriptVariationAllele.pm:866-868), over
        // whose coding span it first writes the reference allele on a model that carries an
        // `_rna_edit` attribute (:872-874). No such write is made here.
        let codon: Option<Vec<u8>> = if is_ref || self.reference_row {
            perl_substr(self.cds, codon_cds_start - 1, codon_len).map(<[u8]>::to_vec)
        } else {
            let alt_cds = self.alternate_cds()?;
            let allele_len = self.seq_length?;
            perl_substr(
                &alt_cds,
                codon_cds_start - 1,
                codon_len + (allele_len - vf_nt_len),
            )
            .map(<[u8]>::to_vec)
        };
        let value = match codon {
            Some(c) if !c.is_empty() => c,
            _ => {
                self.peptide_memo[slot] = Some(Some(b"-".to_vec()));
                b"-".to_vec()
            }
        };
        self.codon_memo[slot] = Some(Some(value.clone()));
        Some(value)
    }

    /// `TranscriptVariationAllele::peptide`.
    fn peptide(&mut self, is_ref: bool) -> Option<Vec<u8>> {
        let slot = usize::from(!is_ref);
        if let Some(m) = &self.peptide_memo[slot] {
            return m.clone();
        }
        self.peptide_memo[slot] = Some(None);
        // An alternate allele that holds unknown bases is translated where
        // its known bases settle a stop codon before the first of them, a
        // gained one (`gains_a_stop_codon_before_unknown_bases`) or the
        // transcript's own, kept
        // (`keeps_a_stop_codon_before_unknown_bases`); a codon that holds an
        // unknown base then reads as `perl_translate_codon` reads it.
        if !is_perl_unambiguous_dna(self.allele(is_ref))
            && (is_ref
                || self.reference_row
                || !(self.gains_a_stop_codon_before_unknown_bases()
                    || self.keeps_a_stop_codon_before_unknown_bases()))
        {
            return None;
        }
        let codon = self.codon(is_ref)?;
        if let Some(Some(p)) = &self.peptide_memo[slot] {
            return Some(p.clone());
        }
        let whole_len = codon.len() / 3 * 3;
        let mut pep = perl_translate(&codon[..whole_len], self.table);
        // SeqEdits apply to the reference peptide only, before the partial-codon `X`.
        if is_ref || self.reference_row {
            if let (Some(a), Some(b)) = (self.span.tl_start, self.span.tl_end) {
                if !edit_reference_peptide(
                    &mut pep,
                    self.seq_edits,
                    &self.pep_seq,
                    self.cds,
                    self.table,
                    a.min(b),
                    a.max(b),
                ) {
                    return None;
                }
            }
        }
        if whole_len < codon.len() && pep != b"*" {
            pep.push(b'X');
        }
        if pep.is_empty() {
            pep.push(b'-');
        }
        self.peptide_memo[slot] = Some(Some(pep.clone()));
        Some(pep)
    }

    /// [`stop_codon_before_unknown_bases`] of the alternate allele's codon
    /// string, where the row has one reading: the reference codons it stands
    /// on hold no stop codon, and the reference peptide either opens the
    /// alternate peptide or is held against last codons of known bases
    /// alone. Asked only of an evaluation that reads unknown bases
    /// ([`Self::reading_unknown_bases`]).
    fn gains_a_stop_codon_before_unknown_bases(&mut self) -> bool {
        if !self.reads_unknown_bases {
            return false;
        }
        let table = self.table;
        let Some(codon) = self.codon(false) else {
            return false;
        };
        if !stop_codon_before_unknown_bases(&codon, table) {
            return false;
        }
        let reference = match self.peptide(true) {
            Some(reference) if reference == b"-" => Vec::new(),
            Some(reference) => reference,
            None => return false,
        };
        // The predicates and the protein notation hold the two peptides
        // against each other from both ends. The first residues are settled
        // through the stop. Where the reference peptide is not the alternate
        // peptide's opening, `protein_altering_variant` asks whether the
        // alternate peptide ends with it (VariationEffect.pm:375-393) and the
        // notation leaves out last residues the two share: the alternate
        // peptide's last codons, as many as the reference peptide has
        // residues, hold known bases alone, or the row's reading rests on the
        // unknown ones.
        let whole = codon.len() / 3 * 3;
        let alternate = perl_translate(&codon[..whole], table);
        let opens_the_alternate = alternate.starts_with(&reference);
        let compared = if opens_the_alternate {
            0
        } else {
            (reference.len() * 3).min(whole)
        };
        !reference.contains(&b'*') && is_unambiguous_dna(&codon[whole - compared..whole])
    }

    /// [`stop_codon_before_unknown_bases`] of the alternate allele's codon
    /// string on a row whose stop codon the stop codon test finds kept, where
    /// the reference peptide opens the alternate peptide. The residues
    /// through that stop are then the same under every completion of the
    /// unknown bases, and `inframe_insertion` reads no residue past it
    /// (VariationEffect.pm:1126-1127). Asked only of an evaluation that reads
    /// such a row ([`Self::reading_unknown_bases`]).
    fn keeps_a_stop_codon_before_unknown_bases(&mut self) -> bool {
        if !self.reads_a_kept_stop_codon {
            return false;
        }
        let table = self.table;
        let Some(codon) = self.codon(false) else {
            return false;
        };
        if !stop_codon_before_unknown_bases(&codon, table) {
            return false;
        }
        let reference = match self.peptide(true) {
            Some(reference) if reference == b"-" => Vec::new(),
            Some(reference) => reference,
            None => return false,
        };
        let whole = codon.len() / 3 * 3;
        perl_translate(&codon[..whole], table).starts_with(&reference)
    }

    /// The evaluation of a row whose alternate allele holds unknown bases
    /// where the row reads a peptide for it: its known bases settle a stop
    /// codon before the first of them. The stop is a gained one
    /// ([`Self::gains_a_stop_codon_before_unknown_bases`]) on a row no test of
    /// the start codon or of the stop codon decides. Those tests read the
    /// allele's letters themselves, so a row they give `start_lost`,
    /// `start_retained_variant` or `stop_lost` keeps their reading. So does a
    /// row whose settled stop ends the protein where the reference protein
    /// ends, with every residue before it the reference's
    /// ([`Self::ref_eq_alt_sequence`]), unless the stop codon test gives it
    /// `stop_retained_variant`: an insertion that test finds to keep the stop
    /// codon reads its peptide where the stop it spells stands before its
    /// first unknown base ([`Self::keeps_a_stop_codon_before_unknown_bases`])
    /// and the alternate protein, held against the reference protein, is
    /// unchanged ([`Self::stop_codon_change`]), and takes `inframe_insertion`
    /// as every completion of the unknown bases does. A row the peptide gives
    /// no `inframe_insertion` stays as the stop codon test gives it. The span
    /// is one the known bases settle
    /// ([`shared_ends_cut_at_known_bases`], of the alleles as the record
    /// writes them where the engines cut them;
    /// [`alleles_open_with_known_bases`] on a record of several alternate
    /// alleles).
    /// `None` for an allele of known bases and for every other row, whose
    /// evaluation is [`Self::new`]'s.
    fn reading_unknown_bases(
        variant: &'a InputVariant,
        transcript: &'a Transcript,
        span_start: u64,
        span_end: u64,
        fasta: Option<&'a vep_fasta::IndexedFasta>,
    ) -> Option<Self> {
        let alternate = variant.alt_allele();
        if alternate.is_empty() || is_perl_unambiguous_dna(alternate) {
            return None;
        }
        let mut unread = Self::new(variant, transcript, span_start, span_end, fasta)?;
        if is_perl_unambiguous_dna(unread.alt_allele) {
            return None;
        }
        let terms = unread.terms();
        let decided_by_a_codon_test = terms.iter().any(|term| {
            matches!(
                term,
                Consequence::StartLost | Consequence::StartRetainedVariant | Consequence::StopLost
            )
        });
        if decided_by_a_codon_test {
            return None;
        }
        // The stop codon test reads the edited sequence itself: where it
        // finds the stop codon kept, known bases spell a stop codon at the
        // old stop position behind the reference's residues. A span that
        // replaces reference bases is left to that test alone.
        let keeps_the_stop_codon = terms.contains(&Consequence::StopRetainedVariant);
        if keeps_the_stop_codon && unread.vf_end >= unread.vf_start {
            return None;
        }
        let written = variant
            .original_allele_string
            .as_deref()
            .map(|alleles| alleles.split('/').collect::<Vec<_>>());
        let alleles = std::iter::once(&variant.ref_allele)
            .chain(&variant.alt_alleles)
            .map(Vec::as_slice);
        if variant.alt_alleles.len() > 1 && !alleles_open_with_known_bases(alleles) {
            return None;
        }
        let settled = match written.as_deref() {
            None => shared_ends_cut_at_known_bases(unread.ref_allele, unread.alt_allele),
            Some([reference, alternate]) => {
                shared_ends_cut_at_known_bases(reference.as_bytes(), alternate.as_bytes())
            }
            Some(_) => false,
        };
        if !settled {
            return None;
        }
        let mut ev = Self::new(variant, transcript, span_start, span_end, fasta)?;
        if keeps_the_stop_codon {
            ev.reads_a_kept_stop_codon = true;
            ev.peptide(false)?;
            // `stop_retained` reads the peptides for a completion whose
            // peptide holds no `X` and the edited sequence for any other
            // (VariationEffect.pm:1313-1318), so the row is read only where
            // the two proteins themselves show the stop codon kept, and only
            // for the term the peptide adds.
            let protein_unchanged =
                matches!(ev.stop_codon_change(), Some(StopCodonChange::Unchanged));
            let terms = ev.terms();
            let reads_its_insertion = terms.contains(&Consequence::StopRetainedVariant)
                && terms.contains(&Consequence::InframeInsertion);
            return (protein_unchanged && reads_its_insertion).then_some(ev);
        }
        ev.reads_unknown_bases = true;
        ev.peptide(false)?;
        (!ev.ref_eq_alt_sequence()).then_some(ev)
    }

    /// [`Self::reading_unknown_bases`] where it applies, [`Self::new`]
    /// otherwise.
    fn for_row(
        variant: &'a InputVariant,
        transcript: &'a Transcript,
        span_start: u64,
        span_end: u64,
        fasta: Option<&'a vep_fasta::IndexedFasta>,
    ) -> Option<Self> {
        Self::reading_unknown_bases(variant, transcript, span_start, span_end, fasta)
            .or_else(|| Self::new(variant, transcript, span_start, span_end, fasta))
    }

    /// `_get_peptide_alleles`: `None` is Perl's empty list.
    fn peptide_alleles(&mut self) -> Option<(Vec<u8>, Vec<u8>)> {
        if let Some(m) = &self.peps_memo {
            return m.clone();
        }
        let alt_pep = self.peptide(false);
        let ref_pep = if alt_pep.is_some() {
            self.peptide(true)
        } else {
            None
        };
        let value = match (ref_pep, alt_pep) {
            (Some(r), Some(a)) if !r.is_empty() && !a.is_empty() => {
                let strip = |p: Vec<u8>| if p == b"-" { Vec::new() } else { p };
                Some((strip(r), strip(a)))
            }
            _ => None,
        };
        self.peps_memo = Some(value.clone());
        value
    }

    /// `_get_codon_alleles`: empty for a frameshift.
    fn codon_alleles(&mut self) -> Option<(Vec<u8>, Vec<u8>)> {
        if self.frameshift() {
            return None;
        }
        let alt_codon = self.codon(false)?;
        let ref_codon = self.codon(true)?;
        let strip = |c: Vec<u8>| if c == b"-" { Vec::new() } else { c };
        Some((strip(ref_codon), strip(alt_codon)))
    }

    fn partial_codon(&self) -> bool {
        let translation_start = self.span.tl_start.and_then(|t| u64::try_from(t).ok());
        crate::mapper::partial_codon(self.transcript, translation_start)
    }

    /// Whether a stop codon stands at the coding end for the stop terms to
    /// read: a whole one, or one that the coding sequence's last one or two
    /// bases begin.
    fn stop_codon_stands(&self) -> bool {
        self.stop_codon_positions() > 0
    }

    /// How many of the last coding positions are bases of the stop codon at
    /// the coding end: three of a whole stop codon, the one or two bases
    /// after the last whole codon of one that polyadenylation completes, none
    /// where no stop codon stands. A span that touches the last whole codon
    /// alone overlaps no stop codon. A record that shows its model complete
    /// reaches the evaluation as its completed copy
    /// ([`completed_by_polyadenylation`]), whose stop codon is whole.
    fn stop_codon_positions(&self) -> i64 {
        match self.stop_codon {
            StopCodonAtTheCodingEnd::None => 0,
            StopCodonAtTheCodingEnd::Whole => 3,
            StopCodonAtTheCodingEnd::BegunByTheTail { bases } => i64::from(bases),
        }
    }

    /// Whether a frameshift whose first changed codon the transcript's
    /// sequence ends inside is named at its first residue, `<Ref><n>fs`, in
    /// place of the deletion of its reference residues that Ensembl VEP's
    /// frameshift routine words (TranscriptVariationAllele.pm:2342-2346): no
    /// stop codon stands at the coding end, and Ensembl VEP reads one lost
    /// there ([`Self::overlaps_the_last_three_coding_positions`]), so it
    /// adds an extension of unknown length to that deletion (1916-1936).
    /// The frame shifts and one or two bases are left of the first codon
    /// the change reaches, so the sequence shows no stop lost, no residue
    /// deleted and none in its place. The HGVS nomenclature describes a
    /// deletion that shifts the frame as a frameshift, and its short form,
    /// a residue and `fs`, marks where the frame breaks. Where Ensembl VEP
    /// reads no lost stop, on a transcript that carries `cds_end_NF` and
    /// for an allele written before the last three coding positions that
    /// the 3' shift carries onto them, the deletion stands as it writes it.
    fn cut_short_frameshift_names_its_first_residue(&self) -> bool {
        !self.stop_codon_stands() && self.overlaps_the_last_three_coding_positions()
    }

    /// `_overlaps_stop_codon` as Ensembl VEP reads it
    /// (VariationEffect.pm:1373-1381): the span's cDNA coordinates against
    /// the last three coding positions, whatever they hold, on a transcript
    /// without `cds_end_NF`.
    fn overlaps_the_last_three_coding_positions(&self) -> bool {
        if self
            .transcript
            .flags
            .iter()
            .any(|flag| &**flag == "cds_end_NF")
        {
            return false;
        }
        let (Some(cs), Some(ce)) = (truthy(self.span.cdna_start), truthy(self.span.cdna_end))
        else {
            return false;
        };
        perl_overlap(cs, ce, self.cdna_coding_end - 2, self.cdna_coding_end)
    }

    /// `_overlaps_stop_codon` (VariationEffect.pm:1363-1385): the span's cDNA
    /// coordinates against the coding positions of the stop codon
    /// ([`Self::stop_codon_positions`]), which Ensembl VEP reads as the last
    /// three on every transcript. Ensembl VEP refuses the overlap for
    /// `cds_end_NF` alone (:1373); here no span overlaps a stop codon where
    /// none stands at the coding end.
    fn overlaps_stop_codon(&mut self) -> bool {
        if let Some(v) = self.overlaps_stop_memo {
            return v;
        }
        self.overlaps_stop_memo = Some(false);
        let positions = self.stop_codon_positions();
        if positions == 0 {
            return false;
        }
        let (Some(cs), Some(ce)) = (truthy(self.span.cdna_start), truthy(self.span.cdna_end))
        else {
            return false;
        };
        let v = perl_overlap(
            cs,
            ce,
            self.cdna_coding_end - (positions - 1),
            self.cdna_coding_end,
        );
        self.overlaps_stop_memo = Some(v);
        v
    }

    /// The stop codon overlap behind the second arm of `stop_retained`
    /// (VariationEffect.pm:1318). A span with reference bases is read at the
    /// stop codon's own genomic bases ([`overlaps_stop_codon_bases`]), an
    /// insertion on the coding sequence
    /// ([`Self::insertion_reaches_stop_codon`]). Ensembl VEP's
    /// `_overlaps_stop_codon_cil` (:1389-1431) moves an insertion's start by
    /// the inserted length, forward on the forward strand and back on the
    /// reverse (:1408-1415), so one coding position overlaps on the reverse
    /// strand from `L - 1` bases before the codon whatever is inserted, and on
    /// the forward strand only between the codon's first and second base.
    fn overlaps_stop_codon_cil(&mut self) -> bool {
        if let Some(v) = self.overlaps_stop_cil_memo {
            return v;
        }
        self.overlaps_stop_cil_memo = Some(false);
        let Some(span) = self.stop_codon_cil_span() else {
            return false;
        };
        if self.coding_region.is_none() {
            return false;
        }
        let Some(mapper) = self
            .transcript
            .vefc
            .as_ref()
            .and_then(|v| v.mapper.as_ref())
        else {
            return false;
        };
        let positions = self.stop_codon_positions();
        let v = if self.vf_end < self.vf_start {
            // Before a stop codon that the last one or two coding bases
            // begin, an insertion inside the last whole codon reaches no
            // base of the stop codon.
            let at_the_stop_codon = positions == 3
                || truthy(self.span.cds_start)
                    .is_some_and(|follows| follows > self.cds.len() as i64 - positions);
            at_the_stop_codon && self.insertion_reaches_stop_codon()
        } else {
            overlaps_stop_codon_bases(
                &mapper.exon_coord_mapper.pairs,
                (self.cdna_coding_start, self.cdna_coding_end),
                positions,
                span,
            )
        };
        self.overlaps_stop_cil_memo = Some(v);
        v
    }

    /// Whether an insertion lies before or inside the stop codon and leaves
    /// every residue before the codon as the reference has it. `f` is the
    /// coding position that follows the insertion (the CDS start, where the
    /// edit of VariationEffect.pm:1521 is anchored) and `n` the coding length.
    /// The test is `f <= n` and, in the sequence that edit gives
    /// ([`Self::edited_stop_codon`]), the reference's residue at every codon
    /// before the old stop position, the same for one coding position on
    /// either strand. Inserted bases that repeat the coding bases from `f` up
    /// to the codon leave every base before it in place; bases that differ
    /// keep a codon's residue where both codons translate to one residue.
    fn insertion_reaches_stop_codon(&self) -> bool {
        let Some(f) = truthy(self.span.cds_start) else {
            return false;
        };
        let n = self.cds.len() as i64;
        if f < 1 || f > n {
            return false;
        }
        if self.feature_seq_bases().is_empty() {
            return f >= n - 2;
        }
        matches!(
            self.edited_stop_codon(),
            EditedStopCodon::Read {
                residues_kept: true,
                ..
            }
        )
    }

    /// The span `_overlaps_stop_codon_cil` tests (VariationEffect.pm:1400-1415):
    /// the variant's, an inserted allele's length added to its start on the
    /// forward strand and subtracted on the reverse; `None` for a span without
    /// coordinates and where no stop codon stands at the coding end, which
    /// Ensembl VEP reads from `cds_end_NF` alone (:1400).
    fn stop_codon_cil_span(&self) -> Option<(i64, i64)> {
        if !self.stop_codon_stands() {
            return None;
        }
        let (mut v_start, v_end) = (self.vf_start, self.vf_end);
        if v_start == 0 || v_end == 0 {
            return None;
        }
        let inserted = &self.feature_seq;
        if v_end < v_start
            && !inserted.is_empty()
            && inserted
                .iter()
                .all(|b| matches!(b, b'A' | b'C' | b'T' | b'G' | b'N'))
        {
            let len = inserted.len() as i64;
            if self.strand == 1 {
                v_start += len;
            } else {
                v_start -= len;
            }
        }
        Some((v_start, v_end))
    }

    /// Whether Ensembl VEP's `stop_retained` holds on this allele for want of
    /// a test alone. An insertion or deletion whose alternate peptide is
    /// undefined or partial, whose span overlaps the last three genomic
    /// positions of the coding region (VariationEffect.pm:1418-1427) and
    /// lacks a cDNA coordinate at one end or a CDS start, leaves the stop
    /// codon edit test before it edits anything (:1508), and :1318 reads a
    /// stop codon kept. `false` where the sequence the record gives is read
    /// instead ([`Self::removed_stop_codon`]). The overlap is the one Ensembl
    /// VEP reads, not the codon's own bases: on a stop codon an intron splits,
    /// a span that takes the stop bases before the intron and none of those
    /// three positions is no stop codon kept there, so no row takes a
    /// notation that neither engine writes. Under `--shift_3prime 1` every
    /// term and every notation is read at the site the flag moves the edit to.
    fn untested_stop_codon_kept_in_ensembl_vep(&mut self) -> bool {
        if !(self.insertion || self.deletion)
            || self.partial_codon()
            || self.stop_lost()
            || self.removed_stop_codon().is_some()
        {
            return false;
        }
        if matches!(
            self.peptide_alleles(),
            Some((_, alt_pep)) if !alt_pep.is_empty() && !alt_pep.contains(&b'X')
        ) {
            return false;
        }
        if truthy(self.span.cdna_start).is_some()
            && truthy(self.span.cdna_end).is_some()
            && truthy(self.span.cds_start).is_some()
        {
            return false;
        }
        let (Some((v_start, v_end)), Some((crs, cre))) =
            (self.stop_codon_cil_span(), self.coding_region)
        else {
            return false;
        };
        let (stop_start, stop_end) = if self.strand == -1 {
            (crs, crs + 2)
        } else {
            (cre - 2, cre)
        };
        perl_overlap(v_start, v_end, stop_start, stop_end)
    }

    /// `_overlaps_start_codon` (VariationEffect.pm:965-990): the span's cDNA
    /// coordinates against the first three coding positions. Ensembl VEP
    /// refuses the overlap for `cds_start_NF` alone (:975); here no span
    /// overlaps a start codon where none stands at the coding start.
    fn overlaps_start_codon(&mut self) -> bool {
        if let Some(v) = self.overlaps_start_memo {
            return v;
        }
        self.overlaps_start_memo = Some(false);
        if !self.start_codon_stands {
            return false;
        }
        let (Some(cs), Some(ce)) = (truthy(self.span.cdna_start), truthy(self.span.cdna_end))
        else {
            return false;
        };
        let v = perl_overlap(cs, ce, self.cdna_coding_start, self.cdna_coding_start + 2);
        self.overlaps_start_memo = Some(v);
        v
    }

    fn feature_seq_bases(&self) -> &[u8] {
        if self.feature_seq == b"-" {
            b""
        } else {
            &self.feature_seq
        }
    }

    /// `_ins_del_stop_altered`: the cDNA-coordinate overlap gate, then the edit test.
    fn ins_del_stop_altered(&mut self) -> bool {
        if let Some(v) = self.ins_del_stop_altered_memo {
            return v;
        }
        self.ins_del_stop_altered_memo = Some(false);
        if !is_perl_unambiguous_dna(self.alt_allele) {
            return false;
        }
        if !self.overlaps_stop_codon() {
            return false;
        }
        let v = self.edited_stop_codon_altered();
        self.ins_del_stop_altered_memo = Some(v);
        v
    }

    /// The edit test of `_ins_del_stop_altered_cil` (VariationEffect.pm:1488-1539)
    /// behind the genomic-coordinate gate of [`Self::overlaps_stop_codon_cil`],
    /// as `stop_retained` reads it: `true` only where the test is made, the
    /// three bases the edited sequence holds at the old stop position are a
    /// stop codon and every codon before them reads the reference's residue.
    /// Ensembl VEP's routine returns 0 before it edits anything for an allele
    /// holding a letter that names no base (:1498) and for a span without
    /// cDNA or CDS coordinates (:1508), and :1318 reads that 0 as a stop codon
    /// kept; here a letter counts only among the bases read, and a test that
    /// cannot be made keeps nothing. Ensembl VEP reads the three bases alone
    /// (:1528-1535), so a deletion that takes coding bases with the stop codon
    /// and leaves a stop codon of the 3' UTR at the old position is kept there.
    fn edited_stop_codon_kept_cil(&mut self) -> bool {
        if let Some(v) = self.edited_stop_codon_kept_cil_memo {
            return v;
        }
        self.edited_stop_codon_kept_cil_memo = Some(false);
        if !self.overlaps_stop_codon_cil() {
            return false;
        }
        let v = match self.edited_stop_codon() {
            EditedStopCodon::Read {
                codon,
                residues_kept,
            } => {
                residues_kept
                    && codon
                        .iter()
                        .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
                    && perl_translate_codon(&codon, self.table) == b'*'
            }
            EditedStopCodon::NotMade | EditedStopCodon::Shorter => false,
        };
        self.edited_stop_codon_kept_cil_memo = Some(v);
        v
    }

    /// The verdict of `_ins_del_stop_altered` after its overlap gate
    /// (VariationEffect.pm:1453-1484): `true` when the edited sequence is
    /// shorter than the CDS or the codon left at the old stop position is not a
    /// stop, `false` where the test is not made, which `stop_lost` (:1263)
    /// reads as a stop codon that is not lost.
    fn edited_stop_codon_altered(&self) -> bool {
        match self.edited_stop_codon() {
            EditedStopCodon::NotMade => false,
            EditedStopCodon::Shorter => true,
            EditedStopCodon::Read { codon, .. } => perl_translate_codon(&codon, self.table) != b'*',
        }
    }

    /// The edit `_ins_del_stop_altered` and its `_cil` twin share after their
    /// overlap gates (VariationEffect.pm:1521-1529): the allele written into
    /// the translateable sequence plus 3' UTR at the CDS start, over the cDNA
    /// length of the span, then the three bases left at the old stop position.
    fn edited_stop_codon(&self) -> EditedStopCodon {
        if !(self.insertion || self.deletion) {
            return EditedStopCodon::NotMade;
        }
        let (Some(cdna_start), Some(cdna_end), Some(cds_start)) = (
            truthy(self.span.cdna_start),
            truthy(self.span.cdna_end),
            truthy(self.span.cds_start),
        ) else {
            return EditedStopCodon::NotMade;
        };
        let mut s = self.cds.to_vec();
        if let Some(utr) = self.utr3.as_ref() {
            s.extend_from_slice(utr);
        }
        let fs = self.feature_seq_bases();
        perl_substr_assign(&mut s, cds_start - 1, cdna_end - cdna_start + 1, fs);
        if s.len() < self.cds.len() {
            return EditedStopCodon::Shorter;
        }
        if self.cds.len() < 3 {
            return EditedStopCodon::NotMade;
        }
        let stop_idx = self.cds.len() - 3;
        EditedStopCodon::Read {
            codon: [s[stop_idx], s[stop_idx + 1], s[stop_idx + 2]],
            residues_kept: residues_before_stop_kept(self.cds, &s, cds_start, self.table),
        }
    }

    /// What the sequence the record gives shows at the old stop position, for
    /// an insertion's or deletion's span that holds all three bases of the stop
    /// codon and lacks a cDNA coordinate at one end or a CDS start, so that
    /// the stop codon edit test cannot be made; `None` for any other span.
    /// The condition is Ensembl VEP's own, `return 0 unless $cdna_start &&
    /// $cdna_end && $cds_start;` (VariationEffect.pm:1508, and :1453 behind
    /// `stop_lost`): a CDS end alone does not make the test. Ensembl VEP
    /// writes `stop_retained_variant` on every such row (:1508 read at
    /// :1318). Here the stop codon is lost unless it is shown kept, and
    /// nothing is read of bases no record gives:
    ///
    /// 1. A deletion that can be written where it takes no coding base and no
    ///    stop base keeps the codon.
    /// 2. Where, at one of the places the edit can be written, an exon base
    ///    stands at or beside each of its ends, the transcript's sequence
    ///    across the edit is given: the exon bases before it, the alternate
    ///    bases and the exon bases after it. The codon is kept where that
    ///    sequence holds a stop codon at the old stop position behind the
    ///    reference's residues, lost where it holds none or ends before it.
    /// 3. Otherwise the edit takes a splice site or the transcript's end with
    ///    exon bases, and only the exon bases before it and, behind an exon
    ///    base, the alternate bases are given; they are read the same way.
    /// 4. A record whose two alleles share every reference base removes
    ///    nothing and is read as the bases it adds
    ///    ([`Self::inserted_alone_stop_codon`]).
    ///
    /// On a model whose stop codon an `_rna_edit` or polyadenylation
    /// completes, one or two of the codon's positions have no genomic base.
    /// A span holds the codon there where it holds every base the genome has
    /// of it, and the bases that complete the codon stand in the transcript's
    /// sequence after its last exon base.
    ///
    /// On the transcript `completed_by_polyadenylation` makes, the adenines
    /// of the poly(A) tail follow the completed codon, and no span removes
    /// them: they follow what an edit that takes the transcript's end leaves.
    /// So a span that starts inside the codon there and holds every base the
    /// genome has of it from there on is read too, on the codon's first bases
    /// and the adenines after them. Only that transcript holds bases past its
    /// exons without an `_rna_edit`; a record whose own edit completes the
    /// codon is read on the bases its exons hold.
    ///
    /// Three bounds are vep-rs's own. The places a deletion can be written at
    /// are read from the reference sequence up to [`EQUIVALENT_SITES_LIMIT`]
    /// positions each way; a place further off is not read, so a deletion
    /// that takes no coding base only there is read at the places within
    /// reach and is `stop_lost` unless one of them shows the codon kept.
    /// Without a reference sequence the record's own place alone is read.
    /// A base that is none of `A`, `C`, `G` and `T` ends the search
    /// ([`Self::equivalent_sites`]).
    fn removed_stop_codon(&mut self) -> Option<GivenStopCodon> {
        if let Some(v) = self.removed_stop_codon_memo {
            return v;
        }
        let v = self.given_stop_codon();
        self.removed_stop_codon_memo = Some(v);
        v
    }

    fn given_stop_codon(&self) -> Option<GivenStopCodon> {
        if !(self.insertion || self.deletion) || !self.stop_codon_stands() || self.vf_start < 1 {
            return None;
        }
        if truthy(self.span.cdna_start).is_some()
            && truthy(self.span.cdna_end).is_some()
            && truthy(self.span.cds_start).is_some()
        {
            return None;
        }
        let mapper = self.transcript.vefc.as_ref()?.mapper.as_ref()?;
        let strand = self.transcript.strand.as_i8();
        let (coding_start, coding_end) = (self.cdna_coding_start, self.cdna_coding_end);
        let axis = ExonAxis::new(
            &mapper.exon_coord_mapper.pairs,
            strand,
            (coding_start, coding_end),
        )?;
        let (mut s, mut e) = if strand == 1 {
            (self.vf_start, self.vf_end)
        } else {
            (-self.vf_end, -self.vf_start)
        };
        if coding_start < 1 || coding_end - 2 < coding_start {
            return None;
        }
        // A stop codon that an `_rna_edit` or polyadenylation completes has
        // positions without a genomic base, which no span removes: the span
        // holds the codon where it holds every base the genome has of it.
        let on_the_genome: Vec<i64> = (coding_end - 2..=coding_end)
            .filter_map(|cdna| axis.position(cdna))
            .collect();
        // The adenines of a poly(A) tail follow the completed codon on the
        // transcript `completed_by_polyadenylation` makes, which alone holds
        // bases past its exons without an `_rna_edit`. A span that starts
        // inside the codon there holds it from that base on.
        let poly_a_tail_follows = on_the_genome.len() < 3
            && self.utr3.as_deref().is_some_and(|tail| !tail.is_empty())
            && !self
                .transcript
                .attributes
                .iter()
                .any(|attribute| &*attribute.code == "_rna_edit");
        let holds_the_codon = on_the_genome.iter().any(|&base| s <= base && base <= e)
            && on_the_genome
                .iter()
                .all(|&base| base <= e && (s <= base || poly_a_tail_follows));
        if !holds_the_codon {
            return None;
        }
        // The transcript's sequence: the 5' UTR, the coding sequence without
        // the filler a start phase puts before it, and the 3' UTR.
        let filler =
            usize::try_from(self.cds.len() as i64 - (coding_end - coding_start + 1)).ok()?;
        let utr5 = self.utr5.as_deref().unwrap_or_default();
        if utr5.len() as i64 != coding_start - 1 {
            return None;
        }
        let mut cdna_seq = utr5.to_vec();
        cdna_seq.extend_from_slice(&self.cds[filler..]);
        cdna_seq.extend_from_slice(self.utr3.as_deref().unwrap_or_default());
        // The edit, with the bases its two alleles share taken off both ends.
        let reference = if self.opposite_strand {
            reverse_complement(self.ref_allele)
        } else {
            self.ref_allele.to_vec()
        };
        let (mut removed, mut alt): (&[u8], &[u8]) = (&reference, self.feature_seq_bases());
        if removed.len() as i64 == e - s + 1 {
            while let (Some(r), Some(a)) = (removed.first(), alt.first()) {
                if !r.eq_ignore_ascii_case(a) {
                    break;
                }
                (removed, alt, s) = (&removed[1..], &alt[1..], s + 1);
            }
            while let (Some(r), Some(a)) = (removed.last(), alt.last()) {
                if !r.eq_ignore_ascii_case(a) {
                    break;
                }
                (removed, alt, e) = (&removed[..removed.len() - 1], &alt[..alt.len() - 1], e - 1);
            }
        }
        if e < s {
            return (!alt.is_empty())
                .then(|| self.inserted_alone_stop_codon(&axis, &cdna_seq, s, alt, filler));
        }
        let (along, back) = if alt.is_empty() {
            self.equivalent_sites(s, e, strand)
        } else {
            (0, 0)
        };
        // A place takes the fewest coding bases, then the fewest exon bases,
        // then lies furthest along the transcript.
        let key = |j: &i64| {
            let (coding, exonic) = axis.taken(s + j, e + j);
            (coding, exonic, -j)
        };
        let best = (-back..=along).min_by_key(key).unwrap_or(0);
        if axis.taken(s + best, e + best).0 == 0 {
            return Some(GivenStopCodon::Retained);
        }
        let before = |pos: i64| axis.exonic_before(pos).min(cdna_seq.len());
        let across = (-back..=along)
            .filter(|j| {
                let (first, last) = (s + j, e + j);
                (axis.cdna(first - 1).is_some() || axis.cdna(first).is_some())
                    && (axis.cdna(last).is_some() || axis.cdna(last + 1).is_some())
            })
            .min_by_key(key);
        let given = match across {
            Some(j) => [
                &cdna_seq[..before(s + j)],
                alt,
                &cdna_seq[before(e + j + 1)..],
            ]
            .concat(),
            None => {
                let mut given = cdna_seq[..before(s + best)].to_vec();
                if axis.cdna(s + best - 1).is_some() {
                    given.extend_from_slice(alt);
                }
                // The completed codon's last bases and the poly(A) tail have
                // no genomic position: they follow what an edit that takes
                // the transcript's end leaves.
                let after = before(e + best + 1);
                if poly_a_tail_follows && after == before(i64::MAX) {
                    given.extend_from_slice(&cdna_seq[after..]);
                }
                given
            }
        };
        let (coding_start, coding_end) = (coding_start as usize, coding_end as usize);
        if given.len() < coding_end || cdna_seq.len() < coding_end {
            return Some(GivenStopCodon::Lost);
        }
        let (given, reference) = (&given[..coding_end], &cdna_seq[..coding_end]);
        if given.eq_ignore_ascii_case(reference) {
            return Some(GivenStopCodon::Retained);
        }
        if perl_translate_codon(&given[coding_end - 3..], self.table) != b'*' {
            return Some(GivenStopCodon::Lost);
        }
        let kept = given[..coding_start - 1].eq_ignore_ascii_case(&reference[..coding_start - 1])
            && same_residues(
                &given[coding_start - 1..coding_end - 3],
                &reference[coding_start - 1..coding_end - 3],
                filler,
                self.table,
            );
        Some(if kept {
            GivenStopCodon::Retained
        } else {
            GivenStopCodon::Changed
        })
    }

    /// The reading of [`Self::removed_stop_codon`] for a record whose two
    /// alleles share every reference base: no base is removed, and `alt`, the
    /// alternate's further bases, stands before the axis position `s`. Away
    /// from the exons, in the 5' UTR and after the stop codon those bases
    /// leave the coding sequence and the codon as they were. Inside the
    /// coding sequence the transcript's sequence is given whole: the codon is
    /// kept where a stop codon stands at the old stop position behind the
    /// reference's residues and lost where bases written between the codon's
    /// own leave none there; bases written before the codon that leave none
    /// there change the frame or a residue, which is neither.
    fn inserted_alone_stop_codon(
        &self,
        axis: &ExonAxis,
        cdna_seq: &[u8],
        s: i64,
        alt: &[u8],
        filler: usize,
    ) -> GivenStopCodon {
        let (coding_start, coding_end) = (
            self.cdna_coding_start as usize,
            self.cdna_coding_end as usize,
        );
        let point = axis.exonic_before(s).min(cdna_seq.len());
        let beside_exon = axis.cdna(s - 1).is_some() || axis.cdna(s).is_some();
        if !beside_exon
            || point < coding_start
            || point >= coding_end
            || cdna_seq.len() < coding_end
        {
            return GivenStopCodon::Retained;
        }
        let given = [&cdna_seq[..point], alt, &cdna_seq[point..]].concat();
        let (given, reference) = (&given[..coding_end], &cdna_seq[..coding_end]);
        if given.eq_ignore_ascii_case(reference) {
            return GivenStopCodon::Retained;
        }
        if perl_translate_codon(&given[coding_end - 3..], self.table) != b'*' {
            return if point > coding_end - 3 {
                GivenStopCodon::Lost
            } else {
                GivenStopCodon::Changed
            };
        }
        let kept = same_residues(
            &given[coding_start - 1..coding_end - 3],
            &reference[coding_start - 1..coding_end - 3],
            filler,
            self.table,
        );
        if kept {
            GivenStopCodon::Retained
        } else {
            GivenStopCodon::Changed
        }
    }

    /// How far a deletion of the axis positions `s..=e` can be written along
    /// the transcript and back against it without changing the edited genome,
    /// read from the reference sequence up to [`EQUIVALENT_SITES_LIMIT`]
    /// positions each way; `(0, 0)` without a reference sequence. A base that
    /// is none of `A`, `C`, `G` and `T` ends the walk: an unknown base is no
    /// evidence that the bases repeat, so no place is read beyond one.
    fn equivalent_sites(&self, s: i64, e: i64, strand: i8) -> (i64, i64) {
        let Some(fasta) = self.fasta else {
            return (0, 0);
        };
        let base = |pos: i64| {
            let genomic = if strand == 1 { pos } else { -pos };
            u64::try_from(genomic)
                .ok()
                .filter(|&g| g > 0)
                .and_then(|g| fasta.base(self.chr, g))
        };
        let repeats = |from: i64, onto: i64, step: i64| {
            (0..EQUIVALENT_SITES_LIMIT)
                .take_while(|i| {
                    matches!(
                        (base(from + step * i), base(onto + step * i)),
                        (Some(x), Some(y))
                            if x == y && matches!(x.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
                    )
                })
                .count() as i64
        };
        (repeats(s, e + 1, 1), repeats(e, s - 1, -1))
    }

    /// `_ins_del_start_altered` (VariationEffect.pm:1028-1075): the allele
    /// replaces its span of the 5' UTR and the coding sequence (:1058), and
    /// the start is kept where the UTR is unchanged and the three bases at
    /// the start position read `ATG` (:1062-1065). Ensembl VEP makes that
    /// test only behind a 5' UTR (:1061) and reads every deletion on a
    /// transcript without one as an altered start (:1069). Here the test is
    /// made on a transcript whose coding sequence opens on its first base as
    /// well, at that base, where no edit changes the UTR.
    fn ins_del_start_altered(&mut self) -> bool {
        if let Some(v) = self.ins_del_start_altered_memo {
            return v;
        }
        self.ins_del_start_altered_memo = Some(false);
        if !is_perl_unambiguous_dna(self.alt_allele) {
            return false;
        }
        if !self.overlaps_start_codon() {
            return false;
        }
        if !(self.insertion || self.deletion) {
            return false;
        }
        let (Some(cdna_start), Some(cdna_end)) =
            (truthy(self.span.cdna_start), truthy(self.span.cdna_end))
        else {
            return false;
        };
        let utr = self.utr5.clone().filter(|u| !u.is_empty());
        let mut s: Vec<u8> = utr.clone().unwrap_or_default();
        s.extend_from_slice(self.cds);
        let fs = self.feature_seq_bases().to_vec();
        perl_substr_assign(&mut s, cdna_start - 1, cdna_end - cdna_start + 1, &fs);
        // The empty UTR of a transcript whose coding sequence opens on its
        // first base is tested as a UTR is; Ensembl VEP skips the test there.
        let tested_utr: Option<&[u8]> = match utr.as_deref() {
            None if self.cdna_coding_start == 1 => Some(&[]),
            other => other,
        };
        if let Some(utr) = tested_utr {
            let atg_start = utr.len();
            let new_sc = s
                .get(atg_start..(atg_start + 3).min(s.len()))
                .unwrap_or(b"");
            let new_utr = &s[..atg_start.min(s.len())];
            if new_utr == utr && new_sc == b"ATG" {
                return false;
            }
        }
        if s.len() < self.cds.len() {
            self.ins_del_start_altered_memo = Some(true);
            return true;
        }
        let v = self.cds != &s[s.len() - self.cds.len()..];
        self.ins_del_start_altered_memo = Some(v);
        v
    }

    fn inv_start_altered(&mut self) -> bool {
        if let Some(v) = self.inv_start_altered_memo {
            return v;
        }
        self.inv_start_altered_memo = Some(false);
        if !is_perl_unambiguous_dna(self.alt_allele) {
            return false;
        }
        if !self.overlaps_start_codon() {
            return false;
        }
        let (Some(cdna_start), Some(cdna_end)) =
            (truthy(self.span.cdna_start), truthy(self.span.cdna_end))
        else {
            return false;
        };
        let Some(utr) = self.utr5.clone().filter(|u| !u.is_empty()) else {
            return false;
        };
        let mut s = utr.clone();
        s.extend_from_slice(self.cds);
        if cdna_end > s.len() as i64 {
            return false;
        }
        let fs = self.feature_seq_bases().to_vec();
        let atg_start = utr.len();
        perl_substr_assign(&mut s, cdna_start - 1, cdna_end - cdna_start + 1, &fs);
        let new_sc = s
            .get(atg_start..(atg_start + 3).min(s.len()))
            .unwrap_or(b"");
        let v = new_sc != b"ATG";
        self.inv_start_altered_memo = Some(v);
        v
    }

    fn snp_start_altered(&mut self) -> bool {
        if let Some(v) = self.snp_start_altered_memo {
            return v;
        }
        self.snp_start_altered_memo = Some(true);
        if !is_perl_unambiguous_dna(self.alt_allele) {
            return false;
        }
        let (Some(cdna_start), Some(cdna_end)) =
            (truthy(self.span.cdna_start), truthy(self.span.cdna_end))
        else {
            return false;
        };
        let mut s: Vec<u8> = self.utr5.clone().unwrap_or_default();
        s.extend_from_slice(self.cds);
        let fs = self.feature_seq_bases().to_vec();
        perl_substr_assign(&mut s, cdna_start - 1, cdna_end - cdna_start + 1, &fs);
        let tail_start = s.len().saturating_sub(self.cds.len());
        let aa = s
            .get(tail_start..(tail_start + 3).min(s.len()))
            .unwrap_or(b"");
        if aa == b"ATG" {
            self.snp_start_altered_memo = Some(false);
        }
        self.snp_start_altered_memo.unwrap_or(true)
    }

    fn start_retained_variant(&mut self) -> bool {
        if !self.overlaps_start_codon() {
            return false;
        }
        if self.snp {
            !self.snp_start_altered()
        } else {
            !self.ins_del_start_altered()
        }
    }

    /// Whether an insertion or deletion at the start codon leaves the whole
    /// coding sequence standing, its first codon moved off the coding
    /// start: the edited 5' UTR and coding sequence end in the
    /// coding sequence (the last test of `_ins_del_start_altered`,
    /// VariationEffect.pm:1071), and the transcript's first codon does not
    /// stand after an unchanged 5' UTR (its first test, 1061-1065, asks that
    /// of the triplet `ATG`). The edited transcript is the transcript with a
    /// longer or shorter 5' UTR, and translation from its start codon gives
    /// the reference protein. The allele is read in place of the span's
    /// bases, as `_ins_del_start_altered` writes it (1058), on a span that
    /// holds bases of an edited model that the genome lacks too.
    fn keeps_the_coding_sequence_whole_behind_the_allele(&mut self) -> bool {
        if !(self.insertion || self.deletion)
            || !is_perl_unambiguous_dna(self.alt_allele)
            || !self.start_retained_variant()
        {
            return false;
        }
        let (Some(cdna_start), Some(cdna_end)) =
            (truthy(self.span.cdna_start), truthy(self.span.cdna_end))
        else {
            return false;
        };
        let utr = self.utr5.clone().unwrap_or_default();
        let mut edited = utr.clone();
        edited.extend_from_slice(self.cds);
        let allele = self.feature_seq_bases().to_vec();
        perl_substr_assign(
            &mut edited,
            cdna_start - 1,
            cdna_end - cdna_start + 1,
            &allele,
        );
        let first_codon_at_the_coding_start = self.cds.len() >= 3
            && edited.get(..utr.len()) == Some(utr.as_slice())
            && edited.get(utr.len()..utr.len() + 3) == self.cds.get(..3);
        !first_codon_at_the_coding_start && edited.ends_with(self.cds)
    }

    fn start_lost(&mut self) -> bool {
        if let Some(v) = self.start_lost_memo {
            return v;
        }
        self.start_lost_memo = Some(false);
        if !self.overlaps_start_codon() {
            return false;
        }
        if self.ins_del_start_altered() && !(self.inframe_insertion() || self.inframe_deletion()) {
            self.start_lost_memo = Some(true);
            return true;
        }
        if self.inv_start_altered() {
            self.start_lost_memo = Some(true);
            return true;
        }
        let Some((ref_pep, alt_pep)) = self.peptide_alleles() else {
            return false;
        };
        if ref_pep.is_empty() {
            return false;
        }
        if alt_pep.is_empty() {
            let v = self.start_codon_removed_without_five_prime_utr();
            self.start_lost_memo = Some(v);
            return v;
        }
        if alt_pep == b"X" {
            return false;
        }
        let v = self.span.tl_start == Some(1)
            && !alt_pep.ends_with(&ref_pep)
            && !alt_pep.starts_with(&ref_pep);
        self.start_lost_memo = Some(v);
        v
    }

    /// Whether whole codons removed from the first coding base take the start
    /// codon of a transcript that has no 5' UTR, asked where the span overlaps
    /// the start codon and the alternate peptide is empty.
    ///
    /// Ensembl VEP writes no start term on such a row: its length test, which
    /// reads the start as altered (VariationEffect.pm:1069), is set aside for
    /// an in-frame deletion (:870), `_inv_start_altered` returns before its
    /// edit on a transcript without a 5' UTR (:926) and the peptide test
    /// returns on the empty alternate peptide (:875). Here the coding sequence
    /// is edited as `_inv_start_altered` edits it behind a 5' UTR (:936-941),
    /// and the start codon is lost unless the three bases at its start code
    /// methionine under the transcript's genetic code: `ATG`, the triplet
    /// that routine asks for, or a triplet that keeps the first residue as a
    /// substitution that both engines read as no lost start does (:880-881).
    /// A triplet that holds a base other than `A`, `C`, `G` and `T` settles
    /// nothing. The test is made where the first codon is a start codon of
    /// that code that codes methionine and the removal starts on the first
    /// coding base. A coding sequence that opens on any other codon is left
    /// as Ensembl VEP writes it: one cut off at its start opens on any
    /// codon. So is a removal that begins inside the start codon and leaves
    /// no whole codon of a coding sequence that ends inside one.
    fn start_codon_removed_without_five_prime_utr(&self) -> bool {
        if self.utr5.as_ref().is_some_and(|utr| !utr.is_empty()) {
            return false;
        }
        let residue = |codon: &[u8]| translate_codon_with_table(codon, self.table);
        if !self
            .cds
            .get(..3)
            .is_some_and(|first| is_start_codon(first, self.table) && residue(first) == b'M')
        {
            return false;
        }
        let (Some(1), Some(cdna_end)) = (self.span.cdna_start, truthy(self.span.cdna_end)) else {
            return false;
        };
        let mut edited = self.cds.to_vec();
        perl_substr_assign(&mut edited, 0, cdna_end, self.feature_seq_bases());
        // `X` is the residue of a triplet that holds a base that is none of
        // the four.
        !matches!(edited.get(..3).map(residue), Some(b'M' | b'X'))
    }

    /// `stop_lost` (VariationEffect.pm:1234-1288): never on a partial codon
    /// (:1240, ahead of the cache); the peptide comparison decides unless the
    /// alternate peptide carries an `X` (:1259), when the edited stop codon of
    /// `_ins_del_stop_altered` decides instead. That routine returns 0 for a
    /// span without coordinates (:1376, :1453), which :1263 reads as a stop
    /// codon that is not lost although the span may hold the whole codon;
    /// such a span is lost unless the sequence the record gives shows the
    /// stop codon kept ([`Self::removed_stop_codon`]).
    fn stop_lost(&mut self) -> bool {
        if self.partial_codon() {
            return false;
        }
        if let Some(v) = self.stop_lost_memo {
            return v;
        }
        self.stop_lost_memo = Some(false);
        let v = match self.peptide_alleles() {
            Some((ref_pep, alt_pep)) if !alt_pep.contains(&b'X') => {
                !alt_pep.contains(&b'*') && ref_pep.contains(&b'*')
            }
            _ => {
                self.ins_del_stop_altered()
                    || self.removed_stop_codon() == Some(GivenStopCodon::Lost)
            }
        };
        self.stop_lost_memo = Some(v);
        v
    }

    fn ref_eq_alt_sequence(&mut self) -> bool {
        let Some((ref_pep, alt_pep)) = self.peptide_alleles() else {
            return false;
        };
        let (Some(tl_start), Some(tl_end)) = (self.span.tl_start, self.span.tl_end) else {
            return false;
        };
        if ref_pep == b"X" && alt_pep == b"X" {
            return false;
        }
        let ref_seq: &[u8] = &self.pep_seq;
        // The tests of VariationEffect.pm:1344 and the first of :1357 see no
        // stop codon in the reference: each holds only where one stands at
        // the coding end.
        let stop_codon_stands = self.stop_codon_stands();
        if stop_codon_stands && tl_start > ref_seq.len() as i64 && alt_pep.first() == Some(&b'*') {
            return true;
        }
        let mut mut_seq = ref_seq.to_vec();
        perl_substr_assign(&mut mut_seq, tl_start - 1, tl_end - tl_start + 1, &alt_pep);
        let mut_substring = &mut_seq[..ref_seq.len().min(mut_seq.len())];
        // VariationEffect.pm:1356-1359: the protein is unchanged up to its old
        // length and the residue right after it is a stop, or the stop sits at
        // the same offset of both peptides.
        let final_stop = (ref_seq.len() < mut_seq.len()).then(|| &mut_seq[ref_seq.len()..]);
        if stop_codon_stands
            && ref_seq == mut_substring
            && final_stop.is_some_and(|t| t.first() == Some(&b'*'))
        {
            return true;
        }
        if let Some(ri) = ref_pep.iter().position(|&b| b == b'*') {
            if alt_pep.iter().position(|&b| b == b'*') == Some(ri) {
                return true;
            }
        }
        false
    }

    /// `stop_retained` (VariationEffect.pm:1290-1324): `ref_eq_alt_sequence`
    /// for an alternate peptide without an `X` (:1313); otherwise an insertion
    /// or deletion that overlaps the stop codon in genomic coordinates and is
    /// tested and found to leave a stop codon in place (:1318, the `_cil`
    /// helpers; [`Self::edited_stop_codon_kept_cil`] for the untested span),
    /// or an untested span holding the whole stop codon that the sequence the
    /// record gives shows kept ([`Self::removed_stop_codon`]).
    fn stop_retained(&mut self) -> bool {
        if self.partial_codon() {
            return false;
        }
        if self.stop_lost() {
            return false;
        }
        if let Some(v) = self.stop_retained_memo {
            return v;
        }
        self.stop_retained_memo = Some(false);
        let v = match self.peptide_alleles() {
            Some((_, alt_pep)) if !alt_pep.is_empty() && !alt_pep.contains(&b'X') => {
                self.ref_eq_alt_sequence()
            }
            _ => {
                ((self.insertion || self.deletion)
                    && self.overlaps_stop_codon_cil()
                    && self.edited_stop_codon_kept_cil())
                    || self.removed_stop_codon() == Some(GivenStopCodon::Retained)
            }
        };
        self.stop_retained_memo = Some(v);
        v
    }

    /// `frameshift` (VariationEffect.pm:1541-1567): the length arithmetic behind
    /// three guards, `partial_codon`, `stop_retained`, and a reference peptide
    /// that begins with the stop codon (:1554, an edit whose first affected
    /// residue is the stop leaves the reading frame intact), which holds for
    /// the stop no residue follows ([`Self::begins_at_the_stop_no_residue_follows`]).
    fn frameshift(&mut self) -> bool {
        if self.partial_codon() {
            return false;
        }
        if self.stop_retained() {
            return false;
        }
        if self.begins_at_the_stop_no_residue_follows() {
            return false;
        }
        self.frameshift_by_length()
    }

    /// Whether `_get_peptide_alleles` yields a reference peptide beginning `*`.
    fn ref_peptide_starts_with_stop(&mut self) -> bool {
        self.peptide_alleles()
            .is_some_and(|(ref_pep, _)| ref_pep.first() == Some(&b'*'))
    }

    /// The third guard of `frameshift`: the reference peptide begins with a
    /// stop after which the transcript's peptide names no residue, so no
    /// residue's frame can shift. Ensembl VEP tests the first letter alone
    /// (VariationEffect.pm:1554), which a stop the peptide holds and reads
    /// on past passes as well ([`Self::begins_at_a_stop_the_peptide_reads_past`]).
    fn begins_at_the_stop_no_residue_follows(&mut self) -> bool {
        self.ref_peptide_starts_with_stop() && !self.begins_at_a_stop_the_peptide_reads_past()
    }

    /// Whether the first residue of the allele's codons is a stop that the
    /// transcript's peptide holds as a stop with a residue named after it (a
    /// polymorphic pseudogene's). An allele whose length change is no
    /// multiple of three shifts the frame of those residues, whatever its
    /// first codon is. An insertion between two codons stands on no codon,
    /// its translation start lying past its translation end, and begins at
    /// no residue.
    fn begins_at_a_stop_the_peptide_reads_past(&self) -> bool {
        let Some(first) = self
            .span
            .tl_start
            .zip(self.span.tl_end)
            .and_then(|(start, end)| (start <= end).then_some(start))
        else {
            return false;
        };
        usize::try_from(first - 1)
            .is_ok_and(|at| self.pep_seq.get(at) == Some(&b'*') && self.peptide_reads_on_past(at))
    }

    /// Whether the transcript's peptide names a residue after the one at the
    /// zero-based place `at`: a stop it holds there is one it reads on past,
    /// and the stop no residue follows is not.
    fn peptide_reads_on_past(&self, at: usize) -> bool {
        at + 1 < self.pep_seq.len()
    }

    /// The length arithmetic of Ensembl VEP's `frameshift` (VariationEffect.pm:1556-1564)
    /// without its guards: the allele length against the CDS span the variant
    /// covers, `false` when either end of that span is undefined.
    fn frameshift_by_length(&self) -> bool {
        let (Some(cds_start), Some(cds_end)) = (self.span.cds_start, self.span.cds_end) else {
            return false;
        };
        let var_len = cds_end - cds_start + 1;
        let Some(allele_len) = self.seq_length else {
            return false;
        };
        (allele_len - var_len).abs() % 3 != 0
    }

    fn stop_gained(&mut self) -> bool {
        if let Some(v) = self.stop_gained_memo {
            return v;
        }
        self.stop_gained_memo = Some(false);
        // VariationEffect.pm:1221-1222: never beside stop_retained or stop_lost.
        if self.stop_retained() || self.stop_lost() {
            return false;
        }
        let Some((ref_pep, alt_pep)) = self.peptide_alleles() else {
            return false;
        };
        let v = alt_pep.contains(&b'*') && !ref_pep.contains(&b'*');
        self.stop_gained_memo = Some(v);
        v
    }

    fn inframe_insertion(&mut self) -> bool {
        let codons = self.codon_alleles();
        if self.start_lost() {
            return false;
        }
        let Some((ref_codon, alt_codon)) = codons else {
            return false;
        };
        if alt_codon.len() <= ref_codon.len() {
            return false;
        }
        let Some((ref_pep, mut alt_pep)) = self.peptide_alleles() else {
            return false;
        };
        if self.start_retained_variant() && alt_pep.ends_with(&ref_pep) {
            return false;
        }
        // VariationEffect.pm:1121: a stop codon still read as one lone stop
        // after the insertion (TAG to TAAG) is not an inframe insertion.
        if ref_pep == b"*" && alt_pep == b"*" {
            return false;
        }
        // `$alt_pep =~ s/\*.+/\*/`
        if let Some(i) = alt_pep.iter().position(|&b| b == b'*') {
            if i + 1 < alt_pep.len() {
                alt_pep.truncate(i + 1);
            }
        }
        alt_pep.starts_with(&ref_pep) || alt_pep.ends_with(&ref_pep)
    }

    fn inframe_deletion(&mut self) -> bool {
        if self.partial_codon() {
            return false;
        }
        let codons = self.codon_alleles();
        let peps = self.peptide_alleles();
        let Some((ref_codon, alt_codon)) = codons else {
            return false;
        };
        if alt_codon.len() >= ref_codon.len() {
            return false;
        }
        // VariationEffect.pm:1175: deleting the stop codon is stop_lost alone.
        if peps.is_some_and(|(ref_pep, _)| ref_pep == b"*") {
            return false;
        }
        if ref_codon.starts_with(&alt_codon) || ref_codon.ends_with(&alt_codon) {
            return true;
        }
        let (r, a) = trim_common_flanks(&ref_codon, &alt_codon);
        a.is_empty() && r.len().is_multiple_of(3)
    }

    fn missense_variant(&mut self) -> bool {
        let Some((ref_pep, alt_pep)) = self.peptide_alleles() else {
            return false;
        };
        if self.start_lost()
            || self.stop_lost()
            || self.stop_gained()
            || self.partial_codon()
            || self.stop_retained()
        {
            return false;
        }
        ref_pep != alt_pep && ref_pep.len() == alt_pep.len()
    }

    fn synonymous_variant(&mut self) -> bool {
        let Some((ref_pep, alt_pep)) = self.peptide_alleles() else {
            return false;
        };
        if ref_pep.is_empty() {
            return false;
        }
        alt_pep == ref_pep
            && !self.stop_retained()
            && !self.start_retained_variant()
            && !alt_pep.contains(&b'X')
            && !ref_pep.contains(&b'X')
    }

    fn protein_altering_variant(&mut self) -> bool {
        let Some((ref_pep, alt_pep)) = self.peptide_alleles() else {
            return false;
        };
        if alt_pep.len() == ref_pep.len() {
            return false;
        }
        if ref_pep.first() == Some(&b'*') || alt_pep.first() == Some(&b'*') {
            return false;
        }
        if alt_pep.starts_with(&ref_pep) || alt_pep.ends_with(&ref_pep) {
            return false;
        }
        if self.inframe_deletion() {
            return false;
        }
        if self.start_lost() {
            return false;
        }
        if self.frameshift() {
            return false;
        }
        true
    }

    fn within_cds(&self) -> bool {
        for c in &self.span.cds_coords {
            if let PerlMapSeg::Coord { start, end } = *c {
                if end > 0 && start <= self.cds.len() as i64 {
                    return true;
                }
            }
        }
        if let Some((crs, cre)) = self.coding_region {
            let in_frameshift_intron = self
                .frameshift_introns
                .iter()
                .any(|&(is, ie)| perl_overlap(self.vf_start, self.vf_end, is, ie));
            if in_frameshift_intron {
                return perl_overlap(self.vf_start, self.vf_end, crs, cre);
            }
        }
        false
    }

    fn coding_unknown(&mut self) -> bool {
        if !self.within_cds() {
            return false;
        }
        let alt_peptide = self.peptide(false);
        let peps = self.peptide_alleles();
        // The peptide of an allele that holds unknown bases is read only
        // past a stop codon its known bases settle, so the `X` of a codon
        // that holds an unknown base stands past the protein's new end and
        // leaves the change known. The `X` of bases past the last whole
        // codon counts as for any allele.
        let partial_codon = self
            .codon(false)
            .is_some_and(|codon| !codon.len().is_multiple_of(3));
        let known_bases = is_perl_unambiguous_dna(self.alt_allele);
        let unknown_residue = |p: &Vec<u8>| {
            if known_bases {
                p.contains(&b'X')
            } else {
                partial_codon && p.last() == Some(&b'X')
            }
        };
        let cond = alt_peptide.is_none()
            || peps.is_none()
            || alt_peptide.as_ref().is_some_and(unknown_residue)
            || peps.as_ref().is_some_and(|(r, _)| r.contains(&b'X'));
        if !cond {
            return false;
        }
        !(self.frameshift()
            || self.inframe_deletion()
            || self.protein_altering_variant()
            || self.start_retained_variant()
            || self.start_lost()
            || self.stop_retained()
            || self.stop_lost())
    }

    /// Whether an allele longer than its genomic span is as long as the CDS
    /// span it stands on. Where an edited model holds bases the genome lacks,
    /// the genomic bases on either side map to CDS positions that are not
    /// neighbours, so the span of an insertion between them is the model's
    /// bases (`perl_map_insert`) and the alternate coding sequence carries the
    /// allele in their place: such an allele replaces them base for base.
    /// Ensembl VEP reads the length class from the genomic span alone
    /// (`_bvf_preds`), so `missense_variant`, which excludes
    /// `increase_length` (Utils/Constants.pm), is never asked of the row, and
    /// one whose residues change takes the default `intergenic_variant`
    /// (BaseVariationFeatureOverlapAllele.pm:285).
    fn stands_on_as_many_model_bases(&self) -> bool {
        self.insertion
            && matches!(
                (self.span.cds_start, self.span.cds_end, self.seq_length),
                (Some(start), Some(end), Some(allele_len)) if allele_len == end - start + 1
            )
    }

    /// Whether an allele longer than its genomic span is shorter than the CDS
    /// span it stands on: among the bases an edited model holds and the genome
    /// lacks ([`Self::stands_on_as_many_model_bases`]) it brings fewer bases than
    /// it replaces, so the model loses bases where the genome gains them. Ensembl
    /// VEP asks `inframe_deletion` only of an allele shorter than its genomic span
    /// (`deletion`, Utils/Constants.pm) and withholds `protein_altering_variant`
    /// where that predicate holds (VariationEffect.pm:388): the row reads
    /// `protein_altering_variant` where the reference codons, read with the span
    /// cut out (TranscriptVariationAllele.pm:872-873), fail it, and otherwise the
    /// default `intergenic_variant` (BaseVariationFeatureOverlapAllele.pm:285).
    fn stands_on_more_model_bases(&self) -> bool {
        self.insertion
            && matches!(
                (self.span.cds_start, self.span.cds_end, self.seq_length),
                (Some(start), Some(end), Some(allele_len)) if allele_len < end - start + 1
            )
    }

    /// The coding terms, evaluated in `%OVERLAP_CONSEQUENCES` rank order with each
    /// predicate's `include` filter applied; `missense_variant` is also asked of
    /// an allele that [`Self::stands_on_as_many_model_bases`], and
    /// `inframe_deletion` of one that [`Self::stands_on_more_model_bases`].
    fn terms(&mut self) -> Vec<Consequence> {
        let mut out = Vec::new();
        if self.stop_gained() {
            out.push(Consequence::StopGained);
        }
        if !self.snp && self.frameshift() {
            out.push(Consequence::FrameshiftVariant);
        }
        if self.stop_lost() {
            out.push(Consequence::StopLost);
        }
        if self.start_lost() {
            out.push(Consequence::StartLost);
        }
        if self.insertion && self.inframe_insertion() {
            out.push(Consequence::InframeInsertion);
        }
        if (self.deletion || self.stands_on_more_model_bases()) && self.inframe_deletion() {
            out.push(Consequence::InframeDeletion);
        }
        if (self.snp || self.stands_on_as_many_model_bases()) && self.missense_variant() {
            out.push(Consequence::MissenseVariant);
        }
        if self.protein_altering_variant() {
            out.push(Consequence::ProteinAlteringVariant);
        }
        if self.partial_codon() {
            out.push(Consequence::IncompleteTerminalCodonVariant);
        }
        if self.start_retained_variant() {
            out.push(Consequence::StartRetainedVariant);
        }
        if self.stop_retained() {
            out.push(Consequence::StopRetainedVariant);
        }
        if self.synonymous_variant() {
            out.push(Consequence::SynonymousVariant);
        }
        if self.coding_unknown() {
            out.push(Consequence::CodingSequenceVariant);
        }
        out
    }
}

/// Whether the row of an alternate allele that holds unknown bases reads a
/// peptide for it (`PerlCodingEval::reading_unknown_bases`), on the span the
/// row is annotated at.
pub(crate) fn unknown_bases_allele_is_translated(
    variant: &InputVariant,
    transcript: &Transcript,
) -> bool {
    PerlCodingEval::reading_unknown_bases(variant, transcript, variant.start, variant.end, None)
        .is_some()
}

/// Ensembl VEP's `codon` / `peptide` strings for both alleles of a coding variant, for
/// tests: `(ref_codon, alt_codon, ref_pep, alt_pep)` with `-` for
/// an empty window and `None` where Ensembl VEP has `undef`.
#[allow(clippy::type_complexity)]
#[cfg(test)]
pub(crate) fn perl_codon_peptides(
    variant: &InputVariant,
    transcript: &Transcript,
    span_start: u64,
    span_end: u64,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<(
    Option<Vec<u8>>,
    Option<Vec<u8>>,
    Option<Vec<u8>>,
    Option<Vec<u8>>,
)> {
    let mut ev = PerlCodingEval::new(variant, transcript, span_start, span_end, fasta)?;
    Some((
        ev.codon(true),
        ev.codon(false),
        ev.peptide(true),
        ev.peptide(false),
    ))
}

/// Ensembl VEP `TranscriptVariationAllele::hgvs_protein` (TranscriptVariationAllele.pm
/// 1657), the text after `p.`.
///
/// `variant` is the allele as the input placed it; `shifted`, when the caller
/// moved an insertion or deletion to its most 3' position (Ensembl VEP's
/// `_return_3prime(1)`), is the shifted allele with its span. Ensembl VEP reads the
/// peptides, the translation coordinates and the alternate CDS from the shifted
/// allele. For an insertion whose run of placements is longer than
/// `perform_shift` walks (303-306), `shifted` is the allele at the run's last
/// placement, the span its HGVSc is described at (`hgvs::hgvsc_shift`), where
/// Ensembl VEP reads the span its walk reaches (1711, 1715-1721).
/// The `coding` pre-consequence predicate (1731) and the cached
/// predicates `stop_lost`, `start_lost`, `partial_codon` and `stop_retained`
/// were filled while the consequences were computed: on the unshifted allele
/// without `--shift_3prime 1`, so those verdicts are taken from `variant` here
/// whatever the shift; on the shifted allele under the flag
/// (`predicates_on_shifted`), where `hgvs_transcript` also clears the cache once
/// the shift moved the allele (1469), so every verdict is the shifted allele's,
/// and a deletion shifted onto the stop codon whose alternate peptide is the
/// partial-codon `X` reads both peptides as `X` (`_get_hgvs_protein_type`,
/// 2060-2061) and prints `Ter<pos>=`. `frameshift` (VariationEffect.pm 1541) is
/// not cached: its guards read the cache (`partial_codon`, `stop_retained`, and
/// a reference peptide starting with the stop, VariationEffect.pm 1554, so a
/// frameshift that begins in the stop codon no residue follows is typed by its
/// peptides; `PerlCodingEval::begins_at_the_stop_no_residue_follows`), but its
/// length arithmetic runs on the CDS span the transcript variation carries at
/// that point, which is the shifted span, so an indel that shifts fully into the
/// CDS is a frameshift there even where the annotated allele straddles an exon
/// boundary.
///
/// An allele annotated in an intron inside the coding region fails that `coding`
/// test, so where the shift carries it wholly into the coding sequence Ensembl
/// VEP writes the shifted HGVSc and no HGVSp, except on a transcript with a
/// frameshift intron, whose exons the test reads 12 bases wider
/// (BaseTranscriptVariation.pm 861-866). The shifted allele has a protein
/// notation, the one written for the same allele given at the shifted position,
/// so every verdict is then the shifted allele's, as under the flag. An
/// insertion the shift leaves between the intron's last base and the first base
/// of the next exon (`inserted_at_the_end_of_an_intron`) is a change at the
/// junction, for which HGVS predicts no protein consequence: it keeps the
/// verdicts of the allele as annotated, and with them the empty HGVSp. On a
/// transcript whose coding sequence does not end on its stop codon
/// (`coding_sequence_ends_on_its_stop_codon`) the verdicts are the shifted
/// allele's where it changes no residue before the stop codon
/// (`changes_no_residue_before_the_stop_codon`), where the codons the
/// transcript holds give its notation
/// (`held_codons_give_the_protein_notation`), or where the allele as annotated
/// passes Ensembl VEP's own test on the wider exons of a frameshift-intron
/// transcript; any other allele there keeps the verdicts of the allele as
/// annotated.
///
/// An allele annotated on the transcript before its coding region, in the 5'
/// UTR or in an intron of it, reads `utr` (`annotated_before_the_coding_start`)
/// and fails the same `coding` test, so where the shift carries it into the
/// coding sequence Ensembl VEP writes the shifted HGVSc and no HGVSp there
/// too. Its verdicts are the shifted allele's under the same conditions, and
/// its notation the one written for the same allele given at the shifted
/// position. The shifted allele is put to the `coding` test like any other:
/// one the shift leaves before the coding start has no notation, and neither
/// has one whose shifted span has no translation start.
///
/// An allele that has no start codon test where the input placed it takes
/// `start_lost`, and the test of a coding sequence kept whole behind the
/// allele, from the shifted allele: an insertion at the last base before the
/// coding start does not reach the start codon's cDNA positions, and a
/// deletion that begins before the transcript's first base has no cDNA start
/// (VariationEffect.pm 977-981). The notation describes the shifted allele,
/// and Ensembl VEP writes the lost start for the same allele given at the
/// shifted position. Without `--shift_3prime 1` it keeps the verdict cached
/// for the allele as placed (VariationEffect.pm 855-862;
/// TranscriptVariationAllele.pm 1796 has dropped the shift before 2161 asks)
/// beside the frameshift of the shifted allele, and names a frameshift at
/// residue 1 with the unknown count, or a deletion of residues from the
/// first.
///
/// On a span whose stop codon edit test is not made, Ensembl VEP's
/// `stop_retained` holds for want of the test (VariationEffect.pm:1508 read at
/// :1318) and guards `frameshift` (:1549). That guard stands here where the
/// allele, as the notation places it, runs across exons: the notation is
/// Ensembl VEP's. It is dropped where the placed allele lies inside one exon,
/// whose alternate CDS shows the altered stop codon.
///
/// `None` where Ensembl VEP returns `undef`: the allele the verdicts are read from does
/// not overlap the coding sequence, the shift leaves the allele outside the
/// transcript, the (shifted) span has no translation start or end, or its
/// reference peptide is undefined. `None` as well for an insertion the input
/// writes at the end of an intron, in a form
/// `intron_end_insertion_has_no_protein_notation` names, and for one the 3'
/// shift leaves at the start of an intron
/// (`insertion_left_at_the_start_of_an_intron`).
///
/// Two more intended divergences. The translations inside follow the transcript's
/// codon table (`hgvsp_fs_peptides`), where BioPerl's argument-less
/// `translate()` of the alternate CDS at 2333, 2492 and 2555 and of the
/// reference CDS at 2450 reads table 1 whatever the transcript's table, against
/// a reference peptide of that table (Ensembl `Transcript::translate`). And
/// Ensembl VEP prints `Met1?` whenever its `start_lost` predicate holds (2161),
/// including the start co-emission pairs on which this engine keeps
/// `start_retained_variant` and drops `start_lost`; the port fires that
/// short-circuit only when it emits `start_lost` itself. On such a pair
/// whose insertion or deletion leaves the whole coding sequence standing
/// behind a start codon the 5' UTR's bases form again
/// (`PerlCodingEval::keeps_the_coding_sequence_whole_behind_the_allele`), the
/// notation is the first residue, unchanged: `Met1=`.
pub fn perl_hgvs_protein(
    variant: &InputVariant,
    shifted: Option<(&InputVariant, u64, u64)>,
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
    predicates_on_shifted: bool,
) -> Option<String> {
    hgvs_protein(
        variant,
        shifted,
        transcript,
        fasta,
        predicates_on_shifted,
        false,
    )
}

/// [`perl_hgvs_protein`] of an insertion that is described where its record
/// is written although the genome's 3' shift moves it: one standing on bases
/// an edited model holds and the genome lacks. The notation is read from the
/// row's own span. Where those bases open the coding sequence that span has
/// no translation start ([`span_on_opening_coding_bases`]), and the notation
/// reads the coding span of the bases themselves, under the verdicts of the
/// row. Ensembl VEP reads the shifted span there (1711, 1724-1733), past
/// the model's bases, and adds the residues of the whole allele after them.
pub(crate) fn hgvs_protein_where_written(
    variant: &InputVariant,
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<String> {
    hgvs_protein(variant, None, transcript, fasta, false, true)
}

/// The span `perl_span` gives an insertion between two exonic bases, read
/// on the bases the model holds between them where they open the coding
/// sequence: the first of them is the model's first coding base. `map_insert`
/// ends the insertion's first half on the base before them (Mapper.pm:522-542),
/// a base of the 5' UTR, which `genomic2cds` makes a gap
/// (TranscriptMapper.pm:441-443): the span has no CDS start and no
/// translation start, although every base it stands on is coding. The span
/// returned starts on the first coding base. `None` for any other span.
fn span_on_opening_coding_bases(transcript: &Transcript, span: &PerlSpan) -> Option<PerlSpan> {
    let mapper = transcript.vefc.as_ref()?.mapper.as_ref()?;
    let coding_start = mapper.cdna_coding_start as i64;
    let [PerlMapSeg::Gap { .. }, last @ PerlMapSeg::Coord { .. }] = span.cds_coords.as_slice()
    else {
        return None;
    };
    let stands_on_opening_bases = span.cdna_start == Some(coding_start)
        && span.cdna_end.is_some_and(|end| end >= coding_start);
    if !stands_on_opening_bases {
        return None;
    }
    let phase = i64::from(mapper.start_phase.max(0));
    Some(PerlSpan {
        cds_start: Some(1 + phase),
        tl_start: Some((3 + phase).div_euclid(3)),
        cds_coords: vec![PerlMapSeg::Coord { start: 1, end: 0 }, *last],
        ..span.clone()
    })
}

/// [`perl_hgvs_protein`], and [`hgvs_protein_where_written`] with
/// `where_written`.
fn hgvs_protein(
    variant: &InputVariant,
    shifted: Option<(&InputVariant, u64, u64)>,
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
    predicates_on_shifted: bool,
    where_written: bool,
) -> Option<String> {
    let unknown_bases =
        !variant.alt_allele().is_empty() && !is_perl_unambiguous_dna(variant.alt_allele());
    // The span the notation's predicates are read at: the shifted span under
    // `--shift_3prime 1` and for an allele annotated in a coding intron or
    // before the coding start (`annotated_before_the_coding_start`) that the
    // 3' shift carries into the coding sequence, the span as written
    // otherwise. An allele of unknown bases keeps the span as written where
    // the shift passes one of them (`shift_passes_an_unknown_base`): no
    // completion of the allele is carried there.
    let predicates_on_shifted = predicates_on_shifted
        || shifted.is_some_and(|(allele, start, end)| {
            let in_a_coding_intron =
                annotated_in_coding_intron(transcript, variant.start, variant.end);
            (in_a_coding_intron
                || annotated_before_the_coding_start(transcript, variant.start, variant.end))
                && !inserted_at_the_end_of_an_intron(allele, transcript)
                && !(unknown_bases
                    && shift_passes_an_unknown_base(variant, start, transcript.strand))
                && (coding_sequence_ends_on_its_stop_codon(transcript)
                    || (in_a_coding_intron
                        && crate::consequences::intronic_allele_reads_as_coding(
                            transcript, variant,
                        ))
                    || changes_no_residue_before_the_stop_codon(
                        allele, start, end, transcript, fasta,
                    )
                    || held_codons_give_the_protein_notation(allele, start, end, transcript, fasta))
        });
    let (predicate_variant, predicate_start, predicate_end) = match shifted {
        Some((shifted_variant, start, end)) if predicates_on_shifted => {
            (shifted_variant, start, end)
        }
        _ => (variant, variant.start, variant.end),
    };
    // An allele that holds unknown bases has a notation only where it reads
    // its peptide at that span (`PerlCodingEval::reading_unknown_bases`): its
    // known bases settle a stop codon there before the first unknown base. A
    // shift that an unknown base of the allele stops leaves that base first
    // among the allele's codons, so no stop codon stands before it. Each
    // evaluation below reads the allele the same way at its own span.
    if unknown_bases {
        let row = PerlCodingEval::reading_unknown_bases(
            predicate_variant,
            transcript,
            predicate_start,
            predicate_end,
            fasta,
        )?;
        // At an exon's edge the notation turns on the inserted bases beside
        // the edge (`intron_end_insertion`,
        // `insertion_left_at_the_start_of_an_intron`), which are not read
        // for an allele that holds unknown bases.
        if !row.inside_one_exon() {
            return None;
        }
        // The engines leave the bases a record of two or more alternate
        // alleles shares with its reference allele uncut (Parser.pm:964-1019
        // cuts a record of one alone), and a notation is read on the allele
        // with those bases clipped, at a span this evaluation does not read.
        if alleles_share_an_end_base(row.ref_allele, row.alt_allele) {
            return None;
        }
    }
    let evaluation = |allele, start, end| {
        if unknown_bases {
            PerlCodingEval::reading_unknown_bases(allele, transcript, start, end, fasta)
        } else {
            PerlCodingEval::new(allele, transcript, start, end, fasta)
        }
    };
    let mut ev = evaluation(predicate_variant, predicate_start, predicate_end)?;
    if !ev.coding_pred(transcript) {
        return None;
    }
    // Tested before every description built below, the stop codon's included:
    // an insertion before an exon that opens inside the stop codon is described
    // at the end of the intron, as before any other coding exon.
    if intron_end_insertion(variant, transcript, fasta)
        .is_some_and(intron_end_insertion_has_no_protein_notation)
    {
        return None;
    }
    if insertion_left_at_the_start_of_an_intron(variant, shifted, transcript, fasta) {
        return None;
    }
    let frameshift_guard =
        ev.partial_codon() || ev.stop_retained() || ev.begins_at_the_stop_no_residue_follows();
    let begins_at_a_held_stop = ev.begins_at_a_stop_the_peptide_reads_past();
    let untested_stop_kept = ev.untested_stop_codon_kept_in_ensembl_vep();
    let stop_lost = ev.stop_lost();
    let cut_short_names_its_first_residue = ev.cut_short_frameshift_names_its_first_residue();
    let mut start_lost = ev.start_lost() && !ev.start_retained_variant();
    // Read on the pair alone, the allele on which Ensembl VEP's `start_lost`
    // holds beside `start_retained_variant`.
    let mut coding_sequence_kept_whole =
        ev.start_lost() && ev.keeps_the_coding_sequence_whole_behind_the_allele();
    let start_codon_tested_where_written = ev.overlaps_start_codon();
    if let Some((shifted_variant, span_start, span_end)) = shifted {
        if !predicates_on_shifted {
            // The shift follows the genome, and a place it reaches outside the
            // transcript has no translation coordinates (1730-1737): so too the
            // place after the last genomic base of a transcript whose coding
            // sequence runs past that base, on a stop codon an edit completes.
            if !perl_overlap(
                span_start as i64,
                span_end as i64,
                transcript.start as i64,
                transcript.end as i64,
            ) {
                return None;
            }
            ev = evaluation(shifted_variant, span_start, span_end)?;
            // The start codon tests of an allele that has none where the
            // input placed it are the shifted allele's; Ensembl VEP keeps the
            // verdict cached for the placed allele (VariationEffect.pm:855-862).
            if !start_codon_tested_where_written {
                start_lost = ev.start_lost() && !ev.start_retained_variant();
                coding_sequence_kept_whole =
                    ev.start_lost() && ev.keeps_the_coding_sequence_whole_behind_the_allele();
            }
        }
    }
    if let Some(span) = where_written
        .then(|| span_on_opening_coding_bases(transcript, &ev.span))
        .flatten()
    {
        ev = evaluation(variant, variant.start, variant.end)?;
        ev.span = span;
    }
    let preds = HgvsPredicates {
        frameshift: !frameshift_guard
            && !(untested_stop_kept && !ev.inside_one_exon())
            && ev.frameshift_by_length(),
        stop_lost,
        start_lost,
        begins_at_a_held_stop,
        cut_short_names_its_first_residue,
    };
    if unknown_bases && !ev.inside_one_exon() {
        return None;
    }
    let tl_start = truthy(ev.span.tl_start)?;
    let tl_end = truthy(ev.span.tl_end)?;
    let alt_pep = ev.peptide(false);
    // An allele that holds unknown bases has a protein notation only where
    // its peptide is defined at the span the notation is read at: where its
    // known bases settle a stop codon there, a gained one
    // (`gains_a_stop_codon_before_unknown_bases`) or the transcript's own,
    // kept (`keeps_a_stop_codon_before_unknown_bases`). Ensembl VEP writes
    // none for any such allele (TranscriptVariationAllele.pm:1683-1684).
    if alt_pep.is_none() && !is_perl_unambiguous_dna(ev.alt_allele) {
        return None;
    }
    // A frameshift's notation names the first changed residue and the place
    // of the new stop and holds no inserted residue against the reference.
    if unknown_bases && !preds.frameshift {
        let reach = alt_pep.as_ref().map_or(0, Vec::len);
        let reference = perl_translate(ev.cds, ev.table);
        if unknown_residue_within_reach(
            &reference,
            tl_start.min(tl_end),
            tl_start.max(tl_end),
            reach,
        ) {
            return None;
        }
    }
    let ref_pep = ev.peptide(true).filter(|p| !p.is_empty())?;
    // The alternate coding sequence holds no base of the 5' UTR, so the frame
    // it shows from residue 1 is one the edited transcript does not read:
    // the start codon the 5' UTR's bases form again opens the coding
    // sequence, whole, and the first residue stands unchanged.
    if coding_sequence_kept_whole {
        if let Some(&first) = ev.pep_seq.first() {
            let description = format!("{}1=", String::from_utf8_lossy(&seq3(&[first])));
            return (!names_a_stop_the_reference_lacks(&ev, &description)).then_some(description);
        }
    }
    let mut n = HgvsProteinNotation {
        ref_pep: Some(ref_pep),
        alt_pep,
        start: tl_start,
        end: tl_end,
        kind: HgvspKind::Unset,
        original_ref: Vec::new(),
        preseq: Vec::new(),
        cut_short: false,
    };
    if n.alt_pep.is_some() && n.alt_pep != n.ref_pep {
        hgvsp_clip_alleles(&mut n);
    }
    hgvsp_protein_type(&ev, &preds, &mut n);
    if n.kind == HgvspKind::Unset {
        return None;
    }
    let description = match hgvsp_peptides(&ev, &preds, &mut n) {
        Some(()) => hgvsp_format(&ev, &preds, &n),
        None => hgvsp_of_an_insertion_no_residue_follows(&ev, &n)?,
    };
    (!names_a_stop_the_reference_lacks(&ev, &description)).then_some(description)
}

/// The forms of an insertion that the input writes between an intron's last
/// base and the first base of the next exon, told apart by the inserted bases
/// as the transcript's strand reads them and the reference bases before the
/// gap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum IntronEndInsertion {
    /// One inserted base, equal to the intron's last base.
    CopyOfTheLastBase,
    /// Two or more inserted bases, equal to the bases of the intron that end
    /// on its last base.
    DuplicationOfTheLastBases,
    /// Bases that are no copy of the bases before the gap and cannot be the
    /// intron's new end: one base, or bases that do not end in the two bases
    /// that end an intron.
    CannotEndTheIntron,
    /// Bases that are no copy of the bases before the gap and end in the two
    /// bases that end an intron.
    EndsInAnAcceptorPair,
}

/// The insertion an allele makes, as the coding notation reads it: the inserted
/// bases on the transcript's strand and the two positions they lie between,
/// the lower first. `variant` is an allele on the forward strand. A pure
/// insertion gives its own bases and gap. An allele the record holds with
/// reference bases, as a record that carries another allele beside an
/// insertion holds it with the base they share, is clipped as `_clip_alleles`
/// (TranscriptVariationAllele.pm 2188) clips it for `hgvs_transcript`: the
/// shared bases from the front and then from the back, on the transcript's
/// strand. It is an insertion where no reference base remains, lying after
/// the bases clipped from the front. `None` for every other allele.
fn clipped_insertion(variant: &InputVariant, reverse: bool) -> Option<(Vec<u8>, u64, u64)> {
    let on_transcript = |allele: &[u8]| {
        let mut bases = allele.to_ascii_uppercase();
        bases.retain(|&base| base != b'-');
        if reverse {
            reverse_complement(&bases)
        } else {
            bases
        }
    };
    let reference = on_transcript(&variant.ref_allele);
    let alternate = on_transcript(variant.alt_allele());
    let front = reference
        .iter()
        .zip(&alternate)
        .take_while(|(reference, alternate)| reference == alternate)
        .count();
    let (reference, alternate) = (&reference[front..], &alternate[front..]);
    let back = reference
        .iter()
        .rev()
        .zip(alternate.iter().rev())
        .take_while(|(reference, alternate)| reference == alternate)
        .count();
    if back != reference.len() || back == alternate.len() {
        return None;
    }
    let inserted = alternate[..alternate.len() - back].to_vec();
    // The record spans `start..=end`, no base for a pure insertion, and the
    // transcript reads the span from `end` down on the reverse strand.
    let front = front as u64;
    let (low, high) = if reverse {
        (variant.end.checked_sub(front)?, variant.end - front + 1)
    } else {
        (
            (variant.start + front).checked_sub(1)?,
            variant.start + front,
        )
    };
    Some((inserted, low, high))
}

/// The form of an insertion written at the end of an intron: `variant`, the
/// allele as the input placed it, on the forward strand, lies between an
/// intron's last base and the first base of the next exon in the transcript's
/// direction, that base is inside the coding region, and the 3' rule of the
/// HGVS nomenclature leaves the insertion in that gap, its first base
/// differing from the exon's first base. The insertion is read as the coding
/// notation reads the allele ([`clipped_insertion`]), so a record that holds
/// it with a base it shares with the reference has the form of the same
/// insertion written alone.
///
/// `None` for every other allele: one the input places anywhere else, wherever
/// the 3' shift carries it; one whose first base equals the exon's first base,
/// which the 3' rule describes inside the exon; one at the end of an intron of
/// at most 13 bases, which Ensembl VEP reads with the exons beside it
/// (`_frameshift`, `BaseTranscriptVariation::_create_intron_trees`); a copy of
/// bases that reach back past the intron's first base; and any allele where no
/// reference gives the bases beside the gap.
fn intron_end_insertion(
    variant: &InputVariant,
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<IntronEndInsertion> {
    let fasta = fasta?;
    let reverse = transcript.strand == Strand::Reverse;
    let (inserted, gap_low, gap_high) = clipped_insertion(variant, reverse)?;
    // The gap lies between `gap_low` and `gap_high`, and the intron's last base
    // is the one of the two the transcript reads first.
    let intron = transcript.introns.iter().find(|intron| {
        intron.start <= intron.end
            && if reverse {
                intron.start == gap_high
            } else {
                intron.end == gap_low
            }
    })?;
    if intron.end - intron.start <= 12 {
        return None;
    }
    let exon_first = if reverse { gap_low } else { gap_high };
    let opens_an_exon = transcript.exons.iter().any(|exon| {
        if reverse {
            exon.end == exon_first
        } else {
            exon.start == exon_first
        }
    });
    let coding = matches!(
        (transcript.coding_region_start, transcript.coding_region_end),
        (Some(start), Some(end)) if start <= exon_first && exon_first <= end
    );
    if !opens_an_exon || !coding {
        return None;
    }

    let on_transcript = |forward: Vec<u8>| {
        if reverse {
            reverse_complement(&forward)
        } else {
            forward
        }
    };
    if !is_perl_unambiguous_dna(&inserted) {
        return None;
    }
    let exon_first_base = on_transcript(fasta.sequence(&variant.chr, exon_first, exon_first)?);
    if inserted.first() == exon_first_base.first() {
        return None;
    }
    let length = inserted.len() as u64;
    let (before_start, before_end) = if reverse {
        (intron.start, intron.start.checked_add(length - 1)?)
    } else {
        ((intron.end + 1).checked_sub(length)?, intron.end)
    };
    let before = on_transcript(fasta.sequence(&variant.chr, before_start, before_end)?);
    if inserted == before {
        return match length {
            _ if length > intron.end - intron.start + 1 => None,
            1 => Some(IntronEndInsertion::CopyOfTheLastBase),
            _ => Some(IntronEndInsertion::DuplicationOfTheLastBases),
        };
    }
    // The two bases that end an intron of the major kind and of the minor kind.
    let ends_an_intron = inserted.ends_with(b"AG") || inserted.ends_with(b"AC");
    Some(if length == 1 || !ends_an_intron {
        IntronEndInsertion::CannotEndTheIntron
    } else {
        IntronEndInsertion::EndsInAnAcceptorPair
    })
}

/// Whether an insertion written at the end of an intron carries no protein
/// notation, where the codon window would read the inserted bases as the
/// exon's first. Each arm is one form of [`IntronEndInsertion`].
fn intron_end_insertion_has_no_protein_notation(form: IntronEndInsertion) -> bool {
    match form {
        // A duplication of the intron's last position: the description names
        // no coding position, and for a change of an intron's last two
        // positions the HGVS nomenclature gives no predicted protein
        // consequence, only the unknown-effect form `p.?`.
        IntronEndInsertion::CopyOfTheLastBase => true,

        // A duplication of two or more intron positions: the description
        // names no coding position, the two bases before the exon stand as
        // they were, and a change that brings a second pair of an intron's
        // end bases beside the first may affect splicing, for which the HGVS
        // nomenclature gives no predicted protein consequence.
        IntronEndInsertion::DuplicationOfTheLastBases => true,

        // An insertion between the intron's last position and the exon's
        // first that leaves the intron without its two end bases: read as the
        // intron's end it changes the intron's last two positions, read as
        // the exon's start it changes the exon's first base, and the HGVS
        // nomenclature gives neither reading a predicted protein consequence.
        IntronEndInsertion::CannotEndTheIntron => true,

        // An insertion between the same two positions whose last two bases
        // can stand as the intron's end: read so, it puts a second pair of
        // end bases before the exon and moves the intron's own upstream;
        // read as the exon's start, it changes the exon's first bases. The
        // HGVS nomenclature gives neither reading a predicted protein
        // consequence.
        IntronEndInsertion::EndsInAnAcceptorPair => true,
    }
}

/// Whether an insertion lies, where the 3' shift of the HGVS notations leaves
/// it, between a coding exon's last base and the first base of the intron
/// after it, the 3' rule of the HGVS nomenclature carrying it no further, and
/// its bases are no copy of the reference bases before that gap. Such an
/// insertion is described between the exon's last position and the intron's
/// first, and carries no protein notation: read as the intron's start it
/// changes the intron's first two positions, read as the exon's end it changes
/// the exon's last base, and the nomenclature gives neither reading a
/// predicted protein consequence.
///
/// `variant` is the allele as the input placed it and `shifted` the allele the
/// 3' shift moved; an allele written in the exon that the shift carries to the
/// gap is read at the gap. An allele the record holds with a base it shares
/// with the reference, which no shift moves, is read as the coding notation
/// clips it ([`clipped_insertion`]). A copy of the bases before the gap is a
/// duplication of coding positions and keeps its notation, also where it is
/// longer than its exon: the copy is told from the reference, base by base, and
/// not from the exon's bounds. `false` for an insertion whose first base
/// equals the intron's first base, which the 3' rule describes inside the
/// intron, and where no reference gives the bases beside the gap.
fn insertion_left_at_the_start_of_an_intron(
    variant: &InputVariant,
    shifted: Option<(&InputVariant, u64, u64)>,
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> bool {
    let Some(fasta) = fasta else {
        return false;
    };
    let reverse = transcript.strand == Strand::Reverse;
    let placed = shifted.map_or(variant, |(moved, ..)| moved);
    let Some((inserted, end, start)) = clipped_insertion(placed, reverse) else {
        return false;
    };
    // The gap lies between `end` and `start`, and the intron's first base is
    // the one of the two the transcript reads last.
    let opens_an_intron = transcript.introns.iter().any(|intron| {
        intron.start <= intron.end
            && if reverse {
                intron.end == end
            } else {
                intron.start == start
            }
    });
    let (exon_last, intron_first) = if reverse { (start, end) } else { (end, start) };
    let closes_an_exon = transcript.exons.iter().any(|exon| {
        if reverse {
            exon.start == exon_last
        } else {
            exon.end == exon_last
        }
    });
    let coding = matches!(
        (transcript.coding_region_start, transcript.coding_region_end),
        (Some(first), Some(last)) if first <= exon_last && exon_last <= last
    );
    if !opens_an_intron || !closes_an_exon || !coding {
        return false;
    }

    let on_transcript = |forward: Vec<u8>| {
        if reverse {
            reverse_complement(&forward)
        } else {
            forward
        }
    };
    if !is_perl_unambiguous_dna(&inserted) {
        return false;
    }
    let Some(intron_first_base) = fasta.sequence(&variant.chr, intron_first, intron_first) else {
        return false;
    };
    if inserted.first() == on_transcript(intron_first_base).first() {
        return false;
    }
    let length = inserted.len() as u64;
    let before = if reverse {
        exon_last
            .checked_add(length - 1)
            .and_then(|last| fasta.sequence(&variant.chr, exon_last, last))
    } else {
        (exon_last + 1)
            .checked_sub(length)
            .filter(|&first| first >= 1)
            .and_then(|first| fasta.sequence(&variant.chr, first, exon_last))
    };
    before.is_some_and(|before| inserted != on_transcript(before))
}

/// The three predicate verdicts `hgvs_protein` reads: `stop_lost` and
/// `start_lost` from Ensembl VEP's `_predicate_cache` (the allele as the input placed
/// it, or the shifted allele under `--shift_3prime 1` and for an allele the
/// input placed in an intron; `start_lost` from the shifted allele as well
/// where the placed allele has no start codon test), `frameshift`
/// recomputed on the shifted span behind the cached `partial_codon`,
/// `stop_retained` and leading-stop guards. `start_lost` already carries the
/// start co-emission gate of this engine. An allele the input placed before
/// the coding start is read as one it placed in an intron.
struct HgvsPredicates {
    frameshift: bool,
    stop_lost: bool,
    start_lost: bool,
    /// Whether the codons of the allele as the predicates read it begin at
    /// a stop the peptide holds and reads on past
    /// (`PerlCodingEval::begins_at_a_stop_the_peptide_reads_past`).
    begins_at_a_held_stop: bool,
    /// Whether a frameshift that the transcript's last base cuts short is
    /// named at its first residue
    /// (`PerlCodingEval::cut_short_frameshift_names_its_first_residue`).
    cut_short_names_its_first_residue: bool,
}

/// Ensembl VEP's `$hgvs_notation` hash. `ref_pep` and `alt_pep` are `None` where the
/// peptide is `undef`; `original_ref` and `preseq` are what `_clip_alleles`
/// records for the stop and duplication checks. `cut_short` is whether the
/// frameshift walk found that the alternate sequence ends one or two bases
/// into the codon of the residue the walk starts at, before the alternate
/// translation reaches it (`_get_fs_peptides`, 2342-2346).
struct HgvsProteinNotation {
    ref_pep: Option<Vec<u8>>,
    alt_pep: Option<Vec<u8>>,
    start: i64,
    end: i64,
    kind: HgvspKind,
    original_ref: Vec<u8>,
    preseq: Vec<u8>,
    cut_short: bool,
}

/// Ensembl VEP's `type` field: `""` is `Bare`, `"="` is `Eq`, `">"` is `Sub`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HgvspKind {
    Unset,
    Fs,
    Ins,
    Del,
    Sub,
    Delins,
    Dup,
    Eq,
    Bare,
}

impl HgvspKind {
    fn as_str(self) -> &'static str {
        match self {
            HgvspKind::Unset | HgvspKind::Bare => "",
            HgvspKind::Fs => "fs",
            HgvspKind::Ins => "ins",
            HgvspKind::Del => "del",
            HgvspKind::Sub => ">",
            HgvspKind::Delins => "delins",
            HgvspKind::Dup => "dup",
            HgvspKind::Eq => "=",
        }
    }
}

/// `Bio::SeqUtils::seq3`: three-letter codes, `Xaa` for anything unknown.
fn seq3(pep: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(pep.len() * 3);
    for &aa in pep {
        out.extend_from_slice(
            crate::hgvs::amino_acid_three_letter(aa.to_ascii_uppercase()).as_bytes(),
        );
    }
    out
}

/// Ensembl VEP's `$s =~ s/Xaa/Ter/g`.
fn replace_xaa_with_ter(pep: &mut [u8]) {
    let mut i = 0;
    while i + 3 <= pep.len() {
        if &pep[i..i + 3] == b"Xaa" {
            pep[i..i + 3].copy_from_slice(b"Ter");
        }
        i += 3;
    }
}

/// `_clip_alleles` (2188) with `numbering` `p`: trims the residues the alleles
/// share from the front, then from the back, and records `original_ref` and
/// `preseq`. Ensembl VEP trims into local copies of `start` and `end` and writes them
/// back only after both loops, so a leading stop on both sides, which returns
/// `Eq` at once, leaves the notation's positions and peptides as they were,
/// residues already trimmed included. The type re-set block runs as written
/// except its `dup` case, which applies only to nucleotide numbering.
fn hgvsp_clip_alleles(n: &mut HgvsProteinNotation) {
    let ref_pep: &[u8] = n.ref_pep.as_deref().unwrap_or_default();
    let alt_pep: &[u8] = n.alt_pep.as_deref().unwrap_or_default();
    n.original_ref = ref_pep.to_vec();
    let mut check_ref = ref_pep;
    let mut check_alt = alt_pep;
    let (mut start, mut end) = (n.start, n.end);
    let mut preseq = Vec::new();
    for _ in 0..ref_pep.len() {
        let next_ref = check_ref.first().copied();
        let next_alt = check_alt.first().copied();
        if next_ref == Some(b'*') && next_alt == Some(b'*') {
            n.kind = HgvspKind::Eq;
            return;
        }
        match next_ref {
            Some(shared) if next_alt == Some(shared) => {
                start += 1;
                check_ref = &check_ref[1..];
                check_alt = &check_alt[1..];
                preseq.push(shared);
            }
            _ => break,
        }
    }
    for _ in 0..check_ref.len() {
        if check_ref.last().is_some() && check_ref.last() == check_alt.last() {
            check_ref = &check_ref[..check_ref.len() - 1];
            check_alt = &check_alt[..check_alt.len() - 1];
            end -= 1;
        } else {
            break;
        }
    }
    n.start = start;
    n.end = end;
    let kind = if check_ref == check_alt {
        Some(HgvspKind::Eq)
    } else if check_ref != b"-" && check_ref.len() == 1 && check_alt.len() == 1 {
        Some(HgvspKind::Sub)
    } else if check_ref.is_empty() && !check_alt.is_empty() {
        Some(HgvspKind::Ins)
    } else if !check_ref.is_empty() && check_alt.is_empty() {
        Some(HgvspKind::Del)
    } else {
        None
    };
    let (check_ref, check_alt) = (check_ref.to_vec(), check_alt.to_vec());
    n.ref_pep = Some(check_ref);
    n.alt_pep = Some(check_alt);
    n.preseq = preseq;
    if let Some(kind) = kind {
        n.kind = kind;
    }
}

/// `_get_hgvs_protein_type` (2047): `fs` from the frameshift predicate; else the
/// first stop of each peptide becomes `X` and the lengths decide; without both
/// peptides the allele lengths less `-` decide.
fn hgvsp_protein_type(
    ev: &PerlCodingEval<'_>,
    preds: &HgvsPredicates,
    n: &mut HgvsProteinNotation,
) {
    if preds.frameshift {
        n.kind = HgvspKind::Fs;
        return;
    }
    if let (Some(r), Some(a)) = (n.ref_pep.as_mut(), n.alt_pep.as_mut()) {
        if let Some(p) = r.iter().position(|&b| b == b'*') {
            r[p] = b'X';
        }
        if let Some(p) = a.iter().position(|&b| b == b'*') {
            a[p] = b'X';
        }
        n.kind = if r.as_slice() == b"-" || r.is_empty() {
            HgvspKind::Ins
        } else if a.is_empty() || a.as_slice() == b"-" {
            HgvspKind::Del
        } else if r.len() == 1 && a.len() == 1 {
            HgvspKind::Sub
        } else if (!a.is_empty() && !r.is_empty() && a.len() != r.len())
            || (a.len() > 1 && r.len() > 1)
        {
            HgvspKind::Delins
        } else {
            HgvspKind::Sub
        };
        return;
    }
    // `_get_allele_length`: `s/\-//` strips the first dash only.
    let len_less_dash = |s: &[u8]| (s.len() - usize::from(s.contains(&b'-'))) as i64;
    let (ref_length, alt_length) = (len_less_dash(ev.ref_allele), len_less_dash(ev.alt_allele));
    if alt_length > 1 {
        n.kind = if n.start == n.end + 1 {
            HgvspKind::Ins
        } else if n.start != n.end {
            HgvspKind::Delins
        } else {
            HgvspKind::Sub
        };
    } else if ref_length > 1 {
        n.kind = HgvspKind::Del;
    }
}

/// `_get_hgvs_peptides` (2114) with three-letter conversion on. `None` where Ensembl VEP
/// returns `undef`: an insertion with no flanking residue to name.
fn hgvsp_peptides(
    ev: &PerlCodingEval<'_>,
    preds: &HgvsPredicates,
    n: &mut HgvsProteinNotation,
) -> Option<()> {
    match n.kind {
        HgvspKind::Fs => hgvsp_fs_peptides(ev, n, preds)?,
        HgvspKind::Ins => {
            hgvsp_post_var_shift(ev, n);
            if !n.alt_pep.as_deref().unwrap_or(b"").contains(&b'*') {
                hgvsp_check_duplication(ev, n);
            }
            if n.kind == HgvspKind::Dup {
                return Some(());
            }
            let min = n.start.min(n.end);
            n.ref_pep = Some(hgvsp_surrounding(ev, min, &n.original_ref, Some(2))?);
        }
        HgvspKind::Del => hgvsp_post_var_shift(ev, n),
        _ => {}
    }
    if let Some(r) = n.ref_pep.as_mut() {
        if r.as_slice() != b"-" {
            *r = seq3(r);
        }
    }
    if let Some(a) = n.alt_pep.as_mut() {
        if a.as_slice() != b"-" {
            *a = seq3(a);
        }
    }
    if n.alt_pep.as_deref() == Some(b"-") {
        n.alt_pep = Some(b"del".to_vec());
    }
    if preds.start_lost {
        n.alt_pep = Some(b"?".to_vec());
        n.kind = HgvspKind::Bare;
    } else if n.kind == HgvspKind::Del {
        let has_word = n
            .ref_pep
            .as_deref()
            .unwrap_or(b"")
            .iter()
            .any(|b| b.is_ascii_alphanumeric() || *b == b'_');
        if has_word {
            n.alt_pep = Some(b"del".to_vec());
        } else {
            hgvsp_del_peptides(ev, n)?;
        }
    } else if n.kind == HgvspKind::Fs {
        if let Some(r) = n.ref_pep.as_mut() {
            r.truncate(3);
        }
    }
    if let Some(r) = n.ref_pep.as_mut() {
        replace_xaa_with_ter(r);
    }
    if let Some(a) = n.alt_pep.as_mut() {
        replace_xaa_with_ter(a);
    }
    Some(())
}

/// `_get_fs_peptides` (2320): the first residue at which the translation of the
/// alternate CDS (3' UTR appended) differs from the reference peptide plus its
/// stop, from `translation_start`. `Del` when the alternate translation ends
/// before that position; `Eq` when both sides reach a stop together.
///
/// Such a `Del` is written as a frameshift named at its first residue where
/// the alternate sequence ends inside that residue's codon, no stop codon
/// stands at the coding end and Ensembl VEP reads one lost there
/// ([`PerlCodingEval::cut_short_frameshift_names_its_first_residue`]), and
/// as a deletion everywhere else.
///
/// `None` where the alternate translation ends before a residue differs: the
/// whole codons the transcript holds give the frameshift no changed residue
/// to be named by. Ensembl VEP's walk leaves by its bound there (2348) and
/// returns the pair it read last, two equal residues, at the position after
/// them (2351-2352, 2361, 2364), so it names a residue that does not change
/// by the name of the one before it (`docs/intended-divergences.md`). A lost
/// start keeps the walk's position, which its own description prints.
///
/// Ensembl VEP translates the alternate CDS here, and in the three paths below,
/// with `translate()` bare, the standard code, while its reference peptide
/// follows the transcript's codon table; a mitochondrial frameshift so compares
/// a table-2 reference with a table-1 alternate. The transcript's table is used
/// on both sides (`docs/intended-divergences.md`).
///
/// Under `start_lost` the description is the lost start, `Met1?` (2161-2165),
/// and residue 1 is compared as Ensembl VEP compares it, under the standard
/// code, so the walk passes residue 1 only where Ensembl VEP's does: a table-2
/// `ATA` first codon, Met under the table and Ile under the standard code,
/// stops it there. The residues after it follow the transcript's table like the
/// rest of the walk.
fn hgvsp_fs_peptides(
    ev: &PerlCodingEval<'_>,
    n: &mut HgvsProteinNotation,
    preds: &HgvsPredicates,
) -> Option<()> {
    let start_lost = preds.start_lost;
    let alt_cds = ev.alternate_cds()?;
    if !alt_cds
        .iter()
        .any(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'-'))
    {
        return None;
    }
    let alt_trans = perl_translate(&alt_cds, ev.table);
    let mut ref_trans = ev.pep_seq.to_vec();
    ref_trans.push(b'*');
    let standard_code_differs_at_one = || {
        alt_cds.get(..3).map(|codon| perl_translate_codon(codon, 1)) != ref_trans.first().copied()
    };
    n.start = truthy(ev.span.tl_start)?;
    if n.start > alt_trans.len() as i64 {
        n.alt_pep = Some(b"del".to_vec());
        n.kind = HgvspKind::Del;
        // One or two bases are left of the codon the walk starts at; none
        // where the change takes whole codons to the transcript's last base.
        n.cut_short = alt_cds.len() > 3 * alt_trans.len();
        return Some(());
    }
    // Ensembl VEP ends the walk on a stop both translations hold and writes
    // it unchanged (2354-2357), which is the change where no residue follows
    // that stop. For an allele whose codons, as the predicates read them,
    // begin at a stop the peptide holds and reads on past, the walk goes on
    // past such a stop to the first residue that differs.
    while n.start <= alt_trans.len() as i64 {
        let i = (n.start - 1) as usize;
        let r = ref_trans.get(i).copied();
        let a = alt_trans.get(i).copied();
        n.ref_pep = Some(r.into_iter().collect());
        n.alt_pep = Some(a.into_iter().collect());
        let read_past = preds.begins_at_a_held_stop && ev.peptide_reads_on_past(i);
        if r == Some(b'*') && a == Some(b'*') && !read_past {
            n.kind = HgvspKind::Eq;
            return Some(());
        }
        let differs = if start_lost && n.start == 1 {
            standard_code_differs_at_one()
        } else {
            r != a
        };
        if differs {
            return Some(());
        }
        n.start += 1;
    }
    start_lost.then_some(())
}

/// `_get_surrounding_peptides` (2368): `length` residues of the reference peptide
/// (plus `original_ref` when it starts with a stop) from 1-based `ref_pos`, or to
/// the end without a length; `None` when the peptide ends at or before `ref_pos`.
/// `ref_pos == 0` reads Perl's `substr(..., -1)`, the final residue.
fn hgvsp_surrounding(
    ev: &PerlCodingEval<'_>,
    ref_pos: i64,
    original_ref: &[u8],
    length: Option<usize>,
) -> Option<Vec<u8>> {
    let mut ref_trans = ev.pep_seq.to_vec();
    if original_ref.first() == Some(&b'*') {
        ref_trans.extend_from_slice(original_ref);
    }
    if ref_trans.len() as i64 <= ref_pos {
        return None;
    }
    let off = if ref_pos >= 1 {
        (ref_pos - 1) as usize
    } else if ref_pos == 0 {
        ref_trans.len() - 1
    } else {
        return None;
    };
    let stop = match length {
        Some(l) => (off + l).min(ref_trans.len()),
        None => ref_trans.len(),
    };
    Some(ref_trans[off..stop].to_vec())
}

/// `_check_for_peptide_duplication` (2442): an inserted peptide equal to the
/// residues just before it (the reference CDS translated with the transcript's
/// table, plus `preseq`) becomes a `Dup` of those residues, three-lettered here
/// and not again.
fn hgvsp_check_duplication(ev: &PerlCodingEval<'_>, n: &mut HgvsProteinNotation) {
    let alt: &[u8] = n.alt_pep.as_deref().unwrap_or_default();
    let mut upstream = perl_translate(ev.cds, ev.table);
    upstream.truncate((n.start - 1).max(0) as usize);
    upstream.extend_from_slice(&n.preseq);
    let test_new_start = n.start - alt.len() as i64 - 1;
    if test_new_start >= 0 && upstream.len() as i64 >= test_new_start + alt.len() as i64 {
        let s = test_new_start as usize;
        if &upstream[s..s + alt.len()] == alt {
            n.kind = HgvspKind::Dup;
            n.end = n.start - 1;
            n.start -= alt.len() as i64;
            n.alt_pep = Some(seq3(alt));
        }
    }
}

/// `_stop_loss_extra_AA` (2477): residues from the variant to the first stop of
/// the alternate translation under the transcript's table; counted from
/// `ref_var_pos` for a frameshift, else past the reference peptide's end. `None`
/// unless positive.
fn hgvsp_stop_loss_extra_aa(ev: &PerlCodingEval<'_>, ref_var_pos: i64, fs: bool) -> Option<i64> {
    if ref_var_pos == 0 {
        return None;
    }
    let alt_cds = ev.alternate_cds()?;
    let alt_trans = perl_translate(&alt_cds, ev.table);
    let stop_end = alt_trans.iter().position(|&b| b == b'*')? as i64 + 1;
    let extra = if fs {
        stop_end - ref_var_pos
    } else {
        stop_end - 1 - ev.pep_seq.len() as i64
    };
    (extra > 0).then_some(extra)
}

/// The count of a frameshift named at residue 1 where no start codon stands
/// at the coding start ([`start_codon_stands_at_the_coding_start`]: a
/// transcript flagged `cds_start_NF`, or one whose record shows no start
/// codon): the place of the first stop of the alternate translation under the
/// transcript's table, as a frameshift named at a later residue counts from
/// its own. The first residue is not shown to be an initiator there, so the
/// new frame is the record's own. `None` for a frameshift named at a later
/// residue, where the new frame meets no stop, and where a start codon
/// stands: the alternate translation read from residue 1 is then not shown to
/// be the variant's protein, since the bases before the coding start can form
/// the start codon again in front of a coding sequence that stands whole.
///
/// Ensembl VEP makes no count for a frameshift named at residue 1 on any
/// transcript: `_stop_loss_extra_AA` returns `undef` for the position 0 it is
/// handed (2483). [`hgvsp_stop_loss_extra_aa`] keeps that reading for every
/// count it makes, the count of an extension named at residue 1 among them.
fn hgvsp_frameshift_count_from_residue_1(
    ev: &PerlCodingEval<'_>,
    first_changed: i64,
) -> Option<i64> {
    if first_changed != 1 || ev.start_codon_stands {
        return None;
    }
    let alt_trans = perl_translate(&ev.alternate_cds()?, ev.table);
    let stop = alt_trans.iter().position(|&b| b == b'*')?;
    Some(stop as i64 + 1)
}

/// `_get_del_peptides` (2544), Ensembl VEP's path for a deletion whose reference peptide
/// window is empty: both peptides from `translation_start` to the end (the alternate
/// side, translated under the transcript's table, cut at its first stop), clipped,
/// three-lettered.
fn hgvsp_del_peptides(ev: &PerlCodingEval<'_>, n: &mut HgvsProteinNotation) -> Option<()> {
    let alt_cds = ev.alternate_cds()?;
    let tl_start = truthy(ev.span.tl_start)?;
    let start0 = ((tl_start - 1).max(0)) as usize;
    let alt_trans = perl_translate(&alt_cds, ev.table);
    let alt_tail = alt_trans.get(start0..).unwrap_or(b"");
    let alt: Vec<u8> = alt_tail
        .split(|&b| b == b'*')
        .next()
        .unwrap_or(b"")
        .to_vec();
    let ref_tail = ev.pep_seq.get(start0..).unwrap_or(b"").to_vec();
    n.alt_pep = Some(alt);
    n.ref_pep = Some(ref_tail);
    n.start = tl_start;
    hgvsp_clip_alleles(n);
    n.alt_pep = Some(seq3(n.alt_pep.as_deref().unwrap_or(b"")));
    n.ref_pep = Some(seq3(n.ref_pep.as_deref().unwrap_or(b"")));
    Some(())
}

/// `_check_peptides_post_var` (2573) plus `_shift_3prime` (2595): rotates an
/// inserted or deleted peptide along the residues after `end` while its first
/// residue matches, moving `start` and `end` with it.
fn hgvsp_post_var_shift(ev: &PerlCodingEval<'_>, n: &mut HgvsProteinNotation) {
    let Some(post_seq) = hgvsp_surrounding(ev, n.end + 1, &n.original_ref, None) else {
        return;
    };
    let seq_to_check = match n.kind {
        HgvspKind::Ins => n.alt_pep.get_or_insert_default(),
        HgvspKind::Del => n.ref_pep.get_or_insert_default(),
        _ => return,
    };
    let deleted_length = seq_to_check.len() as i64;
    let mut i = 0i64;
    while i <= post_seq.len() as i64 - deleted_length {
        let next_del = seq_to_check.first().copied();
        let next_post = post_seq.get(i as usize).copied();
        if next_del.is_some() && next_del == next_post {
            n.start += 1;
            n.end += 1;
            seq_to_check.rotate_left(1);
        } else {
            break;
        }
        i += 1;
    }
}

/// What a change does to the protein at the stop codon, read from the reference
/// protein and the alternate protein alone, as the HGVS nomenclature describes a
/// protein change: by comparing the two sequences, whatever the change is on the
/// DNA level.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StopCodonChange {
    /// The alternate protein is the reference protein, its stop at the same
    /// residue.
    Unchanged,
    /// The whole reference protein is kept and residues follow it: `first`
    /// stands in the stop's place and `added` counts the residues up to the new
    /// stop, the one in the stop's place included, `None` where the transcript
    /// holds no new stop. `stop_inserted` is whether the new stop codon lies
    /// wholly inside the inserted bases.
    Extended {
        first: u8,
        added: Option<usize>,
        stop_inserted: bool,
    },
}

/// The forms in which the description built from the codon window names the
/// stop of the reference protein at its own residue, inserts a stop that
/// stands at that residue, or reads an insertion after the last residue and
/// describes nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WrittenAtStop {
    /// A deletion-insertion of the stop alone.
    Replaced,
    /// An insertion between the last residue and the stop.
    InsertedBefore,
    /// An insertion between two neighbouring residues, described in the
    /// coding sequence's frame, that names one residue or more and then a
    /// stop as inserted. The first stop so named stands at the stop's
    /// residue: the number of the residue before the insertion, the
    /// residues inserted before that stop, and one.
    InsertedWithItsStop,
    /// An insertion between two bases of the last residue's codon, of any
    /// length, that keeps that residue: the codon window reads an insertion
    /// after the last residue, and [`hgvsp_peptides`] names no residue after
    /// it and gives no description.
    InsertedInTheLastCodon,
    /// The stop written as unchanged; `stop_lost` is the allele's own
    /// `stop_lost` predicate, the verdict the row's terms carry.
    Unchanged { stop_lost: bool },
    /// The stop written as deleted, with the count of an extension after it
    /// and no residue in the stop's place.
    DeletedWithExtension,
}

impl PerlCodingEval<'_> {
    /// The half-open range of the alternate coding sequence that the inserted
    /// bases fill, less the bases the two alleles share at either end.
    fn inserted_cds_bases(&self, cds_start: i64, cds_end: i64) -> std::ops::Range<usize> {
        let up_end = (cds_start - 1) as usize;
        let deleted = &self.cds[up_end..cds_end as usize];
        let inserted = self.alt_allele_on_transcript();
        let front = deleted
            .iter()
            .zip(&inserted)
            .take_while(|(d, i)| d == i)
            .count();
        let back = deleted[front..]
            .iter()
            .rev()
            .zip(inserted[front..].iter().rev())
            .take_while(|(d, i)| d == i)
            .count();
        up_end + front..up_end + inserted.len() - back
    }

    /// Whether the allele lies inside one exon: every base it replaces, or for
    /// an insertion both bases it lies between.
    fn inside_one_exon(&self) -> bool {
        let (first, last) = (
            self.vf_start.min(self.vf_end),
            self.vf_start.max(self.vf_end),
        );
        let mapper = self
            .transcript
            .vefc
            .as_ref()
            .and_then(|v| v.mapper.as_ref());
        mapper.is_some_and(|mapper| {
            mapper
                .exon_coord_mapper
                .pairs
                .iter()
                .any(|exon| exon.to_start as i64 <= first && last <= exon.to_end as i64)
        })
    }

    /// The allele's span in the coding sequence, where the alternate protein
    /// can be held against the reference protein: the transcript's peptide is
    /// the translation of a coding sequence that one stop codon closes, and the
    /// allele lies inside one exon and inside the coding sequence (at an exon's
    /// edge the spliced sequence is not known). An allele that holds unknown
    /// bases is compared only where its known bases keep the stop codon
    /// before the first of them
    /// ([`Self::keeps_a_stop_codon_before_unknown_bases`]): the alternate
    /// protein then ends before any unknown base. `None` anywhere else.
    fn span_where_the_proteins_compare(&self) -> Option<(i64, i64)> {
        let residues = self.pep_seq.len();
        let reference = perl_translate(self.cds, self.table);
        if self.cds.len() != 3 * (residues + 1)
            || reference[..residues] != self.pep_seq[..]
            || reference[residues] != b'*'
            || self.pep_seq.contains(&b'*')
            || !(is_perl_unambiguous_dna(self.alt_allele) || self.reads_a_kept_stop_codon)
            || !self.inside_one_exon()
        {
            return None;
        }
        let (cds_start, cds_end) = (self.span.cds_start?, self.span.cds_end?);
        if cds_start < 1 || cds_end > self.cds.len() as i64 || cds_start > cds_end + 1 {
            return None;
        }
        Some((cds_start, cds_end))
    }

    /// The frameshift that a change reaching from before the stop codon into it
    /// is, where the HGVS nomenclature decides it: the first changed residue
    /// (1-based), the residues the alternate protein holds from it up to its
    /// new stop, and whether the transcript holds that stop. The nomenclature
    /// describes deletions "starting 5' of and including the translation
    /// termination (stop) codon" as a frameshift, says the same of
    /// deletion-insertions, and defines a frameshift by translation that
    /// "shifts to another reading frame", so the change, less the bases its two
    /// alleles share at either end, either runs to the stop codon's last base
    /// or is not a whole number of codons. `None` where the proteins are not
    /// compared ([`Self::span_where_the_proteins_compare`]), where no residue
    /// before the stop changes, where the alternate protein keeps the
    /// reference's length and its stop's own residue, and where the first
    /// changed residue is itself a stop.
    fn frameshift_through_the_stop(&self) -> Option<(usize, Vec<u8>, bool)> {
        let residues = self.pep_seq.len();
        let (cds_start, cds_end) = self.span_where_the_proteins_compare()?;
        let stop_codon = self.cds.len() as i64 - 3;
        let inserted = self.alt_allele_on_transcript().len() as i64;
        // The change itself, however the record writes it: the span less the
        // bases the two alleles share at either end.
        let kept = self.inserted_cds_bases(cds_start, cds_end);
        let shared_before = kept.start as i64 - (cds_start - 1);
        let shared_after = cds_start - 1 + inserted - kept.end as i64;
        let (first, last) = (cds_start + shared_before, cds_end - shared_after);
        let whole_stop_codon = last == self.cds.len() as i64;
        let out_of_frame = (cds_end - cds_start + 1 - inserted) % 3 != 0;
        if first > stop_codon || last <= stop_codon || !(whole_stop_codon || out_of_frame) {
            return None;
        }
        let alternate = perl_translate(&self.alternate_cds()?, self.table);
        let same = alternate
            .iter()
            .zip(self.pep_seq.iter())
            .take_while(|(alt, reference)| alt == reference)
            .count();
        let new_stop = alternate.iter().position(|&aa| aa == b'*');
        let end = new_stop.unwrap_or(alternate.len());
        if same >= residues || new_stop == Some(residues) || end <= same {
            return None;
        }
        Some((same + 1, alternate[same..end].to_vec(), new_stop.is_some()))
    }

    /// Whether the allele, less the bases its two alleles share at either
    /// end, is an insertion between two bases of the last residue's codon, of
    /// any length.
    fn inserts_in_the_last_residues_codon(&self) -> bool {
        let (Some(cds_start), Some(cds_end)) = (self.span.cds_start, self.span.cds_end) else {
            return false;
        };
        if cds_start < 1 || cds_end > self.cds.len() as i64 || cds_start > cds_end + 1 {
            return false;
        }
        let inserted = self.alt_allele_on_transcript().len() as i64;
        let kept = self.inserted_cds_bases(cds_start, cds_end);
        let shared_before = kept.start as i64 - (cds_start - 1);
        let shared_after = cds_start - 1 + inserted - kept.end as i64;
        let (first, last) = (cds_start + shared_before, cds_end - shared_after);
        let last_codon_end = 3 * self.pep_seq.len() as i64;
        first == last + 1 && last_codon_end - 2 <= last && first <= last_codon_end
    }

    /// Whether a stop codon closes the coding sequence: it is a whole number
    /// of codons, the last of them is a stop codon, and either the transcript's
    /// peptide holds one residue for every codon before it or no other codon
    /// is a stop codon. The codons are the triplets from the first base, which
    /// the start exon's phase pads into frame, and the translation leaves out
    /// the one or two bases after the last whole codon (`Transcript.pm`
    /// `translateable_seq` and `translate`). A coding sequence that ends inside
    /// a codon is so closed by none, whatever its last three bases spell.
    fn stop_codon_closes_the_coding_sequence(&self) -> bool {
        if !self.cds.len().is_multiple_of(3) {
            return false;
        }
        let Some(last_codon) = self.cds.len().checked_sub(3) else {
            return false;
        };
        if perl_translate_codon(&self.cds[last_codon..], self.table) != b'*' {
            return false;
        }
        self.cds.len() == 3 * (self.pep_seq.len() + 1)
            || !self.cds[..last_codon]
                .chunks_exact(3)
                .any(|codon| perl_translate_codon(codon, self.table) == b'*')
    }

    /// The count of residues that follow the whole reference protein up to a
    /// new stop, where the coding sequence holds stop codons that the
    /// transcript's peptide reads as residues and the allele replaces bases of
    /// the closing stop codon by as many. The peptide is the reference protein,
    /// and a stop codon inside the coding sequence reads as the peptide reads
    /// it, where a translation of the bases ends at the first of them. `None`
    /// where the coding sequence is not of that kind, where the transcript holds
    /// no new stop, and where the new stop is a codon that the peptide reads as
    /// a residue inside the coding sequence: whether such a codon ends the
    /// protein in the 3' UTR is not known from the model.
    fn residues_added_past_stops_read_as_residues(&self) -> Option<i64> {
        let residues = self.pep_seq.len();
        let reference = perl_translate(self.cds, self.table);
        if self.cds.len() != 3 * (residues + 1)
            || reference[residues] != b'*'
            || self.pep_seq.contains(&b'*')
            || !is_perl_unambiguous_dna(self.alt_allele)
            || !self.inside_one_exon()
        {
            return None;
        }
        let read_as_residues: Vec<usize> =
            (0..residues).filter(|&at| reference[at] == b'*').collect();
        let fits =
            (0..residues).all(|at| reference[at] == b'*' || reference[at] == self.pep_seq[at]);
        if read_as_residues.is_empty() || !fits {
            return None;
        }
        let (cds_start, cds_end) = (self.span.cds_start?, self.span.cds_end?);
        let stop_codon = self.cds.len() as i64 - 3;
        let replaced = cds_end - cds_start + 1;
        if cds_start <= stop_codon
            || cds_end > self.cds.len() as i64
            || replaced != self.alt_allele_on_transcript().len() as i64
        {
            return None;
        }
        let alternate_cds = self.alternate_cds()?;
        let alternate = perl_translate(&alternate_cds, self.table);
        let added = alternate
            .get(residues..)?
            .iter()
            .position(|&aa| aa == b'*')?;
        let new_stop = alternate_cds.get(3 * (residues + added)..3 * (residues + added) + 3)?;
        let read_as_a_residue_inside = read_as_residues
            .iter()
            .any(|&at| self.cds[3 * at..3 * at + 3].eq_ignore_ascii_case(new_stop));
        (added > 0 && !read_as_a_residue_inside).then_some(added as i64)
    }

    /// The count of residues from a stop that a substitution changes to the
    /// next stop of the frame, on a transcript whose peptide holds a stop and
    /// names the residues after it (a polymorphic pseudogene). `lost` is the
    /// changed residue's number: a stop the peptide holds, or the stop codon
    /// that closes the coding sequence. `_stop_loss_extra_AA` subtracts the
    /// peptide's length from the place of the first stop of the whole
    /// alternate translation, and returns no count where the changed residue
    /// is the first (TranscriptVariationAllele.pm:2503, 2513 and 2483). On
    /// such a peptide that stop is not the new one where a held stop stands
    /// before the change, and the peptide's last residue is not the one
    /// before a held stop. The next stop is sought from the changed residue
    /// through the alternate coding sequence and the transcript's 3' UTR, and
    /// no farther. The outer `None` where the count is not read here: the
    /// peptide holds no stop or is not the translation of the coding sequence
    /// at each of its codons (the first residue a methionine whatever its
    /// codon), the allele holds an ambiguous base, the residue `lost` is
    /// neither kind of stop, a residue before it changes, or the alternate
    /// bases hold a stop at it. The inner `None` where no stop follows.
    fn residues_from_a_lost_stop_past_stops_the_peptide_holds(
        &self,
        lost: i64,
    ) -> Option<Option<i64>> {
        if !self.pep_seq.contains(&b'*') || !is_perl_unambiguous_dna(self.alt_allele) {
            return None;
        }
        let at = usize::try_from(lost.checked_sub(1)?).ok()?;
        let reference = perl_translate(self.cds, self.table);
        let residues = self.pep_seq.len().min(reference.len());
        let fits = (0..residues).all(|residue| {
            reference[residue] == self.pep_seq[residue] || (residue == 0 && self.pep_seq[0] == b'M')
        });
        let held = at < residues && self.pep_seq[at] == b'*';
        let closing = at == self.pep_seq.len()
            && self.cds.len() == 3 * (at + 1)
            && reference.get(at) == Some(&b'*');
        if !fits || !(held || closing) {
            return None;
        }
        let alternate = perl_translate(&self.alternate_cds()?, self.table);
        if alternate.get(..at)? != &reference[..at] || alternate.get(at)? == &b'*' {
            return None;
        }
        Some(
            alternate[at..]
                .iter()
                .position(|&aa| aa == b'*')
                .map(|residues_on| residues_on as i64),
        )
    }

    /// The count of residues by which a deletion of whole residues that takes
    /// a stop moves the end of the protein, on a transcript whose peptide
    /// holds a stop and names the residues after it (a polymorphic
    /// pseudogene): the place of the first stop the alternate bases meet at
    /// or after the first changed residue, less the place of the stop that
    /// is lost, the first stop of the reference at or after that residue.
    /// `first_changed` is the number of the first residue the description
    /// names. `_stop_loss_extra_AA` subtracts one and the peptide's length
    /// from the place of the first stop of the whole alternate translation
    /// (TranscriptVariationAllele.pm:2503 and 2513), which is this figure
    /// where the lost stop is the codon after the peptide's last residue and
    /// no stop stands before it. On such a peptide a held stop stands before
    /// the change, or the lost stop is not that codon. The next stop is
    /// sought through the alternate coding sequence and the transcript's 3'
    /// UTR, and no farther. `None` where the count is not read here: the
    /// peptide holds no stop or is not the translation of the coding
    /// sequence at each of its codons (the first residue a methionine
    /// whatever its codon), the allele holds an ambiguous base, the first
    /// residue that differs is not `first_changed`, the lost stop is neither
    /// a stop the peptide holds nor the stop codon that closes the coding
    /// sequence, no stop follows, or the figure is not above nought, which
    /// `_stop_loss_extra_AA` leaves without a count on any transcript (2519).
    fn residues_a_deletion_adds_past_a_lost_stop_the_peptide_holds(
        &self,
        first_changed: i64,
    ) -> Option<i64> {
        if !self.pep_seq.contains(&b'*') || !is_perl_unambiguous_dna(self.alt_allele) {
            return None;
        }
        let reference = perl_translate(self.cds, self.table);
        let residues = self.pep_seq.len().min(reference.len());
        let fits = (0..residues).all(|residue| {
            reference[residue] == self.pep_seq[residue] || (residue == 0 && self.pep_seq[0] == b'M')
        });
        let alternate = perl_translate(&self.alternate_cds()?, self.table);
        let first = alternate
            .iter()
            .zip(&reference)
            .position(|(alternate, reference)| alternate != reference)?;
        if !fits || i64::try_from(first).ok()? + 1 != first_changed {
            return None;
        }
        let lost = first + reference[first..].iter().position(|&aa| aa == b'*')?;
        let held = lost < residues && self.pep_seq[lost] == b'*';
        let closing = lost == self.pep_seq.len() && self.cds.len() == 3 * (lost + 1);
        let new = first + alternate[first..].iter().position(|&aa| aa == b'*')?;
        ((held || closing) && new > lost).then(|| (new - lost) as i64)
    }

    /// The count from the first changed residue of a frameshift to the new
    /// stop, where a stop codon inside the coding sequence stands before
    /// that residue and the transcript's peptide goes on past it: the
    /// peptide reads the codon as a residue (a selenocysteine codon) or
    /// holds it as a stop and names the residues after it (a polymorphic
    /// pseudogene). A translation of the alternate bases ends at the first
    /// such codon, before the change (`_stop_loss_extra_AA`,
    /// TranscriptVariationAllele.pm:2503), and gives no count; the peptide
    /// is the reference protein, and the allele changes no residue before
    /// `first_changed`, the first changed residue's number. The new stop is
    /// sought from that residue through the alternate coding sequence and
    /// the transcript's 3' UTR, and no farther. `None` where no such codon
    /// stands before that residue, where the peptide is not the translation
    /// of the coding sequence at each of its other codons (the first residue
    /// a methionine whatever its codon: `Transcript::translate` writes one
    /// over a start codon of the table and an `initial_met` edit over any
    /// other), where the allele holds an ambiguous base, where the first
    /// changed residue is not one of the peptide's or is a stop it holds,
    /// where the new frame meets no stop, and where the new stop is a codon
    /// that the peptide reads as a residue inside the coding sequence:
    /// whether such a codon ends the protein there is not known from the
    /// model.
    fn residues_to_a_new_stop_past_stops_inside_the_coding_sequence(
        &self,
        first_changed: i64,
    ) -> Option<i64> {
        let residues = self.pep_seq.len();
        let first = usize::try_from(first_changed.checked_sub(1)?).ok()?;
        let reference = perl_translate(self.cds, self.table);
        let residues = residues.min(reference.len());
        if first >= residues
            || self.pep_seq[first] == b'*'
            || !is_perl_unambiguous_dna(self.alt_allele)
        {
            return None;
        }
        let fits = (0..residues).all(|at| {
            reference[at] == self.pep_seq[at]
                || (reference[at] == b'*' && self.pep_seq[at] != b'*')
                || (at == 0 && self.pep_seq[0] == b'M')
        });
        let stands_before = reference[..first].contains(&b'*');
        if !fits || !stands_before {
            return None;
        }
        let read_as_residues: Vec<usize> = (0..residues)
            .filter(|&at| reference[at] == b'*' && self.pep_seq[at] != b'*')
            .collect();
        let alternate_cds = self.alternate_cds()?;
        let alternate = perl_translate(&alternate_cds, self.table);
        if alternate.get(..first)? != &reference[..first] {
            return None;
        }
        let to_the_stop = alternate[first..].iter().position(|&aa| aa == b'*')?;
        let at = first + to_the_stop;
        let new_stop = alternate_cds.get(3 * at..3 * at + 3)?;
        let read_as_a_residue_inside = read_as_residues
            .iter()
            .any(|&inside| self.cds[3 * inside..3 * inside + 3].eq_ignore_ascii_case(new_stop));
        (to_the_stop > 0 && !read_as_a_residue_inside).then_some(to_the_stop as i64 + 1)
    }

    /// The count of residues from a stop that is the first changed residue of
    /// a frameshift to the new frame's first stop, on a transcript whose
    /// peptide holds a stop and names the residues after it (a polymorphic
    /// pseudogene). `lost` is that residue's number: a stop the peptide
    /// holds, or the stop codon that closes the coding sequence. Both
    /// engines write such a frameshift as an extension named at the stop
    /// (TranscriptVariationAllele.pm:2007 to 2010). `_stop_loss_extra_AA`
    /// counts it from the first stop of the whole alternate translation
    /// (2503 and 2508), which is a stop the peptide holds before the change
    /// wherever one stands there, and then gives no count. The new stop is
    /// sought from the changed residue through the alternate coding sequence
    /// and the transcript's 3' UTR, and no farther. `None` where the peptide
    /// holds no stop or is not the translation of the coding sequence at
    /// each of its codons (the first residue a methionine whatever its
    /// codon), where the allele holds an ambiguous base, where the residue
    /// `lost` is neither kind of stop, where a residue before it changes or
    /// the alternate bases hold a stop at it, and where the new frame meets
    /// no stop.
    fn residues_from_a_stop_a_frameshift_changes_first(&self, lost: i64) -> Option<i64> {
        if !self.pep_seq.contains(&b'*') || !is_perl_unambiguous_dna(self.alt_allele) {
            return None;
        }
        let at = usize::try_from(lost.checked_sub(1)?).ok()?;
        let reference = perl_translate(self.cds, self.table);
        let residues = self.pep_seq.len().min(reference.len());
        let fits = (0..residues).all(|residue| {
            reference[residue] == self.pep_seq[residue] || (residue == 0 && self.pep_seq[0] == b'M')
        });
        let held = at < residues && self.pep_seq[at] == b'*';
        let closing = at == self.pep_seq.len()
            && self.cds.len() == 3 * (at + 1)
            && reference.get(at) == Some(&b'*');
        if !fits || !(held || closing) {
            return None;
        }
        let alternate = perl_translate(&self.alternate_cds()?, self.table);
        if alternate.get(..at)? != &reference[..at] || alternate.get(at)? == &b'*' {
            return None;
        }
        alternate[at..]
            .iter()
            .position(|&aa| aa == b'*')
            .map(|residues_on| residues_on as i64)
    }

    /// The change at the stop codon that comparing the alternate protein with
    /// the reference protein shows. The reference protein is the transcript's
    /// peptide where it is the translation of a coding sequence that one stop
    /// codon closes; the alternate protein is the translation of the alternate
    /// coding sequence, read on through the 3' UTR to its first stop. `None`
    /// where the proteins are not compared
    /// ([`Self::span_where_the_proteins_compare`]), where a residue of the
    /// reference protein changes, and where no whole codon follows the
    /// reference protein.
    fn stop_codon_change(&self) -> Option<StopCodonChange> {
        let residues = self.pep_seq.len();
        let (cds_start, cds_end) = self.span_where_the_proteins_compare()?;
        let alternate = perl_translate(&self.alternate_cds()?, self.table);
        let tail = alternate.strip_prefix(&self.pep_seq[..])?;
        let first = *tail.first()?;
        let added = tail.iter().position(|&aa| aa == b'*');
        if added == Some(0) {
            return Some(StopCodonChange::Unchanged);
        }
        let stop_inserted = added.is_some_and(|added| {
            let new_stop = 3 * (residues + added);
            let inserted = self.inserted_cds_bases(cds_start, cds_end);
            inserted.start <= new_stop && new_stop + 3 <= inserted.end
        });
        Some(StopCodonChange::Extended {
            first,
            added,
            stop_inserted,
        })
    }
}

/// Whether a form written at the stop codon gives way to the description the
/// two proteins bear out. Each arm is one form: the description built from the
/// codon window, or its lack of one, on the left, what the comparison of the
/// proteins shows on the right.
fn written_at_stop_gives_way(written: WrittenAtStop, change: StopCodonChange) -> bool {
    match (written, change) {
        // A deletion-insertion that replaces the stop by a stop, where the
        // protein does not change.
        (WrittenAtStop::Replaced, StopCodonChange::Unchanged) => true,

        // A deletion-insertion of the stop alone, where the whole protein is
        // kept and residues follow it up to a stop of the 3' UTR or to the
        // transcript's end.
        (
            WrittenAtStop::Replaced,
            StopCodonChange::Extended {
                stop_inserted: false,
                ..
            },
        ) => true,

        // The same where the new stop lies wholly inside the inserted bases.
        (
            WrittenAtStop::Replaced,
            StopCodonChange::Extended {
                stop_inserted: true,
                ..
            },
        ) => true,

        // An insertion between the last residue and the stop, where the whole
        // protein is kept and residues follow it.
        (WrittenAtStop::InsertedBefore, StopCodonChange::Extended { .. }) => true,

        // An insertion whose description names a stop, after one inserted
        // residue or more, at the stop's residue, where the protein does
        // not change: the nomenclature describes an insertion only for a
        // protein that gains residues.
        (WrittenAtStop::InsertedWithItsStop, StopCodonChange::Unchanged) => true,

        // An insertion inside the last residue's codon that the codon window
        // leaves without a description, where the protein does not change:
        // the nomenclature describes the two proteins, and they are equal.
        (WrittenAtStop::InsertedInTheLastCodon, StopCodonChange::Unchanged) => true,

        // The stop written as unchanged on an allele whose own verdict is
        // `stop_lost`, where residues follow the whole protein.
        (WrittenAtStop::Unchanged { stop_lost: true }, StopCodonChange::Extended { .. }) => true,

        // The stop written as deleted with an extension count, where the whole
        // protein is kept and residues follow it: the extension names the
        // residue that stands in the stop's place.
        (WrittenAtStop::DeletedWithExtension, StopCodonChange::Extended { .. }) => true,

        _ => false,
    }
}

/// The residue number at which the stop of an insertion's own residues
/// stands: the number of the residue before the insertion, the residues
/// inserted before their first stop, and one. `None` for a description that
/// is no insertion between two neighbouring residues and for inserted
/// residues that hold no stop. The inserted residues are read in their
/// three-letter names, as [`hgvsp_peptides`] leaves them.
fn residue_of_the_inserted_stop(n: &HgvsProteinNotation) -> Option<i64> {
    let (first, last) = (n.start.min(n.end), n.start.max(n.end));
    if n.kind != HgvspKind::Ins || last != first + 1 {
        return None;
    }
    let inserted = n.alt_pep.as_deref()?;
    let before_the_stop = inserted.chunks(3).position(|name| name == b"Ter")?;
    Some(first + before_the_stop as i64 + 1)
}

/// The notation of a change at the stop codon as the HGVS nomenclature gives it
/// for the two proteins, in place of a form of [`WrittenAtStop`] that the
/// codon window's description takes and
/// [`written_at_stop_gives_way`] names: the stop at its residue number as
/// unchanged where the protein does not change, and an extension named at the
/// stop (the stop's residue number, the first new residue, the count of
/// residues added or `?` where the transcript holds no new stop) where the
/// whole reference protein is kept and residues are added after it.
/// `unchanged` is whether the codon window's description reads the two
/// peptides as equal. `None` leaves that description as it is.
fn hgvsp_at_stop_codon(
    ev: &PerlCodingEval<'_>,
    preds: &HgvsPredicates,
    n: &HgvsProteinNotation,
    unchanged: bool,
) -> Option<String> {
    let stop = ev.pep_seq.len() as i64 + 1;
    let (first, last) = (n.start.min(n.end), n.start.max(n.end));
    let written = if unchanged && n.start == stop {
        WrittenAtStop::Unchanged {
            stop_lost: preds.stop_lost,
        }
    } else if unchanged {
        return None;
    } else if n.kind == HgvspKind::Delins && first == stop && last == stop {
        WrittenAtStop::Replaced
    } else if n.kind == HgvspKind::Ins && first + 1 == stop && last == stop {
        WrittenAtStop::InsertedBefore
    } else if residue_of_the_inserted_stop(n) == Some(stop) {
        WrittenAtStop::InsertedWithItsStop
    } else if n.kind == HgvspKind::Del && first == stop && last == stop && preds.stop_lost {
        WrittenAtStop::DeletedWithExtension
    } else {
        return None;
    };
    hgvsp_the_proteins_bear_out(ev, written)
}

/// The notation of an insertion that [`hgvsp_peptides`] leaves without one,
/// naming no residue after the last, where it is
/// [`WrittenAtStop::InsertedInTheLastCodon`] and the protein does not change:
/// the stop at its residue number as unchanged. `None` for every other
/// description [`hgvsp_peptides`] gives up, and for that form where the
/// proteins are not compared or differ.
fn hgvsp_of_an_insertion_no_residue_follows(
    ev: &PerlCodingEval<'_>,
    n: &HgvsProteinNotation,
) -> Option<String> {
    if n.kind != HgvspKind::Ins || !ev.inserts_in_the_last_residues_codon() {
        return None;
    }
    hgvsp_the_proteins_bear_out(ev, WrittenAtStop::InsertedInTheLastCodon)
}

/// The notation that the two proteins bear out in place of the form `written`,
/// where [`written_at_stop_gives_way`] names the pair. `None` where the
/// proteins are not compared and where the form stands.
fn hgvsp_the_proteins_bear_out(ev: &PerlCodingEval<'_>, written: WrittenAtStop) -> Option<String> {
    let stop = ev.pep_seq.len() as i64 + 1;
    let change = ev.stop_codon_change()?;
    if !written_at_stop_gives_way(written, change) {
        return None;
    }
    Some(match change {
        StopCodonChange::Unchanged => format!("Ter{stop}="),
        StopCodonChange::Extended { first, added, .. } => {
            let residue = crate::hgvs::amino_acid_three_letter(first.to_ascii_uppercase());
            let count = added.map_or("?".to_string(), |added| added.to_string());
            format!("Ter{stop}{residue}extTer{count}")
        }
    })
}

/// The frameshift notation of a change that the codon window writes as a
/// deletion-insertion of a range ending at the stop, where
/// [`PerlCodingEval::frameshift_through_the_stop`] names one: the first changed
/// residue with its number, the new residue, and the new stop's position
/// counted from the first changed residue as 1, `?` where the transcript holds
/// no new stop. `None` leaves the codon window's description as it is.
fn hgvsp_frameshift_through_the_stop(
    ev: &PerlCodingEval<'_>,
    n: &HgvsProteinNotation,
) -> Option<String> {
    let stop = ev.pep_seq.len() as i64 + 1;
    if n.kind != HgvspKind::Delins || n.start >= n.end || n.end != stop {
        return None;
    }
    let (first, new_residues, stop_met) = ev.frameshift_through_the_stop()?;
    let three = |aa: u8| crate::hgvs::amino_acid_three_letter(aa.to_ascii_uppercase());
    let count = if stop_met {
        (new_residues.len() + 1).to_string()
    } else {
        "?".to_string()
    };
    Some(format!(
        "{}{first}{}fsTer{count}",
        three(ev.pep_seq[first - 1]),
        three(new_residues[0])
    ))
}

/// The residue number of a description that names the stop alone, and whether
/// it writes the stop as unchanged: `Ter<N>=`, the stop deleted (`Ter<N>del`),
/// replaced by one residue (`Ter<N><Xaa>`), extended
/// (`Ter<N><Xaa>extTer<count>`) or named as the first residue of a frameshift
/// (`Ter<N>fs`). `None` for every other description.
fn stop_named_alone(description: &str) -> Option<(i64, bool)> {
    let rest = description.strip_prefix("Ter")?;
    let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
    let residue = rest[..digits].parse::<i64>().ok()?;
    let change = &rest[digits..];
    let residue_name = |text: &str| {
        let name = text.as_bytes();
        name.len() == 3
            && name[0].is_ascii_uppercase()
            && name[1..].iter().all(u8::is_ascii_lowercase)
    };
    let extension = change.split_once("extTer").is_some_and(|(first, count)| {
        residue_name(first)
            && (count == "?" || (!count.is_empty() && count.bytes().all(|b| b.is_ascii_digit())))
    });
    (change == "=" || change == "del" || change == "fs" || residue_name(change) || extension)
        .then_some((residue, change == "="))
}

/// Whether a description names, alone, a stop that the reference protein does
/// not hold. Two paths write such a stop. The codon window writes a partial
/// codon and an unknown residue as the stop, so it names a stop on a coding
/// sequence that no stop codon closes and at residue 1 of one that opens
/// inside its first codon. The frameshift walk, [`hgvsp_fs_peptides`], closes
/// the reference translation with a stop whatever the coding sequence's end
/// (`_get_fs_peptides`, 2338), so a frameshift that reaches past the
/// reference protein's last residue names a stop there on a coding sequence
/// that no stop codon closes. The HGVS nomenclature takes as a reference only
/// a transcript with "a complete CDS", whose "translation termination codon"
/// and "first three nucleotides of the CDS" are annotated, so it gives no
/// description of such a stop. That is: one of the five descriptions of
/// [`stop_named_alone`] on a coding sequence that no stop codon closes, at a
/// residue where the transcript's peptide carries no stop; and the stop
/// written as unchanged at residue 1 of a protein that holds a residue there.
fn names_a_stop_the_reference_lacks(ev: &PerlCodingEval<'_>, description: &str) -> bool {
    let Some((residue, unchanged)) = stop_named_alone(description) else {
        return false;
    };
    if ev.stop_codon_closes_the_coding_sequence() {
        return unchanged && residue == 1 && !ev.pep_seq.is_empty();
    }
    !usize::try_from(residue - 1)
        .ok()
        .and_then(|at| ev.pep_seq.get(at))
        .is_some_and(|&aa| aa == b'*')
}

/// `_get_hgvs_protein_format` (1898) with three-letter conversion on and no
/// prediction parentheses, behind [`hgvsp_at_stop_codon`] and
/// [`hgvsp_frameshift_through_the_stop`].
fn hgvsp_format(
    ev: &PerlCodingEval<'_>,
    preds: &HgvsPredicates,
    n: &HgvsProteinNotation,
) -> String {
    let ref_pep: &[u8] = n.ref_pep.as_deref().unwrap_or_default();
    let alt: &[u8] = n.alt_pep.as_deref().unwrap_or_default();
    let (start, end, kind) = (n.start, n.end, n.kind);
    let text = |s: &[u8]| String::from_utf8_lossy(s).into_owned();
    let first3 = |s: &[u8]| text(&s[..3.min(s.len())]);
    let last3 = |s: &[u8]| text(&s[s.len().saturating_sub(3)..]);
    let unchanged = ref_pep == alt && kind != HgvspKind::Fs && kind != HgvspKind::Ins;
    if let Some(at_stop) = hgvsp_at_stop_codon(ev, preds, n, unchanged) {
        return at_stop;
    }
    if let Some(frameshift) = hgvsp_frameshift_through_the_stop(ev, n) {
        return frameshift;
    }
    if unchanged {
        format!("{}{start}=", text(ref_pep))
    } else if preds.stop_lost && (kind == HgvspKind::Del || kind == HgvspKind::Sub) {
        // `_stop_loss_extra_AA` ends at the first stop codon of the bases and
        // counts from the peptide's length, so a stop codon inside the coding
        // sequence that the peptide goes on past, reading it as a residue or
        // holding it as a stop, leaves the count unknown or short.
        let from_the_lost_stop = if kind == HgvspKind::Sub {
            ev.residues_from_a_lost_stop_past_stops_the_peptide_holds(start)
        } else {
            ev.residues_a_deletion_adds_past_a_lost_stop_the_peptide_holds(start)
                .map(Some)
        };
        let extra = from_the_lost_stop.unwrap_or_else(|| {
            hgvsp_stop_loss_extra_aa(ev, start - 1, false).or_else(|| {
                let at_the_closing_stop =
                    kind == HgvspKind::Sub && start == ev.pep_seq.len() as i64 + 1;
                at_the_closing_stop
                    .then(|| ev.residues_added_past_stops_read_as_residues())
                    .flatten()
            })
        });
        let aa_til_stop = extra.map_or("?".to_string(), |extra| extra.to_string());
        let alt = format!("{}extTer{aa_til_stop}", first3(alt));
        if ref_pep.len() > 3 && kind == HgvspKind::Del {
            format!("{}{start}_{}{end}{alt}", first3(ref_pep), last3(ref_pep))
        } else {
            format!("{}{start}{alt}", text(ref_pep))
        }
    } else if kind == HgvspKind::Dup {
        if start < end {
            format!("{}{start}_{}{end}dup", first3(alt), last3(alt))
        } else {
            format!("{}{start}dup", text(alt))
        }
    } else if kind == HgvspKind::Sub {
        format!("{}{start}{}", text(ref_pep), text(alt))
    } else if kind == HgvspKind::Delins || kind == HgvspKind::Ins {
        // `s/Ter\w+/Ter/`: nothing after the first stop is reported.
        let mut alt = alt.to_vec();
        if let Some(p) = alt.windows(3).position(|w| w == b"Ter") {
            let word_len = alt[p + 3..]
                .iter()
                .take_while(|b| b.is_ascii_alphanumeric() || **b == b'_')
                .count();
            alt.drain(p + 3..p + 3 + word_len);
        }
        let mut alt = text(&alt);
        let ref_first = first3(ref_pep);
        let ref_ends_in_x = ref_pep.last() == Some(&b'X');
        let ref_last = if ref_ends_in_x {
            "Ter".to_string()
        } else {
            last3(ref_pep)
        };
        if ref_ends_in_x {
            if let Some(extra) = hgvsp_stop_loss_extra_aa(ev, start - 1, false) {
                alt.push_str(&format!("extTer{extra}"));
            }
        }
        if start == end && kind == HgvspKind::Delins {
            format!("{ref_first}{start}delins{alt}")
        } else {
            let (s, e) = (start.min(end), start.max(end));
            format!("{ref_first}{s}_{ref_last}{e}{}{alt}", kind.as_str())
        }
    } else if kind == HgvspKind::Fs {
        if alt == b"Ter" {
            format!("{}{start}{}", text(ref_pep), text(alt))
        } else {
            let aa_til_stop = hgvsp_stop_loss_extra_aa(ev, start - 1, true);
            if ref_pep == b"Ter" {
                // A frameshift whose first changed residue is the stop codon is an
                // extension (TranscriptVariationAllele.pm 2007-2010), and the count
                // to the new stop leaves out the replaced stop itself. A stop that
                // the peptide holds before the change ends that translation before
                // the lost stop; the count is then read from the lost stop on.
                let count = aa_til_stop
                    .map(|extra| extra - 1)
                    .or_else(|| ev.residues_from_a_stop_a_frameshift_changes_first(start))
                    .map_or("?".to_string(), |count| count.to_string());
                format!("{}{start}{}extTer{count}", text(ref_pep), text(alt))
            } else {
                // A stop codon before the change that the peptide reads as a
                // residue, or holds as a stop and reads on past, ends that
                // translation before the first changed residue; the count is
                // then read from that residue on.
                let count = aa_til_stop
                    .or_else(|| {
                        ev.residues_to_a_new_stop_past_stops_inside_the_coding_sequence(start)
                    })
                    .or_else(|| hgvsp_frameshift_count_from_residue_1(ev, start))
                    .map_or("?".to_string(), |extra| extra.to_string());
                format!("{}{start}{}fsTer{count}", text(ref_pep), text(alt))
            }
        }
    } else if kind == HgvspKind::Del {
        if n.cut_short && preds.cut_short_names_its_first_residue {
            format!("{}{start}fs", first3(ref_pep))
        } else if ref_pep.len() > 3 {
            format!("{}{start}_{}{end}del", first3(ref_pep), last3(ref_pep))
        } else {
            format!("{}{start}del", text(ref_pep))
        }
    } else if start != end {
        format!("{}{start}_{}{end}", text(ref_pep), text(alt))
    } else {
        format!("{}{start}{}", text(ref_pep), text(alt))
    }
}

impl<'a> PerlCodingEval<'a> {
    fn new(
        variant: &'a InputVariant,
        transcript: &'a Transcript,
        span_start: u64,
        span_end: u64,
        fasta: Option<&'a vep_fasta::IndexedFasta>,
    ) -> Option<Self> {
        let vefc = transcript.vefc.as_ref()?;
        let cds = vefc.translateable_seq.as_ref()?.as_bytes();
        if cds.is_empty() {
            return None;
        }
        let mapper = vefc.mapper.as_ref()?;
        let span = perl_span(transcript, span_start, span_end)?;
        let strand = transcript.strand.as_i8();
        let opposite_strand = variant.strand != transcript.strand;
        let ref_allele: &[u8] = if variant.ref_allele.is_empty() {
            b"-"
        } else {
            &variant.ref_allele
        };
        let alt_allele: &[u8] = {
            let a = variant.alt_allele();
            if a.is_empty() {
                b"-"
            } else {
                a
            }
        };
        // `feature_seq` flips only an allele `seq_is_dna` admits
        // (VariationFeatureOverlapAllele.pm:253, 315-321).
        let feature_seq = if opposite_strand && seq_is_dna(alt_allele) {
            reverse_complement(alt_allele)
        } else {
            alt_allele.to_vec()
        };
        let seq_length = if seq_is_dna(alt_allele) {
            Some(if alt_allele == b"-" {
                0
            } else {
                alt_allele.len() as i64
            })
        } else {
            None
        };
        let (vf_start, vf_end) = (span_start as i64, span_end as i64);
        let ref_length = vf_end - vf_start + 1;
        let alt_length = if alt_allele == b"-" {
            0
        } else {
            alt_allele.len() as i64
        };
        let coding_region = match (transcript.coding_region_start, transcript.coding_region_end) {
            (Some(s), Some(e)) if s > 0 && e > 0 => Some((s as i64, e as i64)),
            _ => None,
        };
        let frameshift_introns: Vec<(i64, i64)> = transcript
            .introns
            .iter()
            .filter(|i| (i.end as i64 - i.start as i64).abs() <= 12)
            .map(|i| (i.start as i64, i.end as i64))
            .collect();
        let table = codon_table_for(transcript);
        let pep_seq = transcript_peptide(vefc.peptide.as_deref(), cds, table);
        let start_codon_stands = start_codon_opens(transcript, &pep_seq);
        // A cache without UTR strings still has the UTR lengths: `N` placeholders
        // keep Ensembl VEP's edit anchors (`cdna_start - 1` into 5'UTR + CDS) and window
        // extents in place, translating to `X` where Ensembl VEP reads real bases. The
        // lengths come from the cache, whose loader refuses a transcript whose
        // mapper asks for more than the transcript can hold.
        let padding = |bases: u64| vec![b'N'; bases as usize];
        let utr5 = compute_five_prime_utr_sequence(transcript, fasta).or_else(|| {
            (mapper.cdna_coding_start > 1).then(|| padding(mapper.cdna_coding_start - 1))
        });
        let utr3 = compute_three_prime_utr_sequence(transcript, fasta).or_else(|| {
            let cdna_len = mapper
                .exon_coord_mapper
                .pairs
                .iter()
                .map(|p| p.from_end)
                .max()
                .unwrap_or(0);
            (cdna_len > mapper.cdna_coding_end).then(|| padding(cdna_len - mapper.cdna_coding_end))
        });
        Some(Self {
            transcript,
            cds,
            pep_seq,
            seq_edits: &vefc.seq_edits,
            utr5,
            utr3,
            table,
            span,
            strand,
            opposite_strand,
            reference_row: variant.reference_allele_row,
            start_codon_stands,
            stop_codon: stop_codon_at_the_coding_end(transcript),
            cdna_coding_start: mapper.cdna_coding_start as i64,
            cdna_coding_end: mapper.cdna_coding_end as i64,
            coding_region,
            frameshift_introns,
            vf_start,
            vf_end,
            ref_allele,
            alt_allele,
            feature_seq,
            chr: &variant.chr,
            fasta,
            seq_length,
            snp: alt_length == ref_length,
            insertion: ref_length < alt_length,
            deletion: ref_length > alt_length,
            codon_memo: [None, None],
            peptide_memo: [None, None],
            peps_memo: None,
            overlaps_stop_memo: None,
            overlaps_stop_cil_memo: None,
            overlaps_start_memo: None,
            ins_del_stop_altered_memo: None,
            edited_stop_codon_kept_cil_memo: None,
            removed_stop_codon_memo: None,
            ins_del_start_altered_memo: None,
            inv_start_altered_memo: None,
            snp_start_altered_memo: None,
            start_lost_memo: None,
            stop_lost_memo: None,
            stop_retained_memo: None,
            stop_gained_memo: None,
            reads_unknown_bases: false,
            reads_a_kept_stop_codon: false,
        })
    }

    /// `_bvfo_preds` `coding`: [`coding_span_pred`] on the evaluation's span, and
    /// the span projects to at least one CDS segment (a lone gap counts, as in
    /// Ensembl VEP).
    fn coding_pred(&self, transcript: &Transcript) -> bool {
        coding_span_pred(transcript, self.vf_start, self.vf_end) && !self.span.cds_coords.is_empty()
    }
}

/// The tests of `_bvfo_preds` `coding` that read the span and the transcript
/// alone: the span is within the feature and overlaps the coding region and an
/// exon. `within_feature` is read from the coordinates as given, so an insertion
/// between the transcript's first base and the base before it (`end < start`)
/// is outside the feature; the tests after it read the ordered span.
fn coding_span_pred(transcript: &Transcript, vf_start: i64, vf_end: i64) -> bool {
    let (crs, cre) = match (transcript.coding_region_start, transcript.coding_region_end) {
        (Some(s), Some(e)) if s > 0 && e > 0 => (s as i64, e as i64),
        _ => return false,
    };
    if !perl_overlap(
        vf_start,
        vf_end,
        transcript.start as i64,
        transcript.end as i64,
    ) {
        return false;
    }
    let (lo, hi) = (vf_start.min(vf_end), vf_start.max(vf_end));
    if !perl_overlap(lo, hi, crs, cre) {
        return false;
    }
    let Some(mapper) = transcript.vefc.as_ref().and_then(|v| v.mapper.as_ref()) else {
        return false;
    };
    mapper
        .exon_coord_mapper
        .pairs
        .iter()
        .any(|p| perl_overlap(lo, hi, p.to_start as i64, p.to_end as i64))
}

/// Whether an allele annotated at `start..=end` lies on the transcript, inside
/// its coding region, and overlaps none of its exons, so `_bvfo_preds` reads it
/// `non_coding` (BaseVariationFeatureOverlapAllele.pm:535-537); outside the
/// coding region it reads `utr` (:523-525). `within_feature` is read from the
/// coordinates as given, the coding region and the exons from the ordered span,
/// which for an insertion holds both flanking bases.
pub(crate) fn annotated_in_coding_intron(transcript: &Transcript, start: u64, end: u64) -> bool {
    let (crs, cre) = match (transcript.coding_region_start, transcript.coding_region_end) {
        (Some(s), Some(e)) if s > 0 && e > 0 => (s as i64, e as i64),
        _ => return false,
    };
    let (vf_start, vf_end) = (start as i64, end as i64);
    if !perl_overlap(
        vf_start,
        vf_end,
        transcript.start as i64,
        transcript.end as i64,
    ) {
        return false;
    }
    let (lo, hi) = (vf_start.min(vf_end), vf_start.max(vf_end));
    if !perl_overlap(lo, hi, crs, cre) {
        return false;
    }
    let Some(mapper) = transcript.vefc.as_ref().and_then(|v| v.mapper.as_ref()) else {
        return false;
    };
    !mapper
        .exon_coord_mapper
        .pairs
        .iter()
        .any(|p| perl_overlap(lo, hi, p.to_start as i64, p.to_end as i64))
}

/// Whether an allele annotated at `start..=end` lies on the transcript and
/// wholly before its coding region in the transcript's direction, in the 5'
/// UTR or in an intron of it, so `_bvfo_preds` reads it `utr` and not `coding`
/// (BaseVariationFeatureOverlapAllele.pm:523-525). `within_feature` is read
/// from the coordinates as given and the coding region from the ordered span,
/// which for an insertion holds both flanking bases: an insertion between the
/// last base before the coding start and the first coding base overlaps the
/// coding region and is not such an allele.
pub(crate) fn annotated_before_the_coding_start(
    transcript: &Transcript,
    start: u64,
    end: u64,
) -> bool {
    let (crs, cre) = match (transcript.coding_region_start, transcript.coding_region_end) {
        (Some(s), Some(e)) if s > 0 && e > 0 => (s as i64, e as i64),
        _ => return false,
    };
    let (vf_start, vf_end) = (start as i64, end as i64);
    if !perl_overlap(
        vf_start,
        vf_end,
        transcript.start as i64,
        transcript.end as i64,
    ) {
        return false;
    }
    let (lo, hi) = (vf_start.min(vf_end), vf_start.max(vf_end));
    match transcript.strand {
        Strand::Forward => hi < crs,
        Strand::Reverse => lo > cre,
    }
}

/// Whether `allele`, an insertion as the 3' shift left it, lies between an
/// intron's last base and the first base of the exon after it, in the
/// transcript's direction.
///
/// HGVS describes a change there by the intron's last positions (a duplication
/// ending on them) or by the two positions that flank the junction (any other
/// insertion), and its recommendations on splicing give such a change no
/// definite protein consequence. A record written between those two bases is
/// annotated on the exon: in a form
/// [`intron_end_insertion_has_no_protein_notation`] names it has no protein
/// notation either, and otherwise it keeps the notation read from the exon.
fn inserted_at_the_end_of_an_intron(allele: &InputVariant, transcript: &Transcript) -> bool {
    if allele.start != allele.end + 1 {
        return false;
    }
    let Some(mapper) = transcript.vefc.as_ref().and_then(|v| v.mapper.as_ref()) else {
        return false;
    };
    let pairs = &mapper.exon_coord_mapper.pairs;
    let exonic = |pos: u64| pairs.iter().any(|p| p.to_start <= pos && pos <= p.to_end);
    match transcript.strand {
        Strand::Forward => pairs.iter().any(|p| p.to_start == allele.start) && !exonic(allele.end),
        Strand::Reverse => pairs.iter().any(|p| p.to_end == allele.end) && !exonic(allele.start),
    }
}

/// Whether the 3' shift that carries `written`, an insertion, to a span that
/// starts at `shifted_start` passes a base of the allele outside A, C, G and
/// T. The shift moves the allele one place while the base it passes, the
/// allele's next on the transcript's strand, equals the reference base there
/// (`shift_indel_3prime` of the `consequences` module). A base outside A, C,
/// G and T that it passes equals such a base of the reference, which no
/// completion of the allele equals: the shift of each completion ends before
/// that base.
fn shift_passes_an_unknown_base(
    written: &InputVariant,
    shifted_start: u64,
    strand: Strand,
) -> bool {
    let Some(inserted) = crate::consequences::shift_motif(written) else {
        return false;
    };
    let moved = usize::try_from(shifted_start.abs_diff(written.start)).unwrap_or(usize::MAX);
    let passed = moved.min(inserted.len());
    let passed = match strand {
        Strand::Forward => &inserted[..passed],
        Strand::Reverse => &inserted[inserted.len() - passed..],
    };
    !passed.is_empty() && !is_perl_unambiguous_dna(passed)
}

/// Whether the transcript's coding sequence ends on its stop codon: the
/// transcript does not carry the `cds_end_NF` flag, and it is a whole number of
/// codons, the last of them a stop codon of the transcript's codon table. It is
/// the answer `StopCodonAtTheCodingEnd::Whole` of
/// `stop_codon_at_the_coding_end`, the rule the stop codon's terms are written
/// by. The codons are the triplets from the first base, which the start exon's
/// phase pads into frame (`Transcript.pm` `translateable_seq`), so a coding
/// sequence that ends inside a codon does not end on a stop codon, whatever its
/// last three bases spell.
///
/// The `cds_start_NF` flag does not enter the test: a coding sequence that
/// lacks its start is read from its first base, and its residues are numbered
/// on those codons, as for an allele given in an exon (`TranscriptMapper.pm`
/// `genomic2pep` adds the start phase).
///
/// On a coding sequence that lacks its end or its closing stop codon, an allele
/// the 3' shift carries out of an intron has the protein notation its whole
/// codons give (`held_codons_give_the_protein_notation`). Where they give none
/// it keeps the empty HGVSp Ensembl VEP writes for it, unless it changes no
/// residue before the stop codon (`changes_no_residue_before_the_stop_codon`).
///
/// A row of `crate::consequences::calculate_consequences` on a transcript
/// whose record shows that polyadenylation completes its stop codon asks this
/// of the completed transcript (`completed_by_polyadenylation`), whose coding
/// sequence ends on that stop codon; the cached record ends inside the codon.
fn coding_sequence_ends_on_its_stop_codon(transcript: &Transcript) -> bool {
    stop_codon_at_the_coding_end(transcript) == StopCodonAtTheCodingEnd::Whole
}

/// Whether `allele`, read at `start..=end`, changes no residue before the stop
/// codon: the transcript's coding sequence is a whole number of codons, the last
/// of them a stop codon, and translates to the same residues before it with the
/// allele applied. The first residue such an allele changes is the stop codon or
/// lies past it. `false` on a coding sequence that no stop codon closes, and on
/// one whose length is not a multiple of three, whatever its last three bases
/// spell: the codons are the triplets from the first base, which the start
/// exon's phase pads into frame (`"N" x $start_phase`, Transcript.pm
/// `translateable_seq`).
fn changes_no_residue_before_the_stop_codon(
    allele: &InputVariant,
    start: u64,
    end: u64,
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> bool {
    let Some(ev) = PerlCodingEval::new(allele, transcript, start, end, fasta) else {
        return false;
    };
    let Some(alternate) = ev.alternate_cds() else {
        return false;
    };
    if ev.cds.len() % 3 != 0 {
        return false;
    }
    let reference = perl_translate(ev.cds, ev.table);
    let Some((&b'*', before_the_stop)) = reference.split_last() else {
        return false;
    };
    perl_translate(&alternate, ev.table).starts_with(before_the_stop)
}

/// Whether the codons the transcript holds give the protein notation of
/// `allele`, read at `start..=end`, whatever follows the coding sequence's last
/// whole codon: a residue of those codons differs once the allele is applied.
/// The residues are those of the whole codons from the first base, which the
/// start exon's phase pads into frame.
///
/// A change of other than whole codons is described from the first residue that
/// differs, read on through the 3' UTR: both that residue and the one that
/// replaces it are whole codons of known bases here. The count to the new stop
/// codon is read on the same bases, and a frame that meets none before the
/// transcript's last base is written with the unknown count.
///
/// A change of whole codons is read on the coding sequence alone. Residues it
/// replaces differ, and so does the residue that follows deleted, duplicated
/// or inserted ones. Where no residue of the sequence follows them none
/// differs: the residues past the sequence's end decide where HGVS describes
/// such a change.
///
/// `false` where no residue of the coding sequence's whole codons differs: the
/// change begins in a codon the sequence does not complete, or past it.
fn held_codons_give_the_protein_notation(
    allele: &InputVariant,
    start: u64,
    end: u64,
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> bool {
    let Some(ev) = PerlCodingEval::new(allele, transcript, start, end, fasta) else {
        return false;
    };
    let Some(alternate) = ev.alternate_cds() else {
        return false;
    };
    let after_the_coding_sequence = ev.utr3.as_ref().map_or(0, Vec::len);
    let Some(alternate_coding) = alternate.len().checked_sub(after_the_coding_sequence) else {
        return false;
    };
    let whole_codons_changed = alternate_coding.abs_diff(ev.cds.len()).is_multiple_of(3);
    let read = if whole_codons_changed {
        &alternate[..alternate_coding]
    } else {
        &alternate[..]
    };
    let residues = perl_translate(read, ev.table);
    perl_translate(ev.cds, ev.table)
        .iter()
        .zip(&residues)
        .position(|(before, after)| before != after)
        .is_some_and(|first| residues[first] != b'X')
}

/// Whether the genomic position of an SNV maps into the transcript's first
/// codon: cDNA `cdna_coding_start ..= cdna_coding_start + 2`.
///
/// This is `_overlaps_start_codon` without its `cds_start_NF` gate, so that a
/// start-codon SNV on a 5'-incomplete transcript is also evaluated by the port,
/// which then applies the gate and reads `missense_variant` where a start-codon
/// verdict would otherwise be assumed.
fn snv_in_first_codon(transcript: &Transcript, pos: u64) -> bool {
    match (transcript.coding_region_start, transcript.coding_region_end) {
        (Some(crs), Some(cre)) if crs <= pos && pos <= cre => {}
        _ => return false,
    }
    let Some(mapper) = transcript.vefc.as_ref().and_then(|v| v.mapper.as_ref()) else {
        return false;
    };
    let ccs = mapper.cdna_coding_start;
    if ccs == 0 {
        return false;
    }
    for p in &mapper.exon_coord_mapper.pairs {
        if p.to_start <= pos && pos <= p.to_end {
            let cdna = if p.ori == 1 {
                p.from_start + (pos - p.to_start)
            } else {
                p.from_start + (p.to_end - pos)
            };
            return ccs <= cdna && cdna <= ccs + 2;
        }
    }
    false
}

/// The coding SO terms Ensembl VEP assigns to an allele on `transcript`: every
/// allele but a one-base pair over a one-base span, and such a pair whose
/// position lies in the first codon.
///
/// `None` when the allele is a one-base pair over a one-base span outside the
/// first codon, the transcript has no coding model, or Ensembl VEP's `coding`
/// pre-predicate is false (the caller keeps its own terms). `Some` is the
/// complete coding-term set, possibly empty, and replaces every term for which
/// `is_perl_coding_term` holds.
///
/// A start-codon SNV takes Ensembl VEP's routes as written: `_inv_start_altered`
/// reads `start_lost` from the edited codon not being `ATG` whenever the
/// transcript has a 5' UTR, the peptide route from `translation_start == 1`
/// and a changed residue otherwise, and `start_retained_variant` from the
/// edited codon being the literal `ATG`, on every codon table.
///
/// Where Ensembl VEP emits `start_lost` together with `start_retained_variant`, only
/// `start_retained_variant` is returned. Ensembl VEP reaches that pair only when
/// `_ins_del_start_altered` is false, which certifies that the edited 5'UTR and
/// CDS still begin `ATG` at the CDS start or that the CDS is unchanged, so the
/// start codon is intact and the `start_lost` of the codon-window peptide route
/// is the erroneous member. vep-rs keeps that divergence.
pub fn perl_coding_terms(
    variant: &InputVariant,
    transcript: &Transcript,
    span_start: u64,
    span_end: u64,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<Vec<Consequence>> {
    let ref_len = if variant.ref_allele == b"-" || variant.ref_allele.is_empty() {
        0
    } else {
        variant.ref_allele.len()
    };
    let alt = variant.alt_allele();
    let alt_len = if alt == b"-" || alt.is_empty() {
        0
    } else {
        alt.len()
    };
    // Ensembl VEP reads the frameshift and the stop-codon guard from the span's CDS
    // coordinates (VariationEffect.pm:1554-1564), so a one-base pair over a span
    // of another length, which `--dont_skip` keeps after the reference-length
    // test, is evaluated here and not as the SNV the per-endpoint pass read.
    let snv = ref_len == 1 && alt_len == 1 && span_start == span_end;
    if snv && !snv_in_first_codon(transcript, span_start) {
        return None;
    }
    // A span the `coding` pre-predicate refuses yields `None` whether or not the
    // evaluation state could be built, so its span tests run before the
    // construction, which copies both UTRs.
    if !coding_span_pred(transcript, span_start as i64, span_end as i64) {
        return None;
    }
    let mut ev = PerlCodingEval::for_row(variant, transcript, span_start, span_end, fasta)?;
    if ev.span.cds_coords.is_empty() {
        return None;
    }
    let mut terms = ev.terms();
    if terms.contains(&Consequence::StartRetainedVariant) {
        terms.retain(|c| *c != Consequence::StartLost);
    }
    Some(terms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::{
        cds_ending_in, deletion_of_cds, insertion_after_cds, make_single_exon_transcript,
        make_test_transcript,
    };
    use vep_core::coordinate::Strand;

    #[test]
    fn test_overlay_seq_edits_skips_position_one_and_the_partial_terminal_codon() {
        // Cached translation: `M` written over a CTG start, `U` at codon 3, and a
        // resolved `A` for the 2-base tail of an 11 bp CDS.
        let pep_seq = b"MAUA";
        let mut pep = b"LA*".to_vec();
        overlay_seq_edits(&mut pep, pep_seq, 1, 11);
        assert_eq!(pep, b"LAU");
        let mut tail = b"X".to_vec();
        overlay_seq_edits(&mut tail, pep_seq, 4, 11);
        assert_eq!(tail, b"X");
        let mut past_end = b"*".to_vec();
        overlay_seq_edits(&mut past_end, b"MA", 3, 9);
        assert_eq!(past_end, b"*");
    }

    #[test]
    fn test_trim_common_peptide_flanks() {
        let (r, a) = trim_common_peptide_flanks(b"ABCDE", b"ABXYCDE");
        assert_eq!(r, b"");
        assert_eq!(a, b"XY");

        let (r, a) = trim_common_peptide_flanks(b"ABXCDE", b"ABCDE");
        assert_eq!(r, b"X");
        assert_eq!(a, b"");

        let (r, a) = trim_common_peptide_flanks(b"ABC", b"ABC");
        assert_eq!(r, b"");
        assert_eq!(a, b"");
    }

    #[test]
    fn test_format_codon_display() {
        assert_eq!(format_codon_display(b"GCT", 0), "Gct");
        assert_eq!(format_codon_display(b"GCT", 1), "gCt");
        assert_eq!(format_codon_display(b"GCT", 2), "gcT");
    }

    #[test]
    fn test_snv_synonymous() {
        let tx = make_test_transcript();
        // CDS starts at cDNA pos 51 = genomic 25_000_050
        // The translateable_seq starts with "ATGGCT..." (Met, Ala, ...)
        // CDS pos 3 (codon 1, pos 2) = 'G' in ATG
        // GCT -> GCC = Ala -> Ala (synonymous)
        let variant = InputVariant::new(
            "21".into(),
            25_000_055,
            25_000_055,
            b"T".to_vec(),
            b"C".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 6);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.codon_number, 2);
        assert_eq!(cc.ref_amino_acid, b'A'); // GCT = Ala
        assert_eq!(cc.alt_amino_acid, b'A'); // GCC = Ala
        assert!(!cc.is_frameshift);
    }

    #[test]
    fn test_snv_missense() {
        let tx = make_test_transcript();
        // CDS pos 5 = codon 2, pos 1 = GCT[1] = C
        // GCT -> GAT = Ala -> Asp (missense)
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 5);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.codon_number, 2);
        assert_eq!(cc.ref_amino_acid, b'A'); // GCT = Ala
        assert_eq!(cc.alt_amino_acid, b'D'); // GAT = Asp
        assert!(!cc.is_frameshift);
    }

    #[test]
    fn test_snv_stop_gained() {
        let tx = make_test_transcript();
        // A position where changing one base creates a stop:
        // CDS pos 7 = codon 3, pos 0 = GGA[0] = G
        // GGA -> TGA = Gly -> Stop
        let variant = InputVariant::new(
            "21".into(),
            25_000_056,
            25_000_056,
            b"G".to_vec(),
            b"T".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 7);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.codon_number, 3);
        assert_eq!(cc.ref_amino_acid, b'G'); // GGA = Gly
        assert_eq!(cc.alt_amino_acid, b'*'); // TGA = Stop
    }

    #[test]
    fn test_snv_reverse_strand_complement() {
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;
        // CDS pos 5 = codon 2, pos 1 = GCT[1] = C
        // Alt A on reverse strand should complement to T: GCT -> GTT = Ala -> Val
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_054,
            b"C".to_vec(),
            b"A".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 5);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.ref_amino_acid, b'A');
        assert_eq!(cc.alt_amino_acid, b'V');
        assert_eq!(cc.alt_codon, "gTt");
    }

    #[test]
    fn test_frameshift_insertion() {
        let tx = make_test_transcript();
        let variant = InputVariant::new(
            "21".into(),
            25_000_052,
            25_000_052,
            b"-".to_vec(),
            b"AA".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 3);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert!(cc.is_frameshift); // 2bp insertion -> frameshift
    }

    #[test]
    fn test_inframe_insertion() {
        let tx = make_test_transcript();
        let variant = InputVariant::new(
            "21".into(),
            25_000_052,
            25_000_052,
            b"-".to_vec(),
            b"AAA".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 3);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert!(!cc.is_frameshift); // 3bp insertion -> inframe
    }

    #[test]
    fn test_codon_window_pure_insertion_boundary_uses_inserted_peptide_only() {
        let tx = make_test_transcript();
        let variant = InputVariant::new(
            "21".into(),
            25_000_060,
            25_000_059,
            b"-".to_vec(),
            b"ATT".to_vec(),
        );
        let bounds = CdsSpanBounds {
            cds_start: 11,
            cds_end: 10,
            translation_start: 4,
            translation_end: 3,
        };

        let peps = compute_codon_window_peptide_alleles(&variant, &tx, &bounds)
            .expect("pure insertion boundary should yield codon-window peptides");
        assert!(
            peps.0.is_empty(),
            "Pure insertion boundary should have empty ref peptide window, got {:?}",
            peps.0
        );
        assert_eq!(peps.1, vec![b'I']);
        assert_eq!(
            classify_inframe_indel_by_peptide(&variant, &tx, 11, Some(&bounds)),
            Some(InframeIndelKind::Insertion)
        );
    }

    #[test]
    fn test_codon_window_boundary_insertion_with_ascending_peptide_bounds_uses_inserted_peptide_only(
    ) {
        let tx = make_test_transcript();
        let variant = InputVariant::new(
            "21".into(),
            25_000_060,
            25_000_059,
            b"-".to_vec(),
            b"ATT".to_vec(),
        );
        let bounds = CdsSpanBounds {
            cds_start: 11,
            cds_end: 10,
            translation_start: 3,
            translation_end: 4,
        };

        let peps = compute_codon_window_peptide_alleles(&variant, &tx, &bounds)
            .expect("boundary insertion with ascending peptide bounds should yield peptides");
        assert!(
            peps.0.is_empty(),
            "Boundary insertion should have empty ref peptide window, got {:?}",
            peps.0
        );
        assert_eq!(peps.1, vec![b'I']);
        assert_eq!(
            classify_inframe_indel_by_peptide(&variant, &tx, 10, Some(&bounds)),
            Some(InframeIndelKind::Insertion)
        );
    }

    #[test]
    fn test_codon_window_does_not_extend_past_cds_end_even_with_cached_utr() {
        let mut tx = make_test_transcript();
        // A cached 3' UTR is not appended to this codon window; the whole port in
        // `PerlCodingEval::alternate_cds` is the path that reads it.
        if let Some(vefc) = tx.vefc.as_mut() {
            vefc.three_prime_utr = Some("AA".into());
        } else {
            panic!("test transcript is expected to have VEFC");
        }

        // Last coding base in the test transcript is genomic 25_004_299 (cDNA 900 / CDS 850).
        // The codon window for this position extends past CDS end, so without UTR
        // the function should return None.
        let variant = InputVariant::new(
            "21".into(),
            25_004_299,
            25_004_299,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        let bounds = crate::mapper::map_genomic_span_to_cds_bounds(variant.start, variant.end, &tx)
            .expect("span bounds should resolve at CDS end");
        assert_eq!(bounds.cds_start, 850);
        assert_eq!(bounds.cds_end, 850);

        // The codon window needs to extend past CDS len (850 bases), but without
        // UTR append the alt_cds is only 850 bases, so the window should fail.
        let result = compute_codon_window_peptide_alleles(&variant, &tx, &bounds);
        assert!(
            result.is_none(),
            "Codon window does not extend past the CDS end (this window appends no UTR)"
        );
    }

    #[test]
    fn test_mnv_synonymous_within_single_codon() {
        let tx = make_test_transcript();
        // CDS: ATG GCT GGA AAA TTC GAT ...
        // Codon 2: GCT (Ala). CDS pos 5-6 = CT.
        // Change CT→CC: GCT→GCC = Ala→Ala (synonymous)
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_055,
            b"CT".to_vec(),
            b"CC".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 5);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.codon_number, 2);
        assert_eq!(cc.ref_amino_acid, b'A');
        assert_eq!(cc.alt_amino_acid, b'A');
        assert!(!cc.is_frameshift);
    }

    #[test]
    fn test_mnv_missense_within_single_codon() {
        let tx = make_test_transcript();
        // CDS: ATG GCT GGA AAA TTC GAT ...
        // Codon 2: GCT (Ala). CDS pos 5-6 = CT.
        // Change CT→AA: GCT→GAA = Ala→Glu (missense)
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_055,
            b"CT".to_vec(),
            b"AA".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 5);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.codon_number, 2);
        assert_eq!(cc.ref_amino_acid, b'A'); // GCT = Ala
        assert_eq!(cc.alt_amino_acid, b'E'); // GAA = Glu
        assert!(!cc.is_frameshift);
    }

    #[test]
    fn test_mnv_stop_gained_within_single_codon() {
        let tx = make_test_transcript();
        // CDS: ATG GCT GGA AAA TTC GAT ...
        // Codon 3: GGA (Gly). CDS pos 7-8 = GG.
        // Change GG→TA: GGA→TAA = Gly→Stop (stop_gained)
        let variant = InputVariant::new(
            "21".into(),
            25_000_056,
            25_000_057,
            b"GG".to_vec(),
            b"TA".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 7);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.codon_number, 3);
        assert_eq!(cc.ref_amino_acid, b'G'); // GGA = Gly
        assert_eq!(cc.alt_amino_acid, b'*'); // TAA = Stop
        assert!(!cc.is_frameshift);
    }

    #[test]
    fn test_mnv_cross_codon_stop_gained() {
        let tx = make_test_transcript();
        // CDS: ATG GCT GGA AAA TTC GAT ...
        // CDS pos 6 = T (end of codon 2: GCT), pos 7 = G (start of codon 3: GGA)
        // Codon window: codons 2-3 = GCT GGA
        // Change TG→AT at positions 6-7: GCT→GCA (Ala), GGA→TGA (Stop)
        // First diff at codon 3: G→* (stop_gained)
        let variant = InputVariant::new(
            "21".into(),
            25_000_055,
            25_000_056,
            b"TG".to_vec(),
            b"AT".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 6);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.ref_amino_acid, b'G'); // Gly from codon 3
        assert_eq!(cc.alt_amino_acid, b'*'); // TGA = Stop
        assert!(!cc.is_frameshift);
    }

    #[test]
    fn test_mnv_cross_codon_missense() {
        let tx = make_test_transcript();
        // CDS: ATG GCT GGA AAA TTC GAT ...
        // CDS pos 3-4: G (end of ATG) and G (start of GCT)
        // Codon window: codons 1-2 = ATG GCT
        // Change GG→AA: ATG→ATA (Met→Ile), GCT→ACT (Ala→Thr)
        // First diff at codon 1: M→I (missense)
        let variant = InputVariant::new(
            "21".into(),
            25_000_052,
            25_000_053,
            b"GG".to_vec(),
            b"AA".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 3);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.ref_amino_acid, b'M'); // ATG = Met
        assert_eq!(cc.alt_amino_acid, b'I'); // ATA = Ile
        assert!(!cc.is_frameshift);
    }

    #[test]
    fn test_mnv_reverse_strand() {
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;
        // CDS: ATG GCT GGA ... (still in mRNA sense)
        // Codon 2: GCT (Ala). CDS pos 5-6 = CT.
        // Alt allele on genome TT → revcomp = AA → GCT→GAA = Ala→Glu
        //
        // For reverse-strand transcripts, the mapper maps variant.start (lower
        // genomic position) to the highest CDS position in the variant span.
        // So for a 2-base MNV covering CDS positions 5-6, the mapper returns
        // cds_pos=6.  get_mnv_codon_change adjusts backward by ref_len-1 to
        // find the first affected CDS position.
        let variant = InputVariant::new(
            "21".into(),
            25_000_054,
            25_000_055,
            b"CT".to_vec(),
            b"TT".to_vec(),
        );
        let result = get_codon_change(&variant, &tx, 6);
        assert!(result.is_some());
        let cc = result.unwrap();
        assert_eq!(cc.ref_amino_acid, b'A'); // GCT = Ala
        assert_eq!(cc.alt_amino_acid, b'E'); // GAA = Glu (TT→revcomp→AA)
        assert!(!cc.is_frameshift);
    }

    /// Complex inframe indel where alt peptide does not start/end with ref.
    ///
    /// Ensembl VEP's protein_altering_variant evaluates before inframe_insertion on raw
    /// (un-truncated) peptides. If the alt peptide doesn't contain the ref peptide
    /// as a prefix or suffix, Ensembl VEP returns protein_altering, not inframe_insertion.
    #[test]
    fn test_inframe_indel_raw_peptide_mismatch_returns_protein_altering() {
        let tx = make_test_transcript();
        // CDS layout: ATG GCT GGA AAA TTC GAT GCT GCT GCT ...
        // Complex indel at CDS pos 7-9 (codon 3 = GGA): replace 3bp with 9bp (net +6, inframe).
        //
        // ref = "GGA" (G), alt = "TTTGGAGCT" → peptide [F, G, A]
        // ref_pep = [G]. alt starts with F (not G), ends with A (not G) → raw mismatch.
        let variant = InputVariant::new(
            "21".into(),
            25_000_056,
            25_000_058,
            b"GGA".to_vec(),
            b"TTTGGAGCT".to_vec(), // 9bp alt, net +6bp inframe
        );

        let result = classify_inframe_indel_by_peptide(&variant, &tx, 7, None);
        assert_eq!(
            result,
            Some(InframeIndelKind::ProteinAltering),
            "Complex inframe indel with raw peptide mismatch should be ProteinAltering"
        );
    }

    /// The two-codon and one-codon windows of a frameshift insertion agree on
    /// stop_gained.
    ///
    /// The insert `TAGCC` splices in at 0-based CDS index 6 (bounds `cds_start: 7`),
    /// so both windows read `GCT TAG ...` and both carry the gained stop: the
    /// two-codon window (translation 2..3) is `GCTTAGCCGGA` -> `A*PX` against ref
    /// `GCTGGA` -> `AG`; the one-codon window (translation 2..2) is `GCTTAGCC` ->
    /// `A*X` against ref `GCT` -> `A`.
    #[test]
    fn test_frameshift_stop_gained_two_codon_and_one_codon_windows_agree() {
        let tx = make_test_transcript();
        // Insert 5bp (frameshift) between codons 2 and 3: CDS pos 6 = genomic
        // 25_000_055, so the VEP insertion is start=25_000_056, end=25_000_055.
        let variant = InputVariant::new(
            "21".into(),
            25_000_056, // start
            25_000_055, // end (insertion)
            b"-".to_vec(),
            b"TAGCC".to_vec(), // 5bp frameshift whose TAG lands in frame
        );

        // Bounds spanning 2 codons.
        let wide_bounds = CdsSpanBounds {
            cds_start: 7, // between codons 2 and 3
            cds_end: 6,
            translation_start: 2,
            translation_end: 3, // 2 codons
        };

        // Bounds spanning 1 codon.
        let narrow_bounds = CdsSpanBounds {
            translation_start: 2,
            translation_end: 2, // 1 codon
            ..wide_bounds
        };

        // `frameshift_stop_gained_in_codon_window` reads the bounds it is given.
        let result = frameshift_stop_gained_in_codon_window(&variant, &tx, &wide_bounds);
        let narrow_peps = compute_codon_window_peptide_alleles(&variant, &tx, &narrow_bounds);
        if let Some((ref_pep, alt_pep)) = &narrow_peps {
            let narrow_has_gained_stop = alt_pep.contains(&b'*') && !ref_pep.contains(&b'*');
            assert_eq!(
                result, narrow_has_gained_stop,
                "the two-codon window's stop_gained verdict must equal the one-codon window's. \
                 narrow ref={:?} alt={:?}",
                ref_pep, alt_pep
            );
        }
    }

    /// A 40bp frameshift insertion far from the CDS end is not stop_retained.
    ///
    /// Condition 1 of `is_stop_retained_insertion` scans `alt_pep` only up to the
    /// capped scan limit, `ref_pep.len() + net_insertion_nt / 3` (1 + 13 = 14 here),
    /// and condition 2's `near_stop` gate needs `cds_start` within 15 bases of the
    /// CDS end (21 against 850).
    #[test]
    fn test_stop_retained_not_fired_for_frameshift_insertion_far_from_cds_end() {
        use crate::consequences::is_stop_retained_insertion;
        let tx = make_test_transcript();
        // CDS is 850bp. The 40bp frameshift at CDS position 20 (genomic 25_000_069;
        // VEP insertion: start=25_000_070, end=25_000_069) splices in two bases
        // into codon 7 (bounds `cds_start: 21`), so the insert's TGA is out of
        // frame and the 43-byte alt window translates to `ALLLELLLLLLLLLX`, with
        // no stop inside or beyond the scan limit.
        let variant = InputVariant::new(
            "21".into(),
            25_000_070, // start
            25_000_069, // end (insertion)
            b"-".to_vec(),
            b"GCTGCTGCTTGAGCTGCTGCTGCTGCTGCTGCTGCTGCTG".to_vec(), // 40bp, TGA out of frame
        );

        let bounds = CdsSpanBounds {
            cds_start: 21,
            cds_end: 20,
            translation_start: 7,
            translation_end: 7,
        };

        let result = is_stop_retained_insertion(&variant, &tx, &bounds, true);
        assert!(
            !result,
            "stop_retained must not fire for a frameshift insertion at CDS pos 20 of an 850bp CDS"
        );
    }

    /// Verify that stop_retained condition 1 does fire near the terminal stop.
    ///
    /// CDS is 850bp. Last complete codon (283) at CDS 847-849 (1-based).
    /// CDS 848 in exon 3: cDNA 898, exon 3 offset 297, genomic 25_004_297.
    ///
    /// Splice mechanics: bounds cds_start=848, cds_end=847 (insertion).
    /// splice_start = 847 (0-based). alt_cds = cds[0..847] + alt + cds[847..].
    /// Window starts at 0-based 846 (codon_cds_start=847 1-based → 846 0-based).
    /// cds[846] = G (first byte of codon 283 = GCT).
    /// To form GCT(Ala) as first codon, insert must start with "CT".
    /// Insert "CTGCTTGA" → alt window: G+CTGCTTGA+CT = GCT,GCT,TGA,CT → [A,A,*,X],
    /// and the capped scan `ref_pep.len() + 8 / 3 = 3` reaches the `*`.
    #[test]
    fn test_stop_retained_condition1_fires_near_terminal_stop() {
        use crate::consequences::is_stop_retained_insertion;
        let tx = make_test_transcript();
        let variant = InputVariant::new(
            "21".into(),
            25_004_298, // start (insertion in exon 3)
            25_004_297, // end
            b"-".to_vec(),
            b"CTGCTTGA".to_vec(), // completes Ala, then Ala, then stop
        );

        let bounds = CdsSpanBounds {
            cds_start: 848,
            cds_end: 847,
            translation_start: 283,
            translation_end: 283,
        };

        let result = is_stop_retained_insertion(&variant, &tx, &bounds, true);
        assert!(
            result,
            "stop_retained condition 1 should fire near terminal stop (CDS 847, CDS len 850)"
        );
    }

    /// Pure inframe insertion with stop in codon window peptides.
    ///
    /// Ensembl VEP's stop_gained fires for any coding variant (including pure inframe
    /// insertions) where alt_pep has '*' and ref_pep doesn't. This test verifies
    /// that the codon-window peptide alleles show the gained stop, which the
    /// InframeInsertion path in consequences.rs then reads.
    #[test]
    fn test_inframe_insertion_codon_window_shows_stop_gained() {
        let tx = make_test_transcript();
        // Pure insertion of 9bp at CDS position 10 (codon 4: AAA = K).
        // CDS pos 10 = genomic 25_000_059. VEP insertion: start=25_000_060, end=25_000_059.
        //
        // Ref codon window (codon 4): AAA → [K]
        // Bounds `cds_start: 11` splice the insert in at 0-based CDS index 10, one
        // base into codon 4, so the alt window is A + AAATAAGCT + AA =
        // AAA ATA AGC TAA → [K, I, S, *].
        let variant = InputVariant::new(
            "21".into(),
            25_000_060,
            25_000_059,
            b"-".to_vec(),
            b"AAATAAGCT".to_vec(), // 9bp inframe with TAA stop
        );

        let bounds = CdsSpanBounds {
            cds_start: 11,
            cds_end: 10,
            translation_start: 4,
            translation_end: 4,
        };

        let peps = compute_codon_window_peptide_alleles(&variant, &tx, &bounds);
        assert!(peps.is_some(), "codon window peptides should resolve");
        let (ref_pep, alt_pep) = peps.unwrap();
        assert!(
            !ref_pep.contains(&b'*'),
            "ref_pep should not contain stop: {:?}",
            ref_pep
        );
        assert!(
            alt_pep.contains(&b'*'),
            "alt_pep should contain stop from TAA: {:?}",
            alt_pep
        );
    }

    /// Complex inframe indel where alt peptide preserves ref as prefix, which
    /// classifies as Insertion, not ProteinAltering.
    ///
    /// Verifies the peptide containment check passes when ref is a prefix of alt,
    /// so the gate does not over-reject valid inframe insertions.
    #[test]
    fn test_inframe_indel_raw_peptide_match_returns_insertion() {
        let tx = make_test_transcript();
        // Complex indel at CDS pos 7-9 (codon 3 = GGA): replace 3bp with 9bp.
        // ref = "GGA" (G), alt = "GGAGCTGCT" → peptide [G, A, A]
        // ref_pep = [G]. alt starts with G → the peptide containment check passes
        // → Insertion.
        let variant = InputVariant::new(
            "21".into(),
            25_000_056,
            25_000_058,
            b"GGA".to_vec(),
            b"GGAGCTGCT".to_vec(), // 9bp, net +6bp inframe
        );

        let result = classify_inframe_indel_by_peptide(&variant, &tx, 7, None);
        assert_eq!(
            result,
            Some(InframeIndelKind::Insertion),
            "Complex inframe indel where alt starts with ref should be Insertion"
        );
    }

    /// Exercises the capped scan limit in `is_stop_retained_insertion` condition 1.
    /// For a 3bp insertion (net 1 AA) the capped scan covers `ref_pep.len() + 1 = 2`
    /// positions, `alt_pep[..2]`.
    ///
    /// Uses codon 3 (GGA = G, CDS 7-9) with bounds `cds_start: 9, cds_end: 8`, which
    /// splice the insert in at 0-based CDS index 8: the 6-byte alt window
    /// `alt_cds[6..12]` is GG + GCT + A = "GGGCTA" → GGG(G), CTA(L), so no stop is
    /// in reach of the scan or anywhere in the window. The genomic coordinates play
    /// no part in the window.
    #[test]
    fn test_stop_retained_condition1_capped_scan_finds_no_stop_in_three_bp_insertion() {
        use crate::consequences::is_stop_retained_insertion;
        let tx = make_test_transcript();

        // 3bp insertion two bases into codon 3 (GGA = G). Inserted "GCT" spells no
        // stop in any frame the window reads.
        let variant = InputVariant::new(
            "21".into(),
            25_000_008,
            25_000_007,
            b"-".to_vec(),
            b"GCT".to_vec(), // 3bp = 1 codon (Ala), no stop
        );

        let bounds = CdsSpanBounds {
            cds_start: 9,
            cds_end: 8,
            translation_start: 3,
            translation_end: 3,
        };

        // Verify codon window peptides and that is_stop_retained returns false
        // (no stop codon in the insertion contribution).
        let peps = compute_codon_window_peptide_alleles(&variant, &tx, &bounds);
        assert!(peps.is_some(), "codon window peptides should resolve");
        let (ref_pep, alt_pep) = peps.unwrap();
        assert_eq!(ref_pep.len(), 1, "ref_pep should be single AA");

        let result = is_stop_retained_insertion(&variant, &tx, &bounds, true);
        // The 3bp insertion of GCT (Ala) has no stop codon, so condition 1
        // should not fire regardless of whether ref_pep[0] == alt_pep[0].
        assert!(
            !result,
            "stop_retained should NOT fire for 3bp insertion with no stop. \
             ref_pep={:?}, alt_pep={:?}",
            ref_pep, alt_pep
        );
    }

    /// A 9bp insertion two bases into codon 3 (GGA = G) resolves a one-residue ref
    /// window.
    ///
    /// Bounds `cds_start: 9, cds_end: 8` splice "AAATAAGCT" in at 0-based CDS index
    /// 8, so the 12-byte alt window is GG + AAATAAGCT + A = "GGAAATAAGCTA" →
    /// [G, N, K, L]: the TAA is out of frame and the window has no stop.
    #[test]
    fn test_stop_retained_condition1_nine_bp_insertion_resolves_single_residue_ref_window() {
        let tx = make_test_transcript();

        let variant = InputVariant::new(
            "21".into(),
            25_000_008,
            25_000_007,
            b"-".to_vec(),
            b"AAATAAGCT".to_vec(), // 9bp containing TAA stop
        );

        let bounds = CdsSpanBounds {
            cds_start: 9,
            cds_end: 8,
            translation_start: 3,
            translation_end: 3,
        };

        let peps = compute_codon_window_peptide_alleles(&variant, &tx, &bounds);
        assert!(peps.is_some());
        let (ref_pep, _) = peps.unwrap();
        assert_eq!(ref_pep.len(), 1);
    }

    #[test]
    fn test_is_frameshift_cds_aware_pure_insertion_mod3() {
        // Pure insertion: cds_start > cds_end (vf_nt_len=0), 3bp alt → not frameshift
        let variant = vep_core::variant::InputVariant::new(
            "21".into(),
            100,
            99,
            b"-".to_vec(),
            b"ACG".to_vec(),
        );
        let bounds = CdsSpanBounds {
            cds_start: 11,
            cds_end: 10,
            translation_start: 4,
            translation_end: 4,
        };
        assert!(!is_frameshift_cds_aware(&variant, &bounds));
    }

    #[test]
    fn test_is_frameshift_cds_aware_cds_span_non_mod3() {
        // The CDS span covers 2 bases (partial exon overlap), so the frame test
        // reads abs(allele_len - 2) % 3: a 5 bp alt is in frame, a 4 bp alt is not.
        let variant = vep_core::variant::InputVariant::new(
            "21".into(),
            100,
            101,
            b"AT".to_vec(),
            b"ATCGC".to_vec(),
        );
        let bounds = CdsSpanBounds {
            cds_start: 10,
            cds_end: 11,
            translation_start: 4,
            translation_end: 4,
        };
        // vf_nt_len = 11 - 10 + 1 = 2, allele_len = 5, abs(5-2)%3 = 0 → not frameshift
        assert!(!is_frameshift_cds_aware(&variant, &bounds));

        // Different case: vf_nt_len=2, allele_len=4 → abs(4-2)%3=2 → frameshift
        let variant2 = vep_core::variant::InputVariant::new(
            "21".into(),
            100,
            101,
            b"AT".to_vec(),
            b"ATCG".to_vec(),
        );
        assert!(is_frameshift_cds_aware(&variant2, &bounds));
    }

    #[test]
    fn test_is_frameshift_cds_aware_deletion() {
        // Deletion: vf_nt_len=4, alt="-" (len=0) → abs(0-4)%3=1 → frameshift
        let variant = vep_core::variant::InputVariant::new(
            "21".into(),
            100,
            103,
            b"ATCG".to_vec(),
            b"-".to_vec(),
        );
        let bounds = CdsSpanBounds {
            cds_start: 10,
            cds_end: 13,
            translation_start: 4,
            translation_end: 5,
        };
        assert!(is_frameshift_cds_aware(&variant, &bounds));
    }

    /// `stop_gained` returns 0 when `stop_lost` holds (VariationEffect.pm:1221-1222).
    /// On a coding sequence that ends in a sense codon without `cds_end_NF`,
    /// Ensembl VEP reads the last three coding positions as the stop codon
    /// (:1373-1381): five bases inserted inside the last codon read `K/S*X`, the
    /// codon left at the coding end is no stop codon and makes `stop_lost` hold
    /// through `_ins_del_stop_altered` (:1473-1480), and `stop_gained` yields,
    /// `frameshift_variant,stop_lost`. No stop codon stands at that coding
    /// end, so none is lost here and the stop codon the insertion makes is
    /// gained. Before a stop codon that closes the coding sequence the guard is
    /// not reached: the `X` branch of `stop_lost` needs the allele on the stop
    /// codon (`_overlaps_stop_codon`, :1444), which puts `*` into the reference
    /// peptide and fails `stop_gained`'s own test.
    #[test]
    fn stop_gained_is_written_where_no_stop_codon_stands_to_be_lost() {
        let mut cds = String::from("ATG");
        cds.push_str(&"GCT".repeat(20));
        cds.push_str("AAA");
        let tx = make_single_exon_transcript(Strand::Forward, &cds, "GGCGCAGCA");
        let variant = insertion_after_cds(&tx, 64, "GCTAA");
        let mut ev = PerlCodingEval::new(&variant, &tx, variant.start, variant.end, None)
            .expect("the insertion lies in the coding sequence");
        let (ref_pep, alt_pep) = ev.peptide_alleles().expect("both peptides translate");
        assert_eq!(
            (ref_pep.as_slice(), alt_pep.as_slice()),
            (&b"K"[..], &b"S*X"[..])
        );
        assert!(alt_pep.contains(&b'*') && !ref_pep.contains(&b'*'));
        assert!(!ev.stop_lost());
        assert!(ev.stop_gained());
        assert_eq!(
            ev.terms(),
            vec![Consequence::StopGained, Consequence::FrameshiftVariant]
        );
    }

    /// A record given on the reverse strand (an Ensembl-format or region line
    /// whose strand column reads `-1`, Parser/VEP_input.pm:131 and
    /// Parser/Region.pm:183) writes its alleles on that strand, and against a
    /// forward-strand transcript Ensembl VEP reads each of them
    /// reverse-complemented, for every predicate (`feature_seq`,
    /// VariationFeatureOverlapAllele.pm:253-255). Every case below is the
    /// reverse-strand spelling of a forward record whose terms follow from the
    /// coding sequence of `make_test_transcript`, `ATG GCT GGA AAA TTC GAT` then
    /// `GCT` repeated, with CDS position p at genomic 25_000_049 + p and the
    /// exon ending at 25_000_299, under the standard codon table; read as
    /// written, the alleles would give another amino acid or another term.
    #[test]
    fn reverse_strand_input_alleles_revcomp_for_every_predicate() {
        use crate::consequences::{calculate_consequences, EffectsConfig};
        use vep_core::consequence::Consequence::*;

        let tx = std::sync::Arc::new(make_test_transcript());
        let config = EffectsConfig::default();
        struct Case {
            span: (u64, u64),
            /// The alleles as the reverse-strand record writes them.
            reverse: &'static str,
            /// The same alleles on the forward strand.
            forward: &'static str,
            terms: &'static [Consequence],
            amino_acids: &'static str,
        }
        let cases = [
            // codon 2 GCT: GAT is Asp; as written, GTT is Val
            Case {
                span: (25_000_054, 25_000_054),
                reverse: "G/T",
                forward: "C/A",
                terms: &[MissenseVariant],
                amino_acids: "A/D",
            },
            // codon 4 AAA: AAG is Lys; as written, AAC is Asn
            Case {
                span: (25_000_061, 25_000_061),
                reverse: "T/C",
                forward: "A/G",
                terms: &[SynonymousVariant],
                amino_acids: "K",
            },
            // codon 3 GGA: TGA is a stop; as written, AGA is Arg
            Case {
                span: (25_000_056, 25_000_056),
                reverse: "C/A",
                forward: "G/T",
                terms: &[StopGained],
                amino_acids: "G/*",
            },
            // one base inserted between codons 1 and 2
            Case {
                span: (25_000_053, 25_000_052),
                reverse: "-/T",
                forward: "-/A",
                terms: &[FrameshiftVariant],
                amino_acids: "-/X",
            },
            // a Lys codon inserted between codons 2 and 3; as written, Phe
            Case {
                span: (25_000_056, 25_000_055),
                reverse: "-/TTT",
                forward: "-/AAA",
                terms: &[InframeInsertion],
                amino_acids: "-/K",
            },
            // codon 3 deleted
            Case {
                span: (25_000_056, 25_000_058),
                reverse: "TCC/-",
                forward: "GGA/-",
                terms: &[InframeDeletion],
                amino_acids: "G/-",
            },
            // the exon's last base, codon 84 GCT: TCT is Ser; as written, ACT is Thr
            Case {
                span: (25_000_299, 25_000_299),
                reverse: "C/A",
                forward: "G/T",
                terms: &[MissenseVariant, SpliceRegionVariant],
                amino_acids: "A/S",
            },
            // the start codon ATG: TTG is Leu; as written, the codon is unchanged
            Case {
                span: (25_000_050, 25_000_050),
                reverse: "T/A",
                forward: "A/T",
                terms: &[StartLost],
                amino_acids: "M/L",
            },
        ];
        let record = |(start, end): (u64, u64), alleles: &str, strand: Strand| {
            let (reference, alt) = alleles.split_once('/').unwrap();
            let mut v = InputVariant::new(
                "21".into(),
                start,
                end,
                reference.as_bytes().to_vec(),
                alt.as_bytes().to_vec(),
            );
            v.strand = strand;
            v
        };
        for case in &cases {
            let forward = record(case.span, case.forward, Strand::Forward);
            let reverse = record(case.span, case.reverse, Strand::Reverse);
            let f = calculate_consequences(&forward, &tx, &config).unwrap();
            let r = calculate_consequences(&reverse, &tx, &config).unwrap();
            let label = format!("{} -1 at {}", case.reverse, case.span.0);
            assert_eq!(r.consequences.as_slice(), case.terms, "{label}");
            assert_eq!(r.amino_acids.as_deref(), Some(case.amino_acids), "{label}");
            assert_eq!(f.consequences, r.consequences, "{label}");
            assert_eq!(
                (f.codons, f.amino_acids),
                (r.codons, r.amino_acids),
                "{label}"
            );
        }

        // The coding predicates read the record's own strand, not only the copy
        // `calculate_consequences` places on the forward strand: against a
        // reverse-strand transcript the same reverse-strand record is read as
        // written (`feature_seq` flips nothing when the strands agree).
        let missense = record((25_000_054, 25_000_054), "G/T", Strand::Reverse);
        let cc = get_codon_change(&missense, &tx, 5).unwrap();
        assert_eq!((cc.alt_codon.as_str(), cc.alt_amino_acid), ("gAt", b'D'));
        let mut reverse_tx = make_test_transcript();
        reverse_tx.strand = Strand::Reverse;
        let cc = get_codon_change(&missense, &reverse_tx, 5).unwrap();
        assert_eq!((cc.alt_codon.as_str(), cc.alt_amino_acid), ("gTt", b'V'));
        let start_lost = record((25_000_050, 25_000_050), "T/A", Strand::Reverse);
        assert_eq!(
            perl_coding_terms(&start_lost, &tx, 25_000_050, 25_000_050, None),
            Some(vec![StartLost])
        );
        // Both edits start at CDS position 7, the first base of codon 3.
        let deletion = record((25_000_056, 25_000_058), "TCC/-", Strand::Reverse);
        assert_eq!(
            classify_inframe_indel_by_peptide(&deletion, &tx, 7, None),
            Some(InframeIndelKind::Deletion)
        );
        let insertion = record((25_000_056, 25_000_055), "-/TTT", Strand::Reverse);
        assert_eq!(
            compute_peptide_alleles(&insertion, &tx, 7),
            Some((b"".to_vec(), b"K".to_vec()))
        );
    }

    /// `feature_seq` (VariationFeatureOverlapAllele.pm:246-264): no bases for
    /// `-`, the allele as given on the transcript's strand, its reverse
    /// complement on the other strand when `seq_is_dna` admits it (`:253`,
    /// `:315-321`), and an allele the letter set rejects as given.
    #[test]
    fn feature_seq_follows_the_strands_and_the_dna_gate() {
        use Strand::{Forward, Reverse};
        assert_eq!(&*feature_seq(b"-", Reverse, Forward), b"");
        assert_eq!(&*feature_seq(b"ACG", Forward, Forward), b"ACG");
        assert_eq!(&*feature_seq(b"ACG", Reverse, Reverse), b"ACG");
        assert_eq!(&*feature_seq(b"ACG", Reverse, Forward), b"CGT");
        assert_eq!(&*feature_seq(b"ar", Forward, Reverse), b"yt");
        assert_eq!(&*feature_seq(b"<DEL>", Reverse, Forward), b"<DEL>");
    }

    /// [`rotate_shifted_allele`] on `GATC` written on `allele_strand`, returning
    /// the rotated bases and the left-rotation steps applied.
    fn rotated(shift: u64, allele_strand: Strand, transcript_strand: Strand) -> (String, usize) {
        let mut allele = b"GATC".to_vec();
        let steps = rotate_shifted_allele(&mut allele, shift, allele_strand, transcript_strand);
        (String::from_utf8(allele).unwrap(), steps)
    }

    /// The four strand pairs: the allele written on the transcript's strand
    /// (forward on forward, reverse on reverse) and on the other strand.
    const STRAND_PAIRS: [(Strand, Strand); 4] = [
        (Strand::Forward, Strand::Forward),
        (Strand::Reverse, Strand::Reverse),
        (Strand::Forward, Strand::Reverse),
        (Strand::Reverse, Strand::Forward),
    ];

    /// A shift of zero leaves the allele as written on either strand, as
    /// `shift_feature_seqs` runs its loop no times (TranscriptVariationAllele.pm:1339),
    /// and so do an empty allele and a single base, which no rotation changes.
    #[test]
    fn a_zero_shift_rotates_nothing_on_either_strand() {
        use Strand::{Forward, Reverse};
        for (allele_strand, transcript_strand) in STRAND_PAIRS {
            assert_eq!(
                rotated(0, allele_strand, transcript_strand),
                ("GATC".to_string(), 0)
            );
        }
        assert_eq!(rotate_shifted_allele(&mut [], 3, Forward, Forward), 0);
        let mut base = b"A".to_vec();
        assert_eq!(rotate_shifted_allele(&mut base, 3, Forward, Reverse), 0);
        assert_eq!(base, b"A");
    }

    /// A shift of the allele's length brings every base back to its place on
    /// either strand: along a forward-strand transcript the loop runs the length
    /// (TranscriptVariationAllele.pm:1339-1342), along a reverse-strand one
    /// `seq_length - shift_length` is zero (`:1338`).
    #[test]
    fn a_shift_of_the_alleles_length_leaves_it_as_written() {
        for (allele_strand, transcript_strand) in STRAND_PAIRS {
            assert_eq!(
                rotated(4, allele_strand, transcript_strand),
                ("GATC".to_string(), 0)
            );
        }
    }

    /// A shift of a multiple of the allele's length leaves it as written on
    /// either strand: the forward-strand loop wraps the allele whole
    /// (TranscriptVariationAllele.pm:1339-1342) and the reverse-strand count is
    /// negative (`:1338`), a rotation of nothing that here is the right one.
    #[test]
    fn a_shift_of_a_multiple_of_the_alleles_length_leaves_it_as_written() {
        for shift in [8, 12] {
            for (allele_strand, transcript_strand) in STRAND_PAIRS {
                assert_eq!(
                    rotated(shift, allele_strand, transcript_strand),
                    ("GATC".to_string(), 0)
                );
            }
        }
    }

    /// A shift longer than the allele and not a multiple of its length rotates
    /// by the remainder, as a shift of the remainder does, read on the
    /// transcript's strand: `GATC` moved five bases is `ATCG` when written on
    /// the transcript's strand and `CGAT` when written on the other, the bases a
    /// one-base shift gives. Ensembl VEP 116.3 counts the places from the
    /// transcript's strand alone (TranscriptVariationAllele.pm:1337-1338): it
    /// agrees on a forward-written allele along a forward-strand transcript,
    /// where the loop wraps, leaves a forward-written allele along a
    /// reverse-strand transcript as written, its count `seq_length -
    /// shift_length` being negative (`:1339`), and turns an allele written on
    /// the reverse strand the wrong way at any shift; each splices bases the
    /// variant never produces.
    #[test]
    fn a_shift_longer_than_the_allele_rotates_by_the_remainder_on_both_strands() {
        use Strand::{Forward, Reverse};
        for (allele_strand, transcript_strand) in [(Forward, Forward), (Reverse, Reverse)] {
            assert_eq!(
                rotated(1, allele_strand, transcript_strand),
                ("ATCG".to_string(), 1)
            );
            assert_eq!(
                rotated(5, allele_strand, transcript_strand),
                ("ATCG".to_string(), 1)
            );
            assert_eq!(
                rotated(6, allele_strand, transcript_strand),
                ("TCGA".to_string(), 2)
            );
            assert_eq!(
                rotated(11, allele_strand, transcript_strand),
                ("CGAT".to_string(), 3)
            );
        }
        for (allele_strand, transcript_strand) in [(Forward, Reverse), (Reverse, Forward)] {
            assert_eq!(
                rotated(1, allele_strand, transcript_strand),
                ("CGAT".to_string(), 3)
            );
            assert_eq!(
                rotated(5, allele_strand, transcript_strand),
                ("CGAT".to_string(), 3)
            );
            assert_eq!(
                rotated(6, allele_strand, transcript_strand),
                ("TCGA".to_string(), 2)
            );
            assert_eq!(
                rotated(11, allele_strand, transcript_strand),
                ("ATCG".to_string(), 1)
            );
        }
    }

    /// The per-transcript evaluation reads an opposite-strand allele through
    /// Ensembl's `reverse_comp`, IUPAC letters included (`_get_alternate_cds`,
    /// TranscriptVariationAllele.pm:2413-2417): `A/R` written on the reverse
    /// strand at the T of codon 2, GCT, gives the codon `GCY`, and no peptide for
    /// an ambiguous allele (`peptide`, :701).
    #[test]
    fn perl_coding_eval_reads_an_opposite_strand_allele_through_the_iupac_complement() {
        let tx = make_test_transcript();
        let mut given = InputVariant::new(
            "21".into(),
            25_000_055,
            25_000_055,
            b"A".to_vec(),
            b"R".to_vec(),
        );
        given.strand = Strand::Reverse;
        let (ref_codon, alt_codon, _, alt_pep) =
            perl_codon_peptides(&given, &tx, 25_000_055, 25_000_055, None).unwrap();
        assert_eq!(
            (ref_codon.as_deref(), alt_codon.as_deref()),
            (Some(&b"GCT"[..]), Some(&b"GCY"[..]))
        );
        assert_eq!(alt_pep, None);
    }

    /// A span the `coding` pre-predicate refuses yields `None` although the
    /// evaluation state builds for it: a deletion inside intron 1 (inside the
    /// genomic coding region, in no exon), one inside the 5' UTR and one before
    /// the transcript. A deletion inside the coding sequence yields the terms of
    /// the state built for it.
    #[test]
    fn perl_coding_terms_is_none_outside_the_coding_region() {
        let tx = make_test_transcript();
        let deletion = |start: u64, end: u64, bases: &str| {
            InputVariant::new(
                "21".into(),
                start,
                end,
                bases.as_bytes().to_vec(),
                b"-".to_vec(),
            )
        };
        for variant in [
            deletion(25_001_000, 25_001_003, "ACGT"),
            deletion(25_000_010, 25_000_012, "ACG"),
            deletion(24_999_990, 24_999_992, "ACG"),
        ] {
            let state = PerlCodingEval::new(&variant, &tx, variant.start, variant.end, None)
                .expect("the state builds for a span the pre-predicate refuses");
            assert!(!state.coding_pred(&tx), "{}", variant.allele_string);
            assert_eq!(
                perl_coding_terms(&variant, &tx, variant.start, variant.end, None),
                None,
                "{}",
                variant.allele_string
            );
        }

        let coding = deletion(25_000_100, 25_000_102, "TGC");
        let mut state = PerlCodingEval::new(&coding, &tx, coding.start, coding.end, None)
            .expect("the deletion lies in the coding sequence");
        assert!(state.coding_pred(&tx));
        let expected = state.terms();
        assert!(expected.contains(&Consequence::InframeDeletion));
        assert_eq!(
            perl_coding_terms(&coding, &tx, coding.start, coding.end, None),
            Some(expected)
        );
    }

    /// The reference row of `--process_ref_homs` is Ensembl VEP's reference
    /// allele object, whose codon and peptide come from the translateable
    /// sequence whatever bases the record wrote (TranscriptVariationAllele.pm:868,
    /// 877): over codon 2 (`GCT`, Ala) a two- or three-base reference row written
    /// `TA` or `TAA` is `synonymous_variant`; the same bases as an alternate
    /// allele are spliced into the sequence and read `missense_variant` and
    /// `stop_gained`.
    #[test]
    fn a_reference_row_reads_its_codon_from_the_transcript() {
        let tx = make_test_transcript();
        let terms = |bases: &str, reference_row: bool| {
            let end = 25_000_053 + bases.len() as u64 - 1;
            let mut v = InputVariant::new(
                "21".into(),
                25_000_053,
                end,
                bases.as_bytes().to_vec(),
                bases.as_bytes().to_vec(),
            );
            v.reference_allele_row = reference_row;
            perl_coding_terms(&v, &tx, v.start, v.end, None)
        };
        assert_eq!(
            terms("TA", true),
            Some(vec![Consequence::SynonymousVariant])
        );
        assert_eq!(
            terms("TAA", true),
            Some(vec![Consequence::SynonymousVariant])
        );
        assert_eq!(terms("TA", false), Some(vec![Consequence::MissenseVariant]));
        assert_eq!(terms("TAA", false), Some(vec![Consequence::StopGained]));
    }

    /// A reference row over a stop codon is `stop_retained_variant`: the
    /// reference allele's codon is the stop read from the translateable
    /// sequence (TranscriptVariationAllele.pm:868, 877), so `stop_lost` is
    /// false (Utils/VariationEffect.pm:1259, both peptides `*`) and
    /// `stop_retained` is `ref_eq_alt_sequence` (`:1315`), which holds whatever
    /// bases the record wrote. Ensembl VEP 116.3 writes the row of
    /// `21 30282168 30282169 CA/CG` on the reverse-strand ENST00000340345 under
    /// `--process_ref_homs` as `stop_retained_variant` with `*` and the `CG`
    /// allele's as `stop_lost` with `*/S`; the transcript here ends its coding
    /// sequence in the same `TGA`, read on the reverse strand over the two
    /// genomic bases `CA` that spell its `TG`, so the `CG` allele is `*/R`.
    #[test]
    fn a_reference_row_over_the_stop_codon_is_stop_retained() {
        let mut cds = String::from("ATG");
        cds.push_str(&"GCT".repeat(20));
        cds.push_str("TGA");
        let tx = make_single_exon_transcript(Strand::Reverse, &cds, "GGCGCAGCA");
        let (start, end) = (40_000_011, 40_000_012);
        let terms = |reference: &str, alt: &str, reference_row: bool| {
            let mut v = InputVariant::new(
                "21".into(),
                start,
                end,
                reference.as_bytes().to_vec(),
                alt.as_bytes().to_vec(),
            );
            v.reference_allele_row = reference_row;
            perl_coding_terms(&v, &tx, start, end, None)
        };
        assert_eq!(
            terms("CA", "CA", true),
            Some(vec![Consequence::StopRetainedVariant])
        );
        assert_eq!(
            terms("CG", "CG", true),
            Some(vec![Consequence::StopRetainedVariant])
        );
        assert_eq!(terms("CA", "CG", false), Some(vec![Consequence::StopLost]));
    }

    /// The start and stop predicates read the transcript's own table, as
    /// `Transcript::translate` does through `is_start_codon` and `is_ter_codon`:
    /// `TTG` opens a peptide with `M` under tables 1, 5 and 11 and with `L`
    /// under table 3, where it is no start; `ATT` opens with `M` under tables 2,
    /// 5 and 11 and with `I` under table 1; a terminal `TGA` is dropped under
    /// table 1 and kept as `W` under tables 3 and 5, a terminal `AGA` dropped
    /// under table 2 and kept as `S` under table 5; an ambiguous terminal `TAR`
    /// is dropped under table 1 and `TRA` kept as `X` under table 2, and an `X`
    /// base reads as any base, as BioPerl's `iupac_iub` lists it, so `GCX` is
    /// `A`. The cached-peptide readers follow: `GTG` under table 11 is a start,
    /// so a cached `M` over it is the dump's rewrite and no edit, and `ATT`
    /// under table 11 likewise.
    ///
    /// Under every table `Bio::Tools::CodonTable` holds, for each of the 64
    /// codons in its order: a peptide opens with `M` exactly where its first
    /// codon is a start codon of the table or codes `M` (Transcript.pm:2163)
    /// and loses its last residue exactly where its last codon is a stop codon
    /// of the table (Transcript.pm:2158), and the first residue of a cached
    /// peptide is read as `M` exactly where the first codon is a start codon
    /// of the table.
    #[test]
    fn start_codon_predicates_follow_the_transcripts_table() {
        let starts = |codon: &[u8]| -> Vec<u8> {
            [1u8, 2, 3, 5, 11, 12]
                .into_iter()
                .filter(|&table| is_start_codon(codon, table))
                .collect()
        };
        assert_eq!(starts(b"TTG"), [1, 5, 11]);
        assert_eq!(starts(b"CTG"), [1, 11, 12]);
        assert_eq!(starts(b"ATT"), [2, 5, 11]);
        assert_eq!(starts(b"ATA"), [2, 3, 5, 11]);
        assert_eq!(starts(b"GTG"), [2, 5, 11]);
        assert_eq!(starts(b"ATG"), [1, 2, 3, 5, 11, 12]);
        assert_eq!(starts(b"GCT"), []);

        let peptide =
            |cds: &[u8], table: u8| String::from_utf8(perl_transcript_peptide(cds, table)).unwrap();
        assert_eq!(peptide(b"TTGGCTTGA", 1), "MA");
        assert_eq!(peptide(b"TTGGCTTGA", 3), "LAW");
        assert_eq!(peptide(b"TTGGCTTGA", 5), "MAW");
        assert_eq!(peptide(b"TTGGCTTGA", 11), "MA");
        assert_eq!(peptide(b"ATTGCTAGA", 1), "IAR");
        assert_eq!(peptide(b"ATTGCTAGA", 2), "MA");
        assert_eq!(peptide(b"ATTGCTAGA", 5), "MAS");
        assert_eq!(peptide(b"ATTGCTAGA", 11), "MAR");
        assert_eq!(peptide(b"CTTGCTTAA", 1), "LA");
        assert_eq!(peptide(b"CTTGCTTAA", 3), "TA");
        assert_eq!(peptide(b"ATGGCTTAR", 1), "MA");
        assert_eq!(peptide(b"ATGGCTTRA", 1), "MA");
        assert_eq!(peptide(b"ATGGCTTRA", 2), "MAX");
        assert_eq!(peptide(b"ATGGCXTAA", 1), "MA");
        assert_eq!(perl_translate_codon(b"GCX", 1), b'A');
        assert_eq!(perl_translate_codon(b"TAX", 1), b'X');

        assert_eq!(with_initial_met(b"VAW", b"GTGGCTTGG", 11).as_ref(), b"MAW");
        assert_eq!(with_initial_met(b"VAW", b"GTGGCTTGG", 1).as_ref(), b"VAW");
        assert_eq!(initial_met_edit(b"ATTAACGAA", b"MNE", 1), Some(b'M'));
        assert_eq!(initial_met_edit(b"ATTAACGAA", b"MNE", 11), None);
        assert_eq!(initial_met_edit(b"ATTAACGAA", b"MNE", 2), None);

        let bases = *b"TCAG";
        let mut codons: Vec<[u8; 3]> = Vec::with_capacity(64);
        for first in bases {
            for second in bases {
                for third in bases {
                    codons.push([first, second, third]);
                }
            }
        }
        let tables: Vec<u8> = (0..=u8::MAX)
            .filter(|&table| vep_core::codon::known_table(table))
            .collect();
        let filler = codons
            .iter()
            .copied()
            .find(|codon| {
                translate_codon(codon) != b'M'
                    && tables.iter().all(|&table| {
                        translate_codon_with_table(codon, table) == translate_codon(codon)
                            && !is_start_codon(codon, table)
                            && !is_ter_codon(codon, table)
                    })
            })
            .expect("a codon every table reads as one residue");
        let residue = translate_codon(&filler);
        for &table in &tables {
            let stop = codons
                .iter()
                .copied()
                .find(|codon| is_ter_codon(codon, table))
                .expect("a stop codon of the table");
            for codon in &codons {
                let own = translate_codon_with_table(codon, table);
                let start = is_start_codon(codon, table);
                let opened = [*codon, filler, stop].concat();
                let peptide = perl_transcript_peptide(&opened, table);
                assert_eq!(
                    peptide,
                    [if start { b'M' } else { own }, residue],
                    "table {table}"
                );
                assert_eq!(peptide[0] == b'M', start || own == b'M', "table {table}");

                let mut kept = vec![residue, residue];
                if !is_ter_codon(codon, table) {
                    kept.push(own);
                }
                let closed = [filler, filler, *codon].concat();
                assert_eq!(
                    perl_transcript_peptide(&closed, table),
                    kept,
                    "table {table}"
                );

                let cached = [residue, residue];
                assert_eq!(
                    with_initial_met(&cached, &opened, table).as_ref(),
                    [if start { b'M' } else { residue }, residue],
                    "table {table}"
                );
            }
        }
    }

    /// `hgvsp_del_peptides` reads the alternate sequence under the transcript's
    /// table: coding sequence `ATG GCT AAA TAA` with 3' UTR `TGA GCA AGA CCC TAA`
    /// and the stop codon deleted translates on from residue 4 through `TGA`
    /// (Trp under code 2) and `GCA` to the code-2 stop `AGA`, so the alternate
    /// peptide is `TrpAla`; the standard code stops at the `TGA` at once and
    /// leaves it empty. The reference peptide ends before residue 4 either way.
    #[test]
    fn del_peptides_follow_the_transcripts_table() {
        for (table, alt) in [(2u8, &b"TrpAla"[..]), (1, &b""[..])] {
            let mut tx =
                make_single_exon_transcript(Strand::Forward, "ATGGCTAAATAA", "TGAGCAAGACCCTAA");
            tx.vefc.as_mut().unwrap().codon_table = table;
            let variant = deletion_of_cds(&tx, 10, 12, "TAA");
            let ev = PerlCodingEval::new(&variant, &tx, variant.start, variant.end, None).unwrap();
            let mut n = HgvsProteinNotation {
                ref_pep: None,
                alt_pep: None,
                start: 0,
                end: 0,
                kind: HgvspKind::Del,
                original_ref: Vec::new(),
                preseq: Vec::new(),
                cut_short: false,
            };
            assert_eq!(hgvsp_del_peptides(&ev, &mut n), Some(()));
            assert_eq!(n.alt_pep.as_deref(), Some(alt), "table {table}");
            assert_eq!(n.ref_pep.as_deref(), Some(&b""[..]), "table {table}");
            assert_eq!(n.start, 4, "table {table}");
        }
    }

    fn pair(from_start: u64, from_end: u64, to_start: u64, to_end: u64, ori: i8) -> MapperPair {
        MapperPair {
            from_start,
            from_end,
            to_start,
            to_end,
            ori,
        }
    }

    /// The stop codon overlap reads the codon's own genomic bases
    /// (`_overlaps_stop_codon_cil`, VariationEffect.pm:1418-1422, reads the last
    /// three positions of the genomic coding region): on a codon an intron
    /// splits, the stop bases that end the exon before overlap and the intron
    /// bases beside the codon's last do not, on either strand; a whole codon
    /// reads as its three positions, and an insertion's range overlaps it as
    /// `overlap` reads a range (:1424-1427); and of a codon an edit completes,
    /// only the bases a pair holds overlap, not the sense bases before them.
    #[test]
    fn the_stop_codon_overlap_reads_the_codons_own_bases() {
        // Forward: cDNA 1..=100 on 100..=199 and 101..=200 on 300..=399.
        let forward = [pair(1, 100, 100, 199, 1), pair(101, 200, 300, 399, 1)];
        let overlaps = |coding_end: i64, v: (i64, i64)| {
            overlaps_stop_codon_bases(&forward, (51, coding_end), 3, v)
        };
        // The codon on 198, 199 and 300.
        for (v, expected) in [
            ((198, 198), true),
            ((199, 199), true),
            ((300, 300), true),
            ((298, 299), false),
            ((197, 197), false),
            ((301, 301), false),
        ] {
            assert_eq!(overlaps(101, v), expected, "forward 2+1 {v:?}");
        }
        // The codon on 199, 300 and 301.
        for (v, expected) in [((199, 199), true), ((299, 299), false), ((301, 301), true)] {
            assert_eq!(overlaps(102, v), expected, "forward 1+2 {v:?}");
        }
        // A whole codon on 300..=302.
        for (v, expected) in [((300, 302), true), ((299, 299), false), ((199, 199), false)] {
            assert_eq!(overlaps(103, v), expected, "forward whole {v:?}");
        }
        // An insertion's range, its start past its end: between the whole
        // codon's first two bases, as given (301, 300) and with the start moved
        // by one inserted base (302, 300); between its last two, moved (303, 301).
        for (v, expected) in [((301, 300), true), ((302, 300), true), ((303, 301), false)] {
            assert_eq!(overlaps(103, v), expected, "forward insertion {v:?}");
        }
        // The same range between the two stop bases that end the exon, and
        // between the exon's last base and the intron.
        for (v, expected) in [((199, 198), true), ((200, 199), false)] {
            assert_eq!(overlaps(101, v), expected, "forward 2+1 insertion {v:?}");
        }
        // Reverse: cDNA 1..=100 on 399..=300 and 101..=200 on 199..=100.
        let reverse = [pair(1, 100, 300, 399, -1), pair(101, 200, 100, 199, -1)];
        let overlaps = |coding_end: i64, v: (i64, i64)| {
            overlaps_stop_codon_bases(&reverse, (51, coding_end), 3, v)
        };
        // The codon's first two bases on 301 and 300, its last on 199.
        for (v, expected) in [
            ((199, 199), true),
            ((200, 201), false),
            ((300, 300), true),
            ((301, 301), true),
            ((302, 302), false),
        ] {
            assert_eq!(overlaps(101, v), expected, "reverse 2+1 {v:?}");
        }
        // A whole codon on 150..=152.
        for (v, expected) in [((150, 152), true), ((153, 153), false), ((149, 149), false)] {
            assert_eq!(overlaps(150, v), expected, "reverse whole {v:?}");
        }
        // A one-base exon inside the codon: 199, 250 and 300.
        let three = [
            pair(1, 100, 100, 199, 1),
            pair(101, 101, 250, 250, 1),
            pair(102, 200, 300, 398, 1),
        ];
        for (v, expected) in [((250, 250), true), ((199, 199), true), ((198, 198), false)] {
            assert_eq!(
                overlaps_stop_codon_bases(&three, (51, 102), 3, v),
                expected,
                "three exons {v:?}"
            );
        }
        // A coding sequence of two bases reads no base before its first.
        assert!(!overlaps_stop_codon_bases(
            &forward,
            (99, 100),
            3,
            (197, 197)
        ));
        // A stop codon that the last two coding bases begin, on 198 and 199:
        // the last base of the whole codon before them, on 197, is none of it.
        for (v, expected) in [((198, 198), true), ((199, 199), true), ((197, 197), false)] {
            assert_eq!(
                overlaps_stop_codon_bases(&forward, (51, 100), 2, v),
                expected,
                "begun by two bases {v:?}"
            );
        }
        // A codon an edit completes after its first base, the transcript ending
        // on 199: the pair holds cDNA 100 alone of 100..=102.
        let edited = [pair(1, 100, 100, 199, 1)];
        for (v, expected) in [((199, 199), true), ((197, 198), false), ((200, 201), false)] {
            assert_eq!(
                overlaps_stop_codon_bases(&edited, (51, 102), 3, v),
                expected,
                "edited {v:?}"
            );
        }
    }

    /// [`make_test_transcript`] with a 549-base coding sequence whose stop codon
    /// `stop` has its first two bases on the last two of exon 2 and its third on
    /// the first of exon 3, followed by a 3' UTR beginning `utr3_start`.
    fn make_transcript_with_split_stop_codon(stop: &str, utr3_start: &str) -> Transcript {
        let mut tx = make_test_transcript();
        let vefc = tx.vefc.as_mut().unwrap();
        vefc.translateable_seq = Some(format!("ATG{}{stop}", "GCT".repeat(181)));
        vefc.three_prime_utr = Some(format!("{utr3_start}{}", "GCA".repeat(667)));
        let mapper = vefc.mapper.as_mut().unwrap();
        mapper.cdna_coding_start = 53;
        mapper.cdna_coding_end = 601;
        tx.cdna_coding_start = Some(53);
        tx.cdna_coding_end = Some(601);
        tx.coding_region_start = Some(25_000_052);
        tx.coding_region_end = Some(25_004_000);
        tx.translation_start = tx.coding_region_start;
        tx.translation_end = tx.coding_region_end;
        tx
    }

    /// A deletion of the second base of a stop codon an intron splits after it
    /// is tested on the edited sequence: `TAG` before a 3' UTR beginning `A`
    /// leaves `TGA` and is `stop_retained_variant`, before one beginning `C` it
    /// leaves `TGC` and is `stop_lost`. Ensembl VEP's genomic overlap holds the
    /// intron's last two bases in place of the deleted base, so it writes
    /// `coding_sequence_variant` for the first. Intended divergence.
    #[test]
    fn a_deletion_in_a_split_stop_codon_is_tested_on_the_codons_own_bases() {
        let deleted = |tx: &Transcript| {
            let variant = InputVariant::new(
                "21".into(),
                25_002_299,
                25_002_299,
                b"A".to_vec(),
                b"-".to_vec(),
            );
            perl_coding_terms(&variant, tx, variant.start, variant.end, None)
        };
        assert_eq!(
            deleted(&make_transcript_with_split_stop_codon("TAG", "A")),
            Some(vec![Consequence::StopRetainedVariant])
        );
        assert_eq!(
            deleted(&make_transcript_with_split_stop_codon("TAG", "C")),
            Some(vec![Consequence::StopLost])
        );
    }

    /// A single-exon model whose stop codon an edit completes after its first
    /// `held` bases, beside its twin with the whole codon on the genome: the
    /// twin is [`make_single_exon_transcript`] without a 3' UTR, and the edited
    /// model has the same coding sequence and ends `3 - held` genomic bases
    /// earlier, its mapper pair holding the bases the genome has.
    fn make_stop_completing_edit_twins(
        strand: Strand,
        cds: &str,
        held: u64,
    ) -> (Transcript, Transcript) {
        let twin = make_single_exon_transcript(strand, cds, "");
        let mut edited = twin.clone();
        let cut = 3 - held;
        if strand == Strand::Reverse {
            edited.start += cut;
            edited.exons[0].start += cut;
            edited.coding_region_start = edited.coding_region_start.map(|s| s + cut);
            edited.translation_start = edited.coding_region_start;
        } else {
            edited.end -= cut;
            edited.exons[0].end -= cut;
            edited.coding_region_end = edited.coding_region_end.map(|e| e - cut);
            edited.translation_end = edited.coding_region_end;
        }
        let exons = edited.exons.clone();
        let vefc = edited.vefc.as_mut().unwrap();
        vefc.sorted_exons = exons;
        let mapper = vefc.mapper.as_mut().unwrap();
        let mut pairs = mapper.exon_coord_mapper.pairs.clone();
        pairs[0].from_end -= cut;
        if strand == Strand::Reverse {
            pairs[0].to_start += cut;
        } else {
            pairs[0].to_end -= cut;
        }
        mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
        (edited, twin)
    }

    /// An insertion before a stop codon an edit completes reads as on the twin
    /// whose whole codon is on the genome. Ensembl VEP's genomic overlap
    /// (VariationEffect.pm:1418-1422) holds the stop bases the genome has and
    /// the sense bases before them, so on the reverse strand it accepts an
    /// insertion `3 - held` bases further from the codon than on the twin, and
    /// on the forward strand one base inserted before the codon, and writes
    /// `stop_retained_variant` where the twin is a frameshift. Intended
    /// divergence.
    #[test]
    fn an_insertion_before_a_stop_codon_an_edit_completes_reads_as_on_its_whole_twin() {
        let terms = |variant: &InputVariant, tx: &Transcript| {
            perl_coding_terms(variant, tx, variant.start, variant.end, None)
        };
        // Four bases inserted with `TAA` left at the old stop position: the twin
        // accepts the sites up to three bases before the codon.
        let cds = cds_ending_in("AAG", "TAA");
        for held in [1, 2] {
            let (edited, twin) = make_stop_completing_edit_twins(Strand::Reverse, &cds, held);
            for after in 58..=66 {
                assert_eq!(
                    terms(&insertion_after_cds(&edited, after, "CCTC"), &edited),
                    terms(&insertion_after_cds(&twin, after, "CCTC"), &twin),
                    "reverse, held {held}, after {after}"
                );
            }
            assert_eq!(
                terms(&insertion_after_cds(&twin, 62, "CCTC"), &twin),
                Some(vec![Consequence::FrameshiftVariant]),
                "reverse twin, held {held}"
            );
        }
        // The vertebrate mitochondrial code, the stop codon `AGA` held for two
        // bases: `T` inserted before it leaves `TAG` at the old stop position.
        let cds = cds_ending_in("GCT", "AGA");
        let (mut edited, mut twin) = make_stop_completing_edit_twins(Strand::Forward, &cds, 2);
        for tx in [&mut edited, &mut twin] {
            tx.vefc.as_mut().unwrap().codon_table = 2;
        }
        for after in 62..=66 {
            assert_eq!(
                terms(&insertion_after_cds(&edited, after, "T"), &edited),
                terms(&insertion_after_cds(&twin, after, "T"), &twin),
                "forward, held 2, after {after}"
            );
        }
    }

    /// `coding_sequence_ends_on_its_stop_codon` is the answer `Whole` of
    /// `stop_codon_at_the_coding_end` at every coding end: whole codons whose
    /// last is a stop codon of the transcript's table, with or without
    /// `cds_start_NF`, and no other. A last codon with an ambiguous base is a
    /// stop codon where every reading of it is one
    /// (`Bio::Tools::CodonTable::is_ter_codon`).
    #[test]
    fn coding_sequence_ends_on_its_stop_codon_where_the_rule_reads_a_whole_stop_codon() {
        let check = |case: &str, codons: &[&str], table: u8, flags: &[&str], ends: bool| {
            for strand in [Strand::Forward, Strand::Reverse] {
                let mut tx = make_single_exon_transcript(strand, &codons.concat(), "GGCGCAGCA");
                tx.vefc.as_mut().unwrap().codon_table = table;
                tx.flags = flags
                    .iter()
                    .map(|flag| (*flag).to_string())
                    .collect::<Vec<_>>()
                    .into();
                assert_eq!(
                    stop_codon_at_the_coding_end(&tx) == StopCodonAtTheCodingEnd::Whole,
                    ends,
                    "{strand:?}, {case}: the rule"
                );
                assert_eq!(
                    coding_sequence_ends_on_its_stop_codon(&tx),
                    ends,
                    "{strand:?}, {case}: the notation's test"
                );
            }
        };
        let closed = ["ATG", "GCT", "AAA", "TAA"];
        check("closed by a stop codon", &closed, 1, &[], true);
        check(
            "closed by a codon that is a stop codon with either purine",
            &["ATG", "GCT", "AAA", "TAR"],
            1,
            &[],
            true,
        );
        check(
            "closed, without its start",
            &closed,
            1,
            &["cds_start_NF"],
            true,
        );
        check(
            "closed, without its end",
            &closed,
            1,
            &["cds_end_NF"],
            false,
        );
        check(
            "a last codon that names Pro",
            &["ATG", "GCT", "AAA", "CCA"],
            1,
            &[],
            false,
        );
        let closed_under_table_1 = ["ATG", "GCT", "AAA", "TGA"];
        check(
            "a stop codon of table 1",
            &closed_under_table_1,
            1,
            &[],
            true,
        );
        check(
            "the same codon, Trp under table 2",
            &closed_under_table_1,
            2,
            &[],
            false,
        );
        check(
            "one base past whole codons, the last three a stop codon",
            &["ATG", "GCT", "AAA", "C", "TAA"],
            1,
            &[],
            false,
        );
    }

    /// Insertions written in intron 1 of the fixture, three bases before
    /// exon 2, which opens with `CT`, the last two bases of codon 84. Each
    /// is read with the allele a 3' shift of five bases leaves between
    /// codons 84 and 85, where its bases after the first five read `GAA` and
    /// `TAG`. The shift of the first passes known bases, and its unknown
    /// bases follow the stop codon: it has the notation of its completion.
    /// The shift of the last passes its unknown bases, which only unknown
    /// bases of the reference equal: it has none.
    #[test]
    fn an_unknown_base_the_shift_passes_leaves_the_intron_row_without_a_notation() {
        let mut tx = make_test_transcript();
        let cds = tx
            .vefc
            .as_mut()
            .and_then(|vefc| vefc.translateable_seq.as_mut())
            .expect("coding model");
        cds.truncate(cds.len() - cds.len() % 3);
        cds.push_str("TAA");
        let notation = |inserted: &str| {
            let written = InputVariant::new(
                "21".into(),
                25_001_997,
                25_001_996,
                b"-".to_vec(),
                inserted.as_bytes().to_vec(),
            );
            let carried = crate::hgvs::shifted_variant(
                &written,
                Strand::Forward,
                &tx,
                (25_002_002, 25_002_001),
                None,
            )
            .expect("an insertion");
            perl_hgvs_protein(
                &written,
                Some((&carried, carried.start, carried.end)),
                &tx,
                None,
                false,
            )
        };
        assert_eq!(
            notation("CAGCTGAATAGNNNN").as_deref(),
            Some("Ala84_Ala85insGluTer")
        );
        assert_eq!(
            notation("CAGCTGAATAGAAAA").as_deref(),
            Some("Ala84_Ala85insGluTer")
        );
        assert_eq!(notation("NNNCTGAATAGCCAT"), None);
    }

    /// The bases a shift passes are the allele's first on the transcript's
    /// strand: its first as the forward strand writes it along a
    /// forward-strand transcript, its last along a reverse-strand one. A
    /// shift of the allele's length passes every base.
    #[test]
    fn the_shift_passes_the_bases_the_transcript_reads_first() {
        let passes = |inserted: &[u8], shifted_start: u64, strand: Strand| {
            let written =
                InputVariant::new("21".into(), 1_001, 1_000, b"-".to_vec(), inserted.to_vec());
            shift_passes_an_unknown_base(&written, shifted_start, strand)
        };
        assert!(!passes(b"CAGN", 1_004, Strand::Forward));
        assert!(passes(b"CAGN", 1_005, Strand::Forward));
        assert!(passes(b"NCAG", 1_002, Strand::Forward));
        assert!(!passes(b"NCAG", 998, Strand::Reverse));
        assert!(passes(b"NCAG", 997, Strand::Reverse));
        assert!(passes(b"CAGN", 1_000, Strand::Reverse));
        assert!(passes(b"CNAG", 1_009, Strand::Forward));
        assert!(!passes(b"CAGT", 1_009, Strand::Forward));
    }
}
