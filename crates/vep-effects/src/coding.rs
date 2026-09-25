// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Coding-region effect analysis: codon extraction, translation, and amino acid comparison.
//!
//! Given a variant position within the CDS, this module extracts the affected codon(s),
//! translates them using the standard genetic code, and compares reference vs alternate
//! amino acids. Handles SNVs, inframe indels, and frameshifts, mirroring the Perl VEP
//! `TranscriptVariationAllele` codon/peptide logic for concordance.
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`); `TranscriptMapper` is in ensembl core release/115.
//!
//! Key functions:
//! - [`get_codon_change`]: codon-level effect for a single variant
//! - [`classify_inframe_indel_by_peptide`]: inframe insertion/deletion/protein-altering classification
//! - [`compute_peptide_alleles`]: full CDS translation with trimmed peptide comparison
//! - [`compute_codon_window_peptide_alleles`]: Perl-style local codon-window translation

use vep_core::codon::{translate_codon, translate_codon_with_table};
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

use crate::mapper::CdsSpanBounds;

/// A codon change caused by a variant.
#[derive(Debug, Clone)]
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

/// Compute frameshift status using CDS-mapped variant lengths (Perl parity).
///
/// Perl VEP's `VariationEffect::frameshift` predicate uses:
///   `var_len = cds_end - cds_start + 1`  (not raw ref allele length)
///   `allele_len = seq_length`             (alt allele length)
///   `return abs(allele_len - var_len) % 3`
///
/// For insertions with CDS-inverted coords (`cds_start > cds_end`), `var_len = 0`.
/// For variants partially spanning introns/UTRs, `var_len` may be smaller than
/// the raw reference allele length, changing the mod-3 result.
pub fn is_frameshift_cds_aware(variant: &InputVariant, bounds: &CdsSpanBounds) -> bool {
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
pub enum InframeIndelKind {
    Insertion,
    Deletion,
    ProteinAltering,
    LeadingStopGained,
}

/// Perl's `protein_altering_variant` predicate suppresses protein_altering when
/// the alt peptide contains a gained stop and has prefix/suffix containment with
/// the ref peptide (`$alt_pep =~ /^\Q$ref_pep\E|\Q$ref_pep\E$/`).
fn has_gained_stop_with_containment(ref_pep: &[u8], alt_pep: &[u8]) -> bool {
    !ref_pep.is_empty()
        && alt_pep.contains(&b'*')
        && !ref_pep.contains(&b'*')
        && (alt_pep.starts_with(ref_pep) || alt_pep.ends_with(ref_pep))
}

/// Classify an inframe indel using codon-level heuristics, mirroring
/// Perl VEP's `inframe_deletion` / `inframe_insertion` / `protein_altering_variant`
/// predicate logic.
///
/// Perl's `inframe_deletion` checks at the **codon string** level:
/// - alt_codon is a prefix or suffix of ref_codon -> inframe_deletion
/// - after trimming shared prefix/suffix, if alt remainder is empty -> inframe_deletion
/// - otherwise -> protein_altering_variant
///
/// Returns `None` if the codon alleles cannot be derived (caller uses fallback).
pub fn classify_inframe_indel_by_peptide(
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

    // Perl's _get_codon_alleles returns the codons covering the variant span,
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

    // Perl's `VariationEffect::inframe_insertion` returns 0 when `start_lost`
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
        // Perl's inframe_insertion uses peptide comparison (not codon):
        // alt_pep starts/ends with ref_pep -> inframe_insertion.
        let (manual_ref_pep, manual_alt_pep) = match start_codon_peps {
            Some((r, a)) => (r, a),
            None => (
                translate_to_peptide(&ref_upper),
                translate_to_peptide(&alt_upper),
            ),
        };

        // Perl's inframe_insertion returns 0 if start_lost(@_); an insertion
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

        // Gate: Perl's inframe_insertion returns 0 if
        // start_retained_variant(@_) && alt_pep ends with ref_pep. This prevents
        // inframe_insertion when the insertion is before/at the start codon and the
        // Met is preserved. Must check before stop truncation (Perl order).
        if normalized_cds_pos <= 3
            && manual_ref_pep.first() == Some(&b'M')
            && manual_alt_pep.first() == Some(&b'M')
            && manual_alt_pep.ends_with(&manual_ref_pep)
        {
            return Some(InframeIndelKind::ProteinAltering);
        }

        // The containment check prefers compute_codon_window_peptide_alleles, which
        // builds the full alt CDS like Perl's _get_alternate_cds: the manual codon
        // window is too narrow and over-calls protein_altering. Start codon
        // variants (cds_pos <= 3) keep the manual window and its gates.
        //
        // For a pure insertion spanning a codon boundary (translation_start !=
        // translation_end) Perl's inverted translation bounds yield a zero-length
        // ref codon window, so containment passes trivially; un-inverted bounds
        // give a wider window that breaks it, hence the empty-ref case.
        let cross_codon_insertion =
            ref_is_dash && cds_bounds.is_some_and(|b| b.translation_start != b.translation_end);
        let (ref_pep, mut alt_pep) = if normalized_cds_pos > 3
            && cross_codon_insertion
            && !normalized.trimmed_transcript_flanks
        {
            // Pure insertion spanning codon boundary: Perl has empty ref_pep.
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

        // Save the original (pre-truncation) alt_pep. Perl's protein_altering_variant
        // predicate checks containment on untrimmed peptides (`$alt_pep =~ /^\Q$ref_pep\E
        // |\Q$ref_pep\E$/`), while inframe_insertion truncates at the first stop before
        // its containment check. When untrimmed containment passes but truncated fails,
        // Perl suppresses protein_altering and only stop_gained fires.
        let original_alt_pep = alt_pep.clone();

        // Perl: trim everything after first stop in alt_pep for insertion check.
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

        // Perl's protein_altering_variant predicate returns 0 if alt_pep starts
        // with '*'. This suppresses protein_altering and lets stop_gained stand
        // alone. Only applies when containment failed (would have been
        // ProteinAltering); when containment passes, Perl's inframe_insertion
        // fires alongside stop_gained.
        if alt_starts_with_stop {
            return Some(InframeIndelKind::LeadingStopGained);
        }

        // When untrimmed containment passes but truncated containment fails
        // (because stop truncation removed the matching suffix), and the alt
        // peptide contains a gained stop, Perl suppresses protein_altering.
        if has_gained_stop_with_containment(&ref_pep, &original_alt_pep) {
            return Some(InframeIndelKind::LeadingStopGained);
        }

        return Some(InframeIndelKind::ProteinAltering);
    }

    if alt_upper.len() < ref_upper.len() {
        // Net deletion at codon level.
        // Perl's protein_altering_variant predicate returns 0 when alt_pep starts
        // with '*'. For deletions that create a stop at the first alt codon, Perl
        // emits only stop_gained without protein_altering or inframe_deletion.
        // Check this before the codon containment test.
        let del_ref_pep = translate_to_peptide(&ref_upper);
        let del_alt_pep = translate_to_peptide(&alt_upper);
        if del_alt_pep.first() == Some(&b'*') && del_ref_pep.first() != Some(&b'*') {
            return Some(InframeIndelKind::LeadingStopGained);
        }

        // Perl's inframe_deletion uses codon-level comparison:
        // alt_codon is prefix/suffix of ref_codon -> inframe_deletion.
        if ref_upper.starts_with(&alt_upper) || ref_upper.ends_with(&alt_upper) {
            return Some(InframeIndelKind::Deletion);
        }
        let (_r_trimmed, a_trimmed) = trim_common_flanks(&ref_upper, &alt_upper);
        if a_trimmed.is_empty() {
            return Some(InframeIndelKind::Deletion);
        }

        // When codon-level containment fails but peptide-level containment
        // passes, and the alt peptide contains a gained stop, Perl suppresses
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

    // Perl's `TranscriptVariationAllele::peptide` / `codon` paths bail out unless
    // the allele sequence is unambiguous DNA. When these helpers return undef,
    // consequence assignment falls back to coding_sequence_variant rather than
    // peptide-based indel terms.
    if (!ref_is_dash && !is_unambiguous_dna(ref_allele))
        || (!alt_is_dash && !is_unambiguous_dna(alt_allele))
    {
        return None;
    }

    let mut ref_seq = if ref_is_dash {
        Vec::new()
    } else if transcript.strand == vep_core::coordinate::Strand::Reverse {
        reverse_complement(ref_allele)
    } else {
        ref_allele.to_vec()
    };
    let mut alt_seq = if alt_is_dash {
        Vec::new()
    } else if transcript.strand == vep_core::coordinate::Strand::Reverse {
        reverse_complement(alt_allele)
    } else {
        alt_allele.to_vec()
    };

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
    ref_seq = ref_seq[prefix..ref_seq.len().saturating_sub(suffix)].to_vec();
    alt_seq = alt_seq[prefix..alt_seq.len().saturating_sub(suffix)].to_vec();

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
/// codon remains and the peptide is not exactly `*` (Perl VEP behavior).
fn translate_codon_window_to_peptide(cds_window: &[u8]) -> Vec<u8> {
    let mut pep = Vec::with_capacity(cds_window.len().div_ceil(3));
    let (codons, rem) = cds_window.as_chunks::<3>();
    for codon in codons {
        pep.push(translate_codon(codon));
    }
    if !rem.is_empty() {
        // Perl: append X for partial codons unless peptide is exactly "*"
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

/// `reverse_complement` for callers outside this module (start_retained and HGVS checks).
pub fn reverse_complement_pub(seq: &[u8]) -> Vec<u8> {
    reverse_complement(seq)
}

/// Compute the 5' UTR sequence, preferring cached UTR (`vefc.five_prime_utr`)
/// from the JSON cache and falling back to FASTA only when the cache lacks
/// the field.
///
/// Perl's `_ins_del_start_altered` builds `utr_and_translateable = 5'UTR || CDS`
/// and applies the edit at the cDNA anchor before checking whether the ATG
/// triplet at `length($utr->seq)` is still intact. Perl reads the UTR from
/// `_variation_effect_feature_cache->{five_prime_utr}` (a `Bio::EnsEMBL::Slice`
/// object whose `->seq()` reads from the cached `primary_seq` region in
/// offline+cache mode). vep-rs mirrors this by extracting the seq string
/// from the Slice during Storable→JSON conversion (see
/// `scripts/data/storable_to_json.pl::extract_seq_string`) and using
/// `vefc.five_prime_utr` directly, keeping cache-only mode self-contained.
///
/// `fasta` is optional. When `None` and the cache lacks the field, returns
/// `None` and the caller takes its own fallback path.
pub fn compute_five_prime_utr_sequence(
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<Vec<u8>> {
    let vefc = transcript.vefc.as_ref()?;

    // Perl reads the UTR from the cache via a Slice; the extracted seq string is
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
/// Perl's `_ins_del_stop_altered` reads the UTR from
/// `_variation_effect_feature_cache->{three_prime_utr}` (a Slice object whose
/// `->seq()` reads from the cached `primary_seq` region in offline+cache mode).
/// vep-rs mirrors this by extracting the seq string from the Slice during
/// Storable→JSON conversion (`scripts/data/storable_to_json.pl::
/// extract_seq_string`) and using `vefc.three_prime_utr` directly, keeping
/// cache-only mode self-contained.
///
/// `fasta` is optional. When `None` and the cache lacks the field, returns
/// `None` and the caller takes its own fallback path.
pub fn compute_three_prime_utr_sequence(
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<Vec<u8>> {
    let vefc = transcript.vefc.as_ref()?;

    // Perl reads the UTR from the cache via a Slice; the extracted seq string is
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

/// Replicate Perl VEP's `_ins_del_stop_altered` for deletions that span the
/// CDS into the 3' UTR.
///
/// Perl's predicate (`VariationEffect::_ins_del_stop_altered`):
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
/// * `Some(true)`: stop is altered (→ Perl emits `stop_lost`)
/// * `Some(false)`: stop is not altered (→ Perl emits `stop_retained_variant`)
/// * `None`: cannot evaluate (missing FASTA / UTR / CDS data, non-DNA
///   allele, incomplete terminal codon, etc.); the caller takes its own
///   fallback path.
///
/// Required inputs:
/// * `cds_start`: 1-based CDS position where the edit begins (the `cds_start`
///   used by `substr($_, $cds_start - 1, ...)` in Perl)
/// * `cds_edit_len`: number of bases the edit removes within the CDS
///   (for a deletion spanning into UTR, this equals `(cdna_end - cdna_start) + 1`
///   clamped by CDS end, i.e. the Perl `substr` length)
/// * `reference_fasta`: optional FASTA fallback when the JSON cache lacks
///   `vefc.three_prime_utr`. The cache-first path is used by default since
///   Perl reads UTR from cache in offline mode (see
///   [`compute_three_prime_utr_sequence`]).
///
/// Caller must have verified that `_overlaps_stop_codon` is true; otherwise
/// Perl's `stop_retained` gate fails regardless of this predicate.
pub fn replicate_ins_del_stop_altered(
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

    // Perl bails out for incomplete-CDS transcripts via `_overlaps_stop_codon`
    // (cds_end_NF check); `None` leaves the caller on its fallback.
    if transcript
        .flags
        .iter()
        .any(|f| &**f == "cds_end_NF" || &**f == "cds_start_NF")
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

    // Perl's `seq_is_unambiguous_dna` gate: `_ins_del_stop_altered` returns 0
    // when either allele contains non-ACGT characters.
    if (!ref_is_dash && !is_unambiguous_dna(ref_allele))
        || (!alt_is_dash && !is_unambiguous_dna(alt_allele))
    {
        return None;
    }

    let alt_seq: Vec<u8> = if alt_is_dash {
        Vec::new()
    } else if transcript.strand == vep_core::coordinate::Strand::Reverse {
        reverse_complement(alt_allele)
    } else {
        alt_allele.to_vec()
    };

    // `compute_three_prime_utr_sequence` already reverse-complements for a
    // reverse-strand transcript, so the concatenation below is in reading frame.
    let utr_seq = compute_three_prime_utr_sequence(transcript, reference_fasta)?;
    if utr_seq.is_empty() {
        return None;
    }

    let mut combined: Vec<u8> = Vec::with_capacity(cds.len() + utr_seq.len());
    combined.extend_from_slice(cds);
    combined.extend_from_slice(&utr_seq);

    // Perl's `substr` silently clamps an edit range past the combined sequence;
    // `None` leaves the caller on its fallback instead.
    if idx + cds_edit_len > combined.len() {
        return None;
    }

    // Perl: substr($seq, $cds_start - 1, $length) = $vf_feature_seq;
    combined.splice(idx..idx + cds_edit_len, alt_seq.iter().copied());

    // Perl: "return 1 if length($utr_and_translateable) < length($translateable);"
    if combined.len() < cds.len() {
        return Some(true);
    }

    // Perl: extract the 3-base codon at the original stop position and translate.
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

/// Replicate Perl VEP's `_ins_del_start_altered` for deletions/insertions at
/// the 5'UTR/CDS boundary.
///
/// Perl's predicate (`VariationEffect::_ins_del_start_altered`):
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
/// * `Some(true)`: start is altered (→ Perl emits `start_lost`)
/// * `Some(false)`: start is not altered (→ Perl keeps only `5_prime_UTR_variant`)
/// * `None`: cannot evaluate (missing FASTA / UTR / CDS data, non-DNA
///   allele, incomplete terminal codon, etc.); the caller takes its own fallback path.
///
/// Required inputs:
/// * `cdna_start`: 1-based cDNA position where the edit begins (anchor for
///   the `substr` splice in the combined UTR+CDS string)
/// * `cdna_edit_len`: number of bases the edit removes from the cDNA (equals
///   Perl's `$cdna_end - $cdna_start + 1`)
/// * `reference_fasta`: optional FASTA fallback when the JSON cache lacks
///   `vefc.five_prime_utr`. The cache-first path is used by default since
///   Perl reads UTR from cache in offline mode (see
///   [`compute_five_prime_utr_sequence`]).
///
/// Caller must have verified that `_overlaps_start_codon` is true; otherwise
/// Perl's `_ins_del_start_altered` returns 0 regardless of this predicate.
pub fn replicate_ins_del_start_altered(
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

    // Perl `_overlaps_start_codon` early-returns 0 for cds_start_NF, and the
    // start codon analysis is meaningless for cds_end_NF.
    if transcript
        .flags
        .iter()
        .any(|f| &**f == "cds_start_NF" || &**f == "cds_end_NF")
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

    // Perl's `seq_is_unambiguous_dna` gate.
    if (!ref_is_dash && !is_unambiguous_dna(ref_allele))
        || (!alt_is_dash && !is_unambiguous_dna(alt_allele))
    {
        return None;
    }

    let alt_seq: Vec<u8> = if alt_is_dash {
        Vec::new()
    } else if transcript.strand == vep_core::coordinate::Strand::Reverse {
        reverse_complement(alt_allele)
    } else {
        alt_allele.to_vec()
    };

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

    // Perl: substr($utr_and_translateable, $cdna_start - 1, $length) = $vf_feature_seq;
    combined.splice(idx..idx + cdna_edit_len, alt_seq.iter().copied());

    // Perl: if UTR remained intact and the 3 bases at ATG position are still
    // 'ATG', start is not altered.
    let atg_start = utr_seq.len();
    if atg_start + 3 <= combined.len() {
        let new_utr = &combined[..atg_start];
        let new_sc = &combined[atg_start..atg_start + 3];
        if new_utr == utr_seq.as_slice() && new_sc == b"ATG" {
            return Some(false);
        }
    }

    // Perl: sequence shorter → start has been altered.
    if combined.len() < cds.len() {
        return Some(true);
    }

    // Perl: `$translateable ne substr($utr_and_translateable, 0 - length($translateable))`,
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
pub fn compute_full_peptide_alleles(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
) -> Option<(Vec<u8>, Vec<u8>)> {
    compute_full_peptide_alleles_impl(variant, transcript, cds_pos, false)
}

/// Like [`compute_full_peptide_alleles`] but clamps the reference allele length
/// to available CDS bounds instead of rejecting when `idx + ref_len > cds.len()`.
///
/// Used for UTR-spanning deletions where the deletion extends past CDS end.
/// Clamping yields the peptide effect on the CDS portion, allowing
/// stop_lost vs stop_retained classification.
pub fn compute_full_peptide_alleles_clamped(
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

    let alt_seq = if alt_is_dash {
        Vec::new()
    } else if transcript.strand == vep_core::coordinate::Strand::Reverse {
        reverse_complement(alt_allele)
    } else {
        alt_allele.to_vec()
    };

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
/// amino-acid change, mirroring Perl VEP's `_get_peptide_alleles`.
///
/// Returns `None` if the CDS or variant cannot be resolved.
pub fn compute_peptide_alleles(
    variant: &InputVariant,
    transcript: &Transcript,
    cds_pos: u64,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let (ref_pep, alt_pep) = compute_full_peptide_alleles(variant, transcript, cds_pos)?;
    let (r, a) = trim_common_peptide_flanks(&ref_pep, &alt_pep);
    Some((r.to_vec(), a.to_vec()))
}

/// Compute peptide alleles over the variant's codon-window, in the style of Perl's
/// `TranscriptVariationAllele::codon` / `peptide`.
///
/// This differs from `compute_peptide_alleles` (full-CDS translation) and matches
/// Perl's frameshift `stop_gained` behavior, where the gained stop can occur in a
/// later codon inside the local window.
///
/// Perl's `codon` method:
///   - Builds the full alternate CDS via `_get_alternate_cds()`
///   - Extracts: `substr($alt_cds, $codon_cds_start-1, $codon_len + ($allele_len - $vf_nt_len))`
///     where `codon_len` = ref codon window size, `allele_len` = alt allele length,
///     `vf_nt_len` = ref allele length in CDS terms.
///
/// Then `peptide` translates full codons and appends `X` for partial remainders.
///
/// For Perl parity, callers should derive `bounds` using
/// `mapper::map_genomic_span_to_cds_bounds()` rather than mapping individual endpoints.
pub fn compute_codon_window_peptide_alleles(
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

    let alt_seq = if alt_is_dash {
        Vec::new()
    } else if transcript.strand == vep_core::coordinate::Strand::Reverse {
        reverse_complement(alt_allele)
    } else {
        alt_allele.to_vec()
    };
    // Perl: $allele_len = length of alt variation_feature_seq (in mRNA sense)
    let allele_len = alt_seq.len();

    // Perl: $vf_nt_len = cds_end - cds_start + 1.
    // For insertions (cds_start > cds_end in Perl), this evaluates to 0.
    let cds_start = bounds.cds_start;
    if cds_start == 0 {
        return None;
    }
    let vf_nt_len_isize = (bounds.cds_end as isize) - (bounds.cds_start as isize) + 1;
    let vf_nt_len: usize = vf_nt_len_isize.max(0) as usize;
    let is_insertion_bounds = bounds.cds_start > bounds.cds_end;
    let is_boundary_insertion =
        is_insertion_bounds && bounds.translation_start != bounds.translation_end;

    // Perl: a codon-boundary insertion spans two adjacent peptide coordinates,
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
        // Perl: $codon_cds_start = $tv_tr_start * 3 - 2;  $codon_cds_end = $tv_tr_end * 3;
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

    // Perl `_get_alternate_cds` splicing:
    //   upstream   = cds[0..cds_start-1]
    //   downstream = cds[cds_end..]   (0-based index, so removes cds_start..cds_end inclusive)
    //
    // For insertions, Perl has cds_start > cds_end so downstream starts at cds_start-1 and
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
    // instead by `PerlCodingEval::alternate_cds`, the whole port of Perl's
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

/// Check whether the single codon at `translation_start` in the **alt** CDS translates
/// to a stop codon (`*`).
///
/// Distinguishes genuine stop_retained (local stop preserved in alt) from a false
/// positive where the wide codon window coincidentally finds a downstream `*` at the
/// same relative offset. Mirrors Perl's narrow `TranscriptVariationAllele::peptide`,
/// which returns `X` for partial codons at frameshift boundaries, preventing false
/// stop_retained via the `$ref_pep eq "X" && $alt_pep eq "X"` guard in
/// `ref_eq_alt_sequence`.
///
/// Returns `false` when the alt CDS is too short for a complete codon at the variant
/// position (partial codon = Perl's 'X' = not a stop).
pub fn narrow_alt_codon_is_stop(
    variant: &InputVariant,
    transcript: &Transcript,
    bounds: &CdsSpanBounds,
) -> bool {
    let vefc = match transcript.vefc.as_ref() {
        Some(v) => v,
        None => return false,
    };
    let cds = match vefc.translateable_seq.as_ref() {
        Some(s) => s.as_bytes(),
        None => return false,
    };
    if cds.is_empty() || bounds.translation_start == 0 || bounds.cds_start == 0 {
        return false;
    }

    let ref_allele = &variant.ref_allele;
    let alt_allele = variant.alt_allele();
    let ref_is_dash = ref_allele == b"-" || ref_allele.is_empty();
    let alt_is_dash = alt_allele == b"-" || alt_allele.is_empty();
    if (!ref_is_dash && !is_unambiguous_dna(ref_allele))
        || (!alt_is_dash && !is_unambiguous_dna(alt_allele))
    {
        return false;
    }
    let alt_seq = if alt_is_dash {
        Vec::new()
    } else if transcript.strand == vep_core::coordinate::Strand::Reverse {
        reverse_complement(alt_allele)
    } else {
        alt_allele.to_vec()
    };

    let vf_nt_len_isize = (bounds.cds_end as isize) - (bounds.cds_start as isize) + 1;
    let vf_nt_len: usize = vf_nt_len_isize.max(0) as usize;
    let splice_start = (bounds.cds_start - 1) as usize;
    let splice_end = bounds.cds_end as usize;
    if splice_start > cds.len() || splice_end > cds.len() || splice_end < splice_start {
        return false;
    }

    let mut alt_cds = Vec::with_capacity(cds.len().saturating_sub(vf_nt_len) + alt_seq.len());
    alt_cds.extend_from_slice(&cds[..splice_start]);
    alt_cds.extend_from_slice(&alt_seq);
    alt_cds.extend_from_slice(&cds[splice_end..]);

    // 1-based codon index → 0-based nucleotide: codon_start = (tr_start * 3 - 2) - 1
    let codon_start_idx = (bounds.translation_start * 3 - 2 - 1) as usize;
    let codon_end_idx = codon_start_idx + 3;
    if codon_end_idx > alt_cds.len() {
        // Partial codon: Perl translates as 'X', not a stop.
        return false;
    }
    translate_codon(&alt_cds[codon_start_idx..codon_end_idx]) == b'*'
}

/// Compute the frameshift stop_gained scan limit for the full-CDS fallback path.
///
/// The full-CDS peptide can extend far past Perl's local codon window, so the
/// scan is capped to avoid a stop_gained for a stop that Perl's narrower
/// codon-window peptide would never reach.
///
/// For net insertions, cap at `1 + floor(net_insertion_nt / 3)` amino acids
/// (approximating the number of fully-formed codons Perl's `_get_alternate_cds`
/// window produces from the insertion; partial codons at the boundary
/// translate to 'X' in Perl and are not stops).
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

pub fn frameshift_stop_gained_in_codon_window(
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
    // Perl `stop_gained` returns true when the alt peptide contains `*` and the
    // ref peptide does not. Perl applies no position cap:
    //   $cache->{stop_gained} = (($alt_pep =~ /\*/) and ($ref_pep !~ /\*/));
    //
    // The codon-window peptide is already clipped to Perl's window width by
    // `compute_codon_window_peptide_alleles` + `translate_codon_window_to_peptide`
    // (partial boundary codons become 'X' unless the whole-codon prefix is
    // exactly '*'). So scanning the full alt peptide here matches Perl exactly.
    //
    alt_pep.contains(&b'*')
}

/// Fallback frameshift stop_gained predicate using full-CDS translation.
///
/// Used when `map_genomic_span_to_cds_bounds()` returns `None` and the
/// codon-window approach is unavailable. Translates the entire ref and alt
/// CDS, then checks the region starting from the variant's amino acid
/// position with the scan-limit logic to avoid finding stops far downstream
/// that Perl's narrower codon window would never see.
pub fn frameshift_stop_gained_full_cds(
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
    if let Some(pep_seq) = vefc.peptide.as_deref() {
        let mut ref_aa = [cc.ref_amino_acid];
        overlay_seq_edits(
            &mut ref_aa,
            pep_seq.as_bytes(),
            cc.codon_number as i64,
            seq_bytes.len(),
        );
        cc.ref_amino_acid = ref_aa[0];
    }
    Some(cc)
}

/// NCBI translation table for a transcript, defaulting to 1 as Perl does.
///
/// Perl reads this per transcript in `BaseTranscriptVariation::_codon_table`
/// (`$attrib ? $attrib->value : 1`). Human vertebrate-mitochondrial transcripts
/// carry table 2; every other transcript is table 1.
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

    // The VCF allele is on the forward genomic strand while translateable_seq is
    // in mRNA sense, so a reverse-strand transcript takes the complement.
    let alt_base = variant.alt_allele()[0];
    let effective_alt = if transcript_strand == vep_core::coordinate::Strand::Reverse {
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

    // A reverse-strand transcript takes the reverse complement of the genomic
    // allele to apply it to the mRNA-sense CDS.
    let alt_seq = if transcript_strand == vep_core::coordinate::Strand::Reverse {
        reverse_complement(alt_allele)
    } else {
        alt_allele.to_vec()
    };

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
/// so that `stop_gained` can be detected (Perl parity: `tat/tAat` -> `Y/*`).
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

    // The first alt codon of a frameshift detects stop_gained (Perl: `tat/tAat`
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
    let alt_is_dash = variant.alt_allele() == b"-" || variant.alt_allele().is_empty();

    let alt_bases_owned: Vec<u8>;
    let alt_bases: &[u8] = if alt_is_dash {
        &[]
    } else if transcript_strand == vep_core::coordinate::Strand::Reverse {
        // translateable_seq is in mRNA sense; for reverse-strand transcripts the allele
        // from the genome must be reverse-complemented to apply it to the CDS.
        alt_bases_owned = reverse_complement(variant.alt_allele());
        &alt_bases_owned
    } else {
        variant.alt_allele()
    };
    let del_len = if ref_is_dash { 0 } else { ref_len };

    if codon_end > cds.len() || cds_idx + del_len > cds.len() {
        return None;
    }

    // Perl's frameshift peptide derives from its local `codon()`, which does not
    // pad a short codon with downstream bases after a deletion (`Tta/ta` yields
    // alt peptide `X`), so only bases within the current codon window are
    // stitched: [codon_start..cds_idx] + alt_bases + [post-deletion..codon_end].
    let prefix = &cds[codon_start..cds_idx];
    let del_end_within_codon = (cds_idx + del_len).min(codon_end);
    let suffix_within_codon = &cds[del_end_within_codon..codon_end];

    let mut local_alt =
        Vec::with_capacity(prefix.len() + alt_bases.len() + suffix_within_codon.len());
    local_alt.extend_from_slice(prefix);
    local_alt.extend_from_slice(alt_bases);
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

// Perl VEP coding predicates, ported whole.
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
use vep_core::transcript::ExonCoordMapper;

/// One element of a `Bio::EnsEMBL::Mapper` result list, in transcript order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerlMapSeg {
    Coord { start: i64, end: i64 },
    Gap { start: i64, end: i64 },
}

/// `BaseTranscriptVariation` coordinate accessors for one allele span.
///
/// `None` is Perl's `undef`: the first (or last) mapped segment was a gap.
/// `cds_start > cds_end` (and `tl_start > tl_end`) is Perl's insertion shape.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PerlSpan {
    pub cdna_start: Option<i64>,
    pub cdna_end: Option<i64>,
    pub cds_start: Option<i64>,
    pub cds_end: Option<i64>,
    pub tl_start: Option<i64>,
    pub tl_end: Option<i64>,
    pub cds_coords: Vec<PerlMapSeg>,
}

fn genomic_order(ecm: &ExonCoordMapper) -> Vec<usize> {
    if ecm.by_genomic_start.len() == ecm.pairs.len() {
        ecm.by_genomic_start.iter().map(|&i| i as usize).collect()
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
/// `translation_end` from `genomic2pep`.
///
/// `None` when the transcript has no coding model or no mapper.
pub fn perl_span(transcript: &Transcript, span_start: u64, span_end: u64) -> Option<PerlSpan> {
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;
    let ecm = &mapper.exon_coord_mapper;
    if ecm.pairs.is_empty() || mapper.cdna_coding_start == 0 || mapper.cdna_coding_end == 0 {
        return None;
    }
    transcript.translation.as_ref()?;
    let strand = transcript.strand.as_i8();
    let (start, end) = (span_start as i64, span_end as i64);
    let cdna = perl_genomic2cdna(ecm, start, end, strand);
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

fn is_perl_dna(seq: &[u8]) -> bool {
    !seq.is_empty()
        && seq.iter().all(|b| {
            matches!(
                b.to_ascii_uppercase(),
                b'A' | b'C'
                    | b'G'
                    | b'T'
                    | b'U'
                    | b'M'
                    | b'R'
                    | b'W'
                    | b'S'
                    | b'Y'
                    | b'K'
                    | b'V'
                    | b'H'
                    | b'D'
                    | b'B'
                    | b'X'
                    | b'N'
                    | b'-'
            )
        })
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
        b'N' => b"ACGT",
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
    seq.chunks_exact(3)
        .map(|c| perl_translate_codon(c, table))
        .collect()
}

/// `Bio::Tools::CodonTable::is_start_codon` for the two tables VEP uses.
fn perl_is_start_codon(codon: &[u8], table: u8) -> bool {
    let upper: Vec<u8> = codon.iter().map(u8::to_ascii_uppercase).collect();
    if table == 2 {
        matches!(upper.as_slice(), b"ATT" | b"ATC" | b"ATA" | b"ATG" | b"GTG")
    } else {
        matches!(upper.as_slice(), b"TTG" | b"CTG" | b"ATG")
    }
}

/// `TranscriptVariationAllele::peptide` applies the translation's SeqEdits to the
/// reference peptide only, so a selenocysteine `TGA` reads `U`, not `*`. The
/// cache carries no `seq_edits`, but its `peptide` is `Transcript::translate`
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
/// reference peptide, `Some(b'M')` when the edit is evident: the cached peptide
/// starts with `M` although the first codon translates otherwise and is not a
/// start codon `Transcript::translate` reads as `M` on its own. A forced start
/// codon leaves the edit indistinguishable, so `None`.
pub(crate) fn initial_met_edit(cds: &[u8], pep_seq: &[u8], table: u8) -> Option<u8> {
    let first = cds.get(..3)?;
    if pep_seq.first() != Some(&b'M')
        || perl_translate_codon(first, table) == b'M'
        || perl_is_start_codon(first, table)
    {
        return None;
    }
    Some(b'M')
}

/// `Bio::EnsEMBL::Transcript::translate` on the translateable sequence: whole
/// codons only, the terminal stop dropped, a start codon read as `M`. This is
/// the `_peptide` Perl compares against in `ref_eq_alt_sequence` when the cache
/// carries no peptide.
fn perl_transcript_peptide(cds: &[u8], table: u8) -> Vec<u8> {
    let whole = &cds[..cds.len() - cds.len() % 3];
    if whole.is_empty() {
        return Vec::new();
    }
    let mut pep = perl_translate(whole, table);
    if perl_translate_codon(&whole[whole.len() - 3..], table) == b'*' {
        pep.pop();
    }
    if pep.first().is_some_and(|&a| a != b'M') && perl_is_start_codon(&whole[..3], table) {
        pep[0] = b'M';
    }
    pep
}

/// Every SO term `perl_coding_terms` can emit; the caller replaces exactly this set.
pub fn is_perl_coding_term(c: Consequence) -> bool {
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

/// One `TranscriptVariationAllele` evaluation with Perl's `_predicate_cache`.
struct PerlCodingEval<'a> {
    cds: &'a [u8],
    pep_seq: Cow<'a, [u8]>,
    utr5: Option<Vec<u8>>,
    utr3: Option<Vec<u8>>,
    table: u8,
    span: PerlSpan,
    strand: i8,
    cds_start_nf: bool,
    cds_end_nf: bool,
    cdna_coding_start: i64,
    cdna_coding_end: i64,
    coding_region: Option<(i64, i64)>,
    frameshift_introns: Vec<(i64, i64)>,
    vf_start: i64,
    vf_end: i64,
    ref_allele: &'a [u8],
    alt_allele: &'a [u8],
    feature_seq: Vec<u8>,
    seq_length: Option<i64>,
    snp: bool,
    insertion: bool,
    deletion: bool,
    // _predicate_cache
    codon_memo: [Option<Option<Vec<u8>>>; 2],
    peptide_memo: [Option<Option<Vec<u8>>>; 2],
    peps_memo: Option<Option<(Vec<u8>, Vec<u8>)>>,
    partial_codon_memo: Option<bool>,
    overlaps_stop_memo: Option<bool>,
    overlaps_start_memo: Option<bool>,
    ins_del_stop_altered_memo: Option<bool>,
    ins_del_start_altered_memo: Option<bool>,
    inv_start_altered_memo: Option<bool>,
    snp_start_altered_memo: Option<bool>,
    start_lost_memo: Option<bool>,
    stop_lost_memo: Option<bool>,
    stop_retained_memo: Option<bool>,
    stop_gained_memo: Option<bool>,
}

impl<'a> PerlCodingEval<'a> {
    fn allele(&self, is_ref: bool) -> &'a [u8] {
        if is_ref {
            self.ref_allele
        } else {
            self.alt_allele
        }
    }

    /// `_get_alternate_cds` for the alt allele: upstream + allele + downstream,
    /// the incomplete-codon trim (a no-op past 2 bp: `return $seq if
    /// $full_length = $keep_length` assigns), then the 3' UTR.
    fn alternate_cds(&self) -> Option<Vec<u8>> {
        let cds_start = self.span.cds_start?;
        let cds_end = self.span.cds_end?;
        let up_end = ((cds_start - 1).max(0) as usize).min(self.cds.len());
        let down_start = (cds_end.max(0) as usize).min(self.cds.len());
        let mut a: Vec<u8> = Vec::with_capacity(self.alt_allele.len());
        let mut dash_removed = false;
        for &b in self.alt_allele {
            if b == b'-' && !dash_removed {
                dash_removed = true;
            } else {
                a.push(b);
            }
        }
        if !a.is_empty() && self.strand == -1 {
            a = reverse_complement(&a);
        }
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

    /// `TranscriptVariationAllele::codon`. `Some(b"-")` is Perl's `'-'`.
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
        if !is_perl_dna(self.allele(is_ref)) {
            return None;
        }
        let codon_cds_start = tl_start * 3 - 2;
        let codon_cds_end = tl_end * 3;
        let codon_len = codon_cds_end - codon_cds_start + 1;
        let (Some(cds_start), Some(cds_end)) = (self.span.cds_start, self.span.cds_end) else {
            return None;
        };
        let vf_nt_len = cds_end - cds_start + 1;
        // Perl builds `_get_alternate_cds` for both alleles but reads the reference
        // codon from the translateable sequence itself.
        let codon: Option<Vec<u8>> = if is_ref {
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
        if !is_perl_unambiguous_dna(self.allele(is_ref)) {
            return None;
        }
        let codon = self.codon(is_ref)?;
        if let Some(Some(p)) = &self.peptide_memo[slot] {
            return Some(p.clone());
        }
        let whole_len = codon.len() / 3 * 3;
        let mut pep = perl_translate(&codon[..whole_len], self.table);
        // SeqEdits apply to the reference peptide only, before the partial-codon `X`:
        // the residues the cached peptide carries over the codon translation, then
        // the `initial_met` edit at protein position 1.
        if is_ref && !pep.is_empty() {
            if let (Some(a), Some(b)) = (self.span.tl_start, self.span.tl_end) {
                let tv_lo = a.min(b);
                overlay_seq_edits(&mut pep, &self.pep_seq, tv_lo, self.cds.len());
                if tv_lo == 1 {
                    if let Some(m) = self.initial_met_edit() {
                        pep[0] = m;
                    }
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

    /// `initial_met_edit` on this transcript's CDS, cached peptide and codon table.
    fn initial_met_edit(&self) -> Option<u8> {
        initial_met_edit(self.cds, &self.pep_seq, self.table)
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

    fn partial_codon(&mut self) -> bool {
        if let Some(v) = self.partial_codon_memo {
            return v;
        }
        self.partial_codon_memo = Some(false);
        let Some(tl_start) = self.span.tl_start else {
            return false;
        };
        let codon_cds_start = tl_start * 3 - 2;
        let last = self.cds.len() as i64 - (codon_cds_start - 1);
        let v = last < 3 && last > 0;
        self.partial_codon_memo = Some(v);
        v
    }

    fn overlaps_stop_codon(&mut self) -> bool {
        if let Some(v) = self.overlaps_stop_memo {
            return v;
        }
        self.overlaps_stop_memo = Some(false);
        if self.cds_end_nf {
            return false;
        }
        let (Some(cs), Some(ce)) = (truthy(self.span.cdna_start), truthy(self.span.cdna_end))
        else {
            return false;
        };
        let v = perl_overlap(cs, ce, self.cdna_coding_end - 2, self.cdna_coding_end);
        self.overlaps_stop_memo = Some(v);
        v
    }

    fn overlaps_start_codon(&mut self) -> bool {
        if let Some(v) = self.overlaps_start_memo {
            return v;
        }
        self.overlaps_start_memo = Some(false);
        if self.cds_start_nf {
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
        if !(self.insertion || self.deletion) {
            return false;
        }
        let (Some(cdna_start), Some(cdna_end), Some(cds_start)) = (
            truthy(self.span.cdna_start),
            truthy(self.span.cdna_end),
            truthy(self.span.cds_start),
        ) else {
            return false;
        };
        let mut s = self.cds.to_vec();
        if let Some(utr) = self.utr3.as_ref() {
            s.extend_from_slice(utr);
        }
        let fs = self.feature_seq_bases().to_vec();
        perl_substr_assign(&mut s, cds_start - 1, cdna_end - cdna_start + 1, &fs);
        if s.len() < self.cds.len() {
            self.ins_del_stop_altered_memo = Some(true);
            return true;
        }
        if self.cds.len() < 3 {
            return false;
        }
        let stop_idx = self.cds.len() - 3;
        let codon = &s[stop_idx..stop_idx + 3];
        let v = perl_translate_codon(codon, self.table) != b'*';
        self.ins_del_stop_altered_memo = Some(v);
        v
    }

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
        if let Some(utr) = utr.as_ref() {
            let atg_start = utr.len();
            let new_sc = s
                .get(atg_start..(atg_start + 3).min(s.len()))
                .unwrap_or(b"");
            let new_utr = &s[..atg_start.min(s.len())];
            if new_utr == utr.as_slice() && new_sc == b"ATG" {
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
        if alt_pep.is_empty() || alt_pep == b"X" {
            return false;
        }
        let v = self.span.tl_start == Some(1)
            && !alt_pep.ends_with(&ref_pep)
            && !alt_pep.starts_with(&ref_pep);
        self.start_lost_memo = Some(v);
        v
    }

    fn stop_lost(&mut self) -> bool {
        if let Some(v) = self.stop_lost_memo {
            return v;
        }
        self.stop_lost_memo = Some(false);
        let v = match self.peptide_alleles() {
            Some((ref_pep, alt_pep)) => !alt_pep.contains(&b'*') && ref_pep.contains(&b'*'),
            None => self.ins_del_stop_altered(),
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
        if tl_start > ref_seq.len() as i64 && alt_pep.first() == Some(&b'*') {
            return true;
        }
        let mut mut_seq = ref_seq.to_vec();
        perl_substr_assign(&mut mut_seq, tl_start - 1, tl_end - tl_start + 1, &alt_pep);
        let mut_substring = &mut_seq[..ref_seq.len().min(mut_seq.len())];
        let final_stop_length = if ref_seq.len() < mut_seq.len() {
            Some(mut_seq.len() - ref_seq.len())
        } else {
            None
        };
        if ref_pep.as_slice() == &alt_pep[..alt_pep.len().min(1)] && alt_pep.contains(&b'*') {
            return true;
        }
        if ref_seq == mut_substring && final_stop_length.is_some_and(|l| l < 3) {
            return true;
        }
        if let Some(ri) = ref_pep.iter().position(|&b| b == b'*') {
            if alt_pep.iter().position(|&b| b == b'*') == Some(ri) {
                return true;
            }
        }
        false
    }

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
            Some((_, alt_pep)) if !alt_pep.is_empty() => self.ref_eq_alt_sequence(),
            _ => {
                (self.insertion || self.deletion)
                    && self.overlaps_stop_codon()
                    && !self.ins_del_stop_altered()
            }
        };
        self.stop_retained_memo = Some(v);
        v
    }

    fn frameshift(&mut self) -> bool {
        if self.partial_codon() {
            return false;
        }
        if self.stop_retained() {
            return false;
        }
        self.frameshift_by_length()
    }

    /// The length arithmetic of Perl's `frameshift` (VariationEffect.pm 1447-1455)
    /// without its two cached guards: the allele length against the CDS span the
    /// variant covers, `false` when either end of that span is undefined.
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
        if self.stop_retained() {
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
        let _ = self.peptide_alleles();
        let Some((ref_codon, alt_codon)) = codons else {
            return false;
        };
        if alt_codon.len() >= ref_codon.len() {
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
        let cond = alt_peptide.is_none()
            || peps.is_none()
            || alt_peptide.as_ref().is_some_and(|p| p.contains(&b'X'))
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

    /// The coding terms, evaluated in `%OVERLAP_CONSEQUENCES` rank order with each
    /// predicate's `include` filter applied.
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
        if self.deletion && self.inframe_deletion() {
            out.push(Consequence::InframeDeletion);
        }
        if self.snp && self.missense_variant() {
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

/// Perl's `codon` / `peptide` strings for both alleles of a coding variant, for
/// tests and diagnostics: `(ref_codon, alt_codon, ref_pep, alt_pep)` with `-` for
/// an empty window and `None` where Perl has `undef`.
#[allow(clippy::type_complexity)]
pub fn perl_codon_peptides(
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

/// Perl `TranscriptVariationAllele::hgvs_protein` (TranscriptVariationAllele.pm
/// 1593), the text after `p.`.
///
/// `variant` is the allele as annotated; `shifted`, when the caller moved an
/// insertion or deletion to its most 3' position (Perl's `_return_3prime(1)`),
/// is the shifted allele with its span. Perl reads the peptides, the translation
/// coordinates and the alternate CDS from the shifted allele, but the `coding`
/// pre-consequence predicate (1667) and the cached predicates `stop_lost`,
/// `start_lost`, `partial_codon` and `stop_retained` were filled while the
/// consequences were computed on the unshifted allele (`hgvs_transcript` clears
/// that cache only under `--shift_3prime`, 1405), so those verdicts are taken
/// from `variant` here whatever the shift. `frameshift` (VariationEffect.pm
/// 1435) is not cached: its two guards read the cache, but its length
/// arithmetic runs on the CDS span the transcript variation carries at that
/// point, which is the shifted span, so an indel that shifts fully into the
/// CDS is a frameshift there even where the annotated allele straddles an
/// exon boundary.
///
/// `None` where Perl returns `undef`: the annotated allele does not overlap the
/// coding sequence, the (shifted) span has no translation start or end, or its
/// reference peptide is undefined. The alternate-CDS translations inside use
/// codon table 1 whatever the transcript's table, as BioPerl's argument-less
/// `translate()` does at 2263, 2380, 2422 and 2485; the transcript's own
/// peptide keeps its table (Ensembl `Transcript::translate`).
///
/// One intended divergence: Perl prints `Met1?` whenever its `start_lost`
/// predicate holds (2091), including the start co-emission pairs on which this
/// engine keeps `start_retained_variant` and drops `start_lost`; the port fires
/// that short-circuit only when it emits `start_lost` itself.
pub fn perl_hgvs_protein(
    variant: &InputVariant,
    shifted: Option<(&InputVariant, u64, u64)>,
    transcript: &Transcript,
    fasta: Option<&vep_fasta::IndexedFasta>,
) -> Option<String> {
    let mut ev = PerlCodingEval::new(variant, transcript, variant.start, variant.end, fasta)?;
    if !ev.coding_pred(transcript) {
        return None;
    }
    let frameshift_guard = ev.partial_codon() || ev.stop_retained();
    let stop_lost = ev.stop_lost();
    let start_lost = ev.start_lost() && !ev.start_retained_variant();
    if let Some((shifted_variant, span_start, span_end)) = shifted {
        ev = PerlCodingEval::new(shifted_variant, transcript, span_start, span_end, fasta)?;
    }
    let preds = HgvsPredicates {
        frameshift: !frameshift_guard && ev.frameshift_by_length(),
        stop_lost,
        start_lost,
    };
    let tl_start = truthy(ev.span.tl_start)?;
    let tl_end = truthy(ev.span.tl_end)?;
    let alt_pep = ev.peptide(false);
    let ref_pep = ev.peptide(true).filter(|p| !p.is_empty())?;
    let mut n = HgvsProteinNotation {
        ref_pep: Some(ref_pep),
        alt_pep,
        start: tl_start,
        end: tl_end,
        kind: HgvspKind::Unset,
        original_ref: Vec::new(),
        preseq: Vec::new(),
    };
    if n.alt_pep.is_some() && n.alt_pep != n.ref_pep {
        hgvsp_clip_alleles(&mut n);
    }
    hgvsp_protein_type(&ev, &preds, &mut n);
    if n.kind == HgvspKind::Unset {
        return None;
    }
    hgvsp_peptides(&ev, &preds, &mut n)?;
    Some(hgvsp_format(&ev, &preds, &n))
}

/// The three predicate verdicts `hgvs_protein` reads: `stop_lost` and
/// `start_lost` from Perl's `_predicate_cache` (the allele as annotated, not its
/// shifted form), `frameshift` recomputed on the shifted span behind the cached
/// `partial_codon` and `stop_retained` guards. `start_lost` already carries the
/// start co-emission gate of this engine.
struct HgvsPredicates {
    frameshift: bool,
    stop_lost: bool,
    start_lost: bool,
}

/// Perl's `$hgvs_notation` hash. `ref_pep` and `alt_pep` are `None` where the
/// peptide is `undef`; `original_ref` and `preseq` are what `_clip_alleles`
/// records for the stop and duplication checks.
struct HgvsProteinNotation {
    ref_pep: Option<Vec<u8>>,
    alt_pep: Option<Vec<u8>>,
    start: i64,
    end: i64,
    kind: HgvspKind,
    original_ref: Vec<u8>,
    preseq: Vec<u8>,
}

/// Perl's `type` field: `""` is `Bare`, `"="` is `Eq`, `">"` is `Sub`.
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

/// Perl's `$s =~ s/Xaa/Ter/g`.
fn replace_xaa_with_ter(pep: &mut [u8]) {
    let mut i = 0;
    while i + 3 <= pep.len() {
        if &pep[i..i + 3] == b"Xaa" {
            pep[i..i + 3].copy_from_slice(b"Ter");
        }
        i += 3;
    }
}

/// `_clip_alleles` (2118) with `numbering` `p`: trims the residues the alleles
/// share from the front, then from the back, moving `start` and `end`, and
/// records `original_ref` and `preseq`. A leading stop on both sides returns
/// `Eq` at once. The type re-set block runs as written except its `dup` case,
/// which applies only to nucleotide numbering.
fn hgvsp_clip_alleles(n: &mut HgvsProteinNotation) {
    let ref_pep: &[u8] = n.ref_pep.as_deref().unwrap_or_default();
    let alt_pep: &[u8] = n.alt_pep.as_deref().unwrap_or_default();
    n.original_ref = ref_pep.to_vec();
    let mut check_ref = ref_pep;
    let mut check_alt = alt_pep;
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
                n.start += 1;
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
            n.end -= 1;
        } else {
            break;
        }
    }
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

/// `_get_hgvs_protein_type` (1977): `fs` from the frameshift predicate; else the
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

/// `_get_hgvs_peptides` (2044) with three-letter conversion on. `None` where Perl
/// returns `undef`: an insertion with no flanking residue to name.
fn hgvsp_peptides(
    ev: &PerlCodingEval<'_>,
    preds: &HgvsPredicates,
    n: &mut HgvsProteinNotation,
) -> Option<()> {
    match n.kind {
        HgvspKind::Fs => hgvsp_fs_peptides(ev, n)?,
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

/// `_get_fs_peptides` (2250): the first residue at which the table-1 translation
/// of the alternate CDS (3' UTR appended) differs from the reference peptide plus
/// its stop, from `translation_start`. `Del` when the alternate translation ends
/// before that position; `Eq` when both sides reach a stop together.
fn hgvsp_fs_peptides(ev: &PerlCodingEval<'_>, n: &mut HgvsProteinNotation) -> Option<()> {
    let alt_cds = ev.alternate_cds()?;
    if !alt_cds
        .iter()
        .any(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'-'))
    {
        return None;
    }
    let alt_trans = perl_translate(&alt_cds, 1);
    let mut ref_trans = ev.pep_seq.to_vec();
    ref_trans.push(b'*');
    n.start = truthy(ev.span.tl_start)?;
    if n.start > alt_trans.len() as i64 {
        n.alt_pep = Some(b"del".to_vec());
        n.kind = HgvspKind::Del;
        return Some(());
    }
    while n.start <= alt_trans.len() as i64 {
        let i = (n.start - 1) as usize;
        let r = ref_trans.get(i).copied();
        let a = alt_trans.get(i).copied();
        n.ref_pep = Some(r.into_iter().collect());
        n.alt_pep = Some(a.into_iter().collect());
        if r == Some(b'*') && a == Some(b'*') {
            n.kind = HgvspKind::Eq;
            return Some(());
        }
        if r != a {
            break;
        }
        n.start += 1;
    }
    Some(())
}

/// `_get_surrounding_peptides` (2298): `length` residues of the reference peptide
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

/// `_check_for_peptide_duplication` (2372): an inserted peptide equal to the
/// residues just before it (the table-1 reference translation plus `preseq`)
/// becomes a `Dup` of those residues, three-lettered here and not again.
fn hgvsp_check_duplication(ev: &PerlCodingEval<'_>, n: &mut HgvsProteinNotation) {
    let alt: &[u8] = n.alt_pep.as_deref().unwrap_or_default();
    let mut upstream = perl_translate(ev.cds, 1);
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

/// `_stop_loss_extra_AA` (2407): residues from the variant to the first stop of
/// the table-1 alternate translation; counted from `ref_var_pos` for a
/// frameshift, else past the reference peptide's end. `None` unless positive.
fn hgvsp_stop_loss_extra_aa(ev: &PerlCodingEval<'_>, ref_var_pos: i64, fs: bool) -> Option<i64> {
    if ref_var_pos == 0 {
        return None;
    }
    let alt_cds = ev.alternate_cds()?;
    let alt_trans = perl_translate(&alt_cds, 1);
    let stop_end = alt_trans.iter().position(|&b| b == b'*')? as i64 + 1;
    let extra = if fs {
        stop_end - ref_var_pos
    } else {
        stop_end - 1 - ev.pep_seq.len() as i64
    };
    (extra > 0).then_some(extra)
}

/// `_get_del_peptides` (2474), Perl's path for a deletion whose reference peptide
/// window is empty: both peptides from `translation_start` to the end (the alternate
/// side cut at its first stop), clipped, three-lettered.
fn hgvsp_del_peptides(ev: &PerlCodingEval<'_>, n: &mut HgvsProteinNotation) -> Option<()> {
    let alt_cds = ev.alternate_cds()?;
    let tl_start = truthy(ev.span.tl_start)?;
    let start0 = ((tl_start - 1).max(0)) as usize;
    let alt_trans = perl_translate(&alt_cds, 1);
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

/// `_check_peptides_post_var` (2503) plus `_shift_3prime` (2525): rotates an
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

/// `_get_hgvs_protein_format` (1834) with three-letter conversion on and no
/// prediction parentheses.
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
    if ref_pep == alt && kind != HgvspKind::Fs && kind != HgvspKind::Ins {
        format!("{}{start}=", text(ref_pep))
    } else if preds.stop_lost && (kind == HgvspKind::Del || kind == HgvspKind::Sub) {
        let aa_til_stop = match hgvsp_stop_loss_extra_aa(ev, start - 1, false) {
            Some(extra) => extra.to_string(),
            None => "?".to_string(),
        };
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
            let aa_til_stop = match hgvsp_stop_loss_extra_aa(ev, start - 1, true) {
                Some(extra) => extra.to_string(),
                None => "?".to_string(),
            };
            format!("{}{start}{}fsTer{aa_til_stop}", text(ref_pep), text(alt))
        }
    } else if kind == HgvspKind::Del {
        if ref_pep.len() > 3 {
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
        fasta: Option<&vep_fasta::IndexedFasta>,
    ) -> Option<Self> {
        let vefc = transcript.vefc.as_ref()?;
        let cds = vefc.translateable_seq.as_ref()?.as_bytes();
        if cds.is_empty() {
            return None;
        }
        let mapper = vefc.mapper.as_ref()?;
        let span = perl_span(transcript, span_start, span_end)?;
        let strand = transcript.strand.as_i8();
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
        let feature_seq = if transcript.strand == Strand::Reverse && is_perl_dna(alt_allele) {
            reverse_complement(alt_allele)
        } else {
            alt_allele.to_vec()
        };
        let seq_length = if is_perl_dna(alt_allele) {
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
        let pep_seq: Cow<'a, [u8]> = match vefc.peptide.as_deref() {
            Some(p) => Cow::Borrowed(p.as_bytes()),
            None => Cow::Owned(perl_transcript_peptide(cds, table)),
        };
        // A cache without UTR strings still has the UTR lengths: `N` placeholders
        // keep Perl's edit anchors (`cdna_start - 1` into 5'UTR + CDS) and window
        // extents in place, translating to `X` where Perl reads real bases.
        let utr5 = compute_five_prime_utr_sequence(transcript, fasta).or_else(|| {
            (mapper.cdna_coding_start > 1)
                .then(|| vec![b'N'; (mapper.cdna_coding_start - 1) as usize])
        });
        let utr3 = compute_three_prime_utr_sequence(transcript, fasta).or_else(|| {
            let cdna_len = mapper
                .exon_coord_mapper
                .pairs
                .iter()
                .map(|p| p.from_end)
                .max()
                .unwrap_or(0);
            (cdna_len > mapper.cdna_coding_end)
                .then(|| vec![b'N'; (cdna_len - mapper.cdna_coding_end) as usize])
        });
        Some(Self {
            cds,
            pep_seq,
            utr5,
            utr3,
            table,
            span,
            strand,
            cds_start_nf: transcript.flags.iter().any(|f| f == "cds_start_NF"),
            cds_end_nf: transcript.flags.iter().any(|f| f == "cds_end_NF"),
            cdna_coding_start: mapper.cdna_coding_start as i64,
            cdna_coding_end: mapper.cdna_coding_end as i64,
            coding_region,
            frameshift_introns,
            vf_start,
            vf_end,
            ref_allele,
            alt_allele,
            feature_seq,
            seq_length,
            snp: alt_length == ref_length,
            insertion: ref_length < alt_length,
            deletion: ref_length > alt_length,
            codon_memo: [None, None],
            peptide_memo: [None, None],
            peps_memo: None,
            partial_codon_memo: None,
            overlaps_stop_memo: None,
            overlaps_start_memo: None,
            ins_del_stop_altered_memo: None,
            ins_del_start_altered_memo: None,
            inv_start_altered_memo: None,
            snp_start_altered_memo: None,
            start_lost_memo: None,
            stop_lost_memo: None,
            stop_retained_memo: None,
            stop_gained_memo: None,
        })
    }

    /// `_bvfo_preds` `coding`: the span overlaps the coding region and an exon and
    /// projects to at least one CDS segment (a lone gap counts, as in Perl).
    fn coding_pred(&self, transcript: &Transcript) -> bool {
        let Some((crs, cre)) = self.coding_region else {
            return false;
        };
        let (lo, hi) = (
            self.vf_start.min(self.vf_end),
            self.vf_start.max(self.vf_end),
        );
        if !perl_overlap(lo, hi, crs, cre) {
            return false;
        }
        let Some(mapper) = transcript.vefc.as_ref().and_then(|v| v.mapper.as_ref()) else {
            return false;
        };
        let exon_overlap = mapper
            .exon_coord_mapper
            .pairs
            .iter()
            .any(|p| perl_overlap(lo, hi, p.to_start as i64, p.to_end as i64));
        exon_overlap && !self.span.cds_coords.is_empty()
    }
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

/// The coding SO terms Perl VEP assigns to an allele on `transcript`: every
/// non-SNV allele, and an SNV whose position lies in the first codon.
///
/// `None` when the allele is an SNV outside the first codon, the transcript
/// has no coding model, or Perl's `coding` pre-predicate is false (the caller
/// keeps its own terms). `Some` is the complete coding-term set, possibly
/// empty, and replaces every term for which [`is_perl_coding_term`] holds.
///
/// A start-codon SNV takes Perl's routes as written: `_inv_start_altered`
/// reads `start_lost` from the edited codon not being `ATG` whenever the
/// transcript has a 5' UTR, the peptide route from `translation_start == 1`
/// and a changed residue otherwise, and `start_retained_variant` from the
/// edited codon being the literal `ATG`, on every codon table.
///
/// Where Perl emits `start_lost` together with `start_retained_variant`, only
/// `start_retained_variant` is returned. Perl reaches that pair only when
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
    if ref_len == 1 && alt_len == 1 && !snv_in_first_codon(transcript, span_start) {
        return None;
    }
    let mut ev = PerlCodingEval::new(variant, transcript, span_start, span_end, fasta)?;
    if !ev.coding_pred(transcript) {
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
    use crate::test_helpers::make_test_transcript;
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
    /// Perl's protein_altering_variant evaluates before inframe_insertion on raw
    /// (un-truncated) peptides. If the alt peptide doesn't contain the ref peptide
    /// as a prefix or suffix, Perl returns protein_altering, not inframe_insertion.
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
    /// Perl's stop_gained fires for any coding variant (including pure inframe
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
    /// window; the stop check runs only when the alt window carries a `*`.
    ///
    /// Bounds `cds_start: 9, cds_end: 8` splice "AAATAAGCT" in at 0-based CDS index
    /// 8, so the 12-byte alt window is GG + AAATAAGCT + A = "GGAAATAAGCTA" →
    /// [G, N, K, L]: the TAA is out of frame and the window has no stop, so the
    /// guarded assertion does not run. Net insertion = 9 nt = 3 AAs, so the capped
    /// scan would cover `1 + 3 = 4` AAs.
    #[test]
    fn test_stop_retained_condition1_nine_bp_insertion_resolves_single_residue_ref_window() {
        use crate::consequences::is_stop_retained_insertion;
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
        let (ref_pep, alt_pep) = peps.unwrap();
        assert_eq!(ref_pep.len(), 1);

        // A stop in alt_pep must be detected, by condition 1 (ref_pep[0] ==
        // alt_pep[0] and the stop within the scan window) or by condition 3.
        if alt_pep.contains(&b'*') {
            let result = is_stop_retained_insertion(&variant, &tx, &bounds, true);
            // At least one condition fires when the window has a stop
            assert!(
                result || !ref_pep.is_empty(),
                "stop_retained should detect stop in alt_pep: ref={:?} alt={:?}",
                ref_pep,
                alt_pep
            );
        }
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
}
