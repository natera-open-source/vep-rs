// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Structural variant consequence calculation.
//!
//! This module handles consequence assignment for structural variants (symbolic
//! alleles like `<DEL>`, `<DUP>`, `<INV>`, BND notation, etc.). The approach mirrors
//! Perl VEP's `StructuralVariationOverlapAllele` logic.
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`).

pub mod breakend;
pub mod cnv;
pub mod deletion;
pub mod duplication;
pub mod insertion;
pub mod inversion;

use smallvec::SmallVec;
use vep_core::consequence::{Consequence, ConsequenceList, Impact, TranscriptConsequence};
use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;
use vep_core::variant::{InputVariant, VariantClass};

/// Perl's default `--max_sv_size` threshold (10M). A single-breakend BND whose
/// span exceeds it is routed to `calculate_derived_single_breakend`, which emits
/// `feature_truncation`, a biotype context term and `intron_variant` from two interval tests.
const PERL_DEFAULT_MAX_SV_SIZE: u64 = 10_000_000;

/// Check whether an SV overlaps the mature miRNA subfeature of a miRNA transcript.
///
/// Perl VEP emits `mature_miRNA_variant` when the SV overlaps the genomic region
/// corresponding to the cDNA range in the transcript's `miRNA` attribute. This
/// mirrors `within_mature_miRNA` in `StructuralVariationOverlapAllele.pm`.
pub(crate) fn is_mature_mirna_sv(transcript: &Transcript, sv_start: u64, sv_end: u64) -> bool {
    if &*transcript.biotype != "miRNA" {
        return false;
    }

    let vefc = match transcript.vefc.as_ref() {
        Some(v) => v,
        None => return false,
    };
    let mapper = match vefc.mapper.as_ref() {
        Some(m) => m,
        None => return false,
    };

    for attr in &transcript.attributes {
        if attr.code != "miRNA" {
            continue;
        }
        if let Some((cdna_start, cdna_end)) = parse_mirna_cdna_range(&attr.value) {
            for pair in &mapper.exon_coord_mapper.pairs {
                if cdna_start > pair.from_end || cdna_end < pair.from_start {
                    continue;
                }
                let clamped_cdna_lo = cdna_start.max(pair.from_start);
                let clamped_cdna_hi = cdna_end.min(pair.from_end);

                // Map cDNA to genomic. For forward strand (ori >= 0), from_start
                // maps to to_start. For reverse strand (ori < 0), from_start maps
                // to to_end.
                let (genomic_lo, genomic_hi) = if pair.ori >= 0 {
                    let lo = pair.to_start + (clamped_cdna_lo - pair.from_start);
                    let hi = pair.to_start + (clamped_cdna_hi - pair.from_start);
                    (lo, hi)
                } else {
                    let lo = pair.to_end - (clamped_cdna_hi - pair.from_start);
                    let hi = pair.to_end - (clamped_cdna_lo - pair.from_start);
                    (lo, hi)
                };

                if overlap_bp(sv_start, sv_end, genomic_lo, genomic_hi) > 0 {
                    return true;
                }
            }
        }
    }
    false
}

/// Parse a "N-M" pattern from a miRNA attribute value, returning (start, end) as cDNA positions.
fn parse_mirna_cdna_range(value: &str) -> Option<(u64, u64)> {
    let bytes = value.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        while i < bytes.len() && !bytes[i].is_ascii_digit() {
            i += 1;
        }
        if i >= bytes.len() {
            break;
        }
        let start_i = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        let left: u64 = std::str::from_utf8(&bytes[start_i..i]).ok()?.parse().ok()?;

        if i >= bytes.len() || bytes[i] != b'-' {
            continue;
        }
        i += 1;
        if i >= bytes.len() || !bytes[i].is_ascii_digit() {
            continue;
        }
        let end_i = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        let right: u64 = std::str::from_utf8(&bytes[end_i..i]).ok()?.parse().ok()?;
        return Some((left.min(right), left.max(right)));
    }
    None
}

/// Compute the number of overlapping base pairs between two intervals.
pub(crate) fn overlap_bp(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> u64 {
    let lo = a_start.max(b_start);
    let hi = a_end.min(b_end);
    if lo <= hi {
        hi - lo + 1
    } else {
        0
    }
}

pub(crate) fn overlaps_point(start: u64, end: u64, point: u64) -> bool {
    let lo = start.min(end);
    let hi = start.max(end);
    lo <= point && hi >= point
}

pub(crate) fn overlaps_any_exon(transcript: &Transcript, sv_start: u64, sv_end: u64) -> bool {
    transcript
        .exons
        .iter()
        .any(|exon| overlap_bp(sv_start, sv_end, exon.start, exon.end) > 0)
}

fn transcript_introns(transcript: &Transcript) -> &[vep_core::transcript::Intron] {
    if let Some(ref vefc) = transcript.vefc {
        if !vefc.introns.is_empty() {
            return &vefc.introns;
        }
    }
    &transcript.introns
}

/// Check whether an SV overlaps any intron using Perl VEP's *trimmed* boundaries.
///
/// Perl's `_intron_effects` (BaseTranscriptVariationAllele.pm:143-150) sets the
/// `intronic` flag with `overlap($r_start, $r_end, $intron_start + 2, $intron_end - 2)`
/// since the 2 bases at each intron end are the invariant splice donor (GT) and
/// acceptor (AG) dinucleotides, which Perl attributes to the splice-site terms
/// rather than to `intron_variant`. `within_intron` (VariationEffect.pm:629)
/// returns exactly that flag, so an SV that touches an intron only within those
/// 4 boundary bases gets no `intron_variant` from Perl.
///
/// An intron of <= 4 bp has no interior at all under this rule and can never be
/// `intronic`.
pub(crate) fn overlaps_any_intron_trimmed(
    transcript: &Transcript,
    sv_start: u64,
    sv_end: u64,
) -> bool {
    transcript_introns(transcript).iter().any(|intron| {
        let trimmed_start = intron.start.saturating_add(2);
        let trimmed_end = intron.end.saturating_sub(2);
        if trimmed_start > trimmed_end {
            return false;
        }
        overlap_bp(sv_start, sv_end, trimmed_start, trimmed_end) > 0
    })
}

/// Check if an SV span overlaps the polypyrimidine tract of any intron.
///
/// PPT window: positions -16 to -2 from the acceptor splice site (a 15bp window).
/// On forward strand the acceptor is at the intron end, on reverse at the intron start.
/// Coordinates use the same 1-based inclusive convention as intron boundaries.
///
/// Not called from the SV consequence paths: Perl adds no PPT term to structural
/// deletions, and a large span always overlaps the window.
#[allow(dead_code)]
pub(crate) fn overlaps_polypyrimidine_tract(
    transcript: &Transcript,
    sv_start: u64,
    sv_end: u64,
) -> bool {
    let lo = sv_start.min(sv_end);
    let hi = sv_start.max(sv_end);
    let introns = transcript_introns(transcript);

    introns.iter().any(|intron| {
        // PPT is near the acceptor end of the intron.
        // Forward strand: acceptor at intron_end, PPT = [intron_end - 16, intron_end - 2]
        // Reverse strand: acceptor at intron_start, PPT = [intron_start + 2, intron_start + 16]
        let (ppt_start, ppt_end) = match transcript.strand {
            Strand::Forward => (intron.end.saturating_sub(16), intron.end.saturating_sub(2)),
            Strand::Reverse => (
                intron.start.saturating_add(2),
                intron.start.saturating_add(16),
            ),
        };
        overlap_bp(lo, hi, ppt_start, ppt_end) > 0
    })
}

/// Whether one of a BND's breakends lies close enough to a transcript for Ensembl VEP to
/// annotate that transcript at all.
///
/// This is feature selection, and it is a separate question from which coordinates the
/// consequence terms are then computed over. Ensembl keeps the two apart:
/// `StructuralVariationOverlap::_close_to_feature` tests each entry of
/// `($vf, @{$vf->get_breakends})` against the transcript slice expanded by
/// `MAX_DISTANCE_FROM_TRANSCRIPT` (5000, `Utils/VariationEffect.pm:60`, the same value as
/// the default up/downstream distance), while every positional predicate afterwards reads
/// the local variation feature's span.
fn bnd_breakend_close_to_transcript(
    variant: &InputVariant,
    transcript: &Transcript,
    max_distance: u64,
) -> bool {
    let lo = transcript.start.saturating_sub(max_distance);
    let hi = transcript.end.saturating_add(max_distance);
    let near = |p: u64| p >= lo && p <= hi;

    // A bracket allele's breakend is one point, and it is the coordinate inside the
    // brackets. `StructuralVariationFeature::_parse_breakends` (`:1300-1316`) builds each
    // breakend from the parsed `$alt_chr`/`$alt_pos` alone:
    //
    //   my $slice = Bio::EnsEMBL::Slice->new_fast({
    //     seq_region_name => $alt_chr, start => $alt_pos, end => $alt_pos });
    //   push @$breakends, { string => $alt_string, ..., slice => $slice,
    //     chr => $slice->{seq_region_name}, start => $slice->{start}, end => $slice->{end} };
    //
    // so neither the local POS nor an `INFO/END` span end is ever the breakend of a paired bracket allele.
    if let Some(point) = paired_bnd_breakend_point(variant) {
        return near(point);
    }

    // Without a same-chromosome bracket coordinate there is no per-allele breakend to test.
    // A cross-chromosome mate is deliberately left on the span/POS tests below rather than
    // deselected: `_close_to_feature` rejects its bracket allele on
    // `_compare_seq_region_names`, and both engines annotate the cross-chromosome
    // class through the span/POS tests.
    //
    // A derived single breakend reaches here by construction, because `get_breakends`
    // returns `[]` for the `N.` form (`_parse_breakends` matches it on the
    // `^[\.]?([ACGTN]+)[\.]?$` arm and `next`s without pushing). Its single allele carries
    // no breakend at all, so `within_feature` falls back to `$bvfoa->base_variation_feature`
    // and the test is the whole span. Both spellings of its far end are checked.
    near(variant.start) || near(variant.sv_end.unwrap_or(variant.end))
}

/// The single coordinate Ensembl attaches to a paired BND's bracket allele, when that
/// coordinate is on the same chromosome as the transcript under test.
///
/// `within_feature` is called with `$match_seq_region_names = 1`, and
/// `_close_to_feature` opens with the same comparison, so a bracket allele whose
/// chromosome differs from the feature's is rejected outright rather than tested
/// positionally. Returning `None` in that case keeps both callers on their whole-span
/// test instead of silently deselecting cross-chromosome mates.
fn paired_bnd_breakend_point(variant: &InputVariant) -> Option<u64> {
    if variant.is_single_breakend {
        return None;
    }
    let same_chr = variant
        .mate_chr
        .as_deref()
        .is_some_and(|mc| mc == &*variant.chr);
    variant.mate_pos.filter(|_| same_chr)
}

/// Whether one of a BND's breakends lies within the transcript, which is the whole of
/// Ensembl's `feature_truncation` predicate for a breakend.
///
/// `Utils/VariationEffect.pm` has two arms and only the first can fire for a BND:
///
/// ```text
/// if(chromosome_breakpoint(@_)) {
///     return 1 if within_feature($bvfoa, $feat, $bvfo, $bvfoa->breakend, 1);
/// }
/// return 0 if $feat->isa('Bio::EnsEMBL::Transcript') and not within_cdna(@_);
/// return ( (partial_overlap_feature or complete_within_feature)
///          and (copy_number_loss(@_) or deletion(@_)) );
/// ```
///
/// `copy_number_loss` and `deletion` are both false for a `chromosome_breakpoint`, so the
/// second arm can never grant it to a BND. That leaves `within_feature` against the
/// allele's own breakend, with `$match_seq_region_names = 1` requiring the same
/// chromosome, and `within_feature` expands no slice.
///
/// The distinction from `bnd_breakend_close_to_transcript` is only the 5 kb expansion, so
/// this delegates rather than restating which coordinate is the breakend. Selection uses
/// the expanded form; `feature_truncation` uses this one.
/// VEP emits no `feature_truncation` for a mate one base before a transcript start, so the breakend is not the local POS.
pub(crate) fn bnd_breakend_within_transcript(
    variant: &InputVariant,
    transcript: &Transcript,
) -> bool {
    bnd_breakend_close_to_transcript(variant, transcript, 0)
}

pub(crate) fn overlaps_cds_exon(transcript: &Transcript, sv_start: u64, sv_end: u64) -> bool {
    let Some((cds_start, cds_end)) = crate::consequences::genomic_coding_bounds(transcript) else {
        return false;
    };

    transcript.exons.iter().any(|exon| {
        let coding_lo = exon.start.max(cds_start);
        let coding_hi = exon.end.min(cds_end);
        coding_lo <= coding_hi && overlap_bp(sv_start, sv_end, coding_lo, coding_hi) > 0
    })
}

pub(crate) fn transcript_has_incomplete_cds(transcript: &Transcript) -> bool {
    let facts = transcript.facts();
    facts.cds_start_nf || facts.cds_end_nf
}

pub(crate) fn overlaps_five_prime_utr(transcript: &Transcript, sv_start: u64, sv_end: u64) -> bool {
    let Some((cds_start, cds_end)) = crate::consequences::genomic_coding_bounds(transcript) else {
        return false;
    };

    let within_cdna = overlaps_any_exon(transcript, sv_start, sv_end);
    if !within_cdna {
        return false;
    }

    let before_coding = sv_end >= transcript.start && sv_start < cds_start;
    let after_coding = sv_end > cds_end && sv_start <= transcript.end;
    let start_boundary_is_nf = transcript.facts().cds_start_nf
        && cds_start == transcript.start
        && overlaps_point(sv_start, sv_end, transcript.start);
    let end_boundary_is_nf = transcript.facts().cds_end_nf
        && cds_end == transcript.end
        && overlaps_point(sv_start, sv_end, transcript.end);

    match transcript.strand {
        Strand::Forward => before_coding || start_boundary_is_nf,
        Strand::Reverse => after_coding || end_boundary_is_nf,
    }
}

pub(crate) fn overlaps_three_prime_utr(
    transcript: &Transcript,
    sv_start: u64,
    sv_end: u64,
) -> bool {
    let Some((cds_start, cds_end)) = crate::consequences::genomic_coding_bounds(transcript) else {
        return false;
    };

    let within_cdna = overlaps_any_exon(transcript, sv_start, sv_end);
    if !within_cdna {
        return false;
    }

    let before_coding = sv_end >= transcript.start && sv_start < cds_start;
    let after_coding = sv_end > cds_end && sv_start <= transcript.end;
    let start_boundary_is_nf = transcript.facts().cds_start_nf
        && cds_start == transcript.start
        && overlaps_point(sv_start, sv_end, transcript.start);
    let end_boundary_is_nf = transcript.facts().cds_end_nf
        && cds_end == transcript.end
        && overlaps_point(sv_start, sv_end, transcript.end);

    match transcript.strand {
        Strand::Forward => after_coding || end_boundary_is_nf,
        Strand::Reverse => before_coding || start_boundary_is_nf,
    }
}

/// A BND with a same-chromosome span (a derived single breakend up to `max_sv_size`,
/// or a paired bracket BND of any span): full positional annotation via the deletion
/// module, matching Perl's full StructuralVariationOverlap, which produces the richer
/// (coding_sequence_variant, intron_variant, UTR).
fn calculate_small_bnd_single_breakend(
    variant: &InputVariant,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    // Derived single-breakend BNDs are structural variants; `is_structural` is
    // preserved for deletion.rs.
    let bnd_as_del = variant.clone();
    let mut tc = deletion::calculate(
        &bnd_as_del,
        transcript,
        upstream_distance,
        downstream_distance,
    );

    if let Some(ref mut tc) = tc {
        // Perl's StructuralVariationOverlap assigns no stop_lost or start_lost to a
        // BND-format SV; the deletion module adds them for structural deletions.
        tc.consequences.retain(|c| {
            !matches!(
                c,
                Consequence::StopLost | Consequence::StartLost | Consequence::StartRetainedVariant
            )
        });

        // Post-process: Perl uses feature_truncation (not transcript_ablation) for BND-format SVs.
        if tc.consequences.contains(&Consequence::TranscriptAblation) {
            tc.consequences
                .retain(|c| *c != Consequence::TranscriptAblation);
            // Ensembl builds one allele per entry of `($vf, @{$vf->get_breakends})`,
            // each with `-breakend => $_`, and `feature_truncation`'s only live arm for
            // a `chromosome_breakpoint` is `within_feature($bvfoa->breakend)`. For the
            // `$vf` entry (`N.`) the breakend is the whole feature, so an engulfing
            // span passes; for a bracket entry it is one coordinate, so an engulfing
            // span with both endpoints outside the transcript fails and Ensembl emits
            // the context term alone.
            //
            let paired_breakend_outside =
                !variant.is_single_breakend && !bnd_breakend_within_transcript(variant, transcript);
            if !paired_breakend_outside {
                tc.consequences.push(Consequence::FeatureTruncation);
            }
            let has_context = tc.consequences.iter().any(|c| {
                matches!(
                    c,
                    Consequence::CodingTranscriptVariant
                        | Consequence::NonCodingTranscriptVariant
                        | Consequence::NmdTranscriptVariant
                        | Consequence::MatureMirnaVariant
                )
            });
            if !has_context {
                if transcript.has_cds() {
                    if transcript.is_nmd_transcript() {
                        tc.consequences.push(Consequence::NmdTranscriptVariant);
                    } else {
                        tc.consequences.push(Consequence::CodingTranscriptVariant);
                    }
                } else if is_mature_mirna_sv(transcript, variant.start, variant.end) {
                    tc.consequences.push(Consequence::MatureMirnaVariant);
                } else {
                    // Reached only when the variant engulfed the transcript (this block
                    // is gated on transcript_ablation). Perl's `non_coding_exon_variant`
                    // opens with `return 0 if complete_overlap_feature(@_)`
                    // (Utils/VariationEffect.pm:502-519) and `within_non_coding_gene`
                    // (`:495-500`) is its complement, so an engulfing variant receives
                    // the generic term.
                    tc.consequences
                        .push(Consequence::NonCodingTranscriptVariant);
                }
            }
            tc.consequences.sort_by_key(|c| c.rank());
        }
        // Perl does not emit feature_truncation alongside mature_miRNA_variant.
        if tc.consequences.contains(&Consequence::MatureMirnaVariant)
            && tc.consequences.contains(&Consequence::FeatureTruncation)
        {
            tc.consequences
                .retain(|c| *c != Consequence::FeatureTruncation);
            tc.consequences.sort_by_key(|c| c.rank());
        }

        // Perl emits feature_truncation for a BND-format SV that overlaps the
        // transcript body even on an intron-only overlap completely within it,
        // which the deletion module's exon-overlap requirement misses.
        let is_upstream_downstream_only = tc.consequences.iter().all(|c| {
            matches!(
                c,
                Consequence::UpstreamGeneVariant | Consequence::DownstreamGeneVariant
            )
        });
        // The same breakend condition as the ablation conversion above: for a
        // bracket allele Ensembl's only live `feature_truncation` arm is
        // `within_feature($bvfoa->breakend)` against that single coordinate, and
        // without it this push would undo the suppression above.
        let paired_breakend_outside_body =
            !variant.is_single_breakend && !bnd_breakend_within_transcript(variant, transcript);
        if !tc.consequences.contains(&Consequence::FeatureTruncation)
            && !tc.consequences.contains(&Consequence::TranscriptAblation)
            && !is_upstream_downstream_only
            && !tc.consequences.contains(&Consequence::MatureMirnaVariant)
            && !paired_breakend_outside_body
        {
            tc.consequences.push(Consequence::FeatureTruncation);
            tc.consequences.sort_by_key(|c| c.rank());
        }
        // Every edit above changes the term list; VEP's IMPACT is the most severe
        // impact of the terms a row finally carries, so it is derived once here.
        tc.impact = most_severe_impact(&tc.consequences);
    }
    tc
}

/// The most severe impact among a row's consequence terms.
pub(crate) fn most_severe_impact(consequences: &[Consequence]) -> Impact {
    consequences
        .iter()
        .map(|c| c.impact())
        .min_by_key(|i| match i {
            Impact::HIGH => 0,
            Impact::MODERATE => 1,
            Impact::LOW => 2,
            Impact::MODIFIER => 3,
        })
        .unwrap_or(Impact::MODIFIER)
}

/// Lightweight consequence calculator for derived single-breakend (N.) forms
/// from symbolic `<BND>` variants with SVLEN-based span.
///
/// For the N. allele this path emits `feature_truncation`, one biotype context term
/// (or the exon term for a partially covered non-coding transcript) and
/// `intron_variant`; it does not derive the deletion module's positional sub-terms
/// (coding_sequence_variant, stop_lost, UTR terms, splice).
///
/// The Perl predicate for BND feature_truncation (VariationEffect.pm) has a
/// special early-return: `if(chromosome_breakpoint(@_)) { return 1 if
/// within_feature($bvfoa, $feat, $bvfo, $bvfoa->breakend, 1); }`, checking
/// only that the breakend/span overlaps the transcript body.
fn calculate_derived_single_breakend(
    variant: &InputVariant,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    let sv_start = variant.start;
    let sv_end = variant.sv_end.unwrap_or(variant.end);
    let tx_start = transcript.start;
    let tx_end = transcript.end;

    let max_dist = upstream_distance.max(downstream_distance);
    if sv_end < tx_start.saturating_sub(max_dist) || sv_start > tx_end + max_dist {
        return None;
    }

    let overlaps_transcript = sv_end >= tx_start && sv_start <= tx_end;

    if !overlaps_transcript {
        let mut consequences: ConsequenceList = SmallVec::new();
        let distance;
        match transcript.strand {
            Strand::Forward => {
                if sv_end < tx_start {
                    let dist = tx_start - sv_end;
                    if dist <= upstream_distance {
                        consequences.push(Consequence::UpstreamGeneVariant);
                        distance = Some(dist);
                    } else {
                        return None;
                    }
                } else {
                    let dist = sv_start - tx_end;
                    if dist <= downstream_distance {
                        consequences.push(Consequence::DownstreamGeneVariant);
                        distance = Some(dist);
                    } else {
                        return None;
                    }
                }
            }
            Strand::Reverse => {
                if sv_start > tx_end {
                    let dist = sv_start - tx_end;
                    if dist <= upstream_distance {
                        consequences.push(Consequence::UpstreamGeneVariant);
                        distance = Some(dist);
                    } else {
                        return None;
                    }
                } else {
                    let dist = tx_start - sv_end;
                    if dist <= downstream_distance {
                        consequences.push(Consequence::DownstreamGeneVariant);
                        distance = Some(dist);
                    } else {
                        return None;
                    }
                }
            }
        }

        return Some(build_sv_consequence(transcript, consequences, distance));
    }

    // Perl's BND feature_truncation is within_feature(breakend), simple body
    // overlap, then one context term by biotype. Ensembl short-circuits the
    // specific sub-terms when the variant engulfs the transcript:
    // `non_coding_exon_variant` opens `return 0 if complete_overlap_feature(@_)`,
    // so an engulfing span gets the generic term and no positional sub-terms.
    let engulfs_transcript = sv_start <= transcript.start && sv_end >= transcript.end;

    let mut consequences: ConsequenceList = SmallVec::new();

    // Mature miRNA first: Perl emits only mature_miRNA_variant, no
    // feature_truncation alongside, for BND-format SVs.
    if is_mature_mirna_sv(transcript, sv_start, sv_end) {
        consequences.push(Consequence::MatureMirnaVariant);
    } else {
        // Body overlap is correct here and only here. This path serves the derived single
        // breakend, whose tuple carries the `N.` allele, and that is Ensembl's `$vf` entry
        // of `($vf, @{$vf->get_breakends})`: the allele is constructed with `-breakend => $_`,
        // so for the `$vf` entry the "breakend" passed to `within_feature` IS the whole
        // variation feature, spanning start to end. A paired bracket allele is the other
        // case and is handled in `deletion.rs`, where the breakend is a single coordinate.
        consequences.push(Consequence::FeatureTruncation);

        if transcript.has_cds() {
            if transcript.is_nmd_transcript() {
                consequences.push(Consequence::NmdTranscriptVariant);
            } else {
                consequences.push(Consequence::CodingTranscriptVariant);
            }
        } else {
            // The two non-coding forms are complements chosen by exon overlap, never
            // a hierarchy: `within_non_coding_gene` is the complement of
            // `non_coding_exon_variant` over `within_transcript`
            // (`Utils/VariationEffect.pm`).
            if !engulfs_transcript && overlaps_any_exon(transcript, sv_start, sv_end) {
                consequences.push(Consequence::NonCodingTranscriptExonVariant);
            } else {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
        }
    }

    // `intron_variant` is independent of the biotype context term and of the
    // exon/generic choice above.
    if !engulfs_transcript
        && overlaps_any_intron_trimmed(transcript, sv_start, sv_end)
        && !consequences.contains(&Consequence::IntronVariant)
    {
        consequences.push(Consequence::IntronVariant);
    }

    consequences.sort_by_key(|c| c.rank());

    Some(build_sv_consequence(transcript, consequences, None))
}

/// Build a `TranscriptConsequence` for an SV.
fn build_sv_consequence(
    transcript: &Transcript,
    consequences: ConsequenceList,
    distance: Option<u64>,
) -> TranscriptConsequence {
    use vep_core::consequence::FeatureType;

    let impact = consequences
        .first()
        .map(|c| c.impact())
        .unwrap_or(Impact::MODIFIER);

    TranscriptConsequence {
        transcript_id: transcript.stable_id.clone(),
        feature_start: transcript.start,
        feature_end: transcript.end,
        gene_id: transcript.gene_stable_id.clone(),
        gene_symbol: transcript.gene_symbol.clone(),
        gene_symbol_source: transcript.gene_symbol_source.clone(),
        hgnc_id: transcript.hgnc_id.clone(),
        consequences,
        impact,
        biotype: Some(transcript.biotype.clone()),
        canonical: transcript.canonical,
        strand: match transcript.strand {
            Strand::Forward => 1,
            Strand::Reverse => -1,
        },
        distance,
        feature_type: FeatureType::Transcript,
        flags: transcript.flags.clone(),
        tsl: transcript.tsl,
        ccds: transcript.ccds.clone(),
        swissprot: transcript.swissprot.clone(),
        trembl: transcript.trembl.clone(),
        refseq: transcript.refseq.clone(),
        ..TranscriptConsequence::default()
    }
}

/// Calculate consequences for a structural variant against a transcript.
///
/// Returns `None` if the SV doesn't overlap the transcript region at all
/// (or the SV type is unsupported: stubs return `None` to prevent garbage
/// output).
pub fn calculate_sv_consequences(
    variant: &InputVariant,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    match variant.variant_class {
        VariantClass::StructuralDeletion => {
            deletion::calculate(variant, transcript, upstream_distance, downstream_distance)
        }
        VariantClass::StructuralInsertion => {
            insertion::calculate(variant, transcript, upstream_distance, downstream_distance)
        }
        VariantClass::Duplication | VariantClass::TandemDuplication => {
            duplication::calculate(variant, transcript, upstream_distance, downstream_distance)
        }
        VariantClass::Inversion => {
            inversion::calculate(variant, transcript, upstream_distance, downstream_distance)
        }
        VariantClass::Translocation => {
            // Any BND carrying a span needs span-based annotation, because that is what
            // Ensembl VEP evaluates. `StructuralVariationOverlap::new` builds one allele
            // per entry of `($vf, @{$vf->get_breakends})`, but every positional predicate
            // then takes `$bvf = $bvfoa->base_variation_feature`, the local variation
            // feature, so `non_coding_exon_variant` tests
            // `overlap($bvf->{start}, $bvf->{end}, exon)` over the whole span for every
            // allele. Only `feature_truncation` consults the per-allele breakend. A
            // bracket-notation paired BND has `is_single_breakend = false`
            // (`vcf_parser.rs`) and must still enter here, or its mid-span
            // transcripts would be annotated by nothing and its mate-side ones only
            // at the mate point by `runner.rs::append_bnd_mate_consequences`.
            //
            // An inter-chromosomal mate leaves `end == start`, so it does not enter here
            // and keeps point evaluation, which is also what Ensembl does: both
            // `_close_to_feature` and `within_feature` require matching chromosome names.
            // Feature selection differs between the two BND shapes. A derived single
            // breakend (symbolic `<BND>` with `SVLEN`) has one breakend and a genuine
            // interval, so every transcript the interval covers is selected. A paired
            // bracket BND has two breakends and no deleted interval between them, so
            // Ensembl selects only transcripts a breakend is close to; the span still
            // governs the terms.
            if variant.end > variant.start
                && variant.mate_id.is_some()
                && (variant.is_single_breakend
                    || bnd_breakend_close_to_transcript(
                        variant,
                        transcript,
                        upstream_distance.max(downstream_distance),
                    ))
            {
                // Ensembl's `--max_sv_size` decides whether an SV is annotated at all,
                // never which terms it gets. The cheap path below is an order of
                // magnitude faster than `deletion::calculate` on giant spans and emits
                // `feature_truncation`, a biotype context term and `intron_variant` from
                // two interval tests, the set the full path yields for an engulfing
                // span or a non-coding transcript.
                let sv_span = variant
                    .sv_end
                    .unwrap_or(variant.end)
                    .saturating_sub(variant.start);
                if sv_span > PERL_DEFAULT_MAX_SV_SIZE && variant.is_single_breakend {
                    calculate_derived_single_breakend(
                        variant,
                        transcript,
                        upstream_distance,
                        downstream_distance,
                    )
                } else {
                    calculate_small_bnd_single_breakend(
                        variant,
                        transcript,
                        upstream_distance,
                        downstream_distance,
                    )
                }
            } else {
                breakend::calculate(variant, transcript, upstream_distance, downstream_distance)
            }
        }
        VariantClass::CopyNumberVariation => {
            cnv::calculate(variant, transcript, upstream_distance, downstream_distance)
        }
        VariantClass::MobileElementInsertion => {
            insertion::calculate(variant, transcript, upstream_distance, downstream_distance)
        }
        VariantClass::MobileElementDeletion => {
            deletion::calculate(variant, transcript, upstream_distance, downstream_distance)
        }
        // Complex SVs (<CPX>): Despite not being in Perl's %SO_TERMS, CPX variants
        // are still annotated in default VEP tab output (vep_skip flag is only checked
        // in File-based annotation sources, not in Cache annotation). Use CNV generic
        // path for overlap-only consequences matching Perl's effective annotation.
        VariantClass::ComplexStructural => {
            cnv::calculate(variant, transcript, upstream_distance, downstream_distance)
        }
        // A `<CNV:TR>` record states its alternate allele's length (INFO/RB, or RUC times
        // the unit length) against the reference run SVLEN spans. Perl reads a symbolic
        // tandem repeat as a copy-number gain (`copy_number_gain` in VariationEffect.pm),
        // which is the insertion path here: `feature_elongation`, and
        // `transcript_amplification` when the run contains the transcript. An alternate
        // allele shorter than the run is a loss, which Perl's `copy_number_loss` reads
        // through `feature_truncation`, the deletion path. A record without the fields
        // states no direction and takes the gain path.
        VariantClass::TandemRepeat => {
            let run_bases = variant
                .sv_end
                .unwrap_or(variant.end)
                .saturating_sub(variant.start)
                .saturating_add(1);
            match variant.tr_alt_bases {
                Some(alt_bases) if alt_bases < run_bases => {
                    deletion::calculate(variant, transcript, upstream_distance, downstream_distance)
                }
                _ => insertion::calculate(
                    variant,
                    transcript,
                    upstream_distance,
                    downstream_distance,
                ),
            }
        }
        // Non-structural variants should never reach here.
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::transcript::{
        Attribute, ExonCoordMapper, MapperPair, TranscriptMapper, TranscriptVEFC,
    };
    use vep_core::variant::InputVariant;

    /// Helper: build a minimal structural variant for testing dispatch.
    fn make_sv(variant_class: VariantClass, sv_end: Option<u64>) -> InputVariant {
        let mut v = InputVariant::new(
            "21".into(),
            25_000_000,
            25_006_000,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = variant_class;
        v.is_structural = true;
        v.sv_end = sv_end;
        v
    }

    #[test]
    fn test_dispatch_structural_deletion() {
        let v = make_sv(VariantClass::StructuralDeletion, Some(25_100_000));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "Deletion overlapping transcript should produce consequences"
        );
    }

    #[test]
    fn test_dispatch_structural_insertion() {
        let v = make_sv(VariantClass::StructuralInsertion, Some(25_000_500));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "Insertion within transcript should produce consequences"
        );
    }

    #[test]
    fn test_dispatch_duplication() {
        let v = make_sv(VariantClass::Duplication, Some(25_200_000));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "Duplication overlapping transcript should produce consequences"
        );
    }

    #[test]
    fn test_dispatch_tandem_duplication() {
        let v = make_sv(VariantClass::TandemDuplication, Some(25_200_000));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "Tandem dup overlapping transcript should produce consequences"
        );
    }

    #[test]
    fn test_dispatch_inversion() {
        let v = make_sv(VariantClass::Inversion, Some(25_200_000));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "Inversion overlapping transcript should produce consequences"
        );
    }

    #[test]
    fn test_dispatch_translocation() {
        let v = make_sv(VariantClass::Translocation, None);
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "BND at transcript start should produce consequences"
        );
    }

    #[test]
    fn test_dispatch_cnv() {
        let v = make_sv(VariantClass::CopyNumberVariation, Some(25_200_000));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "CNV overlapping transcript should produce consequences"
        );
    }

    #[test]
    fn test_dispatch_mobile_element_insertion() {
        let v = make_sv(VariantClass::MobileElementInsertion, Some(25_000_500));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "ME insertion within transcript should produce consequences"
        );
    }

    #[test]
    fn test_dispatch_mobile_element_deletion() {
        let v = make_sv(VariantClass::MobileElementDeletion, Some(25_200_000));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "ME deletion overlapping transcript should produce consequences"
        );
    }

    #[test]
    fn test_dispatch_complex_structural() {
        // Perl VEP annotates <CPX> with overlap-only consequences (despite vep_skip).
        // Route through CNV generic path for regional overlap consequences.
        let v = make_sv(VariantClass::ComplexStructural, Some(25_002_100));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "Complex SV (<CPX>) should produce consequences (Perl VEP annotates these)"
        );
    }

    #[test]
    fn test_dispatch_tandem_repeat() {
        let v = make_sv(VariantClass::TandemRepeat, Some(25_000_500));
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "Tandem repeat within transcript should produce consequences"
        );
    }

    /// A `<CNV:TR>` over 300 reference bases inside the first exon and intron of the test
    /// transcript, with the alternate allele's length set by the caller.
    fn make_tandem_repeat(tr_alt_bases: Option<u64>) -> InputVariant {
        let mut v = InputVariant::new(
            "21".into(),
            25_000_101,
            25_000_400,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::TandemRepeat;
        v.is_structural = true;
        v.sv_end = Some(25_000_400);
        v.tr_alt_bases = tr_alt_bases;
        v
    }

    #[test]
    fn test_tandem_repeat_gain_takes_the_insertion_path() {
        // RB above the 300-base run: a copy-number gain, feature_elongation.
        let tx = crate::test_helpers::make_test_transcript();
        let tc = calculate_sv_consequences(&make_tandem_repeat(Some(360)), &tx, 5000, 5000)
            .expect("gain overlapping the transcript annotates");
        assert!(
            tc.consequences.contains(&Consequence::FeatureElongation),
            "gain should carry feature_elongation, got {:?}",
            tc.consequences
        );
        assert!(!tc.consequences.contains(&Consequence::FeatureTruncation));
    }

    #[test]
    fn test_tandem_repeat_loss_takes_the_deletion_path() {
        // RB below the run: a contraction, feature_truncation through the deletion path.
        let tx = crate::test_helpers::make_test_transcript();
        let tc = calculate_sv_consequences(&make_tandem_repeat(Some(240)), &tx, 5000, 5000)
            .expect("loss overlapping the transcript annotates");
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "loss should carry feature_truncation, got {:?}",
            tc.consequences
        );
        assert!(!tc.consequences.contains(&Consequence::FeatureElongation));
    }

    #[test]
    fn test_tandem_repeat_without_alt_length_or_equal_keeps_the_gain_path() {
        // No RB, RUC or RUL, or an alternate allele the length of the run: no direction is
        // stated, and the record reads as Perl's symbolic tandem repeat, a gain.
        let tx = crate::test_helpers::make_test_transcript();
        for alt in [None, Some(300)] {
            let tc = calculate_sv_consequences(&make_tandem_repeat(alt), &tx, 5000, 5000)
                .expect("tandem repeat overlapping the transcript annotates");
            assert!(
                tc.consequences.contains(&Consequence::FeatureElongation),
                "tr_alt_bases {alt:?} should take the gain path, got {:?}",
                tc.consequences
            );
        }
    }

    #[test]
    fn test_noncoding_del_keeps_specific_terms() {
        // Non-coding structural deletion that extends beyond the transcript should
        // get feature_truncation + specific sub-terms. The non-coding normalization
        // only applies to BND-derived paths, not direct deletions.
        // Note: A deletion that fully spans the transcript produces transcript_ablation,
        // not feature_truncation. Use a partial-spanning DEL (extends_beyond case).
        let mut tx = crate::test_helpers::make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;

        // DEL from before transcript start to mid-intron1 (extends beyond, not ablation)
        let mut v = InputVariant::new(
            "21".into(),
            24_999_000,
            25_001_500,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::StructuralDeletion;
        v.is_structural = true;
        v.sv_end = Some(25_001_500);

        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        // Direct deletion: Perl uses specific sub-terms.
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "Should have feature_truncation, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_non_structural_returns_none() {
        let mut v = InputVariant::new(
            "21".into(),
            25_000_100,
            25_000_100,
            b"A".to_vec(),
            b"G".to_vec(),
        );
        // Force it through the SV path even though it is an SNV: it hits the _ => None arm
        v.variant_class = VariantClass::Snv;
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_none(),
            "Non-structural variant should return None from SV dispatch"
        );
    }

    /// Build a non-coding miRNA transcript for testing is_mature_mirna_sv.
    ///
    /// Single exon: genomic 1000-1100, cDNA 1-101.
    /// miRNA attribute "50-70" = cDNA positions 50-70 → genomic 1049-1069 (forward strand).
    fn make_mirna_transcript() -> Transcript {
        use std::sync::Arc;
        use vep_core::transcript::Exon;

        let exons = vec![Exon {
            stable_id: Some("ENSE_MIRNA_001".into()),
            start: 1000,
            end: 1100,
            rank: 1,
            phase: -1,
            end_phase: -1,
        }];

        let mapper_pairs = vec![MapperPair {
            from_start: 1,
            from_end: 101,
            to_start: 1000,
            to_end: 1100,
            ori: 1,
        }];

        let vefc = TranscriptVEFC {
            codon_table: 1,
            five_prime_utr: None,
            three_prime_utr: None,
            translateable_seq: None,
            peptide: None,
            introns: vec![],
            sorted_exons: exons.clone(),
            mapper: Some(TranscriptMapper {
                start_phase: 0,
                cdna_coding_start: 0,
                cdna_coding_end: 0,
                exon_coord_mapper: ExonCoordMapper::new(mapper_pairs),
            }),
            protein_features: vec![],
            protein_function_predictions: None,
            seq_edits: vec![],
        };

        Transcript {
            stable_id: "ENST_MIRNA_001".into(),
            version: Some(1),
            db_id: None,
            gene_stable_id: "ENSG_MIRNA_001".into(),
            chr: "21".into(),
            start: 1000,
            end: 1100,
            strand: Strand::Forward,
            biotype: "miRNA".into(),
            source: "Ensembl".into(),
            description: None,
            gene_symbol: Some("MIR21".into()),
            gene_symbol_source: Some("HGNC".into()),
            hgnc_id: None,
            gene_phenotype: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            exons,
            introns: vec![],
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            translation: None,
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from([]),
            gencode_primary: false,
            attributes: vec![Attribute {
                code: "miRNA".into(),
                value: "50-70".into(),
            }],
            vefc: Some(vefc),
            derived: Default::default(),
        }
    }

    #[test]
    fn test_is_mature_mirna_sv_with_overlap() {
        let tx = make_mirna_transcript();
        // miRNA attribute "50-70" → cDNA 50-70 → genomic 1049-1069 (forward strand).
        // SV 1050-1060 overlaps the mature region.
        assert!(
            is_mature_mirna_sv(&tx, 1050, 1060),
            "SV overlapping mature miRNA region should return true"
        );
    }

    #[test]
    fn test_is_mature_mirna_sv_no_overlap() {
        let tx = make_mirna_transcript();
        // miRNA attribute "50-70" → genomic 1049-1069.
        // SV 1000-1040 is entirely before the mature region.
        assert!(
            !is_mature_mirna_sv(&tx, 1000, 1040),
            "SV not overlapping mature miRNA region should return false"
        );
    }

    #[test]
    fn test_is_mature_mirna_sv_wrong_biotype() {
        let mut tx = make_mirna_transcript();
        tx.biotype = "lncRNA".into();
        // Even with miRNA attribute, non-miRNA biotype should return false.
        assert!(
            !is_mature_mirna_sv(&tx, 1050, 1060),
            "Non-miRNA biotype should return false"
        );
    }

    #[test]
    fn test_parse_mirna_cdna_range() {
        assert_eq!(parse_mirna_cdna_range("50-70"), Some((50, 70)));
        assert_eq!(parse_mirna_cdna_range("  50-70 "), Some((50, 70)));
        assert_eq!(parse_mirna_cdna_range("abc"), None);
        assert_eq!(parse_mirna_cdna_range(""), None);
    }

    /// Helper: build a derived single-breakend variant for testing.
    /// Giant BND: span > 10M, uses lightweight calculator.
    fn make_giant_derived_sb(start: u64, sv_end: u64) -> InputVariant {
        let mut v = InputVariant::new("21".into(), start, sv_end, b"N".to_vec(), b"-".to_vec());
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = true;
        v.sv_end = Some(sv_end);
        v.mate_id = Some("mate_001".into());
        v
    }

    #[test]
    fn test_giant_derived_sb_spanning_transcript_coding() {
        // Giant BND spanning entire coding transcript → feature_truncation + coding_transcript_variant
        let v = make_giant_derived_sb(24_000_000, 45_000_000); // 21M span, above 10M
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "Giant BND spanning transcript should annotate"
        );
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "Should have feature_truncation"
        );
        assert!(
            tc.consequences
                .contains(&Consequence::CodingTranscriptVariant),
            "Should have coding_transcript_variant"
        );
        // Must not have deletion-like sub-terms
        assert!(
            !tc.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Giant BND should not have coding_sequence_variant"
        );
        assert!(
            !tc.consequences.contains(&Consequence::StopLost),
            "Giant BND should not have stop_lost"
        );
        assert!(
            !tc.consequences.contains(&Consequence::IntronVariant),
            "Giant BND should not have intron_variant"
        );
        assert!(
            !tc.consequences.contains(&Consequence::TranscriptAblation),
            "Giant BND should not have transcript_ablation"
        );
    }

    #[test]
    fn test_giant_derived_sb_upstream() {
        // Giant BND span upstream of transcript → upstream_gene_variant
        let tx = crate::test_helpers::make_test_transcript();
        let v = make_giant_derived_sb(10_000_000, 24_993_000); // ends before tx_start - 5000
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        // Span ends at 24_993_000, tx starts at ~25_000_000, distance = 7000 > upstream_distance
        assert!(
            result.is_none(),
            "Giant BND ending >5000bp before transcript should return None"
        );
    }

    #[test]
    fn test_giant_derived_sb_mirna_no_feature_truncation() {
        // Giant BND overlapping miRNA transcript → mature_miRNA_variant only
        let tx = make_mirna_transcript(); // start=1000, end=1100
        let v = make_giant_derived_sb(500, 10_001_500); // 10M+ span, overlaps miRNA
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(
            result.is_some(),
            "Giant BND overlapping miRNA should annotate"
        );
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::MatureMirnaVariant),
            "Should have mature_miRNA_variant"
        );
        assert!(
            !tc.consequences.contains(&Consequence::FeatureTruncation),
            "miRNA should not have feature_truncation alongside"
        );
    }

    #[test]
    fn test_small_derived_sb_uses_deletion_path() {
        // Small BND (< 10M): should still produce deletion-like consequences
        let mut v = InputVariant::new(
            "21".into(),
            25_000_000,
            25_005_000,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = true;
        v.sv_end = Some(25_005_000); // 5K span, below 10M
        v.mate_id = Some("mate_002".into());
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(result.is_some(), "Small BND should annotate");
        let tc = result.unwrap();
        // Small BNDs go through deletion path, so they can have sub-terms
        assert!(
            !tc.consequences.contains(&Consequence::TranscriptAblation),
            "Small BND should not have transcript_ablation (mapped to feature_truncation)"
        );
    }

    #[test]
    fn test_paired_bnd_without_a_span_uses_the_breakend_path() {
        // A paired BND carrying NO span still takes the point-based path. What routes it
        // there is `end == start` and an unset `mate_id`, not `is_single_breakend`: an
        // inter-chromosomal mate leaves the Location a bare position, which is also what
        // Ensembl VEP reports for one.
        let mut v = InputVariant::new(
            "21".into(),
            25_000_000,
            25_000_000,
            b"N".to_vec(),
            b"N[21:30000000[".as_slice().to_vec(),
        );
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = false; // paired, not single
        let tx = crate::test_helpers::make_test_transcript();
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        // Paired BND at transcript start should produce breakend consequences
        assert!(
            result.is_some(),
            "Paired BND at transcript position should produce consequences"
        );
    }

    #[test]
    fn test_paired_bnd_engulfing_a_transcript_with_both_ends_outside_gets_no_truncation() {
        // Ensembl's `feature_truncation` has one live arm for a `chromosome_breakpoint`,
        // `within_feature($bvfoa->breakend)`, because its second arm also needs
        // `copy_number_loss or deletion` and neither holds for a BND. For a bracket allele
        // the breakend is one coordinate, so an engulfing span whose endpoints both fall
        // outside the transcript earns the context term and no truncation.
        //
        // The fixture has to pass selection and fail truncation, which are different
        // predicates: selection expands the transcript by 5 kb, truncation does not. POS at
        // 24,996,000 is inside the expanded region and outside the transcript
        // (25,000,000-25,006,000), and neither the mate at 24,995,999 nor the span end at
        // 25,100,000 is inside it either.
        //
        let tx = make_non_coding_transcript();
        let mut v = InputVariant::new(
            "21".into(),
            24_996_000,
            25_100_000,
            b"N".to_vec(),
            b"N[21:24995999[".as_slice().to_vec(),
        );
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = false;
        v.mate_id = Some("mate_engulf_001".into());
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(24_995_999);

        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("an engulfing BND still annotates the transcript it engulfs");
        assert!(
            !tc.consequences.contains(&Consequence::FeatureTruncation),
            "no breakend lies inside the transcript, so feature_truncation is not owed: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "the context term is still owed via within_non_coding_gene: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_intron_touched_only_in_its_boundary_bases_gets_no_intron_variant() {
        // Ensembl sets its `intronic` flag with
        // `overlap($r_start, $r_end, $intron_start + 2, $intron_end - 2)`, attributing the
        // four invariant donor and acceptor bases to the splice terms instead, and
        // `within_intron` returns exactly that flag.
        //
        // The span covers exon 1 (25,000,000-25,000,299) and only the first two bases of
        // intron 1 (25,000,300-25,001,999), so the exon term is owed and `intron_variant` is
        // not; under the untrimmed test this fixture emits `intron_variant`.
        let tx = make_non_coding_transcript();
        let mut v = InputVariant::new(
            "21".into(),
            25_000_290,
            25_000_301,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::StructuralDeletion;
        v.is_structural = true;
        v.sv_end = Some(25_000_301);

        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("a deletion overlapping exon 1 must annotate the transcript");
        assert!(
            !tc.consequences.contains(&Consequence::IntronVariant),
            "only the two boundary bases of intron 1 are touched, so intron_variant is not \
             owed: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "exon 1 is overlapped, so the exon term is owed: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_paired_bnd_with_a_same_chromosome_span_evaluates_the_span_not_the_mate_point() {
        // The fixture is built so span and point evaluations disagree. Span
        // 24,990,000-25,001,000 covers exon 1 (25,000,000-25,000,299) and reaches into
        // intron 1 (25,000,300-25,001,999), so a span evaluation owes the exon term; the
        // mate point 25,001,000 lies inside intron 1 and inside NO exon, so a point
        // evaluation yields the generic parent instead.
        //
        // `is_single_breakend` is false here, as `vcf_parser.rs` sets it for
        // bracket-notation paired BNDs, so the routing gate must not require it.
        let tx = make_non_coding_transcript();
        let mut v = InputVariant::new(
            "21".into(),
            24_990_000,
            25_001_000,
            b"N".to_vec(),
            b"N[21:25001000[".as_slice().to_vec(),
        );
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = false;
        v.mate_id = Some("mate_paired_span_001".into());
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(25_001_000);

        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("a paired BND whose span overlaps a transcript must annotate it");
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "the span covers exon 1, so the exon term is owed; a mate-point evaluation \
             would emit the generic parent instead: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "the generic parent is the complement of the exon term, never its companion: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "the span also reaches intron 1: {:?}",
            tc.consequences
        );
    }

    /// Builds a paired BND whose bracket mate is on chr21 at `mate_pos`, spanning
    /// `start..=end`. The mate coordinate is the only breakend Ensembl attaches to the
    /// bracket allele (`StructuralVariationFeature::_parse_breakends:1300-1316` sets the
    /// breakend's `chr`/`start`/`end` from the parsed bracket alone), so a fixture that
    /// varies it independently of the span is what separates the two predicates.
    fn make_paired_bnd(start: u64, end: u64, mate_pos: u64) -> InputVariant {
        let mut v = InputVariant::new(
            "21".into(),
            start,
            end,
            b"N".to_vec(),
            format!("N[21:{mate_pos}[").into_bytes(),
        );
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = false;
        v.mate_id = Some("mate_fixture".into());
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(mate_pos);
        v
    }

    /// `feature_truncation` for a paired BND tests the bracket coordinate, never the span
    /// end: a union of the local POS, the span end and the bracket coordinate over-calls
    /// wherever the span end is inside the transcript while the bracket is outside.
    #[test]
    fn test_paired_bnd_span_end_inside_transcript_earns_no_truncation() {
        let tx = make_non_coding_transcript();
        // Transcript 25,000,000-25,006,000. Bracket at 24,996,000 is outside it but within
        // the 5 kb selection window; the span end 25,000,500 is inside it.
        let v = make_paired_bnd(24_996_001, 25_000_500, 24_996_000);
        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("the bracket coordinate is within 5 kb, so the transcript is selected");
        assert!(
            !tc.consequences.contains(&Consequence::FeatureTruncation),
            "the span end being inside the transcript is not a breakend being inside it: {:?}",
            tc.consequences
        );
    }

    /// IMPACT follows the terms the row finally carries. A paired BND whose span
    /// crosses a coding transcript's stop codon collects `stop_lost` from the deletion
    /// module, which the BND rules then drop, and whose bracket lies outside the
    /// transcript earns no `feature_truncation`; what remains is MODIFIER, as VEP prints.
    #[test]
    fn test_paired_bnd_impact_follows_the_terms_it_keeps() {
        let tx = crate::test_helpers::make_test_transcript();
        assert_eq!(tx.translation_end, Some(25_004_299), "fixture assumption");
        // Span 25,001,000-25,005,000 crosses the CDS end inside the transcript
        // (25,000,000-25,006,000); the bracket at 24,996,000 is outside it.
        let v = make_paired_bnd(25_001_000, 25_005_000, 24_996_000);
        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("the span overlaps the transcript, so it is annotated");
        assert!(
            !tc.consequences.contains(&Consequence::FeatureTruncation)
                && !tc.consequences.contains(&Consequence::StopLost),
            "{:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "{:?}",
            tc.consequences
        );
        assert_eq!(
            tc.impact,
            crate::sv::most_severe_impact(&tc.consequences),
            "IMPACT must be derived from the terms kept: {:?}",
            tc.consequences
        );
        assert_eq!(tc.impact, Impact::MODIFIER);
    }

    /// The decisive discriminator: with the bracket mate at POS-1 and `transcript.start ==
    /// POS`, the local POS is inside the transcript and the breakend is one base before it.
    ///
    /// A union predicate emits `feature_truncation` here; Ensembl VEP does not, which
    /// is what proves the breakend is the bracket coordinate rather than the local
    /// POS. gnomAD SV writes the bracket mate at POS-1, so the fixture has the shape of
    /// real records rather than a contrived one.
    #[test]
    fn test_paired_bnd_mate_one_base_before_the_transcript_earns_no_truncation() {
        let tx = make_non_coding_transcript();
        assert_eq!(tx.start, 25_000_000, "fixture assumption");
        let v = make_paired_bnd(25_000_000, 25_010_000, 24_999_999);
        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("the bracket coordinate is one base away, so the transcript is selected");
        assert!(
            !tc.consequences.contains(&Consequence::FeatureTruncation),
            "the breakend is one base BEFORE the transcript, so truncation is not owed even \
             though POS is inside it: {:?}",
            tc.consequences
        );
    }

    /// A paired BND still earns `feature_truncation` when the bracket coordinate IS inside
    /// the transcript, so the narrowing above is not a blanket suppression.
    #[test]
    fn test_paired_bnd_mate_inside_the_transcript_still_earns_truncation() {
        let tx = make_non_coding_transcript();
        let v = make_paired_bnd(24_996_001, 25_010_000, 25_003_000);
        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("the bracket coordinate is inside the transcript");
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "a breakend inside the transcript is the whole of Ensembl's live arm: {:?}",
            tc.consequences
        );
    }

    /// selection for a paired BND is by the bracket coordinate plus
    /// `MAX_DISTANCE_FROM_TRANSCRIPT`, not by span overlap.
    ///
    /// `_close_to_feature` is called per allele on each entry of `($vf, @{$vf->get_breakends})`,
    /// and the bracket allele's entry is one point, so selecting on the span or the
    /// local POS as well over-selects transcripts whose bracket coordinate is farther
    /// than 5 kb away.
    #[test]
    fn test_paired_bnd_whose_mate_is_beyond_5kb_is_not_annotated_at_all() {
        let tx = make_non_coding_transcript();
        // The span reaches deep into the transcript, but the bracket coordinate is 20 kb
        // upstream of it, well outside the 5 kb window.
        let v = make_paired_bnd(24_980_001, 25_000_500, 24_980_000);
        assert!(
            calculate_sv_consequences(&v, &tx, 5000, 5000).is_none(),
            "no breakend of this allele is close to the transcript, so Ensembl VEP emits no \
             row for it at all"
        );
        // With the bracket coordinate moved inside the 5 kb window and nothing else changed,
        // the transcript IS selected, so the assertion above is about the mate and not about
        // some other property of the fixture.
        let near = make_paired_bnd(24_980_001, 25_000_500, 24_996_000);
        assert!(
            calculate_sv_consequences(&near, &tx, 5000, 5000).is_some(),
            "a bracket coordinate 4 kb upstream is within MAX_DISTANCE_FROM_TRANSCRIPT"
        );
    }

    /// Ensembl's tier-2 short-circuit: a mature-miRNA hit is emitted alone.
    ///
    /// `@SORTED_OVERLAP_CONSEQUENCES` is sorted by tier alone
    /// (`BaseVariationFeatureOverlapAllele.pm:69`), and `get_all_OverlapConsequences`
    /// (`:266,275-278`) stops at the first consequence of a higher tier once a tier-2 match is
    /// recorded. `mature_miRNA_variant` is tier 2, `feature_truncation` is tier 3, so rank
    /// does not save the latter despite being lower.
    ///
    /// Without the `is_mature_mirna_sv` early return in `deletion::calculate` this
    /// fails with `[FeatureTruncation, MatureMirnaVariant]`.
    #[test]
    fn test_deletion_partially_overlapping_a_mature_mirna_gives_that_term_alone() {
        let tx = make_mature_mirna_transcript();
        // Transcript 25,000,000-25,000,100 with the mature window at 25,000,049-25,000,069.
        // The deletion starts before the transcript and ends inside the mature window, so it
        // partially overlaps: `transcript_ablation` (tier 1) does not apply and every tier-3
        // term would otherwise be evaluated.
        let mut v = InputVariant::new(
            "21".into(),
            24_999_500,
            25_000_060,
            b"N".to_vec(),
            b"<DEL>".to_vec(),
        );
        v.variant_class = VariantClass::StructuralDeletion;
        v.is_structural = true;
        v.sv_end = Some(25_000_060);

        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("a deletion reaching the mature window must annotate the transcript");
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::MatureMirnaVariant],
            "a mature-miRNA hit stops Ensembl's loop before every tier-3 term: {:?}",
            tc.consequences
        );
    }

    // overlaps_polypyrimidine_tract tests
    //
    // Forward-strand transcript layout:
    //   Intron 1: 25_000_300 - 25_001_999
    //     PPT (fwd) = [intron_end - 16, intron_end - 2] = [25_001_983, 25_001_997]
    //   Intron 2: 25_002_300 - 25_003_999
    //     PPT (fwd) = [25_003_983, 25_003_997]

    #[test]
    fn test_ppt_helper_overlapping_window() {
        let tx = crate::test_helpers::make_test_transcript();
        assert!(overlaps_polypyrimidine_tract(&tx, 25_001_990, 25_001_995));
    }

    #[test]
    fn test_ppt_helper_no_overlap_intron_body() {
        let tx = crate::test_helpers::make_test_transcript();
        assert!(!overlaps_polypyrimidine_tract(&tx, 25_000_500, 25_000_700));
    }

    #[test]
    fn test_ppt_helper_boundary_just_before() {
        let tx = crate::test_helpers::make_test_transcript();
        // One base before PPT start of intron 1
        assert!(!overlaps_polypyrimidine_tract(&tx, 25_001_970, 25_001_982));
    }

    #[test]
    fn test_ppt_helper_boundary_at_start() {
        let tx = crate::test_helpers::make_test_transcript();
        // Exactly at PPT start of intron 1
        assert!(overlaps_polypyrimidine_tract(&tx, 25_001_983, 25_001_983));
    }

    #[test]
    fn test_ppt_helper_reverse_strand() {
        // For reverse strand, acceptor is at intron_start.
        // PPT = [intron_start + 2, intron_start + 16]
        // Intron 1 start = 25_000_300, so reverse PPT = [25_000_302, 25_000_316]
        let mut tx = crate::test_helpers::make_test_transcript();
        tx.strand = Strand::Reverse;
        assert!(overlaps_polypyrimidine_tract(&tx, 25_000_310, 25_000_315));
        // Must not overlap: the forward PPT window of intron 1 is near the intron end
        // but the reverse PPT is near the intron start
        assert!(!overlaps_polypyrimidine_tract(&tx, 25_001_990, 25_001_995));
    }

    #[test]
    fn test_ppt_helper_intron2_overlap() {
        let tx = crate::test_helpers::make_test_transcript();
        // PPT of intron 2 = [25_003_983, 25_003_997]
        assert!(overlaps_polypyrimidine_tract(&tx, 25_003_990, 25_003_995));
    }

    // Non-coding exon variant tests for giant/small BND paths
    //
    // Test transcript (non-coding, same exon/intron layout as standard):
    //   Exon 1: 25_000_000 - 25_000_299
    //   Intron 1: 25_000_300 - 25_001_999
    //   Exon 2: 25_002_000 - 25_002_299
    //   Intron 2: 25_002_300 - 25_003_999
    //   Exon 3: 25_004_000 - 25_006_000

    /// Build a non-coding transcript with the standard 3-exon layout.
    /// Same genomic coordinates as `make_test_transcript()` but biotype=lncRNA,
    /// no CDS, no translation.
    fn make_non_coding_transcript() -> Transcript {
        use std::sync::Arc;
        use vep_core::transcript::{Exon, Intron};

        let exons = vec![
            Exon {
                stable_id: Some("ENSE_NC_001".into()),
                start: 25_000_000,
                end: 25_000_299,
                rank: 1,
                phase: -1,
                end_phase: -1,
            },
            Exon {
                stable_id: Some("ENSE_NC_002".into()),
                start: 25_002_000,
                end: 25_002_299,
                rank: 2,
                phase: -1,
                end_phase: -1,
            },
            Exon {
                stable_id: Some("ENSE_NC_003".into()),
                start: 25_004_000,
                end: 25_006_000,
                rank: 3,
                phase: -1,
                end_phase: -1,
            },
        ];

        let introns = vec![
            Intron {
                start: 25_000_300,
                end: 25_001_999,
                rank: 1,
            },
            Intron {
                start: 25_002_300,
                end: 25_003_999,
                rank: 2,
            },
        ];

        let mapper_pairs = vec![
            MapperPair {
                from_start: 1,
                from_end: 300,
                to_start: 25_000_000,
                to_end: 25_000_299,
                ori: 1,
            },
            MapperPair {
                from_start: 301,
                from_end: 600,
                to_start: 25_002_000,
                to_end: 25_002_299,
                ori: 1,
            },
            MapperPair {
                from_start: 601,
                from_end: 2601,
                to_start: 25_004_000,
                to_end: 25_006_000,
                ori: 1,
            },
        ];

        let vefc = TranscriptVEFC {
            codon_table: 1,
            five_prime_utr: None,
            three_prime_utr: None,
            translateable_seq: None,
            peptide: None,
            introns: introns.clone(),
            sorted_exons: exons.clone(),
            mapper: Some(TranscriptMapper {
                start_phase: 0,
                cdna_coding_start: 0,
                cdna_coding_end: 0,
                exon_coord_mapper: ExonCoordMapper::new(mapper_pairs),
            }),
            protein_features: vec![],
            protein_function_predictions: None,
            seq_edits: vec![],
        };

        Transcript {
            stable_id: "ENST_NC_001".into(),
            version: Some(1),
            db_id: None,
            gene_stable_id: "ENSG_NC_001".into(),
            chr: "21".into(),
            start: 25_000_000,
            end: 25_006_000,
            strand: Strand::Forward,
            biotype: "lncRNA".into(),
            source: "Ensembl".into(),
            description: None,
            gene_symbol: Some("TESTNC1".into()),
            gene_symbol_source: Some("HGNC".into()),
            hgnc_id: None,
            gene_phenotype: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            exons: exons.clone(),
            introns,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            translation: None,
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from([]),
            gencode_primary: false,
            attributes: vec![],
            vefc: Some(vefc),
            derived: Default::default(),
        }
    }

    #[test]
    fn test_giant_bnd_non_coding_exon_overlap() {
        // An engulfing giant BND gets feature_truncation + the generic
        // non_coding_transcript_variant: Perl emits no exon/intron sub-terms for an
        // engulfing span.
        let tx = make_non_coding_transcript();
        let v = make_giant_derived_sb(24_000_000, 45_000_000); // 21M span
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(result.is_some(), "Giant BND should annotate non-coding tx");
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "Should have feature_truncation: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "Should have non_coding_transcript_variant (Perl parity): {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "Should not have non_coding_transcript_exon_variant (Perl uses generic term): {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences.contains(&Consequence::IntronVariant),
            "Should not have intron_variant (Perl uses generic term): {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_giant_bnd_non_coding_intron_only() {
        // Giant BND that only overlaps an intron of a non-coding transcript.
        // Exon 1 at 100-200, intron at 201-20M+300, exon 2 at 20M+301-20M+400.
        // BND from 1000 to 10_002_000 overlaps intron but not either exon.
        use std::sync::Arc;
        use vep_core::transcript::{Exon, Intron};

        let exons = vec![
            Exon {
                stable_id: Some("ENSE_NC_I1".into()),
                start: 100,
                end: 200,
                rank: 1,
                phase: -1,
                end_phase: -1,
            },
            Exon {
                stable_id: Some("ENSE_NC_I2".into()),
                start: 20_000_301,
                end: 20_000_400,
                rank: 2,
                phase: -1,
                end_phase: -1,
            },
        ];
        let introns = vec![Intron {
            start: 201,
            end: 20_000_300,
            rank: 1,
        }];
        let vefc = TranscriptVEFC {
            codon_table: 1,
            five_prime_utr: None,
            three_prime_utr: None,
            translateable_seq: None,
            peptide: None,
            introns: introns.clone(),
            sorted_exons: exons.clone(),
            mapper: Some(TranscriptMapper {
                start_phase: 0,
                cdna_coding_start: 0,
                cdna_coding_end: 0,
                exon_coord_mapper: ExonCoordMapper::new(vec![
                    MapperPair {
                        from_start: 1,
                        from_end: 101,
                        to_start: 100,
                        to_end: 200,
                        ori: 1,
                    },
                    MapperPair {
                        from_start: 102,
                        from_end: 201,
                        to_start: 20_000_301,
                        to_end: 20_000_400,
                        ori: 1,
                    },
                ]),
            }),
            protein_features: vec![],
            protein_function_predictions: None,
            seq_edits: vec![],
        };
        let tx = Transcript {
            stable_id: "ENST_NC_INTRON".into(),
            version: Some(1),
            db_id: None,
            gene_stable_id: "ENSG_NC_INTRON".into(),
            chr: "21".into(),
            start: 100,
            end: 20_000_400,
            strand: Strand::Forward,
            biotype: "lncRNA".into(),
            source: "Ensembl".into(),
            description: None,
            gene_symbol: Some("TESTNCI".into()),
            gene_symbol_source: None,
            hgnc_id: None,
            gene_phenotype: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            exons: exons.clone(),
            introns,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            translation: None,
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from([]),
            gencode_primary: false,
            attributes: vec![],
            vefc: Some(vefc),
            derived: Default::default(),
        };

        // BND: 1000 to 10_002_000 (>10M span). Overlaps intron but not exons.
        let v = make_giant_derived_sb(1000, 10_002_000);
        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(result.is_some(), "Giant BND should annotate");
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "Should have feature_truncation: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "Should have non_coding_transcript_variant (Perl parity): {:?}",
            tc.consequences
        );
        // `intron_variant` IS owed on the three-term intron-only shape
        // `feature_truncation, intron_variant, non_coding_transcript_variant`:
        // `--max_sv_size` decides whether an SV is annotated, never which terms it
        // gets.
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "an intron-only overlap owes intron_variant regardless of span size: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "Should not have non_coding_transcript_exon_variant: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_cds_predicates_fall_back_to_the_mapper_when_fields_are_absent() {
        // The SV CDS predicates must fall back to the mapper when
        // `coding_region_start` / `_end` are absent: Perl derives the bounds through
        // its own mapper, so a transcript whose cache entry carries its CDS only as
        // cDNA offsets still gets `coding_sequence_variant` and both UTR terms.
        //
        // Built by clearing the fields on a transcript that already carries a mapper, so
        // the two configurations differ in exactly one thing and the assertion cannot pass
        // for an unrelated reason.
        let with_fields = crate::test_helpers::make_test_transcript();
        assert!(
            with_fields.coding_region_start.is_some() && with_fields.vefc.is_some(),
            "fixture must carry both the fields and a mapper for this to isolate anything"
        );
        let (lo, hi) = (
            with_fields.coding_region_start.unwrap(),
            with_fields.coding_region_end.unwrap(),
        );
        assert!(
            overlaps_cds_exon(&with_fields, lo, lo),
            "sanity: the CDS start must read as a coding exon overlap"
        );

        let mut cleared = with_fields.clone();
        cleared.coding_region_start = None;
        cleared.coding_region_end = None;
        assert_eq!(
            crate::consequences::genomic_coding_bounds(&cleared),
            Some((lo, hi)),
            "the mapper must reproduce the bounds the cleared fields carried"
        );
        assert!(
            overlaps_cds_exon(&cleared, lo, lo),
            "with the fields cleared the predicate must still see the coding exon"
        );
    }

    #[test]
    fn test_small_bnd_non_coding_partial_exon_overlap_keeps_specific_terms() {
        // This BND extends beyond the transcript's 5' end but stops inside exon 3, so
        // it partially overlaps rather than engulfing, `complete_overlap_feature` is
        // false, and Perl's `non_coding_exon_variant` applies; `deletion::calculate`
        // derives the specific terms and nothing may strip them to the generic parent.
        //
        // All three positional terms are MODIFIER and the retained `feature_truncation` is
        // HIGH, so IMPACT is identical either way: only a term-level assertion can see
        // this, which is why it is asserted per term rather than through the impact field.
        let tx = make_non_coding_transcript(); // 25.000M - 25.006M, 3 exons
        let mut v = InputVariant::new(
            "21".into(),
            24_990_000,
            25_005_000, // inside exon 3 (25.004M - 25.006M), so the transcript survives
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = true;
        v.sv_end = Some(25_005_000);
        v.mate_id = Some("mate_nc_002".into());

        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("partial-overlap BND must annotate a non-coding transcript");
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "partial overlap truncates the feature: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences.contains(&Consequence::TranscriptAblation),
            "the transcript is not engulfed, so nothing is ablated: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "the specific exon term must survive; nothing may strip it to the generic \
             parent: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "Perl's two non-coding forms are mutually exclusive (non_coding_exon_variant \
             at VariationEffect.pm:502-519 against within_non_coding_gene at :495-500), \
             so the generic term must NOT accompany the exon term: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_small_bnd_non_coding_exon_and_intron_gives_the_three_term_form() {
        // Perl's three-term shape,
        // `feature_truncation,non_coding_transcript_exon_variant,intron_variant`: the
        // span extends beyond the 5' end, covers exon 1 and reaches into intron 1, so
        // both positional terms apply and the generic parent does not, Perl's two
        // non-coding forms being complements rather than a hierarchy.
        let tx = make_non_coding_transcript();
        let mut v = InputVariant::new(
            "21".into(),
            24_990_000,
            25_001_000, // covers exon 1 (25.000000-25.000299), ends inside intron 1
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = true;
        v.sv_end = Some(25_001_000);
        v.mate_id = Some("mate_nc_003".into());

        let tc = calculate_sv_consequences(&v, &tx, 5000, 5000)
            .expect("intron-reaching BND must annotate a non-coding transcript");
        for expected in [
            Consequence::FeatureTruncation,
            Consequence::NonCodingTranscriptExonVariant,
            Consequence::IntronVariant,
        ] {
            assert!(
                tc.consequences.contains(&expected),
                "the three-term form requires {expected:?}: {:?}",
                tc.consequences
            );
        }
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "the generic parent must not accompany the exon term: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_small_bnd_non_coding_containment_uses_generic_term() {
        // Small BND that engulfs a non-coding transcript, which is the ablation path.
        // For an engulfing variant the generic term is correct and the exon term is
        // wrong: Perl's `non_coding_exon_variant` returns 0 under
        // `complete_overlap_feature` (Utils/VariationEffect.pm:502-519 testing
        // `:169-178`).
        let tx = make_non_coding_transcript(); // 25M - 25.006M
        let mut v = InputVariant::new(
            "21".into(),
            24_990_000,
            25_010_000,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = true;
        v.sv_end = Some(25_010_000); // 20K span, below 10M
        v.mate_id = Some("mate_nc_001".into());

        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(result.is_some(), "Small BND should annotate non-coding tx");
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "Should have feature_truncation (ablation converted): {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences.contains(&Consequence::TranscriptAblation),
            "Should not have transcript_ablation (converted to truncation): {:?}",
            tc.consequences
        );
        // The generic term comes from the ablation branch's own push: this variant's
        // 20 kb span contains the whole 6 kb transcript, so `complete_overlap_feature`
        // holds and no specific term is applicable.
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "an engulfing BND must not carry the specific exon term: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "BND non-coding should have generic non_coding_transcript_variant: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_small_bnd_coding_intron_keeps_intron_variant() {
        // Perl keeps intron_variant for a coding BND-format SV alongside
        // feature_truncation.
        let tx = crate::test_helpers::make_test_transcript(); // coding, 25M - 25.006M
        let mut v = InputVariant::new(
            "21".into(),
            24_998_000,
            25_001_500,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.is_single_breakend = true;
        v.sv_end = Some(25_001_500); // Spans exon 1 + into intron 1, below 10M
        v.mate_id = Some("mate_coding_001".into());

        let result = calculate_sv_consequences(&v, &tx, 5000, 5000);
        assert!(result.is_some(), "Small BND should annotate coding tx");
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "Should have feature_truncation: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "Small BND should keep intron_variant (Perl keeps it for real-world SVs): {:?}",
            tc.consequences
        );
    }

    /// Build a non-coding **miRNA** transcript whose mature-miRNA subfeature is
    /// resolvable via the cDNA->genomic mapper.
    ///
    /// Layout (forward strand, chr21):
    ///   Exon 1: 25_000_000 - 25_000_100  (cDNA 1-101)
    ///   miRNA attribute "50-70" => cDNA 50-70 => genomic 25_000_049 - 25_000_069
    ///
    /// The transcript is intentionally single-exon and non-coding so that the only
    /// consequence Perl can reach is `mature_miRNA_variant` (tier 2) or the tier-3
    /// context term `non_coding_transcript_variant`.
    pub(crate) fn make_mature_mirna_transcript() -> Transcript {
        use std::sync::Arc;
        use vep_core::transcript::{
            Attribute, Exon, ExonCoordMapper, MapperPair, TranscriptMapper, TranscriptVEFC,
        };

        let exons = vec![Exon {
            stable_id: Some("ENSE_MIR_1".into()),
            start: 25_000_000,
            end: 25_000_100,
            rank: 1,
            phase: -1,
            end_phase: -1,
        }];

        let vefc = TranscriptVEFC {
            codon_table: 1,
            five_prime_utr: None,
            three_prime_utr: None,
            translateable_seq: None,
            peptide: None,
            introns: vec![],
            sorted_exons: exons.clone(),
            mapper: Some(TranscriptMapper {
                start_phase: 0,
                cdna_coding_start: 0,
                cdna_coding_end: 0,
                exon_coord_mapper: ExonCoordMapper::new(vec![MapperPair {
                    from_start: 1,
                    from_end: 101,
                    to_start: 25_000_000,
                    to_end: 25_000_100,
                    ori: 1,
                }]),
            }),
            protein_features: vec![],
            protein_function_predictions: None,
            seq_edits: vec![],
        };

        Transcript {
            stable_id: "ENST_MIR_1".into(),
            version: Some(1),
            db_id: None,
            gene_stable_id: "ENSG_MIR_1".into(),
            chr: "21".into(),
            start: 25_000_000,
            end: 25_000_100,
            strand: Strand::Forward,
            biotype: "miRNA".into(),
            source: "Ensembl".into(),
            description: None,
            gene_symbol: Some("MIRTEST".into()),
            gene_symbol_source: Some("HGNC".into()),
            hgnc_id: None,
            gene_phenotype: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            refseq: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            exons,
            introns: vec![],
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            translation_start: None,
            translation_end: None,
            translation: None,
            cdna_sequence: None,
            protein_sequence: None,
            flags: Arc::from([]),
            gencode_primary: false,
            attributes: vec![Attribute {
                code: "miRNA".into(),
                value: "50-70".into(),
            }],
            vefc: Some(vefc),
            derived: Default::default(),
        }
    }
}
