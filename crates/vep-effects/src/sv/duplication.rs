// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Structural duplication (`<DUP>`, `<DUP:TANDEM>`) consequence calculation.
//!
//! Mirrors Ensembl VEP's `StructuralVariationOverlapAllele` logic for duplications.
//! Consequence assignment follows the Ensembl VEP rules:
//!
//! - **transcript_amplification**: DUP fully contains the transcript (OverlapPC == 100%).
//! - **feature_elongation**: DUP lies entirely within the transcript and overlaps
//!   exonic (cDNA) sequence or lies wholly inside a frameshift intron; a DUP that
//!   extends beyond either transcript boundary does not get it.
//! - Regional sub-consequences are added when the DUP overlaps specific transcript
//!   regions (CDS, UTRs, introns, non-coding exons).
//! - **start_lost**: the DUP span covers a start-codon base with both ends in exons
//!   (`insertion::structural_start_lost`).
//! - **upstream_gene_variant / downstream_gene_variant**: DUP is entirely within the
//!   upstream or downstream region (no transcript body overlap).
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`).

use super::insertion::structural_start_lost;
use super::{
    frameshift_intron_span, is_mature_mirna_sv, overlaps_any_exon, overlaps_any_intron_trimmed,
    overlaps_cds_exon, overlaps_five_prime_utr, overlaps_polypyrimidine_tract,
    overlaps_three_prime_utr, FrameshiftIntronSpan,
};
use smallvec::SmallVec;
use std::sync::Arc;
use vep_core::consequence::{Consequence, ConsequenceList, Impact, TranscriptConsequence};
use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

/// Calculate consequences for a structural duplication against a transcript.
pub fn calculate(
    variant: &InputVariant,
    transcript: &Arc<Transcript>,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    let sv_start = variant.start;
    let sv_end = variant.sv_end.unwrap_or(variant.end);

    let sv_lo = sv_start.min(sv_end);
    let sv_hi = sv_start.max(sv_end);

    let tx_start = transcript.start;
    let tx_end = transcript.end;

    if sv_hi < tx_start.saturating_sub(upstream_distance.max(downstream_distance))
        || sv_lo > tx_end.saturating_add(upstream_distance.max(downstream_distance))
    {
        return None;
    }

    let overlaps_transcript = sv_hi >= tx_start && sv_lo <= tx_end;

    let (upstream_start, upstream_end, downstream_start, downstream_end) = match transcript.strand {
        Strand::Forward => {
            let up_start = tx_start.saturating_sub(upstream_distance);
            let up_end = tx_start.saturating_sub(1);
            let down_start = tx_end.saturating_add(1);
            let down_end = tx_end.saturating_add(downstream_distance);
            (up_start, up_end, down_start, down_end)
        }
        Strand::Reverse => {
            // For reverse strand, "upstream" is after the transcript end,
            // "downstream" is before the transcript start.
            let up_start = tx_end.saturating_add(1);
            let up_end = tx_end.saturating_add(upstream_distance);
            let down_start = tx_start.saturating_sub(downstream_distance);
            let down_end = tx_start.saturating_sub(1);
            (up_start, up_end, down_start, down_end)
        }
    };

    let overlaps_upstream =
        upstream_end >= upstream_start && sv_lo <= upstream_end && sv_hi >= upstream_start;
    let overlaps_downstream =
        downstream_end >= downstream_start && sv_lo <= downstream_end && sv_hi >= downstream_start;

    if !overlaps_transcript && !overlaps_upstream && !overlaps_downstream {
        return None;
    }

    let mut consequences: ConsequenceList = SmallVec::new();

    if overlaps_transcript {
        let fully_contains = sv_lo <= tx_start && sv_hi >= tx_end;

        if fully_contains {
            // CN2 inputs are normalized to Duplication upstream and take this path.
            consequences.push(Consequence::TranscriptAmplification);
        } else if is_mature_mirna_sv(transcript, sv_lo, sv_hi) {
            // Ensembl's tier-2 short-circuit outranks everything in the `else` below:
            // `@SORTED_OVERLAP_CONSEQUENCES = sort {$a->tier <=> $b->tier}`
            // (`BaseVariationFeatureOverlapAllele.pm:69`) orders predicates by tier
            // alone and `get_all_OverlapConsequences` (`:266,275-278`) stops at
            // `last if $assigned_tier && $oc->{tier} > $assigned_tier`.
            // `mature_miRNA_variant` is tier 2; `feature_elongation` and every
            // regional term are tier 3. Only tier-1 `transcript_amplification`,
            // in the arm above, can co-occur.
            consequences.push(Consequence::MatureMirnaVariant);
        } else {
            // Ensembl VEP `feature_elongation` predicate:
            //   within_cdna(@_) AND complete_within_feature(@_) AND copy_number_gain(@_)
            // so it fires when the DUP lies entirely inside the transcript and
            // overlaps exonic sequence, not when it extends beyond. `within_cdna`
            // also holds for a DUP wholly inside a frameshift intron, which takes
            // the terms of `frameshift_intron_span`.
            let complete_within = sv_lo >= tx_start && sv_hi <= tx_end;
            let frameshift_intron = frameshift_intron_span(transcript, sv_lo, sv_hi);
            let overlaps_cdna =
                overlaps_any_exon(transcript, sv_lo, sv_hi) || frameshift_intron.is_some();
            if complete_within && overlaps_cdna {
                consequences.push(Consequence::FeatureElongation);
            }

            add_regional_consequences(
                transcript,
                sv_lo,
                sv_hi,
                frameshift_intron,
                &mut consequences,
            );
        }

        // A row no predicate claims carries Ensembl VEP's default term:
        // `$cons = [$DEFAULT_OVERLAP_CONSEQUENCE] unless @$cons`
        // (`BaseVariationFeatureOverlapAllele.pm:285`), `intergenic_variant`.
        if consequences.is_empty() {
            consequences.push(Consequence::IntergenicVariant);
        }

        let impact = determine_impact(&consequences);
        Some(build_base_consequence(transcript, &consequences, impact))
    } else if overlaps_upstream && !overlaps_downstream {
        consequences.push(Consequence::UpstreamGeneVariant);
        let distance = compute_distance_upstream(transcript, sv_lo, sv_hi);
        let mut tc = build_base_consequence(
            transcript,
            &consequences,
            Consequence::UpstreamGeneVariant.impact(),
        );
        tc.distance = Some(distance);
        Some(tc)
    } else if overlaps_downstream && !overlaps_upstream {
        consequences.push(Consequence::DownstreamGeneVariant);
        let distance = compute_distance_downstream(transcript, sv_lo, sv_hi);
        let mut tc = build_base_consequence(
            transcript,
            &consequences,
            Consequence::DownstreamGeneVariant.impact(),
        );
        tc.distance = Some(distance);
        Some(tc)
    } else {
        // Both flanks but not the body cannot occur (both flanks imply the body);
        // handled rather than panicking.
        None
    }
}

/// Add regional sub-consequences based on which transcript regions the DUP overlaps;
/// `frameshift_intron` is the reading of a DUP wholly inside a frameshift intron.
fn add_regional_consequences(
    transcript: &Transcript,
    sv_lo: u64,
    sv_hi: u64,
    frameshift_intron: Option<FrameshiftIntronSpan>,
    consequences: &mut ConsequenceList,
) {
    let is_protein_coding = transcript.has_cds();

    if is_protein_coding {
        if structural_start_lost(transcript, sv_lo, sv_hi) {
            push_unique(consequences, Consequence::StartLost);
        }
        if overlaps_cds_exon(transcript, sv_lo, sv_hi)
            || frameshift_intron.is_some_and(|span| span.coding)
        {
            push_unique(consequences, Consequence::CodingSequenceVariant);
        }
        if overlaps_five_prime_utr(transcript, sv_lo, sv_hi)
            || frameshift_intron.is_some_and(|span| span.five_prime_utr)
        {
            push_unique(consequences, Consequence::FivePrimeUtrVariant);
        }
        if overlaps_three_prime_utr(transcript, sv_lo, sv_hi)
            || frameshift_intron.is_some_and(|span| span.three_prime_utr)
        {
            push_unique(consequences, Consequence::ThreePrimeUtrVariant);
        }
    } else {
        // No mature-miRNA arm here: the sole caller returns on `is_mature_mirna_sv`
        // before reaching this function, and a second copy of that test could drift.
        if overlaps_any_exon(transcript, sv_lo, sv_hi) {
            push_unique(consequences, Consequence::NonCodingTranscriptExonVariant);
        }
    }

    // Intron overlap on Ensembl's trimmed boundaries: `_intron_effects` sets
    // `intronic` with `overlap($r_start, $r_end, $intron_start + 2, $intron_end - 2)`
    // (`BaseTranscriptVariationAllele.pm:143-150`), so an SV touching only the four
    // invariant donor/acceptor bases gets no `intron_variant`.
    if overlaps_any_intron_trimmed(transcript, sv_lo, sv_hi) {
        push_unique(consequences, Consequence::IntronVariant);
    }
    if overlaps_polypyrimidine_tract(transcript, sv_lo, sv_hi) {
        push_unique(consequences, Consequence::SplicePolypyrimidineTractVariant);
    }

    // Ensembl VEP: for non-coding transcripts, non_coding_transcript_variant is added
    // when there is NO exon overlap (within_non_coding_gene requires not non_coding_exon_variant
    // and not mature_miRNA).
    if !is_protein_coding
        && !transcript.facts().is_nmd
        && !consequences.contains(&Consequence::NonCodingTranscriptExonVariant)
        && !consequences.contains(&Consequence::MatureMirnaVariant)
        && !consequences.contains(&Consequence::NonCodingTranscriptVariant)
    {
        push_unique(consequences, Consequence::NonCodingTranscriptVariant);
    }

    // NMD transcript context: always added alongside other consequences.
    if transcript.facts().is_nmd && !consequences.contains(&Consequence::NmdTranscriptVariant) {
        push_unique(consequences, Consequence::NmdTranscriptVariant);
    }
}

/// Determine the most severe impact from a set of consequences.
fn determine_impact(consequences: &[Consequence]) -> Impact {
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

/// Compute distance from the DUP to the transcript for upstream variants.
fn compute_distance_upstream(transcript: &Transcript, sv_lo: u64, sv_hi: u64) -> u64 {
    match transcript.strand {
        Strand::Forward => {
            // Upstream is before transcript start on forward strand.
            // Distance is from the closest edge of the DUP to transcript start.
            transcript.start.saturating_sub(sv_hi)
        }
        Strand::Reverse => {
            // Upstream is after transcript end on reverse strand.
            sv_lo.saturating_sub(transcript.end)
        }
    }
}

/// Compute distance from the DUP to the transcript for downstream variants.
fn compute_distance_downstream(transcript: &Transcript, sv_lo: u64, sv_hi: u64) -> u64 {
    match transcript.strand {
        Strand::Forward => {
            // Downstream is after transcript end on forward strand.
            sv_lo.saturating_sub(transcript.end)
        }
        Strand::Reverse => {
            // Downstream is before transcript start on reverse strand.
            transcript.start.saturating_sub(sv_hi)
        }
    }
}

/// Build a base `TranscriptConsequence` with transcript metadata populated.
fn build_base_consequence(
    transcript: &Arc<Transcript>,
    consequences: &[Consequence],
    impact: Impact,
) -> TranscriptConsequence {
    let mut sorted_csqs: ConsequenceList = SmallVec::from_slice(consequences);
    sorted_csqs.sort_by_key(|c| c.rank());

    TranscriptConsequence {
        consequences: sorted_csqs,
        impact,
        ..TranscriptConsequence::for_transcript(transcript)
    }
}

/// Push a consequence if not already present.
fn push_unique(consequences: &mut ConsequenceList, consequence: Consequence) {
    if !consequences.contains(&consequence) {
        consequences.push(consequence);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::{make_test_transcript, make_test_transcript_with_flags};
    use vep_core::consequence::FeatureType;
    use vep_core::variant::VariantClass;

    /// [`super::calculate`] on a transcript the test owns: the row shares a copy of
    /// it, as the index shares its handle.
    fn calculate(
        variant: &InputVariant,
        transcript: &Transcript,
        upstream_distance: u64,
        downstream_distance: u64,
    ) -> Option<TranscriptConsequence> {
        super::calculate(
            variant,
            &Arc::new(transcript.clone()),
            upstream_distance,
            downstream_distance,
        )
    }

    /// Helper: build a structural duplication variant.
    fn make_dup(start: u64, sv_end: u64, class: VariantClass) -> InputVariant {
        let mut v = InputVariant::new("21".into(), start, sv_end, b"N".to_vec(), b"-".to_vec());
        v.variant_class = class;
        v.is_structural = true;
        v.sv_end = Some(sv_end);
        v
    }

    // Test transcript layout (from make_test_transcript):
    //   TX:     25_000_000 - 25_006_000  (forward strand)
    //   Exon 1: 25_000_000 - 25_000_299
    //   Intron 1: 25_000_300 - 25_001_999
    //   Exon 2: 25_002_000 - 25_002_299
    //   Intron 2: 25_002_300 - 25_003_999
    //   Exon 3: 25_004_000 - 25_006_000
    //   CDS:    25_000_050 - 25_004_299  (coding_region_start/end)
    //   5' UTR: 25_000_000 - 25_000_049
    //   3' UTR: 25_004_300 - 25_006_000

    #[test]
    fn test_dup_transcript_amplification() {
        let tx = make_test_transcript();
        let v = make_dup(24_990_000, 25_010_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::TranscriptAmplification));
        assert_eq!(tc.impact, Impact::HIGH);
        assert_eq!(
            tc.feature_overlap(v.start, v.end)
                .map(|(_, pc)| format!("{pc:.2}")),
            Some("100.00".to_string())
        );
    }

    #[test]
    fn test_dup_transcript_amplification_exact_boundaries() {
        let tx = make_test_transcript();
        let v = make_dup(25_000_000, 25_006_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::TranscriptAmplification));
    }

    #[test]
    fn test_dup_cn2_full_containment_transcript_amplification() {
        let tx = make_test_transcript();
        let mut v = InputVariant::new(
            "21".into(),
            24_990_000,
            25_010_000,
            b"N".to_vec(),
            b"<CN2>".to_vec(),
        );
        v.variant_class = VariantClass::Duplication;
        v.is_structural = true;
        v.sv_end = Some(25_010_000);

        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::TranscriptAmplification));
    }

    #[test]
    fn test_dup_extends_left_no_feature_elongation() {
        let tx = make_test_transcript();
        // Ensembl VEP: feature_elongation requires complete_within_feature, and this
        // DUP extends beyond the transcript.
        let v = make_dup(24_999_000, 25_003_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            !tc.consequences.contains(&Consequence::FeatureElongation),
            "DUP extending beyond transcript should not get feature_elongation"
        );
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_dup_extends_right_no_feature_elongation() {
        let tx = make_test_transcript();
        // Extends past the transcript end, so not complete_within_feature.
        let v = make_dup(25_003_000, 25_010_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            !tc.consequences.contains(&Consequence::FeatureElongation),
            "DUP extending beyond transcript should not get feature_elongation"
        );
    }

    #[test]
    fn test_dup_internal_cds_overlap_has_feature_elongation() {
        let tx = make_test_transcript();
        // DUP entirely within the transcript, overlapping CDS and introns.
        // Ensembl VEP: feature_elongation fires when complete_within_feature
        // and within_cdna (overlaps exonic sequence).
        let v = make_dup(25_000_100, 25_002_100, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        // feature_elongation: the DUP lies inside the transcript and overlaps exons.
        assert!(
            tc.consequences.contains(&Consequence::FeatureElongation),
            "Internal DUP overlapping exons should get feature_elongation, got: {:?}",
            tc.consequences
        );
        assert!(!tc
            .consequences
            .contains(&Consequence::TranscriptAmplification));
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_dup_internal_intron_only() {
        let tx = make_test_transcript();
        let v = make_dup(25_000_500, 25_001_500, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(!tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert_eq!(tc.impact, Impact::MODIFIER);
    }

    #[test]
    fn test_dup_five_prime_utr_overlap() {
        let tx = make_test_transcript();
        let v = make_dup(25_000_010, 25_000_040, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::FivePrimeUtrVariant));
    }

    #[test]
    fn test_dup_three_prime_utr_overlap() {
        let tx = make_test_transcript();
        let v = make_dup(25_004_500, 25_005_500, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::ThreePrimeUtrVariant));
    }

    /// A reverse-strand transcript whose coding region reaches its first genomic
    /// base has no 3' UTR, and a span starting on that base gets no UTR term:
    /// Perl's `_before_coding` is `overlap(start, end, transcript start,
    /// coding_region_start - 1)`, an inverted window that holds only for a span
    /// starting before the transcript. A span one base further out does hold it.
    #[test]
    fn test_dup_reverse_strand_cds_start_nf_has_no_three_prime_utr() {
        let mut tx = make_test_transcript_with_flags(&["cds_start_NF", "cds_end_NF"]);
        tx.strand = Strand::Reverse;
        tx.coding_region_start = Some(tx.start);

        let v = make_dup(tx.start, 25_001_000, VariantClass::Duplication);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();

        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(!tc.consequences.contains(&Consequence::ThreePrimeUtrVariant));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));

        let v = make_dup(tx.start - 1, 25_001_000, VariantClass::Duplication);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::ThreePrimeUtrVariant));
    }

    #[test]
    fn test_dup_upstream_gene_variant() {
        let tx = make_test_transcript();
        let v = make_dup(24_996_000, 24_998_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::UpstreamGeneVariant));
        assert!(tc.distance.is_some());
        assert_eq!(tc.distance.unwrap(), 25_000_000 - 24_998_000);
    }

    #[test]
    fn test_dup_downstream_gene_variant() {
        let tx = make_test_transcript();
        let v = make_dup(25_007_000, 25_009_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert!(tc.distance.is_some());
        assert_eq!(tc.distance.unwrap(), 25_007_000 - 25_006_000);
    }

    #[test]
    fn test_dup_no_overlap() {
        let tx = make_test_transcript();
        let v = make_dup(26_000_000, 26_100_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_none());
    }

    #[test]
    fn test_dup_tandem_duplication_class() {
        let tx = make_test_transcript();
        // TandemDuplication uses the same logic as Duplication.
        let v = make_dup(24_990_000, 25_010_000, VariantClass::TandemDuplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::TranscriptAmplification));
    }

    #[test]
    fn test_dup_overlap_bp_and_pct_values() {
        let tx = make_test_transcript();
        // DUP over the first 3001 bases of the transcript's 6001-base span.
        let v = make_dup(25_000_000, 25_003_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        let (obp, opct) = tc.feature_overlap(v.start, v.end).unwrap();
        // Overlap: 25_000_000 to 25_003_000 = 3001 bp
        assert_eq!(obp, 3001);
        // TX length = 25_006_000 - 25_000_000 + 1 = 6001
        let expected_pct = 100.0 * 3001.0 / 6001.0;
        assert!((opct - expected_pct).abs() < 0.01);
    }

    #[test]
    fn test_dup_multiple_regions() {
        let tx = make_test_transcript();
        // Contained within the transcript and overlapping exons, so Ensembl VEP
        // fires feature_elongation (complete_within_feature and within_cdna).
        let v = make_dup(25_000_010, 25_005_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(tc.consequences.contains(&Consequence::FivePrimeUtrVariant));
        assert!(tc.consequences.contains(&Consequence::ThreePrimeUtrVariant));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        // A DUP inside the transcript overlapping exons has feature_elongation.
        assert!(
            tc.consequences.contains(&Consequence::FeatureElongation),
            "Internal DUP overlapping exons should get feature_elongation"
        );
        // Both ends exonic (exon 1 and exon 3) and the span covers the start codon.
        assert!(tc.consequences.contains(&Consequence::StartLost));
        assert!(!tc.consequences.contains(&Consequence::StartRetainedVariant));
    }

    /// A DUP inside exon 1 over the start codon (25_000_050-25_000_052) is exactly
    /// `feature_elongation, start_lost, coding_sequence_variant, 5_prime_UTR_variant`,
    /// Perl's structural `start_lost` arm without its co-fired `start_retained_variant`.
    #[test]
    fn test_dup_over_start_codon_within_exon_is_start_lost() {
        let tx = make_test_transcript();
        let v = make_dup(25_000_040, 25_000_121, VariantClass::Duplication);
        let tc = calculate(&v, &tx, 5000, 5000).expect("DUP must annotate the transcript");
        assert_eq!(
            tc.consequences.to_vec(),
            vec![
                Consequence::StartLost,
                Consequence::FeatureElongation,
                Consequence::CodingSequenceVariant,
                Consequence::FivePrimeUtrVariant,
            ],
            "got: {:?}",
            tc.consequences
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    /// The same DUP ending in intron 1 has no `cdna_end`, so no `start_lost`.
    #[test]
    fn test_dup_over_start_codon_ending_in_intron_has_no_start_lost() {
        let tx = make_test_transcript();
        let v = make_dup(25_000_040, 25_001_000, VariantClass::Duplication);
        let tc = calculate(&v, &tx, 5000, 5000).expect("DUP must annotate the transcript");
        assert!(!tc.consequences.contains(&Consequence::StartLost));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_dup_transcript_metadata() {
        let tx = make_test_transcript();
        let v = make_dup(24_990_000, 25_010_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(tc.transcript_id(), "ENST00000000001");
        assert_eq!(tc.gene_id(), "ENSG00000000001");
        assert_eq!(tc.gene_symbol(), Some("TEST1"));
        assert!(tc.canonical);
        assert_eq!(tc.strand, 1);
        assert_eq!(tc.feature_type, FeatureType::Transcript);
        assert_eq!(tc.biotype(), Some("protein_coding"));
    }

    #[test]
    fn test_dup_no_sv_end_fallback() {
        let tx = make_test_transcript();
        // Without sv_end, variant.end is the fallback.
        let mut v = InputVariant::new(
            "21".into(),
            25_000_100,
            25_002_100,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::Duplication;
        v.is_structural = true;
        v.sv_end = None;
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
    }

    #[test]
    fn test_dup_consequence_order_by_severity() {
        let tx = make_test_transcript();
        let v = make_dup(24_999_000, 25_003_000, VariantClass::Duplication);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        // Verify consequences are sorted by rank (lower rank = more severe = first).
        for window in tc.consequences.windows(2) {
            assert!(
                window[0].rank() <= window[1].rank(),
                "Consequences not sorted by rank: {:?} (rank {}) before {:?} (rank {})",
                window[0],
                window[0].rank(),
                window[1],
                window[1].rank()
            );
        }
    }

    /// A duplication wholly inside a frameshift intron. Ensembl VEP 116.2 writes
    /// `feature_elongation,coding_sequence_variant` for `<DUP>` and `<DUP:TANDEM>`
    /// over GRCh37 `21:38600094-38600097` on ENST00000288304 (a 4-base intron inside
    /// its coding region), `feature_elongation,non_coding_transcript_variant` for
    /// GRCh38 `21:6358575-6358576` on the lncRNA ENST00000624847, and
    /// `intron_variant` in an ordinary intron.
    #[test]
    fn a_duplication_inside_a_frameshift_intron_is_an_elongation() {
        let mut tx = crate::sv::deletion::frameshift_intron_test_transcript();
        let terms = |tx: &Transcript, class: VariantClass, start: u64, end: u64| {
            let tc = calculate(&make_dup(start, end, class), tx, 5000, 5000).unwrap();
            let mut terms: Vec<String> = tc.consequences.iter().map(|c| c.to_string()).collect();
            terms.sort();
            terms
        };
        for class in [VariantClass::Duplication, VariantClass::TandemDuplication] {
            // The 9-base intron of the coding region, whole and in part.
            for (start, end) in [(25_002_101, 25_002_109), (25_002_103, 25_002_107)] {
                assert_eq!(
                    terms(&tx, class, start, end),
                    ["coding_sequence_variant", "feature_elongation"],
                    "{class:?} {start}-{end}"
                );
            }
            // The introns of the two UTRs.
            assert_eq!(
                terms(&tx, class, 25_000_021, 25_000_025),
                ["5_prime_UTR_variant", "feature_elongation"]
            );
            assert_eq!(
                terms(&tx, class, 25_005_001, 25_005_013),
                ["3_prime_UTR_variant", "feature_elongation"]
            );
            // Nine bases inside an ordinary intron keep the intronic reading.
            assert_eq!(
                terms(&tx, class, 25_001_000, 25_001_008),
                ["intron_variant"]
            );
        }
        // Without a translation the context term stands beside the elongation.
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        if let Some(mapper) = tx.vefc.as_mut().and_then(|v| v.mapper.as_mut()) {
            mapper.cdna_coding_start = 0;
            mapper.cdna_coding_end = 0;
        }
        assert_eq!(
            terms(&tx, VariantClass::Duplication, 25_002_101, 25_002_109),
            ["feature_elongation", "non_coding_transcript_variant"]
        );
    }

    /// A duplication on only the first or last two bases of an intron, which no
    /// predicate claims: Ensembl VEP 116.2 writes `intergenic_variant` for `<DUP>`,
    /// `<DUP:TANDEM>` and `<CN2>` over GRCh37 `21:38600108-38600109` and
    /// `21:38600522-38600523` on ENST00000288304 (the ends of its 416-base intron),
    /// and `intron_variant` once the span reaches the intron's third base.
    #[test]
    fn a_duplication_on_the_splice_bases_of_an_intron_takes_the_default_term() {
        let tx = make_test_transcript();
        let terms = |class: VariantClass, start: u64, end: u64| {
            let tc = calculate(&make_dup(start, end, class), &tx, 5000, 5000).unwrap();
            tc.consequences
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
        };
        for class in [VariantClass::Duplication, VariantClass::TandemDuplication] {
            // Intron 1 is 25_000_300-25_001_999.
            for (start, end) in [
                (25_000_300, 25_000_301),
                (25_000_300, 25_000_300),
                (25_000_301, 25_000_301),
                (25_001_998, 25_001_999),
            ] {
                assert_eq!(
                    terms(class, start, end),
                    ["intergenic_variant"],
                    "{class:?} {start}-{end}"
                );
            }
            assert_eq!(terms(class, 25_000_301, 25_000_302), ["intron_variant"]);
        }
    }
}
