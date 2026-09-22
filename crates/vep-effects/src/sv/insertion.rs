// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Structural insertion (`<INS>`, `<INS:ME:*>`, `<CNV:TR>`) consequence calculation.
//!
//! Mirrors Perl VEP's `StructuralVariationOverlapAllele` logic for insertions.
//!
//! ## Coordinate model
//!
//! Structural insertions come in two flavours:
//!
//! 1. **Ranged insertions** (`<INS>` with SVLEN): The VCF parser sets
//!    `variant.start` = POS and `variant.end` = POS + abs(SVLEN). Perl VEP
//!    treats the inserted material as occupying this range for overlap
//!    purposes, even though no genomic bases are deleted.
//!
//! 2. **Point insertions** (mobile elements, tandem repeats without SVLEN):
//!    `variant.start == variant.end` (single position). OverlapBP = 1.
//!
//! ## Consequence assignment
//!
//! For each transcript region that the SV range overlaps:
//! - `feature_elongation`: the insertion lies entirely inside the transcript and
//!   overlaps exonic (cDNA) sequence, coding or non-coding (Perl:
//!   `within_cdna AND complete_within_feature AND insertion`); HIGH impact,
//!   because the transcript gets longer.
//! - `coding_sequence_variant`: Overlap with any CDS exon.
//! - `5_prime_UTR_variant` / `3_prime_UTR_variant`: Overlap with UTR regions.
//! - `intron_variant`: Overlap with an intron.
//! - `non_coding_transcript_exon_variant`: Exon overlap in a non-coding transcript.
//! - `coding_transcript_variant` / `non_coding_transcript_variant`: Context term
//!   when the overlap is entirely intronic or the transcript is fully contained.
//! - `upstream_gene_variant` / `downstream_gene_variant`: SV is near but does
//!   not overlap the transcript body.
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`).

use super::{
    is_mature_mirna_sv, overlaps_any_exon, overlaps_any_intron_trimmed, overlaps_cds_exon,
    overlaps_five_prime_utr, overlaps_three_prime_utr,
};
use smallvec::{smallvec, SmallVec};
use vep_core::consequence::{
    Consequence, ConsequenceList, FeatureType, Impact, TranscriptConsequence,
};
use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;
use vep_core::variant::{InputVariant, VariantClass};

/// Calculate consequences for a structural insertion against a transcript.
///
/// Returns `None` when the SV region has no overlap with the transcript's
/// extended region (transcript body + upstream/downstream distance).
pub fn calculate(
    variant: &InputVariant,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    let sv_start = variant.start;
    let sv_end = variant.sv_end.unwrap_or(variant.end);

    let tx_start = transcript.start;
    let tx_end = transcript.end;

    if sv_end < tx_start.saturating_sub(upstream_distance.max(downstream_distance))
        || sv_start > tx_end + upstream_distance.max(downstream_distance)
    {
        return None;
    }

    let (up_boundary, down_boundary) = match transcript.strand {
        Strand::Forward => (
            tx_start.saturating_sub(upstream_distance),
            tx_end.saturating_add(downstream_distance),
        ),
        Strand::Reverse => (
            tx_end.saturating_add(upstream_distance),
            tx_start.saturating_sub(downstream_distance),
        ),
    };

    let mut tc = TranscriptConsequence {
        transcript_id: transcript.stable_id.clone(),
        feature_start: transcript.start,
        feature_end: transcript.end,
        gene_id: transcript.gene_stable_id.clone(),
        gene_symbol: transcript.gene_symbol.clone(),
        gene_symbol_source: transcript.gene_symbol_source.clone(),
        hgnc_id: transcript.hgnc_id.clone(),
        biotype: Some(transcript.biotype.clone()),
        canonical: transcript.canonical,
        strand: match transcript.strand {
            Strand::Forward => 1,
            Strand::Reverse => -1,
        },
        feature_type: FeatureType::Transcript,
        flags: transcript.flags.clone(),
        tsl: transcript.tsl,
        mane_select: transcript.mane_select.clone(),
        mane_plus_clinical: transcript.mane_plus_clinical.clone(),
        appris: transcript.appris.clone(),
        ccds: transcript.ccds.clone(),
        swissprot: transcript.swissprot.clone(),
        trembl: transcript.trembl.clone(),
        refseq: transcript.refseq.clone(),
        ..TranscriptConsequence::default()
    };

    let overlaps_transcript = sv_end >= tx_start && sv_start <= tx_end;

    if !overlaps_transcript {
        let csq =
            classify_upstream_downstream(sv_start, sv_end, transcript, up_boundary, down_boundary);
        if let Some((consequence, distance)) = csq {
            tc.consequences = smallvec![consequence];
            tc.impact = consequence.impact();
            tc.distance = Some(distance);
            return Some(tc);
        }
        return None;
    }

    // Perl VEP `coding_transcript_variant` / `non_coding_transcript_variant`:
    //   complete_overlap_feature and within_coding_gene (or non-coding equivalent)
    // A ranged insertion that fully contains the transcript gets the context term
    // as the sole consequence, like inversions and neutral CNVs. Exception:
    // TandemRepeat (<CNV:TR>) is a copy-number gain that Perl treats like a
    // duplication, so full containment yields `transcript_amplification`.
    let fully_contains = sv_start <= tx_start && sv_end >= tx_end;
    if fully_contains {
        if variant.variant_class == VariantClass::TandemRepeat {
            let consequence = Consequence::TranscriptAmplification;
            tc.consequences = smallvec![consequence];
            tc.impact = consequence.impact();
            return Some(tc);
        }

        let consequence = if transcript.has_cds() {
            Consequence::CodingTranscriptVariant
        } else if transcript.is_nmd_transcript() {
            Consequence::NmdTranscriptVariant
        } else {
            Consequence::NonCodingTranscriptVariant
        };
        tc.consequences = smallvec![consequence];
        tc.impact = consequence.impact();
        return Some(tc);
    }

    let mut consequences: ConsequenceList = SmallVec::new();

    let hits_exon = overlaps_any_exon(transcript, sv_start, sv_end);

    // Perl VEP `within_intron` reads `_intron_effects`, which trims 2bp from each
    // intron end before setting the "intronic" flag, attributing the invariant
    // GT/AG dinucleotides to the splice terms. The trim is unconditional in
    // `BaseTranscriptVariationAllele.pm:143-150`: one `overlap($r_start, $r_end,
    // $intron_start + 2, $intron_end - 2)` with no variant-class or span case.
    let hits_intron = overlaps_any_intron_trimmed(transcript, sv_start, sv_end);

    // Once the trim excludes a point mobile-element insertion's position, Ensembl's
    // predicates all return false and the fallthrough yields `intergenic_variant`
    // rather than a transcript context term.
    let is_point_me =
        sv_start == sv_end && variant.variant_class == VariantClass::MobileElementInsertion;

    let is_protein_coding = transcript.has_cds();

    // Perl VEP `feature_elongation` predicate for insertions/gains:
    //   within_cdna(@_) AND complete_within_feature(@_) AND insertion(@_)
    // complete_within_feature = SV is entirely inside the transcript
    let complete_within = sv_start >= tx_start && sv_end <= tx_end;

    if is_protein_coding {
        let hits_cds = overlaps_cds_exon(transcript, sv_start, sv_end);
        let hits_5utr = overlaps_five_prime_utr(transcript, sv_start, sv_end);
        let hits_3utr = overlaps_three_prime_utr(transcript, sv_start, sv_end);

        if complete_within && hits_exon {
            consequences.push(Consequence::FeatureElongation);
        }

        if hits_cds {
            consequences.push(Consequence::CodingSequenceVariant);
        }
        if hits_5utr {
            consequences.push(Consequence::FivePrimeUtrVariant);
        }
        if hits_3utr {
            consequences.push(Consequence::ThreePrimeUtrVariant);
        }
        if hits_intron {
            consequences.push(Consequence::IntronVariant);

            if sv_start == sv_end {
                add_point_insertion_splice_consequences(&mut consequences, transcript, sv_start);
            } else if variant.variant_class == VariantClass::TandemRepeat {
                // Only a TandemRepeat gets the polypyrimidine term among ranged
                // insertions in Perl VEP, and only by endpoint: the full span would
                // over-call on large SVs.
                add_ranged_endpoint_splice_consequences(
                    &mut consequences,
                    transcript,
                    sv_start,
                    sv_end,
                );
            }
        }

        // Intron-only overlap gets the context term instead of feature_elongation.
        if consequences.is_empty() && hits_intron {
            consequences.push(Consequence::IntronVariant);
        }
    } else {
        // feature_elongation only when INS is within the transcript and
        // overlaps exonic (cDNA) sequence (Perl: within_cdna and complete_within_feature).
        if hits_exon && complete_within {
            consequences.push(Consequence::FeatureElongation);
        }

        if hits_exon {
            if is_mature_mirna_sv(transcript, sv_start, sv_end) {
                consequences.push(Consequence::MatureMirnaVariant);
            } else {
                consequences.push(Consequence::NonCodingTranscriptExonVariant);
            }
        }
        if hits_intron {
            consequences.push(Consequence::IntronVariant);

            if sv_start == sv_end {
                add_point_insertion_splice_consequences(&mut consequences, transcript, sv_start);
            } else if variant.variant_class == VariantClass::TandemRepeat {
                add_ranged_endpoint_splice_consequences(
                    &mut consequences,
                    transcript,
                    sv_start,
                    sv_end,
                );
            }
        }
    }

    // Perl VEP `within_non_coding_gene` predicate:
    //   within_transcript AND NOT translation AND NOT mature_miRNA AND NOT non_coding_exon_variant
    // so an intron-only overlap on a non-coding transcript produces both
    // intron_variant and non_coding_transcript_variant. A point ME insertion at
    // a trimmed splice boundary gets neither: Perl falls through to
    // `intergenic_variant`.
    let me_at_splice_boundary = is_point_me && consequences.is_empty();
    if !is_protein_coding
        && !transcript.is_nmd_transcript()
        && !me_at_splice_boundary
        && !consequences.contains(&Consequence::NonCodingTranscriptExonVariant)
        && !consequences.contains(&Consequence::MatureMirnaVariant)
    {
        consequences.push(Consequence::NonCodingTranscriptVariant);
    }

    // A point ME insertion in the first or last 2bp of an intron (the splice
    // dinucleotides the trimmed check excludes) matches no Perl predicate, so
    // Perl VEP produces `intergenic_variant`.
    if consequences.is_empty() {
        if is_point_me {
            consequences.push(Consequence::IntergenicVariant);
        } else if is_protein_coding {
            consequences.push(Consequence::CodingTranscriptVariant);
        } else {
            consequences.push(Consequence::NonCodingTranscriptVariant);
        }
    }

    if transcript.is_nmd_transcript() && !consequences.contains(&Consequence::NmdTranscriptVariant)
    {
        consequences.push(Consequence::NmdTranscriptVariant);
    }

    consequences.sort_by_key(|c| c.rank());
    consequences.dedup();

    let impact = consequences
        .first()
        .map(|c| c.impact())
        .unwrap_or(Impact::MODIFIER);

    tc.consequences = consequences;
    tc.impact = impact;
    Some(tc)
}

/// Add splice_polypyrimidine_tract_variant for a point insertion within an intron.
///
/// This is the SV equivalent of the small-variant polypyrimidine tract check.
/// Only fires for point insertions (sv_start == sv_end) when the insertion point
/// falls 2-16 bases from the intron acceptor (strand-aware).
fn add_point_insertion_splice_consequences(
    consequences: &mut ConsequenceList,
    transcript: &Transcript,
    insertion_pos: u64,
) {
    let exons = transcript
        .vefc
        .as_ref()
        .map(|v| v.sorted_exons.as_slice())
        .unwrap_or(transcript.exons.as_slice());

    for pair in exons.windows(2) {
        let intron_start = pair[0].end + 1;
        let intron_end = pair[1].start.saturating_sub(1);
        if intron_end < intron_start {
            continue;
        }

        if insertion_pos < intron_start || insertion_pos > intron_end {
            continue;
        }

        let dist_to_genomic_start = insertion_pos - intron_start;
        let dist_to_genomic_end = intron_end - insertion_pos;

        // Forward strand: donor at intron start, acceptor at intron end.
        // Reverse strand: donor at intron end, acceptor at intron start.
        let dist_to_acceptor = match transcript.strand {
            Strand::Forward => dist_to_genomic_end,
            Strand::Reverse => dist_to_genomic_start,
        };

        // Polypyrimidine tract window: positions 2-16 from acceptor.
        if (2..=16).contains(&dist_to_acceptor)
            && !consequences.contains(&Consequence::SplicePolypyrimidineTractVariant)
        {
            consequences.push(Consequence::SplicePolypyrimidineTractVariant);
        }

        break;
    }
}

/// Add splice_polypyrimidine_tract_variant for a ranged insertion when an
/// endpoint falls in the polypyrimidine tract.
///
/// For ranged SVs (sv_start != sv_end), check each endpoint individually
/// against every intron. Only fires when an endpoint (not the full span)
/// falls 2-16bp from the intron acceptor (strand-aware); checking the full span
/// over-calls on large SVs.
fn add_ranged_endpoint_splice_consequences(
    consequences: &mut ConsequenceList,
    transcript: &Transcript,
    sv_start: u64,
    sv_end: u64,
) {
    let exons = transcript
        .vefc
        .as_ref()
        .map(|v| v.sorted_exons.as_slice())
        .unwrap_or(transcript.exons.as_slice());

    for pair in exons.windows(2) {
        let intron_start = pair[0].end + 1;
        let intron_end = pair[1].start.saturating_sub(1);
        if intron_end < intron_start {
            continue;
        }

        for &endpoint in &[sv_start, sv_end] {
            if endpoint < intron_start || endpoint > intron_end {
                continue;
            }

            let dist_to_genomic_start = endpoint - intron_start;
            let dist_to_genomic_end = intron_end - endpoint;

            let dist_to_acceptor = match transcript.strand {
                Strand::Forward => dist_to_genomic_end,
                Strand::Reverse => dist_to_genomic_start,
            };

            if (2..=16).contains(&dist_to_acceptor)
                && !consequences.contains(&Consequence::SplicePolypyrimidineTractVariant)
            {
                consequences.push(Consequence::SplicePolypyrimidineTractVariant);
                return;
            }
        }
    }
}

/// Classify whether the SV is upstream or downstream of the transcript.
///
/// Returns the consequence and the distance in bases.
fn classify_upstream_downstream(
    sv_start: u64,
    sv_end: u64,
    transcript: &Transcript,
    _up_boundary: u64,
    _down_boundary: u64,
) -> Option<(Consequence, u64)> {
    let tx_start = transcript.start;
    let tx_end = transcript.end;

    match transcript.strand {
        Strand::Forward => {
            if sv_end < tx_start {
                // SV is entirely before the transcript on forward strand = upstream.
                let distance = tx_start - sv_end;
                Some((Consequence::UpstreamGeneVariant, distance))
            } else if sv_start > tx_end {
                // SV is entirely after the transcript on forward strand = downstream.
                let distance = sv_start - tx_end;
                Some((Consequence::DownstreamGeneVariant, distance))
            } else {
                None
            }
        }
        Strand::Reverse => {
            if sv_start > tx_end {
                // SV is after the transcript on reverse strand = upstream.
                let distance = sv_start - tx_end;
                Some((Consequence::UpstreamGeneVariant, distance))
            } else if sv_end < tx_start {
                // SV is before the transcript on reverse strand = downstream.
                let distance = tx_start - sv_end;
                Some((Consequence::DownstreamGeneVariant, distance))
            } else {
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::make_test_transcript;
    use vep_core::variant::InputVariant;

    /// Helper to build a structural insertion variant.
    fn make_ins(start: u64, sv_end: u64) -> InputVariant {
        let mut v = InputVariant::new("21".into(), start, sv_end, b"N".to_vec(), b"-".to_vec());
        v.variant_class = vep_core::variant::VariantClass::StructuralInsertion;
        v.is_structural = true;
        v.sv_end = Some(sv_end);
        v.sv_len = Some((sv_end as i64) - (start as i64));
        v
    }

    /// Helper to build a point insertion (mobile element style, no SVLEN range).
    fn make_point_ins(pos: u64) -> InputVariant {
        let mut v = InputVariant::new("21".into(), pos, pos, b"N".to_vec(), b"-".to_vec());
        v.variant_class = vep_core::variant::VariantClass::MobileElementInsertion;
        v.is_structural = true;
        v.sv_end = Some(pos);
        v
    }

    // Test transcript layout:
    // Exon 1:   25_000_000 - 25_000_299  (5'UTR: 25_000_000-25_000_049, CDS: 25_000_050-25_000_299)
    // Intron 1: 25_000_300 - 25_001_999
    // Exon 2:   25_002_000 - 25_002_299  (all CDS)
    // Intron 2: 25_002_300 - 25_003_999
    // Exon 3:   25_004_000 - 25_006_000  (CDS: 25_004_000-25_004_299, 3'UTR: 25_004_300-25_006_000)
    // CDS: 25_000_050 - 25_004_299

    #[test]
    fn test_insertion_completely_outside_returns_none() {
        let tx = make_test_transcript();
        let v = make_ins(24_990_000, 24_994_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_none());
    }

    #[test]
    fn test_insertion_upstream() {
        let tx = make_test_transcript();
        // Insertion point is 1000bp before transcript start (within 5000bp upstream).
        let v = make_ins(24_999_000, 24_999_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::UpstreamGeneVariant));
        assert_eq!(tc.distance, Some(1000));
        assert_eq!(tc.impact, Impact::MODIFIER);
    }

    #[test]
    fn test_insertion_downstream() {
        let tx = make_test_transcript();
        // Insertion point is 2000bp after transcript end (within 5000bp downstream).
        let v = make_ins(25_008_000, 25_008_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert_eq!(tc.distance, Some(2000));
    }

    #[test]
    fn test_insertion_in_5prime_utr() {
        let tx = make_test_transcript();
        let v = make_point_ins(25_000_020);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FivePrimeUtrVariant),
            "Expected 5_prime_UTR_variant, got: {:?}",
            tc.consequences
        );
        assert!(tc.consequences.contains(&Consequence::FeatureElongation));
    }

    #[test]
    fn test_insertion_in_cds_exon() {
        let tx = make_test_transcript();
        let v = make_point_ins(25_002_100);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Expected coding_sequence_variant, got: {:?}",
            tc.consequences
        );
        assert!(tc.consequences.contains(&Consequence::FeatureElongation));
    }

    #[test]
    fn test_insertion_in_3prime_utr() {
        let tx = make_test_transcript();
        let v = make_point_ins(25_005_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
            "Expected 3_prime_UTR_variant, got: {:?}",
            tc.consequences
        );
        assert!(tc.consequences.contains(&Consequence::FeatureElongation));
    }

    #[test]
    fn test_insertion_in_intron() {
        let tx = make_test_transcript();
        let v = make_point_ins(25_001_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "Expected intron_variant, got: {:?}",
            tc.consequences
        );
        // Intronic-only insertions should not get feature_elongation.
        assert!(
            !tc.consequences.contains(&Consequence::FeatureElongation),
            "Intronic-only insertion should not get feature_elongation"
        );
    }

    #[test]
    fn test_large_insertion_spanning_cds_and_intron() {
        let tx = make_test_transcript();
        let v = make_ins(25_000_100, 25_002_200);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureElongation),
            "Expected feature_elongation, got: {:?}",
            tc.consequences
        );
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_large_insertion_fully_contains_coding_transcript() {
        // Ranged insertion spanning the entire transcript and beyond.
        // Perl VEP: when SV fully contains a coding transcript, produce
        // coding_transcript_variant as the sole consequence.
        let tx = make_test_transcript();
        let v = make_ins(24_999_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::CodingTranscriptVariant],
            "INS fully containing coding transcript should produce coding_transcript_variant only, got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences.contains(&Consequence::FeatureElongation),
            "INS extending beyond transcript should not get feature_elongation"
        );
    }

    #[test]
    fn test_tandem_repeat_fully_contains_transcript_gives_transcript_amplification() {
        // TandemRepeat (<CNV:TR>) is a copy-number gain: when it fully contains
        // a transcript, Perl VEP produces `transcript_amplification` (not a context
        // term like coding_transcript_variant). This matches duplication.rs behavior.
        let tx = make_test_transcript();
        let mut v = make_ins(24_999_000, 25_007_000);
        v.variant_class = VariantClass::TandemRepeat;
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::TranscriptAmplification],
            "TandemRepeat fully containing transcript should produce transcript_amplification, got: {:?}",
            tc.consequences
        );
        assert_eq!(tc.impact, Impact::HIGH);
        assert_eq!(
            tc.feature_overlap(v.start, v.end)
                .map(|(_, pc)| format!("{pc:.2}")),
            Some("100.00".to_string()),
            "Fully-containing TandemRepeat should show OverlapPC=100.00"
        );
    }

    #[test]
    fn test_tandem_repeat_fully_contains_noncoding_transcript_gives_transcript_amplification() {
        // Non-coding transcript fully contained by TandemRepeat should also get
        // transcript_amplification, matching Perl VEP behavior.
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;

        let mut v = make_ins(24_999_000, 25_007_000);
        v.variant_class = VariantClass::TandemRepeat;
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::TranscriptAmplification],
            "TandemRepeat fully containing non-coding transcript should produce transcript_amplification, got: {:?}",
            tc.consequences
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_regular_insertion_fully_contains_transcript_gives_context_term() {
        // Regular structural insertions (not TandemRepeat) should still get the
        // context term, not transcript_amplification.
        let tx = make_test_transcript();
        let v = make_ins(24_999_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::CodingTranscriptVariant],
            "Regular INS fully containing coding transcript should produce coding_transcript_variant, got: {:?}",
            tc.consequences
        );
        assert!(!tc
            .consequences
            .contains(&Consequence::TranscriptAmplification));
    }

    #[test]
    fn test_insertion_overlap_bp_and_pct() {
        let tx = make_test_transcript();
        // Insertion fully within transcript: 25_002_000 - 25_002_299 (300bp overlap).
        let v = make_ins(25_002_000, 25_002_299);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.feature_overlap(v.start, v.end).map(|(bp, _)| bp),
            Some(300),
            "OverlapBP should be 300"
        );
        // tx_len = 25_006_000 - 25_000_000 + 1 = 6001
        // overlap_pct = 300 / 6001 * 100 = 4.999...
        let pct = tc.feature_overlap(v.start, v.end).unwrap().1;
        assert!(
            (pct - 5.0).abs() < 0.1,
            "OverlapPC should be ~5.0, got {}",
            pct
        );
    }

    #[test]
    fn test_point_insertion_overlap_bp() {
        let tx = make_test_transcript();
        let v = make_point_ins(25_002_100);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.feature_overlap(v.start, v.end).map(|(bp, _)| bp),
            Some(1),
            "Point insertion OverlapBP should be 1"
        );
    }

    #[test]
    fn test_insertion_impact_with_feature_elongation() {
        let tx = make_test_transcript();
        // Insertion in CDS exon => feature_elongation (HIGH).
        let v = make_point_ins(25_002_100);
        let result = calculate(&v, &tx, 5000, 5000);
        let tc = result.unwrap();
        assert_eq!(
            tc.impact,
            Impact::HIGH,
            "feature_elongation should be HIGH impact"
        );
    }

    #[test]
    fn test_insertion_impact_intron_only() {
        let tx = make_test_transcript();
        // Insertion in intron only => intron_variant (MODIFIER).
        let v = make_point_ins(25_001_000);
        let result = calculate(&v, &tx, 5000, 5000);
        let tc = result.unwrap();
        assert_eq!(
            tc.impact,
            Impact::MODIFIER,
            "Intron-only insertion should be MODIFIER"
        );
    }

    #[test]
    fn test_insertion_transcript_metadata() {
        let tx = make_test_transcript();
        let v = make_point_ins(25_002_100);
        let result = calculate(&v, &tx, 5000, 5000);
        let tc = result.unwrap();
        assert_eq!(&*tc.transcript_id, "ENST00000000001");
        assert_eq!(&*tc.gene_id, "ENSG00000000001");
        assert_eq!(tc.gene_symbol.as_deref(), Some("TEST1"));
        assert_eq!(tc.biotype.as_deref(), Some("protein_coding"));
        assert!(tc.canonical);
        assert_eq!(tc.strand, 1);
    }

    #[test]
    fn test_insertion_spanning_5utr_and_cds() {
        let tx = make_test_transcript();
        let v = make_ins(25_000_020, 25_000_100);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::FivePrimeUtrVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(tc.consequences.contains(&Consequence::FeatureElongation));
    }

    #[test]
    fn test_non_coding_transcript_exon_overlap() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;

        let v = make_point_ins(25_002_100);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "Expected non_coding_transcript_exon_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_non_coding_transcript_intron_overlap() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;

        let v = make_point_ins(25_001_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_reverse_strand_upstream() {
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;

        // On reverse strand, upstream is after the transcript end.
        let v = make_ins(25_007_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::UpstreamGeneVariant),
            "Expected upstream_gene_variant for reverse strand, got: {:?}",
            tc.consequences
        );
        assert_eq!(tc.distance, Some(1000));
        assert_eq!(tc.strand, -1);
    }

    #[test]
    fn test_reverse_strand_downstream() {
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;

        // On reverse strand, downstream is before the transcript start.
        let v = make_ins(24_998_000, 24_998_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences
                .contains(&Consequence::DownstreamGeneVariant),
            "Expected downstream_gene_variant for reverse strand, got: {:?}",
            tc.consequences
        );
        assert_eq!(tc.distance, Some(2000));
    }

    #[test]
    fn test_point_insertion_in_polypyrimidine_tract() {
        let tx = make_test_transcript();
        // Point insertion in intron 1 (25_000_300 - 25_001_999), forward strand.
        // Acceptor is at intron end (25_001_999).
        // Position 25_001_990 is 9bp from acceptor (within 2..=16 range).
        let v = make_point_ins(25_001_990);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences
                .contains(&Consequence::SplicePolypyrimidineTractVariant),
            "Point insertion 9bp from acceptor should get splice_polypyrimidine_tract_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "Should also have intron_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_point_insertion_outside_polypyrimidine_tract() {
        let tx = make_test_transcript();
        // Position 25_001_000 is ~999bp from acceptor (way outside 2..=16 range).
        let v = make_point_ins(25_001_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            !tc.consequences
                .contains(&Consequence::SplicePolypyrimidineTractVariant),
            "Point insertion far from acceptor should not get polypyrimidine, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_ranged_insertion_no_polypyrimidine_endpoints_outside_tract() {
        let tx = make_test_transcript();
        // Ranged insertion (sv_start != sv_end) near the acceptor but both endpoints
        // outside the 2-16bp polypyrimidine window.
        // Intron 1: 25_000_300 - 25_001_999, acceptor at 25_001_999 (forward strand).
        // sv_start=25_001_980: dist_to_acceptor = 19 (>16, outside)
        // sv_end=25_001_998: dist_to_acceptor = 1 (<2, outside)
        let v = make_ins(25_001_980, 25_001_998);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            !tc.consequences
                .contains(&Consequence::SplicePolypyrimidineTractVariant),
            "Ranged insertion with both endpoints outside tract should not get polypyrimidine, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_ranged_insertion_endpoint_in_polypyrimidine_tract() {
        let tx = make_test_transcript();
        // Ranged tandem repeat where sv_end falls in the polypyrimidine tract.
        // Intron 1: 25_000_300 - 25_001_999, acceptor at 25_001_999 (forward strand).
        // sv_start=25_001_500: dist_to_acceptor = 499 (far outside tract)
        // sv_end=25_001_990: dist_to_acceptor = 9 (within 2..=16)
        let mut v = make_ins(25_001_500, 25_001_990);
        v.variant_class = vep_core::variant::VariantClass::TandemRepeat;
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences
                .contains(&Consequence::SplicePolypyrimidineTractVariant),
            "Ranged insertion with endpoint 9bp from acceptor should get polypyrimidine, got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "Should also have intron_variant"
        );
    }

    #[test]
    fn test_ranged_insertion_start_endpoint_in_polypyrimidine_tract() {
        let tx = make_test_transcript();
        // Ranged tandem repeat where sv_start falls in the polypyrimidine tract.
        // Intron 1: 25_000_300 - 25_001_999, acceptor at 25_001_999 (forward strand).
        // sv_start=25_001_993: dist_to_acceptor = 6 (within 2..=16)
        // sv_end=25_002_100: in exon 2 (not in intron)
        let mut v = make_ins(25_001_993, 25_002_100);
        v.variant_class = vep_core::variant::VariantClass::TandemRepeat;
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences
                .contains(&Consequence::SplicePolypyrimidineTractVariant),
            "Ranged insertion with sv_start 6bp from acceptor should get polypyrimidine, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_ranged_insertion_large_span_no_polypyrimidine() {
        let tx = make_test_transcript();
        // Ranged insertion spanning the entire intron. Neither endpoint falls in the
        // polypyrimidine tract (endpoints are at intron boundaries, not 2-16bp from acceptor).
        // Intron 1: 25_000_300 - 25_001_999
        // sv_start=25_000_200 (before intron), sv_end=25_002_100 (after intron)
        let v = make_ins(25_000_200, 25_002_100);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            !tc.consequences
                .contains(&Consequence::SplicePolypyrimidineTractVariant),
            "Large ranged insertion spanning entire intron should not get polypyrimidine, got: {:?}",
            tc.consequences
        );
    }

    // Perl VEP's `_intron_effects` trims 2bp from each intron end before the
    // "intronic" flag check, so a point ME insertion in the first or last 2bp of
    // an intron is `intergenic_variant`.
    //
    // Intron 1: 25_000_300 - 25_001_999
    //   Trimmed: 25_000_302 - 25_001_997
    // Intron 2: 25_002_300 - 25_003_999
    //   Trimmed: 25_002_302 - 25_003_997

    #[test]
    fn test_me_point_ins_at_intron_start_boundary_gives_intergenic() {
        let tx = make_test_transcript();
        // Position 25_000_300 = first base of intron 1 (donor GT dinucleotide).
        // Perl's trimmed check excludes this -> intergenic_variant.
        let v = make_point_ins(25_000_300);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntergenicVariant),
            "ME point insertion at intron start boundary should get intergenic_variant, got: {:?}",
            tc.consequences
        );
        assert!(!tc.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_me_point_ins_at_intron_start_boundary_second_base() {
        let tx = make_test_transcript();
        // Position 25_000_301 = second base of intron 1 (still in GT dinucleotide).
        let v = make_point_ins(25_000_301);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntergenicVariant),
            "ME point insertion at intron start+1 should get intergenic_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_me_point_ins_at_intron_trimmed_start_gives_intron_variant() {
        let tx = make_test_transcript();
        // Position 25_000_302 = first base inside trimmed intron range.
        // Perl considers this intronic -> intron_variant.
        let v = make_point_ins(25_000_302);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "ME point insertion at trimmed intron start should get intron_variant, got: {:?}",
            tc.consequences
        );
        assert!(!tc.consequences.contains(&Consequence::IntergenicVariant));
    }

    #[test]
    fn test_me_point_ins_at_intron_end_boundary_gives_intergenic() {
        let tx = make_test_transcript();
        // Position 25_001_999 = last base of intron 1 (acceptor AG dinucleotide).
        let v = make_point_ins(25_001_999);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntergenicVariant),
            "ME point insertion at intron end boundary should get intergenic_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_me_point_ins_at_intron_end_boundary_penultimate() {
        let tx = make_test_transcript();
        // Position 25_001_998 = second-to-last base of intron 1.
        let v = make_point_ins(25_001_998);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntergenicVariant),
            "ME point insertion at intron end-1 should get intergenic_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_me_point_ins_at_intron_trimmed_end_gives_intron_variant() {
        let tx = make_test_transcript();
        // Position 25_001_997 = last base inside trimmed intron range.
        let v = make_point_ins(25_001_997);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "ME point insertion at trimmed intron end should get intron_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_me_point_ins_deep_in_intron_gives_intron_variant() {
        let tx = make_test_transcript();
        // Position 25_001_000 = deep inside intron 1, well within trimmed range.
        let v = make_point_ins(25_001_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "ME point insertion deep in intron should get intron_variant, got: {:?}",
            tc.consequences
        );
    }

    /// The trim applies to ranged insertions too, not only to point mobile elements.
    ///
    /// Ensembl's trim in `BaseTranscriptVariationAllele.pm:143-150` is one
    /// unconditional `overlap($r_start, $r_end, $intron_start + 2, $intron_end - 2)`
    /// with no variant-class case, so it applies to ranged insertions too.
    ///
    /// The second half keeps the trim from being over-applied: reaching only the
    /// boundary base is not intronic, reaching past it is. Only the `intron_variant`
    /// term is asserted: this test does not fix what else Ensembl VEP emits for a
    /// ranged insertion whose sole transcript overlap is one boundary base.
    #[test]
    fn test_ranged_ins_intron_trim_applies_but_only_to_the_boundary_bases() {
        let tx = make_test_transcript();
        // Intron 1: 25_000_300 - 25_001_999.
        let boundary_only = calculate(&make_ins(25_000_300, 25_000_301), &tx, 5000, 5000).unwrap();
        assert!(
            !boundary_only.consequences.contains(&Consequence::IntronVariant),
            "a ranged insertion touching only the two invariant donor bases is not intronic, got: {:?}",
            boundary_only.consequences
        );
        let reaches_interior =
            calculate(&make_ins(25_000_300, 25_000_302), &tx, 5000, 5000).unwrap();
        assert!(
            reaches_interior
                .consequences
                .contains(&Consequence::IntronVariant),
            "a ranged insertion reaching the first interior base IS intronic, got: {:?}",
            reaches_interior.consequences
        );
    }

    #[test]
    fn test_me_point_ins_at_intron2_start_boundary_gives_intergenic() {
        let tx = make_test_transcript();
        // Intron 2: 25_002_300 - 25_003_999, trimmed: 25_002_302 - 25_003_997
        let v = make_point_ins(25_002_300);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntergenicVariant),
            "ME point insertion at intron 2 start boundary should get intergenic_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_me_point_ins_non_coding_intron_boundary_gives_intergenic() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;

        let v = make_point_ins(25_000_300);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntergenicVariant),
            "ME point insertion at non-coding intron boundary should get intergenic_variant, got: {:?}",
            tc.consequences
        );
        assert!(!tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
    }

    #[test]
    fn test_insertion_consequences_sorted_by_rank() {
        let tx = make_test_transcript();
        let v = make_ins(24_999_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        let tc = result.unwrap();
        // Verify consequences are sorted by rank (lower rank = more severe = first).
        for window in tc.consequences.windows(2) {
            assert!(
                window[0].rank() <= window[1].rank(),
                "Consequences not sorted: {:?} (rank {}) before {:?} (rank {})",
                window[0],
                window[0].rank(),
                window[1],
                window[1].rank()
            );
        }
    }
}
