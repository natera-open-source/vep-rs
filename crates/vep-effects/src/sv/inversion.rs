// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Structural inversion (`<INV>`) consequence calculation.
//!
//! Mirrors Perl VEP's `StructuralVariationOverlapAllele` logic for inversions.
//! Inversions disrupt the reading frame by reversing the strand of the affected
//! region but do not remove or add sequence. Perl VEP consequently never assigns
//! `transcript_ablation` for inversions, even when the INV fully contains the
//! transcript, and reports `coding_transcript_variant` for fully-contained
//! protein-coding transcripts instead.
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`, `Config.pm` being `Utils/Config.pm`).

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
use vep_core::variant::InputVariant;

/// Calculate consequences for a structural inversion against a transcript.
///
/// Inversions produce region-based overlap consequences identical to Perl VEP:
/// - `coding_sequence_variant`, `5_prime_UTR_variant`, `3_prime_UTR_variant`,
///   `intron_variant` when the INV overlaps those regions
/// - `coding_transcript_variant` for protein-coding transcripts fully contained
///   within the INV (100% overlap) where the CDS/UTR/intron overlap is not
///   individually resolved by Perl
/// - `non_coding_transcript_exon_variant`, `non_coding_transcript_variant`
///   for non-coding transcripts
/// - `upstream_gene_variant`, `downstream_gene_variant` when within distance
///   but not overlapping the transcript body
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

    let overlaps_transcript = sv_end >= tx_start && sv_start <= tx_end;

    if !overlaps_transcript {
        return calculate_distance_consequence(
            transcript,
            sv_start,
            sv_end,
            upstream_distance,
            downstream_distance,
        );
    }

    // Overlap metrics (OverlapBP/OverlapPC) use nominal coordinates.
    let overlap_start = sv_start.max(tx_start);
    let overlap_end = sv_end.min(tx_end);
    let overlap_bp = overlap_end - overlap_start + 1;
    let tx_length = tx_end - tx_start + 1;
    let overlap_pc = if tx_length > 0 {
        (overlap_bp as f64 / tx_length as f64) * 100.0
    } else {
        0.0
    };

    let mut consequences: ConsequenceList = SmallVec::new();

    // Perl VEP `coding_transcript_variant` predicate:
    //   (not coding_unknown) AND complete_overlap_feature AND within_coding_gene
    // An inversion (neither gain nor loss) fully containing the transcript gets
    // the context term as the sole consequence, with no regional breakdown.
    let fully_contains = sv_start <= tx_start && sv_end >= tx_end;

    // mature_miRNA_variant is tier 2 in Perl's OverlapConsequence table
    // (Config.pm rank 24). For a complete-overlap SV `_bvfo_preds`
    // (BaseVariationFeatureOverlapAllele.pm:465-477) sets only
    // `complete_overlap`/`within_feature`/<biotype>, so the miRNA predicate is
    // never skipped, and `get_all_OverlapConsequences` (267/277-279) pins
    // `assigned_tier = 2` before every tier-3 context term: Perl emits
    // mature_miRNA_variant alone. Tested before the `fully_contains` branch,
    // which would otherwise short-circuit it.
    if is_mature_mirna_sv(transcript, sv_start, sv_end) {
        return Some(build_transcript_consequence(
            transcript,
            smallvec![Consequence::MatureMirnaVariant],
            Impact::MODIFIER,
            overlap_bp,
            overlap_pc,
            None,
        ));
    }

    if fully_contains {
        if transcript.has_cds() {
            consequences.push(Consequence::CodingTranscriptVariant);
        } else if transcript.is_nmd_transcript() {
            consequences.push(Consequence::NmdTranscriptVariant);
        } else {
            consequences.push(Consequence::NonCodingTranscriptVariant);
        }
    } else {
        if transcript.has_cds() {
            calculate_protein_coding_consequences(&mut consequences, transcript, sv_start, sv_end);
        } else {
            calculate_non_coding_consequences(&mut consequences, transcript, sv_start, sv_end);
        }

        // NMD context: always add alongside other consequences.
        if transcript.is_nmd_transcript()
            && !consequences.contains(&Consequence::NmdTranscriptVariant)
        {
            consequences.push(Consequence::NmdTranscriptVariant);
        }

        if consequences.is_empty() {
            if transcript.has_cds() {
                consequences.push(Consequence::CodingTranscriptVariant);
            } else {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
        }
    }

    consequences.sort_by_key(|c| c.rank());
    consequences.dedup();

    let impact = consequences
        .iter()
        .map(|c| c.impact())
        .min()
        .unwrap_or(Impact::MODIFIER);

    Some(build_transcript_consequence(
        transcript,
        consequences,
        impact,
        overlap_bp,
        overlap_pc,
        None,
    ))
}

/// Calculate consequences for a protein-coding transcript.
fn calculate_protein_coding_consequences(
    consequences: &mut ConsequenceList,
    transcript: &Transcript,
    sv_start: u64,
    sv_end: u64,
) {
    if overlaps_cds_exon(transcript, sv_start, sv_end) {
        push_unique(consequences, Consequence::CodingSequenceVariant);
    }

    if overlaps_five_prime_utr(transcript, sv_start, sv_end) {
        push_unique(consequences, Consequence::FivePrimeUtrVariant);
    }

    if overlaps_three_prime_utr(transcript, sv_start, sv_end) {
        push_unique(consequences, Consequence::ThreePrimeUtrVariant);
    }

    // Perl's `within_intron` reads the `intronic` flag from `_intron_effects`,
    // set with the intron span trimmed by 2bp at each end
    // (BaseTranscriptVariationAllele.pm:143-150), so an inversion touching only
    // the donor/acceptor dinucleotides is not `intron_variant`.
    if overlaps_any_intron_trimmed(transcript, sv_start, sv_end) {
        push_unique(consequences, Consequence::IntronVariant);
    }
}

/// Calculate consequences for a non-coding transcript.
fn calculate_non_coding_consequences(
    consequences: &mut ConsequenceList,
    transcript: &Transcript,
    sv_start: u64,
    sv_end: u64,
) {
    let overlaps_exon = overlaps_any_exon(transcript, sv_start, sv_end);
    // Trimmed intron span, as in `calculate_protein_coding_consequences`.
    let overlaps_intron = overlaps_any_intron_trimmed(transcript, sv_start, sv_end);

    // mature_miRNA overlap already returned above (tier-2 short-circuit), so any
    // transcript reaching here is not a mature-miRNA hit.
    if overlaps_exon {
        push_unique(consequences, Consequence::NonCodingTranscriptExonVariant);
    }
    if overlaps_intron {
        push_unique(consequences, Consequence::IntronVariant);
    }
    // Perl VEP adds non_coding_transcript_variant to a non-coding transcript when
    // there is NO exon overlap and not mature_miRNA.
    if !overlaps_exon && !transcript.is_nmd_transcript() {
        push_unique(consequences, Consequence::NonCodingTranscriptVariant);
    }
}

/// Calculate upstream or downstream consequence for an SV that doesn't overlap
/// the transcript body.
fn calculate_distance_consequence(
    transcript: &Transcript,
    sv_start: u64,
    sv_end: u64,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    let tx_start = transcript.start;
    let tx_end = transcript.end;

    let distance = if sv_end < tx_start {
        tx_start - sv_end
    } else if sv_start > tx_end {
        sv_start - tx_end
    } else {
        return None;
    };

    let is_upstream = match transcript.strand {
        Strand::Forward => sv_end < tx_start,
        Strand::Reverse => sv_start > tx_end,
    };

    if is_upstream && distance <= upstream_distance {
        Some(build_transcript_consequence(
            transcript,
            smallvec![Consequence::UpstreamGeneVariant],
            Impact::MODIFIER,
            0,
            0.0,
            Some(distance),
        ))
    } else if !is_upstream && distance <= downstream_distance {
        Some(build_transcript_consequence(
            transcript,
            smallvec![Consequence::DownstreamGeneVariant],
            Impact::MODIFIER,
            0,
            0.0,
            Some(distance),
        ))
    } else {
        None
    }
}

/// Build a `TranscriptConsequence` with SV-specific fields.
fn build_transcript_consequence(
    transcript: &Transcript,
    consequences: ConsequenceList,
    impact: Impact,
    _overlap_bp: u64,
    _overlap_pc: f64,
    distance: Option<u64>,
) -> TranscriptConsequence {
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
        mane_select: transcript.mane_select.clone(),
        mane_plus_clinical: transcript.mane_plus_clinical.clone(),
        appris: transcript.appris.clone(),
        ccds: transcript.ccds.clone(),
        swissprot: transcript.swissprot.clone(),
        trembl: transcript.trembl.clone(),
        refseq: transcript.refseq.clone(),
        ..TranscriptConsequence::default()
    }
}

fn push_unique(consequences: &mut ConsequenceList, consequence: Consequence) {
    if !consequences.contains(&consequence) {
        consequences.push(consequence);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Shared with `sv::tests` so the tier-ordering fixture has one definition.
    //
    use crate::sv::tests::make_mature_mirna_transcript;
    use crate::test_helpers::{make_test_transcript, make_test_transcript_with_flags};

    use vep_core::variant::VariantClass;

    /// Helper: build a structural inversion variant.
    fn make_inversion(start: u64, sv_end: u64) -> InputVariant {
        let mut v = InputVariant::new("21".into(), start, sv_end, b"N".to_vec(), b"-".to_vec());
        v.variant_class = VariantClass::Inversion;
        v.is_structural = true;
        v.sv_end = Some(sv_end);
        v.allele_string = format!("{}/{}", "N", "inversion");
        v
    }

    // Test transcript layout (from make_test_transcript):
    //   Exon 1:   25_000_000 - 25_000_299  (300bp)
    //   Intron 1: 25_000_300 - 25_001_999  (1700bp)
    //   Exon 2:   25_002_000 - 25_002_299  (300bp)
    //   Intron 2: 25_002_300 - 25_003_999  (1700bp)
    //   Exon 3:   25_004_000 - 25_006_000  (2001bp)
    //
    //   CDS: coding_region_start=25_000_050, coding_region_end=25_004_299
    //   5' UTR: 25_000_000 - 25_000_049 (exon 1 portion)
    //   3' UTR: 25_004_300 - 25_006_000 (exon 3 portion)

    /// Perl's 2bp-trimmed intron span.
    ///
    /// Perl's `within_intron` (VariationEffect.pm:629) returns the `intronic` flag
    /// from `_intron_effects`, which is set with the intron span trimmed 2bp at each
    /// end (BaseTranscriptVariationAllele.pm:143-150):
    ///   `overlap($r_start, $r_end, $intron_start + 2, $intron_end - 2)`
    /// An inversion whose only intron contact is the 2bp splice donor/acceptor
    /// dinucleotides must therefore not get `intron_variant`.
    ///
    /// Intron 1 of `make_test_transcript` is 25_000_300 - 25_001_999, so the trimmed
    /// interior is 25_000_302 - 25_001_997. An INV ending at 25_000_301 touches only
    /// the donor dinucleotide.
    ///
    /// An untrimmed intron test would report `intron_variant` here.
    #[test]
    fn test_inv_touching_only_donor_dinucleotide_has_no_intron_variant() {
        let tx = make_test_transcript();
        // Exon 1 = 25_000_000-25_000_299; INV reaches 2bp into intron 1.
        let v = make_inversion(25_000_200, 25_000_301);
        let csq = calculate(&v, &tx, 5000, 5000).expect("INV must annotate the transcript");
        assert!(
            !csq.consequences.contains(&Consequence::IntronVariant),
            "INV touching only the 2bp donor dinucleotide must NOT be intron_variant \
             (Perl trims intron_start+2), got: {:?}",
            csq.consequences
        );
        // The exonic overlap is still reported.
        assert!(
            csq.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "exonic CDS overlap must still be reported, got: {:?}",
            csq.consequences
        );
    }

    /// Non-coding counterpart of the donor-dinucleotide test above.
    ///
    /// `overlaps_any_intron_trimmed` is called at two sites: once in
    /// `calculate_protein_coding_consequences` and once in
    /// `calculate_non_coding_consequences`. The coding-site test alone does not
    /// cover the non-coding site: with only that site untrimmed the coding tests
    /// stay green while the non-coding path yields `[IntronVariant,
    /// NonCodingTranscriptVariant]`, so each site has its own test.
    #[test]
    fn test_inv_non_coding_touching_only_donor_dinucleotide_has_no_intron_variant() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.protein_id = None;

        // Exon 1 = 25_000_000-25_000_299; INV reaches exactly 2bp into intron 1,
        // i.e. only the donor dinucleotide, which Perl's `intronic` flag trims.
        let v = make_inversion(25_000_200, 25_000_301);
        let csq = calculate(&v, &tx, 5000, 5000).expect("INV must annotate the transcript");
        assert!(
            !csq.consequences.contains(&Consequence::IntronVariant),
            "non-coding INV touching only the 2bp donor dinucleotide must NOT be \
             intron_variant (Perl trims intron_start+2), got: {:?}",
            csq.consequences
        );
        // The exonic overlap on the non-coding path is still reported.
        assert!(
            csq.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "non-coding exonic overlap must still be reported, got: {:?}",
            csq.consequences
        );
    }

    /// Companion to the above: one base further in (25_000_302, the first trimmed
    /// interior base) does yield `intron_variant`. Proves the trim narrows the
    /// window by exactly 2bp rather than disabling the predicate.
    #[test]
    fn test_inv_reaching_trimmed_intron_interior_has_intron_variant() {
        let tx = make_test_transcript();
        let v = make_inversion(25_000_200, 25_000_302);
        let csq = calculate(&v, &tx, 5000, 5000).expect("INV must annotate the transcript");
        assert!(
            csq.consequences.contains(&Consequence::IntronVariant),
            "INV reaching intron_start+2 must be intron_variant, got: {:?}",
            csq.consequences
        );
    }

    /// Acceptor side of the trimmed intron boundary: intron 1 ends at 25_001_999, so the trimmed
    /// interior ends at 25_001_997. An INV starting at 25_001_998 touches only the
    /// acceptor dinucleotide and must not be `intron_variant`.
    #[test]
    fn test_inv_touching_only_acceptor_dinucleotide_has_no_intron_variant() {
        let tx = make_test_transcript();
        let v = make_inversion(25_001_998, 25_002_100);
        let csq = calculate(&v, &tx, 5000, 5000).expect("INV must annotate the transcript");
        assert!(
            !csq.consequences.contains(&Consequence::IntronVariant),
            "INV touching only the 2bp acceptor dinucleotide must NOT be intron_variant, \
             got: {:?}",
            csq.consequences
        );
    }

    /// The miRNA tier-2 short-circuit on the inversion path.
    ///
    /// An inversion fully containing a mature miRNA must produce
    /// `mature_miRNA_variant` alone (Perl tier-2 short-circuit), not the tier-3
    /// `non_coding_transcript_variant`.
    ///
    /// Without the `is_mature_mirna_sv` early return above the `fully_contains`
    /// branch this assert fails with `[NonCodingTranscriptVariant]`.
    #[test]
    fn test_inv_fully_containing_mirna_gives_mature_mirna_only() {
        let tx = make_mature_mirna_transcript();
        let v = make_inversion(24_999_000, 25_001_000);
        let csq = calculate(&v, &tx, 5000, 5000).expect("INV must annotate the miRNA transcript");
        assert_eq!(
            csq.consequences.to_vec(),
            vec![Consequence::MatureMirnaVariant],
            "INV fully containing a mature miRNA must be mature_miRNA_variant alone, got: {:?}",
            csq.consequences
        );
        assert!(
            !csq.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "must not emit the tier-3 non_coding_transcript_variant context term"
        );
    }

    #[test]
    fn test_inv_spanning_entire_transcript() {
        // INV fully contains the transcript: Perl VEP assigns coding_transcript_variant
        // as the sole consequence (not transcript_ablation, no regional breakdown).
        let tx = make_test_transcript();
        let v = make_inversion(24_999_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert_eq!(
            csq.consequences.to_vec(),
            vec![Consequence::CodingTranscriptVariant],
            "INV fully containing coding transcript should produce coding_transcript_variant only, got: {:?}",
            csq.consequences
        );
        assert_eq!(csq.impact, Impact::MODIFIER);
        assert!(
            !csq.consequences.contains(&Consequence::TranscriptAblation),
            "Inversions must not produce transcript_ablation"
        );
    }

    #[test]
    fn test_inv_overlapping_cds_and_utr() {
        let tx = make_test_transcript();
        let v = make_inversion(25_002_100, 25_005_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Expected coding_sequence_variant, got: {:?}",
            csq.consequences
        );
        assert!(
            csq.consequences
                .contains(&Consequence::ThreePrimeUtrVariant),
            "Expected 3_prime_UTR_variant, got: {:?}",
            csq.consequences
        );
        assert!(
            csq.consequences.contains(&Consequence::IntronVariant),
            "Expected intron_variant, got: {:?}",
            csq.consequences
        );
    }

    #[test]
    fn test_inv_intron_only() {
        let tx = make_test_transcript();
        let v = make_inversion(25_000_500, 25_001_500);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences.contains(&Consequence::IntronVariant),
            "Expected intron_variant, got: {:?}",
            csq.consequences
        );
        assert!(
            !csq.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Should not have coding_sequence_variant for intron-only overlap"
        );
    }

    #[test]
    fn test_inv_five_prime_utr_only() {
        let tx = make_test_transcript();
        let v = make_inversion(25_000_010, 25_000_040);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences.contains(&Consequence::FivePrimeUtrVariant),
            "Expected 5_prime_UTR_variant, got: {:?}",
            csq.consequences
        );
        assert!(
            !csq.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Should not have coding_sequence_variant for UTR-only overlap"
        );
    }

    #[test]
    fn test_inv_three_prime_utr_only() {
        let tx = make_test_transcript();
        let v = make_inversion(25_004_500, 25_005_500);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences
                .contains(&Consequence::ThreePrimeUtrVariant),
            "Expected 3_prime_UTR_variant, got: {:?}",
            csq.consequences
        );
        assert!(
            !csq.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Should not have coding_sequence_variant for UTR-only overlap"
        );
    }

    #[test]
    fn test_inv_reverse_strand_cds_start_nf_adds_three_prime_utr() {
        let mut tx = make_test_transcript_with_flags(&["cds_start_NF", "cds_end_NF"]);
        tx.strand = Strand::Reverse;
        tx.coding_region_start = Some(tx.start);

        let v = make_inversion(tx.start, 25_001_000);
        let csq = calculate(&v, &tx, 5000, 5000).unwrap();

        assert!(csq
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(csq
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant));
        assert!(csq.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_inv_upstream() {
        let tx = make_test_transcript();
        let v = make_inversion(24_996_000, 24_998_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences.contains(&Consequence::UpstreamGeneVariant),
            "Expected upstream_gene_variant, got: {:?}",
            csq.consequences
        );
        assert!(csq.distance.is_some());
        assert_eq!(csq.distance.unwrap(), 2000); // 25_000_000 - 24_998_000
    }

    #[test]
    fn test_inv_downstream() {
        let tx = make_test_transcript();
        let v = make_inversion(25_007_000, 25_009_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences
                .contains(&Consequence::DownstreamGeneVariant),
            "Expected downstream_gene_variant, got: {:?}",
            csq.consequences
        );
        assert!(csq.distance.is_some());
        assert_eq!(csq.distance.unwrap(), 1000); // 25_007_000 - 25_006_000
    }

    #[test]
    fn test_inv_too_far_away() {
        let tx = make_test_transcript();
        let v = make_inversion(24_990_000, 24_993_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_none(), "INV beyond distance should return None");
    }

    #[test]
    fn test_inv_non_coding_transcript_exon_overlap() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.protein_id = None;

        let v = make_inversion(25_000_100, 25_000_200);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "Expected non_coding_transcript_exon_variant, got: {:?}",
            csq.consequences
        );
    }

    #[test]
    fn test_inv_non_coding_intron_only() {
        // Non-coding transcript with intron-only overlap.
        // Perl VEP adds both intron_variant and non_coding_transcript_variant.
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.protein_id = None;

        let v = make_inversion(25_000_500, 25_001_500);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences.contains(&Consequence::IntronVariant),
            "Expected intron_variant, got: {:?}",
            csq.consequences
        );
        assert!(
            csq.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "Expected non_coding_transcript_variant for intron-only non-coding overlap, got: {:?}",
            csq.consequences
        );
    }

    #[test]
    fn test_inv_non_coding_fully_contained() {
        // Non-coding transcript fully within INV: Perl VEP assigns
        // non_coding_transcript_variant as the sole consequence.
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.protein_id = None;
        tx.exons = vec![tx.exons[0].clone()];
        tx.introns.clear();
        if let Some(ref mut vefc) = tx.vefc {
            vefc.introns.clear();
        }
        tx.end = 25_000_299;

        let v = make_inversion(24_999_000, 25_001_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert_eq!(
            csq.consequences.to_vec(),
            vec![Consequence::NonCodingTranscriptVariant],
            "INV fully containing non-coding transcript should produce non_coding_transcript_variant only, got: {:?}",
            csq.consequences
        );
    }

    #[test]
    fn test_inv_cds_overlap_only() {
        let tx = make_test_transcript();
        let v = make_inversion(25_002_050, 25_002_250);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Expected coding_sequence_variant, got: {:?}",
            csq.consequences
        );
        assert!(
            !csq.consequences.contains(&Consequence::IntronVariant),
            "Should not have intron_variant for exon-only overlap"
        );
    }

    #[test]
    fn test_inv_transcript_metadata() {
        let tx = make_test_transcript();
        let v = make_inversion(25_000_100, 25_000_200);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert_eq!(&*csq.transcript_id, "ENST00000000001");
        assert_eq!(&*csq.gene_id, "ENSG00000000001");
        assert_eq!(csq.gene_symbol.as_deref(), Some("TEST1"));
        assert_eq!(csq.biotype.as_deref(), Some("protein_coding"));
        assert!(csq.canonical);
        assert_eq!(csq.strand, 1);
        assert_eq!(csq.feature_type, FeatureType::Transcript);
    }

    #[test]
    fn test_inv_consequences_sorted_by_rank() {
        let tx = make_test_transcript();
        let v = make_inversion(24_999_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        for window in csq.consequences.windows(2) {
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

    #[test]
    fn test_inv_reverse_strand_upstream_downstream() {
        // Reverse-strand transcript: upstream is at higher genomic coords.
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;

        // SV beyond transcript end (higher genomic coords) = upstream for reverse strand.
        let v = make_inversion(25_007_000, 25_009_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences.contains(&Consequence::UpstreamGeneVariant),
            "Higher genomic coords should be upstream for reverse strand, got: {:?}",
            csq.consequences
        );

        // SV before transcript start (lower genomic coords) = downstream for reverse strand.
        let v2 = make_inversion(24_996_000, 24_998_000);
        let result2 = calculate(&v2, &tx, 5000, 5000);
        assert!(result2.is_some());
        let csq2 = result2.unwrap();
        assert!(
            csq2.consequences
                .contains(&Consequence::DownstreamGeneVariant),
            "Lower genomic coords should be downstream for reverse strand, got: {:?}",
            csq2.consequences
        );
    }

    #[test]
    fn test_inv_no_sv_end_falls_back_to_end() {
        // If sv_end is None, should use variant.end.
        let mut v = InputVariant::new(
            "21".into(),
            25_000_100,
            25_000_200,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::Inversion;
        v.is_structural = true;
        v.sv_end = None;

        let tx = make_test_transcript();
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
    }

    #[test]
    fn test_inv_coding_transcript_variant_fallback() {
        // A protein-coding transcript with only intron overlap outside the CDS
        // falls back to coding_transcript_variant.
        let mut tx = make_test_transcript();
        tx.coding_region_start = None;
        tx.coding_region_end = None;

        let v = make_inversion(25_000_500, 25_001_500);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences.contains(&Consequence::IntronVariant),
            "Expected intron_variant, got: {:?}",
            csq.consequences
        );
    }

    #[test]
    fn test_inv_nmd_transcript() {
        // NMD transcript with overlap should get NMD_transcript_variant as context.
        let mut tx = make_test_transcript();
        tx.biotype = "nonsense_mediated_decay".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.protein_id = None;
        tx.exons = vec![];
        tx.introns.clear();
        if let Some(ref mut vefc) = tx.vefc {
            vefc.introns.clear();
            vefc.sorted_exons.clear();
        }

        let v = make_inversion(25_000_100, 25_000_200);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let csq = result.unwrap();
        assert!(
            csq.consequences
                .contains(&Consequence::NmdTranscriptVariant),
            "Expected NMD_transcript_variant, got: {:?}",
            csq.consequences
        );
    }
}
