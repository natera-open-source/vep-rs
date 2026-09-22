// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Copy number variation (`<CNV>`, `<CN0>`–`<CN9>`, `<CPX>`) consequence calculation.
//!
//! Perl VEP's `VCF.pm` classifies CN alleles by their structural effect:
//!
//! - **CN0** → deletion class: transcript_ablation if fully contained, otherwise
//!   feature_truncation + regional consequences (like `<DEL>`).
//! - **CN2** → duplication class: transcript_amplification if fully contained,
//!   otherwise feature_elongation + regional consequences (like `<DUP>`).
//! - **CN1, CN3–CN9, generic `<CNV>`** → copy_number_variation class: regional
//!   overlap consequences only. No feature_truncation/elongation/amplification/ablation.
//!
//! `VCF.pm` is `Bio/EnsEMBL/VEP/Parser/VCF.pm` in ensembl-vep release/115; the other
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`, `Config.pm` being `Utils/Config.pm`).

use super::{
    is_mature_mirna_sv, overlaps_any_exon, overlaps_any_intron_trimmed, overlaps_cds_exon,
    overlaps_five_prime_utr, overlaps_three_prime_utr, transcript_has_incomplete_cds,
};
use smallvec::{smallvec, SmallVec};
use vep_core::consequence::{
    Consequence, ConsequenceList, FeatureType, Impact, TranscriptConsequence,
};
use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

/// Structural effect class inferred from the ALT allele, matching Perl VEP's
/// `VCF.pm` classification of CN alleles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CopyNumberClass {
    /// CN0: deletion (0 copies), same behavior as `<DEL>`
    Deletion,
    /// CN2: duplication (2 copies), same behavior as `<DUP>`
    Duplication,
    /// CN1, CN3–CN9, generic `<CNV>`: ambiguous copy number variation
    Generic,
}

/// Parse the copy number from the ALT allele string.
///
/// Matches Perl VEP's `VCF.pm` classification:
///   CN0 → Deletion, CN2 → Duplication, everything else → Generic.
fn parse_copy_number(variant: &InputVariant) -> CopyNumberClass {
    // Use the raw alt allele bytes (not display_allele() which returns the SO term for SVs).
    let alt = String::from_utf8_lossy(variant.alt_allele());
    let upper = alt.to_ascii_uppercase();

    if upper.starts_with("<CN") && upper.ends_with('>') {
        let inner = &upper[3..upper.len() - 1]; // strip "<CN" and ">"
        let num_str = inner.strip_prefix('=').unwrap_or(inner);
        if let Ok(cn) = num_str.parse::<u32>() {
            return match cn {
                0 => CopyNumberClass::Deletion,
                2 => CopyNumberClass::Duplication,
                _ => CopyNumberClass::Generic,
            };
        }
    }

    CopyNumberClass::Generic
}

/// Calculate consequences for a copy number variant against a transcript.
///
/// CNVs use span-based overlap (start to sv_end) like deletions/duplications.
/// The specific consequences depend on the inferred copy number class.
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

    let cn_class = parse_copy_number(variant);

    let fully_contains = sv_start <= tx_start && sv_end >= tx_end;

    // mature_miRNA_variant is tier 2 in Perl's OverlapConsequence table
    // (Config.pm rank 24). For a complete-overlap SV `_bvfo_preds`
    // (BaseVariationFeatureOverlapAllele.pm:465-477) sets only
    // `complete_overlap`/`within_feature`/<biotype>, so the miRNA predicate is
    // never skipped, and `get_all_OverlapConsequences` (lines 267/277-279)
    // records `assigned_tier = 2` and stops before every tier-3 term: Perl
    // emits mature_miRNA_variant alone for a CNV containing a miRNA transcript.
    if is_mature_mirna_sv(transcript, sv_start, sv_end) {
        return Some(build_consequence(
            transcript,
            smallvec![Consequence::MatureMirnaVariant],
            Impact::MODIFIER,
            None,
        ));
    }

    if cn_class == CopyNumberClass::Deletion && fully_contains {
        return Some(build_consequence(
            transcript,
            smallvec![Consequence::TranscriptAblation],
            Impact::HIGH,
            None,
        ));
    }

    if cn_class == CopyNumberClass::Duplication && fully_contains {
        return Some(build_consequence(
            transcript,
            smallvec![Consequence::TranscriptAmplification],
            Impact::HIGH,
            None,
        ));
    }

    // For generic CNVs (CN1, CN3–CN9, <CNV>) that fully contain the transcript,
    // Perl uses only the context term (coding_transcript_variant, etc.).
    if fully_contains {
        let consequence = if transcript.is_nmd_transcript() {
            Consequence::NmdTranscriptVariant
        } else if transcript.has_cds() {
            Consequence::CodingTranscriptVariant
        } else {
            Consequence::NonCodingTranscriptVariant
        };
        let impact = consequence.impact();
        return Some(build_consequence(
            transcript,
            smallvec![consequence],
            impact,
            None,
        ));
    }

    if sv_end < tx_start || sv_start > tx_end {
        return compute_upstream_downstream(
            transcript,
            sv_start,
            sv_end,
            upstream_distance,
            downstream_distance,
        );
    }

    let mut consequences: ConsequenceList = SmallVec::new();
    let is_deletion = cn_class == CopyNumberClass::Deletion;
    let is_duplication = cn_class == CopyNumberClass::Duplication;

    // Perl VEP predicates (only for CN0=deletion and CN2=duplication classes):
    // feature_truncation = within_cdna AND (partial_overlap OR complete_within) AND deletion
    // feature_elongation = within_cdna AND complete_within_feature AND duplication
    let within_cdna = overlaps_any_exon(transcript, sv_start, sv_end);
    let complete_within = sv_start >= tx_start && sv_end <= tx_end;
    let extends_beyond = sv_start < tx_start || sv_end > tx_end;

    if is_deletion && within_cdna && (extends_beyond || complete_within) {
        consequences.push(Consequence::FeatureTruncation);
    }
    if is_duplication && within_cdna && complete_within {
        consequences.push(Consequence::FeatureElongation);
    }

    if is_deletion {
        if let (Some(cds_start), Some(cds_end)) =
            (transcript.coding_region_start, transcript.coding_region_end)
        {
            let has_complete_cds =
                transcript.has_cds() && !transcript_has_incomplete_cds(transcript);
            let (stop_start, stop_end) = match transcript.strand {
                Strand::Forward => (cds_end.saturating_sub(2), cds_end),
                Strand::Reverse => (cds_start, cds_start + 2),
            };
            if sv_start <= stop_end && sv_end >= stop_start && has_complete_cds {
                consequences.push(Consequence::StopLost);
            }
        }
    }

    let is_coding = transcript.has_cds();

    if is_coding {
        // A complete and an incomplete CDS get the same term:
        // `CodingTranscriptVariant` is not emitted on the CNV path.
        if overlaps_cds_exon(transcript, sv_start, sv_end) {
            consequences.push(Consequence::CodingSequenceVariant);
        }

        if overlaps_five_prime_utr(transcript, sv_start, sv_end) {
            consequences.push(Consequence::FivePrimeUtrVariant);
        }

        if overlaps_three_prime_utr(transcript, sv_start, sv_end) {
            consequences.push(Consequence::ThreePrimeUtrVariant);
        }
    }

    // Exon overlap for non-coding transcripts.
    // mature_miRNA overlap already returned above (tier-2 short-circuit), so any
    // transcript reaching here is not a mature-miRNA hit.
    if !is_coding && overlaps_any_exon(transcript, sv_start, sv_end) {
        consequences.push(Consequence::NonCodingTranscriptExonVariant);
    }

    // Intron overlap on Ensembl's trimmed boundaries: `_intron_effects` sets
    // `intronic` with `overlap($r_start, $r_end, $intron_start + 2, $intron_end - 2)`
    // (`BaseTranscriptVariationAllele.pm:143-150`), so an SV touching only the four
    // invariant donor/acceptor bases gets no `intron_variant`.
    if overlaps_any_intron_trimmed(transcript, sv_start, sv_end) {
        consequences.push(Consequence::IntronVariant);
    }

    // Context terms for non-coding transcripts: always added alongside other consequences.
    // Perl's SV annotation adds non_coding_transcript_variant for non-coding transcripts
    // even when other consequences (intron_variant, etc.) are already present.
    // Perl's within_non_coding_gene excludes mature_miRNA.
    if !is_coding
        && !transcript.is_nmd_transcript()
        && !consequences.iter().any(|c| {
            matches!(
                c,
                Consequence::NonCodingTranscriptExonVariant | Consequence::MatureMirnaVariant
            )
        })
    {
        consequences.push(Consequence::NonCodingTranscriptVariant);
    }

    // NMD transcript context is always added alongside other consequences.
    if transcript.is_nmd_transcript() {
        consequences.push(Consequence::NmdTranscriptVariant);
    }

    // Coding transcript context: only as fallback when no other consequences assigned.
    // Perl's SV annotation does not add coding_transcript_variant alongside
    // intron_variant or coding_sequence_variant; it is only a sole-consequence fallback.
    if consequences.is_empty() {
        if is_coding {
            consequences.push(Consequence::CodingTranscriptVariant);
        } else {
            consequences.push(Consequence::NonCodingTranscriptVariant);
        }
    }

    consequences.sort_by_key(|c| c.rank());
    consequences.dedup();

    let impact = consequences
        .first()
        .map(|c| c.impact())
        .unwrap_or(Impact::MODIFIER);

    Some(build_consequence(transcript, consequences, impact, None))
}

/// Compute upstream/downstream gene variant consequence.
fn compute_upstream_downstream(
    transcript: &Transcript,
    sv_start: u64,
    sv_end: u64,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    let tx_start = transcript.start;
    let tx_end = transcript.end;

    let consequence;
    let distance;

    match transcript.strand {
        Strand::Forward => {
            if sv_end < tx_start && tx_start.saturating_sub(sv_end) <= upstream_distance {
                consequence = Consequence::UpstreamGeneVariant;
                distance = tx_start - sv_end;
            } else if sv_start > tx_end && sv_start.saturating_sub(tx_end) <= downstream_distance {
                consequence = Consequence::DownstreamGeneVariant;
                distance = sv_start - tx_end;
            } else {
                return None;
            }
        }
        Strand::Reverse => {
            if sv_start > tx_end && sv_start.saturating_sub(tx_end) <= upstream_distance {
                consequence = Consequence::UpstreamGeneVariant;
                distance = sv_start - tx_end;
            } else if sv_end < tx_start && tx_start.saturating_sub(sv_end) <= downstream_distance {
                consequence = Consequence::DownstreamGeneVariant;
                distance = tx_start - sv_end;
            } else {
                return None;
            }
        }
    }

    Some(build_consequence(
        transcript,
        smallvec![consequence],
        Impact::MODIFIER,
        Some(distance),
    ))
}

/// Build a `TranscriptConsequence` from common transcript fields.
fn build_consequence(
    transcript: &Transcript,
    consequences: ConsequenceList,
    impact: Impact,
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
        cdna_position: None,
        cds_position: None,
        protein_position: None,
        amino_acids: None,
        codons: None,
        protein_id: transcript.protein_id.clone(),
        distance,
        strand: transcript.strand.as_i8(),
        exon: None,
        intron: None,
        hgvsc: None,
        hgvsp: None,
        sift: None,
        polyphen: None,
        domains: Vec::new(),
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
        plugin_data: indexmap::IndexMap::new(),
        loftee_ctx: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Shared with `sv::tests` so the tier-ordering fixture has one definition.
    //
    use crate::sv::tests::make_mature_mirna_transcript;
    use crate::test_helpers::make_test_transcript;

    use vep_core::variant::VariantClass;

    /// Helper: build a CNV variant with specified alt allele and sv_end.
    fn make_cnv(alt: &[u8], start: u64, sv_end: u64) -> InputVariant {
        let mut v = InputVariant::new("21".into(), start, sv_end, b"N".to_vec(), alt.to_vec());
        v.variant_class = VariantClass::CopyNumberVariation;
        v.is_structural = true;
        v.sv_end = Some(sv_end);
        v
    }

    // --- Test transcript layout (from make_test_transcript) ---
    // Transcript: chr21:25_000_000 - 25_006_000, forward strand, protein_coding
    // Exon 1:   25_000_000 - 25_000_299
    // Intron 1: 25_000_300 - 25_001_999
    // Exon 2:   25_002_000 - 25_002_299
    // Intron 2: 25_002_300 - 25_003_999
    // Exon 3:   25_004_000 - 25_006_000
    // CDS: 25_000_050 - 25_004_299
    // 5' UTR: 25_000_000 - 25_000_049
    // 3' UTR: 25_004_300 - 25_006_000

    /// The miRNA tier-2 short-circuit.
    ///
    /// A generic `<CNV>` that fully contains a miRNA transcript must produce
    /// `mature_miRNA_variant` alone, not `non_coding_transcript_variant`.
    ///
    /// Perl mechanism: `mature_miRNA_variant` is tier 2 (Config.pm rank 24). For a
    /// complete-overlap SV, `_bvfo_preds`
    /// (BaseVariationFeatureOverlapAllele.pm:465-477) emits only
    /// `complete_overlap`/`within_feature`/<biotype>, so the miRNA predicate is not
    /// skipped by `_skip_oc`; `get_all_OverlapConsequences` then sets
    /// `assigned_tier = 2` and stops before the tier-3 context terms.
    ///
    /// With `is_mature_mirna_sv` tested below the `fully_contains` short-circuit
    /// this assert fails with `[NonCodingTranscriptVariant]`.
    #[test]
    fn test_cnv_fully_containing_mirna_gives_mature_mirna_only() {
        let tx = make_mature_mirna_transcript();
        // CNV strictly contains the whole transcript (25_000_000-25_000_100).
        let v = make_cnv(b"<CNV>", 24_999_000, 25_001_000);
        let tc = calculate(&v, &tx, 5000, 5000).expect("CNV must annotate the miRNA transcript");
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::MatureMirnaVariant],
            "CNV fully containing a mature miRNA must be mature_miRNA_variant ALONE \
             (Perl tier-2 short-circuit), got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "must not emit the tier-3 non_coding_transcript_variant context term"
        );
    }

    /// The miRNA tier-2 short-circuit, CN0 (deletion-class) variant.
    ///
    /// The tier-2 miRNA short-circuit outranks even `transcript_ablation` (tier 1
    /// but gated on `deletion => 1, complete_overlap => 1`): Perl's output
    /// for a CN0 fully covering a mature miRNA is `mature_miRNA_variant`, because
    /// the miRNA OC is reached and pins the tier before the ablation candidate is
    /// accepted for a non-coding feature.
    ///
    /// With ablation checked first this assert fails with `[TranscriptAblation]`.
    #[test]
    fn test_cn0_fully_containing_mirna_gives_mature_mirna_not_ablation() {
        let tx = make_mature_mirna_transcript();
        let v = make_cnv(b"<CN0>", 24_999_000, 25_001_000);
        let tc = calculate(&v, &tx, 5000, 5000).expect("CN0 must annotate the miRNA transcript");
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::MatureMirnaVariant],
            "CN0 fully containing a mature miRNA must be mature_miRNA_variant alone, got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences.contains(&Consequence::TranscriptAblation),
            "must not emit transcript_ablation for a mature-miRNA hit"
        );
    }

    /// Guard: a CNV fully containing a miRNA transcript but missing the mature
    /// subfeature window still uses the tier-3 context term. This pins the rule to
    /// the mature-miRNA window rather than to the miRNA biotype.
    #[test]
    fn test_cnv_containing_mirna_tx_without_mature_overlap_uses_context_term() {
        let mut tx = make_mature_mirna_transcript();
        // Move the mature window off the transcript so is_mature_mirna_sv is false.
        tx.attributes[0].value = "5000-5010".into();
        let v = make_cnv(b"<CNV>", 24_999_000, 25_001_000);
        let tc = calculate(&v, &tx, 5000, 5000).expect("CNV must annotate the transcript");
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::NonCodingTranscriptVariant],
            "no mature-miRNA overlap => tier-3 non_coding_transcript_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_parse_copy_number_cn0() {
        let v = make_cnv(b"<CN0>", 100, 200);
        assert_eq!(parse_copy_number(&v), CopyNumberClass::Deletion);
    }

    #[test]
    fn test_parse_copy_number_cn1() {
        let v = make_cnv(b"<CN1>", 100, 200);
        assert_eq!(parse_copy_number(&v), CopyNumberClass::Generic);
    }

    #[test]
    fn test_parse_copy_number_cn2() {
        let v = make_cnv(b"<CN2>", 100, 200);
        assert_eq!(parse_copy_number(&v), CopyNumberClass::Duplication);
    }

    #[test]
    fn test_parse_copy_number_cn3_gain() {
        let v = make_cnv(b"<CN3>", 100, 200);
        assert_eq!(parse_copy_number(&v), CopyNumberClass::Generic);
    }

    #[test]
    fn test_parse_copy_number_cn9_gain() {
        let v = make_cnv(b"<CN9>", 100, 200);
        assert_eq!(parse_copy_number(&v), CopyNumberClass::Generic);
    }

    #[test]
    fn test_parse_copy_number_generic_cnv() {
        let v = make_cnv(b"<CNV>", 100, 200);
        assert_eq!(parse_copy_number(&v), CopyNumberClass::Generic);
    }

    #[test]
    fn test_cn0_full_containment_transcript_ablation() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CN0>", 24_999_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::TranscriptAblation]
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_cn0_exact_boundaries_transcript_ablation() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CN0>", 25_000_000, 25_006_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::TranscriptAblation]
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_cn0_partial_cds_overlap_feature_truncation() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CN0>", 25_001_500, 25_003_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "CN0 partial CDS overlap should produce feature_truncation, got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "Should have coding_sequence_variant"
        );
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "Should overlap introns too"
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_cn2_full_containment_transcript_amplification() {
        // CN2 is classified as duplication in Perl VEP, so it gets
        // transcript_amplification when fully containing.
        let tx = make_test_transcript();
        let v = make_cnv(b"<CN2>", 24_999_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::TranscriptAmplification]
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_cn3_full_containment_context_term_only() {
        // CN3 is generic copy_number_variation in Perl VEP, so it gets the context
        // term only when fully containing (no transcript_amplification).
        let tx = make_test_transcript();
        let v = make_cnv(b"<CN3>", 24_999_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::CodingTranscriptVariant]
        );
    }

    #[test]
    fn test_cn4_partial_cds_overlap_no_feature_elongation() {
        // CN4 is generic: no feature_elongation (only CN2 gets that).
        let tx = make_test_transcript();
        let v = make_cnv(b"<CN4>", 25_001_500, 25_003_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            !tc.consequences.contains(&Consequence::FeatureElongation),
            "CN4 should not produce feature_elongation, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_generic_cnv_fully_contains_coding_transcript() {
        // Generic CNV (unknown copy number) that fully contains a coding transcript.
        // Perl VEP: coding_transcript_variant as the sole consequence.
        let tx = make_test_transcript();
        let v = make_cnv(b"<CNV>", 24_999_000, 25_007_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::CodingTranscriptVariant],
            "Generic CNV fully containing coding transcript should produce coding_transcript_variant only, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_generic_cnv_partial_cds_and_intron() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CNV>", 25_000_100, 25_001_500);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert_eq!(tc.impact, Impact::MODIFIER);
    }

    #[test]
    fn test_generic_cnv_intron_only() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CNV>", 25_000_500, 25_001_500);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "Intron-only CNV should have intron_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_generic_cnv_5prime_utr() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CNV>", 25_000_010, 25_000_040);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FivePrimeUtrVariant),
            "Should have 5_prime_UTR_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_generic_cnv_3prime_utr() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CNV>", 25_004_500, 25_005_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::ThreePrimeUtrVariant),
            "Should have 3_prime_UTR_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_cnv_upstream() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CNV>", 24_996_000, 24_998_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::UpstreamGeneVariant]
        );
        assert_eq!(tc.impact, Impact::MODIFIER);
        assert_eq!(tc.distance, Some(2000));
    }

    #[test]
    fn test_cnv_downstream() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CNV>", 25_007_000, 25_009_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::DownstreamGeneVariant]
        );
        assert_eq!(tc.impact, Impact::MODIFIER);
        assert_eq!(tc.distance, Some(1000));
    }

    #[test]
    fn test_cnv_beyond_upstream_distance_returns_none() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CNV>", 24_990_000, 24_993_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_none());
    }

    #[test]
    fn test_cnv_non_coding_transcript() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.protein_id = None;

        let v = make_cnv(b"<CNV>", 25_000_050, 25_000_200);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "Non-coding exon overlap should have non_coding_transcript_exon_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_cnv_non_coding_intron_only() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.translation = None;
        tx.protein_id = None;

        let v = make_cnv(b"<CNV>", 25_000_500, 25_001_500);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
    }

    #[test]
    fn test_cnv_consequence_fields() {
        let tx = make_test_transcript();
        let v = make_cnv(b"<CNV>", 25_000_100, 25_001_500);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert_eq!(&*tc.transcript_id, "ENST00000000001");
        assert_eq!(&*tc.gene_id, "ENSG00000000001");
        assert_eq!(tc.gene_symbol, Some("TEST1".into()));
        assert_eq!(tc.strand, 1);
        assert!(tc.canonical);
        assert_eq!(tc.biotype, Some("protein_coding".into()));
        assert_eq!(tc.feature_type, FeatureType::Transcript);
        assert!(tc.cdna_position.is_none());
        assert!(tc.cds_position.is_none());
        assert!(tc.protein_position.is_none());
        assert!(tc.amino_acids.is_none());
        assert!(tc.codons.is_none());
    }

    #[test]
    fn test_cn0_overlapping_stop_codon_has_stop_lost() {
        let tx = make_test_transcript();
        // CDS end at 25_004_299 => stop codon ~25_004_297-25_004_299
        let v = make_cnv(b"<CN0>", 25_004_200, 25_005_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences.contains(&Consequence::StopLost),
            "CN0 overlapping stop codon should have stop_lost, got: {:?}",
            tc.consequences
        );
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
    }

    #[test]
    fn test_cnv_no_sv_end_uses_variant_end() {
        let tx = make_test_transcript();
        let mut v = InputVariant::new(
            "21".into(),
            25_000_100,
            25_001_500,
            b"N".to_vec(),
            b"<CNV>".to_vec(),
        );
        v.variant_class = VariantClass::CopyNumberVariation;
        v.is_structural = true;
        v.sv_end = None;
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
    }

    #[test]
    fn test_consequences_sorted_by_rank() {
        let tx = make_test_transcript();
        // CN0 partial CDS overlap => feature_truncation + coding_sequence_variant + intron_variant
        let v = make_cnv(b"<CN0>", 25_001_500, 25_003_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        // feature_truncation (rank 10) < coding_sequence_variant (rank 23) < intron_variant (rank 28)
        let ranks: Vec<u32> = tc.consequences.iter().map(|c| c.rank()).collect();
        let mut sorted_ranks = ranks.clone();
        sorted_ranks.sort();
        assert_eq!(ranks, sorted_ranks, "Consequences should be sorted by rank");
    }

    #[test]
    fn test_cnv_reverse_strand_upstream() {
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;
        // For reverse strand, upstream is past the genomic end
        let v = make_cnv(b"<CNV>", 25_007_000, 25_009_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::UpstreamGeneVariant]
        );
        assert_eq!(tc.distance, Some(1000));
    }

    #[test]
    fn test_cnv_reverse_strand_downstream() {
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;
        // For reverse strand, downstream is before the genomic start
        let v = make_cnv(b"<CNV>", 24_996_000, 24_998_000);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::DownstreamGeneVariant]
        );
        assert_eq!(tc.distance, Some(2000));
    }
}
