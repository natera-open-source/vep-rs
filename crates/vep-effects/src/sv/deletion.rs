// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Structural deletion (`<DEL>`) consequence calculation.
//!
//! Mirrors Perl VEP's `StructuralVariationOverlapAllele` logic for deletions.
//! Assigns consequence terms based on how the deletion span overlaps the
//! transcript structure (exons, introns, CDS, UTRs).
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`).

use super::{
    bnd_breakend_within_transcript, is_mature_mirna_sv, overlaps_any_exon as sv_overlaps_any_exon,
    overlaps_any_intron_trimmed as sv_overlaps_any_intron_trimmed, overlaps_cds_exon,
    overlaps_five_prime_utr, overlaps_three_prime_utr,
};
use smallvec::SmallVec;
use vep_core::consequence::{
    Consequence, ConsequenceList, FeatureType, Impact, TranscriptConsequence,
};
use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;
use vep_core::variant::{InputVariant, VariantClass};

/// Calculate consequences for a structural deletion against a transcript.
///
/// Structural deletions can produce: transcript_ablation, feature_truncation,
/// coding_sequence_variant, 5_prime_UTR_variant, 3_prime_UTR_variant,
/// intron_variant, stop_lost, non_coding_transcript_exon_variant,
/// upstream_gene_variant, downstream_gene_variant, and others according to
/// overlap extent. `start_lost` needs both deletion ends in exons and no `cds_start_NF`.
pub fn calculate(
    variant: &InputVariant,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    let del_start = variant.start;
    let del_end = variant.sv_end.unwrap_or(variant.end);
    let tx_start = transcript.start;
    let tx_end = transcript.end;

    if del_end < tx_start.saturating_sub(upstream_distance.max(downstream_distance))
        || del_start > tx_end + upstream_distance.max(downstream_distance)
    {
        return None;
    }

    let mut consequences: ConsequenceList = SmallVec::new();
    let mut distance: Option<u64> = None;

    let overlaps_transcript = del_end >= tx_start && del_start <= tx_end;

    if !overlaps_transcript {
        match transcript.strand {
            Strand::Forward => {
                if del_end < tx_start {
                    let dist = tx_start - del_end;
                    if dist <= upstream_distance {
                        consequences.push(Consequence::UpstreamGeneVariant);
                        distance = Some(dist);
                    }
                } else {
                    let dist = del_start - tx_end;
                    if dist <= downstream_distance {
                        consequences.push(Consequence::DownstreamGeneVariant);
                        distance = Some(dist);
                    }
                }
            }
            Strand::Reverse => {
                if del_start > tx_end {
                    let dist = del_start - tx_end;
                    if dist <= upstream_distance {
                        consequences.push(Consequence::UpstreamGeneVariant);
                        distance = Some(dist);
                    }
                } else {
                    let dist = tx_start - del_end;
                    if dist <= downstream_distance {
                        consequences.push(Consequence::DownstreamGeneVariant);
                        distance = Some(dist);
                    }
                }
            }
        }

        if consequences.is_empty() {
            return None;
        }

        return Some(build_consequence(transcript, consequences, distance));
    }

    if del_start <= tx_start && del_end >= tx_end {
        consequences.push(Consequence::TranscriptAblation);
        return Some(build_consequence(transcript, consequences, None));
    }

    // Mature miRNA is Ensembl's tier-2 short-circuit and outranks every term
    // below: `@SORTED_OVERLAP_CONSEQUENCES = sort {$a->tier <=> $b->tier}`
    // (`BaseVariationFeatureOverlapAllele.pm:69`) orders predicates by tier alone
    // and `get_all_OverlapConsequences` (`:266,275-278`) stops at
    // `last if $assigned_tier && $oc->{tier} > $assigned_tier`. `mature_miRNA_variant`
    // is tier 2 (`Utils/Config.pm`, rank 24); `feature_truncation` (rank 10) and
    // every regional term are tier 3, so rank never puts them earlier. Only
    // tier-1 `transcript_ablation`, handled above, can co-occur.
    if is_mature_mirna_sv(transcript, del_start, del_end) {
        consequences.push(Consequence::MatureMirnaVariant);
        return Some(build_consequence(transcript, consequences, None));
    }

    // has_cds() includes NMD transcripts, which have CDS/UTR like protein_coding.
    let is_coding = transcript.has_cds();
    let coding_start = transcript.coding_region_start;
    let coding_end = transcript.coding_region_end;

    // Perl VEP `feature_truncation` predicate for deletions/losses:
    //   within_cdna(@_) AND (partial_overlap_feature(@_) OR complete_within_feature(@_))
    //   AND (deletion(@_) OR copy_number_loss(@_))
    //
    // within_cdna = overlaps exonic sequence; partial_overlap_feature = extends
    // beyond a boundary; complete_within_feature = entirely inside the
    // transcript. transcript_ablation (complete_overlap + deletion) is already
    // excluded above.
    let extends_beyond = del_start < tx_start || del_end > tx_end;
    let complete_within = del_start >= tx_start && del_end <= tx_end;

    // Perl's `feature_truncation` (Utils/VariationEffect.pm:350-373) has two arms.
    //
    // The first is BND-only and bypasses the cDNA requirement outright:
    //   if (chromosome_breakpoint(@_)) { return 1 if within_feature(..., $bvfoa->breakend, 1) }
    // The second requires `within_cdna(@_)`, and `within_cdna` (`:670-693`) returns 1
    // only when the variant maps to at least one cDNA coordinate, which for a transcript
    // means it overlaps an EXON.
    //
    // `chromosome_breakpoint` (`:311-321`) is true only when the structural variant's
    // class SO term is or contains `chromosome_breakpoint`, which is
    // `VariantClass::Translocation` here (the `VariantClass` SO-term mapping). A symbolic
    // `<DEL>`, `<CNV>` or `<DEL:ME:*>` carries `deletion` or `copy_number_loss` instead
    // and therefore takes the second arm, where an intron-only overlap fails
    // `within_cdna` and gets NO `feature_truncation`.
    //
    // `within_cdna`'s second arm also returns `within_transcript` for a variant in
    // a frameshift intron (12 bp or less, `:687-691`); there is no
    // frameshift-intron fallback here, so that narrow case is under-called.
    let is_chromosome_breakpoint = variant.variant_class == VariantClass::Translocation;
    // Which Ensembl allele this tuple represents decides how `feature_truncation` is tested,
    // and the allele string in the output is what distinguishes them. Ensembl builds one
    // allele per entry of `($vf, @{$vf->get_breakends})`, each constructed with
    // `-breakend => $_`, so `within_feature($bvfoa->breakend)` means:
    //
    //   - the `$vf` entry, displayed `N.` or symbolically: "breakend" is the whole variation
    //     feature, so the test is span-versus-transcript, i.e. body overlap.
    //   - a bracket entry, displayed `N[chr:pos[`: "breakend" is that single coordinate, so
    //     the test is whether that coordinate lies inside the transcript.
    //
    let is_paired_breakend = is_chromosome_breakpoint && !variant.is_single_breakend;

    if is_coding {
        if let (Some(cds_start), Some(cds_end)) = (coding_start, coding_end) {
            let overlaps_cds = del_end >= cds_start && del_start <= cds_end;

            let has_exon_overlap = sv_overlaps_any_exon(transcript, del_start, del_end);
            // For a BND the only live arm of Ensembl's `feature_truncation` is
            // `within_feature($bvfoa->breakend)`; its second arm also requires
            // `copy_number_loss or deletion`, which a `chromosome_breakpoint` never
            // satisfies, so body overlap alone must not grant it.
            let truncates = if is_paired_breakend {
                bnd_breakend_within_transcript(variant, transcript)
            } else {
                has_exon_overlap && (extends_beyond || complete_within)
            };

            if truncates {
                consequences.push(Consequence::FeatureTruncation);
                // NO context term here. `coding_transcript_variant` is gated on
                // engulfment (Utils/VariationEffect.pm:491-493: `complete_overlap_feature`
                // at `:169-178` is the variant covering the whole transcript), unlike its
                // siblings `within_nmd_transcript` (`:477-482`) and
                // `within_non_coding_gene` (`:495-500`), which need only
                // `within_transcript`. `feature_truncation` having held means
                // `partial_overlap_feature or complete_within_feature` (`:364-372`),
                // each of which excludes `complete_overlap_feature`.
            }

            // Perl's structural-variant arms of `frameshift`
            // (Utils/VariationEffect.pm:1459-1471) and `inframe_deletion`
            // (:1186-1201) read the deletion's own length once it lies completely
            // within one exon (`complete_within_feature` against each exon, which
            // is the variant inside the exon, not the exon inside the variant).
            // `inframe_deletion` further needs `cds_coords` to be a single
            // Coordinate, i.e. the deletion inside that exon's coding part;
            // `frameshift` needs only the `coding` pre-predicate, an overlap with
            // the CDS. Either one excludes `coding_unknown` (:1535). Both, like
            // `stop_lost` (:1262), hold only for a `deletion` class: a BND routed
            // through here is a `chromosome_breakpoint` and keeps
            // coding_sequence_variant.
            let within_cds_exon = overlaps_cds && overlaps_cds_exon(transcript, del_start, del_end);
            let sv_len = del_end - del_start + 1;
            let within_one_exon = !is_chromosome_breakpoint
                && transcript
                    .exons
                    .iter()
                    .any(|e| del_start >= e.start && del_end <= e.end);
            let within_one_cds_segment = !is_chromosome_breakpoint
                && transcript
                    .exons
                    .iter()
                    .any(|e| del_start >= e.start.max(cds_start) && del_end <= e.end.min(cds_end));
            let frameshift = within_one_exon && within_cds_exon && !sv_len.is_multiple_of(3);
            let inframe_deletion = within_one_cds_segment && sv_len.is_multiple_of(3);
            if frameshift {
                consequences.push(Consequence::FrameshiftVariant);
            }
            if inframe_deletion {
                consequences.push(Consequence::InframeDeletion);
            }
            if within_cds_exon && !frameshift && !inframe_deletion {
                consequences.push(Consequence::CodingSequenceVariant);
            }

            // `stop_lost` (:1262-1274) and `start_lost` (:886-896) for a structural
            // deletion are genomic overlaps with the three codon bases.
            // `start_lost` also passes `_overlaps_start_codon` (:965-990), which
            // needs `cdna_start` and `cdna_end` defined: both ends of the deletion
            // in exons, and no `cds_start_NF`. Perl co-emits
            // `start_retained_variant` there because `_ins_del_start_altered`
            // returns 0 for every structural allele without reading the sequence,
            // so on a structural deletion the retained term is the erroneous
            // member and vep-rs keeps `start_lost` alone: the reverse of the
            // sequence-variant case in `coding::perl_coding_terms`, where the same
            // predicate has certified the start codon intact.
            let (stop_lo, stop_hi, start_lo, start_hi) = match transcript.strand {
                Strand::Forward => (cds_end.saturating_sub(2), cds_end, cds_start, cds_start + 2),
                Strand::Reverse => (cds_start, cds_start + 2, cds_end.saturating_sub(2), cds_end),
            };
            if !is_chromosome_breakpoint && del_end >= stop_lo && del_start <= stop_hi {
                push_unique(&mut consequences, Consequence::StopLost);
            }
            let ends_in_exons = |pos: u64| {
                transcript
                    .exons
                    .iter()
                    .any(|e| pos >= e.start && pos <= e.end)
            };
            if del_end >= start_lo
                && del_start <= start_hi
                && !transcript.facts().cds_start_nf
                && ends_in_exons(del_start)
                && ends_in_exons(del_end)
            {
                push_unique(&mut consequences, Consequence::StartLost);
            }

            if overlaps_five_prime_utr(transcript, del_start, del_end) {
                push_unique(&mut consequences, Consequence::FivePrimeUtrVariant);
            }
            if overlaps_three_prime_utr(transcript, del_start, del_end) {
                push_unique(&mut consequences, Consequence::ThreePrimeUtrVariant);
            }
        } else {
            // Protein-coding but no CDS coords: treat exon overlaps as coding_sequence_variant.
            let has_exon_overlap_noc = sv_overlaps_any_exon(transcript, del_start, del_end);
            let truncates_noc = if is_paired_breakend {
                bnd_breakend_within_transcript(variant, transcript)
            } else {
                has_exon_overlap_noc && (extends_beyond || complete_within)
            };
            if truncates_noc {
                consequences.push(Consequence::FeatureTruncation);
                // No context term, for the reason the CDS-coords arm above documents:
                // `coding_transcript_variant` needs `complete_overlap_feature`, which
                // `feature_truncation` having held already excludes.
            }
            if has_exon_overlap_noc {
                consequences.push(Consequence::CodingSequenceVariant);
            }
        }
    } else {
        let has_exon_overlap_nc = sv_overlaps_any_exon(transcript, del_start, del_end);
        let truncates_nc = if is_paired_breakend {
            bnd_breakend_within_transcript(variant, transcript)
        } else {
            has_exon_overlap_nc && (extends_beyond || complete_within)
        };
        if truncates_nc {
            consequences.push(Consequence::FeatureTruncation);
            if !has_exon_overlap_nc {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
        }
        // No mature-miRNA arm here: the caller returns on `is_mature_mirna_sv` before
        // any tier-3 term is reached, and a second copy of that test could drift.
        if has_exon_overlap_nc {
            push_unique(
                &mut consequences,
                Consequence::NonCodingTranscriptExonVariant,
            );
        }
    }

    // Intron overlap on Ensembl's trimmed boundaries: `_intron_effects` sets
    // `intronic` with `overlap($r_start, $r_end, $intron_start + 2, $intron_end - 2)`,
    // attributing the four invariant donor/acceptor bases to the splice terms.
    if sv_overlaps_any_intron_trimmed(transcript, del_start, del_end) {
        push_unique(&mut consequences, Consequence::IntronVariant);
    }

    // No splice_polypyrimidine_tract_variant for structural deletions: Perl adds
    // PPT only for small-variant types, and a large deletion spanning an intron
    // always overlaps the PPT window.

    // Context terms for non-coding transcripts (excluding NMD, handled separately).
    // Perl's within_non_coding_gene: not translation and not mature_miRNA and not non_coding_exon_variant
    if !is_coding
        && !transcript.is_nmd_transcript()
        && !consequences.contains(&Consequence::NonCodingTranscriptExonVariant)
        && !consequences.contains(&Consequence::MatureMirnaVariant)
    {
        push_unique(&mut consequences, Consequence::NonCodingTranscriptVariant);
    }

    // NMD transcript context: Perl adds NMD_transcript_variant alongside other
    // consequences (not just as a fallback), e.g., "feature_truncation,coding_sequence_variant,...,NMD_transcript_variant"
    if transcript.is_nmd_transcript() {
        push_unique(&mut consequences, Consequence::NmdTranscriptVariant);
    }

    if consequences.is_empty() {
        if is_coding {
            consequences.push(Consequence::CodingTranscriptVariant);
        } else {
            consequences.push(Consequence::NonCodingTranscriptVariant);
        }
    }

    consequences.sort();
    consequences.dedup();

    Some(build_consequence(transcript, consequences, distance))
}

/// Push a consequence only if not already present.
fn push_unique(consequences: &mut ConsequenceList, consequence: Consequence) {
    if !consequences.contains(&consequence) {
        consequences.push(consequence);
    }
}

/// Build the `TranscriptConsequence` from accumulated data, with the impact of
/// the most severe consequence.
fn build_consequence(
    transcript: &Transcript,
    consequences: ConsequenceList,
    distance: Option<u64>,
) -> TranscriptConsequence {
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

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::variant::{InputVariant, VariantClass};

    /// Helper to build a structural deletion variant.
    fn make_del(chr: &str, start: u64, end: u64) -> InputVariant {
        let mut v = InputVariant::new(chr.into(), start, end, b"N".to_vec(), b"-".to_vec());
        v.variant_class = VariantClass::StructuralDeletion;
        v.is_structural = true;
        v.sv_end = Some(end);
        v
    }

    fn tx() -> Transcript {
        crate::test_helpers::make_test_transcript()
    }

    // Transcript layout reminder (forward strand, chr21):
    //   Exon 1:   25_000_000 - 25_000_299  (CDS starts at 25_000_050)
    //   Intron 1: 25_000_300 - 25_001_999
    //   Exon 2:   25_002_000 - 25_002_299
    //   Intron 2: 25_002_300 - 25_003_999
    //   Exon 3:   25_004_000 - 25_006_000  (CDS ends at 25_004_299)
    //   CDS:      25_000_050 - 25_004_299
    //   5' UTR:   25_000_000 - 25_000_049
    //   3' UTR:   25_004_300 - 25_006_000

    #[test]
    fn test_no_overlap_returns_none() {
        let v = make_del("21", 24_990_000, 24_994_000);
        assert!(calculate(&v, &tx(), 5000, 5000).is_none());
    }

    #[test]
    fn test_transcript_ablation() {
        let v = make_del("21", 24_999_000, 25_007_000);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert_eq!(
            result.consequences.to_vec(),
            vec![Consequence::TranscriptAblation]
        );
        assert_eq!(result.impact, Impact::HIGH);
        assert_eq!(
            format!("{:.2}", result.feature_overlap(v.start, v.end).unwrap().1),
            "100.00"
        );
    }

    #[test]
    fn test_upstream_gene_variant_forward() {
        let v = make_del("21", 24_998_000, 24_999_500);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert_eq!(
            result.consequences.to_vec(),
            vec![Consequence::UpstreamGeneVariant]
        );
        assert_eq!(result.impact, Impact::MODIFIER);
        assert_eq!(result.distance, Some(500));
    }

    #[test]
    fn test_downstream_gene_variant_forward() {
        let v = make_del("21", 25_006_100, 25_006_500);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert_eq!(
            result.consequences.to_vec(),
            vec![Consequence::DownstreamGeneVariant]
        );
        assert_eq!(result.impact, Impact::MODIFIER);
        assert_eq!(result.distance, Some(100));
    }

    #[test]
    fn test_intron_only_structural_within() {
        // Structural deletion entirely within intron 1 (25_000_300 - 25_001_999):
        // no exon overlap, so no feature_truncation.
        let v = make_del("21", 25_000_500, 25_001_500);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(result.consequences.contains(&Consequence::IntronVariant));
        assert!(
            !result
                .consequences
                .contains(&Consequence::FeatureTruncation),
            "Structural intron-only DEL within transcript should not get feature_truncation, got: {:?}",
            result.consequences
        );
        assert!(!result
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
    }

    #[test]
    fn test_intron_only_symbolic_del_extending_beyond_gets_no_truncation() {
        use vep_core::transcript::Exon;

        // A symbolic `<DEL>` reaching past a transcript edge with an intron-only
        // overlap gets no `feature_truncation`: Perl's predicate requires
        // `within_cdna` outside its BND arm (Utils/VariationEffect.pm:361), which
        // needs an exon overlap, and `<DEL>`'s class SO term is `deletion`, not
        // `chromosome_breakpoint`. The fixture's exon 1 does not start at tx_start,
        // leaving an intron-only gap at the start edge.
        let mut t = tx();
        t.exons = vec![
            Exon {
                stable_id: Some("ENSE00000000001".into()),
                start: 25_000_500,
                end: 25_000_799,
                rank: 1,
                phase: -1,
                end_phase: 0,
            },
            Exon {
                stable_id: Some("ENSE00000000002".into()),
                start: 25_002_000,
                end: 25_002_299,
                rank: 2,
                phase: 0,
                end_phase: 0,
            },
            Exon {
                stable_id: Some("ENSE00000000003".into()),
                start: 25_004_000,
                end: 25_006_000,
                rank: 3,
                phase: 0,
                end_phase: -1,
            },
        ];

        // Structural deletion that starts before transcript and ends in the
        // intron-only gap at the start of the transcript body (before exon 1).
        // Del: 24_999_000 - 25_000_400 extends before tx_start, and the overlap is
        // only in the 25_000_000-25_000_400 region which has no exon.
        let v = make_del("21", 24_999_000, 25_000_400);
        let result = calculate(&v, &t, 5000, 5000).unwrap();
        assert!(
            !result
                .consequences
                .contains(&Consequence::FeatureTruncation),
            "a symbolic DEL with no exon overlap fails Perl's within_cdna, whichever \
             transcript edge it reaches past: {:?}",
            result.consequences
        );
        assert!(
            !result
                .consequences
                .contains(&Consequence::CodingSequenceVariant),
            "No exon overlap means no coding_sequence_variant, got: {:?}",
            result.consequences
        );
    }

    #[test]
    fn test_intron_only_breakend_extending_beyond_still_gets_truncation() {
        use vep_core::transcript::Exon;

        // The complement: Perl's `feature_truncation` opens with a BND-only early
        // return that bypasses `within_cdna`:
        //   if (chromosome_breakpoint(@_)) { return 1 if within_feature(..., breakend, 1) }
        // and `chromosome_breakpoint` (`:311-321`) is the class SO term of
        // `VariantClass::Translocation`, so the same intron-only geometry the test
        // above denies a truncation to must receive one here.
        let mut t = tx();
        t.exons = vec![
            Exon {
                stable_id: Some("ENSE00000000001".into()),
                start: 25_000_500,
                end: 25_000_799,
                rank: 1,
                phase: -1,
                end_phase: 0,
            },
            Exon {
                stable_id: Some("ENSE00000000002".into()),
                start: 25_002_000,
                end: 25_002_299,
                rank: 2,
                phase: 0,
                end_phase: 0,
            },
            Exon {
                stable_id: Some("ENSE00000000003".into()),
                start: 25_004_000,
                end: 25_006_000,
                rank: 3,
                phase: 0,
                end_phase: -1,
            },
        ];

        let mut v = make_del("21", 24_999_000, 25_000_400);
        v.variant_class = VariantClass::Translocation;
        let result = calculate(&v, &t, 5000, 5000).unwrap();
        assert!(
            result
                .consequences
                .contains(&Consequence::FeatureTruncation),
            "a chromosome breakpoint takes Perl's BND arm and needs no exon overlap: {:?}",
            result.consequences
        );
        // The truncation must arrive without a coding context term: Ensembl gates
        // `coding_transcript_variant` on `complete_overlap_feature`
        // (Utils/VariationEffect.pm:491-493), which a partial overlap never satisfies.
        assert!(
            !result
                .consequences
                .contains(&Consequence::CodingTranscriptVariant),
            "coding_transcript_variant needs complete_overlap_feature, which a truncating \
             partial overlap excludes: {:?}",
            result.consequences
        );
    }

    #[test]
    fn test_intron_only_explicit_sequence() {
        // An explicit-sequence deletion entirely within intron 1 has no exon overlap,
        // so no feature_truncation.
        let mut v = InputVariant::new(
            "21".into(),
            25_000_500,
            25_001_500,
            b"NNNNN".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::StructuralDeletion;
        v.is_structural = false; // explicit-sequence deletion
        v.sv_end = Some(25_001_500);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(result.consequences.contains(&Consequence::IntronVariant));
        assert!(
            !result
                .consequences
                .contains(&Consequence::FeatureTruncation),
            "Explicit-sequence intron-only DEL should not get feature_truncation, got: {:?}",
            result.consequences
        );
    }

    #[test]
    fn test_coding_sequence_variant() {
        // 101 bp inside exon 2's CDS: Perl's SV `frameshift` arm takes the term.
        let v = make_del("21", 25_002_050, 25_002_150);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(result
            .consequences
            .contains(&Consequence::FrameshiftVariant));
        assert!(!result
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        // Across intron 2 into exon 3 the deletion is inside no single exon, and
        // `coding_unknown` holds.
        let v = make_del("21", 25_002_050, 25_004_050);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(result
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
    }

    #[test]
    fn test_five_prime_utr_variant() {
        let v = make_del("21", 25_000_010, 25_000_040);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(result
            .consequences
            .contains(&Consequence::FivePrimeUtrVariant));
        assert!(!result
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
    }

    #[test]
    fn test_three_prime_utr_variant() {
        let v = make_del("21", 25_004_500, 25_004_800);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(result
            .consequences
            .contains(&Consequence::ThreePrimeUtrVariant));
        assert!(!result
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
    }

    #[test]
    fn test_feature_truncation_with_coding() {
        // Del: 25_002_050 - 25_007_000 (overlaps exon 2 CDS + intron 2 + exon 3 + beyond tx end)
        let v = make_del("21", 25_002_050, 25_007_000);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(result
            .consequences
            .contains(&Consequence::FeatureTruncation));
        assert!(result
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert_eq!(result.impact, Impact::HIGH);
    }

    #[test]
    fn test_stop_lost_forward_strand() {
        // Deletion covers the CDS end (stop codon region) at 25_004_299.
        let v = make_del("21", 25_004_200, 25_004_400);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(
            result.consequences.contains(&Consequence::StopLost),
            "Expected stop_lost, got: {:?}",
            result.consequences
        );
    }

    #[test]
    fn test_start_lost_for_structural_deletion_with_both_ends_exonic() {
        // Perl's SV arm of `start_lost` is a genomic overlap with the start codon
        // once `_overlaps_start_codon` holds (both deletion ends in exons). The
        // 21 bp span from 5' UTR into the CDS is not inside the coding part of the
        // exon, so it is neither an inframe_deletion nor a frameshift and keeps
        // coding_sequence_variant.
        let v = make_del("21", 25_000_040, 25_000_060);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        crate::test_helpers::assert_consequence_set_eq(
            &result,
            &[
                "5_prime_UTR_variant",
                "coding_sequence_variant",
                "feature_truncation",
                "start_lost",
            ],
        );
    }

    #[test]
    fn test_start_lost_withheld_when_a_deletion_end_is_intronic() {
        // Starting in intron 1 leaves Perl's `cdna_start` undefined, so
        // `_overlaps_start_codon` is 0 even though the span covers the ATG.
        let v = make_del("21", 25_000_100, 25_002_010);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(
            !result.consequences.contains(&Consequence::StartLost),
            "got {:?}",
            result.consequences
        );
        let v = make_del("21", 24_999_990, 25_000_060);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(
            !result.consequences.contains(&Consequence::StartLost),
            "got {:?}",
            result.consequences
        );
    }

    #[test]
    fn test_frameshift_for_structural_deletion_within_one_cds_exon() {
        // 68 bp inside exon 2 (all CDS): frameshift_variant replaces
        // coding_sequence_variant, as for the 1000 Genomes chr21 deletions Perl
        // calls `feature_truncation,frameshift_variant`.
        let v = make_del("21", 25_002_010, 25_002_077);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        crate::test_helpers::assert_consequence_set_eq(
            &result,
            &["feature_truncation", "frameshift_variant"],
        );
    }

    #[test]
    fn test_breakend_within_one_cds_exon_keeps_coding_sequence_variant() {
        // A `chromosome_breakpoint` is neither `deletion` nor `copy_number_loss`,
        // so Perl's SV frameshift, inframe_deletion and stop_lost arms are 0 and
        // `coding_unknown` holds.
        let mut v = make_del("21", 25_004_230, 25_004_297);
        v.variant_class = VariantClass::Translocation;
        v.is_single_breakend = true;
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        crate::test_helpers::assert_consequence_set_eq(
            &result,
            &["coding_sequence_variant", "feature_truncation"],
        );
    }

    #[test]
    fn test_inframe_deletion_for_structural_deletion_within_one_cds_exon() {
        let v = make_del("21", 25_002_010, 25_002_075);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        crate::test_helpers::assert_consequence_set_eq(
            &result,
            &["feature_truncation", "inframe_deletion"],
        );
    }

    #[test]
    fn test_structural_deletion_spanning_an_intron_keeps_coding_sequence_variant() {
        // 66 bp from exon 2 into intron 2 is inside no single exon: neither
        // frameshift nor inframe_deletion, and Perl's `coding_unknown` holds.
        let v = make_del("21", 25_002_250, 25_002_315);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        crate::test_helpers::assert_consequence_set_eq(
            &result,
            &[
                "coding_sequence_variant",
                "feature_truncation",
                "intron_variant",
            ],
        );
    }

    #[test]
    fn test_stop_lost_for_structural_deletion_ending_inside_the_stop_codon() {
        // Perl tests `overlap(cre - 2, cre, start, end)`: a deletion that ends on
        // the first stop-codon base still loses the stop. 68 bp inside exon 3.
        let v = make_del("21", 25_004_230, 25_004_297);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        crate::test_helpers::assert_consequence_set_eq(
            &result,
            &["feature_truncation", "frameshift_variant", "stop_lost"],
        );
    }

    #[test]
    fn test_overlap_bp_and_pc() {
        // Deletion partially overlaps transcript.
        // Del: 25_003_000 - 25_007_000 (overlaps 25_003_000-25_006_000 = 3001bp of tx)
        // Tx len: 25_006_000 - 25_000_000 + 1 = 6001
        // OverlapBP = 3001, OverlapPC = 3001/6001*100 = 50.01
        let v = make_del("21", 25_003_000, 25_007_000);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert_eq!(result.feature_overlap(v.start, v.end).unwrap().0, 3001);
        assert_eq!(
            format!("{:.2}", result.feature_overlap(v.start, v.end).unwrap().1),
            "50.01"
        );
    }

    #[test]
    fn test_multiple_consequences_sorted_by_rank() {
        // Del: 25_000_030 - 25_004_500 (covers 5'UTR + CDS + introns + 3'UTR)
        let v = make_del("21", 25_000_030, 25_004_500);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();

        assert!(result.consequences.len() >= 3);
        for pair in result.consequences.windows(2) {
            assert!(
                pair[0].rank() <= pair[1].rank(),
                "Consequences not sorted: {:?} before {:?}",
                pair[0],
                pair[1]
            );
        }
    }

    #[test]
    fn test_non_coding_transcript() {
        let mut t = tx();
        t.biotype = "lncRNA".into();
        t.coding_region_start = None;
        t.coding_region_end = None;
        t.translation = None;

        let v = make_del("21", 25_000_050, 25_000_150);
        let result = calculate(&v, &t, 5000, 5000).unwrap();
        assert!(result
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
        assert!(!result
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
    }

    #[test]
    fn test_non_coding_feature_truncation() {
        let mut t = tx();
        t.biotype = "lncRNA".into();
        t.coding_region_start = None;
        t.coding_region_end = None;
        t.translation = None;

        let v = make_del("21", 24_999_500, 25_000_150);
        let result = calculate(&v, &t, 5000, 5000).unwrap();
        assert!(result
            .consequences
            .contains(&Consequence::FeatureTruncation));
        assert!(result
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
    }

    #[test]
    fn test_too_far_upstream_returns_none() {
        let v = make_del("21", 24_990_000, 24_994_000);
        assert!(calculate(&v, &tx(), 5000, 5000).is_none());
    }

    #[test]
    fn test_too_far_downstream_returns_none() {
        let v = make_del("21", 25_012_000, 25_015_000);
        assert!(calculate(&v, &tx(), 5000, 5000).is_none());
    }

    #[test]
    fn test_transcript_metadata_populated() {
        let v = make_del("21", 24_999_000, 25_007_000);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert_eq!(&*result.transcript_id, "ENST00000000001");
        assert_eq!(&*result.gene_id, "ENSG00000000001");
        assert_eq!(result.gene_symbol, Some("TEST1".into()));
        assert_eq!(result.gene_symbol_source, Some("HGNC".into()));
        assert_eq!(result.hgnc_id, Some("HGNC:0001".into()));
        assert_eq!(result.biotype, Some("protein_coding".into()));
        assert!(result.canonical);
        assert_eq!(result.strand, 1);
        assert_eq!(result.feature_type, FeatureType::Transcript);
    }

    #[test]
    fn test_different_chromosome_is_not_checked_here() {
        let v = make_del("1", 25_000_000, 25_006_000);
        // The function does not check chromosome (the caller does, via
        // TranscriptBinIndex), so positional overlap still yields a result.
        let result = calculate(&v, &tx(), 5000, 5000);
        assert!(result.is_some());
    }

    #[test]
    fn test_exact_transcript_boundary_ablation() {
        let v = make_del("21", 25_000_000, 25_006_000);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert_eq!(
            result.consequences.to_vec(),
            vec![Consequence::TranscriptAblation]
        );
        assert_eq!(
            format!("{:.2}", result.feature_overlap(v.start, v.end).unwrap().1),
            "100.00"
        );
    }

    #[test]
    fn test_spanning_exon_and_intron() {
        let v = make_del("21", 25_000_200, 25_000_500);
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        assert!(result
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert!(result.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_sv_end_fallback_to_variant_end() {
        // When sv_end is None, should use variant.end.
        let mut v = InputVariant::new(
            "21".into(),
            25_000_050,
            25_000_150,
            b"N".to_vec(),
            b"-".to_vec(),
        );
        v.variant_class = VariantClass::StructuralDeletion;
        v.is_structural = true;
        v.sv_end = None; // No sv_end: variant.end = 25_000_150 is used
        let result = calculate(&v, &tx(), 5000, 5000).unwrap();
        // 101 bp from the ATG inside exon 1's CDS: frameshift_variant (with
        // start_lost) rather than coding_sequence_variant; a 5 kb sv_end would
        // have reached exon 3 and read differently.
        assert!(result
            .consequences
            .contains(&Consequence::FrameshiftVariant));
        assert!(result.consequences.contains(&Consequence::StartLost));
        assert!(!result
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
    }
}
