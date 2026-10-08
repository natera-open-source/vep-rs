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
//! 2. **Point insertions** (a symbolic record with neither SVLEN nor END, which
//!    ensembl-io reads as the single base at POS + 1): `variant.start == variant.end`.
//!    OverlapBP = 1.
//!
//! 3. **Insertion pairs** (a record whose END is its position or whose SVLEN is 0):
//!    `variant.start == variant.end + 1`, the insertion between the two bases, as
//!    Ensembl VEP carries it. Its one differing region is the pair itself
//!    (`StructuralVariationOverlapAllele::_get_differing_regions`), so every predicate
//!    reads it as `_intron_effects` reads a sequence insertion; `calculate_insertion_pair`
//!    below is that reading, and the row carries no OverlapBP (overlap length 0).
//!
//! ## Consequence assignment
//!
//! For each transcript region that the SV range overlaps:
//! - `feature_elongation`: the insertion lies entirely inside the transcript and
//!   overlaps exonic (cDNA) sequence, coding or non-coding (Perl:
//!   `within_cdna AND complete_within_feature AND insertion`); HIGH impact,
//!   because the transcript gets longer.
//! - `coding_sequence_variant`: Overlap with any CDS exon.
//! - `start_lost`: the span covers a start-codon base with both ends in exons
//!   (`structural_start_lost`).
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
    engulfed_transcript_row_terms, is_mature_mirna_sv, overlaps_any_exon,
    overlaps_any_intron_trimmed, overlaps_cds_exon, overlaps_five_prime_utr,
    overlaps_polypyrimidine_tract, overlaps_three_prime_utr,
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

    if sv_start == sv_end + 1 {
        return calculate_insertion_pair(
            variant,
            transcript,
            upstream_distance,
            downstream_distance,
        );
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
        protein_id: transcript.protein_id.clone(),
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

        tc.consequences = engulfed_transcript_row_terms(transcript, sv_start, sv_end);
        tc.impact = tc
            .consequences
            .first()
            .map(|c| c.impact())
            .unwrap_or(Impact::MODIFIER);
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

        if structural_start_lost(transcript, sv_start, sv_end) {
            consequences.push(Consequence::StartLost);
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
            add_polypyrimidine_tract(&mut consequences, variant, transcript, sv_start, sv_end);
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
            add_polypyrimidine_tract(&mut consequences, variant, transcript, sv_start, sv_end);
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

/// Perl's reading of a symbolic insertion carried as the pair `start = POS + 1, end = POS`.
///
/// `overlap(start, end, a, b)` on the pair holds only when both flanks lie inside
/// `[a, b]`, so `within_feature`, `complete_within_feature`, the structural `start_lost`
/// arm and the intronic test of `_intron_effects` (`BaseTranscriptVariationAllele.pm`,
/// with its two insertion special cases, `start == intron_start + 2` and
/// `end == intron_end - 2`) all read both flanks. The `exon`, `intron` and `utr`
/// pre-predicates of `_bvfo_preds` read the sorted flank pair, the exon test stretched
/// by 12 bases on a transcript with a frameshift intron. `within_cdna` holds whenever a
/// flank is exonic: `Mapper::map_insert` turns the flank's coordinate into the insert
/// coordinate and drops the gap of an intronic flank; inside a frameshift intron it
/// holds through `within_transcript`. `within_cds` reads that insert coordinate, which
/// is `(c + 1, c)` after an exonic 5' flank at CDS position `c` and `(c, c - 1)` before
/// an exonic 3' flank at `c`, and needs `end > 0` and `start <= length`: the insertion
/// right before the coding region's first base or right after its last is not within
/// the CDS and takes the UTR term through `_before_coding` / `_after_coding`'s insertion
/// special cases instead. `coding_unknown` is `within_cds` for a structural insertion
/// (`inframe_insertion` returns 0 for every structural allele and `frameshift` needs a
/// deletion), `feature_elongation` is `within_cdna and complete_within_feature`, and an
/// empty term list takes Perl's default, `intergenic_variant`.
fn calculate_insertion_pair(
    variant: &InputVariant,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    let s = variant.start;
    let e = variant.sv_end.unwrap_or(variant.end);
    debug_assert_eq!(s, e + 1);
    let (lo, hi) = (e, s);
    let tx_start = transcript.start;
    let tx_end = transcript.end;

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
        protein_id: transcript.protein_id.clone(),
        swissprot: transcript.swissprot.clone(),
        trembl: transcript.trembl.clone(),
        refseq: transcript.refseq.clone(),
        ..TranscriptConsequence::default()
    };

    // `within_feature`: both flanks inside the transcript. Otherwise `_before_start`
    // reads `end` and `_after_end` reads `start`, which is what
    // `classify_upstream_downstream` does with the pair.
    let within_feature = e >= tx_start && s <= tx_end;
    if !within_feature {
        let csq = classify_upstream_downstream(s, e, transcript, up_boundary, down_boundary);
        if let Some((consequence, distance)) = csq {
            tc.consequences = smallvec![consequence];
            tc.impact = consequence.impact();
            tc.distance = Some(distance);
            return Some(tc);
        }
        return None;
    }

    let facts = transcript.facts();
    let stretch = if facts.vefc_has_frameshift_intron {
        12
    } else {
        0
    };
    let exon_pred = overlaps_any_exon(
        transcript,
        lo.saturating_sub(stretch),
        hi.saturating_add(stretch),
    );
    let flank_exonic = overlaps_any_exon(transcript, lo, hi);
    let introns = super::transcript_introns(transcript);
    let in_frameshift_intron = introns
        .iter()
        .any(|i| i.end.saturating_sub(i.start) <= 12 && e >= i.start && s <= i.end);
    let within_cdna = flank_exonic || in_frameshift_intron;
    let intronic = introns.iter().any(|i| {
        let frameshift = i.end.saturating_sub(i.start) <= 12;
        if frameshift && e >= i.start && s <= i.end {
            return false;
        }
        (e >= i.start + 2 && s + 2 <= i.end) || s == i.start + 2 || e + 2 == i.end
    });

    let mut consequences: ConsequenceList = SmallVec::new();

    if transcript.has_cds() {
        if let Some((cds_lo, cds_hi)) = crate::consequences::genomic_coding_bounds(transcript) {
            // `_bvfo_preds`: `coding` for a sorted pair overlapping the coding region on an
            // exon (`cds_coords` of such a pair is never empty); `utr` where the sorted pair
            // reaches past either end of the coding region.
            let coding_pred = hi >= cds_lo && lo <= cds_hi && exon_pred;
            let utr_pred = exon_pred && (lo < cds_lo || hi > cds_hi);

            if coding_pred && insert_within_cds(transcript, s, e) {
                consequences.push(Consequence::CodingSequenceVariant);
            }
            if coding_pred && structural_start_lost(transcript, s, e) {
                consequences.push(Consequence::StartLost);
            }
            if utr_pred && within_cdna {
                let before_coding = s == cds_lo || (e >= tx_start && s < cds_lo);
                let after_coding = e == cds_hi || (e > cds_hi && s <= tx_end);
                let (five_prime, three_prime) = match transcript.strand {
                    Strand::Forward => (before_coding, after_coding),
                    Strand::Reverse => (after_coding, before_coding),
                };
                if five_prime {
                    consequences.push(Consequence::FivePrimeUtrVariant);
                }
                if three_prime {
                    consequences.push(Consequence::ThreePrimeUtrVariant);
                }
            }
        }
        if within_cdna {
            consequences.push(Consequence::FeatureElongation);
        }
    } else {
        if within_cdna {
            consequences.push(Consequence::FeatureElongation);
        }
        let mirna = mature_mirna_contains_pair(transcript, s, e);
        if mirna {
            consequences.push(Consequence::MatureMirnaVariant);
        }
        // `non_coding_exon_variant`: both flanks inside one exon of the overlapped list.
        let exon_variant = !mirna && transcript.exons.iter().any(|x| e >= x.start && s <= x.end);
        if exon_variant {
            consequences.push(Consequence::NonCodingTranscriptExonVariant);
        } else if !mirna && !transcript.is_nmd_transcript() {
            consequences.push(Consequence::NonCodingTranscriptVariant);
        }
    }

    if intronic {
        consequences.push(Consequence::IntronVariant);
    }
    if overlaps_polypyrimidine_tract(transcript, lo, hi) {
        consequences.push(Consequence::SplicePolypyrimidineTractVariant);
    }
    if transcript.is_nmd_transcript() {
        consequences.push(Consequence::NmdTranscriptVariant);
    }
    if consequences.is_empty() {
        consequences.push(Consequence::IntergenicVariant);
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

/// `within_cds` on the insert coordinate of a pair: `(c + 1, c)` after an exonic 5'
/// flank at CDS position `c`, `(c, c - 1)` before an exonic 3' flank at `c`, each
/// needing `end > 0` and `start <= length`; one coordinate `(c + 1, c)` when both flanks
/// share an exon. The 5' and 3' flanks follow the transcript's strand.
fn insert_within_cds(transcript: &Transcript, s: u64, e: u64) -> bool {
    let cds_len = transcript
        .vefc
        .as_ref()
        .and_then(|v| v.translateable_seq.as_ref())
        .map(|t| t.len() as u64)
        .or_else(|| {
            Some(
                transcript
                    .cdna_coding_end?
                    .checked_sub(transcript.cdna_coding_start?)?
                    + 1,
            )
        })
        .unwrap_or(0);
    if cds_len == 0 {
        return false;
    }
    let cds_pos = |pos: u64| match crate::mapper::map_genomic_to_transcript(pos, transcript, 0, 0) {
        Some(crate::mapper::TranscriptPosition::Coding { cds_pos, .. }) => Some(cds_pos),
        _ => None,
    };
    let (five_flank, three_flank) = match transcript.strand {
        Strand::Forward => (e, s),
        Strand::Reverse => (s, e),
    };
    match (cds_pos(five_flank), cds_pos(three_flank)) {
        (Some(_), Some(_)) => true,
        (Some(c), None) => c < cds_len,
        (None, Some(c)) => c >= 2,
        (None, None) => false,
    }
}

/// `within_mature_miRNA` on a pair: both flanks inside a genomic range of the
/// transcript's `miRNA` attribute.
fn mature_mirna_contains_pair(transcript: &Transcript, s: u64, e: u64) -> bool {
    super::mature_mirna_genomic_ranges(transcript)
        .iter()
        .any(|&(lo, hi)| e >= lo && s <= hi)
}

/// Add `splice_polypyrimidine_tract_variant` for an insertion that overlaps an intron.
///
/// Ensembl reads a `<CNV:TR>` as a `VariationFeature` carrying the literal allele
/// string (ensembl-vep `Parser/VCF.pm:530-543`), so the structural rule in
/// `overlaps_polypyrimidine_tract` is not its rule; a tandem repeat is tested by its
/// endpoints. Every other insertion, point or ranged, is a
/// `StructuralVariationFeature` and takes the structural rule over its whole span.
fn add_polypyrimidine_tract(
    consequences: &mut ConsequenceList,
    variant: &InputVariant,
    transcript: &Transcript,
    sv_start: u64,
    sv_end: u64,
) {
    if variant.variant_class == VariantClass::TandemRepeat {
        add_ranged_endpoint_splice_consequences(consequences, transcript, sv_start, sv_end);
    } else if overlaps_polypyrimidine_tract(transcript, sv_start, sv_end)
        && !consequences.contains(&Consequence::SplicePolypyrimidineTractVariant)
    {
        consequences.push(Consequence::SplicePolypyrimidineTractVariant);
    }
}

/// Add splice_polypyrimidine_tract_variant for a tandem repeat when an
/// endpoint falls in the polypyrimidine tract.
///
/// Each endpoint is checked individually against every intron (a point repeat
/// has one). Fires when an endpoint, not the full span, falls 2-16bp from the
/// intron acceptor (strand-aware).
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

/// Perl's structural-variant arm of `start_lost` (Utils/VariationEffect.pm:886-896):
/// a genomic overlap between the allele's span and the three start-codon bases,
/// `coding_region_start..+2` on the forward strand and `coding_region_end-2..` on
/// the reverse, for every structural class. `_overlaps_start_codon` (:965-990) runs
/// first and needs `cdna_start` and `cdna_end` defined, which `cdna_start_unshifted`
/// (BaseTranscriptVariation.pm:194-208) leaves undefined when either end of the span
/// maps to a Gap, so both span ends must lie in exons; it also returns 0 on a
/// `cds_start_NF` transcript. `start_retained_variant` co-fires in Perl on every such
/// allele because `_ins_del_start_altered` returns 0 for a
/// `TranscriptStructuralVariationAllele` without reading sequence; that term is not
/// emitted here.
pub(super) fn structural_start_lost(transcript: &Transcript, sv_start: u64, sv_end: u64) -> bool {
    let (Some(cds_start), Some(cds_end)) =
        (transcript.coding_region_start, transcript.coding_region_end)
    else {
        return false;
    };
    if transcript.facts().cds_start_nf {
        return false;
    }
    let (start_lo, start_hi) = match transcript.strand {
        Strand::Forward => (cds_start, cds_start + 2),
        Strand::Reverse => (cds_end.saturating_sub(2), cds_end),
    };
    let in_exon = |pos: u64| {
        transcript
            .exons
            .iter()
            .any(|e| pos >= e.start && pos <= e.end)
    };
    sv_end >= start_lo && sv_start <= start_hi && in_exon(sv_start) && in_exon(sv_end)
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
        // The span covers the start codon at 25_000_050-25_000_052 with both ends in
        // exon 1, so Perl's structural `start_lost` arm holds; its co-fired
        // `start_retained_variant` is not emitted.
        assert!(tc.consequences.contains(&Consequence::StartLost));
        assert!(!tc.consequences.contains(&Consequence::StartRetainedVariant));
        assert_eq!(tc.impact, Impact::HIGH);
    }

    /// `structural_start_lost` on the reverse strand reads the start codon from
    /// `coding_region_end`: 25_004_297-25_004_299 inside exon 3 (25_004_000-25_006_000).
    #[test]
    fn test_structural_start_lost_reverse_strand_reads_coding_region_end() {
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;
        assert!(structural_start_lost(&tx, 25_004_250, 25_004_350));
        assert!(structural_start_lost(&tx, 25_004_299, 25_004_299));
        assert!(!structural_start_lost(&tx, 25_004_250, 25_004_296));
        assert!(!structural_start_lost(&tx, 25_004_300, 25_004_400));
        // The forward-strand codon (25_000_050-25_000_052) is not a start codon here.
        assert!(!structural_start_lost(&tx, 25_000_040, 25_000_100));
    }

    /// Both span ends must map to exons: a span from exon 1 into intron 1 covers the
    /// forward-strand start codon but has no `cdna_end`, so `_overlaps_start_codon`
    /// returns 0.
    #[test]
    fn test_structural_start_lost_needs_both_ends_in_exons() {
        let tx = make_test_transcript();
        assert!(structural_start_lost(&tx, 25_000_040, 25_000_100));
        assert!(!structural_start_lost(&tx, 25_000_040, 25_001_000));
        assert!(!structural_start_lost(&tx, 24_999_990, 25_000_100));
        // Exon 1 to exon 2: both ends exonic, the codon inside the span.
        assert!(structural_start_lost(&tx, 25_000_040, 25_002_100));
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
    fn test_tandem_repeat_no_polypyrimidine_endpoints_outside_tract() {
        let tx = make_test_transcript();
        // Tandem repeat (sv_start != sv_end) near the acceptor but both endpoints
        // outside the 2-16bp polypyrimidine window.
        // Intron 1: 25_000_300 - 25_001_999, acceptor at 25_001_999 (forward strand).
        // sv_start=25_001_980: dist_to_acceptor = 19 (>16, outside)
        // sv_end=25_001_998: dist_to_acceptor = 1 (<2, outside)
        let mut v = make_ins(25_001_980, 25_001_998);
        v.variant_class = vep_core::variant::VariantClass::TandemRepeat;
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            !tc.consequences
                .contains(&Consequence::SplicePolypyrimidineTractVariant),
            "Tandem repeat with both endpoints outside tract should not get polypyrimidine, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_ranged_insertion_span_over_tract_is_polypyrimidine() {
        let tx = make_test_transcript();
        // The same span as a structural insertion: Perl's `_intron_effects` tests
        // the whole span, `overlap(25_001_980, 25_001_998, 25_001_983, 25_001_997)`,
        // and the span touches no exon, so the term holds.
        let v = make_ins(25_001_980, 25_001_998);
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_some());
        let tc = result.unwrap();
        assert!(
            tc.consequences
                .contains(&Consequence::SplicePolypyrimidineTractVariant),
            "Structural insertion whose span covers the tract should get polypyrimidine, got: {:?}",
            tc.consequences
        );
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
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
        // Ranged insertion spanning the entire intron and the flanking exon ends:
        // Perl's `exon` pre-predicate is 1 and `_skip_oc` drops the term.
        // Intron 1: 25_000_300 - 25_001_999
        // sv_start=25_000_200 (in exon 1), sv_end=25_002_100 (in exon 2)
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

    /// A symbolic insertion carried as the pair `start = POS + 1, end = POS` (a record
    /// whose END is its position or whose SVLEN is 0, ensembl-io `get_end`).
    fn make_pair_ins(pos: u64, class: vep_core::variant::VariantClass) -> InputVariant {
        let mut v = InputVariant::new("21".into(), pos + 1, pos, b"N".to_vec(), b"<INS>".to_vec());
        v.variant_class = class;
        v.is_structural = true;
        v.sv_end = Some(pos);
        v
    }

    fn pair_set(pos: u64, tx: &Transcript) -> Vec<String> {
        let v = make_pair_ins(pos, vep_core::variant::VariantClass::StructuralInsertion);
        let tc = calculate(&v, tx, 5000, 5000).expect("the pair lies within the neighbourhood");
        let mut got: Vec<String> = tc.consequences.iter().map(|c| c.to_string()).collect();
        got.sort();
        got
    }

    fn sorted(terms: &[&str]) -> Vec<String> {
        let mut v: Vec<String> = terms.iter().map(|t| (*t).to_string()).collect();
        v.sort();
        v
    }

    // The pair model, position by position, each set derived from Perl's predicates on
    // the pair (`overlap` holds only when both flanks lie inside the region; the `exon`
    // and `utr` pre-predicates read the sorted flank pair; `within_cdna` holds whenever a
    // flank is exonic). Layout as above; the intron 1 edges are 25_000_300 and
    // 25_001_999, the coding region 25_000_050..25_004_299 (CDS length 850).

    /// `11:7309883-7309884 insertion` on ENST00000318881 (1000 Genomes Phase 3, `<INS:MT>`
    /// with `END=POS`): inside an intron, `intron_variant` alone. Both flanks lie inside
    /// `[intron_start + 2, intron_end - 2]`, no flank is exonic, so no cDNA coordinate.
    #[test]
    fn test_pair_inside_intron_is_intron_variant() {
        let tx = make_test_transcript();
        assert_eq!(pair_set(25_000_500, &tx), sorted(&["intron_variant"]));
    }

    /// The pair between the intron's first and second bases: neither inside the trimmed
    /// intron nor the special cases (`start == intron_start + 2` is false: start is
    /// `intron_start + 1`), no exonic flank, so the term list is empty and the row takes
    /// Perl's default, `intergenic_variant`.
    #[test]
    fn test_pair_between_intron_bases_one_and_two_is_the_default() {
        let tx = make_test_transcript();
        assert_eq!(pair_set(25_000_300, &tx), sorted(&["intergenic_variant"]));
    }

    /// One base further in, `start == intron_start + 2`: the insertion special case of
    /// `_intron_effects` sets `intronic`.
    #[test]
    fn test_pair_at_intron_start_plus_two_is_intronic() {
        let tx = make_test_transcript();
        assert_eq!(pair_set(25_000_301, &tx), sorted(&["intron_variant"]));
    }

    /// `end == intron_end - 2`, the other special case; the pair `(intron_end - 2,
    /// intron_end - 1)` also lies in the polypyrimidine window `[intron_end - 16,
    /// intron_end - 2]` with no exonic flank.
    #[test]
    fn test_pair_at_intron_end_minus_two_is_intronic_with_polypyrimidine() {
        let tx = make_test_transcript();
        assert_eq!(
            pair_set(25_001_997, &tx),
            sorted(&["intron_variant", "splice_polypyrimidine_tract_variant"])
        );
    }

    /// Between the last base of exon 1 (CDS position 250) and the intron: the 5' flank is
    /// exonic, so `within_cdna` holds and `feature_elongation` with it; the insert
    /// coordinate is `(251, 250)`, within the CDS; neither intronic test holds.
    #[test]
    fn test_pair_at_donor_boundary_inside_cds_is_coding_and_elongation() {
        let tx = make_test_transcript();
        assert_eq!(
            pair_set(25_000_299, &tx),
            sorted(&["coding_sequence_variant", "feature_elongation"])
        );
    }

    /// Between the intron and the first base of exon 2 (CDS position 251): the 3' flank
    /// is exonic and its insert coordinate `(251, 250)` is within the CDS.
    #[test]
    fn test_pair_at_acceptor_boundary_inside_cds_is_coding_and_elongation() {
        let tx = make_test_transcript();
        assert_eq!(
            pair_set(25_001_999, &tx),
            sorted(&["coding_sequence_variant", "feature_elongation"])
        );
    }

    /// Between the last UTR base and the first coding base of exon 1: `start ==
    /// coding_region_start`, so `_before_coding` holds and the 5' UTR term with it; the
    /// insert coordinate in CDS space is `(1, 0)`, whose `end` is 0, so not within the
    /// CDS and no coding term.
    #[test]
    fn test_pair_before_first_coding_base_is_5_prime_utr_not_coding() {
        let tx = make_test_transcript();
        assert_eq!(
            pair_set(25_000_049, &tx),
            sorted(&["5_prime_UTR_variant", "feature_elongation"])
        );
    }

    /// After the last coding base of exon 3: `end == coding_region_end`, so
    /// `_after_coding` holds; the insert coordinate `(851, 850)` has `start > length`.
    #[test]
    fn test_pair_after_last_coding_base_is_3_prime_utr_not_coding() {
        let tx = make_test_transcript();
        assert_eq!(
            pair_set(25_004_299, &tx),
            sorted(&["3_prime_UTR_variant", "feature_elongation"])
        );
    }

    /// Inside the 5' UTR: the UTR term and `feature_elongation`.
    #[test]
    fn test_pair_inside_utr() {
        let tx = make_test_transcript();
        assert_eq!(
            pair_set(25_000_010, &tx),
            sorted(&["5_prime_UTR_variant", "feature_elongation"])
        );
    }

    /// Between the first and second bases of the start codon: both flanks inside the
    /// codon, so the structural `start_lost` arm holds beside `coding_sequence_variant`
    /// (`coding_unknown` excludes no start term for a structural allele).
    #[test]
    fn test_pair_inside_start_codon_is_start_lost_and_coding() {
        let tx = make_test_transcript();
        assert_eq!(
            pair_set(25_000_050, &tx),
            sorted(&[
                "coding_sequence_variant",
                "feature_elongation",
                "start_lost"
            ])
        );
    }

    /// One base before the transcript: `end < transcript start`, so not within the
    /// feature; `_before_start` reads `end`, distance 1.
    #[test]
    fn test_pair_before_transcript_start_is_upstream_at_distance_one() {
        let tx = make_test_transcript();
        let v = make_pair_ins(
            24_999_999,
            vep_core::variant::VariantClass::StructuralInsertion,
        );
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::UpstreamGeneVariant]
        );
        assert_eq!(tc.distance, Some(1));
    }

    /// One base after the transcript: `start > transcript end`, `_after_end` reads
    /// `start`, distance 1.
    #[test]
    fn test_pair_after_transcript_end_is_downstream_at_distance_one() {
        let tx = make_test_transcript();
        let v = make_pair_ins(
            25_006_000,
            vep_core::variant::VariantClass::StructuralInsertion,
        );
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::DownstreamGeneVariant]
        );
        assert_eq!(tc.distance, Some(1));
    }

    /// `7:124879111-124879112 Alu_insertion` on ENST00000420224 (1000 Genomes Phase 3,
    /// `<INS:ME:ALU>` with `SVLEN=0`): an intronic pair on a non-coding transcript is
    /// `intron_variant,non_coding_transcript_variant` (`within_non_coding_gene` holds
    /// through `within_transcript`), and in an exon `non_coding_transcript_exon_variant`
    /// with `feature_elongation`.
    #[test]
    fn test_pair_on_non_coding_transcript() {
        let tx = crate::sv::tests::make_non_coding_transcript();
        let v = make_pair_ins(
            25_000_500,
            vep_core::variant::VariantClass::MobileElementInsertion,
        );
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        let mut got: Vec<String> = tc.consequences.iter().map(|c| c.to_string()).collect();
        got.sort();
        assert_eq!(
            got,
            sorted(&["intron_variant", "non_coding_transcript_variant"])
        );
        let v = make_pair_ins(
            25_000_100,
            vep_core::variant::VariantClass::MobileElementInsertion,
        );
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        let mut got: Vec<String> = tc.consequences.iter().map(|c| c.to_string()).collect();
        got.sort();
        assert_eq!(
            got,
            sorted(&["feature_elongation", "non_coding_transcript_exon_variant"])
        );
    }

    /// The pair carries no overlap: VEP computes the overlap length from the allele's own
    /// start and end, 0 for a pair, and writes neither OverlapBP nor OverlapPC.
    #[test]
    fn test_pair_has_no_feature_overlap() {
        let tx = make_test_transcript();
        let v = make_pair_ins(
            25_000_500,
            vep_core::variant::VariantClass::StructuralInsertion,
        );
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert_eq!(tc.feature_overlap(v.start, v.end), None);
    }

    /// The single-base form stays what it was: a symbolic record with neither END nor
    /// SVLEN is the base at POS + 1 in Perl too (`get_end` returns `get_start`), with
    /// OverlapBP 1.
    #[test]
    fn test_point_insertion_unchanged_beside_the_pair() {
        let tx = make_test_transcript();
        let v = make_point_ins(25_000_500);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert_eq!(
            tc.feature_overlap(v.start, v.end).map(|(bp, _)| bp),
            Some(1)
        );
    }
}
