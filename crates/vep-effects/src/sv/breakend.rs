// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Breakend (BND / translocation) consequence calculation.
//!
//! Breakends represent chromosome-level rearrangements encoded in BND notation
//! (bracket/dot syntax). They are point-like: one breakpoint per record. The
//! mate breakpoint is on a different chromosome or distant position
//! and is handled as a separate variant record.
//!
//! Consequence logic (mirrors Perl VEP `StructuralVariationOverlapAllele`):
//!
//! - **Within transcript body (exon, intron, or UTR)**: always
//!   `feature_truncation` (the breakpoint disrupts the transcript), combined
//!   with the positional consequence (coding_sequence_variant, intron_variant,
//!   5_prime_UTR_variant, 3_prime_UTR_variant, etc.).
//! - **Upstream/downstream**: standard distance-based consequences.
//! - **No overlap at all**: None (caller handles intergenic).
//!
//! Perl citations name modules of ensembl-variation release/115
//! (`Bio/EnsEMBL/Variation/...`).

use super::{is_mature_mirna_sv, overlaps_any_exon};
use smallvec::SmallVec;
use vep_core::consequence::{
    Consequence, ConsequenceList, FeatureType, Impact, TranscriptConsequence,
};
use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

/// Calculate consequences for a breakend/translocation against a transcript.
///
/// The breakpoint is at `variant.start` (a single position); it is mapped into
/// transcript coordinates and consequences follow from where it falls.
pub fn calculate(
    variant: &InputVariant,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
) -> Option<TranscriptConsequence> {
    let bp = variant.start;
    let is_native_single_breakend = variant.is_single_breakend && variant.mate_id.is_none();

    if bp
        < transcript
            .start
            .saturating_sub(upstream_distance.max(downstream_distance))
        || bp > transcript.end + upstream_distance.max(downstream_distance)
    {
        return None;
    }

    let mut consequences: ConsequenceList = SmallVec::new();
    let mut distance: Option<u64> = None;
    let mut cdna_position: Option<String> = None;
    let mut cds_position: Option<String> = None;
    let mut protein_position: Option<String> = None;

    // has_cds() includes NMD transcripts, which have CDS/UTR like protein_coding.
    let is_coding = transcript.has_cds();

    let overlaps_transcript = bp >= transcript.start && bp <= transcript.end;

    // Perl VEP adds feature_truncation to an upstream/downstream BND whose mate
    // falls inside the transcript body.
    let same_chr_mate = variant
        .mate_chr
        .as_deref()
        .is_some_and(|mc| mc == &*variant.chr);
    let mate_inside_transcript = variant
        .mate_pos
        .filter(|_| same_chr_mate)
        .map(|mp| mp >= transcript.start && mp <= transcript.end)
        .unwrap_or(false);

    // Perl VEP emits upstream/downstream for a paired allele outside the body only
    // when the mate is also within the extended region; otherwise only the
    // single-breakend form gets those consequences.
    let is_paired_bnd = !variant.is_single_breakend && variant.mate_chr.is_some();
    let max_extend = upstream_distance.max(downstream_distance);
    let mate_near_transcript = variant
        .mate_pos
        .filter(|_| same_chr_mate)
        .map(|mp| {
            mp >= transcript.start.saturating_sub(max_extend) && mp <= transcript.end + max_extend
        })
        .unwrap_or(false);

    // A native single-breakend with sv_end and SVTYPE=DEL is a BND-format
    // deletion (gnomAD `N.` alleles with END set) and takes DEL-style
    // feature_truncation over its span.
    let is_bnd_format_deletion = is_native_single_breakend
        && variant.sv_end.is_some()
        && variant
            .sv_type
            .as_deref()
            .is_some_and(|st| st.eq_ignore_ascii_case("DEL"));

    // A single-breakend derived from a symbolic <BND> with an SVLEN span (end >
    // start, mate_id set) also takes span-based feature_truncation.
    let is_symbolic_bnd_span =
        !is_native_single_breakend && variant.is_single_breakend && variant.end > variant.start;

    let has_span_for_truncation = is_bnd_format_deletion || is_symbolic_bnd_span;
    let bnd_span_overlaps_transcript = if has_span_for_truncation {
        variant.sv_end.map(|sv_end| {
            let sv_start = variant.start;
            let tx_start = transcript.start;
            let tx_end = transcript.end;
            sv_end >= tx_start && sv_start <= tx_end
        })
    } else {
        None
    };

    if !overlaps_transcript && bp < transcript.start {
        let dist = transcript.start - bp;

        // A native single-breakend whose sv_end spans the transcript truncates.
        if bnd_span_overlaps_transcript == Some(true) {
            add_bnd_span_feature_truncation(variant, transcript, is_coding, &mut consequences);
        }

        match transcript.strand {
            Strand::Forward => {
                if dist <= upstream_distance {
                    if mate_inside_transcript {
                        consequences.push(Consequence::FeatureTruncation);
                    }
                    if !is_paired_bnd || mate_near_transcript {
                        consequences.push(Consequence::UpstreamGeneVariant);
                        distance = Some(dist);
                    }
                }
            }
            Strand::Reverse => {
                if dist <= downstream_distance {
                    if mate_inside_transcript {
                        consequences.push(Consequence::FeatureTruncation);
                    }
                    if !is_paired_bnd || mate_near_transcript {
                        consequences.push(Consequence::DownstreamGeneVariant);
                        distance = Some(dist);
                    }
                }
            }
        }
    } else if !overlaps_transcript && bp > transcript.end {
        let dist = bp - transcript.end;

        if bnd_span_overlaps_transcript == Some(true) {
            add_bnd_span_feature_truncation(variant, transcript, is_coding, &mut consequences);
        }

        match transcript.strand {
            Strand::Forward => {
                if dist <= downstream_distance {
                    if mate_inside_transcript {
                        consequences.push(Consequence::FeatureTruncation);
                    }
                    if !is_paired_bnd || mate_near_transcript {
                        consequences.push(Consequence::DownstreamGeneVariant);
                        distance = Some(dist);
                    }
                }
            }
            Strand::Reverse => {
                if dist <= upstream_distance {
                    if mate_inside_transcript {
                        consequences.push(Consequence::FeatureTruncation);
                    }
                    if !is_paired_bnd || mate_near_transcript {
                        consequences.push(Consequence::UpstreamGeneVariant);
                        distance = Some(dist);
                    }
                }
            }
        }
    } else {
        // feature_truncation logic (Perl VEP parity):
        // - Native single-breakend without sv_end: NO feature_truncation
        // - Native single-breakend with sv_end (BND-format DEL): DEL-style feature_truncation
        // - Paired BND with mate inside the same transcript: yes feature_truncation
        // - Paired BND with mate outside the transcript (or on different chr): NO feature_truncation
        // - Derived single-breakend from paired BND: feature_truncation unless in a mature miRNA
        let add_truncation = if is_native_single_breakend && !is_bnd_format_deletion {
            // No SVTYPE=DEL: no feature_truncation.
            false
        } else if has_span_for_truncation {
            // Span-based truncation is added below by add_bnd_span_feature_truncation.
            false
        } else if variant.is_single_breakend {
            // Perl emits only mature_miRNA_variant, no feature_truncation, when the
            // breakpoint falls within a mature miRNA exon.
            !is_mature_mirna_sv(transcript, bp, bp)
        } else if let Some(mate_pos) = variant.mate_pos {
            let same_chr = variant
                .mate_chr
                .as_deref()
                .map(|mc| mc == variant.chr)
                .unwrap_or(false);
            same_chr && mate_pos >= transcript.start && mate_pos <= transcript.end
        } else {
            // A paired BND always carries mate info; truncation is the safe default.
            true
        };
        if add_truncation {
            consequences.push(Consequence::FeatureTruncation);
        }

        if has_span_for_truncation {
            add_bnd_span_feature_truncation(variant, transcript, is_coding, &mut consequences);
        }

        let positional = classify_position_in_transcript(bp, transcript, is_coding);
        match positional {
            BreakendPosition::CodingExon {
                cdna_pos,
                cds_pos,
                prot_pos,
            } => {
                consequences.push(Consequence::CodingSequenceVariant);
                cdna_position = Some(cdna_pos.to_string());
                cds_position = Some(cds_pos.to_string());
                protein_position = Some(prot_pos.to_string());
            }
            BreakendPosition::FivePrimeUtr { cdna_pos } => {
                consequences.push(Consequence::FivePrimeUtrVariant);
                cdna_position = Some(cdna_pos.to_string());
            }
            BreakendPosition::ThreePrimeUtr { cdna_pos } => {
                consequences.push(Consequence::ThreePrimeUtrVariant);
                cdna_position = Some(cdna_pos.to_string());
            }
            BreakendPosition::Intron {
                dist_to_donor,
                dist_to_acceptor,
            } => {
                // trimmed intron boundaries: Ensembl sets its `intronic` flag with
                // `overlap($r_start, $r_end, $intron_start + 2, $intron_end - 2)`
                // (`BaseTranscriptVariationAllele.pm:143-150`). {`dist_to_donor`,
                // `dist_to_acceptor`} is a strand permutation of {`bp - intron.start`,
                // `intron.end - bp`}, so requiring both to be at least 2 is that same
                // interval with no strand case; an intron of 4 bp or less is never
                // `intronic`.
                if dist_to_donor >= 2 && dist_to_acceptor >= 2 {
                    consequences.push(Consequence::IntronVariant);
                }
                add_bnd_splice_consequences(&mut consequences, dist_to_donor, dist_to_acceptor);
                if !is_coding && !transcript.is_nmd_transcript() {
                    consequences.push(Consequence::NonCodingTranscriptVariant);
                }
            }
            BreakendPosition::NonCodingExon { cdna_pos } => {
                if is_mature_mirna_sv(transcript, bp, bp) {
                    consequences.push(Consequence::MatureMirnaVariant);
                } else {
                    consequences.push(Consequence::NonCodingTranscriptExonVariant);
                }
                // Perl: within_non_coding_gene requires not non_coding_exon_variant
                // and not mature_miRNA. So non_coding_transcript_variant is not added
                // when there IS exon overlap (either type).
                cdna_position = Some(cdna_pos.to_string());
            }
        }
    }

    // Perl's StructuralVariationOverlapAllele adds NMD_transcript_variant only
    // when the variant overlaps the transcript body, not for upstream/downstream.
    if transcript.is_nmd_transcript()
        && !consequences.is_empty()
        && !consequences.contains(&Consequence::NmdTranscriptVariant)
        && !consequences.iter().all(|c| {
            matches!(
                c,
                Consequence::UpstreamGeneVariant | Consequence::DownstreamGeneVariant
            )
        })
    {
        consequences.push(Consequence::NmdTranscriptVariant);
    }

    if consequences.is_empty() {
        return None;
    }

    let impact = consequences
        .iter()
        .map(|c| c.impact())
        .min_by_key(|i| match i {
            Impact::HIGH => 0,
            Impact::MODERATE => 1,
            Impact::LOW => 2,
            Impact::MODIFIER => 3,
        })
        .unwrap_or(Impact::MODIFIER);

    consequences.sort_by_key(|c| c.rank());

    Some(TranscriptConsequence {
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
        cdna_position,
        cds_position,
        protein_position,
        amino_acids: None,
        codons: None,
        protein_id: if is_coding {
            transcript.protein_id.clone()
        } else {
            None
        },
        distance,
        strand: match transcript.strand {
            Strand::Forward => 1,
            Strand::Reverse => -1,
        },
        exon: None,
        intron: None,
        hgvsc: None,
        hgvsp: None,
        hgvs_offset: None,
        sift: None,
        polyphen: None,
        domains: vec![],
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
    })
}

/// Calculate consequences for the paired mate-side breakpoint annotation.
///
/// Perl VEP emits a simpler consequence model for the remote mate breakpoint:
/// when the mate position falls within a transcript body it reports only
/// `feature_truncation`; upstream/downstream distance terms are preserved.
///
/// `local_span_is_ranged` says whether the original local variant span covers
/// more than one base, i.e. `end > start` before the caller repointed the
/// variant at the mate. It gates the context terms, and the asymmetry is Perl's:
/// `feature_truncation` (`VariationEffect.pm:350-359`) passes
/// `$bvfoa->breakend`, the mate coordinates, into `within_feature`, while every
/// context predicate (`within_intron:629`, `within_cds:643`) takes the default
/// `$bvf`, the local variation feature. A point breakend therefore satisfies
/// none of them at the mate position and Perl reports bare
/// `feature_truncation`. A fully ranged local span can satisfy them, which is
/// why the gate is on ranged-ness rather than a blanket removal.
pub fn calculate_paired_mate(
    variant: &InputVariant,
    transcript: &Transcript,
    upstream_distance: u64,
    downstream_distance: u64,
    local_span_is_ranged: bool,
) -> Option<TranscriptConsequence> {
    let bp = variant.start;

    if bp
        < transcript
            .start
            .saturating_sub(upstream_distance.max(downstream_distance))
        || bp > transcript.end + upstream_distance.max(downstream_distance)
    {
        return None;
    }

    let mut consequences: ConsequenceList = SmallVec::new();
    let distance: Option<u64> = None;

    if bp >= transcript.start && bp <= transcript.end {
        consequences.push(Consequence::FeatureTruncation);

        // Context terms only for same-chromosome paired BNDs with a ranged local
        // span: Perl's on-demand loading annotates fewer mate-side transcripts on
        // an inter-chromosomal BND, and its context predicates read the local
        // feature, not the mate (see the function doc).
        let same_chr = variant
            .mate_chr
            .as_deref()
            .is_some_and(|mc| mc == &*variant.chr);
        if same_chr && local_span_is_ranged {
            let is_coding = transcript.has_cds();
            if is_coding {
                if let (Some(cds_start), Some(cds_end)) =
                    (transcript.coding_region_start, transcript.coding_region_end)
                {
                    if bp >= cds_start && bp <= cds_end {
                        consequences.push(Consequence::CodingSequenceVariant);
                    }
                }
                if super::overlaps_five_prime_utr(transcript, bp, bp) {
                    consequences.push(Consequence::FivePrimeUtrVariant);
                }
                if super::overlaps_three_prime_utr(transcript, bp, bp) {
                    consequences.push(Consequence::ThreePrimeUtrVariant);
                }
            } else if !transcript.is_nmd_transcript() {
                if super::overlaps_any_exon(transcript, bp, bp) {
                    consequences.push(Consequence::NonCodingTranscriptExonVariant);
                } else {
                    consequences.push(Consequence::NonCodingTranscriptVariant);
                }
            }
            // Trimmed boundaries, as in the `BreakendPosition::Intron` arm above.
            if super::overlaps_any_intron_trimmed(transcript, bp, bp) {
                consequences.push(Consequence::IntronVariant);
            }
            if transcript.is_nmd_transcript() {
                consequences.push(Consequence::NmdTranscriptVariant);
            }
        }
    } else {
        consequences.push(Consequence::IntergenicVariant);
    }

    if consequences.is_empty() {
        return None;
    }

    let impact = consequences
        .iter()
        .map(|c| c.impact())
        .min_by_key(|i| match i {
            Impact::HIGH => 0,
            Impact::MODERATE => 1,
            Impact::LOW => 2,
            Impact::MODIFIER => 3,
        })
        .unwrap_or(Impact::MODIFIER);
    consequences.sort_by_key(|c| c.rank());

    Some(TranscriptConsequence {
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
        protein_id: if transcript.is_protein_coding() {
            transcript.protein_id.clone()
        } else {
            None
        },
        distance,
        strand: match transcript.strand {
            Strand::Forward => 1,
            Strand::Reverse => -1,
        },
        exon: None,
        intron: None,
        hgvsc: None,
        hgvsp: None,
        hgvs_offset: None,
        sift: None,
        polyphen: None,
        domains: vec![],
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
    })
}

/// Apply DEL-style feature_truncation for a native single-breakend with `sv_end`.
///
/// BND-format deletions in gnomAD use single-breakend `N.` notation with SVTYPE=DEL
/// and END set. The SV span is `(variant.start, sv_end)`. This function mirrors the
/// deletion.rs feature_truncation predicate using that span.
fn add_bnd_span_feature_truncation(
    variant: &InputVariant,
    transcript: &Transcript,
    is_coding: bool,
    consequences: &mut ConsequenceList,
) {
    let sv_end = match variant.sv_end {
        Some(e) => e,
        None => return,
    };
    let sv_start = variant.start;
    let tx_start = transcript.start;
    let tx_end = transcript.end;

    if sv_end < tx_start || sv_start > tx_end {
        return;
    }

    let extends_beyond = sv_start < tx_start || sv_end > tx_end;
    let complete_within = sv_start >= tx_start && sv_end <= tx_end;

    let has_exon_overlap = overlaps_any_exon(transcript, sv_start, sv_end);
    let within_cdna = if variant.is_structural && extends_beyond {
        true
    } else {
        has_exon_overlap
    };

    if within_cdna && (extends_beyond || complete_within) {
        consequences.push(Consequence::FeatureTruncation);
        // Intron-only overlap gets the context term; an NMD transcript gets
        // NMD_transcript_variant, not coding_transcript_variant.
        if !has_exon_overlap {
            if transcript.is_nmd_transcript() {
                consequences.push(Consequence::NmdTranscriptVariant);
            } else if is_coding {
                consequences.push(Consequence::CodingTranscriptVariant);
            } else {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
        }
    }
}

/// Position classification for a breakpoint within a transcript body.
enum BreakendPosition {
    CodingExon {
        cdna_pos: u64,
        cds_pos: u64,
        prot_pos: u64,
    },
    FivePrimeUtr {
        cdna_pos: u64,
    },
    ThreePrimeUtr {
        cdna_pos: u64,
    },
    Intron {
        dist_to_donor: u64,
        dist_to_acceptor: u64,
    },
    NonCodingExon {
        cdna_pos: u64,
    },
}

/// Classify where a genomic breakpoint falls within a transcript.
///
/// The breakpoint `bp` is guaranteed to be within `transcript.start..=transcript.end`.
fn classify_position_in_transcript(
    bp: u64,
    transcript: &Transcript,
    is_coding: bool,
) -> BreakendPosition {
    // Check introns first (they occupy more genomic space).
    for intron in &transcript.introns {
        if bp >= intron.start && bp <= intron.end {
            let dist_to_genomic_start = bp - intron.start;
            let dist_to_genomic_end = intron.end - bp;
            // Donor/acceptor assignment is strand-dependent:
            // Forward strand: donor at intron start, acceptor at intron end
            // Reverse strand: donor at intron end, acceptor at intron start
            let (dist_to_donor, dist_to_acceptor) = match transcript.strand {
                Strand::Forward => (dist_to_genomic_start, dist_to_genomic_end),
                Strand::Reverse => (dist_to_genomic_end, dist_to_genomic_start),
            };
            return BreakendPosition::Intron {
                dist_to_donor,
                dist_to_acceptor,
            };
        }
    }

    let cdna_pos = genomic_to_cdna(bp, transcript);

    if !is_coding {
        return BreakendPosition::NonCodingExon {
            cdna_pos: cdna_pos.unwrap_or(1),
        };
    }

    let cdna_coding_start = transcript.cdna_coding_start.unwrap_or(1);
    let cdna_coding_end = transcript.cdna_coding_end.unwrap_or(u64::MAX);

    if let Some(cp) = cdna_pos {
        if cp < cdna_coding_start {
            // cDNA positions are already transcript-oriented, so "before CDS start"
            // is always 5' UTR regardless of genomic strand.
            BreakendPosition::FivePrimeUtr { cdna_pos: cp }
        } else if cp > cdna_coding_end {
            BreakendPosition::ThreePrimeUtr { cdna_pos: cp }
        } else {
            let cds_pos = cp - cdna_coding_start + 1;
            let prot_pos = (cds_pos - 1) / 3 + 1;
            BreakendPosition::CodingExon {
                cdna_pos: cp,
                cds_pos,
                prot_pos,
            }
        }
    } else {
        // Exonic, but the mapper resolved no cDNA position: a coding exon with
        // unknown positions.
        BreakendPosition::CodingExon {
            cdna_pos: 0,
            cds_pos: 0,
            prot_pos: 0,
        }
    }
}

/// Map a genomic position to cDNA position using the transcript mapper pairs.
///
/// Returns `None` if the position doesn't map to any exon (shouldn't happen if
/// the caller already excluded introns).
fn genomic_to_cdna(bp: u64, transcript: &Transcript) -> Option<u64> {
    let vefc = transcript.vefc.as_ref()?;
    let mapper = vefc.mapper.as_ref()?;

    for pair in &mapper.exon_coord_mapper.pairs {
        if bp >= pair.to_start && bp <= pair.to_end {
            // Forward strand: offset from start of exonic block.
            // Reverse strand: offset from end of exonic block.
            let cdna = match transcript.strand {
                Strand::Forward => pair.from_start + (bp - pair.to_start),
                Strand::Reverse => pair.from_start + (pair.to_end - bp),
            };
            return Some(cdna);
        }
    }
    None
}

/// Add splice-region consequences for a breakpoint in an intron.
///
/// Perl VEP uses a simplified splice model for structural variants:
/// only `splice_polypyrimidine_tract_variant` is added when the breakpoint
/// falls within the polypyrimidine tract (2-16 bases from acceptor).
/// Detailed splice sub-terms (splice_region_variant, splice_donor_5th_base,
/// splice_donor_region) are not added for SVs; those apply only to small
/// variants in the main consequence engine.
fn add_bnd_splice_consequences(
    consequences: &mut ConsequenceList,
    _dist_to_donor: u64,
    dist_to_acceptor: u64,
) {
    // Polypyrimidine tract (positions 2..=16 from acceptor).
    if (2..=16).contains(&dist_to_acceptor) {
        consequences.push(Consequence::SplicePolypyrimidineTractVariant);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::make_test_transcript;
    use vep_core::variant::VariantClass;

    /// Build a BND variant at a given position.
    fn make_bnd(pos: u64) -> InputVariant {
        let mut v = InputVariant::new("21".into(), pos, pos, b"N".to_vec(), b"-".to_vec());
        v.variant_class = VariantClass::Translocation;
        v.is_structural = true;
        v.mate_id = Some("bnd_mate".into());
        v
    }

    fn make_native_single_breakend(pos: u64) -> InputVariant {
        let mut v = make_bnd(pos);
        v.alt_alleles = vec![b"A.".to_vec()];
        v.allele_string = "N/A.".into();
        v.is_single_breakend = true;
        v.mate_id = None;
        v
    }

    // ---- Test transcript layout (from make_test_transcript) ----
    //   Exon 1:   25_000_000 - 25_000_299  (cDNA 1-300)
    //   Intron 1: 25_000_300 - 25_001_999
    //   Exon 2:   25_002_000 - 25_002_299  (cDNA 301-600)
    //   Intron 2: 25_002_300 - 25_003_999
    //   Exon 3:   25_004_000 - 25_006_000  (cDNA 601-2601)
    //
    //   CDS: cDNA 51-900  (genomic 25_000_050 - 25_004_299)
    //   5' UTR: cDNA 1-50  (genomic 25_000_000 - 25_000_049)
    //   3' UTR: cDNA 901-2601 (genomic 25_004_300 - 25_006_000)

    #[test]
    fn test_bnd_upstream_forward_strand() {
        let tx = make_test_transcript();
        let v = make_bnd(24_999_000); // 1000bp upstream
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::UpstreamGeneVariant));
        assert_eq!(tc.impact, Impact::MODIFIER);
        assert_eq!(tc.distance, Some(1000));
    }

    #[test]
    fn test_bnd_upstream_paired_with_mate_inside_gets_truncation_and_upstream() {
        // Paired BND upstream of transcript, mate falls inside transcript body.
        // Perl VEP: emits feature_truncation + upstream_gene_variant together.
        let tx = make_test_transcript(); // 25_000_000 - 25_006_000
        let mut v = make_bnd(24_999_000); // 1000bp upstream
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(25_003_000); // Inside transcript body
        v.is_single_breakend = false;
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "Paired BND upstream with mate inside transcript should have feature_truncation"
        );
        assert!(
            tc.consequences.contains(&Consequence::UpstreamGeneVariant),
            "Paired BND with mate inside should also get upstream_gene_variant (Perl VEP parity)"
        );
        assert_eq!(tc.impact, Impact::HIGH);
        assert_eq!(tc.distance, Some(1000));
    }

    #[test]
    fn test_bnd_upstream_paired_with_mate_outside_returns_none() {
        // Paired BND upstream of transcript, mate also outside → no consequences at all.
        // Perl VEP: paired alleles outside the transcript body with no mate inside
        // produce no transcript-level consequence.
        let tx = make_test_transcript();
        let mut v = make_bnd(24_999_000); // 1000bp upstream
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(24_990_000); // Also outside transcript
        v.is_single_breakend = false;
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(
            result.is_none(),
            "Paired BND upstream with mate outside should return None (no consequences)"
        );
    }

    #[test]
    fn test_bnd_downstream_forward_strand() {
        let tx = make_test_transcript();
        let v = make_bnd(25_007_000); // 1000bp downstream
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert_eq!(tc.impact, Impact::MODIFIER);
        assert_eq!(tc.distance, Some(1000));
    }

    #[test]
    fn test_bnd_too_far_returns_none() {
        let tx = make_test_transcript();
        let v = make_bnd(24_990_000); // 10000bp upstream, beyond 5000 limit
        let result = calculate(&v, &tx, 5000, 5000);
        assert!(result.is_none());
    }

    #[test]
    fn test_bnd_in_5prime_utr() {
        let tx = make_test_transcript();
        // 5' UTR is genomic 25_000_000 - 25_000_049 (cDNA 1-50).
        let v = make_bnd(25_000_025);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FivePrimeUtrVariant));
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert_eq!(tc.impact, Impact::HIGH);
        assert_eq!(tc.cdna_position, Some("26".to_string()));
    }

    #[test]
    fn test_bnd_in_3prime_utr() {
        let tx = make_test_transcript();
        // 3' UTR is genomic 25_004_300 - 25_006_000 (cDNA 901-2601).
        let v = make_bnd(25_004_500);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::ThreePrimeUtrVariant));
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_bnd_in_coding_exon() {
        let tx = make_test_transcript();
        // CDS in exon 1: genomic 25_000_050 - 25_000_299 (cDNA 51-300, CDS 1-250).
        let v = make_bnd(25_000_100);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert_eq!(tc.impact, Impact::HIGH);
        // cDNA pos: 25_000_100 - 25_000_000 + 1 = 101
        assert_eq!(tc.cdna_position, Some("101".to_string()));
        // CDS pos: 101 - 51 + 1 = 51
        assert_eq!(tc.cds_position, Some("51".to_string()));
        // Protein pos: (51 - 1) / 3 + 1 = 17
        assert_eq!(tc.protein_position, Some("17".to_string()));
    }

    #[test]
    fn test_bnd_in_coding_exon2() {
        let tx = make_test_transcript();
        // CDS in exon 2: genomic 25_002_000 - 25_002_299 (cDNA 301-600, CDS 251-550).
        let v = make_bnd(25_002_150);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        assert_eq!(tc.impact, Impact::HIGH);
        // cDNA pos: 301 + (25_002_150 - 25_002_000) = 451
        assert_eq!(tc.cdna_position, Some("451".to_string()));
        // CDS pos: 451 - 51 + 1 = 401
        assert_eq!(tc.cds_position, Some("401".to_string()));
    }

    #[test]
    fn test_bnd_in_intron() {
        let tx = make_test_transcript();
        let v = make_bnd(25_001_000); // Middle of intron 1
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        // Perl VEP adds feature_truncation for all breakpoints within the
        // transcript body, including introns.
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert_eq!(tc.impact, Impact::HIGH);
        assert!(tc.cdna_position.is_none());
    }

    #[test]
    fn test_bnd_intron_near_donor_splice_site() {
        let tx = make_test_transcript();
        // Intron 1 starts at 25_000_300. Position +0 from donor = splice donor site.
        // Perl VEP does not add splice_donor_variant for SVs, only polypyrimidine tract.
        let v = make_bnd(25_000_300);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        // No `intron_variant` either: this position is one of the two invariant donor
        // bases Ensembl's `intronic` flag excludes via `$intron_start + 2`.
        assert!(!tc.consequences.contains(&Consequence::IntronVariant));
        assert!(
            !tc.consequences.contains(&Consequence::SpliceDonorVariant),
            "BNDs should not get splice_donor_variant (Perl VEP parity)"
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_bnd_intron_near_acceptor_splice_site() {
        let tx = make_test_transcript();
        // Intron 1 ends at 25_001_999. Position at end = splice acceptor site.
        // Perl VEP does not add splice_acceptor_variant for SVs.
        let v = make_bnd(25_001_999);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        // No `intron_variant`: an invariant acceptor base, excluded by `$intron_end - 2`.
        assert!(!tc.consequences.contains(&Consequence::IntronVariant));
        assert!(
            !tc.consequences
                .contains(&Consequence::SpliceAcceptorVariant),
            "BNDs should not get splice_acceptor_variant (Perl VEP parity)"
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    /// Pins the exact trimmed boundary: offsets 0 and 1 from either intron end are
    /// not `intronic`, offset 2 is the first base that is, per
    /// `overlap($r_start, $r_end, $intron_start + 2, $intron_end - 2)`
    /// (`BaseTranscriptVariationAllele.pm:143-150`).
    #[test]
    fn test_bnd_intron_variant_starts_two_bases_inside_the_intron() {
        let tx = make_test_transcript();
        // Intron 1: 25_000_300 - 25_001_999.
        for (offset_from_start, want_intronic) in [(0u64, false), (1, false), (2, true), (3, true)]
        {
            let v = make_bnd(25_000_300 + offset_from_start);
            let tc = calculate(&v, &tx, 5000, 5000).unwrap();
            assert_eq!(
                tc.consequences.contains(&Consequence::IntronVariant),
                want_intronic,
                "donor side, offset {offset_from_start}: got {:?}",
                tc.consequences
            );
        }
        for (offset_from_end, want_intronic) in [(0u64, false), (1, false), (2, true), (3, true)] {
            let v = make_bnd(25_001_999 - offset_from_end);
            let tc = calculate(&v, &tx, 5000, 5000).unwrap();
            assert_eq!(
                tc.consequences.contains(&Consequence::IntronVariant),
                want_intronic,
                "acceptor side, offset {offset_from_end}: got {:?}",
                tc.consequences
            );
        }
    }

    #[test]
    fn test_bnd_intron_polypyrimidine_tract() {
        let tx = make_test_transcript();
        // Intron 1 ends at 25_001_999. Position 10 bases from acceptor.
        // dist_to_acceptor = 25_001_999 - 25_001_989 = 10 (in polypyrimidine range 2..=16).
        let v = make_bnd(25_001_989);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::SplicePolypyrimidineTractVariant));
    }

    #[test]
    fn test_bnd_intron_polypyrimidine_tract_reverse_strand() {
        // Reverse strand: donor is at intron END, acceptor at intron start.
        // Intron 1: 25_000_300 - 25_001_999.
        // Position 25_000_310 is 10bp from genomic start (=acceptor on reverse strand).
        // dist_to_acceptor = 25_000_310 - 25_000_300 = 10 (in polypyrimidine range 2..=16).
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;
        if let Some(ref mut vefc) = tx.vefc {
            if let Some(ref mut mapper) = vefc.mapper {
                for pair in &mut mapper.exon_coord_mapper.pairs {
                    pair.ori = -1;
                }
            }
        }
        let v = make_bnd(25_000_310);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(
            tc.consequences
                .contains(&Consequence::SplicePolypyrimidineTractVariant),
            "Reverse-strand BND 10bp from acceptor (genomic start) should have polypyrimidine tract"
        );
    }

    #[test]
    fn test_bnd_intron_no_splice_donor_reverse_strand() {
        // Reverse strand: donor is at intron END (25_001_999).
        // Perl VEP does not add splice_donor_variant for SVs.
        let mut tx = make_test_transcript();
        tx.strand = Strand::Reverse;
        if let Some(ref mut vefc) = tx.vefc {
            if let Some(ref mut mapper) = vefc.mapper {
                for pair in &mut mapper.exon_coord_mapper.pairs {
                    pair.ori = -1;
                }
            }
        }
        let v = make_bnd(25_001_999);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        // No `intron_variant`: the trimmed exclusion is symmetric in genomic
        // coordinates, so it holds on the reverse strand at the same two bases.
        assert!(!tc.consequences.contains(&Consequence::IntronVariant));
        assert!(
            !tc.consequences.contains(&Consequence::SpliceDonorVariant),
            "BNDs should not get splice_donor_variant (Perl VEP parity)"
        );
    }

    #[test]
    fn test_bnd_feature_truncation_mate_inside_transcript() {
        // Paired BND with mate inside the transcript → feature_truncation added.
        let tx = make_test_transcript();
        let mut v = make_bnd(25_001_000); // Intron 1
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(25_003_000); // Inside transcript (intron 2)
        v.is_single_breakend = false;
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "BND with mate inside transcript should have feature_truncation"
        );
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_bnd_no_feature_truncation_mate_outside_transcript() {
        // Paired BND with mate outside the transcript → NO feature_truncation.
        // Perl VEP parity: only positional consequence (intron_variant).
        let tx = make_test_transcript();
        let mut v = make_bnd(25_001_000); // Intron 1 (inside transcript)
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(24_990_000); // Outside transcript (before start)
        v.is_single_breakend = false;
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            !tc.consequences.contains(&Consequence::FeatureTruncation),
            "BND with mate outside transcript should not have feature_truncation"
        );
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert_eq!(tc.impact, Impact::MODIFIER); // No HIGH impact without truncation
    }

    #[test]
    fn test_bnd_no_feature_truncation_mate_different_chromosome() {
        // Paired BND with mate on different chromosome → NO feature_truncation.
        let tx = make_test_transcript();
        let mut v = make_bnd(25_001_000); // Intron 1 (inside transcript)
        v.mate_chr = Some("1".into()); // Different chromosome
        v.mate_pos = Some(12345678);
        v.is_single_breakend = false;
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            !tc.consequences.contains(&Consequence::FeatureTruncation),
            "BND with mate on different chromosome should not have feature_truncation"
        );
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_bnd_non_coding_intron() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;
        tx.translation = None;
        tx.protein_id = None;

        let v = make_bnd(25_001_000); // Intron 1
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        // Perl VEP adds feature_truncation for all BND positions within the
        // transcript body, including introns of non-coding transcripts.
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert!(tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptVariant));
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_bnd_non_coding_exon() {
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;
        tx.translation = None;
        tx.protein_id = None;

        let v = make_bnd(25_000_100); // Exon 1
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert!(tc
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
        // Perl: within_non_coding_gene requires not non_coding_exon_variant.
        // So non_coding_transcript_variant is not added when there IS exon overlap.
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "BND in non-coding exon should not have non_coding_transcript_variant"
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_bnd_at_transcript_start_boundary() {
        let tx = make_test_transcript();
        // Exactly at transcript start (25_000_000), which is in the 5' UTR.
        let v = make_bnd(25_000_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FivePrimeUtrVariant));
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
    }

    #[test]
    fn test_bnd_at_transcript_end_boundary() {
        let tx = make_test_transcript();
        // Exactly at transcript end (25_006_000), which is in the 3' UTR.
        let v = make_bnd(25_006_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::ThreePrimeUtrVariant));
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
    }

    #[test]
    fn test_bnd_no_splice_donor_sub_terms() {
        let tx = make_test_transcript();
        // Intron 1 starts at 25_000_300. Position +4 = donor 5th base.
        // Perl VEP does not add splice_donor_5th_base or splice_region for SVs.
        let v = make_bnd(25_000_304);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            !tc.consequences
                .contains(&Consequence::SpliceDonor5thBaseVariant),
            "BNDs should not get splice sub-terms (Perl VEP parity)"
        );
        assert!(!tc
            .consequences
            .contains(&Consequence::SpliceDonorRegionVariant));
        assert!(!tc.consequences.contains(&Consequence::SpliceRegionVariant));
    }

    #[test]
    fn test_bnd_no_splice_donor_region() {
        let tx = make_test_transcript();
        // Intron 1 starts at 25_000_300. Position +3 = donor region.
        // Perl VEP does not add splice_donor_region or splice_region for SVs.
        let v = make_bnd(25_000_303);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            !tc.consequences
                .contains(&Consequence::SpliceDonorRegionVariant),
            "BNDs should not get splice_donor_region (Perl VEP parity)"
        );
        assert!(
            !tc.consequences.contains(&Consequence::SpliceRegionVariant),
            "BNDs should not get splice_region_variant (Perl VEP parity)"
        );
    }

    #[test]
    fn test_bnd_strand_reported_correctly() {
        let tx = make_test_transcript();
        let v = make_bnd(25_000_100);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert_eq!(tc.strand, 1);
    }

    #[test]
    fn test_bnd_metadata_fields() {
        let tx = make_test_transcript();
        let v = make_bnd(25_000_100);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert_eq!(&*tc.transcript_id, "ENST00000000001");
        assert_eq!(&*tc.gene_id, "ENSG00000000001");
        assert_eq!(tc.gene_symbol.as_deref(), Some("TEST1"));
        assert_eq!(tc.biotype.as_deref(), Some("protein_coding"));
        assert!(tc.canonical);
        assert_eq!(tc.feature_type, FeatureType::Transcript);
    }

    #[test]
    fn test_bnd_coding_exon_in_exon3() {
        let tx = make_test_transcript();
        // CDS in exon 3: genomic 25_004_000 - 25_004_299 (cDNA 601-900, CDS 551-850).
        let v = make_bnd(25_004_100);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert!(tc
            .consequences
            .contains(&Consequence::CodingSequenceVariant));
        // cDNA pos: 601 + (25_004_100 - 25_004_000) = 701
        assert_eq!(tc.cdna_position, Some("701".to_string()));
        // CDS pos: 701 - 51 + 1 = 651
        assert_eq!(tc.cds_position, Some("651".to_string()));
    }

    #[test]
    fn test_bnd_intron2() {
        let tx = make_test_transcript();
        // Intron 2: 25_002_300 - 25_003_999.
        let v = make_bnd(25_003_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(tc.consequences.contains(&Consequence::FeatureTruncation));
        assert!(tc.consequences.contains(&Consequence::IntronVariant));
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_paired_mate_intron_is_feature_truncation_only() {
        let tx = make_test_transcript();
        let v = make_bnd(25_001_000);
        let tc = calculate_paired_mate(&v, &tx, 5000, 5000, false).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::FeatureTruncation]
        );
        assert_eq!(tc.impact, Impact::HIGH);
        assert!(tc.cdna_position.is_none());
    }

    #[test]
    fn test_paired_mate_nearby_non_overlap_is_intergenic() {
        let tx = make_test_transcript();
        let v = make_bnd(24_999_000);
        let tc = calculate_paired_mate(&v, &tx, 5000, 5000, false).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::IntergenicVariant]
        );
        assert_eq!(tc.impact, Impact::MODIFIER);
        assert!(tc.distance.is_none());
    }

    /// A point local breakend gets bare `feature_truncation` at the mate,
    /// even when the mate lands mid-transcript on the same chromosome.
    ///
    /// Perl's `feature_truncation` (`VariationEffect.pm:350-359`) passes
    /// `$bvfoa->breakend`, the mate coordinates, into `within_feature`. Every
    /// context predicate (`within_intron:629`, `within_cds:643`) instead takes the
    /// default `$bvf`, the local variation feature. A point local span satisfies
    /// none of them, so Perl reports only `feature_truncation`.
    #[test]
    fn test_paired_mate_point_local_span_gets_no_context_terms() {
        let tx = make_test_transcript();
        // Mate lands in intron 1 (25_000_300 - 25_001_999), same chromosome.
        let mut v = make_bnd(25_001_000);
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(25_001_000);
        let tc = calculate_paired_mate(&v, &tx, 5000, 5000, false).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::FeatureTruncation],
            "Perl's context predicates read the LOCAL feature, which for a point \
             breakend satisfies none of them"
        );
    }

    /// Discriminator: a ranged local span does still get context terms.
    ///
    /// This is why the gate is on ranged-ness rather than removing the context
    /// block: a ranged local span can satisfy Perl's context predicates.
    #[test]
    fn test_paired_mate_ranged_local_span_keeps_context_terms() {
        let tx = make_test_transcript();
        let mut v = make_bnd(25_001_000);
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(25_001_000);
        let tc = calculate_paired_mate(&v, &tx, 5000, 5000, true).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "a ranged local span can satisfy Perl's within_intron; got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_native_single_breakend_in_intron_has_no_feature_truncation() {
        let tx = make_test_transcript();
        let v = make_native_single_breakend(25_001_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert_eq!(tc.consequences.to_vec(), vec![Consequence::IntronVariant]);
        assert_eq!(tc.impact, Impact::MODIFIER);
    }

    #[test]
    fn test_native_single_breakend_in_three_prime_utr_has_no_feature_truncation() {
        let tx = make_test_transcript();
        let v = make_native_single_breakend(25_004_500);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert_eq!(
            tc.consequences.to_vec(),
            vec![Consequence::ThreePrimeUtrVariant]
        );
        assert_eq!(tc.impact, Impact::MODIFIER);
    }

    #[test]
    fn test_bnd_consequences_sorted_by_rank() {
        let tx = make_test_transcript();
        let v = make_bnd(25_000_300);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();

        // Verify consequences are sorted by rank (lower rank = more severe = first).
        for window in tc.consequences.windows(2) {
            assert!(
                window[0].rank() <= window[1].rank(),
                "Consequences not sorted: {:?} (rank {}) should come before {:?} (rank {})",
                window[0],
                window[0].rank(),
                window[1],
                window[1].rank(),
            );
        }
    }

    /// Helper: build a native single-breakend with sv_end and SVTYPE=DEL (BND-format deletion).
    fn make_native_single_breakend_with_sv_end(pos: u64, sv_end: u64) -> InputVariant {
        let mut v = make_native_single_breakend(pos);
        v.sv_end = Some(sv_end);
        v.is_structural = true;
        v.sv_type = Some("DEL".into()); // Required: marks this as a BND-format deletion
        v
    }

    #[test]
    fn test_native_single_breakend_with_sv_end_extending_beyond_transcript() {
        // SV span extends beyond both transcript boundaries → feature_truncation.
        // Breakpoint upstream (within 5000bp distance), sv_end downstream of transcript.
        let tx = make_test_transcript(); // 25_000_000 - 25_006_000
        let v = make_native_single_breakend_with_sv_end(24_998_000, 25_010_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "BND-format DEL spanning beyond transcript should have feature_truncation"
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_native_single_breakend_with_sv_end_within_transcript() {
        // SV span entirely within transcript → feature_truncation (complete_within).
        // Breakpoint at exon 1, sv_end in exon 2 → has exon overlap.
        let tx = make_test_transcript();
        let v = make_native_single_breakend_with_sv_end(25_000_100, 25_002_200);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "BND-format DEL within transcript (exon overlap) should have feature_truncation"
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_native_single_breakend_with_sv_end_upstream_spanning_into_transcript() {
        // Breakpoint upstream (outside transcript), sv_end inside transcript body.
        // The SV span overlaps the transcript and extends beyond its start.
        let tx = make_test_transcript(); // 25_000_000 - 25_006_000
        let v = make_native_single_breakend_with_sv_end(24_998_000, 25_002_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "BND-format DEL from upstream into transcript should have feature_truncation"
        );
        assert!(
            tc.consequences.contains(&Consequence::UpstreamGeneVariant),
            "Breakpoint is upstream so upstream_gene_variant should be present"
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_native_single_breakend_with_sv_end_crossing_one_boundary_truncates() {
        // A span that crosses one transcript boundary and closes inside the body
        // gets feature_truncation whichever end it crosses. The first span opens
        // upstream and closes in intron 1 (covering exon 1); the second opens in
        // intron 1 and runs past the transcript end (covering exons 2 and 3).
        let tx = make_test_transcript(); // coding
        let v = make_native_single_breakend_with_sv_end(24_999_000, 25_000_500);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "BND-format DEL opening upstream of the transcript should have feature_truncation"
        );
        let v2 = make_native_single_breakend_with_sv_end(25_000_500, 25_010_000);
        let tc2 = calculate(&v2, &tx, 5000, 5000).unwrap();
        assert!(
            tc2.consequences.contains(&Consequence::FeatureTruncation),
            "BND-format DEL running past the transcript end should have feature_truncation"
        );
    }

    #[test]
    fn test_native_single_breakend_without_sv_end_gets_no_feature_truncation() {
        // A native single-breakend without sv_end gets no feature_truncation, even
        // with the breakpoint inside the transcript.
        let tx = make_test_transcript();
        let v = make_native_single_breakend(25_001_000); // Intron 1
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            !tc.consequences.contains(&Consequence::FeatureTruncation),
            "Native single-breakend without sv_end should not have feature_truncation"
        );
        assert_eq!(tc.consequences.to_vec(), vec![Consequence::IntronVariant]);
    }

    #[test]
    fn test_native_single_breakend_with_sv_end_non_coding_context_term() {
        // Non-coding transcript: intron-only overlap should get
        // non_coding_transcript_variant as context term.
        let mut tx = make_test_transcript();
        tx.biotype = "lncRNA".into();
        tx.coding_region_start = None;
        tx.coding_region_end = None;
        tx.cdna_coding_start = None;
        tx.cdna_coding_end = None;
        tx.translation = None;
        tx.protein_id = None;

        // The span from 25_000_500 (intron 1) past the transcript end also overlaps
        // exon 2, so this asserts only that feature_truncation is added for a
        // non-coding span overlap.
        let v = make_native_single_breakend_with_sv_end(24_998_000, 25_010_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "BND-format DEL on non-coding transcript should have feature_truncation"
        );
        assert_eq!(tc.impact, Impact::HIGH);
    }

    #[test]
    fn test_native_single_breakend_with_sv_end_no_transcript_overlap() {
        // SV span does not overlap transcript at all → no feature_truncation.
        let tx = make_test_transcript(); // 25_000_000 - 25_006_000
        let v = make_native_single_breakend_with_sv_end(24_990_000, 24_995_000);
        let tc = calculate(&v, &tx, 5000, 5000);
        // Breakpoint at 24_990_000 is beyond the 5000bp upstream distance, so None.
        assert!(
            tc.is_none(),
            "BND-format DEL far from transcript should return None"
        );
    }

    /// Build a miRNA transcript matching the one in mod.rs tests.
    /// Single exon: genomic 1000-1100, cDNA 1-101.
    /// miRNA attribute "50-70" = cDNA positions 50-70 → genomic 1049-1069 (forward strand).
    fn make_mirna_transcript() -> vep_core::transcript::Transcript {
        use std::sync::Arc;
        use vep_core::transcript::{
            Attribute, Exon, ExonCoordMapper, MapperPair, TranscriptMapper, TranscriptVEFC,
        };

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

        vep_core::transcript::Transcript {
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

    /// Build a derived single-breakend variant (from a paired BND).
    /// Has is_single_breakend=true and mate_id set (unlike native single-breakend).
    fn make_derived_single_breakend(pos: u64) -> InputVariant {
        let mut v = make_bnd(pos);
        v.is_single_breakend = true;
        // mate_id stays set from make_bnd, which marks it as "derived" from a paired BND.
        v
    }

    #[test]
    fn test_derived_single_breakend_mirna_no_feature_truncation() {
        // Derived single-breakend at a miRNA exon within the mature miRNA region.
        // Perl emits only mature_miRNA_variant (no feature_truncation).
        let tx = make_mirna_transcript(); // exon 1000-1100, miRNA 1049-1069
        let v = make_derived_single_breakend(1060); // Inside mature miRNA region
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::MatureMirnaVariant),
            "Derived SB in miRNA exon should have mature_miRNA_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences.contains(&Consequence::FeatureTruncation),
            "Derived SB in miRNA exon should not have feature_truncation, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_derived_single_breakend_non_mirna_keeps_feature_truncation() {
        // Derived single-breakend at a non-miRNA non-coding exon.
        // Feature_truncation should still be present.
        let mut tx = make_mirna_transcript();
        tx.biotype = "lncRNA".into(); // Not miRNA → is_mature_mirna_sv returns false
        let v = make_derived_single_breakend(1060); // Inside exon
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "Derived SB in non-miRNA exon should keep feature_truncation, got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::NonCodingTranscriptExonVariant),
            "Should have non_coding_transcript_exon_variant, got: {:?}",
            tc.consequences
        );
    }

    /// Build an NMD variant of the standard test transcript.
    /// NMD transcripts have CDS (has_cds()=true) but biotype="nonsense_mediated_decay".
    fn make_nmd_transcript() -> Transcript {
        let mut tx = crate::test_helpers::make_test_transcript();
        tx.biotype = "nonsense_mediated_decay".into();
        tx
    }

    #[test]
    fn test_nmd_intron_gets_nmd_transcript_variant_not_non_coding() {
        // BND breakpoint in intron of NMD transcript.
        // Perl: feature_truncation,intron_variant,NMD_transcript_variant
        let tx = make_nmd_transcript();
        // Intron 1: 25_000_300 - 25_001_999
        let v = make_derived_single_breakend(25_001_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::NmdTranscriptVariant),
            "NMD intron should have NMD_transcript_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "NMD intron should not have non_coding_transcript_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "Should still have intron_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_nmd_coding_exon_gets_nmd_transcript_variant() {
        // BND breakpoint in coding exon of NMD transcript.
        // Perl: feature_truncation,coding_sequence_variant,NMD_transcript_variant
        let tx = make_nmd_transcript();
        // Exon 1 coding region: 25_000_050 - 25_000_299
        let v = make_derived_single_breakend(25_000_100);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::NmdTranscriptVariant),
            "NMD coding exon should have NMD_transcript_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences
                .contains(&Consequence::CodingSequenceVariant),
            "NMD coding exon should have coding_sequence_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "NMD should not have non_coding_transcript_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_nmd_utr_gets_nmd_transcript_variant() {
        // BND breakpoint in 5' UTR of NMD transcript.
        // NMD transcripts have UTR because they have CDS.
        let tx = make_nmd_transcript();
        // 5' UTR: 25_000_000 - 25_000_049 (before coding_region_start)
        let v = make_derived_single_breakend(25_000_020);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::NmdTranscriptVariant),
            "NMD UTR should have NMD_transcript_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences.contains(&Consequence::FivePrimeUtrVariant),
            "Should have 5_prime_UTR_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_paired_bnd_nmd_intron_gets_nmd_context() {
        // Paired BND with mate inside transcript, breakpoint in intron of NMD transcript.
        let tx = make_nmd_transcript();
        let mut v = make_bnd(25_001_000); // bp in intron
        v.mate_chr = Some("21".into());
        v.mate_pos = Some(25_002_100); // mate in exon 2
        v.is_single_breakend = false;
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::NmdTranscriptVariant),
            "Paired BND in NMD intron should have NMD_transcript_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences.contains(&Consequence::IntronVariant),
            "Should have intron_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_bnd_span_feature_truncation_nmd_gets_nmd_context() {
        // BND-format DEL spanning NMD transcript with intron-only overlap.
        // Should get feature_truncation + NMD_transcript_variant (not non_coding).
        let tx = make_nmd_transcript();
        // SV span extends beyond transcript with intron-only overlap
        let mut v = make_native_single_breakend(24_999_000);
        v.sv_end = Some(25_001_500); // Spans into intron 1 only (no exon overlap)
        v.sv_type = Some("DEL".into());
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::FeatureTruncation),
            "Should have feature_truncation, got: {:?}",
            tc.consequences
        );
        assert!(
            tc.consequences.contains(&Consequence::NmdTranscriptVariant),
            "NMD BND-format DEL should have NMD_transcript_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences
                .contains(&Consequence::NonCodingTranscriptVariant),
            "NMD should not have non_coding_transcript_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_bnd_upstream_nmd_no_nmd_context() {
        // BND breakpoint upstream of an NMD transcript should not get NMD_transcript_variant.
        // Perl only adds NMD context for body/splice overlaps, not distance-only annotations.
        let tx = make_nmd_transcript();
        // Position within upstream distance of transcript start (25_000_000 - 5000 = 24_995_000)
        let v = make_derived_single_breakend(24_996_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences.contains(&Consequence::UpstreamGeneVariant),
            "Should have upstream_gene_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences.contains(&Consequence::NmdTranscriptVariant),
            "Upstream BND should not have NMD_transcript_variant, got: {:?}",
            tc.consequences
        );
    }

    #[test]
    fn test_bnd_downstream_nmd_no_nmd_context() {
        // BND breakpoint downstream of an NMD transcript should not get NMD_transcript_variant.
        let tx = make_nmd_transcript();
        // Position well downstream of transcript end (25_006_000 + 5000 = 25_011_000 boundary)
        let v = make_derived_single_breakend(25_008_000);
        let tc = calculate(&v, &tx, 5000, 5000).unwrap();
        assert!(
            tc.consequences
                .contains(&Consequence::DownstreamGeneVariant),
            "Should have downstream_gene_variant, got: {:?}",
            tc.consequences
        );
        assert!(
            !tc.consequences.contains(&Consequence::NmdTranscriptVariant),
            "Downstream BND should not have NMD_transcript_variant, got: {:?}",
            tc.consequences
        );
    }
}
