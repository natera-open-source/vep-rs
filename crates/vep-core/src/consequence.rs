// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Consequence types and rankings from Ensembl Variation.
//!
//! All 41 Sequence Ontology (SO) consequence terms used by VEP, with their
//! exact ranks, impact classifications, and SO accessions matching the Ensembl VEP
//! implementation in `Bio::EnsEMBL::Variation::Utils::Constants`
//! (ensembl-variation). `OutputFactory.pm` is
//! `Bio/EnsEMBL/VEP/OutputFactory.pm` in ensembl-vep.

use std::borrow::Cow;
use std::fmt;
use std::sync::Arc;

use smallvec::SmallVec;

use crate::prediction::PredictionTerm;
pub use crate::transcript::SourceCache;
use crate::transcript::{file_order_key_suffix, BamEditStatus, Transcript};

/// Per-transcript consequence term list.
///
/// Inline capacity 8 covers the typical 1-8 SO terms a single transcript
/// consequence carries without touching the heap; it spills to a heap buffer
/// beyond that. `Consequence` is a 1-byte `Copy` enum, so the inline array is
/// cheap and the type derefs to `&[Consequence]` for all reads.
pub type ConsequenceList = SmallVec<[Consequence; 8]>;

/// Impact classification for consequences (HIGH, MODERATE, LOW, MODIFIER).
///
/// Matches Ensembl Variation's four-tier impact system. Used for filtering
/// and prioritizing consequences in downstream analysis.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize,
)]
pub enum Impact {
    HIGH,
    MODERATE,
    LOW,
    MODIFIER,
}

impl Impact {
    /// Static string form, matching the [`fmt::Display`] impl, allocation-free.
    ///
    /// Lets the output path write the IMPACT Extra field without materializing
    /// an owned `String` per transcript consequence.
    pub fn as_str(self) -> &'static str {
        match self {
            Impact::HIGH => "HIGH",
            Impact::MODERATE => "MODERATE",
            Impact::LOW => "LOW",
            Impact::MODIFIER => "MODIFIER",
        }
    }
}

impl fmt::Display for Impact {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// All 41 Sequence Ontology (SO) consequence terms used by VEP.
///
/// Each variant represents a specific molecular consequence type, ordered by
/// severity rank (lower rank = more severe). Ranks, impacts, and SO accessions
/// exactly match `Bio::EnsEMBL::Variation::Utils::Constants` from Ensembl
/// Variation.
///
/// See <https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html>
/// for the full consequence hierarchy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum Consequence {
    TranscriptAblation,
    SpliceAcceptorVariant,
    SpliceDonorVariant,
    StopGained,
    FrameshiftVariant,
    StopLost,
    StartLost,
    TranscriptAmplification,
    FeatureElongation,
    FeatureTruncation,
    InframeInsertion,
    InframeDeletion,
    MissenseVariant,
    ProteinAlteringVariant,
    SpliceDonor5thBaseVariant,
    SpliceRegionVariant,
    SpliceDonorRegionVariant,
    SplicePolypyrimidineTractVariant,
    IncompleteTerminalCodonVariant,
    StartRetainedVariant,
    StopRetainedVariant,
    SynonymousVariant,
    CodingSequenceVariant,
    MatureMirnaVariant,
    FivePrimeUtrVariant,
    ThreePrimeUtrVariant,
    NonCodingTranscriptExonVariant,
    IntronVariant,
    NmdTranscriptVariant,
    NonCodingTranscriptVariant,
    CodingTranscriptVariant,
    UpstreamGeneVariant,
    DownstreamGeneVariant,
    TfbsAblation,
    TfbsAmplification,
    TfBindingSiteVariant,
    RegulatoryRegionAblation,
    RegulatoryRegionAmplification,
    RegulatoryRegionVariant,
    IntergenicVariant,
    SequenceVariant,
}

impl Consequence {
    /// Severity rank (1 = most severe). Matches Ensembl Variation exactly.
    pub fn rank(self) -> u32 {
        match self {
            Consequence::TranscriptAblation => 1,
            Consequence::SpliceAcceptorVariant => 2,
            Consequence::SpliceDonorVariant => 3,
            Consequence::StopGained => 4,
            Consequence::FrameshiftVariant => 5,
            Consequence::StopLost => 6,
            Consequence::StartLost => 7,
            Consequence::TranscriptAmplification => 8,
            Consequence::FeatureElongation => 9,
            Consequence::FeatureTruncation => 10,
            Consequence::InframeInsertion => 11,
            Consequence::InframeDeletion => 12,
            Consequence::MissenseVariant => 13,
            Consequence::ProteinAlteringVariant => 14,
            Consequence::SpliceDonor5thBaseVariant => 15,
            Consequence::SpliceRegionVariant => 16,
            Consequence::SpliceDonorRegionVariant => 17,
            Consequence::SplicePolypyrimidineTractVariant => 18,
            Consequence::IncompleteTerminalCodonVariant => 19,
            Consequence::StartRetainedVariant => 20,
            Consequence::StopRetainedVariant => 21,
            Consequence::SynonymousVariant => 22,
            Consequence::CodingSequenceVariant => 23,
            Consequence::MatureMirnaVariant => 24,
            Consequence::FivePrimeUtrVariant => 25,
            Consequence::ThreePrimeUtrVariant => 26,
            Consequence::NonCodingTranscriptExonVariant => 27,
            Consequence::IntronVariant => 28,
            Consequence::NmdTranscriptVariant => 29,
            Consequence::NonCodingTranscriptVariant => 30,
            Consequence::CodingTranscriptVariant => 31,
            Consequence::UpstreamGeneVariant => 32,
            Consequence::DownstreamGeneVariant => 33,
            Consequence::TfbsAblation => 34,
            Consequence::TfbsAmplification => 35,
            Consequence::TfBindingSiteVariant => 36,
            Consequence::RegulatoryRegionAblation => 37,
            Consequence::RegulatoryRegionAmplification => 38,
            Consequence::RegulatoryRegionVariant => 39,
            Consequence::IntergenicVariant => 40,
            Consequence::SequenceVariant => 41,
        }
    }

    /// Impact classification for this consequence.
    pub fn impact(self) -> Impact {
        match self {
            Consequence::TranscriptAblation
            | Consequence::SpliceAcceptorVariant
            | Consequence::SpliceDonorVariant
            | Consequence::StopGained
            | Consequence::FrameshiftVariant
            | Consequence::StopLost
            | Consequence::StartLost
            | Consequence::TranscriptAmplification
            | Consequence::FeatureElongation
            | Consequence::FeatureTruncation => Impact::HIGH,

            Consequence::InframeInsertion
            | Consequence::InframeDeletion
            | Consequence::MissenseVariant
            | Consequence::ProteinAlteringVariant
            | Consequence::TfbsAblation => Impact::MODERATE,

            Consequence::SpliceDonor5thBaseVariant
            | Consequence::SpliceRegionVariant
            | Consequence::SpliceDonorRegionVariant
            | Consequence::SplicePolypyrimidineTractVariant
            | Consequence::IncompleteTerminalCodonVariant
            | Consequence::StartRetainedVariant
            | Consequence::StopRetainedVariant
            | Consequence::SynonymousVariant => Impact::LOW,

            Consequence::CodingSequenceVariant
            | Consequence::MatureMirnaVariant
            | Consequence::FivePrimeUtrVariant
            | Consequence::ThreePrimeUtrVariant
            | Consequence::NonCodingTranscriptExonVariant
            | Consequence::IntronVariant
            | Consequence::NmdTranscriptVariant
            | Consequence::NonCodingTranscriptVariant
            | Consequence::CodingTranscriptVariant
            | Consequence::UpstreamGeneVariant
            | Consequence::DownstreamGeneVariant
            | Consequence::TfbsAmplification
            | Consequence::TfBindingSiteVariant
            | Consequence::RegulatoryRegionAblation
            | Consequence::RegulatoryRegionAmplification
            | Consequence::RegulatoryRegionVariant
            | Consequence::IntergenicVariant
            | Consequence::SequenceVariant => Impact::MODIFIER,
        }
    }

    /// SO term string exactly matching VEP output format.
    pub fn so_term(self) -> &'static str {
        match self {
            Consequence::TranscriptAblation => "transcript_ablation",
            Consequence::SpliceAcceptorVariant => "splice_acceptor_variant",
            Consequence::SpliceDonorVariant => "splice_donor_variant",
            Consequence::StopGained => "stop_gained",
            Consequence::FrameshiftVariant => "frameshift_variant",
            Consequence::StopLost => "stop_lost",
            Consequence::StartLost => "start_lost",
            Consequence::TranscriptAmplification => "transcript_amplification",
            Consequence::FeatureElongation => "feature_elongation",
            Consequence::FeatureTruncation => "feature_truncation",
            Consequence::InframeInsertion => "inframe_insertion",
            Consequence::InframeDeletion => "inframe_deletion",
            Consequence::MissenseVariant => "missense_variant",
            Consequence::ProteinAlteringVariant => "protein_altering_variant",
            Consequence::SpliceDonor5thBaseVariant => "splice_donor_5th_base_variant",
            Consequence::SpliceRegionVariant => "splice_region_variant",
            Consequence::SpliceDonorRegionVariant => "splice_donor_region_variant",
            Consequence::SplicePolypyrimidineTractVariant => "splice_polypyrimidine_tract_variant",
            Consequence::IncompleteTerminalCodonVariant => "incomplete_terminal_codon_variant",
            Consequence::StartRetainedVariant => "start_retained_variant",
            Consequence::StopRetainedVariant => "stop_retained_variant",
            Consequence::SynonymousVariant => "synonymous_variant",
            Consequence::CodingSequenceVariant => "coding_sequence_variant",
            Consequence::MatureMirnaVariant => "mature_miRNA_variant",
            Consequence::FivePrimeUtrVariant => "5_prime_UTR_variant",
            Consequence::ThreePrimeUtrVariant => "3_prime_UTR_variant",
            Consequence::NonCodingTranscriptExonVariant => "non_coding_transcript_exon_variant",
            Consequence::IntronVariant => "intron_variant",
            Consequence::NmdTranscriptVariant => "NMD_transcript_variant",
            Consequence::NonCodingTranscriptVariant => "non_coding_transcript_variant",
            Consequence::CodingTranscriptVariant => "coding_transcript_variant",
            Consequence::UpstreamGeneVariant => "upstream_gene_variant",
            Consequence::DownstreamGeneVariant => "downstream_gene_variant",
            Consequence::TfbsAblation => "TFBS_ablation",
            Consequence::TfbsAmplification => "TFBS_amplification",
            Consequence::TfBindingSiteVariant => "TF_binding_site_variant",
            Consequence::RegulatoryRegionAblation => "regulatory_region_ablation",
            Consequence::RegulatoryRegionAmplification => "regulatory_region_amplification",
            Consequence::RegulatoryRegionVariant => "regulatory_region_variant",
            Consequence::IntergenicVariant => "intergenic_variant",
            Consequence::SequenceVariant => "sequence_variant",
        }
    }

    /// SO accession identifier.
    pub fn so_accession(self) -> &'static str {
        match self {
            Consequence::TranscriptAblation => "SO:0001893",
            Consequence::SpliceAcceptorVariant => "SO:0001574",
            Consequence::SpliceDonorVariant => "SO:0001575",
            Consequence::StopGained => "SO:0001587",
            Consequence::FrameshiftVariant => "SO:0001589",
            Consequence::StopLost => "SO:0001578",
            Consequence::StartLost => "SO:0002012",
            Consequence::TranscriptAmplification => "SO:0001889",
            Consequence::FeatureElongation => "SO:0001907",
            Consequence::FeatureTruncation => "SO:0001906",
            Consequence::InframeInsertion => "SO:0001821",
            Consequence::InframeDeletion => "SO:0001822",
            Consequence::MissenseVariant => "SO:0001583",
            Consequence::ProteinAlteringVariant => "SO:0001818",
            Consequence::SpliceDonor5thBaseVariant => "SO:0001787",
            Consequence::SpliceRegionVariant => "SO:0001630",
            Consequence::SpliceDonorRegionVariant => "SO:0002170",
            Consequence::SplicePolypyrimidineTractVariant => "SO:0002169",
            Consequence::IncompleteTerminalCodonVariant => "SO:0001626",
            Consequence::StartRetainedVariant => "SO:0002019",
            Consequence::StopRetainedVariant => "SO:0001567",
            Consequence::SynonymousVariant => "SO:0001819",
            Consequence::CodingSequenceVariant => "SO:0001580",
            Consequence::MatureMirnaVariant => "SO:0001620",
            Consequence::FivePrimeUtrVariant => "SO:0001623",
            Consequence::ThreePrimeUtrVariant => "SO:0001624",
            Consequence::NonCodingTranscriptExonVariant => "SO:0001792",
            Consequence::IntronVariant => "SO:0001627",
            Consequence::NmdTranscriptVariant => "SO:0001621",
            Consequence::NonCodingTranscriptVariant => "SO:0001619",
            Consequence::CodingTranscriptVariant => "SO:0001968",
            Consequence::UpstreamGeneVariant => "SO:0001631",
            Consequence::DownstreamGeneVariant => "SO:0001632",
            Consequence::TfbsAblation => "SO:0001895",
            Consequence::TfbsAmplification => "SO:0001892",
            Consequence::TfBindingSiteVariant => "SO:0001782",
            Consequence::RegulatoryRegionAblation => "SO:0001894",
            Consequence::RegulatoryRegionAmplification => "SO:0001891",
            Consequence::RegulatoryRegionVariant => "SO:0001566",
            Consequence::IntergenicVariant => "SO:0001628",
            Consequence::SequenceVariant => "SO:0001060",
        }
    }

    /// Parse a consequence from its SO term string.
    pub fn from_so_term(term: &str) -> Option<Consequence> {
        match term {
            "transcript_ablation" => Some(Consequence::TranscriptAblation),
            "splice_acceptor_variant" => Some(Consequence::SpliceAcceptorVariant),
            "splice_donor_variant" => Some(Consequence::SpliceDonorVariant),
            "stop_gained" => Some(Consequence::StopGained),
            "frameshift_variant" => Some(Consequence::FrameshiftVariant),
            "stop_lost" => Some(Consequence::StopLost),
            "start_lost" => Some(Consequence::StartLost),
            "transcript_amplification" => Some(Consequence::TranscriptAmplification),
            "feature_elongation" => Some(Consequence::FeatureElongation),
            "feature_truncation" => Some(Consequence::FeatureTruncation),
            "inframe_insertion" => Some(Consequence::InframeInsertion),
            "inframe_deletion" => Some(Consequence::InframeDeletion),
            "missense_variant" => Some(Consequence::MissenseVariant),
            "protein_altering_variant" => Some(Consequence::ProteinAlteringVariant),
            "splice_donor_5th_base_variant" => Some(Consequence::SpliceDonor5thBaseVariant),
            "splice_region_variant" => Some(Consequence::SpliceRegionVariant),
            "splice_donor_region_variant" => Some(Consequence::SpliceDonorRegionVariant),
            "splice_polypyrimidine_tract_variant" => {
                Some(Consequence::SplicePolypyrimidineTractVariant)
            }
            "incomplete_terminal_codon_variant" => {
                Some(Consequence::IncompleteTerminalCodonVariant)
            }
            "start_retained_variant" => Some(Consequence::StartRetainedVariant),
            "stop_retained_variant" => Some(Consequence::StopRetainedVariant),
            "synonymous_variant" => Some(Consequence::SynonymousVariant),
            "coding_sequence_variant" => Some(Consequence::CodingSequenceVariant),
            "mature_miRNA_variant" => Some(Consequence::MatureMirnaVariant),
            "5_prime_UTR_variant" => Some(Consequence::FivePrimeUtrVariant),
            "3_prime_UTR_variant" => Some(Consequence::ThreePrimeUtrVariant),
            "non_coding_transcript_exon_variant" => {
                Some(Consequence::NonCodingTranscriptExonVariant)
            }
            "intron_variant" => Some(Consequence::IntronVariant),
            "NMD_transcript_variant" => Some(Consequence::NmdTranscriptVariant),
            "non_coding_transcript_variant" => Some(Consequence::NonCodingTranscriptVariant),
            "coding_transcript_variant" => Some(Consequence::CodingTranscriptVariant),
            "upstream_gene_variant" => Some(Consequence::UpstreamGeneVariant),
            "downstream_gene_variant" => Some(Consequence::DownstreamGeneVariant),
            "TFBS_ablation" => Some(Consequence::TfbsAblation),
            "TFBS_amplification" => Some(Consequence::TfbsAmplification),
            "TF_binding_site_variant" => Some(Consequence::TfBindingSiteVariant),
            "regulatory_region_ablation" => Some(Consequence::RegulatoryRegionAblation),
            "regulatory_region_amplification" => Some(Consequence::RegulatoryRegionAmplification),
            "regulatory_region_variant" => Some(Consequence::RegulatoryRegionVariant),
            "intergenic_variant" => Some(Consequence::IntergenicVariant),
            "sequence_variant" => Some(Consequence::SequenceVariant),
            _ => None,
        }
    }

    /// Return the most severe consequence from a slice.
    pub fn most_severe(consequences: &[Consequence]) -> Option<Consequence> {
        consequences.iter().copied().min_by_key(|c| c.rank())
    }

    /// Every term in rank order, so `ALL[i].rank() == i + 1` and a term's
    /// position in the list is `term as usize`.
    pub const ALL: &'static [Consequence] = &[
        Consequence::TranscriptAblation,
        Consequence::SpliceAcceptorVariant,
        Consequence::SpliceDonorVariant,
        Consequence::StopGained,
        Consequence::FrameshiftVariant,
        Consequence::StopLost,
        Consequence::StartLost,
        Consequence::TranscriptAmplification,
        Consequence::FeatureElongation,
        Consequence::FeatureTruncation,
        Consequence::InframeInsertion,
        Consequence::InframeDeletion,
        Consequence::MissenseVariant,
        Consequence::ProteinAlteringVariant,
        Consequence::SpliceDonor5thBaseVariant,
        Consequence::SpliceRegionVariant,
        Consequence::SpliceDonorRegionVariant,
        Consequence::SplicePolypyrimidineTractVariant,
        Consequence::IncompleteTerminalCodonVariant,
        Consequence::StartRetainedVariant,
        Consequence::StopRetainedVariant,
        Consequence::SynonymousVariant,
        Consequence::CodingSequenceVariant,
        Consequence::MatureMirnaVariant,
        Consequence::FivePrimeUtrVariant,
        Consequence::ThreePrimeUtrVariant,
        Consequence::NonCodingTranscriptExonVariant,
        Consequence::IntronVariant,
        Consequence::NmdTranscriptVariant,
        Consequence::NonCodingTranscriptVariant,
        Consequence::CodingTranscriptVariant,
        Consequence::UpstreamGeneVariant,
        Consequence::DownstreamGeneVariant,
        Consequence::TfbsAblation,
        Consequence::TfbsAmplification,
        Consequence::TfBindingSiteVariant,
        Consequence::RegulatoryRegionAblation,
        Consequence::RegulatoryRegionAmplification,
        Consequence::RegulatoryRegionVariant,
        Consequence::IntergenicVariant,
        Consequence::SequenceVariant,
    ];

    /// Whether the term's `OVERLAP_CONSEQUENCES` entry carries
    /// `include => {coding => 1}` (`Utils/Constants.pm`): the thirteen terms
    /// the `Coding consequences` chart of the run statistics counts.
    /// `coding_transcript_variant` is not one.
    pub fn is_coding_include(self) -> bool {
        matches!(
            self,
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
}

impl fmt::Display for Consequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.so_term())
    }
}

impl PartialOrd for Consequence {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Consequence {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.rank().cmp(&other.rank())
    }
}

/// Transcript structural context for the built-in plugins that read the
/// transcript model in the plugin phase (LoFTEE among them), which the
/// summarized [`TranscriptConsequence`] fields do not expose.
///
/// The per-consequence `exon`/`intron` fields are only `"n/total"` summary
/// strings; this context holds what LoFTEE computes SMALL_INTRON,
/// NON_CAN_SPLICE, GC_TO_GT_DONOR and the GERP-weighted END_TRUNC distance
/// from: the transcript's exons and introns in rank order (the introns derived
/// from the exons when the transcript lists none), its strand, coding bounds
/// and `cds_*_NF` flags.
///
/// The code that builds a [`TranscriptConsequence`] attaches this struct while
/// it holds the transcript, on the rows
/// `vep_effects::EffectsConfig::loftee_context` names, which the `vep` binary
/// chooses from the loaded plugins. On every other row the field is `None` and
/// absent from the serialization.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct LofteeContext {
    /// All exons of the transcript, in rank order.
    pub exons: Vec<crate::transcript::Exon>,
    /// All introns of the transcript, in rank order. May be derived from
    /// `exons` when the transcript's own intron list is not populated.
    pub introns: Vec<crate::transcript::Intron>,
    /// Transcript strand (1 forward, -1 reverse).
    pub strand: i8,
    /// Genomic coding-region start (1-based), if coding.
    pub coding_region_start: Option<u64>,
    /// Genomic coding-region end (1-based), if coding.
    pub coding_region_end: Option<u64>,
    /// Genomic translation start (1-based), if coding.
    pub translation_start: Option<u64>,
    /// Genomic translation end (1-based), if coding.
    pub translation_end: Option<u64>,
    /// cDNA coding start (1-based offset within cDNA), if coding.
    pub cdna_coding_start: Option<u64>,
    /// cDNA coding end (1-based offset within cDNA), if coding.
    pub cdna_coding_end: Option<u64>,
    /// `cds_start_NF` annotation present (start codon not found / incomplete CDS).
    pub cds_start_nf: bool,
    /// `cds_end_NF` annotation present (stop codon not found / incomplete CDS).
    pub cds_end_nf: bool,
}

impl LofteeContext {
    /// Build the LoFTEE structural context from a source transcript.
    ///
    /// Clones the exon list (cheap relative to the per-variant work, and done
    /// only for a row that is to carry the context). When the
    /// transcript's own intron list is empty, derives intron coordinates from
    /// the rank-sorted exons (the gaps between consecutive exons), matching
    /// LoFTEE's `transcript->get_all_Introns`.
    pub fn from_transcript(transcript: &crate::transcript::Transcript) -> Self {
        let mut exons = transcript.exons.clone();
        exons.sort_by_key(|e| e.rank);

        let introns = if !transcript.introns.is_empty() {
            let mut introns = transcript.introns.clone();
            introns.sort_by_key(|i| i.rank);
            introns
        } else {
            derive_introns_from_exons(&exons)
        };

        let has_flag = |name: &str| transcript.flags.iter().any(|f| f.as_str() == name);

        Self {
            exons,
            introns,
            strand: transcript.strand.as_i8(),
            coding_region_start: transcript.coding_region_start,
            coding_region_end: transcript.coding_region_end,
            translation_start: transcript.translation_start,
            translation_end: transcript.translation_end,
            cdna_coding_start: transcript.cdna_coding_start,
            cdna_coding_end: transcript.cdna_coding_end,
            cds_start_nf: has_flag("cds_start_NF"),
            cds_end_nf: has_flag("cds_end_NF"),
        }
    }
}

/// Derive intron coordinates from rank-sorted exons.
///
/// Exons are ordered by genomic coordinate (independent of strand); the intron
/// between two genomically adjacent exons spans `[prev.end + 1, next.start - 1]`.
/// Intron rank follows exon rank order. Matches the gaps Ensembl's
/// `get_all_Introns` reports.
fn derive_introns_from_exons(
    sorted_by_rank: &[crate::transcript::Exon],
) -> Vec<crate::transcript::Intron> {
    if sorted_by_rank.len() < 2 {
        return Vec::new();
    }
    let mut by_pos: Vec<&crate::transcript::Exon> = sorted_by_rank.iter().collect();
    by_pos.sort_by_key(|e| e.start);

    let mut introns = Vec::with_capacity(by_pos.len() - 1);
    for (idx, pair) in by_pos.windows(2).enumerate() {
        let (left, right) = (pair[0], pair[1]);
        if right.start > left.end + 1 {
            introns.push(crate::transcript::Intron {
                start: left.end + 1,
                end: right.start - 1,
                rank: (idx + 1) as u32,
            });
        }
    }
    introns
}

/// Consequence annotation for a specific transcript-variant overlap.
///
/// Produced by `vep_effects::calculate_consequences()` for each transcript
/// that a variant overlaps. Contains the assigned SO consequence terms,
/// positional information (cDNA, CDS, protein coordinates), codon/amino acid
/// changes, and a shared handle on the source transcript, whose identifiers,
/// symbol, biotype, flags, cross-references and pick criteria the accessors of
/// the same names read ([`Self::transcript_id`], [`Self::gene_symbol`],
/// [`Self::biotype`], [`Self::mane_select`], ...). Multiple instances are
/// attached to each [`InputVariant`](crate::variant::InputVariant).
///
/// The serde document is flat: the transcript's values appear under their
/// accessor names, in the order [`TRANSCRIPT_CONSEQUENCE_KEYS`] lists, among
/// the row's own fields. `transcript`, `feature_id`, `feature_biotype`,
/// `picked` and `protein_id_display` have no key of their own (`picked` is the
/// `PICK` entry of `flags`, `protein_id_display` shapes `protein_id`), nor
/// does a field whose doc marks it never serialized, and the keys the list
/// marks optional appear only when populated.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct TranscriptConsequence {
    /// The source transcript, shared with the index that holds it: every
    /// per-transcript value the row prints is read from it. `None` on a row
    /// against a regulatory or motif feature and on a row built without a
    /// feature, whose transcript accessors read `feature_id`,
    /// `feature_biotype` or nothing.
    pub transcript: Option<Arc<Transcript>>,
    /// The stable id of a regulatory or motif feature's row, which has no
    /// transcript to read one from ([`Self::transcript_id`]); `None` on a
    /// transcript's row.
    pub feature_id: Option<Arc<str>>,
    /// The biotype of a regulatory feature's row ([`Self::biotype`]); `None`
    /// on a transcript's row and on a structural variant's regulatory row,
    /// which prints none.
    pub feature_biotype: Option<Arc<str>>,
    /// Consequence terms assigned.
    pub consequences: ConsequenceList,
    /// Impact of the most severe consequence.
    pub impact: Impact,
    /// Whether this is the canonical transcript.
    pub canonical: bool,
    /// cDNA position (e.g., "428").
    pub cdna_position: Option<String>,
    /// CDS position (e.g., "301").
    pub cds_position: Option<String>,
    /// Protein position (e.g., "101").
    pub protein_position: Option<String>,
    /// The totals `--total_length` appends to the three positions: the
    /// transcript's cDNA length where the variant overlaps an exon, the coding
    /// sequence length and the peptide length where it lies in the coding
    /// region and the cache carries the sequence. Set whenever the column's
    /// rule holds, position or not, so a `-` column prints `-/<total>`.
    pub cdna_total: Option<u32>,
    pub cds_total: Option<u32>,
    pub protein_total: Option<u32>,
    /// The three positions of the allele's own span where the columns print
    /// another span's: for an allele of a `--minimal` split the columns print the
    /// record's span, and these are the trimmed allele's own, which the VCF and
    /// JSON formats print (`OutputFactory/VCF.pm:305-308` and
    /// `OutputFactory/JSON.pm:217-220` clear the cached positions before the
    /// rejoin, so they are recomputed on the trimmed feature and kept). `None`
    /// when the columns are the allele's own.
    pub own_span_positions: Option<Box<OwnSpanPositions>>,
    /// Amino acid change (e.g., "A/V").
    pub amino_acids: Option<String>,
    /// Codon change (e.g., "gCc/gTc").
    pub codons: Option<String>,
    /// How the row prints the transcript's protein id ([`Self::protein_id`]).
    pub protein_id_display: ProteinIdDisplay,
    /// Distance to transcript (for upstream/downstream).
    pub distance: Option<u64>,
    /// Strand of the transcript (1 or -1).
    pub strand: i8,
    /// Exon number / total (e.g., "3/5").
    pub exon: Option<String>,
    /// Intron number / total (e.g., "2/4").
    pub intron: Option<String>,
    /// HGVS coding notation.
    pub hgvsc: Option<String>,
    /// HGVS protein notation.
    pub hgvsp: Option<String>,
    /// Bases the HGVS notations were shifted 3' along the transcript, signed by
    /// the transcript strand; `None` when nothing shifted.
    pub hgvs_offset: Option<i64>,
    /// Whether `hgvs_offset` counts a move of the HGVSc alone, on a record
    /// that does not shift: every other field of the row, and what a plugin
    /// reads of it, stays where the record is written. Never serialized.
    pub hgvs_offset_of_hgvsc_alone: bool,
    /// The intron offsets of the HGVSc start and end coordinates as
    /// `hgvs_transcript` keeps them (`_hgvs_intron_start_offset`,
    /// `_hgvs_intron_end_offset`, TranscriptVariationAllele.pm:1530-1545): the
    /// signed distance after the exon coordinate of each cDNA position, read
    /// before the two positions are ordered, 0 for an exonic position. `None`
    /// when the row has no HGVSc. Never serialized.
    pub hgvs_intron_offsets: Option<(i32, i32)>,
    /// Bases the insertion or deletion was moved 3' along the transcript before
    /// its terms and positions were computed (`--shift_3prime 1`); 0 when nothing
    /// moved, every other row included.
    pub shift_length: u64,
    /// The `Allele` the row prints under `--shift_3prime 1` when it is not the
    /// allele as the record wrote it: the inserted bases rotated along the
    /// shift on a row whose shifted span has both translation coordinates,
    /// since `codon` rotates the allele's sequence there (`shift_feature_seqs`,
    /// TranscriptVariationAllele.pm:820, 1330-1347) and the row prints that
    /// sequence (OutputFactory.pm:1347); `None` for every other row.
    pub shifted_allele: Option<Arc<str>>,
    /// The reference allele the row was analysed against under
    /// `--use_transcript_ref`, relative to the variant's strand: the transcript's
    /// where the variant maps inside its exons, the input's elsewhere (`USED_REF`,
    /// OutputFactory.pm:1783-1788); `None` without the flag.
    pub used_ref: Option<Arc<str>>,
    /// SIFT prediction and score.
    pub sift: Option<String>,
    /// PolyPhen prediction and score.
    pub polyphen: Option<String>,
    /// The qualitative term behind `sift`, kept for the run statistics, which
    /// count terms where the column may print a score alone. Never serialized.
    pub sift_term: Option<PredictionTerm>,
    /// The qualitative term behind `polyphen`; see `sift_term`. Never serialized.
    pub polyphen_term: Option<PredictionTerm>,
    /// The first residue the variant maps to, Ensembl VEP's `translation_start`
    /// before the `Protein_position` column orders the span's ends (an
    /// insertion's start is one past its end, so the column opens with the
    /// end); `None` where no residue is mapped. Never serialized.
    pub protein_start: Option<u32>,
    /// Overlapping protein domains.
    pub domains: Vec<(String, String)>,
    /// Feature type (Transcript, RegulatoryFeature, MotifFeature).
    pub feature_type: FeatureType,
    /// Whether `--flag_pick`, `--flag_pick_allele` or `--flag_pick_allele_gene`
    /// picked the row: the `PICK` column, and the `PICK` entry after the
    /// transcript's flags in the serialized `flags`.
    pub picked: bool,
    /// The miRNA secondary-structure elements the variant overlaps on a
    /// `miRNA` transcript, as the SO terms `miRNA_loop` and `miRNA_stem`: the
    /// `miRNA` column under `--mirna`. The terms come sorted, one per distinct
    /// structure character, so a span over both bracket kinds lists
    /// `miRNA_stem` twice, and a span reaching one position past the
    /// structure's end lists an empty term. Empty on every other transcript
    /// and when the variant lies outside the structure.
    pub mirna: Vec<String>,
    /// Whether the transcript's model was edited against the cache's alignment
    /// file ([`Transcript::bam_edit_status`]): the `BAM_EDIT` column; `None` on
    /// a transcript no alignment touched.
    pub bam_edit_status: Option<BamEditStatus>,
    /// The length change of the edits an `NM_` or `XM_` model's alignment
    /// inserted or deleted upstream of the allele, Ensembl VEP's
    /// `refseq_misalignment_offset` as `hgvs_transcript` sets it
    /// (`REFSEQ_OFFSET`); the HGVSc is numbered against the edited model and
    /// moves by nothing. `None` when no offset was computed or it is 0, which
    /// prints nothing (OutputFactory.pm:1778).
    pub refseq_offset: Option<i64>,
    /// The alternate allele equals the bases of the transcript's first
    /// `_rna_edit` on the codon path (`invalid_alleles`,
    /// TranscriptVariationAllele.pm:826-837): the row the
    /// `Transcript-assembly mismatch` warning names. Never serialized.
    pub transcript_assembly_mismatch: bool,
    /// Plugin annotation data (transcript-consequence-level).
    /// Keyed by field name (e.g., "REVEL"), value is the annotation string.
    pub plugin_data: indexmap::IndexMap<String, String>,
    /// The keys of `plugin_data` whose empty value on this row is Ensembl VEP's
    /// empty string, an empty tab cell and `""` in JSON, rather than its
    /// undefined value, `-` and no key: a plugin reading a table names the
    /// columns a row has and leaves empty (TranscriptAnnotator.pm:245-250).
    /// Empty on a row with no such key, and never serialized.
    pub plugin_empty_strings: Vec<String>,
    /// Genomic span of the feature (1-based, inclusive), for the structural
    /// variant overlap fields; both zero when unknown.
    pub feature_start: u64,
    pub feature_end: u64,
    /// Transcript structural context for the built-in plugins that read the
    /// transcript model (LoFTEE among them). Attached by the code that builds
    /// the row, on the rows `vep_effects::EffectsConfig::loftee_context` names;
    /// `None` otherwise, and absent from the serialization unless set.
    /// See [`LofteeContext`].
    pub loftee_ctx: Option<Box<LofteeContext>>,
    /// The transcript's exon stable ids in rank order, one per exon of
    /// `loftee_ctx`, for the built-in plugins that print them (NearestExonJB);
    /// attached with the context when the cache was loaded with the ids
    /// (`Transcript::exon_ids`), `None` otherwise. Never serialized.
    pub exon_ids: Option<Arc<[Arc<str>]>>,
    /// The transcript's coding, 3' UTR and protein sequences for the built-in
    /// plugins that translate a changed coding sequence (Downstream). Attached
    /// by the code that builds the row, on the rows
    /// `vep_effects::EffectsConfig::context_sequences` names, `None` on every
    /// other row. Never serialized.
    pub transcript_sequences: Option<Box<TranscriptSequences>>,
    /// The motif fields of a `MotifFeature` row of a sequence variant
    /// (`MotifFeatureVariationAllele_to_output_hash`, OutputFactory.pm:1919-1951):
    /// `None` on every other row, a structural variant's motif row included.
    pub motif: Option<Box<MotifFields>>,
    /// The `CELL_TYPE` list of a `RegulatoryFeature` or `MotifFeature` row under
    /// `--cell_type` (`get_cell_types`, OutputFactory.pm:1968-1978): for each
    /// requested epigenome, in the order requested, the feature's activity in
    /// it, an epigenome the feature's map lacks skipped. Empty on every motif
    /// row, whose features carry no activity map; `None` on a transcript or
    /// intergenic row and on every row of a run without the flag.
    pub cell_type: Option<Box<[CellTypeActivity]>>,
    /// The source transcript's [`Transcript::id_rank`]: the position of its
    /// Ensembl VEP sort key ([`Self::ensembl_vep_key`]) in the byte order of
    /// its chromosome's keys, so two
    /// rows of one chromosome compare by rank as Ensembl VEP orders them; rows
    /// of different chromosomes, an interchromosomal breakend's mate rows among
    /// them, do not. A `MotifFeature` row carries the rank of its motif's
    /// database id, as a string, among its chromosome's motifs: Ensembl VEP
    /// stores motif overlaps under that id and sorts the keys as strings
    /// (`VariationFeature.pm:937`, `:750`). A structural variant's
    /// `RegulatoryFeature` or `MotifFeature` row carries instead the rank of
    /// its feature in the order the features were met, the order Ensembl VEP
    /// lists such rows in (`StructuralVariationFeature.pm:495-509`): source by
    /// source then in file order from the files of `--regulatory_gff` or a
    /// cache baked from them, the cache's listing order from a converted cache.
    /// Never serialized; 0 for a consequence built without a loaded
    /// transcript or a motif without an id.
    pub id_rank: u32,
    /// The database id of a `MotifFeature` row's motif, the key Ensembl VEP
    /// stores the motif's overlap under (`MotifFeatureVariation.pm:70-76`,
    /// `VariationFeature.pm:937`) and the part of [`Self::ensembl_vep_key`]
    /// that is not a stable id. Never serialized; `None` on every other row
    /// and on a motif from a source without one.
    pub motif_db_id: Option<u64>,
    /// Whether the record's span, for an allele of a `--minimal` split
    /// ([`InputVariant::original_start`](crate::variant::InputVariant::original_start)
    /// to `original_end`), has a CDS coordinate on the transcript: the
    /// `within_cds` the `--coding_only` filter asks of a rejoined row in the
    /// default and tab formats, whose factories clear the coordinates cached on
    /// the allele's own span after the rejoin and map the record's span in their
    /// place (OutputFactory.pm:297-323 and 490-502; BaseTranscriptVariation.pm:505-526
    /// and 662-665). `false` for every other row. Never serialized.
    pub record_span_in_cds: bool,
}

/// How a row prints the protein id of its transcript
/// ([`TranscriptConsequence::protein_id`]).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum ProteinIdDisplay {
    /// The id as the cache carries it.
    #[default]
    AsCached,
    /// No id: a breakend row against a transcript without a coding region.
    Hidden,
    /// The id with `.<version>` appended, the translation's version under
    /// `--protein_version` on an id that carries none.
    Versioned(u32),
}

/// The keys of a [`TranscriptConsequence`] document in the order they are
/// written, each with whether it is written only when populated.
pub const TRANSCRIPT_CONSEQUENCE_KEYS: &[(&str, bool)] = &[
    ("transcript_id", false),
    ("gene_id", false),
    ("gene_symbol", false),
    ("gene_symbol_source", false),
    ("hgnc_id", false),
    ("consequences", false),
    ("impact", false),
    ("biotype", false),
    ("canonical", false),
    ("cdna_position", false),
    ("cds_position", false),
    ("protein_position", false),
    ("cdna_total", true),
    ("cds_total", true),
    ("protein_total", true),
    ("own_span_positions", true),
    ("amino_acids", false),
    ("codons", false),
    ("protein_id", false),
    ("distance", false),
    ("strand", false),
    ("exon", false),
    ("intron", false),
    ("hgvsc", false),
    ("hgvsp", false),
    ("hgvs_offset", false),
    ("shift_length", false),
    ("shifted_allele", true),
    ("used_ref", true),
    ("sift", false),
    ("polyphen", false),
    ("domains", false),
    ("feature_type", false),
    ("flags", false),
    ("tsl", false),
    ("mane_select", false),
    ("mane_plus_clinical", false),
    ("appris", false),
    ("ccds", false),
    ("gencode_primary", false),
    ("swissprot", false),
    ("trembl", false),
    ("uniparc", false),
    ("uniprot_isoform", false),
    ("refseq", false),
    ("transcript_version", false),
    ("gene_version", false),
    ("gene_phenotype", false),
    ("mirna", false),
    ("cds_length", true),
    ("cdna_length", false),
    ("source_cache", true),
    ("bam_edit_status", true),
    ("refseq_match", true),
    ("refseq_offset", true),
    ("plugin_data", true),
    ("feature_start", false),
    ("feature_end", false),
    ("loftee_ctx", true),
    ("motif", true),
    ("cell_type", true),
];

impl serde::Serialize for TranscriptConsequence {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;

        let cds_length = self.cds_length();
        let source_cache = self.source_cache();
        let refseq_match = self.refseq_match();
        let written = TRANSCRIPT_CONSEQUENCE_KEYS
            .iter()
            .filter(|(_, optional)| !optional)
            .count()
            + usize::from(self.cdna_total.is_some())
            + usize::from(self.cds_total.is_some())
            + usize::from(self.protein_total.is_some())
            + usize::from(self.own_span_positions.is_some())
            + usize::from(self.shifted_allele.is_some())
            + usize::from(self.used_ref.is_some())
            + usize::from(cds_length.is_some())
            + usize::from(source_cache.is_some())
            + usize::from(self.bam_edit_status.is_some())
            + usize::from(refseq_match.is_some())
            + usize::from(self.refseq_offset.is_some())
            + usize::from(!self.plugin_data.is_empty())
            + usize::from(self.loftee_ctx.is_some())
            + usize::from(self.motif.is_some())
            + usize::from(self.cell_type.is_some());
        let mut state = serializer.serialize_struct("TranscriptConsequence", written)?;
        state.serialize_field("transcript_id", self.transcript_id())?;
        state.serialize_field("gene_id", self.gene_id())?;
        state.serialize_field("gene_symbol", &self.gene_symbol())?;
        state.serialize_field("gene_symbol_source", &self.gene_symbol_source())?;
        state.serialize_field("hgnc_id", &self.hgnc_id())?;
        state.serialize_field("consequences", &self.consequences)?;
        state.serialize_field("impact", &self.impact)?;
        state.serialize_field("biotype", &self.biotype())?;
        state.serialize_field("canonical", &self.canonical)?;
        state.serialize_field("cdna_position", &self.cdna_position)?;
        state.serialize_field("cds_position", &self.cds_position)?;
        state.serialize_field("protein_position", &self.protein_position)?;
        optional_field(&mut state, "cdna_total", &self.cdna_total)?;
        optional_field(&mut state, "cds_total", &self.cds_total)?;
        optional_field(&mut state, "protein_total", &self.protein_total)?;
        optional_field(&mut state, "own_span_positions", &self.own_span_positions)?;
        state.serialize_field("amino_acids", &self.amino_acids)?;
        state.serialize_field("codons", &self.codons)?;
        state.serialize_field("protein_id", &self.protein_id())?;
        state.serialize_field("distance", &self.distance)?;
        state.serialize_field("strand", &self.strand)?;
        state.serialize_field("exon", &self.exon)?;
        state.serialize_field("intron", &self.intron)?;
        state.serialize_field("hgvsc", &self.hgvsc)?;
        state.serialize_field("hgvsp", &self.hgvsp)?;
        state.serialize_field("hgvs_offset", &self.hgvs_offset)?;
        state.serialize_field("shift_length", &self.shift_length)?;
        optional_field(&mut state, "shifted_allele", &self.shifted_allele)?;
        optional_field(&mut state, "used_ref", &self.used_ref)?;
        state.serialize_field("sift", &self.sift)?;
        state.serialize_field("polyphen", &self.polyphen)?;
        state.serialize_field("domains", &self.domains)?;
        state.serialize_field("feature_type", &self.feature_type)?;
        state.serialize_field("flags", &FlagsDocument(self))?;
        state.serialize_field("tsl", &TslDocument(self.tsl()))?;
        state.serialize_field("mane_select", &self.mane_select())?;
        state.serialize_field("mane_plus_clinical", &self.mane_plus_clinical())?;
        state.serialize_field("appris", &self.appris())?;
        state.serialize_field("ccds", &self.ccds())?;
        state.serialize_field("gencode_primary", &self.gencode_primary())?;
        state.serialize_field("swissprot", &self.swissprot())?;
        state.serialize_field("trembl", &self.trembl())?;
        state.serialize_field("uniparc", &self.uniparc())?;
        state.serialize_field("uniprot_isoform", &self.uniprot_isoform())?;
        state.serialize_field("refseq", &self.refseq())?;
        state.serialize_field("transcript_version", &self.transcript_version())?;
        state.serialize_field("gene_version", &self.gene_version())?;
        state.serialize_field("gene_phenotype", &self.gene_phenotype())?;
        state.serialize_field("mirna", &self.mirna)?;
        optional_field(&mut state, "cds_length", &cds_length)?;
        state.serialize_field("cdna_length", &self.cdna_length())?;
        optional_field(&mut state, "source_cache", &source_cache)?;
        optional_field(&mut state, "bam_edit_status", &self.bam_edit_status)?;
        optional_field(&mut state, "refseq_match", &refseq_match)?;
        optional_field(&mut state, "refseq_offset", &self.refseq_offset)?;
        if self.plugin_data.is_empty() {
            state.skip_field("plugin_data")?;
        } else {
            state.serialize_field("plugin_data", &self.plugin_data)?;
        }
        state.serialize_field("feature_start", &self.feature_start)?;
        state.serialize_field("feature_end", &self.feature_end)?;
        optional_field(&mut state, "loftee_ctx", &self.loftee_ctx)?;
        optional_field(&mut state, "motif", &self.motif)?;
        optional_field(&mut state, "cell_type", &self.cell_type)?;
        state.end()
    }
}

/// Writes an `Option` field as a derive writes one marked
/// `skip_serializing_if = "Option::is_none"`: the value when `Some`, else the
/// field is reported to the serializer as skipped.
fn optional_field<S: serde::ser::SerializeStruct, T: serde::Serialize>(
    state: &mut S,
    key: &'static str,
    value: &Option<T>,
) -> Result<(), S::Error> {
    match value {
        Some(_) => state.serialize_field(key, value),
        None => state.skip_field(key),
    }
}

/// The `tsl` entry of a consequence document: the level as the number its
/// digit run spells, the reading Ensembl VEP's JSON format gives the text
/// (`numberify`, `Utils.pm:231`), exact within 64 bits and a double past
/// them. A level whose double is not finite is the text itself, the form the
/// JSON writer gives every value no JSON number holds; `null` without a level.
struct TslDocument<'a>(Option<&'a str>);

impl serde::Serialize for TslDocument<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let Some(digits) = self.0 else {
            return serializer.serialize_none();
        };
        if let Ok(level) = digits.parse::<u64>() {
            return serializer.serialize_u64(level);
        }
        match digits.parse::<f64>() {
            Ok(level) if level.is_finite() => serializer.serialize_f64(level),
            _ => serializer.serialize_str(digits),
        }
    }
}

/// The `flags` entry of a consequence document: the transcript's flags, then
/// `PICK` on a picked row.
struct FlagsDocument<'a>(&'a TranscriptConsequence);

impl serde::Serialize for FlagsDocument<'_> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeSeq;

        let flags = self.0.flags();
        let mut seq = serializer.serialize_seq(Some(flags.len() + usize::from(self.0.picked)))?;
        for flag in flags {
            seq.serialize_element(flag)?;
        }
        if self.0.picked {
            seq.serialize_element("PICK")?;
        }
        seq.end()
    }
}

/// The sequences of a transcript a plugin translates from: the translateable
/// sequence as the cache stores it (`_variation_effect_feature_cache`, the
/// start exon's phase padded as `N`), the 3' UTR, the protein and the codon
/// table. The protein is `None` when the cache lacks it; the UTR when neither
/// the cache nor the reference FASTA yields one ([`Self::from_transcript`]).
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct TranscriptSequences {
    pub cds: Box<str>,
    pub three_prime_utr: Option<Box<str>>,
    pub peptide: Option<Box<str>>,
    pub codon_table: u8,
}

impl TranscriptSequences {
    /// The sequences a plugin translates from, from each of the four fields.
    pub fn new(
        cds: Box<str>,
        three_prime_utr: Option<Box<str>>,
        peptide: Option<Box<str>>,
        codon_table: u8,
    ) -> Self {
        Self {
            cds,
            three_prime_utr,
            peptide,
            codon_table,
        }
    }

    /// The sequences of `transcript`; `None` for a transcript without a
    /// translateable sequence. `three_prime_utr` is kept as the caller resolved
    /// it, the cache's when it carries one, else the UTR read from the reference
    /// FASTA (`vep_effects::coding::compute_three_prime_utr_sequence`).
    pub fn from_transcript(
        transcript: &Transcript,
        three_prime_utr: Option<Vec<u8>>,
    ) -> Option<Self> {
        let vefc = transcript.vefc.as_ref()?;
        let cds = vefc.translateable_seq.as_deref()?;
        Some(Self {
            cds: Box::from(cds),
            three_prime_utr: three_prime_utr.map(|utr| Box::from(String::from_utf8_lossy(&utr))),
            peptide: vefc.peptide.as_deref().map(Box::from),
            codon_table: vefc.codon_table,
        })
    }
}

/// The fields a `MotifFeature` row of a sequence variant carries beside its
/// terms (`MotifFeatureVariationAllele_to_output_hash`, OutputFactory.pm:1919-1951).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct MotifFields {
    /// `MOTIF_NAME`: the binding matrix's stable id.
    pub name: Arc<str>,
    /// `MOTIF_POS`: the variant's first base as a position in the matrix,
    /// counted along the motif's strand; `None` when it falls outside the
    /// matrix.
    pub pos: Option<u32>,
    /// `HIGH_INF_POS`: whether the variant changes one base at a position of
    /// the matrix carrying more than
    /// [`INFORMATIVE_BITS`](crate::binding_matrix::INFORMATIVE_BITS) of
    /// information (`Y` or `N`); `None` for a matrix without counts, which
    /// omits the field.
    pub high_inf_pos: Option<bool>,
    /// `MOTIF_SCORE_CHANGE`: the motif's relative similarity score with the
    /// allele spliced in, minus its score as the cache stores it; `None` when
    /// the allele changes the motif's length, lies outside the matrix, holds a
    /// base other than A, C, G or T, or the matrix has no counts.
    pub score_change: Option<f64>,
    /// `TRANSCRIPTION_FACTORS`: the matrix's factor names, which print as an
    /// empty value when there are none.
    pub factors: Arc<[Arc<str>]>,
}

impl MotifFields {
    /// The fields of a `MotifFeature` row, from each of the five.
    pub fn new(
        name: Arc<str>,
        pos: Option<u32>,
        high_inf_pos: Option<bool>,
        score_change: Option<f64>,
        factors: Arc<[Arc<str>]>,
    ) -> Self {
        Self {
            name,
            pos,
            high_inf_pos,
            score_change,
            factors,
        }
    }
}

/// The activity of a regulatory feature in one epigenome: the five values of
/// the funcgen `regulatory_activity` table, printed as Ensembl VEP prints the
/// value of the feature's activity map (`ACTIVE`, `INACTIVE`, `POISED`,
/// `REPRESSED`, `NA`), and serialized as that spelling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
#[non_exhaustive]
pub enum Activity {
    Active,
    Inactive,
    Poised,
    Repressed,
    Na,
}

impl Activity {
    pub fn as_str(self) -> &'static str {
        match self {
            Activity::Active => "ACTIVE",
            Activity::Inactive => "INACTIVE",
            Activity::Poised => "POISED",
            Activity::Repressed => "REPRESSED",
            Activity::Na => "NA",
        }
    }
}

/// One entry of a row's `CELL_TYPE` list: an epigenome `--cell_type` named and
/// the regulatory feature's activity in it, printed `<name>:<activity>`
/// (OutputFactory.pm:1973-1975).
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct CellTypeActivity {
    pub name: Arc<str>,
    pub activity: Activity,
}

impl fmt::Display for CellTypeActivity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.name, self.activity.as_str())
    }
}

/// Type of genomic feature that a consequence is annotated against.
///
/// Most consequences are against `Transcript` features; `RegulatoryFeature`
/// and `MotifFeature` are used for regulatory consequence annotation. The
/// variants are declared in the order Ensembl VEP lists a variant's feature
/// overlaps (`VariationFeature.pm:855-864`: the transcript rows, then the
/// regulatory rows, then the motif rows), which `Ord` follows: the writers and
/// the picks order a record's rows by it.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, serde::Serialize, serde::Deserialize,
)]
pub enum FeatureType {
    Transcript,
    RegulatoryFeature,
    MotifFeature,
}

impl FeatureType {
    /// Static output string for this feature type.
    ///
    /// Byte-identical to the [`fmt::Display`] impl, but usable from an
    /// allocation-free writer path that emits bytes directly instead of routing
    /// through `core::fmt`.
    pub fn as_str(self) -> &'static str {
        match self {
            FeatureType::Transcript => "Transcript",
            FeatureType::RegulatoryFeature => "RegulatoryFeature",
            FeatureType::MotifFeature => "MotifFeature",
        }
    }
}

impl fmt::Display for FeatureType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// The kind of row a plugin runs on: a consequence against a feature of one of
/// the three [`FeatureType`] kinds, or the intergenic row of a variant that
/// overlaps no transcript.
///
/// The names are the values Ensembl VEP's `feature_types` takes
/// (`Bio::EnsEMBL::Variation::Utils::BaseVepPlugin`), which a shared library
/// lists in `vep_plugin_feature_types` and a built-in plugin returns from
/// `feature_types`; a plugin runs on a row only when its list names the row's
/// kind (`OutputFactory.pm:2232`, `check_feature_type(ref($bvfoa->feature) ||
/// 'Intergenic')`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum RowKind {
    Transcript,
    RegulatoryFeature,
    MotifFeature,
    Intergenic,
}

impl RowKind {
    /// Every kind, in the order Ensembl VEP's overlap list holds the rows.
    pub const ALL: &'static [RowKind] = &[
        RowKind::Transcript,
        RowKind::RegulatoryFeature,
        RowKind::MotifFeature,
        RowKind::Intergenic,
    ];

    /// The one kind a plugin runs on when it names none: Ensembl VEP's default
    /// `feature_types` (`BaseVepPlugin.pm:99`).
    pub const DEFAULT: [RowKind; 1] = [RowKind::Transcript];

    /// The name a plugin lists, `None` for any other string.
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "Transcript" => Some(Self::Transcript),
            "RegulatoryFeature" => Some(Self::RegulatoryFeature),
            "MotifFeature" => Some(Self::MotifFeature),
            "Intergenic" => Some(Self::Intergenic),
            _ => None,
        }
    }

    /// The name a plugin lists for this kind of row.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Transcript => "Transcript",
            Self::RegulatoryFeature => "RegulatoryFeature",
            Self::MotifFeature => "MotifFeature",
            Self::Intergenic => "Intergenic",
        }
    }

    /// The kind of a row: the feature's kind for a consequence row, the
    /// intergenic kind for a row against no feature.
    pub fn of_row(feature: Option<FeatureType>) -> Self {
        feature.map_or(Self::Intergenic, Self::from)
    }
}

impl From<FeatureType> for RowKind {
    fn from(feature: FeatureType) -> Self {
        match feature {
            FeatureType::Transcript => Self::Transcript,
            FeatureType::RegulatoryFeature => Self::RegulatoryFeature,
            FeatureType::MotifFeature => Self::MotifFeature,
        }
    }
}

/// `cDNA_position`, `CDS_position` and `Protein_position` of an allele's own span
/// ([`TranscriptConsequence::own_span_positions`]).
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct OwnSpanPositions {
    pub cdna: Option<String>,
    pub cds: Option<String>,
    pub protein: Option<String>,
}

impl Default for TranscriptConsequence {
    fn default() -> Self {
        Self {
            transcript: None,
            feature_id: None,
            feature_biotype: None,
            consequences: SmallVec::new(),
            impact: Impact::MODIFIER,
            canonical: false,
            cdna_position: None,
            cds_position: None,
            protein_position: None,
            cdna_total: None,
            cds_total: None,
            protein_total: None,
            own_span_positions: None,
            amino_acids: None,
            codons: None,
            protein_id_display: ProteinIdDisplay::AsCached,
            distance: None,
            strand: 1,
            exon: None,
            intron: None,
            hgvsc: None,
            hgvsp: None,
            hgvs_offset: None,
            hgvs_offset_of_hgvsc_alone: false,
            hgvs_intron_offsets: None,
            shift_length: 0,
            shifted_allele: None,
            used_ref: None,
            sift: None,
            polyphen: None,
            sift_term: None,
            polyphen_term: None,
            protein_start: None,
            domains: Vec::new(),
            feature_type: FeatureType::Transcript,
            picked: false,
            mirna: Vec::new(),
            bam_edit_status: None,
            refseq_offset: None,
            transcript_assembly_mismatch: false,
            plugin_data: indexmap::IndexMap::new(),
            plugin_empty_strings: Vec::new(),
            feature_start: 0,
            feature_end: 0,
            loftee_ctx: None,
            exon_ids: None,
            transcript_sequences: None,
            motif: None,
            cell_type: None,
            id_rank: 0,
            motif_db_id: None,
            record_span_in_cds: false,
        }
    }
}

impl TranscriptConsequence {
    /// Stores `value` under `key` as a plugin's value on this row, over any
    /// value an earlier plugin left there, as Ensembl VEP's `run_plugins`
    /// assigns each key of a plugin's result over the line's hash in load order
    /// (`OutputFactory.pm:2239-2241`). The key's mark in
    /// [`Self::plugin_empty_strings`] goes with the earlier value; `empty_string`
    /// names an empty `value` there as Ensembl VEP's empty string.
    pub fn set_plugin_value(&mut self, key: String, value: String, empty_string: bool) {
        self.plugin_empty_strings.retain(|marked| *marked != key);
        if empty_string && value.is_empty() {
            self.plugin_empty_strings.push(key.clone());
        }
        self.plugin_data.insert(key, value);
    }

    /// Removes the plugin value under `key` and its mark in
    /// [`Self::plugin_empty_strings`], keeping the order of the other values.
    pub fn remove_plugin_value(&mut self, key: &str) {
        if self.plugin_data.shift_remove(key).is_some() {
            self.plugin_empty_strings.retain(|marked| marked != key);
        }
    }

    /// A consequence against `transcript`, holding a shared handle on it, with
    /// every field the transcript determines filled from it (canonical, strand,
    /// span, edit status, rank) and every field the variant determines empty:
    /// no terms, `MODIFIER` impact, no positions, notations, predictions or
    /// plugin data. Built once per transcript a variant overlaps, so every field
    /// is spelled rather than taken from `Self::default()`: the handle's
    /// reference count is the row's whole cost, whatever a default value costs.
    pub fn for_transcript(transcript: &Arc<Transcript>) -> Self {
        Self {
            transcript: Some(Arc::clone(transcript)),
            feature_id: None,
            feature_biotype: None,
            consequences: SmallVec::new(),
            impact: Impact::MODIFIER,
            canonical: transcript.canonical,
            cdna_position: None,
            cds_position: None,
            protein_position: None,
            cdna_total: None,
            cds_total: None,
            protein_total: None,
            own_span_positions: None,
            amino_acids: None,
            codons: None,
            protein_id_display: ProteinIdDisplay::AsCached,
            distance: None,
            strand: transcript.strand.as_i8(),
            exon: None,
            intron: None,
            hgvsc: None,
            hgvsp: None,
            hgvs_offset: None,
            hgvs_offset_of_hgvsc_alone: false,
            hgvs_intron_offsets: None,
            shift_length: 0,
            shifted_allele: None,
            used_ref: None,
            sift: None,
            polyphen: None,
            sift_term: None,
            polyphen_term: None,
            protein_start: None,
            domains: Vec::new(),
            feature_type: FeatureType::Transcript,
            picked: false,
            mirna: Vec::new(),
            bam_edit_status: transcript.bam_edit_status,
            refseq_offset: None,
            transcript_assembly_mismatch: false,
            plugin_data: indexmap::IndexMap::new(),
            plugin_empty_strings: Vec::new(),
            feature_start: transcript.start,
            feature_end: transcript.end,
            loftee_ctx: None,
            exon_ids: None,
            transcript_sequences: None,
            motif: None,
            cell_type: None,
            id_rank: transcript.id_rank,
            motif_db_id: None,
            record_span_in_cds: false,
        }
    }

    /// The source transcript, `None` on a row against no transcript.
    pub fn transcript(&self) -> Option<&Transcript> {
        self.transcript.as_deref()
    }

    /// The key Ensembl VEP stores the row's feature variation under and sorts
    /// a record's rows of one kind by, as bytes, without the prefix the rejoin
    /// of a `--minimal` split puts before it (`OutputFactory.pm:2325`): a
    /// transcript's stable id, then `_` and its database id when the cache
    /// carries a non-zero one (`BaseVariationFeature.pm:310-320`, sorted at
    /// `VariationFeature.pm:666`), `_` alone when the cache's record holds no
    /// database id, as a cache built from a GTF or GFF3 does (Ensembl VEP
    /// numbers such a file's transcripts as it builds them, `BaseGXF.pm:489`),
    /// or the place of the gene-model file it was read from
    /// ([`file_order_key_suffix`]); a regulatory feature's stable id (`:918`,
    /// `:710`); a motif's database id in decimal (`:937`, `:750`), or its stable
    /// id when it has none. The byte order of these keys is the order
    /// [`Self::id_rank`] ranks a sequence variant's rows in; a structural
    /// variant's regulatory and motif rows carry the met rank instead
    /// ([`Self::structural_feature_order`]).
    pub fn ensembl_vep_key(&self) -> impl Iterator<Item = u8> + '_ {
        // A transcript of a gene-model file has no database id: its key ends
        // in the file's place instead ([`file_order_key_suffix`]).
        let file_suffix = match self.feature_type {
            FeatureType::Transcript => self
                .source_cache()
                .and_then(SourceCache::file_order)
                .map(file_order_key_suffix),
            _ => None,
        };
        let (stable_id, db_id, separator) = match self.feature_type {
            FeatureType::Transcript => {
                let db_id = match (&file_suffix, self.transcript()) {
                    (None, Some(transcript)) => Some(transcript.db_id),
                    _ => None,
                };
                (
                    self.transcript_id(),
                    db_id.flatten().filter(|&db_id| db_id != 0),
                    matches!(db_id, Some(db_id) if db_id != Some(0)),
                )
            }
            FeatureType::RegulatoryFeature => (self.transcript_id(), None, false),
            FeatureType::MotifFeature => match self.motif_db_id {
                Some(db_id) => ("", Some(db_id), false),
                None => (self.transcript_id(), None, false),
            },
        };
        let separator = separator.then_some(b'_');
        stable_id
            .bytes()
            .chain(separator)
            .chain(db_id.into_iter().flat_map(decimal_digits))
            .chain(file_suffix.into_iter().flat_map(String::into_bytes))
    }

    /// The order of two of one structural variant's `RegulatoryFeature` or
    /// `MotifFeature` rows of one kind, the order Ensembl VEP lists them in
    /// (`StructuralVariationFeature.pm:495-509` reads the lists
    /// `AnnotationType/RegFeat.pm:135` fills as the features are met): by
    /// [`Self::id_rank`], the rank of the order the source lists the features
    /// in, the record order of the files of `--regulatory_gff` or a cache baked
    /// from them and the listing order of a converted cache, then by the
    /// feature's start and its stable id for rows that share a rank.
    pub fn structural_feature_order(&self, other: &Self) -> std::cmp::Ordering {
        self.id_rank
            .cmp(&other.id_rank)
            .then_with(|| self.feature_start.cmp(&other.feature_start))
            .then_with(|| self.transcript_id().cmp(other.transcript_id()))
    }

    /// The feature's stable id: the transcript's (`ENST00000366667`), a
    /// regulatory or motif feature's (`ENSR00000000001`), or `""` on a row
    /// against no feature.
    pub fn transcript_id(&self) -> &str {
        self.transcript_id_handle().map_or("", |id| id)
    }

    /// The shared string behind [`Self::transcript_id`]: every row of one
    /// loaded feature returns the same allocation, so a reader may key on its
    /// address. `None` on a row against no feature.
    pub fn transcript_id_handle(&self) -> Option<&Arc<str>> {
        match &self.transcript {
            Some(transcript) => Some(&transcript.stable_id),
            None => self.feature_id.as_ref(),
        }
    }

    /// The gene's stable id (`ENSG00000142168`); `""` without a transcript.
    pub fn gene_id(&self) -> &str {
        self.gene_id_handle().map_or("", |id| id)
    }

    /// The shared string behind [`Self::gene_id`]; `None` without a transcript.
    pub fn gene_id_handle(&self) -> Option<&Arc<str>> {
        self.transcript().map(|t| &t.gene_stable_id)
    }

    /// The gene symbol (`SOD1`), `None` without one.
    pub fn gene_symbol(&self) -> Option<&str> {
        self.transcript()?.gene_symbol.as_deref()
    }

    /// The gene symbol's source (`HGNC`), `None` without one.
    pub fn gene_symbol_source(&self) -> Option<&str> {
        self.transcript()?.gene_symbol_source.as_deref()
    }

    /// The HGNC id, `None` without one.
    pub fn hgnc_id(&self) -> Option<&str> {
        self.transcript()?.hgnc_id.as_deref()
    }

    /// The feature's biotype: the transcript's (`protein_coding`), a regulatory
    /// feature's (`promoter`), or `None` on a row that prints none.
    pub fn biotype(&self) -> Option<&str> {
        match &self.transcript {
            Some(transcript) => Some(&transcript.biotype),
            None => self.feature_biotype.as_deref(),
        }
    }

    /// The protein id the row prints: the transcript's, with its translation
    /// version appended under `--protein_version`, or `None` on a row that
    /// prints none ([`ProteinIdDisplay`]).
    pub fn protein_id(&self) -> Option<Cow<'_, str>> {
        let id = self.transcript()?.protein_id.as_deref()?;
        match self.protein_id_display {
            ProteinIdDisplay::AsCached => Some(Cow::Borrowed(id)),
            ProteinIdDisplay::Hidden => None,
            ProteinIdDisplay::Versioned(version) => Some(Cow::Owned(format!("{id}.{version}"))),
        }
    }

    /// The transcript's flags (`cds_start_NF`, `cds_end_NF`), empty without a
    /// transcript; a picked row's `PICK` is [`Self::picked`].
    pub fn flags(&self) -> &[String] {
        self.transcript().map_or(&[], |t| &t.flags)
    }

    /// Transcript Support Level: the digit run of the transcript's `TSL`
    /// attribute as text ([`Transcript::tsl`]), `None` without one.
    pub fn tsl(&self) -> Option<&str> {
        self.transcript()?.tsl.as_deref()
    }

    /// The `MANE_Select` attribute's value, `Some` when the transcript carries
    /// the attribute, empty when it carries no value. `--pick`, `--pick_allele`
    /// and `--flag_pick` prefer MANE-curated transcripts over Ensembl-canonical
    /// ones (Ensembl VEP `pick_order` slot 0).
    pub fn mane_select(&self) -> Option<&str> {
        self.transcript()?.mane_select.as_deref()
    }

    /// The `MANE_Plus_Clinical` attribute's value, `Some` on presence as above
    /// (Ensembl VEP `pick_order` slot 1).
    pub fn mane_plus_clinical(&self) -> Option<&str> {
        self.transcript()?.mane_plus_clinical.as_deref()
    }

    /// The APPRIS annotation as the cache spells it (`principal1`,
    /// `alternative2`; Ensembl VEP `pick_order` slot 3).
    pub fn appris(&self) -> Option<&str> {
        self.transcript()?.appris.as_deref()
    }

    /// The CCDS id, `None` without one.
    pub fn ccds(&self) -> Option<&str> {
        self.transcript()?.ccds.as_deref()
    }

    /// The transcript carries a `gencode_primary` attribute: the
    /// `GENCODE_PRIMARY` column under `--flag_gencode_primary`.
    pub fn gencode_primary(&self) -> bool {
        self.transcript().is_some_and(|t| t.gencode_primary)
    }

    /// UniProtKB/Swiss-Prot accession(s), comma-separated as the cache
    /// carries them.
    pub fn swissprot(&self) -> Option<&str> {
        self.transcript()?.swissprot.as_deref()
    }

    /// UniProtKB/TrEMBL accession(s), comma-separated as the cache carries them.
    pub fn trembl(&self) -> Option<&str> {
        self.transcript()?.trembl.as_deref()
    }

    /// UniParc accession(s), comma-separated as the cache carries them.
    pub fn uniparc(&self) -> Option<&str> {
        self.transcript()?.uniparc.as_deref()
    }

    /// UniProtKB isoform accession(s), comma-separated as the cache carries them.
    pub fn uniprot_isoform(&self) -> Option<&str> {
        self.transcript()?.uniprot_isoform.as_deref()
    }

    /// RefSeq transcript id(s), comma-separated as the cache carries them.
    pub fn refseq(&self) -> Option<&str> {
        self.transcript()?.refseq.as_deref()
    }

    /// The transcript's stable-id version: the `.<version>` suffix of the
    /// `Feature` column under `--transcript_version`.
    pub fn transcript_version(&self) -> Option<u32> {
        self.transcript()?.version
    }

    /// The gene's stable-id version: the `.<version>` suffix of the `Gene`
    /// column under `--gene_version`.
    pub fn gene_version(&self) -> Option<u32> {
        self.transcript()?.gene_version
    }

    /// The gene has at least one phenotype association: the `GENE_PHENO`
    /// column under `--gene_phenotype`. Perl truthiness on the cached count
    /// (`$hash->{GENE_PHENO} = 1 if $self->{gene_phenotype} && $tr->{_gene_phenotype}`,
    /// OutputFactory.pm:1595 at release/116.3): a count of zero reads as no
    /// association.
    pub fn gene_phenotype(&self) -> bool {
        self.transcript()
            .and_then(|t| t.gene_phenotype)
            .is_some_and(|count| count > 0)
    }

    /// Length of the transcript's translateable sequence as the cache stores it,
    /// `None` for a transcript without a translation. The `length` pick
    /// criterion reads it (OutputFactory.pm:743-748 at release/116.3).
    pub fn cds_length(&self) -> Option<u32> {
        self.transcript()?.cds_length()
    }

    /// The transcript's spliced length, the sum of its exon lengths
    /// (`Bio::EnsEMBL::Transcript::length`), which the `length` pick criterion
    /// reads for a transcript without a translation; 0 without a transcript.
    pub fn cdna_length(&self) -> u32 {
        self.transcript().map_or(0, Transcript::spliced_length)
    }

    /// The cache the transcript came from in a merged cache (`_source_cache`),
    /// which the `ensembl` and `refseq` pick criteria read; `None` on every
    /// single-source cache, where both criteria score alike.
    pub fn source_cache(&self) -> Option<&SourceCache> {
        self.transcript()?.source_cache.as_ref()
    }

    /// The transcript's `rseq` match codes
    /// ([`RefSeqAlignment::match_codes`](crate::transcript::RefSeqAlignment::match_codes)):
    /// the `REFSEQ_MATCH` list under `--refseq` or `--merged`; `None` on a
    /// transcript without one.
    pub fn refseq_match(&self) -> Option<&[Box<str>]> {
        let codes = self.transcript()?.rseq_codes();
        (!codes.is_empty()).then_some(codes)
    }

    /// Base pairs of a structural variant's span that fall inside the feature
    /// and the percentage of the feature they cover, as VEP's `OverlapBP` and
    /// `OverlapPC` report them (OutputFactory.pm `StructuralVariationOverlap`
    /// branch); `None` when the spans do not meet or the feature span is unknown.
    pub fn feature_overlap(&self, variant_start: u64, variant_end: u64) -> Option<(u64, f64)> {
        if self.feature_end < self.feature_start || self.feature_end == 0 {
            return None;
        }
        let start = variant_start.max(self.feature_start);
        let end = variant_end.min(self.feature_end);
        if end < start {
            return None;
        }
        let bp = end - start + 1;
        let feature_len = self.feature_end - self.feature_start + 1;
        // The quotient first, then times 100, as OutputFactory.pm:2100 evaluates
        // `100 * ($overlap_length / $feature_length)`: the two orders round a
        // two-decimal tie differently (23 of 160 is 14.374999999999998 this way and
        // 14.375 the other, printed 14.37 and 14.38).
        Some((bp, 100.0 * (bp as f64 / feature_len as f64)))
    }

    /// The reference and alternate peptides of the row's amino-acid change, the
    /// two alleles Ensembl VEP's `peptide` reports for the reference and the
    /// variant `TranscriptVariationAllele`: `amino_acids` is `ref/alt` when the
    /// two differ and the one shared residue when they agree, so a synonymous
    /// row returns that residue twice. `None` when the row carries no
    /// amino-acid change.
    pub fn peptides(&self) -> Option<(&str, &str)> {
        let amino_acids = self.amino_acids.as_deref()?;
        Some(match amino_acids.split_once('/') {
            Some((reference, alternate)) => (reference, alternate),
            None => (amino_acids, amino_acids),
        })
    }
}

/// The decimal digits of `n`, most significant first: the bytes of the string
/// Ensembl VEP sorts a database id as ([`TranscriptConsequence::ensembl_vep_key`]).
fn decimal_digits(n: u64) -> impl Iterator<Item = u8> {
    let count = n.checked_ilog10().map_or(1, |magnitude| magnitude + 1);
    (0..count)
        .rev()
        .map(move |place| b'0' + ((n / 10u64.pow(place)) % 10) as u8)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinate::Strand;
    use crate::transcript::{Attribute, Exon, TranscriptVEFC, Translation};

    fn exon(start: u64, end: u64, rank: u32) -> Exon {
        Exon {
            start,
            end,
            rank,
            phase: -1,
            end_phase: -1,
        }
    }

    /// A one-exon forward transcript carrying the given phenotype count,
    /// stable-id versions and UniProt cross-references; every other field is
    /// empty.
    fn annotated_transcript(gene_phenotype: Option<u32>) -> Transcript {
        Transcript {
            stable_id: "ENST00000000001".into(),
            version: Some(7),
            db_id: None,
            gene_stable_id: "ENSG00000000001".into(),
            gene_version: Some(12),
            chr: "1".into(),
            start: 100,
            end: 200,
            strand: Strand::Forward,
            biotype: "protein_coding".into(),
            source: "Ensembl".into(),
            gene_symbol: None,
            gene_symbol_source: None,
            gene_symbol_id: None,
            hgnc_id: None,
            gene_phenotype,
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
            uniparc: Some("UPI0000000001".to_string()),
            uniprot_isoform: Some("P00001-2,P00001-3".to_string()),
            exons: vec![exon(100, 200, 1)],
            exon_ids: None,
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
            gencode_basic: false,
            gencode_primary: false,
            attributes: vec![],
            vefc: None,
            source_cache: None,
            bam_edit_status: None,
            circular_length: None,
            region_files: None,
            derived: Default::default(),
            id_rank: 0,
        }
    }

    /// A consequence reads the transcript's versions, UniProt cross-references
    /// and phenotype flag through its handle, as the cache carries them, the
    /// miRNA terms wait for the variant, and `GENE_PHENO` follows Perl truthiness
    /// on the cached phenotype count
    /// (`$hash->{GENE_PHENO} = 1 if $self->{gene_phenotype} && $tr->{_gene_phenotype}`,
    /// OutputFactory.pm:1595 at release/116.3): a count of zero reads as no
    /// association.
    #[test]
    fn accessors_read_the_borrowed_transcript() {
        let transcript = Arc::new(annotated_transcript(Some(3)));
        let tc = TranscriptConsequence::for_transcript(&transcript);
        assert!(Arc::ptr_eq(tc.transcript.as_ref().unwrap(), &transcript));
        assert_eq!(tc.transcript_id(), "ENST00000000001");
        assert_eq!(tc.gene_id(), "ENSG00000000001");
        assert_eq!(tc.biotype(), Some("protein_coding"));
        assert_eq!(tc.transcript_version(), Some(7));
        assert_eq!(tc.gene_version(), Some(12));
        assert_eq!(tc.uniparc(), Some("UPI0000000001"));
        assert_eq!(tc.uniprot_isoform(), Some("P00001-2,P00001-3"));
        assert_eq!(tc.cdna_length(), 101);
        assert_eq!(tc.cds_length(), None);
        assert!(tc.gene_phenotype());
        assert!(tc.mirna.is_empty());
        assert!(tc.flags().is_empty());
        assert_eq!(tc.protein_id(), None);

        let zero = TranscriptConsequence::for_transcript(&Arc::new(annotated_transcript(Some(0))));
        assert!(!zero.gene_phenotype());
        let unknown = TranscriptConsequence::for_transcript(&Arc::new(annotated_transcript(None)));
        assert!(!unknown.gene_phenotype());

        let none = TranscriptConsequence::default();
        assert_eq!(none.transcript_id(), "");
        assert_eq!(none.gene_id(), "");
        assert_eq!(none.biotype(), None);
        assert_eq!(none.cdna_length(), 0);
        let regulatory = TranscriptConsequence {
            feature_id: Some(Arc::from("ENSR00000000001")),
            feature_biotype: Some(Arc::from("promoter")),
            feature_type: FeatureType::RegulatoryFeature,
            ..TranscriptConsequence::default()
        };
        assert_eq!(regulatory.transcript_id(), "ENSR00000000001");
        assert_eq!(regulatory.biotype(), Some("promoter"));
        assert_eq!(regulatory.gene_id(), "");
    }

    /// The protein id prints as the cache carries it, with the translation's
    /// version appended under `--protein_version`, or not at all on a row that
    /// hides it.
    #[test]
    fn protein_id_follows_its_display() {
        let mut transcript = annotated_transcript(None);
        transcript.protein_id = Some("ENSP00000000001".into());
        let mut tc = TranscriptConsequence::for_transcript(&Arc::new(transcript));
        assert_eq!(tc.protein_id().as_deref(), Some("ENSP00000000001"));
        tc.protein_id_display = ProteinIdDisplay::Versioned(3);
        assert_eq!(tc.protein_id().as_deref(), Some("ENSP00000000001.3"));
        tc.protein_id_display = ProteinIdDisplay::Hidden;
        assert_eq!(tc.protein_id(), None);
    }

    /// A transcript carrying every attribute the document prints: the versions,
    /// symbol, cross-references, pick criteria, flags, a translation of 300
    /// bases, one exon of 1,500 bases and an `rseq` match code.
    fn fully_annotated_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST00000000001".into(),
            version: Some(7),
            gene_stable_id: "ENSG00000000001".into(),
            gene_version: Some(12),
            chr: "1".into(),
            start: 100,
            end: 1599,
            strand: Strand::Reverse,
            biotype: "protein_coding".into(),
            source: "Ensembl".into(),
            gene_symbol: Some("GENE1".into()),
            gene_symbol_source: Some("HGNC".into()),
            hgnc_id: Some("HGNC:1".into()),
            gene_phenotype: Some(3),
            canonical: true,
            mane_select: Some("NM_000001.1".to_string()),
            mane_plus_clinical: Some("NM_000002.1".to_string()),
            tsl: Some("1".to_string()),
            appris: Some("principal1".to_string()),
            ccds: Some("CCDS1.1".to_string()),
            protein_id: Some("ENSP00000000001".into()),
            refseq: Some("NM_000001.1".to_string()),
            swissprot: Some("P00001".to_string()),
            trembl: Some("Q00001".to_string()),
            uniparc: Some("UPI0000000001".to_string()),
            uniprot_isoform: Some("P00001-2".to_string()),
            exons: vec![Exon {
                start: 100,
                end: 1599,
                rank: 1,
                phase: 0,
                end_phase: 0,
            }],
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(300),
            coding_region_start: Some(100),
            coding_region_end: Some(399),
            translation_start: Some(100),
            translation_end: Some(399),
            translation: Some(Translation {
                stable_id: "ENSP00000000001".into(),
                version: Some(2),
                db_id: None,
                start: 1,
                end: 300,
                seq: None,
            }),
            flags: Arc::from(["cds_start_NF".to_string()]),
            gencode_basic: true,
            gencode_primary: true,
            attributes: vec![Attribute {
                code: "rseq_mrna_match".into(),
                value: String::new(),
            }],
            vefc: Some(TranscriptVEFC {
                translateable_seq: Some("A".repeat(300)),
                ..Default::default()
            }),
            source_cache: Some(SourceCache::Ensembl),
            bam_edit_status: Some(BamEditStatus::Ok),
            id_rank: 5,
            ..Transcript::default()
        }
    }

    /// A consequence with every optional key populated, against
    /// [`fully_annotated_transcript`].
    fn fully_populated_consequence() -> TranscriptConsequence {
        let transcript = Arc::new(fully_annotated_transcript());
        let mut tc = TranscriptConsequence::for_transcript(&transcript);
        tc.consequences = smallvec::smallvec![
            Consequence::MissenseVariant,
            Consequence::SpliceRegionVariant
        ];
        tc.impact = Impact::MODERATE;
        tc.cdna_position = Some("428".into());
        tc.cds_position = Some("301".into());
        tc.protein_position = Some("101".into());
        tc.cdna_total = Some(1500);
        tc.cds_total = Some(300);
        tc.protein_total = Some(99);
        tc.own_span_positions = Some(Box::new(OwnSpanPositions {
            cdna: Some("439".into()),
            cds: Some("87".into()),
            protein: Some("29".into()),
        }));
        tc.amino_acids = Some("A/V".into());
        tc.codons = Some("gCc/gTc".into());
        tc.distance = Some(5);
        tc.exon = Some("3/5".into());
        tc.intron = Some("2/4".into());
        tc.hgvsc = Some("ENST00000000001.7:c.301G>A".into());
        tc.hgvsp = Some("ENSP00000000001.2:p.Ala101Val".into());
        tc.hgvs_offset = Some(2);
        tc.shift_length = 1;
        tc.shifted_allele = Some(Arc::from("GAA"));
        tc.used_ref = Some(Arc::from("A"));
        tc.sift = Some("deleterious(0.01)".into());
        tc.polyphen = Some("probably_damaging(0.998)".into());
        tc.domains = vec![("Pfam".to_string(), "PF00080".to_string())];
        tc.mirna = vec!["miRNA_stem".to_string()];
        tc.refseq_offset = Some(3);
        tc.transcript_assembly_mismatch = true;
        tc.plugin_data
            .insert("REVEL_score".to_string(), "0.9".to_string());
        tc.loftee_ctx = Some(Box::new(LofteeContext::from_transcript(&transcript)));
        tc.motif = Some(Box::new(MotifFields {
            name: Arc::from("ENSPFM0290"),
            pos: Some(7),
            high_inf_pos: Some(true),
            score_change: Some(-0.5),
            factors: Arc::from([Arc::from("HNF4A")]),
        }));
        tc.cell_type = Some(Box::from([CellTypeActivity {
            name: Arc::from("HeLa-S3"),
            activity: Activity::Active,
        }]));
        tc.record_span_in_cds = true;
        tc
    }

    /// The top-level keys of a JSON object in document order.
    struct KeyOrder(Vec<String>);

    impl<'de> serde::Deserialize<'de> for KeyOrder {
        fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
            struct Keys;
            impl<'de> serde::de::Visitor<'de> for Keys {
                type Value = KeyOrder;
                fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                    f.write_str("a JSON object")
                }
                fn visit_map<A: serde::de::MapAccess<'de>>(
                    self,
                    mut map: A,
                ) -> Result<KeyOrder, A::Error> {
                    let mut keys = Vec::new();
                    while let Some(key) = map.next_key::<String>()? {
                        map.next_value::<serde::de::IgnoredAny>()?;
                        keys.push(key);
                    }
                    Ok(KeyOrder(keys))
                }
            }
            deserializer.deserialize_map(Keys)
        }
    }

    fn key_order(tc: &TranscriptConsequence) -> Vec<String> {
        let text = serde_json::to_string(tc).unwrap();
        serde_json::from_str::<KeyOrder>(&text).unwrap().0
    }

    /// The document of a consequence with every optional key populated carries
    /// the keys of the fixture, in its order: the key list of a serde derive
    /// over a row that holds the transcript's values as fields of its own. An
    /// empty consequence carries the unconditional keys alone, in the same
    /// order.
    #[test]
    fn serialize_emits_the_derive_key_set_in_order() {
        let expected: Vec<&str> = include_str!("../tests/fixtures/transcript_consequence_keys.txt")
            .lines()
            .collect();
        assert_eq!(key_order(&fully_populated_consequence()), expected);
        let listed: Vec<&str> = TRANSCRIPT_CONSEQUENCE_KEYS
            .iter()
            .map(|(key, _)| *key)
            .collect();
        assert_eq!(listed, expected);
        let unconditional: Vec<&str> = TRANSCRIPT_CONSEQUENCE_KEYS
            .iter()
            .filter(|(_, optional)| !optional)
            .map(|(key, _)| *key)
            .collect();
        assert_eq!(key_order(&TranscriptConsequence::default()), unconditional);
    }

    /// The whole document text of the fully populated consequence and of the
    /// empty one equals the fixture documents of that derive: every value read
    /// through an accessor is encoded as the field of the same name is.
    #[test]
    fn serialize_writes_the_derive_document() {
        assert_eq!(
            serde_json::to_string(&fully_populated_consequence()).unwrap(),
            include_str!("../tests/fixtures/transcript_consequence_document.json").trim_end()
        );
        assert_eq!(
            serde_json::to_string(&TranscriptConsequence::default()).unwrap(),
            include_str!("../tests/fixtures/transcript_consequence_empty.json").trim_end()
        );
        let mut picked =
            TranscriptConsequence::for_transcript(&Arc::new(fully_annotated_transcript()));
        picked.picked = true;
        let value = serde_json::to_value(&picked).unwrap();
        assert_eq!(value["flags"], serde_json::json!(["cds_start_NF", "PICK"]));
    }

    /// The document's `tsl` is the number the level's digit run spells: `007`
    /// is 7, a run past 64 bits a double, one past the double's range the text
    /// itself, and `null` without a level.
    #[test]
    fn serialize_writes_the_level_as_its_number() {
        let tsl = |level: Option<&str>| {
            let mut transcript = fully_annotated_transcript();
            transcript.tsl = level.map(str::to_string);
            let row = TranscriptConsequence::for_transcript(&Arc::new(transcript));
            serde_json::to_value(&row).unwrap()["tsl"].clone()
        };
        assert_eq!(tsl(Some("007")), serde_json::json!(7));
        assert_eq!(tsl(Some("00")), serde_json::json!(0));
        assert_eq!(tsl(Some("999")), serde_json::json!(999));
        assert_eq!(
            tsl(Some("18446744073709551616")),
            serde_json::json!(18446744073709551616.0)
        );
        let past_the_double = "9".repeat(400);
        assert_eq!(
            tsl(Some(&past_the_double)),
            serde_json::Value::String(past_the_double)
        );
        assert_eq!(tsl(None), serde_json::Value::Null);
    }

    /// The coding sequence, protein and codon table come from
    /// `_variation_effect_feature_cache` as the cache stores them, the 3' UTR
    /// from the caller; a transcript whose cache lacks a translateable sequence,
    /// or that has no cache, has no sequences.
    #[test]
    fn transcript_sequences_come_from_the_cache_and_the_callers_utr() {
        let mut transcript = annotated_transcript(None);
        transcript.vefc = Some(TranscriptVEFC {
            codon_table: 2,
            translateable_seq: Some("NNATGCCCTAA".to_string()),
            peptide: Some("XP".to_string()),
            ..Default::default()
        });
        let with_utr =
            TranscriptSequences::from_transcript(&transcript, Some(b"GGGTAG".to_vec())).unwrap();
        assert_eq!(
            with_utr,
            TranscriptSequences {
                cds: Box::from("NNATGCCCTAA"),
                three_prime_utr: Some(Box::from("GGGTAG")),
                peptide: Some(Box::from("XP")),
                codon_table: 2,
            }
        );
        let without_utr = TranscriptSequences::from_transcript(&transcript, None).unwrap();
        assert_eq!(
            without_utr,
            TranscriptSequences {
                three_prime_utr: None,
                ..with_utr
            }
        );

        transcript.vefc.as_mut().unwrap().peptide = None;
        assert_eq!(
            TranscriptSequences::from_transcript(&transcript, None)
                .unwrap()
                .peptide,
            None
        );
        transcript.vefc.as_mut().unwrap().translateable_seq = None;
        assert_eq!(
            TranscriptSequences::from_transcript(&transcript, Some(b"GGG".to_vec())),
            None
        );
        transcript.vefc = None;
        assert_eq!(
            TranscriptSequences::from_transcript(&transcript, None),
            None
        );
    }

    /// `OVERLAP_CONSEQUENCES` (ensembl-variation `Utils/Constants.pm`) marks
    /// thirteen terms `include => {coding => 1}`; `coding_transcript_variant`
    /// is not among them.
    #[test]
    fn coding_include_terms_are_the_thirteen_of_overlap_consequences() {
        let coding: Vec<&str> = Consequence::ALL
            .iter()
            .filter(|c| c.is_coding_include())
            .map(|c| c.so_term())
            .collect();
        assert_eq!(
            coding,
            [
                "stop_gained",
                "frameshift_variant",
                "stop_lost",
                "start_lost",
                "inframe_insertion",
                "inframe_deletion",
                "missense_variant",
                "protein_altering_variant",
                "incomplete_terminal_codon_variant",
                "start_retained_variant",
                "stop_retained_variant",
                "synonymous_variant",
                "coding_sequence_variant",
            ]
        );
        assert!(!Consequence::CodingTranscriptVariant.is_coding_include());
    }

    /// `ALL` lists the terms in rank order, so a term's position is its
    /// discriminant and one less than its rank.
    #[test]
    fn all_lists_every_term_in_rank_order() {
        for (index, term) in Consequence::ALL.iter().enumerate() {
            assert_eq!(*term as usize, index);
            assert_eq!(term.rank() as usize, index + 1);
        }
    }

    #[test]
    fn derive_introns_three_exon_forward() {
        // Exons at 100-200, 300-400, 500-600 -> introns 201-299, 401-499.
        let exons = vec![exon(100, 200, 1), exon(300, 400, 2), exon(500, 600, 3)];
        let introns = derive_introns_from_exons(&exons);
        assert_eq!(introns.len(), 2);
        assert_eq!(
            (introns[0].start, introns[0].end, introns[0].rank),
            (201, 299, 1)
        );
        assert_eq!(
            (introns[1].start, introns[1].end, introns[1].rank),
            (401, 499, 2)
        );
    }

    #[test]
    fn derive_introns_single_exon_is_empty() {
        let exons = vec![exon(100, 200, 1)];
        assert!(derive_introns_from_exons(&exons).is_empty());
    }

    #[test]
    fn derive_introns_adjacent_exons_no_gap() {
        // Exons touching with no intronic gap (200-201) produce no intron.
        let exons = vec![exon(100, 200, 1), exon(201, 300, 2)];
        assert!(derive_introns_from_exons(&exons).is_empty());
    }

    #[test]
    fn derive_introns_orders_by_position_not_input_order() {
        // Reverse-strand-style: rank order opposite to genomic order. The helper
        // sorts by genomic start, so introns come out in coordinate order.
        let exons = vec![exon(500, 600, 1), exon(300, 400, 2), exon(100, 200, 3)];
        let introns = derive_introns_from_exons(&exons);
        assert_eq!(introns.len(), 2);
        assert_eq!((introns[0].start, introns[0].end), (201, 299));
        assert_eq!((introns[1].start, introns[1].end), (401, 499));
    }

    #[test]
    fn test_consequence_ranks_are_ordered() {
        assert!(Consequence::TranscriptAblation.rank() < Consequence::MissenseVariant.rank());
        assert!(Consequence::MissenseVariant.rank() < Consequence::SynonymousVariant.rank());
        assert!(Consequence::SynonymousVariant.rank() < Consequence::IntronVariant.rank());
        assert!(Consequence::IntronVariant.rank() < Consequence::IntergenicVariant.rank());
    }

    #[test]
    fn test_consequence_impacts() {
        assert_eq!(Consequence::StopGained.impact(), Impact::HIGH);
        assert_eq!(Consequence::MissenseVariant.impact(), Impact::MODERATE);
        assert_eq!(Consequence::SynonymousVariant.impact(), Impact::LOW);
        assert_eq!(Consequence::IntronVariant.impact(), Impact::MODIFIER);
    }

    #[test]
    fn test_consequence_so_terms() {
        assert_eq!(Consequence::StopGained.so_term(), "stop_gained");
        assert_eq!(
            Consequence::FivePrimeUtrVariant.so_term(),
            "5_prime_UTR_variant"
        );
        assert_eq!(
            Consequence::NmdTranscriptVariant.so_term(),
            "NMD_transcript_variant"
        );
        assert_eq!(Consequence::TfbsAblation.so_term(), "TFBS_ablation");
    }

    #[test]
    fn test_consequence_roundtrip() {
        for csq in [
            Consequence::TranscriptAblation,
            Consequence::MissenseVariant,
            Consequence::IntergenicVariant,
            Consequence::MatureMirnaVariant,
            Consequence::FivePrimeUtrVariant,
            Consequence::NmdTranscriptVariant,
        ] {
            assert_eq!(Consequence::from_so_term(csq.so_term()), Some(csq));
        }
    }

    #[test]
    fn test_most_severe() {
        let csqs = vec![
            Consequence::IntronVariant,
            Consequence::MissenseVariant,
            Consequence::SynonymousVariant,
        ];
        assert_eq!(
            Consequence::most_severe(&csqs),
            Some(Consequence::MissenseVariant)
        );
    }

    #[test]
    fn test_consequence_ordering() {
        let mut csqs = [
            Consequence::IntergenicVariant,
            Consequence::StopGained,
            Consequence::MissenseVariant,
        ];
        csqs.sort();
        assert_eq!(csqs[0], Consequence::StopGained);
        assert_eq!(csqs[1], Consequence::MissenseVariant);
        assert_eq!(csqs[2], Consequence::IntergenicVariant);
    }

    /// The feature kinds order as Ensembl VEP lists a variant's overlaps
    /// (`VariationFeature.pm:855-864`): transcript, regulatory, motif.
    #[test]
    fn feature_types_order_as_the_overlap_list() {
        let mut kinds = [
            FeatureType::MotifFeature,
            FeatureType::Transcript,
            FeatureType::RegulatoryFeature,
        ];
        kinds.sort();
        assert_eq!(
            kinds,
            [
                FeatureType::Transcript,
                FeatureType::RegulatoryFeature,
                FeatureType::MotifFeature
            ]
        );
    }

    /// A transcript read from a gene-model file keys after the cache's
    /// transcript of the same stable id, whatever that one's database id, and
    /// the files' transcripts of one id follow in argument order: the key ends
    /// in `_:` and the file's place, where Ensembl VEP's ends in a count of the
    /// transcripts it built before (`BaseGXF.pm:489`), which the input's other
    /// records and Perl's hash order decide.
    #[test]
    fn a_gene_model_file_transcript_keys_after_the_cache_transcript_of_its_id() {
        let key = |row: &TranscriptConsequence| {
            String::from_utf8(row.ensembl_vep_key().collect()).unwrap()
        };
        let row = |db_id: Option<u64>, source_cache: Option<SourceCache>| {
            TranscriptConsequence::for_transcript(&Arc::new(Transcript {
                stable_id: "ENST00000000001".into(),
                db_id,
                source_cache,
                ..annotated_transcript(None)
            }))
        };
        let file = |order: u32| {
            row(
                None,
                Some(SourceCache::GeneModelFile {
                    short_name: "genes.gtf.gz".into(),
                    order,
                }),
            )
        };
        assert_eq!(key(&file(1)), "ENST00000000001_:0000000001");
        for cache in [
            row(Some(7), None),
            row(Some(99_999_999_999), None),
            row(None, None),
            row(Some(5), Some(SourceCache::RefSeq)),
        ] {
            assert!(cache.ensembl_vep_key().lt(file(1).ensembl_vep_key()));
        }
        assert!(file(1).ensembl_vep_key().lt(file(2).ensembl_vep_key()));
        assert!(file(9).ensembl_vep_key().lt(file(10).ensembl_vep_key()));
        // Among other stable ids the key stands where a numbered key stands:
        // before an id that continues `_` with a letter.
        let longer = TranscriptConsequence::for_transcript(&Arc::new(Transcript {
            stable_id: "ENST00000000001_PAR_Y".into(),
            db_id: Some(3),
            ..annotated_transcript(None)
        }));
        assert!(file(1).ensembl_vep_key().lt(longer.ensembl_vep_key()));
        assert!(row(Some(7), None)
            .ensembl_vep_key()
            .lt(longer.ensembl_vep_key()));
    }

    /// A cache transcript whose record holds no database id keys as its
    /// stable id and `_`, where Ensembl VEP's key of a GTF or GFF3 transcript
    /// holds `_` and a counter for every transcript but the first it builds
    /// (`BaseGXF.pm:489`, `BaseVariationFeature.pm:316-318`): an id that
    /// continues the stable id with `-` keys first, as it does among numbered
    /// keys. A database id of zero leaves the key bare.
    #[test]
    fn a_transcript_without_a_database_id_keys_as_its_id_and_an_underscore() {
        let key = |row: &TranscriptConsequence| {
            String::from_utf8(row.ensembl_vep_key().collect()).unwrap()
        };
        let row = |stable_id: &str, db_id: Option<u64>| {
            TranscriptConsequence::for_transcript(&Arc::new(Transcript {
                stable_id: stable_id.into(),
                db_id,
                ..annotated_transcript(None)
            }))
        };
        assert_eq!(key(&row("MODEL1", None)), "MODEL1_");
        assert_eq!(key(&row("MODEL1", Some(12))), "MODEL1_12");
        assert_eq!(key(&row("MODEL1", Some(0))), "MODEL1");
        assert!(row("MODEL1-A", None)
            .ensembl_vep_key()
            .lt(row("MODEL1", None).ensembl_vep_key()));
        assert!(row("MODEL1-A", Some(4))
            .ensembl_vep_key()
            .lt(row("MODEL1", Some(3)).ensembl_vep_key()));
        assert!(row("MODEL1", Some(0))
            .ensembl_vep_key()
            .lt(row("MODEL1-A", Some(4)).ensembl_vep_key()));
    }

    /// The key Ensembl VEP 116.3 stores a row's feature variation under: a
    /// transcript's stable id with `_` and its database id, which `sort keys`
    /// orders before the stable id alone when that id is a prefix of the other
    /// (`AT4G12020.10_118783` before `AT4G12020.1_118785`, `_` sorting after
    /// every digit; BaseVariationFeature.pm:310-320, VariationFeature.pm:666); a
    /// transcript with a database id of zero keeps its stable id alone, as
    /// `if(my $dbID = $tr->dbID)` skips it, and one whose cache record holds
    /// none ends in `_` alone; a regulatory feature's stable
    /// id (VariationFeature.pm:918); a motif's database id as the digits Ensembl VEP
    /// sorts as a string (MotifFeatureVariation.pm:70-76, VariationFeature.pm:937),
    /// and its stable id when it has none.
    #[test]
    fn ensembl_vep_key_is_the_stable_id_with_the_database_id() {
        let key = |row: &TranscriptConsequence| {
            String::from_utf8(row.ensembl_vep_key().collect()).unwrap()
        };
        let transcript = |stable_id: &str, db_id: Option<u64>| {
            TranscriptConsequence::for_transcript(&Arc::new(Transcript {
                stable_id: stable_id.into(),
                db_id,
                ..annotated_transcript(None)
            }))
        };
        let ten = transcript("AT4G12020.10", Some(118783));
        let one = transcript("AT4G12020.1", Some(118785));
        assert_eq!(key(&ten), "AT4G12020.10_118783");
        assert_eq!(key(&one), "AT4G12020.1_118785");
        assert!(ten.ensembl_vep_key().lt(one.ensembl_vep_key()));
        assert!(one.transcript_id() < ten.transcript_id());
        assert_eq!(
            key(&transcript("ENST00000000001", None)),
            "ENST00000000001_"
        );
        assert_eq!(
            key(&transcript("ENST00000000001", Some(0))),
            "ENST00000000001"
        );

        let feature = |stable_id: &str, feature_type: FeatureType, motif_db_id: Option<u64>| {
            TranscriptConsequence {
                feature_id: Some(stable_id.into()),
                feature_type,
                motif_db_id,
                ..TranscriptConsequence::default()
            }
        };
        assert_eq!(
            key(&feature(
                "ENSR21_BBBWM",
                FeatureType::RegulatoryFeature,
                Some(9)
            )),
            "ENSR21_BBBWM"
        );
        let motif = feature("ENSM00000575394", FeatureType::MotifFeature, Some(575394));
        assert_eq!(key(&motif), "575394");
        assert_eq!(
            key(&feature("ENSM00000575394", FeatureType::MotifFeature, None)),
            "ENSM00000575394"
        );
        assert_eq!(
            key(&feature("ENSM_TEN", FeatureType::MotifFeature, Some(10))),
            "10"
        );
        assert_eq!(
            key(&feature(
                "ENSM_MAX",
                FeatureType::MotifFeature,
                Some(u64::MAX)
            )),
            u64::MAX.to_string()
        );
    }

    /// Database ids on both sides of six boundaries of the digit count, from
    /// one digit to two up to nineteen to twenty, with a prefix pair, zero and
    /// the largest value, sort by their `decimal_digits` as their decimal
    /// strings do, the order of Perl's string sort. The regulatory index of
    /// `vep-cli` holds its `decimal_string_order` to the same list and the same
    /// expected order.
    #[test]
    fn database_ids_sort_as_their_decimal_strings() {
        const IDS: [u64; 18] = [
            0,
            1,
            9,
            10,
            11,
            57,
            99,
            100,
            570,
            999,
            1_000,
            999_999,
            1_000_000,
            9_999_999_999,
            10_000_000_000,
            9_999_999_999_999_999_999,
            10_000_000_000_000_000_000,
            u64::MAX,
        ];
        const BY_DECIMAL_STRING: [u64; 18] = [
            0,
            1,
            10,
            100,
            1_000,
            1_000_000,
            10_000_000_000,
            10_000_000_000_000_000_000,
            11,
            u64::MAX,
            57,
            570,
            9,
            99,
            999,
            999_999,
            9_999_999_999,
            9_999_999_999_999_999_999,
        ];
        let mut as_strings = IDS;
        as_strings.sort_by_cached_key(|id| id.to_string());
        assert_eq!(as_strings, BY_DECIMAL_STRING);
        let mut sorted = IDS;
        sorted.sort_by(|a, b| decimal_digits(*a).cmp(decimal_digits(*b)));
        assert_eq!(sorted, BY_DECIMAL_STRING);
    }

    /// `amino_acids` prints `ref/alt` for a change and the one residue for a
    /// synonymous row, so the peptides are the two sides or that residue twice.
    #[test]
    fn peptides_split_the_amino_acid_change() {
        let row = |amino_acids: Option<&str>| TranscriptConsequence {
            amino_acids: amino_acids.map(str::to_string),
            ..TranscriptConsequence::default()
        };
        assert_eq!(row(Some("A/V")).peptides(), Some(("A", "V")));
        assert_eq!(row(Some("K")).peptides(), Some(("K", "K")));
        assert_eq!(row(Some("-/H")).peptides(), Some(("-", "H")));
        assert_eq!(row(Some("A/VX")).peptides(), Some(("A", "VX")));
        assert_eq!(row(None).peptides(), None);
    }

    /// A plugin's value written on a row replaces the earlier value under its
    /// key with that value's mark; an empty value is marked only where the
    /// writer names it Ensembl VEP's empty string, and a removal takes the mark
    /// and keeps the order of the other values.
    #[test]
    fn a_plugin_value_written_on_a_row_replaces_the_earlier_value_and_its_mark() {
        let mut row = TranscriptConsequence::default();
        row.set_plugin_value("A".to_string(), String::new(), true);
        row.set_plugin_value("B".to_string(), "b".to_string(), true);
        row.set_plugin_value("C".to_string(), String::new(), false);
        assert_eq!(row.plugin_empty_strings, ["A"]);

        row.set_plugin_value("A".to_string(), String::new(), false);
        assert_eq!(row.plugin_data["A"], "");
        assert!(row.plugin_empty_strings.is_empty());
        row.set_plugin_value("A".to_string(), String::new(), true);
        row.set_plugin_value("A".to_string(), "a".to_string(), true);
        assert_eq!(row.plugin_data["A"], "a");
        assert!(row.plugin_empty_strings.is_empty());

        row.set_plugin_value("C".to_string(), String::new(), true);
        row.set_plugin_value("D".to_string(), "d".to_string(), false);
        row.remove_plugin_value("B");
        row.remove_plugin_value("Absent");
        assert_eq!(row.plugin_data.keys().collect::<Vec<_>>(), ["A", "C", "D"]);
        assert_eq!(row.plugin_empty_strings, ["C"]);
        row.remove_plugin_value("C");
        assert_eq!(row.plugin_data.keys().collect::<Vec<_>>(), ["A", "D"]);
        assert!(row.plugin_empty_strings.is_empty());
    }

    /// The plugin empty-string mark is a three-word list on the row and on the
    /// record, so a row measures its fields' 704 bytes plus the list and a
    /// record its 1288 plus the list on a 64-bit target.
    #[test]
    #[cfg(target_pointer_width = "64")]
    fn plugin_empty_mark_costs_a_row_and_a_record_a_list() {
        assert_eq!(size_of::<Vec<String>>(), 3 * size_of::<usize>());
        assert_eq!(size_of::<TranscriptConsequence>(), 728);
        assert_eq!(size_of::<crate::variant::InputVariant>(), 1312);
    }
}
