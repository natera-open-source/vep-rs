// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Consequence types and rankings from Ensembl Variation.
//!
//! All 41 Sequence Ontology (SO) consequence terms used by VEP, with their
//! exact ranks, impact classifications, and SO accessions matching the Perl
//! implementation in `Bio::EnsEMBL::Variation::Utils::Constants` (ensembl-variation
//! release/115). `OutputFactory.pm` is `Bio/EnsEMBL/VEP/OutputFactory.pm` in
//! ensembl-vep release/115.

use std::fmt;
use std::sync::Arc;

use smallvec::SmallVec;

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
/// Variation release/115.
///
/// See <https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html>
/// for the full consequence hierarchy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
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

/// Transcript structural context required by consequence-filtering plugins
/// (LoFTEE) that the summarized [`TranscriptConsequence`] fields do
/// not expose.
///
/// The per-consequence `exon`/`intron` fields are only `"n/total"` summary
/// strings, and the per-transcript exon/intron coordinate model, genomic CDS
/// bounds, and `cds_*_NF` flags never reach the plugin phase (plugins run with
/// `&mut [InputVariant]`, after the source [`crate::transcript::Transcript`] is
/// out of scope). LoFTEE needs the full structure to compute SMALL_INTRON,
/// NON_CAN_SPLICE, GC_TO_GT_DONOR, the GERP-weighted END_TRUNC distance, etc.
///
/// This struct is attached to a [`TranscriptConsequence`] at construction time
/// **only when a consequence-filtering plugin is active** (see
/// `vep_effects::EffectsConfig`'s plugin gate). For all other runs the field stays
/// `None` and is skipped during serialization, so non-LoFTEE output is
/// byte-identical.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
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
    /// Clones the exon list (cheap relative to the per-variant work, and only
    /// done when a consequence-filtering plugin is active). When the
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
/// changes, and metadata from the source transcript. Multiple instances are
/// attached to each [`InputVariant`](crate::variant::InputVariant).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TranscriptConsequence {
    /// Transcript stable ID (e.g., "ENST00000366667").
    pub transcript_id: Arc<str>,
    /// Gene stable ID (e.g., "ENSG00000142168").
    pub gene_id: Arc<str>,
    /// Gene symbol (e.g., "SOD1").
    pub gene_symbol: Option<Arc<str>>,
    /// Gene symbol source (e.g., "HGNC").
    pub gene_symbol_source: Option<Arc<str>>,
    /// HGNC ID.
    pub hgnc_id: Option<Arc<str>>,
    /// Consequence terms assigned.
    pub consequences: ConsequenceList,
    /// Impact of the most severe consequence.
    pub impact: Impact,
    /// Biotype of the transcript.
    pub biotype: Option<Arc<str>>,
    /// Whether this is the canonical transcript.
    pub canonical: bool,
    /// cDNA position (e.g., "428").
    pub cdna_position: Option<String>,
    /// CDS position (e.g., "301").
    pub cds_position: Option<String>,
    /// Protein position (e.g., "101").
    pub protein_position: Option<String>,
    /// Amino acid change (e.g., "A/V").
    pub amino_acids: Option<String>,
    /// Codon change (e.g., "gCc/gTc").
    pub codons: Option<String>,
    /// Protein (Ensembl) ID.
    pub protein_id: Option<Arc<str>>,
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
    /// SIFT prediction and score.
    pub sift: Option<String>,
    /// PolyPhen prediction and score.
    pub polyphen: Option<String>,
    /// Overlapping protein domains.
    pub domains: Vec<(String, String)>,
    /// Feature type (Transcript, RegulatoryFeature, MotifFeature).
    pub feature_type: FeatureType,
    /// Transcript flags (e.g., "cds_start_NF").
    pub flags: Arc<[String]>,
    /// Transcript Support Level (1-5).
    pub tsl: Option<u8>,
    /// MANE Select transcript ID, if this transcript carries the MANE_Select tag.
    /// Used by `--pick`/`--pick_allele`/`flag_pick` to prefer MANE-curated transcripts
    /// over Ensembl-canonical (Perl VEP `pick_order` slot 0).
    pub mane_select: Option<String>,
    /// MANE Plus Clinical transcript ID, if this transcript carries the
    /// MANE_Plus_Clinical tag. Used by `--pick` (Perl VEP `pick_order` slot 1).
    pub mane_plus_clinical: Option<String>,
    /// APPRIS annotation (e.g. "principal1", "alternative2").
    /// Used by `--pick` to prefer principal isoforms (Perl VEP `pick_order` slot 3).
    pub appris: Option<String>,
    /// CCDS ID.
    pub ccds: Option<String>,
    /// UniProtKB/Swiss-Prot accession(s).
    pub swissprot: Option<String>,
    /// UniProtKB/TrEMBL accession(s).
    pub trembl: Option<String>,
    /// RefSeq transcript ID.
    pub refseq: Option<String>,
    /// Plugin annotation data (transcript-consequence-level).
    /// Keyed by field name (e.g., "REVEL_score"), value is the annotation string.
    #[serde(default, skip_serializing_if = "indexmap::IndexMap::is_empty")]
    pub plugin_data: indexmap::IndexMap<String, String>,
    /// Genomic span of the feature (1-based, inclusive), for the structural
    /// variant overlap fields; both zero when unknown.
    #[serde(default)]
    pub feature_start: u64,
    #[serde(default)]
    pub feature_end: u64,
    /// Transcript structural context for consequence-filtering plugins (LoFTEE).
    /// Populated at consequence-construction time only when such a plugin is
    /// active; `None` (and skipped during serialization) otherwise, so default
    /// output is byte-identical. See [`LofteeContext`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub loftee_ctx: Option<Box<LofteeContext>>,
}

/// Type of genomic feature that a consequence is annotated against.
///
/// Most consequences are against `Transcript` features; `RegulatoryFeature`
/// and `MotifFeature` are used for regulatory consequence annotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
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

impl Default for TranscriptConsequence {
    fn default() -> Self {
        Self {
            transcript_id: Arc::from(""),
            gene_id: Arc::from(""),
            gene_symbol: None,
            gene_symbol_source: None,
            hgnc_id: None,
            consequences: SmallVec::new(),
            impact: Impact::MODIFIER,
            biotype: None,
            canonical: false,
            cdna_position: None,
            cds_position: None,
            protein_position: None,
            amino_acids: None,
            codons: None,
            protein_id: None,
            distance: None,
            strand: 1,
            exon: None,
            intron: None,
            hgvsc: None,
            hgvsp: None,
            sift: None,
            polyphen: None,
            domains: Vec::new(),
            feature_type: FeatureType::Transcript,
            flags: Arc::from([]),
            tsl: None,
            mane_select: None,
            mane_plus_clinical: None,
            appris: None,
            ccds: None,
            swissprot: None,
            trembl: None,
            refseq: None,
            plugin_data: indexmap::IndexMap::new(),
            feature_start: 0,
            feature_end: 0,
            loftee_ctx: None,
        }
    }
}

impl TranscriptConsequence {
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
        Some((bp, 100.0 * bp as f64 / feature_len as f64))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript::Exon;

    fn exon(start: u64, end: u64, rank: u32) -> Exon {
        Exon {
            stable_id: None,
            start,
            end,
            rank,
            phase: -1,
            end_phase: -1,
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
}
