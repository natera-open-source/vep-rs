// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The vocabularies a consequence term is printed in.
//!
//! Ensembl VEP's `--terms` chooses the name of each consequence in the output
//! (`OutputFactory.pm`, `VariationFeatureOverlapAllele_to_output_hash`, where
//! `$self->{terms}.'_term'` names the method called on each overlap
//! consequence): `SO` prints the Sequence Ontology term, `display` the
//! Ensembl display label and `NCBI` the NCBI label, each falling back to the
//! SO term where `%OVERLAP_CONSEQUENCES` (ensembl-variation
//! `Utils/Constants.pm`) carries no label (`OverlapConsequence.pm`,
//! `display_term` and `NCBI_term`). Severity ranks, impacts and the
//! most-severe choice stay on the SO term whatever the vocabulary.

use crate::consequence::Consequence;

/// The vocabulary `--terms` names.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Vocabulary {
    /// Sequence Ontology terms, the default.
    #[default]
    So,
    /// Ensembl display labels (`ESSENTIAL_SPLICE_SITE`, `NON_SYNONYMOUS_CODING`).
    Display,
    /// NCBI labels (`splice-3`, `missense`).
    Ncbi,
}

impl Vocabulary {
    /// The vocabulary a `--terms` value names, matched as Ensembl VEP's
    /// `%VALID` list spells it: `SO`, `display` or `NCBI`.
    pub fn from_flag(value: &str) -> Option<Vocabulary> {
        match value {
            "SO" => Some(Vocabulary::So),
            "display" => Some(Vocabulary::Display),
            "NCBI" => Some(Vocabulary::Ncbi),
            _ => None,
        }
    }

    /// The name of `consequence` in this vocabulary.
    pub fn term(self, consequence: Consequence) -> &'static str {
        match self {
            Vocabulary::So => consequence.so_term(),
            Vocabulary::Display => consequence.display_term(),
            Vocabulary::Ncbi => consequence.ncbi_term(),
        }
    }
}

impl Consequence {
    /// The Ensembl display label, or the SO term for a consequence without one.
    pub fn display_term(self) -> &'static str {
        match self {
            Consequence::SpliceAcceptorVariant | Consequence::SpliceDonorVariant => {
                "ESSENTIAL_SPLICE_SITE"
            }
            Consequence::StopGained => "STOP_GAINED",
            Consequence::FrameshiftVariant => "FRAMESHIFT_CODING",
            Consequence::StopLost => "STOP_LOST",
            Consequence::StartLost
            | Consequence::InframeInsertion
            | Consequence::InframeDeletion
            | Consequence::MissenseVariant => "NON_SYNONYMOUS_CODING",
            Consequence::SpliceDonor5thBaseVariant
            | Consequence::SpliceRegionVariant
            | Consequence::SpliceDonorRegionVariant
            | Consequence::SplicePolypyrimidineTractVariant => "SPLICE_SITE",
            Consequence::IncompleteTerminalCodonVariant => "PARTIAL_CODON",
            Consequence::StartRetainedVariant
            | Consequence::StopRetainedVariant
            | Consequence::SynonymousVariant => "SYNONYMOUS_CODING",
            Consequence::CodingSequenceVariant => "CODING_UNKNOWN",
            Consequence::MatureMirnaVariant => "WITHIN_MATURE_miRNA",
            Consequence::FivePrimeUtrVariant => "5PRIME_UTR",
            Consequence::ThreePrimeUtrVariant => "3PRIME_UTR",
            Consequence::NonCodingTranscriptExonVariant
            | Consequence::NonCodingTranscriptVariant => "WITHIN_NON_CODING_GENE",
            Consequence::IntronVariant => "INTRONIC",
            Consequence::NmdTranscriptVariant => "NMD_TRANSCRIPT",
            Consequence::CodingTranscriptVariant => "CODING_TRANSCRIPT_VARIANT",
            Consequence::UpstreamGeneVariant => "UPSTREAM",
            Consequence::DownstreamGeneVariant => "DOWNSTREAM",
            Consequence::TfBindingSiteVariant | Consequence::RegulatoryRegionVariant => {
                "REGULATORY_REGION"
            }
            Consequence::IntergenicVariant => "INTERGENIC",
            Consequence::SequenceVariant => "SEQUENCE_VARIANT",
            Consequence::TranscriptAblation
            | Consequence::TranscriptAmplification
            | Consequence::FeatureElongation
            | Consequence::FeatureTruncation
            | Consequence::ProteinAlteringVariant
            | Consequence::TfbsAblation
            | Consequence::TfbsAmplification
            | Consequence::RegulatoryRegionAblation
            | Consequence::RegulatoryRegionAmplification => self.so_term(),
        }
    }

    /// The NCBI label, or the SO term for a consequence without one.
    pub fn ncbi_term(self) -> &'static str {
        match self {
            Consequence::SpliceAcceptorVariant => "splice-3",
            Consequence::SpliceDonorVariant => "splice-5",
            Consequence::StopGained => "nonsense",
            Consequence::FrameshiftVariant => "frameshift",
            Consequence::MissenseVariant => "missense",
            Consequence::SynonymousVariant => "cds-synon",
            Consequence::FivePrimeUtrVariant => "untranslated_5",
            Consequence::ThreePrimeUtrVariant => "untranslated_3",
            Consequence::IntronVariant => "intron",
            _ => self.so_term(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `%OVERLAP_CONSEQUENCES` of ensembl-variation `Utils/Constants.pm` at
    /// release/116 (lines 653-1390, plus `$DEFAULT_OVERLAP_CONSEQUENCE` at
    /// 636-650 for `intergenic_variant`): every entry's `SO_term`, its
    /// `display_term` and its `NCBI_term`, an absent key written as `-`.
    const OVERLAP_CONSEQUENCES: &str = "\
transcript_ablation - -
splice_acceptor_variant ESSENTIAL_SPLICE_SITE splice-3
splice_donor_variant ESSENTIAL_SPLICE_SITE splice-5
stop_gained STOP_GAINED nonsense
frameshift_variant FRAMESHIFT_CODING frameshift
stop_lost STOP_LOST -
start_lost NON_SYNONYMOUS_CODING -
transcript_amplification - -
feature_elongation - -
feature_truncation - -
inframe_insertion NON_SYNONYMOUS_CODING -
inframe_deletion NON_SYNONYMOUS_CODING -
missense_variant NON_SYNONYMOUS_CODING missense
protein_altering_variant - -
splice_donor_5th_base_variant SPLICE_SITE -
splice_region_variant SPLICE_SITE -
splice_donor_region_variant SPLICE_SITE -
splice_polypyrimidine_tract_variant SPLICE_SITE -
incomplete_terminal_codon_variant PARTIAL_CODON -
start_retained_variant SYNONYMOUS_CODING -
stop_retained_variant SYNONYMOUS_CODING -
synonymous_variant SYNONYMOUS_CODING cds-synon
coding_sequence_variant CODING_UNKNOWN -
mature_miRNA_variant WITHIN_MATURE_miRNA -
5_prime_UTR_variant 5PRIME_UTR untranslated_5
3_prime_UTR_variant 3PRIME_UTR untranslated_3
non_coding_transcript_exon_variant WITHIN_NON_CODING_GENE -
intron_variant INTRONIC intron
NMD_transcript_variant NMD_TRANSCRIPT -
non_coding_transcript_variant WITHIN_NON_CODING_GENE -
coding_transcript_variant CODING_TRANSCRIPT_VARIANT -
upstream_gene_variant UPSTREAM -
downstream_gene_variant DOWNSTREAM -
TFBS_ablation - -
TFBS_amplification - -
TF_binding_site_variant REGULATORY_REGION -
regulatory_region_ablation - -
regulatory_region_amplification - -
regulatory_region_variant REGULATORY_REGION -
intergenic_variant INTERGENIC -
sequence_variant SEQUENCE_VARIANT -
";

    /// Every consequence renders in each vocabulary as the Perl table says,
    /// the SO term standing in where the table has no label: 32 display
    /// labels and 9 NCBI labels over the 41 terms.
    #[test]
    fn display_and_ncbi_vocabularies_match_ensembl_variation() {
        let mut seen = 0;
        let (mut display_labels, mut ncbi_labels) = (0, 0);
        for line in OVERLAP_CONSEQUENCES.lines() {
            let mut cells = line.split(' ');
            let (so, display, ncbi) = (
                cells.next().unwrap(),
                cells.next().unwrap(),
                cells.next().unwrap(),
            );
            let c = Consequence::from_so_term(so).unwrap_or_else(|| panic!("{so}"));
            let expected = |label: &'static str| if label == "-" { so } else { label };
            assert_eq!(c.display_term(), expected(display), "{so} display_term");
            assert_eq!(c.ncbi_term(), expected(ncbi), "{so} NCBI_term");
            assert_eq!(Vocabulary::So.term(c), so);
            assert_eq!(Vocabulary::Display.term(c), c.display_term());
            assert_eq!(Vocabulary::Ncbi.term(c), c.ncbi_term());
            display_labels += usize::from(display != "-");
            ncbi_labels += usize::from(ncbi != "-");
            seen += 1;
        }
        assert_eq!((seen, display_labels, ncbi_labels), (41, 32, 9));
    }

    /// The `--terms` values are Ensembl VEP's `%VALID` spellings, case included;
    /// the default is the SO vocabulary.
    #[test]
    fn flag_values_name_their_vocabulary() {
        assert_eq!(Vocabulary::from_flag("SO"), Some(Vocabulary::So));
        assert_eq!(Vocabulary::from_flag("display"), Some(Vocabulary::Display));
        assert_eq!(Vocabulary::from_flag("NCBI"), Some(Vocabulary::Ncbi));
        for other in ["so", "ensembl", "Display", "ncbi", ""] {
            assert_eq!(Vocabulary::from_flag(other), None, "{other:?}");
        }
        assert_eq!(Vocabulary::default(), Vocabulary::So);
    }
}
