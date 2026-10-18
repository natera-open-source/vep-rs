// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row filters `--no_intergenic` and `--coding_only`, applied to one input
//! record once every allele of it is annotated and before the pick family.
//!
//! Ensembl VEP applies both inside `get_all_VariationFeatureOverlapAlleles`
//! (`OutputFactory.pm:477-516` at release/116.2) and its structural-variant
//! twin (`:534-561`), before the pick family runs over what is left (`:514`,
//! `:560`): a record that overlaps no transcript has no row at all under
//! `--no_intergenic` (`:482`, `:540`), and under `--coding_only` every overlap
//! that cannot answer `affects_cds`, or answers no, is dropped (`:490-502`,
//! `:547-558`), so the intergenic, regulatory and motif overlaps go with the
//! non-coding transcript ones. A record left with no row is written as the
//! writers write an unannotated record.
//!
//! Perl citations name modules of ensembl-vep release/116.2 (`Bio/EnsEMBL/VEP/...`)
//! and of ensembl-variation release/116 (`Bio/EnsEMBL/Variation/...`).

use vep_core::consequence::{Consequence, FeatureType, TranscriptConsequence};
use vep_core::variant::InputVariant;

/// The two row filters, each a command-line flag.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RowFilters {
    pub coding_only: bool,
    pub no_intergenic: bool,
}

impl RowFilters {
    /// Whether either filter changes a record.
    pub fn any_active(self) -> bool {
        self.coding_only || self.no_intergenic
    }
}

/// Applies the row filters to the alleles of one input record, in Ensembl VEP's
/// order: `--no_intergenic` (`OutputFactory.pm:482`), then `--coding_only`
/// (`:490-502`).
pub fn apply_to_record(record: &mut [InputVariant], filters: RowFilters) {
    if !filters.any_active() {
        return;
    }
    for variant in record.iter_mut() {
        if filters.no_intergenic && variant.transcript_consequences.is_empty() {
            variant.most_severe_consequence = None;
        }
        if filters.coding_only {
            let (start, end) = crate::runner::structural_variant_bounds(variant);
            variant
                .transcript_consequences
                .retain(|tc| affects_cds(start, end, tc));
            // An allele left with no row loses its verdict, as the intergenic
            // overlap cannot answer `affects_cds` and goes too (`:493-500`).
            settle_verdict(variant);
        }
    }
}

/// Sets the allele's `most_severe_consequence` from the transcript rows it has
/// left: the most severe of their terms, or `None` for an allele with no row,
/// which no writer then prints.
pub(crate) fn settle_verdict(variant: &mut InputVariant) {
    variant.most_severe_consequence = Consequence::most_severe(
        &variant
            .transcript_consequences
            .iter()
            .flat_map(|tc| tc.consequences.iter().copied())
            .collect::<Vec<_>>(),
    );
}

/// The consequence terms whose predicates Ensembl VEP evaluates only on an
/// allele inside the coding sequence (`Utils/VariationEffect.pm`: each reads the
/// peptide or codon alleles, the start or stop codon, or the CDS coordinates),
/// so a row carrying one is `within_cds` whatever its position columns print;
/// and `within_cds` puts one on every row, `coding_sequence_variant` where no
/// specific term applies (`coding_unknown`, `:1616-1649`).
const CODING_TERMS: [Consequence; 13] = [
    Consequence::StartLost,
    Consequence::StartRetainedVariant,
    Consequence::StopGained,
    Consequence::StopLost,
    Consequence::StopRetainedVariant,
    Consequence::FrameshiftVariant,
    Consequence::InframeInsertion,
    Consequence::InframeDeletion,
    Consequence::MissenseVariant,
    Consequence::ProteinAlteringVariant,
    Consequence::SynonymousVariant,
    Consequence::IncompleteTerminalCodonVariant,
    Consequence::CodingSequenceVariant,
];

/// Whether the row of the allele spanning `variant_start..=variant_end` against
/// `tc` affects the transcript's CDS: `affects_cds`
/// (`BaseTranscriptVariation.pm:662-665`), which is `within_cds`
/// (`Utils/VariationEffect.pm:636-666`) of the alternate allele. Only a
/// transcript answers it; a regulatory or motif row does not.
///
/// `within_cds` holds when the allele has a CDS coordinate (`:645-651`) or lies
/// in a frameshift intron inside the coding region (`:654-665`). The row shows
/// it as a CDS_position (its first or last CDS coordinate) or as a coding term
/// (`CODING_TERMS`): an allele whose ends both fall outside the CDS, across a
/// whole coding exon or from one untranslated exon segment to another, prints no
/// position and carries the term its inner coordinates earn. An allele covering
/// the whole of a translated transcript (`complete_overlap_feature`, `:169-178`)
/// maps onto every CDS exon, so it has CDS coordinates while its columns print
/// none and its terms are the whole-feature ones.
pub fn affects_cds(variant_start: u64, variant_end: u64, tc: &TranscriptConsequence) -> bool {
    if tc.feature_type != FeatureType::Transcript {
        return false;
    }
    let has_coding_term = tc.consequences.iter().any(|c| CODING_TERMS.contains(c));
    if tc.cds_position.is_some() || has_coding_term {
        return true;
    }
    let translated = tc.cds_length.is_some();
    let covers_feature = tc.feature_end > 0
        && tc.feature_start <= tc.feature_end
        && variant_start <= tc.feature_start
        && variant_end >= tc.feature_end;
    translated && covers_feature
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;
    use vep_core::consequence::{ConsequenceList, Impact};
    use vep_io::output::fields::{self, FieldOptions};
    use vep_io::output::json::JsonOutputFormatter;
    use vep_io::output::vcf_output::VcfOutputFormatter;
    use vep_io::output::OutputFormatter;

    /// A transcript row of `ENST<id>` with `consequences`, spanning
    /// `feature_start..=feature_end`, translated when `cds_length` is set.
    fn tc(
        id: &str,
        consequences: ConsequenceList,
        cds_position: Option<&str>,
        cds_length: Option<u32>,
    ) -> TranscriptConsequence {
        TranscriptConsequence {
            transcript_id: id.into(),
            gene_id: "ENSG1".into(),
            consequences,
            impact: Impact::MODIFIER,
            cds_position: cds_position.map(str::to_string),
            cds_length,
            cdna_length: 900,
            feature_start: 1000,
            feature_end: 5000,
            ..Default::default()
        }
    }

    /// An allele at 21:3000 with the rows `tcs`, intergenic when `tcs` is empty.
    fn allele(alt: &[u8], index: usize, tcs: Vec<TranscriptConsequence>) -> InputVariant {
        let mut v = InputVariant::new("21".into(), 3000, 3000, b"A".to_vec(), alt.to_vec());
        v.allele_index = index;
        v.input_record = 1;
        v.raw_input = Some("21\t3000\trec\tA\tG\t.\t.\t.".into());
        v.transcript_consequences = tcs;
        v.most_severe_consequence = Some(Consequence::IntergenicVariant);
        v
    }

    /// `(allele index, transcript id)` of every row left, `-` for an intergenic row.
    fn kept(record: &[InputVariant]) -> Vec<(usize, String)> {
        let refs: Vec<&InputVariant> = record.iter().collect();
        fields::record_rows(&refs)
            .into_iter()
            .map(|(v, tc)| {
                (
                    v.allele_index,
                    tc.map_or("-".to_string(), |t| t.transcript_id.to_string()),
                )
            })
            .collect()
    }

    const CODING: RowFilters = RowFilters {
        coding_only: true,
        no_intergenic: false,
    };
    const NO_INTERGENIC: RowFilters = RowFilters {
        coding_only: false,
        no_intergenic: true,
    };

    /// `--coding_only` keeps the rows with a CDS coordinate and drops the
    /// UTR, intronic and intergenic ones (`OutputFactory.pm:490-502`:
    /// `affects_cds` on the transcript overlaps, nothing for the intergenic
    /// overlap); an allele left with no row has none, where Ensembl VEP's
    /// default run prints its intergenic row.
    #[test]
    fn coding_only_keeps_cds_rows_and_writes_no_intergenic_row() {
        let mut record = vec![
            allele(
                b"G",
                0,
                vec![
                    tc(
                        "ENST1",
                        smallvec![Consequence::MissenseVariant],
                        Some("301"),
                        Some(900),
                    ),
                    tc(
                        "ENST2",
                        smallvec![Consequence::FivePrimeUtrVariant],
                        None,
                        Some(600),
                    ),
                    tc("ENST3", smallvec![Consequence::IntronVariant], None, None),
                ],
            ),
            allele(b"T", 1, vec![]),
        ];
        apply_to_record(&mut record, CODING);
        assert_eq!(kept(&record), [(0, "ENST1".to_string())]);
        assert_eq!(
            record[0].most_severe_consequence,
            Some(Consequence::MissenseVariant)
        );
        assert_eq!(record[1].most_severe_consequence, None);

        let mut utr_only = vec![allele(
            b"G",
            0,
            vec![tc(
                "ENST2",
                smallvec![Consequence::ThreePrimeUtrVariant],
                None,
                Some(600),
            )],
        )];
        apply_to_record(&mut utr_only, CODING);
        assert!(kept(&utr_only).is_empty());
        assert_eq!(utr_only[0].most_severe_consequence, None);
    }

    /// `--no_intergenic` removes the intergenic row of every allele that
    /// overlaps no transcript (`OutputFactory.pm:482`: an empty overlap list for
    /// a record with an `intergenic_variation`), so the record prints as an
    /// unannotated one: no default or tab row, the VCF line with its INFO column
    /// as it came, and in JSON no `intergenic_consequences` key and
    /// `most_severe_consequence` of `?`.
    #[test]
    fn no_intergenic_suppresses_the_intergenic_row_in_every_format() {
        let mut record = vec![allele(b"G", 0, vec![]), allele(b"T", 1, vec![])];
        record[1].raw_input = record[0].raw_input.clone();
        apply_to_record(&mut record, NO_INTERGENIC);
        assert!(kept(&record).is_empty());
        let refs: Vec<&InputVariant> = record.iter().collect();

        let vcf = VcfOutputFormatter::default().format_record(&refs).unwrap();
        assert_eq!(vcf, ["21\t3000\trec\tA\tG\t.\t.\t."]);

        let json = JsonOutputFormatter::new(FieldOptions::default(), Vec::new(), None)
            .format_record_value(&refs);
        assert_eq!(json["most_severe_consequence"], "?");
        assert!(json.get("intergenic_consequences").is_none());

        // A record with a transcript row keeps it: only the intergenic overlap
        // is in question.
        let mut genic = vec![allele(
            b"G",
            0,
            vec![tc(
                "ENST3",
                smallvec![Consequence::IntronVariant],
                None,
                None,
            )],
        )];
        apply_to_record(&mut genic, NO_INTERGENIC);
        assert_eq!(kept(&genic), [(0, "ENST3".to_string())]);
    }

    /// The shapes `within_cds` accepts without a first or last CDS coordinate
    /// (`Utils/VariationEffect.pm:645-651` for an inner coordinate, `:654-665`
    /// for a frameshift intron inside the coding region) reach the row with no
    /// CDS_position and a coding term: `coding_sequence_variant` where no
    /// specific predicate fires, else the specific term alone, as Ensembl VEP
    /// 116.2 writes `start_lost,…,stop_retained_variant,5_prime_UTR_variant,
    /// 3_prime_UTR_variant,intron_variant` with CDS_position `-` for a deletion
    /// from one untranslated exon segment to another, and keeps the row under
    /// `--coding_only`.
    #[test]
    fn coding_only_keeps_a_coding_term_row_without_a_cds_position() {
        let mut record = vec![allele(
            b"-",
            0,
            vec![
                tc(
                    "ENST1",
                    smallvec![
                        Consequence::SpliceAcceptorVariant,
                        Consequence::CodingSequenceVariant,
                        Consequence::IntronVariant
                    ],
                    None,
                    Some(900),
                ),
                tc(
                    "ENST2",
                    smallvec![Consequence::IntronVariant],
                    None,
                    Some(900),
                ),
            ],
        )];
        apply_to_record(&mut record, CODING);
        assert_eq!(kept(&record), [(0, "ENST1".to_string())]);

        let mut across_the_cds = vec![allele(
            b"-",
            0,
            vec![
                tc(
                    "ENST1",
                    smallvec![
                        Consequence::StartLost,
                        Consequence::StopRetainedVariant,
                        Consequence::FivePrimeUtrVariant,
                        Consequence::ThreePrimeUtrVariant,
                        Consequence::IntronVariant
                    ],
                    None,
                    Some(900),
                ),
                tc(
                    "ENST2",
                    smallvec![
                        Consequence::SpliceAcceptorVariant,
                        Consequence::StopRetainedVariant,
                        Consequence::ThreePrimeUtrVariant,
                        Consequence::IntronVariant
                    ],
                    None,
                    Some(600),
                ),
                tc(
                    "ENST3",
                    smallvec![
                        Consequence::ThreePrimeUtrVariant,
                        Consequence::IntronVariant
                    ],
                    None,
                    Some(600),
                ),
            ],
        )];
        apply_to_record(&mut across_the_cds, CODING);
        assert_eq!(
            kept(&across_the_cds),
            [(0, "ENST1".to_string()), (0, "ENST2".to_string())]
        );
    }

    /// An allele covering a transcript whole (`complete_overlap_feature`,
    /// `Utils/VariationEffect.pm:169-178`) has a CDS coordinate for every CDS
    /// exon when the transcript is translated, and none when it is not
    /// (`TranscriptMapper::genomic2cds` returns one gap for a transcript with no
    /// coding start), so Ensembl VEP 116.2 keeps `transcript_ablation`,
    /// `transcript_amplification` and `coding_transcript_variant` rows of a
    /// protein-coding transcript under `--coding_only` and drops those of a
    /// lncRNA.
    #[test]
    fn coding_only_keeps_a_variant_covering_a_translated_transcript_whole() {
        let mut sv = allele(
            b"<DEL>",
            0,
            vec![
                tc(
                    "ENST1",
                    smallvec![Consequence::TranscriptAblation],
                    None,
                    Some(900),
                ),
                tc(
                    "ENST2",
                    smallvec![Consequence::TranscriptAblation],
                    None,
                    None,
                ),
            ],
        );
        sv.start = 900;
        sv.end = 900;
        sv.sv_end = Some(6000);
        sv.is_structural = true;
        let mut record = vec![sv];
        apply_to_record(&mut record, CODING);
        assert_eq!(kept(&record), [(0, "ENST1".to_string())]);

        // A span that stops short of the transcript's end covers no whole
        // feature; with no CDS coordinate and no coding term the row goes.
        let mut partial = allele(
            b"<DEL>",
            0,
            vec![tc(
                "ENST1",
                smallvec![Consequence::FeatureTruncation, Consequence::IntronVariant],
                None,
                Some(900),
            )],
        );
        partial.start = 900;
        partial.end = 900;
        partial.sv_end = Some(4000);
        partial.is_structural = true;
        let mut record = vec![partial];
        apply_to_record(&mut record, CODING);
        assert!(kept(&record).is_empty());
    }

    /// Every overlap is a `BaseVariationFeatureOverlap` (`OutputFactory.pm:494`)
    /// and only a transcript overlap `can('affects_cds')` (`:495`), so a
    /// regulatory or motif row is dropped too.
    #[test]
    fn coding_only_drops_regulatory_and_motif_rows() {
        let mut regulatory = tc(
            "ENSR1",
            smallvec![Consequence::RegulatoryRegionVariant],
            None,
            None,
        );
        regulatory.feature_type = FeatureType::RegulatoryFeature;
        let mut motif = tc(
            "ENSM1",
            smallvec![Consequence::TfBindingSiteVariant],
            None,
            None,
        );
        motif.feature_type = FeatureType::MotifFeature;
        let mut record = vec![allele(
            b"G",
            0,
            vec![
                tc(
                    "ENST1",
                    smallvec![Consequence::SynonymousVariant],
                    Some("12"),
                    Some(300),
                ),
                regulatory,
                motif,
            ],
        )];
        apply_to_record(&mut record, CODING);
        assert_eq!(kept(&record), [(0, "ENST1".to_string())]);
    }

    /// With neither flag the record is untouched.
    #[test]
    fn inactive_filters_leave_the_record_alone() {
        let mut record = vec![
            allele(
                b"G",
                0,
                vec![tc(
                    "ENST3",
                    smallvec![Consequence::IntronVariant],
                    None,
                    None,
                )],
            ),
            allele(b"T", 1, vec![]),
        ];
        let before = kept(&record);
        apply_to_record(&mut record, RowFilters::default());
        assert_eq!(kept(&record), before);
        assert_eq!(
            record[1].most_severe_consequence,
            Some(Consequence::IntergenicVariant)
        );
    }
}
