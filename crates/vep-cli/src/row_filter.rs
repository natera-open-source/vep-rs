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
//! A `--minimal` split record passes the filters twice: once per split feature
//! over its own span, as the rejoin builds each feature's rows before it merges
//! them (`:2289`, `:2303`), and again as the rejoined record. The second pass
//! decides the rows and reads the CDS coordinates cached on each transcript
//! overlap, which the default and tab factories clear after the rejoin
//! (`:297-302`), so they are mapped again over the record's restored span, and
//! which the VCF and JSON factories clear before it (`OutputFactory/VCF.pm:305-308`,
//! `OutputFactory/JSON.pm:217-220`), so the first pass's own-span coordinates
//! stand. The first pass decides whether the first feature wrote a row and so
//! cached the `REF_ALLELE` every row of the record prints (`:1357`,
//! VariationFeature.pm:1533-1539).
//!
//! Perl citations name modules of ensembl-vep release/116.2 (`Bio/EnsEMBL/VEP/...`)
//! and of ensembl-variation release/116 (`Bio/EnsEMBL/Variation/...`).

use vep_core::consequence::{Consequence, FeatureType, TranscriptConsequence};
use vep_core::variant::InputVariant;

/// The two row filters, each a command-line flag, and the span a row of a
/// `--minimal` split answers `affects_cds` over, which follows the output format.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RowFilters {
    pub coding_only: bool,
    pub no_intergenic: bool,
    /// A row of a `--minimal` split is judged over the allele's own span, as
    /// the VCF and JSON factories judge it, rather than over the record's, as
    /// the default and tab factories do.
    pub cds_of_own_span: bool,
}

impl RowFilters {
    /// Whether either filter changes a record.
    pub fn any_active(self) -> bool {
        self.coding_only || self.no_intergenic
    }
}

/// Applies the row filters to the alleles of one input record, in Ensembl VEP's
/// order: `--no_intergenic` (`OutputFactory.pm:482`), then `--coding_only`
/// (`:490-502`). A `--minimal` split record also learns whether its first
/// feature wrote a row before the rejoin
/// ([`InputVariant::split_first_feature_rowless`]).
pub fn apply_to_record(record: &mut [InputVariant], filters: RowFilters) {
    if !filters.any_active() {
        return;
    }
    // `affects_cds` asks the overlap's alternate alleles alone
    // (BaseTranscriptVariation.pm:662-665), so the reference row of
    // `--process_ref_homs` stays only where its feature has an alternate allele:
    // a reference-homozygous copy and a non-variant record lose every overlap,
    // and so does the feature of a `--minimal` split whose allele equals the
    // reference once trimmed, one feature among the record's (the rows of one
    // feature share a rejoin key prefix).
    let feature_has_alternate: Vec<bool> = record
        .iter()
        .map(|row| {
            record.iter().any(|v| {
                !v.reference_allele_row
                    && !v.no_allele
                    && v.rejoin_key_prefix == row.rejoin_key_prefix
            })
        })
        .collect();
    let split = record.iter().any(|v| v.original_start.is_some());
    let first_feature_rowless = split
        && !record
            .iter()
            .zip(&feature_has_alternate)
            .filter(|(v, _)| v.original_start.is_some() && v.rejoin_key_prefix.is_none())
            .any(|(v, has_alternate)| writes_a_row(v, *has_alternate, filters, true));
    for (variant, has_alternate) in record.iter_mut().zip(feature_has_alternate) {
        variant.split_first_feature_rowless = first_feature_rowless;
        if filters.no_intergenic && !overlaps_a_transcript(variant) {
            variant.transcript_consequences.clear();
            variant.intergenic_overlap = false;
            variant.most_severe_consequence = None;
        }
        if filters.coding_only {
            variant.intergenic_overlap = false;
            if variant.reference_allele_row && !has_alternate {
                variant.transcript_consequences.clear();
            }
            let (start, end) = crate::runner::structural_variant_bounds(variant);
            variant
                .transcript_consequences
                .retain(|tc| affects_cds(start, end, tc, filters.cds_of_own_span));
            // An allele left with no row loses its verdict, as the intergenic
            // overlap cannot answer `affects_cds` and goes too (`:493-500`).
            settle_verdict(variant);
        }
    }
}

/// Whether the allele keeps a row under `filters`, `own_span` naming the span a
/// split row answers `affects_cds` over: an intergenic allele writes none under
/// either filter (`--no_intergenic`, `:482`; `--coding_only` drops the one
/// overlap that cannot answer `affects_cds`, `:494-495`); under `--coding_only` a
/// reference row of a feature with no alternate writes none and any other allele
/// needs a transcript row that `affects_cds`. The first pass over a split
/// feature (`OutputFactory.pm:2289`) is this with `own_span` set.
fn writes_a_row(
    variant: &InputVariant,
    has_alternate: bool,
    filters: RowFilters,
    own_span: bool,
) -> bool {
    if variant.no_allele {
        return false;
    }
    if filters.no_intergenic && !overlaps_a_transcript(variant) {
        return false;
    }
    if variant.transcript_consequences.is_empty() {
        return !filters.no_intergenic
            && !filters.coding_only
            && variant.most_severe_consequence.is_some();
    }
    if !filters.coding_only {
        return true;
    }
    if variant.reference_allele_row && !has_alternate {
        return false;
    }
    let (start, end) = crate::runner::structural_variant_bounds(variant);
    variant
        .transcript_consequences
        .iter()
        .any(|tc| affects_cds(start, end, tc, own_span))
}

/// Whether the allele has a row against a transcript, which is what decides
/// Ensembl VEP's intergenic variation (VariationFeature.pm:826-839) and so
/// `--no_intergenic`. An allele with regulatory and motif rows alone overlaps
/// none whether or not it carries the intergenic row
/// ([`InputVariant::intergenic_overlap`]): of a `--minimal` split record, a
/// later feature's reference row keeps those rows and loses the intergenic one
/// (`crate::runner::apply_record_rules`).
fn overlaps_a_transcript(variant: &InputVariant) -> bool {
    variant
        .transcript_consequences
        .iter()
        .any(|tc| tc.feature_type == FeatureType::Transcript)
}

/// Sets the allele's `most_severe_consequence` from the rows it has left: the
/// most severe of their terms, `intergenic_variant` among them for the
/// intergenic row of an allele that overlaps no transcript
/// ([`InputVariant::intergenic_overlap`]), or `None` for an allele with no row,
/// which no writer then prints.
pub(crate) fn settle_verdict(variant: &mut InputVariant) {
    variant.most_severe_consequence = Consequence::most_severe(
        &variant
            .transcript_consequences
            .iter()
            .flat_map(|tc| tc.consequences.iter().copied())
            .chain(
                variant
                    .intergenic_overlap
                    .then_some(Consequence::IntergenicVariant),
            )
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
/// (`Utils/VariationEffect.pm:636-668`) of the alternate allele. Only a
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
///
/// A row of a `--minimal` split has two spans. Over the record's, the CDS
/// coordinate is [`TranscriptConsequence::record_span_in_cds`], which
/// `cds_coords` maps after the rejoin (`BaseTranscriptVariation.pm:505-526`);
/// over its own (`own_span`), it is the own-span CDS_position
/// ([`TranscriptConsequence::own_span_positions`]). The frameshift-intron test
/// reads the intron effects cached on the allele's own span either way
/// (`:656-665`), so the coding term stands for it on both.
pub fn affects_cds(
    variant_start: u64,
    variant_end: u64,
    tc: &TranscriptConsequence,
    own_span: bool,
) -> bool {
    if tc.feature_type != FeatureType::Transcript {
        return false;
    }
    let has_cds_coordinate = match &tc.own_span_positions {
        Some(own) if own_span => own.cds.is_some(),
        Some(_) => tc.record_span_in_cds,
        None => tc.cds_position.is_some(),
    };
    let has_coding_term = tc.consequences.iter().any(|c| CODING_TERMS.contains(c));
    if has_cds_coordinate || has_coding_term {
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
    use vep_core::consequence::{ConsequenceList, Impact, OwnSpanPositions};
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
        cds_of_own_span: false,
    };
    const NO_INTERGENIC: RowFilters = RowFilters {
        coding_only: false,
        no_intergenic: true,
        cds_of_own_span: false,
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

    /// Under `--coding_only` a reference row (`--process_ref_homs`) stays with the
    /// rows of its own feature alone: of a `--minimal` split record `A A,G`, the
    /// feature `-/-` of the allele equal to the reference has no alternate, so its
    /// reference row loses every overlap (`affects_cds` asks the alternates,
    /// BaseTranscriptVariation.pm:662-665), while the feature `A/G`, whose rows
    /// share the key prefix `A/G_`, keeps its reference row beside `G`.
    #[test]
    fn coding_only_judges_a_reference_row_by_its_own_feature() {
        let coding = || {
            vec![tc(
                "ENST1",
                smallvec![Consequence::SynonymousVariant],
                Some("301"),
                Some(900),
            )]
        };
        let mut equal = allele(b"-", 0, coding());
        equal.reference_allele_row = true;
        equal.original_start = Some(5);
        let mut reference = allele(b"A", 0, coding());
        reference.reference_allele_row = true;
        reference.original_start = Some(5);
        reference.rejoin_key_prefix = Some("A/G_".into());
        let mut alternate = allele(b"G", 1, coding());
        alternate.original_start = Some(5);
        alternate.rejoin_key_prefix = Some("A/G_".into());
        let mut record = vec![equal, reference, alternate];
        apply_to_record(&mut record, CODING);
        assert_eq!(
            kept(&record),
            [(0, "ENST1".to_string()), (1, "ENST1".to_string())]
        );
        assert!(record[1].reference_allele_row && !record[1].transcript_consequences.is_empty());
        assert_eq!(record[0].most_severe_consequence, None);
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

    /// A row of a `--minimal` split answers `affects_cds` over the record's span
    /// for the default and tab factories and over its own for the VCF and JSON
    /// ones: of a record whose first allele is an intronic SNV two bases before a
    /// coding exon and whose second is a deletion reaching the CDS, as Ensembl
    /// VEP 116.2 writes `21 10542393 GTACG CTACG,G` under `--minimal
    /// --coding_only`, both rows stay in the default output and the deletion's
    /// alone in the VCF. The SNV's feature, judged on its own span, wrote no row
    /// before the rejoin, so every row of the record is marked for the record's
    /// own `REF_ALLELE`; with the deletion first, that feature wrote one and no
    /// row is marked.
    #[test]
    fn coding_only_judges_a_split_row_over_the_span_the_format_rejoins_on() {
        let split =
            |alt: &[u8], index: usize, mut row: TranscriptConsequence, own_cds: Option<&str>| {
                row.own_span_positions = Some(Box::new(OwnSpanPositions {
                    cdna: None,
                    cds: own_cds.map(str::to_string),
                    protein: None,
                }));
                row.record_span_in_cds = true;
                let mut v = allele(alt, index, vec![row]);
                v.original_start = Some(2999);
                v.original_end = Some(3003);
                if index > 0 {
                    v.rejoin_key_prefix = Some(format!("{}_", String::from_utf8_lossy(alt)).into());
                }
                v
            };
        let snv = || {
            split(
                b"C",
                0,
                tc(
                    "ENST1",
                    smallvec![Consequence::SpliceAcceptorVariant],
                    None,
                    Some(900),
                ),
                None,
            )
        };
        let deletion = |index: usize| {
            split(
                b"-",
                index,
                tc(
                    "ENST1",
                    smallvec![
                        Consequence::SpliceAcceptorVariant,
                        Consequence::CodingSequenceVariant
                    ],
                    Some("?-68"),
                    Some(900),
                ),
                Some("?-68"),
            )
        };
        let sorted = |record: &[InputVariant]| {
            let mut rows = kept(record);
            rows.sort();
            rows
        };

        let mut record = vec![snv(), deletion(1)];
        apply_to_record(&mut record, CODING);
        assert_eq!(
            sorted(&record),
            [(0, "ENST1".to_string()), (1, "ENST1".to_string())]
        );
        assert!(record.iter().all(|v| v.split_first_feature_rowless));

        let mut record = vec![snv(), deletion(1)];
        apply_to_record(
            &mut record,
            RowFilters {
                cds_of_own_span: true,
                ..CODING
            },
        );
        assert_eq!(sorted(&record), [(1, "ENST1".to_string())]);
        assert!(record.iter().all(|v| v.split_first_feature_rowless));

        let mut deletion_first = deletion(0);
        deletion_first.rejoin_key_prefix = None;
        let mut second = snv();
        second.allele_index = 1;
        second.rejoin_key_prefix = Some("G/C_".into());
        let mut record = vec![deletion_first, second];
        apply_to_record(&mut record, CODING);
        assert_eq!(
            sorted(&record),
            [(0, "ENST1".to_string()), (1, "ENST1".to_string())]
        );
        assert!(record.iter().all(|v| !v.split_first_feature_rowless));
    }

    /// `--no_intergenic` judges an allele by its transcript rows alone: Ensembl
    /// VEP's intergenic variation exists where the transcript variations are
    /// none (VariationFeature.pm:826-839), whatever regulatory and motif
    /// variations the feature carries, and `:482` then drops every row. Of a
    /// `--minimal` split record under `--process_ref_homs --regulatory`, the
    /// record rules leave a later feature's reference row with regulatory rows
    /// alone and no intergenic row (`crate::runner::apply_record_rules`); it
    /// loses its regulatory rows as its alternate sibling does, and the
    /// transcript rows of another feature stay.
    #[test]
    fn no_intergenic_drops_the_regulatory_rows_of_an_allele_without_a_transcript_row() {
        let regulatory = || {
            let mut row = tc(
                "ENSR1",
                smallvec![Consequence::RegulatoryRegionVariant],
                None,
                None,
            );
            row.feature_type = FeatureType::RegulatoryFeature;
            row
        };
        let genic = || {
            vec![tc(
                "ENST1",
                smallvec![Consequence::IntronVariant],
                None,
                Some(900),
            )]
        };
        // The reference row and the alternate row of one split feature, as
        // annotation leaves them: with `rows` of a regulatory feature alone,
        // the intergenic row beside the regulatory row on both.
        let feature =
            |alt: &[u8], index: usize, prefix: Option<&str>, rows: Vec<TranscriptConsequence>| {
                let mut alternate = allele(alt, index, rows);
                alternate.original_start = Some(2999);
                alternate.original_end = Some(3001);
                alternate.rejoin_key_prefix = prefix.map(Into::into);
                alternate.intergenic_overlap = alternate
                    .transcript_consequences
                    .iter()
                    .any(|r| r.feature_type != FeatureType::Transcript);
                if !alternate.transcript_consequences.is_empty() {
                    settle_verdict(&mut alternate);
                }
                let mut reference = alternate.clone();
                reference.reference_allele_row = true;
                reference.allele_index = 0;
                [reference, alternate]
            };

        let mut record = feature(b"A", 0, None, Vec::new()).to_vec();
        record.extend(feature(b"AA", 1, Some("-/AA_"), vec![regulatory()]));
        crate::runner::apply_record_rules(&mut record, None);
        assert!(
            record[2].reference_allele_row
                && !record[2].intergenic_overlap
                && record[2].transcript_consequences.len() == 1
                && record[3].intergenic_overlap
        );
        apply_to_record(&mut record, NO_INTERGENIC);
        assert_eq!(kept(&record), []);
        assert!(record.iter().all(|v| v.most_severe_consequence.is_none()));
        assert!(record.iter().all(|v| v.split_first_feature_rowless));

        let mut record = feature(b"A", 0, None, genic()).to_vec();
        record.extend(feature(b"AA", 1, Some("-/AA_"), vec![regulatory()]));
        crate::runner::apply_record_rules(&mut record, None);
        assert!(record[2].intergenic_overlap && record[3].intergenic_overlap);
        apply_to_record(&mut record, NO_INTERGENIC);
        assert_eq!(
            kept(&record),
            [(0, "ENST1".to_string()), (0, "ENST1".to_string())]
        );
        assert_eq!(
            record
                .iter()
                .map(|v| v.most_severe_consequence)
                .collect::<Vec<_>>(),
            [
                Some(Consequence::IntronVariant),
                Some(Consequence::IntronVariant),
                None,
                None
            ]
        );
        assert!(record.iter().all(|v| !v.split_first_feature_rowless));
    }

    /// A `--minimal` split whose first feature overlaps no transcript wrote no
    /// row before the rejoin under either filter: `--coding_only` drops the
    /// intergenic overlap, which cannot answer `affects_cds`
    /// (`OutputFactory.pm:494-495`), as `--no_intergenic` drops it (`:482`), so
    /// the rows of the record's other allele print the record's own `REF_ALLELE`.
    #[test]
    fn a_split_whose_first_feature_is_intergenic_cached_no_reference() {
        let record = || {
            let mut intergenic = allele(b"C", 0, vec![]);
            intergenic.original_start = Some(2999);
            intergenic.original_end = Some(9003);
            let mut deletion = allele(
                b"-",
                1,
                vec![tc(
                    "ENST1",
                    smallvec![Consequence::CodingSequenceVariant],
                    Some("?-68"),
                    Some(900),
                )],
            );
            deletion.original_start = Some(2999);
            deletion.original_end = Some(9003);
            deletion.rejoin_key_prefix = Some("TACG/-_".into());
            vec![intergenic, deletion]
        };
        for filters in [CODING, NO_INTERGENIC] {
            let mut record = record();
            apply_to_record(&mut record, filters);
            assert_eq!(kept(&record), [(1, "ENST1".to_string())], "{filters:?}");
            assert!(
                record.iter().all(|v| v.split_first_feature_rowless),
                "{filters:?}"
            );
        }
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
