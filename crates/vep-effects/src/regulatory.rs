// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The consequence terms of a variant against a regulatory or a motif feature
//! and the motif fields of a sequence variant's `MotifFeature` row.
//!
//! The terms are Ensembl VEP's regulatory `OverlapConsequence` entries
//! (`Utils/Constants.pm:1271-1372`), six tier-2 terms whose predicates read the
//! variant's span against the feature's and the variant's class. The motif
//! fields (`MotifFeatureVariationAllele.pm`) are the variant's position in the
//! motif, whether that position is informative and how much the allele changes
//! the motif's binding score, read from the binding matrix the motif names
//! (`vep_core::binding_matrix`).
//!
//! Perl citations name modules of ensembl-variation (`Bio/EnsEMBL/Variation/...`).

use std::sync::Arc;

use smallvec::SmallVec;
use vep_core::binding_matrix::BindingMatrix;
use vep_core::codon::{reverse_complement, seq_is_dna};
use vep_core::consequence::{Consequence, ConsequenceList, FeatureType, MotifFields};
use vep_core::variant::InputVariant;

use crate::sv::{perl_deletion_class, perl_gain_class};

/// The terms of `variant` against a feature of `feature_type` spanning
/// `feature_start..=feature_end`, in rank order, where `start..=end` is the span
/// of the variant the feature was matched on (its own span, or one breakend's
/// position). The overlap term always holds: the feature was found by overlap,
/// which is the `within_feature` predicate (`Utils/VariationEffect.pm:126-144`,
/// `:1653-1656`, `:1684-1687`). A variant covering the whole feature
/// (`complete_overlap_feature`, `:169-178`) adds the ablation term when it is a
/// deletion (`feature_ablation`, `:323-328`) or the amplification term when it
/// is a copy-number gain (`feature_amplification`, `:330-335`). All six terms
/// are tier 2 (`Utils/Constants.pm`), so the tier loop
/// (`BaseVariationFeatureOverlapAllele.pm:265-283`) writes both terms of a
/// covering deletion or gain and reaches no tier-3 term for any regulatory row.
pub fn regulatory_terms(
    variant: &InputVariant,
    start: u64,
    end: u64,
    feature_type: FeatureType,
    feature_start: u64,
    feature_end: u64,
) -> ConsequenceList {
    let (ablation, amplification, overlap) = if feature_type == FeatureType::MotifFeature {
        (
            Consequence::TfbsAblation,
            Consequence::TfbsAmplification,
            Consequence::TfBindingSiteVariant,
        )
    } else {
        (
            Consequence::RegulatoryRegionAblation,
            Consequence::RegulatoryRegionAmplification,
            Consequence::RegulatoryRegionVariant,
        )
    };
    let mut terms: ConsequenceList = SmallVec::new();
    if start <= feature_start && end >= feature_end {
        if perl_deletion_class(variant) {
            terms.push(ablation);
        } else if perl_gain_class(variant) {
            terms.push(amplification);
        }
    }
    terms.push(overlap);
    terms
}

/// A motif feature as the arithmetic reads it.
#[derive(Debug, Clone, Copy)]
pub struct Motif<'a> {
    /// The feature's genomic span, 1-based and inclusive.
    pub start: u64,
    pub end: u64,
    /// The feature's strand, `1` or `-1`.
    pub strand: i8,
    /// The motif's sequence on its own strand, one base per position of the
    /// span (`_motif_feature_seq`, `MotifFeatureVariation.pm:147-152`);
    /// `None` for a source that carries none.
    pub seq: Option<&'a [u8]>,
    pub matrix: &'a BindingMatrix,
}

/// A sequence variant's allele as the arithmetic reads it: the span and the
/// two alleles on the forward strand of the genome, `-` for an empty allele.
#[derive(Debug, Clone, Copy)]
pub struct Allele<'a> {
    pub start: u64,
    pub end: u64,
    pub reference: &'a [u8],
    pub alternate: &'a [u8],
}

impl Motif<'_> {
    /// The number of positions of the feature (`Bio::EnsEMBL::Feature::length`).
    fn feature_length(&self) -> i64 {
        self.end as i64 - self.start as i64 + 1
    }

    /// A genomic position as a position in the motif, 1-based along the motif's
    /// strand: its offset from the feature's start, flipped by the matrix length
    /// on the reverse strand (`motif_start`, `MotifFeatureVariationAllele.pm:99-103`).
    fn motif_position(&self, genomic: u64) -> i64 {
        let relative = genomic as i64 - self.start as i64 + 1;
        if self.strand < 0 {
            self.matrix.length as i64 - relative + 1
        } else {
            relative
        }
    }

    /// `seq`, a forward-strand allele, on the motif's own strand: as given when
    /// the motif is on the forward strand, its reverse complement on the
    /// reverse strand, and a symbolic allele as given (`feature_seq`,
    /// `VariationFeatureOverlapAllele.pm:246-264`, which complements when the
    /// record's strand differs from the motif's; the caller's complement of a
    /// reverse-strand record onto the forward strand composes with this one to
    /// that test).
    fn feature_seq<'s>(&self, seq: &'s [u8]) -> std::borrow::Cow<'s, [u8]> {
        if self.strand < 0 && seq_is_dna(seq) {
            std::borrow::Cow::Owned(reverse_complement(seq))
        } else {
            std::borrow::Cow::Borrowed(seq)
        }
    }
}

/// `motif_start` (`MotifFeatureVariationAllele.pm:90-112`): the variant's
/// genomic start as a position in the motif, `None` outside 1 to the feature's
/// length.
pub fn motif_start(motif: &Motif<'_>, allele: &Allele<'_>) -> Option<u32> {
    let position = motif.motif_position(allele.start);
    (position >= 1 && position <= motif.feature_length()).then_some(position as u32)
}

/// `motif_end` (`MotifFeatureVariationAllele.pm:122-143`): the variant's
/// genomic end as a position in the motif, `None` below 1; a position past the
/// motif's end is kept, for `motif_score_delta` to clip against.
pub fn motif_end(motif: &Motif<'_>, allele: &Allele<'_>) -> Option<i64> {
    let position = motif.motif_position(allele.end);
    (position >= 1).then_some(position)
}

/// `in_informative_position` (`MotifFeatureVariationAllele.pm:153-173`): whether
/// a one-base variant whose allele is not `-` falls on an informative position
/// of the matrix; `false` for any other variant, for a start outside the
/// feature, and for a matrix without counts.
pub fn in_informative_position(motif: &Motif<'_>, allele: &Allele<'_>) -> bool {
    if allele.start != allele.end || allele.alternate == b"-" {
        return false;
    }
    motif_start(motif, allele).is_some_and(|start| motif.matrix.is_position_informative(start))
}

/// `motif_score_delta` (`MotifFeatureVariationAllele.pm:183-265`): the motif's
/// relative similarity score with the allele spliced in, minus the score of the
/// motif's own sequence. `None` when either allele is `-` or the two differ in
/// length (the motif would change length), when the variant's start lies
/// outside the feature or its end before the motif's first position, when the
/// source carries no motif sequence or the matrix no counts, or when the
/// spliced sequence is not the motif's length. An allele running past the
/// motif's end is clipped to the positions inside it (`:222-225`).
///
/// The allele is spliced in at its first position in the motif, the smaller of
/// `motif_start` and `motif_end`. On a reverse-strand motif `motif_start` is
/// the position of the variant's genomic start, the allele's last base along
/// the motif, and Ensembl VEP splices there (`:245`), placing a multi-base
/// allele one or more positions too far along the motif and, at the motif's
/// end, lengthening the sequence so that no score is written; a one-base
/// allele is unaffected. vep-rs scores the allele where it lies.
pub fn motif_score_delta(motif: &Motif<'_>, allele: &Allele<'_>) -> Option<f64> {
    let alternate = motif.feature_seq(allele.alternate);
    let reference = motif.feature_seq(allele.reference);
    if alternate.as_ref() == b"-"
        || reference.as_ref() == b"-"
        || alternate.len() != reference.len()
    {
        return None;
    }
    let start = i64::from(motif_start(motif, allele)?);
    let end = motif_end(motif, allele)?;
    let motif_seq = motif.seq?;
    let motif_len = motif_seq.len() as i64;
    let mut alternate: &[u8] = &alternate;
    if end > motif_len {
        let kept = (motif_len - start + 1).max(0) as usize;
        alternate = &alternate[..kept.min(alternate.len())];
    }
    if alternate.len() as i64 > motif.feature_length() {
        return None;
    }
    let reference_affinity = motif.matrix.relative_similarity(motif_seq)?;
    let first = start.min(end.min(motif_len)) as usize - 1;
    if first + alternate.len() > motif_seq.len() {
        return None;
    }
    let mut spliced = motif_seq.to_vec();
    spliced[first..first + alternate.len()].copy_from_slice(alternate);
    let variant_affinity = motif.matrix.relative_similarity(&spliced)?;
    Some(variant_affinity - reference_affinity)
}

/// The motif fields of a `MotifFeature` row
/// (`MotifFeatureVariationAllele_to_output_hash`, `OutputFactory.pm:1910-1942`):
/// the matrix's id and factors, `MOTIF_POS`, and, for a matrix with counts,
/// `HIGH_INF_POS` and `MOTIF_SCORE_CHANGE`, the score only for an alternate
/// allele of A, C, G and T bases (`:1935-1939`).
pub fn motif_fields(motif: &Motif<'_>, allele: &Allele<'_>) -> MotifFields {
    let scored = motif.matrix.has_frequencies();
    let score_change = if scored && is_plain_dna(allele.alternate) {
        motif_score_delta(motif, allele)
    } else {
        None
    };
    MotifFields {
        name: Arc::clone(&motif.matrix.stable_id),
        pos: motif_start(motif, allele),
        high_inf_pos: scored.then(|| in_informative_position(motif, allele)),
        score_change,
        factors: Arc::clone(&motif.matrix.transcription_factors),
    }
}

/// `/^[ACGT]+$/` on the allele as written (`OutputFactory.pm:1938`).
fn is_plain_dna(allele: &[u8]) -> bool {
    !allele.is_empty() && allele.iter().all(|base| b"ACGT".contains(base))
}

#[cfg(test)]
mod tests {
    use super::*;
    use vep_core::consequence::Impact;
    use vep_core::variant::VariantClass;

    fn sequence(start: u64, end: u64, reference: &[u8], alternate: &[u8]) -> InputVariant {
        InputVariant::new(
            "21".into(),
            start,
            end,
            reference.to_vec(),
            alternate.to_vec(),
        )
    }

    /// A structural record over 100..=130 as the VCF reader builds it: one class,
    /// the record's term, and the alternate allele as written.
    fn structural(class: VariantClass, alt: &str, record_term: Option<&str>) -> InputVariant {
        let mut v = sequence(100, 130, b"N", alt.as_bytes());
        v.is_structural = true;
        v.variant_class = class;
        v.sv_end = Some(130);
        v.record_so_term = record_term.map(str::to_string);
        v
    }

    fn terms(
        variant: &InputVariant,
        feature_type: FeatureType,
        from: u64,
        to: u64,
    ) -> Vec<Consequence> {
        regulatory_terms(variant, variant.start, variant.end, feature_type, from, to).to_vec()
    }

    /// The deletion and gain predicates on the feature-covering spans
    /// (`Utils/VariationEffect.pm:192-216`, `:246-309`) and the class a
    /// copy-number record reads as (ensembl-vep's `Parser.pm:702-706`): a 14-base
    /// deletion over a 12-base motif ablates it, a `<DEL>` inside a regulatory
    /// feature is the overlap term alone, a `<DUP>` over a motif amplifies it, a
    /// lone `<CN0>` is a deletion and a lone `<CN2>` a duplication while
    /// `<CN0>,<CN2>` is a `copy_number_variation` that is neither, and `CAG` to
    /// `CAGCAG` spanning a motif is a tandem repeat, a copy-number gain.
    #[test]
    fn ablation_and_amplification_follow_variation_effect() {
        let deletion = sequence(100, 113, b"ACGTACGTACGTAC", b"-");
        let motif = terms(&deletion, FeatureType::MotifFeature, 101, 112);
        assert_eq!(
            motif,
            [Consequence::TfbsAblation, Consequence::TfBindingSiteVariant]
        );
        assert_eq!(motif[0].impact(), Impact::MODERATE);

        let inside = structural(VariantClass::StructuralDeletion, "<DEL>", Some("deletion"));
        assert_eq!(
            terms(&inside, FeatureType::RegulatoryFeature, 90, 140),
            [Consequence::RegulatoryRegionVariant]
        );
        assert_eq!(
            terms(&inside, FeatureType::RegulatoryFeature, 100, 130),
            [
                Consequence::RegulatoryRegionAblation,
                Consequence::RegulatoryRegionVariant
            ]
        );

        let duplication = structural(VariantClass::Duplication, "<DUP>", Some("duplication"));
        assert_eq!(
            terms(&duplication, FeatureType::MotifFeature, 110, 121),
            [
                Consequence::TfbsAmplification,
                Consequence::TfBindingSiteVariant
            ]
        );

        // The reader classes a lone `<CN0>` as a deletion and a lone `<CN2>` as a
        // duplication; a record read without a term keeps the copy-number class
        // and the allele decides.
        for (class, term) in [
            (VariantClass::StructuralDeletion, Some("deletion")),
            (VariantClass::CopyNumberVariation, None),
        ] {
            let loss = structural(class, "<CN0>", term);
            assert!(perl_deletion_class(&loss) && !perl_gain_class(&loss));
            assert_eq!(
                terms(&loss, FeatureType::RegulatoryFeature, 110, 121)[0],
                Consequence::RegulatoryRegionAblation
            );
        }
        for (class, term) in [
            (VariantClass::Duplication, Some("duplication")),
            (VariantClass::CopyNumberVariation, None),
        ] {
            let gain = structural(class, "<CN2>", term);
            assert!(perl_gain_class(&gain) && !perl_deletion_class(&gain));
            assert_eq!(
                terms(&gain, FeatureType::MotifFeature, 110, 121)[0],
                Consequence::TfbsAmplification
            );
        }
        for alt in ["<CN0>", "<CN2>"] {
            let mixed = structural(
                VariantClass::CopyNumberVariation,
                alt,
                Some("copy_number_variation"),
            );
            assert!(!perl_deletion_class(&mixed) && !perl_gain_class(&mixed));
            assert_eq!(
                terms(&mixed, FeatureType::MotifFeature, 110, 121),
                [Consequence::TfBindingSiteVariant]
            );
        }

        let repeat = sequence(100, 102, b"CAG", b"CAGCAG");
        assert_eq!(
            terms(&repeat, FeatureType::MotifFeature, 100, 102),
            [
                Consequence::TfbsAmplification,
                Consequence::TfBindingSiteVariant
            ]
        );
        let longer = sequence(100, 102, b"CAG", b"CAGCAT");
        assert_eq!(
            terms(&longer, FeatureType::MotifFeature, 100, 102),
            [Consequence::TfBindingSiteVariant]
        );
        let insertion = sequence(101, 100, b"-", b"CAGCAG");
        assert!(!perl_gain_class(&insertion) && !perl_deletion_class(&insertion));
    }

    /// The classes of the other structural terms: an insertion, an inversion, a
    /// tandem duplication, a tandem repeat and a breakpoint
    /// (`Utils/VariationEffect.pm:205-213`, `:268-309`).
    #[test]
    fn structural_classes_read_as_their_terms() {
        let gains = [
            (VariantClass::TandemDuplication, "tandem_duplication"),
            (VariantClass::TandemRepeat, "tandem_repeat"),
        ];
        for (class, term) in gains {
            let v = structural(class, "<DUP:TANDEM>", Some(term));
            assert!(perl_gain_class(&v) && !perl_deletion_class(&v), "{term}");
        }
        let deletion = structural(
            VariantClass::MobileElementDeletion,
            "<DEL:ME:ALU>",
            Some("Alu_deletion"),
        );
        assert!(perl_deletion_class(&deletion) && !perl_gain_class(&deletion));
        let neither = [
            (VariantClass::StructuralInsertion, "insertion"),
            (VariantClass::Inversion, "inversion"),
            (VariantClass::Translocation, "chromosome_breakpoint"),
            (VariantClass::MobileElementInsertion, "Alu_insertion"),
        ];
        for (class, term) in neither {
            let v = structural(class, "<X>", Some(term));
            assert!(!perl_gain_class(&v) && !perl_deletion_class(&v), "{term}");
        }
    }

    fn matrix(id: &str, counts: &[[u32; 4]], factors: &[&str]) -> BindingMatrix {
        BindingMatrix::new(
            id,
            counts.len() as u64,
            counts.iter().map(|row| row.map(f64::from)).collect(),
            factors.iter().map(|f| Arc::from(*f)).collect(),
        )
        .unwrap()
    }

    /// ENSPFM0042 (CTCF, 17 positions) as the Ensembl VEP test cache stores it,
    /// with the motif ENSM00000230511 at 21:25039939-25039955 on the forward
    /// strand, sequence `AGTGCCATCTTCTGTCT`.
    fn ctcf() -> BindingMatrix {
        matrix(
            "ENSPFM0042",
            &[
                [5461, 545, 1745, 4231],
                [2613, 450, 11439, 632],
                [0, 8101, 0, 668],
                [64, 0, 12567, 86],
                [639, 9651, 8, 46],
                [45, 10170, 0, 0],
                [5877, 5220, 207, 241],
                [87, 7119, 0, 3957],
                [0, 8928, 1, 0],
                [39, 48, 23, 14241],
                [10634, 290, 1632, 1371],
                [91, 1625, 9894, 774],
                [186, 5181, 0, 8026],
                [2, 0, 8521, 121],
                [409, 113, 7515, 2252],
                [2682, 4990, 153, 9808],
                [4226, 3598, 2490, 1308],
            ],
            &["CTCF"],
        )
    }

    const MOTIF_SEQ: &[u8] = b"AGTGCCATCTTCTGTCT";
    const START: u64 = 25_039_939;
    const END: u64 = 25_039_955;

    fn motif<'a>(matrix: &'a BindingMatrix, strand: i8, seq: &'a [u8]) -> Motif<'a> {
        Motif {
            start: START,
            end: END,
            strand,
            seq: Some(seq),
            matrix,
        }
    }

    fn allele<'a>(start: u64, end: u64, reference: &'a [u8], alternate: &'a [u8]) -> Allele<'a> {
        Allele {
            start,
            end,
            reference,
            alternate,
        }
    }

    fn printed(delta: Option<f64>) -> Option<String> {
        delta.map(|d| format!("{d:.3}"))
    }

    type Fields = (Option<u32>, Option<bool>, Option<String>);

    fn fields(motif: &Motif<'_>, allele: &Allele<'_>) -> Fields {
        let fields = motif_fields(motif, allele);
        (
            fields.pos,
            fields.high_inf_pos,
            printed(fields.score_change),
        )
    }

    /// The score of `seq` with `bases` written over its positions from `first`
    /// (1-based), minus its own score.
    fn delta_of(matrix: &BindingMatrix, seq: &[u8], first: usize, bases: &[u8]) -> Option<String> {
        let mut spliced = seq.to_vec();
        spliced[first - 1..first - 1 + bases.len()].copy_from_slice(bases);
        printed(Some(
            matrix.relative_similarity(&spliced).unwrap()
                - matrix.relative_similarity(seq).unwrap(),
        ))
    }

    /// The case table of `motif_score_delta` and its companions, on the matrix
    /// and motif Ensembl VEP's own test cache carries. Position 2 G>A prints
    /// `-0.012`, Ensembl VEP's value; position 3 (1.611 bits) is informative.
    /// A multi-base allele inside the motif is scored; one running past the
    /// motif's end is clipped and keeps its MOTIF_POS (`:222-225`); one starting
    /// before the motif has neither MOTIF_POS nor a score (`:107-108`, `:208-211`);
    /// a deletion is `N` with no score; a one-base reference replaced by two bases
    /// is `Y` with no score (`:160` tests the span, `:196-204` the lengths). On a
    /// reverse-strand motif the allele is reverse-complemented and a one-base
    /// change scores at the genomic base's motif position; a multi-base allele is
    /// spliced at its first motif position, where Ensembl VEP splices at its last.
    #[test]
    fn motif_score_delta_cases() {
        let ctcf = ctcf();
        let forward = motif(&ctcf, 1, MOTIF_SEQ);
        // SNV inside: position 2 G>A.
        assert_eq!(
            fields(&forward, &allele(START + 1, START + 1, b"G", b"A")),
            (Some(2), Some(false), Some("-0.012".into()))
        );
        // SNV at the informative position 3.
        let (pos, informative, delta) = fields(&forward, &allele(START + 2, START + 2, b"T", b"C"));
        assert_eq!((pos, informative), (Some(3), Some(true)));
        assert_eq!(delta, delta_of(&ctcf, MOTIF_SEQ, 3, b"C"));
        // MNV inside: positions 2 and 3, `GT` to `AC`.
        assert_eq!(
            fields(&forward, &allele(START + 1, START + 2, b"GT", b"AC")),
            (Some(2), Some(false), delta_of(&ctcf, MOTIF_SEQ, 2, b"AC"))
        );
        // A 3' overhang: positions 16 to 18, the allele clipped to `GG` at 16 and 17.
        assert_eq!(
            fields(&forward, &allele(END - 1, END + 1, b"CTA", b"GGC")),
            (Some(16), Some(false), Some("-0.024".into()))
        );
        // A 5' overhang: the start lies before the motif.
        assert_eq!(
            fields(&forward, &allele(START - 1, START, b"CA", b"TT")),
            (None, Some(false), None)
        );
        // A one-base deletion.
        assert_eq!(
            fields(&forward, &allele(START + 2, START + 2, b"T", b"-")),
            (Some(3), Some(false), None)
        );
        // One reference base replaced by two: informative, no score.
        assert_eq!(
            fields(&forward, &allele(START + 2, START + 2, b"T", b"CA")),
            (Some(3), Some(true), None)
        );
        // An insertion lies between two positions: no base changes.
        assert_eq!(
            fields(&forward, &allele(START + 3, START + 2, b"-", b"A")),
            (Some(4), Some(false), None)
        );

        // The motif on the reverse strand: its stranded sequence is the reverse
        // complement of the forward one, and position 1 is the genomic end.
        let reversed = reverse_complement(MOTIF_SEQ);
        let reverse = motif(&ctcf, -1, &reversed);
        // Position 2 is the genomic base END - 1, `G` on the motif strand and
        // `C` on the genome: G>A on the motif is C>T on the genome.
        assert_eq!(&reversed[1..3], b"GA");
        assert_eq!(
            fields(&reverse, &allele(END - 1, END - 1, b"C", b"T")),
            (Some(2), Some(false), delta_of(&ctcf, &reversed, 2, b"A"))
        );
        // A reverse-strand MNV over motif positions 2 and 3 (genomic END - 2 to
        // END - 1): `GA` to `AC` on the motif is `TC` to `GT` on the genome;
        // MOTIF_POS is the genomic start's position, 3.
        assert_eq!(
            fields(&reverse, &allele(END - 2, END - 1, b"TC", b"GT")),
            (Some(3), Some(false), delta_of(&ctcf, &reversed, 2, b"AC"))
        );
        // A reverse-strand MNV at the motif's end (positions 16 and 17, the
        // genomic start): scored where the allele lies, `GG` over `CT`.
        assert_eq!(&reversed[15..17], b"CT");
        assert_eq!(
            fields(&reverse, &allele(START, START + 1, b"AG", b"CC")),
            (Some(17), Some(false), delta_of(&ctcf, &reversed, 16, b"GG"))
        );
        // Past the genomic end of a reverse-strand motif the end position falls
        // below 1: MOTIF_POS stays, no score.
        assert_eq!(
            fields(&reverse, &allele(END - 1, END + 1, b"CTA", b"GGC")),
            (Some(2), Some(false), None)
        );
        // Before the genomic start of a reverse-strand motif the start position
        // exceeds the motif's length: neither MOTIF_POS nor a score.
        assert_eq!(
            fields(&reverse, &allele(START - 1, START, b"CA", b"TT")),
            (None, Some(false), None)
        );
    }

    /// Two motifs of the GRCh38 cache over one site, 21:5090724 A>T, with the
    /// values Ensembl VEP 116.2 writes for them: the forward motif
    /// ENSM00000574814 (ENSPFM0290, 5090718-5090731) is at an informative
    /// position with the score change `-0.089`, the reverse motif
    /// ENSM00000574815 (ENSPFM0484, 5090723-5090738) is not, with `0.039`.
    #[test]
    fn cache_motifs_reproduce_ensembl_vep_values() {
        let hnf4a = matrix(
            "ENSPFM0290",
            &[
                [3910, 2113, 5352, 1410],
                [5089, 194, 4956, 62],
                [169, 43, 11449, 463],
                [170, 53, 10535, 2299],
                [57, 222, 734, 9359],
                [96, 10935, 834, 2011],
                [9723, 6, 405, 10],
                [9037, 198, 1376, 1038],
                [9604, 4, 1365, 7],
                [52, 1, 11817, 96],
                [54, 8, 10851, 2556],
                [11, 142, 811, 8999],
                [188, 10995, 700, 2137],
                [9604, 148, 1810, 105],
            ],
            &["HNF4A", "NR2F1", "RXRA"],
        );
        let runx3 = matrix(
            "ENSPFM0484",
            &[
                [603, 281, 184, 774],
                [1170, 180, 625, 262],
                [1386, 113, 235, 52],
                [20, 1538, 64, 28],
                [17, 1537, 51, 26],
                [459, 41, 2242, 55],
                [20, 1546, 16, 39],
                [1396, 11, 57, 179],
                [1344, 72, 365, 44],
                [1384, 90, 166, 64],
                [103, 1565, 122, 33],
                [77, 1550, 73, 92],
                [488, 114, 1592, 75],
                [105, 1518, 118, 156],
                [1297, 168, 374, 173],
                [726, 212, 586, 295],
            ],
            &["RUNX3"],
        );
        let site = allele(5_090_724, 5_090_724, b"A", b"T");
        let forward = Motif {
            start: 5_090_718,
            end: 5_090_731,
            strand: 1,
            seq: Some(b"GAGGGCAGAGGTCT"),
            matrix: &hnf4a,
        };
        assert_eq!(
            fields(&forward, &site),
            (Some(7), Some(true), Some("-0.089".into()))
        );
        let reverse = Motif {
            start: 5_090_723,
            end: 5_090_738,
            strand: -1,
            seq: Some(b"TGACCCCAGACCTCTG"),
            matrix: &runx3,
        };
        assert_eq!(
            fields(&reverse, &site),
            (Some(15), Some(false), Some("0.039".into()))
        );
        let fields = motif_fields(&forward, &site);
        assert_eq!(fields.name.as_ref(), "ENSPFM0290");
        assert_eq!(
            fields.factors.iter().map(AsRef::as_ref).collect::<Vec<_>>(),
            ["HNF4A", "NR2F1", "RXRA"]
        );
    }

    /// The fields Ensembl VEP writes only for a matrix with counts
    /// (`OutputFactory.pm:1935-1939`): a matrix known by its length alone
    /// omits HIGH_INF_POS and MOTIF_SCORE_CHANGE and keeps MOTIF_POS; an allele
    /// with a base outside A, C, G and T has no score (`:1938`) but a
    /// HIGH_INF_POS; a motif without a sequence has no score; `motif_end` keeps
    /// a position past the motif and drops one before it.
    #[test]
    fn fields_follow_the_matrix_and_the_allele() {
        let stub = BindingMatrix::new(
            "ENSPFM_STUB",
            17,
            Vec::new(),
            Arc::from([Arc::from("CTCF")]),
        )
        .unwrap();
        let fields_of_stub = motif_fields(
            &motif(&stub, 1, MOTIF_SEQ),
            &allele(START + 2, START + 2, b"T", b"C"),
        );
        assert_eq!(fields_of_stub.name.as_ref(), "ENSPFM_STUB");
        assert_eq!(fields_of_stub.pos, Some(3));
        assert_eq!(fields_of_stub.high_inf_pos, None);
        assert_eq!(fields_of_stub.score_change, None);
        assert_eq!(fields_of_stub.factors.len(), 1);

        let ctcf = ctcf();
        let forward = motif(&ctcf, 1, MOTIF_SEQ);
        assert_eq!(
            fields(&forward, &allele(START + 2, START + 2, b"T", b"N")),
            (Some(3), Some(true), None)
        );
        assert!(motif_score_delta(&forward, &allele(START + 2, START + 2, b"T", b"N")).is_none());
        let unsequenced = Motif {
            seq: None,
            ..forward
        };
        assert_eq!(
            fields(&unsequenced, &allele(START + 2, START + 2, b"T", b"C")),
            (Some(3), Some(true), None)
        );
        assert_eq!(
            motif_end(&forward, &allele(END + 3, END + 3, b"A", b"C")),
            Some(20)
        );
        assert_eq!(
            motif_end(&forward, &allele(START - 1, START - 1, b"A", b"C")),
            None
        );
    }
}
