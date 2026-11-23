// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The start term of an insertion or deletion that overlaps the start codon
//! and leaves a start codon at the start of the coding sequence, on a
//! transcript with and without a 5' UTR.
//!
//! The Sequence Ontology names `start_retained_variant` a change of at least
//! one base of the start codon after which the start remains (SO:0002019).
//! Behind a 5' UTR Ensembl VEP writes the allele over its span of the UTR and
//! the coding sequence (VariationEffect.pm:1058) and reads the start as kept
//! where the UTR is unchanged and the three bases at the start position read
//! `ATG` (:1061-1065). On a transcript whose coding sequence opens on its
//! first base that test is not made (:1061), every deletion reads as an
//! altered start (:1069) and the row carries `start_lost` where the allele
//! shifts the frame (:870) or the peptide test holds (:878-882): the
//! documented difference `start_codon_kept_without_five_prime_utr`. An
//! insertion there is read by the comparison of the edited sequence's tail
//! with the coding sequence (:1071), which bases added inside the start codon
//! pass only where their last bases repeat the codon's bases before them.
//! vep-rs makes the kept-start test on such a transcript as well, at the
//! transcript's first base.
//!
//! Each test lays one coding sequence on a single-exon transcript with its
//! 5' UTR and without it, on both strands, and holds the terms and the IMPACT
//! of the row [`crate::consequences::calculate_consequences`] gives.

use crate::consequences::{calculate_consequences, EffectsConfig};
use crate::test_helpers::{
    assert_consequence_set_eq, deletion_of_cds, insertion_after_cds, make_single_exon_transcript,
    SINGLE_EXON_UTR5,
};
use std::sync::Arc;
use vep_core::consequence::{Impact, TranscriptConsequence};
use vep_core::coordinate::Strand;
use vep_core::transcript::{ExonCoordMapper, MapperPair, Transcript};
use vep_core::variant::InputVariant;

const STRANDS: [Strand; 2] = [Strand::Forward, Strand::Reverse];

/// A start codon, five more codons and a stop codon, with a second `ATG` out
/// of frame at coding bases 5 to 7.
const SECOND_START_CODON_AT_BASE_FIVE: &str =
    concat!("ATG", "CAT", "GGC", "CAA", "GTT", "TGG", "TAA");

/// A start codon whose third base is repeated by the fourth coding base,
/// four more codons and a stop codon.
const THIRD_BASE_REPEATED: &str = concat!("ATG", "GCC", "AAA", "GTT", "TGG", "TAA");

/// A start codon, a codon that ends on the start codon's third base, a codon
/// that opens with it, two more codons and a stop codon.
const THIRD_BASE_AFTER_ONE_CODON: &str = concat!("ATG", "GCG", "GTT", "TGG", "AAA", "TAA");

/// [`make_single_exon_transcript`] without its 5' UTR: the coding sequence
/// opens on the transcript's first base, so Ensembl VEP's `_five_prime_utr` is
/// undefined (Transcript.pm:1755). The 3' UTR stays.
fn without_five_prime_utr(mut tx: Transcript) -> Transcript {
    let coding_len = tx.cdna_coding_end.expect("coding") - SINGLE_EXON_UTR5;
    let ori = if tx.strand == Strand::Reverse {
        tx.end -= SINGLE_EXON_UTR5;
        -1
    } else {
        tx.start += SINGLE_EXON_UTR5;
        1
    };
    let (start, end) = (tx.start, tx.end);
    let vefc = tx.vefc.as_mut().expect("coding model");
    for exon in tx.exons.iter_mut().chain(vefc.sorted_exons.iter_mut()) {
        exon.start = start;
        exon.end = end;
        exon.phase = 0;
    }
    let mapper = vefc.mapper.as_mut().expect("mapper");
    mapper.cdna_coding_start = 1;
    mapper.cdna_coding_end = coding_len;
    mapper.exon_coord_mapper = ExonCoordMapper::new(vec![MapperPair {
        from_start: 1,
        from_end: end - start + 1,
        to_start: start,
        to_end: end,
        ori,
    }]);
    tx.cdna_coding_start = Some(1);
    tx.cdna_coding_end = Some(coding_len);
    let translation = tx.translation.as_mut().expect("translation");
    translation.start = 1;
    translation.end = coding_len;
    tx
}

fn row(variant: &InputVariant, tx: &Transcript) -> TranscriptConsequence {
    calculate_consequences(variant, &Arc::new(tx.clone()), &EffectsConfig::default())
        .expect("a row")
}

/// The row of the deletion of coding bases `from..=to` of `cds`, laid with a
/// 5' UTR and without one: the genomic bases are the same on both, since the
/// transcript without a UTR keeps the coding region where it stands.
fn rows(strand: Strand, cds: &str, from: usize, to: usize) -> [TranscriptConsequence; 2] {
    let with_utr = make_single_exon_transcript(strand, cds, "");
    let variant = deletion_of_cds(&with_utr, from as u64, to as u64, &cds[from - 1..to]);
    let without_utr = without_five_prime_utr(with_utr.clone());
    [row(&variant, &with_utr), row(&variant, &without_utr)]
}

/// The row of the insertion of `seq` after coding base `after` of `cds`, laid
/// with a 5' UTR and without one.
fn rows_of_insertion(
    strand: Strand,
    cds: &str,
    after: u64,
    seq: &str,
) -> [TranscriptConsequence; 2] {
    let with_utr = make_single_exon_transcript(strand, cds, "");
    let variant = insertion_after_cds(&with_utr, after, seq);
    let without_utr = without_five_prime_utr(with_utr.clone());
    [row(&variant, &with_utr), row(&variant, &without_utr)]
}

/// Four bases removed from the first coding base, the next three reading
/// `ATG`: the row Ensembl VEP writes behind a 5' UTR,
/// `frameshift_variant,start_retained_variant`, on both transcripts.
#[test]
fn a_deletion_from_the_first_coding_base_that_leaves_a_start_codon_keeps_the_start() {
    for strand in STRANDS {
        for tc in rows(strand, SECOND_START_CODON_AT_BASE_FIVE, 1, 4) {
            assert_consequence_set_eq(&tc, &["frameshift_variant", "start_retained_variant"]);
            assert_eq!(tc.impact, Impact::HIGH, "{strand:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("1-4"), "{strand:?}");
            assert_eq!(tc.protein_position.as_deref(), Some("1-2"), "{strand:?}");
        }
    }
}

/// The third base of the start codon removed, the fourth coding base being
/// the same base: the start of the coding sequence still reads `ATG`.
#[test]
fn a_deletion_of_the_start_codons_third_base_before_its_repeat_keeps_the_start() {
    for strand in STRANDS {
        for tc in rows(strand, THIRD_BASE_REPEATED, 3, 3) {
            assert_consequence_set_eq(&tc, &["frameshift_variant", "start_retained_variant"]);
            assert_eq!(tc.impact, Impact::HIGH, "{strand:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("3"), "{strand:?}");
            assert_eq!(tc.protein_position.as_deref(), Some("1"), "{strand:?}");
        }
    }
}

/// Four bases removed from the first coding base, the next three reading no
/// start codon: `start_lost` on both transcripts, as Ensembl VEP writes it.
#[test]
fn a_deletion_from_the_first_coding_base_that_leaves_no_start_codon_loses_the_start() {
    for strand in STRANDS {
        for tc in rows(strand, THIRD_BASE_REPEATED, 1, 4) {
            assert_consequence_set_eq(&tc, &["frameshift_variant", "start_lost"]);
            assert_eq!(tc.impact, Impact::HIGH, "{strand:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("1-4"), "{strand:?}");
        }
    }
}

/// Four bases added after the start codon's first base, opening with the
/// codon's two remaining bases: the codon stands whole before the rest of
/// them, and the tail of the edited sequence is not the coding sequence.
#[test]
fn bases_added_inside_the_start_codon_that_leave_it_whole_keep_the_start() {
    for strand in STRANDS {
        for tc in rows_of_insertion(strand, THIRD_BASE_REPEATED, 1, "TGCC") {
            assert_consequence_set_eq(&tc, &["frameshift_variant", "start_retained_variant"]);
            assert_eq!(tc.impact, Impact::HIGH, "{strand:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("1-2"), "{strand:?}");
        }
    }
}

/// One codon's worth of bases added at the same place: `inframe_insertion`
/// beside the kept-start term, as Ensembl VEP writes it behind a 5' UTR,
/// where it writes `inframe_insertion` alone on a transcript without one.
#[test]
fn a_codon_added_inside_the_start_codon_that_leaves_it_whole_keeps_the_start() {
    for strand in STRANDS {
        for tc in rows_of_insertion(strand, THIRD_BASE_REPEATED, 1, "TGC") {
            assert_consequence_set_eq(&tc, &["inframe_insertion", "start_retained_variant"]);
            assert_eq!(tc.impact, Impact::MODERATE, "{strand:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("1-2"), "{strand:?}");
        }
    }
}

/// One codon's worth of bases removed from the start codon's third base, the
/// base after them repeating it. Ensembl VEP's peptide test names
/// `start_lost` for such a row (:878-882), beside `start_retained_variant`
/// behind a 5' UTR and alone on a transcript without one; the start codon
/// stands, so the row keeps `start_retained_variant` alone on both, with the
/// IMPACT of `inframe_deletion`.
#[test]
fn a_codon_removed_from_the_start_codons_third_base_before_its_repeat_keeps_the_start() {
    for strand in STRANDS {
        for tc in rows(strand, THIRD_BASE_AFTER_ONE_CODON, 3, 5) {
            assert_consequence_set_eq(&tc, &["inframe_deletion", "start_retained_variant"]);
            assert_eq!(tc.impact, Impact::MODERATE, "{strand:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("3-5"), "{strand:?}");
            assert_eq!(tc.protein_position.as_deref(), Some("1-2"), "{strand:?}");
        }
    }
}
