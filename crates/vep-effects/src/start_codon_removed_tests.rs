// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The start term of a deletion that removes whole codons from the first base
//! of a coding sequence, on a transcript with and without a 5' UTR.
//!
//! The Sequence Ontology names `start_lost` a change of at least one base of
//! the start codon (SO:0002012). A deletion of whole codons from the first
//! coding base leaves the triplet after them at the start of the coding
//! sequence, and where that triplet codes no methionine the start codon is
//! gone. Behind a 5' UTR Ensembl VEP reads the edited sequence:
//! `_inv_start_altered` edits the UTR and the coding sequence and tests the
//! three bases at the old start position against `ATG`
//! (VariationEffect.pm:923-941). On a transcript whose coding sequence opens
//! on its first base the routine returns before the edit (:926), the length
//! test, which reads the start as altered (:1069), is set aside for a
//! deletion of whole codons (:870) and the peptide test returns on the empty
//! alternate peptide (:875), so Ensembl VEP writes `inframe_deletion` alone
//! for the same removal of the same codon: the documented difference
//! `start_codon_removed_without_five_prime_utr`. vep-rs reads the edited
//! sequence on such a transcript as well, where the coding sequence opens on
//! a start codon that codes methionine.
//!
//! Each test lays one coding sequence on a single-exon transcript with its
//! 5' UTR and without it, on both strands, removes whole codons from the
//! first coding base and holds the terms and the IMPACT of the row
//! [`crate::consequences::calculate_consequences`] gives.

use crate::consequences::{calculate_consequences, EffectsConfig};
use crate::test_helpers::{
    assert_consequence_set_eq, deletion_of_cds, make_single_exon_transcript, SINGLE_EXON_UTR5,
};
use std::sync::Arc;
use vep_core::consequence::{Impact, TranscriptConsequence};
use vep_core::coordinate::Strand;
use vep_core::transcript::{ExonCoordMapper, MapperPair, Transcript};
use vep_core::variant::InputVariant;

const STRANDS: [Strand; 2] = [Strand::Forward, Strand::Reverse];

/// A start codon, eight sense codons of which none is `ATG`, and a stop codon.
const CODING_SEQUENCE: &str = "ATGGCCCTGGAAAAGCGTTCCGACGTGTAA";

/// [`CODING_SEQUENCE`] with a second `ATG` as its second codon.
const CODING_SEQUENCE_WITH_TWO_START_CODONS: &str = "ATGATGCTGGAAAAGCGTTCCGACGTGTAA";

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

/// The deletion of the first `codons` codons of `cds`, laid on the transcript
/// with its 5' UTR: the genomic bases are the same on the transcript without
/// one, which keeps the coding region where it stands.
fn removal_of_leading_codons(with_utr: &Transcript, cds: &str, codons: usize) -> InputVariant {
    let bases = codons * 3;
    deletion_of_cds(with_utr, 1, bases as u64, &cds[..bases])
}

fn row(variant: &InputVariant, tx: &Transcript) -> TranscriptConsequence {
    calculate_consequences(variant, &Arc::new(tx.clone()), &EffectsConfig::default())
        .expect("a row")
}

/// The row of the removal on `cds` laid with a 5' UTR and without one.
fn rows(strand: Strand, cds: &str, codons: usize) -> [TranscriptConsequence; 2] {
    let with_utr = make_single_exon_transcript(strand, cds, "");
    let variant = removal_of_leading_codons(&with_utr, cds, codons);
    let without_utr = without_five_prime_utr(with_utr.clone());
    [row(&variant, &with_utr), row(&variant, &without_utr)]
}

/// The start codon removed whole, the next triplet no start codon: the row
/// Ensembl VEP writes behind a 5' UTR, `start_lost,inframe_deletion` with
/// IMPACT `HIGH`, on both transcripts.
#[test]
fn start_codon_removed_whole_is_start_lost_with_and_without_a_five_prime_utr() {
    for strand in STRANDS {
        for tc in rows(strand, CODING_SEQUENCE, 1) {
            assert_consequence_set_eq(&tc, &["start_lost", "inframe_deletion"]);
            assert_eq!(tc.impact, Impact::HIGH, "{strand:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("1-3"), "{strand:?}");
            assert_eq!(tc.protein_position.as_deref(), Some("1"), "{strand:?}");
            assert_eq!(tc.amino_acids.as_deref(), Some("M/-"), "{strand:?}");
            assert_eq!(tc.codons.as_deref(), Some("ATG/-"), "{strand:?}");
        }
    }
}

/// The first two codons removed whole: the same terms, the triplet left at
/// the start of the coding sequence being the third codon.
#[test]
fn two_leading_codons_removed_whole_are_start_lost_without_a_five_prime_utr() {
    for strand in STRANDS {
        for tc in rows(strand, CODING_SEQUENCE, 2) {
            assert_consequence_set_eq(&tc, &["start_lost", "inframe_deletion"]);
            assert_eq!(tc.impact, Impact::HIGH, "{strand:?}");
            assert_eq!(tc.cds_position.as_deref(), Some("1-6"), "{strand:?}");
        }
    }
}

/// Every codon but the last removed from the first coding base: the stop
/// codon is left at the start of the coding sequence.
#[test]
fn leading_codons_removed_up_to_the_last_are_start_lost_without_a_five_prime_utr() {
    let codons = CODING_SEQUENCE.len() / 3 - 1;
    for strand in STRANDS {
        let [_, without_utr] = rows(strand, CODING_SEQUENCE, codons);
        assert!(
            without_utr
                .consequences
                .iter()
                .any(|term| term.to_string() == "start_lost"),
            "{strand:?}: {:?}",
            without_utr.consequences
        );
    }
}

/// A second `ATG` behind the removed start codon stands at the start of the
/// edited coding sequence: no `start_lost` on either transcript. Behind a
/// 5' UTR `_ins_del_start_altered` sees the start codon in place and the row
/// is `start_retained_variant` (VariationEffect.pm:1061-1066). Without one
/// Ensembl VEP's routine reads the shorter sequence as altered (:1069) and
/// writes no start term, where the same test is made at the transcript's
/// first base and the row is `start_retained_variant` as well, the
/// documented difference `start_codon_kept_without_five_prime_utr`.
#[test]
fn start_codon_removed_before_a_second_start_codon_is_not_start_lost() {
    for strand in STRANDS {
        let [with_utr, without_utr] = rows(strand, CODING_SEQUENCE_WITH_TWO_START_CODONS, 1);
        assert_consequence_set_eq(&with_utr, &["start_retained_variant", "inframe_deletion"]);
        assert_consequence_set_eq(
            &without_utr,
            &["start_retained_variant", "inframe_deletion"],
        );
        assert_eq!(without_utr.impact, Impact::MODERATE, "{strand:?}");
    }
}

/// `_overlaps_start_codon` is 0 on a transcript flagged `cds_start_NF`
/// (VariationEffect.pm:975), whose first triplet is no start codon by its
/// annotation: no start term.
#[test]
fn leading_codon_removed_on_a_transcript_incomplete_at_its_start_has_no_start_term() {
    for strand in STRANDS {
        let with_utr = make_single_exon_transcript(strand, CODING_SEQUENCE, "");
        let variant = removal_of_leading_codons(&with_utr, CODING_SEQUENCE, 1);
        let mut incomplete = without_five_prime_utr(with_utr);
        incomplete.flags = vec!["cds_start_NF".to_string()].into();
        assert_consequence_set_eq(&row(&variant, &incomplete), &["inframe_deletion"]);
    }
}

/// One base of the start codon removed: `_ins_del_start_altered` reads the
/// shorter sequence as altered and nothing sets the length test aside
/// (VariationEffect.pm:870, :1069), so both engines write `start_lost` with
/// and without a 5' UTR.
#[test]
fn one_base_of_the_start_codon_removed_is_start_lost_with_and_without_a_five_prime_utr() {
    for strand in STRANDS {
        let with_utr = make_single_exon_transcript(strand, CODING_SEQUENCE, "");
        let variant = deletion_of_cds(&with_utr, 1, 1, "A");
        let without_utr = without_five_prime_utr(with_utr.clone());
        for tx in [&with_utr, &without_utr] {
            assert_consequence_set_eq(&row(&variant, tx), &["frameshift_variant", "start_lost"]);
        }
    }
}

/// [`CODING_SEQUENCE`] opening on `first` and `second` in place of its first
/// two codons.
fn coding_sequence_opening_on(first: &str, second: &str) -> String {
    CODING_SEQUENCE.replacen("ATGGCC", &format!("{first}{second}"), 1)
}

/// The row of the first codon's removal on a transcript without a 5' UTR
/// whose coding sequence opens on `first` and `second`, under the genetic
/// code `table`.
fn row_without_utr_opening_on(
    strand: Strand,
    (first, second): (&str, &str),
    table: u8,
) -> TranscriptConsequence {
    let cds = coding_sequence_opening_on(first, second);
    let with_utr = make_single_exon_transcript(strand, &cds, "");
    let variant = removal_of_leading_codons(&with_utr, &cds, 1);
    let mut without_utr = without_five_prime_utr(with_utr);
    without_utr.vefc.as_mut().expect("coding model").codon_table = table;
    row(&variant, &without_utr)
}

/// The vertebrate mitochondrial code and the invertebrate one, whose
/// transcripts have no 5' UTR: `ATA` is a start codon that codes methionine
/// under both.
const MITOCHONDRIAL_CODES: [u8; 2] = [2, 5];

/// A methionine start codon other than `ATG` removed whole, the triplet left
/// at the head coding another residue: `start_lost`.
#[test]
fn a_mitochondrial_methionine_start_codon_removed_whole_is_start_lost() {
    for strand in STRANDS {
        for table in MITOCHONDRIAL_CODES {
            for first in ["ATA", "ATG"] {
                let tc = row_without_utr_opening_on(strand, (first, "GCC"), table);
                assert_consequence_set_eq(&tc, &["start_lost", "inframe_deletion"]);
                assert_eq!(tc.impact, Impact::HIGH, "{strand:?} {first} {table}");
            }
        }
    }
}

/// The triplet left where the start codon stood codes methionine, as the
/// removed codon did: the removed codon again, `ATG` after `ATA`, or `ATA`
/// after `ATG`. Both engines read a substitution that keeps the first residue
/// as no lost start (VariationEffect.pm:880-881), and the row carries no
/// `start_lost`. Where the triplet left is `ATA` it keeps `inframe_deletion`
/// alone, as Ensembl VEP writes it; where it is `ATG` the kept-start test
/// reads the start codon in place and the row carries
/// `start_retained_variant` as well, the documented difference
/// `start_codon_kept_without_five_prime_utr`. Under the ascidian and the
/// trematode mitochondrial codes `ATA` codes methionine without being a
/// start codon of the code.
#[test]
fn a_methionine_codon_left_where_the_start_codon_stood_is_not_start_lost() {
    for strand in STRANDS {
        for table in MITOCHONDRIAL_CODES {
            for codons in [("ATA", "ATA"), ("ATG", "ATA")] {
                let tc = row_without_utr_opening_on(strand, codons, table);
                assert_consequence_set_eq(&tc, &["inframe_deletion"]);
                assert_eq!(tc.impact, Impact::MODERATE, "{strand:?} {codons:?} {table}");
            }
            let tc = row_without_utr_opening_on(strand, ("ATA", "ATG"), table);
            assert_consequence_set_eq(&tc, &["inframe_deletion", "start_retained_variant"]);
            assert_eq!(tc.impact, Impact::MODERATE, "{strand:?} {table}");
        }
        for table in [13, 21] {
            let tc = row_without_utr_opening_on(strand, ("ATG", "ATA"), table);
            assert_consequence_set_eq(&tc, &["inframe_deletion"]);
            assert_eq!(tc.impact, Impact::MODERATE, "{strand:?} {table}");
        }
    }
}

/// The triplet left where the start codon stood holds a base that is none
/// of `A`, `C`, `G` and `T`: its residue is not known, so no lost start can
/// be read from it.
#[test]
fn an_unknown_base_left_where_the_start_codon_stood_is_not_start_lost() {
    for strand in STRANDS {
        let tc = row_without_utr_opening_on(strand, ("ATG", "NCC"), 1);
        assert!(
            !tc.consequences
                .iter()
                .any(|term| term.to_string() == "start_lost"),
            "{strand:?}: {:?}",
            tc.consequences
        );
    }
}

/// A coding sequence that opens on a codon that is no methionine start codon
/// of its genetic code, on a transcript that carries no `cds_start_NF`:
/// `inframe_deletion` alone, as Ensembl VEP writes it. `GCA` is no start
/// codon; `ATA` is one under the mitochondrial codes and none under the
/// standard code, and under the ascidian and the trematode mitochondrial
/// codes it codes methionine without being a start codon; `CTG` under the
/// standard code and `GTG` under the vertebrate mitochondrial code are start
/// codons that code another residue, which a coding sequence cut off at its
/// start opens on as often as on any other codon.
#[test]
fn a_first_codon_that_is_no_methionine_start_codon_removed_whole_is_not_start_lost() {
    for strand in STRANDS {
        for (first, table) in [
            ("GCA", 1),
            ("ATA", 1),
            ("ATA", 13),
            ("ATA", 21),
            ("CTG", 1),
            ("GTG", 2),
            ("CTG", 2),
        ] {
            let tc = row_without_utr_opening_on(strand, (first, "GCC"), table);
            assert_consequence_set_eq(&tc, &["inframe_deletion"]);
            assert_eq!(tc.impact, Impact::MODERATE, "{strand:?} {first} {table}");
        }
    }
}

/// The whole coding sequence removed on a transcript that holds a 3' UTR: no
/// coding base is left where the start codon stood, and the row carries
/// `start_lost` beside the same terms with a 5' UTR and without one.
#[test]
fn the_whole_coding_sequence_removed_is_start_lost_with_and_without_a_five_prime_utr() {
    let codons = CODING_SEQUENCE.len() / 3;
    for strand in STRANDS {
        let with_utr = make_single_exon_transcript(strand, CODING_SEQUENCE, "CTGATTCC");
        let variant = removal_of_leading_codons(&with_utr, CODING_SEQUENCE, codons);
        let without_utr = without_five_prime_utr(with_utr.clone());
        let terms = |tx: &Transcript| {
            let mut terms: Vec<String> = row(&variant, tx)
                .consequences
                .iter()
                .map(|term| term.to_string())
                .collect();
            terms.sort();
            terms
        };
        let behind_a_utr = terms(&with_utr);
        assert!(
            behind_a_utr.iter().any(|term| term == "start_lost"),
            "{strand:?}: {behind_a_utr:?}"
        );
        assert_eq!(terms(&without_utr), behind_a_utr, "{strand:?}");
    }
}

/// A removal that begins inside the start codon and leaves fewer than three
/// coding bases, on a coding sequence that is no whole number of codons: the
/// alternate peptide is empty there as well, and the span does not start on
/// the first coding base. The rule reads removals of whole codons from that
/// base alone, so the row without a 5' UTR stays as Ensembl VEP writes it,
/// without `start_lost`, where the row with one carries the term.
#[test]
fn a_removal_that_begins_inside_the_start_codon_is_left_as_ensembl_vep_writes_it() {
    let cds = &CODING_SEQUENCE[..10];
    for strand in STRANDS {
        let with_utr = make_single_exon_transcript(strand, cds, "");
        let variant = deletion_of_cds(&with_utr, 2, 10, &cds[1..]);
        let without_utr = without_five_prime_utr(with_utr.clone());
        let names_a_lost_start = |tx: &Transcript| {
            row(&variant, tx)
                .consequences
                .iter()
                .any(|term| term.to_string() == "start_lost")
        };
        assert!(names_a_lost_start(&with_utr), "{strand:?}");
        assert!(!names_a_lost_start(&without_utr), "{strand:?}");
    }
}
