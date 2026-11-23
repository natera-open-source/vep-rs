// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of an insertion or deletion at the start codon that
//! leaves the whole coding sequence standing behind it.
//!
//! Each test lays a 5' UTR of 50 bases before the coding sequence of the test
//! transcript, which opens `ATG GCT GGA AAA`, and holds the row's terms, its
//! HGVSp and, on the forward strand, its HGVSc to values derived from those
//! bases. Where the bases that end the 5' UTR form the start codon again, the
//! edited transcript is the transcript with a longer or shorter 5' UTR: every
//! codon stands and the notation is the first residue, unchanged. Ensembl VEP
//! writes `Met1?` on each such row, its `start_lost` holding beside
//! `start_retained_variant` (VariationEffect.pm:871 and 961,
//! TranscriptVariationAllele.pm:2161-2164).

use std::sync::Arc;

use crate::concordance_tests::{
    bases_from, hgvsc_fasta_fixture, make_compact_test_transcript,
    make_descending_reverse_strand_transcript,
};
use crate::consequences::{calculate_consequences, EffectsConfig};
use crate::hgvs::{generate_hgvs, HgvsNotation, HgvsOptions};
use crate::test_helpers::{
    cds_ending_in, insertion_after_cds, make_single_exon_transcript, SINGLE_EXON_UTR5,
};
use vep_core::coordinate::Strand;
use vep_core::transcript::{ExonCoordMapper, MapperPair, Transcript};
use vep_core::variant::InputVariant;

/// Genomic position of coding base 1 on the compact forward transcript, whose
/// 5' UTR is genomic 100 to 149.
const FIRST_CODING_BASE: u64 = 150;

/// A 5' UTR of 50 bases that ends in `tail`, eight bases.
fn utr_ending(tail: &str) -> String {
    assert_eq!(tail.len(), 8);
    format!("{}{tail}", "GCCACC".repeat(7))
}

/// The compact forward transcript with a cached 5' UTR that ends in `tail`.
fn forward(tail: &str) -> Transcript {
    let mut tx = make_compact_test_transcript();
    tx.vefc.as_mut().expect("a coding model").five_prime_utr = Some(utr_ending(tail));
    tx
}

/// The reverse-strand transcript with the same cached 5' UTR. Its coding base
/// 1 is genomic 25,005,950 and coding base 2 the base below it.
fn reverse(tail: &str) -> Transcript {
    let mut tx = make_descending_reverse_strand_transcript();
    tx.vefc.as_mut().expect("a coding model").five_prime_utr = Some(utr_ending(tail));
    tx
}

/// The single-exon transcript without its 5' UTR: the transcript opens on
/// coding base 1, genomic 40,000,031.
fn without_a_utr() -> Transcript {
    let mut tx = make_single_exon_transcript(Strand::Forward, &cds_ending_in("GCT", "TAA"), "");
    let coding = cds_ending_in("GCT", "TAA").len() as u64;
    tx.start += SINGLE_EXON_UTR5;
    tx.exons[0].start += SINGLE_EXON_UTR5;
    tx.cdna_coding_start = Some(1);
    tx.cdna_coding_end = Some(coding);
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.sorted_exons = tx.exons.clone();
    let mapper = vefc.mapper.as_mut().expect("a mapper");
    mapper.cdna_coding_start = 1;
    mapper.cdna_coding_end = coding;
    mapper.exon_coord_mapper = ExonCoordMapper::new(vec![MapperPair {
        from_start: 1,
        from_end: coding,
        to_start: tx.start,
        to_end: tx.end,
        ori: 1,
    }]);
    tx.derived = Default::default();
    tx
}

/// The Sequence Ontology terms of the row, sorted.
fn terms(variant: &InputVariant, tx: &Transcript) -> Vec<&'static str> {
    let row = calculate_consequences(variant, &Arc::new(tx.clone()), &EffectsConfig::default())
        .expect("a row on the transcript");
    let mut terms: Vec<_> = row.consequences.iter().map(|term| term.so_term()).collect();
    terms.sort_unstable();
    terms
}

/// The notations of `variant` on the forward transcript, over a reference
/// that holds its 5' UTR and the coding bases of its first exon.
fn notations(variant: &InputVariant, tx: &Transcript) -> HgvsNotation {
    let vefc = tx.vefc.as_ref().expect("a coding model");
    let utr = vefc.five_prime_utr.as_deref().expect("a cached 5' UTR");
    let cds = vefc
        .translateable_seq
        .as_deref()
        .expect("a coding sequence");
    let exon = format!("{utr}{}", &cds[..250]);
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(100, &exon));
    generate_hgvs(variant, tx, Some(&fasta), HgvsOptions::default())
}

/// An insertion after the coding base `after` of the forward transcript.
fn insertion_after(after: u64, bases: &str) -> InputVariant {
    let before = FIRST_CODING_BASE + after - 1;
    InputVariant::new(
        "1".into(),
        before + 1,
        before,
        b"-".to_vec(),
        bases.as_bytes().to_vec(),
    )
}

/// A deletion of the coding bases `first` to `last` of the forward transcript.
fn deletion_of(first: u64, last: u64, bases: &str) -> InputVariant {
    InputVariant::new(
        "1".into(),
        FIRST_CODING_BASE + first - 1,
        FIRST_CODING_BASE + last - 1,
        bases.as_bytes().to_vec(),
        b"-".to_vec(),
    )
}

const FRAMESHIFT_KEPT_START: [&str; 2] = ["frameshift_variant", "start_retained_variant"];
const UNCHANGED: Option<&str> = Some("ENSP00000000001.1:p.Met1=");

/// Four bases written after coding base 1 that repeat the last three bases
/// of the 5' UTR and that base: the copy's last base and the two bases after
/// it are the start codon, four bases later, and the coding sequence follows
/// it whole.
#[test]
fn an_insertion_that_repeats_the_utr_end_and_coding_base_1_changes_no_residue() {
    let tx = forward("GCCACCGC");
    let variant = insertion_after(1, "CGCA");
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = notations(&variant, &tx);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.-3_1dup"));
    assert_eq!(row.hgvsp.as_deref(), UNCHANGED);
}

/// One base written after coding base 1 that repeats it: the copy opens the
/// start codon, one base later.
#[test]
fn a_copy_of_coding_base_1_written_after_it_changes_no_residue() {
    let tx = forward("GCCACCGC");
    let variant = insertion_after(1, "A");
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = notations(&variant, &tx);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.1dup"));
    assert_eq!(row.hgvsp.as_deref(), UNCHANGED);
}

/// Two bases written between coding bases 2 and 3 that repeat bases 1 and 2:
/// the copy and base 3 are the start codon, two bases later.
#[test]
fn a_copy_of_coding_bases_1_and_2_written_after_them_changes_no_residue() {
    let tx = forward("GCCACCGC");
    let variant = insertion_after(2, "AT");
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = notations(&variant, &tx);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.1_2dup"));
    assert_eq!(row.hgvsp.as_deref(), UNCHANGED);
}

/// Coding base 1 deleted behind a 5' UTR that ends in the same base: that
/// base and coding bases 2 and 3 are the start codon, one base earlier.
#[test]
fn a_deletion_of_coding_base_1_behind_a_utr_that_ends_in_it_changes_no_residue() {
    let tx = forward("GCCACCGA");
    let variant = deletion_of(1, 1, "A");
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = notations(&variant, &tx);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.1del"));
    assert_eq!(row.hgvsp.as_deref(), UNCHANGED);
}

/// Coding bases 1 and 2 deleted behind a 5' UTR that ends in the same two
/// bases: those two and coding base 3 are the start codon, two bases earlier.
#[test]
fn a_deletion_of_coding_bases_1_and_2_behind_a_utr_that_ends_in_them_changes_no_residue() {
    let tx = forward("GCCACCAT");
    let variant = deletion_of(1, 2, "AT");
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = notations(&variant, &tx);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.1_2del"));
    assert_eq!(row.hgvsp.as_deref(), UNCHANGED);
}

/// The whole start codon deleted behind a 5' UTR that ends in the start
/// codon's bases: a length of whole codons, and the same reading. The 5' UTR's
/// last three bases are the start codon and the coding sequence follows them.
#[test]
fn a_deletion_of_the_start_codon_behind_a_utr_that_ends_in_it_changes_no_residue() {
    let tx = forward("GCCACATG");
    let variant = deletion_of(1, 3, "ATG");
    assert_eq!(
        terms(&variant, &tx),
        ["inframe_deletion", "start_retained_variant"]
    );
    let row = notations(&variant, &tx);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.1_3del"));
    assert_eq!(row.hgvsp.as_deref(), UNCHANGED);
}

/// The reverse-strand twin of the four-base insertion: the bases the record
/// writes are the reverse complement of the copy, between coding bases 1
/// and 2.
#[test]
fn an_insertion_on_the_reverse_strand_that_repeats_the_utr_end_changes_no_residue() {
    let tx = reverse("GCCACCGC");
    let variant = InputVariant::new(
        "21".into(),
        25_005_950,
        25_005_949,
        b"-".to_vec(),
        b"TGCG".to_vec(),
    );
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = generate_hgvs(&variant, &tx, None, HgvsOptions::default());
    assert_eq!(row.hgvsp.as_deref(), UNCHANGED);
}

/// Four bases written after coding base 1 that open with coding bases 2 and
/// 3 and end in base 1: the start codon stands at the coding start behind
/// the unchanged 5' UTR, and the frame read from it changes from residue 2
/// (`ATG CAT GGC TGG`, then the alanine codons read one base on, which hold
/// no stop codon). Both engines write that frameshift.
#[test]
fn an_insertion_that_leaves_a_start_codon_at_the_coding_start_stays_a_frameshift() {
    let tx = forward("GCCACCGC");
    let variant = insertion_after(1, "TGCA");
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = notations(&variant, &tx);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala2HisfsTer?")
    );
}

/// Coding base 1 deleted behind a 5' UTR that ends in another base: no start
/// codon is formed again, the row carries `start_lost`, and the notation is
/// the lost start's.
#[test]
fn a_deletion_of_coding_base_1_that_forms_no_start_codon_keeps_the_lost_start_notation() {
    let tx = forward("GCCACCGC");
    let variant = deletion_of(1, 1, "A");
    assert_eq!(terms(&variant, &tx), ["frameshift_variant", "start_lost"]);
    let row = notations(&variant, &tx);
    assert_eq!(row.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Met1?"));
}

/// The same reading on a transcript whose coding sequence opens on another
/// codon, `CTG`: four bases written after coding base 1 that open with
/// coding bases 2 and 3 and end in base 1 leave that codon at the coding
/// start behind the unchanged 5' UTR. The frame read from it changes from
/// residue 2 (`CTG ACT GGC TGG`), and the 3' rule writes the same bases after
/// the first codon.
#[test]
fn an_insertion_that_leaves_the_first_codon_of_a_transcript_that_opens_on_another_codon() {
    let mut tx = forward("GCCACCGC");
    tx.vefc
        .as_mut()
        .expect("a coding model")
        .translateable_seq
        .as_mut()
        .expect("a coding sequence")
        .replace_range(..1, "C");
    let variant = insertion_after(1, "TGAC");
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = notations(&variant, &tx);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.3_4insACTG"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala2ThrfsTer?")
    );
}

/// An insertion is read so whatever the 5' UTR holds: the bases it writes end
/// in coding base 1, so the coding sequence follows them whole. The
/// single-exon transcript carries no 5' UTR bases and the run names no
/// reference.
#[test]
fn an_insertion_that_ends_in_coding_base_1_needs_no_base_of_the_utr() {
    let tx = make_single_exon_transcript(Strand::Forward, &cds_ending_in("GCT", "TAA"), "");
    let variant = insertion_after_cds(&tx, 1, "CCCA");
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = generate_hgvs(&variant, &tx, None, HgvsOptions::default());
    assert_eq!(row.hgvsp.as_deref(), UNCHANGED);
}

/// A deletion of the last base of the 5' UTR and coding base 1, behind a 5'
/// UTR whose base before them is coding base 1 again: the coding sequence
/// stands whole, and the row keeps `start_retained_variant`. Its span opens
/// in the 5' UTR and has no translation start, so it carries no protein
/// notation, as in Ensembl VEP (TranscriptVariationAllele.pm:1730-1737).
#[test]
fn a_deletion_across_the_utr_end_and_coding_base_1_has_no_protein_notation() {
    let tx = forward("GCCACCAC");
    let variant = InputVariant::new("1".into(), 149, 150, b"CA".to_vec(), b"-".to_vec());
    assert_eq!(
        terms(&variant, &tx),
        ["5_prime_UTR_variant", "start_retained_variant"]
    );
    let row = notations(&variant, &tx);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.-1_1del"));
    assert_eq!(row.hgvsp, None);
}

/// The coding base after which the two insertions below are written, and the
/// base after it, on the transcript without a 5' UTR.
const NO_UTR_CODING_BASE_1: u64 = 40_000_031;

/// An insertion after coding base 1 of the transcript without a 5' UTR.
fn insertion_after_base_1_without_a_utr(bases: &str) -> InputVariant {
    InputVariant::new(
        "21".into(),
        NO_UTR_CODING_BASE_1 + 1,
        NO_UTR_CODING_BASE_1,
        b"-".to_vec(),
        bases.as_bytes().to_vec(),
    )
}

/// Four bases that end in coding base 1, written after it on a transcript
/// with no 5' UTR: the transcript gains four bases before a coding sequence
/// that stands whole. Ensembl VEP's `start_lost` holds by its peptide test
/// (VariationEffect.pm:872-882), so the row is one of the pair.
#[test]
fn an_insertion_on_a_transcript_with_no_utr_that_keeps_the_coding_sequence_changes_no_residue() {
    let tx = without_a_utr();
    let variant = insertion_after_base_1_without_a_utr("CCCA");
    assert_eq!(terms(&variant, &tx), FRAMESHIFT_KEPT_START);
    let row = generate_hgvs(&variant, &tx, None, HgvsOptions::default());
    assert_eq!(row.hgvsp.as_deref(), UNCHANGED);
}

/// Three such bases on the same transcript are a whole codon before the
/// start codon. Ensembl VEP's `start_lost` does not hold there: with no 5'
/// UTR its test of the old coding start is not made (VariationEffect.pm:926)
/// and the alternate peptide ends in the reference residue (880-882). The
/// row is not one of the pair, the reading of the pair is not made on it,
/// and its notation is left as the insertion's own form gives it: the
/// alternate peptide `TM` is clipped of the `M` it shares with the
/// reference, which leaves `Thr` inserted before residue 1, between a
/// residue 0 and residue 1 that the formatter both names by the residue it
/// reads at position 0, the protein's last, `Ala`.
#[test]
fn a_whole_codon_before_the_start_of_a_transcript_with_no_utr_is_not_read_as_the_pair() {
    let tx = without_a_utr();
    let variant = insertion_after_base_1_without_a_utr("CCA");
    assert_eq!(terms(&variant, &tx), ["start_retained_variant"]);
    let row = generate_hgvs(&variant, &tx, None, HgvsOptions::default());
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala0_Ala1insThr")
    );
}
