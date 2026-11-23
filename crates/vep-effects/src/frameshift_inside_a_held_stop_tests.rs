// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The terms and the protein notation of an insertion or a deletion whose
//! length change is no multiple of three and whose first base lies inside the
//! codon of a stop that the transcript's peptide holds and reads on past (a
//! polymorphic pseudogene).
//!
//! Each test lays a made-up coding sequence of 24 codons on the first exon of
//! the compact test transcript: a start codon, alanine codons, a stop codon
//! at residue 10 that the peptide holds as a stop, ten alanine codons, three
//! tail codons and the closing stop codon, with a 3' UTR after it. The
//! expected terms and notations are derived from those bases. Ensembl VEP's
//! `frameshift` returns 0 wherever the reference peptide of the allele's
//! codons begins with a stop (VariationEffect.pm:1554), which holds for the
//! stop codon no residue follows and not for a stop the peptide reads on
//! past.

use crate::coding::perl_coding_terms;
use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture, make_compact_test_transcript};
use crate::hgvs::{generate_hgvs, HgvsOptions};
use vep_core::consequence::Consequence;
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

/// The cDNA position of the first coding base; exon 1 lies at genomic 100 to
/// 399, so coding base `b` is genomic `149 + b`.
const FIRST_CODING_CDNA: u64 = 51;

/// Codons 21 to 23: threonine, leucine and serine.
const TAIL: [&str; 3] = ["ACC", "CTA", "AGC"];

/// The stop codon that closes the coding sequence, residue 24.
const CLOSING_STOP: &str = "TAG";

/// The 3' UTR: its third to fifth bases are a stop codon in the frame that
/// one added base opens.
const UTR: [&str; 3] = ["CA", "TGA", "TTTT"];

/// The coding sequence: a start codon, eight alanine codons, a stop codon at
/// residue 10, ten alanine codons, the codons of `TAIL` and the closing stop.
fn coding_sequence() -> String {
    format!(
        "ATG{}TGA{}{}{CLOSING_STOP}",
        "GCT".repeat(8),
        "GCT".repeat(10),
        TAIL.concat()
    )
}

/// The compact test transcript with `cds` and `utr` laid from its first
/// coding base and `peptide` as its record's peptide, with no sequence edit.
fn model(cds: &str, utr: Option<&str>, peptide: &str) -> Transcript {
    let mut tx = make_compact_test_transcript();
    let last = FIRST_CODING_CDNA + cds.len() as u64 - 1;
    let genomic_last = 99 + last;
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.translateable_seq = Some(cds.to_string());
    vefc.three_prime_utr = utr.map(str::to_string);
    vefc.peptide = Some(peptide.to_string());
    vefc.seq_edits.clear();
    let mapper = vefc.mapper.as_mut().expect("a mapper");
    mapper.cdna_coding_end = last;
    tx.cdna_coding_end = Some(last);
    tx.coding_region_end = Some(genomic_last);
    tx.translation_end = Some(genomic_last);
    tx
}

/// The model of a polymorphic pseudogene: its peptide holds the stop of
/// residue 10 as a stop and names the 13 residues after it.
fn pseudogene() -> Transcript {
    model(
        &coding_sequence(),
        Some(&UTR.concat()),
        &format!("M{}*{}TLS", "A".repeat(8), "A".repeat(10)),
    )
}

/// The terms and the HGVSp of `variant` on `tx`, on a reference that holds
/// the transcript's coding sequence and 3' UTR.
fn row(tx: &Transcript, variant: &InputVariant) -> (Vec<Consequence>, Option<String>) {
    let vefc = tx.vefc.as_ref().expect("a coding model");
    let on_genome = format!(
        "{}{}",
        vefc.translateable_seq
            .as_deref()
            .expect("a coding sequence"),
        vefc.three_prime_utr.as_deref().unwrap_or_default()
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(150, &on_genome));
    let terms = perl_coding_terms(variant, tx, variant.start, variant.end, Some(&fasta))
        .expect("a row inside the coding sequence");
    let hgvsp = generate_hgvs(variant, tx, Some(&fasta), HgvsOptions::default()).hgvsp;
    (terms, hgvsp)
}

/// The deletion of the coding bases `first` to `last` of `tx`.
fn deletion(tx: &Transcript, first: usize, last: usize) -> InputVariant {
    let vefc = tx.vefc.as_ref().expect("a coding model");
    let cds = vefc
        .translateable_seq
        .as_deref()
        .expect("a coding sequence");
    InputVariant::new(
        "1".into(),
        149 + first as u64,
        149 + last as u64,
        cds.as_bytes()[first - 1..last].to_vec(),
        b"-".to_vec(),
    )
}

/// The insertion of `inserted` between the coding bases `after` and
/// `after + 1`.
fn insertion(after: usize, inserted: &str) -> InputVariant {
    InputVariant::new(
        "1".into(),
        150 + after as u64,
        149 + after as u64,
        b"-".to_vec(),
        inserted.as_bytes().to_vec(),
    )
}

/// Coding base 30 is the last base of the stop codon of residue 10. Without
/// it the codon's two bases and the first base of codon 11 spell a
/// tryptophan codon, so the stop is the first residue that changes, and the
/// frame of every residue the peptide names after it shifts: nine leucine
/// codons, a leucine codon over the end of codon 20, a proline codon and a
/// stop codon at residue 22, twelve residues from the lost stop. Ensembl VEP
/// writes `coding_sequence_variant` and the stop unchanged.
#[test]
fn a_base_lost_inside_a_stop_the_peptide_holds_shifts_the_frame() {
    let tx = pseudogene();
    let (terms, hgvsp) = row(&tx, &deletion(&tx, 30, 30));
    assert_eq!(terms, vec![Consequence::FrameshiftVariant]);
    assert_eq!(
        hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10TrpextTer12")
    );
}

/// One base added between the first and second bases of the stop codon of
/// residue 10 makes the codon a serine codon and shifts every codon after
/// it. The new frame runs past the closing stop and meets its first stop at
/// residue 26, the third to fifth bases of the 3' UTR: sixteen residues from
/// the lost stop. Ensembl VEP writes `coding_sequence_variant` and a
/// deletion-insertion that ends in a stop.
#[test]
fn a_base_added_inside_a_stop_the_peptide_holds_shifts_the_frame() {
    let tx = pseudogene();
    let (terms, hgvsp) = row(&tx, &insertion(28, "C"));
    assert_eq!(terms, vec![Consequence::FrameshiftVariant]);
    assert_eq!(
        hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10SerextTer16")
    );
}

/// The last base of codon 9 and the first base of the stop codon of residue 10
/// are the same base, so a copy of it added between the two, or between the
/// stop codon's first and second bases, leaves one and the same sequence:
/// residue 10 a leucine and the new frame's first stop at residue 26. Ensembl
/// VEP writes the first record as this frameshift and the second, whose
/// first base lies inside the stop codon, as `coding_sequence_variant` with
/// a deletion-insertion.
#[test]
fn one_change_written_on_either_side_of_the_stop_codons_first_base_gets_one_row() {
    let tx = pseudogene();
    let before_the_codon = row(&tx, &insertion(27, "T"));
    let inside_the_codon = row(&tx, &insertion(28, "T"));
    assert_eq!(before_the_codon.0, vec![Consequence::FrameshiftVariant]);
    assert_eq!(
        before_the_codon.1.as_deref(),
        Some("ENSP00000000001.1:p.Ter10LeuextTer16")
    );
    assert_eq!(inside_the_codon, before_the_codon);
}

/// Four bases lost from the second base of the stop codon of residue 10: the
/// codon's first base and two bases of codon 11 spell a leucine codon, and
/// the new frame meets a stop at residue 21. A deletion whose length is no
/// multiple of three is no in-frame deletion.
#[test]
fn four_bases_lost_from_inside_the_stop_are_no_inframe_deletion() {
    let tx = pseudogene();
    let (terms, hgvsp) = row(&tx, &deletion(&tx, 29, 32));
    assert_eq!(terms, vec![Consequence::FrameshiftVariant]);
    assert_eq!(
        hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10LeuextTer11")
    );
}

/// Coding base 29 is the middle base of the stop codon of residue 10.
/// Without it the codon's two other bases and the first base of codon 11
/// spell another stop codon, so residue 10 stays a stop that the peptide
/// reads on past, and the frame of the residues after it shifts all the
/// same: residue 11, an alanine, is the first that changes, to a leucine,
/// and the new frame's stop is residue 22, the twelfth residue counted from
/// the changed one.
#[test]
fn a_stop_the_new_frame_keeps_still_shifts_the_residues_after_it() {
    let tx = pseudogene();
    let (terms, hgvsp) = row(&tx, &deletion(&tx, 29, 29));
    assert_eq!(terms, vec![Consequence::FrameshiftVariant]);
    assert_eq!(
        hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala11LeufsTer12")
    );
}

/// Four bases added between codon 9 and the stop codon of residue 10 spell
/// a stop codon and one base more: residue 10 stays a stop and the frame of
/// the residues after it shifts. The record's first base lies in codon 9,
/// so the row is a frameshift in both engines, and Ensembl VEP's walk for
/// the first changed residue ends at the stop both translations hold
/// (TranscriptVariationAllele.pm:2354-2357). The walk goes on past such a
/// stop only for an allele whose codons begin at it, so this row keeps the
/// notation Ensembl VEP writes.
#[test]
fn a_frameshift_written_before_the_stop_codon_keeps_the_walk_as_it_is() {
    let tx = pseudogene();
    let (terms, hgvsp) = row(&tx, &insertion(27, "TAAC"));
    assert_eq!(
        terms,
        vec![Consequence::StopGained, Consequence::FrameshiftVariant]
    );
    assert_eq!(hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter10="));
}

/// The same four bases added between two stop codons that the peptide
/// holds, of residues 9 and 10. The record's first base lies between the
/// two codons, so its codons hold no residue and begin at no stop: the row
/// is a frameshift in both engines and keeps the notation Ensembl VEP
/// writes, as the row before a lone stop codon does.
#[test]
fn a_frameshift_written_between_two_stop_codons_keeps_the_walk_as_it_is() {
    let cds = format!(
        "ATG{}TGATGA{}{}{CLOSING_STOP}",
        "GCT".repeat(7),
        "GCT".repeat(10),
        TAIL.concat()
    );
    let tx = model(
        &cds,
        Some(&UTR.concat()),
        &format!("M{}**{}TLS", "A".repeat(7), "A".repeat(10)),
    );
    let (terms, hgvsp) = row(&tx, &insertion(27, "TAAC"));
    assert_eq!(
        terms,
        vec![Consequence::StopGained, Consequence::FrameshiftVariant]
    );
    assert_eq!(hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter10="));
}

/// The stop codon that closes the coding sequence has no residue after it,
/// so a base lost inside it shifts no residue's frame: the row stays
/// `stop_lost`, as Ensembl VEP writes it.
#[test]
fn a_base_lost_inside_the_closing_stop_codon_is_no_frameshift() {
    let tx = pseudogene();
    let (terms, _) = row(&tx, &deletion(&tx, 72, 72));
    assert_eq!(terms, vec![Consequence::StopLost]);
}

/// A peptide whose last letter is the stop of residue 10 names no residue
/// after it: the coding sequence ends at that codon, and a base lost inside
/// it shifts no residue's frame.
#[test]
fn a_stop_that_ends_the_peptide_shifts_no_frame() {
    let cds = format!("ATG{}TGA", "GCT".repeat(8));
    let tx = model(&cds, Some(&UTR.concat()), &format!("M{}*", "A".repeat(8)));
    let (terms, _) = row(&tx, &deletion(&tx, 30, 30));
    assert!(!terms.contains(&Consequence::FrameshiftVariant));
}
