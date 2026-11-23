// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The count of an extension in the protein notation of a substitution that
//! changes a stop, on a transcript whose peptide holds a stop and names the
//! residues after it (a polymorphic pseudogene).
//!
//! Each test lays a made-up coding sequence of 23 codons, with or without a
//! stop codon after them, on the first exon of the compact test transcript
//! and holds the row's HGVSp, in some tests its HGVSc too, to values derived
//! from those bases.
//! Ensembl VEP subtracts the peptide's length from the place of the first
//! stop of the alternate sequence's translation
//! (TranscriptVariationAllele.pm:2503 and 2513), and returns no count where
//! the changed residue is the first (2483).

use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture, make_compact_test_transcript};
use crate::hgvs::{generate_hgvs, HgvsNotation, HgvsOptions};
use vep_core::transcript::{SeqEdit, Transcript};
use vep_core::variant::InputVariant;

/// The cDNA position of the first coding base; exon 1 lies at genomic 100 to
/// 399, so coding base `b` is genomic `149 + b`.
const FIRST_CODING_CDNA: u64 = 51;

/// The number of residues the peptide names.
const RESIDUES: usize = 23;

/// The stop codon that closes a coding sequence.
const CLOSING_STOP: &str = "TAG";

/// The 3' UTR after a closing stop codon: two alanine codons, a stop codon
/// and two more bases.
const UTR: [&str; 4] = ["GCT", "GCT", "TGA", "TT"];

/// A coding sequence of 23 codons and the peptide that is its translation: a
/// start codon, alanine codons, and a stop codon at each residue number of
/// `stops`, which the peptide holds as a stop.
fn coding(stops: &[usize]) -> (String, String) {
    let mut bases = String::new();
    let mut peptide = String::new();
    for residue in 1..=RESIDUES {
        let (codon, letter) = if stops.contains(&residue) {
            ("TGA", '*')
        } else if residue == 1 {
            ("ATG", 'M')
        } else {
            ("GCT", 'A')
        };
        bases.push_str(codon);
        peptide.push(letter);
    }
    (bases, peptide)
}

/// The compact test transcript with `cds` laid from its first coding base and
/// `utr` after it, as a polymorphic pseudogene's record gives it: `peptide`
/// and no sequence edit.
fn pseudogene(cds: &str, peptide: &str, utr: Option<&str>) -> Transcript {
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

/// The pseudogene whose peptide holds a stop at each residue number of
/// `stops`, whose coding sequence a stop codon closes, with the 3' UTR after
/// it.
fn closed(stops: &[usize]) -> Transcript {
    let (bases, peptide) = coding(stops);
    pseudogene(
        &format!("{bases}{CLOSING_STOP}"),
        &peptide,
        Some(&UTR.concat()),
    )
}

/// The row of a substitution of coding base `at` of `tx` by `base`, on a
/// reference that holds the transcript's coding sequence and 3' UTR.
fn substitution(tx: &Transcript, at: usize, base: u8) -> HgvsNotation {
    let vefc = tx.vefc.as_ref().expect("a coding model");
    let on_genome = format!(
        "{}{}",
        vefc.translateable_seq
            .as_deref()
            .expect("a coding sequence"),
        vefc.three_prime_utr.as_deref().unwrap_or_default()
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(150, &on_genome));
    let reference = on_genome.as_bytes()[at - 1];
    let variant = InputVariant::new(
        "1".into(),
        149 + at as u64,
        149 + at as u64,
        vec![reference],
        vec![base],
    );
    generate_hgvs(&variant, tx, Some(&fasta), HgvsOptions::default())
}

/// The peptide holds one stop, at residue 10, and names the 13 residues after
/// it. Coding base 30 is the third base of codon 10; with a `G` there the
/// codon is a tryptophan codon and the frame reads on through the alanine
/// codons 11 to 23 to the stop codon that closes the coding sequence, residue
/// 24: 14 residues on from the lost stop. Ensembl VEP finds that stop,
/// subtracts one and the peptide's 23 residues, and writes the unknown count
/// for the nought that leaves.
#[test]
fn a_lost_stop_the_peptide_holds_counts_to_the_stop_codon_that_closes_the_coding_sequence() {
    let row = substitution(&closed(&[10]), 30, b'G');
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.30A>G"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10TrpextTer14")
    );
}

/// The peptide holds stops at residues 10 and 15. With the stop of residue 10
/// lost, the next stop of the frame is the one the peptide holds at residue
/// 15, five residues on. The first stop of Ensembl VEP's alternate
/// translation is that same residue 15, and 15 - 1 - 23 is not positive.
#[test]
fn a_lost_stop_counts_to_the_next_stop_the_peptide_holds() {
    let row = substitution(&closed(&[10, 15]), 30, b'G');
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10TrpextTer5")
    );
}

/// Coding base 45 is the third base of codon 15. With the stop of residue 15
/// lost, the frame reads on to the closing stop codon, residue 24, nine
/// residues on. The first stop of Ensembl VEP's alternate translation is
/// residue 10, before the change.
#[test]
fn a_later_lost_stop_counts_past_the_stop_the_peptide_holds_before_it() {
    let row = substitution(&closed(&[10, 15]), 45, b'G');
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.45A>G"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter15TrpextTer9")
    );
}

/// Coding base 70 is the first base of the closing stop codon, residue 24.
/// With a `C` there the codon is a glutamine codon and the frame reads two
/// alanine codons of the 3' UTR and then its stop codon, residue 27: three
/// residues on. The first stop of Ensembl VEP's alternate translation is the
/// stop the peptide holds at residue 10.
#[test]
fn a_lost_closing_stop_behind_a_stop_the_peptide_holds_counts_into_the_three_prime_utr() {
    let row = substitution(&closed(&[10]), 70, b'C');
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.70T>C"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter24GlnextTer3")
    );
}

/// A coding sequence that no stop codon closes: 23 codons and one base, which
/// the first two bases of the 3' UTR complete to an alanine codon; an alanine
/// codon and a stop codon follow. With the stop of residue 10 lost, the frame
/// meets that stop at residue 26, 16 residues on. Ensembl VEP finds the same
/// stop and counts from the peptide's last residue, 26 - 1 - 23: the two
/// residues past the peptide's end.
#[test]
fn the_count_runs_from_the_lost_stop_and_not_from_the_peptides_last_residue() {
    let (bases, peptide) = coding(&[10]);
    let tx = pseudogene(&format!("{bases}G"), &peptide, Some("CTGCTTAATT"));
    let row = substitution(&tx, 30, b'G');
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10TrpextTer16")
    );
}

/// The stop is sought through the coding sequence and the transcript's 3'
/// UTR and no farther. No stop codon closes this model's coding sequence and
/// the bases of its 3' UTR, read from the reference, are `T` alone, so the
/// frame meets no stop after residue 10 and the count stays unknown.
#[test]
fn the_count_stays_unknown_where_no_stop_follows_in_the_transcript() {
    let (bases, peptide) = coding(&[10]);
    let mut tx = pseudogene(&bases, &peptide, None);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let row = substitution(&tx, 30, b'G');
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10TrpextTer?")
    );
}

/// A peptide that opens with a stop, as a short immunoglobulin segment's
/// does: residue 1 is a stop and so is residue 6. With the stop of residue 1
/// lost the frame meets the stop of residue 6, five residues on. Ensembl VEP
/// hands its count the number of the residue before the changed one, nought,
/// which the routine reads as no argument.
#[test]
fn a_lost_stop_at_the_peptides_first_residue_counts_as_any_other() {
    let mut tx = closed(&[1, 6]);
    tx.flags = vec!["cds_start_NF".to_string()].into();
    let row = substitution(&tx, 3, b'G');
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter1TrpextTer5")
    );
}

/// A transcript whose peptide holds no stop keeps the count both engines
/// write: with the closing stop of residue 24 lost, the stop codon of the 3'
/// UTR stands three residues on.
#[test]
fn a_lost_closing_stop_on_a_peptide_that_holds_no_stop_keeps_its_count() {
    let row = substitution(&closed(&[]), 70, b'C');
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter24GlnextTer3")
    );
}

/// A peptide that holds a stop over a codon that is no stop codon is not the
/// translation of these bases, and the row keeps the unknown count.
#[test]
fn the_count_stays_unknown_where_the_peptide_holds_a_stop_over_another_codon() {
    let (bases, peptide) = coding(&[10]);
    let mut residues = peptide.into_bytes();
    residues[2] = b'*';
    let tx = pseudogene(
        &format!("{bases}{CLOSING_STOP}"),
        &String::from_utf8(residues).expect("ASCII residues"),
        Some(&UTR.concat()),
    );
    let row = substitution(&tx, 30, b'G');
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10TrpextTer?")
    );
}

/// A model that reads the stop codon of residue 5 as a selenocysteine and
/// holds the stop codon of residue 10 as a stop: its peptide is not the
/// translation of its coding sequence at residue 5, so the lost stop of
/// residue 10 keeps the unknown count.
#[test]
fn the_count_stays_unknown_where_the_peptide_also_reads_a_stop_codon_as_a_residue() {
    let (bases, peptide) = coding(&[5, 10]);
    let mut residues = peptide.into_bytes();
    residues[4] = b'U';
    let mut tx = pseudogene(
        &format!("{bases}{CLOSING_STOP}"),
        &String::from_utf8(residues).expect("ASCII residues"),
        Some(&UTR.concat()),
    );
    tx.vefc.as_mut().expect("a coding model").seq_edits = vec![SeqEdit {
        start: 5,
        end: 5,
        alt_seq: "U".to_string(),
    }];
    let row = substitution(&tx, 30, b'G');
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10TrpextTer?")
    );
}
