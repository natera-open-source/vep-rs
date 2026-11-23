// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The count of an extension in the protein notation of a deletion of whole
//! residues that takes a stop, on a transcript whose peptide holds a stop and
//! names the residues after it (a polymorphic pseudogene).
//!
//! Each test lays a made-up coding sequence of 23 codons, with or without a
//! stop codon after them, on the first exon of the compact test transcript,
//! deletes whole codons of it and holds the row's HGVSp to a value derived
//! from those bases. The count is the place of the first stop the alternate
//! bases meet at or after the first changed residue, less the place of the
//! stop that is lost. Ensembl VEP subtracts one and the peptide's length from
//! the place of the first stop of the whole alternate translation and writes
//! the unknown count for a figure that is not above nought
//! (TranscriptVariationAllele.pm:2503, 2513 and 2519).

use crate::concordance_tests::{bases_from, hgvsc_fasta_fixture, make_compact_test_transcript};
use crate::hgvs::{generate_hgvs, HgvsNotation, HgvsOptions};
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

/// The cDNA position of the first coding base; exon 1 lies at genomic 100 to
/// 399, so coding base `b` is genomic `149 + b`.
const FIRST_CODING_CDNA: u64 = 51;

/// The number of residues the peptide names.
const RESIDUES: usize = 23;

/// The stop codon that closes a coding sequence.
const CLOSING_STOP: &str = "TAG";

/// The 3' UTR after a closing stop codon: two proline codons, a stop codon
/// and two more bases.
const UTR: [&str; 4] = ["CCA", "CCA", "TGA", "TT"];

/// A coding sequence of 23 codons and the peptide that is its translation: a
/// start codon, alanine codons, an aspartate codon at each residue number of
/// `aspartates`, and a stop codon at each residue number of `stops`, which
/// the peptide holds as a stop.
fn coding(stops: &[usize], aspartates: &[usize]) -> (String, String) {
    let mut bases = String::new();
    let mut peptide = String::new();
    for residue in 1..=RESIDUES {
        let (codon, letter) = if stops.contains(&residue) {
            ("TGA", '*')
        } else if residue == 1 {
            ("ATG", 'M')
        } else if aspartates.contains(&residue) {
            ("GAT", 'D')
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
/// `stops` and an aspartate at each of `aspartates`, whose coding sequence a
/// stop codon closes, with the 3' UTR after it.
fn closed(stops: &[usize], aspartates: &[usize]) -> Transcript {
    let (bases, peptide) = coding(stops, aspartates);
    pseudogene(
        &format!("{bases}{CLOSING_STOP}"),
        &peptide,
        Some(&UTR.concat()),
    )
}

/// The row of a deletion of the codons of residues `first` to `last` of `tx`,
/// on a reference that holds the transcript's coding sequence and 3' UTR.
fn deletion(tx: &Transcript, first: usize, last: usize) -> HgvsNotation {
    let vefc = tx.vefc.as_ref().expect("a coding model");
    let on_genome = format!(
        "{}{}",
        vefc.translateable_seq
            .as_deref()
            .expect("a coding sequence"),
        vefc.three_prime_utr.as_deref().unwrap_or_default()
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(150, &on_genome));
    let (from, to) = (3 * (first - 1) + 1, 3 * last);
    let variant = InputVariant::new(
        "1".into(),
        149 + from as u64,
        149 + to as u64,
        on_genome.as_bytes()[from - 1..to].to_vec(),
        b"-".to_vec(),
    );
    generate_hgvs(&variant, tx, Some(&fasta), HgvsOptions::default())
}

/// The peptide holds one stop, at residue 10, after an aspartate. With the
/// codons of residues 9 and 10 deleted the alanine of residue 11 stands at
/// residue 9, the first changed residue, and the frame reads on to the stop
/// codon that closes the coding sequence, which stood at residue 24 and
/// stands at residue 22: 12 residues on from the lost stop. Ensembl VEP
/// finds that stop, subtracts one and the peptide's 23 residues, and writes
/// the unknown count for the figure below nought.
#[test]
fn a_deletion_through_a_stop_the_peptide_holds_counts_to_the_stop_codon_that_closes_the_coding_sequence(
) {
    let row = deletion(&closed(&[10], &[9]), 9, 10);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp9_Ter10delextTer12")
    );
}

/// A first residue written as a methionine over a first codon that codes
/// another residue, here a leucine codon, leaves the peptide the
/// translation of its coding sequence: the count is the one of the test
/// above.
#[test]
fn a_methionine_written_over_another_first_codon_keeps_the_count() {
    let (bases, peptide) = coding(&[10], &[9]);
    let tx = pseudogene(
        &format!("CTG{}{CLOSING_STOP}", &bases[3..]),
        &peptide,
        Some(&UTR.concat()),
    );
    let row = deletion(&tx, 9, 10);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp9_Ter10delextTer12")
    );
}

/// The peptide holds stops at residues 10 and 15. With the codons of
/// residues 10 and 11 deleted the first changed residue is the stop itself,
/// and the stop of residue 15 stands at residue 13: three residues on. The
/// first stop of Ensembl VEP's alternate translation is that same stop, and
/// 13 - 1 - 23 is not above nought.
#[test]
fn a_deletion_that_opens_at_a_stop_counts_to_the_next_stop_the_peptide_holds() {
    let row = deletion(&closed(&[10, 15], &[]), 10, 11);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10_Ala11delextTer3")
    );
}

/// The stop alone deleted: the alanine of residue 11 stands in its place and
/// the stop of residue 15 at residue 14, four residues on. The description
/// keeps the form the codon window gives it, the stop written as deleted.
#[test]
fn a_deletion_of_the_stop_alone_counts_to_the_next_stop_the_peptide_holds() {
    let row = deletion(&closed(&[10, 15], &[]), 10, 10);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter10delextTer4")
    );
}

/// The peptide holds stops at residues 5, 10 and 15. With the codons of
/// residues 9 and 10 deleted the stop of residue 15 stands at residue 13,
/// three residues on from the lost stop. The first stop of Ensembl VEP's
/// alternate translation is residue 5, before the change.
#[test]
fn a_deletion_counts_past_the_stop_the_peptide_holds_before_it() {
    let row = deletion(&closed(&[5, 10, 15], &[9]), 9, 10);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp9_Ter10delextTer3")
    );
}

/// The stop codon that closes the coding sequence, residue 24, behind a stop
/// the peptide holds at residue 10. With the codons of residues 23 and 24
/// deleted the first proline of the 3' UTR stands at residue 23 and the stop
/// codon of the 3' UTR at residue 25, one residue on from the lost stop. The
/// first stop of Ensembl VEP's alternate translation is residue 10.
#[test]
fn a_deletion_through_the_closing_stop_behind_a_stop_the_peptide_holds_counts_into_the_three_prime_utr(
) {
    let row = deletion(&closed(&[10], &[23]), 23, 24);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp23_Ter24delextTer1")
    );
}

/// A coding sequence of 23 codons that no stop codon closes; four alanine
/// codons and a stop codon follow in the 3' UTR. With the codons of residues
/// 19 and 20 deleted the stop of the 3' UTR stands at residue 26, six
/// residues on from the lost stop of residue 20. Ensembl VEP finds the same
/// stop and counts from the peptide's last residue, 26 - 1 - 23: two.
#[test]
fn the_count_runs_from_the_lost_stop_and_not_from_the_peptides_last_residue() {
    let (bases, peptide) = coding(&[20], &[19]);
    let tx = pseudogene(&bases, &peptide, Some("GCTGCTGCTGCTTGATT"));
    let row = deletion(&tx, 19, 20);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp19_Ter20delextTer6")
    );
}

/// A peptide that opens with a stop, as a short immunoglobulin segment's
/// does: residue 1 is a stop and so is residue 6. With the first codon
/// deleted the alanine of residue 2 stands at residue 1 and the stop of
/// residue 6 at residue 5, four residues on. Ensembl VEP hands its count the
/// number of the residue before the description's first one, nought, which
/// the routine reads as no argument (TranscriptVariationAllele.pm:2483).
#[test]
fn a_deletion_that_opens_at_the_peptides_first_residue_counts_as_any_other() {
    let mut tx = closed(&[1, 6], &[]);
    tx.flags = vec!["cds_start_NF".to_string()].into();
    let row = deletion(&tx, 1, 1);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter1delextTer4")
    );
}

/// The peptide holds stops at residues 10 and 12. With the codons of
/// residues 9 and 10 deleted the stop of residue 12 stands at residue 10,
/// the lost stop's own residue: the protein ends where it ended, and the
/// row keeps the unknown count both engines write.
#[test]
fn the_count_stays_unknown_where_the_next_stop_stands_at_the_lost_stops_residue() {
    let row = deletion(&closed(&[10, 12], &[9]), 9, 10);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp9_Ter10delextTer?")
    );
}

/// The peptide holds stops at residues 10 and 11. With the codons of
/// residues 9 and 10 deleted the stop of residue 11 stands at residue 9, the
/// first changed residue, before the lost stop's residue: the row keeps the
/// unknown count both engines write.
#[test]
fn the_count_stays_unknown_where_a_stop_stands_at_the_first_changed_residue() {
    let row = deletion(&closed(&[10, 11], &[9]), 9, 10);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp9_Ter10delextTer?")
    );
}

/// The stop is sought through the coding sequence and the transcript's 3'
/// UTR and no farther. No stop codon closes this model's coding sequence and
/// the bases of its 3' UTR, read from the reference, are `T` alone, so the
/// frame meets no stop after residue 19 and the count stays unknown.
#[test]
fn the_count_stays_unknown_where_no_stop_follows_in_the_transcript() {
    let (bases, peptide) = coding(&[20], &[19]);
    let mut tx = pseudogene(&bases, &peptide, None);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let row = deletion(&tx, 19, 20);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp19_Ter20delextTer?")
    );
}

/// A transcript whose peptide holds no stop keeps the count both engines
/// write: with the codons of residue 23 and of the closing stop deleted, the
/// stop codon of the 3' UTR stands at residue 25, one residue on from the
/// lost stop, and 25 - 1 - 23 is the same figure.
#[test]
fn a_deletion_through_the_closing_stop_of_a_peptide_that_holds_no_stop_keeps_its_count() {
    let row = deletion(&closed(&[], &[23]), 23, 24);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp23_Ter24delextTer1")
    );
}

/// A peptide that holds a stop over a codon that is no stop codon is not the
/// translation of these bases, and the row keeps the unknown count.
#[test]
fn the_count_stays_unknown_where_the_peptide_holds_a_stop_over_another_codon() {
    let (bases, peptide) = coding(&[10], &[9]);
    let mut residues = peptide.into_bytes();
    residues[2] = b'*';
    let tx = pseudogene(
        &format!("{bases}{CLOSING_STOP}"),
        &String::from_utf8(residues).expect("ASCII residues"),
        Some(&UTR.concat()),
    );
    let row = deletion(&tx, 9, 10);
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp9_Ter10delextTer?")
    );
}
