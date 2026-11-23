// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of a frameshift whose alternate translation ends
//! before a residue differs from the reference protein: the documented
//! difference `hgvsp_frameshift_with_no_changed_residue`.
//!
//! Each test lays a made-up coding sequence on a transcript of one exon that
//! ends on its last coding base, with no 3' UTR, and holds the row's HGVSc
//! and HGVSp to values derived from those bases. Ensembl VEP walks the two
//! translations from the row's residue for the first one that differs
//! (TranscriptVariationAllele.pm:2348). Where the alternate translation ends
//! first, the walk leaves by its bound with the pair it read last, two equal
//! residues (2351-2352), and the position raised once more (2361), and
//! returns them as the first changed residue (2364): a residue named at the
//! number after it, changed to itself.

use crate::hgvs::{generate_hgvs, HgvsNotation, HgvsOptions};
use crate::test_helpers::{make_single_exon_transcript, SINGLE_EXON_UTR5};
use vep_core::coordinate::Strand;
use vep_core::transcript::Transcript;
use vep_core::variant::InputVariant;

/// The chromosome of the single-exon test transcript.
const CHROMOSOME: &str = "21";

/// The genomic position of the transcript's first base.
const MODEL_START: u64 = 40_000_001;

/// The genomic position of coding base 1, after the 5' UTR.
const FIRST_CODING: u64 = MODEL_START + SINGLE_EXON_UTR5;

/// Codons 1 to 6: a start codon and five alanine codons, 18 bases.
fn opening() -> String {
    format!("ATG{}", "GCT".repeat(5))
}

/// A coding sequence that no stop codon closes: the six opening codons, a
/// proline codon, an aspartate codon and one base of a ninth codon. Eight
/// whole codons, 25 bases.
fn open_coding_sequence() -> String {
    format!("{}CCGGACT", opening())
}

/// The protein of [`open_coding_sequence`] and of [`closed_coding_sequence`]:
/// eight residues.
const PEPTIDE: &str = "MAAAAAPD";

/// The same eight codons closed by a stop codon, 27 bases.
fn closed_coding_sequence() -> String {
    format!("{}CCGGACTAA", opening())
}

/// A forward-strand transcript of one exon that ends on the last base of
/// `cds`: no 3' UTR, and `peptide` as its protein.
fn model(cds: &str, peptide: &str) -> Transcript {
    let mut tx = make_single_exon_transcript(Strand::Forward, cds, "");
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.three_prime_utr = None;
    vefc.peptide = Some(peptide.to_string());
    tx
}

/// The transcript of [`open_coding_sequence`], flagged as its kind is.
fn open() -> Transcript {
    let mut tx = model(&open_coding_sequence(), PEPTIDE);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    tx
}

/// The coding sequence of `tx`.
fn coding(tx: &Transcript) -> String {
    tx.vefc
        .as_ref()
        .and_then(|vefc| vefc.translateable_seq.clone())
        .expect("a coding sequence")
}

/// A reference that holds the coding sequence of `tx` at its place and one
/// filler base everywhere else, from 1,000 bases before the transcript to
/// 1,000 after it. Nothing is written before the bases held, so the file
/// takes the room of the bases it holds.
fn reference_of(tx: &Transcript) -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
    use std::io::{Seek, SeekFrom, Write};

    let (start, end) = (tx.start - 1000, tx.end + 1000);
    let mut held = vec![b'T'; (end - start + 1) as usize];
    let coding = coding(tx);
    let at = (FIRST_CODING - start) as usize;
    held[at..at + coding.len()].copy_from_slice(coding.as_bytes());
    held.push(b'\n');
    let dir = tempfile::tempdir().expect("tempdir");
    let fa = dir.path().join("ref.fa");
    let header = format!(">{CHROMOSOME}\n");
    let mut file = std::fs::File::create(&fa).expect("create FASTA");
    file.write_all(header.as_bytes())
        .expect("write FASTA header");
    file.seek(SeekFrom::Start(header.len() as u64 + start - 1))
        .expect("seek to the first base held");
    file.write_all(&held).expect("write FASTA bases");
    drop(file);
    std::fs::write(
        dir.path().join("ref.fa.fai"),
        format!(
            "{CHROMOSOME}\t{end}\t{}\t{end}\t{}\n",
            header.len(),
            end + 1
        ),
    )
    .expect("write FASTA index");
    let fasta = vep_fasta::IndexedFasta::from_path_mmap(&fa).expect("load indexed FASTA");
    (dir, fasta)
}

/// The row of a deletion of the coding bases `first` to `last` of `tx`.
fn deletion(tx: &Transcript, first: usize, last: usize) -> HgvsNotation {
    let (_dir, fasta) = reference_of(tx);
    let deleted = coding(tx).as_bytes()[first - 1..last].to_vec();
    let variant = InputVariant::new(
        CHROMOSOME.into(),
        FIRST_CODING + first as u64 - 1,
        FIRST_CODING + last as u64 - 1,
        deleted,
        b"-".to_vec(),
    );
    generate_hgvs(&variant, tx, Some(&fasta), HgvsOptions::default())
}

/// The row of `inserted` written between the coding bases `after` and
/// `after + 1` of `tx`.
fn insertion(tx: &Transcript, after: usize, inserted: &[u8]) -> HgvsNotation {
    let (_dir, fasta) = reference_of(tx);
    let variant = InputVariant::new(
        CHROMOSOME.into(),
        FIRST_CODING + after as u64,
        FIRST_CODING + after as u64 - 1,
        b"-".to_vec(),
        inserted.to_vec(),
    );
    generate_hgvs(&variant, tx, Some(&fasta), HgvsOptions::default())
}

/// Coding base 24 is the third base of codon 8, the last whole codon. Without
/// it the codon closes on the base after it and still codes aspartate, so the
/// eight whole codons the changed sequence holds are the reference's eight
/// residues and no base is left over. Ensembl VEP's walk reads the pair at
/// residue 8, finds it equal, and leaves by its bound at 9: it writes
/// `p.Asp9AspfsTer?`, on a protein of eight residues that no stop closes. No
/// residue the transcript's bases give changes, and no notation is written.
#[test]
fn a_deletion_that_keeps_every_residue_of_an_open_coding_sequence_has_no_protein_notation() {
    let row = deletion(&open(), 24, 24);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.24del"));
    assert_eq!(row.offset, None);
    assert_eq!(row.hgvsp, None);
}

/// One base written between coding bases 23 and 24 closes codon 8 as an
/// aspartate codon, and the last two bases of the coding sequence follow it.
/// They spell a leucine codon whatever base comes third, at residue 9, where
/// the reference protein holds no residue to compare it with. Ensembl VEP
/// writes `p.Asp9AspfsTer?`; no notation is written.
#[test]
fn an_insertion_that_keeps_every_residue_of_an_open_coding_sequence_has_no_protein_notation() {
    let row = insertion(&open(), 23, b"T");
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.23_24insT"));
    assert_eq!(row.offset, None);
    assert_eq!(row.hgvsp, None);
}

/// Without coding bases 21 to 24, the third base of codon 7 and the whole of
/// codon 8, codon 7 closes on the last base of the coding sequence and still
/// codes proline: seven whole codons, the reference's first seven residues,
/// and no base left over. Ensembl VEP's walk leaves by its bound at 8 and
/// writes `p.Pro8ProfsTer?`, naming the reference's aspartate by the proline
/// before it. The changed sequence holds no base of residue 8, so whether it
/// changes is not known, and no notation is written.
#[test]
fn a_deletion_that_leaves_no_base_of_the_next_residue_has_no_protein_notation() {
    let row = deletion(&open(), 21, 24);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.21_24del"));
    assert_eq!(row.offset, None);
    assert_eq!(row.hgvsp, None);
}

/// Without coding bases 21 and 22, the third base of codon 7 and the first of
/// codon 8, codon 7 closes on the second base of codon 8 and still codes
/// proline, and the last two bases of the coding sequence follow it. They
/// spell a leucine codon whatever base comes third, where the reference
/// protein holds an aspartate. A protein is read over whole codons alone, as
/// the reference protein of a coding sequence that ends on such a pair is
/// (`Bio/Tools/CodonTable.pm:448` reads the pair only when complete codons
/// are asked for), so no residue is read from them. Ensembl VEP's walk leaves
/// by its bound at 8 and writes `p.Pro8ProfsTer?`; no notation is written.
#[test]
fn two_bases_after_the_last_whole_codon_that_spell_another_residue_have_no_protein_notation() {
    let row = deletion(&open(), 21, 22);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.21_22del"));
    assert_eq!(row.offset, None);
    assert_eq!(row.hgvsp, None);
}

/// A coding sequence of eight whole codons that no stop codon closes: the
/// start codon, four alanine codons, a proline, a leucine and a threonine
/// codon. Without coding bases 18 to 21, the third base of codon 6 and the
/// whole of codon 7, codon 6 closes on the first base of codon 8 and still
/// codes proline, and the last two bases follow it. They spell a leucine
/// codon whatever base comes third, and the reference protein holds a leucine
/// at residue 7. Ensembl VEP writes `p.Pro7ProfsTer?`; no notation is written.
#[test]
fn two_bases_after_the_last_whole_codon_that_spell_the_same_residue_have_no_protein_notation() {
    let cds = format!("ATG{}CCGCTAACT", "GCT".repeat(4));
    let mut tx = model(&cds, "MAAAAPLT");
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let row = deletion(&tx, 18, 21);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.18_21del"));
    assert_eq!(row.offset, None);
    assert_eq!(row.hgvsp, None);
}

/// A coding sequence that its stop codon closes, with no 3' UTR in the model.
/// Without coding base 24 codon 8 closes on the stop codon's first base and
/// still codes aspartate; the stop codon's last two bases follow it and spell
/// no one residue. Every residue is kept and the stop is broken, with no
/// codon in the model to name what stands in its place. Ensembl VEP's walk
/// leaves by its bound at 9 and writes `p.Asp9AspfsTer?`, naming the stop by
/// the residue before it; no notation is written.
#[test]
fn a_deletion_that_keeps_every_residue_and_breaks_the_stop_codon_has_no_protein_notation() {
    let tx = model(&closed_coding_sequence(), PEPTIDE);
    let row = deletion(&tx, 24, 24);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.24del"));
    assert_eq!(row.offset, None);
    assert_eq!(row.hgvsp, None);
}

/// Coding base 21 is the third base of codon 7 and equals the base after it,
/// so the 3' rule names the deletion at base 22. Codon 7 still codes proline,
/// and codon 8 reads the last three bases of the coding sequence, a threonine
/// codon where the reference holds an aspartate. The walk finds that residue
/// and leaves by it: both engines write the frameshift at residue 8.
#[test]
fn a_frameshift_whose_walk_finds_a_changed_residue_keeps_its_notation() {
    let row = deletion(&open(), 21, 21);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.22del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp8ThrfsTer?")
    );
}
