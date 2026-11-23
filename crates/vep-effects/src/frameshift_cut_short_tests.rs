// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of a frameshift that the end of the transcript cuts
//! short: the `HGVSp` pair of the documented difference
//! `stop_codon_overlap_without_stop_codon`.
//!
//! Each test lays a made-up coding sequence on a transcript of one exon that
//! ends on its last coding base, with no 3' UTR, and holds the row's HGVSc
//! and HGVSp to values derived from those bases. Where a deletion shifts the
//! frame and fewer than three bases are left from the first changed codon on,
//! the alternate translation ends before the first changed residue. Ensembl
//! VEP's frameshift routine then words the row's reference residues as
//! deleted (TranscriptVariationAllele.pm:2342-2346) and, where its lost-stop
//! test holds, on every transcript without `cds_end_NF`
//! (VariationEffect.pm:1373-1381), adds an extension of unknown length
//! (TranscriptVariationAllele.pm:1916-1936).
//!
//! No stop codon stands on the coding sequences here, so none is lost; and
//! the sequence shows the first changed codon cut short, not its residue
//! absent with the rest of the protein standing. The row is a frameshift,
//! and the HGVS nomenclature's short form, a residue and `fs`, marks where
//! the frame breaks. Where Ensembl VEP reads no lost stop it writes the
//! deletion alone, and the row keeps it.

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

/// Eight whole codons, 24 bases: a start codon, five alanine codons, a
/// proline codon at 19-21 and an aspartate codon at 22-24.
fn whole_codons() -> String {
    format!("ATG{}CCGGAC", "GCT".repeat(5))
}

/// The protein of [`whole_codons`]: eight residues.
const PEPTIDE: &str = "MAAAAAPD";

/// A forward-strand transcript of one exon that ends on the last base of
/// `cds`: no 3' UTR, and `peptide` as its protein.
fn model_of(cds: &str, peptide: &str) -> Transcript {
    let mut tx = make_single_exon_transcript(Strand::Forward, cds, "");
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.three_prime_utr = None;
    vefc.peptide = Some(peptide.to_string());
    tx
}

/// [`model_of`] with [`PEPTIDE`] as the protein.
fn model(cds: &str) -> Transcript {
    model_of(cds, PEPTIDE)
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
/// 1,000 after it.
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

/// The protein notation of `row`, after the protein's id.
fn protein_notation(row: &HgvsNotation) -> Option<&str> {
    row.hgvsp
        .as_deref()
        .map(|hgvsp| hgvsp.trim_start_matches("ENSP00000000001.1:"))
}

/// Codon 8, the last, names aspartate, and nothing follows it. Without one
/// of its bases, or without its first two, the changed sequence holds seven
/// whole codons, the reference's first seven residues, and two bases or one
/// of an eighth. Ensembl VEP writes `p.Asp8delextTer?` on each: the residue
/// deleted, and a stop lost where none stands. The frame shifts at residue
/// 8, and the row names that residue as the first one affected.
#[test]
fn a_deletion_inside_the_last_codon_of_an_unclosed_coding_sequence_is_a_frameshift_at_its_residue()
{
    let tx = model(&whole_codons());
    for (first, last, hgvsc) in [
        (22, 22, "c.22del"),
        (23, 23, "c.23del"),
        (24, 24, "c.24del"),
        (22, 23, "c.22_23del"),
    ] {
        let row = deletion(&tx, first, last);
        assert_eq!(
            row.hgvsc.as_deref(),
            Some(format!("ENST00000000001.1:{hgvsc}").as_str())
        );
        assert_eq!(protein_notation(&row), Some("p.Asp8fs"), "{hgvsc}");
    }
}

/// Without coding bases 19 to 22, the whole of codon 7 and the first base of
/// codon 8, the changed sequence holds six whole codons and two bases of a
/// seventh. Ensembl VEP writes `p.Pro7_Asp8delextTer?`, the two residues
/// deleted and a stop lost. The frame shifts at residue 7, the first of the
/// two, and the row names it alone.
#[test]
fn a_deletion_from_an_earlier_codon_into_the_last_is_a_frameshift_at_the_first_of_its_residues() {
    let row = deletion(&model(&whole_codons()), 19, 22);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.19_22del"));
    assert_eq!(protein_notation(&row), Some("p.Pro7fs"));
}

/// A coding sequence that ends one base into a ninth codon, a base that
/// adenines complete to no stop codon. Without the first two bases of codon 8
/// the changed sequence holds seven whole codons and two bases: the same
/// frameshift at residue 8.
#[test]
fn a_frameshift_cut_short_before_one_base_of_a_next_codon_is_named_at_its_first_residue() {
    let tx = model(&format!("{}C", whole_codons()));
    let row = deletion(&tx, 22, 23);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.22_23del"));
    assert_eq!(protein_notation(&row), Some("p.Asp8fs"));
}

/// A coding sequence that ends two bases into a ninth codon, bases that
/// adenines complete to a leucine codon. Without the last two bases of codon
/// 8 and both bases after it, the changed sequence holds seven whole codons
/// and one base: the frameshift is named at residue 8. A deletion that
/// begins in the two bases after the last whole codon has no residue of the
/// protein to name and keeps no notation: one written there, and one written
/// on the last base of codon 8 that the 3' shift names at the equal base
/// after it, where the frameshift walk starts at a ninth residue.
#[test]
fn a_frameshift_cut_short_before_two_bases_of_a_next_codon_is_named_at_a_residue_of_the_protein() {
    let tx = model(&format!("{}CT", whole_codons()));
    let row = deletion(&tx, 23, 26);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.23_26del"));
    assert_eq!(protein_notation(&row), Some("p.Asp8fs"));
    let in_the_tail = deletion(&tx, 25, 25);
    assert_eq!(
        in_the_tail.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.25del")
    );
    assert_eq!(protein_notation(&in_the_tail), None);
    let carried_into_the_tail = deletion(&tx, 24, 24);
    assert_eq!(
        carried_into_the_tail.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.25del")
    );
    assert_eq!(carried_into_the_tail.offset, Some(1));
    assert_eq!(protein_notation(&carried_into_the_tail), None);
}

/// The same eight codons on a transcript that carries `cds_end_NF`. Ensembl
/// VEP's lost-stop test does not hold under the flag (VariationEffect.pm:1373),
/// so it writes the deletion alone, `p.Asp8del`, and the row keeps that form.
#[test]
fn a_frameshift_cut_short_on_a_transcript_flagged_incomplete_at_its_end_keeps_the_deletion() {
    let mut tx = model(&whole_codons());
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let row = deletion(&tx, 22, 22);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.22del"));
    assert_eq!(protein_notation(&row), Some("p.Asp8del"));
}

/// A coding sequence that ends one base into a ninth codon. Without the whole
/// of codon 8 and the base after it, four bases, the changed sequence ends on
/// its seventh whole codon: no base of a codon is left for the frame to break
/// in. Ensembl VEP reads no lost stop where the changed allele holds no base
/// (VariationEffect.pm:1259-1260) and writes the range deleted
/// (TranscriptVariationAllele.pm:2018-2028), and the row keeps that form.
#[test]
fn a_deletion_of_whole_codons_to_the_last_base_is_no_frameshift_cut_short() {
    let tx = model(&format!("{}C", whole_codons()));
    let row = deletion(&tx, 22, 25);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.22_25del"));
    assert_eq!(protein_notation(&row), Some("p.Asp8_Ter9del"));
}

/// Eight whole codons whose seventh is three adenines and whose eighth
/// opens on two more: a run of five adenines over coding bases 19 to 23. A
/// record that deletes the first of them is written before the last three
/// coding positions, and the 3' shift names the deletion at base 23, inside
/// the last codon, where two bases of it are left. Ensembl VEP reads its
/// lost-stop test where the record is written
/// (VariationEffect.pm:1375-1381), finds no overlap and writes the deletion
/// alone; the row keeps that form.
#[test]
fn a_deletion_the_shift_carries_into_the_last_codon_keeps_the_deletion() {
    let cds = format!("ATG{}AAAAAC", "GCT".repeat(5));
    let row = deletion(&model_of(&cds, "MAAAAAKN"), 19, 19);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.23del"));
    assert_eq!(row.offset, Some(4));
    assert_eq!(protein_notation(&row), Some("p.Asn8del"));
}
