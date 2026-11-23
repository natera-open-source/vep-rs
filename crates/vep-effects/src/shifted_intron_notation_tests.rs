// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of an indel written in a coding intron and shifted into
//! the coding sequence, on five shapes: a change of the stop codon that keeps
//! every residue before it, a stop codon that stays a stop, a shift longer
//! than the coding coordinate it ends at, a coding sequence that lacks its
//! start and ends on its stop codon, and a coding sequence that does not end
//! on its stop codon.
//!
//! Each test lays a made-up coding sequence and 3' UTR on the three-exon test
//! transcript, on either strand, and holds the row's HGVSc, its offset and its
//! HGVSp to values derived from those bases. Ensembl VEP writes no HGVSp on any
//! of these rows: `hgvs_protein` tests the `coding` predicate of the allele as
//! written, in the intron (TranscriptVariationAllele.pm:1731). It writes the
//! row's HGVSc and offset, except for the insertion of 700 bases, which it
//! describes where its bounded walk stops (`perform_shift`, 303-306).

use crate::concordance_tests::{
    bases_from, compact_on_chr1, hgvsc_fasta_fixture, make_compact_test_transcript,
    make_descending_reverse_strand_transcript,
};
use crate::hgvs::{generate_hgvs, HgvsNotation, HgvsOptions};
use vep_core::transcript::{ExonCoordMapper, Transcript};
use vep_core::variant::InputVariant;

/// The genomic base a cDNA position maps to through the transcript's mapper pairs.
fn genomic_of(tx: &Transcript, cdna: u64) -> u64 {
    let pair = tx
        .vefc
        .as_ref()
        .and_then(|vefc| vefc.mapper.as_ref())
        .expect("a transcript with a mapper")
        .exon_coord_mapper
        .pairs
        .iter()
        .find(|pair| pair.from_start <= cdna && cdna <= pair.from_end)
        .expect("a cDNA position inside an exon");
    let along = cdna - pair.from_start;
    if pair.ori == 1 {
        pair.to_start + along
    } else {
        pair.to_end - along
    }
}

/// `tx` with the coding sequence `cds` from cDNA `first` and the 3' UTR `utr`
/// after it: the mapper, the transcript's coding fields and both sequences.
fn coding(mut tx: Transcript, first: u64, cds: &str, utr: &str) -> Transcript {
    let last = first + cds.len() as u64 - 1;
    let (genomic_first, genomic_last) = (genomic_of(&tx, first), genomic_of(&tx, last));
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.translateable_seq = Some(cds.to_string());
    vefc.three_prime_utr = Some(utr.to_string());
    let mapper = vefc.mapper.as_mut().expect("a mapper");
    mapper.cdna_coding_start = first;
    mapper.cdna_coding_end = last;
    tx.cdna_coding_start = Some(first);
    tx.cdna_coding_end = Some(last);
    tx.coding_region_start = Some(genomic_first.min(genomic_last));
    tx.coding_region_end = Some(genomic_first.max(genomic_last));
    tx.translation_start = Some(genomic_first);
    tx.translation_end = Some(genomic_last);
    tx
}

/// A start codon, 82 alanine codons and `tail`: residues 1 to 83 and what
/// follows them.
fn alanines_then(tail: &str) -> String {
    format!("ATG{}{tail}", "GCT".repeat(82))
}

/// The row of `variant` on `tx` beside the row of `at_the_shift`, the same
/// change given where the 3' shift carries it.
fn rows(
    tx: &Transcript,
    fasta: &vep_fasta::IndexedFasta,
    variant: &InputVariant,
    at_the_shift: &InputVariant,
) -> (HgvsNotation, HgvsNotation) {
    (
        generate_hgvs(variant, tx, Some(fasta), HgvsOptions::default()),
        generate_hgvs(at_the_shift, tx, Some(fasta), HgvsOptions::default()),
    )
}

fn deletion(position: u64, base: &[u8]) -> InputVariant {
    InputVariant::new("1".into(), position, position, base.to_vec(), b"-".to_vec())
}

/// An insertion of `base` between `position - 1` and `position`.
fn insertion(position: u64, base: &[u8]) -> InputVariant {
    InputVariant::new(
        "1".into(),
        position,
        position - 1,
        b"-".to_vec(),
        base.to_vec(),
    )
}

/// The forward transcript's stop codon `TAA` opens exon 2 at 2100..2102, coding
/// bases 250..252, behind a `T` that ends intron 1; the 3' UTR reads `GCATTAG`.
/// Deleting that last intron base shifts one base onto the stop codon's first
/// base. The reading frame runs on as `AAG CAT TAG`: every residue before the
/// stop is kept, the stop's residue becomes a lysine, a histidine follows and
/// the next codon is a stop, an extension of two residues named at the stop.
/// A transcript flagged `cds_start_NF` carries the same notation, since the
/// change leaves every residue before the stop as it is.
#[test]
fn a_deletion_shifted_from_an_intron_onto_the_stop_codon_is_an_extension_named_at_the_stop() {
    let unflagged = coding(
        make_compact_test_transcript(),
        52,
        &alanines_then("TAA"),
        "GCATTAGTT",
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "TTAAGCATTAG"));
    for flags in [vec![], vec!["cds_start_NF"]] {
        let mut tx = unflagged.clone();
        tx.flags = flags.iter().map(ToString::to_string).collect();
        let (row, alone) = rows(&tx, &fasta, &deletion(2099, b"T"), &deletion(2100, b"T"));
        assert_eq!(
            row.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.250del"),
            "{flags:?}"
        );
        assert_eq!(row.offset, Some(1), "{flags:?}");
        assert_eq!(
            row.hgvsp.as_deref(),
            Some("ENSP00000000001.1:p.Ter84LysextTer2"),
            "{flags:?}"
        );
        assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
    }
}

/// The same on the reverse strand, where the shift runs towards lower
/// coordinates: the stop codon opens exon 2 at 4100..4098 and intron 1 ends at
/// 4101, so the forward strand reads, upwards from 4091, the reverse complement
/// of the bases the forward fixture lays down from 2099. Deleting the last
/// intron base shifts one base onto the stop codon, offset -1, and the frame
/// runs on as on the forward strand.
#[test]
fn a_deletion_shifted_onto_the_stop_codon_on_the_reverse_strand_is_the_same_extension() {
    let tx = coding(
        compact_on_chr1(make_descending_reverse_strand_transcript()),
        52,
        &alanines_then("TAA"),
        "GCATTAGTT",
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(4091, "CTAATGCTTAA"));
    let (row, alone) = rows(&tx, &fasta, &deletion(4101, b"A"), &deletion(4100, b"A"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(-1));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter84LysextTer2")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// Intron 1 splits the stop codon `TAG`: its `T` is the last base of exon 1,
/// coding base 250, and `AG` opens exon 2 at 2100..2101 behind an `A` that ends
/// the intron; the 3' UTR opens with an `A`. Deleting the last intron base
/// shifts one base onto coding base 251, and the codon reads `TGA` with the
/// UTR's first base: a stop still stands at residue 84 and the protein is kept
/// whole, so the notation names the stop as unchanged.
#[test]
fn a_deletion_shifted_into_a_stop_codon_that_stays_a_stop_names_the_stop_unchanged() {
    let tx = coding(
        make_compact_test_transcript(),
        51,
        &alanines_then("TAG"),
        "ACCC",
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "AAGACCC"));
    let (row, alone) = rows(&tx, &fasta, &deletion(2099, b"A"), &deletion(2100, b"A"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.251del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(row.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter84="));
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// A shift longer than the coding coordinate it ends at. Intron 1 of the
/// forward transcript is a run of `T` up to 2099 and exon 2 opens with two
/// more, coding bases 251 and 252, the end of the valine codon `GTT`; `GCT AAG`
/// follow. A `T` inserted between 1840 and 1841 shifts 261 bases, past the
/// junction, to stand after coding base 252. Residue 84 is kept, residue 85
/// reads `TGC` and the next codon `TAA`: a frameshift at the alanine of residue
/// 85 with a new stop two residues on.
#[test]
fn an_insertion_shifted_further_than_its_coding_coordinate_keeps_the_frameshift_notation() {
    let tx = coding(
        make_compact_test_transcript(),
        51,
        &alanines_then("GTTGCTAAGGCTGCTTAA"),
        "CC",
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2102, "GCTAAG"));
    let (row, alone) = rows(&tx, &fasta, &insertion(1841, b"T"), &insertion(2102, b"T"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.252dup"));
    assert_eq!(row.offset, Some(261));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala85CysfsTer2")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// The same on the reverse strand: intron 1 reads a run of `A` on the
/// transcript down to 4101 and exon 2 opens with two more at 4100 and 4099, the
/// end of the glutamate codon `GAA`; `GCT AAG` follow, the complement of
/// `CTTAGC` upwards from 4093. A `T` of the forward strand inserted between 4360
/// and 4361 shifts 262 bases to stand after coding base 252. Residue 85 reads
/// `AGC` and the next codon `TAA`.
#[test]
fn an_insertion_shifted_further_than_its_coding_coordinate_on_the_reverse_strand() {
    let tx = coding(
        compact_on_chr1(make_descending_reverse_strand_transcript()),
        51,
        &alanines_then("GAAGCTAAGGCTGCTTAA"),
        "CC",
    );
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(4093, "CTTAGC"));
    let (row, alone) = rows(&tx, &fasta, &insertion(4361, b"T"), &insertion(4099, b"T"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.252dup"));
    assert_eq!(row.offset, Some(-262));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala85SerfsTer2")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// The compact forward transcript flagged `cds_start_NF`, its coding sequence
/// read from cDNA 52: a start codon, 82 alanine codons, then `GAC AAG TAA` from
/// the first base of exon 2 at 2100, coding bases 250 to 258, and the 3' UTR
/// `GCATTAGTT`. The sequence is 86 whole codons, the last of them the stop.
fn forward_transcript_that_lacks_its_coding_start() -> Transcript {
    let mut tx = coding(
        make_compact_test_transcript(),
        52,
        &alanines_then("GACAAGTAA"),
        "GCATTAGTT",
    );
    tx.flags = vec!["cds_start_NF".to_string()].into();
    tx
}

/// A coding sequence that lacks its start and ends on its stop codon. With a
/// `G` ending intron 1 at 2099 before the `G` of the aspartate codon, deleting
/// that last intron base shifts one base onto coding base 250. The frame runs on
/// as `ACA AGT AAG CAT TAG` through the 3' UTR: a frameshift at the aspartate of
/// residue 84, a threonine in its place and the new stop five residues on.
/// Ensembl VEP writes that notation for the deletion given at 2100 and none for
/// the one written in the intron, on a flagged transcript as on any other: no
/// test before TranscriptVariationAllele.pm:1731 reads a transcript flag.
#[test]
fn a_deletion_shifted_out_of_an_intron_keeps_its_frameshift_where_the_coding_start_is_not_found() {
    let tx = forward_transcript_that_lacks_its_coding_start();
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "GGACAAGTAAGCATTAG"));
    let (row, alone) = rows(&tx, &fasta, &deletion(2099, b"G"), &deletion(2100, b"G"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp84ThrfsTer5")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// The same transcript flagged `cds_end_NF` as well carries the same notation:
/// the first residue that differs and the one in its place are whole codons of
/// the sequence, and the new stop is counted on the bases after it.
#[test]
fn a_deletion_shifted_out_of_an_intron_keeps_its_frameshift_where_the_transcript_carries_both_flags(
) {
    let mut tx = forward_transcript_that_lacks_its_coding_start();
    tx.flags = vec!["cds_start_NF".to_string(), "cds_end_NF".to_string()].into();
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "GGACAAGTAAGCATTAG"));
    let (row, alone) = rows(&tx, &fasta, &deletion(2099, b"G"), &deletion(2100, b"G"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp84ThrfsTer5")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// The same on the reverse strand: intron 1 ends at 4101 and exon 2 opens at
/// 4100 with the aspartate codon, so the forward strand reads, upwards from
/// 4085, the reverse complement of the bases the forward fixture lays down from
/// 2099. Deleting the last intron base shifts one base onto coding base 250,
/// offset -1, and the frame runs on as on the forward strand.
#[test]
fn a_deletion_shifted_out_of_an_intron_keeps_its_frameshift_on_the_reverse_strand_too() {
    let mut tx = coding(
        compact_on_chr1(make_descending_reverse_strand_transcript()),
        52,
        &alanines_then("GACAAGTAA"),
        "GCATTAGTT",
    );
    tx.flags = vec!["cds_start_NF".to_string()].into();
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(4085, "CTAATGCTTACTTGTCC"));
    let (row, alone) = rows(&tx, &fasta, &deletion(4101, b"C"), &deletion(4100, b"C"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(-1));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp84ThrfsTer5")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// A coding sequence that opens inside a codon. The start exon's phase is 1, so
/// the cache pads the sequence with one unknown base (`Transcript.pm`
/// `translateable_seq`): codon 1 is that base and coding bases 1 and 2, codons
/// 2 to 83 are alanines, and the glutamate codon `GAG` of residue 84 is coding
/// bases 249 to 251, its `A` the first base of exon 2 at 2100 behind an `A` that
/// ends intron 1. `CAT AAG TAA` close the sequence, 87 whole codons with the
/// pad. Deleting the last intron base shifts one base onto coding base 250, and
/// the frame runs on as `GGC ATA AGT AAG CAT TAG`: a frameshift at the glutamate
/// of residue 84, a glycine in its place and the new stop six residues on. The
/// residue is numbered on the padded codons, as for an allele given in an exon
/// (TranscriptMapper.pm `genomic2pep` adds the start phase).
#[test]
fn a_deletion_shifted_out_of_an_intron_is_numbered_on_the_codons_the_start_phase_pads() {
    let known = format!("AT{}GAGCATAAGTAA", "GCT".repeat(82));
    let mut tx = coding(make_compact_test_transcript(), 52, &known, "GCATTAGTT");
    tx.flags = vec!["cds_start_NF".to_string()].into();
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.translateable_seq = Some(format!("N{known}"));
    vefc.mapper.as_mut().expect("a mapper").start_phase = 1;
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "AAGCATAAGTAAGCATTAG"));
    let (row, alone) = rows(&tx, &fasta, &deletion(2099, b"A"), &deletion(2100, b"A"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Glu84GlyfsTer6")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// An in-frame insertion on a coding sequence that lacks its start. Intron 1 of
/// the forward transcript is a run of `T` up to 2099 and exon 2 opens with two
/// more, coding bases 251 and 252, the end of the valine codon `GTT`; `GCT AAG
/// TAA` follow, 87 whole codons in all. Three `T` inserted between 2098 and
/// 2099 shift three bases, past the junction, to stand after coding base 252:
/// a duplication of the intron's last base and the exon's first two. The frame
/// is kept, and a phenylalanine codon stands between the valine of residue 84
/// and the alanine of residue 85.
#[test]
fn an_in_frame_insertion_shifted_out_of_an_intron_keeps_its_notation_where_the_coding_start_is_not_found(
) {
    let mut tx = coding(
        make_compact_test_transcript(),
        51,
        &alanines_then("GTTGCTAAGTAA"),
        "CC",
    );
    tx.flags = vec!["cds_start_NF".to_string()].into();
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2102, "GCTAAGTAACC"));
    let (row, alone) = rows(
        &tx,
        &fasta,
        &insertion(2099, b"TTT"),
        &insertion(2102, b"TTT"),
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.251-1_252dup")
    );
    assert_eq!(row.offset, Some(3));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Val84_Ala85insPhe")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// The compact forward transcript flagged `cds_start_NF` with a coding sequence
/// of 402 bases read from cDNA 51: a start codon, 82 alanine codons, a valine
/// codon across the first junction, 49 alanine codons and the stop codon, which
/// ends at coding base 402, inside exon 2.
fn transcript_whose_coding_sequence_ends_inside_exon_2() -> Transcript {
    let cds = alanines_then(&format!("GTT{}TAA", "GCT".repeat(49)));
    assert_eq!(cds.len(), 402);
    let mut tx = coding(make_compact_test_transcript(), 51, &cds, "CC");
    tx.flags = vec!["cds_start_NF".to_string()].into();
    tx
}

/// An insertion whose run of placements is longer than Ensembl VEP's shift
/// walks and ends in an intron, on a coding sequence that lacks its start. The
/// reference is a run of `T` through intron 1, exon 2 and the first 51 bases of
/// intron 2, up to 2450. Seven hundred `T` inserted between 1899 and 1900 stand
/// at the end of their run 551 bases on, in intron 2, where the HGVSc describes
/// them as a duplication from intron 1 to intron 2. Both notations are read at
/// that placement, which has no translation coordinates, so the row has no
/// HGVSp.
#[test]
fn an_insertion_whose_run_ends_in_an_intron_has_no_notation_where_the_coding_start_is_not_found() {
    let tx = transcript_whose_coding_sequence_ends_inside_exon_2();
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2451, "G"));
    let row = generate_hgvs(
        &insertion(1900, &[b'T'; 700]),
        &tx,
        Some(&fasta),
        HgvsOptions::default(),
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.251-349_*148+51dup")
    );
    assert_eq!(row.offset, Some(551));
    assert_eq!(row.hgvsp, None);
}

/// The same insertion where the run ends inside exon 2: the reference is a run
/// of `T` up to 2230, coding base 381. The 700 `T` stand at the end of their
/// run 331 bases on, after coding base 381, as a duplication that opens in the
/// first intron. Ensembl VEP's walk stops 301 bases on and its row has no
/// HGVSp. The notation is read at the end of the run: one base more than 233
/// codons is inserted before the alanine codon of residue 128, which reads as
/// a phenylalanine, and the shifted frame meets no stop codon in the coding
/// sequence or the two bases of the 3' UTR.
#[test]
fn an_insertion_whose_run_ends_in_the_coding_sequence_keeps_its_notation_where_the_coding_start_is_not_found(
) {
    let tx = transcript_whose_coding_sequence_ends_inside_exon_2();
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2231, "G"));
    let row = generate_hgvs(
        &insertion(1900, &[b'T'; 700]),
        &tx,
        Some(&fasta),
        HgvsOptions::default(),
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.251-569_381dup")
    );
    assert_eq!(row.offset, Some(331));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala128PhefsTer?")
    );
}

/// `tx` cut off after its coding sequence, as a transcript whose coding end is
/// not found: the exon that holds the last coding base ends on it and no exon
/// follows.
fn cut_off_after_its_coding_sequence(mut tx: Transcript) -> Transcript {
    let last = tx.cdna_coding_end.expect("a coding sequence");
    let end = genomic_of(&tx, last);
    tx.end = end;
    tx.exons.retain(|exon| exon.start <= end);
    tx.introns.retain(|intron| intron.end < end);
    for exon in &mut tx.exons {
        exon.end = exon.end.min(end);
    }
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.three_prime_utr = None;
    vefc.sorted_exons.retain(|exon| exon.start <= end);
    vefc.introns.retain(|intron| intron.end < end);
    for exon in &mut vefc.sorted_exons {
        exon.end = exon.end.min(end);
    }
    let mapper = &mut vefc.mapper.as_mut().expect("a mapper").exon_coord_mapper;
    mapper.pairs.retain(|pair| pair.from_start <= last);
    for pair in &mut mapper.pairs {
        pair.from_end = pair.from_end.min(last);
        pair.to_end = pair.to_end.min(end);
    }
    *mapper = ExonCoordMapper::new(std::mem::take(&mut mapper.pairs));
    tx
}

/// The forward transcript flagged `cds_end_NF`: 83 residues in exon 1 and
/// `tail` from the first base of exon 2 at 2100, coding base 250, to the
/// transcript's last base. No stop codon closes the coding sequence and no base
/// follows it.
fn forward_transcript_whose_coding_end_is_not_found(tail: &str) -> Transcript {
    let mut tx = cut_off_after_its_coding_sequence(coding(
        make_compact_test_transcript(),
        52,
        &alanines_then(tail),
        "",
    ));
    tx.flags = vec!["cds_end_NF".to_string()].into();
    tx
}

/// A frameshift whose new stop codon lies inside the coding sequence. Exon 2
/// reads `GAT AAG CAT TGG CTA GC`, residues 84 to 88 and two bases of a codon
/// the sequence does not complete. With a `G` ending intron 1 at 2099, deleting
/// that base shifts one base onto coding base 250, and the frame runs on as
/// `ATA AGC ATT GGC TAG`: a frameshift at the aspartate of residue 84, an
/// isoleucine in its place and the new stop five residues on, every base of it
/// a coding base the transcript holds. Ensembl VEP writes that notation for the
/// deletion given at 2100 and none for the one written in the intron: no test
/// at or before TranscriptVariationAllele.pm:1731 reads a transcript flag.
#[test]
fn a_frameshift_shifted_out_of_an_intron_keeps_its_notation_where_the_coding_end_is_not_found() {
    let tx = forward_transcript_whose_coding_end_is_not_found("GATAAGCATTGGCTAGC");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "GGATAAGCATTGGCTAGC"));
    let (row, alone) = rows(&tx, &fasta, &deletion(2099, b"G"), &deletion(2100, b"G"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp84IlefsTer5")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// A residue replaced by a stop. Exon 2 reads `GTA GCA TTG GCT AGC`, residues
/// 84 to 88. Deleting the `G` that ends intron 1 shifts one base onto coding
/// base 250, and the valine codon of residue 84 becomes `TAG`.
#[test]
fn a_stop_shifted_out_of_an_intron_keeps_its_notation_where_the_coding_end_is_not_found() {
    let tx = forward_transcript_whose_coding_end_is_not_found("GTAGCATTGGCTAGC");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "GGTAGCATTGGCTAGC"));
    let (row, alone) = rows(&tx, &fasta, &deletion(2099, b"G"), &deletion(2100, b"G"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(row.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Val84Ter"));
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// A frameshift whose frame meets no stop codon before the transcript's last
/// base. Exon 2 reads `GAT GCT GCT`, residues 84 to 86. Deleting the `G` that
/// ends intron 1 shifts one base onto coding base 250, and the frame runs on as
/// `ATG CTG CT`: a methionine in place of the aspartate of residue 84, and the
/// count to the new stop is the unknown one.
#[test]
fn a_frameshift_that_meets_no_stop_codon_is_counted_as_unknown_where_the_coding_end_is_not_found() {
    let tx = forward_transcript_whose_coding_end_is_not_found("GATGCTGCT");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "GGATGCTGCT"));
    let (row, alone) = rows(&tx, &fasta, &deletion(2099, b"G"), &deletion(2100, b"G"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp84MetfsTer?")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// A deleted residue that a residue of the coding sequence follows.
/// Intron 1 ends `GAC` at 2097 to 2099 and exon 2 reads `GAC AAG CAT TGG CT`,
/// residues 84 to 87 and two bases more. Deleting the three intron bases shifts
/// three bases, onto coding bases 250 to 252: the aspartate of residue 84 is
/// deleted, and the lysine after it stands unchanged.
#[test]
fn a_deleted_residue_shifted_out_of_an_intron_keeps_its_notation_where_the_coding_end_is_not_found()
{
    let tx = forward_transcript_whose_coding_end_is_not_found("GACAAGCATTGGCT");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2097, "GACGACAAGCATTGGCT"));
    let written = InputVariant::new("1".into(), 2097, 2099, b"GAC".to_vec(), b"-".to_vec());
    let shifted = InputVariant::new("1".into(), 2100, 2102, b"GAC".to_vec(), b"-".to_vec());
    let (row, alone) = rows(&tx, &fasta, &written, &shifted);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250_252del"));
    assert_eq!(row.offset, Some(3));
    assert_eq!(row.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Asp84del"));
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// A change of whole codons that replaces the sequence's last whole codon. The
/// coding sequence opens one base earlier, so codon 84 is the last base of exon
/// 1 and the first two of exon 2: `G AT`, an aspartate, and `CGT`, an arginine,
/// closes the sequence as its last whole codon. Intron 1 ends `ATC` at 2097 to
/// 2099, the bases exon 2 opens with, so deleting them shifts three bases, onto
/// coding bases 251 to 253. Codon 84 is then `G GT`: a glycine stands where the
/// aspartate and the arginine stood. No residue of the sequence follows, and
/// none is needed: HGVS moves no deletion-insertion, so the bases the
/// transcript holds give the description whole.
#[test]
fn residues_replaced_in_the_last_whole_codon_keep_their_notation_where_the_coding_end_is_not_found()
{
    let mut tx = cut_off_after_its_coding_sequence(coding(
        make_compact_test_transcript(),
        51,
        &alanines_then("GATCGT"),
        "",
    ));
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2097, "ATCATCGT"));
    let written = InputVariant::new("1".into(), 2097, 2099, b"ATC".to_vec(), b"-".to_vec());
    let shifted = InputVariant::new("1".into(), 2100, 2102, b"ATC".to_vec(), b"-".to_vec());
    let (row, alone) = rows(&tx, &fasta, &written, &shifted);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.251_253del"));
    assert_eq!(row.offset, Some(3));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp84_Arg85delinsGly")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// The same deletion where the aspartate is the sequence's last whole codon:
/// exon 2 reads `GAC` and the transcript ends on it. No residue of the coding
/// sequence follows the deleted one, so the residues past the sequence's end
/// decide where the description stands, and the row keeps the empty HGVSp.
#[test]
fn a_deletion_of_the_last_whole_codon_has_no_notation_where_the_coding_end_is_not_found() {
    let tx = forward_transcript_whose_coding_end_is_not_found("GAC");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2097, "GACGAC"));
    let written = InputVariant::new("1".into(), 2097, 2099, b"GAC".to_vec(), b"-".to_vec());
    let row = generate_hgvs(&written, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250_252del"));
    assert_eq!(row.offset, Some(3));
    assert_eq!(row.hgvsp, None);
}

/// A frameshift in the sequence's last whole codon. Exon 2 reads `GAC` and the
/// transcript ends on it. Deleting the `G` that ends intron 1 shifts one base
/// onto coding base 250 and leaves two bases of residue 84's codon: no whole
/// codon gives the residue that replaces the aspartate, and the row keeps the
/// empty HGVSp.
#[test]
fn a_frameshift_in_the_last_whole_codon_has_no_notation_where_the_coding_end_is_not_found() {
    let tx = forward_transcript_whose_coding_end_is_not_found("GAC");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "GGAC"));
    let row = generate_hgvs(
        &deletion(2099, b"G"),
        &tx,
        Some(&fasta),
        HgvsOptions::default(),
    );
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(row.hgvsp, None);
}

/// A deleted residue in a run of the same residue that reaches the sequence's
/// end. Intron 1 ends `GAC` and exon 2 reads `GAC GAC`, the aspartates of
/// residues 84 and 85, on which the transcript ends. Deleting the three intron
/// bases shifts six bases, onto the last whole codon. No residue of the
/// sequence's whole codons differs: an aspartate past its end would carry the
/// description further, and the row keeps the empty HGVSp.
#[test]
fn a_deletion_in_a_run_that_reaches_the_end_has_no_notation_where_the_coding_end_is_not_found() {
    let tx = forward_transcript_whose_coding_end_is_not_found("GACGAC");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2097, "GACGACGAC"));
    let written = InputVariant::new("1".into(), 2097, 2099, b"GAC".to_vec(), b"-".to_vec());
    let row = generate_hgvs(&written, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.253_255del"));
    assert_eq!(row.offset, Some(6));
    assert_eq!(row.hgvsp, None);
}

/// A residue duplicated where the whole codons end. Intron 1 ends `GAC` and
/// exon 2 reads `GAC T`, the aspartate of residue 84 and one base of a codon
/// the sequence does not complete. Three bases `GAC` inserted before those
/// intron bases shift six bases, to stand between coding bases 252 and 253. The
/// second aspartate follows the last whole codon, where no residue of the
/// sequence stands to differ, and the row keeps the empty HGVSp.
#[test]
fn a_duplication_after_the_last_whole_codon_has_no_notation_where_the_coding_end_is_not_found() {
    let tx = forward_transcript_whose_coding_end_is_not_found("GACT");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2097, "GACGACT"));
    let row = generate_hgvs(
        &insertion(2097, b"GAC"),
        &tx,
        Some(&fasta),
        HgvsOptions::default(),
    );
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250_252dup"));
    assert_eq!(row.offset, Some(6));
    assert_eq!(row.hgvsp, None);
}

/// A frameshift whose new frame keeps every residue to the sequence's end.
/// Exon 2 reads `CTT A`, the leucine of residue 84 and one base more. Deleting
/// the `C` that ends intron 1 shifts one base onto coding base 250, and the
/// last whole codon reads `TTA`, a leucine still: no residue of the sequence's
/// whole codons differs, and the row keeps the empty HGVSp.
#[test]
fn a_frameshift_that_keeps_every_residue_has_no_notation_where_the_coding_end_is_not_found() {
    let tx = forward_transcript_whose_coding_end_is_not_found("CTTA");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2099, "CCTTA"));
    let row = generate_hgvs(
        &deletion(2099, b"C"),
        &tx,
        Some(&fasta),
        HgvsOptions::default(),
    );
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.250del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!(row.hgvsp, None);
}
