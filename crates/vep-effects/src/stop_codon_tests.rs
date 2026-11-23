// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation at the stop codon, on the shapes of recorded rows.
//!
//! Each test lays a transcript model at its own coordinates (its exons, strand
//! and coding bounds) over a reference that holds synthetic bases around the
//! transcript alone, gives the allele with the lengths and the place of the
//! row's record, and asserts the notation the HGVS nomenclature gives for the
//! two proteins: the stop written as unchanged where the protein does not
//! change, an extension named at the stop where the whole reference protein is
//! kept and residues are added after it, with their count where the peptide
//! reads a stop codon inside the coding sequence as a residue, a frameshift
//! named at the first changed residue where the nomenclature gives one for a
//! change that reaches from before the stop codon into it, and no notation
//! where the description would name a stop that the reference protein does
//! not hold. An insertion inside the codon of the last residue, which the
//! codon window leaves without a description, takes the stop written as
//! unchanged where the protein does not change, at any length; one there that
//! adds residues before the stop, and one between that codon and the stop
//! codon, are asserted as written, without a notation.
//!
//! The tests after them assert the stop codon's terms, `stop_lost` and
//! `stop_retained_variant`, which are written only where a stop codon stands
//! at the coding end. Each lays a made-up transcript on both strands and
//! asserts the whole term set of a row at a coding end that a stop codon
//! closes, that polyadenylation completes or that holds none. Ensembl VEP
//! reads the last three coding positions as the stop codon on every
//! transcript without `cds_end_NF`, whatever they hold
//! (VariationEffect.pm:1373-1381).
//!
//! A coding sequence that ends one or two bases into a stop codon on the
//! transcript's last base is read by what its record shows. Where the record
//! shows the model complete, by a mitochondrial codon table or a consensus
//! coding sequence id, every term and cell is that of the coding sequence
//! with the adenines that complete the codon and the adenines of the poly(A)
//! tail after it, which a change that removes bases reads on into. Where it
//! does not, no stop codon stands and no row takes a stop term from the
//! coding end. A structural deletion is read at the bases the genome has of
//! a stop codon the record shows completed.

use crate::consequences::{calculate_consequences, EffectsConfig};
use crate::hgvs::{generate_hgvs, HgvsNotation, HgvsOptions};
use crate::test_helpers::{assert_consequence_set_eq, make_single_exon_transcript};
use std::sync::Arc;
use vep_core::codon::reverse_complement;
use vep_core::consequence::TranscriptConsequence;
use vep_core::coordinate::Strand;
use vep_core::transcript::{Attribute, Exon, ExonCoordMapper, Intron, MapperPair, Transcript};
use vep_core::variant::InputVariant;

/// A transcript model at its own coordinates.
struct RowModel<'a> {
    chr: &'static str,
    strand: Strand,
    /// The exons in transcript order, each by its first and last genomic base.
    exons: &'a [(u64, u64)],
    cdna_coding_start: u64,
    cdna_coding_end: u64,
    /// The mapper's start phase: the bases the first codon lacks where the
    /// coding sequence opens inside it, `-1` or `0` where it opens on a codon.
    start_phase: i8,
}

impl RowModel<'_> {
    /// The transcript's first and last genomic base.
    fn span(&self) -> (u64, u64) {
        let first = self.exons.iter().map(|e| e.0).min().expect("an exon");
        let last = self.exons.iter().map(|e| e.1).max().expect("an exon");
        (first, last)
    }
}

/// A reference for the model's chromosome on which the transcript is read at
/// its own coordinates: the record ends 1,000 bases after the transcript and
/// the file holds its bases from 1,000 before the transcript on alone, one
/// filler base but for `edits`, each a run of forward-strand bases laid from a
/// genomic position. Nothing is written before the bases held, so the file
/// takes the room of the bases it holds.
fn reference_of(
    model: &RowModel<'_>,
    edits: &[(u64, &str)],
) -> (tempfile::TempDir, vep_fasta::IndexedFasta) {
    use std::io::{Seek, SeekFrom, Write};

    let (first, last) = model.span();
    let (start, end) = (first - 1000, last + 1000);
    let mut held = vec![b'T'; (end - start + 1) as usize];
    for &(pos, bases) in edits {
        let at = (pos - start) as usize;
        held[at..at + bases.len()].copy_from_slice(bases.as_bytes());
    }
    held.push(b'\n');
    let dir = tempfile::tempdir().expect("tempdir");
    let fa = dir.path().join("ref.fa");
    let header = format!(">{}\n", model.chr);
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
            "{}\t{end}\t{}\t{end}\t{}\n",
            model.chr,
            header.len(),
            end + 1
        ),
    )
    .expect("write FASTA index");
    let fasta = vep_fasta::IndexedFasta::from_path_mmap(&fa).expect("load indexed FASTA");
    (dir, fasta)
}

/// The model's transcript over `fasta`, with the coding sequence and the 3' UTR
/// a cache carries read from the reference: the coding sequence opens with one
/// `N` per base its first codon lacks.
fn transcript_of(model: &RowModel<'_>, fasta: &vep_fasta::IndexedFasta) -> Transcript {
    let reverse = model.strand == Strand::Reverse;
    let mut cdna: Vec<u8> = Vec::new();
    let mut pairs = Vec::new();
    for &(start, end) in model.exons {
        let bases = fasta
            .sequence(model.chr, start, end)
            .expect("the exon's bases");
        let from_start = cdna.len() as u64 + 1;
        cdna.extend(if reverse {
            reverse_complement(&bases)
        } else {
            bases
        });
        pairs.push(MapperPair {
            from_start,
            from_end: cdna.len() as u64,
            to_start: start,
            to_end: end,
            ori: if reverse { -1 } else { 1 },
        });
    }
    let exon_of = |cdna_pos: u64| {
        pairs
            .iter()
            .position(|p| p.from_start <= cdna_pos && cdna_pos <= p.from_end)
            .expect("a cDNA position on an exon")
    };
    let genomic = |cdna_pos: u64| {
        let p = &pairs[exon_of(cdna_pos)];
        if reverse {
            p.to_end - (cdna_pos - p.from_start)
        } else {
            p.to_start + (cdna_pos - p.from_start)
        }
    };
    let (coding_start, coding_end) = (model.cdna_coding_start, model.cdna_coding_end);
    let coding = &cdna[coding_start as usize - 1..coding_end as usize];
    let mut cds = "N".repeat(model.start_phase.max(0) as usize);
    cds.push_str(std::str::from_utf8(coding).expect("ASCII bases"));
    let utr3 = String::from_utf8(cdna[coding_end as usize..].to_vec()).expect("ASCII bases");

    let mut tx = make_single_exon_transcript(model.strand, &cds, &utr3);
    let (first, last) = model.span();
    tx.chr = model.chr.into();
    tx.start = first;
    tx.end = last;
    tx.exons = model
        .exons
        .iter()
        .enumerate()
        .map(|(i, &(start, end))| Exon {
            start,
            end,
            rank: i as u32 + 1,
            phase: -1,
            end_phase: -1,
        })
        .collect();
    tx.introns = model
        .exons
        .windows(2)
        .enumerate()
        .map(|(i, pair)| {
            let (start, end) = if reverse {
                (pair[1].1 + 1, pair[0].0 - 1)
            } else {
                (pair[0].1 + 1, pair[1].0 - 1)
            };
            Intron {
                start,
                end,
                rank: i as u32 + 1,
            }
        })
        .collect();
    let (at_start, at_end) = (genomic(coding_start), genomic(coding_end));
    tx.cdna_coding_start = Some(coding_start);
    tx.cdna_coding_end = Some(coding_end);
    tx.coding_region_start = Some(at_start.min(at_end));
    tx.coding_region_end = Some(at_start.max(at_end));
    tx.translation_start = Some(at_start);
    tx.translation_end = Some(at_end);
    if model.start_phase > 0 {
        tx.flags = vec!["cds_start_NF".to_string()].into();
    }
    let (start_exon, end_exon) = (exon_of(coding_start), exon_of(coding_end));
    let translation = tx.translation.as_mut().expect("a translation");
    translation.start = coding_start - pairs[start_exon].from_start + 1;
    translation.end = coding_end - pairs[end_exon].from_start + 1;
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.sorted_exons = tx.exons.clone();
    vefc.introns = tx.introns.clone();
    let mapper = vefc.mapper.as_mut().expect("a mapper");
    mapper.start_phase = model.start_phase;
    mapper.cdna_coding_start = coding_start;
    mapper.cdna_coding_end = coding_end;
    mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
    tx
}

/// The shape of the row at GRCh37 1:90470803 on ENST00000361911: three bases
/// inserted after CDS base 209, a copy of the three after them, which the 3'
/// rule carries 3 bases on, to between the second and third base of the stop
/// codon (residue 71). The first inserted base closes a stop codon with the
/// stop codon's first two bases, at the stop's own residue, so the protein is
/// the reference protein and the stop is written as unchanged, where the codon
/// window gives a deletion-insertion that replaces the stop by a stop.
#[test]
fn insertion_in_the_stop_codon_that_leaves_the_protein_unchanged_is_the_stop_unchanged() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Forward,
        exons: &[
            (90_460_693, 90_460_832),
            (90_463_669, 90_463_713),
            (90_470_562, 90_470_597),
            (90_470_692, 90_472_008),
        ],
        cdna_coding_start: 125,
        cdna_coding_end: 337,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[(90_460_817, "ATG"), (90_470_787, "ACATACACGTCGGTAGAGTAA")],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "1".into(),
        90_470_804,
        90_470_803,
        b"-".to_vec(),
        b"GTA".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, Some(3));
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter71="));
}

/// The shape of the row at GRCh37 3:182804570 on ENST00000466650, a
/// reverse-strand transcript: four bases inserted after the first base of the
/// stop codon (CDS base 142, residue 48). The first two close a stop codon with
/// that base, at the stop's own residue, so the protein is the reference
/// protein and the stop is written as unchanged.
#[test]
fn insertion_in_the_stop_codon_that_leaves_the_protein_unchanged_on_the_reverse_strand() {
    let model = RowModel {
        chr: "3",
        strand: Strand::Reverse,
        exons: &[
            (182_817_140, 182_817_318),
            (182_812_347, 182_812_393),
            (182_804_481, 182_804_576),
            (182_788_998, 182_789_145),
            (182_788_836, 182_788_932),
        ],
        cdna_coding_start: 91,
        cdna_coding_end: 234,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (182_804_569, "CTAGCTAG"),
            (182_812_347, "AGCGTCGCGAGTC"),
            (182_817_226, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "3".into(),
        182_804_571,
        182_804_570,
        b"-".to_vec(),
        b"AGTC".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter48="));
}

/// The notation of bases inserted after the first base of the codon of
/// residue 10 on a made-up forward-strand transcript of one exon. Its coding
/// sequence opens at 61,031 and holds eleven codons, those of residues 10
/// (Lys) and 11 (Val) at 61,058 and 61,061. Where `closed`, the stop codon
/// follows them as residue 12; where not, the transcript ends on the codon of
/// residue 11 and no stop codon closes its coding sequence. The first
/// inserted base equals the codon's second base in every test below, so the 3'
/// rule carries the insertion one base on, inside the same codon.
fn insertion_in_the_codon_of_residue_10(closed: bool, inserted: &[u8]) -> HgvsNotation {
    insertion_after_base(closed, 61_058, inserted)
}

/// The notation of `inserted` written after the base `after` of the made-up
/// transcript of `insertion_in_the_codon_of_residue_10`, whose codons of
/// residues 9 (Ile), 10 (Lys) and 11 (Val) open at 61,055, 61,058 and 61,061.
fn insertion_after_base(closed: bool, after: u64, inserted: &[u8]) -> HgvsNotation {
    let model = RowModel {
        chr: "7",
        strand: Strand::Forward,
        exons: if closed {
            &[(61_001, 61_120)]
        } else {
            &[(61_001, 61_063)]
        },
        cdna_coding_start: 31,
        cdna_coding_end: if closed { 66 } else { 63 },
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (61_031, "ATGGCACGTGAC"),
            (61_043, "AACTGCGAGCAT"),
            (61_055, "ATCAAGGTCTAA"),
        ],
    );
    let mut tx = transcript_of(&model, &fasta);
    if !closed {
        tx.flags = vec!["cds_end_NF".to_string()].into();
    }
    let variant = InputVariant::new(
        "7".into(),
        after + 1,
        after,
        b"-".to_vec(),
        inserted.to_vec(),
    );
    generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default())
}

/// Eight bases inserted before the stop codon that leave the protein as it
/// was: they complete the codon of residue 10 to the same residue, then read
/// as a codon of residue 11's residue and a stop codon, which stands at the
/// stop's own residue, 12. The alternate protein is the reference protein, so
/// the stop is written as unchanged, where the codon window gives an insertion
/// of that residue and a stop between residues 10 and 11.
#[test]
fn insertion_before_the_stop_codon_whose_own_stop_stands_at_the_stops_residue_is_the_stop_unchanged(
) {
    let hgvs = insertion_in_the_codon_of_residue_10(true, b"AAGTGTGA");
    assert_eq!(hgvs.offset, Some(1));
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
}

/// The same with one base more, a whole number of codons: the protein is the
/// reference protein and the stop is written as unchanged.
#[test]
fn insertion_of_whole_codons_whose_own_stop_stands_at_the_stops_residue_is_the_stop_unchanged() {
    let hgvs = insertion_in_the_codon_of_residue_10(true, b"AAGTGTGAC");
    assert_eq!(hgvs.offset, Some(1));
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
}

/// Bases inserted two residues before the last that leave the protein as it
/// was: written after the first base of the codon of residue 9, they complete
/// that codon to the same residue, then read as codons of the residues at 10
/// and 11 and a stop codon, which stands at the stop's own residue, 12. The
/// first inserted base equals the codon's second base, so the 3' rule carries
/// the insertion one base on. The alternate protein is the reference protein,
/// so the stop is written as unchanged, where the codon window gives an
/// insertion of those two residues and a stop between residues 9 and 10. One
/// base more, a whole number of codons, reads the same.
#[test]
fn insertion_of_residues_whose_own_stop_stands_at_the_stops_residue_is_the_stop_unchanged() {
    for inserted in [&b"TTAAAGTGTGA"[..], &b"TTAAAGTGTGAC"[..]] {
        let hgvs = insertion_after_base(true, 61_055, inserted);
        assert_eq!(hgvs.offset, Some(1));
        assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
    }
}

/// An insertion that changes the protein: the nine bases of the
/// test above with another codon in the place of residue 11's. Residue 11
/// changes before the inserted stop, the two proteins differ, and the
/// insertion of the new residue and a stop keeps its notation.
#[test]
fn insertion_of_whole_codons_that_changes_the_residue_before_its_stop_keeps_its_notation() {
    let hgvs = insertion_in_the_codon_of_residue_10(true, b"AACTGTGAC");
    assert_eq!(hgvs.offset, Some(1));
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Lys10_Val11insLeuTer")
    );
}

/// The same change in eight bases, which shift the frame: residue
/// 11 changes and the notation is the frameshift named there, its new stop
/// one residue on.
#[test]
fn insertion_that_shifts_the_frame_and_changes_the_residue_before_its_stop_is_a_frameshift() {
    let hgvs = insertion_in_the_codon_of_residue_10(true, b"AACTGTGA");
    assert_eq!(hgvs.offset, Some(1));
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Val11LeufsTer2")
    );
}

/// Where the proteins are not compared: the nine bases that keep
/// both residues, on the transcript whose coding sequence no stop codon
/// closes. The reference protein holds no stop after residue 11, so the
/// inserted stop is no stop left in place, and the insertion of the residue
/// and a stop keeps its notation.
#[test]
fn insertion_ending_in_a_stop_on_a_coding_sequence_no_stop_codon_closes_keeps_its_notation() {
    let hgvs = insertion_in_the_codon_of_residue_10(false, b"AAGTGTGAC");
    assert_eq!(hgvs.offset, Some(1));
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Lys10_Val11insValTer")
    );
}

/// The notation and the row of bases inserted after coding base `after` of a
/// made-up transcript of one exon, on either strand, with `inserted` as the
/// transcript reads it. The coding sequence opens at transcript base 31 and
/// holds twelve codons: those of residue 10 (Lys) and of the last residue, 11
/// (Val), at coding bases 28 and 31, and the stop codon (residue 12) at coding
/// base 34. The 3' UTR opens with a stop codon in the frame that two more
/// bases set.
fn insertion_after_coding_base(
    strand: Strand,
    after: u64,
    inserted: &[u8],
    options: HgvsOptions,
) -> (HgvsNotation, TranscriptConsequence) {
    insertion_after_coding_base_as(strand, after, inserted, options, false)
}

/// [`insertion_after_coding_base`], with the record written one of two ways:
/// as the inserted bases alone between two reference bases, or, under
/// `shared_base`, as the reference base that follows the gap on the forward
/// strand replaced by the inserted bases and that base.
fn insertion_after_coding_base_as(
    strand: Strand,
    after: u64,
    inserted: &[u8],
    options: HgvsOptions,
    shared_base: bool,
) -> (HgvsNotation, TranscriptConsequence) {
    let reverse = strand == Strand::Reverse;
    let model = RowModel {
        chr: "7",
        strand,
        exons: &[(62_001, 62_120)],
        cdna_coding_start: 31,
        cdna_coding_end: 66,
        start_phase: -1,
    };
    // The first coding base lies at 62,031 on the forward strand and at 62,090
    // on the reverse strand, where the bases are laid as their complement.
    let laid: Vec<(u64, Vec<u8>)> = [
        (0, "ATGGCACGTGAC"),
        (12, "AACTGCGAGCAT"),
        (24, "ATCAAGGTCTAA"),
        (36, "CTGATT"),
    ]
    .into_iter()
    .map(|(offset, bases): (u64, &str)| {
        if reverse {
            let at = 62_090 - offset - (bases.len() as u64 - 1);
            (at, reverse_complement(bases.as_bytes()))
        } else {
            (62_031 + offset, bases.as_bytes().to_vec())
        }
    })
    .collect();
    let edits: Vec<(u64, &str)> = laid
        .iter()
        .map(|(at, bases)| (*at, std::str::from_utf8(bases).expect("ASCII bases")))
        .collect();
    let (_dir, fasta) = reference_of(&model, &edits);
    let tx = transcript_of(&model, &fasta);
    let (start, alt) = if reverse {
        (62_091 - after, reverse_complement(inserted))
    } else {
        (62_031 + after, inserted.to_vec())
    };
    let variant = if shared_base {
        let base = laid
            .iter()
            .find_map(|(at, bases)| bases.get(start.checked_sub(*at)? as usize).copied())
            .expect("a laid base after the gap");
        let with_base = [alt.as_slice(), &[base]].concat();
        InputVariant::new("7".into(), start, start, vec![base], with_base)
    } else {
        InputVariant::new("7".into(), start, start - 1, b"-".to_vec(), alt)
    };
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), options);
    let config = EffectsConfig {
        reference_fasta: Some(Arc::new(fasta)),
        hgvs: options,
        ..EffectsConfig::default()
    };
    let row =
        calculate_consequences(&variant, &Arc::new(tx), &config).expect("a row of the transcript");
    (hgvs, row)
}

/// Five bases inserted after the second base of the codon of the last
/// residue, 11, that leave the protein as it was: the first completes that
/// codon to the same residue and the next three are a stop codon, which stands
/// at the stop's own residue, 12. The codon window reads an insertion after
/// the last residue, names no residue after it and gives no description. The
/// alternate protein is the reference protein, so the stop is written as
/// unchanged.
#[test]
fn insertion_of_five_bases_in_the_codon_of_the_last_residue_that_keeps_the_protein_is_the_stop_unchanged(
) {
    let (hgvs, _) =
        insertion_after_coding_base(Strand::Forward, 32, b"ATGAG", HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
}

/// The same with three bases more, eight in all: the protein is the reference
/// protein and the stop is written as unchanged.
#[test]
fn insertion_of_eight_bases_in_the_codon_of_the_last_residue_that_keeps_the_protein_is_the_stop_unchanged(
) {
    let (hgvs, _) =
        insertion_after_coding_base(Strand::Forward, 32, b"ATGAGGAG", HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
}

/// The five bases on the reverse-strand transcript: the stop is written as
/// unchanged.
#[test]
fn insertion_of_five_bases_in_the_codon_of_the_last_residue_that_keeps_the_protein_on_the_reverse_strand(
) {
    let (hgvs, _) =
        insertion_after_coding_base(Strand::Reverse, 32, b"ATGAG", HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
}

/// The eight bases on the reverse-strand transcript: the stop is written as
/// unchanged.
#[test]
fn insertion_of_eight_bases_in_the_codon_of_the_last_residue_that_keeps_the_protein_on_the_reverse_strand(
) {
    let (hgvs, _) =
        insertion_after_coding_base(Strand::Reverse, 32, b"ATGAGGAG", HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
}

/// Five bases inserted after the first base of the codon of the last residue,
/// a copy of the five bases after them, described where the input places them
/// (`--shift_hgvs 0`): the protein is the reference protein and the stop is
/// written as unchanged. The 3' rule would carry the same bases past the stop
/// codon, where no protein notation is written.
#[test]
fn copy_of_the_bases_after_it_inserted_in_the_codon_of_the_last_residue_and_left_in_place_is_the_stop_unchanged(
) {
    let in_place = HgvsOptions {
        shift_hgvs: false,
        ..HgvsOptions::default()
    };
    let (hgvs, _) = insertion_after_coding_base(Strand::Forward, 31, b"TCTAA", in_place);
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
    let (carried, _) =
        insertion_after_coding_base(Strand::Forward, 31, b"TCTAA", HgvsOptions::default());
    assert_eq!(carried.offset, Some(5));
    assert_eq!(carried.hgvsp, None);
}

/// The prediction form (`--hgvsp_use_prediction`): the stop written as
/// unchanged stands inside the parentheses, for the five bases in the codon of
/// the last residue as for the eight bases that end in a stop of their own one
/// codon before it.
#[test]
fn insertion_in_the_codon_of_the_last_residue_that_keeps_the_protein_is_the_stop_unchanged_in_the_prediction_form(
) {
    let prediction = HgvsOptions {
        hgvsp_use_prediction: true,
        ..HgvsOptions::default()
    };
    let (one_codon_before, _) =
        insertion_after_coding_base(Strand::Forward, 28, b"AAGTGTGA", prediction);
    assert_eq!(
        one_codon_before.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.(Ter12=)")
    );
    let (in_the_last_codon, _) =
        insertion_after_coding_base(Strand::Forward, 32, b"ATGAG", prediction);
    assert_eq!(
        in_the_last_codon.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.(Ter12=)")
    );
}

/// A whole number of codons in the same place, six bases and nine: the first
/// completes the codon of the last residue to the same residue and the next
/// three are a stop codon at the stop's own residue, 12. The codon window
/// gives no description of an insertion after the last residue. The alternate
/// protein is the reference protein, so the stop is written as unchanged.
#[test]
fn insertion_of_whole_codons_in_the_codon_of_the_last_residue_that_keeps_the_protein_is_the_stop_unchanged(
) {
    for inserted in [&b"ATGAGG"[..], &b"ATGAGGAGT"[..]] {
        let (hgvs, _) =
            insertion_after_coding_base(Strand::Forward, 32, inserted, HgvsOptions::default());
        assert_eq!(hgvs.offset, None);
        assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
    }
}

/// The six bases and the nine on the reverse-strand transcript: the stop is
/// written as unchanged.
#[test]
fn insertion_of_whole_codons_in_the_codon_of_the_last_residue_that_keeps_the_protein_on_the_reverse_strand(
) {
    for inserted in [&b"ATGAGG"[..], &b"ATGAGGAGT"[..]] {
        let (hgvs, _) =
            insertion_after_coding_base(Strand::Reverse, 32, inserted, HgvsOptions::default());
        assert_eq!(hgvs.offset, None);
        assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
    }
}

/// The six bases in the prediction form (`--hgvsp_use_prediction`): the stop
/// written as unchanged stands inside the parentheses on either strand.
#[test]
fn insertion_of_whole_codons_in_the_codon_of_the_last_residue_that_keeps_the_protein_is_the_stop_unchanged_in_the_prediction_form(
) {
    let prediction = HgvsOptions {
        hgvsp_use_prediction: true,
        ..HgvsOptions::default()
    };
    for strand in [Strand::Forward, Strand::Reverse] {
        let (hgvs, _) = insertion_after_coding_base(strand, 32, b"ATGAGG", prediction);
        assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.(Ter12=)"));
    }
}

/// The same insertions in a record that holds a base its two alleles share:
/// the reference base after the gap on the forward strand, replaced by the
/// inserted bases and that base. The two proteins are those of the insertion
/// alone, at five, six, eight and nine bases on either strand, and the stop
/// is written as unchanged.
#[test]
fn insertion_in_the_codon_of_the_last_residue_held_with_a_shared_base_is_read_as_the_pure_insertion(
) {
    for strand in [Strand::Forward, Strand::Reverse] {
        for inserted in [
            &b"ATGAG"[..],
            &b"ATGAGG"[..],
            &b"ATGAGGAG"[..],
            &b"ATGAGGAGT"[..],
        ] {
            let (pure, _) =
                insertion_after_coding_base(strand, 32, inserted, HgvsOptions::default());
            let (held, _) =
                insertion_after_coding_base_as(strand, 32, inserted, HgvsOptions::default(), true);
            assert_eq!(held.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter12="));
            assert_eq!(held.hgvsp, pure.hgvsp);
        }
    }
}

/// Six bases inserted after the second base of the codon of the last residue
/// that keep that residue and put two residues before the stop codon: the
/// protein gains residues 12 and 13 and its stop moves to residue 14. The two
/// proteins differ, and the insertion keeps what the codon window writes, no
/// notation.
#[test]
fn insertion_of_whole_codons_in_the_codon_of_the_last_residue_that_adds_residues_keeps_what_it_writes(
) {
    let (hgvs, _) =
        insertion_after_coding_base(Strand::Forward, 32, b"AAAGGG", HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp, None);
}

/// Five bases inserted after the first base of the codon of the last residue
/// that change that residue and end in a stop codon at the stop's own residue:
/// the two proteins differ and the change keeps its notation.
#[test]
fn insertion_in_the_codon_of_the_last_residue_that_changes_the_residue_keeps_its_notation() {
    let (hgvs, _) =
        insertion_after_coding_base(Strand::Forward, 31, b"CATGA", HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Val11AlafsTer2")
    );
}

/// Five bases inserted after the second base of the codon of the last residue
/// that keep that residue and put a residue in the stop's place: the protein
/// runs on to a stop codon of the 3' UTR, three residues later, and the change
/// keeps its notation.
#[test]
fn insertion_in_the_codon_of_the_last_residue_that_moves_the_stop_keeps_its_notation() {
    let (hgvs, _) =
        insertion_after_coding_base(Strand::Forward, 32, b"AAAGG", HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter12LysextTer3")
    );
}

/// Five bases inserted between the codon of the last residue and the stop
/// codon, described where the input places them (`--shift_hgvs 0`), that open
/// with a stop codon and leave the protein as it was: the two bases they lie
/// between are not both bases of the codon of the last residue, and no
/// notation is written.
#[test]
fn insertion_between_the_codon_of_the_last_residue_and_the_stop_codon_keeps_what_it_writes() {
    let in_place = HgvsOptions {
        shift_hgvs: false,
        ..HgvsOptions::default()
    };
    let (hgvs, _) = insertion_after_coding_base(Strand::Forward, 33, b"TGAGG", in_place);
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp, None);
}

/// The consequence terms of the insertions that keep the protein, five, six,
/// eight and nine bases on either strand: the notation of the stop as
/// unchanged changes no term of the row.
#[test]
fn insertion_in_the_codon_of_the_last_residue_that_keeps_the_protein_keeps_its_consequence_terms() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for inserted in [
            &b"ATGAG"[..],
            &b"ATGAGG"[..],
            &b"ATGAGGAG"[..],
            &b"ATGAGGAGT"[..],
        ] {
            let (_, row) =
                insertion_after_coding_base(strand, 32, inserted, HgvsOptions::default());
            assert_consequence_set_eq(&row, &["inframe_insertion", "stop_retained_variant"]);
        }
    }
}

/// The stop written as unchanged where it stays: the shape of the
/// row at GRCh37 1:879532 on ENST00000341065, a substitution of CDS base 1768
/// that turns the stop codon (residue 590) into another stop codon. The protein
/// is the reference protein.
#[test]
fn stop_to_stop_substitution_keeps_the_stop_written_as_unchanged() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Forward,
        exons: &[
            (865_692, 865_716),
            (866_419, 866_469),
            (871_152, 871_276),
            (874_420, 874_509),
            (874_655, 874_792),
            (876_524, 876_686),
            (877_516, 877_631),
            (877_790, 877_868),
            (877_939, 878_438),
            (878_633, 878_757),
            (879_078, 879_188),
            (879_288, 879_955),
        ],
        cdna_coding_start: 1,
        cdna_coding_end: 1_769,
        start_phase: 1,
    };
    let (_dir, fasta) = reference_of(&model, &[(879_513, "CCGGAGGTACGAGCCCTATAA")]);
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new("1".into(), 879_532, 879_532, b"A".to_vec(), b"G".to_vec());
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter590="));
}

/// The shape of the row at GRCh37 1:17345377 on ENST00000375499, a
/// reverse-strand transcript: one base inserted after the first base of the
/// stop codon (CDS base 841, residue 281). The whole protein is kept, a residue
/// stands in the stop's place and 22 residues in all follow up to a stop of the
/// 3' UTR: an extension named at the stop, where the codon window gives a
/// deletion-insertion of the stop alone.
#[test]
fn stop_replaced_by_residues_up_to_a_stop_of_the_utr_is_an_extension() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Reverse,
        exons: &[
            (17_380_443, 17_380_665),
            (17_371_256, 17_371_383),
            (17_359_555, 17_359_640),
            (17_355_095, 17_355_231),
            (17_354_244, 17_354_360),
            (17_350_468, 17_350_569),
            (17_349_103, 17_349_225),
            (17_345_217, 17_345_453),
        ],
        cdna_coding_start: 152,
        cdna_coding_end: 994,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (17_345_311, "CTA"),
            (17_345_376, "TTACACACTGCTCTCAAAGGT"),
            (17_380_512, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "1".into(),
        17_345_378,
        17_345_377,
        b"-".to_vec(),
        b"A".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter281LeuextTer22")
    );
}

/// The shape of the row at GRCh37 1:151318416 on ENST00000421986, a
/// reverse-strand transcript: one base inserted after the first base of the
/// stop codon (CDS base 133, residue 45). The whole protein is kept and 80
/// residues follow to the transcript's end without a stop, so the extension's
/// count is unknown.
#[test]
fn stop_replaced_by_residues_without_a_new_stop_is_an_extension_of_unknown_length() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Reverse,
        exons: &[
            (151_319_627, 151_319_682),
            (151_318_965, 151_319_082),
            (151_318_681, 151_318_809),
            (151_318_404, 151_318_433),
            (151_317_582, 151_317_664),
            (151_317_204, 151_317_323),
            (151_316_991, 151_317_012),
        ],
        cdna_coding_start: 188,
        cdna_coding_end: 322,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (151_318_415, "CTAAATAACGACCTCTCCT"),
            (151_318_681, "CG"),
            (151_318_794, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "1".into(),
        151_318_417,
        151_318_416,
        b"-".to_vec(),
        b"G".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter45SerextTer?")
    );
}

/// The shape of the row at GRCh38 9:130500954 on ENST00001100851, given where
/// the 3' rule carries it: 23 bases between the second and third base of the
/// stop codon (after CDS base 842, residue 281). The whole protein is kept and
/// 73 residues follow up to a stop of the 3' UTR.
#[test]
fn insertion_carried_into_the_stop_codon_is_an_extension_at_its_shifted_place() {
    let model = RowModel {
        chr: "9",
        strand: Strand::Forward,
        exons: &[
            (130_444_888, 130_444_995),
            (130_452_224, 130_452_333),
            (130_454_305, 130_454_373),
            (130_458_401, 130_458_589),
            (130_464_111, 130_464_167),
            (130_466_725, 130_466_799),
            (130_470_834, 130_470_904),
            (130_471_485, 130_471_515),
            (130_476_871, 130_476_961),
            (130_479_716, 130_479_800),
            (130_480_385, 130_480_449),
            (130_500_976, 130_502_361),
        ],
        cdna_coding_start: 114,
        cdna_coding_end: 956,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (130_452_229, "ATG"),
            (130_480_434, "CAGCTGAGTCCGCCCG"),
            (130_500_976, "ATTAA"),
            (130_501_174, "TAA"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "9".into(),
        130_500_980,
        130_500_979,
        b"-".to_vec(),
        b"CTTGAGGGGCGTCACCCATCGTT".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter281TyrextTer73")
    );
}

/// An extension named at the stop: the shape of the row at GRCh37
/// 1:1153052 on ENST00000263741, a substitution of CDS base 1045, the first
/// base of the stop codon (residue 349), after which 6 residues follow the
/// whole protein up to a stop of the 3' UTR.
#[test]
fn stop_lost_by_a_substitution_keeps_its_extension() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Reverse,
        exons: &[
            (1_167_272, 1_167_411),
            (1_163_848, 1_164_326),
            (1_159_212, 1_159_348),
            (1_158_624, 1_158_737),
            (1_154_149, 1_154_307),
            (1_153_838, 1_154_013),
            (1_152_311, 1_153_184),
        ],
        cdna_coding_start: 294,
        cdna_coding_end: 1_340,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (1_153_032, "CTA"),
            (1_153_050, "CTAGGCAGTTTTTGCGCAACC"),
            (1_164_171, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "1".into(),
        1_153_052,
        1_153_052,
        b"A".to_vec(),
        b"G".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter349GlnextTer6")
    );
}

/// The shape of the row at GRCh38 3:81591146 on ENST00001019268, a
/// reverse-strand transcript: 41 bases inserted after the second base of the
/// stop codon (CDS base 1130, residue 377). The whole protein is kept and 10
/// residues follow up to a stop codon that lies wholly inside the inserted
/// bases: an extension named at the stop, where the codon window gives a
/// deletion-insertion of the stop that names every added residue and that stop.
#[test]
fn inserted_bases_that_hold_the_new_stop_extend_the_protein() {
    let model = RowModel {
        chr: "3",
        strand: Strand::Reverse,
        exons: &[
            (81_761_375, 81_761_642),
            (81_705_444, 81_705_613),
            (81_670_838, 81_670_953),
            (81_649_796, 81_649_921),
            (81_648_856, 81_648_991),
            (81_646_392, 81_646_482),
            (81_642_781, 81_642_990),
            (81_593_904, 81_594_023),
            (81_591_037, 81_591_164),
            (81_586_092, 81_586_190),
            (81_581_165, 81_581_275),
            (81_577_925, 81_578_096),
            (81_536_911, 81_537_095),
            (81_535_195, 81_535_325),
            (81_499_110, 81_499_227),
            (81_489_707, 81_490_463),
        ],
        cdna_coding_start: 126,
        cdna_coding_end: 1_256,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (81_591_146, "CTAATTTTCTGGACACCCG"),
            (81_593_904, "GC"),
            (81_761_515, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "3".into(),
        81_591_147,
        81_591_146,
        b"-".to_vec(),
        b"GGTGTGGTGTTCATCCGCGTACCCACCACCGGCGAAGTCTA".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter377TyrextTer10")
    );
}

/// The shape of the row at GRCh38 7:5997387 on ENST00000699757, a
/// reverse-strand transcript: ten bases inserted after the first base of the
/// stop codon (CDS base 190, residue 64). The whole protein is kept and 2
/// residues follow up to a stop codon wholly inside the inserted bases.
#[test]
fn inserted_bases_that_hold_the_new_stop_two_residues_on_extend_the_protein() {
    let model = RowModel {
        chr: "7",
        strand: Strand::Reverse,
        exons: &[
            (6_008_997, 6_009_084),
            (6_005_892, 6_006_026),
            (6_003_972, 6_004_058),
            (5_999_108, 5_999_226),
            (5_997_326, 5_997_423),
            (5_995_534, 5_995_633),
            (5_991_973, 5_992_057),
            (5_989_800, 5_989_955),
            (5_986_759, 5_987_620),
            (5_982_824, 5_982_991),
            (5_978_596, 5_978_696),
            (5_977_588, 5_977_757),
            (5_970_955, 5_973_542),
        ],
        cdna_coding_start: 276,
        cdna_coding_end: 467,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[(5_997_386, "TTATAAGTTAGACGGCTGACC"), (6_004_004, "CAT")],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "7".into(),
        5_997_388,
        5_997_387,
        b"-".to_vec(),
        b"ATTTAACACA".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter64LeuextTer2")
    );
}

/// The shape of the row at GRCh37 1:865692 on ENST00000341065, whose coding
/// sequence opens on the second base of its first codon and is closed by a stop
/// codon. A substitution of CDS base 1 changes that codon alone; the codon
/// window reads the unknown residue as the stop and writes it as unchanged at
/// residue 1. The reference protein holds no stop there, so no protein notation
/// is written.
#[test]
fn stop_unchanged_at_residue_1_of_an_incomplete_first_codon_has_no_protein_notation() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Forward,
        exons: &[
            (865_692, 865_716),
            (866_419, 866_469),
            (871_152, 871_276),
            (874_420, 874_509),
            (874_655, 874_792),
            (876_524, 876_686),
            (877_516, 877_631),
            (877_790, 877_868),
            (877_939, 878_438),
            (878_633, 878_757),
            (879_078, 879_188),
            (879_288, 879_955),
        ],
        cdna_coding_start: 1,
        cdna_coding_end: 1_769,
        start_phase: 1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[(865_692, "C"), (879_513, "CTGGTAGGACTGTCGATGTGA")],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new("1".into(), 865_692, 865_692, b"C".to_vec(), b"A".to_vec());
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp, None);
}

/// The shape of the row at GRCh37 16:89815072 on ENST00000561660, a
/// reverse-strand transcript whose coding sequence opens inside its first
/// codon: three bases inserted after the first base of the stop codon (CDS base
/// 546, residue 183). The whole protein is kept, 1 residue follows it, and the
/// last inserted base re-forms a stop codon with the stop codon's last two
/// bases: an extension named at the stop, where the codon window gives an
/// insertion between the last residue and the stop.
#[test]
fn residue_inserted_before_a_stop_re_formed_with_one_inserted_base_is_an_extension() {
    let model = RowModel {
        chr: "16",
        strand: Strand::Reverse,
        exons: &[
            (89_831_298, 89_831_452),
            (89_828_357, 89_828_430),
            (89_824_985, 89_825_113),
            (89_818_546, 89_818_630),
            (89_815_067, 89_815_175),
            (89_813_239, 89_813_298),
            (89_812_992, 89_813_096),
            (89_811_435, 89_811_479),
        ],
        cdna_coding_start: 1,
        cdna_coding_end: 548,
        start_phase: 1,
    };
    let (_dir, fasta) = reference_of(&model, &[(89_815_071, "CTAAGGAGGGGAACAATCACC")]);
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "16".into(),
        89_815_073,
        89_815_072,
        b"-".to_vec(),
        b"AGG".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter183SerextTer1")
    );
}

/// The shape of the row at GRCh37 16:79245647 on ENST00000406884: 45 bases
/// inserted after CDS base 659, a copy of the 45 after them, which the 3' rule
/// carries 45 bases on, to between the second and third base of the stop codon
/// (residue 235). The whole protein is kept, 15 residues follow it, and the
/// last two inserted bases re-form a stop codon with the stop codon's last
/// base.
#[test]
fn residues_inserted_before_a_stop_re_formed_with_two_inserted_bases_are_an_extension() {
    let model = RowModel {
        chr: "16",
        strand: Strand::Forward,
        exons: &[
            (78_133_551, 78_133_782),
            (78_142_320, 78_142_384),
            (78_143_675, 78_143_732),
            (78_148_873, 78_149_051),
            (78_198_080, 78_198_186),
            (79_245_505, 79_246_564),
        ],
        cdna_coding_start: 126,
        cdna_coding_end: 830,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (78_133_676, "ATG"),
            (
                79_245_637,
                "GTCCCTCGTTGTGTCGTAGACCAACTCGTTTTATCTACGGCAGATACTGCGCTCTAG",
            ),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "16".into(),
        79_245_648,
        79_245_647,
        b"-".to_vec(),
        b"TGTCGTAGACCAACTCGTTTTATCTACGGCAGATACTGCGCTCTA".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, Some(45));
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter235TyrextTer15")
    );
}

/// The shape of the row at GRCh37 3:37061788 on ENST00000447829, whose coding
/// sequence opens inside its first codon, given where the 3' rule carries it:
/// twelve bases after the first base of the stop codon (CDS base 174, residue
/// 59). The whole protein is kept, 4 residues follow it, and the last inserted
/// base re-forms a stop codon with the stop codon's last two bases.
#[test]
fn insertion_carried_to_the_stop_codon_before_a_re_formed_stop_is_an_extension() {
    let model = RowModel {
        chr: "3",
        strand: Strand::Forward,
        exons: &[
            (37_053_533, 37_053_590),
            (37_055_923, 37_056_035),
            (37_061_801, 37_061_954),
            (37_067_128, 37_067_410),
        ],
        cdna_coding_start: 1,
        cdna_coding_end: 176,
        start_phase: 1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[(37_056_020, "TTGCTATGTGACGCAT"), (37_061_801, "ACTAG")],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "3".into(),
        37_061_804,
        37_061_803,
        b"-".to_vec(),
        b"TGGAGGATCGAT".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter59LeuextTer4")
    );
}

/// A lost stop on a transcript that ends before a whole codon
/// follows: the shape of the row at GRCh37 22:29083885 on ENST00000403642,
/// whose last base is the stop codon's last (residue 453). Deleting CDS base
/// 1358 leaves two bases after the last residue, no residue is added, and the
/// notation is left as the codon window gives it.
#[test]
fn lost_stop_with_no_whole_codon_after_it_keeps_its_notation() {
    let model = RowModel {
        chr: "22",
        strand: Strand::Reverse,
        exons: &[
            (29_130_391, 29_130_709),
            (29_115_383, 29_115_473),
            (29_107_897, 29_108_005),
            (29_105_994, 29_106_047),
            (29_099_493, 29_099_554),
            (29_095_826, 29_095_925),
            (29_092_889, 29_092_975),
            (29_091_698, 29_091_861),
            (29_091_115, 29_091_230),
            (29_090_020, 29_090_105),
            (29_085_123, 29_085_203),
            (29_083_885, 29_083_974),
        ],
        cdna_coding_start: 1,
        cdna_coding_end: 1_359,
        start_phase: 0,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[(29_083_885, "TCATCGCGGGACGGTTCGTAG"), (29_130_707, "CAT")],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "22".into(),
        29_083_886,
        29_083_886,
        b"C".to_vec(),
        b"-".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ter453="));
}

/// The shape of the row at GRCh37 1:2238201 on ENST00000378536: the first base
/// of the stop codon deleted (CDS base 2185, residue 729). The stop is lost,
/// the whole protein is kept and 47 residues follow it up to a stop of the 3'
/// UTR: an extension named at the stop, where the codon window writes the stop
/// as unchanged.
#[test]
fn deletion_in_the_stop_codon_that_loses_the_stop_is_an_extension() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Forward,
        exons: &[
            (2_160_134, 2_161_174),
            (2_234_417, 2_234_542),
            (2_234_724, 2_234_839),
            (2_235_279, 2_235_541),
            (2_235_732, 2_236_024),
            (2_237_459, 2_237_689),
            (2_238_016, 2_241_558),
        ],
        cdna_coding_start: 73,
        cdna_coding_end: 2_259,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (2_160_206, "ATG"),
            (2_238_184, "ATAATTCTAACTACTCCCTGA"),
            (2_238_344, "TAG"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "1".into(),
        2_238_202,
        2_238_202,
        b"T".to_vec(),
        b"-".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter729AspextTer47")
    );
}

/// The shape of the row at GRCh37 11:31811483 on ENST00000241001, a
/// reverse-strand transcript: the first base of the stop codon deleted (CDS
/// base 1267, residue 423). The stop is lost and 7 residues follow the whole
/// protein to the transcript's end without a stop, so the extension's count is
/// unknown.
#[test]
fn deletion_in_the_stop_codon_without_a_new_stop_is_an_extension_of_unknown_length() {
    let model = RowModel {
        chr: "11",
        strand: Strand::Reverse,
        exons: &[
            (31_832_761, 31_832_785),
            (31_832_376, 31_832_563),
            (31_828_397, 31_828_473),
            (31_827_950, 31_828_010),
            (31_824_252, 31_824_382),
            (31_823_109, 31_823_324),
            (31_822_239, 31_822_404),
            (31_816_178, 31_816_336),
            (31_815_580, 31_815_662),
            (31_815_200, 31_815_350),
            (31_814_986, 31_815_101),
            (31_812_258, 31_812_408),
            (31_811_461, 31_811_567),
        ],
        cdna_coding_start: 342,
        cdna_coding_end: 1_610,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[(31_811_482, "TCATATGGAAGGTGAGAGGAG"), (31_827_957, "CAT")],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "11".into(),
        31_811_484,
        31_811_484,
        b"A".to_vec(),
        b"-".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter423GluextTer?")
    );
}

/// The shape of the row at GRCh38 12:112450359 on ENST00001096176: CDS bases
/// 174 to 177, the last base of the last residue's codon and the whole stop
/// codon (residue 59), replaced by one base. Residue 58 changes, the stop is
/// lost and the alternate protein meets a stop of the 3' UTR at position 4
/// counted from the changed residue: a frameshift named at residue 58, where
/// the codon window gives a deletion-insertion of the range that ends at the
/// stop.
#[test]
fn change_from_the_last_codon_to_the_stop_codons_last_base_is_a_frameshift() {
    let model = RowModel {
        chr: "12",
        strand: Strand::Forward,
        exons: &[
            (112_418_904, 112_419_125),
            (112_446_276, 112_446_393),
            (112_450_318, 112_450_512),
            (112_453_195, 112_453_249),
        ],
        cdna_coding_start: 209,
        cdna_coding_end: 385,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (112_419_112, "ATG"),
            (112_450_342, "ATAAACTCGATGCCCAACTGACGGTCCTAA"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "12".into(),
        112_450_359,
        112_450_362,
        b"CTGA".to_vec(),
        b"G".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asn58LysfsTer4")
    );
}

/// The shape of the row at GRCh38 3:128908990 on ENST00000681585: CDS bases
/// 1426 to 1431, the last residue's codon and the stop codon (residue 477),
/// replaced by three bases. Residue 476 changes, the stop is lost and the
/// alternate protein meets a stop of the 3' UTR at position 24 counted from the
/// changed residue: a frameshift named at residue 476.
#[test]
fn last_codon_and_stop_codon_replaced_by_one_codon_are_a_frameshift() {
    let model = RowModel {
        chr: "3",
        strand: Strand::Forward,
        exons: &[
            (128_879_614, 128_879_841),
            (128_884_653, 128_884_746),
            (128_893_555, 128_893_656),
            (128_895_310, 128_895_416),
            (128_896_436, 128_896_536),
            (128_897_632, 128_897_710),
            (128_899_287, 128_899_461),
            (128_901_276, 128_901_349),
            (128_902_553, 128_902_628),
            (128_904_062, 128_904_132),
            (128_904_386, 128_904_505),
            (128_906_121, 128_906_249),
            (128_908_185, 128_908_264),
            (128_908_923, 128_909_099),
            (128_909_344, 128_909_421),
            (128_910_021, 128_910_149),
            (128_910_741, 128_910_813),
            (128_912_507, 128_913_099),
        ],
        cdna_coding_start: 79,
        cdna_coding_end: 1_509,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (128_879_692, "ATG"),
            (128_908_975, "CGACCTAGCCCCTCGAAGTGATCCCGTCGGGCCCTCTTTCCCTGGCCCAGAAAGAACGACTATTCAGAGCGTCGGCCCAAGCAAGATTAG"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "3".into(),
        128_908_990,
        128_908_995,
        b"AAGTGA".to_vec(),
        b"GAT".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Lys476AspfsTer24")
    );
}

/// A change into the stop codon that shifts the frame, after which the protein
/// keeps its length: the shape of the row at GRCh38 9:127818308 on
/// ENST00001020284, a reverse-strand transcript. CDS bases 1483 to 1493, up to
/// the second base of the stop codon (residue 498), are replaced by ten bases.
/// The residues from 495 on are read from other bases and a stop codon forms at
/// the stop's own residue; the change shifts the frame all the same, so it
/// reads as a frameshift named at the first changed residue.
#[test]
fn change_into_the_stop_codon_that_keeps_the_proteins_length_is_a_frameshift() {
    let model = RowModel {
        chr: "9",
        strand: Strand::Reverse,
        exons: &[
            (127_854_289, 127_854_719),
            (127_843_094, 127_843_245),
            (127_829_687, 127_829_827),
            (127_826_510, 127_826_672),
            (127_825_695, 127_825_860),
            (127_825_231, 127_825_357),
            (127_824_800, 127_824_974),
            (127_824_304, 127_824_446),
            (127_819_900, 127_820_037),
            (127_819_622, 127_819_660),
            (127_818_716, 127_818_832),
            (127_818_120, 127_818_372),
            (127_817_149, 127_817_203),
            (127_815_943, 127_816_053),
            (127_815_016, 127_815_806),
        ],
        cdna_coding_start: 365,
        cdna_coding_end: 1_858,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (127_818_306, "TCTAGAGCCTATCCCGGCAGGC"),
            (127_854_353, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "9".into(),
        127_818_308,
        127_818_318,
        b"TAGAGCCTATC".to_vec(),
        b"AAAATGTGTA".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Asp495TyrfsTer4")
    );
}

/// A deletion through the stop codon that shifts the frame, after which the
/// protein keeps its length: the shape of the row at GRCh38 X:49076647 on
/// ENST00001057692, a reverse-strand transcript. CDS bases 264 to 270, from the
/// last base of residue 88's codon to the last base of the stop codon (residue
/// 90), are deleted. Residue 88 is read as before, residue 89 changes and a
/// stop codon of the 3' UTR stands at the stop's own residue; the deletion
/// shifts the frame all the same, so it reads as a frameshift named at residue
/// 89, the first changed residue.
#[test]
fn deletion_through_the_stop_codon_that_keeps_the_proteins_length_is_a_frameshift() {
    let model = RowModel {
        chr: "X",
        strand: Strand::Reverse,
        exons: &[
            (49_079_751, 49_079_882),
            (49_078_041, 49_078_112),
            (49_077_837, 49_077_911),
            (49_077_575, 49_077_611),
            (49_076_617, 49_076_750),
            (49_076_430, 49_076_524),
            (49_075_866, 49_075_945),
            (49_075_545, 49_075_753),
            (49_075_364, 49_075_465),
            (49_075_136, 49_075_281),
            (49_074_445, 49_074_912),
        ],
        cdna_coding_start: 150,
        cdna_coding_end: 419,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (49_076_641, "TCATCCCTTACGCTCTGACCAACCGCTG"),
            (49_078_093, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "X".into(),
        49_076_648,
        49_076_654,
        b"TTACGCT".to_vec(),
        b"-".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala89GlyfsTer2")
    );
}

/// A change of a whole number of codons that ends before the stop
/// codon's last base: the shape of the row at GRCh37 9:35077384 on
/// ENST00000425676, a reverse-strand transcript. CDS bases 317 to 320, the last
/// two bases of the last residue's codon and the first two of the stop codon
/// (residue 107), are replaced by one base. The nomenclature's frameshift
/// sentences do not reach it, and the notation is left as the codon window
/// gives it.
#[test]
fn in_frame_change_short_of_the_stop_codons_last_base_keeps_its_notation() {
    let model = RowModel {
        chr: "9",
        strand: Strand::Reverse,
        exons: &[
            (35_079_438, 35_079_690),
            (35_079_148, 35_079_238),
            (35_078_602, 35_078_733),
            (35_077_261, 35_077_396),
            (35_076_968, 35_077_098),
            (35_076_721, 35_076_867),
            (35_076_429, 35_076_580),
            (35_075_959, 35_076_025),
            (35_075_462, 35_075_751),
            (35_075_276, 35_075_322),
            (35_074_924, 35_075_079),
            (35_074_368, 35_074_491),
            (35_073_832, 35_074_213),
        ],
        cdna_coding_start: 170,
        cdna_coding_end: 490,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (
                35_077_341,
                "TCAGAACCAATGCGCGACTAGCTTGAACGGCACCTTAGAAACTCAGAAGGTACTAG",
            ),
            (35_078_602, "TGAGCCG"),
            (35_079_519, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "9".into(),
        35_077_384,
        35_077_387,
        b"CAGA".to_vec(),
        b"G".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Phe106_Ter107delinsSer")
    );
}

/// A record whose two alleles share their last base: a record at
/// the place of the row at GRCh38 12:112450359 on ENST00001096176, CDS bases
/// 174 to 177 replaced by one base that equals base 177. The change itself is
/// the deletion of bases 174 to 176, a whole number of codons that ends before
/// the last base of the stop codon (residue 59), so no frameshift is written
/// and the notation is left as the codon window gives it.
#[test]
fn record_whose_alleles_share_their_last_base_is_read_by_the_change_itself() {
    let model = RowModel {
        chr: "12",
        strand: Strand::Forward,
        exons: &[
            (112_418_904, 112_419_125),
            (112_446_276, 112_446_393),
            (112_450_318, 112_450_512),
            (112_453_195, 112_453_249),
        ],
        cdna_coding_start: 209,
        cdna_coding_end: 385,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (112_419_112, "ATG"),
            (112_450_342, "CATAAGAGAGGCAAGTGTTAGTGGGCTTGA"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "12".into(),
        112_450_359,
        112_450_362,
        b"TTAG".to_vec(),
        b"G".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Cys58_Ter59delinsTrp")
    );
}

/// The shape of the row at GRCh37 1:3548777 on ENST00000424367, a
/// reverse-strand transcript whose coding sequence ends on the first base of a
/// codon and is closed by no stop codon: a substitution of that base, CDS base
/// 913. The codon window reads the partial codon as the stop and writes it as
/// unchanged at residue 305. The reference protein holds no stop there, so no
/// protein notation is written.
#[test]
fn stop_unchanged_on_a_coding_sequence_no_stop_codon_closes_has_no_protein_notation() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Reverse,
        exons: &[
            (3_566_495, 3_566_590),
            (3_563_972, 3_564_124),
            (3_563_230, 3_563_346),
            (3_555_314, 3_555_386),
            (3_553_559, 3_553_662),
            (3_552_508, 3_552_594),
            (3_551_561, 3_551_638),
            (3_549_962, 3_550_067),
            (3_548_777, 3_548_902),
        ],
        cdna_coding_start: 28,
        cdna_coding_end: 940,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[(3_548_777, "CTTTACGGATTGGCCCCGCGATGGT"), (3_566_561, "CAT")],
    );
    let mut tx = transcript_of(&model, &fasta);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let variant = InputVariant::new(
        "1".into(),
        3_548_777,
        3_548_777,
        b"C".to_vec(),
        b"T".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp, None);
}

/// The shape of the row at GRCh37 1:224492436 on ENST00000488718, a
/// reverse-strand transcript whose coding sequence ends on the first base of a
/// codon and is closed by no stop codon: that base, CDS base 553, deleted. The
/// codon window writes the stop at residue 185 as deleted. The reference
/// protein holds no stop there, so no protein notation is written.
#[test]
fn stop_deleted_on_a_coding_sequence_no_stop_codon_closes_has_no_protein_notation() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Reverse,
        exons: &[
            (224_517_773, 224_517_840),
            (224_514_093, 224_514_166),
            (224_505_580, 224_505_632),
            (224_505_396, 224_505_495),
            (224_499_502, 224_499_559),
            (224_492_736, 224_492_868),
            (224_492_437, 224_492_514),
        ],
        cdna_coding_start: 12,
        cdna_coding_end: 564,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (224_492_437, "CTTGTAAGAATTGACTATACCACGA"),
            (224_517_827, "CAT"),
        ],
    );
    let mut tx = transcript_of(&model, &fasta);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let variant = InputVariant::new(
        "1".into(),
        224_492_437,
        224_492_437,
        b"C".to_vec(),
        b"-".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp, None);
}

/// The shape of the row at GRCh37 11:62459848 on ENST00000448568, a
/// reverse-strand transcript whose coding sequence ends on the second base of a
/// codon and is closed by no stop codon: one base inserted between those two
/// bases, after CDS base 670, which makes the codon whole. The codon window
/// writes the stop at residue 224 as replaced by one residue. The reference
/// protein holds no stop there, so no protein notation is written.
#[test]
fn stop_replaced_by_one_residue_on_a_coding_sequence_no_stop_codon_closes_has_no_protein_notation()
{
    let model = RowModel {
        chr: "11",
        strand: Strand::Reverse,
        exons: &[
            (62_473_475, 62_473_690),
            (62_472_773, 62_473_089),
            (62_469_940, 62_470_021),
            (62_462_040, 62_462_183),
            (62_460_135, 62_460_269),
            (62_459_848, 62_459_945),
        ],
        cdna_coding_start: 322,
        cdna_coding_end: 992,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (62_459_848, "CACTGGCGAGGGATTCCTCCGTGGGG"),
            (62_472_982, "CAT"),
        ],
    );
    let mut tx = transcript_of(&model, &fasta);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let variant = InputVariant::new(
        "11".into(),
        62_459_849,
        62_459_848,
        b"-".to_vec(),
        b"A".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp, None);
}

/// The shape of the row at GRCh37 16:28498985 on ENST00000568443, a
/// reverse-strand transcript whose coding sequence ends on the second base of a
/// codon and is closed by no stop codon: one base inserted after CDS base 524,
/// inside the last whole codon, whose residue stays while the bases after it
/// make a whole codon. The frameshift walk closes the reference translation
/// with a stop, so the description is an extension named at a stop at residue
/// 176. The reference protein holds no stop there, so no protein notation is
/// written.
#[test]
fn extension_on_a_coding_sequence_no_stop_codon_closes_has_no_protein_notation() {
    let model = RowModel {
        chr: "16",
        strand: Strand::Reverse,
        exons: &[
            (28_503_341, 28_503_623),
            (28_502_803, 28_503_156),
            (28_500_611, 28_500_707),
            (28_499_912, 28_499_983),
            (28_498_983, 28_499_062),
        ],
        cdna_coding_start: 360,
        cdna_coding_end: 886,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (28_498_983, "ATGGCTAGTAGCCACCAGGAGATCAT"),
            (28_503_078, "CAT"),
        ],
    );
    let mut tx = transcript_of(&model, &fasta);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let variant = InputVariant::new(
        "16".into(),
        28_498_986,
        28_498_985,
        b"-".to_vec(),
        b"C".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(hgvs.hgvsp, None);
}

/// A made-up forward-strand transcript that no flag marks incomplete: its
/// coding sequence is three whole codons and two more bases, and its last three
/// coding bases spell a stop codon in no frame of the protein. One base inserted
/// after CDS base 7, inside the last whole codon, is a copy of the two bases
/// after it, so the 3' rule carries it 2 bases on, to between CDS bases 9 and
/// 10, where it closes a stop codon with the two bases the coding sequence ends
/// on. `_get_fs_peptides` closes the reference translation with a stop whatever
/// the coding sequence's end (`TranscriptVariationAllele.pm:2338`), so the
/// frameshift walk writes the stop as unchanged at residue 4. The reference
/// protein holds three residues and no stop, so no protein notation is written.
#[test]
fn frameshift_past_the_last_whole_codon_of_an_unflagged_coding_sequence_has_no_protein_notation() {
    let model = RowModel {
        chr: "3",
        strand: Strand::Forward,
        exons: &[(10_001, 10_040)],
        cdna_coding_start: 6,
        cdna_coding_end: 16,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(&model, &[(10_001, "GCGCCATGGGCTTTAAGCC")]);
    let tx = transcript_of(&model, &fasta);
    assert!(tx.flags.is_empty());
    let variant = InputVariant::new("3".into(), 10_013, 10_012, b"-".to_vec(), b"T".to_vec());
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert!(hgvs
        .hgvsc
        .as_deref()
        .is_some_and(|c| c.ends_with(":c.9dup")));
    assert_eq!(hgvs.offset, Some(2));
    assert_eq!(hgvs.hgvsp, None);
}

/// A stop that the transcript's own peptide carries: the shape of
/// the row at GRCh37 12:6940463 on ENST00000606935, whose coding sequence is
/// closed by no stop codon and holds a stop codon at residue 356. A
/// substitution of CDS base 1068, that codon's last base, loses it. The
/// extension names a stop the reference protein holds, and is left as written.
#[test]
fn extension_at_a_stop_the_peptide_carries_keeps_its_notation() {
    let model = RowModel {
        chr: "12",
        strand: Strand::Forward,
        exons: &[
            (6_937_606, 6_938_102),
            (6_939_023, 6_939_175),
            (6_939_517, 6_939_718),
            (6_939_803, 6_939_934),
            (6_940_380, 6_940_516),
            (6_940_989, 6_941_078),
            (6_942_756, 6_942_808),
            (6_942_906, 6_942_973),
            (6_943_089, 6_943_213),
            (6_946_149, 6_946_250),
            (6_946_587, 6_946_737),
            (6_946_896, 6_947_013),
            (6_947_122, 6_947_197),
            (6_948_164, 6_948_304),
            (6_948_461, 6_948_625),
        ],
        cdna_coding_start: 1,
        cdna_coding_end: 2_210,
        start_phase: 0,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (6_937_606, "ATG"),
            (6_940_449, "TGCGGTCCACAATAATATGATCCT"),
            (6_948_600, "AGAAAGCCCTCTCGGTACTCCATGCG"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "12".into(),
        6_940_463,
        6_940_463,
        b"A".to_vec(),
        b"T".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter356TyrextTer?")
    );
}

/// A coding sequence that ends inside a codon while its last three bases
/// spell a stop codon: the shape of the rows at GRCh37 22:23235998 on
/// ENST00000390320, whose coding sequence is one base longer than whole
/// codons from its first base, as the exon's phase reads it, and holds no
/// other stop codon in the frame of its last three bases. A substitution of
/// its last base, CDS base 127, changes the lone first base of the codon of
/// residue 43. The codon window writes that codon's unknown residue as the
/// stop for the reference and for every allele alike
/// (`TranscriptVariationAllele.pm:774-776`, `:2181-2182`): the stop unchanged
/// at residue 43, whether the allele leaves a stop codon in the last three
/// bases or turns them into a residue's codon. The reference protein holds 42
/// residues and no stop, so no protein notation is written.
#[test]
fn stop_unchanged_where_the_last_three_bases_spell_a_stop_codon_out_of_frame_has_no_protein_notation(
) {
    let model = RowModel {
        chr: "22",
        strand: Strand::Forward,
        exons: &[(23_235_872, 23_235_998)],
        cdna_coding_start: 1,
        cdna_coding_end: 127,
        start_phase: 0,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (23_235_872, "GCATCGTTGGCTCCGAATTCGTATATGGCGCCCCGCTTTCGGCATTCGGAGCTGTGCTTCAAGGTGGCACATTCTTGGCGTGGCATTAGGACGCAAAGTTATGCCGACAGGAGCAATTCCAAACTAA"),
        ],
    );
    let mut tx = transcript_of(&model, &fasta);
    tx.flags = vec!["cds_start_NF".to_string()].into();
    // The first allele leaves a stop codon in the last three bases; the other
    // two turn them into a residue's codon.
    for alternate in [b"G", b"C", b"T"] {
        let variant = InputVariant::new(
            "22".into(),
            23_235_998,
            23_235_998,
            b"A".to_vec(),
            alternate.to_vec(),
        );
        let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
        assert_eq!(hgvs.offset, None);
        assert_eq!(hgvs.hgvsp, None, "alternate {}", alternate[0] as char);
    }
}

/// Another description on a coding sequence no stop codon closes:
/// the shape of the row at GRCh37 20:58452518 on ENST00000446834, a
/// reverse-strand transcript whose coding sequence ends on the second base of a
/// codon. CDS bases 3067 to 3071, its last whole codon and the two bases after
/// it, are deleted. The frameshift walk writes a deletion of the range the
/// codon window gives, from residue 1023 to a stop at residue 1024, a
/// description that does not name the stop alone, and it is left as written.
#[test]
fn range_deleted_up_to_the_end_of_an_unclosed_coding_sequence_keeps_its_notation() {
    let model = RowModel {
        chr: "20",
        strand: Strand::Reverse,
        exons: &[
            (58_508_286, 58_508_416),
            (58_507_117, 58_507_209),
            (58_497_445, 58_497_514),
            (58_496_365, 58_496_508),
            (58_495_414, 58_495_542),
            (58_494_548, 58_494_652),
            (58_494_213, 58_494_237),
            (58_491_541, 58_491_626),
            (58_490_520, 58_490_605),
            (58_490_220, 58_490_254),
            (58_489_164, 58_489_306),
            (58_489_030, 58_489_082),
            (58_487_428, 58_487_473),
            (58_486_805, 58_486_900),
            (58_482_419, 58_482_479),
            (58_476_717, 58_476_865),
            (58_475_799, 58_475_873),
            (58_475_193, 58_475_339),
            (58_471_484, 58_471_583),
            (58_470_492, 58_470_652),
            (58_468_579, 58_468_648),
            (58_468_197, 58_468_254),
            (58_467_623, 58_467_756),
            (58_467_045, 58_467_481),
            (58_461_762, 58_461_889),
            (58_460_983, 58_461_019),
            (58_457_150, 58_457_220),
            (58_457_000, 58_457_057),
            (58_456_492, 58_456_547),
            (58_456_021, 58_456_078),
            (58_455_358, 58_455_526),
            (58_453_063, 58_453_100),
            (58_452_519, 58_452_610),
        ],
        cdna_coding_start: 271,
        cdna_coding_end: 3_341,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (58_452_519, "ACAGGGGCGGCGATCTCGATTTCATT"),
            (58_497_466, "CAT"),
        ],
    );
    let mut tx = transcript_of(&model, &fasta);
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let variant = InputVariant::new(
        "20".into(),
        58_452_519,
        58_452_523,
        b"ACAGG".to_vec(),
        b"-".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert!(hgvs
        .hgvsp
        .as_deref()
        .is_some_and(|notation| notation.ends_with("1023_Ter1024del")));
}

/// The shape of the row at GRCh37 1:26142208 on ENST00000374315, whose coding
/// sequence holds a stop codon at residue 428 that the transcript's peptide
/// reads as a residue: a substitution of CDS base 1670, in the closing stop
/// codon (residue 557), loses the stop, and 35 residues follow the whole
/// protein up to a stop of the 3' UTR that is another codon than the one
/// inside. The count ported from Ensembl VEP ends at the stop codon inside and
/// is unknown; the extension carries the count of residues added.
#[test]
fn stop_lost_past_a_stop_codon_read_as_a_residue_counts_the_residues_added() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Forward,
        exons: &[
            (26_126_684, 26_126_904),
            (26_127_534, 26_127_651),
            (26_131_633, 26_131_766),
            (26_135_071, 26_135_280),
            (26_135_517, 26_135_641),
            (26_136_174, 26_136_311),
            (26_137_945, 26_138_026),
            (26_138_182, 26_138_370),
            (26_139_178, 26_139_283),
            (26_140_372, 26_140_484),
            (26_140_568, 26_140_669),
            (26_142_039, 26_144_712),
        ],
        cdna_coding_start: 39,
        cdna_coding_end: 1_709,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (26_126_722, "ATG"),
            (26_139_271, "TTGATTGTCTAGC"),
            (26_140_372, "TCTGGGTA"),
            (26_142_183, "GATGGGAGAAGGTTAACAGCTGTTTAGAGTCAAAGTTGTCTTTCTTCGGCTTGTGTACATAAGGCGTATTATCATCTCAGGTCCTCTAGCGCTCAATTGGCACATCTTCGTACCCTTAGGAGCATGTTCTGA"),
        ],
    );
    let mut tx = transcript_of(&model, &fasta);
    // The model's peptide reads the stop codon inside the coding sequence as a residue.
    let vefc = tx.vefc.as_mut().expect("a coding model");
    let coding = vefc.translateable_seq.clone().expect("a coding sequence");
    let mut peptide = vep_core::codon::translate_sequence(coding.as_bytes(), 1);
    peptide.pop();
    for residue in peptide.iter_mut().filter(|residue| **residue == b'*') {
        *residue = b'U';
    }
    vefc.peptide = Some(String::from_utf8(peptide).expect("ASCII residues"));
    let variant = InputVariant::new(
        "1".into(),
        26_142_208,
        26_142_208,
        b"A".to_vec(),
        b"C".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter557SerextTer35")
    );
}

/// A new stop of the codon that the peptide reads as a residue: the
/// shape of the row at GRCh37 22:19864628 on ENST00000400519, a reverse-strand
/// transcript whose coding sequence holds a stop codon at residue 522 that its
/// peptide reads as a residue. A substitution of CDS base 1572, in the closing
/// stop codon (residue 524), loses the stop, and the first stop codon of the 3'
/// UTR, 31 residues on, is that same codon. Whether it ends the protein there
/// is not known from the model, and the count stays unknown.
#[test]
fn stop_lost_where_the_new_stop_is_a_codon_read_as_a_residue_keeps_the_unknown_count() {
    let model = RowModel {
        chr: "22",
        strand: Strand::Reverse,
        exons: &[
            (19_929_224, 19_929_326),
            (19_918_553, 19_918_618),
            (19_907_066, 19_907_122),
            (19_906_383, 19_906_527),
            (19_905_666, 19_905_740),
            (19_903_288, 19_903_366),
            (19_902_737, 19_902_799),
            (19_898_900, 19_898_970),
            (19_886_572, 19_886_591),
            (19_885_562, 19_885_653),
            (19_882_930, 19_883_104),
            (19_870_848, 19_870_984),
            (19_868_145, 19_868_240),
            (19_867_702, 19_867_794),
            (19_865_889, 19_865_960),
            (19_865_613, 19_865_710),
            (19_863_070, 19_864_757),
        ],
        cdna_coding_start: 1,
        cdna_coding_end: 1_572,
        start_phase: 0,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (19_864_535, "TCAGAGGGTCTTACTAGTACCATGGTAAAGTTCCTCGAGTCCTAATACTCTGGCTTGCACCAGTCCACTAGCCCGATATTGACACGAGATTAGTTAATTTCAGGCATCCACTCCTCTGCC"),
            (19_929_324, "CAT"),
        ],
    );
    let mut tx = transcript_of(&model, &fasta);
    // The model's peptide reads the stop codon inside the coding sequence as a residue.
    let vefc = tx.vefc.as_mut().expect("a coding model");
    let coding = vefc.translateable_seq.clone().expect("a coding sequence");
    let mut peptide = vep_core::codon::translate_sequence(coding.as_bytes(), 1);
    peptide.pop();
    for residue in peptide.iter_mut().filter(|residue| **residue == b'*') {
        *residue = b'U';
    }
    vefc.peptide = Some(String::from_utf8(peptide).expect("ASCII residues"));
    let variant = InputVariant::new(
        "22".into(),
        19_864_628,
        19_864_628,
        b"T".to_vec(),
        b"G".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter524TyrextTer?")
    );
}

/// The shape of the row at GRCh37 7:151292534 on ENST00000488258, a
/// reverse-strand transcript: the three bases of the stop codon (residue 230),
/// CDS bases 688 to 690, deleted. The whole protein is kept and 11 residues
/// follow it up to a stop of the 3' UTR. The codon window writes the stop as
/// deleted, with the extension's count and no residue in the stop's place; the
/// extension is named at the stop with the first new residue.
#[test]
fn stop_codon_deleted_before_residues_of_the_utr_is_an_extension_with_its_first_residue() {
    let model = RowModel {
        chr: "7",
        strand: Strand::Reverse,
        exons: &[
            (151_573_592, 151_573_933),
            (151_483_556, 151_483_627),
            (151_478_238, 151_478_517),
            (151_372_506, 151_372_723),
            (151_292_431, 151_292_540),
            (151_273_457, 151_273_538),
            (151_271_977, 151_272_035),
            (151_269_750, 151_269_795),
            (151_267_257, 151_267_311),
            (151_265_780, 151_265_928),
        ],
        cdna_coding_start: 229,
        cdna_coding_end: 918,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (151_292_499, "TTAAATAGCCTGTTGCTTCGAATATGCCTTTCTCCCTCAACC"),
            (151_372_506, "TATAGCGTTGAGAATCGGATT"),
            (151_573_703, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "7".into(),
        151_292_535,
        151_292_537,
        b"TCA".to_vec(),
        b"-".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter230GlyextTer11")
    );
}

/// The shape of the row at GRCh37 20:43251710 on ENST00000539235, a
/// reverse-strand transcript: CDS bases 225 to 227, the last base of the last
/// residue's codon and the first two of the stop codon (residue 76), deleted.
/// The last residue is read as before, the stop is lost and 17 residues follow
/// the whole protein up to a stop of the 3' UTR: an extension named at the stop
/// with the first new residue, where the codon window writes the stop as
/// deleted with the extension's count.
#[test]
fn deletion_from_the_last_codon_into_the_stop_codon_that_keeps_the_last_residue_is_an_extension() {
    let model = RowModel {
        chr: "20",
        strand: Strand::Reverse,
        exons: &[
            (43_280_216, 43_280_332),
            (43_264_868, 43_264_929),
            (43_257_688, 43_257_810),
            (43_251_648, 43_251_719),
            (43_251_470, 43_251_571),
            (43_251_229, 43_251_293),
            (43_249_659, 43_249_788),
            (43_248_940, 43_249_042),
            (43_248_451, 43_248_488),
        ],
        cdna_coding_start: 85,
        cdna_coding_end: 312,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (
                43_251_656,
                "CTATAAGCGCGGGGCATTAATTGGCACGCTAACGCACGTGCATGCTACCATTTCTTACAACCAC",
            ),
            (43_257_688, "CATTTCAAACAGGCAAG"),
            (43_280_246, "CAT"),
        ],
    );
    let tx = transcript_of(&model, &fasta);
    let variant = InputVariant::new(
        "20".into(),
        43_251_711,
        43_251_713,
        b"TAC".to_vec(),
        b"-".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter76GluextTer17")
    );
}

/// Where the proteins are not compared: on the model of
/// ENST00000374315, whose coding sequence holds a stop codon at residue 428
/// that the transcript's peptide reads as a residue, the three bases of the
/// closing stop codon (residue 557), CDS bases 1669 to 1671, are deleted. The
/// notation is left as the codon window gives it: the stop written as deleted,
/// with an unknown count.
#[test]
fn stop_codon_deleted_where_the_peptide_reads_a_stop_codon_inside_keeps_its_notation() {
    let model = RowModel {
        chr: "1",
        strand: Strand::Forward,
        exons: &[
            (26_126_684, 26_126_904),
            (26_127_534, 26_127_651),
            (26_131_633, 26_131_766),
            (26_135_071, 26_135_280),
            (26_135_517, 26_135_641),
            (26_136_174, 26_136_311),
            (26_137_945, 26_138_026),
            (26_138_182, 26_138_370),
            (26_139_178, 26_139_283),
            (26_140_372, 26_140_484),
            (26_140_568, 26_140_669),
            (26_142_039, 26_144_712),
        ],
        cdna_coding_start: 39,
        cdna_coding_end: 1_709,
        start_phase: -1,
    };
    let (_dir, fasta) = reference_of(
        &model,
        &[
            (26_126_722, "ATG"),
            (26_139_271, "AACATCCGCTGAT"),
            (26_140_372, "GCTTATGG"),
            (
                26_142_183,
                "AGAGAAAAGTTCTTAGAGCACTTATAAATAGGGAAGTCCAACGACATTTCAGCGGAGAGCACATGA",
            ),
        ],
    );
    let mut tx = transcript_of(&model, &fasta);
    // The model's peptide reads the stop codon inside the coding sequence as a residue.
    let vefc = tx.vefc.as_mut().expect("a coding model");
    let coding = vefc.translateable_seq.clone().expect("a coding sequence");
    let mut peptide = vep_core::codon::translate_sequence(coding.as_bytes(), 1);
    peptide.pop();
    for residue in peptide.iter_mut().filter(|residue| **residue == b'*') {
        *residue = b'U';
    }
    vefc.peptide = Some(String::from_utf8(peptide).expect("ASCII residues"));
    let variant = InputVariant::new(
        "1".into(),
        26_142_207,
        26_142_209,
        b"TAA".to_vec(),
        b"-".to_vec(),
    );
    let hgvs = generate_hgvs(&variant, &tx, Some(&fasta), HgvsOptions::default());
    assert_eq!(hgvs.offset, None);
    assert_eq!(
        hgvs.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ter557delextTer?")
    );
}

/// A made-up transcript on chromosome 11, as the transcript reads it: `utr5`,
/// then `coding`, then `utr3`. It is one exon, or two where
/// `utr3_before_the_intron` gives the count of 3' UTR bases its first exon
/// holds before an intron. `start_phase` is the count of bases the first
/// codon lacks, `0` where the coding sequence opens on a codon.
pub(crate) struct MadeUp<'a> {
    pub(crate) strand: Strand,
    pub(crate) utr5: &'a str,
    pub(crate) coding: &'a str,
    pub(crate) utr3: &'a str,
    pub(crate) utr3_before_the_intron: Option<usize>,
    pub(crate) start_phase: i8,
}

impl MadeUp<'_> {
    /// The transcript over a reference that holds its bases, a 100-base
    /// intron and 60 bases after its last base, the same on either strand as
    /// the transcript reads them.
    pub(crate) fn laid(&self) -> Laid {
        const FIRST_BASE_HELD: u64 = 63_001;
        let reverse = self.strand == Strand::Reverse;
        let in_first_exon = self.utr3_before_the_intron.unwrap_or(self.utr3.len());
        let mut read = format!(
            "{}{}{}",
            self.utr5,
            self.coding,
            &self.utr3[..in_first_exon]
        );
        let mut exons_read = vec![(0, read.len())];
        if self.utr3_before_the_intron.is_some() {
            read.push_str(&"GTAC".repeat(25));
            let second = read.len();
            read.push_str(&self.utr3[in_first_exon..]);
            exons_read.push((second, read.len()));
        }
        read.push_str(&"CTGA".repeat(15));
        let length = read.len() as u64;
        let genomic = |offset: usize| {
            if reverse {
                FIRST_BASE_HELD + length - 1 - offset as u64
            } else {
                FIRST_BASE_HELD + offset as u64
            }
        };
        let exons: Vec<(u64, u64)> = exons_read
            .iter()
            .map(|&(from, to)| {
                let (a, b) = (genomic(from), genomic(to - 1));
                (a.min(b), a.max(b))
            })
            .collect();
        let held = if reverse {
            reverse_complement(read.as_bytes())
        } else {
            read.into_bytes()
        };
        let held = String::from_utf8(held).expect("ASCII bases");
        let model = RowModel {
            chr: "11",
            strand: self.strand,
            exons: &exons,
            cdna_coding_start: self.utr5.len() as u64 + 1,
            cdna_coding_end: (self.utr5.len() + self.coding.len()) as u64,
            start_phase: if self.start_phase > 0 {
                self.start_phase
            } else {
                -1
            },
        };
        let (dir, fasta) = reference_of(&model, &[(FIRST_BASE_HELD, &held)]);
        let mut transcript = transcript_of(&model, &fasta);
        transcript.flags = Vec::<String>::new().into();
        Laid {
            _dir: Arc::new(dir),
            fasta: Arc::new(fasta),
            transcript,
            reverse,
            first_coding_base: genomic(self.utr5.len()),
        }
    }
}

/// A [`MadeUp`] transcript and its reference. A place on it is counted along
/// the genome in the transcript's direction from the first coding base, which
/// is `1`: `0` and below lie before the coding sequence, and a place past the
/// last exon base lies in the intron or after the transcript.
pub(crate) struct Laid {
    _dir: Arc<tempfile::TempDir>,
    fasta: Arc<vep_fasta::IndexedFasta>,
    /// The transcript as [`transcript_of`] gives it, without a flag: codon
    /// table 1 and no cached translation.
    pub(crate) transcript: Transcript,
    reverse: bool,
    first_coding_base: u64,
}

impl Laid {
    fn genomic(&self, place: i64) -> u64 {
        let first = self.first_coding_base as i64;
        let at = if self.reverse {
            first - (place - 1)
        } else {
            first + (place - 1)
        };
        at as u64
    }

    fn on_forward_strand(&self, bases: &str) -> Vec<u8> {
        if self.reverse {
            reverse_complement(bases.as_bytes())
        } else {
            bases.as_bytes().to_vec()
        }
    }

    /// The reference sequence the transcript lies on.
    pub(crate) fn reference(&self) -> &vep_fasta::IndexedFasta {
        &self.fasta
    }

    /// The row of `variant` on the transcript, with the reference sequence.
    pub(crate) fn row(&self, variant: &InputVariant) -> TranscriptConsequence {
        let config = EffectsConfig {
            reference_fasta: Some(Arc::clone(&self.fasta)),
            ..EffectsConfig::default()
        };
        calculate_consequences(variant, &Arc::new(self.transcript.clone()), &config)
            .expect("a row of the transcript")
    }

    /// The row of `variant` on the transcript with its HGVS notations.
    pub(crate) fn notated_row(&self, variant: &InputVariant) -> TranscriptConsequence {
        let config = EffectsConfig {
            reference_fasta: Some(Arc::clone(&self.fasta)),
            compute_hgvs: true,
            ..EffectsConfig::default()
        };
        calculate_consequences(variant, &Arc::new(self.transcript.clone()), &config)
            .expect("a row of the transcript")
    }

    /// `transcript` laid at the same places over the same reference.
    fn with_transcript(&self, transcript: Transcript) -> Laid {
        Laid {
            _dir: Arc::clone(&self._dir),
            fasta: Arc::clone(&self.fasta),
            transcript,
            reverse: self.reverse,
            first_coding_base: self.first_coding_base,
        }
    }

    /// The bases at the places `from..=to` replaced by `bases`, which the
    /// transcript reads; no bases is a deletion.
    pub(crate) fn replacement(&self, from: i64, to: i64, bases: &str) -> InputVariant {
        let (a, b) = (self.genomic(from), self.genomic(to));
        let (start, end) = (a.min(b), a.max(b));
        let reference = self
            .fasta
            .sequence("11", start, end)
            .expect("the replaced bases");
        let alternate = if bases.is_empty() {
            b"-".to_vec()
        } else {
            self.on_forward_strand(bases)
        };
        InputVariant::new("11".into(), start, end, reference, alternate)
    }

    /// The bases at the places `from..=to` deleted.
    pub(crate) fn deletion(&self, from: i64, to: i64) -> InputVariant {
        self.replacement(from, to, "")
    }

    /// `bases`, which the transcript reads, inserted after `place`.
    pub(crate) fn insertion_after(&self, place: i64, bases: &str) -> InputVariant {
        let (a, b) = (self.genomic(place), self.genomic(place + 1));
        InputVariant::new(
            "11".into(),
            a.max(b),
            a.min(b),
            b"-".to_vec(),
            self.on_forward_strand(bases),
        )
    }
}

/// The 5' UTR of the made-up transcripts of the tests below.
fn made_up_utr5() -> String {
    "GCA".repeat(10)
}

/// A 3' UTR of 18 bases that holds no stop codon in any frame.
fn made_up_utr3() -> String {
    ["GGC", "TTC", "AGG", "CAC", "GTC", "ACG"].concat()
}

/// A coding sequence that opens with the codons of Met, Ala and Arg and goes
/// on with `rest`: a fourth codon at coding bases 10 to 12 and the bases
/// after it.
fn coding_of(rest: &[&str]) -> String {
    let mut coding = ["ATG", "GCC", "CGT"].concat();
    for bases in rest {
        coding.push_str(bases);
    }
    coding
}

/// The made-up transcript of one exon that reads [`made_up_utr5`], `coding`
/// and `utr3`.
fn coding_end(strand: Strand, coding: &str, utr3: &str) -> Laid {
    MadeUp {
        strand,
        utr5: &made_up_utr5(),
        coding,
        utr3,
        utr3_before_the_intron: None,
        start_phase: 0,
    }
    .laid()
}

/// The same transcript with an intron after the first
/// `utr3_before_the_intron` bases of its 3' UTR.
fn coding_end_before_an_intron(
    strand: Strand,
    coding: &str,
    utr3: &str,
    utr3_before_the_intron: usize,
) -> Laid {
    MadeUp {
        strand,
        utr5: &made_up_utr5(),
        coding,
        utr3,
        utr3_before_the_intron: Some(utr3_before_the_intron),
        start_phase: 0,
    }
    .laid()
}

/// Asserts that the row of `variant` on `laid` carries exactly `expected`, in
/// any order.
pub(crate) fn assert_row(laid: &Laid, variant: &InputVariant, expected: &[&str], case: &str) {
    let row = laid.row(variant);
    let mut terms: Vec<String> = row.consequences.iter().map(|c| c.to_string()).collect();
    terms.sort();
    let mut expected: Vec<String> = expected.iter().map(|term| (*term).to_string()).collect();
    expected.sort();
    assert_eq!(terms, expected, "{case}");
}

/// The coding ends at which no stop codon stands, each as the bases after the
/// opening three codons and whether a 3' UTR follows them. The last whole
/// codon names Pro. A tail of one or two bases opens a codon that the 3' UTR
/// goes on with or that ends the transcript, where adenines after it would
/// spell the codon of Lys.
const ENDS_WITHOUT_A_STOP_CODON: [(&str, &[&str], bool); 6] = [
    ("whole codons before a 3' UTR", &["CCA"], true),
    ("whole codons at the transcript's end", &["CCA"], false),
    ("a tail of one base before a 3' UTR", &["CCA", "A"], true),
    ("a tail of two bases before a 3' UTR", &["CCA", "AA"], true),
    (
        "a tail of one base at the transcript's end",
        &["CCA", "A"],
        false,
    ),
    (
        "a tail of two bases at the transcript's end",
        &["CCA", "AA"],
        false,
    ),
];

/// The transcript of one row of [`ENDS_WITHOUT_A_STOP_CODON`] and its coding
/// length.
fn end_without_a_stop_codon(strand: Strand, rest: &[&str], utr3: bool) -> (Laid, i64) {
    let coding = coding_of(rest);
    let utr3 = if utr3 { made_up_utr3() } else { String::new() };
    (coding_end(strand, &coding, &utr3), coding.len() as i64)
}

/// The shape of the row at GRCh38 11:111908773 on ENST00000526167: one base
/// deleted from the last whole codon of a coding sequence that ends one or
/// two bases into the next codon. The last three coding positions are no
/// codon of the frame and no stop codon, so the deletion shifts the frame and
/// loses no stop codon. Ensembl VEP reads those three positions as the stop
/// codon on every transcript without `cds_end_NF` (VariationEffect.pm:1373-1381)
/// and writes `stop_lost` beside the frameshift, for a sequence left shorter
/// (:1469) or for three bases left there that are no stop codon (:1473-1480).
#[test]
fn deletion_in_the_last_whole_codon_of_a_coding_sequence_that_ends_inside_a_codon_loses_no_stop() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (end, rest, utr3) in ENDS_WITHOUT_A_STOP_CODON {
            let tail = rest.len() - 1;
            if tail == 0 {
                continue;
            }
            let (laid, length) = end_without_a_stop_codon(strand, rest, utr3);
            // The last three coding positions that lie in the last whole codon.
            for deleted in length - 2..=12 {
                let case = format!("{strand:?}, {end}, coding base {deleted} deleted");
                assert_row(
                    &laid,
                    &laid.deletion(deleted, deleted),
                    &["frameshift_variant"],
                    &case,
                );
            }
        }
    }
}

/// The same deletion on a coding sequence of whole codons whose last codon
/// names a residue: no stop codon stands at its last three positions.
#[test]
fn deletion_in_the_last_codon_of_a_coding_sequence_whose_last_codon_names_a_residue_loses_no_stop()
{
    for strand in [Strand::Forward, Strand::Reverse] {
        for (end, rest, utr3) in ENDS_WITHOUT_A_STOP_CODON {
            if rest.len() > 1 {
                continue;
            }
            let (laid, _) = end_without_a_stop_codon(strand, rest, utr3);
            for deleted in [10, 11, 12] {
                let case = format!("{strand:?}, {end}, coding base {deleted} deleted");
                assert_row(
                    &laid,
                    &laid.deletion(deleted, deleted),
                    &["frameshift_variant"],
                    &case,
                );
            }
        }
    }
}

/// One base inserted between two of the last three coding positions, inside
/// the last whole codon, where no stop codon stands: a frameshift alone.
/// Ensembl VEP writes `stop_lost` beside it, for the three bases the insertion
/// leaves at the old coding end (VariationEffect.pm:1473-1480).
#[test]
fn insertion_in_the_last_whole_codon_of_a_coding_sequence_no_stop_codon_closes_loses_no_stop() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (end, rest, utr3) in ENDS_WITHOUT_A_STOP_CODON {
            let (laid, length) = end_without_a_stop_codon(strand, rest, utr3);
            // The places between two of the last three coding positions that
            // lie inside the last whole codon: none behind a tail of two.
            for after in length - 2..12 {
                let case = format!("{strand:?}, {end}, a base inserted after coding base {after}");
                assert_row(
                    &laid,
                    &laid.insertion_after(after, "G"),
                    &["frameshift_variant"],
                    &case,
                );
            }
        }
    }
}

/// A one-base deletion in the last whole codon after which the three bases
/// at the old coding end are a stop codon made with the first base of the 3'
/// UTR. No stop codon stood there, so none is kept: the deletion is a
/// frameshift. Ensembl VEP reads a stop codon kept (VariationEffect.pm:1318)
/// and withholds the frameshift behind it (:1549).
#[test]
fn deletion_that_leaves_a_stop_codon_of_the_3_prime_utr_at_the_old_coding_end_keeps_no_stop() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (end, rest, first_utr_base) in [
            ("whole codons, the last of Cys", &["TGC"][..], "A"),
            ("a tail of one base", &["CTG", "A"][..], "G"),
            ("a tail of two bases", &["CCA", "TA"][..], "G"),
        ] {
            let utr3 = format!("{first_utr_base}{}", made_up_utr3());
            let laid = coding_end(strand, &coding_of(rest), &utr3);
            let case = format!("{strand:?}, {end}, coding base 12 deleted");
            assert_row(
                &laid,
                &laid.deletion(12, 12),
                &["frameshift_variant"],
                &case,
            );
        }
    }
}

/// The three bases of the last codon deleted where the first three bases of
/// the 3' UTR are a stop codon: whole codons are removed and no stop codon
/// stood at the coding end, so the row is `inframe_deletion` alone. Ensembl
/// VEP writes `stop_retained_variant` beside it for the 3' UTR bases that
/// move to the old coding end (VariationEffect.pm:1318).
#[test]
fn deletion_of_the_last_codon_before_a_stop_codon_of_the_3_prime_utr_is_an_inframe_deletion_alone()
{
    for strand in [Strand::Forward, Strand::Reverse] {
        let utr3 = format!("TAA{}", made_up_utr3());
        let laid = coding_end(strand, &coding_of(&["CCA"]), &utr3);
        let case = format!("{strand:?}");
        assert_row(&laid, &laid.deletion(10, 12), &["inframe_deletion"], &case);
    }
}

/// A deletion from inside a coding sequence that no stop codon closes into
/// the 3' UTR, both ends in the exon: it removes coding bases and 3' UTR bases
/// and no stop codon, whether or not the 3' UTR bases left at the old coding
/// end are a stop codon. Ensembl VEP writes `stop_lost` or, where they are,
/// `stop_retained_variant` (VariationEffect.pm:1263, 1378-1381, 1466-1480,
/// 1318).
#[test]
fn deletion_from_a_coding_sequence_no_stop_codon_closes_into_the_3_prime_utr_is_a_coding_sequence_variant(
) {
    let expected = ["coding_sequence_variant", "3_prime_UTR_variant"];
    for strand in [Strand::Forward, Strand::Reverse] {
        // The 3' UTR holds a stop codon at its second to fourth base.
        let utr3 = format!("GTGA{}", made_up_utr3());
        for (end, rest, no_stop_codon_left, stop_codon_left) in [
            ("whole codons", &["CCA"][..], (11, 14), (10, 13)),
            ("a tail of one base", &["CCA", "A"][..], (12, 14), (11, 14)),
            (
                "a tail of two bases",
                &["CCA", "AA"][..],
                (11, 15),
                (12, 15),
            ),
        ] {
            let laid = coding_end(strand, &coding_of(rest), &utr3);
            for (left, (from, to)) in [
                ("no stop codon", no_stop_codon_left),
                ("a stop codon", stop_codon_left),
            ] {
                let case = format!(
                    "{strand:?}, {end}, coding bases {from} to {to} deleted, {left} left at the old end"
                );
                assert_row(&laid, &laid.deletion(from, to), &expected, &case);
            }
        }
    }
}

/// A deletion from inside a coding sequence that no stop codon closes past
/// the transcript's last base, over the last three coding positions: no stop
/// codon is removed. Ensembl VEP writes `stop_retained_variant` for such a
/// span, whose stop codon edit test it cannot make (VariationEffect.pm:1508
/// read at :1318).
#[test]
fn deletion_from_a_coding_sequence_no_stop_codon_closes_past_the_transcripts_end_is_a_coding_sequence_variant(
) {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (end, rest, utr3) in ENDS_WITHOUT_A_STOP_CODON {
            let (laid, length) = end_without_a_stop_codon(strand, rest, utr3);
            let utr3_length = if utr3 { made_up_utr3().len() as i64 } else { 0 };
            let (from, to) = (length - 3, length + utr3_length + 10);
            let case = format!("{strand:?}, {end}, coding base {from} to 10 bases past the end");
            assert_row(
                &laid,
                &laid.deletion(from, to),
                &["coding_sequence_variant", "3_prime_UTR_variant"],
                &case,
            );
        }
    }
}

/// The same deletion with its 3' end in an intron, which follows the coding
/// end at once or after six bases of the 3' UTR.
#[test]
fn deletion_from_a_coding_sequence_no_stop_codon_closes_into_an_intron_is_a_coding_sequence_variant(
) {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (end, rest, utr3_before_the_intron) in [
            ("whole codons, six bases of 3' UTR", &["CCA"][..], 6),
            ("whole codons at the exon's end", &["CCA"][..], 0),
            (
                "a tail of two bases, six bases of 3' UTR",
                &["CCA", "AA"][..],
                6,
            ),
            ("a tail of one base at the exon's end", &["CCA", "A"][..], 0),
        ] {
            let coding = coding_of(rest);
            let laid = coding_end_before_an_intron(
                strand,
                &coding,
                &made_up_utr3(),
                utr3_before_the_intron,
            );
            let length = coding.len() as i64;
            let (from, to) = (length - 3, length + utr3_before_the_intron as i64 + 10);
            let case = format!("{strand:?}, {end}, coding base {from} to 10 bases into the intron");
            assert_row(
                &laid,
                &laid.deletion(from, to),
                &[
                    "coding_sequence_variant",
                    "3_prime_UTR_variant",
                    "splice_donor_variant",
                    "splice_donor_5th_base_variant",
                    "intron_variant",
                ],
                &case,
            );
        }
    }
}

/// One base inserted just before the first of the last three coding
/// positions that leaves a stop codon at the old coding end, where none
/// stood: a frameshift on both strands. Ensembl VEP's genomic overlap reaches
/// the insertion on the reverse strand alone (VariationEffect.pm:1408-1427),
/// where it writes `stop_retained_variant` beside an in-frame term.
#[test]
fn insertion_before_the_last_coding_bases_that_leaves_a_stop_codon_there_is_a_frameshift() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for (end, rest, utr3, after) in [
            ("whole codons before a 3' UTR", &["AAC"][..], true, 9),
            (
                "whole codons at the transcript's end",
                &["AAC"][..],
                false,
                9,
            ),
            (
                "a tail of one base before a 3' UTR",
                &["CAG", "C"][..],
                true,
                10,
            ),
            (
                "a tail of two bases before a 3' UTR",
                &["CCA", "AC"][..],
                true,
                11,
            ),
            (
                "a tail of two bases at the transcript's end",
                &["CCA", "AC"][..],
                false,
                11,
            ),
        ] {
            let (laid, _) = end_without_a_stop_codon(strand, rest, utr3);
            let case = format!("{strand:?}, {end}, a base inserted after coding base {after}");
            assert_row(
                &laid,
                &laid.insertion_after(after, "T"),
                &["frameshift_variant"],
                &case,
            );
        }
    }
}

/// Three bases inserted in the last whole codon that keep its residue and
/// bring a stop codon behind it, on a coding sequence no stop codon closes:
/// a stop codon is gained where none stood. Ensembl VEP reads the protein
/// unchanged up to its old length with a stop after it as a stop codon kept
/// (VariationEffect.pm:1357) and withholds `stop_gained` behind that (:1221),
/// with `cds_end_NF` as without it.
#[test]
fn insertion_that_brings_a_stop_codon_behind_the_last_residue_gains_a_stop_where_none_stood() {
    let expected = ["stop_gained", "inframe_insertion"];
    for strand in [Strand::Forward, Strand::Reverse] {
        for (end, rest, utr3) in [
            ("whole codons before a 3' UTR", &["GCA"][..], true),
            ("whole codons at the transcript's end", &["GCA"][..], false),
            (
                "a tail of one base before a 3' UTR",
                &["GCA", "A"][..],
                true,
            ),
            (
                "a tail of two bases before a 3' UTR",
                &["GCA", "AA"][..],
                true,
            ),
            (
                "a tail of two bases at the transcript's end",
                &["GCA", "AA"][..],
                false,
            ),
        ] {
            let (mut laid, _) = end_without_a_stop_codon(strand, rest, utr3);
            let inserted = laid.insertion_after(11, "ATA");
            assert_row(&laid, &inserted, &expected, &format!("{strand:?}, {end}"));
            laid.transcript.flags = vec!["cds_end_NF".to_string()].into();
            let case = format!("{strand:?}, {end}, with cds_end_NF");
            assert_row(&laid, &inserted, &expected, &case);
        }
    }
}

/// Where no stop codon stands, a substitution in the last whole codon is read
/// from its codons (VariationEffect.pm:1082-1098, 1228), and an allele that
/// starts in the bases after the last whole codon is read as lying in an
/// incomplete terminal codon (:1601-1610), in Ensembl VEP and here.
#[test]
fn substitution_in_the_last_whole_codon_and_change_in_the_bases_after_it_read_as_ensembl_vep_reads_them(
) {
    let partial = [
        "incomplete_terminal_codon_variant",
        "coding_sequence_variant",
    ];
    for strand in [Strand::Forward, Strand::Reverse] {
        let case = format!("{strand:?}");
        // The last whole codon names Tyr.
        let laid = coding_end(strand, &coding_of(&["TAC"]), &made_up_utr3());
        assert_row(
            &laid,
            &laid.replacement(11, 11, "T"),
            &["missense_variant"],
            &case,
        );
        assert_row(
            &laid,
            &laid.replacement(12, 12, "T"),
            &["synonymous_variant"],
            &case,
        );
        assert_row(
            &laid,
            &laid.replacement(12, 12, "A"),
            &["stop_gained"],
            &case,
        );

        let laid = coding_end(strand, &coding_of(&["CCA", "AA"]), &made_up_utr3());
        assert_row(&laid, &laid.replacement(13, 13, "C"), &partial, &case);
        assert_row(&laid, &laid.deletion(14, 14), &partial, &case);
        assert_row(&laid, &laid.deletion(13, 14), &partial, &case);
        assert_row(
            &laid,
            &laid.insertion_after(12, "G"),
            &[
                "inframe_insertion",
                "incomplete_terminal_codon_variant",
                "coding_sequence_variant",
            ],
            &case,
        );
        assert_row(
            &laid,
            &laid.deletion(14, 16),
            &[
                "incomplete_terminal_codon_variant",
                "coding_sequence_variant",
                "3_prime_UTR_variant",
            ],
            &case,
        );
    }
}

/// A one-base deletion or insertion in a stop codon that closes the coding
/// sequence keeps its stop term: `stop_lost` where the three bases left at
/// the stop codon's place are no stop codon, `stop_retained_variant` where
/// the first base of the 3' UTR makes one with the two bases left.
#[test]
fn one_base_change_in_a_stop_codon_that_closes_the_coding_sequence_keeps_its_stop_term() {
    for strand in [Strand::Forward, Strand::Reverse] {
        let case = format!("{strand:?}");
        let laid = coding_end(strand, &coding_of(&["CCA", "TAA"]), &made_up_utr3());
        assert_row(&laid, &laid.deletion(13, 13), &["stop_lost"], &case);
        assert_row(
            &laid,
            &laid.deletion(14, 14),
            &["stop_retained_variant"],
            &case,
        );
        assert_row(
            &laid,
            &laid.deletion(15, 15),
            &["stop_retained_variant"],
            &case,
        );
        assert_row(&laid, &laid.insertion_after(14, "C"), &["stop_lost"], &case);
        assert_row(
            &laid,
            &laid.insertion_after(13, "G"),
            &["stop_retained_variant"],
            &case,
        );

        // The transcript ends on its stop codon.
        let laid = coding_end(strand, &coding_of(&["CCA", "TAA"]), "");
        assert_row(&laid, &laid.deletion(13, 13), &["stop_lost"], &case);
    }
}

/// The three bases of a stop codon that closes the coding sequence deleted
/// where the next three bases of the 3' UTR are a stop codon: `stop_lost`,
/// from the peptides (VariationEffect.pm:1259-1260).
#[test]
fn deletion_of_a_stop_codon_that_closes_the_coding_sequence_is_stop_lost() {
    for strand in [Strand::Forward, Strand::Reverse] {
        let utr3 = format!("TAG{}", made_up_utr3());
        let laid = coding_end(strand, &coding_of(&["CCA", "TAA"]), &utr3);
        assert_row(
            &laid,
            &laid.deletion(13, 15),
            &["stop_lost"],
            &format!("{strand:?}"),
        );
    }
}

/// A deletion from a stop codon that closes the coding sequence into the 3'
/// UTR keeps its stop term: `stop_retained_variant` where the 3' UTR bases
/// left at the stop codon's place are a stop codon, `stop_lost` where they
/// are not.
#[test]
fn deletion_from_a_stop_codon_that_closes_the_coding_sequence_into_the_3_prime_utr_keeps_its_stop_term(
) {
    for strand in [Strand::Forward, Strand::Reverse] {
        let case = format!("{strand:?}");
        let utr3 = format!("GTGA{}", made_up_utr3());
        let laid = coding_end(strand, &coding_of(&["CCA", "TAA"]), &utr3);
        assert_row(
            &laid,
            &laid.deletion(14, 17),
            &["stop_retained_variant", "3_prime_UTR_variant"],
            &case,
        );
        assert_row(
            &laid,
            &laid.deletion(14, 16),
            &["stop_lost", "3_prime_UTR_variant"],
            &case,
        );
    }
}

/// A deletion from before a stop codon that closes the coding sequence past
/// the transcript's last base, or into an intron after six bases of the 3'
/// UTR, removes the stop codon: `stop_lost`.
#[test]
fn deletion_over_a_stop_codon_that_closes_the_coding_sequence_past_the_exons_end_is_stop_lost() {
    for strand in [Strand::Forward, Strand::Reverse] {
        let case = format!("{strand:?}");
        let coding = coding_of(&["CCA", "TAA"]);
        for utr3 in [made_up_utr3(), String::new()] {
            let laid = coding_end(strand, &coding, &utr3);
            let to = (coding.len() + utr3.len() + 10) as i64;
            assert_row(
                &laid,
                &laid.deletion(12, to),
                &["stop_lost", "3_prime_UTR_variant"],
                &format!("{case}, a 3' UTR of {} bases", utr3.len()),
            );
        }
        let laid = coding_end_before_an_intron(strand, &coding, &made_up_utr3(), 6);
        assert_row(
            &laid,
            &laid.deletion(12, coding.len() as i64 + 6 + 10),
            &[
                "stop_lost",
                "3_prime_UTR_variant",
                "splice_donor_variant",
                "splice_donor_5th_base_variant",
                "intron_variant",
            ],
            &format!("{case}, into the intron"),
        );
    }
}

/// Four bases inserted just before a stop codon that closes the coding
/// sequence, the first three of them a stop codon: a stop codon stands at the
/// old place and the row keeps `stop_retained_variant`.
#[test]
fn insertion_before_a_stop_codon_that_closes_the_coding_sequence_that_leaves_a_stop_codon_there_is_stop_retained(
) {
    for strand in [Strand::Forward, Strand::Reverse] {
        let laid = coding_end(strand, &coding_of(&["CCA", "TAA"]), &made_up_utr3());
        assert_row(
            &laid,
            &laid.insertion_after(12, "TGAC"),
            &["inframe_insertion", "stop_retained_variant"],
            &format!("{strand:?}"),
        );
    }
}

/// Three bases inserted in the codon of the last residue that keep the
/// residue and bring a stop codon behind it, before a stop codon that closes
/// the coding sequence: the stop codon is kept.
#[test]
fn insertion_that_brings_a_stop_codon_behind_the_last_residue_before_a_stop_codon_is_stop_retained()
{
    for strand in [Strand::Forward, Strand::Reverse] {
        let laid = coding_end(strand, &coding_of(&["GCA", "TAA"]), &made_up_utr3());
        assert_row(
            &laid,
            &laid.insertion_after(11, "ATA"),
            &["inframe_insertion", "stop_retained_variant"],
            &format!("{strand:?}"),
        );
    }
}

/// A coding sequence that ends one or two bases into a codon at the
/// transcript's last base, where those bases followed by adenines are a stop
/// codon, on a record that does not show the model complete: the standard
/// code and no consensus coding sequence id. No stop codon stands there, and
/// no row takes a stop term from the coding end: a deletion or insertion
/// inside the last whole codon, a deletion that runs from it into the bases
/// after it, and a deletion past the transcript's end, where Ensembl VEP
/// reads the last three coding positions as the stop codon
/// (VariationEffect.pm:1373-1381, 1418-1427). An allele that starts in the
/// bases after the last whole codon lies in an incomplete terminal codon, as
/// Ensembl VEP reads it (:1601-1610).
#[test]
fn change_in_the_last_whole_codon_takes_no_stop_term_where_the_record_does_not_show_a_completed_stop_codon(
) {
    let frameshift = ["frameshift_variant"];
    let partial = [
        "incomplete_terminal_codon_variant",
        "coding_sequence_variant",
    ];
    for strand in [Strand::Forward, Strand::Reverse] {
        let case = format!("{strand:?}, a tail of one base");
        let laid = coding_end(strand, &coding_of(&["CCA", "T"]), "");
        // Inside the last whole codon.
        assert_row(&laid, &laid.deletion(11, 11), &frameshift, &case);
        assert_row(&laid, &laid.deletion(12, 12), &frameshift, &case);
        assert_row(&laid, &laid.deletion(10, 11), &frameshift, &case);
        assert_row(&laid, &laid.deletion(11, 12), &frameshift, &case);
        assert_row(&laid, &laid.insertion_after(11, "G"), &frameshift, &case);
        assert_row(
            &laid,
            &laid.replacement(11, 11, "T"),
            &["missense_variant"],
            &case,
        );
        // From the last whole codon into the base after it.
        assert_row(&laid, &laid.deletion(12, 13), &frameshift, &case);
        assert_row(&laid, &laid.deletion(11, 13), &["inframe_deletion"], &case);
        // The base after the last whole codon.
        assert_row(&laid, &laid.replacement(13, 13, "C"), &partial, &case);
        assert_row(
            &laid,
            &laid.deletion(13, 13),
            &[
                "incomplete_terminal_codon_variant",
                "protein_altering_variant",
            ],
            &case,
        );
        // Past the transcript's end over the last three coding positions.
        assert_row(
            &laid,
            &laid.deletion(10, 23),
            &["coding_sequence_variant", "3_prime_UTR_variant"],
            &case,
        );
        // Three bases that keep the last residue and bring a stop codon
        // behind it gain a stop codon where none stood.
        let laid = coding_end(strand, &coding_of(&["GCA", "T"]), "");
        assert_row(
            &laid,
            &laid.insertion_after(11, "ATA"),
            &["inframe_insertion", "stop_gained"],
            &case,
        );
        // One base inside the last whole codon that leaves a stop codon at
        // the last three coding positions.
        let laid = coding_end(strand, &coding_of(&["CAG", "T"]), "");
        assert_row(&laid, &laid.insertion_after(10, "T"), &frameshift, &case);

        let case = format!("{strand:?}, a tail of two bases");
        let laid = coding_end(strand, &coding_of(&["CCA", "TA"]), "");
        // Inside the last whole codon.
        assert_row(&laid, &laid.deletion(12, 12), &frameshift, &case);
        assert_row(&laid, &laid.deletion(11, 12), &frameshift, &case);
        assert_row(&laid, &laid.insertion_after(11, "TGAC"), &frameshift, &case);
        // From the last whole codon into the bases after it.
        assert_row(&laid, &laid.deletion(12, 13), &frameshift, &case);
        assert_row(&laid, &laid.deletion(10, 13), &frameshift, &case);
        assert_row(&laid, &laid.deletion(12, 14), &["inframe_deletion"], &case);
        // The bases after the last whole codon.
        assert_row(&laid, &laid.replacement(14, 14, "C"), &partial, &case);
        assert_row(&laid, &laid.deletion(13, 13), &partial, &case);
        assert_row(
            &laid,
            &laid.deletion(13, 14),
            &[
                "incomplete_terminal_codon_variant",
                "protein_altering_variant",
            ],
            &case,
        );
        // Past the transcript's end over the last three coding positions.
        assert_row(
            &laid,
            &laid.deletion(11, 24),
            &["coding_sequence_variant", "3_prime_UTR_variant"],
            &case,
        );
    }
}

/// What in a record shows that polyadenylation completes the stop codon its
/// coding sequence ends inside.
#[derive(Clone, Copy, Debug)]
enum ShownBy {
    /// A mitochondrial codon table: 2, the vertebrate code, or 5, the
    /// invertebrate code.
    CodonTable(u8),
    /// A consensus coding sequence id, on a transcript of the standard code.
    ConsensusId,
}

/// The records each test below lays its coding ends on.
const RECORDS_THAT_SHOW_THE_STOP_CODON: [ShownBy; 3] = [
    ShownBy::CodonTable(2),
    ShownBy::CodonTable(5),
    ShownBy::ConsensusId,
];

/// The made-up transcript that reads [`coding_of`] `rest` and ends on its
/// last coding base, on a record that shows the model complete.
fn end_that_adenines_complete(strand: Strand, rest: &[&str], shown_by: ShownBy) -> Laid {
    shown(coding_end(strand, &coding_of(rest), ""), shown_by)
}

/// `laid` with what `shown_by` names on its record.
fn shown(mut laid: Laid, shown_by: ShownBy) -> Laid {
    match shown_by {
        ShownBy::CodonTable(table) => {
            laid.transcript
                .vefc
                .as_mut()
                .expect("a coding model")
                .codon_table = table;
        }
        ShownBy::ConsensusId => laid.transcript.ccds = Some("CCDS1.1".into()),
    }
    laid
}

/// The transcript of [`end_that_adenines_complete`] on a genome that holds
/// what polyadenylation adds: its coding sequence ends on the completed stop
/// codon and a 3' UTR of 30 adenines follows, each a base of the reference.
fn end_with_the_adenines_on_the_genome(strand: Strand, rest: &[&str], shown_by: ShownBy) -> Laid {
    let mut coding = coding_of(rest);
    let adenines = 3 - coding.len() % 3;
    coding.push_str(&"A".repeat(adenines));
    shown(coding_end(strand, &coding, &"A".repeat(30)), shown_by)
}

/// The count of adenines the made-up records below carry after the stop
/// codon as their 3' UTR sequence: more than a made-up transcript has bases.
const POLY_A_TAIL_IN_THE_RECORD: usize = 90;

/// The protein notation of `row` after its `p.`.
pub(crate) fn protein_notation(row: &TranscriptConsequence) -> Option<&str> {
    let written = row.hgvsp.as_deref()?;
    written
        .rsplit_once(":p.")
        .map(|(_, description)| description)
}

/// The terms of `row` in the order it carries them.
pub(crate) fn terms_of(row: &TranscriptConsequence) -> Vec<String> {
    row.consequences.iter().map(|c| c.to_string()).collect()
}

/// A change at a made-up coding end, by the places it takes, so that one
/// change is laid on each record it is read on.
#[derive(Clone, Copy, Debug)]
enum ChangeAtTheEnd {
    Replaced(i64, &'static str),
    Deleted(i64, i64),
    InsertedAfter(i64, &'static str),
}

impl ChangeAtTheEnd {
    fn on(self, laid: &Laid) -> InputVariant {
        match self {
            Self::Replaced(place, bases) => laid.replacement(place, place, bases),
            Self::Deleted(from, to) => laid.deletion(from, to),
            Self::InsertedAfter(place, bases) => laid.insertion_after(place, bases),
        }
    }
}

/// The transcript of `laid` as a record that carries its stop codon whole:
/// the adenines that complete the codon follow the coding sequence, which
/// ends that many cDNA positions later, and no mapper pair holds them, as on
/// a model whose `_rna_edit` appends them. `poly_a_tail` more adenines stand
/// after the stop codon as the record's 3' UTR sequence, without a mapper
/// pair either; with none the record ends on its stop codon.
fn with_the_stop_codon_in_the_record(laid: &Laid, poly_a_tail: usize) -> Laid {
    let mut transcript = laid.transcript.clone();
    let vefc = transcript.vefc.as_mut().expect("a coding model");
    let coding = vefc.translateable_seq.as_mut().expect("a coding sequence");
    let adenines = 3 - coding.len() % 3;
    coding.push_str(&"A".repeat(adenines));
    vefc.three_prime_utr = (poly_a_tail > 0).then(|| "A".repeat(poly_a_tail));
    vefc.mapper.as_mut().expect("a mapper").cdna_coding_end += adenines as u64;
    transcript.cdna_coding_end = transcript.cdna_coding_end.map(|end| end + adenines as u64);
    laid.with_transcript(transcript)
}

/// A substitution in the one or two bases of a stop codon that
/// polyadenylation completes, on a record that shows the model complete. The
/// codon is those bases and the adenines after them, residue 5 of the made-up
/// protein. A change that leaves a codon naming a residue loses the stop
/// codon, and the protein goes on with no stop codon in sight; a change that
/// leaves a stop codon of the table keeps it. The codon cell shows the
/// adenines as unchanged bases, and the positions are the record's own.
/// Ensembl VEP reads an allele that starts in those bases as lying in an
/// incomplete terminal codon (VariationEffect.pm:1601-1610).
#[test]
fn substitution_in_a_stop_codon_that_adenines_complete_loses_or_keeps_the_stop_codon() {
    use ShownBy::{CodonTable, ConsensusId};
    // The record, the bases after the last whole codon, the places replaced
    // and the bases put there; then the term, `Amino_acids`, `Codons` and
    // the protein notation.
    let rows = [
        // One base: the codon of Gln under every code.
        (
            CodonTable(2),
            "T",
            (13, 13),
            "C",
            "stop_lost",
            "*/Q",
            "Taa/Caa",
            "Ter5GlnextTer?",
        ),
        (
            CodonTable(5),
            "T",
            (13, 13),
            "C",
            "stop_lost",
            "*/Q",
            "Taa/Caa",
            "Ter5GlnextTer?",
        ),
        (
            ConsensusId,
            "T",
            (13, 13),
            "C",
            "stop_lost",
            "*/Q",
            "Taa/Caa",
            "Ter5GlnextTer?",
        ),
        // The second of two bases: the codon of Trp under the mitochondrial
        // codes, a stop codon of the standard code.
        (
            CodonTable(2),
            "TA",
            (14, 14),
            "G",
            "stop_lost",
            "*/W",
            "tAa/tGa",
            "Ter5TrpextTer?",
        ),
        (
            CodonTable(5),
            "TA",
            (14, 14),
            "G",
            "stop_lost",
            "*/W",
            "tAa/tGa",
            "Ter5TrpextTer?",
        ),
        (
            ConsensusId,
            "TA",
            (14, 14),
            "G",
            "stop_retained_variant",
            "*",
            "tAa/tGa",
            "Ter5=",
        ),
        // Both bases: a stop codon of the vertebrate mitochondrial code, the
        // codon of Ser under the invertebrate one and of Arg under the
        // standard code.
        (
            CodonTable(2),
            "TA",
            (13, 14),
            "AG",
            "stop_retained_variant",
            "*",
            "TAa/AGa",
            "Ter5=",
        ),
        (
            CodonTable(5),
            "TA",
            (13, 14),
            "AG",
            "stop_lost",
            "*/S",
            "TAa/AGa",
            "Ter5SerextTer?",
        ),
        (
            ConsensusId,
            "TA",
            (13, 14),
            "AG",
            "stop_lost",
            "*/R",
            "TAa/AGa",
            "Ter5ArgextTer?",
        ),
    ];
    for strand in [Strand::Forward, Strand::Reverse] {
        for (shown_by, tail, (from, to), bases, term, amino_acids, codons, notation) in rows {
            let case = format!("{strand:?}, {shown_by:?}, places {from} to {to}");
            let laid = end_that_adenines_complete(strand, &["CCA", tail], shown_by);
            let row = laid.notated_row(&laid.replacement(from, to, bases));
            let terms: Vec<String> = row.consequences.iter().map(|c| c.to_string()).collect();
            assert_eq!(terms, [term], "{case}");
            assert_eq!(row.amino_acids.as_deref(), Some(amino_acids), "{case}");
            assert_eq!(row.codons.as_deref(), Some(codons), "{case}");
            let written = row.hgvsp.as_deref().and_then(|p| p.rsplit_once(":p."));
            assert_eq!(
                written.map(|(_, description)| description),
                Some(notation),
                "{case}"
            );
            let places = if from == to {
                from.to_string()
            } else {
                format!("{from}-{to}")
            };
            assert_eq!(row.cds_position.as_deref(), Some(places.as_str()), "{case}");
            assert_eq!(row.protein_position.as_deref(), Some("5"), "{case}");
        }
    }
}

/// A deletion of one base or of both bases of a stop codon that
/// polyadenylation completes. The adenines that follow spell, with what is
/// left, the codon of Lys or of Glu: the stop codon is lost. Ensembl VEP
/// reads the allele as lying in an incomplete terminal codon
/// (VariationEffect.pm:1601-1610).
#[test]
fn deletion_in_a_stop_codon_that_adenines_complete_is_stop_lost() {
    let lost = ["stop_lost"];
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            let case = format!("{strand:?}, {shown_by:?}");
            let laid = end_that_adenines_complete(strand, &["CCA", "T"], shown_by);
            assert_row(&laid, &laid.deletion(13, 13), &lost, &case);
            let laid = end_that_adenines_complete(strand, &["CCA", "TA"], shown_by);
            assert_row(&laid, &laid.deletion(13, 13), &lost, &case);
            assert_row(&laid, &laid.deletion(13, 14), &lost, &case);
        }
        // The other stop codon of the vertebrate mitochondrial code that two
        // bases begin.
        let case = format!("{strand:?}, two bases of the other stop codon");
        let laid = end_that_adenines_complete(strand, &["CCA", "AG"], ShownBy::CodonTable(2));
        assert_row(&laid, &laid.deletion(13, 13), &lost, &case);
        assert_row(&laid, &laid.deletion(14, 14), &lost, &case);
        assert_row(&laid, &laid.deletion(13, 14), &lost, &case);
    }
}

/// An insertion at a stop codon that polyadenylation completes. Between the
/// last whole codon and the codon's first base, one base shifts the frame
/// and three bases add a residue before the stop codon. Between the codon's
/// two genomic bases, a base that leaves a stop codon at the old place keeps
/// it, and three bases that read a residue and then the stop codon are an
/// in-frame insertion. Ensembl VEP writes `inframe_insertion` beside the
/// incomplete terminal codon for every length (VariationEffect.pm:1100-1128).
#[test]
fn insertion_at_a_stop_codon_that_adenines_complete_is_read_on_the_whole_codon() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            for tail in ["T", "TA"] {
                let case = format!("{strand:?}, {shown_by:?}, a tail of {}", tail.len());
                let laid = end_that_adenines_complete(strand, &["CCA", tail], shown_by);
                assert_row(
                    &laid,
                    &laid.insertion_after(12, "C"),
                    &["frameshift_variant"],
                    &case,
                );
                assert_row(
                    &laid,
                    &laid.insertion_after(12, "GCT"),
                    &["inframe_insertion"],
                    &case,
                );
            }
            let case = format!("{strand:?}, {shown_by:?}, between the two bases");
            let laid = end_that_adenines_complete(strand, &["CCA", "TA"], shown_by);
            assert_row(
                &laid,
                &laid.insertion_after(13, "A"),
                &["stop_retained_variant"],
                &case,
            );
            assert_row(
                &laid,
                &laid.insertion_after(13, "GCT"),
                &["inframe_insertion"],
                &case,
            );
        }
    }
}

/// A deletion or insertion of one base inside the last whole codon before a
/// stop codon that polyadenylation completes touches no base of the stop
/// codon: a frameshift alone. Ensembl VEP reads the last three coding
/// positions as the stop codon and writes `stop_lost` beside the frameshift
/// (VariationEffect.pm:1378-1381, 1469).
#[test]
fn one_base_change_in_the_last_whole_codon_before_a_stop_codon_that_adenines_complete_is_a_frameshift_alone(
) {
    let frameshift = ["frameshift_variant"];
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            for tail in ["T", "TA"] {
                let case = format!("{strand:?}, {shown_by:?}, a tail of {}", tail.len());
                let laid = end_that_adenines_complete(strand, &["CCA", tail], shown_by);
                for deleted in 10..=12 {
                    assert_row(&laid, &laid.deletion(deleted, deleted), &frameshift, &case);
                }
                assert_row(&laid, &laid.insertion_after(10, "T"), &frameshift, &case);
                assert_row(&laid, &laid.insertion_after(11, "G"), &frameshift, &case);
            }
        }
    }
}

/// A deletion of the last bases of the last whole codon with the bases of a
/// stop codon that polyadenylation completes. Two bases shift the frame and
/// lose the stop codon. Three bases keep the frame: where the bases left and
/// the adenines spell a stop codon the row is an in-frame deletion alone,
/// and where they name a residue the stop codon is lost beside it.
#[test]
fn deletion_from_the_last_whole_codon_through_a_stop_codon_that_adenines_complete_reads_the_bases_left(
) {
    let lost_in_frame = ["stop_lost", "inframe_deletion"];
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            let case = format!("{strand:?}, {shown_by:?}");
            let laid = end_that_adenines_complete(strand, &["CCA", "T"], shown_by);
            assert_row(
                &laid,
                &laid.deletion(12, 13),
                &["frameshift_variant", "stop_lost"],
                &case,
            );
            // The base left and the adenines are the codon of Gln.
            assert_row(&laid, &laid.deletion(11, 13), &lost_in_frame, &case);
            // The two bases left and the adenine are the codon of Pro.
            let laid = end_that_adenines_complete(strand, &["CCA", "TA"], shown_by);
            assert_row(&laid, &laid.deletion(12, 14), &lost_in_frame, &case);
            // The codon of Ser, whose first base the adenines complete to a
            // stop codon.
            let laid = end_that_adenines_complete(strand, &["TCA", "T"], shown_by);
            assert_row(&laid, &laid.deletion(11, 13), &["inframe_deletion"], &case);
            // The codon of Tyr, whose first two bases the adenine completes
            // to a stop codon.
            let laid = end_that_adenines_complete(strand, &["TAC", "TA"], shown_by);
            assert_row(&laid, &laid.deletion(12, 14), &["inframe_deletion"], &case);
        }
    }
}

/// A deletion that runs past the last base of a transcript whose stop codon
/// polyadenylation completes, and holds every base the genome has of that
/// codon, removes the stop codon: `stop_lost`, whether the deletion starts
/// inside the coding sequence before the codon or on the codon's first base.
/// The `3_prime_UTR_variant` beside it is the term every span across the
/// last base of a transcript without a 3' UTR takes. Ensembl VEP writes
/// `stop_retained_variant` for the first, a stop codon edit test it does not
/// make (VariationEffect.pm:1318 with 1508), and the incomplete terminal
/// codon for the second (:1601-1610).
#[test]
fn deletion_past_the_end_of_a_transcript_whose_stop_codon_adenines_complete_is_stop_lost() {
    let lost = ["stop_lost", "3_prime_UTR_variant"];
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            for tail in ["T", "TA"] {
                let case = format!("{strand:?}, {shown_by:?}, a tail of {}", tail.len());
                let laid = end_that_adenines_complete(strand, &["CCA", tail], shown_by);
                let past_the_end = 12 + tail.len() as i64 + 10;
                for from in [10, 12, 13] {
                    assert_row(
                        &laid,
                        &laid.deletion(from, past_the_end),
                        &lost,
                        &format!("{case}, from coding base {from}"),
                    );
                }
            }
        }
    }
}

/// A deletion that starts on the second of two bases that begin a stop codon
/// polyadenylation completes, and runs past the transcript's last base,
/// leaves the codon's first base with the adenines of the poly(A) tail after
/// it: the stop codon is kept where that base and two adenines are a stop
/// codon of the table, and lost where they name a residue. A deletion of both
/// bases past the end leaves no base of the codon and is `stop_lost`, as the
/// deletion of the one base of a stop codon that one base begins, unless the
/// reference repeats the codon's first base after it: it is then the deletion
/// from the second base, written one place back. Ensembl VEP reads each as
/// lying in an incomplete terminal codon (VariationEffect.pm:1601-1610).
#[test]
fn deletion_from_the_second_base_of_a_stop_codon_that_adenines_complete_past_the_end_reads_the_base_left(
) {
    let kept: &[&str] = &["stop_retained_variant", "3_prime_UTR_variant"];
    let lost: &[&str] = &["stop_lost", "3_prime_UTR_variant"];
    // Each two-base tail with the row of a deletion from its second base.
    // The first base of the first three, with two adenines, is a stop codon
    // of the standard code and of table 2; the first base of the last,
    // whose stop codon is one of table 2 alone, names Lys with them.
    let ends = [
        (ShownBy::CodonTable(2), "an adenine second", "TA", kept),
        (ShownBy::ConsensusId, "an adenine second", "TA", kept),
        (ShownBy::ConsensusId, "another base second", "TG", kept),
        (ShownBy::CodonTable(2), "an adenine first", "AG", lost),
    ];
    for strand in [Strand::Forward, Strand::Reverse] {
        for (shown_by, tail_with, tail, from_the_second_base) in ends {
            let laid = end_that_adenines_complete(strand, &["CCA", tail], shown_by);
            let case = |deleted: &str, past_the_end: i64| {
                format!("{strand:?}, {shown_by:?}, {tail_with}, {deleted} and {past_the_end} past the end")
            };
            // The base after each of these deletions differs from the
            // deletion's first base, so none can be written further along.
            for past_the_end in [1, 4] {
                assert_row(
                    &laid,
                    &laid.deletion(14, 14 + past_the_end),
                    from_the_second_base,
                    &case("the second tail base", past_the_end),
                );
            }
            for past_the_end in [4, 10] {
                assert_row(
                    &laid,
                    &laid.deletion(13, 14 + past_the_end),
                    lost,
                    &case("both tail bases", past_the_end),
                );
            }
            // The second base past the end is the first base of the three
            // tails that open on it: both tail bases and one base more is
            // there the deletion from the second tail base.
            assert_row(
                &laid,
                &laid.deletion(13, 15),
                from_the_second_base,
                &case("both tail bases", 1),
            );
        }
    }
}

/// `begun`, whose coding sequence ends inside its stop codon, as a record
/// whose `_rna_edit` replaces the bases after the last whole codon by the
/// whole stop codon, as a cache holds a model whose alignment completes it:
/// the record ends on its stop codon and carries the edit's attribute.
fn with_the_stop_codon_from_an_edit(begun: &Laid) -> Laid {
    let mut laid = with_the_stop_codon_in_the_record(begun, 0);
    let coding = begun
        .transcript
        .vefc
        .as_ref()
        .and_then(|vefc| vefc.translateable_seq.as_deref())
        .expect("a coding sequence");
    let tail = &coding[coding.len() - coding.len() % 3..];
    let end = begun.transcript.cdna_coding_end.expect("a coding end");
    laid.transcript.attributes.push(Attribute {
        code: "_rna_edit".into(),
        value: format!(
            "{} {end} {tail}{}",
            end + 1 - tail.len() as u64,
            "A".repeat(3 - tail.len())
        ),
    });
    laid
}

/// A record whose own `_rna_edit` completes its stop codon ends on that
/// codon, with no poly(A) tail in its sequence, and a deletion past its last
/// base is read on the bases the genome has of the codon alone: one that
/// holds every such base is `stop_lost`, and one that starts on the second
/// of two takes no stop term.
#[test]
fn deletion_past_the_end_of_a_record_whose_edit_completes_its_stop_codon_reads_the_genome_bases() {
    let lost = ["stop_lost", "3_prime_UTR_variant"];
    let no_stop_term = ["coding_sequence_variant", "3_prime_UTR_variant"];
    let table_2 = ShownBy::CodonTable(2);
    for strand in [Strand::Forward, Strand::Reverse] {
        let two = with_the_stop_codon_from_an_edit(&end_that_adenines_complete(
            strand,
            &["CCA", "TA"],
            table_2,
        ));
        for past_the_end in [1, 2, 3, 4, 10] {
            let case =
                format!("{strand:?}, two stop bases on the genome, {past_the_end} past the end");
            assert_row(
                &two,
                &two.deletion(14, 14 + past_the_end),
                &no_stop_term,
                &format!("{case}, from the second"),
            );
            assert_row(
                &two,
                &two.deletion(13, 14 + past_the_end),
                &lost,
                &format!("{case}, both"),
            );
        }
        let one = with_the_stop_codon_from_an_edit(&end_that_adenines_complete(
            strand,
            &["CCA", "T"],
            table_2,
        ));
        // One base past the end the deletion can be written off the
        // transcript, which keeps the codon.
        for past_the_end in [2, 3, 4, 10] {
            let case =
                format!("{strand:?}, one stop base on the genome, {past_the_end} past the end");
            assert_row(&one, &one.deletion(13, 13 + past_the_end), &lost, &case);
        }
    }
}

/// The same deletions on a model whose record carries the stop codon whole,
/// two or one of its bases without a genomic position: the genome's bases of
/// the stop codon are removed and the row is `stop_lost`, as on the record
/// that leaves the adenines to polyadenylation.
#[test]
fn deletion_past_the_end_of_a_model_whose_record_completes_its_stop_codon_is_stop_lost() {
    let lost = ["stop_lost", "3_prime_UTR_variant"];
    for strand in [Strand::Forward, Strand::Reverse] {
        for tail in ["T", "TA"] {
            let case = format!("{strand:?}, a tail of {}", tail.len());
            let begun = coding_end(strand, &coding_of(&["CCA", tail]), "");
            let laid = with_the_stop_codon_in_the_record(&begun, 0);
            let past_the_end = 12 + tail.len() as i64 + 10;
            for from in [10, 12, 13] {
                assert_row(
                    &laid,
                    &laid.deletion(from, past_the_end),
                    &lost,
                    &format!("{case}, from coding base {from}"),
                );
            }
            // A deletion inside the transcript reads the codon as before.
            assert_row(&laid, &laid.deletion(13, 13), &["stop_lost"], &case);
        }
    }
}

/// A deletion of bases of a stop codon that polyadenylation completes is
/// read on into the poly(A) tail: what is left of the codon and the adenines
/// after it name a residue, and no stop codon follows. The stop codon is
/// lost and the notation is the extension that meets no new stop codon,
/// never the stop written as unchanged.
#[test]
fn deletion_in_a_stop_codon_that_adenines_complete_is_an_extension_that_meets_no_stop_codon() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            // The one base, the first of two and both: adenines alone are
            // left at the codon's place, the codon of Lys.
            for (tail, from, to) in [("T", 13, 13), ("TA", 13, 13), ("TA", 13, 14)] {
                let case = format!("{strand:?}, {shown_by:?}, {tail:?}, {from} to {to}");
                let laid = end_that_adenines_complete(strand, &["CCA", tail], shown_by);
                let row = laid.notated_row(&laid.deletion(from, to));
                assert_eq!(terms_of(&row), ["stop_lost"], "{case}");
                assert_eq!(protein_notation(&row), Some("Ter5LysextTer?"), "{case}");
                assert_eq!(row.protein_position.as_deref(), Some("5"), "{case}");
            }
        }
        // The other stop codon of the vertebrate mitochondrial code that two
        // bases begin: its second base and the adenines are the codon of
        // Glu, its first base and the adenines the codon of Lys.
        let laid = end_that_adenines_complete(strand, &["CCA", "AG"], ShownBy::CodonTable(2));
        for (from, to, notation) in [
            (13, 13, "Ter5GluextTer?"),
            (14, 14, "Ter5LysextTer?"),
            (13, 14, "Ter5LysextTer?"),
        ] {
            let case = format!("{strand:?}, the other stop codon, {from} to {to}");
            let row = laid.notated_row(&laid.deletion(from, to));
            assert_eq!(terms_of(&row), ["stop_lost"], "{case}");
            assert_eq!(protein_notation(&row), Some(notation), "{case}");
        }
    }
}

/// A deletion of the second base of a stop codon that polyadenylation
/// completes, where that base is an adenine: the adenine after it takes its
/// place and the message reads the same stop codon. The stop codon is kept
/// and the notation writes it as unchanged.
#[test]
fn deletion_of_an_adenine_of_a_stop_codon_that_adenines_complete_keeps_the_stop_codon() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            let case = format!("{strand:?}, {shown_by:?}");
            let laid = end_that_adenines_complete(strand, &["CCA", "TA"], shown_by);
            let row = laid.notated_row(&laid.deletion(14, 14));
            assert_eq!(terms_of(&row), ["stop_retained_variant"], "{case}");
            assert_eq!(protein_notation(&row), Some("Ter5="), "{case}");
            assert_eq!(row.cds_position.as_deref(), Some("14"), "{case}");
            assert_eq!(row.protein_position.as_deref(), Some("5"), "{case}");
        }
    }
}

/// An adenine inserted where the adenines of a stop codon that
/// polyadenylation completes stand leaves the message as it was. After the
/// transcript's last base it lies in the poly(A) tail, outside the
/// transcript, and the row is the one an insertion after the last base of
/// any transcript takes, with no coding term. Before the adenine that is the
/// second base of the codon it lies inside the transcript and keeps the stop
/// codon.
#[test]
fn adenine_inserted_among_the_adenines_of_a_stop_codon_that_adenines_complete_changes_no_codon() {
    let after_the_transcript = ["downstream_gene_variant"];
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            let case = format!("{strand:?}, {shown_by:?}");
            let laid = end_that_adenines_complete(strand, &["CCA", "T"], shown_by);
            let inserted = laid.insertion_after(13, "A");
            assert_row(&laid, &inserted, &after_the_transcript, &case);
            let laid = end_that_adenines_complete(strand, &["CCA", "TA"], shown_by);
            let inserted = laid.insertion_after(14, "A");
            assert_row(&laid, &inserted, &after_the_transcript, &case);
            let row = laid.notated_row(&laid.insertion_after(13, "A"));
            assert_eq!(terms_of(&row), ["stop_retained_variant"], "{case}");
        }
    }
}

/// A change before a stop codon that polyadenylation completes whose new
/// reading frame runs through the codon's place is notated from the poly(A)
/// tail too: the residues after the change are read from the adenines, and
/// no stop codon follows them.
#[test]
fn frame_that_runs_through_a_stop_codon_that_adenines_complete_reads_the_poly_a_tail() {
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            let case = format!("{strand:?}, {shown_by:?}");
            let laid = end_that_adenines_complete(strand, &["CCA", "T"], shown_by);
            // The last base of the last whole codon: the codon still names
            // Pro, and the codon of Lys stands at the stop codon's place.
            let row = laid.notated_row(&laid.deletion(12, 12));
            assert_eq!(terms_of(&row), ["frameshift_variant"], "{case}");
            assert_eq!(protein_notation(&row), Some("Ter5LysextTer?"), "{case}");
            // That base and the stop codon's own.
            let row = laid.notated_row(&laid.deletion(12, 13));
            assert_eq!(
                terms_of(&row),
                ["frameshift_variant", "stop_lost"],
                "{case}"
            );
            assert_eq!(protein_notation(&row), Some("Ter5LysextTer?"), "{case}");
            // Every coding base after the first codon: the codon of Lys
            // follows Met.
            let row = laid.notated_row(&laid.deletion(4, 13));
            assert_eq!(
                terms_of(&row),
                ["frameshift_variant", "stop_lost"],
                "{case}"
            );
            assert_eq!(protein_notation(&row), Some("Ala2LysfsTer?"), "{case}");
            // Three bases through both bases of the stop codon: the last
            // whole codon is completed to the codon of Pro and the codon of
            // Lys follows.
            let laid = end_that_adenines_complete(strand, &["CCA", "TA"], shown_by);
            let row = laid.notated_row(&laid.deletion(12, 14));
            assert_eq!(terms_of(&row), ["stop_lost", "inframe_deletion"], "{case}");
            assert_eq!(protein_notation(&row), Some("Ter5LysextTer?"), "{case}");
        }
    }
}

/// Every cell of a row on a record that shows its stop codon completed is
/// the cell of the same change on a record that carries the completed codon
/// and, after it, the poly(A) tail as its 3' UTR sequence: the terms, the
/// positions, the residues, the codons and both notations, over each
/// substitution, deletion and insertion at the last codons and each deletion
/// past the transcript's end.
///
/// A change inside the transcript reads the same on a genome that holds the
/// adenines, where the completed stop codon and a 3' UTR of adenines are
/// bases of the reference. That transcript's coordinates add three things
/// the record without them does not have: its cDNA length; the place the 3'
/// rule gives a notation, which moves an allele along the genome's adenines
/// there and, on the record that ends on its coding sequence, along the
/// bases after the transcript; and the insertion after the last base, which
/// lies inside the transcript only there. So the notations are held to each
/// other where the 3' rule moves neither allele.
#[test]
fn row_on_a_record_that_shows_its_stop_codon_completed_is_the_row_of_the_completed_record() {
    use ChangeAtTheEnd::{Deleted, InsertedAfter, Replaced};
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            for tail in ["T", "TA", "TG"] {
                let begun = end_that_adenines_complete(strand, &["CCA", tail], shown_by);
                if crate::coding::stop_codon_at_the_coding_end(&begun.transcript)
                    == crate::coding::StopCodonAtTheCodingEnd::None
                {
                    continue;
                }
                let in_the_record =
                    with_the_stop_codon_in_the_record(&begun, POLY_A_TAIL_IN_THE_RECORD);
                let on_the_genome =
                    end_with_the_adenines_on_the_genome(strand, &["CCA", tail], shown_by);
                let last = 12 + tail.len() as i64;
                let mut inside: Vec<ChangeAtTheEnd> = Vec::new();
                let mut past_the_end: Vec<ChangeAtTheEnd> = Vec::new();
                for place in 1..=last {
                    for bases in ["A", "C", "G", "T"] {
                        let change = Replaced(place, bases).on(&begun);
                        if change.ref_allele != change.alt_alleles[0] {
                            inside.push(Replaced(place, bases));
                        }
                    }
                    for to in place..=last {
                        inside.push(Deleted(place, to));
                    }
                    past_the_end.push(Deleted(place, last + 10));
                }
                for after in 1..last {
                    for bases in ["A", "C", "G", "T", "GCT", "TAA", "TGAC"] {
                        inside.push(InsertedAfter(after, bases));
                    }
                }
                for change in inside.iter().chain(&past_the_end) {
                    let case = format!(
                        "{strand:?}, {shown_by:?}, a tail of {}, {change:?}",
                        tail.len()
                    );
                    let laid = change.on(&begun);
                    let (row, expected) =
                        (begun.notated_row(&laid), in_the_record.notated_row(&laid));
                    assert_eq!(row.consequences, expected.consequences, "{case}");
                    assert_eq!(row.impact, expected.impact, "{case}");
                    assert_eq!(row.cdna_position, expected.cdna_position, "{case}");
                    assert_eq!(row.cds_position, expected.cds_position, "{case}");
                    assert_eq!(row.protein_position, expected.protein_position, "{case}");
                    assert_eq!(row.amino_acids, expected.amino_acids, "{case}");
                    assert_eq!(row.codons, expected.codons, "{case}");
                    assert_eq!(row.hgvsc, expected.hgvsc, "{case}");
                    assert_eq!(row.hgvsp, expected.hgvsp, "{case}");
                    assert_eq!(row.hgvs_offset, expected.hgvs_offset, "{case}");
                }
                for change in &inside {
                    let case = format!(
                        "{strand:?}, {shown_by:?}, a tail of {}, {change:?}, on the genome",
                        tail.len()
                    );
                    let row = begun.notated_row(&change.on(&begun));
                    let expected = on_the_genome.notated_row(&change.on(&on_the_genome));
                    assert_eq!(row.consequences, expected.consequences, "{case}");
                    assert_eq!(row.impact, expected.impact, "{case}");
                    assert_eq!(row.cdna_position, expected.cdna_position, "{case}");
                    assert_eq!(row.cds_position, expected.cds_position, "{case}");
                    assert_eq!(row.protein_position, expected.protein_position, "{case}");
                    assert_eq!(row.amino_acids, expected.amino_acids, "{case}");
                    assert_eq!(row.codons, expected.codons, "{case}");
                    if row.hgvs_offset.is_none() && expected.hgvs_offset.is_none() {
                        assert_eq!(row.hgvsc, expected.hgvsc, "{case}");
                        assert_eq!(row.hgvsp, expected.hgvsp, "{case}");
                    }
                }
            }
        }
    }
}

/// A row on a record that shows its stop codon completed belongs to that
/// record: it holds the transcript it was asked for, and the lengths beside
/// its positions are the record's own, the coding sequence without the
/// adenines.
#[test]
fn row_on_a_record_that_shows_its_stop_codon_completed_keeps_the_records_own_lengths() {
    for strand in [Strand::Forward, Strand::Reverse] {
        let laid = end_that_adenines_complete(strand, &["CCA", "T"], ShownBy::CodonTable(2));
        let transcript = Arc::new(laid.transcript.clone());
        let config = EffectsConfig {
            reference_fasta: Some(Arc::clone(&laid.fasta)),
            ..EffectsConfig::default()
        };
        let row = calculate_consequences(&laid.replacement(13, 13, "C"), &transcript, &config)
            .expect("a row of the transcript");
        assert_eq!(row.consequences.len(), 1, "{strand:?}");
        assert_eq!(row.consequences[0].to_string(), "stop_lost", "{strand:?}");
        assert!(
            row.transcript
                .as_ref()
                .is_some_and(|held| Arc::ptr_eq(held, &transcript)),
            "{strand:?}"
        );
        assert_eq!(
            (row.cdna_total, row.cds_total, row.protein_total),
            (Some(43), Some(13), None),
            "{strand:?}"
        );
    }
}

/// A structural variant is read on the record as the cache holds it, where
/// its own arm asks the one rule: a structural deletion that takes the base
/// of a stop codon that polyadenylation completes loses the stop codon,
/// beside the terms of a deletion that runs past the transcript's end. A
/// structural deletion of the last base of the last whole codon takes no
/// base of the stop codon. Ensembl VEP reads the first span as starting in
/// an incomplete terminal codon, with no lost stop codon
/// (VariationEffect.pm:1240), and holds the second against the last three
/// coding positions (:1268-1280).
#[test]
fn structural_deletion_that_takes_a_stop_codon_that_adenines_complete_is_stop_lost() {
    let structural = |laid: &Laid, from: i64, to: i64| {
        let (a, b) = (laid.genomic(from), laid.genomic(to));
        let (start, end) = (a.min(b), a.max(b));
        let mut deletion = InputVariant::new("11".into(), start, end, b"N".to_vec(), b"-".to_vec());
        deletion.variant_class = vep_core::variant::VariantClass::StructuralDeletion;
        deletion.is_structural = true;
        deletion.sv_end = Some(end);
        deletion
    };
    for strand in [Strand::Forward, Strand::Reverse] {
        for shown_by in RECORDS_THAT_SHOW_THE_STOP_CODON {
            let case = format!("{strand:?}, {shown_by:?}");
            let laid = end_that_adenines_complete(strand, &["CCA", "T"], shown_by);
            let deletion = structural(&laid, 13, 23);
            let transcript = Arc::new(laid.transcript.clone());
            let row = calculate_consequences(&deletion, &transcript, &EffectsConfig::default())
                .expect("a row of the transcript");
            let cached = crate::sv::calculate_sv_consequences(&deletion, &transcript, 5000, 5000)
                .expect("a row of the transcript");
            assert_eq!(row.consequences, cached.consequences, "{case}");
            assert_row(
                &laid,
                &deletion,
                &[
                    "feature_truncation",
                    "stop_lost",
                    "coding_sequence_variant",
                    "3_prime_UTR_variant",
                ],
                &case,
            );
            assert_row(
                &laid,
                &structural(&laid, 12, 12),
                &["feature_truncation", "frameshift_variant"],
                &case,
            );
        }
    }
}

/// What the one rule of the stop terms reads at each made-up coding end
/// (`stop_codon_at_the_coding_end`).
#[test]
fn stop_codon_at_a_coding_end_is_whole_begun_by_its_last_bases_or_none() {
    use crate::coding::{stop_codon_at_the_coding_end, StopCodonAtTheCodingEnd};
    let (none, whole) = (
        StopCodonAtTheCodingEnd::None,
        StopCodonAtTheCodingEnd::Whole,
    );
    let utr3 = made_up_utr3();
    for strand in [Strand::Forward, Strand::Reverse] {
        let at = |rest: &[&str], utr3: &str| {
            stop_codon_at_the_coding_end(&coding_end(strand, &coding_of(rest), utr3).transcript)
        };
        for (end, rest, has_utr3) in ENDS_WITHOUT_A_STOP_CODON {
            let (laid, _) = end_without_a_stop_codon(strand, rest, has_utr3);
            assert_eq!(
                stop_codon_at_the_coding_end(&laid.transcript),
                none,
                "{end}"
            );
        }
        // Whole codons, the last a stop codon, with a 3' UTR and without one.
        assert_eq!(at(&["CCA", "TAA"], &utr3), whole);
        assert_eq!(at(&["CCA", "TAA"], ""), whole);
        // Another stop codon inside the coding sequence leaves the last one
        // closing it.
        assert_eq!(at(&["TGA", "CCA", "TAA"], &utr3), whole);
        // One or two bases that adenines complete to a stop codon, on the
        // transcript's last base: the standard code and no consensus coding
        // sequence id do not show the model complete, and no stop codon
        // stands.
        let begun = |bases| StopCodonAtTheCodingEnd::BegunByTheTail { bases };
        assert_eq!(at(&["CCA", "T"], ""), none);
        assert_eq!(at(&["CCA", "TA"], ""), none);
        assert_eq!(at(&["CCA", "TG"], ""), none);
        // A consensus coding sequence id shows it; an empty one and the
        // placeholder of a record without one do not.
        for (id, shown) in [("CCDS1.1", true), ("", false), ("-", false)] {
            let mut laid = coding_end(strand, &coding_of(&["CCA", "T"]), "");
            laid.transcript.ccds = Some(id.into());
            assert_eq!(
                stop_codon_at_the_coding_end(&laid.transcript),
                if shown { begun(1) } else { none },
                "id {id:?}"
            );
        }
        // The vertebrate, invertebrate, echinoderm and ascidian mitochondrial
        // codes show it; the mold and protozoan one, under which the same
        // base begins a stop codon, does not.
        for (table, shown) in [(2, true), (5, true), (9, true), (13, true), (4, false)] {
            let mut laid = coding_end(strand, &coding_of(&["CCA", "T"]), "");
            laid.transcript
                .vefc
                .as_mut()
                .expect("a coding model")
                .codon_table = table;
            assert_eq!(
                stop_codon_at_the_coding_end(&laid.transcript),
                if shown { begun(1) } else { none },
                "table {table}"
            );
        }
        // The same bases before a 3' UTR open a codon the 3' UTR goes on with.
        assert_eq!(at(&["CCA", "T"], &utr3), none);
        assert_eq!(at(&["CCA", "TA"], &utr3), none);

        let mut flagged = coding_end(strand, &coding_of(&["CCA", "TAA"]), &utr3);
        flagged.transcript.flags = vec!["cds_end_NF".to_string()].into();
        assert_eq!(stop_codon_at_the_coding_end(&flagged.transcript), none);

        // A last codon, and two last bases that adenines complete to it, that
        // are a stop codon of table 2 and name Arg under table 1.
        for (rest, utr3, under_table_2) in [
            (&["CCA", "AGA"][..], utr3.as_str(), whole),
            (&["CCA", "AG"][..], "", begun(2)),
        ] {
            let mut laid = coding_end(strand, &coding_of(rest), utr3);
            assert_eq!(stop_codon_at_the_coding_end(&laid.transcript), none);
            laid.transcript
                .vefc
                .as_mut()
                .expect("a coding model")
                .codon_table = 2;
            assert_eq!(
                stop_codon_at_the_coding_end(&laid.transcript),
                under_table_2
            );
        }
    }
}

/// `replicate_ins_del_stop_altered` makes its test only where a stop codon
/// stands at the coding end: a deletion from the last codon into the 3' UTR
/// is not evaluated on a coding sequence whose last codon names a residue,
/// and is on one that a stop codon closes, where the bases left at the stop
/// codon's place decide.
#[test]
fn stop_codon_edit_test_of_a_deletion_into_the_3_prime_utr_is_made_only_where_a_stop_codon_stands()
{
    use crate::coding::replicate_ins_del_stop_altered;
    // The 3' UTR holds a stop codon at its second to fourth base.
    let utr3 = format!("GTGA{}", made_up_utr3());
    for strand in [Strand::Forward, Strand::Reverse] {
        let altered = |laid: &Laid, from: i64, to: i64| {
            replicate_ins_del_stop_altered(
                &laid.deletion(from, to),
                &laid.transcript,
                from as u64,
                (to - from + 1) as usize,
                Some(laid.reference()),
            )
        };
        let open = coding_end(strand, &coding_of(&["CCA"]), &utr3);
        assert_eq!(altered(&open, 11, 14), None, "{strand:?}");
        assert_eq!(altered(&open, 10, 13), None, "{strand:?}");
        let closed = coding_end(strand, &coding_of(&["CCA", "TAA"]), &utr3);
        assert_eq!(altered(&closed, 14, 16), Some(true), "{strand:?}");
        assert_eq!(altered(&closed, 14, 17), Some(false), "{strand:?}");
    }
}
