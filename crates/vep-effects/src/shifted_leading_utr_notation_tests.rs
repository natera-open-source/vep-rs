// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The protein notation of an indel written before the coding start, in the
//! 5' UTR or in an intron of it, and shifted into the coding sequence.
//!
//! Each test lays a made-up 5' UTR, coding sequence and 3' UTR on the
//! three-exon test transcript, on either strand, and holds the row's HGVSc,
//! its offset and its HGVSp to values derived from those bases, beside the row
//! of the same change given where the 3' shift carries it. Ensembl VEP writes
//! the HGVSc and the offset of each row and no HGVSp: `hgvs_protein` tests the
//! `coding` predicate of the allele as written
//! (TranscriptVariationAllele.pm:1731), which a span that does not overlap the
//! coding region never gets (BaseVariationFeatureOverlapAllele.pm:523-525).
//! The last test holds no row: it holds the test of such a span,
//! `coding::annotated_before_the_coding_start`, at its bounds on both strands.

use std::sync::Arc;

use crate::concordance_tests::{
    bases_from, compact_on_chr1, hgvsc_fasta_fixture, make_compact_test_transcript,
    make_descending_reverse_strand_transcript,
};
use crate::consequences::{calculate_consequences, EffectsConfig};
use crate::hgvs::{generate_hgvs, HgvsNotation, HgvsOptions};
use crate::test_helpers::make_test_transcript;
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

/// `tx` with the coding sequence `cds` from cDNA `first`, the 5' UTR `utr5`
/// before it and the 3' UTR `utr3` after it: the mapper, the transcript's
/// coding fields and the three sequences.
fn coding(mut tx: Transcript, first: u64, utr5: &str, cds: &str, utr3: &str) -> Transcript {
    assert_eq!(utr5.len() as u64, first - 1);
    let last = first + cds.len() as u64 - 1;
    let (genomic_first, genomic_last) = (genomic_of(&tx, first), genomic_of(&tx, last));
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.five_prime_utr = Some(utr5.to_string());
    vefc.translateable_seq = Some(cds.to_string());
    vefc.three_prime_utr = Some(utr3.to_string());
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

/// A start codon and 82 alanine codons, then `tail`.
fn alanines_then(tail: &str) -> String {
    format!("ATG{}{tail}", "GCT".repeat(82))
}

/// A 5' UTR of `length` bases that ends in `tail`.
fn utr_ending(length: usize, tail: &str) -> String {
    let filler: String = "GCCACC".chars().cycle().take(length - tail.len()).collect();
    format!("{filler}{tail}")
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

/// An insertion of `bases` between `position - 1` and `position`.
fn insertion(position: u64, bases: &[u8]) -> InputVariant {
    InputVariant::new(
        "1".into(),
        position,
        position - 1,
        b"-".to_vec(),
        bases.to_vec(),
    )
}

/// The Sequence Ontology terms of the row, sorted.
fn terms(variant: &InputVariant, tx: &Transcript) -> Vec<&'static str> {
    let row = calculate_consequences(variant, &Arc::new(tx.clone()), &EffectsConfig::default())
        .expect("a row on the transcript");
    let mut terms: Vec<_> = row.consequences.iter().map(|term| term.so_term()).collect();
    terms.sort_unstable();
    terms
}

/// The forward transcript's 5' UTR is 100..149 and ends `AC`; coding base 1
/// is 150 and the coding sequence opens `ATG GCT GCT`. Seven bases written
/// between 147 and 148, in the 5' UTR, repeat the seven after them (`AC` and
/// coding bases 1 to 5), so the 3' shift carries them seven bases on, to stand
/// after coding base 5: a duplication from two bases before the coding start
/// to coding base 5. The frame then reads `ATG GCA CAT GGC TGC ...`: residue 2
/// is kept, residue 3 is a histidine, and the frame runs through the coding
/// sequence and into the 3' UTR `GCTAGGGTT`, where `AGC TAG` end it 86
/// residues on. The row's term stays that of the 5' UTR.
#[test]
fn a_duplication_written_in_the_5_prime_utr_is_named_as_the_frameshift_at_its_3_prime_place() {
    let utr5 = utr_ending(50, "AC");
    let cds = alanines_then("TAA");
    let tx = coding(make_compact_test_transcript(), 51, &utr5, &cds, "GCTAGGGTT");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(100, &format!("{utr5}{}", &cds[..250])));
    let written = insertion(148, b"ACATGGC");
    let (row, alone) = rows(&tx, &fasta, &written, &insertion(155, b"ACATGGC"));
    assert_eq!(terms(&written, &tx), ["5_prime_UTR_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.-2_5dup"));
    assert_eq!(row.offset, Some(7));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala3HisfsTer86")
    );
    assert_eq!(alone.offset, None);
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// The same on the reverse strand, where the shift runs towards lower
/// coordinates: exon 1 is 5801..6100, the 5' UTR 6100..6051 and coding base 1
/// 6050, so the forward strand reads, upwards from 5801, the reverse
/// complement of the exon's bases. The seven bases are written between 6053
/// and 6052 and stand, shifted, between 6045 and 6046.
#[test]
fn a_duplication_written_in_the_5_prime_utr_of_a_reverse_strand_transcript_is_the_same_frameshift()
{
    let utr5 = utr_ending(50, "AC");
    let cds = alanines_then("TAA");
    let tx = coding(
        compact_on_chr1(make_descending_reverse_strand_transcript()),
        51,
        &utr5,
        &cds,
        "GCTAGGGTT",
    );
    let exon = format!("{utr5}{}", &cds[..250]);
    let forward_strand = String::from_utf8(vep_core::codon::reverse_complement(exon.as_bytes()))
        .expect("ASCII bases");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(5801, &forward_strand));
    let written = insertion(6053, b"GCCATGT");
    let (row, alone) = rows(&tx, &fasta, &written, &insertion(6046, b"GCCATGT"));
    assert_eq!(terms(&written, &tx), ["5_prime_UTR_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.-2_5dup"));
    assert_eq!(row.offset, Some(-7));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala3HisfsTer86")
    );
    assert_eq!(alone.offset, None);
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// The 5' UTR ends in three `A` behind a `C`, and the start codon's `A` is a
/// fourth: deleting the first of the run, at 147, is deleting coding base 1,
/// three bases on. The row carries the protein notation of the deletion given
/// at coding base 1, whatever that notation is: the start codon stands again
/// on the 5' UTR's last base, and the notation of such a change is the start
/// codon's own subject.
#[test]
fn a_deletion_in_a_run_that_ends_on_coding_base_1_takes_the_notation_of_the_deleted_start_base() {
    let utr5 = utr_ending(50, "CAAA");
    let cds = alanines_then("TAA");
    let tx = coding(make_compact_test_transcript(), 51, &utr5, &cds, "GCTAGGGTT");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(100, &format!("{utr5}{}", &cds[..250])));
    let written = deletion(147, b"A");
    let (row, alone) = rows(&tx, &fasta, &written, &deletion(150, b"A"));
    assert_eq!(terms(&written, &tx), ["5_prime_UTR_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.1del"));
    assert_eq!(row.offset, Some(3));
    assert_eq!(alone.offset, None);
    assert!(alone.hgvsp.is_some());
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// The 5' UTR ends in a copy of the start codon: deleting that copy, at
/// 147..149, is deleting coding bases 1 to 3, three bases on, and the row
/// carries the notation of the start codon deleted where it stands.
#[test]
fn a_start_codon_copy_deleted_in_the_5_prime_utr_takes_the_notation_of_the_deleted_start_codon() {
    let utr5 = utr_ending(50, "CCATG");
    let cds = alanines_then("TAA");
    let tx = coding(make_compact_test_transcript(), 51, &utr5, &cds, "GCTAGGGTT");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(100, &format!("{utr5}{}", &cds[..250])));
    let written = InputVariant::new("1".into(), 147, 149, b"ATG".to_vec(), b"-".to_vec());
    let there = InputVariant::new("1".into(), 150, 152, b"ATG".to_vec(), b"-".to_vec());
    let (row, alone) = rows(&tx, &fasta, &written, &there);
    assert_eq!(terms(&written, &tx), ["5_prime_UTR_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.1_3del"));
    assert_eq!(row.offset, Some(3));
    assert_eq!(alone.offset, None);
    assert!(alone.hgvsp.is_some());
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// An allele written in an intron of the 5' UTR. The coding sequence starts
/// at cDNA 303, so exon 2 opens at 2100 with the 5' UTR bases `CC` and coding
/// base 1 is 2102; intron 1 ends in a `G` at 2099. Eight bases written between
/// 2098 and 2099 repeat the eight after them (the intron's last base, `CC` and
/// coding bases 1 to 5) and shift eight bases on, to stand after coding base
/// 5. The frame then reads `ATG GCG CCA TGG CTG ...`: residue 3 is a proline,
/// and the 3' UTR `CTAGGG` ends the frame 86 residues on with `AAC TAG`.
#[test]
fn a_duplication_written_in_an_intron_of_the_5_prime_utr_is_named_at_its_3_prime_place() {
    let utr5 = utr_ending(302, "CC");
    let cds = alanines_then("TAA");
    let tx = coding(make_compact_test_transcript(), 303, &utr5, &cds, "CTAGGG");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(2098, &format!("AGCC{}", &cds[..40])));
    let written = insertion(2099, b"GCCATGGC");
    let (row, alone) = rows(&tx, &fasta, &written, &insertion(2107, b"GCCATGGC"));
    assert_eq!(terms(&written, &tx), ["splice_acceptor_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.-2-1_5dup"));
    assert_eq!(row.offset, Some(8));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala3ProfsTer86")
    );
    assert_eq!(alone.offset, None);
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// On a coding sequence whose end is not found, the whole codons the
/// transcript holds give the notation: the seven bases of the first test on a
/// coding sequence of a start codon and 82 alanine codons, with no base after
/// it, change residue 3 to a histidine and meet no stop codon.
#[test]
fn the_held_codons_give_the_notation_where_the_coding_end_is_not_found() {
    let utr5 = utr_ending(50, "AC");
    let cds = alanines_then("");
    let mut tx = coding(make_compact_test_transcript(), 51, &utr5, &cds, "");
    tx.flags = vec!["cds_end_NF".to_string()].into();
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(100, &format!("{utr5}{cds}")));
    let (row, alone) = rows(
        &tx,
        &fasta,
        &insertion(148, b"ACATGGC"),
        &insertion(155, b"ACATGGC"),
    );
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.-2_5dup"));
    assert_eq!(row.offset, Some(7));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ala3HisfsTer?")
    );
    assert_eq!((&row.hgvsc, &row.hgvsp), (&alone.hgvsc, &alone.hgvsp));
}

/// Three alleles written in the 5' UTR that keep an empty HGVSp, as the same
/// change given at its 3' place does. A deletion in a run that ends in the 5'
/// UTR stays there. Two bases that repeat the 5' UTR's last two are left
/// between its last base and coding base 1. And a deletion of three bases
/// whose run ends on coding base 1 spans, shifted, the 5' UTR's last two bases
/// and that base, where no translation start is read.
#[test]
fn an_allele_the_shift_leaves_at_or_before_the_coding_start_keeps_an_empty_hgvsp() {
    let cds = alanines_then("TAA");
    let lay = |utr5: &str| {
        let tx = coding(make_compact_test_transcript(), 51, utr5, &cds, "GCTAGGGTT");
        let fixture = hgvsc_fasta_fixture(&bases_from(100, &format!("{utr5}{}", &cds[..250])));
        (tx, fixture)
    };

    let (tx, (_dir, fasta)) = lay(&utr_ending(50, "CGGTC"));
    let (row, alone) = rows(&tx, &fasta, &deletion(146, b"G"), &deletion(147, b"G"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.-3del"));
    assert_eq!(row.offset, Some(1));
    assert_eq!((row.hgvsp, alone.hgvsp), (None, None));

    let (tx, (_dir, fasta)) = lay(&utr_ending(50, "CGGTC"));
    let (row, alone) = rows(&tx, &fasta, &insertion(148, b"TC"), &insertion(150, b"TC"));
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.-2_-1dup"));
    assert_eq!(row.offset, Some(2));
    assert_eq!((row.hgvsp, alone.hgvsp), (None, None));

    let (tx, (_dir, fasta)) = lay(&utr_ending(50, "CACGACG"));
    let written = InputVariant::new("1".into(), 144, 146, b"ACG".to_vec(), b"-".to_vec());
    let there = InputVariant::new("1".into(), 148, 150, b"CGA".to_vec(), b"-".to_vec());
    let (row, alone) = rows(&tx, &fasta, &written, &there);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.-2_1del"));
    assert_eq!(row.offset, Some(4));
    assert_eq!((row.hgvsp, alone.hgvsp), (None, None));
}

/// The three-exon test transcript with intron 1 cut to `length` bases and exon
/// 2 moved up to it: the exons, the introns and the mapper agree.
fn with_first_intron_of(length: u64) -> Transcript {
    let mut tx = make_test_transcript();
    let intron_end = tx.introns[0].start + length - 1;
    let (exon_start, exon_end) = (intron_end + 1, intron_end + 300);
    tx.exons[1].start = exon_start;
    tx.exons[1].end = exon_end;
    tx.introns[0].end = intron_end;
    tx.introns[1].start = exon_end + 1;
    let vefc = tx.vefc.as_mut().expect("a coding model");
    vefc.sorted_exons[1].start = exon_start;
    vefc.sorted_exons[1].end = exon_end;
    vefc.introns[0].end = intron_end;
    vefc.introns[1].start = exon_end + 1;
    let mapper = vefc.mapper.as_mut().expect("a mapper");
    let mut pairs = mapper.exon_coord_mapper.pairs.clone();
    pairs[1].to_start = exon_start;
    pairs[1].to_end = exon_end;
    mapper.exon_coord_mapper = ExonCoordMapper::new(pairs);
    tx
}

/// An insertion the shift leaves between an intron's last base and the exon
/// after it keeps an empty HGVSp, wherever its record writes it. The coding
/// sequence starts at cDNA 299, so exon 1 ends on coding bases 1 and 2, `AT`,
/// at 398 and 399; intron 1 is 14 bases, 400..413, and repeats `AT`; exon 2
/// opens at 414 with the start codon's `G`. The 5' UTR ends `ATATAT`. `AT`
/// written between 393 and 394, in the 5' UTR, shifts 20 bases, through both
/// coding bases and the whole intron, to stand between 413 and 414: a
/// duplication of the intron's last two bases. The same two bases written in
/// the intron, between 403 and 404, shift 10 bases to the same place, and
/// written there they do not move.
#[test]
fn an_insertion_the_shift_leaves_at_the_end_of_an_intron_keeps_an_empty_hgvsp() {
    let utr5 = utr_ending(298, "ATATAT");
    let cds = alanines_then("TAA");
    let tx = coding(
        compact_on_chr1(with_first_intron_of(14)),
        299,
        &utr5,
        &cds,
        "GCTAGGGTT",
    );
    let mut reference = bases_from(100, &format!("{utr5}AT"));
    reference.extend(bases_from(400, "ATATATATATATAT"));
    reference.extend(bases_from(414, &cds[2..42]));
    let (_dir, fasta) = hgvsc_fasta_fixture(&reference);
    let in_the_utr = insertion(394, b"AT");
    let (row, at_the_junction) = rows(&tx, &fasta, &in_the_utr, &insertion(414, b"AT"));
    let in_the_intron = generate_hgvs(
        &insertion(404, b"AT"),
        &tx,
        Some(&fasta),
        HgvsOptions::default(),
    );
    assert_eq!(terms(&in_the_utr, &tx), ["5_prime_UTR_variant"]);
    for (notation, offset) in [
        (&row, Some(20)),
        (&in_the_intron, Some(10)),
        (&at_the_junction, None),
    ] {
        assert_eq!(
            notation.hgvsc.as_deref(),
            Some("ENST00000000001.1:c.3-2_3-1dup")
        );
        assert_eq!(notation.offset, offset);
        assert_eq!(notation.hgvsp, None);
    }
}

/// `USED_REF` stays the bases the record wrote. Ensembl VEP rewrites the
/// reference allele with the shifted bases only past the test the allele as
/// written fails (TranscriptVariationAllele.pm:1730-1737 before 1771-1786).
/// The 5' UTR ends `ATGGCT` and the coding sequence opens `ATG GCT ATG GAA`:
/// deleting the 5' UTR's last six bases is deleting coding bases 5 to 10, ten
/// bases on, which read `CTATGG`. The row keeps `ATGGCT` beside the protein
/// notation.
#[test]
fn used_ref_of_an_allele_written_in_the_5_prime_utr_is_not_rotated_by_the_notation() {
    let utr5 = utr_ending(50, "CCATGGCT");
    let cds = format!("ATGGCTATGGAA{}TAA", "GCT".repeat(80));
    let tx = coding(make_compact_test_transcript(), 51, &utr5, &cds, "GCTAGGGTT");
    let (_dir, fasta) = hgvsc_fasta_fixture(&bases_from(100, &format!("{utr5}{}", &cds[..250])));
    let written = InputVariant::new("1".into(), 144, 149, b"ATGGCT".to_vec(), b"-".to_vec());
    let config = EffectsConfig {
        reference_fasta: Some(Arc::new(fasta)),
        compute_hgvs: true,
        use_transcript_ref: true,
        ..Default::default()
    };
    let row =
        calculate_consequences(&written, &Arc::new(tx), &config).expect("a row on the transcript");
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.5_10del"));
    assert_eq!(row.hgvs_offset, Some(10));
    assert!(row.hgvsp.is_some());
    assert_eq!(row.used_ref.as_deref(), Some("ATGGCT"));
}

/// The bounds of `annotated_before_the_coding_start`. Ensembl VEP reads `utr`
/// alone where the ordered span does not overlap the coding region
/// (BaseVariationFeatureOverlapAllele.pm:463, 523-525) and reads the
/// transcript's own span against the ends as given (:480). The forward
/// transcript runs 100..6100 with coding bases 150..2101; the reverse one
/// starts at 6100 with coding base 1 at 6050 and its last coding base at 4099.
/// A span that holds coding base 1, as its last base or as the base after an
/// insertion, overlaps the coding region.
#[test]
fn a_span_is_before_the_coding_start_only_where_it_ends_before_coding_base_1() {
    let cds = alanines_then("TAA");
    let utr5 = utr_ending(50, "AC");

    let forward = coding(make_compact_test_transcript(), 51, &utr5, &cds, "GCTAGGGTT");
    assert_eq!(
        (forward.coding_region_start, forward.coding_region_end),
        (Some(150), Some(2101))
    );
    for (span, start, end, before) in [
        ("an insertion between 148 and 149", 149, 148, true),
        ("a deletion that ends on 149", 140, 149, true),
        (
            "an insertion between 149 and coding base 1",
            150,
            149,
            false,
        ),
        ("a deletion that ends on coding base 1", 149, 150, false),
        ("coding base 1", 150, 150, false),
        (
            "an insertion before the transcript's first base",
            100,
            99,
            false,
        ),
        ("two bases before the transcript", 98, 99, false),
        ("a base of the 3' UTR", 2110, 2110, false),
    ] {
        assert_eq!(
            crate::coding::annotated_before_the_coding_start(&forward, start, end),
            before,
            "forward strand: {span}"
        );
    }

    let reverse = coding(
        compact_on_chr1(make_descending_reverse_strand_transcript()),
        51,
        &utr5,
        &cds,
        "GCTAGGGTT",
    );
    assert_eq!(
        (reverse.coding_region_start, reverse.coding_region_end),
        (Some(4099), Some(6050))
    );
    for (span, start, end, before) in [
        ("an insertion between 6051 and 6052", 6052, 6051, true),
        ("a deletion that ends on 6051", 6051, 6060, true),
        (
            "an insertion between coding base 1 and 6051",
            6051,
            6050,
            false,
        ),
        ("a deletion that ends on coding base 1", 6050, 6051, false),
        ("coding base 1", 6050, 6050, false),
        (
            "an insertion before the transcript's first base",
            6101,
            6100,
            false,
        ),
        ("two bases before the transcript", 6101, 6102, false),
        ("a base of the 3' UTR", 4090, 4090, false),
    ] {
        assert_eq!(
            crate::coding::annotated_before_the_coding_start(&reverse, start, end),
            before,
            "reverse strand: {span}"
        );
    }
}
