// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! An insertion the 3' shift leaves between a coding exon's last base and the
//! first base of the intron after it, whose bases are no copy of the bases
//! before that gap: an insertion flanked by the exon's last position and the
//! intron's first, which carries no protein notation.

use super::{terms, JunctionModel, ENST00000227471_EXONS_7_8_9, ENST00000312957_EXONS_4_5};
use vep_core::coordinate::Strand;

/// GRCh37 6, ENST00000348925, exons 9 and 10 of 11: exon 9 ends at 24528394
/// on coding base 1382, the second base of codon 461, before an intron of
/// 3,952 bases. The exon's last two bases differ from each other and from
/// the intron's first base.
const ENST00000348925_EXONS_9_10: JunctionModel = JunctionModel {
    chr: "6",
    strand: Strand::Forward,
    exons: &[(24528225, 24528394), (24532347, 24532405)],
    first_cdna: 1241,
    cdna_coding_start: 29,
    bases: &[(24528393, "ACGT")],
};

/// GRCh37 1, ENST00000368147, a reverse-strand transcript, exons 31 and 32
/// of 52: exon 31 ends at 158613112 on coding base 4442, the second base of
/// codon 1481, before an intron of 345 bases whose first base is 158613111.
const ENST00000368147_EXONS_31_32: JunctionModel = JunctionModel {
    chr: "1",
    strand: Strand::Reverse,
    exons: &[(158613112, 158613215), (158612604, 158612766)],
    first_cdna: 4520,
    cdna_coding_start: 182,
    bases: &[(158613110, "AC")],
};

/// GRCh37 3, ENST00000231790, exons 4 and 5 of 19: exon 4 ends at 37045965
/// on coding base 380, the second base of codon 127, before an intron of
/// 2,516 bases.
const ENST00000231790_EXONS_4_5: JunctionModel = JunctionModel {
    chr: "3",
    strand: Strand::Forward,
    exons: &[(37045892, 37045965), (37048482, 37048554)],
    first_cdna: 523,
    cdna_coding_start: 217,
    bases: &[(37045966, "GT")],
};

/// GRCh37 22, ENST00000417588, a reverse-strand transcript, exons 2 to 4 of
/// 13: exon 3, of 148 bases, ends at 29120965 on coding base 592, the first
/// base of codon 198, before an intron of 12,959 bases whose first base is
/// 29120964; exon 2 holds 125 bases and the intron between them 118. The
/// coding sequence opens in exon 1 and ends at cDNA position 612, in exon 4,
/// on a stop codon. Exon 4 opens on bases that differ from the first bases
/// of exon 2, and bases 25 to 27 of the intron between exons 2 and 3 are a
/// stop codon.
const ENST00000417588_EXONS_2_3_4: JunctionModel = JunctionModel {
    chr: "22",
    strand: Strand::Reverse,
    exons: &[
        (29121231, 29121355),
        (29120965, 29121112),
        (29107897, 29108005),
    ],
    first_cdna: 320,
    cdna_coding_start: 1,
    bases: &[
        (29107986, "TTA"),
        (29107997, "GCGATGCTGCTGGAC"),
        (29120963, "AC"),
        (29121104, "GCGATGCTGCTGGAC"),
        (29121204, "TTA"),
    ],
};

/// The shape of the row at GRCh37 3:37045965 on ENST00000231790: three
/// bases that are no copy of the bases before them, written between the
/// last base of exon 4 and the intron's first base, where the 3' rule
/// leaves them. The coding notation is an insertion between coding base 380
/// and the intron's first position; the codon window reads the inserted
/// bases as the exon's last, an insertion of one residue between residues 127 and 128. One base written after the
/// last base of exon 4 of ENST00000312957 is the same form.
#[test]
fn insertion_left_after_an_exons_last_base_has_no_protein_notation() {
    let row = ENST00000231790_EXONS_4_5.row(37045966, 37045965, "-", "CAC");
    assert_eq!(terms(&row), ["inframe_insertion", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.380_380+1insCAC")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("380-381"));
    assert_eq!(row.protein_position.as_deref(), Some("127"));
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000312957_EXONS_4_5.row(27066221, 27066220, "-", "C");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.394_394+1insC")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("394-395"));
    assert_eq!(row.protein_position.as_deref(), Some("132"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// The shape of the row at GRCh37 1:158613111 on ENST00000368147, a
/// reverse-strand transcript: three bases that are no copy of the bases
/// before them, written between the last base of exon 31, 158613112, and
/// the intron's first base. The codon window gives a replacement at residue 1481 that brings a stop.
#[test]
fn insertion_left_after_an_exons_last_base_on_the_reverse_strand_has_no_protein_notation() {
    let row = ENST00000368147_EXONS_31_32.row(158613112, 158613111, "-", "TAG");
    assert_eq!(
        terms(&row),
        [
            "protein_altering_variant",
            "splice_region_variant",
            "stop_gained"
        ]
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.4442_4442+1insCTA")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("4442-4443"));
    assert_eq!(row.protein_position.as_deref(), Some("1481"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// The shape of the row at GRCh37 6:24528392 on ENST00000348925: three
/// bases written two bases before the last base of exon 9, which the 3'
/// rule carries 2 bases, to between that base and the intron's first. The
/// last inserted base is then the exon's own last base, and the three are
/// no copy of the three before them. The codon window gives
/// a duplication of residue 461.
#[test]
fn insertion_the_shift_carries_to_an_exons_last_base_has_no_protein_notation() {
    let row = ENST00000348925_EXONS_9_10.row(24528393, 24528392, "-", "ACA");
    assert_eq!(terms(&row), ["inframe_insertion", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.1382_1382+1insAAC")
    );
    assert_eq!(row.hgvs_offset, Some(2));
    assert_eq!(row.cds_position.as_deref(), Some("1380-1381"));
    assert_eq!(row.protein_position.as_deref(), Some("460-461"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// One base that is no copy, written between the last base of exon 7 of
/// ENST00000227471 and the first base of an intron of 2 bases: the form
/// holds before an intron of any length.
#[test]
fn insertion_left_before_an_intron_of_two_bases_has_no_protein_notation() {
    let row = ENST00000227471_EXONS_7_8_9.row(67765166, 67765165, "-", "G");
    assert_eq!(
        terms(&row),
        [
            "frameshift_variant",
            "intron_variant",
            "splice_region_variant"
        ]
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.885_885+1insC")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("885-886"));
    assert_eq!(row.protein_position.as_deref(), Some("295-296"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// A copy of the bases before the gap is a duplication of coding positions
/// and keeps the notation its bases give: one base equal to the last base of
/// exon 4 of ENST00000231790, three equal to its last three, and the shape of
/// the row at GRCh37 22:29120964 on ENST00000417588, 391 bases equal to the
/// 391 before them, written after the last base of an exon of 148 bases.
/// That copy holds the whole exon, the intron before it and the whole exon
/// before that; the reference shows it a copy all the same, and the frame
/// it opens ends on the stop codon the copied intron holds.
///
/// Exons 4 and 5 of the first model hold the filler base alone, so codons
/// 103 to 151, the model's last, read one residue: the frame one copied
/// base opens reads the reference's residues to the last of them and has
/// no protein notation, where Ensembl VEP's walk of the two translations
/// leaves by its bound one position on, at 152
/// (`TranscriptVariationAllele.pm:2348-2364`), and the residue three copied
/// bases add is named at the last of that run, 151, both past the row's
/// protein position, 127. On the third row the copy's first two bases close
/// codon 198 with the exon's last base, so the first changed residue is the
/// row's own, and the copied intron's stop codon is the 51st codon from it.
#[test]
fn copy_of_the_bases_before_an_introns_first_base_keeps_its_notation() {
    let row = ENST00000231790_EXONS_4_5.row(37045966, 37045965, "-", "T");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.380dup"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("380-381"));
    assert_eq!(row.protein_position.as_deref(), Some("127"));
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000231790_EXONS_4_5.row(37045966, 37045965, "-", "TTT");
    assert_eq!(terms(&row), ["inframe_insertion", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.378_380dup"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("380-381"));
    assert_eq!(row.protein_position.as_deref(), Some("127"));
    assert_eq!(row.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Phe151dup"));

    let row = ENST00000417588_EXONS_2_3_4.row_ending(Some(612), 29120965, 29120964, "-", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGCGATGCTGCTGGACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    assert_eq!(
        terms(&row),
        ["frameshift_variant", "splice_region_variant", "stop_gained"]
    );
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.320_592dup"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("592-593"));
    assert_eq!(row.protein_position.as_deref(), Some("198"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Thr198LysfsTer51")
    );
}

/// One base equal to the intron's first base, written between the last
/// base of exon 4 of ENST00000231790 and that base: the 3' rule carries it
/// 1 base into the intron, where it is a duplication of the intron's first
/// position.
#[test]
fn insertion_the_three_prime_rule_carries_into_the_intron_keeps_its_value() {
    let row = ENST00000231790_EXONS_4_5.row(37045966, 37045965, "-", "G");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.380+1dup"));
    assert_eq!(row.hgvs_offset, Some(1));
    assert_eq!(row.cds_position.as_deref(), Some("380-381"));
    assert_eq!(row.protein_position.as_deref(), Some("127"));
    assert_eq!(row.hgvsp.as_deref(), None);
}
