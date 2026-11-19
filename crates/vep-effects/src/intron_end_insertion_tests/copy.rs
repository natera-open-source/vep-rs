// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! One inserted base equal to the intron's last base, written between that
//! base and the exon's first base: a duplication of the intron's last position,
//! which carries no protein notation.

use super::{terms, JunctionModel, ENST00000312957_EXONS_4_5};
use vep_core::coordinate::Strand;

/// GRCh37 18, ENST00000389658, a reverse-strand transcript, exons 45 and 46
/// of 63: exon 46 opens at 6975035 on coding base 6490, the first base of
/// codon 2164, behind an intron of 900 bases whose last base is 6975036.
const ENST00000389658_EXONS_45_46: JunctionModel = JunctionModel {
    chr: "18",
    strand: Strand::Reverse,
    exons: &[(6975936, 6976079), (6974902, 6975035)],
    first_cdna: 6440,
    cdna_coding_start: 95,
    bases: &[(6975027, "GCGATGCTGCTGGAC")],
};

/// GRCh37 7, exons 30 to 32 of ENST00000538115, which ENST00000542585 and
/// ENST00000546292 hold at the same cDNA positions behind the same coding
/// start: exon 30 ends at 100371475, an intron of 2 bases follows, and exon
/// 31 opens at 100371478 on coding base 5767, the first base of codon 1923.
/// Exon 30's last base and the intron's two bases are one base three times.
const ENST00000538115_EXONS_30_31_32: JunctionModel = JunctionModel {
    chr: "7",
    strand: Strand::Forward,
    exons: &[
        (100371344, 100371475),
        (100371478, 100371487),
        (100372950, 100373197),
    ],
    first_cdna: 5783,
    cdna_coding_start: 149,
    bases: &[(100371475, "GGGCAGCATCGC")],
};

/// GRCh37 7, ENST00000546213, exons 30 to 32 of 46, the same junction as
/// [`ENST00000538115_EXONS_30_31_32`] behind a later coding start: exon 31
/// opens on coding base 1300, the first base of codon 434.
const ENST00000546213_EXONS_30_31_32: JunctionModel = JunctionModel {
    chr: "7",
    strand: Strand::Forward,
    exons: &[
        (100371404, 100371475),
        (100371478, 100371487),
        (100372950, 100373197),
    ],
    first_cdna: 5645,
    cdna_coding_start: 4418,
    bases: &[(100371475, "GGGCAGCATCGC")],
};

/// The shape of the row at GRCh37 21:27070987 on ENST00000312957, written
/// at the junction: one base equal to the intron's last base, between that
/// base and the first base of exon 5. The coding notation is a duplication
/// of the intron's last position; the codon window reads the inserted base
/// as the exon's first, a frameshift at residue 132.
#[test]
fn copy_of_an_introns_last_base_written_before_the_exon_has_no_protein_notation() {
    let row = ENST00000312957_EXONS_4_5.row(27070989, 27070988, "-", "G");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.395-1dup"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("394-395"));
    assert_eq!(row.protein_position.as_deref(), Some("132"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// The shape of the row at GRCh37 18:6975035 on ENST00000389658, a
/// reverse-strand transcript: one base equal to the intron's last base,
/// 6975036, written between it and the first base of exon 46. The codon
/// window gives a frameshift at residue 2164.
#[test]
fn copy_of_an_introns_last_base_on_the_reverse_strand_has_no_protein_notation() {
    let row = ENST00000389658_EXONS_45_46.row(6975036, 6975035, "-", "C");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.6490-1dup"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("6489-6490"));
    assert_eq!(row.protein_position.as_deref(), Some("2163-2164"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// The shape of the four rows of the record at GRCh37 7:100371474, on
/// ENST00000538115, ENST00000542585, ENST00000546292 and ENST00000546213:
/// one base equal to the three after it, written before the last base of
/// exon 30, which the 3' rule carries 3 bases, over that base and an
/// intron of 2 bases, to before the first base of exon 31. The input
/// writes the allele inside exon 30 and not at the end of an intron, so
/// the rows keep the protein notation Ensembl VEP writes there too,
/// a frameshift at residue 1923 and at residue 434.
#[test]
fn copy_written_in_the_exon_before_an_intron_of_two_bases_keeps_its_notation() {
    let row = ENST00000538115_EXONS_30_31_32.row(100371475, 100371474, "-", "G");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.5767-1dup"));
    assert_eq!(row.hgvs_offset, Some(3));
    assert_eq!(row.cds_position.as_deref(), Some("5765-5766"));
    assert_eq!(row.protein_position.as_deref(), Some("1922"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Gln1923AlafsTer?")
    );

    let row = ENST00000546213_EXONS_30_31_32.row(100371475, 100371474, "-", "G");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.1300-1dup"));
    assert_eq!(row.hgvs_offset, Some(3));
    assert_eq!(row.cds_position.as_deref(), Some("1298-1299"));
    assert_eq!(row.protein_position.as_deref(), Some("433"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Gln434AlafsTer?")
    );
}
