// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Inserted bases that are no copy of the bases before them and end in the two
//! bases that end an intron, written between the intron's last base and the
//! exon's first base: an insertion flanked by the intron's last position and
//! the exon's first, which carries no protein notation.

use super::{terms, JunctionModel, ENST00000227471_EXONS_7_8_9, ENST00000321219_EXONS_8_9};
use vep_core::coordinate::Strand;

/// GRCh37 10, ENST00000353479, a reverse-strand transcript, exons 38 and 39
/// of 56: exon 39 opens at 105800876 on coding base 2648, the second base of
/// codon 883, behind an intron of 184 bases whose last base is 105800877.
pub(super) const ENST00000353479_EXONS_38_39: JunctionModel = JunctionModel {
    chr: "10",
    strand: Strand::Reverse,
    exons: &[(105801061, 105801102), (105800823, 105800876)],
    first_cdna: 2897,
    cdna_coding_start: 292,
    bases: &[(105800868, "GCGATGCTGCTGGAC")],
};

/// Three bases that are no copy of the bases before them and end in the two
/// bases that end an intron, of the major kind and then of the minor kind,
/// written between the last base of intron 8 of ENST00000321219 and the
/// first base of exon 9. The coding notation is an insertion between the
/// intron's last position and coding base 1240; the codon window reads the
/// inserted bases as the exon's first, an insertion of one residue between residues 413 and 414.
#[test]
fn insertion_ending_in_the_bases_that_end_an_intron_written_before_the_exon_has_no_protein_notation(
) {
    let row = ENST00000321219_EXONS_8_9.row(38877586, 38877585, "-", "GAG");
    assert_eq!(terms(&row), ["inframe_insertion", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.1240-1_1240insGAG")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("1239-1240"));
    assert_eq!(row.protein_position.as_deref(), Some("413-414"));
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000321219_EXONS_8_9.row(38877586, 38877585, "-", "GAC");
    assert_eq!(terms(&row), ["inframe_insertion", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.1240-1_1240insGAC")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("1239-1240"));
    assert_eq!(row.protein_position.as_deref(), Some("413-414"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// The shape of the row at GRCh37 10:105800876 on ENST00000353479, a
/// reverse-strand transcript: 16 bases that are no copy of the bases before
/// them and end in the two bases that end an intron, written between the
/// intron's last base, 105800877, and the first base of exon 39. The codon
/// window gives a frameshift at residue 883.
#[test]
fn insertion_ending_in_the_bases_that_end_an_intron_on_the_reverse_strand_has_no_protein_notation()
{
    let row = ENST00000353479_EXONS_38_39.row(105800877, 105800876, "-", "CTAAGCTAGCTAGCTA");
    assert_eq!(
        terms(&row),
        ["frameshift_variant", "splice_region_variant", "stop_gained"]
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.2648-1_2648insTAGCTAGCTAGCTTAG")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("2647-2648"));
    assert_eq!(row.protein_position.as_deref(), Some("883"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// Three bases ending in the two bases that end an intron, written between
/// the last base of an intron of 2 bases of ENST00000227471 and the exon's
/// first base: Ensembl VEP reads an intron of at most 13 bases with the
/// exons beside it, and the row keeps its protein notation.
#[test]
fn insertion_ending_in_the_bases_that_end_an_intron_after_an_intron_of_two_bases_keeps_its_notation(
) {
    let row = ENST00000227471_EXONS_7_8_9.row(67765164, 67765163, "-", "CTA");
    assert_eq!(
        terms(&row),
        [
            "inframe_insertion",
            "intron_variant",
            "splice_region_variant",
            "stop_gained"
        ]
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.886-1_886insTAG")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("885-886"));
    assert_eq!(row.protein_position.as_deref(), Some("295-296"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Lys295_Gln296insTer")
    );
}
