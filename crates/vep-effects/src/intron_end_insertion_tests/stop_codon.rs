// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! An insertion written between an intron's last base and the first base of an
//! exon that opens inside the stop codon: it carries no protein notation, as
//! before any other coding exon.

use super::{terms, JunctionModel};
use vep_core::coordinate::Strand;

/// GRCh37 21, ENST00000399326, exons 23 and 24 of 37: exon 24 opens on coding
/// base 3048, the last base of the stop codon, behind an intron of 2,760
/// bases; exon 23 ends on the stop codon's first two bases. The three bases
/// are a stop codon, and the intron's last base differs from the exon's
/// first.
const ENST00000399326_EXONS_23_24: JunctionModel = JunctionModel {
    chr: "21",
    strand: Strand::Forward,
    exons: &[(35199122, 35199167), (35201928, 35202049)],
    first_cdna: 3003,
    cdna_coding_start: 2,
    bases: &[(35199167, "AGT"), (35201922, "GTCCAGA")],
};

/// Exon 24 of ENST00000399326 opens on the stop codon's last base. One base
/// equal to the intron's last base, six bases equal to the intron's last
/// six, one base that is no copy, and three bases ending in the two bases
/// that end an intron, each written between the intron's last base and
/// that exon's first base: the change is described at the end of the
/// intron as it is before any other exon, and the row carries no protein
/// notation, where the codon window gives the stop written as unchanged or a replacement of the stop, at residue 1016.
#[test]
fn insertion_before_an_exon_that_opens_inside_the_stop_codon_has_no_protein_notation() {
    let row = ENST00000399326_EXONS_23_24.row_ending(Some(3049), 35201928, 35201927, "-", "G");
    assert_eq!(
        terms(&row),
        ["splice_region_variant", "stop_retained_variant"]
    );
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.*-1dup"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("3047-3048"));
    assert_eq!(row.protein_position.as_deref(), Some("1016"));
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000399326_EXONS_23_24.row_ending(Some(3049), 35201928, 35201927, "-", "GTCCAG");
    assert_eq!(
        terms(&row),
        [
            "inframe_insertion",
            "splice_region_variant",
            "stop_retained_variant"
        ]
    );
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.*-6_*-1dup"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("3047-3048"));
    assert_eq!(row.protein_position.as_deref(), Some("1016"));
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000399326_EXONS_23_24.row_ending(Some(3049), 35201928, 35201927, "-", "C");
    assert_eq!(terms(&row), ["splice_region_variant", "stop_lost"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.*-1_3048insC")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("3047-3048"));
    assert_eq!(row.protein_position.as_deref(), Some("1016"));
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000399326_EXONS_23_24.row_ending(Some(3049), 35201928, 35201927, "-", "TAG");
    assert_eq!(terms(&row), ["splice_region_variant", "stop_lost"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.*-1_3048insTAG")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("3047-3048"));
    assert_eq!(row.protein_position.as_deref(), Some("1016"));
    assert_eq!(row.hgvsp.as_deref(), None);
}
