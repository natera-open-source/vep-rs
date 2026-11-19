// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Inserted bases that are no copy of the bases before them and cannot be the
//! intron's new end, written between the intron's last base and the exon's
//! first base: an insertion flanked by the intron's last position and the
//! exon's first, which carries no protein notation.

use super::{terms, JunctionModel, ENST00000321219_EXONS_8_9};
use vep_core::coordinate::Strand;

/// GRCh37 15, ENST00000396659, a reverse-strand transcript, exons 6 and 7 of
/// 9: exon 7 opens at 45657058 on coding base 979, the first base of codon
/// 327, behind an intron of 1,185 bases whose last base is 45657059.
const ENST00000396659_EXONS_6_7: JunctionModel = JunctionModel {
    chr: "15",
    strand: Strand::Reverse,
    exons: &[(45658244, 45658408), (45656995, 45657058)],
    first_cdna: 1154,
    cdna_coding_start: 341,
    bases: &[(45657050, "GCGATGCTGCTGGAC")],
};

/// The shape of the row at GRCh37 21:38877585 on ENST00000321219: three
/// bases that are no copy of the bases before them and do not end in the
/// two bases that end an intron, written between the intron's last base and
/// the first base of exon 9. The coding notation is an insertion between
/// the intron's last position and coding base 1240; the codon window reads
/// the inserted bases as the exon's first, an insertion of a stop between residues 413 and 414.
#[test]
fn insertion_that_cannot_end_the_intron_written_before_the_exon_has_no_protein_notation() {
    let row = ENST00000321219_EXONS_8_9.row(38877586, 38877585, "-", "TAA");
    assert_eq!(
        terms(&row),
        ["inframe_insertion", "splice_region_variant", "stop_gained"]
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.1240-1_1240insTAA")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("1239-1240"));
    assert_eq!(row.protein_position.as_deref(), Some("413-414"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// The shape of the row at GRCh37 15:45657058 on ENST00000396659, a
/// reverse-strand transcript: one base that is no copy of the intron's last
/// base, 45657059, written between it and the first base of exon 7. The
/// codon window gives a frameshift at residue 327.
#[test]
fn insertion_that_cannot_end_the_intron_on_the_reverse_strand_has_no_protein_notation() {
    let row = ENST00000396659_EXONS_6_7.row(45657059, 45657058, "-", "A");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.979-1_979insT")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("978-979"));
    assert_eq!(row.protein_position.as_deref(), Some("326-327"));
    assert_eq!(row.hgvsp.as_deref(), None);
}
