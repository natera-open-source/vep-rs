// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Two or more inserted bases equal to the bases of the intron that end on its
//! last base, written between that base and the exon's first base: a
//! duplication of intron positions alone, which carries no protein notation.

use super::{terms, JunctionModel, ENST00000297268_EXONS_36_37};
use vep_core::coordinate::Strand;

/// GRCh37 14, ENST00000206765, a reverse-strand transcript, exons 14 and 15
/// of 15: exon 15 opens at 24718747 on coding base 2226, the third base of
/// codon 742, behind an intron of 4,610 bases whose last base is 24718748.
const ENST00000206765_EXONS_14_15: JunctionModel = JunctionModel {
    chr: "14",
    strand: Strand::Reverse,
    exons: &[(24723358, 24723494), (24718320, 24718747)],
    first_cdna: 2213,
    cdna_coding_start: 125,
    bases: &[(24718739, "GCGATGCTGCTGGAC")],
};

/// The shape of the row at GRCh37 7:94049846 on ENST00000297268, written at
/// the junction: six bases equal to the intron's last six, between the
/// intron's last base and the first base of exon 37. The coding notation is
/// a duplication of six intron positions; the codon window reads the
/// inserted bases as the exon's first, an insertion of two residues between residues 729 and 730.
#[test]
fn duplication_ending_on_an_introns_last_base_written_before_the_exon_has_no_protein_notation() {
    let row = ENST00000297268_EXONS_36_37.row(94049853, 94049852, "-", "GTCCAG");
    assert_eq!(terms(&row), ["inframe_insertion", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.2188-6_2188-1dup")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("2187-2188"));
    assert_eq!(row.protein_position.as_deref(), Some("729-730"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// The shape of the row at GRCh37 14:24718747 on ENST00000206765, a
/// reverse-strand transcript: two bases equal to the intron's last two,
/// 24718749 and 24718748, written between the intron's last base and the
/// first base of exon 15. The codon window gives a frameshift at residue 742.
#[test]
fn duplication_ending_on_an_introns_last_base_on_the_reverse_strand_has_no_protein_notation() {
    let row = ENST00000206765_EXONS_14_15.row(24718748, 24718747, "-", "CT");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.2226-2_2226-1dup")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("2225-2226"));
    assert_eq!(row.protein_position.as_deref(), Some("742"));
    assert_eq!(row.hgvsp.as_deref(), None);
}
