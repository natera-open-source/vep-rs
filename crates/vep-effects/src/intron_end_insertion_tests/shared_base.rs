// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! An insertion at either end of an intron that the record holds with a
//! reference base it shares, as a record that carries another allele beside the
//! insertion holds it: the change is the insertion the coding notation
//! describes, and it carries no protein notation where the same insertion
//! written alone carries none.

use super::acceptor_pair::ENST00000353479_EXONS_38_39;
use super::{
    terms, ENST00000297268_EXONS_36_37, ENST00000312957_EXONS_4_5, ENST00000321219_EXONS_8_9,
};

/// The insertion of the row at GRCh37 10:105800876 on ENST00000353479 as a
/// record that carries another allele writes it: the reference allele is
/// the base at 105800876, the first base of exon 39 on this reverse-strand
/// transcript, and the alternate allele that base and the 16 inserted
/// bases after it. The coding notation is the insertion between the
/// intron's last position and the exon's first coding base that the same
/// insertion written alone has, and the row carries no protein notation as
/// that one does.
#[test]
fn insertion_written_with_the_base_before_it_before_an_exon_has_no_protein_notation() {
    let row = ENST00000353479_EXONS_38_39.row(105800876, 105800876, "G", "GCTAAGCTAGCTAGCTA");
    assert_eq!(
        terms(&row),
        ["frameshift_variant", "splice_region_variant", "stop_gained"]
    );
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.2648-1_2648insTAGCTAGCTAGCTTAG")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("2648"));
    assert_eq!(row.protein_position.as_deref(), Some("883"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// One base that is no copy, written after the last base of exon 4 of
/// ENST00000312957 with that base before it as the reference allele: the
/// insertion between the exon's last coding base and the intron's first
/// position, which carries no protein notation.
#[test]
fn insertion_written_with_the_base_before_it_after_an_exons_last_base_has_no_protein_notation() {
    let row = ENST00000312957_EXONS_4_5.row(27066220, 27066220, "T", "TC");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.394_394+1insC")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("394"));
    assert_eq!(row.protein_position.as_deref(), Some("132"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// The four forms of an insertion at the end of an intron, each written
/// with the exon's first base after it as the reference allele, on
/// forward-strand transcripts: one base equal to the intron's last base,
/// six equal to its last six, three that cannot end the intron, and three
/// ending in the two bases that end an intron. The shared base is clipped
/// from the back as the coding notation clips it, and the rows carry no
/// protein notation.
#[test]
fn insertion_written_with_the_base_after_it_before_an_exon_has_no_protein_notation() {
    let row = ENST00000312957_EXONS_4_5.row(27070989, 27070989, "C", "GC");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.395-1_395insG")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("395"));
    assert_eq!(row.protein_position.as_deref(), Some("132"));
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000297268_EXONS_36_37.row(94049853, 94049853, "C", "GTCCAGC");
    assert_eq!(terms(&row), ["inframe_insertion", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.2188-1_2188insGTCCAG")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("2188"));
    assert_eq!(row.protein_position.as_deref(), Some("730"));
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000321219_EXONS_8_9.row(38877586, 38877586, "C", "TAAC");
    assert_eq!(terms(&row), ["splice_region_variant", "stop_gained"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.1240-1_1240insTAA")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("1240"));
    assert_eq!(row.protein_position.as_deref(), Some("414"));
    assert_eq!(row.hgvsp.as_deref(), None);

    let row = ENST00000321219_EXONS_8_9.row(38877586, 38877586, "C", "GAGC");
    assert_eq!(terms(&row), ["inframe_insertion", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.1240-1_1240insGAG")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("1240"));
    assert_eq!(row.protein_position.as_deref(), Some("414"));
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// Three bases written before the first base of exon 9 of ENST00000321219
/// with the intron's last base before them as the reference allele: the
/// allele as the record holds it lies in the intron, and the row carries no
/// protein notation whatever its form. Its splice term is the one of the
/// insertion between the intron's last base and the exon, where Ensembl VEP
/// reads the intron's last base as changed
/// (VariationFeatureOverlapAllele.pm:452-454) and writes the acceptor term.
#[test]
fn insertion_written_with_an_intron_base_before_it_keeps_its_empty_notation() {
    let row = ENST00000321219_EXONS_8_9.row(38877585, 38877585, "G", "GTAA");
    assert_eq!(terms(&row), ["splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.1240-1_1240insTAA")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), None);
    assert_eq!(row.protein_position.as_deref(), None);
    assert_eq!(row.hgvsp.as_deref(), None);
}

/// Beside the last base of exon 4 of ENST00000312957: a copy of
/// that base written with it as the reference allele, a duplication of a
/// coding position; two bases replaced by three that share the first, of
/// which a base is replaced; and a substitution of the exon's last base. The
/// third base of the three stands on the intron's first base, past the span,
/// which gives no donor term: the replaced exonic base gives
/// `splice_region_variant`.
#[test]
fn allele_sharing_a_base_that_is_no_insertion_at_a_junction_keeps_its_notation() {
    let row = ENST00000312957_EXONS_4_5.row(27066220, 27066220, "T", "TT");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.394dup"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("394"));
    assert_eq!(row.protein_position.as_deref(), Some("132"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ser132PhefsTer?")
    );

    let row = ENST00000312957_EXONS_4_5.row(27066219, 27066220, "TT", "TAA");
    assert_eq!(terms(&row), ["frameshift_variant", "splice_region_variant"]);
    assert_eq!(
        row.hgvsc.as_deref(),
        Some("ENST00000000001.1:c.394delinsAA")
    );
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("393-394"));
    assert_eq!(row.protein_position.as_deref(), Some("131-132"));
    assert_eq!(
        row.hgvsp.as_deref(),
        Some("ENSP00000000001.1:p.Ser132AsnfsTer?")
    );

    let row = ENST00000312957_EXONS_4_5.row(27066220, 27066220, "T", "A");
    assert_eq!(terms(&row), ["missense_variant", "splice_region_variant"]);
    assert_eq!(row.hgvsc.as_deref(), Some("ENST00000000001.1:c.394T>A"));
    assert_eq!(row.hgvs_offset, None);
    assert_eq!(row.cds_position.as_deref(), Some("394"));
    assert_eq!(row.protein_position.as_deref(), Some("132"));
    assert_eq!(row.hgvsp.as_deref(), Some("ENSP00000000001.1:p.Ser132Thr"));
}
