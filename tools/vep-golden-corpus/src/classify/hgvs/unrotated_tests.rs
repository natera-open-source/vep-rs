// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the class of an insertion left unrotated: the two readings of a moved
//! allele, and the rows the predicate takes and leaves.

use super::*;

const REFERENCE: &str = "ENSP1.1:p.Trp365PhefsTer42";
const CANDIDATE: &str = "ENSP1.1:p.Trp365SerfsTer42";

fn key_of(allele: &str, feature_type: &str) -> FieldKey {
    FieldKey {
        location: "21:100-101".to_string(),
        allele: allele.to_string(),
        feature: "ENST1".to_string(),
        feature_type: feature_type.to_string(),
        consequences: "frameshift_variant".to_string(),
    }
}

/// A row of the coding change with the offset given, and `changes` over its fields.
fn row(offset: &str, changes: &[(&str, &str)]) -> FieldRow {
    let mut row: FieldRow = [
        ("cDNA_position", "1200-1201"),
        ("CDS_position", "1090-1091"),
        ("HGVSc", "ENST1.1:c.1092_1093dup"),
        ("HGVS_OFFSET", offset),
    ]
    .iter()
    .map(|(name, value)| (name.to_string(), value.to_string()))
    .collect();
    for (name, value) in changes {
        row.insert(name.to_string(), value.to_string());
    }
    row
}

/// The predicate on an `HGVSp` of a reverse-strand transcript's row of `allele`.
fn taken(allele: &str, perl_row: &FieldRow, rust_row: &FieldRow) -> bool {
    let key = key_of(allele, "Transcript");
    unrotated_shifted_insertion(
        "HGVSp",
        &key,
        REFERENCE,
        CANDIDATE,
        "-1",
        Some(perl_row),
        Some(rust_row),
    )
}

#[test]
fn the_two_readings_of_a_moved_allele_part_past_its_length_on_a_reverse_strand() {
    // Ensembl VEP's count of turns on a reverse-strand transcript is the allele's length
    // less the shift, none once the shift is the longer
    // (TranscriptVariationAllele.pm:1337-1339); the shift's own count is taken modulo the
    // allele's length on either strand.
    let readings = |allele: &str, offset: i128, reverse: bool| {
        (
            hgvs_shifted_allele(allele, offset, reverse),
            shift_rotated_allele(allele, offset, reverse),
        )
    };
    let both = |text: &str| (text.to_string(), text.to_string());
    assert_eq!(
        readings("GA", 3, true),
        ("GA".to_string(), "AG".to_string())
    );
    assert_eq!(
        readings("GAT", 4, true),
        ("GAT".to_string(), "TGA".to_string())
    );
    assert_eq!(readings("GAT", 2, true), both("ATG"));
    assert_eq!(readings("GAT", 3, true), both("GAT"));
    assert_eq!(readings("GAT", 6, true), both("GAT"));
    assert_eq!(readings("GAT", 4, false), both("ATG"));
    // The sign of the offset changes neither reading, and no bases turn to none.
    assert_eq!(readings("GAT", -4, true), readings("GAT", 4, true));
    assert_eq!(shift_rotated_allele("", 3, true), "");
    assert_eq!(shift_rotated_allele("GAT", i128::MIN, true), "ATG");
    assert_eq!(shift_rotated_allele("GAT", i128::MIN, false), "TGA");
}

#[test]
fn the_class_takes_a_reverse_strand_insertion_shifted_past_its_length() {
    let shifted = row("-3", &[]);
    assert!(taken("GA", &shifted, &shifted));
    // A positive offset of the same size, and a longer allele one or two bases past a turn.
    for (allele, offset) in [("GA", "3"), ("GAT", "-4"), ("GAT", "-5")] {
        let rows = row(offset, &[]);
        assert!(taken(allele, &rows, &rows), "{allele} {offset}");
    }
    // A shift within the allele's length or a multiple of it leaves the two readings equal,
    // as does an allele of one base repeated; an allele of one base and a deletion are not
    // two or more bases.
    for (allele, offset) in [
        ("GA", "-1"),
        ("GA", "-2"),
        ("GA", "-4"),
        ("GAT", "-2"),
        ("GAT", "-6"),
        ("GG", "-3"),
        ("G", "-3"),
        ("-", "-3"),
        ("GN", "-3"),
    ] {
        let rows = row(offset, &[]);
        assert!(!taken(allele, &rows, &rows), "{allele} {offset}");
    }
}

#[test]
fn the_class_reads_one_coding_change_from_both_rows() {
    let shifted = row("-3", &[]);
    // A run that does not compare the shift carries none on its rows, and names nothing;
    // neither does a zero, a text that is no integer or one 64 bits do not hold.
    for offset in ["-", "", "0", "-0", "3.0", "99999999999999999999"] {
        let rows = row(offset, &[]);
        assert!(!taken("GA", &rows, &rows), "{offset}");
    }
    // Rows that differ in their shift, a position or their `HGVSc` state two coding
    // changes, and rows without an `HGVSc` none.
    for change in [
        ("HGVS_OFFSET", "-5"),
        ("cDNA_position", "1201-1202"),
        ("CDS_position", "1091-1092"),
        ("HGVSc", "ENST1.1:c.1093_1094dup"),
    ] {
        assert!(!taken("GA", &shifted, &row("-3", &[change])), "{change:?}");
    }
    for unnamed in ["-", ""] {
        let rows = row("-3", &[("HGVSc", unnamed)]);
        assert!(!taken("GA", &rows, &rows), "{unnamed:?}");
    }
    assert!(!unrotated_shifted_insertion(
        "HGVSp",
        &key_of("GA", "Transcript"),
        REFERENCE,
        CANDIDATE,
        "-1",
        None,
        None
    ));
}

#[test]
fn the_class_takes_an_hgvsp_of_one_protein_on_a_reverse_strand_transcript_alone() {
    let shifted = row("-3", &[]);
    let rows = (Some(&shifted), Some(&shifted));
    let transcript = key_of("GA", "Transcript");
    let asked = |name: &str, key: &FieldKey, vep: &str, vep_rs: &str, strand: &str| {
        unrotated_shifted_insertion(name, key, vep, vep_rs, strand, rows.0, rows.1)
    };
    assert!(asked("HGVSp", &transcript, REFERENCE, CANDIDATE, "-1"));
    // Another field, another strand or none, and a row of another feature type.
    assert!(!asked("HGVSc", &transcript, REFERENCE, CANDIDATE, "-1"));
    for strand in ["1", "-", ""] {
        assert!(!asked("HGVSp", &transcript, REFERENCE, CANDIDATE, strand));
    }
    let regulatory = key_of("GA", "RegulatoryFeature");
    assert!(!asked("HGVSp", &regulatory, REFERENCE, CANDIDATE, "-1"));
    // Two protein names, a side without a notation, and a value without a name or a
    // description.
    for (vep, vep_rs) in [
        (REFERENCE, "ENSP2.1:p.Trp365SerfsTer42"),
        (REFERENCE, "-"),
        ("-", CANDIDATE),
        (REFERENCE, ""),
        (":p.Trp365PhefsTer42", ":p.Trp365SerfsTer42"),
        ("ENSP1.1:p.", CANDIDATE),
        (REFERENCE, "ENSP1.1:p."),
    ] {
        assert!(
            !asked("HGVSp", &transcript, vep, vep_rs, "-1"),
            "{vep} {vep_rs}"
        );
    }
}
