// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cells of a line, what a column's rules read of each part, and the header's mapped
//! columns with every refusal of a header that cannot be placed.

use super::*;
use crate::map::read_field_map;

const MAP: &str = "\
fastvep_column\tpart\tvep_field\tsource\tkind\tnote\tlist_separator\treference_set
FV_A\t1\tA_ONE\tS\ttext\tn\t-\t-
FV_A\t3-4\tA_RANGE\tS\ttext\tn\t-\t-
FV_B\t1.score\tB_SCORE\tS\tnumber\tn\t-\t-
FV_B\t2.word\tB_WORD\tS\ttext\tn\t-\t-
FV_B\t2.score\tB_TWO\tS\tnumber\tn\t-\t-
FV_B\t3.word\tB_HALF\tS\ttext\tn\t-\t-
FV_B\t3\tB_WHOLE\tS\ttext\tn\t-\t-
FV_B\t70.word\tB_FAR\tS\ttext\tn\t-\t-
FV_B\t80-99999999999999999999\tB_REST\tS\ttext\tn\t-\t-
";

fn map(kept: Option<&[&str]>) -> FieldMap {
    let kept: Option<Vec<Vec<u8>>> =
        kept.map(|names| names.iter().map(|name| name.as_bytes().to_vec()).collect());
    read_field_map("map.tsv", MAP.as_bytes(), kept.as_deref()).unwrap()
}

fn header(names: &[&str]) -> Vec<Vec<u8>> {
    names.iter().map(|name| name.as_bytes().to_vec()).collect()
}

const BUILT_IN: [&str; 5] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
];

fn with(extra: &[&str]) -> Vec<Vec<u8>> {
    header(&[&BUILT_IN[..], extra].concat())
}

#[test]
fn a_line_ends_at_its_line_feed_and_any_other_carriage_return_is_encoded() {
    assert_eq!(&*line_text(b"a\tb\n"), b"a\tb");
    assert_eq!(&*line_text(b"a\tb\r\n"), b"a\tb");
    assert_eq!(&*line_text(b"a\tb"), b"a\tb");
    assert_eq!(&*line_text(b"a\tb\r"), b"a\tb");
    assert_eq!(&*line_text(b"a\r\tb\r\r\n"), b"a%0D\tb%0D");
    assert_eq!(&*line_text(b"\n"), b"");
    assert_eq!(&*line_text(b"\r\n"), b"");
    assert_eq!(&*line_text(b""), b"");
    assert_eq!(&*line_text(b"a\n\n"), b"a\n");
    assert_eq!(
        header_cells(b"#Uploaded_variation\tLocation\tFV_\rA\r\n"),
        [
            b"Uploaded_variation".to_vec(),
            b"Location".to_vec(),
            b"FV_%0DA".to_vec()
        ]
    );
    assert_eq!(header_cells(b"#\n"), [b"".to_vec()]);
    assert_eq!(header_cells(b"#a"), [b"a".to_vec()]);
    assert_eq!(header_cells(b"# a\t#b\n"), [b" a".to_vec(), b"#b".to_vec()]);
}

#[test]
fn a_part_is_read_whole_by_half_or_not_at_all() {
    let placed = place_columns(Some(&with(&["FV_A", "FV_B"])), &map(None), false).unwrap();
    let (a, b) = (&placed.columns[0], &placed.columns[1]);
    let halves = |word, score| Reading::Halves { word, score };
    for (index, expected) in [
        (0, Reading::Unread),
        (1, Reading::Whole),
        (2, Reading::Unread),
        (3, Reading::Whole),
        (4, Reading::Whole),
        (5, Reading::Unread),
        (64, Reading::Unread),
        (u64::MAX, Reading::Unread),
    ] {
        assert_eq!(a.reading(index), expected, "FV_A part {index}");
    }
    for (index, expected) in [
        (1, halves(false, true)),
        (2, halves(true, true)),
        (3, Reading::Whole),
        (4, Reading::Unread),
        (63, Reading::Unread),
        (69, Reading::Unread),
        (70, halves(true, false)),
        (79, Reading::Unread),
        (80, Reading::Whole),
        (1 << 40, Reading::Whole),
        (u64::MAX, Reading::Whole),
    ] {
        assert_eq!(b.reading(index), expected, "FV_B part {index}");
    }
}

#[test]
fn the_mapped_columns_are_placed_with_their_kept_fields_in_column_order() {
    let placed = place_columns(
        Some(&with(&["FV_PHYLOP", "FV_B", "Other", "FV_A", "FV_PHYLOP"])),
        &map(None),
        false,
    )
    .unwrap();
    let places: Vec<(usize, &str)> = placed
        .columns
        .iter()
        .map(|column| (column.index, column.name.as_str()))
        .collect();
    assert_eq!(places, [(6, "FV_B"), (8, "FV_A")]);
    assert_eq!(
        placed.fields,
        [
            "B_SCORE", "B_WORD", "B_TWO", "B_HALF", "B_WHOLE", "B_FAR", "B_REST", "A_ONE",
            "A_RANGE"
        ]
    );
    assert_eq!(placed.columns[1].slots, [Some(7), Some(8)]);
    assert_eq!(
        placed.unmapped,
        [b"FV_PHYLOP".to_vec(), b"FV_PHYLOP".to_vec()]
    );
    assert_eq!(placed.width, 9);
    // A mapped column among the first five still leaves a row five cells wide.
    let placed = place_columns(
        Some(&header(&["FV_A", "Location", "Allele", "Gene", "Feature"])),
        &map(None),
        false,
    )
    .unwrap();
    assert_eq!((placed.columns[0].index, placed.width), (0, 5));
    // A run that names its fields keeps those alone, and their columns must be there.
    let kept = map(Some(&["A_RANGE", "B_TWO"]));
    let placed = place_columns(Some(&with(&["FV_A", "FV_B"])), &kept, true).unwrap();
    assert_eq!(placed.fields, ["A_RANGE", "B_TWO"]);
    assert_eq!(placed.columns[0].slots, [None, Some(0)]);
    assert_eq!(
        placed.columns[1].slots,
        [None, None, Some(1), None, None, None, None]
    );
}

#[test]
fn a_header_that_cannot_be_placed_is_refused() {
    let refusal = |header: Option<&[Vec<u8>]>, map: &FieldMap, named: bool| {
        place_columns(header, map, named).unwrap_err()
    };
    let all = map(None);
    assert_eq!(
        refusal(None, &all, false),
        "the input has no column header line"
    );
    assert_eq!(
        refusal(Some(&with(&["FV_PHYLOP"])), &all, false),
        "the field map names none of the header's columns: Uploaded_variation, Location, \
         Allele, Gene, Feature, FV_PHYLOP"
    );
    assert_eq!(
        refusal(Some(&[b"".to_vec()]), &all, false),
        "the field map names none of the header's columns: "
    );
    assert_eq!(
        refusal(Some(&[b"FV_\xffA".to_vec(), b"x".to_vec()]), &all, false),
        "the field map names none of the header's columns: FV_\\udcffA, x"
    );
    assert_eq!(
        refusal(
            Some(&with(&[
                "FV_B",
                "FV_A",
                "FV_B",
                "FV_PHYLOP",
                "FV_A",
                "FV_B"
            ])),
            &all,
            false
        ),
        "the header names a mapped column more than once: FV_B (columns 6, 8 and 11), FV_A \
         (columns 7 and 10)"
    );
    // The repeated column is refused before the absent column of a named field.
    let kept = map(Some(&["A_ONE", "B_WORD", "B_FAR"]));
    assert_eq!(
        refusal(Some(&with(&["FV_A", "FV_A"])), &kept, true),
        "the header names a mapped column more than once: FV_A (columns 6 and 7)"
    );
    assert_eq!(
        refusal(Some(&with(&["FV_A"])), &kept, true),
        "the header lacks the column of a field --fields names: B_FAR (FV_B), B_WORD (FV_B)"
    );
    // A run that names no field asks for no column.
    assert!(place_columns(Some(&with(&["FV_A"])), &all, false).is_ok());
}
