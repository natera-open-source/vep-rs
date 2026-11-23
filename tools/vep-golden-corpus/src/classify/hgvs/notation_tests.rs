// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the HGVSc notations: the order of cDNA coordinates, the two body patterns
//! and the counts read off a notation, held to the tables of
//! `tests/fixtures/hgvs_notations/`.

use super::*;

/// The rows of a table of `tests/fixtures/hgvs_notations/` after its header, by column.
fn table(name: &str) -> Vec<Vec<String>> {
    let path = format!(
        "{}/tests/fixtures/hgvs_notations/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{path}: {error}"));
    let rows: Vec<Vec<String>> = text
        .strip_suffix('\n')
        .expect("a table ends with a line end")
        .split('\n')
        .skip(1)
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    assert!(!rows.is_empty(), "{path} holds no row");
    rows
}

/// A cell of a table: `None` for the text `None`.
fn cell(text: &str) -> Option<&str> {
    (text != "None").then_some(text)
}

/// A cell that holds a number.
fn number(text: &str) -> i128 {
    text.parse()
        .unwrap_or_else(|_| panic!("{text} is no number"))
}

#[test]
fn cdna_coordinate_order_walks_the_transcript() {
    // The 5' UTR counts down to c.-1, the coding sequence up, the 3' UTR up from c.*1; an
    // intron offset follows its exon base, `-M` before it and `+M` after; `*-M`, the stop
    // codon's last base plus a negative intron distance printed without that base's
    // coordinate, sits after every coding coordinate and before `*1`.
    let coordinates = [
        "-5-2", "-5", "-5+2", "-1", "1", "3-1", "3", "3+1", "300", "*-3", "*1", "*12-3", "*12",
    ];
    let mut sorted = coordinates;
    sorted.reverse();
    sorted.sort_by_key(|coordinate| cdna_coordinate_order(coordinate).expect("a coordinate"));
    assert_eq!(sorted, coordinates);
    assert_eq!(cdna_coordinate_order("*-3"), Some((1, 0, -3)));
    assert_eq!(cdna_coordinate_order("-228+19"), Some((-1, -228, 19)));
}

#[test]
fn coordinates_are_placed_as_the_table_places_them() {
    let rows = table("coordinates.tsv");
    assert_eq!(rows.len(), 41);
    for row in rows {
        let coordinate = row[0].as_str();
        let order = (
            i8::try_from(number(&row[1])).unwrap(),
            number(&row[2]),
            number(&row[3]),
        );
        assert_eq!(
            cdna_coordinate_order(coordinate),
            Some(order),
            "{coordinate}"
        );
        assert_eq!(
            cdna_coordinate_index(coordinate),
            Some(number(&row[4])),
            "{coordinate}"
        );
        assert_eq!(
            intron_offset_text(coordinate).as_deref(),
            Some(row[3].as_str()),
            "{coordinate}"
        );
    }
}

#[test]
fn a_coordinate_is_read_as_int_reads_its_numbers() {
    // A number that is absent reads as 0, a sign with digits at the end as the intron
    // offset, before one line end as before none, and the rest as `int()` reads a text.
    let placed = [
        ("", Some((0, 0, 0))),
        ("*", Some((1, 0, 0))),
        ("-", Some((-1, 0, 0))),
        ("+5", Some((0, 0, 5))),
        ("--5", Some((-1, 0, -5))),
        ("5\n", Some((0, 5, 0))),
        ("5+3\n", Some((0, 5, 3))),
        (" 5 ", Some((0, 5, 0))),
        ("1_0+2", Some((0, 10, 2))),
        ("x", None),
        ("5x", None),
        ("5+", None),
        ("*x-3", None),
        ("5+3\n\n", None),
        ("5-+3", None),
    ];
    for (text, expected) in placed {
        assert_eq!(cdna_coordinate_order(text), expected, "{text:?}");
    }
    // A number 64 bits do not hold is no coordinate, and a body that carries one no notation.
    assert_eq!(cdna_coordinate_order("9223372036854775808"), None);
    assert_eq!(clipped_body("c.9223372036854775808del"), None);
    assert_eq!(indel_body("c.9223372036854775808del"), None);
    assert!(clipped_body("c.9223372036854775807del").is_some());
    assert!(indel_body("c.9223372036854775807del").is_some());
}

#[test]
fn clipped_body_is_the_comparators_clipped_notation_pattern() {
    // The pattern is the one the format comparator reads a clipped notation with, so a row
    // that comparator masks under the clipped class is one `classify` names; a change to
    // either alphabet is made in both. The table holds what that pattern reads of each body.
    assert_eq!(
        HGVSC_CLIPPED_BODY,
        concat!(
            r"(?P<numbering>[cn])\.(?P<start>\*?-?\d+(?:[+-]\d+)?)(?:_(?P<end>\*?-?\d+(?:[+-]\d+)?))?",
            r"(?:[ACGTN]>(?P<alternate_base>[ACGT])|(?P<kind>delins|ins)(?P<bases>[ACGT]+)|(?P<plain>inv|del|dup))",
        )
    );
    let rows = table("clipped_bodies.tsv");
    assert_eq!(rows.len(), 127);
    let (mut read, mut without_offsets) = (0, 0);
    for row in &rows {
        let body = row[0].replace("\\n", "\n");
        let expected = cell(&row[1]).map(|numbering| ClippedBody {
            numbering,
            start: &row[2],
            end: cell(&row[3]),
            alternate_base: cell(&row[4]),
            kind: cell(&row[5]),
            bases: cell(&row[6]),
            plain: cell(&row[7]),
        });
        assert_eq!(clipped_body(&body), expected, "{body:?}");
        let offsets = expected.as_ref().and_then(notation_intron_offsets);
        let expected_offsets = cell(&row[8]).map(|start| (start.to_string(), row[9].clone()));
        assert_eq!(offsets, expected_offsets, "{body:?}");
        read += usize::from(expected.is_some());
        without_offsets += usize::from(expected.is_some() && offsets.is_none());
    }
    assert_eq!((read, without_offsets), (111, 1));
}

#[test]
fn indel_bodies_are_read_as_the_table_reads_them() {
    let rows = table("indel_bodies.tsv");
    let mut read = 0;
    for row in &rows {
        let body = row[0].replace("\\n", "\n");
        let expected = cell(&row[1]).map(|numbering| IndelBody {
            numbering,
            start: &row[2],
            end: cell(&row[3]),
            kind: &row[4],
        });
        assert_eq!(indel_body(&body), expected, "{body:?}");
        read += usize::from(expected.is_some());
    }
    assert_eq!((rows.len(), read), (127, 76));
}

#[test]
fn exonic_run_counts_a_notations_bases_as_the_comparator_does() {
    // The bases from one cDNA coordinate to another, none where a coordinate carries an
    // intron offset or the two lie in different regions of the transcript. The table holds
    // every pair of fourteen coordinates as the format comparator counts it, so a notation
    // that comparator counts for the class of added bases outside the transcript is one
    // `classify` counts.
    assert_eq!(exonic_run("-103", "-95"), Some(9));
    assert_eq!(exonic_run("14", "14"), Some(1));
    assert_eq!(exonic_run("*5", "*9"), Some(5));
    assert_eq!(exonic_run("20+6", "20+6"), None);
    assert_eq!(exonic_run("-3", "2"), None);
    assert_eq!(exonic_run("20", "*1"), None);
    let rows = table("exonic_run.tsv");
    assert_eq!(rows.len(), 14 * 14);
    for row in rows {
        let expected = cell(&row[2]).map(number);
        assert_eq!(
            exonic_run(&row[0], &row[1]),
            expected,
            "{} {}",
            row[0],
            row[1]
        );
    }
}

#[test]
fn a_notation_names_the_rows_transcript_with_or_without_a_version() {
    let named = |hgvsc: &str, feature: &str| names_transcript(hgvsc, feature).is_some();
    assert!(named("ENST1.1:n.15dup", "ENST1"));
    assert!(named("ENST1.1:n.15dup", "ENST1.4"));
    assert!(named("ENST1:n.15dup", "ENST1.1"));
    assert!(!named("ENST2.1:n.15dup", "ENST1"));
    assert!(!named("ENST1.1:n.16A>T", "ENST1"));
    assert!(!named("ENST1.1", "ENST1"));
    assert!(clipped_notation("ENST1.1:n.16A>T", "ENST1").is_some());
    assert!(clipped_notation("-", "-").is_none());
    assert!(clipped_notation("", "").is_none());
    assert!(clipped_notation("ENST2.1:n.16A>T", "ENST1").is_none());
    assert_eq!(notation_name("ENST1.1:n.15dup"), "ENST1.1");
    assert_eq!(notation_name("-"), "-");
}

#[test]
fn a_location_and_a_cdna_position_are_cut_at_their_first_separators() {
    assert_eq!(span_texts("21:100-106"), ("100", "106"));
    assert_eq!(span_texts("21:100"), ("100", ""));
    assert_eq!(span_texts("100-106"), ("100", "106"));
    assert_eq!(span_texts("HLA:A:5-7-9"), ("A:5", "7-9"));
    assert_eq!(span_positions("21:100-106"), Some(7));
    assert_eq!(span_positions("21:100"), Some(1));
    assert_eq!(span_positions("21:100-99"), Some(0));
    assert_eq!(span_positions("21:100-"), Some(1));
    assert_eq!(span_positions("21:-100"), None);
    assert_eq!(span_positions("21:a-100"), None);
    assert_eq!(span_positions("21:100-b"), None);
    assert_eq!(position_ends("?-6/2000"), ("?", "6"));
    assert_eq!(position_ends("8485-?"), ("8485", "?"));
    assert_eq!(position_ends("147"), ("147", ""));
    assert_eq!(position_ends("-"), ("", ""));
    assert_eq!(position_ends(""), ("", ""));
}

#[test]
fn a_rows_fields_read_empty_where_the_row_or_the_field_is_absent() {
    let row: FieldRow = [("HGVSc", "ENST1.1:n.15dup"), ("HGVS_OFFSET", "-")]
        .into_iter()
        .map(|(name, value)| (name.to_string(), value.to_string()))
        .collect();
    assert_eq!(field(Some(&row), "HGVSc"), "ENST1.1:n.15dup");
    assert_eq!(field(Some(&row), "STRAND"), "");
    assert_eq!(field(None, "HGVSc"), "");
    assert!(holds(Some(&row), "HGVS_OFFSET", "-"));
    assert!(!holds(Some(&row), "STRAND", ""));
    assert!(!holds(None, "HGVSc", ""));
    assert_eq!(intron_offset_pair(Some(&row)), ("", ""));
    assert!(blank("") && blank("-") && !blank("0") && !blank(" "));
    assert!(non_zero_integer("-7") && non_zero_integer("21"));
    assert!(!non_zero_integer("0") && !non_zero_integer("07") && !non_zero_integer("+7"));
    assert!(!non_zero_integer("7\n") && !non_zero_integer(""));
}
