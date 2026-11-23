// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of a field map as rules, the parts a row may name, and the names `--fields`
//! keeps.

use super::*;

const HEADER: &str =
    "fastvep_column\tpart\tvep_field\tsource\tkind\tnote\tlist_separator\treference_set";

/// A map of one row per `(column, part, field, separator)`.
fn map_text(rows: &[(&str, &str, &str, &str)]) -> String {
    let mut text = format!("{HEADER}\n");
    for (column, part, field, separator) in rows {
        text.push_str(&format!(
            "{column}\t{part}\t{field}\tS\ttext\tnote\t{separator}\t-\n"
        ));
    }
    text
}

fn read(text: &str, kept: Option<&[Vec<u8>]>) -> Result<FieldMap, Failure> {
    read_field_map("map.tsv", text.as_bytes(), kept)
}

fn refusal(text: &str) -> String {
    match read(text, None) {
        Err(Failure::Refused(message)) => message,
        other => panic!("not a refusal: {other:?}"),
    }
}

/// The rule of a map of one row whose part is `part`.
fn rule_of(part: &str) -> Rule {
    let map = read(&map_text(&[("FV_X", part, "X_ONE", "-")]), None).unwrap();
    map.columns[0].1[0].clone()
}

#[test]
fn the_shipped_map_reads_as_one_rule_a_row_under_its_columns() {
    let shipped = include_str!("../data/fastvep_field_map.tsv");
    let map = read(shipped, None).unwrap();
    let columns: Vec<(&str, usize)> = map
        .columns
        .iter()
        .map(|(name, rules)| (name.as_str(), rules.len()))
        .collect();
    assert_eq!(
        columns,
        [
            ("FV_ALPHAMISSENSE", 2),
            ("FV_REVEL", 1),
            ("FV_GNOMAD_GENE", 2),
            ("SpliceAI", 1),
            ("FV_PRIMATEAI", 1),
            ("FV_DBNSFP", 2),
            ("FV_CLINVAR", 5),
            ("FV_GNOMAD", 14),
        ]
    );
    let splice = &map.rules(b"SpliceAI").unwrap()[0];
    assert_eq!(
        (splice.first, splice.last, splice.half),
        (1, 9, Half::Whole)
    );
    let sift = &map.rules(b"FV_DBNSFP").unwrap()[0];
    assert_eq!((sift.first, sift.last, sift.half), (1, 1, Half::Score));
    assert_eq!(sift.field, "SIFT_score");
    let significance = &map.rules(b"FV_CLINVAR").unwrap()[0];
    assert_eq!(significance.separator, "|");
    assert_eq!(map.rules(b"FV_CLINVAR").unwrap()[1].separator, "");
    assert!(map
        .columns
        .iter()
        .all(|(_, rules)| rules.iter().all(|rule| rule.selected)));
    assert!(map.rules(b"FV_PHYLOP").is_none());
}

#[test]
fn a_part_is_a_number_a_rising_range_or_one_half_in_the_digits_of_any_script() {
    for (part, first, last, half) in [
        ("1", 1, 1, Half::Whole),
        ("007", 7, 7, Half::Whole),
        ("2-5", 2, 5, Half::Whole),
        ("01-2", 1, 2, Half::Whole),
        ("3.word", 3, 3, Half::Word),
        ("4.score", 4, 4, Half::Score),
        ("\u{661}", 1, 1, Half::Whole),
        ("1-\u{662}", 1, 2, Half::Whole),
        ("\u{ff11}\u{ff12}", 12, 12, Half::Whole),
        ("18446744073709551615", u64::MAX, u64::MAX, Half::Whole),
        ("99999999999999999999999", u64::MAX, u64::MAX, Half::Whole),
        ("5-99999999999999999999999", 5, u64::MAX, Half::Whole),
        (
            "18446744073709551616-18446744073709551617",
            u64::MAX,
            u64::MAX,
            Half::Whole,
        ),
    ] {
        let rule = rule_of(part);
        assert_eq!(
            (rule.first, rule.last, rule.half),
            (first, last, half),
            "{part:?}"
        );
    }
}

#[test]
fn a_part_of_another_shape_is_refused_with_its_repr() {
    for (part, shown) in [
        ("first", "'first'"),
        ("", "''"),
        ("1.", "'1.'"),
        ("1.Word", "'1.Word'"),
        ("1.word.score", "'1.word.score'"),
        ("-1", "'-1'"),
        ("1-", "'1-'"),
        ("1-2-3", "'1-2-3'"),
        (" 1", "' 1'"),
        ("1 ", "'1 '"),
        ("1_0", "'1_0'"),
        ("\u{b2}", "'\u{b2}'"),
        ("'", "\"'\""),
        ("a\"b'c", "'a\"b\\'c'"),
        ("\u{7f}", "'\\x7f'"),
        ("\u{200b}1", "'\\u200b1'"),
    ] {
        assert_eq!(
            refusal(&map_text(&[("FV_X", part, "X_ONE", "-")])),
            format!("map.tsv:2: part {shown} is not N, N-M, N.word or N.score"),
            "{part:?}"
        );
    }
    for part in [
        "0",
        "00",
        "0-1",
        "3-2",
        "2-2",
        "2-02",
        "\u{662}-\u{661}",
        "99999999999999999999999-1",
        "18446744073709551617-18446744073709551616",
        "18446744073709551616-18446744073709551616",
    ] {
        assert_eq!(
            refusal(&map_text(&[("FV_X", part, "X_ONE", "-")])),
            format!(
                "map.tsv:2: part {} must count from 1, and a range must rise",
                py_repr_str(part)
            ),
            "{part:?}"
        );
    }
}

#[test]
fn a_part_past_the_digit_limit_of_int_ends_in_the_error_python_raises() {
    let limit = "9".repeat(4300);
    assert_eq!(rule_of(&limit).first, u64::MAX);
    let past = "0".repeat(4301);
    let line = "ValueError: Exceeds the limit (4300 digits) for integer string conversion: value \
                has 4301 digits; use sys.set_int_max_str_digits() to increase the limit";
    assert_eq!(
        read(&map_text(&[("FV_X", &past, "X_ONE", "-")]), None),
        Err(Failure::Uncaught(line.to_string()))
    );
    let range = format!("1-{past}");
    assert_eq!(
        read(&map_text(&[("FV_X", &range, "X_ONE", "-")]), None),
        Err(Failure::Uncaught(line.to_string()))
    );
}

#[test]
fn a_malformed_map_is_refused_with_its_line_in_the_order_the_checks_run() {
    assert_eq!(
        read_field_map("m", b"\xff\xfe\n", None),
        Err(Failure::Refused(
            "the field map m is not UTF-8 text".to_string()
        ))
    );
    let header = "map.tsv:1: the header is not fastvep_column part vep_field source kind note \
                  list_separator reference_set";
    for text in ["", "\n", "fastvep_column\tpart\n", "\u{feff}fastvep_column"] {
        assert_eq!(refusal(text), header, "{text:?}");
    }
    assert_eq!(refusal(&format!("\u{feff}{HEADER}\n")), header);
    assert_eq!(refusal(&format!("{HEADER}\t\n")), header);
    assert_eq!(
        refusal(&format!("{HEADER}\n\n")),
        "map.tsv:2: 1 columns, not 8"
    );
    assert_eq!(
        refusal(&format!("{HEADER}\nFV_X\t1\tX\tS\ttext\tn\t-\n")),
        "map.tsv:2: 7 columns, not 8"
    );
    assert_eq!(
        refusal(&format!("{HEADER}\nFV_X\t1\tX\tS\ttext\tn\t-\t-\tmore\n")),
        "map.tsv:2: 9 columns, not 8"
    );
    // A row that fails several checks is refused by the first: the part, the kind, an
    // empty cell, a field mapped twice, the separator.
    let row = |part: &str, field: &str, kind: &str, separator: &str| {
        format!(
            "{HEADER}\nFV_X\t1\tX\tS\ttext\tn\t-\t-\n\
             FV_X\t{part}\t{field}\tS\t{kind}\tn\t{separator}\t-\n"
        )
    };
    assert_eq!(
        refusal(&row("0", "X", "float", "")),
        "map.tsv:3: part '0' must count from 1, and a range must rise"
    );
    assert_eq!(
        refusal(&row("2", "X", "float", "")),
        "map.tsv:3: kind 'float' is not text or number"
    );
    assert_eq!(
        refusal(&row("2", "X", "text", "")),
        "map.tsv:3: an empty fastvep_column, vep_field, list_separator or reference_set"
    );
    assert_eq!(
        refusal(&row("2-3", "X", "text", "|")),
        "map.tsv:3: X is mapped twice"
    );
    assert_eq!(
        refusal(&row("2-3", "Y", "number", "|")),
        "map.tsv:3: a list separator needs a part written N, not '2-3'"
    );
    assert_eq!(
        refusal(&row("2.score", "Y", "number", ",")),
        "map.tsv:3: a list separator needs a part written N, not '2.score'"
    );
}

#[test]
fn a_map_line_ends_where_splitlines_ends_it() {
    let rows = [("FV_X", "1", "X_ONE", "-"), ("FV_Y", "2", "Y_TWO", ",")];
    let text = map_text(&rows);
    let expected = read(&text, None).unwrap();
    assert_eq!(expected.columns.len(), 2);
    for end in [
        "\r\n", "\r", "\u{0b}", "\u{0c}", "\u{1c}", "\u{1d}", "\u{1e}", "\u{85}", "\u{2028}",
    ] {
        assert_eq!(
            read(&text.replace('\n', end), None).unwrap(),
            expected,
            "{end:?}"
        );
    }
    assert_eq!(read(text.trim_end(), None).unwrap(), expected);
    // A line end inside a cell cuts the row in two.
    assert_eq!(
        refusal(&text.replace("text\tnote", "text\tno\u{85}te")),
        "map.tsv:2: 6 columns, not 8"
    );
}

#[test]
fn fields_keeps_the_names_it_lists_and_refuses_a_name_the_map_lacks() {
    let names = |value: &str| kept_fields(Some(value.as_bytes())).unwrap().unwrap();
    assert_eq!(kept_fields(None), Ok(None));
    assert_eq!(names("A"), [b"A".to_vec()]);
    assert_eq!(
        names("A,B,A"),
        [b"A".to_vec(), b"B".to_vec(), b"A".to_vec()]
    );
    for (value, message) in [
        ("", "--fields names no field"),
        (",", "--fields names no field beside a comma: ,"),
        ("A,,B", "--fields names no field beside a comma: A,,B"),
        ("A,", "--fields names no field beside a comma: A,"),
    ] {
        assert_eq!(
            kept_fields(Some(value.as_bytes())),
            Err(Failure::Refused(message.to_string())),
            "{value:?}"
        );
    }
    assert_eq!(
        kept_fields(Some(b"A,\xff,")),
        Err(Failure::Refused(
            "--fields names no field beside a comma: A,\\udcff,".to_string()
        ))
    );
    let text = map_text(&[("FV_X", "1", "X_ONE", "-"), ("FV_X", "2", "X_TWO", "-")]);
    let map = read(&text, Some(&names("X_TWO"))).unwrap();
    let selected: Vec<bool> = map.columns[0].1.iter().map(|rule| rule.selected).collect();
    assert_eq!(selected, [false, true]);
    // The names the map lacks are listed once each in code point order, a byte outside
    // UTF-8 as its lone surrogate: after U+D7FF and before U+E000.
    let kept = [
        b"zz".to_vec(),
        b"\xee\x80\x80".to_vec(),
        b"X_ONE".to_vec(),
        b"\xff".to_vec(),
        b"\xed\x9f\xbf".to_vec(),
        b"zz".to_vec(),
        b"Zeta".to_vec(),
    ];
    assert_eq!(
        read(&text, Some(&kept)),
        Err(Failure::Refused(
            "--fields names no row of map.tsv: Zeta, zz, \u{d7ff}, \\udcff, \u{e000}".to_string()
        ))
    );
}
