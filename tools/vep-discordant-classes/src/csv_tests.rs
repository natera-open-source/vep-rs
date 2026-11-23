// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reader and the writer held to what `csv.reader` and `csv.writer` make of the same
//! text.

use super::*;

/// The records of `text`, fed in runs of `run` bytes.
fn records(text: &str, run: usize) -> Result<Vec<Vec<String>>, String> {
    let mut reader = Reader::new();
    let mut all = Vec::new();
    let mut sink = |record: Record<'_>| -> Result<(), String> {
        let fields = record.fields().count();
        assert_eq!(record.is_empty(), fields == 0);
        assert_eq!(record.get(fields), None);
        all.push(
            record
                .fields()
                .map(|field| String::from_utf8(field.to_vec()).unwrap())
                .collect(),
        );
        Ok(())
    };
    for piece in text.as_bytes().chunks(run.max(1)) {
        reader.feed(piece, &mut sink)?;
    }
    reader.finish(&mut sink)?;
    Ok(all)
}

/// Holds `text` to `expected` at every run length.
fn assert_reads(text: &str, expected: &[&[&str]]) {
    let expected: Vec<Vec<String>> = expected
        .iter()
        .map(|record| record.iter().map(|field| field.to_string()).collect())
        .collect();
    for run in 1..=text.len() + 1 {
        assert_eq!(
            records(text, run).unwrap(),
            expected,
            "{text:?} in runs of {run}"
        );
    }
}

#[test]
fn lines_without_a_quote_are_cut_at_their_tabs() {
    assert_reads("", &[]);
    assert_reads("a\tb\n", &[&["a", "b"]]);
    assert_reads("a\tb", &[&["a", "b"]]);
    assert_reads("a\t\n", &[&["a", ""]]);
    assert_reads("a\t", &[&["a", ""]]);
    assert_reads("\t\n", &[&["", ""]]);
    assert_reads(" \n", &[&[" "]]);
    assert_reads("a\0b\tc\n", &[&["a\0b", "c"]]);
    assert_reads(
        "a\x0bb\x0cc\x1cd\u{2028}e\u{85}f\tz\n",
        &[&["a\x0bb\x0cc\x1cd\u{2028}e\u{85}f", "z"]],
    );
}

#[test]
fn an_empty_line_is_a_record_of_no_field() {
    assert_reads("\n", &[&[]]);
    assert_reads("\n\n", &[&[], &[]]);
    assert_reads("a\n\nb\n", &[&["a"], &[], &["b"]]);
}

#[test]
fn a_field_that_opens_with_a_quote_runs_to_its_closing_quote() {
    assert_reads("\"a\tb\"\tc\n", &[&["a\tb", "c"]]);
    assert_reads("\"a\nb\"\tc\n", &[&["a\nb", "c"]]);
    assert_reads("\"a\"\"b\"\tc\n", &[&["a\"b", "c"]]);
    assert_reads("\"\"\n", &[&[""]]);
    assert_reads("\"\n\n\"\n", &[&["\n\n"]]);
    assert_reads("x\t\"\"\ty\n", &[&["x", "", "y"]]);
}

#[test]
fn text_behind_a_closing_quote_is_appended_to_the_field() {
    assert_reads("\"a\"b\tc\n", &[&["ab", "c"]]);
    assert_reads("\"a\" \tb\n", &[&["a ", "b"]]);
    assert_reads("\"a\"b\"c\"\td\n", &[&["ab\"c\"", "d"]]);
}

#[test]
fn a_quote_inside_a_field_is_a_character() {
    assert_reads("a\"b\"\tc\n", &[&["a\"b\"", "c"]]);
    assert_reads("a\t b\"\n", &[&["a", " b\""]]);
}

#[test]
fn a_quote_left_open_ends_with_the_input() {
    assert_reads("a\t\"bc", &[&["a", "bc"]]);
    assert_reads("a\t\"bc\n", &[&["a", "bc\n"]]);
    assert_reads("a\t\"", &[&["a", ""]]);
    assert_reads("a\t\"\"", &[&["a", ""]]);
    assert_reads("a\n\"b\tc\nd\te\n", &[&["a"], &["b\tc\nd\te\n"]]);
}

#[test]
fn a_record_behind_a_quoted_one_is_read_afresh() {
    assert_reads(
        "\"q\"\tr\ns\tt\n\"u\nv\"\nw\n",
        &[&["q", "r"], &["s", "t"], &["u\nv"], &["w"]],
    );
}

#[test]
fn a_field_holds_at_most_the_limit_in_code_points() {
    let lengths = |text: &str| -> Result<Vec<usize>, String> {
        Ok(records(text, 4096)?[0]
            .iter()
            .map(|field| field.chars().count())
            .collect())
    };
    let too_large = Err(FIELD_TOO_LARGE.to_string());
    for (unit, name) in [("x", "one byte"), ("\u{e9}", "two bytes")] {
        let at_limit = unit.repeat(FIELD_LIMIT);
        assert_eq!(
            lengths(&format!("{at_limit}\ty\n")),
            Ok(vec![FIELD_LIMIT, 1]),
            "{name}"
        );
        assert_eq!(
            lengths(&format!("{at_limit}{unit}\ty\n")),
            too_large,
            "{name}"
        );
        assert_eq!(
            lengths(&format!("\"{at_limit}\"\ty\n")),
            Ok(vec![FIELD_LIMIT, 1]),
            "{name}"
        );
        assert_eq!(
            lengths(&format!("\"{at_limit}{unit}\"\ty\n")),
            too_large,
            "{name}"
        );
        assert_eq!(
            lengths(&format!("y\t{at_limit}{unit}")),
            too_large,
            "{name}"
        );
    }
    let across_lines = format!("\"{}\n\"\ty\n", "x".repeat(FIELD_LIMIT - 1));
    assert_eq!(lengths(&across_lines), Ok(vec![FIELD_LIMIT, 1]));
    let across_lines = format!("\"{}\n\"\ty\n", "x".repeat(FIELD_LIMIT));
    assert_eq!(lengths(&across_lines), too_large);
}

#[test]
fn a_field_is_written_between_quotes_only_where_it_holds_a_tab_a_quote_or_a_line_ending() {
    let written = |field: &str| {
        let mut out = Vec::new();
        push_field(&mut out, field.as_bytes());
        String::from_utf8(out).unwrap()
    };
    for plain in [
        "",
        "a",
        "a b",
        " lead",
        "trail ",
        "a\0b",
        "a\x0bb\x0c\x1c\u{2028}\u{85}",
        "\u{e9}",
    ] {
        assert_eq!(written(plain), plain);
    }
    assert_eq!(written("a\tb"), "\"a\tb\"");
    assert_eq!(written("a\"b"), "\"a\"\"b\"");
    assert_eq!(written("a\nb"), "\"a\nb\"");
    assert_eq!(written("a\rb"), "\"a\rb\"");
    assert_eq!(written("\""), "\"\"\"\"");
    assert_eq!(written("\"\""), "\"\"\"\"\"\"");
}
