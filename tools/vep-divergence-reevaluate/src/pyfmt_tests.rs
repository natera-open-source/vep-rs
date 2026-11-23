// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::testing::{cases, texts};
use crate::value::Value;

fn bytes_of(hex: &str) -> Vec<u8> {
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).expect("hex digits"))
        .collect()
}

#[test]
fn a_str_is_written_stripped_split_and_quoted_as_cpython_does() {
    for case in cases("str") {
        let s = case.get_str("s").expect("a string");
        let expect = |key: &str| case.get_str(key).expect("an expected text");
        assert_eq!(py_repr_str(s), expect("repr"), "repr({s:?})");
        assert_eq!(py_strip(s), expect("strip"), "{s:?}.strip()");
        assert_eq!(py_rstrip(s), expect("rstrip"), "{s:?}.rstrip()");
        assert_eq!(py_lstrip(s), expect("lstrip"), "{s:?}.lstrip()");
        assert_eq!(shlex_quote(s), expect("quote"), "shlex.quote({s:?})");
        assert_eq!(
            py_splitlines(s),
            texts(case.get("splitlines").expect("lines")),
            "{s:?}.splitlines()"
        );
        if let Some(Value::Bool(space)) = case.get("isspace") {
            let c = s.chars().next().expect("one character");
            assert_eq!(is_py_space(c), *space, "{s:?}.isspace()");
        }
    }
}

#[test]
fn a_list_of_str_is_written_and_joined_as_cpython_does() {
    for case in cases("list") {
        let items = texts(case.get("items").expect("items"));
        assert_eq!(py_repr_str_list(&items), case.get_str("repr").unwrap());
        assert_eq!(shlex_join(&items), case.get_str("join").unwrap());
    }
}

#[test]
fn a_command_is_split_as_shlex_splits_it() {
    for case in cases("split") {
        let text = case.get_str("text").expect("a text");
        match case.get_str("error") {
            Some(error) => assert_eq!(shlex_split(text), Err(error), "shlex.split({text:?})"),
            None => assert_eq!(
                shlex_split(text),
                Ok(texts(case.get("words").expect("words"))),
                "shlex.split({text:?})"
            ),
        }
    }
}

#[test]
fn bytes_are_decoded_as_the_utf8_codec_decodes_them() {
    for case in cases("utf8") {
        let hex = case.get_str("hex").expect("hex");
        let bytes = bytes_of(hex);
        assert_eq!(
            decode_backslashreplace(&bytes),
            case.get_str("backslashreplace").unwrap(),
            "{hex}"
        );
        match case.get_str("error") {
            Some(error) => assert_eq!(decode_utf8_strict(&bytes), Err(error.to_string()), "{hex}"),
            None => assert_eq!(
                decode_utf8_strict(&bytes),
                Ok(case.get_str("strict").unwrap())
            ),
        }
    }
}

#[test]
fn the_last_line_of_a_text_is_its_last_after_stripping() {
    assert_eq!(
        last_line("WARNING: x\nERROR: no cache\n"),
        "ERROR: no cache"
    );
    assert_eq!(last_line("one\r\ntwo \n\n"), "two");
    assert_eq!(last_line(" \n\t"), "");
    assert_eq!(last_line(""), "");
}

#[test]
fn a_count_is_read_as_int_reads_it() {
    for (text, expected) in [
        ("20", 20),
        ("-0", 0),
        ("+2", 2),
        ("007", 7),
        ("1_000", 1000),
    ] {
        assert_eq!(py_int_clamped(text), Some(expected), "int({text:?})");
    }
    for (text, expected) in [(" 4 ", 4), ("\t3\n", 3), ("-3", -3)] {
        assert_eq!(py_int_clamped(text), Some(expected), "int({text:?})");
    }
    for text in [
        "", "x", "1.5", "0x10", "_1", "1_", "1__0", "- 3", "\u{ff13}",
    ] {
        assert_eq!(py_int_clamped(text), None, "int({text:?})");
    }
}

#[test]
fn a_count_is_written_as_format_groups_it() {
    for (n, expected) in [
        (0, "0"),
        (7, "7"),
        (999, "999"),
        (1000, "1,000"),
        (12345, "12,345"),
        (123456, "123,456"),
        (1234567, "1,234,567"),
    ] {
        assert_eq!(py_int_grouped(n), expected, "format({n}, \",\")");
    }
}

#[test]
fn a_text_mode_read_ends_lines_at_universal_newlines() {
    assert_eq!(universal_newlines("a\r\nb\rc\n"), "a\nb\nc\n");
    assert_eq!(lines_keepends("a\nb\n\nc"), ["a\n", "b\n", "\n", "c"]);
    assert!(lines_keepends("").is_empty());
}
