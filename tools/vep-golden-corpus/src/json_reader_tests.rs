// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `loads` on the texts `json.loads` reads and on the ones it refuses: each value is held to
//! the text `json.dumps` writes for it and each refusal to the text of the
//! `json.JSONDecodeError`.

use super::*;

fn read(text: &str) -> Result<String, String> {
    loads(text)
        .map(|value| dumps(&value, &Style::PLAIN))
        .map_err(|error| error.to_string())
}

#[test]
fn values_read_as_json_loads_reads_them() {
    let cases = [
        ("NaN", "NaN"),
        ("Infinity", "Infinity"),
        ("-Infinity", "-Infinity"),
        (
            "[1e5, 1E-2, 1.5e+3, -0, -0.0, 0e0, 1e400, -1e400, 1.0, 100, 1e16, 1.50]",
            "[100000.0, 0.01, 1500.0, 0, -0.0, 0.0, Infinity, -Infinity, 1.0, 100, 1e+16, 1.5]",
        ),
        (r#"{"a":1,"b":2,"a":3}"#, r#"{"a": 3, "b": 2}"#),
        ("\"\\ud83d\\ude00\"", "\"\\ud83d\\ude00\""),
        (
            "\"\\u00e9\\n\\/\\b\\f\\r\\t\\\"\\\\\"",
            "\"\\u00e9\\n/\\b\\f\\r\\t\\\"\\\\\"",
        ),
        (" \t\r\n [ \t\r\n ] \t\r\n ", "[]"),
        (
            r#"{"a":{"b":[{"c":null}]}}"#,
            r#"{"a": {"b": [{"c": null}]}}"#,
        ),
        ("\"\u{7f}\"", r#""\u007f""#),
        ("[true,false,null]", "[true, false, null]"),
        (
            "123456789012345678901234567890",
            "123456789012345678901234567890",
        ),
        (
            "-123456789012345678901234567890",
            "-123456789012345678901234567890",
        ),
        ("\"\\u0041\\u00e9\"", "\"A\\u00e9\""),
        (r#"{"a":"b" , "c" : [ ] }"#, r#"{"a": "b", "c": []}"#),
        ("\"\\u1234\"", "\"\\u1234\""),
    ];
    for (text, expected) in cases {
        assert_eq!(read(text), Ok(expected.to_string()), "{text:?}");
    }
}

#[test]
fn refused_texts_carry_the_decode_error_of_json_loads() {
    let cases = [
        ("", "Expecting value: line 1 column 1 (char 0)"),
        (" ", "Expecting value: line 1 column 2 (char 1)"),
        (
            "{",
            "Expecting property name enclosed in double quotes: line 1 column 2 (char 1)",
        ),
        ("[", "Expecting value: line 1 column 2 (char 1)"),
        ("[1,", "Expecting value: line 1 column 4 (char 3)"),
        (
            "[1,]",
            "Illegal trailing comma before end of array: line 1 column 3 (char 2)",
        ),
        (
            r#"{"a":1,}"#,
            "Illegal trailing comma before end of object: line 1 column 7 (char 6)",
        ),
        (
            r#"{"a""#,
            "Expecting ':' delimiter: line 1 column 5 (char 4)",
        ),
        (r#"{"a":"#, "Expecting value: line 1 column 6 (char 5)"),
        (
            r#"{"a":1"#,
            "Expecting ',' delimiter: line 1 column 7 (char 6)",
        ),
        (
            r#"{"a" 1}"#,
            "Expecting ':' delimiter: line 1 column 6 (char 5)",
        ),
        (
            "{a:1}",
            "Expecting property name enclosed in double quotes: line 1 column 2 (char 1)",
        ),
        ("[1 2]", "Expecting ',' delimiter: line 1 column 4 (char 3)"),
        (
            r#"{"a":1 "b":2}"#,
            "Expecting ',' delimiter: line 1 column 8 (char 7)",
        ),
        (
            "\"abc",
            "Unterminated string starting at: line 1 column 1 (char 0)",
        ),
        (
            "\"a\\",
            "Unterminated string starting at: line 1 column 1 (char 0)",
        ),
        (r#""a\x""#, "Invalid \\escape: line 1 column 3 (char 2)"),
        (
            "\"\\u1234",
            "Unterminated string starting at: line 1 column 1 (char 0)",
        ),
        (
            "\"\\ud800\\udc00",
            "Unterminated string starting at: line 1 column 1 (char 0)",
        ),
        (
            "\"a\tb\"",
            "Invalid control character at: line 1 column 3 (char 2)",
        ),
        (
            "\"a\nb\"",
            "Invalid control character at: line 1 column 3 (char 2)",
        ),
        ("01", "Extra data: line 1 column 2 (char 1)"),
        ("-", "Expecting value: line 1 column 1 (char 0)"),
        ("-x", "Expecting value: line 1 column 1 (char 0)"),
        ("1.", "Extra data: line 1 column 2 (char 1)"),
        ("1e", "Extra data: line 1 column 2 (char 1)"),
        ("1e+", "Extra data: line 1 column 2 (char 1)"),
        (".5", "Expecting value: line 1 column 1 (char 0)"),
        ("+1", "Expecting value: line 1 column 1 (char 0)"),
        ("nul", "Expecting value: line 1 column 1 (char 0)"),
        ("-Inf", "Expecting value: line 1 column 1 (char 0)"),
        ("nan", "Expecting value: line 1 column 1 (char 0)"),
        ("1 2", "Extra data: line 1 column 3 (char 2)"),
        ("[]]", "Extra data: line 1 column 3 (char 2)"),
        (
            "\u{feff}{}",
            "Unexpected UTF-8 BOM (decode using utf-8-sig): line 1 column 1 (char 0)",
        ),
        ("{\"a\":1}\n\n x", "Extra data: line 3 column 2 (char 10)"),
        (
            "{\n \"a\": 1,\n \"b\": \n}",
            "Expecting value: line 4 column 1 (char 18)",
        ),
        (
            "[\n1,\n2,\n\n  ]",
            "Illegal trailing comma before end of array: line 3 column 2 (char 6)",
        ),
        (
            "\"\u{e9}\\q\"",
            "Invalid \\escape: line 1 column 3 (char 2)",
        ),
        (
            "[\"\u{e9}\", \u{e9}]",
            "Expecting value: line 1 column 7 (char 6)",
        ),
        (
            "\n\n  [1, 2,, 3]",
            "Expecting value: line 3 column 9 (char 10)",
        ),
        (
            r#"{"a":1,,}"#,
            "Expecting property name enclosed in double quotes: line 1 column 8 (char 7)",
        ),
        ("[,]", "Expecting value: line 1 column 2 (char 1)"),
        (
            "[1\u{0b}]",
            "Expecting ',' delimiter: line 1 column 3 (char 2)",
        ),
        ("\u{0b}1", "Expecting value: line 1 column 1 (char 0)"),
        (r#""\U0041""#, "Invalid \\escape: line 1 column 2 (char 1)"),
        ("1.5.5", "Extra data: line 1 column 4 (char 3)"),
        ("0x10", "Extra data: line 1 column 2 (char 1)"),
    ];
    for (text, expected) in cases {
        assert_eq!(read(text), Err(expected.to_string()), "{text:?}");
    }
}

#[test]
fn an_escape_without_four_hex_digits_is_refused_at_its_u() {
    let cases = [
        (r#""\u12""#, "line 1 column 3 (char 2)"),
        (r#""\u12G4""#, "line 1 column 3 (char 2)"),
        (r#""\u123"#, "line 1 column 3 (char 2)"),
        (r#""\ud800\u12""#, "line 1 column 9 (char 8)"),
        (r#""\ud800\uZZZZ" "#, "line 1 column 9 (char 8)"),
        (r#""\ud800\u"#, "line 1 column 9 (char 8)"),
        ("\"\\ud800\\u12\u{e9}4\"", "line 1 column 9 (char 8)"),
    ];
    for (text, place) in cases {
        let expected = format!("{BAD_UNICODE_ESCAPE}: {place}");
        assert_eq!(read(text), Err(expected), "{text:?}");
    }
}

#[test]
fn texts_outside_the_reader_are_named_after_any_decode_error() {
    assert_eq!(
        loads(r#"["a", "\ud800x"]"#),
        Err(JsonError::Unsupported(
            "a lone surrogate escape at char 7".to_string()
        ))
    );
    assert_eq!(
        read(r#"["\udc00\ud800", 1"#),
        Err("Expecting ',' delimiter: line 1 column 19 (char 18)".to_string())
    );
    let deep = format!("{}{}", "[".repeat(MAX_DEPTH), "]".repeat(MAX_DEPTH));
    assert!(loads(&deep).is_ok());
    assert_eq!(
        loads(&format!("[{deep}]")),
        Err(JsonError::Unsupported(
            "nesting past 1000 levels at char 1000".to_string()
        ))
    );
    assert!(loads(&"1".repeat(MAX_INT_DIGITS)).is_ok());
    assert_eq!(
        loads(&format!("-{}", "1".repeat(MAX_INT_DIGITS + 1))),
        Err(JsonError::Unsupported(
            "an integer of more than 4300 digits at char 0".to_string()
        ))
    );
}

#[test]
fn an_object_behaves_as_a_dict() {
    let mut object = Object::new();
    object.insert("b", 1i64);
    object.insert("a", "x");
    assert_eq!(object.insert("b", true), Some(Value::from(1i64)));
    assert_eq!(object.keys().collect::<Vec<_>>(), ["b", "a"]);
    assert_eq!(object.get("b"), Some(&Value::Bool(true)));
    assert_eq!(object.remove("b"), Some(Value::Bool(true)));
    assert_eq!(object.remove("b"), None);
    object.get_or_insert("list", Vec::new());
    object
        .get_or_insert("list", "unread")
        .as_list_mut()
        .unwrap()
        .push(Value::Null);
    assert_eq!(
        dumps_object(&object, &Style::PLAIN),
        r#"{"a": "x", "list": [null]}"#
    );
    let other: Object = [("list", Value::strings(&["y"])), ("z", Value::Null)]
        .into_iter()
        .collect();
    object.update(other);
    object.retain(|key, _| key != "a");
    assert_eq!(
        dumps_object(&object, &Style::PLAIN),
        r#"{"list": ["y"], "z": null}"#
    );
    // Two objects of the same entries are equal in any order of their keys.
    assert_eq!(
        loads(r#"{"a": [1, {"b": 2, "c": 3}]}"#),
        loads(r#"{"a": [1, {"c": 3, "b": 2}]}"#)
    );
    assert_ne!(loads(r#"{"a": 1}"#), loads(r#"{"a": 1.0}"#));
}

#[test]
fn an_object_of_many_keys_keeps_its_order_and_finds_each() {
    let mut object = Object::new();
    for number in 0..200usize {
        object.insert(format!("key{number}"), number);
    }
    object.insert("key7", "seven");
    assert_eq!(object.remove("key150"), Some(Value::from(150usize)));
    object.insert("key150", "back");
    assert_eq!(object.len(), 200);
    assert_eq!(object.get("key7").and_then(Value::as_str), Some("seven"));
    assert_eq!(object.get("key199").and_then(Value::as_i64), Some(199));
    let keys: Vec<&str> = object.keys().collect();
    assert_eq!(keys[..3], ["key0", "key1", "key2"]);
    assert_eq!(keys[149..151], ["key149", "key151"]);
    assert_eq!(keys[199], "key150");
    object.retain(|key, _| key.ends_with('0'));
    assert_eq!(object.len(), 20);
    assert!(object.contains_key("key190") && !object.contains_key("key191"));
}

#[test]
fn the_accessors_read_one_kind_each() {
    let value =
        loads(r#"{"s": "t", "n": 7, "big": 99999999999999999999, "f": 0.0, "l": [], "o": {}, "y": true, "z": null}"#)
            .unwrap();
    assert_eq!(value.get("s").and_then(Value::as_str), Some("t"));
    assert_eq!(value.get("n").and_then(Value::as_i64), Some(7));
    assert_eq!(value.get("big").and_then(Value::as_i64), None);
    assert_eq!(value.get("y").and_then(Value::as_bool), Some(true));
    assert_eq!(value.get("y").and_then(Value::as_i64), None);
    assert!(value.get("l").and_then(Value::as_list).is_some());
    assert!(value.get("o").and_then(Value::as_object).is_some());
    assert_eq!(value.get("absent"), None);
    assert_eq!(Value::Null.get("s"), None);
    let truths: Vec<bool> = ["s", "n", "big", "f", "l", "o", "y", "z"]
        .iter()
        .map(|key| value.get(key).unwrap().truthy())
        .collect();
    assert_eq!(truths, [true, true, true, false, false, false, true, false]);
    let Some(Value::Int(big)) = value.get("big") else {
        panic!("an integer past 64 bits stays an integer");
    };
    assert_eq!(big.text(), "99999999999999999999");
    assert!(!big.is_negative() && !big.is_zero());
}
