// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `dumps` under the three styles, held to the text `json.dumps` writes for the same value.

use super::*;

fn sample() -> Value {
    let mut numbers =
        loads("[1, -2, 1.5, 1e16, 1e-7, NaN, Infinity, -Infinity, true, false, null]").unwrap();
    let mut object = Object::new();
    object.insert("release", "116");
    object.insert(
        "covers",
        Value::strings(&["symbol", "\u{e9}\u{1F600}\u{7f}\u{1f}\"\\/"]),
    );
    object.insert("empty_list", Vec::new());
    object.insert("empty", Object::new());
    object.insert("n", std::mem::replace(&mut numbers, Value::Null));
    object.insert(
        "nested",
        loads(r#"{"a": {"b": [[], {}, [1, [2]]]}}"#).unwrap(),
    );
    Value::Object(object)
}

#[test]
fn the_plain_style_is_json_dumps() {
    assert_eq!(
        dumps(&sample(), &Style::PLAIN),
        concat!(
            "{\"release\": \"116\", \"covers\": [\"symbol\", ",
            "\"\\u00e9\\ud83d\\ude00\\u007f\\u001f\\\"\\\\/\"], \"empty_list\": [], \"empty\": {}, ",
            "\"n\": [1, -2, 1.5, 1e+16, 1e-07, NaN, Infinity, -Infinity, true, false, null], ",
            "\"nested\": {\"a\": {\"b\": [[], {}, [1, [2]]]}}}"
        )
    );
}

#[test]
fn the_manifest_style_is_json_dumps_with_an_indent_of_one() {
    assert_eq!(
        dumps(&sample(), &Style::INDENT_1),
        concat!(
            "{\n \"release\": \"116\",\n \"covers\": [\n  \"symbol\",\n  ",
            "\"\\u00e9\\ud83d\\ude00\\u007f\\u001f\\\"\\\\/\"\n ],\n \"empty_list\": [],\n ",
            "\"empty\": {},\n \"n\": [\n  1,\n  -2,\n  1.5,\n  1e+16,\n  1e-07,\n  NaN,\n  ",
            "Infinity,\n  -Infinity,\n  true,\n  false,\n  null\n ],\n \"nested\": {\n  ",
            "\"a\": {\n   \"b\": [\n    [],\n    {},\n    [\n     1,\n     [\n      2\n     ]\n    ]\n",
            "   ]\n  }\n }\n}"
        )
    );
    assert_eq!(dumps(&Value::List(Vec::new()), &Style::INDENT_1), "[]");
    assert_eq!(dumps(&Value::Object(Object::new()), &Style::INDENT_1), "{}");
    assert_eq!(dumps(&Value::from("x"), &Style::INDENT_1), "\"x\"");
    assert_eq!(
        dumps(&loads("[[]]").unwrap(), &Style::INDENT_1),
        "[\n []\n]"
    );
}

#[test]
fn the_compact_style_keeps_characters_outside_ascii() {
    assert_eq!(
        dumps(&sample(), &Style::COMPACT_UNICODE),
        concat!(
            "{\"release\":\"116\",\"covers\":[\"symbol\",",
            "\"\u{e9}\u{1F600}\u{7f}\\u001f\\\"\\\\/\"],\"empty_list\":[],\"empty\":{},",
            "\"n\":[1,-2,1.5,1e+16,1e-07,NaN,Infinity,-Infinity,true,false,null],",
            "\"nested\":{\"a\":{\"b\":[[],{},[1,[2]]]}}}"
        )
    );
}
