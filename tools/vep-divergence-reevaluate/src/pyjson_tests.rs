// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::testing::cases;

#[test]
fn a_document_is_read_and_written_as_cpython_json_does() {
    for case in cases("json") {
        let text = case.get_str("text").expect("a document");
        let value = loads(text);
        if case.has("refused") {
            assert_eq!(value, None, "json.loads({text:?})");
            continue;
        }
        let value = value.unwrap_or_else(|| panic!("json.loads({text:?}) refused"));
        let expect = |key: &str| case.get_str(key).expect("an expected text");
        assert_eq!(value.py_repr(), expect("repr"), "repr of {text}");
        assert_eq!(value.py_str(), expect("str"), "str of {text}");
        assert_eq!(dumps_compact_sorted(&value), expect("compact"), "{text}");
        assert_eq!(dumps_indented(&value, 1), expect("indented"), "{text}");
        assert_eq!(
            Value::Bool(value.truthy()),
            *case.get("truthy").unwrap(),
            "truth of {text}"
        );
    }
}

#[test]
fn half_a_surrogate_pair_reads_as_the_replacement_character() {
    assert_eq!(
        loads(r#""\ud83d""#),
        Some(Value::Str("\u{fffd}".to_string()))
    );
    assert_eq!(
        loads(r#""\ud83dx\ude00""#),
        Some(Value::Str("\u{fffd}x\u{fffd}".to_string()))
    );
}
