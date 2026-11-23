// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::record::parse_toml;
use crate::testing::{cases, texts};

#[test]
fn a_toml_value_is_written_as_cpython_writes_what_tomllib_reads() {
    for case in cases("toml") {
        let text = case.get_str("text").expect("a document");
        let document = parse_toml(text.as_bytes()).unwrap_or_else(|e| panic!("{text}: {e}"));
        let value = document.get("v").expect("the value");
        assert_eq!(
            value.py_repr(),
            case.get_str("repr").unwrap(),
            "repr of {text}"
        );
        assert_eq!(
            value.py_str(),
            case.get_str("str").unwrap(),
            "str of {text}"
        );
        assert_eq!(
            Value::Bool(value.truthy()),
            *case.get("truthy").unwrap(),
            "truth of {text}"
        );
        assert_eq!(
            as_list(value),
            texts(case.get("as_list").unwrap()),
            "as_list of {text}"
        );
    }
}

#[test]
fn a_table_is_read_by_key_with_python_defaults() {
    let record = parse_toml(b"flag = false\nzero = 0\n[run]\nstats = true\nname = ''\n").unwrap();
    assert!(record.has("flag") && !record.has("absent"));
    assert!(!record.truthy_at("flag") && !record.truthy_at("zero") && !record.truthy_at("absent"));
    assert!(record.truthy_or_true("absent") && !record.truthy_or_true("flag"));
    assert!(record.section("run").truthy_at("stats"));
    assert!(!record.section("run").truthy_at("name"));
    assert!(record.section("absent").items().is_empty());
    assert_eq!(record.section("run").get_str("name"), Some(""));
    assert_eq!(record.get_str("flag"), None);
    assert!(record.get("flag").unwrap().is_int() && record.get("zero").unwrap().is_int());
    assert_eq!(record.get("flag").unwrap().as_i128(), Some(0));
    assert!(!record.section("run").is_int());
}
