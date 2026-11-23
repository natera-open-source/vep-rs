// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison of the two caches' `info.json`: whether they are equal and, where they
//! are not, each key whose values differ with `str()` of both values cut to 200
//! characters.

use crate::json::{py_eq, py_ne, Object, Value};
use crate::pyfmt::py_str;
use crate::report::Json;

/// The characters of a value's `str()` the report keeps.
const SHOWN: usize = 200;

fn shown(value: Option<&Value>) -> Json {
    let mut text = py_str(value.unwrap_or(&Value::Null));
    text.truncate(SHOWN);
    Json::Text(text)
}

/// The `info_json` block of the report; `None` is a cache without the file.
pub fn compare(baseline: Option<&Object>, builder: Option<&Object>) -> Json {
    let (b, c) = match (baseline, builder) {
        (None, None) => {
            return Json::object([
                ("matches", Json::Bool(true)),
                ("diffs", Json::List(Vec::new())),
                ("note", Json::text("neither cache has info.json")),
            ])
        }
        (None, Some(_)) => return missing("info.json missing in baseline"),
        (Some(_), None) => return missing("info.json missing in builder"),
        (Some(b), Some(c)) => (b, c),
    };
    let same = b.len() == c.len()
        && b.iter()
            .all(|(key, value)| c.get(key).is_some_and(|other| py_eq(value, other)));
    if same {
        return Json::object([
            ("matches", Json::Bool(true)),
            ("diffs", Json::List(Vec::new())),
        ]);
    }
    let mut keys: Vec<&[u8]> = b.iter().chain(c.iter()).map(|(key, _)| key).collect();
    keys.sort_unstable();
    keys.dedup();
    let diffs = keys
        .into_iter()
        .filter_map(|key| {
            let (bval, cval) = (b.get(key), c.get(key));
            let differs = py_ne(bval.unwrap_or(&Value::Null), cval.unwrap_or(&Value::Null));
            differs.then(|| {
                Json::object([
                    ("field", Json::Text(crate::json::code_points(key))),
                    ("baseline", shown(bval)),
                    ("builder", shown(cval)),
                ])
            })
        })
        .collect();
    Json::object([("matches", Json::Bool(false)), ("diffs", Json::List(diffs))])
}

fn missing(side: &str) -> Json {
    Json::object([
        ("matches", Json::Bool(false)),
        ("diffs", Json::List(vec![Json::text(side)])),
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::load;
    use crate::report::dumps;

    fn object(text: &str) -> Object {
        match load(text.as_bytes()).unwrap() {
            Value::Object(object) => object,
            other => panic!("{other:?}"),
        }
    }

    fn compared(baseline: &str, builder: &str) -> String {
        dumps(&compare(Some(&object(baseline)), Some(&object(builder)))).replace("\n", "")
    }

    #[test]
    fn a_missing_file_is_named() {
        let info = object("{}");
        let line = |b: Option<&Object>, c: Option<&Object>| dumps(&compare(b, c)).replace("\n", "");
        assert_eq!(
            line(None, None),
            r#"{  "matches": true,  "diffs": [],  "note": "neither cache has info.json"}"#
        );
        assert_eq!(
            line(None, Some(&info)),
            r#"{  "matches": false,  "diffs": [    "info.json missing in baseline"  ]}"#
        );
        assert_eq!(
            line(Some(&info), None),
            r#"{  "matches": false,  "diffs": [    "info.json missing in builder"  ]}"#
        );
    }

    #[test]
    fn equal_documents_match_in_any_key_order() {
        let matched = r#"{  "matches": true,  "diffs": []}"#;
        assert_eq!(compared("{}", "{}"), matched);
        assert_eq!(
            compared(
                r#"{"a": 1, "b": [1, {"c": NaN}]}"#,
                r#"{"b": [1.0, {"c": NaN}], "a": true}"#
            ),
            matched
        );
    }

    #[test]
    fn differing_keys_are_listed_in_sorted_order_with_both_values() {
        assert_eq!(
            compared(
                r#"{"species": "homo_sapiens", "b": 1, "only_b": [1, "x"], "same": {"k": 2}}"#,
                r#"{"same": {"k": 2.0}, "only_c": null, "species": "mus", "b": "1"}"#,
            ),
            concat!(
                r#"{  "matches": false,  "diffs": [    {"#,
                r#"      "field": "b",      "baseline": "1",      "builder": "1"    },    {"#,
                r#"      "field": "only_b",      "baseline": "[1, 'x']",      "builder": "None"    },    {"#,
                r#"      "field": "species",      "baseline": "homo_sapiens",      "builder": "mus"    }  ]}"#,
            )
        );
        // A key one side lacks and the other holds as null is no difference, though the
        // documents differ; a NaN both hold is one.
        assert_eq!(
            compared(r#"{"a": null}"#, "{}"),
            r#"{  "matches": false,  "diffs": []}"#
        );
        assert_eq!(
            compared(r#"{"a": NaN, "b": 1}"#, r#"{"a": NaN, "b": 2}"#),
            concat!(
                r#"{  "matches": false,  "diffs": [    {"#,
                r#"      "field": "a",      "baseline": "nan",      "builder": "nan"    },    {"#,
                r#"      "field": "b",      "baseline": "1",      "builder": "2"    }  ]}"#,
            )
        );
    }

    #[test]
    fn a_value_is_shown_to_200_characters() {
        let long = "x".repeat(300);
        let text = compared(
            &format!(r#"{{"k": "{long}"}}"#),
            &format!(r#"{{"k": ["{long}"]}}"#),
        );
        let list = format!("['{}", "x".repeat(198));
        assert_eq!(
            text,
            format!(
                r#"{{  "matches": false,  "diffs": [    {{      "field": "k",      "baseline": "{}",      "builder": "{list}"    }}  ]}}"#,
                "x".repeat(200)
            )
        );
    }
}
