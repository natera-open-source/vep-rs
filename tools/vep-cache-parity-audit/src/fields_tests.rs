// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::{json_fields, Fields};
use crate::jsonpp::{decode, Value};

fn fields_of(text: &str) -> Fields {
    let Ok(Value::Dict(entries)) = decode(text.as_bytes()) else {
        panic!("an object: {text}");
    };
    json_fields(text.as_bytes(), &entries)
}

fn shown(fields: &Fields) -> Vec<(String, String)> {
    fields
        .iter()
        .map(|(path, kind)| {
            (
                String::from_utf8_lossy(path).into_owned(),
                String::from_utf8_lossy(kind).into_owned(),
            )
        })
        .collect()
}

fn pairs(expected: &[(&str, &str)]) -> Vec<(String, String)> {
    expected
        .iter()
        .map(|(path, kind)| ((*path).to_string(), (*kind).to_string()))
        .collect()
}

/// Each expectation is what the same walk in Perl gives for `decode_json` of the same text.
#[test]
fn a_json_transcript_is_its_paths_and_kinds() {
    let fields = fields_of(
        r#"{"stable_id":"ENST1","start":5,"flag":true,"none":null,"exons":[{"start":1,"end":2},{"start":3,"phase":null}],"cache":{"pairs":[[1,2],[3,4]],"empty":[],"object":{}}}"#,
    );
    assert_eq!(
        shown(&fields),
        pairs(&[
            ("cache", "HASH"),
            ("cache.empty", "ARRAY"),
            ("cache.object", "HASH"),
            ("cache.pairs", "ARRAY"),
            ("cache.pairs[]", "ARRAY"),
            ("cache.pairs[][]", "scalar"),
            ("exons", "ARRAY"),
            ("exons[]", "HASH"),
            ("exons[].end", "scalar"),
            ("exons[].phase", "undef"),
            ("exons[].start", "scalar"),
            ("flag", "JSON::PP::Boolean(SCALAR)"),
            ("none", "undef"),
            ("stable_id", "scalar"),
            ("start", "scalar"),
        ])
    );
}

#[test]
fn a_key_that_is_empty_or_zero_opens_no_path_of_its_own() {
    // `$prefix ? "$prefix.$k" : $k`: a false prefix is left out.
    let fields = fields_of(r#"{"0":{"a":1},"":{"b":2},"x":{"0":{"c":3}},"k":1,"k":[2]}"#);
    assert_eq!(
        shown(&fields),
        pairs(&[
            ("", "HASH"),
            ("0", "HASH"),
            ("a", "scalar"),
            ("b", "scalar"),
            ("k", "ARRAY"),
            ("k[]", "scalar"),
            ("x", "HASH"),
            ("x.0", "HASH"),
            ("x.0.c", "scalar"),
        ])
    );
}
