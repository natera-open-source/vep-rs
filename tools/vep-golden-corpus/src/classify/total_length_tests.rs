// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Tests of the JSON end key class under `--total_length`.

use super::*;
use crate::json::loads;

fn run(text: &str) -> Object {
    loads(text)
        .expect("the text is JSON")
        .as_object()
        .expect("an object")
        .clone()
}

fn flags(texts: &[&str]) -> Vec<String> {
    texts.iter().map(|text| text.to_string()).collect()
}

#[test]
fn a_range_with_a_total_gives_both_engines_end_values() {
    let pair = |reference: &str, candidate: &str| Some((reference.into(), candidate.into()));
    assert_eq!(range_end("484-487/2674"), pair("484", "487/2674"));
    assert_eq!(range_end("39-40/533"), pair("39", "40/533"));
    // A range whose first position is unknown: Ensembl VEP writes no end key.
    assert_eq!(range_end("?-487/2674"), pair("-", "487/2674"));
}

#[test]
fn any_other_column_text_gives_none() {
    for column in [
        "484/2674",
        "484-487",
        "484-?/2674",
        "?-487",
        "-/2674",
        "-",
        "",
        "484-487/",
        "484-/2674",
        "x-487/2674",
        "484-487/x",
    ] {
        assert_eq!(range_end(column), None, "{column}");
    }
}

#[test]
fn a_run_documents_the_class_under_the_flag_where_it_compares_json() {
    let total = flags(&["--hgvs", "--total_length"]);
    assert!(documents_range_ends(&total, &run(r#"{"name": "a"}"#)));
    assert!(documents_range_ends(
        &total,
        &run(r#"{"formats": ["default", "json"]}"#)
    ));
    assert!(!documents_range_ends(
        &total,
        &run(r#"{"formats": ["default", "tab", "vcf"]}"#)
    ));
    assert!(!documents_range_ends(&flags(&["--hgvs"]), &run("{}")));
    assert!(!documents_range_ends(&[], &run("{}")));
}
