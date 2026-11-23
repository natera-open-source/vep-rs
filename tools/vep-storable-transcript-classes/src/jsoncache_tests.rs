// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Each expectation is what CPython 3.14's `json.loads` and the walk over its result
//! give for the same text.

use std::collections::HashSet;

use super::{scan, Event, Peptide};

fn wanted(ids: &[&str]) -> HashSet<String> {
    ids.iter().map(|id| (*id).to_string()).collect()
}

fn found(id: &str, peptide: Peptide) -> Event {
    Event::Found {
        id: id.to_string(),
        peptide,
    }
}

fn text(peptide: &str) -> Peptide {
    Peptide::Text(peptide.to_string())
}

#[test]
fn the_first_transcript_of_a_wanted_id_gives_its_peptide() {
    let shard = br#"{"1": [
        {"stable_id": "ENST1.4", "variation_effect_feature_cache": {"peptide": "MKG"}},
        {"stable_id": "ENST2", "variation_effect_feature_cache": {"peptide": "MK*"}},
        {"stable_id": "ENST1.9", "variation_effect_feature_cache": {"peptide": "other"}},
        {"stable_id": "ENST3"},
        {"stable_id": "ENST4", "variation_effect_feature_cache": null},
        {"stable_id": "ENST5", "variation_effect_feature_cache": {"peptide": null}},
        {"stable_id": "ENST6", "variation_effect_feature_cache": {"peptide": 5}},
        {"stable_id": "ENST7", "variation_effect_feature_cache": {"peptide": "a", "peptide": "b"}},
        "no object", 7, null, ["a list"]
    ], "2": "no list", "3": [{"stable_id": "ENST8", "variation_effect_feature_cache": {}}]}"#;
    let ids = wanted(&[
        "ENST1", "ENST3", "ENST4", "ENST5", "ENST6", "ENST7", "ENST8",
    ]);
    assert_eq!(
        scan(shard, &ids),
        Ok(vec![
            found("ENST1", text("MKG")),
            found("ENST3", Peptide::Absent),
            found("ENST4", Peptide::Absent),
            found("ENST5", Peptide::Absent),
            found("ENST6", Peptide::Other),
            found("ENST7", text("b")),
            found("ENST8", Peptide::Absent),
        ])
    );
}

#[test]
fn a_shard_that_is_a_list_holds_its_transcripts_itself() {
    let shard = br#" [ {"stable_id": "A.1.2", "variation_effect_feature_cache": {"peptide": "\u004d\n\ud83d\ude00"}} ] "#;
    assert_eq!(
        scan(shard, &wanted(&["A"])),
        Ok(vec![found("A", text("M\n\u{1f600}"))])
    );
}

#[test]
fn a_key_written_twice_holds_its_last_value_at_its_first_place() {
    let shard = br#"{"1": [{"stable_id": "A", "stable_id": "B"}], "2": [{"stable_id": "C"}], "1": [{"stable_id": "B"}, {"stable_id": "A"}]}"#;
    assert_eq!(
        scan(shard, &wanted(&["A", "B", "C"])),
        Ok(vec![
            found("B", Peptide::Absent),
            found("A", Peptide::Absent),
            found("C", Peptide::Absent),
        ])
    );
}

#[test]
fn an_id_that_is_a_false_value_is_the_empty_id() {
    for id in [
        "null", "false", "0", "0.0", "-0", "0e7", "1e-999", "\"\"", "[]", "{}",
    ] {
        let shard = format!(
            r#"[{{"stable_id": {id}, "variation_effect_feature_cache": {{"peptide": "M"}}}}]"#
        );
        assert_eq!(
            scan(shard.as_bytes(), &wanted(&[""])),
            Ok(vec![found("", text("M"))]),
            "{id}"
        );
    }
    // A true value that is no string has no `split`.
    for id in ["1", "true", "0.5", "[0]", "{\"a\": 1}", "NaN", "-Infinity"] {
        let shard = format!(r#"[{{"stable_id": {id}}}]"#);
        assert!(scan(shard.as_bytes(), &wanted(&[])).is_err(), "{id}");
    }
}

#[test]
fn a_cache_that_is_a_true_value_and_no_object_is_unreadable_for_a_wanted_id() {
    let shard = br#"[{"stable_id": "A", "variation_effect_feature_cache": "text"},
        {"stable_id": "B", "variation_effect_feature_cache": [1]},
        {"stable_id": "C", "variation_effect_feature_cache": ""},
        {"stable_id": "D", "variation_effect_feature_cache": 3},
        {"stable_id": "A", "variation_effect_feature_cache": {"peptide": "M"}},
        {"stable_id": "A", "variation_effect_feature_cache": true}]"#;
    assert_eq!(
        scan(shard, &wanted(&["A", "C"])),
        Ok(vec![
            Event::Unreadable {
                id: "A".to_string()
            },
            found("C", Peptide::Absent),
            found("A", text("M")),
        ])
    );
}

#[test]
fn text_that_json_load_refuses_is_refused() {
    for shard in [
        "",
        " ",
        "{",
        "[1,]",
        "{\"a\": 1,}",
        "{'a': 1}",
        "{\"a\" 1}",
        "[1 2]",
        "[] []",
        "\"tab\there\"",
        "\"bad \\x escape\"",
        "\"short \\u12\"",
        "\"open",
        "01",
        "-",
        "+1",
        ".5",
        "nul",
        "\u{feff}[]",
        "7",
        "\"a string\"",
        "null",
        "true",
    ] {
        assert!(scan(shard.as_bytes(), &wanted(&[])).is_err(), "{shard:?}");
    }
    assert!(scan(b"[\xff]", &wanted(&[])).is_err());
    assert!(scan(format!("[{}]", "9".repeat(4301)).as_bytes(), &wanted(&[])).is_err());
    // What it takes.
    for shard in [
        "[]",
        "{}",
        " \t\r\n[ ] \n",
        "[1.5e+3, -0.0, 1E5, NaN, Infinity, -Infinity, 1.0e-2]",
        "[\"\\ud800 alone\", \"\\u00e9\\/\\b\\f\\\"\\\\\"]",
        "{\"a\": {\"b\": [null, true, false]}}",
    ] {
        assert_eq!(
            scan(shard.as_bytes(), &wanted(&[])),
            Ok(vec![]),
            "{shard:?}"
        );
    }
    assert!(scan(format!("[{}]", "9".repeat(4300)).as_bytes(), &wanted(&[])).is_ok());
}

#[test]
fn nesting_past_the_limit_is_refused_whatever_its_depth() {
    let lists = |depth: usize| [vec![b'['; depth], vec![b']'; depth]].concat();
    assert_eq!(scan(&lists(1000), &wanted(&[])), Ok(vec![]));
    for depth in [1002, 200_000] {
        let refused = scan(&lists(depth), &wanted(&[])).expect_err("a refusal");
        assert!(
            refused.starts_with("nesting deeper than 1,000 levels"),
            "{depth} lists: {refused}"
        );
    }
    // Objects count as lists do.
    let objects = |depth: usize| {
        let mut text = "{\"k\":".repeat(depth).into_bytes();
        text.extend_from_slice(b"null");
        text.extend(std::iter::repeat_n(b'}', depth));
        text
    };
    assert_eq!(scan(&objects(1000), &wanted(&[])), Ok(vec![]));
    assert!(scan(&objects(200_000), &wanted(&[])).is_err());
}
