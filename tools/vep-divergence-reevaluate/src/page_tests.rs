// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::dump::{dump_page, section_value};
use crate::pyjson::{dumps_compact_sorted, loads};
use crate::value::Value;

const FIXTURES: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures");

/// `tests/fixtures/expected_page_dump.jsonl` holds the sections, rows and named keys
/// `check_divergence_order.py` reads of `tests/fixtures/page.md`.
#[test]
fn the_fixture_page_is_read_as_the_order_checker_reads_it() {
    let text = std::fs::read_to_string(format!("{FIXTURES}/page.md")).expect("the page");
    let sections = parse_page(&text, &PagePatterns::new()).expect("a page that parses");
    let expected = std::fs::read_to_string(format!("{FIXTURES}/expected_page_dump.jsonl"))
        .expect("the expected dump");
    assert_eq!(dump_page(&sections), expected);
    let keys: Vec<&str> = leaves(&sections)
        .iter()
        .flat_map(|section| section.rows.iter().map(|row| row.key.as_str()))
        .collect();
    assert_eq!(
        keys,
        [
            "first_entry",
            "second_entry",
            "third_entry",
            "fourth_entry",
            "plugin_entry",
            "cadd_entry",
            "broken_entry",
            "broken_cadd_entry"
        ]
    );
}

/// `tests/fixtures/page_cases.jsonl` holds small pages, each with the sections the order
/// checker reads of it or the line it refuses the page with.
#[test]
fn small_pages_are_read_or_refused_as_the_order_checker_does() {
    let table =
        std::fs::read_to_string(format!("{FIXTURES}/page_cases.jsonl")).expect("the case table");
    let patterns = PagePatterns::new();
    let mut cases = 0;
    for line in table.lines() {
        let case = loads(line).expect("a JSON object a line");
        let text = case.get_str("text").expect("a page text");
        let parsed = parse_page(text, &patterns);
        match case.get_str("error") {
            Some(error) => assert_eq!(parsed, Err(error.to_string()), "{text:?}"),
            None => {
                let sections = parsed.unwrap_or_else(|error| panic!("{text:?}: {error}"));
                let got = Value::List(sections.iter().map(section_value).collect());
                assert_eq!(
                    dumps_compact_sorted(&got),
                    dumps_compact_sorted(case.get("sections").expect("sections")),
                    "{text:?}"
                );
            }
        }
        cases += 1;
    }
    assert!(cases > 20, "the table holds {cases} cases");
}

#[test]
fn class_counts_name_the_three_classes_in_order() {
    let classes = [
        "observed",
        "malformed",
        "observed",
        "unobserved",
        "observed",
    ];
    assert_eq!(
        class_counts(classes.iter().copied()),
        "3 observed, 1 unobserved, 1 malformed"
    );
    assert_eq!(
        class_counts(std::iter::empty()),
        "0 observed, 0 unobserved, 0 malformed"
    );
}
