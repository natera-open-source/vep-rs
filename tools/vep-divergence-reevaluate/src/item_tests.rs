// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;
use crate::pyjson::loads;
use crate::value::Value;

/// `tests/fixtures/item_span_cases.jsonl` holds Rust source texts, each with an item's
/// name and the first and the last line of its span, `null` where the text defines no
/// such item. `scripts/ci/test_admission.py` holds the admission script's reader to the
/// same table, so the two readers give one answer.
#[test]
fn an_item_spans_the_lines_the_case_table_states() {
    let table = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/item_span_cases.jsonl"
    ))
    .expect("the case table");
    let mut cases = 0;
    let mut undefined = 0;
    for line in table.lines() {
        let case = loads(line).expect("a JSON object a line");
        let text = case.get_str("text").expect("a source text");
        let item = case.get_str("item").expect("an item's name");
        let line_at = |key: &str| {
            case.get(key)
                .and_then(Value::as_i128)
                .map(|number| usize::try_from(number).expect("a line number"))
        };
        let expected = line_at("first").zip(line_at("last"));
        assert_eq!(item_span(text, item), expected, "{item} in {text:?}");
        cases += 1;
        undefined += usize::from(expected.is_none());
    }
    assert!(cases > 100, "the table holds {cases} cases");
    assert!(
        undefined > 10,
        "{undefined} cases name an item no text defines"
    );
}

#[test]
fn items_are_named_by_their_chain_and_listed_as_they_end() {
    let text =
        "struct Thing;\nimpl Thing {\n    fn method() {}\n}\nmod tests {\n    fn pins() {}\n}\n";
    let named: Vec<(&str, usize, usize)> = vec![
        ("Thing", 1, 1),
        ("Thing::method", 3, 3),
        ("Thing", 2, 4),
        ("tests::pins", 6, 6),
        ("tests", 5, 7),
    ];
    let found = items(text);
    let found: Vec<(&str, usize, usize)> = found
        .iter()
        .map(|(name, first, last)| (name.as_str(), *first, *last))
        .collect();
    assert_eq!(found, named);
    // The first item of a name is the one a row means.
    assert_eq!(item_span(text, "Thing"), Some((1, 1)));
    assert_eq!(item_span(text, "pins"), None);
}

/// The fixture source the fixture records name defines the items their tables name.
#[test]
fn the_fixture_source_defines_the_items_the_fixture_records_name() {
    let source = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/implementation/recorded.rs"
    ))
    .expect("the fixture source");
    for item in [
        "downstream_row",
        "PLAIN_EQUALS_SIGN",
        "is_escaped_equals_sign",
        "tests::the_downstream_row_is_written",
        "tests::the_equals_sign_is_written_plain",
    ] {
        assert!(item_span(&source, item).is_some(), "{item}");
    }
}
