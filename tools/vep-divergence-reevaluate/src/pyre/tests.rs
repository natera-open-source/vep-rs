// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::Regex;
use crate::pyjson::loads;
use crate::value::Value;

fn span(value: &Value) -> (usize, usize) {
    let Value::List(pair) = value else {
        panic!("a span is a pair: {value:?}");
    };
    let number = |v: &Value| v.py_str().parse::<usize>().expect("a character index");
    (number(&pair[0]), number(&pair[1]))
}

/// `tests/fixtures/regex_cases.jsonl` holds what CPython's `re` gives a pattern on a text:
/// the message of a pattern it refuses, or every match of `finditer`, the groups of
/// `search`, the span of `match`, whether `fullmatch` matches, and the texts `sub` and
/// `split` give.
#[test]
fn patterns_read_and_match_as_cpython_re() {
    let table = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/regex_cases.jsonl"
    ))
    .expect("the case table");
    let mut cases = 0;
    let mut failures = Vec::new();
    for line in table.lines() {
        let case = loads(line).expect("a JSON object a line");
        let pattern = case.get_str("pattern").expect("a pattern");
        cases += 1;
        let compiled = Regex::new(pattern);
        if let Some(error) = case.get("error") {
            let got = compiled.as_ref().err().map(|e| e.to_string());
            let expected = error.as_str().map(str::to_string);
            if got != expected {
                failures.push(format!("{pattern:?}: error {got:?}, expected {expected:?}"));
            }
            continue;
        }
        let regex = match compiled {
            Ok(regex) => regex,
            Err(error) => {
                failures.push(format!("{pattern:?}: refused with {error}"));
                continue;
            }
        };
        let text = case.get_str("text").expect("a text");
        let label = format!("{pattern:?} on {text:?}");
        let Some(Value::List(spans)) = case.get("spans") else {
            panic!("{label}: no spans");
        };
        let expected: Vec<(usize, usize)> = spans.iter().map(span).collect();
        let found = regex.find_iter(text);
        let got: Vec<(usize, usize)> = found.iter().map(|m| (m.start(), m.end())).collect();
        if got != expected {
            failures.push(format!("{label}: finditer {got:?}, expected {expected:?}"));
            continue;
        }
        let groups = regex.search(text).map(|m| {
            let mut groups = Vec::new();
            let mut index = 0;
            while index < m.spans.len() {
                groups.push(m.group(index).map_or(Value::None, Value::Str));
                index += 1;
            }
            Value::List(groups)
        });
        if groups.as_ref().unwrap_or(&Value::None) != case.get("groups").expect("groups") {
            failures.push(format!("{label}: groups {groups:?}"));
        }
        if regex.is_match(text) == expected.is_empty() {
            failures.push(format!("{label}: is_match"));
        }
        let matched = regex.match_start(text).map(|m| (m.start(), m.end()));
        let expected_match = match case.get("match").expect("match") {
            Value::None => None,
            pair => Some(span(pair)),
        };
        if matched != expected_match {
            failures.push(format!(
                "{label}: match {matched:?}, expected {expected_match:?}"
            ));
        }
        if Value::Bool(regex.is_full_match(text)) != *case.get("fullmatch").expect("fullmatch") {
            failures.push(format!("{label}: fullmatch"));
        }
        let sub = regex.sub("<masked>", text);
        if Some(sub.as_str()) != case.get_str("sub") {
            failures.push(format!("{label}: sub {sub:?}"));
        }
        if let Some(split) = case.get("split") {
            let got = Value::List(regex.split(text).into_iter().map(Value::Str).collect());
            if got != *split {
                failures.push(format!("{label}: split {got:?}"));
            }
        }
    }
    assert!(cases > 1000, "the table holds {cases} cases");
    assert!(
        failures.is_empty(),
        "{} of {cases} cases differ:\n{}",
        failures.len(),
        failures
            .iter()
            .take(60)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn a_named_group_is_read_by_its_name() {
    let regex = Regex::new(r"^(?P<tier>[1-3]) · (?P<figure>[^·]+?) · (?P<tail>.+)$").unwrap();
    let found = regex
        .match_start("2 · 156.8 per million · 0 flag sets, 28 corpora")
        .unwrap();
    assert_eq!(found.name("tier").as_deref(), Some("2"));
    assert_eq!(found.name("figure").as_deref(), Some("156.8 per million"));
    assert_eq!(
        found.name("tail").as_deref(),
        Some("0 flag sets, 28 corpora")
    );
    assert_eq!(found.name("absent"), None);
    assert!(regex.match_start("tier 2").is_none());
}
