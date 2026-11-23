// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The matcher against CPython 3.14's `difflib`.
//!
//! `tests/fixtures/difflib/cases.txt` holds pairs of sequences with the places
//! `difflib.unified_diff(a, b, n=0)` marks: four lines a case, `a` and the items of the
//! first sequence, `b` and those of the second, `-` and the places of `a` that are
//! removed, `+` and the places of `b` that are added. The sequences run from none to
//! 320 items over 2 to 600 distinct ones, on both sides of the 200 items from which
//! an item that fills more than one place in a hundred takes no part in the search.

use super::changed_lines;

fn items<'a>(line: &'a str, mark: &str) -> Vec<&'a str> {
    let rest = line
        .strip_prefix(mark)
        .unwrap_or_else(|| panic!("a line that opens with {mark:?}: {line:?}"));
    rest.split_whitespace().collect()
}

#[test]
fn the_marked_lines_are_difflibs() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/difflib/cases.txt"
    );
    let text = std::fs::read_to_string(path).expect("the cases");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len() % 4, 0);
    let places = |line: &str, mark: &str| -> Vec<usize> {
        items(line, mark)
            .iter()
            .map(|place| place.parse().expect("a place"))
            .collect()
    };
    let mut cases = 0;
    for case in lines.chunks_exact(4) {
        let (a, b) = (items(case[0], "a"), items(case[1], "b"));
        let (removed, added) = changed_lines(&a, &b);
        assert_eq!(removed, places(case[2], "-"), "removed, case {cases}");
        assert_eq!(added, places(case[3], "+"), "added, case {cases}");
        cases += 1;
    }
    assert_eq!(cases, 154);
}

#[test]
fn equal_sequences_have_no_marked_line() {
    let lines = ["{", "  'a' => 1", "}"];
    assert_eq!(changed_lines(&lines, &lines), (vec![], vec![]));
    assert_eq!(changed_lines(&[], &[]), (vec![], vec![]));
}
