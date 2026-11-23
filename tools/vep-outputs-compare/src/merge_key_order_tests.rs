// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The merge of two key files walks each forward once and takes each to ascend with every
//! key once. Past a key below the key before it the walk would step over a key both files
//! hold and count it as missing and as extra at once, so it stops there, carries the key and
//! the key before it, and names the side and neither key.

use super::*;
use crate::TempDir;

/// The merge of two key files that hold `reference` and `scored` in the order given, with no
/// discordant row written and no consequence set counted.
fn merge_of(reference: &[&str], scored: &[&str]) -> io::Result<MergeCounts> {
    let scratch = TempDir::create(&std::env::temp_dir(), "key_order_test_")?;
    let (perl, rust) = (
        scratch.path.join("perl.keys"),
        scratch.path.join("rust.keys"),
    );
    let lines = |keys: &[&str]| {
        keys.iter()
            .map(|key| format!("{key}\n"))
            .collect::<String>()
    };
    std::fs::write(&perl, lines(reference))?;
    std::fs::write(&rust, lines(scored))?;
    compare_sorted_key_files(&perl, &rust, "f", None, false, Collector::None)
}

/// The stop of the merge: its message, the key it is on and the key before it.
fn stop_of(reference: &[&str], scored: &[&str]) -> (String, String, String) {
    match merge_of(reference, scored) {
        Err(error) => {
            let order = key_order(&error).unwrap_or_else(|| panic!("{error}"));
            let text = |key: &[u8]| String::from_utf8(key.to_vec()).unwrap();
            (error.to_string(), text(order.key()), text(order.before()))
        }
        Ok(counts) => panic!(
            "the merge went on and counted {} keys on both sides",
            counts.intersection
        ),
    }
}

/// The whole message of a stop on the side `option` names, before any line of an output
/// file is named for its two keys.
fn message_of_a_stop_on(option: &str) -> String {
    format!(
        "the sorted keys of the {option} side do not ascend: a key of the side does not sort \
         after the key before it; the comparison needs each side's keys in order with every \
         key once"
    )
}

#[test]
fn a_key_below_the_key_before_it_on_the_scored_side_names_that_side_and_no_key() {
    let (message, key, before) = stop_of(
        &["a1", "b2", "c3", "d4", "e5"],
        &["a1", "b2", "c3", "e5", "d4"],
    );
    assert_eq!(message, message_of_a_stop_on("--rust-dir"));
    assert_eq!((key.as_str(), before.as_str()), ("d4", "e5"));
}

#[test]
fn a_key_out_of_order_after_a_key_the_reference_alone_holds_names_that_side_and_no_key() {
    let (message, key, before) = stop_of(&["a1", "b2", "d4", "c3"], &["a1", "x9"]);
    assert_eq!(message, message_of_a_stop_on("--perl-dir"));
    assert_eq!((key.as_str(), before.as_str()), ("c3", "d4"));
}

#[test]
fn a_key_out_of_order_after_a_key_the_scored_side_alone_holds_names_that_side_and_no_key() {
    let (message, key, before) = stop_of(&["a1", "x9"], &["a1", "b2", "d4", "c3"]);
    assert_eq!(message, message_of_a_stop_on("--rust-dir"));
    assert_eq!((key.as_str(), before.as_str()), ("c3", "d4"));
}

#[test]
fn keys_that_ascend_are_counted_to_the_last_key_of_each_side() {
    let counts = merge_of(&["a1", "b2", "d4", "e5", "f6"], &["a1", "c3", "d4"]).unwrap();
    assert_eq!(
        (
            counts.perl_count,
            counts.rust_count,
            counts.intersection,
            counts.missing,
            counts.extra
        ),
        (5, 3, 2, 3, 1)
    );
}
