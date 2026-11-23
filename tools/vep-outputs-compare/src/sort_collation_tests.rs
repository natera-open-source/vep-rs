// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The sorts order keys by byte. The merge compares keys by byte, so a sort in any other
//! order, a locale's collation among them, would walk it past keys both sides hold and lose
//! them from the intersection with no error raised.

use super::*;
use crate::TempDir;

/// Two keys whose byte order and numeric order differ: `0` sorts before `:`.
const KEYS: [&str; 2] = ["1:117329", "10:110755"];

/// The lines a [`KeySorter`] writes for `keys`, spilling a run at `memory_bound` bytes.
fn sorted_to_a_file(keys: &[&str], memory_bound: usize) -> Vec<String> {
    let scratch = TempDir::create(&std::env::temp_dir(), "collation_test_").unwrap();
    let mut sorter = KeySorter::new(&scratch.path, "keys", memory_bound);
    for key in keys {
        sorter.push(key.as_bytes()).unwrap();
    }
    let out = scratch.path.join("keys.sorted");
    sorter.finish(&out).unwrap();
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

/// The keys a partition's sort holds for `keys`.
fn sorted_in_memory(keys: &[&str]) -> Vec<String> {
    let lines: String = keys.iter().map(|key| format!("{key}\n")).collect();
    let sorted = SortedKeys::from_lines(lines.into_bytes());
    (0..sorted.len())
        .map(|index| String::from_utf8(sorted.get(index).to_vec()).unwrap())
        .collect()
}

#[test]
fn keys_are_written_in_byte_order_where_it_differs_from_numeric_order() {
    assert_eq!(
        sorted_to_a_file(&KEYS, DEFAULT_MEMORY_BOUND),
        ["10:110755", "1:117329"]
    );
}

/// Any of the three may sort a side's keys: the sort held in memory, the merge of spilled
/// runs past the memory bound, and the sort of one partition.
#[test]
fn every_sort_writes_the_order_of_a_bytewise_sort() {
    let mut bytewise = KEYS.to_vec();
    bytewise.sort_unstable_by(|a, b| a.as_bytes().cmp(b.as_bytes()));
    assert_eq!(sorted_to_a_file(&KEYS, DEFAULT_MEMORY_BOUND), bytewise);
    assert_eq!(sorted_to_a_file(&KEYS, 1), bytewise);
    assert_eq!(sorted_in_memory(&KEYS), bytewise);
}
