// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A Location is keyed by the numbers of its coordinates.

use super::*;
use crate::TempDir;

fn canonical(location: &str) -> String {
    let mut written = Vec::new();
    push_canonical_location(location, &mut written);
    String::from_utf8(written).unwrap()
}

#[test]
fn a_location_is_keyed_by_the_numbers_of_its_coordinates() {
    for (written, keyed) in [
        ("21:30425651-030425949", "21:30425651-30425949"),
        ("21:30425651-+30425949", "21:30425651-30425949"),
        ("21:+30425650", "21:30425650"),
        ("21:30425651-30425949.5", "21:30425651-30425949.5"),
        ("21:30425651", "21:30425651"),
        ("MT", "MT"),
    ] {
        assert_eq!(canonical(written), keyed, "{written}");
    }
    let scratch = TempDir::create(&std::env::temp_dir(), "location_test_").unwrap();
    let output = scratch.path.join("perl.txt");
    std::fs::write(
        &output,
        "x\t21:30425651-030425949\tdeletion\tG\tT1\tTranscript\tupstream_gene_variant\n",
    )
    .unwrap();
    let mut prep = KeyPrep::default();
    let mut keys = Vec::new();
    for_each_record(File::open(&output).unwrap(), |record| {
        let mut key = Vec::new();
        if prep.prepare(record, false, None, &mut key) {
            keys.push(String::from_utf8(key).unwrap());
        }
        Ok(())
    })
    .unwrap();
    assert_eq!(
        keys,
        ["21:30425651-30425949\tdeletion\tT1\tTranscript\tupstream_gene_variant"],
        "an end written with a leading zero keys the row at its number"
    );
}
