// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

fn tally_of(data: &str) -> Vec<(String, u64)> {
    Tally::of(data.as_bytes())
        .contigs
        .into_iter()
        .map(|(contig, count)| (String::from_utf8(contig).unwrap(), count))
        .collect()
}

#[test]
fn contigs_are_tallied_in_the_order_first_met() {
    assert_eq!(
        tally_of("##h\n2\t5\n1\t6\n2\t7\n#c\n\nno tab\n2"),
        [
            ("2".to_string(), 3),
            ("1".to_string(), 1),
            (String::new(), 1),
            ("no tab".to_string(), 1)
        ]
    );
    assert_eq!(tally_of(""), []);
}

#[test]
fn a_tally_merged_from_runs_is_the_tally_of_the_file() {
    let data = "1\ta\n2\tb\n1\tc\n3\td\n2\te\n3\tf\n4\tg\n";
    let whole = tally_of(data);
    for cut in memchr::memchr_iter(b'\n', data.as_bytes()) {
        let mut merged = Tally::of(&data.as_bytes()[..=cut]);
        merged.merge(Tally::of(&data.as_bytes()[cut + 1..]));
        let merged: Vec<(String, u64)> = merged
            .contigs
            .into_iter()
            .map(|(contig, count)| (String::from_utf8(contig).unwrap(), count))
            .collect();
        assert_eq!(merged, whole);
    }
}

fn slice_of(records: u64, tally: &str, data: &str) -> (String, String) {
    let mut slice = Slice::new(records, tally.as_bytes());
    let mut out = Vec::new();
    // Line by line: a slice counts across the runs it is handed.
    for line in data.split_inclusive('\n') {
        slice.take(line.as_bytes(), &mut out);
    }
    (
        String::from_utf8(out).unwrap(),
        String::from_utf8(slice.table(b"label")).unwrap(),
    )
}

#[test]
fn every_contig_gives_its_share_evenly_spaced_and_no_more_than_it_holds() {
    let data = "#h\na\t1\na\t2\na\t3\na\t4\nb\t1\n#mid\nc\t1\nc\t2\nc\t3\na\t5\na\t6";
    let tally = "a\t6\nb\t1\nc\t3\n";
    // Seven over three contigs: three for the first, two for the others, b holds one.
    let (kept, table) = slice_of(7, tally, data);
    assert_eq!(kept, "#h\na\t2\na\t4\nb\t1\n#mid\nc\t2\nc\t3\na\t6\n");
    assert_eq!(table, "label\ta\t6\t3\nlabel\tb\t1\t1\nlabel\tc\t3\t2\n");
    // More than the file holds: every record.
    let (kept, table) = slice_of(100, tally, data);
    assert_eq!(kept, format!("{data}\n"));
    assert_eq!(table, "label\ta\t6\t6\nlabel\tb\t1\t1\nlabel\tc\t3\t3\n");
    // Fewer than the contigs: the first contigs of the tally give one each.
    let (kept, table) = slice_of(2, tally, data);
    assert_eq!(kept, "#h\nb\t1\n#mid\na\t6\n");
    assert_eq!(table, "label\ta\t6\t1\nlabel\tb\t1\t1\nlabel\tc\t3\t0\n");
}

#[test]
fn a_contig_the_tally_does_not_hold_gives_no_record() {
    let (kept, table) = slice_of(4, "a\t2\n", "a\t1\nz\t1\na\t2\n");
    assert_eq!(kept, "a\t1\na\t2\n");
    assert_eq!(table, "label\ta\t2\t2\nlabel\tz\t\t\n");
}
