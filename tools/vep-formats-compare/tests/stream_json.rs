// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reads planted pairs of JSON Lines outputs and their input, every line written here,
//! through the streamed record groups at one worker, two and the host's default and
//! through blocks of several sizes, and holds the readings to one remainder and to the
//! remainder each shape must leave.

mod stream_records_common;

use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use std::path::Path;
use stream_records_common::{
    block as counts, counter, line_runs, pairs, plan_of, remainder_at, texts_of, Remainder, Scratch,
};
use vep_formats_compare::stream::{self, Format};

/// Two items of a consequence list: an intron row and a row with an exon number.
const TX1: &str = r#"{"variant_allele":"G","transcript_id":"tx1","consequence_terms":["intron_variant"],"impact":"MODIFIER"}"#;
const TX2: &str = r#"{"variant_allele":"G","transcript_id":"tx2","consequence_terms":["missense_variant"],"impact":"MODERATE","exon":"2/5"}"#;

/// Record `record`'s line of the input.
fn input_line(record: usize) -> String {
    format!("1\t{}\trs{record}\tA\tG\t.\t.\t.", 1000 + record)
}

/// The members of record `record`'s object before its consequence lists.
fn head(record: usize) -> String {
    format!(
        r#""input":"{}","id":"rs{record}","seq_region_name":"1","start":{},"end":{},"allele_string":"A/G","strand":1"#,
        input_line(record).replace('\t', "\\t"),
        1000 + record,
        1000 + record
    )
}

/// Record `record`'s object, its transcript consequences the `items`, with `more` after
/// them: nothing, or further members each after a comma.
fn object_with(record: usize, items: &[&str], more: &str) -> String {
    format!(
        r#"{{{},"transcript_consequences":[{}]{more}}}"#,
        head(record),
        items.join(",")
    )
}

fn object(record: usize, items: &[&str]) -> String {
    object_with(record, items, "")
}

fn output(rows: &[String]) -> String {
    rows.iter().map(|row| format!("{row}\n")).collect()
}

/// Ten records, two transcript consequences each.
fn ten_records() -> Vec<String> {
    (0..10).map(|r| object(r, &[TX1, TX2])).collect()
}

fn ten_inputs() -> Vec<String> {
    (0..10).map(input_line).collect()
}

fn remainder_of(expected: &Path, actual: &Path, input: &Path, more: &[&str]) -> Remainder {
    stream_records_common::remainder_of(Format::Json, expected, actual, input, more)
}

fn remainder_of_texts(expected: &str, actual: &str, inputs: &[String], more: &[&str]) -> Remainder {
    let scratch = Scratch::new();
    let e = scratch.file("expected.json", expected.as_bytes());
    let a = scratch.file("actual.json", actual.as_bytes());
    let i = scratch.file("input.vcf", output(inputs).as_bytes());
    remainder_of(&e, &a, &i, more)
}

fn remainder(
    expected: &[String],
    actual: &[String],
    inputs: &[String],
    more: &[&str],
) -> Remainder {
    remainder_of_texts(&output(expected), &output(actual), inputs, more)
}

/// The lines of a remainder text.
fn data(text: &str) -> Vec<String> {
    text.lines().map(str::to_string).collect()
}

#[test]
fn outputs_written_alike_leave_no_line() {
    let rows = ten_records();
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert!(r.expected.is_empty() && r.actual.is_empty());
    assert_eq!(counter(&r.counters, "whole"), "null");
    // A line is its record entry and one entry a consequence item.
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
    let every = pairs(&[
        ("allele_string", 10),
        ("consequence_terms", 20),
        ("end", 10),
        ("exon", 10),
        ("id", 10),
        ("impact", 20),
        ("seq_region_name", 10),
        ("start", 10),
        ("strand", 10),
        ("transcript_id", 20),
        ("variant_allele", 20),
    ]);
    assert_eq!(counts(&r.counters, "kept_populated"), every);
    assert_eq!(counts(&r.counters, "kept_present"), every);
    assert_eq!(counter(&r.counters, "kept_under"), "{}");
    assert_eq!(counter(&r.counters, "kept_before"), "[]");
    assert_eq!(counter(&r.counters, "line_runs"), line_runs("[]", "[]"));
    assert_eq!(counter(&r.counters, "groups"), "10");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "0");
    assert_eq!(
        counter(&r.counters, "expected"),
        "{\"lines\": 10, \"lines_handed_over\": 0, \"bytes_handed_over\": 0}"
    );
}

#[test]
fn lines_of_the_same_entries_in_another_text_are_kept_back() {
    let expected = ten_records();
    // The members of an object in another order, white space between them, an escape for
    // a solidus, and the identity of the engine that wrote the line.
    let mut actual = expected.clone();
    actual[2] = format!(
        r#"{{ "transcript_consequences" : [{TX1}, {}], "vep_rs_version":"1.2.3", {} }}"#,
        r#"{"exon":"2\/5","impact":"MODERATE","consequence_terms":["missense_variant"],"transcript_id":"tx2","variant_allele":"G"}"#,
        head(2)
    );
    let r = remainder(&expected, &actual, &ten_inputs(), &[]);
    assert!(r.expected.is_empty() && r.actual.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
    assert_eq!(counter(&r.counters, "transcript_id"), "20");
    // The items of a list in another order are other entries in order, and a whole number
    // spelt with a fraction is another value.
    let mut actual = expected.clone();
    actual[2] = object(2, &[TX2, TX1]);
    actual[5] = actual[5].replace("\"strand\":1", "\"strand\":1.0");
    let r = remainder(&expected, &actual, &ten_inputs(), &[]);
    assert_eq!(
        data(&r.expected),
        [&expected[2..3], &expected[5..6]].concat()
    );
    assert_eq!(data(&r.actual), [&actual[2..3], &actual[5..6]].concat());
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 6], [3, 12]]");
}

#[test]
fn a_differing_member_hands_its_line_over() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual[3] = actual[3].replace("MODERATE", "HIGH");
    let r = remainder(&expected, &actual, &ten_inputs(), &[]);
    assert_eq!(data(&r.expected), expected[3..4]);
    assert_eq!(data(&r.actual), actual[3..4]);
    assert_eq!(counter(&r.counters, "kept_entries"), "27");
    // The line's three entries stand after the nine entries of the three lines before it.
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 9]]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[1, 3]]", "[[1, 3]]")
    );
    assert_eq!(counter(&r.counters, "impact"), "18");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "1");
    assert_eq!(
        counter(&r.counters, "expected"),
        format!(
            "{{\"lines\": 10, \"lines_handed_over\": 1, \"bytes_handed_over\": {}}}",
            expected[3].len() + 1
        )
    );
    // A member one side alone writes, at the top level and in an item.
    let mut actual = expected.clone();
    actual[0] = object_with(
        0,
        &[TX1, TX2],
        r#","most_severe_consequence":"missense_variant""#,
    );
    actual[9] = actual[9].replace(r#","exon":"2/5""#, "");
    let r = remainder(&expected, &actual, &ten_inputs(), &[]);
    assert_eq!(data(&r.expected), [&expected[..1], &expected[9..]].concat());
    assert_eq!(data(&r.actual), [&actual[..1], &actual[9..]].concat());
}

#[test]
fn a_blank_line_is_no_line_of_a_group_and_a_line_of_the_output() {
    let rows = ten_records();
    // Blank lines at other places on the two sides: empty, and white space alone.
    let mut expected = format!("\n{}", output(&rows[..4]));
    expected.push_str(" \t\n");
    expected.push_str(&output(&rows[4..]));
    let mut actual = output(&rows[..7]);
    actual.push_str("\n\n");
    actual.push_str(&output(&rows[7..]));
    let r = remainder_of_texts(&expected, &actual, &ten_inputs(), &[]);
    assert!(r.expected.is_empty() && r.actual.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
    assert_eq!(
        counter(&r.counters, "expected"),
        "{\"lines\": 10, \"lines_handed_over\": 0, \"bytes_handed_over\": 0}"
    );
    // A line handed over stands in the output after the blank lines before it too.
    let other = actual.replace("\"rs8\"", "\"rs80\"");
    let r = remainder_of_texts(&expected, &other, &ten_inputs(), &[]);
    assert_eq!(data(&r.expected), rows[8..9]);
    assert_eq!(
        r.actual,
        format!("{}\n", rows[8].replace("\"rs8\"", "\"rs80\""))
    );
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 24]]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[1, 10]]", "[[1, 10]]")
    );
    // With everything handed over the blank lines are left out all the same.
    let mut broken = expected.clone();
    broken.push_str("{\"input\":\n");
    let r = remainder_of_texts(&broken, &actual, &ten_inputs(), &[]);
    assert_eq!(
        counter(&r.counters, "whole"),
        "\"a line that is not one JSON object\""
    );
    assert_eq!(r.expected, format!("{}{{\"input\":\n", output(&rows)));
    assert_eq!(r.actual, output(&rows));
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[4, 1], [7, 2]]", "[[7, 0], [3, 2]]")
    );
    assert_eq!(counter(&r.counters, "kept_before"), "[[31, 0]]");
}

#[test]
fn a_key_of_more_lines_than_the_input_holds_is_handed_over() {
    // Record 2's line twice on both sides, the same text.
    let mut rows = ten_records();
    rows.insert(3, rows[2].clone());
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert_eq!(data(&r.expected), rows[2..4]);
    assert_eq!(data(&r.actual), rows[2..4]);
    assert_eq!(counter(&r.counters, "kept_entries"), "27");
    assert_eq!(counter(&r.counters, "kept_before"), "[[6, 6]]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[2, 2]]", "[[2, 2]]")
    );
    // The input holding the line twice, each output line names one of the two.
    let mut inputs = ten_inputs();
    inputs.push(input_line(2));
    let r = remainder(&rows, &rows, &inputs, &[]);
    assert!(r.expected.is_empty() && r.actual.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "33");
    assert_eq!(counter(&r.counters, "groups"), "10");
}

#[test]
fn a_line_of_a_key_the_input_lacks_is_handed_over() {
    let rows = ten_records();
    let mut inputs = ten_inputs();
    inputs.remove(7);
    let r = remainder(&rows, &rows, &inputs, &[]);
    assert_eq!(data(&r.expected), rows[7..8]);
    assert_eq!(data(&r.actual), rows[7..8]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 21]]");
    // The key is the first five columns of the `input` text: its later columns are no
    // part of it, and an object without the text has the empty key.
    let mut later = rows.clone();
    later[7] = later[7].replace("G\\t.\\t.\\t.\"", "G\\t50\\tq10\\tAC=2\"");
    let r = remainder(&later, &later, &ten_inputs(), &[]);
    assert!(r.expected.is_empty());
    let mut bare = rows.clone();
    bare[7] = r#"{"id":"rs7","transcript_consequences":[]}"#.to_string();
    let r = remainder(&bare, &bare, &ten_inputs(), &[]);
    assert_eq!(data(&r.expected), bare[7..8]);
    let mut inputs = ten_inputs();
    inputs[7] = String::new();
    let r = remainder(&bare, &bare, &inputs, &[]);
    assert!(r.expected.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "28");
}

#[test]
fn a_listed_record_is_handed_over_with_equal_lines() {
    let rows = ten_records();
    let scratch = Scratch::new();
    let list = scratch.file("records.txt", b"4\n99\n");
    let more = ["--force-records", list.to_str().unwrap()];
    let r = remainder(&rows, &rows, &ten_inputs(), &more);
    assert_eq!(data(&r.expected), rows[4..5]);
    assert_eq!(data(&r.actual), rows[4..5]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 12]]");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "1");
}

#[test]
fn a_reference_list_of_a_kind_off_the_kept_list_hands_its_line_over() {
    let regulatory = r#","regulatory_feature_consequences":[{"variant_allele":"G","regulatory_feature_id":"reg1","consequence_terms":["regulatory_region_variant"]}]"#;
    let mut rows = ten_records();
    rows[5] = object_with(5, &[TX1], regulatory);
    // A list without an item is no entry of its kind.
    rows[6] = object_with(6, &[TX1], r#","motif_feature_consequences":[]"#);
    let scratch = Scratch::new();
    let list = scratch.file("kept.txt", b"transcript\nintergenic\n");
    let more = ["--keep-values", "kind", list.to_str().unwrap()];
    let r = remainder(&rows, &rows, &ten_inputs(), &more);
    assert_eq!(data(&r.expected), rows[5..6]);
    assert_eq!(data(&r.actual), rows[5..6]);
    assert_eq!(counter(&r.counters, "kept_entries"), "26");
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert!(r.expected.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "29");
    // A list that leaves the transcript kind out hands over every line with an item of it.
    let list = scratch.file("regulatory.txt", b"regulatory_feature\nmotif_feature\n");
    let more = ["--keep-values", "kind", list.to_str().unwrap()];
    let r = remainder(&rows, &rows, &ten_inputs(), &more);
    assert_eq!(counter(&r.counters, "groups_handed_over"), "10");
    let bare = vec![object(0, &[])];
    let r = remainder(&bare, &bare, &ten_inputs(), &more);
    assert_eq!(counter(&r.counters, "groups_handed_over"), "0");
}

#[test]
fn a_reference_line_holding_a_listed_text_hands_its_line_over() {
    let mut rows = ten_records();
    rows[6] = rows[6].replace("\"2/5\"", "\"misspeltKey7\"");
    let scratch = Scratch::new();
    let list = scratch.file("texts.txt", b"misspeltKey\n");
    let more = ["--force-text", list.to_str().unwrap()];
    let r = remainder(&rows, &rows, &ten_inputs(), &more);
    assert_eq!(data(&r.expected), rows[6..7]);
    assert_eq!(data(&r.actual), rows[6..7]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "27");
    // In the key of a reference line, everything is handed over.
    let mut named = rows.clone();
    named[1] = named[1].replace("rs1", "misspeltKey1");
    let mut inputs = ten_inputs();
    inputs[1] = inputs[1].replace("rs1", "misspeltKey1");
    let r = remainder(&named, &named, &inputs, &more);
    assert_eq!(
        counter(&r.counters, "whole"),
        "\"a listed text stands in a header line of the reference or in a line's key\""
    );
    assert_eq!(data(&r.expected), named);
    assert_eq!(counter(&r.counters, "kept_entries"), "0");
    assert_eq!(counter(&r.counters, "kept_before"), "[[30, 0]]");
    let r = remainder(&named, &named, &inputs, &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
}

#[test]
fn lines_in_another_order_are_handed_over_unless_the_order_is_unread() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual.swap(2, 3);
    let r = remainder(&expected, &actual, &ten_inputs(), &[]);
    assert_eq!(data(&r.expected), expected[2..4]);
    assert_eq!(data(&r.actual), actual[2..4]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[6, 6]]");
    let r = remainder(&expected, &actual, &ten_inputs(), &["--unordered"]);
    assert!(r.expected.is_empty() && r.actual.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
}

#[test]
fn lines_labelled_by_their_place_are_one_group_in_the_order_of_the_shared_keys() {
    // Records 11 to 13 are no records of the input: the comparison labels their lines by
    // their place in the output, and two such lines at one place share their keys.
    let inputs = vec![input_line(0)];
    let known = object(0, &[TX1]);
    let (first, second, third) = (object(11, &[TX1]), object(12, &[TX1]), object(13, &[TX1]));
    let expected = vec![known.clone(), first.clone()];
    let actual = vec![second.clone(), third.clone(), known.clone()];
    let r = remainder(&expected, &actual, &inputs, &[]);
    assert_eq!(data(&r.expected), expected);
    assert_eq!(data(&r.actual), actual);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "0");
    let r = remainder(&expected, &actual, &inputs, &["--unordered"]);
    assert_eq!(data(&r.expected), expected[1..]);
    assert_eq!(data(&r.actual), actual[..2]);
    // On both sides after the same kept-back lines, they move nothing.
    let actual = vec![known.clone(), second.clone()];
    let r = remainder(&expected, &actual, &inputs, &[]);
    assert_eq!(data(&r.expected), expected[1..]);
    assert_eq!(data(&r.actual), actual[1..]);
    assert_eq!(counter(&r.counters, "kept_entries"), "2");
}

#[test]
fn a_line_the_comparison_may_fail_to_convert_hands_its_line_over() {
    let plain = ten_records();
    for exon in [
        // A digit outside ASCII, a run of 641 digits, 62 arrays inside one another below
        // the three levels of an item's member, and an escape that may spell any of them.
        "\"\u{b2}\"".to_string(),
        "7".repeat(641),
        format!("{}1{}", "[".repeat(62), "]".repeat(62)),
        r#""2\u002f5""#.to_string(),
    ] {
        let mut rows = plain.clone();
        rows[4] = rows[4].replace("\"2/5\"", &exon);
        let r = remainder(&rows, &rows, &ten_inputs(), &[]);
        assert_eq!(data(&r.expected), rows[4..5]);
        assert_eq!(data(&r.actual), rows[4..5]);
        assert_eq!(counter(&r.counters, "kept_entries"), "27");
    }
    // At the two bounds the line is kept back: 640 digits, 64 levels deep. A reverse
    // solidus before the letter is no escape of that kind.
    let mut rows = plain.clone();
    let exon = format!("{}{}{}", "[".repeat(61), "7".repeat(640), "]".repeat(61));
    rows[4] = rows[4].replace("\"2/5\"", &exon);
    rows[5] = rows[5].replace("\"2/5\"", r#""2\\u002f5""#);
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert!(r.expected.is_empty() && r.actual.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
}

/// A value of `levels` arrays and objects inside one another, an array the outermost.
fn nested(levels: usize) -> String {
    let mut text = "1".to_string();
    for level in (0..levels).rev() {
        text = match level % 2 {
            0 => format!("[{text}]"),
            _ => format!("{{\"a\":{text}}}"),
        };
    }
    text
}

#[test]
fn a_line_is_handed_over_by_the_depth_of_its_nesting_and_not_by_its_brackets() {
    // Forty consequence items stand side by side: 82 opening brackets, four levels.
    let items: Vec<String> = (0..40)
        .map(|t| {
            format!(
                r#"{{"variant_allele":"G","transcript_id":"tx{t}","consequence_terms":["intron_variant"]}}"#
            )
        })
        .collect();
    let items: Vec<&str> = items.iter().map(String::as_str).collect();
    let mut rows = ten_records();
    rows[4] = object(4, &items);
    assert_eq!(rows[4].matches(['[', '{']).count(), 82);
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert!(r.expected.is_empty() && r.actual.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "68");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "0");
    // A member of an item stands three levels deep: 61 levels below it are 64, kept back,
    // and 62 are 65, handed over.
    let in_an_item = |levels: usize| {
        let item = format!(
            r#"{{"variant_allele":"G","transcript_id":"tx3","deep":{}}}"#,
            nested(levels)
        );
        let mut rows = ten_records();
        rows[6] = object(6, &[TX1, &item]);
        (remainder(&rows, &rows, &ten_inputs(), &[]), rows)
    };
    let (r, _) = in_an_item(61);
    assert!(r.expected.is_empty() && r.actual.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
    let (r, rows) = in_an_item(62);
    assert_eq!(data(&r.expected), rows[6..7]);
    assert_eq!(data(&r.actual), rows[6..7]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "27");
    // A member of the object stands one level deep, and an array without an item is a
    // level: 63 arrays inside one another are 64 levels, 64 are 65.
    let at_the_top = |arrays: usize| {
        let value = format!("{}{}", "[".repeat(arrays), "]".repeat(arrays));
        let mut rows = ten_records();
        rows[6] = object_with(6, &[TX1, TX2], &format!(r#","deep":{value}"#));
        remainder(&rows, &rows, &ten_inputs(), &[])
    };
    let r = at_the_top(63);
    assert!(r.expected.is_empty() && r.actual.is_empty());
    let r = at_the_top(64);
    assert_eq!(data(&r.expected).len(), 1);
    assert_eq!(counter(&r.counters, "groups_handed_over"), "1");
    // A text of brackets is no nesting.
    let mut rows = ten_records();
    rows[6] = rows[6].replace("\"2/5\"", &format!("\"{}\"", "[{".repeat(40)));
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert!(r.expected.is_empty() && r.actual.is_empty());
}

#[test]
fn every_member_of_a_record_and_of_an_item_is_counted_by_its_presence_and_its_value() {
    // An empty text, an absent value, an empty list, a list of one absent value, a list of
    // two, an object and a null; the engine's identity and the input are no fields.
    let rows = vec![
        r#"{"input":"1\t1000\trs0\tA\tG","vep_rs_version":"1.2.3","id":"","start":1000,"colocated_variants":[{"id":"rs0"}],"flag":null,"transcript_consequences":[{"variant_allele":"-","transcript_id":"tx1","flags":[],"domains":["-"],"consequence_terms":["a","b"]}],"intergenic_consequences":[{"variant_allele":"G"}],"regulatory_feature_consequences":"none"}"#.to_string(),
        object(1, &[TX1, TX2]),
        object(2, &[TX1]),
    ];
    let inputs: Vec<String> = (0..3).map(input_line).collect();
    let r = remainder(&rows, &rows, &inputs, &[]);
    assert!(r.expected.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "8");
    let mut actual = rows.clone();
    actual[1] = object(1, &[TX1]);
    let r = remainder(&rows, &actual, &inputs, &[]);
    assert_eq!(data(&r.expected), rows[1..2]);
    assert_eq!(counter(&r.counters, "kept_entries"), "5");
    // Three entries after the three of the first line.
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 3]]");
    assert_eq!(
        counts(&r.counters, "kept_present"),
        pairs(&[
            ("allele_string", 1),
            ("colocated_variants", 1),
            ("consequence_terms", 2),
            ("domains", 1),
            ("end", 1),
            ("flag", 1),
            ("flags", 1),
            ("id", 2),
            ("impact", 1),
            ("seq_region_name", 1),
            ("start", 2),
            ("strand", 1),
            ("transcript_id", 2),
            ("variant_allele", 3),
        ])
    );
    assert_eq!(
        counts(&r.counters, "kept_populated"),
        pairs(&[
            ("allele_string", 1),
            ("colocated_variants", 1),
            ("consequence_terms", 2),
            ("end", 1),
            ("flag", 1),
            ("id", 1),
            ("impact", 1),
            ("seq_region_name", 1),
            ("start", 2),
            ("strand", 1),
            ("transcript_id", 2),
            ("variant_allele", 2),
        ])
    );
}

#[test]
fn a_floor_is_counted_under_its_fields_in_the_consequence_entries() {
    let mut rows = ten_records();
    // A record's own member of the name is no consequence entry, and `-` is no value.
    rows[0] = object_with(0, &[TX1], r#","exon":"1/1""#);
    rows[1] = rows[1].replace("\"2/5\"", "\"-\"");
    let scratch = Scratch::new();
    let floors = scratch.file(
        "floors.txt",
        b"EXON\texon\nIMPACT\tseverity\timpact\nSYMBOL\tgene_symbol\tid\n",
    );
    let more = ["--floor-fields", floors.to_str().unwrap()];
    let r = remainder(&rows, &rows, &ten_inputs(), &more);
    assert!(r.expected.is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "29");
    assert_eq!(
        counts(&r.counters, "kept_under"),
        pairs(&[("EXON", 8), ("IMPACT", 19)])
    );
    assert_eq!(counter(&r.counters, "exon"), "9");
    let mut actual = rows.clone();
    actual[9] = object(9, &[TX2]);
    let r = remainder(&rows, &actual, &ten_inputs(), &more);
    assert_eq!(data(&r.expected), rows[9..]);
    assert_eq!(
        counts(&r.counters, "kept_under"),
        pairs(&[("EXON", 7), ("IMPACT", 17)])
    );
}

#[test]
fn a_floor_that_finds_no_field_hands_over_the_lines_it_may_be_counted_in() {
    let mut rows = ten_records();
    // An item's member that is a list of objects, an item's member named with a dot, and
    // an object with custom annotations.
    let domains =
        r#"{"variant_allele":"G","transcript_id":"tx3","domains":[{"db":"Pfam","name":"PF1"}]}"#;
    let dotted = r#"{"variant_allele":"G","transcript_id":"tx3","plugin.score":"0.5"}"#;
    rows[3] = object(3, &[TX1, domains]);
    rows[4] = object(4, &[TX1, dotted]);
    rows[5] = object_with(
        5,
        &[TX1],
        r#","custom_annotations":{"source":[{"name":"x"}]}"#,
    );
    let scratch = Scratch::new();
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert!(r.expected.is_empty());
    // Every floor finds a field in some entry.
    let found = scratch.file("found.txt", b"EXON\texon\nID\tid\n");
    let r = remainder(
        &rows,
        &rows,
        &ten_inputs(),
        &["--floor-fields", found.to_str().unwrap()],
    );
    assert!(r.expected.is_empty());
    assert_eq!(counts(&r.counters, "kept_under"), pairs(&[("EXON", 7)]));
    // One floor finds none.
    let lost = scratch.file("lost.txt", b"EXON\texon\nSCORE\tscore\tprediction\n");
    let more = ["--floor-fields", lost.to_str().unwrap()];
    let r = remainder(&rows, &rows, &ten_inputs(), &more);
    assert_eq!(data(&r.expected), rows[3..6]);
    assert_eq!(data(&r.actual), rows[3..6]);
    assert_eq!(counts(&r.counters, "kept_under"), pairs(&[("EXON", 7)]));
    assert_eq!(counter(&r.counters, "whole"), "null");
    // A field one line of the candidate alone carries is a field found.
    let mut actual = rows.clone();
    actual[9] = actual[9].replace("\"exon\":", "\"score\":");
    let r = remainder(&rows, &actual, &ten_inputs(), &more);
    assert_eq!(data(&r.expected), rows[9..]);
    assert_eq!(data(&r.actual), actual[9..]);
}

#[test]
fn a_reference_entry_of_a_kind_off_the_kept_list_carries_no_field_of_a_floor() {
    let motif = r#","motif_feature_consequences":[{"variant_allele":"G","motif_feature_id":"m1","motif_name":"M1"}]"#;
    let domains =
        r#"{"variant_allele":"G","transcript_id":"tx3","domains":[{"db":"Pfam","name":"PF1"}]}"#;
    let mut expected = ten_records();
    expected[2] = object(2, &[TX1, domains]);
    let mut actual = expected.clone();
    expected[8] = object_with(8, &[TX1], motif);
    actual[8] = object(8, &[TX1]);
    let scratch = Scratch::new();
    let floors = scratch.file("floors.txt", b"MOTIF_NAME\tmotif_name\n");
    let kinds = scratch.file("kinds.txt", b"transcript\nintergenic\n");
    let floor = ["--floor-fields", floors.to_str().unwrap()];
    // The reference's motif entry carries the field: the floor is found.
    let r = remainder(&expected, &actual, &ten_inputs(), &floor);
    assert_eq!(data(&r.expected), expected[8..9]);
    // Left out of the comparison, it carries none, and the candidate's entry of the kind
    // would.
    let more = [
        floor[0],
        floor[1],
        "--keep-values",
        "kind",
        kinds.to_str().unwrap(),
    ];
    let r = remainder(&expected, &actual, &ten_inputs(), &more);
    assert_eq!(
        data(&r.expected),
        [&expected[2..3], &expected[8..9]].concat()
    );
    let r = remainder(&actual, &expected, &ten_inputs(), &more);
    assert_eq!(data(&r.expected), actual[8..9]);
}

#[test]
fn shapes_a_line_s_text_does_not_settle_hand_everything_over() {
    let rows = ten_records();
    let whole = |line: &str, reason: &str| {
        let text = format!("{}{line}\n", output(&rows));
        let r = remainder_of_texts(&text, &text, &ten_inputs(), &[]);
        assert_eq!(r.expected, text, "{line}");
        assert_eq!(r.actual, text);
        assert_eq!(counter(&r.counters, "whole"), format!("\"{reason}\""));
        assert_eq!(counter(&r.counters, "kept_entries"), "0");
        assert_eq!(counter(&r.counters, "kept_populated"), "{}");
        assert_eq!(counter(&r.counters, "kept_present"), "{}");
        assert_eq!(counter(&r.counters, "kept_under"), "{}");
        assert_eq!(counter(&r.counters, "groups"), "0");
    };
    for line in [
        "## a note",
        "[1, 2]",
        "{\"id\": 1} {\"id\": 2}",
        "{\"id\": 01}",
        "{\"id\": 'rs1'}",
        "{\"id\": \"a\tb\"}",
        "{\"id\": NaN}",
        "{\"id\": 1,}",
    ] {
        whole(line, "a line that is not one JSON object");
    }
    whole(
        r#"{"id":"rs1","start":1,"id":"rs2"}"#,
        "a JSON object naming a member twice",
    );
    whole(
        r#"{"colocated_variants":[{"id":"rs1","id":"rs2"}]}"#,
        "a JSON object naming a member twice",
    );
    whole(r#"{"id":"\ud800"}"#, "a JSON escape of a lone surrogate");
    whole(
        r#"{"transcript_consequences":[{"variant_allele":"G"},"tx"]}"#,
        "a consequence list item that is not an object",
    );
    for item in [
        r#"{"variant_allele":"G|T","transcript_id":"tx1"}"#,
        r#"{"variant_allele":"G","transcript_id":"tx1#2"}"#,
    ] {
        whole(
            &object(3, &[item]),
            "a variant_allele or feature id value holds | or #",
        );
    }
}

#[test]
fn a_mark_outside_the_key_members_settles_its_entry() {
    let mut rows = ten_records();
    rows[1] = rows[1].replace("\"2/5\"", "\"2#5|6\"");
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
}

#[test]
fn every_line_end_and_gzipped_files_read_as_the_plain_text() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual[6] = object(6, &[TX1]);
    let plain = remainder(&expected, &actual, &ten_inputs(), &[]);
    assert_eq!(data(&plain.expected), expected[6..7]);
    let gz = |text: &str| {
        let mut gz = GzEncoder::new(Vec::new(), Compression::default());
        gz.write_all(text.as_bytes()).unwrap();
        gz.finish().unwrap()
    };
    let scratch = Scratch::new();
    let e = scratch.file("e.json.gz", &gz(&output(&expected).replace('\n', "\r\n")));
    let a = scratch.file("a.json", output(&actual).replace('\n', "\r").as_bytes());
    // The input's header lines are no records.
    let input = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\n{}",
        output(&ten_inputs())
    );
    let i = scratch.file("i.vcf.gz", &gz(&input.replace('\n', "\r\n")));
    let r = remainder_of(&e, &a, &i, &[]);
    assert_eq!(r.expected, plain.expected);
    assert_eq!(r.actual, plain.actual);
    assert_eq!(r.counters, plain.counters);
}

#[test]
fn many_lines_split_alike_at_every_worker_count() {
    // 3,000 records of up to four items; every 97th record differs in one item, every
    // 211th is on the reference side alone, every 389th record's key returns 40 lines
    // later, the input holding its line twice, and a blank line follows every 50th line of
    // the candidate.
    let item = |t: usize, impact: &str| {
        format!(r#"{{"variant_allele":"G","transcript_id":"tx{t}","impact":"{impact}"}}"#)
    };
    let entries = |row: &String| 1 + row.matches("\"transcript_id\"").count();
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    let mut inputs = Vec::new();
    let mut late: Vec<(usize, usize)> = Vec::new();
    let mut differing = Vec::new();
    for r in 0..3000usize {
        while late.first().is_some_and(|(at, _)| *at <= r) {
            let (_, record) = late.remove(0);
            let row = object(record, &[&item(9, "MODIFIER")]);
            inputs.push(input_line(record));
            expected.push(row.clone());
            actual.push(row);
        }
        inputs.push(input_line(r));
        let items: Vec<String> = (0..1 + r % 4).map(|t| item(t, "MODIFIER")).collect();
        let items: Vec<&str> = items.iter().map(String::as_str).collect();
        let row_e = object(r, &items);
        let row_a = if r % 97 == 5 {
            row_e.replacen("MODIFIER", "LOW", 1)
        } else {
            row_e.clone()
        };
        if r % 97 == 5 || r % 211 == 9 {
            differing.push(row_e.clone());
        }
        expected.push(row_e);
        if r % 211 != 9 {
            actual.push(row_a);
        }
        if r % 50 == 49 {
            actual.push(String::new());
        }
        if r % 389 == 11 {
            late.push((r + 40, r));
        }
    }
    let scratch = Scratch::new();
    let e = scratch.file("expected.json", output(&expected).as_bytes());
    let a = scratch.file("actual.json", output(&actual).as_bytes());
    let i = scratch.file("input.vcf", output(&inputs).as_bytes());
    let three = remainder_of(&e, &a, &i, &[]);
    let seven = remainder_at(Format::Json, &e, &a, &i, &[], 7, Some(4096));
    assert_eq!(texts_of(Format::Json, &seven), three);
    assert_eq!(counter(&three.counters, "whole"), "null");
    assert_eq!(counter(&three.counters, "groups"), "3000");
    let handed: usize = counter(&three.counters, "groups_handed_over")
        .parse()
        .unwrap();
    assert!(
        handed >= differing.len() && handed < 3000 / 4,
        "{handed} groups handed over"
    );
    // Every entry is kept back or on the remainder, on each side.
    let kept: usize = counter(&three.counters, "kept_entries").parse().unwrap();
    let on = |rows: &[String]| rows.iter().map(entries).sum::<usize>();
    assert_eq!(kept + on(&data(&three.expected)), on(&expected));
    let written: Vec<String> = actual
        .iter()
        .filter(|row| !row.is_empty())
        .cloned()
        .collect();
    assert_eq!(kept + on(&data(&three.actual)), on(&written));
    let handed_rows = data(&three.expected);
    for row in &differing {
        assert!(handed_rows.contains(row), "a differing line was kept back");
    }
    // The runs of the candidate's lines give each remainder line its place in the output,
    // blank lines counted.
    let runs = counter(&three.counters, "line_runs");
    let runs = &runs[runs.find("\"actual\": [[").unwrap() + "\"actual\": [[".len()..];
    let runs = &runs[..runs.find("]]}").unwrap()];
    let handed_rows = data(&three.actual);
    let mut at = 0;
    for run in runs.split("], [") {
        let (n, before) = run.split_once(", ").unwrap();
        let (n, before): (usize, usize) = (n.parse().unwrap(), before.parse().unwrap());
        for row in &handed_rows[at..at + n] {
            assert_eq!(&actual[at + before], row);
            at += 1;
        }
    }
    assert_eq!(at, handed_rows.len());
}

#[test]
fn a_plan_that_names_what_its_format_does_not_read_is_refused() {
    let scratch = Scratch::new();
    let rows = ten_records();
    let e = scratch.file("e.json", output(&rows).as_bytes());
    let i = scratch.file("i.vcf", output(&ten_inputs()).as_bytes());
    let refused = |plan: &stream::Plan| {
        let error = stream::remainder(&e, &e, plan).err().unwrap();
        error.to_string()
    };
    let floor = || vec![("EXON".to_string(), vec!["exon".to_string()])];
    let mut lacking = plan_of(Format::Json, &i, &[], 2);
    lacking.input = None;
    assert_eq!(
        refused(&lacking),
        "the json format needs the input the outputs were written from"
    );
    let mut located = plan_of(Format::Json, &i, &[], 2);
    located.force_locations = vec!["1:1000".to_string()];
    assert_eq!(
        refused(&located),
        "the json format reads no listed Locations"
    );
    let mut kept = plan_of(Format::Json, &i, &[], 2);
    kept.keep_values = Some(("Feature_type".to_string(), vec!["Transcript".to_string()]));
    assert_eq!(
        refused(&kept),
        "the column of the kept values of the json format is kind"
    );
    let mut vcf = plan_of(Format::Vcf, &i, &[], 2);
    vcf.floor_fields = floor();
    assert_eq!(refused(&vcf), "the vcf format reads no floors");
    for (format, name) in [(Format::Default, "default"), (Format::Tab, "tab")] {
        let mut floored = plan_of(format, &i, &[], 2);
        floored.input = None;
        floored.floor_fields = floor();
        assert_eq!(
            refused(&floored),
            format!("the {name} format reads no floors")
        );
    }
}
