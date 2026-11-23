// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reads planted pairs of VCF outputs and their input, every line written here, through
//! the streamed record groups at one worker, two and the host's default and through blocks
//! of several sizes, and holds the readings to one remainder and to the remainder each
//! shape must leave.

mod stream_records_common;

use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use std::path::Path;
use stream_records_common::{
    counter, line_runs, pairs, plan_of, populated, remainder_at, texts_of, Remainder, Scratch,
};
use vep_formats_compare::stream::{self, Format};

const VERSION: &str = "##fileformat=VCFv4.2";
const CSQ: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|Feature_type|Feature|EXON\">";
const NAMES: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
/// Two chunks of the Format list of `CSQ`: an intron row and a row with an exon number.
const TX1: &str = "G|intron_variant|MODIFIER|Transcript|tx1|";
const TX2: &str = "G|missense_variant|MODERATE|Transcript|tx2|2/5";

/// The first five columns of record `record`'s line.
fn key(record: usize) -> String {
    format!("1\t{}\trs{record}\tA\tG", 1000 + record)
}

/// Record `record`'s line of the input.
fn input_line(record: usize) -> String {
    format!("{}\t.\t.\t.", key(record))
}

/// Record `record`'s line of an output, its CSQ the `chunks`.
fn line(record: usize, chunks: &[&str]) -> String {
    format!("{}\t.\tPASS\tCSQ={}", key(record), chunks.join(","))
}

fn lines(rows: &[String]) -> String {
    rows.iter().map(|row| format!("{row}\n")).collect()
}

fn output(rows: &[String]) -> String {
    format!("{VERSION}\n{CSQ}\n{NAMES}\n{}", lines(rows))
}

/// Ten records, two chunks each.
fn ten_records() -> Vec<String> {
    (0..10).map(|r| line(r, &[TX1, TX2])).collect()
}

fn ten_inputs() -> Vec<String> {
    (0..10).map(input_line).collect()
}

fn remainder_of(expected: &Path, actual: &Path, input: &Path, more: &[&str]) -> Remainder {
    stream_records_common::remainder_of(Format::Vcf, expected, actual, input, more)
}

fn remainder_of_texts(expected: &str, actual: &str, inputs: &[String], more: &[&str]) -> Remainder {
    let scratch = Scratch::new();
    let e = scratch.file("expected.vcf", expected.as_bytes());
    let a = scratch.file("actual.vcf", actual.as_bytes());
    let i = scratch.file("input.vcf", lines(inputs).as_bytes());
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

/// The data lines of a remainder text, its header lines checked and dropped.
fn data(text: &str) -> Vec<String> {
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some(VERSION));
    assert_eq!(lines.next(), Some(CSQ));
    assert_eq!(lines.next(), Some(NAMES));
    lines.map(str::to_string).collect()
}

#[test]
fn outputs_written_alike_leave_their_header_lines_alone() {
    let rows = ten_records();
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert!(data(&r.expected).is_empty() && data(&r.actual).is_empty());
    assert_eq!(counter(&r.counters, "whole"), "null");
    // A line is its record entry and one entry a chunk.
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
    assert_eq!(
        populated(&r.counters),
        pairs(&[
            ("ALT", 10),
            ("Allele", 20),
            ("CHROM", 10),
            ("CSQ_present", 10),
            ("Consequence", 20),
            ("EXON", 10),
            ("FILTER", 10),
            ("Feature", 20),
            ("Feature_type", 20),
            ("ID", 10),
            ("IMPACT", 20),
            ("POS", 10),
            ("QUAL", 10),
            ("REF", 10),
        ])
    );
    assert_eq!(counter(&r.counters, "kept_before"), "[]");
    assert_eq!(counter(&r.counters, "line_runs"), line_runs("[]", "[]"));
    assert_eq!(counter(&r.counters, "groups"), "10");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "0");
    assert_eq!(
        counter(&r.counters, "expected"),
        "{\"lines\": 10, \"lines_handed_over\": 0, \"bytes_handed_over\": 0}"
    );
    // A record of the input that neither output writes is no group.
    let r = remainder(&rows[..4], &rows[1..5], &ten_inputs(), &[]);
    assert_eq!(counter(&r.counters, "groups"), "5");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "2");
}

#[test]
fn a_differing_chunk_hands_its_line_over() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual[3] = line(3, &[TX1, "G|missense_variant|HIGH|Transcript|tx2|2/5"]);
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
    assert_eq!(counter(&r.counters, "IMPACT"), "18");
    assert_eq!(counter(&r.counters, "CSQ_present"), "9");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "1");
    assert_eq!(
        counter(&r.counters, "expected"),
        format!(
            "{{\"lines\": 10, \"lines_handed_over\": 1, \"bytes_handed_over\": {}}}",
            expected[3].len() + 1
        )
    );
    // A column the comparison reads no entry from is text of the line all the same.
    let mut actual = expected.clone();
    actual[8].push_str("\tGT\t0/1");
    let r = remainder(&expected, &actual, &ten_inputs(), &[]);
    assert_eq!(data(&r.expected), expected[8..9]);
    assert_eq!(data(&r.actual), actual[8..9]);
}

#[test]
fn a_line_one_side_alone_writes_is_handed_over_from_either_side() {
    let both = ten_records();
    let mut inputs = ten_inputs();
    inputs.push(input_line(10));
    let mut longer = both.clone();
    longer.insert(4, line(10, &[TX1]));
    let r = remainder(&longer, &both, &inputs, &[]);
    assert_eq!(data(&r.expected), longer[4..5]);
    assert!(data(&r.actual).is_empty());
    assert_eq!(counter(&r.counters, "kept_before"), "[[2, 12]]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[1, 4]]", "[]")
    );
    let r = remainder(&both, &longer, &inputs, &[]);
    assert!(data(&r.expected).is_empty());
    assert_eq!(data(&r.actual), longer[4..5]);
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
    assert_eq!(counter(&r.counters, "kept_before"), "[]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[]", "[[1, 4]]")
    );
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
    assert!(data(&r.expected).is_empty() && data(&r.actual).is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "33");
    assert_eq!(counter(&r.counters, "groups"), "10");
    // One side with the line once is another count of lines.
    let r = remainder(&rows, &ten_records(), &inputs, &[]);
    assert_eq!(data(&r.expected), rows[2..4]);
    assert_eq!(data(&r.actual), rows[2..3]);
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
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[1, 7]]", "[[1, 7]]")
    );
    // The key is the first five columns as text: another ID is another key.
    let mut inputs = ten_inputs();
    inputs[7] = inputs[7].replace("rs7", "rs70");
    let r = remainder(&rows, &rows, &inputs, &[]);
    assert_eq!(data(&r.expected), rows[7..8]);
    // The input's columns after the fifth are no part of it.
    let mut inputs = ten_inputs();
    inputs[7] = format!("{}\t50\tq10\tAC=2\tGT\t0/1", key(7));
    let r = remainder(&rows, &rows, &inputs, &[]);
    assert!(data(&r.expected).is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
}

#[test]
fn a_listed_record_is_handed_over_with_equal_lines() {
    let rows = ten_records();
    let scratch = Scratch::new();
    // The ordinals count the input's data lines from 0; one past the last names no line.
    let list = scratch.file("records.txt", b"4\n99\n4\n");
    let more = ["--force-records", list.to_str().unwrap()];
    let mut inputs = vec![VERSION.to_string(), NAMES.to_string()];
    inputs.extend(ten_inputs());
    let r = remainder(&rows, &rows, &inputs, &more);
    assert_eq!(data(&r.expected), rows[4..5]);
    assert_eq!(data(&r.actual), rows[4..5]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 12]]");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "1");
    // A listed line of a key the input holds twice hands over both lines of the key.
    let mut inputs = ten_inputs();
    inputs.push(input_line(2));
    let mut twice = rows.clone();
    twice.push(rows[2].clone());
    let list = scratch.file("second.txt", b"10\n");
    let more = ["--force-records", list.to_str().unwrap(), "--unordered"];
    let r = remainder(&twice, &twice, &inputs, &more);
    assert_eq!(data(&r.expected), [&twice[2..3], &twice[10..]].concat());
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 6], [3, 27]]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[1, 2], [1, 9]]", "[[1, 2], [1, 9]]")
    );
    // Read with its order, what stands between the two lines is handed over with them.
    let r = remainder(&twice, &twice, &inputs, &more[..2]);
    assert_eq!(data(&r.expected), twice[2..]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[27, 6]]");
    let r = remainder(&twice, &twice, &inputs, &[]);
    assert!(data(&r.expected).is_empty());
}

#[test]
fn a_reference_chunk_of_a_value_off_the_kept_list_hands_its_line_over() {
    let mut rows = ten_records();
    rows[5] = line(
        5,
        &[
            TX1,
            "G|regulatory_region_variant|MODIFIER|RegulatoryFeature|reg1|",
        ],
    );
    // A line without a CSQ has no entry of the column.
    rows[6] = format!("{}\t.\tPASS\t.", key(6));
    let scratch = Scratch::new();
    let list = scratch.file("kept.txt", b"Transcript\n-\n\n");
    let more = ["--keep-values", "Feature_type", list.to_str().unwrap()];
    let r = remainder(&rows, &rows, &ten_inputs(), &more);
    assert_eq!(data(&r.expected), rows[5..6]);
    assert_eq!(data(&r.actual), rows[5..6]);
    assert_eq!(counter(&r.counters, "kept_entries"), "25");
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert!(data(&r.expected).is_empty());
    // A Format list without the column reads every chunk's value as the empty text.
    let narrow = CSQ.replace("|Feature_type|", "|Kind|");
    let text = output(&rows).replace(CSQ, &narrow);
    let r = remainder_of_texts(&text, &text, &ten_inputs(), &more);
    assert_eq!(counter(&r.counters, "kept_entries"), "28");
    let list = scratch.file("valued.txt", b"Transcript\n-\n");
    let more = ["--keep-values", "Feature_type", list.to_str().unwrap()];
    let r = remainder_of_texts(&text, &text, &ten_inputs(), &more);
    assert_eq!(counter(&r.counters, "kept_entries"), "1");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "9");
}

#[test]
fn a_reference_line_holding_a_listed_text_hands_its_line_over() {
    let mut rows = ten_records();
    rows[6] = line(
        6,
        &[TX1, "G|intron_variant|MODIFIER|Transcript|tx2|misspeltKey7"],
    );
    let scratch = Scratch::new();
    let list = scratch.file("texts.txt", b"misspeltKey\n");
    let more = ["--force-text", list.to_str().unwrap()];
    let r = remainder(&rows, &rows, &ten_inputs(), &more);
    assert_eq!(data(&r.expected), rows[6..7]);
    assert_eq!(data(&r.actual), rows[6..7]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "27");
    // The candidate's lines and header lines are not rewritten.
    let plain = ten_records();
    let noted = output(&plain).replace(VERSION, "##source=misspeltKey");
    let r = remainder_of_texts(&output(&plain), &noted, &ten_inputs(), &more);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
}

#[test]
fn a_listed_text_in_a_reference_header_line_or_a_key_hands_everything_over() {
    let reason = "\"a listed text stands in a header line of the reference or in a line's key\"";
    let scratch = Scratch::new();
    let list = scratch.file("texts.txt", b"misspeltKey\n");
    let more = ["--force-text", list.to_str().unwrap()];
    let rows = ten_records();
    // In the Format list of both outputs.
    let text = output(&rows).replace("|EXON\">", "|misspeltKey100way\">");
    let r = remainder_of_texts(&text, &text, &ten_inputs(), &more);
    assert_eq!(counter(&r.counters, "whole"), reason);
    assert_eq!(r.expected, text);
    assert_eq!(r.actual, text);
    assert_eq!(counter(&r.counters, "kept_entries"), "0");
    assert_eq!(counter(&r.counters, "kept_before"), "[[30, 0]]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[10, 0]]", "[[10, 0]]")
    );
    let r = remainder_of_texts(&text, &text, &ten_inputs(), &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "misspeltKey100way"), "10");
    // In a line's ID, a column of its key.
    let mut named = rows.clone();
    named[1] = named[1].replace("rs1", "misspeltKey1");
    let mut inputs = ten_inputs();
    inputs[1] = inputs[1].replace("rs1", "misspeltKey1");
    let r = remainder(&named, &named, &inputs, &more);
    assert_eq!(counter(&r.counters, "whole"), reason);
    assert_eq!(data(&r.expected), named);
}

#[test]
fn lines_of_one_key_that_stand_apart_are_one_group() {
    // Record 2's key returns after record 4's line, and the input holds its line twice.
    let mut expected = ten_records();
    expected.insert(5, line(2, &[TX1]));
    let mut inputs = ten_inputs();
    inputs.insert(5, input_line(2));
    let r = remainder(&expected, &expected, &inputs, &[]);
    assert!(data(&r.expected).is_empty());
    assert_eq!(counter(&r.counters, "groups"), "10");
    assert_eq!(counter(&r.counters, "kept_entries"), "32");
    // A difference in the late line hands over its group and what stands inside its span:
    // the lines of records 3 and 4.
    let mut actual = expected.clone();
    actual[5] = line(2, &[TX2]);
    let r = remainder(&expected, &actual, &inputs, &[]);
    assert_eq!(data(&r.expected), expected[2..6]);
    assert_eq!(data(&r.actual), actual[2..6]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[11, 6]]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[4, 2]]", "[[4, 2]]")
    );
    assert_eq!(counter(&r.counters, "groups_handed_over"), "3");
    // Read without order, the group alone is handed over.
    let r = remainder(&expected, &actual, &inputs, &["--unordered"]);
    assert_eq!(
        data(&r.expected),
        [&expected[2..3], &expected[5..6]].concat()
    );
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 6], [2, 12]]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[1, 2], [1, 4]]", "[[1, 2], [1, 4]]")
    );
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
    assert!(data(&r.expected).is_empty() && data(&r.actual).is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
    // Two lines of one key in another order are another text in order.
    let mut inputs = ten_inputs();
    inputs.push(input_line(9));
    let mut expected = ten_records();
    expected.push(line(9, &[TX1]));
    let mut actual = expected.clone();
    actual.swap(9, 10);
    let r = remainder(&expected, &actual, &inputs, &["--unordered"]);
    assert_eq!(data(&r.expected), expected[9..]);
    assert_eq!(data(&r.actual), actual[9..]);
}

#[test]
fn lines_labelled_by_their_place_are_one_group_in_the_order_of_the_shared_keys() {
    // Records 11 to 13 are no records of the input: the comparison labels their lines by
    // their place in the output, and two such lines at one place share their keys.
    let inputs = vec![input_line(0)];
    let known = line(0, &[TX1]);
    let (first, second, third) = (line(11, &[TX1]), line(12, &[TX1]), line(13, &[TX1]));
    // The reference's line at place 1 follows the known line, the candidate's precedes
    // it: the known line is no line to keep back.
    let expected = vec![known.clone(), first.clone()];
    let actual = vec![second.clone(), third.clone(), known.clone()];
    let r = remainder(&expected, &actual, &inputs, &[]);
    assert_eq!(data(&r.expected), expected);
    assert_eq!(data(&r.actual), actual);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "0");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "4");
    let r = remainder(&expected, &actual, &inputs, &["--unordered"]);
    assert_eq!(data(&r.expected), expected[1..]);
    assert_eq!(data(&r.actual), actual[..2]);
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[1, 1]]", "[[2, 0]]")
    );
    // Such lines on one side alone share no key with the other side.
    let r = remainder(&expected, &expected[..1], &inputs, &[]);
    assert_eq!(data(&r.expected), expected[1..]);
    assert!(data(&r.actual).is_empty());
    // On both sides after the same kept-back lines, they move nothing.
    let actual = vec![known.clone(), second.clone()];
    let r = remainder(&expected, &actual, &inputs, &[]);
    assert_eq!(data(&r.expected), expected[1..]);
    assert_eq!(data(&r.actual), actual[1..]);
    assert_eq!(counter(&r.counters, "kept_entries"), "2");
    // A line past the input's lines of its key is such a line: the second line of
    // record 0 stands after record 1's line on the reference side and before it on the
    // candidate's.
    let inputs = vec![input_line(0), input_line(1)];
    let other = line(1, &[TX2]);
    let expected = vec![known.clone(), other.clone(), known.clone()];
    let actual = vec![known.clone(), known.clone(), other.clone()];
    let r = remainder(&expected, &actual, &inputs, &[]);
    assert_eq!(data(&r.expected), expected);
    assert_eq!(data(&r.actual), actual);
    let r = remainder(&expected, &actual, &inputs, &["--unordered"]);
    assert_eq!(data(&r.expected), [&expected[..1], &expected[2..]].concat());
    assert_eq!(data(&r.actual), actual[..2]);
}

#[test]
fn a_line_the_comparison_may_fail_to_convert_hands_its_line_over() {
    let plain = ten_records();
    for exon in [
        // A digit outside ASCII, a run of 641 digits and 65 opening brackets.
        "\u{b2}".to_string(),
        "7".repeat(641),
        format!("{}1", "{\"a\":".repeat(65)),
    ] {
        let mut rows = plain.clone();
        let chunk = format!("G|intron_variant|MODIFIER|Transcript|tx2|{exon}");
        rows[4] = line(4, &[TX1, &chunk]);
        let r = remainder(&rows, &rows, &ten_inputs(), &[]);
        assert_eq!(data(&r.expected), rows[4..5]);
        assert_eq!(data(&r.actual), rows[4..5]);
        assert_eq!(counter(&r.counters, "kept_entries"), "27");
    }
    // At the two bounds the line is kept back.
    let mut rows = plain.clone();
    let chunk = format!(
        "G|intron_variant|MODIFIER|Transcript|tx2|{}{}",
        "7".repeat(640),
        "[".repeat(64)
    );
    rows[4] = line(4, &[TX1, &chunk]);
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert!(data(&r.expected).is_empty() && data(&r.actual).is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
}

#[test]
fn every_field_of_a_record_entry_and_of_a_chunk_is_counted_by_its_value() {
    // No CSQ and two other INFO items; one chunk between two other items; three chunks;
    // and a record of an absent ID, an empty QUAL, an absent FILTER and one chunk of an
    // allele alone.
    let rows = vec![
        format!("{}\t.\t.\tAC=2;AN=4", key(0)),
        format!("{}\t.\tPASS\tAC=1;CSQ={TX2};AN=2", key(1)),
        line(2, &[TX1, TX2, TX1]),
        "1\t1003\t-\tA\tG\t\t-\tCSQ=G|||||".to_string(),
    ];
    let mut inputs: Vec<String> = (0..3).map(input_line).collect();
    inputs.push("1\t1003\t-\tA\tG\t.\t.\t.".to_string());
    let r = remainder(&rows, &rows, &inputs, &[]);
    assert!(data(&r.expected).is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "9");
    let mut actual = rows.clone();
    actual[2] = line(2, &[TX1, TX2]);
    let r = remainder(&rows, &actual, &inputs, &[]);
    assert_eq!(data(&r.expected), rows[2..3]);
    assert_eq!(data(&r.actual), actual[2..3]);
    assert_eq!(counter(&r.counters, "kept_entries"), "5");
    // Four entries after the three of the first two lines.
    assert_eq!(counter(&r.counters, "kept_before"), "[[4, 3]]");
    assert_eq!(
        counter(&r.counters, "line_runs"),
        line_runs("[[1, 2]]", "[[1, 2]]")
    );
    assert_eq!(
        populated(&r.counters),
        pairs(&[
            ("ALT", 3),
            ("Allele", 2),
            ("CHROM", 3),
            ("CSQ_present", 3),
            ("Consequence", 1),
            ("EXON", 1),
            ("FILTER", 2),
            ("Feature", 1),
            ("Feature_type", 1),
            ("ID", 2),
            ("IMPACT", 1),
            ("INFO_without_CSQ", 2),
            ("POS", 3),
            ("QUAL", 2),
            ("REF", 3),
        ])
    );
}

#[test]
fn the_last_csq_item_of_a_line_is_its_csq() {
    let mut rows = ten_records();
    rows[0] = format!("{}\t.\tPASS\tCSQ=first|chunk;CSQ={TX1},{TX2}", key(0));
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
    assert_eq!(counter(&r.counters, "Allele"), "20");
}

#[test]
fn shapes_a_line_s_text_does_not_settle_hand_everything_over() {
    let rows = ten_records();
    let whole = |expected: String, actual: String, inputs: &[String], reason: &str| {
        let r = remainder_of_texts(&expected, &actual, inputs, &[]);
        assert_eq!(r.expected, expected);
        assert_eq!(r.actual, actual);
        assert_eq!(counter(&r.counters, "whole"), format!("\"{reason}\""));
        assert_eq!(counter(&r.counters, "kept_entries"), "0");
        assert_eq!(counter(&r.counters, "kept_populated"), "{}");
        assert_eq!(counter(&r.counters, "groups"), "0");
        r.counters
    };
    let same = |text: String, reason: &str| whole(text.clone(), text, &ten_inputs(), reason);
    let counters = same(
        output(&rows) + NAMES + "\n",
        "a header line follows a data line",
    );
    assert_eq!(counter(&counters, "kept_before"), "[[30, 0]]");
    assert_eq!(
        counter(&counters, "line_runs"),
        line_runs("[[10, 0]]", "[[10, 0]]")
    );
    whole(
        output(&rows),
        output(&rows).replace("|EXON\">", "|INTRON\">"),
        &ten_inputs(),
        "the two outputs' CSQ Format lists differ",
    );
    whole(
        output(&rows),
        output(&rows).replace(&format!("{CSQ}\n"), ""),
        &ten_inputs(),
        "the two outputs' CSQ Format lists differ",
    );
    same(
        output(&rows).replace(" Format: ", " Fields: "),
        "a CSQ header line has no Format list",
    );
    same(
        output(&rows).replace("|EXON\">", "|Allele\">"),
        "the CSQ Format list names a field twice",
    );
    same(
        output(&rows) + &key(3) + "\t.\tPASS\n",
        "a data line has fewer than 8 columns",
    );
    same(output(&rows) + "\n", "a data line has fewer than 8 columns");
    let narrow = line(4, &[TX1, "G|intron_variant|MODIFIER|Transcript|tx2"]);
    same(
        output(&[rows[..4].to_vec(), vec![narrow]].concat()),
        "a CSQ chunk's width is not the Format list's",
    );
    // No Format list at all is a list of no field, which no chunk is as wide as.
    same(
        output(&rows).replace(&format!("{CSQ}\n"), ""),
        "a CSQ chunk's width is not the Format list's",
    );
    for chunk in [
        "G#2|intron_variant|MODIFIER|Transcript|tx2|",
        "G|intron_variant|MODIFIER|Transcript#2|tx2|",
        "G|intron_variant|MODIFIER|Transcript|tx2#2|",
    ] {
        same(
            output(&[line(0, &[TX1, chunk])]),
            "an Allele, Feature or Feature_type value holds #",
        );
    }
}

#[test]
fn a_mark_outside_the_key_fields_and_a_line_without_a_csq_settle_their_entries() {
    let mut rows = ten_records();
    rows[1] = line(1, &[TX1, "G|intron_variant#2|MODIFIER|Transcript|tx2|2#5"]);
    rows[2] = format!("{}\t.\tPASS\t.", key(2));
    let r = remainder(&rows, &rows, &ten_inputs(), &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "28");
    // Without a Format list a line without a CSQ is its record entry.
    let bare: Vec<String> = (0..10).map(|r| format!("{}\t.\tPASS\t.", key(r))).collect();
    let text = output(&bare).replace(&format!("{CSQ}\n"), "");
    let r = remainder_of_texts(&text, &text, &ten_inputs(), &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "10");
    assert_eq!(counter(&r.counters, "INFO_without_CSQ"), "10");
}

#[test]
fn header_lines_may_differ_where_the_format_list_is_one() {
    let rows = ten_records();
    let other = output(&rows)
        .replace(VERSION, "##fileformat=VCFv4.3\n##source=the candidate")
        .replace(
            "from Ensembl VEP. Format: ",
            "of another description. Format: ",
        );
    let r = remainder_of_texts(&output(&rows), &other, &ten_inputs(), &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert!(data(&r.expected).is_empty());
    assert_eq!(r.actual, other[..other.find("1\t1000").unwrap()]);
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
    // The last CSQ header line of each side is the one in force.
    let earlier = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: Allele|Gene\">";
    let twice = output(&rows).replace(CSQ, &format!("{earlier}\n{CSQ}"));
    let r = remainder_of_texts(&output(&rows), &twice, &ten_inputs(), &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "30");
}

#[test]
fn every_line_end_and_gzipped_files_read_as_the_plain_text() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual[6] = line(6, &[TX1]);
    let plain = remainder(&expected, &actual, &ten_inputs(), &[]);
    assert_eq!(data(&plain.expected), expected[6..7]);
    let gz = |text: &str| {
        let mut gz = GzEncoder::new(Vec::new(), Compression::default());
        gz.write_all(text.as_bytes()).unwrap();
        gz.finish().unwrap()
    };
    let scratch = Scratch::new();
    let e = scratch.file("e.vcf.gz", &gz(&output(&expected).replace('\n', "\r\n")));
    let a = scratch.file("a.vcf", output(&actual).replace('\n', "\r").as_bytes());
    // The input's header lines are no records.
    let input = format!("{VERSION}\n{NAMES}\n{}", lines(&ten_inputs()));
    let i = scratch.file("i.vcf.gz", &gz(&input.replace('\n', "\r\n")));
    let r = remainder_of(&e, &a, &i, &[]);
    assert_eq!(r.expected, plain.expected);
    assert_eq!(r.actual, plain.actual);
    assert_eq!(r.counters, plain.counters);
}

#[test]
fn many_lines_split_alike_at_every_worker_count() {
    // 3,000 records of up to four chunks; every 97th record differs in one chunk, every
    // 211th is on the reference side alone, and every 389th record's key returns 40 lines
    // later, the input holding its line twice.
    let chunk = |t: usize, impact: &str| format!("G|intron_variant|{impact}|Transcript|tx{t}|");
    let entries = |row: &String| 1 + row.matches("|Transcript|").count();
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    let mut inputs = Vec::new();
    let mut late: Vec<(usize, String)> = Vec::new();
    let mut differing = Vec::new();
    for r in 0..3000usize {
        while late.first().is_some_and(|(at, _)| *at <= r) {
            let (_, row) = late.remove(0);
            inputs.push(row[..row.find("\t.\tPASS").unwrap()].to_string() + "\t.\t.\t.");
            expected.push(row.clone());
            actual.push(row);
        }
        inputs.push(input_line(r));
        let chunks: Vec<String> = (0..1 + r % 4).map(|t| chunk(t, "MODIFIER")).collect();
        let chunks: Vec<&str> = chunks.iter().map(String::as_str).collect();
        let row_e = line(r, &chunks);
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
        if r % 389 == 11 {
            late.push((r + 40, line(r, &[&chunk(9, "MODIFIER")])));
        }
    }
    let scratch = Scratch::new();
    let e = scratch.file("expected.vcf", output(&expected).as_bytes());
    let a = scratch.file("actual.vcf", output(&actual).as_bytes());
    let i = scratch.file("input.vcf", lines(&inputs).as_bytes());
    let three = remainder_of(&e, &a, &i, &[]);
    let seven = remainder_at(Format::Vcf, &e, &a, &i, &[], 7, Some(4096));
    assert_eq!(texts_of(Format::Vcf, &seven), three);
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
    assert_eq!(kept + on(&data(&three.actual)), on(&actual));
    // A kept-back line is on neither remainder, so every line that differs is on the
    // reference's.
    let handed_rows = data(&three.expected);
    for row in &differing {
        assert!(handed_rows.contains(row), "a differing line was kept back");
    }
    // The runs of the reference's lines add up to its remainder and to what stands
    // before each of its lines.
    let runs = counter(&three.counters, "line_runs");
    let runs = &runs["{\"expected\": [[".len()..runs.find("]], \"actual\"").unwrap()];
    let mut at = 0;
    for run in runs.split("], [") {
        let (n, before) = run.split_once(", ").unwrap();
        let (n, before): (usize, usize) = (n.parse().unwrap(), before.parse().unwrap());
        for row in &handed_rows[at..at + n] {
            assert_eq!(&expected[at + before], row);
            at += 1;
        }
    }
    assert_eq!(at, handed_rows.len());
}

#[test]
fn a_plan_that_names_what_its_format_does_not_read_is_refused() {
    let scratch = Scratch::new();
    let rows = ten_records();
    let e = scratch.file("e.vcf", output(&rows).as_bytes());
    let i = scratch.file("i.vcf", lines(&ten_inputs()).as_bytes());
    let refused = |plan: &stream::Plan| {
        let error = stream::remainder(&e, &e, plan).err().unwrap();
        error.to_string()
    };
    let mut lacking = plan_of(Format::Vcf, &i, &[], 2);
    lacking.input = None;
    assert_eq!(
        refused(&lacking),
        "the vcf format needs the input the outputs were written from"
    );
    let mut located = plan_of(Format::Vcf, &i, &[], 2);
    located.force_locations = vec!["1:1000".to_string()];
    assert_eq!(
        refused(&located),
        "the vcf format reads no listed Locations"
    );
    let mut floored = plan_of(Format::Vcf, &i, &[], 2);
    floored.floor_fields = vec![("IMPACT".to_string(), vec!["impact".to_string()])];
    assert_eq!(refused(&floored), "the vcf format reads no floors");
    for (format, name) in [(Format::Default, "default"), (Format::Tab, "tab")] {
        let with_input = plan_of(format, &i, &[], 2);
        assert_eq!(
            refused(&with_input),
            format!("the {name} format reads no input")
        );
        let mut listed = plan_of(format, &i, &[], 2);
        listed.input = None;
        listed.force_records = vec![1];
        assert_eq!(
            refused(&listed),
            format!("the {name} format reads no listed records")
        );
    }
    let absent = scratch.0.join("absent.vcf");
    let unread = plan_of(Format::Vcf, &absent, &[], 2);
    assert_eq!(
        refused(&unread),
        format!(
            "cannot read {}: No such file or directory",
            absent.display()
        )
    );
}

#[test]
fn a_default_or_tab_remainder_has_its_line_runs_as_well() {
    let text = "## a note\n#Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tExtra\n\
                rec1\t1:7\tG\ttx1\tTranscript\tIMPACT=LOW\n\
                rec2\t1:8\tG\ttx1\tTranscript\tIMPACT=LOW\n";
    let other = text.replace(
        "1:8\tG\ttx1\tTranscript\tIMPACT=LOW",
        "1:8\tG\ttx1\tTranscript\t-",
    );
    let scratch = Scratch::new();
    let e = scratch.file("e.txt", text.as_bytes());
    let a = scratch.file("a.txt", other.as_bytes());
    let mut plan = plan_of(Format::Default, &e, &[], 2);
    plan.input = None;
    let default = stream::remainder(&e, &a, &plan).unwrap();
    plan.format = Format::Tab;
    assert_eq!(stream::remainder(&e, &a, &plan).unwrap(), default);
    assert_eq!(default.kept_before, [(1, 1)]);
    assert_eq!(default.expected_line_runs, [(1, 1)]);
    assert_eq!(default.actual_line_runs, [(1, 1)]);
    assert_eq!((default.actual_lines, default.actual_handed_over), (2, 1));
    assert_eq!(default.actual_bytes_handed_over, 28);
}
