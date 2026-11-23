// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reads planted pairs of outputs, every row written here, through the streamed record
//! groups at one worker, two and the host's default and through blocks of several sizes,
//! and holds the readings to one remainder and to the remainder each shape must leave.

use flate2::write::{DeflateEncoder, GzEncoder};
use flate2::Compression;
use std::fmt::Write as _;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};
use vep_formats_compare::error::Error;
use vep_formats_compare::stream::{self, Format, Plan};

const COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tExtra";
const NOTE: &str = "## written by the reference engine";

static SCRATCH: AtomicU32 = AtomicU32::new(0);

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-compare-stream-test-{}-{}",
            std::process::id(),
            SCRATCH.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    fn file(&self, name: &str, bytes: &[u8]) -> PathBuf {
        let path = self.0.join(name);
        fs::write(&path, bytes).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// One data row: a transcript row of record `record` at `location`, its allele one of the
/// four bases by `record`, with `extra` as its Extra column.
fn row(record: usize, location: &str, feature: &str, extra: &str) -> String {
    let allele = ["A", "C", "G", "T"][record % 4];
    format!("rec{record}\t{location}\t{allele}\tgene{record}\t{feature}\tTranscript\tintron_variant\t{extra}")
}

fn output(rows: &[String]) -> String {
    let mut text = format!("{NOTE}\n{COLUMNS}\n");
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    text
}

/// A remainder as the three texts a comparison reads it from: each side's lines, each
/// ended by `\n`, and the counts.
#[derive(Debug, PartialEq)]
struct Remainder {
    expected: String,
    actual: String,
    counters: String,
}

fn worker_count() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
}

/// The lines of a list file, an empty line being the empty value.
fn list_lines(path: &str) -> Vec<String> {
    let text = fs::read_to_string(path).unwrap();
    let mut lines: Vec<String> = text.split('\n').map(str::to_string).collect();
    if lines.last().is_some_and(String::is_empty) {
        lines.pop();
    }
    lines
}

/// The plan of the options `more`: `--force-locations FILE`, `--force-text FILE`,
/// `--keep-values COLUMN FILE` and `--unordered`.
fn plan_of(more: &[&str], threads: usize) -> Plan {
    let mut plan = Plan {
        format: Format::Default,
        force_locations: Vec::new(),
        force_text: Vec::new(),
        keep_values: None,
        unordered: false,
        threads,
        input: None,
        force_records: Vec::new(),
        floor_fields: Vec::new(),
    };
    let mut more = more.iter();
    while let Some(option) = more.next() {
        let mut value = || *more.next().unwrap();
        match *option {
            "--force-locations" => plan.force_locations = list_lines(value()),
            "--force-text" => plan.force_text = list_lines(value()),
            "--keep-values" => {
                let column = value().to_string();
                plan.keep_values = Some((column, list_lines(value())));
            }
            "--unordered" => plan.unordered = true,
            other => panic!("{other} is no option of a plan"),
        }
    }
    plan
}

fn json_string(text: &str) -> String {
    let mut out = String::from("\"");
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// The counts of a remainder as one JSON text.
fn counters_json(r: &stream::Remainder) -> String {
    let whole = r.whole.as_deref().map_or("null".to_string(), json_string);
    let populated: Vec<String> = r
        .kept_populated
        .iter()
        .map(|(name, n)| format!("  {}: {n}", json_string(name)))
        .collect();
    let runs: Vec<String> = r
        .kept_before
        .iter()
        .map(|(n, b)| format!("[{n}, {b}]"))
        .collect();
    let side = |lines: u64, handed_over: u64, bytes: u64| {
        format!(
            "{{\"lines\": {lines}, \"handed_over\": {handed_over}, \"bytes_handed_over\": {bytes}}}"
        )
    };
    let block = |items: &[String]| {
        if items.is_empty() {
            String::new()
        } else {
            format!("\n{}\n ", items.join(",\n"))
        }
    };
    format!(
        "{{\n \"whole\": {whole},\n \"kept_entries\": {},\n \"kept_populated\": {{{}}},\n \"kept_before\": [{}],\n \"groups\": {},\n \"groups_handed_over\": {},\n \"expected\": {},\n \"actual\": {}\n}}\n",
        r.kept_entries,
        block(&populated),
        runs.join(", "),
        r.groups,
        r.groups_handed_over,
        side(r.expected_lines, r.expected_handed_over, r.expected_bytes_handed_over),
        side(r.actual_lines, r.actual_handed_over, r.actual_bytes_handed_over),
    )
}

fn file_text(lines: &[String]) -> String {
    lines.iter().fold(String::new(), |mut text, line| {
        text.push_str(line);
        text.push('\n');
        text
    })
}

/// The remainder of one pair, the same at one worker, at two and at the default, and
/// through blocks shorter than a line, of a few lines and of the default size.
fn remainder_of(expected: &Path, actual: &Path, more: &[&str]) -> Remainder {
    let mut runs = Vec::new();
    for (threads, block) in [
        (1, None),
        (2, None),
        (worker_count(), None),
        (2, Some(48)),
        (3, Some(700)),
    ] {
        let plan = plan_of(more, threads);
        let done = match block {
            Some(block) => stream::remainder_in_blocks(expected, actual, &plan, block),
            None => stream::remainder(expected, actual, &plan),
        };
        let done = done.unwrap_or_else(|error| panic!("{error}"));
        assert_eq!(done.restarts, 0);
        assert_eq!(done.candidate_error, None);
        runs.push(done);
    }
    for run in &runs[1..] {
        assert_eq!(run, &runs[0]);
    }
    let run = &runs[0];
    Remainder {
        expected: file_text(&run.expected),
        actual: file_text(&run.actual),
        counters: counters_json(run),
    }
}

fn remainder(expected: &[String], actual: &[String], more: &[&str]) -> Remainder {
    let scratch = Scratch::new();
    let e = scratch.file("expected.txt", output(expected).as_bytes());
    let a = scratch.file("actual.txt", output(actual).as_bytes());
    remainder_of(&e, &a, more)
}

/// The data rows of a remainder text, its header lines checked and dropped.
fn data(text: &str) -> Vec<String> {
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some(NOTE));
    assert_eq!(lines.next(), Some(COLUMNS));
    lines.map(str::to_string).collect()
}

fn counter(counters: &str, name: &str) -> String {
    let key = format!("\"{name}\": ");
    let at = counters
        .find(&key)
        .unwrap_or_else(|| panic!("{name} in {counters}"))
        + key.len();
    counters[at..]
        .lines()
        .next()
        .unwrap()
        .trim_end_matches(',')
        .to_string()
}

/// Ten records at ten Locations, two transcript rows each.
fn ten_records() -> Vec<String> {
    (0..10)
        .flat_map(|r| {
            let location = format!("1:{}", 1000 + r);
            [
                row(r, &location, "tx1", "IMPACT=MODIFIER;STRAND=1"),
                row(r, &location, "tx2", "IMPACT=MODIFIER"),
            ]
        })
        .collect()
}

#[test]
fn outputs_written_alike_leave_their_header_lines_alone() {
    let rows = ten_records();
    let r = remainder(&rows, &rows, &[]);
    assert!(data(&r.expected).is_empty() && data(&r.actual).is_empty());
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "20");
    assert_eq!(counter(&r.counters, "IMPACT"), "20");
    assert_eq!(counter(&r.counters, "STRAND"), "10");
    assert_eq!(counter(&r.counters, "Extra"), "20");
    assert_eq!(counter(&r.counters, "kept_before"), "[]");
    assert_eq!(counter(&r.counters, "groups"), "10");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "0");
}

#[test]
fn a_differing_field_hands_its_group_over_whole() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual[7] = row(3, "1:1003", "tx2", "IMPACT=LOW");
    let r = remainder(&expected, &actual, &[]);
    assert_eq!(data(&r.expected), expected[6..8]);
    assert_eq!(data(&r.actual), actual[6..8]);
    assert_eq!(counter(&r.counters, "kept_entries"), "18");
    assert_eq!(counter(&r.counters, "kept_before"), "[[2, 6]]");
    assert_eq!(counter(&r.counters, "IMPACT"), "18");
    assert_eq!(
        counter(&r.counters, "expected"),
        format!(
            "{{\"lines\": 20, \"handed_over\": 2, \"bytes_handed_over\": {}}}",
            expected[6].len() + expected[7].len() + 2
        )
    );
}

#[test]
fn a_row_one_side_alone_writes_is_handed_over_from_either_side() {
    let both = ten_records();
    let mut longer = both.clone();
    longer.insert(4, row(1, "1:1001", "tx3", "IMPACT=MODIFIER"));
    longer.push(row(10, "1:2000", "tx1", "IMPACT=MODIFIER"));
    let r = remainder(&longer, &both, &[]);
    assert_eq!(data(&r.expected), [&longer[2..5], &longer[21..]].concat());
    assert_eq!(data(&r.actual), both[2..4]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[3, 2], [1, 18]]");
    let r = remainder(&both, &longer, &[]);
    assert_eq!(data(&r.expected), both[2..4]);
    assert_eq!(data(&r.actual), [&longer[2..5], &longer[21..]].concat());
    assert_eq!(counter(&r.counters, "kept_before"), "[[2, 2]]");
}

#[test]
fn a_repeated_row_is_kept_back_only_with_its_repeat_on_both_sides() {
    let mut expected = ten_records();
    expected.insert(1, expected[0].clone());
    let r = remainder(&expected, &expected, &[]);
    assert!(data(&r.expected).is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "21");
    let r = remainder(&expected, &ten_records(), &[]);
    assert_eq!(data(&r.expected), expected[..3]);
    assert_eq!(data(&r.actual), expected[1..3]);
}

#[test]
fn a_listed_location_is_handed_over_with_equal_rows() {
    let rows = ten_records();
    let scratch = Scratch::new();
    let list = scratch.file("locations.txt", b"1:1004\n9:9\n");
    let r = remainder(&rows, &rows, &["--force-locations", list.to_str().unwrap()]);
    assert_eq!(data(&r.expected), rows[8..10]);
    assert_eq!(data(&r.actual), rows[8..10]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[2, 8]]");
}

#[test]
fn a_reference_row_of_a_value_off_the_kept_list_hands_its_group_over() {
    let mut rows = ten_records();
    rows[5] = rows[5].replace("\tTranscript\t", "\tRegulatoryFeature\t");
    let scratch = Scratch::new();
    let list = scratch.file("kept.txt", b"Transcript\n-\n\n");
    let more = ["--keep-values", "Feature_type", list.to_str().unwrap()];
    let r = remainder(&rows, &rows, &more);
    assert_eq!(data(&r.expected), rows[4..6]);
    assert_eq!(data(&r.actual), rows[4..6]);
    let r = remainder(&rows, &rows, &[]);
    assert!(data(&r.expected).is_empty());
}

#[test]
fn a_reference_row_holding_a_listed_text_hands_its_group_over() {
    let mut rows = ten_records();
    rows[12] = row(6, "1:1006", "tx1", "IMPACT=MODIFIER;misspeltKey100way=0.5");
    let scratch = Scratch::new();
    let list = scratch.file("texts.txt", b"misspeltKey\n");
    let r = remainder(&rows, &rows, &["--force-text", list.to_str().unwrap()]);
    assert_eq!(data(&r.expected), rows[12..14]);
    assert_eq!(data(&r.actual), rows[12..14]);
    assert_eq!(counter(&r.counters, "kept_entries"), "18");
}

#[test]
fn a_listed_text_in_the_column_line_or_a_key_column_hands_everything_over() {
    let reason = "\"a listed text stands in the reference's column line or a key column\"";
    let scratch = Scratch::new();
    let list = scratch.file("texts.txt", b"misspeltKey\n");
    let more = ["--force-text", list.to_str().unwrap()];
    let mut rows = ten_records();
    rows[3] = row(1, "1:1001", "misspeltKey7", "IMPACT=MODIFIER");
    let r = remainder(&rows, &rows, &more);
    assert_eq!(counter(&r.counters, "whole"), reason);
    assert_eq!(data(&r.expected), rows);
    // The same text in a column that is no key column hands its group over, no more.
    let mut rows = ten_records();
    rows[3] = rows[3].replace("\tgene1\t", "\tmisspeltKey7\t");
    let r = remainder(&rows, &rows, &more);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(data(&r.expected), rows[2..4]);
    // In the column line of both outputs.
    let text = output(&ten_records()).replace("\tGene\t", "\tmisspeltKey100way\t");
    let e = scratch.file("e.txt", text.as_bytes());
    let a = scratch.file("a.txt", text.as_bytes());
    let r = remainder_of(&e, &a, &more);
    assert_eq!(counter(&r.counters, "whole"), reason);
    assert_eq!(r.expected, text);
    let r = remainder_of(&e, &a, &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "20");
}

#[test]
fn a_location_whose_rows_stand_apart_is_one_group() {
    // Record 2's Location is also the Location of a row after record 4's rows.
    let mut expected = ten_records();
    expected.insert(10, row(2, "1:1002", "tx9", "IMPACT=MODIFIER"));
    let r = remainder(&expected, &expected, &[]);
    assert!(data(&r.expected).is_empty());
    assert_eq!(counter(&r.counters, "groups"), "10");
    // A difference in the late row hands over its group and what stands inside its span:
    // the rows of records 3 and 4.
    let mut actual = expected.clone();
    actual[10] = row(2, "1:1002", "tx9", "IMPACT=LOW");
    let r = remainder(&expected, &actual, &[]);
    assert_eq!(data(&r.expected), expected[4..11]);
    assert_eq!(data(&r.actual), actual[4..11]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[7, 4]]");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "3");
    // Read without order, the group alone is handed over.
    let r = remainder(&expected, &actual, &["--unordered"]);
    assert_eq!(
        data(&r.expected),
        [&expected[4..6], &expected[10..11]].concat()
    );
    assert_eq!(counter(&r.counters, "kept_before"), "[[2, 4], [1, 8]]");
}

#[test]
fn a_group_with_its_rows_in_another_order_is_handed_over() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual.swap(8, 9);
    let r = remainder(&expected, &actual, &[]);
    assert_eq!(data(&r.expected), expected[8..10]);
    assert_eq!(data(&r.actual), actual[8..10]);
}

#[test]
fn groups_in_another_order_are_handed_over_unless_the_order_is_unread() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual.swap(2, 4);
    actual.swap(3, 5);
    let r = remainder(&expected, &actual, &[]);
    assert_eq!(data(&r.expected), expected[2..6]);
    assert_eq!(data(&r.actual), actual[2..6]);
    assert_eq!(counter(&r.counters, "kept_before"), "[[4, 2]]");
    let r = remainder(&expected, &actual, &["--unordered"]);
    assert!(data(&r.expected).is_empty() && data(&r.actual).is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "20");
}

#[test]
fn note_lines_may_differ_and_a_differing_column_line_hands_everything_over() {
    let rows = ten_records();
    let scratch = Scratch::new();
    let e = scratch.file("expected.txt", output(&rows).as_bytes());
    let other = output(&rows).replace(NOTE, "## written by the candidate\n## with a second note");
    let a = scratch.file("actual.txt", other.as_bytes());
    let r = remainder_of(&e, &a, &[]);
    assert!(data(&r.expected).is_empty());
    assert_eq!(
        r.actual,
        format!("## written by the candidate\n## with a second note\n{COLUMNS}\n")
    );
    assert_eq!(counter(&r.counters, "kept_entries"), "20");

    let renamed = output(&rows).replace("\tGene\t", "\tGENE\t");
    let a = scratch.file("renamed.txt", renamed.as_bytes());
    let r = remainder_of(&e, &a, &[]);
    assert_eq!(r.expected, output(&rows));
    assert_eq!(r.actual, renamed);
    assert_eq!(
        counter(&r.counters, "whole"),
        "\"the two outputs' column lines differ\""
    );
    assert_eq!(counter(&r.counters, "kept_entries"), "0");
    assert_eq!(counter(&r.counters, "kept_populated"), "{}");
}

#[test]
fn shapes_a_line_s_text_does_not_settle_hand_everything_over() {
    let rows = ten_records();
    let whole = |text: String, reason: &str| {
        let scratch = Scratch::new();
        let e = scratch.file("expected.txt", text.as_bytes());
        let a = scratch.file("actual.txt", text.as_bytes());
        let r = remainder_of(&e, &a, &[]);
        assert_eq!(r.expected, text);
        assert_eq!(r.actual, text);
        assert_eq!(counter(&r.counters, "whole"), format!("\"{reason}\""));
    };
    whole(
        output(&rows) + COLUMNS + "\n",
        "a header line follows a data line",
    );
    whole(
        format!("{NOTE}\n{}\n", rows[0]),
        "a data line precedes every column line",
    );
    whole(
        output(&rows).replace("\tGene\t", "\tAllele\t"),
        "the column line names a column twice",
    );
    whole(
        output(&rows) + "short\trow\n",
        "a data line's width is not the column line's",
    );
    whole(
        output(&rows) + "\n",
        "a data line's width is not the column line's",
    );
    whole(
        output(&[row(0, "2:7|8", "tx1", "IMPACT=MODIFIER")]),
        "a key column holds | or #",
    );
    whole(
        output(&[row(0, "2:7", "tx1#2", "IMPACT=MODIFIER")]),
        "a key column holds | or #",
    );
    whole(
        output(&[row(0, "2:7 8", "tx1", "IMPACT=MODIFIER")]),
        "a Location holds a space",
    );
    whole(
        output(&[row(0, "2:7", "tx1", "IMPACT=MODIFIER;Feature=tx2")]),
        "an Extra pair is named as a key column",
    );
}

#[test]
fn every_line_end_and_a_gzipped_output_read_as_the_plain_text() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual[7] = row(3, "1:1003", "tx2", "IMPACT=LOW");
    let plain = remainder(&expected, &actual, &[]);

    let scratch = Scratch::new();
    let crlf = output(&expected).replace('\n', "\r\n");
    let e = scratch.file("expected.txt", crlf.trim_end().as_bytes());
    let mut gz = GzEncoder::new(Vec::new(), Compression::default());
    let text = output(&actual);
    let (first, second) = text.split_at(text.len() / 2);
    gz.write_all(first.as_bytes()).unwrap();
    let mut bytes = gz.finish().unwrap();
    let mut gz = GzEncoder::new(Vec::new(), Compression::default());
    gz.write_all(second.replace('\n', "\r").as_bytes()).unwrap();
    bytes.extend(gz.finish().unwrap());
    // Zero bytes after the last member are padding.
    bytes.extend_from_slice(&[0; 64]);
    let a = scratch.file("actual.txt.gz", &bytes);
    let r = remainder_of(&e, &a, &[]);
    assert_eq!(r, plain);
}

/// One gzip member holding `text` under a header of `flags`, `fields` standing between the
/// fixed header and `stream`, the raw deflate stream.
fn gzip_member(flags: u8, fields: &[u8], stream: &[u8], text: &[u8]) -> Vec<u8> {
    let mut crc = flate2::Crc::new();
    crc.update(text);
    let mut member = vec![0x1f, 0x8b, 8, flags, 0, 0, 0, 0, 0, 0xff];
    member.extend_from_slice(fields);
    member.extend_from_slice(stream);
    member.extend_from_slice(&crc.sum().to_le_bytes());
    member.extend_from_slice(&(text.len() as u32).to_le_bytes());
    member
}

fn deflate(text: &[u8]) -> Vec<u8> {
    let mut stream = DeflateEncoder::new(Vec::new(), Compression::default());
    stream.write_all(text).unwrap();
    stream.finish().unwrap()
}

/// One final fixed-code deflate block that writes `before`, copies three bytes from five
/// bytes back, and writes `after`.
fn deflate_with_a_copy(before: &[u8], after: &[u8]) -> Vec<u8> {
    fn code(bits: &mut Vec<bool>, value: u16, width: u8) {
        bits.extend((0..width).rev().map(|i| value >> i & 1 == 1));
    }
    fn literals(bits: &mut Vec<bool>, text: &[u8]) {
        for &byte in text {
            assert!(byte < 144, "a literal of the eight-bit codes");
            code(bits, 0x30 + u16::from(byte), 8);
        }
    }
    // The last block, of the fixed codes.
    let mut bits = vec![true, true, false];
    literals(&mut bits, before);
    // Length 3, then the distance code of 5 and 6 with its extra bit clear.
    code(&mut bits, 1, 7);
    code(&mut bits, 4, 5);
    bits.push(false);
    literals(&mut bits, after);
    code(&mut bits, 0, 7);
    bits.chunks(8)
        .map(|byte| {
            byte.iter()
                .enumerate()
                .fold(0u8, |acc, (k, &bit)| acc | u8::from(bit) << k)
        })
        .collect()
}

#[test]
fn a_gzip_header_is_read_as_the_comparator_reads_it() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual[7] = row(3, "1:1003", "tx2", "IMPACT=LOW");
    let plain = remainder(&expected, &actual, &[]);

    let scratch = Scratch::new();
    let (e_text, a_text) = (output(&expected), output(&actual));
    // Every flag bit the format reserves beside a header checksum that is not the header's.
    let e = scratch.file(
        "expected.txt.gz",
        &gzip_member(
            0xe2,
            &[0xff, 0xff],
            &deflate(e_text.as_bytes()),
            e_text.as_bytes(),
        ),
    );
    for (flags, fields) in [
        (0x02u8, &b"\x12\x34"[..]),
        (0x20, &b""[..]),
        (0x1e, &b"\x02\x00xyname\0note\0\0\0"[..]),
    ] {
        let a = scratch.file(
            "actual.txt.gz",
            &gzip_member(
                flags,
                fields,
                &deflate(a_text.as_bytes()),
                a_text.as_bytes(),
            ),
        );
        let r = remainder_of(&e, &a, &[]);
        assert_eq!(r, plain);
    }
}

#[test]
fn a_row_the_comparator_may_fail_to_convert_hands_its_group_over() {
    let plain = ten_records();
    for extra in [
        // A digit outside ASCII, a run of 641 digits and 65 opening brackets.
        "IMPACT=MODIFIER;SHIFT_LENGTH=\u{b2}".to_string(),
        format!("IMPACT=MODIFIER;SHIFT_LENGTH={}", "7".repeat(641)),
        format!("IMPACT=MODIFIER;OBJECT={}1", "{\"a\":".repeat(65)),
    ] {
        let mut rows = plain.clone();
        rows[9] = row(4, "1:1004", "tx2", &extra);
        let r = remainder(&rows, &rows, &[]);
        assert_eq!(data(&r.expected), rows[8..10]);
        assert_eq!(data(&r.actual), rows[8..10]);
        assert_eq!(counter(&r.counters, "kept_entries"), "18");
        assert_eq!(counter(&r.counters, "kept_before"), "[[2, 8]]");
    }
    // At the two bounds the row is kept back.
    let mut rows = plain.clone();
    let extra = format!(
        "IMPACT=MODIFIER;SHIFT_LENGTH={};OBJECT={}",
        "7".repeat(640),
        "[".repeat(64)
    );
    rows[9] = row(4, "1:1004", "tx2", &extra);
    let r = remainder(&rows, &rows, &[]);
    assert!(data(&r.expected).is_empty() && data(&r.actual).is_empty());
    assert_eq!(counter(&r.counters, "kept_entries"), "20");
}

/// 3,000 records of up to four rows; every 97th record differs in one row, every 211th
/// has a row on one side alone, and every 389th record's Location returns 40 rows later.
/// With the two sides, the number of records that differ.
fn many_rows() -> (Vec<String>, Vec<String>, usize) {
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    let mut late: Vec<(usize, String)> = Vec::new();
    let mut differing_groups = 0;
    for r in 0..3000usize {
        let location = format!("{}:{}", 1 + r / 1000, 5000 + 3 * r);
        while late.first().is_some_and(|(at, _)| *at <= r) {
            let (_, row) = late.remove(0);
            expected.push(row.clone());
            actual.push(row);
        }
        for t in 0..1 + r % 4 {
            let feature = format!("tx{t}");
            let row_e = row(r, &location, &feature, "IMPACT=MODIFIER;DISTANCE=7");
            let row_a = if r % 97 == 5 && t == 0 {
                row(r, &location, &feature, "IMPACT=MODIFIER;DISTANCE=8")
            } else {
                row_e.clone()
            };
            expected.push(row_e);
            if !(r % 211 == 9 && t == 0) {
                actual.push(row_a);
            }
        }
        if r % 97 == 5 || r % 211 == 9 {
            differing_groups += 1;
        }
        if r % 389 == 11 {
            late.push((r + 40, row(r, &location, "late", "IMPACT=MODIFIER")));
        }
    }
    (expected, actual, differing_groups)
}

#[test]
fn many_rows_split_alike_at_every_worker_count() {
    let (expected, actual, differing_groups) = many_rows();
    let scratch = Scratch::new();
    let e = scratch.file("expected.txt", output(&expected).as_bytes());
    let a = scratch.file("actual.txt", output(&actual).as_bytes());
    let three = remainder_of(&e, &a, &[]);
    let seven = counters_json(&stream::remainder(&e, &a, &plan_of(&[], 7)).unwrap());
    assert_eq!(seven, three.counters);
    assert_eq!(counter(&three.counters, "whole"), "null");
    assert_eq!(counter(&three.counters, "groups"), "3000");
    let handed: usize = counter(&three.counters, "groups_handed_over")
        .parse()
        .unwrap();
    assert!(
        handed >= differing_groups && handed < 3000 / 4,
        "{handed} groups handed over"
    );
    let kept: usize = counter(&three.counters, "kept_entries").parse().unwrap();
    assert_eq!(kept + data(&three.expected).len(), expected.len());
    assert_eq!(kept + data(&three.actual).len(), actual.len());
    // A kept-back row is on neither remainder, so every row that differs is on both.
    let handed_rows = data(&three.expected);
    for (r, row_e) in expected
        .iter()
        .enumerate()
        .filter(|(_, row)| row.contains("DISTANCE=7"))
    {
        let record: usize = row_e[3..row_e.find('\t').unwrap()].parse().unwrap();
        if record % 97 == 5 || record % 211 == 9 {
            assert!(
                handed_rows.contains(row_e),
                "row {r} of a differing record was kept back"
            );
        }
    }
}

#[test]
fn an_output_named_as_the_help_option_is_read() {
    let rows = ten_records();
    let scratch = Scratch::new();
    let e = scratch.file("-h", output(&rows).as_bytes());
    let a = scratch.file("--help", output(&rows).as_bytes());
    let r = remainder_of(&e, &a, &[]);
    assert_eq!(counter(&r.counters, "kept_entries"), "20");
}

fn refused(expected: &Path, actual: &Path) -> String {
    let mut texts = Vec::new();
    for (threads, block) in [(1, None), (2, None), (3, Some(40))] {
        let plan = plan_of(&[], threads);
        let done = match block {
            Some(block) => stream::remainder_in_blocks(expected, actual, &plan, block),
            None => stream::remainder(expected, actual, &plan),
        };
        match done {
            Err(Error::Comparator(text)) => texts.push(text),
            Ok(_) => panic!("{} and {} were read", expected.display(), actual.display()),
        }
    }
    assert!(texts.iter().all(|text| *text == texts[0]), "{texts:?}");
    texts.remove(0)
}

#[test]
fn an_input_that_cannot_be_read_is_refused_in_the_comparator_s_words() {
    let scratch = Scratch::new();
    let text = output(&ten_records());
    let good = scratch.file("good.txt", text.as_bytes());
    let mut bytes = text.clone().into_bytes();
    bytes.extend_from_slice(b"rec\t2:7\t\xff\n");
    let bad = scratch.file("bad.txt", &bytes);
    // A deflate stream that copies from before its start, which zlib refuses, under the
    // trailer of the text a decoder taking the copy writes: zero bytes for the two places
    // before the start, then the first byte.
    let whole = text.clone().into_bytes();
    let (before, after) = whole.split_at(3);
    let taken = [before, &[0, 0, before[0]], after].concat();
    let copied = scratch.file(
        "copied.txt.gz",
        &gzip_member(0, b"", &deflate_with_a_copy(before, after), &taken),
    );
    let absent = scratch.0.join("absent.txt");
    let not_utf8 = format!(
        "'utf-8' codec can't decode byte 0xff in position {}: invalid start byte",
        text.len() + 8
    );
    for (expected, actual, unread, detail) in [
        (&good, &bad, &bad, not_utf8.as_str()),
        (&bad, &good, &bad, not_utf8.as_str()),
        (&bad, &absent, &bad, not_utf8.as_str()),
        (&good, &absent, &absent, "No such file or directory"),
        (&absent, &bad, &absent, "No such file or directory"),
        (
            &good,
            &copied,
            &copied,
            "Error -3 while decompressing data: invalid distance too far back",
        ),
        (&good, &scratch.0, &scratch.0, "Is a directory"),
    ] {
        assert_eq!(
            refused(expected, actual),
            format!("cannot read {}: {detail}", unread.display())
        );
    }
}

#[test]
fn a_gzip_file_gzip_open_refuses_is_refused_in_its_words() {
    let scratch = Scratch::new();
    let text = output(&ten_records());
    let good = scratch.file("good.txt", text.as_bytes());
    let member = gzip_member(0, b"", &deflate(text.as_bytes()), text.as_bytes());
    let last = member.len() - 1;
    let mut wrong_length = member.clone();
    wrong_length[last] ^= 1;
    let mut wrong_checksum = member.clone();
    wrong_checksum[last - 4] ^= 1;
    let mut crc = flate2::Crc::new();
    crc.update(text.as_bytes());
    let checksum = format!(
        "CRC check failed {:#x} != {:#x}",
        crc.sum() ^ 1 << 24,
        crc.sum()
    );
    let mut other_method = member.clone();
    other_method[2] = 7;
    let cut_short = "Compressed file ended before the end-of-stream marker was reached";
    for (bytes, detail) in [
        (wrong_length, "Incorrect length of data produced"),
        (wrong_checksum, checksum.as_str()),
        (other_method, "Unknown compression method"),
        (member[..last].to_vec(), cut_short),
        (member[..40].to_vec(), cut_short),
        (member[..5].to_vec(), cut_short),
        (
            [&member[..], b"more"].concat(),
            "Not a gzipped file (b'mo')",
        ),
        (
            [&member[..], &[0, 0, 0, 0x1f]].concat(),
            "Not a gzipped file (b'\\x1f')",
        ),
        (b"plain text".to_vec(), "Not a gzipped file (b'pl')"),
    ] {
        let bad = scratch.file("bad.txt.gz", &bytes);
        let expected = format!("cannot read {}: {detail}", bad.display());
        assert_eq!(refused(&good, &bad), expected);
        assert_eq!(refused(&bad, &good), expected);
    }
    // A gzip file that fails after text that is not UTF-8 fails as a gzip file.
    let mut bytes = text.clone().into_bytes();
    bytes.extend_from_slice(b"rec\t2:7\t\xff\n");
    let mut member = gzip_member(0, b"", &deflate(&bytes), &bytes);
    let last = member.len() - 1;
    member[last] ^= 1;
    let bad = scratch.file("both.txt.gz", &member);
    assert_eq!(
        refused(&bad, &good),
        format!(
            "cannot read {}: Incorrect length of data produced",
            bad.display()
        )
    );
    // An empty gzip file is an empty text.
    let empty = scratch.file("empty.txt.gz", b"");
    let r = stream::remainder(&empty, &good, &plan_of(&[], 2)).unwrap();
    assert!(r.expected.is_empty());
    assert_eq!(r.actual.len(), 22);
    assert_eq!(
        r.whole.as_deref(),
        Some("the two outputs' column lines differ")
    );
}

#[test]
fn a_candidate_that_cannot_be_read_comes_after_a_reference_that_may_not_parse() {
    let scratch = Scratch::new();
    let absent = scratch.0.join("absent.txt");
    let unread = format!(
        "cannot read {}: No such file or directory",
        absent.display()
    );
    // A reference of rows the comparison reads without an error of its own: the
    // candidate's error ends the call.
    let good = scratch.file("good.txt", output(&ten_records()).as_bytes());
    assert_eq!(refused(&good, &absent), unread);
    // A reference row of another width, a row past a limit of the comparison's
    // interpreter and a reference without a column line may end the comparison first: the
    // whole reference comes back beside the candidate's error.
    let mut wide = ten_records();
    wide.push("short\trow".to_string());
    let mut digits = ten_records();
    digits[4] = row(2, "1:1002", "tx1", &format!("N={}", "7".repeat(641)));
    for text in [
        output(&wide),
        output(&digits),
        format!("{NOTE}\n{}\n", ten_records()[0]),
    ] {
        let fragile = scratch.file("fragile.txt", text.as_bytes());
        for (threads, block) in [(1, None), (2, Some(64))] {
            let plan = plan_of(&[], threads);
            let r = match block {
                Some(block) => stream::remainder_in_blocks(&fragile, &absent, &plan, block),
                None => stream::remainder(&fragile, &absent, &plan),
            }
            .unwrap();
            assert_eq!(r.candidate_error, Some(Error::Comparator(unread.clone())));
            assert_eq!(file_text(&r.expected), text);
            assert!(r.actual.is_empty());
            assert_eq!(r.kept_entries, 0);
            assert_eq!(r.whole, None);
        }
    }
}

#[test]
fn a_pipe_is_read_once_and_gives_the_remainder_of_the_file() {
    let expected = ten_records();
    let mut actual = expected.clone();
    actual[7] = row(3, "1:1003", "tx2", "IMPACT=LOW");
    let plain = remainder(&expected, &actual, &[]);

    let scratch = Scratch::new();
    let e = scratch.file("expected.txt", output(&expected).as_bytes());
    let pipe = scratch.0.join("actual.pipe");
    let made = std::process::Command::new("mkfifo")
        .arg(&pipe)
        .status()
        .unwrap();
    assert!(made.success());
    let text = output(&actual);
    let r = std::thread::scope(|scope| {
        scope.spawn(|| fs::write(&pipe, text.as_bytes()).unwrap());
        stream::remainder_in_blocks(&e, &pipe, &plan_of(&[], 2), 64).unwrap()
    });
    assert_eq!(file_text(&r.expected), plain.expected);
    assert_eq!(file_text(&r.actual), plain.actual);
    assert_eq!(counters_json(&r), plain.counters);
}

#[test]
fn a_location_that_returns_blocks_later_is_one_group() {
    // Locations of two rows; every 50th returns 200 records later with a third row, which
    // every other time differs on the candidate side.
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    let mut late: Vec<(usize, String, String)> = Vec::new();
    for r in 0..=2200usize {
        while late.first().is_some_and(|(at, _, _)| *at <= r) {
            let (_, row_e, row_a) = late.remove(0);
            expected.push(row_e);
            actual.push(row_a);
        }
        if r >= 2000 {
            continue;
        }
        let location = format!("3:{}", 100 + r);
        for feature in ["tx1", "tx2"] {
            expected.push(row(r, &location, feature, "IMPACT=MODIFIER"));
            actual.push(row(r, &location, feature, "IMPACT=MODIFIER"));
        }
        if r % 50 == 7 {
            let row_e = row(r, &location, "late", "IMPACT=MODIFIER");
            let row_a = if r % 100 == 7 {
                row(r, &location, "late", "IMPACT=LOW")
            } else {
                row_e.clone()
            };
            late.push((r + 200, row_e, row_a));
        }
    }
    assert_eq!(expected.len(), 4040);
    // Read without order, the 20 groups that differ are handed over, three rows each.
    let r = remainder(&expected, &actual, &["--unordered"]);
    assert_eq!(counter(&r.counters, "groups"), "2000");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "20");
    assert_eq!(counter(&r.counters, "kept_entries"), "3980");
    let handed_rows = data(&r.expected);
    assert_eq!(handed_rows.len(), 60);
    for row_e in expected.iter().filter(|row| row.contains("\tlate\t")) {
        let record: usize = row_e[3..row_e.find('\t').unwrap()].parse().unwrap();
        assert_eq!(handed_rows.contains(row_e), record % 100 == 7);
    }
    // Read with order, a group that differs also hands over the records inside its span,
    // and the spans of the 20 reach from the eighth record to the last.
    let r = remainder(&expected, &actual, &[]);
    assert_eq!(counter(&r.counters, "groups"), "2000");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "1993");
    assert_eq!(counter(&r.counters, "kept_entries"), "14");
    assert_eq!(counter(&r.counters, "kept_before"), "[[4026, 14]]");
}

#[test]
fn a_list_holds_no_value_with_a_tab_or_a_line_end() {
    assert_eq!(
        stream::listed(["1:5", "", "a\tb", "a\nb", "a\rb", "5\u{2032}UTR"]),
        ["1:5", "", "5\u{2032}UTR"]
    );
}

/// Records 0 to `last`, two rows each, the rows of record k standing around the first row
/// of record k + 1, and one more record after them; the first row of record 0 differs on
/// the candidate side.
fn chained(last: usize) -> (Vec<String>, Vec<String>) {
    let mut order = vec![0, 1, 0];
    for record in 2..=last {
        order.extend([record, record - 1]);
    }
    order.push(last);
    let mut seen = vec![0; last + 1];
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    for record in order {
        seen[record] += 1;
        let location = format!("5:{}", 100 + record);
        let feature = format!("tx{}", seen[record]);
        expected.push(row(record, &location, &feature, "IMPACT=MODIFIER"));
        let impact = if (record, seen[record]) == (0, 1) {
            "IMPACT=LOW"
        } else {
            "IMPACT=MODIFIER"
        };
        actual.push(row(record, &location, &feature, impact));
    }
    for feature in ["tx1", "tx2"] {
        expected.push(row(900, "5:900", feature, "IMPACT=MODIFIER"));
        actual.push(row(900, "5:900", feature, "IMPACT=MODIFIER"));
    }
    (expected, actual)
}

#[test]
fn spans_that_open_one_another_past_the_rounds_hand_everything_over() {
    // Each round of the second order condition hands over the next record of the chain:
    // 63 rounds settle the order, and the record after the chain is kept back.
    let (expected, actual) = chained(63);
    let r = remainder(&expected, &actual, &[]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "kept_entries"), "2");
    assert_eq!(counter(&r.counters, "kept_before"), "[[128, 0]]");
    assert_eq!(counter(&r.counters, "groups"), "65");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "64");
    // One record more and the rounds end with the order unsettled.
    let (expected, actual) = chained(64);
    let r = remainder(&expected, &actual, &[]);
    assert_eq!(
        counter(&r.counters, "whole"),
        "\"the kept-back lines found no common order\""
    );
    assert_eq!(data(&r.expected), expected);
    assert_eq!(data(&r.actual), actual);
    assert_eq!(counter(&r.counters, "kept_entries"), "0");
    assert_eq!(counter(&r.counters, "kept_before"), "[[132, 0]]");
    assert_eq!(counter(&r.counters, "groups"), "66");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "66");
    // Read without order, the record that differs alone is handed over.
    let r = remainder(&expected, &actual, &["--unordered"]);
    assert_eq!(counter(&r.counters, "whole"), "null");
    assert_eq!(counter(&r.counters, "groups_handed_over"), "1");
    assert_eq!(counter(&r.counters, "kept_entries"), "130");
}
