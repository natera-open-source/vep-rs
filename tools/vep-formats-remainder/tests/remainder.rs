// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over planted pairs of outputs, every row written here, at one worker,
//! two and the host's default, and holds the three runs' files to the same bytes and to
//! the remainder each shape must leave.

use flate2::write::{DeflateEncoder, GzEncoder};
use flate2::Compression;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};

const BIN: &str = env!("CARGO_BIN_EXE_vep-formats-remainder");
const COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tExtra";
const NOTE: &str = "## written by the reference engine";

static SCRATCH: AtomicU32 = AtomicU32::new(0);

struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-remainder-test-{}-{}",
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

struct Remainder {
    expected: String,
    actual: String,
    counters: String,
}

fn invoke(
    dir: &Path,
    expected: &Path,
    actual: &Path,
    more: &[&str],
    threads: Option<&str>,
) -> Output {
    let mut command = Command::new(BIN);
    command
        .args(["--expected", expected.to_str().unwrap()])
        .args(["--actual", actual.to_str().unwrap()])
        .args(["--out-dir", dir.to_str().unwrap()])
        .args(more);
    match threads {
        Some(n) => command.env("VEP_FORMATS_REMAINDER_THREADS", n),
        None => command.env_remove("VEP_FORMATS_REMAINDER_THREADS"),
    };
    command.output().unwrap()
}

/// The remainder of one pair, the same bytes at one worker, at two and at the default.
fn remainder_of(scratch: &Scratch, expected: &Path, actual: &Path, more: &[&str]) -> Remainder {
    let mut runs = Vec::new();
    for (tag, threads) in [("one", Some("1")), ("two", Some("2")), ("default", None)] {
        let dir = scratch.0.join(tag);
        let done = invoke(&dir, expected, actual, more, threads);
        assert!(
            done.status.success(),
            "{}",
            String::from_utf8_lossy(&done.stderr)
        );
        assert!(done.stdout.is_empty() && done.stderr.is_empty());
        runs.push(Remainder {
            expected: fs::read_to_string(dir.join("expected.txt")).unwrap(),
            actual: fs::read_to_string(dir.join("actual.txt")).unwrap(),
            counters: fs::read_to_string(dir.join("counters.json")).unwrap(),
        });
    }
    for run in &runs[1..] {
        assert_eq!(run.expected, runs[0].expected);
        assert_eq!(run.actual, runs[0].actual);
        assert_eq!(run.counters, runs[0].counters);
    }
    runs.remove(0)
}

fn remainder(expected: &[String], actual: &[String], more: &[&str]) -> Remainder {
    let scratch = Scratch::new();
    let e = scratch.file("expected.txt", output(expected).as_bytes());
    let a = scratch.file("actual.txt", output(actual).as_bytes());
    remainder_of(&scratch, &e, &a, more)
}

/// The data rows of a remainder file, its header lines checked and dropped.
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
    let r = remainder_of(&scratch, &e, &a, &more);
    assert_eq!(counter(&r.counters, "whole"), reason);
    assert_eq!(r.expected, text);
    let r = remainder_of(&scratch, &e, &a, &[]);
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
    let r = remainder_of(&scratch, &e, &a, &[]);
    assert!(data(&r.expected).is_empty());
    assert_eq!(
        r.actual,
        format!("## written by the candidate\n## with a second note\n{COLUMNS}\n")
    );
    assert_eq!(counter(&r.counters, "kept_entries"), "20");

    let renamed = output(&rows).replace("\tGene\t", "\tGENE\t");
    let a = scratch.file("renamed.txt", renamed.as_bytes());
    let r = remainder_of(&scratch, &e, &a, &[]);
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
        let r = remainder_of(&scratch, &e, &a, &[]);
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
    let r = remainder_of(&scratch, &e, &a, &[]);
    assert_eq!(r.expected, plain.expected);
    assert_eq!(r.actual, plain.actual);
    assert_eq!(r.counters, plain.counters);
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
        let r = remainder_of(&scratch, &e, &a, &[]);
        assert_eq!(r.expected, plain.expected);
        assert_eq!(r.actual, plain.actual);
        assert_eq!(r.counters, plain.counters);
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

#[test]
fn many_rows_split_alike_at_every_worker_count() {
    // 3,000 records of up to four rows; every 97th record differs in one row, every 211th
    // has a row on one side alone, and every 389th record's Location returns 40 rows later.
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
    let scratch = Scratch::new();
    let e = scratch.file("expected.txt", output(&expected).as_bytes());
    let a = scratch.file("actual.txt", output(&actual).as_bytes());
    let three = remainder_of(&scratch, &e, &a, &[]);
    let seven = {
        let dir = scratch.0.join("seven");
        assert!(invoke(&dir, &e, &a, &[], Some("7")).status.success());
        fs::read_to_string(dir.join("counters.json")).unwrap()
    };
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
    scratch.file("-h", output(&rows).as_bytes());
    scratch.file("--help", output(&rows).as_bytes());
    let done = Command::new(BIN)
        .current_dir(&scratch.0)
        .args(["--expected", "-h", "--actual", "--help", "--out-dir", "out"])
        .env("VEP_FORMATS_REMAINDER_THREADS", "2")
        .output()
        .unwrap();
    assert!(done.status.success() && done.stdout.is_empty() && done.stderr.is_empty());
    let counters = fs::read_to_string(scratch.0.join("out").join("counters.json")).unwrap();
    assert_eq!(counter(&counters, "kept_entries"), "20");
}

#[test]
fn an_input_this_tool_does_not_read_is_refused_in_one_line() {
    let scratch = Scratch::new();
    let good = scratch.file("good.txt", output(&ten_records()).as_bytes());
    let mut bytes = output(&ten_records()).into_bytes();
    bytes.extend_from_slice(b"rec\t2:7\t\xff\n");
    let bad = scratch.file("bad.txt", &bytes);
    // A deflate stream that copies from before its start, which zlib refuses, under the
    // trailer of the text a decoder taking the copy writes: zero bytes for the two places
    // before the start, then the first byte.
    let text = output(&ten_records()).into_bytes();
    let (before, after) = text.split_at(3);
    let taken = [before, &[0, 0, before[0]], after].concat();
    let copied = scratch.file(
        "copied.txt.gz",
        &gzip_member(0, b"", &deflate_with_a_copy(before, after), &taken),
    );
    for (expected, actual, detail) in [
        (&good, &bad, "the text is not UTF-8"),
        (&bad, &good, "the text is not UTF-8"),
        (&good, &scratch.0.join("absent.txt"), "cannot read"),
        (&good, &copied, "cannot read"),
    ] {
        let done = invoke(&scratch.0.join("out"), expected, actual, &[], Some("2"));
        assert_eq!(done.status.code(), Some(1));
        let stderr = String::from_utf8(done.stderr).unwrap();
        assert_eq!(stderr.lines().count(), 1, "{stderr}");
        assert!(
            stderr.starts_with("ERROR: [vep-formats-remainder] ") && stderr.contains(detail),
            "{stderr}"
        );
    }
    let done = invoke(
        &scratch.0.join("out"),
        &good,
        &good,
        &["--threads", "2"],
        None,
    );
    assert_eq!(done.status.code(), Some(2));
    let done = invoke(&scratch.0.join("out"), &good, &good, &[], Some("none"));
    assert_eq!(done.status.code(), Some(1));
}
