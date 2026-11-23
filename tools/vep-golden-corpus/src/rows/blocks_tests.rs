// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The blocks of a file hold the lines `TextLines` gives, whatever the line ends and the
//! count of threads, and a file that cannot be read is refused as `iter_rows` refuses it.

use super::*;
use crate::rows::{iter_rows, substitution_alts, Row, TextLines};
use crate::TempDir;
use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-blocks-").unwrap()
}

fn file(directory: &TempDir, name: &str, bytes: &[u8]) -> String {
    let path = directory.path.join(name);
    fs::write(&path, bytes).unwrap();
    path.to_str().unwrap().to_string()
}

/// The lines of a file as the blocks give them.
fn block_lines(path: &str, threads: usize) -> Vec<String> {
    let mut all = Vec::new();
    let read = |block: &str| lines(block).map(str::to_string).collect::<Vec<String>>();
    let take = |_: &str, found: Vec<String>| {
        all.extend(found);
        Ok(())
    };
    read_blocks(path, threads, read, take).unwrap();
    all
}

/// The lines of a file as `TextLines` gives them.
fn text_lines(path: &str) -> Vec<String> {
    let mut reader = TextLines::open(path).unwrap();
    let mut all = Vec::new();
    while let Some(line) = reader.next_line().unwrap() {
        all.push(line.to_string());
    }
    all
}

#[test]
fn the_lines_of_a_block_end_where_universal_newlines_end_them() {
    let split = |text: &'static str| lines(text).collect::<Vec<&str>>();
    assert_eq!(split(""), Vec::<&str>::new());
    assert_eq!(split("a"), ["a"]);
    assert_eq!(split("a\n"), ["a"]);
    assert_eq!(split("a\n\nb"), ["a", "", "b"]);
    assert_eq!(split("a\rb\r\nc\n\rd\r"), ["a", "b", "c", "", "d"]);
    assert_eq!(split("\n"), [""]);
    assert_eq!(split("\r\n\r\n"), ["", ""]);
}

#[test]
fn a_block_is_cut_behind_a_line_end_it_can_tell() {
    assert_eq!(line_end_before(b"ab\ncd\nef", 0), Some(6));
    assert_eq!(line_end_before(b"ab\ncd\nef", 6), Some(6));
    assert_eq!(line_end_before(b"ab\ncd\nef", 7), None);
    // A `\r` is a line end once the byte behind it is read and is no `\n`.
    assert_eq!(line_end_before(b"ab\rcd\ref", 0), Some(6));
    assert_eq!(line_end_before(b"ab\rcd\r", 0), Some(3));
    assert_eq!(line_end_before(b"abcd\r", 0), None);
    assert_eq!(line_end_before(b"abcd\rx", 5), Some(5));
    assert_eq!(line_end_before(b"ab\r\ncd", 0), Some(4));
    assert_eq!(line_end_before(b"abcd", 0), None);
}

#[test]
fn blocks_hold_the_lines_of_the_file_on_every_count_of_threads() {
    let scratch = scratch();
    // Lines of every length up to past a block, with each line end and bytes outside UTF-8.
    let mut state = 99u32;
    let mut next = || {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12345);
        (state >> 16) as usize
    };
    let mut body: Vec<u8> = Vec::new();
    for line in 0..40_000 {
        let length = match next() % 50 {
            0 => 0,
            1 => 5000,
            _ => next() % 300,
        };
        for place in 0..length {
            body.push(b'a' + ((line + place) % 26) as u8);
        }
        match next() % 40 {
            0 => body.extend(b"\xff\xfe"),
            1 => body.extend("\u{e9}\u{4e2d}".as_bytes()),
            2 => body.extend(b"\xe2\x82"),
            _ => {}
        }
        body.extend(match next() % 6 {
            0 => &b"\r\n"[..],
            1 => &b"\r"[..],
            _ => &b"\n"[..],
        });
    }
    let long = vec![b'x'; 3 * BLOCK + 17];
    let cases: [(&str, Vec<u8>); 7] = [
        ("mixed.txt", body.clone()),
        ("no-last-end.txt", [&body[..], b"last line"].concat()),
        ("ends-in-return.txt", [&body[..], b"last line\r"].concat()),
        (
            "long-line.txt",
            [&long[..], b"\n", &long[..], b"\r", &long[..]].concat(),
        ),
        (
            "returns-only.txt",
            body.iter()
                .map(|byte| if *byte == b'\n' { b'\r' } else { *byte })
                .collect(),
        ),
        ("empty.txt", Vec::new()),
        ("one-end.txt", b"\n".to_vec()),
    ];
    for (name, bytes) in cases {
        let path = file(&scratch, name, &bytes);
        let expected = text_lines(&path);
        for threads in [1, 2, 7] {
            assert!(
                block_lines(&path, threads) == expected,
                "{name} on {threads} threads"
            );
        }
    }
    assert!(body.len() > 3 * BLOCK);
}

#[test]
fn rows_are_the_rows_iter_rows_gives() {
    let scratch = scratch();
    let mut text = "## header\n#Uploaded_variation\tLocation\n".to_string();
    for place in 0..2000 {
        let position = 1000 + place * 10;
        text += &format!(
            "rs{place}\t21:{position}-{}\tTG\tENSG1\tENST1\tTranscript\tintron_variant,missense_variant\t1\t2\t3\t-\t-\t-\tIMPACT=LOW;STRAND=-1;FLAGS=cds_end_NF\n",
            position + 1
        );
        text += "short\tline\n\n";
    }
    let reference = file(&scratch, "ref.txt", text.as_bytes());
    let input = file(
        &scratch,
        "in.vcf",
        b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t1000\trs0\tAC\tTG\t.\t.\t.\n",
    );
    let substitutions = substitution_alts(&input).unwrap();
    let specs = ["IMPACT=", "HGVSc"];
    let expected: Vec<Row> = iter_rows(&reference, "s", &specs, Some(&substitutions))
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(expected.len(), 2000);
    assert!(expected[0].substitution && !expected[1].substitution);
    let reader = RowReader::new("s", &specs, Some(&substitutions));
    let mut found: Vec<Row> = Vec::new();
    let read = |block: &str| {
        let rows = lines(block).filter_map(|line| reader.view(line));
        rows.map(|view| view.to_row()).collect::<Vec<Row>>()
    };
    let take = |_: &str, rows: Vec<Row>| {
        found.extend(rows);
        Ok(())
    };
    read_blocks(&reference, 3, read, take).unwrap();
    assert!(found == expected);
    assert!(reader.view("#a\tb\tc\td\te\tf\tg").is_none());
    assert!(reader.view("a\tb\tc\td\te\tf").is_none());
    assert_eq!(
        reader.view("a\t1:2\tc\td\te\tf\tg").unwrap().positions,
        ["-"; 3]
    );
}

#[test]
fn a_file_that_cannot_be_read_is_refused_once_the_blocks_before_are_taken() {
    let scratch = scratch();
    let absent = scratch.path.join("absent.txt");
    let absent = absent.to_str().unwrap();
    for threads in [1, 4] {
        let refusal = read_blocks(absent, threads, |_| (), |_, ()| Ok(())).unwrap_err();
        assert_eq!(
            refusal.line,
            format!(
                "ERROR: [build_golden_corpus] {absent} is not a readable text file ([Errno 2] No such file or directory: '{absent}')"
            )
        );
    }
    // A gzip file cut short: the lines of the blocks before the cut are taken, then the
    // refusal comes.
    let mut gzip = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    let line = [&[b'q'; 99][..], b"\n"].concat();
    let stored = line.repeat(600);
    for _ in 0..80 {
        gzip.push(0);
        gzip.extend((stored.len() as u16).to_le_bytes());
        gzip.extend((!(stored.len() as u16)).to_le_bytes());
        gzip.extend(&stored);
    }
    let cut = file(&scratch, "cut.txt.gz", &gzip);
    for threads in [1, 4] {
        let taken = AtomicUsize::new(0);
        let read = |block: &str| lines(block).count();
        let take = |_: &str, count: usize| {
            taken.fetch_add(count, Ordering::Relaxed);
            Ok(())
        };
        let refusal = read_blocks(&cut, threads, read, take).unwrap_err();
        assert_eq!(
            refusal.line,
            format!(
                "ERROR: [build_golden_corpus] {cut} is not a readable text file (Compressed file ended before the end-of-stream marker was reached)"
            )
        );
        let taken = taken.load(Ordering::Relaxed);
        assert!(taken > 0 && taken <= 48_000, "{taken}");
    }
    // A failure of the pass ends the reading, whatever is read ahead.
    let many = file(&scratch, "many.txt", &line.repeat(200_000));
    for threads in [1, 3] {
        let mut blocks = 0;
        let take = |_: &str, ()| {
            blocks += 1;
            if blocks == 2 {
                Err(Failure::own("the second block"))
            } else {
                Ok(())
            }
        };
        let refusal = read_blocks(&many, threads, |_| (), take).unwrap_err();
        assert_eq!(refusal, Failure::own("the second block"));
        assert_eq!(blocks, 2);
    }
}

#[test]
fn the_thread_bound_is_a_whole_number_of_one_or_more() {
    let read = |value: &str| thread_bound(Some(OsStr::new(value)));
    assert_eq!(read("3"), Ok(3));
    assert_eq!(read(" 12 "), Ok(12));
    assert!(read("").is_ok_and(|threads| threads >= 1));
    assert!(thread_bound(None).is_ok_and(|threads| threads >= 1));
    for value in ["0", "-1", "two", "1.5"] {
        assert_eq!(
            read(value),
            Err(Failure::own(format!(
                "VEP_GOLDEN_CORPUS_THREADS={value}: expected a whole number of threads, 1 or more"
            ))),
        );
    }
}
