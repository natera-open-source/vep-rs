// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of `rows::text`: the lines of plain and gzip files, and the reason a file
//! that cannot be read gives, each as CPython's `open` and `gzip.open` give them.

use super::*;
use crate::TempDir;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write as _;

const TWO_LINES: &[u8] = b"a\tb\nc\td\n";

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-text-").unwrap()
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(bytes).unwrap();
    encoder.finish().unwrap()
}

/// Writes `bytes` to `<directory>/<name>` and reads the file's lines: the lines given, and
/// the error that ended the reading.
fn read(directory: &TempDir, name: &str, bytes: &[u8]) -> (Vec<String>, Option<ReadError>) {
    let path = directory.path.join(name);
    std::fs::write(&path, bytes).unwrap();
    lines_of(path.to_str().unwrap())
}

fn lines_of(path: &str) -> (Vec<String>, Option<ReadError>) {
    let mut lines = Vec::new();
    let mut reader = match TextLines::open(path) {
        Ok(reader) => reader,
        Err(error) => return (lines, Some(error)),
    };
    loop {
        match reader.next_line() {
            Ok(Some(line)) => lines.push(line.to_string()),
            Ok(None) => return (lines, None),
            Err(error) => return (lines, Some(error)),
        }
    }
}

fn cpython(reason: &str) -> Option<ReadError> {
    Some(ReadError::Cpython(reason.to_string()))
}

#[test]
fn lines_end_at_universal_newlines_and_bytes_outside_utf8_are_replaced() {
    let directory = scratch();
    assert_eq!(
        read(&directory, "cr.txt", b"a\rb\r\nc\n\rd"),
        (["a", "b", "c", "", "d"].map(str::to_string).to_vec(), None)
    );
    assert_eq!(
        read(
            &directory,
            "bad_utf8.txt",
            b"a\xffb\xe2\x82\n\xe2\x82\xacok\xf0\x9f\n"
        ),
        (
            ["a\u{fffd}b\u{fffd}", "\u{20ac}ok\u{fffd}"]
                .map(str::to_string)
                .to_vec(),
            None
        )
    );
    assert_eq!(
        read(&directory, "bom.txt", b"\xef\xbb\xbf#x\n"),
        (vec!["\u{feff}#x".to_string()], None)
    );
    assert_eq!(read(&directory, "empty.txt", b""), (Vec::new(), None));
    assert_eq!(
        read(&directory, "blank.txt", b"\n\r\n\r"),
        (vec![String::new(); 3], None)
    );
    assert_eq!(
        read(
            &directory,
            "nul.txt",
            b"a\0b\x0b\x0cc\xc2\x85d\xe2\x80\xa8e"
        ),
        (vec!["a\0b\u{0b}\u{0c}c\u{85}d\u{2028}e".to_string()], None)
    );
}

#[test]
fn a_line_longer_than_the_bytes_read_at_a_time_is_one_line() {
    let directory = scratch();
    let long = "x".repeat(3 * CHUNK + 17);
    let text = format!("{long}\r\nshort\n{long}");
    let (lines, error) = read(&directory, "long.txt", text.as_bytes());
    assert_eq!(error, None);
    assert_eq!(lines, [long.as_str(), "short", long.as_str()]);
    let (lines, error) = read(&directory, "long.txt.gz", &gzip(text.as_bytes()));
    assert_eq!(error, None);
    assert_eq!(lines, [long.as_str(), "short", long.as_str()]);
}

#[test]
fn a_carriage_return_at_the_end_of_a_read_waits_for_its_line_feed() {
    let directory = scratch();
    // A line end that falls on either side of the bytes read at a time.
    for width in [2 * CHUNK - 2, 2 * CHUNK - 1, 2 * CHUNK, CHUNK - 1, CHUNK] {
        for ending in ["\r\n", "\r", "\n"] {
            let first = "y".repeat(width);
            let text = format!("{first}{ending}z{ending}");
            let (lines, error) = read(&directory, "edge.txt", text.as_bytes());
            assert_eq!(error, None);
            assert_eq!(lines, [first.as_str(), "z"], "{width} {ending:?}");
            let (lines, error) = read(&directory, "edge.txt.gz", &gzip(text.as_bytes()));
            assert_eq!(error, None);
            assert_eq!(lines, [first.as_str(), "z"], "{width} {ending:?} gzip");
        }
    }
}

#[test]
fn a_file_that_cannot_be_opened_gives_the_os_error() {
    let directory = scratch();
    let root = directory.path.to_str().unwrap();
    for name in ["absent.txt", "absent.gz"] {
        assert_eq!(
            lines_of(&format!("{root}/{name}")),
            (
                Vec::new(),
                cpython(&format!(
                    "[Errno 2] No such file or directory: '{root}/{name}'"
                ))
            )
        );
    }
    for name in ["dir.txt", "dir.gz"] {
        std::fs::create_dir(directory.path.join(name)).unwrap();
        assert_eq!(
            lines_of(&format!("{root}/{name}")),
            (
                Vec::new(),
                cpython(&format!("[Errno 21] Is a directory: '{root}/{name}'"))
            )
        );
    }
}

#[test]
fn gzip_members_are_read_in_turn_with_the_zero_bytes_between_them() {
    let directory = scratch();
    let good = gzip(TWO_LINES);
    let two_lines = ["a\tb", "c\td"].map(str::to_string).to_vec();
    let four_lines = [two_lines.clone(), two_lines.clone()].concat();
    assert_eq!(read(&directory, "empty.gz", b""), (Vec::new(), None));
    assert_eq!(
        read(&directory, "one.vcf.gz", &good),
        (two_lines.clone(), None)
    );
    assert_eq!(
        read(&directory, "one.vcf.bgz", &good),
        (two_lines.clone(), None)
    );
    assert_eq!(
        read(&directory, "two.gz", &[good.clone(), good.clone()].concat()),
        (four_lines.clone(), None)
    );
    assert_eq!(
        read(&directory, "padded.gz", &[&good[..], &[0, 0, 0]].concat()),
        (two_lines.clone(), None)
    );
    assert_eq!(
        read(
            &directory,
            "padded_two.gz",
            &[&good[..], &[0, 0], &good[..]].concat()
        ),
        (four_lines, None)
    );
    // A member whose header carries an extra field, a name, a comment and a header checksum.
    let mut flagged = vec![0x1f, 0x8b, 8, 4 | 8 | 16 | 2, 0, 0, 0, 0, 0, 3, 3, 0];
    flagged.extend(b"abcname\0comment\0\x12\x34");
    flagged.extend(&good[10..]);
    assert_eq!(
        read(&directory, "flags.gz", &flagged),
        (two_lines.clone(), None)
    );
    // A plain name is read as it is, whatever its bytes.
    let (lines, error) = read(&directory, "gzip.txt", &good);
    assert_eq!(error, None);
    assert!(lines[0].starts_with("\u{1f}\u{fffd}"));
    let path = directory.path.join("gzip.txt");
    let mut forced = TextLines::open_gzip(path.to_str().unwrap()).unwrap();
    assert_eq!(forced.next_line(), Ok(Some("a\tb")));
    assert!(forced.advance().unwrap());
    assert_eq!(forced.line(), "c\td");
    assert_eq!(forced.next_line(), Ok(None));
    let bytes = ByteSource::open(directory.path.join("two.gz").to_str().unwrap())
        .and_then(ByteSource::read_to_end)
        .unwrap();
    assert_eq!(bytes, [TWO_LINES, TWO_LINES].concat());
    let bytes = ByteSource::open_plain(directory.path.join("one.vcf.gz").to_str().unwrap())
        .and_then(ByteSource::read_to_end)
        .unwrap();
    assert_eq!(bytes, good);
}

#[test]
fn a_file_that_is_not_gzip_is_named_by_its_first_bytes() {
    let directory = scratch();
    let good = gzip(TWO_LINES);
    for (name, bytes, reason) in [
        ("one.gz", &b"\x1f"[..], "Not a gzipped file (b'\\x1f')"),
        ("plain.gz", b"##fileformat\n", "Not a gzipped file (b'##')"),
        ("quote.gz", b"'\"x", "Not a gzipped file (b'\\'\"')"),
        ("tab.gz", b"\t\\", "Not a gzipped file (b'\\t\\\\')"),
        ("hi.gz", b"\xff\x00", "Not a gzipped file (b'\\xff\\x00')"),
        (
            "method.gz",
            b"\x1f\x8b\x07\x00\x00\x00\x00\x00\x00\x03",
            "Unknown compression method",
        ),
        ("magic_only.gz", b"\x1f\x8b", CUT_SHORT),
        ("header_only.gz", &good[..10], CUT_SHORT),
        ("cut.gz", &good[..good.len() / 2], CUT_SHORT),
    ] {
        assert_eq!(
            read(&directory, name, bytes),
            (Vec::new(), cpython(reason)),
            "{name}"
        );
    }
    assert_eq!(
        CUT_SHORT,
        "Compressed file ended before the end-of-stream marker was reached"
    );
    assert_eq!(py_repr_bytes(b"it's"), "b\"it's\"");
    assert_eq!(py_repr_bytes(b"\r\n \x7f~"), "b'\\r\\n \\x7f~'");
    assert_eq!(py_repr_bytes(b""), "b''");
}

#[test]
fn a_member_that_ends_wrongly_gives_its_lines_and_then_the_reason() {
    let directory = scratch();
    let good = gzip(TWO_LINES);
    let two_lines = ["a\tb", "c\td"].map(str::to_string).to_vec();
    let mut wrong_checksum = good.clone();
    let at = wrong_checksum.len() - 8;
    wrong_checksum[at] ^= 0xff;
    let mut wrong_length = good.clone();
    let at = wrong_length.len() - 1;
    wrong_length[at] ^= 0x01;
    for (name, bytes, reason) in [
        ("cut_trailer.gz", good[..good.len() - 3].to_vec(), CUT_SHORT),
        (
            "garbage_after.gz",
            [&good[..], b"xy"].concat(),
            "Not a gzipped file (b'xy')",
        ),
        (
            "garbage_one.gz",
            [&good[..], b"x"].concat(),
            "Not a gzipped file (b'x')",
        ),
        (
            "crc.gz",
            wrong_checksum,
            "CRC check failed 0x2471a84f != 0x2471a8b0",
        ),
        ("len.gz", wrong_length, "Incorrect length of data produced"),
    ] {
        assert_eq!(
            read(&directory, name, &bytes),
            (two_lines.clone(), cpython(reason)),
            "{name}"
        );
    }
    // Without a last line end the last line waits for bytes that never come.
    let open_ended = gzip(b"a\tb\nc\td");
    assert_eq!(
        read(&directory, "open.gz", &open_ended[..open_ended.len() - 3]),
        (vec!["a\tb".to_string()], cpython(CUT_SHORT))
    );
    let carriage = gzip(b"a\tb\rc\td\r");
    assert_eq!(
        read(&directory, "carriage.gz", &carriage[..carriage.len() - 3]),
        (vec!["a\tb".to_string()], cpython(CUT_SHORT))
    );
}

#[test]
fn a_cut_gzip_file_gives_every_whole_line_of_what_inflates() {
    let directory = scratch();
    let mut text = Vec::new();
    for number in 0..200_000 {
        writeln!(text, "line {number}\tx").unwrap();
    }
    let whole = gzip(&text);
    let (lines, error) = read(&directory, "bigcut.gz", &whole[..whole.len() / 2]);
    assert_eq!(error, cpython(CUT_SHORT));
    assert!(lines.len() > 1000 && lines.len() < 200_000);
    for (number, line) in lines.iter().enumerate() {
        assert_eq!(*line, format!("line {number}\tx"));
    }
    let (lines, error) = read(&directory, "big.gz", &whole);
    assert_eq!((lines.len(), error), (200_000, None));
}

#[test]
fn corrupt_deflate_data_is_refused_in_zlibs_words() {
    let directory = scratch();
    let mut corrupt = gzip(TWO_LINES);
    corrupt[10] = 0xff;
    let (lines, read_error) = read(&directory, "corrupt.gz", &corrupt);
    assert_eq!(lines, Vec::<String>::new());
    // The reason is zlib's, on CPython's line.
    let named =
        ReadError::Cpython("Error -3 while decompressing data: invalid block type".to_string());
    assert_eq!(read_error, Some(named));
    // The tool's own words stand where the library names no reason.
    let error = ReadError::Tool(CORRUPT.to_string());
    assert_eq!(error.reason(), "corrupt deflate data in a gzip member");
    assert_eq!(
        not_readable_text("in.vcf.gz", &error),
        Failure::own(
            "in.vcf.gz is not a readable text file (corrupt deflate data in a gzip member)"
        )
    );
    assert_eq!(
        not_readable_text("in.vcf.gz", &ReadError::Cpython(CUT_SHORT.to_string())),
        failure(format!(
            "in.vcf.gz is not a readable text file ({CUT_SHORT})"
        ))
    );
    assert_eq!(
        ReadError::Cpython("x".to_string()).failure(|reason| format!("{reason} y")),
        failure("x y")
    );
}
