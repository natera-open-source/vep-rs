// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reads one statistics report as a stream: the file's bytes (the gzip members of a `.gz`
//! name inflated), strict UTF-8, the lines of `str.splitlines()`, the sections.
//!
//! The file is never held whole. What a whole read would report first is still reported
//! first: a file that cannot be read to its end outranks text that is not UTF-8, which
//! outranks text that is not a statistics report.

use crate::pyfmt::{os_error_text, py_bytes_repr};
use crate::report::{Parser, Section};
use flate2::{Crc, Decompress, FlushDecompress, Status};
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;

/// The line boundaries of `str.splitlines()` over text read with universal newlines: `\n`,
/// `\r`, `\r\n` as one, `\v`, `\f`, U+001C to U+001E, U+0085, U+2028 and U+2029. Text
/// arrives in pieces; a line is handed on without its ending, an empty last line never.
#[derive(Default)]
pub struct LineSplitter {
    pending: String,
    after_carriage_return: bool,
}

impl LineSplitter {
    /// Takes the next piece of text and hands on every line it completes.
    pub fn push(&mut self, text: &str, line: &mut impl FnMut(&str)) {
        let mut text = text;
        if self.after_carriage_return && !text.is_empty() {
            text = text.strip_prefix('\n').unwrap_or(text);
            self.after_carriage_return = false;
        }
        let bytes = text.as_bytes();
        let (mut start, mut at) = (0, 0);
        while at < bytes.len() {
            let width = match bytes[at] {
                b'\n' | 0x0b | 0x0c | 0x1c..=0x1e => 1,
                b'\r' => match bytes.get(at + 1) {
                    Some(b'\n') => 2,
                    Some(_) => 1,
                    None => {
                        self.after_carriage_return = true;
                        1
                    }
                },
                0xc2 if bytes[at + 1] == 0x85 => 2,
                0xe2 if bytes[at + 1] == 0x80 && matches!(bytes[at + 2], 0xa8 | 0xa9) => 3,
                _ => {
                    at += 1;
                    continue;
                }
            };
            if self.pending.is_empty() {
                line(&text[start..at]);
            } else {
                self.pending.push_str(&text[start..at]);
                line(&self.pending);
                self.pending.clear();
            }
            at += width;
            start = at;
        }
        self.pending.push_str(&text[start..]);
    }

    /// Hands on the last line where the text did not end in a line boundary.
    pub fn finish(&mut self, line: &mut impl FnMut(&str)) {
        if !self.pending.is_empty() {
            line(&self.pending);
            self.pending.clear();
        }
    }
}

/// The message of the `UnicodeDecodeError` CPython's strict UTF-8 decoder raises for the
/// undecodable bytes `bad` found `start` bytes into the text.
fn decode_error(bad: &[u8], start: u64, reason: &str) -> String {
    match bad {
        [byte] => {
            format!("'utf-8' codec can't decode byte 0x{byte:02x} in position {start}: {reason}")
        }
        _ => format!(
            "'utf-8' codec can't decode bytes in position {start}-{}: {reason}",
            start + bad.len() as u64 - 1
        ),
    }
}

/// Decodes a report's bytes as they arrive and parses the lines they complete.
#[derive(Default)]
struct Decoder {
    /// The bytes of a character the last piece left unfinished.
    carry: Vec<u8>,
    /// Bytes decoded so far, the carry not counted.
    decoded: u64,
    error: Option<String>,
    splitter: LineSplitter,
    parser: Parser,
}

impl Decoder {
    fn push(&mut self, bytes: &[u8]) {
        if self.error.is_some() {
            return;
        }
        let joined;
        let data = if self.carry.is_empty() {
            bytes
        } else {
            self.carry.extend_from_slice(bytes);
            joined = std::mem::take(&mut self.carry);
            &joined
        };
        let (valid, failure) = match std::str::from_utf8(data) {
            Ok(_) => (data.len(), None),
            Err(error) => (error.valid_up_to(), Some(error.error_len())),
        };
        let text = std::str::from_utf8(&data[..valid]).expect("the decoded prefix is UTF-8");
        let parser = &mut self.parser;
        self.splitter.push(text, &mut |line| parser.line(line));
        self.decoded += valid as u64;
        match failure {
            None => {}
            Some(None) => self.carry = data[valid..].to_vec(),
            Some(Some(length)) => {
                let bad = &data[valid..valid + length];
                let reason = match bad[0] {
                    0x80..=0xc1 | 0xf5..=0xff => "invalid start byte",
                    _ => "invalid continuation byte",
                };
                self.error = Some(decode_error(bad, self.decoded, reason));
            }
        }
    }

    fn finish(mut self) -> Result<Result<Vec<Section>, String>, String> {
        if let Some(error) = self.error {
            return Err(error);
        }
        if !self.carry.is_empty() {
            return Err(decode_error(
                &self.carry,
                self.decoded,
                "unexpected end of data",
            ));
        }
        let parser = &mut self.parser;
        self.splitter.finish(&mut |line| parser.line(line));
        Ok(self.parser.finish())
    }
}

/// The message of the `EOFError` Python's `gzip` raises for a file that ends inside a
/// member.
const ENDED_EARLY: &str = "Compressed file ended before the end-of-stream marker was reached";

const FHCRC: u8 = 2;
const FEXTRA: u8 = 4;
const FNAME: u8 = 8;
const FCOMMENT: u8 = 16;

fn read_error(error: io::Error) -> String {
    os_error_text(&error, None)
}

/// Up to `count` bytes: fewer only where the file ends.
fn read_up_to(reader: &mut impl Read, count: usize) -> Result<Vec<u8>, String> {
    let mut bytes = Vec::with_capacity(count);
    reader
        .take(count as u64)
        .read_to_end(&mut bytes)
        .map_err(read_error)?;
    Ok(bytes)
}

fn read_exact(reader: &mut impl Read, count: usize) -> Result<Vec<u8>, String> {
    let bytes = read_up_to(reader, count)?;
    if bytes.len() < count {
        return Err(ENDED_EARLY.to_string());
    }
    Ok(bytes)
}

/// Skips a zero-terminated header field; a file that ends inside it is left to the next
/// read to refuse.
fn skip_past_zero(reader: &mut impl BufRead) -> Result<(), String> {
    loop {
        let buffered = reader.fill_buf().map_err(read_error)?;
        let (skip, found) = match buffered.iter().position(|byte| *byte == 0) {
            Some(at) => (at + 1, true),
            None => (buffered.len(), buffered.is_empty()),
        };
        reader.consume(skip);
        if found {
            return Ok(());
        }
    }
}

/// One member's header, as Python's `gzip` reads it. False where the file ends before the
/// member's first byte.
fn read_member_header(reader: &mut impl BufRead) -> Result<bool, String> {
    let magic = read_up_to(reader, 2)?;
    if magic.is_empty() {
        return Ok(false);
    }
    if magic != [0x1f, 0x8b] {
        return Err(format!("Not a gzipped file ({})", py_bytes_repr(&magic)));
    }
    let fixed = read_exact(reader, 8)?;
    let (method, flags) = (fixed[0], fixed[1]);
    if method != 8 {
        return Err("Unknown compression method".to_string());
    }
    if flags & FEXTRA != 0 {
        let length = read_exact(reader, 2)?;
        read_exact(
            reader,
            usize::from(u16::from_le_bytes([length[0], length[1]])),
        )?;
    }
    if flags & FNAME != 0 {
        skip_past_zero(reader)?;
    }
    if flags & FCOMMENT != 0 {
        skip_past_zero(reader)?;
    }
    if flags & FHCRC != 0 {
        read_exact(reader, 2)?;
    }
    Ok(true)
}

/// Inflates every member of a gzip file into `sink`, as Python's `gzip` reads one: members
/// one after another, each held to its CRC and its length, zero bytes between and after
/// them skipped, an empty file an empty text. The error is the text of the exception
/// Python raises for the same file.
fn read_gzip(mut reader: impl BufRead, sink: &mut impl FnMut(&[u8])) -> Result<(), String> {
    // Under 260 bytes of room the inflater stays out of its fast loop, which names every
    // fault it meets `repeated call with bad state` where zlib names the fault itself.
    let mut output = [0u8; 256];
    while read_member_header(&mut reader)? {
        let mut inflater = Decompress::new(false);
        let mut crc = Crc::new();
        loop {
            let input = reader.fill_buf().map_err(read_error)?;
            let at_end = input.is_empty();
            let (before_in, before_out) = (inflater.total_in(), inflater.total_out());
            let status = inflater
                .decompress(input, &mut output, FlushDecompress::None)
                .map_err(|error| {
                    let fault = error.message().unwrap_or("invalid input data");
                    format!("Error -3 while decompressing data: {fault}")
                })?;
            let consumed = (inflater.total_in() - before_in) as usize;
            let produced = (inflater.total_out() - before_out) as usize;
            reader.consume(consumed);
            crc.update(&output[..produced]);
            sink(&output[..produced]);
            if status == Status::StreamEnd {
                break;
            }
            if at_end && produced == 0 {
                return Err(ENDED_EARLY.to_string());
            }
        }
        let trailer = read_exact(&mut reader, 8)?;
        let stored_crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
        let stored_length = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
        if stored_crc != crc.sum() {
            return Err(format!(
                "CRC check failed {stored_crc:#x} != {:#x}",
                crc.sum()
            ));
        }
        if stored_length != crc.amount() {
            return Err("Incorrect length of data produced".to_string());
        }
        loop {
            let buffered = reader.fill_buf().map_err(read_error)?;
            let zeros = buffered.iter().take_while(|byte| **byte == 0).count();
            let more = zeros > 0 && zeros == buffered.len();
            reader.consume(zeros);
            if !more {
                break;
            }
        }
    }
    Ok(())
}

fn read_plain(mut file: impl Read, sink: &mut impl FnMut(&[u8])) -> Result<(), String> {
    let mut buffer = vec![0u8; 1 << 16];
    loop {
        match file.read(&mut buffer) {
            Ok(0) => return Ok(()),
            Ok(count) => sink(&buffer[..count]),
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(read_error(error)),
        }
    }
}

/// Reads the statistics report at `path`, gunzipped when the name ends in `.gz`.
///
/// The outer error is why the file cannot be read, in the words of the exception Python
/// raises for it (an `OSError`, a `UnicodeDecodeError`, an `EOFError`, a `zlib.error`); the
/// inner one is why its text is not a statistics report.
pub fn read_report(path: &[u8]) -> Result<Result<Vec<Section>, String>, String> {
    let mut file = File::open(Path::new(OsStr::from_bytes(path)))
        .map_err(|error| os_error_text(&error, Some(path)))?;
    // Python refuses a directory as it opens it; the kernel only as it is read.
    if file.metadata().is_ok_and(|metadata| metadata.is_dir()) {
        let error = file
            .read(&mut [0u8; 1])
            .err()
            .unwrap_or_else(|| io::ErrorKind::IsADirectory.into());
        return Err(os_error_text(&error, Some(path)));
    }
    let mut decoder = Decoder::default();
    let mut sink = |bytes: &[u8]| decoder.push(bytes);
    if path.ends_with(b".gz") {
        read_gzip(BufReader::new(file), &mut sink)?;
    } else {
        read_plain(file, &mut sink)?;
    }
    decoder.finish()
}

/// The sections of a report held as text.
#[cfg(test)]
pub fn parse_text(text: &str) -> Result<Vec<Section>, String> {
    let mut parser = Parser::default();
    let mut splitter = LineSplitter::default();
    splitter.push(text, &mut |line| parser.line(line));
    splitter.finish(&mut |line| parser.line(line));
    parser.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lines_of(pieces: &[&str]) -> Vec<String> {
        let mut lines = Vec::new();
        let mut splitter = LineSplitter::default();
        for piece in pieces {
            splitter.push(piece, &mut |line: &str| lines.push(line.to_string()));
        }
        splitter.finish(&mut |line: &str| lines.push(line.to_string()));
        lines
    }

    #[test]
    fn lines_end_where_splitlines_ends_them() {
        let text = "[A]\r\nk\t1\rk2\t2\r\n\r\n[B]\u{0b}q\t1\u{0c}z\t2\u{1c}y\t3\u{1d}x\t4\u{1e}w\t5\u{85}v\t6\u{2028}u\t7\u{2029}t\t8";
        let expected = [
            "[A]", "k\t1", "k2\t2", "", "[B]", "q\t1", "z\t2", "y\t3", "x\t4", "w\t5", "v\t6",
            "u\t7", "t\t8",
        ];
        assert_eq!(lines_of(&[text]), expected);
        // The same lines however the text is cut, a `\r\n` cut in two included.
        for cut in 0..=text.len() {
            if text.is_char_boundary(cut) {
                assert_eq!(
                    lines_of(&[&text[..cut], &text[cut..]]),
                    expected,
                    "cut at {cut}"
                );
            }
        }
        assert_eq!(lines_of(&[""]), Vec::<String>::new());
        assert_eq!(lines_of(&["a\n"]), ["a"]);
        assert_eq!(lines_of(&["a\n\n"]), ["a", ""]);
        assert_eq!(lines_of(&["a\r", "\r", "\nb"]), ["a", "", "b"]);
        assert_eq!(lines_of(&["a\r", "", "\nb\r"]), ["a", "b"]);
        assert_eq!(
            lines_of(&["\u{1f}\n\u{0b}\u{1c}  k\tv\n"]),
            ["\u{1f}", "", "", "  k\tv"]
        );
    }

    fn decode(pieces: &[&[u8]]) -> Result<Result<Vec<Section>, String>, String> {
        let mut decoder = Decoder::default();
        for piece in pieces {
            decoder.push(piece);
        }
        decoder.finish()
    }

    #[test]
    fn undecodable_bytes_are_named_as_python_names_them() {
        let cases: [(&[u8], &str); 17] = [
            (
                b"[A]\n\xff\n",
                "byte 0xff in position 4: invalid start byte",
            ),
            (
                b"[A]\nab\xc3",
                "byte 0xc3 in position 6: unexpected end of data",
            ),
            (
                b"[A]\n\xe2\x82",
                "bytes in position 4-5: unexpected end of data",
            ),
            (
                b"[A]\n\xe2\x82A",
                "bytes in position 4-5: invalid continuation byte",
            ),
            (
                b"\xe0\x80\x80",
                "byte 0xe0 in position 0: invalid continuation byte",
            ),
            (
                b"\xf0\x90\x80",
                "bytes in position 0-2: unexpected end of data",
            ),
            (
                b"\xf0\x90A",
                "bytes in position 0-1: invalid continuation byte",
            ),
            (
                b"\xf0\x90\x80A",
                "bytes in position 0-2: invalid continuation byte",
            ),
            (
                b"\xf4\x90\x80\x80",
                "byte 0xf4 in position 0: invalid continuation byte",
            ),
            (
                b"\xed\xa0\x80",
                "byte 0xed in position 0: invalid continuation byte",
            ),
            (b"\xc0\x80", "byte 0xc0 in position 0: invalid start byte"),
            (b"\xf5\x80", "byte 0xf5 in position 0: invalid start byte"),
            (b"a\x80", "byte 0x80 in position 1: invalid start byte"),
            (b"\xf0\x90", "bytes in position 0-1: unexpected end of data"),
            (b"\xe2", "byte 0xe2 in position 0: unexpected end of data"),
            (b"\xf0", "byte 0xf0 in position 0: unexpected end of data"),
            (
                b"\xf0\x80",
                "byte 0xf0 in position 0: invalid continuation byte",
            ),
        ];
        for (bytes, expected) in cases {
            let expected = format!("'utf-8' codec can't decode {expected}");
            assert_eq!(decode(&[bytes]), Err(expected.clone()), "{bytes:?}");
            // The position counts from the start of the text however the bytes arrive.
            for cut in 0..=bytes.len() {
                let pieces = [&bytes[..cut], &bytes[cut..]];
                assert_eq!(
                    decode(&pieces),
                    Err(expected.clone()),
                    "{bytes:?} cut at {cut}"
                );
            }
        }
    }

    #[test]
    fn a_character_cut_between_two_pieces_is_one_character() {
        let text = "[caf\u{e9} \u{1f600}]\nk\u{2028}\u{4e00}\t1\n";
        let whole = parse_text(text).unwrap();
        assert_eq!(whole.len(), 1);
        assert_eq!(whole[0].title, "caf\u{e9} \u{1f600}");
        assert_eq!(
            whole[0].rows,
            [
                ("k".to_string(), None),
                ("\u{4e00}".to_string(), Some("1".to_string()))
            ]
        );
        let bytes = text.as_bytes();
        for cut in 0..=bytes.len() {
            let sections = decode(&[&bytes[..cut], &bytes[cut..]]).unwrap().unwrap();
            assert_eq!(sections, whole, "cut at {cut}");
        }
    }

    #[test]
    fn text_that_is_not_utf8_outranks_text_that_is_not_a_report() {
        assert_eq!(
            decode(&[b"k\t1\n[A]\n"]),
            Ok(Err(
                "line 1 precedes the first [section] header: 'k\\t1'".to_string()
            ))
        );
        assert_eq!(
            decode(&[b"k\t1\n[A]\n\xff"]),
            Err(
                "'utf-8' codec can't decode byte 0xff in position 8: invalid start byte"
                    .to_string()
            )
        );
    }
}
