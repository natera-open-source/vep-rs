// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The input as the text Python reads from a file opened in text mode: strict UTF-8 with
//! universal newlines, handed on in runs so that the input is never held whole.

use memchr::memchr;
use std::io::{self, Read};

/// The bytes read from the input at a time.
pub const RUN_BYTES: usize = 1 << 20;

/// Why the text of an input ends early.
#[derive(Debug)]
pub enum Failure {
    Read(io::Error),
    /// The offset, counted from 0, of the first byte that is not UTF-8.
    NotUtf8(u64),
}

/// The text of one input, a run at a time.
pub struct Text<R> {
    source: R,
    run_bytes: usize,
    /// The bytes held between two runs: at most three of a sequence that the end of the
    /// last run cut, then the run being read.
    raw: Vec<u8>,
    /// The offset in the input of `raw[0]`.
    offset: u64,
    /// The last byte handed on stood for a CR: an LF that follows belongs to it.
    after_cr: bool,
    ended: bool,
}

impl<R: Read> Text<R> {
    pub fn new(source: R) -> Text<R> {
        Text::with_run_bytes(source, RUN_BYTES)
    }

    pub fn with_run_bytes(source: R, run_bytes: usize) -> Text<R> {
        Text {
            source,
            run_bytes: run_bytes.max(1),
            raw: Vec::new(),
            offset: 0,
            after_cr: false,
            ended: false,
        }
    }

    /// Replaces `run` with the next run of text, every line ending in it written as LF.
    /// Returns false, and leaves `run` empty, at the end of the input.
    pub fn next_run(&mut self, run: &mut Vec<u8>) -> Result<bool, Failure> {
        run.clear();
        while run.is_empty() && !self.ended {
            let carried = self.raw.len();
            self.raw.resize(carried + self.run_bytes, 0);
            let read =
                read_full(&mut self.source, &mut self.raw[carried..]).map_err(Failure::Read)?;
            self.raw.truncate(carried + read);
            self.ended = read < self.run_bytes;
            let valid = match std::str::from_utf8(&self.raw) {
                Ok(_) => self.raw.len(),
                // A sequence that the end of the run cut is whole in the next run.
                Err(error) if error.error_len().is_none() && !self.ended => error.valid_up_to(),
                Err(error) => {
                    return Err(Failure::NotUtf8(self.offset + error.valid_up_to() as u64))
                }
            };
            translate(&self.raw[..valid], &mut self.after_cr, run);
            self.raw.drain(..valid);
            self.offset += valid as u64;
        }
        Ok(!run.is_empty())
    }
}

/// Reads until `buffer` is full or the source ends; returns the bytes read.
fn read_full(source: &mut impl Read, buffer: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buffer.len() {
        match source.read(&mut buffer[filled..]) {
            Ok(0) => break,
            Ok(read) => filled += read,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
    Ok(filled)
}

/// Appends `text` with CR LF and a lone CR each written as LF. `after_cr` carries a CR at
/// the end of one call to the LF that may open the next.
fn translate(text: &[u8], after_cr: &mut bool, out: &mut Vec<u8>) {
    let mut rest = text;
    if rest.is_empty() {
        return;
    }
    if std::mem::take(after_cr) {
        if let [b'\n', tail @ ..] = rest {
            rest = tail;
        }
    }
    while let Some(at) = memchr(b'\r', rest) {
        out.extend_from_slice(&rest[..at]);
        out.push(b'\n');
        rest = &rest[at + 1..];
        match rest {
            [b'\n', tail @ ..] => rest = tail,
            [] => *after_cr = true,
            _ => {}
        }
    }
    out.extend_from_slice(rest);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A source that hands out at most `step` bytes a read.
    struct Trickle<'a>(&'a [u8], usize);

    impl Read for Trickle<'_> {
        fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
            let taken = self.0.len().min(self.1).min(buffer.len());
            buffer[..taken].copy_from_slice(&self.0[..taken]);
            self.0 = &self.0[taken..];
            Ok(taken)
        }
    }

    fn text_of(input: &[u8], run_bytes: usize) -> Result<Vec<u8>, Failure> {
        let mut text = Text::with_run_bytes(Trickle(input, 3), run_bytes);
        let mut all = Vec::new();
        let mut run = Vec::new();
        while text.next_run(&mut run)? {
            all.extend_from_slice(&run);
        }
        Ok(all)
    }

    #[test]
    fn every_line_ending_is_written_as_a_line_feed_at_any_run_length() {
        let input = "a\r\nb\rc\n\r\n\r\r\nd\u{e9}\u{20ac}\u{1f600}\r".as_bytes();
        let expected = "a\nb\nc\n\n\n\nd\u{e9}\u{20ac}\u{1f600}\n".as_bytes();
        for run_bytes in 1..=input.len() + 1 {
            assert_eq!(
                text_of(input, run_bytes).unwrap(),
                expected,
                "runs of {run_bytes}"
            );
        }
    }

    #[test]
    fn text_without_a_carriage_return_is_handed_on_as_read() {
        let input = "file\tsource\n\n\u{2028}x\u{85}\x0b\x0c\x1c\0\n".as_bytes();
        for run_bytes in [1, 2, 5, 1 << 20] {
            assert_eq!(text_of(input, run_bytes).unwrap(), input);
        }
        assert_eq!(text_of(b"", 4).unwrap(), b"");
    }

    #[test]
    fn a_byte_outside_utf_8_is_named_by_its_offset_at_any_run_length() {
        let cases: [(&[u8], u64); 6] = [
            (b"ab\xffcd", 2),
            (b"\x1f\x8b\x08", 1),
            (b"abc\xe2\x82", 3),
            (b"abc\xe2\x82x", 3),
            (b"a\xed\xa0\x80", 1),
            (b"ok\n\xc0\xaf", 3),
        ];
        for (input, offset) in cases {
            for run_bytes in 1..=input.len() + 1 {
                match text_of(input, run_bytes) {
                    Err(Failure::NotUtf8(at)) => assert_eq!(at, offset, "{input:?}"),
                    other => panic!("{input:?}: {other:?}"),
                }
            }
        }
    }
}
