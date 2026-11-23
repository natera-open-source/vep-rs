// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The lines of a text as `read_text` and `text_lines` of `crate::input` give them
//! together: `read_text` reads in text mode, which turns `\r\n` and a lone `\r` into
//! `\n` before `text_lines` splits at `\n`, so each of the three ends a line, and a last
//! line without its end is a line.

use memchr::{memchr2, memrchr2};
use std::ops::Range;

/// The places of the lines of `bytes`, each without its end.
pub(crate) struct Lines<'a> {
    bytes: &'a [u8],
    at: usize,
}

pub(crate) fn lines(bytes: &[u8]) -> Lines<'_> {
    Lines { bytes, at: 0 }
}

impl Iterator for Lines<'_> {
    type Item = Range<usize>;

    fn next(&mut self) -> Option<Range<usize>> {
        if self.at >= self.bytes.len() {
            return None;
        }
        let start = self.at;
        let end = match memchr2(b'\n', b'\r', &self.bytes[start..]) {
            Some(i) => {
                let end = start + i;
                let crlf = self.bytes[end] == b'\r' && self.bytes.get(end + 1) == Some(&b'\n');
                self.at = end + 1 + usize::from(crlf);
                end
            }
            None => {
                self.at = self.bytes.len();
                self.bytes.len()
            }
        };
        Some(start..end)
    }
}

/// Where `bytes`, the start of a longer text, can be cut after a line end: the offset
/// after the last `\n`, or after the last `\r` that a byte other than `\n` follows.
/// `None` where `bytes` holds no such line end.
pub(crate) fn cut(bytes: &[u8]) -> Option<usize> {
    let last = memrchr2(b'\n', b'\r', bytes)?;
    if bytes[last] == b'\n' || last + 1 < bytes.len() {
        return Some(last + 1);
    }
    // A `\r` that ends `bytes` may be the first half of a `\r\n`.
    memrchr2(b'\n', b'\r', &bytes[..last]).map(|before| before + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texts(bytes: &[u8]) -> Vec<String> {
        lines(bytes)
            .map(|line| String::from_utf8(bytes[line].to_vec()).unwrap())
            .collect()
    }

    #[test]
    fn every_line_end_of_text_mode_ends_a_line() {
        assert_eq!(
            texts(b"#h\tx\r\nr1\t1\rr2\t2\nr3\t3"),
            ["#h\tx", "r1\t1", "r2\t2", "r3\t3"]
        );
        assert_eq!(
            texts(b"a\rb\r\nc\n\rd\r\r\ne"),
            ["a", "b", "c", "", "d", "", "e"]
        );
    }

    #[test]
    fn an_empty_line_before_the_end_stays_and_the_last_end_adds_none() {
        assert_eq!(texts(b"a\n\n"), ["a", ""]);
        assert_eq!(texts(b"a\r"), ["a"]);
        assert_eq!(texts(b""), Vec::<String>::new());
        assert_eq!(texts(b"\n"), [""]);
        assert_eq!(texts(b"\r\n\r"), ["", ""]);
    }

    #[test]
    fn a_text_is_cut_after_a_line_end_and_never_inside_one() {
        assert_eq!(cut(b"ab\ncd"), Some(3));
        assert_eq!(cut(b"ab\ncd\n"), Some(6));
        assert_eq!(cut(b"ab\r\ncd"), Some(4));
        assert_eq!(cut(b"ab\rcd"), Some(3));
        // The last byte may be half of a line end.
        assert_eq!(cut(b"ab\ncd\r"), Some(3));
        assert_eq!(cut(b"ab\rcd\r"), Some(3));
        assert_eq!(cut(b"ab\r\r"), Some(3));
        assert_eq!(cut(b"abcd\r"), None);
        assert_eq!(cut(b"abcd"), None);
        assert_eq!(cut(b""), None);
    }

    #[test]
    fn pieces_cut_at_line_ends_hold_the_lines_of_the_whole() {
        let mut bytes = Vec::new();
        for i in 0..200 {
            let end = ["\n", "\r\n", "\r"][i % 3];
            bytes.extend_from_slice(format!("row{i}\tvalue{end}").as_bytes());
        }
        bytes.extend_from_slice(b"last");
        let whole = texts(&bytes);
        assert_eq!(whole.len(), 201);
        for size in [1, 2, 7, 64, 1000] {
            let mut pieces = Vec::new();
            let mut from = 0;
            while from < bytes.len() {
                let mut to = (from + size).min(bytes.len());
                let at = loop {
                    if to == bytes.len() {
                        break to;
                    }
                    match cut(&bytes[from..to]) {
                        Some(at) => break from + at,
                        None => to = (to + size).min(bytes.len()),
                    }
                };
                pieces.extend(texts(&bytes[from..at]));
                from = at;
            }
            assert_eq!(pieces, whole, "pieces of {size} bytes");
        }
    }
}
