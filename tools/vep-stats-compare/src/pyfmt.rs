// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Text the outputs must reproduce character for character from CPython: `str.strip()`,
//! `int()` of a `str`, `repr()` of a `str` and of `bytes`, `json.dumps(ensure_ascii=True)`
//! string escapes, `str(pathlib.Path(arg))`, an argument decoded under `surrogateescape`, a
//! line written under `backslashreplace` and the text of an `OSError`.

use crate::unicode::{decimal_value, is_printable};
use std::cmp::Ordering;
use std::fmt::Write;
use std::io;

/// `str.isspace()`: the Unicode White_Space set plus U+001C to U+001F, which Python counts
/// through their bidirectional class and `char::is_whitespace` does not. The `\s` of `re`
/// is the same set.
pub fn is_py_space(c: char) -> bool {
    matches!(
        c,
        '\u{09}'..='\u{0d}'
            | '\u{1c}'..='\u{20}'
            | '\u{85}'
            | '\u{a0}'
            | '\u{1680}'
            | '\u{2000}'..='\u{200a}'
            | '\u{2028}'
            | '\u{2029}'
            | '\u{202f}'
            | '\u{205f}'
            | '\u{3000}'
    )
}

/// `str.strip()` with no argument.
pub fn py_strip(s: &str) -> &str {
    s.trim_matches(is_py_space)
}

/// The code points of the `str` Python decodes an argument's bytes to under
/// `surrogateescape`: UTF-8 where the bytes are UTF-8, U+DC80 to U+DCFF for each byte that
/// is not.
pub fn surrogateescape(bytes: &[u8]) -> Vec<u32> {
    let mut units = Vec::with_capacity(bytes.len());
    for chunk in bytes.utf8_chunks() {
        units.extend(chunk.valid().chars().map(u32::from));
        units.extend(chunk.invalid().iter().map(|b| 0xdc00 + u32::from(*b)));
    }
    units
}

/// The code points of a `str` as written to a stream whose error handler is
/// `backslashreplace`, as Python's standard error is: a lone surrogate as `\udcxx`, every
/// other code point as itself.
pub fn backslashreplace(units: &[u32]) -> String {
    let mut out = String::with_capacity(units.len());
    for &unit in units {
        match char::from_u32(unit) {
            Some(c) => out.push(c),
            None => write!(out, "\\u{unit:04x}").unwrap(),
        }
    }
    out
}

/// `repr()` of a `str`: single quotes unless the text holds a single quote and no double
/// one, the quote and the backslash escaped, `\t`, `\n` and `\r` by name, and every other
/// code point `str.isprintable()` refuses as `\xNN`, `\uNNNN` or `\UNNNNNNNN`.
pub fn py_repr(units: &[u32]) -> String {
    let quote = if units.contains(&0x27) && !units.contains(&0x22) {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(units.len() + 2);
    out.push(quote);
    for &unit in units {
        match unit {
            0x5c => out.push_str("\\\\"),
            0x09 => out.push_str("\\t"),
            0x0a => out.push_str("\\n"),
            0x0d => out.push_str("\\r"),
            _ if unit == u32::from(quote) => {
                out.push('\\');
                out.push(quote);
            }
            _ if is_printable(unit) => {
                out.push(char::from_u32(unit).expect("a printable code point is a scalar value"));
            }
            0..=0xff => write!(out, "\\x{unit:02x}").unwrap(),
            0x100..=0xffff => write!(out, "\\u{unit:04x}").unwrap(),
            _ => write!(out, "\\U{unit:08x}").unwrap(),
        }
    }
    out.push(quote);
    out
}

/// `repr()` of a `str` held as a Rust string.
pub fn py_repr_str(text: &str) -> String {
    let units: Vec<u32> = text.chars().map(u32::from).collect();
    py_repr(&units)
}

/// `repr()` of `bytes`: `b'...'`, double-quoted where the bytes hold a single quote and no
/// double one, the quote and the backslash escaped, `\t`, `\n` and `\r` by name, every other
/// byte outside `0x20..=0x7e` as `\xNN`.
pub fn py_bytes_repr(bytes: &[u8]) -> String {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(bytes.len() + 3);
    out.push('b');
    out.push(quote);
    for &byte in bytes {
        match byte {
            b'\\' => out.push_str("\\\\"),
            b'\t' => out.push_str("\\t"),
            b'\n' => out.push_str("\\n"),
            b'\r' => out.push_str("\\r"),
            _ if char::from(byte) == quote => {
                out.push('\\');
                out.push(quote);
            }
            b' '..=b'~' => out.push(char::from(byte)),
            _ => write!(out, "\\x{byte:02x}").unwrap(),
        }
    }
    out.push(quote);
    out
}

/// Appends `json.dumps(s)` with `ensure_ascii=True` for the code points of a `str`: the
/// seven short escapes, every other code point outside `0x20..=0x7e` as `\u` and four
/// lowercase hex digits, a surrogate pair above U+FFFF, a lone surrogate as itself.
pub fn push_json_units(out: &mut String, units: impl IntoIterator<Item = u32>) {
    out.push('"');
    for unit in units {
        match unit {
            0x22 => out.push_str("\\\""),
            0x5c => out.push_str("\\\\"),
            0x0a => out.push_str("\\n"),
            0x0d => out.push_str("\\r"),
            0x09 => out.push_str("\\t"),
            0x08 => out.push_str("\\b"),
            0x0c => out.push_str("\\f"),
            0x20..=0x7e => out.push(unit as u8 as char),
            0x10000.. => {
                let offset = unit - 0x10000;
                let (high, low) = (0xd800 + (offset >> 10), 0xdc00 + (offset & 0x3ff));
                write!(out, "\\u{high:04x}\\u{low:04x}").unwrap();
            }
            _ => write!(out, "\\u{unit:04x}").unwrap(),
        }
    }
    out.push('"');
}

/// `str(pathlib.PurePosixPath(arg))` on the argument's bytes, as Python applies it to the
/// surrogate-escaped `str` of any `argv` entry: repeated separators and `.` components
/// dropped, a trailing separator dropped, exactly two leading slashes kept, `..` left alone,
/// `.` for an empty result.
pub fn py_path_str(arg: &[u8]) -> Vec<u8> {
    let leading_slashes = arg.iter().take_while(|b| **b == b'/').count();
    let root: &[u8] = match leading_slashes {
        0 => b"",
        2 => b"//",
        _ => b"/",
    };
    let parts: Vec<&[u8]> = arg
        .split(|b| *b == b'/')
        .filter(|part| !part.is_empty() && *part != b".")
        .collect();
    let joined = [root, &parts.join(&b'/')].concat();
    if joined.is_empty() {
        b".".to_vec()
    } else {
        joined
    }
}

/// `os.strerror(errno)` of an operating-system error, or the error's own text where it
/// carries no error number.
pub fn strerror(error: &io::Error) -> String {
    let text = error.to_string();
    match error.raw_os_error() {
        Some(code) => text
            .strip_suffix(&format!(" (os error {code})"))
            .unwrap_or(&text)
            .to_string(),
        None => text,
    }
}

/// `str()` of the `OSError` Python raises for an operating-system error: `[Errno N] text`,
/// followed by `: ` and the `repr()` of the file name where the failing call named a file.
pub fn os_error_text(error: &io::Error, filename: Option<&[u8]>) -> String {
    let mut text = match error.raw_os_error() {
        Some(code) => format!("[Errno {code}] {}", strerror(error)),
        None => error.to_string(),
    };
    if let Some(filename) = filename {
        text.push_str(": ");
        text.push_str(&py_repr(&surrogateescape(filename)));
    }
    text
}

/// The digits `int()` reads before it refuses a text: `sys.get_int_max_str_digits()` as
/// CPython starts with it.
const INT_MAX_STR_DIGITS: usize = 4300;

/// A Python `int` of any size: a sign and the decimal digits of its magnitude without
/// leading zeros, zero being `0` and never negative.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PyInt {
    negative: bool,
    magnitude: String,
}

impl PyInt {
    /// The integer of a small count.
    pub fn from_u64(value: u64) -> PyInt {
        PyInt {
            negative: false,
            magnitude: value.to_string(),
        }
    }

    /// `int(text)` in base 10: the blanks of `str.isspace()` dropped at both ends except
    /// U+001C to U+001F, which `int()` refuses, one optional sign, decimal digits of any
    /// script with single underscores between them, at most 4300 digits. None where Python
    /// raises `ValueError`.
    pub fn parse(text: &str) -> Option<PyInt> {
        // CPython first rewrites the text to ASCII: a blank outside ASCII becomes a space,
        // a digit of another script its ASCII digit, and the text ends at the first other
        // character outside ASCII, which no rule below accepts.
        let mut ascii = Vec::with_capacity(text.len());
        for c in text.chars() {
            if u32::from(c) < 127 {
                ascii.push(c as u8);
            } else if is_py_space(c) {
                ascii.push(b' ');
            } else {
                ascii.push(b'0' + decimal_value(c)?);
            }
        }
        let is_blank = |b: &u8| matches!(*b, b'\t'..=b'\r' | b' ');
        let start = ascii.iter().position(|b| !is_blank(b))?;
        let end = ascii.iter().rposition(|b| !is_blank(b))? + 1;
        let body = &ascii[start..end];
        let (negative, digits) = match body.first()? {
            b'-' => (true, &body[1..]),
            b'+' => (false, &body[1..]),
            _ => (false, body),
        };
        if digits.first()? == &b'_' || digits.last()? == &b'_' {
            return None;
        }
        let mut magnitude = String::with_capacity(digits.len());
        let mut previous = 0u8;
        let mut count = 0usize;
        for &b in digits {
            match b {
                b'0'..=b'9' => {
                    count += 1;
                    if magnitude != "0" {
                        magnitude.push(char::from(b));
                    } else if b != b'0' {
                        magnitude.clear();
                        magnitude.push(char::from(b));
                    }
                }
                b'_' if previous != b'_' => {}
                _ => return None,
            }
            previous = b;
        }
        if count > INT_MAX_STR_DIGITS {
            return None;
        }
        Some(PyInt {
            negative: negative && magnitude != "0",
            magnitude,
        })
    }

    /// Whether the integer is below zero.
    pub fn is_negative(&self) -> bool {
        self.negative
    }

    /// Whether the integer is at least one.
    pub fn is_positive(&self) -> bool {
        !self.negative && self.magnitude != "0"
    }

    /// `str(n)`.
    pub fn decimal(&self) -> String {
        if self.negative {
            format!("-{}", self.magnitude)
        } else {
            self.magnitude.clone()
        }
    }

    /// `self <= other`.
    pub fn at_most(&self, other: &PyInt) -> bool {
        let order = match (self.negative, other.negative) {
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            (false, false) => cmp_magnitude(&self.magnitude, &other.magnitude),
            (true, true) => cmp_magnitude(&other.magnitude, &self.magnitude),
        };
        order != Ordering::Greater
    }

    /// `abs(self - other) <= bound` for a bound that is not negative.
    pub fn within(&self, other: &PyInt, bound: &PyInt) -> bool {
        let distance = if self.negative == other.negative {
            sub_magnitude(&self.magnitude, &other.magnitude)
        } else {
            add_magnitude(&self.magnitude, &other.magnitude)
        };
        cmp_magnitude(&distance, &bound.magnitude) != Ordering::Greater
    }
}

fn cmp_magnitude(a: &str, b: &str) -> Ordering {
    a.len().cmp(&b.len()).then_with(|| a.cmp(b))
}

fn add_magnitude(a: &str, b: &str) -> String {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    let mut digits = Vec::with_capacity(a.len().max(b.len()) + 1);
    let mut carry = 0u8;
    for i in 0..a.len().max(b.len()) {
        let digit = |n: &[u8]| n.len().checked_sub(i + 1).map_or(0, |at| n[at] - b'0');
        let sum = digit(a) + digit(b) + carry;
        digits.push(b'0' + sum % 10);
        carry = sum / 10;
    }
    if carry > 0 {
        digits.push(b'0' + carry);
    }
    digits.reverse();
    String::from_utf8(digits).expect("decimal digits")
}

/// The larger magnitude less the smaller.
fn sub_magnitude(a: &str, b: &str) -> String {
    let (large, small) = match cmp_magnitude(a, b) {
        Ordering::Less => (b.as_bytes(), a.as_bytes()),
        _ => (a.as_bytes(), b.as_bytes()),
    };
    let mut digits = Vec::with_capacity(large.len());
    let mut borrow = 0u8;
    for i in 0..large.len() {
        let top = large[large.len() - 1 - i] - b'0';
        let bottom = small
            .len()
            .checked_sub(i + 1)
            .map_or(0, |at| small[at] - b'0')
            + borrow;
        let (digit, next) = if top >= bottom {
            (top - bottom, 0)
        } else {
            (top + 10 - bottom, 1)
        };
        digits.push(b'0' + digit);
        borrow = next;
    }
    while digits.len() > 1 && digits.last() == Some(&b'0') {
        digits.pop();
    }
    digits.reverse();
    String::from_utf8(digits).expect("decimal digits")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn units(text: &str) -> Vec<u32> {
        text.chars().map(u32::from).collect()
    }

    #[test]
    fn strip_set_is_pythons_isspace() {
        let expected: Vec<u32> = [
            0x9, 0xa, 0xb, 0xc, 0xd, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x85, 0xa0, 0x1680, 0x2000,
            0x2001, 0x2002, 0x2003, 0x2004, 0x2005, 0x2006, 0x2007, 0x2008, 0x2009, 0x200a, 0x2028,
            0x2029, 0x202f, 0x205f, 0x3000,
        ]
        .to_vec();
        let actual: Vec<u32> = (0..=0x10ffffu32)
            .filter_map(char::from_u32)
            .filter(|c| is_py_space(*c))
            .map(u32::from)
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(py_strip("\u{1f} k\tv\u{a0}\u{3000}"), "k\tv");
        assert_eq!(py_strip(" \t x \u{200b}"), "x \u{200b}");
        assert_eq!(py_strip("\u{1c}\u{2028} "), "");
    }

    #[test]
    fn surrogateescape_maps_each_undecodable_byte_to_one_surrogate() {
        assert_eq!(surrogateescape(b"a\xc3\xa9"), [0x61, 0xe9]);
        assert_eq!(surrogateescape(b"a\xffb"), [0x61, 0xdcff, 0x62]);
        assert_eq!(surrogateescape(b"\xe2\x82A"), [0xdce2, 0xdc82, 0x41]);
        assert_eq!(surrogateescape(b"\xed\xb2\x80"), [0xdced, 0xdcb2, 0xdc80]);
        assert_eq!(surrogateescape(b"\xf0\x90"), [0xdcf0, 0xdc90]);
        assert_eq!(
            backslashreplace(&surrogateescape(b"out/\xffcaf\xc3\xa9.json")),
            "out/\\udcffcaf\u{e9}.json"
        );
    }

    #[test]
    fn str_repr_matches_python() {
        let cases = [
            ("x", "'x'"),
            ("", "''"),
            ("it's", "\"it's\""),
            ("it's \"x\"", "'it\\'s \"x\"'"),
            ("say \"x\"", "'say \"x\"'"),
            (" \u{feff}[A]", "' \\ufeff[A]'"),
            ("  k\tv", "'  k\\tv'"),
            ("caf\u{e9}\u{a0}\u{200b}", "'caf\u{e9}\\xa0\\u200b'"),
            (
                "a\\b\tc\u{7f}\0\u{a0}\u{1f600}\u{e0001}\u{d7ff}\u{e000}",
                "'a\\\\b\\tc\\x7f\\x00\\xa0\u{1f600}\\U000e0001\\ud7ff\\ue000'",
            ),
            ("\u{378}\u{61c}\u{ad}", "'\\u0378\\u061c\\xad'"),
            ("\r\n\u{1f}\u{85}", "'\\r\\n\\x1f\\x85'"),
            ("\u{4e00}\u{3000}\u{fffd}", "'\u{4e00}\\u3000\u{fffd}'"),
        ];
        for (text, expected) in cases {
            assert_eq!(py_repr(&units(text)), expected, "repr({text:?})");
            assert_eq!(py_repr_str(text), expected, "repr({text:?})");
        }
        assert_eq!(
            py_repr(&surrogateescape(b"/data/\xff's.txt")),
            "\"/data/\\udcff's.txt\""
        );
    }

    #[test]
    fn bytes_repr_matches_python() {
        let cases: [(&[u8], &str); 9] = [
            (b"ab", "b'ab'"),
            (b"\x1f", "b'\\x1f'"),
            (b"'", "b\"'\""),
            (b"'\"", "b'\\'\"'"),
            (b"\"\"", "b'\"\"'"),
            (b"\t\n", "b'\\t\\n'"),
            (b"\r\\", "b'\\r\\\\'"),
            (b"\x00\xff", "b'\\x00\\xff'"),
            (b"\x7f~", "b'\\x7f~'"),
        ];
        for (bytes, expected) in cases {
            assert_eq!(py_bytes_repr(bytes), expected, "repr({bytes:?})");
        }
    }

    #[test]
    fn json_string_escapes_match_ensure_ascii() {
        let mut out = String::new();
        push_json_units(
            &mut out,
            units("\u{1f}\u{7f}\u{e9}\u{1F600}\u{fffd}/\\\"\n\t\u{08}\u{0c}\r\0"),
        );
        assert_eq!(
            out,
            r#""\u001f\u007f\u00e9\ud83d\ude00\ufffd/\\\"\n\t\b\f\r\u0000""#
        );
        let mut out = String::new();
        push_json_units(&mut out, surrogateescape(b"row \xff\xc3\xa9"));
        assert_eq!(out, r#""row \udcff\u00e9""#);
    }

    #[test]
    fn path_str_matches_pathlib() {
        let cases = [
            ("./a//b/", "a/b"),
            ("/data/work/out_summary.txt", "/data/work/out_summary.txt"),
            ("//x", "//x"),
            ("///x", "/x"),
            ("", "."),
            (".", "."),
            ("a/./b", "a/b"),
            ("a/b/.", "a/b"),
            ("./", "."),
            ("x/../y", "x/../y"),
            ("report.gz/", "report.gz"),
        ];
        for (arg, expected) in cases {
            assert_eq!(
                py_path_str(arg.as_bytes()),
                expected.as_bytes(),
                "Path({arg:?})"
            );
        }
        assert_eq!(py_path_str(b"./ref\xff.txt/"), b"ref\xff.txt");
    }

    #[test]
    fn os_error_text_matches_python() {
        let missing = io::Error::from_raw_os_error(2);
        assert_eq!(strerror(&missing), "No such file or directory");
        assert_eq!(
            os_error_text(&missing, Some(b"out/none.txt")),
            "[Errno 2] No such file or directory: 'out/none.txt'"
        );
        assert_eq!(
            os_error_text(&io::Error::from_raw_os_error(21), Some(b"it's")),
            "[Errno 21] Is a directory: \"it's\""
        );
        assert_eq!(
            os_error_text(&io::Error::from_raw_os_error(5), None),
            "[Errno 5] Input/output error"
        );
        assert_eq!(
            strerror(&io::Error::from_raw_os_error(20)),
            "Not a directory"
        );
    }

    fn int(text: &str) -> Option<String> {
        PyInt::parse(text).map(|n| n.decimal())
    }

    #[test]
    fn int_reads_what_python_reads() {
        let read = [
            ("7", "7"),
            (" 7", "7"),
            ("7 ", "7"),
            ("+7", "7"),
            ("-7", "-7"),
            ("007", "7"),
            ("7_000", "7000"),
            ("0_0", "0"),
            ("\u{667}", "7"),
            ("\u{ff17}", "7"),
            ("\t7\n", "7"),
            ("\u{0b}7\u{0c}", "7"),
            ("\u{a0}7", "7"),
            ("\u{2000}7\u{3000}", "7"),
            ("\u{85}7", "7"),
            ("\u{2028}7", "7"),
            ("\u{1680} 7", "7"),
            ("-0", "0"),
            ("\u{1810}", "0"),
            ("\u{1d7d7}", "9"),
            ("\u{966}\u{967}\u{968}", "12"),
            ("7\u{667}", "77"),
            ("-\u{663}", "-3"),
            (
                "999999999999999999999999999999",
                "999999999999999999999999999999",
            ),
            (
                "-999999999999999999999999999999",
                "-999999999999999999999999999999",
            ),
            ("000000000000000000000000000000000000000012", "12"),
        ];
        for (text, expected) in read {
            assert_eq!(int(text).as_deref(), Some(expected), "int({text:?})");
        }
        let refused = [
            "",
            " ",
            "-",
            "+",
            "7.0",
            "1e3",
            "0x7",
            "\u{1c}7",
            "7\u{1f}",
            "\u{1e}",
            "--7",
            "+-7",
            "- 1",
            "+ 1",
            "1 2",
            "1_",
            "_1",
            "1__2",
            "-_1",
            "\u{b2}",
            "\u{bd}",
            "\u{4e00}",
            "7\0",
            "\x007",
            "\u{feff}7",
            "7\u{200b}",
            "7\u{7f}",
        ];
        for text in refused {
            assert_eq!(int(text), None, "int({text:?})");
        }
    }

    #[test]
    fn int_refuses_more_than_4300_digits() {
        assert_eq!(int(&"1".repeat(4300)), Some("1".repeat(4300)));
        assert_eq!(int(&"1".repeat(4301)), None);
        assert_eq!(int(&format!("{}1", "0".repeat(5000))), None);
        assert_eq!(int(&format!("{}1", "0".repeat(4299))).as_deref(), Some("1"));
        assert_eq!(
            int(&format!("{}1", "1_".repeat(4299))),
            Some("1".repeat(4300))
        );
        assert_eq!(int(&format!("{}1", "1_".repeat(4300))), None);
    }

    fn n(text: &str) -> PyInt {
        PyInt::parse(text).unwrap()
    }

    #[test]
    fn int_arithmetic_holds_beyond_machine_words() {
        assert!(n("963").within(&n("964"), &n("1")));
        assert!(!n("963").within(&n("965"), &n("1")));
        assert!(n("965").within(&n("963"), &n("2")));
        assert!(n("-1").within(&n("1"), &n("2")));
        assert!(!n("-1").within(&n("1"), &n("1")));
        assert!(n("-5").within(&n("-7"), &n("2")));
        assert!(n("1000").within(&n("999"), &n("1")));
        assert!(!n("1000").within(&n("998"), &n("1")));
        let big = n("340282366920938463463374607431768211456");
        let next = n("340282366920938463463374607431768211457");
        assert!(big.within(&next, &n("1")));
        assert!(!big.within(&next, &n("0")));
        assert!(!big.within(&n("-1"), &big));
        assert!(big.within(&n("-1"), &next));
        assert!(n("0").at_most(&n("0")));
        assert!(n("-3").at_most(&n("0")));
        assert!(n("-3").at_most(&n("-2")));
        assert!(!n("-2").at_most(&n("-3")));
        assert!(!n("1").at_most(&n("0")));
        assert!(n("1").at_most(&big));
        assert!(!next.at_most(&big));
        assert!(n("-1").is_negative() && !n("-0").is_negative());
        assert!(n("1").is_positive() && !n("0").is_positive() && !n("-1").is_positive());
        assert_eq!(n("-0"), n("0"));
        assert_eq!(n("+007"), PyInt::from_u64(7));
    }
}
