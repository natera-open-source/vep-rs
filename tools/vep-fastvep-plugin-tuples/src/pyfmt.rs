// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Text the outputs must reproduce character for character from CPython: an argument
//! decoded under `surrogateescape`, a line written under `backslashreplace`, `repr()` of a
//! `str`, `str(pathlib.Path(arg))`, `str.splitlines()`, the size of an `int()` and the text
//! and the class name of an `OSError`.

use crate::unicode::is_printable;
use std::cmp::Ordering;
use std::fmt::Write;
use std::io;

/// The code points of the `str` Python decodes bytes to under `surrogateescape`: UTF-8
/// where the bytes are UTF-8, U+DC80 to U+DCFF for each byte that is not.
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

/// Bytes of an argument or of the input as Python writes them to standard error.
pub fn shown(bytes: &[u8]) -> String {
    backslashreplace(&surrogateescape(bytes))
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

/// `str(pathlib.PurePosixPath(arg))` on the argument's bytes: repeated separators and `.`
/// components dropped, a trailing separator dropped, exactly two leading slashes kept, `..`
/// left alone, `.` for an empty result.
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

/// `str.splitlines()`: a line ends at `\n`, `\r`, `\r\n`, `\v`, `\f`, U+001C to U+001E,
/// U+0085, U+2028 or U+2029, and a text that ends at a line end holds no empty last line.
/// Text read with universal newlines splits the same way.
pub fn py_splitlines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = text.char_indices().peekable();
    while let Some((at, c)) = chars.next() {
        let ends = matches!(
            c,
            '\n' | '\r' | '\u{0b}' | '\u{0c}' | '\u{1c}'
                ..='\u{1e}' | '\u{85}' | '\u{2028}' | '\u{2029}'
        );
        if !ends {
            continue;
        }
        lines.push(&text[start..at]);
        start = at + c.len_utf8();
        if c == '\r' && chars.peek().is_some_and(|(_, next)| *next == '\n') {
            chars.next();
            start += 1;
        }
    }
    if start < text.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// The digits `int()` reads before it refuses a text: `sys.get_int_max_str_digits()` as
/// CPython starts with it.
pub const INT_MAX_STR_DIGITS: usize = 4300;

/// The order of two whole numbers written in ASCII decimal digits, leading zeros allowed.
pub fn cmp_decimal(a: &str, b: &str) -> Ordering {
    let (a, b) = (a.trim_start_matches('0'), b.trim_start_matches('0'));
    a.len().cmp(&b.len()).then_with(|| a.cmp(b))
}

/// A whole number written in ASCII decimal digits as a `u64`, the largest `u64` where it
/// does not fit.
pub fn saturating_u64(digits: &str) -> u64 {
    digits.bytes().fold(0u64, |value, digit| {
        value
            .saturating_mul(10)
            .saturating_add(u64::from(digit - b'0'))
    })
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

/// The last line of the traceback Python prints for an `OSError` no handler catches: the
/// class Python raises for the error number, then `[Errno N] text`.
pub fn uncaught_os_error(error: &io::Error) -> String {
    let class = match error.kind() {
        io::ErrorKind::BrokenPipe => "BrokenPipeError",
        io::ErrorKind::PermissionDenied => "PermissionError",
        io::ErrorKind::IsADirectory => "IsADirectoryError",
        io::ErrorKind::NotFound => "FileNotFoundError",
        io::ErrorKind::TimedOut => "TimeoutError",
        _ => "OSError",
    };
    match error.raw_os_error() {
        Some(code) => format!("{class}: [Errno {code}] {}", strerror(error)),
        None => format!("{class}: {error}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_byte_outside_utf8_is_shown_as_its_escaped_surrogate() {
        assert_eq!(shown(b"caf\xc3\xa9"), "caf\u{e9}");
        assert_eq!(shown(b"abs\xffent"), "abs\\udcffent");
        assert_eq!(shown(b"\xe2\x80x\xc3"), "\\udce2\\udc80x\\udcc3");
        assert_eq!(shown(b"\xed\xb2\x80"), "\\udced\\udcb2\\udc80");
    }

    #[test]
    fn repr_picks_its_quote_and_escapes_what_is_not_printable() {
        assert_eq!(py_repr_str("first"), "'first'");
        assert_eq!(py_repr_str("it's"), "\"it's\"");
        assert_eq!(py_repr_str("it's \"x\""), "'it\\'s \"x\"'");
        assert_eq!(py_repr_str("a\\b\tc\rd"), "'a\\\\b\\tc\\rd'");
        assert_eq!(py_repr_str("\u{7f}\u{a0}\u{e9}"), "'\\x7f\\xa0\u{e9}'");
        assert_eq!(py_repr_str("\u{200b}\u{1f600}"), "'\\u200b\u{1f600}'");
        assert_eq!(py_repr_str("\u{e0001}"), "'\\U000e0001'");
        assert_eq!(py_repr_str(""), "''");
    }

    #[test]
    fn a_path_is_normalized_as_pathlib_writes_it() {
        for (arg, expected) in [
            ("", "."),
            (".", "."),
            ("./", "."),
            ("a", "a"),
            ("./a//b/./c/", "a/b/c"),
            ("/", "/"),
            ("//", "//"),
            ("///", "/"),
            ("//a", "//a"),
            ("///a/", "/a"),
            ("a/../b", "a/../b"),
            ("..", ".."),
        ] {
            assert_eq!(py_path_str(arg.as_bytes()), expected.as_bytes(), "{arg:?}");
        }
        assert_eq!(py_path_str(b"./\xff/"), b"\xff");
    }

    #[test]
    fn lines_end_where_splitlines_ends_them() {
        assert_eq!(py_splitlines(""), Vec::<&str>::new());
        assert_eq!(py_splitlines("\n"), [""]);
        assert_eq!(py_splitlines("a"), ["a"]);
        assert_eq!(py_splitlines("a\nb\n"), ["a", "b"]);
        assert_eq!(py_splitlines("a\r\n\rb\n\n"), ["a", "", "b", ""]);
        assert_eq!(
            py_splitlines("a\u{0b}b\u{0c}c\u{1c}d\u{1d}e\u{1e}f\u{85}g\u{2028}h\u{2029}i"),
            ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
        );
        assert_eq!(py_splitlines("a\u{1f}b\tc d"), ["a\u{1f}b\tc d"]);
        assert_eq!(py_splitlines("\r"), [""]);
        assert_eq!(py_splitlines("x\r\r\ny"), ["x", "", "y"]);
    }

    #[test]
    fn decimal_texts_compare_and_saturate_as_whole_numbers() {
        assert_eq!(cmp_decimal("2", "10"), Ordering::Less);
        assert_eq!(cmp_decimal("010", "9"), Ordering::Greater);
        assert_eq!(cmp_decimal("000", "0"), Ordering::Equal);
        assert_eq!(cmp_decimal("7", "07"), Ordering::Equal);
        assert_eq!(
            cmp_decimal("18446744073709551616", "18446744073709551617"),
            Ordering::Less
        );
        assert_eq!(saturating_u64("0"), 0);
        assert_eq!(saturating_u64("0012"), 12);
        assert_eq!(saturating_u64("18446744073709551615"), u64::MAX);
        assert_eq!(saturating_u64("18446744073709551616"), u64::MAX);
        assert_eq!(saturating_u64("99999999999999999999999999"), u64::MAX);
    }

    #[test]
    fn an_operating_system_error_reads_as_python_writes_it() {
        let missing = io::Error::from_raw_os_error(2);
        assert_eq!(strerror(&missing), "No such file or directory");
        assert_eq!(
            uncaught_os_error(&missing),
            "FileNotFoundError: [Errno 2] No such file or directory"
        );
        let full = io::Error::from_raw_os_error(28);
        assert_eq!(
            uncaught_os_error(&full),
            "OSError: [Errno 28] No space left on device"
        );
        let pipe = io::Error::from_raw_os_error(32);
        assert_eq!(
            uncaught_os_error(&pipe),
            "BrokenPipeError: [Errno 32] Broken pipe"
        );
        let other = io::Error::other("cut short");
        assert_eq!(strerror(&other), "cut short");
        assert_eq!(uncaught_os_error(&other), "OSError: cut short");
    }
}
