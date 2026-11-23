// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Text the tool reproduces character for character from CPython: `str(pathlib.Path(arg))`,
//! `repr()` of a `str` and of a `bytes`, the `OSError` and `UnicodeDecodeError` lines, the
//! decimal digits `\d` matches and `float()` reads, and what the error stream prints of a
//! byte that is not UTF-8.

use std::io;
use std::str::Utf8Error;

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

/// The zero of every run of ten decimal digits (general category Nd) in Unicode 16.0, in
/// code point order; each run holds the digits 0 to 9 in order.
const DECIMAL_DIGIT_ZEROS: [u32; 76] = [
    0x30, 0x660, 0x6F0, 0x7C0, 0x966, 0x9E6, 0xA66, 0xAE6, 0xB66, 0xBE6, 0xC66, 0xCE6, 0xD66,
    0xDE6, 0xE50, 0xED0, 0xF20, 0x1040, 0x1090, 0x17E0, 0x1810, 0x1946, 0x19D0, 0x1A80, 0x1A90,
    0x1B50, 0x1BB0, 0x1C40, 0x1C50, 0xA620, 0xA8D0, 0xA900, 0xA9D0, 0xA9F0, 0xAA50, 0xABF0, 0xFF10,
    0x104A0, 0x10D30, 0x10D40, 0x11066, 0x110F0, 0x11136, 0x111D0, 0x112F0, 0x11450, 0x114D0,
    0x11650, 0x116C0, 0x116D0, 0x116DA, 0x11730, 0x118E0, 0x11950, 0x11BF0, 0x11C50, 0x11D50,
    0x11DA0, 0x11F50, 0x16130, 0x16A60, 0x16AC0, 0x16B50, 0x16D70, 0x1CCF0, 0x1D7CE, 0x1D7D8,
    0x1D7E2, 0x1D7EC, 0x1D7F6, 0x1E140, 0x1E2F0, 0x1E4F0, 0x1E5F1, 0x1E950, 0x1FBF0,
];

/// The value of a character `\d` matches in a `str` pattern, which is the digit `float()`
/// reads it as; `None` for any other character.
pub fn decimal_digit(c: char) -> Option<u8> {
    let code = c as u32;
    let run = DECIMAL_DIGIT_ZEROS.partition_point(|zero| *zero <= code);
    let value = code - DECIMAL_DIGIT_ZEROS[run.checked_sub(1)?];
    (value < 10).then_some(value as u8)
}

/// Whether `repr()` escapes a character above U+007F. `str.isprintable()` is false for the
/// control, format, separator, surrogate, private-use and unassigned categories; the ones
/// a path or a cell can plausibly hold are listed, and any other character prints as
/// written.
fn repr_escapes(c: char) -> bool {
    c.is_control()
        || c.is_whitespace()
        || matches!(
            c as u32,
            0xAD | 0x200B..=0x200F | 0x202A..=0x202E | 0x2060..=0x206F | 0xFEFF
        )
}

/// Appends `repr(s)` for the `str` Python decodes `bytes` to under `surrogateescape`, as an
/// `OSError` line names a file and a `ValueError` line a token: single quotes unless the
/// text holds one and no double quote, the backslash, the quote, tab, newline and carriage
/// return as their two-character escapes, any other control character and a byte that is
/// not UTF-8 as a hex escape.
pub fn push_py_repr(out: &mut Vec<u8>, bytes: &[u8]) {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        b'"'
    } else {
        b'\''
    };
    out.push(quote);
    for chunk in bytes.utf8_chunks() {
        for c in chunk.valid().chars() {
            match c {
                '\\' => out.extend_from_slice(b"\\\\"),
                '\t' => out.extend_from_slice(b"\\t"),
                '\n' => out.extend_from_slice(b"\\n"),
                '\r' => out.extend_from_slice(b"\\r"),
                c if c as u32 == u32::from(quote) => out.extend_from_slice(&[b'\\', quote]),
                ' '..='~' => out.push(c as u8),
                c if (c as u32) < 0x80 || repr_escapes(c) => {
                    let code = c as u32;
                    let escape = match code {
                        0..=0xFF => format!("\\x{code:02x}"),
                        0x100..=0xFFFF => format!("\\u{code:04x}"),
                        _ => format!("\\U{code:08x}"),
                    };
                    out.extend_from_slice(escape.as_bytes());
                }
                c => out.extend_from_slice(c.encode_utf8(&mut [0; 4]).as_bytes()),
            }
        }
        for byte in chunk.invalid() {
            out.extend_from_slice(format!("\\udc{byte:02x}").as_bytes());
        }
    }
    out.push(quote);
}

/// `repr(b)` of a `bytes`: `b'...'` with the same quote choice and two-character escapes
/// as a `str`, every other byte outside the printable ASCII range as `\xNN`.
pub fn py_bytes_repr(bytes: &[u8]) -> String {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::from("b");
    out.push(quote);
    for &byte in bytes {
        match byte {
            b'\\' => out.push_str("\\\\"),
            b'\t' => out.push_str("\\t"),
            b'\n' => out.push_str("\\n"),
            b'\r' => out.push_str("\\r"),
            byte if byte == quote as u8 => {
                out.push('\\');
                out.push(quote);
            }
            b' '..=b'~' => out.push(byte as char),
            byte => out.push_str(&format!("\\x{byte:02x}")),
        }
    }
    out.push(quote);
    out
}

/// Text as Python's error stream prints it: UTF-8 as written, a byte that is not UTF-8 as
/// the `\udcNN` escape of the surrogate it was decoded to.
pub fn push_error_stream_text(out: &mut Vec<u8>, bytes: &[u8]) {
    for chunk in bytes.utf8_chunks() {
        out.extend_from_slice(chunk.valid().as_bytes());
        for byte in chunk.invalid() {
            out.extend_from_slice(format!("\\udc{byte:02x}").as_bytes());
        }
    }
}

/// `str(OSError)` without a file name: `[Errno N] <strerror>`. An error that carries no
/// errno prints its own text.
pub fn os_error_text(error: &io::Error) -> String {
    match error.raw_os_error() {
        Some(errno) => {
            let text = error.to_string();
            let strerror = text
                .strip_suffix(&format!(" (os error {errno})"))
                .unwrap_or(&text);
            format!("[Errno {errno}] {strerror}")
        }
        None => error.to_string(),
    }
}

/// `str(OSError)` of a failed `open()`: `[Errno N] <strerror>: <repr of the name>`.
pub fn os_error_line(error: &io::Error, name: &[u8]) -> Vec<u8> {
    let mut line = os_error_text(error).into_bytes();
    line.extend_from_slice(b": ");
    push_py_repr(&mut line, name);
    line
}

/// `str(UnicodeDecodeError)` of the strict UTF-8 codec for the first invalid sequence of
/// `bytes`, which sit `offset` bytes into the decoded data: an invalid start byte, an
/// invalid continuation byte (the sequence up to the byte that broke it) or an unexpected
/// end of data (the truncated sequence).
pub fn utf8_error_text(bytes: &[u8], offset: u64, error: &Utf8Error) -> String {
    let start = error.valid_up_to();
    let first = bytes[start];
    let (length, reason) = match error.error_len() {
        None => (bytes.len() - start, "unexpected end of data"),
        Some(1) if !(0xC2..=0xF4).contains(&first) => (1, "invalid start byte"),
        Some(length) => (length, "invalid continuation byte"),
    };
    let position = offset + start as u64;
    if length == 1 {
        format!("'utf-8' codec can't decode byte 0x{first:02x} in position {position}: {reason}")
    } else {
        let last = position + length as u64 - 1;
        format!("'utf-8' codec can't decode bytes in position {position}-{last}: {reason}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path(arg: &str) -> String {
        String::from_utf8(py_path_str(arg.as_bytes())).unwrap()
    }

    #[test]
    fn a_path_is_normalised_as_pathlib_prints_it() {
        assert_eq!(path("out.parquet"), "out.parquet");
        assert_eq!(path("./out.parquet/"), "out.parquet");
        assert_eq!(path("a//b/./c/"), "a/b/c");
        assert_eq!(path("/a/../b"), "/a/../b");
        assert_eq!(path("//a"), "//a");
        assert_eq!(path("///a"), "/a");
        assert_eq!(path("."), ".");
        assert_eq!(path("./"), ".");
        assert_eq!(path("/"), "/");
        assert_eq!(py_path_str(b"./p\xff/"), b"p\xff");
    }

    #[test]
    fn decimal_digits_are_the_runs_of_ten() {
        assert_eq!(decimal_digit('0'), Some(0));
        assert_eq!(decimal_digit('9'), Some(9));
        assert_eq!(decimal_digit('/'), None);
        assert_eq!(decimal_digit(':'), None);
        assert_eq!(decimal_digit('\u{661}'), Some(1));
        assert_eq!(decimal_digit('\u{669}'), Some(9));
        assert_eq!(decimal_digit('\u{66A}'), None);
        assert_eq!(decimal_digit('\u{FF17}'), Some(7));
        assert_eq!(decimal_digit('\u{1D7CE}'), Some(0));
        assert_eq!(decimal_digit('\u{1D7FF}'), Some(9));
        assert_eq!(decimal_digit('\u{1FBF9}'), Some(9));
        assert_eq!(decimal_digit('\u{1FBFA}'), None);
        assert_eq!(decimal_digit('\u{BD}'), None);
        assert_eq!(decimal_digit('a'), None);
        assert_eq!(DECIMAL_DIGIT_ZEROS.len() * 10, 760);
        assert!(DECIMAL_DIGIT_ZEROS.windows(2).all(|w| w[0] + 10 <= w[1]));
    }

    fn repr(bytes: &[u8]) -> String {
        let mut out = Vec::new();
        push_py_repr(&mut out, bytes);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn a_str_repr_quotes_and_escapes_as_python_does() {
        assert_eq!(repr(b"out/x.tab"), "'out/x.tab'");
        assert_eq!(repr(b""), "''");
        assert_eq!(repr(b"it's"), "\"it's\"");
        assert_eq!(repr(b"it's \"x\""), "'it\\'s \"x\"'");
        assert_eq!(repr(b"a\\b\tc\nd\re"), "'a\\\\b\\tc\\nd\\re'");
        assert_eq!(repr(b"\x00\x1f\x7f"), "'\\x00\\x1f\\x7f'");
        assert_eq!(repr("\u{131}nf\n".as_bytes()), "'\u{131}nf\\n'");
        assert_eq!(repr("\u{a0}\u{ad}".as_bytes()), "'\\xa0\\xad'");
        assert_eq!(repr("\u{200b}\u{2028}".as_bytes()), "'\\u200b\\u2028'");
        assert_eq!(repr(b"p\xff"), "'p\\udcff'");
    }

    #[test]
    fn a_bytes_repr_quotes_and_escapes_as_python_does() {
        assert_eq!(py_bytes_repr(b"ju"), "b'ju'");
        assert_eq!(py_bytes_repr(b"\x1f"), "b'\\x1f'");
        assert_eq!(py_bytes_repr(b"'\x00"), "b\"'\\x00\"");
        assert_eq!(py_bytes_repr(b"\"'"), "b'\"\\''");
        assert_eq!(py_bytes_repr(b"\\\n"), "b'\\\\\\n'");
        assert_eq!(py_bytes_repr(b"\xff~"), "b'\\xff~'");
    }

    #[test]
    fn an_os_error_line_names_the_errno_and_the_file() {
        let missing = io::Error::from_raw_os_error(2);
        assert_eq!(
            os_error_text(&missing),
            "[Errno 2] No such file or directory"
        );
        assert_eq!(
            os_error_line(&missing, b"nodir/out.tab"),
            b"[Errno 2] No such file or directory: 'nodir/out.tab'"
        );
        let other = io::Error::other("no errno");
        assert_eq!(os_error_text(&other), "no errno");
    }

    #[test]
    fn the_error_stream_escapes_a_byte_that_is_not_utf8() {
        let mut out = Vec::new();
        push_error_stream_text(&mut out, b"a\xffb\xc3\xa9");
        assert_eq!(out, b"a\\udcffb\xc3\xa9");
    }

    fn decode_error(bytes: &[u8], offset: u64) -> String {
        let error = std::str::from_utf8(bytes).unwrap_err();
        utf8_error_text(bytes, offset, &error)
    }

    #[test]
    fn a_utf8_error_reads_as_the_codec_reports_it() {
        let cases: &[(&[u8], &str)] = &[
            (
                b"\xff\xfe\x00",
                "byte 0xff in position 0: invalid start byte",
            ),
            (b"ab\x80", "byte 0x80 in position 2: invalid start byte"),
            (b"\xc0\x80", "byte 0xc0 in position 0: invalid start byte"),
            (
                b"\xf5\x80\x80\x80",
                "byte 0xf5 in position 0: invalid start byte",
            ),
            (
                b"\xc3(",
                "byte 0xc3 in position 0: invalid continuation byte",
            ),
            (
                b"\xe0\x80\x80",
                "byte 0xe0 in position 0: invalid continuation byte",
            ),
            (
                b"\xed\xa0\x80",
                "byte 0xed in position 0: invalid continuation byte",
            ),
            (
                b"\xf4\x90\x80\x80",
                "byte 0xf4 in position 0: invalid continuation byte",
            ),
            (
                b"a\xe1\x80(",
                "bytes in position 1-2: invalid continuation byte",
            ),
            (
                b"\xf1\x80\x80(",
                "bytes in position 0-2: invalid continuation byte",
            ),
            (
                b"\xf1\x80(",
                "bytes in position 0-1: invalid continuation byte",
            ),
            (b"ab\xc3", "byte 0xc3 in position 2: unexpected end of data"),
            (b"\xe1\x80", "bytes in position 0-1: unexpected end of data"),
            (
                b"\xf1\x80\x80",
                "bytes in position 0-2: unexpected end of data",
            ),
            (
                b"\xe1(",
                "byte 0xe1 in position 0: invalid continuation byte",
            ),
        ];
        for (bytes, tail) in cases {
            assert_eq!(
                decode_error(bytes, 0),
                format!("'utf-8' codec can't decode {tail}"),
                "{bytes:?}"
            );
        }
        assert_eq!(
            decode_error(b"ab\x80", 1 << 33),
            "'utf-8' codec can't decode byte 0x80 in position 8589934594: invalid start byte"
        );
    }
}
