// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The CPython rules the tool's lines follow beside those of `pyfmt`: `repr()` and `str()` of
//! a string, a list and a JSON value, the text of a `UnicodeDecodeError` and of an `OSError`,
//! a file read as `Path.read_text(encoding="utf-8")` reads it, `str.splitlines()`, `int()`
//! for an option's count, `str(pathlib.Path(arg))` with its `/` join and
//! `hashlib.sha256(data).hexdigest()`.

use crate::json::Value;
use crate::pyfmt::{is_py_space, py_float_repr, py_path_str};
use sha2::{Digest, Sha256};
use std::fmt::Write;

/// The assigned code points from U+0080 that `str.isprintable()` refuses under Unicode 16:
/// the control, format, private-use and separator classes, as inclusive ranges.
const UNPRINTABLE: [(u32, u32); 27] = [
    (0x80, 0xA0),
    (0xAD, 0xAD),
    (0x600, 0x605),
    (0x61C, 0x61C),
    (0x6DD, 0x6DD),
    (0x70F, 0x70F),
    (0x890, 0x891),
    (0x8E2, 0x8E2),
    (0x1680, 0x1680),
    (0x180E, 0x180E),
    (0x2000, 0x200F),
    (0x2028, 0x202F),
    (0x205F, 0x2064),
    (0x2066, 0x206F),
    (0x3000, 0x3000),
    (0xE000, 0xF8FF),
    (0xFEFF, 0xFEFF),
    (0xFFF9, 0xFFFB),
    (0x110BD, 0x110BD),
    (0x110CD, 0x110CD),
    (0x13430, 0x1343F),
    (0x1BCA0, 0x1BCA3),
    (0x1D173, 0x1D17A),
    (0xE0001, 0xE0001),
    (0xE0020, 0xE007F),
    (0xF0000, 0xFFFFD),
    (0x100000, 0x10FFFD),
];

fn printable(c: char) -> bool {
    let code = c as u32;
    !UNPRINTABLE
        .iter()
        .any(|(first, last)| (*first..=*last).contains(&code))
}

/// `repr(s)` of a `str`: single quotes unless the text holds one and no double quote, the
/// quote and the backslash escaped, `\t`, `\n` and `\r` by name, every other control
/// character and every unprintable assigned code point as `\x`, `\u` or `\U` and lowercase
/// hex digits. A code point Unicode leaves unassigned is written as it is.
pub fn py_repr_str(s: &str) -> String {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            _ if c == quote => {
                out.push('\\');
                out.push(c);
            }
            ' '..='~' => out.push(c),
            _ if c > '\u{7f}' && printable(c) => out.push(c),
            _ => {
                let code = c as u32;
                if code <= 0xff {
                    write!(out, "\\x{code:02x}").unwrap();
                } else if code <= 0xffff {
                    write!(out, "\\u{code:04x}").unwrap();
                } else {
                    write!(out, "\\U{code:08x}").unwrap();
                }
            }
        }
    }
    out.push(quote);
    out
}

/// `repr(items)` of a list of `str`.
pub fn py_repr_strs<S: AsRef<str>>(items: &[S]) -> String {
    let reprs: Vec<String> = items
        .iter()
        .map(|item| py_repr_str(item.as_ref()))
        .collect();
    format!("[{}]", reprs.join(", "))
}

/// `repr(value)` of what `json.loads` returns: `None`, `True` and `False`, an integer's
/// digits, `float.__repr__` with `nan` and `inf`, a string as `py_repr_str` writes it, and a
/// list or a dict of those.
pub fn py_repr(value: &Value) -> String {
    match value {
        Value::Null => "None".to_string(),
        Value::Bool(true) => "True".to_string(),
        Value::Bool(false) => "False".to_string(),
        Value::Int(int) => int.text().to_string(),
        Value::Float(x) if x.is_nan() => "nan".to_string(),
        Value::Float(x) if x.is_infinite() => if *x > 0.0 { "inf" } else { "-inf" }.to_string(),
        Value::Float(x) => py_float_repr(*x),
        Value::Str(s) => py_repr_str(s),
        Value::List(items) => {
            let reprs: Vec<String> = items.iter().map(py_repr).collect();
            format!("[{}]", reprs.join(", "))
        }
        Value::Object(object) => {
            let reprs: Vec<String> = object
                .iter()
                .map(|(key, value)| format!("{}: {}", py_repr_str(key), py_repr(value)))
                .collect();
            format!("{{{}}}", reprs.join(", "))
        }
    }
}

/// `str(value)` of what `json.loads` returns: a string as it is, anything else its `repr`.
pub fn py_str(value: &Value) -> String {
    match value {
        Value::Str(s) => s.clone(),
        other => py_repr(other),
    }
}

/// `str(exc)` of the `UnicodeDecodeError` that `bytes.decode("utf-8")` raises, `None` for
/// bytes that are UTF-8: the first byte no sequence can start with, the first sequence a
/// byte breaks off, or the sequence the data ends in.
pub fn utf8_decode_error(bytes: &[u8]) -> Option<String> {
    let start = std::str::from_utf8(bytes).err()?.valid_up_to();
    let rest = &bytes[start..];
    let continuation = |i: usize| rest.get(i).is_some_and(|b| b & 0xC0 == 0x80);
    let lead = rest[0];
    // The count of bytes the error covers, zero for data that ends inside a sequence.
    let broken = match lead {
        0x80..=0xC1 | 0xF5..=0xFF => 1,
        0xC2..=0xDF => usize::from(rest.len() >= 2),
        0xE0..=0xEF => {
            let second_refused = rest.len() >= 2
                && (!continuation(1)
                    || (lead == 0xE0 && rest[1] < 0xA0)
                    || (lead == 0xED && rest[1] >= 0xA0));
            if second_refused {
                1
            } else if rest.len() >= 3 {
                2
            } else {
                0
            }
        }
        _ => {
            let second_refused = rest.len() >= 2
                && (!continuation(1)
                    || (lead == 0xF0 && rest[1] < 0x90)
                    || (lead == 0xF4 && rest[1] >= 0x90));
            if second_refused {
                1
            } else if rest.len() >= 3 && !continuation(2) {
                2
            } else if rest.len() >= 4 {
                3
            } else {
                0
            }
        }
    };
    let (end, reason) = match (broken, lead) {
        (0, _) => (bytes.len(), "unexpected end of data"),
        (1, 0x80..=0xC1 | 0xF5..=0xFF) => (start + 1, "invalid start byte"),
        (n, _) => (start + n, "invalid continuation byte"),
    };
    Some(if end == start + 1 {
        format!("'utf-8' codec can't decode byte 0x{lead:02x} in position {start}: {reason}")
    } else {
        format!(
            "'utf-8' codec can't decode bytes in position {start}-{}: {reason}",
            end - 1
        )
    })
}

/// `str(exc)` of the `OSError` a file operation on `path` raises:
/// `[Errno <n>] <strerror>: <repr(path)>`.
pub fn os_error_text(error: &std::io::Error, path: &str) -> String {
    match error.raw_os_error() {
        Some(code) => {
            let described = std::io::Error::from_raw_os_error(code).to_string();
            let strerror = described
                .strip_suffix(&format!(" (os error {code})"))
                .unwrap_or(&described);
            format!("[Errno {code}] {strerror}: {}", py_repr_str(path))
        }
        None => format!("{error}: {}", py_repr_str(path)),
    }
}

/// The text of a file as `Path(path).read_text(encoding="utf-8")` returns it, `\r\n` and a
/// lone `\r` read as `\n`; the error is `str(exc)` of the `OSError` or `UnicodeDecodeError`
/// the read raises.
pub fn read_text(path: &str) -> Result<String, String> {
    let bytes = std::fs::read(path).map_err(|error| os_error_text(&error, path))?;
    match String::from_utf8(bytes) {
        Ok(text) => Ok(universal_newlines(text)),
        Err(error) => Err(utf8_decode_error(error.as_bytes()).expect("bytes outside UTF-8")),
    }
}

/// Universal newlines on read: `\r\n` and a lone `\r` become `\n`.
pub fn universal_newlines(text: String) -> String {
    if !text.contains('\r') {
        return text;
    }
    text.replace("\r\n", "\n").replace('\r', "\n")
}

/// `text.splitlines()`: the lines between the boundaries `\n`, `\r`, `\r\n`, U+000B, U+000C,
/// U+001C to U+001E, U+0085, U+2028 and U+2029, with no empty line after a last boundary.
pub fn splitlines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = text.char_indices().peekable();
    while let Some((at, c)) = chars.next() {
        let boundary = matches!(
            c,
            '\n' | '\r' | '\u{0b}' | '\u{0c}' | '\u{1c}'
                ..='\u{1e}' | '\u{85}' | '\u{2028}' | '\u{2029}'
        );
        if !boundary {
            continue;
        }
        lines.push(&text[start..at]);
        let mut end = at + c.len_utf8();
        if c == '\r' && chars.peek().is_some_and(|(_, next)| *next == '\n') {
            chars.next();
            end += 1;
        }
        start = end;
    }
    if start < text.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// `int(text)` within 64 bits: the blanks `int()` drops around the number (those of
/// `str.isspace()` but U+001C to U+001F), an optional sign, then decimal digits with single
/// underscores between them. `None` where `int()` raises `ValueError`, for a digit outside
/// ASCII, which `int()` reads, and for a number 64 bits do not hold.
pub fn py_int(text: &str) -> Option<i64> {
    let body = text.trim_matches(|c| is_py_space(c) && !('\u{1c}'..='\u{1f}').contains(&c));
    let (negative, digits) = match body.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, body.strip_prefix('+').unwrap_or(body)),
    };
    if digits.is_empty()
        || digits.starts_with('_')
        || digits.ends_with('_')
        || digits.contains("__")
    {
        return None;
    }
    let mut magnitude = 0i128;
    for byte in digits.bytes() {
        match byte {
            b'_' => {}
            b'0'..=b'9' => {
                magnitude = magnitude * 10 + i128::from(byte - b'0');
                if magnitude > i128::from(i64::MAX) + 1 {
                    return None;
                }
            }
            _ => return None,
        }
    }
    i64::try_from(if negative { -magnitude } else { magnitude }).ok()
}

/// `str(pathlib.Path(arg))`: repeated separators and `.` components dropped, a trailing
/// separator dropped, exactly two leading slashes kept, `..` left alone, `.` for an empty
/// result.
pub fn path_str(arg: &str) -> String {
    String::from_utf8(py_path_str(arg.as_bytes())).expect("a path cut at its separators")
}

/// `str(pathlib.Path(base) / name)`: `name` alone when it is absolute, else the two joined.
pub fn path_join(base: &str, name: &str) -> String {
    if name.starts_with('/') {
        path_str(name)
    } else if base.ends_with('/') {
        path_str(&format!("{base}{name}"))
    } else {
        path_str(&format!("{base}/{name}"))
    }
}

/// `pathlib.Path(path).is_file()`: a regular file, or a link to one; false where the path
/// cannot be read.
pub fn is_file(path: &str) -> bool {
    std::fs::metadata(path).is_ok_and(|metadata| metadata.is_file())
}

/// `pathlib.Path(path).is_dir()`: a directory, or a link to one; false where the path
/// cannot be read.
pub fn is_dir(path: &str) -> bool {
    std::fs::metadata(path).is_ok_and(|metadata| metadata.is_dir())
}

/// `hashlib.sha256(data).hexdigest()`: the digest of `data` as 64 lowercase hexadecimal
/// digits, the form a corpus's `provenance.json` records a file under.
pub fn sha256_hex(data: &[u8]) -> String {
    let mut hex = String::with_capacity(64);
    for byte in Sha256::digest(data) {
        let _ = write!(hex, "{byte:02x}");
    }
    hex
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::loads;

    #[test]
    fn the_digests_of_the_standard_examples() {
        assert_eq!(
            sha256_hex(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(
            sha256_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
    }

    #[test]
    fn repr_of_a_string_follows_cpython() {
        let cases = [
            ("x", "'x'"),
            ("it's", "\"it's\""),
            ("say \"hi\"", "'say \"hi\"'"),
            ("both ' and \"", "'both \\' and \"'"),
            ("tab\there", "'tab\\there'"),
            ("nl\n", "'nl\\n'"),
            ("cr\r", "'cr\\r'"),
            ("bs\\", "'bs\\\\'"),
            ("\0\u{1f}\u{7f}", "'\\x00\\x1f\\x7f'"),
            ("\u{e9}", "'\u{e9}'"),
            ("\u{a0}", "'\\xa0'"),
            ("\u{ad}", "'\\xad'"),
            ("\u{85}", "'\\x85'"),
            ("\u{200b}", "'\\u200b'"),
            ("\u{2028}", "'\\u2028'"),
            ("\u{1F600}", "'\u{1F600}'"),
            ("\u{e000}", "'\\ue000'"),
            ("\u{e0001}", "'\\U000e0001'"),
            ("", "''"),
        ];
        for (text, expected) in cases {
            assert_eq!(py_repr_str(text), expected, "repr({text:?})");
        }
        assert_eq!(py_repr_strs(&["a", "b'c"]), "['a', \"b'c\"]");
        assert_eq!(py_repr_strs::<&str>(&[]), "[]");
    }

    #[test]
    fn repr_and_str_of_a_value_follow_cpython() {
        let value = loads(
            r#"{"a": 1, "b": [true, null, 1.5, "x", NaN, Infinity, -0.0, 1e16], "c": {"d": "e'f"}}"#,
        )
        .unwrap();
        assert_eq!(
            py_repr(&value),
            "{'a': 1, 'b': [True, None, 1.5, 'x', nan, inf, -0.0, 1e+16], 'c': {'d': \"e'f\"}}"
        );
        assert_eq!(py_str(&Value::Str("21:100".to_string())), "21:100");
        assert_eq!(py_str(&loads("[\"a\"]").unwrap()), "['a']");
        assert_eq!(py_str(&loads("-Infinity").unwrap()), "-inf");
    }

    #[test]
    fn a_decode_error_reads_as_cpython_words_it() {
        let cases: [(&[u8], &str); 25] = [
            (b"ab\xff", "byte 0xff in position 2: invalid start byte"),
            (b"\xc3", "byte 0xc3 in position 0: unexpected end of data"),
            (
                b"a\xc3(",
                "byte 0xc3 in position 1: invalid continuation byte",
            ),
            (b"\xe2\x82", "bytes in position 0-1: unexpected end of data"),
            (
                b"\xe2\x28\xa1",
                "byte 0xe2 in position 0: invalid continuation byte",
            ),
            (
                b"\xe2\x82\x28",
                "bytes in position 0-1: invalid continuation byte",
            ),
            (
                b"\xf0\x9f\x98",
                "bytes in position 0-2: unexpected end of data",
            ),
            (
                b"\xf0\x28\x8c\xbc",
                "byte 0xf0 in position 0: invalid continuation byte",
            ),
            (
                b"\xf0\x9f\x28\xbc",
                "bytes in position 0-1: invalid continuation byte",
            ),
            (
                b"\xf0\x9f\x98\x28",
                "bytes in position 0-2: invalid continuation byte",
            ),
            (
                b"\xed\xa0\x80",
                "byte 0xed in position 0: invalid continuation byte",
            ),
            (
                b"\xed\xa0",
                "byte 0xed in position 0: invalid continuation byte",
            ),
            (
                b"\xe0\x80\x80",
                "byte 0xe0 in position 0: invalid continuation byte",
            ),
            (
                b"\xe0\x9f",
                "byte 0xe0 in position 0: invalid continuation byte",
            ),
            (
                b"\xf4\x90\x80\x80",
                "byte 0xf4 in position 0: invalid continuation byte",
            ),
            (
                b"\xf4\x90",
                "byte 0xf4 in position 0: invalid continuation byte",
            ),
            (
                b"\xf0\x80",
                "byte 0xf0 in position 0: invalid continuation byte",
            ),
            (b"\xf5", "byte 0xf5 in position 0: invalid start byte"),
            (b"\xc0\x80", "byte 0xc0 in position 0: invalid start byte"),
            (b"\x80", "byte 0x80 in position 0: invalid start byte"),
            (b"ok\xe2", "byte 0xe2 in position 2: unexpected end of data"),
            (b"\xf0\x9f", "bytes in position 0-1: unexpected end of data"),
            (b"\xf0", "byte 0xf0 in position 0: unexpected end of data"),
            (b"\xed\x9f", "bytes in position 0-1: unexpected end of data"),
            (b"\xee\x80", "bytes in position 0-1: unexpected end of data"),
        ];
        for (bytes, expected) in cases {
            assert_eq!(
                utf8_decode_error(bytes),
                Some(format!("'utf-8' codec can't decode {expected}")),
                "{bytes:?}"
            );
        }
        assert_eq!(utf8_decode_error("caf\u{e9}".as_bytes()), None);
    }

    #[test]
    fn lines_split_where_splitlines_splits() {
        assert_eq!(
            splitlines(
                "a\rb\r\nc\n\nd\u{0b}e\u{0c}f\u{1c}g\u{1d}h\u{1e}i\u{85}j\u{2028}k\u{2029}l\u{1f}m"
            ),
            ["a", "b", "c", "", "d", "e", "f", "g", "h", "i", "j", "k", "l\u{1f}m"]
        );
        assert_eq!(splitlines(""), Vec::<&str>::new());
        assert_eq!(splitlines("\n"), [""]);
        assert_eq!(splitlines("a\n"), ["a"]);
        assert_eq!(universal_newlines("a\r\nb\rc\n".to_string()), "a\nb\nc\n");
    }

    #[test]
    fn a_count_is_read_as_int_reads_it_within_64_bits() {
        let read = [
            ("20", 20),
            ("-0", 0),
            ("+2", 2),
            ("007", 7),
            ("1_000", 1000),
            (" 3 ", 3),
            ("\t3\n", 3),
            ("\u{0b}3\u{0c}", 3),
            ("\u{a0}3", 3),
            ("\u{85}3", 3),
            ("\u{2028}3", 3),
            ("9223372036854775807", i64::MAX),
            ("-9223372036854775808", i64::MIN),
        ];
        for (text, expected) in read {
            assert_eq!(py_int(text), Some(expected), "int({text:?})");
        }
        let refused = [
            "",
            " ",
            "x",
            "1.5",
            "0x10",
            "_1",
            "1_",
            "1__0",
            "- 3",
            "--3",
            "\u{1c}3",
            "3\u{1f}",
            "\u{663}",
            "9223372036854775808",
            "-9223372036854775809",
        ];
        for text in refused {
            assert_eq!(py_int(text), None, "int({text:?})");
        }
    }

    #[test]
    fn paths_join_as_pathlib_joins_them() {
        assert_eq!(path_str(""), ".");
        assert_eq!(path_str("corpus/"), "corpus");
        assert_eq!(path_join(".", "manifest.json"), "manifest.json");
        assert_eq!(path_join("corpus", "manifest.json"), "corpus/manifest.json");
        assert_eq!(path_join("/", "manifest.json"), "/manifest.json");
        assert_eq!(path_join("a/b", "../c"), "a/b/../c");
        assert_eq!(path_join("a", "/abs/file"), "/abs/file");
    }

    #[test]
    fn an_os_error_names_its_number_its_text_and_the_path() {
        let error = std::fs::read("/nonexistent/x").unwrap_err();
        assert_eq!(
            os_error_text(&error, "/nonexistent/x"),
            "[Errno 2] No such file or directory: '/nonexistent/x'"
        );
        assert_eq!(
            read_text("/nonexistent/it's"),
            Err("[Errno 2] No such file or directory: \"/nonexistent/it's\"".to_string())
        );
    }
}
