// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! CPython's text rules the comparison's report depends on: `str.strip()`, `float()` and
//! `int()` over text, `float.__repr__`, `repr()` of a `str`, `json.dumps` string escapes,
//! `urllib.parse.unquote` and the name and stem of a `pathlib` path.

use regex::Regex;
use std::borrow::Cow;
use std::fmt::Write;
use std::sync::OnceLock;

/// `str.isspace()`: the Unicode White_Space set plus U+001C to U+001F.
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

/// The code point of the digit zero of every run of ten decimal digits (Unicode 16).
const DECIMAL_ZEROS: [u32; 76] = [
    0x30, 0x660, 0x6f0, 0x7c0, 0x966, 0x9e6, 0xa66, 0xae6, 0xb66, 0xbe6, 0xc66, 0xce6, 0xd66,
    0xde6, 0xe50, 0xed0, 0xf20, 0x1040, 0x1090, 0x17e0, 0x1810, 0x1946, 0x19d0, 0x1a80, 0x1a90,
    0x1b50, 0x1bb0, 0x1c40, 0x1c50, 0xa620, 0xa8d0, 0xa900, 0xa9d0, 0xa9f0, 0xaa50, 0xabf0, 0xff10,
    0x104a0, 0x10d30, 0x10d40, 0x11066, 0x110f0, 0x11136, 0x111d0, 0x112f0, 0x11450, 0x114d0,
    0x11650, 0x116c0, 0x116d0, 0x116da, 0x11730, 0x118e0, 0x11950, 0x11bf0, 0x11c50, 0x11d50,
    0x11da0, 0x11f50, 0x16130, 0x16a60, 0x16ac0, 0x16b50, 0x16d70, 0x1ccf0, 0x1d7ce, 0x1d7d8,
    0x1d7e2, 0x1d7ec, 0x1d7f6, 0x1e140, 0x1e2f0, 0x1e4f0, 0x1e5f1, 0x1e950, 0x1fbf0,
];

/// The value of a decimal digit of any script.
pub fn decimal_digit(c: char) -> Option<u8> {
    let code = c as u32;
    let run = DECIMAL_ZEROS.partition_point(|zero| *zero <= code);
    let zero = *DECIMAL_ZEROS.get(run.checked_sub(1)?)?;
    (code - zero < 10).then(|| (code - zero) as u8)
}

/// The ASCII text `float()` and `int()` parse: a decimal digit of any script as its ASCII
/// digit and a space outside ASCII as a blank, the ends trimmed of blanks and of the ASCII
/// white space (U+001C to U+001F are not trimmed), an underscore dropped where it stands
/// between two digits. `None` where a character outside ASCII is neither, where an
/// underscore stands anywhere else, and for a text with nothing but white space.
fn number_text(text: &str) -> Option<String> {
    let mut ascii = String::with_capacity(text.len());
    for c in text.chars() {
        if (c as u32) < 127 {
            ascii.push(c);
        } else if is_py_space(c) {
            ascii.push(' ');
        } else {
            ascii.push((b'0' + decimal_digit(c)?) as char);
        }
    }
    let blank = |c: char| matches!(c, ' ' | '\t' | '\n' | '\r' | '\u{0b}' | '\u{0c}');
    let trimmed = ascii.trim_matches(blank);
    if trimmed.is_empty() {
        return None;
    }
    if !ascii.contains('_') {
        return Some(trimmed.to_string());
    }
    // The underscores are judged over the text before its ends are trimmed.
    let mut without = String::with_capacity(ascii.len());
    let mut previous = '\0';
    for c in ascii.chars() {
        if c == '_' {
            if !previous.is_ascii_digit() {
                return None;
            }
        } else {
            if previous == '_' && !c.is_ascii_digit() {
                return None;
            }
            without.push(c);
        }
        previous = c;
    }
    if previous == '_' {
        return None;
    }
    Some(without.trim_matches(blank).to_string())
}

/// `float(text)`, `None` where CPython raises `ValueError`.
pub fn py_float(text: &str) -> Option<f64> {
    let number = number_text(text)?;
    let bytes = number.as_bytes();
    // Rust's grammar is CPython's once the text is ASCII without blanks and underscores:
    // a sign, then `inf`, `infinity` or `nan` in any case, or digits with an optional point
    // and exponent. Anything else Rust accepts is excluded by the characters allowed here.
    let allowed = bytes.iter().all(|byte| {
        byte.is_ascii_digit() || matches!(byte, b'+' | b'-' | b'.') || byte.is_ascii_alphabetic()
    });
    if !allowed {
        return None;
    }
    number.parse::<f64>().ok()
}

/// An `int` of any size: its sign and its decimal digits without leading zeros.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PyInt {
    negative: bool,
    digits: String,
}

impl PyInt {
    pub fn from_i64(value: i64) -> PyInt {
        PyInt {
            negative: value < 0,
            digits: value.unsigned_abs().to_string(),
        }
    }

    pub fn is_negative(&self) -> bool {
        self.negative
    }

    /// The value where it fits, the nearest `i128` where it does not.
    pub fn saturating_i128(&self) -> i128 {
        match self.digits.parse::<i128>() {
            Ok(magnitude) if self.negative => -magnitude,
            Ok(magnitude) => magnitude,
            Err(_) if self.negative => i128::MIN,
            Err(_) => i128::MAX,
        }
    }
}

impl std::fmt::Display for PyInt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.negative {
            f.write_str("-")?;
        }
        f.write_str(&self.digits)
    }
}

/// The most decimal digits `int()` converts from text.
const INT_MAX_STR_DIGITS: usize = 4300;

/// `int(text)`, `None` where CPython raises `ValueError`.
pub fn py_int(text: &str) -> Option<PyInt> {
    let number = number_text(text)?;
    let (negative, digits) = match number.as_bytes()[0] {
        b'-' => (true, &number[1..]),
        b'+' => (false, &number[1..]),
        _ => (false, &number[..]),
    };
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    if digits.len() > INT_MAX_STR_DIGITS {
        return None;
    }
    let significant = digits.trim_start_matches('0');
    Some(if significant.is_empty() {
        PyInt {
            negative: false,
            digits: "0".to_string(),
        }
    } else {
        PyInt {
            negative,
            digits: significant.to_string(),
        }
    })
}

pub use crate::pytext::float_repr as py_float_repr;

/// `json.dumps(x)` of a float: its `repr`, and `NaN`, `Infinity` and `-Infinity`.
pub fn json_float(x: f64) -> String {
    if x.is_nan() {
        "NaN".to_string()
    } else if x.is_infinite() {
        if x > 0.0 { "Infinity" } else { "-Infinity" }.to_string()
    } else {
        py_float_repr(x)
    }
}

/// `f"{x:.2f}"`: two decimals of the binary value, a tie rounded to the even digit.
pub fn fixed_two(x: f64) -> String {
    format!("{x:.2}")
}

/// Appends `json.dumps(s)` with `ensure_ascii=True`: the seven short escapes, every other
/// code point outside `0x20..=0x7e` as `\u` and four lowercase hex digits, a surrogate pair
/// above U+FFFF.
pub fn push_json_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            ' '..='~' => out.push(c),
            _ => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units) {
                    write!(out, "\\u{unit:04x}").unwrap();
                }
            }
        }
    }
    out.push('"');
}

/// `str.isprintable()` of one character outside ASCII: false for the Other and Separator
/// general categories.
fn printable(c: char) -> bool {
    static UNPRINTABLE: OnceLock<Regex> = OnceLock::new();
    let unprintable = UNPRINTABLE
        .get_or_init(|| Regex::new(r"^[\p{C}\p{Z}]$").expect("the category class compiles"));
    !unprintable.is_match(c.encode_utf8(&mut [0; 4]))
}

/// Appends `repr(s)` of a `str`: single quotes unless the text holds one and no double
/// quote, the quote and the backslash escaped, `\t`, `\n` and `\r` by name, any other
/// unprintable character as `\x`, `\u` or `\U` with two, four or eight hex digits.
pub fn push_str_repr(out: &mut String, s: &str) {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
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
            _ if (c as u32) < 0x100 && ((c as u32) < 0x80 || !printable(c)) => {
                write!(out, "\\x{:02x}", c as u32).unwrap();
            }
            _ if printable(c) => out.push(c),
            _ if (c as u32) < 0x10000 => write!(out, "\\u{:04x}", c as u32).unwrap(),
            _ => write!(out, "\\U{:08x}", c as u32).unwrap(),
        }
    }
    out.push(quote);
}

/// `urllib.parse.unquote(s)`: within each run of ASCII characters every `%` and two hex
/// digits is the byte they name, and the run's bytes are read as UTF-8 with replacement;
/// characters outside ASCII stand as they are, and a `%` without two hex digits stays.
pub fn unquote(s: &str) -> Cow<'_, str> {
    if !s.contains('%') {
        return Cow::Borrowed(s);
    }
    let hex = |byte: u8| (byte as char).to_digit(16).map(|digit| digit as u8);
    let mut out = String::with_capacity(s.len());
    let mut run: Vec<u8> = Vec::new();
    let bytes = s.as_bytes();
    let mut at = 0;
    while at < bytes.len() {
        let byte = bytes[at];
        if byte >= 0x80 {
            // A character outside ASCII ends the run before it.
            if !run.is_empty() {
                out.push_str(&String::from_utf8_lossy(&run));
                run.clear();
            }
            let width = s[at..].chars().next().map_or(1, char::len_utf8);
            out.push_str(&s[at..at + width]);
            at += width;
            continue;
        }
        if byte == b'%' && at + 2 < bytes.len() {
            if let (Some(high), Some(low)) = (hex(bytes[at + 1]), hex(bytes[at + 2])) {
                run.push(high << 4 | low);
                at += 3;
                continue;
            }
        }
        run.push(byte);
        at += 1;
    }
    if !run.is_empty() {
        out.push_str(&String::from_utf8_lossy(&run));
    }
    Cow::Owned(out)
}

/// `str(pathlib.PurePosixPath(arg))`: repeated separators and `.` components dropped, a
/// trailing separator dropped, exactly two leading slashes kept, `.` for an empty result.
pub fn path_str(arg: &str) -> String {
    let leading_slashes = arg.bytes().take_while(|byte| *byte == b'/').count();
    let root = match leading_slashes {
        0 => "",
        2 => "//",
        _ => "/",
    };
    let parts: Vec<&str> = arg
        .split('/')
        .filter(|part| !part.is_empty() && *part != ".")
        .collect();
    let joined = format!("{root}{}", parts.join("/"));
    if joined.is_empty() {
        ".".to_string()
    } else {
        joined
    }
}

/// `pathlib.PurePosixPath(arg).name`: the last component that is neither empty nor `.`,
/// empty where there is none.
pub fn path_name(arg: &str) -> &str {
    arg.split('/')
        .rfind(|part| !part.is_empty() && *part != ".")
        .unwrap_or("")
}

/// `pathlib.PurePosixPath(arg).stem`: the name up to its last dot, where what stands before
/// that dot holds a character other than a dot; the whole name otherwise.
pub fn path_stem(arg: &str) -> &str {
    let name = path_name(arg);
    match name.rfind('.') {
        Some(dot) if !name[..dot].trim_start_matches('.').is_empty() => &name[..dot],
        _ => name,
    }
}

/// `str(pathlib.PurePosixPath(arg).parent)`.
pub fn path_parent(arg: &str) -> String {
    let whole = path_str(arg);
    let root_len = whole.bytes().take_while(|byte| *byte == b'/').count();
    let (root, tail) = whole.split_at(root_len);
    if tail.is_empty() || tail == "." {
        return whole;
    }
    match tail.rfind('/') {
        Some(slash) => format!("{root}{}", &tail[..slash]),
        None if root.is_empty() => ".".to_string(),
        None => root.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_reads_text_as_cpython_does() {
        let read = [
            ("\u{a0}1", 1.0),
            ("\u{3000}1\u{2003}", 1.0),
            ("1\u{0b}", 1.0),
            ("\u{85}1", 1.0),
            ("1_0", 10.0),
            (" 1_0 ", 10.0),
            ("1e1_0", 1e10),
            ("\u{661}\u{662}\u{663}", 123.0),
            ("\u{ff11}\u{ff12}", 12.0),
            ("1.", 1.0),
            (".5", 0.5),
            ("\t1\n", 1.0),
            ("1.5e3", 1500.0),
            ("+.5e-2", 0.005),
            ("1E2", 100.0),
            ("inf", f64::INFINITY),
            ("-InFinity", f64::NEG_INFINITY),
            ("1e999", f64::INFINITY),
            ("-1e999", f64::NEG_INFINITY),
            ("0.30000000000000004", 0.1 + 0.2),
        ];
        for (text, expected) in read {
            assert_eq!(py_float(text), Some(expected), "{text:?}");
        }
        for text in ["nan", "+nan", "NAN", "-nan"] {
            assert!(py_float(text).is_some_and(f64::is_nan), "{text:?}");
        }
        let refused = [
            "1\u{1c}", "\u{1c}1", "1__0", "_1", "1_", "1_.5", "1._5", "1e_1", ".", "1e", "1e+",
            "infi", "0x10", "", " ", "1 2", "1\0", "\u{b2}", "1,5", "- 1", "+-1", "1\u{a0}2", "e5",
            "in", "1f", "1d5", "--1", "1.2.3", "infinit", "nane",
        ];
        for text in refused {
            assert_eq!(py_float(text), None, "{text:?}");
        }
    }

    #[test]
    fn int_reads_text_as_cpython_does() {
        let read = [
            ("\u{a0}1", "1"),
            ("1\u{0b}", "1"),
            ("1_0", "10"),
            (" 1_0 ", "10"),
            ("\u{661}\u{662}\u{663}", "123"),
            ("\t1\n", "1"),
            ("-1", "-1"),
            ("+7", "7"),
            ("-0", "0"),
            ("007", "7"),
            ("-000", "0"),
            (
                "340282366920938463463374607431768211456",
                "340282366920938463463374607431768211456",
            ),
        ];
        for (text, expected) in read {
            assert_eq!(
                py_int(text).map(|n| n.to_string()).as_deref(),
                Some(expected),
                "{text:?}"
            );
        }
        let refused = [
            "1\u{1c}", "1__0", "_1", "1_", "1e1_0", "1.", ".5", "inf", "0x10", "", " ", "1 2",
            "1\0", "\u{b2}", "1.5", "- 1", "+-1", "+", "-",
        ];
        for text in refused {
            assert_eq!(py_int(text), None, "{text:?}");
        }
        assert!(py_int(&"9".repeat(4300)).is_some());
        assert!(py_int(&"9".repeat(4301)).is_none());
        assert!(py_int("-5").unwrap().is_negative());
        assert_eq!(py_int("-5").unwrap().saturating_i128(), -5);
        assert_eq!(
            py_int(&"9".repeat(60)).unwrap().saturating_i128(),
            i128::MAX
        );
        assert_eq!(PyInt::from_i64(-3).to_string(), "-3");
    }

    #[test]
    fn float_repr_matches_cpython() {
        let cases = [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (0.5, "0.5"),
            (0.0001, "0.0001"),
            (0.00005, "5e-05"),
            (0.000001, "1e-06"),
            (1e15, "1000000000000000.0"),
            (1e16, "1e+16"),
            (123.456, "123.456"),
            (0.1 + 0.2, "0.30000000000000004"),
            (100.0, "100.0"),
            (1.0 / 3.0 * 100.0, "33.33333333333333"),
            (2.0 / 3.0 * 100.0, "66.66666666666666"),
            (1e-300, "1e-300"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
        ];
        for (x, expected) in cases {
            assert_eq!(py_float_repr(x), expected, "repr({x:?})");
        }
        assert_eq!(py_float_repr(f64::NAN), "nan");
        assert_eq!(json_float(f64::NAN), "NaN");
        assert_eq!(json_float(f64::INFINITY), "Infinity");
        assert_eq!(json_float(f64::NEG_INFINITY), "-Infinity");
        assert_eq!(json_float(1e-6), "1e-06");
        // 26,217 of 262,144 records as a percentage lies half way between two texts.
        assert_eq!(json_float(26217.0 / 262144.0 * 100.0), "10.000991821289062");
    }

    #[test]
    fn two_decimals_round_a_tie_of_the_binary_value_to_even() {
        let cases = [
            (0.125, "0.12"),
            (0.375, "0.38"),
            (2.675, "2.67"),
            (1.005, "1.00"),
            (0.005, "0.01"),
            (0.015, "0.01"),
            (0.025, "0.03"),
            (99.995, "100.00"),
            (99.994999, "99.99"),
            (1.0 / 800.0 * 100.0, "0.12"),
            (3.0 / 800.0 * 100.0, "0.38"),
            (1.0 / 3.0 * 100.0, "33.33"),
            (2.0 / 3.0 * 100.0, "66.67"),
            (100.0, "100.00"),
            (0.0, "0.00"),
            (12.5, "12.50"),
        ];
        for (x, expected) in cases {
            assert_eq!(fixed_two(x), expected, "{x:?}");
        }
    }

    #[test]
    fn str_repr_matches_cpython() {
        let cases = [
            ("0.100", "'0.100'"),
            ("", "''"),
            ("it's", "\"it's\""),
            ("say \"hi\"", "'say \"hi\"'"),
            ("it's \"both\"", "'it\\'s \"both\"'"),
            ("a\\b", "'a\\\\b'"),
            ("a\tb\nc\rd", "'a\\tb\\nc\\rd'"),
            ("\u{0}\u{1f}\u{7f}", "'\\x00\\x1f\\x7f'"),
            ("\u{80}\u{a0}\u{ad}", "'\\x80\\xa0\\xad'"),
            ("\u{e9}\u{fffd}\u{4e2d}", "'\u{e9}\u{fffd}\u{4e2d}'"),
            ("\u{2028}\u{200b}\u{3000}", "'\\u2028\\u200b\\u3000'"),
            ("\u{1f600}", "'\u{1f600}'"),
            ("\u{e0001}\u{10ffff}", "'\\U000e0001\\U0010ffff'"),
            ("\u{378}", "'\\u0378'"),
        ];
        for (text, expected) in cases {
            let mut out = String::new();
            push_str_repr(&mut out, text);
            assert_eq!(out, expected, "{text:?}");
        }
    }

    #[test]
    fn json_string_escapes_match_ensure_ascii() {
        let mut out = String::new();
        push_json_str(
            &mut out,
            "\u{1f}\u{7f}\u{e9}\u{1F600}\u{fffd}/\\\"\n\t\u{08}\u{0c}\r\0",
        );
        assert_eq!(
            out,
            r#""\u001f\u007f\u00e9\ud83d\ude00\ufffd/\\\"\n\t\b\f\r\u0000""#
        );
    }

    #[test]
    fn unquote_matches_urllib() {
        let cases = [
            ("abc", "abc"),
            ("abc%20def", "abc def"),
            ("%3B%3d%7C%2c", ";=|,"),
            ("100%", "100%"),
            ("%", "%"),
            ("%2", "%2"),
            ("%zz", "%zz"),
            ("%%41", "%A"),
            ("%4%41", "%4A"),
            ("%C3%A9", "\u{e9}"),
            ("%c3%a9t%C3%A9", "\u{e9}t\u{e9}"),
            ("%C3", "\u{fffd}"),
            ("%C3\u{e9}%A9", "\u{fffd}\u{e9}\u{fffd}"),
            ("%e2%82", "\u{fffd}"),
            ("%e2%82%ac", "\u{20ac}"),
            ("%FF%FE", "\u{fffd}\u{fffd}"),
            ("%ED%A0%80", "\u{fffd}\u{fffd}\u{fffd}"),
            ("%F0%9F%98%80", "\u{1f600}"),
            ("%F0%9F%98", "\u{fffd}"),
            ("%C0%AF", "\u{fffd}\u{fffd}"),
            ("%F4%90%80%80", "\u{fffd}\u{fffd}\u{fffd}\u{fffd}"),
            ("\u{e9}%41\u{e9}", "\u{e9}A\u{e9}"),
            ("%00%09", "\0\t"),
            ("a%2", "a%2"),
            ("a%2g", "a%2g"),
            ("%41%", "A%"),
        ];
        for (text, expected) in cases {
            assert_eq!(unquote(text), expected, "{text:?}");
        }
        // The hexadecimal digits are read in either case.
        assert_eq!(unquote(&"%e2%82%ac".to_uppercase()), "\u{20ac}");
    }

    #[test]
    fn path_parts_match_pathlib() {
        let cases = [
            ("a/b.txt", "b.txt", "b", "a/b.txt", "a"),
            ("a/b/", "b", "b", "a/b", "a"),
            ("/", "", "", "/", "/"),
            (".", "", "", ".", "."),
            ("..", "..", "..", "..", "."),
            ("a/.", "a", "a", "a", "."),
            ("", "", "", ".", "."),
            ("foo.", "foo.", "foo", "foo.", "."),
            ("..a", "..a", "..a", "..a", "."),
            (".bashrc", ".bashrc", ".bashrc", ".bashrc", "."),
            ("a.b.c", "a.b.c", "a.b", "a.b.c", "."),
            ("x.tar.gz", "x.tar.gz", "x.tar", "x.tar.gz", "."),
            ("./perl.txt", "perl.txt", "perl", "perl.txt", "."),
            ("a//b//c.d", "c.d", "c", "a/b/c.d", "a/b"),
            ("...", "...", "...", "...", "."),
            ("a/..", "..", "..", "a/..", "a"),
            (".x.", ".x.", ".x", ".x.", "."),
            ("x..", "x..", "x.", "x..", "."),
            ("x. ", "x. ", "x", "x. ", "."),
            ("/x", "x", "x", "/x", "/"),
            ("//x/y", "y", "y", "//x/y", "//x"),
            ("///x", "x", "x", "/x", "/"),
        ];
        for (arg, name, stem, whole, parent) in cases {
            assert_eq!(path_name(arg), name, "name of {arg:?}");
            assert_eq!(path_stem(arg), stem, "stem of {arg:?}");
            assert_eq!(path_str(arg), whole, "str of {arg:?}");
            assert_eq!(path_parent(arg), parent, "parent of {arg:?}");
        }
    }

    #[test]
    fn strip_set_is_cpythons_isspace() {
        assert_eq!(py_strip("\u{1f}A/T\u{a0}\u{3000}"), "A/T");
        assert_eq!(py_strip(" \t x \u{200b}"), "x \u{200b}");
        let count = (0..=0x10ffffu32)
            .filter_map(char::from_u32)
            .filter(|c| is_py_space(*c))
            .count();
        assert_eq!(count, 29);
    }
}
