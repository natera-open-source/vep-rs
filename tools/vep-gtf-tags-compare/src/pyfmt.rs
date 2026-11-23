// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Text the tool reads and writes character for character as CPython does: `str.strip()`,
//! `int()`, `float.__repr__`, `json.dumps(ensure_ascii=True)` string escapes, `f"{n:,}"`
//! digit grouping and the order `sorted()` gives file names.

use std::cmp::Ordering;
use std::fmt::Write;

/// `str.isspace()`: the Unicode White_Space set plus U+001C to U+001F, which Python counts
/// through their bidirectional class and `char::is_whitespace` does not.
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

/// The integer digits `int()` converts in one call.
const INT_MAX_STR_DIGITS: usize = 4300;

/// `int(text)`, saturated to an `i64`: surrounding whitespace, one sign, ASCII digits with
/// single underscores between them. `None` where `int()` raises: U+001C to U+001F are not
/// whitespace to it, and more than 4,300 digits are refused. A digit outside ASCII, which
/// `int()` converts, is refused as well.
pub fn py_int(text: &str) -> Option<i64> {
    let is_int_space = |c: char| is_py_space(c) && !('\u{1c}'..='\u{1f}').contains(&c);
    let body = text.trim_matches(is_int_space);
    let (negative, digits) = match body.as_bytes().first()? {
        b'-' => (true, &body[1..]),
        b'+' => (false, &body[1..]),
        _ => (false, body),
    };
    let bytes = digits.as_bytes();
    if bytes.first().is_none_or(|b| *b == b'_') || bytes.last() == Some(&b'_') {
        return None;
    }
    let mut magnitude = 0i64;
    let mut count = 0usize;
    let mut after_underscore = false;
    for &b in bytes {
        match b {
            b'0'..=b'9' => {
                magnitude = magnitude
                    .saturating_mul(10)
                    .saturating_add(i64::from(b - b'0'));
                count += 1;
                after_underscore = false;
            }
            b'_' if !after_underscore => after_underscore = true,
            _ => return None,
        }
    }
    if count > INT_MAX_STR_DIGITS {
        return None;
    }
    Some(if negative { -magnitude } else { magnitude })
}

/// `repr(x)` for an f64: the shortest round-trip digits, `.0` on an integral value, exponent
/// form below 1e-4 and from 1e16 with a signed two-digit exponent, `nan`, `inf` and `-inf`.
pub fn py_float_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_sign_negative() {
        return format!("-{}", py_float_repr(-x));
    }
    if x.is_infinite() {
        return "inf".to_string();
    }
    let exp_form = format!("{x:e}");
    let (mantissa, exponent) = exp_form.split_once('e').expect("LowerExp output");
    let exponent: i32 = exponent.parse().expect("LowerExp exponent");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let decpt = exponent + 1;
    if decpt <= -4 || decpt > 16 {
        // LowerExp's mantissa is already Python's `d.ddd`; only the exponent's sign and
        // width differ.
        let sign = if exponent < 0 { '-' } else { '+' };
        let magnitude = exponent.abs();
        format!("{mantissa}e{sign}{magnitude:02}")
    } else if decpt <= 0 {
        format!("0.{}{}", "0".repeat(decpt.unsigned_abs() as usize), digits)
    } else if decpt as usize >= digits.len() {
        format!("{}{}.0", digits, "0".repeat(decpt as usize - digits.len()))
    } else {
        format!(
            "{}.{}",
            &digits[..decpt as usize],
            &digits[decpt as usize..]
        )
    }
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

/// `f"{n:,}"`.
pub fn py_grouped(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// The order `sorted()` gives two file names: by code point of the name as Python decodes
/// it, each byte outside UTF-8 standing as the lone surrogate U+DC00 plus the byte. Two
/// UTF-8 names compare as their bytes do.
pub fn py_name_order(a: &[u8], b: &[u8]) -> Ordering {
    match (std::str::from_utf8(a), std::str::from_utf8(b)) {
        (Ok(_), Ok(_)) => a.cmp(b),
        _ => surrogate_escaped(a).cmp(&surrogate_escaped(b)),
    }
}

fn surrogate_escaped(mut name: &[u8]) -> Vec<u32> {
    let mut out = Vec::with_capacity(name.len());
    loop {
        match std::str::from_utf8(name) {
            Ok(text) => {
                out.extend(text.chars().map(u32::from));
                return out;
            }
            Err(error) => {
                let (valid, rest) = name.split_at(error.valid_up_to());
                let text = std::str::from_utf8(valid).expect("the valid prefix");
                out.extend(text.chars().map(u32::from));
                let invalid = error.error_len().unwrap_or(rest.len());
                out.extend(rest[..invalid].iter().map(|b| 0xdc00 + u32::from(*b)));
                name = &rest[invalid..];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
            .map(|c| c as u32)
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(py_strip("\u{1f}tag\u{a0}\u{3000}"), "tag");
        assert_eq!(py_strip(" \t x \u{200b}"), "x \u{200b}");
    }

    #[test]
    fn int_reads_what_pythons_int_reads() {
        let cases: [(&str, Option<i64>); 24] = [
            ("5", Some(5)),
            (" 5 ", Some(5)),
            ("\t5\n", Some(5)),
            ("\u{0b}5\u{0c}", Some(5)),
            ("\u{a0}5", Some(5)),
            ("\u{3000}5\u{2003}", Some(5)),
            ("+5", Some(5)),
            ("-5", Some(-5)),
            ("-0", Some(0)),
            ("007", Some(7)),
            ("1_000", Some(1000)),
            ("- 5", None),
            ("+-5", None),
            ("_1", None),
            ("1_", None),
            ("1__0", None),
            ("", None),
            (" ", None),
            ("0x10", None),
            ("5.0", None),
            ("1e3", None),
            ("\u{1c}5", None),
            ("5\u{1c}", None),
            ("5\0", None),
        ];
        for (text, expected) in cases {
            assert_eq!(py_int(text), expected, "int({text:?})");
        }
    }

    #[test]
    fn int_saturates_and_stops_at_4300_digits() {
        assert_eq!(py_int(&"9".repeat(4300)), Some(i64::MAX));
        assert_eq!(py_int(&format!("-{}", "9".repeat(4300))), Some(-i64::MAX));
        assert_eq!(py_int(&"9".repeat(4301)), None);
        assert_eq!(py_int(&"0".repeat(4301)), None);
        assert_eq!(py_int(&format!("{}1", "1_".repeat(2150))), Some(i64::MAX));
    }

    #[test]
    fn float_repr_matches_python() {
        let cases = [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (1.5, "1.5"),
            (0.5, "0.5"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (0.000015, "1.5e-05"),
            (1e15, "1000000000000000.0"),
            (1e16, "1e+16"),
            (1e22, "1e+22"),
            (123456789012345680.0, "1.2345678901234568e+17"),
            (100000.0, "100000.0"),
            (-0.01, "-0.01"),
            (0.1 + 0.2, "0.30000000000000004"),
            (f64::NAN, "nan"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
        ];
        for (x, expected) in cases {
            assert_eq!(py_float_repr(x), expected, "repr({x:?})");
        }
    }

    #[test]
    fn json_string_escapes_match_ensure_ascii() {
        let mut out = String::new();
        push_json_str(
            &mut out,
            "\u{1f}\u{7f}\u{e9}\u{1F600}\u{fffd}/\\\"\n\t\u{08}\u{0c}\r\0",
        );
        let escape = |hex: &str| format!("{}u{hex}", '\\');
        let expected = format!(
            r#""{}{}{}{}{}{}/\\\"\n\t\b\f\r{}""#,
            escape("001f"),
            escape("007f"),
            escape("00e9"),
            escape("d83d"),
            escape("de00"),
            escape("fffd"),
            escape("0000"),
        );
        assert_eq!(out, expected);
    }

    #[test]
    fn grouping_matches_format_comma() {
        assert_eq!(py_grouped(0), "0");
        assert_eq!(py_grouped(999), "999");
        assert_eq!(py_grouped(1000), "1,000");
        assert_eq!(py_grouped(1234567), "1,234,567");
        assert_eq!(py_grouped(125_572_619), "125,572,619");
    }

    #[test]
    fn names_sort_by_code_point_with_stray_bytes_as_lone_surrogates() {
        assert_eq!(py_name_order(b"1-2.json", b"10-11.json"), Ordering::Less);
        assert_eq!(py_name_order(b"a", b"a"), Ordering::Equal);
        assert_eq!(
            py_name_order("\u{e9}".as_bytes(), "\u{4e2d}".as_bytes()),
            Ordering::Less
        );
        // U+DCFF, the stray byte 0xff, sorts below U+E000 although 0xff is the larger byte.
        assert_eq!(
            py_name_order(b"\xff", "\u{e000}".as_bytes()),
            Ordering::Less
        );
        assert_eq!(
            py_name_order(b"\xff", "\u{d7ff}".as_bytes()),
            Ordering::Greater
        );
        assert_eq!(
            surrogate_escaped(b"a\xe2\x82z\xc3\xa9\x80"),
            [0x61, 0xdce2, 0xdc82, 0x7a, 0xe9, 0xdc80]
        );
    }
}
