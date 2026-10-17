// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Text the report must reproduce character for character from CPython: `str.strip()`,
//! `float.__repr__`, `round(x, 6)`, `json.dumps(ensure_ascii=True)` string escapes,
//! `f"{n:,}"` digit grouping and `str(pathlib.Path(arg))`.

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

/// `repr(x)` for a finite f64: the shortest round-trip digits, `.0` on an integral value,
/// exponent form below 1e-4 and from 1e16 with a signed two-digit exponent.
pub fn py_float_repr(x: f64) -> String {
    if x.is_sign_negative() {
        return format!("-{}", py_float_repr(-x));
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

/// `round(x, 6)`: the f64 nearest to the half-to-even six-decimal rounding of x's exact
/// binary value.
pub fn py_round6(x: f64) -> f64 {
    format!("{x:.6}").parse().expect("fixed-point float text")
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
pub fn py_grouped(n: u64) -> String {
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
        assert_eq!(py_strip("\u{1f}A/T\u{a0}\u{3000}"), "A/T");
        assert_eq!(py_strip(" \t x \u{200b}"), "x \u{200b}");
    }

    #[test]
    fn float_repr_matches_python() {
        let cases = [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (0.5, "0.5"),
            (0.666667, "0.666667"),
            (0.007812, "0.007812"),
            (0.0001, "0.0001"),
            (0.00005, "5e-05"),
            (0.000001, "1e-06"),
            (0.000015, "1.5e-05"),
            (0.000099, "9.9e-05"),
            (1e15, "1000000000000000.0"),
            (1e16, "1e+16"),
            (123.456, "123.456"),
            (0.1 + 0.2, "0.30000000000000004"),
        ];
        for (x, expected) in cases {
            assert_eq!(py_float_repr(x), expected, "repr({x:?})");
        }
    }

    #[test]
    fn six_decimal_rounding_is_half_even_on_the_binary_value() {
        assert_eq!(py_float_repr(py_round6(1.0 / 128.0)), "0.007812");
        assert_eq!(py_float_repr(py_round6(2.0 / 3.0)), "0.666667");
        assert_eq!(py_float_repr(py_round6(0.0000005)), "0.0");
        assert_eq!(py_float_repr(py_round6(0.0000015)), "2e-06");
        assert_eq!(py_float_repr(py_round6(0.0000025)), "3e-06");
        assert_eq!(py_float_repr(py_round6(0.9999995)), "1.0");
        assert_eq!(py_float_repr(py_round6(1e-7)), "0.0");
        assert_eq!(py_float_repr(py_round6(1.0)), "1.0");
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
    fn grouping_matches_format_comma() {
        assert_eq!(py_grouped(0), "0");
        assert_eq!(py_grouped(999), "999");
        assert_eq!(py_grouped(1000), "1,000");
        assert_eq!(py_grouped(1234567), "1,234,567");
        assert_eq!(py_grouped(125_572_619), "125,572,619");
    }

    #[test]
    fn path_str_matches_pathlib() {
        let cases = [
            ("./a//b/", "a/b"),
            ("/data/work/perl/output.txt", "/data/work/perl/output.txt"),
            ("//x", "//x"),
            ("///x", "/x"),
            ("", "."),
            (".", "."),
            ("a/./b", "a/b"),
            ("a/b/.", "a/b"),
            ("./", "."),
            ("x/../y", "x/../y"),
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
}
