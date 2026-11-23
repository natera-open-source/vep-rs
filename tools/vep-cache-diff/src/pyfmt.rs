// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Text the tool writes character for character as CPython does: `float.__repr__`,
//! `round(x, 6)`, `json.dumps(ensure_ascii=True)` string escapes, `str()` and `repr()` of
//! what `json.load` returns, a file name decoded under `surrogateescape` and the order
//! `sorted()` gives such names.

use crate::json::{code_points, Int, Value};
use crate::unicode::is_printable;
use std::cmp::Ordering;
use std::fmt::Write;

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

/// `round(x, 6)`: the double nearest to the six-place decimal nearest to `x`, a tie between
/// two such decimals going to the even one.
pub fn py_round6(x: f64) -> f64 {
    format!("{x:.6}").parse().expect("a fixed-point literal")
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
                write!(
                    out,
                    "\\u{:04x}\\u{:04x}",
                    0xd800 + (offset >> 10),
                    0xdc00 + (offset & 0x3ff)
                )
                .unwrap();
            }
            _ => write!(out, "\\u{unit:04x}").unwrap(),
        }
    }
    out.push('"');
}

/// The code points of the `str` Python decodes a file name's bytes to under
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

/// The order `sorted()` gives two file names: by code point of the name as Python decodes
/// it. Two UTF-8 names compare as their bytes do.
pub fn py_name_order(a: &[u8], b: &[u8]) -> Ordering {
    match (std::str::from_utf8(a), std::str::from_utf8(b)) {
        (Ok(_), Ok(_)) => a.cmp(b),
        _ => surrogateescape(a).cmp(&surrogateescape(b)),
    }
}

fn push_text(out: &mut Vec<u32>, text: &str) {
    out.extend(text.chars().map(u32::from));
}

/// Appends `repr()` of a `str`: single quotes unless the text holds a single quote and no
/// double one, the quote and the backslash escaped, `\t`, `\n` and `\r` by name, and every
/// other code point `str.isprintable()` refuses as `\xNN`, `\uNNNN` or `\UNNNNNNNN`.
fn push_repr_str(out: &mut Vec<u32>, units: &[u32]) {
    let quote = if units.contains(&0x27) && !units.contains(&0x22) {
        0x22
    } else {
        0x27
    };
    out.push(quote);
    for &unit in units {
        match unit {
            0x5c => push_text(out, "\\\\"),
            0x09 => push_text(out, "\\t"),
            0x0a => push_text(out, "\\n"),
            0x0d => push_text(out, "\\r"),
            _ if unit == quote => out.extend([0x5c, quote]),
            _ if is_printable(unit) => out.push(unit),
            0..=0xff => push_text(out, &format!("\\x{unit:02x}")),
            0x100..=0xffff => push_text(out, &format!("\\u{unit:04x}")),
            _ => push_text(out, &format!("\\U{unit:08x}")),
        }
    }
    out.push(quote);
}

/// Appends `repr(value)` of what `json.load` returns.
fn push_repr(out: &mut Vec<u32>, value: &Value) {
    match value {
        Value::Null => push_text(out, "None"),
        Value::Bool(true) => push_text(out, "True"),
        Value::Bool(false) => push_text(out, "False"),
        Value::Int(Int::Small(n)) => push_text(out, &n.to_string()),
        Value::Int(Int::Big(digits)) => push_text(out, digits),
        Value::Float(x) => push_text(out, &py_float_repr(*x)),
        Value::Str(text) => push_repr_str(out, &code_points(text)),
        Value::List(items) => {
            out.push(0x5b);
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    push_text(out, ", ");
                }
                push_repr(out, item);
            }
            out.push(0x5d);
        }
        Value::Object(object) => {
            out.push(0x7b);
            for (i, (key, item)) in object.iter().enumerate() {
                if i > 0 {
                    push_text(out, ", ");
                }
                push_repr_str(out, &code_points(key));
                push_text(out, ": ");
                push_repr(out, item);
            }
            out.push(0x7d);
        }
    }
}

/// The code points of `str(value)`: a string as it stands, anything else as `repr()`.
pub fn py_str(value: &Value) -> Vec<u32> {
    let mut out = Vec::new();
    match value {
        Value::Str(text) => out = code_points(text),
        other => push_repr(&mut out, other),
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::load;

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
            (0.000001, "1e-06"),
            (0.333333, "0.333333"),
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
    fn round_and_fixed_point_match_python() {
        // `repr(round(x, 6))` and `f"{round(x, 6):.4f}"` of CPython for each `x`.
        let cases = [
            (1.0 / 128.0, "0.007812", "0.0078"),
            (3.0 / 128.0, "0.023438", "0.0234"),
            (1.0 / 32.0, "0.03125", "0.0312"),
            (3.0 / 32.0, "0.09375", "0.0938"),
            (2.0 / 3.0, "0.666667", "0.6667"),
            (1.0 / 3.0, "0.333333", "0.3333"),
            (1.0 / 1_000_000.0, "1e-06", "0.0000"),
            (3.0 / 2_000_000.0, "2e-06", "0.0000"),
            (5e-7, "0.0", "0.0000"),
            (1e-7, "0.0", "0.0000"),
            (0.00005, "5e-05", "0.0001"),
            (0.99999949, "0.999999", "1.0000"),
            (0.9999995, "1.0", "1.0000"),
            (644_426.0 / 644_427.0, "0.999998", "1.0000"),
            (0.5, "0.5", "0.5000"),
            (0.125, "0.125", "0.1250"),
            (1.0, "1.0", "1.0000"),
            (0.0, "0.0", "0.0000"),
        ];
        for (x, repr, fixed) in cases {
            let rounded = py_round6(x);
            assert_eq!(py_float_repr(rounded), repr, "round({x:?}, 6)");
            assert_eq!(format!("{rounded:.4}"), fixed, "{x:?} to four places");
        }
    }

    #[test]
    fn json_string_escapes_match_ensure_ascii() {
        let mut out = String::new();
        let text = "\u{1f}\u{7f}\u{e9}\u{1F600}\u{fffd}/\\\"\n\t\u{08}\u{0c}\r\0";
        let mut units: Vec<u32> = text.chars().map(u32::from).collect();
        units.extend([0xd800, 0xdcff]);
        push_json_units(&mut out, units);
        let escape = |hex: &str| format!("{}u{hex}", '\\');
        let expected = format!(
            r#""{}{}{}{}{}{}/\\\"\n\t\b\f\r{}{}{}""#,
            escape("001f"),
            escape("007f"),
            escape("00e9"),
            escape("d83d"),
            escape("de00"),
            escape("fffd"),
            escape("0000"),
            escape("d800"),
            escape("dcff"),
        );
        assert_eq!(out, expected);
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
            surrogateescape(b"a\xe2\x82z\xc3\xa9\x80"),
            [0x61, 0xdce2, 0xdc82, 0x7a, 0xe9, 0xdc80]
        );
    }

    fn written(text: &str) -> String {
        let units = py_str(&load(text.as_bytes()).unwrap());
        let mut out = String::new();
        push_json_units(&mut out, units);
        out
    }

    #[test]
    fn str_of_a_value_matches_python() {
        // `json.dumps(str(json.loads(text)))` of CPython for each text.
        let cases = [
            ("null", r#""None""#),
            ("true", r#""True""#),
            ("false", r#""False""#),
            ("1", r#""1""#),
            ("-0", r#""0""#),
            ("1.5", r#""1.5""#),
            ("1e16", r#""1e+16""#),
            ("1e-5", r#""1e-05""#),
            ("-0.0", r#""-0.0""#),
            ("NaN", r#""nan""#),
            ("-Infinity", r#""-inf""#),
            (
                "10000000000000000000000000",
                r#""10000000000000000000000000""#,
            ),
            (r#""st'r""#, r#""st'r""#),
            ("[]", r#""[]""#),
            ("{}", r#""{}""#),
            (
                r#"["a", "b\"", "c'd\"", 1, null, {"k": [1.0, NaN, Infinity]}]"#,
                r#""['a', 'b\"', 'c\\'d\"', 1, None, {'k': [1.0, nan, inf]}]""#,
            ),
            (
                r#"{"a": "x\n\t\\", "b": true, "a": false}"#,
                r#""{'a': False, 'b': True}""#,
            ),
        ];
        for (text, expected) in cases {
            assert_eq!(written(text), expected, "str({text})");
        }
        // Characters `repr()` escapes by number, and those it leaves as they are.
        let text = format!(
            "[\"\u{e9}\u{200b}{}\u{1F600}\u{7f}\u{80}\u{a0}\u{ad}\u{e0001}\"]",
            "\\ud800"
        );
        let expected = format!(
            "\"['{}{}{}{}{}{}{}{}{}']\"",
            "\\u00e9",
            "\\\\u200b",
            "\\\\ud800",
            "\\ud83d\\ude00",
            "\\\\x7f",
            "\\\\x80",
            "\\\\xa0",
            "\\\\xad",
            "\\\\U000e0001",
        );
        assert_eq!(written(&text), expected);
    }
}
