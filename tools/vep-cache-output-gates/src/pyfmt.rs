// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Text the tool reads and writes character for character as CPython does: `float()` of a
//! command-line word, `float.__repr__`, `round(x, 6)`, `f"{x:.6f}"`,
//! `json.dump(indent=2)` under `ensure_ascii=True`, `str()` of what `json.load` returns and
//! `str(pathlib.Path(arg))`.

use crate::json::{code_points, Value};
use crate::unicode::is_printable;
use std::fmt::Write;

/// `float(word)` for a command-line word in ASCII: blanks at both ends, a sign, digits with
/// single underscores between digits, a fraction, an exponent, or `inf`, `infinity` or
/// `nan` in any case. `None` for any other word.
pub fn py_float(word: &[u8]) -> Option<f64> {
    if !word.is_ascii() {
        return None;
    }
    for (i, byte) in word.iter().enumerate() {
        if *byte == b'_' {
            let digit_at = |at: Option<usize>| {
                at.and_then(|at| word.get(at))
                    .is_some_and(|b| b.is_ascii_digit())
            };
            if !digit_at(i.checked_sub(1)) || !digit_at(Some(i + 1)) {
                return None;
            }
        }
    }
    // The blanks `float()` strips from a word in ASCII: tab to carriage return and the
    // space. U+001C to U+001F are blanks to it only in a word that holds text outside ASCII.
    let blank = |b: &u8| matches!(*b, 0x09..=0x0d | b' ');
    let start = word.iter().position(|b| !blank(b))?;
    let end = word.iter().rposition(|b| !blank(b))? + 1;
    let text: String = word[start..end]
        .iter()
        .filter(|b| **b != b'_')
        .map(|b| char::from(*b))
        .collect();
    let unsigned = text.strip_prefix(['+', '-']).unwrap_or(&text);
    let lower = unsigned.to_ascii_lowercase();
    if lower == "inf" || lower == "infinity" || lower == "nan" {
        return text.parse().ok();
    }
    let (mantissa, exponent) = match unsigned.find(['e', 'E']) {
        Some(at) => (&unsigned[..at], Some(&unsigned[at + 1..])),
        None => (unsigned, None),
    };
    let digits = |part: &str| part.bytes().all(|b| b.is_ascii_digit());
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    if !digits(whole) || !digits(fraction) || whole.len() + fraction.len() == 0 {
        return None;
    }
    if let Some(exponent) = exponent {
        let exponent = exponent.strip_prefix(['+', '-']).unwrap_or(exponent);
        if exponent.is_empty() || !digits(exponent) {
            return None;
        }
    }
    text.parse().ok()
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

/// `round(x, 6)`: the double nearest to the six-place decimal nearest to `x`, a tie between
/// two such decimals going to the even one; `nan` and the infinities as they are.
pub fn py_round6(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    format!("{x:.6}").parse().expect("a fixed-point literal")
}

/// `f"{x:.6f}"`.
pub fn py_fixed6(x: f64) -> String {
    if x.is_nan() {
        "nan".to_string()
    } else if x.is_infinite() {
        if x > 0.0 { "inf" } else { "-inf" }.to_string()
    } else {
        format!("{x:.6}")
    }
}

/// A float as `json.dump` writes one: `float.__repr__`, and `NaN`, `Infinity` and
/// `-Infinity` for the values JSON has no number for.
fn json_float(x: f64) -> String {
    if x.is_nan() {
        "NaN".to_string()
    } else if x.is_infinite() {
        if x > 0.0 { "Infinity" } else { "-Infinity" }.to_string()
    } else {
        py_float_repr(x)
    }
}

/// Appends `json.dumps(s)` with `ensure_ascii=True` for the code points of a `str`: the
/// seven short escapes, every other code point outside `0x20..=0x7e` as `\u` and four
/// lowercase hex digits, a surrogate pair above U+FFFF, a lone surrogate as itself.
fn push_json_units(out: &mut String, units: impl IntoIterator<Item = u32>) {
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

fn push_indent(out: &mut String, level: usize) {
    out.push('\n');
    for _ in 0..level {
        out.push_str("  ");
    }
}

/// Appends `json.dumps(value, indent=2)`, `level` being the nesting the value stands at.
pub fn push_json(out: &mut String, value: &Value, level: usize) {
    match value {
        Value::Null => out.push_str("null"),
        Value::Bool(flag) => out.push_str(if *flag { "true" } else { "false" }),
        Value::Int(n) => write!(out, "{n}").unwrap(),
        Value::BigInt(digits) => out.push_str(digits),
        Value::Float(x) => out.push_str(&json_float(*x)),
        Value::Str(text) => push_json_units(out, code_points(text)),
        Value::List(items) if items.is_empty() => out.push_str("[]"),
        Value::List(items) => {
            out.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                push_indent(out, level + 1);
                push_json(out, item, level + 1);
            }
            push_indent(out, level);
            out.push(']');
        }
        Value::Object(object) if object.is_empty() => out.push_str("{}"),
        Value::Object(object) => {
            out.push('{');
            for (i, (key, item)) in object.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                push_indent(out, level + 1);
                push_json_units(out, code_points(key));
                out.push_str(": ");
                push_json(out, item, level + 1);
            }
            push_indent(out, level);
            out.push('}');
        }
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
        Value::Int(n) => push_text(out, &n.to_string()),
        Value::BigInt(digits) => push_text(out, digits),
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

/// `str(value)` as `print` writes it: a string as it stands, anything else as `repr()`.
/// `None` where the text holds a surrogate, which UTF-8 does not encode.
pub fn py_str(value: &Value) -> Option<String> {
    let mut units = Vec::new();
    match value {
        Value::Str(text) => units = code_points(text),
        other => push_repr(&mut units, other),
    }
    units.into_iter().map(char::from_u32).collect()
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::load;

    #[test]
    fn float_of_a_word_matches_python() {
        // `repr(float(word))` of CPython for each word.
        let cases = [
            ("0.999", "0.999"),
            ("1", "1.0"),
            ("1.", "1.0"),
            (".5", "0.5"),
            ("-.5e1", "-5.0"),
            ("+1E-3", "0.001"),
            ("1e-05", "1e-05"),
            ("0.00001", "1e-05"),
            ("1e16", "1e+16"),
            ("1_000.000_1", "1000.0001"),
            ("1e1_0", "10000000000.0"),
            (" 0.01\n", "0.01"),
            ("\t\x0b\x0c\r\n 7 \n", "7.0"),
            ("1.e5", "100000.0"),
            ("inf", "inf"),
            ("-INF", "-inf"),
            ("+Infinity", "inf"),
            ("nan", "nan"),
            ("-NaN", "nan"),
            ("1e400", "inf"),
            ("-1e400", "-inf"),
            ("1e-400", "0.0"),
            ("-0", "-0.0"),
            ("0000.5", "0.5"),
        ];
        for (word, expected) in cases {
            let read = py_float(word.as_bytes()).map(py_float_repr);
            assert_eq!(read.as_deref(), Some(expected), "float({word:?})");
        }
        // Words `float()` raises `ValueError` for.
        for word in [
            "",
            " ",
            ".",
            "+",
            "e5",
            "1e",
            "1e+",
            "1.5.2",
            "1 .5",
            "- 1",
            "0x10",
            "1_",
            "_1",
            "1__0",
            "1_.5",
            "1._5",
            "1e_5",
            "1_e5",
            "in_f",
            "infinit",
            "infinityx",
            "nan(1)",
            "1,5",
            "1d5",
            "--1",
            "\x1c7",
            "7\x1f",
        ] {
            assert_eq!(py_float(word.as_bytes()), None, "float({word:?})");
        }
        // Digits and blanks outside ASCII, which `float()` reads, are refused.
        assert_eq!(py_float("\u{a0}1".as_bytes()), None);
        assert_eq!(py_float("\u{661}".as_bytes()), None);
        assert_eq!(py_float(b"1\xff"), None);
    }

    #[test]
    fn float_repr_matches_python() {
        let cases = [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (1.5, "1.5"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (0.000015, "1.5e-05"),
            (0.333333, "0.333333"),
            (1e15, "1000000000000000.0"),
            (1e16, "1e+16"),
            (1e22, "1e+22"),
            (123456789012345680.0, "1.2345678901234568e+17"),
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
        // `repr(round(x, 6))` and `f"{x:.6f}"` of CPython for each `x`.
        let cases = [
            (1.0 / 1000.0, "0.001", "0.001000"),
            (2.0 / 3.0, "0.666667", "0.666667"),
            (1.0 / 3.0, "0.333333", "0.333333"),
            (1.0 / 1_000_000.0, "1e-06", "0.000001"),
            (3.0 / 2_000_000.0, "2e-06", "0.000002"),
            (5e-7, "0.0", "0.000000"),
            (2.5e-6, "3e-06", "0.000003"),
            (1e-7, "0.0", "0.000000"),
            (0.9999995, "1.0", "1.000000"),
            (0.99999949, "0.999999", "0.999999"),
            (644_426.0 / 644_427.0, "0.999998", "0.999998"),
            (0.9993331, "0.999333", "0.999333"),
            (1.0, "1.0", "1.000000"),
            (0.0, "0.0", "0.000000"),
            (-0.0, "-0.0", "-0.000000"),
            (123456.7890125, "123456.789012", "123456.789012"),
            (1e22, "1e+22", "10000000000000000000000.000000"),
            (f64::NAN, "nan", "nan"),
            (f64::INFINITY, "inf", "inf"),
            (f64::NEG_INFINITY, "-inf", "-inf"),
        ];
        for (x, repr, fixed) in cases {
            assert_eq!(py_float_repr(py_round6(x)), repr, "round({x:?}, 6)");
            assert_eq!(py_fixed6(x), fixed, "{x:?} to six places");
        }
    }

    fn dumped(text: &str) -> String {
        let mut out = String::new();
        push_json(&mut out, &load(text.as_bytes()).unwrap(), 0);
        out
    }

    #[test]
    fn json_matches_dump_with_an_indent_of_two() {
        // `json.dumps(json.loads(text), indent=2)` of CPython for each text.
        let cases = [
            ("null", "null"),
            ("[true, false]", "[\n  true,\n  false\n]"),
            ("[]", "[]"),
            ("{}", "{}"),
            ("-0", "0"),
            ("-0.0", "-0.0"),
            ("1.0e0", "1.0"),
            ("1E5", "100000.0"),
            ("1e16", "1e+16"),
            ("0.999998", "0.999998"),
            ("[NaN, Infinity, -Infinity, 1e400]", "[\n  NaN,\n  Infinity,\n  -Infinity,\n  Infinity\n]"),
            (
                "170141183460469231731687303715884105728",
                "170141183460469231731687303715884105728",
            ),
            (
                r#"{"a": [], "b": {}, "c": [1, [2, {}]], "d": {"e": null}, "a": "x"}"#,
                "{\n  \"a\": \"x\",\n  \"b\": {},\n  \"c\": [\n    1,\n    [\n      2,\n      {}\n    ]\n  ],\n  \"d\": {\n    \"e\": null\n  }\n}",
            ),
        ];
        for (text, expected) in cases {
            assert_eq!(dumped(text), expected, "{text}");
        }
        let text = format!(
            "\"\u{1f}\u{7f}\u{e9}\u{1F600}\u{fffd}/\\\\\\\"\\n\\t\\b\\f\\r\\u0000{}\"",
            "\\ud800\\udcff\\ud83d\\ude00"
        );
        let escape = |hex: &str| format!("{}u{hex}", '\\');
        let expected = format!(
            r#""{}{}{}{}{}{}/\\\"\n\t\b\f\r{}{}{}{}{}""#,
            escape("001f"),
            escape("007f"),
            escape("00e9"),
            escape("d83d"),
            escape("de00"),
            escape("fffd"),
            escape("0000"),
            escape("d800"),
            escape("dcff"),
            escape("d83d"),
            escape("de00"),
        );
        assert_eq!(dumped(&text.replace('\u{1f}', "\\u001f")), expected);
    }

    #[test]
    fn str_of_a_value_matches_python() {
        // `str(json.loads(text))` of CPython for each text.
        let cases = [
            ("null", "None"),
            ("true", "True"),
            ("false", "False"),
            ("1", "1"),
            ("-0", "0"),
            ("1.5", "1.5"),
            ("1e16", "1e+16"),
            ("-0.0", "-0.0"),
            ("NaN", "nan"),
            ("-Infinity", "-inf"),
            (
                "10000000000000000000000000000000000000000",
                "10000000000000000000000000000000000000000",
            ),
            (r#""st'r""#, "st'r"),
            ("\"\"", ""),
            ("[]", "[]"),
            ("{}", "{}"),
            (
                r#"["a", "b\"", "c'd\"", 1, null, {"k": [1.0, NaN, Infinity]}]"#,
                r#"['a', 'b"', 'c\'d"', 1, None, {'k': [1.0, nan, inf]}]"#,
            ),
            (
                r#"{"a": "x\n\t\\", "b": true, "a": false}"#,
                "{'a': False, 'b': True}",
            ),
            (
                "[\"\u{e9}\u{200b}\u{1F600}\u{7f}\u{80}\u{a0}\u{ad}\u{e0001}\"]",
                "['\u{e9}\\u200b\u{1F600}\\x7f\\x80\\xa0\\xad\\U000e0001']",
            ),
        ];
        for (text, expected) in cases {
            let value = load(text.as_bytes()).unwrap();
            assert_eq!(py_str(&value).as_deref(), Some(expected), "str({text})");
        }
        // A lone surrogate inside a list is written as an escape; one that is the text
        // itself has no UTF-8.
        let inside = load(br#"["\ud800"]"#).unwrap();
        assert_eq!(py_str(&inside).as_deref(), Some("['\\ud800']"));
        let alone = load(br#""a\ud800""#).unwrap();
        assert_eq!(py_str(&alone), None);
    }

    #[test]
    fn path_str_matches_pathlib() {
        let cases = [
            ("./a//b/", "a/b"),
            ("/data/work/verdict.json", "/data/work/verdict.json"),
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
        assert_eq!(py_path_str(b"./ref\xff.json/"), b"ref\xff.json");
    }
}
