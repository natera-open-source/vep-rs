// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The text rules of CPython the command line and its messages read by, beside those of
//! `pytext`: `" ".join(text.split())`, `repr` of `bytes`, `str(int(text))` of a number of
//! any length, and the bytes a stream of a given encoding takes for a line.

use crate::pytext::{decimal_value, is_space, split};
pub use crate::pytext::{isdecimal as is_decimal, repr, strip};

/// `" ".join(text.split())`: every run of whitespace as one space, none at either end.
pub fn collapsed(text: &str) -> String {
    split(text).join(" ")
}

/// A character as the `backslashreplace` error handler escapes it.
fn escaped(character: char) -> String {
    match character as u32 {
        code @ 0..=0xff => format!("\\x{code:02x}"),
        code @ 0x100..=0xffff => format!("\\u{code:04x}"),
        code => format!("\\U{code:08x}"),
    }
}

/// `repr(value)` of `bytes`: `b'...'` with the quote rule of `repr(str)`, `\t`, `\n`, `\r`,
/// and every byte outside printable ASCII as `\xNN`.
pub fn repr_bytes(value: &[u8]) -> String {
    let quote = if value.contains(&b'\'') && !value.contains(&b'"') {
        '"'
    } else {
        '\''
    };
    let mut written = String::from("b");
    written.push(quote);
    for &byte in value {
        match byte {
            b'\\' => written.push_str("\\\\"),
            b'\t' => written.push_str("\\t"),
            b'\n' => written.push_str("\\n"),
            b'\r' => written.push_str("\\r"),
            other if other as char == quote => {
                written.push('\\');
                written.push(quote);
            }
            b' '..=b'~' => written.push(byte as char),
            other => written.push_str(&format!("\\x{other:02x}")),
        }
    }
    written.push(quote);
    written
}

/// Whether the character is a decimal digit of some script.
pub fn is_decimal_digit(character: char) -> bool {
    decimal_value(character).is_some()
}

/// The digits a text conversion of `int` takes at most.
const MOST_INTEGER_DIGITS: usize = 4300;

/// `str(int(text))`: the number without its leading zeros, with what the conversion admits
/// around it: whitespace at either end, a sign, the digits of any script and single
/// underscores between digits. `None` where `int` raises.
pub fn integer_text(text: &str) -> Option<String> {
    // `int` drops the whitespace `str.strip()` drops but for U+001C to U+001F.
    let body = text.trim_matches(|c: char| is_space(c) && !matches!(c, '\u{1c}'..='\u{1f}'));
    let (negative, written) = match body.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, body.strip_prefix('+').unwrap_or(body)),
    };
    if written.is_empty()
        || written.starts_with('_')
        || written.ends_with('_')
        || written.contains("__")
    {
        return None;
    }
    let mut digits = String::new();
    for character in written.chars().filter(|&character| character != '_') {
        digits.push(char::from_digit(decimal_value(character)?, 10)?);
    }
    if digits.len() > MOST_INTEGER_DIGITS {
        return None;
    }
    let significant = digits.trim_start_matches('0');
    Some(match (significant.is_empty(), negative) {
        (true, _) => "0".to_string(),
        (false, true) => format!("-{significant}"),
        (false, false) => significant.to_string(),
    })
}

/// What a stream encodes text in.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Encoding {
    Utf8,
    Ascii,
    Latin1,
}

impl Encoding {
    /// The bytes of `text` under the `backslashreplace` error handler: a character the
    /// encoding lacks is written as `\xNN`, `\uNNNN` or `\UNNNNNNNN`.
    pub fn encode(self, text: &str) -> Vec<u8> {
        let most = match self {
            Encoding::Utf8 => return text.as_bytes().to_vec(),
            Encoding::Ascii => 0x7f,
            Encoding::Latin1 => 0xff,
        };
        let mut bytes = Vec::with_capacity(text.len());
        for character in text.chars() {
            if character as u32 <= most {
                bytes.push(character as u32 as u8);
            } else {
                bytes.extend_from_slice(escaped(character).as_bytes());
            }
        }
        bytes
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bytes_are_quoted_as_repr_quotes_them() {
        assert_eq!(repr_bytes(b"xy"), "b'xy'");
        assert_eq!(repr_bytes(b"\x1f"), "b'\\x1f'");
        assert_eq!(repr_bytes(b"'\""), "b'\\'\"'");
        assert_eq!(repr_bytes(b"'"), "b\"'\"");
        assert_eq!(repr_bytes(b"\t\n\\"), "b'\\t\\n\\\\'");
    }

    #[test]
    fn whitespace_is_pythons_set() {
        assert_eq!(strip("\u{1c} a b\u{a0}\u{85}"), "a b");
        assert_eq!(collapsed(" a \t\n b\u{1f}c "), "a b c");
        assert_eq!(strip("\u{200b}a"), "\u{200b}a");
    }

    #[test]
    fn decimal_digits_of_any_script_are_numbers() {
        assert!(is_decimal("12") && is_decimal("٣٤") && is_decimal("１２"));
        assert!(!is_decimal("") && !is_decimal("²") && !is_decimal("1_0") && !is_decimal(" 5"));
        assert_eq!(integer_text("٣٤").as_deref(), Some("34"));
        assert_eq!(integer_text(" 1_0 ").as_deref(), Some("10"));
        assert_eq!(integer_text("+8").as_deref(), Some("8"));
        assert_eq!(integer_text("-5").as_deref(), Some("-5"));
        assert_eq!(integer_text("-000").as_deref(), Some("0"));
        assert_eq!(integer_text("007").as_deref(), Some("7"));
        for refused in ["", "²", "1__0", "_1", "1_", "x", "1.0", "+", "- 5", "0x10"] {
            assert_eq!(integer_text(refused), None, "{refused}");
        }
        let wide = "9".repeat(4300);
        assert_eq!(integer_text(&wide).as_deref(), Some(wide.as_str()));
        assert_eq!(integer_text(&format!("0{wide}")), None);
    }

    #[test]
    fn whitespace_around_a_number_is_what_int_drops() {
        for kept in [
            " 1 ",
            "\t1\n",
            "\u{0b}1\u{0c}\r",
            "\u{a0}1\u{3000}",
            "\u{85}1\u{2028}",
        ] {
            assert_eq!(integer_text(kept).as_deref(), Some("1"), "{kept:?}");
        }
        for refused in ["1\u{1f}", "\u{1c}1", " \u{1d}1", "1\u{1e} ", "\u{200b}1"] {
            assert_eq!(integer_text(refused), None, "{refused:?}");
        }
    }

    #[test]
    fn stream_that_lacks_a_character_takes_its_escape() {
        let line = "keys: Kéy (1) \u{2028} \u{1f600}";
        assert_eq!(Encoding::Utf8.encode(line), line.as_bytes());
        assert_eq!(
            Encoding::Ascii.encode(line),
            b"keys: K\\xe9y (1) \\u2028 \\U0001f600"
        );
        assert_eq!(
            Encoding::Latin1.encode(line),
            b"keys: K\xe9y (1) \\u2028 \\U0001f600"
        );
    }
}
