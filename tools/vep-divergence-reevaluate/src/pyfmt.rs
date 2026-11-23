// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The CPython rules the tool's bytes follow, each under its CPython name: `repr()` of a
//! `str` and of a `list` of them, `float.__repr__`, `int()`, `str.strip()` and
//! `str.splitlines()`, `shlex.quote`, `shlex.join` and `shlex.split`, the `utf-8` codec's
//! `backslashreplace` and strict readings, universal newlines, and the `json.dumps` string
//! escapes.

use std::fmt::Write;
use unicode_general_category::{get_general_category, GeneralCategory};

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

/// `str.rstrip()` with no argument.
pub fn py_rstrip(s: &str) -> &str {
    s.trim_end_matches(is_py_space)
}

/// `str.lstrip()` with no argument.
pub fn py_lstrip(s: &str) -> &str {
    s.trim_start_matches(is_py_space)
}

/// `str.splitlines()`: a line ends at `\n`, `\r`, `\r\n`, `\v`, `\f`, U+001C to U+001E,
/// U+0085, U+2028 or U+2029, and a text that ends on a line end has no empty last line.
pub fn py_splitlines(s: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = s.char_indices().peekable();
    while let Some((i, c)) = chars.next() {
        let is_break = matches!(
            c,
            '\n' | '\r'
                | '\u{0b}'
                | '\u{0c}'
                | '\u{1c}'
                | '\u{1d}'
                | '\u{1e}'
                | '\u{85}'
                | '\u{2028}'
                | '\u{2029}'
        );
        if !is_break {
            continue;
        }
        lines.push(&s[start..i]);
        let mut end = i + c.len_utf8();
        if c == '\r' {
            if let Some((j, '\n')) = chars.peek().copied() {
                chars.next();
                end = j + 1;
            }
        }
        start = end;
    }
    if start < s.len() {
        lines.push(&s[start..]);
    }
    lines
}

/// The last line of a text with its outer whitespace removed, the empty string for a text of
/// whitespace alone: `text.strip().splitlines()[-1]`.
pub fn last_line(text: &str) -> &str {
    py_splitlines(py_strip(text)).last().copied().unwrap_or("")
}

/// `str.isprintable()` of one character: every category but Cc, Cf, Cs, Co, Cn, Zl, Zp and
/// Zs, and the ASCII space.
pub fn is_py_printable(c: char) -> bool {
    if c == ' ' {
        return true;
    }
    !matches!(
        get_general_category(c),
        GeneralCategory::Control
            | GeneralCategory::Format
            | GeneralCategory::Surrogate
            | GeneralCategory::PrivateUse
            | GeneralCategory::Unassigned
            | GeneralCategory::LineSeparator
            | GeneralCategory::ParagraphSeparator
            | GeneralCategory::SpaceSeparator
    )
}

/// `str.isalpha()` of one character: a letter of any of the five letter categories.
pub fn is_py_alpha(c: char) -> bool {
    matches!(
        get_general_category(c),
        GeneralCategory::UppercaseLetter
            | GeneralCategory::LowercaseLetter
            | GeneralCategory::TitlecaseLetter
            | GeneralCategory::ModifierLetter
            | GeneralCategory::OtherLetter
    )
}

/// `str.isdecimal()` of one character: a decimal digit of any script.
pub fn is_py_decimal(c: char) -> bool {
    get_general_category(c) == GeneralCategory::DecimalNumber
}

/// `str.isalnum()` of one character: a letter, a decimal digit, or a character with a digit
/// or numeric value, which are the letter numbers and the other numbers.
pub fn is_py_alnum(c: char) -> bool {
    matches!(
        get_general_category(c),
        GeneralCategory::UppercaseLetter
            | GeneralCategory::LowercaseLetter
            | GeneralCategory::TitlecaseLetter
            | GeneralCategory::ModifierLetter
            | GeneralCategory::OtherLetter
            | GeneralCategory::DecimalNumber
            | GeneralCategory::LetterNumber
            | GeneralCategory::OtherNumber
    )
}

/// `repr(s)` of a `str`: single quotes unless the text holds one and no double quote; the
/// quote, the backslash, `\t`, `\n` and `\r` escaped; every other character that is not
/// printable as `\xNN`, `\uNNNN` or `\UNNNNNNNN`.
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
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            ' '..='~' => out.push(c),
            c if is_py_printable(c) => out.push(c),
            c => {
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

/// `repr(values)` of a `list` of `str`.
pub fn py_repr_str_list<S: AsRef<str>>(values: &[S]) -> String {
    let items: Vec<String> = values.iter().map(|v| py_repr_str(v.as_ref())).collect();
    format!("[{}]", items.join(", "))
}

/// `repr(x)` of a float: the shortest round-trip digits, `.0` on an integral value, exponent
/// form below 1e-4 and from 1e16 with a signed two-digit exponent, and `inf`, `-inf`, `nan`.
pub fn py_float_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    if x.is_sign_negative() {
        return format!("-{}", py_float_repr(-x));
    }
    let exp_form = format!("{x:e}");
    let (mantissa, exponent) = exp_form.split_once('e').expect("LowerExp output");
    let exponent: i32 = exponent.parse().expect("LowerExp exponent");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let decpt = exponent + 1;
    if decpt <= -4 || decpt > 16 {
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

/// `int(text)` for an option that takes a count, clamped to `i64`: the blanks `int()` drops
/// around the number (those of `str.isspace()` but U+001C to U+001F), an optional sign,
/// then decimal digits with single underscores between them, at any length. None where
/// `int()` raises `ValueError`, and for a digit outside ASCII, which `int()` reads.
pub fn py_int_clamped(text: &str) -> Option<i64> {
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
    let mut value = 0i64;
    for byte in digits.bytes() {
        match byte {
            b'_' => {}
            b'0'..=b'9' => {
                value = value
                    .saturating_mul(10)
                    .saturating_add(i64::from(byte - b'0'));
            }
            _ => return None,
        }
    }
    Some(if negative { -value } else { value })
}

/// `format(n, ",")` for a count: its decimal digits in groups of three from the right, a
/// comma between the groups.
pub fn py_int_grouped(n: usize) -> String {
    let digits = n.to_string();
    let mut grouped = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(digit);
    }
    grouped
}

/// `shlex.quote(s)`: the word as it is when it holds only ASCII letters, digits and
/// `@%+=:,./-_`, else in single quotes with each single quote written `'"'"'`.
pub fn shlex_quote(s: &str) -> String {
    if s.is_empty() {
        return "''".to_string();
    }
    let safe = |b: u8| b.is_ascii_alphanumeric() || b"@%+=:,./-_".contains(&b);
    if s.bytes().all(safe) {
        return s.to_string();
    }
    format!("'{}'", s.replace('\'', "'\"'\"'"))
}

/// `shlex.join(words)`.
pub fn shlex_join<S: AsRef<str>>(words: &[S]) -> String {
    let quoted: Vec<String> = words.iter().map(|w| shlex_quote(w.as_ref())).collect();
    quoted.join(" ")
}

/// `shlex.split(s)`: POSIX word splitting on blanks, tabs and line ends, with single quotes,
/// double quotes and the backslash; no comments. The error is the `ValueError` text.
pub fn shlex_split(s: &str) -> Result<Vec<String>, &'static str> {
    #[derive(Clone, Copy, PartialEq)]
    enum State {
        Blank,
        Word,
        Quote(char),
        Escape,
    }
    let mut tokens = Vec::new();
    let mut token = String::new();
    let mut quoted = false;
    let mut state = State::Blank;
    // The state an escape returns to.
    let mut escaped = State::Word;
    let is_blank = |c: char| matches!(c, ' ' | '\t' | '\r' | '\n');
    for c in s.chars() {
        match state {
            State::Blank => {
                if is_blank(c) {
                    // A token cannot be open in this state under POSIX rules.
                } else if c == '\\' {
                    escaped = State::Word;
                    state = State::Escape;
                } else if c == '\'' || c == '"' {
                    state = State::Quote(c);
                } else {
                    token.push(c);
                    state = State::Word;
                }
            }
            State::Quote(q) => {
                quoted = true;
                if c == q {
                    state = State::Word;
                } else if c == '\\' && q == '"' {
                    escaped = state;
                    state = State::Escape;
                } else {
                    token.push(c);
                }
            }
            State::Escape => {
                if let State::Quote(q) = escaped {
                    if c != '\\' && c != q {
                        token.push('\\');
                    }
                }
                token.push(c);
                state = escaped;
            }
            State::Word => {
                if is_blank(c) {
                    state = State::Blank;
                    if !token.is_empty() || quoted {
                        tokens.push(std::mem::take(&mut token));
                        quoted = false;
                    }
                } else if c == '\'' || c == '"' {
                    state = State::Quote(c);
                } else if c == '\\' {
                    escaped = State::Word;
                    state = State::Escape;
                } else {
                    token.push(c);
                }
            }
        }
    }
    match state {
        State::Quote(_) => return Err("No closing quotation"),
        State::Escape => return Err("No escaped character"),
        State::Blank | State::Word => {}
    }
    if !token.is_empty() || quoted {
        tokens.push(token);
    }
    Ok(tokens)
}

/// `bytes.decode("utf-8", "backslashreplace")`: the UTF-8 as it is and every byte outside
/// it as the four characters `\xNN`.
pub fn decode_backslashreplace(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len());
    for chunk in bytes.utf8_chunks() {
        out.push_str(chunk.valid());
        for byte in chunk.invalid() {
            write!(out, "\\x{byte:02x}").unwrap();
        }
    }
    out
}

/// `bytes.decode("utf-8")`: the text, or the `UnicodeDecodeError` message of the first
/// sequence the codec refuses.
pub fn decode_utf8_strict(bytes: &[u8]) -> Result<&str, String> {
    let error = match std::str::from_utf8(bytes) {
        Ok(text) => return Ok(text),
        Err(error) => error,
    };
    let start = error.valid_up_to();
    let rest = &bytes[start..];
    let lead = rest[0];
    let continuation = |b: u8| (0x80..0xc0).contains(&b);
    // (bytes in the refused range, reason), as CPython's decoder states them.
    let (length, reason) = if !(0xc2..0xf5).contains(&lead) {
        (1, "invalid start byte")
    } else {
        let needed = if lead < 0xe0 {
            2
        } else if lead < 0xf0 {
            3
        } else {
            4
        };
        // The second byte's own range for the leads that narrow it.
        let second_ok = |b: u8| match lead {
            0xe0 => (0xa0..0xc0).contains(&b),
            0xed => (0x80..0xa0).contains(&b),
            0xf0 => (0x90..0xc0).contains(&b),
            0xf4 => (0x80..0x90).contains(&b),
            _ => continuation(b),
        };
        if rest.len() < 2 {
            (rest.len(), "unexpected end of data")
        } else if !second_ok(rest[1]) {
            (1, "invalid continuation byte")
        } else if needed == 2 {
            unreachable!("a two-byte sequence with a continuation byte is valid")
        } else if rest.len() < 3 {
            (rest.len(), "unexpected end of data")
        } else if !continuation(rest[2]) {
            (2, "invalid continuation byte")
        } else if needed == 3 {
            unreachable!("a three-byte sequence with two continuation bytes is valid")
        } else if rest.len() < 4 {
            (rest.len(), "unexpected end of data")
        } else {
            (3, "invalid continuation byte")
        }
    };
    Err(if length == 1 {
        format!("'utf-8' codec can't decode byte 0x{lead:02x} in position {start}: {reason}")
    } else {
        format!(
            "'utf-8' codec can't decode bytes in position {start}-{}: {reason}",
            start + length - 1
        )
    })
}

/// Universal newlines as a text-mode read applies them: `\r\n` and a lone `\r` become `\n`.
pub fn universal_newlines(text: &str) -> String {
    if !text.contains('\r') {
        return text.to_string();
    }
    text.replace("\r\n", "\n").replace('\r', "\n")
}

/// The lines a text-mode file yields, each with its `\n` when it has one.
pub fn lines_keepends(text: &str) -> Vec<&str> {
    text.split_inclusive('\n').collect()
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

#[cfg(test)]
#[path = "pyfmt_tests.rs"]
mod tests;
