// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Text as CPython writes and reads it: `str(pathlib.Path(...))`, a file read in text
//! mode, `repr()` of a `str`, `json.dumps` of a string, `{:,}` of an integer.

use std::fmt::Write;

use unicode_general_category::{get_general_category, GeneralCategory};

/// `str(pathlib.PurePosixPath(arg))`: repeated separators and `.` components dropped, a
/// trailing separator dropped, exactly two leading slashes kept, `..` left alone, `.`
/// for an empty result.
pub fn path_str(arg: &[u8]) -> Vec<u8> {
    let leading_slashes = arg.iter().take_while(|byte| **byte == b'/').count();
    let root: &[u8] = match leading_slashes {
        0 => b"",
        2 => b"//",
        _ => b"/",
    };
    let parts: Vec<&[u8]> = arg
        .split(|byte| *byte == b'/')
        .filter(|part| !part.is_empty() && *part != b".")
        .collect();
    let joined = [root, &parts.join(&b'/')].concat();
    if joined.is_empty() {
        b".".to_vec()
    } else {
        joined
    }
}

/// `path / name` for a normalised path and a name without a separator.
pub fn path_join(path: &[u8], name: &[u8]) -> Vec<u8> {
    let mut joined = if path == b"." {
        Vec::new()
    } else {
        path.to_vec()
    };
    if !joined.is_empty() && joined.last() != Some(&b'/') {
        joined.push(b'/');
    }
    joined.extend_from_slice(name);
    joined
}

/// A file read in text mode with universal newlines: `\r\n` and a lone `\r` are `\n`.
pub fn universal_newlines(text: &str) -> String {
    if !text.contains('\r') {
        return text.to_string();
    }
    let mut out = String::with_capacity(text.len());
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\r' {
            if chars.peek() == Some(&'\n') {
                chars.next();
            }
            out.push('\n');
        } else {
            out.push(c);
        }
    }
    out
}

/// `str.isprintable()` of one character: every category but Cc, Cf, Cs, Co, Cn, Zl, Zp
/// and Zs, and the ASCII space.
pub fn is_printable(c: char) -> bool {
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

/// What `\d` matches in a `str` pattern: a decimal digit of any script.
pub fn is_decimal(c: char) -> bool {
    get_general_category(c) == GeneralCategory::DecimalNumber
}

/// `repr(s)` of a `str`: single quotes unless the text holds one and no double quote;
/// the quote, the backslash, `\t`, `\n` and `\r` escaped; every other character that is
/// not printable as `\xNN`, `\uNNNN` or `\UNNNNNNNN`.
pub fn repr_str(s: &str) -> String {
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
            c if is_printable(c) => out.push(c),
            c => {
                let code = c as u32;
                if code <= 0xff {
                    write!(out, "\\x{code:02x}").expect("a string takes text");
                } else if code <= 0xffff {
                    write!(out, "\\u{code:04x}").expect("a string takes text");
                } else {
                    write!(out, "\\U{code:08x}").expect("a string takes text");
                }
            }
        }
    }
    out.push(quote);
    out
}

/// `repr(values)` of a `list` of `str`.
pub fn repr_str_list<S: AsRef<str>>(values: &[S]) -> String {
    let items: Vec<String> = values
        .iter()
        .map(|value| repr_str(value.as_ref()))
        .collect();
    format!("[{}]", items.join(", "))
}

/// A string as `json.dumps` writes it under `ensure_ascii`: `"` and `\` escaped, and
/// every character outside space to `~` as `\n`, `\r`, `\t`, `\b`, `\f` or `\uNNNN`, in
/// two halves past U+FFFF.
pub fn push_json_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            ' '..='~' => out.push(c),
            c => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units) {
                    write!(out, "\\u{unit:04x}").expect("a string takes text");
                }
            }
        }
    }
    out.push('"');
}

/// `f"{n:,}"`: the digits in groups of three.
pub fn grouped(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            out.push(',');
        }
        out.push(digit);
    }
    out
}

/// The text of an operating system error without Rust's ` (os error <n>)`.
pub fn os_text(error: &std::io::Error) -> String {
    let text = error.to_string();
    match text.rfind(" (os error ") {
        Some(at) if text.ends_with(')') => text[..at].to_string(),
        _ => text,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path(arg: &str) -> String {
        String::from_utf8(path_str(arg.as_bytes())).expect("UTF-8")
    }

    /// Each expectation is `str(pathlib.Path(arg))` of CPython 3.14.
    #[test]
    fn a_path_is_normalised_as_pathlib_writes_it() {
        assert_eq!(path("a/b"), "a/b");
        assert_eq!(path("a//b/"), "a/b");
        assert_eq!(path("./a/./b/."), "a/b");
        assert_eq!(path(""), ".");
        assert_eq!(path("."), ".");
        assert_eq!(path("/"), "/");
        assert_eq!(path("//a"), "//a");
        assert_eq!(path("///a"), "/a");
        assert_eq!(path("a/../b"), "a/../b");
    }

    #[test]
    fn a_name_joins_a_path_with_one_separator() {
        let join = |path: &str, name: &str| {
            String::from_utf8(path_join(path.as_bytes(), name.as_bytes())).expect("UTF-8")
        };
        assert_eq!(join("out", "x.tsv"), "out/x.tsv");
        assert_eq!(join(".", "x.tsv"), "x.tsv");
        assert_eq!(join("/", "x.tsv"), "/x.tsv");
    }

    /// Each expectation is `repr(s)` of CPython 3.14 (Unicode 16.0.0).
    #[test]
    fn repr_quotes_and_escapes_as_cpython_does() {
        assert_eq!(repr_str("a\tb\n"), "'a\\tb\\n'");
        assert_eq!(repr_str("it's"), "\"it's\"");
        assert_eq!(repr_str("it's \"x\""), "'it\\'s \"x\"'");
        assert_eq!(repr_str("back\\slash"), "'back\\\\slash'");
        assert_eq!(
            repr_str("\u{0}\u{1f}\u{7f}\u{a0}\u{ad}"),
            "'\\x00\\x1f\\x7f\\xa0\\xad'"
        );
        assert_eq!(repr_str("\u{e9}\u{263a}"), "'\u{e9}\u{263a}'");
        assert_eq!(
            repr_str("\u{2028}\u{e000}\u{1f600}\u{e0001}"),
            "'\\u2028\\ue000\u{1f600}\\U000e0001'"
        );
        assert_eq!(repr_str_list(&["a", "b'c"]), "['a', \"b'c\"]");
    }

    /// Each expectation is `json.dumps(s)` of CPython 3.14.
    #[test]
    fn a_json_string_is_ascii() {
        let json = |s: &str| {
            let mut out = String::new();
            push_json_str(&mut out, s);
            out
        };
        assert_eq!(json("a\"b\\c"), "\"a\\\"b\\\\c\"");
        assert_eq!(
            json("\n\r\t\u{8}\u{c}\u{1}\u{7f}"),
            "\"\\n\\r\\t\\b\\f\\u0001\\u007f\""
        );
        assert_eq!(json("\u{e9}\u{1f600}"), "\"\\u00e9\\ud83d\\ude00\"");
    }

    #[test]
    fn digits_are_grouped_by_three() {
        assert_eq!(grouped(0), "0");
        assert_eq!(grouped(999), "999");
        assert_eq!(grouped(1000), "1,000");
        assert_eq!(grouped(195_232), "195,232");
        assert_eq!(grouped(1_234_567), "1,234,567");
    }

    #[test]
    fn universal_newlines_make_every_line_end_a_line_feed() {
        assert_eq!(universal_newlines("a\r\nb\rc\n\r\n"), "a\nb\nc\n\n");
    }
}
