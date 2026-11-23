// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The CPython rules the tool's bytes follow, character for character: `str.isspace()`,
//! `str.splitlines()`, `int()`, `repr()` of a `str`, `json.dumps(ensure_ascii=True)` string
//! escapes, `str(pathlib.Path(arg))`, the code point order of two `str` values and the
//! `surrogateescape` reading of a command-line argument or a file name that is not UTF-8.

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

/// The lines of `str.splitlines()`: a line ends at `\n`, `\r`, `\r\n`, U+000B, U+000C, U+001C
/// to U+001E, U+0085, U+2028 or U+2029, and no empty line follows the last line end.
pub fn py_splitlines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = text.char_indices().peekable();
    while let Some((at, c)) = chars.next() {
        let ends = matches!(
            c,
            '\n' | '\r' | '\u{0b}' | '\u{0c}' | '\u{1c}'
                ..='\u{1e}' | '\u{85}' | '\u{2028}' | '\u{2029}'
        );
        if !ends {
            continue;
        }
        lines.push(&text[start..at]);
        start = at + c.len_utf8();
        if c == '\r' && chars.peek().is_some_and(|(_, next)| *next == '\n') {
            chars.next();
            start += 1;
        }
    }
    if start < text.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// An `int()` value: its decimal text as `str(int)` writes it, at any length, and the value
/// as a count, 0 for a negative one and `u64::MAX` for one past it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PyInt {
    pub decimal: String,
    pub count: u64,
}

/// The most digits `int()` reads from a `str` before it refuses the conversion.
const INT_MAX_STR_DIGITS: usize = 4300;

/// The first code point of each run of ten decimal digits outside ASCII, as Unicode 16.0
/// lists them: the digits `int()` reads beside `0` to `9`.
const DECIMAL_ZEROS: [u32; 75] = [
    0x0660, 0x06f0, 0x07c0, 0x0966, 0x09e6, 0x0a66, 0x0ae6, 0x0b66, 0x0be6, 0x0c66, 0x0ce6, 0x0d66,
    0x0de6, 0x0e50, 0x0ed0, 0x0f20, 0x1040, 0x1090, 0x17e0, 0x1810, 0x1946, 0x19d0, 0x1a80, 0x1a90,
    0x1b50, 0x1bb0, 0x1c40, 0x1c50, 0xa620, 0xa8d0, 0xa900, 0xa9d0, 0xa9f0, 0xaa50, 0xabf0, 0xff10,
    0x104a0, 0x10d30, 0x10d40, 0x11066, 0x110f0, 0x11136, 0x111d0, 0x112f0, 0x11450, 0x114d0,
    0x11650, 0x116c0, 0x116d0, 0x116da, 0x11730, 0x118e0, 0x11950, 0x11bf0, 0x11c50, 0x11d50,
    0x11da0, 0x11f50, 0x16130, 0x16a60, 0x16ac0, 0x16b50, 0x16d70, 0x1ccf0, 0x1d7ce, 0x1d7d8,
    0x1d7e2, 0x1d7ec, 0x1d7f6, 0x1e140, 0x1e2f0, 0x1e4f0, 0x1e5f1, 0x1e950, 0x1fbf0,
];

/// The ASCII character `int()` reads `c` as before it parses: a decimal digit of any script
/// as its ASCII digit, a blank of `str.isspace()` outside ASCII as a space, an ASCII
/// character as itself. None for any other character.
fn int_char(c: char) -> Option<char> {
    if c.is_ascii() {
        return Some(c);
    }
    if is_py_space(c) {
        return Some(' ');
    }
    let code = u32::from(c);
    let run = DECIMAL_ZEROS
        .partition_point(|zero| *zero <= code)
        .checked_sub(1)?;
    let digit = code - DECIMAL_ZEROS[run];
    (digit < 10).then(|| char::from(b'0' + digit as u8))
}

/// Whether `c` is a decimal digit of any script: the `\d` of `re` on a `str`.
pub fn is_py_decimal(c: char) -> bool {
    int_char(c).is_some_and(|read| read.is_ascii_digit())
}

/// `int(text)`: blanks around the number (those of `str.isspace()` but U+001C to U+001F),
/// an optional sign, then decimal digits of any script with single underscores between
/// them, up to 4,300 digits. None where `int()` raises `ValueError`.
pub fn py_int(text: &str) -> Option<PyInt> {
    let ascii: String = text.chars().map(int_char).collect::<Option<String>>()?;
    let body = ascii.trim_matches(|c| matches!(c, '\t'..='\r' | ' '));
    let (negative, digits) = match body.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, body.strip_prefix('+').unwrap_or(body)),
    };
    if digits.is_empty()
        || digits.starts_with('_')
        || digits.ends_with('_')
        || digits.contains("__")
        || !digits.bytes().all(|b| b == b'_' || b.is_ascii_digit())
    {
        return None;
    }
    let plain: String = digits.chars().filter(|c| *c != '_').collect();
    if plain.len() > INT_MAX_STR_DIGITS {
        return None;
    }
    let significant = plain.trim_start_matches('0');
    if significant.is_empty() {
        return Some(PyInt {
            decimal: "0".to_string(),
            count: 0,
        });
    }
    let (decimal, count) = if negative {
        (format!("-{significant}"), 0)
    } else {
        (
            significant.to_string(),
            significant.parse().unwrap_or(u64::MAX),
        )
    };
    Some(PyInt { decimal, count })
}

/// The code points of the `str` Python holds for these bytes: the UTF-8 in them, and each
/// byte outside UTF-8 as the lone surrogate U+DC00 plus the byte.
fn code_points(bytes: &[u8]) -> impl Iterator<Item = u32> + '_ {
    bytes.utf8_chunks().flat_map(|chunk| {
        chunk
            .valid()
            .chars()
            .map(u32::from)
            .chain(chunk.invalid().iter().map(|byte| 0xdc00 + u32::from(*byte)))
    })
}

/// The order of the two `str` values Python holds for these bytes, by code point.
pub fn py_str_cmp(a: &[u8], b: &[u8]) -> Ordering {
    code_points(a).cmp(code_points(b))
}

/// `repr()` of the `str` Python holds for these bytes: double quotes where the text holds a
/// single quote and no double one, `\\`, `\t`, `\n`, `\r`, `\xNN` below U+0020, at U+007F
/// and from U+0080 to U+00A0 and U+00AD, `\uNNNN` for the other blanks of `str.isspace()`,
/// and `\udcNN` for a byte outside UTF-8. Every other character is written as it is; of
/// those, `repr()` escapes the ones Unicode does not class as printable.
pub fn py_repr(bytes: &[u8]) -> String {
    let quote = if bytes.contains(&b'\'') && !bytes.contains(&b'"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(bytes.len() + 2);
    out.push(quote);
    for chunk in bytes.utf8_chunks() {
        for c in chunk.valid().chars() {
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
                '\0'..='\u{1f}' | '\u{7f}'..='\u{a0}' | '\u{ad}' => {
                    write!(out, "\\x{:02x}", u32::from(c)).unwrap();
                }
                c if is_py_space(c) => write!(out, "\\u{:04x}", u32::from(c)).unwrap(),
                c => out.push(c),
            }
        }
        for byte in chunk.invalid() {
            write!(out, "\\udc{byte:02x}").unwrap();
        }
    }
    out.push(quote);
    out
}

/// Appends `json.dumps(s)` with `ensure_ascii=True`: the seven short escapes, every other
/// code point outside `0x20..=0x7e` as `\u` and four lowercase hex digits, a surrogate pair
/// above U+FFFF.
pub fn push_json_str(out: &mut String, s: &str) {
    out.push('"');
    push_json_chars(out, s);
    out.push('"');
}

/// Appends `json.dumps(s)` for the `str` Python holds for these bytes: the UTF-8 in them as
/// `push_json_str` writes it and each byte outside UTF-8 as the lone surrogate `\udcXX`.
pub fn push_json_os_str(out: &mut String, bytes: &[u8]) {
    out.push('"');
    for chunk in bytes.utf8_chunks() {
        push_json_chars(out, chunk.valid());
        for byte in chunk.invalid() {
            write!(out, "\\udc{byte:02x}").unwrap();
        }
    }
    out.push('"');
}

fn push_json_chars(out: &mut String, s: &str) {
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
}

/// The text Python's stderr writes for the `str` it holds for these bytes, under its
/// `backslashreplace` handler: the UTF-8 as it is and each byte outside UTF-8 as the six
/// characters `\udcXX`.
pub fn py_stderr_os_str(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len());
    for chunk in bytes.utf8_chunks() {
        out.push_str(chunk.valid());
        for byte in chunk.invalid() {
            write!(out, "\\udc{byte:02x}").unwrap();
        }
    }
    out
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

    #[test]
    fn space_set_is_pythons_isspace() {
        let expected: Vec<u32> = [
            0x9, 0xa, 0xb, 0xc, 0xd, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x85, 0xa0, 0x1680, 0x2000,
            0x2001, 0x2002, 0x2003, 0x2004, 0x2005, 0x2006, 0x2007, 0x2008, 0x2009, 0x200a, 0x2028,
            0x2029, 0x202f, 0x205f, 0x3000,
        ]
        .to_vec();
        let actual: Vec<u32> = (0..=0x10ffffu32)
            .filter_map(char::from_u32)
            .filter(|c| is_py_space(*c))
            .map(u32::from)
            .collect();
        assert_eq!(actual, expected);
    }

    #[test]
    fn lines_end_where_splitlines_ends_them() {
        assert_eq!(py_splitlines(""), Vec::<&str>::new());
        assert_eq!(py_splitlines("a"), ["a"]);
        assert_eq!(py_splitlines("a\n"), ["a"]);
        assert_eq!(py_splitlines("a\n\nb"), ["a", "", "b"]);
        assert_eq!(py_splitlines("a\r\nb\rc\n"), ["a", "b", "c"]);
        assert_eq!(py_splitlines("\r\n\r"), ["", ""]);
        assert_eq!(
            py_splitlines("a\u{0b}b\u{0c}c\u{1c}d\u{1d}e\u{1e}f\u{85}g\u{2028}h\u{2029}i\u{1f}j"),
            ["a", "b", "c", "d", "e", "f", "g", "h", "i\u{1f}j"]
        );
    }

    #[test]
    fn an_integer_is_read_as_int_reads_it() {
        let read = [
            ("20", "20", 20),
            ("0", "0", 0),
            ("-0", "0", 0),
            ("+2", "2", 2),
            ("-1", "-1", 0),
            ("007", "7", 7),
            ("1_000", "1000", 1000),
            (" 3 ", "3", 3),
            ("\t\n\u{0b}\u{0c}\r 3\n", "3", 3),
            ("\u{a0}3\u{3000}", "3", 3),
            ("\u{ff13}", "3", 3),
            ("\u{669}\u{660}_\u{6f1}", "901", 901),
            ("-\u{1d7ce}\u{1d7ff}\u{1fbf9}", "-99", 0),
            ("\u{85}+\u{967}0\u{2028}", "10", 10),
            ("18446744073709551615", "18446744073709551615", u64::MAX),
            (
                "99999999999999999999999999",
                "99999999999999999999999999",
                u64::MAX,
            ),
            (
                "-99999999999999999999999999",
                "-99999999999999999999999999",
                0,
            ),
        ];
        for (text, decimal, count) in read {
            let value = py_int(text).unwrap_or_else(|| panic!("int({text:?})"));
            assert_eq!((value.decimal.as_str(), value.count), (decimal, count));
        }
        let refused = [
            "",
            " ",
            "x",
            "1.5",
            "1e2",
            "0x10",
            "_1",
            "1_",
            "1__0",
            "+_1",
            "- 3",
            "+-3",
            "--3",
            "-",
            "+",
            "1 0",
            "\u{1c}3",
            "3\u{1f}",
            "\u{b2}",
            "\u{2460}",
            "\u{65f}",
            "\u{66a}",
            "\u{1d7cd}",
            "\u{1d800}",
            "3\u{a0}4",
            "\u{feff}5",
            "5\0",
        ];
        for text in refused {
            assert_eq!(py_int(text), None, "int({text:?})");
        }
    }

    #[test]
    fn a_decimal_digit_of_any_script_is_a_digit_and_no_other_character_is() {
        for digit in ['0', '9', '\u{663}', '\u{ff15}', '\u{116e3}', '\u{1fbf9}'] {
            assert!(is_py_decimal(digit), "{digit:?}");
        }
        for other in [
            'a',
            '-',
            '.',
            ' ',
            '\u{a0}',
            '\u{b2}',
            '\u{bd}',
            '\u{2460}',
            '\u{1d7cd}',
        ] {
            assert!(!is_py_decimal(other), "{other:?}");
        }
    }

    #[test]
    fn an_integer_is_read_up_to_the_digit_limit() {
        let most = "1".repeat(4300);
        assert_eq!(py_int(&most).unwrap().decimal, most);
        assert_eq!(
            py_int(&format!("-{most}")).unwrap().decimal,
            format!("-{most}")
        );
        assert_eq!(py_int(&"0".repeat(4300)).unwrap().decimal, "0");
        assert_eq!(py_int(&"1".repeat(4301)), None);
        assert_eq!(py_int(&"0".repeat(4301)), None);
        let underscored = format!("{}1", "1_".repeat(4299));
        assert_eq!(py_int(&underscored).unwrap().decimal, most);
    }

    #[test]
    fn strings_order_by_code_point_with_surrogates_for_stray_bytes() {
        assert_eq!(py_str_cmp(b"a", b"b"), Ordering::Less);
        assert_eq!(py_str_cmp(b"a", b"a"), Ordering::Equal);
        assert_eq!(py_str_cmp(b"ab", b"a"), Ordering::Greater);
        assert_eq!(py_str_cmp(b"Z", b"a"), Ordering::Less);
        // U+E000 and U+1F600 sort after the lone surrogate a stray byte is read as, though
        // their first UTF-8 bytes are smaller than 0xff.
        assert_eq!(py_str_cmp(b"\xff", "\u{e000}".as_bytes()), Ordering::Less);
        assert_eq!(py_str_cmp(b"\xff", "\u{1F600}".as_bytes()), Ordering::Less);
        assert_eq!(
            py_str_cmp(b"\xff", "\u{d7ff}".as_bytes()),
            Ordering::Greater
        );
        assert_eq!(py_str_cmp(b"\x80", b"\xff"), Ordering::Less);
    }

    #[test]
    fn repr_quotes_and_escapes_as_python_does() {
        let cases: [(&[u8], &str); 12] = [
            (b"", "''"),
            (b",", "','"),
            (b" , ", "' , '"),
            (b"it's", r#""it's""#),
            (b"\"q\"", r#"'"q"'"#),
            (b"a'b\"c", r#"'a\'b"c'"#),
            (b"a\\b", r"'a\\b'"),
            (
                b"\t\n\r\x0b\x0c\x1c\x1f\x00\x7f",
                r"'\t\n\r\x0b\x0c\x1c\x1f\x00\x7f'",
            ),
            (
                "\u{85}\u{a0}\u{ad}\u{1680}\u{2003}\u{2028}\u{3000}".as_bytes(),
                concat!(
                    r"'\x85\xa0\xad",
                    r"\",
                    "u1680",
                    r"\",
                    "u2003",
                    r"\",
                    "u2028",
                    r"\",
                    "u3000'"
                ),
            ),
            (
                "d\u{e9}j\u{e0}/\u{4e2d}".as_bytes(),
                "'d\u{e9}j\u{e0}/\u{4e2d}'",
            ),
            (b"set\xff\xe2\x82.txt", r"'set\udcff\udce2\udc82.txt'"),
            (b"plain/set.txt", "'plain/set.txt'"),
        ];
        for (bytes, expected) in cases {
            assert_eq!(py_repr(bytes), expected, "{bytes:?}");
        }
    }

    #[test]
    fn json_string_escapes_match_ensure_ascii() {
        let mut out = String::new();
        push_json_str(
            &mut out,
            "\u{1f}\u{7f}\u{e9}\u{1F600}\u{fffd}/\\\"\n\t\u{08}\u{0c}\r\0",
        );
        let expected = concat!(
            r#""\u001f\u007f"#,
            r"\",
            "u00e9",
            r"\",
            "ud83d",
            r"\",
            "ude00",
            r"\",
            "ufffd",
            r#"/\\\"\n\t\b\f\r\u0000""#
        );
        assert_eq!(out, expected);
    }

    #[test]
    fn bytes_outside_utf8_read_as_surrogateescape_reads_them() {
        let cases: [(&[u8], &str, &str); 4] = [
            (
                b"in\xff\xe2\x82.vcf",
                r#""in\udcff\udce2\udc82.vcf""#,
                r"in\udcff\udce2\udc82.vcf",
            ),
            (
                b"d\xc0\x80r/x\xc3\xa9.vcf",
                concat!(r#""d\udcc0\udc80r/x"#, r"\", r#"u00e9.vcf""#),
                "d\\udcc0\\udc80r/x\u{e9}.vcf",
            ),
            (b"plain/name.vcf", r#""plain/name.vcf""#, "plain/name.vcf"),
            (
                b"\xf0\x9f\x98\x80\x80",
                concat!(r#"""#, r"\", "ud83d", r"\", r#"ude00\udc80""#),
                "\u{1F600}\\udc80",
            ),
        ];
        for (bytes, json, stderr) in cases {
            let mut out = String::new();
            push_json_os_str(&mut out, bytes);
            assert_eq!(out, json, "{bytes:?}");
            assert_eq!(py_stderr_os_str(bytes), stderr, "{bytes:?}");
        }
    }

    #[test]
    fn path_str_matches_pathlib() {
        let cases = [
            ("./a//b/", "a/b"),
            ("/data/work/inputs", "/data/work/inputs"),
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
        assert_eq!(py_path_str(b"./in\xff/"), b"in\xff");
    }
}
