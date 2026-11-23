// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A JSON document read as a stream, by the grammar of Python's `json.load`: strict JSON,
//! plus the constants `NaN`, `Infinity` and `-Infinity`, plus `\u` escapes that spell a lone
//! surrogate; no byte order mark, no control character inside a string, only space, tab,
//! line feed and carriage return between tokens, one value a document, and no integer of
//! more than 4,300 digits. The text is UTF-8 and refused where it is not.
//!
//! The scanner holds a window of the input, not the document: a caller walks the values it
//! reads and skips the rest, and every skipped value is still held to the grammar. A string
//! is unescaped into bytes that are UTF-8 but for a lone surrogate, which is written in
//! three bytes as UTF-8 would write its code point; such bytes keep the order and the
//! equality of the code points.

use crate::input::READ_BUF;
use crate::pyfmt::py_float_repr;
use memchr::memchr2;
use std::io::{self, Read};

/// The integer digits Python's `int` converts from text.
const INT_MAX_STR_DIGITS: usize = 4300;

/// The type of a value that was read past.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Null,
    Bool,
    Number,
    Str,
    Array,
    Object,
}

/// A value that was read past: its type and whether Python counts it true (a non-empty
/// string, array or object, `true`, a number other than zero, `NaN` included).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Skipped {
    pub kind: Kind,
    pub truthy: bool,
}

fn refuse(what: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, what)
}

pub struct Scanner<R> {
    reader: R,
    buf: Vec<u8>,
    pos: usize,
    end: usize,
    eof: bool,
    nesting: Vec<bool>,
}

impl<R: Read> Scanner<R> {
    /// A scanner that reads through `window`: the buffer an earlier scanner gave back, or
    /// an empty one.
    pub fn new(reader: R, mut window: Vec<u8>) -> Scanner<R> {
        if window.len() < READ_BUF {
            window.resize(READ_BUF, 0);
        }
        Scanner::over(reader, window)
    }

    /// The window, at the size the longest token grew it to, for the thread's next scanner.
    pub fn into_window(self) -> Vec<u8> {
        self.buf
    }

    fn over(reader: R, window: Vec<u8>) -> Scanner<R> {
        Scanner {
            reader,
            buf: window,
            pos: 0,
            end: 0,
            eof: false,
            nesting: Vec::new(),
        }
    }

    /// Reads more input behind the unread bytes, which move to the front of the window and
    /// keep their offsets from `pos`; false at the end of the input.
    fn fill(&mut self) -> io::Result<bool> {
        if self.eof {
            return Ok(false);
        }
        if self.pos > 0 {
            self.buf.copy_within(self.pos..self.end, 0);
            self.end -= self.pos;
            self.pos = 0;
        }
        if self.end == self.buf.len() {
            self.buf.resize(self.buf.len() * 2, 0);
        }
        loop {
            match self.reader.read(&mut self.buf[self.end..]) {
                Ok(0) => {
                    self.eof = true;
                    return Ok(false);
                }
                Ok(n) => {
                    self.end += n;
                    return Ok(true);
                }
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
    }

    /// The unread byte `offset` bytes ahead, `None` past the end of the input.
    fn at(&mut self, offset: usize) -> io::Result<Option<u8>> {
        while self.pos + offset >= self.end {
            if !self.fill()? {
                return Ok(None);
            }
        }
        Ok(Some(self.buf[self.pos + offset]))
    }

    fn peek(&mut self) -> io::Result<Option<u8>> {
        self.at(0)
    }

    /// Whether the unread input starts with `text`.
    fn starts_with(&mut self, text: &[u8]) -> io::Result<bool> {
        if text.is_empty() {
            return Ok(true);
        }
        Ok(
            self.at(text.len() - 1)?.is_some()
                && self.buf[self.pos..self.pos + text.len()] == *text,
        )
    }

    fn skip_whitespace(&mut self) -> io::Result<()> {
        loop {
            while self.pos < self.end {
                match self.buf[self.pos] {
                    b' ' | b'\t' | b'\n' | b'\r' => self.pos += 1,
                    _ => return Ok(()),
                }
            }
            if !self.fill()? {
                return Ok(());
            }
        }
    }

    /// Starts the document: refuses a byte order mark and stops at the value's first byte.
    pub fn begin(&mut self) -> io::Result<()> {
        if self.starts_with(b"\xef\xbb\xbf")? {
            return Err(refuse("byte order mark"));
        }
        self.skip_whitespace()
    }

    /// Ends the document: nothing but whitespace may follow its value.
    pub fn finish(&mut self) -> io::Result<()> {
        self.skip_whitespace()?;
        match self.peek()? {
            None => Ok(()),
            Some(_) => Err(refuse("extra data after the value")),
        }
    }

    /// The first byte of the value the scanner stands at.
    pub fn value_start(&mut self) -> io::Result<u8> {
        self.peek()?.ok_or_else(|| refuse("expecting a value"))
    }

    /// Enters the array or object whose opening byte the scanner stands at; false when it
    /// is empty, and then it is closed as well.
    pub fn enter(&mut self, closer: u8) -> io::Result<bool> {
        self.pos += 1;
        self.skip_whitespace()?;
        if self.peek()? == Some(closer) {
            self.pos += 1;
            return Ok(false);
        }
        Ok(true)
    }

    /// After a member of an array or object: true at the next member, false when the
    /// closing byte ended the container.
    pub fn another(&mut self, closer: u8) -> io::Result<bool> {
        self.skip_whitespace()?;
        match self.peek()? {
            Some(b',') => {
                self.pos += 1;
                self.skip_whitespace()?;
                Ok(true)
            }
            Some(c) if c == closer => {
                self.pos += 1;
                Ok(false)
            }
            _ => Err(refuse("expecting ',' or the end of the container")),
        }
    }

    /// Reads an object member's name and the `:` behind it, and stops at the member's
    /// value. The name is appended to `name`.
    pub fn key(&mut self, name: &mut Vec<u8>) -> io::Result<()> {
        self.member_name(Some(name))
    }

    fn member_name(&mut self, name: Option<&mut Vec<u8>>) -> io::Result<()> {
        if self.peek()? != Some(b'"') {
            return Err(refuse("expecting a member name in double quotes"));
        }
        self.string(name)?;
        self.skip_whitespace()?;
        if self.peek()? != Some(b':') {
            return Err(refuse("expecting ':' after a member name"));
        }
        self.pos += 1;
        self.skip_whitespace()
    }

    /// Reads the string whose opening quote the scanner stands at, appending its unescaped
    /// content to `out`; true when the string is not empty.
    pub fn string(&mut self, mut out: Option<&mut Vec<u8>>) -> io::Result<bool> {
        self.pos += 1;
        let mut filled = false;
        let mut clean = 0;
        loop {
            let Some(found) = memchr2(b'"', b'\\', &self.buf[self.pos + clean..self.end]) else {
                clean = self.end - self.pos;
                if !self.fill()? {
                    return Err(refuse("unterminated string"));
                }
                continue;
            };
            let stop = self.pos + clean + found;
            let run = &self.buf[self.pos..stop];
            // A run ends at an ASCII byte, so it is UTF-8 exactly where the document is.
            let control = run.iter().fold(false, |any, b| any | (*b < 0x20));
            if control || std::str::from_utf8(run).is_err() {
                return Err(refuse("control character or invalid UTF-8 in a string"));
            }
            filled |= !run.is_empty();
            if let Some(out) = out.as_deref_mut() {
                out.extend_from_slice(run);
            }
            let quote = self.buf[stop] == b'"';
            self.pos = stop + 1;
            clean = 0;
            if quote {
                return Ok(filled);
            }
            filled = true;
            let code_point = self.escape()?;
            if let Some(out) = out.as_deref_mut() {
                push_code_point(out, code_point);
            }
        }
    }

    /// The code point of the escape whose backslash was read. A high surrogate joins a low
    /// surrogate escape right behind it; any other surrogate stands alone.
    fn escape(&mut self) -> io::Result<u32> {
        let Some(letter) = self.peek()? else {
            return Err(refuse("unterminated string"));
        };
        self.pos += 1;
        Ok(match letter {
            b'"' | b'\\' | b'/' => u32::from(letter),
            b'b' => 0x08,
            b'f' => 0x0c,
            b'n' => 0x0a,
            b'r' => 0x0d,
            b't' => 0x09,
            b'u' => {
                let high = self.hex4(0)?;
                self.pos += 4;
                if (0xd800..0xdc00).contains(&high) && self.starts_with(b"\\u")? {
                    let low = self.hex4(2)?;
                    if (0xdc00..0xe000).contains(&low) {
                        self.pos += 6;
                        return Ok(0x10000 + ((high - 0xd800) << 10) + (low - 0xdc00));
                    }
                }
                high
            }
            _ => return Err(refuse("invalid escape in a string")),
        })
    }

    /// The four hex digits `offset` bytes ahead.
    fn hex4(&mut self, offset: usize) -> io::Result<u32> {
        let mut value = 0;
        for i in offset..offset + 4 {
            let digit = match self.at(i)? {
                Some(b @ b'0'..=b'9') => b - b'0',
                Some(b @ b'a'..=b'f') => b - b'a' + 10,
                Some(b @ b'A'..=b'F') => b - b'A' + 10,
                _ => return Err(refuse("invalid \\u escape in a string")),
            };
            value = value << 4 | u32::from(digit);
        }
        Ok(value)
    }

    /// Reads a value that is no string, array or object: `null`, `true`, `false`, a number
    /// or one of the three constants. With `text`, appends what Python's `str()` writes
    /// for the value.
    pub fn scalar(&mut self, text: Option<&mut Vec<u8>>) -> io::Result<Skipped> {
        let (word, written, kind, truthy): (&[u8], &[u8], Kind, bool) = match self.value_start()? {
            b'n' => (b"null", b"None", Kind::Null, false),
            b't' => (b"true", b"True", Kind::Bool, true),
            b'f' => (b"false", b"False", Kind::Bool, false),
            b'N' => (b"NaN", b"nan", Kind::Number, true),
            b'I' => (b"Infinity", b"inf", Kind::Number, true),
            b'-' if self.at(1)? == Some(b'I') => (b"-Infinity", b"-inf", Kind::Number, true),
            _ => return self.number(text),
        };
        if !self.starts_with(word)? {
            return Err(refuse("expecting a value"));
        }
        self.pos += word.len();
        if let Some(text) = text {
            text.extend_from_slice(written);
        }
        Ok(Skipped { kind, truthy })
    }

    /// Reads a number: an optional `-`, `0` or digits that do not start with `0`, an
    /// optional fraction, an optional exponent. A fraction or an exponent makes it a float;
    /// an integer is zero only as `0` or `-0`.
    fn number(&mut self, text: Option<&mut Vec<u8>>) -> io::Result<Skipped> {
        let mut len = usize::from(self.at(0)? == Some(b'-'));
        let sign = len;
        match self.at(len)? {
            Some(b'0') => len += 1,
            Some(b'1'..=b'9') => {
                len += 1;
                while matches!(self.at(len)?, Some(b'0'..=b'9')) {
                    len += 1;
                }
            }
            _ => return Err(refuse("expecting a value")),
        }
        let integer_digits = len - sign;
        let mut float = false;
        if self.at(len)? == Some(b'.') && matches!(self.at(len + 1)?, Some(b'0'..=b'9')) {
            float = true;
            len += 2;
            while matches!(self.at(len)?, Some(b'0'..=b'9')) {
                len += 1;
            }
        }
        if matches!(self.at(len)?, Some(b'e' | b'E')) {
            let mut exponent = len + 1;
            if matches!(self.at(exponent)?, Some(b'+' | b'-')) {
                exponent += 1;
            }
            if matches!(self.at(exponent)?, Some(b'0'..=b'9')) {
                while matches!(self.at(exponent)?, Some(b'0'..=b'9')) {
                    exponent += 1;
                }
                float = true;
                len = exponent;
            }
        }
        let token = &self.buf[self.pos..self.pos + len];
        let truthy = if float {
            let value: f64 = std::str::from_utf8(token)
                .ok()
                .and_then(|digits| digits.parse().ok())
                .ok_or_else(|| refuse("unreadable number"))?;
            if let Some(text) = text {
                text.extend_from_slice(py_float_repr(value).as_bytes());
            }
            value != 0.0
        } else {
            if integer_digits > INT_MAX_STR_DIGITS {
                return Err(refuse("integer of more than 4300 digits"));
            }
            let zero = token[sign] == b'0';
            if let Some(text) = text {
                text.extend_from_slice(if zero { b"0" } else { token });
            }
            !zero
        };
        self.pos += len;
        Ok(Skipped {
            kind: Kind::Number,
            truthy,
        })
    }

    /// Reads past the value the scanner stands at, at any depth of nesting.
    pub fn skip_value(&mut self) -> io::Result<Skipped> {
        match self.value_start()? {
            b'"' => Ok(Skipped {
                kind: Kind::Str,
                truthy: self.string(None)?,
            }),
            b'{' | b'[' => self.skip_container(),
            _ => self.scalar(None),
        }
    }

    fn skip_container(&mut self) -> io::Result<Skipped> {
        let kind = if self.buf[self.pos] == b'{' {
            Kind::Object
        } else {
            Kind::Array
        };
        // `nesting` holds, innermost last, whether each open container is an object.
        self.nesting.clear();
        loop {
            // The scanner stands at a value.
            match self.value_start()? {
                opener @ (b'{' | b'[') => {
                    let object = opener == b'{';
                    if self.enter(if object { b'}' } else { b']' })? {
                        self.nesting.push(object);
                        if object {
                            self.member_name(None)?;
                        }
                        continue;
                    }
                    if self.nesting.is_empty() {
                        return Ok(Skipped {
                            kind,
                            truthy: false,
                        });
                    }
                }
                b'"' => {
                    self.string(None)?;
                }
                _ => {
                    self.scalar(None)?;
                }
            }
            // The scanner stands behind a value inside the innermost container.
            loop {
                let object = *self.nesting.last().expect("an open container");
                if self.another(if object { b'}' } else { b']' })? {
                    if object {
                        self.member_name(None)?;
                    }
                    break;
                }
                self.nesting.pop();
                if self.nesting.is_empty() {
                    return Ok(Skipped { kind, truthy: true });
                }
            }
        }
    }
}

/// Appends the UTF-8 of `code_point`; a surrogate takes the three bytes its code point
/// would take.
fn push_code_point(out: &mut Vec<u8>, code_point: u32) {
    let c = code_point;
    match c {
        0..0x80 => out.push(c as u8),
        0x80..0x800 => out.extend_from_slice(&[0xc0 | (c >> 6) as u8, 0x80 | (c & 0x3f) as u8]),
        0x800..0x10000 => out.extend_from_slice(&[
            0xe0 | (c >> 12) as u8,
            0x80 | (c >> 6 & 0x3f) as u8,
            0x80 | (c & 0x3f) as u8,
        ]),
        _ => out.extend_from_slice(&[
            0xf0 | (c >> 18) as u8,
            0x80 | (c >> 12 & 0x3f) as u8,
            0x80 | (c >> 6 & 0x3f) as u8,
            0x80 | (c & 0x3f) as u8,
        ]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reads `text` as one document that is skipped whole, through a window of `capacity`
    /// bytes at the start.
    fn skip(text: &[u8], capacity: usize) -> io::Result<Skipped> {
        let mut scanner = Scanner::over(text, vec![0; capacity]);
        scanner.begin()?;
        let skipped = scanner.skip_value()?;
        scanner.finish()?;
        Ok(skipped)
    }

    fn accepted(text: &str) -> Skipped {
        let whole = skip(text.as_bytes(), 4096).unwrap_or_else(|e| panic!("{text:?}: {e}"));
        for capacity in [1, 2, 3, 5] {
            assert_eq!(
                skip(text.as_bytes(), capacity).unwrap(),
                whole,
                "{text:?} through {capacity} bytes"
            );
        }
        whole
    }

    fn refused(text: &[u8]) {
        for capacity in [1, 3, 4096] {
            assert!(
                skip(text, capacity).is_err(),
                "{:?} through {capacity} bytes",
                String::from_utf8_lossy(text)
            );
        }
    }

    /// A backslash and `u` before four hex digits.
    fn u(hex: &str) -> String {
        format!("{}u{hex}", '\\')
    }

    #[test]
    fn values_python_reads_are_read_with_their_truth() {
        let cases = [
            ("null", Kind::Null, false),
            ("true", Kind::Bool, true),
            ("false", Kind::Bool, false),
            ("NaN", Kind::Number, true),
            ("Infinity", Kind::Number, true),
            ("-Infinity", Kind::Number, true),
            ("0", Kind::Number, false),
            ("-0", Kind::Number, false),
            ("7", Kind::Number, true),
            ("-12", Kind::Number, true),
            ("0.0", Kind::Number, false),
            ("-0.0", Kind::Number, false),
            ("0e0", Kind::Number, false),
            ("0.0e-0", Kind::Number, false),
            ("1e-400", Kind::Number, false),
            ("1e999", Kind::Number, true),
            ("1.5e+3", Kind::Number, true),
            ("1E5", Kind::Number, true),
            ("-1E-2", Kind::Number, true),
            (r#""""#, Kind::Str, false),
            (r#""x""#, Kind::Str, true),
            (r#""\n""#, Kind::Str, true),
            ("\"\u{7f}\u{e9}\u{2028}\u{1F600}\"", Kind::Str, true),
            ("[]", Kind::Array, false),
            ("[ ]", Kind::Array, false),
            ("[[]]", Kind::Array, true),
            ("[null]", Kind::Array, true),
            ("{}", Kind::Object, false),
            ("{\t}", Kind::Object, false),
            (r#"{"a":{}}"#, Kind::Object, true),
            (
                " \t\r\n[1 , 2.5,\"x\" , {\"k\" : [ ] , \"l\":{ }}]\n",
                Kind::Array,
                true,
            ),
            (
                r#"{"a":1,"b":[true,false,null,NaN,-Infinity],"a":{"c":"d"}}"#,
                Kind::Object,
                true,
            ),
        ];
        for (text, kind, truthy) in cases {
            assert_eq!(accepted(text), Skipped { kind, truthy }, "{text:?}");
        }
        accepted(&"9".repeat(4300));
        accepted(&format!("-{}", "9".repeat(4300)));
        accepted(&format!("1{}.0", "0".repeat(5000)));
        accepted(&format!("{}1{}", "[".repeat(5000), "]".repeat(5000)));
    }

    #[test]
    fn text_python_refuses_is_refused() {
        let cases = [
            "",
            "   ",
            "-",
            "-NaN",
            "+1",
            "01",
            "1.",
            ".5",
            "1e",
            "1e+",
            "nul",
            "nullx",
            "truefalse",
            "-Inf",
            "[1,]",
            "[,1]",
            "[1 2]",
            "[1",
            "[",
            r#"{"a":1,}"#,
            r#"{"a" 1}"#,
            r#"{"a":}"#,
            "{a:1}",
            "{'a':1}",
            r#"{"a":1"#,
            "{1:2}",
            "[]x",
            "[] []",
            "[]\u{0c}",
            "\u{0b}[]",
            "\u{feff}[]",
            r#""abc"#,
            "\"a\tb\"",
            "\"a\nb\"",
            r#""\x41""#,
            r#""\"#,
            r#""\u12G4""#,
            r#""\u 123""#,
            r#""\u12""#,
            "[1]]",
        ];
        for text in cases {
            refused(text.as_bytes());
        }
        refused(&"9".repeat(4301).into_bytes());
        refused(&format!("-{}", "9".repeat(4301)).into_bytes());
        refused(b"\"\xff\"");
        refused(b"\"\xe2\x82\"");
        refused(b"\"\xed\xa0\x80\"");
        refused(b"[\"ok\", \"\xc3\"]");
        refused(format!("\"{}{}\"", u("d83d"), u("12G4")).as_bytes());
    }

    fn unescaped(text: &str) -> Vec<u8> {
        let mut whole = None;
        for capacity in [1, 2, 3, 7, 4096] {
            let mut scanner = Scanner::over(text.as_bytes(), vec![0; capacity]);
            scanner.begin().unwrap();
            let mut out = Vec::new();
            scanner.string(Some(&mut out)).unwrap();
            scanner.finish().unwrap();
            assert_eq!(*whole.get_or_insert(out.clone()), out);
        }
        whole.unwrap()
    }

    #[test]
    fn strings_are_unescaped_and_surrogates_paired_as_python_pairs_them() {
        assert_eq!(unescaped(r#""ENST00000457565""#), b"ENST00000457565");
        assert_eq!(
            unescaped(r#""a\"b\\c\/d\be\ff\ng\rh\ti""#),
            b"a\"b\\c/d\x08e\x0cf\ng\rh\ti"
        );
        assert_eq!(
            unescaped(&format!("\"{}{}{}\"", u("0041"), u("00e9"), u("4E2d"))),
            "A\u{e9}\u{4e2d}".as_bytes()
        );
        assert_eq!(
            unescaped(&format!("\"{}{}\"", u("d83d"), u("de00"))),
            "\u{1F600}".as_bytes()
        );
        assert_eq!(
            unescaped(&format!("\"\u{1F600}{}{}\"", u("D83D"), u("DE00"))),
            "\u{1F600}\u{1F600}".as_bytes()
        );
        // A lone high surrogate, a high one before a letter escape, and a low one first.
        assert_eq!(unescaped(&format!("\"{}\"", u("d83d"))), b"\xed\xa0\xbd");
        assert_eq!(
            unescaped(&format!("\"{}{}\"", u("d83d"), u("0041"))),
            b"\xed\xa0\xbdA"
        );
        assert_eq!(
            unescaped(&format!("\"{}{}\"", u("de00"), u("d83d"))),
            b"\xed\xb8\x80\xed\xa0\xbd"
        );
        assert_eq!(
            unescaped(&format!("\"{}{}{}\"", u("d83d"), u("d83d"), u("de00"))),
            b"\xed\xa0\xbd\xf0\x9f\x98\x80"
        );
    }

    fn written(text: &str) -> String {
        let mut scanner = Scanner::over(text.as_bytes(), vec![0; 2]);
        scanner.begin().unwrap();
        let mut out = Vec::new();
        scanner.scalar(Some(&mut out)).unwrap();
        scanner.finish().unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn scalars_are_written_as_pythons_str_writes_them() {
        let cases = [
            ("null", "None"),
            ("true", "True"),
            ("false", "False"),
            ("NaN", "nan"),
            ("Infinity", "inf"),
            ("-Infinity", "-inf"),
            ("12", "12"),
            ("-12", "-12"),
            ("0", "0"),
            ("-0", "0"),
            ("1.5", "1.5"),
            ("1.50", "1.5"),
            ("-0.0", "-0.0"),
            ("100.0", "100.0"),
            ("1e5", "100000.0"),
            ("1e22", "1e+22"),
            ("1.0e-5", "1e-05"),
            ("1e999", "inf"),
            ("-1e999", "-inf"),
            ("1e-400", "0.0"),
            (
                "123456789012345678901234567890",
                "123456789012345678901234567890",
            ),
        ];
        for (text, expected) in cases {
            assert_eq!(written(text), expected, "str({text})");
        }
    }

    #[test]
    fn members_are_walked_and_the_rest_skipped() {
        let text = br#" { "skip" : [1, {"deep": [[], {}]}], "keep": "v", "n": 3 } "#;
        for capacity in [1, 4, 4096] {
            let mut scanner = Scanner::over(&text[..], vec![0; capacity]);
            scanner.begin().unwrap();
            assert_eq!(scanner.value_start().unwrap(), b'{');
            assert!(scanner.enter(b'}').unwrap());
            let mut seen = Vec::new();
            loop {
                let mut name = Vec::new();
                scanner.key(&mut name).unwrap();
                let mut value = Vec::new();
                if scanner.value_start().unwrap() == b'"' {
                    scanner.string(Some(&mut value)).unwrap();
                } else {
                    scanner.skip_value().unwrap();
                }
                seen.push((name, value));
                if !scanner.another(b'}').unwrap() {
                    break;
                }
            }
            scanner.finish().unwrap();
            assert_eq!(
                seen,
                [
                    (b"skip".to_vec(), Vec::new()),
                    (b"keep".to_vec(), b"v".to_vec()),
                    (b"n".to_vec(), Vec::new()),
                ]
            );
        }
    }
}
