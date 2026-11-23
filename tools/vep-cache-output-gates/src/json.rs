// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A JSON document read by the grammar of Python's `json.load`: strict JSON, plus the
//! constants `NaN`, `Infinity` and `-Infinity`, plus `\u` escapes that spell a lone
//! surrogate; no byte order mark, no control character inside a string, only space, tab,
//! line feed and carriage return between tokens, one value a document, and no integer of
//! more than 4,300 digits. The text is UTF-8 and refused where it is not.
//!
//! A string is unescaped into bytes that are UTF-8 but for a lone surrogate, which is
//! written in three bytes as UTF-8 would write its code point; such bytes keep the order
//! and the equality of the code points.

use std::collections::HashMap;

/// The nesting of a document: an array or object inside this many others is refused.
pub const MAX_DEPTH: usize = 1000;

/// The integer digits Python's `int` converts from text.
const INT_MAX_STR_DIGITS: usize = 4300;

pub type Res<T> = Result<T, String>;

/// A JSON value as `json.load` returns one.
#[derive(Clone, Debug)]
pub enum Value {
    Null,
    Bool(bool),
    /// An integer that 128 bits hold.
    Int(i128),
    /// A longer integer: its digits as the document writes them, a `-` before a negative
    /// number.
    BigInt(Box<str>),
    Float(f64),
    Str(Box<[u8]>),
    List(Vec<Value>),
    Object(Object),
}

/// A JSON object with a `dict`'s behaviour: a repeated key keeps its first place and takes
/// its last value.
#[derive(Clone, Debug, Default)]
pub struct Object {
    entries: Vec<(Box<[u8]>, Value)>,
    index: HashMap<Box<[u8]>, usize>,
}

impl Object {
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// `d.get(key)`.
    pub fn get(&self, key: &[u8]) -> Option<&Value> {
        self.index.get(key).map(|place| &self.entries[*place].1)
    }

    /// `d[key] = value`.
    pub fn insert(&mut self, key: &[u8], value: Value) {
        if let Some(place) = self.index.get(key) {
            self.entries[*place].1 = value;
            return;
        }
        self.index.insert(key.into(), self.entries.len());
        self.entries.push((key.into(), value));
    }

    /// The entries in the object's order.
    pub fn iter(&self) -> impl Iterator<Item = (&[u8], &Value)> {
        self.entries.iter().map(|(key, value)| (&**key, value))
    }

    /// The keys in the object's order.
    pub fn keys(&self) -> impl Iterator<Item = &[u8]> {
        self.entries.iter().map(|(key, _)| &**key)
    }
}

/// An object of `entries` in their order, for a report that is written.
pub fn object(entries: Vec<(&str, Value)>) -> Value {
    let mut out = Object::default();
    for (key, value) in entries {
        out.insert(key.as_bytes(), value);
    }
    Value::Object(out)
}

/// A string value of `text`.
pub fn string(text: &str) -> Value {
    Value::Str(text.as_bytes().into())
}

/// `bool(value)`.
pub fn truthy(value: &Value) -> bool {
    match value {
        Value::Null => false,
        Value::Bool(flag) => *flag,
        Value::Int(n) => *n != 0,
        Value::BigInt(_) => true,
        Value::Float(x) => *x != 0.0,
        Value::Str(text) => !text.is_empty(),
        Value::List(items) => !items.is_empty(),
        Value::Object(object) => !object.is_empty(),
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Reader<'_> {
    fn refuse<T>(&self, what: &str) -> Res<T> {
        Err(format!("{what} at byte {}", self.pos))
    }

    fn at(&self, offset: usize) -> Option<u8> {
        self.bytes.get(self.pos + offset).copied()
    }

    fn peek(&self) -> Option<u8> {
        self.at(0)
    }

    fn starts_with(&self, text: &[u8]) -> bool {
        self.bytes[self.pos..].starts_with(text)
    }

    fn skip_whitespace(&mut self) {
        while let Some(b' ' | b'\t' | b'\n' | b'\r') = self.peek() {
            self.pos += 1;
        }
    }

    fn value(&mut self, depth: usize) -> Res<Value> {
        match self.peek() {
            Some(b'"') => Ok(Value::Str(self.string()?.into())),
            Some(b'{') => self.object(depth),
            Some(b'[') => self.list(depth),
            Some(b'n') => self.word(b"null", Value::Null),
            Some(b't') => self.word(b"true", Value::Bool(true)),
            Some(b'f') => self.word(b"false", Value::Bool(false)),
            Some(b'N') => self.word(b"NaN", Value::Float(f64::NAN)),
            Some(b'I') => self.word(b"Infinity", Value::Float(f64::INFINITY)),
            Some(b'-') if self.at(1) == Some(b'I') => {
                self.word(b"-Infinity", Value::Float(f64::NEG_INFINITY))
            }
            _ => self.number(),
        }
    }

    fn word(&mut self, word: &[u8], value: Value) -> Res<Value> {
        if !self.starts_with(word) {
            return self.refuse("expecting a value");
        }
        self.pos += word.len();
        Ok(value)
    }

    fn deeper(&self, depth: usize) -> Res<usize> {
        if depth >= MAX_DEPTH {
            return self.refuse("a value nested past 1000 levels");
        }
        Ok(depth + 1)
    }

    fn list(&mut self, depth: usize) -> Res<Value> {
        let inside = self.deeper(depth)?;
        self.pos += 1;
        self.skip_whitespace();
        let mut items = Vec::new();
        if self.peek() == Some(b']') {
            self.pos += 1;
            return Ok(Value::List(items));
        }
        loop {
            items.push(self.value(inside)?);
            self.skip_whitespace();
            match self.peek() {
                Some(b',') => {
                    self.pos += 1;
                    self.skip_whitespace();
                }
                Some(b']') => {
                    self.pos += 1;
                    return Ok(Value::List(items));
                }
                _ => return self.refuse("expecting ',' or ']'"),
            }
        }
    }

    fn object(&mut self, depth: usize) -> Res<Value> {
        let inside = self.deeper(depth)?;
        self.pos += 1;
        self.skip_whitespace();
        let mut object = Object::default();
        if self.peek() == Some(b'}') {
            self.pos += 1;
            return Ok(Value::Object(object));
        }
        loop {
            if self.peek() != Some(b'"') {
                return self.refuse("expecting a member name in double quotes");
            }
            let name = self.string()?;
            self.skip_whitespace();
            if self.peek() != Some(b':') {
                return self.refuse("expecting ':' after a member name");
            }
            self.pos += 1;
            self.skip_whitespace();
            let value = self.value(inside)?;
            object.insert(&name, value);
            self.skip_whitespace();
            match self.peek() {
                Some(b',') => {
                    self.pos += 1;
                    self.skip_whitespace();
                }
                Some(b'}') => {
                    self.pos += 1;
                    return Ok(Value::Object(object));
                }
                _ => return self.refuse("expecting ',' or '}'"),
            }
        }
    }

    /// Reads the string whose opening quote the reader stands at, unescaped.
    fn string(&mut self) -> Res<Vec<u8>> {
        self.pos += 1;
        let mut out = Vec::new();
        loop {
            let rest = &self.bytes[self.pos..];
            let Some(found) = rest.iter().position(|b| *b == b'"' || *b == b'\\') else {
                return self.refuse("unterminated string");
            };
            let run = &rest[..found];
            if run.iter().any(|b| *b < 0x20) {
                return self.refuse("control character in a string");
            }
            out.extend_from_slice(run);
            let quote = rest[found] == b'"';
            self.pos += found + 1;
            if quote {
                return Ok(out);
            }
            let code_point = self.escape()?;
            push_code_point(&mut out, code_point);
        }
    }

    /// The code point of the escape whose backslash was read. A high surrogate joins a low
    /// surrogate escape right behind it; any other surrogate stands alone.
    fn escape(&mut self) -> Res<u32> {
        let Some(letter) = self.peek() else {
            return self.refuse("unterminated string");
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
                if (0xd800..0xdc00).contains(&high) && self.starts_with(b"\\u") {
                    let low = self.hex4(2)?;
                    if (0xdc00..0xe000).contains(&low) {
                        self.pos += 6;
                        return Ok(0x10000 + ((high - 0xd800) << 10) + (low - 0xdc00));
                    }
                }
                high
            }
            _ => return self.refuse("invalid escape in a string"),
        })
    }

    /// The four hex digits `offset` bytes ahead.
    fn hex4(&self, offset: usize) -> Res<u32> {
        let mut value = 0;
        for i in offset..offset + 4 {
            let digit = match self.at(i) {
                Some(b @ b'0'..=b'9') => b - b'0',
                Some(b @ b'a'..=b'f') => b - b'a' + 10,
                Some(b @ b'A'..=b'F') => b - b'A' + 10,
                _ => return self.refuse("invalid \\u escape in a string"),
            };
            value = value << 4 | u32::from(digit);
        }
        Ok(value)
    }

    /// Reads the number the reader stands at: an optional `-`, `0` or digits that do not
    /// start with `0`, an optional fraction, an optional exponent. A fraction or an
    /// exponent makes it a float.
    fn number(&mut self) -> Res<Value> {
        let digit = |at: usize| matches!(self.at(at), Some(b'0'..=b'9'));
        let mut len = usize::from(self.peek() == Some(b'-'));
        let sign = len;
        match self.at(len) {
            Some(b'0') => len += 1,
            Some(b'1'..=b'9') => {
                len += 1;
                while digit(len) {
                    len += 1;
                }
            }
            _ => return self.refuse("expecting a value"),
        }
        let digits = len - sign;
        let mut float = false;
        if self.at(len) == Some(b'.') && digit(len + 1) {
            float = true;
            len += 2;
            while digit(len) {
                len += 1;
            }
        }
        if matches!(self.at(len), Some(b'e' | b'E')) {
            let mut exponent = len + 1;
            if matches!(self.at(exponent), Some(b'+' | b'-')) {
                exponent += 1;
            }
            if digit(exponent) {
                while digit(exponent) {
                    exponent += 1;
                }
                float = true;
                len = exponent;
            }
        }
        if !float && digits > INT_MAX_STR_DIGITS {
            return self.refuse("integer of more than 4300 digits");
        }
        let text = std::str::from_utf8(&self.bytes[self.pos..self.pos + len])
            .expect("a number token is ASCII");
        let value = if float {
            Value::Float(text.parse().expect("a float literal"))
        } else {
            match text.parse::<i128>() {
                Ok(n) => Value::Int(n),
                Err(_) => Value::BigInt(text.into()),
            }
        };
        self.pos += len;
        Ok(value)
    }
}

/// Reads `bytes` as one document.
pub fn load(bytes: &[u8]) -> Res<Value> {
    if let Err(fault) = std::str::from_utf8(bytes) {
        return Err(format!("not UTF-8 at byte {}", fault.valid_up_to()));
    }
    let mut reader = Reader { bytes, pos: 0 };
    if reader.starts_with(b"\xef\xbb\xbf") {
        return reader.refuse("byte order mark");
    }
    reader.skip_whitespace();
    let value = reader.value(0)?;
    reader.skip_whitespace();
    match reader.peek() {
        None => Ok(value),
        Some(_) => reader.refuse("extra data after the value"),
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

/// The code points of a string as [`load`] unescapes one.
pub fn code_points(text: &[u8]) -> Vec<u32> {
    let mut out = Vec::with_capacity(text.len());
    let mut at = 0;
    while at < text.len() {
        let lead = u32::from(text[at]);
        let (len, mut code_point) = match lead {
            0..0x80 => (1, lead),
            0x80..0xe0 => (2, lead & 0x1f),
            0xe0..0xf0 => (3, lead & 0x0f),
            _ => (4, lead & 0x07),
        };
        for byte in &text[at + 1..at + len] {
            code_point = code_point << 6 | u32::from(byte & 0x3f);
        }
        out.push(code_point);
        at += len;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text_of(value: &Value) -> String {
        match value {
            Value::Null => "null".to_string(),
            Value::Bool(flag) => flag.to_string(),
            Value::Int(n) => format!("int {n}"),
            Value::BigInt(digits) => format!("big {digits}"),
            Value::Float(x) => format!("float {x:?}"),
            Value::Str(text) => format!("str {:x?}", code_points(text)),
            Value::List(items) => {
                let inner: Vec<String> = items.iter().map(text_of).collect();
                format!("[{}]", inner.join(", "))
            }
            Value::Object(object) => {
                let inner: Vec<String> = object
                    .iter()
                    .map(|(key, value)| format!("{:x?}: {}", code_points(key), text_of(value)))
                    .collect();
                format!("{{{}}}", inner.join(", "))
            }
        }
    }

    fn read(text: &str) -> String {
        match load(text.as_bytes()) {
            Ok(value) => text_of(&value),
            Err(detail) => format!("refused: {detail}"),
        }
    }

    #[test]
    fn scalars_are_read_as_json_load_reads_them() {
        for (text, expected) in [
            ("null", "null"),
            (" true\n", "true"),
            ("false", "false"),
            ("0", "int 0"),
            ("-0", "int 0"),
            ("-17", "int -17"),
            (
                "170141183460469231731687303715884105727",
                "int 170141183460469231731687303715884105727",
            ),
            (
                "170141183460469231731687303715884105728",
                "big 170141183460469231731687303715884105728",
            ),
            (
                "-170141183460469231731687303715884105729",
                "big -170141183460469231731687303715884105729",
            ),
            ("1.5", "float 1.5"),
            ("-0.0", "float -0.0"),
            ("1e3", "float 1000.0"),
            ("1E-2", "float 0.01"),
            ("1e400", "float inf"),
            ("NaN", "float NaN"),
            ("Infinity", "float inf"),
            ("-Infinity", "float -inf"),
            ("\"\"", "str []"),
        ] {
            assert_eq!(read(text), expected, "{text}");
        }
    }

    #[test]
    fn what_json_load_refuses_is_refused() {
        for text in [
            "",
            " ",
            "01",
            "-",
            "1.",
            ".5",
            "1e",
            "+1",
            "nul",
            "True",
            "nan",
            "-NaN",
            "[1,]",
            "[1 2]",
            "{\"a\":1,}",
            "{a: 1}",
            "{\"a\" 1}",
            "\"a",
            "\"a\tb\"",
            "\"\\x41\"",
            "\"\\u12g4\"",
            "\"\\ud800\\u12\"",
            "1 2",
            "\u{feff}1",
            "'a'",
        ] {
            assert!(read(text).starts_with("refused: "), "{text:?}");
        }
        assert!(load(b"\"\xff\"").is_err());
        assert!(load(b"1 \xff").is_err());
        let long_integer = "9".repeat(4301);
        assert!(load(long_integer.as_bytes()).is_err());
        assert!(load("9".repeat(4300).as_bytes()).is_ok());
        assert!(load(format!("1.{}", "9".repeat(5000)).as_bytes()).is_ok());
        let nested = |levels: usize| format!("{}{}", "[".repeat(levels), "]".repeat(levels));
        assert!(load(nested(1000).as_bytes()).is_ok());
        assert!(load(nested(1001).as_bytes()).is_err());
    }

    #[test]
    fn strings_are_unescaped_with_surrogates_as_python_keeps_them() {
        assert_eq!(
            read(r#""a\"\\\/\b\f\n\r\tz""#),
            "str [61, 22, 5c, 2f, 8, c, a, d, 9, 7a]"
        );
        assert_eq!(read("\"\u{e9}\u{1f600}\""), "str [e9, 1f600]");
        assert_eq!(read(r#""\u00e9\ud83d\ude00""#), "str [e9, 1f600]");
        assert_eq!(read(r#""\ud83d""#), "str [d83d]");
        assert_eq!(read(r#""\ud83dx""#), "str [d83d, 78]");
        assert_eq!(read(r#""\ude00\ud83d""#), "str [de00, d83d]");
        assert_eq!(read(r#""\ud83d\u0041""#), "str [d83d, 41]");
    }

    #[test]
    fn a_repeated_key_keeps_its_first_place_and_its_last_value() {
        assert_eq!(
            read(r#"{"b": 1, "a": [], "b": {"c": null}}"#),
            "{[62]: {[63]: null}, [61]: []}"
        );
        let Value::Object(object) = load(br#"{"b": 1, "a": 2, "b": 3}"#).unwrap() else {
            panic!("an object");
        };
        assert!(matches!(object.get(b"b"), Some(Value::Int(3))));
        assert!(object.get(b"c").is_none());
        assert_eq!(object.keys().collect::<Vec<_>>(), [b"b", b"a"]);
    }

    #[test]
    fn truth_is_pythons() {
        for (text, expected) in [
            ("null", false),
            ("false", false),
            ("true", true),
            ("0", false),
            ("-0.0", false),
            ("0.0", false),
            ("2", true),
            ("NaN", true),
            ("\"\"", false),
            ("\"a\"", true),
            ("[]", false),
            ("[0]", true),
            ("{}", false),
            ("{\"a\": 0}", true),
            ("170141183460469231731687303715884105728", true),
        ] {
            assert_eq!(truthy(&load(text.as_bytes()).unwrap()), expected, "{text}");
        }
    }
}
