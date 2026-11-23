// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! JSON text as `JSON::PP` 4.06 reads it under `decode_json` and writes a string under
//! `encode_json`.
//!
//! The reader takes exactly the texts `decode_json` takes: UTF-8, white space of
//! blank, tab, line feed and carriage return, the literals `true`, `false` and `null`,
//! numbers without a leading zero before another digit, strings without a control
//! character and with every surrogate escape in a pair, no trailing comma, a nesting of
//! at most 512, and nothing after the value. A text it refuses is one `decode_json`
//! dies on.

/// The deepest nesting `decode_json` takes.
const MAX_DEPTH: usize = 512;

/// A value of a JSON text. A string and a key are the place of their text between the
/// quotes; a number is the place of its text.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Null,
    Bool(bool),
    Number(usize, usize),
    Str(usize, usize),
    List(Vec<Value>),
    /// The entries in the order of the text; a key written twice stands twice.
    Dict(Vec<((usize, usize), Value)>),
}

struct Parser<'a> {
    text: &'a [u8],
    at: usize,
    depth: usize,
}

type Parsed<T> = Result<T, &'static str>;

impl Parser<'_> {
    fn peek(&self) -> Option<u8> {
        self.text.get(self.at).copied()
    }

    fn white(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\r' | b'\n')) {
            self.at += 1;
        }
    }

    fn value(&mut self) -> Parsed<Value> {
        self.white();
        match self.peek() {
            None => Err("the text ends where a value is expected"),
            Some(b'{') => self.object(),
            Some(b'[') => self.array(),
            Some(b'"') => {
                let (from, to) = self.string()?;
                Ok(Value::Str(from, to))
            }
            Some(b'0'..=b'9' | b'-') => self.number(),
            Some(_) => self.word(),
        }
    }

    fn word(&mut self) -> Parsed<Value> {
        let rest = &self.text[self.at..];
        let (value, length) = if rest.starts_with(b"true") {
            (Value::Bool(true), 4)
        } else if rest.starts_with(b"null") {
            (Value::Null, 4)
        } else if rest.starts_with(b"false") {
            (Value::Bool(false), 5)
        } else {
            return Err("neither array, object, number, string or atom");
        };
        self.at += length;
        Ok(value)
    }

    fn digits(&mut self) -> usize {
        let from = self.at;
        while matches!(self.peek(), Some(b'0'..=b'9')) {
            self.at += 1;
        }
        self.at - from
    }

    fn number(&mut self) -> Parsed<Value> {
        let from = self.at;
        if self.peek() == Some(b'-') {
            self.at += 1;
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                return Err("no digits after initial minus");
            }
        }
        if self.peek() == Some(b'0')
            && matches!(
                self.text.get(self.at + 1),
                Some(b'0'..=b'9' | b'a'..=b'd' | b'f' | b'A'..=b'D' | b'F')
            )
        {
            return Err("leading zero must not be followed by another digit");
        }
        self.digits();
        if self.peek() == Some(b'.') {
            self.at += 1;
            if self.digits() == 0 {
                return Err("no digits after decimal point");
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.at += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.at += 1;
            }
            if self.digits() == 0 {
                return Err("no digits after exp sign");
            }
        }
        Ok(Value::Number(from, self.at))
    }

    /// A string's text between its quotes, checked.
    fn string(&mut self) -> Parsed<(usize, usize)> {
        if self.peek() != Some(b'"') {
            return Err("unexpected end of string while parsing JSON string");
        }
        self.at += 1;
        let from = self.at;
        // The high half of a surrogate pair that waits for its low half.
        let mut high = false;
        loop {
            match self.peek() {
                None => return Err("unexpected end of string while parsing JSON string"),
                Some(b'"') => {
                    if high {
                        return Err("missing low surrogate character in surrogate pair");
                    }
                    self.at += 1;
                    return Ok((from, self.at - 1));
                }
                Some(b'\\') => match self.text.get(self.at + 1) {
                    Some(b'b' | b't' | b'n' | b'f' | b'r' | b'\\' | b'"' | b'/') => self.at += 2,
                    Some(b'u') => {
                        let Some(digits) = self.text.get(self.at + 2..self.at + 6) else {
                            return Err("unexpected end of string while parsing JSON string");
                        };
                        if !digits.iter().all(u8::is_ascii_hexdigit) {
                            return Err("unexpected end of string while parsing JSON string");
                        }
                        let code = hex(digits);
                        match code {
                            0xd800..=0xdbff => high = true,
                            0xdc00..=0xdfff => {
                                if !high {
                                    return Err(
                                        "missing high surrogate character in surrogate pair",
                                    );
                                }
                                high = false;
                            }
                            _ if high => return Err("surrogate pair expected"),
                            _ => {}
                        }
                        self.at += 6;
                    }
                    _ => return Err("illegal backslash escape sequence in string"),
                },
                Some(byte) if byte < 0x20 => {
                    return Err("invalid character encountered while parsing JSON string")
                }
                Some(_) => self.at += 1,
            }
        }
    }

    fn array(&mut self) -> Parsed<Value> {
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            return Err("json text or perl structure exceeds maximum nesting level");
        }
        self.at += 1;
        self.white();
        let mut items = Vec::new();
        if self.peek() == Some(b']') {
            self.at += 1;
            self.depth -= 1;
            return Ok(Value::List(items));
        }
        loop {
            items.push(self.value()?);
            self.white();
            match self.peek() {
                Some(b']') => {
                    self.at += 1;
                    self.depth -= 1;
                    return Ok(Value::List(items));
                }
                Some(b',') => self.at += 1,
                _ => return Err(", or ] expected while parsing array"),
            }
        }
    }

    fn object(&mut self) -> Parsed<Value> {
        self.depth += 1;
        if self.depth > MAX_DEPTH {
            return Err("json text or perl structure exceeds maximum nesting level");
        }
        self.at += 1;
        self.white();
        let mut entries = Vec::new();
        if self.peek() == Some(b'}') {
            self.at += 1;
            self.depth -= 1;
            return Ok(Value::Dict(entries));
        }
        loop {
            let key = self.string()?;
            self.white();
            if self.peek() != Some(b':') {
                return Err("':' expected");
            }
            self.at += 1;
            entries.push((key, self.value()?));
            self.white();
            match self.peek() {
                Some(b'}') => {
                    self.at += 1;
                    self.depth -= 1;
                    return Ok(Value::Dict(entries));
                }
                Some(b',') => {
                    self.at += 1;
                    self.white();
                }
                _ => return Err(", or } expected while parsing object/hash"),
            }
        }
    }
}

fn hex(digits: &[u8]) -> u32 {
    digits.iter().fold(0, |value, digit| {
        value * 16 + (*digit as char).to_digit(16).unwrap_or(0)
    })
}

/// The value of a whole text; the error is why `decode_json` dies on it.
pub fn decode(text: &[u8]) -> Result<Value, &'static str> {
    if std::str::from_utf8(text).is_err() {
        return Err("malformed UTF-8 character in JSON string");
    }
    let mut parser = Parser {
        text,
        at: 0,
        depth: 0,
    };
    let value = parser.value()?;
    parser.white();
    if parser.at != text.len() {
        return Err("garbage after JSON object");
    }
    Ok(value)
}

/// The characters of a string between its quotes, as UTF-8: its escapes read.
pub fn string_text(raw: &[u8]) -> Vec<u8> {
    if !raw.contains(&b'\\') {
        return raw.to_vec();
    }
    let mut out = Vec::with_capacity(raw.len());
    let mut at = 0;
    while at < raw.len() {
        if raw[at] != b'\\' {
            out.push(raw[at]);
            at += 1;
            continue;
        }
        let simple = match raw[at + 1] {
            b'b' => Some(8u8),
            b't' => Some(9),
            b'n' => Some(10),
            b'f' => Some(12),
            b'r' => Some(13),
            b'u' => None,
            other => Some(other),
        };
        if let Some(byte) = simple {
            out.push(byte);
            at += 2;
            continue;
        }
        let mut code = hex(&raw[at + 2..at + 6]);
        at += 6;
        if (0xd800..0xdc00).contains(&code) {
            // A high half pairs with the low half behind it; one that another high
            // half follows is dropped, as the reader keeps the later of the two.
            let low = match raw.get(at..at + 6) {
                Some([b'\\', b'u', digits @ ..]) => hex(digits),
                _ => 0,
            };
            if !(0xdc00..0xe000).contains(&low) {
                continue;
            }
            code = 0x10000 + ((code - 0xd800) << 10) + (low - 0xdc00);
            at += 6;
        }
        let character = char::from_u32(code).unwrap_or(char::REPLACEMENT_CHARACTER);
        out.extend_from_slice(character.encode_utf8(&mut [0u8; 4]).as_bytes());
    }
    out
}

/// A string as `encode_json` writes it: `"`, `\`, line feed, carriage return, tab, form
/// feed and backspace escaped by letter, every other character below 32 as `\u00xx`,
/// and all else, `/` and DEL included, as itself.
pub fn write_string(characters: &[u8], out: &mut Vec<u8>) {
    out.push(b'"');
    for &byte in characters {
        match byte {
            b'"' => out.extend_from_slice(b"\\\""),
            b'\\' => out.extend_from_slice(b"\\\\"),
            b'\n' => out.extend_from_slice(b"\\n"),
            b'\r' => out.extend_from_slice(b"\\r"),
            b'\t' => out.extend_from_slice(b"\\t"),
            0x0c => out.extend_from_slice(b"\\f"),
            0x08 => out.extend_from_slice(b"\\b"),
            0x00..=0x1f => out.extend_from_slice(format!("\\u00{byte:02x}").as_bytes()),
            _ => out.push(byte),
        }
    }
    out.push(b'"');
}

#[cfg(test)]
#[path = "jsonpp_tests.rs"]
mod tests;
