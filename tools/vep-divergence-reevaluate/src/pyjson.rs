// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! JSON as CPython's `json` module reads and writes it: `json.loads` of one document, with
//! `NaN`, `Infinity` and `-Infinity` read as numbers, integers at any length and a repeated
//! key keeping its first place and its last value; and the two `json.dumps` layouts the
//! tool writes.

use crate::pyfmt::{push_json_str, py_float_repr};
use crate::value::Value;

/// A document deeper than this is refused; CPython stops near a thousand levels.
const MAX_DEPTH: usize = 900;

struct Reader<'a> {
    chars: &'a [char],
    pos: usize,
}

impl Reader<'_> {
    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn skip_whitespace(&mut self) {
        while matches!(self.peek(), Some(' ' | '\t' | '\n' | '\r')) {
            self.pos += 1;
        }
    }

    fn literal(&mut self, word: &str) -> bool {
        let end = self.pos + word.len();
        if end <= self.chars.len() && self.chars[self.pos..end].iter().copied().eq(word.chars()) {
            self.pos = end;
            return true;
        }
        false
    }

    fn value(&mut self, depth: usize) -> Option<Value> {
        if depth > MAX_DEPTH {
            return None;
        }
        match self.peek()? {
            '"' => {
                self.pos += 1;
                self.string().map(Value::Str)
            }
            '{' => {
                self.pos += 1;
                self.object(depth)
            }
            '[' => {
                self.pos += 1;
                self.array(depth)
            }
            'n' => self.literal("null").then_some(Value::None),
            't' => self.literal("true").then_some(Value::Bool(true)),
            'f' => self.literal("false").then_some(Value::Bool(false)),
            'N' => self.literal("NaN").then_some(Value::Float(f64::NAN)),
            'I' => self
                .literal("Infinity")
                .then_some(Value::Float(f64::INFINITY)),
            '-' if self.literal("-Infinity") => Some(Value::Float(f64::NEG_INFINITY)),
            _ => self.number(),
        }
    }

    fn digits(&mut self) -> usize {
        let start = self.pos;
        while matches!(self.peek(), Some('0'..='9')) {
            self.pos += 1;
        }
        self.pos - start
    }

    fn number(&mut self) -> Option<Value> {
        let start = self.pos;
        if self.peek() == Some('-') {
            self.pos += 1;
        }
        match self.peek() {
            Some('1'..='9') => {
                self.digits();
            }
            Some('0') => self.pos += 1,
            _ => {
                self.pos = start;
                return None;
            }
        }
        let mut is_float = false;
        if self.peek() == Some('.') && matches!(self.chars.get(self.pos + 1), Some('0'..='9')) {
            is_float = true;
            self.pos += 1;
            self.digits();
        }
        if matches!(self.peek(), Some('e' | 'E')) {
            let exponent = self.pos;
            self.pos += 1;
            if matches!(self.peek(), Some('-' | '+')) {
                self.pos += 1;
            }
            if self.digits() > 0 {
                is_float = true;
            } else {
                self.pos = exponent;
            }
        }
        let text: String = self.chars[start..self.pos].iter().collect();
        if is_float {
            return text.parse().ok().map(Value::Float);
        }
        Some(Value::Int(if text == "-0" {
            "0".to_string()
        } else {
            text
        }))
    }

    fn hex4(&mut self) -> Option<u32> {
        let end = self.pos + 4;
        let digits = self.chars.get(self.pos..end)?;
        let mut code = 0;
        for digit in digits {
            code = code * 16 + digit.to_digit(16)?;
        }
        self.pos = end;
        Some(code)
    }

    /// The string after its opening quote.
    fn string(&mut self) -> Option<String> {
        let mut out = String::new();
        loop {
            let c = self.peek()?;
            self.pos += 1;
            match c {
                '"' => return Some(out),
                '\\' => {
                    let escape = self.peek()?;
                    self.pos += 1;
                    match escape {
                        '"' => out.push('"'),
                        '\\' => out.push('\\'),
                        '/' => out.push('/'),
                        'b' => out.push('\u{08}'),
                        'f' => out.push('\u{0c}'),
                        'n' => out.push('\n'),
                        'r' => out.push('\r'),
                        't' => out.push('\t'),
                        'u' => {
                            let mut code = self.hex4()?;
                            if (0xd800..0xdc00).contains(&code)
                                && self.chars.get(self.pos) == Some(&'\\')
                                && self.chars.get(self.pos + 1) == Some(&'u')
                            {
                                let back = self.pos;
                                self.pos += 2;
                                match self.hex4() {
                                    Some(low) if (0xdc00..0xe000).contains(&low) => {
                                        code = 0x10000 + ((code - 0xd800) << 10) + (low - 0xdc00);
                                    }
                                    _ => self.pos = back,
                                }
                            }
                            // Half of a surrogate pair on its own is no character of a
                            // text this tool holds: it reads as U+FFFD.
                            out.push(char::from_u32(code).unwrap_or('\u{fffd}'));
                        }
                        _ => return None,
                    }
                }
                c if (c as u32) < 0x20 => return None,
                c => out.push(c),
            }
        }
    }

    /// The object after its opening brace.
    fn object(&mut self, depth: usize) -> Option<Value> {
        let mut items: Vec<(String, Value)> = Vec::new();
        self.skip_whitespace();
        if self.peek() == Some('}') {
            self.pos += 1;
            return Some(Value::Dict(items));
        }
        loop {
            if self.peek() != Some('"') {
                return None;
            }
            self.pos += 1;
            let key = self.string()?;
            self.skip_whitespace();
            if self.peek() != Some(':') {
                return None;
            }
            self.pos += 1;
            self.skip_whitespace();
            let value = self.value(depth + 1)?;
            match items.iter_mut().find(|(k, _)| *k == key) {
                Some(item) => item.1 = value,
                None => items.push((key, value)),
            }
            self.skip_whitespace();
            match self.peek()? {
                '}' => {
                    self.pos += 1;
                    return Some(Value::Dict(items));
                }
                ',' => {
                    self.pos += 1;
                    self.skip_whitespace();
                }
                _ => return None,
            }
        }
    }

    /// The array after its opening bracket.
    fn array(&mut self, depth: usize) -> Option<Value> {
        let mut items = Vec::new();
        self.skip_whitespace();
        if self.peek() == Some(']') {
            self.pos += 1;
            return Some(Value::List(items));
        }
        loop {
            items.push(self.value(depth + 1)?);
            self.skip_whitespace();
            match self.peek()? {
                ']' => {
                    self.pos += 1;
                    return Some(Value::List(items));
                }
                ',' => {
                    self.pos += 1;
                    self.skip_whitespace();
                }
                _ => return None,
            }
        }
    }
}

/// `json.loads(text)`; None where it raises `JSONDecodeError`.
pub fn loads(text: &str) -> Option<Value> {
    let chars: Vec<char> = text.chars().collect();
    let mut reader = Reader {
        chars: &chars,
        pos: 0,
    };
    reader.skip_whitespace();
    let value = reader.value(0)?;
    reader.skip_whitespace();
    (reader.pos == chars.len()).then_some(value)
}

fn push_float(out: &mut String, x: f64) {
    if x.is_nan() {
        out.push_str("NaN");
    } else if x.is_infinite() {
        out.push_str(if x > 0.0 { "Infinity" } else { "-Infinity" });
    } else {
        out.push_str(&py_float_repr(x));
    }
}

/// `json.dumps(value, sort_keys=True, separators=(",", ":"))`.
pub fn dumps_compact_sorted(value: &Value) -> String {
    let mut out = String::new();
    push_compact(&mut out, value);
    out
}

fn push_compact(out: &mut String, value: &Value) {
    match value {
        Value::None => out.push_str("null"),
        Value::Bool(true) => out.push_str("true"),
        Value::Bool(false) => out.push_str("false"),
        Value::Int(digits) => out.push_str(digits),
        Value::Float(x) => push_float(out, *x),
        Value::Str(s) => push_json_str(out, s),
        Value::List(items) => {
            out.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                push_compact(out, item);
            }
            out.push(']');
        }
        Value::Dict(items) => {
            let mut sorted: Vec<&(String, Value)> = items.iter().collect();
            sorted.sort_by(|a, b| a.0.cmp(&b.0));
            out.push('{');
            for (i, (key, item)) in sorted.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                push_json_str(out, key);
                out.push(':');
                push_compact(out, item);
            }
            out.push('}');
        }
        // No JSON document holds a date or a time.
        Value::DateTime(_) | Value::Date(_) | Value::Time(_) => push_json_str(out, &value.py_str()),
    }
}

/// `json.dumps(value, indent=indent)`: an item a line, each level `indent` blanks deeper,
/// an empty list as `[]` and an empty object as `{}`.
pub fn dumps_indented(value: &Value, indent: usize) -> String {
    let mut out = String::new();
    push_indented(&mut out, value, indent, 0);
    out
}

fn push_indented(out: &mut String, value: &Value, indent: usize, level: usize) {
    let newline = |out: &mut String, level: usize| {
        out.push('\n');
        out.push_str(&" ".repeat(indent * level));
    };
    match value {
        Value::List(items) if !items.is_empty() => {
            out.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                newline(out, level + 1);
                push_indented(out, item, indent, level + 1);
            }
            newline(out, level);
            out.push(']');
        }
        Value::Dict(items) if !items.is_empty() => {
            out.push('{');
            for (i, (key, item)) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                newline(out, level + 1);
                push_json_str(out, key);
                out.push_str(": ");
                push_indented(out, item, indent, level + 1);
            }
            newline(out, level);
            out.push('}');
        }
        Value::List(_) => out.push_str("[]"),
        Value::Dict(_) => out.push_str("{}"),
        scalar => push_compact(out, scalar),
    }
}

#[cfg(test)]
#[path = "pyjson_tests.rs"]
mod tests;
