// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The JSON value a manifest and a cache record are held in: an object keeps its keys in the
//! order they were written and an integer its digits. `loads` reads a text as `json.loads`
//! reads it and `dumps` writes a value as `json.dumps` writes it.

use crate::pyfmt::{push_json_str, py_float_repr};
use std::cell::Cell;
use std::collections::HashMap;
use std::fmt;
use std::fmt::Write;

/// The nesting `loads` reads: an array or object inside this many others is refused.
pub const MAX_DEPTH: usize = 1000;

/// The digits of the longest integer `loads` reads, the bound past which `int()` raises
/// `ValueError`.
pub const MAX_INT_DIGITS: usize = 4300;

/// The message of the decode error for a `\u` escape without four hex digits, in the words
/// of `json.JSONDecodeError`.
const BAD_UNICODE_ESCAPE: &str = "Invalid \\uXXXX escape";

/// An integer of any length, held as its decimal digits with a `-` before a negative one.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Int(String);

impl Int {
    /// The digits, with a leading `-` for a negative number.
    pub fn text(&self) -> &str {
        &self.0
    }

    /// The number, where 64 bits hold it.
    pub fn as_i64(&self) -> Option<i64> {
        self.0.parse().ok()
    }

    pub fn is_negative(&self) -> bool {
        self.0.starts_with('-')
    }

    pub fn is_zero(&self) -> bool {
        self.0 == "0"
    }
}

impl From<i64> for Int {
    fn from(number: i64) -> Int {
        Int(number.to_string())
    }
}

impl From<usize> for Int {
    fn from(number: usize) -> Int {
        Int(number.to_string())
    }
}

/// A JSON value as `json.loads` returns one: `None`, a `bool`, an `int`, a `float`, a `str`,
/// a `list` or a `dict`. Two values are equal as Python's are, but for numbers: an integer
/// equals an integer of the same digits and never a float or a boolean.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Null,
    Bool(bool),
    Int(Int),
    Float(f64),
    Str(String),
    List(Vec<Value>),
    Object(Object),
}

impl Value {
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(text) => Some(text),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Value::Bool(flag) => Some(*flag),
            _ => None,
        }
    }

    /// The number of an integer that 64 bits hold; `None` for a float and a boolean.
    pub fn as_i64(&self) -> Option<i64> {
        match self {
            Value::Int(int) => int.as_i64(),
            _ => None,
        }
    }

    pub fn as_list(&self) -> Option<&[Value]> {
        match self {
            Value::List(items) => Some(items),
            _ => None,
        }
    }

    pub fn as_list_mut(&mut self) -> Option<&mut Vec<Value>> {
        match self {
            Value::List(items) => Some(items),
            _ => None,
        }
    }

    pub fn as_object(&self) -> Option<&Object> {
        match self {
            Value::Object(object) => Some(object),
            _ => None,
        }
    }

    pub fn as_object_mut(&mut self) -> Option<&mut Object> {
        match self {
            Value::Object(object) => Some(object),
            _ => None,
        }
    }

    /// The value under `key` of an object; `None` for a key it lacks and for any other value.
    pub fn get(&self, key: &str) -> Option<&Value> {
        self.as_object().and_then(|object| object.get(key))
    }

    /// Python's truth of the value: false for `None`, `False`, zero, the empty string, the
    /// empty list and the empty dict.
    pub fn truthy(&self) -> bool {
        match self {
            Value::Null => false,
            Value::Bool(flag) => *flag,
            Value::Int(int) => !int.is_zero(),
            Value::Float(x) => *x != 0.0,
            Value::Str(text) => !text.is_empty(),
            Value::List(items) => !items.is_empty(),
            Value::Object(object) => !object.is_empty(),
        }
    }

    /// A list of strings from `names`.
    pub fn strings<S: AsRef<str>>(names: &[S]) -> Value {
        Value::List(
            names
                .iter()
                .map(|name| Value::from(name.as_ref()))
                .collect(),
        )
    }
}

impl From<&str> for Value {
    fn from(text: &str) -> Value {
        Value::Str(text.to_string())
    }
}

impl From<String> for Value {
    fn from(text: String) -> Value {
        Value::Str(text)
    }
}

impl From<bool> for Value {
    fn from(flag: bool) -> Value {
        Value::Bool(flag)
    }
}

impl From<i64> for Value {
    fn from(number: i64) -> Value {
        Value::Int(Int::from(number))
    }
}

impl From<usize> for Value {
    fn from(number: usize) -> Value {
        Value::Int(Int::from(number))
    }
}

impl From<Vec<Value>> for Value {
    fn from(items: Vec<Value>) -> Value {
        Value::List(items)
    }
}

impl From<Object> for Value {
    fn from(object: Object) -> Value {
        Value::Object(object)
    }
}

/// The entries past which an object keeps an index of its keys.
const INDEXED_FROM: usize = 16;

/// A JSON object with a `dict`'s behaviour: the keys stay in the order they were first set,
/// setting a key the object holds replaces its value in place, and two objects are equal
/// when they hold the same keys with equal values, in any order.
#[derive(Clone, Debug, Default)]
pub struct Object {
    entries: Vec<(String, Value)>,
    index: Option<HashMap<String, usize>>,
}

impl Object {
    pub fn new() -> Object {
        Object::default()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    fn position(&self, key: &str) -> Option<usize> {
        match &self.index {
            Some(index) => index.get(key).copied(),
            None => self.entries.iter().position(|(held, _)| held == key),
        }
    }

    fn reindex(&mut self) {
        self.index = (self.entries.len() >= INDEXED_FROM).then(|| {
            self.entries
                .iter()
                .enumerate()
                .map(|(place, (key, _))| (key.clone(), place))
                .collect()
        });
    }

    pub fn contains_key(&self, key: &str) -> bool {
        self.position(key).is_some()
    }

    /// `d.get(key)`.
    pub fn get(&self, key: &str) -> Option<&Value> {
        self.position(key).map(|place| &self.entries[place].1)
    }

    pub fn get_mut(&mut self, key: &str) -> Option<&mut Value> {
        self.position(key).map(|place| &mut self.entries[place].1)
    }

    /// `d[key] = value`: a key the object holds keeps its place and takes the value, whose
    /// predecessor is returned; another key goes last.
    pub fn insert(&mut self, key: impl Into<String>, value: impl Into<Value>) -> Option<Value> {
        let key = key.into();
        let value = value.into();
        if let Some(place) = self.position(&key) {
            return Some(std::mem::replace(&mut self.entries[place].1, value));
        }
        if let Some(index) = &mut self.index {
            index.insert(key.clone(), self.entries.len());
        }
        self.entries.push((key, value));
        if self.index.is_none() && self.entries.len() >= INDEXED_FROM {
            self.reindex();
        }
        None
    }

    /// `d.setdefault(key, default)`: the value under `key`, `default` set there first when
    /// the object lacks the key.
    pub fn get_or_insert(&mut self, key: &str, default: impl Into<Value>) -> &mut Value {
        if !self.contains_key(key) {
            self.insert(key, default);
        }
        self.get_mut(key).expect("the key was set")
    }

    /// `d.pop(key, None)`: the value leaves and the other keys keep their order.
    pub fn remove(&mut self, key: &str) -> Option<Value> {
        let place = self.position(key)?;
        let (_, value) = self.entries.remove(place);
        self.reindex();
        Some(value)
    }

    /// Keeps the entries `keep` answers true for, in their order.
    pub fn retain(&mut self, mut keep: impl FnMut(&str, &Value) -> bool) {
        self.entries.retain(|(key, value)| keep(key, value));
        self.reindex();
    }

    /// `d.update(other)`: every entry of `other` set in its order.
    pub fn update(&mut self, other: Object) {
        for (key, value) in other.entries {
            self.insert(key, value);
        }
    }

    /// The entries in the object's order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &Value)> {
        self.entries
            .iter()
            .map(|(key, value)| (key.as_str(), value))
    }

    /// The entries in the object's order, each value open to change.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (&str, &mut Value)> {
        self.entries
            .iter_mut()
            .map(|(key, value)| (key.as_str(), value))
    }

    pub fn keys(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|(key, _)| key.as_str())
    }

    pub fn values(&self) -> impl Iterator<Item = &Value> {
        self.entries.iter().map(|(_, value)| value)
    }
}

impl PartialEq for Object {
    fn eq(&self, other: &Object) -> bool {
        self.len() == other.len()
            && self
                .iter()
                .all(|(key, value)| other.get(key) == Some(value))
    }
}

impl<K: Into<String>, V: Into<Value>> FromIterator<(K, V)> for Object {
    fn from_iter<I: IntoIterator<Item = (K, V)>>(entries: I) -> Object {
        let mut object = Object::new();
        for (key, value) in entries {
            object.insert(key, value);
        }
        object
    }
}

/// Why `loads` refuses a text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum JsonError {
    /// The `json.JSONDecodeError` `json.loads` raises: its message, and the place it names as
    /// a line, a column and a character index, each counted in code points.
    Decode {
        message: &'static str,
        line: usize,
        column: usize,
        position: usize,
    },
    /// A text outside what this reader holds, which `json.loads` reads or ends on another
    /// exception for: a lone surrogate escape, nesting past `MAX_DEPTH`, an integer of more
    /// than `MAX_INT_DIGITS` digits.
    Unsupported(String),
}

impl fmt::Display for JsonError {
    /// `str(exc)` of a `json.JSONDecodeError`; the detail of an unsupported text.
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            JsonError::Decode {
                message,
                line,
                column,
                position,
            } => write!(
                out,
                "{message}: line {line} column {column} (char {position})"
            ),
            JsonError::Unsupported(detail) => out.write_str(detail),
        }
    }
}

/// `json.loads(text)`: any JSON value between blanks, `NaN`, `Infinity` and `-Infinity`
/// among the numbers, a repeated key keeping its first place and its last value, a control
/// character inside a string refused.
pub fn loads(text: &str) -> Result<Value, JsonError> {
    let reader = Reader {
        text,
        bytes: text.as_bytes(),
        lone_surrogate: Cell::new(None),
    };
    if text.starts_with('\u{feff}') {
        return Err(reader.error("Unexpected UTF-8 BOM (decode using utf-8-sig)", 0));
    }
    let (value, end) = reader.value(reader.blanks(0))?;
    let end = reader.blanks(end);
    if end != text.len() {
        return Err(reader.error("Extra data", end));
    }
    match reader.lone_surrogate.get() {
        Some(at) => Err(reader.unsupported("a lone surrogate escape", at)),
        None => Ok(value),
    }
}

struct Reader<'a> {
    text: &'a str,
    bytes: &'a [u8],
    /// The offset of the first surrogate escape that pairs with no other.
    lone_surrogate: Cell<Option<usize>>,
}

impl Reader<'_> {
    /// The decode error `message` at the byte offset `at`, which is the start of a character.
    fn error(&self, message: &'static str, at: usize) -> JsonError {
        let before = &self.text[..at];
        let position = before.chars().count();
        let line = before.bytes().filter(|byte| *byte == b'\n').count() + 1;
        let column = match before.rfind('\n') {
            Some(newline) => before[newline..].chars().count(),
            None => position + 1,
        };
        JsonError::Decode {
            message,
            line,
            column,
            position,
        }
    }

    fn unsupported(&self, what: &str, at: usize) -> JsonError {
        let position = self.text[..at].chars().count();
        JsonError::Unsupported(format!("{what} at char {position}"))
    }

    /// The offset of the first byte from `at` that is not a space, tab, newline or return.
    fn blanks(&self, mut at: usize) -> usize {
        while matches!(self.bytes.get(at), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            at += 1;
        }
        at
    }

    /// The value that starts at `at` and the offset after it. The arrays and objects that
    /// are open are held in a list, so that the nesting of a text costs no stack.
    fn value(&self, mut at: usize) -> Result<(Value, usize), JsonError> {
        let bytes = self.bytes;
        let mut open: Vec<Open> = Vec::new();
        loop {
            let rest = bytes.get(at..).unwrap_or_default();
            let (mut value, mut end) = match rest.first() {
                Some(b'{' | b'[') if open.len() >= MAX_DEPTH => {
                    return Err(self.unsupported(&format!("nesting past {MAX_DEPTH} levels"), at))
                }
                Some(b'[') => {
                    at = self.blanks(at + 1);
                    if bytes.get(at) != Some(&b']') {
                        open.push(Open::List(Vec::new()));
                        continue;
                    }
                    (Value::List(Vec::new()), at + 1)
                }
                Some(b'{') => {
                    at = self.blanks(at + 1);
                    if bytes.get(at) != Some(&b'}') {
                        let (key, value_at) = self.key(at)?;
                        open.push(Open::Object(Object::new(), key));
                        at = value_at;
                        continue;
                    }
                    (Value::Object(Object::new()), at + 1)
                }
                _ => self.scalar(at)?,
            };
            // The value goes to the array or object it stands in, which closes here or
            // takes a further value after a comma.
            loop {
                let Some(parent) = open.last_mut() else {
                    return Ok((value, end));
                };
                let (close, trailing) = match parent {
                    Open::List(items) => {
                        items.push(value);
                        (b']', "Illegal trailing comma before end of array")
                    }
                    Open::Object(object, key) => {
                        object.insert(std::mem::take(key), value);
                        (b'}', "Illegal trailing comma before end of object")
                    }
                };
                let after = self.blanks(end);
                match bytes.get(after) {
                    Some(byte) if *byte == close => {
                        value = match open.pop().expect("the parent") {
                            Open::List(items) => Value::List(items),
                            Open::Object(object, _) => Value::Object(object),
                        };
                        end = after + 1;
                    }
                    Some(b',') => {
                        at = self.blanks(after + 1);
                        if bytes.get(at) == Some(&close) {
                            return Err(self.error(trailing, after));
                        }
                        if let Some(Open::Object(_, key)) = open.last_mut() {
                            (*key, at) = self.key(at)?;
                        }
                        break;
                    }
                    _ => return Err(self.error("Expecting ',' delimiter", after)),
                }
            }
        }
    }

    /// The key at `at` and the offset of its value, past the `:` between them.
    fn key(&self, at: usize) -> Result<(String, usize), JsonError> {
        if self.bytes.get(at) != Some(&b'"') {
            return Err(self.error("Expecting property name enclosed in double quotes", at));
        }
        let (key, after) = self.string(at + 1)?;
        let colon = self.blanks(after);
        if self.bytes.get(colon) != Some(&b':') {
            return Err(self.error("Expecting ':' delimiter", colon));
        }
        Ok((key, self.blanks(colon + 1)))
    }

    /// The string, constant or number that starts at `at` and the offset after it.
    fn scalar(&self, at: usize) -> Result<(Value, usize), JsonError> {
        let rest = self.bytes.get(at..).unwrap_or_default();
        match rest.first() {
            None => Err(self.error("Expecting value", at)),
            Some(b'"') => {
                let (text, end) = self.string(at + 1)?;
                Ok((Value::Str(text), end))
            }
            Some(b'n') if rest.starts_with(b"null") => Ok((Value::Null, at + 4)),
            Some(b't') if rest.starts_with(b"true") => Ok((Value::Bool(true), at + 4)),
            Some(b'f') if rest.starts_with(b"false") => Ok((Value::Bool(false), at + 5)),
            Some(b'N') if rest.starts_with(b"NaN") => Ok((Value::Float(f64::NAN), at + 3)),
            Some(b'I') if rest.starts_with(b"Infinity") => {
                Ok((Value::Float(f64::INFINITY), at + 8))
            }
            Some(b'-') if rest.starts_with(b"-Infinity") => {
                Ok((Value::Float(f64::NEG_INFINITY), at + 9))
            }
            Some(_) => self.number(at),
        }
    }

    /// A number as the scanner matches one: an optional `-`, `0` or digits that start with
    /// no `0`, a fraction when a digit follows the point, and an exponent when a digit ends
    /// it; an exponent without one is left unread.
    fn number(&self, start: usize) -> Result<(Value, usize), JsonError> {
        let bytes = self.bytes;
        let digit = |at: usize| bytes.get(at).is_some_and(u8::is_ascii_digit);
        let mut at = start;
        if bytes[at] == b'-' {
            at += 1;
        }
        match bytes.get(at) {
            Some(b'1'..=b'9') => {
                while digit(at) {
                    at += 1;
                }
            }
            Some(b'0') => at += 1,
            _ => return Err(self.error("Expecting value", start)),
        }
        let mut float = false;
        if bytes.get(at) == Some(&b'.') && digit(at + 1) {
            float = true;
            at += 2;
            while digit(at) {
                at += 1;
            }
        }
        if matches!(bytes.get(at), Some(b'e' | b'E')) && at + 1 < bytes.len() {
            let exponent = at;
            at += 1;
            if matches!(bytes.get(at), Some(b'-' | b'+')) && at + 1 < bytes.len() {
                at += 1;
            }
            while digit(at) {
                at += 1;
            }
            if bytes[at - 1].is_ascii_digit() {
                float = true;
            } else {
                at = exponent;
            }
        }
        let literal = &self.text[start..at];
        if float {
            let number = literal.parse().expect("a float literal");
            return Ok((Value::Float(number), at));
        }
        if literal.trim_start_matches('-').len() > MAX_INT_DIGITS {
            return Err(self.unsupported(
                &format!("an integer of more than {MAX_INT_DIGITS} digits"),
                start,
            ));
        }
        let digits = if literal == "-0" { "0" } else { literal };
        Ok((Value::Int(Int(digits.to_string())), at))
    }

    /// The string whose opening quote is before `start`, and the offset after its closing
    /// quote.
    fn string(&self, start: usize) -> Result<(String, usize), JsonError> {
        let bytes = self.bytes;
        let begin = start - 1;
        let mut text = String::new();
        let mut end = start;
        loop {
            let mut next = end;
            let stop = loop {
                match bytes.get(next) {
                    None => return Err(self.error("Unterminated string starting at", begin)),
                    Some(byte @ (b'"' | b'\\')) => break *byte,
                    Some(byte) if *byte <= 0x1f => {
                        return Err(self.error("Invalid control character at", next))
                    }
                    Some(_) => next += 1,
                }
            };
            text.push_str(&self.text[end..next]);
            next += 1;
            if stop == b'"' {
                return Ok((text, next));
            }
            let Some(escape) = bytes.get(next) else {
                return Err(self.error("Unterminated string starting at", begin));
            };
            if *escape != b'u' {
                end = next + 1;
                text.push(match escape {
                    b'"' => '"',
                    b'\\' => '\\',
                    b'/' => '/',
                    b'b' => '\u{08}',
                    b'f' => '\u{0c}',
                    b'n' => '\n',
                    b'r' => '\r',
                    b't' => '\t',
                    _ => return Err(self.error("Invalid \\escape", next - 1)),
                });
                continue;
            }
            // Four hex digits, or the escape is refused at its `u`.
            let mut unit = self
                .hex4(next + 1)
                .ok_or_else(|| self.error(BAD_UNICODE_ESCAPE, next))?;
            end = next + 5;
            // A high surrogate takes the low surrogate escape behind it; a second escape that
            // is no low surrogate is read on its own.
            let paired = (0xD800..=0xDBFF).contains(&unit)
                && bytes.get(end) == Some(&b'\\')
                && bytes.get(end + 1) == Some(&b'u');
            if paired {
                let low = self
                    .hex4(end + 2)
                    .ok_or_else(|| self.error(BAD_UNICODE_ESCAPE, end + 1))?;
                if (0xDC00..=0xDFFF).contains(&low) {
                    unit = 0x10000 + (((unit - 0xD800) << 10) | (low - 0xDC00));
                    end += 6;
                }
            }
            match char::from_u32(unit) {
                Some(character) => text.push(character),
                None => {
                    // A decode error later in the text comes before this refusal.
                    if self.lone_surrogate.get().is_none() {
                        self.lone_surrogate.set(Some(next - 1));
                    }
                    text.push(char::REPLACEMENT_CHARACTER);
                }
            }
        }
    }

    /// The number four hex digits at `at` spell.
    fn hex4(&self, at: usize) -> Option<u32> {
        let digits = self.bytes.get(at..at + 4)?;
        digits.iter().try_fold(0u32, |number, digit| {
            Some(number << 4 | char::from(*digit).to_digit(16)?)
        })
    }
}

/// An array or an object the reader is inside: its items, or its entries and the key
/// whose value is being read.
enum Open {
    List(Vec<Value>),
    Object(Object, String),
}

/// How `dumps` lays a value out: the arguments of `json.dumps` that change its text.
#[derive(Clone, Copy, Debug)]
pub struct Style {
    /// The spaces of one indent level; `None` writes the value on one line.
    pub indent: Option<usize>,
    pub item_separator: &'static str,
    pub key_separator: &'static str,
    /// Whether every character outside ASCII is written as a `\u` escape.
    pub ensure_ascii: bool,
}

impl Style {
    /// `json.dumps(value)`.
    pub const PLAIN: Style = Style {
        indent: None,
        item_separator: ", ",
        key_separator: ": ",
        ensure_ascii: true,
    };

    /// `json.dumps(value, indent=1)`, the manifest's layout.
    pub const INDENT_1: Style = Style {
        indent: Some(1),
        item_separator: ",",
        key_separator: ": ",
        ensure_ascii: true,
    };

    /// `json.dumps(value, separators=(",", ":"), ensure_ascii=False)`.
    pub const COMPACT_UNICODE: Style = Style {
        indent: None,
        item_separator: ",",
        key_separator: ":",
        ensure_ascii: false,
    };
}

/// `json.dumps(value, ...)` under `style`: an integer's digits, a float as
/// `float.__repr__` with `NaN`, `Infinity` and `-Infinity`, an empty list as `[]` and an
/// empty dict as `{}`.
pub fn dumps(value: &Value, style: &Style) -> String {
    let mut out = String::new();
    write_value(&mut out, value, style, 0);
    out
}

/// `dumps` of an object, without a `Value` around it.
pub fn dumps_object(object: &Object, style: &Style) -> String {
    let mut out = String::new();
    write_object(&mut out, object, style, 0);
    out
}

fn write_string(out: &mut String, text: &str, style: &Style) {
    if style.ensure_ascii {
        push_json_str(out, text);
        return;
    }
    out.push('"');
    for character in text.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            '\0'..='\u{1f}' => write!(out, "\\u{:04x}", character as u32).unwrap(),
            _ => out.push(character),
        }
    }
    out.push('"');
}

fn write_break(out: &mut String, style: &Style, level: usize) {
    if let Some(indent) = style.indent {
        out.push('\n');
        out.extend(std::iter::repeat_n(' ', indent * level));
    }
}

fn write_value(out: &mut String, value: &Value, style: &Style, level: usize) {
    match value {
        Value::Null => out.push_str("null"),
        Value::Bool(true) => out.push_str("true"),
        Value::Bool(false) => out.push_str("false"),
        Value::Int(int) => out.push_str(int.text()),
        Value::Float(x) if x.is_nan() => out.push_str("NaN"),
        Value::Float(x) if x.is_infinite() => {
            out.push_str(if *x > 0.0 { "Infinity" } else { "-Infinity" })
        }
        Value::Float(x) => out.push_str(&py_float_repr(*x)),
        Value::Str(text) => write_string(out, text, style),
        Value::List(items) if items.is_empty() => out.push_str("[]"),
        Value::List(items) => {
            out.push('[');
            for (place, item) in items.iter().enumerate() {
                if place > 0 {
                    out.push_str(style.item_separator);
                }
                write_break(out, style, level + 1);
                write_value(out, item, style, level + 1);
            }
            write_break(out, style, level);
            out.push(']');
        }
        Value::Object(object) => write_object(out, object, style, level),
    }
}

fn write_object(out: &mut String, object: &Object, style: &Style, level: usize) {
    if object.is_empty() {
        out.push_str("{}");
        return;
    }
    out.push('{');
    for (place, (key, item)) in object.iter().enumerate() {
        if place > 0 {
            out.push_str(style.item_separator);
        }
        write_break(out, style, level + 1);
        write_string(out, key, style);
        out.push_str(style.key_separator);
        write_value(out, item, style, level + 1);
    }
    write_break(out, style, level);
    out.push('}');
}

#[cfg(test)]
#[path = "json_reader_tests.rs"]
mod reader_tests;

#[cfg(test)]
#[path = "json_writer_tests.rs"]
mod writer_tests;
