// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A JSON document read by the grammar of Python's `json.load`: strict JSON, plus the
//! constants `NaN`, `Infinity` and `-Infinity`, plus `\u` escapes that spell a lone
//! surrogate; no byte order mark, no control character inside a string, only space, tab,
//! line feed and carriage return between tokens, one value a document, and no integer of
//! more than 4,300 digits. The text is UTF-8 and refused where it is not.
//!
//! A caller walks the values it reads and skips the rest, and every skipped value is still
//! held to the grammar. A value that is kept is a [`Value`], compared and counted true as
//! Python compares and counts what `json.load` returns. A string is unescaped into bytes
//! that are UTF-8 but for a lone surrogate, which is written in three bytes as UTF-8 would
//! write its code point; such bytes keep the order and the equality of the code points.

use memchr::memchr2;
use std::collections::HashMap;

/// The nesting of a value that is kept: an array or object inside this many others is
/// refused. A value that is skipped nests to any depth.
pub const MAX_DEPTH: usize = 1000;

/// The integer digits Python's `int` converts from text.
const INT_MAX_STR_DIGITS: usize = 4300;

/// The entries past which an object keeps an index of its keys.
const INDEXED_FROM: usize = 32;

pub type Res<T> = Result<T, String>;

/// An integer of any length: the number where 64 bits hold it, its digits where not.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Int {
    Small(i64),
    /// The digits as the document writes them, a `-` before a negative number.
    Big(Box<str>),
}

/// A JSON value as `json.load` returns one.
#[derive(Clone, Debug)]
pub enum Value {
    Null,
    Bool(bool),
    Int(Int),
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
    index: Option<HashMap<Box<[u8]>, usize>>,
}

impl Object {
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    fn position(&self, key: &[u8]) -> Option<usize> {
        match &self.index {
            Some(index) => index.get(key).copied(),
            None => self.entries.iter().position(|(held, _)| **held == *key),
        }
    }

    /// `d.get(key)`.
    pub fn get(&self, key: &[u8]) -> Option<&Value> {
        self.position(key).map(|place| &self.entries[place].1)
    }

    /// `d[key] = value`.
    pub fn insert(&mut self, key: &[u8], value: Value) {
        if let Some(place) = self.position(key) {
            self.entries[place].1 = value;
            return;
        }
        if let Some(index) = &mut self.index {
            index.insert(key.into(), self.entries.len());
        }
        self.entries.push((key.into(), value));
        if self.index.is_none() && self.entries.len() >= INDEXED_FROM {
            self.index = Some(
                self.entries
                    .iter()
                    .enumerate()
                    .map(|(place, (key, _))| (key.clone(), place))
                    .collect(),
            );
        }
    }

    /// The entries in the object's order.
    pub fn iter(&self) -> impl Iterator<Item = (&[u8], &Value)> {
        self.entries.iter().map(|(key, value)| (&**key, value))
    }
}

/// A number as Python compares it: `True` is 1 and `False` is 0.
enum Number<'a> {
    Small(i64),
    Big(&'a str),
    Float(f64),
}

fn number(value: &Value) -> Option<Number<'_>> {
    match value {
        Value::Bool(flag) => Some(Number::Small(i64::from(*flag))),
        Value::Int(Int::Small(n)) => Some(Number::Small(*n)),
        Value::Int(Int::Big(digits)) => Some(Number::Big(digits)),
        Value::Float(x) => Some(Number::Float(*x)),
        _ => None,
    }
}

/// Whether an integer equals a float, compared exactly as Python compares them.
fn small_equals_float(n: i64, x: f64) -> bool {
    x.fract() == 0.0 && (x as i128) == i128::from(n)
}

fn big_equals_float(digits: &str, x: f64) -> bool {
    x.fract() == 0.0 && format!("{x:.0}") == digits
}

fn numbers_equal(a: &Number, b: &Number) -> bool {
    match (a, b) {
        (Number::Small(x), Number::Small(y)) => x == y,
        (Number::Big(x), Number::Big(y)) => x == y,
        (Number::Small(_), Number::Big(_)) | (Number::Big(_), Number::Small(_)) => false,
        // Every NaN `json.load` returns is one object, and a list or dict takes an object
        // as equal to itself.
        (Number::Float(x), Number::Float(y)) => x == y || (x.is_nan() && y.is_nan()),
        (Number::Small(n), Number::Float(x)) | (Number::Float(x), Number::Small(n)) => {
            small_equals_float(*n, *x)
        }
        (Number::Big(digits), Number::Float(x)) | (Number::Float(x), Number::Big(digits)) => {
            big_equals_float(digits, *x)
        }
    }
}

/// `a == b` as a list or a dict asks it of two of its items: numbers by value across
/// `bool`, `int` and `float`, a dict by its keys in any order, `NaN` equal to `NaN`.
pub fn py_eq(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Null, Value::Null) => true,
        (Value::Str(x), Value::Str(y)) => x == y,
        (Value::List(x), Value::List(y)) => {
            x.len() == y.len() && x.iter().zip(y).all(|(p, q)| py_eq(p, q))
        }
        (Value::Object(x), Value::Object(y)) => {
            x.len() == y.len()
                && x.iter()
                    .all(|(key, value)| y.get(key).is_some_and(|other| py_eq(value, other)))
        }
        _ => match (number(a), number(b)) {
            (Some(x), Some(y)) => numbers_equal(&x, &y),
            _ => false,
        },
    }
}

/// `a != b` as the operator answers for two values: [`py_eq`] negated, but `NaN` differs
/// from `NaN` where the two are the operands themselves.
pub fn py_ne(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Float(x), Value::Float(y)) if x.is_nan() && y.is_nan() => true,
        _ => !py_eq(a, b),
    }
}

/// Python's truth of a value: false for `None`, `False`, zero, the empty string, the empty
/// list and the empty dict.
pub fn truthy(value: &Value) -> bool {
    match value {
        Value::Null => false,
        Value::Bool(flag) => *flag,
        Value::Int(Int::Small(n)) => *n != 0,
        Value::Int(Int::Big(_)) => true,
        Value::Float(x) => *x != 0.0,
        Value::Str(text) => !text.is_empty(),
        Value::List(items) => !items.is_empty(),
        Value::Object(object) => !object.is_empty(),
    }
}

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

/// What is known of a value that was read past: its type, whether Python counts it true,
/// and `len()` of an array (its items) or an object (its distinct keys).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Shape {
    pub kind: Kind,
    pub truthy: bool,
    pub len: usize,
}

impl Shape {
    /// The shape of `None`, which `dict.get` returns for a key the dict lacks.
    pub const ABSENT: Shape = Shape {
        kind: Kind::Null,
        truthy: false,
        len: 0,
    };

    /// `len(value) if isinstance(value, list) else 0`.
    pub fn list_len(&self) -> usize {
        match self.kind {
            Kind::Array => self.len,
            _ => 0,
        }
    }

    /// `len(value) if isinstance(value, (list, dict)) else 0`.
    pub fn container_len(&self) -> usize {
        match self.kind {
            Kind::Array | Kind::Object => self.len,
            _ => 0,
        }
    }
}

/// A number token: its length and whether a fraction or an exponent makes it a float.
struct Token {
    len: usize,
    float: bool,
}

pub struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
    nesting: Vec<bool>,
}

impl<'a> Reader<'a> {
    pub fn new(bytes: &'a [u8]) -> Reader<'a> {
        Reader {
            bytes,
            pos: 0,
            nesting: Vec::new(),
        }
    }

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

    /// Starts the document: refuses a byte order mark and stops at the value's first byte.
    pub fn begin(&mut self) -> Res<()> {
        if self.starts_with(b"\xef\xbb\xbf") {
            return self.refuse("byte order mark");
        }
        self.skip_whitespace();
        Ok(())
    }

    /// Ends the document: nothing but whitespace may follow its value.
    pub fn finish(&mut self) -> Res<()> {
        self.skip_whitespace();
        match self.peek() {
            None => Ok(()),
            Some(_) => self.refuse("extra data after the value"),
        }
    }

    /// The first byte of the value the reader stands at.
    pub fn value_start(&self) -> Res<u8> {
        match self.peek() {
            Some(byte) => Ok(byte),
            None => self.refuse("expecting a value"),
        }
    }

    /// Enters the array or object whose opening byte the reader stands at; false when it
    /// is empty, and then it is closed as well.
    pub fn enter(&mut self, closer: u8) -> Res<bool> {
        self.pos += 1;
        self.skip_whitespace();
        if self.peek() == Some(closer) {
            self.pos += 1;
            return Ok(false);
        }
        Ok(true)
    }

    /// After a member of an array or object: true at the next member, false when the
    /// closing byte ended the container.
    pub fn another(&mut self, closer: u8) -> Res<bool> {
        self.skip_whitespace();
        match self.peek() {
            Some(b',') => {
                self.pos += 1;
                self.skip_whitespace();
                Ok(true)
            }
            Some(c) if c == closer => {
                self.pos += 1;
                Ok(false)
            }
            _ => self.refuse("expecting ',' or the end of the container"),
        }
    }

    /// Reads an object member's name and the `:` behind it, and stops at the member's
    /// value. The name is appended to `name`.
    pub fn key(&mut self, name: &mut Vec<u8>) -> Res<()> {
        self.member_name(Some(name))
    }

    fn member_name(&mut self, name: Option<&mut Vec<u8>>) -> Res<()> {
        if self.peek() != Some(b'"') {
            return self.refuse("expecting a member name in double quotes");
        }
        self.string(name)?;
        self.skip_whitespace();
        if self.peek() != Some(b':') {
            return self.refuse("expecting ':' after a member name");
        }
        self.pos += 1;
        self.skip_whitespace();
        Ok(())
    }

    /// Reads the string whose opening quote the reader stands at, appending its unescaped
    /// content to `out`; true when the string is not empty.
    pub fn string(&mut self, mut out: Option<&mut Vec<u8>>) -> Res<bool> {
        self.pos += 1;
        let mut filled = false;
        loop {
            let rest = &self.bytes[self.pos..];
            let Some(found) = memchr2(b'"', b'\\', rest) else {
                return self.refuse("unterminated string");
            };
            let run = &rest[..found];
            // A run ends at an ASCII byte, so it is UTF-8 exactly where the document is.
            let control = run.iter().fold(false, |any, b| any | (*b < 0x20));
            if control || (!run.is_ascii() && std::str::from_utf8(run).is_err()) {
                return self.refuse("control character or invalid UTF-8 in a string");
            }
            filled |= !run.is_empty();
            if let Some(out) = out.as_deref_mut() {
                out.extend_from_slice(run);
            }
            let quote = rest[found] == b'"';
            self.pos += found + 1;
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

    /// The constant the reader stands at: its text and its value.
    fn constant(&self) -> Option<(&'static [u8], Value)> {
        Some(match self.peek()? {
            b'n' => (b"null", Value::Null),
            b't' => (b"true", Value::Bool(true)),
            b'f' => (b"false", Value::Bool(false)),
            b'N' => (b"NaN", Value::Float(f64::NAN)),
            b'I' => (b"Infinity", Value::Float(f64::INFINITY)),
            b'-' if self.at(1) == Some(b'I') => (b"-Infinity", Value::Float(f64::NEG_INFINITY)),
            _ => return None,
        })
    }

    /// Measures the number the reader stands at: an optional `-`, `0` or digits that do
    /// not start with `0`, an optional fraction, an optional exponent.
    fn number(&self) -> Res<Token> {
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
        Ok(Token { len, float })
    }

    /// Reads a value that is no string, array or object: `null`, `true`, `false`, a number
    /// or one of the three constants.
    fn scalar(&mut self) -> Res<Value> {
        if let Some((word, value)) = self.constant() {
            if !self.starts_with(word) {
                return self.refuse("expecting a value");
            }
            self.pos += word.len();
            return Ok(value);
        }
        let token = self.number()?;
        let text = std::str::from_utf8(&self.bytes[self.pos..self.pos + token.len])
            .expect("a number token is ASCII");
        let value = if token.float {
            Value::Float(text.parse().expect("a float literal"))
        } else {
            Value::Int(match text.parse::<i64>() {
                Ok(n) => Int::Small(n),
                Err(_) => Int::Big(text.into()),
            })
        };
        self.pos += token.len;
        Ok(value)
    }

    /// Reads past a scalar without its value.
    fn skip_scalar(&mut self) -> Res<()> {
        if let Some((word, _)) = self.constant() {
            if !self.starts_with(word) {
                return self.refuse("expecting a value");
            }
            self.pos += word.len();
            return Ok(());
        }
        self.pos += self.number()?.len;
        Ok(())
    }

    /// Reads past the value the reader stands at, at any depth of nesting.
    pub fn skip_value(&mut self) -> Res<()> {
        match self.value_start()? {
            b'"' => self.string(None).map(|_| ()),
            b'{' | b'[' => self.skip_container(),
            _ => self.skip_scalar(),
        }
    }

    fn skip_container(&mut self) -> Res<()> {
        // `nesting` holds, innermost last, whether each open container is an object.
        self.nesting.clear();
        loop {
            // The reader stands at a value.
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
                        return Ok(());
                    }
                }
                b'"' => {
                    self.string(None)?;
                }
                _ => self.skip_scalar()?,
            }
            // The reader stands behind a value inside the innermost container.
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
                    return Ok(());
                }
            }
        }
    }

    /// Reads past the value the reader stands at and returns its shape.
    pub fn shape(&mut self) -> Res<Shape> {
        match self.value_start()? {
            b'"' => {
                let truthy = self.string(None)?;
                Ok(Shape {
                    kind: Kind::Str,
                    truthy,
                    len: 0,
                })
            }
            b'[' => {
                let mut len = 0;
                if self.enter(b']')? {
                    loop {
                        self.skip_value()?;
                        len += 1;
                        if !self.another(b']')? {
                            break;
                        }
                    }
                }
                Ok(Shape {
                    kind: Kind::Array,
                    truthy: len > 0,
                    len,
                })
            }
            b'{' => {
                let mut keys: Vec<Vec<u8>> = Vec::new();
                if self.enter(b'}')? {
                    loop {
                        let mut name = Vec::new();
                        self.key(&mut name)?;
                        self.skip_value()?;
                        keys.push(name);
                        if !self.another(b'}')? {
                            break;
                        }
                    }
                }
                keys.sort_unstable();
                keys.dedup();
                Ok(Shape {
                    kind: Kind::Object,
                    truthy: !keys.is_empty(),
                    len: keys.len(),
                })
            }
            _ => {
                let value = self.scalar()?;
                let kind = match value {
                    Value::Null => Kind::Null,
                    Value::Bool(_) => Kind::Bool,
                    _ => Kind::Number,
                };
                Ok(Shape {
                    kind,
                    truthy: truthy(&value),
                    len: 0,
                })
            }
        }
    }

    /// Reads the value the reader stands at and keeps it. The arrays and objects that are
    /// open are held in a list, so that the nesting of a document costs no stack.
    pub fn value(&mut self) -> Res<Value> {
        let mut open: Vec<Open> = Vec::new();
        loop {
            // The reader stands at a value.
            let mut value = match self.value_start()? {
                b'"' => {
                    let mut text = Vec::new();
                    self.string(Some(&mut text))?;
                    Value::Str(text.into())
                }
                opener @ (b'[' | b'{') => {
                    if open.len() >= MAX_DEPTH {
                        return self.refuse("a kept value nested past 1000 levels");
                    }
                    if opener == b'[' {
                        if self.enter(b']')? {
                            open.push(Open::List(Vec::new()));
                            continue;
                        }
                        Value::List(Vec::new())
                    } else {
                        if self.enter(b'}')? {
                            let mut name = Vec::new();
                            self.key(&mut name)?;
                            open.push(Open::Object(Object::default(), name));
                            continue;
                        }
                        Value::Object(Object::default())
                    }
                }
                _ => self.scalar()?,
            };
            // The value goes to the array or object it stands in, which closes here or
            // takes a further value after a comma.
            loop {
                let closer = match open.last_mut() {
                    None => return Ok(value),
                    Some(Open::List(items)) => {
                        items.push(value);
                        b']'
                    }
                    Some(Open::Object(object, name)) => {
                        object.insert(name, value);
                        b'}'
                    }
                };
                if self.another(closer)? {
                    if let Some(Open::Object(_, name)) = open.last_mut() {
                        name.clear();
                        self.key(name)?;
                    }
                    break;
                }
                value = match open.pop().expect("the open container") {
                    Open::List(items) => Value::List(items),
                    Open::Object(object, _) => Value::Object(object),
                };
            }
        }
    }
}

/// An array or an object the reader is inside: its items, or its members and the name
/// whose value is being read.
enum Open {
    List(Vec<Value>),
    Object(Object, Vec<u8>),
}

/// Reads `bytes` as one document and keeps its value.
pub fn load(bytes: &[u8]) -> Res<Value> {
    let mut reader = Reader::new(bytes);
    reader.begin()?;
    let value = reader.value()?;
    reader.finish()?;
    Ok(value)
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

/// The code points of a string as [`Reader::string`] unescapes one.
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

/// The bytes [`Reader::string`] would hold for `code_points`.
pub fn from_code_points(code_points: &[u32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(code_points.len());
    for code_point in code_points {
        push_code_point(&mut out, *code_point);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn skip(text: &[u8]) -> Res<()> {
        let mut reader = Reader::new(text);
        reader.begin()?;
        reader.skip_value()?;
        reader.finish()
    }

    fn shape(text: &str) -> Shape {
        let mut reader = Reader::new(text.as_bytes());
        reader.begin().unwrap();
        let shape = reader.shape().unwrap_or_else(|e| panic!("{text:?}: {e}"));
        reader.finish().unwrap();
        skip(text.as_bytes()).unwrap_or_else(|e| panic!("{text:?}: {e}"));
        let kept = load(text.as_bytes()).unwrap_or_else(|e| panic!("{text:?}: {e}"));
        assert_eq!(truthy(&kept), shape.truthy, "{text:?}");
        shape
    }

    /// A backslash and `u` before four hex digits.
    fn u(hex: &str) -> String {
        format!("{}u{hex}", '\\')
    }

    #[test]
    fn values_python_reads_are_read_with_their_truth_and_length() {
        let cases = [
            ("null", Kind::Null, false, 0),
            ("true", Kind::Bool, true, 0),
            ("false", Kind::Bool, false, 0),
            ("NaN", Kind::Number, true, 0),
            ("Infinity", Kind::Number, true, 0),
            ("-Infinity", Kind::Number, true, 0),
            ("0", Kind::Number, false, 0),
            ("-0", Kind::Number, false, 0),
            ("7", Kind::Number, true, 0),
            ("-12", Kind::Number, true, 0),
            ("0.0", Kind::Number, false, 0),
            ("-0.0", Kind::Number, false, 0),
            ("0e0", Kind::Number, false, 0),
            ("0.0e-0", Kind::Number, false, 0),
            ("1e-400", Kind::Number, false, 0),
            ("1e999", Kind::Number, true, 0),
            ("1.5e+3", Kind::Number, true, 0),
            ("1E5", Kind::Number, true, 0),
            ("-1E-2", Kind::Number, true, 0),
            ("123456789012345678901234567890", Kind::Number, true, 0),
            (r#""""#, Kind::Str, false, 0),
            (r#""x""#, Kind::Str, true, 0),
            (r#""\n""#, Kind::Str, true, 0),
            ("\"\u{7f}\u{e9}\u{2028}\u{1F600}\"", Kind::Str, true, 0),
            ("[]", Kind::Array, false, 0),
            ("[ ]", Kind::Array, false, 0),
            ("[[]]", Kind::Array, true, 1),
            ("[null, 2, [3, 4]]", Kind::Array, true, 3),
            ("{}", Kind::Object, false, 0),
            ("{\t}", Kind::Object, false, 0),
            (r#"{"a":{}}"#, Kind::Object, true, 1),
            (
                r#"{"a":1,"b":[true,false,null,NaN,-Infinity],"a":{"c":"d"}}"#,
                Kind::Object,
                true,
                2,
            ),
            (
                " \t\r\n[1 , 2.5,\"x\" , {\"k\" : [ ] , \"l\":{ }}]\n",
                Kind::Array,
                true,
                4,
            ),
        ];
        for (text, kind, truthy, len) in cases {
            assert_eq!(shape(text), Shape { kind, truthy, len }, "{text:?}");
        }
        shape(&"9".repeat(4300));
        shape(&format!("-{}", "9".repeat(4300)));
        shape(&format!("1{}.0", "0".repeat(5000)));
        skip(format!("{}1{}", "[".repeat(5000), "]".repeat(5000)).as_bytes()).unwrap();
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
            assert!(skip(text.as_bytes()).is_err(), "skipping {text:?}");
            assert!(load(text.as_bytes()).is_err(), "keeping {text:?}");
        }
        let long = "9".repeat(4301);
        for text in [long.clone(), format!("-{long}"), format!("[{long}]")] {
            assert!(skip(text.as_bytes()).is_err());
            assert!(load(text.as_bytes()).is_err());
        }
        for bytes in [
            &b"\"\xff\""[..],
            b"\"\xe2\x82\"",
            b"\"\xed\xa0\x80\"",
            b"[\"ok\", \"\xc3\"]",
        ] {
            assert!(skip(bytes).is_err());
            assert!(load(bytes).is_err());
        }
        assert!(skip(format!("\"{}{}\"", u("d83d"), u("12G4")).as_bytes()).is_err());
    }

    #[test]
    fn a_kept_value_nests_to_a_thousand_levels_and_a_skipped_one_deeper() {
        let nested = |depth: usize| format!("{}1{}", "[".repeat(depth), "]".repeat(depth));
        assert!(load(nested(MAX_DEPTH).as_bytes()).is_ok());
        assert!(load(nested(MAX_DEPTH + 1).as_bytes()).is_err());
        assert!(skip(nested(MAX_DEPTH + 1).as_bytes()).is_ok());
    }

    fn unescaped(text: &str) -> Vec<u8> {
        match load(text.as_bytes()).unwrap() {
            Value::Str(bytes) => bytes.into_vec(),
            other => panic!("{other:?}"),
        }
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
        let points = [0x41, 0xe9, 0x4e2d, 0xd83d, 0x1f600, 0x10ffff];
        assert_eq!(code_points(&from_code_points(&points)), points);
    }

    fn eq(a: &str, b: &str) -> bool {
        let (a, b) = (load(a.as_bytes()).unwrap(), load(b.as_bytes()).unwrap());
        let equal = py_eq(&a, &b);
        assert_eq!(py_eq(&b, &a), equal);
        equal
    }

    fn ne(a: &str, b: &str) -> bool {
        py_ne(&load(a.as_bytes()).unwrap(), &load(b.as_bytes()).unwrap())
    }

    #[test]
    fn values_compare_as_python_compares_them() {
        for (a, b) in [
            ("null", "null"),
            ("1", "1.0"),
            ("true", "1"),
            ("false", "0.0"),
            ("-0", "0"),
            ("0.0", "-0.0"),
            ("[1, true]", "[1.0, 1]"),
            (r#"{"a": 1, "b": [2]}"#, r#"{"b": [2.0], "a": true}"#),
            (r#"{"a": 1, "a": 2}"#, r#"{"a": 2}"#),
            ("18014398509481984", "18014398509481984.0"),
            ("1e22", "10000000000000000000000"),
            ("-1e22", "-10000000000000000000000"),
            (
                "123456789012345678901234567890",
                "123456789012345678901234567890",
            ),
            ("1e400", "Infinity"),
            ("[NaN]", "[NaN]"),
            (r#"{"a": NaN}"#, r#"{"a": NaN}"#),
            (r#""é""#, "\"\u{e9}\""),
        ] {
            assert!(eq(a, b), "{a} == {b}");
            assert!(!ne(a, b), "not {a} != {b}");
        }
        for (a, b) in [
            ("null", "0"),
            ("null", "false"),
            ("null", "\"\""),
            ("1", "\"1\""),
            ("[]", "{}"),
            ("[1]", "[1, 2]"),
            ("[1, 2]", "[2, 1]"),
            (r#"{"a": 1}"#, r#"{"a": 1, "b": 2}"#),
            (r#"{"a": 1}"#, r#"{"b": 1}"#),
            ("18014398509481985", "18014398509481984.0"),
            ("1e22", "10000000000000000000001"),
            ("1e23", "100000000000000000000000"),
            (
                "123456789012345678901234567890",
                "123456789012345678901234567891",
            ),
            ("9223372036854775807", "9223372036854775808"),
            ("0.5", "0"),
            ("NaN", "1"),
            ("NaN", "Infinity"),
        ] {
            assert!(!eq(a, b), "{a} != {b}");
            assert!(ne(a, b), "{a} != {b}");
        }
        // The operands themselves: NaN differs from NaN, though a list holding it does not.
        assert!(eq("NaN", "NaN"));
        assert!(ne("NaN", "NaN"));
        assert!(!ne("[NaN]", "[NaN]"));
    }

    #[test]
    fn an_object_keeps_first_places_and_last_values() {
        let text = r#"{"b": 1, "a": 2, "b": 3}"#;
        let Value::Object(object) = load(text.as_bytes()).unwrap() else {
            panic!()
        };
        let keys: Vec<&[u8]> = object.iter().map(|(key, _)| key).collect();
        assert_eq!(keys, [b"b", b"a"]);
        assert!(py_eq(object.get(b"b").unwrap(), &Value::Int(Int::Small(3))));
        assert!(object.get(b"c").is_none());
        // Past the size where an index is kept.
        let many: String = (0..100)
            .map(|i| format!("\"k{}\": {i}", i % 60))
            .collect::<Vec<_>>()
            .join(", ");
        let Value::Object(object) = load(format!("{{{many}}}").as_bytes()).unwrap() else {
            panic!()
        };
        assert_eq!(object.len(), 60);
        assert!(py_eq(
            object.get(b"k5").unwrap(),
            &Value::Int(Int::Small(65))
        ));
        assert!(py_eq(
            object.get(b"k59").unwrap(),
            &Value::Int(Int::Small(59))
        ));
    }
}
