// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! JSON values read and written as the reference implementation's interpreter does.

use crate::pyfmt;
use indexmap::IndexMap;
use std::fmt;
use std::fmt::Write as _;

/// The deepest nesting of lists and objects `loads` reads: one level more is the
/// interpreter's `RecursionError`. The interpreter's own bound is what is left of its C
/// stack, so it moves with the stack limit and, by a few dozen levels, from one process
/// to the next: CPython 3.14 under a stack limit of 10240 kB first refuses between
/// 81,564 and 81,594 levels, of lists and of objects alike, and this is the middle.
pub const MAX_DEPTH: usize = 81_578;

/// The frames the interpreter allows on its call stack: `json_compact` is a recursive
/// function there (`compact_frames`).
pub const RECURSION_LIMIT: usize = 1000;

/// A JSON value as `json.loads` returns it.
///
/// A value drops its lists and objects without recursion, so a pattern cannot move a
/// field out of an owned value: match on a reference, or take the field with `into_str`,
/// `into_list` or `into_dict`. `Clone` and `==` recurse once a level of nesting.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Null,
    Bool(bool),
    Int(i64),
    /// A whole number beyond `i64`, as its decimal text.
    BigInt(String),
    Float(f64),
    Str(String),
    List(Vec<Value>),
    /// Members in insertion order; a member named twice keeps its first place and takes
    /// the last value.
    Dict(IndexMap<String, Value>),
}

impl Value {
    /// The whole number `decimal` spells, as `int_of_big` of `pyfmt` returns it: `Int`
    /// where it fits `i64`, else `BigInt`.
    pub fn of_int(decimal: &str) -> Value {
        match decimal.parse() {
            Ok(whole) => Value::Int(whole),
            Err(_) => Value::BigInt(decimal.to_string()),
        }
    }

    /// The text of a `str`.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(text) => Some(text),
            _ => None,
        }
    }

    /// The items of a `list`.
    pub fn as_list(&self) -> Option<&[Value]> {
        match self {
            Value::List(items) => Some(items),
            _ => None,
        }
    }

    /// The members of a `dict`.
    pub fn as_dict(&self) -> Option<&IndexMap<String, Value>> {
        match self {
            Value::Dict(members) => Some(members),
            _ => None,
        }
    }

    /// `value.get(name)` of a `dict`; `None` for any other value and for an absent member.
    pub fn get(&self, name: &str) -> Option<&Value> {
        self.as_dict()?.get(name)
    }

    /// The text of a `str`, or the value back.
    pub fn into_str(mut self) -> Result<String, Value> {
        if let Value::Str(text) = &mut self {
            return Ok(std::mem::take(text));
        }
        Err(self)
    }

    /// The items of a `list`, or the value back.
    pub fn into_list(mut self) -> Result<Vec<Value>, Value> {
        if let Value::List(items) = &mut self {
            return Ok(std::mem::take(items));
        }
        Err(self)
    }

    /// The members of a `dict`, or the value back.
    pub fn into_dict(mut self) -> Result<IndexMap<String, Value>, Value> {
        if let Value::Dict(members) = &mut self {
            return Ok(std::mem::take(members));
        }
        Err(self)
    }

    fn holds_values(&self) -> bool {
        match self {
            Value::List(items) => !items.is_empty(),
            Value::Dict(members) => !members.is_empty(),
            _ => false,
        }
    }

    /// Moves the lists and objects directly inside `self` to the values still to free.
    fn release_nested(&mut self, pending: &mut Vec<Value>) {
        let mut release = |item: &mut Value| {
            if item.holds_values() {
                pending.push(std::mem::replace(item, Value::Null));
            }
        };
        match self {
            Value::List(items) => items.iter_mut().for_each(&mut release),
            Value::Dict(members) => members.values_mut().for_each(&mut release),
            _ => {}
        }
    }
}

impl Drop for Value {
    fn drop(&mut self) {
        let mut pending = Vec::new();
        self.release_nested(&mut pending);
        while let Some(mut value) = pending.pop() {
            value.release_nested(&mut pending);
        }
    }
}

/// What `json.loads` raises, and the one text CPython reads that a `String` cannot hold.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JsonErrorKind {
    /// `json.JSONDecodeError`, a `ValueError`.
    Decode,
    /// The `ValueError` of a whole number of more digits than the interpreter converts.
    IntDigits,
    /// `RecursionError`: the text nests deeper than the interpreter reads.
    TooDeep,
    /// The text is JSON, and one of its strings holds a `\u` escape of half a surrogate
    /// pair: CPython keeps the half in the `str`.
    LoneSurrogate,
}

/// A text `loads` does not return a value for; `Display` is `str()` of the exception.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JsonError {
    pub kind: JsonErrorKind,
    pub message: String,
}

impl JsonError {
    /// Whether an `except ValueError` catches it.
    pub fn is_value_error(&self) -> bool {
        matches!(self.kind, JsonErrorKind::Decode | JsonErrorKind::IntDigits)
    }

    /// Whether an `except RecursionError` catches it.
    pub fn is_recursion_error(&self) -> bool {
        self.kind == JsonErrorKind::TooDeep
    }
}

impl fmt::Display for JsonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for JsonError {}

/// The arguments of `json.dumps` the comparator passes; `default()` is `json.dumps(v)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DumpOptions {
    pub indent: Option<usize>,
    pub sort_keys: bool,
    pub ensure_ascii: bool,
}

impl Default for DumpOptions {
    fn default() -> DumpOptions {
        DumpOptions {
            indent: None,
            sort_keys: false,
            ensure_ascii: true,
        }
    }
}

/// A list or an object whose end the reader has not met.
enum Open {
    List(Vec<Value>),
    /// The members read, and the name of the member whose value is being read.
    Dict(IndexMap<String, Value>, String),
}

/// The interpreter's message for a `\u` escape that four hex digits do not follow.
const INVALID_U_ESCAPE: &str = concat!("Invalid \\uXX", "XX escape");

struct Reader<'a> {
    text: &'a str,
    bytes: &'a [u8],
    /// The byte offset of the first escape of half a surrogate pair.
    lone_surrogate: Option<usize>,
    /// The character a string takes for half a surrogate pair.
    surrogate_mark: char,
}

impl Reader<'_> {
    /// The count of code points before byte `at`.
    fn chars_before(&self, at: usize) -> usize {
        self.bytes[..at]
            .iter()
            .filter(|byte| (**byte & 0xc0) != 0x80)
            .count()
    }

    fn error(&self, what: &str, at: usize) -> JsonError {
        let position = self.chars_before(at);
        let before = &self.bytes[..at];
        let line = 1 + before.iter().filter(|byte| **byte == b'\n').count();
        let column = match before.iter().rposition(|byte| *byte == b'\n') {
            Some(newline) => position - self.chars_before(newline),
            None => position + 1,
        };
        JsonError {
            kind: JsonErrorKind::Decode,
            message: format!("{what}: line {line} column {column} (char {position})"),
        }
    }

    fn skip_space(&self, mut at: usize) -> usize {
        while matches!(self.bytes.get(at), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            at += 1;
        }
        at
    }

    /// The value that starts at byte `at`, and the offset after it.
    fn value(&mut self, mut at: usize) -> Result<(Value, usize), JsonError> {
        let bytes = self.bytes;
        let mut open: Vec<Open> = Vec::new();
        loop {
            let rest = &bytes[at.min(bytes.len())..];
            let mut value = match rest.first() {
                None => return Err(self.error("Expecting value", at)),
                Some(b'"') => {
                    let (text, end) = self.string(at + 1)?;
                    at = end;
                    Value::Str(text)
                }
                Some(b'{') => {
                    if open.len() >= MAX_DEPTH {
                        return Err(too_deep("object"));
                    }
                    at = self.skip_space(at + 1);
                    if bytes.get(at) == Some(&b'}') {
                        at += 1;
                        Value::Dict(IndexMap::new())
                    } else {
                        let (name, after) = self.member_name(at)?;
                        at = after;
                        open.push(Open::Dict(IndexMap::new(), name));
                        continue;
                    }
                }
                Some(b'[') => {
                    if open.len() >= MAX_DEPTH {
                        return Err(too_deep("array"));
                    }
                    at = self.skip_space(at + 1);
                    if bytes.get(at) == Some(&b']') {
                        at += 1;
                        Value::List(Vec::new())
                    } else {
                        open.push(Open::List(Vec::new()));
                        continue;
                    }
                }
                Some(b'n') if rest.starts_with(b"null") => {
                    at += 4;
                    Value::Null
                }
                Some(b't') if rest.starts_with(b"true") => {
                    at += 4;
                    Value::Bool(true)
                }
                Some(b'f') if rest.starts_with(b"false") => {
                    at += 5;
                    Value::Bool(false)
                }
                Some(b'N') if rest.starts_with(b"NaN") => {
                    at += 3;
                    Value::Float(f64::NAN)
                }
                Some(b'I') if rest.starts_with(b"Infinity") => {
                    at += 8;
                    Value::Float(f64::INFINITY)
                }
                Some(b'-') if rest.starts_with(b"-Infinity") => {
                    at += 9;
                    Value::Float(f64::NEG_INFINITY)
                }
                Some(_) => {
                    let (number, end) = self.number(at)?;
                    at = end;
                    number
                }
            };
            // The value ends every list and object that closes right after it.
            loop {
                match open.last_mut() {
                    None => return Ok((value, at)),
                    Some(Open::List(items)) => {
                        items.push(value);
                        at = self.skip_space(at);
                        match bytes.get(at) {
                            Some(b']') => {
                                at += 1;
                                value = Value::List(std::mem::take(items));
                                open.pop();
                            }
                            Some(b',') => {
                                let comma = at;
                                at = self.skip_space(at + 1);
                                if bytes.get(at) == Some(&b']') {
                                    return Err(self.error(
                                        "Illegal trailing comma before end of array",
                                        comma,
                                    ));
                                }
                                break;
                            }
                            _ => return Err(self.error("Expecting ',' delimiter", at)),
                        }
                    }
                    Some(Open::Dict(members, name)) => {
                        members.insert(std::mem::take(name), value);
                        at = self.skip_space(at);
                        match bytes.get(at) {
                            Some(b'}') => {
                                at += 1;
                                value = Value::Dict(std::mem::take(members));
                                open.pop();
                            }
                            Some(b',') => {
                                let comma = at;
                                at = self.skip_space(at + 1);
                                if bytes.get(at) == Some(&b'}') {
                                    return Err(self.error(
                                        "Illegal trailing comma before end of object",
                                        comma,
                                    ));
                                }
                                let (next, after) = self.member_name(at)?;
                                let Some(Open::Dict(_, name)) = open.last_mut() else {
                                    unreachable!("the innermost open value is this object");
                                };
                                *name = next;
                                at = after;
                                break;
                            }
                            _ => return Err(self.error("Expecting ',' delimiter", at)),
                        }
                    }
                }
            }
        }
    }

    /// The member name at byte `at`, and the offset of the member's value.
    fn member_name(&mut self, at: usize) -> Result<(String, usize), JsonError> {
        if self.bytes.get(at) != Some(&b'"') {
            return Err(self.error("Expecting property name enclosed in double quotes", at));
        }
        let (name, end) = self.string(at + 1)?;
        let colon = self.skip_space(end);
        if self.bytes.get(colon) != Some(&b':') {
            return Err(self.error("Expecting ':' delimiter", colon));
        }
        Ok((name, self.skip_space(colon + 1)))
    }

    /// The number at byte `start`, and the offset after it.
    fn number(&self, start: usize) -> Result<(Value, usize), JsonError> {
        let bytes = self.bytes;
        let digit = |at: usize| bytes.get(at).is_some_and(u8::is_ascii_digit);
        let mut at = start;
        if bytes[at] == b'-' {
            at += 1;
        }
        let whole = at;
        match bytes.get(at) {
            Some(b'1'..=b'9') => {
                while digit(at) {
                    at += 1;
                }
            }
            Some(b'0') => at += 1,
            _ => return Err(self.error("Expecting value", start)),
        }
        let whole_digits = at - whole;
        let mut float = false;
        if bytes.get(at) == Some(&b'.') && digit(at + 1) {
            float = true;
            at += 2;
            while digit(at) {
                at += 1;
            }
        }
        // An exponent needs a digit; without one the number ends before the `e`.
        if matches!(bytes.get(at), Some(b'e' | b'E')) {
            let mut exponent = at + 1;
            // A sign counts only with a character after it.
            if matches!(bytes.get(exponent), Some(b'+' | b'-')) && exponent + 1 < bytes.len() {
                exponent += 1;
            }
            if digit(exponent) {
                float = true;
                at = exponent;
                while digit(at) {
                    at += 1;
                }
            }
        }
        let literal = &self.text[start..at];
        let value = if float {
            Value::Float(literal.parse().unwrap_or(f64::NAN))
        } else {
            match literal.parse::<i64>() {
                Ok(whole) => Value::Int(whole),
                Err(_) if whole_digits > pyfmt::INT_MAX_STR_DIGITS => {
                    return Err(JsonError {
                        kind: JsonErrorKind::IntDigits,
                        message: pyfmt::IntError::TooManyDigits(whole_digits).to_string(),
                    });
                }
                Err(_) => Value::BigInt(literal.to_string()),
            }
        };
        Ok((value, at))
    }

    /// The four hexadecimal digits at byte `at` as a number.
    fn hex4(&self, at: usize) -> Option<u32> {
        let digits = self.bytes.get(at..at + 4)?;
        let mut code = 0;
        for digit in digits {
            code = code * 16 + (*digit as char).to_digit(16)?;
        }
        Some(code)
    }

    /// The string whose opening quote stands before byte `start`, and the offset after
    /// its closing quote.
    fn string(&mut self, start: usize) -> Result<(String, usize), JsonError> {
        let bytes = self.bytes;
        let unterminated =
            |reader: &Self| reader.error("Unterminated string starting at", start - 1);
        let mut out = String::new();
        let mut chunk = start;
        loop {
            let mut at = chunk;
            let stop = loop {
                match bytes.get(at) {
                    None => return Err(unterminated(self)),
                    Some(b'"') => break b'"',
                    Some(b'\\') => break b'\\',
                    Some(control) if *control <= 0x1f => {
                        return Err(self.error("Invalid control character at", at));
                    }
                    Some(_) => at += 1,
                }
            };
            if stop == b'"' {
                if chunk == start {
                    return Ok((self.text[start..at].to_string(), at + 1));
                }
                out.push_str(&self.text[chunk..at]);
                return Ok((out, at + 1));
            }
            out.push_str(&self.text[chunk..at]);
            let backslash = at;
            let Some(&escape) = bytes.get(at + 1) else {
                return Err(unterminated(self));
            };
            at += 2;
            let c = match escape {
                b'"' => '"',
                b'\\' => '\\',
                b'/' => '/',
                b'b' => '\x08',
                b'f' => '\x0c',
                b'n' => '\n',
                b'r' => '\r',
                b't' => '\t',
                b'u' => {
                    let invalid = |reader: &Self, u: usize| reader.error(INVALID_U_ESCAPE, u);
                    let Some(mut code) = self.hex4(at) else {
                        return Err(invalid(self, at - 1));
                    };
                    at += 4;
                    // The escape of a pair's second half right after its first half makes
                    // one code point with it; any other escape there is read by itself.
                    if (0xd800..0xdc00).contains(&code) && bytes[at..].starts_with(b"\\u") {
                        let Some(second) = self.hex4(at + 2) else {
                            return Err(invalid(self, at + 1));
                        };
                        if (0xdc00..0xe000).contains(&second) {
                            code = 0x10000 + ((code - 0xd800) << 10) + (second - 0xdc00);
                            at += 6;
                        }
                    }
                    match char::from_u32(code) {
                        Some(c) => c,
                        None => {
                            self.lone_surrogate.get_or_insert(backslash);
                            self.surrogate_mark
                        }
                    }
                }
                _ => return Err(self.error("Invalid \\escape", backslash)),
            };
            out.push(c);
            chunk = at;
        }
    }
}

fn too_deep(what: &str) -> JsonError {
    JsonError {
        kind: JsonErrorKind::TooDeep,
        message: format!(
            "Stack overflow (used 10200 kB) while decoding a JSON {what} from a unicode string"
        ),
    }
}

/// The whole of `text` as one value, and the byte offset of its first escape of half a
/// surrogate pair, each such half read as `surrogate_mark`.
fn read_whole(text: &str, surrogate_mark: char) -> Result<(Value, Option<usize>), JsonError> {
    let mut reader = Reader {
        text,
        bytes: text.as_bytes(),
        lone_surrogate: None,
        surrogate_mark,
    };
    if text.starts_with('\u{feff}') {
        return Err(reader.error("Unexpected UTF-8 BOM (decode using utf-8-sig)", 0));
    }
    let start = reader.skip_space(0);
    let (value, end) = reader.value(start)?;
    let end = reader.skip_space(end);
    if end != text.len() {
        return Err(reader.error("Extra data", end));
    }
    Ok((value, reader.lone_surrogate))
}

/// A private-use character that `text` neither holds nor spells as an escaped pair.
fn unused_mark(text: &str) -> Option<char> {
    const PRIVATE: std::ops::RangeInclusive<u32> = 0xf0000..=0x10fffd;
    let mut used: std::collections::HashSet<u32> = text
        .chars()
        .map(u32::from)
        .filter(|code| PRIVATE.contains(code))
        .collect();
    let bytes = text.as_bytes();
    let hex4 = |at: usize| {
        let digits = std::str::from_utf8(bytes.get(at..at + 4)?).ok()?;
        u32::from_str_radix(digits, 16).ok()
    };
    for (at, _) in text.match_indices("\\u") {
        if let (Some(high @ 0xd800..0xdc00), Some(b"\\u"), Some(low @ 0xdc00..0xe000)) =
            (hex4(at + 2), bytes.get(at + 6..at + 8), hex4(at + 8))
        {
            used.insert(0x10000 + ((high - 0xd800) << 10) + (low - 0xdc00));
        }
    }
    PRIVATE
        .filter(|code| !used.contains(code))
        .find_map(char::from_u32)
}

/// Whether a string of `value`, a member name included, holds `mark`.
fn holds_char(value: &Value, mark: char) -> bool {
    let mut found = false;
    walk_depths(value, |value, _| match value {
        Value::Str(text) => found |= text.contains(mark),
        Value::Dict(members) => found |= members.keys().any(|name| name.contains(mark)),
        _ => {}
    });
    found
}

/// `json.loads(text)`.
///
/// A text the interpreter reads into a value that holds half a surrogate pair in one of
/// its strings is `JsonErrorKind::LoneSurrogate`. A text it refuses for another reason is
/// that refusal, whatever escapes stand before it, and a text whose every such string is
/// replaced by a later member of the same name is read as the interpreter reads it.
pub fn loads(text: &str) -> Result<Value, JsonError> {
    let (value, lone_surrogate) = read_whole(text, '\u{fffd}')?;
    let Some(at) = lone_surrogate else {
        return Ok(value);
    };
    // The halves are read again as a character the text has nowhere else: a string of
    // the value holds it exactly when the interpreter's string holds half a pair.
    if let Some(mark) = unused_mark(text) {
        let (value, _) = read_whole(text, mark)?;
        if !holds_char(&value, mark) {
            return Ok(value);
        }
    }
    let before = text.as_bytes()[..at]
        .iter()
        .filter(|byte| (**byte & 0xc0) != 0x80)
        .count();
    Err(JsonError {
        kind: JsonErrorKind::LoneSurrogate,
        message: format!("a string holds half a surrogate pair (char {before})"),
    })
}

/// Appends `text` as `json.dumps` writes a `str`.
fn write_string(out: &mut String, text: &str, ensure_ascii: bool) {
    out.push('"');
    let plain =
        |byte: u8| byte >= 0x20 && byte != b'"' && byte != b'\\' && (!ensure_ascii || byte < 0x7f);
    if text.bytes().all(plain) {
        out.push_str(text);
        out.push('"');
        return;
    }
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\x08' => out.push_str("\\b"),
            '\x0c' => out.push_str("\\f"),
            c if (c as u32) < 0x20 || (ensure_ascii && (c as u32) > 0x7e) => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units) {
                    let _ = write!(out, "\\u{unit:04x}");
                }
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

/// How a writer lays a value out.
struct Layout {
    sort_keys: bool,
    ensure_ascii: bool,
    indent: Option<usize>,
    item_separator: &'static str,
    key_separator: &'static str,
    /// `json_compact`: a float through `_json_float`.
    compact_floats: bool,
}

/// A list or an object the writer is inside, with its members still to write.
enum Walk<'a> {
    List(std::slice::Iter<'a, Value>),
    Dict(std::vec::IntoIter<(&'a String, &'a Value)>),
}

fn write_json(out: &mut String, value: &Value, layout: &Layout) {
    let new_line = |out: &mut String, level: usize| {
        if let Some(indent) = layout.indent {
            out.push('\n');
            out.extend(std::iter::repeat_n(' ', indent * level));
        }
    };
    let mut walks: Vec<Walk> = Vec::new();
    let mut current = value;
    loop {
        let mut opened = false;
        match current {
            Value::Null => out.push_str("null"),
            Value::Bool(true) => out.push_str("true"),
            Value::Bool(false) => out.push_str("false"),
            Value::Int(whole) => {
                let _ = write!(out, "{whole}");
            }
            Value::BigInt(whole) => out.push_str(whole),
            Value::Float(float) if layout.compact_floats => write_json_float(out, *float),
            Value::Float(float) => {
                if float.is_nan() {
                    out.push_str("NaN");
                } else if float.is_infinite() {
                    out.push_str(if *float > 0.0 {
                        "Infinity"
                    } else {
                        "-Infinity"
                    });
                } else {
                    pyfmt::write_repr_f64(out, *float);
                }
            }
            Value::Str(text) => write_string(out, text, layout.ensure_ascii),
            Value::List(items) if items.is_empty() => out.push_str("[]"),
            Value::Dict(members) if members.is_empty() => out.push_str("{}"),
            Value::List(items) => {
                out.push('[');
                walks.push(Walk::List(items.iter()));
                opened = true;
            }
            Value::Dict(members) => {
                out.push('{');
                let mut pairs: Vec<(&String, &Value)> = members.iter().collect();
                if layout.sort_keys {
                    pairs.sort_unstable_by(|a, b| a.0.cmp(b.0));
                }
                walks.push(Walk::Dict(pairs.into_iter()));
                opened = true;
            }
        }
        // The next value to write: the next member of the innermost list or object that
        // has one, every one that has none closed on the way.
        loop {
            let level = walks.len();
            let next = match walks.last_mut() {
                None => return,
                Some(Walk::List(items)) => items.next().map(|item| (None, item)),
                Some(Walk::Dict(pairs)) => pairs.next().map(|(name, item)| (Some(name), item)),
            };
            match next {
                Some((name, item)) => {
                    if !opened {
                        out.push_str(layout.item_separator);
                    }
                    new_line(out, level);
                    if let Some(name) = name {
                        write_string(out, name, layout.ensure_ascii);
                        out.push_str(layout.key_separator);
                    }
                    current = item;
                    break;
                }
                None => {
                    new_line(out, level - 1);
                    out.push(match walks.pop() {
                        Some(Walk::Dict(_)) => '}',
                        _ => ']',
                    });
                    opened = false;
                }
            }
        }
    }
}

/// `json.dumps(value, indent=..., sort_keys=..., ensure_ascii=...)`.
pub fn dumps(value: &Value, options: &DumpOptions) -> String {
    let mut out = String::new();
    write_json(
        &mut out,
        value,
        &Layout {
            sort_keys: options.sort_keys,
            ensure_ascii: options.ensure_ascii,
            indent: options.indent,
            item_separator: if options.indent.is_some() { "," } else { ", " },
            key_separator: ": ",
            compact_floats: false,
        },
    );
    out
}

/// `json.dumps(text, ensure_ascii=...)` of one `str`.
pub fn dumps_str(text: &str, ensure_ascii: bool) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    write_string(&mut out, text, ensure_ascii);
    out
}

/// Compact text, the members of an object sorted by name.
///
/// Any value is written, however deep it nests; `compact_frames` tells a caller whether
/// the value nests past the call depth at which the comparison refuses it.
pub fn dumps_compact(value: &Value) -> String {
    let mut out = String::new();
    write_compact(&mut out, value);
    out
}

/// Appends `dumps_compact(value)` to `out`.
pub fn write_compact(out: &mut String, value: &Value) {
    write_json(
        out,
        value,
        &Layout {
            sort_keys: true,
            ensure_ascii: false,
            indent: None,
            item_separator: ",",
            key_separator: ":",
            compact_floats: true,
        },
    );
}

fn write_json_float(out: &mut String, value: f64) {
    if !value.is_finite() {
        out.push_str("null");
        return;
    }
    let from = out.len();
    pyfmt::write_repr_f64(out, value);
    // `1e-05` and `1e+16` lose the padding and the plus sign of their exponent.
    if let Some(e) = out[from..].find('e') {
        let exponent = from + e + 1;
        let negative = out[exponent..].starts_with('-');
        let digits = out[exponent + 1..].trim_start_matches('0').to_string();
        out.truncate(exponent);
        if negative {
            out.push('-');
        }
        out.push_str(if digits.is_empty() { "0" } else { &digits });
    }
}

/// A float as compact JSON writes it: `null` for a value that is not finite, else the
/// shortest text that reads back, its exponent without padding or plus sign.
pub fn json_float(value: f64) -> String {
    let mut out = String::with_capacity(24);
    write_json_float(&mut out, value);
    out
}

/// A value as one field's text: a string as it is, a list as its items joined by commas,
/// any other value as compact JSON.
pub fn json_scalar(value: &Value) -> String {
    let mut out = String::new();
    // The items of a list, and of the lists inside it, joined with commas in order.
    let mut lists: Vec<std::slice::Iter<Value>> = Vec::new();
    let mut current = value;
    loop {
        let mut opened = false;
        match current {
            Value::Str(text) => out.push_str(text),
            Value::List(items) => {
                lists.push(items.iter());
                opened = true;
            }
            other => out.push_str(&dumps_compact(other)),
        }
        loop {
            let Some(items) = lists.last_mut() else {
                return out;
            };
            match items.next() {
                Some(item) => {
                    if !opened {
                        out.push(',');
                    }
                    current = item;
                    break;
                }
                None => {
                    lists.pop();
                    opened = false;
                }
            }
        }
    }
}

/// `repr(value)` of the object `json.loads` returns.
pub fn py_repr(value: &Value) -> String {
    let mut out = String::new();
    let mut walks: Vec<Walk> = Vec::new();
    let mut current = value;
    loop {
        let mut opened = false;
        match current {
            Value::Null => out.push_str("None"),
            Value::Bool(true) => out.push_str("True"),
            Value::Bool(false) => out.push_str("False"),
            Value::Int(whole) => {
                let _ = write!(out, "{whole}");
            }
            Value::BigInt(whole) => out.push_str(whole),
            Value::Float(float) => pyfmt::write_repr_f64(&mut out, *float),
            Value::Str(text) => out.push_str(&pyfmt::repr_str(text)),
            Value::List(items) => {
                out.push('[');
                walks.push(Walk::List(items.iter()));
                opened = true;
            }
            Value::Dict(members) => {
                out.push('{');
                walks.push(Walk::Dict(members.iter().collect::<Vec<_>>().into_iter()));
                opened = true;
            }
        }
        loop {
            let next = match walks.last_mut() {
                None => return out,
                Some(Walk::List(items)) => items.next().map(|item| (None, item)),
                Some(Walk::Dict(pairs)) => pairs.next().map(|(name, item)| (Some(name), item)),
            };
            match next {
                Some((name, item)) => {
                    if !opened {
                        out.push_str(", ");
                    }
                    if let Some(name) = name {
                        out.push_str(&pyfmt::repr_str(name));
                        out.push_str(": ");
                    }
                    current = item;
                    break;
                }
                None => {
                    out.push(match walks.pop() {
                        Some(Walk::Dict(_)) => '}',
                        _ => ']',
                    });
                    opened = false;
                }
            }
        }
    }
}

/// `str(value)` of the object `json.loads` returns.
pub fn py_str(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        other => py_repr(other),
    }
}

/// A number as the interpreter compares it: a `bool` is 1 or 0.
enum Number<'a> {
    Int(i64),
    Big(&'a str),
    Float(f64),
}

fn number_of(value: &Value) -> Option<Number<'_>> {
    match value {
        Value::Bool(flag) => Some(Number::Int(i64::from(*flag))),
        Value::Int(whole) => Some(Number::Int(*whole)),
        Value::BigInt(whole) => Some(Number::Big(whole)),
        Value::Float(float) => Some(Number::Float(*float)),
        _ => None,
    }
}

/// `a == b` of two objects `json.loads` returns.
///
/// Numbers compare by value whatever their type (`1 == 1.0 == True`), objects whatever
/// the order of their members. Two NaN compare unequal by themselves and equal as items
/// of two lists or objects: `json.loads` returns one NaN object for every `NaN` it
/// reads, and a list compares an item with itself as equal.
pub fn py_eq(a: &Value, b: &Value) -> bool {
    eq_inside(a, b, false)
}

fn eq_inside(a: &Value, b: &Value, item: bool) -> bool {
    match (a, b) {
        (Value::Null, Value::Null) => true,
        (Value::Str(a), Value::Str(b)) => a == b,
        (Value::List(a), Value::List(b)) => {
            a.len() == b.len() && a.iter().zip(b).all(|(a, b)| eq_inside(a, b, true))
        }
        (Value::Dict(a), Value::Dict(b)) => {
            a.len() == b.len()
                && a.iter()
                    .all(|(name, a)| b.get(name).is_some_and(|b| eq_inside(a, b, true)))
        }
        _ => match (number_of(a), number_of(b)) {
            (Some(Number::Int(a)), Some(Number::Int(b))) => a == b,
            (Some(Number::Big(a)), Some(Number::Big(b))) => a == b,
            (Some(Number::Float(a)), Some(Number::Float(b))) => {
                a == b || (item && a.is_nan() && b.is_nan())
            }
            (Some(Number::Float(float)), Some(Number::Int(whole)))
            | (Some(Number::Int(whole)), Some(Number::Float(float))) => pyfmt::int_of_f64(float)
                .is_ok_and(|text| float.fract() == 0.0 && text == whole.to_string()),
            (Some(Number::Float(float)), Some(Number::Big(whole)))
            | (Some(Number::Big(whole)), Some(Number::Float(float))) => {
                pyfmt::int_of_f64(float).is_ok_and(|text| float.fract() == 0.0 && text == whole)
            }
            _ => false,
        },
    }
}

/// Why `float(value)` raises; each holds the exception's text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FloatError {
    /// `TypeError`: the value is `None`, a list or an object.
    Type(String),
    /// `ValueError`: the value is a `str` that spells no number.
    Value(String),
    /// `OverflowError`: the value is a whole number beyond the largest float.
    Overflow(String),
}

impl fmt::Display for FloatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FloatError::Type(text) | FloatError::Value(text) | FloatError::Overflow(text) => {
                f.write_str(text)
            }
        }
    }
}

/// `float(value)` of the object `json.loads` returns.
pub fn py_float(value: &Value) -> Result<f64, FloatError> {
    let no_number = |class: &str| {
        FloatError::Type(format!(
            "float() argument must be a string or a real number, not '{class}'"
        ))
    };
    match value {
        Value::Null => Err(no_number("NoneType")),
        Value::List(_) => Err(no_number("list")),
        Value::Dict(_) => Err(no_number("dict")),
        Value::Bool(flag) => Ok(f64::from(u8::from(*flag))),
        Value::Int(whole) => Ok(*whole as f64),
        Value::BigInt(whole) => pyfmt::f64_of_int(whole)
            .ok_or_else(|| FloatError::Overflow("int too large to convert to float".to_string())),
        Value::Float(float) => Ok(*float),
        Value::Str(text) => {
            pyfmt::float_of(text).ok_or_else(|| FloatError::Value(pyfmt::float_error(text)))
        }
    }
}

/// `bool(value)` of the object `json.loads` returns.
pub fn truthy(value: &Value) -> bool {
    match value {
        Value::Null => false,
        Value::Bool(flag) => *flag,
        Value::Int(whole) => *whole != 0,
        Value::BigInt(_) => true,
        Value::Float(float) => *float != 0.0,
        Value::Str(text) => !text.is_empty(),
        Value::List(items) => !items.is_empty(),
        Value::Dict(members) => !members.is_empty(),
    }
}

/// Calls `visit` with every value inside `value`, itself included, and the count of
/// lists and objects around it.
fn walk_depths(value: &Value, mut visit: impl FnMut(&Value, usize)) {
    let mut pending = vec![(value, 0usize)];
    while let Some((value, around)) = pending.pop() {
        visit(value, around);
        match value {
            Value::List(items) => pending.extend(items.iter().map(|item| (item, around + 1))),
            Value::Dict(members) => pending.extend(members.values().map(|item| (item, around + 1))),
            _ => {}
        }
    }
}

/// The count of lists and objects around the deepest value: 0 for a scalar, 1 for an
/// empty list.
pub fn depth(value: &Value) -> usize {
    let mut deepest = 0;
    walk_depths(value, |value, around| {
        let own = usize::from(matches!(value, Value::List(_) | Value::Dict(_)));
        deepest = deepest.max(around + own);
    });
    deepest
}

/// The call depth `json_compact(value)` reaches, counted as CPython counts the frames of
/// the recursive writing, its own first frame included.
///
/// Called by a function that is itself `n` frames deep, `json_compact` raises
/// `RecursionError` exactly when `n + compact_frames(value)` is more than
/// `RECURSION_LIMIT`. A list or an object costs two frames a level (the function and the
/// generator that joins its items); the innermost call costs what it calls: four for
/// `None`, a `bool` or a member name (`json.dumps` and its encoder), three for a `str`,
/// two for a float, one for a whole number.
pub fn compact_frames(value: &Value) -> usize {
    let mut most = 0;
    walk_depths(value, |value, around| {
        let own = match value {
            Value::Null | Value::Bool(_) => 4,
            Value::Str(_) => 3,
            Value::Float(_) => 2,
            Value::Int(_) | Value::BigInt(_) => 1,
            Value::List(_) => 2,
            Value::Dict(members) if members.is_empty() => 2,
            Value::Dict(_) => 4,
        };
        most = most.max(2 * around + own);
    });
    most
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(text: &str) -> Value {
        loads(text).unwrap_or_else(|error| panic!("{text}: {error}"))
    }

    fn refusal(text: &str) -> String {
        match loads(text) {
            Err(error) if error.kind == JsonErrorKind::Decode => error.message,
            other => panic!("{text}: {other:?}"),
        }
    }

    #[test]
    fn json_scalar_renders_like_the_rust_comparator() {
        assert_eq!(json_scalar(&Value::Str("x".into())), "x");
        assert_eq!(json_scalar(&read(r#"["a", "b"]"#)), "a,b");
        assert_eq!(json_scalar(&Value::Int(1)), "1");
        assert_eq!(json_scalar(&Value::Bool(true)), "true");
        assert_eq!(json_scalar(&Value::Null), "null");
        assert_eq!(json_scalar(&Value::Float(0.00001)), "1e-5");
        assert_eq!(json_scalar(&Value::Float(1e16)), "1e16");
        assert_eq!(json_scalar(&Value::Float(0.5)), "0.5");
        assert_eq!(
            json_scalar(&read(r#"{"name": "PF1", "db": "Pfam"}"#)),
            r#"{"db":"Pfam","name":"PF1"}"#
        );
        assert_eq!(
            json_scalar(&read(r#"[{"id": "rs1", "af": 0.25}]"#)),
            r#"{"af":0.25,"id":"rs1"}"#
        );
    }

    #[test]
    fn a_text_is_read_into_the_values_the_interpreter_makes_of_it() {
        assert_eq!(read(" null "), Value::Null);
        assert_eq!(read("-0"), Value::Int(0));
        assert_eq!(read("-0.0"), Value::Float(-0.0));
        assert!(matches!(read("-0.0"), Value::Float(zero) if zero.is_sign_negative()));
        assert_eq!(read("1E400"), Value::Float(f64::INFINITY));
        assert_eq!(read("-Infinity"), Value::Float(f64::NEG_INFINITY));
        assert!(matches!(read("NaN"), Value::Float(nan) if nan.is_nan()));
        assert_eq!(read("9223372036854775807"), Value::Int(i64::MAX));
        assert_eq!(
            read("-9223372036854775809"),
            Value::BigInt("-9223372036854775809".into())
        );
        assert_eq!(
            read(r#""a\u00e9\ud83d\ude00\/\b""#),
            Value::Str("aé\u{1f600}/\x08".into())
        );
        let object = read(r#"{"b":1,"a":2,"b":[true,{}]}"#);
        let members = object.as_dict().unwrap();
        assert_eq!(members.keys().collect::<Vec<_>>(), ["b", "a"]);
        assert_eq!(members["b"], read("[true,{}]"));
        assert_eq!(object.get("a"), Some(&Value::Int(2)));
    }

    #[test]
    fn a_refused_text_carries_the_interpreters_message_with_its_place_in_code_points() {
        for (text, message) in [
            ("", "Expecting value: line 1 column 1 (char 0)"),
            (
                "\u{feff}1",
                "Unexpected UTF-8 BOM (decode using utf-8-sig): line 1 column 1 (char 0)",
            ),
            ("1.", "Extra data: line 1 column 2 (char 1)"),
            ("1e+", "Extra data: line 1 column 2 (char 1)"),
            ("-Infinit", "Expecting value: line 1 column 1 (char 0)"),
            (
                "\"abc",
                "Unterminated string starting at: line 1 column 1 (char 0)",
            ),
            (
                "\"\\u1234",
                "Unterminated string starting at: line 1 column 1 (char 0)",
            ),
            ("\"a\\x\"", "Invalid \\escape: line 1 column 3 (char 2)"),
            (
                "\"\\u12\"",
                concat!("Invalid \\uXX", "XX escape: line 1 column 3 (char 2)"),
            ),
            (
                "\"\\ud83d\\u12\"",
                concat!("Invalid \\uXX", "XX escape: line 1 column 9 (char 8)"),
            ),
            (
                "\"a\tb\"",
                "Invalid control character at: line 1 column 3 (char 2)",
            ),
            (
                "[1,]",
                "Illegal trailing comma before end of array: line 1 column 3 (char 2)",
            ),
            (
                "{\"a\":1,}",
                "Illegal trailing comma before end of object: line 1 column 7 (char 6)",
            ),
            (
                "{a:1}",
                "Expecting property name enclosed in double quotes: line 1 column 2 (char 1)",
            ),
            (
                "{\"a\" 1}",
                "Expecting ':' delimiter: line 1 column 6 (char 5)",
            ),
            ("[1 2]", "Expecting ',' delimiter: line 1 column 4 (char 3)"),
            (
                "\n\n  [1,\n 2,\n x]",
                "Expecting value: line 5 column 2 (char 13)",
            ),
            ("[\"é\", é]", "Expecting value: line 1 column 7 (char 6)"),
        ] {
            assert_eq!(refusal(text), message, "{text:?}");
        }
    }

    #[test]
    fn half_a_surrogate_pair_counts_only_in_a_string_the_value_keeps() {
        assert_eq!(loads(r#"{"a":"\ud83d","a":1}"#), Ok(read(r#"{"a":1}"#)));
        assert_eq!(
            loads(r#"{"a":["\ud83d",{"\ude00":2}],"b":3,"a":"\ufffd"}"#),
            Ok(read(r#"{"a":"\ufffd","b":3}"#))
        );
        let private = "{\"a\":\"\\ud83d\",\"a\":\"\u{f0000}\\udb80\\udc01\"}";
        assert_eq!(loads(private), Ok(read("{\"a\":\"\u{f0000}\u{f0001}\"}")));
        for kept in [
            r#"{"a":"\ud83d","b":1}"#,
            r#"{"\ud83d":1,"\ud83e":2}"#,
            r#"{"\ud83d":1,"\ufffd":1}"#,
            r#"{"a":1,"a":"\ud83d"}"#,
            r#"{"a":{"b":"\ud83d","b":1,"c":["\ude00"]}}"#,
        ] {
            assert_eq!(
                loads(kept).unwrap_err().kind,
                JsonErrorKind::LoneSurrogate,
                "{kept}"
            );
        }
        assert_eq!(
            loads(r#"{"a":"\ud83d","a":1,}"#).unwrap_err().kind,
            JsonErrorKind::Decode
        );
    }

    #[test]
    fn half_a_surrogate_pair_too_many_digits_and_too_deep_a_nesting_are_told_apart() {
        assert_eq!(
            loads(r#""\ud83d""#).unwrap_err().kind,
            JsonErrorKind::LoneSurrogate
        );
        assert_eq!(
            loads(r#"["\ude00x"]"#).unwrap_err().kind,
            JsonErrorKind::LoneSurrogate
        );
        assert_eq!(
            loads(r#""\ud83d\u0041"#).unwrap_err().kind,
            JsonErrorKind::Decode
        );
        let digits = loads(&"9".repeat(4301)).unwrap_err();
        assert_eq!(digits.kind, JsonErrorKind::IntDigits);
        assert!(digits.is_value_error() && digits.message.contains("value has 4301 digits"));
        assert!(matches!(read(&"9".repeat(4300)), Value::BigInt(_)));
        for (open, close) in [("[", "]"), ("{\"a\":", "}")] {
            let deepest = format!("{}1{}", open.repeat(MAX_DEPTH), close.repeat(MAX_DEPTH));
            assert_eq!(depth(&read(&deepest)), MAX_DEPTH);
            let deeper = format!(
                "{}1{}",
                open.repeat(MAX_DEPTH + 1),
                close.repeat(MAX_DEPTH + 1)
            );
            let error = loads(&deeper).unwrap_err();
            assert!(error.is_recursion_error() && !error.is_value_error());
        }
        assert_eq!(
            loads(&"[".repeat(60_000)).unwrap_err().kind,
            JsonErrorKind::Decode
        );
        assert_eq!(
            loads(&"[".repeat(100_000)).unwrap_err().kind,
            JsonErrorKind::TooDeep
        );
    }

    #[test]
    fn dumps_writes_the_separators_escapes_and_floats_of_the_interpreter() {
        let value = read(
            r#"{"b": [1, 2.5, {"é": "a\"\n\u007f\ud83d\ude00"}], "a": {}, "c": [], "n": [NaN, -Infinity, 1e16, null, true]}"#,
        );
        assert_eq!(
            dumps(&value, &DumpOptions::default()),
            r#"{"b": [1, 2.5, {"\u00e9": "a\"\n\u007f\ud83d\ude00"}], "a": {}, "c": [], "n": [NaN, -Infinity, 1e+16, null, true]}"#
        );
        let sorted = DumpOptions {
            sort_keys: true,
            ensure_ascii: false,
            ..DumpOptions::default()
        };
        assert_eq!(
            dumps(&value, &sorted),
            "{\"a\": {}, \"b\": [1, 2.5, {\"é\": \"a\\\"\\n\u{7f}\u{1f600}\"}], \"c\": [], \"n\": [NaN, -Infinity, 1e+16, null, true]}"
        );
        let indented = DumpOptions {
            indent: Some(1),
            ..DumpOptions::default()
        };
        assert_eq!(
            dumps(&read(r#"{"a": [1, {"b": []}], "c": {}}"#), &indented),
            "{\n \"a\": [\n  1,\n  {\n   \"b\": []\n  }\n ],\n \"c\": {}\n}"
        );
        assert_eq!(dumps_str("a\u{1}é", true), "\"a\\u0001\\u00e9\"");
        assert_eq!(dumps_str("a\u{1}é", false), "\"a\\u0001é\"");
    }

    #[test]
    fn the_compact_form_sorts_members_and_writes_floats_without_exponent_padding() {
        let value =
            read(r#"{"b": [1, 2.50, 1e-05, 1E+22, NaN], "a": {"z": null, "é": "é\"", "A": -0.0}}"#);
        assert_eq!(
            dumps_compact(&value),
            r#"{"a":{"A":-0.0,"z":null,"é":"é\""},"b":[1,2.5,1e-5,1e22,null]}"#
        );
        let mut line = String::from("x=");
        write_compact(&mut line, &value);
        assert_eq!(line[2..], dumps_compact(&value));
        assert_eq!(json_float(1.5e300), "1.5e300");
        assert_eq!(json_float(5e-324), "5e-324");
        assert_eq!(json_float(f64::INFINITY), "null");
        assert_eq!(
            json_scalar(&read(r#"["a", ["b", 1], [], {"k": [2]}]"#)),
            r#"a,b,1,,{"k":[2]}"#
        );
    }

    #[test]
    fn a_value_prints_and_compares_as_the_interpreters_object() {
        let value = read(r#"{"a": [1, 2.0, "x'y", null, true], "b": {}}"#);
        assert_eq!(
            py_repr(&value),
            r#"{'a': [1, 2.0, "x'y", None, True], 'b': {}}"#
        );
        assert_eq!(py_str(&value), py_repr(&value));
        assert_eq!(py_str(&Value::Str("x".into())), "x");
        assert_eq!(py_repr(&Value::Str("x".into())), "'x'");
        assert!(py_eq(&read("1"), &read("1.0")) && py_eq(&read("true"), &read("1")));
        assert!(py_eq(
            &read(r#"{"a":1,"b":[2]}"#),
            &read(r#"{"b":[2.0],"a":true}"#)
        ));
        assert!(!py_eq(&read("NaN"), &read("NaN")) && py_eq(&read("[NaN]"), &read("[NaN]")));
        assert!(!py_eq(&read("1"), &read("\"1\"")) && !py_eq(&read("[]"), &read("{}")));
        assert!(py_eq(
            &read("1e30"),
            &read("1000000000000000019884624838656")
        ));
        assert_eq!(Value::of_int("-12"), Value::Int(-12));
        assert_eq!(
            Value::of_int("9223372036854775808"),
            read("9223372036854775808")
        );
        assert_eq!(py_float(&read("true")), Ok(1.0));
        assert_eq!(py_float(&read("\" 1_0 \"")), Ok(10.0));
        assert_eq!(
            py_float(&read("9223372036854775808")),
            Ok(9223372036854775808.0)
        );
        assert_eq!(
            py_float(&Value::Null).unwrap_err().to_string(),
            "float() argument must be a string or a real number, not 'NoneType'"
        );
        assert!(matches!(
            py_float(&read("\"x\"")),
            Err(FloatError::Value(_))
        ));
        assert!(matches!(
            py_float(&Value::of_int(&"9".repeat(400))),
            Err(FloatError::Overflow(_))
        ));
        assert!(
            !truthy(&read("0.0"))
                && !truthy(&read("[]"))
                && truthy(&read("NaN"))
                && truthy(&read("[0]"))
        );
    }

    #[test]
    fn the_frames_of_the_compact_form_count_two_a_level_and_the_innermost_call() {
        assert_eq!(compact_frames(&read("1")), 1);
        assert_eq!(compact_frames(&read("[]")), 2);
        assert_eq!(compact_frames(&read("[[1]]")), 5);
        assert_eq!(compact_frames(&read(r#"{"a":1}"#)), 4);
        assert_eq!(compact_frames(&read(r#"{"a":{"a":1}}"#)), 6);
        assert_eq!(compact_frames(&read(r#"[1.5, "x", null]"#)), 6);
        assert_eq!(depth(&read(r#"{"a":[1,{"b":[]}]}"#)), 4);
        assert_eq!(depth(&read("1")), 0);
    }

    #[test]
    fn a_value_nested_as_deep_as_the_reader_goes_is_written_and_dropped() {
        let text = format!("{}{}", "[".repeat(MAX_DEPTH), "]".repeat(MAX_DEPTH));
        let value = read(&text);
        assert_eq!(dumps_compact(&value), text);
        assert_eq!(json_scalar(&value), "");
        assert_eq!(py_repr(&value), text);
        assert_eq!(compact_frames(&value), 2 * MAX_DEPTH);
        assert_eq!(value.into_list().map(|items| items.len()), Ok(1));
    }
}
