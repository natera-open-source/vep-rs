// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The stored translation of named transcripts in a vep-rs JSON cache: the shards
//! `transcripts/<contig>/*.json` and `*.json.gz` in the order `sorted()` gives their
//! `pathlib` paths, each read as `json.load` reads it.

use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::io::Read;
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

/// The deepest nesting a shard may have.
const MAX_DEPTH: usize = 1000;
/// The most digits an integer may have (`sys.get_int_max_str_digits()`).
const MAX_INT_DIGITS: usize = 4300;

/// The value a transcript's cache holds under `peptide`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Peptide {
    /// `None`: the key or the cache is absent, or the value is `null`.
    Absent,
    Text(String),
    /// A value that is no string, or a string that holds half a surrogate pair: equal
    /// to no peptide of a dump.
    Other,
}

/// What one transcript of a shard says, in the shard's order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event {
    Found {
        id: String,
        peptide: Peptide,
    },
    /// The transcript's cache is a value `.get` cannot be called on: an error where
    /// the id has no peptide yet.
    Unreadable {
        id: String,
    },
}

/// A parsed JSON value: strings and keys as the place of their text between the quotes.
#[derive(Debug)]
enum Value {
    Null,
    Bool(bool),
    /// A number, with whether it is zero.
    Number {
        zero: bool,
    },
    Str(usize, usize),
    List(Vec<Value>),
    Dict(Vec<((usize, usize), Value)>),
}

impl Value {
    /// Python's truth of the value.
    fn truthy(&self) -> bool {
        match self {
            Value::Null => false,
            Value::Bool(value) => *value,
            Value::Number { zero } => !zero,
            Value::Str(from, to) => from != to,
            Value::List(items) => !items.is_empty(),
            Value::Dict(items) => !items.is_empty(),
        }
    }
}

struct Parser<'a> {
    text: &'a [u8],
    at: usize,
}

type Parsed<T> = Result<T, String>;

impl<'a> Parser<'a> {
    fn fail<T>(&self, what: &str) -> Parsed<T> {
        Err(format!("{what} at byte {}", self.at))
    }

    fn skip_whitespace(&mut self) {
        while matches!(self.text.get(self.at), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.at += 1;
        }
    }

    fn literal(&mut self, word: &str) -> bool {
        if self.text[self.at..].starts_with(word.as_bytes()) {
            self.at += word.len();
            true
        } else {
            false
        }
    }

    fn value(&mut self, depth: usize) -> Parsed<Value> {
        if depth > MAX_DEPTH {
            return self.fail("nesting deeper than 1,000 levels");
        }
        match self.text.get(self.at) {
            Some(b'"') => {
                let (from, to) = self.string()?;
                Ok(Value::Str(from, to))
            }
            Some(b'{') => self.dict(depth),
            Some(b'[') => self.list(depth),
            Some(b'n') if self.literal("null") => Ok(Value::Null),
            Some(b't') if self.literal("true") => Ok(Value::Bool(true)),
            Some(b'f') if self.literal("false") => Ok(Value::Bool(false)),
            Some(b'N') if self.literal("NaN") => Ok(Value::Number { zero: false }),
            Some(b'I') if self.literal("Infinity") => Ok(Value::Number { zero: false }),
            Some(b'-') if self.literal("-Infinity") => Ok(Value::Number { zero: false }),
            Some(b'-' | b'0'..=b'9') => self.number(),
            _ => self.fail("expecting a value"),
        }
    }

    /// A string's text between its quotes, checked: no control character, every
    /// escape one JSON defines.
    fn string(&mut self) -> Parsed<(usize, usize)> {
        self.at += 1;
        let from = self.at;
        loop {
            match self.text.get(self.at) {
                None => return self.fail("an unterminated string"),
                Some(b'"') => {
                    self.at += 1;
                    return Ok((from, self.at - 1));
                }
                Some(b'\\') => match self.text.get(self.at + 1) {
                    Some(b'"' | b'\\' | b'/' | b'b' | b'f' | b'n' | b'r' | b't') => {
                        self.at += 2;
                    }
                    Some(b'u') => {
                        let digits = self.text.get(self.at + 2..self.at + 6);
                        if !digits.is_some_and(|digits| digits.iter().all(u8::is_ascii_hexdigit)) {
                            return self.fail("a \\u escape without four hexadecimal digits");
                        }
                        self.at += 6;
                    }
                    _ => return self.fail("an invalid \\escape"),
                },
                Some(byte) if *byte < 0x20 => return self.fail("an invalid control character"),
                Some(_) => self.at += 1,
            }
        }
    }

    fn digits(&mut self) -> usize {
        let from = self.at;
        while matches!(self.text.get(self.at), Some(b'0'..=b'9')) {
            self.at += 1;
        }
        self.at - from
    }

    /// `-?(0|[1-9]\d*)(\.\d+)?([eE][-+]?\d+)?`, as much of it as stands there.
    fn number(&mut self) -> Parsed<Value> {
        let from = self.at;
        if self.text.get(self.at) == Some(&b'-') {
            self.at += 1;
        }
        match self.text.get(self.at) {
            Some(b'0') => self.at += 1,
            Some(b'1'..=b'9') => {
                self.digits();
            }
            _ => return self.fail("expecting a value"),
        }
        let integer_end = self.at;
        let mut float = false;
        if self.text.get(self.at) == Some(&b'.')
            && matches!(self.text.get(self.at + 1), Some(b'0'..=b'9'))
        {
            self.at += 1;
            self.digits();
            float = true;
        }
        if matches!(self.text.get(self.at), Some(b'e' | b'E')) {
            let mark = self.at;
            self.at += 1;
            if matches!(self.text.get(self.at), Some(b'-' | b'+')) {
                self.at += 1;
            }
            if self.digits() == 0 {
                self.at = mark;
            } else {
                float = true;
            }
        }
        let text = &self.text[from..self.at];
        if float {
            let number: f64 = std::str::from_utf8(text)
                .ok()
                .and_then(|text| text.parse().ok())
                .unwrap_or(0.0);
            return Ok(Value::Number {
                zero: number == 0.0,
            });
        }
        let digits = &text[usize::from(text[0] == b'-')..integer_end - from];
        if digits.len() > MAX_INT_DIGITS {
            return self.fail("an integer of more than 4,300 digits");
        }
        Ok(Value::Number {
            zero: digits == b"0",
        })
    }

    fn list(&mut self, depth: usize) -> Parsed<Value> {
        self.at += 1;
        let mut items = Vec::new();
        self.skip_whitespace();
        if self.text.get(self.at) == Some(&b']') {
            self.at += 1;
            return Ok(Value::List(items));
        }
        loop {
            items.push(self.value(depth + 1)?);
            self.skip_whitespace();
            match self.text.get(self.at) {
                Some(b']') => {
                    self.at += 1;
                    return Ok(Value::List(items));
                }
                Some(b',') => {
                    self.at += 1;
                    self.skip_whitespace();
                }
                _ => return self.fail("expecting ',' or ']'"),
            }
        }
    }

    fn dict(&mut self, depth: usize) -> Parsed<Value> {
        self.at += 1;
        let mut items = Vec::new();
        self.skip_whitespace();
        if self.text.get(self.at) == Some(&b'}') {
            self.at += 1;
            return Ok(Value::Dict(items));
        }
        loop {
            if self.text.get(self.at) != Some(&b'"') {
                return self.fail("expecting a key in double quotes");
            }
            let key = self.string()?;
            self.skip_whitespace();
            if self.text.get(self.at) != Some(&b':') {
                return self.fail("expecting ':'");
            }
            self.at += 1;
            self.skip_whitespace();
            items.push((key, self.value(depth + 1)?));
            self.skip_whitespace();
            match self.text.get(self.at) {
                Some(b'}') => {
                    self.at += 1;
                    return Ok(Value::Dict(items));
                }
                Some(b',') => {
                    self.at += 1;
                    self.skip_whitespace();
                }
                _ => return self.fail("expecting ',' or '}'"),
            }
        }
    }
}

/// The text of a string between its quotes, its escapes read; `None` where it holds
/// half a surrogate pair.
fn decode(raw: &[u8]) -> Option<String> {
    if !raw.contains(&b'\\') {
        return std::str::from_utf8(raw).ok().map(str::to_string);
    }
    let hex = |at: usize| -> u32 {
        std::str::from_utf8(&raw[at..at + 4])
            .ok()
            .and_then(|digits| u32::from_str_radix(digits, 16).ok())
            .unwrap_or(0)
    };
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
            b'f' => Some(12),
            b'n' => Some(b'\n'),
            b'r' => Some(b'\r'),
            b't' => Some(b'\t'),
            b'u' => None,
            other => Some(other),
        };
        if let Some(byte) = simple {
            out.push(byte);
            at += 2;
            continue;
        }
        let mut code = hex(at + 2);
        at += 6;
        if (0xd800..0xdc00).contains(&code) && raw[at..].starts_with(b"\\u") && raw.len() >= at + 6
        {
            let low = hex(at + 2);
            if (0xdc00..0xe000).contains(&low) {
                code = 0x10000 + ((code - 0xd800) << 10) + (low - 0xdc00);
                at += 6;
            }
        }
        let character = char::from_u32(code)?;
        out.extend_from_slice(character.encode_utf8(&mut [0u8; 4]).as_bytes());
    }
    String::from_utf8(out).ok()
}

struct Shard<'a> {
    text: &'a [u8],
}

impl Shard<'_> {
    fn string(&self, value: &Value) -> Option<Option<String>> {
        match value {
            Value::Str(from, to) => Some(decode(&self.text[*from..*to])),
            _ => None,
        }
    }

    /// `dict.get(key)`: the value of the last entry under the key.
    fn get<'v>(&self, items: &'v [((usize, usize), Value)], key: &str) -> Option<&'v Value> {
        items
            .iter()
            .rev()
            .find(|((from, to), _)| {
                let raw = &self.text[*from..*to];
                if raw.contains(&b'\\') {
                    decode(raw).as_deref() == Some(key)
                } else {
                    raw == key.as_bytes()
                }
            })
            .map(|(_, value)| value)
    }
}

/// What the transcripts of one shard say about the wanted ids, in the shard's order.
/// The error is why `json.load` or the walk over its result ends on an exception.
pub fn scan(text: &[u8], wanted: &HashSet<String>) -> Result<Vec<Event>, String> {
    let checked = std::str::from_utf8(text).map_err(|_| "the shard is not UTF-8".to_string())?;
    if checked.starts_with('\u{feff}') {
        return Err("the shard opens with a byte order mark".to_string());
    }
    let mut parser = Parser { text, at: 0 };
    parser.skip_whitespace();
    let data = parser.value(0)?;
    parser.skip_whitespace();
    if parser.at != text.len() {
        return parser.fail("extra data");
    }
    let shard = Shard { text };
    // The lists of a dict in the order of their keys' first writing, each key with its
    // last value.
    let lists: Vec<&Vec<Value>> = match &data {
        Value::List(items) => vec![items],
        Value::Dict(items) => {
            let mut last: HashMap<Option<String>, usize> = HashMap::new();
            let mut order = Vec::new();
            for (index, ((from, to), _)) in items.iter().enumerate() {
                let key = decode(&text[*from..*to]);
                if last.insert(key.clone(), index).is_none() {
                    order.push(key);
                }
            }
            order
                .iter()
                .filter_map(|key| match &items[last[key]].1 {
                    Value::List(list) => Some(list),
                    _ => None,
                })
                .collect()
        }
        _ => return Err("the shard is neither a list nor an object".to_string()),
    };
    let mut events = Vec::new();
    let mut met: HashSet<String> = HashSet::new();
    for transcript in lists.into_iter().flatten() {
        let Value::Dict(fields) = transcript else {
            continue;
        };
        // `(t.get("stable_id") or "").split(".")[0]`.
        let id = match shard.get(fields, "stable_id") {
            Some(value) if value.truthy() => match shard.string(value) {
                Some(Some(id)) => id.split('.').next().unwrap_or("").to_string(),
                // Half a surrogate pair: the part before the first dot may still be text.
                Some(None) => match value {
                    Value::Str(from, to) => {
                        let raw = &text[*from..*to];
                        let head = raw.split(|byte| *byte == b'.').next().unwrap_or(b"");
                        match decode(head) {
                            Some(id) => id,
                            None => continue,
                        }
                    }
                    _ => continue,
                },
                None => return Err("a stable_id that is no string".to_string()),
            },
            _ => String::new(),
        };
        if !wanted.contains(&id) || met.contains(&id) {
            continue;
        }
        let event = match shard.get(fields, "variation_effect_feature_cache") {
            Some(cache) if cache.truthy() => match cache {
                Value::Dict(items) => {
                    let peptide = match shard.get(items, "peptide") {
                        None | Some(Value::Null) => Peptide::Absent,
                        Some(value) => match shard.string(value) {
                            Some(Some(text)) => Peptide::Text(text),
                            _ => Peptide::Other,
                        },
                    };
                    Event::Found {
                        id: id.clone(),
                        peptide,
                    }
                }
                _ => Event::Unreadable { id: id.clone() },
            },
            _ => Event::Found {
                id: id.clone(),
                peptide: Peptide::Absent,
            },
        };
        // An id whose cache could be read is settled for the rest of the shard.
        if matches!(event, Event::Found { .. }) {
            met.insert(id);
        }
        events.push(event);
    }
    Ok(events)
}

fn os_path(bytes: &[u8]) -> PathBuf {
    PathBuf::from(OsString::from_vec(bytes.to_vec()))
}

/// The names of a directory, whatever they name. A directory that cannot be listed
/// has none.
fn names(directory: &Path) -> Vec<Vec<u8>> {
    std::fs::read_dir(directory)
        .map(|listing| {
            listing
                .filter_map(Result::ok)
                .map(|entry| entry.file_name().as_bytes().to_vec())
                .collect()
        })
        .unwrap_or_default()
}

/// The shards of a cache: `transcripts/*/*.json` and `transcripts/*/*.json.gz`, where
/// the first `*` is any directory and the second any name, in the order of their paths'
/// components.
pub fn shards(cache: &[u8]) -> Vec<Vec<u8>> {
    let transcripts = crate::py::path_join(cache, b"transcripts");
    let mut found: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    for contig in names(&os_path(&transcripts)) {
        let directory = crate::py::path_join(&transcripts, &contig);
        if !os_path(&directory).is_dir() {
            continue;
        }
        for name in names(&os_path(&directory)) {
            if name.ends_with(b".json") || name.ends_with(b".json.gz") {
                found.push((contig.clone(), name));
            }
        }
    }
    found.sort();
    found
        .into_iter()
        .map(|(contig, name)| {
            crate::py::path_join(&crate::py::path_join(&transcripts, &contig), &name)
        })
        .collect()
}

/// The bytes of a shard: gunzipped, every member in turn, when its name ends in `.gz`.
pub fn read(path: &[u8]) -> Result<Vec<u8>, String> {
    let bytes = std::fs::read(os_path(path)).map_err(|error| crate::py::os_text(&error))?;
    if !path.ends_with(b".gz") {
        return Ok(bytes);
    }
    let mut out = Vec::new();
    MultiGzDecoder::new(bytes.as_slice())
        .read_to_end(&mut out)
        .map_err(|_| "the gzip stream is damaged or short".to_string())?;
    Ok(out)
}

#[cfg(test)]
#[path = "jsoncache_tests.rs"]
mod tests;
