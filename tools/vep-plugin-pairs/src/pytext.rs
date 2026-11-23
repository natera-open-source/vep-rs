// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! CPython's text rules, in which the comparison's rows, rules and reports are defined:
//! the white space of `str.strip()` and `str.split()`, `str.splitlines()`, `str.isdigit()`
//! and `str.isdecimal()`, `str.partition()`, `repr()` of a `str`, `int()` and `float()` of a
//! text, `bytes.decode("utf-8", "replace")`, `float.__repr__`, `str(pathlib.Path(arg))`, and
//! the `json` module's reading and writing of a value.
//!
//! Text is counted in code points, as a Python `str` is.

use std::fmt::Write;

/// `str.isspace()` of one character: the Unicode White_Space set plus U+001C to U+001F,
/// which Python counts through their bidirectional class and `char::is_whitespace` does not.
pub fn is_space(c: char) -> bool {
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

/// The body of a regular-expression class that matches what `\s` matches in a Python `str`
/// pattern, for use between `[` and `]` or after `[^`.
pub const RE_SPACE: &str =
    r"\t-\r\x1c-\x20\x{85}\x{a0}\x{1680}\x{2000}-\x{200a}\x{2028}\x{2029}\x{202f}\x{205f}\x{3000}";

/// `str.strip()` with no argument.
pub fn strip(s: &str) -> &str {
    s.trim_matches(is_space)
}

/// `str.lstrip()` with no argument.
pub fn lstrip(s: &str) -> &str {
    s.trim_start_matches(is_space)
}

/// `str.rstrip()` with no argument.
pub fn rstrip(s: &str) -> &str {
    s.trim_end_matches(is_space)
}

/// `str.split()` with no argument: the runs between white space, none of them empty.
pub fn split(s: &str) -> Vec<&str> {
    s.split(is_space).filter(|part| !part.is_empty()).collect()
}

/// `str.splitlines()`: the lines without their terminators, a `\r\n` pair one terminator,
/// and no empty last line after a final terminator.
pub fn splitlines(s: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = s.char_indices().peekable();
    while let Some((at, c)) = chars.next() {
        let breaks = matches!(
            c,
            '\n' | '\r' | '\u{0b}' | '\u{0c}' | '\u{1c}'
                ..='\u{1e}' | '\u{85}' | '\u{2028}' | '\u{2029}'
        );
        if !breaks {
            continue;
        }
        lines.push(&s[start..at]);
        start = at + c.len_utf8();
        if c == '\r' && chars.peek().is_some_and(|&(_, next)| next == '\n') {
            chars.next();
            start += 1;
        }
    }
    if start < s.len() {
        lines.push(&s[start..]);
    }
    lines
}

/// `str.partition(sep)`: the text before the first `sep`, `sep` where found (else empty)
/// and the text after it.
pub fn partition<'a>(s: &'a str, sep: &str) -> (&'a str, &'a str, &'a str) {
    match s.find(sep) {
        Some(at) => (&s[..at], &s[at..at + sep.len()], &s[at + sep.len()..]),
        None => (s, "", ""),
    }
}

/// `str.rpartition(sep)`: the text before the last `sep`, `sep` where found (else empty)
/// and the text after it; the whole text is the last part where `sep` is not found.
pub fn rpartition<'a>(s: &'a str, sep: &str) -> (&'a str, &'a str, &'a str) {
    match s.rfind(sep) {
        Some(at) => (&s[..at], &s[at..at + sep.len()], &s[at + sep.len()..]),
        None => ("", "", s),
    }
}

/// The value of a character `str.isdecimal()` takes, a decimal digit of any script.
pub fn decimal_value(c: char) -> Option<u32> {
    let code = c as u32;
    if code < 0x80 {
        return c.to_digit(10);
    }
    let run = DECIMAL_ZEROS.partition_point(|&zero| zero <= code);
    let zero = *DECIMAL_ZEROS.get(run.checked_sub(1)?)?;
    (code - zero < 10).then_some(code - zero)
}

/// `str.isdecimal()`: every character a decimal digit of some script, and at least one.
pub fn isdecimal(s: &str) -> bool {
    !s.is_empty() && s.chars().all(|c| decimal_value(c).is_some())
}

/// `str.isdigit()`: every character a decimal digit or another digit (a superscript, a
/// circled digit), and at least one.
pub fn isdigit(s: &str) -> bool {
    let other = |c: char| {
        let code = c as u32;
        let at = OTHER_DIGITS.partition_point(|&(first, _)| first <= code);
        at > 0 && code <= OTHER_DIGITS[at - 1].1
    };
    !s.is_empty() && s.chars().all(|c| decimal_value(c).is_some() || other(c))
}

/// `str.isascii() and str.isdigit()`: written in the digits `0` to `9` alone, at least one.
pub fn is_ascii_digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// The text `int()` and `float()` parse: every decimal digit as its ASCII digit, white space
/// above U+007E as a blank and any other character above it as `?`, then the ASCII blanks
/// around it dropped. U+001C to U+001F stay, so a number beside one is refused.
fn number_text(text: &str) -> String {
    let mapped: String = text
        .chars()
        .map(|c| match c {
            _ if (c as u32) < 127 => c,
            _ if is_space(c) => ' ',
            _ => decimal_value(c).map_or('?', |digit| char::from(b'0' + digit as u8)),
        })
        .collect();
    mapped
        .trim_matches(|c| matches!(c, ' ' | '\t'..='\r'))
        .to_string()
}

/// `text` without the underscores Python allows inside a number, each between two digits;
/// `None` where one stands elsewhere.
fn without_underscores(text: &str) -> Option<String> {
    let mut out = String::with_capacity(text.len());
    let mut previous = '\0';
    for c in text.chars() {
        if c == '_' {
            if !previous.is_ascii_digit() {
                return None;
            }
        } else {
            if previous == '_' && !c.is_ascii_digit() {
                return None;
            }
            out.push(c);
        }
        previous = c;
    }
    (previous != '_').then_some(out)
}

/// `int(text)`: white space around the number dropped, an optional sign, decimal digits of
/// any script with single underscores between them. `None` where Python raises
/// `ValueError`, which it does for a number of more than 4,300 digits whatever its value,
/// and for a value outside `i128`, which Python reads.
pub fn int(text: &str) -> Option<i128> {
    if is_ascii_digits(text) && text.len() <= 18 {
        return text.parse().ok();
    }
    let body = number_text(text);
    let (negative, digits) = match body.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, body.strip_prefix('+').unwrap_or(&body)),
    };
    let digits = without_underscores(digits)?;
    if !is_ascii_digits(&digits) || digits.len() > INT_MAX_DIGITS {
        return None;
    }
    let mut value: i128 = 0;
    for b in digits.bytes() {
        value = value.checked_mul(10)?.checked_sub(i128::from(b - b'0'))?;
    }
    if negative {
        Some(value)
    } else {
        value.checked_neg()
    }
}

/// `float(text)`: white space around the number dropped, an optional sign, then `inf`,
/// `infinity` or `nan` in any case, or decimal digits of any script with an optional
/// fraction and exponent and single underscores between digits; a number past the largest
/// double is an infinity. `None` where Python raises `ValueError`.
pub fn float(text: &str) -> Option<f64> {
    let body = without_underscores(&number_text(text))?;
    let unsigned = body.strip_prefix(['+', '-']).unwrap_or(&body);
    let word = unsigned.to_ascii_lowercase();
    if !matches!(word.as_str(), "inf" | "infinity" | "nan") {
        let (mantissa, exponent) = match unsigned.find(['e', 'E']) {
            Some(at) => (&unsigned[..at], Some(&unsigned[at + 1..])),
            None => (unsigned, None),
        };
        let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
        let digits = |part: &str| part.bytes().all(|b| b.is_ascii_digit());
        if !digits(whole) || !digits(fraction) || whole.len() + fraction.len() == 0 {
            return None;
        }
        if let Some(exponent) = exponent {
            if !is_ascii_digits(exponent.strip_prefix(['+', '-']).unwrap_or(exponent)) {
                return None;
            }
        }
    }
    body.parse().ok()
}

/// The text one lower than `digits` in its last place, where `x` lies exactly half way
/// between the two, the lower one ends in an even digit and it reads back as `x`. `digits`
/// is Rust's shortest text of `x` at the power of ten `exponent`, which is the upper of two
/// texts equally near.
fn even_text_of_a_tie(x: f64, digits: &str, exponent: i32) -> Option<String> {
    let last = *digits.as_bytes().last()?;
    // A text of fifteen digits is the only one of its length that reads back as its
    // double, so two texts of one length tie from sixteen digits on.
    if digits.len() < 16 || (last - b'0').is_multiple_of(2) {
        return None;
    }
    let lower = format!("{}{}", &digits[..digits.len() - 1], char::from(last - 1));
    // 767 significant digits write any double exactly.
    let exact = format!("{x:.766e}");
    let (exact_mantissa, _) = exact.split_once('e')?;
    let exact_digits: String = exact_mantissa.chars().filter(|c| *c != '.').collect();
    let half_way = exact_digits.trim_end_matches('0').strip_suffix('5') == Some(lower.as_str());
    let last_place = exponent - (digits.len() as i32 - 1);
    (half_way && format!("{lower}e{last_place}").parse() == Ok(x)).then_some(lower)
}

/// The significant digits `float.__repr__` writes for a finite `x` that is not negative, and
/// the power of ten of the first of them.
fn repr_digits(x: f64) -> (String, i32) {
    let shortest = format!("{x:e}");
    let (mantissa, exponent) = shortest.split_once('e').expect("LowerExp output");
    let exponent: i32 = exponent.parse().expect("LowerExp exponent");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let digits = even_text_of_a_tie(x, &digits, exponent).unwrap_or(digits);
    (digits, exponent)
}

/// `repr(x)` of a float: the shortest round-trip digits, of two such texts equally near
/// `x` the one that ends in an even digit, `.0` on an integral value, exponent form below
/// 1e-4 and from 1e16 with a signed two-digit exponent, and `nan`, `inf` and `-inf`.
pub fn float_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    if x.is_sign_negative() {
        return format!("-{}", float_repr(-x));
    }
    let (digits, exponent) = repr_digits(x);
    let decpt = exponent + 1;
    if decpt <= -4 || decpt > 16 {
        let (first, rest) = digits.split_at(1);
        let point = if rest.is_empty() { "" } else { "." };
        let sign = if exponent < 0 { '-' } else { '+' };
        let magnitude = exponent.abs();
        format!("{first}{point}{rest}e{sign}{magnitude:02}")
    } else if decpt <= 0 {
        format!("0.{}{}", "0".repeat(decpt.unsigned_abs() as usize), digits)
    } else if decpt as usize >= digits.len() {
        format!("{}{}.0", digits, "0".repeat(decpt as usize - digits.len()))
    } else {
        format!(
            "{}.{}",
            &digits[..decpt as usize],
            &digits[decpt as usize..]
        )
    }
}

/// `str.isprintable()` of one character.
fn is_printable(c: char) -> bool {
    let code = c as u32;
    if code < 0x80 {
        return (0x20..0x7f).contains(&code);
    }
    let at = UNPRINTABLE.partition_point(|&(first, _)| first <= code);
    at == 0 || code > UNPRINTABLE[at - 1].1
}

/// `repr()` of a `str`: single quotes, or double quotes where the text holds a single quote
/// and no double one; a backslash, the quote, `\n`, `\r` and `\t` escaped, and every
/// character `str.isprintable()` refuses as `\x`, `\u` or `\U` with its code in hexadecimal.
pub fn repr(s: &str) -> String {
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
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ if c == quote => {
                out.push('\\');
                out.push(c);
            }
            _ if is_printable(c) => out.push(c),
            _ => {
                let code = c as u32;
                if code < 0x100 {
                    write!(out, "\\x{code:02x}").unwrap();
                } else if code < 0x10000 {
                    write!(out, "\\u{code:04x}").unwrap();
                } else {
                    write!(out, "\\U{code:08x}").unwrap();
                }
            }
        }
    }
    out.push(quote);
    out
}

/// `bytes.decode("utf-8", "replace")`: one U+FFFD for each maximal part of an ill-formed
/// sequence, which is what `String::from_utf8_lossy` writes.
pub fn decode(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// `str(pathlib.PurePosixPath(arg))` on the argument's bytes, as Python applies it to the
/// surrogate-escaped `str` of any `argv` entry: repeated separators and `.` components
/// dropped, a trailing separator dropped, exactly two leading slashes kept, `..` left alone,
/// `.` for an empty result.
pub fn path_str(arg: &[u8]) -> Vec<u8> {
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

/// A value of the `json` module. An integer keeps its decimal text, since Python's has no
/// bound; an object keeps its members in the order a Python `dict` holds them.
#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Int(String),
    Float(f64),
    Str(String),
    List(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    /// An integer value.
    pub fn int(value: impl Into<i128>) -> Json {
        Json::Int(value.into().to_string())
    }

    /// A text value.
    pub fn str(value: impl Into<String>) -> Json {
        Json::Str(value.into())
    }

    /// An object of `members`, a later member of a name replacing the earlier one's value
    /// in its place, as assignments to a `dict` do.
    pub fn object<K: Into<String>>(members: impl IntoIterator<Item = (K, Json)>) -> Json {
        let mut object = Vec::new();
        for (name, value) in members {
            put(&mut object, name.into(), value);
        }
        Json::Object(object)
    }

    /// The member `name` of an object; `None` for another value or a missing member.
    pub fn get(&self, name: &str) -> Option<&Json> {
        match self {
            Json::Object(members) => members
                .iter()
                .find_map(|(key, value)| (key == name).then_some(value)),
            _ => None,
        }
    }
}

fn put(object: &mut Vec<(String, Json)>, name: String, value: Json) {
    match object.iter_mut().find(|(key, _)| *key == name) {
        Some(member) => member.1 = value,
        None => object.push((name, value)),
    }
}

/// Why `json_loads` reads no value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JsonError {
    /// The text is not JSON as `json.loads` reads it (its `ValueError`).
    Invalid,
    /// The value is nested deeper than `JSON_NESTING`.
    TooDeep,
}

/// The nesting past which a JSON text is refused: the depth at which CPython's default
/// recursion limit ends the walk that writes a parsed object back in canonical form.
pub const JSON_NESTING: usize = 497;

/// The most digits `int()` converts (`sys.get_int_max_str_digits()`): a text of more, an
/// integer of a JSON text among them, is a `ValueError`. The sign, the white space around
/// the number and the underscores between its digits are no digits.
const INT_MAX_DIGITS: usize = 4300;

/// `json.loads(text)`: the value of a JSON text with blanks around it, read as CPython
/// reads it. `NaN`, `Infinity` and `-Infinity` are numbers, a number past the largest
/// double is an infinity, an object's later member of a name replaces the earlier one's
/// value in its place, a control character inside a string is refused, and a `\u` escape
/// of half a surrogate pair is read as U+FFFD, which a Python `str` holds as the half.
pub fn json_loads(text: &str) -> Result<Json, JsonError> {
    let mut reader = JsonReader { text, at: 0 };
    reader.blanks();
    let value = reader.value(0)?;
    reader.blanks();
    if reader.at == text.len() {
        Ok(value)
    } else {
        Err(JsonError::Invalid)
    }
}

/// A reader over the text's bytes: every byte the grammar names is ASCII, so `at` stands
/// on a character boundary wherever a slice is taken.
struct JsonReader<'a> {
    text: &'a str,
    at: usize,
}

impl JsonReader<'_> {
    fn peek(&self) -> Option<u8> {
        self.text.as_bytes().get(self.at).copied()
    }

    fn blanks(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.at += 1;
        }
    }

    fn word(&mut self, word: &str) -> bool {
        let found = self.text.as_bytes()[self.at..].starts_with(word.as_bytes());
        if found {
            self.at += word.len();
        }
        found
    }

    fn expect(&mut self, byte: u8) -> Result<(), JsonError> {
        if self.peek() == Some(byte) {
            self.at += 1;
            Ok(())
        } else {
            Err(JsonError::Invalid)
        }
    }

    fn value(&mut self, depth: usize) -> Result<Json, JsonError> {
        match self.peek().ok_or(JsonError::Invalid)? {
            b'"' => {
                self.at += 1;
                self.string().map(Json::Str)
            }
            b'{' | b'[' if depth >= JSON_NESTING => Err(JsonError::TooDeep),
            b'{' => {
                self.at += 1;
                self.object(depth + 1)
            }
            b'[' => {
                self.at += 1;
                self.list(depth + 1)
            }
            b'n' if self.word("null") => Ok(Json::Null),
            b't' if self.word("true") => Ok(Json::Bool(true)),
            b'f' if self.word("false") => Ok(Json::Bool(false)),
            b'N' if self.word("NaN") => Ok(Json::Float(f64::NAN)),
            b'I' if self.word("Infinity") => Ok(Json::Float(f64::INFINITY)),
            b'-' if self.word("-Infinity") => Ok(Json::Float(f64::NEG_INFINITY)),
            _ => self.number(),
        }
    }

    fn object(&mut self, depth: usize) -> Result<Json, JsonError> {
        let mut members = Vec::new();
        self.blanks();
        if self.peek() == Some(b'}') {
            self.at += 1;
            return Ok(Json::Object(members));
        }
        loop {
            self.expect(b'"')?;
            let name = self.string()?;
            self.blanks();
            self.expect(b':')?;
            self.blanks();
            let value = self.value(depth)?;
            put(&mut members, name, value);
            self.blanks();
            if self.peek() == Some(b'}') {
                self.at += 1;
                return Ok(Json::Object(members));
            }
            self.expect(b',')?;
            self.blanks();
        }
    }

    fn list(&mut self, depth: usize) -> Result<Json, JsonError> {
        let mut items = Vec::new();
        self.blanks();
        if self.peek() == Some(b']') {
            self.at += 1;
            return Ok(Json::List(items));
        }
        loop {
            items.push(self.value(depth)?);
            self.blanks();
            if self.peek() == Some(b']') {
                self.at += 1;
                return Ok(Json::List(items));
            }
            self.expect(b',')?;
            self.blanks();
        }
    }

    fn hex4(&mut self) -> Result<u32, JsonError> {
        let mut code = 0;
        for _ in 0..4 {
            let digit = self
                .peek()
                .and_then(|byte| char::from(byte).to_digit(16))
                .ok_or(JsonError::Invalid)?;
            code = code * 16 + digit;
            self.at += 1;
        }
        Ok(code)
    }

    /// The text of a string whose opening quote is read.
    fn string(&mut self) -> Result<String, JsonError> {
        let mut out = String::new();
        loop {
            let run = self.at;
            while self
                .peek()
                .is_some_and(|byte| byte != b'"' && byte != b'\\' && byte >= 0x20)
            {
                self.at += 1;
            }
            out.push_str(&self.text[run..self.at]);
            let byte = self.peek().ok_or(JsonError::Invalid)?;
            self.at += 1;
            match byte {
                b'"' => return Ok(out),
                b'\\' => {}
                _ => return Err(JsonError::Invalid),
            }
            let escaped = self.peek().ok_or(JsonError::Invalid)?;
            self.at += 1;
            match escaped {
                b'"' | b'\\' | b'/' => out.push(char::from(escaped)),
                b'b' => out.push('\u{08}'),
                b'f' => out.push('\u{0c}'),
                b'n' => out.push('\n'),
                b'r' => out.push('\r'),
                b't' => out.push('\t'),
                b'u' => {
                    let mut code = self.hex4()?;
                    if (0xd800..0xdc00).contains(&code) && self.word("\\u") {
                        let back = self.at - 2;
                        match self.hex4() {
                            Ok(low) if (0xdc00..0xe000).contains(&low) => {
                                code = 0x10000 + ((code - 0xd800) << 10) + (low - 0xdc00);
                            }
                            _ => self.at = back,
                        }
                    }
                    out.push(char::from_u32(code).unwrap_or('\u{fffd}'));
                }
                _ => return Err(JsonError::Invalid),
            }
        }
    }

    fn digits(&mut self) -> usize {
        let start = self.at;
        while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
            self.at += 1;
        }
        self.at - start
    }

    fn number(&mut self) -> Result<Json, JsonError> {
        let start = self.at;
        if self.peek() == Some(b'-') {
            self.at += 1;
        }
        match self.peek() {
            Some(b'0') => self.at += 1,
            Some(b'1'..=b'9') => {
                self.digits();
            }
            _ => return Err(JsonError::Invalid),
        }
        let mut float = false;
        let bytes = self.text.as_bytes();
        if self.peek() == Some(b'.') && bytes.get(self.at + 1).is_some_and(u8::is_ascii_digit) {
            self.at += 1;
            self.digits();
            float = true;
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            let back = self.at;
            self.at += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.at += 1;
            }
            if self.digits() == 0 {
                self.at = back;
            } else {
                float = true;
            }
        }
        let written = &self.text[start..self.at];
        if float {
            return written.parse().map(Json::Float).or(Err(JsonError::Invalid));
        }
        let digits = written.trim_start_matches('-');
        if digits.len() > INT_MAX_DIGITS {
            return Err(JsonError::Invalid);
        }
        Ok(Json::Int(
            if digits == "0" { digits } else { written }.to_string(),
        ))
    }
}

/// Appends `json.dumps(s)` of a text: the quote, the backslash and the control characters
/// below U+0020 escaped, and under `ensure_ascii` every other character outside
/// `0x20..=0x7e` as `\u` and four lowercase hex digits, a surrogate pair above U+FFFF.
pub fn push_json_str(out: &mut String, s: &str, ensure_ascii: bool) {
    out.push('"');
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
            _ if !ensure_ascii && c >= ' ' => out.push(c),
            _ => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units) {
                    write!(out, "\\u{unit:04x}").unwrap();
                }
            }
        }
    }
    out.push('"');
}

fn push_dumped(out: &mut String, value: &Json, indent: Option<usize>, level: usize) {
    let line = |out: &mut String, level: usize| {
        if let Some(width) = indent {
            out.push('\n');
            out.push_str(&" ".repeat(width * level));
        }
    };
    let between = if indent.is_some() { "," } else { ", " };
    match value {
        Json::Null => out.push_str("null"),
        Json::Bool(true) => out.push_str("true"),
        Json::Bool(false) => out.push_str("false"),
        Json::Int(digits) => out.push_str(digits),
        Json::Float(x) if x.is_nan() => out.push_str("NaN"),
        Json::Float(x) if x.is_infinite() => {
            out.push_str(if *x < 0.0 { "-Infinity" } else { "Infinity" })
        }
        Json::Float(x) => out.push_str(&float_repr(*x)),
        Json::Str(text) => push_json_str(out, text, true),
        Json::List(items) if items.is_empty() => out.push_str("[]"),
        Json::Object(members) if members.is_empty() => out.push_str("{}"),
        Json::List(items) => {
            out.push('[');
            for (number, item) in items.iter().enumerate() {
                if number > 0 {
                    out.push_str(between);
                }
                line(out, level + 1);
                push_dumped(out, item, indent, level + 1);
            }
            line(out, level);
            out.push(']');
        }
        Json::Object(members) => {
            out.push('{');
            for (number, (name, member)) in members.iter().enumerate() {
                if number > 0 {
                    out.push_str(between);
                }
                line(out, level + 1);
                push_json_str(out, name, true);
                out.push_str(": ");
                push_dumped(out, member, indent, level + 1);
            }
            line(out, level);
            out.push('}');
        }
    }
}

/// `json.dumps(value)`: one line, `, ` between items and `: ` after a name, every character
/// outside ASCII escaped, the members in the object's order.
pub fn json_dumps(value: &Json) -> String {
    let mut out = String::new();
    push_dumped(&mut out, value, None, 0);
    out
}

/// `json.dumps(value, indent=2)`: an item or member per line, two blanks per level, every
/// character outside ASCII escaped, the members in the object's order.
pub fn json_dumps_indented(value: &Json) -> String {
    let mut out = String::new();
    push_dumped(&mut out, value, Some(2), 0);
    out
}

fn push_compact(out: &mut String, value: &Json) {
    match value {
        Json::Null => out.push_str("null"),
        Json::Bool(true) => out.push_str("true"),
        Json::Bool(false) => out.push_str("false"),
        Json::Int(digits) => out.push_str(digits),
        Json::Float(x) if !x.is_finite() => out.push_str("null"),
        Json::Float(x) => {
            let written = float_repr(*x);
            match written.split_once('e') {
                Some((mantissa, exponent)) => {
                    let exponent: i32 = exponent.parse().expect("repr exponent");
                    write!(out, "{mantissa}e{exponent}").unwrap();
                }
                None => out.push_str(&written),
            }
        }
        Json::Str(text) => push_json_str(out, text, false),
        Json::List(items) => {
            out.push('[');
            for (number, item) in items.iter().enumerate() {
                if number > 0 {
                    out.push(',');
                }
                push_compact(out, item);
            }
            out.push(']');
        }
        Json::Object(members) => {
            let mut sorted: Vec<&(String, Json)> = members.iter().collect();
            sorted.sort_by(|a, b| a.0.cmp(&b.0));
            out.push('{');
            for (number, (name, member)) in sorted.into_iter().enumerate() {
                if number > 0 {
                    out.push(',');
                }
                push_json_str(out, name, false);
                out.push(':');
                push_compact(out, member);
            }
            out.push('}');
        }
    }
}

/// A JSON value as compact text with object keys sorted by code point: no blank, a text as
/// it is but for the quote, the backslash and the control characters, a float as its
/// shortest round trip with the exponent unpadded, and `null` for one that is not finite.
pub fn json_compact(value: &Json) -> String {
    let mut out = String::new();
    push_compact(&mut out, value);
    out
}

// The three tables hold Unicode 16.0.0, the version CPython 3.14 reads text by.

/// The digit zero of every script's run of ten decimal digits, ascending.
const DECIMAL_ZEROS: [u32; 76] = [
    0x30, 0x660, 0x6f0, 0x7c0, 0x966, 0x9e6, 0xa66, 0xae6, 0xb66, 0xbe6, 0xc66, 0xce6, 0xd66,
    0xde6, 0xe50, 0xed0, 0xf20, 0x1040, 0x1090, 0x17e0, 0x1810, 0x1946, 0x19d0, 0x1a80, 0x1a90,
    0x1b50, 0x1bb0, 0x1c40, 0x1c50, 0xa620, 0xa8d0, 0xa900, 0xa9d0, 0xa9f0, 0xaa50, 0xabf0, 0xff10,
    0x104a0, 0x10d30, 0x10d40, 0x11066, 0x110f0, 0x11136, 0x111d0, 0x112f0, 0x11450, 0x114d0,
    0x11650, 0x116c0, 0x116d0, 0x116da, 0x11730, 0x118e0, 0x11950, 0x11bf0, 0x11c50, 0x11d50,
    0x11da0, 0x11f50, 0x16130, 0x16a60, 0x16ac0, 0x16b50, 0x16d70, 0x1ccf0, 0x1d7ce, 0x1d7d8,
    0x1d7e2, 0x1d7ec, 0x1d7f6, 0x1e140, 0x1e2f0, 0x1e4f0, 0x1e5f1, 0x1e950, 0x1fbf0,
];

/// The characters `str.isdigit()` takes and `str.isdecimal()` refuses, as ascending ranges.
#[rustfmt::skip]
const OTHER_DIGITS: [(u32, u32); 20] = [
    (0xb2, 0xb3), (0xb9, 0xb9), (0x1369, 0x1371), (0x19da, 0x19da),
    (0x2070, 0x2070), (0x2074, 0x2079), (0x2080, 0x2089), (0x2460, 0x2468),
    (0x2474, 0x247c), (0x2488, 0x2490), (0x24ea, 0x24ea), (0x24f5, 0x24fd),
    (0x24ff, 0x24ff), (0x2776, 0x277e), (0x2780, 0x2788), (0x278a, 0x2792),
    (0x10a40, 0x10a43), (0x10e60, 0x10e68), (0x11052, 0x1105a), (0x1f100, 0x1f10a),
];

/// The characters above U+007F that `str.isprintable()` refuses, as ascending ranges.
#[rustfmt::skip]
const UNPRINTABLE: [(u32, u32); 736] = [
    (0x80, 0xa0), (0xad, 0xad), (0x378, 0x379), (0x380, 0x383),
    (0x38b, 0x38b), (0x38d, 0x38d), (0x3a2, 0x3a2), (0x530, 0x530),
    (0x557, 0x558), (0x58b, 0x58c), (0x590, 0x590), (0x5c8, 0x5cf),
    (0x5eb, 0x5ee), (0x5f5, 0x605), (0x61c, 0x61c), (0x6dd, 0x6dd),
    (0x70e, 0x70f), (0x74b, 0x74c), (0x7b2, 0x7bf), (0x7fb, 0x7fc),
    (0x82e, 0x82f), (0x83f, 0x83f), (0x85c, 0x85d), (0x85f, 0x85f),
    (0x86b, 0x86f), (0x88f, 0x896), (0x8e2, 0x8e2), (0x984, 0x984),
    (0x98d, 0x98e), (0x991, 0x992), (0x9a9, 0x9a9), (0x9b1, 0x9b1),
    (0x9b3, 0x9b5), (0x9ba, 0x9bb), (0x9c5, 0x9c6), (0x9c9, 0x9ca),
    (0x9cf, 0x9d6), (0x9d8, 0x9db), (0x9de, 0x9de), (0x9e4, 0x9e5),
    (0x9ff, 0xa00), (0xa04, 0xa04), (0xa0b, 0xa0e), (0xa11, 0xa12),
    (0xa29, 0xa29), (0xa31, 0xa31), (0xa34, 0xa34), (0xa37, 0xa37),
    (0xa3a, 0xa3b), (0xa3d, 0xa3d), (0xa43, 0xa46), (0xa49, 0xa4a),
    (0xa4e, 0xa50), (0xa52, 0xa58), (0xa5d, 0xa5d), (0xa5f, 0xa65),
    (0xa77, 0xa80), (0xa84, 0xa84), (0xa8e, 0xa8e), (0xa92, 0xa92),
    (0xaa9, 0xaa9), (0xab1, 0xab1), (0xab4, 0xab4), (0xaba, 0xabb),
    (0xac6, 0xac6), (0xaca, 0xaca), (0xace, 0xacf), (0xad1, 0xadf),
    (0xae4, 0xae5), (0xaf2, 0xaf8), (0xb00, 0xb00), (0xb04, 0xb04),
    (0xb0d, 0xb0e), (0xb11, 0xb12), (0xb29, 0xb29), (0xb31, 0xb31),
    (0xb34, 0xb34), (0xb3a, 0xb3b), (0xb45, 0xb46), (0xb49, 0xb4a),
    (0xb4e, 0xb54), (0xb58, 0xb5b), (0xb5e, 0xb5e), (0xb64, 0xb65),
    (0xb78, 0xb81), (0xb84, 0xb84), (0xb8b, 0xb8d), (0xb91, 0xb91),
    (0xb96, 0xb98), (0xb9b, 0xb9b), (0xb9d, 0xb9d), (0xba0, 0xba2),
    (0xba5, 0xba7), (0xbab, 0xbad), (0xbba, 0xbbd), (0xbc3, 0xbc5),
    (0xbc9, 0xbc9), (0xbce, 0xbcf), (0xbd1, 0xbd6), (0xbd8, 0xbe5),
    (0xbfb, 0xbff), (0xc0d, 0xc0d), (0xc11, 0xc11), (0xc29, 0xc29),
    (0xc3a, 0xc3b), (0xc45, 0xc45), (0xc49, 0xc49), (0xc4e, 0xc54),
    (0xc57, 0xc57), (0xc5b, 0xc5c), (0xc5e, 0xc5f), (0xc64, 0xc65),
    (0xc70, 0xc76), (0xc8d, 0xc8d), (0xc91, 0xc91), (0xca9, 0xca9),
    (0xcb4, 0xcb4), (0xcba, 0xcbb), (0xcc5, 0xcc5), (0xcc9, 0xcc9),
    (0xcce, 0xcd4), (0xcd7, 0xcdc), (0xcdf, 0xcdf), (0xce4, 0xce5),
    (0xcf0, 0xcf0), (0xcf4, 0xcff), (0xd0d, 0xd0d), (0xd11, 0xd11),
    (0xd45, 0xd45), (0xd49, 0xd49), (0xd50, 0xd53), (0xd64, 0xd65),
    (0xd80, 0xd80), (0xd84, 0xd84), (0xd97, 0xd99), (0xdb2, 0xdb2),
    (0xdbc, 0xdbc), (0xdbe, 0xdbf), (0xdc7, 0xdc9), (0xdcb, 0xdce),
    (0xdd5, 0xdd5), (0xdd7, 0xdd7), (0xde0, 0xde5), (0xdf0, 0xdf1),
    (0xdf5, 0xe00), (0xe3b, 0xe3e), (0xe5c, 0xe80), (0xe83, 0xe83),
    (0xe85, 0xe85), (0xe8b, 0xe8b), (0xea4, 0xea4), (0xea6, 0xea6),
    (0xebe, 0xebf), (0xec5, 0xec5), (0xec7, 0xec7), (0xecf, 0xecf),
    (0xeda, 0xedb), (0xee0, 0xeff), (0xf48, 0xf48), (0xf6d, 0xf70),
    (0xf98, 0xf98), (0xfbd, 0xfbd), (0xfcd, 0xfcd), (0xfdb, 0xfff),
    (0x10c6, 0x10c6), (0x10c8, 0x10cc), (0x10ce, 0x10cf), (0x1249, 0x1249),
    (0x124e, 0x124f), (0x1257, 0x1257), (0x1259, 0x1259), (0x125e, 0x125f),
    (0x1289, 0x1289), (0x128e, 0x128f), (0x12b1, 0x12b1), (0x12b6, 0x12b7),
    (0x12bf, 0x12bf), (0x12c1, 0x12c1), (0x12c6, 0x12c7), (0x12d7, 0x12d7),
    (0x1311, 0x1311), (0x1316, 0x1317), (0x135b, 0x135c), (0x137d, 0x137f),
    (0x139a, 0x139f), (0x13f6, 0x13f7), (0x13fe, 0x13ff), (0x1680, 0x1680),
    (0x169d, 0x169f), (0x16f9, 0x16ff), (0x1716, 0x171e), (0x1737, 0x173f),
    (0x1754, 0x175f), (0x176d, 0x176d), (0x1771, 0x1771), (0x1774, 0x177f),
    (0x17de, 0x17df), (0x17ea, 0x17ef), (0x17fa, 0x17ff), (0x180e, 0x180e),
    (0x181a, 0x181f), (0x1879, 0x187f), (0x18ab, 0x18af), (0x18f6, 0x18ff),
    (0x191f, 0x191f), (0x192c, 0x192f), (0x193c, 0x193f), (0x1941, 0x1943),
    (0x196e, 0x196f), (0x1975, 0x197f), (0x19ac, 0x19af), (0x19ca, 0x19cf),
    (0x19db, 0x19dd), (0x1a1c, 0x1a1d), (0x1a5f, 0x1a5f), (0x1a7d, 0x1a7e),
    (0x1a8a, 0x1a8f), (0x1a9a, 0x1a9f), (0x1aae, 0x1aaf), (0x1acf, 0x1aff),
    (0x1b4d, 0x1b4d), (0x1bf4, 0x1bfb), (0x1c38, 0x1c3a), (0x1c4a, 0x1c4c),
    (0x1c8b, 0x1c8f), (0x1cbb, 0x1cbc), (0x1cc8, 0x1ccf), (0x1cfb, 0x1cff),
    (0x1f16, 0x1f17), (0x1f1e, 0x1f1f), (0x1f46, 0x1f47), (0x1f4e, 0x1f4f),
    (0x1f58, 0x1f58), (0x1f5a, 0x1f5a), (0x1f5c, 0x1f5c), (0x1f5e, 0x1f5e),
    (0x1f7e, 0x1f7f), (0x1fb5, 0x1fb5), (0x1fc5, 0x1fc5), (0x1fd4, 0x1fd5),
    (0x1fdc, 0x1fdc), (0x1ff0, 0x1ff1), (0x1ff5, 0x1ff5), (0x1fff, 0x200f),
    (0x2028, 0x202f), (0x205f, 0x206f), (0x2072, 0x2073), (0x208f, 0x208f),
    (0x209d, 0x209f), (0x20c1, 0x20cf), (0x20f1, 0x20ff), (0x218c, 0x218f),
    (0x242a, 0x243f), (0x244b, 0x245f), (0x2b74, 0x2b75), (0x2b96, 0x2b96),
    (0x2cf4, 0x2cf8), (0x2d26, 0x2d26), (0x2d28, 0x2d2c), (0x2d2e, 0x2d2f),
    (0x2d68, 0x2d6e), (0x2d71, 0x2d7e), (0x2d97, 0x2d9f), (0x2da7, 0x2da7),
    (0x2daf, 0x2daf), (0x2db7, 0x2db7), (0x2dbf, 0x2dbf), (0x2dc7, 0x2dc7),
    (0x2dcf, 0x2dcf), (0x2dd7, 0x2dd7), (0x2ddf, 0x2ddf), (0x2e5e, 0x2e7f),
    (0x2e9a, 0x2e9a), (0x2ef4, 0x2eff), (0x2fd6, 0x2fef), (0x3000, 0x3000),
    (0x3040, 0x3040), (0x3097, 0x3098), (0x3100, 0x3104), (0x3130, 0x3130),
    (0x318f, 0x318f), (0x31e6, 0x31ee), (0x321f, 0x321f), (0xa48d, 0xa48f),
    (0xa4c7, 0xa4cf), (0xa62c, 0xa63f), (0xa6f8, 0xa6ff), (0xa7ce, 0xa7cf),
    (0xa7d2, 0xa7d2), (0xa7d4, 0xa7d4), (0xa7dd, 0xa7f1), (0xa82d, 0xa82f),
    (0xa83a, 0xa83f), (0xa878, 0xa87f), (0xa8c6, 0xa8cd), (0xa8da, 0xa8df),
    (0xa954, 0xa95e), (0xa97d, 0xa97f), (0xa9ce, 0xa9ce), (0xa9da, 0xa9dd),
    (0xa9ff, 0xa9ff), (0xaa37, 0xaa3f), (0xaa4e, 0xaa4f), (0xaa5a, 0xaa5b),
    (0xaac3, 0xaada), (0xaaf7, 0xab00), (0xab07, 0xab08), (0xab0f, 0xab10),
    (0xab17, 0xab1f), (0xab27, 0xab27), (0xab2f, 0xab2f), (0xab6c, 0xab6f),
    (0xabee, 0xabef), (0xabfa, 0xabff), (0xd7a4, 0xd7af), (0xd7c7, 0xd7ca),
    (0xd7fc, 0xf8ff), (0xfa6e, 0xfa6f), (0xfada, 0xfaff), (0xfb07, 0xfb12),
    (0xfb18, 0xfb1c), (0xfb37, 0xfb37), (0xfb3d, 0xfb3d), (0xfb3f, 0xfb3f),
    (0xfb42, 0xfb42), (0xfb45, 0xfb45), (0xfbc3, 0xfbd2), (0xfd90, 0xfd91),
    (0xfdc8, 0xfdce), (0xfdd0, 0xfdef), (0xfe1a, 0xfe1f), (0xfe53, 0xfe53),
    (0xfe67, 0xfe67), (0xfe6c, 0xfe6f), (0xfe75, 0xfe75), (0xfefd, 0xff00),
    (0xffbf, 0xffc1), (0xffc8, 0xffc9), (0xffd0, 0xffd1), (0xffd8, 0xffd9),
    (0xffdd, 0xffdf), (0xffe7, 0xffe7), (0xffef, 0xfffb), (0xfffe, 0xffff),
    (0x1000c, 0x1000c), (0x10027, 0x10027), (0x1003b, 0x1003b), (0x1003e, 0x1003e),
    (0x1004e, 0x1004f), (0x1005e, 0x1007f), (0x100fb, 0x100ff), (0x10103, 0x10106),
    (0x10134, 0x10136), (0x1018f, 0x1018f), (0x1019d, 0x1019f), (0x101a1, 0x101cf),
    (0x101fe, 0x1027f), (0x1029d, 0x1029f), (0x102d1, 0x102df), (0x102fc, 0x102ff),
    (0x10324, 0x1032c), (0x1034b, 0x1034f), (0x1037b, 0x1037f), (0x1039e, 0x1039e),
    (0x103c4, 0x103c7), (0x103d6, 0x103ff), (0x1049e, 0x1049f), (0x104aa, 0x104af),
    (0x104d4, 0x104d7), (0x104fc, 0x104ff), (0x10528, 0x1052f), (0x10564, 0x1056e),
    (0x1057b, 0x1057b), (0x1058b, 0x1058b), (0x10593, 0x10593), (0x10596, 0x10596),
    (0x105a2, 0x105a2), (0x105b2, 0x105b2), (0x105ba, 0x105ba), (0x105bd, 0x105bf),
    (0x105f4, 0x105ff), (0x10737, 0x1073f), (0x10756, 0x1075f), (0x10768, 0x1077f),
    (0x10786, 0x10786), (0x107b1, 0x107b1), (0x107bb, 0x107ff), (0x10806, 0x10807),
    (0x10809, 0x10809), (0x10836, 0x10836), (0x10839, 0x1083b), (0x1083d, 0x1083e),
    (0x10856, 0x10856), (0x1089f, 0x108a6), (0x108b0, 0x108df), (0x108f3, 0x108f3),
    (0x108f6, 0x108fa), (0x1091c, 0x1091e), (0x1093a, 0x1093e), (0x10940, 0x1097f),
    (0x109b8, 0x109bb), (0x109d0, 0x109d1), (0x10a04, 0x10a04), (0x10a07, 0x10a0b),
    (0x10a14, 0x10a14), (0x10a18, 0x10a18), (0x10a36, 0x10a37), (0x10a3b, 0x10a3e),
    (0x10a49, 0x10a4f), (0x10a59, 0x10a5f), (0x10aa0, 0x10abf), (0x10ae7, 0x10aea),
    (0x10af7, 0x10aff), (0x10b36, 0x10b38), (0x10b56, 0x10b57), (0x10b73, 0x10b77),
    (0x10b92, 0x10b98), (0x10b9d, 0x10ba8), (0x10bb0, 0x10bff), (0x10c49, 0x10c7f),
    (0x10cb3, 0x10cbf), (0x10cf3, 0x10cf9), (0x10d28, 0x10d2f), (0x10d3a, 0x10d3f),
    (0x10d66, 0x10d68), (0x10d86, 0x10d8d), (0x10d90, 0x10e5f), (0x10e7f, 0x10e7f),
    (0x10eaa, 0x10eaa), (0x10eae, 0x10eaf), (0x10eb2, 0x10ec1), (0x10ec5, 0x10efb),
    (0x10f28, 0x10f2f), (0x10f5a, 0x10f6f), (0x10f8a, 0x10faf), (0x10fcc, 0x10fdf),
    (0x10ff7, 0x10fff), (0x1104e, 0x11051), (0x11076, 0x1107e), (0x110bd, 0x110bd),
    (0x110c3, 0x110cf), (0x110e9, 0x110ef), (0x110fa, 0x110ff), (0x11135, 0x11135),
    (0x11148, 0x1114f), (0x11177, 0x1117f), (0x111e0, 0x111e0), (0x111f5, 0x111ff),
    (0x11212, 0x11212), (0x11242, 0x1127f), (0x11287, 0x11287), (0x11289, 0x11289),
    (0x1128e, 0x1128e), (0x1129e, 0x1129e), (0x112aa, 0x112af), (0x112eb, 0x112ef),
    (0x112fa, 0x112ff), (0x11304, 0x11304), (0x1130d, 0x1130e), (0x11311, 0x11312),
    (0x11329, 0x11329), (0x11331, 0x11331), (0x11334, 0x11334), (0x1133a, 0x1133a),
    (0x11345, 0x11346), (0x11349, 0x1134a), (0x1134e, 0x1134f), (0x11351, 0x11356),
    (0x11358, 0x1135c), (0x11364, 0x11365), (0x1136d, 0x1136f), (0x11375, 0x1137f),
    (0x1138a, 0x1138a), (0x1138c, 0x1138d), (0x1138f, 0x1138f), (0x113b6, 0x113b6),
    (0x113c1, 0x113c1), (0x113c3, 0x113c4), (0x113c6, 0x113c6), (0x113cb, 0x113cb),
    (0x113d6, 0x113d6), (0x113d9, 0x113e0), (0x113e3, 0x113ff), (0x1145c, 0x1145c),
    (0x11462, 0x1147f), (0x114c8, 0x114cf), (0x114da, 0x1157f), (0x115b6, 0x115b7),
    (0x115de, 0x115ff), (0x11645, 0x1164f), (0x1165a, 0x1165f), (0x1166d, 0x1167f),
    (0x116ba, 0x116bf), (0x116ca, 0x116cf), (0x116e4, 0x116ff), (0x1171b, 0x1171c),
    (0x1172c, 0x1172f), (0x11747, 0x117ff), (0x1183c, 0x1189f), (0x118f3, 0x118fe),
    (0x11907, 0x11908), (0x1190a, 0x1190b), (0x11914, 0x11914), (0x11917, 0x11917),
    (0x11936, 0x11936), (0x11939, 0x1193a), (0x11947, 0x1194f), (0x1195a, 0x1199f),
    (0x119a8, 0x119a9), (0x119d8, 0x119d9), (0x119e5, 0x119ff), (0x11a48, 0x11a4f),
    (0x11aa3, 0x11aaf), (0x11af9, 0x11aff), (0x11b0a, 0x11bbf), (0x11be2, 0x11bef),
    (0x11bfa, 0x11bff), (0x11c09, 0x11c09), (0x11c37, 0x11c37), (0x11c46, 0x11c4f),
    (0x11c6d, 0x11c6f), (0x11c90, 0x11c91), (0x11ca8, 0x11ca8), (0x11cb7, 0x11cff),
    (0x11d07, 0x11d07), (0x11d0a, 0x11d0a), (0x11d37, 0x11d39), (0x11d3b, 0x11d3b),
    (0x11d3e, 0x11d3e), (0x11d48, 0x11d4f), (0x11d5a, 0x11d5f), (0x11d66, 0x11d66),
    (0x11d69, 0x11d69), (0x11d8f, 0x11d8f), (0x11d92, 0x11d92), (0x11d99, 0x11d9f),
    (0x11daa, 0x11edf), (0x11ef9, 0x11eff), (0x11f11, 0x11f11), (0x11f3b, 0x11f3d),
    (0x11f5b, 0x11faf), (0x11fb1, 0x11fbf), (0x11ff2, 0x11ffe), (0x1239a, 0x123ff),
    (0x1246f, 0x1246f), (0x12475, 0x1247f), (0x12544, 0x12f8f), (0x12ff3, 0x12fff),
    (0x13430, 0x1343f), (0x13456, 0x1345f), (0x143fb, 0x143ff), (0x14647, 0x160ff),
    (0x1613a, 0x167ff), (0x16a39, 0x16a3f), (0x16a5f, 0x16a5f), (0x16a6a, 0x16a6d),
    (0x16abf, 0x16abf), (0x16aca, 0x16acf), (0x16aee, 0x16aef), (0x16af6, 0x16aff),
    (0x16b46, 0x16b4f), (0x16b5a, 0x16b5a), (0x16b62, 0x16b62), (0x16b78, 0x16b7c),
    (0x16b90, 0x16d3f), (0x16d7a, 0x16e3f), (0x16e9b, 0x16eff), (0x16f4b, 0x16f4e),
    (0x16f88, 0x16f8e), (0x16fa0, 0x16fdf), (0x16fe5, 0x16fef), (0x16ff2, 0x16fff),
    (0x187f8, 0x187ff), (0x18cd6, 0x18cfe), (0x18d09, 0x1afef), (0x1aff4, 0x1aff4),
    (0x1affc, 0x1affc), (0x1afff, 0x1afff), (0x1b123, 0x1b131), (0x1b133, 0x1b14f),
    (0x1b153, 0x1b154), (0x1b156, 0x1b163), (0x1b168, 0x1b16f), (0x1b2fc, 0x1bbff),
    (0x1bc6b, 0x1bc6f), (0x1bc7d, 0x1bc7f), (0x1bc89, 0x1bc8f), (0x1bc9a, 0x1bc9b),
    (0x1bca0, 0x1cbff), (0x1ccfa, 0x1ccff), (0x1ceb4, 0x1ceff), (0x1cf2e, 0x1cf2f),
    (0x1cf47, 0x1cf4f), (0x1cfc4, 0x1cfff), (0x1d0f6, 0x1d0ff), (0x1d127, 0x1d128),
    (0x1d173, 0x1d17a), (0x1d1eb, 0x1d1ff), (0x1d246, 0x1d2bf), (0x1d2d4, 0x1d2df),
    (0x1d2f4, 0x1d2ff), (0x1d357, 0x1d35f), (0x1d379, 0x1d3ff), (0x1d455, 0x1d455),
    (0x1d49d, 0x1d49d), (0x1d4a0, 0x1d4a1), (0x1d4a3, 0x1d4a4), (0x1d4a7, 0x1d4a8),
    (0x1d4ad, 0x1d4ad), (0x1d4ba, 0x1d4ba), (0x1d4bc, 0x1d4bc), (0x1d4c4, 0x1d4c4),
    (0x1d506, 0x1d506), (0x1d50b, 0x1d50c), (0x1d515, 0x1d515), (0x1d51d, 0x1d51d),
    (0x1d53a, 0x1d53a), (0x1d53f, 0x1d53f), (0x1d545, 0x1d545), (0x1d547, 0x1d549),
    (0x1d551, 0x1d551), (0x1d6a6, 0x1d6a7), (0x1d7cc, 0x1d7cd), (0x1da8c, 0x1da9a),
    (0x1daa0, 0x1daa0), (0x1dab0, 0x1deff), (0x1df1f, 0x1df24), (0x1df2b, 0x1dfff),
    (0x1e007, 0x1e007), (0x1e019, 0x1e01a), (0x1e022, 0x1e022), (0x1e025, 0x1e025),
    (0x1e02b, 0x1e02f), (0x1e06e, 0x1e08e), (0x1e090, 0x1e0ff), (0x1e12d, 0x1e12f),
    (0x1e13e, 0x1e13f), (0x1e14a, 0x1e14d), (0x1e150, 0x1e28f), (0x1e2af, 0x1e2bf),
    (0x1e2fa, 0x1e2fe), (0x1e300, 0x1e4cf), (0x1e4fa, 0x1e5cf), (0x1e5fb, 0x1e5fe),
    (0x1e600, 0x1e7df), (0x1e7e7, 0x1e7e7), (0x1e7ec, 0x1e7ec), (0x1e7ef, 0x1e7ef),
    (0x1e7ff, 0x1e7ff), (0x1e8c5, 0x1e8c6), (0x1e8d7, 0x1e8ff), (0x1e94c, 0x1e94f),
    (0x1e95a, 0x1e95d), (0x1e960, 0x1ec70), (0x1ecb5, 0x1ed00), (0x1ed3e, 0x1edff),
    (0x1ee04, 0x1ee04), (0x1ee20, 0x1ee20), (0x1ee23, 0x1ee23), (0x1ee25, 0x1ee26),
    (0x1ee28, 0x1ee28), (0x1ee33, 0x1ee33), (0x1ee38, 0x1ee38), (0x1ee3a, 0x1ee3a),
    (0x1ee3c, 0x1ee41), (0x1ee43, 0x1ee46), (0x1ee48, 0x1ee48), (0x1ee4a, 0x1ee4a),
    (0x1ee4c, 0x1ee4c), (0x1ee50, 0x1ee50), (0x1ee53, 0x1ee53), (0x1ee55, 0x1ee56),
    (0x1ee58, 0x1ee58), (0x1ee5a, 0x1ee5a), (0x1ee5c, 0x1ee5c), (0x1ee5e, 0x1ee5e),
    (0x1ee60, 0x1ee60), (0x1ee63, 0x1ee63), (0x1ee65, 0x1ee66), (0x1ee6b, 0x1ee6b),
    (0x1ee73, 0x1ee73), (0x1ee78, 0x1ee78), (0x1ee7d, 0x1ee7d), (0x1ee7f, 0x1ee7f),
    (0x1ee8a, 0x1ee8a), (0x1ee9c, 0x1eea0), (0x1eea4, 0x1eea4), (0x1eeaa, 0x1eeaa),
    (0x1eebc, 0x1eeef), (0x1eef2, 0x1efff), (0x1f02c, 0x1f02f), (0x1f094, 0x1f09f),
    (0x1f0af, 0x1f0b0), (0x1f0c0, 0x1f0c0), (0x1f0d0, 0x1f0d0), (0x1f0f6, 0x1f0ff),
    (0x1f1ae, 0x1f1e5), (0x1f203, 0x1f20f), (0x1f23c, 0x1f23f), (0x1f249, 0x1f24f),
    (0x1f252, 0x1f25f), (0x1f266, 0x1f2ff), (0x1f6d8, 0x1f6db), (0x1f6ed, 0x1f6ef),
    (0x1f6fd, 0x1f6ff), (0x1f777, 0x1f77a), (0x1f7da, 0x1f7df), (0x1f7ec, 0x1f7ef),
    (0x1f7f1, 0x1f7ff), (0x1f80c, 0x1f80f), (0x1f848, 0x1f84f), (0x1f85a, 0x1f85f),
    (0x1f888, 0x1f88f), (0x1f8ae, 0x1f8af), (0x1f8bc, 0x1f8bf), (0x1f8c2, 0x1f8ff),
    (0x1fa54, 0x1fa5f), (0x1fa6e, 0x1fa6f), (0x1fa7d, 0x1fa7f), (0x1fa8a, 0x1fa8e),
    (0x1fac7, 0x1facd), (0x1fadd, 0x1fade), (0x1faea, 0x1faef), (0x1faf9, 0x1faff),
    (0x1fb93, 0x1fb93), (0x1fbfa, 0x1ffff), (0x2a6e0, 0x2a6ff), (0x2b73a, 0x2b73f),
    (0x2b81e, 0x2b81f), (0x2cea2, 0x2ceaf), (0x2ebe1, 0x2ebef), (0x2ee5e, 0x2f7ff),
    (0x2fa1e, 0x2ffff), (0x3134b, 0x3134f), (0x323b0, 0xe00ff), (0xe01f0, 0x10ffff),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn white_space_is_pythons_isspace() {
        let expected: Vec<u32> = [
            0x9, 0xa, 0xb, 0xc, 0xd, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x85, 0xa0, 0x1680, 0x2000,
            0x2001, 0x2002, 0x2003, 0x2004, 0x2005, 0x2006, 0x2007, 0x2008, 0x2009, 0x200a, 0x2028,
            0x2029, 0x202f, 0x205f, 0x3000,
        ]
        .to_vec();
        let actual: Vec<u32> = (0..=0x10ffffu32)
            .filter_map(char::from_u32)
            .filter(|c| is_space(*c))
            .map(|c| c as u32)
            .collect();
        assert_eq!(actual, expected);
        assert_eq!(strip("\u{1f}A/T\u{a0}\u{3000}"), "A/T");
        assert_eq!(strip(" \t x \u{200b}"), "x \u{200b}");
        assert_eq!(lstrip("\u{85} x "), "x ");
        assert_eq!(rstrip(" x\u{2028}\n"), " x");
    }

    #[test]
    fn split_without_a_separator_keeps_the_runs_between_white_space() {
        assert_eq!(split("  a  b "), ["a", "b"]);
        assert_eq!(split("a\u{1c}b\u{a0}c\u{3000}d"), ["a", "b", "c", "d"]);
        assert!(split("").is_empty() && split("   ").is_empty());
        assert_eq!(split("\ta\nb\u{0b}c\u{0c}d\re"), ["a", "b", "c", "d", "e"]);
    }

    #[test]
    fn splitlines_breaks_where_python_breaks() {
        assert_eq!(splitlines("a\nb\r\nc\rd"), ["a", "b", "c", "d"]);
        assert_eq!(splitlines("a\n"), ["a"]);
        assert_eq!(splitlines("\n"), [""]);
        assert!(splitlines("").is_empty());
        assert_eq!(splitlines("a\u{0b}b\u{0c}c"), ["a", "b", "c"]);
        assert_eq!(
            splitlines("a\u{1c}b\u{1d}b\u{1e}b\u{1f}b"),
            ["a", "b", "b", "b\u{1f}b"]
        );
        assert_eq!(
            splitlines("a\u{85}b\u{2028}c\u{2029}d"),
            ["a", "b", "c", "d"]
        );
        assert_eq!(splitlines("a\r\n\r\nb"), ["a", "", "b"]);
        assert_eq!(splitlines("a\n\rb"), ["a", "", "b"]);
    }

    #[test]
    fn partition_and_rpartition_part_at_the_first_and_the_last_separator() {
        assert_eq!(partition("a=b=c", "="), ("a", "=", "b=c"));
        assert_eq!(partition("abc", "="), ("abc", "", ""));
        assert_eq!(partition("=a", "="), ("", "=", "a"));
        assert_eq!(rpartition("a=b=c", "="), ("a=b", "=", "c"));
        assert_eq!(rpartition("abc", "="), ("", "", "abc"));
        assert_eq!(partition("a::b", "::"), ("a", "::", "b"));
    }

    #[test]
    fn digits_of_every_script_are_decimal_and_superscripts_are_digits_only() {
        assert!(isdecimal("0123456789") && isdigit("0123456789"));
        assert!(isdecimal("\u{661}\u{662}") && isdecimal("\u{ff11}") && isdecimal("\u{1fbf9}"));
        assert!(!isdecimal("\u{b2}") && isdigit("\u{b2}") && isdigit("\u{2460}"));
        assert!(!isdecimal("") && !isdigit("") && !isdigit("1a") && !isdigit("-1"));
        assert!(!isdigit("\u{bd}") && !isdecimal("\u{2180}"));
        assert_eq!(decimal_value('\u{666}'), Some(6));
        assert_eq!(decimal_value('\u{1d7ff}'), Some(9));
        assert_eq!(decimal_value('\u{1d800}'), None);
        assert!(is_ascii_digits("042") && !is_ascii_digits("\u{661}") && !is_ascii_digits(""));
        let decimal = (0..=0x10ffffu32)
            .filter_map(char::from_u32)
            .filter(|c| decimal_value(*c).is_some())
            .count();
        assert_eq!(decimal, 760);
    }

    #[test]
    fn int_reads_what_python_reads() {
        assert_eq!(int(" 1 "), Some(1));
        assert_eq!(int("\u{a0}1\u{3000}"), Some(1));
        assert_eq!(int("1_000"), Some(1000));
        assert_eq!(int("+5"), Some(5));
        assert_eq!(int("-0"), Some(0));
        assert_eq!(int("\u{661}\u{662}"), Some(12));
        assert_eq!(int("\u{ff11}"), Some(1));
        assert_eq!(
            int("-170141183460469231731687303715884105728"),
            Some(i128::MIN)
        );
        assert_eq!(
            int("170141183460469231731687303715884105727"),
            Some(i128::MAX)
        );
        for refused in [
            "\u{1c}1",
            "1\u{1f}",
            "_1",
            "1_",
            "1__0",
            "0x10",
            "\u{b2}",
            "",
            " ",
            "1.0",
            "1e3",
            "+ 5",
            "--5",
            "170141183460469231731687303715884105728",
        ] {
            assert_eq!(int(refused), None, "{refused:?}");
        }
    }

    #[test]
    fn int_refuses_a_number_of_more_digits_than_python_converts() {
        let zeros = |count: usize| "0".repeat(count);
        assert_eq!(int(&zeros(4300)), Some(0));
        assert_eq!(int(&zeros(4301)), None);
        assert_eq!(int(&format!("{}7", zeros(4299))), Some(7));
        assert_eq!(int(&format!("{}7", zeros(4300))), None);
        // 4,300 digits with the 39 of the largest value, and one more.
        assert_eq!(
            int(&format!("{}{}", zeros(4261), i128::MAX)),
            Some(i128::MAX)
        );
        assert_eq!(int(&format!("{}{}", zeros(4262), i128::MAX)), None);
        assert_eq!(
            int(&format!("-{}{}", zeros(4261), i128::MIN.unsigned_abs())),
            Some(i128::MIN)
        );
        // The sign, the white space and the underscores are not counted.
        assert_eq!(int(&format!(" -{}7\n", zeros(4299))), Some(-7));
        assert_eq!(int(&format!("+{}7", zeros(4300))), None);
        assert_eq!(int(&format!("{}7", "0_".repeat(4299))), Some(7));
        assert_eq!(int(&format!("{}7", "0_".repeat(4300))), None);
        // The digits of another script are.
        assert_eq!(int(&format!("{}7", "\u{660}".repeat(4299))), Some(7));
        assert_eq!(int(&"\u{ff10}".repeat(4301)), None);
    }

    #[test]
    fn float_reads_what_python_reads() {
        let read = |text: &str| float(text).map(float_repr);
        for (text, expected) in [
            (" 1 ", "1.0"),
            ("\u{a0}1.5\u{3000}", "1.5"),
            ("1_000.5", "1000.5"),
            ("+5", "5.0"),
            ("-0", "-0.0"),
            ("\u{661}.\u{662}", "1.2"),
            ("inf", "inf"),
            ("-Infinity", "-inf"),
            ("nan", "nan"),
            ("NaN", "nan"),
            ("1e999", "inf"),
            ("1.", "1.0"),
            (".5", "0.5"),
            ("1.e5", "100000.0"),
        ] {
            assert_eq!(read(text).as_deref(), Some(expected), "{text:?}");
        }
        for refused in [
            "\u{1c}1", "1\u{1f}", "_1", "1_", "1__0", "1_.5", "1._5", "0x10", ".", "1e", "1e+",
            ".e5", "", "abc", "1,5", "infinit", "1 2", "+-1", "\u{b2}",
        ] {
            assert_eq!(read(refused), None, "{refused:?}");
        }
    }

    #[test]
    fn float_repr_matches_python() {
        let cases = [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (0.5, "0.5"),
            (0.666667, "0.666667"),
            (0.007812, "0.007812"),
            (0.0001, "0.0001"),
            (0.00005, "5e-05"),
            (0.000001, "1e-06"),
            (0.000015, "1.5e-05"),
            (0.000099, "9.9e-05"),
            (1e15, "1000000000000000.0"),
            (1e16, "1e+16"),
            (123.456, "123.456"),
            (0.1 + 0.2, "0.30000000000000004"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
        ];
        for (x, expected) in cases {
            assert_eq!(float_repr(x), expected, "repr({x:?})");
        }
        assert_eq!(float_repr(f64::NAN), "nan");
    }

    /// Each value is given by its exact decimal expansion, which ends in a 5 one place past
    /// the digits `repr` writes.
    #[test]
    fn float_repr_writes_the_even_digit_of_two_texts_equally_near() {
        let lower_is_even = [
            ("1125899906842624.25", "1125899906842624.2"),
            ("10.0009918212890625", "10.000991821289062"),
            ("64.002227783203125", "64.00222778320312"),
            ("724325075565378.25", "724325075565378.2"),
            ("151.749176025390625", "151.74917602539062"),
            ("0.73834991455078125", "0.7383499145507812"),
            ("0.00238513946533203125", "0.0023851394653320312"),
            ("0.0000088214874267578125", "8.821487426757812e-06"),
            ("0.0000000298023223876953125", "2.9802322387695312e-08"),
            ("-1370970376886769.25", "-1370970376886769.2"),
        ];
        let lower_is_odd = [
            ("893978950817812.75", "893978950817812.8"),
            ("1175353826352732.75", "1175353826352732.8"),
            ("2.25878143310546875", "2.2587814331054688"),
            ("64.000701904296875", "64.00070190429688"),
            ("0.00000083446502685546875", "8.344650268554688e-07"),
            ("0.000000178813934326171875", "1.7881393432617188e-07"),
        ];
        for (exact, expected) in lower_is_even.into_iter().chain(lower_is_odd) {
            let x: f64 = exact.parse().unwrap();
            assert_eq!(float_repr(x), expected, "repr({exact})");
        }
        assert_eq!(float_repr(26217.0 / 262144.0 * 100.0), "10.000991821289062");
        assert_eq!(float_repr(83889.0 / 131072.0 * 100.0), "64.00222778320312");
        // Below a power of two the doubles lie twice as close, so the lower text of this
        // tie reads back as another double and the upper one is written, odd digit and all.
        assert_eq!(float_repr(2f64.powi(-24)), "5.960464477539063e-08");
        assert_eq!(float_repr(-(2f64.powi(-24))), "-5.960464477539063e-08");
    }

    #[test]
    fn repr_of_a_text_matches_python() {
        for (text, expected) in [
            ("a b", "'a b'"),
            ("", "''"),
            ("it's", "\"it's\""),
            ("say \"hi\"", "'say \"hi\"'"),
            ("both ' and \"", "'both \\' and \"'"),
            ("tab\there", "'tab\\there'"),
            ("nl\nx", "'nl\\nx'"),
            ("cr\rx", "'cr\\rx'"),
            ("back\\slash", "'back\\\\slash'"),
            ("\0\u{1f}\u{7f}", "'\\x00\\x1f\\x7f'"),
            ("\u{80}\u{a0}\u{ad}", "'\\x80\\xa0\\xad'"),
            ("caf\u{e9}", "'caf\u{e9}'"),
            ("\u{378}", "'\\u0378'"),
            ("\u{2028}", "'\\u2028'"),
            ("\u{1f600}", "'\u{1f600}'"),
            ("\u{e0001}", "'\\U000e0001'"),
            ("\u{fffd}", "'\u{fffd}'"),
            ("\u{feff}", "'\\ufeff'"),
        ] {
            assert_eq!(repr(text), expected, "{text:?}");
        }
    }

    #[test]
    fn ill_formed_utf8_decodes_with_one_replacement_per_maximal_part() {
        let cases: [(&[u8], &str); 20] = [
            (b"abc", "abc"),
            (b"\xff", "\u{fffd}"),
            (b"\xc0\x80", "\u{fffd}\u{fffd}"),
            (b"\xe0\x80\x80", "\u{fffd}\u{fffd}\u{fffd}"),
            (b"\xed\xa0\x80", "\u{fffd}\u{fffd}\u{fffd}"),
            (b"\xf4\x90\x80\x80", "\u{fffd}\u{fffd}\u{fffd}\u{fffd}"),
            (b"\xe2\x82", "\u{fffd}"),
            (b"\xe2\x82A", "\u{fffd}A"),
            (b"\xf0\x90\x80A", "\u{fffd}A"),
            (b"\xf0\x90", "\u{fffd}"),
            (b"A\xc3", "A\u{fffd}"),
            (b"\xc3(", "\u{fffd}("),
            (
                b"\xf8\x88\x80\x80\x80",
                "\u{fffd}\u{fffd}\u{fffd}\u{fffd}\u{fffd}",
            ),
            (b"\xe2(\xa1", "\u{fffd}(\u{fffd}"),
            (b"\xf0(\x8c\xbc", "\u{fffd}(\u{fffd}\u{fffd}"),
            (b"\xf0\x9f\x98\x80", "\u{1f600}"),
            (b"\xf0\x9f\x98", "\u{fffd}"),
            (b"\xc2\xa9\xc2", "\u{a9}\u{fffd}"),
            (b"\x80\x80", "\u{fffd}\u{fffd}"),
            (b"\xf1\x80\x80\xe2\x82\xac", "\u{fffd}\u{20ac}"),
        ];
        for (bytes, expected) in cases {
            assert_eq!(decode(bytes), expected, "{bytes:?}");
        }
    }

    #[test]
    fn path_str_matches_pathlib() {
        let cases = [
            ("./a//b/", "a/b"),
            (
                "/data/work/reference/output.txt",
                "/data/work/reference/output.txt",
            ),
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
                path_str(arg.as_bytes()),
                expected.as_bytes(),
                "Path({arg:?})"
            );
        }
        assert_eq!(path_str(b"./ref\xff.txt/"), b"ref\xff.txt");
    }

    fn sample() -> Json {
        Json::object([
            ("a", Json::int(1)),
            (
                "b",
                Json::List(vec![
                    Json::int(1),
                    Json::Float(2.5),
                    Json::str("x\u{e9}\n"),
                    Json::Null,
                    Json::Bool(true),
                    Json::Bool(false),
                    Json::Object(Vec::new()),
                    Json::List(Vec::new()),
                ]),
            ),
            (
                "c",
                Json::object([
                    ("d", Json::Float(1e16)),
                    ("e", Json::Float(1e-5)),
                    ("f", Json::Float(-0.0)),
                    ("g", Json::Float(f64::INFINITY)),
                    ("h", Json::Float(f64::NAN)),
                    (
                        "i",
                        Json::Int("1000000000000000000000000000000".to_string()),
                    ),
                ]),
            ),
            ("\u{e9}", Json::str("\u{1f600}")),
        ])
    }

    #[test]
    fn dumps_writes_what_json_dumps_writes() {
        let indented = "{\n  \"a\": 1,\n  \"b\": [\n    1,\n    2.5,\n    \"x\\u00e9\\n\",\n    null,\n    true,\n    false,\n    {},\n    []\n  ],\n  \"c\": {\n    \"d\": 1e+16,\n    \"e\": 1e-05,\n    \"f\": -0.0,\n    \"g\": Infinity,\n    \"h\": NaN,\n    \"i\": 1000000000000000000000000000000\n  },\n  \"\\u00e9\": \"\\ud83d\\ude00\"\n}";
        assert_eq!(json_dumps_indented(&sample()), indented);
        let one_line = "{\"a\": 1, \"b\": [1, 2.5, \"x\\u00e9\\n\", null, true, false, {}, []], \"c\": {\"d\": 1e+16, \"e\": 1e-05, \"f\": -0.0, \"g\": Infinity, \"h\": NaN, \"i\": 1000000000000000000000000000000}, \"\\u00e9\": \"\\ud83d\\ude00\"}";
        assert_eq!(json_dumps(&sample()), one_line);
        assert_eq!(json_dumps_indented(&Json::List(Vec::new())), "[]");
        assert_eq!(json_dumps_indented(&Json::str("\u{7f}")), "\"\\u007f\"");
    }

    #[test]
    fn compact_text_sorts_keys_and_unpads_exponents() {
        let compact = "{\"a\":1,\"b\":[1,2.5,\"x\u{e9}\\n\",null,true,false,{},[]],\"c\":{\"d\":1e16,\"e\":1e-5,\"f\":-0.0,\"g\":null,\"h\":null,\"i\":1000000000000000000000000000000},\"\u{e9}\":\"\u{1f600}\"}";
        assert_eq!(json_compact(&sample()), compact);
        assert_eq!(
            json_compact(&Json::str("\u{7f}\u{1}\u{2028}")),
            "\"\u{7f}\\u0001\u{2028}\""
        );
    }

    /// An object's text in compact form, any other text as it is.
    fn canonical(text: &str) -> String {
        match json_loads(text) {
            Ok(value @ Json::Object(_)) => json_compact(&value),
            _ => text.to_string(),
        }
    }

    #[test]
    fn loads_reads_what_json_loads_reads() {
        for (text, expected) in [
            (
                r#"{"b":1,"a":[1,2.0,"x",null,true,false,{"z":1e5,"y":-0.0}],"a":2}"#,
                r#"{"a":2,"b":1}"#,
            ),
            (
                r#"{"a":1e999,"b":-1e999,"c":NaN,"d":Infinity,"e":-Infinity,"f":-0,"g":1E2,"h":1.50,"i":0.1e-6,"j":123456789012345678901234567890}"#,
                r#"{"a":null,"b":null,"c":null,"d":null,"e":null,"f":0,"g":100.0,"h":1.5,"i":1e-7,"j":123456789012345678901234567890}"#,
            ),
            (
                " {\"a\" : \"\\u00e9\\ud83d\\ude00\\n\\/\" } ",
                "{\"a\":\"\u{e9}\u{1f600}\\n/\"}",
            ),
            (r#"{"a":{"b":{"c":[]}}}"#, r#"{"a":{"b":{"c":[]}}}"#),
            (
                "{\"k\":\"caf\u{e9} \u{2028} \u{7f} \\u0000\"}",
                "{\"k\":\"caf\u{e9} \u{2028} \u{7f} \\u0000\"}",
            ),
            (r#"{"a":2.5e-3}"#, r#"{"a":0.0025}"#),
            (r#"{"\u0041":1,"A":2}"#, r#"{"A":2}"#),
            ("{}", "{}"),
        ] {
            assert_eq!(canonical(text), expected, "{text}");
        }
        for unread in [
            r#"{"a":01}"#,
            r#"{"a":1,}"#,
            "[1,]",
            "{\"a\":\"x\ty\"}",
            r#"{"a":1}x"#,
            "\u{feff}{\"a\":1}",
            r#"{"a":.5}"#,
            r#"{"a":1.}"#,
            r#"{"a":+1}"#,
            r#"{"a":"\x"}"#,
            "{'a':1}",
            r#"{"a":tru}"#,
            r#"{"a":1 "b":2}"#,
            r#"{"a" 1}"#,
            r#"{"a":-}"#,
            r#"{"a":-Infinit}"#,
            r#"{"a":2e}"#,
            r#"{"a":2e+}"#,
            "",
            "{",
            r#"{"a":"\u12"}"#,
            r#"{"a":"x"#,
        ] {
            assert_eq!(json_loads(unread), Err(JsonError::Invalid), "{unread}");
        }
        assert_eq!(json_loads("[]"), Ok(Json::List(Vec::new())));
        assert_eq!(json_loads("\"x\""), Ok(Json::str("x")));
        assert_eq!(json_loads("1"), Ok(Json::int(1)));
        assert_eq!(json_loads("1e5"), Ok(Json::Float(100000.0)));
        assert_eq!(json_loads("-0"), Ok(Json::int(0)));
        assert_eq!(
            json_loads(r#"{"a":"\ud800"}"#),
            Ok(Json::object([("a", Json::str("\u{fffd}"))]))
        );
        assert_eq!(
            json_loads(&format!("{{\"a\":{}}}", "9".repeat(4300))).map(|_| ()),
            Ok(())
        );
        assert_eq!(
            json_loads(&format!("{{\"a\":{}}}", "9".repeat(4301))),
            Err(JsonError::Invalid)
        );
    }

    #[test]
    fn a_text_nested_past_the_bound_is_too_deep() {
        let nested = |depth: usize| format!("{}1{}", "{\"a\":".repeat(depth), "}".repeat(depth));
        assert!(json_loads(&nested(JSON_NESTING)).is_ok());
        assert_eq!(
            json_loads(&nested(JSON_NESTING + 1)),
            Err(JsonError::TooDeep)
        );
        let lists = format!(
            "{}{}",
            "[".repeat(JSON_NESTING + 1),
            "]".repeat(JSON_NESTING + 1)
        );
        assert_eq!(json_loads(&lists), Err(JsonError::TooDeep));
    }
}
