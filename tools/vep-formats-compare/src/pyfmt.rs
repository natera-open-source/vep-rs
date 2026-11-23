// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The text forms and text tests of the reference implementation's interpreter.
//!
//! Every function here answers as CPython 3 answers for a `str`: lengths, indexes and
//! slices count code points, never bytes. `sorted()` of `str` values is code-point order,
//! which is the byte order of UTF-8, so `Ord` of `str` and of `String` is that order and
//! no function stands for it.

mod tables;

use std::fmt;
use std::fmt::Write as _;

/// The most digits `int(text)` converts.
pub const INT_MAX_STR_DIGITS: usize = 4300;

/// Why `int(text)` raises `ValueError`; `Display` is the exception's text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntError {
    /// The text is not a whole number; holds `repr(text)`.
    Invalid(String),
    /// The text is a whole number of more digits than the interpreter converts; holds the
    /// count of its digits.
    TooManyDigits(usize),
}

impl fmt::Display for IntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            // The interpreter cuts the repr at 200 code points.
            IntError::Invalid(repr) => write!(
                f,
                "invalid literal for int() with base 10: {}",
                prefix(repr, 200)
            ),
            IntError::TooManyDigits(digits) => write!(
                f,
                "Exceeds the limit ({INT_MAX_STR_DIGITS} digits) for integer string conversion: \
                 value has {digits} digits; use sys.set_int_max_str_digits() to increase the limit"
            ),
        }
    }
}

/// Why `int(value)` raises for a float; `Display` is the exception's text.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FloatToIntError {
    /// A NaN: `ValueError`.
    Nan,
    /// An infinity: `OverflowError`.
    Infinity,
}

impl fmt::Display for FloatToIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            FloatToIntError::Nan => "cannot convert float NaN to integer",
            FloatToIntError::Infinity => "cannot convert float infinity to integer",
        })
    }
}

fn in_set(table: &[u32], c: char) -> bool {
    table.partition_point(|&bound| bound <= c as u32) & 1 == 1
}

/// `c.isprintable()` of one character.
pub fn char_is_printable(c: char) -> bool {
    if c.is_ascii() {
        (' '..='~').contains(&c)
    } else {
        in_set(tables::PRINTABLE, c)
    }
}

/// `repr(text)`.
pub fn repr_str(text: &str) -> String {
    let quote = if text.contains('\'') && !text.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::with_capacity(text.len() + 2);
    out.push(quote);
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if char_is_printable(c) => out.push(c),
            c if (c as u32) < 0x100 => {
                let _ = write!(out, "\\x{:02x}", c as u32);
            }
            c if (c as u32) < 0x10000 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => {
                let _ = write!(out, "\\U{:08x}", c as u32);
            }
        }
    }
    out.push(quote);
    out
}

/// `repr(text)` of a `str | None`: `None` for an absent text.
pub fn repr_opt_str(text: Option<&str>) -> String {
    match text {
        Some(text) => repr_str(text),
        None => "None".to_string(),
    }
}

/// `repr(value)` and `str(value)` of a float.
pub fn repr_f64(value: f64) -> String {
    let mut out = String::with_capacity(24);
    write_repr_f64(&mut out, value);
    out
}

/// Appends `repr(value)` of a float to `out`.
pub fn write_repr_f64(out: &mut String, value: f64) {
    if value.is_nan() {
        out.push_str("nan");
        return;
    }
    if value.is_infinite() {
        out.push_str(if value > 0.0 { "inf" } else { "-inf" });
        return;
    }
    // The fewest digits that read back as the value, as `d[.ddd]e<exponent>`. Where two
    // texts of that many digits read back as the value, the interpreter prints the
    // nearer one and, of two equally near, the one that ends in an even digit, which is
    // the value rounded to that many digits; the formatter's own choice stands only
    // where that rounding reads back as another value.
    let magnitude = value.abs();
    let mut shortest = ShortBuffer::new();
    let _ = write!(shortest, "{magnitude:e}");
    let places = shortest
        .as_str()
        .bytes()
        .take_while(|byte| *byte != b'e')
        .filter(u8::is_ascii_digit)
        .count()
        .saturating_sub(1);
    let mut rounded = ShortBuffer::new();
    let _ = write!(rounded, "{magnitude:.places$e}");
    let chosen = if rounded.as_str().parse() == Ok(magnitude) {
        rounded.as_str()
    } else {
        shortest.as_str()
    };
    let (mantissa, exponent) = chosen.split_once('e').unwrap_or((chosen, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let (lead, tail) = match mantissa.split_once('.') {
        Some((lead, tail)) => (lead, tail.trim_end_matches('0')),
        None => (mantissa, ""),
    };
    if value.is_sign_negative() {
        out.push('-');
    }
    // The digits stand for 0.<lead><tail> times ten to the `point`.
    let point = exponent + 1;
    if point <= -4 || point > 16 {
        out.push_str(lead);
        if !tail.is_empty() {
            out.push('.');
            out.push_str(tail);
        }
        out.push('e');
        out.push(if exponent < 0 { '-' } else { '+' });
        let _ = write!(out, "{:02}", exponent.abs());
    } else if point <= 0 {
        out.push_str("0.");
        for _ in 0..-point {
            out.push('0');
        }
        out.push_str(lead);
        out.push_str(tail);
    } else {
        let point = point as usize;
        let digits = lead.len() + tail.len();
        out.push_str(lead);
        if digits <= point {
            out.push_str(tail);
            for _ in digits..point {
                out.push('0');
            }
            out.push_str(".0");
        } else {
            let split = point - lead.len();
            out.push_str(&tail[..split]);
            out.push('.');
            out.push_str(&tail[split..]);
        }
    }
}

/// A float's shortest exponent form without a heap allocation: at most 17 digits, a
/// point, `e`, a sign and three exponent digits.
struct ShortBuffer {
    bytes: [u8; 32],
    len: usize,
}

impl ShortBuffer {
    fn new() -> ShortBuffer {
        ShortBuffer {
            bytes: [0; 32],
            len: 0,
        }
    }

    fn as_str(&self) -> &str {
        std::str::from_utf8(&self.bytes[..self.len]).unwrap_or("")
    }
}

impl fmt::Write for ShortBuffer {
    fn write_str(&mut self, text: &str) -> fmt::Result {
        let end = self.len + text.len();
        if end > self.bytes.len() {
            return Err(fmt::Error);
        }
        self.bytes[self.len..end].copy_from_slice(text.as_bytes());
        self.len = end;
        Ok(())
    }
}

/// `f"{value:.{digits}f}"`.
pub fn format_fixed(value: f64, digits: usize) -> String {
    if value.is_nan() {
        "nan".to_string()
    } else if value.is_infinite() {
        if value > 0.0 { "inf" } else { "-inf" }.to_string()
    } else {
        format!("{value:.digits$}")
    }
}

/// The white space of C's `isspace`, which `int` and `float` skip around a number.
fn is_c_space(byte: u8) -> bool {
    matches!(byte, b'\t' | b'\n' | 0x0b | 0x0c | b'\r' | b' ')
}

/// The value of a decimal digit of any script.
fn decimal_value(c: char) -> Option<u8> {
    if c.is_ascii_digit() {
        return Some(c as u8 - b'0');
    }
    if c.is_ascii() {
        return None;
    }
    let at = tables::DECIMAL.partition_point(|&bound| bound <= c as u32);
    (at & 1 == 1).then(|| ((c as u32 - tables::DECIMAL[at - 1]) % 10) as u8)
}

/// The text as `int` and `float` read it: every decimal digit as its ASCII digit and
/// every white space beyond ASCII as a space. The text ends in a `?`, which no number
/// holds, at the first character beyond ASCII that is neither.
fn digits_and_spaces_to_ascii(text: &str) -> std::borrow::Cow<'_, [u8]> {
    if text.is_ascii() {
        return std::borrow::Cow::Borrowed(text.as_bytes());
    }
    let mut ascii = Vec::with_capacity(text.len());
    for c in text.chars() {
        if (c as u32) < 127 {
            ascii.push(c as u8);
        } else if char_is_space(c) {
            ascii.push(b' ');
        } else if let Some(digit) = decimal_value(c) {
            ascii.push(b'0' + digit);
        } else {
            ascii.push(b'?');
            break;
        }
    }
    std::borrow::Cow::Owned(ascii)
}

/// `int(text)` where the number fits `i64`; `None` where `int` raises or the number is
/// beyond `i64` (`int_of_big` tells the two apart).
pub fn int_of(text: &str) -> Option<i64> {
    let bytes = text.as_bytes();
    if !bytes.is_empty() && bytes.len() <= 18 && bytes.iter().all(u8::is_ascii_digit) {
        return Some(
            bytes
                .iter()
                .fold(0i64, |n, digit| n * 10 + i64::from(digit - b'0')),
        );
    }
    int_checked(text).ok()?.parse().ok()
}

/// `str(int(text))`: the number's decimal text, of any length the interpreter converts;
/// `None` where `int` raises.
pub fn int_of_big(text: &str) -> Option<String> {
    int_checked(text).ok()
}

/// `str(int(text))`, or why `int` raises.
pub fn int_checked(text: &str) -> Result<String, IntError> {
    let invalid = || IntError::Invalid(repr_str(text));
    let ascii = digits_and_spaces_to_ascii(text);
    let bytes: &[u8] = &ascii;
    let mut at = 0;
    while at < bytes.len() && is_c_space(bytes[at]) {
        at += 1;
    }
    let mut negative = false;
    match bytes.get(at) {
        Some(b'+') => at += 1,
        Some(b'-') => {
            negative = true;
            at += 1;
        }
        _ => {}
    }
    if bytes.get(at) == Some(&b'_') {
        return Err(invalid());
    }
    let first = at;
    let mut digits = 0usize;
    let mut previous = 0u8;
    while at < bytes.len() && (bytes[at].is_ascii_digit() || bytes[at] == b'_') {
        if bytes[at] == b'_' {
            if previous == b'_' {
                return Err(invalid());
            }
        } else {
            digits += 1;
        }
        previous = bytes[at];
        at += 1;
    }
    if previous == b'_' {
        return Err(invalid());
    }
    if at == first {
        return Err(invalid());
    }
    let last = at;
    while at < bytes.len() && is_c_space(bytes[at]) {
        at += 1;
    }
    // The interpreter reads the text as a C string: a NUL ends what it judges, and makes
    // the text invalid only once the number before it has passed the count of digits.
    if at != bytes.len() && bytes[at] != 0 {
        return Err(invalid());
    }
    if digits > INT_MAX_STR_DIGITS {
        return Err(IntError::TooManyDigits(digits));
    }
    if at != bytes.len() {
        return Err(invalid());
    }
    let mut number = String::with_capacity(digits + 1);
    let significant = bytes[first..last]
        .iter()
        .filter(|byte| **byte != b'_')
        .skip_while(|byte| **byte == b'0');
    for digit in significant {
        if number.is_empty() && negative {
            number.push('-');
        }
        number.push(*digit as char);
    }
    if number.is_empty() {
        number.push('0');
    }
    Ok(number)
}

/// `float(text)`; `None` where `float` raises `ValueError`.
pub fn float_of(text: &str) -> Option<f64> {
    let ascii = digits_and_spaces_to_ascii(text);
    let mut bytes: &[u8] = &ascii;
    // An underscore stands between two digits, judged on the text as given.
    let without_underscores: Vec<u8>;
    if bytes.contains(&b'_') {
        let mut kept = Vec::with_capacity(bytes.len());
        let mut previous = 0u8;
        for &byte in bytes {
            if byte == b'_' {
                if !previous.is_ascii_digit() {
                    return None;
                }
            } else {
                if previous == b'_' && !byte.is_ascii_digit() {
                    return None;
                }
                kept.push(byte);
            }
            previous = byte;
        }
        if previous == b'_' {
            return None;
        }
        without_underscores = kept;
        bytes = &without_underscores;
    }
    while let [first, rest @ ..] = bytes {
        if !is_c_space(*first) {
            break;
        }
        bytes = rest;
    }
    while let [rest @ .., last] = bytes {
        if !is_c_space(*last) || rest.is_empty() {
            break;
        }
        bytes = rest;
    }
    if bytes.is_empty() {
        return None;
    }
    let (negative, unsigned) = match bytes {
        [b'-', rest @ ..] => (true, rest),
        [b'+', rest @ ..] => (false, rest),
        _ => (false, bytes),
    };
    let named =
        if unsigned.eq_ignore_ascii_case(b"inf") || unsigned.eq_ignore_ascii_case(b"infinity") {
            Some(f64::INFINITY)
        } else if unsigned.eq_ignore_ascii_case(b"nan") {
            Some(f64::NAN)
        } else {
            None
        };
    if let Some(value) = named {
        return Some(if negative { -value } else { value });
    }
    // digits [. digits] [e [sign] digits], with a digit before or after the point
    let mut at = 0;
    let digits_from = |mut at: usize| {
        while at < unsigned.len() && unsigned[at].is_ascii_digit() {
            at += 1;
        }
        at
    };
    let whole_end = digits_from(at);
    let mut mantissa_digits = whole_end - at;
    at = whole_end;
    if unsigned.get(at) == Some(&b'.') {
        let fraction_end = digits_from(at + 1);
        mantissa_digits += fraction_end - (at + 1);
        at = fraction_end;
    }
    if mantissa_digits == 0 {
        return None;
    }
    if matches!(unsigned.get(at), Some(b'e' | b'E')) {
        let mut exponent = at + 1;
        if matches!(unsigned.get(exponent), Some(b'+' | b'-')) {
            exponent += 1;
        }
        let exponent_end = digits_from(exponent);
        if exponent_end == exponent {
            return None;
        }
        at = exponent_end;
    }
    if at != unsigned.len() {
        return None;
    }
    std::str::from_utf8(bytes).ok()?.parse().ok()
}

/// The text of the `ValueError` `float(text)` raises.
pub fn float_error(text: &str) -> String {
    format!("could not convert string to float: {}", repr_str(text))
}

/// `str(int(value))` of a float: the decimal text of the value cut toward zero.
pub fn int_of_f64(value: f64) -> Result<String, FloatToIntError> {
    if value.is_nan() {
        return Err(FloatToIntError::Nan);
    }
    if value.is_infinite() {
        return Err(FloatToIntError::Infinity);
    }
    let whole = value.trunc();
    if whole == 0.0 {
        return Ok("0".to_string());
    }
    Ok(format!("{whole:.0}"))
}

/// `float(n)` of a whole number given as its decimal text; `None` where `float` raises
/// `OverflowError`.
pub fn f64_of_int(decimal: &str) -> Option<f64> {
    let value: f64 = decimal.parse().ok()?;
    value.is_finite().then_some(value)
}

/// `c.isdecimal()` of one character.
pub fn char_is_decimal(c: char) -> bool {
    c.is_ascii_digit() || (!c.is_ascii() && in_set(tables::DECIMAL, c))
}

/// `c.isdigit()` of one character.
pub fn char_is_digit(c: char) -> bool {
    c.is_ascii_digit() || (!c.is_ascii() && in_set(tables::DIGIT, c))
}

/// `c.isnumeric()` of one character.
pub fn char_is_numeric(c: char) -> bool {
    c.is_ascii_digit() || (!c.is_ascii() && in_set(tables::NUMERIC, c))
}

/// `c.isalpha()` of one character.
pub fn char_is_alpha(c: char) -> bool {
    c.is_ascii_alphabetic() || (!c.is_ascii() && in_set(tables::ALPHA, c))
}

/// `c.isspace()` of one character.
pub fn char_is_space(c: char) -> bool {
    matches!(
        c,
        '\t'..='\r'
            | '\x1c'..='\x1f'
            | ' '
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

fn all_chars(text: &str, test: impl Fn(char) -> bool) -> bool {
    !text.is_empty() && text.chars().all(test)
}

/// `text.isdigit()`.
pub fn is_digit(text: &str) -> bool {
    all_chars(text, char_is_digit)
}

/// `text.isdecimal()`.
pub fn is_decimal(text: &str) -> bool {
    all_chars(text, char_is_decimal)
}

/// `text.isnumeric()`.
pub fn is_numeric(text: &str) -> bool {
    all_chars(text, char_is_numeric)
}

/// `text.isalpha()`.
pub fn is_alpha(text: &str) -> bool {
    all_chars(text, char_is_alpha)
}

/// `text.isalnum()`.
pub fn is_alnum(text: &str) -> bool {
    all_chars(text, |c| {
        char_is_alpha(c) || char_is_decimal(c) || char_is_digit(c) || char_is_numeric(c)
    })
}

/// `text.isspace()`.
pub fn is_space(text: &str) -> bool {
    all_chars(text, char_is_space)
}

/// `text.isprintable()`: true of the empty text.
pub fn is_printable(text: &str) -> bool {
    text.chars().all(char_is_printable)
}

fn char_is_lower(c: char) -> bool {
    c.is_ascii_lowercase() || (!c.is_ascii() && in_set(tables::LOWER, c))
}

fn char_is_upper(c: char) -> bool {
    c.is_ascii_uppercase() || (!c.is_ascii() && in_set(tables::UPPER, c))
}

fn char_is_title(c: char) -> bool {
    !c.is_ascii() && in_set(tables::TITLE, c)
}

/// `text.islower()`.
pub fn is_lower(text: &str) -> bool {
    let mut cased = false;
    for c in text.chars() {
        if char_is_upper(c) || char_is_title(c) {
            return false;
        }
        cased = cased || char_is_lower(c);
    }
    cased
}

/// `text.isupper()`.
pub fn is_upper(text: &str) -> bool {
    let mut cased = false;
    for c in text.chars() {
        if char_is_lower(c) || char_is_title(c) {
            return false;
        }
        cased = cased || char_is_upper(c);
    }
    cased
}

/// `text.isascii()`: true of the empty text.
pub fn is_ascii(text: &str) -> bool {
    text.is_ascii()
}

/// `text.strip()`.
pub fn strip(text: &str) -> &str {
    text.trim_matches(char_is_space)
}

/// `text.lstrip()`.
pub fn lstrip(text: &str) -> &str {
    text.trim_start_matches(char_is_space)
}

/// `text.rstrip()`.
pub fn rstrip(text: &str) -> &str {
    text.trim_end_matches(char_is_space)
}

/// `text.strip(chars)`.
pub fn strip_chars<'a>(text: &'a str, chars: &str) -> &'a str {
    text.trim_matches(|c| chars.contains(c))
}

/// `text.lstrip(chars)`.
pub fn lstrip_chars<'a>(text: &'a str, chars: &str) -> &'a str {
    text.trim_start_matches(|c| chars.contains(c))
}

/// `text.rstrip(chars)`.
pub fn rstrip_chars<'a>(text: &'a str, chars: &str) -> &'a str {
    text.trim_end_matches(|c| chars.contains(c))
}

/// `text.split()`.
pub fn split_whitespace(text: &str) -> Vec<&str> {
    text.split(char_is_space)
        .filter(|word| !word.is_empty())
        .collect()
}

/// `text.split(sep)` for a non-empty `sep`.
pub fn split<'a>(text: &'a str, sep: &str) -> Vec<&'a str> {
    text.split(sep).collect()
}

/// `text.split(sep, maxsplit)` for a non-empty `sep`.
pub fn split_n<'a>(text: &'a str, sep: &str, maxsplit: usize) -> Vec<&'a str> {
    text.splitn(maxsplit.saturating_add(1), sep).collect()
}

/// `text.rsplit(sep, maxsplit)` for a non-empty `sep`.
pub fn rsplit_n<'a>(text: &'a str, sep: &str, maxsplit: usize) -> Vec<&'a str> {
    let mut parts: Vec<&str> = text.rsplitn(maxsplit.saturating_add(1), sep).collect();
    parts.reverse();
    parts
}

/// `text.splitlines()`.
pub fn split_lines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut start = 0;
    let mut chars = text.char_indices().peekable();
    while let Some((at, c)) = chars.next() {
        if !matches!(
            c,
            '\n' | '\r'
                | '\x0b'
                | '\x0c'
                | '\x1c'
                | '\x1d'
                | '\x1e'
                | '\u{85}'
                | '\u{2028}'
                | '\u{2029}'
        ) {
            continue;
        }
        lines.push(&text[start..at]);
        start = at + c.len_utf8();
        if c == '\r' && chars.next_if(|(_, next)| *next == '\n').is_some() {
            start += 1;
        }
    }
    if start < text.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// `text.partition(sep)` for a non-empty `sep`.
pub fn partition<'a>(text: &'a str, sep: &str) -> (&'a str, &'a str, &'a str) {
    match text.find(sep) {
        Some(at) => (
            &text[..at],
            &text[at..at + sep.len()],
            &text[at + sep.len()..],
        ),
        None => (text, "", ""),
    }
}

/// `text.rpartition(sep)` for a non-empty `sep`.
pub fn rpartition<'a>(text: &'a str, sep: &str) -> (&'a str, &'a str, &'a str) {
    match text.rfind(sep) {
        Some(at) => (
            &text[..at],
            &text[at..at + sep.len()],
            &text[at + sep.len()..],
        ),
        None => ("", "", text),
    }
}

/// Appends the image of `c` under a case map and its list of longer images.
fn push_mapped(out: &mut String, c: char, map: &[i32], special: &[u32]) {
    let code = c as u32 as i32;
    let runs = map.len() / 4;
    let (mut low, mut high) = (0, runs);
    while low < high {
        let middle = (low + high) / 2;
        if map[middle * 4] <= code {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    if low > 0 {
        let run = &map[(low - 1) * 4..low * 4];
        if code <= run[1] && (code - run[0]) % run[2] == 0 {
            out.push(char::from_u32((code + run[3]) as u32).unwrap_or(c));
            return;
        }
    }
    let (mut low, mut high) = (0, special.len() / 4);
    while low < high {
        let middle = (low + high) / 2;
        if special[middle * 4] < c as u32 {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    if low < special.len() / 4 && special[low * 4] == c as u32 {
        for &image in &special[low * 4 + 1..low * 4 + 4] {
            if image != 0 {
                out.push(char::from_u32(image).unwrap_or(c));
            }
        }
        return;
    }
    out.push(c);
}

/// Whether the capital sigma at byte `at` of `text` ends a word: a cased character
/// stands before it and none after it, the characters `tables::CASE_IGNORABLE` lists
/// skipped on both sides.
fn sigma_is_final(text: &str, at: usize) -> bool {
    let counts = |c: &char| !in_set(tables::CASE_IGNORABLE, *c);
    let Some(before) = text[..at].chars().rev().find(counts) else {
        return false;
    };
    if !in_set(tables::CASED, before) {
        return false;
    }
    match text[at + 'Σ'.len_utf8()..].chars().find(counts) {
        Some(after) => !in_set(tables::CASED, after),
        None => true,
    }
}

/// `text.lower()`.
pub fn lower(text: &str) -> String {
    if text.is_ascii() {
        return text.to_ascii_lowercase();
    }
    let mut out = String::with_capacity(text.len());
    for (at, c) in text.char_indices() {
        if c.is_ascii() {
            out.push(c.to_ascii_lowercase());
        } else if c == 'Σ' {
            out.push(if sigma_is_final(text, at) { 'ς' } else { 'σ' });
        } else {
            push_mapped(&mut out, c, tables::LOWER_MAP, tables::LOWER_SPECIAL);
        }
    }
    out
}

/// `text.upper()`.
pub fn upper(text: &str) -> String {
    if text.is_ascii() {
        return text.to_ascii_uppercase();
    }
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        if c.is_ascii() {
            out.push(c.to_ascii_uppercase());
        } else {
            push_mapped(&mut out, c, tables::UPPER_MAP, tables::UPPER_SPECIAL);
        }
    }
    out
}

/// `text.removeprefix(prefix)`.
pub fn removeprefix<'a>(text: &'a str, prefix: &str) -> &'a str {
    text.strip_prefix(prefix).unwrap_or(text)
}

/// `text.removesuffix(suffix)`.
pub fn removesuffix<'a>(text: &'a str, suffix: &str) -> &'a str {
    text.strip_suffix(suffix).unwrap_or(text)
}

/// `len(text)`: the count of code points.
pub fn len(text: &str) -> usize {
    if text.is_ascii() {
        text.len()
    } else {
        text.chars().count()
    }
}

/// The byte offset of the code point at index `chars`, the text's length past its end.
fn byte_offset(text: &str, chars: usize) -> usize {
    if chars >= text.len() {
        return text.len();
    }
    if text.is_ascii() {
        return chars;
    }
    text.char_indices()
        .nth(chars)
        .map_or(text.len(), |(at, _)| at)
}

/// `text[:chars]` for a count of code points.
pub fn prefix(text: &str, chars: usize) -> &str {
    &text[..byte_offset(text, chars)]
}

/// `text[start:stop]`: either bound may be absent or negative, as in a slice.
pub fn slice(text: &str, start: Option<i64>, stop: Option<i64>) -> &str {
    let length = len(text) as i64;
    let clamp = |bound: Option<i64>, absent: i64| match bound {
        None => absent,
        Some(bound) if bound < 0 => (bound + length).max(0),
        Some(bound) => bound.min(length),
    };
    let (start, stop) = (clamp(start, 0), clamp(stop, length));
    if start >= stop {
        return "";
    }
    let from = byte_offset(text, start as usize);
    let to = from + byte_offset(&text[from..], (stop - start) as usize);
    &text[from..to]
}

/// The end of the bytes the UTF-8 decoder refuses at `start`, and its reason.
fn utf8_refusal(bytes: &[u8], start: usize) -> (usize, &'static str) {
    const START: &str = "invalid start byte";
    const CONTINUATION: &str = "invalid continuation byte";
    const END: &str = "unexpected end of data";
    let rest = &bytes[start..];
    let lead = rest[0];
    let continues = |byte: u8| (0x80..0xc0).contains(&byte);
    if !(0xc2..0xf5).contains(&lead) {
        return (start + 1, START);
    }
    let Some(&second) = rest.get(1) else {
        return (bytes.len(), END);
    };
    // A second byte that would spell a shorter form, a surrogate or a code point past
    // the last one is refused as a continuation.
    let second_fits = continues(second)
        && match lead {
            0xe0 => second >= 0xa0,
            0xed => second < 0xa0,
            0xf0 => second >= 0x90,
            0xf4 => second < 0x90,
            _ => true,
        };
    if !second_fits {
        return (start + 1, CONTINUATION);
    }
    // A two-byte form with a fitting second byte is well formed, so the first byte asks for more.
    let Some(&third) = rest.get(2) else {
        return (bytes.len(), END);
    };
    if !continues(third) {
        return (start + 2, CONTINUATION);
    }
    if rest.len() < 4 {
        return (bytes.len(), END);
    }
    (start + 3, CONTINUATION)
}

/// `bytes.decode("utf-8")`; the error is the text of the `UnicodeDecodeError`.
pub fn decode_utf8(bytes: &[u8]) -> Result<&str, String> {
    let error = match std::str::from_utf8(bytes) {
        Ok(text) => return Ok(text),
        Err(error) => error,
    };
    let start = error.valid_up_to();
    let (end, reason) = utf8_refusal(bytes, start);
    Err(if end == start + 1 {
        format!(
            "'utf-8' codec can't decode byte 0x{:02x} in position {start}: {reason}",
            bytes[start]
        )
    } else {
        format!(
            "'utf-8' codec can't decode bytes in position {start}-{}: {reason}",
            end - 1
        )
    })
}

/// `bytes.decode("utf-8", "replace")`.
pub fn decode_utf8_replace(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len());
    let mut at = 0;
    loop {
        match std::str::from_utf8(&bytes[at..]) {
            Ok(text) => {
                out.push_str(text);
                return out;
            }
            Err(error) => {
                let start = at + error.valid_up_to();
                // The bytes before the refused ones are UTF-8, as the error says.
                out.push_str(std::str::from_utf8(&bytes[at..start]).unwrap_or(""));
                out.push('\u{fffd}');
                at = utf8_refusal(bytes, start).0;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_repr_picks_the_quote_the_text_does_not_hold_and_escapes_what_does_not_print() {
        assert_eq!(repr_str("MODIFIER"), "'MODIFIER'");
        assert_eq!(repr_str("it's"), "\"it's\"");
        assert_eq!(repr_str("a \"b\" c's"), "'a \"b\" c\\'s'");
        assert_eq!(repr_str("a\tb\nc\rd\\e"), "'a\\tb\\nc\\rd\\\\e'");
        assert_eq!(
            repr_str("\0\x1f\x7f\u{a0}\u{ad}"),
            "'\\x00\\x1f\\x7f\\xa0\\xad'"
        );
        assert_eq!(repr_str("é\u{2028}\u{e000}"), "'é\\u2028\\ue000'");
        assert_eq!(repr_str("\u{1f600}\u{e0001}"), "'\u{1f600}\\U000e0001'");
        assert_eq!(repr_opt_str(None), "None");
        assert_eq!(repr_opt_str(Some("")), "''");
    }

    #[test]
    fn a_float_prints_its_shortest_digits_with_an_exponent_outside_sixteen_places() {
        for (value, text) in [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (0.5, "0.5"),
            (0.1, "0.1"),
            (1e-5, "1e-05"),
            (0.0001, "0.0001"),
            (0.00012345, "0.00012345"),
            (1e15, "1000000000000000.0"),
            (1e16, "1e+16"),
            (1.5e16, "1.5e+16"),
            (123456789.125, "123456789.125"),
            (1e22, "1e+22"),
            (1.7976931348623157e308, "1.7976931348623157e+308"),
            (5e-324, "5e-324"),
            (2.2250738585072014e-308, "2.2250738585072014e-308"),
            (-1234.5, "-1234.5"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            (f64::NAN, "nan"),
        ] {
            assert_eq!(repr_f64(value), text);
        }
    }

    #[test]
    fn a_fixed_form_rounds_the_exact_binary_value_half_to_even() {
        assert_eq!(format_fixed(0.5, 0), "0");
        assert_eq!(format_fixed(1.5, 0), "2");
        assert_eq!(format_fixed(2.5, 0), "2");
        assert_eq!(format_fixed(0.125, 2), "0.12");
        assert_eq!(format_fixed(2.675, 2), "2.67");
        assert_eq!(format_fixed(-0.001, 2), "-0.00");
        assert_eq!(format_fixed(1e21, 1), "1000000000000000000000.0");
        assert_eq!(format_fixed(f64::NAN, 3), "nan");
        assert_eq!(format_fixed(f64::NEG_INFINITY, 3), "-inf");
    }

    #[test]
    fn a_whole_number_is_read_through_white_space_a_sign_underscores_and_any_digits() {
        assert_eq!(int_of("42"), Some(42));
        assert_eq!(int_of(" \t+1_000\n"), Some(1000));
        assert_eq!(int_of("-0"), Some(0));
        assert_eq!(int_of("\u{a0}١٢٣\u{3000}"), Some(123));
        assert_eq!(int_of("9223372036854775807"), Some(i64::MAX));
        assert_eq!(int_of("-9223372036854775808"), Some(i64::MIN));
        assert_eq!(int_of("9223372036854775808"), None);
        assert_eq!(
            int_of_big("-0009223372036854775808000").as_deref(),
            Some("-9223372036854775808000")
        );
        for refused in [
            "", " ", "+", "_1", "1_", "1__0", "1.0", "1e3", "0x10", "+ 1", "²", "\x1c1", "1 2",
        ] {
            assert_eq!(int_of_big(refused), None, "{refused:?}");
        }
        assert_eq!(
            int_checked("x'y").unwrap_err().to_string(),
            "invalid literal for int() with base 10: \"x'y\""
        );
    }

    #[test]
    fn a_whole_number_of_more_digits_than_the_interpreter_converts_is_told_apart() {
        let long = "9".repeat(4301);
        assert_eq!(int_checked(&long), Err(IntError::TooManyDigits(4301)));
        assert_eq!(
            int_checked(&long).unwrap_err().to_string(),
            "Exceeds the limit (4300 digits) for integer string conversion: value has 4301 digits; \
             use sys.set_int_max_str_digits() to increase the limit"
        );
        assert_eq!(int_checked(&"9".repeat(4300)).map(|n| n.len()), Ok(4300));
        assert_eq!(int_of_big(&long), None);
    }

    #[test]
    fn a_float_is_read_with_its_names_underscores_and_any_digits() {
        assert_eq!(float_of("1.5"), Some(1.5));
        assert_eq!(float_of(" 1_0.5e1_0 "), Some(10.5e10));
        assert_eq!(float_of(".5"), Some(0.5));
        assert_eq!(float_of("5."), Some(5.0));
        assert_eq!(float_of("-Infinity"), Some(f64::NEG_INFINITY));
        assert_eq!(float_of("INF"), Some(f64::INFINITY));
        assert!(float_of("+nan").is_some_and(f64::is_nan));
        assert_eq!(float_of("1e400"), Some(f64::INFINITY));
        assert_eq!(float_of("１２.５"), Some(12.5));
        for refused in [
            "", " ", ".", "e5", "1e", "1_", "_1", "1__0", "1_.5", "0x10", "infinit", "nan_", "1 2",
            "\x1c1",
        ] {
            assert_eq!(float_of(refused), None, "{refused:?}");
        }
        assert_eq!(
            float_error("a'b"),
            "could not convert string to float: \"a'b\""
        );
    }

    #[test]
    fn a_float_cut_to_a_whole_number_keeps_every_digit() {
        assert_eq!(int_of_f64(1.9).as_deref(), Ok("1"));
        assert_eq!(int_of_f64(-1.9).as_deref(), Ok("-1"));
        assert_eq!(int_of_f64(-0.5).as_deref(), Ok("0"));
        assert_eq!(int_of_f64(1e23).as_deref(), Ok("99999999999999991611392"));
        assert_eq!(int_of_f64(f64::NAN), Err(FloatToIntError::Nan));
        assert_eq!(int_of_f64(f64::INFINITY), Err(FloatToIntError::Infinity));
        let past_53_bits = (1_u64 << 53) + 1;
        assert_eq!(
            f64_of_int(&past_53_bits.to_string()),
            Some(9007199254740992.0)
        );
        assert_eq!(f64_of_int(&"9".repeat(400)), None);
    }

    #[test]
    fn the_text_tests_read_the_characters_as_the_interpreter_classes_them() {
        assert!(is_digit("123") && is_digit("²") && !is_digit("") && !is_digit("1.0"));
        assert!(is_decimal("١٢٣") && !is_decimal("²") && !is_decimal(""));
        assert!(is_numeric("½") && is_numeric("四") && !is_numeric("a"));
        assert!(is_alpha("aé四") && !is_alpha("a1") && !is_alpha(""));
        assert!(is_alnum("a1½") && !is_alnum("a-1"));
        assert!(is_space(" \t\x1c\u{a0}\u{3000}") && !is_space("") && !is_space("\u{200b}"));
        assert!(is_printable("") && is_printable("a b") && !is_printable("a\u{a0}"));
        assert!(is_lower("a1") && !is_lower("1") && !is_lower("aB") && !is_lower("aǅ"));
        assert!(is_upper("A1") && !is_upper("") && !is_upper("Ab") && !is_upper("Aǅ"));
        assert!(is_ascii("") && is_ascii("a") && !is_ascii("é"));
    }

    #[test]
    fn stripping_and_splitting_read_white_space_and_line_ends_as_the_interpreter_does() {
        assert_eq!(strip("\x1f a b\u{a0}\n"), "a b");
        assert_eq!(lstrip("  a  "), "a  ");
        assert_eq!(rstrip("  a  "), "  a");
        assert_eq!(strip_chars("--a-b--", "-"), "a-b");
        assert_eq!(lstrip_chars("xyax", "xy"), "ax");
        assert_eq!(rstrip_chars("1.500", "0"), "1.5");
        assert_eq!(split_whitespace("  a\tb\u{2003}c \n"), ["a", "b", "c"]);
        assert!(split_whitespace(" \n").is_empty());
        assert_eq!(split("a,,b", ","), ["a", "", "b"]);
        assert_eq!(split("", ","), [""]);
        assert_eq!(split_n("a: b: c", ": ", 1), ["a", "b: c"]);
        assert_eq!(rsplit_n("a|b|c", "|", 1), ["a|b", "c"]);
        assert_eq!(rsplit_n("aaa", "aa", 1), ["a", ""]);
        assert_eq!(
            split_lines("a\nb\r\nc\rd\x0be\x0cf\x1cg\u{85}h\u{2028}i\n"),
            ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
        );
        assert_eq!(split_lines("\n\n"), ["", ""]);
        assert!(split_lines("").is_empty());
        assert_eq!(partition("a=b=c", "="), ("a", "=", "b=c"));
        assert_eq!(partition("abc", "="), ("abc", "", ""));
        assert_eq!(rpartition("a=b=c", "="), ("a=b", "=", "c"));
        assert_eq!(rpartition("abc", "="), ("", "", "abc"));
        assert_eq!(removeprefix("chr21", "chr"), "21");
        assert_eq!(removesuffix("a.gz", ".gz"), "a");
        assert_eq!(removesuffix("a.gz", ".bz"), "a.gz");
    }

    #[test]
    fn case_changes_follow_the_interpreters_maps() {
        assert_eq!(lower("ACGT"), "acgt");
        assert_eq!(upper("acgt"), "ACGT");
        assert_eq!(upper("straße"), "STRASSE");
        assert_eq!(lower("İ"), "i\u{307}");
        assert_eq!(lower("ΟΔΟΣ ΣΟΦΟΣ"), "οδος σοφος");
        assert_eq!(lower("Σ"), "σ");
        assert_eq!(upper("ǆ"), "Ǆ");
        assert_eq!(upper("ﬁ"), "FI");
    }

    #[test]
    fn lengths_and_slices_count_code_points() {
        assert_eq!(len("aé四"), 3);
        assert_eq!(prefix("aé四b", 3), "aé四");
        assert_eq!(prefix("ab", 5), "ab");
        assert_eq!(slice("aé四b", Some(1), Some(3)), "é四");
        assert_eq!(slice("aé四b", Some(-2), None), "四b");
        assert_eq!(slice("aé四b", None, Some(-1)), "aé四");
        assert_eq!(slice("aé四b", Some(3), Some(1)), "");
        assert_eq!(slice("aé四b", Some(-9), Some(9)), "aé四b");
    }

    #[test]
    fn bytes_that_are_not_utf8_are_refused_or_replaced_as_the_interpreter_does() {
        assert_eq!(decode_utf8(b"abc"), Ok("abc"));
        assert_eq!(
            decode_utf8(b"a\xff").unwrap_err(),
            "'utf-8' codec can't decode byte 0xff in position 1: invalid start byte"
        );
        assert_eq!(
            decode_utf8(b"\xe2\x82").unwrap_err(),
            "'utf-8' codec can't decode bytes in position 0-1: unexpected end of data"
        );
        assert_eq!(
            decode_utf8(b"\xf0\x9f\x98A").unwrap_err(),
            "'utf-8' codec can't decode bytes in position 0-2: invalid continuation byte"
        );
        assert_eq!(decode_utf8_replace(b"a\xffb\xe2\x82"), "a\u{fffd}b\u{fffd}");
        assert_eq!(
            decode_utf8_replace(b"\xed\xa0\x80"),
            "\u{fffd}\u{fffd}\u{fffd}"
        );
    }
}
