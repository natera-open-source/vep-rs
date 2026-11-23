// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Perl's numbers as the conversion makes them: what a string is as a number,
//! what `+` and `-` return, and the text JSON::XS writes for the result.
//!
//! A Perl number is an integer (signed or unsigned 64-bit) or a double. Arithmetic on
//! two values that are exact integers stays an integer unless it leaves the 64-bit
//! range; anything else is a double. JSON::XS writes an integer in decimal and a
//! double with C's `%.15g`.

/// A Perl numeric value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Num {
    /// An integer in Perl's range: `i64::MIN ..= u64::MAX`.
    Int(i128),
    /// A double: one stored as a double, or the result of arithmetic. Arithmetic takes
    /// it for an integer when it is whole and below 2^53.
    Float(f64),
    /// A double read from a string that is not an integer's text: one with a decimal
    /// point, or with text after the number, or no number at all. Arithmetic never
    /// takes it for an integer, whole or not.
    Decimal(f64),
}

const INT_MIN: i128 = i64::MIN as i128;
const INT_MAX: i128 = u64::MAX as i128;
/// Doubles hold every integer below 2^53 exactly; Perl treats a double as an exact
/// integer only under that bound.
const EXACT_DOUBLE: f64 = 9_007_199_254_740_992.0;

impl Num {
    pub fn from_i64(value: i64) -> Num {
        Num::Int(i128::from(value))
    }

    /// A double as it enters arithmetic.
    pub fn from_f64(value: f64) -> Num {
        Num::Float(value)
    }

    pub fn as_f64(self) -> f64 {
        match self {
            Num::Int(value) => value as f64,
            Num::Float(value) | Num::Decimal(value) => value,
        }
    }

    /// The value as the exact integer Perl's integer-preserving arithmetic sees: an
    /// integer as it is, a double when it is whole and below 2^53.
    fn exact_int(self) -> Option<i128> {
        match self {
            Num::Int(value) => Some(value),
            Num::Float(value) => {
                (value.fract() == 0.0 && value.abs() < EXACT_DOUBLE).then_some(value as i128)
            }
            Num::Decimal(_) => None,
        }
    }

    fn combine(self, other: Num, int: fn(i128, i128) -> i128, float: fn(f64, f64) -> f64) -> Num {
        if let (Some(left), Some(right)) = (self.exact_int(), other.exact_int()) {
            let result = int(left, right);
            if (INT_MIN..=INT_MAX).contains(&result) {
                return Num::Int(result);
            }
        }
        Num::Float(float(self.as_f64(), other.as_f64()))
    }

    /// Perl's `+`.
    pub fn plus(self, other: Num) -> Num {
        self.combine(
            other,
            |left, right| left + right,
            |left, right| left + right,
        )
    }

    /// Perl's `-`.
    pub fn minus(self, other: Num) -> Num {
        self.combine(
            other,
            |left, right| left - right,
            |left, right| left - right,
        )
    }

    /// Perl's `<=>` for values that are not NaN; `None` when either is.
    pub fn compare(self, other: Num) -> Option<std::cmp::Ordering> {
        match (self, other) {
            (Num::Int(left), Num::Int(right)) => Some(left.cmp(&right)),
            _ => self.as_f64().partial_cmp(&other.as_f64()),
        }
    }

    /// Perl's own string for the number (`"$x"`).
    pub fn to_perl_string(self) -> Vec<u8> {
        match self {
            Num::Int(value) => {
                let mut out = Vec::new();
                write_int(value, &mut out);
                out
            }
            Num::Float(value) | Num::Decimal(value) => double_to_perl_string(value),
        }
    }

    /// The JSON text JSON::XS writes: decimal for an integer, `%.15g` for a double.
    pub fn write_json(self, out: &mut Vec<u8>) {
        match self {
            Num::Int(value) => write_int(value, out),
            Num::Float(value) | Num::Decimal(value) => write_g15(value, out),
        }
    }
}

/// Perl's own string for a double (`"$x"`): `%.15g`, with `Inf`, `-Inf` and `NaN` for
/// the values C prints in lower case.
pub fn double_to_perl_string(value: f64) -> Vec<u8> {
    if value.is_nan() {
        return b"NaN".to_vec();
    }
    if value.is_infinite() {
        return if value < 0.0 {
            b"-Inf".to_vec()
        } else {
            b"Inf".to_vec()
        };
    }
    let mut out = Vec::new();
    write_g15(value, &mut out);
    out
}

pub fn write_int(value: i128, out: &mut Vec<u8>) {
    let mut digits = [0u8; 40];
    let mut at = digits.len();
    let mut rest = value.unsigned_abs();
    loop {
        at -= 1;
        digits[at] = b'0' + (rest % 10) as u8;
        rest /= 10;
        if rest == 0 {
            break;
        }
    }
    if value < 0 {
        out.push(b'-');
    }
    out.extend_from_slice(&digits[at..]);
}

/// C's `%.15g` of a double, as glibc prints it.
pub fn write_g15(value: f64, out: &mut Vec<u8>) {
    const PRECISION: i32 = 15;
    if value.is_nan() {
        out.extend_from_slice(if value.is_sign_negative() {
            b"-nan"
        } else {
            b"nan"
        });
        return;
    }
    if value.is_sign_negative() {
        out.push(b'-');
    }
    let magnitude = value.abs();
    if magnitude.is_infinite() {
        out.extend_from_slice(b"inf");
        return;
    }
    if magnitude == 0.0 {
        out.push(b'0');
        return;
    }
    // Fifteen significant digits, correctly rounded, as `d.dddddddddddddde<exponent>`.
    let scientific = format!("{magnitude:.14e}");
    let (mantissa, exponent) = scientific.split_once('e').expect("an exponent");
    let exponent: i32 = exponent.parse().expect("a decimal exponent");
    let mut digits: Vec<u8> = mantissa.bytes().filter(u8::is_ascii_digit).collect();
    while digits.len() > 1 && digits.last() == Some(&b'0') {
        digits.pop();
    }
    if !(-4..PRECISION).contains(&exponent) {
        out.push(digits[0]);
        if digits.len() > 1 {
            out.push(b'.');
            out.extend_from_slice(&digits[1..]);
        }
        out.push(b'e');
        out.push(if exponent < 0 { b'-' } else { b'+' });
        let exponent = exponent.unsigned_abs();
        if exponent < 10 {
            out.push(b'0');
        }
        out.extend_from_slice(exponent.to_string().as_bytes());
    } else if exponent < 0 {
        out.extend_from_slice(b"0.");
        out.resize(out.len() + (-exponent - 1) as usize, b'0');
        out.extend_from_slice(&digits);
    } else {
        let whole = exponent as usize + 1;
        if digits.len() <= whole {
            out.extend_from_slice(&digits);
            out.resize(out.len() + whole - digits.len(), b'0');
        } else {
            out.extend_from_slice(&digits[..whole]);
            out.push(b'.');
            out.extend_from_slice(&digits[whole..]);
        }
    }
}

fn is_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

/// A string as a number, the way Perl reads one in numeric context: white space, a
/// sign, then decimal digits with an optional fraction and exponent, or an infinity or
/// NaN word; the longest such prefix counts and a string with none is 0.
///
/// What arithmetic then makes of it follows Perl's own cases. Digits alone that fit 64
/// bits are an integer. A number written with an exponent is an integer when its value
/// is whole and in the 64-bit range. A number with a decimal point and no exponent, a
/// number followed by other text, digits past 64 bits and a string with no number are
/// doubles that arithmetic does not take for integers.
pub fn numify(text: &[u8]) -> Num {
    let mut at = 0;
    while at < text.len() && is_space(text[at]) {
        at += 1;
    }
    let start = at;
    let negative = match text.get(at) {
        Some(b'-') => {
            at += 1;
            true
        }
        Some(b'+') => {
            at += 1;
            false
        }
        _ => false,
    };
    let digits_start = at;
    while at < text.len() && text[at].is_ascii_digit() {
        at += 1;
    }
    let whole_digits = at - digits_start;
    let mut end = at;
    let mut has_point = false;
    if text.get(at) == Some(&b'.') {
        let mut fraction = at + 1;
        while fraction < text.len() && text[fraction].is_ascii_digit() {
            fraction += 1;
        }
        // `1.` and `.5` are numbers; a point with digits on neither side is not.
        if whole_digits > 0 || fraction > at + 1 {
            has_point = true;
            end = fraction;
        }
    }
    if end == digits_start {
        return special(&text[digits_start..], negative);
    }
    let mut has_exponent = false;
    if matches!(text.get(end), Some(b'e' | b'E')) {
        let mut exponent = end + 1;
        if matches!(text.get(exponent), Some(b'+' | b'-')) {
            exponent += 1;
        }
        let exponent_digits = exponent;
        while exponent < text.len() && text[exponent].is_ascii_digit() {
            exponent += 1;
        }
        if exponent > exponent_digits {
            has_exponent = true;
            end = exponent;
        }
    }
    // White space may follow a number; anything else after it makes the string one
    // Perl reads the number out of without taking the string for a number.
    let clean = text[end..].iter().all(|byte| is_space(*byte));
    if clean && !has_point && !has_exponent {
        let mut value: i128 = 0;
        let mut fits = true;
        for digit in &text[digits_start..end] {
            value = value * 10 + i128::from(digit - b'0');
            if value > INT_MAX {
                fits = false;
                break;
            }
        }
        let value = if negative { -value } else { value };
        if fits && value >= INT_MIN {
            return Num::Int(value);
        }
    }
    let prefix = std::str::from_utf8(&text[start..end]).expect("ASCII number text");
    // Rust reads the same decimal grammar and rounds to the nearest double.
    let value = prefix.parse::<f64>().unwrap_or(0.0);
    if clean && has_exponent && value.fract() == 0.0 {
        let whole = value as i128;
        if (INT_MIN..=INT_MAX).contains(&whole) && whole as f64 == value {
            return Num::Int(whole);
        }
    }
    Num::Decimal(value)
}

/// The infinity and NaN words Perl reads, case-insensitively, after the sign.
fn special(rest: &[u8], negative: bool) -> Num {
    let lower: Vec<u8> = rest.iter().take(8).map(u8::to_ascii_lowercase).collect();
    if lower.starts_with(b"inf") {
        return Num::Decimal(if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    if lower.starts_with(b"nan") {
        return Num::Decimal(f64::NAN);
    }
    Num::Decimal(0.0)
}

#[cfg(test)]
#[path = "number_tests.rs"]
mod tests;
