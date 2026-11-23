// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Numbers printed as Perl prints them, for the plugins whose Ensembl VEP
//! modules write one: `sprintf` with `%g` (LoFtool.pm:89), `%.2f` and `%.3g`
//! (pLI.pm:110) or `%.4f` (gnomADc.pm:175), and a number joined into a string
//! (LoFTEE's `'GERP_DIST:' . $gerp_dist`), which Perl writes with fifteen
//! significant digits.
//!
//! Perl's `sprintf` is C's: `%.Nf` prints `N` digits after the point, and `%.Pg`
//! prints `P` significant digits in the fixed form when the decimal exponent `X`
//! of the rounded value satisfies `-4 <= X < P`, else in the exponent form with
//! a sign and at least two exponent digits, and drops the trailing zeros of the
//! fraction either way (`6.89e-05`, `0.000138`, `1e+06`). Rust's `{:e}` writes
//! the exponent as `e-5`, so the exponent form is assembled here.
//!
//! The text those modules hand `sprintf` is a table cell, which Perl numifies
//! first: `numify` reads a cell as Perl does. A parameter value met in
//! arithmetic or a numeric comparison is read by `numify_arithmetic`, which
//! differs on a few `inf` and `nan` spellings.
//!
//! Perl multiplies and divides two integers as integers where the result
//! allows it, so a product of whole numbers past `2^53` is exact where the
//! same product in doubles has rounded: `PerlNumber` is an operand of those
//! two operators, read from a text or carried out of an earlier operation.

/// `sprintf("%g", value)`: six significant digits, trailing zeros dropped.
pub fn sprintf_g(value: f64) -> String {
    percent_g(value, 6)
}

/// The number Perl makes of a cell handed to `sprintf` (LoFtool.pm:89,
/// pLI.pm:145), with the text left after it, empty when the cell was read
/// whole. Perl's `my_atof3` reads leading blanks and one sign, then runs its
/// own `inf` and `nan` grammar ([`grok_infnan`]) from that sign when it is a
/// `-`, or after a `+` with a sign of its own, and keeps that reading only
/// where the grammar consumed text: `inf`, `infinity` and Visual C's `1.#INF`
/// and `1.#IND` whatever follows them, and a `nan`, with its `q` or `s` marks
/// and Visual C's `1.#` in front, only when text follows it. A `0x` or `0b`
/// form after the sign is `0` whatever the sign, as Perl refuses `strtod`
/// those. The rest goes to C's `strtod`, which skips blanks of its own, then
/// reads a sign of its own, `inf` or `nan` again, a hexadecimal number, or the
/// longest decimal prefix: digits with at most one point, at least one digit,
/// and an exponent when it carries digits. So `1.5abc` is `1.5`, `1,5` is `1`,
/// `1.5e` is `1.5`, `--1` is `1`, `- 5` is `-5`, `0abc` is `0`, `-0x1A` is `0`,
/// `+ 0x1A` is `26` (the `0x` refusal looks at the character after the sign,
/// here the blank), `+-0x1A` is `-26`, `1e999` is infinite,
/// `1.#QNAN` is `1` where `1.#QNAN(123)` is not a number, and `qnan` is
/// nothing numeric. `None` where nothing numeric leads the cell (`NA`, `-`,
/// `.`), which Perl numifies to `0`.
pub(crate) fn numify(text: &str) -> Option<(f64, &str)> {
    /// Perl's grammar where it consumed text.
    fn consumed(unsigned: &str) -> Option<(f64, &str)> {
        grok_infnan(unsigned).and_then(|(value, rest, consumed)| consumed.then_some((value, rest)))
    }
    let text = text.trim_start_matches(perl_space);
    let (negative, body) = sign(text);
    let perl_grammar = if negative {
        consumed(body).map(|(value, rest)| (-value, rest))
    } else {
        signed(body, consumed)
    };
    if let Some(found) = perl_grammar {
        return Some(found);
    }
    if let Some(rest) = radix_mark(body) {
        return Some((0.0, rest));
    }
    let (value, rest) = strtod(body)?;
    Some((if negative { -value } else { value }, rest))
}

/// A blank of Perl's `isSPACE` set, skipped by its number readers before the
/// number: space, tab, line feed, vertical tab, form feed, carriage return.
fn perl_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\x0b' | '\x0c' | '\r')
}

/// The number Perl makes of `text` as the right operand of an arithmetic
/// operator or a numeric comparison, or as the left operand beside an
/// integer (`0.5 + $v`, `$v + 0`, `$v == 1`, `$v < 0`): the infinity or
/// not-a-number Perl's `grok_number` flags from the text after one sign,
/// whatever trails the word and whether or not the grammar consumed it, else
/// the number [`numify`] reads; so `qnan`, `1.#QNAN` and `-SNaN` are not a
/// number and `1.#INFI` is infinite here, where `sprintf` reads `0`, `1`, `0`
/// and `1`. A left operand beside a right operand Perl does not read as an
/// integer, a double such as `0.5`, is converted as `sprintf` converts it:
/// `$v + 0.5` and `$v < 0.5` take [`numify`]'s reading of those spellings.
/// `None` where nothing numeric leads `text`, which Perl numifies to `0`.
/// A zero keeps [`numify`]'s sign: `-0` is negative zero here, and the
/// integer zero to Perl on the right of an operator or beside an integer.
pub(crate) fn numify_arithmetic(text: &str) -> Option<f64> {
    let text = text.trim_start_matches(perl_space);
    let (negative, body) = sign(text);
    if let Some((value, _, _)) = grok_infnan(body) {
        // `grok_number` keeps the sign for an infinity and clears it for a nan.
        return Some(if negative && value.is_infinite() {
            -value
        } else {
            value
        });
    }
    numify(text).map(|(value, _)| value)
}

/// `2^53`: every whole number below it is a double, and Perl's operators read a
/// double out of arithmetic as an integer only below it
/// (`NV_PRESERVES_UV_BITS`, sv.c `S_sv_2iuv_common`).
const DOUBLE_INTEGERS: u64 = 1 << 53;

/// The magnitude of the lowest signed 64-bit integer, `2^63`.
const IV_MIN_MAGNITUDE: u64 = 1 << 63;

/// An operand of Perl's `*` and `/` as Perl 5.34 holds it on a 64-bit build
/// (`pp_multiply` and `pp_divide`, pp.c, under `PERL_PRESERVE_IVUV`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum PerlNumber {
    /// An integer of the signed or the unsigned 64-bit range, as its sign and
    /// its magnitude; a negative one reaches `-2^63`. A zero is negative only
    /// as the operand of an exponent spelling of negative zero (`-0e0`),
    /// whose double keeps the sign.
    Integer { negative: bool, magnitude: u64 },
    /// A double out of arithmetic, or of a text Perl reads as one. The
    /// operators read it as an integer where it is a whole number below
    /// `2^53`.
    Double(f64),
    /// The double of a text with a decimal point and no exponent, which the
    /// operators multiply and divide as a double whatever its value
    /// (`IS_NUMBER_IN_UV | IS_NUMBER_NOT_INT` sets no public integer flag).
    Point(f64),
    /// A text Perl's two converters read differently, such as `qnan`,
    /// `1.#QNAN` and `-0`: `value` is `grok_number`'s reading, taken for the
    /// right operand and for the left operand of an integer (`sv_2iv_flags`),
    /// and `left_of_double` is `my_atof3`'s, taken for the left operand of a
    /// double (`sv_2nv_flags`).
    Spelled { value: f64, left_of_double: f64 },
}

impl PerlNumber {
    /// The operand Perl makes of a text it reads whole, raising no `isn't
    /// numeric` warning (`grok_number`, numeric.c): blanks, one sign, then
    /// digits with at most one point and at least one digit, an exponent that
    /// carries digits, and blanks to the end; or the ten bytes `0 but true`;
    /// or an infinity or a not-a-number in a spelling [`infnan_whole`] takes;
    /// or a `-` and blanks to the end, which `grok_number` returns as a sign
    /// alone and Perl reads as zero. Digits alone are an integer where they
    /// fit 64 bits, and so is an exponent spelling whose double is a whole
    /// number that fits them (`1e3`, `1.5e1`); `None` for a text Perl
    /// converts only under the warning.
    pub(crate) fn from_text(text: &str) -> Option<Self> {
        if text == "0 but true" {
            return Some(Self::integer(false, 0));
        }
        let (negative, body) = sign(text.trim_start_matches(perl_space));
        if let Some(number) = Self::decimal_whole(negative, body) {
            return Some(number);
        }
        if negative && !body.is_empty() && body.chars().all(perl_space) {
            return Some(Self::integer(false, 0));
        }
        // `grok_number` keeps the sign for an infinity and clears it for a
        // not-a-number.
        let value = infnan_whole(body)?;
        let value = if negative && value.is_infinite() {
            -value
        } else {
            value
        };
        // `sv_2nv_flags` hands `my_atof3` the text as a C string, which ends
        // at its first NUL byte.
        let c_string = text.split('\0').next().unwrap_or(text);
        let left_of_double = numify(c_string).map_or(0.0, |(converted, _)| converted);
        Some(
            if left_of_double == value || (left_of_double.is_nan() && value.is_nan()) {
                Self::Double(value)
            } else {
                Self::Spelled {
                    value,
                    left_of_double,
                }
            },
        )
    }

    /// The operand of a decimal number that fills `body`, the text after the
    /// sign, up to trailing blanks; `None` where `body` is not one.
    fn decimal_whole(negative: bool, body: &str) -> Option<Self> {
        let bytes = body.as_bytes();
        let digits_end = |from: usize| {
            from + bytes[from..]
                .iter()
                .take_while(|byte| byte.is_ascii_digit())
                .count()
        };
        let integer_end = digits_end(0);
        let mut end = integer_end;
        let point = bytes.get(end) == Some(&b'.');
        if point {
            end = digits_end(end + 1);
            // `5.` and `.5` are numbers; a lone `.` is not.
            if end == 1 {
                return None;
            }
        }
        if end == 0 {
            return None;
        }
        let exponent = matches!(bytes.get(end), Some(b'e' | b'E'));
        if exponent {
            let exponent_digits = match bytes.get(end + 1) {
                Some(b'+' | b'-') => end + 2,
                _ => end + 1,
            };
            end = digits_end(exponent_digits);
            if end == exponent_digits {
                return None;
            }
        }
        if !body[end..].chars().all(perl_space) {
            return None;
        }
        let signed = |value: f64| if negative { -value } else { value };
        if !point && !exponent {
            return Some(match body[..end].parse::<u64>() {
                // `-0` is the integer zero, and negative zero where
                // `my_atof3` converts the text.
                Ok(0) if negative => Self::Spelled {
                    value: 0.0,
                    left_of_double: -0.0,
                },
                Ok(magnitude) if !negative || magnitude <= IV_MIN_MAGNITUDE => {
                    Self::integer(negative, magnitude)
                }
                Ok(magnitude) => Self::Double(-(magnitude as f64)),
                Err(_) => Self::Double(signed(body[..end].parse().ok()?)),
            });
        }
        let value = signed(body[..end].parse().ok()?);
        if !exponent {
            return Some(Self::Point(value));
        }
        // An exponent spelling is an integer where its double is a whole
        // number of the signed or the unsigned range (`S_sv_2iuv_common` and
        // `S_sv_2iuv_non_preserve`, sv.c).
        let fits = value >= -(IV_MIN_MAGNITUDE as f64) && value < 2.0 * IV_MIN_MAGNITUDE as f64;
        Some(if fits && value == value.trunc() {
            // The double stays beside the integer, so a negative zero keeps
            // its sign where the operand is multiplied as a double.
            Self::Integer {
                negative: value.is_sign_negative(),
                magnitude: value.abs() as u64,
            }
        } else {
            Self::Double(value)
        })
    }

    /// An integer of the given sign and magnitude, zero never negative.
    fn integer(negative: bool, magnitude: u64) -> Self {
        Self::Integer {
            negative: negative && magnitude != 0,
            magnitude,
        }
    }

    /// The operand as the sign and magnitude of the integer Perl's operators
    /// read (`SvIV_please_nomg`): an integer, or a double out of arithmetic,
    /// or `grok_number`'s reading of a text, that is a whole number below
    /// `2^53`.
    fn whole(self) -> Option<(bool, u64)> {
        match self {
            Self::Integer {
                negative,
                magnitude,
            } => Some((negative && magnitude != 0, magnitude)),
            Self::Double(value) | Self::Spelled { value, .. }
                if value.abs() < DOUBLE_INTEGERS as f64 && value == value.trunc() =>
            {
                Some((value < 0.0, value.abs() as u64))
            }
            Self::Double(_) | Self::Point(_) | Self::Spelled { .. } => None,
        }
    }

    /// The operand as a double: an integer rounds to the nearest one.
    pub(crate) fn double(self) -> f64 {
        match self {
            Self::Integer {
                negative: false,
                magnitude,
            } => magnitude as f64,
            Self::Integer {
                negative: true,
                magnitude,
            } => -(magnitude as f64),
            Self::Double(value) | Self::Point(value) | Self::Spelled { value, .. } => value,
        }
    }

    /// The operand as the double the left side of an operator is where the
    /// right side is a double: Perl tries the left side as an integer only
    /// beside an integer, and converts it as a double otherwise.
    fn double_left_of(self, right: Self) -> f64 {
        match self {
            Self::Spelled { left_of_double, .. } if right.whole().is_none() => left_of_double,
            _ => self.double(),
        }
    }

    /// `self * right` (`pp_multiply`): the exact product of two integers
    /// where it fits the unsigned range, or the signed range for a negative
    /// product, else the product of the two doubles.
    pub(crate) fn multiply(self, right: Self) -> Self {
        if let (Some((right_negative, right)), Some((left_negative, left))) =
            (right.whole(), self.whole())
        {
            if let Some(product) = left.checked_mul(right) {
                if left_negative == right_negative {
                    return Self::integer(false, product);
                }
                if product <= IV_MIN_MAGNITUDE {
                    return Self::integer(true, product);
                }
            }
        }
        Self::Double(self.double_left_of(right) * right.double())
    }

    /// `self / right` (`pp_divide`): the exact quotient of two integers where
    /// the dividend is past `2^53` and the divisor divides it, else the
    /// quotient of the two doubles; `None` for a zero divisor, where Perl
    /// dies with `Illegal division by zero`.
    pub(crate) fn divide(self, right: Self) -> Option<Self> {
        if let (Some((right_negative, divisor)), Some((left_negative, dividend))) =
            (right.whole(), self.whole())
        {
            if divisor == 0 {
                return None;
            }
            if dividend >= divisor && dividend > DOUBLE_INTEGERS && dividend % divisor == 0 {
                let quotient = dividend / divisor;
                return Some(if left_negative == right_negative {
                    Self::integer(false, quotient)
                } else if quotient <= IV_MIN_MAGNITUDE {
                    Self::integer(true, quotient)
                } else {
                    Self::Double(-(quotient as f64))
                });
            }
        }
        let divisor = right.double();
        if divisor == 0.0 {
            return None;
        }
        Some(Self::Double(self.double_left_of(right) / divisor))
    }
}

/// The infinity or not-a-number that fills `text`, the text after one sign,
/// with nothing Perl's own grammar counts as trailing and no warning of its
/// digit readers (`grok_infnan`, numeric.c, returning no
/// `IS_NUMBER_TRAILING`). After Visual C's `1.#` or `1#`, or with no prefix:
/// `inf` or `infinity`, then blanks, zeros standing before them where the
/// word is Visual C's `inf` (`1.#INF00`); or, after Visual C's prefix alone,
/// `ind`, then zeros; or `nan` with an optional `q` or `s` before it, then
/// blanks, or a `q` or `s` and blanks, or a payload in parentheses with
/// anything after them. A payload is decimal digits, or `0x` and hexadecimal
/// digits or `0b` and binary digits that name at most `0xffffffff`, one
/// underscore allowed between two digits. Blanks and zeros run to the end of
/// the text or to a NUL byte. A second sign is not read.
fn infnan_whole(text: &str) -> Option<f64> {
    /// Whether nothing but `skipped` characters stands before the end of
    /// `rest` or a NUL byte.
    fn ends_after(rest: &str, skipped: impl Fn(char) -> bool) -> bool {
        let rest = rest.trim_start_matches(skipped);
        rest.is_empty() || rest.starts_with('\0')
    }
    let (visual_c, word) = match text.strip_prefix('1') {
        Some(hashed) => (
            true,
            hashed
                .strip_prefix(".#")
                .or_else(|| hashed.strip_prefix('#'))?,
        ),
        None => (false, text),
    };
    if let Some(rest) = after(word, "inf") {
        let rest = match rest.as_bytes().first() {
            Some(b'i' | b'I') => after(rest, "inity")?,
            _ if visual_c => rest.trim_start_matches('0'),
            _ => rest,
        };
        return ends_after(rest, perl_space).then_some(f64::INFINITY);
    }
    if visual_c {
        if let Some(rest) = after(word, "ind") {
            return ends_after(rest, |c| c == '0').then_some(f64::NAN);
        }
    }
    let quiet_or_signalling = ['q', 's', 'Q', 'S'];
    let unmarked = word.strip_prefix(quiet_or_signalling).unwrap_or(word);
    let rest = after(unmarked, "nan")?;
    let tail = rest.strip_prefix(quiet_or_signalling).unwrap_or(rest);
    let Some(payload) = tail.strip_prefix('(') else {
        return ends_after(tail, perl_space).then_some(f64::NAN);
    };
    // Perl looks for the closing parenthesis from the payload's second
    // character on, and no further than a NUL byte.
    let close = payload
        .char_indices()
        .skip(1)
        .take_while(|(_, c)| *c != '\0')
        .find_map(|(at, c)| (c == ')').then_some(at))?;
    nan_payload(&payload[..close]).then_some(f64::NAN)
}

/// Whether Perl reads a not-a-number's payload whole and without a warning:
/// decimal digits; or `0x` or `0b` and digits of that radix, one underscore
/// allowed between two of them, naming at most `0xffffffff`, past which
/// `grok_hex` and `grok_bin` warn that the number is not portable.
fn nan_payload(payload: &str) -> bool {
    let radix = match payload.as_bytes() {
        [b'0', b'x' | b'X', third, ..] if third.is_ascii_hexdigit() => 16,
        [b'0', b'b' | b'B', b'0' | b'1', ..] => 2,
        _ => return payload.bytes().all(|byte| byte.is_ascii_digit()),
    };
    let digits = &payload[2..];
    let separated = digits.split('_').all(|run| !run.is_empty());
    separated
        && u64::from_str_radix(&digits.replace('_', ""), radix)
            .is_ok_and(|value| value <= u64::from(u32::MAX))
}

/// What C's `strtod` reads from the start of `text`: blanks (C's `isspace` in
/// the C locale, Perl's `isSPACE` set), a sign, then `inf`, `infinity` or
/// `nan` in any case, a hexadecimal number or the longest decimal prefix, with
/// the text after it; `None` when nothing numeric leads `text`.
fn strtod(text: &str) -> Option<(f64, &str)> {
    let text = text.trim_start_matches(perl_space);
    signed(text, |body| {
        inf_or_nan(body)
            .or_else(|| hex_prefix(body))
            .or_else(|| decimal_prefix(body))
    })
}

/// `read` applied after an optional sign, the value negated by a `-`.
fn signed(text: &str, read: impl Fn(&str) -> Option<(f64, &str)>) -> Option<(f64, &str)> {
    let (negative, body) = sign(text);
    let (value, rest) = read(body)?;
    Some((if negative { -value } else { value }, rest))
}

/// The text after the `0` of a `0x` or `0b` form in either case, which Perl
/// numifies to `0` without reading further.
fn radix_mark(text: &str) -> Option<&str> {
    match text.as_bytes() {
        [b'0', b'x' | b'X' | b'b' | b'B', ..] => Some(&text[1..]),
        _ => None,
    }
}

/// `text` after a leading `+` or `-`, and whether it was `-`.
fn sign(text: &str) -> (bool, &str) {
    match text.as_bytes().first() {
        Some(b'-') => (true, &text[1..]),
        Some(b'+') => (false, &text[1..]),
        _ => (false, text),
    }
}

/// `text` after `word`, when `word` leads it in any case.
fn after<'a>(text: &'a str, word: &str) -> Option<&'a str> {
    let head = text.as_bytes().get(..word.len())?;
    head.eq_ignore_ascii_case(word.as_bytes())
        .then(|| &text[word.len()..])
}

/// `inf`, `infinity` or `nan` leading `text` in any case, as C's `strtod`
/// reads them, with the text after the word.
fn inf_or_nan(text: &str) -> Option<(f64, &str)> {
    if let Some(rest) = after(text, "infinity").or_else(|| after(text, "inf")) {
        return Some((f64::INFINITY, rest));
    }
    after(text, "nan").map(|rest| (f64::NAN, rest))
}

/// Perl's own `inf` and `nan` grammar (`grok_infnan`) on `text`, its sign
/// already taken: the infinity or not-a-number it flags, the text after the
/// word, and whether it consumed the text, which `my_atof3` demands before
/// keeping the reading and `grok_number` does not. The words, in any case:
/// `inf` and `infinity`, `nan` with an optional `q` or `s` mark on either
/// side and an optional parenthesised payload, and Visual C's `1.#` (the
/// point optional) before `inf`, `ind` or the `nan` forms. Nothing is
/// consumed by an `inf` followed by an `i` short of `infinity`, or by a `nan`
/// that ends the text after its marks or opens a payload without a closing
/// parenthesis or with an empty one; a Perl `nan` consumes every character
/// after it otherwise, so the text returned after the word is the diagnostic's,
/// not what Perl leaves.
fn grok_infnan(text: &str) -> Option<(f64, &str, bool)> {
    let (visual_c, word) = match text.strip_prefix('1') {
        Some(hashed) => (
            true,
            hashed
                .strip_prefix(".#")
                .or_else(|| hashed.strip_prefix('#'))?,
        ),
        None => (false, text),
    };
    if let Some(rest) = after(word, "inf") {
        return Some(match (after(word, "infinity"), rest.as_bytes().first()) {
            (Some(rest), _) => (f64::INFINITY, rest, true),
            (None, Some(b'i' | b'I')) => (f64::INFINITY, rest, false),
            (None, _) => (f64::INFINITY, rest, true),
        });
    }
    if visual_c {
        if let Some(rest) = after(word, "ind") {
            return Some((f64::NAN, rest, true));
        }
    }
    let quiet_or_signalling = ['q', 's', 'Q', 'S'];
    let unmarked = word.strip_prefix(quiet_or_signalling).unwrap_or(word);
    let rest = after(unmarked, "nan")?;
    let tail = rest.strip_prefix(quiet_or_signalling).unwrap_or(rest);
    let consumed = match tail.strip_prefix('(') {
        None => !tail.is_empty(),
        // Perl looks for the closing parenthesis from the payload's second
        // character on.
        Some(payload) => payload.match_indices(')').any(|(close, _)| close > 0),
    };
    Some((f64::NAN, rest, consumed))
}

/// A hexadecimal number leading `text`, as `strtod` reads one: `0x` in either
/// case, hexadecimal digits with at most one point and at least one digit,
/// and a binary exponent (`p`, a sign, decimal digits) when it carries digits,
/// with the text after it, the value rounded once to a double. `None` when
/// no digit follows the `0x`, where `strtod` reads the `0` alone.
fn hex_prefix(text: &str) -> Option<(f64, &str)> {
    let bytes = text.as_bytes();
    if !matches!(bytes, [b'0', b'x' | b'X', ..]) {
        return None;
    }
    let digits_end = |from: usize| {
        from + bytes[from..]
            .iter()
            .take_while(|byte| byte.is_ascii_hexdigit())
            .count()
    };
    let integer_end = digits_end(2);
    let mut end = integer_end;
    let mut fraction: &[u8] = &[];
    if bytes.get(end) == Some(&b'.') {
        let fraction_end = digits_end(end + 1);
        if fraction_end > end + 1 || integer_end > 2 {
            fraction = &bytes[end + 1..fraction_end];
            end = fraction_end;
        }
    }
    if end == 2 {
        return None;
    }
    // The first 32 significant digits fill a `u128`, which a double rounds
    // once; a dropped digit sets the lowest bit when it is not zero, so the
    // rounding sees whether anything lay below the kept digits.
    const MANTISSA_DIGITS: usize = 32;
    let integer = &bytes[2..integer_end];
    let mut mantissa: u128 = 0;
    let mut kept = 0;
    let mut dropped: i32 = 0;
    let mut dropped_nonzero = false;
    for &byte in integer.iter().chain(fraction) {
        let digit = (byte as char).to_digit(16).expect("a hexadecimal digit");
        if kept == 0 && digit == 0 {
            continue;
        }
        if kept < MANTISSA_DIGITS {
            mantissa = mantissa * 16 + u128::from(digit);
            kept += 1;
        } else {
            dropped += 1;
            dropped_nonzero |= digit != 0;
        }
    }
    if dropped_nonzero {
        mantissa |= 1;
    }
    // Four bits per digit: the dropped digits stood below the kept ones, and
    // the fraction digits, kept or not, stood right of the point.
    let mut exponent = 4 * (dropped - fraction.len() as i32);
    if matches!(bytes.get(end), Some(b'p' | b'P')) {
        let (negative, power_start) = match bytes.get(end + 1) {
            Some(b'-') => (true, end + 2),
            Some(b'+') => (false, end + 2),
            _ => (false, end + 1),
        };
        let power_end = power_start
            + bytes[power_start..]
                .iter()
                .take_while(|byte| byte.is_ascii_digit())
                .count();
        if power_end > power_start {
            let power: i32 = text[power_start..power_end].parse().unwrap_or(i32::MAX);
            exponent = exponent.saturating_add(if negative { -power } else { power });
            end = power_end;
        }
    }
    Some((hex_value(mantissa, exponent), &text[end..]))
}

/// The power of two of the smallest subnormal double, the grid every
/// subnormal value sits on.
const SUBNORMAL_GRID: i32 = -1074;

/// `mantissa * 2^exponent` rounded once to a double. A normal result rounds
/// when the mantissa becomes a double and is then scaled exactly; a result
/// below the smallest normal double is rounded onto the subnormal grid
/// directly, half to even, since a double's 53 bits would round it a second
/// time when scaled.
fn hex_value(mantissa: u128, exponent: i32) -> f64 {
    let bits_below_grid = SUBNORMAL_GRID.saturating_sub(exponent);
    if bits_below_grid > 0 {
        let (on_grid, under) = if bits_below_grid < u128::BITS as i32 {
            (
                mantissa >> bits_below_grid,
                mantissa & ((1u128 << bits_below_grid) - 1),
            )
        } else {
            (0, mantissa)
        };
        if on_grid < 1 << 52 {
            let round_up = if bits_below_grid > u128::BITS as i32 {
                // The whole mantissa lies under half a grid step.
                false
            } else {
                let half = 1u128 << (bits_below_grid - 1);
                under > half || (under == half && on_grid & 1 == 1)
            };
            return scale_by_two((on_grid + u128::from(round_up)) as f64, SUBNORMAL_GRID);
        }
    }
    scale_by_two(mantissa as f64, exponent)
}

/// `value * 2^exponent`, scaled in steps a double holds exactly, so a value
/// beyond the double's range reaches infinity or zero as the product does and
/// a zero stays a zero whatever the exponent.
fn scale_by_two(mut value: f64, mut exponent: i32) -> f64 {
    const STEP: i32 = 1000;
    while exponent != 0 && value != 0.0 && value.is_finite() {
        let step = exponent.clamp(-STEP, STEP);
        value *= 2f64.powi(step);
        exponent -= step;
    }
    value
}

/// The longest decimal prefix of `text`, as `strtod` reads one, parsed, with
/// the text after it.
fn decimal_prefix(text: &str) -> Option<(f64, &str)> {
    let bytes = text.as_bytes();
    let digits_end = |from: usize| {
        from + bytes[from..]
            .iter()
            .take_while(|byte| byte.is_ascii_digit())
            .count()
    };
    let integer_end = digits_end(0);
    let mut end = integer_end;
    if bytes.get(end) == Some(&b'.') {
        let fraction_end = digits_end(end + 1);
        // `5.` and `.5` are numbers; a lone `.` is not.
        if fraction_end > end + 1 || integer_end > 0 {
            end = fraction_end;
        }
    }
    if end == 0 {
        return None;
    }
    if matches!(bytes.get(end), Some(b'e' | b'E')) {
        let exponent_digits = match bytes.get(end + 1) {
            Some(b'+' | b'-') => end + 2,
            _ => end + 1,
        };
        let exponent_end = digits_end(exponent_digits);
        if exponent_end > exponent_digits {
            end = exponent_end;
        }
    }
    let value = text[..end].parse().ok()?;
    Some((value, &text[end..]))
}

/// A floating-point number in Perl's string context: `%.15g`, the fifteen
/// significant digits of a double (`27/113` is `0.238938053097345`, `1e15` is
/// `1e+15`), and `0` for a zero of either sign, which `sprintf("%.15g")` writes
/// as `-0`.
pub(crate) fn perl_number(value: f64) -> String {
    vep_core::variation::perl_number_string(value)
}

/// `sprintf("%.3g", value)`: three significant digits, trailing zeros dropped.
pub fn sig3(value: f64) -> String {
    percent_g(value, 3)
}

/// `sprintf("%.<decimals>f", value)`: a fixed number of digits after the point.
pub fn fixed(value: f64, decimals: usize) -> String {
    if let Some(text) = non_finite(value) {
        return text;
    }
    fixed_by_integers(value, decimals).unwrap_or_else(|| format!("{value:.decimals$}"))
}

/// Decimal places [`fixed_by_integers`] handles: `10^decimals` times a
/// 53-bit mantissa stays under `2^83`, so the scaled value fits a `u128`
/// beside a binary exponent of up to [`INTEGER_PATH_MAX_EXPONENT`].
const INTEGER_PATH_MAX_DECIMALS: usize = 9;

/// Binary exponents [`fixed_by_integers`] handles: the scaled value, under
/// `2^83`, shifted left by this many bits stays under `2^123`.
const INTEGER_PATH_MAX_EXPONENT: i32 = 40;

/// `%.<decimals>f` of a finite double through integer arithmetic on its exact
/// binary value, the digits C's `printf` prints: the value is `m * 2^e` with
/// `m` an integer, so `m * 10^decimals * 2^e` is rounded to an integer with a
/// tie going to the even neighbour, and that integer is written with the point
/// `decimals` digits from its end. A negative value keeps its sign even when
/// every digit is zero, as `-0.0000` prints for `-0.0`. `None` for a decimal
/// count of zero or above [`INTEGER_PATH_MAX_DECIMALS`], and for a magnitude
/// of `2^93` or more, whose exponent `e` exceeds [`INTEGER_PATH_MAX_EXPONENT`];
/// the caller formats those the general way.
fn fixed_by_integers(value: f64, decimals: usize) -> Option<String> {
    if decimals == 0 || decimals > INTEGER_PATH_MAX_DECIMALS {
        return None;
    }
    let bits = value.to_bits();
    let exponent_bits = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1u64 << 52) - 1);
    // A subnormal has no implicit leading bit and the smallest exponent.
    let (mantissa, exponent) = if exponent_bits == 0 {
        (fraction, -1074)
    } else {
        (fraction | (1u64 << 52), exponent_bits - 1075)
    };
    if exponent > INTEGER_PATH_MAX_EXPONENT {
        return None;
    }
    let scaled = u128::from(mantissa) * 10u128.pow(decimals as u32);
    let rounded = if exponent >= 0 {
        scaled << exponent
    } else if exponent <= -84 {
        // `scaled` is below `2^83`, so the value is under half a unit.
        0
    } else {
        let shift = exponent.unsigned_abs();
        let quotient = scaled >> shift;
        let remainder = scaled & ((1u128 << shift) - 1);
        let half = 1u128 << (shift - 1);
        let round_up = remainder > half || (remainder == half && quotient & 1 == 1);
        quotient + u128::from(round_up)
    };
    // The digits are written backwards into a stack buffer, the fraction first
    // and zero-padded to `decimals`, so the returned string is allocated once at
    // its exact length (the rounded value is below `2^123`, at most 38 digits).
    let mut digits = [b'0'; 64];
    let mut start = digits.len();
    let mut rest = rounded;
    for _ in 0..decimals {
        start -= 1;
        digits[start] = b'0' + (rest % 10) as u8;
        rest /= 10;
    }
    start -= 1;
    digits[start] = b'.';
    loop {
        start -= 1;
        digits[start] = b'0' + (rest % 10) as u8;
        rest /= 10;
        if rest == 0 {
            break;
        }
    }
    if value.is_sign_negative() {
        start -= 1;
        digits[start] = b'-';
    }
    let text = std::str::from_utf8(&digits[start..]).expect("ASCII digits, a point and a sign");
    Some(text.to_owned())
}

/// `sprintf("%.<precision>g", value)` for a precision of at least one.
fn percent_g(value: f64, precision: usize) -> String {
    if let Some(text) = non_finite(value) {
        return text;
    }
    // Rust's exponent form rounds to `precision` significant digits, so its
    // exponent is the `X` of C's rule for the rounded value.
    let exponent_form = format!("{value:.*e}", precision - 1);
    let (mantissa, exponent) = exponent_form
        .split_once('e')
        .expect("Rust's exponent form carries an e");
    let exponent: i32 = exponent.parse().expect("an integer exponent");
    if -4 <= exponent && exponent < precision as i32 {
        let decimals = (precision as i32 - 1 - exponent) as usize;
        return trim_fraction(&format!("{value:.decimals$}")).to_string();
    }
    let sign = if exponent < 0 { '-' } else { '+' };
    format!(
        "{}e{sign}{:02}",
        trim_fraction(mantissa),
        exponent.unsigned_abs()
    )
}

/// Perl prints `Inf`, `-Inf` and `NaN` where C prints `inf` and `nan`.
fn non_finite(value: f64) -> Option<String> {
    if value.is_nan() {
        Some("NaN".to_string())
    } else if value.is_infinite() {
        Some(if value > 0.0 { "Inf" } else { "-Inf" }.to_string())
    } else {
        None
    }
}

/// A decimal text without the trailing zeros of its fraction, and without the
/// point when nothing of the fraction remains.
fn trim_fraction(text: &str) -> &str {
    if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        text
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Perl's `%g` on the LoFtool scores file's exponent and fixed forms
    /// (`6.89E-05` and `1.38E-04` are the two spellings the file carries) and on
    /// the whole numbers and halves the other plugins print.
    #[test]
    fn percent_g_matches_perl_for_exponent_and_fixed_forms() {
        assert_eq!(sprintf_g(6.89e-05), "6.89e-05");
        assert_eq!(sprintf_g(1.38e-04), "0.000138");
        assert_eq!(sprintf_g(1.0), "1");
        assert_eq!(sprintf_g(0.5), "0.5");
        assert_eq!(sprintf_g(0.0), "0");
        assert_eq!(sprintf_g(100000.0), "100000");
        assert_eq!(sprintf_g(1_000_000.0), "1e+06");
        assert_eq!(sprintf_g(999_999.5), "1e+06");
        assert_eq!(sprintf_g(0.0001), "0.0001");
        assert_eq!(sprintf_g(0.00009999), "9.999e-05");
        assert_eq!(sprintf_g(123.456789), "123.457");
        assert_eq!(sprintf_g(-2.5e-7), "-2.5e-07");
        assert_eq!(sprintf_g(1.5e300), "1.5e+300");
    }

    /// Perl's numification of a cell, each case with the value its `sprintf`
    /// prints: a number read whole, a number before other text, a second sign
    /// read by `strtod` or by Perl's own grammar after a `+`, a blank after the
    /// sign, which `strtod` skips before its own reading, `inf` and `nan` by
    /// their first letters and in Visual C's spelling, a `0x` or `0b` form as
    /// `0` and as `strtod`'s hexadecimal after a second sign or a blank, and
    /// cells with nothing numeric in front.
    #[test]
    fn numify_reads_a_cell_as_perl_does() {
        fn read(text: &str) -> Option<(String, &str)> {
            numify(text).map(|(value, rest)| (sprintf_g(value), rest))
        }
        for (cell, value) in [
            ("1.26E-01", "0.126"),
            ("6.89E-05", "6.89e-05"),
            ("+1.5", "1.5"),
            ("-.5", "-0.5"),
            ("5.", "5"),
            ("00012", "12"),
            ("1.e5", "100000"),
            ("\x0b\x0c20", "20"),
            ("- 5", "-5"),
            ("-\x0b5", "-5"),
            ("- -5", "5"),
            ("+ 0x1A", "26"),
            ("- inf", "-Inf"),
            ("Inf", "Inf"),
            ("-Infinity", "-Inf"),
            ("NaN", "NaN"),
            ("-0", "-0"),
            ("1e999", "Inf"),
            ("1e-400", "0"),
        ] {
            assert_eq!(read(cell), Some((value.to_string(), "")), "{cell}");
        }
        for (cell, value, rest) in [
            ("1.5abc", "1.5", "abc"),
            ("12.5%", "12.5", "%"),
            ("1,5", "1", ",5"),
            ("1_000", "1", "_000"),
            ("1.2.3", "1.2", ".3"),
            ("1.5e", "1.5", "e"),
            ("1.5e+", "1.5", "e+"),
            ("1.5E-3x", "0.0015", "x"),
            ("0x1A", "0", "x1A"),
            ("-0x1A", "0", "x1A"),
            ("-0X1A", "0", "X1A"),
            ("-0b11", "0", "b11"),
            ("0x1p3", "0", "x1p3"),
            ("+-0x1A", "-26", ""),
            ("--0x1A", "26", ""),
            ("-+0x1A", "-26", ""),
            ("+-0X1A", "-26", ""),
            ("+-0x1p3", "-8", ""),
            ("+-0x1p+2", "-4", ""),
            ("+-0x1Ap-2", "-6.5", ""),
            ("+-0x1.8p1", "-3", ""),
            ("+-0x.8", "-0.5", ""),
            ("+-0x1.", "-1", ""),
            ("+-0x1pz", "-1", "pz"),
            ("+-0xg", "-0", "xg"),
            ("+-0x", "-0", "x"),
            ("+-0x.", "-0", "x."),
            ("+-0b11", "-0", "b11"),
            ("+-0xABCDEF0123456789ABCDEF0123456789A", "-3.65388e+39", ""),
            ("+-0x000000000000000000000000000000001", "-1", ""),
            (
                "+-0x0.000000000000000000000000000000001",
                "-1.83671e-40",
                "",
            ),
            ("+-0x0p1030", "-0", ""),
            ("+-0x0p99999999999", "-0", ""),
            ("+-0x10p-1078", "-4.94066e-324", ""),
            ("+-0x2p-1024", "-1.11254e-308", ""),
            (
                "+-0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFp-1100",
                "-2.50521e-293",
                "",
            ),
            ("+-0x1p99999", "-Inf", ""),
            ("+-0x1.8p2147483647", "-Inf", ""),
            ("+-0x1p-99999999999", "-0", ""),
            ("+-0x.8p-99999999999", "-0", ""),
            ("+-0x1.0p-2147483647", "-0", ""),
            ("+-0xA0000000000001p-1128", "-1.4822e-323", ""),
            ("+-0x1p-1074", "-4.94066e-324", ""),
            ("+-0x1p-1075", "-0", ""),
            ("+-0x3p-1075", "-9.88131e-324", ""),
            ("+-0x5p-1075", "-9.88131e-324", ""),
            ("+-0x7p-1075", "-1.97626e-323", ""),
            ("+-0xFFFFFFFFFFFFF8p-1130", "-4.94066e-324", ""),
            ("+-0x1.FFFFFFFFFFFFFp-1023", "-2.22507e-308", ""),
            ("+-0x1p-1022", "-2.22507e-308", ""),
            ("+-0x10000000000000000000000000000001p-1200", "-0", ""),
            ("+-1.#INF", "-Inf", ""),
            ("+-1.#IND", "NaN", ""),
            ("++1.#INF", "Inf", ""),
            ("-+1.#INF", "-1", "#INF"),
            ("0abc", "0", "abc"),
            ("-0abc", "-0", "abc"),
            ("- 5abc", "-5", "abc"),
            ("- 1.#INF", "-1", "#INF"),
            ("--1", "1", ""),
            ("+-1", "-1", ""),
            ("-+1", "-1", ""),
            ("++1", "1", ""),
            ("+-0", "-0", ""),
            ("--0", "0", ""),
            ("--.5", "0.5", ""),
            ("+-inf", "-Inf", ""),
            ("--inf", "Inf", ""),
            ("infx", "Inf", "x"),
            ("Infinit", "Inf", "init"),
            ("1.#Infinity", "Inf", ""),
            ("1.#INFI", "1", "#INFI"),
            ("1.#INFINIT", "1", "#INFINIT"),
            ("nan(123)", "NaN", "(123)"),
            ("NaN123", "NaN", "123"),
            ("nanq", "NaN", "q"),
            ("nan(", "NaN", "("),
            ("qnan5", "NaN", "5"),
            ("qnan(1)", "NaN", "(1)"),
            ("1.#QNAN(123)", "NaN", "(123)"),
            ("1.#QNAN0", "NaN", "0"),
            ("1.#nan00", "NaN", "00"),
            ("1.#INF", "Inf", ""),
            ("1#INF", "Inf", ""),
            ("-1.#INF", "-Inf", ""),
            ("1.#INFabc", "Inf", "abc"),
            ("1.#IND", "NaN", ""),
            ("1.#IND5", "NaN", "5"),
            ("1.#QNAN", "1", "#QNAN"),
            ("1#SNAN", "1", "#SNAN"),
            ("1.#NAN", "1", "#NAN"),
            ("1.#nanq", "1", "#nanq"),
            ("-1.#QNAN", "-1", "#QNAN"),
            ("1.#QNAN(", "1", "#QNAN("),
            ("1.#QNAN()", "1", "#QNAN()"),
            ("2.#INF", "2", "#INF"),
            ("--1.#INF", "1", "#INF"),
            ("1e5.5", "100000", ".5"),
            ("1ß", "1", "ß"),
        ] {
            assert_eq!(read(cell), Some((value.to_string(), rest)), "{cell}");
        }
        for cell in [
            "NA", "n/a", "-", "+", ".", "e5", ".e5", "abc", "None", "---1", "+-+1", "..5", ",5",
            "_1", "in", "na", "١٢", "ß1", "", "qnan", "snan", "snanq", "-SNaN", "qnan()", "sqnan",
        ] {
            assert_eq!(read(cell), None, "{cell}");
        }
        assert_eq!(numify(" 1.5"), Some((1.5, "")));
        // The 33rd digit decides the tie the 32 kept digits leave.
        let (tie, rest) = numify("+-0x800000000000040000000000000000001").unwrap();
        assert_eq!(
            (fixed(tie, 2), rest),
            (
                "-2722258935367508312169906666768733044736.00".to_string(),
                ""
            )
        );
    }

    /// On the right of an arithmetic operator, or on its left beside an
    /// integer, Perl takes `grok_number`'s infinity or not-a-number whether
    /// or not its grammar consumed the text, so the `nan` spellings `sprintf`
    /// reads as `0` or `1` are not a number here, `1.#INFI` is infinite, and
    /// every other text reads as the cell does.
    #[test]
    fn numify_arithmetic_takes_grok_numbers_inf_and_nan() {
        for text in [
            "qnan", "snan", "snanq", "-SNaN", "1.#QNAN", "1#SNAN", "1.#NAN", "1.#nanq", "-1.#QNAN",
            "1.#QNAN(", "nan", "nanq", "nan(123)", "nan(abc)", "1.#IND", "+-1.#IND", "- nan",
        ] {
            assert!(numify_arithmetic(text).unwrap().is_nan(), "{text}");
        }
        for (text, value) in [
            ("1.#INFI", f64::INFINITY),
            ("-1.#INFINIT", f64::NEG_INFINITY),
            ("inf", f64::INFINITY),
            ("+-inf", f64::NEG_INFINITY),
            ("- inf", f64::NEG_INFINITY),
            ("--1.#INF", 1.0),
            ("- 1.#INF", -1.0),
            ("1.5abc", 1.5),
            ("--1", 1.0),
            ("- 5", -5.0),
            ("0x1A", 0.0),
            ("+-0x1A", -26.0),
            ("1e999", f64::INFINITY),
        ] {
            assert_eq!(numify_arithmetic(text), Some(value), "{text}");
        }
        for text in ["", "NA", "sqnan", "---1", "."] {
            assert_eq!(numify_arithmetic(text), None, "{text}");
        }
    }

    /// A `PerlNumber` as a `perl` one-liner shows the same operand: `I` and
    /// the integer, or `N` and the bits of the double in hexadecimal, `nan`
    /// for a not-a-number.
    fn shown(number: PerlNumber) -> String {
        match number {
            PerlNumber::Integer {
                negative,
                magnitude,
            } => format!("I {}{magnitude}", if negative { "-" } else { "" }),
            other if other.double().is_nan() => "N nan".to_string(),
            other => format!("N {:#018x}", other.double().to_bits()),
        }
    }

    /// The product of the `*`-joined texts of `factors`, left to right.
    fn product_of(factors: &str) -> PerlNumber {
        factors
            .split('*')
            .map(|text| PerlNumber::from_text(text).unwrap())
            .reduce(PerlNumber::multiply)
            .unwrap()
    }

    /// The operand of a text, each as Perl 5.34 on a 64-bit build holds the
    /// product of the text and `"1"` under `use warnings`: an integer for
    /// digits alone and for an exponent spelling of a whole number, a double
    /// for a text with a point, for a whole number past 64 bits and for every
    /// spelling of an infinity or a not-a-number Perl takes whole; and `None`
    /// for each text the same product raises `isn't numeric` on.
    #[test]
    fn perl_number_reads_a_text_as_perls_arithmetic_does() {
        let one = PerlNumber::from_text("1").unwrap();
        for (text, expected) in [
            ("94906267", "I 94906267"),
            (" 94906267 ", "I 94906267"),
            ("+94906267", "I 94906267"),
            ("094906267", "I 94906267"),
            ("-0", "I 0"),
            ("18446744073709551615", "I 18446744073709551615"),
            ("-9223372036854775808", "I -9223372036854775808"),
            ("3e0", "I 3"),
            ("+.5e1", "I 5"),
            ("1.e5", "I 100000"),
            ("1e17", "I 100000000000000000"),
            ("1e19", "I 10000000000000000000"),
            ("72057594037927937e0", "I 72057594037927936"),
            ("1e-400", "I 0"),
            ("-1e-400", "I 0"),
            ("0 but true", "I 0"),
            ("-\t", "I 0"),
            ("1.0", "N 0x3ff0000000000000"),
            ("94906267.", "N 0x4196a09e6c000000"),
            (".5", "N 0x3fe0000000000000"),
            ("-0.0", "N 0x8000000000000000"),
            ("1e20", "N 0x4415af1d78b58c40"),
            ("18446744073709551616", "N 0x43f0000000000000"),
            ("-9223372036854775809", "N 0xc3e0000000000000"),
            ("1.5e0", "N 0x3ff8000000000000"),
            ("1e400", "N 0x7ff0000000000000"),
            ("-1e400", "N 0xfff0000000000000"),
            ("inf", "N 0x7ff0000000000000"),
            ("-Infinity", "N 0xfff0000000000000"),
            ("+INF", "N 0x7ff0000000000000"),
            ("infinity ", "N 0x7ff0000000000000"),
            ("1.#INF", "N 0x7ff0000000000000"),
            ("-1.#INF", "N 0xfff0000000000000"),
            ("1#INF", "N 0x7ff0000000000000"),
            ("1.#infinity", "N 0x7ff0000000000000"),
            ("1.#INF0", "N 0x7ff0000000000000"),
            ("-1#INF00 ", "N 0xfff0000000000000"),
            ("nan", "N nan"),
            ("-NaN", "N nan"),
            ("nan ", "N nan"),
            ("nanq", "N nan"),
            ("NaNS", "N nan"),
            ("qnan", "N nan"),
            ("sNaN", "N nan"),
            ("nan(1)", "N nan"),
            ("nan(0x7ff)", "N nan"),
            ("nan(0b101)", "N nan"),
            ("nan(1)x", "N nan"),
            ("nan(1)(2)", "N nan"),
            ("1.#IND", "N nan"),
            ("1.#IND0", "N nan"),
            ("1.#QNAN", "N nan"),
            ("1#NAN", "N nan"),
            ("nan(0xffffffff)", "N nan"),
            ("nan(0x1_0)", "N nan"),
        ] {
            let read = PerlNumber::from_text(text).map(|number| shown(number.multiply(one)));
            assert_eq!(read.as_deref(), Some(expected), "{text:?}");
        }
        for text in [
            "",
            " ",
            "high",
            "0.5x",
            "0.5e",
            "0x10",
            "1_000",
            "0,5",
            "infinite",
            "nanx",
            "infx",
            "infinit",
            "infinityx",
            "NaN(",
            "nan()",
            "nan(1",
            "nan(a)",
            "nan(1.5)",
            "nan(0x100000000)",
            "nan(0x1__0)",
            "1.#INFx",
            "1.#INF 0",
            "1.#INF0x",
            "1.#INFINITY0",
            "inf0",
            "1.#INDx",
            "1.#IND ",
            "1.#IN",
            "0 but true ",
            " 0 but true",
            "0 BUT TRUE",
            "+-inf",
            "--nan",
            ".",
            "+",
            "-",
            "+ ",
            "1e",
            ".e1",
            "++1",
            "- 1",
            "1e5.5",
            "0.5\0",
            "\u{a0}0.5",
        ] {
            assert_eq!(PerlNumber::from_text(text), None, "{text:?}");
        }
        // A point without an exponent stays a double whatever its value, and
        // an exponent spelling of a whole number is an integer.
        assert_eq!(PerlNumber::from_text("3.0"), Some(PerlNumber::Point(3.0)));
        assert_eq!(
            PerlNumber::from_text("3e0"),
            Some(PerlNumber::Integer {
                negative: false,
                magnitude: 3
            })
        );
        assert_eq!(
            PerlNumber::from_text("1.5e0"),
            Some(PerlNumber::Double(1.5))
        );
    }

    /// Products and quotients of texts, each as Perl 5.34 on a 64-bit build
    /// computes it over fresh strings: two integers multiply exactly while
    /// the product fits the unsigned range, or the signed range when it is
    /// negative, and as doubles past it; a double out of arithmetic that is a
    /// whole number below `2^53` is an integer again; a text with a point is
    /// never one; a quotient is an integer where the dividend is past `2^53`
    /// and the divisor divides it; and a zero divisor is where Perl dies. A
    /// spelling Perl's two converters read differently is a not-a-number on
    /// the right and beside an integer, and `0` (`qnan`) or `1` (`1.#QNAN`)
    /// on the left of a double, where the text ends at a NUL byte; `-0` is
    /// negative zero there alone, and an exponent spelling of negative zero
    /// beside any double.
    #[test]
    fn perl_number_multiplies_and_divides_as_perl_does() {
        for (factors, expected) in [
            ("300000007*300000007", "I 90000004200000049"),
            ("3037000500*3037000500", "I 9223372037000250000"),
            ("4294967295*4294967297", "I 18446744073709551615"),
            ("4294967296*4294967296", "N 0x43f0000000000000"),
            ("4294967297*4294967297", "N 0x43f0000000200000"),
            ("-3037000500*3037000500", "N 0xc3e0000000011578"),
            ("-4294967296*2147483648", "I -9223372036854775808"),
            ("-4294967296*2147483649", "N 0xc3e0000000200000"),
            ("-4294967295*-4294967297", "I 18446744073709551615"),
            ("300000007*300000007*200", "I 18000000840000009800"),
            ("300000007*300000007*210", "N 0x43f064a4aaa6a2e2"),
            ("300000007.0*300000007", "N 0x4373fbe86e81fea3"),
            ("300000007*300000007.0", "N 0x4373fbe86e81fea3"),
            ("3.00000007e8*300000007", "I 90000004200000049"),
            ("300000007e0*300000007e0", "I 90000004200000049"),
            ("0.5*4*300000007*300000007", "I 180000008400000098"),
            ("1.5*2*300000007*300000007", "I 270000012600000147"),
            ("72057594037927937*3", "I 216172782113783811"),
            ("72057594037927937*1.5", "N 0x4378000000000000"),
            ("72057594037927937e0*3", "I 216172782113783808"),
            ("18446744073709551615*1", "I 18446744073709551615"),
            ("18446744073709551616*1", "N 0x43f0000000000000"),
            ("-0.0*5", "N 0x8000000000000000"),
            ("0*-5", "I 0"),
            ("inf*0", "N nan"),
            ("nan*2", "N nan"),
            ("-inf*-1", "N 0x7ff0000000000000"),
            ("qnan*2", "N nan"),
            ("qnan*0.5", "N 0x0000000000000000"),
            ("0.5*qnan", "N nan"),
            ("snan*0.25*4", "I 0"),
            ("qnan(1)*0.5", "N nan"),
            ("1.#QNAN*2", "N nan"),
            ("1.#QNAN*0.5", "N 0x3fe0000000000000"),
            ("0.5*1.#QNAN", "N nan"),
            ("-1.#QNAN*0.5", "N 0xbfe0000000000000"),
            ("1.#QNAN *0.5", "N nan"),
            ("1.#NAN*1.5*2", "N 0x4008000000000000"),
            ("qnan\0*0.5", "N 0x0000000000000000"),
            ("1.#QNAN\0*0.5", "N 0x3fe0000000000000"),
            ("qnan\0*2", "N nan"),
            ("-0e0*5", "I 0"),
            ("-0e0*0.5", "N 0x8000000000000000"),
            ("0.5*-0e0", "N 0x8000000000000000"),
            ("-1e-400*0.5", "N 0x8000000000000000"),
            ("-0*5", "I 0"),
            ("0.5*-0", "N 0x0000000000000000"),
            ("-0*0.5", "N 0x8000000000000000"),
            ("-5*-0", "I 0"),
            ("-5*-0*0.5", "N 0x0000000000000000"),
            ("-00*0.5*2", "I 0"),
        ] {
            assert_eq!(shown(product_of(factors)), expected, "{factors}");
        }
        for (dividend, divisor, expected) in [
            ("72057594037927937*3", "3*1", Some("I 72057594037927937")),
            ("72057594037927937*9", "3*3", Some("I 72057594037927937")),
            ("72057594037927937*3", "7*1", Some("N 0x435b6db6db6db6db")),
            ("300000007*300000007", "300000007*1", Some("I 300000007")),
            ("67108864*134217728", "2*1", Some("N 0x4330000000000000")),
            ("134217728*67108865", "134217728*1", Some("I 67108865")),
            ("2*3", "3*1", Some("N 0x4000000000000000")),
            ("3*1", "72057594037927937*3", Some("N 0x3c70000000000000")),
            ("-72057594037927937*3", "3*1", Some("I -72057594037927937")),
            (
                "18446744073709551615*1",
                "-1*1",
                Some("N 0xc3f0000000000000"),
            ),
            ("72057594037927937*3", "3.0*1", Some("I 72057594037927937")),
            ("72057594037927937*3", "1.5*2", Some("I 72057594037927937")),
            ("4294967296*4294967296", "4*1", Some("N 0x43d0000000000000")),
            (
                "-9223372036854775808*1",
                "-1*1",
                Some("I 9223372036854775808"),
            ),
            (
                "-9223372036854775808*1",
                "1*1",
                Some("I -9223372036854775808"),
            ),
            ("5*1", "0*7", None),
            ("5*1", "0.0*7", None),
            ("5*1", "-0.0*7", None),
            ("5*1", "inf*1", Some("N 0x0000000000000000")),
            ("5*1", "nan*1", Some("N nan")),
            ("0*1", "5*1", Some("N 0x0000000000000000")),
            ("inf*1", "inf*1", Some("N nan")),
            ("8*1", "qnan*0.5", None),
            ("8*1", "qnan*2", Some("N nan")),
            ("8*1", "1.#QNAN*0.5", Some("N 0x4030000000000000")),
        ] {
            let quotient = product_of(dividend).divide(product_of(divisor)).map(shown);
            assert_eq!(quotient.as_deref(), expected, "{dividend} / {divisor}");
        }
    }

    /// pLI's `%.2f` and `%.3g` on two values of its file, and gnomADc's `%.4f`
    /// on an integer depth.
    #[test]
    fn two_decimals_and_three_significant_figures() {
        assert_eq!(fixed(0.17335234048116, 2), "0.17");
        assert_eq!(sig3(0.17335234048116), "0.173");
        assert_eq!(fixed(1.33038194561114e-19, 2), "0.00");
        assert_eq!(sig3(1.33038194561114e-19), "1.33e-19");
        assert_eq!(fixed(12345.0, 4), "12345.0000");
        assert_eq!(fixed(30.7, 4), "30.7000");
        assert_eq!(sig3(1.0), "1");
        assert_eq!(sig3(0.999_9), "1");
        assert_eq!(sig3(1234.0), "1.23e+03");
    }

    /// `fixed` writes what C's `printf` writes, through the integer path and,
    /// at zero or ten decimals, the general formatter: a tie on the exact
    /// binary value goes to the even neighbour, a negative value keeps its sign
    /// when every digit rounds to zero, and the values a coverage file carries
    /// (exponent forms, integer depths, fractions) round to four decimals.
    #[test]
    fn fixed_rounds_the_exact_binary_value_half_to_even() {
        assert_eq!(fixed(0.03125, 4), "0.0312");
        assert_eq!(fixed(0.09375, 4), "0.0938");
        assert_eq!(fixed(0.15625, 4), "0.1562");
        assert_eq!(fixed(30.03125, 4), "30.0312");
        assert_eq!(fixed(2.5, 0), "2");
        assert_eq!(fixed(1.5, 0), "2");
        assert_eq!(fixed(-0.0, 4), "-0.0000");
        assert_eq!(fixed(-1e-9, 4), "-0.0000");
        assert_eq!(fixed(-1e-300, 2), "-0.00");
        assert_eq!(fixed(5e-324, 4), "0.0000");
        assert_eq!(fixed(1.1177e-01, 4), "0.1118");
        assert_eq!(fixed(1.3681e-06, 4), "0.0000");
        assert_eq!(fixed(5.3355e-05, 4), "0.0001");
        assert_eq!(fixed(7_862_476.0, 4), "7862476.0000");
        assert_eq!(fixed(26_238_668.666_666_668, 4), "26238668.6667");
        assert_eq!(fixed(1e15, 4), "1000000000000000.0000");
        assert_eq!(fixed(1e20, 4), "100000000000000000000.0000");
        assert_eq!(fixed(123_456_789.000_05, 4), "123456789.0000");
        assert_eq!(fixed(0.99995, 4), "1.0000");
        assert_eq!(fixed(1.00005, 4), "1.0001");
        assert_eq!(fixed(-27.0 / 113.0, 3), "-0.239");
        assert_eq!(
            fixed(2.0f64.powi(70), 9),
            "1180591620717411303424.000000000"
        );
        assert_eq!(fixed(0.1, 10), "0.1000000000");
        // A magnitude of `2^93` or more takes the general formatter, which
        // prints the exact binary value of the double.
        assert!(fixed_by_integers(2.0f64.powi(93), 4).is_none());
        assert_eq!(fixed(2.0f64.powi(93), 4), format!("{:.4}", 2.0f64.powi(93)));
        assert_eq!(fixed(1e30, 4), "1000000000000000019884624838656.0000");
    }

    /// Every double of a deterministic sweep over the exponent range, at every
    /// decimal count the integer path takes, prints as the general formatter
    /// prints it; the sweep is the proof the two agree, since the integer path
    /// never consults the general one.
    #[test]
    fn fixed_agrees_with_the_general_formatter_across_the_exponent_range() {
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        // The exponent range reaches the integer path's zero, shift-right and
        // shift-left branches. The exponent is drawn from the bits above the
        // word's low byte rather than from the top byte alone, which carries
        // the sign and would confine each sign to part of the range.
        const LOWEST_EXPONENT: i32 = -144;
        let exponent_count = (INTEGER_PATH_MAX_EXPONENT - LOWEST_EXPONENT + 1) as u64;
        for _ in 0..40_000 {
            let word = next();
            let exponent = LOWEST_EXPONENT + ((word >> 8) % exponent_count) as i32;
            let exponent_bits = (1075 + exponent) as u64;
            let bits = (word & 0x800F_FFFF_FFFF_FFFF) | (exponent_bits << 52);
            let value = f64::from_bits(bits);
            for decimals in 1..=INTEGER_PATH_MAX_DECIMALS {
                assert_eq!(
                    fixed_by_integers(value, decimals).unwrap(),
                    format!("{value:.decimals$}"),
                    "{value:e} with {decimals} decimals"
                );
            }
        }
        for subnormal in [
            5e-324,
            -5e-324,
            2.2250738585072014e-308,
            f64::MIN_POSITIVE / 3.0,
        ] {
            assert_eq!(
                fixed_by_integers(subnormal, 4).unwrap(),
                format!("{subnormal:.4}")
            );
        }
    }

    /// Perl's string context writes a double with fifteen significant digits,
    /// the fixed form up to `1e14` and the exponent form from `1e15`.
    #[test]
    fn perl_number_writes_fifteen_significant_digits() {
        assert_eq!(perl_number(27.0 / 113.0), "0.238938053097345");
        assert_eq!(perl_number(0.1 + 0.2), "0.3");
        assert_eq!(perl_number(2.0 / 3.0), "0.666666666666667");
        assert_eq!(perl_number(15.2040980383754), "15.2040980383754");
        assert_eq!(perl_number(-1767.31707369443), "-1767.31707369443");
        assert_eq!(perl_number(0.5), "0.5");
        assert_eq!(perl_number(-0.0), "0");
        assert_eq!(perl_number(100.0), "100");
        assert_eq!(perl_number(1e-5), "1e-05");
        assert_eq!(perl_number(1e14), "100000000000000");
        assert_eq!(perl_number(1e15), "1e+15");
        assert_eq!(perl_number(1e21), "1e+21");
    }

    #[test]
    fn non_finite_values_print_as_perl_spells_them() {
        assert_eq!(sprintf_g(f64::INFINITY), "Inf");
        assert_eq!(sprintf_g(f64::NEG_INFINITY), "-Inf");
        assert_eq!(sprintf_g(f64::NAN), "NaN");
        assert_eq!(fixed(f64::NAN, 2), "NaN");
    }
}
