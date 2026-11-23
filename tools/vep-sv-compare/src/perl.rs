// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Ensembl VEP's readings of a text: truth, a number, a POS column, the INFO column.
//!
//! An integer is held in 128 bits: a reading or a sum past that word stays at its bound.

use std::cmp::Ordering;
use std::ops::{Add, Mul, Sub};

use memchr::{memchr, memchr_iter, memmem};

use crate::pyfmt::py_strip;

/// The magnitude an infinite reading is held at: the largest signed 64-bit integer, a span
/// over every `--max_sv_size`.
pub const SATURATED_MAGNITUDE: i128 = (1 << 63) - 1;

/// 2^127, the first double past every integer of the word.
const TWO_POW_127: f64 = 170_141_183_460_469_231_731_687_303_715_884_105_728.0;

/// Whether Perl reads an INFO value as true: defined, not empty and not `0`.
pub fn perl_true(value: Option<&str>) -> bool {
    !matches!(value, None | Some("") | Some("0"))
}

/// The fields Perl's list-context `split` keeps: the trailing empty fields dropped and the
/// middle ones kept, so `G,,T,,` is `G`, the empty field, `T`.
pub fn perl_fields<S: AsRef<str>>(mut fields: Vec<S>) -> Vec<S> {
    while fields.last().is_some_and(|field| field.as_ref().is_empty()) {
        fields.pop();
    }
    fields
}

/// The value of `key=value` in an INFO column, Perl's `defined $info->{key}`: `None` for an
/// absent key or a bare flag, the last occurrence winning as in the hash Ensembl VEP's parser
/// builds (BaseVCF4.pm:581-583), the text up to a second `=` kept (`END2=30425700=5` is
/// `30425700`, `SVTYPE==DEL` the empty value).
pub fn info_written<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    if is_plain_key(key) {
        return last_item_named(info, key).flatten().map(to_second_equals);
    }
    // The name of an item ends at its first `=`: a key that holds one names no item.
    if key.contains('=') {
        return None;
    }
    let mut value = None;
    for item in info_items(info) {
        if let Some(rest) = item.strip_prefix(key) {
            match rest.strip_prefix('=') {
                Some(text) => value = Some(to_second_equals(text)),
                None if rest.is_empty() => value = None,
                None => {}
            }
        }
    }
    value
}

/// The `;`-separated items of an INFO column, an empty column being one empty item.
fn info_items(info: &str) -> impl Iterator<Item = &str> {
    let mut start = 0;
    memchr_iter(b';', info.as_bytes())
        .chain(std::iter::once(info.len()))
        .map(move |end| {
            let item = &info[start..end];
            start = end + 1;
            item
        })
}

/// The text of an item after its name and `=`, up to a second `=`.
fn to_second_equals(text: &str) -> &str {
    text.split_once('=').map_or(text, |(kept, _)| kept)
}

/// Whether a key is one [`last_item_named`] reads: not empty, not `.`, without `=` or `;`.
fn is_plain_key(key: &str) -> bool {
    !key.is_empty() && key != "." && !key.contains(['=', ';'])
}

/// The last item of an INFO column whose name is `key`, as what stands after the name:
/// `Some(text)` for `key=text`, the text running to the end of the item, and `None` for
/// the key alone. The column is searched for the key, not cut into its items: an item
/// starts at the start of the column or after a `;`, and no match that starts elsewhere
/// can cover the start of an item, a key holding no `;`.
fn last_item_named<'a>(info: &'a str, key: &str) -> Option<Option<&'a str>> {
    let bytes = info.as_bytes();
    let mut last = None;
    for at in memmem::find_iter(bytes, key.as_bytes()) {
        if at > 0 && bytes[at - 1] != b';' {
            continue;
        }
        let after = at + key.len();
        match bytes.get(after) {
            None | Some(b';') => last = Some(None),
            Some(b'=') => {
                let text = after + 1;
                let end = memchr(b';', &bytes[text..]).map_or(bytes.len(), |n| text + n);
                last = Some(Some(&info[text..end]));
            }
            Some(_) => {}
        }
    }
    last
}

/// The INFO column as the hash Ensembl VEP's parser builds (BaseVCF4.pm:581-583): one entry
/// a key, in the order the keys first appear, each with the value of its last occurrence.
/// `KEY=` and a bare `KEY` have the empty value, a value ends at a second `=`, and an empty
/// item or `.` is no entry.
pub fn info_fields(info: &str) -> Vec<(&str, &str)> {
    let mut fields: Vec<(&str, &str)> = Vec::new();
    for item in info.split(';') {
        if item.is_empty() || item == "." {
            continue;
        }
        let (key, rest) = item.split_once('=').unwrap_or((item, ""));
        let value = rest.split_once('=').map_or(rest, |(kept, _)| kept);
        match fields.iter_mut().find(|(held, _)| *held == key) {
            Some(entry) => entry.1 = value,
            None => fields.push((key, value)),
        }
    }
    fields
}

/// The value [`info_fields`] holds under a key; `None` for a key the column does not name.
pub fn info_value<'a>(fields: &[(&'a str, &'a str)], key: &str) -> Option<&'a str> {
    fields
        .iter()
        .rev()
        .find(|(held, _)| *held == key)
        .map(|(_, value)| *value)
}

/// The value the hash of [`info_fields`] holds under one key, read without building the
/// hash: `info_value(&info_fields(info), key)`.
pub fn info_field<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    if is_plain_key(key) {
        return last_item_named(info, key).map(|after| after.map_or("", to_second_equals));
    }
    if key.contains('=') {
        return None;
    }
    let mut value = None;
    for item in info_items(info) {
        if item.is_empty() || item == "." {
            continue;
        }
        if let Some(rest) = item.strip_prefix(key) {
            match rest.strip_prefix('=') {
                Some(text) => value = Some(to_second_equals(text)),
                None if rest.is_empty() => value = Some(""),
                None => {}
            }
        }
    }
    value
}

/// Perl's numeric reading of a text as the coordinate arithmetic holds it: a double, or the
/// integer an infinite reading saturates to and the sums over it. The two kinds compare by
/// value, exactly: no integer is rounded to a double to be compared.
#[derive(Clone, Copy, Debug)]
pub enum PerlNumber {
    Float(f64),
    Integer(i128),
}

impl PerlNumber {
    /// The magnitude, of the same kind.
    pub fn abs(self) -> PerlNumber {
        match self {
            PerlNumber::Float(value) => PerlNumber::Float(value.abs()),
            PerlNumber::Integer(value) => PerlNumber::Integer(value.saturating_abs()),
        }
    }

    /// The integer part, the fraction dropped toward zero; a double past the word is the
    /// word's bound.
    pub fn int(self) -> i128 {
        match self {
            PerlNumber::Float(value) => value as i128,
            PerlNumber::Integer(value) => value,
        }
    }

    /// The quotient by an integer, always a double: of a double, by the divisor's nearest
    /// double; of an integer, the double nearest the exact quotient. `None` for a zero
    /// divisor.
    pub fn div_int(self, divisor: i128) -> Option<PerlNumber> {
        if divisor == 0 {
            return None;
        }
        Some(PerlNumber::Float(match self {
            PerlNumber::Float(value) => value / divisor as f64,
            PerlNumber::Integer(value) => integer_quotient(value, divisor),
        }))
    }
}

/// The double nearest the exact quotient of two integers, ties to even; `divisor` is not 0.
fn integer_quotient(dividend: i128, divisor: i128) -> f64 {
    const EXACT_AS_DOUBLE: u128 = 1 << 53;
    const ROUNDS_ONCE: u128 = 1 << 55;
    let (mut quotient, by) = (dividend.unsigned_abs(), divisor.unsigned_abs());
    let magnitude = if quotient <= EXACT_AS_DOUBLE && by <= EXACT_AS_DOUBLE {
        quotient as f64 / by as f64
    } else if quotient == 0 {
        0.0
    } else {
        let mut remainder = quotient % by;
        quotient /= by;
        let mut fraction_bits = 0;
        while quotient < ROUNDS_ONCE {
            remainder <<= 1;
            quotient <<= 1;
            if remainder >= by {
                remainder -= by;
                quotient |= 1;
            }
            fraction_bits += 1;
        }
        // The lowest bit lies below the rounding position and stands for every bit of the
        // quotient not written out.
        if remainder != 0 {
            quotient |= 1;
        }
        quotient as f64 * 2f64.powi(-fraction_bits)
    };
    if (dividend < 0) == (divisor < 0) {
        magnitude
    } else {
        -magnitude
    }
}

/// How a double orders against an integer by exact value; `None` for NaN.
fn float_cmp_integer(float: f64, integer: i128) -> Option<Ordering> {
    if float.is_nan() {
        return None;
    }
    if float >= TWO_POW_127 {
        return Some(Ordering::Greater);
    }
    if float < -TWO_POW_127 {
        return Some(Ordering::Less);
    }
    let whole = float.trunc();
    let by_whole = (whole as i128).cmp(&integer);
    Some(by_whole.then_with(|| (float - whole).partial_cmp(&0.0).expect("finite")))
}

impl PartialOrd for PerlNumber {
    fn partial_cmp(&self, other: &PerlNumber) -> Option<Ordering> {
        use PerlNumber::{Float, Integer};
        match (*self, *other) {
            (Float(a), Float(b)) => a.partial_cmp(&b),
            (Integer(a), Integer(b)) => Some(a.cmp(&b)),
            (Float(a), Integer(b)) => float_cmp_integer(a, b),
            (Integer(a), Float(b)) => float_cmp_integer(b, a).map(Ordering::reverse),
        }
    }
}

impl PartialEq for PerlNumber {
    fn eq(&self, other: &PerlNumber) -> bool {
        self.partial_cmp(other) == Some(Ordering::Equal)
    }
}

impl PartialEq<i128> for PerlNumber {
    fn eq(&self, other: &i128) -> bool {
        *self == PerlNumber::Integer(*other)
    }
}

impl PartialOrd<i128> for PerlNumber {
    fn partial_cmp(&self, other: &i128) -> Option<Ordering> {
        self.partial_cmp(&PerlNumber::Integer(*other))
    }
}

/// The sum with an integer: a double takes the integer as its nearest double, an integer
/// stays one.
impl Add<i128> for PerlNumber {
    type Output = PerlNumber;

    fn add(self, other: i128) -> PerlNumber {
        match self {
            PerlNumber::Float(value) => PerlNumber::Float(value + other as f64),
            PerlNumber::Integer(value) => PerlNumber::Integer(value.saturating_add(other)),
        }
    }
}

/// The difference from an integer, kinds as for the sum.
impl Sub<i128> for PerlNumber {
    type Output = PerlNumber;

    fn sub(self, other: i128) -> PerlNumber {
        match self {
            PerlNumber::Float(value) => PerlNumber::Float(value - other as f64),
            PerlNumber::Integer(value) => PerlNumber::Integer(value.saturating_sub(other)),
        }
    }
}

/// The product with an integer, kinds as for the sum.
impl Mul<i128> for PerlNumber {
    type Output = PerlNumber;

    fn mul(self, other: i128) -> PerlNumber {
        match self {
            PerlNumber::Float(value) => PerlNumber::Float(value * other as f64),
            PerlNumber::Integer(value) => PerlNumber::Integer(value.saturating_mul(other)),
        }
    }
}

fn is_ascii_blank(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\r' | '\x0b' | '\x0c')
}

/// The length of the decimal a text opens with: the digits `0` to `9` with an optional
/// fraction, or a fraction alone, then an exponent when it is complete; `None` when the text
/// opens with neither.
fn decimal_len(text: &[u8]) -> Option<usize> {
    let digits = |from: usize| {
        text.get(from..).map_or(0, |rest| {
            rest.iter().take_while(|b| b.is_ascii_digit()).count()
        })
    };
    let whole = digits(0);
    let mut len = if whole > 0 {
        match text.get(whole) {
            Some(b'.') => whole + 1 + digits(whole + 1),
            _ => whole,
        }
    } else if text.first() == Some(&b'.') && digits(1) > 0 {
        1 + digits(1)
    } else {
        return None;
    };
    if matches!(text.get(len), Some(b'e' | b'E')) {
        let signed = usize::from(matches!(text.get(len + 1), Some(b'+' | b'-')));
        let exponent = digits(len + 1 + signed);
        if exponent > 0 {
            len += 1 + signed + exponent;
        }
    }
    Some(len)
}

/// The number Perl reads off the front of a text, as `(number, rest)`: blanks skipped, then
/// an optional sign and a decimal; `None` when the text opens with no number.
fn leading_number(text: &str) -> Option<(&str, &str)> {
    let body = text.trim_start_matches(is_ascii_blank);
    let sign = usize::from(body.starts_with(['+', '-']));
    let len = sign + decimal_len(&body.as_bytes()[sign..])?;
    Some(body.split_at(len))
}

/// Whether the whole text is a number Perl reads without a remainder: blanks in front, an
/// optional sign, the digits `0` to `9` with an optional fraction and exponent.
pub fn is_perl_number(text: &str) -> bool {
    matches!(leading_number(text), Some((_, "")))
}

/// Whether the text is one of `words` under an optional sign with blanks of any script
/// around, letters in either case. A dotted capital and a dotless small `i` stand for `i`,
/// as the case-insensitive match this reading is defined by takes them.
fn is_spelling_of(text: &str, words: &[&str]) -> bool {
    let core = py_strip(text);
    let core = core.strip_prefix(['+', '-']).unwrap_or(core);
    let letter = |c: char| match c {
        '\u{130}' | '\u{131}' => 'i',
        other => other.to_ascii_lowercase(),
    };
    words
        .iter()
        .any(|word| core.chars().map(letter).eq(word.chars()))
}

/// Whether the text spells Perl's NaN: `nan` in any case, an optional sign, blanks around.
pub fn is_perl_nan(text: &str) -> bool {
    is_spelling_of(text, &["nan"])
}

/// Whether the text spells Perl's infinity: `inf` or `infinity` in any case, an optional
/// sign, blanks around.
pub fn is_perl_infinity(text: &str) -> bool {
    is_spelling_of(text, &["inf", "infinity"])
}

/// Perl's numeric reading of a text, as `abs($svlen)` and `max(@svlen)` apply it: leading
/// blanks, an optional sign, digits with an optional fraction and exponent; anything else,
/// `.` and the empty text included, is 0 (`300abc` is 300). A reading beyond a double
/// (`1e400`, Perl's `Inf`) is the integer [`SATURATED_MAGNITUDE`] with its sign, so it stays
/// a number the coordinate arithmetic can hold. `NaN`, `Inf` and `Infinity`, Perl's NaN and
/// infinity, read as 0 here, the length vep-rs gives the spellings, so a Location derived
/// from them is the one vep-rs prints; `perl_skip_route` reads the spellings themselves.
pub fn perl_numeric_value(text: &str) -> PerlNumber {
    let Some((number, _)) = leading_number(text) else {
        return PerlNumber::Float(0.0);
    };
    let value: f64 = number.parse().expect("a decimal literal");
    if value == f64::INFINITY {
        PerlNumber::Integer(SATURATED_MAGNITUDE)
    } else if value == f64::NEG_INFINITY {
        PerlNumber::Integer(-SATURATED_MAGNITUDE)
    } else {
        PerlNumber::Float(value)
    }
}

/// The POS column as Ensembl VEP numifies it (BaseVCF4.pm:263 `$start++` on the structural
/// path): a decimal integer exactly (`4611686018427387904`), any other spelling by
/// [`perl_numeric_value`] with its fraction dropped (`3e7` is 30000000, `30425949.0` and
/// ` 30425949` are 30425949, `.`, `0x1D04A1D` and fullwidth digits are 0, `-5` is -5).
pub fn perl_pos(text: &str) -> i128 {
    let stripped = py_strip(text);
    let digits = stripped.strip_prefix(['+', '-']).unwrap_or(stripped);
    if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
        let bound = if stripped.starts_with('-') {
            i128::MIN
        } else {
            i128::MAX
        };
        return stripped.parse().unwrap_or(bound);
    }
    perl_numeric_value(text).int()
}

/// Whether a POS column denotes a non-negative number (`30425949`, `+30425949`, `030425949`,
/// `3e7`, `30425949.0`, ` 30425949`); `.`, the empty text, `0x1D04A1D`, fullwidth digits,
/// `30_425_949`, `-5` and `NaN` do not, and Ensembl VEP numifies them ([`perl_pos`]) where
/// vep-rs refuses the line. `None`, the decimal spelling of a parsed integer, is a number.
pub fn perl_pos_is_a_number(text: Option<&str>) -> bool {
    let Some(text) = text else {
        return true;
    };
    let core = text.trim_matches(is_ascii_blank);
    let core = core.strip_prefix('+').unwrap_or(core);
    decimal_len(core.as_bytes()) == Some(core.len())
}

/// Whether Perl reads the text as a number that is not an integer: a fraction (`9.5`,
/// `300.7`), a value past a double (`1e400`), `NaN` or an infinity (`Inf`, `Infinity`), with
/// an optional sign and blanks around; `9`, `9.0`, `3e2` are integers, and `abc`, `.`, the
/// empty text and a fullwidth digit are 0.
pub fn perl_number_is_not_an_integer(text: &str) -> bool {
    let stripped = py_strip(text);
    if is_perl_nan(stripped) || is_perl_infinity(stripped) {
        return true;
    }
    let Some((number, "")) = leading_number(stripped) else {
        return false;
    };
    let value: f64 = number.parse().expect("a decimal literal");
    !(value.is_finite() && value.fract() == 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TWO_POW_63: i128 = 1 << 63;
    const FULLWIDTH_POS: &str = "\u{ff13}\u{ff10}\u{ff14}\u{ff12}\u{ff15}\u{ff19}\u{ff14}\u{ff19}";

    fn float(value: f64) -> PerlNumber {
        PerlNumber::Float(value)
    }

    /// Whether the reading is a double of exactly this value, the sign of zero included.
    fn is_float(number: PerlNumber, value: f64) -> bool {
        matches!(number, PerlNumber::Float(read) if read.to_bits() == value.to_bits())
    }

    #[test]
    fn truth_is_defined_not_empty_and_not_zero() {
        assert!(perl_true(Some("DEL")));
        assert!(perl_true(Some("00")));
        assert!(perl_true(Some(" ")));
        assert!(!perl_true(Some("0")));
        assert!(!perl_true(Some("")));
        assert!(!perl_true(None));
    }

    #[test]
    fn trailing_empty_fields_are_no_fields() {
        assert_eq!(
            perl_fields("G,,T,,".split(',').collect()),
            vec!["G", "", "T"]
        );
        assert_eq!(perl_fields("".split(',').collect()), Vec::<&str>::new());
        assert_eq!(
            perl_fields(vec![String::from("<DEL>"), String::new()]),
            vec![String::from("<DEL>")]
        );
    }

    #[test]
    fn an_info_value_is_read_to_a_second_equals_sign() {
        assert_eq!(
            info_fields("SVTYPE==DEL;END=30425949"),
            vec![("SVTYPE", ""), ("END", "30425949")]
        );
        let fields = info_fields("SVTYPE=BND;CHR2=21;END2=30425700=5");
        assert_eq!(info_value(&fields, "END2"), Some("30425700"));
        assert_eq!(
            info_fields("END2=30425700=5=7;CHR2"),
            vec![("END2", "30425700"), ("CHR2", "")]
        );
        assert_eq!(info_written("SVTYPE==DEL;END=30425949", "SVTYPE"), Some(""));
        assert_eq!(
            info_written("SVTYPE=BND;CHR2=21;END2=30425700=5", "END2"),
            Some("30425700")
        );
    }

    #[test]
    fn the_last_occurrence_of_a_key_holds_its_first_place() {
        let fields = info_fields("END=5;SVTYPE=DEL;;.;END=7;FLAG;END=9");
        assert_eq!(fields, vec![("END", "9"), ("SVTYPE", "DEL"), ("FLAG", "")]);
        assert_eq!(info_value(&fields, "END"), Some("9"));
        assert_eq!(info_value(&fields, "FLAG"), Some(""));
        assert_eq!(info_value(&fields, "SVLEN"), None);
        assert_eq!(info_fields("."), Vec::new());
        assert_eq!(info_fields(""), Vec::new());
        assert_eq!(info_written("END=5;END", "END"), None);
        assert_eq!(info_written("END;END=5", "END"), Some("5"));
        assert_eq!(info_written("ENDX=5", "END"), None);
        assert_eq!(info_written("END=", "END"), Some(""));
    }

    /// Every text of up to `longest` characters over `alphabet`.
    fn texts(alphabet: &[char], longest: usize) -> Vec<String> {
        let mut all = vec![String::new()];
        let mut last = all.clone();
        for _ in 0..longest {
            last = last
                .iter()
                .flat_map(|text| alphabet.iter().map(move |c| format!("{text}{c}")))
                .collect();
            all.extend(last.iter().cloned());
        }
        all
    }

    /// [`info_written`] read item by item through `str::split`.
    fn written_item_by_item<'a>(info: &'a str, key: &str) -> Option<&'a str> {
        let mut value = None;
        for item in info.split(';') {
            match item.split_once('=') {
                Some((name, text)) if name == key => {
                    value = Some(text.split_once('=').map_or(text, |(kept, _)| kept));
                }
                None if item == key => value = None,
                _ => {}
            }
        }
        value
    }

    #[test]
    fn a_key_is_read_alike_item_by_item_and_through_the_hash() {
        let keys = [
            "A", "B", "AB", "AA", "ABA", "A.", "", ".", "A=", "=", "A;", ";", "\u{e9}",
        ];
        for info in texts(&['A', 'B', '=', ';', '.', '\u{e9}'], 6) {
            for key in keys {
                assert_eq!(
                    info_written(&info, key),
                    written_item_by_item(&info, key),
                    "{key:?} written in {info:?}"
                );
                assert_eq!(
                    info_field(&info, key),
                    info_value(&info_fields(&info), key),
                    "{key:?} held in {info:?}"
                );
            }
        }
    }

    #[test]
    fn perl_numification() {
        let cases: [(&str, f64); 8] = [
            ("300", 300.0),
            ("-300", -300.0),
            ("300abc", 300.0),
            (" 3e2", 300.0),
            ("abc", 0.0),
            ("", 0.0),
            (".", 0.0),
            ("300.7", 300.7),
        ];
        for (text, value) in cases {
            assert!(is_float(perl_numeric_value(text), value), "{text:?}");
        }
        assert!(matches!(
            perl_numeric_value("1e400"),
            PerlNumber::Integer(SATURATED_MAGNITUDE)
        ));
        assert_eq!(perl_numeric_value("1e400"), TWO_POW_63 - 1);
        assert_eq!(perl_numeric_value("-1e400"), -(TWO_POW_63 - 1));
        assert!(is_float(perl_numeric_value("1.7e308"), 1.7e308));
    }

    #[test]
    fn a_number_is_read_off_the_front_in_the_digits_0_to_9() {
        let cases: [(&str, f64); 25] = [
            ("+3", 3.0),
            ("3.", 3.0),
            (".5", 0.5),
            ("5.e3", 5000.0),
            (".5E-1", 0.05),
            ("3e", 3.0),
            ("3e+", 3.0),
            ("3e5.2", 300000.0),
            ("1_000", 1.0),
            ("\t\n\x0b\x0c\r 7", 7.0),
            ("\u{a0}7", 0.0),
            ("\u{ff14}", 0.0),
            ("7\u{ff14}", 7.0),
            ("+-3", 0.0),
            ("- 3", 0.0),
            ("0x1D", 0.0),
            ("NaN", 0.0),
            ("Inf", 0.0),
            ("-Infinity", 0.0),
            ("e5", 0.0),
            ("-.5", -0.5),
            ("-0", -0.0),
            ("\x1c7", 0.0),
            ("1e-400", 0.0),
            ("0.1e1", 1.0),
        ];
        for (text, value) in cases {
            assert!(is_float(perl_numeric_value(text), value), "{text:?}");
        }
        assert!(is_perl_number("300"));
        assert!(is_perl_number(" -300.5e+2"));
        assert!(!is_perl_number("300 "));
        assert!(!is_perl_number("300abc"));
        assert!(!is_perl_number(""));
        assert!(!is_perl_number("."));
        assert!(!is_perl_number("\u{ff14}"));
    }

    #[test]
    fn sums_keep_an_integer_exact_and_a_double_a_double() {
        let saturated = perl_numeric_value("1e400").abs() + 1000;
        assert!(matches!(saturated, PerlNumber::Integer(_)));
        assert_eq!(saturated, 1000 + TWO_POW_63 - 1);
        assert_eq!(saturated.int(), 1000 + TWO_POW_63 - 1);
        let fractional = perl_numeric_value("-300.7").abs() + 30_425_650;
        assert!(is_float(fractional, 30_425_950.7));
        assert_eq!(fractional.int(), 30_425_950);
        assert_eq!(perl_numeric_value("-300.7").int(), -300);
        assert!(is_float(perl_numeric_value("2000.5") - 1001, 999.5));
        assert!(perl_numeric_value("2000.5") - 1001 > 999);
        assert_eq!(perl_numeric_value("-1e400").abs(), TWO_POW_63 - 1);
        assert_eq!(
            (perl_numeric_value("1e400") - 1001).int(),
            TWO_POW_63 - 1002
        );
        assert!(is_float(perl_numeric_value("4.5") * 3, 13.5));
        assert_eq!(
            (perl_numeric_value("1e400") * 3).int(),
            3 * (TWO_POW_63 - 1)
        );
        assert_eq!(perl_numeric_value("1e300").int(), i128::MAX);
        assert_eq!(perl_numeric_value("-1e300").int(), i128::MIN);
        assert_eq!(
            PerlNumber::Integer(i128::MAX) + 1,
            PerlNumber::Integer(i128::MAX)
        );
    }

    #[test]
    fn an_integer_and_a_double_compare_by_exact_value() {
        let two_pow_63 = float(9_223_372_036_854_775_808.0);
        assert!(PerlNumber::Integer(TWO_POW_63 - 1) < two_pow_63);
        assert!(two_pow_63 > TWO_POW_63 - 1);
        assert_eq!(two_pow_63, TWO_POW_63);
        assert!(two_pow_63 < TWO_POW_63 + 1);
        assert!(float(0.5) > 0 && float(0.5) < 1);
        assert!(float(-0.5) < 0 && float(-0.5) > -1);
        assert_eq!(float(-0.0), 0);
        assert!(float(1e300) > i128::MAX);
        assert!(float(-1e300) < i128::MIN);
        assert_eq!(float(-TWO_POW_127), i128::MIN);
        assert!(float(TWO_POW_127) > i128::MAX);
        assert_eq!(float(f64::NAN).partial_cmp(&0), None);
        assert_ne!(float(f64::NAN), float(f64::NAN));
        assert_eq!(float(300.0), PerlNumber::Integer(300));
        assert!(PerlNumber::Integer(301) > float(300.5));
    }

    #[test]
    fn a_quotient_is_the_double_nearest_the_exact_one() {
        let quotient = |number: PerlNumber, by: i128| match number.div_int(by) {
            Some(PerlNumber::Float(value)) => value,
            other => panic!("{other:?}"),
        };
        assert_eq!(quotient(float(9.0), 3), 3.0);
        assert_eq!(quotient(float(10.0), 4), 2.5);
        assert_eq!(quotient(PerlNumber::Integer(7), 2), 3.5);
        assert_eq!(quotient(PerlNumber::Integer(-7), 2), -3.5);
        assert_eq!(quotient(PerlNumber::Integer(1), 3), 1.0 / 3.0);
        assert!(quotient(PerlNumber::Integer(0), -5).is_sign_negative());
        assert_eq!(float(1.0).div_int(0), None);
        assert_eq!(PerlNumber::Integer(1).div_int(0), None);
        let saturated = PerlNumber::Integer(SATURATED_MAGNITUDE);
        let expected: [(i128, f64); 10] = [
            (1, 9.223372036854776e18),
            (2, 4.611686018427388e18),
            (3, 3.0744573456182584e18),
            (7, 1.3176245766935393e18),
            (1023, 9016003946094600.0),
            (1025, 8998411743272952.0),
            (4097, 2251250192056328.0),
            (-3, -3.0744573456182584e18),
            (12_345_677, 747_093_256_761.437_7),
            (99_999_999_977, 92233720.38976151),
        ];
        for (by, value) in expected {
            assert_eq!(quotient(saturated, by), value, "{by}");
        }
        let two_pow_100: i128 = 1 << 100;
        assert_eq!(quotient(PerlNumber::Integer(i128::MAX), i128::MAX - 1), 1.0);
        assert_eq!(
            quotient(PerlNumber::Integer(1), two_pow_100),
            2f64.powi(-100)
        );
        assert_eq!(
            quotient(PerlNumber::Integer(two_pow_100 + 1), two_pow_100),
            1.0
        );
        // A tie goes to the even double; one more unit below the rounding position breaks it.
        assert_eq!(
            quotient(PerlNumber::Integer(two_pow_100 + (1 << 47)), two_pow_100),
            1.0
        );
        assert_eq!(
            quotient(
                PerlNumber::Integer(two_pow_100 + (1 << 47) + 1),
                two_pow_100
            ),
            1.0000000000000002
        );
        assert_eq!(
            quotient(PerlNumber::Integer(3i128.pow(70)), 7i128.pow(20)),
            3.1370903737172776e16
        );
    }

    #[test]
    fn test_a_pos_in_another_spelling_is_read_as_its_number() {
        let deletion = "SVTYPE=DEL;END=30425949";
        let records = crate::record::testing::vcf(&[
            ["21", "3e7", "exp", "N", "<DEL>", ".", ".", deletion],
            ["21", "30425949.0", "float0", "A", "G", ".", ".", "."],
            ["21", " 30425949", "lead", "A", "G", ".", ".", "."],
            ["21", "30425949 ", "trail", "A", "G", ".", ".", "."],
            ["21", "+30425949", "plus", "N", "<DEL>", ".", ".", deletion],
            [
                "21",
                "030425949",
                "zeropad",
                "N",
                "<DEL>",
                ".",
                ".",
                deletion,
            ],
            ["21", "abc", "abc", "A", "G", ".", ".", "."],
            ["21", ".", "dot", "N", "<DEL>", ".", ".", deletion],
            [
                "21",
                FULLWIDTH_POS,
                "fullwidth",
                "N",
                "<DEL>",
                ".",
                ".",
                deletion,
            ],
            [
                "21",
                "18446744073709551616",
                "huge",
                "N",
                "<DEL>",
                ".",
                ".",
                "SVTYPE=DEL;SVLEN=-300",
            ],
        ]);
        assert_eq!(
            records.iter().map(|v| v.pos).collect::<Vec<_>>(),
            [
                30_000_000,
                30_425_949,
                30_425_949,
                30_425_949,
                30_425_949,
                30_425_949,
                0,
                0,
                0,
                1 << 64
            ]
        );
        assert_eq!(
            records[..3]
                .iter()
                .map(|v| v.pos_text.as_deref())
                .collect::<Vec<_>>(),
            [Some("3e7"), Some("30425949.0"), Some(" 30425949")]
        );
        assert_eq!(
            (
                records[0].location.as_str(),
                records[0].variant_type.as_str()
            ),
            ("21:30000000", "Symbolic_DEL")
        );
        let spellings: [(&str, i128); 10] = [
            ("3e7", 30_000_000),
            ("30425949.", 30_425_949),
            ("3.0425949E7", 30_425_949),
            ("1e3", 1000),
            ("-5", -5),
            ("30_425_949", 30),
            ("0x1D04A1D", 0),
            ("NaN", 0),
            ("1e400", TWO_POW_63 - 1),
            ("4611686018427387904", 1 << 62),
        ];
        for (text, value) in spellings {
            assert_eq!(perl_pos(text), value, "{text:?}");
        }
        for text in [
            "30425949",
            "+30425949",
            "030425949",
            "3e7",
            "30425949.0",
            "30425949.",
            " 30425949",
            "30425949 ",
            "3.0425949E7",
            "1e3",
        ] {
            assert!(perl_pos_is_a_number(Some(text)), "{text:?}");
        }
        assert!(perl_pos_is_a_number(None));
        for text in [
            ".",
            "",
            "0x1D04A1D",
            FULLWIDTH_POS,
            "30_425_949",
            "-5",
            "-2147483648",
            "NaN",
            "Inf",
            "e5",
            "30425949,30425950",
        ] {
            assert!(!perl_pos_is_a_number(Some(text)), "{text:?}");
        }
    }

    #[test]
    fn a_pos_keeps_its_sign_its_blanks_and_the_word() {
        assert_eq!(perl_pos("\u{a0}30425949\u{2003}"), 30_425_949);
        assert_eq!(perl_pos("\u{a0}3e7"), 0);
        assert_eq!(perl_pos("-0"), 0);
        assert_eq!(perl_pos("+"), 0);
        assert_eq!(perl_pos("-30425949.9"), -30_425_949);
        assert_eq!(perl_pos(&"9".repeat(60)), i128::MAX);
        assert_eq!(perl_pos(&format!("-{}", "9".repeat(60))), i128::MIN);
        assert_eq!(perl_pos("1e300"), i128::MAX);
        assert!(!perl_pos_is_a_number(Some("\u{a0}30425949")));
        assert!(perl_pos_is_a_number(Some("\x0b.5e+3\x0c")));
        assert!(perl_pos_is_a_number(Some("5.")));
        assert!(perl_pos_is_a_number(Some(" +5 ")));
        assert!(!perl_pos_is_a_number(Some("+")));
        assert!(!perl_pos_is_a_number(Some("++5")));
        assert!(!perl_pos_is_a_number(Some("5 5")));
        assert!(!perl_pos_is_a_number(Some("5e")));
        assert!(!perl_pos_is_a_number(Some("\x1c5")));
    }

    #[test]
    fn a_number_that_is_not_an_integer() {
        for text in ["9.5", "300.7", "-8.5", "NaN", "-Inf", "Infinity", "1e400"] {
            assert!(perl_number_is_not_an_integer(text), "{text:?}");
        }
        for text in ["9", "9.0", "3e2", "+30425949", "abc", ".", "", "\u{ff14}"] {
            assert!(!perl_number_is_not_an_integer(text), "{text:?}");
        }
        for text in [
            " 9.5 ",
            "\u{a0}-nan\u{3000}",
            "+INFINITY",
            "-1e400",
            "1e-5",
            "\u{130}nf",
        ] {
            assert!(perl_number_is_not_an_integer(text), "{text:?}");
        }
        for text in [
            "9.5x", "infinit", "nan0", "+-nan", "9 .5", "1e5", "1e-400", " 9",
        ] {
            assert!(!perl_number_is_not_an_integer(text), "{text:?}");
        }
    }

    #[test]
    fn nan_and_infinity_are_spelled_in_any_case() {
        for text in ["NaN", "nan", " +NAN\t", "-nan", "\u{1f}nAn\u{85}"] {
            assert!(is_perl_nan(text), "{text:?}");
        }
        for text in ["", "na", "nans", "n an", "+ nan", "--nan"] {
            assert!(!is_perl_nan(text), "{text:?}");
        }
        for text in [
            "Inf",
            "-Infinity",
            " +INF ",
            "iNfInItY",
            "\u{130}nf",
            "\u{131}nf\u{131}n\u{131}ty",
        ] {
            assert!(is_perl_infinity(text), "{text:?}");
        }
        for text in ["", "in", "infinit", "infinity1", "\u{ff49}nf", "1e400"] {
            assert!(!is_perl_infinity(text), "{text:?}");
        }
    }
}
