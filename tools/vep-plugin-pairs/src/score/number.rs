// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The value rule: two values are equal when their text is equal, or when both are decimal
//! numbers that differ, in exact decimal arithmetic, by at most the tolerance of the larger
//! of the two and by at most the tolerance outright.

use std::cmp::Ordering;
use std::collections::HashMap;

/// The bytes of a value's text a decimal number has at most.
pub const DECIMAL_MAX_TEXT: usize = 64;
/// The tolerance is ten to the power of minus this.
const TOLERANCE_PLACES: i32 = 6;
/// The tolerance as the report states it.
pub const TOLERANCE: &str = "0.000001";
/// The digits of the widest whole number the short arithmetic holds with room to shift it
/// by the tolerance's places.
const SHORT_DIGITS: usize = 31;

/// A decimal number as its text spells it: a sign, the digits without leading zeros and
/// the power of ten of the last digit. Zero has no digit.
#[derive(Debug, Clone, Copy)]
pub struct Decimal {
    negative: bool,
    digits: [u8; DECIMAL_MAX_TEXT],
    count: usize,
    exponent: i32,
}

impl Decimal {
    fn digits(&self) -> &[u8] {
        &self.digits[..self.count]
    }

    fn is_zero(&self) -> bool {
        self.count == 0
    }

    /// The order of two numbers' sizes, their signs set aside.
    fn cmp_size(&self, other: &Decimal) -> Ordering {
        match (self.is_zero(), other.is_zero()) {
            (true, true) => return Ordering::Equal,
            (true, false) => return Ordering::Less,
            (false, true) => return Ordering::Greater,
            (false, false) => {}
        }
        // The power of ten above the first digit decides unless both share it.
        let (top, other_top) = (
            self.exponent + self.count as i32,
            other.exponent + other.count as i32,
        );
        top.cmp(&other_top).then_with(|| {
            let longest = self.count.max(other.count);
            let digit = |number: &Decimal, at: usize| number.digits().get(at).copied().unwrap_or(0);
            (0..longest)
                .map(|at| digit(self, at).cmp(&digit(other, at)))
                .find(|order| order.is_ne())
                .unwrap_or(Ordering::Equal)
        })
    }

    /// The order of two numbers by value: a zero of either sign equals the other, and two
    /// spellings of one number are equal.
    pub fn cmp_value(&self, other: &Decimal) -> Ordering {
        let sign = |number: &Decimal| match (number.is_zero(), number.negative) {
            (true, _) => 0,
            (false, true) => -1,
            (false, false) => 1,
        };
        sign(self).cmp(&sign(other)).then_with(|| {
            if sign(self) < 0 {
                other.cmp_size(self)
            } else {
                self.cmp_size(other)
            }
        })
    }

    /// The digits as a whole number at the power of ten `exponent`, which is at most the
    /// number's own; `None` where the short arithmetic cannot hold it.
    fn short(&self, exponent: i32) -> Option<u128> {
        let shift = usize::try_from(self.exponent - exponent).ok()?;
        if self.count + shift > SHORT_DIGITS {
            return None;
        }
        let whole = self
            .digits()
            .iter()
            .fold(0u128, |whole, &digit| whole * 10 + u128::from(digit));
        Some(whole * 10u128.pow(shift as u32))
    }

    /// The digits as a whole number at the power of ten `exponent`, which is at most the
    /// number's own, the first digit first.
    fn long(&self, exponent: i32) -> Vec<u8> {
        let shift = (self.exponent - exponent).max(0) as usize;
        let mut whole = self.digits().to_vec();
        whole.resize(self.count + shift, 0);
        whole
    }
}

/// The value as a decimal number: an optional sign, digits with an optional fraction and an
/// optional exponent of up to three digits, within `DECIMAL_MAX_TEXT` bytes; `None` for any
/// other text.
pub fn as_decimal(value: &[u8]) -> Option<Decimal> {
    if value.len() > DECIMAL_MAX_TEXT {
        return None;
    }
    let (negative, rest) = match value.split_first() {
        Some((b'-', rest)) => (true, rest),
        Some((b'+', rest)) => (false, rest),
        _ => (false, value),
    };
    let run = |text: &[u8]| text.iter().take_while(|byte| byte.is_ascii_digit()).count();
    let whole = run(rest);
    let mut at = whole;
    let mut fraction = 0;
    if rest.get(at) == Some(&b'.') {
        fraction = run(&rest[at + 1..]);
        if whole == 0 && fraction == 0 {
            return None;
        }
        at += 1 + fraction;
    } else if whole == 0 {
        return None;
    }
    let mut exponent = 0i32;
    if matches!(rest.get(at), Some(b'e' | b'E')) {
        let mut power = &rest[at + 1..];
        let lowers = match power.split_first() {
            Some((b'-', after)) => {
                power = after;
                true
            }
            Some((b'+', after)) => {
                power = after;
                false
            }
            _ => false,
        };
        if power.is_empty() || power.len() > 3 || run(power) != power.len() {
            return None;
        }
        let size = power
            .iter()
            .fold(0i32, |size, &digit| size * 10 + i32::from(digit - b'0'));
        exponent = if lowers { -size } else { size };
    } else if at != rest.len() {
        return None;
    }
    let mut number = Decimal {
        negative,
        digits: [0; DECIMAL_MAX_TEXT],
        count: 0,
        exponent: exponent - fraction as i32,
    };
    let written = rest[..whole]
        .iter()
        .chain(&rest[whole + usize::from(at > whole)..][..fraction]);
    for &digit in written {
        if number.count > 0 || digit != b'0' {
            number.digits[number.count] = digit - b'0';
            number.count += 1;
        }
    }
    Some(number)
}

fn significant(whole: &[u8]) -> &[u8] {
    let first = whole.iter().position(|&digit| digit != 0);
    first.map_or(&[], |first| &whole[first..])
}

fn cmp_whole(one: &[u8], other: &[u8]) -> Ordering {
    let (one, other) = (significant(one), significant(other));
    one.len().cmp(&other.len()).then_with(|| one.cmp(other))
}

/// The sum of two whole numbers written first digit first.
fn sum(one: &[u8], other: &[u8]) -> Vec<u8> {
    let mut total = Vec::with_capacity(one.len().max(other.len()) + 1);
    let (mut ones, mut others) = (one.iter().rev(), other.iter().rev());
    let mut carried = 0;
    loop {
        let (digit, other_digit) = (ones.next(), others.next());
        if digit.is_none() && other_digit.is_none() {
            break;
        }
        let column = digit.copied().unwrap_or(0) + other_digit.copied().unwrap_or(0) + carried;
        total.push(column % 10);
        carried = column / 10;
    }
    total.push(carried);
    total.reverse();
    total
}

/// The larger of two whole numbers less the smaller, written first digit first.
fn difference(one: &[u8], other: &[u8]) -> Vec<u8> {
    let (larger, smaller) = if cmp_whole(one, other).is_lt() {
        (other, one)
    } else {
        (one, other)
    };
    let mut apart = Vec::with_capacity(larger.len());
    let mut smaller = smaller.iter().rev();
    let mut borrowed = 0;
    for &digit in larger.iter().rev() {
        let taken = smaller.next().copied().unwrap_or(0) + borrowed;
        borrowed = u8::from(digit < taken);
        apart.push(digit + 10 * borrowed - taken);
    }
    apart.reverse();
    apart
}

/// Whether `apart`, a whole number at the power of ten `exponent`, is at most the
/// tolerance: at most one at the power of ten minus `TOLERANCE_PLACES`.
fn within_the_tolerance(apart: &[u8], exponent: i32) -> bool {
    let apart = significant(apart);
    let places = -TOLERANCE_PLACES - exponent;
    if apart.is_empty() {
        return true;
    }
    if places < 0 {
        return false;
    }
    match apart.len().cmp(&(places as usize + 1)) {
        Ordering::Less => true,
        Ordering::Greater => false,
        Ordering::Equal => apart[0] == 1 && apart[1..].iter().all(|&digit| digit == 0),
    }
}

fn agree_long(first: &Decimal, second: &Decimal, exponent: i32) -> bool {
    let (one, other) = (first.long(exponent), second.long(exponent));
    let apart = if first.negative == second.negative {
        difference(&one, &other)
    } else {
        sum(&one, &other)
    };
    if !within_the_tolerance(&apart, exponent) {
        return false;
    }
    let larger = if cmp_whole(&one, &other).is_lt() {
        &other
    } else {
        &one
    };
    let mut share = apart;
    share.resize(share.len() + TOLERANCE_PLACES as usize, 0);
    cmp_whole(&share, larger).is_le()
}

fn agree_short(one: u128, other: u128, same_sign: bool, exponent: i32) -> bool {
    let apart = if same_sign {
        one.abs_diff(other)
    } else {
        one + other
    };
    let places = -TOLERANCE_PLACES - exponent;
    let outright =
        apart == 0 || (places >= 0 && (places >= 38 || apart <= 10u128.pow(places as u32)));
    outright && apart * 10u128.pow(TOLERANCE_PLACES as u32) <= one.max(other)
}

/// Whether two decimal numbers differ by at most the tolerance of the larger of the two in
/// size and by at most the tolerance outright.
pub fn decimals_agree(first: &Decimal, second: &Decimal) -> bool {
    let exponent = first.exponent.min(second.exponent);
    match (first.short(exponent), second.short(exponent)) {
        (Some(one), Some(other)) => {
            let same_sign = first.negative == second.negative;
            agree_short(one, other, same_sign, exponent)
        }
        _ => agree_long(first, second, exponent),
    }
}

/// Whether both values are decimal numbers that agree (`decimals_agree`).
pub fn numbers_agree(one: &[u8], other: &[u8]) -> bool {
    match (as_decimal(one), as_decimal(other)) {
        (Some(first), Some(second)) => decimals_agree(&first, &second),
        _ => false,
    }
}

/// `values` in their order without the ones `paired` counts, each as often as counted.
pub fn left_over<'a>(values: Vec<&'a [u8]>, paired: &HashMap<&[u8], usize>) -> Vec<&'a [u8]> {
    if paired.is_empty() {
        return values;
    }
    let mut owed = paired.clone();
    let mut left = Vec::with_capacity(values.len());
    for value in values {
        match owed.get_mut(value) {
            Some(times) if *times > 0 => *times -= 1,
            _ => left.push(value),
        }
    }
    left
}

/// What pairing two sides' values gives: the pairs equal as text, the pairs equal as
/// numbers alone, and the values each side has left, in its order.
pub type Pairs<'a> = (usize, usize, Vec<&'a [u8]>, Vec<&'a [u8]>);

fn counted<'a>(values: &[&'a [u8]]) -> HashMap<&'a [u8], usize> {
    let mut counts = HashMap::with_capacity(values.len());
    for &value in values {
        *counts.entry(value).or_default() += 1;
    }
    counts
}

fn numbers_in_order<'a>(values: &[&'a [u8]]) -> Vec<(Decimal, &'a [u8])> {
    let mut numbers: Vec<(Decimal, &[u8])> = values
        .iter()
        .filter_map(|&value| as_decimal(value).map(|number| (number, value)))
        .collect();
    numbers.sort_by(|one, other| one.0.cmp_value(&other.0).then_with(|| one.1.cmp(other.1)));
    numbers
}

/// The product of two sides' counts of values up to which equal text is paired by scanning.
const FEW_VALUES: usize = 64;

/// The pairs equal as text of two short sides and what each side has left: each value of
/// the reference pairs with the engine's first unpaired copy of it, which leaves
/// each side without the first copies of a value that both sides hold.
fn pair_few_as_text<'a>(
    ref_values: &[&'a [u8]],
    eng_values: &[&'a [u8]],
) -> (usize, Vec<&'a [u8]>, Vec<&'a [u8]>) {
    let mut taken = [false; FEW_VALUES];
    let mut ref_left = Vec::new();
    for &value in ref_values {
        let same = (0..eng_values.len()).find(|&at| !taken[at] && eng_values[at] == value);
        match same {
            Some(at) => taken[at] = true,
            None => ref_left.push(value),
        }
    }
    let eng_left: Vec<&[u8]> = eng_values
        .iter()
        .zip(taken)
        .filter(|(_, taken)| !taken)
        .map(|(&value, _)| value)
        .collect();
    (ref_values.len() - ref_left.len(), ref_left, eng_left)
}

/// Pairs two sides' values one to one, equal text first, then as many equal numbers as the
/// two sides allow, a value as often as its side holds it.
pub fn pair_values<'a>(ref_values: &[&'a [u8]], eng_values: &[&'a [u8]]) -> Pairs<'a> {
    if ref_values.is_empty() || eng_values.is_empty() {
        return (0, 0, ref_values.to_vec(), eng_values.to_vec());
    }
    let (as_text, mut ref_left, mut eng_left) = if ref_values.len() * eng_values.len() <= FEW_VALUES
    {
        pair_few_as_text(ref_values, eng_values)
    } else {
        let eng_counts = counted(eng_values);
        let as_text: HashMap<&[u8], usize> = counted(ref_values)
            .into_iter()
            .filter_map(|(value, times)| {
                let shared = times.min(eng_counts.get(value).copied().unwrap_or(0));
                (shared > 0).then_some((value, shared))
            })
            .collect();
        (
            as_text.values().sum(),
            left_over(ref_values.to_vec(), &as_text),
            left_over(eng_values.to_vec(), &as_text),
        )
    };
    let mut as_numbers = 0;
    if !ref_left.is_empty() && !eng_left.is_empty() {
        // The numbers one number agrees with are a run of the other side's order that
        // moves up with it, so walking both orders and dropping the smaller of two that do
        // not agree pairs the most numbers the rule allows.
        let ref_numbers = numbers_in_order(&ref_left);
        let eng_numbers = numbers_in_order(&eng_left);
        let mut ref_numeric: HashMap<&[u8], usize> = HashMap::new();
        let mut eng_numeric: HashMap<&[u8], usize> = HashMap::new();
        let (mut i, mut j) = (0, 0);
        while i < ref_numbers.len() && j < eng_numbers.len() {
            let (reference, engine) = (&ref_numbers[i], &eng_numbers[j]);
            if decimals_agree(&reference.0, &engine.0) {
                *ref_numeric.entry(reference.1).or_default() += 1;
                *eng_numeric.entry(engine.1).or_default() += 1;
                as_numbers += 1;
                i += 1;
                j += 1;
            } else if reference.0.cmp_value(&engine.0).is_lt() {
                i += 1;
            } else {
                j += 1;
            }
        }
        ref_left = left_over(ref_left, &ref_numeric);
        eng_left = left_over(eng_left, &eng_numeric);
    }
    (as_text, as_numbers, ref_left, eng_left)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The rule on the long arithmetic alone, whatever the numbers' sizes.
    fn agree_on_long_arithmetic(one: &[u8], other: &[u8]) -> bool {
        match (as_decimal(one), as_decimal(other)) {
            (Some(first), Some(second)) => {
                agree_long(&first, &second, first.exponent.min(second.exponent))
            }
            _ => false,
        }
    }

    #[test]
    fn numbers_agree_within_the_tolerance_in_decimal_arithmetic() {
        let cases: &[(&[u8], &[u8], bool)] = &[
            (b"1.5", b"1.500001", true),
            (b"1.5", b"1.5000011", false),
            (b"1.5", b"1.499999", true),
            (b"1.5", b"1.4999989", false),
            (b"0.5", b"0.5000005", true),
            (b"0.5", b"0.500001", false),
            (b"0.5", b"0.4999995", true),
            (b"0.5", b"0.49999949", false),
            (b"-0.5", b"-0.5000005", true),
            (b"-0.5", b"0.5", false),
            (b"4E-24", b"4.000001E-24", true),
            (b"4E-24", b"4.00001E-24", false),
            (b"-0.000001", b"0", false),
            (b"0", b"1e-30", false),
            (b"0", b"0.0", true),
            (b"-0", b"0", true),
            (b"1000000", b"1000001", false),
            (b"1000000", b"1000000.000001", true),
            (b"1", b"2", false),
            (b"100", b"1e2", true),
            (b"1e-3", b"0.001", true),
            (b"5E-8", b"0.00000005", true),
            (b"1e999", b"1e999", true),
            (b"1e999", b"-1e999", false),
            (
                b"12345678901234567890.0000001",
                b"12345678901234567890.0000012",
                false,
            ),
            (
                b"12345678901234567890.0000001",
                b"12345678901234567890.0000011",
                true,
            ),
            (b"0.0000010000000000000000000000000001", b"0", false),
            (b"1.", b"1", true),
            (b".5", b"0.5", true),
            (b"1e1000", b"1e1000", false),
            (b"nan", b"nan", false),
            (b"inf", b"inf", false),
            (b"1_0", b"10", false),
            (b" 1", b"1", false),
            (b"0x10", b"16", false),
            (b"1,2", b"1,2.0", false),
            (b"abc", b"abc", false),
        ];
        for &(one, other, agree) in cases {
            let shown = (String::from_utf8_lossy(one), String::from_utf8_lossy(other));
            assert_eq!(numbers_agree(one, other), agree, "{shown:?}");
            assert_eq!(numbers_agree(other, one), agree, "{shown:?}");
            assert_eq!(agree_on_long_arithmetic(one, other), agree, "{shown:?}");
            assert_eq!(agree_on_long_arithmetic(other, one), agree, "{shown:?}");
        }
    }

    #[test]
    fn tolerance_is_the_published_one() {
        assert_eq!(TOLERANCE, "0.000001");
        assert_eq!(
            TOLERANCE,
            format!("0.{}1", "0".repeat(TOLERANCE_PLACES as usize - 1))
        );
    }

    #[test]
    fn numbers_close_in_size_alone_or_one_apart_are_not_equal() {
        let small_numbers_apart = [
            ("3E-24", "6E-9"),
            ("7E-9", "8E-10"),
            ("3E-11", "4E-18"),
            ("3E-18", "3E-8"),
            ("4E-7", "7E-7"),
            ("2E-16", "4E-8"),
            ("3E-9", "9E-7"),
        ];
        let whole_numbers_one_apart = [
            ("1", "2"),
            ("0", "1"),
            ("1000000", "1000001"),
            ("250", "251"),
            ("-7", "-8"),
        ];
        for (one, other) in small_numbers_apart.iter().chain(&whole_numbers_one_apart) {
            assert!(
                !numbers_agree(one.as_bytes(), other.as_bytes()),
                "{one} {other}"
            );
            assert!(
                !agree_on_long_arithmetic(one.as_bytes(), other.as_bytes()),
                "{one} {other}"
            );
        }
    }

    #[test]
    fn two_spellings_of_one_number_are_equal_as_numbers() {
        let spellings_of_one_number = [
            ("5E-9", "5e-09"),
            ("5E-9", "0.000000005"),
            ("4E-24", "4.0E-24"),
            ("1.0", "1"),
            ("100", "1e2"),
            ("0.10", "0.1"),
        ];
        for (one, other) in spellings_of_one_number {
            assert!(
                numbers_agree(one.as_bytes(), other.as_bytes()),
                "{one} {other}"
            );
            assert!(
                agree_on_long_arithmetic(one.as_bytes(), other.as_bytes()),
                "{one} {other}"
            );
            let (first, second) = (
                as_decimal(one.as_bytes()).unwrap(),
                as_decimal(other.as_bytes()).unwrap(),
            );
            assert_eq!(first.cmp_value(&second), Ordering::Equal, "{one} {other}");
        }
    }

    #[test]
    fn decimal_number_is_a_sign_digits_a_fraction_and_a_short_exponent() {
        for text in [
            "0", "-0", "+1", "1.", ".5", "-.5e3", "1e5", "1E+5", "1e-005", "001.100", "1e999",
        ] {
            assert!(as_decimal(text.as_bytes()).is_some(), "{text}");
        }
        let longest = format!("{}1", "0".repeat(DECIMAL_MAX_TEXT - 1));
        assert!(as_decimal(longest.as_bytes()).is_some());
        assert!(as_decimal(format!("0{longest}").as_bytes()).is_none());
        for text in [
            "", "+", "-", ".", "e5", "1e", "1e+", "1e1000", "1e-1000", "1.2.3", "1 ", " 1", "--1",
            "1e5.0", "0x1", "1d3", "١", "Infinity", "NaN",
        ] {
            assert!(as_decimal(text.as_bytes()).is_none(), "{text}");
        }
    }

    #[test]
    fn numbers_order_by_value_whatever_their_spelling() {
        let order = |one: &str, other: &str| {
            as_decimal(one.as_bytes())
                .unwrap()
                .cmp_value(&as_decimal(other.as_bytes()).unwrap())
        };
        for (smaller, larger) in [
            ("1", "2"),
            ("-2", "-1"),
            ("-1", "0"),
            ("-1e-999", "0"),
            ("0", "1e-999"),
            ("9", "10"),
            ("0.9", "1"),
            ("1.05", "1.5"),
            ("99e-1", "1e1"),
            ("1e10", "1.0000000001e10"),
            ("-1e999", "1e-999"),
            ("123456789012345678901234567890123456789", "1e39"),
        ] {
            assert_eq!(order(smaller, larger), Ordering::Less, "{smaller} {larger}");
            assert_eq!(
                order(larger, smaller),
                Ordering::Greater,
                "{smaller} {larger}"
            );
        }
        for (one, other) in [
            ("0", "-0"),
            ("0.0", "0e5"),
            ("1.0", "1"),
            ("10", "1e1"),
            ("-5e-1", "-0.50"),
        ] {
            assert_eq!(order(one, other), Ordering::Equal, "{one} {other}");
        }
    }

    #[test]
    fn short_and_long_arithmetic_give_one_answer() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = |bound: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % bound
        };
        let mut agreed = 0;
        for _ in 0..20_000 {
            let digits = 1 + next(12);
            let whole: u64 = (0..digits).fold(0, |whole, _| whole * 10 + next(10));
            let nudged = whole + next(4) - next(4).min(whole);
            let exponent = next(40) as i64 - 20;
            let other_exponent = exponent - next(3) as i64;
            let scale = 10u64.pow((exponent - other_exponent) as u32);
            let signs = [["", ""], ["-", "-"], ["", "-"]][next(3) as usize];
            let one = format!("{}{whole}e{exponent}", signs[0]);
            let other = format!("{}{}e{other_exponent}", signs[1], nudged * scale + next(2));
            let short = numbers_agree(one.as_bytes(), other.as_bytes());
            assert_eq!(
                short,
                agree_on_long_arithmetic(one.as_bytes(), other.as_bytes()),
                "{one} {other}"
            );
            agreed += usize::from(short);
        }
        assert!(agreed > 500, "{agreed}");
    }

    fn values<'a>(values: &[&'a str]) -> Vec<&'a [u8]> {
        values.iter().map(|value| value.as_bytes()).collect()
    }

    fn paired(reference: &[&str], engine: &[&str]) -> (usize, usize, Vec<String>, Vec<String>) {
        let (reference, engine) = (values(reference), values(engine));
        let (as_text, as_numbers, ref_left, eng_left) = pair_values(&reference, &engine);
        let shown = |left: Vec<&[u8]>| {
            left.into_iter()
                .map(|value| String::from_utf8_lossy(value).into_owned())
                .collect()
        };
        (as_text, as_numbers, shown(ref_left), shown(eng_left))
    }

    fn expected(
        as_text: usize,
        as_numbers: usize,
        ref_left: &[&str],
        eng_left: &[&str],
    ) -> (usize, usize, Vec<String>, Vec<String>) {
        let owned = |left: &[&str]| left.iter().map(|value| value.to_string()).collect();
        (as_text, as_numbers, owned(ref_left), owned(eng_left))
    }

    #[test]
    fn numbers_pair_as_many_as_the_tolerance_allows() {
        assert_eq!(
            paired(&["1.5", "1.500001"], &["1.5000005", "14.999992e-1"]),
            expected(0, 2, &[], &[])
        );
        assert_eq!(
            paired(&["0.5", "0.500001"], &["0.5000005", "4.999992e-1"]),
            expected(0, 1, &["0.500001"], &["4.999992e-1"])
        );
        assert_eq!(
            paired(&["a", "1", "3"], &["3.0000001", "b", "0.9999995", "7"]),
            expected(0, 2, &["a"], &["b", "7"])
        );
        assert_eq!(
            paired(&["1", "x"], &["1", "1.0", "x"]),
            expected(2, 0, &[], &["1.0"])
        );
        assert_eq!(
            paired(&["1", "3"], &["2.9999999"]),
            expected(0, 1, &["1"], &[])
        );
        assert_eq!(
            paired(&["2.9999999"], &["1", "3"]),
            expected(0, 1, &[], &["1"])
        );
        assert_eq!(
            paired(&["1.0", "1.0"], &["1", "1.00"]),
            expected(0, 2, &[], &[])
        );
    }

    #[test]
    fn value_held_more_than_once_pairs_as_often_as_both_sides_hold_it() {
        assert_eq!(
            paired(&["a", "a", "b"], &["a", "b", "b"]),
            expected(2, 0, &["a"], &["b"])
        );
        assert_eq!(
            paired(&["1.0", "1.0"], &["1"]),
            expected(0, 1, &["1.0"], &[])
        );
        assert_eq!(paired(&["1", "1"], &["1", "1.0"]), expected(1, 1, &[], &[]));
        assert_eq!(
            paired(&["x", "1", "x", "1"], &["1.00", "x"]),
            expected(1, 1, &["x", "1"], &[])
        );
        assert_eq!(paired(&["a", "a"], &[]), expected(0, 0, &["a", "a"], &[]));
    }

    #[test]
    fn equal_numbers_pair_in_the_order_of_their_bytes() {
        // Of two spellings of one number the one whose bytes sort first pairs first.
        assert_eq!(
            paired(&["1.0", "1"], &["1.00"]),
            expected(0, 1, &["1.0"], &[])
        );
        assert_eq!(
            paired(&["1.00"], &["1.0", "1"]),
            expected(0, 1, &[], &["1.0"])
        );
        assert_eq!(paired(&["0", "-0"], &["0.0"]), expected(0, 1, &["0"], &[]));
    }

    #[test]
    fn short_and_long_sides_pair_equal_text_alike() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = |bound: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % bound) as usize
        };
        let pool: [&[u8]; 5] = [b"a", b"b", b"", b"1", b"1.0"];
        for _ in 0..5_000 {
            let reference: Vec<&[u8]> = (0..next(9)).map(|_| pool[next(5)]).collect();
            let engine: Vec<&[u8]> = (0..next(9)).map(|_| pool[next(5)]).collect();
            if reference.is_empty() || engine.is_empty() {
                continue;
            }
            let counts = counted(&engine);
            let as_text: HashMap<&[u8], usize> = counted(&reference)
                .into_iter()
                .filter_map(|(value, times)| {
                    let shared = times.min(counts.get(value).copied().unwrap_or(0));
                    (shared > 0).then_some((value, shared))
                })
                .collect();
            let by_counting = (
                as_text.values().sum::<usize>(),
                left_over(reference.clone(), &as_text),
                left_over(engine.clone(), &as_text),
            );
            assert_eq!(pair_few_as_text(&reference, &engine), by_counting);
        }
    }

    #[test]
    fn values_left_over_keep_their_order_less_the_ones_counted() {
        let held = values(&["a", "b", "a", "c", "a"]);
        let none = HashMap::new();
        assert_eq!(left_over(held.clone(), &none), held);
        let counts: HashMap<&[u8], usize> = HashMap::from([
            (b"a".as_slice(), 2),
            (b"c".as_slice(), 1),
            (b"z".as_slice(), 4),
        ]);
        assert_eq!(left_over(held, &counts), values(&["b", "a"]));
    }
}
