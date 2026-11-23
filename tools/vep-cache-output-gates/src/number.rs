// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The arithmetic of the gates as Python does it on what `json.load` returns: an `int` of
//! any sign (`True` and `False` among them, as 1 and 0) or a `float`. Two integers
//! subtract and divide exactly, the quotient rounded once to the nearest double; an
//! integer beside a float becomes the nearest double first; an integer is compared with a
//! float by value, with no conversion.

use crate::json::Value;

/// A number of the gates' arithmetic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Number {
    Int(i128),
    Float(f64),
}

/// The number a value is, or why it is none: the words follow the value's name in an
/// error line.
pub fn number(value: &Value) -> Result<Number, String> {
    let what = match value {
        Value::Bool(flag) => return Ok(Number::Int(i128::from(*flag))),
        Value::Int(n) => return Ok(Number::Int(*n)),
        Value::Float(x) => return Ok(Number::Float(*x)),
        Value::BigInt(_) => "an integer of more than 128 bits",
        Value::Null => "null, not a number",
        Value::Str(_) => "a string, not a number",
        Value::List(_) => "an array, not a number",
        Value::Object(_) => "an object, not a number",
    };
    Err(format!("is {what}"))
}

impl Number {
    /// `self > 0`.
    pub fn positive(self) -> bool {
        match self {
            Number::Int(n) => n > 0,
            Number::Float(x) => x > 0.0,
        }
    }

    /// `float(self)`.
    pub fn to_float(self) -> f64 {
        match self {
            Number::Int(n) => n as f64,
            Number::Float(x) => x,
        }
    }

    /// `self >= x`.
    pub fn at_least(self, x: f64) -> bool {
        match self {
            Number::Int(n) => int_at_least(n, x),
            Number::Float(own) => own >= x,
        }
    }
}

/// `value == 0`: true of a zero of either kind and of `False`, false of anything else,
/// whatever its type.
pub fn equals_zero(value: &Value) -> bool {
    match value {
        Value::Bool(flag) => !flag,
        Value::Int(n) => *n == 0,
        Value::Float(x) => *x == 0.0,
        _ => false,
    }
}

/// `n >= x` for an integer and a float, by value.
fn int_at_least(n: i128, x: f64) -> bool {
    // 2 to the 127th, one past the largest `i128`.
    const BOUND: f64 = 170_141_183_460_469_231_731_687_303_715_884_105_728.0;
    if x.is_nan() || x >= BOUND {
        return false;
    }
    if x < -BOUND {
        return true;
    }
    let floor = x.floor();
    let whole = floor as i128;
    n > whole || (n == whole && x == floor)
}

/// `abs(c - b) / b` for a `b` above zero.
pub fn relative_delta(c: Number, b: Number) -> Result<f64, String> {
    match (c, b) {
        (Number::Int(c), Number::Int(b)) => {
            let difference = c
                .checked_sub(b)
                .ok_or("differs from the baseline by more than 128 bits")?;
            Ok(divide(difference.unsigned_abs(), b.unsigned_abs()))
        }
        _ => Ok((c.to_float() - b.to_float()).abs() / b.to_float()),
    }
}

/// `n / d` of two integers, `d` above zero: the double nearest to the quotient, a tie
/// going to the even one.
fn divide(n: u128, d: u128) -> f64 {
    // Below 2 to the 53rd both are doubles as they stand, and one division rounds once.
    const EXACT: u128 = 1 << 53;
    if n < EXACT && d < EXACT {
        return n as f64 / d as f64;
    }
    // Long division, a bit of the quotient a step, until the quotient holds 64 bits or
    // the remainder is spent: `n / d` is then `(quotient + remainder / d) / 2^scale`.
    let mut quotient = n / d;
    let mut remainder = n % d;
    let mut scale = 0u64;
    while remainder != 0 && quotient < 1 << 63 {
        let carry = remainder >> 127 != 0;
        remainder <<= 1;
        quotient <<= 1;
        if carry || remainder >= d {
            remainder = remainder.wrapping_sub(d);
            quotient |= 1;
        }
        scale += 1;
    }
    // A remainder left stands for a fraction below one: as the quotient's lowest bit it
    // sits under the 53 bits a double keeps and under the bit that decides a tie, so the
    // one rounding of the conversion is the rounding of the whole quotient.
    let rounded = (quotient | u128::from(remainder != 0)) as f64;
    rounded * f64::from_bits((1023 - scale) << 52)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::load;
    use crate::pyfmt::py_float_repr;

    #[test]
    fn integer_division_matches_python() {
        // `repr(n / d)` of CPython for each pair.
        let cases: [(u128, u128, &str); 29] = [
            (1, 1000, "0.001"),
            (2, 3, "0.6666666666666666"),
            (1, 3, "0.3333333333333333"),
            (46_609, 125_572_915, "0.0003711708054240837"),
            (
                9_007_199_254_740_993,
                9_007_199_254_740_995,
                "0.9999999999999998",
            ),
            (9_007_199_254_740_993, 3, "3002399751580331.0"),
            (9_007_199_254_740_993, 1, "9007199254740992.0"),
            (9_007_199_254_740_995, 1, "9007199254740996.0"),
            (9_007_199_254_740_997, 2, "4503599627370498.0"),
            (1, 9_007_199_254_740_993, "1.1102230246251564e-16"),
            (3, 9_007_199_254_740_993, "3.330669073875469e-16"),
            (18_446_744_073_709_551_615, 3, "6.148914691236517e+18"),
            (
                18_446_744_073_709_551_615,
                18_446_744_073_709_551_614,
                "1.0",
            ),
            (
                36_893_488_147_419_103_231,
                18_446_744_073_709_551_617,
                "2.0",
            ),
            (
                170_141_183_460_469_231_731_687_303_715_884_105_727,
                3,
                "5.671372782015641e+37",
            ),
            (
                1,
                170_141_183_460_469_231_731_687_303_715_884_105_727,
                "5.877471754111438e-39",
            ),
            (
                340_282_366_920_938_463_463_374_607_431_768_211_455,
                340_282_366_920_938_463_463_374_607_431_768_211_453,
                "1.0",
            ),
            (
                340_282_366_920_938_463_463_374_607_431_768_211_455,
                7,
                "4.861176670299121e+37",
            ),
            (
                7,
                340_282_366_920_938_463_463_374_607_431_768_211_455,
                "2.0571151139390031e-38",
            ),
            (0, 9_007_199_254_740_993, "0.0"),
            (
                340_282_366_920_938_463_463_374_607_431_768_211_455,
                1,
                "3.402823669209385e+38",
            ),
            (
                170_141_183_460_469_250_621_153_235_194_464_960_512,
                1,
                "1.7014118346046923e+38",
            ),
            (
                170_141_183_460_469_250_621_153_235_194_464_960_513,
                1,
                "1.7014118346046927e+38",
            ),
            (
                170_141_183_460_469_250_621_153_235_194_464_960_512,
                3_145_728,
                "5.408642560973779e+31",
            ),
            (
                100_000_000_000_000_000_000_000_000_000_000_000_007,
                10_000_000_000_000_000_009,
                "1e+19",
            ),
            (
                3_802_951_800_684_688_204_490_109_616_129,
                2_535_301_200_456_458_802_993_406_410_752,
                "1.5",
            ),
            (18_446_744_073_709_551_617, 2048, "9007199254740992.0"),
            (18_446_744_073_709_553_665, 4096, "4503599627370497.0"),
            (23_058_430_092_136_939_521, 1024, "2.251799813685248e+16"),
        ];
        for (n, d, expected) in cases {
            assert_eq!(py_float_repr(divide(n, d)), expected, "{n} / {d}");
        }
    }

    #[test]
    fn a_delta_follows_the_kinds_of_its_two_numbers() {
        let delta = |c: Number, b: Number| relative_delta(c, b).map(py_float_repr);
        assert_eq!(
            delta(Number::Int(1001), Number::Int(1000)).unwrap(),
            "0.001"
        );
        assert_eq!(delta(Number::Int(999), Number::Int(1000)).unwrap(), "0.001");
        assert_eq!(delta(Number::Int(-5), Number::Int(2)).unwrap(), "3.5");
        assert_eq!(delta(Number::Float(1.5), Number::Int(1)).unwrap(), "0.5");
        assert_eq!(delta(Number::Int(3), Number::Float(1.5)).unwrap(), "1.0");
        assert_eq!(
            delta(Number::Float(0.1), Number::Float(0.3)).unwrap(),
            "0.6666666666666666"
        );
        assert_eq!(
            delta(Number::Float(f64::INFINITY), Number::Int(1)).unwrap(),
            "inf"
        );
        assert_eq!(
            delta(Number::Int(1), Number::Float(f64::INFINITY)).unwrap(),
            "nan"
        );
        assert!(delta(Number::Int(i128::MIN), Number::Int(1)).is_err());
        assert_eq!(
            delta(Number::Int(i128::MIN + 1), Number::Int(1)).unwrap(),
            "1.7014118346046923e+38"
        );
    }

    #[test]
    fn an_integer_is_compared_with_a_float_by_value() {
        let two_53 = 9_007_199_254_740_992i128;
        for (n, x, expected) in [
            (1, 0.999, true),
            (0, 0.999, false),
            (1, 1.0, true),
            (1, 1.0000000000000002, false),
            (-1, -0.5, false),
            (-1, -1.0, true),
            (-1, -1.5, true),
            (two_53 + 1, 9_007_199_254_740_992.0, true),
            (two_53 + 1, 9_007_199_254_740_994.0, false),
            (two_53 - 1, 9_007_199_254_740_992.0, false),
            (i128::MAX, 1.7014118346046923e38, false),
            (i128::MAX, 1.7014118346046921e38, true),
            (i128::MIN, -1.7014118346046923e38, true),
            (i128::MIN, -1.7014118346046921e38, false),
            (i128::MIN, f64::NEG_INFINITY, true),
            (i128::MAX, f64::INFINITY, false),
            (0, f64::NAN, false),
        ] {
            assert_eq!(Number::Int(n).at_least(x), expected, "{n} >= {x:?}");
        }
        assert!(Number::Float(0.999).at_least(0.999));
        assert!(!Number::Float(f64::NAN).at_least(0.0));
    }

    #[test]
    fn values_become_numbers_as_python_counts_them() {
        let read = |text: &str| number(&load(text.as_bytes()).unwrap());
        assert_eq!(read("true"), Ok(Number::Int(1)));
        assert_eq!(read("false"), Ok(Number::Int(0)));
        assert_eq!(read("-7"), Ok(Number::Int(-7)));
        assert_eq!(read("0.5"), Ok(Number::Float(0.5)));
        for text in [
            "null",
            "\"1\"",
            "[]",
            "{}",
            "170141183460469231731687303715884105728",
        ] {
            assert!(read(text).is_err(), "{text}");
        }
        assert!(Number::Int(1).positive() && !Number::Int(0).positive());
        assert!(Number::Float(0.5).positive() && !Number::Float(f64::NAN).positive());
        for (text, expected) in [
            ("0", true),
            ("0.0", true),
            ("-0.0", true),
            ("false", true),
            ("true", false),
            ("1", false),
            ("null", false),
            ("\"\"", false),
            ("[]", false),
            ("NaN", false),
        ] {
            assert_eq!(
                equals_zero(&load(text.as_bytes()).unwrap()),
                expected,
                "{text}"
            );
        }
    }
}
