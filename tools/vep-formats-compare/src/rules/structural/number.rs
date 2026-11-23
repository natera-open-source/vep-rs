// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The numbers of the structural-variant comparator: a whole number, or a float where the
//! interpreter's arithmetic gives one.
//!
//! A coordinate is a whole number of `Int`. A length or an end read from INFO is a float
//! unless it is past a double, and the sum of a whole number and a float is a float, so a
//! position beyond 2^53 loses its last digits there exactly as it does in the interpreter.

use crate::pyfmt;
use std::cmp::Ordering;

/// A whole number of the comparator. A number of 39 digits or more is held as the nearest
/// number `Int` holds.
pub type Int = i128;

/// A whole number or a float, as the interpreter's arithmetic keeps the two apart.
#[derive(Debug, Clone, Copy)]
pub enum Number {
    Int(Int),
    Float(f64),
}

/// `float(whole)`.
fn float_of(whole: Int) -> f64 {
    whole as f64
}

/// How a whole number stands to a finite float, compared by value as the interpreter compares
/// an `int` with a `float`.
fn int_to_float(whole: Int, float: f64) -> Ordering {
    const TWO_127: f64 = 170_141_183_460_469_231_731_687_303_715_884_105_728.0;
    if float >= TWO_127 {
        return Ordering::Less;
    }
    if float < -TWO_127 {
        return Ordering::Greater;
    }
    let cut = float.trunc();
    match whole.cmp(&(cut as Int)) {
        Ordering::Equal if float > cut => Ordering::Less,
        Ordering::Equal if float < cut => Ordering::Greater,
        ordering => ordering,
    }
}

impl Number {
    /// `abs(self)`.
    pub fn abs(self) -> Number {
        match self {
            Number::Int(whole) => Number::Int(whole.saturating_abs()),
            Number::Float(float) => Number::Float(float.abs()),
        }
    }

    /// `whole + self`.
    pub fn added_to(self, whole: Int) -> Number {
        match self {
            Number::Int(own) => Number::Int(whole.saturating_add(own)),
            Number::Float(own) => Number::Float(float_of(whole) + own),
        }
    }

    /// `self + 1`.
    pub fn plus_one(self) -> Number {
        match self {
            Number::Int(own) => Number::Int(own.saturating_add(1)),
            Number::Float(own) => Number::Float(own + 1.0),
        }
    }

    /// `self - whole`.
    pub fn minus(self, whole: Int) -> Number {
        match self {
            Number::Int(own) => Number::Int(own.saturating_sub(whole)),
            Number::Float(own) => Number::Float(own - float_of(whole)),
        }
    }

    /// How `whole` stands to `self`; `None` where `self` is a NaN, to which nothing compares.
    pub fn compared_with(self, whole: Int) -> Option<Ordering> {
        match self {
            Number::Int(own) => Some(whole.cmp(&own)),
            Number::Float(own) if own.is_nan() => None,
            Number::Float(own) => Some(int_to_float(whole, own)),
        }
    }

    /// `whole > self`.
    pub fn is_below(self, whole: Int) -> bool {
        self.compared_with(whole) == Some(Ordering::Greater)
    }

    /// `self > whole`.
    pub fn is_above(self, whole: Int) -> bool {
        self.compared_with(whole) == Some(Ordering::Less)
    }

    /// `int(self)`: a float cut toward zero.
    pub fn int(self) -> Int {
        match self {
            Number::Int(own) => own,
            Number::Float(own) => own.trunc() as Int,
        }
    }
}

/// `a + b`, the nearest number `Int` holds where the sum is past it.
pub fn plus(a: Int, b: Int) -> Int {
    a.saturating_add(b)
}

/// `a - b`, the nearest number `Int` holds where the difference is past it.
pub fn minus(a: Int, b: Int) -> Int {
    a.saturating_sub(b)
}

/// The whole number a decimal text spells, the nearest one `Int` holds where the text
/// spells one past it.
pub fn nearest(decimal: &str) -> Int {
    decimal.parse().unwrap_or(if decimal.starts_with('-') {
        Int::MIN
    } else {
        Int::MAX
    })
}

/// The whole number a run of digits with an optional sign spells, in the digits of any
/// script as `int()` reads them. A run of more digits than the interpreter converts ends
/// the run as the interpreter's `ValueError`.
pub fn whole_number(text: &str) -> Int {
    match pyfmt::int_checked(text) {
        Ok(decimal) => nearest(&decimal),
        Err(error) => crate::uncaught!("{error}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_whole_number_and_a_float_compare_by_value() {
        let two_53: Int = 1 << 53;
        assert_eq!(
            Number::Float(9007199254740992.0).compared_with(two_53 + 1),
            Some(Ordering::Greater)
        );
        assert_eq!(
            Number::Float(9007199254740992.0).compared_with(two_53),
            Some(Ordering::Equal)
        );
        assert_eq!(
            Number::Float(100.5).compared_with(100),
            Some(Ordering::Less)
        );
        assert_eq!(
            Number::Float(-100.5).compared_with(-100),
            Some(Ordering::Greater)
        );
        assert_eq!(
            Number::Float(1e300).compared_with(Int::MAX),
            Some(Ordering::Less)
        );
        assert_eq!(
            Number::Float(-1e300).compared_with(Int::MIN),
            Some(Ordering::Greater)
        );
        assert_eq!(Number::Float(f64::NAN).compared_with(0), None);
        assert!(Number::Float(300.7).is_below(301));
        assert!(Number::Float(300.7).is_above(300));
        assert!(!Number::Float(300.0).is_above(300));
    }

    #[test]
    fn a_sum_with_a_float_is_a_float_and_with_a_whole_number_a_whole_number() {
        let pos: Int = (1 << 62) + 1;
        assert_eq!(Number::Float(0.0).added_to(pos).int(), 1 << 62);
        assert_eq!(
            Number::Int(9_223_372_036_854_775_807).added_to(pos).int(),
            pos + 9_223_372_036_854_775_807
        );
        assert_eq!(Number::Float(300.7).added_to(30_425_649).int(), 30_425_949);
        assert_eq!(Number::Float(-2.5).abs().int(), 2);
        assert_eq!(Number::Float(-2.5).int(), -2);
    }

    #[test]
    fn a_whole_number_is_read_in_the_digits_of_any_script() {
        assert_eq!(whole_number("030425949"), 30_425_949);
        assert_eq!(whole_number("-5"), -5);
        assert_eq!(whole_number("٣٠"), 30);
        assert_eq!(
            whole_number("4611686018427387904000"),
            4_611_686_018_427_387_904_000
        );
        assert_eq!(whole_number(&"9".repeat(60)), Int::MAX);
        assert_eq!(whole_number(&format!("-{}", "9".repeat(60))), Int::MIN);
    }
}
