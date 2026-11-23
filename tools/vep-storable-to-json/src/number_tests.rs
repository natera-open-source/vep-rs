// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Expected values are perl 5.32.1's: `printf "%.15g"` for the texts, and for the
//! strings `my $v = $s + 0` with the result's integer flag read through `B`.

use super::{double_to_perl_string, numify, write_g15, Num};

/// Ten to the fifteenth, and two to the fifty-third plus one, in digits.
const TEN_TO_15: &str = concat!("1", "00000", "00000", "00000");
const TWO_TO_53_PLUS_1: &str = concat!("90071", "99254", "740993");

fn g15(value: f64) -> String {
    let mut out = Vec::new();
    write_g15(value, &mut out);
    String::from_utf8(out).expect("ASCII")
}

fn json(value: Num) -> String {
    let mut out = Vec::new();
    value.write_json(&mut out);
    String::from_utf8(out).expect("ASCII")
}

#[test]
fn a_double_is_written_as_c_prints_fifteen_significant_digits() {
    for (value, expected) in [
        (0.1, "0.1"),
        (1.5, "1.5"),
        (100.0, "100"),
        (1e15, "1e+15"),
        (1e14, "100000000000000"),
        (1.234_567_890_123_456_8e17, "1.23456789012346e+17"),
        (0.0001, "0.0001"),
        (0.00001, "1e-05"),
        (5.461_942_896_429_71, "5.46194289642971"),
        (10.678_436_016_133_5, "10.6784360161335"),
        (1.0 / 3.0, "0.333333333333333"),
        (2.0 / 3.0, "0.666666666666667"),
        (1e100, "1e+100"),
        (-2.5, "-2.5"),
        (4.1, "4.1"),
        (4.4, "4.4"),
        (0.1 + 0.2, "0.3"),
        (999_999_999_999_999.9, "1e+15"),
        (99_999_999_999_999.9, "99999999999999.9"),
        (1e21, "1e+21"),
        (5e-324, "4.94065645841247e-324"),
        (f64::MAX, "1.79769313486232e+308"),
        (123_456.789, "123456.789"),
        (0.5, "0.5"),
        (12_345_678_901_234_567_890.0, "1.23456789012346e+19"),
        (2.5e-7, "2.5e-07"),
        (1024.0, "1024"),
        (0.0, "0"),
        (-0.0, "-0"),
        (f64::INFINITY, "inf"),
        (f64::NEG_INFINITY, "-inf"),
        (f64::NAN, "nan"),
    ] {
        assert_eq!(g15(value), expected, "{value:e}");
    }
}

#[test]
fn perls_own_string_for_a_double_names_infinity_and_nan_in_capitals() {
    assert_eq!(double_to_perl_string(f64::INFINITY), b"Inf");
    assert_eq!(double_to_perl_string(f64::NEG_INFINITY), b"-Inf");
    assert_eq!(double_to_perl_string(f64::NAN), b"NaN");
    assert_eq!(double_to_perl_string(1e15), b"1e+15");
    assert_eq!(double_to_perl_string(0.25), b"0.25");
}

#[test]
fn a_string_plus_zero_is_the_integer_or_double_perl_makes() {
    // (the string, whether `$s + 0` is an integer, its text)
    for (text, integer, expected) in [
        ("12", true, "12"),
        (" 12", true, "12"),
        ("12abc", false, "12"),
        ("-7", true, "-7"),
        ("+7", true, "7"),
        ("1.5", false, "1.5"),
        (".5", false, "0.5"),
        ("5.", false, "5"),
        ("1e3", true, "1000"),
        ("1E3", true, "1000"),
        ("1e", false, "1"),
        ("abc", false, "0"),
        ("", false, "0"),
        ("0x10", false, "0"),
        ("1_000", false, "1"),
        ("18446744073709551615", true, "18446744073709551615"),
        ("18446744073709551616", false, "1.84467440737096e+19"),
        ("-9223372036854775808", true, "-9223372036854775808"),
        ("-9223372036854775809", false, "-9.22337203685478e+18"),
        ("inf", false, "inf"),
        ("-Infinity", false, "-inf"),
        ("nan", false, "nan"),
        ("0.25xyz", false, "0.25"),
        ("1.5e-3", false, "0.0015"),
        ("  \t42  ", true, "42"),
        ("-0", true, "0"),
        ("-0.0", false, "0"),
        ("45000787", true, "45000787"),
        ("5.46194289642971", false, "5.46194289642971"),
        ("4.1", false, "4.1"),
        ("1e15", true, TEN_TO_15),
        ("1e16", true, "10000000000000000"),
        (TWO_TO_53_PLUS_1, true, TWO_TO_53_PLUS_1),
    ] {
        let sum = numify(text.as_bytes()).plus(Num::Int(0));
        assert_eq!(matches!(sum, Num::Int(_)), integer, "{text:?} gave {sum:?}");
        assert_eq!(json(sum), expected, "{text:?}");
    }
}

#[test]
fn a_whole_double_from_a_string_with_a_point_stays_a_double() {
    // Ten to the fifteenth with `.0` after it, plus 0, is a double, which C prints with
    // an exponent; the same digits without the point are an integer.
    let with_point = format!("{TEN_TO_15}.0");
    assert_eq!(
        json(numify(with_point.as_bytes()).plus(Num::Int(0))),
        "1e+15"
    );
    assert_eq!(
        json(numify(TEN_TO_15.as_bytes()).plus(Num::Int(0))),
        TEN_TO_15
    );
}

#[test]
fn arithmetic_on_exact_integers_stays_an_integer_until_it_leaves_64_bits() {
    assert_eq!(Num::Int(5).minus(Num::Int(7)), Num::Int(-2));
    // A stored double that is whole enters arithmetic as an integer.
    assert_eq!(Num::Float(1000.0).plus(Num::Int(0)), Num::Int(1000));
    assert_eq!(Num::Float(1.5).plus(Num::Int(1)), Num::Float(2.5));
    assert_eq!(
        Num::Int(u64::MAX as i128).plus(Num::Int(1)),
        Num::Float(18_446_744_073_709_551_616.0)
    );
    assert_eq!(
        Num::Int(i64::MIN as i128).minus(Num::Int(1)),
        Num::Float(-9_223_372_036_854_775_809.0)
    );
    // 2^53 and above, a double is not taken for an integer.
    assert_eq!(
        Num::Float(9_007_199_254_740_992.0).plus(Num::Int(0)),
        Num::Float(9_007_199_254_740_992.0)
    );
}

#[test]
fn comparison_orders_integers_exactly_and_knows_no_order_for_nan() {
    use std::cmp::Ordering;
    assert_eq!(
        Num::Int(9_007_199_254_740_993).compare(Num::Int(9_007_199_254_740_992)),
        Some(Ordering::Greater)
    );
    assert_eq!(Num::Int(2).compare(Num::Float(2.5)), Some(Ordering::Less));
    assert_eq!(Num::Decimal(f64::NAN).compare(Num::Int(1)), None);
}
