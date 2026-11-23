// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `py_round_half_up` on the values a reader of a per-class CSV rounds: the string is the one
//! `Decimal(str(x)).quantize(..., rounding=ROUND_HALF_UP)` gives, an exact half included.

use super::*;

#[test]
fn a_half_in_the_seventh_decimal_rounds_up_where_fixed_formatting_rounds_down() {
    assert_eq!(py_round_half_up(0.1234565, 6), "0.123457");
    assert_eq!(format!("{:.6}", 0.1234565), "0.123456");
}

#[test]
fn an_exact_half_rounds_up_and_not_to_even() {
    assert_eq!(py_round_half_up(0.125, 2), "0.13");
    assert_ne!(py_round_half_up(0.125, 2), format!("{:.2}", 0.125));
}

#[test]
fn a_value_near_one_keeps_its_six_decimals() {
    assert_eq!(py_round_half_up(0.9999284, 6), "0.999928");
}
