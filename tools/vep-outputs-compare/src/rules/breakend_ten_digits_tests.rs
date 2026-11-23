// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Every reader of a coordinate takes the digits `0` to `9` alone: a position, a span's end
//! or a mate's position written in the decimal digits of another script is no number to it.

use super::*;

/// `number` written in the decimal digits of another script than `0` to `9`.
fn other_digits(number: u64) -> String {
    let digit = |digit: char| char::from_u32(0x0660 + digit.to_digit(10).unwrap()).unwrap();
    number.to_string().chars().map(digit).collect()
}

#[test]
fn the_location_span_reader_takes_the_digits_0_to_9_alone() {
    assert_eq!(
        location_span("21:5-9"),
        Some(("21", Coordinate::Word(5), Coordinate::Word(9)))
    );
    let written = format!("21:{}-{}", other_digits(5), other_digits(9));
    assert_eq!(location_span(&written), None);
}

#[test]
fn the_reader_of_a_coordinate_that_is_no_integer_takes_the_digits_0_to_9_alone() {
    assert!(!non_integer_end("21:5-9"));
    let written = format!("21:{}-{}", other_digits(5), other_digits(9));
    assert!(non_integer_end(&written));
}

/// An end between blanks an integer's reading drops is read, and one after a blank it
/// refuses is not.
#[test]
fn the_span_size_reader_takes_the_digits_0_to_9_alone() {
    assert!(span_exceeds_max_sv_size("21:5-2147478648"));
    let written = format!("21:5-{}", other_digits(2147478648));
    assert!(!span_exceeds_max_sv_size(&written));
    assert!(span_exceeds_max_sv_size("21:5-\u{a0}2147478648\u{2003}"));
    assert!(!span_exceeds_max_sv_size("21:5-\u{1c}2147478648"));
}

#[test]
fn the_mate_position_reader_takes_the_digits_0_to_9_alone() {
    assert_eq!(
        bracket_mate("N[5:30000000["),
        Some(("5", Coordinate::Word(30000000)))
    );
    let written = format!("N[5:{}[", other_digits(30000000));
    assert_eq!(bracket_mate(&written), None);
}

#[test]
fn the_reader_of_vep_rs_s_position_takes_the_digits_0_to_9_alone() {
    assert_eq!(vep_rs_pos("30425650"), Some(30425650));
    assert_eq!(vep_rs_pos(&other_digits(30425650)), None);
}

/// A bracket breakend allele whose mate position, the buffer bound, is written in other
/// digits on a chromosome no cache names.
#[test]
fn a_mate_position_in_other_digits_at_the_buffer_bound_names_no_breakend_mate() {
    let allele = format!("N[Chr21:{}[", other_digits(2147478648));
    assert_eq!(bracket_mate(&allele), None);
    assert_eq!(reference_mate_chromosome(&allele), None);
}
