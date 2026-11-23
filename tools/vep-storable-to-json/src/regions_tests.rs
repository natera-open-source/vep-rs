// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Expected values are perl 5.34.0's for the pattern and the `print` the comments quote.

use super::{matches_length, SequenceRegion, DECIMAL_DIGITS};
use crate::json::PerlStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Answer {
    Match,
    No,
    Dies,
}
use Answer::{Dies, Match, No};

/// The bytes a string of hexadecimal digit pairs names.
fn bytes(hex: &str) -> Vec<u8> {
    (0..hex.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&hex[at..at + 2], 16).expect("two hexadecimal digits"))
        .collect()
}

fn answer(hex: &str, utf8: bool) -> Answer {
    let bytes = bytes(hex);
    match matches_length(PerlStr {
        bytes: &bytes,
        utf8,
    }) {
        Ok(true) => Match,
        Ok(false) => No,
        Err(die) => {
            assert_eq!(die.0, b"Malformed UTF-8 character (fatal)", "{hex}");
            Dies
        }
    }
}

// eval { $length =~ /^[1-9]\d*$/ } on a string without the UTF-8 flag.
#[test]
fn a_byte_string_is_a_length_when_it_is_ascii_digits_without_a_leading_zero() {
    for (hex, expected) in [
        ("35", Match),
        ("35303030", Match),
        ("353030300a", Match),
        ("353030300a0a", No),
        ("30", No),
        ("3035", No),
        ("", No),
        ("0a", No),
        ("3541", No),
        ("3520", No),
        ("2035", No),
        ("2d35", No),
        ("2b35", No),
        ("316533", No),
        ("352e30", No),
        ("35b2", No),
        ("3580", No),
        ("35d9a1", No),
        ("8035", No),
        (
            "393939393939393939393939393939393939393939393939393939393939",
            Match,
        ),
        ("350a80", No),
        ("35ff0a", No),
    ] {
        assert_eq!(answer(hex, false), expected, "{hex}");
    }
}

// The same match on the same bytes under the UTF-8 flag (Encode::_utf8_on).
#[test]
fn a_flagged_string_is_matched_by_character_and_dies_as_perl_does_on_damage() {
    for (hex, expected) in [
        ("35", Match),
        ("35303030", Match),
        ("353030300a", Match),
        ("353030300a0a", No),
        ("370a", Match),
        ("30", No),
        ("3035", No),
        ("", No),
        ("0a", No),
        ("0a0a", No),
        ("41", No),
        ("3541", No),
        ("3520", No),
        ("2035", No),
        ("2d35", No),
        ("316533", No),
        ("35d9a1", Match),
        ("35d9a10a", Match),
        ("35e0a5a6", Match),
        ("35efbc90", Match),
        ("35f09d9f8e", Match),
        ("35f09d9fbf", Match),
        ("35f09fafb9", Match),
        ("35d9a0d9a9", Match),
        ("35d99f", No),
        ("35d9aa", No),
        ("d9a135", No),
        ("35c2b2", No),
        ("35dfbf", No),
        ("35e0a080", No),
        ("35eda080", No),
        ("35efbfbf", No),
        ("35f4908080", No),
        ("35f7bfbfbf", No),
        ("35f888808080", No),
        ("35fc8480808080", No),
        ("35fe828080808080", No),
        ("3580", Dies),
        ("35bf", Dies),
        ("358037", No),
        ("35c0b5", Dies),
        ("35c1bf", Dies),
        ("35c0", Dies),
        ("35c2", No),
        ("35c3", No),
        ("35c241", No),
        ("35c237", No),
        ("35d9", No),
        ("35d941", No),
        ("35d937", No),
        ("35d90a", No),
        ("35e08080", No),
        ("35e09fbf", No),
        ("35e081a0", No),
        ("35e0", No),
        ("35e0a5", No),
        ("35e041a6", Dies),
        ("35e0a541", No),
        ("35f0808080", No),
        ("35f08fbfbf", No),
        ("35f09d9f", No),
        ("35f880808080", No),
        ("35f888", No),
        ("35fe", No),
        ("35ff", No),
        ("8035", No),
        ("c2b235", No),
        ("d9", No),
        ("c0b5", Dies),
        ("e0a5", No),
        ("35d9a180", No),
        ("35d9a1d9", No),
        ("350a80", Dies),
        ("350ad9", No),
        ("354180", Dies),
        ("35800a", Dies),
        ("80", Dies),
        ("bf", Dies),
        ("800a", Dies),
        ("800a0a", No),
        ("c0", Dies),
        ("c1", Dies),
        ("c2", No),
        ("ff", No),
        ("fe", No),
        ("f04180", Dies),
        ("e080", No),
        ("e041", No),
        ("e08080", No),
        ("f080", No),
        ("f880", No),
        ("e00a80", Dies),
        ("e00a", No),
        ("ff80", No),
        ("ff4180", Dies),
        ("d9a180", No),
        ("d98080", No),
        ("4180", Dies),
        ("3080", Dies),
        ("358080", Dies),
        ("35c28080", No),
        ("35f48fbfbf80", No),
        ("35d9a1800a", No),
        ("35c1", Dies),
        ("35c10a", Dies),
        ("35c20a", No),
        ("353580", Dies),
        ("350a0a", No),
        ("35800a0a", No),
        ("80350a", No),
        ("35eda08080", No),
        ("35fdbfbfbfbfbf", No),
    ] {
        assert_eq!(answer(hex, true), expected, "{hex}");
    }
}

// The same match under the flag where the first byte announces a character of the
// string's whole length, which hides the newline that ends the string, and on strings a
// byte more or less than that.
#[test]
fn a_first_byte_that_announces_the_whole_string_hides_the_final_newline() {
    for (hex, expected) in [
        ("c00a", No),
        ("c10a", No),
        ("c1800a", Dies),
        ("c2c00a", Dies),
        ("dfc00a", Dies),
        ("e0c0", Dies),
        ("e0c00a", No),
        ("e0c00a0a", No),
        ("edc00a", No),
        ("efc10a", No),
        ("e041c00a", Dies),
        ("f041c00a", No),
        ("f4c1c00a", No),
        ("f7c0c00a", No),
        ("f04141c00a", Dies),
        ("f84141c00a", No),
        ("fb4141c10a", No),
        ("f8414141c00a", Dies),
        ("fc414141c00a", No),
        ("fd414141c00a", No),
        ("fc0a9fbcaa0a", No),
        ("fc41414141c00a", Dies),
        ("fe41414141c00a", No),
        ("fe4141414141c00a", Dies),
    ] {
        assert_eq!(answer(hex, true), expected, "{hex}");
        assert_eq!(answer(hex, false), No, "{hex}");
    }
    // The byte 0xff announces 13 bytes: itself, the filler, the two bytes at the end.
    for (filler, expected) in [(9, Dies), (10, No), (11, Dies)] {
        let hex = format!("ff{}c00a", "41".repeat(filler));
        assert_eq!(answer(&hex, true), expected, "{hex}");
    }
}

// The characters chr($_) =~ /^\d$/ holds for, 0 to 0x10FFFF: 650 in 61 ranges, the ASCII
// digits among them.
#[test]
fn the_digit_table_holds_the_decimal_digits_past_ascii_in_ascending_ranges() {
    let mut digits = 0;
    let mut after = 0x39;
    for (first, last) in DECIMAL_DIGITS {
        assert!(after < first && first <= last, "{first:x}-{last:x}");
        digits += last - first + 1;
        after = last;
    }
    assert_eq!(DECIMAL_DIGITS.len() + 1, 61);
    assert_eq!(digits + 10, 650);
}

fn row(name: &[u8], length: &[u8], utf8: bool, circular: Option<bool>) -> (Vec<u8>, bool) {
    let mut out = Vec::new();
    let wide = SequenceRegion {
        length: length.to_vec(),
        utf8,
        circular,
    }
    .write_row(name, &mut out);
    (out, wide)
}

// print $ofh join("\t", $chr, $length, defined $circular ? ($circular ? 1 : 0) : '-'), "\n";
#[test]
fn a_row_is_the_name_the_length_and_the_flag_between_tabs() {
    assert_eq!(
        row(b"1", b"5000", false, Some(false)),
        (b"1\t5000\t0\n".to_vec(), false)
    );
    assert_eq!(
        row(b"MT", b"16000", false, Some(true)),
        (b"MT\t16000\t1\n".to_vec(), false)
    );
    assert_eq!(
        row(b"scaffold_a", b"700", false, None),
        (b"scaffold_a\t700\t-\n".to_vec(), false)
    );
    // A length that ends in a newline is printed with it.
    assert_eq!(
        row(b"1", b"5000\n", false, Some(true)),
        (b"1\t5000\n\t1\n".to_vec(), false)
    );
}

// The same print for the chromosome directory "caf\xe9": beside the length 77, beside
// "7000" under the UTF-8 flag, and beside "5\x{661}", for which Perl warns
// `Wide character in print`.
#[test]
fn a_digit_past_ascii_writes_the_whole_row_as_utf8_and_warns() {
    assert_eq!(
        row(b"caf\xe9", b"77", false, Some(true)),
        (b"caf\xe9\t77\t1\n".to_vec(), false)
    );
    assert_eq!(
        row(b"caf\xe9", b"7000", true, Some(true)),
        (b"caf\xe9\t7000\t1\n".to_vec(), false)
    );
    assert_eq!(
        row(b"caf\xe9", b"5\xd9\xa1", true, Some(true)),
        (b"caf\xc3\xa9\t5\xd9\xa1\t1\n".to_vec(), true)
    );
}
