// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rewriting of a Location takes the digits `0` to `9` alone.

use super::*;

fn canonical(location: &str) -> String {
    let mut written = Vec::new();
    push_canonical_location(location, &mut written);
    String::from_utf8(written).expect("UTF-8")
}

/// A coordinate in the decimal digits of another script, a `+` before it kept, and one that
/// ends in a superscript digit, stay as written.
#[test]
fn a_location_is_rewritten_where_its_coordinates_are_written_in_the_digits_0_to_9() {
    assert_eq!(canonical("21:005-+9"), "21:5-9");
    let other_digits = "21:\u{665}-\u{669}";
    assert_eq!(canonical(other_digits), other_digits);
    let signed = "21:+\u{665}";
    assert_eq!(canonical(signed), signed);
    assert_eq!(canonical("21:10\u{b2}"), "21:10\u{b2}");
}
