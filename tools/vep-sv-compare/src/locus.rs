// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What several masks read off a Location or a breakend allele: the span, the chromosome
//! and the mate, with the constants the masks share.
//!
//! A coordinate is read up to 128 bits; a longer one is no coordinate here.

use std::sync::LazyLock;

use regex::Regex;

/// The reference's breakend-to-feature distance gate, in bases.
pub const BREAKEND_FEATURE_DISTANCE: i128 = 5000;

/// The Allele the scored engine writes for a tandem-repeat record.
pub const CNV_TR_RUST_ALLELE: &str = "tandem_repeat";

/// The last three columns of the row a record with no overlapping feature gets.
pub const INTERGENIC_TUPLE_TAIL: [&str; 3] = ["-", "-", "intergenic_variant"];

static LOCATION_SPAN: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^([^:]+):(-?[0-9]+)(?:-(-?[0-9]+))?$").expect("regex"));
static BRACKET_MATE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]").expect("regex"));

/// `chrom:start-end` or `chrom:pos` as `(chrom, start, end)`, a coordinate below zero
/// included; `None` for any other text.
pub fn location_span(loc: &str) -> Option<(&str, i128, i128)> {
    let m = LOCATION_SPAN.captures(loc)?;
    let start: i128 = m[2].parse().ok()?;
    let end = match m.get(3) {
        Some(end) => end.as_str().parse().ok()?,
        None => start,
    };
    Some((m.get(1).expect("chrom group").as_str(), start, end))
}

/// A chromosome name without a leading `chr` in any case.
pub fn strip_chr(chrom: &str) -> &str {
    if chrom.len() >= 3 && chrom.is_char_boundary(3) && chrom[..3].eq_ignore_ascii_case("chr") {
        &chrom[3..]
    } else {
        chrom
    }
}

/// The chromosome a Location names, without a leading `chr`; `None` for a text with no
/// colon.
pub fn chromosome_of_location(loc: &str) -> Option<&str> {
    let (chrom, _) = loc.split_once(':')?;
    Some(strip_chr(chrom))
}

/// The value of a text of decimal digits of any script, as Python's `int` reads it. The
/// decimal digits of a script are ten consecutive characters from its zero, and the runs of
/// neighbouring scripts follow each other in whole tens, so a digit is its distance from the
/// start of its run, modulo ten.
fn decimal_value(digits: &str) -> Option<i128> {
    let mut value: i128 = 0;
    for c in digits.chars() {
        let digit = match c.to_digit(10) {
            Some(digit) => digit,
            None => {
                let mut zero = c as u32;
                while zero > 0 && char::from_u32(zero - 1).is_some_and(is_decimal_digit) {
                    zero -= 1;
                }
                (c as u32 - zero) % 10
            }
        };
        value = value.checked_mul(10)?.checked_add(i128::from(digit))?;
    }
    Some(value)
}

fn is_decimal_digit(c: char) -> bool {
    static DIGIT: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^\d$").expect("regex"));
    let mut buffer = [0u8; 4];
    DIGIT.is_match(c.encode_utf8(&mut buffer))
}

/// The mate `(chromosome, position)` a bracket-notation breakend allele names, the
/// chromosome without a leading `chr`; `None` for any other allele.
pub fn bracket_mate(allele: &str) -> Option<(&str, i128)> {
    let m = BRACKET_MATE.captures(allele)?;
    let chrom = m.get(1).expect("chrom group").as_str();
    Some((strip_chr(chrom), decimal_value(&m[2])?))
}

/// Whether the allele is bracket-notation breakend text, whose mate may be on another
/// chromosome.
pub fn allele_names_another_chromosome(allele: &str) -> bool {
    allele.contains('[') || allele.contains(']')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_location_is_a_span_or_a_point() {
        assert_eq!(location_span("21:5-9"), Some(("21", 5, 9)));
        assert_eq!(location_span("21:7"), Some(("21", 7, 7)));
        assert_eq!(location_span("21:-4-295"), Some(("21", -4, 295)));
        assert_eq!(location_span("21:5-9.5"), None);
        assert_eq!(location_span("21"), None);
        assert_eq!(location_span("21:5-9\n"), None);
    }

    #[test]
    fn a_chromosome_loses_its_prefix_in_any_case() {
        assert_eq!(chromosome_of_location("chr21:5"), Some("21"));
        assert_eq!(chromosome_of_location("CHR21:5"), Some("21"));
        assert_eq!(chromosome_of_location("21:5"), Some("21"));
        assert_eq!(chromosome_of_location("21"), None);
    }

    #[test]
    fn a_bracket_allele_names_its_mate() {
        assert_eq!(bracket_mate("N[chr22:44739371["), Some(("22", 44_739_371)));
        assert_eq!(bracket_mate("]Chr5:10]N"), Some(("5", 10)));
        assert_eq!(bracket_mate("N."), None);
        assert_eq!(bracket_mate("N[22:\u{ff11}\u{ff12}["), Some(("22", 12)));
        assert_eq!(bracket_mate("N[22:\u{0663}\u{0660}["), Some(("22", 30)));
        assert!(allele_names_another_chromosome("N[22:5["));
        assert!(!allele_names_another_chromosome("deletion"));
    }
}
