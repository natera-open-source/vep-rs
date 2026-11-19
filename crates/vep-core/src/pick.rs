// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The category sort Ensembl VEP picks one transcript row with: the
//! `pick_worst_VariationFeatureOverlapAllele` loop of `OutputFactory.pm:773-808`
//! (`--pick` and its family) and the `_pick_likely_transcript` loop of
//! `VariationFeatureAdaptor.pm:2200-2225` (the transcript a gene name resolves
//! to), which rank the candidates by one category at a time.
//!
//! Each candidate holds one value per category in an info hash. Every value is
//! a number the code assigns (`0` or `1`, a negated length, a rank) except the
//! TSL, which is the digit run the attribute's value matched (`$1`) kept as a
//! string, or the number `100` without one. The loop sorts by the value as a
//! number (`<=>`) and tests a tie as text (`eq`), and the two differ only on a
//! TSL string that is not the canonical spelling of its number: `007` sorts
//! with `7` and ties with neither it nor the default.
//!
//! Perl citations name modules of ensembl-vep release/116.2 (`Bio/EnsEMBL/VEP/...`)
//! and ensembl-variation release/116 (`Bio/EnsEMBL/Variation/...`).

use std::cmp::Ordering;

/// One category's value for one candidate: a number the code assigned, or the
/// TSL attribute's digit run as text, the one value the loops keep as a string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rank<'a> {
    Number(i64),
    Text(&'a str),
}

/// A value as Perl's `<=>` reads it: an integer within 64 bits exactly, a
/// longer digit run as a double.
#[derive(Clone, Copy)]
enum Numeric {
    Int(i128),
    Float(f64),
}

impl Rank<'_> {
    /// The TSL of a transcript from its level ([`crate::transcript::tsl_level`]),
    /// the number 100 without one (`OutputFactory.pm:726`,
    /// `VariationFeatureAdaptor.pm:2163`).
    pub fn tsl(level: Option<&str>) -> Rank<'_> {
        match level {
            Some(level) => Rank::Text(level),
            None => Rank::Number(100),
        }
    }

    fn numeric(self) -> Numeric {
        match self {
            Rank::Number(n) => Numeric::Int(i128::from(n)),
            Rank::Text(digits) => match digits.parse::<u64>() {
                Ok(n) => Numeric::Int(i128::from(n)),
                // Past 64 bits Perl converts the string to a double, infinite
                // past the double's range.
                Err(_) => Numeric::Float(digits.parse::<f64>().unwrap_or(f64::INFINITY)),
            },
        }
    }

    /// Perl's `<=>` between two values of one category
    /// (`OutputFactory.pm:789`, `VariationFeatureAdaptor.pm:2207`), as a total
    /// order: two integers compare exactly, two doubles as doubles, and a
    /// double, which only a digit run past 64 bits yields, sorts after every
    /// integer. Perl rounds the integer to a double for that last pair, so an
    /// integer within one rounding step of the double ties with it there; the
    /// sort needs one consistent answer, and the exact one keeps the order
    /// transitive.
    pub fn numeric_cmp(self, other: Self) -> Ordering {
        match (self.numeric(), other.numeric()) {
            (Numeric::Int(a), Numeric::Int(b)) => a.cmp(&b),
            (Numeric::Int(_), Numeric::Float(_)) => Ordering::Less,
            (Numeric::Float(_), Numeric::Int(_)) => Ordering::Greater,
            (Numeric::Float(a), Numeric::Float(b)) => a.total_cmp(&b),
        }
    }

    /// Perl's `eq` between two values of one category
    /// (`OutputFactory.pm:796`, `VariationFeatureAdaptor.pm:2214`): a number's
    /// text is its decimal, so a digit run equals a number only when it is that
    /// decimal, with no leading zero.
    pub fn text_eq(self, other: Self) -> bool {
        match (self, other) {
            (Rank::Number(a), Rank::Number(b)) => a == b,
            (Rank::Text(a), Rank::Text(b)) => a == b,
            (Rank::Number(n), Rank::Text(digits)) | (Rank::Text(digits), Rank::Number(n)) => {
                let canonical = digits == "0" || !digits.starts_with('0');
                canonical && digits.parse::<i64>() == Ok(n)
            }
        }
    }
}

/// The candidate the loop returns from `candidates`, `None` for none
/// (`OutputFactory.pm:773-811`, `VariationFeatureAdaptor.pm:2204-2225`): for
/// each category in order, the candidates left are sorted by the category's
/// value as a number (a stable sort, so equal values keep their order), the
/// first is taken, and those following it whose value is the same as text stay
/// for the next category; one alone is the pick, as is the first left after the
/// last category. `value` is one candidate's value in one category.
pub fn pick_by_categories<'a, C: Copy>(
    candidates: &[usize],
    categories: &[C],
    value: impl Fn(usize, C) -> Rank<'a>,
) -> Option<usize> {
    let mut left: Vec<usize> = candidates.to_vec();
    if left.is_empty() {
        return None;
    }
    for &category in categories {
        left.sort_by(|&a, &b| value(a, category).numeric_cmp(value(b, category)));
        let picked = left[0];
        let picked_value = value(picked, category);
        let tied = left
            .iter()
            .take_while(|&&candidate| value(candidate, category).text_eq(picked_value))
            .count();
        if tied == 1 {
            return Some(picked);
        }
        left.truncate(tied);
    }
    Some(left[0])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A TSL digit run compares as its number, exactly to 64 bits and as a
    /// double past them (`"300" <=> "2"` is 1), and the default is the number
    /// 100, so `tsl300` ranks after a transcript without a level and `tsl00`
    /// before every level.
    #[test]
    fn a_digit_run_compares_as_its_number() {
        let tsl = |level: &'static str| Rank::tsl(Some(level));
        assert_eq!(tsl("300").numeric_cmp(tsl("2")), Ordering::Greater);
        assert_eq!(tsl("300").numeric_cmp(Rank::tsl(None)), Ordering::Greater);
        assert_eq!(tsl("00").numeric_cmp(tsl("1")), Ordering::Less);
        assert_eq!(tsl("007").numeric_cmp(tsl("7")), Ordering::Equal);
        assert_eq!(tsl("0100").numeric_cmp(Rank::tsl(None)), Ordering::Equal);
        assert_eq!(
            tsl("18446744073709551615").numeric_cmp(tsl("18446744073709551614")),
            Ordering::Greater
        );
        assert_eq!(
            tsl("184467440737095516150").numeric_cmp(tsl("18446744073709551615")),
            Ordering::Greater
        );
        assert_eq!(
            Rank::Number(-5).numeric_cmp(Rank::Number(-4)),
            Ordering::Less
        );
    }

    /// The order is total: a run past 64 bits sorts after every integer, so
    /// two 64-bit integers one apart, which both round to the same double,
    /// still compare exactly against each other and consistently against it,
    /// and a sort over many such values is well defined.
    #[test]
    fn a_run_past_sixty_four_bits_sorts_after_every_integer() {
        let past = Rank::Text("18446744073709551616");
        let high = Rank::Text("18446744073709551615");
        let low = Rank::Text("18446744073709551614");
        assert_eq!(high.numeric_cmp(past), Ordering::Less);
        assert_eq!(low.numeric_cmp(past), Ordering::Less);
        assert_eq!(past.numeric_cmp(low), Ordering::Greater);
        assert_eq!(high.numeric_cmp(low), Ordering::Greater);
        let further = Rank::Text("200000000000000000000");
        assert_eq!(past.numeric_cmp(further), Ordering::Less);
        assert_eq!(further.numeric_cmp(further), Ordering::Equal);
        let mut many: Vec<Rank> = [past, high, low, Rank::Number(100)]
            .into_iter()
            .cycle()
            .take(30)
            .collect();
        many.sort_by(|a, b| a.numeric_cmp(*b));
        assert_eq!(many[0], Rank::Number(100));
        assert_eq!(many[29], past);
    }

    /// The tie test reads the text: `007` is not `7`, `0100` is not the default
    /// `100`, and `100` is.
    #[test]
    fn a_tie_is_the_same_text() {
        let tsl = |level: &'static str| Rank::tsl(Some(level));
        assert!(tsl("7").text_eq(tsl("7")));
        assert!(!tsl("007").text_eq(tsl("7")));
        assert!(tsl("100").text_eq(Rank::tsl(None)));
        assert!(!tsl("0100").text_eq(Rank::tsl(None)));
        assert!(Rank::Number(0).text_eq(tsl("0")));
        assert!(!Rank::Number(0).text_eq(tsl("00")));
        assert!(Rank::Number(3).text_eq(Rank::Number(3)));
    }

    /// The loop over two categories: the numeric minimum wins a category; a
    /// candidate tied as a number but not as text is dropped, and the tied
    /// prefix goes on to the next category.
    #[test]
    fn the_loop_shrinks_to_the_text_ties_of_the_numeric_minimum() {
        // Candidates: (tsl, length).
        let infos: [(Rank, i64); 4] = [
            (Rank::Text("7"), -10),
            (Rank::Text("007"), -50),
            (Rank::Text("7"), -40),
            (Rank::Number(100), -90),
        ];
        let value = |c: usize, category: u8| match category {
            0 => infos[c].0,
            _ => Rank::Number(infos[c].1),
        };
        // `007` sorts with the two `7`s but is not their tie, so the two `7`s
        // go on alone and the longer wins the length category.
        assert_eq!(pick_by_categories(&[0, 2, 1, 3], &[0u8, 1], value), Some(2));
        // Between the two `7`s, `007` ends the tied prefix at the first `7`,
        // which is the pick; listed first, `007` is the pick itself.
        assert_eq!(pick_by_categories(&[0, 1, 2, 3], &[0u8, 1], value), Some(0));
        assert_eq!(pick_by_categories(&[1, 0, 2, 3], &[0u8, 1], value), Some(1));
        // A single category with a full tie keeps the first candidate.
        assert_eq!(pick_by_categories(&[2, 0], &[0u8], value), Some(2));
        assert_eq!(pick_by_categories(&[], &[0u8], value), None);
    }
}
