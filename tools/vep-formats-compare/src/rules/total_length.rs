// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The end key of a position range under `--total_length`, in the JSON format.

use crate::entry::Entry;
use crate::shared::blank;

/// The JSON end keys the class reads, each beside the start key of its coordinate.
pub const RANGE_END_FIELDS: [(&str, &str); 3] = [
    ("cdna_end", "cdna_start"),
    ("cds_end", "cds_start"),
    ("protein_end", "protein_start"),
];

/// A run of ASCII digits as a number; `None` for an empty text, another character or a
/// number `u64` does not hold.
fn position(text: &str) -> Option<u64> {
    if text.is_empty() || !text.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    text.parse().ok()
}

/// Whether the end key `name` of the reference's entry `e` and the candidate's entry `a` has
/// the shape of `total_length_json_range_end`, `ev` and `av` being its two values, `-` for
/// none. The candidate's value is a last position and a total (`487/2674`). For a range
/// whose start is known both entries carry the start key with one value, the reference's
/// end repeats it and the candidate's last position lies above it. For a range whose start
/// is unknown neither entry carries the start key and the reference carries no end.
pub fn total_length_range_end_pair(e: &Entry, a: &Entry, name: &str, ev: &str, av: &str) -> bool {
    let Some((_, start_key)) = RANGE_END_FIELDS.iter().find(|(end, _)| *end == name) else {
        return false;
    };
    let Some((last, total)) = av.split_once('/') else {
        return false;
    };
    let (Some(last), Some(_)) = (position(last), position(total)) else {
        return false;
    };
    let start_of = |entry: &'_ Entry| -> Option<String> {
        entry
            .fields
            .get(*start_key)
            .filter(|value| !blank(value))
            .cloned()
    };
    match (start_of(e), start_of(a)) {
        (Some(first), Some(candidate_first)) => {
            first == candidate_first
                && ev == first
                && position(&first).is_some_and(|first| first < last)
        }
        (None, None) => blank(ev),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(fields: &[(&str, &str)]) -> Entry {
        Entry {
            fields: fields
                .iter()
                .map(|(name, value)| (name.to_string(), value.to_string()))
                .collect(),
            ..Entry::default()
        }
    }

    #[test]
    fn a_range_end_the_reference_writes_as_its_start_is_taken() {
        let reference = entry(&[("cdna_start", "484"), ("cdna_end", "484")]);
        let candidate = entry(&[("cdna_start", "484"), ("cdna_end", "487/2674")]);
        assert!(total_length_range_end_pair(
            &reference, &candidate, "cdna_end", "484", "487/2674"
        ));
        for (end, start) in RANGE_END_FIELDS {
            let reference = entry(&[(start, "39"), (end, "39")]);
            let candidate = entry(&[(start, "39"), (end, "40/533")]);
            assert!(
                total_length_range_end_pair(&reference, &candidate, end, "39", "40/533"),
                "{end}"
            );
        }
    }

    #[test]
    fn a_range_of_unknown_start_is_taken_where_the_reference_writes_neither_key() {
        let reference = entry(&[("variant_allele", "-")]);
        let candidate = entry(&[("cdna_end", "487/2674")]);
        assert!(total_length_range_end_pair(
            &reference, &candidate, "cdna_end", "-", "487/2674"
        ));
        // An end the reference writes beside no start is another row.
        let reference = entry(&[("cdna_end", "487")]);
        assert!(!total_length_range_end_pair(
            &reference, &candidate, "cdna_end", "487", "487/2674"
        ));
    }

    #[test]
    fn every_other_pair_stays_a_difference() {
        let reference = entry(&[("cdna_start", "484"), ("cdna_end", "484")]);
        let pair = |candidate: &Entry, name: &str, ev: &str, av: &str| {
            total_length_range_end_pair(&reference, candidate, name, ev, av)
        };
        // An end without a total, an end at or below the start, and a total that is no number.
        let candidate = entry(&[("cdna_start", "484"), ("cdna_end", "487")]);
        assert!(!pair(&candidate, "cdna_end", "484", "487"));
        let candidate = entry(&[("cdna_start", "484"), ("cdna_end", "484/2674")]);
        assert!(!pair(&candidate, "cdna_end", "484", "484/2674"));
        let candidate = entry(&[("cdna_start", "484"), ("cdna_end", "480/2674")]);
        assert!(!pair(&candidate, "cdna_end", "484", "480/2674"));
        let candidate = entry(&[("cdna_start", "484"), ("cdna_end", "487/x")]);
        assert!(!pair(&candidate, "cdna_end", "484", "487/x"));
        // A start the two engines write differently, or one engine alone writes.
        let candidate = entry(&[("cdna_start", "485"), ("cdna_end", "487/2674")]);
        assert!(!pair(&candidate, "cdna_end", "484", "487/2674"));
        let candidate = entry(&[("cdna_end", "487/2674")]);
        assert!(!pair(&candidate, "cdna_end", "484", "487/2674"));
        // A reference end that is not its start, and a key the class does not read.
        let candidate = entry(&[("cdna_start", "484"), ("cdna_end", "487/2674")]);
        assert!(!pair(&candidate, "cdna_end", "486", "487/2674"));
        assert!(!pair(&candidate, "cdna_start", "484", "487/2674"));
        // A single position with its total is a string at both ends in both engines.
        let single = entry(&[("cdna_start", "484/2674"), ("cdna_end", "484/2674")]);
        assert!(!total_length_range_end_pair(
            &single, &candidate, "cdna_end", "484/2674", "487/2674"
        ));
    }
}
