// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The counts of one projection of one engine.

use super::api::CountedProjection;
use super::error::{Scored, ScorerError};
use super::sources::text;
use std::collections::{BTreeMap, BTreeSet};

/// Per field name: the tuples paired, those paired as numbers alone and those each side has
/// left. Per source: the tuples of each side set aside, per class the keys under which it
/// set a tuple aside (`masked`) or would if it were applied (`not_applied`), and the keys
/// under which a side writes one field more than once (`shared`).
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Tally {
    pub paired: BTreeMap<Vec<u8>, u64>,
    pub text_only: BTreeMap<Vec<u8>, u64>,
    pub ref_left: BTreeMap<Vec<u8>, u64>,
    pub eng_left: BTreeMap<Vec<u8>, u64>,
    pub ref_aside: BTreeMap<String, u64>,
    pub eng_aside: BTreeMap<String, u64>,
    pub masked: BTreeMap<String, BTreeMap<String, u64>>,
    pub not_applied: BTreeMap<String, BTreeMap<String, u64>>,
    pub shared: BTreeMap<String, u64>,
}

fn held(counts: &BTreeMap<Vec<u8>, u64>, field: &[u8]) -> u64 {
    counts.get(field).copied().unwrap_or(0)
}

impl Tally {
    /// Every field name the tally counts a tuple of, paired or left, in byte order.
    pub fn fields(&self) -> BTreeSet<Vec<u8>> {
        self.paired
            .keys()
            .chain(self.ref_left.keys())
            .chain(self.eng_left.keys())
            .cloned()
            .collect()
    }

    /// Counts what the settling step found of one field under one key: the pairs, those
    /// equal as numbers alone among them, and the values each side has left.
    pub fn count(
        &mut self,
        name: &[u8],
        paired: u64,
        as_numbers: u64,
        ref_left: u64,
        eng_left: u64,
    ) {
        for (counts, more) in [
            (&mut self.paired, paired),
            (&mut self.text_only, as_numbers),
            (&mut self.ref_left, ref_left),
            (&mut self.eng_left, eng_left),
        ] {
            match counts.get_mut(name) {
                Some(count) => *count += more,
                None => {
                    counts.insert(name.to_vec(), more);
                }
            }
        }
    }

    /// Adds the counts of `other`, the tally of other keys of the same projection. A count
    /// of zero keeps its name.
    pub fn add(&mut self, other: Tally) {
        fn sum<K: Ord>(counts: &mut BTreeMap<K, u64>, more: BTreeMap<K, u64>) {
            for (name, count) in more {
                *counts.entry(name).or_default() += count;
            }
        }
        sum(&mut self.paired, other.paired);
        sum(&mut self.text_only, other.text_only);
        sum(&mut self.ref_left, other.ref_left);
        sum(&mut self.eng_left, other.eng_left);
        sum(&mut self.ref_aside, other.ref_aside);
        sum(&mut self.eng_aside, other.eng_aside);
        sum(&mut self.shared, other.shared);
        for (counts, more) in [
            (&mut self.masked, other.masked),
            (&mut self.not_applied, other.not_applied),
        ] {
            for (source, classes) in more {
                sum(counts.entry(source).or_default(), classes);
            }
        }
    }

    /// Adds what the passes counted of one projection: per field the pairs equal as text,
    /// and per source the keys under which a side writes one field more than once. The
    /// passes state each side's values of a field as well: their pairs and the values the
    /// settling step was handed of the field, or the counts and the remainder disagree,
    /// which is a refusal naming the first such field in byte order.
    pub fn take(&mut self, counted: &CountedProjection) -> Scored<()> {
        let stated: BTreeMap<&[u8], [u64; 3]> = counted
            .values
            .iter()
            .map(|(name, counts)| (name.as_slice(), *counts))
            .collect();
        let mut fields = self.fields();
        fields.extend(stated.keys().map(|name| name.to_vec()));
        for field in fields {
            let [reference, engine, paired] =
                stated.get(field.as_slice()).copied().unwrap_or([0, 0, 0]);
            let settled = i128::from(held(&self.paired, &field));
            let handed = (
                settled + i128::from(held(&self.ref_left, &field)),
                settled + i128::from(held(&self.eng_left, &field)),
            );
            let owed = (
                i128::from(reference) - i128::from(paired),
                i128::from(engine) - i128::from(paired),
            );
            if handed != owed {
                return Err(ScorerError(format!(
                    "the pairing tool's counts of {} and its remainder disagree",
                    text(&field)
                )));
            }
            *self.paired.entry(field).or_default() += paired;
        }
        for (source, keys) in &counted.shared_keys {
            *self.shared.entry(source.clone()).or_default() += keys;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counted(values: &[(&str, [u64; 3])], shared: &[(&str, u64)]) -> CountedProjection {
        CountedProjection {
            values: values
                .iter()
                .map(|(name, counts)| (name.as_bytes().to_vec(), *counts))
                .collect(),
            shared_keys: shared
                .iter()
                .map(|(source, keys)| (source.to_string(), *keys))
                .collect(),
        }
    }

    #[test]
    fn pairs_and_shared_keys_the_passes_counted_are_added() {
        let mut tally = Tally::default();
        tally.count(b"PLANT_score", 1, 1, 2, 0);
        tally
            .take(&counted(
                &[("PLANT_score", [10, 8, 7]), ("PLANT_label", [4, 4, 4])],
                &[("Planted", 3)],
            ))
            .unwrap();
        assert_eq!(tally.paired[b"PLANT_score".as_slice()], 8);
        assert_eq!(tally.paired[b"PLANT_label".as_slice()], 4);
        assert_eq!(tally.text_only[b"PLANT_score".as_slice()], 1);
        assert_eq!(tally.shared["Planted"], 3);
        let fields: Vec<Vec<u8>> = tally.fields().into_iter().collect();
        assert_eq!(fields, [b"PLANT_label".to_vec(), b"PLANT_score".to_vec()]);
    }

    #[test]
    fn counts_that_disagree_with_the_remainder_are_refused() {
        let mut short = Tally::default();
        short.count(b"PLANT_score", 0, 0, 1, 0);
        assert_eq!(
            short.take(&counted(&[("PLANT_score", [3, 1, 1])], &[])),
            Err(ScorerError(
                "the pairing tool's counts of PLANT_score and its remainder disagree".to_string()
            ))
        );
        let mut unstated = Tally::default();
        unstated.count(b"PLANT\xff", 0, 0, 0, 1);
        assert_eq!(
            unstated.take(&counted(&[], &[])),
            Err(ScorerError(
                "the pairing tool's counts of PLANT\u{fffd} and its remainder disagree".to_string()
            ))
        );
        let mut settled = Tally::default();
        settled.count(b"PLANT_score", 0, 0, 0, 0);
        assert_eq!(settled.take(&counted(&[], &[])), Ok(()));
        assert_eq!(settled.paired[b"PLANT_score".as_slice()], 0);
    }
}
