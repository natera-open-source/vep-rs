// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The bounded reservoir a category keeps its candidates in. Every candidate has a
//! priority, a hash of the seed and of what the candidate is, and the reservoir keeps the
//! lowest: what it holds at the end depends on the candidates alone, never on the order
//! they came in.

use std::collections::BTreeMap;

/// What a candidate is: the sequence by its number in reading order, a span on it and a
/// name. It breaks a tie between equal priorities and tells two candidates apart.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Key {
    pub sequence: usize,
    pub start: u64,
    pub end: u64,
    pub id: String,
}

/// The group a candidate is spread by: its sequence and its transcript's biotype.
pub type Stratum = (usize, String);

/// A reservoir of at most `want` candidates, kept twice over: the `want` strata whose
/// best candidate has the lowest priority, one candidate each, and the `want` candidates
/// of the lowest priority whatever their stratum.
#[derive(Debug, Clone)]
pub struct Reservoir<T> {
    want: usize,
    spread: BTreeMap<Stratum, (u64, Key, T)>,
    lowest: BTreeMap<(u64, Key), T>,
}

impl<T: Clone> Reservoir<T> {
    /// An empty reservoir that keeps `want` candidates.
    pub fn new(want: usize) -> Reservoir<T> {
        Reservoir {
            want,
            spread: BTreeMap::new(),
            lowest: BTreeMap::new(),
        }
    }

    /// Whether the reservoir holds no candidate.
    pub fn is_empty(&self) -> bool {
        self.lowest.is_empty()
    }

    /// The stratum whose kept candidate has the highest priority.
    fn last_stratum(&self) -> Option<(&Stratum, (u64, &Key))> {
        self.spread
            .iter()
            .map(|(stratum, (priority, key, _))| (stratum, (*priority, key)))
            .max_by(|a, b| a.1.cmp(&b.1))
    }

    fn spread_admits(&self, priority: u64, key: &Key, stratum: &Stratum) -> bool {
        if let Some((kept, kept_key, _)) = self.spread.get(stratum) {
            return (priority, key) < (*kept, kept_key);
        }
        self.spread.len() < self.want
            || self
                .last_stratum()
                .is_some_and(|(_, last)| (priority, key) < last)
    }

    fn lowest_admits(&self, priority: u64, key: &Key) -> bool {
        self.lowest.len() < self.want
            || self
                .lowest
                .last_key_value()
                .is_some_and(|((kept, kept_key), _)| (priority, key) < (*kept, kept_key))
    }

    /// Whether a candidate of this priority, key and stratum would be kept: asked before
    /// the candidate is made, so that only a candidate with a place is worked out.
    pub fn admits(&self, priority: u64, key: &Key, stratum: &Stratum) -> bool {
        self.want > 0
            && (self.spread_admits(priority, key, stratum) || self.lowest_admits(priority, key))
    }

    /// Keeps `item` where it has a place, and drops what it pushes out.
    pub fn insert(&mut self, priority: u64, key: Key, stratum: Stratum, item: T) {
        if self.want == 0 {
            return;
        }
        if self.spread_admits(priority, &key, &stratum) {
            self.spread
                .insert(stratum, (priority, key.clone(), item.clone()));
            if self.spread.len() > self.want {
                if let Some(last) = self.last_stratum().map(|(stratum, _)| stratum.clone()) {
                    self.spread.remove(&last);
                }
            }
        }
        if self.lowest_admits(priority, &key) {
            self.lowest.insert((priority, key), item);
            if self.lowest.len() > self.want {
                self.lowest.pop_last();
            }
        }
    }

    /// The kept candidates, `want` at most: one a stratum by priority, then the lowest
    /// priorities that are not among them.
    pub fn picks(&self) -> Vec<T> {
        let mut spread: Vec<&(u64, Key, T)> = self.spread.values().collect();
        spread.sort_by(|a, b| (a.0, &a.1).cmp(&(b.0, &b.1)));
        let mut picks: Vec<(&Key, &T)> = spread.iter().map(|(_, key, item)| (key, item)).collect();
        for ((_, key), item) in &self.lowest {
            if picks.len() >= self.want {
                break;
            }
            if !picks.iter().any(|(picked, _)| *picked == key) {
                picks.push((key, item));
            }
        }
        picks.truncate(self.want);
        picks.into_iter().map(|(_, item)| item.clone()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(id: &str) -> Key {
        Key {
            sequence: 0,
            start: 1,
            end: 2,
            id: id.to_owned(),
        }
    }

    /// Candidates of four strata; the priority is the number, the stratum its tens.
    fn candidates() -> Vec<(u64, Key, Stratum)> {
        [3, 11, 12, 14, 25, 31, 38, 40, 41, 47]
            .iter()
            .map(|n: &u64| (*n, key(&n.to_string()), (0, (n / 10).to_string())))
            .collect()
    }

    fn fill(want: usize, order: &[usize]) -> Vec<u64> {
        let all = candidates();
        let mut reservoir = Reservoir::new(want);
        for index in order {
            let (priority, key, stratum) = all[*index].clone();
            if reservoir.admits(priority, &key, &stratum) {
                reservoir.insert(priority, key, stratum, priority);
            }
        }
        reservoir.picks()
    }

    #[test]
    fn one_candidate_is_the_lowest_priority() {
        assert_eq!(fill(1, &[5, 0, 9, 3]), [3]);
        assert_eq!(fill(0, &[5, 0, 9, 3]), Vec::<u64>::new());
    }

    #[test]
    fn the_picks_spread_over_strata_then_fill_from_the_lowest() {
        let forward: Vec<usize> = (0..10).collect();
        assert_eq!(fill(3, &forward), [3, 11, 25]);
        assert_eq!(fill(4, &forward), [3, 11, 25, 31]);
        assert_eq!(fill(6, &forward), [3, 11, 25, 31, 40, 12]);
    }

    #[test]
    fn the_order_candidates_come_in_changes_nothing() {
        let forward: Vec<usize> = (0..10).collect();
        let backward: Vec<usize> = (0..10).rev().collect();
        let mixed = [4, 9, 0, 7, 2, 5, 1, 8, 3, 6];
        for want in 1..=10 {
            let expected = fill(want, &forward);
            assert_eq!(fill(want, &backward), expected, "want {want}");
            assert_eq!(fill(want, &mixed), expected, "want {want}");
        }
    }

    #[test]
    fn a_candidate_offered_twice_is_kept_once() {
        assert_eq!(fill(4, &[1, 1, 2, 1]), [11, 12]);
    }
}
