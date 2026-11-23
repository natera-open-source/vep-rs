// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Counts by name, in the order the names were first counted.

use indexmap::IndexMap;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Counter(pub IndexMap<String, i64>);

impl Counter {
    pub fn new() -> Counter {
        Counter(IndexMap::new())
    }

    /// The count of `name`, 0 for a name never counted.
    pub fn get(&self, name: &str) -> i64 {
        self.0.get(name).copied().unwrap_or(0)
    }

    pub fn add(&mut self, name: &str, count: i64) {
        match self.0.get_mut(name) {
            Some(held) => *held += count,
            None => {
                self.0.insert(name.to_string(), count);
            }
        }
    }

    pub fn contains(&self, name: &str) -> bool {
        self.0.contains_key(name)
    }

    /// Adds every count of `other`, a name new to `self` placed after the names it holds.
    pub fn update(&mut self, other: &Counter) {
        for (name, count) in &other.0 {
            self.add(name, *count);
        }
    }

    /// Removes `name`, the names after it keeping their order.
    pub fn remove(&mut self, name: &str) {
        self.0.shift_remove(name);
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    pub fn iter(&self) -> impl Iterator<Item = (&str, i64)> {
        self.0.iter().map(|(name, count)| (name.as_str(), *count))
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_name_never_counted_reads_zero_and_adding_inserts_it() {
        let mut counts = Counter::new();
        assert_eq!(counts.get("a"), 0);
        assert!(!counts.contains("a"));
        counts.add("b", 2);
        counts.add("a", 1);
        counts.add("b", 3);
        assert_eq!(counts.get("b"), 5);
        assert_eq!(counts.iter().collect::<Vec<_>>(), [("b", 5), ("a", 1)]);
        counts.remove("b");
        assert_eq!(counts.iter().collect::<Vec<_>>(), [("a", 1)]);
    }
}
