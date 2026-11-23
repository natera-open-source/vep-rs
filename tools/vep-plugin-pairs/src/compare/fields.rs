// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The plugins of one comparison and the registry entries each is compared on.

use crate::registry::{field_sets, FieldMatcher};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

/// The names met so far: each one's number, and the names by their numbers.
#[derive(Default)]
struct Numbered {
    numbers: HashMap<Arc<str>, u32>,
    names: Vec<Arc<str>>,
}

/// The field names a file pair's outputs carry for the run's plugins, each under one
/// number, so that a record in the scratch names its fields by number.
#[derive(Default)]
pub struct FieldNames {
    known: Mutex<Numbered>,
}

impl FieldNames {
    /// The number of a name, the next free one for a name not met before.
    pub fn number(&self, name: &str) -> u32 {
        let mut known = self.known.lock().expect("no reader panics");
        if let Some(number) = known.numbers.get(name) {
            return *number;
        }
        let (name, number) = (Arc::<str>::from(name), known.names.len() as u32);
        known.names.push(name.clone());
        known.numbers.insert(name, number);
        number
    }

    /// The names by their numbers.
    pub fn names(&self) -> Vec<Arc<str>> {
        self.known.lock().expect("no reader panics").names.clone()
    }
}

/// One plugin the run names.
pub struct Requested {
    pub name: String,
    /// The plugin's entries that no plugin named with it and standing before it in the
    /// registry's order has an entry for. Each compared record needs one owner, so a name
    /// two named plugins write counts for the first of them in that order.
    pub registered: Vec<String>,
    /// Matches `registered`: the fields read for the plugin.
    pub claims: FieldMatcher,
    /// Matches every entry of the plugin.
    pub writes: FieldMatcher,
}

/// The plugins of a run.
pub struct Run {
    /// The names the command line gives, in its order, a repeated name as often as given.
    pub plugins: Vec<String>,
    /// Each named plugin once, in the order it is first named.
    pub distinct: Vec<Requested>,
    /// Per name of `plugins`, its place in `distinct`.
    pub slots: Vec<usize>,
}

impl Run {
    /// The run over `plugins`, every one a field set of the registry. `fields`, the
    /// command line's override, stands in for the entries of the run's only plugin.
    pub fn new(plugins: Vec<String>, fields: Option<Vec<String>>) -> Run {
        let mut in_registry_order: Vec<Requested> = Vec::new();
        for set in field_sets() {
            if !plugins.contains(&set.name) {
                continue;
            }
            let entries = match &fields {
                Some(given) if plugins.len() == 1 => given,
                _ => &set.entries,
            };
            let registered: Vec<String> = entries
                .iter()
                .filter(|entry| {
                    !in_registry_order
                        .iter()
                        .any(|earlier| earlier.writes.claims(entry))
                })
                .cloned()
                .collect();
            in_registry_order.push(Requested {
                name: set.name.clone(),
                claims: FieldMatcher::new(&registered),
                registered,
                writes: FieldMatcher::new(entries),
            });
        }
        let mut distinct: Vec<Requested> = Vec::new();
        let mut slots = Vec::with_capacity(plugins.len());
        for name in &plugins {
            let slot = match distinct.iter().position(|known| &known.name == name) {
                Some(slot) => slot,
                None => {
                    let at = in_registry_order
                        .iter()
                        .position(|known| &known.name == name)
                        .expect("a named plugin is a field set of the registry");
                    distinct.push(in_registry_order.remove(at));
                    distinct.len() - 1
                }
            };
            slots.push(slot);
        }
        Run {
            plugins,
            distinct,
            slots,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(names: &[&str]) -> Vec<String> {
        names.iter().map(|name| name.to_string()).collect()
    }

    #[test]
    fn a_name_two_named_plugins_write_is_registered_for_the_first_in_registry_order() {
        let run = Run::new(names(&["TranscriptAnnotator", "dbNSFP"]), None);
        assert_eq!(run.distinct[0].name, "TranscriptAnnotator");
        assert_eq!(
            run.distinct[0].registered,
            ["transcript_annotator_<column>", "Comment"]
        );
        assert!(run.distinct[1]
            .registered
            .iter()
            .any(|entry| entry == "SIFT_pred"));
        assert!(run.distinct[0].writes.claims("SIFT_pred"));
        assert!(!run.distinct[0].claims.claims("SIFT_pred"));
        let alone = Run::new(names(&["TranscriptAnnotator"]), None);
        assert_eq!(alone.distinct[0].registered.len(), 4);
    }

    #[test]
    fn a_repeated_plugin_is_one_plugin_of_the_run() {
        let run = Run::new(names(&["REVEL", "CADD", "REVEL"]), None);
        assert_eq!(run.plugins, ["REVEL", "CADD", "REVEL"]);
        assert_eq!(run.slots, [0, 1, 0]);
        assert_eq!(run.distinct.len(), 2);
    }

    #[test]
    fn the_override_stands_in_for_the_entries_of_an_only_plugin() {
        let fields = Some(names(&["LoF", "x_<N>"]));
        let run = Run::new(names(&["LoFTEE"]), fields.clone());
        assert_eq!(run.distinct[0].registered, ["LoF", "x_<N>"]);
        assert!(
            run.distinct[0].claims.claims("x_7") && !run.distinct[0].claims.claims("LoF_filter")
        );
        let two = Run::new(names(&["LoFTEE", "REVEL"]), fields);
        assert!(two.distinct[0].claims.claims("LoF_filter"));
    }
}
