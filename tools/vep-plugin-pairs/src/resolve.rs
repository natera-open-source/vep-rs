// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The fields a run has met and what settles a new one: the source of a field name and the
//! canonical form of a Location the passes leave unsettled, each field asked once.

use crate::order::cmp_column;
use crate::rules::Res;
use std::collections::HashMap;
use std::sync::{Mutex, RwLock};

/// What a field name is to the scored run.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Source {
    /// A column of Ensembl VEP's own.
    Core,
    /// A key no source claims.
    Unclaimed,
    /// A field of a source, by its number in the order the fields are first met.
    Field(u32),
}

/// What a resolver says of a field name.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Answer {
    /// A column of Ensembl VEP's own.
    Core,
    /// A key no source claims.
    Unclaimed,
    /// A field of the named source.
    Source(String),
}

/// What knows the sources of a run. The passes call it from any thread.
pub trait Resolver: Sync {
    /// The source of the field `name`; a refusal ends the run in its message.
    fn field(&self, name: &[u8]) -> Res<Answer>;

    /// The canonical form of a Location that `order::canonical_location` leaves unsettled.
    fn location(&self, written: &[u8]) -> Res<Vec<u8>>;
}

#[derive(Default)]
struct Claimed {
    names: Vec<Vec<u8>>,
    source: Vec<u32>,
    sources: Vec<String>,
}

/// The fields met so far and the resolver that settles a new one.
pub struct Registry<'a> {
    resolved: RwLock<HashMap<Vec<u8>, Source>>,
    claimed: Mutex<Claimed>,
    resolver: &'a (dyn Resolver + 'a),
    /// Held while a new field is asked, so that each is asked once.
    asking: Mutex<()>,
}

/// The claimed fields at one moment: per field its name, its source and its place among the
/// names in the order tuple lines sort them.
pub struct Table {
    pub names: Vec<Vec<u8>>,
    pub source: Vec<u32>,
    pub rank: Vec<u32>,
    pub sources: Vec<String>,
}

impl<'a> Registry<'a> {
    pub fn new(resolver: &'a (dyn Resolver + 'a)) -> Registry<'a> {
        Registry {
            resolved: RwLock::default(),
            claimed: Mutex::default(),
            resolver,
            asking: Mutex::default(),
        }
    }

    /// The source of `name`, asked of the resolver the first time the name is met.
    pub fn source_of(&self, name: &[u8]) -> Res<Source> {
        if let Some(&source) = self.resolved.read().unwrap().get(name) {
            return Ok(source);
        }
        let _asking = self.asking.lock().unwrap();
        if let Some(&source) = self.resolved.read().unwrap().get(name) {
            return Ok(source);
        }
        let source = match self.resolver.field(name)? {
            Answer::Core => Source::Core,
            Answer::Unclaimed => Source::Unclaimed,
            Answer::Source(source) => {
                let mut claimed = self.claimed.lock().unwrap();
                let known = claimed.sources.iter().position(|other| *other == source);
                let number = known.unwrap_or_else(|| {
                    claimed.sources.push(source);
                    claimed.sources.len() - 1
                });
                claimed.names.push(name.to_vec());
                claimed.source.push(number as u32);
                Source::Field(claimed.names.len() as u32 - 1)
            }
        };
        self.resolved.write().unwrap().insert(name.to_vec(), source);
        Ok(source)
    }

    /// The canonical form of a Location that `canonical_location` leaves to the resolver.
    pub fn location(&self, written: &[u8]) -> Res<Vec<u8>> {
        self.resolver.location(written)
    }

    pub fn table(&self) -> Table {
        let claimed = self.claimed.lock().unwrap();
        let mut order: Vec<u32> = (0..claimed.names.len() as u32).collect();
        order.sort_by(|&a, &b| cmp_column(&claimed.names[a as usize], &claimed.names[b as usize]));
        let mut rank = vec![0; order.len()];
        for (place, &field) in order.iter().enumerate() {
            rank[field as usize] = place as u32;
        }
        Table {
            names: claimed.names.clone(),
            source: claimed.source.clone(),
            rank,
            sources: claimed.sources.clone(),
        }
    }
}
