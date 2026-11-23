// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The field paths of a transcript and the kind of value under each, on both sides of
//! the conversion.
//!
//! A path joins the keys from the transcript down with `.`; the items of an array share
//! the path `<array>[]`, so a field any item holds is one path. A kind is `scalar`,
//! `undef`, `HASH` or `ARRAY`, with the class before it for a blessed value
//! (`Bio::EnsEMBL::Exon(HASH)`), and `(HASH,cycle)` for a hash that is being walked
//! further up the path, which is not walked again.

use std::collections::{BTreeMap, HashSet};

use crate::jsonpp::{string_text, Value as Json};
use crate::perl::Sv;
use crate::shards::Shard;
use crate::storable::{NodeId, Value};

/// The kinds of a transcript's paths. Paths are the UTF-8 of their characters.
pub type Fields = BTreeMap<Vec<u8>, Vec<u8>>;

/// An array met again on its own path. Where arrays alone stand between, the same walk
/// in Perl has no end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArrayCycle;

/// `$prefix ? "$prefix.$k" : $k`: a prefix that is empty or `0` is a false string.
fn child_path(prefix: &[u8], key: &[u8]) -> Vec<u8> {
    if prefix.is_empty() || prefix == b"0" {
        return key.to_vec();
    }
    let mut path = Vec::with_capacity(prefix.len() + 1 + key.len());
    path.extend_from_slice(prefix);
    path.push(b'.');
    path.extend_from_slice(key);
    path
}

fn item_path(prefix: &[u8]) -> Vec<u8> {
    [prefix, b"[]"].concat()
}

/// What a scalar of a Storable image is to the walk.
enum Target {
    Hash(NodeId),
    Array(NodeId),
    /// A reference to anything else: `SCALAR`, `REF`, `REGEXP` or `VSTRING`.
    Other(NodeId, &'static str),
    Defined,
    Undef,
}

struct StorableWalk<'a> {
    shard: &'a Shard,
    out: Fields,
    /// The hashes on the path from the transcript down.
    on_path: HashSet<NodeId>,
    /// The arrays on the path: the walk's own rule stops at a hash on its path alone, and
    /// an array met again is an `ArrayCycle`.
    arrays_on_path: HashSet<NodeId>,
}

impl StorableWalk<'_> {
    fn target(&self, sv: Sv) -> Target {
        let image = &self.shard.image;
        let Sv::Node(id) = self.shard.read(sv) else {
            return Target::Undef;
        };
        match image.value(id) {
            Value::Ref { target, .. } => match image.value(target) {
                Value::Hash { .. } => Target::Hash(target),
                Value::Array(_) => Target::Array(target),
                Value::Ref { .. } => Target::Other(target, "REF"),
                Value::Regexp(_) => Target::Other(target, "REGEXP"),
                _ if image.vstring(target).is_some() => Target::Other(target, "VSTRING"),
                _ => Target::Other(target, "SCALAR"),
            },
            Value::Undef | Value::ImmortalUndef | Value::ElementUndef => Target::Undef,
            _ => Target::Defined,
        }
    }

    /// `_kind`: the class and the kind of what a reference points at.
    fn reference_kind(&self, target: NodeId, kind: &str) -> Vec<u8> {
        // A class whose name is a false string counts as none.
        match self.shard.image.class(target) {
            Some(b"0") | None => kind.as_bytes().to_vec(),
            Some(b"") => format!("main({kind})").into_bytes(),
            Some(class) => [class, b"(", kind.as_bytes(), b")"].concat(),
        }
    }

    fn field(&mut self, child: Sv, path: Vec<u8>) -> Result<(), ArrayCycle> {
        match self.target(child) {
            Target::Hash(hash) => {
                let mut kind = self.reference_kind(hash, "HASH");
                if self.on_path.contains(&hash) {
                    // A blessed kind ends in `)` and takes the marker; `HASH` alone stays.
                    if kind.ends_with(b")") {
                        kind.pop();
                        kind.extend_from_slice(b",cycle)");
                    }
                    self.out.insert(path, kind);
                    return Ok(());
                }
                self.out.insert(path.clone(), kind);
                self.hash(hash, &path)
            }
            Target::Array(array) => {
                let kind = self.reference_kind(array, "ARRAY");
                self.out.insert(path.clone(), kind);
                self.array(array, &path)
            }
            Target::Other(target, kind) => {
                let kind = self.reference_kind(target, kind);
                self.out.insert(path, kind);
                Ok(())
            }
            Target::Defined => {
                self.out.insert(path, b"scalar".to_vec());
                Ok(())
            }
            Target::Undef => {
                self.out.insert(path, b"undef".to_vec());
                Ok(())
            }
        }
    }

    fn hash(&mut self, hash: NodeId, prefix: &[u8]) -> Result<(), ArrayCycle> {
        if !self.on_path.insert(hash) {
            return Ok(());
        }
        for (key, value) in crate::dumper::hash_entries(&self.shard.image, hash) {
            let path = child_path(prefix, &key.characters());
            self.field(Sv::Node(value), path)?;
        }
        self.on_path.remove(&hash);
        Ok(())
    }

    fn array(&mut self, array: NodeId, prefix: &[u8]) -> Result<(), ArrayCycle> {
        if !self.arrays_on_path.insert(array) {
            return Err(ArrayCycle);
        }
        let items = self.shard.image.items(array).unwrap_or(&[]).to_vec();
        for item in items {
            let element = match self.shard.image.value(item) {
                Value::ImmortalUndef => Sv::Missing,
                _ => Sv::Node(item),
            };
            self.field(element, item_path(prefix))?;
        }
        self.arrays_on_path.remove(&array);
        Ok(())
    }
}

/// The fields of a transcript of a Storable image: `transcript` is the hash.
pub fn storable_fields(shard: &Shard, transcript: NodeId) -> Result<Fields, ArrayCycle> {
    let mut walk = StorableWalk {
        shard,
        out: Fields::new(),
        on_path: HashSet::new(),
        arrays_on_path: HashSet::new(),
    };
    walk.hash(transcript, b"")?;
    Ok(walk.out)
}

/// The value an object holds under a key: its last entry of that key.
pub fn json_get<'v>(
    text: &[u8],
    entries: &'v [((usize, usize), Json)],
    key: &[u8],
) -> Option<&'v Json> {
    entries
        .iter()
        .rev()
        .find(|((from, to), _)| {
            let raw = &text[*from..*to];
            if raw.contains(&b'\\') {
                string_text(raw) == key
            } else {
                raw == key
            }
        })
        .map(|(_, value)| value)
}

fn json_field(text: &[u8], value: &Json, path: Vec<u8>, out: &mut Fields) {
    match value {
        Json::Dict(entries) => {
            out.insert(path.clone(), b"HASH".to_vec());
            json_object(text, entries, &path, out);
        }
        Json::List(items) => {
            out.insert(path.clone(), b"ARRAY".to_vec());
            for item in items {
                json_field(text, item, item_path(&path), out);
            }
        }
        Json::Null => {
            out.insert(path, b"undef".to_vec());
        }
        Json::Bool(_) => {
            out.insert(path, b"JSON::PP::Boolean(SCALAR)".to_vec());
        }
        Json::Number(..) | Json::Str(..) => {
            out.insert(path, b"scalar".to_vec());
        }
    }
}

fn json_object(text: &[u8], entries: &[((usize, usize), Json)], prefix: &[u8], out: &mut Fields) {
    // A key written twice holds its last value: the keys in the order of their names,
    // each once.
    let mut keys: Vec<(Vec<u8>, &Json)> = entries
        .iter()
        .map(|((from, to), value)| (string_text(&text[*from..*to]), value))
        .collect();
    keys.reverse();
    keys.sort_by(|left, right| left.0.cmp(&right.0));
    keys.dedup_by(|later, earlier| later.0 == earlier.0);
    for (key, value) in keys {
        json_field(text, value, child_path(prefix, &key), out);
    }
}

/// The fields of a transcript of a JSON shard: `entries` are the object's.
pub fn json_fields(text: &[u8], entries: &[((usize, usize), Json)]) -> Fields {
    let mut out = Fields::new();
    json_object(text, entries, b"", &mut out);
    out
}

#[cfg(test)]
#[path = "fields_tests.rs"]
mod tests;
