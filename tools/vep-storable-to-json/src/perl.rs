// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What Perl does with a retrieved value, for the steps of the conversion: `defined`, `ref`,
//! `reftype`, truth, a hash fetch, an array walk, a number, a string.
//!
//! A scalar the conversion holds is an [`Sv`]: a node of the image, or the undef Perl
//! yields for a key or element that is not there. A step that makes Perl die is a
//! [`Die`] whose text is Perl's message without the ` at <file> line <n>.` it appends.

use std::borrow::Cow;

use crate::json::PerlStr;
use crate::number::{double_to_perl_string, numify, write_int, Num};
use crate::storable::{Image, NodeId, Value};

/// The text of a Perl `die`, without a final newline.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Die(pub Vec<u8>);

impl Die {
    pub fn new(text: impl Into<Vec<u8>>) -> Self {
        Die(text.into())
    }
}

/// A scalar in the conversion's hands.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sv {
    Node(NodeId),
    /// Undef from a hash key that is absent or an array element that does not exist.
    Missing,
}

/// The kind of thing a reference points at, as `Scalar::Util`'s `reftype` names it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefType {
    Hash,
    Array,
    /// A scalar that is not itself a reference.
    Scalar,
    /// A scalar that is a reference.
    Ref,
    Regexp,
}

impl RefType {
    fn name(self) -> &'static [u8] {
        match self {
            RefType::Hash => b"HASH",
            RefType::Array => b"ARRAY",
            RefType::Scalar => b"SCALAR",
            RefType::Ref => b"REF",
            RefType::Regexp => b"Regexp",
        }
    }
}

/// The image with Perl's operations on it.
#[derive(Clone, Copy)]
pub struct Perl<'a> {
    pub image: &'a Image,
}

impl<'a> Perl<'a> {
    pub fn new(image: &'a Image) -> Self {
        Self { image }
    }

    pub fn value(&self, sv: Sv) -> Option<Value> {
        match sv {
            Sv::Node(id) => Some(self.image.value(id)),
            Sv::Missing => None,
        }
    }

    /// Perl's `defined`.
    pub fn defined(&self, sv: Sv) -> bool {
        match sv {
            Sv::Node(id) => self.image.is_defined(id),
            Sv::Missing => false,
        }
    }

    /// The node a reference points at; `None` for a scalar that is not a reference.
    pub fn target(&self, sv: Sv) -> Option<NodeId> {
        match sv {
            Sv::Node(id) => self.image.target(id),
            Sv::Missing => None,
        }
    }

    pub fn reftype_of_node(&self, target: NodeId) -> RefType {
        match self.image.value(target) {
            Value::Hash { .. } => RefType::Hash,
            Value::Array(_) => RefType::Array,
            Value::Ref { .. } => RefType::Ref,
            Value::Regexp(_) => RefType::Regexp,
            _ => RefType::Scalar,
        }
    }

    /// `Scalar::Util`'s `reftype`: `None` for a scalar that is not a reference.
    pub fn reftype(&self, sv: Sv) -> Option<RefType> {
        self.target(sv).map(|target| self.reftype_of_node(target))
    }

    /// Perl's `ref`: the class of a blessed target, else the name of its kind; `None`
    /// for a scalar that is not a reference.
    pub fn ref_name(&self, sv: Sv) -> Option<&'a [u8]> {
        let target = self.target(sv)?;
        Some(match self.image.class(target) {
            // Blessing into the empty name blesses into `main`.
            Some(b"") => b"main",
            Some(class) => class,
            None => self.reftype_of_node(target).name(),
        })
    }

    /// `ref $x` in a condition: false for a scalar that is not a reference, and for a
    /// reference blessed into the class `0`, whose name is a false string.
    pub fn ref_is_true(&self, sv: Sv) -> bool {
        self.ref_name(sv).is_some_and(|name| name != b"0")
    }

    /// `ref $x eq 'ARRAY'` and its like.
    pub fn ref_is(&self, sv: Sv, name: &[u8]) -> bool {
        self.ref_name(sv) == Some(name)
    }

    /// `(reftype($x) // '') eq 'HASH'`: the hash a reference points at, blessed or not.
    pub fn hash_target(&self, sv: Sv) -> Option<NodeId> {
        self.target(sv)
            .filter(|target| matches!(self.image.value(*target), Value::Hash { .. }))
    }

    /// `(reftype($x) // '') eq 'ARRAY'`: the array a reference points at.
    pub fn array_target(&self, sv: Sv) -> Option<NodeId> {
        self.target(sv)
            .filter(|target| matches!(self.image.value(*target), Value::Array(_)))
    }

    /// A fetch from a hash node.
    pub fn fetch(&self, hash: NodeId, key: &[u8]) -> Sv {
        match self.image.fetch(hash, key) {
            Some(id) => Sv::Node(id),
            None => Sv::Missing,
        }
    }

    /// The hash `$sv->{...}` reads from. An undef scalar comes to life as an empty
    /// hash, which reads as `None`; a reference to anything but a hash and a scalar
    /// that is no reference make Perl die.
    pub fn hash_of(&self, sv: Sv) -> Result<Option<NodeId>, Die> {
        if !self.defined(sv) {
            return Ok(None);
        }
        match self.target(sv) {
            Some(target) => match self.image.value(target) {
                Value::Hash { .. } => Ok(Some(target)),
                _ => Err(Die::new("Not a HASH reference")),
            },
            None => Err(self.strict_refs(sv, "a HASH")),
        }
    }

    /// `$sv->{key}`.
    pub fn get(&self, sv: Sv, key: &[u8]) -> Result<Sv, Die> {
        Ok(match self.hash_of(sv)? {
            Some(hash) => self.fetch(hash, key),
            None => Sv::Missing,
        })
    }

    /// The items `@{$sv}` walks: an element that does not exist reads as undef.
    pub fn items_of(&self, sv: Sv) -> Result<Vec<Sv>, Die> {
        let Some(target) = self.target(sv) else {
            if !self.defined(sv) {
                return Err(Die::new(
                    "Can't use an undefined value as an ARRAY reference",
                ));
            }
            return Err(self.strict_refs(sv, "an ARRAY"));
        };
        match self.image.items(target) {
            Some(items) => Ok(items.iter().map(|item| self.element(*item)).collect()),
            None => Err(Die::new("Not an ARRAY reference")),
        }
    }

    /// An array item as the scalar Perl reads: undef for an element that does not exist.
    pub fn element(&self, item: NodeId) -> Sv {
        match self.image.value(item) {
            Value::ImmortalUndef => Sv::Missing,
            _ => Sv::Node(item),
        }
    }

    /// Perl's message for a string used as a reference under `strict refs`.
    fn strict_refs(&self, sv: Sv, kind: &str) -> Die {
        let text = self
            .string(sv)
            .map_or_else(|_| Vec::new(), |text| text.bytes.into_owned());
        let mut message = b"Can't use string (\"".to_vec();
        message.extend_from_slice(&text[..text.len().min(32)]);
        message.extend_from_slice(if text.len() > 32 { b"\"..." } else { b"\"" });
        message
            .extend_from_slice(format!(") as {kind} ref while \"strict refs\" in use").as_bytes());
        Die(message)
    }

    /// Perl's truth of a scalar: false for undef, the empty string, the string `0` and
    /// a number that is zero; a reference is true.
    pub fn is_true(&self, sv: Sv) -> bool {
        match self.value(sv) {
            None | Some(Value::Undef | Value::ImmortalUndef | Value::ElementUndef | Value::No) => {
                false
            }
            Some(Value::Yes | Value::Ref { .. }) => true,
            Some(Value::Bool(value)) => value,
            Some(Value::Integer(value)) => value != 0,
            Some(Value::Double(value)) => value != 0.0,
            Some(Value::String { bytes, .. }) => {
                let bytes = self.image.bytes(bytes);
                !(bytes.is_empty() || bytes == b"0")
            }
            // A container or pattern is never a scalar's value; a reference to one is.
            Some(Value::Array(_) | Value::Hash { .. } | Value::Regexp(_)) => true,
        }
    }

    /// The scalar in numeric context. Undef is 0. A reference would be its address in
    /// memory, which no run repeats: the caller's shard fails instead.
    pub fn number(&self, sv: Sv) -> Result<Num, Die> {
        Ok(match self.value(sv) {
            None | Some(Value::Undef | Value::ImmortalUndef | Value::ElementUndef) => Num::Int(0),
            Some(Value::Yes | Value::Bool(true)) => Num::Int(1),
            Some(Value::No | Value::Bool(false)) => Num::Int(0),
            Some(Value::Integer(value)) => Num::from_i64(value),
            Some(Value::Double(value)) => Num::from_f64(value),
            Some(Value::String { bytes, .. }) => numify(self.image.bytes(bytes)),
            Some(Value::Ref { .. } | Value::Array(_) | Value::Hash { .. } | Value::Regexp(_)) => {
                return Err(reference_as_value("a number"))
            }
        })
    }

    /// The scalar in string context (`"$x"`). Undef is the empty string.
    pub fn string(&self, sv: Sv) -> Result<OwnedStr<'a>, Die> {
        Ok(match self.value(sv) {
            None | Some(Value::Undef | Value::ImmortalUndef | Value::ElementUndef) => {
                OwnedStr::ascii(Vec::new())
            }
            Some(Value::Yes | Value::Bool(true)) => OwnedStr::ascii(b"1".to_vec()),
            Some(Value::No | Value::Bool(false)) => OwnedStr::ascii(Vec::new()),
            Some(Value::Integer(value)) => {
                let mut text = Vec::new();
                write_int(i128::from(value), &mut text);
                OwnedStr::ascii(text)
            }
            Some(Value::Double(value)) => OwnedStr::ascii(double_to_perl_string(value)),
            Some(Value::String { bytes, utf8 }) => OwnedStr {
                bytes: Cow::Borrowed(self.image.bytes(bytes)),
                utf8,
            },
            Some(Value::Ref { .. } | Value::Array(_) | Value::Hash { .. } | Value::Regexp(_)) => {
                return Err(reference_as_value("a string"))
            }
        })
    }
}

fn reference_as_value(what: &str) -> Die {
    Die::new(format!(
        "a reference is used as {what}: Perl writes its address in memory, which no other run repeats"
    ))
}

/// A Perl string the conversion made or borrowed from the image.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedStr<'a> {
    pub bytes: Cow<'a, [u8]>,
    pub utf8: bool,
}

impl<'a> OwnedStr<'a> {
    pub fn ascii(bytes: Vec<u8>) -> Self {
        Self {
            bytes: Cow::Owned(bytes),
            utf8: false,
        }
    }

    pub fn as_perl(&self) -> PerlStr<'_> {
        PerlStr {
            bytes: &self.bytes,
            utf8: self.utf8,
        }
    }

    /// The string's characters as UTF-8 (see [`PerlStr::characters`]).
    pub fn characters(&self) -> Vec<u8> {
        self.as_perl().characters().into_owned()
    }
}
