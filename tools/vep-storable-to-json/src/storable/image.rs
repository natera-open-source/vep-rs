// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A decoded image: every stored value as a node under its Storable tag, with the bytes
//! of strings and hash keys left in the image's own buffer.

use super::header::Header;

/// One stored value. The number is the value's Storable tag: the order in which the
/// writer met it, so the root is tag 0 and an `SX_OBJECT` back-reference names its
/// target by it. Two places of an image that hold the same `NodeId` hold the same Perl
/// value: a shared reference is a shared `NodeId`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId(pub(super) u32);

impl NodeId {
    /// The Storable tag.
    pub fn tag(self) -> u32 {
        self.0
    }
}

/// A run of bytes in the image buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub(super) start: usize,
    pub(super) len: u32,
}

/// A run of array items or hash entries in the image's tables.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Run {
    pub(super) start: usize,
    pub(super) len: u32,
}

/// What a node is, as Perl sees it after `retrieve`. Blessing is not part of the value:
/// [`Image::class`] gives the class a node is blessed into.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Value {
    /// `SX_UNDEF`: an undefined scalar.
    Undef,
    /// `SX_SV_UNDEF`: Perl's immortal undef. As an array item it marks an element that
    /// does not exist.
    ImmortalUndef,
    /// `SX_SVUNDEF_ELEM`: an array element that exists and is the immortal undef.
    ElementUndef,
    /// `SX_SV_YES`: Perl's immortal true, a scalar that is at once the string `1` and
    /// the number 1.
    Yes,
    /// `SX_SV_NO`: Perl's immortal false, at once the empty string and the number 0.
    No,
    /// `SX_BOOLEAN_TRUE` or `SX_BOOLEAN_FALSE`: a boolean scalar that is not an immortal.
    Bool(bool),
    /// `SX_BYTE`, `SX_INTEGER` or `SX_NETINT`: a scalar stored as an integer.
    Integer(i64),
    /// `SX_DOUBLE`: a scalar stored as a double.
    Double(f64),
    /// `SX_SCALAR`, `SX_LSCALAR`, `SX_UTF8STR` or `SX_LUTF8STR`: a scalar stored as a
    /// string. `utf8` is Perl's UTF-8 flag: set, the bytes are UTF-8 for characters;
    /// clear, each byte is one character.
    String { bytes: Span, utf8: bool },
    /// `SX_REF`, `SX_WEAKREF`, `SX_OVERLOAD` or `SX_WEAKOVERLOAD`: a scalar that is a
    /// reference to `target`.
    Ref {
        target: NodeId,
        weak: bool,
        overloaded: bool,
    },
    /// `SX_ARRAY`: an array; [`Image::items`] gives its items in index order.
    Array(Run),
    /// `SX_HASH` or `SX_FLAG_HASH`: a hash; [`Image::entries`] gives its entries in
    /// stored order. `restricted` is the `SHV_RESTRICTED` hash flag.
    Hash { entries: Run, restricted: bool },
    /// `SX_REGEXP`: a compiled pattern; [`Image::regexp`] gives its source and flags.
    Regexp(u32),
}

/// The per-key flags of a flagged hash; none is set in a plain hash.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KeyFlags(pub u8);

impl KeyFlags {
    pub const UTF8: u8 = 0x01;
    pub const WAS_UTF8: u8 = 0x02;
    pub const LOCKED: u8 = 0x04;
    pub const PLACEHOLDER: u8 = 0x10;

    /// The key bytes are UTF-8 for characters.
    pub fn utf8(self) -> bool {
        self.0 & Self::UTF8 != 0
    }
    /// The key was UTF-8 when stored into the hash and is held as one byte a character.
    pub fn was_utf8(self) -> bool {
        self.0 & Self::WAS_UTF8 != 0
    }
    /// The value is read-only in a restricted hash.
    pub fn locked(self) -> bool {
        self.0 & Self::LOCKED != 0
    }
    /// A deleted key of a restricted hash: Perl's `keys`, `exists` and fetch do not see it.
    pub fn placeholder(self) -> bool {
        self.0 & Self::PLACEHOLDER != 0
    }
}

/// One stored key and value of a hash.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HashEntry {
    pub(super) key: Span,
    pub flags: KeyFlags,
    pub value: NodeId,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Node {
    pub(super) value: Value,
    /// Index into the class table; `NO_CLASS` for an unblessed value.
    pub(super) class: u32,
}

pub(super) const NO_CLASS: u32 = u32::MAX;

/// A whole decoded image.
#[derive(Debug)]
pub struct Image {
    pub(super) data: Vec<u8>,
    pub(super) header: Header,
    pub(super) nodes: Vec<Node>,
    pub(super) items: Vec<NodeId>,
    pub(super) entries: Vec<HashEntry>,
    pub(super) classes: Vec<Span>,
    pub(super) regexps: Vec<(Span, Span)>,
    pub(super) vstrings: Vec<(NodeId, Span)>,
    pub(super) root: NodeId,
    pub(super) end: usize,
}

impl Image {
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// The value `store` was given a reference to. Perl's `retrieve` returns a
    /// reference to it.
    pub fn root(&self) -> NodeId {
        self.root
    }

    /// The number of tagged values.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Every tagged value, in tag order.
    pub fn node_ids(&self) -> impl Iterator<Item = NodeId> {
        (0..self.nodes.len() as u32).map(NodeId)
    }

    /// Bytes of the buffer the image occupies, header included. Perl's `fd_retrieve`
    /// reads this many and leaves what follows unread.
    pub fn image_len(&self) -> usize {
        self.end
    }

    /// Bytes of the buffer after the image.
    pub fn trailing_len(&self) -> usize {
        self.data.len() - self.end
    }

    pub fn value(&self, id: NodeId) -> Value {
        self.nodes[id.0 as usize].value
    }

    /// The class the value is blessed into, as the bytes of its name.
    pub fn class(&self, id: NodeId) -> Option<&[u8]> {
        match self.nodes[id.0 as usize].class {
            NO_CLASS => None,
            index => Some(self.bytes(self.classes[index as usize])),
        }
    }

    /// The bytes of a [`Span`] a value or an entry carries.
    pub fn bytes(&self, span: Span) -> &[u8] {
        &self.data[span.start..span.start + span.len as usize]
    }

    /// The bytes and UTF-8 flag of a string scalar.
    pub fn string(&self, id: NodeId) -> Option<(&[u8], bool)> {
        match self.value(id) {
            Value::String { bytes, utf8 } => Some((self.bytes(bytes), utf8)),
            _ => None,
        }
    }

    /// The node a reference points at.
    pub fn target(&self, id: NodeId) -> Option<NodeId> {
        match self.value(id) {
            Value::Ref { target, .. } => Some(target),
            _ => None,
        }
    }

    /// The items of an array, in index order. An item whose value is
    /// [`Value::ImmortalUndef`] is an element that does not exist.
    pub fn items(&self, id: NodeId) -> Option<&[NodeId]> {
        match self.value(id) {
            Value::Array(run) => Some(&self.items[run.start..run.start + run.len as usize]),
            _ => None,
        }
    }

    /// The entries of a hash, in stored order, placeholders included. A key stored
    /// twice appears twice; Perl keeps the later value.
    pub fn entries(&self, id: NodeId) -> Option<&[HashEntry]> {
        match self.value(id) {
            Value::Hash { entries, .. } => {
                Some(&self.entries[entries.start..entries.start + entries.len as usize])
            }
            _ => None,
        }
    }

    /// The bytes of an entry's key.
    pub fn key(&self, entry: &HashEntry) -> &[u8] {
        self.bytes(entry.key)
    }

    /// The value a hash holds under `key`, as a Perl fetch sees it: the last entry whose
    /// key bytes equal `key`, placeholders skipped. For a key of ASCII bytes this is
    /// Perl's own equality whatever the entry's UTF-8 flag.
    pub fn fetch(&self, id: NodeId, key: &[u8]) -> Option<NodeId> {
        self.entries(id)?
            .iter()
            .rev()
            .find(|entry| !entry.flags.placeholder() && self.key(entry) == key)
            .map(|entry| entry.value)
    }

    /// The source and the flag letters of a compiled pattern.
    pub fn regexp(&self, id: NodeId) -> Option<(&[u8], &[u8])> {
        match self.value(id) {
            Value::Regexp(index) => {
                let (pattern, flags) = self.regexps[index as usize];
                Some((self.bytes(pattern), self.bytes(flags)))
            }
            _ => None,
        }
    }

    /// The version-string magic a scalar carries (`v1.2.3` holds bytes 1, 2, 3 and this
    /// is the literal's text).
    pub fn vstring(&self, id: NodeId) -> Option<&[u8]> {
        self.vstrings
            .iter()
            .find(|(node, _)| *node == id)
            .map(|(_, span)| self.bytes(*span))
    }

    /// Perl's `defined` on the value.
    pub fn is_defined(&self, id: NodeId) -> bool {
        !matches!(
            self.value(id),
            Value::Undef | Value::ImmortalUndef | Value::ElementUndef
        )
    }
}
