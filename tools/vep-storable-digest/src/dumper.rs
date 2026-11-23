// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The text `Data::Dumper::Dumper` returns for one value, as the module's compiled
//! implementation (`Dumpxs` of Data::Dumper 2.179) builds it under
//! `$Data::Dumper::Sortkeys = 1` and `$Data::Dumper::Terse = 1`, with
//! `$Data::Dumper::Indent` 0 or 1 and every other setting at its default.
//!
//! The rules, each from `DD_dump`:
//!
//! - A hash is `{'key' => value,...}` with its keys in the order of Perl's `cmp`; an
//!   array is `[item,...]`. Under `Indent = 1` every key and item opens a line indented
//!   two spaces a level, the closing bracket has a line of its own, and the text ends
//!   with a newline; under `Indent = 0` there is no white space but the ` => `.
//! - A blessed value is `bless( <value>, '<class>' )`; a compiled pattern blessed into
//!   `Regexp` is `qr/<source>/<flags>` alone.
//! - A scalar stored as an integer is its digits, in quotes when they are more than ten
//!   characters; a string is in single quotes with `'` and `\` escaped, even when it is
//!   all digits; a double is the quoted text Perl gives it; undef is `undef`. A string
//!   with Perl's UTF-8 flag that holds a character above 127 is in double quotes with
//!   such characters as `\x{<hex>}` and `"`, `\`, `$` and `@` escaped.
//! - A reference to a scalar is `\<scalar>`, a blessed one `do{\(my $o = <scalar>)}`.
//! - A value met a second time is written as the path of its first place
//!   (`$VAR1->{'exons'}[0]`): a hash or an array by the path of the reference that led
//!   to it, a scalar that a reference also points at as `\<path>` or `${<path>}`.
//!
//! The value is a [`Val`]: a scalar of the image, or one of an [`Overlay`], the values a
//! caller builds beside the image (a copy without some keys, a summary in the place of
//! a field).

use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::HashMap;

use crate::number::{double_to_perl_string, write_int};
use crate::shards::Shard;
use crate::storable::{HashEntry, Image, NodeId, Value};

/// `$Data::Dumper::Maxrecurse`: the depth at which the module gives up.
const MAX_RECURSE: usize = 1000;

/// A value of an [`Overlay`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VId(u32);

/// A scalar handed to the dumper.
#[derive(Debug, Clone)]
pub enum Val<'a> {
    /// The image's own scalar. Met twice, it is the same scalar.
    Node(NodeId),
    /// A scalar of its own that holds a copy of the image's scalar: the same number or
    /// string, or a reference to the same target.
    CopyOf(NodeId),
    Undef,
    Int(i64),
    /// An integer above the signed range.
    UInt(u64),
    Str {
        bytes: Cow<'a, [u8]>,
        utf8: bool,
    },
    /// A reference to a value of the overlay.
    Ref(VId),
    /// A reference of its own to a value of the image: what `retrieve` returns for the
    /// image's root.
    RefTo(NodeId),
}

/// A hash key as the scalar Perl's hash iterator returns for it.
#[derive(Debug, Clone)]
pub struct Key<'a> {
    pub pv: Cow<'a, [u8]>,
    pub utf8: bool,
}

impl<'a> Key<'a> {
    pub fn ascii(name: &'a [u8]) -> Self {
        Key {
            pv: Cow::Borrowed(name),
            utf8: false,
        }
    }

    /// The key of a stored entry. A key that was UTF-8 when it went into the hash comes
    /// back with the flag, whatever its characters.
    pub fn of_entry(image: &'a Image, entry: &HashEntry) -> Self {
        let bytes = image.key(entry);
        if entry.flags.utf8() {
            Key {
                pv: Cow::Borrowed(bytes),
                utf8: true,
            }
        } else if entry.flags.was_utf8() {
            Key {
                pv: latin1_as_utf8(bytes),
                utf8: true,
            }
        } else {
            Key {
                pv: Cow::Borrowed(bytes),
                utf8: false,
            }
        }
    }

    /// The key's characters as UTF-8: what Perl's `cmp` and a hash lookup compare.
    pub fn characters(&self) -> Cow<'_, [u8]> {
        if self.utf8 {
            Cow::Borrowed(&self.pv)
        } else {
            latin1_as_utf8(&self.pv)
        }
    }
}

/// Bytes that are one character each, as the UTF-8 of those characters.
fn latin1_as_utf8(bytes: &[u8]) -> Cow<'_, [u8]> {
    if bytes.is_ascii() {
        return Cow::Borrowed(bytes);
    }
    let mut out = Vec::with_capacity(bytes.len() + 8);
    for &byte in bytes {
        if byte < 0x80 {
            out.push(byte);
        } else {
            out.push(0xc0 | (byte >> 6));
            out.push(0x80 | (byte & 0x3f));
        }
    }
    Cow::Owned(out)
}

/// What a reference of the overlay points at.
#[derive(Debug, Clone)]
pub enum VNode<'a> {
    /// An unblessed hash.
    Hash(Vec<(Key<'a>, Val<'a>)>),
    /// An unblessed array.
    Array(Vec<Val<'a>>),
    /// A scalar that is not a reference.
    Scalar(Val<'a>),
}

/// Values a caller builds beside the image.
#[derive(Debug, Default)]
pub struct Overlay<'a> {
    nodes: Vec<VNode<'a>>,
}

impl<'a> Overlay<'a> {
    pub fn push(&mut self, node: VNode<'a>) -> VId {
        self.nodes.push(node);
        VId(self.nodes.len() as u32 - 1)
    }

    pub fn set(&mut self, id: VId, node: VNode<'a>) {
        self.nodes[id.0 as usize] = node;
    }

    pub fn clear(&mut self) {
        self.nodes.clear();
    }
}

/// The entries of an image's hash as Perl's `keys` and a fetch see them, in the order
/// the image stores them: placeholders skipped, and a key stored twice once, at its
/// first place and with its later value. Perl's own order of a hash is not this one
/// and changes from run to run.
pub fn hash_entries_in_stored_order(image: &Image, hash: NodeId) -> Vec<(Key<'_>, NodeId)> {
    let entries = image.entries(hash).unwrap_or(&[]);
    let mut live: Vec<(Key<'_>, NodeId)> = Vec::with_capacity(entries.len());
    let mut place: HashMap<Vec<u8>, usize> = HashMap::new();
    for entry in entries.iter().filter(|entry| !entry.flags.placeholder()) {
        let key = Key::of_entry(image, entry);
        match place.entry(key.characters().into_owned()) {
            std::collections::hash_map::Entry::Occupied(at) => live[*at.get()].1 = entry.value,
            std::collections::hash_map::Entry::Vacant(at) => {
                at.insert(live.len());
                live.push((key, entry.value));
            }
        }
    }
    live
}

/// The same entries in the order of Perl's `cmp` on their keys: the order
/// `$Data::Dumper::Sortkeys = 1` writes them in.
pub fn hash_entries(image: &Image, hash: NodeId) -> Vec<(Key<'_>, NodeId)> {
    let entries = image.entries(hash).unwrap_or(&[]);
    let mut live: Vec<(Key<'_>, NodeId)> = entries
        .iter()
        .filter(|entry| !entry.flags.placeholder())
        .map(|entry| (Key::of_entry(image, entry), entry.value))
        .collect();
    live.sort_by(|left, right| key_order(&left.0, &right.0));
    // The sort is stable, so equal keys stand in stored order: the first keeps its
    // place and takes the last one's value.
    live.dedup_by(|later, earlier| {
        let same = key_order(&later.0, &earlier.0) == Ordering::Equal;
        if same {
            earlier.1 = later.1;
        }
        same
    });
    live
}

/// Why the module would not return a text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DumpError {
    /// `Recursion limit of 1000 exceeded`.
    TooDeep,
    /// The image holds an array or a hash where a scalar belongs.
    ContainerAsScalar,
}

impl std::fmt::Display for DumpError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DumpError::TooDeep => write!(formatter, "Recursion limit of {MAX_RECURSE} exceeded"),
            DumpError::ContainerAsScalar => {
                formatter.write_str("the image holds an array or a hash where a scalar belongs")
            }
        }
    }
}

/// What the module's table of seen values is keyed on: the address of a value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Addr {
    Node(NodeId),
    Virtual(VId),
    /// Perl's three immortal scalars: every place that holds one holds the same scalar.
    Undef,
    Yes,
    No,
}

struct Seen {
    name: Vec<u8>,
    /// Set once the value's text is complete; a scalar is named only after that.
    marked: bool,
}

struct Dumper<'i, 'v, 'o> {
    shard: &'i Shard,
    overlay: &'v Overlay<'i>,
    pretty: bool,
    seen: HashMap<Addr, Seen>,
    out: &'o mut Vec<u8>,
    failure: Option<DumpError>,
}

/// `$Data::Dumper::Indent`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Indent {
    /// 0: one line, no white space.
    None,
    /// 1: a line a key or item, two spaces a level.
    Fixed,
}

/// Appends the module's text for `val` to `out`.
pub fn dump<'i>(
    shard: &'i Shard,
    overlay: &Overlay<'i>,
    val: &Val<'i>,
    indent: Indent,
    out: &mut Vec<u8>,
) -> Result<(), DumpError> {
    let mut dumper = Dumper {
        shard,
        overlay,
        pretty: indent == Indent::Fixed,
        seen: HashMap::new(),
        out,
        failure: None,
    };
    let mut name = b"$VAR1".to_vec();
    dumper.value(val, &mut name, 0);
    if dumper.pretty {
        dumper.out.push(b'\n');
    }
    match dumper.failure {
        Some(failure) => Err(failure),
        None => Ok(()),
    }
}

/// Whether a child's path joins its parent's with `->`: after the top name, and after a
/// dereferenced scalar (`${...}`).
fn needs_arrow(name: &[u8]) -> bool {
    let last = name.last().copied();
    (last.is_some() && last != Some(b']') && last != Some(b'}'))
        || (name.len() > 4 && (name[1] == b'{' || (name[0] == b'\\' && name[2] == b'{')))
}

impl<'i> Dumper<'i, '_, '_> {
    fn image(&self) -> &'i Image {
        &self.shard.image
    }

    fn put(&mut self, text: &[u8]) {
        self.out.extend_from_slice(text);
    }

    /// The separator and the indentation that open a line at `level`.
    fn line(&mut self, level: usize) {
        if self.pretty {
            self.out.push(b'\n');
            for _ in 0..level {
                self.out.extend_from_slice(b"  ");
            }
        }
    }

    fn node_addr(&self, id: NodeId) -> Addr {
        match self.image().value(id) {
            Value::ImmortalUndef | Value::ElementUndef => Addr::Undef,
            Value::Yes => Addr::Yes,
            Value::No => Addr::No,
            _ => Addr::Node(id),
        }
    }

    fn value(&mut self, val: &Val<'i>, name: &mut Vec<u8>, level: usize) {
        if self.failure.is_some() {
            return;
        }
        match val {
            Val::Node(id) => self.node(*id, true, name, level),
            Val::CopyOf(id) => self.node(*id, false, name, level),
            Val::Undef => self.put(b"undef"),
            Val::Int(number) => self.integer(i128::from(*number)),
            Val::UInt(number) => self.integer(i128::from(*number)),
            Val::Str { bytes, utf8 } => self.string(bytes, *utf8),
            Val::Ref(id) => self.reference(Addr::Virtual(*id), name, level),
            Val::RefTo(id) => {
                let addr = self.node_addr(*id);
                self.reference(addr, name, level);
            }
        }
    }

    /// A scalar of the image: itself (`own`) or a copy of it.
    fn node(&mut self, id: NodeId, own: bool, name: &mut Vec<u8>, level: usize) {
        let value = self.image().value(id);
        if let Value::Ref { target, .. } = value {
            if self.shard.cleared(id) {
                // A weak reference whose target is gone is an undefined scalar.
                return self.plain(own.then_some(Addr::Node(id)), name, |dumper| {
                    dumper.put(b"undef");
                });
            }
            let addr = self.node_addr(target);
            return self.reference(addr, name, level);
        }
        if matches!(
            value,
            Value::Array(_) | Value::Hash { .. } | Value::Regexp(_)
        ) {
            self.failure = Some(DumpError::ContainerAsScalar);
            return;
        }
        // An immortal is the same scalar wherever it stands, copies excepted.
        let addr = match self.node_addr(id) {
            Addr::Node(_) if !own => None,
            Addr::Yes | Addr::No if !own => None,
            addr => Some(addr),
        };
        self.plain(addr, name, |dumper| dumper.scalar(id, value));
    }

    /// A scalar that is no reference, with the module's bookkeeping around its text.
    fn plain(&mut self, addr: Option<Addr>, name: &[u8], text: impl FnOnce(&mut Self)) {
        if let Some(addr) = addr {
            match self.seen.get(&addr) {
                Some(seen) if seen.marked => {
                    let path = seen.name.clone();
                    self.put(b"${");
                    self.put(&path);
                    self.put(b"}");
                    return;
                }
                Some(_) => {}
                // The immortal undef is looked up and never entered.
                None if addr == Addr::Undef => {}
                None => {
                    let mut path = Vec::with_capacity(name.len() + 1);
                    path.push(b'\\');
                    path.extend_from_slice(name);
                    self.seen.insert(
                        addr,
                        Seen {
                            name: path,
                            marked: false,
                        },
                    );
                }
            }
        }
        text(self);
        if let Some(seen) = addr.and_then(|addr| self.seen.get_mut(&addr)) {
            seen.marked = true;
        }
    }

    fn scalar(&mut self, id: NodeId, value: Value) {
        match value {
            Value::Undef | Value::ImmortalUndef | Value::ElementUndef => self.put(b"undef"),
            Value::Yes | Value::Bool(true) => self.put(b"1"),
            Value::No | Value::Bool(false) => self.put(b"''"),
            Value::Integer(number) => self.integer(i128::from(number)),
            Value::Double(number) => {
                let text = double_to_perl_string(number);
                self.string(&text, false);
            }
            Value::String { bytes, utf8 } => {
                let image = self.image();
                match image.vstring(id) {
                    // A version string is written as its literal.
                    Some(literal) => self.put(literal),
                    None => self.string(image.bytes(bytes), utf8),
                }
            }
            Value::Ref { .. } | Value::Array(_) | Value::Hash { .. } | Value::Regexp(_) => {}
        }
    }

    /// An integer: bare up to ten characters, in quotes past that.
    fn integer(&mut self, number: i128) {
        let mut digits = Vec::with_capacity(24);
        write_int(number, &mut digits);
        if digits.len() > 10 {
            self.out.push(b'\'');
            self.put(&digits);
            self.out.push(b'\'');
        } else {
            self.put(&digits);
        }
    }

    fn string(&mut self, bytes: &[u8], utf8: bool) {
        quote(bytes, utf8, self.out);
    }

    /// A reference whose target has the address `addr`.
    fn reference(&mut self, addr: Addr, name: &mut Vec<u8>, level: usize) {
        if let Some(seen) = self.seen.get(&addr) {
            let path = seen.name.clone();
            self.put(&path);
            return;
        }
        self.seen.insert(
            addr,
            Seen {
                name: name.clone(),
                marked: false,
            },
        );
        if level >= MAX_RECURSE {
            self.failure = Some(DumpError::TooDeep);
            return;
        }
        match addr {
            Addr::Virtual(id) => self.virtual_target(id, name, level),
            Addr::Node(target) => self.node_target(target, name, level),
            // A reference to an immortal: the target is found through any node that
            // holds it.
            Addr::Undef => {
                self.put(b"\\");
                self.put(b"undef");
            }
            Addr::Yes => self.put(b"\\1"),
            Addr::No => self.put(b"\\''"),
        }
        if let Some(seen) = self.seen.get_mut(&addr) {
            seen.marked = true;
        }
    }

    fn virtual_target(&mut self, id: VId, name: &mut Vec<u8>, level: usize) {
        let overlay = self.overlay;
        match &overlay.nodes[id.0 as usize] {
            VNode::Hash(entries) => {
                let mut order: Vec<&(Key<'i>, Val<'i>)> = entries.iter().collect();
                order.sort_by(|left, right| key_order(&left.0, &right.0));
                self.hash(order.iter().map(|entry| (&entry.0, &entry.1)), name, level);
            }
            VNode::Array(items) => self.array(items.iter(), items.len(), name, level),
            VNode::Scalar(val) => {
                let mut inner = dereferenced(name);
                self.put(b"\\");
                self.value(val, &mut inner, level + 1);
            }
        }
    }

    fn node_target(&mut self, target: NodeId, name: &mut Vec<u8>, level: usize) {
        let image = self.image();
        let value = image.value(target);
        let pattern = matches!(value, Value::Regexp(_));
        let class = match image.class(target) {
            // Blessing into the empty name blesses into `main`.
            Some(b"") => Some(&b"main"[..]),
            Some(class) => Some(class),
            // A compiled pattern is an object of `Regexp` unless it is blessed elsewhere.
            None if pattern => Some(&b"Regexp"[..]),
            None => None,
        };
        let blessed = class.filter(|class| !(pattern && *class == b"Regexp"));
        if blessed.is_some() {
            self.put(b"bless( ");
        }
        match value {
            Value::Regexp(_) => self.pattern(target),
            Value::Array(_) => {
                let items: Vec<Val<'i>> = image
                    .items(target)
                    .unwrap_or(&[])
                    .iter()
                    .map(|item| Val::Node(*item))
                    .collect();
                self.array(items.iter(), items.len(), name, level);
            }
            Value::Hash { .. } => {
                let entries: Vec<(Key<'i>, Val<'i>)> = hash_entries(image, target)
                    .into_iter()
                    .map(|(key, value)| (key, Val::Node(value)))
                    .collect();
                self.hash(
                    entries.iter().map(|entry| (&entry.0, &entry.1)),
                    name,
                    level,
                );
            }
            _ => {
                // A reference to a scalar, which may itself be a reference.
                let mut inner = dereferenced(name);
                if class.is_some() {
                    self.put(b"do{\\(my $o = ");
                    self.value(&Val::Node(target), &mut inner, level + 1);
                    self.put(b")}");
                } else {
                    self.put(b"\\");
                    self.value(&Val::Node(target), &mut inner, level + 1);
                }
            }
        }
        if let Some(class) = blessed {
            self.put(b", '");
            escape_single(class, self.out);
            self.put(b"' )");
        }
    }

    fn array<'x>(
        &mut self,
        items: impl Iterator<Item = &'x Val<'i>>,
        count: usize,
        name: &mut Vec<u8>,
        level: usize,
    ) where
        'i: 'x,
    {
        self.put(b"[");
        let base = name.len();
        if needs_arrow(name) {
            name.extend_from_slice(b"->");
        }
        name.push(b'[');
        let stem = name.len();
        for (index, item) in items.enumerate() {
            name.truncate(stem);
            write_int(index as i128, name);
            name.push(b']');
            self.line(level + 1);
            self.value(item, name, level + 1);
            if index + 1 < count {
                self.put(b",");
            }
        }
        if count > 0 {
            self.line(level);
        }
        name.truncate(base);
        self.put(b"]");
    }

    /// A hash whose entries come in the order of their keys.
    fn hash<'x>(
        &mut self,
        entries: impl Iterator<Item = (&'x Key<'i>, &'x Val<'i>)>,
        name: &mut Vec<u8>,
        level: usize,
    ) where
        'i: 'x,
    {
        self.put(b"{");
        let base = name.len();
        if needs_arrow(name) {
            name.extend_from_slice(b"->");
        }
        name.push(b'{');
        let stem = name.len();
        let mut any = false;
        for (key, val) in entries {
            if any {
                self.put(b",");
            }
            any = true;
            self.line(level + 1);
            // The key's quoted text is also its step in the path.
            let at = self.out.len();
            quote(&key.pv, key.utf8, self.out);
            name.truncate(stem);
            name.extend_from_slice(&self.out[at..]);
            name.push(b'}');
            self.put(b" => ");
            self.value(val, name, level + 1);
        }
        if any {
            self.line(level);
        }
        name.truncate(base);
        self.put(b"}");
    }

    /// `qr/<source>/<flags>`, with every `/` of the source that no backslash precedes
    /// escaped.
    fn pattern(&mut self, id: NodeId) {
        let (source, flags) = self.image().regexp(id).unwrap_or((b"", b""));
        self.put(b"qr/");
        let mut from = 0;
        let mut at = 0;
        while at < source.len() {
            match source[at] {
                b'\\' => at += 1,
                b'/' => {
                    self.put(&source[from..at]);
                    self.put(b"\\/");
                    from = at + 1;
                }
                _ => {}
            }
            at += 1;
        }
        self.put(&source[from.min(source.len())..]);
        self.put(b"/");
        self.put(flags);
    }
}

/// The path of the scalar a reference at `name` points at: `${<name>}`.
fn dereferenced(name: &[u8]) -> Vec<u8> {
    let mut inner = Vec::with_capacity(name.len() + 3);
    inner.extend_from_slice(b"${");
    inner.extend_from_slice(name);
    inner.push(b'}');
    inner
}

/// The bytes with `'` and `\` escaped by a backslash.
fn escape_single(bytes: &[u8], out: &mut Vec<u8>) {
    for &byte in bytes {
        if byte == b'\'' || byte == b'\\' {
            out.push(b'\\');
        }
        out.push(byte);
    }
}

/// A string or a key in quotes. Without Perl's UTF-8 flag, and with it when every
/// byte is ASCII, single quotes; else double quotes with each character above 127 as
/// `\x{<hex>}`.
pub fn quote(bytes: &[u8], utf8: bool, out: &mut Vec<u8>) {
    if !utf8 || bytes.is_ascii() {
        out.push(b'\'');
        escape_single(bytes, out);
        out.push(b'\'');
        return;
    }
    out.push(b'"');
    let mut at = 0;
    while at < bytes.len() {
        let byte = bytes[at];
        if byte >= 0x80 {
            // A sequence that is not well formed is one byte long and the character 0.
            let (code, length) = decode_character(&bytes[at..]).unwrap_or((0, 1));
            out.extend_from_slice(format!("\\x{{{code:x}}}").as_bytes());
            at += length;
            continue;
        }
        if (0x20..0x7f).contains(&byte) && matches!(byte, b'"' | b'\\' | b'$' | b'@') {
            out.push(b'\\');
        }
        out.push(byte);
        at += 1;
    }
    out.push(b'"');
}

/// One character of Perl's UTF-8 at the start of `bytes`, which opens with a byte above
/// 127: the code point and the bytes it takes. Perl's form runs past Unicode's (up to
/// seven bytes here) and takes surrogates; a sequence that is short, overlong or opens
/// with a continuation byte is `None`.
fn decode_character(bytes: &[u8]) -> Option<(u64, usize)> {
    let lead = bytes[0];
    let (length, mut code) = match lead {
        0xc0..=0xdf => (2, u64::from(lead & 0x1f)),
        0xe0..=0xef => (3, u64::from(lead & 0x0f)),
        0xf0..=0xf7 => (4, u64::from(lead & 0x07)),
        0xf8..=0xfb => (5, u64::from(lead & 0x03)),
        0xfc..=0xfd => (6, u64::from(lead & 0x01)),
        0xfe => (7, 0),
        _ => return None,
    };
    if bytes.len() < length {
        return None;
    }
    for &byte in &bytes[1..length] {
        if byte & 0xc0 != 0x80 {
            return None;
        }
        code = (code << 6) | u64::from(byte & 0x3f);
    }
    // The shortest form alone is well formed.
    let shortest = match code {
        0..=0x7f => 1,
        0x80..=0x7ff => 2,
        0x800..=0xffff => 3,
        0x1_0000..=0x1f_ffff => 4,
        0x20_0000..=0x3ff_ffff => 5,
        0x400_0000..=0x7fff_ffff => 6,
        _ => 7,
    };
    (shortest == length).then_some((code, length))
}

/// The order of Perl's `cmp` on two keys: by character.
pub fn key_order(left: &Key<'_>, right: &Key<'_>) -> Ordering {
    if left.utf8 == right.utf8 {
        // The same encoding on both sides: bytes order as characters do.
        return left.pv.cmp(&right.pv);
    }
    left.characters().cmp(&right.characters())
}
