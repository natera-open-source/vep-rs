// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The decoder: one pass over the bytes, one node a tagged value.
//!
//! Storable gives a tag to every value in the order it meets it and writes containers
//! depth first: a reference before its target, an array before its items, a hash before
//! its entries, each entry as the value and then the key. `SX_BLESS` and `SX_IX_BLESS`
//! carry a class for the next tagged value alone. The decoder walks the same order with
//! its own stack of open containers, so an image nested to any depth reads in constant
//! call depth.

use super::error::{ErrorKind, StorableError};
use super::header::{ByteOrder, Header};
use super::image::{HashEntry, Image, KeyFlags, Node, NodeId, Run, Span, Value, NO_CLASS};

/// The opcodes of binary major 2, by number, as Storable.xs names them.
pub mod op {
    pub const SX_OBJECT: u8 = 0;
    pub const SX_LSCALAR: u8 = 1;
    pub const SX_ARRAY: u8 = 2;
    pub const SX_HASH: u8 = 3;
    pub const SX_REF: u8 = 4;
    pub const SX_UNDEF: u8 = 5;
    pub const SX_INTEGER: u8 = 6;
    pub const SX_DOUBLE: u8 = 7;
    pub const SX_BYTE: u8 = 8;
    pub const SX_NETINT: u8 = 9;
    pub const SX_SCALAR: u8 = 10;
    pub const SX_TIED_ARRAY: u8 = 11;
    pub const SX_TIED_HASH: u8 = 12;
    pub const SX_TIED_SCALAR: u8 = 13;
    pub const SX_SV_UNDEF: u8 = 14;
    pub const SX_SV_YES: u8 = 15;
    pub const SX_SV_NO: u8 = 16;
    pub const SX_BLESS: u8 = 17;
    pub const SX_IX_BLESS: u8 = 18;
    pub const SX_HOOK: u8 = 19;
    pub const SX_OVERLOAD: u8 = 20;
    pub const SX_TIED_KEY: u8 = 21;
    pub const SX_TIED_IDX: u8 = 22;
    pub const SX_UTF8STR: u8 = 23;
    pub const SX_LUTF8STR: u8 = 24;
    pub const SX_FLAG_HASH: u8 = 25;
    pub const SX_CODE: u8 = 26;
    pub const SX_WEAKREF: u8 = 27;
    pub const SX_WEAKOVERLOAD: u8 = 28;
    pub const SX_VSTRING: u8 = 29;
    pub const SX_LVSTRING: u8 = 30;
    pub const SX_SVUNDEF_ELEM: u8 = 31;
    pub const SX_REGEXP: u8 = 32;
    pub const SX_LOBJECT: u8 = 33;
    pub const SX_BOOLEAN_TRUE: u8 = 34;
    pub const SX_BOOLEAN_FALSE: u8 = 35;

    /// The name of an opcode, `unknown` past the last one defined.
    pub fn name(opcode: u8) -> &'static str {
        const NAMES: [&str; 36] = [
            "SX_OBJECT",
            "SX_LSCALAR",
            "SX_ARRAY",
            "SX_HASH",
            "SX_REF",
            "SX_UNDEF",
            "SX_INTEGER",
            "SX_DOUBLE",
            "SX_BYTE",
            "SX_NETINT",
            "SX_SCALAR",
            "SX_TIED_ARRAY",
            "SX_TIED_HASH",
            "SX_TIED_SCALAR",
            "SX_SV_UNDEF",
            "SX_SV_YES",
            "SX_SV_NO",
            "SX_BLESS",
            "SX_IX_BLESS",
            "SX_HOOK",
            "SX_OVERLOAD",
            "SX_TIED_KEY",
            "SX_TIED_IDX",
            "SX_UTF8STR",
            "SX_LUTF8STR",
            "SX_FLAG_HASH",
            "SX_CODE",
            "SX_WEAKREF",
            "SX_WEAKOVERLOAD",
            "SX_VSTRING",
            "SX_LVSTRING",
            "SX_SVUNDEF_ELEM",
            "SX_REGEXP",
            "SX_LOBJECT",
            "SX_BOOLEAN_TRUE",
            "SX_BOOLEAN_FALSE",
        ];
        NAMES.get(opcode as usize).copied().unwrap_or("unknown")
    }
}

/// `SHV_RESTRICTED`, the one hash flag.
const SHV_RESTRICTED: u8 = 0x01;
/// `SHV_K_ISSV`: the key is a stored scalar, not bytes.
const SHV_K_ISSV: u8 = 0x08;
/// `SHR_U32_RE_LEN`: the pattern length of an `SX_REGEXP` is four bytes, not one.
const SHR_U32_RE_LEN: u8 = 0x01;

/// A container whose target, items or entries are still being read.
enum Open {
    Ref {
        node: u32,
    },
    Array {
        node: u32,
        remaining: u32,
        base: usize,
    },
    Hash {
        node: u32,
        remaining: u32,
        base: usize,
        flagged: bool,
        restricted: bool,
    },
}

struct Decoder<'a> {
    data: &'a [u8],
    pos: usize,
    order: ByteOrder,
    integer_size: usize,
    double_size: usize,
    nodes: Vec<Node>,
    items: Vec<NodeId>,
    entries: Vec<HashEntry>,
    classes: Vec<Span>,
    regexps: Vec<(Span, Span)>,
    vstrings: Vec<(NodeId, Span)>,
    open: Vec<Open>,
    open_items: Vec<NodeId>,
    open_entries: Vec<HashEntry>,
    /// The class `SX_BLESS` or `SX_IX_BLESS` named for the next tagged value.
    pending_class: u32,
    /// The version string `SX_VSTRING` or `SX_LVSTRING` carried for the next value.
    pending_vstring: Option<Span>,
}

pub(super) fn decode(data: Vec<u8>) -> Result<Image, StorableError> {
    let header = Header::parse(&data)?;
    let mut decoder = Decoder {
        data: &data,
        pos: header.body_offset,
        order: header.byte_order,
        integer_size: header.integer_size as usize,
        double_size: header
            .native
            .as_ref()
            .and_then(|native| native.double_size)
            .unwrap_or(8) as usize,
        nodes: Vec::new(),
        items: Vec::new(),
        entries: Vec::new(),
        classes: Vec::new(),
        regexps: Vec::new(),
        vstrings: Vec::new(),
        open: Vec::new(),
        open_items: Vec::new(),
        open_entries: Vec::new(),
        pending_class: NO_CLASS,
        pending_vstring: None,
    };
    let root = decoder.run()?;
    let Decoder {
        pos,
        nodes,
        items,
        entries,
        classes,
        regexps,
        vstrings,
        ..
    } = decoder;
    Ok(Image {
        data,
        header,
        nodes,
        items,
        entries,
        classes,
        regexps,
        vstrings,
        root,
        end: pos,
    })
}

impl Decoder<'_> {
    fn error(&self, kind: ErrorKind, offset: usize) -> StorableError {
        StorableError::new(kind, offset)
    }

    fn truncated(&self) -> StorableError {
        StorableError::new(ErrorKind::Truncated, self.data.len())
    }

    fn byte(&mut self) -> Result<u8, StorableError> {
        let byte = *self.data.get(self.pos).ok_or_else(|| self.truncated())?;
        self.pos += 1;
        Ok(byte)
    }

    fn take(&mut self, len: usize) -> Result<&[u8], StorableError> {
        let end = self
            .pos
            .checked_add(len)
            .filter(|end| *end <= self.data.len())
            .ok_or_else(|| self.truncated())?;
        let bytes = &self.data[self.pos..end];
        self.pos = end;
        Ok(bytes)
    }

    fn span(&mut self, len: usize) -> Result<Span, StorableError> {
        let start = self.pos;
        self.take(len)?;
        Ok(Span {
            start,
            len: len as u32,
        })
    }

    /// Four bytes, big-endian in every image: an object tag or a network integer.
    fn be_u32(&mut self) -> Result<u32, StorableError> {
        let bytes: [u8; 4] = self.take(4)?.try_into().expect("four bytes");
        Ok(u32::from_be_bytes(bytes))
    }

    /// A four-byte length in the image's own order (`RLEN`).
    fn length(&mut self) -> Result<u32, StorableError> {
        let bytes: [u8; 4] = self.take(4)?.try_into().expect("four bytes");
        Ok(match self.order {
            ByteOrder::Big => u32::from_be_bytes(bytes),
            ByteOrder::Little => u32::from_le_bytes(bytes),
        })
    }

    /// An eight-byte length in the image's own order (`READ_U64`).
    fn long_length(&mut self) -> Result<u64, StorableError> {
        let bytes: [u8; 8] = self.take(8)?.try_into().expect("eight bytes");
        Ok(match self.order {
            ByteOrder::Big => u64::from_be_bytes(bytes),
            ByteOrder::Little => u64::from_le_bytes(bytes),
        })
    }

    /// Tags a value: Storable's `SEEN`. The class named for it blesses it and the
    /// version string attaches to it.
    fn tag(&mut self, value: Value, offset: usize) -> Result<u32, StorableError> {
        let index = u32::try_from(self.nodes.len())
            .ok()
            .filter(|index| *index != u32::MAX)
            .ok_or_else(|| {
                self.error(
                    ErrorKind::BadLength("more values than the reader holds"),
                    offset,
                )
            })?;
        self.nodes.push(Node {
            value,
            class: std::mem::replace(&mut self.pending_class, NO_CLASS),
        });
        if let Some(vstring) = self.pending_vstring.take() {
            self.vstrings.push((NodeId(index), vstring));
        }
        Ok(index)
    }

    fn string(&mut self, len: u64, utf8: bool, offset: usize) -> Result<u32, StorableError> {
        let len = u32::try_from(len)
            .map_err(|_| self.error(ErrorKind::BadLength("string of 4 GiB or more"), offset))?;
        let bytes = self.span(len as usize)?;
        self.tag(Value::String { bytes, utf8 }, offset)
    }

    /// The count of an array or hash. Storable reads it as a signed 32-bit number and
    /// a negative one yields an empty container.
    fn count(&mut self) -> Result<u32, StorableError> {
        let count = self.length()? as i32;
        let count = count.max(0) as u32;
        // Every item or entry takes at least one byte: a count past the bytes left is
        // an image that ends early, refused before anything is reserved for it.
        if count as usize > self.data.len() - self.pos {
            return Err(self.truncated());
        }
        Ok(count)
    }

    fn open_reference(
        &mut self,
        weak: bool,
        overloaded: bool,
        offset: usize,
    ) -> Result<(), StorableError> {
        // The reference is tagged before its target, which may point back at it.
        let node = self.tag(
            Value::Ref {
                target: NodeId(u32::MAX),
                weak,
                overloaded,
            },
            offset,
        )?;
        self.open.push(Open::Ref { node });
        Ok(())
    }

    fn open_array(&mut self, offset: usize) -> Result<Option<u32>, StorableError> {
        let count = self.count()?;
        let node = self.tag(Value::Array(Run { start: 0, len: 0 }), offset)?;
        if count == 0 {
            return Ok(Some(node));
        }
        self.open.push(Open::Array {
            node,
            remaining: count,
            base: self.open_items.len(),
        });
        Ok(None)
    }

    fn open_hash(&mut self, flagged: bool, offset: usize) -> Result<Option<u32>, StorableError> {
        let restricted = flagged && self.byte()? & SHV_RESTRICTED != 0;
        let count = self.count()?;
        let node = self.tag(
            Value::Hash {
                entries: Run { start: 0, len: 0 },
                restricted,
            },
            offset,
        )?;
        if count == 0 {
            return Ok(Some(node));
        }
        self.open.push(Open::Hash {
            node,
            remaining: count,
            base: self.open_entries.len(),
            flagged,
            restricted,
        });
        Ok(None)
    }

    /// Reads up to the next value that is whole once its opcode's own bytes are read
    /// (a scalar, an empty container, a back-reference) and returns its tag, or `None`
    /// when the opcode opened a container.
    fn start_value(&mut self) -> Result<Option<u32>, StorableError> {
        loop {
            let offset = self.pos;
            let opcode = self.byte()?;
            let unsupported = |decoder: &Self| {
                decoder.error(
                    ErrorKind::UnsupportedOpcode {
                        opcode,
                        name: op::name(opcode),
                    },
                    offset,
                )
            };
            let tagged = match opcode {
                op::SX_OBJECT => {
                    let tag = self.be_u32()?;
                    return self.back_reference(u64::from(tag), offset).map(Some);
                }
                op::SX_BLESS => {
                    let mut len = u32::from(self.byte()?);
                    if len & 0x80 != 0 {
                        len = self.length()?;
                    }
                    let name = self.span(len as usize)?;
                    self.pending_class = self.classes.len() as u32;
                    self.classes.push(name);
                    continue;
                }
                op::SX_IX_BLESS => {
                    let mut index = u32::from(self.byte()?);
                    if index & 0x80 != 0 {
                        index = self.length()?;
                    }
                    if index as usize >= self.classes.len() {
                        return Err(self.error(ErrorKind::UnknownClass(u64::from(index)), offset));
                    }
                    self.pending_class = index;
                    continue;
                }
                op::SX_VSTRING | op::SX_LVSTRING => {
                    let len = if opcode == op::SX_VSTRING {
                        u32::from(self.byte()?)
                    } else {
                        self.length()?
                    };
                    self.pending_vstring = Some(self.span(len as usize)?);
                    continue;
                }
                op::SX_UNDEF => self.tag(Value::Undef, offset)?,
                op::SX_SV_UNDEF => self.tag(Value::ImmortalUndef, offset)?,
                op::SX_SVUNDEF_ELEM => self.tag(Value::ElementUndef, offset)?,
                op::SX_SV_YES => self.tag(Value::Yes, offset)?,
                op::SX_SV_NO => self.tag(Value::No, offset)?,
                op::SX_BOOLEAN_TRUE => self.tag(Value::Bool(true), offset)?,
                op::SX_BOOLEAN_FALSE => self.tag(Value::Bool(false), offset)?,
                op::SX_BYTE => {
                    let byte = self.byte()?;
                    self.tag(Value::Integer(i64::from(byte) - 128), offset)?
                }
                op::SX_NETINT => {
                    let value = self.be_u32()? as i32;
                    self.tag(Value::Integer(i64::from(value)), offset)?
                }
                op::SX_INTEGER => {
                    let little = self.order == ByteOrder::Little;
                    let value = match self.integer_size {
                        8 => {
                            let bytes: [u8; 8] = self.take(8)?.try_into().expect("eight bytes");
                            if little {
                                i64::from_le_bytes(bytes)
                            } else {
                                i64::from_be_bytes(bytes)
                            }
                        }
                        _ => {
                            let bytes: [u8; 4] = self.take(4)?.try_into().expect("four bytes");
                            i64::from(if little {
                                i32::from_le_bytes(bytes)
                            } else {
                                i32::from_be_bytes(bytes)
                            })
                        }
                    };
                    self.tag(Value::Integer(value), offset)?
                }
                op::SX_DOUBLE => {
                    if self.double_size != 8 {
                        return Err(
                            self.error(ErrorKind::BadLength("double that is not 8 bytes"), offset)
                        );
                    }
                    let bytes: [u8; 8] = self.take(8)?.try_into().expect("eight bytes");
                    let value = match self.order {
                        ByteOrder::Little => f64::from_le_bytes(bytes),
                        ByteOrder::Big => f64::from_be_bytes(bytes),
                    };
                    self.tag(Value::Double(value), offset)?
                }
                op::SX_SCALAR | op::SX_UTF8STR => {
                    let len = u64::from(self.byte()?);
                    self.string(len, opcode == op::SX_UTF8STR, offset)?
                }
                op::SX_LSCALAR | op::SX_LUTF8STR => {
                    let len = u64::from(self.length()?);
                    self.string(len, opcode == op::SX_LUTF8STR, offset)?
                }
                op::SX_REF | op::SX_WEAKREF | op::SX_OVERLOAD | op::SX_WEAKOVERLOAD => {
                    self.open_reference(
                        matches!(opcode, op::SX_WEAKREF | op::SX_WEAKOVERLOAD),
                        matches!(opcode, op::SX_OVERLOAD | op::SX_WEAKOVERLOAD),
                        offset,
                    )?;
                    return Ok(None);
                }
                op::SX_ARRAY => return self.open_array(offset),
                op::SX_HASH => return self.open_hash(false, offset),
                op::SX_FLAG_HASH => return self.open_hash(true, offset),
                op::SX_REGEXP => {
                    let flags = self.byte()?;
                    let pattern_len = if flags & SHR_U32_RE_LEN != 0 {
                        self.length()?
                    } else {
                        u32::from(self.byte()?)
                    };
                    let pattern = self.span(pattern_len as usize)?;
                    let flags_len = self.byte()?;
                    let flags = self.span(flags_len as usize)?;
                    self.regexps.push((pattern, flags));
                    self.tag(Value::Regexp(self.regexps.len() as u32 - 1), offset)?
                }
                op::SX_LOBJECT => {
                    // A string of 2 GiB or more, or a tag past 32 bits. The large array
                    // and hash forms hold more items than the reader's tables.
                    let kind = self.byte()?;
                    match kind {
                        op::SX_OBJECT => {
                            let tag = self.long_length()?;
                            return self.back_reference(tag, offset).map(Some);
                        }
                        op::SX_LSCALAR | op::SX_LUTF8STR => {
                            let len = self.long_length()?;
                            self.string(len, kind == op::SX_LUTF8STR, offset)?
                        }
                        _ => return Err(unsupported(self)),
                    }
                }
                _ => return Err(unsupported(self)),
            };
            return Ok(Some(tagged));
        }
    }

    /// `SX_OBJECT`: the value an earlier tag names. A class or version string named
    /// for it is dropped, as Storable drops the class.
    fn back_reference(&mut self, tag: u64, offset: usize) -> Result<u32, StorableError> {
        self.pending_class = NO_CLASS;
        self.pending_vstring = None;
        if tag >= self.nodes.len() as u64 {
            return Err(self.error(ErrorKind::UnknownTag(tag), offset));
        }
        Ok(tag as u32)
    }

    /// The key that follows a hash value: its flags in a flagged hash, a four-byte
    /// length and the bytes.
    fn key(&mut self, flagged: bool) -> Result<(Span, KeyFlags), StorableError> {
        let offset = self.pos;
        let flags = if flagged { self.byte()? } else { 0 };
        if flags & SHV_K_ISSV != 0 {
            return Err(self.error(ErrorKind::ScalarHashKey, offset));
        }
        let len = self.length()?;
        if len > i32::MAX as u32 {
            return Err(self.error(ErrorKind::BadLength("negative hash key length"), offset));
        }
        Ok((self.span(len as usize)?, KeyFlags(flags)))
    }

    fn run(&mut self) -> Result<NodeId, StorableError> {
        loop {
            let Some(mut done) = self.start_value()? else {
                continue;
            };
            // Hand the finished value to the container that waits for it, and a
            // container that this completes to its own, up to the root.
            loop {
                let Some(top) = self.open.last_mut() else {
                    return Ok(NodeId(done));
                };
                match top {
                    Open::Ref { node } => {
                        let node = *node;
                        if let Value::Ref { target, .. } = &mut self.nodes[node as usize].value {
                            *target = NodeId(done);
                        }
                        self.open.pop();
                        done = node;
                    }
                    Open::Array {
                        node,
                        remaining,
                        base,
                    } => {
                        self.open_items.push(NodeId(done));
                        *remaining -= 1;
                        if *remaining > 0 {
                            break;
                        }
                        let (node, base) = (*node, *base);
                        let start = self.items.len();
                        self.items.extend(self.open_items.drain(base..));
                        self.nodes[node as usize].value = Value::Array(Run {
                            start,
                            len: (self.items.len() - start) as u32,
                        });
                        self.open.pop();
                        done = node;
                    }
                    Open::Hash {
                        node,
                        remaining,
                        base,
                        flagged,
                        restricted,
                    } => {
                        let (node, base, flagged, restricted) =
                            (*node, *base, *flagged, *restricted);
                        *remaining -= 1;
                        let last = *remaining == 0;
                        let (key, flags) = self.key(flagged)?;
                        self.open_entries.push(HashEntry {
                            key,
                            flags,
                            value: NodeId(done),
                        });
                        if !last {
                            break;
                        }
                        let start = self.entries.len();
                        self.entries.extend(self.open_entries.drain(base..));
                        self.nodes[node as usize].value = Value::Hash {
                            entries: Run {
                                start,
                                len: (self.entries.len() - start) as u32,
                            },
                            restricted,
                        };
                        self.open.pop();
                        done = node;
                    }
                }
            }
        }
    }
}
