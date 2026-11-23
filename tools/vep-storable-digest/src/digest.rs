// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The two digests of a transcript and the line a shard's transcript becomes.
//!
//! `md5_full` is the MD5 of the transcript's `Data::Dumper` text. `md5_content` is the
//! MD5 of the same text of a copy the program makes first: the copy leaves out the keys
//! a cache build changes at every depth, holds every scalar that is all digits as a
//! number, and is unblessed throughout.

use std::borrow::Cow;
use std::collections::HashMap;

use crate::dumper::{dump, hash_entries, Indent, Key, Overlay, VId, VNode, Val};
use crate::md5::md5_hex;
use crate::number::double_to_perl_string;
use crate::perl::{Die, Perl, Sv};
use crate::shards::{transcripts, Shard};
use crate::storable::{NodeId, Value};
use crate::text::{strip_version, Text};

/// What the copy behind `md5_content` leaves out and turns into numbers.
#[derive(Debug, Clone, Copy)]
pub struct ContentRule {
    /// The keys left out at every depth.
    pub strip: &'static [&'static [u8]],
    /// Whether a scalar that is all digits becomes a number.
    pub integers: bool,
}

/// The rule of the digest: database identifiers, dates and the gene-version label go,
/// and an integer is a number however the cache stores it.
pub const CONTENT_RULE: ContentRule = ContentRule {
    strip: &[
        b"dbID",
        b"adaptor",
        b"created_date",
        b"modified_date",
        b"_gene_version",
    ],
    integers: true,
};

/// One transcript of a shard.
#[derive(Debug, Clone)]
pub struct Record {
    /// The stable id without its version.
    pub id: Text,
    /// The line's bytes, final newline included.
    pub line: Vec<u8>,
    /// `md5_content`.
    pub content: [u8; 32],
}

/// The transcripts of one shard, and what ended the shard early if anything did.
#[derive(Debug, Default)]
pub struct ShardDigest {
    pub records: Vec<Record>,
    pub failure: Option<Die>,
}

/// `/^-?\d+\z/` on bytes, with ASCII digits.
fn is_integer_text(bytes: &[u8]) -> bool {
    let digits = bytes.strip_prefix(b"-").unwrap_or(bytes);
    !digits.is_empty() && digits.iter().all(u8::is_ascii_digit)
}

/// `$text + 0` for a text of an optional minus and digits: an integer while it fits
/// 64 bits, signed or unsigned, and past that the double of the same value, which
/// Data::Dumper writes as a string.
fn integer_value<'a>(text: &[u8]) -> Val<'a> {
    let (negative, digits) = match text.strip_prefix(b"-") {
        Some(digits) => (true, digits),
        None => (false, text),
    };
    let first = digits
        .iter()
        .position(|byte| *byte != b'0')
        .unwrap_or(digits.len());
    let significant = &digits[first..];
    // Twenty digits hold every 64-bit value; a longer run is past them.
    let magnitude: Option<u128> = (significant.len() <= 20).then(|| {
        significant
            .iter()
            .fold(0u128, |sum, digit| sum * 10 + u128::from(digit - b'0'))
    });
    match magnitude {
        Some(value) if negative && value <= 1 << 63 => Val::Int((-(value as i128)) as i64),
        Some(value) if !negative && value <= i64::MAX as u128 => Val::Int(value as i64),
        Some(value) if !negative && value <= u128::from(u64::MAX) => Val::UInt(value as u64),
        _ => {
            let number: f64 = std::str::from_utf8(text)
                .ok()
                .and_then(|text| text.parse().ok())
                .unwrap_or(0.0);
            Val::Str {
                bytes: Cow::Owned(double_to_perl_string(number)),
                utf8: false,
            }
        }
    }
}

/// A hash or an array of the copy whose entries or items are still being copied.
enum Open<'i> {
    Hash {
        copy: VId,
        rest: std::vec::IntoIter<(Key<'i>, NodeId)>,
        entries: Vec<(Key<'i>, Val<'i>)>,
        /// The key whose value is being copied.
        key: Option<Key<'i>>,
    },
    Array {
        copy: VId,
        rest: std::slice::Iter<'i, NodeId>,
        items: Vec<Val<'i>>,
    },
}

/// Builds the copy `strip_copy` returns.
struct Copier<'i, 'o> {
    shard: &'i Shard,
    overlay: &'o mut Overlay<'i>,
    /// The copy of each hash and array met, so that a value reached twice is copied once.
    copies: HashMap<NodeId, VId>,
    rule: ContentRule,
    /// The hashes and arrays around the value in hand, outermost first.
    open: Vec<Open<'i>>,
}

impl<'i> Copier<'i, '_> {
    /// The copy of the scalar `id`. A hash or an array on the way waits in `open` while
    /// its entries or items are copied, so a value nested to any depth is copied in
    /// constant call depth.
    fn copy(&mut self, id: NodeId) -> Result<Val<'i>, Die> {
        let mut done = self.start(id)?;
        loop {
            let Some(top) = self.open.last_mut() else {
                return Ok(done.expect("a copy once no hash or array is open"));
            };
            // The innermost open hash or array takes the value copied for it and names
            // the next one it holds; with none left its copy is whole.
            let next = match top {
                Open::Hash {
                    copy,
                    rest,
                    entries,
                    key,
                } => {
                    if let (Some(key), Some(value)) = (key.take(), done.take()) {
                        entries.push((key, value));
                    }
                    let strip = self.rule.strip;
                    let kept = rest.find(|(name, _)| {
                        !strip
                            .iter()
                            .any(|left_out| name.characters().as_ref() == *left_out)
                    });
                    match kept {
                        Some((name, value)) => {
                            *key = Some(name);
                            Some(value)
                        }
                        None => {
                            self.overlay
                                .set(*copy, VNode::Hash(std::mem::take(entries)));
                            done = Some(Val::Ref(*copy));
                            None
                        }
                    }
                }
                Open::Array { copy, rest, items } => {
                    if let Some(value) = done.take() {
                        items.push(value);
                    }
                    match rest.next() {
                        Some(item) => Some(*item),
                        None => {
                            self.overlay.set(*copy, VNode::Array(std::mem::take(items)));
                            done = Some(Val::Ref(*copy));
                            None
                        }
                    }
                }
            };
            match next {
                Some(value) => done = self.start(value)?,
                None => {
                    self.open.pop();
                }
            }
        }
    }

    /// The copy of the scalar `id` where it is whole at once; `None` where `id` is a
    /// reference to a hash or an array met for the first time, which is put on `open`.
    fn start(&mut self, id: NodeId) -> Result<Option<Val<'i>>, Die> {
        let image = &self.shard.image;
        let value = image.value(id);
        Ok(Some(match value {
            Value::Ref { .. } if self.shard.cleared(id) => Val::Undef,
            Value::Ref { target, .. } => return self.reference(id, target),
            Value::Undef | Value::ImmortalUndef | Value::ElementUndef => Val::Undef,
            _ if !self.rule.integers => Val::CopyOf(id),
            Value::Integer(number) => Val::Int(number),
            Value::Yes | Value::Bool(true) => Val::Int(1),
            Value::Double(number) => {
                // A whole double below 10^15 reads as digits and becomes that integer.
                if is_integer_text(&double_to_perl_string(number)) {
                    Val::Int(number as i64)
                } else {
                    Val::CopyOf(id)
                }
            }
            Value::String { bytes, .. } => {
                let bytes = image.bytes(bytes);
                if is_integer_text(bytes) {
                    integer_value(bytes)
                } else {
                    Val::CopyOf(id)
                }
            }
            _ => Val::CopyOf(id),
        }))
    }

    /// The copy of a reference `id` to `target`, as [`Copier::start`] gives it.
    fn reference(&mut self, id: NodeId, target: NodeId) -> Result<Option<Val<'i>>, Die> {
        if let Some(copy) = self.copies.get(&target) {
            return Ok(Some(Val::Ref(*copy)));
        }
        let image = &self.shard.image;
        let value = image.value(target);
        let class = match image.class(target) {
            Some(b"") => Some(&b"main"[..]),
            class => class,
        };
        // `blessed` gives the class name, which is a false value for the class `0`.
        let blessed = class.is_some_and(|class| class != b"0");
        let is_hash = matches!(value, Value::Hash { .. });
        let is_array = matches!(value, Value::Array(_));
        let is_scalar = !is_hash
            && !is_array
            && !matches!(value, Value::Ref { .. } | Value::Regexp(_))
            && image.vstring(target).is_none();
        // `ref` is the class of a blessed value and the kind's name of any other.
        let named = |name: &[u8], kind: bool| match class {
            Some(class) => class == name,
            None => kind,
        };
        if named(b"HASH", is_hash) || (blessed && is_hash) {
            if !is_hash {
                return Err(Die::new("Not a HASH reference"));
            }
            let copy = self.overlay.push(VNode::Hash(Vec::new()));
            self.copies.insert(target, copy);
            self.open.push(Open::Hash {
                copy,
                rest: hash_entries(image, target).into_iter(),
                entries: Vec::new(),
                key: None,
            });
            return Ok(None);
        }
        if named(b"ARRAY", is_array) || (blessed && is_array) {
            let Some(items) = image.items(target) else {
                return Err(Die::new("Not an ARRAY reference"));
            };
            let copy = self.overlay.push(VNode::Array(Vec::new()));
            self.copies.insert(target, copy);
            self.open.push(Open::Array {
                copy,
                rest: items.iter(),
                items: Vec::with_capacity(items.len()),
            });
            return Ok(None);
        }
        if named(b"SCALAR", is_scalar) || (blessed && is_scalar) {
            // A reference to a copy of the scalar; no two share one.
            return Ok(Some(Val::Ref(
                self.overlay.push(VNode::Scalar(Val::CopyOf(target))),
            )));
        }
        // Any other reference stays the reference it is, to the value it points at.
        Ok(Some(Val::CopyOf(id)))
    }
}

/// The digest lines of one shard's transcripts under `contig`, in the shard's order.
pub fn digest_shard(shard: Option<&Shard>, contig: &[u8], rule: ContentRule) -> ShardDigest {
    let mut digest = ShardDigest::default();
    digest.failure = digest_into(shard, contig, rule, &mut digest.records).err();
    digest
}

fn digest_into(
    shard: Option<&Shard>,
    contig: &[u8],
    rule: ContentRule,
    records: &mut Vec<Record>,
) -> Result<(), Die> {
    let list = transcripts(shard)?;
    let Some(shard) = shard else {
        return Ok(());
    };
    let perl = Perl::new(&shard.image);
    let fetch = |transcript: Sv, key: &[u8]| -> Result<Sv, Die> {
        Ok(shard.read(perl.get(transcript, key)?))
    };
    let mut text = Vec::new();
    for transcript in list {
        // `next unless ref $t`.
        let Sv::Node(node) = transcript else {
            continue;
        };
        if !perl.ref_is_true(transcript) {
            continue;
        }
        // `$t->{stable_id} // $t->{_stable_id}`.
        let mut id = fetch(transcript, b"stable_id")?;
        if !perl.defined(id) {
            id = fetch(transcript, b"_stable_id")?;
        }
        if !perl.defined(id) {
            continue;
        }
        let id = perl.string(id)?;
        let mut id = Text {
            bytes: id.bytes.into_owned(),
            utf8: id.utf8,
        };
        strip_version(&mut id.bytes);

        let overlay = Overlay::default();
        text.clear();
        dump(shard, &overlay, &Val::Node(node), Indent::None, &mut text)
            .map_err(|error| Die::new(error.to_string()))?;
        let full = md5_hex(&text);

        let mut overlay = Overlay::default();
        let copy = Copier {
            shard,
            overlay: &mut overlay,
            copies: HashMap::new(),
            rule,
            open: Vec::new(),
        }
        .copy(node)?;
        text.clear();
        dump(shard, &overlay, &copy, Indent::None, &mut text)
            .map_err(|error| Die::new(error.to_string()))?;
        let content = md5_hex(&text);

        // `$t->{start} // ''`, `$t->{end} // ''`.
        let start = perl.string(fetch(transcript, b"start")?)?;
        let end = perl.string(fetch(transcript, b"end")?)?;
        let mut line = id.clone();
        line.push(b"\t", false);
        line.push(contig, false);
        line.push(b"\t", false);
        line.push(&start.bytes, start.utf8);
        line.push(b"\t", false);
        line.push(&end.bytes, end.utf8);
        line.push(b"\t", false);
        line.push(&full, false);
        line.push(b"\t", false);
        line.push(&content, false);
        let mut line = line.printed();
        line.push(b'\n');
        records.push(Record { id, line, content });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{integer_value, is_integer_text};
    use crate::dumper::Val;

    fn shown(text: &str) -> String {
        match integer_value(text.as_bytes()) {
            Val::Int(number) => format!("int {number}"),
            Val::UInt(number) => format!("uint {number}"),
            Val::Str { bytes, .. } => format!("str {}", String::from_utf8_lossy(&bytes)),
            other => format!("{other:?}"),
        }
    }

    /// `perl -e 'print $ARGV[0] + 0'` on each text, with the kind Perl's scalar has.
    #[test]
    fn digits_become_the_number_perl_adds_zero_to() {
        assert_eq!(shown("0"), "int 0");
        assert_eq!(shown("-0"), "int 0");
        assert_eq!(shown("007"), "int 7");
        assert_eq!(shown("-123456789"), "int -123456789");
        assert_eq!(shown("9223372036854775807"), "int 9223372036854775807");
        assert_eq!(shown("9223372036854775808"), "uint 9223372036854775808");
        assert_eq!(shown("18446744073709551615"), "uint 18446744073709551615");
        assert_eq!(shown("18446744073709551616"), "str 1.84467440737096e+19");
        assert_eq!(shown("-9223372036854775808"), "int -9223372036854775808");
        assert_eq!(shown("-9223372036854775809"), "str -9.22337203685478e+18");
        assert_eq!(
            shown("000000000000000000000000000000000000000000012"),
            "int 12"
        );
    }

    #[test]
    fn digits_alone_with_one_leading_minus_are_an_integer_text() {
        for text in ["0", "-0", "12", "-12", "0012"] {
            assert!(is_integer_text(text.as_bytes()), "{text}");
        }
        for text in ["", "-", "+1", "1.0", " 1", "1 ", "1e3", "--1", "1\n"] {
            assert!(!is_integer_text(text.as_bytes()), "{text:?}");
        }
    }
}
