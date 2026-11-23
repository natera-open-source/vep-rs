// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The JSON format as `parse_json` (`crate::parse::json`) reads it: the key of a
//! line's input record, its entries, whether two lines parse into the same entries, and
//! the counts the comparison adds for a kept-back line.
//!
//! A line is one object. Its record entry is its top-level members other than `input`,
//! `vep_rs_version` and the four consequence lists; each item of a consequence list that is
//! an array is one more entry. The comparison reads a member's value as one text
//! (`json_scalar`): a string itself, an array its items' texts joined by `,`, anything else
//! its compact JSON with an object's members sorted.
//!
//! A record group is every object of a side whose `input` text opens with the same five
//! tab-separated columns, the key of a VCF line (`super::vcf`): the lines of one key name
//! the input's lines of that key in order, and a line past them is labelled by its place
//! in the output. A line of white space alone is no line of a group.
//!
//! The canonical form of a line is its key, the record entry's members sorted by name, and
//! each consequence list's items in order with their members sorted by name; a value's form
//! is an object's members sorted by name, an array's items in order, a string's unescaped
//! text, a number's value, or the literal. A number's value is what the comparison's reader
//! (`pyjson::loads`) makes of its text: a whole number, kept as its digits, where the text
//! has neither a fraction nor an exponent, and a 64-bit float otherwise, so `2.5e-05` and
//! `0.000025` are one form and `1` and `1.0` are two. Every part is tagged and every text
//! carries its length, so two lines share a form only when the comparison parses them into
//! the same entries. Lines of two forms may still parse into the same entries (`"a,b"` and
//! `["a","b"]`, `1e400` and `-1e400`): the form is the finer comparison.

use super::grouping::{Grouping, Scanned, Side, Tally};
use super::vcf::{self, InputKeys, LISTED_TEXT_SHAPE, OPENING_BRACKETS};
use super::Plan;
use crate::error::Error;
use memchr::{memchr2, memchr2_iter, memchr_iter, memmem};
use std::borrow::Cow;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::ops::Range;
use std::sync::Arc;

const MEMBER_MARK: &str = "a variant_allele or feature id value holds | or #";

/// What the comparison names as the groups it must see, beside the input's keys.
struct Force {
    input: InputKeys,
    floors: Floors,
    texts: Vec<String>,
    /// The kinds of consequence list that are not kept, as bits of `Object::kinds`.
    dropped: u8,
}

pub(crate) struct Json {
    force: Arc<Force>,
}

impl Json {
    pub(crate) fn new(plan: &Plan, block: usize) -> crate::error::Result<Json> {
        let dropped = match &plan.keep_values {
            Some((column, _)) if column != "kind" => {
                return Err(Error::comparator(
                    "the column of the kept values of the json format is kind",
                ));
            }
            Some((_, kept)) => {
                let kinds = KINDS.iter().enumerate();
                kinds
                    .filter(|(_, kind)| !kept.iter().any(|value| value == *kind))
                    .fold(0, |bits, (i, _)| bits | 1 << i)
            }
            None => 0,
        };
        let floors = plan.floor_fields.iter().map(|(column, fields)| Floor {
            column: column.clone(),
            fields: fields.clone(),
        });
        Ok(Json {
            force: Arc::new(Force {
                input: vcf::input_keys(plan, "json", block)?,
                floors: Floors::new(floors.collect()),
                texts: plan.force_text.clone(),
                dropped,
            }),
        })
    }
}

/// The shapes of a line under which every line is handed over, in the order one is
/// reported.
const LINE_SHAPES: [&str; 7] = [
    NOT_ONE_OBJECT,
    NAMED_TWICE,
    LONE_SURROGATE,
    ITEM_NOT_AN_OBJECT,
    TOO_LONG,
    MEMBER_MARK,
    LISTED_TEXT_SHAPE,
];

/// The bit of the shape `shape` among `LINE_SHAPES`, the first's for a shape not listed.
fn shape_bit(shape: &str) -> u32 {
    1 << LINE_SHAPES
        .iter()
        .position(|listed| *listed == shape)
        .unwrap_or(0)
}

impl Grouping for Json {
    type Header = ();
    type Side = JsonSide;

    const SHAPES: &'static [&'static str] = &LINE_SHAPES;

    const PLACED: Option<u64> = Some(vcf::PLACED);

    fn is_header(_line: &[u8]) -> bool {
        false
    }

    fn is_blank(line: &[u8]) -> bool {
        line.iter().all(|&c| is_space(c))
    }

    fn header_line(_header: &mut (), _line: &[u8]) {}

    fn side(&self, _header: (), _data: bool, reference: bool) -> JsonSide {
        JsonSide {
            force: Arc::clone(&self.force),
            // The lists name reference lines alone.
            dropped: if reference { self.force.dropped } else { 0 },
            texts: vcf::listed_texts(&self.force.texts, reference),
        }
    }

    fn between(&self, _expected: &JsonSide, _actual: &JsonSide) -> Option<u32> {
        None
    }

    fn limit(&self, key: u64) -> Option<u64> {
        Some(self.force.input.limit(key))
    }

    /// A floor that finds none of its fields as a member of an entry is looked for
    /// inside the members of the consequence items.
    fn strict(&self, carried: u64) -> bool {
        let floors = &self.force.floors;
        let listed = floors.floors().len();
        let all = if listed < 64 {
            (1u64 << listed) - 1
        } else {
            u64::MAX
        };
        listed > 64 || carried & all != all || floors.reads_object_members()
    }
}

pub(crate) struct JsonSide {
    force: Arc<Force>,
    /// On the reference side, the kinds of consequence list that are not kept.
    dropped: u8,
    /// On the reference side, the listed texts.
    texts: Vec<memmem::Finder<'static>>,
}

impl Side for JsonSide {
    type Memo = ();

    fn shape(&self, _header_after_data: bool) -> Option<u32> {
        None
    }

    fn scan(&self, _memo: &mut (), text: &[u8]) -> Scanned {
        let object = match scan(text, &self.force.floors, self.dropped) {
            Line::Object(object) => object,
            Line::Malformed(shape) => return unread(shape),
            Line::Blank => return unread(NOT_ONE_OBJECT),
        };
        let key = self.force.input.token(&object.key);
        // A line of many consequence items is flat: what the comparison's interpreter
        // may stop at is how deep the line's objects and arrays nest, not how many they
        // are. A `\u` escape can spell a character outside ASCII, a digit of a long run
        // or a bracket on a line whose own bytes show none.
        let fragile = vcf::unconverted(text)
            || object.depth as usize > OPENING_BRACKETS
            || object.unicode_escape;
        let mut scanned = Scanned {
            key,
            hash: object.form,
            forced: fragile
                || object.kinds & self.dropped != 0
                || self.force.input.forced(key)
                || self.texts.iter().any(|listed| listed.find(text).is_some()),
            fragile,
            shapes: 0,
            strict: object.opens_object || object.dotted_member || object.custom_annotations,
            carried: object.carried,
        };
        if object.key_mark {
            scanned.shapes |= shape_bit(MEMBER_MARK);
        }
        // A rewrite of the key's text would move the line to another input record.
        if self
            .texts
            .iter()
            .any(|listed| listed.find(&object.key).is_some())
        {
            scanned.shapes |= shape_bit(LISTED_TEXT_SHAPE);
        }
        scanned
    }

    fn entries(&self, text: &[u8]) -> u32 {
        match scan(text, &self.force.floors, 0) {
            Line::Object(object) => object.entries,
            Line::Malformed(_) | Line::Blank => 1,
        }
    }

    /// A kept-back line is held against the other side's by its canonical form: two
    /// lines of two texts may be the same entries.
    fn kept_text(&self, line: &[u8], text: &mut Vec<u8>) {
        form_text(line, text);
    }

    fn populated(
        &self,
        block: &[u8],
        kept: &[Range<usize>],
        tally: &mut Tally,
        entries: &mut Vec<u32>,
    ) {
        let floors = &self.force.floors;
        let mut counts = Counts::new();
        let mut under = vec![0u64; floors.floors().len()];
        for line in kept {
            entries.push(count(&block[line.clone()], &mut counts, floors, &mut under));
        }
        let columns = floors.floors().iter().map(|floor| floor.column.as_bytes());
        for (counted, more) in [
            (&mut tally.populated, counts.populated()),
            (&mut tally.present, counts.present()),
            (
                &mut tally.under,
                columns
                    .zip(under)
                    .filter(|(_, count)| *count > 0)
                    .map(|(column, count)| (column.to_vec(), count))
                    .collect(),
            ),
        ] {
            for (name, count) in more {
                *counted.entry(name).or_insert(0) += count;
            }
        }
    }
}

/// The scan of a line that is not read as one object.
fn unread(shape: &str) -> Scanned {
    Scanned {
        key: 0,
        hash: 0,
        forced: true,
        fragile: true,
        shapes: shape_bit(shape),
        strict: false,
        carried: 0,
    }
}

/// The kinds of the four consequence lists, in the order of the bits of `Object::kinds`.
pub(crate) const KINDS: [&str; 4] = [
    "transcript",
    "intergenic",
    "regulatory_feature",
    "motif_feature",
];
const LISTS: [&[u8]; 4] = [
    b"transcript_consequences",
    b"intergenic_consequences",
    b"regulatory_feature_consequences",
    b"motif_feature_consequences",
];
const INPUT: &[u8] = b"input";
const IDENTITY: &[u8] = b"vep_rs_version";
const CUSTOM_ANNOTATIONS: &[u8] = b"custom_annotations";
/// The members of a consequence item the comparison builds the entry's key from.
const KEY_MEMBERS: [&[u8]; 4] = [
    b"variant_allele",
    b"transcript_id",
    b"regulatory_feature_id",
    b"motif_feature_id",
];

/// The shapes of a line that is not read: `Line::Malformed` carries one of them.
pub(crate) const NOT_ONE_OBJECT: &str = "a line that is not one JSON object";
pub(crate) const NAMED_TWICE: &str = "a JSON object naming a member twice";
pub(crate) const LONE_SURROGATE: &str = "a JSON escape of a lone surrogate";
pub(crate) const ITEM_NOT_AN_OBJECT: &str = "a consequence list item that is not an object";
pub(crate) const TOO_LONG: &str = "a line of 4 GiB or more";

pub(crate) enum Line<'a> {
    /// Empty, or white space alone: no line of a group.
    Blank,
    /// Not one object the comparison reads as this module does; the text names the shape.
    Malformed(&'static str),
    Object(Object<'a>),
}

pub(crate) struct Object<'a> {
    /// The first five tab-separated columns of the unescaped `input` text; empty without a
    /// string `input`.
    pub(crate) key: Cow<'a, [u8]>,
    /// 1 and the items of the consequence lists that are arrays.
    pub(crate) entries: u32,
    /// A digest of the canonical form: equal forms, equal digests.
    pub(crate) form: u128,
    /// Bit `i`: the consequence list of `KINDS[i]` is a non-empty array.
    pub(crate) kinds: u8,
    /// Bit `i`: an entry of the line carries a field of floor `i` as a member: the record's
    /// entry, or a consequence item of a list `scan` was not told to leave out. A floor
    /// past the 64th is never reported.
    pub(crate) carried: u64,
    /// A consequence item has a member whose value's text opens with `{`.
    pub(crate) opens_object: bool,
    /// A consequence item has a member whose name holds `.`.
    pub(crate) dotted_member: bool,
    /// The object has a `custom_annotations` member.
    pub(crate) custom_annotations: bool,
    /// A consequence item's `variant_allele`, `transcript_id`, `regulatory_feature_id` or
    /// `motif_feature_id` text holds `|` or `#`.
    pub(crate) key_mark: bool,
    /// The line holds a `\u` escape, in a name or a value.
    pub(crate) unicode_escape: bool,
    /// The objects and arrays open at once at the line's deepest place, the line's own
    /// object the first: what the comparison's interpreter nests its calls by.
    pub(crate) depth: u32,
}

/// A floor of the comparison: the column it is asked on and the member names it is counted
/// under.
pub(crate) struct Floor {
    pub(crate) column: String,
    pub(crate) fields: Vec<String>,
}

/// Names by their bytes, each with the index it was first given.
#[derive(Default)]
struct Names {
    names: Vec<Box<[u8]>>,
    /// Open addressing over a power of two: 0 for a free slot, else the name's index and 1.
    slots: Vec<u32>,
}

fn name_hash(name: &[u8]) -> u64 {
    const K: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut h = K ^ name.len() as u64;
    let mut words = name.chunks_exact(8);
    for word in &mut words {
        let word = u64::from_le_bytes(word.try_into().expect("eight bytes"));
        h = (h ^ word).wrapping_mul(K).rotate_left(27);
    }
    let rest = words.remainder();
    if !rest.is_empty() {
        let mut word = [0u8; 8];
        word[..rest.len()].copy_from_slice(rest);
        h = (h ^ u64::from_le_bytes(word)).wrapping_mul(K);
    }
    h ^ (h >> 32)
}

impl Names {
    fn find(&self, name: &[u8]) -> Option<usize> {
        if self.slots.is_empty() {
            return None;
        }
        let mask = self.slots.len() - 1;
        let mut slot = name_hash(name) as usize & mask;
        loop {
            match self.slots[slot] {
                0 => return None,
                n if *self.names[n as usize - 1] == *name => return Some(n as usize - 1),
                _ => slot = (slot + 1) & mask,
            }
        }
    }

    /// The index of `name`, the next index for a name that has none.
    fn index(&mut self, name: &[u8]) -> usize {
        if let Some(index) = self.find(name) {
            return index;
        }
        self.names.push(name.into());
        if self.names.len() * 2 > self.slots.len() {
            self.slots = vec![0; (self.names.len() * 4).next_power_of_two().max(64)];
            for index in 0..self.names.len() {
                self.place(index);
            }
        } else {
            self.place(self.names.len() - 1);
        }
        self.names.len() - 1
    }

    fn place(&mut self, index: usize) {
        let mask = self.slots.len() - 1;
        let mut slot = name_hash(&self.names[index]) as usize & mask;
        while self.slots[slot] != 0 {
            slot = (slot + 1) & mask;
        }
        self.slots[slot] = index as u32 + 1;
    }
}

/// The floors of a comparison, indexed for the member names of a line.
#[derive(Default)]
pub(crate) struct Floors {
    floors: Vec<Floor>,
    fields: Names,
    /// Per field, the floors that list it: as bits for the first 64, and all of them.
    bits: Vec<u64>,
    listing: Vec<Vec<u32>>,
    /// Bit `n`: a field is `n` bytes long, modulo 64.
    lengths: u64,
}

impl Floors {
    pub(crate) fn new(floors: Vec<Floor>) -> Floors {
        let mut indexed = Floors::default();
        for (n, floor) in floors.iter().enumerate() {
            for field in &floor.fields {
                let field = field.as_bytes();
                let index = indexed.fields.index(field);
                if index == indexed.bits.len() {
                    indexed.bits.push(0);
                    indexed.listing.push(Vec::new());
                }
                if n < 64 {
                    indexed.bits[index] |= 1 << n;
                }
                if indexed.listing[index].last() != Some(&(n as u32)) {
                    indexed.listing[index].push(n as u32);
                }
                indexed.lengths |= 1 << (field.len() & 63);
            }
        }
        indexed.floors = floors;
        indexed
    }

    pub(crate) fn floors(&self) -> &[Floor] {
        &self.floors
    }

    /// Whether a field holds `.`: the comparison reads such a field, where no member has
    /// that name, as a member of an object under the name before the first `.`, which
    /// `count` does not read. `count` is then exact on the lines whose `opens_object` is
    /// false alone, where the comparison finds no object either.
    pub(crate) fn reads_object_members(&self) -> bool {
        self.fields.names.iter().any(|name| name.contains(&b'.'))
    }

    fn field(&self, name: &[u8]) -> Option<usize> {
        if self.lengths & (1 << (name.len() & 63)) == 0 {
            return None;
        }
        self.fields.find(name)
    }
}

/// The forms' tags. A value opens with one of the first eight, so the end of an array and
/// the member and end of an object are told from a value.
const NULL: u8 = b'n';
const TRUE: u8 = b't';
const FALSE: u8 = b'f';
/// A whole number: its length and its digits follow.
const WHOLE: u8 = b'#';
/// A float: the eight bytes of its bits follow.
const FLOAT: u8 = b'.';
const STRING: u8 = b'"';
const ARRAY: u8 = b'[';
const OBJECT: u8 = b'{';
const ARRAY_END: u8 = b']';
const MEMBER: u8 = b'm';
const OBJECT_END: u8 = b'}';
/// The bytes of a tag and a length before a text.
const HEAD: usize = 5;

/// What an object or an array is to the line.
const OTHER: u8 = 0;
const TOP: u8 = 1;
const LIST: u8 = 2;
const ITEM: u8 = 3;

/// What a top-level member is to the record entry.
const REGULAR: u8 = 0;
const INPUT_MEMBER: u8 = 1;
const IDENTITY_MEMBER: u8 = 2;
const LIST_MEMBER: u8 = 3;

fn top_role(name: &[u8]) -> u8 {
    if name == INPUT {
        INPUT_MEMBER
    } else if name == IDENTITY {
        IDENTITY_MEMBER
    } else {
        LISTS
            .iter()
            .position(|list| *list == name)
            .map_or(REGULAR, |i| LIST_MEMBER + i as u8)
    }
}

#[derive(Clone, Copy)]
struct Member {
    /// The member's form in `out`, and its unescaped name there.
    start: usize,
    end: usize,
    name: usize,
    name_len: usize,
    /// The name's first eight bytes: two names stand as these numbers do, and as their
    /// bytes do where the numbers are one.
    head: u64,
    role: u8,
}

fn name_head(name: &[u8]) -> u64 {
    let mut head = [0u8; 8];
    let len = name.len().min(8);
    head[..len].copy_from_slice(&name[..len]);
    u64::from_be_bytes(head)
}

#[derive(Clone, Copy)]
struct Frame {
    object: bool,
    role: u8,
    /// The consequence list of a `LIST` array.
    list: u8,
    /// Whether the members read so far stand in the order of their names.
    sorted: bool,
    /// The first of the object's members in `members`.
    base: usize,
}

#[derive(Default)]
struct Work {
    /// The forms in the order the line writes its members.
    out: Vec<u8>,
    /// An object's forms while its members are put in order.
    tmp: Vec<u8>,
    members: Vec<Member>,
    frames: Vec<Frame>,
}

#[derive(Default)]
struct Scratch {
    work: Work,
    form: Vec<u8>,
}

thread_local! {
    static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch::default());
}

#[derive(Default)]
struct Facts {
    /// The `input` text in the line, where it is a string without an escape.
    input_raw: Option<(usize, usize)>,
    entries: u32,
    kinds: u8,
    carried: u64,
    opens_object: bool,
    dotted_member: bool,
    custom_annotations: bool,
    key_mark: bool,
    unicode_escape: bool,
    depth: u32,
}

/// A text's first five tab-separated columns: the comparison's `vcf_line_key`.
fn vcf_line_key(text: &[u8]) -> &[u8] {
    match memchr_iter(b'\t', text).nth(4) {
        Some(tab) => &text[..tab],
        None => text,
    }
}

fn skip_ws(b: &[u8], mut p: usize) -> usize {
    while let Some(b' ' | b'\t' | b'\n' | b'\r') = b.get(p) {
        p += 1;
    }
    p
}

/// The first byte at or after `i` that ends a run of a string's plain text: `"`, `\` or a
/// control character; the text's length without one.
fn stop(b: &[u8], mut i: usize) -> usize {
    const LOW: u64 = 0x0101_0101_0101_0101;
    const HIGH: u64 = 0x8080_8080_8080_8080;
    while i + 8 <= b.len() {
        let word = u64::from_le_bytes(b[i..i + 8].try_into().expect("eight bytes"));
        let quote = word ^ (LOW * u64::from(b'"'));
        let slash = word ^ (LOW * u64::from(b'\\'));
        // Each term marks its first such byte exactly; a mark above it may be spurious.
        let marks = (quote.wrapping_sub(LOW) & !quote)
            | (slash.wrapping_sub(LOW) & !slash)
            | (word.wrapping_sub(LOW * 0x20) & !word);
        if marks & HIGH != 0 {
            return i + ((marks & HIGH).trailing_zeros() / 8) as usize;
        }
        i += 8;
    }
    while i < b.len() && b[i] != b'"' && b[i] != b'\\' && b[i] >= 0x20 {
        i += 1;
    }
    i
}

fn hex4(b: &[u8], at: usize) -> Option<u32> {
    let digits = b.get(at..at + 4)?;
    let mut unit = 0;
    for &d in digits {
        unit = unit * 16 + char::from(d).to_digit(16)?;
    }
    Some(unit)
}

/// The `\u` escapes read so far: how many, and how many of them spell `|` or `#`.
#[derive(Clone, Copy, Default)]
struct Escapes {
    unicode: u32,
    marks: u32,
}

/// Appends the text of the escape whose `\` stands at `i` and returns the place after it.
fn escape(
    b: &[u8],
    i: usize,
    out: &mut Vec<u8>,
    seen: &mut Escapes,
) -> Result<usize, &'static str> {
    let simple = match b.get(i + 1) {
        Some(b'"') => b'"',
        Some(b'\\') => b'\\',
        Some(b'/') => b'/',
        Some(b'b') => 0x08,
        Some(b'f') => 0x0c,
        Some(b'n') => b'\n',
        Some(b'r') => b'\r',
        Some(b't') => b'\t',
        Some(b'u') => {
            let unit = hex4(b, i + 2).ok_or(NOT_ONE_OBJECT)?;
            let mut next = i + 6;
            let code = match unit {
                0xd800..=0xdbff => {
                    if b.get(next..next + 2) != Some(b"\\u") {
                        return Err(LONE_SURROGATE);
                    }
                    let low = hex4(b, next + 2).ok_or(NOT_ONE_OBJECT)?;
                    if !(0xdc00..=0xdfff).contains(&low) {
                        return Err(LONE_SURROGATE);
                    }
                    next += 6;
                    0x10000 + ((unit - 0xd800) << 10) + (low - 0xdc00)
                }
                0xdc00..=0xdfff => return Err(LONE_SURROGATE),
                _ => unit,
            };
            let c = char::from_u32(code).ok_or(NOT_ONE_OBJECT)?;
            seen.unicode += 1;
            seen.marks += u32::from(c == '|' || c == '#');
            out.extend_from_slice(c.encode_utf8(&mut [0; 4]).as_bytes());
            return Ok(next);
        }
        _ => return Err(NOT_ONE_OBJECT),
    };
    out.push(simple);
    Ok(i + 2)
}

/// Appends the unescaped text of the string whose opening quote stands at `quote` and
/// returns the place of its closing quote and whether the text is the line's own bytes.
fn string(
    b: &[u8],
    quote: usize,
    out: &mut Vec<u8>,
    seen: &mut Escapes,
) -> Result<(usize, bool), &'static str> {
    let mut i = quote + 1;
    let mut run = i;
    let mut plain = true;
    loop {
        i = stop(b, i);
        match b.get(i) {
            Some(b'"') => {
                out.extend_from_slice(&b[run..i]);
                return Ok((i, plain));
            }
            Some(b'\\') => {
                out.extend_from_slice(&b[run..i]);
                plain = false;
                i = escape(b, i, out, seen)?;
                run = i;
            }
            _ => return Err(NOT_ONE_OBJECT),
        }
    }
}

/// The end of the number that starts at `p`, by the grammar of RFC 8259.
fn number(b: &[u8], p: usize) -> Option<usize> {
    let digits = |mut i: usize| {
        let start = i;
        while b.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
        }
        (i > start).then_some(i)
    };
    let mut i = p + usize::from(b.get(p) == Some(&b'-'));
    i = match b.get(i) {
        Some(b'0') => i + 1,
        _ => digits(i)?,
    };
    if b.get(i) == Some(&b'.') {
        i = digits(i + 1)?;
    }
    if let Some(b'e' | b'E') = b.get(i) {
        i += 1;
        i = digits(i + usize::from(matches!(b.get(i), Some(b'+' | b'-'))))?;
    }
    Some(i)
}

struct Parser<'a, 'w> {
    b: &'a [u8],
    p: usize,
    floors: &'w Floors,
    /// The consequence lists, as bits of `Object::kinds`, whose items carry no floor.
    uncounted: u8,
    out: &'w mut Vec<u8>,
    tmp: &'w mut Vec<u8>,
    members: &'w mut Vec<Member>,
    frames: &'w mut Vec<Frame>,
    form: &'w mut Vec<u8>,
    facts: Facts,
    /// The form of each consequence list's value in `out`.
    lists: [Option<(usize, usize)>; 4],
    /// The unescaped text of a string `input` in `out`.
    input: Option<(usize, usize)>,
    /// The last string read, in the line, where it holds no escape.
    plain: Option<(usize, usize)>,
    escapes: Escapes,
    /// Whether the value about to be read opens the text of a consequence item's member.
    opens: bool,
    /// The consequence item's member being read: where its value starts, the marks before
    /// it, and whether it is one the entry's key is built from.
    item_value: usize,
    item_marks: u32,
    item_key: bool,
}

impl Parser<'_, '_> {
    fn ws(&mut self) {
        self.p = skip_ws(self.b, self.p);
    }

    /// Reads the string at `p` as a text of the form: `tag`, its length, its bytes.
    /// Returns where the bytes start in `out`.
    fn text(&mut self, tag: u8) -> Result<usize, &'static str> {
        let head = self.out.len();
        self.out.push(tag);
        self.out.extend_from_slice(&[0; 4]);
        let (end, plain) = string(self.b, self.p, self.out, &mut self.escapes)?;
        self.plain = plain.then_some((self.p + 1, end));
        self.p = end + 1;
        let at = head + HEAD;
        let len = (self.out.len() - at) as u32;
        self.out[head + 1..at].copy_from_slice(&len.to_le_bytes());
        Ok(at)
    }

    /// Reads a member's name and its `:` at `p`, and opens the member.
    fn name(&mut self) -> Result<(), &'static str> {
        let start = self.out.len();
        let at = self.text(MEMBER)?;
        self.ws();
        if self.b.get(self.p) != Some(&b':') {
            return Err(NOT_ONE_OBJECT);
        }
        self.p += 1;
        self.ws();
        let frame = self.frames.last_mut().expect("the open object");
        let name = &self.out[at..];
        let head = name_head(name);
        if let Some(last) = self.members[frame.base..].last() {
            if last.head > head
                || (last.head == head && self.out[last.name..last.name + last.name_len] >= *name)
            {
                frame.sorted = false;
            }
        }
        let mut role = REGULAR;
        match frame.role {
            TOP => {
                role = top_role(name);
                if role == REGULAR {
                    self.facts.custom_annotations |= name == CUSTOM_ANNOTATIONS;
                    if let Some(field) = self.floors.field(name) {
                        self.facts.carried |= self.floors.bits[field];
                    }
                }
            }
            ITEM => {
                self.facts.dotted_member |= name.contains(&b'.');
                if self.uncounted & (1 << frame.list) == 0 {
                    if let Some(field) = self.floors.field(name) {
                        self.facts.carried |= self.floors.bits[field];
                    }
                }
                self.item_key = KEY_MEMBERS.contains(&name);
                self.item_value = self.p;
                self.item_marks = self.escapes.marks;
                self.opens = true;
            }
            _ => {}
        }
        self.members.push(Member {
            start,
            end: start,
            name: at,
            name_len: self.out.len() - at,
            head,
            role,
        });
        Ok(())
    }

    /// Closes the member whose value ended before `end` in the line.
    fn close_member(&mut self, role: u8, end: usize) {
        let member = self.members.last_mut().expect("the open member");
        member.end = self.out.len();
        let value = member.name + member.name_len;
        match role {
            ITEM => {
                if self.item_key
                    && (self.escapes.marks != self.item_marks
                        || memchr2(b'|', b'#', &self.b[self.item_value..end]).is_some())
                {
                    self.facts.key_mark = true;
                }
            }
            TOP if member.role == INPUT_MEMBER => {
                if self.out[value] == STRING {
                    self.input = Some((value + HEAD, member.end));
                    self.facts.input_raw = self.plain;
                }
            }
            TOP if member.role >= LIST_MEMBER => {
                self.lists[usize::from(member.role - LIST_MEMBER)] = Some((value, member.end));
            }
            _ => {}
        }
    }

    /// Closes the open object: its members in the order of their names, none named twice.
    fn close_object(&mut self) -> Result<(), &'static str> {
        let frame = self.frames.pop().expect("the open object");
        let first = self.members.get(frame.base).map_or(0, |m| m.start);
        if !frame.sorted {
            let out = &*self.out;
            let name = |m: &Member| &out[m.name..m.name + m.name_len];
            let members = &mut self.members[frame.base..];
            members.sort_unstable_by(|x, y| x.head.cmp(&y.head).then_with(|| name(x).cmp(name(y))));
            let twice = |w: &[Member]| w[0].head == w[1].head && name(&w[0]) == name(&w[1]);
            if members.windows(2).any(twice) {
                return Err(NAMED_TWICE);
            }
        }
        let members = &self.members[frame.base..];
        if frame.role == TOP {
            let key = self
                .input
                .map_or(&[][..], |(at, end)| vcf_line_key(&self.out[at..end]));
            self.form.clear();
            self.form
                .extend_from_slice(&(key.len() as u32).to_le_bytes());
            self.form.extend_from_slice(key);
            for m in members.iter().filter(|m| m.role == REGULAR) {
                self.form.extend_from_slice(&self.out[m.start..m.end]);
            }
            self.form.push(OBJECT_END);
            for list in self.lists {
                match list {
                    Some((at, end)) if self.out[at] == ARRAY => {
                        self.form.extend_from_slice(&self.out[at..end]);
                    }
                    _ => self.form.extend_from_slice(&[ARRAY, ARRAY_END]),
                }
            }
        } else {
            if !frame.sorted {
                self.tmp.clear();
                self.tmp.extend_from_slice(&self.out[first..]);
                self.out.truncate(first);
                for m in members {
                    self.out
                        .extend_from_slice(&self.tmp[m.start - first..m.end - first]);
                }
            }
            self.out.push(OBJECT_END);
        }
        self.members.truncate(frame.base);
        Ok(())
    }

    fn literal(&mut self, word: &[u8], tag: u8) -> Result<(), &'static str> {
        if !self.b[self.p..].starts_with(word) {
            return Err(NOT_ONE_OBJECT);
        }
        self.p += word.len();
        self.out.push(tag);
        Ok(())
    }

    /// Appends the form of the number that stands from `p` to before `end`: the value the
    /// comparison's reader makes of the token. A token with neither a fraction nor an
    /// exponent is a whole number of any length, whose form is its digits with its sign,
    /// `-0` being `0`; any other token is the float its text parses to, a token past the
    /// range of a float an infinity, and its form is that float's bits, so two tokens the
    /// reader rounds to one float are one form. Every such token of the grammar parses; one
    /// that did not would keep its text.
    fn number_form(&mut self, end: usize) {
        let token = &self.b[self.p..end];
        let float = token
            .iter()
            .any(|c| matches!(c, b'.' | b'e' | b'E'))
            .then(|| std::str::from_utf8(token).ok()?.parse::<f64>().ok())
            .flatten();
        match float {
            Some(float) => {
                self.out.push(FLOAT);
                self.out.extend_from_slice(&float.to_bits().to_le_bytes());
            }
            None => {
                let digits: &[u8] = if token == b"-0" { b"0" } else { token };
                self.out.push(WHOLE);
                self.out
                    .extend_from_slice(&(digits.len() as u32).to_le_bytes());
                self.out.extend_from_slice(digits);
            }
        }
    }

    /// Reads the line as one object and leaves its canonical form in `form`.
    fn object(mut self) -> Result<Facts, &'static str> {
        self.ws();
        if self.b.get(self.p) != Some(&b'{') {
            return Err(NOT_ONE_OBJECT);
        }
        self.facts.entries = 1;
        loop {
            // A value starts at `p`.
            let c = *self.b.get(self.p).ok_or(NOT_ONE_OBJECT)?;
            let parent = self.frames.last().copied();
            if parent.is_some_and(|f| f.role == LIST) && c != b'{' {
                return Err(ITEM_NOT_AN_OBJECT);
            }
            let mut opening_string = false;
            if self.opens && c != b'[' {
                self.opens = false;
                self.facts.opens_object |= c == b'{';
                opening_string = c == b'"';
            }
            match c {
                b'{' => {
                    let role = match parent {
                        None => TOP,
                        Some(f) if f.role == LIST => {
                            self.facts.entries += 1;
                            self.facts.kinds |= 1 << f.list;
                            ITEM
                        }
                        Some(_) => OTHER,
                    };
                    self.out.push(OBJECT);
                    self.frames.push(Frame {
                        object: true,
                        role,
                        list: parent.map_or(0, |f| f.list),
                        sorted: true,
                        base: self.members.len(),
                    });
                    self.facts.depth = self.facts.depth.max(self.frames.len() as u32);
                    self.p += 1;
                    self.ws();
                    match self.b.get(self.p) {
                        Some(b'"') => {
                            self.name()?;
                            continue;
                        }
                        Some(b'}') => {
                            self.p += 1;
                            self.close_object()?;
                        }
                        _ => return Err(NOT_ONE_OBJECT),
                    }
                }
                b'[' => {
                    let list = match (parent, self.members.last()) {
                        (Some(f), Some(m)) if f.role == TOP && m.role >= LIST_MEMBER => {
                            Some(m.role - LIST_MEMBER)
                        }
                        _ => None,
                    };
                    self.out.push(ARRAY);
                    // An array without an item is a level all the same.
                    self.facts.depth = self.facts.depth.max(self.frames.len() as u32 + 1);
                    self.p += 1;
                    self.ws();
                    if self.b.get(self.p) == Some(&b']') {
                        self.p += 1;
                        self.out.push(ARRAY_END);
                        self.opens = false;
                    } else {
                        self.frames.push(Frame {
                            object: false,
                            role: if list.is_some() { LIST } else { OTHER },
                            list: list.unwrap_or(0),
                            sorted: true,
                            base: 0,
                        });
                        continue;
                    }
                }
                b'"' => {
                    let at = self.text(STRING)?;
                    if opening_string && self.out.get(at) == Some(&b'{') {
                        self.facts.opens_object = true;
                    }
                }
                b'-' | b'0'..=b'9' => {
                    let end = number(self.b, self.p).ok_or(NOT_ONE_OBJECT)?;
                    self.number_form(end);
                    self.p = end;
                }
                b't' => self.literal(b"true", TRUE)?,
                b'f' => self.literal(b"false", FALSE)?,
                b'n' => self.literal(b"null", NULL)?,
                _ => return Err(NOT_ONE_OBJECT),
            }
            // A value ended at `p`: close what ends with it.
            loop {
                let Some(frame) = self.frames.last().copied() else {
                    self.ws();
                    self.facts.unicode_escape = self.escapes.unicode > 0;
                    return if self.p == self.b.len() {
                        Ok(self.facts)
                    } else {
                        Err(NOT_ONE_OBJECT)
                    };
                };
                let end = self.p;
                self.ws();
                let c = *self.b.get(self.p).ok_or(NOT_ONE_OBJECT)?;
                self.p += 1;
                if frame.object {
                    self.close_member(frame.role, end);
                    match c {
                        b',' => {
                            self.ws();
                            if self.b.get(self.p) != Some(&b'"') {
                                return Err(NOT_ONE_OBJECT);
                            }
                            self.name()?;
                            break;
                        }
                        b'}' => self.close_object()?,
                        _ => return Err(NOT_ONE_OBJECT),
                    }
                } else {
                    match c {
                        b',' => {
                            self.ws();
                            break;
                        }
                        b']' => {
                            self.out.push(ARRAY_END);
                            self.frames.pop();
                        }
                        _ => return Err(NOT_ONE_OBJECT),
                    }
                }
            }
        }
    }
}

/// Reads `line` as one object, leaves its canonical form in `form` and returns the rest.
fn canonical(
    line: &[u8],
    floors: &Floors,
    uncounted: u8,
    work: &mut Work,
    form: &mut Vec<u8>,
) -> Result<Facts, &'static str> {
    if u32::try_from(line.len()).is_err() {
        return Err(TOO_LONG);
    }
    work.out.clear();
    work.members.clear();
    work.frames.clear();
    Parser {
        b: line,
        p: 0,
        floors,
        uncounted,
        out: &mut work.out,
        tmp: &mut work.tmp,
        members: &mut work.members,
        frames: &mut work.frames,
        form,
        facts: Facts::default(),
        lists: [None; 4],
        input: None,
        plain: None,
        escapes: Escapes::default(),
        opens: false,
        item_value: 0,
        item_marks: 0,
        item_key: false,
    }
    .object()
}

/// A digest of a form's bytes: 16 bytes a step through two multiplications, each folded.
fn digest(bytes: &[u8]) -> u128 {
    const K: [u64; 4] = [
        0xa076_1d64_78bd_642f,
        0xe703_7ed1_a0b4_28db,
        0x8ebc_6af0_9c88_c6e3,
        0x5899_65cc_7537_4cc3,
    ];
    let fold = |a: u64, b: u64| {
        let product = u128::from(a) * u128::from(b);
        (product as u64) ^ ((product >> 64) as u64)
    };
    let word = |bytes: &[u8]| u64::from_le_bytes(bytes.try_into().expect("eight bytes"));
    let len = bytes.len() as u64;
    let (mut x, mut y) = (K[0] ^ len, K[1] ^ len);
    let mut steps = bytes.chunks_exact(16);
    for step in &mut steps {
        let (a, b) = (word(&step[..8]), word(&step[8..]));
        x = fold(a ^ K[2], b ^ x);
        y = fold(b ^ K[3], a ^ y);
    }
    let mut last = [0u8; 16];
    last[..steps.remainder().len()].copy_from_slice(steps.remainder());
    let (a, b) = (word(&last[..8]), word(&last[8..]));
    x = fold(a ^ K[2], b ^ x);
    y = fold(b ^ K[3], a ^ y);
    (u128::from(fold(x ^ K[1], y ^ K[2])) << 64) | u128::from(fold(y ^ K[0], x ^ K[3]))
}

/// Whether the comparison's `str.strip` takes the byte off a line's ends.
fn is_space(c: u8) -> bool {
    matches!(c, b'\t'..=b'\r' | 0x1c..=0x1f | b' ')
}

/// Reads one line. `floors` may be empty; `uncounted` names, as bits of `Object::kinds`,
/// the consequence lists whose items are left out of `Object::carried`.
///
/// The line is valid UTF-8. A line of white space outside ASCII alone, which the comparison
/// takes for a blank one, is `Malformed` here.
pub(crate) fn scan<'a>(line: &'a [u8], floors: &Floors, uncounted: u8) -> Line<'a> {
    if line.first() != Some(&b'{') && line.iter().all(|&c| is_space(c)) {
        return Line::Blank;
    }
    SCRATCH.with_borrow_mut(|scratch| {
        let facts = match canonical(
            line,
            floors,
            uncounted,
            &mut scratch.work,
            &mut scratch.form,
        ) {
            Ok(facts) => facts,
            Err(shape) => return Line::Malformed(shape),
        };
        let form = &scratch.form;
        let key_len = u32::from_le_bytes(form[..4].try_into().expect("four bytes")) as usize;
        let key = match facts.input_raw {
            // A string without an escape holds no tab: the text is the key.
            Some((start, end)) => Cow::Borrowed(&line[start..end]),
            None => Cow::Owned(form[4..4 + key_len].to_vec()),
        };
        Line::Object(Object {
            key,
            entries: facts.entries,
            form: digest(form),
            kinds: facts.kinds,
            carried: facts.carried,
            opens_object: facts.opens_object,
            dotted_member: facts.dotted_member,
            custom_annotations: facts.custom_annotations,
            key_mark: facts.key_mark,
            unicode_escape: facts.unicode_escape,
            depth: facts.depth,
        })
    })
}

/// What stands before a line's canonical form, and before the text of a line that has
/// none, in `form_text`.
const FORM: u8 = b'F';
const NO_FORM: u8 = b'R';
/// The byte `form_text` writes before the stand-in of a byte it may not write.
const STUFFED: u8 = 0x01;

/// Appends one line of text, ended by `\n`, that two lines `scan` reads as objects share
/// exactly when their canonical forms are the same: the form with every line end and
/// every `STUFFED` byte written as two bytes. A line that is not read as an object is
/// its own text, written the same way.
pub(crate) fn form_text(line: &[u8], text: &mut Vec<u8>) {
    SCRATCH.with_borrow_mut(|scratch| {
        let Scratch { work, form } = scratch;
        let floors = Floors::default();
        let (tag, bytes) = match canonical(line, &floors, 0, work, form) {
            Ok(_) => (FORM, form.as_slice()),
            Err(_) => (NO_FORM, line),
        };
        text.push(tag);
        let mut start = 0;
        for at in memchr2_iter(b'\n', STUFFED, bytes) {
            text.extend_from_slice(&bytes[start..at]);
            text.extend_from_slice(&[STUFFED, if bytes[at] == b'\n' { 2 } else { 3 }]);
            start = at + 1;
        }
        text.extend_from_slice(&bytes[start..]);
        text.push(b'\n');
    });
}

/// Whether two lines that `scan` reads as objects parse into the same entries: their
/// canonical forms compared in full, as `form_text` writes them.
#[cfg(test)]
fn alike(expected: &[u8], actual: &[u8]) -> bool {
    let text = |line: &[u8]| {
        let mut text = Vec::new();
        form_text(line, &mut text);
        text
    };
    text(expected) == text(actual)
}

/// Per member name, the entries that carry the member and those that carry it with a value
/// other than the empty text and `-`.
#[derive(Default)]
pub(crate) struct Counts {
    names: Names,
    present: Vec<u64>,
    populated: Vec<u64>,
    /// Per name, the name counted after it the last time: the next lookup's first guess,
    /// right wherever the entries write their members in one order.
    follows: Vec<u32>,
    last: usize,
}

impl Counts {
    pub(crate) fn new() -> Counts {
        Counts::default()
    }

    fn index(&mut self, name: &[u8]) -> usize {
        let index = self.names.index(name);
        if index == self.present.len() {
            self.present.push(0);
            self.populated.push(0);
            self.follows.push(0);
        }
        index
    }

    fn add(&mut self, name: &[u8], valued: bool) {
        let guess = self.follows.get(self.last).map_or(0, |&i| i as usize);
        let index = if self.names.names.get(guess).is_some_and(|n| **n == *name) {
            guess
        } else {
            let index = self.index(name);
            self.follows[self.last] = index as u32;
            index
        };
        self.present[index] += 1;
        self.populated[index] += u64::from(valued);
        self.last = index;
    }

    fn by_name(&self, counts: &[u64]) -> BTreeMap<Vec<u8>, u64> {
        let names = self.names.names.iter().zip(counts);
        names
            .filter(|(_, n)| **n > 0)
            .map(|(name, n)| (name.to_vec(), *n))
            .collect()
    }

    /// Per name, the entries carrying a value: the comparison's `populated_counts`. A name
    /// no entry carries with a value is left out.
    pub(crate) fn populated(&self) -> BTreeMap<Vec<u8>, u64> {
        self.by_name(&self.populated)
    }

    /// Per name, the entries carrying the member, whatever its value.
    pub(crate) fn present(&self) -> BTreeMap<Vec<u8>, u64> {
        self.by_name(&self.present)
    }
}

/// What a value's text is to the counts: the empty text, `-`, or any other.
#[derive(Clone, Copy, PartialEq)]
enum Text {
    Empty,
    Dash,
    Other,
}

/// A reader of a line `scan` has read as an object: it tells the values apart and checks
/// nothing.
struct Cursor<'a> {
    b: &'a [u8],
    p: usize,
}

impl<'a> Cursor<'a> {
    fn peek(&self) -> u8 {
        self.b.get(self.p).copied().unwrap_or(0)
    }

    fn ws(&mut self) {
        self.p = skip_ws(self.b, self.p);
    }

    /// Past the string that opens at `p`: its text in the line and whether it holds an
    /// escape.
    fn string(&mut self) -> (&'a [u8], bool) {
        let start = (self.p + 1).min(self.b.len());
        let mut i = start;
        let mut escaped = false;
        loop {
            i = stop(self.b, i);
            match self.b.get(i) {
                Some(b'\\') => {
                    escaped = true;
                    i += 2;
                }
                Some(b'"') | None => break,
                Some(_) => i += 1,
            }
        }
        let end = i.min(self.b.len());
        self.p = (end + 1).min(self.b.len());
        (&self.b[start..end], escaped)
    }

    /// Past what is left of the array or object entered `depth` deep.
    fn close(&mut self, mut depth: usize) {
        while self.p < self.b.len() {
            match self.b[self.p] {
                b'"' => {
                    self.string();
                    continue;
                }
                b'{' | b'[' => depth += 1,
                b'}' | b']' => {
                    depth -= 1;
                    if depth == 0 {
                        self.p += 1;
                        return;
                    }
                }
                _ => {}
            }
            self.p += 1;
        }
    }

    /// Past a value.
    fn skip(&mut self) {
        match self.peek() {
            b'"' => {
                self.string();
            }
            b'{' | b'[' => {
                self.p += 1;
                self.close(1);
            }
            _ => {
                while self.p < self.b.len()
                    && !matches!(
                        self.b[self.p],
                        b',' | b'}' | b']' | b' ' | b'\t' | b'\n' | b'\r'
                    )
                {
                    self.p += 1;
                }
            }
        }
    }

    /// Past a value, and what its text is: a string's own, the empty text for an empty
    /// array, the item's text for an array of one item, any other text otherwise.
    fn text(&mut self) -> Text {
        let mut arrays = 0usize;
        let mut text = loop {
            match self.peek() {
                b'[' => {
                    self.p += 1;
                    self.ws();
                    if self.peek() == b']' {
                        self.p += 1;
                        break Text::Empty;
                    }
                    arrays += 1;
                }
                b'"' => {
                    let (raw, escaped) = self.string();
                    break match raw {
                        b"" => Text::Empty,
                        b"-" => Text::Dash,
                        _ if escaped && raw.eq_ignore_ascii_case(b"\\u002d") => Text::Dash,
                        _ => Text::Other,
                    };
                }
                _ => {
                    self.skip();
                    break Text::Other;
                }
            }
        };
        for _ in 0..arrays {
            self.ws();
            if self.peek() == b']' {
                self.p += 1;
            } else {
                text = Text::Other;
                self.close(1);
            }
        }
        text
    }

    /// Past a member's name and its `:`; the name, unescaped into `unescaped` where the
    /// line escapes it. `None` where no member stands at `p`.
    fn name<'n>(&mut self, unescaped: &'n mut Vec<u8>) -> Option<&'n [u8]>
    where
        'a: 'n,
    {
        if self.peek() != b'"' {
            return None;
        }
        let quote = self.p;
        let (raw, escaped) = self.string();
        self.ws();
        if self.peek() != b':' {
            return None;
        }
        self.p += 1;
        self.ws();
        if !escaped {
            return Some(raw);
        }
        unescaped.clear();
        string(self.b, quote, unescaped, &mut Escapes::default()).ok()?;
        Some(unescaped)
    }
}

/// Adds the entries of a line that `scan` reads as an object to `counts`, and to `under[i]`
/// the consequence entries that carry a field of floor `i` with a text other than `-`.
///
/// Returns the entries of the line: 1 and the items of its consequence lists.
///
/// Two lines that are alike have the same counts. A field holding `.` is counted where a
/// member has that name (`Floors::reads_object_members`).
pub(crate) fn count(line: &[u8], counts: &mut Counts, floors: &Floors, under: &mut [u64]) -> u32 {
    let mut c = Cursor { b: line, p: 0 };
    let mut unescaped = Vec::new();
    let mut hits: Vec<u32> = Vec::new();
    let mut entries = 1u32;
    c.ws();
    if c.peek() != b'{' {
        return entries;
    }
    c.p += 1;
    loop {
        c.ws();
        if c.peek() == b',' {
            c.p += 1;
            c.ws();
        }
        let Some(name) = c.name(&mut unescaped) else {
            return entries;
        };
        let role = top_role(name);
        if role == REGULAR {
            let valued = c.text() == Text::Other;
            counts.add(name, valued);
            continue;
        }
        if role < LIST_MEMBER || c.peek() != b'[' {
            c.skip();
            continue;
        }
        c.p += 1;
        loop {
            c.ws();
            match c.peek() {
                b',' => c.p += 1,
                b'{' => {
                    c.p += 1;
                    entries = entries.saturating_add(1);
                    hits.clear();
                    loop {
                        c.ws();
                        if c.peek() == b',' {
                            c.p += 1;
                            c.ws();
                        }
                        let Some(name) = c.name(&mut unescaped) else {
                            break;
                        };
                        let text = c.text();
                        counts.add(name, text == Text::Other);
                        if text != Text::Dash {
                            if let Some(field) = floors.field(name) {
                                hits.extend_from_slice(&floors.listing[field]);
                            }
                        }
                    }
                    if c.peek() == b'}' {
                        c.p += 1;
                    }
                    hits.sort_unstable();
                    hits.dedup();
                    for &floor in &hits {
                        if let Some(n) = under.get_mut(floor as usize) {
                            *n += 1;
                        }
                    }
                }
                _ => break,
            }
        }
        if c.peek() == b']' {
            c.p += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(line: &[u8]) -> Line<'_> {
        scan(line, &Floors::default(), 0)
    }

    fn object(line: &str) -> Object<'_> {
        match read(line.as_bytes()) {
            Line::Object(object) => object,
            Line::Blank => panic!("blank: {line}"),
            Line::Malformed(shape) => panic!("{shape}: {line}"),
        }
    }

    /// The shape a line is refused for, `blank` or `object`.
    fn shape(line: &str) -> &'static str {
        match read(line.as_bytes()) {
            Line::Malformed(shape) => shape,
            Line::Blank => "blank",
            Line::Object(_) => "object",
        }
    }

    fn same(a: &str, b: &str) -> bool {
        let alike = alike(a.as_bytes(), b.as_bytes());
        assert_eq!(alike, object(a).form == object(b).form, "{a} and {b}");
        alike
    }

    fn floors(lines: &[&str]) -> Floors {
        let floor = |line: &&str| {
            let mut words = line.split(' ').map(str::to_string);
            Floor {
                column: words.next().expect("a column"),
                fields: words.collect(),
            }
        };
        Floors::new(lines.iter().map(floor).collect())
    }

    fn named(counts: BTreeMap<Vec<u8>, u64>) -> Vec<(String, u64)> {
        let name = |(name, n)| (String::from_utf8(name).expect("a name"), n);
        counts.into_iter().map(name).collect()
    }

    fn pairs(expected: &[(&str, u64)]) -> Vec<(String, u64)> {
        expected
            .iter()
            .map(|(n, c)| ((*n).to_string(), *c))
            .collect()
    }

    #[test]
    fn the_key_is_the_first_five_columns_of_the_unescaped_input() {
        let key = |line: &'static str| object(line).key;
        assert_eq!(
            *key(r#"{"input":"21\t7\trs1\tA\tC\t.\t.\t."}"#),
            *b"21\t7\trs1\tA\tC"
        );
        assert_eq!(
            *key(r#"{"input":"21\u00097\t.\tA\tC"}"#),
            *b"21\t7\t.\tA\tC"
        );
        assert_eq!(*key(r#"{"input":"21\t7"}"#), *b"21\t7");
        assert_eq!(*key(r#"{"a":1,"input":"21 7 7 A/C"}"#), *b"21 7 7 A/C");
        assert_eq!(*key(r#"{"input":"\t\t\t\t\t"}"#), *b"\t\t\t\t");
        assert_eq!(
            *key(concat!(r#"{"input":"café\\t"#, "\u{1f600}", r#""}"#)),
            *"caf\u{e9}\\t\u{1f600}".as_bytes()
        );
        for line in [
            r#"{}"#,
            r#"{"input":7}"#,
            r#"{"input":["1\t2"]}"#,
            r#"{"input":""}"#,
        ] {
            assert_eq!(*key(line), *b"", "{line}");
        }
        assert!(matches!(key(r#"{"input":"21 7"}"#), Cow::Borrowed(_)));
        assert!(matches!(key(r#"{"input":"21\t7"}"#), Cow::Owned(_)));
    }

    #[test]
    fn a_blank_line_is_one_the_comparison_strips_to_nothing() {
        for line in ["", " ", "\t", " \t ", "\x0b\x0c", "\x1c\x1d\x1e\x1f", " \r"] {
            assert_eq!(shape(line), "blank", "{line:?}");
        }
        for line in ["\0", "\x1b", "\u{a0}", "\u{3000} ", " x"] {
            assert_eq!(shape(line), NOT_ONE_OBJECT, "{line:?}");
        }
    }

    #[test]
    fn a_line_that_is_not_one_object_of_the_grammar_is_refused() {
        let refused = [
            "{",
            "[]",
            "\"x\"",
            "5",
            "null",
            "{]",
            r#"{"a"}"#,
            r#"{"a":}"#,
            r#"{"a":1,}"#,
            r#"{,"a":1}"#,
            r#"{"a":1 "b":2}"#,
            "{'a':1}",
            "{a:1}",
            r#"{"a":1}}"#,
            r#"{"a":1}{"b":2}"#,
            r#"{"a":1} x"#,
            "{\"a\":1}\x0c",
            "\x0c{\"a\":1}",
            "\u{feff}{\"a\":1}",
            "{\"a\":1}\u{a0}",
            r#"{"a":01}"#,
            r#"{"a":+1}"#,
            r#"{"a":.5}"#,
            r#"{"a":1.}"#,
            r#"{"a":1e}"#,
            r#"{"a":1e+}"#,
            r#"{"a":0x10}"#,
            r#"{"a":-}"#,
            r#"{"a":tru}"#,
            r#"{"a":True}"#,
            r#"{"a":NaN}"#,
            r#"{"a":-Infinity}"#,
            r#"{"a":"x}"#,
            "{\"a\":\"x\ty\"}",
            "{\"a\":\"x\x1fy\"}",
            "{\"a\tb\":1}",
            r#"{"a":"\x41"}"#,
            r#"{"a":"\u12"}"#,
            r#"{"a":"\u12G4"}"#,
            r#"{"a":"\U0041"}"#,
            r#"{"a":"\"}"#,
            r#"{"a":[1,]}"#,
            r#"{"a":[,1]}"#,
            r#"{"a":[1 2]}"#,
            r#"{"a":[1}"#,
            r#"{"a" 1}"#,
            r#"{"a":1,,"b":2}"#,
            r#"{"a":"\ud83d\uZZZZ"}"#,
            "{\"a\":1}\0",
        ];
        for line in refused {
            assert_eq!(shape(line), NOT_ONE_OBJECT, "{line:?}");
        }
        let read = [
            "{}",
            " {} ",
            "\t{\"a\":1}\r",
            "{ \"a\" : [ 1 , 2 ] }",
            r#"{"":{"":[]}}"#,
            r#"{"a":-0,"b":0.5e-7,"c":1E+2,"d":"\"\\\/\b\f\n\r\t\u0000\u007f"}"#,
            "{\"a\":\"\u{7f}\u{e9}\u{1f600}\"}",
        ];
        for line in read {
            assert_eq!(shape(line), "object", "{line:?}");
        }
        let deep = format!("{{\"a\":{}{}}}", "[".repeat(100_000), "]".repeat(100_000));
        assert_eq!(shape(&deep), "object");
    }

    #[test]
    fn a_member_named_twice_and_a_lone_surrogate_are_refused() {
        for line in [
            r#"{"a":1,"a":2}"#,
            r#"{"a":1,"\u0061":2}"#,
            r#"{"b":1,"a":2,"b":3}"#,
            r#"{"x":{"b":1,"c":2,"b":3}}"#,
            r#"{"input":"1","input":"2"}"#,
            r#"{"transcript_consequences":[{"a":1,"a":1}]}"#,
            r#"{"":1,"":2}"#,
        ] {
            assert_eq!(shape(line), NAMED_TWICE, "{line}");
        }
        for line in [
            r#"{"a":"\ud800"}"#,
            r#"{"a":"\udc00"}"#,
            r#"{"a":"\ud800A"}"#,
            r#"{"a":"\ud800x"}"#,
            r#"{"a":"\uDBFF\uD800"}"#,
            r#"{"\udfff":1}"#,
        ] {
            assert_eq!(shape(line), LONE_SURROGATE, "{line}");
        }
        assert_eq!(shape(r#"{"a":1,"A":2,"b":{"a":3}}"#), "object");
        assert_eq!(shape("{\"a\":\"\u{1f600}\u{1f600}\"}"), "object");
    }

    #[test]
    fn an_entry_is_the_object_and_each_item_of_a_list_that_is_an_array() {
        let line = r#"{"transcript_consequences":[{},{"a":1}],"intergenic_consequences":[],"regulatory_feature_consequences":{"a":[1]},"motif_feature_consequences":[{}]}"#;
        let read = object(line);
        assert_eq!((read.entries, read.kinds), (4, 0b1001));
        assert_eq!((object("{}").entries, object("{}").kinds), (1, 0));
        // A list under another name, and a list's name inside an item, are members.
        let line = r#"{"x":{"transcript_consequences":[{}]},"intergenic_consequences":[{"motif_feature_consequences":[{}]}]}"#;
        assert_eq!((object(line).entries, object(line).kinds), (2, 0b0010));
        for line in [
            r#"{"transcript_consequences":[1]}"#,
            r#"{"motif_feature_consequences":[{},"x"]}"#,
            r#"{"intergenic_consequences":[[{}]]}"#,
            r#"{"regulatory_feature_consequences":[null]}"#,
        ] {
            assert_eq!(shape(line), ITEM_NOT_AN_OBJECT, "{line}");
        }
    }

    #[test]
    fn the_order_of_members_is_nothing_and_the_order_of_items_is_something() {
        let a = r#"{"input":"1\t2\t.\tA\tC\t.","id":".","start":2,"colocated_variants":[{"id":"rs1","frequencies":{"C":{"af":0.1,"eur":0.2}}}],"transcript_consequences":[{"variant_allele":"C","impact":"LOW"},{"variant_allele":"C","impact":"HIGH"}]}"#;
        let b = r#"{ "transcript_consequences" : [ {"impact":"LOW","variant_allele":"C"}, {"impact":"HIGH","variant_allele":"C"} ], "colocated_variants":[{"frequencies":{"C":{"eur":0.2,"af":0.1}},"id":"rs1"}], "start":2, "id":".", "input":"1\t2\t.\tA\tC\tPASS" }"#;
        assert!(same(a, b));
        let items_swapped = r#"{"input":"1\t2\t.\tA\tC\t.","id":".","start":2,"colocated_variants":[{"id":"rs1","frequencies":{"C":{"af":0.1,"eur":0.2}}}],"transcript_consequences":[{"variant_allele":"C","impact":"HIGH"},{"variant_allele":"C","impact":"LOW"}]}"#;
        assert!(!same(a, items_swapped));
        assert!(!same(r#"{"a":[1,2]}"#, r#"{"a":[2,1]}"#));
        assert!(!same(
            r#"{"input":"1\t2\t.\tA\tC"}"#,
            r#"{"input":"1\t2\t.\tA\tG"}"#
        ));
    }

    #[test]
    fn the_members_set_aside_and_a_list_that_is_no_array_are_nothing() {
        let bare = r#"{"id":"x"}"#;
        for line in [
            r#"{"id":"x","vep_rs_version":"1.2.3"}"#,
            r#"{"vep_rs_version":{"a":[1]},"id":"x"}"#,
            r#"{"id":"x","transcript_consequences":[]}"#,
            r#"{"id":"x","intergenic_consequences":{"a":1},"motif_feature_consequences":"x"}"#,
            r#"{"id":"x","input":5}"#,
            r#"{"id":"x","input":""}"#,
        ] {
            assert!(same(bare, line), "{line}");
        }
        assert!(!same(bare, r#"{"id":"x","input":"1"}"#));
        assert!(!same(bare, r#"{"id":"x","transcript_consequences":[{}]}"#));
        // An item is an entry of its list's kind.
        assert!(!same(
            r#"{"transcript_consequences":[{"a":1}]}"#,
            r#"{"intergenic_consequences":[{"a":1}]}"#
        ));
        assert!(!same(
            r#"{"transcript_consequences":[{"a":1},{"b":2}]}"#,
            r#"{"transcript_consequences":[{"a":1}],"intergenic_consequences":[{"b":2}]}"#
        ));
        // Inside an item they are members like any other.
        assert!(!same(
            r#"{"transcript_consequences":[{"input":"1"}]}"#,
            r#"{"transcript_consequences":[{"input":"2"}]}"#
        ));
    }

    #[test]
    fn a_value_is_compared_by_its_kind_and_a_number_by_the_value_it_reads_as() {
        for (a, b) in [
            ("1", "1.0"),
            ("0", "0.0"),
            ("-0", "-0.0"),
            ("-0.0", "0.0"),
            ("100", "1e2"),
            ("0.1", "0.10000000000000002"),
            ("2.568e-05", "2.569e-05"),
            ("18014398509481985", "18014398509481984"),
            ("9223372036854775807", "9223372036854775808"),
            (
                "123456789012345678901234567890",
                "123456789012345678901234567891",
            ),
            (
                "123456789012345678901234567890",
                "123456789012345678901234567890.0",
            ),
            ("123456789012345678901234567890", "1.2345678901234568e29"),
            ("1e400", "-1e400"),
            ("1e400", "null"),
            ("1.5", "\"1.5\""),
            ("1", "\"1\""),
            ("\"a,b\"", "[\"a\",\"b\"]"),
            ("[\"a\"]", "\"a\""),
            ("[]", "\"\""),
            ("null", "\"null\""),
            ("true", "false"),
            ("{}", "[]"),
            ("{\"a\":1}", "\"{\\\"a\\\":1}\""),
            ("[[1],2]", "[1,[2]]"),
            ("[1,2]", "[12]"),
            ("\"ab\"", "\"a\",\"k\":\"b\""),
        ] {
            let (a, b) = (format!("{{\"v\":{a}}}"), format!("{{\"v\":{b}}}"));
            assert!(!same(&a, &b), "{a} and {b}");
        }
        assert!(same(r#"{"v":"a\/bé"}"#, "{\"v\":\"a/b\u{e9}\"}"));
        assert!(same(r#"{"v":1.50e+3,"w":[ ]}"#, r#"{"w":[],"v":1.50e+3}"#));
        for (a, b) in [
            ("1.0", "1.00"),
            ("1e2", "100.0"),
            ("1E2", "1e+2"),
            ("-0", "0"),
            ("2.568e-05", "0.00002568"),
            ("9.391e-06", "9.391e-6"),
            ("1.50", "15e-1"),
            ("0.1", "0.10000000000000001"),
            ("0.1", "0.1000000000000000055511151231257827"),
            ("18014398509481985.0", "18014398509481984.0"),
            (
                "123456789012345678901234567890.0",
                "123456789012345678901234567891.0",
            ),
            ("123456789012345678901234567890.0", "1.2345678901234568e29"),
            ("0e5", "0.0"),
            ("-0e0", "-0.0"),
            ("1e-999", "0.0"),
            ("1e400", "1e999"),
            ("-1e400", "-17e308"),
        ] {
            let (a, b) = (format!("{{\"v\":{a}}}"), format!("{{\"v\":{b}}}"));
            assert!(same(&a, &b), "{a} and {b}");
        }
    }

    /// A line of the shape of an output with co-located variants: `number` is a member of
    /// the record, of an object and of an array inside an array of the record, and of a
    /// consequence item and an object of it.
    fn numbered(number: &str) -> String {
        format!(
            r#"{{"input":"21\t7\trs1\tA\tC","start":7,"score":{number},"colocated_variants":[{{"id":"rs1","frequencies":{{"C":{{"af":{number},"eur":[{number}]}}}}}}],"transcript_consequences":[{{"variant_allele":"C","cadd":{number},"plug":{{"x":{number}}}}}]}}"#
        )
    }

    #[test]
    fn two_spellings_of_a_number_are_alike_exactly_where_the_comparison_reads_one_value() {
        let entries = |line: &str| {
            crate::parse::json::parse_json(line, &[]).expect("the comparison reads the line")
        };
        let numbers = [
            "0",
            "-0",
            "1",
            "-1",
            "10",
            "100",
            "1.0",
            "1.00",
            "-1.0",
            "1e0",
            "1E0",
            "1e+0",
            "10e-1",
            "0.1e1",
            "1e2",
            "1.0e2",
            "100.0",
            "0.0",
            "0.00",
            "0e0",
            "0e9",
            "-0.0",
            "-0e0",
            "0.00002568",
            "2.568e-05",
            "2.568e-5",
            "2.568E-5",
            "25.68e-6",
            "2.569e-05",
            "9.391e-06",
            "9.391e-6",
            "0.000009391",
            "0.1",
            "0.10",
            "0.10000000000000001",
            "0.1000000000000000055511151231257827",
            "0.10000000000000002",
            "1.5",
            "1.50",
            "15e-1",
            "18014398509481984",
            "18014398509481985",
            "18014398509481984.0",
            "18014398509481985.0",
            "9223372036854775807",
            "9223372036854775808",
            "9223372036854775808.0",
            "-9223372036854775808",
            "-9223372036854775809",
            "123456789012345678901234567890",
            "123456789012345678901234567891",
            "123456789012345678901234567890.0",
            "123456789012345678901234567891.0",
            "1.2345678901234568e29",
            "1e22",
            "10000000000000000000000.0",
            "1e23",
            "1e-7",
            "0.0000001",
            "5e-324",
            "4.9e-324",
            "1e-999",
            "-1e-999",
            "1.7976931348623157e308",
            "17976931348623157e292",
        ];
        let mut alike_pairs = 0;
        for a in numbers {
            for b in numbers {
                let (x, y) = (numbered(a), numbered(b));
                let alike = same(&x, &y);
                assert_eq!(alike, entries(&x) == entries(&y), "{a} and {b}");
                alike_pairs += usize::from(alike && a != b);
            }
        }
        assert_eq!(alike_pairs, 136);
        // A number past the range of a float is no number to the comparison, which reads
        // it as it reads `null`: one such number is alike another of its sign alone.
        for (a, b, alike) in [
            ("1e400", "1e999", true),
            ("-1e400", "-17e308", true),
            ("1e400", "-1e400", false),
            ("1e400", "null", false),
        ] {
            let (x, y) = (numbered(a), numbered(b));
            assert_eq!(entries(&x), entries(&y), "{a} and {b}");
            assert_eq!(same(&x, &y), alike, "{a} and {b}");
        }
    }

    #[test]
    fn the_shapes_that_hand_a_group_or_everything_over_are_seen() {
        let flags = |line: &str| {
            let read = object(line);
            (
                read.opens_object,
                read.dotted_member,
                read.custom_annotations,
                read.key_mark,
            )
        };
        let item = |member: &str| format!("{{\"transcript_consequences\":[{{{member}}}]}}");
        for member in [
            r#""x":{}"#,
            r#""x":"{""#,
            r#""x":"{x""#,
            r#""x":"\u007bx""#,
            r#""x":[{"a":1},2]"#,
            r#""x":[[["{"]],1]"#,
        ] {
            assert_eq!(
                flags(&item(member)),
                (true, false, false, false),
                "{member}"
            );
        }
        for member in [
            r#""x":" {""#,
            r#""x":[[],{}]"#,
            r#""x":["","{"]"#,
            r#""x":[[[]],"{"]"#,
            r#""x":[1,{}]"#,
            r#""x":"""#,
            r#""variant_allele":"-","transcript_id":"T1.2","gene_id":"a|b#c""#,
        ] {
            assert_eq!(
                flags(&item(member)),
                (false, false, false, false),
                "{member}"
            );
        }
        assert_eq!(flags(&item(r#""a.b":1"#)), (false, true, false, false));
        assert_eq!(flags(&item(r#""a\u002eb":1"#)), (false, true, false, false));
        assert_eq!(
            flags(&item(r#""x":{"a.b":1}"#)),
            (true, false, false, false)
        );
        for member in [
            r#""variant_allele":"A|C""#,
            r#""variant_allele":"|""#,
            r#""transcript_id":"T#2""#,
            r##""regulatory_feature_id":["a","#"]"##,
            r#""motif_feature_id":"x|""#,
            r#""variant_allele":"\u007C""#,
            r#""transcript_id":"\u0023""#,
        ] {
            assert_eq!(
                flags(&item(member)),
                (false, false, false, true),
                "{member}"
            );
        }
        assert_eq!(
            flags(&item(r#""variant_allele":{"a|b":1}"#)),
            (true, false, false, true)
        );
        // The record's own members hand nothing over but the custom annotations.
        assert_eq!(
            flags(r#"{"x":{"a":1},"a.b":"{","variant_allele":"|","custom_annotations":null}"#),
            (false, false, true, false)
        );
        assert_eq!(
            flags(&item(r#""custom_annotations":1"#)),
            (false, false, false, false)
        );
    }

    #[test]
    fn a_unicode_escape_is_seen_wherever_it_stands() {
        for line in [
            r#"{"shift_length":"\u00b2"}"#,
            r#"{"input":"1\u00092"}"#,
            r#"{"vep_rs_version":["\u0031"]}"#,
            r#"{"a\u0062":1}"#,
            r#"{"transcript_consequences":[{"x":{"y":"\\\u005C"}}]}"#,
            r#"{"a":"\ud83d\ude00"}"#,
        ] {
            assert!(object(line).unicode_escape, "{line}");
        }
        for line in [
            r#"{"shift_length":"2"}"#,
            r#"{"input":"1\t2","a":"\"\\\/\b\f\n\r"}"#,
            r#"{"a":"\\u00b2","b":"\\\\u0041"}"#,
            "{\"shift_length\":\"\u{b2}\"}",
        ] {
            assert!(!object(line).unicode_escape, "{line}");
        }
    }

    #[test]
    fn a_line_s_depth_is_its_objects_and_arrays_open_at_once() {
        for (line, depth) in [
            (r#"{}"#, 1),
            (r#"{"a":1,"b":"[[[{"}"#, 1),
            (r#"{"a":[]}"#, 2),
            (r#"{"a":{},"b":[1,2]}"#, 2),
            (r#"{"a":[[],[[]],{"b":1}]}"#, 4),
            (
                r#"{"transcript_consequences":[{"x":1},{"x":[1]},{"x":2}]}"#,
                4,
            ),
            (r#"{"a":[{"b":[{"c":{}}]}],"d":[[1]]}"#, 6),
        ] {
            assert_eq!(object(line).depth, depth, "{line}");
        }
        // Forty items of a list stand side by side, and 70 arrays inside one another.
        let items = vec![r#"{"x":["y"]}"#; 40].join(",");
        let flat = format!(r#"{{"transcript_consequences":[{items}]}}"#);
        assert_eq!(object(&flat).depth, 4);
        let nested = format!(r#"{{"a":{}{}}}"#, "[".repeat(70), "]".repeat(70));
        assert_eq!(object(&nested).depth, 71);
    }

    #[test]
    fn a_floor_is_carried_where_an_entry_has_one_of_its_fields_as_a_member() {
        let floors = floors(&[
            "IMPACT impact",
            "Feature transcript_id motif_feature_id",
            "id id",
            "x never",
        ]);
        let carried = |line: &str, uncounted: u8| match scan(line.as_bytes(), &floors, uncounted) {
            Line::Object(read) => read.carried,
            _ => panic!("not read: {line}"),
        };
        let line = r#"{"id":"x","transcript_consequences":[{"impact":"-"}],"motif_feature_consequences":[{"motif_feature_id":"m","impact":1}],"colocated_variants":[{"never":1}]}"#;
        assert_eq!(carried(line, 0), 0b111);
        assert_eq!(carried(line, 0b0001), 0b111);
        assert_eq!(carried(line, 0b1000), 0b101);
        assert_eq!(carried(line, 0b1001), 0b100);
        assert_eq!(carried(r#"{"impact":1,"input":"id"}"#, 0b1111), 0b001);
        assert_eq!(
            carried(
                r#"{"x":{"impact":1},"transcript_consequences":{"impact":1}}"#,
                0
            ),
            0
        );
        assert_eq!(object(line).carried, 0);
        // The 65th floor and those after it are never reported.
        let many: Vec<String> = (0..70).map(|i| format!("c{i} f{i}")).collect();
        let many = self::floors(&many.iter().map(String::as_str).collect::<Vec<_>>());
        let Line::Object(read) = scan(br#"{"f0":1,"f63":1,"f64":1,"f69":1}"#, &many, 0) else {
            panic!("not read");
        };
        assert_eq!(read.carried, 1 | 1 << 63);
    }

    #[test]
    fn the_counts_are_the_comparison_s_of_the_line_s_entries() {
        let floors = floors(&[
            "IMPACT impact",
            "both impact biotype",
            "Feature transcript_id",
            "o a.b",
        ]);
        assert!(floors.reads_object_members());
        assert!(!self::floors(&["IMPACT impact"]).reads_object_members());
        let mut counts = Counts::new();
        let mut under = [0u64; 4];
        let line = r#"{"input":"1\t2","vep_rs_version":"1","id":"rs1","empty":"","dash":"-","dashes":["-"],"none":[[]],"two":["-","-"],"zero":0,"transcript_consequences":[{"impact":"HIGH","biotype":"-","id":[""]},{"impact":"-","biotype":["-"],"a.b":""},{"impact":[],"transcript_id":"T"}],"intergenic_consequences":"x","motif_feature_consequences":[{"none":null}]}"#;
        count(line.as_bytes(), &mut counts, &floors, &mut under);
        assert_eq!(
            named(counts.populated()),
            pairs(&[
                ("id", 1),
                ("impact", 1),
                ("none", 1),
                ("transcript_id", 1),
                ("two", 1),
                ("zero", 1)
            ])
        );
        assert_eq!(
            named(counts.present()),
            pairs(&[
                ("a.b", 1),
                ("biotype", 2),
                ("dash", 1),
                ("dashes", 1),
                ("empty", 1),
                ("id", 2),
                ("impact", 3),
                ("none", 2),
                ("transcript_id", 1),
                ("two", 1),
                ("zero", 1)
            ])
        );
        // An entry counts once under a floor, with any text but `-`, the empty one too.
        assert_eq!(under, [2, 2, 1, 1]);
        // A second line adds to the counts of the first.
        let mut again = vec![0u64; under.len()];
        count(line.as_bytes(), &mut counts, &floors, &mut again);
        count(br#"{"new":1}"#, &mut counts, &floors, &mut again);
        assert_eq!(counts.populated()[&b"impact"[..]], 2);
        assert_eq!(counts.present()[&b"biotype"[..]], 4);
        assert_eq!(counts.present()[&b"new"[..]], 1);
        assert_eq!(again, under);
    }

    #[test]
    fn many_names_keep_their_counts() {
        let mut counts = Counts::new();
        for round in 0..3 {
            for i in (0..500).rev() {
                let line = format!("{{\"n{i}\":{round},\"m{}\":\"-\"}}", i % 7);
                count(line.as_bytes(), &mut counts, &Floors::default(), &mut []);
            }
        }
        let (populated, present) = (counts.populated(), counts.present());
        assert_eq!((populated.len(), present.len()), (500, 507));
        assert!(populated.values().all(|&n| n == 3));
        assert_eq!(present[&b"m0"[..]], 3 * 72);
    }

    /// What a cases file held.
    #[derive(Default)]
    struct Tally {
        lines: u64,
        objects: u64,
        /// Per shape, the lines the comparison's interpreter reads and this module refuses.
        refused_more: BTreeMap<String, u64>,
        /// Per relation of a pair: the pairs, those of equal entries, those alike.
        pairs: BTreeMap<String, [u64; 3]>,
        /// The floors with a field holding `.` on lines with an object's text: those `count`
        /// answers as the comparison does, and the others.
        object_member_floors: [u64; 2],
    }

    fn unquoted(text: &str) -> Vec<u8> {
        let mut bytes = text.bytes();
        let mut out = Vec::new();
        while let Some(b) = bytes.next() {
            if b == b'%' {
                let hex = [
                    bytes.next().expect("a digit"),
                    bytes.next().expect("a digit"),
                ];
                let hex = std::str::from_utf8(&hex).expect("two digits");
                out.push(u8::from_str_radix(hex, 16).expect("two hex digits"));
            } else {
                out.push(b);
            }
        }
        out
    }

    fn table(text: &str) -> BTreeMap<Vec<u8>, u64> {
        let pair = |pair: &str| {
            let (name, n) = pair.split_once(':').expect("a count");
            (unquoted(name), n.parse().expect("a number"))
        };
        text.split(',')
            .filter(|p| !p.is_empty())
            .map(pair)
            .collect()
    }

    fn numbers(text: &str) -> Vec<u64> {
        text.split(',')
            .map(|n| n.parse().expect("a number"))
            .collect()
    }

    fn shown(line: &[u8]) -> String {
        String::from_utf8_lossy(&line[..line.len().min(400)]).into_owned()
    }

    /// Holds one line to the comparison's answers: `kind`, the status and, for an object,
    /// its `name=value` answers.
    fn line_case(words: &[&str], line: &[u8], floors: &Floors, tally: &mut Tally) {
        let (kind, status) = (words[0], words[1]);
        let answer: BTreeMap<&str, &str> = words[2..]
            .iter()
            .map(|w| w.split_once('=').expect("an answer"))
            .collect();
        let more: u8 = answer.get("more").map_or(0, |n| n.parse().expect("bits"));
        tally.lines += 1;
        let mut refused =
            |shape: &str| *tally.refused_more.entry(shape.to_string()).or_default() += 1;
        match (status, scan(line, floors, 0)) {
            ("blank", Line::Blank) | ("error", Line::Malformed(_)) => {}
            ("blank", Line::Malformed(_)) if !line.is_ascii() => {
                refused("white space outside ASCII alone");
            }
            ("object", Line::Malformed(shape)) => {
                let bit = match shape {
                    NOT_ONE_OBJECT => 1,
                    NAMED_TWICE => 2,
                    LONE_SURROGATE => 4,
                    _ => 0,
                };
                assert!(more & bit != 0, "{shape}, more={more}: {}", shown(line));
                refused(if bit == 1 {
                    "a NaN or Infinity constant"
                } else {
                    shape
                });
            }
            ("object", Line::Object(object)) => {
                assert_eq!(more, 0, "read: {}", shown(line));
                tally.objects += 1;
                object_case(&answer, line, &object, floors, tally);
            }
            (status, _) => panic!("the comparison reads {status} ({kind}): {}", shown(line)),
        }
    }

    fn object_case(
        answer: &BTreeMap<&str, &str>,
        line: &[u8],
        object: &Object,
        floors: &Floors,
        tally: &mut Tally,
    ) {
        let at = shown(line);
        assert_eq!(*object.key, *unquoted(answer["key"]), "key: {at}");
        assert_eq!(
            object.entries.to_string(),
            answer["entries"],
            "entries: {at}"
        );
        assert_eq!(object.kinds.to_string(), answer["kinds"], "kinds: {at}");
        let flags = [
            object.opens_object,
            object.dotted_member,
            object.custom_annotations,
            object.key_mark,
        ];
        let flags: String = flags.iter().map(|&f| if f { '1' } else { '0' }).collect();
        assert_eq!(flags, answer["flags"], "flags: {at}");
        assert_eq!(
            object.carried.to_string(),
            answer["carried"],
            "carried: {at}"
        );
        let escape = if object.unicode_escape { "1" } else { "0" };
        assert_eq!(escape, answer["uescape"], "a unicode escape: {at}");
        let by = numbers(answer["by"]);
        for uncounted in [0b1111, (line.len() % 16) as u8] {
            let Line::Object(part) = scan(line, floors, uncounted) else {
                panic!("not read again: {at}");
            };
            let counted = (0..4).filter(|i| uncounted & (1 << i) == 0);
            let carried = counted.fold(by[0], |all, i| all | by[i + 1]);
            assert_eq!(part.carried, carried, "carried without {uncounted}: {at}");
            assert_eq!(part.form, object.form, "form: {at}");
        }
        let mut counts = Counts::new();
        let mut under = vec![0u64; floors.floors().len()];
        assert_eq!(
            count(line, &mut counts, floors, &mut under),
            object.entries,
            "entries counted: {at}"
        );
        assert_eq!(
            counts.populated(),
            table(answer["populated"]),
            "populated: {at}"
        );
        assert_eq!(counts.present(), table(answer["present"]), "present: {at}");
        let expected = numbers(answer["under"]);
        for (i, floor) in floors.floors().iter().enumerate() {
            if object.opens_object && floor.fields.iter().any(|f| f.contains('.')) {
                tally.object_member_floors[usize::from(under[i] != expected[i])] += 1;
            } else {
                assert_eq!(under[i], expected[i], "under {}: {at}", floor.column);
            }
        }
    }

    /// Holds a pair of lines to whether the comparison parses them into equal entries.
    fn pair_case(relation: &str, equal: bool, a: &[u8], b: &[u8], tally: &mut Tally) {
        let (Line::Object(x), Line::Object(y)) = (read(a), read(b)) else {
            assert!(!equal && relation == "unread", "not read: {}", shown(a));
            return;
        };
        let alike_ab = alike(a, b);
        let at = || format!("{}\n{}", shown(a), shown(b));
        assert_eq!(alike_ab, alike(b, a), "one way alone: {}", at());
        assert_eq!(alike_ab, x.form == y.form, "digests: {}", at());
        assert!(equal || !alike_ab, "alike, of other entries: {}", at());
        let one_form = ["shuffle", "respell_escape", "list:empty_absent_or_no_array"];
        if one_form.contains(&relation) || relation.starts_with("set_aside:") {
            assert!(equal && alike_ab, "{relation}, not alike: {}", at());
        }
        let slot = tally.pairs.entry(relation.to_string()).or_default();
        slot[0] += 1;
        slot[1] += u64::from(equal);
        slot[2] += u64::from(alike_ab);
    }

    /// Holds every case of a cases text: `F column field...` lines name the floors; a line
    /// case is `L kind answers` and the line; a pair is `P relation 1|0`, whether the
    /// comparison parses the two into equal entries, and the two lines.
    fn check(cases: &[u8]) -> Tally {
        let mut lines = cases.split(|&b| b == b'\n');
        let mut listed = Vec::new();
        let mut floors = None;
        let mut tally = Tally::default();
        while let Some(head) = lines.next() {
            let head = std::str::from_utf8(head).expect("a case's first line");
            let words: Vec<&str> = head.split(' ').collect();
            match words[0] {
                "" => {}
                "F" => listed.push(Floor {
                    column: words[1].to_string(),
                    fields: words[2..].iter().map(|f| (*f).to_string()).collect(),
                }),
                "L" => {
                    let floors =
                        floors.get_or_insert_with(|| Floors::new(std::mem::take(&mut listed)));
                    let line = lines.next().expect("a case's line");
                    line_case(&words[1..], line, floors, &mut tally);
                }
                "P" => {
                    let (a, b) = (lines.next().expect("a line"), lines.next().expect("a line"));
                    pair_case(words[1], words[2] == "1", a, b, &mut tally);
                }
                _ => panic!("no case: {head}"),
            }
        }
        tally
    }

    /// Lines and pairs of every kind with the comparison's own answers, `parse_json` and its
    /// counts run over each line alone: 106 of the 207 lines are objects to both, 21 are
    /// read by the comparison's interpreter and refused here; of the 64 pairs both read, 36
    /// have equal entries and 22 of those are alike, the two that spell one float two ways
    /// among them.
    #[test]
    fn the_committed_cases_hold() {
        let tally = check(include_bytes!("json_cases.txt"));
        let sum = |all: [u64; 3], p: &[u64; 3]| [all[0] + p[0], all[1] + p[1], all[2] + p[2]];
        assert_eq!((tally.lines, tally.objects), (207, 106));
        assert_eq!(tally.refused_more.values().sum::<u64>(), 21);
        assert_eq!(tally.pairs.values().fold([0; 3], sum), [64, 36, 22]);
        assert_eq!(tally.pairs["number"], [4, 2, 2]);
    }
}

#[cfg(test)]
mod format_tests {
    use super::*;
    use crate::stream::Format;
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicU32, Ordering};

    static SCRATCH: AtomicU32 = AtomicU32::new(0);

    fn scratch(text: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-compare-json-{}-{}",
            std::process::id(),
            SCRATCH.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, text).unwrap();
        path
    }

    fn plan(input: &Path) -> Plan {
        Plan {
            format: Format::Json,
            force_locations: Vec::new(),
            force_text: Vec::new(),
            keep_values: None,
            unordered: false,
            threads: 1,
            input: Some(input.to_path_buf()),
            force_records: Vec::new(),
            floor_fields: Vec::new(),
        }
    }

    #[test]
    fn a_json_line_is_a_blank_line_or_a_line_of_a_group() {
        for line in ["", " ", "\t \x0b\x0c", "\x1c\x1d\x1e\x1f"] {
            assert!(Json::is_blank(line.as_bytes()), "{line:?}");
        }
        for line in ["{}", "#x", " x", "\u{a0}", "\0"] {
            assert!(!Json::is_blank(line.as_bytes()), "{line:?}");
            assert!(!Json::is_header(line.as_bytes()), "{line:?}");
        }
    }

    #[test]
    fn a_line_s_scan_reads_its_key_its_form_and_what_hands_it_over() {
        let input = scratch("1\t7\trs1\tA\tG\n1\t8\trs2\tA\tG\n");
        let mut plan = plan(&input);
        plan.force_records = vec![0];
        plan.force_text = vec!["misspelt".to_string()];
        plan.keep_values = Some(("kind".to_string(), vec!["transcript".to_string()]));
        plan.floor_fields = vec![
            ("IMPACT".to_string(), vec!["impact".to_string()]),
            ("am_class".to_string(), vec!["am_class".to_string()]),
        ];
        let format = Json::new(&plan, 64).unwrap();
        std::fs::remove_file(input).unwrap();
        let reference = format.side((), true, true);
        let candidate = format.side((), true, false);
        let scan = |side: &JsonSide, line: &str| side.scan(&mut (), line.as_bytes());
        let listed = scan(&reference, r#"{"input":"1\t7\trs1\tA\tG\t.","id":"rs1"}"#);
        assert!(listed.forced && !listed.fragile);
        assert_eq!((listed.key, listed.shapes), (0, 0));
        let plain = r#"{"input":"1\t8\trs2\tA\tG","transcript_consequences":[{"impact":"LOW"}]}"#;
        let turned =
            r#"{ "transcript_consequences": [{"impact": "LOW"}], "input": "1\t8\trs2\tA\tG\tGT" }"#;
        let kept = scan(&reference, plain);
        assert!(!kept.forced && !kept.strict);
        assert_eq!((kept.key, kept.shapes, kept.carried), (1, 0, 1));
        assert_eq!(kept.hash, scan(&candidate, turned).hash);
        let kept_text = |side: &JsonSide, line: &str| {
            let mut text = Vec::new();
            side.kept_text(line.as_bytes(), &mut text);
            text
        };
        assert_eq!(kept_text(&reference, plain), kept_text(&candidate, turned));
        let other = r#"{"input":"1\t8\trs2\tA\tG","transcript_consequences":[{"impact":"HIGH"}]}"#;
        assert_ne!(kept_text(&reference, plain), kept_text(&reference, other));
        // A frequency spelt two ways is one form and one text; a whole number and the
        // float of its value are two.
        let frequency = |number: &str| {
            format!(
                r#"{{"input":"1\t8\trs2\tA\tG","colocated_variants":[{{"frequencies":{{"G":{{"af":{number}}}}}}}]}}"#
            )
        };
        let (padded, positional) = (frequency("2.568e-05"), frequency("0.00002568"));
        assert_eq!(
            scan(&reference, &padded).hash,
            scan(&candidate, &positional).hash
        );
        assert_eq!(
            kept_text(&reference, &padded),
            kept_text(&candidate, &positional)
        );
        assert_ne!(
            kept_text(&reference, &frequency("1")),
            kept_text(&candidate, &frequency("1.0"))
        );
        // A line end and the byte that stands before a stand-in are written as two bytes:
        // the text is one line, and two forms that differ in such a byte are two texts.
        let ended = kept_text(&reference, r#"{"a":"x\ny"}"#);
        assert_eq!(ended.iter().filter(|byte| **byte == b'\n').count(), 1);
        assert_eq!(ended.last(), Some(&b'\n'));
        assert_ne!(ended, kept_text(&reference, r#"{"a":"x\u0001\u0002y"}"#));
        assert_ne!(ended, kept_text(&reference, r#"{"a":"x\u0001\u0003y"}"#));
        // A line that is no object is its own text.
        assert_eq!(kept_text(&reference, "[1]"), kept_text(&candidate, "[1]"));
        assert_ne!(kept_text(&reference, "[1]"), kept_text(&candidate, "[2]"));
        assert_eq!(reference.entries(plain.as_bytes()), 2);
        // A list of a kind off the kept list, and a listed text, name reference lines
        // alone; an item of such a list carries no field of a floor there.
        let regulatory =
            r#"{"input":"1\t8\trs2\tA\tG","regulatory_feature_consequences":[{"am_class":"x"}]}"#;
        let (dropped, written) = (scan(&reference, regulatory), scan(&candidate, regulatory));
        assert!(dropped.forced && !written.forced);
        assert_eq!((dropped.carried, written.carried), (0, 2));
        let named = r#"{"input":"1\t8\trs2\tA\tG","note":"misspelt7"}"#;
        assert!(scan(&reference, named).forced && !scan(&candidate, named).forced);
        // What a floor without a field may be looked for in.
        for line in [
            r#"{"input":"x","transcript_consequences":[{"plug":{"am_class":"x"}}]}"#,
            r#"{"input":"x","transcript_consequences":[{"plug.x":1}]}"#,
            r#"{"input":"x","custom_annotations":{}}"#,
        ] {
            assert!(scan(&candidate, line).strict, "{line}");
        }
        assert!(format.strict(0) && format.strict(1) && !format.strict(3));
        let shape = |line: &str| {
            let shapes = scan(&reference, line).shapes;
            Json::SHAPES[shapes.trailing_zeros() as usize]
        };
        assert_eq!(shape("[]"), NOT_ONE_OBJECT);
        assert_eq!(shape(r#"{"a":1,"a":2}"#), NAMED_TWICE);
        assert_eq!(shape(r#"{"a":"\ud800"}"#), LONE_SURROGATE);
        assert_eq!(
            shape(r#"{"transcript_consequences":[1]}"#),
            ITEM_NOT_AN_OBJECT
        );
        assert_eq!(
            shape(r#"{"transcript_consequences":[{"variant_allele":"A|C"}]}"#),
            MEMBER_MARK
        );
        assert_eq!(
            shape(r#"{"input":"1\t8\tmisspelt2\tA\tG"}"#),
            LISTED_TEXT_SHAPE
        );
        assert_eq!(reference.entries(b"[]"), 1);
        // A line past a reader's limit: a byte outside ASCII, a `\u` escape, a nesting.
        for line in [
            "{\"a\":\"\u{e9}\"}".to_string(),
            r#"{"a":"\u0041"}"#.to_string(),
            format!("{{\"a\":{}1{}}}", "[".repeat(64), "]".repeat(64)),
        ] {
            let scanned = scan(&candidate, &line);
            assert!(scanned.forced && scanned.fragile, "{line}");
        }
        let nested = format!("{{\"a\":{}1{}}}", "[".repeat(63), "]".repeat(63));
        assert!(!scan(&candidate, &nested).fragile);
    }

    #[test]
    fn a_plan_names_what_the_format_reads() {
        let input = scratch("1\t7\trs1\tA\tG\n");
        let refused = |plan: &Plan| Json::new(plan, 64).err().unwrap().to_string();
        let mut kept = plan(&input);
        kept.keep_values = Some(("Feature_type".to_string(), Vec::new()));
        assert_eq!(
            refused(&kept),
            "the column of the kept values of the json format is kind"
        );
        let mut lacking = plan(&input);
        lacking.input = None;
        assert_eq!(
            refused(&lacking),
            "the json format needs the input the outputs were written from"
        );
        let without_floors = Json::new(&plan(&input), 64).unwrap();
        assert!(!without_floors.strict(0));
        std::fs::remove_file(input).unwrap();
    }
}
