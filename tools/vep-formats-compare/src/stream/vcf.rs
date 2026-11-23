// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The VCF format as `parse_vcf` (`crate::parse::vcf`) reads it: a data line is one
//! record entry and one entry a chunk of its CSQ, the last `CSQ=` item of its INFO column,
//! each chunk cut by the Format list of the last `##INFO=<ID=CSQ` header line. A line names
//! its input record by its first five columns.
//!
//! A record group is every data line of a side that shares those five columns. The lines
//! of one key name the input's lines of that key in order, and a line past them is
//! labelled by its place in the output: a group of more lines than the input holds of its
//! key is handed over, and its lines past that count are lines of one more group, which
//! the order of the shared keys reads beside the groups by key.

use super::columns::hash;
use super::digest;
use super::grouping::{Grouping, Scanned, Side, Tally};
use super::lines::lines;
use super::source::{Block, Buffers, Source};
use super::work::{in_order, joined, Closing, Collecting, Jobs, Outcome, Permits};
use super::{utf8, Format, Plan};
use crate::error::{Error, Result};
use memchr::{memchr, memchr2_iter, memchr_iter, memmem};
use std::any::Any;
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::iter::once;
use std::ops::Range;
use std::panic::resume_unwind;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::channel;
use std::sync::Arc;

const CSQ_HEADER: &[u8] = b"##INFO=<ID=CSQ";
const FORMAT: &[u8] = b"Format: ";
const CSQ_ITEM: &[u8] = b"CSQ=";
/// The CSQ fields the comparison joins into an entry's key.
const KEY_FIELDS: [&[u8]; 3] = [b"Allele", b"Feature", b"Feature_type"];
/// The fields of a line's record entry: its first seven columns, its INFO items other
/// than the CSQ joined by `;`, and whether it has a CSQ.
const RECORD_FIELDS: [&[u8]; 9] = [
    b"CHROM",
    b"POS",
    b"ID",
    b"REF",
    b"ALT",
    b"QUAL",
    b"FILTER",
    b"INFO_without_CSQ",
    b"CSQ_present",
];

const HEADER_AFTER_DATA: u32 = 0;
const NO_FORMAT_LIST: u32 = 1;
const REPEATED_FIELD: u32 = 2;
const FORMAT_LISTS_DIFFER: u32 = 3;
const SHORT_LINE: u32 = 4;
const CHUNK_WIDTH: u32 = 5;
const KEY_MARK: u32 = 6;
const LISTED_TEXT: u32 = 7;

/// The shape of a listed text in a header line of the reference or in a line's key, in
/// the words the VCF and JSON formats share.
pub(crate) const LISTED_TEXT_SHAPE: &str =
    "a listed text stands in a header line of the reference or in a line's key";

/// Set on the key of a line whose key the input lacks: the rest is the hash of its key.
const UNKNOWN: u64 = 1 << 63;
/// The key of the group of the lines the comparison labels by their place in the output.
/// The key of an input's line is below 2^32.
pub(crate) const PLACED: u64 = 1 << 62;

/// A data line's identity, as `vcf_line_key` reads it: its first five columns as text.
pub(crate) fn line_key(text: &[u8]) -> &[u8] {
    memchr_iter(b'\t', text)
        .nth(4)
        .map_or(text, |end| &text[..end])
}

/// The longest run of digits every build of the comparison's interpreter converts to an
/// integer: a longer one may be refused.
const DIGIT_RUN: usize = 640;
/// The opening brackets under which an object's nesting stays well inside the
/// interpreter's call depth: the brackets of a line of text, the levels of a JSON line.
pub(crate) const OPENING_BRACKETS: usize = 64;
/// The digits of a POS the comparison reads as a whole number of 64 bits.
const POSITION_DIGITS: usize = 18;

/// Whether the line holds a byte outside ASCII or a run of digits longer than every build
/// of the comparison's interpreter converts.
pub(crate) fn unconverted(text: &[u8]) -> bool {
    !text.is_ascii()
        || (text.len() > DIGIT_RUN
            && text
                .split(|b| !b.is_ascii_digit())
                .any(|run| run.len() > DIGIT_RUN))
}

/// Whether the comparison may fail to read a value of the line: it converts some values
/// (a whole number, an object's text) before it compares them, takes a digit outside
/// ASCII for a digit without converting it, refuses a number past a length, and stops at
/// an object nested past its call depth. A line of
/// such a shape hands its group over, so the comparison meets it whatever is kept back.
fn past_a_reader_s_limit(text: &[u8]) -> bool {
    unconverted(text)
        || memchr2_iter(b'[', b'{', text)
            .nth(OPENING_BRACKETS)
            .is_some()
}

/// What a worker reads off one block of the input: the keys of its data lines, one key
/// after another, and where each key ends; or the message of text that is not UTF-8.
#[derive(Default)]
struct BlockKeys {
    bytes: Vec<u8>,
    ends: Vec<usize>,
    failure: Option<String>,
}

impl BlockKeys {
    fn of(block: &Block, skip: &AtomicBool) -> BlockKeys {
        let bytes = block.bytes();
        let mut keys = BlockKeys {
            failure: utf8::failure(bytes, block.offset),
            ..BlockKeys::default()
        };
        if keys.failure.is_some() || skip.load(Ordering::Relaxed) {
            return keys;
        }
        for line in lines(bytes) {
            let line = &bytes[line];
            if line.first() != Some(&b'#') {
                keys.bytes.extend_from_slice(line_key(line));
                keys.ends.push(keys.bytes.len());
            }
        }
        keys
    }
}

/// The input's data lines by key: an output line names its input record by its key, the
/// lines of one key taking the input's lines of that key in order. Each key of the input
/// is a group, the groups numbered as the input first names their keys.
pub(crate) struct InputKeys {
    hash: fn(&[u8]) -> u64,
    /// The keys' bytes, one key after another, and where each key ends.
    bytes: Vec<u8>,
    ends: Vec<u64>,
    /// Key hash to the group of the input's first key of that hash.
    ids: HashMap<u64, u32>,
    /// The groups of the keys whose hash an earlier key of the input has.
    others: HashMap<Vec<u8>, u32>,
    /// Per group, the input's data lines of its key.
    lines: Vec<u32>,
    /// Per group, whether a listed record is a line of its key.
    forced: Vec<bool>,
}

impl InputKeys {
    fn new(hash: fn(&[u8]) -> u64) -> InputKeys {
        InputKeys {
            hash,
            bytes: Vec::new(),
            ends: Vec::new(),
            ids: HashMap::new(),
            others: HashMap::new(),
            lines: Vec::new(),
            forced: Vec::new(),
        }
    }

    /// The keys of the data lines of the input at `path`, every line that does not start
    /// with `#`, read once through blocks of `block` bytes: one thread reads (and
    /// inflates), `threads` workers cut the blocks into keys. `records` are ordinals of
    /// those lines, counted from 0: the group of a listed line's key is handed over.
    pub(crate) fn read(
        path: &Path,
        records: &[u64],
        block: usize,
        threads: usize,
    ) -> Result<InputKeys> {
        // Two keys of one hash that the input does not tell apart form one group, which is
        // kept back only when the lines of both are the same on the two sides.
        InputKeys::hashed(path, records, block, threads, hash)
    }

    fn hashed(
        path: &Path,
        records: &[u64],
        block: usize,
        threads: usize,
        hash: fn(&[u8]) -> u64,
    ) -> Result<InputKeys> {
        let mut records = records.to_vec();
        records.sort_unstable();
        records.dedup();
        let mut listed = records.into_iter().peekable();
        let mut keys = InputKeys::new(hash);
        let source = Source::open(path)?;
        let buffers = Buffers::new(block.max(1));
        let threads = threads.max(1);
        let permits = Permits::new(threads + 2);
        let (stop, skip) = (AtomicBool::new(false), AtomicBool::new(false));
        // The jobs borrow what stands before them.
        let jobs = Jobs::new();
        let (results, found) = channel::<(u64, Outcome<BlockKeys>)>();
        let mut ordinal = 0u64;
        // Text that is not UTF-8 is reported once the file has been read to its end,
        // where it may fail to be read at all.
        let mut failure: Option<String> = None;
        let mut refused: Option<Error> = None;
        let mut panicked: Option<Box<dyn Any + Send>> = None;
        let read = std::thread::scope(|scope| {
            let _closing = Closing(&jobs);
            for _ in 0..threads {
                scope.spawn(|| jobs.work());
            }
            let reader = {
                let (jobs, buffers, source) = (&jobs, &buffers, &source);
                let (permits, stop, skip) = (&permits, &stop, &skip);
                scope.spawn(move || -> Result<()> {
                    let mut blocks = source.blocks(buffers)?;
                    while !stop.load(Ordering::Relaxed) {
                        permits.take();
                        let next = blocks.next().map_err(|detail| source.failed(&detail))?;
                        let Some(block) = next else {
                            break;
                        };
                        jobs.queue(block.index, results.clone(), move || {
                            let keys = BlockKeys::of(&block, skip);
                            buffers.give(block.into_buffer());
                            keys
                        });
                    }
                    Ok(())
                })
            };
            let collecting = Collecting {
                permits: &permits,
                stop: &stop,
            };
            in_order(found, |outcome: Outcome<BlockKeys>| {
                permits.give();
                let found = match outcome {
                    Ok(found) => found,
                    Err(panic) => {
                        panicked.get_or_insert(panic);
                        stop.store(true, Ordering::Relaxed);
                        return;
                    }
                };
                if failure.is_some() || refused.is_some() {
                    return;
                }
                if found.failure.is_some() {
                    failure = found.failure;
                    skip.store(true, Ordering::Relaxed);
                    return;
                }
                let mut start = 0;
                for &end in &found.ends {
                    let forced = listed.next_if_eq(&ordinal).is_some();
                    if let Err(error) = keys.line(&found.bytes[start..end], forced) {
                        refused = Some(error);
                        stop.store(true, Ordering::Relaxed);
                        return;
                    }
                    ordinal += 1;
                    start = end;
                }
            });
            let read = joined(reader);
            drop(collecting);
            read
        });
        if let Some(panic) = panicked {
            resume_unwind(panic);
        }
        read?;
        match (failure, refused) {
            (Some(detail), _) => Err(source.failed(&detail)),
            (None, Some(error)) => Err(error),
            (None, None) => Ok(keys),
        }
    }

    fn key(&self, group: u32) -> &[u8] {
        let group = group as usize;
        let start = group.checked_sub(1).map_or(0, |before| self.ends[before]);
        &self.bytes[start as usize..self.ends[group] as usize]
    }

    /// The group of `key`, of hash `hashed`, where it is a key of the input.
    fn group(&self, hashed: u64, key: &[u8]) -> Option<u32> {
        let group = *self.ids.get(&hashed)?;
        if self.key(group) == key {
            return Some(group);
        }
        self.others.get(key).copied()
    }

    /// Takes one data line of the input, of the key `key`.
    fn line(&mut self, key: &[u8], forced: bool) -> Result<()> {
        let hashed = (self.hash)(key);
        let group = match self.group(hashed, key) {
            Some(group) => group,
            None => {
                let next = u32::try_from(self.lines.len()).map_err(|_| {
                    Error::comparator("the input holds more keys than can be numbered")
                })?;
                match self.ids.entry(hashed) {
                    Entry::Vacant(slot) => {
                        slot.insert(next);
                    }
                    Entry::Occupied(_) => {
                        self.others.insert(key.to_vec(), next);
                    }
                }
                self.bytes.extend_from_slice(key);
                self.ends.push(self.bytes.len() as u64);
                self.lines.push(0);
                self.forced.push(false);
                next
            }
        };
        let lines = &mut self.lines[group as usize];
        *lines = lines.saturating_add(1);
        self.forced[group as usize] |= forced;
        Ok(())
    }

    /// What tells the group of a line of `key`: the group itself for a key of the input,
    /// the key's hash under `UNKNOWN` for any other.
    pub(crate) fn token(&self, key: &[u8]) -> u64 {
        let hashed = (self.hash)(key);
        match self.group(hashed, key) {
            Some(group) => u64::from(group),
            None => hashed | UNKNOWN,
        }
    }

    /// The input's data lines of the key a line's `token` stands for.
    pub(crate) fn limit(&self, token: u64) -> u64 {
        usize::try_from(token)
            .ok()
            .and_then(|group| self.lines.get(group))
            .map_or(0, |lines| u64::from(*lines))
    }

    /// Whether a listed record is a line of the key a line's `token` stands for.
    pub(crate) fn forced(&self, token: u64) -> bool {
        usize::try_from(token)
            .ok()
            .and_then(|group| self.forced.get(group))
            .is_some_and(|forced| *forced)
    }
}

/// The error of a plan that names what the format `format` does not read.
pub(crate) fn unread(format: &str, what: &str) -> Error {
    Error::comparator(format!("the {format} format reads no {what}"))
}

/// The keys of the input a plan of the VCF or JSON format names.
pub(crate) fn input_keys(plan: &Plan, format: &str, block: usize) -> Result<InputKeys> {
    if !plan.force_locations.is_empty() {
        return Err(unread(format, "listed Locations"));
    }
    if plan.force_text.iter().any(String::is_empty) {
        return Err(Error::comparator(
            "the texts that hand a group over: an empty text matches every line",
        ));
    }
    let Some(input) = &plan.input else {
        return Err(Error::comparator(format!(
            "the {format} format needs the input the outputs were written from"
        )));
    };
    InputKeys::read(input, &plan.force_records, block, plan.threads)
}

/// The error of a plan of the default or tab format that names what the VCF and JSON
/// formats alone read.
pub(crate) fn columns_alone(plan: &Plan) -> Result<()> {
    let format = match plan.format {
        Format::Tab => "tab",
        _ => "default",
    };
    if plan.input.is_some() {
        return Err(unread(format, "input"));
    }
    if !plan.force_records.is_empty() {
        return Err(unread(format, "listed records"));
    }
    if !plan.floor_fields.is_empty() {
        return Err(unread(format, "floors"));
    }
    Ok(())
}

/// The texts of `texts` as finders, for a side whose lines the lists name.
pub(crate) fn listed_texts(texts: &[String], reference: bool) -> Vec<memmem::Finder<'static>> {
    let texts = texts.iter().filter(|_| reference);
    texts
        .map(|text| memmem::Finder::new(text.as_bytes()).into_owned())
        .collect()
}

/// What the comparison names as the groups it must see, beside the input's keys.
struct Force {
    input: InputKeys,
    texts: Vec<String>,
    keep: Option<(Vec<u8>, HashSet<Vec<u8>>)>,
}

pub(crate) struct Vcf {
    force: Arc<Force>,
}

impl Vcf {
    pub(crate) fn new(plan: &Plan, block: usize) -> Result<Vcf> {
        if !plan.floor_fields.is_empty() {
            return Err(unread("vcf", "floors"));
        }
        let input = input_keys(plan, "vcf", block)?;
        let keep = plan.keep_values.as_ref().map(|(column, values)| {
            let values = values.iter().map(|value| value.as_bytes().to_vec());
            (column.as_bytes().to_vec(), values.collect())
        });
        Ok(Vcf {
            force: Arc::new(Force {
                input,
                texts: plan.force_text.clone(),
                keep,
            }),
        })
    }
}

/// A side's header lines before its first data line: every line, and what the
/// `##INFO=<ID=CSQ` lines among them give.
#[derive(Default)]
pub(crate) struct CsqHeader {
    /// The header lines, each ended by a line end.
    text: Vec<u8>,
    /// The Format list of the last CSQ header line that has one, as the line writes it.
    list: Option<Vec<u8>>,
    /// Whether a CSQ header line has no Format list.
    unlisted: bool,
}

impl Grouping for Vcf {
    type Header = CsqHeader;
    type Side = VcfSide;

    const SHAPES: &'static [&'static str] = &[
        "a header line follows a data line",
        "a CSQ header line has no Format list",
        "the CSQ Format list names a field twice",
        "the two outputs' CSQ Format lists differ",
        "a data line has fewer than 8 columns",
        "a CSQ chunk's width is not the Format list's",
        "an Allele, Feature or Feature_type value holds #",
        LISTED_TEXT_SHAPE,
    ];

    const PLACED: Option<u64> = Some(PLACED);

    fn is_header(line: &[u8]) -> bool {
        line.first() == Some(&b'#')
    }

    fn header_line(header: &mut CsqHeader, line: &[u8]) {
        header.text.extend_from_slice(line);
        header.text.push(b'\n');
        if !line.starts_with(CSQ_HEADER) {
            return;
        }
        let Some(at) = memmem::find(line, FORMAT) else {
            header.unlisted = true;
            return;
        };
        let mut list = &line[at + FORMAT.len()..];
        while let Some(shorter) = list.strip_suffix(b"\">") {
            list = shorter;
        }
        header.list = Some(list.to_vec());
    }

    fn side(&self, header: CsqHeader, _data: bool, reference: bool) -> VcfSide {
        let fields: Vec<Vec<u8>> = match &header.list {
            Some(list) => list.split(|&b| b == b'|').map(<[u8]>::to_vec).collect(),
            None => Vec::new(),
        };
        let place = |name: &[u8]| fields.iter().position(|field| field == name);
        let texts = listed_texts(&self.force.texts, reference);
        // The lists name reference lines alone.
        let keep = match &self.force.keep {
            Some((column, _)) if reference => Some(place(column)),
            _ => None,
        };
        VcfSide {
            force: Arc::clone(&self.force),
            // A rewrite of the reference's Format list would rename a field on one side
            // alone.
            listed: texts.iter().any(|text| text.find(&header.text).is_some()),
            texts,
            keep,
            key: KEY_FIELDS.map(place),
            // The comparison reads a chunk into one value a name, the last of a name twice.
            repeated: (1..fields.len()).any(|i| fields[..i].contains(&fields[i])),
            unlisted: header.unlisted,
            fields,
        }
    }

    fn between(&self, expected: &VcfSide, actual: &VcfSide) -> Option<u32> {
        (expected.fields != actual.fields).then_some(FORMAT_LISTS_DIFFER)
    }

    fn trailing(&self, expected: &VcfSide, _actual: &VcfSide) -> Option<u32> {
        expected.listed.then_some(LISTED_TEXT)
    }

    fn limit(&self, key: u64) -> Option<u64> {
        Some(self.force.input.limit(key))
    }
}

pub(crate) struct VcfSide {
    force: Arc<Force>,
    /// On the reference side, the listed texts.
    texts: Vec<memmem::Finder<'static>>,
    /// Whether a header line holds a listed text.
    listed: bool,
    /// On the reference side under a kept column, the place of the column in a chunk:
    /// `None` within where the Format list lacks it, every chunk's value then being the
    /// empty text.
    keep: Option<Option<usize>>,
    /// The Format list in force for every data line: the one of the last CSQ header line,
    /// empty without such a line.
    fields: Vec<Vec<u8>>,
    /// The places of `KEY_FIELDS`, in that order; `None` for a field the list lacks.
    key: [Option<usize>; 3],
    repeated: bool,
    unlisted: bool,
}

/// The ends of a line's first eight columns, a column the line lacks ending where the
/// line does, and how many of the eight it has.
fn columns(text: &[u8]) -> ([usize; 8], usize) {
    let mut ends = [text.len(); 8];
    let mut columns = 1;
    for (end, tab) in ends.iter_mut().zip(memchr_iter(b'\t', text)) {
        *end = tab;
        columns += 1;
    }
    (ends, columns.min(8))
}

fn valued(value: &[u8]) -> bool {
    !value.is_empty() && value != b"-"
}

/// An INFO column's CSQ, the value of its last `CSQ=` item, and whether its other items
/// joined by `;` are a value.
fn read_info(info: &[u8]) -> (Option<&[u8]>, bool) {
    let mut csq = None;
    let mut others = 0usize;
    let mut first: &[u8] = b"";
    let mut start = 0;
    for end in memchr_iter(b';', info).chain(once(info.len())) {
        let item = &info[start..end];
        match item.strip_prefix(CSQ_ITEM) {
            Some(value) => csq = Some(value),
            None => {
                if others == 0 {
                    first = item;
                }
                others += 1;
            }
        }
        start = end + 1;
    }
    (csq, others > 1 || (others == 1 && valued(first)))
}

impl VcfSide {
    /// The chunks of a CSQ: how many, whether one hands the line's group over, and the
    /// shapes they show.
    fn chunks(&self, csq: &[u8]) -> (u32, bool, u32) {
        let width = self.fields.len();
        let marked = memchr(b'#', csq).is_some();
        let kept = self.force.keep.as_ref().zip(self.keep);
        let (mut chunks, mut forced, mut shapes) = (0u32, false, 0u32);
        let mut start = 0;
        for end in memchr_iter(b',', csq).chain(once(csq.len())) {
            let chunk = &csq[start..end];
            start = end + 1;
            chunks = chunks.saturating_add(1);
            if !marked && kept.is_none() {
                if memchr_iter(b'|', chunk).count() + 1 != width {
                    shapes |= 1 << CHUNK_WIDTH;
                }
                continue;
            }
            let mut kept_value: &[u8] = b"";
            let (mut index, mut from) = (0usize, 0usize);
            for bar in memchr_iter(b'|', chunk).chain(once(chunk.len())) {
                let value = &chunk[from..bar];
                // The comparison numbers a repeated key with `#` and a count.
                if marked && self.key.contains(&Some(index)) && value.contains(&b'#') {
                    shapes |= 1 << KEY_MARK;
                }
                if kept.is_some_and(|(_, place)| place == Some(index)) {
                    kept_value = value;
                }
                index += 1;
                from = bar + 1;
            }
            if index != width {
                shapes |= 1 << CHUNK_WIDTH;
            }
            if let Some(((_, values), _)) = kept {
                forced |= !values.contains(kept_value);
            }
        }
        (chunks, forced, shapes)
    }
}

impl Side for VcfSide {
    type Memo = ();

    fn shape(&self, header_after_data: bool) -> Option<u32> {
        if header_after_data {
            Some(HEADER_AFTER_DATA)
        } else if self.unlisted {
            Some(NO_FORMAT_LIST)
        } else if self.repeated {
            Some(REPEATED_FIELD)
        } else {
            None
        }
    }

    fn scan(&self, _memo: &mut (), text: &[u8]) -> Scanned {
        let (ends, columns) = columns(text);
        let key = &text[..ends[4]];
        let token = self.force.input.token(key);
        // A POS of more digits is one the comparison refuses, and a POS is read where
        // the line has a CSQ.
        let position = ends[1] - (ends[0] + 1).min(ends[1]);
        let fragile = past_a_reader_s_limit(text) || position > POSITION_DIGITS;
        let mut scanned = Scanned {
            key: token,
            hash: digest::of_line(text),
            forced: fragile || self.force.input.forced(token),
            fragile,
            shapes: 0,
            strict: false,
            carried: 0,
        };
        // A rewrite of the key's text would move the line to another input record.
        if self.texts.iter().any(|text| text.find(key).is_some()) {
            scanned.shapes |= 1 << LISTED_TEXT;
        }
        scanned.forced |= self.texts.iter().any(|listed| listed.find(text).is_some());
        if columns < 8 {
            scanned.shapes |= 1 << SHORT_LINE;
            return scanned;
        }
        if let (Some(csq), _) = read_info(&text[ends[6] + 1..ends[7]]) {
            let (_, forced, shapes) = self.chunks(csq);
            scanned.forced |= forced;
            scanned.shapes |= shapes;
        }
        scanned
    }

    fn entries(&self, text: &[u8]) -> u32 {
        let (ends, columns) = columns(text);
        if columns < 8 {
            return 1;
        }
        match read_info(&text[ends[6] + 1..ends[7]]).0 {
            Some(csq) => 1 + memchr_iter(b',', csq).count() as u32 + 1,
            None => 1,
        }
    }

    /// The populated counts of kept-back lines, as `populated_counts` counts their
    /// entries: the record entry's fields, and every chunk's value under its field of
    /// the Format list.
    fn populated(
        &self,
        block: &[u8],
        kept: &[Range<usize>],
        tally: &mut Tally,
        entries: &mut Vec<u32>,
    ) {
        let mut record = [0u64; 9];
        let mut fields = vec![0u64; self.fields.len()];
        for line in kept {
            let text = &block[line.clone()];
            let (ends, columns) = columns(text);
            if columns < 8 {
                entries.push(1);
                continue;
            }
            let mut start = 0;
            for (count, &end) in record.iter_mut().zip(&ends[..7]) {
                *count += u64::from(valued(&text[start..end]));
                start = end + 1;
            }
            let (csq, others) = read_info(&text[start..ends[7]]);
            record[7] += u64::from(others);
            // `yes` or `no`, a value either way.
            record[8] += 1;
            let Some(csq) = csq else {
                entries.push(1);
                continue;
            };
            let mut chunks = 1u32;
            let (mut index, mut from) = (0usize, 0usize);
            for end in memchr2_iter(b'|', b',', csq).chain(once(csq.len())) {
                if valued(&csq[from..end]) {
                    if let Some(count) = fields.get_mut(index) {
                        *count += 1;
                    }
                }
                if csq.get(end) == Some(&b',') {
                    index = 0;
                    chunks += 1;
                } else {
                    index += 1;
                }
                from = end + 1;
            }
            entries.push(1 + chunks);
        }
        let record = RECORD_FIELDS.iter().copied().zip(record);
        let fields = self.fields.iter().map(Vec::as_slice).zip(fields);
        for (name, count) in record.chain(fields).filter(|(_, count)| *count > 0) {
            *tally.populated.entry(name.to_vec()).or_insert(0) += count;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::sync::atomic::AtomicU32;

    const CSQ: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. Format: Allele|Consequence|Feature_type|Feature\">";

    static SCRATCH: AtomicU32 = AtomicU32::new(0);

    fn scratch(name: &str, text: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-compare-vcf-{}-{}-{name}",
            std::process::id(),
            SCRATCH.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, text).unwrap();
        path
    }

    fn plan(input: &Path) -> Plan {
        Plan {
            format: Format::Vcf,
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

    /// The side of `text`'s header lines under `plan`, and the text's data lines.
    fn side(plan: &Plan, text: &str, reference: bool) -> (VcfSide, Vec<String>) {
        let format = Vcf::new(plan, 64).unwrap();
        let mut header = CsqHeader::default();
        let mut data = Vec::new();
        for line in text.lines() {
            if data.is_empty() && Vcf::is_header(line.as_bytes()) {
                Vcf::header_line(&mut header, line.as_bytes());
            } else {
                data.push(line.to_string());
            }
        }
        (format.side(header, !data.is_empty(), reference), data)
    }

    fn header_side(text: &str) -> VcfSide {
        let input = scratch("header-input", "");
        let (side, _) = side(&plan(&input), text, true);
        std::fs::remove_file(input).unwrap();
        side
    }

    #[test]
    fn a_key_is_the_first_five_columns_as_text() {
        assert_eq!(line_key(b"1\t7\t.\tA\tG\t.\t.\tx"), b"1\t7\t.\tA\tG");
        assert_eq!(line_key(b"1\t7\t.\tA\tG"), b"1\t7\t.\tA\tG");
        assert_eq!(line_key(b"1\t7\t.\tA\tG\t"), b"1\t7\t.\tA\tG");
        assert_eq!(line_key(b"1\t7"), b"1\t7");
        assert_eq!(line_key(b""), b"");
    }

    #[test]
    fn the_last_csq_header_line_gives_the_format_list() {
        let text = format!("##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Description=\"Format: Old|List\">\n{CSQ}\n#CHROM\n");
        let side = header_side(&text);
        assert_eq!(side.shape(false), None);
        let names: Vec<&[u8]> = vec![b"Allele", b"Consequence", b"Feature_type", b"Feature"];
        assert_eq!(side.fields, names);
        assert_eq!(side.key, [Some(0), Some(3), Some(2)]);
        // Every closing `">` is cut off the list's end, and the first `Format: ` starts it.
        let side = header_side("##INFO=<ID=CSQ,Description=\"Format: A|Format: B\">\">\n");
        let names: Vec<&[u8]> = vec![b"A", b"Format: B"];
        assert_eq!(side.fields, names);
        let side = header_side("##INFO=<ID=CSQ,Description=\"Format: \">\n");
        assert_eq!(side.fields, vec![b""]);
        let side = header_side("##fileformat=VCFv4.2\n");
        assert!(side.fields.is_empty() && side.shape(false).is_none());
    }

    #[test]
    fn header_lines_the_comparison_reads_otherwise_are_named() {
        let shape = |text: &str, header_after_data: bool| {
            let shape = header_side(text).shape(header_after_data);
            shape.map(|shape| Vcf::SHAPES[shape as usize])
        };
        assert_eq!(
            shape("##INFO=<ID=CSQ,Description=\"none\">\n", false),
            Some("a CSQ header line has no Format list")
        );
        assert_eq!(
            shape("##INFO=<ID=CSQ,Description=\"Format: A|B|A\">\n", false),
            Some("the CSQ Format list names a field twice")
        );
        assert_eq!(
            shape(&format!("{CSQ}\n"), true),
            Some("a header line follows a data line")
        );
        // A list an earlier line repeats a name in is not the list in force.
        assert_eq!(
            shape(
                &format!("##INFO=<ID=CSQ,Description=\"Format: A|A\">\n{CSQ}\n"),
                false
            ),
            None
        );
    }

    #[test]
    fn the_last_csq_item_is_read_and_the_other_items_are_one_field() {
        assert_eq!(read_info(b"CSQ=a|b"), (Some(&b"a|b"[..]), false));
        assert_eq!(read_info(b"AC=1;CSQ=a;AN=2"), (Some(&b"a"[..]), true));
        assert_eq!(read_info(b"CSQ=a;CSQ=b"), (Some(&b"b"[..]), false));
        assert_eq!(read_info(b"CSQ="), (Some(&b""[..]), false));
        assert_eq!(read_info(b"."), (None, true));
        assert_eq!(read_info(b"-"), (None, false));
        assert_eq!(read_info(b""), (None, false));
        assert_eq!(read_info(b";CSQ=a"), (Some(&b"a"[..]), false));
        assert_eq!(read_info(b"-;CSQ=a"), (Some(&b"a"[..]), false));
        // Two empty items join to `;`, a value.
        assert_eq!(read_info(b";"), (None, true));
        assert_eq!(read_info(b"csq=a"), (None, true));
    }

    /// The populated counts and the entries of the data lines of `text`, counted in two
    /// calls.
    fn counts(text: &str) -> (Vec<(String, u64)>, Vec<u32>) {
        let input = scratch("counts-input", "");
        let (side, _) = side(&plan(&input), text, true);
        std::fs::remove_file(input).unwrap();
        let kept: Vec<Range<usize>> = lines(text.as_bytes())
            .filter(|line| !Vcf::is_header(&text.as_bytes()[line.clone()]))
            .collect();
        let mut tally = Tally::default();
        let mut entries = Vec::new();
        let (first, second) = kept.split_at(kept.len() / 2);
        side.populated(text.as_bytes(), first, &mut tally, &mut entries);
        side.populated(text.as_bytes(), second, &mut tally, &mut entries);
        assert!(tally.present.is_empty() && tally.under.is_empty());
        let of_lines: Vec<u32> = kept
            .iter()
            .map(|line| side.entries(&text.as_bytes()[line.clone()]))
            .collect();
        assert_eq!(entries, of_lines);
        let populated = tally.populated.into_iter();
        let named = populated.map(|(name, n)| (String::from_utf8(name).unwrap(), n));
        (named.collect(), entries)
    }

    fn pairs(expected: &[(&str, u64)]) -> Vec<(String, u64)> {
        expected
            .iter()
            .map(|(n, c)| ((*n).to_string(), *c))
            .collect()
    }

    #[test]
    fn a_line_counts_its_record_entry_and_every_chunk() {
        let text = format!(
            "{CSQ}\n\
             1\t7\trs1\tA\tG\t.\t-\tAC=1;CSQ=G|missense|Transcript|tx1,G||Transcript|tx2\n\
             1\t8\t\tA\tG\t\t.\tCSQ=G|-|-|\n\
             1\t9\t.\tA\tG\t.\t.\t.\tGT\t0/1\n\
             1\t9\t.\tA\n"
        );
        let (populated, entries) = counts(&text);
        assert_eq!(
            populated,
            pairs(&[
                ("ALT", 3),
                ("Allele", 3),
                ("CHROM", 3),
                ("CSQ_present", 3),
                ("Consequence", 1),
                ("FILTER", 2),
                ("Feature", 2),
                ("Feature_type", 2),
                ("ID", 2),
                ("INFO_without_CSQ", 2),
                ("POS", 3),
                ("QUAL", 2),
                ("REF", 3),
            ])
        );
        assert_eq!(entries, [3, 2, 1, 1]);
    }

    #[test]
    fn a_csq_field_named_as_a_record_field_adds_to_one_count() {
        let text = "##INFO=<ID=CSQ,Description=\"Format: Allele|ID\">\n\
                    1\t7\trs1\tA\tG\t.\t.\tCSQ=G|x,G|y\n";
        assert_eq!(
            counts(text).0,
            pairs(&[
                ("ALT", 1),
                ("Allele", 2),
                ("CHROM", 1),
                ("CSQ_present", 1),
                ("FILTER", 1),
                ("ID", 3),
                ("POS", 1),
                ("QUAL", 1),
                ("REF", 1),
            ])
        );
    }

    /// The keys of the data lines of `text`, a file read through blocks of `block` bytes,
    /// under `hashed`.
    fn keys_of(text: &str, records: &[u64], block: usize, hashed: fn(&[u8]) -> u64) -> InputKeys {
        let path = scratch("keys", text);
        let threads = 1 + block % 3;
        let keys = InputKeys::hashed(&path, records, block, threads, hashed).unwrap();
        std::fs::remove_file(path).unwrap();
        keys
    }

    #[test]
    fn keys_of_one_hash_are_groups_of_their_own() {
        let text = "#CHROM\n1\t7\ta\tA\tG\tq\n1\t8\tb\tA\tG\n1\t7\ta\tA\tG\tr\n2\t70\tc\tA\tG\n";
        // A hash by length: the first three keys share one, the last has its own.
        let by_length: fn(&[u8]) -> u64 = |key| key.len() as u64;
        for hashed in [hash as fn(&[u8]) -> u64, by_length] {
            for block in [1, 16, 1 << 20] {
                let keys = keys_of(text, &[2, 3], block, hashed);
                assert_eq!(keys.lines, [2, 1, 1]);
                assert_eq!(keys.forced, [true, false, true]);
                assert_eq!(keys.token(b"1\t7\ta\tA\tG"), 0);
                assert_eq!(keys.token(b"1\t8\tb\tA\tG"), 1);
                assert_eq!(keys.token(b"2\t70\tc\tA\tG"), 2);
                assert_eq!((keys.limit(0), keys.limit(1), keys.limit(2)), (2, 1, 1));
                assert!(keys.forced(0) && !keys.forced(1) && keys.forced(2));
                for absent in [&b"1\t9\tb\tA\tG"[..], b"1\t7\ta\tA", b""] {
                    let token = keys.token(absent);
                    assert_eq!(token, hashed(absent) | UNKNOWN);
                    assert_eq!((keys.limit(token), keys.forced(token)), (0, false));
                }
                assert_eq!((keys.limit(PLACED), keys.forced(PLACED)), (0, false));
            }
        }
    }

    #[test]
    fn the_keys_of_a_long_input_are_numbered_as_it_first_names_them_at_any_worker_count() {
        // 5,000 lines of 613 keys: line `n` has the key of position `7 n` modulo 613.
        let position = |line: usize| (line * 7) % 613;
        let mut text = String::from("##a header line\n");
        for line in 0..5000 {
            text.push_str(&format!("1\t{}\trs\tA\tG\tq{line}\n", position(line)));
        }
        let records: [u64; 4] = [0, 612, 613, 4999];
        let path = scratch("long", &text);
        for (block, threads) in [(16, 1), (300, 4), (1 << 20, 8)] {
            let keys = InputKeys::read(&path, &records, block, threads).unwrap();
            assert_eq!(keys.lines.len(), 613);
            assert_eq!(keys.lines.iter().map(|n| u64::from(*n)).sum::<u64>(), 5000);
            for group in [0, 1, 100, 612] {
                let key = format!("1\t{}\trs\tA\tG", position(group));
                assert_eq!(keys.key(group as u32), key.as_bytes());
                assert_eq!(keys.token(key.as_bytes()), group as u64);
            }
            // The group of a listed line is the first line of its key.
            let listed: Vec<usize> = records
                .iter()
                .map(|line| position(*line as usize))
                .map(|listed| (0..613).find(|first| position(*first) == listed).unwrap())
                .collect();
            for group in 0..613 {
                assert_eq!(keys.forced[group], listed.contains(&group), "group {group}");
            }
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn an_input_is_its_data_lines_whatever_ends_them() {
        let keys = keys_of("##a\r\n1\t7\r\n\r\n#b\r1\t7\n", &[1], 4, hash);
        // The empty line is a data line, of the empty key.
        assert_eq!(keys.lines, [2, 1]);
        assert_eq!(keys.forced, [false, true]);
        assert_eq!(keys.token(b""), 1);
        let absent = Path::new("/nonexistent/input.vcf");
        let error = InputKeys::read(absent, &[], 64, 2).err().unwrap();
        assert_eq!(
            error.to_string(),
            "cannot read /nonexistent/input.vcf: No such file or directory"
        );
        let path = scratch("not-utf8", "1\t7\n");
        std::fs::write(&path, b"1\t7\n\xff\n").unwrap();
        let error = InputKeys::read(&path, &[], 2, 2).err().unwrap();
        std::fs::remove_file(&path).unwrap();
        assert!(error
            .to_string()
            .ends_with("'utf-8' codec can't decode byte 0xff in position 4: invalid start byte"));
    }

    #[test]
    fn a_line_s_scan_reads_its_key_its_chunks_and_what_hands_it_over() {
        let input = scratch(
            "scan-input",
            "1\t7\trs1\tA\tG\n1\t8\trs2\tA\tG\n1\t8\trs2\tA\tG\n",
        );
        let mut plan = plan(&input);
        plan.force_records = vec![0];
        plan.force_text = vec!["misspelt".to_string()];
        plan.keep_values = Some((
            "Feature_type".to_string(),
            vec!["Transcript".to_string(), String::new()],
        ));
        let text = format!("{CSQ}\n");
        let (reference, _) = side(&plan, &text, true);
        let (candidate, _) = side(&plan, &text, false);
        std::fs::remove_file(input).unwrap();
        let scan = |side: &VcfSide, line: &str| side.scan(&mut (), line.as_bytes());
        let listed = scan(&reference, "1\t7\trs1\tA\tG\t.\t.\tCSQ=G|x|Transcript|tx1");
        assert!(listed.forced && !listed.fragile);
        assert_eq!((listed.key, listed.shapes), (0, 0));
        let plain = "1\t8\trs2\tA\tG\t.\t.\tCSQ=G|x|Transcript|tx1,G|x||tx2";
        let kept = scan(&reference, plain);
        assert!(!kept.forced);
        assert_eq!((kept.key, kept.shapes), (1, 0));
        assert_eq!(kept.hash, scan(&candidate, plain).hash);
        // A value off the kept list, and a listed text, name reference lines alone.
        let regulatory = "1\t8\trs2\tA\tG\t.\t.\tCSQ=G|x|RegulatoryFeature|reg1";
        assert!(scan(&reference, regulatory).forced && !scan(&candidate, regulatory).forced);
        let named = "1\t8\trs2\tA\tG\t.\t.\tCSQ=G|misspelt7|Transcript|tx1";
        assert!(scan(&reference, named).forced && !scan(&candidate, named).forced);
        // A key the input lacks.
        let absent = scan(&reference, "1\t9\trs3\tA\tG\t.\t.\t.");
        assert_eq!(absent.key, hash(b"1\t9\trs3\tA\tG") | UNKNOWN);
        assert!(!absent.forced);
        let shapes = |line: &str| scan(&reference, line).shapes;
        assert_eq!(shapes("1\t8\trs2\tA\tG\t.\t."), 1 << SHORT_LINE);
        assert_eq!(shapes(""), 1 << SHORT_LINE);
        assert_eq!(
            shapes("1\t8\trs2\tA\tG\t.\t.\tCSQ=G|x|Transcript"),
            1 << CHUNK_WIDTH
        );
        assert_eq!(
            shapes("1\t8\trs2\tA\tG\t.\t.\tCSQ=G|x|Transcript|tx#2"),
            1 << KEY_MARK
        );
        assert_eq!(shapes("1\t8\trs2\tA\tG\t.\t.\tCSQ=G|x#2|Transcript|tx"), 0);
        assert_eq!(shapes("1\t8\tmisspelt2\tA\tG\t.\t.\t."), 1 << LISTED_TEXT);
        // A line past a reader's limit, and a POS past the digits of a whole number.
        let fragile = scan(&candidate, "1\t8\trs2\tA\tG\t.\t.\tAC=\u{b2}");
        assert!(fragile.forced && fragile.fragile);
        let long = format!("1\t{}\trs2\tA\tG\t.\t.\t.", "7".repeat(19));
        assert!(scan(&candidate, &long).fragile);
        let short = format!("1\t{}\trs2\tA\tG\t.\t.\t.", "7".repeat(18));
        assert!(!scan(&candidate, &short).fragile);
    }

    #[test]
    fn a_plan_names_what_the_format_reads() {
        let input = scratch("plan-input", "1\t7\trs1\tA\tG\n");
        let refused = |plan: &Plan| Vcf::new(plan, 64).err().unwrap().to_string();
        let mut lacking = plan(&input);
        lacking.input = None;
        assert_eq!(
            refused(&lacking),
            "the vcf format needs the input the outputs were written from"
        );
        let mut located = plan(&input);
        located.force_locations = vec!["1:7".to_string()];
        assert_eq!(
            refused(&located),
            "the vcf format reads no listed Locations"
        );
        let mut floored = plan(&input);
        floored.floor_fields = vec![("IMPACT".to_string(), vec!["impact".to_string()])];
        assert_eq!(refused(&floored), "the vcf format reads no floors");
        let mut empty = plan(&input);
        empty.force_text = vec![String::new()];
        assert_eq!(
            refused(&empty),
            "the texts that hand a group over: an empty text matches every line"
        );
        let mut columns = plan(&input);
        columns.format = Format::Tab;
        assert_eq!(
            columns_alone(&columns).err().unwrap().to_string(),
            "the tab format reads no input"
        );
        std::fs::remove_file(input).unwrap();
    }
}
