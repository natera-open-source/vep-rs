// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The transcripts of a candidate JSON cache as the row rules of a coding end ask them: for
//! the transcript of a row, the stop codon and the start codon vep-rs reads at the two ends
//! of its coding sequence (`CodingEnds`).
//!
//! A cache files its transcripts under `transcripts/<sequence>/`, one shard a file named
//! `*.json` or, gzipped, `*.json.gz`; a shard is an array of transcript objects or an object
//! of such arrays, a `null` among them standing for no transcript. The shards of one
//! sequence are read when a row on that sequence first asks, and no other. A record is read
//! into the facts `coding_ends` asks (`RecordFacts`) as vep-rs loads it, and a record that
//! states one of them in a form vep-rs reads otherwise or refuses is not read
//! (`StopEnd::Unread`). `read_shard` says which file under a shard's name is passed over and
//! which ends the reading with an error that names it.
//!
//! The reader is written for the three tool roots that compile this file, so it asks of
//! them nothing but the gzip reader and the byte search they share.

use super::coding_ends::{coding_ends, CodingEnds, MapperFacts, RecordFacts};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

/// One value of a shard's JSON text. A number is kept as the text writes it.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Json {
    Null,
    Bool(bool),
    Number(String),
    Text(String),
    List(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    /// The member `key` of an object, the last where the text repeats the key; `None` for
    /// an object without it and for any other value.
    pub(crate) fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(members) => members
                .iter()
                .rev()
                .find(|(name, _)| name == key)
                .map(|(_, value)| value),
            _ => None,
        }
    }

    fn as_text(&self) -> Option<&str> {
        match self {
            Json::Text(text) => Some(text),
            _ => None,
        }
    }

    fn as_list(&self) -> Option<&[Json]> {
        match self {
            Json::List(items) => Some(items),
            _ => None,
        }
    }

    /// An integer of a cache transcript, a number in a native cache and its text in a
    /// converted one: an optional `-` and the digits `0` to `9`, within 64 bits. `None` for
    /// any other value.
    pub(crate) fn integer(&self) -> Option<i64> {
        let digits = match self {
            Json::Number(digits) | Json::Text(digits) => digits,
            _ => return None,
        };
        let unsigned = digits.strip_prefix('-').unwrap_or(digits);
        if unsigned.is_empty() || !unsigned.bytes().all(|byte| byte.is_ascii_digit()) {
            return None;
        }
        digits.parse().ok()
    }
}

/// The members of a transcript object the reading asks; every other member is passed over
/// without being kept.
const RECORD_MEMBERS: [&str; 12] = [
    "stable_id",
    "start",
    "end",
    "strand",
    "ccds",
    "cdna_coding_start",
    "cdna_coding_end",
    "attributes",
    "exons",
    "dbID",
    "version",
    "source_cache",
];
const FEATURE_CACHE: &str = "variation_effect_feature_cache";
const FEATURE_CACHE_MEMBERS: [&str; 5] = [
    "translateable_seq",
    "peptide",
    "codon_table",
    "mapper",
    "introns",
];

/// Which members of an object are kept.
#[derive(Clone, Copy, PartialEq)]
enum Keep {
    All,
    Record,
    FeatureCache,
}

/// The deepest nesting of arrays and objects a shard is read to.
const DEEPEST: usize = 64;

struct Parser<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Parser<'a> {
    fn refused(&self, what: &str) -> String {
        format!("{what} at byte {}", self.at)
    }

    fn blanks(&mut self) {
        while matches!(self.bytes.get(self.at), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.at += 1;
        }
    }

    fn expect(&mut self, byte: u8) -> Result<(), String> {
        if self.bytes.get(self.at) == Some(&byte) {
            self.at += 1;
            Ok(())
        } else {
            Err(self.refused(&format!("no `{}`", char::from(byte))))
        }
    }

    fn word(&mut self, word: &str, value: Json) -> Result<Json, String> {
        if self.bytes[self.at..].starts_with(word.as_bytes()) {
            self.at += word.len();
            Ok(value)
        } else {
            Err(self.refused("no JSON value"))
        }
    }

    /// The four hexadecimal digits of a `\u` escape.
    fn code_unit(&mut self) -> Result<u32, String> {
        let digits = self
            .bytes
            .get(self.at..self.at + 4)
            .and_then(|digits| std::str::from_utf8(digits).ok())
            .and_then(|digits| u32::from_str_radix(digits, 16).ok())
            .ok_or_else(|| self.refused("no four hexadecimal digits"))?;
        self.at += 4;
        Ok(digits)
    }

    /// A string, the opening quote at the cursor; `None` where `keep` is false, the string
    /// passed over.
    fn string(&mut self, keep: bool) -> Result<Option<String>, String> {
        self.expect(b'"')?;
        let mut text: Vec<u8> = Vec::new();
        loop {
            let rest = &self.bytes[self.at..];
            let Some(stop) = memchr::memchr2(b'"', b'\\', rest) else {
                return Err(self.refused("a string without its closing quote"));
            };
            if keep {
                text.extend_from_slice(&rest[..stop]);
            }
            self.at += stop + 1;
            if rest[stop] == b'"' {
                break;
            }
            let Some(&escape) = self.bytes.get(self.at) else {
                return Err(self.refused("a string without its closing quote"));
            };
            self.at += 1;
            let unit = match escape {
                b'"' | b'\\' | b'/' => u32::from(escape),
                b'b' => 0x08,
                b'f' => 0x0c,
                b'n' => 0x0a,
                b'r' => 0x0d,
                b't' => 0x09,
                b'u' => self.code_unit()?,
                _ => return Err(self.refused("no escape of a string")),
            };
            let scalar =
                if (0xd800..0xdc00).contains(&unit) && self.bytes[self.at..].starts_with(b"\\u") {
                    self.at += 2;
                    let low = self.code_unit()?;
                    if !(0xdc00..0xe000).contains(&low) {
                        return Err(self.refused("half a surrogate pair"));
                    }
                    0x10000 + ((unit - 0xd800) << 10) + (low - 0xdc00)
                } else {
                    unit
                };
            let Some(scalar) = char::from_u32(scalar) else {
                return Err(self.refused("half a surrogate pair"));
            };
            if keep {
                text.extend_from_slice(scalar.encode_utf8(&mut [0; 4]).as_bytes());
            }
        }
        if !keep {
            return Ok(None);
        }
        String::from_utf8(text)
            .map(Some)
            .map_err(|_| self.refused("a string that is not UTF-8"))
    }

    fn number(&mut self) -> Result<Json, String> {
        let start = self.at;
        while matches!(
            self.bytes.get(self.at),
            Some(b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9')
        ) {
            self.at += 1;
        }
        let digits = &self.bytes[start..self.at];
        if !digits.iter().any(u8::is_ascii_digit) {
            return Err(self.refused("no JSON value"));
        }
        Ok(Json::Number(String::from_utf8_lossy(digits).into_owned()))
    }

    /// The value at the cursor, its objects keeping the members `keep` names; `None` where
    /// `kept` is false, the value passed over.
    fn value(&mut self, keep: Keep, kept: bool, depth: usize) -> Result<Option<Json>, String> {
        if depth > DEEPEST {
            return Err(self.refused("nesting past the depth a shard is read to"));
        }
        self.blanks();
        let value = match self.bytes.get(self.at) {
            Some(b'"') => self.string(kept)?.map(Json::Text),
            Some(b'{') => {
                self.at += 1;
                let mut members = Vec::new();
                self.blanks();
                if self.bytes.get(self.at) == Some(&b'}') {
                    self.at += 1;
                } else {
                    loop {
                        self.blanks();
                        let name = self.string(true)?.unwrap_or_default();
                        self.blanks();
                        self.expect(b':')?;
                        let (inner, wanted) = match keep {
                            Keep::All => (Keep::All, kept),
                            Keep::Record if name == FEATURE_CACHE => (Keep::FeatureCache, kept),
                            Keep::Record => (Keep::All, kept && RECORD_MEMBERS.contains(&&*name)),
                            Keep::FeatureCache => {
                                (Keep::All, kept && FEATURE_CACHE_MEMBERS.contains(&&*name))
                            }
                        };
                        if let Some(value) = self.value(inner, wanted, depth + 1)? {
                            members.push((name, value));
                        }
                        self.blanks();
                        match self.bytes.get(self.at) {
                            Some(b',') => self.at += 1,
                            Some(b'}') => {
                                self.at += 1;
                                break;
                            }
                            _ => return Err(self.refused("an object without its end")),
                        }
                    }
                }
                kept.then_some(Json::Object(members))
            }
            Some(b'[') => {
                self.at += 1;
                let mut items = Vec::new();
                self.blanks();
                if self.bytes.get(self.at) == Some(&b']') {
                    self.at += 1;
                } else {
                    loop {
                        // An array hands its own kind of keeping to its items: the records
                        // of a shard are the items of its arrays.
                        if let Some(item) = self.value(keep, kept, depth + 1)? {
                            items.push(item);
                        }
                        self.blanks();
                        match self.bytes.get(self.at) {
                            Some(b',') => self.at += 1,
                            Some(b']') => {
                                self.at += 1;
                                break;
                            }
                            _ => return Err(self.refused("an array without its end")),
                        }
                    }
                }
                kept.then_some(Json::List(items))
            }
            Some(b'n') => Some(self.word("null", Json::Null)?),
            Some(b't') => Some(self.word("true", Json::Bool(true))?),
            Some(b'f') => Some(self.word("false", Json::Bool(false))?),
            Some(_) => Some(self.number()?),
            None => return Err(self.refused("no JSON value")),
        };
        Ok(value.filter(|_| kept))
    }
}

/// The transcript objects of a shard's text, each reduced to the members the reading asks:
/// the items of a bare array, or of every array of an object of arrays, a `null` dropped.
pub(crate) fn shard_records(text: &str) -> Result<Vec<Json>, String> {
    let mut parser = Parser {
        bytes: text.as_bytes(),
        at: 0,
    };
    parser.blanks();
    let records = match parser.bytes.get(parser.at) {
        Some(b'[') => parser.value(Keep::Record, true, 0)?,
        Some(b'{') => {
            // The members of the object are the arrays; none of them is passed over.
            parser.at += 1;
            let mut items = Vec::new();
            parser.blanks();
            if parser.bytes.get(parser.at) == Some(&b'}') {
                parser.at += 1;
            } else {
                loop {
                    parser.blanks();
                    parser.string(false)?;
                    parser.blanks();
                    parser.expect(b':')?;
                    parser.blanks();
                    if parser.bytes.get(parser.at) != Some(&b'[') {
                        return Err("a member of the shard's object is no array".to_string());
                    }
                    match parser.value(Keep::Record, true, 1)? {
                        Some(Json::List(array)) => items.extend(array),
                        _ => return Err("a member of the shard's object is no array".to_string()),
                    }
                    parser.blanks();
                    match parser.bytes.get(parser.at) {
                        Some(b',') => parser.at += 1,
                        Some(b'}') => {
                            parser.at += 1;
                            break;
                        }
                        _ => return Err(parser.refused("an object without its end")),
                    }
                }
            }
            Some(Json::List(items))
        }
        _ => {
            return Err("neither an array of transcripts nor an object of such arrays".to_string())
        }
    };
    parser.blanks();
    if parser.at != parser.bytes.len() {
        return Err(parser.refused("text after the shard's value"));
    }
    let Some(Json::List(items)) = records else {
        return Err("neither an array of transcripts nor an object of such arrays".to_string());
    };
    let records: Vec<Json> = items
        .into_iter()
        .filter(|record| *record != Json::Null)
        .collect();
    if !records
        .iter()
        .all(|record| matches!(record, Json::Object(_)))
    {
        return Err("a transcript is not an object".to_string());
    }
    Ok(records)
}

/// A member that is text, `Some(None)` for one the record lacks or holds as `null`; `None`
/// for any other value, which vep-rs refuses.
fn text(member: Option<&Json>) -> Option<Option<&str>> {
    match member {
        None | Some(Json::Null) => Some(None),
        Some(Json::Text(text)) => Some(Some(text)),
        Some(_) => None,
    }
}

/// The record's `start` and `end`, `None` where it states either in no form an integer is
/// read from.
fn record_span(record: &Json) -> Option<(i64, i64)> {
    Some((
        record.get("start").and_then(Json::integer)?,
        record.get("end").and_then(Json::integer)?,
    ))
}

/// The mapper pairs and the exons of a record that states its strand, its two coding bounds
/// on the cDNA, at least one pair of five whole numbers and at least one exon of two; `None`
/// for any other record. The exons are in genomic order.
type Coordinates = (Vec<(i64, i64, i64, i64, i64)>, Vec<(i64, i64)>);
fn coordinates(record: &Json, mapper: &Json) -> Option<Coordinates> {
    let integer = |item: &Json, key: &str| item.get(key).and_then(Json::integer);
    let strand = integer(record, "strand")?;
    integer(record, "cdna_coding_start")?;
    integer(record, "cdna_coding_end")?;
    if !matches!(strand, 1 | -1) {
        return None;
    }
    let pairs = mapper.get("pairs").and_then(Json::as_list)?;
    let exons = record.get("exons").and_then(Json::as_list)?;
    if pairs.is_empty() || exons.is_empty() {
        return None;
    }
    let pairs = pairs
        .iter()
        .map(|pair| {
            Some((
                integer(pair, "from_start")?,
                integer(pair, "from_end")?,
                integer(pair, "to_start")?,
                integer(pair, "to_end")?,
                integer(pair, "ori")?,
            ))
        })
        .collect::<Option<Vec<_>>>()?;
    let mut exons = exons
        .iter()
        .map(|exon| Some((integer(exon, "start")?, integer(exon, "end")?)))
        .collect::<Option<Vec<_>>>()?;
    exons.sort_unstable();
    Some((pairs, exons))
}

/// The `codon_table` of a record's feature cache: `Some(None)` where it states none,
/// `Some(Some(id))` for a number or a text of the digits `0` to `9` alone, `None` for any
/// other value, which names no genetic code.
fn codon_table(member: Option<&Json>) -> Option<Option<i64>> {
    match member {
        None | Some(Json::Null) => Some(None),
        Some(Json::Number(digits) | Json::Text(digits))
            if !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit()) =>
        {
            digits.parse().ok().map(Some)
        }
        Some(_) => None,
    }
}

/// The facts of a record whose coding sequence is `coding`, read into its coding ends;
/// `None` for a record that states one of them in a form vep-rs reads otherwise or refuses:
/// an attribute whose `code` or `value` is no text, a `ccds` or a `peptide` that is neither
/// text nor `null`, a `codon_table` that names no genetic code, a `start_phase` or a coding
/// bound of the mapper that is no small whole number, a mapper beside which the record
/// states no coordinates (`coordinates`), or `introns` that are no list of two whole numbers
/// each. The introns are the cached ones, or where the record caches none the stretches
/// between its exons. The record's 3' UTR is not handed to the reading: this comparison
/// derives no protein notation, so no transcript's sequence is kept for one.
fn read_coding_ends(record: &Json, cached: &Json, coding: &str) -> Option<CodingEnds> {
    let attribute_codes: Vec<&str> = match record.get("attributes") {
        None => Vec::new(),
        Some(attributes) => attributes
            .as_list()?
            .iter()
            .map(|attribute| {
                attribute.get("value")?.as_text()?;
                attribute.get("code")?.as_text()
            })
            .collect::<Option<_>>()?,
    };
    let codon_table = codon_table(cached.get("codon_table"))?;
    let mut introns: Vec<(i64, i64)> = match cached.get("introns") {
        None | Some(Json::Null) => Vec::new(),
        Some(introns) => introns
            .as_list()?
            .iter()
            .map(|intron| {
                Some((
                    intron.get("start").and_then(Json::integer)?,
                    intron.get("end").and_then(Json::integer)?,
                ))
            })
            .collect::<Option<_>>()?,
    };
    let pairs;
    let mapper = match cached.get("mapper") {
        None | Some(Json::Null) => None,
        Some(mapper) => {
            if !matches!(mapper, Json::Object(_)) {
                return None;
            }
            let bound = |name: &str| match mapper.get(name) {
                None | Some(Json::Null) => Some(0),
                Some(bound) => bound.integer().filter(|bound| *bound >= 0),
            };
            let start_phase = match mapper.get("start_phase") {
                None => 0,
                Some(phase) => phase
                    .integer()
                    .filter(|phase| i8::try_from(*phase).is_ok())?,
            };
            let exons;
            (pairs, exons) = coordinates(record, mapper)?;
            if introns.is_empty() {
                introns = exons
                    .windows(2)
                    .filter_map(|exons| {
                        Some((exons[0].1.checked_add(1)?, exons[1].0.checked_sub(1)?))
                    })
                    .collect();
            }
            Some(MapperFacts {
                start_phase,
                coding_start: bound("cdna_coding_start")?,
                coding_end: bound("cdna_coding_end")?,
                pairs: &pairs,
            })
        }
    };
    Some(coding_ends(&RecordFacts {
        attribute_codes: &attribute_codes,
        consensus_id: text(record.get("ccds"))?,
        strand: record
            .get("strand")
            .and_then(Json::integer)
            .filter(|strand| matches!(strand, 1 | -1)),
        codon_table,
        coding: Some(coding),
        peptide: text(cached.get("peptide"))?,
        mapper,
        span: record_span(record),
        introns: &introns,
    }))
}

/// The stable id of one transcript object and its coding ends, read as vep-rs loads the
/// record; `None` for a transcript without a stable id that is text. A transcript without a
/// coding sequence that is text is not read, nor is one whose record `read_coding_ends`
/// refuses; a record that is not read keeps its span.
pub(crate) fn coding_ends_of(record: &Json) -> Option<(String, CodingEnds)> {
    let stable_id = record.get("stable_id")?.as_text()?.to_string();
    let cached = record
        .get(FEATURE_CACHE)
        .filter(|cached| matches!(cached, Json::Object(_)));
    let coding = cached
        .and_then(|cached| cached.get("translateable_seq"))
        .and_then(Json::as_text);
    let unread = |tail_length: u8| CodingEnds {
        span: record_span(record),
        ..CodingEnds::unread(tail_length)
    };
    let ends = match (cached, coding) {
        (Some(cached), Some(coding)) => read_coding_ends(record, cached, coding)
            .unwrap_or_else(|| unread((coding.chars().count() % 3) as u8)),
        _ => unread(0),
    };
    Some((stable_id, ends))
}

/// The bytes of a shard, a name that ends in `.gz` read through gzip.
fn shard_bytes(shard: &Path) -> io::Result<Vec<u8>> {
    let mut bytes = Vec::new();
    let file = std::fs::File::open(shard)?;
    if shard.extension().is_some_and(|extension| extension == "gz") {
        MultiGzDecoder::new(file).read_to_end(&mut bytes)?;
    } else {
        let mut file = file;
        file.read_to_end(&mut bytes)?;
    }
    Ok(bytes)
}

/// Whether `text` is one JSON array or object and nothing after it.
fn is_array_or_object(text: &str) -> bool {
    let mut parser = Parser {
        bytes: text.as_bytes(),
        at: 0,
    };
    parser.blanks();
    let opens = matches!(parser.bytes.get(parser.at), Some(b'[' | b'{'));
    let read = opens && parser.value(Keep::All, false, 0).is_ok();
    parser.blanks();
    read && parser.at == parser.bytes.len()
}

/// The transcript objects of the file at `shard`. A file under a shard's name that is no
/// shard holds none and is passed over, since a copied cache can hold one beside its
/// shards: a file that cannot be read whole, bytes that are not UTF-8 and a text that is no
/// JSON array or object. A JSON array or object that is neither an array of transcript
/// objects nor an object of such arrays is a failure that names the file.
pub(crate) fn read_shard(shard: &Path) -> io::Result<Vec<Json>> {
    let Ok(bytes) = shard_bytes(shard) else {
        return Ok(Vec::new());
    };
    let Ok(text) = std::str::from_utf8(&bytes) else {
        return Ok(Vec::new());
    };
    match shard_records(text) {
        Ok(records) => Ok(records),
        Err(reason) if is_array_or_object(text) => Err(io::Error::other(format!(
            "{}: is no shard of transcripts: {reason}",
            shard.display()
        ))),
        Err(_) => Ok(Vec::new()),
    }
}

/// The shards of the directory `directory` in the order of their names, the `*.json` before
/// the `*.json.gz`; a directory that cannot be listed has none.
pub(crate) fn shards_of(directory: &Path) -> Vec<PathBuf> {
    let mut names: Vec<std::ffi::OsString> = std::fs::read_dir(directory)
        .map(|entries| entries.flatten().map(|entry| entry.file_name()).collect())
        .unwrap_or_default();
    names.sort();
    let suffixed = |suffix: &'static str| {
        names
            .iter()
            .filter(move |name| name.to_string_lossy().ends_with(suffix))
            .map(|name| directory.join(name))
    };
    suffixed(".json").chain(suffixed(".json.gz")).collect()
}

/// The coding ends of each record of one stable id, a reading that two records share once:
/// a record that two shards repeat is one, a copy of the transcript on another span or
/// another transcript under the same id is a second.
type Records = HashMap<String, Vec<Arc<CodingEnds>>>;

/// What reading one shard came to, with the shard's place in the order of the shards: the
/// stable id and the coding ends of each of its transcripts, or the failure that names it.
type ShardRead = (usize, io::Result<Vec<(String, CodingEnds)>>);

/// The lowest and the highest genomic base a Location names after its sequence name, an
/// insertion's being the two bases it stands between; `None` for a text that names no
/// position there.
pub fn location_span(location: &str) -> Option<(i64, i64)> {
    let (_, position) = location.rsplit_once(':')?;
    let digits = |text: &str| -> Option<i64> {
        (!text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit()))
            .then(|| text.parse().ok())
            .flatten()
    };
    let (first, last) = match position.split_once('-') {
        Some((first, last)) => (digits(first)?, digits(last)?),
        None => {
            let base = digits(position)?;
            (base, base)
        }
    };
    Some((first.min(last), first.max(last)))
}

/// A Feature without the version a row writes after its stable id (`--transcript_version`):
/// the text before a last dot that only digits follow.
fn unversioned(feature: &str) -> &str {
    match feature.rsplit_once('.') {
        Some((stable_id, version))
            if !version.is_empty() && version.bytes().all(|byte| byte.is_ascii_digit()) =>
        {
            stable_id
        }
        _ => feature,
    }
}

/// The coding ends of the transcripts of a candidate JSON cache, read a sequence at a time
/// when a row first asks for a transcript on it.
#[derive(Debug)]
pub struct CacheTranscripts {
    cache: PathBuf,
    threads: usize,
    read: Mutex<HashMap<String, Arc<Records>>>,
    /// The names the cache's `chr_synonyms.txt` pairs each name with, read once.
    synonyms: OnceLock<HashMap<String, Vec<String>>>,
}

/// The names a synonyms file pairs each name with, as vep-rs reads the file
/// (`crates/vep-core/src/synonyms.rs`): the first word of a line and each later word are
/// synonyms of each other, and two later words of one line are not.
fn synonyms(text: &str) -> HashMap<String, Vec<String>> {
    let mut paired: HashMap<String, Vec<String>> = HashMap::new();
    for line in text.lines() {
        let mut words = line.split_whitespace();
        let Some(reference) = words.next() else {
            continue;
        };
        for synonym in words {
            for (name, other) in [(reference, synonym), (synonym, reference)] {
                let others = paired.entry(name.to_string()).or_default();
                if !others.iter().any(|held| held == other) {
                    others.push(other.to_string());
                }
            }
        }
    }
    paired
}

/// The spellings under which a sequence an output names may be filed: the name itself,
/// the name without a leading `chr` or with one, and for the mitochondrial names `MT` and
/// `M` each other (vep-rs writes a record's `M` as `MT`).
fn spellings(sequence: &str) -> Vec<String> {
    let mut names = vec![sequence.to_string()];
    let bare = match sequence.get(..3) {
        Some(prefix) if prefix.eq_ignore_ascii_case("chr") => &sequence[3..],
        _ => sequence,
    };
    for mitochondrial in [("MT", "M"), ("M", "MT")] {
        if bare == mitochondrial.0 {
            names.push(mitochondrial.1.to_string());
        }
    }
    for name in names.clone() {
        let other = match name.get(..3) {
            Some(prefix) if prefix.eq_ignore_ascii_case("chr") => name[3..].to_string(),
            _ => format!("chr{name}"),
        };
        names.push(other);
    }
    names.dedup();
    names
}

impl CacheTranscripts {
    /// The transcripts of the cache directory `cache`, its shards read on `threads` threads
    /// at most. Nothing is read here.
    pub fn new(cache: &Path, threads: usize) -> CacheTranscripts {
        CacheTranscripts {
            cache: cache.to_path_buf(),
            threads: threads.max(1),
            read: Mutex::new(HashMap::new()),
            synonyms: OnceLock::new(),
        }
    }

    /// The directory under `transcripts/` that files the sequence an output names
    /// `sequence`: the directory of that name, or of one of its spellings (`spellings`), or
    /// of the first name the cache's `chr_synonyms.txt` pairs with a spelling that the cache
    /// holds, as vep-rs finds a record's sequence in a cache
    /// (`crates/vep-core/src/synonyms.rs`, `source_chr_name`). `None` where no directory
    /// answers, or two spellings lead to two.
    fn directory_of(&self, sequence: &str) -> Option<String> {
        let root = self.cache.join("transcripts");
        let held = |name: &str| {
            let named =
                !name.is_empty() && name != "." && name != ".." && !name.contains(['/', '\\']);
            named && root.join(name).is_dir()
        };
        if held(sequence) {
            return Some(sequence.to_string());
        }
        let synonyms = self.synonyms.get_or_init(|| {
            std::fs::read_to_string(self.cache.join("chr_synonyms.txt"))
                .map(|text| synonyms(&text))
                .unwrap_or_default()
        });
        let mut found: Vec<String> = Vec::new();
        for spelling in spellings(sequence) {
            let paired = synonyms.get(&spelling).map_or(&[][..], Vec::as_slice);
            let directory = match held(&spelling) {
                true => Some(spelling.as_str()),
                false => paired.iter().map(String::as_str).find(|name| held(name)),
            };
            if let Some(directory) = directory.filter(|name| !found.iter().any(|f| f == name)) {
                found.push(directory.to_string());
            }
        }
        match found.len() {
            1 => found.pop(),
            _ => None,
        }
    }

    /// The records of the sequence the directory `sequence` files, by stable id. The first
    /// failure in the order of the shards ends the reading.
    fn records(&self, sequence: &str) -> io::Result<Arc<Records>> {
        let mut read = self
            .read
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if let Some(records) = read.get(sequence) {
            return Ok(Arc::clone(records));
        }
        let root = self.cache.join("transcripts");
        let shards = shards_of(&root.join(sequence));
        let next = AtomicUsize::new(0);
        let read_some = || -> Vec<ShardRead> {
            let mut done = Vec::new();
            loop {
                let index = next.fetch_add(1, Ordering::Relaxed);
                let Some(shard) = shards.get(index) else {
                    return done;
                };
                let records = read_shard(shard)
                    .map(|records| records.iter().filter_map(coding_ends_of).collect());
                done.push((index, records));
            }
        };
        let workers = self.threads.min(shards.len()).max(1);
        let mut by_shard: Vec<ShardRead> = std::thread::scope(|scope| {
            let spawned: Vec<_> = (1..workers).map(|_| scope.spawn(read_some)).collect();
            let mut done = read_some();
            for worker in spawned {
                match worker.join() {
                    Ok(more) => done.extend(more),
                    Err(panic) => std::panic::resume_unwind(panic),
                }
            }
            done
        });
        by_shard.sort_unstable_by_key(|(index, _)| *index);
        let mut records = Records::new();
        for (_, of_shard) in by_shard {
            for (stable_id, ends) in of_shard? {
                let readings = records.entry(stable_id).or_default();
                if !readings.iter().any(|held| **held == ends) {
                    readings.push(Arc::new(ends));
                }
            }
        }
        let records = Arc::new(records);
        read.insert(sequence.to_string(), Arc::clone(&records));
        Ok(records)
    }

    /// The coding ends of the transcript `feature` of a row at `location`: of the record of
    /// that stable id, on the sequence the Location names, whose span shares a base with the
    /// Location, an insertion's being the two bases it stands between. The sequence is the
    /// directory of the Location's own name, or the one the cache's synonyms give it
    /// (`directory_of`). A Feature carrying a version the cache's stable id lacks is looked
    /// up without it. `None` for a Location that names no position or no directory of the
    /// cache, for a transcript the cache does not hold, and where no record of the stable id
    /// holds the Location or more than one does: one id stands on several records of a
    /// sequence, and a row names one of them by where it lies. A cache without a
    /// `transcripts/` directory is a failure.
    pub fn coding_ends(
        &self,
        location: &str,
        feature: &str,
    ) -> io::Result<Option<Arc<CodingEnds>>> {
        if !self.cache.join("transcripts").is_dir() {
            return Err(io::Error::other(format!(
                "--json-cache: {} has no transcripts/ directory",
                self.cache.display()
            )));
        }
        let Some((first, last)) = location_span(location) else {
            return Ok(None);
        };
        let named = location
            .rsplit_once(':')
            .map_or("", |(sequence, _)| sequence);
        let Some(sequence) = self.directory_of(named) else {
            return Ok(None);
        };
        let records = self.records(&sequence)?;
        let Some(readings) = records
            .get(feature)
            .or_else(|| records.get(unversioned(feature)))
        else {
            return Ok(None);
        };
        let mut holding = readings.iter().filter(|ends| {
            ends.span
                .is_some_and(|(start, end)| first <= end && start <= last)
        });
        Ok(match (holding.next(), holding.next()) {
            (Some(ends), None) => Some(Arc::clone(ends)),
            _ => None,
        })
    }
}

#[cfg(test)]
#[path = "cache_transcripts_tests.rs"]
mod tests;
