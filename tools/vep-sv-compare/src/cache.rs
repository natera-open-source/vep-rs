// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The vep-rs JSON cache as the authority for which chromosome a transcript is on, where
//! it lies and how it is built.
//!
//! A cache files its transcripts under `transcripts/<chromosome>/`, one shard a file named
//! `*.json` or, gzipped, `*.json.gz`. A shard is an array of transcript records or an
//! object of such arrays (`{"21": [record, null, ...]}`), read in file order; where an
//! object repeats a key, the last array stands at the first one's place. The directories
//! are read in name order and the shards of each in name order, and a later record of one
//! chromosome and id takes the place of an earlier one, whatever the number of threads.
//!
//! A shard is read as Python's `json` module reads a text: beside JSON it takes the numbers
//! `NaN`, `Infinity` and `-Infinity`, an integer of more than 4,300 digits is no number,
//! and the last of a repeated key stands. A shard that cannot be opened, is not UTF-8 or
//! is not such a text holds nothing, as does a gzipped shard that is cut short or fails
//! its checksum.
//!
//! A record is read as Python reads the object: a value that it takes as false stands for
//! an absent one, and a shape on which its reading ends in an exception ends the loading in
//! a [`Failure`] with that exception's text. Those shapes are an id that is an array or an
//! object, an infinite number where an integer is read, a number where a list is read, a
//! value that is no object where one is asked for a key (an element of the exons, the pairs
//! or the attributes, the feature cache, its mapper), and a mature miRNA range of more than
//! 4,300 digits. The first such record in the order of the shards and of their records
//! decides; within a record the id comes first, then start, end and strand, and, for a
//! record that gives those three as integers, the exons, the mapper (its pairs under two
//! coding bounds above zero) and the attributes. A gzipped shard whose deflate stream is
//! corrupt and a `transcripts` directory that cannot be listed end the loading as well.
//!
//! A record whose id is a number other than zero, or `true`, makes its chromosome known
//! and enters no view, since no row of an output can name it. The transcript total counts
//! such ids as Python's set of a chromosome's ids does: numbers that are equal are one id.
//!
//! An integer is a JSON number or a text that Python's `int` reads: a fraction is cut
//! towards zero, `true` is 1, and a text may carry blanks around it, a sign, single
//! underscores between digits and the digits of any script. It is read up to 128 bits; a
//! longer one is no integer here.

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fs;
use std::io;
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};

use flate2::{Crc, Decompress, FlushDecompress, Status};
use memchr::memchr2;

use crate::failure::Failure;
use crate::mate_context::{
    decimal_digit, mature_range_failure, model_from_cache_record, CacheRecord, MapperPair,
    TranscriptModel,
};
use crate::rules::py_int;
use crate::text::{os_failure, py_join, py_path, INT_MAX_STR_DIGITS};
use crate::threads::run_in_order;

/// Chromosome to the transcript ids the cache holds for it.
pub type TranscriptsByChromosome = HashMap<String, HashSet<String>>;

/// Chromosome, then transcript id, to the transcript's `(start, end)`.
pub type TranscriptSpans = HashMap<String, HashMap<String, (i128, i128)>>;

/// Chromosome, then transcript id, to the transcript's structure.
pub type TranscriptModels = HashMap<String, HashMap<String, TranscriptModel>>;

/// Three views of one cache, each keyed by chromosome first: the pseudoautosomal
/// transcripts are filed under both X and Y with their own coordinates on each.
///
/// A view holds a chromosome as soon as one record of that chromosome has something to
/// file in it, so a view is empty exactly when the cache gave it nothing: no record with an
/// id, no record with a span, no record with a structure.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct CacheAuthority {
    /// Chromosome to the transcript ids the cache holds for it.
    pub by_chr: TranscriptsByChromosome,
    /// Chromosome, then transcript id, to the transcript's `(start, end)`: the records that
    /// give both as integers.
    pub spans: TranscriptSpans,
    /// Chromosome, then transcript id, to the transcript's structure: the records that give
    /// start, end and strand as integers.
    pub models: TranscriptModels,
    /// Chromosome to the ids of its records that are a number or `true`, each number once.
    /// No row of an output can name such an id and no view above holds one; they are
    /// counted in [`CacheAuthority::transcript_total`].
    pub number_ids: HashMap<String, HashSet<String>>,
}

impl CacheAuthority {
    /// The transcript ids filed under a chromosome, or `None` when the cache files none.
    pub fn transcripts(&self, chrom: &str) -> Option<&HashSet<String>> {
        self.by_chr.get(chrom)
    }

    pub fn span(&self, chrom: &str, transcript: &str) -> Option<(i128, i128)> {
        self.spans.get(chrom)?.get(transcript).copied()
    }

    pub fn model(&self, chrom: &str, transcript: &str) -> Option<&TranscriptModel> {
        self.models.get(chrom)?.get(transcript)
    }

    /// How many transcripts the cache holds, as the first line of stdout states it: the
    /// distinct ids of each chromosome, a number or `true` among them, as the set of ids
    /// Python files under a chromosome counts them.
    pub fn transcript_total(&self) -> u64 {
        let texts: usize = self.by_chr.values().map(HashSet::len).sum();
        let numbers: usize = self.number_ids.values().map(HashSet::len).sum();
        (texts + numbers) as u64
    }

    /// Files one record of `chrom`. A record whose id is no text makes its chromosome known
    /// to each view it would be filed in and adds no id to a view.
    fn file(&mut self, chrom: &str, record: Filed) {
        let Filed {
            id,
            number_id,
            span,
            model,
        } = record;
        if let Some(number) = number_id {
            filed_under(&mut self.number_ids, chrom).insert(number);
        }
        let ids = filed_under(&mut self.by_chr, chrom);
        if let Some(id) = &id {
            if !ids.contains(id) {
                ids.insert(id.clone());
            }
        }
        if let Some(span) = span {
            let spans = filed_under(&mut self.spans, chrom);
            if let Some(id) = &id {
                spans.insert(id.clone(), span);
            }
        }
        if let Some(model) = model {
            let models = filed_under(&mut self.models, chrom);
            if let Some(id) = id {
                models.insert(id, model);
            }
        }
    }
}

fn filed_under<'a, V: Default>(view: &'a mut HashMap<String, V>, chrom: &str) -> &'a mut V {
    if !view.contains_key(chrom) {
        view.insert(chrom.to_string(), V::default());
    }
    view.get_mut(chrom).expect("the chromosome was just filed")
}

/// What one record files: its id when that is a text, the text that stands for an id that
/// is a number, its span and its structure.
struct Filed {
    id: Option<String>,
    number_id: Option<String>,
    span: Option<(i128, i128)>,
    model: Option<TranscriptModel>,
}

/// The shards that may be read past the last one filed: the bound on what is held between
/// the reading and the filing.
const SHARDS_A_ROUND: usize = 256;

/// Reads the cache under `dir` once, its shards across `threads` workers, and returns the
/// three views of it; all three are empty when the directory holds no `transcripts`
/// directory. The records of a shard are filed, in the order of the shards, while later
/// shards are read. A record or a shard on which Python's reading ends in an exception is a
/// failure with that exception's text, the first in the order of the shards deciding.
pub fn load_cache_authority(dir: &Path, threads: usize) -> Result<CacheAuthority, Failure> {
    let mut authority = CacheAuthority::default();
    let shards = cache_shards(dir)?;
    let read = |index: usize| -> Result<Vec<Filed>, Failure> {
        let (chrom, path) = &shards[index];
        let filed = shard_transcripts(path)?
            .into_iter()
            .map(|record| Filed {
                span: record.start.zip(record.end),
                model: model_from_cache_record(chrom, &record),
                id: record.stable_id,
                number_id: record.number_id,
            })
            .collect();
        Ok(filed)
    };
    let mut failure = None;
    let file = |index: usize, records: Result<Vec<Filed>, Failure>| match records {
        Ok(records) => {
            for record in records {
                authority.file(&shards[index].0, record);
            }
            ControlFlow::Continue(())
        }
        Err(first) => {
            failure = Some(first);
            ControlFlow::Break(())
        }
    };
    run_in_order(threads, shards.len(), SHARDS_A_ROUND, read, file);
    match failure {
        Some(failure) => Err(failure),
        None => Ok(authority),
    }
}

/// The chromosome and span views of the cache; see [`load_cache_authority`].
pub fn load_transcript_index(
    dir: &Path,
    threads: usize,
) -> Result<(TranscriptsByChromosome, TranscriptSpans), Failure> {
    let authority = load_cache_authority(dir, threads)?;
    Ok((authority.by_chr, authority.spans))
}

/// Each chromosome's transcript ids; see [`load_cache_authority`].
pub fn load_transcripts_by_chromosome(
    dir: &Path,
    threads: usize,
) -> Result<TranscriptsByChromosome, Failure> {
    Ok(load_transcript_index(dir, threads)?.0)
}

/// The entries of a directory in the order of their names, as bytes.
fn entries_by_name(dir: &Path) -> io::Result<Vec<PathBuf>> {
    let mut paths: Vec<PathBuf> = fs::read_dir(dir)?
        .flatten()
        .map(|entry| entry.path())
        .collect();
    paths.sort_unstable_by(|a, b| a.file_name().cmp(&b.file_name()));
    Ok(paths)
}

/// The shards of the cache under `dir`, each with the chromosome its directory names: the
/// directories under `transcripts` in name order and, in each, the entries named `*.json`
/// or `*.json.gz` in name order. A `transcripts` directory that cannot be listed is a
/// failure; a chromosome's directory that cannot be listed holds no shard.
pub fn cache_shards(dir: &Path) -> Result<Vec<(String, PathBuf)>, Failure> {
    let root = dir.join("transcripts");
    if !root.is_dir() {
        return Ok(Vec::new());
    }
    let chromosomes = entries_by_name(&root).map_err(|error| {
        let spelled = py_join(&py_path(dir.as_os_str()), OsStr::new("transcripts"));
        os_failure(&error, Some(&spelled))
    })?;
    let mut shards = Vec::new();
    for chr_dir in chromosomes {
        if !chr_dir.is_dir() {
            continue;
        }
        let Some(chrom) = chr_dir.file_name().map(|name| name.to_string_lossy()) else {
            continue;
        };
        for shard in entries_by_name(&chr_dir).unwrap_or_default() {
            let name = shard.file_name().map(|name| name.to_string_lossy());
            if name.is_some_and(|name| name.ends_with(".json") || name.ends_with(".json.gz")) {
                shards.push((chrom.to_string(), shard));
            }
        }
    }
    Ok(shards)
}

/// The transcript records of a shard file, in file order: those with a `stable_id` that
/// Python reads as true. None for a file that cannot be read whole; a failure where
/// Python's reading of the file ends in an exception.
pub fn shard_transcripts(shard: &Path) -> Result<Vec<CacheRecord>, Failure> {
    let Ok(mut bytes) = fs::read(shard) else {
        return Ok(Vec::new());
    };
    if shard.extension().is_some_and(|extension| extension == "gz") {
        match gunzip(&bytes)? {
            Some(inflated) => bytes = inflated,
            None => return Ok(Vec::new()),
        }
    }
    match std::str::from_utf8(&bytes) {
        Ok(text) => Ok(parse_shard(text)?.unwrap_or_default()),
        Err(_) => Ok(Vec::new()),
    }
}

/// The length of the gzip member header at the start of `bytes`, where the member's
/// deflate stream begins; `None` for bytes that begin no member of the deflate method.
fn gzip_header(bytes: &[u8]) -> Option<usize> {
    let [0x1f, 0x8b, 8, flags, ..] = bytes else {
        return None;
    };
    let mut at = 10;
    if flags & 4 != 0 {
        let length = bytes.get(at..at + 2)?;
        at += 2 + usize::from(u16::from_le_bytes([length[0], length[1]]));
    }
    // A name and a comment each run to a zero byte, or to the end of the file.
    for text in [8, 16] {
        if flags & text != 0 {
            let rest = bytes.get(at..)?;
            at += rest
                .iter()
                .position(|byte| *byte == 0)
                .map_or(rest.len(), |zero| zero + 1);
        }
    }
    if flags & 2 != 0 {
        at += 2;
    }
    (at <= bytes.len()).then_some(at)
}

/// The bytes a gzipped shard holds, as Python's `gzip` module reads the file: every member
/// in turn, zero bytes after a member passed over. `None` for a file that module refuses or
/// finds cut short (no member at its start, a checksum or a length that does not match,
/// other bytes after a member); a deflate stream that is itself corrupt is a failure.
fn gunzip(bytes: &[u8]) -> Result<Option<Vec<u8>>, Failure> {
    let mut text = Vec::new();
    let mut rest = bytes;
    while !rest.is_empty() {
        let Some(deflate) = gzip_header(rest) else {
            return Ok(None);
        };
        rest = &rest[deflate..];
        let member = text.len();
        let mut inflater = Decompress::new(false);
        loop {
            text.reserve((2 * rest.len()).max(1 << 15));
            let (read, written) = (inflater.total_in(), inflater.total_out());
            let status = inflater
                .decompress_vec(rest, &mut text, FlushDecompress::None)
                .map_err(|_| Failure::value("Error -3 while decompressing data"))?;
            let consumed = (inflater.total_in() - read) as usize;
            rest = &rest[consumed..];
            if status == Status::StreamEnd {
                break;
            }
            if consumed == 0 && inflater.total_out() == written {
                return Ok(None);
            }
        }
        let Some((trailer, after)) = rest.split_at_checked(8) else {
            return Ok(None);
        };
        let mut crc = Crc::new();
        crc.update(&text[member..]);
        let length = (text.len() - member) as u32;
        if trailer[..4] != crc.sum().to_le_bytes() || trailer[4..] != length.to_le_bytes() {
            return Ok(None);
        }
        rest = after;
        while let [0, after @ ..] = rest {
            rest = after;
        }
    }
    Ok(Some(text))
}

/// The transcript records of a shard's text, in file order; `None` for a text Python's
/// `json` module does not read; a failure where Python's reading of a record ends in an
/// exception.
pub fn parse_shard(text: &str) -> Result<Option<Vec<CacheRecord>>, Failure> {
    match Json::new(text).shard() {
        Ok(read) => read.map(Some),
        Err(NotJson) => Ok(None),
    }
}

/// `AttributeError`: a value that is no object was asked for a key.
fn no_get(kind: &str) -> Failure {
    Failure::value(format!("'{kind}' object has no attribute 'get'"))
}

/// `TypeError`: a number or `true` where a list is run through.
fn not_iterable(kind: &str) -> Failure {
    Failure::value(format!("'{kind}' object is not iterable"))
}

/// `TypeError`: an id that is an array or an object.
fn unhashable(kind: &str) -> Failure {
    Failure::value(format!(
        "cannot use '{kind}' as a set element (unhashable type: '{kind}')"
    ))
}

/// `OverflowError`: an infinite number where an integer is read.
fn infinite() -> Failure {
    Failure::value("cannot convert float infinity to integer")
}

/// What reading one field came to: its value, or the exception Python's reading of it ends
/// in.
type Read<T> = Result<T, Failure>;

/// A reading that counts only where Python reaches the field: elsewhere an exception in it
/// is never met, and the field is as absent.
fn reached<T: Default>(read: Read<T>, reached: bool) -> Read<T> {
    if reached {
        read
    } else {
        Ok(read.unwrap_or_default())
    }
}

/// 2^127, the first double past every integer of the word.
const TWO_POW_127: f64 = 170_141_183_460_469_231_731_687_303_715_884_105_728.0;

/// What Python's `int` makes of a value.
#[derive(Clone, Copy)]
enum Integer {
    None,
    Some(i128),
    /// An infinite number: `int` raises on it.
    Infinite,
}

impl Integer {
    fn value(self) -> Option<i128> {
        match self {
            Integer::Some(value) => Some(value),
            Integer::None | Integer::Infinite => None,
        }
    }

    fn finite(self) -> Read<Option<i128>> {
        match self {
            Integer::Infinite => Err(infinite()),
            Integer::None | Integer::Some(_) => Ok(self.value()),
        }
    }
}

/// The integers of an object when it gives every one of them; `None` when one is missing;
/// the failure of an infinite one.
fn integers<const N: usize>(numbers: [Integer; N]) -> Read<Option<[i128; N]>> {
    let mut values = [0; N];
    let mut whole = true;
    for (value, number) in values.iter_mut().zip(numbers) {
        match number.finite()? {
            Some(integer) => *value = integer,
            None => whole = false,
        }
    }
    Ok(whole.then_some(values))
}

/// `int()` of a text: [`py_int`] with the digits of any script.
fn text_integer(text: &str) -> Option<i128> {
    if text.is_ascii() {
        return py_int(text);
    }
    let ascii: String = text
        .chars()
        .map(|c| {
            decimal_digit(c)
                .and_then(|digit| char::from_digit(digit, 10))
                .unwrap_or(c)
        })
        .collect();
    py_int(&ascii)
}

/// A record's `stable_id` as Python's truth test and the views read it.
enum StableId {
    /// No id: null, false, zero, an empty text, an empty array or object.
    Absent,
    Text(String),
    /// A number other than zero, or `true`: an id no row of an output can name, as the
    /// text of [`number_id`].
    Number(String),
    /// An array or an object that holds something, by the name of its Python type.
    Unhashable(&'static str),
}

/// What a record's `variation_effect_feature_cache.mapper` gives.
struct Mapper {
    coding_start: Integer,
    coding_end: Integer,
    pairs: Read<Vec<MapperPair>>,
}

impl Default for Mapper {
    fn default() -> Mapper {
        Mapper {
            coding_start: Integer::None,
            coding_end: Integer::None,
            pairs: Ok(Vec::new()),
        }
    }
}

/// A record as its text gives it, before the order in which Python meets its parts decides
/// what it comes to.
struct RecordReading {
    id: StableId,
    start: Integer,
    end: Integer,
    strand: Integer,
    biotype: String,
    exons: Read<Vec<(i128, i128)>>,
    mapper: Read<Mapper>,
    mirna: Read<Vec<String>>,
}

impl RecordReading {
    /// The record, `None` for one without an id, or the exception Python's reading of it
    /// ends in. Python files the id, converts start, end and strand, and only for a record
    /// that gives all three as integers reads the exons, then the mapper, whose pairs it
    /// reads under two coding bounds above zero, then the attributes.
    fn settle(self) -> Read<Option<CacheRecord>> {
        let (stable_id, number_id) = match self.id {
            StableId::Absent => return Ok(None),
            StableId::Unhashable(kind) => return Err(unhashable(kind)),
            StableId::Text(id) => (Some(id), None),
            StableId::Number(id) => (None, Some(id)),
        };
        let (start, end, strand) = (
            self.start.finite()?,
            self.end.finite()?,
            self.strand.finite()?,
        );
        let modelled = start.is_some() && end.is_some() && strand.is_some();
        let exons = reached(self.exons, modelled)?;
        let mapper = reached(self.mapper, modelled)?;
        let (cdna_coding_start, cdna_coding_end) = if modelled {
            (mapper.coding_start.finite()?, mapper.coding_end.finite()?)
        } else {
            (mapper.coding_start.value(), mapper.coding_end.value())
        };
        let coding = |bound: Option<i128>| bound.is_some_and(|bound| bound > 0);
        let paired = modelled && coding(cdna_coding_start) && coding(cdna_coding_end);
        let pairs = reached(mapper.pairs, paired)?;
        let mirna = reached(self.mirna, modelled)?;
        Ok(Some(CacheRecord {
            stable_id,
            number_id,
            start,
            end,
            strand,
            biotype: self.biotype,
            exons,
            cdna_coding_start,
            cdna_coding_end,
            pairs,
            mirna,
        }))
    }
}

/// The text that stands for an id that is a number or `true`: two ids have one text
/// exactly when Python holds them as one element of a set. `true` is 1, a number without
/// a fraction is its integer whatever its spelling (`7.0` and `7`, `1e23` and
/// `99999999999999991611392`), and every `NaN` is one id. `literal` is the id as the
/// shard writes it.
fn number_id(scalar: &Scalar, literal: &str) -> String {
    match scalar {
        Scalar::Float(value) if value.is_nan() => "nan".to_string(),
        Scalar::Float(value) if value.is_infinite() => {
            if *value > 0.0 { "inf" } else { "-inf" }.to_string()
        }
        // Every digit of a double without a fraction is written, so the text is that of
        // the integer equal to it.
        Scalar::Float(value) if value.fract() == 0.0 => format!("{value:.0}"),
        // The shortest digits that give the double back, with an exponent no integer has.
        Scalar::Float(value) => format!("{value:e}"),
        Scalar::Bool(_) => "1".to_string(),
        // An integer is written without a leading zero and without a plus sign.
        Scalar::Int(_) | Scalar::LongInt | Scalar::Null | Scalar::Text(_) => literal.to_string(),
    }
}

/// A text that Python's `json` module does not read.
struct NotJson;

type Parsed<T> = Result<T, NotJson>;

/// A value that is neither an array nor an object, as Python holds it.
enum Scalar<'a> {
    Null,
    Bool(bool),
    Int(i128),
    /// An integer of more than 128 bits.
    LongInt,
    Float(f64),
    Text(Cow<'a, str>),
}

impl Scalar<'_> {
    /// Whether Python reads the value as true.
    fn is_true(&self) -> bool {
        match self {
            Scalar::Null => false,
            Scalar::Bool(value) => *value,
            Scalar::Int(value) => *value != 0,
            Scalar::LongInt => true,
            Scalar::Float(value) => *value != 0.0,
            Scalar::Text(text) => !text.is_empty(),
        }
    }

    /// The name of the value's Python type, as an exception's text gives it.
    fn type_name(&self) -> &'static str {
        match self {
            Scalar::Null => "NoneType",
            Scalar::Bool(_) => "bool",
            Scalar::Int(_) | Scalar::LongInt => "int",
            Scalar::Float(_) => "float",
            Scalar::Text(_) => "str",
        }
    }
}

/// Four hexadecimal digits at `at`, as their number.
fn hex4(bytes: &[u8], at: usize) -> Parsed<u32> {
    let digits = bytes.get(at..at + 4).ok_or(NotJson)?;
    digits
        .iter()
        .try_fold(0, |value, digit| {
            Some(value * 16 + char::from(*digit).to_digit(16)?)
        })
        .ok_or(NotJson)
}

/// The reader of one shard's text: the JSON that Python's `json` module reads.
struct Json<'a> {
    text: &'a str,
    at: usize,
    /// The arrays and objects open while a value is passed over.
    open: Vec<u8>,
}

/// The text itself: blanks, literals, numbers, strings, and passing over a value.
impl<'a> Json<'a> {
    fn new(text: &'a str) -> Json<'a> {
        Json {
            text,
            at: 0,
            open: Vec::new(),
        }
    }

    fn bytes(&self) -> &'a [u8] {
        self.text.as_bytes()
    }

    fn peek(&self) -> Option<u8> {
        self.bytes().get(self.at).copied()
    }

    fn blank(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.at += 1;
        }
    }

    fn eat(&mut self, byte: u8) -> bool {
        let found = self.peek() == Some(byte);
        self.at += usize::from(found);
        found
    }

    fn expect(&mut self, byte: u8) -> Parsed<()> {
        if self.eat(byte) {
            Ok(())
        } else {
            Err(NotJson)
        }
    }

    fn word(&mut self, word: &str) -> Parsed<()> {
        if self.bytes()[self.at..].starts_with(word.as_bytes()) {
            self.at += word.len();
            Ok(())
        } else {
            Err(NotJson)
        }
    }

    fn digits(&mut self) {
        while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
            self.at += 1;
        }
    }

    /// Passes over a number and says whether it has a fraction or an exponent. A dot that
    /// no digit follows and an exponent without a digit are left where they are, no part
    /// of the number, and an integer of more than 4,300 digits is none.
    fn number(&mut self) -> Parsed<bool> {
        let bytes = self.bytes();
        self.eat(b'-');
        let integer = self.at;
        match self.peek() {
            Some(b'0') => self.at += 1,
            Some(b'1'..=b'9') => self.digits(),
            _ => return Err(NotJson),
        }
        let integer_digits = self.at - integer;
        let mut float = false;
        if self.peek() == Some(b'.') && bytes.get(self.at + 1).is_some_and(u8::is_ascii_digit) {
            self.at += 1;
            self.digits();
            float = true;
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            let mut exponent = self.at + 1;
            if matches!(bytes.get(exponent), Some(b'+' | b'-')) {
                exponent += 1;
            }
            if bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
                self.at = exponent;
                self.digits();
                float = true;
            }
        }
        if !float && integer_digits > INT_MAX_STR_DIGITS {
            return Err(NotJson);
        }
        Ok(float)
    }

    /// Passes over a string, the reader at its opening quote, and says whether it holds an
    /// escape.
    fn pass_string(&mut self) -> Parsed<bool> {
        let bytes = self.bytes();
        let mut at = self.at + 1;
        let mut escaped = false;
        loop {
            let run = memchr2(b'"', b'\\', &bytes[at..]).ok_or(NotJson)?;
            if bytes[at..at + run].iter().any(|byte| *byte < 0x20) {
                return Err(NotJson);
            }
            at += run;
            if bytes[at] == b'"' {
                self.at = at + 1;
                return Ok(escaped);
            }
            escaped = true;
            at += match bytes.get(at + 1) {
                Some(b'"' | b'\\' | b'/' | b'b' | b'f' | b'n' | b'r' | b't') => 2,
                Some(b'u') => {
                    hex4(bytes, at + 2)?;
                    6
                }
                _ => return Err(NotJson),
            };
        }
    }

    /// Reads a string, the reader at its opening quote. A surrogate that is no half of a
    /// pair, which Python keeps as it is, is read as U+FFFD.
    fn string(&mut self) -> Parsed<Cow<'a, str>> {
        let start = self.at + 1;
        let escaped = self.pass_string()?;
        let mut rest = &self.text[start..self.at - 1];
        if !escaped {
            return Ok(Cow::Borrowed(rest));
        }
        let mut text = String::with_capacity(rest.len());
        while let Some(at) = rest.find('\\') {
            text.push_str(&rest[..at]);
            let bytes = rest.as_bytes();
            let mut used = 2;
            text.push(match bytes[at + 1] {
                b'b' => '\u{8}',
                b'f' => '\u{c}',
                b'n' => '\n',
                b'r' => '\r',
                b't' => '\t',
                b'u' => {
                    let unit = hex4(bytes, at + 2)?;
                    used = 6;
                    let low = ((0xd800..0xdc00).contains(&unit)
                        && rest[at + 6..].starts_with("\\u"))
                    .then(|| hex4(bytes, at + 8).ok())
                    .flatten()
                    .filter(|low| (0xdc00..0xe000).contains(low));
                    let point = match low {
                        Some(low) => {
                            used = 12;
                            0x10000 + ((unit - 0xd800) << 10) + (low - 0xdc00)
                        }
                        None => unit,
                    };
                    char::from_u32(point).unwrap_or('\u{fffd}')
                }
                quoted => char::from(quoted),
            });
            rest = &rest[at + used..];
        }
        text.push_str(rest);
        Ok(Cow::Owned(text))
    }

    /// Passes over a value that is neither an array nor an object.
    fn pass_scalar(&mut self) -> Parsed<()> {
        match self.peek() {
            Some(b'"') => self.pass_string().map(|_| ()),
            Some(b'n') => self.word("null"),
            Some(b't') => self.word("true"),
            Some(b'f') => self.word("false"),
            Some(b'N') => self.word("NaN"),
            Some(b'I') => self.word("Infinity"),
            Some(b'-') if self.bytes()[self.at..].starts_with(b"-Infinity") => {
                self.word("-Infinity")
            }
            _ => self.number().map(|_| ()),
        }
    }

    /// Reads a value that is neither an array nor an object.
    fn scalar(&mut self) -> Parsed<Scalar<'a>> {
        let start = self.at;
        match self.peek() {
            Some(b'"') => return self.string().map(Scalar::Text),
            Some(b'-' | b'0'..=b'9') if !self.bytes()[start..].starts_with(b"-Infinity") => {
                let float = self.number()?;
                let literal = &self.text[start..self.at];
                return Ok(if float {
                    Scalar::Float(literal.parse().map_err(|_| NotJson)?)
                } else {
                    literal.parse().map_or(Scalar::LongInt, Scalar::Int)
                });
            }
            _ => self.pass_scalar()?,
        }
        Ok(match self.bytes()[start] {
            b'n' => Scalar::Null,
            b't' => Scalar::Bool(true),
            b'f' => Scalar::Bool(false),
            b'N' => Scalar::Float(f64::NAN),
            b'I' => Scalar::Float(f64::INFINITY),
            _ => Scalar::Float(f64::NEG_INFINITY),
        })
    }

    /// Passes over a key and its colon, to the start of the key's value.
    fn pass_key(&mut self) -> Parsed<()> {
        if self.peek() != Some(b'"') {
            return Err(NotJson);
        }
        self.pass_string()?;
        self.blank();
        self.expect(b':')?;
        self.blank();
        Ok(())
    }

    /// Passes over one value of any kind and any depth.
    fn pass(&mut self) -> Parsed<()> {
        loop {
            match self.peek() {
                Some(b'[') => {
                    self.at += 1;
                    self.blank();
                    if !self.eat(b']') {
                        self.open.push(b']');
                        continue;
                    }
                }
                Some(b'{') => {
                    self.at += 1;
                    self.blank();
                    if !self.eat(b'}') {
                        self.open.push(b'}');
                        self.pass_key()?;
                        continue;
                    }
                }
                _ => self.pass_scalar()?,
            }
            // A value has ended: it ends the arrays and objects it is the last value of,
            // or another value follows it.
            loop {
                let Some(&close) = self.open.last() else {
                    return Ok(());
                };
                self.blank();
                if self.eat(close) {
                    self.open.pop();
                    continue;
                }
                self.expect(b',')?;
                self.blank();
                if close == b'}' {
                    self.pass_key()?;
                }
                break;
            }
        }
    }

    /// Passes over the array or object at the reader and says whether it holds nothing.
    fn pass_container(&mut self) -> Parsed<bool> {
        let start = self.at;
        self.pass()?;
        Ok(self.bytes()[start + 1..self.at - 1]
            .iter()
            .all(|byte| matches!(byte, b' ' | b'\t' | b'\n' | b'\r')))
    }

    /// Passes over a value of any kind and says whether Python reads it as true, with the
    /// name of its Python type.
    fn truth(&mut self) -> Parsed<(bool, &'static str)> {
        match self.peek() {
            Some(b'[') => Ok((!self.pass_container()?, "list")),
            Some(b'{') => Ok((!self.pass_container()?, "dict")),
            _ => {
                let scalar = self.scalar()?;
                Ok((scalar.is_true(), scalar.type_name()))
            }
        }
    }

    /// Inside an array: moves to its next element, or past its end. `first` is true before
    /// the first element.
    fn element(&mut self, first: &mut bool) -> Parsed<bool> {
        self.blank();
        if self.eat(b']') {
            return Ok(false);
        }
        if !std::mem::take(first) {
            self.expect(b',')?;
            self.blank();
        }
        Ok(true)
    }

    /// Inside an object: reads its next key and moves to that key's value, or past the
    /// object's end. `first` is true before the first key.
    fn key(&mut self, first: &mut bool) -> Parsed<Option<Cow<'a, str>>> {
        self.blank();
        if self.eat(b'}') {
            return Ok(None);
        }
        if !std::mem::take(first) {
            self.expect(b',')?;
            self.blank();
        }
        if self.peek() != Some(b'"') {
            return Err(NotJson);
        }
        let key = self.string()?;
        self.blank();
        self.expect(b':')?;
        self.blank();
        Ok(Some(key))
    }
}

/// What the comparison reads of a shard.
impl<'a> Json<'a> {
    /// A value as the integer Python's `int` makes of it.
    fn integer(&mut self) -> Parsed<Integer> {
        if matches!(self.peek(), Some(b'[' | b'{')) {
            self.pass()?;
            return Ok(Integer::None);
        }
        let value = match self.scalar()? {
            Scalar::Null | Scalar::LongInt => None,
            Scalar::Bool(value) => Some(i128::from(value)),
            Scalar::Int(value) => Some(value),
            Scalar::Float(value) if value.is_infinite() => return Ok(Integer::Infinite),
            Scalar::Float(value) => (value.abs() < TWO_POW_127).then(|| value.trunc() as i128),
            Scalar::Text(text) => text_integer(&text),
        };
        Ok(value.map_or(Integer::None, Integer::Some))
    }

    /// A value as its text; `None` for a value of any other kind.
    fn text(&mut self) -> Parsed<Option<String>> {
        if matches!(self.peek(), Some(b'[' | b'{')) {
            self.pass()?;
            return Ok(None);
        }
        Ok(match self.scalar()? {
            Scalar::Text(text) => Some(text.into_owned()),
            _ => None,
        })
    }

    fn stable_id(&mut self) -> Parsed<StableId> {
        let held = |kind, empty| {
            if empty {
                StableId::Absent
            } else {
                StableId::Unhashable(kind)
            }
        };
        let start = self.at;
        Ok(match self.peek() {
            Some(b'[') => held("list", self.pass_container()?),
            Some(b'{') => held("dict", self.pass_container()?),
            _ => match self.scalar()? {
                Scalar::Text(text) if !text.is_empty() => StableId::Text(text.into_owned()),
                other if other.is_true() => {
                    StableId::Number(number_id(&other, &self.text[start..self.at]))
                }
                _ => StableId::Absent,
            },
        })
    }

    /// Whether the value at the reader is an object to look into. A value Python reads as
    /// false is passed over as absent; any other value that is no object has no key to be
    /// asked for.
    fn object_or_absent(&mut self) -> Parsed<Read<bool>> {
        if self.peek() == Some(b'{') {
            return Ok(Ok(true));
        }
        let (given, kind) = self.truth()?;
        Ok(if given { Err(no_get(kind)) } else { Ok(false) })
    }

    /// The objects of a field Python runs through asking each element for its keys, each
    /// read by `object`, or the exception that ends in: a number or `true` is not iterable,
    /// and an element that is no object has no key to be asked for, as have the characters
    /// of a text and the keys of an object. The field is read to its end either way.
    fn objects<T>(
        &mut self,
        mut object: impl FnMut(&mut Json<'a>) -> Parsed<Read<Option<T>>>,
    ) -> Parsed<Read<Vec<T>>> {
        if self.peek() != Some(b'[') {
            return Ok(match self.truth()? {
                (false, _) => Ok(Vec::new()),
                (true, "str" | "dict") => Err(no_get("str")),
                (true, kind) => Err(not_iterable(kind)),
            });
        }
        self.at += 1;
        let mut objects = Vec::new();
        let mut failure = None;
        let mut first = true;
        while self.element(&mut first)? {
            let element = if self.peek() == Some(b'{') {
                object(self)?
            } else {
                Err(no_get(self.truth()?.1))
            };
            match element {
                Ok(found) if failure.is_none() => objects.extend(found),
                Err(met) if failure.is_none() => failure = Some(met),
                _ => {}
            }
        }
        Ok(failure.map_or(Ok(objects), Err))
    }

    /// An exon: its `start` and `end`, when both are integers.
    fn exon(&mut self) -> Parsed<Read<Option<(i128, i128)>>> {
        let mut bounds = [Integer::None; 2];
        self.at += 1;
        let mut first = true;
        while let Some(key) = self.key(&mut first)? {
            match &*key {
                "start" => bounds[0] = self.integer()?,
                "end" => bounds[1] = self.integer()?,
                _ => self.pass()?,
            }
        }
        Ok(integers(bounds).map(|bounds| bounds.map(|[start, end]| (start, end))))
    }

    /// A pair of a mapper, when it gives all five numbers as integers.
    fn pair(&mut self) -> Parsed<Read<Option<MapperPair>>> {
        const KEYS: [&str; 5] = ["from_start", "from_end", "to_start", "to_end", "ori"];
        let mut numbers = [Integer::None; 5];
        self.at += 1;
        let mut first = true;
        while let Some(key) = self.key(&mut first)? {
            match KEYS.iter().position(|name| *name == key) {
                Some(at) => numbers[at] = self.integer()?,
                None => self.pass()?,
            }
        }
        Ok(integers(numbers).map(|numbers| {
            numbers.map(|[from_start, from_end, to_start, to_end, ori]| MapperPair {
                from_start,
                from_end,
                to_start,
                to_end,
                ori,
            })
        }))
    }

    /// An attribute: its text value when its code is `miRNA`.
    fn attribute(&mut self) -> Parsed<Read<Option<String>>> {
        let (mut code, mut value) = (None, None);
        self.at += 1;
        let mut first = true;
        while let Some(key) = self.key(&mut first)? {
            match &*key {
                "code" => code = self.text()?,
                "value" => value = self.text()?,
                _ => self.pass()?,
            }
        }
        let value = value.filter(|_| code.as_deref() == Some("miRNA"));
        Ok(match value.as_deref().and_then(mature_range_failure) {
            Some(failure) => Err(failure),
            None => Ok(value),
        })
    }

    /// A record's `variation_effect_feature_cache`: its `mapper`.
    fn feature_cache(&mut self) -> Parsed<Read<Mapper>> {
        match self.object_or_absent()? {
            Ok(true) => {}
            Ok(false) => return Ok(Ok(Mapper::default())),
            Err(failure) => return Ok(Err(failure)),
        }
        let mut mapper = Ok(Mapper::default());
        self.at += 1;
        let mut first = true;
        while let Some(key) = self.key(&mut first)? {
            match &*key {
                "mapper" => mapper = self.mapper()?,
                _ => self.pass()?,
            }
        }
        Ok(mapper)
    }

    fn mapper(&mut self) -> Parsed<Read<Mapper>> {
        match self.object_or_absent()? {
            Ok(true) => {}
            Ok(false) => return Ok(Ok(Mapper::default())),
            Err(failure) => return Ok(Err(failure)),
        }
        let mut mapper = Mapper::default();
        self.at += 1;
        let mut first = true;
        while let Some(key) = self.key(&mut first)? {
            match &*key {
                "cdna_coding_start" => mapper.coding_start = self.integer()?,
                "cdna_coding_end" => mapper.coding_end = self.integer()?,
                "pairs" => mapper.pairs = self.objects(Json::pair)?,
                _ => self.pass()?,
            }
        }
        Ok(Ok(mapper))
    }

    /// One object of a shard's array, the reader at its opening brace.
    fn record(&mut self) -> Parsed<RecordReading> {
        let mut record = RecordReading {
            id: StableId::Absent,
            start: Integer::None,
            end: Integer::None,
            strand: Integer::None,
            biotype: String::new(),
            exons: Ok(Vec::new()),
            mapper: Ok(Mapper::default()),
            mirna: Ok(Vec::new()),
        };
        self.at += 1;
        let mut first = true;
        while let Some(key) = self.key(&mut first)? {
            match &*key {
                "stable_id" => record.id = self.stable_id()?,
                "start" => record.start = self.integer()?,
                "end" => record.end = self.integer()?,
                "strand" => record.strand = self.integer()?,
                "biotype" => record.biotype = self.text()?.unwrap_or_default(),
                "exons" => record.exons = self.objects(Json::exon)?,
                "variation_effect_feature_cache" => record.mapper = self.feature_cache()?,
                "attributes" => record.mirna = self.objects(Json::attribute)?,
                _ => self.pass()?,
            }
        }
        Ok(record)
    }

    /// The records of an array, the reader at its opening bracket: its objects with an id,
    /// in order, or the first exception Python's reading of one ends in. The array is read
    /// to its end either way.
    fn records(&mut self) -> Parsed<Read<Vec<CacheRecord>>> {
        self.at += 1;
        let mut records = Vec::new();
        let mut failure = None;
        let mut first = true;
        while self.element(&mut first)? {
            if self.peek() != Some(b'{') {
                self.pass()?;
                continue;
            }
            match self.record()?.settle() {
                Ok(record) if failure.is_none() => records.extend(record),
                Err(met) if failure.is_none() => failure = Some(met),
                _ => {}
            }
        }
        Ok(failure.map_or(Ok(records), Err))
    }

    /// The records of an object of arrays, the reader at its opening brace: the arrays in
    /// the order of their keys, the last value of a repeated key at the first one's place.
    /// Python runs through every value before it reads a record, and a value that is a
    /// number other than zero or `true` is not iterable.
    fn arrays(&mut self) -> Parsed<Read<Vec<CacheRecord>>> {
        let mut values: Vec<Result<Read<Vec<CacheRecord>>, &'static str>> = Vec::new();
        let mut place: HashMap<String, usize> = HashMap::new();
        self.at += 1;
        let mut first = true;
        while let Some(key) = self.key(&mut first)? {
            let value = if self.peek() == Some(b'[') {
                Ok(self.records()?)
            } else {
                match self.truth()? {
                    (true, kind @ ("int" | "float" | "bool")) => Err(kind),
                    _ => Ok(Ok(Vec::new())),
                }
            };
            match place.get(&*key) {
                Some(&at) => values[at] = value,
                None => {
                    place.insert(key.into_owned(), values.len());
                    values.push(value);
                }
            }
        }
        if let Some(kind) = values.iter().find_map(|value| value.as_ref().err()) {
            return Ok(Err(not_iterable(kind)));
        }
        let mut records = Vec::new();
        for read in values.into_iter().flatten() {
            match read {
                Ok(found) => records.extend(found),
                Err(failure) => return Ok(Err(failure)),
            }
        }
        Ok(Ok(records))
    }

    /// A whole shard: an array of records, or an object of such arrays; any other value
    /// holds none.
    fn shard(&mut self) -> Parsed<Read<Vec<CacheRecord>>> {
        self.blank();
        let read = match self.peek() {
            Some(b'[') => self.records()?,
            Some(b'{') => self.arrays()?,
            _ => {
                self.pass_scalar()?;
                Ok(Vec::new())
            }
        };
        self.blank();
        if self.at == self.text.len() {
            Ok(read)
        } else {
            Err(NotJson)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TempDir;

    /// A cache directory holding the given `(chromosome, shard name, bytes)`.
    fn cache_of(shards: &[(&str, &str, &[u8])]) -> TempDir {
        let dir = TempDir::new("vep_sv_compare_cache_").expect("a scratch directory");
        for (chrom, name, bytes) in shards {
            let folder = dir.path.join("transcripts").join(chrom);
            fs::create_dir_all(&folder).expect("a chromosome directory");
            fs::write(folder.join(name), bytes).expect("a shard");
        }
        dir
    }

    fn ids(names: &[&str]) -> HashSet<String> {
        names.iter().map(|name| name.to_string()).collect()
    }

    /// The views of a cache Python reads to its end.
    fn loaded(dir: &TempDir, threads: usize) -> CacheAuthority {
        load_cache_authority(&dir.path, threads).expect("a cache Python reads")
    }

    /// The records of a shard's text Python reads to its end.
    fn records(text: &str) -> Vec<CacheRecord> {
        let shard = parse_shard(text).expect("records Python reads");
        shard.unwrap_or_else(|| panic!("no shard: {text}"))
    }

    /// The one record of a shard written as a JSON object with an id.
    fn record(fields: &str) -> CacheRecord {
        let text = format!(r#"[{{"stable_id": "T", {fields}}}]"#);
        let mut records = records(&text);
        assert_eq!(records.len(), 1, "{text}");
        records.remove(0)
    }

    /// The text of the exception Python's reading of a shard ends in.
    fn failure(text: &str) -> String {
        match parse_shard(text) {
            Err(failure) => failure.message,
            Ok(read) => panic!("{text}: no failure, {read:?}"),
        }
    }

    /// The text of the exception the record of [`modelled`] with these fields ends in.
    fn record_failure(fields: &str) -> String {
        failure(&modelled(fields))
    }

    /// A shard of one record that gives its id, start, end and strand, and `fields`.
    fn modelled(fields: &str) -> String {
        format!(r#"[{{"stable_id": "T", "start": 1, "end": 2, "strand": 1, {fields}}}]"#)
    }

    const NOT_AN_INTEGER: &str = "cannot convert float infinity to integer";

    fn no_get_text(kind: &str) -> String {
        format!("'{kind}' object has no attribute 'get'")
    }

    fn not_iterable_text(kind: &str) -> String {
        format!("'{kind}' object is not iterable")
    }

    /// Caches derived from the reference's own files nest transcripts under a chromosome
    /// key; native ones are a bare list.
    #[test]
    fn test_cache_loader_reads_both_shard_layouts() {
        let dir = cache_of(&[
            (
                "21",
                "0-1000.json",
                br#"[{"stable_id": "ENST_A"}, {"stable_id": "ENST_B"}]"#,
            ),
            (
                "22",
                "0-1000.json",
                br#"{"22": [{"stable_id": "ENST_C"}, null]}"#,
            ),
        ]);
        let got = load_transcripts_by_chromosome(&dir.path, 1);
        let expected: TranscriptsByChromosome = [
            ("21".to_string(), ids(&["ENST_A", "ENST_B"])),
            ("22".to_string(), ids(&["ENST_C"])),
        ]
        .into();
        assert_eq!(got, Ok(expected));
    }

    #[test]
    fn test_cache_loader_returns_empty_for_a_missing_directory() {
        let got = load_transcripts_by_chromosome(Path::new("/nonexistent/cache"), 1);
        assert_eq!(got, Ok(TranscriptsByChromosome::new()));
    }

    /// A pseudoautosomal transcript is filed under X and Y with different coordinates; one
    /// key per chromosome keeps both.
    #[test]
    fn test_models_are_keyed_by_chromosome_and_id() {
        let shard = |start: i128| {
            format!(
                r#"[{{"stable_id": "ENST_PAR", "start": {start}, "end": {end}, "strand": 1,
                "biotype": "protein_coding", "exons": [{{"start": {start}, "end": {end}}}],
                "variation_effect_feature_cache": {{"mapper": {{"pairs": [], "pair_count": 0}}}}}}]"#,
                end = start + 999
            )
        };
        let dir = cache_of(&[
            ("X", "1-1000000.json", shard(200000).as_bytes()),
            ("Y", "1-1000000.json", shard(150000).as_bytes()),
        ]);
        let models = loaded(&dir, 1).models;
        assert_eq!(models["X"]["ENST_PAR"].start, 200000);
        assert_eq!(models["Y"]["ENST_PAR"].start, 150000);
    }

    #[test]
    fn test_a_storable_shaped_shard_is_read() {
        let dir = cache_of(&[(
            "21",
            "1-1000000.json",
            br#"{"21": [null, {"stable_id": "ENST_A", "start": 10, "end": 20, "strand": -1,
            "biotype": "lncRNA", "exons": [{"start": 10, "end": 20}],
            "attributes": [{"code": "miRNA", "value": "3-7"}]}]}"#,
        )]);
        let authority = loaded(&dir, 1);
        let model = authority.model("21", "ENST_A").expect("a model");
        assert_eq!(model.mirna_mature, [(3, 7)]);
        assert_eq!(model.cds, None);
    }

    #[test]
    fn shards_are_read_in_name_order_and_the_last_record_of_an_id_stands() {
        let at = |start: i128| {
            format!(r#"[{{"stable_id": "T", "start": {start}, "end": 9, "strand": 1}}]"#)
        };
        let dir = cache_of(&[
            ("21", "2-3.json", at(2).as_bytes()),
            ("21", "10-11.json", at(10).as_bytes()),
            ("21", "A.json", at(7).as_bytes()),
            ("21", "._2-3.json", b"\x00\x05\x16\x07"),
            ("21", "notes.txt", at(99).as_bytes()),
            ("21", "3-4.json.bak", at(98).as_bytes()),
            ("3", "1-2.json", at(3).as_bytes()),
            ("10", "1-2.json", at(1).as_bytes()),
        ]);
        fs::write(dir.path.join("transcripts").join("stray.json"), at(97)).expect("a file");
        fs::create_dir(dir.path.join("transcripts").join("21").join("folder.json"))
            .expect("a directory");
        let names: Vec<String> = cache_shards(&dir.path)
            .expect("a cache that can be listed")
            .iter()
            .map(|(chrom, path)| {
                let name = path.file_name().expect("a name").to_string_lossy();
                format!("{chrom}/{name}")
            })
            .collect();
        assert_eq!(
            names,
            [
                "10/1-2.json",
                "21/._2-3.json",
                "21/10-11.json",
                "21/2-3.json",
                "21/A.json",
                "21/folder.json",
                "3/1-2.json"
            ]
        );
        for threads in [0, 1, 3, 64] {
            let authority = loaded(&dir, threads);
            assert_eq!(authority.span("21", "T"), Some((7, 9)), "{threads} threads");
            assert_eq!(authority.span("10", "T"), Some((1, 9)));
            assert_eq!(authority.model("3", "T").map(|model| model.start), Some(3));
            assert_eq!(authority.by_chr.len(), 3);
        }
    }

    #[test]
    fn the_last_record_of_an_id_stands_across_rounds_of_shards() {
        let count = 2 * SHARDS_A_ROUND + 7;
        let shards: Vec<(String, String)> = (0..count)
            .map(|at| {
                let record = format!(
                    r#"[{{"stable_id": "T", "start": {at}, "end": {count}, "strand": 1}},
                        {{"stable_id": "T{at}", "start": {at}, "end": {count}}}]"#
                );
                (format!("{at:05}.json"), record)
            })
            .collect();
        let named: Vec<(&str, &str, &[u8])> = shards
            .iter()
            .map(|(name, record)| ("21", name.as_str(), record.as_bytes()))
            .collect();
        let dir = cache_of(&named);
        let last = (count - 1) as i128;
        for threads in [1, 5] {
            let authority = loaded(&dir, threads);
            assert_eq!(authority.span("21", "T"), Some((last, count as i128)));
            assert_eq!(
                authority.model("21", "T").map(|model| model.start),
                Some(last)
            );
            assert_eq!(authority.by_chr["21"].len(), count + 1);
            assert_eq!(authority.span("21", "T300"), Some((300, count as i128)));
        }
    }

    #[test]
    fn a_later_record_without_a_span_leaves_the_earlier_span_and_structure() {
        let dir = cache_of(&[(
            "21",
            "1-2.json",
            br#"[{"stable_id": "T", "start": 5, "end": 9, "strand": 1},
                 {"stable_id": "T", "start": "x", "end": 20, "strand": 1},
                 {"stable_id": "U", "start": 1, "end": 2},
                 {"stable_id": "V"}]"#,
        )]);
        let authority = loaded(&dir, 2);
        assert_eq!(authority.span("21", "T"), Some((5, 9)));
        assert_eq!(authority.model("21", "T").map(|model| model.end), Some(9));
        assert_eq!(authority.span("21", "U"), Some((1, 2)));
        assert_eq!(authority.model("21", "U"), None);
        assert_eq!(authority.span("21", "V"), None);
        assert_eq!(authority.transcripts("21"), Some(&ids(&["T", "U", "V"])));
        assert_eq!(authority.transcripts("22"), None);
        let (by_chr, spans) = load_transcript_index(&dir.path, 2).expect("a cache");
        assert_eq!((by_chr, spans), (authority.by_chr, authority.spans));
    }

    /// A gzip member of the deflate method around `deflate`, with the checksum and the
    /// length of `text`.
    fn gzip_member(deflate: &[u8], text: &[u8]) -> Vec<u8> {
        let mut crc = Crc::new();
        crc.update(text);
        let mut member = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3];
        member.extend_from_slice(deflate);
        member.extend_from_slice(&crc.sum().to_le_bytes());
        member.extend_from_slice(&(text.len() as u32).to_le_bytes());
        member
    }

    fn gzip(text: &str) -> Vec<u8> {
        use std::io::Write;
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(text.as_bytes()).expect("gzip");
        encoder.finish().expect("gzip")
    }

    #[test]
    fn a_gzipped_shard_is_read_and_an_unreadable_one_holds_nothing() {
        let good = gzip(r#"[{"stable_id": "ENST_GZ", "start": 1, "end": 2}]"#);
        let padded = [gzip(r#"[{"stable_id": "ENST_PADDED"}]"#), vec![0; 7]].concat();
        let two = [
            gzip(r#"[{"stable_id": "ENST_ONE"},"#),
            vec![0; 3],
            gzip(r#" {"stable_id": "ENST_TWO"}]"#),
        ]
        .concat();
        let named = {
            let plain = gzip(r#"[{"stable_id": "ENST_NAMED"}]"#);
            let mut member = vec![0x1f, 0x8b, 8, 0xfc];
            member.extend_from_slice(&plain[4..10]);
            member.extend_from_slice(b"\x02\x00abname\x00comment\x00");
            member.extend_from_slice(&plain[10..]);
            member
        };
        let wrong_crc = {
            let mut member = good.clone();
            let at = member.len() - 8;
            member[at] ^= 0xff;
            member
        };
        let wrong_length = {
            let mut member = good.clone();
            let at = member.len() - 4;
            member[at] ^= 0xff;
            member
        };
        let other_method = {
            let mut member = good.clone();
            member[2] = 7;
            member
        };
        let dir = cache_of(&[
            ("21", "01.json.gz", &good),
            ("21", "02.json.gz", &padded),
            ("21", "03.json.gz", &two),
            ("21", "04.json.gz", &named),
            ("21", "05.json.gz", &good[..good.len() - 12]),
            ("21", "06.json.gz", &good[..good.len() - 3]),
            ("21", "07.json.gz", &wrong_crc),
            ("21", "08.json.gz", &wrong_length),
            ("21", "09.json.gz", &other_method),
            (
                "21",
                "10.json.gz",
                &[good.clone(), b"xyz".to_vec()].concat(),
            ),
            ("21", "11.json.gz", &[vec![0], good.clone()].concat()),
            (
                "21",
                "12.json.gz",
                br#"[{"stable_id": "ENST_PLAIN_IN_GZ"}]"#,
            ),
            ("21", "13.json.gz", b""),
            ("21", "14.json.gz", &good[..10]),
            ("21", "15.json", &good),
            ("21", "a.json", br#"[{"stable_id": "ENST_CUT"}"#),
            ("21", "b.json", br#"[{"stable_id": "ENST_TWICE"}] []"#),
            (
                "21",
                "c.json",
                b"[{\"stable_id\": \"ENST_LATIN\", \"x\": \"\xe9\"}]",
            ),
            (
                "21",
                "d.json",
                b"\xef\xbb\xbf[{\"stable_id\": \"ENST_BOM\"}]",
            ),
            ("21", "e.json", b""),
            ("21", "f.json", br#" [{"stable_id": "ENST_F"}] "#),
        ]);
        let got = load_transcripts_by_chromosome(&dir.path, 4).expect("a cache");
        assert_eq!(
            got["21"],
            ids(&[
                "ENST_GZ",
                "ENST_PADDED",
                "ENST_ONE",
                "ENST_TWO",
                "ENST_NAMED",
                "ENST_F"
            ])
        );
    }

    #[test]
    fn a_corrupt_deflate_stream_ends_the_loading() {
        let text = br#"[{"stable_id": "T"}]"#;
        let good = gzip(r#"[{"stable_id": "T"}]"#);
        let block_of_no_type = gzip_member(b"\x07\x00\x00", text);
        let stored_block_of_two_lengths = gzip_member(b"\x01\x05\x00\x00\x00hello", b"hello");
        let after_a_member = [good.clone(), block_of_no_type.clone()].concat();
        let after_bytes_that_are_no_text =
            [gzip("["), vec![0xff; 0], block_of_no_type.clone()].concat();
        for corrupt in [
            &block_of_no_type,
            &stored_block_of_two_lengths,
            &after_a_member,
            &after_bytes_that_are_no_text,
        ] {
            let dir = cache_of(&[("21", "1.json", text), ("21", "2.json.gz", corrupt)]);
            let failure = load_cache_authority(&dir.path, 2).expect_err("a corrupt stream");
            assert_eq!(failure.message, "Error -3 while decompressing data");
        }
        // The same bytes under a name that is not gzipped are bytes that are no text.
        let dir = cache_of(&[("21", "1.json", text), ("21", "2.json", &block_of_no_type)]);
        assert_eq!(loaded(&dir, 1).by_chr["21"], ids(&["T"]));
    }

    #[test]
    fn a_transcripts_directory_that_cannot_be_listed_ends_the_loading() {
        use std::os::unix::fs::PermissionsExt;
        let lock = |path: &Path, mode: u32| {
            fs::set_permissions(path, fs::Permissions::from_mode(mode)).expect("a mode");
        };
        let dir = cache_of(&[
            ("21", "1.json", br#"[{"stable_id": "T"}]"#),
            ("22", "1.json", br#"[{"stable_id": "U"}]"#),
        ]);
        let root = dir.path.join("transcripts");
        // A chromosome's directory that cannot be listed, and a shard that cannot be
        // opened, hold nothing.
        lock(&root.join("21"), 0o000);
        lock(&root.join("22").join("1.json"), 0o000);
        let unlisted = fs::read_dir(root.join("21")).is_err();
        let got = load_cache_authority(&dir.path, 1);
        lock(&root.join("21"), 0o755);
        if unlisted {
            let got = got.expect("a cache Python reads");
            assert_eq!(got.by_chr, TranscriptsByChromosome::new());
        }
        lock(&root, 0o000);
        let unlisted = fs::read_dir(&root).is_err();
        let got = load_cache_authority(&dir.path, 1);
        lock(&root, 0o755);
        // A user the mode does not bind lists the directory all the same.
        if unlisted {
            let failure = got.expect_err("a directory that cannot be listed");
            assert_eq!(failure.kind, crate::failure::FailureKind::Io);
            assert_eq!(
                failure.message,
                format!("[Errno 13] Permission denied: '{}'", root.display())
            );
        }
        // The path is quoted as Python's `repr` quotes a text: in double quotes where it
        // holds a single one.
        let quoted = dir.path.join("it's");
        let root = quoted.join("transcripts");
        fs::create_dir_all(root.join("21")).expect("a directory");
        lock(&root, 0o000);
        let unlisted = fs::read_dir(&root).is_err();
        let got = load_cache_authority(&quoted, 1);
        lock(&root, 0o755);
        if unlisted {
            assert_eq!(
                got.expect_err("a directory that cannot be listed").message,
                format!("[Errno 13] Permission denied: \"{}\"", root.display())
            );
        }
        let file = TempDir::new("vep_sv_compare_cache_").expect("a scratch directory");
        fs::write(file.path.join("transcripts"), "x").expect("a file");
        assert_eq!(
            load_cache_authority(&file.path, 1),
            Ok(CacheAuthority::default())
        );
    }

    #[test]
    fn a_shard_is_an_array_or_an_object_of_arrays() {
        let named = |text: &str| -> Vec<String> {
            records(text)
                .into_iter()
                .map(|record| record.stable_id.unwrap_or_else(|| "(no text)".to_string()))
                .collect()
        };
        assert_eq!(
            named(
                r#"[{"stable_id": "A"}, null, 5, "B", [{"stable_id": "C"}], {}, {"stable_id": "D"}]"#
            ),
            ["A", "D"]
        );
        assert_eq!(
            named(
                r#"{"x": [{"stable_id": "A"}], "y": null, "z": {"stable_id": "B"}, "v": "text",
                    "u": 0, "t": false, "w": [{"stable_id": "C"}]}"#
            ),
            ["A", "C"]
        );
        // A repeated key keeps its first place and its last value.
        assert_eq!(
            named(
                r#"{"x": [{"stable_id": "A"}], "y": [{"stable_id": "B"}], "x": [{"stable_id": "C"}]}"#
            ),
            ["C", "B"]
        );
        assert_eq!(named(r#"{"x": 5, "x": [{"stable_id": "A"}]}"#), ["A"]);
        for scalar in [
            "5",
            "null",
            r#""text""#,
            "true",
            "{}",
            "[]",
            "NaN",
            "-Infinity",
            "-0",
            "1e5",
        ] {
            assert_eq!(named(scalar), [""; 0], "{scalar}");
        }
        for broken in [
            "",
            "[",
            "[1,]",
            "{'a': 1}",
            "[] x",
            "\u{feff}[]",
            "[01]",
            "[1.]",
            "[.5]",
            "[1e]",
            "[+1]",
            "[-]",
            "[nul]",
            "[nan]",
            "[-NaN]",
            "[+Infinity]",
            "[1 2]",
            r#"{"a" 1}"#,
            r#"{"a": 1,}"#,
            r#"{1: 2}"#,
            r#"["a\qb"]"#,
            r#"["a\u12g4"]"#,
            "[\"a\tb\"]",
            r#"["a"#,
            "[[]",
            "[]]",
            "{}}",
            r#"{"a": }"#,
            "[,1]",
            "\u{a0}[]",
            "[\u{2003}]",
        ] {
            assert_eq!(parse_shard(broken), Ok(None), "{broken}");
        }
        for read in [
            " \t\r\n[ ] \n",
            "[1.5e-7, -0.0, 1E+2, 0e0, -1]",
            concat!(r#"["\"\\\/\b\f\n\r\té"#, "\u{1f600}", r#""]"#),
            r#"{"a": {"b": [[], {}, [{"c": null}]]}, "d": [true, false]}"#,
            r#"["\ud800", "\udc00x", "\ud800A"]"#,
        ] {
            assert_eq!(named(read), [""; 0], "{read}");
        }
    }

    #[test]
    fn a_shard_is_read_as_pythons_json_reads_it() {
        let id_of = |fields: &str| record(fields).stable_id;
        let held = Some("T".to_string());
        // Three numbers beside JSON's, in a field the comparison does not read.
        assert_eq!(
            id_of(r#""x": [NaN, Infinity, -Infinity, 1e999, -1e999]"#),
            held
        );
        // An integer is a number up to 4,300 digits, a number with a fraction at any length.
        let digits = "1".repeat(INT_MAX_STR_DIGITS);
        assert_eq!(id_of(&format!(r#""x": {digits}"#)), held);
        assert_eq!(id_of(&format!(r#""x": -{digits}"#)), held);
        assert_eq!(id_of(&format!(r#""x": {digits}1.5"#)), held);
        assert_eq!(id_of(&format!(r#""x": {digits}1e5"#)), held);
        let longer = format!(r#"[{{"stable_id": "T", "x": [{{"y": {digits}1}}]}}]"#);
        assert_eq!(parse_shard(&longer), Ok(None));
        // A surrogate that is no half of a pair is kept, in a key and in a value.
        let lone = record(concat!(
            r#""\udc00": 1, "biotype": "a\ud800b"#,
            "\u{1f600}\""
        ));
        assert_eq!(lone.biotype, "a\u{fffd}b\u{1f600}");
        // A key is compared after its escapes are read, and the last of a repeated one stands.
        assert_eq!(
            record(r#""start": 5, "start": 6, "start": 7"#).start,
            Some(7)
        );
        // A field that is not read may be of any depth.
        let deep = format!(r#""x": {}{}"#, "[".repeat(10_000), "]".repeat(10_000));
        assert_eq!(id_of(&deep), held);
        let deep = format!(
            r#""x": {}1{}"#,
            r#"{"k": ["#.repeat(5_000),
            "]}".repeat(5_000)
        );
        assert_eq!(id_of(&deep), held);
    }

    #[test]
    fn a_stable_id_is_what_python_reads_as_true() {
        let kept = |id: &str| -> Option<Option<String>> {
            let text = format!(r#"[{{"stable_id": {id}}}]"#);
            records(&text).pop().map(|record| record.stable_id)
        };
        assert_eq!(kept(r#""ENST1""#), Some(Some("ENST1".to_string())));
        for absent in [
            r#""""#, "null", "false", "0", "0.0", "-0", "-0.0", "[]", "{}", "[ ]", "{\n}",
        ] {
            assert_eq!(kept(absent), None, "{absent}");
        }
        for other in ["5", "-1", "0.5", "true", "NaN", "Infinity", "1e999"] {
            assert_eq!(kept(other), Some(None), "{other}");
        }
        // The last of a repeated key stands.
        let twice = records(
            r#"[{"stable_id": "A", "stable_id": "B"}, {"stable_id": "A", "stable_id": null},
                {"stable_id": ["a"], "stable_id": "C"}]"#,
        );
        let twice: Vec<Option<String>> = twice.into_iter().map(|record| record.stable_id).collect();
        assert_eq!(twice, [Some("B".to_string()), Some("C".to_string())]);
    }

    #[test]
    fn an_id_that_is_no_text_makes_its_chromosome_known_and_files_no_id() {
        let dir = cache_of(&[
            (
                "21",
                "1-2.json",
                br#"[{"stable_id": 7, "start": 1, "end": 2}]"#,
            ),
            ("22", "1-2.json", br#"[{"stable_id": true}]"#),
            (
                "X",
                "1-2.json",
                br#"[{"stable_id": 0, "start": 1, "end": 2}]"#,
            ),
        ]);
        let authority = loaded(&dir, 1);
        assert_eq!(authority.transcripts("21"), Some(&ids(&[])));
        assert_eq!(authority.transcripts("22"), Some(&ids(&[])));
        assert_eq!(authority.transcripts("X"), None);
        assert_eq!(authority.spans.keys().collect::<Vec<_>>(), ["21"]);
        assert!(authority.models.is_empty());
    }

    #[test]
    fn an_id_that_is_a_number_is_counted_once_a_chromosome_and_enters_no_view() {
        // One record under a number beside two under texts.
        let dir = cache_of(&[
            (
                "21",
                "0-1.json",
                br#"[{"stable_id": 7, "start": 500, "end": 3000}]"#,
            ),
            (
                "22",
                "0-1.json",
                br#"[{"stable_id": "ENST22_X", "start": 16150000, "end": 16160000},
                     {"stable_id": "ENST22_M", "start": 26036000, "end": 26040000}]"#,
            ),
        ]);
        let authority = loaded(&dir, 1);
        assert_eq!(authority.by_chr.len(), 2);
        assert_eq!(authority.transcript_total(), 3);
        assert_eq!(authority.transcripts("21"), Some(&ids(&[])));
        assert_eq!(
            authority.transcripts("22"),
            Some(&ids(&["ENST22_X", "ENST22_M"]))
        );
        assert_eq!(authority.spans["21"], HashMap::new());
        assert_eq!(authority.spans["22"].len(), 2);
        assert!(authority.models.is_empty());
    }

    #[test]
    fn ids_that_are_numbers_are_one_id_where_python_holds_them_equal() {
        // Chromosome 21 holds 15 ids to Python: the texts `7` and `T`, and 7, 1, 1e23, NaN,
        // the two infinities, 0.5, -7, 2^128, 2^128 + 1, 1.5, 100 and 8.
        let first = br#"[
            {"stable_id": 7, "start": 500, "end": 3000}, {"stable_id": 7.0}, {"stable_id": true},
            {"stable_id": 1}, {"stable_id": 1.0}, {"stable_id": 1e23},
            {"stable_id": 99999999999999991611392}, {"stable_id": NaN}, {"stable_id": NaN},
            {"stable_id": Infinity}, {"stable_id": 1e999}, {"stable_id": -Infinity},
            {"stable_id": 0.5}, {"stable_id": 5e-1}, {"stable_id": "7"}, {"stable_id": -7},
            {"stable_id": -7.0}, {"stable_id": 340282366920938463463374607431768211456},
            {"stable_id": 340282366920938463463374607431768211456},
            {"stable_id": 340282366920938463463374607431768211457},
            {"stable_id": 3.4028236692093846346e38}, {"stable_id": 0}, {"stable_id": 0.0},
            {"stable_id": -0.0}, {"stable_id": false}, {"stable_id": ""}, {"stable_id": null},
            {"stable_id": 1.5}, {"stable_id": 15e-1}, {"stable_id": -1e999}, {"stable_id": 1E2},
            {"stable_id": 100}]"#;
        let dir = cache_of(&[
            ("21", "0-1.json", first),
            (
                "21",
                "1-2.json",
                br#"[{"stable_id": NaN}, {"stable_id": 7}, {"stable_id": 8}, {"stable_id": "T"}]"#,
            ),
            (
                "22",
                "0-1.json",
                br#"[{"stable_id": "ENST22_X", "start": 16150000, "end": 16160000},
                     {"stable_id": 7}, {"stable_id": true}]"#,
            ),
            // A repeated key of an object keeps its last array.
            (
                "X",
                "0-1.json",
                br#"{"X": [{"stable_id": 2.5}], "Y": [{"stable_id": 3}],
                     "X": [{"stable_id": 2.50}, {"stable_id": 4}]}"#,
            ),
        ]);
        for threads in [1, 4] {
            let authority = loaded(&dir, threads);
            let numbers = |chrom: &str| authority.number_ids[chrom].len();
            assert_eq!(
                (numbers("21"), numbers("22"), numbers("X")),
                (13, 2, 3),
                "{threads}"
            );
            assert_eq!(authority.transcript_total(), 21, "{threads}");
            assert_eq!(authority.transcripts("21"), Some(&ids(&["7", "T"])));
            assert_eq!(authority.transcripts("22"), Some(&ids(&["ENST22_X"])));
            assert_eq!(authority.transcripts("X"), Some(&ids(&[])));
            assert_eq!(authority.spans["21"], HashMap::new());
            assert_eq!(authority.spans["22"].len(), 1);
            assert!(authority.models.is_empty());
        }
    }

    #[test]
    fn the_text_of_an_id_that_is_a_number_is_that_of_its_value() {
        let text = |id: &str| -> Option<String> {
            let shard = format!(r#"[{{"stable_id": {id}}}]"#);
            records(&shard).pop().and_then(|record| record.number_id)
        };
        for (id, expected) in [
            ("7", "7"),
            ("7.0", "7"),
            ("-7.0", "-7"),
            ("true", "1"),
            ("1.0", "1"),
            ("1E2", "100"),
            ("1e23", "99999999999999991611392"),
            ("99999999999999991611392", "99999999999999991611392"),
            (
                "3.4028236692093846346e38",
                "340282366920938463463374607431768211456",
            ),
            (
                "340282366920938463463374607431768211457",
                "340282366920938463463374607431768211457",
            ),
            ("0.5", "5e-1"),
            ("5e-1", "5e-1"),
            ("15e-1", "1.5e0"),
            ("NaN", "nan"),
            ("Infinity", "inf"),
            ("1e999", "inf"),
            ("-Infinity", "-inf"),
            ("-1e999", "-inf"),
        ] {
            assert_eq!(text(id).as_deref(), Some(expected), "{id}");
        }
        for id in [r#""7""#, "0", "0.0", "-0.0", "false", "null", r#""""#] {
            assert_eq!(text(id), None, "{id}");
        }
    }

    #[test]
    fn an_integer_is_a_number_or_a_text_python_reads_as_one() {
        let start = |value: &str| record(&format!(r#""start": {value}"#)).start;
        for (value, expected) in [
            ("12", Some(12)),
            ("-12", Some(-12)),
            ("12.9", Some(12)),
            ("-12.9", Some(-12)),
            ("1e3", Some(1000)),
            ("-0", Some(0)),
            ("true", Some(1)),
            ("false", Some(0)),
            (r#""12""#, Some(12)),
            (r#"" +1_2\n""#, Some(12)),
            (r#""-12""#, Some(-12)),
            (r#""١٢""#, Some(12)),
            (r#""1２""#, Some(12)),
            ("18446744073709551616", Some(1 << 64)),
            ("-9223372036854775809", Some(-(1 << 63) - 1)),
            ("1e30", Some(1_000_000_000_000_000_019_884_624_838_656)),
            ("170141183460469231731687303715884105727", Some(i128::MAX)),
            ("-170141183460469231731687303715884105728", Some(i128::MIN)),
            (
                r#""170141183460469231731687303715884105727""#,
                Some(i128::MAX),
            ),
            ("170141183460469231731687303715884105728", None),
            ("-170141183460469231731687303715884105729", None),
            (r#""170141183460469231731687303715884105728""#, None),
            (
                "1.7e38",
                Some(169_999_999_999_999_998_061_923_293_023_115_935_744),
            ),
            ("1.8e38", None),
            ("-1.8e38", None),
            ("1e308", None),
            ("NaN", None),
            ("null", None),
            (r#""""#, None),
            (r#""12.0""#, None),
            (r#""0x10""#, None),
            (r#""1__2""#, None),
            (r#""²""#, None),
            (r#""Infinity""#, None),
            (r#""NaN""#, None),
            ("[12]", None),
            (r#"{"start": 12}"#, None),
        ] {
            assert_eq!(start(value), expected, "{value}");
        }
        assert_eq!(record(r#""start": 1, "start": 2"#).start, Some(2));
        assert_eq!(record(r#""start": 1, "start": null"#).start, None);
    }

    #[test]
    fn a_record_gives_its_exons_pairs_and_mature_ranges() {
        let got = record(
            r#""biotype": "miRNA", "end": 30, "strand": -1, "start": 10,
            "exons": [{"start": 20, "end": 30, "phase": 0}, {"start": 10}, {},
                      {"end": "15", "start": 10.5}, {"start": NaN, "end": 4}],
            "variation_effect_feature_cache": {"peptide": "x", "mapper": {
              "cdna_coding_end": 9, "cdna_coding_start": 2, "pair_count": 2,
              "pairs": [{"from_start": 1, "from_end": 11, "to_start": 20, "to_end": 30, "ori": -1},
                        {"from_start": 12, "from_end": 17, "to_start": 10, "to_end": 15},
                        {}]}},
            "attributes": [{"code": "miRNA", "value": "3-7"}, {"value": "8-9", "code": "miRNA"},
                           {"code": "gencode_basic", "value": "1-2"}, {"code": "miRNA", "value": 5},
                           {"code": "miRNA"}, {}]"#,
        );
        let expected = CacheRecord {
            stable_id: Some("T".to_string()),
            number_id: None,
            start: Some(10),
            end: Some(30),
            strand: Some(-1),
            biotype: "miRNA".to_string(),
            exons: vec![(20, 30), (10, 15)],
            cdna_coding_start: Some(2),
            cdna_coding_end: Some(9),
            pairs: vec![MapperPair {
                from_start: 1,
                from_end: 11,
                to_start: 20,
                to_end: 30,
                ori: -1,
            }],
            mirna: vec!["3-7".to_string(), "8-9".to_string()],
        };
        assert_eq!(got, expected);
        let model = model_from_cache_record("21", &got).expect("a model");
        assert_eq!(model.exons, [(10, 15), (20, 30)]);
        assert_eq!(model.cds, Some((22, 29)));
        assert_eq!(model.mirna_mature, [(3, 7), (8, 9)]);
    }

    #[test]
    fn a_field_python_reads_as_false_is_absent() {
        let plain = CacheRecord {
            stable_id: Some("T".to_string()),
            start: Some(1),
            end: Some(2),
            strand: Some(1),
            ..CacheRecord::default()
        };
        for fields in [
            r#""exons": null, "variation_effect_feature_cache": null, "attributes": null"#,
            r#""exons": 0, "variation_effect_feature_cache": 0, "attributes": 0.0"#,
            r#""exons": false, "variation_effect_feature_cache": false, "attributes": false"#,
            r#""exons": "", "variation_effect_feature_cache": "", "attributes": """#,
            r#""exons": [], "variation_effect_feature_cache": [], "attributes": []"#,
            r#""exons": {}, "variation_effect_feature_cache": {}, "attributes": {}"#,
            r#""exons": [{}], "attributes": [{}], "biotype": null"#,
            r#""variation_effect_feature_cache": {"mapper": null}, "biotype": 5"#,
            r#""variation_effect_feature_cache": {"mapper": 0}, "biotype": ["miRNA"]"#,
            r#""variation_effect_feature_cache": {"mapper": []}, "biotype": {"a": 1}"#,
            r#""variation_effect_feature_cache": {"mapper": {"pairs": null}}"#,
            r#""variation_effect_feature_cache": {"mapper": {"cdna_coding_start": 1}, "mapper": 0}"#,
        ] {
            assert_eq!(
                records(&modelled(fields)),
                std::slice::from_ref(&plain),
                "{fields}"
            );
        }
    }

    /// Each shape on which Python's reading of a record ends in an exception, with the
    /// text of that exception.
    #[test]
    fn a_shape_python_raises_on_ends_the_loading_with_the_exceptions_text() {
        assert_eq!(
            failure(r#"[{"stable_id": ["a"]}]"#),
            "cannot use 'list' as a set element (unhashable type: 'list')"
        );
        assert_eq!(
            failure(r#"[{"stable_id": {"a": 1}}]"#),
            "cannot use 'dict' as a set element (unhashable type: 'dict')"
        );
        for infinite in ["Infinity", "-Infinity", "1e999", "-1e999", "1.0E+400"] {
            for key in ["start", "end", "strand"] {
                let text = format!(r#"[{{"stable_id": "T", "{key}": {infinite}}}]"#);
                assert_eq!(failure(&text), NOT_AN_INTEGER, "{text}");
            }
            for fields in [
                format!(r#""exons": [{{"start": 1, "end": 2}}, {{"start": {infinite}}}]"#),
                format!(r#""exons": [{{"end": {infinite}, "start": "x"}}]"#),
                format!(
                    r#""variation_effect_feature_cache": {{"mapper": {{"cdna_coding_start": {infinite}}}}}"#
                ),
                format!(
                    r#""variation_effect_feature_cache": {{"mapper": {{"cdna_coding_start": 0, "cdna_coding_end": {infinite}}}}}"#
                ),
                format!(
                    r#""variation_effect_feature_cache": {{"mapper": {{"cdna_coding_start": 1, "cdna_coding_end": 1, "pairs": [{{"ori": {infinite}}}]}}}}"#
                ),
            ] {
                assert_eq!(record_failure(&fields), NOT_AN_INTEGER, "{fields}");
            }
        }
        let bounds = r#""cdna_coding_start": 1, "cdna_coding_end": 1"#;
        let listed = |field: &str, value: &str| match field {
            "pairs" => format!(
                r#""variation_effect_feature_cache": {{"mapper": {{{bounds}, "pairs": {value}}}}}"#
            ),
            _ => format!(r#""{field}": {value}"#),
        };
        for field in ["exons", "pairs", "attributes"] {
            for (value, kind) in [
                ("5", "int"),
                ("-1", "int"),
                ("0.5", "float"),
                ("NaN", "float"),
                ("true", "bool"),
            ] {
                assert_eq!(
                    record_failure(&listed(field, value)),
                    not_iterable_text(kind),
                    "{field} {value}"
                );
            }
            for (value, kind) in [
                (r#""ab""#, "str"),
                (r#"{"a": 1}"#, "str"),
                ("[null]", "NoneType"),
                ("[5]", "int"),
                ("[0]", "int"),
                ("[0.5]", "float"),
                ("[true]", "bool"),
                ("[false]", "bool"),
                (r#"["s"]"#, "str"),
                (r#"[""]"#, "str"),
                ("[[1]]", "list"),
                ("[[]]", "list"),
                ("[{}, null, 5]", "NoneType"),
            ] {
                assert_eq!(
                    record_failure(&listed(field, value)),
                    no_get_text(kind),
                    "{field} {value}"
                );
            }
        }
        for (value, kind) in [
            (r#""x""#, "str"),
            ("[1]", "list"),
            ("5", "int"),
            ("0.5", "float"),
            ("true", "bool"),
        ] {
            let cache = format!(r#""variation_effect_feature_cache": {value}"#);
            assert_eq!(record_failure(&cache), no_get_text(kind), "{cache}");
            let mapper = format!(r#""variation_effect_feature_cache": {{"mapper": {value}}}"#);
            assert_eq!(record_failure(&mapper), no_get_text(kind), "{mapper}");
        }
        let long = "1".repeat(INT_MAX_STR_DIGITS + 1);
        let limit = |digits: usize| crate::text::int_conversion_failure(digits).message;
        for (value, digits) in [
            (format!("{long}-5"), 4301),
            (format!("5-{long}"), 4301),
            (format!(" {long}1-{long} "), 4302),
            (format!("{}-5", "\u{0663}".repeat(4301)), 4301),
        ] {
            let attributes = format!(r#""attributes": [{{"value": "{value}", "code": "miRNA"}}]"#);
            assert_eq!(record_failure(&attributes), limit(digits));
        }
        for kind in ["5", "-1", "true", "0.5", "NaN", "Infinity"] {
            let name = match kind {
                "5" | "-1" => "int",
                "true" => "bool",
                _ => "float",
            };
            let text = format!(r#"{{"21": [{{"stable_id": "T"}}], "22": {kind}}}"#);
            assert_eq!(failure(&text), not_iterable_text(name), "{text}");
        }
    }

    /// Beside each shape Python raises on, the shapes it reads to their end.
    #[test]
    fn a_shape_python_does_not_reach_ends_nothing() {
        let read = |text: &str| records(text).len();
        // A record without an id is not looked into.
        assert_eq!(
            read(r#"[{"start": Infinity, "exons": [null], "attributes": 5}]"#),
            0
        );
        assert_eq!(
            read(r#"[{"stable_id": "", "start": Infinity}, {"stable_id": [], "end": 1e999}]"#),
            0
        );
        // A record without an integer start, end or strand has no structure to read.
        for lacking in [
            r#""end": 2, "strand": 1"#,
            r#""start": 1, "end": 2, "strand": "x""#,
            r#""start": 1, "end": NaN, "strand": 1"#,
            r#""start": 1e308, "end": 2, "strand": 1"#,
        ] {
            let text = format!(
                r#"[{{"stable_id": "T", {lacking}, "exons": [null], "attributes": 5,
                "variation_effect_feature_cache": "x"}}]"#
            );
            assert_eq!(read(&text), 1, "{text}");
        }
        // The pairs are read only under two coding bounds above zero.
        for bounds in [
            "",
            r#""cdna_coding_start": 0, "cdna_coding_end": 5,"#,
            r#""cdna_coding_start": 5, "cdna_coding_end": -1,"#,
            r#""cdna_coding_start": 5, "cdna_coding_end": NaN,"#,
            r#""cdna_coding_start": "x", "cdna_coding_end": 5,"#,
        ] {
            for pairs in ["5", "[null]", r#"[{"ori": Infinity}]"#, r#""ab""#] {
                let fields = format!(
                    r#""variation_effect_feature_cache": {{"mapper": {{{bounds} "pairs": {pairs}}}}}"#
                );
                assert_eq!(read(&modelled(&fields)), 1, "{fields}");
            }
        }
        // The last of a repeated key stands, whatever the earlier one held.
        for fields in [
            r#""exons": [null], "exons": []"#,
            r#""start": Infinity, "start": 5"#,
            r#""attributes": true, "attributes": null"#,
            r#""variation_effect_feature_cache": "x", "variation_effect_feature_cache": {}"#,
            r#""variation_effect_feature_cache": {"mapper": 5, "mapper": {}}"#,
            r#""exons": [{"start": Infinity, "start": 1, "end": 2}]"#,
        ] {
            assert_eq!(read(&modelled(fields)), 1, "{fields}");
        }
        // A mature range is converted for the code `miRNA` alone, up to 4,300 digits a run.
        let long = "1".repeat(INT_MAX_STR_DIGITS + 1);
        for attribute in [
            format!(r#"{{"code": "x", "value": "{long}-5"}}"#),
            format!(r#"{{"value": "{long}-5"}}"#),
            format!(r#"{{"code": "miRNA", "value": "{}-5"}}"#, &long[1..]),
            format!(r#"{{"code": "miRNA", "value": "{long}"}}"#),
            format!(r#"{{"code": "miRNA", "value": "{long}-x"}}"#),
            format!(r#"{{"code": "miRNA", "value": {long}.5}}"#),
        ] {
            let fields = format!(r#""attributes": [{attribute}]"#);
            assert_eq!(read(&modelled(&fields)), 1);
        }
        // A text that is no number, where an integer is read, is no integer.
        assert_eq!(record(&format!(r#""start": "{long}""#)).start, None);
        // A value of an object of arrays that Python reads as false, or runs through to no
        // record, ends nothing.
        for value in [
            "0",
            "0.0",
            "null",
            "false",
            r#""abc""#,
            r#"{"a": [{"stable_id": "U"}]}"#,
        ] {
            let text = format!(r#"{{"21": [{{"stable_id": "T"}}], "22": {value}}}"#);
            assert_eq!(read(&text), 1, "{text}");
        }
    }

    #[test]
    fn the_first_exception_in_pythons_order_decides() {
        // Within a record: the id, the three bounds, the exons, the mapper, the attributes,
        // whatever the order of its keys.
        for (fields, expected) in [
            (
                r#""attributes": true, "exons": [null]"#,
                no_get_text("NoneType"),
            ),
            (
                r#""attributes": [null], "variation_effect_feature_cache": "x""#,
                no_get_text("str"),
            ),
            (
                r#""variation_effect_feature_cache": [1], "exons": [5]"#,
                no_get_text("int"),
            ),
            (
                r#""attributes": 5, "variation_effect_feature_cache": {"mapper": {"cdna_coding_start": 1, "cdna_coding_end": 1, "pairs": true}}"#,
                not_iterable_text("bool"),
            ),
            (
                r#""variation_effect_feature_cache": {"mapper": {"pairs": true, "cdna_coding_end": Infinity, "cdna_coding_start": 1}}"#,
                NOT_AN_INTEGER.to_string(),
            ),
            (
                r#""exons": [{"start": 1}, 5, null, {"end": Infinity}]"#,
                no_get_text("int"),
            ),
            (
                r#""exons": [{"end": Infinity}, 5]"#,
                NOT_AN_INTEGER.to_string(),
            ),
        ] {
            assert_eq!(record_failure(fields), expected, "{fields}");
        }
        assert_eq!(
            failure(
                r#"[{"exons": [null], "strand": Infinity, "start": 1, "end": 2, "stable_id": "T"}]"#
            ),
            NOT_AN_INTEGER
        );
        assert_eq!(
            failure(r#"[{"start": Infinity, "exons": [null], "stable_id": ["a"]}]"#),
            "cannot use 'list' as a set element (unhashable type: 'list')"
        );
        // Across records, the first; a shard that is no text Python reads holds nothing,
        // whatever its records.
        let second = modelled(r#""exons": [5]"#);
        let second = &second[1..second.len() - 1];
        let two = format!(r#"[{{"stable_id": "A"}}, {second}, {{"stable_id": ["c"]}}]"#);
        assert_eq!(failure(&two), no_get_text("int"));
        assert_eq!(
            parse_shard(&format!(r#"[{second}, {{"stable_id": "C""#)),
            Ok(None)
        );
        assert_eq!(parse_shard(&format!("[{second}] x")), Ok(None));
        assert_eq!(
            parse_shard(&format!(r#"[{second}, {{"x": 1e}}]"#)),
            Ok(None)
        );
        // In an object of arrays every value is run through before a record is read.
        assert_eq!(
            failure(r#"{"a": [{"stable_id": ["c"]}], "b": 5, "c": true}"#),
            not_iterable_text("int")
        );
        assert_eq!(
            failure(r#"{"a": [{"stable_id": "T"}], "b": [{"stable_id": ["c"]}], "a": 0.5}"#),
            not_iterable_text("float")
        );
        assert_eq!(
            failure(
                r#"{"a": [{"stable_id": {"x": 1}}], "b": [{"stable_id": ["c"]}], "a": [{"stable_id": ["d"]}], "a": []}"#
            ),
            "cannot use 'list' as a set element (unhashable type: 'list')"
        );
        // Across shards, the first in the order of the names, at any number of threads.
        let good = br#"[{"stable_id": "A"}]"#;
        let dir = cache_of(&[
            ("21", "1.json", good),
            (
                "21",
                "2.json",
                br#"[{"stable_id": "T", "start": Infinity}]"#,
            ),
            ("21", "3.json", b"["),
            ("3", "1.json", br#"[{"stable_id": ["c"]}]"#),
            ("3", "2.json", good),
        ]);
        for threads in [1, 2, 8] {
            let failure = load_cache_authority(&dir.path, threads).expect_err("a failure");
            assert_eq!(failure.message, NOT_AN_INTEGER, "{threads} threads");
            assert_eq!(failure.kind, crate::failure::FailureKind::Value);
        }
        assert!(load_transcript_index(&dir.path, 2).is_err());
        assert!(load_transcripts_by_chromosome(&dir.path, 2).is_err());
    }
}
