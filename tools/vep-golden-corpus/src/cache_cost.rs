// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What a candidate row costs the fixture cache under `select --cache`: the compressed
//! length of the transcript records of the whole cache that overlap the row within the
//! flank.

use crate::error::Failure;
use crate::json::{dumps, loads, Style, Value};
use crate::pytext::{is_dir, os_error_text, path_join, py_int, utf8_decode_error};
use crate::rows::{ByteSource, Row};
use flate2::{Compress, Compression, FlushCompress, Status};
use regex::Regex;
use std::cell::RefCell;
use std::collections::{BTreeSet, HashMap};
use std::sync::atomic::{AtomicUsize, Ordering};

/// The most threads that read shards ahead: each holds the text of one shard, and a shard
/// of a whole cache can be tens of megabytes.
const MOST_THREADS: usize = 16;

/// The length of `zlib.compress(data)`: the zlib stream of stock zlib's deflate at level 6
/// under a 32 KiB window and memory level 8, header and checksum included.
pub fn compressed_length(data: &[u8]) -> usize {
    let mut deflate = Compress::new(Compression::new(6), true);
    let mut out = [0u8; 64 * 1024];
    loop {
        let read = deflate.total_in() as usize;
        let status = deflate
            .compress(&data[read..], &mut out, FlushCompress::Finish)
            .expect("deflate of bytes held in memory");
        if status == Status::StreamEnd {
            return deflate.total_out() as usize;
        }
    }
}

/// The cost of one transcript record: the compressed length of its text as
/// `json.dumps(record, separators=(",", ":"), ensure_ascii=False)` writes it, in UTF-8.
pub fn record_cost(record: &Value) -> usize {
    compressed_length(dumps(record, &Style::COMPACT_UNICODE).as_bytes())
}

/// `int(value)` of a record's `start` or `end`: an integer, a text `int()` reads, a number
/// with a fraction cut towards zero, `True` and `False` as 1 and 0.
fn coordinate(value: Option<&Value>) -> Option<i64> {
    match value? {
        Value::Int(number) => number.as_i64(),
        Value::Str(text) => py_int(text),
        Value::Bool(flag) => Some(i64::from(*flag)),
        Value::Float(number) if number.is_finite() && number.abs() < 9.0e18 => {
            Some(number.trunc() as i64)
        }
        _ => None,
    }
}

/// The transcripts of one shard: start, end and cost of each record.
pub type Transcripts = Vec<(i64, i64, usize)>;

/// The offset of the first byte from `at` that JSON does not read as blank.
fn past_blanks(bytes: &[u8], mut at: usize) -> usize {
    while matches!(bytes.get(at), Some(b' ' | b'\t' | b'\n' | b'\r')) {
        at += 1;
    }
    at
}

/// The offset behind the string whose opening quote stands at `at`; `None` for a string
/// the text ends in.
fn past_string(bytes: &[u8], at: usize) -> Option<usize> {
    let mut at = at + 1;
    loop {
        at += memchr::memchr2(b'"', b'\\', bytes.get(at..)?)?;
        if bytes[at] == b'"' {
            return Some(at + 1);
        }
        at += 2;
    }
}

/// The offset behind the JSON value that starts at `at`, as far as its brackets and
/// strings tell: an object or a list up to the bracket that closes it, a string up to its
/// closing quote, anything else up to a blank, a comma or a closing bracket. Whether the
/// text between is a value is for its reader to say.
fn past_value(bytes: &[u8], at: usize) -> Option<usize> {
    match bytes.get(at)? {
        b'"' => past_string(bytes, at),
        b'{' | b'[' => {
            let mut depth = 0usize;
            let mut at = at;
            loop {
                match bytes.get(at)? {
                    b'"' => {
                        at = past_string(bytes, at)?;
                        continue;
                    }
                    b'{' | b'[' => depth += 1,
                    b'}' | b']' => depth -= 1,
                    _ => {}
                }
                at += 1;
                if depth == 0 {
                    return Some(at);
                }
            }
        }
        _ => {
            let rest = &bytes[at..];
            let end = rest
                .iter()
                .position(|byte| matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | b',' | b']' | b'}'));
            Some(at + end.unwrap_or(rest.len()))
        }
    }
}

/// The texts of the items of a JSON list, in order, so that a shard is read one record at
/// a time; `None` for a text that is no list between blanks.
fn list_items(text: &str) -> Option<Vec<&str>> {
    let bytes = text.as_bytes();
    let mut at = past_blanks(bytes, 0);
    if bytes.get(at) != Some(&b'[') {
        return None;
    }
    at = past_blanks(bytes, at + 1);
    let mut items = Vec::new();
    if bytes.get(at) == Some(&b']') {
        at += 1;
    } else {
        loop {
            let end = past_value(bytes, at).filter(|end| *end > at)?;
            items.push(&text[at..end]);
            at = past_blanks(bytes, end);
            match bytes.get(at)? {
                b',' => at = past_blanks(bytes, at + 1),
                b']' => {
                    at += 1;
                    break;
                }
                _ => return None,
            }
        }
    }
    (past_blanks(bytes, at) == bytes.len()).then_some(items)
}

/// Reads a shard as `json.load(opener(path)(path, "rt", encoding="utf-8"))` reads it and
/// costs every record. A shard that cannot be read as a list of records with a `start` and
/// an `end` is refused in the tool's own words.
pub fn shard_transcripts(path: &str) -> Result<Transcripts, Failure> {
    let unreadable = |reason: &str| Failure::own(format!("{path} is not a cache shard ({reason})"));
    let bytes = ByteSource::open(path)
        .and_then(ByteSource::read_to_end)
        .map_err(|error| unreadable(error.reason()))?;
    let text = String::from_utf8(bytes).map_err(|error| {
        unreadable(&utf8_decode_error(error.as_bytes()).expect("bytes outside UTF-8"))
    })?;
    let items = match list_items(&text) {
        Some(items) => items,
        // An empty object and an empty text hold no record, as iterating them yields none.
        None => match loads(&text) {
            Ok(Value::Object(object)) if object.is_empty() => Vec::new(),
            Ok(Value::Str(held)) if held.is_empty() => Vec::new(),
            _ => return Err(unreadable("no JSON list of transcript records")),
        },
    };
    items
        .into_iter()
        .enumerate()
        .map(|(place, item)| {
            let record = loads(item)
                .map_err(|error| unreadable(&format!("record {place} is not JSON: {error}")))?;
            let bounds = coordinate(record.get("start")).zip(coordinate(record.get("end")));
            let (start, end) = bounds.ok_or_else(|| {
                unreadable(&format!(
                    "record {place} has no start and end that are whole numbers of 64 bits"
                ))
            })?;
            Ok((start, end, record_cost(&record)))
        })
        .collect()
}

/// A shard file of one chromosome: the first and last position of its name, and its path.
type Shard = (u64, u64, String);

/// The compressed bytes a row adds to the fixture cache the pruner cuts for it: every
/// transcript of the source cache whose span, extended by the flank, overlaps the row's
/// location (and a bracket breakend's mate position), costed at its own deflate size. A
/// shard is read on first touch; the source stores a transcript in every shard it reaches,
/// so only the shards overlapping the flanked location are read. A plain file under
/// `transcripts/` names no sequence and is passed over.
pub struct CacheCost {
    flank: i64,
    shards: HashMap<String, Vec<Shard>>,
    loaded: RefCell<HashMap<String, Result<Transcripts, Failure>>>,
    memo: RefCell<HashMap<(String, i64, i64), usize>>,
    mate: Regex,
}

/// The first and last position a shard file is named after: `<first>-<last>.json`, with or
/// without `.gz`, in ASCII digits; a number past 64 bits reads as the largest.
fn shard_bounds(name: &str) -> Option<(u64, u64)> {
    let stem = name.strip_suffix(".gz").unwrap_or(name);
    let (first, last) = stem.strip_suffix(".json")?.split_once('-')?;
    let number = |digits: &str| {
        let numeric = !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit());
        numeric.then(|| digits.parse::<u64>().unwrap_or(u64::MAX))
    };
    Some((number(first)?, number(last)?))
}

/// The names in a directory; a directory that cannot be listed is refused in the tool's
/// own words.
fn entries(directory: &str) -> Result<Vec<String>, Failure> {
    let unlisted = |error: &std::io::Error| {
        Failure::own(format!(
            "{directory} is not a readable directory ({})",
            os_error_text(error, directory)
        ))
    };
    let mut names = Vec::new();
    for entry in std::fs::read_dir(directory).map_err(|error| unlisted(&error))? {
        let name = entry.map_err(|error| unlisted(&error))?.file_name();
        let name = name.into_string().map_err(|name| {
            Failure::own(format!(
                "{directory} holds {name:?}, a name that is not UTF-8"
            ))
        })?;
        names.push(name);
    }
    Ok(names)
}

impl CacheCost {
    /// The cost reader of the cache `cache` (`str(Path(cache))`) under `flank` bases: the
    /// shard files under each directory of `transcripts/`, read when a row first reaches
    /// them. A `transcripts/` that cannot be listed, or that holds a directory that cannot
    /// be listed, is refused in the tool's own words.
    pub fn new(cache: &str, flank: i64) -> Result<CacheCost, Failure> {
        let transcripts = path_join(cache, "transcripts");
        let mut shards = HashMap::new();
        for chromosome in entries(&transcripts)? {
            let directory = path_join(&transcripts, &chromosome);
            if !is_dir(&directory) {
                continue;
            }
            let mut files: Vec<Shard> = entries(&directory)?
                .into_iter()
                .filter_map(|name| {
                    let (first, last) = shard_bounds(&name)?;
                    Some((first, last, path_join(&directory, &name)))
                })
                .collect();
            files.sort();
            shards.insert(chromosome, files);
        }
        Ok(CacheCost {
            flank,
            shards,
            loaded: RefCell::new(HashMap::new()),
            memo: RefCell::new(HashMap::new()),
            mate: Regex::new(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]").expect("the mate pattern"),
        })
    }

    /// The directories of `transcripts/` the reader holds, in name order.
    pub fn sequences(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.shards.keys().map(String::as_str).collect();
        names.sort_unstable();
        names
    }

    fn flanked(&self, start: i64, end: i64) -> (i128, i128) {
        let flank = i128::from(self.flank);
        (i128::from(start) - flank, i128::from(end) + flank)
    }

    /// The shards of `chromosome` that overlap `start` to `end` extended by the flank.
    fn reached<'a>(
        &'a self,
        chromosome: &str,
        start: i64,
        end: i64,
    ) -> impl Iterator<Item = &'a Shard> + 'a {
        let (low, high) = self.flanked(start, end);
        let shards = self.shards.get(chromosome).map_or(&[][..], Vec::as_slice);
        shards
            .iter()
            .filter(move |(first, last, _)| i128::from(*first) <= high && i128::from(*last) >= low)
    }

    /// The summed cost of the transcripts of `chromosome` that overlap `start` to `end`
    /// extended by the flank, a transcript counted once in every shard that holds it.
    pub fn span_cost(&self, chromosome: &str, start: i64, end: i64) -> Result<usize, Failure> {
        let key = (chromosome.to_string(), start, end);
        if let Some(cost) = self.memo.borrow().get(&key) {
            return Ok(*cost);
        }
        let (low, high) = self.flanked(start, end);
        let mut total = 0;
        for (_, _, path) in self.reached(chromosome, start, end) {
            let mut loaded = self.loaded.borrow_mut();
            if !loaded.contains_key(path) {
                loaded.insert(path.clone(), shard_transcripts(path));
            }
            let transcripts = loaded[path].as_ref().map_err(Clone::clone)?;
            total += transcripts
                .iter()
                .filter(|(first, last, _)| i128::from(*first) <= high && i128::from(*last) >= low)
                .map(|(_, _, cost)| cost)
                .sum::<usize>();
        }
        self.memo.borrow_mut().insert(key, total);
        Ok(total)
    }

    /// The spans `row_cost` costs for a row, in the order it costs them, each with whether
    /// its cost stands in place of the ones before it: the row's own location, then for
    /// each bracket mate of the allele the whole span to a mate on the row's chromosome,
    /// or the mate's own position on another. `None` for a location `int()` does not read.
    fn spans(&self, row: &Row) -> Option<Vec<(String, i64, i64, bool)>> {
        let (chromosome, coordinates) = row.location.split_once(':').unwrap_or((&row.location, ""));
        let (a, b) = coordinates.split_once('-').unwrap_or((coordinates, ""));
        let first = py_int(a)?;
        let second = if b.is_empty() { first } else { py_int(b)? };
        let (start, end) = (first.min(second), first.max(second));
        let mut spans = vec![(chromosome.to_string(), start, end, true)];
        for found in self.mate.captures_iter(&row.allele) {
            let name = &found[1];
            let mate_chromosome = name.strip_prefix("chr").unwrap_or(name);
            let mate = found[2].parse::<i64>().unwrap_or(i64::MAX);
            if mate_chromosome == chromosome {
                // A same-chromosome mate is annotated over the whole breakend-pair span.
                spans.push((chromosome.to_string(), start.min(mate), end.max(mate), true));
            } else {
                spans.push((mate_chromosome.to_string(), mate, mate, false));
            }
        }
        Some(spans)
    }

    /// The bytes `row` adds to the fixture cache; 0 for a location `int()` does not read.
    pub fn row_cost(&self, row: &Row) -> Result<usize, Failure> {
        let mut total = 0;
        for (chromosome, start, end, replaces) in self.spans(row).unwrap_or_default() {
            let cost = self.span_cost(&chromosome, start, end)?;
            total = if replaces { cost } else { total + cost };
        }
        Ok(total)
    }

    /// Reads, on at most `threads` threads and never more than sixteen, the shards the
    /// costs of `rows` reach and no cost has read: what `row_cost` would read one after
    /// another. A shard that cannot be read is refused when a cost reaches it.
    pub fn read_ahead<'a>(&self, rows: impl IntoIterator<Item = &'a Row>, threads: usize) {
        let mut paths: BTreeSet<&str> = BTreeSet::new();
        for row in rows {
            for (chromosome, start, end, _) in self.spans(row).unwrap_or_default() {
                paths.extend(
                    self.reached(&chromosome, start, end)
                        .map(|shard| shard.2.as_str()),
                );
            }
        }
        let mut loaded = self.loaded.borrow_mut();
        let paths: Vec<&str> = paths
            .into_iter()
            .filter(|path| !loaded.contains_key(*path))
            .collect();
        let next = AtomicUsize::new(0);
        let read = || {
            let mut read = Vec::new();
            loop {
                let index = next.fetch_add(1, Ordering::Relaxed);
                let Some(path) = paths.get(index) else {
                    return read;
                };
                read.push((path.to_string(), shard_transcripts(path)));
            }
        };
        let read: Vec<(String, Result<Transcripts, Failure>)> = std::thread::scope(|scope| {
            let workers: Vec<_> = (0..threads.min(MOST_THREADS).min(paths.len()).max(1))
                .map(|_| scope.spawn(read))
                .collect();
            let joined = workers.into_iter().flat_map(|worker| {
                worker
                    .join()
                    .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
            });
            joined.collect()
        });
        loaded.extend(read);
    }
}

#[cfg(test)]
#[path = "cache_cost_tests.rs"]
mod tests;
