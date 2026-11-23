// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The passes of one run, in their order: the inputs opened, the Locations split into
//! shards, each input extracted, each engine paired with the reference and the rows of the
//! remainder written.
//!
//! The remainder's files are written into the run's directory,
//! `input<n>.<projection>.<side>` and `input<n>.rows.<side>` for engine `n` counted from 1;
//! the counts come back to the caller. For the `score` program the run's sources answer
//! the passes directly (`Direct`), and what the passes found is handed on as the settling
//! step takes it (`extracted`, `counted`).

use super::api::{
    text as decoded, Counted, CountedProjection, Extracted as Settled, ScorerError, SourceOf,
    Sources,
};
use crate::extract::{extract, Extracted, Run};
use crate::input::{lines, scratch_file, Input, Kind};
use crate::order::{canonical_location, cmp_column, strip};
use crate::pair::{self, Counts};
use crate::resolve::{Answer, Registry, Resolver, Table};
use crate::rows;
use crate::rules::{Refusal, Res};
use crate::slots::Slots;
use std::collections::hash_map::{Entry, HashMap};
use std::fs;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

/// The locations sampled per shard and input to place the shards' bounds.
const SAMPLES_PER_SHARD: u64 = 32;
/// The input bytes a shard is sized for where the caller names no number of shards.
const SHARD_INPUT_BYTES: u64 = 128 << 20;
const MOST_SHARDS: u64 = 4096;
/// The digits of a coordinate beyond which CPython's `int()` refuses the text.
const MOST_COORDINATE_DIGITS: usize = 4300;

/// What a run of the passes is asked for beside its inputs.
pub struct Options {
    /// Also write the rows a class can read for the remainder.
    pub rows: bool,
    /// The threads that work at one time: those of each pass, those that take a plain
    /// input's digest beside the passes and the caller's own.
    pub threads: usize,
    /// The shards the Locations are split into; by the inputs' size where `None`.
    pub shards: Option<usize>,
    /// The distinct keys of one input that no source claims, at most.
    pub distinct_keys: usize,
    /// The bytes of an input one thread reads at a time.
    pub block_bytes: u64,
}

impl Options {
    pub fn new(rows: bool) -> Options {
        Options {
            rows,
            threads: std::thread::available_parallelism().map_or(1, usize::from),
            shards: None,
            distinct_keys: 100_000,
            block_bytes: 8 << 20,
        }
    }
}

/// One input as the passes found it.
pub struct InputFacts {
    pub kind: Kind,
    /// The sha256 of the file as stored.
    pub sha256: String,
    pub extracted: Extracted,
}

/// One engine paired with the reference.
pub struct EngineFacts {
    pub table: Table,
    /// The counts of the feature projection and of the allele projection.
    pub counts: [Counts; 2],
    /// The rows written for the reference and for the engine.
    pub rows: [u64; 2],
}

/// What the passes hand over: the inputs, the reference first, and every engine's pairing.
pub struct Passed {
    pub inputs: Vec<InputFacts>,
    pub engines: Vec<EngineFacts>,
}

/// The Location with each coordinate written as its number, for any Location, those
/// `order::canonical_location` leaves alone included: a coordinate counts as a number when
/// it is ASCII digits after an optional `+`, and a Location one of whose coordinates has
/// more digits than CPython's `int()` converts stays as written.
pub fn canonical(written: &[u8]) -> Vec<u8> {
    let location = strip(written);
    let Some(colon) = memchr::memchr(b':', location) else {
        return location.to_vec();
    };
    let mut canonical = location[..=colon].to_vec();
    for (index, coordinate) in location[colon + 1..]
        .split(|&byte| byte == b'-')
        .enumerate()
    {
        if index > 0 {
            canonical.push(b'-');
        }
        let digits = coordinate.strip_prefix(b"+").unwrap_or(coordinate);
        if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
            canonical.extend_from_slice(coordinate);
            continue;
        }
        if digits.len() > MOST_COORDINATE_DIGITS {
            return location.to_vec();
        }
        let significant = digits.iter().position(|&digit| digit != b'0');
        canonical.extend_from_slice(significant.map_or(b"0", |at| &digits[at..]));
    }
    canonical
}

/// The sources of a scored run as the passes' resolver.
pub struct Direct<'a, 'r> {
    pub sources: &'a Sources<'r>,
}

impl Resolver for Direct<'_, '_> {
    fn field(&self, name: &[u8]) -> Res<Answer> {
        match self.sources.of(name) {
            Ok(SourceOf::Core) => Ok(Answer::Core),
            Ok(SourceOf::Unclaimed) => Ok(Answer::Unclaimed),
            Ok(SourceOf::Named(source)) => Ok(Answer::Source(source)),
            Err(ScorerError(message)) => Err(Refusal(message)),
        }
    }

    fn location(&self, written: &[u8]) -> Res<Vec<u8>> {
        Ok(canonical(written))
    }
}

/// The keys no source claims as text, in their bytes' order, each with the times its input
/// writes it. Keys of different bytes can decode to the same text, where their counts add.
fn listed_keys(unclaimed: &[(Vec<u8>, u64)]) -> Vec<(String, u64)> {
    let mut listed: Vec<(String, u64)> = Vec::new();
    let mut places: HashMap<String, usize> = HashMap::with_capacity(unclaimed.len());
    for (name, times) in unclaimed {
        match places.entry(decoded(name)) {
            Entry::Occupied(place) => listed[*place.get()].1 += times,
            Entry::Vacant(free) => {
                listed.push((free.key().clone(), *times));
                free.insert(listed.len() - 1);
            }
        }
    }
    listed
}

/// One input as the settling step takes it, with the files of its `side` of the remainder
/// against the engine that is input `engine`.
pub fn extracted(
    ident: &str,
    path: &str,
    facts: &InputFacts,
    out: &Path,
    engine: usize,
    side: &str,
) -> Settled {
    let file = |part: &str| out.join(format!("input{engine}.{part}.{side}"));
    Settled {
        ident: ident.to_string(),
        path: PathBuf::from(path),
        kind: facts.kind,
        sha256: facts.sha256.clone(),
        lines: facts.extracted.lines,
        outside: facts.extracted.outside,
        unclaimed: listed_keys(&facts.extracted.unclaimed),
        header: decoded(&facts.extracted.header),
        feature: file("feature"),
        allele: file("allele"),
        rows: file("rows"),
    }
}

/// What the passes counted of one engine, per projection: each field's values on the two
/// sides with the pairs equal as text, and per source the keys a side shares.
pub fn counted(facts: &EngineFacts) -> Counted {
    let projection = |counts: &Counts| {
        let mut values = Vec::new();
        for (field, name) in facts.table.names.iter().enumerate() {
            let (reference, engine) = (counts.values[0][field], counts.values[1][field]);
            if reference > 0 || engine > 0 {
                values.push((name.clone(), [reference, engine, counts.paired[field]]));
            }
        }
        let shared_keys = facts
            .table
            .sources
            .iter()
            .zip(&counts.shared)
            .filter(|(_, &keys)| keys > 0)
            .map(|(source, &keys)| (source.clone(), keys))
            .collect();
        CountedProjection {
            values,
            shared_keys,
        }
    };
    Counted {
        feature: projection(&facts.counts[0]),
        allele: projection(&facts.counts[1]),
    }
}

/// The Locations that start each shard but the first: the inputs' Locations at evenly
/// spaced offsets, in the lines' order, cut into `shards` runs. Any bounds give the same
/// files; these keep the shards near one size.
fn shard_bounds(inputs: &[Input], shards: usize, registry: &Registry) -> Res<Vec<Vec<u8>>> {
    let mut sampled: Vec<Vec<u8>> = Vec::new();
    let mut window = vec![0u8; 1 << 16];
    let points = shards as u64 * SAMPLES_PER_SHARD;
    for input in inputs {
        for point in 0..if shards > 1 { points } else { 0 } {
            let offset =
                (u128::from(input.length()) * u128::from(point) / u128::from(points)) as u64;
            let count = input.file.read_at(&mut window, offset)?;
            let mut read = &window[..count];
            if offset > 0 {
                // The window opens inside a line, whose end starts the first whole one.
                let Some(end) = memchr::memchr(b'\n', read) else {
                    continue;
                };
                read = &read[end + 1..];
            }
            let Some(row) = lines(read).next() else {
                continue;
            };
            let column = match input.kind {
                Kind::Default if row.starts_with(b"#") => continue,
                Kind::Default => row.split(|&byte| byte == b'\t').nth(1),
                Kind::Tuples => row.split(|&byte| byte == b'\t').next(),
            };
            let Some(written) = column.filter(|_| !row.is_empty()) else {
                continue;
            };
            sampled.push(match canonical_location(written) {
                Some(location) => location,
                None => registry.location(written)?,
            });
        }
    }
    sampled.sort_by(|a, b| cmp_column(a, b));
    sampled.dedup();
    let count = sampled.len();
    let mut bounds: Vec<Vec<u8>> = (1..shards)
        .map(|shard| sampled.get(shard * count / shards).cloned())
        .take_while(Option::is_some)
        .flatten()
        .collect();
    bounds.dedup();
    Ok(bounds)
}

/// Raises its flag when it is dropped.
struct Raised<'a>(&'a AtomicBool);

impl Drop for Raised<'_> {
    fn drop(&mut self) {
        self.0.store(true, Ordering::Relaxed);
    }
}

fn passes(
    reference: &Path,
    engines: &[PathBuf],
    out: &Path,
    scratch: &Path,
    options: &Options,
    registry: &Registry,
    told: &dyn Fn(&str),
) -> Res<Passed> {
    let clock = Instant::now();
    let timed = std::env::var_os("VEP_PLUGIN_PAIRS_TIMINGS").is_some();
    let note = |what: &str| {
        if timed {
            told(&format!(
                "vep-plugin-pairs: {what} at {:.1} s",
                clock.elapsed().as_secs_f64()
            ));
        }
    };
    let mut inputs = Vec::new();
    for (number, path) in std::iter::once(reference)
        .chain(engines.iter().map(PathBuf::as_path))
        .enumerate()
    {
        let plain = scratch_file(scratch, number, "plain");
        inputs.push(Input::open(path, number > 0, &plain, options.block_bytes)?);
    }
    let slots = Slots::new(options.threads);
    let digests_end = AtomicBool::new(false);
    let (extracted, paired_engines) = std::thread::scope(|scope| -> Res<_> {
        let mut own = slots.take();
        // Raised on every way out, so that no digest is read on behind a refusal.
        let _raised = Raised(&digests_end);
        let digests: Vec<_> = inputs
            .iter()
            .filter(|input| input.digest_is_due())
            .map(|input| {
                let (slots, end) = (&slots, &digests_end);
                scope.spawn(move || {
                    let _slot = slots.take();
                    input.digest(end);
                })
            })
            .collect();
        let longest = inputs.iter().map(Input::length).max().unwrap_or(0);
        let shards = options
            .shards
            .unwrap_or((longest / SHARD_INPUT_BYTES + 1).min(MOST_SHARDS) as usize);
        let splitters = shard_bounds(&inputs, shards, registry)?;
        let run = Run {
            registry,
            splitters: &splitters,
            scratch,
            threads: options.threads,
            slots: &slots,
            distinct_keys: options.distinct_keys,
        };
        note(&format!(
            "opened {} inputs, {} shards",
            inputs.len(),
            run.shards()
        ));
        let first = extract(&run, &inputs[0], 0, &mut own)?;
        note("extracted the reference");
        if first.lines == 0 {
            return Err(Refusal(format!(
                "{}: the reference has no data row",
                inputs[0].shown
            )));
        }
        let mut extracted = Vec::new();
        let mut paired_engines = Vec::new();
        for number in 1..inputs.len() {
            let engine = extract(&run, &inputs[number], number, &mut own)?;
            note(&format!("extracted engine {number}"));
            let table = registry.table();
            let kind = inputs[number].kind;
            let paired = pair::pair(&run, &table, number, kind, options.rows, out, &mut own)?;
            note(&format!("paired engine {number}"));
            let mut written = [0, paired.tuple_rows];
            if options.rows {
                let wanted = &paired.locations;
                let path = scratch_file(out, number, "rows.reference");
                written[0] = rows::write(&run, &inputs[0], &first, wanted, &path, &mut own)?;
                if kind == Kind::Default {
                    let path = scratch_file(out, number, "rows.engine");
                    written[1] =
                        rows::write(&run, &inputs[number], &engine, wanted, &path, &mut own)?;
                }
                note(&format!("wrote the rows of engine {number}"));
            }
            for shard in 0..run.shards() {
                let _ = fs::remove_file(scratch_file(scratch, number, &format!("shard{shard}")));
            }
            extracted.push(engine);
            paired_engines.push(EngineFacts {
                table,
                counts: paired.counts,
                rows: written,
            });
        }
        let lent = own.lent();
        for digest in digests {
            digest.join().expect("a digest thread panicked");
        }
        drop(lent);
        note("digested the inputs");
        extracted.insert(0, first);
        Ok((extracted, paired_engines))
    })?;
    let mut facts = Vec::new();
    for (input, extracted) in inputs.iter().zip(extracted) {
        facts.push(InputFacts {
            kind: input.kind,
            sha256: input.sha256()?.to_string(),
            extracted,
        });
    }
    Ok(Passed {
        inputs: facts,
        engines: paired_engines,
    })
}

/// Runs the passes over `reference` and each of `engines`, asking `registry` for the source
/// of every field name, and writes the remainder's files into `out`, which it creates. The
/// passes' own scratch, `out/shards`, is removed before the return. With
/// `VEP_PLUGIN_PAIRS_TIMINGS` in the environment the time each pass ended is `told`, a line
/// a call.
pub fn run(
    reference: &Path,
    engines: &[PathBuf],
    out: &Path,
    options: &Options,
    registry: &Registry,
    told: &dyn Fn(&str),
) -> Res<Passed> {
    let scratch = out.join("shards");
    fs::create_dir_all(&scratch)?;
    let outcome = passes(reference, engines, out, &scratch, options, registry, told);
    let _ = fs::remove_dir_all(&scratch);
    outcome
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settled(written: &str) -> String {
        String::from_utf8_lossy(&canonical(written.as_bytes())).into_owned()
    }

    #[test]
    fn location_outside_ascii_keeps_its_bytes_and_numbers_its_ascii_coordinates() {
        assert_eq!(settled("21:²"), "21:²");
        assert_eq!(settled("21:٣"), "21:٣");
        assert_eq!(settled("chré:007-+08"), "chré:7-8");
        assert_eq!(settled(" 21:000 "), "21:0");
        assert_eq!(settled("no_coordinates"), "no_coordinates");
    }

    #[test]
    fn unclaimed_keys_of_one_text_add_their_counts_at_the_place_of_the_first() {
        let keys = |keys: &[(&[u8], u64)]| -> Vec<(Vec<u8>, u64)> {
            keys.iter()
                .map(|(name, times)| (name.to_vec(), *times))
                .collect()
        };
        let listed = listed_keys(&keys(&[
            (b"A", 1),
            (b"\xef\xbf\xbdA", 2),
            (b"\xf0\x9f\x98\x80", 3),
            (b"\xfeA", 4),
            (b"\xffA", 5),
        ]));
        assert_eq!(
            listed,
            [
                ("A".to_string(), 1),
                ("\u{fffd}A".to_string(), 11),
                ("\u{1f600}".to_string(), 3),
            ]
        );
        assert_eq!(listed_keys(&[]), []);
    }

    #[test]
    fn location_with_a_coordinate_beyond_the_integer_conversion_stays_as_written() {
        let long = format!("21:{}", "9".repeat(300));
        assert_eq!(settled(&long), long);
        let zeros = format!("21:{}7", "0".repeat(4299));
        assert_eq!(settled(&zeros), "21:7");
        let refused = format!(" 21:05-{}7 ", "0".repeat(4300));
        assert_eq!(settled(&refused), refused.trim());
    }
}
