// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The heavy path of `scripts/concordance/score_plugin_tuples.py` on every core.
//!
//! The scorer drives this binary and stays the one place the plugin registry, the number
//! rule, the documented classes, the residual files and the report are defined. The binary
//! reads the reference and each engine's output, extracts the populated plugin fields by
//! the scorer's extraction rules, asking the scorer over standard output and input for the
//! source of each field name, splits the values by Location into shards, sorts and pairs
//! the values equal as text, and hands over the counts (`pairs.json`) and everything it does
//! not pair, the remainder, as sorted tuple lines beside the rows a class can read for them.

mod extract;
mod handover;
mod input;
mod order;
mod pair;
mod resolve;
mod rows;

use extract::{extract, Run};
use handover::{pairs_json, EngineFacts, InputFacts};
use input::{lines, scratch_file, Input, Kind};
use order::{canonical_location, cmp_column};
use resolve::Registry;
use std::ffi::OsString;
use std::fs;
use std::io;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

const USAGE: &str = "usage: vep-plugin-pairs --reference PATH --engine PATH [--engine PATH ...] --out DIR [--rows] [--threads N] [--shards N] [--distinct-keys N] [--block-bytes N]";
const HELP: &str = "\
Pairs the plugin values of a reference output with each engine's for score_plugin_tuples.py.

options:
  --reference PATH   Ensembl VEP default-format output, plain or gzip
  --engine PATH      a scored engine's default-format output or tuple file; repeat per engine
  --out DIR          directory for pairs.json, the remainder's files and the scratch
  --rows             also write the rows a class can read for the remainder
  --threads N        threads of each pass (default: every core)
  --shards N         shards the Locations are split into (default: by the inputs' size)
  --distinct-keys N  distinct keys of one input that no source claims, at most (default 100000)
  --block-bytes N    bytes of an input read at a time by one thread (default 8388608)

The source of every field name is asked on standard output and answered on standard input.
VEP_PLUGIN_PAIRS_TIMINGS in the environment prints the time of each pass on standard error.";
/// The locations sampled per shard and input to place the shards' bounds.
const SAMPLES_PER_SHARD: u64 = 32;
/// The input bytes a shard is sized for where `--shards` is not given.
const SHARD_INPUT_BYTES: u64 = 128 << 20;
const MOST_SHARDS: u64 = 4096;

/// A refusal, printed as one `ERROR: [vep-plugin-pairs]` line.
pub struct Refusal(pub String);

pub type Res<T> = Result<T, Refusal>;

impl From<io::Error> for Refusal {
    fn from(error: io::Error) -> Refusal {
        Refusal(error.to_string())
    }
}

struct Args {
    reference: PathBuf,
    engines: Vec<PathBuf>,
    out: PathBuf,
    rows: bool,
    threads: usize,
    shards: Option<usize>,
    distinct_keys: usize,
    block_bytes: u64,
}

/// The run's options, or `None` for the help: `-h` or `--help` where an option's name stands,
/// since after an option that takes a path either is the path.
fn parse_args(argv: Vec<OsString>) -> Result<Option<Args>, String> {
    let mut reference = None;
    let mut engines = Vec::new();
    let mut out = None;
    let mut rows = false;
    let mut threads = std::thread::available_parallelism().map_or(1, usize::from);
    let mut shards = None;
    let mut distinct_keys = 100_000;
    let mut block_bytes = 8 << 20;
    let mut rest = argv.into_iter();
    while let Some(argument) = rest.next() {
        let name = argument.to_string_lossy().into_owned();
        let mut value = || {
            rest.next()
                .ok_or_else(|| format!("argument {name}: expected one argument"))
        };
        let mut number = |least: u64| {
            let text = value()?.to_string_lossy().into_owned();
            match text.parse::<u64>() {
                Ok(number) if number >= least => Ok(number),
                _ => Err(format!("argument {name}: invalid value: '{text}'")),
            }
        };
        match name.as_str() {
            "-h" | "--help" => return Ok(None),
            "--reference" => reference = Some(PathBuf::from(value()?)),
            "--engine" => engines.push(PathBuf::from(value()?)),
            "--out" => out = Some(PathBuf::from(value()?)),
            "--rows" => rows = true,
            "--threads" => threads = number(1)? as usize,
            "--shards" => shards = Some(number(1)? as usize),
            "--distinct-keys" => distinct_keys = number(0)? as usize,
            "--block-bytes" => block_bytes = number(1)?,
            _ => return Err(format!("unrecognized arguments: {name}")),
        }
    }
    let required = |name: &str| format!("the following arguments are required: {name}");
    if engines.is_empty() {
        return Err(required("--engine"));
    }
    Ok(Some(Args {
        reference: reference.ok_or_else(|| required("--reference"))?,
        engines,
        out: out.ok_or_else(|| required("--out"))?,
        rows,
        threads,
        shards,
        distinct_keys,
        block_bytes,
    }))
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

fn pairs(args: &Args, scratch: &Path) -> Res<()> {
    let clock = Instant::now();
    let timed = std::env::var_os("VEP_PLUGIN_PAIRS_TIMINGS").is_some();
    let note = |what: &str| {
        if timed {
            eprintln!(
                "vep-plugin-pairs: {what} at {:.1} s",
                clock.elapsed().as_secs_f64()
            );
        }
    };
    let registry = Registry::over_standard_streams();
    let mut inputs = Vec::new();
    for (number, path) in std::iter::once(&args.reference)
        .chain(&args.engines)
        .enumerate()
    {
        let plain = scratch_file(scratch, number, "plain");
        inputs.push(Input::open(path, number > 0, &plain, args.block_bytes)?);
    }
    let longest = inputs.iter().map(Input::length).max().unwrap_or(0);
    let shards = args
        .shards
        .unwrap_or((longest / SHARD_INPUT_BYTES + 1).min(MOST_SHARDS) as usize);
    let splitters = shard_bounds(&inputs, shards, &registry)?;
    let run = Run {
        registry: &registry,
        splitters: &splitters,
        scratch,
        threads: args.threads,
        distinct_keys: args.distinct_keys,
    };
    note(&format!(
        "opened {} inputs, {} shards",
        inputs.len(),
        run.shards()
    ));
    let reference = extract(&run, &inputs[0], 0)?;
    note("extracted the reference");
    if reference.lines == 0 {
        return Err(Refusal(format!(
            "{}: the reference has no data row",
            inputs[0].shown
        )));
    }
    let mut extracted = Vec::new();
    let mut engines = Vec::new();
    for number in 1..inputs.len() {
        let engine = extract(&run, &inputs[number], number)?;
        note(&format!("extracted engine {number}"));
        let table = registry.table();
        let kind = inputs[number].kind;
        let paired = pair::pair(&run, &table, number, kind, args.rows, &args.out)?;
        note(&format!("paired engine {number}"));
        let mut rows = [0, paired.tuple_rows];
        if args.rows {
            let wanted = &paired.locations;
            let path = scratch_file(&args.out, number, "rows.reference");
            rows[0] = rows::write(&run, &inputs[0], &reference, wanted, &path)?;
            if kind == Kind::Default {
                let path = scratch_file(&args.out, number, "rows.engine");
                rows[1] = rows::write(&run, &inputs[number], &engine, wanted, &path)?;
            }
            note(&format!("wrote the rows of engine {number}"));
        }
        for shard in 0..run.shards() {
            let _ = fs::remove_file(scratch_file(scratch, number, &format!("shard{shard}")));
        }
        extracted.push(engine);
        engines.push(EngineFacts {
            table,
            counts: paired.counts,
            rows,
        });
    }
    let mut facts = Vec::new();
    for (input, extracted) in inputs
        .iter_mut()
        .zip(std::iter::once(reference).chain(extracted))
    {
        facts.push(InputFacts {
            kind: input.kind,
            sha256: input.sha256()?.to_string(),
            extracted,
        });
    }
    note("digested the inputs");
    let staged = args.out.join("pairs.json.tmp");
    fs::write(&staged, pairs_json(&facts, &engines))?;
    fs::rename(&staged, args.out.join("pairs.json"))?;
    Ok(())
}

fn run(args: &Args) -> Res<()> {
    let scratch = args.out.join("shards");
    fs::create_dir_all(&scratch)?;
    let outcome = pairs(args, &scratch);
    let _ = fs::remove_dir_all(&scratch);
    outcome
}

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().skip(1).collect();
    let args = match parse_args(argv) {
        Ok(Some(args)) => args,
        Ok(None) => {
            println!("{USAGE}\n\n{HELP}");
            return ExitCode::SUCCESS;
        }
        Err(message) => {
            eprintln!("{USAGE}\nERROR: [vep-plugin-pairs] {message}");
            return ExitCode::from(2);
        }
    };
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(Refusal(message)) => {
            eprintln!("ERROR: [vep-plugin-pairs] {message}");
            ExitCode::from(2)
        }
    }
}
