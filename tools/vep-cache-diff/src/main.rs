// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Structural comparison of two JSON cache directories, a baseline and a builder's.
//!
//! The transcript shards of the two caches are paired on their path under `transcripts/`
//! and compared record by record on the stable id; the variation files are paired on the
//! chromosome and joined on the variation name; the two `info.json` are compared key by
//! key. The shards and the chromosomes are compared across the host's cores, and the
//! report is one JSON document with the counts, the match rates and the mismatches by
//! field.

mod cli;
mod info;
mod input;
mod json;
mod paths;
mod pyfmt;
mod report;
mod threads;
mod transcripts;
mod unicode;
mod variations;

use json::{Object, Value};
use paths::Pair;
use report::Json;
use std::ffi::OsString;
use std::fs;
use std::io::{self, Write};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;
use std::process::ExitCode;

/// How many shard names with differences the report lists.
const REGIONS_LISTED: usize = 50;

/// The counts both comparisons report, as the members that open their blocks.
fn count_members(names: [u64; 2], shared: u64, matched: u64, differed: u64) -> Vec<(String, Json)> {
    [
        ("total_baseline", names[0]),
        ("total_builder", names[1]),
        ("shared", shared),
        ("only_baseline", names[0] - shared),
        ("only_builder", names[1] - shared),
        ("matched", matched),
        ("differed", differed),
    ]
    .into_iter()
    .map(|(name, count)| (name.to_string(), Json::Count(count)))
    .collect()
}

/// `dict(Counter.most_common())`: the names by count, the largest first, equal counts in
/// the order of `met`.
fn most_common(mut met: Vec<(String, u64)>) -> Json {
    met.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
    Json::Object(
        met.into_iter()
            .map(|(name, count)| (name, Json::Count(count)))
            .collect(),
    )
}

/// One comparison's block of the report and its stdout line.
struct Block {
    report: Json,
    line: String,
}

fn ids(list: &[Box<[u8]>]) -> Json {
    Json::List(
        list.iter()
            .map(|id| Json::Text(json::code_points(id)))
            .collect(),
    )
}

/// Compares the transcript shards of the two caches.
fn diff_transcripts(baseline: &Path, builder: &Path, threads: usize) -> Result<Block, String> {
    let pairs = paths::shard_pairs(baseline, builder);
    let regions = threads::run_indexed(threads, pairs.len(), |i| {
        transcripts::compare_region(&pairs[i])
    })?;
    let mut names = [0u64; 2];
    let (mut shared, mut matched, mut differed) = (0u64, 0u64, 0u64);
    let mut mismatches = [0u64; transcripts::MISMATCH_NAMES];
    // The mismatch names in the order the shards, their ids and a record's comparison
    // first meet them.
    let mut met: Vec<usize> = Vec::new();
    let mut listed = Vec::new();
    for (pair, region) in pairs.iter().zip(&regions) {
        names[0] += region.baseline_count as u64;
        names[1] += region.builder_count as u64;
        shared += region.shared as u64;
        matched += region.matched as u64;
        differed += region.differed as u64;
        let mut first: Vec<(usize, usize)> = (0..transcripts::MISMATCH_NAMES)
            .filter(|name| mismatches[*name] == 0)
            .filter_map(|name| region.first_mismatch[name].map(|id| (id, name)))
            .collect();
        first.sort_unstable();
        met.extend(first.into_iter().map(|(_, name)| name));
        for (total, count) in mismatches.iter_mut().zip(region.mismatches) {
            *total += count as u64;
        }
        if region.has_differences() && listed.len() < REGIONS_LISTED {
            listed.push(Json::object([
                ("region", Json::Text(pyfmt::surrogateescape(&pair.name))),
                ("baseline_count", Json::count(region.baseline_count)),
                ("builder_count", Json::count(region.builder_count)),
                ("shared", Json::count(region.shared)),
                ("only_baseline", ids(&region.first_only_baseline)),
                ("only_builder", ids(&region.first_only_builder)),
                ("matched", Json::count(region.matched)),
                ("differed", Json::count(region.differed)),
            ]));
        }
    }
    let rate = report::match_rate(matched, shared);
    let mut members = count_members(names, shared, matched, differed);
    let met = met
        .into_iter()
        .map(|name| (transcripts::mismatch_name(name), mismatches[name]))
        .collect();
    members.push(("field_mismatches".to_string(), most_common(met)));
    members.push(("match_rate".to_string(), Json::Float(rate)));
    members.push(("regions_with_diffs".to_string(), Json::List(listed)));
    Ok(Block {
        report: Json::Object(members),
        line: report::summary_line(shared, matched, differed, rate),
    })
}

/// A cache's `info.json`: `None` where the cache has none, its object otherwise. A
/// document that is no object is kept as it is, for the comparison that reads it to
/// refuse.
fn read_info(cache: &Path) -> Result<Option<Value>, String> {
    let path = cache.join("info.json");
    if fs::metadata(&path).is_err() {
        return Ok(None);
    }
    let named = |detail: String| format!("{}: {detail}", path.display());
    let text = fs::read(&path).map_err(|e| named(e.to_string()))?;
    json::load(&text).map(Some).map_err(named)
}

fn as_object(info: &Option<Value>) -> Option<&Object> {
    match info {
        Some(Value::Object(object)) => Some(object),
        _ => None,
    }
}

/// Compares the variation files of the two caches chromosome by chromosome.
fn diff_variations(
    baseline: &Path,
    builder: &Path,
    infos: &[Option<Value>; 2],
    (threads, budget): (usize, usize),
) -> Result<Block, String> {
    let columns = |side: usize, cache: &Path| {
        variations::columns(infos[side].as_ref())
            .map_err(|e| format!("{}: {e}", cache.join("info.json").display()))
    };
    let baseline_columns = columns(0, baseline)?;
    let builder_columns = columns(1, builder)?;
    let frequency = variations::frequency_columns(&baseline_columns, &builder_columns);
    let layouts = [
        variations::Layout::new(&baseline_columns, &frequency),
        variations::Layout::new(&builder_columns, &frequency),
    ];
    let pairs = paths::chromosome_pairs(baseline, builder);
    // The largest chromosomes start first, so that the longest of them is not the last.
    let size = |pair: &Pair| -> u64 {
        [&pair.baseline, &pair.builder]
            .into_iter()
            .flatten()
            .map(|path| fs::metadata(path).map_or(0, |meta| meta.len()))
            .sum()
    };
    let mut order: Vec<(u64, usize)> = pairs.iter().map(size).zip(0..).collect();
    order.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    let spills = variations::Spills::default();
    let compared = threads::run_indexed(threads, pairs.len(), |turn| {
        let chromosome = order[turn].1;
        variations::compare_chromosome(&pairs[chromosome], &layouts, budget, &spills, chromosome)
    })?;
    drop(spills);
    let mut by_chromosome: Vec<Option<variations::Counts>> = vec![None; pairs.len()];
    for (turn, counts) in compared.into_iter().enumerate() {
        by_chromosome[order[turn].1] = Some(counts);
    }
    let mut total = variations::Counts::default();
    // The mismatch names in the order the chromosomes first meet them; within one
    // chromosome, in the order a record's comparison meets them.
    let mut met: Vec<usize> = Vec::new();
    for counts in by_chromosome.into_iter().flatten() {
        for side in 0..2 {
            total.names[side] += counts.names[side];
        }
        total.shared += counts.shared;
        total.matched += counts.matched;
        total.differed += counts.differed;
        for name in 0..variations::MISMATCH_NAMES {
            if total.mismatches[name] == 0 && counts.mismatches[name] > 0 {
                met.push(name);
            }
            total.mismatches[name] += counts.mismatches[name];
        }
    }
    let rate = report::match_rate(total.matched, total.shared);
    let mut members = count_members(total.names, total.shared, total.matched, total.differed);
    let met = met
        .into_iter()
        .map(|name| {
            (
                variations::mismatch_name(name).to_string(),
                total.mismatches[name],
            )
        })
        .collect();
    members.push(("field_mismatches".to_string(), most_common(met)));
    let frequency = frequency
        .iter()
        .map(|name| Json::Text(json::code_points(name)))
        .collect();
    members.push((
        "frequency_columns_compared".to_string(),
        Json::List(frequency),
    ));
    members.push(("match_rate".to_string(), Json::Float(rate)));
    Ok(Block {
        report: Json::Object(members),
        line: report::summary_line(total.shared, total.matched, total.differed, rate),
    })
}

/// Writes one line on stdout at once, so that a run that fails has written every line
/// before the failure. A path behind `line` is written as its bytes.
fn say(line: &str, path: Option<&Path>) -> Result<(), String> {
    let mut bytes = line.as_bytes().to_vec();
    if let Some(path) = path {
        bytes.extend_from_slice(path.as_os_str().as_bytes());
    }
    bytes.push(b'\n');
    let mut stdout = io::stdout().lock();
    stdout
        .write_all(&bytes)
        .and_then(|()| stdout.flush())
        .map_err(|e| format!("stdout: {e}"))
}

/// A path as the report holds it: a byte outside UTF-8 as the lone surrogate Python
/// decodes it to.
fn in_report(path: &Path) -> Json {
    Json::Text(pyfmt::surrogateescape(path.as_os_str().as_bytes()))
}

fn write_report(path: &Path, report: &Json) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut text = report::dumps(report);
    text.push('\n');
    fs::write(path, text)
}

/// Runs the comparison: each stdout line as its step ends, then the report.
fn run(baseline: &Path, builder: &Path, output: &Path) -> Result<(), String> {
    let threads = threads::thread_count()?;
    let budget = threads::buffer_bytes()?;
    say("Comparing caches:", None)?;
    say("  Baseline: ", Some(baseline))?;
    say("  Builder:  ", Some(builder))?;
    say("  Comparing transcripts...", None)?;
    let transcripts = diff_transcripts(baseline, builder, threads)?;
    say(&transcripts.line, None)?;
    say("  Comparing variations...", None)?;
    let infos = [read_info(baseline)?, read_info(builder)?];
    let variations = diff_variations(baseline, builder, &infos, (threads, budget))?;
    say(&variations.line, None)?;
    say("  Comparing info.json...", None)?;
    let info_json = info::compare(as_object(&infos[0]), as_object(&infos[1]));
    let matches = matches!(&info_json, Json::Object(members) if members[0].1 == Json::Bool(true));
    say(
        if matches {
            "    matches=True"
        } else {
            "    matches=False"
        },
        None,
    )?;
    let report = Json::object([
        ("baseline_dir", in_report(baseline)),
        ("builder_dir", in_report(builder)),
        ("transcripts", transcripts.report),
        ("variations", variations.report),
        ("info_json", info_json),
    ]);
    write_report(output, &report).map_err(|e| format!("{}: {e}", output.display()))?;
    say("\nReport written to: ", Some(output))
}

/// A path as Python writes it on stderr: a byte outside UTF-8 as `\udcXX`.
fn on_stderr(path: &Path) -> String {
    pyfmt::surrogateescape(path.as_os_str().as_bytes())
        .into_iter()
        .map(|unit| match char::from_u32(unit) {
            Some(c) => c.to_string(),
            None => format!("\\u{unit:04x}"),
        })
        .collect()
}

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().skip(1).collect();
    let args = match cli::parse(&argv) {
        cli::Parsed::Run(args) => args,
        cli::Parsed::Help => {
            println!("{}\n\n{}", cli::USAGE, cli::HELP);
            return ExitCode::SUCCESS;
        }
        cli::Parsed::Usage(e) => {
            eprintln!("{}\nvep-cache-diff: error: {e}", cli::USAGE);
            return ExitCode::from(2);
        }
    };
    let resolved = [&args.baseline_dir, &args.builder_dir, &args.output_report]
        .map(|arg| paths::py_resolve(arg));
    let [baseline, builder, output] = match resolved {
        [Ok(baseline), Ok(builder), Ok(output)] => [baseline, builder, output],
        _ => {
            eprintln!("ERROR: [vep-cache-diff] the working directory cannot be read");
            return ExitCode::FAILURE;
        }
    };
    for (side, directory) in [("baseline", &baseline), ("builder", &builder)] {
        if !directory.is_dir() {
            eprintln!(
                "ERROR: [diff_caches] {side} directory not found: {}",
                on_stderr(directory)
            );
            return ExitCode::FAILURE;
        }
    }
    match run(&baseline, &builder, &output) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("ERROR: [vep-cache-diff] {e}");
            ExitCode::FAILURE
        }
    }
}
