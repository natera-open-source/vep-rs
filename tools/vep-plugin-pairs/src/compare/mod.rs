// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `compare` program: the per-field comparison of two outputs' plugin values.
//!
//! One file pair, or every pair of two directories, is compared on the fields the
//! registry lists for the requested plugins. Each output's data lines are read once, a
//! plain file in ranges on the threads and a gzip file as one stream, and what a row gives a
//! plugin is written to the shard of the row's key in the run's scratch. A thread then
//! reads one shard of both sides, sorts it and walks the two sides together, so that a
//! pair's memory is one shard a thread, whatever the outputs' length. What a pair keeps
//! is its counts, the fields it met and its first mismatches.

mod cli;
mod fields;
mod files;
mod read;
mod report;
mod shards;
mod tally;
mod text;

use crate::registry::{all_plugin_names, field_sets};
use crate::score::{interrupt, paths};
use cli::{Command, Options};
use fields::{FieldNames, Run};
use read::Reading;
use report::FileResults;
use shards::{Scratch, Shard};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;
use tally::{Counting, PluginResult, Terms};
use text::{path_parent, path_stem, path_str, py_float_repr, py_strip};

/// Writes to standard output; a reader that went away is not an error of the comparison.
fn say(text: &str) {
    let mut out = std::io::stdout().lock();
    let _ = out.write_all(text.as_bytes());
    let _ = out.flush();
}

/// Writes a line to standard error; a reader that went away takes nothing, and the run's
/// exit status stands.
fn stderr_line(line: &str) {
    let _ = writeln!(std::io::stderr().lock(), "{line}");
}

/// A failure of the run itself, as one line on standard error.
fn fail(message: &str) -> i32 {
    stderr_line(&format!("ERROR: [vep-plugin-pairs] {message}"));
    1
}

/// The run's scratch directory, removed with everything in it when dropped.
struct Removed(PathBuf);

impl Drop for Removed {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The length a pipe is taken to have where the shards are counted, and the times a gzip
/// file's length is taken to grow when unpacked.
const PIPE_BYTES: u64 = 32 << 30;
const GZIP_GROWTH: u64 = 16;

/// The shards of a file pair: one per `shard_bytes` of the longer output, so that a shard
/// holds no more than about that much of either side.
fn shard_count(paths: [&str; 2], shard_bytes: u64) -> usize {
    let length = |path: &str| match std::fs::metadata(path) {
        Ok(found) if !found.is_file() => PIPE_BYTES,
        Ok(found) if read::is_gzip_path(path) => found.len().saturating_mul(GZIP_GROWTH),
        Ok(found) => found.len(),
        Err(_) => 0,
    };
    let longer = length(paths[0]).max(length(paths[1]));
    longer
        .div_ceil(shard_bytes.max(1))
        .clamp(1, shards::MOST_SHARDS) as usize
}

/// One file pair's result per distinct plugin of the run, or why the pair cannot be read.
type Compared = Result<Vec<PluginResult>, String>;

/// Compares one file pair, the run's pair `number`, on at most `threads` threads: the two
/// outputs are read into the pair's shards under `scratch`, both at once where there are
/// two threads, and the shards are compared, a shard a thread.
fn compare_files(
    (perl_path, rust_path): (&str, &str),
    number: usize,
    run: &Run,
    options: &Options,
    threads: usize,
    scratch: &Path,
) -> Compared {
    let clock = Instant::now();
    let timed = std::env::var_os("VEP_PLUGIN_PAIRS_TIMINGS").is_some();
    let note = |what: &str| {
        if timed {
            stderr_line(&format!(
                "vep-plugin-pairs: {perl_path} {what} at {:.1} s",
                clock.elapsed().as_secs_f64()
            ));
        }
    };
    let shards = shard_count([perl_path, rust_path], options.shard_bytes);
    let in_scratch = |error: std::io::Error| format!("{}: {error}", scratch.display());
    let scratch =
        Scratch::new(scratch.join(format!("pair{number}")), shards).map_err(in_scratch)?;
    let names = FieldNames::default();
    let extract = |path: &str, side: usize, threads: usize| {
        let reading = Reading {
            threads,
            block_bytes: options.block_bytes,
        };
        read::extract(path, run, &names, &scratch, side, reading)
            .map_err(|error| format!("{path}: {error}"))
    };
    let (perl, rust) = if threads >= 2 {
        // A gzip file is one stream on one thread, which leaves the others to the other side.
        let perl_threads = match (read::is_gzip_path(perl_path), read::is_gzip_path(rust_path)) {
            (true, _) => 1,
            (false, true) => threads - 1,
            (false, false) => threads - threads / 2,
        };
        std::thread::scope(|scope| {
            let rust = scope.spawn(|| extract(rust_path, 1, (threads - perl_threads).max(1)));
            let perl = extract(perl_path, 0, perl_threads);
            (perl, rust.join().expect("the reader does not panic"))
        })
    } else {
        (extract(perl_path, 0, 1), extract(rust_path, 1, 1))
    };
    perl?;
    rust?;
    if timed {
        note(&format!(
            "read into {shards} shards of {} bytes",
            scratch.bytes()
        ));
    }

    let names = names.names();
    let terms = Terms {
        names: &names,
        tolerance: options.tolerance,
        max_mismatches: options.max_mismatches,
    };
    let next = AtomicUsize::new(0);
    let tallies: Mutex<Vec<Vec<Counting>>> = Mutex::new(Vec::new());
    let failed: Mutex<Option<std::io::Error>> = Mutex::new(None);
    let work = || loop {
        let shard = next.fetch_add(1, Ordering::Relaxed);
        if shard >= shards
            || interrupt::requested()
            || failed.lock().expect("no worker panics").is_some()
        {
            break;
        }
        let sides = Shard::load(&scratch, 0, shard)
            .and_then(|perl| Ok((perl, Shard::load(&scratch, 1, shard)?)));
        match sides {
            Ok((perl, rust)) => {
                let tally = tally::compare_shard(&perl, &rust, run, &terms);
                tallies.lock().expect("no worker panics").push(tally);
            }
            Err(error) => *failed.lock().expect("no worker panics") = Some(error),
        }
    };
    std::thread::scope(|scope| {
        for _ in 1..threads.min(shards) {
            scope.spawn(work);
        }
        work();
    });
    if let Some(error) = failed.into_inner().expect("no worker panics") {
        return Err(format!("{}: {error}", scratch.directory().display()));
    }
    if interrupt::requested() {
        return Err("interrupted".to_string());
    }
    let tallies = tallies.into_inner().expect("no worker panics");
    let results = tally::join(tallies, run, &terms);
    note("compared");
    Ok(results)
}

/// Every file pair's results in the pairs' order, the pairs compared on as many threads as
/// the run allows. The first pair in order that cannot be read is the run's failure.
fn compare_all(
    pairs: &[(String, String, String)],
    run: &Run,
    options: &Options,
    scratch: &Path,
) -> Result<Vec<FileResults>, String> {
    let workers = options.threads.min(pairs.len()).max(1);
    let threads_of_a_pair = (options.threads / workers).max(1);
    let next = AtomicUsize::new(0);
    let done: Vec<Mutex<Option<Compared>>> = pairs.iter().map(|_| Mutex::new(None)).collect();
    std::thread::scope(|scope| {
        for _ in 0..workers {
            scope.spawn(|| loop {
                let at = next.fetch_add(1, Ordering::Relaxed);
                let Some((perl_path, rust_path, _)) = pairs.get(at) else {
                    break;
                };
                let paths = (perl_path.as_str(), rust_path.as_str());
                let result = compare_files(paths, at, run, options, threads_of_a_pair, scratch);
                *done[at].lock().expect("no worker panics") = Some(result);
            });
        }
    });
    pairs
        .iter()
        .zip(done)
        .map(|((_, _, name), result)| {
            let plugins = result
                .into_inner()
                .expect("no worker panics")
                .expect("every pair is compared")?;
            Ok(FileResults {
                name: name.clone(),
                plugins,
            })
        })
        .collect()
}

/// Runs `compare` over the command line after the word and returns the exit status: 0
/// where every plugin passes, 1 where one fails or the run cannot be made, 2 for a wrong
/// command line or a plugin the registry does not know, 130 where a stop signal ended it.
pub fn main(args: Vec<String>) -> i32 {
    let options = match cli::parse(&args) {
        Command::Help => {
            say(&format!("{}\n\n{}\n", cli::USAGE, cli::HELP));
            return 0;
        }
        Command::Refused(message) => {
            stderr_line(&format!(
                "{}\nvep-plugin-pairs compare: error: {message}",
                cli::USAGE
            ));
            return 2;
        }
        Command::Run(options) => options,
    };

    let plugins: Vec<String> = if options.plugins.eq_ignore_ascii_case("all") {
        all_plugin_names().to_vec()
    } else {
        let named: Vec<String> = options
            .plugins
            .split(',')
            .map(|name| py_strip(name).to_string())
            .collect();
        let known = |name: &String| field_sets().iter().any(|set| &set.name == name);
        let unknown: Vec<&str> = named
            .iter()
            .filter(|name| !known(name))
            .map(String::as_str)
            .collect();
        if !unknown.is_empty() {
            // A misspelt name would otherwise be compared on no field at all.
            let mut sets: Vec<&str> = field_sets().iter().map(|set| set.name.as_str()).collect();
            sets.sort_unstable();
            stderr_line(&format!(
                "ERROR: unknown plugin(s): {}. Known: {}",
                unknown.join(", "),
                sets.join(", ")
            ));
            return 2;
        }
        named
    };
    let fields = options
        .fields
        .as_deref()
        .filter(|fields| !fields.is_empty())
        .map(|fields| {
            fields
                .split(',')
                .map(|field| py_strip(field).to_string())
                .collect()
        });
    let run = Run::new(plugins, fields);

    let given = |value: &Option<String>| value.clone().filter(|value| !value.is_empty());
    let pairs = match (
        given(&options.perl_dir),
        given(&options.rust_dir),
        given(&options.perl_output),
        given(&options.rust_output),
    ) {
        (Some(perl_dir), Some(rust_dir), _, _) => {
            let pairs = files::matching_pairs(&perl_dir, &rust_dir);
            if pairs.is_empty() {
                stderr_line(&format!(
                    "ERROR: no matching file pairs found in {perl_dir} and {rust_dir}"
                ));
                return 1;
            }
            pairs
        }
        (_, _, Some(perl_output), Some(rust_output)) => {
            let name = path_stem(&perl_output).to_string();
            vec![(perl_output, rust_output, name)]
        }
        _ => {
            stderr_line(
                "ERROR: provide either --perl-dir/--rust-dir or --perl-output/--rust-output",
            );
            return 1;
        }
    };

    say(&format!(
        "Plugins: {}\nFile pairs: {}\nTolerance: {}\n",
        run.plugins.join(", "),
        pairs.len(),
        py_float_repr(options.tolerance)
    ));
    // An interrupt or a termination signal ends the run with its scratch removed.
    let mut taken = interrupt::Taken::interrupt();
    taken.termination();
    let compared = paths::make_scratch(options.work_dir.as_deref(), "vep-plugin-pairs-compare_")
        .map_err(|failure| failure.text())
        .and_then(|scratch| {
            let scratch = Removed(scratch);
            compare_all(&pairs, &run, &options, &scratch.0)
        });
    let stopped = interrupt::requested();
    drop(taken);
    if stopped {
        stderr_line("ERROR: [vep-plugin-pairs] interrupted before the report was written");
        return 130;
    }
    let results = match compared {
        Ok(results) => results,
        Err(message) => return fail(&message),
    };
    let (printed, any_failure) = report::printed(&results, &run, &options.min_annotated);
    say(&printed);

    if let Some(json_output) = given(&options.json_output) {
        let json = report::json(&results, &run, options.tolerance, &options.min_annotated);
        let written = std::fs::create_dir_all(path_parent(&json_output))
            .and_then(|()| std::fs::write(path_str(&json_output), json));
        if let Err(error) = written {
            return fail(&format!("{json_output}: {error}"));
        }
        say(&format!("\nJSON results written to: {json_output}\n"));
    }

    if any_failure {
        say("\nOVERALL: FAIL (one or more plugins have discordances)\n");
        1
    } else {
        say("\nOVERALL: PASS\n");
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_pair_has_one_shard_per_so_many_bytes_of_its_longer_output() {
        let directory =
            std::env::temp_dir().join(format!("vep-plugin-pairs-count-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&directory);
        std::fs::create_dir_all(&directory).unwrap();
        let file = |name: &str, bytes: usize| {
            let path = directory.join(name);
            std::fs::write(&path, vec![b'x'; bytes]).unwrap();
            path.to_str().unwrap().to_string()
        };
        let (empty, short, long) = (file("empty", 0), file("short", 1000), file("long", 10_000));
        let packed = file("packed.txt.GZ", 100);
        let missing = directory.join("missing").to_str().unwrap().to_string();
        assert_eq!(shard_count([&empty, &empty], 64), 1);
        assert_eq!(shard_count([&short, &long], 1000), 10);
        assert_eq!(shard_count([&long, &short], 1000), 10);
        assert_eq!(shard_count([&long, &short], 3000), 4);
        assert_eq!(shard_count([&long, &short], 1 << 27), 1);
        assert_eq!(shard_count([&short, &packed], 1000), 2);
        assert_eq!(shard_count([&long, &missing], 1), 4096);
        assert_eq!(shard_count([&missing, &missing], 1), 1);
        std::fs::remove_dir_all(directory).unwrap();
    }
}
