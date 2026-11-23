// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Prepares a directory of benchmark VCFs for a concordance run.
//!
//! Every `.vcf` and `.vcf.gz` file of the input directory is written to the output
//! directory as plain text: header lines as they are, variant lines under one contig
//! naming and, under the canonical-contig filter, on the canonical contigs alone, the first
//! lines only where a cap is given. `manifest.json` records the run and, for each file, its
//! counts and the SHA256 digest of what was written. The prepared files, the manifest and the
//! stdout lines are fixed to the byte; `README.md` states the contract.

mod cli;
mod contigs;
mod gzip;
mod lines;
mod manifest;
mod paths;
mod pipeline;
mod pyfmt;

use cli::{Args, ChromMode, Command, ContigsFlag};
use contigs::{ContigSet, NoSet};
use lines::{Mode, Rules};
use manifest::{FileSummary, Manifest};
use paths::{as_path, Style};
use pipeline::{Pool, Preparation};
use std::collections::HashMap;
use std::ffi::OsString;
use std::io::{self, Write};
use std::num::NonZero;
use std::panic::AssertUnwindSafe;
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex};
use std::time::SystemTime;

/// The environment variable that bounds the worker threads; unset or empty, the host's
/// available cores are the bound. An environment variable and not an option: the command
/// line names the preparation, and the thread count changes no byte of it.
const THREADS_ENV: &str = "VEP_PREPARE_BENCHMARK_VCFS_THREADS";

/// The inputs prepared at once, at most.
const FILES_AT_ONCE: usize = 8;
/// The pieces of one input that are read and still unwritten, at most.
const PIECES_IN_FLIGHT: usize = 12;

/// How a run ends without its manifest.
enum Failure {
    /// A line for stderr, whole.
    Line(String),
    /// What failed, for the tool's own error line.
    Detail(String),
}

fn detail(context: &[u8], error: impl std::fmt::Display) -> Failure {
    Failure::Detail(format!("{}: {error}", pyfmt::py_stderr_os_str(context)))
}

fn thread_count() -> Result<usize, Failure> {
    let available = || std::thread::available_parallelism().map_or(1, NonZero::get);
    match std::env::var_os(THREADS_ENV) {
        None => Ok(available()),
        Some(value) if value.is_empty() => Ok(available()),
        Some(value) => value
            .to_str()
            .and_then(|text| text.trim().parse::<usize>().ok())
            .filter(|threads| *threads >= 1)
            .ok_or_else(|| {
                Failure::Detail(format!(
                    "{THREADS_ENV}={}: expected a whole number of threads, 1 or more",
                    value.to_string_lossy()
                ))
            }),
    }
}

/// The naming a file's variant lines are given: the one asked for, or under `auto` the
/// cache's where the input's first variant line and the cache are named differently.
fn resolve_mode(requested: ChromMode, input: Style, cache: Style) -> Mode {
    match requested {
        ChromMode::Keep => Mode::Keep,
        ChromMode::Strip => Mode::Strip,
        ChromMode::Add => Mode::Add,
        ChromMode::Auto => match (input, cache) {
            (Style::Chr, Style::Plain) => Mode::Strip,
            (Style::Plain, Style::Chr) => Mode::Add,
            _ => Mode::Keep,
        },
    }
}

fn write_stdout(bytes: &[u8]) -> Result<(), Failure> {
    let mut stdout = io::stdout().lock();
    stdout
        .write_all(bytes)
        .and_then(|()| stdout.flush())
        .map_err(|e| Failure::Detail(format!("stdout: {e}")))
}

/// What every input of a run shares.
struct Run<'a> {
    input_dir: &'a [u8],
    output_dir: &'a [u8],
    requested: ChromMode,
    cache_style: Style,
    contigs: Option<Arc<ContigSet>>,
    cap: Option<u64>,
    pool: &'a Pool,
    in_flight: usize,
}

impl Run<'_> {
    /// Prepares the input `name`: reads its first variant line for its naming, then
    /// writes its prepared file.
    fn prepare(&self, name: &[u8]) -> Result<FileSummary, String> {
        let input_file = paths::join(self.input_dir, name);
        let output_file = paths::join(self.output_dir, paths::prepared_name(name));
        let gzip = paths::is_gzip_name(name);
        let failed =
            |path: &[u8], error: io::Error| format!("{}: {error}", pyfmt::py_stderr_os_str(path));
        let input_style = pipeline::first_variant_style(as_path(&input_file), gzip)
            .map_err(|e| failed(&input_file, e))?;
        let mode = resolve_mode(self.requested, input_style, self.cache_style);
        std::fs::create_dir_all(as_path(self.output_dir))
            .map_err(|e| failed(self.output_dir, e))?;
        let rules = Arc::new(Rules {
            mode,
            contigs: self.contigs.clone(),
        });
        let prepared = Preparation {
            rules: &rules,
            cap: self.cap,
            pool: self.pool,
            in_flight: self.in_flight,
        }
        .file(as_path(&input_file), as_path(&output_file), gzip)
        .map_err(|e| failed(&input_file, e))?;
        Ok(FileSummary {
            input_file,
            output_file,
            input_style: input_style.name(),
            cache_style: self.cache_style.name(),
            mode: mode.name(),
            counts: prepared.counts,
            sha256: prepared.sha256,
        })
    }

    /// The stdout line of one prepared file.
    fn line(&self, name: &[u8], file: &FileSummary) -> Vec<u8> {
        let mut line = b"prepared ".to_vec();
        line.extend_from_slice(name);
        line.extend_from_slice(b" -> ");
        line.extend_from_slice(paths::prepared_name(name));
        let counts = &file.counts;
        let mut rest = format!(
            " (mode={}, variants={}, headers={}",
            file.mode, counts.variants, counts.headers
        );
        if self.contigs.is_some() {
            rest.push_str(&format!(", canonical_dropped={}", counts.dropped));
            if counts.chr_stripped {
                rest.push_str(" chr_stripped=true");
            }
        }
        rest.push_str(")\n");
        line.extend_from_slice(rest.as_bytes());
        line
    }

    /// Prepares every input and writes its stdout line, the lines in input order. Inputs
    /// that are written under one name are prepared one after another in input order;
    /// the others at once. After a failure no input behind the failed one is started, the
    /// inputs before it are prepared and reported, and the failure of the first is returned.
    fn prepare_all(
        &self,
        names: &[Vec<u8>],
        files_at_once: usize,
    ) -> Result<Vec<FileSummary>, Failure> {
        let mut by_output: HashMap<&[u8], usize> = HashMap::new();
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for (index, name) in names.iter().enumerate() {
            let group = *by_output
                .entry(paths::prepared_name(name))
                .or_insert(groups.len());
            if group == groups.len() {
                groups.push(Vec::new());
            }
            groups[group].push(index);
        }
        let next_group = AtomicUsize::new(0);
        let first_failed = AtomicUsize::new(usize::MAX);
        let done: Mutex<Vec<Option<Result<FileSummary, String>>>> =
            Mutex::new(names.iter().map(|_| None).collect());
        let changed = Condvar::new();
        let work = || loop {
            let group = next_group.fetch_add(1, Ordering::Relaxed);
            let Some(indexes) = groups.get(group) else {
                return;
            };
            for &index in indexes {
                if index > first_failed.load(Ordering::Relaxed) {
                    break;
                }
                let name = &names[index];
                let result = std::panic::catch_unwind(AssertUnwindSafe(|| self.prepare(name)))
                    .unwrap_or_else(|_| {
                        Err(format!(
                            "{}: a thread of the run ended early",
                            pyfmt::py_stderr_os_str(name)
                        ))
                    });
                let failed = result.is_err();
                if failed {
                    first_failed.fetch_min(index, Ordering::Relaxed);
                }
                done.lock().expect("no panic under the lock")[index] = Some(result);
                changed.notify_all();
                if failed {
                    break;
                }
            }
        };
        std::thread::scope(|scope| {
            for _ in 0..files_at_once.min(groups.len()).max(1) {
                scope.spawn(work);
            }
            let mut files = Vec::with_capacity(names.len());
            for (index, name) in names.iter().enumerate() {
                let mut slots = done.lock().expect("no panic under the lock");
                while slots[index].is_none() {
                    slots = changed.wait(slots).expect("no panic under the lock");
                }
                let result = slots[index].take().expect("waited for above");
                drop(slots);
                let file = result.map_err(Failure::Detail)?;
                write_stdout(&self.line(name, &file))?;
                files.push(file);
            }
            Ok(files)
        })
    }
}

fn run(args: &Args) -> Result<(), Failure> {
    let threads = thread_count()?;
    let resolved = |arg: &[u8]| paths::resolve(arg).map_err(|e| detail(arg, e));
    let input_dir = resolved(&args.input_dir)?;
    let output_dir = resolved(&args.output_dir)?;
    let cache_dir = match &args.cache_dir {
        Some(arg) if !arg.is_empty() => Some(resolved(arg)?),
        _ => None,
    };
    let (contigs, policy) = match &args.canonical_contigs {
        ContigsFlag::Off => (None, "disabled"),
        ContigsFlag::Named(value) => match ContigSet::named(value) {
            Ok(set) => (Some(Arc::new(set)), "named_set"),
            Err(NoSet::Empty(line)) => return Err(Failure::Line(line)),
            Err(NoSet::Unreadable(why)) => return Err(Failure::Detail(why)),
        },
        ContigsFlag::Bare => match args.assembly.as_deref() {
            None | Some(b"GRCh37") | Some(b"GRCh38") => (
                Some(Arc::new(ContigSet::human())),
                "universal_ensembl_style",
            ),
            Some(assembly) => {
                let mut warning =
                    b"WARNING: --canonical-contigs requested but --assembly is ".to_vec();
                warning.extend_from_slice(assembly);
                warning.extend_from_slice(
                    b", which names no canonical set; the filter is OFF for this run. Name the \
                      set with --canonical-contigs <file or list> (a cache's info.json \
                      contig_order) or pass --no-canonical-contigs to silence this warning.\n",
                );
                write_stdout(&warning)?;
                (None, "disabled")
            }
        },
    };
    std::fs::create_dir_all(as_path(&output_dir)).map_err(|e| detail(&output_dir, e))?;
    let names = paths::list_inputs(&input_dir);
    if names.is_empty() {
        return Err(Failure::Detail(format!(
            "No .vcf or .vcf.gz files found in {}",
            pyfmt::py_stderr_os_str(&input_dir)
        )));
    }
    let cache_style = paths::cache_style(cache_dir.as_deref())
        .map_err(|e| detail(cache_dir.as_deref().unwrap_or_default(), e))?;
    let files_at_once = threads.min(FILES_AT_ONCE);
    let in_flight = threads.clamp(2, PIECES_IN_FLIGHT);
    // More workers than pieces can be in flight would stand idle.
    let pool = Pool::new(threads.min(files_at_once * in_flight));
    let run = Run {
        input_dir: &input_dir,
        output_dir: &output_dir,
        requested: args.chrom_mode,
        cache_style,
        contigs,
        cap: args.max_variants.as_ref().map(|max| max.count),
        pool: &pool,
        in_flight,
    };
    let files = run.prepare_all(&names, files_at_once)?;
    let manifest = Manifest {
        generated_at_utc: manifest::utc_isoformat(SystemTime::now()),
        input_dir: &input_dir,
        output_dir: &output_dir,
        requested_chrom_mode: args.chrom_mode.name(),
        cache_dir: cache_dir.as_deref(),
        max_variants: args.max_variants.as_ref().map(|max| max.decimal.as_str()),
        assembly: args.assembly.as_deref(),
        policy,
        contigs: run.contigs.as_ref().map(|set| set.sorted()),
        files: &files,
    };
    let manifest_path = paths::join(&output_dir, b"manifest.json");
    std::fs::write(as_path(&manifest_path), manifest.dumps())
        .map_err(|e| detail(&manifest_path, e))?;
    let mut wrote = b"wrote manifest: ".to_vec();
    wrote.extend_from_slice(&manifest_path);
    wrote.push(b'\n');
    write_stdout(&wrote)
}

fn main() -> ExitCode {
    let argv: Vec<OsString> = std::env::args_os().skip(1).collect();
    let args = match cli::parse(&argv) {
        Ok(Command::Run(args)) => args,
        Ok(Command::Help) => {
            println!("{}\n\n{}", cli::USAGE, cli::HELP);
            return ExitCode::SUCCESS;
        }
        Err(e) => {
            eprintln!("{}\n{}: error: {e}", cli::USAGE, cli::PROG);
            return ExitCode::from(2);
        }
    };
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(Failure::Line(line)) => {
            eprintln!("{line}");
            ExitCode::FAILURE
        }
        Err(Failure::Detail(why)) => {
            eprintln!("ERROR: [{}] {why}", cli::PROG);
            ExitCode::FAILURE
        }
    }
}
