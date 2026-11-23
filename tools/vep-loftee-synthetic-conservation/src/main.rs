// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `vep-loftee-synthetic-conservation` command: reads one converted cache and writes
//! the made-up table of base scores, its tabix index and the made-up database of exon
//! scores into a directory, then prints what it wrote.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use vep_loftee_synthetic_conservation::generate::{generate, Options};
use vep_seeded_variants::cache::Cache;

const COMPONENT: &str = "vep-loftee-synthetic-conservation";

const USAGE: &str = "\
usage:
  vep-loftee-synthetic-conservation --cache <dir> --seed <n> --out-dir <dir>
      [--threads <n>] [--per-sequence <file.tsv>] [--force]
      [--skip-transcripts-without-exons]
  vep-loftee-synthetic-conservation --version";

/// The options that take a value.
const VALUED: [&str; 5] = ["cache", "seed", "out-dir", "threads", "per-sequence"];

/// The options that take none.
const SWITCHES: [&str; 2] = ["force", "skip-transcripts-without-exons"];

/// How a run ends short of success: a command line the tool cannot read (exit 2) or a
/// failure (exit 1).
enum Stop {
    Usage(String),
    Failed(String),
}

/// What the command line asks for.
enum Command {
    Write {
        cache: PathBuf,
        out_dir: PathBuf,
        per_sequence: Option<PathBuf>,
        options: Options,
    },
    Version,
    Help,
}

/// The options of the line with their values, a switch with an empty one. An argument
/// that is no option, an option outside the two lists, an option given twice, a valued
/// option without a value and a switch with one are each a command line the tool cannot
/// read.
fn split(args: &[String]) -> Result<Vec<(String, String)>, Stop> {
    let mut given: Vec<(String, String)> = Vec::new();
    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        let Some(option) = arg.strip_prefix("--") else {
            return Err(Stop::Usage(format!("{arg} is no option")));
        };
        let (name, inline) = match option.split_once('=') {
            Some((name, value)) => (name, Some(value.to_owned())),
            None => (option, None),
        };
        if given.iter().any(|(earlier, _)| earlier == name) {
            return Err(Stop::Usage(format!("--{name} is given twice")));
        }
        let value = if VALUED.contains(&name) {
            inline
                .or_else(|| rest.next().cloned())
                .ok_or_else(|| Stop::Usage(format!("--{name} takes a value")))?
        } else if !SWITCHES.contains(&name) {
            return Err(Stop::Usage(format!("no option --{name}")));
        } else if inline.is_some() {
            return Err(Stop::Usage(format!("--{name} takes no value")));
        } else {
            String::new()
        };
        given.push((name.to_owned(), value));
    }
    Ok(given)
}

fn read(args: &[String]) -> Result<Command, Stop> {
    match args {
        [only] if only == "--version" => return Ok(Command::Version),
        [only] if only == "--help" => return Ok(Command::Help),
        _ => {}
    }
    let given = split(args)?;
    let value = |name: &str| {
        given
            .iter()
            .find(|(given, _)| given == name)
            .map(|(_, value)| value.clone())
    };
    let required =
        |name: &str| value(name).ok_or_else(|| Stop::Usage(format!("--{name} is required")));
    let cache = PathBuf::from(required("cache")?);
    let seed = required("seed")?
        .parse::<u64>()
        .map_err(|_| Stop::Usage("--seed takes a whole number".to_owned()))?;
    let out_dir = PathBuf::from(required("out-dir")?);
    let threads = value("threads")
        .map(|count| match count.parse::<usize>() {
            Ok(count) if count >= 1 => Ok(count),
            _ => Err(Stop::Usage(
                "--threads takes a whole number above zero".to_owned(),
            )),
        })
        .transpose()?;
    Ok(Command::Write {
        cache,
        out_dir,
        per_sequence: value("per-sequence").map(PathBuf::from),
        options: Options {
            seed,
            threads,
            skip_transcripts_without_exons: value("skip-transcripts-without-exons").is_some(),
            force: value("force").is_some(),
        },
    })
}

/// Writes `text` to the file at `path`, through a file beside it that is renamed when
/// it is whole.
fn write_file(path: &Path, text: &str) -> Result<(), Stop> {
    let failed = |error: std::io::Error| Stop::Failed(format!("{}: {error}", path.display()));
    let mut partial = path.as_os_str().to_owned();
    partial.push(format!(".{}.tmp", std::process::id()));
    let partial = PathBuf::from(partial);
    fs::write(&partial, text)
        .and_then(|()| fs::rename(&partial, path))
        .map_err(|error| {
            let _ = fs::remove_file(&partial);
            failed(error)
        })
}

fn run(command: Command) -> Result<(), Stop> {
    match command {
        Command::Version => {
            println!("{COMPONENT} {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        Command::Help => {
            println!("{USAGE}");
            Ok(())
        }
        Command::Write {
            cache,
            out_dir,
            per_sequence,
            options,
        } => {
            if let Some(file) = per_sequence.as_deref() {
                if file.exists() && !options.force {
                    return Err(Stop::Failed(format!(
                        "{} exists; --force overwrites it",
                        file.display()
                    )));
                }
                // The file is renamed into place after the three files of the run, where
                // a rename over a directory fails and leaves them.
                if fs::symlink_metadata(file).is_ok_and(|entry| entry.is_dir()) {
                    return Err(Stop::Failed(format!("{} is a directory", file.display())));
                }
            }
            let cache = Cache::open(&cache).map_err(|error| Stop::Failed(error.to_string()))?;
            let summary = generate(&cache, &out_dir, &options).map_err(Stop::Failed)?;
            if let Some(file) = per_sequence.as_deref() {
                write_file(file, &summary.per_sequence())?;
            }
            print!("{}", summary.lines());
            Ok(())
        }
    }
}

/// The command line after the program's name. An argument that is not Unicode text is a
/// command line the tool cannot read, where `env::args` panics on one.
fn arguments() -> Result<Vec<String>, Stop> {
    env::args_os()
        .skip(1)
        .map(|arg| {
            arg.into_string()
                .map_err(|_| Stop::Usage("an argument is not Unicode text".to_owned()))
        })
        .collect()
}

fn main() -> ExitCode {
    match arguments().and_then(|args| read(&args)).and_then(run) {
        Ok(()) => ExitCode::SUCCESS,
        Err(Stop::Usage(what)) => {
            eprintln!("ERROR: [{COMPONENT}] {what}");
            eprintln!("{USAGE}");
            ExitCode::from(2)
        }
        Err(Stop::Failed(what)) => {
            eprintln!("ERROR: [{COMPONENT}] {what}");
            ExitCode::from(1)
        }
    }
}
