// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `vep-genome-features` command: `extract` reads one converted cache and writes its
//! features as JSON, and `gather` writes one table row a genome from a directory of such
//! files.

use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use vep_genome_features::extract::{extract, Request};
use vep_genome_features::features::FEATURES_VERSION;
use vep_genome_features::gather::gather;

const COMPONENT: &str = "vep-genome-features";

const USAGE: &str = "\
usage:
  vep-genome-features extract <converted cache> --out <features.json>
      [--sequence-regions <file>] [--genome <key>] [--threads <n>]
  vep-genome-features gather <directory> --out <features.tsv>
  vep-genome-features --version";

/// How a run ends short of success: a command line the tool cannot read (exit 2) or a
/// failure (exit 1).
enum Stop {
    Usage(String),
    Failed(String),
}

/// What the command line asks for.
enum Command {
    Extract {
        cache: PathBuf,
        out: PathBuf,
        sequence_regions: Option<PathBuf>,
        genome: Option<String>,
        threads: Option<usize>,
    },
    Gather {
        dir: PathBuf,
        out: PathBuf,
    },
    Version,
    Help,
}

/// The one positional argument and the valued options of a subcommand's line. An option
/// outside `valued`, an option given twice, a second positional argument and a missing
/// one are each a command line the tool cannot read.
fn split(args: &[String], valued: &[&str]) -> Result<(String, Vec<(String, String)>), Stop> {
    let mut positional: Option<String> = None;
    let mut given: Vec<(String, String)> = Vec::new();
    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        let Some(option) = arg.strip_prefix("--") else {
            if positional.is_some() {
                return Err(Stop::Usage(format!("{arg} is one argument too many")));
            }
            positional = Some(arg.clone());
            continue;
        };
        let (name, inline) = match option.split_once('=') {
            Some((name, value)) => (name, Some(value.to_owned())),
            None => (option, None),
        };
        if !valued.contains(&name) {
            return Err(Stop::Usage(format!("no option --{name}")));
        }
        if given.iter().any(|(earlier, _)| earlier == name) {
            return Err(Stop::Usage(format!("--{name} is given twice")));
        }
        let value = inline
            .or_else(|| rest.next().cloned())
            .ok_or_else(|| Stop::Usage(format!("--{name} takes a value")))?;
        given.push((name.to_owned(), value));
    }
    let positional =
        positional.ok_or_else(|| Stop::Usage("the directory to read is missing".to_owned()))?;
    Ok((positional, given))
}

fn read(args: &[String]) -> Result<Command, Stop> {
    let Some((first, rest)) = args.split_first() else {
        return Err(Stop::Usage("no subcommand".to_owned()));
    };
    let value = |given: &[(String, String)], name: &str| {
        given
            .iter()
            .find(|(given, _)| given == name)
            .map(|(_, value)| value.clone())
    };
    let out = |given: &[(String, String)]| {
        value(given, "out")
            .map(PathBuf::from)
            .ok_or_else(|| Stop::Usage("--out is required".to_owned()))
    };
    match first.as_str() {
        "--version" if rest.is_empty() => Ok(Command::Version),
        "--help" if rest.is_empty() => Ok(Command::Help),
        "extract" => {
            let (cache, given) = split(rest, &["out", "sequence-regions", "genome", "threads"])?;
            let threads = value(&given, "threads")
                .map(|count| match count.parse::<usize>() {
                    Ok(count) if count >= 1 => Ok(count),
                    _ => Err(Stop::Usage(
                        "--threads takes a whole number above zero".to_owned(),
                    )),
                })
                .transpose()?;
            Ok(Command::Extract {
                cache: PathBuf::from(cache),
                out: out(&given)?,
                sequence_regions: value(&given, "sequence-regions").map(PathBuf::from),
                genome: value(&given, "genome"),
                threads,
            })
        }
        "gather" => {
            let (dir, given) = split(rest, &["out"])?;
            Ok(Command::Gather {
                dir: PathBuf::from(dir),
                out: out(&given)?,
            })
        }
        other => Err(Stop::Usage(format!("no subcommand {other}"))),
    }
}

/// Writes `text` to the file at `path`.
fn write_file(path: &Path, text: &str) -> Result<(), Stop> {
    let failed = |error: io::Error| Stop::Failed(format!("{}: {error}", path.display()));
    let mut out = BufWriter::new(File::create(path).map_err(failed)?);
    out.write_all(text.as_bytes()).map_err(failed)?;
    out.flush().map_err(failed)
}

fn run(command: Command) -> Result<(), Stop> {
    match command {
        Command::Version => {
            println!(
                "{COMPONENT} {} (features version {FEATURES_VERSION})",
                env!("CARGO_PKG_VERSION")
            );
            Ok(())
        }
        Command::Help => {
            println!("{USAGE}");
            Ok(())
        }
        Command::Extract {
            cache,
            out,
            sequence_regions,
            genome,
            threads,
        } => {
            let request = Request {
                sequence_regions: sequence_regions.as_deref(),
                genome: genome.as_deref(),
                threads,
            };
            let features = extract(&cache, &request).map_err(Stop::Failed)?;
            let mut text = serde_json::to_string_pretty(&features)
                .map_err(|error| Stop::Failed(error.to_string()))?;
            text.push('\n');
            write_file(&out, &text)
        }
        Command::Gather { dir, out } => {
            let table = gather(&dir).map_err(Stop::Failed)?;
            write_file(&out, &table)
        }
    }
}

/// The command line after the program's name. An argument that is not UTF-8 is a command
/// line the tool cannot read, where `env::args` panics on one.
fn arguments() -> Result<Vec<String>, Stop> {
    env::args_os()
        .skip(1)
        .map(|arg| {
            arg.into_string()
                .map_err(|_| Stop::Usage("an argument is not UTF-8".to_owned()))
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
