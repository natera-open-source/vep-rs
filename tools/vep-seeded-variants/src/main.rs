// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `vep-seeded-variants` command: reads one converted cache and writes its records
//! as a VCF, with the summary beside it on request.

use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use vep_seeded_variants::cache::{Cache, SEQUENCE_REGIONS_FILE};
use vep_seeded_variants::fasta::Fasta;
use vep_seeded_variants::generate::{generate, Options, Records, WIDE_RECORDS};
use vep_seeded_variants::vcf::{write_summary, write_vcf, Header};
use vep_seeded_variants::GENERATOR_VERSION;

const COMPONENT: &str = "vep-seeded-variants";

const USAGE: &str = "\
usage:
  vep-seeded-variants --cache <dir> --seed <n> --output <file.vcf>
      [--records light|wide] [--per-category <n>] [--sequence <name>]...
      [--fasta <file>] [--sequence-regions <file>] [--summary <file.tsv>]
  vep-seeded-variants --version";

/// How a run ends short of success: a command line the tool cannot read (exit 2) or a
/// failure (exit 1).
enum Stop {
    Usage(String),
    Failed(String),
}

/// What the command line asks for.
struct Request {
    cache: PathBuf,
    output: PathBuf,
    summary: Option<PathBuf>,
    fasta: Option<PathBuf>,
    sequence_regions: Option<PathBuf>,
    options: Options,
}

/// What the command line is: a run, or a question the tool answers and ends on.
enum Command {
    Run(Request),
    Version,
    Help,
}

/// The options that take a value; `--sequence` alone may be given more than once.
const VALUED: [&str; 9] = [
    "cache",
    "seed",
    "output",
    "records",
    "per-category",
    "sequence",
    "fasta",
    "sequence-regions",
    "summary",
];

fn read(args: &[String]) -> Result<Command, Stop> {
    let mut given: Vec<(String, String)> = Vec::new();
    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        let Some(option) = arg.strip_prefix("--") else {
            return Err(Stop::Usage(format!("{arg} is not an option")));
        };
        let (name, inline) = match option.split_once('=') {
            Some((name, value)) => (name, Some(value.to_owned())),
            None => (option, None),
        };
        match name {
            "version" | "help" if inline.is_some() => {
                return Err(Stop::Usage(format!("--{name} takes no value")));
            }
            "version" => return Ok(Command::Version),
            "help" => return Ok(Command::Help),
            _ => {}
        }
        if !VALUED.contains(&name) {
            return Err(Stop::Usage(format!("no option --{name}")));
        }
        if name != "sequence" && given.iter().any(|(earlier, _)| earlier == name) {
            return Err(Stop::Usage(format!("--{name} is given twice")));
        }
        let value = inline
            .or_else(|| rest.next().cloned())
            .ok_or_else(|| Stop::Usage(format!("--{name} takes a value")))?;
        given.push((name.to_owned(), value));
    }
    let value = |name: &str| {
        given
            .iter()
            .find(|(given, _)| given == name)
            .map(|(_, value)| value.as_str())
    };
    let required =
        |name: &str| value(name).ok_or_else(|| Stop::Usage(format!("--{name} is required")));
    let seed: u64 = required("seed")?
        .parse()
        .map_err(|_| Stop::Usage("--seed takes a whole number from 0 to 2^64 - 1".to_owned()))?;
    let per_category = value("per-category")
        .map(|count| match count.parse::<usize>() {
            Ok(count) if count >= 1 => Ok(count),
            _ => Err(Stop::Usage(
                "--per-category takes a whole number above zero".to_owned(),
            )),
        })
        .transpose()?;
    let records = match (value("records").unwrap_or("light"), per_category) {
        ("light", None) => Records::Light,
        ("light", Some(_)) => {
            return Err(Stop::Usage(
                "--per-category goes with --records wide".to_owned(),
            ))
        }
        ("wide", count) => Records::Wide(count.unwrap_or(WIDE_RECORDS)),
        (other, _) => {
            return Err(Stop::Usage(format!(
                "--records takes light or wide, not {other}"
            )))
        }
    };
    let sequences = given
        .iter()
        .filter(|(name, _)| name == "sequence")
        .map(|(_, value)| value.clone())
        .collect();
    let cache = PathBuf::from(required("cache")?);
    let output = PathBuf::from(required("output")?);
    let summary = value("summary").map(PathBuf::from);
    if summary.as_ref() == Some(&output) {
        return Err(Stop::Usage(
            "--summary names the file --output names".to_owned(),
        ));
    }
    Ok(Command::Run(Request {
        cache,
        output,
        summary,
        fasta: value("fasta").map(PathBuf::from),
        sequence_regions: value("sequence-regions").map(PathBuf::from),
        options: Options {
            seed,
            records,
            sequences,
        },
    }))
}

/// Writes a file through `write`, whole or not at all as far as an error allows.
fn write_file(
    path: &Path,
    write: impl FnOnce(&mut BufWriter<File>) -> io::Result<()>,
) -> Result<(), Stop> {
    let failed = |error: io::Error| Stop::Failed(format!("{}: {error}", path.display()));
    let mut out = BufWriter::new(File::create(path).map_err(failed)?);
    write(&mut out).map_err(failed)?;
    out.flush().map_err(failed)
}

fn run(request: &Request) -> Result<(), Stop> {
    let side_file = match &request.sequence_regions {
        Some(path) => Some(path.clone()),
        None => Some(request.cache.join(SEQUENCE_REGIONS_FILE)).filter(|path| path.is_file()),
    };
    let cache = Cache::open_with_regions(&request.cache, side_file.as_deref())
        .map_err(|error| Stop::Failed(error.to_string()))?;
    let mut fasta = request
        .fasta
        .as_deref()
        .map(Fasta::open)
        .transpose()
        .map_err(Stop::Failed)?;
    let generated = generate(&cache, fasta.as_mut(), &request.options)
        .map_err(|failure| Stop::Failed(failure.to_string()))?;
    let header = Header {
        seed: request.options.seed,
        records: request.options.records.describe(),
        reference_alleles: if fasta.is_some() { "fasta" } else { "cache" }.to_owned(),
    };
    write_file(&request.output, |out| {
        write_vcf(out, &header, &generated.records)
    })?;
    if let Some(summary) = &request.summary {
        write_file(summary, |out| write_summary(out, &generated.summary))?;
    }
    Ok(())
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
    let outcome = arguments().and_then(|args| match read(&args)? {
        Command::Version => {
            println!(
                "{COMPONENT} {} (generator version {GENERATOR_VERSION})",
                env!("CARGO_PKG_VERSION")
            );
            Ok(())
        }
        Command::Help => {
            println!("{USAGE}");
            Ok(())
        }
        Command::Run(request) => run(&request),
    });
    match outcome {
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
