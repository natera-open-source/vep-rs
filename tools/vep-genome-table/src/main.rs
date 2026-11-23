// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `vep-genome-table` command: `check`, `row`, `keys` and `put` over a genome table
//! directory, and `digest` over a default-format output.

use std::env;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use vep_genome_table::{check, digest_file, header, load, put, Failure, COLUMNS};

const COMPONENT: &str = "vep-genome-table";

const USAGE: &str = "\
usage:
  vep-genome-table check <dir> [--complete]
  vep-genome-table row <dir> <genome key> [--field <column>]
  vep-genome-table keys <dir> [--set <set>] [--outcome <prefix>]
  vep-genome-table put <dir> <row file>...
  vep-genome-table digest <output file>";

/// How a run ends short of success: a command line the tool cannot read (exit 2) or a
/// refusal, one line a failure (exit 1).
enum Stop {
    Usage(String),
    Failed(Vec<String>),
}

impl From<Vec<Failure>> for Stop {
    fn from(failures: Vec<Failure>) -> Stop {
        Stop::Failed(failures.iter().map(Failure::to_string).collect())
    }
}

/// A subcommand's arguments: the positional ones in order and each `--name` given.
struct Arguments {
    positional: Vec<String>,
    options: Vec<(String, Option<String>)>,
}

impl Arguments {
    /// Splits `args` into positional arguments and options. An option of `valued` takes
    /// its value from `--name=value` or from the next argument; one of `switches` takes
    /// none; any other `--name`, and an option given twice, is refused.
    fn read(args: &[String], valued: &[&str], switches: &[&str]) -> Result<Arguments, Stop> {
        let mut positional = Vec::new();
        let mut options: Vec<(String, Option<String>)> = Vec::new();
        let mut rest = args.iter();
        while let Some(arg) = rest.next() {
            let Some(option) = arg.strip_prefix("--") else {
                positional.push(arg.clone());
                continue;
            };
            let (name, inline) = match option.split_once('=') {
                Some((name, value)) => (name, Some(value.to_owned())),
                None => (option, None),
            };
            if options.iter().any(|(given, _)| given == name) {
                return Err(Stop::Usage(format!("--{name} is given twice")));
            }
            let value = if valued.contains(&name) {
                let value = inline
                    .or_else(|| rest.next().cloned())
                    .ok_or_else(|| Stop::Usage(format!("--{name} takes a value")))?;
                Some(value)
            } else if switches.contains(&name) {
                if inline.is_some() {
                    return Err(Stop::Usage(format!("--{name} takes no value")));
                }
                None
            } else {
                return Err(Stop::Usage(format!("no option --{name}")));
            };
            options.push((name.to_owned(), value));
        }
        Ok(Arguments {
            positional,
            options,
        })
    }

    fn has(&self, name: &str) -> bool {
        self.options.iter().any(|(given, _)| given == name)
    }

    fn value(&self, name: &str) -> Option<&str> {
        self.options
            .iter()
            .find(|(given, _)| given == name)
            .and_then(|(_, value)| value.as_deref())
    }
}

/// Writes `text` to standard output. A reader that closed the pipe ends the output
/// without a failure.
fn print(text: &str) -> Result<(), Stop> {
    match io::stdout().lock().write_all(text.as_bytes()) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == io::ErrorKind::BrokenPipe => Ok(()),
        Err(error) => Err(Stop::Failed(vec![format!("standard output: {error}")])),
    }
}

fn run(args: &[String]) -> Result<(), Stop> {
    let Some((command, rest)) = args.split_first() else {
        return Err(Stop::Usage("no subcommand".to_owned()));
    };
    match command.as_str() {
        "check" => {
            let arguments = Arguments::read(rest, &[], &["complete"])?;
            let [dir] = arguments.positional.as_slice() else {
                return Err(Stop::Usage("check takes one directory".to_owned()));
            };
            let summary = check(Path::new(dir), arguments.has("complete"))?;
            print(&format!("{summary}: ok\n"))
        }
        "row" => {
            let arguments = Arguments::read(rest, &["field"], &[])?;
            let [dir, key] = arguments.positional.as_slice() else {
                return Err(Stop::Usage(
                    "row takes a directory and a genome key".to_owned(),
                ));
            };
            let field = arguments.value("field");
            if let Some(column) = field {
                if !COLUMNS.contains(&column) {
                    return Err(Stop::Usage(format!(
                        "no column {column}; the columns are {}",
                        COLUMNS.join(", ")
                    )));
                }
            }
            let rows = load(Path::new(dir))?;
            let Some(row) = rows.iter().find(|row| row.key() == *key) else {
                return Err(Stop::Failed(vec![format!(
                    "{dir}: no row has the key {key}"
                )]));
            };
            match field.and_then(|column| row.field(column)) {
                Some(value) => print(&format!("{value}\n")),
                None => print(&format!("{}\n{}\n", header(), row.line())),
            }
        }
        "keys" => {
            let arguments = Arguments::read(rest, &["set", "outcome"], &[])?;
            let [dir] = arguments.positional.as_slice() else {
                return Err(Stop::Usage("keys takes one directory".to_owned()));
            };
            let set = arguments.value("set");
            let outcome = arguments.value("outcome");
            let mut text = String::new();
            for row in load(Path::new(dir))? {
                if set.is_some_and(|set| row.set != set)
                    || outcome.is_some_and(|prefix| !row.outcome.starts_with(prefix))
                {
                    continue;
                }
                text.push_str(&row.key());
                text.push('\n');
            }
            print(&text)
        }
        "put" => {
            let arguments = Arguments::read(rest, &[], &[])?;
            let Some((dir, row_files)) = arguments.positional.split_first() else {
                return Err(Stop::Usage(
                    "put takes a directory and row files".to_owned(),
                ));
            };
            if row_files.is_empty() {
                return Err(Stop::Usage("put takes at least one row file".to_owned()));
            }
            let row_files: Vec<PathBuf> = row_files.iter().map(PathBuf::from).collect();
            let done = put(Path::new(dir), &row_files)?;
            print(&format!(
                "{}: {} added, {} replaced\n",
                done.summary, done.added, done.replaced
            ))
        }
        "digest" => {
            let arguments = Arguments::read(rest, &[], &[])?;
            let [output] = arguments.positional.as_slice() else {
                return Err(Stop::Usage("digest takes one output file".to_owned()));
            };
            let digest = digest_file(Path::new(output))
                .map_err(|error| Stop::Failed(vec![format!("{output}: {error}")]))?;
            print(&format!("{digest}\n"))
        }
        "help" | "--help" | "-h" => print(&format!("{USAGE}\n")),
        other => Err(Stop::Usage(format!("no subcommand {other}"))),
    }
}

/// The command line after the program's name. An argument that is not UTF-8 is a command
/// line the tool cannot read, where `env::args` panics on one.
fn arguments() -> Result<Vec<String>, Stop> {
    env::args_os()
        .skip(1)
        .map(|arg| {
            arg.into_string().map_err(|arg| {
                Stop::Usage(format!(
                    "the argument {} is not UTF-8",
                    arg.to_string_lossy()
                ))
            })
        })
        .collect()
}

fn main() -> ExitCode {
    match arguments().and_then(|args| run(&args)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(Stop::Usage(what)) => {
            eprintln!("ERROR: [{COMPONENT}] {what}");
            eprintln!("{USAGE}");
            ExitCode::from(2)
        }
        Err(Stop::Failed(lines)) => {
            for line in lines {
                eprintln!("ERROR: [{COMPONENT}] {line}");
            }
            ExitCode::from(1)
        }
    }
}
