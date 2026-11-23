// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The plugin comparison's command line: the first word names the command, `score` or
//! `compare`, and the words after it are that command's own command line.

use std::io::Write;
use std::process::ExitCode;
use vep_plugin_pairs::{compare, score};

const USAGE: &str = "usage: vep-plugin-pairs {score,compare} ...";
const HELP: &str = "\
The plugin comparison of two VEP outputs.

commands:
  score    scores the plugin tuples of each engine's output against one Ensembl VEP
           default-format output and writes the counts a raw and an adjusted F1 are derived from
  compare  compares the fields plugins write between an Ensembl VEP output and a vep-rs output

Each command prints its options under --help.";

/// The exit status of a command that returns its own, kept within a byte.
fn status(code: i32) -> ExitCode {
    ExitCode::from(u8::try_from(code).unwrap_or(1))
}

/// Writes `text` and a line feed to a standard stream. A stream whose reader went away
/// takes nothing, and the exit status of the command line stands.
fn say(mut stream: impl Write, text: &str) {
    let _ = writeln!(stream, "{text}");
    let _ = stream.flush();
}

fn main() -> ExitCode {
    let mut words = std::env::args_os()
        .skip(1)
        .map(|word| word.to_string_lossy().into_owned());
    let refused = |reason: &str| {
        let line = format!("{USAGE}\nERROR: [vep-plugin-pairs] {reason}");
        say(std::io::stderr().lock(), &line);
        ExitCode::from(2)
    };
    let command: fn(Vec<String>) -> i32 = match words.next().as_deref() {
        Some("score") => score::main,
        Some("compare") => compare::main,
        Some("-h" | "--help") => {
            say(std::io::stdout().lock(), &format!("{USAGE}\n\n{HELP}"));
            return ExitCode::SUCCESS;
        }
        Some(word) => {
            return refused(&format!(
                "no command {word}: the first word is score or compare"
            ));
        }
        None => return refused("the first word names the command: score or compare"),
    };
    status(command(words.collect()))
}
