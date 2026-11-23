// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The golden corpus builder. A corpus is a small VCF of exemplar records, one or more for
//! every consequence set an Ensembl VEP reference output holds, with a manifest that records
//! how it was selected, the runs the golden tests make over it and the documented
//! divergences of each run. `select` builds a corpus, `classify` writes a run's divergences,
//! `check` holds a corpus against its references, `identifiers` replaces the identifiers
//! under restricted redistribution terms in stored files, `covers` and `flag-set` keep
//! the manifest's lists, and `generated` writes and checks the record of the files a corpus
//! does not commit. Every byte the tool writes is fixed; `README.md` states the contract.

#[path = "../../vep-formats-compare/src/rules/allele_number.rs"]
pub mod allele_number;
#[path = "../../vep-formats-compare/src/rules/coding_end.rs"]
pub mod coding_end;
#[path = "../../vep-formats-compare/src/genetic_codes.rs"]
pub mod genetic_codes;
#[path = "../../vep-outputs-compare/src/input.rs"]
pub mod input;
#[path = "../../vep-outputs-compare/src/pyfmt.rs"]
pub mod pyfmt;
#[path = "../../vep-outputs-compare/src/rules/mod.rs"]
pub mod rules;
#[path = "../../vep-formats-compare/src/rules/shifted_start_test.rs"]
pub mod shifted_start_test;
#[path = "../../vep-formats-compare/src/rules/start_test_past_the_allele.rs"]
pub mod start_test_past_the_allele;

pub mod cache_cost;
pub mod check;
pub mod classes;
pub mod classify;
pub mod cli;
pub mod constants;
pub mod covers;
pub mod error;
pub mod flag_set;
pub mod generated;
pub mod identifiers;
pub mod json;
pub mod manifest;
pub mod mate_context;
pub mod names;
pub mod pytext;
pub mod rows;
pub mod sampling;
pub mod scratch;
pub mod select;

pub use scratch::TempDir;

use error::{finish, say_raw, Failure};
use std::ffi::OsString;
use std::io::Write;

/// Runs the tool on the arguments after the program's name and returns its exit code: 0
/// after the work or the help, 1 after a refusal of the run, 2 after a refused command
/// line.
pub fn run(arguments: Vec<OsString>) -> i32 {
    let outcome = text_arguments(arguments).and_then(|arguments| match cli::parse(&arguments) {
        cli::Parsed::Help(text) => {
            say_raw(text);
            Ok(0)
        }
        cli::Parsed::Refused(text) => {
            let _ = std::io::stderr().lock().write_all(text.as_bytes());
            Ok(2)
        }
        cli::Parsed::Run(command) => match *command {
            cli::Command::Select(args) => select::run(&args),
            // The readers `classify` and `identifiers` open read a bound of no shape as one
            // thread, so the bound is refused here, as `select` and `check` refuse it before
            // their first read.
            cli::Command::Classify(args) => {
                rows::blocks::thread_count().and_then(|_| classify::run(&args))
            }
            cli::Command::Check(args) => check::run(&args),
            cli::Command::Identifiers(args) => {
                rows::blocks::thread_count().and_then(|_| identifiers::run(&args))
            }
            cli::Command::Covers(args) => covers::run(&args),
            cli::Command::FlagSet(args) => flag_set::run(&args),
            cli::Command::Generated(args) => generated::run(&args),
        },
    });
    match outcome {
        Ok(code) => finish(code),
        Err(failure) => {
            let _ = writeln!(std::io::stderr().lock(), "{}", failure.line);
            finish(failure.code)
        }
    }
}

/// The arguments as text; one that is not UTF-8 is refused.
fn text_arguments(arguments: Vec<OsString>) -> Result<Vec<String>, Failure> {
    arguments
        .into_iter()
        .enumerate()
        .map(|(place, argument)| {
            argument
                .into_string()
                .map_err(|_| Failure::own(format!("argument {} is not UTF-8", place + 1)))
        })
        .collect()
}
