// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How a run speaks and ends: the lines of standard output, the error lines of standard
//! error, and the failure that ends a subcommand with one line and an exit code.

use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};

/// The name every error line of the program's own carries in brackets, a fixed part of the
/// line.
pub const PROGRAM: &str = "build_golden_corpus";

/// The name in brackets on a line the tool writes where the contract names no line.
pub const TOOL: &str = "vep-golden-corpus";

/// A run that ends before its work is done: one line on standard error and an exit code. A
/// subcommand returns it where the contract says the run is refused with a message
/// (`Err(failure(...))`), and returns `Ok(code)` after error lines of its own
/// (`complain`) or after its work.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Failure {
    /// The whole line, without its newline.
    pub line: String,
    pub code: i32,
}

impl Failure {
    /// A refusal the tool states in its own name: `ERROR: [vep-golden-corpus] <detail>`,
    /// exit code 1.
    pub fn own(detail: impl AsRef<str>) -> Failure {
        Failure {
            line: format!("ERROR: [{TOOL}] {}", detail.as_ref()),
            code: 1,
        }
    }
}

/// A refusal of the program's own: `ERROR: [build_golden_corpus] <message>`, exit code 1.
pub fn failure(message: impl AsRef<str>) -> Failure {
    Failure {
        line: format!("ERROR: [{PROGRAM}] {}", message.as_ref()),
        code: 1,
    }
}

static OUTPUT_LOST: AtomicBool = AtomicBool::new(false);

/// Writes `text` and a newline to standard output.
pub fn say(text: &str) {
    say_raw(text);
    say_raw("\n");
}

/// Writes `text` to standard output as it is.
pub fn say_raw(text: &str) {
    if std::io::stdout().lock().write_all(text.as_bytes()).is_err() {
        OUTPUT_LOST.store(true, Ordering::Relaxed);
    }
}

/// Writes `ERROR: [build_golden_corpus] <message>` to standard error: an error line after
/// which the subcommand goes on to its next line or returns its exit code.
pub fn complain(message: &str) {
    let _ = writeln!(std::io::stderr().lock(), "ERROR: [{PROGRAM}] {message}");
}

/// Writes `WARNING: [build_golden_corpus] <message>` to standard error: a line about a run
/// that goes on.
pub fn warn(message: &str) {
    let _ = writeln!(std::io::stderr().lock(), "WARNING: [{PROGRAM}] {message}");
}

/// The exit code of a run that ends with `code`: standard output is flushed, and a run
/// whose output could not be written ends with exit code 1 and a line that says so.
pub fn finish(code: i32) -> i32 {
    if std::io::stdout().lock().flush().is_err() || OUTPUT_LOST.load(Ordering::Relaxed) {
        let _ = writeln!(
            std::io::stderr().lock(),
            "ERROR: [{TOOL}] standard output is not writable"
        );
        return 1;
    }
    code
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_two_failures_carry_their_lines_and_codes() {
        assert_eq!(
            failure("covers is not a list of names"),
            Failure {
                line: "ERROR: [build_golden_corpus] covers is not a list of names".to_string(),
                code: 1
            }
        );
        assert_eq!(
            Failure::own("corpus/manifest.json is not writable"),
            Failure {
                line: "ERROR: [vep-golden-corpus] corpus/manifest.json is not writable".to_string(),
                code: 1
            }
        );
    }
}
