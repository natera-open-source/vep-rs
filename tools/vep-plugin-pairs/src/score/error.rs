// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The refusal of the `score` program.

use crate::rules::Refusal;
use std::io;
use std::path::Path;

/// A refusal of the `score` program: the message of its one
/// `ERROR: [score_plugin_tuples] <message>` line, which ends the run with exit status 2.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScorerError(pub String);

pub type Scored<T> = Result<T, ScorerError>;

impl ScorerError {
    /// A failed operating-system call on `path` in the words CPython gives an `OSError`,
    /// `[Errno N] <reason>: '<path>'`, the path as `repr` writes it, every run of whitespace
    /// written as one space.
    pub fn os(error: &io::Error, path: &Path) -> ScorerError {
        let shown = path.to_string_lossy();
        let message = match error.raw_os_error() {
            Some(code) => {
                let told = error.to_string();
                let reason = told
                    .strip_suffix(&format!(" (os error {code})"))
                    .unwrap_or(&told);
                format!("[Errno {code}] {reason}: {}", crate::pytext::repr(&shown))
            }
            None => error.to_string(),
        };
        ScorerError(message.split_whitespace().collect::<Vec<_>>().join(" "))
    }
}

/// A class rule's refusal ends the run in the same line as the program's own.
impl From<Refusal> for ScorerError {
    fn from(refusal: Refusal) -> ScorerError {
        ScorerError(refusal.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn failed_system_call_reads_as_the_interpreter_words_it() {
        let missing = io::Error::from_raw_os_error(2);
        assert_eq!(
            ScorerError::os(&missing, Path::new("/no/such  file")).0,
            "[Errno 2] No such file or directory: '/no/such file'"
        );
        let denied = io::Error::from_raw_os_error(13);
        assert_eq!(
            ScorerError::os(&denied, Path::new("it's")).0,
            "[Errno 13] Permission denied: \"it's\""
        );
    }
}
