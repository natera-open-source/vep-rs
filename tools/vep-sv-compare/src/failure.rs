// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The one error type of the comparator: a failure that ends the run with exit code 1.

use std::fmt;

/// What kind of failure ended the comparison.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FailureKind {
    /// An input that cannot be scored: a truncated output row, an unknown release, a mask
    /// called without the records it reads.
    Value,
    /// A state the comparator refuses to score over: a scored-engine row the cache rules
    /// out, a registry rule without a mask.
    Assertion,
    /// A file that cannot be opened, read or written.
    Io,
    /// Input bytes that are not UTF-8.
    Decode,
}

/// A failure with the detail printed after `ERROR: [vep-sv-compare] ` on stderr.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Failure {
    pub kind: FailureKind,
    pub message: String,
}

impl Failure {
    pub fn value(message: impl Into<String>) -> Failure {
        Failure {
            kind: FailureKind::Value,
            message: message.into(),
        }
    }

    pub fn assertion(message: impl Into<String>) -> Failure {
        Failure {
            kind: FailureKind::Assertion,
            message: message.into(),
        }
    }

    pub fn io(message: impl Into<String>) -> Failure {
        Failure {
            kind: FailureKind::Io,
            message: message.into(),
        }
    }

    pub fn decode(message: impl Into<String>) -> Failure {
        Failure {
            kind: FailureKind::Decode,
            message: message.into(),
        }
    }
}

impl fmt::Display for Failure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for Failure {}
