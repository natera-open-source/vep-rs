// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparator's own error, and the payload of a run that ends on an input no reading of
//! the comparison expects.

use std::fmt;

/// A malformed input or a refused argument: the run ends with
/// `ERROR: [compare_vep_formats] <text>` and exit code 2.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    Comparator(String),
}

impl Error {
    pub fn comparator(text: impl Into<String>) -> Error {
        Error::Comparator(text.into())
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Comparator(text) => f.write_str(text),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

/// The panic payload of a run that meets an input no reading of the comparison expects: what
/// went wrong. `handover::run` ends such a run with `ERROR: [compare_vep_formats] <text>` and
/// exit code 1.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Uncaught(pub String);

/// Ends the run on an input no reading of the comparison expects; the arguments are those of
/// `format!`.
#[macro_export]
macro_rules! uncaught {
    ($($arg:tt)*) => {
        std::panic::panic_any($crate::error::Uncaught(format!($($arg)*)))
    };
}
