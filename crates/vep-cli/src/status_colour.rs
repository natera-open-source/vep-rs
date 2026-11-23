// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Whether the timestamped lines a binary of this crate logs on standard
//! error carry colour codes: on a terminal alone, so the lines a file or a
//! pipe takes are plain text.

use std::ffi::OsStr;
use std::io::IsTerminal;

/// Whether the lines logged on standard error carry colour codes: only where
/// standard error is a terminal and `NO_COLOR` is unset or empty. The value is
/// for `with_ansi`, which replaces the log layer's own reading of `NO_COLOR`,
/// so the variable is read here.
pub fn on_standard_error() -> bool {
    carried(
        std::io::stderr().is_terminal(),
        std::env::var_os("NO_COLOR").as_deref(),
    )
}

/// [`on_standard_error`] for a standard error that is a `terminal` or is not,
/// and for what `NO_COLOR` holds.
fn carried(terminal: bool, no_color: Option<&OsStr>) -> bool {
    terminal && no_color.is_none_or(OsStr::is_empty)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Colour codes go to a terminal alone, and not there once `NO_COLOR` holds
    /// a value; a file or a pipe takes none whatever the variable holds.
    #[test]
    fn colour_is_carried_on_a_terminal_alone_and_not_under_no_color() {
        let set = OsStr::new("1");
        let empty = OsStr::new("");
        assert!(carried(true, None));
        assert!(carried(true, Some(empty)));
        assert!(!carried(true, Some(set)));
        for no_color in [None, Some(empty), Some(set)] {
            assert!(!carried(false, no_color));
        }
    }
}
