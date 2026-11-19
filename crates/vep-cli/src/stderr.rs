// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Standard error as the `vep` run writes it: the banner, the status, progress
//! and summary lines, the echo of a warning and the `ERROR:` line a refusal
//! ends with all pass through this module. The `vep` binary calls
//! [`drop_refused_writes`] before its first line, and from then on a write the
//! stream refuses is dropped and the run goes on, as the `warn` of Ensembl
//! VEP's `warning_msg` (BaseVEP.pm:821) fails without ending its run on a
//! device that refuses every write: a run whose output is whole exits 0, and a
//! refusal keeps its exit status whether or not its line was written. A process
//! that has not called it, a program embedding the library or a test, writes
//! through [`write()`] as `eprint!` does: a test harness captures the text, and
//! a refused stream panics.

use std::fmt;
use std::io::{self, Write};
use std::sync::atomic::{AtomicBool, Ordering};

/// Whether [`write()`] drops a text standard error refuses.
static DROP_REFUSED: AtomicBool = AtomicBool::new(false);

/// Makes every later [`write()`] of the process drop a text that standard
/// error refuses.
pub fn drop_refused_writes() {
    DROP_REFUSED.store(true, Ordering::Relaxed);
}

/// Writes `args` to standard error: as [`write_or_drop`] once the process has
/// called [`drop_refused_writes`], and as `eprint!` until then.
pub fn write(args: fmt::Arguments<'_>) {
    if DROP_REFUSED.load(Ordering::Relaxed) {
        write_or_drop(args);
    } else {
        eprint!("{args}");
    }
}

/// Writes `args` to standard error and drops the text when the stream refuses
/// it, in any process: for a write whose panic would leave the pipeline
/// waiting.
pub fn write_or_drop(args: fmt::Arguments<'_>) {
    let _ = io::stderr().write_fmt(args);
}

/// Standard error for a caller that takes a writer: every write counts as
/// taken whole, a refused one included.
#[derive(Debug, Clone, Copy, Default)]
pub struct Stream;

impl Write for Stream {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let _ = io::stderr().write_all(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
