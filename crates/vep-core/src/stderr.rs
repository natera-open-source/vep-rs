// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Standard error as every program of the workspace writes it. A binary calls
//! [`drop_refused_writes`] before its first line, and from then on a write the
//! stream refuses, on a device that refuses every write or a pipe whose reader
//! has gone, is dropped: the run goes on and ends with its own exit status,
//! whether or not a line was written. A process that has not called it, a
//! program embedding a library crate or a test, writes through [`write()`] as
//! `eprint!` does: a test harness captures the text, and a refused stream
//! panics.
//!
//! A progress line ([`progress`], [`progress_or_drop`]) is written after a
//! carriage return and has no line end, so on a terminal the next progress
//! line overwrites it. Every other text written through this module ends a
//! progress line that is still open with a line end first, and so opens a line
//! of its own.

use std::fmt;
use std::io::{self, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, MutexGuard, PoisonError};

/// Whether [`write()`] drops a text standard error refuses.
static DROP_REFUSED: AtomicBool = AtomicBool::new(false);

/// Whether the last text written through this module was a progress line,
/// which has no line end. A write holds the lock from its read of the flag to
/// the end of its text, so no other thread's text comes between the two.
static PROGRESS_OPEN: Mutex<bool> = Mutex::new(false);

/// What a text on standard error is.
#[derive(Clone, Copy)]
enum Kind {
    /// A progress line, which the next progress line overwrites.
    Progress,
    /// Any other text.
    Other,
}

/// Makes every later [`write()`] of the process drop a text that standard
/// error refuses.
pub fn drop_refused_writes() {
    DROP_REFUSED.store(true, Ordering::Relaxed);
}

/// Writes `args` to standard error: as [`write_or_drop`] once the process has
/// called [`drop_refused_writes`], and as `eprint!` until then. A progress line
/// that is still open is ended first.
pub fn write(args: fmt::Arguments<'_>) {
    emit(Kind::Other, args, DROP_REFUSED.load(Ordering::Relaxed));
}

/// Writes `args` to standard error and drops the text when the stream refuses
/// it, in any process: for a write whose panic would end the run. A progress
/// line that is still open is ended first.
pub fn write_or_drop(args: fmt::Arguments<'_>) {
    emit(Kind::Other, args, true);
}

/// Writes `args` to standard error as a progress line, after a carriage return
/// and with no line end. A text the stream refuses is treated as [`write()`]
/// treats one.
pub fn progress(args: fmt::Arguments<'_>) {
    emit(Kind::Progress, args, DROP_REFUSED.load(Ordering::Relaxed));
}

/// Writes `args` to standard error as [`progress`] does and drops the text
/// when the stream refuses it, in any process.
pub fn progress_or_drop(args: fmt::Arguments<'_>) {
    emit(Kind::Progress, args, true);
}

/// Writes `args` after what [`lead`] puts before a text of `kind`; an empty
/// text writes nothing and leaves an open progress line open. `args` is
/// formatted before the lock is taken, so a value whose formatting writes a
/// line through this module does not wait on this write.
fn emit(kind: Kind, args: fmt::Arguments<'_>, drop_refused: bool) {
    let text = fmt::format(args);
    if text.is_empty() {
        return;
    }
    let mut open = progress_open();
    let lead = lead(&mut open, kind);
    if drop_refused {
        let _ = io::stderr().write_fmt(format_args!("{lead}{text}"));
    } else {
        eprint!("{lead}{text}");
    }
}

/// The lock on [`PROGRESS_OPEN`]; a thread that panicked on a refused write
/// while it held the lock leaves the flag readable.
fn progress_open() -> MutexGuard<'static, bool> {
    PROGRESS_OPEN.lock().unwrap_or_else(PoisonError::into_inner)
}

/// What stands before a text of `kind` when a progress line is `open`, which
/// it is afterwards only where the text is one: a carriage return before a
/// progress line, the line end that closes an open progress line before any
/// other text, and nothing before a text that follows a line end.
fn lead(open: &mut bool, kind: Kind) -> &'static str {
    match kind {
        Kind::Progress => {
            *open = true;
            "\r"
        }
        Kind::Other if std::mem::take(open) => "\n",
        Kind::Other => "",
    }
}

/// Standard error for a caller that takes a writer, a log subscriber among
/// them: every write counts as taken whole, a refused one included, and ends a
/// progress line that is still open first.
#[derive(Debug, Clone, Copy, Default)]
pub struct Stream;

impl Write for Stream {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if !buf.is_empty() {
            let mut open = progress_open();
            let mut stderr = io::stderr().lock();
            let _ = stderr
                .write_all(lead(&mut open, Kind::Other).as_bytes())
                .and_then(|()| stderr.write_all(buf));
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{lead, Kind};

    /// A progress line opens with a carriage return whatever stands before it.
    /// The first other text after one opens with the line end that closes it,
    /// and a text after that opens with nothing.
    #[test]
    fn a_line_end_closes_an_open_progress_line_before_any_other_text() {
        let mut open = false;
        let texts = [
            Kind::Other,
            Kind::Progress,
            Kind::Progress,
            Kind::Other,
            Kind::Other,
            Kind::Progress,
            Kind::Other,
        ];
        let leads = texts.map(|kind| lead(&mut open, kind));
        assert_eq!(leads, ["", "\r", "\r", "\n", "", "\r", "\n"]);
        assert!(!open);
    }
}
