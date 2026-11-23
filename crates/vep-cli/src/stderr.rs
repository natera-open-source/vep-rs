// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Standard error as the `vep` run writes it: the banner, the status, progress
//! and summary lines, the echo of a warning and the `ERROR:` line a refusal
//! ends with all pass through this module, as do the warning of `vep filter`
//! and the lines of `vep recode`; `vep filter` writes its `ERROR:` line
//! itself, as the bytes of the input it quotes, and drops a refused one the
//! same way. Its items are those of [`vep_core::stderr`], whose one switch
//! a built-in plugin's own line and the cache tools read too. The `vep` binary
//! calls [`drop_refused_writes`] before its first line, and from then on a
//! write the stream refuses is dropped and the run goes on, as the `warn` of
//! Ensembl VEP's `warning_msg` (BaseVEP.pm:821) fails without ending its run on
//! a device that refuses every write: a run whose output is whole exits 0, and
//! a refusal keeps its exit status whether or not its line was written. A
//! process that has not called it, a program embedding the library or a test,
//! writes through [`write()`] as `eprint!` does: a test harness captures the
//! text, and a refused stream panics.

pub use vep_core::stderr::{drop_refused_writes, write, write_or_drop, Stream};
pub(crate) use vep_core::stderr::{progress, progress_or_drop};
