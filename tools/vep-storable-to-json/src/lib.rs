// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tool's parts as a library, so that the tests reach them: the reader of Perl
//! Storable images under `storable`, and over it the conversion's rules for a
//! transcript shard, a regulatory shard, `info.json` and `sequence_regions.tsv`, the
//! JSON text of its two encoders, and the run that orders them.

pub mod gzip;
pub mod info;
pub mod json;
pub mod number;
pub mod perl;
pub mod pool;
pub mod regions;
pub mod regulatory;
pub mod run;
pub mod shard;
pub mod storable;
pub mod transcripts;
