// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tool's parts as a library, so that the tests reach them.
//!
//! The reader of Perl Storable images, Perl's handling of a retrieved value, Storable's
//! outcomes for an image and the worker threads are the modules of
//! `tools/vep-storable-to-json`; a retrieved shard and the entries of its hashes are
//! the modules of `tools/vep-storable-digest`. All are compiled here from their own
//! files. Over them: the walk `File::Find` makes of a directory (`find`), Perl's sort
//! and random numbers (`perlsort`), JSON as `JSON::PP` reads and writes it (`jsonpp`),
//! the field paths of a transcript (`fields`) and the audit (`audit`).

#[allow(dead_code, unused_imports)]
#[path = "../../vep-storable-to-json/src/storable/mod.rs"]
pub mod storable;

#[path = "../../vep-storable-to-json/src/json.rs"]
pub mod json;
#[path = "../../vep-storable-to-json/src/number.rs"]
pub mod number;
#[path = "../../vep-storable-to-json/src/perl.rs"]
pub mod perl;
#[path = "../../vep-storable-to-json/src/pool.rs"]
pub mod pool;
#[path = "../../vep-storable-to-json/src/shard.rs"]
pub mod shard;

#[path = "../../vep-storable-digest/src/dumper.rs"]
pub mod dumper;
#[path = "../../vep-storable-digest/src/shards.rs"]
pub mod shards;

pub mod audit;
pub mod fields;
pub mod find;
pub mod jsonpp;
pub mod perlsort;
