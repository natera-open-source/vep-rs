// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tool's parts as a library, so that the tests reach them.
//!
//! The reader of Perl Storable images and what Perl does with a retrieved value are the
//! modules of `tools/vep-storable-to-json`, compiled here from their own files. Over
//! them: the transcript shards of a cache directory, each read as `gunzip` and
//! `fd_retrieve` give it (`shards`), the text `Data::Dumper` returns for a value
//! (`dumper`), Perl strings joined and printed (`text`), MD5 (`md5`), the two digests
//! of a transcript (`digest`) and the run (`run`).

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

pub mod digest;
pub mod dumper;
pub mod md5;
pub mod run;
pub mod shards;
pub mod text;
