// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The tool's parts as a library, so that the tests reach them.
//!
//! The reader of Perl Storable images and what Perl does with a retrieved value are the
//! modules of `tools/vep-storable-to-json`; the transcript shards of a cache directory,
//! the text `Data::Dumper` returns for a value and Perl strings joined and printed are
//! the modules of `tools/vep-storable-digest`. All are compiled here from their own
//! files. Over them: the dump of one transcript and the run (`run`).

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
#[path = "../../vep-storable-digest/src/text.rs"]
pub mod text;

pub mod run;
