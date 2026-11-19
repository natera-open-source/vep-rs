// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The library behind the vep command of vep-rs, a Rust reimplementation of the Ensembl Variant Effect Predictor (VEP).
//!
//! # Public API surface
//!
//! The `pub use`s below are the library surface for downstream consumers,
//! including a long-lived process that holds an `Annotator` and annotates many
//! batches from one load. The contract they carry is stated on
//! `docs/developers/embedding.md`: how the crates are added, what each
//! re-export is for, the path that writes the `vep --json` bytes, and how the
//! surface is versioned; `crates/vep-cli/tests/library_embedding.rs` proves
//! the bytes. Code that pulls from one of these re-exports is opting in to
//! that contract; code that reaches into the submodules below relies on
//! implementation details that may move between releases.
//!
//! The `pub` submodules let the in-tree binaries (`vep`, and `bench_plugins` under the `bench`
//! feature), integration tests, and downstream crates reach into them; the rest are `pub(crate)`.

pub use annotator::{
    build_thread_pool, row_to_input_variant, validate_raw_row, AnnotationError,
    AnnotationResources, AnnotationStats, Annotator, PredictionConfig, RawRow,
};
pub use config::Config;
pub use pick::{apply_to_batch, apply_to_record, FilterConfig, PickCriterion};
pub use row_filter::RowFilters;
pub use transcript_index::{LazyTranscriptIndexes, MaterializedTranscriptIndexes, TranscriptSet};

// The per-variant data type and core consequence types are re-exported so
// callers can construct and inspect annotations without a direct `vep-core`
// dependency.
pub use vep_core::consequence::{Consequence, Impact, TranscriptConsequence};
pub use vep_core::variant::InputVariant;

// The JSON formatter is re-exported from `vep-io` so embedding consumers need
// no separate dependency.
pub use vep_io::output::json::JsonOutputFormatter;

pub mod affinity;
pub mod annotator;
pub mod argv;
pub mod cache_resolve;
pub mod config;
pub(crate) mod delimiter;
pub mod dylib_plugins;
pub mod filter_cli;
pub mod filter_input;
pub(crate) mod filter_regex;
pub mod filter_set;
pub mod frequency_filter;
pub(crate) mod glob;
pub mod individuals;
pub mod input_shape;
pub mod json_cache;
pub mod nearest;
pub mod option_rules;
pub mod output_compress;
pub mod pick;
pub(crate) mod pipeline;
pub mod recode;
pub mod recode_cli;
pub mod regulatory;
pub(crate) mod render;
pub mod row_filter;
pub mod runner;
pub(crate) mod stats;
pub mod stderr;
pub mod sv_source;
pub mod transcript_index;
pub mod transcript_lookup;
pub(crate) mod validate;
pub mod variation_matcher;
pub mod variation_name_index;
pub mod variation_source;
pub mod vcf_parser;
pub mod warnings;

pub mod args;
