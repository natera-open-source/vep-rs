// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Library entry point for vep-cli.
//!
//! # Public API surface
//!
//! The `pub use`s below are the **intended** library surface for downstream
//! consumers, including a long-lived service that holds an `Annotator` directly
//! and serves many requests from one load. Code that pulls from one of these
//! re-exports is opting in to a stable contract; code that reaches into the
//! submodules below relies on implementation details that may move between
//! releases.
//!
//! Submodules remain `pub` so the in-tree binaries (`vep`, `bench_plugins`),
//! integration tests, and downstream crates can keep reaching into them.

pub use annotator::{
    build_thread_pool, row_to_input_variant, validate_raw_row, AnnotationError,
    AnnotationResources, AnnotationStats, Annotator, PredictionConfig, RawRow,
};
pub use config::Config;
pub use pick::{apply_filters, FilterConfig};
pub use transcript_index::{
    LazyTranscriptIndexes, MaterializedTranscriptIndexes, TranscriptIndexImpl,
};

// The per-variant data type and core consequence types are re-exported so
// callers can construct and inspect annotations without a direct `vep-core`
// dependency.
pub use vep_core::consequence::{Consequence, Impact, TranscriptConsequence};
pub use vep_core::variant::InputVariant;

// The JSON formatter is re-exported from `vep-io` so embedding consumers need
// no separate dependency.
pub use vep_io::output::json::JsonOutputFormatter;

pub mod annotator;
pub mod config;
pub mod dylib_plugins;
pub mod json_cache;
pub mod pick;
pub(crate) mod pipeline;
pub mod runner;
pub mod transcript_index;
pub mod variation_matcher;
pub mod vcf_parser;

pub mod args;
