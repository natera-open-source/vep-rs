// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The library behind the vep command of vep-rs, a Rust reimplementation of the Ensembl Variant Effect Predictor (VEP).
//!
//! # Public API surface
//!
//! The interface of this library is what this documentation shows: the items
//! re-exported at the crate root, the modules listed below, and the types
//! their signatures and public fields name. `docs/developers/embedding.md`
//! states the contract a program that links the library relies on: how the
//! crates are added, what each item is for, the path that writes the
//! `vep --json` bytes, and how the surface is versioned;
//! `crates/vep-cli/tests/library_embedding.rs` proves the bytes.
//!
//! - The settings are a [`Config`]. [`Config::from_argv`] resolves a command
//!   line as the `vep` binary does, and `Config::default()` is the
//!   configuration of a command line that gives no flag, with nothing read
//!   from the host. [`prepass`] reads the `--config` file, `vep.ini` and the
//!   `VEP_*` variables into the command line `from_argv` takes, and
//!   [`resolve_fork`] gives the thread count a `--fork` value resolves to.
//! - An [`Annotator`] holds one loaded cache and annotates many batches from
//!   it. [`Annotator::from_config`] builds it from the settings;
//!   [`Annotator::new`] builds it from [`AnnotationResources`] a program
//!   assembles itself, with [`build_thread_pool`], [`PredictionConfig`],
//!   [`LazyTranscriptIndexes`], [`TranscriptSet`], [`VariationSource`] and
//!   [`ColocationRules`].
//! - The input is a [`RawRow`], checked by [`validate_raw_row`] and read by
//!   [`row_to_input_variant`], or a VCF data line read by [`parse_vcf_line`]
//!   and [`parse_vcf_line_with`] into [`InputVariant`] values.
//! - [`apply_to_record`] and [`apply_to_batch`] apply a [`FilterConfig`] to
//!   annotated records: the row filters of [`RowFilters`], then the pick
//!   flags in the order of its [`PickCriterion`] list.
//! - [`JsonOutputFormatter`] writes a record as `vep --json` writes it.
//! - [`Runner`] is the whole run of the `vep` binary, from the input file to
//!   the output file.
//!
//! The five items the pages name by a module path are the same items at that
//! path: `vep_cli::argv::prepass`, `vep_cli::config::resolve_fork`,
//! `vep_cli::runner::Runner`, `vep_cli::vcf_parser::parse_vcf_line` and
//! `vep_cli::vcf_parser::parse_vcf_line_with`.
//!
//! An item this documentation does not show is the binary's own and may
//! change in any release. Where such an item is `pub`, it is so for the `vep`
//! and `bench_plugins` binaries, the tests and the cache tools of this
//! workspace.

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

// The items the pages name by a module path, at the root as well.
pub use argv::prepass;
pub use config::resolve_fork;
pub use runner::Runner;
pub use vcf_parser::{parse_vcf_line, parse_vcf_line_with};

// The types a documented signature or field names whose module is hidden.
pub use cache_resolve::SourceType;
pub use nearest::NearestType;
pub use variation_matcher::ColocationRules;
pub use variation_source::VariationSource;

#[doc(hidden)]
pub mod affinity;
pub mod annotator;
pub mod argv;
pub(crate) mod cache_left_out;
#[doc(hidden)]
pub mod cache_resolve;
pub mod config;
pub(crate) mod delimiter;
pub(crate) mod dylib_plugins;
#[doc(hidden)]
pub mod filter_cli;
pub(crate) mod filter_input;
pub(crate) mod filter_regex;
pub(crate) mod filter_set;
pub(crate) mod frequency_filter;
pub(crate) mod gene_model;
pub(crate) mod glob;
pub(crate) mod individuals;
#[doc(hidden)]
pub mod input_shape;
#[doc(hidden)]
pub mod json_cache;
#[doc(hidden)]
pub mod nearest;
#[doc(hidden)]
pub mod option_rules;
pub(crate) mod output_compress;
pub mod pick;
pub(crate) mod pipeline;
pub(crate) mod recode;
#[doc(hidden)]
pub mod recode_cli;
#[doc(hidden)]
pub mod regulatory;
pub(crate) mod render;
pub mod row_filter;
pub mod runner;
pub(crate) mod stats;
#[doc(hidden)]
pub mod status_colour;
#[doc(hidden)]
pub mod stderr;
pub(crate) mod sv_source;
pub mod transcript_index;
pub(crate) mod transcript_lookup;
pub(crate) mod validate;
#[doc(hidden)]
pub mod variation_matcher;
pub(crate) mod variation_name_index;
#[doc(hidden)]
pub mod variation_source;
pub mod vcf_parser;
#[doc(hidden)]
pub mod warnings;

#[doc(hidden)]
pub mod args;
