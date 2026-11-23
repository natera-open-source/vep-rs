// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Built-in VEP plugins compiled directly into the binary.
//!
//! This crate provides:
//! - `BuiltinPlugin` trait for native plugin implementations
//! - [`TabixAnnotator`] batch-oriented tabix query engine
//! - [`BinaryAnnotator`] mmap-backed binary annotation store
//! - [`AnnotationStore`] trait unifying tabix and binary backends
//! - `BuiltinRegistry` for name-based plugin lookup and lifecycle management
//! - `PluginParams` parser for the `--plugin Name,param1,param2` parameters
//! - [`perl_format`] for scores printed as Perl's `sprintf` prints them
//! - `bigwig`: the opener every bigWig-backed plugin and source goes through
//! - Concrete plugin implementations for common clinical annotations
//! - [`custom`]: the `--custom` argument grammar and the sources it names
//!
//! An item named above without a link is for the workspace's own crates: it
//! is public and hidden from the API documentation.

pub mod annotation_store;
#[doc(hidden)]
pub mod bigwig;
pub mod binary_store;
pub mod custom;
#[doc(hidden)]
pub mod params;
pub mod perl_format;
mod pool_pass;
#[doc(hidden)]
pub mod registry;
mod sqlite;
pub mod tabix;
pub mod traits;

#[cfg(test)]
mod test_fixtures;

mod plugins;

pub use annotation_store::AnnotationStore;
pub use binary_store::BinaryAnnotator;
#[doc(hidden)]
pub use params::{Grammar, PluginParams};
#[doc(hidden)]
pub use registry::{
    unported_reason, BatchOutcome, BuiltinRegistry, PluginFailure, PluginPrefetch, PluginTimings,
    UNPORTED_ENSEMBL_PLUGINS,
};
pub use tabix::TabixAnnotator;
pub use traits::PluginError;
#[doc(hidden)]
pub use traits::{BuiltinPlugin, OutputFormat, PrefetchData, RunContext};
