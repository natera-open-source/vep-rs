// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Built-in VEP plugins compiled directly into the binary.
//!
//! This crate provides:
//! - [`BuiltinPlugin`] trait for native plugin implementations
//! - [`TabixAnnotator`] batch-oriented tabix query engine
//! - [`BinaryAnnotator`] mmap-backed binary annotation store
//! - [`AnnotationStore`] trait unifying tabix and binary backends
//! - [`BuiltinRegistry`] for name-based plugin lookup and lifecycle management
//! - Concrete plugin implementations for common clinical annotations

pub mod annotation_store;
pub mod binary_store;
pub mod registry;
pub mod tabix;
pub mod traits;

#[cfg(test)]
mod test_fixtures;

pub mod plugins;

pub use annotation_store::AnnotationStore;
pub use binary_store::BinaryAnnotator;
pub use registry::{BuiltinRegistry, PluginTimings};
pub use tabix::TabixAnnotator;
pub use traits::{BuiltinPlugin, PluginError, PrefetchData};
