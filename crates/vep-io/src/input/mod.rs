// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Input parsers for variant formats.

pub mod ensembl;
pub mod hgvs;
pub mod region;
pub mod vcf;

use vep_core::variant::InputVariant;

use crate::error::IoError;

/// Trait for parsing input variants from various formats.
pub trait InputParser: Send {
    /// Read the next variant from the input.
    /// Returns `None` at end of input.
    fn next_variant(&mut self) -> Option<Result<InputVariant, IoError>>;
}
