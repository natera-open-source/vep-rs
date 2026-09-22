// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Input parsers and output formatters for VEP.
//!
//! This crate provides:
//! - Input parsers for VCF (via noodles, with coordinate conversion to VEP
//!   format), Ensembl default, HGVS genomic and region formats
//! - Default, tab, VCF, JSON and Parquet-intermediate output formatters

pub mod error;
pub mod input;
pub mod output;

pub use error::IoError;
pub use input::InputParser;
pub use output::OutputFormatter;
