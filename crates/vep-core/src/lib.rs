// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Core types for vep-rs, a Rust implementation of Ensembl VEP.
//!
//! This crate provides the foundational data types used throughout the VEP
//! Rust implementation: variants, transcripts, consequences, genomic coordinates,
//! codon translation, and allele utilities.

pub mod allele;
pub mod assembly;
pub mod cache_info;
pub mod codon;
pub mod consequence;
pub mod coordinate;
pub mod prediction;
pub mod synonyms;
pub mod transcript;
pub mod variant;
pub mod variation;

/// VEP version number (matches Perl VEP release).
pub const VEP_VERSION: u32 = 115;

/// VEP sub-version number.
pub const VEP_SUB_VERSION: u32 = 2;

/// Default upstream/downstream distance in base pairs.
pub const DEFAULT_DISTANCE: u64 = 5000;

/// Default buffer size (number of variants per batch).
pub const DEFAULT_BUFFER_SIZE: usize = 5000;

/// Default cache region size in base pairs (1 Mb).
pub const DEFAULT_CACHE_REGION_SIZE: u64 = 1_000_000;
