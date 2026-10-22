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
pub mod cache_manifest;
pub mod codon;
pub mod consequence;
pub mod coordinate;
pub mod minimal;
pub mod prediction;
pub mod so_ontology;
pub mod spdi;
pub mod sv_record;
pub mod sv_routing;
pub mod sv_span;
pub mod synonyms;
pub mod terms;
pub mod transcript;
pub mod variant;
pub mod variation;
pub mod vcf_info;

/// VEP version number (matches the Ensembl VEP release).
pub const VEP_VERSION: u32 = 116;

/// VEP sub-version number.
pub const VEP_SUB_VERSION: u32 = 2;
