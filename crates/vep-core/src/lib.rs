// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Core types for vep-rs, a Rust implementation of Ensembl VEP.
//!
//! This crate provides the foundational data types used throughout the VEP
//! Rust implementation: variants, transcripts, consequences, genomic coordinates,
//! codon translation, and allele utilities.

pub mod allele;
#[doc(hidden)]
pub mod allele_trim;
pub mod assembly;
pub mod binding_matrix;
pub mod cache_info;
pub mod cache_manifest;
pub mod codon;
pub mod consequence;
pub mod coordinate;
pub mod custom_annotation;
pub mod minimal;
pub mod pick;
pub mod prediction;
pub mod regulatory_gff;
pub mod so_ontology;
pub mod spdi;
pub mod stderr;
pub mod sv_record;
pub mod sv_routing;
pub mod sv_span;
pub mod synonyms;
pub mod terms;
pub mod transcript;
pub mod variant;
pub mod variation;
pub mod vcf_header;
pub mod vcf_info;
pub mod warning_hold;

/// VEP version number (matches the Ensembl VEP release).
pub const VEP_VERSION: u32 = 116;

/// VEP sub-version number.
pub const VEP_SUB_VERSION: u32 = 3;

/// The Ensembl API modules of the Ensembl VEP release that [`VEP_VERSION`] and
/// [`VEP_SUB_VERSION`] name, each with the version that release's container
/// image reports for it, `<release>.<the first seven characters of its commit>`,
/// sorted by module name as Ensembl VEP lists them (`OutputFactory/VCF.pm`
/// lines 244-256, `Utils.pm` lines 506-531). They state the reference release's
/// install: no module of the list runs in vep-rs.
pub const VEP_API_MODULE_VERSIONS: &[(&str, &str)] = &[
    ("ensembl", "116.0d85231"),
    ("ensembl-compara", "116.4038708"),
    ("ensembl-funcgen", "116.90049ea"),
    ("ensembl-io", "116.6afb5dc"),
    ("ensembl-variation", "116.db7394a"),
];
