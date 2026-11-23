// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The structural-variant comparator: tuple concordance between Ensembl VEP and a scored
//! engine across a directory of VCF files, with the divergence masks of a reference release.
//!
//! `pyfmt`, `input` and `rules` are the source files of `tools/vep-outputs-compare`, compiled
//! into this crate under the names those files use for each other, so that the registry of
//! divergence classes, the pair rules and the Perl-reading helpers have one definition for
//! both comparators. `coding_end` and `genetic_codes` are source files of
//! `tools/vep-formats-compare`, which `rules` asks by those names.

#[allow(dead_code)]
#[path = "../../vep-outputs-compare/src/input.rs"]
pub mod input;
#[allow(dead_code)]
#[path = "../../vep-outputs-compare/src/pyfmt.rs"]
pub mod pyfmt;
#[allow(dead_code)]
#[path = "../../vep-outputs-compare/src/rules/mod.rs"]
pub mod rules;

/// The shapes of the class `alternate_codon_past_the_coding_end` and the genetic codes they
/// are translated under: source files of `tools/vep-formats-compare`, which
/// `rules::coding_ends` asks by these names.
#[path = "../../vep-formats-compare/src/rules/coding_end.rs"]
pub mod coding_end;
#[path = "../../vep-formats-compare/src/genetic_codes.rs"]
pub mod genetic_codes;

pub mod cache;
pub mod cli;
pub mod compare;
pub mod failure;
pub mod locus;
pub mod masks;
pub mod mate_context;
pub mod metrics;
pub mod normalize;
pub mod output;
pub mod perl;
pub mod record;
pub mod registry;
pub mod report;
pub mod strings;
pub mod tempdir;
pub mod text;
pub mod threads;
pub mod vcf;

pub use failure::{Failure, FailureKind};
pub use tempdir::TempDir;
