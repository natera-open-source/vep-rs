// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A seeded generator of input records for any converted Ensembl VEP cache: about thirty
//! records a genome, one for each situation of the genome check
//! (`docs/genome-checks.md`), the same bytes for the same cache and seed on any machine.
//! [`cache`] reads a converted cache, [`select`] holds one function a situation,
//! [`generate`] runs them over a cache and [`vcf`] writes the records.

pub mod cache;
pub mod fasta;
pub mod generate;
pub mod random;
pub mod reservoir;
pub mod select;
pub mod vcf;

/// The number of the generating algorithm. It is written in every output's header and
/// printed by `--version`, and it changes whenever the same cache and seed would give
/// other bytes.
pub const GENERATOR_VERSION: u32 = 5;
