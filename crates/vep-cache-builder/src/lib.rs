// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The library of the cache builder serves `vep`, which reads a gene-model
//! file named beside its cache through [`gene_model`]; the other modules are
//! the builder's own and are hidden from the documentation.
//!
//! Among them are the readers of the cache builder: a GFF3 or GTF transcript
//! model read as Ensembl VEP's `--gff` and `--gtf` sources read one
//! (`gff3`, `gtf`) and the transcript records built from it (`builder`),
//! with the Ensembl inputs a build of the cache's content adds
//! (`core_dumps`, `fasta`), the regulatory tree baked from the
//! regulatory GFF files (`regulatory`), the variation files (`variation`)
//! with the cells a variation database's tables carry (`variant_cells`)
//! and the indexed genome a cache directory carries (`fasta_index`).
//!
//! The `vep-cache-builder` binary writes those records as a cache directory.
//! `vep` reads the same records in memory for a gene-model file named beside
//! its cache ([`gene_model`]), so one reader serves both.

#[doc(hidden)]
pub mod builder;
#[doc(hidden)]
pub mod core_dumps;
#[doc(hidden)]
pub mod count;
#[doc(hidden)]
pub mod dump_rows;
#[doc(hidden)]
pub mod fasta;
#[doc(hidden)]
pub mod fasta_index;
pub mod gene_model;
#[doc(hidden)]
pub mod gff3;
#[doc(hidden)]
pub mod gtf;
mod gvf;
#[doc(hidden)]
pub mod info;
mod line;
mod mapper;
#[doc(hidden)]
pub mod regulatory;
mod sorted_runs;
mod translation_row;
#[doc(hidden)]
pub mod variant_cells;
#[doc(hidden)]
pub mod variation;
