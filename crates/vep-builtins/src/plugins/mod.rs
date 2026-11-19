// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Built-in plugin implementations.
//!
//! Each module implements a specific VEP plugin using the [`crate::BuiltinPlugin`] trait.
//! Tabix-based plugins use the shared [`crate::TabixAnnotator`] engine for batched I/O.

pub mod alpha_missense;
pub mod ancestral_allele;
pub mod avi;
pub mod bayesdel;
pub mod blosum62;
pub mod cadd;
pub mod clinpred;
pub mod conservation;
pub mod dbnsfp;
pub mod dbscsnv;
pub mod dosage_sensitivity;
pub mod downstream;
pub mod gene_scores;
pub mod gnomadc;
pub mod gnomadmt;
pub mod gwas;
pub mod hgvs_intron_offset;
pub mod loeuf;
pub mod loftee;
pub mod loftool;
pub mod maxentscan;
pub mod nearest_exon_jb;
pub mod nmd;
pub mod open_targets;
pub mod pli;
pub mod primate_ai;
pub mod revel;
pub mod spliceai;
pub mod splicevault;
pub mod transcript_annotator;
pub mod tss_distance;
pub mod utr_annotator;
