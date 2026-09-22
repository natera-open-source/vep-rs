// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Built-in plugin implementations.
//!
//! Each module implements a specific VEP plugin using the [`BuiltinPlugin`] trait.
//! Tabix-based plugins use the shared [`TabixAnnotator`] engine for batched I/O.

pub mod alpha_missense;
pub mod cadd;
pub mod dbnsfp;
pub mod dbscsnv;
pub mod gnomadc;
pub mod gwas;
pub mod loftee;
pub mod loftool;
pub mod pli;
pub mod revel;
pub mod spliceai;
