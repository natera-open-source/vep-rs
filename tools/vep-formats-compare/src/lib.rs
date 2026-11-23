// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Per-column comparison of a reference and a candidate VEP output.

pub mod argparse;
pub mod cache;
pub mod cache_classes;
pub mod classes;
pub mod cli;
pub mod coding_ends;
pub mod compare;
pub mod compare_files;
pub mod consts;
pub mod counter;
pub mod entry;
pub mod error;
pub mod fasta;
pub mod field_class;
pub mod fields;
pub mod floors;
pub mod genetic_codes;
pub mod handover;
pub mod input;
pub mod lists;
pub mod manifest;
pub mod nearest;
pub mod parse;
pub mod pattern;
pub mod pyfmt;
pub mod pyjson;
pub mod record_rows;
pub mod regulatory;
pub mod remainder;
pub mod report;
pub mod rules;
pub mod self_parity;
pub mod shared;
pub mod stream;
pub mod warnings;

/// The shapes of the class `alternate_codon_past_the_coding_end`, under the name by which
/// `coding_ends::terms` asks them in every tool root that compiles that file.
pub use rules::coding_end;

/// The component the error lines name.
pub const COMPONENT: &str = "compare_vep_formats";
