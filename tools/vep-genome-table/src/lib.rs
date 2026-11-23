// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The genome table: one row a ready-made Ensembl VEP cache, saying whether vep-rs and
//! Ensembl VEP annotate that genome's generated records alike. This library reads, checks
//! and rewrites a table directory and computes the digest a row carries;
//! `docs/genome-checks.md` describes the files.

mod digest;
mod reader;
mod row;
mod table;

pub use digest::{digest_file, digest_reader};
pub use row::{
    genome_key, header, is_bacterial_set, is_set, outcome_kind, OutcomeKind, Row, COLLECTION_FORMS,
    COLUMNS, NAMED_SETS, TRANSCRIPTS,
};
pub use table::{
    check, load, put, Failure, Put, Summary, BACTERIA_FILE, GENOMES_FILE, OUTCOMES_COLUMNS,
    OUTCOMES_FILE, SETS_COLUMNS, SETS_FILE,
};
