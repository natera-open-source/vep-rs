// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Made-up conservation files for the LoFTEE plugin, written from a converted Ensembl VEP
//! cache so that vep-rs and Ensembl VEP are compared on the same files: a tabix-indexed
//! table with one score for each base of every exon, and a SQLite database with one total
//! and one pair of PhyloCSF scores for each exon of every transcript. Every value is made
//! up, a function of a seed and the row's key. [`values`] holds those functions, [`exons`]
//! reads a sequence's transcripts, [`table`] and [`database`] write the two files and
//! [`generate`] runs them over a cache.
//!
//! What the plugin reads of each file is cited as `LoF.pm:line` and `gerp_dist.pl:line`,
//! the lines of the plugin's `master` branch.

pub mod database;
pub mod exons;
pub mod generate;
pub mod table;
pub mod values;

/// The name of the table of base scores.
pub const TABLE_FILE: &str = "GERP_scores.synthetic.txt.gz";

/// The name of the table's tabix index.
pub const INDEX_FILE: &str = "GERP_scores.synthetic.txt.gz.tbi";

/// The name of the database of exon scores.
pub const DATABASE_FILE: &str = "phylocsf_gerp.synthetic.sql";
