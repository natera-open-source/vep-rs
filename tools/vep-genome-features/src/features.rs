// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The features file: what `extract` writes as JSON and `gather` reads back. Every list
//! has one order, so the same cache gives the same bytes on any machine.

use serde::{Deserialize, Serialize};

/// The number of the file's shape, written as `features_version`. It changes whenever a
/// key is renamed, removed or counted another way.
pub const FEATURES_VERSION: u32 = 2;

/// The features of one converted cache.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Features {
    /// The shape of this file, [`FEATURES_VERSION`].
    pub features_version: u32,
    /// The genome's key, where the run named one.
    pub genome: Option<String>,
    /// The species of `info.json`.
    pub species: Option<String>,
    /// The assembly of `info.json`.
    pub assembly: Option<String>,
    /// Every transcript of the cache: each record of a sequence counted once, told apart
    /// by its stable id, start, end and `dbID`.
    pub transcripts: u64,
    /// The transcripts that carry a coding sequence.
    pub coding_transcripts: u64,
    /// The coding transcripts by genetic code, the lowest code first.
    pub genetic_codes: Vec<GeneticCode>,
    /// The transcripts by biotype, in bytewise order of the name.
    pub biotypes: Vec<Biotype>,
    /// The coding transcripts whose coding sequence is cut off at an end.
    pub incomplete_coding: IncompleteCoding,
    /// The sequences of the cache.
    pub sequences: Sequences,
    /// The sequence synonyms of `chr_synonyms.txt`.
    pub synonyms: Synonyms,
    /// The cache holds variation files and `info.json` names their columns.
    pub variant_data: bool,
    /// `info.json` names SIFT predictions.
    pub sift_data: bool,
    /// `info.json` declares regulatory features.
    pub regulatory_data: bool,
    /// The transcripts whose exons, in the cache's order, run against position on the
    /// transcript's strand.
    pub exons_against_position: u64,
}

/// The coding transcripts of one genetic code.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GeneticCode {
    /// The code's number; a transcript that names none, or names 0, counts under 1.
    pub code: u8,
    /// The coding transcripts on the code.
    pub coding_transcripts: u64,
    /// The sequences that carry them, in the cache's reading order.
    pub sequences: Vec<SequenceCount>,
    /// The first codon of the coding sequences that start complete, in bytewise order
    /// of the codon.
    pub opening_codons: Vec<OpeningCodon>,
}

/// The coding transcripts of one genetic code on one sequence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceCount {
    pub name: String,
    pub coding_transcripts: u64,
}

/// The coding transcripts of one genetic code whose coding sequence starts complete and
/// opens with one codon.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OpeningCodon {
    /// The three bases, upper case.
    pub codon: String,
    /// The transcripts that open with it.
    pub transcripts: u64,
    /// Of those, the ones whose cached peptide opens with methionine.
    pub peptide_opens_with_methionine: u64,
    /// Of those, the ones whose cached peptide opens with another residue. The rest
    /// carry no cached peptide.
    pub peptide_opens_otherwise: u64,
}

/// The transcripts of one biotype.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Biotype {
    /// The biotype's name, empty for a transcript that names none.
    pub biotype: String,
    pub transcripts: u64,
}

/// The coding transcripts that carry the `cds_start_NF` or the `cds_end_NF` attribute.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IncompleteCoding {
    /// With `cds_start_NF`.
    pub start: u64,
    /// With `cds_end_NF`.
    pub end: u64,
    /// With either attribute, a transcript with both counted once.
    pub transcripts: u64,
    /// `transcripts` over the coding transcripts, to six decimal places; 0 in a cache
    /// without a coding transcript.
    pub share: f64,
}

/// The sequences of the cache.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sequences {
    /// The sequences the cache holds features for: the directories under `transcripts`.
    pub with_features: u64,
    /// Of those, the ones that hold a transcript.
    pub with_transcripts: u64,
    /// What the side file of sequence lengths states, `None` for a cache read
    /// without one.
    pub side_file: Option<SideFile>,
}

/// What the side file `sequence_regions.tsv` states.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SideFile {
    /// The file's rows.
    pub listed: u64,
    /// The longest sequence, the first of the file among equals; `None` for a file
    /// without a row.
    pub longest: Option<Longest>,
    /// The sequences the file marks circular, in the file's order.
    pub circular: Vec<String>,
}

/// The longest sequence of the side file.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Longest {
    pub name: String,
    pub length: u64,
}

/// The sequence synonyms of the cache.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Synonyms {
    /// The cache names at least one synonym.
    pub present: bool,
    /// The pairs of a sequence's name and another name for it.
    pub pairs: u64,
}
