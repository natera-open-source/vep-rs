// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the `score` program hands the settling step and takes back from it.

pub use super::class_table::{Classes, Declared, Reading, TableFacts, NOT_MET};
pub use super::error::{Scored, ScorerError};
pub use super::residuals::{residual_path, Residuals, Side};
pub use super::settle::{score_engine, source_entries, SourceEntry, Spread};
pub use super::sources::{text, without_chr, SourceOf, Sources, RESOLVED_LIMIT};
pub use super::tally::Tally;
use crate::input::Kind;
pub use crate::registry::Registry;
pub use crate::rules::{BuiltIn, Consult, Reads, Rules};
use std::path::PathBuf;

/// The projections of a run, in the order every per-projection pair is held.
pub const PROJECTIONS: [&str; 2] = ["feature", "allele"];

/// A keyword of a class's consultation that a run hands the classes of a plugin.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Keyword {
    GwasPairs,
    RevelTable,
    UtrannotatorSweep,
    DownstreamSequences,
    /// The applied classes.
    Declared,
    SpliceWindows,
    SpliceaiTable,
    ExonSpans,
}

/// Per registry plugin name, the keywords the classes the table declares for that plugin
/// are handed, each read from the run's `Reads`.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct PluginReads(pub Vec<(String, Vec<Keyword>)>);

/// The scored run's input VCF by identifier.
pub trait NamedLines: Sync {
    /// The input's lines whose ID column has the identifier `name` as its first
    /// `;`-joined element, byte for byte, each as its first five columns in the VCF's
    /// order, tab-separated, with `name` as the third.
    fn named(&self, name: &[u8]) -> Vec<String>;
}

/// One input after the passes: what they counted of it and, against one engine, the files
/// of its side of the remainder.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Extracted {
    pub ident: String,
    /// The input as the command line names it.
    pub path: PathBuf,
    pub kind: Kind,
    pub sha256: String,
    /// The input's rows, or the lines of a tuple file.
    pub lines: u64,
    /// A tuple file's tuples whose field no source owns.
    pub outside: u64,
    /// The keys no source claims as text, in their bytes' order, each with its count.
    pub unclaimed: Vec<(String, u64)>,
    /// The column header a default-format input's rows are read under.
    pub header: String,
    /// The side's feature tuples of the remainder, sorted.
    pub feature: PathBuf,
    /// The side's allele tuples of the remainder, sorted.
    pub allele: PathBuf,
    /// The side's rows at the Locations of the remainder, sorted; a tuple file's tuples
    /// under the remainder's keys.
    pub rows: PathBuf,
}

/// What the passes counted of one projection of one engine.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct CountedProjection {
    /// Per field name that either side wrote, the reference's values, the engine's values
    /// and the pairs equal as text.
    pub values: Vec<(Vec<u8>, [u64; 3])>,
    /// Per source, the keys under which a side writes one field more than once.
    pub shared_keys: Vec<(String, u64)>,
}

/// What the passes counted of one engine against the reference.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Counted {
    pub feature: CountedProjection,
    pub allele: CountedProjection,
}
