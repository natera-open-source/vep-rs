// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The classes whose rules read the candidate's JSON cache (`--json-cache`).
//!
//! `data/json_cache_classes.tsv` names each class for which a comparison builds a side table
//! from the cache where `--masks` names the class, with that table. The command line tests
//! the masks against the names of this file to decide which tables it builds, and the test of
//! a program outside the tool reads the file by its path to hold that program's lists of the
//! classes it hands the cache. For the classes of the `TRANSCRIPT_READING_CLASSES` column of
//! `data/classes.tsv`, listed here under `downstream_sequences`, the command line builds
//! `partial_codon_models` as well.
//!
//! The structural-variant comparator opens a cache it is given whichever of its classes is
//! named, and refuses a run without one for the classes of the `CACHE_READING_CLASSES` column
//! of `data/classes.tsv`, which are the two this file names for it. The score table of
//! `--spliceai-table` reads the cache on that option alone, whatever the masks name, so its
//! class has no line here.
//!
//! `stop_codon_settled_before_unknown_bases` reads `DownstreamSequences` for one kind of row
//! alone, that of an insertion written in an intron and carried into the coding sequence,
//! which holds no `Codons` cell. Without the table the class reads its rows with a `Codons`
//! cell alone, as it reads them beside one.

use std::sync::LazyLock;

const TABLE: &str = include_str!("../data/json_cache_classes.tsv");

/// A side table a comparison builds from the candidate's JSON cache.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CacheTable {
    /// The transcripts of the structural-variant comparator (`StructuralVariants`).
    StructuralVariants,
    /// `SequenceEdges`, which a comparison builds beside `--reference-fasta` alone.
    SequenceEdges,
    /// `PartialCodonModels`.
    PartialCodonModels,
    /// `LengthEdits`.
    LengthEdits,
    /// `OrganelleCodonTables`.
    CodonTables,
    /// `DownstreamSequences`.
    DownstreamSequences,
    /// `TranscriptEdges`.
    TranscriptEdges,
    /// `IntronEnds`.
    IntronEnds,
    /// `SpliceWindows`.
    SpliceWindows,
}

impl CacheTable {
    /// Every table, in the order a comparison builds them.
    pub const ALL: [CacheTable; 9] = [
        CacheTable::StructuralVariants,
        CacheTable::SequenceEdges,
        CacheTable::PartialCodonModels,
        CacheTable::LengthEdits,
        CacheTable::CodonTables,
        CacheTable::DownstreamSequences,
        CacheTable::TranscriptEdges,
        CacheTable::IntronEnds,
        CacheTable::SpliceWindows,
    ];

    /// The table's name in the `table` column of the file.
    pub fn name(self) -> &'static str {
        match self {
            CacheTable::StructuralVariants => "structural_variants",
            CacheTable::SequenceEdges => "sequence_edges",
            CacheTable::PartialCodonModels => "partial_codon_models",
            CacheTable::LengthEdits => "length_edits",
            CacheTable::CodonTables => "codon_tables",
            CacheTable::DownstreamSequences => "downstream_sequences",
            CacheTable::TranscriptEdges => "transcript_edges",
            CacheTable::IntronEnds => "intron_ends",
            CacheTable::SpliceWindows => "splice_windows",
        }
    }

    /// The table the `table` column names as `name`; `None` for any other text.
    fn named(name: &str) -> Option<CacheTable> {
        CacheTable::ALL
            .into_iter()
            .find(|table| table.name() == name)
    }
}

/// The classes of the file with their tables, and the classes of each table.
struct CacheClasses {
    rows: Vec<(&'static str, CacheTable)>,
    by_table: Vec<Vec<&'static str>>,
}

impl CacheClasses {
    fn read() -> CacheClasses {
        let mut lines = TABLE.lines();
        assert_eq!(
            lines.next(),
            Some("class\ttable"),
            "data/json_cache_classes.tsv: the header names the columns class and table"
        );
        let mut rows: Vec<(&'static str, CacheTable)> = Vec::new();
        for line in lines {
            let cells: Vec<&'static str> = line.split('\t').collect();
            let [class, table] = cells[..] else {
                panic!("data/json_cache_classes.tsv: a line of other than two cells: {line}")
            };
            let table = CacheTable::named(table).unwrap_or_else(|| {
                panic!("data/json_cache_classes.tsv: {class} names no table: {table}")
            });
            assert!(
                !class.is_empty() && rows.iter().all(|(held, _)| *held != class),
                "data/json_cache_classes.tsv: a class has one line: {class}"
            );
            rows.push((class, table));
        }
        let by_table = CacheTable::ALL
            .into_iter()
            .map(|table| {
                rows.iter()
                    .filter(|(_, held)| *held == table)
                    .map(|(class, _)| *class)
                    .collect()
            })
            .collect();
        CacheClasses { rows, by_table }
    }
}

fn cache_classes() -> &'static CacheClasses {
    static CLASSES: LazyLock<CacheClasses> = LazyLock::new(CacheClasses::read);
    &CLASSES
}

/// Every class of the file with the table a comparison builds for it, in the file's order.
pub fn classes() -> &'static [(&'static str, CacheTable)] {
    &cache_classes().rows
}

/// The classes for which a comparison builds `table`, in the file's order.
pub fn classes_reading(table: CacheTable) -> &'static [&'static str] {
    let place = CacheTable::ALL
        .into_iter()
        .position(|held| held == table)
        .unwrap_or_default();
    &cache_classes().by_table[place]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classes::tables;

    #[test]
    fn every_line_names_a_class_of_the_class_file_and_every_table_has_a_class() {
        let class_file = include_str!("../data/classes.tsv");
        let known: Vec<&str> = class_file
            .lines()
            .skip(1)
            .filter_map(|line| line.split('\t').next())
            .collect();
        for (class, _) in classes() {
            assert!(
                known.contains(class),
                "{class} is no class of data/classes.tsv"
            );
        }
        for table in CacheTable::ALL {
            assert!(!classes_reading(table).is_empty(), "{}", table.name());
        }
        let listed: usize = CacheTable::ALL
            .into_iter()
            .map(|table| classes_reading(table).len())
            .sum();
        assert_eq!(listed, classes().len());
    }

    #[test]
    fn the_tables_that_are_columns_of_the_class_file_hold_the_same_classes() {
        let tables = tables();
        assert_eq!(
            classes_reading(CacheTable::StructuralVariants),
            tables.cache_reading_classes()
        );
        assert_eq!(
            classes_reading(CacheTable::IntronEnds),
            tables.intron_gap_classes()
        );
        for class in tables.transcript_reading_classes() {
            assert!(
                classes_reading(CacheTable::DownstreamSequences).contains(class),
                "{class}"
            );
        }
    }
}
