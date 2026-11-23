// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The database of exon scores, a SQLite file of two tables with one row of each for
//! every exon of every transcript record. LoFTEE reads `gerp_exons` with
//! `SELECT * FROM gerp_exons where transcript_id = ? AND exon_num = ?`
//! (`gerp_dist.pl:129`) and its column `gerp` (`gerp_dist.pl:133`), and `phylocsf_data`
//! with `SELECT * FROM phylocsf_data WHERE transcript = ? AND exon_number = ?`
//! (`LoF.pm:570`) and its columns `corresponding_orf_score` and `max_score`
//! (`LoF.pm:293-296`). Every value of the file is made up.
//!
//! Each table has an index on the two columns its query names: a query of a genome-wide
//! table without one reads the whole table. An index changes no row a query returns.

use std::path::Path;

use rusqlite::{params, Connection, OpenFlags};
use sha2::{Digest, Sha256};

use crate::exons::Record;
use crate::values::{exon_total, hundredths_real, phylocsf_scores, thousandths_real};

/// The two tables, as `scripts/fixtures/loftee_conservation_database.py` creates them.
/// The key columns carry a declared type: the plugin binds both values of a query as
/// text, and SQLite compares an INTEGER column with a bound text as numbers.
pub const SCHEMA: [&str; 2] = [
    "CREATE TABLE gerp_exons (transcript_id TEXT, exon_num INTEGER, gerp REAL)",
    "CREATE TABLE phylocsf_data (transcript TEXT, exon_number INTEGER, \
     corresponding_orf_score REAL, max_score REAL)",
];

/// One index a table, on the two columns the plugin's query names.
pub const INDEXES: [&str; 2] = [
    "CREATE INDEX gerp_exons_key ON gerp_exons (transcript_id, exon_num)",
    "CREATE INDEX phylocsf_data_key ON phylocsf_data (transcript, exon_number)",
];

/// The rows of each table in the order their digest takes them.
const ORDERED: [(&str, &str); 2] = [
    (
        "gerp_exons",
        "SELECT transcript_id, exon_num, gerp FROM gerp_exons \
         ORDER BY transcript_id, exon_num, gerp",
    ),
    (
        "phylocsf_data",
        "SELECT transcript, exon_number, corresponding_orf_score, max_score FROM phylocsf_data \
         ORDER BY transcript, exon_number, corresponding_orf_score, max_score",
    ),
];

/// The writer of the database: the two tables, their rows inside one transaction, then
/// the two indexes.
pub struct DatabaseWriter {
    connection: Connection,
    seed: u64,
}

fn failed(path: &Path, error: rusqlite::Error) -> String {
    format!("{}: {error}", path.display())
}

impl DatabaseWriter {
    /// Creates the database at `path`, a file that does not exist. The file is one a
    /// failed run throws away, so it is written without a rollback journal.
    pub fn create(path: &Path, seed: u64) -> Result<DatabaseWriter, String> {
        let connection = Connection::open(path).map_err(|error| failed(path, error))?;
        connection
            .execute_batch("PRAGMA journal_mode = OFF; PRAGMA synchronous = OFF; BEGIN")
            .and_then(|()| {
                SCHEMA
                    .iter()
                    .try_for_each(|table| connection.execute_batch(table))
            })
            .map_err(|error| failed(path, error))?;
        Ok(DatabaseWriter { connection, seed })
    }

    /// Adds the rows of `records`: for each record and each of its exons, from 1, one
    /// row of `gerp_exons` and one of `phylocsf_data`.
    pub fn add(&mut self, records: &[Record]) -> Result<(), rusqlite::Error> {
        let mut total = self
            .connection
            .prepare("INSERT INTO gerp_exons VALUES (?1, ?2, ?3)")?;
        let mut scores = self
            .connection
            .prepare("INSERT INTO phylocsf_data VALUES (?1, ?2, ?3, ?4)")?;
        for record in records {
            let id = record.stable_id.as_str();
            for exon in 1..=record.exons {
                let number = exon as i64;
                let gerp = hundredths_real(exon_total(self.seed, id, exon));
                total.execute(params![id, number, gerp])?;
                let phylocsf = phylocsf_scores(self.seed, id, exon);
                scores.execute(params![
                    id,
                    number,
                    thousandths_real(phylocsf.annotated),
                    thousandths_real(phylocsf.highest),
                ])?;
            }
        }
        Ok(())
    }

    /// Ends the transaction, builds the two indexes and closes the file.
    pub fn finish(self) -> Result<(), rusqlite::Error> {
        self.connection.execute_batch("COMMIT")?;
        for index in INDEXES {
            self.connection.execute_batch(index)?;
        }
        self.connection.close().map_err(|(_, error)| error)
    }
}

/// What a written database holds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatabaseRows {
    /// The rows of `gerp_exons`.
    pub gerp_exons: u64,
    /// The rows of `phylocsf_data`.
    pub phylocsf_data: u64,
    /// The sha256 of the rows as lines: the table's name and the row's values, joined
    /// by tabs, first the rows of `gerp_exons` in the order of its three columns, then
    /// those of `phylocsf_data` in the order of its four. A REAL stands as the shortest
    /// decimal that reads back as the same double.
    pub sha256: String,
}

/// Reads the database at `path` back, read-only, into its row counts and digest.
pub fn read_rows(path: &Path) -> Result<DatabaseRows, String> {
    let connection = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)
        .map_err(|error| failed(path, error))?;
    connection
        .execute_batch("PRAGMA temp_store = MEMORY")
        .map_err(|error| failed(path, error))?;
    let mut hasher = Sha256::new();
    let mut counts = [0u64; 2];
    for ((table, select), count) in ORDERED.iter().zip(&mut counts) {
        let mut statement = connection
            .prepare(select)
            .map_err(|error| failed(path, error))?;
        let scores = statement.column_count() - 2;
        let mut rows = statement.query([]).map_err(|error| failed(path, error))?;
        while let Some(row) = rows.next().map_err(|error| failed(path, error))? {
            let id: String = row.get(0).map_err(|error| failed(path, error))?;
            let exon: i64 = row.get(1).map_err(|error| failed(path, error))?;
            let mut line = format!("{table}\t{id}\t{exon}");
            for column in 2..2 + scores {
                let score: f64 = row.get(column).map_err(|error| failed(path, error))?;
                line.push_str(&format!("\t{score}"));
            }
            line.push('\n');
            hasher.update(line.as_bytes());
            *count += 1;
        }
    }
    let mut sha256 = String::with_capacity(64);
    for byte in hasher.finalize() {
        sha256.push_str(&format!("{byte:02x}"));
    }
    Ok(DatabaseRows {
        gerp_exons: counts[0],
        phylocsf_data: counts[1],
        sha256,
    })
}
