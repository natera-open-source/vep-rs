// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests share: made-up caches written into scratch directories, a runner of
//! the built binary, and readers of the three files a run writes.

#![allow(dead_code)]

use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use noodles::bgzf;
use rusqlite::{Connection, OpenFlags};
use serde_json::{json, Value};
use vep_loftee_synthetic_conservation::{DATABASE_FILE, INDEX_FILE, TABLE_FILE};

/// A directory of this test's own, removed when the value is dropped.
pub struct Scratch(PathBuf);

impl Scratch {
    pub fn new(label: &str) -> Scratch {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!(
            "vep-loftee-synthetic-conservation-{label}-{}-{}",
            std::process::id(),
            COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).expect("a scratch directory");
        Scratch(dir)
    }

    pub fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A transcript record of the smallest shape the cache reader takes: a stable id, a
/// strand, a biotype and its exons, with the place the exons span.
pub fn transcript(id: &str, strand: i8, biotype: &str, exons: &[(u64, u64)]) -> Value {
    let start = exons.iter().map(|exon| exon.0).min().expect("an exon");
    let end = exons.iter().map(|exon| exon.1).max().expect("an exon");
    json!({
        "stable_id": id,
        "start": start,
        "end": end,
        "strand": strand,
        "biotype": biotype,
        "exons": exons
            .iter()
            .map(|(start, end)| json!({"start": start, "end": end}))
            .collect::<Vec<_>>(),
    })
}

/// `record` under the database id `id`.
pub fn with_db_id(mut record: Value, id: u64) -> Value {
    record["dbID"] = json!(id);
    record
}

/// Starts a made-up cache under `dir` whose `info.json` lists `sequences` in order.
pub fn new_cache(dir: &Path, sequences: &[&str]) -> PathBuf {
    let cache = dir.join("json_cache");
    fs::create_dir_all(cache.join("transcripts")).expect("a cache directory");
    let info = json!({"species": "made_up", "assembly": "MADE1", "contig_order": sequences});
    fs::write(cache.join("info.json"), info.to_string()).expect("a written info.json");
    cache
}

/// Writes `records` as the region file `file` of `sequence`.
pub fn write_shard(cache: &Path, sequence: &str, file: &str, records: &[Value]) {
    let dir = cache.join("transcripts").join(sequence);
    fs::create_dir_all(&dir).expect("a sequence directory");
    fs::write(
        dir.join(file),
        serde_json::to_vec(records).expect("a shard's bytes"),
    )
    .expect("a written shard");
}

/// A made-up transcript: its sequence, stable id, strand and biotype, and its exons
/// from its 5' end to its 3' end.
pub type MadeUp = (
    &'static str,
    &'static str,
    i8,
    &'static str,
    &'static [(u64, u64)],
);

/// The transcripts of the small made-up cache. `first` holds two coding transcripts
/// whose exons share bases and one lncRNA; `second` holds one coding transcript and one
/// snRNA.
pub const SMALL: [MadeUp; 5] = [
    (
        "first",
        "MADE_T1",
        1,
        "protein_coding",
        &[(101, 200), (301, 450)],
    ),
    (
        "first",
        "MADE_T2",
        -1,
        "protein_coding",
        &[(701, 800), (401, 500), (151, 250)],
    ),
    ("first", "MADE_N1", 1, "lncRNA", &[(1001, 1100)]),
    (
        "second",
        "MADE_T3",
        1,
        "protein_coding",
        &[(11, 30), (51, 60)],
    ),
    ("second", "MADE_N2", -1, "snRNA", &[(100, 110)]),
];

/// The sequences of the small made-up cache, in the order `info.json` lists them.
pub const SMALL_SEQUENCES: [&str; 2] = ["first", "second"];

/// Writes the small made-up cache under `dir`.
pub fn small_cache(dir: &Path) -> PathBuf {
    let cache = new_cache(dir, &SMALL_SEQUENCES);
    for sequence in SMALL_SEQUENCES {
        let records: Vec<Value> = SMALL
            .iter()
            .filter(|(on, ..)| *on == sequence)
            .map(|(_, id, strand, biotype, exons)| transcript(id, *strand, biotype, exons))
            .collect();
        write_shard(&cache, sequence, "1-1000000.json", &records);
    }
    cache
}

/// How a run of the built binary ended.
pub struct Ended {
    pub code: Option<i32>,
    pub stdout: String,
    pub stderr: String,
}

/// Runs the built binary.
pub fn run(args: &[&str]) -> Ended {
    let output = Command::new(env!("CARGO_BIN_EXE_vep-loftee-synthetic-conservation"))
        .args(args)
        .output()
        .expect("the built binary runs");
    Ended {
        code: output.status.code(),
        stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
    }
}

/// Runs the binary on `cache` with `seed` into `out` and gives what it printed.
pub fn write_files(cache: &Path, out: &Path, seed: u64, more: &[&str]) -> String {
    let seed = seed.to_string();
    let mut args = vec![
        "--cache",
        cache.to_str().expect("a path"),
        "--seed",
        &seed,
        "--out-dir",
        out.to_str().expect("a path"),
    ];
    args.extend_from_slice(more);
    let ended = run(&args);
    assert_eq!(ended.code, Some(0), "stderr: {}", ended.stderr);
    assert_eq!(ended.stderr, "");
    ended.stdout
}

/// The fields after `key` of the line a run printed for it.
pub fn printed(stdout: &str, key: &str) -> Vec<String> {
    let line = stdout
        .lines()
        .find(|line| line.split('\t').next() == Some(key))
        .unwrap_or_else(|| panic!("no line for {key} in {stdout}"));
    line.split('\t').skip(1).map(str::to_owned).collect()
}

/// The value of `name=` among the fields of the line a run printed for `key`.
pub fn printed_field(stdout: &str, key: &str, name: &str) -> String {
    printed(stdout, key)
        .iter()
        .find_map(|field| field.strip_prefix(&format!("{name}=")).map(str::to_owned))
        .unwrap_or_else(|| panic!("no {name} on the line of {key}"))
}

pub fn table_path(out: &Path) -> PathBuf {
    out.join(TABLE_FILE)
}

pub fn index_path(out: &Path) -> PathBuf {
    out.join(INDEX_FILE)
}

pub fn database_path(out: &Path) -> PathBuf {
    out.join(DATABASE_FILE)
}

/// The table's rows as text, decompressed.
pub fn table_text(out: &Path) -> String {
    let mut reader = bgzf::io::Reader::new(File::open(table_path(out)).expect("the table"));
    let mut text = String::new();
    reader.read_to_string(&mut text).expect("the table's rows");
    text
}

/// One row of the table: the sequence, the position and the score as written.
pub type Row = (String, u64, String);

/// The row a line of the table holds.
pub fn row_of(line: &str) -> Row {
    let fields: Vec<&str> = line.split('\t').collect();
    let [sequence, position, score] = fields.as_slice() else {
        panic!("{line:?} does not have three fields");
    };
    (
        (*sequence).to_owned(),
        position.parse().expect("a position"),
        (*score).to_owned(),
    )
}

/// The table's rows in file order.
pub fn table_rows(out: &Path) -> Vec<Row> {
    table_text(out).lines().map(row_of).collect()
}

/// A score as written, in thousandths.
pub fn thousandths(score: &str) -> i64 {
    let (sign, size) = match score.strip_prefix('-') {
        Some(size) => (-1, size),
        None => (1, score),
    };
    let (whole, fraction) = size.split_once('.').unwrap_or((size, ""));
    assert!(fraction.len() <= 3, "{score} has more than three decimals");
    let whole: i64 = whole.parse().expect("a whole part");
    let fraction: i64 = format!("{fraction:0<3}").parse().expect("a fraction");
    sign * (whole * 1000 + fraction)
}

/// The database a run wrote, opened read-only.
pub fn database(out: &Path) -> Connection {
    Connection::open_with_flags(database_path(out), OpenFlags::SQLITE_OPEN_READ_ONLY)
        .expect("the database")
}

/// The rows of `gerp_exons` in the order they were inserted: the stable id, the exon
/// number and the total.
pub fn totals(out: &Path) -> Vec<(String, i64, f64)> {
    let connection = database(out);
    let mut statement = connection
        .prepare("SELECT transcript_id, exon_num, gerp FROM gerp_exons ORDER BY rowid")
        .expect("a query");
    let rows = statement
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))
        .expect("the rows");
    rows.map(|row| row.expect("a row")).collect()
}

/// The rows of `phylocsf_data` in the order they were inserted: the stable id, the exon
/// number, the annotated score and the highest score.
pub fn phylocsf(out: &Path) -> Vec<(String, i64, f64, f64)> {
    let connection = database(out);
    let mut statement = connection
        .prepare(
            "SELECT transcript, exon_number, corresponding_orf_score, max_score \
             FROM phylocsf_data ORDER BY rowid",
        )
        .expect("a query");
    let rows = statement
        .query_map([], |row| {
            Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?))
        })
        .expect("the rows");
    rows.map(|row| row.expect("a row")).collect()
}
