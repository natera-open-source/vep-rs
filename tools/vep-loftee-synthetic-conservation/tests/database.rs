// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The database of exon scores: its two tables and their indexes, which exons have a
//! row, what a row holds, and what the plugin's own queries read from it.

mod common;

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use common::{
    database, new_cache, phylocsf, printed, printed_field, small_cache, table_rows, thousandths,
    totals, transcript, with_db_id, write_files, write_shard, Scratch, SMALL,
};
use rusqlite::{Connection, OpenFlags};
use vep_loftee_synthetic_conservation::values::{
    exon_total, hundredths_real, phylocsf_scores, thousandths_real,
};
use vep_loftee_synthetic_conservation::DATABASE_FILE;

/// The query of `gerp_dist.pl:129`.
const TOTAL_QUERY: &str = "SELECT * FROM gerp_exons where transcript_id = ? AND exon_num = ?";

/// The query of `LoF.pm:570`.
const PHYLOCSF_QUERY: &str = "SELECT * FROM phylocsf_data WHERE transcript = ? AND exon_number = ?";

#[test]
fn every_exon_has_a_total_and_its_phylocsf_rows() {
    let scratch = Scratch::new("every-exon");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    let stdout = write_files(&cache, &out, 1, &[]);
    // One row of each table for each exon, a transcript's exons numbered from 1, in the
    // order of the sequences, of the transcripts by start and stable id, and of the
    // exons.
    let order = ["MADE_T1", "MADE_T2", "MADE_N1", "MADE_T3", "MADE_N2"];
    let mut keys = Vec::new();
    for id in order {
        let (.., exons) = SMALL.iter().find(|(_, of, ..)| *of == id).unwrap();
        keys.extend((1..=exons.len() as i64).map(|exon| (id.to_owned(), exon)));
    }
    assert_eq!(keys.len(), 9);

    let totals = totals(&out);
    let total_keys: Vec<(String, i64)> = totals
        .iter()
        .map(|(id, exon, _)| (id.clone(), *exon))
        .collect();
    assert_eq!(total_keys, keys);
    for (id, exon, total) in &totals {
        assert_eq!(*total, hundredths_real(exon_total(1, id, *exon as u64)));
    }

    let scores = phylocsf(&out);
    let score_keys: Vec<(String, i64)> = scores
        .iter()
        .map(|(id, exon, ..)| (id.clone(), *exon))
        .collect();
    assert_eq!(score_keys, keys);
    for (id, exon, annotated, highest) in &scores {
        let drawn = phylocsf_scores(1, id, *exon as u64);
        assert_eq!(*annotated, thousandths_real(drawn.annotated));
        assert_eq!(*highest, thousandths_real(drawn.highest));
        assert!(highest >= annotated, "{id} exon {exon}");
    }

    assert_eq!(printed(&stdout, "exons"), ["9"]);
    assert_eq!(printed_field(&stdout, DATABASE_FILE, "rows"), "18");
    assert_eq!(printed_field(&stdout, DATABASE_FILE, "gerp_exons"), "9");
    assert_eq!(printed_field(&stdout, DATABASE_FILE, "phylocsf_data"), "9");
}

#[test]
fn an_exon_total_is_not_the_sum_of_its_base_scores() {
    let scratch = Scratch::new("own-total");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    write_files(&cache, &out, 1, &[]);
    let scores: BTreeMap<(String, u64), i64> = table_rows(&out)
        .into_iter()
        .map(|(sequence, position, score)| ((sequence, position), thousandths(&score)))
        .collect();
    let totals: BTreeMap<(String, i64), f64> = totals(&out)
        .into_iter()
        .map(|(id, exon, total)| ((id, exon), total))
        .collect();
    let mut exons = 0;
    for (sequence, id, _, _, spans) in SMALL {
        for (number, (start, end)) in spans.iter().enumerate() {
            let summed: i64 = (*start..=*end)
                .map(|position| scores[&(sequence.to_owned(), position)])
                .sum();
            let stored = totals[&(id.to_owned(), number as i64 + 1)];
            assert_ne!(
                thousandths_real(summed),
                stored,
                "{id} exon {}: the total is the sum of the exon's rows",
                number + 1
            );
            exons += 1;
        }
    }
    assert_eq!(exons, 9);
}

/// The columns of a table: each one's place, name and declared type.
type Columns = Vec<(i64, String, String)>;

/// The tables of a database, each with its columns.
fn tables(connection: &Connection) -> Vec<(String, Columns)> {
    let mut statement = connection
        .prepare("SELECT name FROM sqlite_master WHERE type = 'table' ORDER BY name")
        .unwrap();
    let names: Vec<String> = statement
        .query_map([], |row| row.get(0))
        .unwrap()
        .map(|name| name.unwrap())
        .collect();
    names
        .into_iter()
        .map(|name| {
            let mut statement = connection
                .prepare(&format!("PRAGMA table_info({name})"))
                .unwrap();
            let columns = statement
                .query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))
                .unwrap()
                .map(|column| column.unwrap())
                .collect();
            (name, columns)
        })
        .collect()
}

#[test]
fn the_tables_are_the_ones_the_corpus_generator_writes() {
    let committed = Path::new(env!("CARGO_MANIFEST_DIR")).join(
        "../../tests/golden/116/GRCh37-plugin-loftee-fixture/plugin_data/phylocsf_gerp.synthetic.sql",
    );
    let committed = Connection::open_with_flags(committed, OpenFlags::SQLITE_OPEN_READ_ONLY)
        .expect("the committed database");
    let scratch = Scratch::new("tables");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    write_files(&cache, &out, 1, &[]);
    let written = tables(&database(&out));
    assert_eq!(written, tables(&committed));
    let column = |cid: i64, name: &str, declared: &str| (cid, name.to_owned(), declared.to_owned());
    assert_eq!(
        written,
        [
            (
                "gerp_exons".to_owned(),
                vec![
                    column(0, "transcript_id", "TEXT"),
                    column(1, "exon_num", "INTEGER"),
                    column(2, "gerp", "REAL"),
                ]
            ),
            (
                "phylocsf_data".to_owned(),
                vec![
                    column(0, "transcript", "TEXT"),
                    column(1, "exon_number", "INTEGER"),
                    column(2, "corresponding_orf_score", "REAL"),
                    column(3, "max_score", "REAL"),
                ]
            ),
        ]
    );
}

#[test]
fn the_plugin_s_queries_read_one_row_through_an_index() {
    let scratch = Scratch::new("queries");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    write_files(&cache, &out, 1, &[]);
    let connection = database(&out);
    // The plugin binds both values of a query as text.
    let total: f64 = connection
        .query_row(TOTAL_QUERY, ["MADE_T2", "3"], |row| row.get("gerp"))
        .unwrap();
    assert_eq!(total, hundredths_real(exon_total(1, "MADE_T2", 3)));
    let (annotated, highest): (f64, f64) = connection
        .query_row(PHYLOCSF_QUERY, ["MADE_T2", "3"], |row| {
            Ok((row.get("corresponding_orf_score")?, row.get("max_score")?))
        })
        .unwrap();
    let drawn = phylocsf_scores(1, "MADE_T2", 3);
    assert_eq!(annotated, thousandths_real(drawn.annotated));
    assert_eq!(highest, thousandths_real(drawn.highest));
    // An exon number the transcript does not have reads no row.
    let none = connection.query_row(TOTAL_QUERY, ["MADE_T2", "4"], |row| {
        row.get::<_, f64>("gerp")
    });
    assert!(matches!(none, Err(rusqlite::Error::QueryReturnedNoRows)));

    for (query, index) in [
        (TOTAL_QUERY, "gerp_exons_key"),
        (PHYLOCSF_QUERY, "phylocsf_data_key"),
    ] {
        let plan: Vec<String> = connection
            .prepare(&format!("EXPLAIN QUERY PLAN {query}"))
            .unwrap()
            .query_map(["MADE_T2", "3"], |row| row.get("detail"))
            .unwrap()
            .map(|step| step.unwrap())
            .collect();
        assert!(
            plan.iter()
                .any(|step| step.starts_with("SEARCH") && step.contains(index)),
            "{query}: {plan:?}"
        );
    }
    let indexes: Vec<String> = connection
        .prepare("SELECT name FROM sqlite_master WHERE type = 'index' ORDER BY name")
        .unwrap()
        .query_map([], |row| row.get(0))
        .unwrap()
        .map(|name| name.unwrap())
        .collect();
    assert_eq!(indexes, ["gerp_exons_key", "phylocsf_data_key"]);
}

#[test]
fn a_transcript_filed_in_two_region_files_has_one_row_an_exon() {
    let scratch = Scratch::new("two-region-files");
    let cache = new_cache(scratch.path(), &["long"]);
    // A gene that reaches both regions: the cache files both of its transcripts in
    // both region files.
    let low = transcript(
        "MADE_LOW",
        1,
        "protein_coding",
        &[(500_001, 500_050), (500_101, 500_150)],
    );
    let across = transcript(
        "MADE_ACROSS",
        1,
        "protein_coding",
        &[
            (999_901, 1_000_000),
            (1_000_101, 1_000_200),
            (1_000_301, 1_000_310),
        ],
    );
    let gene = [low, across];
    write_shard(&cache, "long", "1-1000000.json", &gene);
    write_shard(&cache, "long", "1000001-2000000.json", &gene);
    let out = scratch.path().join("out");
    let stdout = write_files(&cache, &out, 1, &[]);
    assert_eq!(printed(&stdout, "transcripts"), ["2"]);
    assert_eq!(printed(&stdout, "exons"), ["5"]);
    assert_eq!(printed(&stdout, "exon_bases"), ["310"]);
    let keys: Vec<(String, i64)> = totals(&out)
        .into_iter()
        .map(|(id, exon, _)| (id, exon))
        .collect();
    let expected: Vec<(String, i64)> = [
        ("MADE_LOW", 1),
        ("MADE_LOW", 2),
        ("MADE_ACROSS", 1),
        ("MADE_ACROSS", 2),
        ("MADE_ACROSS", 3),
    ]
    .iter()
    .map(|(id, exon)| ((*id).to_owned(), *exon))
    .collect();
    assert_eq!(keys, expected);
    let scored: Vec<(String, i64)> = phylocsf(&out)
        .into_iter()
        .map(|(id, exon, ..)| (id, exon))
        .collect();
    assert_eq!(scored, expected);
    assert_eq!(table_rows(&out).len(), 310);
}

#[test]
fn one_stable_id_under_two_dbids_keeps_the_exons_of_both_records() {
    let scratch = Scratch::new("two-dbids");
    let cache = new_cache(scratch.path(), &["long"]);
    // One stable id at two places under two database ids, and at the first place under
    // a third: three records, as the cache reader tells them apart.
    let first = with_db_id(
        transcript("MADE_X", 1, "protein_coding", &[(101, 200), (301, 400)]),
        11,
    );
    let elsewhere = with_db_id(
        transcript(
            "MADE_X",
            1,
            "protein_coding",
            &[(1_001, 1_100), (1_201, 1_300), (1_401, 1_500)],
        ),
        12,
    );
    let aligned = with_db_id(
        transcript("MADE_X", 1, "cDNA", &[(101, 200), (301, 400)]),
        13,
    );
    write_shard(
        &cache,
        "long",
        "1-1000000.json",
        &[first, elsewhere, aligned],
    );
    let out = scratch.path().join("out");
    let stdout = write_files(&cache, &out, 1, &[]);
    assert_eq!(printed(&stdout, "transcripts"), ["3"]);
    assert_eq!(printed(&stdout, "exons"), ["7"]);
    assert_eq!(printed(&stdout, "exon_bases"), ["500"]);
    assert_eq!(printed_field(&stdout, DATABASE_FILE, "gerp_exons"), "7");
    assert_eq!(printed_field(&stdout, DATABASE_FILE, "phylocsf_data"), "7");

    // The exons of both places have rows in the table.
    let positions: BTreeSet<u64> = table_rows(&out)
        .into_iter()
        .map(|(_, position, _)| position)
        .collect();
    let expected: BTreeSet<u64> = (101..=200)
        .chain(301..=400)
        .chain(1_001..=1_100)
        .chain(1_201..=1_300)
        .chain(1_401..=1_500)
        .collect();
    assert_eq!(positions, expected);

    // Each record has a row for each of its exons: three rows under exon 1 and under
    // exon 2, one under exon 3, and the rows of one key hold one value.
    let mut by_exon: BTreeMap<i64, Vec<f64>> = BTreeMap::new();
    for (id, exon, total) in totals(&out) {
        assert_eq!(id, "MADE_X");
        by_exon.entry(exon).or_default().push(total);
    }
    let counts: Vec<(i64, usize)> = by_exon
        .iter()
        .map(|(exon, totals)| (*exon, totals.len()))
        .collect();
    assert_eq!(counts, [(1, 3), (2, 3), (3, 1)]);
    for (exon, totals) in &by_exon {
        assert!(
            totals.iter().all(|total| *total == totals[0]),
            "exon {exon}"
        );
        assert_eq!(
            totals[0],
            hundredths_real(exon_total(1, "MADE_X", *exon as u64))
        );
    }
    let mut scored: BTreeMap<i64, Vec<(f64, f64)>> = BTreeMap::new();
    for (_, exon, annotated, highest) in phylocsf(&out) {
        scored.entry(exon).or_default().push((annotated, highest));
    }
    let counts: Vec<(i64, usize)> = scored
        .iter()
        .map(|(exon, scores)| (*exon, scores.len()))
        .collect();
    assert_eq!(counts, [(1, 3), (2, 3), (3, 1)]);
    for scores in scored.values() {
        assert!(scores.iter().all(|pair| *pair == scores[0]));
    }
}

#[test]
fn the_values_stand_on_both_sides_of_what_the_plugin_tests() {
    let scratch = Scratch::new("both-sides");
    let cache = new_cache(scratch.path(), &["many"]);
    // One transcript of 400 exons of two bases.
    let exons: Vec<(u64, u64)> = (0..400)
        .map(|exon| (exon * 10 + 1, exon * 10 + 2))
        .collect();
    write_shard(
        &cache,
        "many",
        "1-1000000.json",
        &[transcript("MADE_MANY", 1, "protein_coding", &exons)],
    );
    let out = scratch.path().join("out");
    write_files(&cache, &out, 1, &[]);
    // `LoF.pm:268` holds the walk's total against 180.
    let totals = totals(&out);
    assert_eq!(totals.len(), 400);
    assert!(totals.iter().any(|(.., total)| *total > 180.0));
    assert!(totals.iter().any(|(.., total)| *total <= 180.0));
    assert!(totals.iter().any(|(.., total)| *total < 0.0));
    assert!(totals
        .iter()
        .all(|(.., total)| (-300.0..=900.0).contains(total)));
    // `LoF.pm:295-296` flags a negative annotated score by the sign of the highest.
    let scores = phylocsf(&out);
    assert_eq!(scores.len(), 400);
    let on = |keep: &dyn Fn(f64, f64) -> bool| {
        scores
            .iter()
            .filter(|(.., annotated, highest)| keep(*annotated, *highest))
            .count()
    };
    assert!(on(&|annotated, _| annotated < 0.0) > 0 && on(&|annotated, _| annotated > 0.0) > 0);
    assert!(on(&|_, highest| highest < 0.0) > 0 && on(&|_, highest| highest > 0.0) > 0);
    assert!(on(&|annotated, _| annotated == 0.0) > 0 && on(&|_, highest| highest == 0.0) > 0);
    let unflagged = on(&|annotated, _| annotated >= 0.0);
    let unlikely = on(&|annotated, highest| annotated < 0.0 && highest > 0.0);
    let weak = on(&|annotated, highest| annotated < 0.0 && highest <= 0.0);
    assert_eq!(unflagged + unlikely + weak, 400);
    for (flag, count) in [("none", unflagged), ("unlikely", unlikely), ("weak", weak)] {
        assert!((100..=170).contains(&count), "{flag}: {count} of 400");
    }
    // A score has at most three decimals.
    for (.., annotated, highest) in &scores {
        for score in [annotated, highest] {
            assert_eq!((score * 1000.0).round() / 1000.0, *score);
        }
    }
}
