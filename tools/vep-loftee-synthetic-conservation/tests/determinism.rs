// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A seed fixes the rows: the same cache and seed write the same rows on any number of
//! threads, and another seed writes other values under the same keys.

mod common;

use std::fs;

use common::{
    index_path, phylocsf, printed_field, small_cache, table_path, table_rows, table_text, totals,
    write_files, Scratch,
};
use vep_loftee_synthetic_conservation::{DATABASE_FILE, TABLE_FILE};

#[test]
fn the_same_seed_writes_the_same_rows() {
    let scratch = Scratch::new("same-seed");
    let cache = small_cache(scratch.path());
    let (one, other) = (scratch.path().join("one"), scratch.path().join("other"));
    let printed_one = write_files(&cache, &one, 7, &["--threads", "1"]);
    let printed_other = write_files(&cache, &other, 7, &["--threads", "4"]);
    assert_eq!(printed_one, printed_other);
    assert_eq!(table_text(&one), table_text(&other));
    assert_eq!(totals(&one), totals(&other));
    assert_eq!(phylocsf(&one), phylocsf(&other));
    // The table and its index are the same bytes on any number of threads.
    assert_eq!(
        fs::read(table_path(&one)).unwrap(),
        fs::read(table_path(&other)).unwrap()
    );
    assert_eq!(
        fs::read(index_path(&one)).unwrap(),
        fs::read(index_path(&other)).unwrap()
    );
}

#[test]
fn another_seed_writes_other_rows() {
    let scratch = Scratch::new("other-seed");
    let cache = small_cache(scratch.path());
    let (one, other) = (scratch.path().join("one"), scratch.path().join("other"));
    let printed_one = write_files(&cache, &one, 1, &[]);
    let printed_other = write_files(&cache, &other, 2, &[]);
    for file in [TABLE_FILE, DATABASE_FILE] {
        assert_eq!(
            printed_field(&printed_one, file, "rows"),
            printed_field(&printed_other, file, "rows")
        );
        let digest = printed_field(&printed_one, file, "rows_sha256");
        assert_eq!(digest.len(), 64);
        assert_ne!(digest, printed_field(&printed_other, file, "rows_sha256"));
    }

    // The same keys in the same order, and other values under most of them.
    let (rows_one, rows_other) = (table_rows(&one), table_rows(&other));
    assert_eq!(rows_one.len(), rows_other.len());
    let mut moved = 0;
    for (row, other) in rows_one.iter().zip(&rows_other) {
        assert_eq!((&row.0, row.1), (&other.0, other.1));
        moved += usize::from(row.2 != other.2);
    }
    assert!(moved > rows_one.len() * 9 / 10, "{moved} scores differ");

    let (totals_one, totals_other) = (totals(&one), totals(&other));
    assert_eq!(totals_one.len(), totals_other.len());
    for (total, other) in totals_one.iter().zip(&totals_other) {
        assert_eq!((&total.0, total.1), (&other.0, other.1));
        assert_ne!(total.2, other.2, "{} exon {}", total.0, total.1);
    }
    let (scores_one, scores_other) = (phylocsf(&one), phylocsf(&other));
    assert_eq!(scores_one.len(), scores_other.len());
    for (scores, other) in scores_one.iter().zip(&scores_other) {
        assert_eq!((&scores.0, scores.1), (&other.0, other.1));
        assert_ne!(
            (scores.2, scores.3),
            (other.2, other.3),
            "{} exon {}",
            scores.0,
            scores.1
        );
    }
}
