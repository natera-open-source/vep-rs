// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The table of base scores and its index: which bases have a row, what a row holds
//! and what a region read through the index gives, by this crate's own reading of the
//! index and by the `tabix` program where the machine has one.

mod common;

use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};
use std::process::Command;

use common::{
    index_path, new_cache, printed, printed_field, row_of, small_cache, table_path, table_rows,
    table_text, totals, transcript, write_files, write_shard, Row, Scratch, SMALL, SMALL_SEQUENCES,
};
use noodles::bgzf;
use noodles::core::region::Interval;
use noodles::core::Position;
use noodles::csi::{self, BinningIndex};
use noodles::tabix;
use serde_json::Value;
use vep_loftee_synthetic_conservation::values::{base_score, thousandths_text};
use vep_loftee_synthetic_conservation::{DATABASE_FILE, INDEX_FILE, TABLE_FILE};

/// The bases of the small made-up cache that lie in an exon, each once: the sequences
/// in the cache's order, a sequence's bases in ascending position.
fn small_bases() -> Vec<(String, u64)> {
    let mut bases = Vec::new();
    for sequence in SMALL_SEQUENCES {
        let on: BTreeSet<u64> = SMALL
            .iter()
            .filter(|(of, ..)| *of == sequence)
            .flat_map(|(.., exons)| exons.iter().flat_map(|(start, end)| *start..=*end))
            .collect();
        bases.extend(
            on.into_iter()
                .map(|position| (sequence.to_owned(), position)),
        );
    }
    bases
}

fn places(rows: &[Row]) -> Vec<(String, u64)> {
    rows.iter()
        .map(|(sequence, position, _)| (sequence.clone(), *position))
        .collect()
}

#[test]
fn every_exon_base_has_one_row() {
    let scratch = Scratch::new("every-base");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    let stdout = write_files(&cache, &out, 1, &[]);
    let rows = table_rows(&out);
    // 150, 200, 100 and 100 bases on the first sequence, 20, 10 and 11 on the second.
    assert_eq!(rows.len(), 591);
    assert_eq!(places(&rows), small_bases());
    for (sequence, position, score) in &rows {
        assert_eq!(
            *score,
            thousandths_text(base_score(1, sequence, *position)),
            "{sequence}:{position}"
        );
    }
    assert_eq!(printed(&stdout, "exon_bases"), ["591"]);
    assert_eq!(printed_field(&stdout, TABLE_FILE, "rows"), "591");
    let text = table_text(&out);
    assert!(text.starts_with("first\t101\t") && text.ends_with('\n'));
}

#[test]
fn a_base_in_two_exons_is_written_once() {
    let scratch = Scratch::new("shared-base");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    write_files(&cache, &out, 1, &[]);
    let rows = table_rows(&out);
    // Bases 151 to 200 and 401 to 450 lie in an exon of MADE_T1 and in one of MADE_T2.
    for position in (151..=200).chain(401..=450) {
        let held = rows
            .iter()
            .filter(|(sequence, at, _)| sequence == "first" && *at == position)
            .count();
        assert_eq!(held, 1, "first:{position}");
    }

    // A further transcript whose exons lie inside those of MADE_T1 adds exons and no row.
    let records: Vec<Value> = SMALL
        .iter()
        .filter(|(of, ..)| *of == "first")
        .map(|(_, id, strand, biotype, exons)| transcript(id, *strand, biotype, exons))
        .chain([transcript(
            "MADE_T4",
            1,
            "protein_coding",
            &[(101, 200), (120, 130), (301, 450)],
        )])
        .collect();
    write_shard(&cache, "first", "1-1000000.json", &records);
    let again = scratch.path().join("again");
    let stdout = write_files(&cache, &again, 1, &[]);
    assert_eq!(table_rows(&again), rows);
    assert_eq!(printed(&stdout, "exons"), ["12"]);
    assert_eq!(printed(&stdout, "exon_bases"), ["591"]);
}

#[test]
fn non_coding_transcripts_get_rows() {
    let scratch = Scratch::new("non-coding");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    write_files(&cache, &out, 1, &[]);
    let rows = places(&table_rows(&out));
    // The lncRNA's one exon on the first sequence and the snRNA's on the second.
    for position in 1001..=1100 {
        assert!(rows.contains(&("first".to_owned(), position)));
    }
    for position in 100..=110 {
        assert!(rows.contains(&("second".to_owned(), position)));
    }
    let keys: Vec<(String, i64)> = totals(&out)
        .into_iter()
        .map(|(id, exon, _)| (id, exon))
        .collect();
    assert!(keys.contains(&("MADE_N1".to_owned(), 1)));
    assert!(keys.contains(&("MADE_N2".to_owned(), 1)));
}

#[test]
fn every_written_file_name_holds_synthetic() {
    let scratch = Scratch::new("names");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    write_files(&cache, &out, 1, &[]);
    let mut names: Vec<String> = fs::read_dir(&out)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    assert_eq!(names, [TABLE_FILE, INDEX_FILE, DATABASE_FILE]);
    assert!(names.iter().all(|name| name.contains("synthetic")));
}

/// The `tabix` program: the one under the home directory's `.local/bin`, or the first
/// on the search path.
fn tabix_program() -> Option<PathBuf> {
    let home = std::env::var_os("HOME").map(|home| PathBuf::from(home).join(".local/bin"));
    let searched = std::env::var_os("PATH")
        .map(|path| std::env::split_paths(&path).collect::<Vec<_>>())
        .unwrap_or_default();
    home.into_iter()
        .chain(searched)
        .map(|dir| dir.join("tabix"))
        .find(|program| program.is_file())
}

/// What the `tabix` program prints for `args`.
fn tabix_lines(program: &Path, args: &[&str]) -> Vec<String> {
    let output = Command::new(program).args(args).output().unwrap();
    assert!(
        output.status.success(),
        "tabix {args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout)
        .lines()
        .map(str::to_owned)
        .collect()
}

/// The sequences the index names, in its order.
fn indexed_sequences(index: &tabix::Index) -> Vec<String> {
    index
        .header()
        .expect("a header")
        .reference_sequence_names()
        .iter()
        .map(ToString::to_string)
        .collect()
}

/// The rows of `sequence` from `start` to `end` that the index at `out` leads to: the
/// rows of the chunks it names for the region, those outside the region left out as a
/// reader of a tabix index leaves them out.
fn region_rows(out: &Path, sequence: &str, start: u64, end: u64) -> Vec<Row> {
    let index = tabix::fs::read(index_path(out)).expect("the index");
    let id = indexed_sequences(&index)
        .iter()
        .position(|name| name == sequence)
        .expect("a sequence of the index");
    let at = |position: u64| Position::try_from(position as usize).expect("a position");
    let interval: Interval = (at(start)..=at(end)).into();
    let chunks = index.query(id, interval).expect("the region's chunks");
    let mut reader = bgzf::io::Reader::new(File::open(table_path(out)).expect("the table"));
    csi::io::Query::new(&mut reader, chunks)
        .lines()
        .map(|line| row_of(&line.expect("a row")))
        .filter(|(of, position, _)| of == sequence && (start..=end).contains(position))
        .collect()
}

#[test]
fn the_index_gives_a_region_exactly_its_rows() {
    let scratch = Scratch::new("index");
    let cache = new_cache(scratch.path(), &["bare", "big", "small"]);
    write_shard(&cache, "bare", "1-1000000.json", &[]);
    // The first transcript's rows fill many blocks of the table; the second ends at the
    // last position a region read through a tabix index names.
    let last = 536_870_911;
    write_shard(
        &cache,
        "big",
        "1-1000000.json",
        &[
            transcript(
                "MADE_B1",
                1,
                "protein_coding",
                &[(5_000, 45_000), (60_000, 60_010), (100_000, 130_000)],
            ),
            transcript("MADE_B2", -1, "lncRNA", &[(last - 12, last)]),
        ],
    );
    write_shard(
        &cache,
        "small",
        "1-1000000.json",
        &[transcript("MADE_S1", 1, "protein_coding", &[(1, 3)])],
    );
    let out = scratch.path().join("out");
    write_files(&cache, &out, 1, &["--threads", "3"]);
    let rows = table_rows(&out);
    assert_eq!(rows.len(), 40_001 + 11 + 30_001 + 13 + 3);
    let table = table_path(&out);
    assert!(fs::metadata(&table).unwrap().len() > 4 * 65_536);

    let index = tabix::fs::read(index_path(&out)).expect("the index");
    assert_eq!(indexed_sequences(&index), ["big", "small"]);
    // The index's first fields, as `tabix -s 1 -b 2 -e 2` states them: the magic
    // string, two sequences, the generic format, the sequence in column 1, the begin
    // and the end in column 2, `#` for a comment line, no skipped line, and the names.
    let mut stated = Vec::new();
    bgzf::io::Reader::new(File::open(index_path(&out)).unwrap())
        .read_to_end(&mut stated)
        .unwrap();
    assert_eq!(&stated[..4], b"TBI\x01");
    let fields: Vec<i32> = stated[4..36]
        .chunks(4)
        .map(|field| i32::from_le_bytes(field.try_into().unwrap()))
        .collect();
    assert_eq!(fields, [2, 0, 1, 2, 2, i32::from(b'#'), 0, 10]);
    assert_eq!(&stated[36..46], b"big\0small\0");

    let regions: [(&str, u64, u64, usize); 13] = [
        ("big", 5_000, 5_000, 1),
        ("big", 1, 4_999, 0),
        ("big", 44_990, 45_010, 11),
        ("big", 45_001, 59_999, 0),
        ("big", 16_380, 16_390, 11),
        ("big", 59_990, 60_020, 11),
        ("big", 100_000, 130_000, 30_001),
        ("big", 129_999, last - 13, 2),
        ("big", last - 5, last, 6),
        ("big", 1, last, 70_026),
        ("small", 1, 1, 1),
        ("small", 2, 10, 2),
        ("small", 4, 100, 0),
    ];
    let expected = |sequence: &str, start: u64, end: u64| -> Vec<Row> {
        rows.iter()
            .filter(|(of, position, _)| of == sequence && (start..=end).contains(position))
            .cloned()
            .collect()
    };
    for (sequence, start, end, count) in regions {
        let wanted = expected(sequence, start, end);
        assert_eq!(wanted.len(), count, "{sequence}:{start}-{end}");
        assert_eq!(
            region_rows(&out, sequence, start, end),
            wanted,
            "{sequence}:{start}-{end}"
        );
    }

    let Some(program) = tabix_program() else {
        eprintln!("no tabix program on this machine: its reading of the index is not held");
        return;
    };
    let table = table.to_str().unwrap();
    assert_eq!(tabix_lines(&program, &["-l", table]), ["big", "small"]);
    for (sequence, start, end, _) in regions {
        let region = format!("{sequence}:{start}-{end}");
        let read: Vec<Row> = tabix_lines(&program, &[table, &region])
            .iter()
            .map(|line| row_of(line))
            .collect();
        assert_eq!(read, expected(sequence, start, end), "{region}");
    }
    // The index the program writes for the sequence in column 1 and the position in
    // column 2 states the header this crate writes.
    let copy = scratch.path().join("copy.txt.gz");
    fs::copy(table, &copy).unwrap();
    let copy = copy.to_str().unwrap();
    tabix_lines(&program, &["-f", "-s", "1", "-b", "2", "-e", "2", copy]);
    let made = tabix::fs::read(format!("{copy}.tbi")).expect("the program's index");
    assert_eq!(made.header(), index.header());
}
