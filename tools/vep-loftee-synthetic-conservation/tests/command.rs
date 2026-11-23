// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The built binary: the lines it prints, the files it refuses to overwrite, its exit
//! codes and its error lines.

mod common;

use std::fs;
use std::path::Path;

use common::{
    database_path, index_path, new_cache, phylocsf, printed, run, small_cache, table_path,
    table_rows, table_text, totals, transcript, write_files, write_shard, Scratch,
};
use serde_json::json;
use sha2::{Digest, Sha256};
use vep_loftee_synthetic_conservation::{DATABASE_FILE, INDEX_FILE, TABLE_FILE};

fn sha256(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

/// The lines the database's digest is taken over, written here from the rows a query
/// of each table gives: the table's name and the row's values, the rows of a table in
/// the order of its columns.
fn database_lines(out: &Path) -> String {
    let mut totals = totals(out);
    totals.sort_by(|a, b| (&a.0, a.1, a.2).partial_cmp(&(&b.0, b.1, b.2)).unwrap());
    let mut scores = phylocsf(out);
    scores.sort_by(|a, b| {
        (&a.0, a.1, a.2, a.3)
            .partial_cmp(&(&b.0, b.1, b.2, b.3))
            .unwrap()
    });
    let mut lines = String::new();
    for (id, exon, total) in totals {
        lines.push_str(&format!("gerp_exons\t{id}\t{exon}\t{total}\n"));
    }
    for (id, exon, annotated, highest) in scores {
        lines.push_str(&format!(
            "phylocsf_data\t{id}\t{exon}\t{annotated}\t{highest}\n"
        ));
    }
    lines
}

#[test]
fn a_run_prints_its_seed_its_counts_and_one_line_a_file() {
    let scratch = Scratch::new("printed");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    let stdout = write_files(&cache, &out, 1, &[]);
    let lines = database_lines(&out);
    assert_eq!(lines.lines().count(), 18);
    assert!(lines.starts_with("gerp_exons\tMADE_N1\t1\t"));
    let expected = format!(
        "seed\t1\n\
         sequences\t2\n\
         transcripts\t5\n\
         exons\t9\n\
         exon_bases\t591\n\
         {TABLE_FILE}\trows=591\trows_sha256={}\n\
         {DATABASE_FILE}\trows=18\tgerp_exons=9\tphylocsf_data=9\trows_sha256={}\n",
        sha256(table_text(&out).as_bytes()),
        sha256(lines.as_bytes())
    );
    assert_eq!(stdout, expected);
}

#[test]
fn the_per_sequence_file_has_one_line_a_sequence() {
    let scratch = Scratch::new("per-sequence");
    let cache = small_cache(scratch.path());
    // A sequence without a transcript has a line and no row.
    fs::create_dir_all(cache.join("transcripts/bare")).unwrap();
    let out = scratch.path().join("out");
    let table = scratch.path().join("per-sequence.tsv");
    let stdout = write_files(
        &cache,
        &out,
        1,
        &["--per-sequence", table.to_str().unwrap()],
    );
    assert_eq!(printed(&stdout, "sequences"), ["3"]);
    assert_eq!(
        fs::read_to_string(&table).unwrap(),
        "sequence\ttranscripts\texons\texon_bases\n\
         first\t3\t6\t550\n\
         second\t2\t3\t41\n\
         bare\t0\t0\t0\n"
    );
    let rows = table_rows(&out);
    for (sequence, bases) in [("first", 550), ("second", 41), ("bare", 0)] {
        let held = rows.iter().filter(|(of, ..)| of == sequence).count();
        assert_eq!(held, bases, "{sequence}");
    }
    // The file is refused where it exists, as the three files are.
    let ended = run(&[
        "--cache",
        cache.to_str().unwrap(),
        "--seed",
        "1",
        "--out-dir",
        scratch.path().join("again").to_str().unwrap(),
        "--per-sequence",
        table.to_str().unwrap(),
    ]);
    assert_eq!(ended.code, Some(1));
    assert_eq!(
        ended.stderr,
        format!(
            "ERROR: [vep-loftee-synthetic-conservation] {} exists; --force overwrites it\n",
            table.display()
        )
    );
    assert!(!scratch.path().join("again").exists());
}

#[test]
fn an_existing_file_is_refused_without_force() {
    let scratch = Scratch::new("existing");
    let cache = small_cache(scratch.path());
    let out = scratch.path().join("out");
    write_files(&cache, &out, 1, &[]);
    let files = [table_path(&out), index_path(&out), database_path(&out)];
    let written: Vec<Vec<u8>> = files.iter().map(|file| fs::read(file).unwrap()).collect();
    let again = |out: &Path, more: &[&str]| {
        let mut args = vec![
            "--cache",
            cache.to_str().unwrap(),
            "--seed",
            "2",
            "--out-dir",
            out.to_str().unwrap(),
        ];
        args.extend_from_slice(more);
        run(&args)
    };

    let ended = again(&out, &[]);
    assert_eq!(ended.code, Some(1));
    assert_eq!(ended.stdout, "");
    assert_eq!(
        ended.stderr,
        format!(
            "ERROR: [vep-loftee-synthetic-conservation] {} exists; --force overwrites it\n",
            files[0].display()
        )
    );
    // Any one of the three files is enough, and the run writes nothing beside it.
    for (name, bytes) in [TABLE_FILE, INDEX_FILE, DATABASE_FILE].iter().zip(&written) {
        let alone = scratch.path().join(format!("alone-{name}"));
        fs::create_dir_all(&alone).unwrap();
        fs::write(alone.join(name), bytes).unwrap();
        let ended = again(&alone, &[]);
        assert_eq!(ended.code, Some(1), "{name}");
        assert!(
            ended.stderr.contains(&format!("{name} exists;")),
            "{}",
            ended.stderr
        );
        assert_eq!(fs::read_dir(&alone).unwrap().count(), 1, "{name}");
    }
    for (file, bytes) in files.iter().zip(&written) {
        assert_eq!(&fs::read(file).unwrap(), bytes);
    }

    let ended = again(&out, &["--force"]);
    assert_eq!(ended.code, Some(0), "{}", ended.stderr);
    assert_eq!(printed(&ended.stdout, "seed"), ["2"]);
    assert_ne!(&fs::read(&files[0]).unwrap(), &written[0]);
    assert_eq!(fs::read_dir(&out).unwrap().count(), 3);
}

#[test]
fn a_directory_under_a_file_s_name_is_refused_with_force() {
    let scratch = Scratch::new("directory");
    let cache = small_cache(scratch.path());
    for name in [TABLE_FILE, INDEX_FILE, DATABASE_FILE] {
        let out = scratch.path().join(format!("out-{name}"));
        fs::create_dir_all(out.join(name)).unwrap();
        let ended = run(&[
            "--cache",
            cache.to_str().unwrap(),
            "--seed",
            "1",
            "--out-dir",
            out.to_str().unwrap(),
            "--force",
        ]);
        assert_eq!(ended.code, Some(1), "{name}");
        assert_eq!(ended.stdout, "", "{name}");
        assert_eq!(
            ended.stderr,
            format!(
                "ERROR: [vep-loftee-synthetic-conservation] {} is a directory\n",
                out.join(name).display()
            )
        );
        // The run writes no file beside the directory.
        assert_eq!(fs::read_dir(&out).unwrap().count(), 1, "{name}");
    }

    // A directory under the name of the per-sequence file is refused before the run
    // writes a file.
    let taken = scratch.path().join("taken");
    fs::create_dir_all(&taken).unwrap();
    let out = scratch.path().join("out");
    let ended = run(&[
        "--cache",
        cache.to_str().unwrap(),
        "--seed",
        "1",
        "--out-dir",
        out.to_str().unwrap(),
        "--per-sequence",
        taken.to_str().unwrap(),
        "--force",
    ]);
    assert_eq!(ended.code, Some(1));
    assert_eq!(ended.stdout, "");
    assert_eq!(
        ended.stderr,
        format!(
            "ERROR: [vep-loftee-synthetic-conservation] {} is a directory\n",
            taken.display()
        )
    );
    assert!(!out.exists());
}

#[test]
fn a_transcript_without_exons_fails_the_run() {
    let scratch = Scratch::new("without-exons");
    let cache = new_cache(scratch.path(), &["first", "ring"]);
    write_shard(
        &cache,
        "first",
        "1-1000000.json",
        &[transcript("MADE_T1", 1, "protein_coding", &[(101, 200)])],
    );
    // A transcript stored across the origin of a circular sequence, one whose exon
    // ends outside the transcript's own span, and one the reader gives an exon for.
    write_shard(
        &cache,
        "ring",
        "1-1000000.json",
        &[
            json!({"stable_id": "MADE_RING", "start": -50, "end": 120, "strand": 1,
                   "biotype": "protein_coding", "exons": [{"start": -50, "end": 120}]}),
            json!({"stable_id": "MADE_OUT", "start": 150, "end": 200, "strand": 1,
                   "biotype": "lncRNA", "exons": [{"start": 150, "end": 260}]}),
            transcript("MADE_R1", -1, "protein_coding", &[(301, 330)]),
        ],
    );
    let out = scratch.path().join("out");
    let args = [
        "--cache",
        cache.to_str().unwrap(),
        "--seed",
        "1",
        "--out-dir",
        out.to_str().unwrap(),
    ];
    let ended = run(&args);
    assert_eq!(ended.code, Some(1));
    assert_eq!(ended.stdout, "");
    assert_eq!(
        ended.stderr,
        "ERROR: [vep-loftee-synthetic-conservation] transcript records the cache reader gives \
         no exon for: 2, the first MADE_RING on sequence ring; \
         --skip-transcripts-without-exons writes the files without them\n"
    );
    // The failed run leaves no file, under a file's own name or another.
    assert_eq!(fs::read_dir(&out).unwrap().count(), 0);

    let mut skipping = args.to_vec();
    skipping.push("--skip-transcripts-without-exons");
    let ended = run(&skipping);
    assert_eq!(ended.code, Some(0), "{}", ended.stderr);
    let lines: Vec<&str> = ended.stdout.lines().collect();
    assert_eq!(
        lines[..6],
        [
            "seed\t1",
            "sequences\t2",
            "transcripts\t4",
            "transcripts_without_exons\t2",
            "exons\t2",
            "exon_bases\t130",
        ]
    );
    let keys: Vec<(String, i64)> = totals(&out)
        .into_iter()
        .map(|(id, exon, _)| (id, exon))
        .collect();
    assert_eq!(keys, [("MADE_T1".to_owned(), 1), ("MADE_R1".to_owned(), 1)]);
    assert_eq!(table_rows(&out).len(), 130);

    // A cache whose every record has an exon prints the line with a zero under the
    // option and no such line without it.
    let whole = small_cache(&scratch.path().join("whole"));
    let whole_out = scratch.path().join("whole-out");
    let stdout = write_files(&whole, &whole_out, 1, &["--skip-transcripts-without-exons"]);
    assert_eq!(printed(&stdout, "transcripts_without_exons"), ["0"]);
    let plain = write_files(&whole, &scratch.path().join("plain-out"), 1, &[]);
    assert!(!plain.contains("transcripts_without_exons"));
}

#[test]
fn a_command_line_the_tool_cannot_read_ends_with_the_usage() {
    let scratch = Scratch::new("usage");
    let out = scratch.path().join("out");
    let out = out.to_str().unwrap();
    for (args, what) in [
        (&[][..], "--cache is required"),
        (
            &["--cache", "c", "--out-dir", out][..],
            "--seed is required",
        ),
        (
            &["--cache", "c", "--seed", "1"][..],
            "--out-dir is required",
        ),
        (
            &["--cache", "c", "--seed", "one", "--out-dir", out][..],
            "--seed takes a whole number",
        ),
        (
            &[
                "--cache",
                "c",
                "--seed",
                "1",
                "--out-dir",
                out,
                "--threads",
                "0",
            ][..],
            "--threads takes a whole number above zero",
        ),
        (
            &["--cache", "c", "--seed", "1", "--out-dir", out, "--records"][..],
            "no option --records",
        ),
        (
            &["--cache", "c", "--seed", "1", "--out-dir", out, "--force=1"][..],
            "--force takes no value",
        ),
        (
            &["--cache", "c", "--seed", "1", "--seed", "2"][..],
            "--seed is given twice",
        ),
        (&["--cache", "c", "--seed"][..], "--seed takes a value"),
        (&["c"][..], "c is no option"),
        (&["--version", "--force"][..], "no option --version"),
    ] {
        let ended = run(args);
        assert_eq!(ended.code, Some(2), "{args:?}");
        assert_eq!(ended.stdout, "", "{args:?}");
        let mut lines = ended.stderr.lines();
        assert_eq!(
            lines.next(),
            Some(format!("ERROR: [vep-loftee-synthetic-conservation] {what}").as_str()),
            "{args:?}"
        );
        assert_eq!(lines.next(), Some("usage:"), "{args:?}");
    }
    assert!(!Path::new(out).exists());

    let ended = run(&["--version"]);
    assert_eq!(ended.code, Some(0));
    assert_eq!(ended.stdout, "vep-loftee-synthetic-conservation 0.1.0\n");
    let ended = run(&["--help"]);
    assert_eq!(ended.code, Some(0));
    assert!(ended.stdout.starts_with("usage:\n"));
    assert!(ended.stdout.contains("--skip-transcripts-without-exons"));
}

#[test]
fn a_cache_that_cannot_be_read_fails_the_run() {
    let scratch = Scratch::new("no-cache");
    let out = scratch.path().join("out");
    let missing = scratch.path().join("missing");
    let ended = run(&[
        "--cache",
        missing.to_str().unwrap(),
        "--seed",
        "1",
        "--out-dir",
        out.to_str().unwrap(),
    ]);
    assert_eq!(ended.code, Some(1));
    assert_eq!(ended.stdout, "");
    assert_eq!(ended.stderr.lines().count(), 1, "{}", ended.stderr);
    assert!(ended
        .stderr
        .starts_with("ERROR: [vep-loftee-synthetic-conservation] "));
    assert!(!out.exists());

    // A region file that is no JSON fails the run and leaves no file.
    let cache = small_cache(scratch.path());
    fs::write(cache.join("transcripts/second/1-1000000.json"), "{").unwrap();
    let ended = run(&[
        "--cache",
        cache.to_str().unwrap(),
        "--seed",
        "1",
        "--out-dir",
        out.to_str().unwrap(),
    ]);
    assert_eq!(ended.code, Some(1));
    assert_eq!(ended.stderr.lines().count(), 1, "{}", ended.stderr);
    assert!(ended.stderr.contains("transcripts/second/1-1000000.json"));
    assert_eq!(fs::read_dir(&out).unwrap().count(), 0);
}

#[test]
fn an_exon_above_the_last_position_of_an_index_fails_the_run() {
    let scratch = Scratch::new("above");
    let cache = new_cache(scratch.path(), &["far"]);
    write_shard(
        &cache,
        "far",
        "1-1000000.json",
        &[transcript(
            "MADE_FAR",
            1,
            "lncRNA",
            &[(536_870_902, 536_870_912)],
        )],
    );
    let out = scratch.path().join("out");
    let ended = run(&[
        "--cache",
        cache.to_str().unwrap(),
        "--seed",
        "1",
        "--out-dir",
        out.to_str().unwrap(),
    ]);
    assert_eq!(ended.code, Some(1));
    assert_eq!(
        ended.stderr,
        "ERROR: [vep-loftee-synthetic-conservation] sequence far has an exon base at 536870912, \
         above the 536870911 a tabix index reaches\n"
    );
    assert_eq!(fs::read_dir(&out).unwrap().count(), 0);
}
