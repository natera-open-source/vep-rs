// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary over tables made at run time, each held to bytes laid out row by row: a
//! table of many runs and batches under every thread bound, quoted fields that cross the
//! runs, and fields at and past the reader's limit.

mod common;

use common::{assert_clean_exit, assert_refused, run_bound, text, Scratch, BOUNDS};
use std::fmt::Write;
use std::fs;

const HEADER: &str = "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set";

/// A table of `variants` variants in blocks of fifty, the reference-only rows of a block
/// before its scored-only rows, and the table and the summary the tool writes for it. The
/// shape of a variant follows its number: a pair that differs in the end of its Location, a
/// pair that differs in its consequence set, a pair under two Features, a row the reference
/// alone wrote and a row the scored engine alone wrote.
fn blocks(variants: usize, note: &str) -> (String, String, String) {
    let mut table = format!("{HEADER}\n");
    let mut classed = format!("{HEADER}\tmismatch_type\r\n");
    let mut counts = [0usize; 5];
    let mut row = |fields: [&str; 4], class: usize| {
        let [source, location, feature, set] = fields;
        let file = if class.is_multiple_of(2) {
            "even.txt"
        } else {
            "odd.txt"
        };
        let line = format!("{file}\t{source}\t{location}\t{note}\t{feature}\tTranscript\t{set}");
        let name = [
            "coordinate_only",
            "consequence_swap",
            "transcript_selection",
            "true_missing",
            "true_extra",
        ][class];
        writeln!(table, "{line}").unwrap();
        write!(classed, "{line}\t{name}\r\n").unwrap();
        counts[class] += 1;
    };
    for block in 0..variants.div_ceil(50) {
        let numbers = block * 50..((block + 1) * 50).min(variants);
        for n in numbers.clone() {
            match n % 5 {
                0 => row(
                    ["missing_in_rust", &format!("1:{n}-{}", n + 1), "T", "x"],
                    0,
                ),
                1 => row(["missing_in_rust", &format!("2:{n}"), "T", "x"], 1),
                2 => row(["missing_in_rust", &format!("3:{n}"), "Ta", "x"], 2),
                3 => row(["missing_in_rust", &format!("4:{n}"), "T", "x"], 3),
                _ => {}
            }
        }
        for n in numbers {
            match n % 5 {
                0 => row(["extra_in_rust", &format!("1:{n}-{}", n + 2), "T", "x"], 0),
                1 => row(["extra_in_rust", &format!("2:{n}"), "T", "y"], 1),
                2 => row(["extra_in_rust", &format!("3:{n}"), "Tb", "x"], 2),
                4 => row(["extra_in_rust", &format!("5:{n}"), "T", "x"], 4),
                _ => {}
            }
        }
    }
    let total: usize = counts.iter().sum();
    let mut summary = format!(
        "{rule}\nDISCORDANT CLASSIFICATION SUMMARY\n{rule}\n\nTotal discordants: {total}\n",
        rule = "=".repeat(70)
    );
    // The variants are a multiple of five: the three paired classes hold a quarter of the
    // rows each and are met in this order, the two unpaired ones an eighth each.
    summary.push_str("\nBy mismatch type:\n");
    for (name, share) in [
        ("coordinate_only", " 25.0"),
        ("consequence_swap", " 25.0"),
        ("transcript_selection", " 25.0"),
        ("true_missing", " 12.5"),
        ("true_extra", " 12.5"),
    ] {
        writeln!(
            summary,
            "  {name:25} {:>8} ({share}%)",
            total / if share == " 25.0" { 4 } else { 8 }
        )
        .unwrap();
    }
    summary.push_str("\nBy mismatch type + direction:\n");
    for (name, source) in [
        ("coordinate_only", "missing_in_rust"),
        ("consequence_swap", "missing_in_rust"),
        ("transcript_selection", "missing_in_rust"),
        ("true_missing", "missing_in_rust"),
        ("coordinate_only", "extra_in_rust"),
        ("consequence_swap", "extra_in_rust"),
        ("transcript_selection", "extra_in_rust"),
        ("true_extra", "extra_in_rust"),
    ] {
        writeln!(summary, "  {name:25} {source:12} {:>8}", total / 8).unwrap();
    }
    summary.push_str("\nBy file + mismatch type (top 30):\n");
    for (file, name, share) in [
        ("even.txt", "coordinate_only", " 25.0"),
        ("odd.txt", "consequence_swap", " 25.0"),
        ("even.txt", "transcript_selection", " 25.0"),
        ("odd.txt", "true_missing", " 12.5"),
        ("even.txt", "true_extra", " 12.5"),
    ] {
        writeln!(
            summary,
            "  {file:35} {name:25} {:>8} ({share}%)",
            total / if share == " 25.0" { 4 } else { 8 }
        )
        .unwrap();
    }
    (table, classed, summary)
}

#[test]
fn a_table_of_many_runs_is_classed_alike_under_every_thread_bound() {
    let scratch = Scratch::new("many-runs");
    let (table, classed, summary) = blocks(100_000, "A");
    assert!(table.len() > 5 << 20, "several runs of text");
    fs::write(scratch.path("in.tsv"), &table).unwrap();
    for bound in BOUNDS {
        let output = run_bound(&scratch.0, &["in.tsv", "-o", "out.tsv"], bound);
        assert_clean_exit(&output);
        assert_eq!(text(&output.stdout), summary, "{bound:?}");
        assert!(
            fs::read(scratch.path("out.tsv")).unwrap() == classed.as_bytes(),
            "{bound:?}"
        );
    }
}

#[test]
fn line_endings_of_every_kind_are_read_across_the_runs() {
    let scratch = Scratch::new("endings");
    let (table, classed, summary) = blocks(40_000, "A");
    // CR LF, then a lone CR, then LF, row after row.
    let mut mixed = String::with_capacity(table.len() + table.len() / 100);
    for (index, line) in table.lines().enumerate() {
        mixed.push_str(line);
        mixed.push_str(["\r\n", "\r", "\n"][index % 3]);
    }
    fs::write(scratch.path("in.tsv"), mixed.trim_end()).unwrap();
    for bound in [Some("1"), Some("3"), None] {
        let output = run_bound(&scratch.0, &["in.tsv", "-o", "out.tsv"], bound);
        assert_clean_exit(&output);
        assert_eq!(text(&output.stdout), summary, "{bound:?}");
        assert!(
            fs::read(scratch.path("out.tsv")).unwrap() == classed.as_bytes(),
            "{bound:?}"
        );
    }
}

#[test]
fn quoted_fields_that_cross_the_runs_are_written_between_quotes_again() {
    let scratch = Scratch::new("quoted-runs");
    // The Allele of every row holds a tab, a line ending and a quote: it is read from
    // between quotes and is written between quotes, its quote doubled.
    let (table, classed, summary) = blocks(30_000, "\"an\tallele\nof \"\"two\"\" lines\"");
    assert!(table.len() > 3 << 20, "several runs of text");
    fs::write(scratch.path("in.tsv"), &table).unwrap();
    for bound in [Some("1"), Some("2"), Some("5"), None] {
        let output = run_bound(&scratch.0, &["in.tsv", "-o", "out.tsv"], bound);
        assert_clean_exit(&output);
        assert_eq!(text(&output.stdout), summary, "{bound:?}");
        assert!(
            fs::read(scratch.path("out.tsv")).unwrap() == classed.as_bytes(),
            "{bound:?}"
        );
    }
}

#[test]
fn a_field_holds_at_most_131072_code_points() {
    let scratch = Scratch::new("limit");
    let row = |feature: &str| {
        format!("{HEADER}\no.txt\textra_in_rust\t1:10\tA\t{feature}\tTranscript\tx\n")
    };
    for (unit, name) in [("x", "one byte"), ("\u{e9}", "two bytes")] {
        let at_limit = unit.repeat(131_072);
        fs::write(scratch.path("in.tsv"), row(&at_limit)).unwrap();
        let output = run_bound(&scratch.0, &["in.tsv", "-o", "out.tsv"], Some("2"));
        assert_clean_exit(&output);
        assert_eq!(
            text(&fs::read(scratch.path("out.tsv")).unwrap()),
            format!(
                "{HEADER}\tmismatch_type\r\no.txt\textra_in_rust\t1:10\tA\t{at_limit}\tTranscript\tx\ttrue_extra\r\n"
            ),
            "{name}"
        );
        fs::remove_file(scratch.path("out.tsv")).unwrap();

        for past_limit in [format!("{at_limit}{unit}"), format!("\"{at_limit}{unit}\"")] {
            fs::write(scratch.path("in.tsv"), row(&past_limit)).unwrap();
            let output = run_bound(&scratch.0, &["in.tsv", "-o", "out.tsv"], Some("2"));
            assert_refused(&output, "in.tsv: field larger than field limit (131072)");
            assert!(output.stdout.is_empty(), "{name}");
            assert!(!scratch.path("out.tsv").exists(), "{name}");
        }
    }
}
