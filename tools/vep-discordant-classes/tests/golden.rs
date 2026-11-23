// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over every fixture table and compares the classified table, the summary
//! and the exit code with the goldens, which hold the tool's contract, under every thread
//! bound.

mod common;

use common::{assert_clean_exit, assert_refused, bytes, fixture, run_bound, text, Scratch, BOUNDS};
use std::fs;

/// The fixtures and the exit code of each, from `tests/fixtures/cases.tsv`.
fn cases() -> Vec<(String, i32)> {
    text(&bytes("cases.tsv"))
        .lines()
        .map(|line| {
            let (file, code) = line.split_once('\t').expect("a fixture and its exit code");
            (file.to_string(), code.parse().expect("an exit code"))
        })
        .collect()
}

/// The name a fixture's goldens carry: its file name up to the first dot.
fn name_of(file: &str) -> &str {
    file.split('.').next().unwrap_or(file)
}

#[test]
fn every_fixture_matches_its_goldens_under_every_thread_bound() {
    let cases = cases();
    assert_eq!(cases.len(), 24);
    let scratch = Scratch::new("goldens");
    for (file, code) in &cases {
        let name = name_of(file);
        let table = scratch.path(&format!("{name}.out.tsv"));
        let expected_table = fixture(&format!("expected_{name}.tsv"));
        for bound in BOUNDS {
            let _ = fs::remove_file(&table);
            let output = run_bound(
                &scratch.0,
                &[fixture(file).as_os_str(), "-o".as_ref(), table.as_os_str()],
                bound,
            );
            let at = format!("{file} under {bound:?}: {}", text(&output.stderr));
            assert_eq!(output.status.code(), Some(*code), "{at}");
            assert_eq!(
                output.stdout,
                bytes(&format!("expected_{name}.stdout")),
                "{at}"
            );
            assert_eq!(table.exists(), expected_table.exists(), "{at}");
            if expected_table.exists() {
                assert_eq!(
                    fs::read(&table).unwrap(),
                    fs::read(&expected_table).unwrap(),
                    "{at}"
                );
            }
            if *code == 0 {
                assert_clean_exit(&output);
            } else {
                let error = text(&output.stderr);
                assert!(
                    error.starts_with("ERROR: [vep-discordant-classes] ")
                        && error.ends_with('\n')
                        && error.lines().count() == 1,
                    "{at}"
                );
            }
        }
    }
}

#[test]
fn every_golden_belongs_to_a_fixture_and_every_whole_run_has_both() {
    let cases = cases();
    let mut goldens: Vec<String> = fs::read_dir(fixture(""))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|file| file.starts_with("expected_"))
        .collect();
    goldens.sort();
    let mut expected = Vec::new();
    for (file, code) in &cases {
        let name = name_of(file);
        let table = format!("expected_{name}.tsv");
        let summary = format!("expected_{name}.stdout");
        if *code == 0 {
            assert!(
                goldens.contains(&table) && goldens.contains(&summary),
                "{file}"
            );
        }
        expected.extend([table, summary].into_iter().filter(|g| goldens.contains(g)));
    }
    expected.sort();
    assert_eq!(goldens, expected);
}

/// Runs the binary over one refused fixture and holds its error line.
fn assert_fixture_refused(file: &str, detail: &str) {
    let scratch = Scratch::new(&format!("refused-{file}"));
    let table = scratch.path("out.tsv");
    let input = fixture(file);
    let output = run_bound(
        &scratch.0,
        &[input.as_os_str(), "-o".as_ref(), table.as_os_str()],
        Some("2"),
    );
    assert_refused(&output, &format!("{}: {detail}", input.display()));
}

#[test]
fn a_table_without_its_columns_is_refused_by_the_names_it_lacks() {
    let schema = "Unsupported discordant TSV schema; missing required columns:";
    assert_fixture_refused("missing_columns.tsv", &format!("{schema} allele, feature"));
    assert_fixture_refused("byte_order_mark.tsv", &format!("{schema} file_name"));
    assert_fixture_refused(
        "blank_first_line.tsv",
        &format!("{schema} allele, consequence_set, feature, file_name, location, source"),
    );
    assert_fixture_refused("empty.tsv", "no header row");
}

#[test]
fn a_table_that_is_not_utf_8_is_refused_by_the_offset_of_its_first_such_byte() {
    assert_fixture_refused("not_utf8.tsv", "the byte at offset 143 is not UTF-8");
    assert_fixture_refused("gzipped.tsv.gz", "the byte at offset 1 is not UTF-8");
}

#[test]
fn a_row_that_ends_before_its_location_is_refused_by_its_number() {
    assert_fixture_refused("no_location.tsv", "row 2 ends before its location column");
}

#[test]
fn a_summary_that_cannot_print_a_row_names_the_column_the_row_ends_before() {
    let ends = "the summary ends at a row that ends before its";
    assert_fixture_refused(
        "classified_short_row.tsv",
        &format!("{ends} mismatch_type column"),
    );
    assert_fixture_refused("summary_no_source.tsv", &format!("{ends} source column"));
    assert_fixture_refused("summary_no_file.tsv", &format!("{ends} file_name column"));
}

#[test]
fn the_goldens_hold_what_the_fixtures_are_made_to_show() {
    let table = |name: &str| text(&bytes(&format!("expected_{name}.tsv")));
    let summary = |name: &str| text(&bytes(&format!("expected_{name}.stdout")));
    // Every class, and rows that end in CR LF.
    let tuple = table("tuple");
    for class in [
        "coordinate_only",
        "consequence_swap",
        "transcript_selection",
        "true_missing",
        "true_extra",
    ] {
        assert!(tuple.contains(&format!("\t{class}\r\n")), "{class}");
    }
    assert_eq!(tuple.matches("\r\n").count(), 26);
    // Starts pair across chromosomes.
    assert!(table("locations").starts_with(
        "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\tmismatch_type\r\n\
         o.txt\tmissing_in_rust\tchr1:100\tA\tT1\tTranscript\tx\tcoordinate_only\r\n\
         o.txt\textra_in_rust\tchr2:100-150\tA\tT1\tTranscript\tx\tcoordinate_only\r\n"
    ));
    // A short row is filled, a long one cut, a field with a tab or a quote is quoted.
    assert!(table("short_rows")
        .contains("o.txt\tmissing_in_rust\t1:100-200\t\t\t\t\tcoordinate_only\r\n"));
    assert!(table("short_rows")
        .contains("o.txt\tmissing_in_rust\t4:1\tG\tT4\tTranscript\ty\ttrue_missing\r\n"));
    assert!(table("quoted").contains("\t\"an \"\"allele\"\"\"\t\"T\n2\"\t"));
    // A table classed before keeps its values, a pair of the second pairing excepted.
    assert!(table("classified").contains(
        "o.txt\tmissing_in_rust\t2:5\tA\tT2\tTranscript\tx\tconsequence_swap\tconsequence_swap\r\n"
    ));
    assert!(table("classified").contains(
        "o.txt\textra_in_rust\t4:7\tC\tT3\tTranscript\tx\ta label of its own\ta label of its own\r\n"
    ));
    // Half-way percentages go to the even tenth, and the third block holds thirty rows.
    let wide = summary("summary");
    assert!(
        wide.contains("              21 (  5.2%)\n")
            && wide.contains("              19 (  4.8%)\n")
    );
    let third_block = wide.split("(top 30):\n").nth(1).unwrap();
    assert_eq!(third_block.lines().count(), 30);
    assert!(summary("header_only").ends_with("Total discordants: 0\n\nBy mismatch type:\n\nBy mismatch type + direction:\n\nBy file + mismatch type (top 30):\n"));
    // A summary that stops keeps the lines before the row it cannot print.
    assert!(summary("summary_no_source")
        .ends_with("  true_missing              missing_in_rust        2\n"));
}
