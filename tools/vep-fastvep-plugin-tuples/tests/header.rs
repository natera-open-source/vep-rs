// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary over column headers: which line is the header, a source column without a
//! map row, and the headers a run refuses, each leaving an earlier tuple file as it was.

mod common;

use common::{at, row, run, run_text, table, tuple, Scratch, ERROR};
use std::fs;

#[test]
fn the_header_is_the_last_single_hash_line_before_the_first_row() {
    let scratch = Scratch::new();
    let text = [
        "#fabricated comment\n".to_string(),
        table(&["FV_PRIMATEAI", "FV_REVEL"], &[]),
        "## FV_REVEL\ta meta line, never the header\n".to_string(),
        row("21:1000", "G", "ENST1", &["-", "G|0.5"]),
        "\n#".to_string(),
        row("21:1500", "A", "ENST9", &["A|0.1", "A|0.9"]),
        "\n".to_string(),
        row("21:2000", "T", "ENST2", &["T|0.7312", "-"]),
        "\n".to_string(),
    ]
    .concat();
    let run = run_text(&scratch, &text, &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(
        run.tuples(),
        [
            tuple(["21:1000", "G", "ENST1", "REVEL", "0.5"]),
            tuple(["21:2000", "T", "ENST2", "PrimateAI", "0.7312"]),
        ]
    );
}

#[test]
fn a_source_column_without_a_map_row_is_named_in_the_summary() {
    let scratch = Scratch::new();
    let rows = [row("21:1000", "G", "ENST1", &["G|1.234", "G|0.5"])];
    let run = run(&scratch, &["FV_PHYLOP", "FV_REVEL"], &rows, &[]);
    assert_eq!(run.tuples(), at(&[("REVEL", "0.5")]));
    assert!(run.stderr.ends_with("unmapped columns: FV_PHYLOP\n"));
}

fn a_header_naming_a_mapped_column_more_than_once_is_refused_with_its_positions(
    columns: &[&str],
    options: &[&str],
    named: &str,
) {
    let scratch = Scratch::new();
    scratch.write("tuples.tsv", b"an earlier run\n");
    let cells: Vec<String> = (0..columns.len())
        .map(|place| format!("G|0.{place}"))
        .collect();
    let cells: Vec<&str> = cells.iter().map(String::as_str).collect();
    let run = run(
        &scratch,
        columns,
        &[row("21:1000", "G", "ENST1", &cells)],
        options,
    );
    assert_eq!(run.code, 1);
    assert_eq!(
        run.stderr,
        format!("{ERROR}the header names a mapped column more than once: {named}\n")
    );
    assert_eq!(run.tuple_bytes(), b"an earlier run\n");
}

#[test]
fn a_header_naming_a_mapped_column_more_than_once_is_refused_with_its_positions_one_column_twice() {
    a_header_naming_a_mapped_column_more_than_once_is_refused_with_its_positions(
        &["FV_REVEL", "FV_REVEL"],
        &[],
        "FV_REVEL (columns 15 and 16)",
    );
}

#[test]
fn a_header_naming_a_mapped_column_more_than_once_is_refused_with_its_positions_two_columns_and_an_unmapped_one(
) {
    a_header_naming_a_mapped_column_more_than_once_is_refused_with_its_positions(
        &[
            "FV_PHYLOP",
            "FV_REVEL",
            "FV_PHYLOP",
            "FV_GNOMAD_GENE",
            "FV_REVEL",
            "FV_GNOMAD_GENE",
            "FV_REVEL",
        ],
        &[],
        "FV_REVEL (columns 16, 19 and 21), FV_GNOMAD_GENE (columns 18 and 20)",
    );
}

#[test]
fn a_header_naming_a_mapped_column_more_than_once_is_refused_with_its_positions_a_column_fields_leaves_out(
) {
    a_header_naming_a_mapped_column_more_than_once_is_refused_with_its_positions(
        &["FV_GNOMAD_GENE", "FV_REVEL", "FV_REVEL"],
        &["--fields", "pLI_gene_value"],
        "FV_REVEL (columns 16 and 17)",
    );
}

fn an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing(
    text: &str,
    message: &str,
) {
    let scratch = Scratch::new();
    scratch.write("tuples.tsv", b"an earlier run\n");
    let run = run_text(&scratch, text, &[]);
    assert_eq!(run.code, 1);
    assert!(
        run.stderr.starts_with(&format!("{ERROR}{message}")),
        "{}",
        run.stderr
    );
    assert_eq!(run.stderr.matches('\n').count(), 1);
    assert_eq!(fs::read(&run.tuple_file).unwrap(), b"an earlier run\n");
}

#[test]
fn an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing_no_mapped_column() {
    an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing(
        &table(&[], &[row("21:1000", "G", "ENST1", &[])]),
        "the field map names none of the header's columns: Uploaded_variation, Location,",
    );
}

#[test]
fn an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing_no_mapped_column_no_rows()
{
    an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing(
        &table(&["FV_PHYLOP"], &[]),
        "the field map names none of the header's columns: ",
    );
}

#[test]
fn an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing_no_header() {
    an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing(
        &(row("21:1000", "G", "ENST1", &["G|0.5"]) + "\n"),
        "the input has no column header line",
    );
}

#[test]
fn an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing_empty_input() {
    an_input_whose_columns_cannot_be_placed_is_refused_and_writes_nothing(
        "",
        "the input has no column header line",
    );
}
