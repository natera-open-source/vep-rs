// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary's options and files: `--fields`, a malformed command line, an input or a
//! tuple file that cannot be opened, a write that fails, a tuple file that is the input,
//! and the output of the normalizer beside the raw output.

mod common;

use common::{
    at, binary, exit_code, one_cell, row, row_of_gene, table, text, Scratch, ERROR, FIELD_MAP,
};
use std::fs;
use std::path::Path;

/// The gene-constraint column and a cell in its layout.
const GENE_CONSTRAINT: (&str, &str) = ("FV_GNOMAD_GENE", "GENE1|0.9981|0.2345|3.21|-0.45");

#[test]
fn fields_keeps_the_named_fields_and_counts_the_rest() {
    let scratch = Scratch::new();
    let (column, cell) = GENE_CONSTRAINT;
    let run = one_cell(&scratch, column, cell, &["--fields", "pLI_gene_value"]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(run.tuples(), at(&[("pLI_gene_value", "0.9981")]));
    assert!(run.stderr.contains("; by field: pLI_gene_value=1; "));
    assert!(run.stderr.contains("values outside --fields: 1;"));
    assert!(run
        .stderr
        .contains("dropped parts with no Ensembl VEP field: 2 "));
}

#[test]
fn fields_refuses_a_name_the_map_lacks() {
    let scratch = Scratch::new();
    let run = one_cell(
        &scratch,
        "FV_REVEL",
        "G|0.5",
        &["--fields", "REVEL,pLI_gene"],
    );
    assert_eq!(run.code, 1);
    assert_eq!(
        run.stderr,
        format!("{ERROR}--fields names no row of {FIELD_MAP}: pLI_gene\n")
    );
    assert!(!run.wrote());
}

fn fields_refuses_an_empty_value_and_an_empty_name_beside_a_comma(fields: &str, message: &str) {
    let scratch = Scratch::new();
    let run = one_cell(&scratch, "FV_REVEL", "G|0.5", &["--fields", fields]);
    assert_eq!(run.code, 1);
    assert_eq!(run.stderr, format!("{ERROR}{message}\n"));
    assert!(!run.wrote());
}

#[test]
fn fields_refuses_an_empty_value_and_an_empty_name_beside_a_comma_empty_value() {
    fields_refuses_an_empty_value_and_an_empty_name_beside_a_comma("", "--fields names no field");
}

#[test]
fn fields_refuses_an_empty_value_and_an_empty_name_beside_a_comma_between_two_commas() {
    fields_refuses_an_empty_value_and_an_empty_name_beside_a_comma(
        "REVEL,,am_class",
        "--fields names no field beside a comma: REVEL,,am_class",
    );
}

#[test]
fn fields_refuses_an_empty_value_and_an_empty_name_beside_a_comma_after_the_last_comma() {
    fields_refuses_an_empty_value_and_an_empty_name_beside_a_comma(
        "REVEL,",
        "--fields names no field beside a comma: REVEL,",
    );
}

fn fields_refuses_a_named_field_whose_column_the_header_lacks(fields: &str, absent: &str) {
    let scratch = Scratch::new();
    let (column, cell) = GENE_CONSTRAINT;
    let run = one_cell(&scratch, column, cell, &["--fields", fields]);
    assert_eq!(run.code, 1);
    assert_eq!(
        run.stderr,
        format!("{ERROR}the header lacks the column of a field --fields names: {absent}\n")
    );
    assert!(!run.wrote());
}

#[test]
fn fields_refuses_a_named_field_whose_column_the_header_lacks_every_named_field() {
    fields_refuses_a_named_field_whose_column_the_header_lacks("REVEL", "REVEL (FV_REVEL)");
}

#[test]
fn fields_refuses_a_named_field_whose_column_the_header_lacks_one_of_several() {
    fields_refuses_a_named_field_whose_column_the_header_lacks(
        "pLI_gene_value,am_class,REVEL",
        "REVEL (FV_REVEL), am_class (FV_ALPHAMISSENSE)",
    );
}

fn a_malformed_command_line_is_a_usage_error_and_writes_nothing(options: &[&str]) {
    let scratch = Scratch::new();
    let run = one_cell(&scratch, "FV_REVEL", "G|0.5", options);
    assert_eq!(run.code, 2);
    assert!(run.stderr.starts_with("usage: "), "{}", run.stderr);
    assert!(!run.wrote());
}

#[test]
fn a_malformed_command_line_is_a_usage_error_and_writes_nothing_fields_without_a_value() {
    a_malformed_command_line_is_a_usage_error_and_writes_nothing(&["--fields"]);
}

#[test]
fn a_malformed_command_line_is_a_usage_error_and_writes_nothing_map_without_a_value() {
    a_malformed_command_line_is_a_usage_error_and_writes_nothing(&["--field-map"]);
}

#[test]
fn a_malformed_command_line_is_a_usage_error_and_writes_nothing_unknown_option() {
    a_malformed_command_line_is_a_usage_error_and_writes_nothing(&[
        "--fields",
        "REVEL",
        "--verbose",
    ]);
}

#[test]
fn a_malformed_command_line_is_a_usage_error_and_writes_nothing_third_path() {
    a_malformed_command_line_is_a_usage_error_and_writes_nothing(&["a-third-path"]);
}

#[test]
fn an_input_that_cannot_be_read_is_refused() {
    let scratch = Scratch::new();
    let absent = text(&scratch.path("absent.tab"));
    let output = binary(&[&absent, &text(&scratch.path("tuples.tsv"))]);
    assert_eq!(exit_code(&output), 1);
    assert_eq!(
        String::from_utf8(output.stderr).unwrap(),
        format!("{ERROR}cannot read the input {absent}: No such file or directory\n")
    );
    assert!(!scratch.path("tuples.tsv").exists());
}

/// One row of the REVEL column under its header.
fn revel_table() -> String {
    table(&["FV_REVEL"], &[row("21:1000", "G", "ENST1", &["G|0.5"])])
}

#[test]
fn a_tuple_file_that_cannot_be_opened_is_refused() {
    let scratch = Scratch::new();
    let raw = scratch.write("raw.tab", revel_table());
    let absent = text(&scratch.path("absent").join("tuples.tsv"));
    let output = binary(&[&raw, &absent]);
    assert_eq!(exit_code(&output), 1);
    assert_eq!(
        String::from_utf8(output.stderr).unwrap(),
        format!("{ERROR}cannot write the tuple file {absent}: No such file or directory\n")
    );
}

#[test]
fn an_earlier_tuple_file_is_replaced_by_a_run_that_is_not_refused() {
    let scratch = Scratch::new();
    scratch.write("tuples.tsv", b"an earlier run\n");
    let run = one_cell(&scratch, "FV_REVEL", "G|0.5", &[]);
    assert_eq!(run.code, 0, "{}", run.stderr);
    assert_eq!(run.tuples(), at(&[("REVEL", "0.5")]));
}

#[test]
fn a_write_that_fails_ends_in_an_error_that_is_not_a_refusal() {
    // The test needs a device that refuses writes.
    if !Path::new("/dev/full").exists() {
        return;
    }
    let scratch = Scratch::new();
    let raw = scratch.write("raw.tab", revel_table());
    let output = binary(&[raw.as_str(), "/dev/full"]);
    assert_ne!(exit_code(&output), 0);
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(
        stderr.contains("OSError") && !stderr.contains(ERROR),
        "{stderr}"
    );
}

#[test]
fn a_tuple_file_that_is_the_input_is_refused_and_the_input_is_kept() {
    let scratch = Scratch::new();
    let written = revel_table();
    let raw = scratch.write("raw.tab", &written);
    let output = binary(&[&raw, &raw]);
    assert_eq!(exit_code(&output), 1);
    assert_eq!(
        String::from_utf8(output.stderr).unwrap(),
        format!("{ERROR}the tuple file {raw} is the input\n")
    );
    assert_eq!(fs::read_to_string(&raw).unwrap(), written);
}

#[test]
fn the_normalizer_output_gives_the_same_tuples_as_the_raw_output() {
    let scratch = Scratch::new();
    let columns = ["FV_REVEL", "FV_GNOMAD_GENE"];
    let rows = [
        row_of_gene("-", "21:100", "G", "-", &["G|0.5", "-"]),
        row(
            "21:200",
            "T",
            "ENST1",
            &["T|0.25", "GENE1|0.9981|0.2345|3.21|-0.45"],
        ),
    ];
    let written = table(&columns, &rows);
    let raw = scratch.write("raw.tab", &written);
    // The normalizer rewrites the Feature_type of a row without a feature, and nothing
    // else of this input.
    let rewritten = written.replace("\t-\t-\tTranscript\t", "\t-\t-\t-\t");
    assert_ne!(rewritten, written);
    let normalized = scratch.write("normalized.tab", &rewritten);
    let from_raw = binary(&[&raw, &text(&scratch.path("raw.tsv"))]);
    let from_normalized = binary(&[&normalized, &text(&scratch.path("normalized.tsv"))]);
    assert_eq!(exit_code(&from_raw), 0);
    assert_eq!(exit_code(&from_normalized), 0);
    assert_eq!(
        fs::read(scratch.path("raw.tsv")).unwrap(),
        fs::read(scratch.path("normalized.tsv")).unwrap()
    );
    assert_eq!(from_raw.stderr, from_normalized.stderr);
    assert!(String::from_utf8(from_raw.stderr)
        .unwrap()
        .starts_with("fastvep_plugin_tuples: wrote 4 tuples from 2 rows "));
}
