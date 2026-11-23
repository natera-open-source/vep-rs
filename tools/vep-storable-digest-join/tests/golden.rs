// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The binary against reference bytes.
//!
//! `tests/fixtures/inputs/` holds digest files and parity tables, and
//! `tests/fixtures/expected/<case>/` one reference run over them, made by CPython 3.14:
//! the arguments (`args`, one a line), standard output, standard error, the exit code,
//! and under `tree/` every file the run wrote or changed. A file the run names with a
//! character outside printable ASCII, a double quote or a backslash lies under `files/`
//! by a plain name, and `files.tsv` beside it spells the path the run wrote it to. Each
//! test runs the binary with the same arguments in a copy of the inputs and asks for the
//! same bytes: the two streams, the exit code, and the directory as the run leaves it.
//! Standard error is held as the tool prints it where the reference is `argparse`'s own
//! text or a traceback (the README's "Outside the contract"), and a sweep id is held as
//! the marker `{sweep-id}`, in whose place the test puts an id of the right shape.

mod common;

use common::{copy_tree, fixtures, run_in, tree, written, Scratch};

fn lines(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// The marker the fixtures hold in the place of a sweep id.
const SWEEP_ID_MARKER: &str = "{sweep-id}";

/// A sweep id of the shape the tool takes: `sweep-`, eight digits, `T`, six digits, `Z`.
fn sweep_id() -> String {
    format!("sweep-{}T{}Z", "12345678", "123456")
}

fn case(id: &str) {
    let expected = fixtures().join("expected").join(id);
    let arguments: Vec<String> = std::fs::read_to_string(expected.join("args"))
        .expect("args")
        .lines()
        .map(|line| line.replace(SWEEP_ID_MARKER, &sweep_id()))
        .collect();
    let argument_refs: Vec<&str> = arguments.iter().map(String::as_str).collect();
    let scratch = Scratch::new(id);
    copy_tree(&fixtures().join("inputs"), scratch.path());
    let run = run_in(scratch.path(), &argument_refs);

    let exit: i32 = std::fs::read_to_string(expected.join("exit"))
        .expect("exit")
        .trim()
        .parse()
        .expect("an exit code");
    assert_eq!(
        run.status.code(),
        Some(exit),
        "{id}: exit code; stderr: {}",
        lines(&run.stderr)
    );
    let stdout = std::fs::read(expected.join("stdout")).expect("stdout");
    assert_eq!(lines(&run.stdout), lines(&stdout), "{id}: standard output");
    assert_eq!(run.stdout, stdout, "{id}: standard output bytes");
    let stderr = std::fs::read(expected.join("stderr")).expect("stderr");
    assert_eq!(lines(&run.stderr), lines(&stderr), "{id}: standard error");
    assert_eq!(run.stderr, stderr, "{id}: standard error bytes");

    // The directory after the run: the inputs, with what the reference run wrote or
    // changed in their place.
    let mut wanted = tree(&fixtures().join("inputs"));
    wanted.extend(written(&expected).into_iter().map(|(name, bytes)| {
        let text = String::from_utf8(bytes).expect("UTF-8");
        (
            name,
            text.replace(SWEEP_ID_MARKER, &sweep_id()).into_bytes(),
        )
    }));
    let left = tree(scratch.path());
    let names = |files: &std::collections::BTreeMap<String, Vec<u8>>| {
        files.keys().cloned().collect::<Vec<_>>()
    };
    assert_eq!(names(&left), names(&wanted), "{id}: files after the run");
    for (name, bytes) in &wanted {
        assert_eq!(
            lines(&left[name]),
            lines(bytes),
            "{id}: {name} after the run"
        );
        assert!(left[name] == *bytes, "{id}: the bytes of {name}");
    }
}

macro_rules! cases {
    ($($name:ident => $id:literal,)*) => {
        $(
            #[test]
            fn $name() {
                case($id);
            }
        )*

        #[test]
        fn every_expected_directory_has_its_test() {
            let mut listed: Vec<&str> = vec![$($id),*];
            listed.sort_unstable();
            let mut present: Vec<String> = std::fs::read_dir(fixtures().join("expected"))
                .expect("expected directory")
                .map(|entry| entry.expect("entry").file_name().to_string_lossy().into_owned())
                .collect();
            present.sort();
            assert_eq!(present, listed);
        }
    };
}

cases! {
    // The four lists, the summary and the line on standard output.
    two_files_are_joined_on_the_id => "plain",
    the_output_directory_is_made_with_its_parents => "nested-out-dir",
    line_ends_of_every_kind_end_a_line => "line-ends",
    one_file_under_two_names_is_identical => "identical",
    one_version_for_both_files_shares_its_keys_and_lists => "same-version",
    two_empty_files_are_content_identical => "empty-files",
    an_empty_file_leaves_every_id_on_one_side => "one-empty",
    versions_and_ids_outside_ascii => "unicode",
    counts_past_a_thousand_are_grouped => "thousands",
    options_shortened_and_with_an_equals_sign => "prefixes",
    // The parity table.
    the_row_of_the_assembly_gains_the_columns => "parity",
    a_sweep_id_is_written_to_the_row => "parity-sweep",
    columns_the_table_holds_are_written_over => "parity-rerun",
    quoted_short_long_and_empty_records_of_a_table => "parity-shapes",
    an_assembly_with_a_comma_and_a_quote_left_open => "parity-last-row",
    two_columns_of_one_name => "parity-duplicate-names",
    a_table_of_one_column_keeps_its_empty_record_quoted => "parity-one-column",
    a_table_path_is_normalised_before_its_suffix_is_replaced => "parity-no-suffix",
    two_rows_of_the_assembly_fail_the_run => "parity-two-rows",
    no_row_of_the_assembly_fails_the_run => "parity-no-row",
    an_empty_table_has_no_row => "parity-empty",
    a_table_without_an_assembly_is_refused => "parity-without-assembly",
    an_assembly_without_a_table_is_refused => "assembly-without-parity",
    an_empty_assembly_counts_as_none => "empty-assembly",
    // The sweep id.
    a_sweep_id_of_another_shape_is_refused => "bad-sweep-id",
    an_empty_sweep_id_counts_as_none => "empty-sweep-id",
    a_sweep_id_without_a_table_is_checked_and_not_written => "sweep-id-alone",
    // Lines that are no digest lines.
    a_line_of_five_columns_is_quoted_and_refused => "five-columns",
    a_line_of_seven_columns_is_quoted_and_refused => "seven-columns",
    an_empty_line_is_quoted_and_refused => "blank-line",
    the_first_file_is_the_one_reported => "both-malformed",
    // The command line.
    options_that_are_required_are_named => "missing-required",
    an_option_that_is_unknown_is_named => "unknown-option",
    a_prefix_of_two_options_is_named => "ambiguous-option",
    an_option_without_its_value_is_named => "no-value",
    // Refusals in the tool's own words.
    a_digest_file_that_is_absent => "tb-missing-dump",
    a_digest_file_that_is_not_utf8 => "tb-not-utf8",
    an_output_directory_that_is_a_file => "tb-out-dir-is-a-file",
    a_table_without_an_assembly_column => "tb-parity-no-assembly",
    a_table_that_is_absent => "tb-parity-missing",
}

#[test]
fn help_goes_to_standard_output() {
    let scratch = Scratch::new("help");
    for arguments in [&["-h"][..], &["--help"], &["--version-a", "1", "--he"]] {
        let run = run_in(scratch.path(), arguments);
        assert_eq!(run.status.code(), Some(0), "{arguments:?}");
        assert!(lines(&run.stdout).starts_with("usage: vep-storable-digest-join [-h] --version-a"));
        assert!(run.stderr.is_empty());
    }
}
