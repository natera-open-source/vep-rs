// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Every long flag the `vep` binary accepts has a row in `tests/parity/matrix.tsv`.
//!
//! clap is the arbiter of the flag set: hidden flags, the flags of flattened structs and the
//! generated `--help` and `--version` all come back from `Command::get_arguments`, where
//! `scripts/parity/check_parity_matrix.py` can only read the attributes of `args.rs`. A flag
//! `--foo-bar` or `--foo_bar` is the row `foo_bar`.

use std::collections::BTreeSet;

use clap::CommandFactory;
use vep_cli::args::Args;

const MATRIX: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/parity/matrix.tsv"
));

fn matrix_items() -> BTreeSet<&'static str> {
    let mut lines = MATRIX.lines();
    let header = lines
        .next()
        .expect("tests/parity/matrix.tsv has a header row");
    assert_eq!(
        header.split('\t').next(),
        Some("item"),
        "the first column of tests/parity/matrix.tsv is item"
    );
    lines
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').next().unwrap_or_default())
        .collect()
}

#[test]
fn every_long_flag_has_a_matrix_row() {
    let items = matrix_items();
    let mut command = Args::command();
    command.build();
    assert!(
        command.get_subcommands().next().is_none(),
        "subcommands are not rows of the matrix"
    );
    let longs: Vec<&str> = command
        .get_arguments()
        .filter_map(|arg| arg.get_long())
        .collect();
    assert!(
        longs.contains(&"input_file"),
        "clap reports no --input_file among {} long flags",
        longs.len()
    );
    let missing: Vec<String> = longs
        .iter()
        .map(|long| long.replace('-', "_"))
        .filter(|row| !items.contains(row.as_str()))
        .collect();
    assert!(
        missing.is_empty(),
        "flags without a row in tests/parity/matrix.tsv: {missing:?}"
    );
}
