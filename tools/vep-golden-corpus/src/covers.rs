// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `covers` subcommand: lists the parity-matrix rows a corpus or one of its flag sets
//! covers, or adds rows of the matrix to that list and rewrites the manifest.

use crate::cli::CoversArgs;
use crate::error::{complain, say, Failure};
use crate::json::Value;
use crate::manifest::{
    covers_names, covers_owner, matrix_items, names_outside_matrix, names_under, read_manifest,
    target_run, target_run_mut, without_repeats, write_manifest,
};
use crate::pytext::path_str;

/// Runs `covers`. A name that is not a row of the matrix is refused with one line a name
/// and exit code 1, and nothing is written; the list is printed when `--list` is given or
/// nothing is added.
pub fn run(args: &CoversArgs) -> Result<i32, Failure> {
    let corpus = path_str(&args.corpus);
    let mut manifest = read_manifest(&corpus)?;
    let flag_set = args.flag_set.as_deref();
    let mut covers = names_under(target_run(&manifest, flag_set)?, "covers");
    if !args.add.is_empty() {
        let added = covers_names(&args.add)?;
        let matrix = path_str(&args.matrix);
        let items = matrix_items(&matrix)?;
        let errors = names_outside_matrix(&added, &items, &matrix, &covers_owner(flag_set));
        for error in &errors {
            complain(error);
        }
        if !errors.is_empty() {
            return Ok(1);
        }
        covers = without_repeats(covers.into_iter().chain(added));
        target_run_mut(&mut manifest, flag_set)?.insert("covers", Value::strings(&covers));
        write_manifest(&corpus, &manifest)?;
    }
    if args.list || args.add.is_empty() {
        for name in &covers {
            say(name);
        }
    }
    Ok(0)
}
