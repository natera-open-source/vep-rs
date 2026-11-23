// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `covers` on a copy of the two-record corpus: a row of the matrix is added to the
//! manifest's list or to a flag set's, a name that is no row is refused before anything is
//! written, and the list is printed one name a line.

mod common;

use common::{
    clean, failed, manifest, manifest_text, matrix, names, run, scratch, two_record_corpus,
};
use vep_golden_corpus::json::Value;

#[test]
fn covers_subcommand_adds_a_row_and_refuses_a_non_row() {
    let scratch = scratch("covers");
    let out = two_record_corpus(&scratch.path);
    let matrix = matrix(
        &scratch.path,
        &["symbol", "max_af", "flag_matrix_concordance"],
    );
    let before = manifest_text(&out);
    assert_eq!(
        run(&[
            "covers",
            &out,
            "--add",
            "no_such_row",
            "--add",
            "max_af",
            "--matrix",
            &matrix
        ]),
        failed(&format!(
            "covers names no_such_row, which is not a row of {matrix}"
        ))
    );
    assert_eq!(manifest_text(&out), before, "a refused run writes nothing");
    assert_eq!(
        run(&[
            "covers",
            &out,
            "--add",
            "flag_matrix_concordance",
            "--add",
            "symbol",
            "--matrix",
            &matrix
        ]),
        clean(0, "")
    );
    assert_eq!(
        run(&["covers", &out, "--list", "--matrix", &matrix]),
        clean(0, "symbol\nflag_matrix_concordance\n")
    );
    let after = manifest_text(&out);
    // Only the list changed: the manifest round-trips byte for byte around it.
    assert_eq!(
        after.replace(
            "\"covers\": [\n  \"symbol\",\n  \"flag_matrix_concordance\"\n ]",
            "\"covers\": [\n  \"symbol\"\n ]"
        ),
        before
    );
    assert_eq!(
        run(&[
            "flag-set",
            &out,
            "--name",
            "af",
            "--flag=--max_af",
            "--matrix",
            &matrix
        ])
        .code,
        0
    );
    let before = manifest_text(&out);
    assert_eq!(
        run(&[
            "covers",
            &out,
            "--flag-set",
            "af",
            "--add",
            "ghost",
            "--add",
            "max_af",
            "--matrix",
            &matrix
        ]),
        failed(&format!(
            "flag set af covers names ghost, which is not a row of {matrix}"
        ))
    );
    assert_eq!(manifest_text(&out), before);
    assert_eq!(
        run(&[
            "covers",
            &out,
            "--flag-set",
            "af",
            "--add",
            "max_af",
            "--matrix",
            &matrix
        ]),
        clean(0, "")
    );
    let written = manifest(&out);
    let sets = written.get("flag_sets").and_then(Value::as_list).unwrap();
    assert_eq!(names(sets[0].as_object().unwrap(), "covers"), ["max_af"]);
    assert_eq!(
        names(&written, "covers"),
        ["symbol", "flag_matrix_concordance"]
    );
    assert_eq!(
        run(&[
            "covers",
            &out,
            "--flag-set",
            "none",
            "--list",
            "--matrix",
            &matrix
        ]),
        failed("the manifest has no flag set 'none'")
    );
}

#[test]
fn covers_lists_with_an_addition_only_on_request_and_keeps_each_name_once() {
    let scratch = scratch("covers-list");
    let out = two_record_corpus(&scratch.path);
    let matrix = matrix(&scratch.path, &["symbol", "max_af", "pick"]);
    assert_eq!(run(&["covers", &out]), clean(0, "symbol\n"));
    assert_eq!(
        run(&[
            "covers", &out, "--add", "max_af", "--add", "max_af", "--add", "symbol", "--list",
            "--matrix", &matrix
        ]),
        clean(0, "symbol\nmax_af\n")
    );
    // The matrix is read only where a name is added.
    assert_eq!(
        run(&["covers", &out, "--list", "--matrix", "no/such/matrix.tsv"]),
        clean(0, "symbol\nmax_af\n")
    );
    // An empty name is refused under the name of the option `select` and `flag-set` take.
    assert_eq!(
        run(&["covers", &out, "--add", "pick", "--add", " ", "--matrix", &matrix]),
        failed("--covers wants a non-empty row name")
    );
    assert_eq!(names(&manifest(&out), "covers"), ["symbol", "max_af"]);
}
