// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use super::*;

fn texts(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| value.to_string()).collect()
}

fn outcome(verdict: &'static str, ensembl_vep: &[&str], vep_rs: &[&str], detail: &str) -> Outcome {
    Outcome {
        entry: "x".to_string(),
        verdict,
        ensembl_vep: Some(texts(ensembl_vep)),
        vep_rs: Some(texts(vep_rs)),
        detail: detail.to_string(),
        implementation: match verdict {
            "still-present" => Vec::new(),
            _ => vec![("code", texts(&["a.rs::f"]))],
        },
    }
}

#[test]
fn one_runs_outcome_is_the_entrys() {
    let own = outcome("closed", &["a"], &["a"], "Ensembl VEP ['a'], vep-rs ['a']");
    assert_eq!(entry_outcome(vec![own.clone()]), own);
}

#[test]
fn runs_that_read_alike_give_the_entry_their_verdict_and_their_values_in_order() {
    let present = entry_outcome(vec![
        outcome("still-present", &["a"], &[], ""),
        outcome("still-present", &["b", "c"], &["d"], ""),
    ]);
    assert_eq!(present.entry, "x");
    assert_eq!(present.verdict, "still-present");
    assert_eq!(present.ensembl_vep, Some(texts(&["a", "b", "c"])));
    assert_eq!(present.vep_rs, Some(texts(&["d"])));
    assert_eq!(present.detail, "");
    assert!(present.implementation.is_empty());

    let closed = entry_outcome(vec![
        outcome("closed", &["a"], &["a"], "one"),
        outcome("closed", &["b"], &["b"], "two"),
    ]);
    assert_eq!(closed.verdict, "closed");
    assert_eq!(closed.detail, "own run: one; also[1]: two");
    assert_eq!(closed.implementation, [("code", texts(&["a.rs::f"]))]);
}

#[test]
fn runs_that_read_differently_give_changed_and_name_the_run_that_moved() {
    let mixed = entry_outcome(vec![
        outcome("still-present", &["a"], &[], ""),
        outcome("still-present", &["b"], &[], ""),
        outcome("closed", &["c"], &["c"], "Ensembl VEP ['c'], vep-rs ['c']"),
    ]);
    assert_eq!(mixed.verdict, "changed");
    assert_eq!(mixed.detail, "also[2]: Ensembl VEP ['c'], vep-rs ['c']");
    assert_eq!(mixed.ensembl_vep, Some(texts(&["a", "b", "c"])));
    assert_eq!(mixed.vep_rs, Some(texts(&["c"])));
    assert_eq!(mixed.implementation, [("code", texts(&["a.rs::f"]))]);
}

#[test]
fn a_run_that_could_not_be_made_makes_the_entry_an_error() {
    let failed = entry_outcome(vec![
        outcome("still-present", &["a"], &[], ""),
        Outcome::error("x", "input /nowhere does not exist".to_string()),
    ]);
    assert_eq!(failed.verdict, "error");
    assert_eq!(failed.detail, "also[1]: input /nowhere does not exist");
    assert_eq!(failed.ensembl_vep, None);
    assert_eq!(failed.vep_rs, None);
    assert!(failed.implementation.is_empty());
}
