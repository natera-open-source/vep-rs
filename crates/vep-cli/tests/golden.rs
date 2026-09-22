// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Golden-corpus tests: every consequence-set combination observed on the
//! datasets listed in `manuscript/data/dataset_inventory.csv`, per Ensembl
//! release and assembly, annotated with the corpus's own pruned cache and
//! compared column by column with the Ensembl VEP output committed beside it.
//!
//! One run per corpus; every mismatch is collected and reported in one failure
//! message, grouped by column and consequence set, so a change is measured
//! against the whole corpus rather than its first failing row.

mod common;

use std::collections::BTreeSet;

use common::{
    compare_entries, compare_tab_headers, corpora, parse_default, read_expected, run_vep,
    Documented,
};

#[test]
fn default_output_matches_vep_on_every_corpus() {
    let corpora = corpora();
    assert!(!corpora.is_empty(), "no corpus found under tests/golden");
    let tmp = tempfile::tempdir().unwrap();
    let mut failures = Vec::new();
    for corpus in &corpora {
        let expected = read_expected(&corpus.dir, "default.txt");
        let actual = run_vep(corpus, "default", tmp.path(), &[]);
        let documented = Documented::from_manifest(&corpus.manifest);
        let header = compare_tab_headers(&expected, &actual);
        let body = compare_entries(
            &parse_default(&expected),
            &parse_default(&actual),
            &documented,
        );
        if !header.is_empty() || !body.is_empty() {
            failures.push(format!(
                "== corpus {}/{} (default format)\n{header}{body}",
                corpus.release, corpus.assembly
            ));
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

/// The manifest is internally consistent with the committed inputs and
/// expected output: every record it lists is in `variants.vcf`, every
/// combination it claims to cover appears in VEP's output for the corpus, and
/// the generator left no combination without an exemplar.
#[test]
fn manifest_covers_every_combination_it_claims() {
    for corpus in corpora() {
        let m = &corpus.manifest;
        let records = m["records"].as_array().unwrap();
        let inputs = common::input_records(&corpus);
        assert_eq!(
            records.len(),
            inputs.len(),
            "{}/{}: manifest records vs variants.vcf data lines",
            corpus.release,
            corpus.assembly
        );
        for (i, (rec, line)) in records.iter().zip(inputs.iter()).enumerate() {
            assert_eq!(
                rec["vcf"].as_str().unwrap(),
                line,
                "{}/{}: record {i} differs between manifest and variants.vcf",
                corpus.release,
                corpus.assembly
            );
        }
        assert!(
            m["uncovered_combinations"].as_array().unwrap().is_empty(),
            "{}/{}: combinations without an exemplar: {:?}",
            corpus.release,
            corpus.assembly,
            m["uncovered_combinations"]
        );
        let expected = read_expected(&corpus.dir, "default.txt");
        let observed: BTreeSet<String> = parse_default(&expected)
            .iter()
            .map(|e| common::consequence_set(&e.consequence))
            .collect();
        let claimed: BTreeSet<String> = m["combinations"]
            .as_object()
            .unwrap()
            .keys()
            .map(|k| common::consequence_set(k))
            .collect();
        let unmet: Vec<&String> = claimed.difference(&observed).collect();
        assert!(
            unmet.is_empty(),
            "{}/{}: {} combination(s) in the manifest never appear in expected/default.txt: {:?}",
            corpus.release,
            corpus.assembly,
            unmet.len(),
            unmet
        );
        let listed: BTreeSet<String> = records
            .iter()
            .flat_map(|r| r["combinations"].as_array().unwrap().iter())
            .map(|c| common::consequence_set(c.as_str().unwrap()))
            .collect();
        assert_eq!(
            listed, claimed,
            "{}/{}: the records' combinations and the combination table disagree",
            corpus.release, corpus.assembly
        );
    }
}

/// Each corpus stays inside its in-repository size budget.
#[test]
fn each_corpus_fits_its_size_budget() {
    const BUDGET_BYTES: u64 = 15 * 1024 * 1024;
    for corpus in corpora() {
        let size = common::dir_size(&corpus.dir);
        assert!(
            size <= BUDGET_BYTES,
            "{}/{} is {} bytes on disk, over the {} byte budget",
            corpus.release,
            corpus.assembly,
            size,
            BUDGET_BYTES
        );
    }
}
