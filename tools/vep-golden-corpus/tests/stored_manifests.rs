// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The committed corpora under `tests/golden/`: every manifest is read and written back to
//! its bytes, passes the shape check, records no unexplained residual, and names only mask
//! classes the golden harness and this tool both hold.

mod common;

use common::{committed_manifests, corpus_label, names, object, repository, run, RECODER_CLASSES};
use regex::Regex;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs;
use vep_golden_corpus::constants::{
    mask_classes, CORRECTED_BEHAVIOUR_CLASSES, HARNESS_ONLY_MASK_CLASSES, ORDER_ONLY,
    REGULATORY_FEATURE_ROWS, UNEXPLAINED,
};
use vep_golden_corpus::json::{Object, Value};
use vep_golden_corpus::manifest::{manifest_shape_error, manifest_text, masks_shape_error};

fn harness() -> String {
    fs::read_to_string(repository().join("crates/vep-cli/tests/common/mod.rs")).unwrap()
}

/// The lines of the list a harness constant opens with `opening`, without the opening line.
fn harness_block(harness: &str, opening: &str) -> Vec<String> {
    let start = harness.find(opening).expect("the constant");
    let end = start + harness[start..].find("];").expect("the list's end");
    let lines = harness[start..end].lines().skip(1);
    lines.map(|line| line.trim().to_string()).collect()
}

/// Every run of a manifest under the name a message gives it, the manifest's own first and
/// then each flag set's.
fn runs(manifest: &Object) -> Vec<(String, &Object)> {
    let mut runs = vec![("the manifest's run".to_string(), manifest)];
    let sets = manifest.get("flag_sets").and_then(Value::as_list);
    for set in sets.unwrap_or_default() {
        let set = set.as_object().expect("a flag set");
        let name = set.get("name").and_then(Value::as_str).unwrap_or("None");
        runs.push((format!("flag set {name}"), set));
    }
    runs
}

/// The places of a manifest, its own run and each flag set's, that record a row or a field
/// under the unexplained class: an entry of a classify list or a key of a summary.
fn unexplained_entries(manifest: &Object) -> Vec<String> {
    let mut found = Vec::new();
    for (owner, run) in runs(manifest) {
        for key in ["divergences", "vep_rs_only_tuples", "field_divergences"] {
            let entries = run.get(key).and_then(Value::as_list).unwrap_or_default();
            let count = entries
                .iter()
                .filter(|entry| {
                    entry.get("expected_divergence").and_then(Value::as_str) == Some(UNEXPLAINED)
                })
                .count();
            if count > 0 {
                found.push(format!("{owner}: {count} in {key}"));
            }
        }
        for key in [
            "divergence_summary",
            "vep_rs_only_summary",
            "field_divergence_summary",
        ] {
            if run
                .get(key)
                .is_some_and(|summary| summary.get(UNEXPLAINED).is_some())
            {
                found.push(format!("{owner}: {key}"));
            }
        }
    }
    found
}

#[test]
fn every_committed_manifest_is_read_and_written_back_to_its_bytes() {
    let mut written_back = 0;
    for path in committed_manifests() {
        let stored = fs::read_to_string(&path).unwrap();
        let manifest = object(&stored);
        let written = manifest_text(&manifest);
        if manifest.contains_key("tool_cases") {
            // A tool corpus is laid out by hand; what the writer makes of it is its own
            // fixed point.
            assert_eq!(manifest_text(&object(&written)), written, "{path:?}");
        } else {
            assert_eq!(written, stored, "{}", corpus_label(&path));
            written_back += 1;
        }
    }
    assert!(written_back >= 109, "{written_back} corpora written back");
}

#[test]
fn every_committed_manifest_passes_the_builder_shape_check() {
    let mut refused = Vec::new();
    for path in committed_manifests() {
        let mut manifest = object(&fs::read_to_string(&path).unwrap());
        let corpus = corpus_label(&path);
        let variants = path.parent().unwrap().join("variants.vcf");
        if manifest.contains_key("tool_cases") {
            assert!(
                !variants.exists(),
                "{corpus}: a tool corpus holds no variants.vcf"
            );
            let unknown: BTreeSet<String> = names(&manifest, "masks")
                .into_iter()
                .filter(|class| !RECODER_CLASSES.contains(&class.as_str()))
                .collect();
            if !unknown.is_empty() {
                refused.push(format!(
                    "{corpus}: masks names classes the recoder comparison lacks {unknown:?}"
                ));
            }
            manifest.remove("masks");
        } else {
            assert!(
                variants.is_file(),
                "{corpus}: a corpus the harness reads holds variants.vcf"
            );
        }
        if let Some(reason) = manifest_shape_error(&manifest) {
            refused.push(format!("{corpus}: {reason}"));
        }
    }
    assert!(
        refused.is_empty(),
        "committed manifests the tool refuses:\n{}",
        refused.join("\n")
    );
}

#[test]
fn no_committed_manifest_carries_an_unexplained_residual() {
    let explained = object(r#"{"divergences": [{"expected_divergence": "splice_family_swap"}]}"#);
    assert_eq!(unexplained_entries(&explained), Vec::<String>::new());
    let entry = format!(r#"{{"location": "21:100", "expected_divergence": "{UNEXPLAINED}"}}"#);
    let unexplained = object(&format!(
        r#"{{"divergences": [{entry}], "divergence_summary": {{"{UNEXPLAINED}": 1}},
            "flag_sets": [{{"name": "hgvs", "vep_rs_only_tuples": [{entry}],
                "field_divergences": [{entry}, {entry}], "vep_rs_only_summary": {{"{UNEXPLAINED}": 1}},
                "field_divergence_summary": {{"{UNEXPLAINED}": 2}}}}]}}"#
    ));
    assert_eq!(
        unexplained_entries(&unexplained),
        [
            "the manifest's run: 1 in divergences",
            "the manifest's run: divergence_summary",
            "flag set hgvs: 1 in vep_rs_only_tuples",
            "flag set hgvs: 2 in field_divergences",
            "flag set hgvs: vep_rs_only_summary",
            "flag set hgvs: field_divergence_summary",
        ]
    );
    let mut found = Vec::new();
    for path in committed_manifests() {
        let manifest = object(&fs::read_to_string(&path).unwrap());
        for place in unexplained_entries(&manifest) {
            found.push(format!("{}: {place}", corpus_label(&path)));
        }
    }
    assert!(
        found.is_empty(),
        "committed manifests that carry {UNEXPLAINED}:\n{}",
        found.join("\n")
    );
}

#[test]
fn corrected_behaviour_classes_agree_across_the_tool_and_the_harness() {
    let corrected: HashSet<&str> = CORRECTED_BEHAVIOUR_CLASSES.into_iter().collect();
    assert_eq!(corrected.len(), CORRECTED_BEHAVIOUR_CLASSES.len());
    let block = harness_block(
        &harness(),
        "const CORRECTED_BEHAVIOUR_CLASSES: &[&str] = &[",
    );
    let harness_names: HashSet<String> = block
        .iter()
        .filter(|line| line.starts_with('"'))
        .map(|line| line.trim_matches(',').trim_matches('"').to_string())
        .collect();
    assert!(!harness_names.is_empty());
    let outside: Vec<&String> = harness_names
        .iter()
        .filter(|name| !corrected.contains(name.as_str()))
        .collect();
    assert!(
        outside.is_empty(),
        "harness classes the tool lacks: {outside:?}"
    );
    // A class the harness lacks is one no committed corpus names in its `masks`, the one
    // place the harness refuses an unknown name.
    let mut named_in_manifests = HashSet::new();
    for path in committed_manifests() {
        let manifest = object(&fs::read_to_string(&path).unwrap());
        for (_, run) in runs(&manifest) {
            named_in_manifests.extend(names(run, "masks"));
        }
    }
    let unheld: Vec<&&str> = corrected
        .iter()
        .filter(|class| !harness_names.contains(**class) && named_in_manifests.contains(**class))
        .collect();
    assert!(
        unheld.is_empty(),
        "classes a corpus names that the harness lacks: {unheld:?}"
    );
    assert_eq!(
        masks_shape_error(&Value::strings(&CORRECTED_BEHAVIOUR_CLASSES), "masks"),
        None
    );
}

#[test]
fn harness_mask_classes_are_the_tool_mask_classes() {
    let harness = harness();
    let declared = Regex::new(r#"pub const ([A-Z_]+): &str = "([a-z0-9_]+)";"#).unwrap();
    let constants: HashMap<&str, &str> = declared
        .captures_iter(&harness)
        .map(|found| {
            let (_, [name, value]) = found.extract();
            (name, value)
        })
        .collect();
    let block = harness_block(&harness, "const MASK_CLASSES: &[&str] = &[");
    let listed: Vec<&str> = block
        .iter()
        .filter(|line| !line.is_empty())
        .map(|line| line.trim_end_matches(','))
        .collect();
    let unresolved: Vec<&&str> = listed
        .iter()
        .filter(|name| !constants.contains_key(**name))
        .collect();
    assert!(
        unresolved.is_empty(),
        "the harness's MASK_CLASSES names constants without a string: {unresolved:?}"
    );
    let harness_classes: Vec<&str> = listed.iter().map(|name| constants[name]).collect();
    let mut tool_classes = vec![REGULATORY_FEATURE_ROWS, ORDER_ONLY];
    tool_classes.extend(HARNESS_ONLY_MASK_CLASSES);
    let all = mask_classes();
    assert!(
        harness_classes.iter().all(|class| all.contains(class)),
        "harness mask classes the tool lacks"
    );
    let as_set = |classes: &[&'static str]| classes.iter().copied().collect::<HashSet<&str>>();
    assert_eq!(
        harness_classes.iter().copied().collect::<HashSet<&str>>(),
        as_set(&tool_classes),
        "harness and tool harness-only mask classes differ"
    );
    assert_eq!(
        harness_classes.len(),
        harness_classes.iter().collect::<HashSet<_>>().len()
    );
    assert_eq!(tool_classes.len(), as_set(&tool_classes).len());
}

#[test]
fn covers_and_flag_set_list_every_committed_corpus() {
    for path in committed_manifests() {
        let manifest = object(&fs::read_to_string(&path).unwrap());
        if manifest.contains_key("masks") && manifest.contains_key("tool_cases") {
            continue;
        }
        let corpus = path.parent().unwrap().to_str().unwrap();
        let listed = run(&["covers", corpus, "--list"]);
        let covers: String = names(&manifest, "covers")
            .iter()
            .map(|name| format!("{name}\n"))
            .collect();
        assert_eq!(
            (listed.code, listed.stdout.as_str(), listed.stderr.as_str()),
            (0, covers.as_str(), ""),
            "{corpus}"
        );
        let sets = manifest.get("flag_sets").and_then(Value::as_list);
        let listed = run(&["flag-set", corpus, "--list"]);
        assert_eq!((listed.code, listed.stderr.as_str()), (0, ""), "{corpus}");
        assert_eq!(
            listed.stdout.lines().count(),
            sets.unwrap_or_default().len(),
            "{corpus}"
        );
        for (line, set) in listed.stdout.lines().zip(sets.unwrap_or_default()) {
            let name = set.get("name").and_then(Value::as_str).unwrap();
            let columns: Vec<&str> = line.split('\t').collect();
            assert_eq!(columns.len(), 8, "{corpus} {name}");
            assert_eq!(columns[0], name, "{corpus}");
            assert_eq!(
                columns[7],
                names(set.as_object().unwrap(), "covers").join(","),
                "{corpus} {name}"
            );
        }
    }
}
