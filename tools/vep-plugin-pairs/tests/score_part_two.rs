// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `score` program over small outputs under a planted registry and planted classes:
//! the value rule in whole runs, a tuple-file engine, the allele projection with the class
//! an allele tuple rests on, and the sources a run names.

mod support;

use serde_json::{json, Value};
use std::collections::BTreeSet;
use support::{
    allele, counts, entry, feature, feature_residuals, raw_counts, reference_rows, refusal, row,
    row_with, sha256, Planted, PlantedRegistry, PlantedRules, Run, Scratch, ERROR, PLANTED,
};
use vep_plugin_pairs::score::api::Sources;
use vep_plugin_pairs::score::number::numbers_agree;

const SMALL_NUMBERS_APART: [(&str, &str); 7] = [
    ("3E-24", "6E-9"),
    ("7E-9", "8E-10"),
    ("3E-11", "4E-18"),
    ("3E-18", "3E-8"),
    ("4E-7", "7E-7"),
    ("2E-16", "4E-8"),
    ("3E-9", "9E-7"),
];
const WHOLE_NUMBERS_ONE_APART: [(&str, &str); 5] = [
    ("1", "2"),
    ("0", "1"),
    ("1000000", "1000001"),
    ("250", "251"),
    ("-7", "-8"),
];
const SPELLINGS_OF_ONE_NUMBER: [(&str, &str); 6] = [
    ("5E-9", "5e-09"),
    ("5E-9", "0.000000005"),
    ("4E-24", "4.0E-24"),
    ("1.0", "1"),
    ("100", "1e2"),
    ("0.10", "0.1"),
];

/// The Planted entry and the residual lines of a run whose rows each write one pair of
/// values, the reference's and the engine's, under one key.
fn scored_pairs(tag: &str, pairs: &[(&str, &str)]) -> (Value, Vec<String>) {
    let scratch = Scratch::new(tag);
    let planted = Planted::new();
    let rows = |values: Vec<&str>| -> Vec<String> {
        values
            .iter()
            .enumerate()
            .map(|(at, value)| {
                row(
                    &format!("1:{}", at + 100),
                    "T",
                    "ENST1",
                    &format!("PLANT_score={value}"),
                )
            })
            .collect()
    };
    let reference = scratch.output(
        "reference.txt",
        &rows(pairs.iter().map(|pair| pair.0).collect()),
    );
    let engine = scratch.output(
        "engine.txt",
        &rows(pairs.iter().map(|pair| pair.1).collect()),
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().residuals(&residuals).no_table(),
    );
    assert_eq!(
        report["value_rule"],
        json!({
            "numeric_tolerance": "0.000001",
            "numeric_tolerance_applies": [
                "to the difference as a share of the larger number",
                "to the difference",
            ],
        })
    );
    (
        feature(&report, "Planted").clone(),
        feature_residuals(&residuals, "Planted"),
    )
}

#[test]
fn numbers_close_in_size_alone_or_one_apart_are_not_equal() {
    for (what, pairs) in [
        ("two small numbers", &SMALL_NUMBERS_APART[..]),
        ("whole numbers one apart", &WHOLE_NUMBERS_ONE_APART[..]),
    ] {
        for (one, other) in pairs {
            assert!(
                !numbers_agree(one.as_bytes(), other.as_bytes()),
                "{one} {other}"
            );
        }
        let (item, lines) = scored_pairs("apart", pairs);
        let held = pairs.len() as i64;
        assert_eq!(counts(&item), [held, held, 0, held, held, 0], "{what}");
        assert_eq!(item["text_only_differences"], 0, "{what}");
        assert_eq!(item["residual_rows"], 2 * held, "{what}");
        assert_eq!(lines.len(), 2 * pairs.len(), "{what}");
    }
}

#[test]
fn two_spellings_of_one_number_are_a_text_only_pair() {
    for (one, other) in SPELLINGS_OF_ONE_NUMBER {
        assert!(
            numbers_agree(one.as_bytes(), other.as_bytes()),
            "{one} {other}"
        );
    }
    let (item, lines) = scored_pairs("spellings", &SPELLINGS_OF_ONE_NUMBER);
    let pairs = SPELLINGS_OF_ONE_NUMBER.len() as i64;
    assert_eq!(counts(&item), [pairs; 6]);
    assert_eq!(item["text_only_differences"], pairs);
    assert_eq!(item["residual_rows"], 0);
    assert!(lines.is_empty(), "{lines:?}");
}

#[test]
fn value_rule_at_its_boundary_and_text_only_differences() {
    let scratch = Scratch::new("boundary");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=1.5"),
            row("1:200", "T", "ENST1", "PLANT_score=1.5"),
            row("1:300", "T", "ENST1", "PLANT_score=0.50"),
            row("1:400", "T", "ENST1", "PLANT_score=0.5;PLANT_label=A"),
            row("1:500", "T", "ENST1", "PLANT_score=0.5"),
            row("1:600", "T", "ENST1", "PLANT_score=0.5"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=1.500001"),
            row("1:200", "T", "ENST1", "PLANT_score=1.5000011"),
            row("1:300", "T", "ENST1", "PLANT_score=0.5"),
            row("1:400", "T", "ENST1", "PLANT_score=0.5;PLANT_label=a"),
            row("1:500", "T", "ENST1", "PLANT_score=0.5000005"),
            row("1:600", "T", "ENST1", "PLANT_score=0.500001"),
        ],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine), ("same", &reference)],
        &Run::new(),
    );
    let item = feature(&report, "Planted");
    assert_eq!(counts(item), [7, 7, 4, 7, 7, 4]);
    assert_eq!(item["text_only_differences"], 3);
    assert_eq!(
        entry(&report, "Planted", "same", "feature")["text_only_differences"],
        0
    );
    let item = allele(&report, "Planted");
    assert_eq!(counts(item), [7, 7, 4, 7, 7, 4]);
    assert_eq!(item["text_only_differences"], 3);
}

#[test]
fn values_under_one_key_pair_one_to_one() {
    let scratch = Scratch::new("one-to-one");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=0.1234561"),
            row("1:100", "T", "ENST1", "PLANT_score=0.1234562"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[row("1:100", "T", "ENST1", "PLANT_score=0.12345615")],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    let item = feature(&report, "Planted");
    assert_eq!(counts(item), [2, 1, 1, 2, 1, 1]);
    assert_eq!(item["text_only_differences"], 1);
}

#[test]
fn tuple_file_engine() {
    let scratch = Scratch::new("tuple-file");
    let planted = Planted::new();
    let reference = scratch.output("reference.txt", &reference_rows());
    let engine = scratch.tuples(
        "engine.tsv",
        &[
            ["1:100", "T", "ENST1", "PLANT_score", "0.5"],
            ["1:100", "T", "ENST1", "PLANT_label", "benign"],
            ["1:100", "T", "ENST2", "PLANT_score", "0.50"],
            ["1:200-201", "-", "ENST1", "PLANT_score", "0.25"],
            ["1:200-201", "-", "ENST1", "PLANT_label", "-"],
            ["1:200-201", "-", "ENST1", "IMPACT", "MODIFIER"],
            ["1:200-201", "-", "ENST1", "PlantedDynamic_a", "x"],
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    let item = feature(&report, "Planted");
    assert_eq!(counts(item), [4; 6]);
    assert_eq!(item["text_only_differences"], 1);
    assert_eq!(
        counts(feature(&report, "PlantedDynamic")),
        [2, 1, 1, 2, 1, 1]
    );
    let described = &report["engines"][0];
    assert_eq!(described["format"], "tuples");
    assert_eq!(described["tuples"], 7);
    assert_eq!(described["tuples_outside_every_source"], 1);
    assert_eq!(
        described["sha256"],
        json!(sha256(&std::fs::read(&engine).unwrap()))
    );
}

#[test]
fn allele_projection_drops_the_feature_and_deduplicates() {
    let scratch = Scratch::new("allele-projection");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=0.9"),
            row("1:100", "T", "ENST2", "PLANT_score=0.9"),
            row("1:100", "T", "ENST3", "IMPACT=LOW"),
            row("1:100", "G", "ENST1", "PLANT_score=0.2"),
        ],
    );
    let engine = scratch.tuples(
        "engine.tsv",
        &[
            ["1:100", "T", "ENST1", "PLANT_score", "0.9"],
            ["1:100", "T", "ENST3", "PLANT_score", "0.9"],
            ["1:100", "T", "OWN_MODEL.1", "PLANT_score", "0.9"],
            ["1:100", "G", "OWN_MODEL.1", "PLANT_score", "0.2"],
            ["1:100", "C", "OWN_MODEL.1", "PLANT_score", "0.7"],
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    assert_eq!(counts(feature(&report, "Planted")), [3, 5, 1, 3, 5, 1]);
    assert_eq!(counts(allele(&report, "Planted")), [2, 3, 2, 2, 3, 2]);
}

#[test]
fn allele_tuple_is_set_aside_when_every_unpaired_feature_tuple_was() {
    let scratch = Scratch::new("allele-aside");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "IMPACT=LOW"),
            row("1:100", "T", "ENST2", "IMPACT=LOW"),
            row("1:200", "C", "ENST1", "IMPACT=LOW"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "LOEUF=0"),
            row("1:100", "T", "ENST2", "LOEUF=0"),
            row("1:200", "C", "ENST1", "LOEUF=0"),
            row("1:200", "C", "ENST7", "LOEUF=0"),
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    let item = feature(&report, "LOEUF");
    assert_eq!(counts(item), [0, 4, 0, 0, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_engine_zero": 3}));
    let item = allele(&report, "LOEUF");
    assert_eq!(counts(item), [0, 2, 0, 0, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_engine_zero": 1}));
    let mut every: Vec<&str> = PLANTED.iter().map(|(name, _)| *name).collect();
    every.sort_unstable();
    assert_eq!(report["not_populated"], json!(every));
}

/// A class table that gives the one field of the plugin GWAS two classes.
const TWO_CLASSES_OF_ONE_FIELD: &str =
    "planted_row_cdna_one\tGWAS\tGWAS_pmid\nplanted_reference_value\tGWAS\tGWAS_pmid\n";

#[test]
fn allele_tuple_under_an_applied_and_a_held_class_rests_on_the_held_one() {
    let scratch = Scratch::new("held-class");
    let planted = Planted::registry(&[("GWAS", &["GWAS_pmid"])]);
    let reference = scratch.output(
        "reference.txt",
        &[
            row_with("1:100", "T", "ENST1", "GWAS_pmid=17", "1", "v"),
            row("1:100", "T", "ENST2", "GWAS_pmid=17"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row_with("1:100", "T", "ENST1", "IMPACT=LOW", "1", "v"),
            row("1:100", "T", "ENST2", "IMPACT=LOW"),
        ],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new()
            .table(TWO_CLASSES_OF_ONE_FIELD)
            .applied(&["planted_row_cdna_one"]),
    );
    let by_feature = feature(&report, "GWAS");
    assert_eq!(raw_counts(by_feature), [2, 0, 0]);
    assert_eq!(by_feature["masked"], json!({"planted_row_cdna_one": 1}));
    assert_eq!(
        by_feature["not_applied"],
        json!({"planted_reference_value": 1})
    );
    assert_eq!(by_feature["residual_rows"], 1);
    let by_allele = allele(&report, "GWAS");
    assert_eq!(raw_counts(by_allele), [1, 0, 0]);
    assert_eq!(by_allele["masked"], json!({"planted_row_cdna_one": 0}));
    assert_eq!(
        by_allele["not_applied"],
        json!({"planted_reference_value": 1})
    );
    assert_eq!(by_allele["residual_rows"], 1);
    for item in [by_feature, by_allele] {
        let keys = item.as_object().unwrap().keys();
        assert!(
            keys.into_iter().all(|key| !key.starts_with("adj_")),
            "{item}"
        );
    }
}

#[test]
fn several_unpaired_values_of_one_allele_each_find_their_class() {
    let scratch = Scratch::new("several-values");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "MaxEntScan_diff=1.5"),
            row("1:100", "T", "ENST2", "MaxEntScan_diff=2.5"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "MaxEntScan_diff=-1.25"),
            row("1:100", "T", "ENST2", "MaxEntScan_diff=-2.25"),
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    for projection in ["feature", "allele"] {
        let item = entry(&report, "MaxEntScan", "engine", projection);
        assert_eq!(counts(item), [2, 2, 0, 0, 0, 0], "{projection}");
        assert_eq!(item["residual_rows"], 0, "{projection}");
    }
}

#[test]
fn allele_tuple_whose_feature_tuples_all_paired_stays_in_both_counts() {
    let scratch = Scratch::new("all-paired");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "LOEUF=0"),
            row("1:100", "T", "ENST2", "IMPACT=LOW"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "LOEUF=0.0"),
            row("1:100", "T", "ENST2", "LOEUF=0"),
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    let item = feature(&report, "LOEUF");
    assert_eq!(counts(item), [1, 2, 1, 1, 1, 1]);
    assert_eq!(item["masked"], json!({"planted_engine_zero": 1}));
    assert_eq!(item["text_only_differences"], 1);
    let item = allele(&report, "LOEUF");
    assert_eq!(counts(item), [1, 2, 1, 1, 2, 1]);
    assert_eq!(item["masked"], json!({"planted_engine_zero": 0}));
    assert_eq!(item["text_only_differences"], 0);
}

#[test]
fn allele_tuple_is_credited_to_the_first_class_in_table_order() {
    let scratch = Scratch::new("first-class");
    let planted = Planted::registry(&[("GWAS", &["GWAS_pmid"])]);
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "GWAS_pmid=17"),
            row_with("1:100", "T", "ENST2", "GWAS_pmid=17", "1", "v"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "IMPACT=LOW"),
            row_with("1:100", "T", "ENST2", "IMPACT=LOW", "1", "v"),
        ],
    );
    let run = Run::new().table(TWO_CLASSES_OF_ONE_FIELD);
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &run);
    assert_eq!(
        feature(&report, "GWAS")["masked"],
        json!({"planted_row_cdna_one": 1, "planted_reference_value": 1})
    );
    let item = allele(&report, "GWAS");
    assert_eq!(counts(item), [1, 0, 0, 0, 0, 0]);
    assert_eq!(
        item["masked"],
        json!({"planted_row_cdna_one": 1, "planted_reference_value": 0})
    );
    let held = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &run.applied(&[]),
    );
    assert_eq!(
        allele(&held, "GWAS")["not_applied"],
        json!({"planted_row_cdna_one": 1})
    );
}

#[test]
fn custom_source_groups_the_fields_under_its_prefix() {
    let scratch = Scratch::new("custom-prefix");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "ClinVar=7105;ClinVar_CLNSIG=Benign;ClinVarious=1;gnomAD_AF=0.01",
        )],
    );
    let engine = scratch.tuples(
        "engine.tsv",
        &[
            ["1:100", "T", "ENST1", "ClinVar", "7105"],
            ["1:100", "T", "ENST1", "ClinVar_CLNSIG", "Pathogenic"],
            ["1:100", "T", "ENST1", "gnomAD_AF", "0.0100"],
        ],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().custom(&[("ClinVar", "ClinVar"), ("gnomAD", "gnomAD_")]),
    );
    assert_eq!(counts(feature(&report, "ClinVar")), [3, 2, 1, 3, 2, 1]);
    assert_eq!(counts(feature(&report, "gnomAD")), [1; 6]);
    assert_eq!(feature(&report, "gnomAD")["text_only_differences"], 1);
    assert_eq!(
        report["sources"],
        json!({
            "ClinVar": ["ClinVar", "ClinVar_CLNSIG", "ClinVarious"],
            "gnomAD": ["gnomAD_AF"],
        })
    );
}

#[test]
fn custom_source_may_not_take_a_registry_plugins_name() {
    let registry = PlantedRegistry::planted();
    let custom = [("Planted".to_string(), "PLANT_".to_string())];
    let refused = Sources::new(&[], &custom, None, &registry)
        .err()
        .map(|error| error.0);
    assert_eq!(
        refused.as_deref(),
        Some("--custom-source Planted names a plugin of the registry")
    );
}

#[test]
fn non_plugin_field_set_of_the_registry_is_no_source() {
    let scratch = Scratch::new("non-plugin");
    let planted = Planted::registry(&[("HGVS", &["HGVSc"]), ("Planted", &["PLANT_score"])]);
    let reference = scratch.output(
        "reference.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "HGVSc=ENST1:c.1A>T;PLANT_score=1",
        )],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new().no_table(),
    );
    assert_eq!(report["plugins"], json!(["Planted"]));
    let sources: BTreeSet<&str> = report["scores"]
        .as_array()
        .unwrap()
        .iter()
        .map(|item| item["source"].as_str().unwrap())
        .collect();
    assert_eq!(sources, BTreeSet::from(["Planted"]));
    assert_eq!(report["not_populated"], json!([]));
    let refused = refusal(planted.try_score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new().no_table().plugins(&["HGVS"]),
    ));
    assert_eq!(refused, "--plugins names HGVS, no plugin of the registry");
}

#[test]
fn source_only_the_engine_writes_is_scored_and_listed_as_not_populated() {
    let scratch = Scratch::new("engine-only");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "T", "ENST1", "PLANT_score=1")],
    );
    let engine = scratch.output(
        "engine.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "PLANT_score=1;PlantedDynamic_a=x",
        )],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().custom(&[("Unused", "Unused_")]),
    );
    assert_eq!(
        counts(feature(&report, "PlantedDynamic")),
        [0, 1, 0, 0, 1, 0]
    );
    let mut unpopulated: Vec<&str> = PLANTED
        .iter()
        .map(|(name, _)| *name)
        .filter(|name| *name != "Planted")
        .chain(["Unused"])
        .collect();
    unpopulated.sort_unstable();
    assert_eq!(report["not_populated"], json!(unpopulated));
}

#[test]
fn plugin_is_named_by_its_registry_name_or_an_alias_the_registry_lists() {
    let scratch = Scratch::new("alias");
    let planted = Planted {
        registry: PlantedRegistry::planted().alias("Loaded", "Planted"),
        rules: PlantedRules::planted(),
    };
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "T", "ENST1", "PLANT_score=1")],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new().no_table().plugins(&["Loaded", "Planted"]),
    );
    assert_eq!(report["plugins"], json!(["Planted"]));
    assert_eq!(counts(feature(&report, "Planted")), [1; 6]);
}

#[test]
fn field_the_registry_gives_a_plugin_is_scored_though_named_like_a_core_key() {
    let scratch = Scratch::new("core-named");
    let mut rows = PLANTED.to_vec();
    rows[0] = ("Planted", &["PLANT_score", "PLANT_label", "IMPACT"]);
    let planted = Planted::registry(&rows);
    let reference = scratch.output(
        "reference.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "IMPACT=LOW;PLANT_score=1;LOEUF=0.1",
        )],
    );
    let with_the_plugin = planted.score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new().plugins(&["Planted", "LOEUF"]),
    );
    assert_eq!(
        with_the_plugin["sources"]["Planted"],
        json!(["IMPACT", "PLANT_score"])
    );
    assert_eq!(counts(feature(&with_the_plugin, "Planted")), [2; 6]);
    assert_eq!(with_the_plugin["reference"]["unclaimed_keys"], json!({}));
    let without_it = planted.score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new().plugins(&["LOEUF"]),
    );
    assert_eq!(without_it["sources"], json!({"LOEUF": ["LOEUF"]}));
    assert_eq!(
        without_it["reference"]["unclaimed_keys"],
        json!({"PLANT_score": 1})
    );
}

#[test]
fn run_command_is_handed_to_the_registry_which_can_refuse_the_run() {
    const REFUSAL: &str = "the run loads Planted and sets --planted: PLANT_label is both the plugin's column and Ensembl VEP's own key";
    let scratch = Scratch::new("run-command");
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "T", "ENST1", "PLANT_score=1;PLANT_label=x")],
    );
    let no_class = scratch.file("no-class.tsv", "");
    let report_path = scratch.path("report.json");
    let run = |planted: &Planted, name: &str, more: &[&str]| {
        let engine = format!("engine={reference}");
        let (work, residuals) = (scratch.path("work"), scratch.path(name));
        let mut arguments = vec![
            "--reference",
            &reference,
            "--engine",
            &engine,
            "--plugins",
            "Planted",
            "--class-table",
            &no_class,
            "--work-dir",
            &work,
            "--sort-memory",
            "16M",
            "--report",
            &report_path,
            "--residuals",
            &residuals,
        ];
        arguments.extend(more);
        planted.run_main(&arguments)
    };
    let stated =
        || -> Value { serde_json::from_slice(&std::fs::read(&report_path).unwrap()).unwrap() };
    let (code, _, error) = run(&Planted::new(), "unasked", &[]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert_eq!(stated()["run_command"], Value::Null);
    let refusing = || Planted {
        registry: PlantedRegistry::planted().resolving(|field, _, flags| {
            let planted = flags.is_some_and(|flags| flags.iter().any(|flag| flag == "--planted"));
            if planted && field == "PLANT_label" {
                return Err(REFUSAL.to_string());
            }
            Ok(field.starts_with("PLANT_").then(|| "Planted".to_string()))
        }),
        rules: PlantedRules::planted(),
    };
    let handed = refusing();
    let line = ["--run-command", "vep --plugin 'Planted,a file' --offline"];
    let (code, _, error) = run(&handed, "handed", &line);
    assert_eq!((code, error.as_str()), (0, ""));
    let words = support::strings(&["vep", "--plugin", "Planted,a file", "--offline"]);
    let asked: BTreeSet<Option<Vec<String>>> = handed
        .registry
        .asked()
        .into_iter()
        .map(|(_, _, flags)| flags)
        .collect();
    assert_eq!(asked, BTreeSet::from([Some(words.clone())]));
    assert_eq!(stated()["run_command"], json!(words));
    let (code, lines, error) = run(&refusing(), "refused", &["--run-command", "vep --planted"]);
    assert_eq!((code, lines.len()), (2, 0));
    assert_eq!(error, format!("{ERROR}{REFUSAL}\n"));
    assert!(scratch.exists("work") && scratch.names("work").is_empty());
    assert!(scratch.exists("refused") && scratch.names("refused").is_empty());
}

#[test]
fn run_without_a_plugin_scores_its_custom_sources_alone() {
    let scratch = Scratch::new("custom-alone");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "T", "ENST1", "PLANT_score=1;Custom=id1")],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new()
            .no_table()
            .plugins(&[])
            .custom(&[("Custom", "Custom")]),
    );
    assert_eq!(counts(feature(&report, "Custom")), [1; 6]);
    assert!(planted.registry.asked().is_empty());
}
