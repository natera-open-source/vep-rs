// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `score` program over small outputs under a planted registry and planted classes:
//! which source a field counts for, which tuples a documented class sets aside, how a
//! class table is read and refused, and what the residual files hold.

mod support;

use serde_json::{json, Map, Value};
use std::collections::BTreeSet;
use std::os::unix::fs::PermissionsExt;
use support::{
    allele, counts, entry, feature, feature_residuals, names, raw_counts, reference_rows, refusal,
    residual_lines, row, sha256, strings, Planted, PlantedRegistry, PlantedRules, Run, Scratch,
    COUNT_KEYS, HEADER,
};

#[test]
fn reference_against_itself_pairs_every_tuple_of_every_source() {
    let scratch = Scratch::new("itself");
    let planted = Planted::new();
    let reference = scratch.output("reference.txt", &reference_rows());
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new().custom(&[("Custom", "Custom")]),
    );
    let scores = report["scores"].as_array().unwrap();
    let sources: BTreeSet<&str> = scores
        .iter()
        .map(|item| item["source"].as_str().unwrap())
        .collect();
    assert_eq!(
        sources,
        BTreeSet::from([
            "Planted",
            "PlantedDynamic",
            "LOEUF",
            "MaxEntScan",
            "GWAS",
            "Custom"
        ])
    );
    for item in scores {
        let held = counts(item);
        assert!(held.iter().all(|count| *count == held[0]), "{item}");
        let masked = item["masked"].as_object().unwrap();
        assert!(masked.values().all(|keys| keys == 0), "{item}");
        assert_eq!(item["not_applied"], json!({}), "{item}");
        assert_eq!(item["text_only_differences"], 0, "{item}");
        assert!(COUNT_KEYS.iter().all(|key| item[*key].is_i64()), "{item}");
    }
    assert_eq!(counts(feature(&report, "Planted")), [4; 6]);
    assert_eq!(counts(allele(&report, "Planted")), [3; 6]);
    assert_eq!(counts(feature(&report, "PlantedDynamic")), [2; 6]);
    assert_eq!(counts(feature(&report, "Custom")), [2; 6]);
    assert_eq!(report["sources"]["Custom"], json!(["Custom", "Custom_SIG"]));
    assert_eq!(
        report["not_populated"],
        json!(["NearestExonJB", "REVEL", "UTRAnnotator"])
    );
    assert_eq!(
        report["tuple"]["allele"],
        json!(["Location", "Allele", "field", "value"])
    );
    assert!(planted.calls().is_empty());
}

#[test]
fn documented_class_difference_moves_the_raw_counts_alone() {
    let scratch = Scratch::new("documented");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let mut changed = rows.clone();
    changed[1] = rows[1].replace("MaxEntScan_diff=1.5", "MaxEntScan_diff=-2.25");
    let engine = scratch.output("engine.txt", &changed);
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    for projection in ["feature", "allele"] {
        let item = entry(&report, "MaxEntScan", "engine", projection);
        assert_eq!(counts(item), [1, 1, 0, 0, 0, 0], "{projection}");
        assert_eq!(
            item["masked"],
            json!({"planted_two_numbers": 1}),
            "{projection}"
        );
    }
    assert_eq!(counts(feature(&report, "Planted")), [4; 6]);
}

#[test]
fn documented_class_sets_aside_a_tuple_one_side_alone_writes() {
    let scratch = Scratch::new("one-side");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let mut filled = rows.clone();
    filled[2] = rows[2].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;LOEUF=0");
    let engine = scratch.output("engine.txt", &filled);
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    let item = feature(&report, "LOEUF");
    assert_eq!(counts(item), [1, 2, 1, 1, 1, 1]);
    assert_eq!(item["masked"], json!({"planted_engine_zero": 1}));
}

#[test]
fn undocumented_difference_moves_the_raw_and_the_adjusted_counts() {
    let scratch = Scratch::new("undocumented");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let mut changed = rows.clone();
    changed[0] = rows[0].replace("PLANT_label=benign", "PLANT_label=pathogenic");
    let engine = scratch.output("engine.txt", &changed);
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    for (projection, expected) in [
        ("feature", [4, 4, 3, 4, 4, 3]),
        ("allele", [3, 3, 2, 3, 3, 2]),
    ] {
        let item = entry(&report, "Planted", "engine", projection);
        assert_eq!(counts(item), expected, "{projection}");
        assert_eq!(item["masked"], json!({}), "{projection}");
        assert_eq!(item["not_applied"], json!({}), "{projection}");
    }
}

#[test]
fn class_is_not_applied_to_another_plugins_field() {
    let scratch = Scratch::new("other-plugin");
    let planted =
        Planted::registry(&[("LOEUF", &["LOEUF", "LOEUF_b"]), ("Other", &["LOEUF_twin"])]);
    let reference = scratch.output("reference.txt", &[row("1:100", "T", "ENST1", "IMPACT=LOW")]);
    let engine = scratch.output(
        "engine.txt",
        &[row("1:100", "T", "ENST1", "LOEUF=0;LOEUF_b=0;LOEUF_twin=0")],
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new()
            .residuals(&residuals)
            .table("planted_engine_zero\tLOEUF\tLOEUF,LOEUF_t*\n"),
    );
    for (projection, columns) in [("feature", "1:100\tT\tENST1"), ("allele", "1:100\tT")] {
        let own = entry(&report, "LOEUF", "engine", projection);
        assert_eq!(counts(own), [0, 2, 0, 0, 1, 0], "{projection}");
        assert_eq!(own["masked"], json!({"planted_engine_zero": 1}));
        assert_eq!(
            residual_lines(&residuals, "LOEUF", "engine", projection),
            [format!("engine\t{columns}\tLOEUF_b\t0")]
        );
        let other = entry(&report, "Other", "engine", projection);
        assert_eq!(counts(other), [0, 1, 0, 0, 1, 0], "{projection}");
        assert_eq!(other["masked"], json!({}));
        assert_eq!(other["not_applied"], json!({}));
        assert_eq!(
            residual_lines(&residuals, "Other", "engine", projection),
            [format!("engine\t{columns}\tLOEUF_twin\t0")]
        );
    }
    let asked: Vec<(String, String)> = planted
        .calls()
        .into_iter()
        .map(|call| (call.klass, call.name))
        .collect();
    assert_eq!(
        asked,
        [("planted_engine_zero".to_string(), "LOEUF".to_string())]
    );
}

/// The raw counts of every source's feature projection, by source.
fn raw_by_source(report: &Value) -> Value {
    let found: Map<String, Value> = report["scores"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|item| item["projection"] == "feature")
        .map(|item| {
            let source = item["source"].as_str().unwrap().to_string();
            (source, json!(raw_counts(item)))
        })
        .collect();
    Value::Object(found)
}

#[test]
fn field_two_plugins_write_counts_once_for_the_runs_own_plugin() {
    let scratch = Scratch::new("two-writers");
    let planted = Planted::registry(&[
        ("dbNSFP", &["SIFT_score", "SIFT_pred"]),
        ("TranscriptAnnotator", &["SIFT_score", "Comment"]),
    ]);
    let reference = scratch.output(
        "reference.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "SIFT_score=0.1;Comment=kept;IMPACT=LOW",
        )],
    );
    let engine = scratch.output(
        "engine.txt",
        &[row("1:100", "T", "ENST1", "SIFT_score=0.9;Comment=kept")],
    );
    let runs = [
        (
            vec!["TranscriptAnnotator"],
            json!({"TranscriptAnnotator": [2, 2, 1]}),
            json!({}),
            json!([]),
        ),
        (
            vec!["dbNSFP"],
            json!({"dbNSFP": [1, 1, 0]}),
            json!({"Comment": 1}),
            json!([]),
        ),
        (
            vec!["TranscriptAnnotator", "dbNSFP"],
            json!({"TranscriptAnnotator": [1, 1, 1]}),
            json!({"SIFT_score": 1}),
            json!(["dbNSFP"]),
        ),
    ];
    for (plugins, expected, unclaimed, not_populated) in runs {
        let residuals = scratch.path(&plugins.join("-"));
        let report = planted.score(
            &scratch,
            &reference,
            &[("engine", &engine)],
            &Run::new()
                .residuals(&residuals)
                .no_table()
                .plugins(&plugins),
        );
        assert_eq!(report["plugins"], json!(plugins));
        assert_eq!(raw_by_source(&report), expected, "{plugins:?}");
        assert_eq!(report["not_populated"], not_populated, "{plugins:?}");
        assert_eq!(
            report["reference"]["unclaimed_keys"], unclaimed,
            "{plugins:?}"
        );
        assert_eq!(
            report["engines"][0]["unclaimed_keys"], unclaimed,
            "{plugins:?}"
        );
        let mut written: Vec<String> = names(&residuals)
            .iter()
            .map(|name| name.split('.').next().unwrap().to_string())
            .collect();
        written.sort();
        let mut sources: Vec<String> = expected
            .as_object()
            .unwrap()
            .keys()
            .flat_map(|source| [source.clone(), source.clone()])
            .collect();
        sources.sort();
        assert_eq!(written, sources, "{plugins:?}");
        if let [plugin] = plugins.as_slice() {
            assert_eq!(
                feature_residuals(&residuals, plugin),
                [
                    "reference\t1:100\tT\tENST1\tSIFT_score\t0.1",
                    "engine\t1:100\tT\tENST1\tSIFT_score\t0.9",
                ]
            );
        }
    }
}

#[test]
fn field_is_asked_of_the_registry_among_the_runs_plugins_alone() {
    let scratch = Scratch::new("asked");
    let planted = Planted {
        registry: PlantedRegistry::planted().resolving(|field, _, _| {
            let depth = field.strip_prefix("depth_").is_some_and(|rest| {
                !rest.is_empty() && rest.chars().all(|digit| digit.is_ascii_digit())
            });
            Ok(depth.then(|| "Planted".to_string()))
        }),
        rules: PlantedRules::planted(),
    };
    let reference = scratch.output(
        "reference.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "depth_10=1;depth_x=2;PLANT_score=3",
        )],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new()
            .no_table()
            .plugins(&["Planted", "LOEUF", "Planted"]),
    );
    assert_eq!(report["plugins"], json!(["Planted", "LOEUF"]));
    assert_eq!(report["sources"], json!({"Planted": ["depth_10"]}));
    assert_eq!(
        report["reference"]["unclaimed_keys"],
        json!({"PLANT_score": 1, "depth_x": 1})
    );
    assert_eq!(report["not_populated"], json!(["LOEUF"]));
    let mut asked = planted.registry.asked();
    asked.sort();
    let loaded = strings(&["Planted", "LOEUF"]);
    let question = |field: &str| (field.to_string(), loaded.clone(), None);
    assert_eq!(
        asked,
        [
            question("PLANT_score"),
            question("depth_10"),
            question("depth_x")
        ]
    );
}

#[test]
fn run_names_its_sources() {
    let scratch = Scratch::new("names-sources");
    let planted = Planted::new();
    let reference = scratch.output("reference.txt", &reference_rows());
    let engines = [("engine", reference.as_str())];
    let refused = refusal(planted.try_score(
        &scratch,
        &reference,
        &engines,
        &Run::new().plugins(&["Planted", "NoSuchPlugin"]),
    ));
    assert!(
        refused.contains("--plugins names NoSuchPlugin, no plugin of the registry"),
        "{refused}"
    );
    let refused = refusal(planted.try_score(
        &scratch,
        &reference,
        &engines,
        &Run::new().plugins(&[]).no_table(),
    ));
    assert!(refused.contains("no source to score"), "{refused}");
    let report = planted.score(
        &scratch,
        &reference,
        &engines,
        &Run::new()
            .custom(&[("Custom", "Custom")])
            .plugins(&[])
            .no_table(),
    );
    let sources: Vec<&str> = report["scores"]
        .as_array()
        .unwrap()
        .iter()
        .map(|item| item["source"].as_str().unwrap())
        .collect();
    assert_eq!(sources, ["Custom", "Custom"]);
    assert_eq!(report["plugins"], json!([]));
    assert_eq!(report["not_populated"], json!([]));
    assert_eq!(
        report["reference"]["unclaimed_keys"],
        json!({
            "GWAS_pmid": 1,
            "LOEUF": 1,
            "MaxEntScan_diff": 1,
            "PLANT_label": 1,
            "PLANT_score": 3,
            "PlantedDynamic_a": 1,
            "PlantedDynamic_b": 1,
        })
    );
}

#[test]
fn custom_source_claims_only_what_no_plugin_of_the_run_claims() {
    let scratch = Scratch::new("custom-claims");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "PLANT_score=1;PLANT_other=2;PLANT_label=x",
        )],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new()
            .custom(&[("Wide", "PLANT_"), ("Wider", "PLANT")])
            .plugins(&["Planted"]),
    );
    assert_eq!(
        report["sources"],
        json!({"Planted": ["PLANT_label", "PLANT_score"], "Wide": ["PLANT_other"]})
    );
    assert_eq!(counts(feature(&report, "Planted")), [2; 6]);
    assert_eq!(counts(feature(&report, "Wide")), [1; 6]);
    assert_eq!(report["not_populated"], json!(["Wider"]));
}

#[test]
fn keys_no_source_claims_are_listed_per_input() {
    let scratch = Scratch::new("unclaimed");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row(
                "1:100",
                "T",
                "ENST1",
                "IMPACT=LOW;Mystery=1;PLANT_score=1;Bare;Blank=-",
            ),
            row("1:200", "T", "ENST1", "impact=LOW;Mystery=2;Other=x"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[row("1:100", "T", "ENST1", "PLANT_score=1;Sidecar=7")],
    );
    let listed = scratch.tuples(
        "engine.tsv",
        &[
            ["1:100", "T", "ENST1", "PLANT_score", "1"],
            ["1:100", "T", "ENST1", "IMPACT", "LOW"],
            ["1:100", "T", "ENST1", "Sidecar", "7"],
            ["1:100", "T", "ENST1", "Sidecar", "-"],
        ],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine), ("listed", &listed)],
        &Run::new(),
    );
    assert_eq!(
        report["reference"]["unclaimed_keys"],
        json!({"Mystery": 2, "Other": 1})
    );
    assert_eq!(
        report["engines"][0]["unclaimed_keys"],
        json!({"Sidecar": 1})
    );
    assert_eq!(
        report["engines"][1]["unclaimed_keys"],
        json!({"Sidecar": 1})
    );
    assert_eq!(report["engines"][1]["tuples_outside_every_source"], 2);
}

#[test]
fn unclaimed_keys_that_read_alike_keep_every_count() {
    let scratch = Scratch::new("read-alike");
    let planted = Planted::new();
    let mut written = HEADER.as_bytes().to_vec();
    written.extend_from_slice(
        b"v\t1:100\tT\tENSG1\tENST1\tTranscript\tmissense_variant\t-\t-\t-\t-\t-\t-\t",
    );
    written.extend_from_slice(b"K\xff=1;K\xfe=1;K\xfe=2;PLANT_score=1\n");
    let reference = scratch.file("reference.txt", written);
    let report = planted.score(&scratch, &reference, &[("engine", &reference)], &Run::new());
    assert_eq!(
        report["reference"]["unclaimed_keys"],
        json!({"K\u{fffd}": 3})
    );
    assert_eq!(counts(feature(&report, "Planted")), [1; 6]);
}

#[test]
fn allele_tuples_find_their_verdicts_whatever_bytes_their_values_carry() {
    let scratch = Scratch::new("value-bytes");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "GWAS_pmid=x"),
            row("1:100", "T", "ENST2", "GWAS_pmid=x\u{1}"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "IMPACT=LOW"),
            row("1:100", "T", "ENST2", "IMPACT=LOW"),
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    for (projection, keys) in [("feature", 2), ("allele", 1)] {
        let item = entry(&report, "GWAS", "engine", projection);
        assert_eq!(counts(item), [2, 0, 0, 0, 0, 0], "{projection}");
        assert_eq!(
            item["masked"],
            json!({"planted_reference_value": keys}),
            "{projection}"
        );
    }
}

/// The reference rows with a changed MaxEntScan value on the second and a `LOEUF=0` the
/// third alone writes.
fn changed_and_filled(rows: &[String]) -> Vec<String> {
    let mut engine = rows.to_vec();
    engine[1] = rows[1].replace("MaxEntScan_diff=1.5", "MaxEntScan_diff=-2.25");
    engine[2] = rows[2].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;LOEUF=0");
    engine
}

#[test]
fn classes_named_on_the_command_line_are_the_ones_applied() {
    let scratch = Scratch::new("named-classes");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let engine = scratch.output("engine.txt", &changed_and_filled(&rows));
    let engines = [("engine", engine.as_str())];
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &engines,
        &Run::new()
            .residuals(&residuals)
            .applied(&["planted_two_numbers"]),
    );
    let applied: Vec<&str> = report["classes"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|item| item["applied"] == true)
        .map(|item| item["class"].as_str().unwrap())
        .collect();
    assert_eq!(applied, ["planted_two_numbers"]);
    assert_eq!(
        report["classes"][1],
        json!({
            "class": "planted_engine_zero",
            "plugins": {"LOEUF": ["LOEUF"]},
            "applied": false,
            "further_columns": ["met", "-"],
        })
    );
    let table = scratch.path("classes.tsv");
    assert_eq!(
        report["class_table"],
        json!({"path": table, "sha256": sha256(&std::fs::read(&table).unwrap())})
    );
    for (projection, columns) in [("feature", "1:100\tT\tENST3"), ("allele", "1:100\tT")] {
        let applied = entry(&report, "MaxEntScan", "engine", projection);
        assert_eq!(counts(applied), [1, 1, 0, 0, 0, 0], "{projection}");
        assert_eq!(applied["masked"], json!({"planted_two_numbers": 1}));
        assert_eq!(applied["not_applied"], json!({}));
        let held = entry(&report, "LOEUF", "engine", projection);
        assert_eq!(raw_counts(held), [1, 2, 1], "{projection}");
        let keys = held.as_object().unwrap().keys();
        assert!(!keys.into_iter().any(|key| key.starts_with("adj_")));
        assert_eq!(held["masked"], json!({}));
        assert_eq!(held["not_applied"], json!({"planted_engine_zero": 1}));
        assert_eq!(held["residual_rows"], 1);
        assert_eq!(
            residual_lines(&residuals, "LOEUF", "engine", projection),
            [format!("engine\t{columns}\tLOEUF\t0")]
        );
    }
    let nothing = planted.score(&scratch, &reference, &engines, &Run::new().applied(&[]));
    let classes = nothing["classes"].as_array().unwrap();
    assert!(!classes.iter().any(|item| item["applied"] == true));
    assert_eq!(
        feature(&nothing, "MaxEntScan")["not_applied"],
        json!({"planted_two_numbers": 1})
    );
    let refused = refusal(planted.try_score(
        &scratch,
        &reference,
        &engines,
        &Run::new().applied(&["planted_missing"]),
    ));
    assert!(
        refused.contains("--classes names planted_missing"),
        "{refused}"
    );
}

#[test]
fn applied_class_is_asked_before_one_that_is_not() {
    let scratch = Scratch::new("applied-first");
    let planted = Planted::registry(&[("GWAS", &["GWAS_pmid"])]);
    let table = "planted_reference_value\tGWAS\tGWAS_*\nplanted_two_numbers\tGWAS\tGWAS_pmid\n";
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "T", "ENST1", "GWAS_pmid=17")],
    );
    let engine = scratch.output("engine.txt", &[row("1:100", "T", "ENST1", "GWAS_pmid=18")]);
    let engines = [("engine", engine.as_str())];
    let first = planted.score(&scratch, &reference, &engines, &Run::new().table(table));
    assert_eq!(
        feature(&first, "GWAS")["masked"],
        json!({"planted_reference_value": 1, "planted_two_numbers": 0})
    );
    let second = planted.score(
        &scratch,
        &reference,
        &engines,
        &Run::new().table(table).applied(&["planted_two_numbers"]),
    );
    let second = feature(&second, "GWAS");
    assert_eq!(second["masked"], json!({"planted_two_numbers": 1}));
    assert_eq!(second["not_applied"], json!({}));
    assert_eq!(counts(second), [1, 1, 0, 0, 0, 0]);
}

#[test]
fn no_class_table_declares_no_class() {
    let scratch = Scratch::new("no-table");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let mut filled = rows.clone();
    filled[2] = rows[2].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;LOEUF=0");
    let engine = scratch.output("engine.txt", &filled);
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().no_table(),
    );
    assert_eq!(report["classes"], json!([]));
    assert!(report["class_table"].is_null());
    let item = feature(&report, "LOEUF");
    assert_eq!(counts(item), [1, 2, 1, 1, 2, 1]);
    assert_eq!(item["masked"], json!({}));
    assert_eq!(item["not_applied"], json!({}));
    assert!(planted.calls().is_empty());
}

#[test]
fn malformed_class_table_is_refused() {
    let cases = [
        (
            "planted_engine_zero\tLOEUF\n",
            "line 1 lacks a class, its plugins and their field names",
        ),
        (
            "#header\nno_such_class\tLOEUF\tLOEUF\n",
            "line 2 names no_such_class, no class of",
        ),
        (
            "planted_engine_zero\tLOEUF\tLOEUF\nplanted_engine_zero\tGWAS\tGWAS_*\n",
            "line 2 declares planted_engine_zero",
        ),
        (
            "planted_engine_zero\tLOEUF;GWAS\tLOEUF\n",
            "line 1 gives 2 plugins and 1 field lists",
        ),
        (
            "planted_engine_zero\tNoSuchPlugin\tLOEUF\n",
            "line 1 names 'NoSuchPlugin', no source",
        ),
        (
            "planted_engine_zero\tLOEUF\t ,\n",
            "line 1 names 'LOEUF', no source, or no field of it",
        ),
        (
            "planted_engine_zero\tLOEUF; LOEUF\tLOEUF;LOEUF_b\n",
            "line 1 names 'LOEUF' twice",
        ),
    ];
    for (table, message) in cases {
        let scratch = Scratch::new("malformed-table");
        let planted = Planted::new();
        let reference = scratch.output("reference.txt", &reference_rows());
        let refused = refusal(planted.try_score(
            &scratch,
            &reference,
            &[("engine", &reference)],
            &Run::new().table(table),
        ));
        assert!(refused.contains(message), "{refused}");
        assert!(!scratch.exists("work"), "{message}");
    }
}

#[test]
fn class_table_that_is_not_text_is_refused() {
    let scratch = Scratch::new("table-bytes");
    let planted = Planted::new();
    let reference = scratch.output("reference.txt", &reference_rows());
    let table = scratch.file("classes.tsv", b"planted_engine_zero\tLOEUF\tLOEUF\xff\n");
    let refused = refusal(planted.try_score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new().table_file(&table).plugins(&["LOEUF"]),
    ));
    assert!(
        refused.contains("classes.tsv: the class table is not UTF-8 text"),
        "{refused}"
    );
    assert!(!scratch.exists("work"));
}

#[test]
fn star_is_the_class_tables_only_wildcard() {
    let scratch = Scratch::new("wildcard");
    let planted = Planted::registry(&[("LOEUF", &["LOEUF", "LOEUF_b", "LOEU?", "L.EUF"])]);
    let reference = scratch.output("reference.txt", &[row("1:100", "T", "ENST1", "IMPACT=LOW")]);
    let engine = scratch.output(
        "engine.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "LOEUF=0;LOEUF_b=0;LOEU?=0;L.EUF=0",
        )],
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new()
            .residuals(&residuals)
            .table("planted_engine_zero\tLOEUF\tLOEU?,[A-Z]OEUF,LO*F_b,L.EUF\n"),
    );
    let item = feature(&report, "LOEUF");
    assert_eq!(counts(item), [0, 4, 0, 0, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_engine_zero": 3}));
    assert_eq!(
        feature_residuals(&residuals, "LOEUF"),
        ["engine\t1:100\tT\tENST1\tLOEUF\t0"]
    );
    assert_eq!(
        report["classes"][0]["plugins"],
        json!({"LOEUF": ["LOEU?", "[A-Z]OEUF", "LO*F_b", "L.EUF"]})
    );
}

#[test]
fn class_of_two_plugins_reads_each_plugins_own_fields() {
    let scratch = Scratch::new("two-plugins");
    let planted = Planted::new();
    let reference = scratch.output("reference.txt", &[row("1:100", "T", "ENST1", "IMPACT=LOW")]);
    let engine = scratch.output(
        "engine.txt",
        &[row("1:100", "T", "ENST1", "LOEUF=0;GWAS_pmid=0")],
    );
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().table("planted_engine_zero\tLOEUF;GWAS\tLOEUF;GWAS_pmid\n"),
    );
    assert_eq!(
        report["classes"],
        json!([{
            "class": "planted_engine_zero",
            "plugins": {"LOEUF": ["LOEUF"], "GWAS": ["GWAS_pmid"]},
            "applied": true,
            "further_columns": [],
        }])
    );
    for source in ["LOEUF", "GWAS"] {
        let item = feature(&report, source);
        assert_eq!(counts(item), [0, 1, 0, 0, 0, 0], "{source}");
        assert_eq!(
            item["masked"],
            json!({"planted_engine_zero": 1}),
            "{source}"
        );
    }
}

#[test]
fn residual_files_are_empty_when_a_reference_scores_against_itself() {
    let scratch = Scratch::new("empty-residuals");
    let planted = Planted::new();
    let reference = scratch.output("reference.txt", &reference_rows());
    let residuals = scratch.path("kept/residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new()
            .custom(&[("Custom", "Custom")])
            .residuals(&residuals),
    );
    let scores = report["scores"].as_array().unwrap();
    let written = names(&residuals);
    let mut stated: Vec<String> = scores
        .iter()
        .map(|item| {
            format!(
                "{}.engine.{}.tsv",
                item["source"].as_str().unwrap(),
                item["projection"].as_str().unwrap()
            )
        })
        .collect();
    stated.sort();
    assert_eq!(written, stated);
    assert_eq!(written.len(), 12);
    for name in &written {
        let content = std::fs::read(format!("{residuals}/{name}")).unwrap();
        assert!(content.is_empty(), "{name}");
    }
    assert!(scores.iter().all(|item| item["residual_rows"] == 0));
}

#[test]
fn residual_file_has_one_line_per_side_for_an_undocumented_difference() {
    let scratch = Scratch::new("residual-lines");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let mut changed = rows.clone();
    changed[0] = rows[0].replace("PLANT_label=benign", "PLANT_label=pathogenic");
    changed[3] = rows[3].replace("PLANT_score=0.25", "PLANT_score=0.125");
    let engine = scratch.output("engine.txt", &changed);
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().residuals(&residuals),
    );
    assert_eq!(
        feature_residuals(&residuals, "Planted"),
        [
            "reference\t1:100\tT\tENST1\tPLANT_label\tbenign",
            "engine\t1:100\tT\tENST1\tPLANT_label\tpathogenic",
            "engine\t1:200-201\t-\tENST1\tPLANT_score\t0.125",
            "reference\t1:200-201\t-\tENST1\tPLANT_score\t0.25",
        ]
    );
    assert_eq!(
        residual_lines(&residuals, "Planted", "engine", "allele"),
        [
            "reference\t1:100\tT\tPLANT_label\tbenign",
            "engine\t1:100\tT\tPLANT_label\tpathogenic",
            "engine\t1:200-201\t-\tPLANT_score\t0.125",
            "reference\t1:200-201\t-\tPLANT_score\t0.25",
        ]
    );
    assert_eq!(feature(&report, "Planted")["residual_rows"], 4);
    assert_eq!(allele(&report, "Planted")["residual_rows"], 4);
    assert!(feature_residuals(&residuals, "LOEUF").is_empty());
    assert_eq!(feature(&report, "LOEUF")["residual_rows"], 0);
}

#[test]
fn residual_file_has_no_line_for_a_documented_difference() {
    let scratch = Scratch::new("documented-residuals");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let engine = scratch.output("engine.txt", &changed_and_filled(&rows));
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().residuals(&residuals),
    );
    for source in ["MaxEntScan", "LOEUF"] {
        for projection in ["feature", "allele"] {
            let item = entry(&report, source, "engine", projection);
            let masked = item["masked"].as_object().unwrap();
            assert!(
                masked.values().any(|keys| keys != 0),
                "{source} {projection}"
            );
            assert!(item["intersection"].as_i64() < item["rust"].as_i64());
            assert_eq!(item["residual_rows"], 0, "{source} {projection}");
            assert!(residual_lines(&residuals, source, "engine", projection).is_empty());
        }
    }
}

#[test]
fn residuals_need_names_a_file_name_can_carry() {
    let scratch = Scratch::new("residual-names");
    let planted = Planted::new();
    let reference = scratch.output("reference.txt", &reference_rows());
    let engines = [("../engine", reference.as_str())];
    let refused = refusal(planted.try_score(
        &scratch,
        &reference,
        &engines,
        &Run::new().residuals(&scratch.path("residuals")),
    ));
    assert!(
        refused.contains("cannot be part of a residual file name"),
        "{refused}"
    );
    assert!(!scratch.exists("residuals"));
    let report = planted.score(&scratch, &reference, &engines, &Run::new());
    assert_eq!(
        counts(entry(&report, "Planted", "../engine", "feature")),
        [4; 6]
    );
}

#[test]
fn residual_file_names_may_not_collide() {
    let scratch = Scratch::new("residual-collision");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "T", "ENST1", "X_a=1;Y_a=1")],
    );
    let engines = [("B.C", reference.as_str()), ("C", reference.as_str())];
    let run = Run::new()
        .no_table()
        .plugins(&[])
        .custom(&[("A", "X_"), ("A.B", "Y_")]);
    let kept = run.clone().residuals(&scratch.path("kept"));
    let refused = refusal(planted.try_score(&scratch, &reference, &engines, &kept));
    assert!(refused.contains("share a residual file name"), "{refused}");
    let twice = [engines[1], engines[1]];
    let refused = refusal(planted.try_score(&scratch, &reference, &twice, &kept));
    assert!(refused.contains("share a residual file name"), "{refused}");
    assert!(!scratch.exists("kept"));
    let report = planted.score(&scratch, &reference, &engines, &run);
    assert_eq!(report["scores"].as_array().unwrap().len(), 8);
}

#[test]
fn residuals_directory_of_an_earlier_run_is_refused() {
    let scratch = Scratch::new("earlier-residuals");
    let planted = Planted::new();
    let reference = scratch.output("reference.txt", &reference_rows());
    let residuals = scratch.path("residuals");
    let run = Run::new().residuals(&residuals);
    planted.score(&scratch, &reference, &[("engine", &reference)], &run);
    let earlier = names(&residuals);
    assert!(!earlier.is_empty());
    let refused = refusal(planted.try_score(&scratch, &reference, &[("later", &reference)], &run));
    assert!(
        refused.contains("holds .tsv files of an earlier run"),
        "{refused}"
    );
    assert_eq!(names(&residuals), earlier);
}

#[test]
fn residuals_directory_that_cannot_be_written_is_refused_before_scoring() {
    // SAFETY: `geteuid` reads the process's user and cannot fail.
    if unsafe { libc::geteuid() } == 0 {
        // The superuser writes in any directory.
        return;
    }
    let scratch = Scratch::new("locked-residuals");
    let planted = Planted::new();
    let reference = scratch.output("reference.txt", &reference_rows());
    let residuals = scratch.path("residuals");
    std::fs::create_dir(&residuals).unwrap();
    std::fs::set_permissions(&residuals, std::fs::Permissions::from_mode(0o555)).unwrap();
    let scored = planted.try_score(
        &scratch,
        &reference,
        &[("engine", &reference)],
        &Run::new().residuals(&residuals),
    );
    std::fs::set_permissions(&residuals, std::fs::Permissions::from_mode(0o755)).unwrap();
    let refused = refusal(scored);
    assert!(refused.contains("cannot write in --residuals"), "{refused}");
    assert!(scratch.exists("work") && scratch.names("work").is_empty());
}
