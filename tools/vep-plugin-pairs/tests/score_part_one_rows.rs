// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows the `score` program consults a class with, under a planted registry and
//! planted classes: the row that wrote a value, the row that wrote none, the rows of the
//! same input record, and how the values left under a key two records share are paired.

mod support;

use serde_json::json;
use std::collections::BTreeMap;
use support::{
    allele, counts, entry, feature, feature_residuals, residual_lines, row, row_with, Call,
    Planted, PlantedRegistry, PlantedRules, Run, Scratch,
};
use vep_plugin_pairs::rules::{Entry, ABSENT};

/// A row of the record `uploaded` that writes no cDNA position.
fn named(uploaded: &str, location: &str, allele: &str, feature: &str, extra: &str) -> String {
    row_with(location, allele, feature, extra, "-", uploaded)
}

fn uploaded(row: &Entry) -> &str {
    row.fields.get("Uploaded_variation").unwrap()
}

#[test]
fn missing_and_wrong_values_on_a_one_allele_row_are_residual_under_a_declared_class() {
    let scratch = Scratch::new("one-allele");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            named("rec1", "1:100", "T", "ENST1", "REVEL=0.51"),
            named("rec2", "1:200", "T", "ENST1", "REVEL=0.51"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            named("rec1", "1:100", "T", "ENST1", "IMPACT=LOW"),
            named("rec2", "1:200", "T", "ENST1", "REVEL=0.99"),
        ],
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().residuals(&residuals),
    );
    let item = feature(&report, "REVEL");
    assert_eq!(counts(item), [2, 1, 0, 2, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_sibling_value": 0}));
    assert_eq!(item["not_applied"], json!({}));
    assert_eq!(item["residual_rows"], 3);
    assert_eq!(
        feature_residuals(&residuals, "REVEL"),
        [
            "reference\t1:100\tT\tENST1\tREVEL\t0.51",
            "reference\t1:200\tT\tENST1\tREVEL\t0.51",
            "engine\t1:200\tT\tENST1\tREVEL\t0.99",
        ]
    );
    let calls = planted.calls();
    let asked: Vec<(&str, &str, &str)> = calls
        .iter()
        .map(|call| {
            (
                call.klass.as_str(),
                call.expected.as_str(),
                call.actual.as_str(),
            )
        })
        .collect();
    assert_eq!(
        asked,
        [
            ("planted_sibling_value", "0.51", ABSENT),
            ("planted_sibling_value", "0.51", "0.99"),
        ]
    );
    assert!(calls.iter().all(|call| call.record_rows.is_none()));
}

#[test]
fn predicate_is_handed_the_reference_rows_of_the_same_record() {
    let scratch = Scratch::new("record-rows");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            named("1_300_C/A/G", "1:300", "A", "ENST1", "REVEL=0.7"),
            named("1_300_C/A/G", "1:300", "A", "ENST2", "REVEL=0.7"),
            named("1_300_C/A/G", "1:300", "G", "ENST1", "REVEL=0.7"),
            named("1_300_G/C", "1:300", "C", "ENST1", "REVEL=0.2"),
            named("1_300_G/T", "1:300", "T", "ENST1", "REVEL=0.7"),
            named("1_400_A/G", "1:400", "G", "ENST1", "REVEL=0.7"),
            named("1_500_C/A/G", "1:500", "A", "ENST1", "REVEL=0.7"),
            named("1_500_C/A/G", "1:500", "A", "ENST1", "IMPACT=LOW"),
            named("1_500_C/A/G", "1:500", "G", "ENST1", "REVEL=0.7"),
            named("rs6", "1:600", "A", "ENST1", "REVEL=0.7"),
            named("rs6", "1:600", "G", "ENST1", "REVEL=0.7"),
        ],
    );
    let other = |location: &str, allele: &str, feature: &str, extra: &str| {
        named("other-name", location, allele, feature, extra)
    };
    let engine = scratch.output(
        "engine.txt",
        &[
            other("1:300", "A", "ENST1", "REVEL=0.7"),
            other("1:300", "A", "ENST2", "REVEL=0.7"),
            other("1:300", "G", "ENST1", "IMPACT=LOW"),
            other("1:300", "C", "ENST1", "REVEL=0.2"),
            other("1:300", "T", "ENST1", "IMPACT=LOW"),
            other("1:400", "G", "ENST1", "IMPACT=LOW"),
            other("1:500", "A", "ENST1", "REVEL=0.7"),
            other("1:500", "G", "ENST1", "IMPACT=LOW"),
            other("1:600", "A", "ENST1", "REVEL=0.7"),
            other("1:600", "G", "ENST1", "IMPACT=LOW"),
        ],
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().residuals(&residuals),
    );
    let item = feature(&report, "REVEL");
    assert_eq!(counts(item), [10, 5, 5, 9, 5, 5]);
    assert_eq!(item["masked"], json!({"planted_sibling_value": 1}));
    assert_eq!(
        feature_residuals(&residuals, "REVEL"),
        [
            "reference\t1:300\tT\tENST1\tREVEL\t0.7",
            "reference\t1:400\tG\tENST1\tREVEL\t0.7",
            "reference\t1:500\tG\tENST1\tREVEL\t0.7",
            "reference\t1:600\tG\tENST1\tREVEL\t0.7",
        ]
    );
    let calls = planted.calls();
    let handed: BTreeMap<(&str, &str), Option<&Vec<Entry>>> = calls
        .iter()
        .map(|call| {
            (
                (call.e.location.as_deref().unwrap(), call.e.allele.as_str()),
                call.record_rows.as_ref(),
            )
        })
        .collect();
    let rows = |location: &'static str, allele: &'static str| handed[&(location, allele)];
    let both: Vec<(&str, &str)> = rows("1:300", "G")
        .unwrap()
        .iter()
        .map(|other| (other.allele.as_str(), other.feature.as_str()))
        .collect();
    assert_eq!(both, [("A", "ENST1"), ("G", "ENST1")]);
    let own: Vec<&str> = rows("1:300", "T").unwrap().iter().map(uploaded).collect();
    assert_eq!(own, ["1_300_G/T"]);
    let alone: Vec<&str> = rows("1:400", "G")
        .unwrap()
        .iter()
        .map(|other| other.allele.as_str())
        .collect();
    assert_eq!(alone, ["G"]);
    assert!(rows("1:500", "G").is_none());
    assert!(rows("1:600", "G").is_none());
    assert!(calls
        .iter()
        .all(|call| uploaded(call.a.as_ref().unwrap()) == "other-name"));
    let by_allele = allele(&report, "REVEL");
    assert_eq!(counts(by_allele), [9, 4, 4, 8, 4, 4]);
    assert_eq!(by_allele["masked"], json!({"planted_sibling_value": 1}));
}

#[test]
fn class_is_not_consulted_without_a_row_that_wrote_the_value_or_none() {
    let scratch = Scratch::new("no-row");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            named("v1", "1:100", "T", "ENST1", "LOEUF=0.31"),
            named("v3", "1:200", "T", "ENST1", "GWAS_pmid=17"),
            named("v4", "1:200", "T", "ENST1", "GWAS_pmid=18"),
            named("v5", "1:300", "T", "ENST1", "LOEUF="),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            named("v1", "1:100", "T", "ENST1", "LOEUF=0.31"),
            named("v2", "1:100", "T", "ENST1", "LOEUF=0"),
            named("v3", "1:200", "T", "ENST1", "GWAS_pmid=17"),
            named("v5", "1:300", "T", "ENST1", "LOEUF=0"),
        ],
    );
    let listed = scratch.tuples(
        "engine.tsv",
        &[
            ["1:100", "T", "ENST1", "LOEUF", "0.31"],
            ["1:100", "T", "ENST1", "LOEUF", "0"],
            ["1:200", "T", "ENST1", "GWAS_pmid", "17"],
            ["1:300", "T", "ENST1", "LOEUF", "0"],
        ],
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine), ("listed", &listed)],
        &Run::new().residuals(&residuals),
    );
    for ident in ["engine", "listed"] {
        let zero = entry(&report, "LOEUF", ident, "feature");
        assert_eq!(counts(zero), [1, 3, 1, 1, 2, 1], "{ident}");
        assert_eq!(zero["masked"], json!({"planted_engine_zero": 1}), "{ident}");
        assert_eq!(
            residual_lines(&residuals, "LOEUF", ident, "feature"),
            ["engine\t1:100\tT\tENST1\tLOEUF\t0"],
            "{ident}"
        );
        let left = entry(&report, "GWAS", ident, "feature");
        assert_eq!(counts(left), [2, 1, 1, 2, 1, 1], "{ident}");
        assert_eq!(
            left["masked"],
            json!({"planted_reference_value": 0}),
            "{ident}"
        );
        assert_eq!(
            residual_lines(&residuals, "GWAS", ident, "feature"),
            ["reference\t1:200\tT\tENST1\tGWAS_pmid\t18"],
            "{ident}"
        );
    }
    let calls = planted.calls();
    // The class, the reference row's Location and its own LOEUF value, and the two values.
    let asked: Vec<[&str; 5]> = calls
        .iter()
        .map(|call| {
            [
                call.klass.as_str(),
                call.e.location.as_deref().unwrap(),
                call.e.fields.get("LOEUF").unwrap(),
                call.expected.as_str(),
                call.actual.as_str(),
            ]
        })
        .collect();
    assert_eq!(
        asked,
        [["planted_engine_zero", "1:300", "", ABSENT, "0"]; 2]
    );
}

#[test]
fn items_in_another_order_are_set_aside_as_order_only() {
    let scratch = Scratch::new("order-only");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row(
                "1:100",
                "T",
                "ENST1",
                "NearestExonJB=ENSE1|5|start|90,ENSE2|5|end|70",
            ),
            row("1:300", "C", "ENST1", "NearestExonJB=ENSE3|1|start|10"),
            row("1:400", "C", "ENST1", "NearestExonJB=ENSE4|1|start|10"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row(
                "1:100",
                "T",
                "ENST1",
                "NearestExonJB=ENSE2|5|end|70,ENSE1|5|start|90",
            ),
            row("1:300", "C", "ENST1", "NearestExonJB=ENSE3|2|start|10"),
            row("1:400", "C", "ENST1", "IMPACT=LOW"),
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    let item = feature(&report, "NearestExonJB");
    assert_eq!(counts(item), [3, 2, 0, 2, 1, 0]);
    assert_eq!(item["masked"], json!({"order_only": 1}));
}

#[test]
fn field_class_is_consulted_only_where_both_sides_wrote_the_row() {
    let scratch = Scratch::new("both-rows");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "GWAS_pmid=17"),
            row("1:100", "T", "ENST2", "GWAS_pmid=18"),
        ],
    );
    let without_the_field =
        scratch.output("engine.txt", &[row("1:100", "T", "ENST1", "IMPACT=LOW")]);
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &without_the_field)],
        &Run::new(),
    );
    let item = feature(&report, "GWAS");
    assert_eq!(counts(item), [2, 0, 0, 1, 0, 0]);
    assert_eq!(item["masked"], json!({"planted_reference_value": 1}));
}

#[test]
fn field_class_reads_the_row_that_wrote_the_value() {
    let scratch = Scratch::new("writing-row");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "5UTR_consequence=first"),
            row_with("1:100", "T", "ENST1", "5UTR_consequence=second", "1", "v"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "5UTR_consequence=third"),
            row("1:100", "T", "ENST1", "IMPACT=LOW"),
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    let item = feature(&report, "UTRAnnotator");
    assert_eq!(counts(item), [2, 1, 0, 1, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_row_cdna_one": 1}));
}

#[test]
fn row_that_wrote_no_value_is_read_for_one_unpaired_tuple_alone() {
    let scratch = Scratch::new("no-value-row");
    let planted = Planted::new();
    let key = |extra: &str| row("1:100", "T", "ENST1", extra);
    let once = scratch.output("once.txt", &[key("IMPACT=LOW")]);
    let twice = scratch.output("twice.txt", &[key("IMPACT=LOW"), key("IMPACT=LOW")]);
    let zeros = scratch.output("zeros.txt", &[key("LOEUF=0"), key("LOEUF=0")]);
    let listed = scratch.tuples("zeros.tsv", &[["1:100", "T", "ENST1", "LOEUF", "0"]; 2]);
    let engines = [("engine", zeros.as_str()), ("listed", listed.as_str())];
    let residuals = scratch.path("residuals");
    let report = planted.score(&scratch, &once, &engines, &Run::new().residuals(&residuals));
    for ident in ["engine", "listed"] {
        let item = entry(&report, "LOEUF", ident, "feature");
        assert_eq!(counts(item), [0, 2, 0, 0, 1, 0], "{ident}");
        assert_eq!(item["masked"], json!({"planted_engine_zero": 1}), "{ident}");
        assert_eq!(
            residual_lines(&residuals, "LOEUF", ident, "feature"),
            ["engine\t1:100\tT\tENST1\tLOEUF\t0"],
            "{ident}"
        );
        let by_allele = entry(&report, "LOEUF", ident, "allele");
        assert_eq!(counts(by_allele), [0, 1, 0, 0, 1, 0], "{ident}");
        assert_eq!(
            by_allele["masked"],
            json!({"planted_engine_zero": 0}),
            "{ident}"
        );
    }
    let both = planted.score(&scratch, &twice, &engines, &Run::new());
    for ident in ["engine", "listed"] {
        assert_eq!(
            counts(entry(&both, "LOEUF", ident, "feature")),
            [0, 2, 0, 0, 0, 0],
            "{ident}"
        );
    }
    let engines = [("once", once.as_str()), ("twice", twice.as_str())];
    let pmids = scratch.output("pmids.txt", &[key("GWAS_pmid=17"), key("GWAS_pmid=17")]);
    let kept = scratch.path("written");
    let written = planted.score(&scratch, &pmids, &engines, &Run::new().residuals(&kept));
    let left = entry(&written, "GWAS", "once", "feature");
    assert_eq!(counts(left), [2, 0, 0, 1, 0, 0]);
    assert_eq!(left["masked"], json!({"planted_reference_value": 1}));
    assert_eq!(
        residual_lines(&kept, "GWAS", "once", "feature"),
        ["reference\t1:100\tT\tENST1\tGWAS_pmid\t17"]
    );
    assert_eq!(
        counts(entry(&written, "GWAS", "twice", "feature")),
        [2, 0, 0, 0, 0, 0]
    );
    let distinct = scratch.output("distinct.txt", &[key("GWAS_pmid=17"), key("GWAS_pmid=18")]);
    let kept = scratch.path("apart");
    let apart = planted.score(&scratch, &distinct, &engines, &Run::new().residuals(&kept));
    assert_eq!(
        residual_lines(&kept, "GWAS", "once", "feature"),
        ["reference\t1:100\tT\tENST1\tGWAS_pmid\t18"]
    );
    assert_eq!(
        counts(entry(&apart, "GWAS", "twice", "feature")),
        [2, 0, 0, 0, 0, 0]
    );
}

#[test]
fn each_unpaired_copy_of_a_value_is_judged_with_a_row_that_wrote_it() {
    for positions in [["1", "-"], ["-", "1"]] {
        let scratch = Scratch::new("copies");
        let planted = Planted::new();
        let reference = scratch.output(
            "reference.txt",
            &positions.map(|position| {
                row_with("1:100", "T", "ENST1", "5UTR_consequence=x", position, "v")
            }),
        );
        let unwritten = row("1:100", "T", "ENST1", "IMPACT=LOW");
        let engine = scratch.output("engine.txt", &[unwritten.clone(), unwritten]);
        let residuals = scratch.path("residuals");
        let report = planted.score(
            &scratch,
            &reference,
            &[("engine", &engine)],
            &Run::new().residuals(&residuals),
        );
        let item = feature(&report, "UTRAnnotator");
        assert_eq!(counts(item), [2, 0, 0, 1, 0, 0], "{positions:?}");
        assert_eq!(
            item["masked"],
            json!({"planted_row_cdna_one": 1}),
            "{positions:?}"
        );
        assert_eq!(
            feature_residuals(&residuals, "UTRAnnotator"),
            ["reference\t1:100\tT\tENST1\t5UTR_consequence\tx"],
            "{positions:?}"
        );
        let calls = planted.calls();
        let read: Vec<&str> = calls
            .iter()
            .map(|call| call.e.fields.get("cDNA_position").unwrap())
            .collect();
        assert_eq!(read, positions);
    }
}

#[test]
fn copies_of_a_value_from_two_records_are_each_handed_their_own_records_rows() {
    let scratch = Scratch::new("two-records");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            named("1_300_C/A/G", "1:300", "A", "ENST1", "REVEL=0.7"),
            named("1_300_C/A/G", "1:300", "G", "ENST1", "REVEL=0.7"),
            named("1_300_C/G", "1:300", "G", "ENST1", "REVEL=0.7"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:300", "A", "ENST1", "REVEL=0.7"),
            row("1:300", "G", "ENST1", "IMPACT=LOW"),
            row("1:300", "G", "ENST1", "IMPACT=LOW"),
        ],
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().residuals(&residuals),
    );
    let item = feature(&report, "REVEL");
    assert_eq!(counts(item), [3, 1, 1, 2, 1, 1]);
    assert_eq!(item["masked"], json!({"planted_sibling_value": 1}));
    assert_eq!(
        feature_residuals(&residuals, "REVEL"),
        ["reference\t1:300\tG\tENST1\tREVEL\t0.7"]
    );
    let calls = planted.calls();
    let handed: Vec<Vec<&str>> = calls
        .iter()
        .map(|call| {
            let rows = call.record_rows.as_ref().unwrap();
            rows.iter().map(uploaded).collect()
        })
        .collect();
    assert_eq!(
        handed,
        [vec!["1_300_C/A/G", "1_300_C/A/G"], vec!["1_300_C/G"]]
    );
}

/// The planted registry under one class, `planted_reversed_value`, which accepts an
/// engine's value that is the reference's read backwards.
fn reversed_value() -> Planted {
    Planted {
        registry: PlantedRegistry::planted(),
        rules: PlantedRules::none().class("planted_reversed_value", &[], |c| {
            c.expected == c.actual.chars().rev().collect::<String>()
        }),
    }
}

/// What each consult was put: the reference's value, the engine's, and the record of the
/// reference's row and of the engine's.
fn judged(calls: &[Call]) -> Vec<(&str, &str, &str, Option<&str>)> {
    calls
        .iter()
        .map(|call| {
            (
                call.expected.as_str(),
                call.actual.as_str(),
                uploaded(&call.e),
                call.a
                    .as_ref()
                    .and_then(|row| row.fields.get("Uploaded_variation")),
            )
        })
        .collect()
}

/// Two input records that print one Location and Allele share a key on a Feature. Their
/// values left pair by record, the rows' `Uploaded_variation`, whatever order the bytes of
/// the values give: a value with the other side's value of its record, else against its
/// record's row that wrote none, never with another record's value. A tuple file states no
/// rows, so its values pair in their order.
#[test]
fn values_left_under_a_key_two_records_share_pair_by_record() {
    let scratch = Scratch::new("pair-by-record");
    let planted = reversed_value();
    let table = "planted_reversed_value\tUTRAnnotator\t5UTR_consequence\tmet\t-\norder_only\tNearestExonJB\tNearestExonJB\tmet\t-\n";
    let both = |uploaded: &str, extra: &str| named(uploaded, "1:100", "T", "ENST1", extra);
    let one = |uploaded: &str, extra: &str| named(uploaded, "1:200", "T", "ENST1", extra);
    let items = |uploaded: &str, extra: &str| named(uploaded, "1:300", "C", "ENST1", extra);
    let reference = scratch.output(
        "reference.txt",
        &[
            both("rs1", "5UTR_consequence=ay"),
            both("rs2", "5UTR_consequence=bx"),
            one("rs3", "IMPACT=LOW"),
            one("rs4", "5UTR_consequence=cz"),
            items("rs5", "NearestExonJB=ENSE1|1|start|10,ENSE4|4|end|40"),
            items("rs6", "NearestExonJB=ENSE2|2|start|20,ENSE3|3|end|30"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            both("rs1", "5UTR_consequence=ya"),
            both("rs2", "5UTR_consequence=xb"),
            one("rs3", "5UTR_consequence=zc"),
            one("rs4", "IMPACT=LOW"),
            items("rs5", "NearestExonJB=ENSE4|4|end|40,ENSE1|1|start|10"),
            items("rs6", "NearestExonJB=ENSE3|3|end|30,ENSE2|2|start|20"),
        ],
    );
    let listed = scratch.tuples(
        "engine.tsv",
        &[
            ["1:100", "T", "ENST1", "5UTR_consequence", "ya"],
            ["1:100", "T", "ENST1", "5UTR_consequence", "xb"],
        ],
    );
    let residuals = scratch.path("residuals");
    let run = Run::new()
        .table(table)
        .plugins(&["UTRAnnotator", "NearestExonJB"]);
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &run.clone().residuals(&residuals),
    );
    let item = feature(&report, "UTRAnnotator");
    assert_eq!(counts(item), [3, 3, 0, 1, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_reversed_value": 1}));
    assert_eq!(
        feature_residuals(&residuals, "UTRAnnotator"),
        [
            "reference\t1:200\tT\tENST1\t5UTR_consequence\tcz",
            "engine\t1:200\tT\tENST1\t5UTR_consequence\tzc",
        ]
    );
    let with_rows = planted.calls();
    let mut put = judged(&with_rows);
    put.sort();
    assert_eq!(
        put,
        [
            (ABSENT, "zc", "rs3", Some("rs3")),
            ("ay", "ya", "rs1", Some("rs1")),
            ("bx", "xb", "rs2", Some("rs2")),
            ("cz", ABSENT, "rs4", Some("rs4")),
        ]
    );
    let in_order = feature(&report, "NearestExonJB");
    assert_eq!(counts(in_order), [2, 2, 0, 0, 0, 0]);
    assert_eq!(in_order["masked"], json!({"order_only": 1}));
    let report = planted.score(&scratch, &reference, &[("listed", &listed)], &run);
    let without_rows = entry(&report, "UTRAnnotator", "listed", "feature");
    assert_eq!(counts(without_rows), [3, 2, 0, 3, 2, 0]);
    assert_eq!(without_rows["masked"], json!({"planted_reversed_value": 0}));
    let calls = planted.calls();
    let values: Vec<(&str, &str)> = judged(&calls[with_rows.len()..])
        .iter()
        .map(|pair| (pair.0, pair.1))
        .collect();
    assert_eq!(values, [("ay", "xb"), ("bx", "ya"), ("cz", ABSENT)]);
}

/// Where two rows of a side under a shared key write no value of the field, a value the
/// other side has left stands against the row of its own record among them, the second as
/// readily as the first. A key that the rows of one side alone share pairs by record as
/// well: the value whose record the other side wrote no row for is put to no class.
#[test]
fn value_left_stands_against_its_records_row_among_the_rows_that_wrote_none() {
    let scratch = Scratch::new("rows-without-value");
    let planted = reversed_value();
    let three = |uploaded: &str, extra: &str| named(uploaded, "1:100", "T", "ENST1", extra);
    let one_side = |uploaded: &str, extra: &str| named(uploaded, "1:200", "T", "ENST1", extra);
    let reference = scratch.output(
        "reference.txt",
        &[
            three("rs1", "IMPACT=LOW"),
            three("rs2", "IMPACT=LOW"),
            three("rs3", "5UTR_consequence=cz"),
            one_side("rs4", "5UTR_consequence=ay"),
            one_side("rs5", "5UTR_consequence=bx"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            three("rs1", "IMPACT=LOW"),
            three("rs2", "5UTR_consequence=yb"),
            three("rs3", "IMPACT=LOW"),
            one_side("rs5", "5UTR_consequence=xb"),
        ],
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new()
            .residuals(&residuals)
            .table("planted_reversed_value\tUTRAnnotator\t5UTR_consequence\tmet\t-\n")
            .plugins(&["UTRAnnotator"]),
    );
    let item = feature(&report, "UTRAnnotator");
    assert_eq!(counts(item), [3, 2, 0, 2, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_reversed_value": 1}));
    assert_eq!(
        feature_residuals(&residuals, "UTRAnnotator"),
        [
            "reference\t1:100\tT\tENST1\t5UTR_consequence\tcz",
            "engine\t1:100\tT\tENST1\t5UTR_consequence\tyb",
            "reference\t1:200\tT\tENST1\t5UTR_consequence\tay",
        ]
    );
    let calls = planted.calls();
    let mut put = judged(&calls);
    put.sort();
    assert_eq!(
        put,
        [
            (ABSENT, "yb", "rs2", Some("rs2")),
            ("bx", "xb", "rs5", Some("rs5")),
            ("cz", ABSENT, "rs3", Some("rs3")),
        ]
    );
}

#[test]
fn tuple_file_states_no_rows_so_every_reference_row_is_read() {
    let scratch = Scratch::new("tuple-file-rows");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "GWAS_pmid=17"),
            row("1:100", "T", "ENST2", "GWAS_pmid=18"),
        ],
    );
    let engine = scratch.tuples(
        "engine.tsv",
        &[
            ["1:100", "T", "ENST2", "PLANT_label", "kept"],
            ["1:900", "A", "ENST9", "GWAS_pmid", "3"],
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    let item = feature(&report, "GWAS");
    assert_eq!(counts(item), [2, 1, 0, 0, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_reference_value": 2}));
    let calls = planted.calls();
    let stated: Vec<Vec<(&str, &str)>> = calls
        .iter()
        .map(|call| call.a.as_ref().unwrap().fields.iter().collect())
        .collect();
    assert_eq!(
        stated,
        [
            vec![("Location", "1:100"), ("Allele", "T"), ("Feature", "ENST1")],
            vec![
                ("Location", "1:100"),
                ("Allele", "T"),
                ("Feature", "ENST2"),
                ("PLANT_label", "kept")
            ],
        ]
    );
}
