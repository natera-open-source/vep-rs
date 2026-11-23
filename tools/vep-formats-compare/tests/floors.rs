// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Floors on plugin columns in the JSON format, read from the committed plugin corpora.

use indexmap::IndexMap;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use vep_formats_compare::compare_files::compare_files;
use vep_formats_compare::floors::{json_floor_location, StrListOrStr};
use vep_formats_compare::input::{input_records, read_text};
use vep_formats_compare::lists::json_object_member;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::Side;

/// Rows of the AlphaMissense corpus's references carrying the plugin's fields: every
/// missense row under the default options and under cols=all, the rows whose transcript
/// version the plugin matched under transcript_match=1.
const ALPHAMISSENSE_ROWS: i64 = 245;
const ALPHAMISSENSE_MATCHED_ROWS: i64 = 16;
/// Rows of the OpenTargets corpus's reference carrying the plugin's fields.
const OPENTARGETS_ROWS: i64 = 383;

fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// A reference compared with itself in `fmt` under `floors`.
fn compare_pair(
    fmt: &str,
    path: &Path,
    floors: &[(&str, i64)],
    inputs: Option<&[String]>,
) -> Value {
    let floors: IndexMap<String, i64> = floors
        .iter()
        .map(|(column, minimum)| (column.to_string(), *minimum))
        .collect();
    compare_files(
        fmt,
        path,
        path,
        inputs,
        &mut Documented::new(),
        &BTreeSet::new(),
        &floors,
        &IndexMap::new(),
        4,
        false,
        None,
        None,
        None,
        &mut Side::default(),
    )
    .unwrap()
}

fn text<'a>(report: &'a Value, key: &str) -> &'a str {
    report.get(key).and_then(Value::as_str).unwrap()
}

fn texts(report: &Value, key: &str) -> Vec<String> {
    report
        .get(key)
        .and_then(Value::as_list)
        .unwrap()
        .iter()
        .map(|item| item.as_str().unwrap().to_string())
        .collect()
}

fn column<'a>(report: &'a Value, name: &str) -> &'a Value {
    report
        .get("columns")
        .and_then(|columns| columns.get(name))
        .unwrap()
}

fn fields(items: &[&str]) -> Option<StrListOrStr> {
    Some(StrListOrStr::StrList(
        items.iter().map(|item| item.to_string()).collect(),
    ))
}

#[test]
fn a_floor_on_a_plugin_column_is_counted_inside_the_plugin_object() {
    let alphamissense = corpus("GRCh38-plugin-alphamissense");
    let inputs = input_records(&alphamissense.join("variants.vcf")).unwrap();
    let cases: [(&str, i64, &[&str]); 3] = [
        (
            "expected",
            ALPHAMISSENSE_ROWS,
            &["am_class", "am_pathogenicity"],
        ),
        (
            "expected_transcript-match",
            ALPHAMISSENSE_MATCHED_ROWS,
            &["am_class"],
        ),
        (
            "expected_cols-all",
            ALPHAMISSENSE_ROWS,
            &["am_genome", "am_transcript_id"],
        ),
    ];
    for (expected_dir, rows, columns) in cases {
        let floors: Vec<(&str, i64)> = columns.iter().map(|name| (*name, rows)).collect();
        let json_path = alphamissense.join(expected_dir).join("json.jsonl.gz");
        let report = compare_pair("json", &json_path, &floors, Some(inputs.as_slice()));
        assert_eq!(text(&report, "verdict"), "green", "{expected_dir}");
        assert!(texts(&report, "floors_not_applicable").is_empty());
        for name in columns {
            let counted = column(&report, name);
            assert_eq!(
                texts(counted, "floor_fields"),
                [format!("alphamissense.{name}")]
            );
            assert_eq!(counted.get("floor_met"), Some(&Value::Bool(true)));
            assert_eq!(counted.get("expected_populated"), Some(&Value::Int(rows)));
            assert_eq!(counted.get("actual_populated"), Some(&Value::Int(rows)));
        }
        let default = compare_pair(
            "default",
            &alphamissense.join(expected_dir).join("default.txt.gz"),
            &floors,
            None,
        );
        assert_eq!(text(&default, "verdict"), "green");
        for name in columns {
            assert_eq!(
                column(&default, name).get("expected_populated"),
                Some(&Value::Int(rows))
            );
        }
    }
    let short = compare_pair(
        "json",
        &alphamissense.join("expected").join("json.jsonl.gz"),
        &[("am_class", ALPHAMISSENSE_ROWS + 1)],
        Some(inputs.as_slice()),
    );
    assert_eq!(text(&short, "verdict"), "red");
    assert_eq!(
        texts(&short, "reasons"),
        [format!(
            "floor am_class={}: expected populated {ALPHAMISSENSE_ROWS}, actual populated {ALPHAMISSENSE_ROWS}",
            ALPHAMISSENSE_ROWS + 1
        )]
    );
}

#[test]
fn a_core_column_stays_out_of_a_plugin_object_and_every_object_counts() {
    let spliceai_corpus = corpus("GRCh38-plugin-spliceai-fixture");
    let inputs = input_records(&spliceai_corpus.join("variants.vcf")).unwrap();
    let json_path = spliceai_corpus.join("expected").join("json.jsonl.gz");
    let spliceai = parse_json(&read_text(&json_path).unwrap(), &inputs).unwrap();
    assert!(!spliceai.iter().any(|e| e.fields.contains_key("symbol")));
    assert!(spliceai.iter().any(|e| {
        json_object_member(e, "spliceai", "SYMBOL").is_some_and(|symbol| !symbol.is_empty())
    }));
    assert_eq!(
        json_floor_location("SYMBOL", &spliceai, &spliceai),
        fields(&["gene_symbol"])
    );
    assert_eq!(
        json_floor_location("DS_AG", &spliceai, &spliceai),
        fields(&["spliceai.ds_ag"])
    );
    let json_report = compare_pair(
        "json",
        &json_path,
        &[("SYMBOL", 1)],
        Some(inputs.as_slice()),
    );
    let default = compare_pair(
        "default",
        &spliceai_corpus.join("expected").join("default.txt.gz"),
        &[("SYMBOL", 1)],
        None,
    );
    for report in [&json_report, &default] {
        assert_eq!(text(report, "verdict"), "red");
        assert_eq!(
            texts(report, "reasons"),
            ["floor SYMBOL=1: expected populated 0, actual populated 0"]
        );
    }
    let opentargets = corpus("GRCh38-plugin-opentargets");
    let inputs = input_records(&opentargets.join("variants.vcf")).unwrap();
    let mut rows = 0;
    for expected_dir in ["expected", "expected_all-columns-every-record"] {
        let default = compare_pair(
            "default",
            &opentargets.join(expected_dir).join("default.txt.gz"),
            &[("OpenTargets_gwasGeneId", OPENTARGETS_ROWS)],
            None,
        );
        rows = match column(&default, "OpenTargets_gwasGeneId").get("expected_populated") {
            Some(Value::Int(rows)) => *rows,
            other => panic!("no count of populated rows: {other:?}"),
        };
        let report = compare_pair(
            "json",
            &opentargets.join(expected_dir).join("json.jsonl.gz"),
            &[("OpenTargets_gwasGeneId", rows)],
            Some(inputs.as_slice()),
        );
        assert_eq!(text(&report, "verdict"), "green", "{expected_dir}");
        let counted = column(&report, "OpenTargets_gwasGeneId");
        assert_eq!(
            texts(counted, "floor_fields"),
            ["opentargets.opentargets_gwasgeneid"]
        );
        assert_eq!(counted.get("expected_populated"), Some(&Value::Int(rows)));
        assert_eq!(counted.get("actual_populated"), Some(&Value::Int(rows)));
    }
    assert!(rows >= OPENTARGETS_ROWS);
}
