// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `fields.json` as CPython's `json.dumps(report, indent=2)` renders it, and the stdout
//! summary lines.

use crate::merge::{FieldStats, MergeResult, PLAIN_FIELDS};
use crate::pyfmt::{push_json_str, py_float_repr, py_grouped, py_round6};
use std::fmt::Write;

pub enum Json {
    Null,
    Int(u64),
    Float(f64),
    Str(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

/// `json.dumps(value, indent=2)` with the default `ensure_ascii=True` and separators.
pub fn dumps(value: &Json) -> String {
    let mut out = String::new();
    write_value(value, 0, &mut out);
    out
}

fn indent(level: usize, out: &mut String) {
    out.push('\n');
    for _ in 0..level {
        out.push_str("  ");
    }
}

fn write_value(value: &Json, level: usize, out: &mut String) {
    match value {
        Json::Null => out.push_str("null"),
        Json::Int(n) => write!(out, "{n}").unwrap(),
        Json::Float(x) => out.push_str(&py_float_repr(*x)),
        Json::Str(s) => push_json_str(out, s),
        Json::Array(items) => {
            if items.is_empty() {
                out.push_str("[]");
                return;
            }
            out.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                indent(level + 1, out);
                write_value(item, level + 1, out);
            }
            indent(level, out);
            out.push(']');
        }
        Json::Object(items) => {
            if items.is_empty() {
                out.push_str("{}");
                return;
            }
            out.push('{');
            for (i, (key, item)) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                indent(level + 1, out);
                push_json_str(out, key);
                out.push_str(": ");
                write_value(item, level + 1, out);
            }
            indent(level, out);
            out.push('}');
        }
    }
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn agreement(stats: &FieldStats) -> Option<f64> {
    if stats.compared == 0 {
        None
    } else {
        Some(py_round6(stats.equal as f64 / stats.compared as f64))
    }
}

fn field_block(stats: &FieldStats, examples: i64) -> Json {
    let top = stats
        .top_mismatches(examples)
        .iter()
        .map(|shape| {
            Json::Object(vec![
                ("reference".into(), Json::Str(text(shape.reference))),
                ("scored".into(), Json::Str(text(shape.scored))),
                ("count".into(), Json::Int(shape.count)),
            ])
        })
        .collect();
    Json::Object(vec![
        ("compared".into(), Json::Int(stats.compared)),
        ("equal".into(), Json::Int(stats.equal)),
        ("differ".into(), Json::Int(stats.compared - stats.equal)),
        (
            "agreement".into(),
            agreement(stats).map_or(Json::Null, Json::Float),
        ),
        (
            "distinct_mismatch_shapes".into(),
            Json::Int(stats.distinct_mismatch_shapes() as u64),
        ),
        ("top_mismatches".into(), Json::Array(top)),
    ])
}

/// The report object in the Python's key order.
pub fn build(
    result: &MergeResult,
    examples: i64,
    reference_input: &str,
    scored_input: &str,
) -> Json {
    let totals = &result.totals;
    let totals_block = Json::Object(vec![
        ("reference_keys".into(), Json::Int(totals.reference_keys)),
        ("scored_keys".into(), Json::Int(totals.scored_keys)),
        ("intersection".into(), Json::Int(totals.intersection)),
        ("reference_only".into(), Json::Int(totals.reference_only)),
        ("scored_only".into(), Json::Int(totals.scored_only)),
        (
            "reference_duplicate_keys".into(),
            Json::Int(totals.reference_duplicate_keys),
        ),
        (
            "scored_duplicate_keys".into(),
            Json::Int(totals.scored_duplicate_keys),
        ),
    ]);
    let fields = PLAIN_FIELDS
        .iter()
        .zip(&result.fields)
        .map(|(name, stats)| (name.to_string(), field_block(stats, examples)))
        .collect();
    let extra = result
        .extra
        .iter()
        .map(|(key, stats)| (text(key), field_block(stats, examples)))
        .collect();
    Json::Object(vec![
        ("totals".into(), totals_block),
        ("fields".into(), Json::Object(fields)),
        ("extra".into(), Json::Object(extra)),
        (
            "inputs".into(),
            Json::Object(vec![
                ("reference".into(), Json::Str(reference_input.to_string())),
                ("scored".into(), Json::Str(scored_input.to_string())),
            ]),
        ),
        (
            "key".into(),
            Json::Array(
                [
                    "Location",
                    "Allele",
                    "Feature",
                    "Feature_type",
                    "consequence_set",
                ]
                .iter()
                .map(|k| Json::Str(k.to_string()))
                .collect(),
            ),
        ),
        (
            "comparison".into(),
            Json::Str("verbatim per field on tuples both sides emit under the same key".into()),
        ),
    ])
}

/// The stdout lines: the key totals, then one row per plain field and per `Extra:<key>`.
pub fn summary(result: &MergeResult) -> String {
    let totals = &result.totals;
    let mut out = format!(
        "keys: reference {} scored {} intersection {} (reference-only {}, scored-only {})\n",
        py_grouped(totals.reference_keys),
        py_grouped(totals.scored_keys),
        py_grouped(totals.intersection),
        py_grouped(totals.reference_only),
        py_grouped(totals.scored_only)
    );
    let plain = PLAIN_FIELDS
        .iter()
        .zip(&result.fields)
        .map(|(name, stats)| (name.to_string(), stats));
    let extra = result
        .extra
        .iter()
        .map(|(key, stats)| (format!("Extra:{}", text(key)), stats));
    for (name, stats) in plain.chain(extra) {
        let agr = agreement(stats).map_or_else(|| "None".to_string(), py_float_repr);
        writeln!(
            out,
            "  {name:<24} equal {:>12} / {:>12}  agreement {agr}",
            py_grouped(stats.equal),
            py_grouped(stats.compared)
        )
        .unwrap();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dumps_matches_json_dumps_indent_2() {
        let value = Json::Object(vec![
            ("a".into(), Json::Str("\u{e9}".into())),
            ("e".into(), Json::Object(vec![])),
            ("l".into(), Json::Array(vec![])),
            ("n".into(), Json::Null),
            ("f".into(), Json::Float(1.0)),
            (
                "list".into(),
                Json::Array(vec![
                    Json::Int(1),
                    Json::Object(vec![("k".into(), Json::Int(2))]),
                ]),
            ),
        ]);
        let expected = "{\n  \"a\": \"\\u00e9\",\n  \"e\": {},\n  \"l\": [],\n  \"n\": null,\n  \"f\": 1.0,\n  \"list\": [\n    1,\n    {\n      \"k\": 2\n    }\n  ]\n}";
        assert_eq!(dumps(&value), expected);
    }
}
