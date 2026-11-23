// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The five gates, scored in order: the consequence F1 of the concordance summary against
//! its floor, then four relative deltas between the baseline's metrics and the builder's
//! against their ceilings. Each gate adds its line as it is scored, so a report that stops
//! a later gate leaves the lines of the gates before it.

use crate::json::{object, string, truthy, Object, Value};
use crate::number::{equals_zero, number, relative_delta, Number};
use crate::pyfmt::{py_fixed6, py_float_repr, py_round6, py_str};
use std::collections::BTreeSet;

/// The four reports of a run; no concordance summary where the operand is empty or names
/// nothing.
pub struct Reports {
    pub baseline: Value,
    pub builder: Value,
    pub cache_diff: Value,
    pub concordance: Option<Value>,
}

/// The five thresholds, in the order of the gates.
pub struct Thresholds {
    pub f1: f64,
    pub total: f64,
    pub colocated: f64,
    pub frequency: f64,
    pub consequence: f64,
}

/// What the gates come to: the verdict report and whether every gate passed.
pub struct Verdict {
    pub report: Value,
    pub passed: bool,
}

fn word(passed: bool) -> &'static str {
    if passed {
        "PASS"
    } else {
        "FAIL"
    }
}

/// A report that must be an object for the gate that reads it.
fn as_object<'a>(value: &'a Value, name: &str) -> Result<&'a Object, String> {
    match value {
        Value::Object(object) => Ok(object),
        _ => Err(format!("{name} is not an object")),
    }
}

/// The `aggregate` block of the concordance summary, as `report and "aggregate" in report`
/// finds it: none in a summary that is absent, empty or without the key, and a summary of
/// another type that holds the word stops the run where Python would index it.
fn aggregate(concordance: Option<&Value>) -> Result<Option<&Value>, String> {
    let name = "the concordance report";
    match concordance {
        Some(report) if truthy(report) => match report {
            Value::Object(object) => Ok(object.get(b"aggregate")),
            Value::List(items) => {
                let named =
                    |item: &Value| matches!(item, Value::Str(text) if **text == *b"aggregate");
                if items.iter().any(named) {
                    return Err(format!("{name} is an array, not an object"));
                }
                Ok(None)
            }
            Value::Str(text) => {
                if text.windows(9).any(|window| window == b"aggregate") {
                    return Err(format!("{name} is a string, not an object"));
                }
                Ok(None)
            }
            _ => Err(format!("{name} is a number or true, not an object")),
        },
        _ => Ok(None),
    }
}

/// Gate 1: the F1 of the concordance summary's aggregate is at or above its floor.
fn consequence_f1(
    concordance: Option<&Value>,
    gate: f64,
    lines: &mut String,
) -> Result<(Value, bool), String> {
    let name = string("consequence_f1");
    let Some(aggregate) = aggregate(concordance)? else {
        lines.push_str("  Consequence F1: N/A (concordance report missing) FAIL\n");
        let scored = object(vec![
            ("gate", name),
            ("threshold", Value::Float(gate)),
            ("actual", Value::Null),
            ("passed", Value::Bool(false)),
            ("note", string("concordance report not available")),
        ]);
        return Ok((scored, false));
    };
    let aggregate = as_object(aggregate, "aggregate of the concordance report")?;
    let f1 = match aggregate.get(b"f1") {
        Some(value) => {
            number(value).map_err(|why| format!("f1 of the concordance report {why}"))?
        }
        None => Number::Float(0.0),
    };
    let passed = f1.at_least(gate);
    // `round()` of an integer is that integer, and of a float the float nearest to six
    // places.
    let actual = match f1 {
        Number::Int(n) => Value::Int(n),
        Number::Float(x) => Value::Float(py_round6(x)),
    };
    lines.push_str(&format!(
        "  Consequence F1: {} (gate >= {}) {}\n",
        py_fixed6(f1.to_float()),
        py_float_repr(gate),
        word(passed)
    ));
    let scored = object(vec![
        ("gate", name),
        ("threshold", Value::Float(gate)),
        ("actual", actual),
        ("passed", Value::Bool(passed)),
    ]);
    Ok((scored, passed))
}

/// The delta of one count: `abs(builder - baseline) / baseline` over a baseline above
/// zero; under any other baseline 0.0 where the builder's count equals zero and 1.0 where
/// it does not, whatever the builder's value is.
fn delta_of(baseline: &Value, builder: &Value) -> Result<f64, String> {
    let b = number(baseline).map_err(|why| format!("the baseline's count {why}"))?;
    if b.positive() {
        let c = number(builder).map_err(|why| format!("the builder's count {why}"))?;
        relative_delta(c, b).map_err(|why| format!("the builder's count {why}"))
    } else if equals_zero(builder) {
        Ok(0.0)
    } else {
        Ok(1.0)
    }
}

/// One of the two gates over a single count: its key in the metrics, its name in the
/// report and its label on the line.
struct CountGate {
    key: &'static str,
    name: &'static str,
    label: &'static str,
}

/// Gates 2 and 3: the relative delta of one count is below its ceiling.
fn count_delta(
    which: &CountGate,
    baseline: &Object,
    builder: &Object,
    gate: f64,
    lines: &mut String,
) -> Result<(Value, bool), String> {
    let zero = Value::Int(0);
    let b = baseline.get(which.key.as_bytes()).unwrap_or(&zero);
    let c = builder.get(which.key.as_bytes()).unwrap_or(&zero);
    let delta = delta_of(b, c).map_err(|why| format!("{}: {why}", which.key))?;
    let passed = delta < gate;
    let shown = |value: &Value| {
        py_str(value).ok_or_else(|| {
            format!(
                "{}: a count holds a surrogate, which no line can carry",
                which.key
            )
        })
    };
    lines.push_str(&format!(
        "  {}: {} ({} vs {}) (gate < {}) {}\n",
        which.label,
        py_fixed6(delta),
        shown(b)?,
        shown(c)?,
        py_float_repr(gate),
        word(passed)
    ));
    let scored = object(vec![
        ("gate", string(which.name)),
        ("threshold", Value::Float(gate)),
        ("baseline", b.clone()),
        ("builder", c.clone()),
        ("actual", Value::Float(py_round6(delta))),
        ("passed", Value::Bool(passed)),
    ]);
    Ok((scored, passed))
}

/// One side's block of counts by name: an object's entries, nothing where the metrics
/// hold no such block, and `Inert` for an empty array or an empty string, which names no
/// count and cannot be asked for one.
enum Side<'a> {
    Counts(&'a Object),
    Absent,
    Inert,
}

impl<'a> Side<'a> {
    fn of(metrics: &'a Object, key: &str, whose: &str) -> Result<Side<'a>, String> {
        match metrics.get(key.as_bytes()) {
            None => Ok(Side::Absent),
            Some(Value::Object(object)) => Ok(Side::Counts(object)),
            Some(Value::List(items)) if items.is_empty() => Ok(Side::Inert),
            Some(Value::Str(text)) if text.is_empty() => Ok(Side::Inert),
            Some(_) => Err(format!("{key} of the {whose} metrics is not an object")),
        }
    }

    fn names(&self) -> impl Iterator<Item = &'a [u8]> {
        let object = match self {
            Side::Counts(object) => Some(*object),
            _ => None,
        };
        object.into_iter().flat_map(Object::keys)
    }

    fn get(&self, name: &[u8]) -> Option<&'a Value> {
        match self {
            Side::Counts(object) => object.get(name),
            _ => None,
        }
    }
}

/// One of the two gates over a block of counts by name.
struct BlockGate {
    key: &'static str,
    name: &'static str,
    label: &'static str,
}

/// Gates 4 and 5: over every name either side counts, the largest relative delta is
/// below its ceiling. The details stand in the code point order of the names.
fn largest_delta(
    which: &BlockGate,
    baseline: &Object,
    builder: &Object,
    gate: f64,
    lines: &mut String,
) -> Result<(Value, bool), String> {
    let b_side = Side::of(baseline, which.key, "baseline")?;
    let c_side = Side::of(builder, which.key, "builder")?;
    let names: BTreeSet<&[u8]> = b_side.names().chain(c_side.names()).collect();
    if !names.is_empty() {
        for (side, whose) in [(&b_side, "baseline"), (&c_side, "builder")] {
            if matches!(side, Side::Inert) {
                return Err(format!(
                    "{} of the {whose} metrics is not an object",
                    which.key
                ));
            }
        }
    }
    let zero = Value::Int(0);
    let mut largest = 0.0;
    let mut details = Object::default();
    for name in names {
        let b = b_side.get(name).unwrap_or(&zero);
        let c = c_side.get(name).unwrap_or(&zero);
        let delta = delta_of(b, c)
            .map_err(|why| format!("{}, {}: {why}", which.key, String::from_utf8_lossy(name)))?;
        details.insert(
            name,
            object(vec![
                ("baseline", b.clone()),
                ("builder", c.clone()),
                ("delta", Value::Float(py_round6(delta))),
            ]),
        );
        // `max(largest, delta)`: the delta takes the place only where it is above.
        if delta > largest {
            largest = delta;
        }
    }
    let passed = largest < gate;
    lines.push_str(&format!(
        "  {}: {} (gate < {}) {}\n",
        which.label,
        py_fixed6(largest),
        py_float_repr(gate),
        word(passed)
    ));
    let scored = object(vec![
        ("gate", string(which.name)),
        ("threshold", Value::Float(gate)),
        ("actual", Value::Float(py_round6(largest))),
        ("details", Value::Object(details)),
        ("passed", Value::Bool(passed)),
    ]);
    Ok((scored, passed))
}

/// The three values the verdict carries over from the cache comparison report, `null`
/// for a section or a key the report does not hold.
fn cache_diff_summary(cache_diff: &Value) -> Result<Value, String> {
    let report = as_object(cache_diff, "the cache comparison report")?;
    let carried = |section: &str, key: &str| -> Result<Value, String> {
        let Some(block) = report.get(section.as_bytes()) else {
            return Ok(Value::Null);
        };
        let block = as_object(block, &format!("{section} of the cache comparison report"))?;
        Ok(block.get(key.as_bytes()).cloned().unwrap_or(Value::Null))
    };
    Ok(object(vec![
        (
            "transcript_match_rate",
            carried("transcripts", "match_rate")?,
        ),
        ("variation_match_rate", carried("variations", "match_rate")?),
        ("info_json_matches", carried("info_json", "matches")?),
    ]))
}

/// Scores the five gates over `reports`, a line a gate onto `lines`, and builds the
/// verdict report.
pub fn score(
    reports: &Reports,
    thresholds: &Thresholds,
    lines: &mut String,
) -> Result<Verdict, String> {
    let mut gates = Vec::new();
    let mut all_pass = true;
    let mut keep = |(scored, passed): (Value, bool)| {
        gates.push(scored);
        all_pass &= passed;
    };

    keep(consequence_f1(
        reports.concordance.as_ref(),
        thresholds.f1,
        lines,
    )?);

    let baseline = as_object(&reports.baseline, "the baseline metrics")?;
    let builder = as_object(&reports.builder, "the builder metrics")?;
    let total = CountGate {
        key: "total_annotation_lines",
        name: "total_annotations_delta",
        label: "Total annotations delta",
    };
    keep(count_delta(
        &total,
        baseline,
        builder,
        thresholds.total,
        lines,
    )?);
    let colocated = CountGate {
        key: "existing_variation_count",
        name: "colocated_variant_delta",
        label: "Co-located variant delta",
    };
    keep(count_delta(
        &colocated,
        baseline,
        builder,
        thresholds.colocated,
        lines,
    )?);
    let frequency = BlockGate {
        key: "frequency_field_counts",
        name: "frequency_field_delta",
        label: "Frequency field max delta",
    };
    keep(largest_delta(
        &frequency,
        baseline,
        builder,
        thresholds.frequency,
        lines,
    )?);
    let consequence = BlockGate {
        key: "consequence_counts",
        name: "per_consequence_type_delta",
        label: "Per-consequence max delta",
    };
    keep(largest_delta(
        &consequence,
        baseline,
        builder,
        thresholds.consequence,
        lines,
    )?);

    let summary = cache_diff_summary(&reports.cache_diff)?;
    let report = object(vec![
        ("verdict", string(word(all_pass))),
        ("gates", Value::List(gates)),
        ("cache_diff_summary", summary),
    ]);
    Ok(Verdict {
        report,
        passed: all_pass,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::json::load;
    use crate::pyfmt::push_json;

    fn thresholds() -> Thresholds {
        Thresholds {
            f1: 0.999,
            total: 0.001,
            colocated: 0.01,
            frequency: 0.01,
            consequence: 0.01,
        }
    }

    fn reports(
        baseline: &str,
        builder: &str,
        cache_diff: &str,
        concordance: Option<&str>,
    ) -> Reports {
        let read = |text: &str| load(text.as_bytes()).unwrap();
        Reports {
            baseline: read(baseline),
            builder: read(builder),
            cache_diff: read(cache_diff),
            concordance: concordance.map(read),
        }
    }

    /// The lines of a run and, where every gate was scored, the verdict word.
    fn run(reports: &Reports) -> (String, Result<&'static str, String>) {
        let mut lines = String::new();
        let outcome = score(reports, &thresholds(), &mut lines).map(|verdict| word(verdict.passed));
        (lines, outcome)
    }

    #[test]
    fn equal_metrics_pass_every_gate() {
        let metrics = r#"{"total_annotation_lines": 10, "existing_variation_count": 4,
            "frequency_field_counts": {"AF": 3}, "consequence_counts": {"intron_variant": 10}}"#;
        let (lines, outcome) = run(&reports(
            metrics,
            metrics,
            "{}",
            Some(r#"{"aggregate": {"f1": 1.0}}"#),
        ));
        assert_eq!(outcome, Ok("PASS"));
        assert_eq!(
            lines,
            "  Consequence F1: 1.000000 (gate >= 0.999) PASS\n\
             \x20 Total annotations delta: 0.000000 (10 vs 10) (gate < 0.001) PASS\n\
             \x20 Co-located variant delta: 0.000000 (4 vs 4) (gate < 0.01) PASS\n\
             \x20 Frequency field max delta: 0.000000 (gate < 0.01) PASS\n\
             \x20 Per-consequence max delta: 0.000000 (gate < 0.01) PASS\n"
        );
    }

    #[test]
    fn a_missing_summary_fails_the_first_gate_alone() {
        for concordance in [
            None,
            Some("{}"),
            Some("null"),
            Some("[\"f1\"]"),
            Some("\"f1\""),
        ] {
            let (lines, outcome) = run(&reports("{}", "{}", "{}", concordance));
            assert_eq!(outcome, Ok("FAIL"), "{concordance:?}");
            assert_eq!(
                lines,
                "  Consequence F1: N/A (concordance report missing) FAIL\n\
                 \x20 Total annotations delta: 0.000000 (0 vs 0) (gate < 0.001) PASS\n\
                 \x20 Co-located variant delta: 0.000000 (0 vs 0) (gate < 0.01) PASS\n\
                 \x20 Frequency field max delta: 0.000000 (gate < 0.01) PASS\n\
                 \x20 Per-consequence max delta: 0.000000 (gate < 0.01) PASS\n"
            );
        }
    }

    #[test]
    fn a_report_of_the_wrong_shape_stops_at_the_gate_that_reads_it() {
        let lines_before = |reports: &Reports| {
            let (lines, outcome) = run(reports);
            assert!(outcome.is_err());
            lines.lines().count()
        };
        let summary = Some(r#"{"aggregate": {"f1": 1.0}}"#);
        // The concordance summary: before any line.
        for concordance in [
            "[\"aggregate\"]",
            "\"an aggregate\"",
            "7",
            "true",
            r#"{"aggregate": []}"#,
            r#"{"aggregate": {"f1": null}}"#,
        ] {
            assert_eq!(
                lines_before(&reports("{}", "{}", "{}", Some(concordance))),
                0
            );
        }
        // The metrics as a whole and the two single counts: after the first line or the
        // second.
        assert_eq!(lines_before(&reports("[]", "{}", "{}", summary)), 1);
        assert_eq!(lines_before(&reports("{}", "null", "{}", summary)), 1);
        let count = |b: &str, c: &str| {
            reports(
                &format!(r#"{{"existing_variation_count": {b}}}"#),
                &format!(r#"{{"existing_variation_count": {c}}}"#),
                "{}",
                summary,
            )
        };
        assert_eq!(lines_before(&count("null", "1")), 2);
        assert_eq!(lines_before(&count("\"1\"", "1")), 2);
        assert_eq!(lines_before(&count("1", "null")), 2);
        assert_eq!(lines_before(&count("1", "[]")), 2);
        assert_eq!(lines_before(&count("0", "\"a\\ud800\"")), 2);
        // The blocks of counts by name: after the third line or the fourth.
        let block = |key: &str, b: &str, c: &str| {
            reports(
                &format!(r#"{{"{key}": {b}}}"#),
                &format!(r#"{{"{key}": {c}}}"#),
                "{}",
                summary,
            )
        };
        for (b, c) in [
            ("7", "{}"),
            ("{}", "null"),
            ("[1]", "{}"),
            ("\"a\"", "{}"),
            ("[]", "{\"AF\": 1}"),
            ("{\"AF\": 1}", "\"\""),
            ("{\"AF\": null}", "{}"),
            ("{\"AF\": 1}", "{\"AF\": \"1\"}"),
        ] {
            assert_eq!(
                lines_before(&block("frequency_field_counts", b, c)),
                3,
                "{b} {c}"
            );
            assert_eq!(
                lines_before(&block("consequence_counts", b, c)),
                4,
                "{b} {c}"
            );
        }
        // The cache comparison report: after all five lines.
        for cache_diff in [
            "[]",
            "null",
            r#"{"transcripts": []}"#,
            r#"{"info_json": 1}"#,
        ] {
            assert_eq!(lines_before(&reports("{}", "{}", cache_diff, summary)), 5);
        }
    }

    #[test]
    fn a_count_of_any_type_scores_beside_a_baseline_of_zero() {
        let scored = |b: &str, c: &str| {
            let reports = reports(
                &format!(
                    r#"{{"total_annotation_lines": {b}, "frequency_field_counts": {{"AF": {b}}}}}"#
                ),
                &format!(
                    r#"{{"total_annotation_lines": {c}, "frequency_field_counts": {{"AF": {c}}}}}"#
                ),
                "{}",
                None,
            );
            let (lines, outcome) = run(&reports);
            assert_eq!(outcome, Ok("FAIL"));
            let wanted: Vec<&str> = lines.lines().collect();
            format!("{} | {}", wanted[1].trim(), wanted[3].trim())
        };
        assert_eq!(
            scored("0", "null"),
            "Total annotations delta: 1.000000 (0 vs None) (gate < 0.001) FAIL | Frequency field max delta: 1.000000 (gate < 0.01) FAIL"
        );
        assert_eq!(
            scored("-3", "[\"a\", {\"b\": 1.0}]"),
            "Total annotations delta: 1.000000 (-3 vs ['a', {'b': 1.0}]) (gate < 0.001) FAIL | Frequency field max delta: 1.000000 (gate < 0.01) FAIL"
        );
        assert_eq!(
            scored("false", "\"text\""),
            "Total annotations delta: 1.000000 (False vs text) (gate < 0.001) FAIL | Frequency field max delta: 1.000000 (gate < 0.01) FAIL"
        );
        assert_eq!(
            scored("0.0", "false"),
            "Total annotations delta: 0.000000 (0.0 vs False) (gate < 0.001) PASS | Frequency field max delta: 0.000000 (gate < 0.01) PASS"
        );
        assert_eq!(
            scored("NaN", "170141183460469231731687303715884105728"),
            "Total annotations delta: 1.000000 (nan vs 170141183460469231731687303715884105728) (gate < 0.001) FAIL | Frequency field max delta: 1.000000 (gate < 0.01) FAIL"
        );
        assert_eq!(
            scored("true", "2"),
            "Total annotations delta: 1.000000 (True vs 2) (gate < 0.001) FAIL | Frequency field max delta: 1.000000 (gate < 0.01) FAIL"
        );
    }

    #[test]
    fn details_stand_in_code_point_order_and_empty_blocks_name_nothing() {
        let (_, outcome) = run(&reports(
            r#"{"frequency_field_counts": []}"#,
            r#"{"frequency_field_counts": ""}"#,
            "{}",
            None,
        ));
        assert_eq!(outcome, Ok("FAIL"));
        let reports = reports(
            r#"{"consequence_counts": {"b": 2, "a": 4, "\u00e9": 1, "B": 1}}"#,
            r#"{"consequence_counts": {"c": 1, "a": 3, "b": 2, "\ud83d\ude00": 0}}"#,
            r#"{"transcripts": {"match_rate": 0.5}, "info_json": {"matches": false}}"#,
            None,
        );
        let mut lines = String::new();
        let verdict = score(&reports, &thresholds(), &mut lines).unwrap();
        let mut text = String::new();
        push_json(&mut text, &verdict.report, 0);
        let keys: Vec<&str> = text
            .lines()
            .filter(|line| line.starts_with("        \""))
            .collect();
        assert_eq!(
            keys,
            [
                "        \"B\": {",
                "        \"a\": {",
                "        \"b\": {",
                "        \"c\": {",
                "        \"\\u00e9\": {",
                "        \"\\ud83d\\ude00\": {",
            ]
        );
        assert!(text.contains("\"actual\": 1.0,\n      \"details\": {"));
        assert!(text.ends_with(
            "  \"cache_diff_summary\": {\n    \"transcript_match_rate\": 0.5,\n    \"variation_match_rate\": null,\n    \"info_json_matches\": false\n  }\n}"
        ));
    }
}
