// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The compare mode's report files: `summary.json` as `json.dumps(summary, indent=2)` writes
//! it and `summary.md`.

use crate::pyfmt::{push_json_str, py_float_repr};
use std::time::{SystemTime, UNIX_EPOCH};

/// A JSON value in the order its keys were inserted.
pub enum Json {
    Null,
    Bool(bool),
    Int(i128),
    Float(f64),
    Str(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    pub fn str(s: &str) -> Json {
        Json::Str(s.to_string())
    }

    pub fn uint(n: u64) -> Json {
        Json::Int(i128::from(n))
    }

    /// `Counter.most_common` as JSON: a list of two-element lists.
    pub fn pairs(pairs: &[(String, u64)]) -> Json {
        Json::Array(
            pairs
                .iter()
                .map(|(k, n)| Json::Array(vec![Json::str(k), Json::uint(*n)]))
                .collect(),
        )
    }

    pub fn strings(items: &[String]) -> Json {
        Json::Array(items.iter().map(|s| Json::str(s)).collect())
    }

    pub fn object(entries: Vec<(&str, Json)>) -> Json {
        Json::Object(
            entries
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect(),
        )
    }

    /// `json.dumps(value, indent=2)`: no trailing newline.
    pub fn dumps(&self) -> String {
        let mut out = String::new();
        self.write(&mut out, 0);
        out
    }

    fn write(&self, out: &mut String, level: usize) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Int(n) => out.push_str(&n.to_string()),
            Json::Float(x) => out.push_str(&py_float_repr(*x)),
            Json::Str(s) => push_json_str(out, s),
            Json::Array(items) if items.is_empty() => out.push_str("[]"),
            Json::Object(items) if items.is_empty() => out.push_str("{}"),
            Json::Array(items) => {
                out.push('[');
                for (i, item) in items.iter().enumerate() {
                    out.push_str(if i == 0 { "\n" } else { ",\n" });
                    out.push_str(&"  ".repeat(level + 1));
                    item.write(out, level + 1);
                }
                out.push('\n');
                out.push_str(&"  ".repeat(level));
                out.push(']');
            }
            Json::Object(items) => {
                out.push('{');
                for (i, (key, value)) in items.iter().enumerate() {
                    out.push_str(if i == 0 { "\n" } else { ",\n" });
                    out.push_str(&"  ".repeat(level + 1));
                    push_json_str(out, key);
                    out.push_str(": ");
                    value.write(out, level + 1);
                }
                out.push('\n');
                out.push_str(&"  ".repeat(level));
                out.push('}');
            }
        }
    }
}

/// `datetime.now(timezone.utc).isoformat()`: microseconds shown unless zero, `+00:00`.
pub fn utc_now_isoformat() -> String {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default();
    let seconds = now.as_secs() as i64;
    let micros = now.subsec_micros();
    let days = seconds.div_euclid(86_400);
    let secs_of_day = seconds.rem_euclid(86_400);
    // Civil date from days since 1970-01-01 (proleptic Gregorian).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = if month <= 2 { year + 1 } else { year };
    let fraction = if micros == 0 {
        String::new()
    } else {
        format!(".{micros:06}")
    };
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}{fraction}+00:00",
        secs_of_day / 3600,
        secs_of_day % 3600 / 60,
        secs_of_day % 60
    )
}

/// `f"{x:.6f}"`.
pub fn fixed6(x: f64) -> String {
    format!("{x:.6}")
}

/// Precision, recall and F1 from the intersection and the two sides' tuple counts.
pub fn prf(intersection: u64, rust_count: u64, perl_count: u64) -> (f64, f64, f64) {
    let precision = if rust_count > 0 {
        intersection as f64 / rust_count as f64
    } else {
        0.0
    };
    let recall = if perl_count > 0 {
        intersection as f64 / perl_count as f64
    } else {
        0.0
    };
    if precision + recall == 0.0 {
        return (precision, recall, 0.0);
    }
    (
        precision,
        recall,
        2.0 * precision * recall / (precision + recall),
    )
}

/// One paired file's metrics, the `FileMetrics` record of the report.
pub struct FileMetrics {
    pub file_name: String,
    pub perl_tuple_count: u64,
    pub rust_tuple_count: u64,
    pub intersection_count: u64,
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
    pub missing_in_rust: u64,
    pub extra_in_rust: u64,
    pub top_missing_consequence_sets: Vec<(String, u64)>,
    pub top_extra_consequence_sets: Vec<(String, u64)>,
    pub context_collapsed_precision: f64,
    pub context_collapsed_recall: f64,
    pub context_collapsed_f1: f64,
    pub context_collapsed_missing_in_rust: u64,
    pub context_collapsed_extra_in_rust: u64,
}

impl FileMetrics {
    pub fn json(&self) -> Json {
        Json::object(vec![
            ("file_name", Json::str(&self.file_name)),
            ("perl_tuple_count", Json::uint(self.perl_tuple_count)),
            ("rust_tuple_count", Json::uint(self.rust_tuple_count)),
            ("intersection_count", Json::uint(self.intersection_count)),
            ("precision", Json::Float(self.precision)),
            ("recall", Json::Float(self.recall)),
            ("f1", Json::Float(self.f1)),
            ("missing_in_rust", Json::uint(self.missing_in_rust)),
            ("extra_in_rust", Json::uint(self.extra_in_rust)),
            (
                "top_missing_consequence_sets",
                Json::pairs(&self.top_missing_consequence_sets),
            ),
            (
                "top_extra_consequence_sets",
                Json::pairs(&self.top_extra_consequence_sets),
            ),
            (
                "context_collapsed_precision",
                Json::Float(self.context_collapsed_precision),
            ),
            (
                "context_collapsed_recall",
                Json::Float(self.context_collapsed_recall),
            ),
            (
                "context_collapsed_f1",
                Json::Float(self.context_collapsed_f1),
            ),
            (
                "context_collapsed_missing_in_rust",
                Json::uint(self.context_collapsed_missing_in_rust),
            ),
            (
                "context_collapsed_extra_in_rust",
                Json::uint(self.context_collapsed_extra_in_rust),
            ),
        ])
    }
}

/// A side's tuple totals with the derived precision, recall and F1.
pub struct Scores {
    pub perl_tuple_count: u64,
    pub rust_tuple_count: u64,
    pub intersection_count: u64,
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
}

impl Scores {
    pub fn new(intersection: u64, rust: u64, perl: u64) -> Scores {
        let (precision, recall, f1) = prf(intersection, rust, perl);
        Scores {
            perl_tuple_count: perl,
            rust_tuple_count: rust,
            intersection_count: intersection,
            precision,
            recall,
            f1,
        }
    }

    fn entries(&self) -> Vec<(&'static str, Json)> {
        vec![
            ("perl_tuple_count", Json::uint(self.perl_tuple_count)),
            ("rust_tuple_count", Json::uint(self.rust_tuple_count)),
            ("intersection_count", Json::uint(self.intersection_count)),
            ("precision", Json::Float(self.precision)),
            ("recall", Json::Float(self.recall)),
            ("f1", Json::Float(self.f1)),
        ]
    }
}

/// The adjusted block of the aggregate.
pub struct Adjusted<'a> {
    pub scores: Scores,
    pub excluded_perl_tuples: u64,
    pub excluded_rust_tuples: u64,
    pub excluded_categories: &'a [(String, crate::rules::filter::BucketStats)],
    /// The reference release and its class definitions; absent under the published 115.2 shape.
    pub reference: Option<(&'a str, &'a [crate::rules::table::DivergenceClass])>,
}

/// The aggregate block of the report.
pub struct Aggregate<'a> {
    pub raw: Scores,
    pub top_missing_consequence_sets: Vec<(String, u64)>,
    pub top_extra_consequence_sets: Vec<(String, u64)>,
    pub context_collapsed: Scores,
    pub adjusted: Adjusted<'a>,
}

impl Aggregate<'_> {
    pub fn json(&self) -> Json {
        let mut entries = self.raw.entries();
        entries.push((
            "top_missing_consequence_sets",
            Json::pairs(&self.top_missing_consequence_sets),
        ));
        entries.push((
            "top_extra_consequence_sets",
            Json::pairs(&self.top_extra_consequence_sets),
        ));
        let mut collapsed = self.context_collapsed.entries();
        collapsed.push((
            "collapsed_terms",
            Json::Array(
                sorted_context_terms()
                    .iter()
                    .map(|t| Json::str(t))
                    .collect(),
            ),
        ));
        entries.push(("context_collapsed", Json::object(collapsed)));
        let adjusted = &self.adjusted;
        let mut adj = adjusted.scores.entries();
        adj.push((
            "excluded_perl_tuples",
            Json::uint(adjusted.excluded_perl_tuples),
        ));
        adj.push((
            "excluded_rust_tuples",
            Json::uint(adjusted.excluded_rust_tuples),
        ));
        adj.push((
            "excluded_categories",
            Json::Object(
                adjusted
                    .excluded_categories
                    .iter()
                    .map(|(bucket, stats)| {
                        (
                            bucket.clone(),
                            Json::Object(
                                stats
                                    .iter()
                                    .map(|(k, v)| (k.to_string(), Json::uint(*v)))
                                    .collect(),
                            ),
                        )
                    })
                    .collect(),
            ),
        ));
        if let Some((release, classes)) = adjusted.reference {
            adj.push(("reference_release", Json::str(release)));
            adj.push((
                "divergence_classes",
                Json::Object(
                    classes
                        .iter()
                        .map(|c| {
                            (
                                c.bucket.to_string(),
                                Json::object(vec![
                                    ("definition", Json::str(c.definition)),
                                    ("excludes", Json::Bool(c.excludes)),
                                    ("taxonomy_class", Json::str(c.taxonomy_class)),
                                    ("perl_citation", Json::str(c.perl_citation)),
                                ]),
                            )
                        })
                        .collect(),
                ),
            ));
        }
        entries.push(("adjusted", Json::object(adj)));
        Json::object(entries)
    }
}

/// `sorted(TRANSCRIPT_CONTEXT_TERMS)`.
pub fn sorted_context_terms() -> Vec<&'static str> {
    let mut terms = crate::input::TRANSCRIPT_CONTEXT_TERMS.to_vec();
    terms.sort_unstable();
    terms
}

/// The text of `summary.md`, a trailing newline included.
pub fn summary_markdown(aggregate: &Aggregate, files: &[FileMetrics]) -> String {
    let mut lines: Vec<String> = vec![
        "# Concordance Summary".to_string(),
        String::new(),
        "## Aggregate".to_string(),
        String::new(),
        format!("- Perl tuples: {}", aggregate.raw.perl_tuple_count),
        format!("- Rust tuples: {}", aggregate.raw.rust_tuple_count),
        format!("- Intersection: {}", aggregate.raw.intersection_count),
        format!("- Precision: {}", fixed6(aggregate.raw.precision)),
        format!("- Recall: {}", fixed6(aggregate.raw.recall)),
        format!("- F1: {}", fixed6(aggregate.raw.f1)),
    ];
    let adjusted = &aggregate.adjusted;
    if let Some((release, _)) = adjusted.reference {
        lines.push(format!("- Reference release: Ensembl VEP {release}"));
    }
    let count_of = |stats: &crate::rules::filter::BucketStats, name: &str| -> u64 {
        stats
            .iter()
            .find(|(k, _)| *k == name)
            .map_or(0, |(_, v)| *v)
    };
    let mut buckets: Vec<&(String, crate::rules::filter::BucketStats)> =
        adjusted.excluded_categories.iter().collect();
    buckets.sort_by_key(|(_, stats)| std::cmp::Reverse(count_of(stats, "count")));
    let cats: Vec<String> = buckets
        .iter()
        .filter(|(_, stats)| count_of(stats, "count") > 0)
        .map(|(bucket, stats)| {
            format!(
                "{bucket}: {} ({} excluded)",
                count_of(stats, "count"),
                count_of(stats, "excluded")
            )
        })
        .collect();
    if !cats.is_empty() {
        lines.push(format!("- Divergence buckets: {}", cats.join(", ")));
    }
    let not_evaluated: Vec<String> = adjusted
        .excluded_categories
        .iter()
        .map(|(bucket, stats)| (bucket, count_of(stats, crate::rules::filter::NOT_EVALUATED)))
        .filter(|(_, rows)| *rows > 0)
        .map(|(bucket, rows)| format!("{bucket}: {rows} row{}", if rows == 1 { "" } else { "s" }))
        .collect();
    if !not_evaluated.is_empty() {
        lines.push(format!(
            "- Not evaluated (the run gives no sequence names of its cache): {}",
            not_evaluated.join(", ")
        ));
    }
    if adjusted.excluded_perl_tuples > 0 {
        lines.push(format!(
            "- **Adjusted F1: {}** (excluded {}P + {}R intended-divergence tuples)",
            fixed6(adjusted.scores.f1),
            adjusted.excluded_perl_tuples,
            adjusted.excluded_rust_tuples
        ));
    }
    lines.push(format!(
        "- Context-collapsed F1: {} (collapsed terms: {})",
        fixed6(aggregate.context_collapsed.f1),
        sorted_context_terms().join(", ")
    ));
    if let Some((_, classes)) = adjusted.reference {
        if !classes.is_empty() {
            lines.push(String::new());
            lines.push("## Divergence classes".to_string());
            lines.push(String::new());
            lines.push("| Class | Treatment | Definition | Ensembl VEP source |".to_string());
            lines.push("|---|---|---|---|".to_string());
            for c in classes {
                let verdict = if c.excludes {
                    "sets aside"
                } else {
                    "counts only"
                };
                lines.push(format!(
                    "| `{}` | {verdict} | {} | {} |",
                    c.bucket, c.definition, c.perl_citation
                ));
            }
        }
    }
    lines.push(String::new());
    lines.push("## Per File".to_string());
    lines.push(String::new());
    lines.push(
        "| File | Perl Tuples | Rust Tuples | Intersection | Precision | Recall | F1 |".to_string(),
    );
    lines.push("|---|---:|---:|---:|---:|---:|---:|".to_string());
    for item in files {
        lines.push(format!(
            "| {} | {} | {} | {} | {} | {} | {} |",
            item.file_name,
            item.perl_tuple_count,
            item.rust_tuple_count,
            item.intersection_count,
            fixed6(item.precision),
            fixed6(item.recall),
            fixed6(item.f1)
        ));
    }
    lines.push(String::new());
    lines.push("## Per File (Context-Collapsed)".to_string());
    lines.push(String::new());
    lines.push("| File | Precision | Recall | F1 | Missing | Extra |".to_string());
    lines.push("|---|---:|---:|---:|---:|---:|".to_string());
    for item in files {
        lines.push(format!(
            "| {} | {} | {} | {} | {} | {} |",
            item.file_name,
            fixed6(item.context_collapsed_precision),
            fixed6(item.context_collapsed_recall),
            fixed6(item.context_collapsed_f1),
            item.context_collapsed_missing_in_rust,
            item.context_collapsed_extra_in_rust
        ));
    }
    lines.push(String::new());
    lines.join("\n")
}

#[cfg(test)]
#[path = "report_ratio_tests.rs"]
mod ratio_tests;
