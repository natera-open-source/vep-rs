// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison's printed report and its JSON report.

use super::fields::Run;
use super::tally::{Mismatch, PluginResult, Stats};
use super::text::{fixed_two, json_float, push_json_str, push_str_repr, PyInt};
use std::fmt::Write;

/// One file pair's results: its name and the result of each distinct plugin of the run.
pub struct FileResults {
    pub name: String,
    pub plugins: Vec<PluginResult>,
}

/// A plugin's counts summed over the file pairs.
#[derive(Default)]
struct Totals {
    matched: u64,
    mismatch: u64,
    perl_only: u64,
    rust_only: u64,
    annotated: u64,
    perl_annotated: u64,
    rust_annotated: u64,
}

impl Totals {
    fn of(results: &[FileResults], slot: usize) -> Totals {
        let mut totals = Totals::default();
        for file in results {
            totals.add(&file.plugins[slot].stats);
        }
        totals
    }

    fn add(&mut self, stats: &Stats) {
        self.matched += stats.matched;
        self.mismatch += stats.mismatch;
        self.perl_only += stats.perl_only;
        self.rust_only += stats.rust_only;
        self.annotated += stats.annotated_keys;
        self.perl_annotated += stats.perl_annotated;
        self.rust_annotated += stats.rust_annotated;
    }

    fn concordance_pct(&self) -> Option<f64> {
        (self.annotated > 0).then(|| self.matched as f64 / self.annotated as f64 * 100.0)
    }
}

/// Whether a count is below the floor of annotated records.
fn below(count: u64, floor: &PyInt) -> bool {
    (count as i128) < floor.saturating_i128()
}

/// Appends a side's fields as CPython prints a `dict` of `str`.
fn push_dict_repr(out: &mut String, data: &[(String, String)]) {
    out.push('{');
    for (at, (field, value)) in data.iter().enumerate() {
        if at > 0 {
            out.push_str(", ");
        }
        push_str_repr(out, field);
        out.push_str(": ");
        push_str_repr(out, value);
    }
    out.push('}');
}

/// The printed report: per requested plugin its section with each file pair's counts and
/// first mismatches, the counts summed over the pairs and the plugin's status. The second
/// value tells whether any plugin fails.
///
/// A plugin fails where either engine annotated fewer than `min_annotated` records, which
/// measures nothing, and on any mismatch or one-sided record.
pub fn printed(results: &[FileResults], run: &Run, min_annotated: &PyInt) -> (String, bool) {
    let mut out = String::new();
    let mut any_failure = false;
    let rule = "=".repeat(60);
    for (name, slot) in run.plugins.iter().zip(&run.slots) {
        let _ = writeln!(out, "\n{rule}\n  Plugin: {name}\n{rule}");
        for file in results {
            let result = &file.plugins[*slot];
            let stats = &result.stats;
            let _ = writeln!(out, "\n  File: {}", file.name);
            let shown = &result.fields[..result.fields.len().min(5)];
            let more = if result.fields.len() > 5 { "..." } else { "" };
            let _ = writeln!(out, "    Fields: {}{more}", shown.join(", "));
            if !result.shared_fields.is_empty() {
                let _ = writeln!(
                    out,
                    "    Also entries of another requested plugin: {}",
                    result.shared_fields.join(", ")
                );
            }
            let _ = writeln!(out, "    Perl records: {}", result.perl_annotated_records);
            let _ = writeln!(out, "    Rust records: {}", result.rust_annotated_records);
            let _ = writeln!(
                out,
                "    Match: {}  Mismatch: {}  Perl-only: {}  Rust-only: {}",
                stats.matched, stats.mismatch, stats.perl_only, stats.rust_only
            );
            match stats.concordance_pct {
                None => {
                    out.push_str("    Concordance: n/a (no records annotated by either engine)\n")
                }
                Some(pct) => {
                    let _ = writeln!(out, "    Concordance: {}%", fixed_two(pct));
                }
            }
            if !result.mismatches_sample.is_empty() {
                out.push_str("    First mismatches:\n");
                for mismatch in result.mismatches_sample.iter().take(3) {
                    let [location, allele, feature] = &mismatch.key;
                    let _ = writeln!(
                        out,
                        "      {location} {allele} {feature}: {}",
                        mismatch.difference.name()
                    );
                    for (side, data) in [("Perl", &mismatch.perl), ("Rust", &mismatch.rust)] {
                        if !data.is_empty() {
                            let _ = write!(out, "        {side}: ");
                            push_dict_repr(&mut out, data);
                            out.push('\n');
                        }
                    }
                }
            }
        }
        let totals = Totals::of(results, *slot);
        let aggregate = match totals.concordance_pct() {
            Some(pct) => format!("{}%", fixed_two(pct)),
            None => "n/a".to_string(),
        };
        let _ = writeln!(
            out,
            "\n  AGGREGATE: {}/{} match ({aggregate})",
            totals.matched, totals.annotated
        );
        let _ = writeln!(
            out,
            "  ANNOTATED: perl={}  rust={}  (floor {min_annotated})",
            totals.perl_annotated, totals.rust_annotated
        );
        let mut reasons = Vec::new();
        if below(totals.perl_annotated, min_annotated) {
            reasons.push(format!(
                "Perl annotated {} records (< {min_annotated})",
                totals.perl_annotated
            ));
        }
        if below(totals.rust_annotated, min_annotated) {
            reasons.push(format!(
                "vep-rs annotated {} records (< {min_annotated})",
                totals.rust_annotated
            ));
        }
        if totals.mismatch > 0 {
            reasons.push(format!("{} value mismatches", totals.mismatch));
        }
        if totals.perl_only > 0 {
            reasons.push(format!("{} Perl-only records", totals.perl_only));
        }
        if totals.rust_only > 0 {
            reasons.push(format!("{} vep-rs-only records", totals.rust_only));
        }
        if reasons.is_empty() {
            out.push_str("  STATUS: PASS\n");
        } else {
            let _ = writeln!(out, "  STATUS: FAIL ({})", reasons.join("; "));
            any_failure = true;
        }
    }
    (out, any_failure)
}

/// A JSON value, written as CPython's `json.dumps(value, indent=2)` writes it.
enum Json {
    Null,
    Bool(bool),
    /// An integer's decimal text.
    Int(String),
    Float(f64),
    Text(String),
    List(Vec<Json>),
    /// Members in their order; a later member of a name already present takes its place.
    Object(Vec<(String, Json)>),
}

impl Json {
    fn count(value: u64) -> Json {
        Json::Int(value.to_string())
    }

    fn texts(values: &[String]) -> Json {
        Json::List(values.iter().cloned().map(Json::Text).collect())
    }

    fn object(members: Vec<(&str, Json)>) -> Json {
        Json::Object(
            members
                .into_iter()
                .map(|(name, value)| (name.to_string(), value))
                .collect(),
        )
    }

    /// Sets a member as assignment to a `dict` key does.
    fn set(members: &mut Vec<(String, Json)>, name: &str, value: Json) {
        match members.iter_mut().find(|(known, _)| known == name) {
            Some(member) => member.1 = value,
            None => members.push((name.to_string(), value)),
        }
    }

    fn write(&self, out: &mut String, depth: usize) {
        let indent = |out: &mut String, depth: usize| {
            out.push('\n');
            out.extend(std::iter::repeat_n("  ", depth));
        };
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(true) => out.push_str("true"),
            Json::Bool(false) => out.push_str("false"),
            Json::Int(text) => out.push_str(text),
            Json::Float(value) => out.push_str(&json_float(*value)),
            Json::Text(text) => push_json_str(out, text),
            Json::List(items) if items.is_empty() => out.push_str("[]"),
            Json::List(items) => {
                out.push('[');
                for (at, item) in items.iter().enumerate() {
                    if at > 0 {
                        out.push(',');
                    }
                    indent(out, depth + 1);
                    item.write(out, depth + 1);
                }
                indent(out, depth);
                out.push(']');
            }
            Json::Object(members) if members.is_empty() => out.push_str("{}"),
            Json::Object(members) => {
                out.push('{');
                for (at, (name, value)) in members.iter().enumerate() {
                    if at > 0 {
                        out.push(',');
                    }
                    indent(out, depth + 1);
                    push_json_str(out, name);
                    out.push_str(": ");
                    value.write(out, depth + 1);
                }
                indent(out, depth);
                out.push('}');
            }
        }
    }
}

fn side(data: &[(String, String)]) -> Json {
    Json::Object(
        data.iter()
            .map(|(field, value)| (field.clone(), Json::Text(value.clone())))
            .collect(),
    )
}

fn mismatch_json(mismatch: &Mismatch) -> Json {
    Json::object(vec![
        ("key", Json::texts(&mismatch.key)),
        ("type", Json::Text(mismatch.difference.name().to_string())),
        ("perl", side(&mismatch.perl)),
        ("rust", side(&mismatch.rust)),
    ])
}

fn result_json(result: &PluginResult) -> Json {
    let stats = &result.stats;
    Json::object(vec![
        ("fields", Json::texts(&result.fields)),
        ("shared_fields", Json::texts(&result.shared_fields)),
        (
            "perl_annotated_records",
            Json::count(result.perl_annotated_records),
        ),
        (
            "rust_annotated_records",
            Json::count(result.rust_annotated_records),
        ),
        (
            "stats",
            Json::object(vec![
                ("total_keys", Json::count(stats.total_keys)),
                ("match", Json::count(stats.matched)),
                ("mismatch", Json::count(stats.mismatch)),
                ("perl_only", Json::count(stats.perl_only)),
                ("rust_only", Json::count(stats.rust_only)),
                ("both_empty", Json::count(stats.both_empty)),
                (
                    "concordance_pct",
                    stats.concordance_pct.map_or(Json::Null, Json::Float),
                ),
                ("annotated_keys", Json::count(stats.annotated_keys)),
                ("perl_annotated", Json::count(stats.perl_annotated)),
                ("rust_annotated", Json::count(stats.rust_annotated)),
            ]),
        ),
        (
            "mismatches_sample",
            Json::List(result.mismatches_sample.iter().map(mismatch_json).collect()),
        ),
        ("total_mismatches", Json::count(result.total_mismatches)),
    ])
}

/// The JSON report: the run's plugins, tolerance and floor, the number of file pairs, every
/// pair's results and, per plugin, the counts summed over the pairs with whether it passes.
pub fn json(results: &[FileResults], run: &Run, tolerance: f64, min_annotated: &PyInt) -> String {
    let mut files: Vec<(String, Json)> = Vec::new();
    for file in results {
        let mut plugins: Vec<(String, Json)> = Vec::new();
        for (plugin, result) in run.distinct.iter().zip(&file.plugins) {
            plugins.push((plugin.name.clone(), result_json(result)));
        }
        Json::set(&mut files, &file.name, Json::Object(plugins));
    }
    let mut aggregates: Vec<(String, Json)> = Vec::new();
    for (slot, plugin) in run.distinct.iter().enumerate() {
        let totals = Totals::of(results, slot);
        let pass = totals.mismatch == 0
            && totals.perl_only == 0
            && totals.rust_only == 0
            && !below(totals.perl_annotated, min_annotated)
            && !below(totals.rust_annotated, min_annotated);
        aggregates.push((
            plugin.name.clone(),
            Json::object(vec![
                ("match", Json::count(totals.matched)),
                ("mismatch", Json::count(totals.mismatch)),
                ("perl_only", Json::count(totals.perl_only)),
                ("rust_only", Json::count(totals.rust_only)),
                ("annotated", Json::count(totals.annotated)),
                ("perl_annotated", Json::count(totals.perl_annotated)),
                ("rust_annotated", Json::count(totals.rust_annotated)),
                (
                    "concordance_pct",
                    totals.concordance_pct().map_or(Json::Null, Json::Float),
                ),
                ("pass", Json::Bool(pass)),
            ]),
        ));
    }
    let report = Json::object(vec![
        ("plugins", Json::texts(&run.plugins)),
        ("tolerance", Json::Float(tolerance)),
        ("min_annotated", Json::Int(min_annotated.to_string())),
        ("file_count", Json::count(results.len() as u64)),
        ("results", Json::Object(files)),
        ("aggregates", Json::Object(aggregates)),
    ]);
    let mut out = String::new();
    report.write(&mut out, 0);
    out
}
