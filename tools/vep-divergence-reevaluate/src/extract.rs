// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The compared value: what each compare kind reads from one engine's run, how the values
//! are normalised, and the verdict over the two engines' values and the recorded ones.

use crate::engine::{read_lines, EngineRun, RunError};
use crate::pyfmt::{last_line, py_rstrip, py_splitlines, py_strip};
use crate::pyjson;
use crate::pypath;
use crate::pyre::Regex;
use crate::record::Reproduction;
use crate::value::{as_list, Value, EMPTY_DICT};
use std::collections::HashMap;

const DEFAULT_COLUMNS: [&str; 4] = ["Location", "Allele", "Feature", "Consequence"];
/// Any absolute path, so a value naming the host tree and one naming a container mount
/// read alike; the slash of an allele string (`-/A`, `C//G`, `*/Q`) follows `-`, `/` or
/// `*` and is not one.
pub const PATH: &str = r"(?<![\w.*/-])/(?:[\w.+@%:,=~-]+/)*[\w.+@%:,=~-]+";
/// The compare kinds that read the engine's failure itself, so a non-zero exit is their
/// subject.
const FAILURE_KINDS: [&str; 3] = ["exit", "stderr", "stdout"];

type Row = HashMap<String, String>;

fn compile(pattern: &str) -> Result<Regex, RunError> {
    Regex::new(pattern).map_err(|error| RunError(error.to_string()))
}

/// `str.partition("=")`: the key, and the value when the item holds an `=`.
fn key_value(item: &str) -> (&str, &str) {
    item.split_once('=').unwrap_or((item, ""))
}

/// The default or tab output's rows, the `Extra` column's items folded in as columns.
pub fn table_rows(lines: &[String]) -> Vec<Row> {
    let mut columns: Vec<&str> = Vec::new();
    let mut rows = Vec::new();
    for line in lines {
        if line.starts_with("##") {
            continue;
        }
        if let Some(header) = line.strip_prefix('#') {
            columns = header.split('\t').collect();
            continue;
        }
        if py_strip(line).is_empty() || columns.is_empty() {
            continue;
        }
        let mut row: Row = columns
            .iter()
            .zip(line.split('\t'))
            .map(|(column, cell)| (column.to_string(), cell.to_string()))
            .collect();
        if let Some(extra) = row.get("Extra").filter(|extra| *extra != "-").cloned() {
            for item in extra.split(';') {
                let (key, value) = key_value(item);
                row.entry(key.to_string())
                    .or_insert_with(|| value.to_string());
            }
        }
        row.insert("_line".to_string(), line.clone());
        rows.push(row);
    }
    rows
}

fn where_column(key: &str) -> Option<&'static str> {
    Some(match key {
        "record" => "Uploaded_variation",
        "location" => "Location",
        "allele" => "Allele",
        "gene" => "Gene",
        "feature" => "Feature",
        "feature_type" => "Feature_type",
        "consequence" => "Consequence",
        _ => return None,
    })
}

fn sorted_terms(text: &str) -> Vec<&str> {
    let mut terms: Vec<&str> = text.split(',').collect();
    terms.sort_unstable();
    terms
}

fn row_matches(row: &Row, where_: &Value, pattern: Option<&Regex>) -> bool {
    for (key, want) in where_.items() {
        if key == "pattern" {
            if !pattern.is_some_and(|pattern| pattern.is_match(&row["_line"])) {
                return false;
            }
            continue;
        }
        let Some(column) = where_column(key) else {
            continue;
        };
        let have = row.get(column).map_or("", String::as_str);
        let wanted = as_list(want);
        let found = if key == "consequence" {
            wanted.iter().any(|w| sorted_terms(have) == sorted_terms(w))
        } else {
            wanted.iter().any(|w| w == have)
        };
        if !found {
            return false;
        }
    }
    true
}

/// A column or `Extra` item of a row, `-` when the row lacks it.
fn cell<'a>(row: &'a Row, name: &str) -> &'a str {
    row.get(name).map_or("-", String::as_str)
}

/// The lines of one `[section]` of a statistics report, up to the next blank line or
/// section.
pub fn section_lines(lines: &[String], section: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut inside = false;
    for line in lines {
        if py_strip(line) == section {
            inside = true;
            continue;
        }
        if inside {
            if py_strip(line).is_empty() || line.starts_with('[') {
                break;
            }
            out.push(line.clone());
        }
    }
    out
}

/// The values at a dotted path: `a.b[].c` iterates a list, `a[key=val].c` keeps the
/// members whose `key` equals `val`.
pub fn json_path<'a>(value: &'a Value, path: &str) -> Vec<&'a Value> {
    let parts = Regex::new(r"[^.\[\]]+|\[[^\]]*\]")
        .expect("the path pattern")
        .find_all(path);
    let mut current = vec![value];
    for part in parts {
        let mut next = Vec::new();
        if part == "[]" {
            for v in current {
                if let Value::List(items) = v {
                    next.extend(items.iter());
                }
            }
        } else if part.starts_with('[') {
            let (key, want) = key_value(&part[1..part.len() - 1]);
            for v in current {
                if let Value::List(items) = v {
                    next.extend(items.iter().filter(|item| {
                        item.is_dict()
                            && item
                                .get(key)
                                .map_or_else(|| "None".to_string(), Value::py_str)
                                == want
                    }));
                }
            }
        } else {
            for v in current {
                if let Some(member) = v.get(&part) {
                    next.push(member);
                }
            }
        }
        current = next;
    }
    current
}

/// A JSON value as a compared string: a string as it is, any other value as compact JSON
/// with its keys sorted.
pub fn json_value(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        other => pyjson::dumps_compact_sorted(other),
    }
}

/// One data line of a VCF output.
struct VcfRecord<'a> {
    chrom: &'a str,
    pos: &'a str,
    id: &'a str,
    info: HashMap<&'a str, &'a str>,
    line: &'a str,
}

/// A VCF output's INFO formats (`CSQ` and the like, from their `##INFO` lines) and its
/// data lines as fields.
fn vcf_records(lines: &[String]) -> (HashMap<String, Vec<String>>, Vec<VcfRecord<'_>>) {
    let info_format =
        Regex::new(r#"##INFO=<ID=([^,]+),.*Format: ([^"]+)""#).expect("the INFO pattern");
    let mut formats = HashMap::new();
    let mut records = Vec::new();
    for line in lines {
        if line.starts_with("##INFO=<ID=") {
            if let Some(found) = info_format.match_start(line) {
                let id = found.group(1).expect("the INFO id");
                let names = found.group(2).expect("the INFO format");
                formats.insert(id, names.split('|').map(str::to_string).collect());
            }
            continue;
        }
        if line.starts_with('#') || py_strip(line).is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let mut info = HashMap::new();
        if let Some(column) = fields.get(7) {
            for item in column.split(';') {
                let (key, value) = key_value(item);
                info.insert(key, value);
            }
        }
        records.push(VcfRecord {
            chrom: fields[0],
            pos: fields.get(1).copied().unwrap_or(""),
            id: fields.get(2).copied().unwrap_or(""),
            info,
            line,
        });
    }
    (formats, records)
}

/// Whether a VCF or JSON record meets every `where` key the record level knows.
fn record_matches(
    id: &str,
    chrom: &str,
    pos: &str,
    line: &str,
    where_: &Value,
    pattern: Option<&Regex>,
) -> bool {
    let location = format!("{chrom}:{pos}");
    for (key, want) in where_.items() {
        let have = match key.as_str() {
            "pattern" => {
                if !pattern.is_some_and(|pattern| pattern.is_match(line)) {
                    return false;
                }
                continue;
            }
            "record" | "id" => id,
            "chrom" => chrom,
            "pos" => pos,
            "location" => &location,
            _ => continue,
        };
        if !as_list(want).iter().any(|w| w == have) {
            return false;
        }
    }
    true
}

/// `str(obj.get(key, ""))` of a JSON record.
fn json_field(record: &Value, key: &str) -> String {
    record.get(key).map_or_else(String::new, Value::py_str)
}

/// The compared value read from one engine's run, as a list of strings.
pub fn extract(repro: &Reproduction, run: &EngineRun) -> Result<Vec<String>, RunError> {
    let c = repro.compare();
    let kind = c.get_str("kind").unwrap_or("");
    let where_ = c.get("where").unwrap_or(&EMPTY_DICT);
    let pattern = || compile(c.get_str("pattern").unwrap_or(""));
    let matching = |lines: Vec<String>, regex: &Regex| -> Vec<String> {
        lines
            .into_iter()
            .filter(|line| regex.is_match(line))
            .collect()
    };
    match kind {
        "exit" => return Ok(vec![run.exit_code.to_string()]),
        "stderr" | "stdout" => {
            let text = if kind == "stderr" {
                &run.stderr
            } else {
                &run.stdout
            };
            let lines = py_splitlines(text)
                .into_iter()
                .map(str::to_string)
                .collect();
            return Ok(matching(lines, &pattern()?));
        }
        "warnings" => return Ok(matching(read_lines(&run.warnings())?, &pattern()?)),
        "stats" => {
            return Ok(section_lines(
                &read_lines(&run.stats())?,
                c.get_str("section").unwrap_or(""),
            ))
        }
        _ => {}
    }
    let lines = read_lines(&run.output)?;
    let where_pattern = match where_.get("pattern") {
        Some(pattern) => Some(compile(&pattern.py_str())?),
        None => None,
    };
    let where_pattern = where_pattern.as_ref();
    match kind {
        "header" => {
            let regex = pattern()?;
            Ok(lines
                .into_iter()
                .filter(|line| line.starts_with('#') && regex.is_match(line))
                .collect())
        }
        "lines" => Ok(matching(lines, &pattern()?)),
        "json" => {
            let mut values = Vec::new();
            for line in &lines {
                if py_strip(line).is_empty() {
                    continue;
                }
                let Some(record) = pyjson::loads(line) else {
                    continue;
                };
                if !record.is_dict() {
                    return Err(RunError(format!(
                        "'{}' object has no attribute 'get'",
                        record.py_type()
                    )));
                }
                let id = json_field(&record, "id");
                let chrom = json_field(&record, "seq_region_name");
                let pos = json_field(&record, "start");
                if !record_matches(&id, &chrom, &pos, line, where_, where_pattern) {
                    continue;
                }
                if let Some(input) = where_.get("input") {
                    let have = record.get_str("input");
                    if !as_list(input).iter().any(|w| Some(w.as_str()) == have) {
                        continue;
                    }
                }
                let path = c.get_str("path").unwrap_or("");
                values.extend(json_path(&record, path).into_iter().map(json_value));
            }
            Ok(values)
        }
        "vcf" => {
            let (formats, records) = vcf_records(&lines);
            let mut values = Vec::new();
            for record in &records {
                if !record_matches(
                    record.id,
                    record.chrom,
                    record.pos,
                    record.line,
                    where_,
                    where_pattern,
                ) {
                    continue;
                }
                if c.truthy_at("line") {
                    values.push(record.line.to_string());
                    continue;
                }
                let key = c.get_str("info").unwrap_or("CSQ");
                let field = c.get_str("field").unwrap_or("");
                if field == key {
                    values.push(record.info.get(key).copied().unwrap_or("-").to_string());
                    continue;
                }
                let (Some(names), Some(entries)) = (formats.get(key), record.info.get(key)) else {
                    values.push("-".to_string());
                    continue;
                };
                for entry in entries.split(',') {
                    let entry_row: HashMap<&str, &str> = names
                        .iter()
                        .map(String::as_str)
                        .zip(entry.split('|'))
                        .collect();
                    let narrowed_out = |where_key: &str, column: &str| {
                        where_.get(where_key).is_some_and(|want| {
                            let have = entry_row.get(column);
                            !as_list(want).iter().any(|w| Some(&w.as_str()) == have)
                        })
                    };
                    if narrowed_out("feature", "Feature") || narrowed_out("allele", "Allele") {
                        continue;
                    }
                    values.push(entry_row.get(field).copied().unwrap_or("-").to_string());
                }
            }
            Ok(values)
        }
        _ => {
            let mut matched = Vec::new();
            for row in table_rows(&lines) {
                if row_matches(&row, where_, where_pattern) {
                    matched.push(row);
                }
            }
            if kind == "count" {
                return Ok(vec![matched.len().to_string()]);
            }
            if kind == "field" {
                let field = c.get_str("field").unwrap_or("");
                return Ok(matched
                    .iter()
                    .map(|row| {
                        if c.truthy_at("with_key") {
                            format!("{}={}", cell(row, "Feature"), cell(row, field))
                        } else {
                            cell(row, field).to_string()
                        }
                    })
                    .collect());
            }
            let names: Vec<String> = match c.get("columns") {
                Some(columns) => as_list(columns),
                None => DEFAULT_COLUMNS.iter().map(|n| n.to_string()).collect(),
            };
            Ok(matched
                .iter()
                .map(|row| {
                    names
                        .iter()
                        .map(|name| cell(row, name))
                        .collect::<Vec<_>>()
                        .join(" | ")
                })
                .collect())
        }
    }
}

/// The values as they are compared: trailing whitespace removed unless `strip = false`,
/// absolute paths replaced by `<path>` unless `mask_paths = false`, each `mask` pattern
/// replaced by `<masked>`, and the list sorted unless `sorted = false`.
pub fn normalise(values: &[String], compare: &Value) -> Result<Vec<String>, RunError> {
    let strip = compare.truthy_or_true("strip");
    let path = if compare.truthy_or_true("mask_paths") {
        Some(compile(PATH)?)
    } else {
        None
    };
    let mut masks = Vec::new();
    if let Some(Value::List(patterns)) = compare.get("mask") {
        for pattern in patterns {
            masks.push(compile(&pattern.py_str())?);
        }
    }
    let mut out = Vec::with_capacity(values.len());
    for value in values {
        let mut value = if strip {
            py_rstrip(value).to_string()
        } else {
            value.clone()
        };
        if let Some(path) = &path {
            value = path.sub("<path>", &value);
        }
        for mask in &masks {
            value = mask.sub("<masked>", &value);
        }
        out.push(value);
    }
    if compare.truthy_or_true("sorted") {
        out.sort();
    }
    Ok(out)
}

/// `closed` when the engines' values are equal, `still-present` when each equals its
/// recorded value, `changed` otherwise.
pub fn verdict(
    expected: &Value,
    compare: &Value,
    ensembl: &[String],
    vep_rs: &[String],
) -> Result<&'static str, RunError> {
    let want = |key: &str| normalise(&as_list(expected.section(key)), compare);
    let want_ensembl = want("ensembl_vep")?;
    let want_vep_rs = want("vep_rs")?;
    Ok(if ensembl == vep_rs {
        "closed"
    } else if ensembl == want_ensembl && vep_rs == want_vep_rs {
        "still-present"
    } else {
        "changed"
    })
}

/// Why an engine's run cannot be compared: an exit other than the one recorded for it, or
/// no output to read. `expected.<engine>_exit` records an exit status the entry's
/// divergence produces; without it, a kind that reads the output takes a non-zero exit or
/// a missing output file as a run that failed, never as the recorded value.
pub fn run_failure(repro: &Reproduction, engine: &str, run: &EngineRun) -> Option<String> {
    let recorded = repro
        .expected()
        .get(&format!("{}_exit", engine.replace('-', "_")));
    let kind = repro.compare().get_str("kind").unwrap_or("");
    let failure_kind = FAILURE_KINDS.contains(&kind);
    match recorded {
        Some(recorded) => {
            if recorded.as_i128() != Some(i128::from(run.exit_code)) {
                return Some(format!(
                    "{engine} exited {}, the entry records {}: {}",
                    run.exit_code,
                    recorded.py_str(),
                    last_line(&run.stderr)
                ));
            }
        }
        None if failure_kind => return None,
        None => {
            if run.exit_code != 0 {
                return Some(format!(
                    "{engine} exited {}: {}",
                    run.exit_code,
                    last_line(&run.stderr)
                ));
            }
        }
    }
    let reads_output = !failure_kind && kind != "warnings" && kind != "stats";
    let zero_or_none = recorded.is_none_or(|r| r.as_i128() == Some(0));
    if reads_output && zero_or_none && !pypath::exists(&run.output) {
        return Some(format!("{engine} wrote no output file"));
    }
    None
}

#[cfg(test)]
#[path = "extract_tests.rs"]
mod tests;
