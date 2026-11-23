// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The records of a run, numbered as the golden tests number them, the names an output
//! finds them under, the records the reference carried with no row and the records it
//! skipped.

use super::outputs::TupleSets;
use crate::constants::DEFAULT_INPUT;
use crate::error::{failure, Failure};
use crate::json::{loads, JsonError, Object, Value};
use crate::names::record_names;
use crate::pyfmt::py_strip;
use crate::pytext::{is_file, path_join, py_repr_str, py_str};
use crate::rows::{not_readable_text, opener, ByteSource, ReadError, TextLines};
use regex::Regex;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::sync::OnceLock;

/// One record of a run, at its ordinal as the golden tests number records.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RunRecord {
    /// The record's line.
    pub vcf: String,
    /// Its place among the manifest's `records`; `None` for a line of a flag set's own input
    /// that twins no manifest record.
    pub manifest_index: Option<usize>,
    /// The reason Ensembl VEP gave for skipping it, where the manifest records one.
    pub reference_warning: Option<String>,
}

/// A VCF line's CHROM, POS, ID, REF and ALT, the columns a twin of the line (another
/// encoding, sample columns added or dropped) shares with it.
pub fn first_five_columns(line: &str) -> &str {
    match line.match_indices('\t').nth(4) {
        Some((at, _)) => &line[..at],
        None => line,
    }
}

/// The manifest's `records` as run records; a manifest without the list, or with a record
/// that is no object with a `vcf` line, is refused in the tool's own words.
fn manifest_records(manifest: &Object) -> Result<Vec<RunRecord>, Failure> {
    let held = manifest
        .get("records")
        .and_then(Value::as_list)
        .ok_or_else(|| Failure::own("the manifest holds no list of records"))?;
    held.iter()
        .enumerate()
        .map(|(index, record)| {
            let vcf = record.get("vcf").and_then(Value::as_str).ok_or_else(|| {
                Failure::own(format!("manifest record {index} holds no vcf line"))
            })?;
            let warning = record.get("reference_warning").and_then(Value::as_str);
            Ok(RunRecord {
                vcf: vcf.to_string(),
                manifest_index: Some(index),
                reference_warning: warning.map(str::to_string),
            })
        })
        .collect()
}

/// Whether a line is a VCF record: five tab-separated columns or more, a POS of decimal
/// digits (`str.isdecimal()`) and a REF of bases or `.`, never the `/`-joined allele string
/// or the class word of an Ensembl-format line.
fn vcf_record(line: &str) -> bool {
    static DECIMAL: OnceLock<Regex> = OnceLock::new();
    let decimal = DECIMAL.get_or_init(|| Regex::new(r"^\p{Nd}+$").expect("a pattern"));
    let columns: Vec<&str> = line.split('\t').collect();
    columns.len() >= 5
        && decimal.is_match(columns[1])
        && !columns[3].is_empty()
        && columns[3]
            .bytes()
            .all(|base| b"ACGTNacgtn.".contains(&base))
}

/// The text of a region line before its allele: `chr:start-end` and the strand where the
/// line writes one. `None` for a line that is no region line: one that holds a blank
/// (Utils/VEP.pm:403 reads a region in a line of one column) or that Ensembl VEP's pattern
/// of a region does not fit (Utils/VEP.pm:391, matched without regard to case).
fn region_position(line: &str) -> Option<&str> {
    static REGION: OnceLock<Regex> = OnceLock::new();
    // Both cases are spelled out: under `(?i)` a letter class also takes the letters outside
    // ASCII that fold into it, which Ensembl VEP's match of bytes refuses.
    let region = REGION.get_or_init(|| {
        Regex::new(r"^[^:]+:[0-9]+-[0-9]+(?::[-+]?1)?[/:]([A-Za-z0-9:]{3,}|[ACGTNacgtn-]+)$")
            .expect("a pattern")
    });
    if line.contains(char::is_whitespace) {
        return None;
    }
    let allele = region.captures(line)?.get(1)?;
    // The allele follows a `/` or a `:` of one byte.
    Some(&line[..allele.start() - 1])
}

/// The region lines among the records of a run, by ordinal: a line of a flag set's own
/// input that holds no VCF record and is a region line (`region_position`) names its
/// record in both engines' outputs, since Ensembl VEP gives a region record its line as
/// its name (Parser/Region.pm:198 and :213) and writes that name in the first column
/// (OutputFactory.pm:895-901). `None` at the ordinal of any other line, whose record its
/// VCF line names; no ordinal for a run over `variants.vcf` or over a file that holds a
/// VCF record. A file that cannot be read is refused with its path.
pub fn region_lines(corpus: &str, target: &Object) -> Result<Vec<Option<String>>, Failure> {
    let input_name = match target.get("input").filter(|input| input.truthy()) {
        Some(Value::Str(name)) if name != DEFAULT_INPUT => name.as_str(),
        _ => return Ok(Vec::new()),
    };
    let path = path_join(corpus, input_name);
    let unreadable = |error: ReadError| not_readable_text(&path, &error);
    let mut reader = opener(&path).map_err(unreadable)?;
    let mut lines: Vec<Option<String>> = Vec::new();
    while reader.advance().map_err(unreadable)? {
        let line = reader.line();
        if line.starts_with('#') {
            continue;
        }
        if vcf_record(line) {
            return Ok(Vec::new());
        }
        lines.push(region_position(line).map(|_| line.to_string()));
    }
    Ok(lines)
}

/// The names the records of a run carry in the first column of an output.
#[derive(Debug)]
pub struct RunNames {
    /// Uploaded_variation to the ordinals of the records it names, in order.
    by_name: HashMap<String, Vec<usize>>,
    /// The region line at each ordinal that has one (`region_lines`).
    region_lines: Vec<Option<String>>,
}

impl RunNames {
    /// The names of `records`: a record a region line names (`region_lines`) carries that
    /// line alone, any other the names of its VCF line (`record_names`).
    pub fn of(
        records: &[RunRecord],
        region_lines: Vec<Option<String>>,
    ) -> Result<RunNames, Failure> {
        let mut by_name: HashMap<String, Vec<usize>> = HashMap::new();
        for (ordinal, record) in records.iter().enumerate() {
            let names = match region_lines.get(ordinal) {
                Some(Some(line)) => vec![line.clone()],
                _ => record_names(&record.vcf)?,
            };
            for name in names {
                by_name.entry(name).or_default().push(ordinal);
            }
        }
        Ok(RunNames {
            by_name,
            region_lines,
        })
    }

    /// The ordinals of the records `name` names, in order.
    pub fn get(&self, name: &str) -> Option<&Vec<usize>> {
        self.by_name.get(name)
    }

    /// The region line that names the record at `ordinal`, where one does.
    fn region_line(&self, ordinal: usize) -> Option<&str> {
        self.region_lines.get(ordinal).and_then(Option::as_deref)
    }
}

/// The records of a run, numbered as the golden harness numbers them (every line of the
/// run's input that is not a header, in order): the manifest's `records` for a run over
/// `variants.vcf`; for a flag set over another `input`, that file's lines. A line sharing
/// its first five columns with a manifest record is that record (a twin of `variants.vcf`
/// in another encoding keeps its `reference_warning` and takes `reference_rows`). A file
/// holding a VCF record (`vcf_record`) is a VCF whose other lines are records of their own,
/// which carry no `reference_warning`, so a record the reference skipped in such an input
/// is refused as unexplained. A file holding none is `variants.vcf` re-encoded in another
/// format or delimiter (`--format hgvs`, an Ensembl-format line, a space-delimited VCF
/// line), record for record, and each line is the manifest record at its ordinal; such a
/// file whose record count is not the manifest's is refused with its path, as is a line of
/// a VCF input that is neither a twin nor a VCF record, and a file that cannot be read.
pub fn run_records(
    corpus: &str,
    manifest: &Object,
    target: &Object,
) -> Result<Vec<RunRecord>, Failure> {
    let input_name = match target.get("input").filter(|input| input.truthy()) {
        None => DEFAULT_INPUT,
        Some(Value::Str(name)) => name.as_str(),
        Some(_) => return Err(Failure::own("the run's input is not a file name")),
    };
    let held = manifest_records(manifest)?;
    if input_name == DEFAULT_INPUT {
        return Ok(held);
    }
    // The record of the manifest that carries a head, the last where several do.
    let by_head: HashMap<&str, usize> = held
        .iter()
        .enumerate()
        .map(|(index, record)| (first_five_columns(&record.vcf), index))
        .collect();
    let path = path_join(corpus, input_name);
    let lines = input_lines(&path)?;
    if lines.iter().any(|line| vcf_record(line)) {
        let mut records = Vec::with_capacity(lines.len());
        for (ordinal, line) in lines.into_iter().enumerate() {
            if let Some(index) = by_head.get(first_five_columns(&line)) {
                records.push(held[*index].clone());
            } else if vcf_record(&line) {
                records.push(RunRecord {
                    vcf: line,
                    manifest_index: None,
                    reference_warning: None,
                });
            } else {
                let opening: String = line.chars().take(40).collect();
                return Err(failure(format!(
                    "{path} record {ordinal} ({}) is no VCF record in a VCF input",
                    py_repr_str(&opening)
                )));
            }
        }
        return Ok(records);
    }
    if lines.len() != held.len() {
        return Err(failure(format!(
            "{path} holds no VCF record and its {} records are not the manifest's {} re-encoded",
            lines.len(),
            held.len()
        )));
    }
    Ok(held)
}

/// The lines of an input file that are no header, in order; a file that cannot be read is
/// refused with its path.
fn input_lines(path: &str) -> Result<Vec<String>, Failure> {
    let unreadable = |error: ReadError| not_readable_text(path, &error);
    let mut reader = opener(path).map_err(unreadable)?;
    let mut lines: Vec<String> = Vec::new();
    while reader.advance().map_err(unreadable)? {
        if !reader.line().starts_with('#') {
            lines.push(reader.line().to_string());
        }
    }
    Ok(lines)
}

/// The lines of the input of a flag set that is `variants.vcf` re-encoded in another format
/// or delimiter, each the record at its ordinal (`run_records`): an output of such a run
/// names a record by its own line, where it names one at all. `None` for a run over
/// `variants.vcf` and for an input that holds a VCF record.
pub fn reencoded_lines(corpus: &str, target: &Object) -> Result<Option<Vec<String>>, Failure> {
    let input_name = match target.get("input").filter(|input| input.truthy()) {
        Some(Value::Str(name)) if name != DEFAULT_INPUT => name.as_str(),
        _ => return Ok(None),
    };
    let lines = input_lines(&path_join(corpus, input_name))?;
    let reencoded = !lines.iter().any(|line| vcf_record(line));
    Ok(reencoded.then_some(lines))
}

/// Whether the text a gzip file inflates to ends in a line end.
fn ends_in_a_line_end(path: &str) -> bool {
    ByteSource::open_gzip(path)
        .and_then(ByteSource::read_to_end)
        .map_or(true, |bytes| matches!(bytes.last(), Some(b'\n' | b'\r')))
}

/// Ordinals of the records the reference carried with no row: the JSON output at
/// `json_path`, when the file exists, has an object for the record whose
/// `most_severe_consequence` is `?`, the object Ensembl VEP writes for a feature with no
/// overlap allele (every ALT a `*` or `<DEL:*>`, equal to the reference or a repeat of an
/// earlier one; VariationFeatureOverlap.pm:459-491) and for a record a row filter
/// (`--coding_only`, `--no_intergenic`) left no consequence of (JSON.pm:387-388). A record
/// is the object whose `input` opens with its first five columns, so a twin input carrying
/// sample columns names it too. A file that is not gzip, or holds a line that is not a JSON
/// object, is refused with its path; the reason a line is not JSON is the text of the
/// `json.JSONDecodeError` for the line with its line end.
pub fn carried_records(records: &[RunRecord], json_path: &str) -> Result<BTreeSet<usize>, Failure> {
    let mut carried = BTreeSet::new();
    if !is_file(json_path) {
        return Ok(carried);
    }
    let mut by_head: HashMap<&str, Vec<usize>> = HashMap::new();
    for (ordinal, record) in records.iter().enumerate() {
        by_head
            .entry(first_five_columns(&record.vcf))
            .or_default()
            .push(ordinal);
    }
    let not_json = |reason: &str| format!("{json_path} is not a readable JSON output ({reason})");
    let unreadable = |error: ReadError| error.failure(not_json);
    let mut lines = TextLines::open_gzip(json_path).map_err(unreadable)?;
    let mut more = lines.advance().map_err(unreadable)?;
    let mut number = 0usize;
    while more {
        number += 1;
        let line = lines.line().to_string();
        // The line after this one is read first: a line end followed this line unless it
        // is the file's last, and a file that fails further on fails after this line.
        let next = lines.advance();
        if !py_strip(&line).is_empty() {
            let object = match loads(&line) {
                Ok(Value::Object(object)) => object,
                Ok(_) => {
                    return Err(failure(not_json(&format!(
                        "line {number} is not a JSON object"
                    ))))
                }
                Err(JsonError::Decode { .. }) => {
                    let ended = !matches!(next, Ok(false)) || ends_in_a_line_end(json_path);
                    let text = if ended { format!("{line}\n") } else { line };
                    let reason = match loads(&text) {
                        Err(error) => error.to_string(),
                        Ok(_) => String::new(),
                    };
                    return Err(failure(not_json(&reason)));
                }
                Err(JsonError::Unsupported(detail)) => {
                    return Err(Failure::own(format!(
                        "{json_path} line {number} holds {detail}"
                    )))
                }
            };
            if object
                .get("most_severe_consequence")
                .and_then(Value::as_str)
                == Some("?")
            {
                let input = object.get("input").map(py_str).unwrap_or_default();
                if let Some(ordinals) = by_head.get(first_five_columns(&input)) {
                    carried.extend(ordinals.iter().copied());
                }
            }
        }
        more = next.map_err(unreadable)?;
    }
    Ok(carried)
}

/// Ordinals of the records no Ensembl VEP row names (Ensembl VEP dropped them before
/// annotation), and the records that contradict their `reference_warning`: a zero-row
/// record without one (its name unresolved, or the reason unrecorded) and a warned record
/// with rows. A record that carries no row by its nature is kept, not skipped, when it has
/// no warning: a non-variant record (ALT `.`, written under `--allow_non_variant` as its VCF
/// line and JSON object and no default row, Parser/VCF.pm:259-267) and a record of
/// `carried`, one the reference's JSON names as a feature with no allele to annotate
/// (`carried_records`), one a row filter or the frequency filter of the run left without a
/// row on either side, or one the frequency filter dropped and vep-rs kept for a zero
/// frequency. A problem names its record by its ordinal and the first five columns of its
/// VCF line, or, for a record a region line names (`region_lines`), by that line's position
/// (`region_position`), which holds no allele.
pub fn reference_skipped_records(
    records: &[RunRecord],
    perl_names: &TupleSets,
    by_name: &RunNames,
    carried: &BTreeSet<usize>,
) -> (BTreeSet<usize>, Vec<String>) {
    let with_rows: HashSet<usize> = perl_names
        .values()
        .flatten()
        .filter_map(|name| by_name.get(name))
        .flatten()
        .copied()
        .collect();
    let mut skipped: BTreeSet<usize> = (0..records.len())
        .filter(|ordinal| !with_rows.contains(ordinal))
        .collect();
    let mut problems = Vec::new();
    for (ordinal, record) in records.iter().enumerate() {
        let columns: Vec<&str> = record.vcf.split('\t').collect();
        let region = by_name.region_line(ordinal);
        let label = match region.and_then(region_position) {
            Some(position) => position.to_string(),
            None => columns[..columns.len().min(5)].join(" "),
        };
        let non_variant = region.is_none() && columns.len() > 4 && columns[4] == ".";
        match &record.reference_warning {
            None if skipped.contains(&ordinal) => {
                if non_variant || carried.contains(&ordinal) {
                    skipped.remove(&ordinal);
                } else {
                    problems.push(format!(
                        "record {ordinal} ({label}) has no VEP row and no reference_warning"
                    ));
                }
            }
            Some(warning) if !skipped.contains(&ordinal) => problems.push(format!(
                "record {ordinal} ({label}) has VEP rows and reference_warning {}",
                py_repr_str(warning)
            )),
            _ => {}
        }
    }
    (skipped, problems)
}

#[cfg(test)]
#[path = "records_tests.rs"]
mod tests;
