// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The SpliceAI field class of a row whose score record holds several entries: the
//! questions `classify` puts to the format comparator, whose rule the class is, the answers
//! it holds, and the entries that state the row's value under the names of the VCF and JSON
//! formats.

use super::comparator::{option, Comparator, Ended, Outcome};
use super::fields::{FieldKey, FieldModels, FieldQuestion, FieldRow, FieldRun, Verdict};
use super::protein::whole;
use crate::error::Failure;
use crate::json::{dumps, loads, Object, Style, Value};
use crate::pytext::{path_join, path_str};
use crate::rows::ByteSource;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::sync::LazyLock;

/// The class of a SpliceAI field on a row whose score record holds several comma-joined
/// entries: Ensembl VEP files the record under its first entry's symbol with every entry as
/// one value (SpliceAI.pm:326-332, :305), where vep-rs writes a row every entry of the
/// record that names the row's gene. The format comparator decides it from the score files `classify
/// --spliceai-table` names (`ask`).
pub const SPLICEAI_MULTI_ENTRY_CELL: &str = "spliceai_multi_entry_cell";

/// The fields of a row the class reads, which ride along with the compared fields: the
/// columns a question carries of each of its two rows, in this order.
pub const SPLICEAI_ROW_FIELDS: [&str; 12] = [
    "SpliceAI_pred",
    "SpliceAI_pred_SYMBOL",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL",
    "SpliceAI_pred_DP_DG",
    "SpliceAI_pred_DP_DL",
    "SpliceAI_cutoff",
    "SYMBOL",
];

/// The nine fields SpliceAI.pm:276-287 cut from a value for the VCF format, whatever the
/// run's options, and the key of the object the JSON format writes them under.
pub const SPLICEAI_MEMBER_FIELDS: [&str; 9] = [
    "SpliceAI_pred_SYMBOL",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_DP_AG",
    "SpliceAI_pred_DP_AL",
    "SpliceAI_pred_DP_DG",
    "SpliceAI_pred_DP_DL",
];
pub const SPLICEAI_JSON_FIELD: &str = "spliceai";

/// A member of a SpliceAI value that the JSON format writes as a number.
pub const SPLICEAI_NUMBER: &str = r"-?\d+(\.\d+)?";

static NUMBER: LazyLock<Regex> = LazyLock::new(|| whole(SPLICEAI_NUMBER));

/// `json.dumps(value, sort_keys=True, separators=(",", ":"))` of an object whose members are
/// set in the order of their names.
const COMPACT: Style = Style {
    indent: None,
    item_separator: ",",
    key_separator: ":",
    ensure_ascii: true,
};

/// The cell of a question for a field its row lacks.
const LACKING: &str = "-";

/// The fields the format comparator accepts for the class, by the key of their row.
pub type Accepted = HashMap<FieldKey, HashSet<String>>;

/// What the class reads of a run: whether the run names score files, and the fields the
/// format comparator accepts under them. The default is a run without score files.
#[derive(Default)]
pub struct Models {
    pub score_files: bool,
    pub accepted: Accepted,
}

impl Models {
    /// The models of a run before the format comparator is asked.
    pub fn of(run: &FieldRun<'_>) -> Models {
        Models {
            score_files: !run.spliceai_tables.is_empty(),
            accepted: Accepted::new(),
        }
    }

    /// The fields of a row the class reads beside the compared ones: none without the
    /// score files.
    pub fn row_fields(&self) -> &'static [&'static str] {
        if self.score_files {
            &SPLICEAI_ROW_FIELDS
        } else {
            &[]
        }
    }
}

/// One question of the class as a line of the command's questions: the key's Location,
/// Allele and Feature, the field, Ensembl VEP's value and vep-rs's, then the fields of
/// `SPLICEAI_ROW_FIELDS` as Ensembl VEP's row writes them and as vep-rs's row does, `-` for
/// a field a row lacks. The key's consequence set is no column: the class does not read it.
pub fn question_line(question: &FieldQuestion<'_>) -> String {
    let key = question.key;
    let mut cells = vec![
        key.location.as_str(),
        key.allele.as_str(),
        key.feature.as_str(),
        question.name,
        question.vep,
        question.vep_rs,
    ];
    for row in [question.perl_row, question.rust_row] {
        cells.extend(SPLICEAI_ROW_FIELDS.iter().map(|name| {
            row.and_then(|row| row.get(*name))
                .map_or(LACKING, String::as_str)
        }));
    }
    cells.join("\t") + "\n"
}

/// A refusal of the score files under the option's name.
fn refusal(detail: impl AsRef<str>) -> Failure {
    Failure::own(format!("--spliceai-table: {}", detail.as_ref()))
}

/// Why an index the command could not read is none, as this tool reads the file: the reason
/// it is no readable gzip file, else that its bytes are no tabix index.
fn index_reason(index: &str) -> String {
    match ByteSource::open_gzip(index).and_then(ByteSource::read_to_end) {
        Err(error) => error.reason().to_string(),
        Ok(_) => "not a tabix index".to_string(),
    }
}

/// The failure of a run on which the command left no answers, `tables` the score files as
/// the command was given them. A refusal of the command (exit code 2 and its own error
/// line) is this run's refusal with the same detail, two details worded by this tool: an
/// index that cannot be read, with the reason `index_reason` gives, and a block of a score
/// file that does not inflate. Any other ending names the build and how it ended.
pub fn failure(comparator: &Comparator, tables: &[String], ended: &Ended) -> Failure {
    let Some(detail) = ended.detail().filter(|_| ended.code == Some(2)) else {
        return comparator.ended(ended);
    };
    for table in tables {
        let index = format!("{table}.tbi");
        let unread = format!("--spliceai-table: cannot read the tabix index {index}: ");
        if detail.starts_with(&unread) {
            let reason = index_reason(&index);
            return refusal(format!("{index} is not a readable tabix index ({reason})"));
        }
        if detail.starts_with(&format!("--spliceai-table: cannot read {table}: Error -")) {
            return refusal(format!("{table} holds a block that does not inflate"));
        }
    }
    Failure::own(detail)
}

/// The fields of `questions` the class accepts: the answers of one run of the format
/// comparator's `field-class` command over every question, under the run's score files and
/// the corpus cache where it holds a `transcripts` directory. A run without a question is
/// asked as well, since the command reads the index and the header of every score file
/// whatever it is asked. A run of the command that leaves no answers is this run's
/// `failure`, and so are answers that are not one `0` or `1` a question.
pub fn ask(
    comparator: &Comparator,
    run: &FieldRun<'_>,
    questions: &[FieldQuestion<'_>],
) -> Result<Accepted, Failure> {
    let tables: Vec<String> = run
        .spliceai_tables
        .iter()
        .map(|table| path_str(table))
        .collect();
    let mut options: Vec<OsString> = tables
        .iter()
        .map(|table| option("--spliceai-table", table))
        .collect();
    let transcripts = path_join(run.cache_dir, "transcripts");
    if std::fs::metadata(transcripts).is_ok_and(|held| held.is_dir()) {
        options.push(option("--json-cache", run.cache_dir));
    }
    let lines: String = questions.iter().map(question_line).collect();
    let answered = match comparator.field_class(SPLICEAI_MULTI_ENTRY_CELL, &options, &lines)? {
        Outcome::Answered(answered) => answered,
        Outcome::Ended(ended) => return Err(failure(comparator, &tables, &ended)),
    };
    let accepts = comparator.accepted(SPLICEAI_MULTI_ENTRY_CELL, questions.len(), &answered)?;
    let mut accepted = Accepted::new();
    for (question, accepts) in questions.iter().zip(accepts) {
        if accepts {
            let fields = accepted.entry(question.key.clone()).or_default();
            fields.insert(question.name.to_string());
        }
    }
    Ok(accepted)
}

/// The class `spliceai_multi_entry_cell` of the dispatch: the field is one the format
/// comparator accepts on the question's key (`Models::accepted`). Nothing in a run without
/// score files, whose models hold no answer.
pub fn spliceai_multi_entry_cell(question: &FieldQuestion<'_>, models: &FieldModels) -> Verdict {
    let accepted = models.spliceai.accepted.get(question.key);
    if accepted.is_some_and(|fields| fields.contains(question.name)) {
        return Verdict::Class(SPLICEAI_MULTI_ENTRY_CELL, None);
    }
    Verdict::Pass
}

/// The nine members of a default-format row's SpliceAI value: the first nine `|`-joined
/// parts of its `SpliceAI_pred`, else its nine `SpliceAI_pred_` fields; `None` for a row
/// that writes neither whole.
pub fn spliceai_members(row: &FieldRow) -> Option<Vec<String>> {
    let field = |name: &str| row.get(name).map_or("-", String::as_str);
    let members: Vec<String> = match field("SpliceAI_pred") {
        "-" => SPLICEAI_MEMBER_FIELDS
            .iter()
            .map(|name| field(name).to_string())
            .collect(),
        value => value.split('|').take(9).map(str::to_string).collect(),
    };
    (members.len() == 9 && !members.iter().any(|member| member == "-")).then_some(members)
}

/// The nine members of vep-rs's default-format row: where its `SpliceAI_pred` holds several
/// comma-joined entries, the first behind its allele and each later one with its allele
/// part, each member is its values over the entries joined by commas, as vep-rs writes the
/// nine fields of a row whose gene several entries of one score record name; any other row
/// is read as `spliceai_members` reads it. `None` for a row that writes no whole value or
/// an entry of fewer than nine parts.
pub fn spliceai_entry_members(row: &FieldRow) -> Option<Vec<String>> {
    let value = row.get("SpliceAI_pred").map_or("-", String::as_str);
    if value == "-" || !value.contains(',') {
        return spliceai_members(row);
    }
    let entries: Vec<Vec<&str>> = value
        .split(',')
        .enumerate()
        .map(|(place, entry)| entry.split('|').skip(usize::from(place > 0)).collect())
        .collect();
    if entries.iter().any(|entry| entry.len() < 9) {
        return None;
    }
    let members = (0..9).map(|at| {
        let values: Vec<&str> = entries.iter().map(|entry| entry[at]).collect();
        values.join(",")
    });
    Some(members.collect())
}

/// The question of the class about a row whose compared fields agree: the first member
/// field whose two values differ, Ensembl VEP's member and vep-rs's.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MemberQuestion {
    pub name: &'static str,
    pub vep: String,
    pub vep_rs: String,
}

/// The member question of a row on which both engines write one `SpliceAI_pred` text: where
/// the nine members Ensembl VEP cuts from it (`spliceai_members`) are not the nine vep-rs
/// writes for it (`spliceai_entry_members`), the VCF and JSON formats differ on a row the
/// default format writes alike, as on a row whose gene alone the several entries of its
/// score record name. `None` for any other pair of rows.
pub fn member_question(perl_row: &FieldRow, rust_row: &FieldRow) -> Option<MemberQuestion> {
    let value = perl_row
        .get("SpliceAI_pred")
        .filter(|value| *value != "-")?;
    if rust_row.get("SpliceAI_pred") != Some(value) {
        return None;
    }
    let (theirs, ours) = (
        spliceai_members(perl_row)?,
        spliceai_entry_members(rust_row)?,
    );
    let at = (0..9).find(|at| theirs[*at] != ours[*at])?;
    Some(MemberQuestion {
        name: SPLICEAI_MEMBER_FIELDS[at],
        vep: theirs[at].clone(),
        vep_rs: ours[at].clone(),
    })
}

/// A member's text as the number both engines write for it: an integer where the value is
/// one, so a score written `0.00` is `0` and `0.10` is `0.1`.
fn number(text: &str) -> Value {
    let value: f64 = text.parse().unwrap_or(f64::NAN);
    if value == 0.0 {
        return Value::from(0i64);
    }
    if value.is_finite() && value.fract() == 0.0 {
        if let Ok(integer) = loads(&format!("{value:.0}")) {
            return integer;
        }
    }
    Value::Float(value)
}

/// A row's JSON `spliceai` object as the golden harness prints an object: no spaces, its
/// members in the byte order of their names, the symbol a text, any other member a number
/// where its text is one and a text otherwise, and the cutoff flag among them where the row
/// carries one; `-` for a row without the value. Both engines write a number by its value.
pub fn spliceai_json_object(members: Option<&[String]>, flag: &str) -> String {
    let Some(members) = members else {
        return "-".to_string();
    };
    let mut written: Vec<(&str, Value)> = SPLICEAI_MEMBER_FIELDS
        .iter()
        .zip(members)
        .map(|(field, member)| {
            let name = field.strip_prefix("SpliceAI_pred_").unwrap_or(field);
            let value = if name != "SYMBOL" && NUMBER.is_match(member) {
                number(member)
            } else {
                Value::from(member.as_str())
            };
            (name, value)
        })
        .collect();
    if flag != "-" {
        written.push(("SpliceAI_cutoff", Value::from(flag)));
    }
    written.sort_by_key(|(name, _)| *name);
    let object: Object = written.into_iter().collect();
    dumps(&Value::Object(object), &COMPACT)
}

/// The entries the VCF and JSON formats need for a row `classify` writes under
/// `spliceai_multi_entry_cell` (`written`, one of its entries), each value derived from the
/// two rows' default-format values: one a member field of `SPLICEAI_MEMBER_FIELDS` whose
/// two members differ, left out where `compared` names the field and the run's default
/// format writes it itself and where `written` is that field's own entry, and one for the
/// JSON object (`spliceai_json_object`). Ensembl VEP's members are the nine it cuts from
/// its value and vep-rs's the nine it writes (`spliceai_entry_members`). The
/// golden harness reads a documented field by the name its format gives it.
pub fn spliceai_other_format_entries(
    written: &Object,
    perl_row: &FieldRow,
    rust_row: &FieldRow,
    compared: &[String],
) -> Vec<Object> {
    let (theirs, ours) = (spliceai_members(perl_row), spliceai_entry_members(rust_row));
    let own_field = written.get("field").and_then(Value::as_str);
    let entry = |field: &str, vep: &str, vep_rs: &str| {
        let mut entry = written.clone();
        entry.insert("field", field);
        entry.insert("vep", vep);
        entry.insert("vep_rs", vep_rs);
        entry
    };
    let member = |members: &Option<Vec<String>>, at: usize| {
        members
            .as_ref()
            .map_or_else(|| "-".to_string(), |members| members[at].clone())
    };
    let mut entries = Vec::new();
    for (at, field) in SPLICEAI_MEMBER_FIELDS.iter().enumerate() {
        let (vep, vep_rs) = (member(&theirs, at), member(&ours, at));
        let written_already =
            compared.iter().any(|name| name == field) || own_field == Some(*field);
        if vep != vep_rs && !written_already {
            entries.push(entry(field, &vep, &vep_rs));
        }
    }
    fn flag(row: &FieldRow) -> &str {
        row.get("SpliceAI_cutoff").map_or("-", String::as_str)
    }
    entries.push(entry(
        SPLICEAI_JSON_FIELD,
        &spliceai_json_object(theirs.as_deref(), flag(perl_row)),
        &spliceai_json_object(ours.as_deref(), flag(rust_row)),
    ));
    entries
}

#[cfg(test)]
#[path = "spliceai_tests.rs"]
mod tests;
