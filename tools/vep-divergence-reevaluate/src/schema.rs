// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The schema of a reproduction record: the keys each table may hold, the type of each
//! value, what each compare kind needs, the form of the input class's basis and the shape
//! of the `implementation` table. A violation is one line; `reproductions/README.md`
//! beside the records states the schema.

use crate::page::{PagePatterns, INPUT_CLASSES};
use crate::pyfmt::{py_repr_str, py_splitlines, py_strip};
use crate::pyre::Regex;
use crate::record::{Reproduction, RUN_TABLES};
use crate::value::{as_list, Value, EMPTY_DICT};

pub const VERDICTS: [&str; 4] = ["still-present", "closed", "changed", "error"];
/// The output formats and the file name each is written under.
pub const FORMATS: [(&str, &str); 4] = [
    ("default", "default.txt"),
    ("tab", "tab.txt"),
    ("vcf", "vcf.vcf"),
    ("json", "json.jsonl"),
];
/// The input formats and the extension an inline input of each is written with.
pub const INPUT_EXTENSIONS: [(&str, &str); 7] = [
    ("vcf", "vcf"),
    ("ensembl", "txt"),
    ("hgvs", "txt"),
    ("spdi", "txt"),
    ("region", "txt"),
    ("id", "txt"),
    ("guess", "txt"),
];
pub const COMPARE_KINDS: [&str; 12] = [
    "rows", "field", "count", "header", "lines", "warnings", "stderr", "stdout", "exit", "stats",
    "json", "vcf",
];
const TABLE_WHERE_KEYS: [&str; 8] = [
    "record",
    "location",
    "allele",
    "gene",
    "feature",
    "feature_type",
    "consequence",
    "pattern",
];
const RECORD_WHERE_KEYS: [&str; 6] = ["record", "id", "chrom", "pos", "location", "pattern"];
const JSON_WHERE_KEYS: [&str; 7] = [
    "record", "id", "chrom", "pos", "location", "pattern", "input",
];
const VCF_WHERE_KEYS: [&str; 8] = [
    "record", "id", "chrom", "pos", "location", "pattern", "feature", "allele",
];
/// The programs a run names and Ensembl VEP's name for each.
pub const TOOLS: [(&str, &str); 3] = [
    ("vep", "vep"),
    ("filter", "filter_vep"),
    ("recode", "variant_recoder"),
];

/// The keys of each table; the empty name is the record's top level.
const KEYS: [(&str, &[&str]); 7] = [
    (
        "",
        &[
            "entry",
            "page",
            "verified_against",
            "image",
            "input_class",
            "input_class_basis",
            "notes",
            "prepare",
            "input",
            "cache",
            "run",
            "compare",
            "expected",
            "also",
            "implementation",
        ],
    ),
    (
        "input",
        &["text", "path", "format", "header", "samples", "name"],
    ),
    (
        "cache",
        &[
            "species",
            "assembly",
            "source",
            "version",
            "json_cache",
            "corpus",
            "fasta",
            "gxf",
            "gxf_kind",
            "perl_cache",
            "chromosomes",
            "writable",
        ],
    ),
    (
        "run",
        &[
            "flags",
            "ensembl_vep_flags",
            "vep_rs_flags",
            "format",
            "format_flag",
            "stats",
            "buffer_size",
            "tool",
            "args",
            "stdin",
            "perl_hash_seed",
            "perl_mode",
            "assembly_flag",
            "tools",
        ],
    ),
    (
        "compare",
        &[
            "kind",
            "where",
            "field",
            "columns",
            "with_key",
            "pattern",
            "path",
            "section",
            "info",
            "line",
            "mask_paths",
            "mask",
            "sorted",
            "strip",
        ],
    ),
    (
        "expected",
        &["ensembl_vep", "vep_rs", "ensembl_vep_exit", "vep_rs_exit"],
    ),
    ("implementation", &["code", "tests", "rules", "none"]),
];

/// The keys whose value is a string.
const STRING_KEYS: [(&str, &str); 26] = [
    ("", "notes"),
    ("input", "text"),
    ("input", "path"),
    ("input", "format"),
    ("input", "name"),
    ("cache", "species"),
    ("cache", "assembly"),
    ("cache", "source"),
    ("cache", "version"),
    ("cache", "json_cache"),
    ("cache", "corpus"),
    ("cache", "fasta"),
    ("cache", "gxf"),
    ("cache", "gxf_kind"),
    ("cache", "perl_cache"),
    ("run", "format"),
    ("run", "tool"),
    ("run", "perl_hash_seed"),
    ("run", "perl_mode"),
    ("compare", "kind"),
    ("compare", "field"),
    ("compare", "pattern"),
    ("compare", "path"),
    ("compare", "section"),
    ("compare", "info"),
    ("implementation", "none"),
];
/// The keys whose value is a list of strings.
const LIST_KEYS: [(&str, &str); 14] = [
    ("", "prepare"),
    ("input", "header"),
    ("input", "samples"),
    ("cache", "chromosomes"),
    ("run", "flags"),
    ("run", "ensembl_vep_flags"),
    ("run", "vep_rs_flags"),
    ("run", "args"),
    ("compare", "columns"),
    ("compare", "mask"),
    ("run", "tools"),
    ("implementation", "code"),
    ("implementation", "tests"),
    ("implementation", "rules"),
];

/// The keys of the `implementation` table that list items: the code that holds the
/// recorded behaviour, the tests that pin vep-rs's row and the comparison rules that
/// recognise the difference.
pub const IMPLEMENTATION_LISTS: [&str; 3] = ["code", "tests", "rules"];

/// The fixed words of an `unobserved` entry's basis, around the format or option contract
/// it names; a basis may go on after the closing words.
const UNOBSERVED_OPENING: &str = "well-formed under ";
const UNOBSERVED_CLOSING: &str = "; not seen in the measured datasets";
/// The basis of an `observed` entry of tier 3: the real data its page names, up to the
/// first colon a space follows, then the figure the page measures there.
pub const SIGHTING: &str = r"^seen in (?P<where>[^\s:](?:(?!: ).)*)(?<=\S): (?P<what>\S.*)$";
pub const SIGHTING_FORM: &str =
    "seen in <the real data its page names>: <the figure its page measures there>";

/// The path and the item of an `implementation` row, `<path>::<item>`: a path relative to
/// the repository that climbs out of no directory, and an item's chain of names. `None`
/// for a row of another shape.
pub fn implementation_row(row: &str) -> Option<(&str, &str)> {
    let (path, item) = row.split_once("::")?;
    let relative = !path.is_empty()
        && !path.starts_with('/')
        && path
            .split('/')
            .all(|part| !part.is_empty() && part != "." && part != "..");
    let named = !item.is_empty()
        && item.split("::").all(|part| {
            !part.is_empty()
                && part
                    .chars()
                    .all(|ch| ch.is_ascii_alphanumeric() || ch == '_')
        });
    (relative && named).then_some((path, item))
}

/// What is wrong with the shape of an `implementation` table whose values have their
/// types: a row that is not `<path>::<item>`, a table with neither `code` nor `none` or
/// with both, a blank `none`, and a table without `tests`.
fn implementation_errors(table: &Value) -> Vec<String> {
    let mut errors = Vec::new();
    let rows = |key: &str| match table.get(key) {
        Some(Value::List(rows)) => rows.as_slice(),
        _ => &[],
    };
    for key in IMPLEMENTATION_LISTS {
        for row in rows(key) {
            let row = row.as_str().unwrap_or_default();
            if implementation_row(row).is_none() {
                errors.push(format!(
                    "implementation.{key} names {}, which is not <path>::<item>",
                    py_repr_str(row)
                ));
            }
        }
    }
    let none = table.get_str("none");
    if none.is_some_and(|text| py_strip(text).is_empty() || py_splitlines(text) != [text]) {
        errors.push("implementation.none is not one line of text".to_string());
    }
    match (rows("code").is_empty(), none.is_some()) {
        (true, false) => errors.push("implementation needs code or none".to_string()),
        (false, true) => errors.push(
            "implementation holds code and none: none is for an entry with no code of its own"
                .to_string(),
        ),
        _ => {}
    }
    if rows("tests").is_empty() {
        errors.push("implementation needs tests".to_string());
    }
    errors
}

fn names<const N: usize>(pairs: &[(&str, &str); N]) -> String {
    pairs
        .iter()
        .map(|(name, _)| *name)
        .collect::<Vec<_>>()
        .join(", ")
}

/// `repr()` of what `dict.get(key)` returns: `None` for an absent key.
fn repr_at(table: &Value, key: &str) -> String {
    table
        .get(key)
        .map_or_else(|| "None".to_string(), Value::py_repr)
}

/// What a schema check holds a record against, compiled once.
pub struct Schema {
    pub page: PagePatterns,
    pub sighting: Regex,
}

impl Schema {
    pub fn new() -> Self {
        Schema {
            page: PagePatterns::new(),
            sighting: Regex::new(SIGHTING).expect("the sighting pattern"),
        }
    }

    /// What is wrong with a record's `input_class` and `input_class_basis`: one of the
    /// three classes, and one line of basis in the form the class asks for.
    pub fn input_class_errors(&self, data: &Value) -> Vec<String> {
        let mut errors = Vec::new();
        for key in ["input_class", "input_class_basis"] {
            if !data
                .get_str(key)
                .is_some_and(|text| !py_strip(text).is_empty())
            {
                errors.push(format!("{key} is missing"));
            }
        }
        if !errors.is_empty() {
            return errors;
        }
        let input_class = data.get_str("input_class").expect("a string");
        let basis = data.get_str("input_class_basis").expect("a string");
        if !INPUT_CLASSES.contains(&input_class) {
            return vec![format!(
                "input_class {} is not one of {}",
                py_repr_str(input_class),
                INPUT_CLASSES.join(", ")
            )];
        }
        if py_splitlines(basis) != [basis] {
            return vec!["input_class_basis is not one line".to_string()];
        }
        let cell_tier = self
            .page
            .cell
            .match_start(basis)
            .and_then(|m| m.name("tier"));
        let is_sighting = self.sighting.match_start(basis).is_some();
        match input_class {
            "observed" => {
                let tier_1_or_2_cell = cell_tier.is_some_and(|tier| tier == "1" || tier == "2");
                if !tier_1_or_2_cell && !is_sighting {
                    errors.push(format!(
                        "input_class_basis of an observed entry is neither a tier 1 or tier 2 `Tier and count` cell nor `{SIGHTING_FORM}`"
                    ));
                }
            }
            "unobserved" => {
                let well_formed = basis.strip_prefix(UNOBSERVED_OPENING).is_some_and(|rest| {
                    rest.split_once(UNOBSERVED_CLOSING)
                        .is_some_and(|(contract, _)| !py_strip(contract).is_empty())
                });
                if !well_formed {
                    errors.push(format!(
                        "input_class_basis of an unobserved entry does not read `{UNOBSERVED_OPENING}<format or option contract>{UNOBSERVED_CLOSING}`"
                    ));
                }
            }
            _ => {
                if cell_tier.is_some() || is_sighting || basis.starts_with(UNOBSERVED_OPENING) {
                    errors.push(
                        "input_class_basis of a malformed entry reads as the basis of another class, not as the rule the input breaks"
                            .to_string(),
                    );
                }
            }
        }
        errors
    }

    /// The schema violations of one record, as one line each.
    pub fn validate(&self, repro: &Reproduction) -> Vec<String> {
        let d = &repro.data;
        let mut errors = Vec::new();
        let prefix = |table: &str| {
            if table.is_empty() {
                String::new()
            } else {
                format!("{table}.")
            }
        };
        let section = |table: &str| -> &Value {
            if table.is_empty() {
                d
            } else {
                d.get(table).unwrap_or(&EMPTY_DICT)
            }
        };
        let mut tables = true;
        for (table, allowed) in KEYS {
            let section = section(table);
            if !section.is_dict() {
                errors.push(format!("{table} is not a table"));
                tables = false;
                continue;
            }
            for (key, _) in section.items() {
                if !allowed.contains(&key.as_str()) {
                    errors.push(format!("unknown key {}{key}", prefix(table)));
                }
            }
        }
        if !tables {
            // A section that is not a table has no keys to read.
            return errors;
        }
        if d.get_str("entry") != Some(repro.name.as_str()) {
            errors.push(format!(
                "entry is {}, the file is {}.toml",
                repr_at(d, "entry"),
                repro.name
            ));
        }
        for key in ["page", "verified_against", "image"] {
            if !d.get_str(key).is_some_and(|text| !text.is_empty()) {
                errors.push(format!("{key} is missing"));
            }
        }
        if d.get_str("image")
            .is_some_and(|image| !image.contains("@sha256:"))
        {
            errors.push("image is not pinned by digest".to_string());
        }
        errors.extend(self.input_class_errors(d));
        let mut typed = true;
        for (table, key) in STRING_KEYS {
            if let Some(value) = section(table).get(key) {
                if value.as_str().is_none() {
                    errors.push(format!("{}{key} is not a string", prefix(table)));
                    typed = false;
                }
            }
        }
        for (table, key) in LIST_KEYS {
            if let Some(value) = section(table).get(key) {
                let strings = matches!(value, Value::List(items) if items.iter().all(|v| v.as_str().is_some()));
                if !strings {
                    errors.push(format!("{}{key} is not a list of strings", prefix(table)));
                    typed = false;
                }
            }
        }
        if !typed {
            // The checks below read the typed values.
            return errors;
        }
        let input = repro.input();
        if input.has("text") == input.has("path") {
            errors.push("input needs exactly one of text and path".to_string());
        }
        let format = input.get_str("format").unwrap_or("vcf");
        if !INPUT_EXTENSIONS.iter().any(|(name, _)| *name == format) {
            errors.push(format!(
                "input.format {} is not one of {}",
                repr_at(input, "format"),
                names(&INPUT_EXTENSIONS)
            ));
        }
        let cache = repro.cache();
        let run = repro.run();
        let tool = run.get_str("tool").unwrap_or("vep");
        if !TOOLS.iter().any(|(name, _)| *name == tool) {
            errors.push(format!(
                "run.tool {} is not one of {}",
                py_repr_str(tool),
                names(&TOOLS)
            ));
        }
        if tool == "vep" {
            if !cache.truthy_at("assembly") {
                errors.push("cache.assembly is missing".to_string());
            }
            if !(cache.truthy_at("json_cache") || cache.truthy_at("gxf")) {
                errors.push("cache needs json_cache or gxf".to_string());
            }
            let source = cache.get_str("source").unwrap_or("ensembl");
            if !["ensembl", "merged", "refseq"].contains(&source) {
                errors.push("cache.source is not ensembl, merged or refseq".to_string());
            }
        } else if !run.truthy_at("args") {
            errors.push("a tool run needs run.args".to_string());
        }
        match d.get("compare") {
            None => errors.push("compare is missing".to_string()),
            Some(compare) => self.validate_compare(compare, &mut errors),
        }
        let expected = d.get("expected");
        if let Some(expected) = expected {
            for key in ["ensembl_vep_exit", "vep_rs_exit"] {
                if expected.get(key).is_some_and(|value| !value.is_int()) {
                    errors.push(format!("expected.{key} is not an integer"));
                }
            }
        }
        match expected.map(|e| (e.get("ensembl_vep"), e.get("vep_rs"))) {
            Some((Some(ensembl_vep), Some(vep_rs))) => {
                if as_list(ensembl_vep) == as_list(vep_rs) {
                    errors.push(
                        "expected.ensembl_vep equals expected.vep_rs: the entry records no divergence"
                            .to_string(),
                    );
                }
            }
            _ => errors.push("expected needs ensembl_vep and vep_rs".to_string()),
        }
        let perl_mode = run.get_str("perl_mode").unwrap_or("offline");
        if perl_mode != "offline" && perl_mode != "database" {
            errors.push("run.perl_mode is not offline or database".to_string());
        }
        if let Some(format) = run.get_str("format") {
            if !FORMATS.iter().any(|(name, _)| *name == format) {
                errors.push(format!(
                    "run.format {} is not one of {}",
                    py_repr_str(format),
                    names(&FORMATS)
                ));
            }
        }
        if let Some(implementation) = d.get("implementation") {
            errors.extend(implementation_errors(implementation));
        }
        if errors.is_empty() {
            // A further run is the record with some tables replaced, so it is held to the
            // schema once the record's own tables pass it.
            errors.extend(self.also_errors(repro));
        }
        errors
    }

    /// What is wrong with the further runs of a record whose own tables pass the schema:
    /// an `also` that is no list of tables, a key a further run does not take, a run
    /// without values of its own, and each violation of the record as that run reads it,
    /// named by the run's place in the list, counted from 1.
    fn also_errors(&self, repro: &Reproduction) -> Vec<String> {
        let tables = match repro.data.get("also") {
            None => return Vec::new(),
            Some(Value::List(tables)) if tables.iter().all(Value::is_dict) => tables,
            Some(_) => return vec!["also is not a list of tables".to_string()],
        };
        let mut errors = Vec::new();
        for (index, table) in tables.iter().enumerate() {
            let place = index + 1;
            for (key, _) in table.items() {
                if key != "notes" && !RUN_TABLES.contains(&key.as_str()) {
                    errors.push(format!("unknown key also[{place}].{key}"));
                }
            }
            if table.get("notes").is_some_and(|v| v.as_str().is_none()) {
                errors.push(format!("also[{place}].notes is not a string"));
            }
            if !table.has("expected") {
                errors.push(format!("also[{place}] needs expected"));
            }
        }
        if !errors.is_empty() {
            return errors;
        }
        for (place, run) in repro.runs().iter().enumerate().skip(1) {
            for problem in self.validate(run) {
                errors.push(format!("also[{place}]: {problem}"));
            }
        }
        errors
    }

    fn validate_compare(&self, compare: &Value, errors: &mut Vec<String>) {
        let kind = compare.get_str("kind");
        // `None` where the record names no kind, as the messages below print it.
        let kind_text = kind.unwrap_or("None");
        if !kind.is_some_and(|kind| COMPARE_KINDS.contains(&kind)) {
            errors.push(format!(
                "compare.kind {} is not one of {}",
                repr_at(compare, "kind"),
                COMPARE_KINDS.join(", ")
            ));
        }
        let whole_line = compare.truthy_at("line");
        let needs_field = kind == Some("field") || (kind == Some("vcf") && !whole_line);
        if needs_field && !compare.truthy_at("field") {
            errors.push(format!("compare.kind {kind_text} needs compare.field"));
        }
        let needs_pattern = matches!(
            kind,
            Some("header" | "lines" | "warnings" | "stderr" | "stdout")
        );
        if needs_pattern && !compare.truthy_at("pattern") {
            errors.push(format!("compare.kind {kind_text} needs compare.pattern"));
        }
        if kind == Some("json") && !compare.truthy_at("path") {
            errors.push("compare.kind json needs compare.path".to_string());
        }
        if kind == Some("stats") && !compare.truthy_at("section") {
            errors.push("compare.kind stats needs compare.section".to_string());
        }
        let mut where_ = compare.section("where");
        if !where_.is_dict() {
            errors.push("compare.where is not a table".to_string());
            where_ = &EMPTY_DICT;
        }
        let allowed: &[&str] = match kind {
            Some("rows" | "field" | "count") => &TABLE_WHERE_KEYS,
            Some("json") => &JSON_WHERE_KEYS,
            // A whole data line is matched at the record level alone.
            Some("vcf") if whole_line => &RECORD_WHERE_KEYS,
            Some("vcf") => &VCF_WHERE_KEYS,
            _ => &[],
        };
        let unknown: Vec<&str> = where_
            .items()
            .iter()
            .map(|(key, _)| key.as_str())
            .filter(|key| !allowed.contains(key))
            .collect();
        if !unknown.is_empty() {
            errors.push(format!(
                "compare.where has keys {} that compare.kind {kind_text} does not read",
                unknown.join(", ")
            ));
        }
        let mut patterns: Vec<(&str, &Value)> = Vec::new();
        if let Some(pattern) = compare.get("pattern") {
            patterns.push(("compare.pattern", pattern));
        }
        if let Some(pattern) = where_.get("pattern") {
            patterns.push(("compare.where.pattern", pattern));
        }
        if let Some(Value::List(mask)) = compare.get("mask") {
            patterns.extend(mask.iter().map(|pattern| ("compare.mask", pattern)));
        }
        for (label, pattern) in patterns {
            let refusal = match pattern {
                Value::Str(pattern) => Regex::new(pattern).err().map(|error| error.to_string()),
                Value::List(_) | Value::Dict(_) => Some(format!(
                    "cannot use 'tuple' as a dict key (unhashable type: '{}')",
                    pattern.py_type()
                )),
                _ => Some("first argument must be string or compiled pattern".to_string()),
            };
            if let Some(refusal) = refusal {
                errors.push(format!("{label} is not a regular expression: {refusal}"));
            }
        }
    }
}

impl Default for Schema {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
#[path = "schema_tests.rs"]
mod tests;
