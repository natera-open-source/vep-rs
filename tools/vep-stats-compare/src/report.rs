// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The sections of a statistics report and the comparison of two reports.
//!
//! The grammar read is `Stats.pm:806-831` at release/116.3: a `[Section]` line opens each
//! section, a row is `key<TAB>value`, sections are separated by one blank line, and the HTML
//! tags of the three tabular sections are stripped before they are written. The three
//! tabular sections come first and in this order: `[VEP run statistics]`
//! (`Stats.pm:624-681`), `[Data version]` (`:695-700`; the one line
//! `No data version information available.` when the cache names no data versions) and
//! `[General statistics]` (`:716-740`, a row whose value is undefined is written as its key
//! alone). Every later section is a chart (`:505-609`): `key<TAB>count` rows in the chart's
//! sort order.
//!
//! The rules:
//!
//! * A tabular section compares row by row: the same keys in the same order, every value
//!   verbatim, except the masked rows, whose values no second run reproduces and which
//!   compare by key alone ([`MASKED_ROWS`] and every `--mask-row` label); in
//!   `Annotation sources` every path (a token carrying `/`) is replaced by `[PATH]` before
//!   the values compare, so `Cache: /a/cache` and `Cache: [PATH]/116_GRCh38` agree.
//! * A chart section compares as a set of `(key, count)` pairs, since `sort_keys`
//!   (`Stats.pm:760-790`) breaks ties in hash order: two keys with one count, or two
//!   non-numeric chromosome names sharing a first character, come out in an order that
//!   differs from run to run. Under a tolerance of at least 1 the two counts of a key may
//!   differ by up to the tolerance, a key absent on one side counting as 0 there; at the
//!   tolerance 0 the two sets must be equal.
//! * The two files must carry the same sections in the same order: a chart Ensembl VEP
//!   writes only under a flag (`[SIFT summary]`), only with a coding row
//!   (`[Position in protein]`) or only with a FASTA (the per-chromosome distributions) is
//!   expected of the candidate under the same run.
//! * One row is held to a documented difference: `Novel / existing variants` of
//!   `[General statistics]`, where the expected report's command line switches the
//!   frequency check on and its row counts every processed variant as existing. The
//!   candidate's row then agrees when it is the row of some count of existing variants
//!   among the expected report's `Variants processed`, and the report names the row under
//!   `documented` with its class, [`EXISTING_VARIANTS_CLASS`].
//! * One kind of key is held to a second documented difference: a key of
//!   `[Consequences (all)]` that the expected report alone holds and that names no
//!   consequence type, where the expected report's command line names a `--custom` source
//!   whose key is `Consequence`. Ensembl VEP counts the names of that source's records
//!   among the consequence types, so the key is set aside and the report names it under
//!   `documented` with its class, [`CONSEQUENCE_KEY_CLASS`]. A key that names a consequence
//!   type, and a key the candidate holds, compare as any other.

use crate::json::Json;
use crate::pyfmt::{is_py_space, py_repr_str, py_strip, surrogateescape, PyInt};
use crate::terms::is_consequence_term;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};

/// The label both `mode` and `format` of the report carry, and the first word on stdout.
pub const MODE: &str = "stats";

pub const RUN_STATISTICS: &str = "VEP run statistics";
pub const DATA_VERSION: &str = "Data version";
pub const GENERAL_STATISTICS: &str = "General statistics";
const TABULAR_SECTIONS: [&str; 3] = [RUN_STATISTICS, DATA_VERSION, GENERAL_STATISTICS];

/// The rows of `[VEP run statistics]` whose values differ between any two runs
/// (`Stats.pm:645-663`): the command line names the engine's own options and paths, the
/// times and the run time are the run's own, the input and output rows are paths.
pub const MASKED_ROWS: [&str; 7] = [
    "Command line options",
    "Start time",
    "End time",
    "Run time",
    "Input file",
    "Input data",
    "Output file",
];

/// `Stats.pm:666-677`: `Cache: <dir>`, `<db> on <host>` and `Custom: <file> (<type>)` joined
/// by `; `; the directory and the files are paths.
const PATH_ROW: &str = "Annotation sources";
const MASKED: &str = "<masked>";
const PATH: &str = "[PATH]";

/// The documented difference of the `Novel / existing variants` row: under a frequency
/// check Ensembl VEP counts a variant as existing where it holds a list of known variants
/// (`Stats.pm:244`), and the check leaves every variant such a list, an empty one too.
pub const EXISTING_VARIANTS_CLASS: &str = "stats_existing_variants_check_frequency";
const OPTIONS_ROW: &str = "Command line options";
const PROCESSED_ROW: &str = "Variants processed";
const EXISTING_ROW: &str = "Novel / existing variants";

/// `--filter_common` sets `check_frequency` (`Config.pm:441-446`).
const FREQUENCY_FLAGS: [&str; 2] = ["--check_frequency", "--filter_common"];

/// The documented difference of the `Consequences (all)` chart: Ensembl VEP counts every
/// item of a row's `Consequence` list (`Stats.pm:282`), and under a `--custom` source whose
/// key is `Consequence` the list holds the names of the source's records on the row's
/// allele after the row's terms (`OutputFactory.pm:1301`, called at `:1426-1435`).
pub const CONSEQUENCE_KEY_CLASS: &str = "stats_consequences_custom_key_named_consequence";
const CONSEQUENCES_ALL: &str = "Consequences (all)";
const CUSTOM_FLAG: &str = "--custom";
const CONSEQUENCE_KEY: &str = "Consequence";

/// A double holds every integer up to here, so a percentage of two counts within it is
/// the one Perl computes.
const MAX_COUNT: u64 = 1 << 53;

/// One `[title]` block: its rows as `(key, value)` pairs in file order, a row without a tab
/// carrying no value.
#[derive(Debug, PartialEq, Eq)]
pub struct Section {
    pub title: String,
    pub rows: Vec<(String, Option<String>)>,
}

impl Section {
    /// Whether the section is a chart: any section but the three tabular ones.
    pub fn is_chart(&self) -> bool {
        !TABULAR_SECTIONS.contains(&self.title.as_str())
    }
}

/// Builds the sections of a statistics report from its lines in file order. A non-blank
/// line before the first section header is refused: the file is not a statistics report.
#[derive(Default)]
pub struct Parser {
    sections: Vec<Section>,
    refused: Option<String>,
    lines: u64,
}

impl Parser {
    /// Takes the next line of the report, without its line ending.
    pub fn line(&mut self, line: &str) {
        self.lines += 1;
        if self.refused.is_some() || py_strip(line).is_empty() {
            return;
        }
        if line.len() >= 2 && line.starts_with('[') && line.ends_with(']') {
            self.sections.push(Section {
                title: line[1..line.len() - 1].to_string(),
                rows: Vec::new(),
            });
            return;
        }
        let Some(section) = self.sections.last_mut() else {
            self.refused = Some(format!(
                "line {} precedes the first [section] header: {}",
                self.lines,
                py_repr_str(line)
            ));
            return;
        };
        let row = match line.split_once('\t') {
            Some((key, value)) => (key.to_string(), Some(value.to_string())),
            None => (line.to_string(), None),
        };
        section.rows.push(row);
    }

    /// The sections in file order, or the reason the text is not a statistics report.
    pub fn finish(self) -> Result<Vec<Section>, String> {
        match self.refused {
            Some(reason) => Err(reason),
            None => Ok(self.sections),
        }
    }
}

/// The labels whose values are set aside: [`MASKED_ROWS`] and the `--mask-row` labels.
pub struct MaskedRows {
    /// Every label once, in the order Python sorts the `str` of each.
    sorted: Vec<Vec<u8>>,
    /// The labels a row key can equal: a key is text, so a label that is not UTF-8 names
    /// no row.
    keys: HashSet<String>,
}

impl MaskedRows {
    pub fn new(mask_rows: &[Vec<u8>]) -> MaskedRows {
        let mut sorted: Vec<Vec<u8>> = MASKED_ROWS
            .iter()
            .map(|label| label.as_bytes().to_vec())
            .chain(mask_rows.iter().cloned())
            .collect();
        sorted.sort_by_cached_key(|label| surrogateescape(label));
        sorted.dedup();
        let keys = sorted
            .iter()
            .filter_map(|label| String::from_utf8(label.clone()).ok())
            .collect();
        MaskedRows { sorted, keys }
    }

    fn contains(&self, key: &str) -> bool {
        self.keys.contains(key)
    }
}

fn text_or_null(value: Option<&str>) -> Json {
    value.map_or(Json::Null, Json::str)
}

/// `Annotation sources` with every path replaced: a token carrying a slash, up to a blank
/// or the `; ` that joins two sources, becomes `[PATH]`.
fn mask_paths(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    let mut token_start = 0;
    let push_token = |out: &mut String, token: &str| {
        out.push_str(if token.contains('/') { PATH } else { token });
    };
    for (at, c) in value.char_indices() {
        if c == ';' || is_py_space(c) {
            push_token(&mut out, &value[token_start..at]);
            out.push(c);
            token_start = at + c.len_utf8();
        }
    }
    push_token(&mut out, &value[token_start..]);
    out
}

/// The value as it compares: `<masked>` for a masked row, paths replaced in
/// `Annotation sources`, verbatim otherwise.
fn mask_value<'a>(key: &str, value: Option<&'a str>, masked: &MaskedRows) -> Option<Cow<'a, str>> {
    if masked.contains(key) {
        return Some(Cow::Borrowed(MASKED));
    }
    match value {
        Some(value) if key == PATH_ROW => Some(Cow::Owned(mask_paths(value))),
        value => value.map(Cow::Borrowed),
    }
}

/// Whether the report's command line switches the frequency check on. `Config.pm:734-749`
/// writes every option the run set by its full name, each as `--<name>` and a value
/// other than 1 after it, so a flag stands as a word of its own and is off before a `0`.
fn frequency_checked(report: &[Section]) -> bool {
    let run = report
        .iter()
        .find(|section| section.title == RUN_STATISTICS);
    let options = run.and_then(|run| run.rows.iter().find(|(key, _)| key == OPTIONS_ROW));
    let Some((_, Some(options))) = options else {
        return false;
    };
    let mut words = options.split(' ').peekable();
    while let Some(word) = words.next() {
        if FREQUENCY_FLAGS.contains(&word) && words.peek() != Some(&"0") {
            return true;
        }
    }
    false
}

/// The key a `--custom` value gives its source. `AnnotationSourceAdaptor.pm:351-363` reads
/// a value that holds `=` as `key=value` items between commas, the last item of a key the
/// one kept and the text of an item up to its second `=` alone, and any other value by
/// position, the key second. A value that names no key is keyed by its file's base name
/// (`AnnotationSource/File.pm:154-158`), which the `.gz` or `.bgz` ending Ensembl VEP asks
/// of a VCF source (ensembl-io `TabixParser.pm:52`) keeps from spelling `Consequence`.
fn custom_key(source: &str) -> Option<&str> {
    if !source.contains('=') {
        return source.split(',').nth(1);
    }
    let mut keys = source.split(',').filter_map(|item| {
        let mut parts = item.split('=');
        (parts.next() == Some("short_name"))
            .then(|| parts.next())
            .flatten()
    });
    keys.next_back()
}

/// Whether the report's command line names a `--custom` source whose key is `Consequence`.
/// `Config.pm:734-749` writes each source as the word `--custom` and its value after it.
/// Its path mask (`:796`) replaces the text from the last `,` or `=` through the directory
/// of a later source given by position, that word included, so the key of such a source is
/// not read, nor a key the replaced text held.
fn custom_source_keyed_consequence(report: &[Section]) -> bool {
    let run = report
        .iter()
        .find(|section| section.title == RUN_STATISTICS);
    let options = run.and_then(|run| run.rows.iter().find(|(key, _)| key == OPTIONS_ROW));
    let Some((_, Some(options))) = options else {
        return false;
    };
    let mut words = options.split(' ');
    while let Some(word) = words.next() {
        if word == CUSTOM_FLAG && words.next().and_then(custom_key) == Some(CONSEQUENCE_KEY) {
            return true;
        }
    }
    false
}

/// A count as `[General statistics]` writes one: decimal digits alone, no leading zero.
fn plain_count(text: &str) -> Option<u64> {
    let digits = !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit());
    if !digits || (text.len() > 1 && text.starts_with('0')) {
        return None;
    }
    text.parse().ok().filter(|count| *count <= MAX_COUNT)
}

/// The `Novel / existing variants` row of a run that counts `existing` of its
/// `processed` variants as existing, as `Stats.pm:727-734` writes it: `-` for a count of
/// none, else both counts, each with `100 * (count / processed)` under `%.1f`.
fn existing_row(processed: u64, existing: u64) -> String {
    if existing == 0 {
        return "-".to_string();
    }
    let novel = processed - existing;
    let percent = |count: u64| 100.0 * (count as f64 / processed as f64);
    format!(
        "{novel} ({:.1}) / {existing} ({:.1})",
        percent(novel),
        percent(existing)
    )
}

/// Whether the two values of `Novel / existing variants` differ only as the documented
/// difference makes them: `expected` counts each of the `processed` variants as existing,
/// and `actual` is the row of a count from none to all of them.
fn existing_rows_agree(
    expected: Option<&str>,
    actual: Option<&str>,
    processed: Option<&str>,
) -> bool {
    let (Some(expected), Some(actual)) = (expected, actual) else {
        return false;
    };
    let Some(processed) = processed.and_then(plain_count) else {
        return false;
    };
    if expected != existing_row(processed, processed) {
        return false;
    }
    // The count of existing variants is the digits after the row's ` / `.
    let existing = match actual.split_once(" / ") {
        Some((_, existing)) => existing.split(' ').next().and_then(plain_count),
        None => Some(0),
    };
    existing.is_some_and(|count| count <= processed && actual == existing_row(processed, count))
}

/// A section's rows as a map that keeps the file's key order.
struct Keyed<'a> {
    keys: Vec<&'a str>,
    values: HashMap<&'a str, Option<&'a str>>,
}

impl<'a> Keyed<'a> {
    /// The section's rows as a map, every key printed twice recorded as a difference; the
    /// first value of a repeated key is the one kept.
    fn new(section: &'a Section, side: &str, differences: &mut Vec<Json>) -> Keyed<'a> {
        let mut keyed = Keyed {
            keys: Vec::with_capacity(section.rows.len()),
            values: HashMap::with_capacity(section.rows.len()),
        };
        for (key, value) in &section.rows {
            if keyed.values.contains_key(key.as_str()) {
                differences.push(Json::object([
                    ("kind", Json::str("duplicate")),
                    ("row", Json::str(key)),
                    ("side", Json::str(side)),
                ]));
                continue;
            }
            keyed.keys.push(key);
            keyed.values.insert(key, value.as_deref());
        }
        keyed
    }

    fn has(&self, key: &str) -> bool {
        self.values.contains_key(key)
    }

    fn value(&self, key: &str) -> Option<&'a str> {
        self.values.get(key).copied().flatten()
    }
}

/// Row by row: the same keys in the same order, the masked values set aside. Under
/// `frequency_checked` a `Novel / existing variants` row of `[General statistics]` that
/// differs as its documented difference does goes to `documented`, not to the differences.
fn compare_tabular(
    expected: &Section,
    actual: &Section,
    masked: &MaskedRows,
    frequency_checked: bool,
    documented: &mut Vec<Documented>,
) -> Vec<Json> {
    let mut differences = Vec::new();
    let e_rows = Keyed::new(expected, "expected", &mut differences);
    let a_rows = Keyed::new(actual, "actual", &mut differences);
    for key in e_rows.keys.iter().filter(|key| !a_rows.has(key)) {
        differences.push(Json::object([
            ("kind", Json::str("missing")),
            ("row", Json::str(key)),
            ("expected", text_or_null(e_rows.value(key))),
        ]));
    }
    for key in a_rows.keys.iter().filter(|key| !e_rows.has(key)) {
        differences.push(Json::object([
            ("kind", Json::str("unexpected")),
            ("row", Json::str(key)),
            ("actual", text_or_null(a_rows.value(key))),
        ]));
    }
    let shared = |rows: &Keyed<'_>, other: &Keyed<'_>| -> Vec<String> {
        let shared = rows.keys.iter().filter(|key| other.has(key));
        shared.map(|key| key.to_string()).collect()
    };
    let shared_e = shared(&e_rows, &a_rows);
    let shared_a = shared(&a_rows, &e_rows);
    if shared_e != shared_a {
        differences.push(Json::object([
            ("kind", Json::str("order")),
            ("expected", Json::strs(shared_e.iter().map(String::as_str))),
            ("actual", Json::strs(shared_a.iter().map(String::as_str))),
        ]));
    }
    for key in &shared_e {
        let e_value = mask_value(key, e_rows.value(key), masked);
        let a_value = mask_value(key, a_rows.value(key), masked);
        if e_value == a_value {
            continue;
        }
        let is_existing_row = expected.title == GENERAL_STATISTICS && key == EXISTING_ROW;
        if frequency_checked
            && is_existing_row
            && existing_rows_agree(
                e_value.as_deref(),
                a_value.as_deref(),
                e_rows.value(PROCESSED_ROW),
            )
        {
            documented.push(Documented {
                class: EXISTING_VARIANTS_CLASS,
                section: expected.title.clone(),
                row: key.clone(),
                expected: e_value.map(Cow::into_owned),
                actual: a_value.map(Cow::into_owned),
            });
            continue;
        }
        differences.push(Json::object([
            ("kind", Json::str("value")),
            ("row", Json::str(key)),
            ("expected", text_or_null(e_value.as_deref())),
            ("actual", text_or_null(a_value.as_deref())),
        ]));
    }
    differences
}

/// A row whose two values differ as a documented difference makes them differ: it is
/// set aside, and the report names it.
struct Documented {
    class: &'static str,
    section: String,
    row: String,
    expected: Option<String>,
    actual: Option<String>,
}

impl Documented {
    fn json(&self) -> Json {
        Json::object([
            ("class", Json::str(self.class)),
            ("section", Json::str(&self.section)),
            ("row", Json::str(&self.row)),
            ("expected", text_or_null(self.expected.as_deref())),
            ("actual", text_or_null(self.actual.as_deref())),
        ])
    }

    /// The row's line on stdout.
    fn line(&self) -> String {
        format!(
            "[{}]: {} set aside as {}",
            self.section, self.row, self.class
        )
    }
}

/// A chart row's count: an integer where `int()` reads the value as one, the text as
/// written where it does not, nothing for a row without a tab or a key the side lacks.
#[derive(PartialEq, Eq)]
enum Count<'a> {
    Absent,
    Int(PyInt),
    Text(&'a str),
}

impl<'a> Count<'a> {
    fn of(value: Option<&'a str>) -> Count<'a> {
        match value {
            None => Count::Absent,
            Some(text) => PyInt::parse(text).map_or(Count::Text(text), Count::Int),
        }
    }

    fn json(&self) -> Json {
        match self {
            Count::Absent => Json::Null,
            Count::Int(value) => Json::Int(value.clone()),
            Count::Text(text) => Json::str(text),
        }
    }
}

/// Whether a key of `[Consequences (all)]` is one the documented difference of a source
/// keyed `Consequence` adds: the expected side counts it, as an integer, the candidate does
/// not hold it, and it names no consequence type.
fn record_name_counted_as_a_consequence(key: &str, expected: &Count<'_>, in_actual: bool) -> bool {
    matches!(expected, Count::Int(_)) && !in_actual && !is_consequence_term(key)
}

/// As a set of `(key, count)` pairs; under a tolerance of at least 1 a key's two counts may
/// differ by up to the tolerance, an absent key counting as 0. Under `keyed_consequence` a
/// key of `[Consequences (all)]` that the expected side alone counts and that names no
/// consequence type goes to `documented`, not to the differences.
fn compare_chart(
    expected: &Section,
    actual: &Section,
    tolerance: &PyInt,
    masked: &MaskedRows,
    keyed_consequence: bool,
    documented: &mut Vec<Documented>,
) -> Vec<Json> {
    let mut differences = Vec::new();
    let e_rows = Keyed::new(expected, "expected", &mut differences);
    let a_rows = Keyed::new(actual, "actual", &mut differences);
    let only_actual = a_rows.keys.iter().filter(|key| !e_rows.has(key));
    for key in e_rows.keys.iter().chain(only_actual) {
        if masked.contains(key) {
            continue;
        }
        let e_count = Count::of(e_rows.value(key));
        let a_count = Count::of(a_rows.value(key));
        if keyed_consequence
            && expected.title == CONSEQUENCES_ALL
            && record_name_counted_as_a_consequence(key, &e_count, a_rows.has(key))
        {
            documented.push(Documented {
                class: CONSEQUENCE_KEY_CLASS,
                section: expected.title.clone(),
                row: key.to_string(),
                expected: e_rows.value(key).map(str::to_string),
                actual: None,
            });
            continue;
        }
        let (kind, within) = if e_rows.has(key) && a_rows.has(key) {
            let close = match (&e_count, &a_count) {
                (Count::Int(e), Count::Int(a)) => e.within(a, tolerance),
                _ => false,
            };
            ("count", e_count == a_count || close)
        } else {
            let (kind, present) = if e_rows.has(key) {
                ("missing", &e_count)
            } else {
                ("unexpected", &a_count)
            };
            let small = matches!(present, Count::Int(count) if count.at_most(tolerance));
            (kind, tolerance.is_positive() && small)
        };
        if !within {
            differences.push(Json::object([
                ("kind", Json::str(kind)),
                ("row", Json::str(key)),
                ("expected", e_count.json()),
                ("actual", a_count.json()),
            ]));
        }
    }
    differences
}

/// The comparison of two statistics reports.
pub struct Report {
    /// One line per section that differs and per section on one side alone; none means
    /// the verdict is green.
    pub reasons: Vec<String>,
    /// One line per row set aside under a documented difference.
    set_aside: Vec<String>,
    /// The JSON report without the two paths.
    comparison: Json,
}

impl Report {
    /// `green` or `red`.
    pub fn verdict(&self) -> &'static str {
        if self.reasons.is_empty() {
            "green"
        } else {
            "red"
        }
    }

    /// The JSON report: the comparison, then the two paths as the command line gave them.
    pub fn json(&self, expected: &[u8], actual: &[u8]) -> Json {
        let Json::Object(mut members) = self.comparison.clone() else {
            unreachable!("the report is an object");
        };
        members.push(("expected".to_string(), Json::Arg(expected.to_vec())));
        members.push(("actual".to_string(), Json::Arg(actual.to_vec())));
        Json::Object(members)
    }

    /// The lines written to stdout: `stats: <verdict>`, then each reason and each row set
    /// aside, indented by two.
    pub fn summary(&self) -> String {
        let mut out = format!("{MODE}: {}\n", self.verdict());
        for line in self.reasons.iter().chain(&self.set_aside) {
            out.push_str("  ");
            out.push_str(line);
            out.push('\n');
        }
        out
    }
}

/// The titles once each, in file order, that the other side carries too.
fn shared_titles<'a>(titles: &[&'a str], other: &[&str]) -> Vec<&'a str> {
    let mut seen = HashSet::new();
    let shared = titles.iter().filter(|title| other.contains(title));
    shared
        .filter(|title| seen.insert(**title))
        .copied()
        .collect()
}

/// Compares two parsed reports. The checks run in this order and the first that fails is
/// the error: a negative tolerance, an expected text that is not a statistics report, an
/// actual text that is not one.
pub fn compare_reports(
    expected: Result<Vec<Section>, String>,
    actual: Result<Vec<Section>, String>,
    tolerance: &PyInt,
    mask_rows: &[Vec<u8>],
) -> Result<Report, String> {
    if tolerance.is_negative() {
        return Err(format!("--tolerance {} is negative", tolerance.decimal()));
    }
    let masked = MaskedRows::new(mask_rows);
    let expected = expected?;
    let actual = actual?;
    let e_titles: Vec<&str> = expected.iter().map(|s| s.title.as_str()).collect();
    let a_titles: Vec<&str> = actual.iter().map(|s| s.title.as_str()).collect();
    let mut differences: Vec<(String, Json)> = Vec::new();
    let mut reasons: Vec<String> = Vec::new();
    for (side, titles) in [("expected", &e_titles), ("actual", &a_titles)] {
        if titles.is_empty() {
            reasons.push(format!("the {side} file carries no [section] header"));
        }
        for title in titles {
            let repeated = format!("section [{title}] appears twice in the {side} file");
            let copies = titles.iter().filter(|t| *t == title).count();
            if copies > 1 && !reasons.contains(&repeated) {
                reasons.push(repeated);
            }
        }
    }
    let mut set_difference = |title: &str, found: Json| {
        match differences.iter_mut().find(|(t, _)| t == title) {
            Some(entry) => entry.1 = found,
            None => differences.push((title.to_string(), found)),
        };
    };
    let one_sided = |kind: &str| Json::List(vec![Json::object([("kind", Json::str(kind))])]);
    for title in e_titles.iter().filter(|t| !a_titles.contains(t)) {
        set_difference(title, one_sided("section_missing"));
        reasons.push(format!("section [{title}] only in the expected file"));
    }
    for title in a_titles.iter().filter(|t| !e_titles.contains(t)) {
        set_difference(title, one_sided("section_unexpected"));
        reasons.push(format!("section [{title}] only in the actual file"));
    }
    let shared_e = shared_titles(&e_titles, &a_titles);
    let shared_a = shared_titles(&a_titles, &e_titles);
    let order_agrees = shared_e == shared_a;
    if !order_agrees {
        reasons.push("the shared sections come in another order in the actual file".to_string());
    }
    // A repeated title is reported above; on either side its first copy is the one compared.
    let mut actual_by_title: HashMap<&str, &Section> = HashMap::new();
    for section in &actual {
        actual_by_title.entry(&section.title).or_insert(section);
    }
    let frequency_checked = frequency_checked(&expected);
    let keyed_consequence = custom_source_keyed_consequence(&expected);
    let mut documented: Vec<Documented> = Vec::new();
    let mut rows_compared = 0u64;
    let mut compared_titles = HashSet::new();
    for section in &expected {
        let Some(other) = actual_by_title.get(section.title.as_str()) else {
            continue;
        };
        if !compared_titles.insert(section.title.as_str()) {
            continue;
        }
        let found = if section.is_chart() {
            compare_chart(
                section,
                other,
                tolerance,
                &masked,
                keyed_consequence,
                &mut documented,
            )
        } else {
            compare_tabular(section, other, &masked, frequency_checked, &mut documented)
        };
        rows_compared += section.rows.len() as u64;
        if !found.is_empty() {
            let noun = if found.len() == 1 {
                "difference"
            } else {
                "differences"
            };
            reasons.push(format!("[{}]: {} {noun}", section.title, found.len()));
            set_difference(&section.title, Json::List(found));
        }
    }
    let verdict = if reasons.is_empty() { "green" } else { "red" };
    let comparison = Json::object([
        ("mode", Json::str(MODE)),
        ("format", Json::str(MODE)),
        ("verdict", Json::str(verdict)),
        ("reasons", Json::strs(reasons.iter().map(String::as_str))),
        ("tolerance", Json::Int(tolerance.clone())),
        (
            "masked_rows",
            Json::List(masked.sorted.iter().cloned().map(Json::Arg).collect()),
        ),
        (
            "sections",
            Json::object([
                ("expected", Json::strs(e_titles.iter().copied())),
                ("actual", Json::strs(a_titles.iter().copied())),
                ("order_agrees", Json::Bool(order_agrees)),
            ]),
        ),
        (
            "compared",
            Json::object([
                (
                    "sections",
                    Json::Int(PyInt::from_u64(shared_e.len() as u64)),
                ),
                ("rows", Json::Int(PyInt::from_u64(rows_compared))),
            ]),
        ),
        (
            "documented",
            Json::List(documented.iter().map(Documented::json).collect()),
        ),
        ("differences", Json::Object(differences)),
    ]);
    Ok(Report {
        reasons,
        set_aside: documented.iter().map(Documented::line).collect(),
        comparison,
    })
}

#[cfg(test)]
#[path = "report_tests.rs"]
mod tests;
