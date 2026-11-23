// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The warnings mode.

use crate::argparse::path_text;
use crate::classes;
use crate::classes::CADD_STRUCTURAL_ALLELE_MATCH;
use crate::classes::FORK_WARNING_PID_PREFIX;
use crate::classes::JSON_MISMATCH_WARNING_NAMES_NO_RECORD;
use crate::classes::REGULATORY_FEATURE_ROWS;
use crate::error::Error;
use crate::error::Result;
use crate::input::read_text;
use crate::input::text_lines;
use crate::pyfmt;
use crate::pyjson::Value;
use indexmap::IndexMap;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::LazyLock;

/// The 0.1.0-to-0.5.0 and the cache-aware warnings-file phrases for --regulatory on a cache
/// without regulatory data.
pub const REGULATORY_NOTICE: &str = "regulatory-feature annotation is not implemented";

pub const MISSING_REGULATORY_NOTICE: &str = "this cache holds no regulatory data";

pub const FORKED_CHILD_RELAY_RE_PATTERN: &str = r"^WARNING:[\s\x1c-\x1f]*\d+ : ";
pub static FORKED_CHILD_RELAY_RE: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(FORKED_CHILD_RELAY_RE_PATTERN).expect("the pattern compiles")
});

pub const WARNING_LINE_RE_PATTERN: &str = r"^WARNING:[\s\x1c-\x1f]*(?:line (\d+) skipped|variant skipped)[\s\x1c-\x1f]*(?:\((.*?)\): )?(.*)(?:\n?\z)";
pub static WARNING_LINE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(WARNING_LINE_RE_PATTERN).expect("the pattern compiles"));

pub const WARNING_OTHER_RE_PATTERN: &str =
    r"^WARNING:[\s\x1c-\x1f]*(.*?)(?:[\s\x1c-\x1f]+(?:on )?line (\d+))?(?:\n?\z)";
pub static WARNING_OTHER_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(WARNING_OTHER_RE_PATTERN).expect("the pattern compiles"));

pub const SKIPPED_LINE_RE_PATTERN: &str =
    r"^Line[\s\x1c-\x1f]+([^\s\x1c-\x1f]+)[\s\x1c-\x1f]*\t(.*)\t(.*)(?:\n?\z)";
pub static SKIPPED_LINE_RE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(SKIPPED_LINE_RE_PATTERN).expect("the pattern compiles"));

/// A path token. `mask_description` reads the same tokens without the pattern.
pub const PATH_RE_PATTERN: &str = r"(?<![^\s\x1c-\x1f])(?:~|\.{1,2})?/[^\s\x1c-\x1f]+";
pub static PATH_RE: LazyLock<fancy_regex::Regex> =
    LazyLock::new(|| fancy_regex::Regex::new(PATH_RE_PATTERN).expect("the pattern compiles"));

/// The Perl run-time warnings a plugin's undefined allele string raises before the plugin fails
/// (Utils/Sequence.pm:1193 and :1252 at release/116), written to the warnings file as they are.
pub const UNDEFINED_ALLELE_WARNING_PATTERN: &str =
    r"Use of uninitialized value \$(?:orig_a_alt|alt) in .*Sequence\.pm line \d+";
pub static UNDEFINED_ALLELE_WARNING: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(UNDEFINED_ALLELE_WARNING_PATTERN).expect("the pattern compiles")
});

/// The description of the transcript-assembly mismatch warning before the record's name
/// (OutputFactory.pm:1727 at release/116.3).
pub const MISMATCH_WARNING: &str = "Transcript-assembly mismatch in";

/// The interpreter's warning of a concatenation that joins an undefined value in
/// `OutputFactory.pm`, bare of a relay prefix; group 1 runs through the statement's line
/// number and so names the statement.
pub const UNDEFINED_CONCATENATION_WARNING_PATTERN: &str = r"^(Use of uninitialized value in concatenation \(\.\) or string at (?:.*/)?OutputFactory\.pm line \d+)[,.]";
pub static UNDEFINED_CONCATENATION_WARNING: LazyLock<regex::Regex> = LazyLock::new(|| {
    regex::Regex::new(UNDEFINED_CONCATENATION_WARNING_PATTERN).expect("the pattern compiles")
});

pub const EXCEPTION_RULE_PATTERN: &str = r"^-{10,}(?: EXCEPTION -{10,})?(?:\n?\z)";
pub static EXCEPTION_RULE: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(EXCEPTION_RULE_PATTERN).expect("the pattern compiles"));

/// The input line a warning names: `int()` of a run of digits, of any length the
/// interpreter converts, held as its decimal text.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct LineNumber(String);

impl LineNumber {
    /// `int(digits)`; a run of digits `int()` refuses ends the run with its `ValueError` text.
    pub fn of(digits: &str) -> LineNumber {
        match pyfmt::int_checked(digits) {
            Ok(decimal) => LineNumber(decimal),
            Err(error) => crate::uncaught!("{error}"),
        }
    }

    /// The number as the report holds it.
    pub fn value(&self) -> Value {
        match self.0.parse() {
            Ok(number) => Value::Int(number),
            Err(_) => Value::BigInt(self.0.clone()),
        }
    }
}

impl Ord for LineNumber {
    /// Numeric order: no number is negative and none has a leading zero.
    fn cmp(&self, other: &LineNumber) -> Ordering {
        (self.0.len(), &self.0).cmp(&(other.0.len(), &other.0))
    }
}

impl PartialOrd for LineNumber {
    fn partial_cmp(&self, other: &LineNumber) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// One entry of a warnings file: the input line it names, if any, and its description.
pub type WarningEntry = (Option<LineNumber>, String);

/// Whether a run of characters without a blank is a path: it opens with `/`, `~/`, `./` or
/// `../` and goes on after that `/`.
fn is_path(token: &str) -> bool {
    ["~", "..", ".", ""].iter().any(|prefix| {
        token
            .strip_prefix(prefix)
            .and_then(|rest| rest.strip_prefix('/'))
            .is_some_and(|rest| !rest.is_empty())
    })
}

/// A warning description with every path token replaced by `<path>`.
///
/// The tokens `PATH_RE` matches are read in one pass over the text: the pattern's engine
/// refuses a text of a million characters without a path.
pub fn mask_description(text: &str) -> String {
    let mut masked = String::with_capacity(text.len());
    let mut rest = text;
    while !rest.is_empty() {
        let blanks = rest
            .find(|c| !pyfmt::char_is_space(c))
            .unwrap_or(rest.len());
        masked.push_str(&rest[..blanks]);
        rest = &rest[blanks..];
        let end = rest.find(pyfmt::char_is_space).unwrap_or(rest.len());
        let token = &rest[..end];
        masked.push_str(if is_path(token) { "<path>" } else { token });
        rest = &rest[end..];
    }
    pyfmt::strip(&masked).to_string()
}

/// The (input line number, description) set of a `_warnings.txt` or a
/// `--skipped_variants_file` output; `BaseVEP.pm:802-863` and
/// `BaseRunner.pm:364-380` are the two line shapes.
pub fn parse_warning_lines(text: &str) -> BTreeSet<WarningEntry> {
    let mut out = BTreeSet::new();
    for line in text_lines(text) {
        if pyfmt::strip(&line).is_empty()
            || line.starts_with("[VEP skipped the following variants from")
        {
            continue;
        }
        if let Some(found) = SKIPPED_LINE_RE.captures(&line) {
            let number = Some(&found[1])
                .filter(|number| pyfmt::is_digit(number))
                .map(LineNumber::of);
            out.insert((number, mask_description(&found[3])));
            continue;
        }
        if let Some(found) = WARNING_LINE_RE.captures(&line) {
            if found.get(1).is_some() || line.starts_with("WARNING: variant skipped") {
                let number = found.get(1).map(|number| LineNumber::of(number.as_str()));
                out.insert((number, mask_description(&found[3])));
                continue;
            }
        }
        if let Some(found) = WARNING_OTHER_RE.captures(&line) {
            let number = found.get(2).map(|number| LineNumber::of(number.as_str()));
            out.insert((number, mask_description(&found[1])));
            continue;
        }
        out.insert((None, mask_description(&line)));
    }
    out
}

/// The reference's warnings text without the lines reporting one plugin's failure: each
/// `WARNING: Plugin '<name>' went wrong` line (`failure` is its text) with the exception
/// block Ensembl VEP writes after it, the dashed `EXCEPTION` header through the closing
/// dashed line (`BaseVEP.pm` writes the exception's own text; a block the next `WARNING`
/// line interrupts ends there), and with `undefined_allele_lines` the Perl warnings about
/// the undefined allele a CADD failure raises, wherever they stand. Returns the text and the
/// number of lines dropped, blocks aside.
pub fn without_plugin_failure(
    text: &str,
    failure: &str,
    undefined_allele_lines: bool,
) -> (String, i64) {
    let rule = |line: &str| EXCEPTION_RULE.is_match(pyfmt::strip(line));
    let lines: Vec<&str> = text.split('\n').collect();
    let mut kept: Vec<&str> = Vec::new();
    let mut dropped = 0;
    let mut i = 0;
    while i < lines.len() {
        let line = lines[i];
        i += 1;
        if line.contains(failure) {
            dropped += 1;
            // An exception block follows a failure whose text the plugin threw: its dashed
            // header, then its lines up to and including the closing dashed rule.
            if i < lines.len() && lines[i].contains("EXCEPTION") && rule(lines[i]) {
                i += 1;
                while i < lines.len() && !lines[i].starts_with("WARNING") {
                    let closing = rule(lines[i]) && !lines[i].contains("EXCEPTION");
                    i += 1;
                    if closing {
                        break;
                    }
                }
            }
            continue;
        }
        if undefined_allele_lines && UNDEFINED_ALLELE_WARNING.is_match(line) {
            dropped += 1;
            continue;
        }
        kept.push(line);
    }
    (kept.join("\n"), dropped)
}

/// The reference's warnings text with the `<pid> : ` a forked parent wrote before a
/// child's first line removed from every line carrying one, and the number of lines changed.
/// A child's own warning then reads as the serial run writes it; a child's interpreter
/// warning, bare after the prefix, stays a line of its own.
pub fn without_forked_child_prefix(text: &str) -> (String, i64) {
    let mut rewritten = Vec::new();
    let mut changed = 0;
    for line in text.split('\n') {
        let stripped = FORKED_CHILD_RELAY_RE.replace(line, "");
        if stripped != line {
            changed += 1;
        }
        rewritten.push(stripped);
    }
    (rewritten.join("\n"), changed)
}

/// The one entry a line of a warnings file reads as, where it reads as one.
fn entry_of(line: &str) -> Option<WarningEntry> {
    let mut entries = parse_warning_lines(line).into_iter();
    let entry = entries.next()?;
    entries.next().is_none().then_some(entry)
}

/// Whether a line is the mismatch warning with no record name after it.
fn names_no_record(line: &str) -> bool {
    entry_of(line)
        .is_some_and(|(number, description)| number.is_none() && description == MISMATCH_WARNING)
}

/// Whether a line is the mismatch warning with a record's name after it.
fn names_a_record(line: &str) -> bool {
    entry_of(line).is_some_and(|(_, description)| {
        description
            .strip_prefix(MISMATCH_WARNING)
            .is_some_and(|name| name.starts_with(pyfmt::char_is_space))
    })
}

/// The statement a bare interpreter line of `UNDEFINED_CONCATENATION_WARNING` names.
fn undefined_concatenation(line: &str) -> Option<&str> {
    let found = UNDEFINED_CONCATENATION_WARNING.captures(line)?;
    found.get(1).map(|statement| statement.as_str())
}

/// The two warnings texts without the lines of `json_mismatch_warning_names_no_record`, and
/// the number of lines left out of the two. The class reads its own shape alone: the
/// reference holds the mismatch warning with no record name and none with a name, and the
/// candidate holds one with a name; any other pair of texts is returned as it is. Left out
/// are the reference's unnamed lines, the candidate's named lines and each bare line of the
/// reference that is the interpreter's warning of the statement that wrote an unnamed line.
/// That statement joins the missing name before it writes the line, so its warning is the
/// line directly before an unnamed line, read here without a relay prefix; an interpreter
/// warning of another statement stays.
pub fn without_unnamed_mismatch_warning(expected: &str, actual: &str) -> (String, String, i64) {
    let reference: Vec<&str> = expected.split('\n').collect();
    let candidate: Vec<&str> = actual.split('\n').collect();
    let unnamed: Vec<usize> = (0..reference.len())
        .filter(|at| names_no_record(reference[*at]))
        .collect();
    let own_shape = !unnamed.is_empty()
        && !reference.iter().any(|line| names_a_record(line))
        && candidate.iter().any(|line| names_a_record(line));
    if !own_shape {
        return (expected.to_string(), actual.to_string(), 0);
    }
    let statements: BTreeSet<String> = unnamed
        .iter()
        .filter_map(|at| at.checked_sub(1))
        .filter_map(|before| {
            let bare = FORKED_CHILD_RELAY_RE.replace(reference[before], "");
            undefined_concatenation(&bare).map(str::to_string)
        })
        .collect();
    let mut left_out = 0;
    let mut kept = |of_the_class: bool| {
        left_out += i64::from(of_the_class);
        !of_the_class
    };
    let reference: Vec<&str> = reference
        .into_iter()
        .filter(|line| {
            let written_by_the_statement = undefined_concatenation(line)
                .is_some_and(|statement| statements.contains(statement));
            kept(names_no_record(line) || written_by_the_statement)
        })
        .collect();
    let candidate: Vec<&str> = candidate
        .into_iter()
        .filter(|line| kept(names_a_record(line)))
        .collect();
    (reference.join("\n"), candidate.join("\n"), left_out)
}

/// The entries as the report lists them: the numbered ones first by number then by
/// description, then the others by description.
fn rows<'a>(entries: impl Iterator<Item = &'a WarningEntry>) -> Value {
    let mut entries: Vec<&WarningEntry> = entries.collect();
    entries.sort_by(|(a, a_text), (b, b_text)| {
        a.is_none()
            .cmp(&b.is_none())
            .then_with(|| a.cmp(b))
            .then_with(|| a_text.cmp(b_text))
    });
    let row = |(number, description): &WarningEntry| {
        Value::Dict(IndexMap::from([
            (
                "line".to_string(),
                number.as_ref().map_or(Value::Null, LineNumber::value),
            ),
            ("description".to_string(), Value::Str(description.clone())),
        ]))
    };
    Value::List(entries.into_iter().map(row).collect())
}

fn count(entries: usize) -> Value {
    Value::Int(i64::try_from(entries).unwrap_or(i64::MAX))
}

/// The two files as one-sided entry sets. Under `regulatory_feature_rows` the candidate's
/// notice that it annotates no regulatory feature is required and left out of the sets; under
/// `fork_warning_pid_prefix` the `<pid> : ` before a relayed line of the reference is
/// removed; under a class of `PLUGIN_WARNING_CLASSES` the reference's lines naming that
/// plugin's failure are left out; under `json_mismatch_warning_names_no_record` the
/// reference's mismatch warning that names no record and the candidate's that name one are
/// left out (`without_unnamed_mismatch_warning`); each line changed or left out is counted
/// under its class in `masks_applied`.
pub fn compare_warning_files(
    expected_path: &Path,
    actual_path: &Path,
    masks: &[&str],
) -> Result<Value> {
    let mut expected_text = read_text(expected_path)?;
    let mut applied: IndexMap<String, Value> = IndexMap::new();
    if masks.contains(&FORK_WARNING_PID_PREFIX) {
        let (text, relayed) = without_forked_child_prefix(&expected_text);
        expected_text = text;
        if relayed != 0 {
            applied.insert(FORK_WARNING_PID_PREFIX.to_string(), Value::Int(relayed));
        }
    }
    for (klass, failure) in classes::tables().plugin_warning_classes() {
        if masks.contains(klass) {
            let (text, dropped) = without_plugin_failure(
                &expected_text,
                failure,
                *klass == CADD_STRUCTURAL_ALLELE_MATCH,
            );
            expected_text = text;
            if dropped != 0 {
                applied.insert(klass.to_string(), Value::Int(dropped));
            }
        }
    }
    let mut actual_text = read_text(actual_path)?;
    if masks.contains(&JSON_MISMATCH_WARNING_NAMES_NO_RECORD) {
        let (reference, candidate, left_out) =
            without_unnamed_mismatch_warning(&expected_text, &actual_text);
        expected_text = reference;
        actual_text = candidate;
        if left_out != 0 {
            applied.insert(
                JSON_MISMATCH_WARNING_NAMES_NO_RECORD.to_string(),
                Value::Int(left_out),
            );
        }
    }
    let expected = parse_warning_lines(&expected_text);
    let mut actual = parse_warning_lines(&actual_text);
    let mut notice_missing = false;
    if masks.contains(&REGULATORY_FEATURE_ROWS) {
        let entries = actual.len();
        actual.retain(|(_, description)| {
            !description.contains(REGULATORY_NOTICE)
                && !description.contains(MISSING_REGULATORY_NOTICE)
        });
        notice_missing = actual.len() == entries;
    }
    Ok(Value::Dict(IndexMap::from([
        ("expected".to_string(), Value::Str(path_text(expected_path))),
        ("actual".to_string(), Value::Str(path_text(actual_path))),
        ("expected_entries".to_string(), count(expected.len())),
        ("actual_entries".to_string(), count(actual.len())),
        (
            "shared".to_string(),
            count(expected.intersection(&actual).count()),
        ),
        (
            "expected_only".to_string(),
            rows(expected.difference(&actual)),
        ),
        (
            "actual_only".to_string(),
            rows(actual.difference(&expected)),
        ),
        (
            "regulatory_notice_missing".to_string(),
            Value::Bool(notice_missing),
        ),
        ("masks_applied".to_string(), Value::Dict(applied)),
    ])))
}

/// The warnings-mode report: the `_warnings.txt` pair and, when given, the
/// skipped-variants pair, each as one-sided sets; `masks` applies to the warnings pair.
pub fn compare_warnings(
    expected: &Path,
    actual: &Path,
    expected_skipped: Option<&Path>,
    actual_skipped: Option<&Path>,
    masks: &[&str],
) -> Result<Value> {
    let mut sections = vec![("warnings", compare_warning_files(expected, actual, masks)?)];
    match (expected_skipped, actual_skipped) {
        (None, None) => {}
        (Some(expected_skipped), Some(actual_skipped)) => sections.push((
            "skipped",
            compare_warning_files(expected_skipped, actual_skipped, &[])?,
        )),
        _ => {
            return Err(Error::comparator(
                "--expected-skipped and --actual-skipped go together",
            ))
        }
    }
    let mut reasons = Vec::new();
    for (name, section) in &sections {
        let listed = |key: &str| {
            section
                .get(key)
                .and_then(Value::as_list)
                .map_or(0, <[_]>::len)
        };
        for (key, side) in [("expected_only", "expected"), ("actual_only", "actual")] {
            if listed(key) != 0 {
                reasons.push(Value::Str(format!(
                    "{name}: {} entries only in the {side} file",
                    listed(key)
                )));
            }
        }
        if section.get("regulatory_notice_missing") == Some(&Value::Bool(true)) {
            reasons.push(Value::Str(format!(
                "{name}: the regulatory notice is missing from the actual file"
            )));
        }
    }
    let verdict = if reasons.is_empty() { "green" } else { "red" };
    let mut report = IndexMap::from([
        ("mode".to_string(), Value::Str("warnings".to_string())),
        ("format".to_string(), Value::Str("warnings".to_string())),
        ("verdict".to_string(), Value::Str(verdict.to_string())),
        ("reasons".to_string(), Value::List(reasons)),
    ]);
    for (name, section) in sections {
        report.insert(name.to_string(), section);
    }
    Ok(Value::Dict(report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Uncaught;

    fn entry(number: Option<&str>, description: &str) -> WarningEntry {
        (number.map(LineNumber::of), description.to_string())
    }

    fn parsed(text: &str) -> Vec<WarningEntry> {
        parse_warning_lines(text).into_iter().collect()
    }

    #[test]
    fn a_path_token_is_masked_whole() {
        for path in ["/x", "//", "..//x", "././x", "~//", "~/a", "../b"] {
            assert_eq!(mask_description(path), "<path>", "{path}");
        }
        for other in [
            "/", "~/", "./", "../", ".../a", "a/b", "(/x)", "at:/x", "~x/y",
        ] {
            assert_eq!(mask_description(other), other, "{other}");
        }
        // A token opens after a blank of the interpreter's set, and after no other character.
        assert_eq!(mask_description("a\u{1c}/x"), "a\u{1c}<path>");
        assert_eq!(mask_description("a\u{a0}/x"), "a\u{a0}<path>");
        assert_eq!(mask_description("a\u{200b}/x"), "a\u{200b}/x");
        assert_eq!(
            mask_description("  see /a/b and ~/c, not ../ or x/y  "),
            "see <path> and <path> not ../ or x/y"
        );
        assert_eq!(mask_description(""), "");
    }

    #[test]
    fn the_path_scan_reads_the_tokens_the_pattern_matches() {
        let prefixes = ["", "~", ".", "..", "...", "a", "~.", "/"];
        let rests = ["", "/", "/x", "//", "/x/y", "x", "/é"];
        for prefix in prefixes {
            for rest in rests {
                let token = format!("{prefix}{rest}");
                let text = format!("a {token}\u{1f}b");
                let by_pattern = PATH_RE.replace_all(&text, "<path>");
                assert_eq!(mask_description(&text), by_pattern, "{token}");
            }
        }
    }

    #[test]
    fn a_line_of_millions_of_characters_is_masked() {
        let long = "word ".repeat(600_000);
        let masked = mask_description(&format!("{long}/a/b"));
        assert!(masked.ends_with("word <path>"));
        assert_eq!(masked.len(), long.len() + "<path>".len());
    }

    #[test]
    fn the_patterns_read_the_interpreters_blanks_and_digits() {
        let blank = regex::Regex::new(r"^[\s\x1c-\x1f]$").expect("the pattern compiles");
        let digit = regex::Regex::new(r"^\d$").expect("the pattern compiles");
        let mut buffer = [0u8; 4];
        for c in (0..=0x10ffffu32).filter_map(char::from_u32) {
            let text: &str = c.encode_utf8(&mut buffer);
            assert_eq!(
                blank.is_match(text),
                pyfmt::char_is_space(c),
                "U+{:04X}",
                c as u32
            );
            assert_eq!(
                digit.is_match(text),
                pyfmt::char_is_decimal(c),
                "U+{:04X}",
                c as u32
            );
        }
    }

    #[test]
    fn line_numbers_order_by_value_at_any_length() {
        let mut numbers: Vec<LineNumber> = ["100", "20", "99999999999999999999", "007", "0"]
            .into_iter()
            .map(LineNumber::of)
            .collect();
        numbers.sort();
        let values: Vec<Value> = numbers.iter().map(LineNumber::value).collect();
        assert_eq!(
            values,
            [
                Value::Int(0),
                Value::Int(7),
                Value::Int(20),
                Value::Int(100),
                Value::BigInt("99999999999999999999".to_string()),
            ]
        );
        // A digit of another script is a digit.
        assert_eq!(LineNumber::of("\u{663}").value(), Value::Int(3));
    }

    #[test]
    fn a_warning_line_gives_its_input_line_and_its_description() {
        for (line, number, description) in [
            ("WARNING: line 13 skipped: text", Some("13"), ": text"),
            (
                "WARNING:variant skipped (x): t",
                None,
                "variant skipped (x): t",
            ),
            ("WARNING: on line 9", Some("9"), "on"),
            ("WARNING: line 8", None, "line 8"),
            ("WARNING: line 007 skipped (a): b): c", Some("7"), "b): c"),
            ("WARNING: line \u{663} skipped (): x", Some("3"), "x"),
            ("Line abc\tx\ty", None, "y"),
            ("Line 5\tonly one tab", None, "Line 5\tonly one tab"),
            (
                "Line 12   \t21\t5\tid\ttext /a/b here",
                Some("12"),
                "text <path> here",
            ),
            (
                "WARNING: variant skipped (21 5 x): no /data/x.fa",
                None,
                "no <path>",
            ),
            (
                "WARNING: Unable to open ./two.fa for reading on line 12345678901234567890123",
                Some("12345678901234567890123"),
                "Unable to open <path> for reading",
            ),
            ("plain ~/x text", None, "plain <path> text"),
            (
                "WARNING: line 99999999999999999999 skipped (): big",
                Some("99999999999999999999"),
                "big",
            ),
        ] {
            assert_eq!(parsed(line), [entry(number, description)], "{line}");
        }
        assert_eq!(
            parsed("[VEP skipped the following variants from /x]\n   \n\n"),
            []
        );
        // Equal entries are one entry, and a final carriage return is no part of a line.
        assert_eq!(
            parsed("WARNING: b\r\nWARNING: a\nWARNING: b\n"),
            [entry(None, "a"), entry(None, "b")]
        );
    }

    #[test]
    fn a_run_of_digits_the_interpreter_refuses_ends_the_run() {
        let what = |text: String| {
            let payload = std::panic::catch_unwind(|| parse_warning_lines(&text))
                .expect_err("the call ends the run");
            let Uncaught(what) = payload
                .downcast_ref::<Uncaught>()
                .expect("the run ends as an uncaught exception");
            what.clone()
        };
        assert_eq!(
            what("Line \u{b2}\tx\ty".to_string()),
            "invalid literal for int() with base 10: '\u{b2}'"
        );
        assert_eq!(
            what(format!("WARNING: line {} skipped (): x", "1".repeat(4301))),
            "Exceeds the limit (4300 digits) for integer string conversion: value has 4301 \
             digits; use sys.set_int_max_str_digits() to increase the limit"
        );
    }

    #[test]
    fn a_relayed_prefix_is_removed_once_a_line() {
        assert_eq!(
            without_forked_child_prefix("WARNING: 11 : WARNING: a\nWARNING: b\n"),
            ("WARNING: a\nWARNING: b\n".to_string(), 1)
        );
        assert_eq!(
            without_forked_child_prefix(
                "WARNING:11 : x\nWARNING: 12 :y\nWARNING: \u{661}\u{661} : z"
            ),
            ("x\nWARNING: 12 :y\nz".to_string(), 2)
        );
    }

    /// The interpreter's warning of the statement at `line` of `file`.
    fn undefined_at(file: &str, line: u32) -> String {
        format!(
            "Use of uninitialized value in concatenation (.) or string at {file} line {line}, \
             <__ANONIO__> line 4."
        )
    }

    #[test]
    fn the_unnamed_mismatch_warning_goes_with_its_statements_interpreter_lines() {
        let factory = "/src/modules/Bio/EnsEMBL/VEP/OutputFactory.pm";
        let unnamed = "WARNING: Transcript-assembly mismatch in ";
        let named = "WARNING: Transcript-assembly mismatch in rs1\nWARNING: other\n";
        let own = undefined_at(factory, 1727);
        let other_line = undefined_at(factory, 1800);
        // One child's block once its relay prefix is off: the statement's warning before the
        // unnamed line and after it, and another statement's warning.
        let block = [&own, unnamed, &own, &other_line, ""].join("\n");
        assert_eq!(
            without_unnamed_mismatch_warning(&block, named),
            (format!("{other_line}\n"), "WARNING: other\n".to_string(), 4)
        );
        // A relayed line names the statement through its prefix and stays; a warning that
        // ends at the statement's line number is the statement's too.
        let ends_there = format!(
            "{} line 1718.",
            own.split(" line ").next().expect("a place")
        );
        let relayed = format!("WARNING: 11 : {ends_there}");
        let block = [&relayed, unnamed, &ends_there].join("\n");
        assert_eq!(
            without_unnamed_mismatch_warning(&block, named),
            (relayed.clone(), "WARNING: other\n".to_string(), 3)
        );
        // The statement is read from the line directly before an unnamed line and from a
        // file of the factory's name alone: any other interpreter line stays.
        let elsewhere = undefined_at("/src/modules/Bio/EnsEMBL/VEP/Other.pm", 5);
        for block in [
            [unnamed, &own].join("\n"),
            [&own, "WARNING: between", unnamed].join("\n"),
            [&elsewhere, unnamed].join("\n"),
        ] {
            let (reference, candidate, left_out) = without_unnamed_mismatch_warning(&block, named);
            assert_eq!(reference.split('\n').count(), block.split('\n').count() - 1);
            assert!(!reference.contains(unnamed), "{block}");
            assert_eq!((candidate.as_str(), left_out), ("WARNING: other\n", 2));
        }
    }

    #[test]
    fn the_unnamed_mismatch_warning_is_left_out_of_its_own_shape_alone() {
        let unnamed = "WARNING: Transcript-assembly mismatch in \n";
        let named = "WARNING: Transcript-assembly mismatch in rs1\n";
        for (reference, candidate) in [
            // No unnamed line in the reference.
            ("", named),
            (named, named),
            // No named line in the candidate.
            (unnamed, ""),
            (unnamed, unnamed),
            (unnamed, "WARNING: Transcript-assembly mismatch inside rs1\n"),
            // A named line in the reference beside the unnamed one.
            ("WARNING: Transcript-assembly mismatch in rs2\nWARNING: Transcript-assembly mismatch in \n", named),
        ] {
            assert_eq!(
                without_unnamed_mismatch_warning(reference, candidate),
                (reference.to_string(), candidate.to_string(), 0),
                "{reference}{candidate}"
            );
        }
        // A line is read as the comparison reads it: without its blanks, and under a line
        // number where one closes it.
        assert!(names_no_record("WARNING:Transcript-assembly mismatch in"));
        assert!(!names_no_record(
            "WARNING: Transcript-assembly mismatch in line 7"
        ));
        assert!(names_a_record(
            "WARNING: Transcript-assembly mismatch in rs1 on line 7"
        ));
        assert!(!names_a_record(unnamed));
        assert_eq!(
            without_unnamed_mismatch_warning(unnamed, named),
            (String::new(), String::new(), 2)
        );
    }

    #[test]
    fn a_failure_line_goes_with_its_exception_block() {
        // The header and the closing rule are read without the blanks around them.
        assert_eq!(
            without_plugin_failure(
                "a\nfail here\n  ---------- EXCEPTION ----------  \nMSG\n-----------\nkept\n\
                 Use of uninitialized value $alt in x Sequence.pm line 3\n",
                "fail",
                true
            ),
            ("a\nkept\n".to_string(), 2)
        );
        // A block without its closing rule runs to the end of the text.
        assert_eq!(
            without_plugin_failure(
                "fail\n---------- EXCEPTION ----------\nMSG\nlast",
                "fail",
                false
            ),
            (String::new(), 1)
        );
        // Without the switch the undefined-allele line stays.
        assert_eq!(
            without_plugin_failure(
                "Use of uninitialized value $alt in Sequence.pm line 3",
                "fail",
                false
            ),
            (
                "Use of uninitialized value $alt in Sequence.pm line 3".to_string(),
                0
            )
        );
    }
}
