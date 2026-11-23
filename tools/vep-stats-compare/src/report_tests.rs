// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The parser on an Ensembl VEP 116.3 statistics report, the masks, the set rule of the
//! chart sections, the tolerance, the verdict and the report's shape. The candidate texts
//! are the reference report with one change each, plus a candidate-shaped report whose
//! volatile rows and tie orders differ as a second engine's would.

use super::*;
use crate::input::parse_text;
use std::collections::HashMap;

/// The statistics report Ensembl VEP 116.3 wrote for the `GRCh38-statistics` corpus.
const REFERENCE: &str = include_str!("../tests/fixtures/reference.summary.log");
/// The same run under the prediction flags, which add two charts.
const PREDICTIONS: &str = include_str!("../tests/fixtures/predictions.summary.log");

/// `Stats.pm:806-831` writes the three tabular sections, then the charts of
/// `generate_chart_data` (`Stats.pm:505-609`) in its order; the fixture run carries no
/// prediction flag and no FASTA, so it has neither the two prediction charts nor a
/// distribution chart.
const FIXTURE_SECTIONS: [&str; 9] = [
    "VEP run statistics",
    "Data version",
    "General statistics",
    "Variant classes",
    "Consequences (most severe)",
    "Consequences (all)",
    "Coding consequences",
    "Variants by chromosome",
    "Position in protein",
];

/// The volatile rows as a second run of another engine writes them: its own command line
/// (no `--database 0`, a `--json_cache` directory), UTC times, its own paths.
const CANDIDATE_VOLATILE_ROWS: [(&str, &str); 7] = [
    (
        "Annotation sources",
        "Cache: /home/user/vep-rs/tests/golden/116/GRCh38-statistics/json_cache",
    ),
    (
        "Command line options",
        "vep --assembly GRCh38 --buffer_size 5000 --force_overwrite --input_file [PATH]/variants.vcf --json_cache [PATH]/json_cache --offline --output_file [PATH]/default.txt --quiet --species homo_sapiens --stats_text",
    ),
    ("Start time", "2030-01-01 00:00:00"),
    ("End time", "2030-01-01 00:00:01"),
    ("Run time", "1 seconds"),
    (
        "Input file",
        "/home/user/vep-rs/tests/golden/116/GRCh38-statistics/variants.vcf",
    ),
    ("Output file", "/data/out/GRCh38-statistics-default.txt"),
];

fn header(line: &str) -> Option<&str> {
    line.strip_prefix('[')?.strip_suffix(']')
}

/// The text with the value of row `key` of `[section]` replaced.
fn with_row(text: &str, section: &str, key: &str, value: &str) -> String {
    let mut out = Vec::new();
    let mut current = None;
    let mut replaced = false;
    for line in text.lines() {
        if let Some(title) = header(line) {
            current = Some(title);
        } else if current == Some(section) && line.split('\t').next() == Some(key) {
            out.push(format!("{key}\t{value}"));
            replaced = true;
            continue;
        }
        out.push(line.to_string());
    }
    assert!(replaced, "no row {key:?} in [{section}]");
    out.join("\n") + "\n"
}

fn section_lines<'a>(text: &'a str, section: &str) -> Vec<&'a str> {
    let mut lines = Vec::new();
    let mut current = None;
    for line in text.lines() {
        if let Some(title) = header(line) {
            current = Some(title);
        } else if current == Some(section) && !line.is_empty() {
            lines.push(line);
        }
    }
    lines
}

/// The text with the rows of `[section]` replaced by `rows`.
fn with_section_lines(text: &str, section: &str, rows: &[&str]) -> String {
    let mut out: Vec<&str> = Vec::new();
    let mut current = None;
    for line in text.lines() {
        if let Some(title) = header(line) {
            current = Some(title);
            out.push(line);
            if current == Some(section) {
                out.extend(rows);
            }
            continue;
        }
        if current == Some(section) && !line.is_empty() {
            continue;
        }
        out.push(line);
    }
    out.join("\n") + "\n"
}

/// The reference as a second engine would write it: the volatile rows its own, the
/// equal-count ties of `[Variant classes]` in key order.
fn candidate_shaped(text: &str) -> String {
    let mut text = text.to_string();
    for (key, value) in CANDIDATE_VOLATILE_ROWS {
        text = with_row(&text, RUN_STATISTICS, key, value);
    }
    let classes = section_lines(&text, "Variant classes");
    assert_eq!(classes[1..3], ["sequence_alteration\t7", "indel\t7"]);
    let mut swapped = classes.clone();
    swapped.swap(1, 2);
    with_section_lines(&text, "Variant classes", &swapped)
}

fn compare_with(
    expected: &str,
    actual: &str,
    tolerance: &str,
    mask_rows: &[&str],
) -> Result<Report, String> {
    let mask_rows: Vec<Vec<u8>> = mask_rows.iter().map(|m| m.as_bytes().to_vec()).collect();
    let tolerance = PyInt::parse(tolerance).expect("an integer tolerance");
    compare_reports(
        parse_text(expected),
        parse_text(actual),
        &tolerance,
        &mask_rows,
    )
}

fn compare(expected: &str, actual: &str) -> Report {
    compare_with(expected, actual, "0", &[]).unwrap()
}

/// One member of the report, by the path of keys to it.
fn member<'a>(report: &'a Report, path: &[&str]) -> &'a Json {
    path.iter().fold(&report.comparison, |value, key| {
        value
            .get(key)
            .unwrap_or_else(|| panic!("no member {key:?}"))
    })
}

fn differences<'a>(report: &'a Report, section: &str) -> &'a Json {
    member(report, &["differences", section])
}

fn int(value: u64) -> Json {
    Json::Int(PyInt::from_u64(value))
}

fn keys(section: &Section) -> Vec<&str> {
    section.rows.iter().map(|(key, _)| key.as_str()).collect()
}

fn rows(section: &Section) -> HashMap<&str, Option<&str>> {
    let rows = section.rows.iter();
    rows.map(|(key, value)| (key.as_str(), value.as_deref()))
        .collect()
}

#[test]
fn sections_parse_in_order() {
    let sections = parse_text(REFERENCE).unwrap();
    let titles: Vec<&str> = sections.iter().map(|s| s.title.as_str()).collect();
    assert_eq!(titles, FIXTURE_SECTIONS);
    let charts: Vec<bool> = sections.iter().map(Section::is_chart).collect();
    assert_eq!(
        charts,
        [false, false, false, true, true, true, true, true, true]
    );
    let run_stats = &sections[0];
    assert_eq!(
        keys(run_stats),
        [
            "VEP version (API)",
            "Annotation sources",
            "Species",
            "Command line options",
            "Start time",
            "End time",
            "Run time",
            "Input file",
            "Output file",
        ]
    );
    assert_eq!(rows(run_stats)["VEP version (API)"], Some(" 116.3 (116)"));
    assert_eq!(rows(&sections[2])["Lines of input read"], Some("340"));
    assert_eq!(
        rows(&sections[3]),
        HashMap::from([
            ("duplication", Some("1")),
            ("sequence_alteration", Some("7")),
            ("indel", Some("7")),
            ("substitution", Some("9")),
            ("deletion", Some("35")),
            ("insertion", Some("37")),
            ("SNV", Some("243")),
        ])
    );
    // A row without a tab (a [General statistics] value undefined) keeps no value; text
    // before the first header is not a statistics report.
    assert_eq!(
        parse_text("[General statistics]\nOverlapped genes\n").unwrap()[0].rows,
        [("Overlapped genes".to_string(), None)]
    );
    assert_eq!(
        parse_text("Lines of input read\t3\n[General statistics]\n"),
        Err("line 1 precedes the first [section] header: 'Lines of input read\\t3'".to_string())
    );
}

#[test]
fn masked_rows_never_differ() {
    let candidate = candidate_shaped(REFERENCE);
    assert_ne!(candidate, REFERENCE);
    let report = compare(REFERENCE, &candidate);
    assert_eq!(report.verdict(), "green", "{:?}", report.reasons);
    assert_eq!(member(&report, &["differences"]), &Json::Object(Vec::new()));
    assert_eq!(
        member(&report, &["compared"]),
        &Json::object([
            ("sections", int(9)),
            ("rows", int(9 + 13 + 7 + 7 + 22 + 26 + 9 + 2 + 10)),
        ])
    );
    let Json::List(masked_rows) = member(&report, &["masked_rows"]) else {
        panic!("masked_rows is a list");
    };
    for label in MASKED_ROWS {
        assert!(masked_rows.contains(&Json::Arg(label.as_bytes().to_vec())));
    }
    // The path mask leaves the structure of Annotation sources to compare: a second source
    // on one side alone is a difference, the directory alone never is.
    let masked = with_row(
        REFERENCE,
        RUN_STATISTICS,
        "Annotation sources",
        "Cache: [PATH]/116_GRCh38",
    );
    assert_eq!(compare(REFERENCE, &masked).verdict(), "green");
    let custom = with_row(
        REFERENCE,
        RUN_STATISTICS,
        "Annotation sources",
        "Cache: /work/cache/homo_sapiens/116_GRCh38; Custom: /work/corpus/custom/sites.vcf.gz (vcf)",
    );
    let report = compare(REFERENCE, &custom);
    assert_eq!(report.verdict(), "red");
    assert_eq!(
        differences(&report, "VEP run statistics"),
        &Json::List(vec![Json::object([
            ("kind", Json::str("value")),
            ("row", Json::str("Annotation sources")),
            ("expected", Json::str("Cache: [PATH]")),
            ("actual", Json::str("Cache: [PATH]; Custom: [PATH] (vcf)")),
        ])])
    );
    // The rows that are not masked compare verbatim.
    let species = with_row(REFERENCE, RUN_STATISTICS, "Species", "mus_musculus");
    assert_eq!(
        differences(&compare(REFERENCE, &species), "VEP run statistics"),
        &Json::List(vec![Json::object([
            ("kind", Json::str("value")),
            ("row", Json::str("Species")),
            ("expected", Json::str("homo_sapiens")),
            ("actual", Json::str("mus_musculus")),
        ])])
    );
    let version = with_row(
        REFERENCE,
        RUN_STATISTICS,
        "VEP version (API)",
        "116.3 (116)",
    );
    assert_eq!(compare(REFERENCE, &version).verdict(), "red");
    // --mask-row sets one more label aside, in whichever section it occurs.
    let genes = with_row(REFERENCE, GENERAL_STATISTICS, "Overlapped genes", "97");
    assert_eq!(compare(REFERENCE, &genes).verdict(), "red");
    let masked = compare_with(REFERENCE, &genes, "0", &["Overlapped genes"]).unwrap();
    assert_eq!(masked.verdict(), "green");
}

#[test]
fn chart_sections_compare_as_sets() {
    let rows = section_lines(REFERENCE, "Consequences (most severe)");
    let mut swapped_rows = rows.clone();
    swapped_rows.swap(0, 1);
    let swapped = with_section_lines(REFERENCE, "Consequences (most severe)", &swapped_rows);
    assert_ne!(swapped, REFERENCE);
    let report = compare(REFERENCE, &swapped);
    assert_eq!(report.verdict(), "green", "{:?}", report.reasons);
    // Rows of a tabular section keep their order: Stats.pm writes them in a fixed sequence.
    let run_rows = section_lines(REFERENCE, RUN_STATISTICS);
    let mut reordered_rows = run_rows.clone();
    reordered_rows.swap(0, 2);
    let reordered = with_section_lines(REFERENCE, RUN_STATISTICS, &reordered_rows);
    let report = compare(REFERENCE, &reordered);
    assert_eq!(report.verdict(), "red");
    let parsed = parse_text(REFERENCE).unwrap();
    let later = &keys(&parsed[0])[3..];
    let in_order = |first: [&str; 3]| Json::strs(first.into_iter().chain(later.iter().copied()));
    assert_eq!(
        differences(&report, "VEP run statistics"),
        &Json::List(vec![Json::object([
            ("kind", Json::str("order")),
            (
                "expected",
                in_order(["VEP version (API)", "Annotation sources", "Species"])
            ),
            (
                "actual",
                in_order(["Species", "Annotation sources", "VEP version (API)"])
            ),
        ])])
    );
    // A key on one side alone is a set difference, whatever its count.
    let classes = section_lines(REFERENCE, "Variant classes");
    let with_class = |row: &'static str| {
        let rows: Vec<&str> = classes.iter().copied().chain([row]).collect();
        with_section_lines(REFERENCE, "Variant classes", &rows)
    };
    let report = compare(REFERENCE, &with_class("inversion\t0"));
    assert_eq!(
        differences(&report, "Variant classes"),
        &Json::List(vec![Json::object([
            ("kind", Json::str("unexpected")),
            ("row", Json::str("inversion")),
            ("expected", Json::Null),
            ("actual", int(0)),
        ])])
    );
    assert_eq!(report.reasons, ["[Variant classes]: 1 difference"]);
    // A key printed twice is reported, never merged.
    let report = compare(REFERENCE, &with_class("SNV\t1"));
    assert_eq!(
        differences(&report, "Variant classes"),
        &Json::List(vec![Json::object([
            ("kind", Json::str("duplicate")),
            ("row", Json::str("SNV")),
            ("side", Json::str("actual")),
        ])])
    );
    // In a tabular section the first value of a repeated key is the one compared, and the
    // repeat adds no order difference.
    let rows = section_lines(REFERENCE, GENERAL_STATISTICS);
    let mut repeated_rows = rows[..5].to_vec();
    repeated_rows.push("Overlapped genes\t1");
    repeated_rows.extend(&rows[5..]);
    let repeated = with_section_lines(REFERENCE, GENERAL_STATISTICS, &repeated_rows);
    let report = compare(REFERENCE, &repeated);
    assert_eq!(
        differences(&report, "General statistics"),
        &Json::List(vec![Json::object([
            ("kind", Json::str("duplicate")),
            ("row", Json::str("Overlapped genes")),
            ("side", Json::str("actual")),
        ])])
    );
}

#[test]
fn tolerance_applies_per_chart_section() {
    let verdict = |actual: &str, tolerance: &str| {
        let report = compare_with(REFERENCE, actual, tolerance, &[]).unwrap();
        report.verdict()
    };
    let off_by_one = with_row(REFERENCE, "Consequences (all)", "intron_variant", "964");
    let red = compare(REFERENCE, &off_by_one);
    assert_eq!(red.verdict(), "red");
    assert_eq!(
        differences(&red, "Consequences (all)"),
        &Json::List(vec![Json::object([
            ("kind", Json::str("count")),
            ("row", Json::str("intron_variant")),
            ("expected", int(963)),
            ("actual", int(964)),
        ])])
    );
    let green = compare_with(REFERENCE, &off_by_one, "1", &[]).unwrap();
    assert_eq!(green.verdict(), "green", "{:?}", green.reasons);
    assert_eq!(member(&green, &["tolerance"]), &int(1));
    // The tolerance is a count bound within each chart section, not a budget across them:
    // two charts each off by one pass at 1, and one chart off by two fails at 1.
    let two_charts = with_row(
        &off_by_one,
        "Coding consequences",
        "missense_variant",
        "362",
    );
    assert_eq!(verdict(&two_charts, "1"), "green");
    let off_by_two = with_row(REFERENCE, "Consequences (all)", "intron_variant", "965");
    assert_eq!(verdict(&off_by_two, "1"), "red");
    assert_eq!(verdict(&off_by_two, "2"), "green");
    // Under a tolerance a key on one side alone counts as 0 on the other.
    let mut classes = section_lines(REFERENCE, "Variant classes");
    classes.push("inversion\t1");
    let rare = with_section_lines(REFERENCE, "Variant classes", &classes);
    assert_eq!(verdict(&rare, "0"), "red");
    assert_eq!(verdict(&rare, "1"), "green");
    // The tabular sections stay exact under any tolerance.
    let lines = with_row(REFERENCE, GENERAL_STATISTICS, "Lines of input read", "341");
    assert_eq!(verdict(&lines, "5"), "red");
    assert_eq!(
        compare_with(REFERENCE, REFERENCE, "-1", &[]).err(),
        Some("--tolerance -1 is negative".to_string())
    );
}

#[test]
fn verdict_is_red_on_any_count_difference() {
    let changed = with_row(
        REFERENCE,
        GENERAL_STATISTICS,
        "Overlapped transcripts",
        "842",
    );
    let report = compare(REFERENCE, &changed);
    assert_eq!(report.verdict(), "red");
    assert_eq!(report.reasons, ["[General statistics]: 1 difference"]);
    assert_eq!(
        member(&report, &["differences"]),
        &Json::object([(
            "General statistics",
            Json::List(vec![Json::object([
                ("kind", Json::str("value")),
                ("row", Json::str("Overlapped transcripts")),
                ("expected", Json::str("843")),
                ("actual", Json::str("842")),
            ])])
        )])
    );
    // A [Data version] row, a missing row and a missing section are each red.
    let data_version = with_row(REFERENCE, DATA_VERSION, "dbSNP", "157");
    assert_eq!(compare(REFERENCE, &data_version).verdict(), "red");
    let rows = section_lines(REFERENCE, GENERAL_STATISTICS);
    let fewer = with_section_lines(REFERENCE, GENERAL_STATISTICS, &rows[..rows.len() - 1]);
    let report = compare(REFERENCE, &fewer);
    assert_eq!(
        differences(&report, "General statistics"),
        &Json::List(vec![Json::object([
            ("kind", Json::str("missing")),
            ("row", Json::str("Overlapped regulatory features")),
            ("expected", Json::str("-")),
        ])])
    );
    let report = compare(REFERENCE, PREDICTIONS);
    assert_eq!(report.verdict(), "red");
    assert_eq!(
        report.reasons,
        [
            "section [SIFT summary] only in the actual file",
            "section [PolyPhen summary] only in the actual file",
        ]
    );
    let Json::List(actual_titles) = member(&report, &["sections", "actual"]) else {
        panic!("the titles are a list");
    };
    assert_eq!(
        actual_titles[7..9],
        [Json::str("SIFT summary"), Json::str("PolyPhen summary")]
    );
    assert_eq!(
        member(&report, &["sections", "order_agrees"]),
        &Json::Bool(true)
    );
    // A run with a FASTA writes the distribution charts between [Variants by chromosome]
    // and [Position in protein]; the same sections in another order are red.
    let moved = REFERENCE
        .replace(
            "\n[Position in protein]\n",
            "\n[Position in protein moved]\n",
        )
        .replace("\n[Variants by chromosome]\n", "\n[Position in protein]\n")
        .replace(
            "\n[Position in protein moved]\n",
            "\n[Variants by chromosome]\n",
        );
    let report = compare(REFERENCE, &moved);
    assert_eq!(
        member(&report, &["sections", "order_agrees"]),
        &Json::Bool(false)
    );
    assert_eq!(
        report.reasons[0],
        "the shared sections come in another order in the actual file"
    );
    assert_eq!(
        member(&report, &["sections", "expected"]),
        &Json::strs(FIXTURE_SECTIONS)
    );
    // A repeated section title is reported once; on either side its first copy is the one
    // compared and the rows compared count it once.
    let twice = format!("{REFERENCE}\n[Variant classes]\nSNV\t243\n");
    let report = compare(REFERENCE, &twice);
    assert_eq!(
        report.reasons,
        ["section [Variant classes] appears twice in the actual file"]
    );
    assert_eq!(
        member(&report, &["sections", "order_agrees"]),
        &Json::Bool(true)
    );
    let report = compare(&twice, REFERENCE);
    assert_eq!(
        report.reasons,
        ["section [Variant classes] appears twice in the expected file"]
    );
    assert_eq!(
        member(&report, &["compared"]),
        member(&compare(REFERENCE, REFERENCE), &["compared"])
    );
    // Two empty reports agree on nothing.
    assert_eq!(
        compare("", "").reasons,
        [
            "the expected file carries no [section] header",
            "the actual file carries no [section] header",
        ]
    );
}

#[test]
fn report_carries_its_keys_in_one_shape() {
    let report = compare(REFERENCE, REFERENCE);
    let Json::Object(members) = &report.comparison else {
        panic!("the report is an object");
    };
    let keys: Vec<&str> = members.iter().map(|(key, _)| key.as_str()).collect();
    assert_eq!(
        keys,
        [
            "mode",
            "format",
            "verdict",
            "reasons",
            "tolerance",
            "masked_rows",
            "sections",
            "compared",
            "documented",
            "differences",
        ]
    );
    assert_eq!(member(&report, &["documented"]), &Json::List(Vec::new()));
    assert_eq!(member(&report, &["mode"]), &Json::str("stats"));
    assert_eq!(member(&report, &["format"]), &Json::str("stats"));
    assert_eq!(member(&report, &["verdict"]), &Json::str("green"));
    assert_eq!(member(&report, &["reasons"]), &Json::List(Vec::new()));
    // The command line adds the two paths, after every key of the comparison.
    let Json::Object(with_paths) = report.json(b"e/ref.txt", b"a/out_summary.txt") else {
        panic!("the report is an object");
    };
    assert_eq!(with_paths[..10], members[..]);
    assert_eq!(
        with_paths[10..],
        [
            ("expected".to_string(), Json::Arg(b"e/ref.txt".to_vec())),
            (
                "actual".to_string(),
                Json::Arg(b"a/out_summary.txt".to_vec())
            ),
        ]
    );
    assert_eq!(report.summary(), "stats: green\n");
    let red = compare(
        REFERENCE,
        &with_row(REFERENCE, GENERAL_STATISTICS, "Variants processed", "338"),
    );
    assert_eq!(
        red.summary(),
        "stats: red\n  [General statistics]: 1 difference\n"
    );
}

/// A made-up report of a run under `--check_frequency` in the shape Ensembl VEP writes
/// it: eleven variants processed and every one of them counted as existing.
const FREQUENCY: &str = include_str!("../tests/fixtures/frequency_expected.summary.log");
/// The candidate's report of the same run, with ten of the eleven counted as existing.
const FREQUENCY_CANDIDATE: &str = include_str!("../tests/fixtures/frequency_actual.summary.log");
const PROCESSED: &str = "Variants processed";
const EXISTING: &str = "Novel / existing variants";
const ALL_OF_ELEVEN: &str = "0 (0.0) / 11 (100.0)";
const TEN_OF_ELEVEN: &str = "1 (9.1) / 10 (90.9)";

/// The two frequency reports with `processed` variants, every one of them existing in
/// the expected report and the candidate's row reading `actual`.
fn frequency_reports(processed: &str, actual: &str) -> (String, String) {
    let all = format!("0 (0.0) / {processed} (100.0)");
    let expected = with_row(FREQUENCY, GENERAL_STATISTICS, PROCESSED, processed);
    let expected = with_row(&expected, GENERAL_STATISTICS, EXISTING, &all);
    let candidate = with_row(
        FREQUENCY_CANDIDATE,
        GENERAL_STATISTICS,
        PROCESSED,
        processed,
    );
    let candidate = with_row(&candidate, GENERAL_STATISTICS, EXISTING, actual);
    (expected, candidate)
}

fn frequency(processed: &str, actual: &str) -> Report {
    let (expected, candidate) = frequency_reports(processed, actual);
    compare(&expected, &candidate)
}

/// The `documented` list of a report whose one row set aside is the existing-variant row.
fn set_aside(expected: &str, actual: &str) -> Json {
    Json::List(vec![Json::object([
        (
            "class",
            Json::str("stats_existing_variants_check_frequency"),
        ),
        ("section", Json::str("General statistics")),
        ("row", Json::str(EXISTING)),
        ("expected", Json::str(expected)),
        ("actual", Json::str(actual)),
    ])])
}

/// The differences of `[General statistics]` when the existing-variant row alone differs.
fn existing_row_differs(expected: &str, actual: Json) -> Json {
    Json::List(vec![Json::object([
        ("kind", Json::str("value")),
        ("row", Json::str(EXISTING)),
        ("expected", Json::str(expected)),
        ("actual", actual),
    ])])
}

fn assert_existing_row_is_red(report: &Report, expected: &str, actual: Json) {
    assert_eq!(report.verdict(), "red", "{actual:?}");
    assert_eq!(
        differences(report, "General statistics"),
        &existing_row_differs(expected, actual)
    );
    assert_eq!(member(report, &["documented"]), &Json::List(Vec::new()));
}

#[test]
fn a_frequency_check_sets_the_all_existing_row_aside() {
    let report = compare(FREQUENCY, FREQUENCY_CANDIDATE);
    assert_eq!(report.verdict(), "green", "{:?}", report.reasons);
    assert_eq!(member(&report, &["differences"]), &Json::Object(Vec::new()));
    assert_eq!(
        member(&report, &["documented"]),
        &set_aside(ALL_OF_ELEVEN, TEN_OF_ELEVEN)
    );
    assert_eq!(
        report.summary(),
        "stats: green\n  [General statistics]: Novel / existing variants set aside as \
         stats_existing_variants_check_frequency\n"
    );
    // `Stats.pm:727-734` writes `-` for a run that counts no variant as existing.
    let dash = frequency("11", "-");
    assert_eq!(dash.verdict(), "green", "{:?}", dash.reasons);
    assert_eq!(
        member(&dash, &["documented"]),
        &set_aside(ALL_OF_ELEVEN, "-")
    );
    // `--filter_common` switches the frequency check on (`Config.pm:441-446`).
    let common = FREQUENCY.replace("--check_frequency", "--filter_common");
    assert_ne!(common, FREQUENCY);
    let report = compare(&common, FREQUENCY_CANDIDATE);
    assert_eq!(report.verdict(), "green", "{:?}", report.reasons);
    assert_eq!(
        member(&report, &["documented"]),
        &set_aside(ALL_OF_ELEVEN, TEN_OF_ELEVEN)
    );
    // The percentages are those `sprintf('%.1f')` prints of `100 * (count / processed)`:
    // a tie goes to the even digit, as 0.25, 6.25 and 93.75 do.
    for (processed, row) in [
        ("400", "399 (99.8) / 1 (0.2)"),
        ("16", "15 (93.8) / 1 (6.2)"),
        ("2000", "1999 (100.0) / 1 (0.1)"),
        ("1000", "995 (99.5) / 5 (0.5)"),
        ("40", "37 (92.5) / 3 (7.5)"),
        ("3", "2 (66.7) / 1 (33.3)"),
        ("2", "1 (50.0) / 1 (50.0)"),
    ] {
        let report = frequency(processed, row);
        assert_eq!(report.verdict(), "green", "{processed} {row}");
        let all = format!("0 (0.0) / {processed} (100.0)");
        assert_eq!(member(&report, &["documented"]), &set_aside(&all, row));
    }
    // Two rows that agree as written are no difference, and nothing is set aside.
    let same = frequency("11", ALL_OF_ELEVEN);
    assert_eq!(same.verdict(), "green", "{:?}", same.reasons);
    assert_eq!(member(&same, &["documented"]), &Json::List(Vec::new()));
    assert_eq!(same.summary(), "stats: green\n");
    // A row set aside stands beside the differences of a red report, after the reasons.
    let (expected, candidate) = frequency_reports("11", TEN_OF_ELEVEN);
    let genes = with_row(&candidate, GENERAL_STATISTICS, "Overlapped genes", "7");
    let report = compare(&expected, &genes);
    assert_eq!(
        member(&report, &["documented"]),
        &set_aside(ALL_OF_ELEVEN, TEN_OF_ELEVEN)
    );
    assert_eq!(
        report.summary(),
        "stats: red\n  [General statistics]: 1 difference\n  [General statistics]: Novel / \
         existing variants set aside as stats_existing_variants_check_frequency\n"
    );
}

#[test]
fn the_row_set_aside_is_held_to_the_processed_count() {
    for actual in [
        // The two counts sum to another total than the eleven processed.
        "2 (18.2) / 10 (90.9)",
        "0 (0.0) / 10 (100.0)",
        "1 (9.1) / 9 (81.8)",
        // The row counts more existing variants than the eleven processed.
        "0 (0.0) / 12 (100.0)",
        // A percentage is not the one of its count in eleven.
        "1 (9.0) / 10 (90.9)",
        "1 (9.1) / 10 (91.0)",
        "1 (9.09) / 10 (90.91)",
        "1 (9) / 10 (91)",
        // A run that counts no existing variant writes the dash, never a split.
        "11 (100.0) / 0 (0.0)",
        // Another spelling of the row is not the row.
        "1 (9.1)/10 (90.9)",
        "01 (9.1) / 10 (90.9)",
        "1 (9.1) / 010 (90.9)",
        "1 (9.1) / +10 (90.9)",
        "1 (9.1) / 10 (90.9) ",
        " 1 (9.1) / 10 (90.9)",
        "1 (9.1) / 10 (90.9) / 0 (0.0)",
        "1 / 10",
        "10",
        "--",
        "",
    ] {
        let report = frequency("11", actual);
        assert_existing_row_is_red(&report, ALL_OF_ELEVEN, Json::str(actual));
    }
    assert_eq!(frequency("16", "15 (93.8) / 1 (6.3)").verdict(), "red");
    assert_eq!(frequency("400", "399 (99.8) / 1 (0.3)").verdict(), "red");
    // A candidate row without a value is no count of existing variants.
    let valueless = FREQUENCY_CANDIDATE.replace(
        "Novel / existing variants\t1 (9.1) / 10 (90.9)\n",
        "Novel / existing variants\n",
    );
    assert_ne!(valueless, FREQUENCY_CANDIDATE);
    let report = compare(FREQUENCY, &valueless);
    assert_existing_row_is_red(&report, ALL_OF_ELEVEN, Json::Null);
    // `--mask-row` keeps its meaning: the label's values are set aside unread, and the
    // report names no documented row.
    let (expected, candidate) = frequency_reports("11", "2 (18.2) / 10 (90.9)");
    let masked = compare_with(&expected, &candidate, "0", &[EXISTING]).unwrap();
    assert_eq!(masked.verdict(), "green", "{:?}", masked.reasons);
    assert_eq!(member(&masked, &["documented"]), &Json::List(Vec::new()));
}

#[test]
fn the_rule_reads_the_frequency_check_from_the_expected_command_line() {
    let red = |expected: &str| {
        let report = compare(expected, FREQUENCY_CANDIDATE);
        assert_existing_row_is_red(&report, ALL_OF_ELEVEN, Json::str(TEN_OF_ELEVEN));
    };
    // The candidate's command line names the flag in every case below: the rule reads the
    // expected report's alone.
    assert!(FREQUENCY_CANDIDATE.contains(" --check_frequency "));
    let without = include_str!("../tests/fixtures/frequency_off_expected.summary.log");
    assert!(!without.contains("--check_frequency") && !without.contains("--filter_common"));
    red(without);
    // `Config.pm:745` writes a value other than 1 after its option, and 0 leaves the
    // check off.
    red(&FREQUENCY.replace("--check_frequency", "--check_frequency 0"));
    red(&FREQUENCY.replace("--check_frequency", "--filter_common 0"));
    // An option is its whole word: a longer name and a value that carries the name are
    // not the flag.
    red(&FREQUENCY.replace("--check_frequency", "--check_frequency_of"));
    red(&FREQUENCY.replace("--check_frequency", "--plugin Name,--check_frequency"));
    red(&FREQUENCY.replace("--check_frequency", "check_frequency"));
    // A report with no command line row, a row without a value or no run section
    // states no frequency check.
    let run_rows = section_lines(FREQUENCY, RUN_STATISTICS);
    let no_options: Vec<&str> = run_rows
        .iter()
        .copied()
        .filter(|row| !row.starts_with("Command line options\t"))
        .collect();
    assert_eq!(no_options.len(), run_rows.len() - 1);
    red(&with_section_lines(FREQUENCY, RUN_STATISTICS, &no_options));
    let valueless: Vec<&str> = run_rows
        .iter()
        .map(|row| match row.strip_prefix("Command line options\t") {
            Some(_) => "Command line options",
            None => row,
        })
        .collect();
    red(&with_section_lines(FREQUENCY, RUN_STATISTICS, &valueless));
    let general_on = FREQUENCY.find("[General statistics]").unwrap();
    red(&FREQUENCY[general_on..]);
    // One of two flags switched off leaves the other on.
    let both = FREQUENCY.replace("--check_frequency", "--check_frequency 0 --filter_common");
    let report = compare(&both, FREQUENCY_CANDIDATE);
    assert_eq!(report.verdict(), "green", "{:?}", report.reasons);
}

#[test]
fn an_expected_row_of_another_shape_compares_as_written() {
    // Under the flag, an expected row that counts a variant as novel is a count like any
    // other.
    let split = with_row(
        FREQUENCY,
        GENERAL_STATISTICS,
        EXISTING,
        "2 (18.2) / 9 (81.8)",
    );
    let report = compare(&split, FREQUENCY_CANDIDATE);
    assert_existing_row_is_red(&report, "2 (18.2) / 9 (81.8)", Json::str(TEN_OF_ELEVEN));
    let dash = with_row(FREQUENCY, GENERAL_STATISTICS, EXISTING, "-");
    let report = compare(&dash, FREQUENCY_CANDIDATE);
    assert_existing_row_is_red(&report, "-", Json::str(TEN_OF_ELEVEN));
    // The expected row counts every variant of another total than the report's
    // `Variants processed`.
    let twelve = with_row(
        FREQUENCY,
        GENERAL_STATISTICS,
        EXISTING,
        "0 (0.0) / 12 (100.0)",
    );
    for actual in [TEN_OF_ELEVEN, "2 (16.7) / 10 (83.3)", "-"] {
        let candidate = with_row(FREQUENCY_CANDIDATE, GENERAL_STATISTICS, EXISTING, actual);
        let report = compare(&twelve, &candidate);
        assert_existing_row_is_red(&report, "0 (0.0) / 12 (100.0)", Json::str(actual));
    }
    // `Variants processed` is no plain count, or the expected report lacks the row.
    for processed in ["11 ", "+11", "011", "1e1", "11.0", "", "0"] {
        let all = format!("0 (0.0) / {processed} (100.0)");
        let report = frequency(processed, "-");
        assert_existing_row_is_red(&report, &all, Json::str("-"));
    }
    let rows = section_lines(FREQUENCY, GENERAL_STATISTICS);
    let unprocessed: Vec<&str> = rows
        .iter()
        .copied()
        .filter(|row| !row.starts_with("Variants processed\t"))
        .collect();
    assert_eq!(unprocessed.len(), rows.len() - 1);
    let expected = with_section_lines(FREQUENCY, GENERAL_STATISTICS, &unprocessed);
    let report = compare(&expected, FREQUENCY_CANDIDATE);
    assert_eq!(report.verdict(), "red");
    assert_eq!(member(&report, &["documented"]), &Json::List(Vec::new()));
    // The row of another section is a row like any other: a chart compares it as a count
    // and another tabular section as written.
    let chart = "[VEP run statistics]\nCommand line options\tvep --check_frequency\n\n\
                 [Variant classes]\nNovel / existing variants\t0 (0.0) / 11 (100.0)\n";
    let report = compare(chart, &chart.replace(ALL_OF_ELEVEN, TEN_OF_ELEVEN));
    assert_eq!(report.verdict(), "red");
    assert_eq!(member(&report, &["documented"]), &Json::List(Vec::new()));
    let tabular = "[VEP run statistics]\nCommand line options\tvep --check_frequency\n\n\
                   [Data version]\nVariants processed\t11\n\
                   Novel / existing variants\t0 (0.0) / 11 (100.0)\n";
    let report = compare(tabular, &tabular.replace(ALL_OF_ELEVEN, TEN_OF_ELEVEN));
    assert_eq!(
        differences(&report, "Data version"),
        &existing_row_differs(ALL_OF_ELEVEN, Json::str(TEN_OF_ELEVEN))
    );
    assert_eq!(member(&report, &["documented"]), &Json::List(Vec::new()));
}

/// Ensembl VEP 116.3's report of one record of the `GRCh38-custom` corpus under an `exact`
/// VCF source whose key is `Consequence`, its paths and times rewritten and its command
/// line without `--dir_plugins`, `--fork` and `--format`. The source's record `exact_3`
/// names one allele of the record, and the report counts its name once for each of that
/// allele's nine rows.
const CONSEQUENCE_KEY: &str =
    include_str!("../tests/fixtures/consequence_key_expected.summary.log");
/// vep-rs's report of the same run, its paths and times rewritten: the rows' terms alone.
const CONSEQUENCE_KEY_CANDIDATE: &str =
    include_str!("../tests/fixtures/consequence_key_actual.summary.log");
const ALL_CONSEQUENCES: &str = "Consequences (all)";
const KEYED_SOURCE: &str =
    "--custom file=[PATH]/sites.vcf.gz,short_name=Consequence,format=vcf,type=exact,fields=AF%CLASS";

/// One row of `documented` for a record name the expected report alone counts.
fn name_set_aside(name: &str, count: &str) -> Json {
    Json::object([
        (
            "class",
            Json::str("stats_consequences_custom_key_named_consequence"),
        ),
        ("section", Json::str(ALL_CONSEQUENCES)),
        ("row", Json::str(name)),
        ("expected", Json::str(count)),
        ("actual", Json::Null),
    ])
}

/// One difference of a chart section.
fn chart_difference(kind: &str, row: &str, expected: Json, actual: Json) -> Json {
    Json::object([
        ("kind", Json::str(kind)),
        ("row", Json::str(row)),
        ("expected", expected),
        ("actual", actual),
    ])
}

/// The expected report with its `--custom` source written as `source`.
fn keyed_as(source: &str) -> String {
    let text = CONSEQUENCE_KEY.replace(KEYED_SOURCE, source);
    assert_ne!(text, CONSEQUENCE_KEY, "{source}");
    text
}

#[test]
fn a_source_keyed_consequence_sets_the_record_name_aside() {
    assert_eq!(
        section_lines(CONSEQUENCE_KEY, ALL_CONSEQUENCES),
        [
            "exact_3\t9",
            "missense_variant\t7",
            "intron_variant\t2",
            "NMD_transcript_variant\t2"
        ]
    );
    let report = compare(CONSEQUENCE_KEY, CONSEQUENCE_KEY_CANDIDATE);
    assert_eq!(report.verdict(), "green", "{:?}", report.reasons);
    assert_eq!(member(&report, &["differences"]), &Json::Object(Vec::new()));
    assert_eq!(
        member(&report, &["documented"]),
        &Json::List(vec![name_set_aside("exact_3", "9")])
    );
    assert_eq!(
        report.summary(),
        "stats: green\n  [Consequences (all)]: exact_3 set aside as \
         stats_consequences_custom_key_named_consequence\n"
    );
    // The source by position, the key second (`AnnotationSourceAdaptor.pm:356`), and the
    // source as the second of two.
    for source in [
        "--custom [PATH]/sites.vcf.gz,Consequence,vcf,exact,0,AF,CLASS",
        "--custom file=[PATH]/regions.bed.gz,short_name=regions,format=bed --custom \
         file=[PATH]/sites.vcf.gz,short_name=Consequence,format=vcf,type=exact",
        "--custom file=[PATH]/sites.vcf.gz,short_name=sites,short_name=Consequence,format=vcf",
    ] {
        let report = compare(&keyed_as(source), CONSEQUENCE_KEY_CANDIDATE);
        assert_eq!(report.verdict(), "green", "{source}: {:?}", report.reasons);
        assert_eq!(
            member(&report, &["documented"]),
            &Json::List(vec![name_set_aside("exact_3", "9")]),
            "{source}"
        );
    }
    // Every name is set aside, in the expected report's order.
    let mut lines = section_lines(CONSEQUENCE_KEY, ALL_CONSEQUENCES);
    lines.insert(1, "21:10542415-10542415\t4");
    let two_names = with_section_lines(CONSEQUENCE_KEY, ALL_CONSEQUENCES, &lines);
    let report = compare(&two_names, CONSEQUENCE_KEY_CANDIDATE);
    assert_eq!(report.verdict(), "green", "{:?}", report.reasons);
    assert_eq!(
        member(&report, &["documented"]),
        &Json::List(vec![
            name_set_aside("exact_3", "9"),
            name_set_aside("21:10542415-10542415", "4"),
        ])
    );
    // A name set aside stands beside the differences of a red report, after the reasons.
    let fewer = with_row(
        CONSEQUENCE_KEY_CANDIDATE,
        ALL_CONSEQUENCES,
        "missense_variant",
        "6",
    );
    let report = compare(CONSEQUENCE_KEY, &fewer);
    assert_eq!(
        differences(&report, ALL_CONSEQUENCES),
        &Json::List(vec![chart_difference(
            "count",
            "missense_variant",
            int(7),
            int(6)
        )])
    );
    assert_eq!(
        report.summary(),
        "stats: red\n  [Consequences (all)]: 1 difference\n  [Consequences (all)]: exact_3 set \
         aside as stats_consequences_custom_key_named_consequence\n"
    );
}

#[test]
fn the_rule_reads_the_key_from_the_expected_command_line() {
    let name_differs = Json::List(vec![chart_difference(
        "missing",
        "exact_3",
        int(9),
        Json::Null,
    )]);
    let unkeyed = [
        // Another key, a key that differs in a letter or in its case, and a later item
        // of the key that replaces `Consequence`.
        keyed_as("--custom file=[PATH]/sites.vcf.gz,short_name=sites,format=vcf,type=exact"),
        keyed_as("--custom file=[PATH]/sites.vcf.gz,short_name=Consequences,format=vcf"),
        keyed_as("--custom file=[PATH]/sites.vcf.gz,short_name=consequence,format=vcf"),
        keyed_as("--custom file=[PATH]/sites.vcf.gz,short_name=Consequence,short_name=sites"),
        // By position the key is the second item alone.
        keyed_as("--custom [PATH]/Consequence,sites,vcf,exact,0,Consequence"),
        // Two sources, the later one by position, as the path mask writes them
        // (`Config.pm:796`): without the later one's `--custom` and the earlier one's last
        // value, which is its format in the first and its key in the second.
        keyed_as(
            "--custom [PATH]/regions.bed.gz,regions,[PATH]/sites.vcf.gz,Consequence,vcf,exact,0",
        ),
        keyed_as(
            "--custom file=[PATH]/sites.vcf.gz,format=vcf,type=exact,\
             short_name=[PATH]/regions.bed.gz,regions,bed",
        ),
        // The word with no value after it, and the key under another option.
        keyed_as("--custom"),
        keyed_as("--plugin file=[PATH]/sites.vcf.gz,short_name=Consequence"),
        // No command line at all.
        CONSEQUENCE_KEY
            .lines()
            .filter(|line| !line.starts_with("Command line options\t"))
            .collect::<Vec<_>>()
            .join("\n"),
    ];
    for expected in &unkeyed {
        let report = compare(expected, CONSEQUENCE_KEY_CANDIDATE);
        assert_eq!(report.verdict(), "red");
        assert_eq!(differences(&report, ALL_CONSEQUENCES), &name_differs);
        assert_eq!(member(&report, &["documented"]), &Json::List(Vec::new()));
    }
    // The candidate's command line is not read: its own report names the key here.
    let report = compare(&unkeyed[0], CONSEQUENCE_KEY_CANDIDATE);
    assert!(CONSEQUENCE_KEY_CANDIDATE.contains(KEYED_SOURCE));
    assert_eq!(differences(&report, ALL_CONSEQUENCES), &name_differs);
}

#[test]
fn a_consequence_term_or_a_key_the_candidate_holds_is_compared() {
    // A term the candidate lacks is a difference under the key too, in each vocabulary.
    for (term, count) in [("intron_variant", 2), ("NMD_transcript_variant", 2)] {
        let lines: Vec<&str> = section_lines(CONSEQUENCE_KEY_CANDIDATE, ALL_CONSEQUENCES)
            .into_iter()
            .filter(|line| line.split('\t').next() != Some(term))
            .collect();
        let candidate = with_section_lines(CONSEQUENCE_KEY_CANDIDATE, ALL_CONSEQUENCES, &lines);
        let report = compare(CONSEQUENCE_KEY, &candidate);
        assert_eq!(report.verdict(), "red", "{term}");
        assert_eq!(
            differences(&report, ALL_CONSEQUENCES),
            &Json::List(vec![chart_difference(
                "missing",
                term,
                int(count),
                Json::Null
            )])
        );
        assert_eq!(
            member(&report, &["documented"]),
            &Json::List(vec![name_set_aside("exact_3", "9")])
        );
    }
    let relabel = |text: &str| {
        text.replace("intron_variant\t2", "INTRONIC\t2")
            .replace("missense_variant\t7", "missense\t7")
    };
    let display = relabel(CONSEQUENCE_KEY);
    let candidate = relabel(CONSEQUENCE_KEY_CANDIDATE).replace("INTRONIC\t2\n", "");
    let report = compare(&display, &candidate);
    assert_eq!(
        differences(&report, ALL_CONSEQUENCES),
        &Json::List(vec![chart_difference(
            "missing",
            "INTRONIC",
            int(2),
            Json::Null
        )])
    );
    // The first `missense` line of a report is that of `[Consequences (all)]`.
    let candidate = relabel(CONSEQUENCE_KEY_CANDIDATE).replacen("missense\t7\n", "", 1);
    let report = compare(&display, &candidate);
    assert_eq!(
        differences(&report, ALL_CONSEQUENCES),
        &Json::List(vec![chart_difference(
            "missing",
            "missense",
            int(7),
            Json::Null
        )])
    );
    // A name the candidate holds compares as any key: equal counts agree and nothing is
    // set aside, other counts differ.
    let mut lines = section_lines(CONSEQUENCE_KEY_CANDIDATE, ALL_CONSEQUENCES);
    lines.insert(0, "exact_3\t9");
    let named = with_section_lines(CONSEQUENCE_KEY_CANDIDATE, ALL_CONSEQUENCES, &lines);
    let report = compare(CONSEQUENCE_KEY, &named);
    assert_eq!(report.verdict(), "green", "{:?}", report.reasons);
    assert_eq!(member(&report, &["documented"]), &Json::List(Vec::new()));
    let report = compare(
        CONSEQUENCE_KEY,
        &with_row(&named, ALL_CONSEQUENCES, "exact_3", "4"),
    );
    assert_eq!(
        differences(&report, ALL_CONSEQUENCES),
        &Json::List(vec![chart_difference("count", "exact_3", int(9), int(4))])
    );
    assert_eq!(member(&report, &["documented"]), &Json::List(Vec::new()));
    // A name the candidate alone holds, a count that is no integer and a name in another
    // chart are differences.
    let report = compare(CONSEQUENCE_KEY_CANDIDATE, &named);
    assert_eq!(
        differences(&report, ALL_CONSEQUENCES),
        &Json::List(vec![chart_difference(
            "unexpected",
            "exact_3",
            Json::Null,
            int(9)
        )])
    );
    let report = compare(
        &with_row(CONSEQUENCE_KEY, ALL_CONSEQUENCES, "exact_3", "nine"),
        CONSEQUENCE_KEY_CANDIDATE,
    );
    assert_eq!(
        differences(&report, ALL_CONSEQUENCES),
        &Json::List(vec![chart_difference(
            "missing",
            "exact_3",
            Json::str("nine"),
            Json::Null
        )])
    );
    let coding = with_section_lines(
        CONSEQUENCE_KEY,
        "Coding consequences",
        &["exact_3\t9", "missense_variant\t7"],
    );
    let report = compare(&coding, CONSEQUENCE_KEY_CANDIDATE);
    assert_eq!(
        differences(&report, "Coding consequences"),
        &Json::List(vec![chart_difference(
            "missing",
            "exact_3",
            int(9),
            Json::Null
        )])
    );
    assert_eq!(
        member(&report, &["documented"]),
        &Json::List(vec![name_set_aside("exact_3", "9")])
    );
}
