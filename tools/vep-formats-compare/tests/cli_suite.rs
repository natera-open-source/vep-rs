// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Command-line cases, each run through the binary on planted files and held to the
//! report, the standard streams and the exit code it must give.

mod support;

use indexmap::IndexMap;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use support::{args_of, compact, golden, member, report, run, scratch, text_of, Ran};
use vep_formats_compare::compare::{compare_entries, CompareOptions, FLOOR_NOT_APPLICABLE};
use vep_formats_compare::compare_files::compare_files;
use vep_formats_compare::consts::EXAMPLES_PER_GROUP;
use vep_formats_compare::floors::json_floor_fields;
use vep_formats_compare::input::input_records;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::nearest::NearestSites;
use vep_formats_compare::parse::recoder::recoder_identifier_indel;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::pyjson::{self, DumpOptions, Value};
use vep_formats_compare::rules::Side;

const ERROR: &str = "ERROR: [compare_vep_formats] ";

fn json(text: &str) -> Value {
    pyjson::loads(text).unwrap_or_else(|error| panic!("{text} is not JSON: {error}"))
}

/// Whether `value` is the value `expected` spells, the order of an object's members
/// aside.
fn is(value: &Value, expected: &str) -> bool {
    pyjson::py_eq(value, &json(expected))
}

fn strings(value: &Value) -> Vec<&str> {
    value
        .as_list()
        .expect("a list")
        .iter()
        .map(|item| item.as_str().expect("a text"))
        .collect()
}

fn path_text(path: &Path) -> &str {
    path.to_str().expect("the path is UTF-8")
}

fn corpus() -> PathBuf {
    golden().join("GRCh38-release116")
}

fn expected(name: &str) -> PathBuf {
    corpus().join("expected").join(name)
}

/// The data rows of an output with `rows` in the place of its own.
fn rebuilt(text: &str, rows: &[&str]) -> String {
    let (header, _) = support::header_and_data_lines(text);
    support::rebuilt(&header, rows)
}

// The fabricated rows of the manifest tests: record `i` is one deletion on one transcript,
// written in every text format with the consequence set its side gives it.
const KEPT_SET: &str = "3_prime_UTR_variant,stop_retained_variant";
const LOST_SET: &str = "3_prime_UTR_variant,stop_lost";
const OTHER_SET: &str = "3_prime_UTR_variant,coding_sequence_variant";
const HEADER_14: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";

fn impact(set: &str) -> &'static str {
    match set {
        KEPT_SET => "LOW",
        LOST_SET => "HIGH",
        _ => "MODIFIER",
    }
}

/// The VCF line, the Location and the transcript of the fabricated record `ordinal`.
fn fabricated_record(ordinal: usize) -> (String, String, String) {
    let position = 100 * (ordinal + 1);
    (
        format!("21\t{}\tdel{ordinal}\tCAT\tC\t.\t.\t.", position - 1),
        format!("21:{position}-{}", position + 1),
        format!("ENST{}", ordinal + 1),
    )
}

fn quoted(text: &str) -> String {
    pyjson::dumps_str(text, true)
}

/// The `divergences` item of the fabricated record `ordinal` under the class `class`, with
/// `more` members after its own.
fn fabricated_item(ordinal: usize, class: &str, more: &str) -> String {
    let (_, location, feature) = fabricated_record(ordinal);
    format!(
        "{{\"location\": {}, \"allele\": \"-\", \"feature\": {}, \"feature_type\": \"Transcript\", \
         \"record_indices\": [{ordinal}], \"vep_consequence_set\": {}, \
         \"vep_rs_consequence_sets\": [{}], \"expected_divergence\": {}{more}}}",
        quoted(&location),
        quoted(&feature),
        quoted(KEPT_SET),
        quoted(LOST_SET),
        quoted(class)
    )
}

fn default_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    extra: &str,
) -> String {
    format!("{name}\t{location}\t{allele}\tENSG1\t{feature}\tTranscript\t{consequence}\t-\t-\t-\t-\t-\t-\t{extra}")
}

fn default_text(rows: &[String]) -> String {
    let mut text = format!("{HEADER_14}\n");
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    text
}

/// One row per fabricated record in the format `fmt`, record `i` carrying `sets[i]`.
fn fabricated_output(fmt: &str, sets: &[&str]) -> String {
    let records: Vec<_> = (0..sets.len()).map(fabricated_record).collect();
    let rows = records.iter().zip(sets);
    match fmt {
        "default" => default_text(
            &rows
                .enumerate()
                .map(|(i, ((_, location, feature), terms))| {
                    let extra = format!("IMPACT={}", impact(terms));
                    default_row(&format!("del{i}"), location, "-", feature, terms, &extra)
                })
                .collect::<Vec<_>>(),
        ),
        "tab" => {
            let mut text = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIMPACT\n".to_string();
            for (i, ((_, location, feature), terms)) in rows.enumerate() {
                text.push_str(&format!(
                    "del{i}\t{location}\t-\tENSG1\t{feature}\tTranscript\t{terms}\t{}\n",
                    impact(terms)
                ));
            }
            text
        }
        "vcf" => {
            let mut text = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|Feature_type|Feature\">\n".to_string();
            for ((line, _, feature), terms) in rows {
                text.push_str(&format!(
                    "{}CSQ=-|{}|{}|Transcript|{feature}\n",
                    &line[..line.len() - 1],
                    terms.replace(',', "&"),
                    impact(terms)
                ));
            }
            text
        }
        "json" => {
            let mut text = String::new();
            for (i, ((line, _, feature), terms)) in rows.enumerate() {
                let term_list: Vec<String> = terms.split(',').map(quoted).collect();
                text.push_str(&format!(
                    "{{\"input\": {}, \"id\": \"del{i}\", \"allele_string\": \"AT/-\", \
                     \"most_severe_consequence\": {}, \"transcript_consequences\": \
                     [{{\"variant_allele\": \"-\", \"transcript_id\": {}, \"consequence_terms\": [{}], \
                     \"impact\": {}}}]}}\n",
                    quoted(line),
                    term_list[0],
                    quoted(feature),
                    term_list.join(", "),
                    quoted(impact(terms))
                ));
            }
            text
        }
        other => panic!("no fabricated output in the format {other}"),
    }
}

/// The exit code, the report where one was written and the run of the comparator on the
/// fabricated pair in the format `fmt`: the reference carrying the sets `expected` and the
/// candidate `actual`, read beside `manifest` under `options`.
fn compare_fabricated(
    dir: &Path,
    fmt: &str,
    expected: &[&str],
    actual: &[&str],
    manifest: Option<&str>,
    options: &[&str],
) -> (i32, Option<Value>, Ran) {
    std::fs::create_dir_all(dir).expect("the directory is made");
    let plant = |name: &str, text: String| {
        let path = dir.join(name);
        std::fs::write(&path, text).expect("a file is planted");
        path
    };
    let reference = plant("vep.out", fabricated_output(fmt, expected));
    let candidate = plant("candidate.out", fabricated_output(fmt, actual));
    let lines: String = (0..expected.len())
        .map(|i| format!("{}\n", fabricated_record(i).0))
        .collect();
    let input = plant("variants.vcf", lines);
    let written = dir.join("report.json");
    let mut argv = args_of(&[
        "--format",
        fmt,
        "--expected",
        path_text(&reference),
        "--actual",
        path_text(&candidate),
        "--report",
        path_text(&written),
        "--input",
        path_text(&input),
    ]);
    if let Some(manifest) = manifest {
        let manifest = plant("manifest.json", manifest.to_string());
        argv.extend(args_of(&["--manifest", path_text(&manifest)]));
    }
    argv.extend(args_of(options));
    let ran = run(dir, &argv);
    let report = written.is_file().then(|| report(&written));
    (ran.code, report, ran)
}

/// The red reason of one entry whose Consequence differs, as the format names the field.
fn consequence_reason(fmt: &str) -> String {
    let field = if fmt == "json" {
        "consequence_terms"
    } else {
        "Consequence"
    };
    format!("field {field}: 1 differing entries")
}

/// A manifest whose own run and one flag set each document one fabricated record over the
/// one input, with a second flag set that carries no list of its own and a third over
/// another input.
fn runs_manifest() -> String {
    format!(
        "{{\"divergences\": [{}], \"flag_sets\": [\
         {{\"name\": \"other\", \"flags\": [\"--other\"], \"divergences\": [{}]}}, \
         {{\"name\": \"plain\", \"flags\": [\"--plain\"]}}, \
         {{\"name\": \"elsewhere\", \"input\": \"elsewhere.vcf\", \"divergences\": [{}]}}]}}",
        fabricated_item(0, "some_class", ""),
        fabricated_item(1, "some_class", ""),
        fabricated_item(1, "some_class", "")
    )
}

const RUNS_MASKS: [&str; 2] = ["--masks", "some_class"];

/// The line a comparison beside `runs_manifest` ends on where it names `some_class` and no
/// run: three of the manifest's runs annotate the compared input.
const UNNAMED_RUNS: &str = "--run is not given, and the manifest lists 3 runs over variants.vcf \
    (default, other, plain) that document some_class, which --masks names: name the pair's run";
const WITH_REFERENCE: &str = "vep_rs_consequence_sets_with_reference_sequence";

/// The text formats a manifest's reader serves; a Parquet candidate is not planted here.
const TEXT_FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

fn with(options: &[&str], more: &[&str]) -> Vec<String> {
    options.iter().chain(more).map(|o| o.to_string()).collect()
}

fn refs(options: &[String]) -> Vec<&str> {
    options.iter().map(String::as_str).collect()
}

#[test]
fn field_pairs_accept_a_value_pair() {
    let text = text_of(&expected("default.txt.gz"));
    let (_, rows) = support::header_and_data_lines(&text);
    let changed = rows[0].replace("IMPACT=MODIFIER", "IMPACT=LOW");
    let mut other_rows = rows.clone();
    other_rows[0] = &changed;
    let reference = parse_default(&text).expect("the reference parses");
    let candidate = parse_default(&rebuilt(&text, &other_rows)).expect("the candidate parses");
    let spec = "IMPACT=MODIFIER=>LOW;Extra=IMPACT=MODIFIER=>IMPACT=LOW".to_string();
    let pairs =
        vep_formats_compare::cli::parse_field_pairs(Some(&[spec])).expect("the pairs parse");
    let side = Side::default();
    let options = |field_pairs| CompareOptions {
        floors: None,
        field_pairs: Some(field_pairs),
        check_order: true,
        examples: EXAMPLES_PER_GROUP,
        order_only: None,
        floor_fields: None,
        input_lines: None,
        kept: None,
        side: &side,
    };
    let accepted = compare_entries(
        &reference,
        &candidate,
        &mut Documented::new(),
        &options(&pairs),
    )
    .expect("the entries compare");
    assert_eq!(member(&accepted, &["verdict"]).as_str(), Some("green"));
    assert!(is(
        member(&accepted, &["columns", "IMPACT", "pairs_accepted"]),
        "1"
    ));
    let mut high = IndexMap::new();
    high.insert(
        "IMPACT".to_string(),
        BTreeSet::from([("MODIFIER".to_string(), "HIGH".to_string())]),
    );
    let other = compare_entries(
        &reference,
        &candidate,
        &mut Documented::new(),
        &options(&high),
    )
    .expect("the entries compare");
    assert_eq!(member(&other, &["verdict"]).as_str(), Some("red"));
}

#[test]
fn unordered_compares_the_shared_entries_without_their_order() {
    let dir = scratch("suite-unordered");
    let corpus = golden().join("GRCh38-pick-summary");
    let manifest = json(&text_of(&corpus.join("manifest.json")));
    let unordered: Vec<&str> = member(&manifest, &["flag_sets"])
        .as_list()
        .expect("the flag sets")
        .iter()
        .filter(|set| set.get("unordered").is_some_and(pyjson::truthy))
        .map(|set| member(set, &["name"]).as_str().expect("a name"))
        .collect();
    assert_eq!(unordered, ["pick-allele", "per-gene", "pick-allele-gene"]);
    let reference = corpus.join("expected_per-gene/default.txt.gz");
    let text = text_of(&reference);
    let (_, rows) = support::header_and_data_lines(&text);
    // The two gene rows of the first record, in the hash order of the reference's own run.
    let first_two: Vec<&str> = rows[..2]
        .iter()
        .map(|row| row.split('\t').nth(4).expect("a Feature"))
        .collect();
    let mut sorted = first_two.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, ["ENST00000815440", "ENST00000815693"]);
    let swapped = dir.join("swapped.txt");
    let mut swapped_rows = rows.clone();
    swapped_rows.swap(0, 1);
    std::fs::write(&swapped, rebuilt(&text, &swapped_rows)).expect("a file is planted");
    let written = dir.join("r.json");
    let args = [
        "--format",
        "default",
        "--expected",
        path_text(&reference),
        "--actual",
        path_text(&swapped),
        "--report",
        path_text(&written),
    ];
    let ran = run(&dir, &args);
    assert_eq!(
        (ran.code, ran.out().as_str()),
        (1, "default: red\n  order differs at shared position 0\n")
    );
    let report_ordered = report(&written);
    assert!(is(member(&report_ordered, &["unordered"]), "false"));
    assert_eq!(
        member(&report_ordered, &["order", "expected"]).as_str(),
        Some(format!("21:8857289|C|{}|Transcript", first_two[0]).as_str())
    );
    let ran = run(&dir, &with(&args, &["--unordered"]));
    assert_eq!((ran.code, ran.out().as_str()), (0, "default: green\n"));
    let report_unordered = report(&written);
    assert!(is(member(&report_unordered, &["unordered"]), "true"));
    assert!(is(member(&report_unordered, &["order"]), "null"));
    assert!(is(
        member(&report_unordered, &["counts", "compared"]),
        &rows.len().to_string()
    ));
    let dropped = dir.join("dropped.txt");
    std::fs::write(&dropped, rebuilt(&text, &rows[1..])).expect("a file is planted");
    let ran = run(
        &dir,
        &[
            "--format",
            "default",
            "--expected",
            path_text(&reference),
            "--actual",
            path_text(&dropped),
            "--unordered",
        ],
    );
    assert_eq!(
        (ran.code, ran.out().as_str()),
        (1, "default: red\n  missing in the candidate: 1 entries\n")
    );
}

#[test]
fn an_item_documents_its_with_reference_sets_where_they_are_not_ensembl_veps() {
    let dir = scratch("suite-with-reference-sets");
    let more = format!(", {}: [{}]", quoted(WITH_REFERENCE), quoted(OTHER_SET));
    let manifest = format!(
        "{{\"divergences\": [{}]}}",
        fabricated_item(0, "some_class", &more)
    );
    let options = ["--masks", "some_class", "--reference-sequence", "present"];
    for fmt in TEXT_FORMATS {
        let (status, report, _) = compare_fabricated(
            &dir.join(format!("observed-{fmt}")),
            fmt,
            &[KEPT_SET],
            &[OTHER_SET],
            Some(&manifest),
            &options,
        );
        let report = report.expect("a report");
        assert_eq!(
            status,
            0,
            "{fmt}: {}",
            compact(member(&report, &["reasons"]))
        );
        assert!(
            is(
                member(&report, &["masks", "applied"]),
                "{\"some_class\": 1}"
            ),
            "{fmt}"
        );
        let (status, report, _) = compare_fabricated(
            &dir.join(format!("crossed-{fmt}")),
            fmt,
            &[KEPT_SET],
            &[LOST_SET],
            Some(&manifest),
            &options,
        );
        let reasons = compact(member(&report.expect("a report"), &["reasons"]));
        assert_eq!(
            (status, reasons.as_str()),
            (1, "[\"documented divergences not observed: 1\"]"),
            "{fmt}"
        );
    }
}

fn a_named_run_reads_its_own_items_alone(fmt: &str) {
    let dir = scratch(&format!("suite-named-run-alone-{fmt}"));
    let manifest = runs_manifest();
    let own_pair: (&[&str], &[&str]) = (&[KEPT_SET, KEPT_SET], &[LOST_SET, KEPT_SET]);
    for name in ["default", "plain"] {
        let options = with(&RUNS_MASKS, &["--run", name]);
        let (status, report, _) = compare_fabricated(
            &dir.join(name),
            fmt,
            own_pair.0,
            own_pair.1,
            Some(&manifest),
            &refs(&options),
        );
        let report = report.expect("a report");
        assert_eq!(
            status,
            0,
            "{name}: {}",
            compact(member(&report, &["reasons"]))
        );
        assert_eq!(member(&report, &["verdict"]).as_str(), Some("green"));
        assert!(is(
            member(&report, &["masks", "applied"]),
            "{\"some_class\": 1}"
        ));
        assert!(is(
            member(&report, &["masks", "available"]),
            "{\"some_class\": 1}"
        ));
        assert_eq!(member(&report, &["manifest_run"]).as_str(), Some(name));
    }
    // Unnamed, the pair would be held to the item of the run `other` as well.
    let (status, report, ran) = compare_fabricated(
        &dir.join("unnamed"),
        fmt,
        own_pair.0,
        own_pair.1,
        Some(&manifest),
        &RUNS_MASKS,
    );
    assert_eq!(status, 2);
    assert!(report.is_none());
    assert_eq!(ran.err(), format!("{ERROR}{UNNAMED_RUNS}\n"));
}

fn a_named_run_still_reports_its_own_item_that_is_not_observed(fmt: &str) {
    let dir = scratch(&format!("suite-named-run-not-observed-{fmt}"));
    let manifest = runs_manifest();
    for (name, record) in [("default", 0), ("other", 1)] {
        let options = with(&RUNS_MASKS, &["--run", name]);
        let (status, report, _) = compare_fabricated(
            &dir.join(name),
            fmt,
            &[KEPT_SET, KEPT_SET],
            &[KEPT_SET, KEPT_SET],
            Some(&manifest),
            &refs(&options),
        );
        let report = report.expect("a report");
        assert_eq!(status, 1, "{name}");
        assert!(
            is(
                member(&report, &["reasons"]),
                "[\"documented divergences not observed: 1\"]"
            ),
            "{name}"
        );
        let stale = strings(member(&report, &["stale_divergences"]));
        assert_eq!(stale.len(), 1);
        assert!(stale[0].contains(&format!("ENST{}", record + 1)));
    }
}

fn a_difference_only_another_run_documents_is_not_set_aside(fmt: &str) {
    let dir = scratch(&format!("suite-other-run-difference-{fmt}"));
    let manifest = runs_manifest();
    let both_differ: (&[&str], &[&str]) = (&[KEPT_SET, KEPT_SET], &[LOST_SET, LOST_SET]);
    let options = with(&RUNS_MASKS, &["--run", "default"]);
    let (status, report, _) = compare_fabricated(
        &dir.join("own"),
        fmt,
        both_differ.0,
        both_differ.1,
        Some(&manifest),
        &refs(&options),
    );
    let report = report.expect("a report");
    assert_eq!(status, 1);
    let reasons = strings(member(&report, &["reasons"]));
    assert!(
        reasons.contains(&consequence_reason(fmt).as_str()),
        "{reasons:?}"
    );
    assert!(is(
        member(&report, &["masks", "applied"]),
        "{\"some_class\": 1}"
    ));
    assert!(is(member(&report, &["stale_divergences"]), "[]"));
    // Unnamed, the item of the run `other` would set the second difference aside.
    let (status, report, ran) = compare_fabricated(
        &dir.join("unnamed"),
        fmt,
        both_differ.0,
        both_differ.1,
        Some(&manifest),
        &RUNS_MASKS,
    );
    assert_eq!(status, 2);
    assert!(report.is_none());
    assert_eq!(ran.err(), format!("{ERROR}{UNNAMED_RUNS}\n"));
}

/// A comparison that names no run reads the lists of every run over its input where they
/// document no class it names, and the one run's where a single run annotates the input.
fn an_unnamed_comparison_reads_what_no_run_would_read_otherwise(fmt: &str) {
    let dir = scratch(&format!("suite-unnamed-run-{fmt}"));
    let manifest = runs_manifest();
    let own_pair: (&[&str], &[&str]) = (&[KEPT_SET, KEPT_SET], &[LOST_SET, KEPT_SET]);
    for (name, options) in [
        ("another-class", &["--masks", "another_class"][..]),
        ("no-class", &[][..]),
    ] {
        let (status, report, _) = compare_fabricated(
            &dir.join(name),
            fmt,
            own_pair.0,
            own_pair.1,
            Some(&manifest),
            options,
        );
        let report = report.expect("a report");
        assert_eq!(status, 1, "{name}");
        let reasons = strings(member(&report, &["reasons"]));
        assert!(
            reasons.contains(&consequence_reason(fmt).as_str()),
            "{name}: {reasons:?}"
        );
        assert!(is(member(&report, &["masks", "applied"]), "{}"), "{name}");
        assert!(
            is(
                member(&report, &["masks", "available"]),
                "{\"some_class\": 2}"
            ),
            "{name}"
        );
        assert!(is(member(&report, &["manifest_run"]), "null"), "{name}");
    }
    // The second run of this manifest annotates another input, so one run is read.
    let one_run = format!(
        "{{\"divergences\": [{}], \"flag_sets\": [\
         {{\"name\": \"elsewhere\", \"input\": \"elsewhere.vcf\", \"divergences\": [{}]}}]}}",
        fabricated_item(0, "some_class", ""),
        fabricated_item(1, "some_class", "")
    );
    let (status, report, _) = compare_fabricated(
        &dir.join("one-run"),
        fmt,
        own_pair.0,
        own_pair.1,
        Some(&one_run),
        &RUNS_MASKS,
    );
    let report = report.expect("a report");
    assert_eq!(status, 0, "{}", compact(member(&report, &["reasons"])));
    assert!(is(
        member(&report, &["masks", "applied"]),
        "{\"some_class\": 1}"
    ));
    assert!(is(member(&report, &["manifest_run"]), "null"));
}

/// One test a text format for each of the three cases of a named run and for the unnamed one.
macro_rules! per_format {
    ($($name:ident: $case:ident($fmt:literal);)*) => {
        $(
            #[test]
            fn $name() {
                $case($fmt);
            }
        )*
    };
}

per_format! {
    a_named_run_reads_its_own_items_alone_in_the_default_format: a_named_run_reads_its_own_items_alone("default");
    a_named_run_reads_its_own_items_alone_in_the_tab_format: a_named_run_reads_its_own_items_alone("tab");
    a_named_run_reads_its_own_items_alone_in_the_vcf_format: a_named_run_reads_its_own_items_alone("vcf");
    a_named_run_reads_its_own_items_alone_in_the_json_format: a_named_run_reads_its_own_items_alone("json");
    a_named_run_still_reports_its_own_item_that_is_not_observed_in_the_default_format: a_named_run_still_reports_its_own_item_that_is_not_observed("default");
    a_named_run_still_reports_its_own_item_that_is_not_observed_in_the_tab_format: a_named_run_still_reports_its_own_item_that_is_not_observed("tab");
    a_named_run_still_reports_its_own_item_that_is_not_observed_in_the_vcf_format: a_named_run_still_reports_its_own_item_that_is_not_observed("vcf");
    a_named_run_still_reports_its_own_item_that_is_not_observed_in_the_json_format: a_named_run_still_reports_its_own_item_that_is_not_observed("json");
    a_difference_only_another_run_documents_is_not_set_aside_in_the_default_format: a_difference_only_another_run_documents_is_not_set_aside("default");
    a_difference_only_another_run_documents_is_not_set_aside_in_the_tab_format: a_difference_only_another_run_documents_is_not_set_aside("tab");
    a_difference_only_another_run_documents_is_not_set_aside_in_the_vcf_format: a_difference_only_another_run_documents_is_not_set_aside("vcf");
    a_difference_only_another_run_documents_is_not_set_aside_in_the_json_format: a_difference_only_another_run_documents_is_not_set_aside("json");
    an_unnamed_comparison_reads_what_no_run_would_read_otherwise_in_the_default_format: an_unnamed_comparison_reads_what_no_run_would_read_otherwise("default");
    an_unnamed_comparison_reads_what_no_run_would_read_otherwise_in_the_tab_format: an_unnamed_comparison_reads_what_no_run_would_read_otherwise("tab");
    an_unnamed_comparison_reads_what_no_run_would_read_otherwise_in_the_vcf_format: an_unnamed_comparison_reads_what_no_run_would_read_otherwise("vcf");
    an_unnamed_comparison_reads_what_no_run_would_read_otherwise_in_the_json_format: an_unnamed_comparison_reads_what_no_run_would_read_otherwise("json");
}

#[test]
fn a_named_runs_record_level_item_has_nothing_to_observe_on_a_record_object() {
    let dir = scratch("suite-record-level-item");
    let (line, location, _) = fabricated_record(0);
    let item = fabricated_item(0, "some_class", "")
        .replace("\"feature\": \"ENST1\"", "\"feature\": \"-\"")
        .replace(
            "\"feature_type\": \"Transcript\"",
            "\"feature_type\": \"-\"",
        );
    let manifest = format!(
        "{{\"divergences\": [], \"flag_sets\": [{{\"name\": \"summary\", \"flags\": [\"--summary\"], \"divergences\": [{item}]}}]}}"
    );
    std::fs::write(dir.join("variants.vcf"), format!("{line}\n")).expect("a file is planted");
    std::fs::write(dir.join("manifest.json"), manifest).expect("a file is planted");
    let compare = |fmt: &str, reference: String, candidate: String, options: &[&str]| {
        std::fs::write(dir.join("vep.out"), reference).expect("a file is planted");
        std::fs::write(dir.join("candidate.out"), candidate).expect("a file is planted");
        let written = dir.join("report.json");
        let mut argv = vec![
            "--format",
            fmt,
            "--expected",
            "vep.out",
            "--actual",
            "candidate.out",
            "--report",
            "report.json",
            "--input",
            "variants.vcf",
            "--manifest",
            "manifest.json",
            "--masks",
            "some_class",
        ];
        argv.extend(options);
        std::fs::remove_file(&written).ok();
        let ran = run(&dir, &argv);
        let report = written.is_file().then(|| report(&written));
        (ran.code, report, ran)
    };
    let record_object = |severe: &str| {
        format!(
            "{{\"input\": {}, \"id\": \"del0\", \"allele_string\": \"AT/-\", \"most_severe_consequence\": {}}}\n",
            quoted(&line),
            quoted(severe)
        )
    };
    let summary_row = |terms: &str| {
        let extra = format!("IMPACT={}", impact(terms));
        // The one row of a `--summary` run names no gene, feature or feature type.
        default_text(&[default_row("del0", &location, "-", "-", terms, &extra)
            .replace("\tENSG1\t-\tTranscript\t", "\t-\t-\t-\t")])
    };
    let (status, report, _) = compare(
        "json",
        record_object("stop_retained_variant"),
        record_object("stop_lost"),
        &["--run", "summary"],
    );
    let report = report.expect("a report");
    assert_eq!(status, 0);
    assert!(is(member(&report, &["stale_divergences"]), "[]"));
    assert!(is(
        member(&report, &["masks", "applied"]),
        "{\"some_class\": 1}"
    ));
    // The manifest's own run annotates the input too, and the item is the run `summary`'s.
    let (status, report, ran) = compare(
        "json",
        record_object("stop_retained_variant"),
        record_object("stop_lost"),
        &[],
    );
    assert_eq!(status, 2);
    assert!(report.is_none());
    assert_eq!(
        ran.err(),
        format!(
            "{ERROR}--run is not given, and the manifest lists 2 runs over variants.vcf \
             (default, summary) that document some_class, which --masks names: name the pair's \
             run\n"
        )
    );
    let not_observed = "[\"documented divergences not observed: 1\"]";
    let (status, report, _) = compare(
        "default",
        summary_row(KEPT_SET),
        summary_row(LOST_SET),
        &["--run", "summary"],
    );
    assert_eq!(status, 0);
    assert!(is(
        member(
            &report.expect("a report"),
            &["counts", "documented_divergences_skipped"]
        ),
        "1"
    ));
    let (status, report, _) = compare(
        "default",
        summary_row(KEPT_SET),
        summary_row(OTHER_SET),
        &["--run", "summary"],
    );
    assert_eq!(status, 1);
    assert!(is(
        member(&report.expect("a report"), &["reasons"]),
        not_observed
    ));
}

#[test]
fn a_run_the_manifest_does_not_name_or_of_another_input_is_refused() {
    let dir = scratch("suite-run-refused");
    let manifest = runs_manifest();
    let refused = [
        (
            Some(manifest.as_str()),
            "missing",
            "--run missing: the manifest has no flag set of that name (default names its own run)",
        ),
        (
            Some(manifest.as_str()),
            "elsewhere",
            "--run elsewhere annotates elsewhere.vcf, and --input names variants.vcf",
        ),
        (
            None,
            "default",
            "--run default names a run of a manifest, and none is given",
        ),
    ];
    for (n, (manifest, name, message)) in refused.into_iter().enumerate() {
        let options = with(&RUNS_MASKS, &["--run", name]);
        let (status, report, ran) = compare_fabricated(
            &dir.join(n.to_string()),
            "default",
            &[KEPT_SET, KEPT_SET],
            &[LOST_SET, KEPT_SET],
            manifest,
            &refs(&options),
        );
        assert_eq!(status, 2);
        assert!(report.is_none());
        assert_eq!(ran.err(), format!("{ERROR}{message}\n"));
    }
}

#[test]
fn cli_writes_the_report_and_exits_by_verdict() {
    let dir = scratch("suite-report-and-verdict");
    let default = expected("default.txt.gz");
    let written = dir.join("r.json");
    let ran = run(
        &dir,
        &[
            "--format",
            "default",
            "--expected",
            path_text(&default),
            "--actual",
            path_text(&default),
            "--report",
            path_text(&written),
        ],
    );
    assert_eq!(ran.code, 0);
    let report_green = report(&written);
    assert_eq!(member(&report_green, &["verdict"]).as_str(), Some("green"));
    assert_eq!(member(&report_green, &["mode"]).as_str(), Some("compare"));
    assert_eq!(member(&report_green, &["format"]).as_str(), Some("default"));
    assert_eq!(ran.out(), "default: green\n");

    let mutated = dir.join("mutated.txt");
    let text = text_of(&default).replacen("IMPACT=MODIFIER", "IMPACT=LOW", 1);
    std::fs::write(&mutated, text).expect("a file is planted");
    let ran = run(
        &dir,
        &[
            "--format",
            "default",
            "--expected",
            path_text(&default),
            "--actual",
            path_text(&mutated),
            "--report",
            path_text(&written),
            "--floors",
            "IMPACT=1",
        ],
    );
    assert_eq!(ran.code, 1);
    assert_eq!(
        member(&report(&written), &["verdict"]).as_str(),
        Some("red")
    );
    assert!(ran.out().starts_with("default: red\n"));
    assert!(ran.out().contains("  field IMPACT: 1 differing entries"));
}

#[test]
fn cli_self_parity_and_warnings_modes() {
    let dir = scratch("suite-two-modes");
    let written = dir.join("sp.json");
    let input = corpus().join("variants.vcf");
    let ran = run(
        &dir,
        &[
            "--self-parity",
            "--default",
            path_text(&expected("default.txt.gz")),
            "--tab",
            path_text(&expected("tab.txt.gz")),
            "--vcf",
            path_text(&expected("vcf.vcf.gz")),
            "--json",
            path_text(&expected("json.jsonl.gz")),
            "--input",
            path_text(&input),
            "--report",
            path_text(&written),
        ],
    );
    assert_eq!(ran.code, 0);
    let self_parity = report(&written);
    assert_eq!(
        member(&self_parity, &["mode"]).as_str(),
        Some("self-parity")
    );
    assert_eq!(member(&self_parity, &["verdict"]).as_str(), Some("green"));
    assert_eq!(
        member(&self_parity, &["files", "input"]).as_str(),
        Some(path_text(&input))
    );
    let empty = dir.join("w.txt");
    std::fs::write(&empty, "").expect("a file is planted");
    let ran = run(
        &dir,
        &[
            "--warnings",
            "--expected",
            path_text(&empty),
            "--actual",
            path_text(&empty),
            "--report",
            path_text(&written),
        ],
    );
    assert_eq!(ran.code, 0);
    assert_eq!(
        member(&report(&written), &["format"]).as_str(),
        Some("warnings")
    );
}

#[test]
fn cli_errors_exit_2_with_one_error_line() {
    let dir = scratch("suite-error-lines");
    let vcf = expected("vcf.vcf.gz");
    let ran = run(
        &dir,
        &[
            "--format",
            "vcf",
            "--expected",
            path_text(&vcf),
            "--actual",
            path_text(&vcf),
        ],
    );
    assert_eq!(ran.code, 2);
    assert_eq!(
        ran.err(),
        format!("{ERROR}--format vcf needs --input, the input VCF\n")
    );
    let ran = run(
        &dir,
        &[
            "--self-parity",
            "--default",
            path_text(&expected("default.txt.gz")),
        ],
    );
    assert_eq!(ran.code, 2);
    assert!(ran.err().contains("needs --tab, --vcf, --json"));
    let none = dir.join("none.txt");
    let ran = run(
        &dir,
        &[
            "--format",
            "default",
            "--expected",
            path_text(&none),
            "--actual",
            path_text(&none),
        ],
    );
    assert_eq!(ran.code, 2);
    assert!(ran.err().starts_with(&format!("{ERROR}cannot read ")));
    assert_eq!(run::<&str>(&dir, &[]).code, 2);
}

fn one_shard_cache(root: &Path, shard: &str) -> PathBuf {
    let chromosome = root.join("transcripts/1");
    std::fs::create_dir_all(&chromosome).expect("the cache is made");
    std::fs::write(chromosome.join("1-1000000.json"), shard).expect("a shard is planted");
    root.to_path_buf()
}

#[test]
fn nearest_sites_from_cache_rejects_a_shard_of_another_shape() {
    let dir = scratch("suite-nearest-sites");
    let other_shapes = [
        "\"x\"",
        "7",
        "null",
        "[7]",
        "[\"x\"]",
        "[[]]",
        "{\"1-1000000\": \"x\"}",
        "{\"1-1000000\": {\"a\": 1}}",
    ];
    for (n, shard) in other_shapes.into_iter().enumerate() {
        let cache = one_shard_cache(&dir.join(format!("shape{n}")), shard);
        match NearestSites::from_cache(&cache) {
            Err(error) => assert!(error.to_string().contains("1-1000000.json"), "{shard}"),
            Ok(_) => panic!("the shard {shard} is read"),
        }
    }
    for (name, shard) in [("empty_array", "[]"), ("empty_map", "{}")] {
        let sites = NearestSites::from_cache(&one_shard_cache(&dir.join(name), shard))
            .expect("an empty shard is read");
        assert!(sites.sites.is_empty(), "{name}");
    }
    let default = expected("default.txt.gz");
    let cache = one_shard_cache(&dir.join("cli"), "\"x\"");
    let ran = run(
        &dir,
        &[
            "--format",
            "default",
            "--expected",
            path_text(&default),
            "--actual",
            path_text(&default),
            "--nearest-sites",
            path_text(&cache),
        ],
    );
    assert_eq!(ran.code, 2);
    let err = ran.err();
    assert!(err.starts_with(ERROR), "{err}");
    assert!(err.contains("1-1000000.json"), "{err}");
    assert_eq!(err.matches('\n').count(), 1);
}

#[test]
fn unreadable_inputs_and_an_unwritable_report_exit_2() {
    let dir = scratch("suite-unreadable");
    let default = expected("default.txt.gz");
    let whole = std::fs::read(&default).expect("the reference is read");
    let truncated = dir.join("truncated.txt.gz");
    std::fs::write(&truncated, &whole[..whole.len() / 2]).expect("a file is planted");
    let corrupt = dir.join("corrupt.txt.gz");
    let mut flipped = whole.clone();
    for byte in &mut flipped[2000..2400] {
        *byte ^= 0xff;
    }
    std::fs::write(&corrupt, flipped).expect("a file is planted");
    let binary = dir.join("binary.txt");
    std::fs::write(&binary, &whole).expect("a file is planted");
    for actual in [&truncated, &corrupt, &binary] {
        let ran = run(
            &dir,
            &[
                "--format",
                "default",
                "--expected",
                path_text(&default),
                "--actual",
                path_text(actual),
            ],
        );
        assert_eq!(ran.code, 2);
        let err = ran.err();
        let opening = format!("{ERROR}cannot read {}: ", path_text(actual));
        assert!(err.starts_with(&opening), "{err}");
        assert_eq!(err.matches('\n').count(), 1);
    }
    let under_a_file = binary.join("sub/report.json");
    let ran = run(
        &dir,
        &[
            "--format",
            "default",
            "--expected",
            path_text(&default),
            "--actual",
            path_text(&default),
            "--report",
            path_text(&under_a_file),
        ],
    );
    assert_eq!(ran.code, 2);
    assert!(ran.err().starts_with(&format!("{ERROR}cannot write ")));
}

#[test]
fn json_floors_count_under_the_json_keys() {
    let fields = |column: &str| json_floor_fields(column);
    let named = |names: &[&str]| {
        Some(
            names
                .iter()
                .map(|name| name.to_string())
                .collect::<Vec<_>>(),
        )
    };
    assert_eq!(fields("DOMAINS"), named(&["domains"]));
    assert_eq!(fields("ENSP"), named(&["protein_id"]));
    assert_eq!(fields("SYMBOL"), named(&["gene_symbol"]));
    assert_eq!(fields("SIFT"), named(&["sift_prediction", "sift_score"]));
    assert_eq!(
        fields("PolyPhen"),
        named(&["polyphen_prediction", "polyphen_score"])
    );
    assert_eq!(
        fields("Protein_position"),
        named(&["protein_start", "protein_end"])
    );
    assert_eq!(
        fields("Feature"),
        named(&["transcript_id", "regulatory_feature_id", "motif_feature_id"])
    );
    for column in [
        "Existing_variation",
        "AF",
        "gnomADe_NFE_AF",
        "MAX_AF",
        "FREQS",
        "CLIN_SIG",
        "Feature_type",
        "Extra",
    ] {
        assert_eq!(fields(column), None, "{column}");
    }

    let protein_features = golden().join("GRCh38-protein-features");
    let inputs = input_records(&protein_features.join("variants.vcf")).expect("the input is read");
    let floors: IndexMap<String, i64> = [
        ("DOMAINS", 157),
        ("SIFT", 219),
        ("PolyPhen", 163),
        ("ENSP", 451),
        ("Existing_variation", 1),
        ("Feature_type", 1),
    ]
    .into_iter()
    .map(|(column, floor)| (column.to_string(), floor))
    .collect();
    let compared =
        |fmt: &str, path: &Path, inputs: Option<&[String]>, floors: &IndexMap<String, i64>| {
            compare_files(
                fmt,
                path,
                path,
                inputs,
                &mut Documented::new(),
                &BTreeSet::new(),
                floors,
                &IndexMap::new(),
                EXAMPLES_PER_GROUP,
                false,
                None,
                None,
                None,
                &mut Side::default(),
            )
            .expect("the files compare")
        };
    let json_path = protein_features.join("expected/json.jsonl.gz");
    let report_json = compared("json", &json_path, Some(&inputs), &floors);
    assert_eq!(member(&report_json, &["verdict"]).as_str(), Some("green"));
    assert!(is(member(&report_json, &["reasons"]), "[]"));
    assert!(is(
        member(&report_json, &["floors_not_applicable"]),
        "[\"Existing_variation\", \"Feature_type\"]"
    ));
    let not_applicable = format!(
        "{{\"expected_populated\": 0, \"actual_populated\": 0, \"differing\": 0, \"pairs_accepted\": 0, \
         \"nearest_accepted\": 0, \"order_only\": 0, \"individual_allele_order\": 0, \
         \"motif_score_reverse_strand_mnv\": 0, \"floor\": 1, \"floor_met\": null, \"floor_note\": {}}}",
        quoted(FLOOR_NOT_APPLICABLE)
    );
    assert!(is(
        member(&report_json, &["columns", "Existing_variation"]),
        &not_applicable
    ));
    assert!(is(
        member(&report_json, &["columns", "Feature_type", "floor_met"]),
        "null"
    ));
    assert_eq!(
        member(&report_json, &["columns", "Feature_type", "floor_note"]).as_str(),
        Some(FLOOR_NOT_APPLICABLE)
    );
    for (column, keys) in [
        ("DOMAINS", "[\"domains\"]"),
        ("SIFT", "[\"sift_prediction\", \"sift_score\"]"),
        ("PolyPhen", "[\"polyphen_prediction\", \"polyphen_score\"]"),
        ("ENSP", "[\"protein_id\"]"),
    ] {
        let col = member(&report_json, &["columns", column]);
        let floor = floors[column].to_string();
        assert!(is(member(col, &["floor_met"]), "true"), "{column}");
        assert!(is(member(col, &["floor_fields"]), keys), "{column}");
        assert!(is(member(col, &["expected_populated"]), &floor), "{column}");
        assert!(is(member(col, &["actual_populated"]), &floor), "{column}");
    }
    // The same counts the default format reaches under the column names themselves.
    let default_path = protein_features.join("expected/default.txt.gz");
    let report_text = compared("default", &default_path, None, &floors);
    assert!(is(member(&report_text, &["floors_not_applicable"]), "[]"));
    assert!(is(
        member(&report_text, &["reasons"]),
        "[\"floor Existing_variation=1: expected populated 0, actual populated 0\"]"
    ));
    assert!(is(
        member(&report_text, &["columns", "Feature_type", "floor_met"]),
        "true"
    ));
    assert!(member(&report_text, &["columns", "DOMAINS"])
        .get("floor_fields")
        .is_none());
    assert!(is(
        member(&report_text, &["columns", "DOMAINS", "expected_populated"]),
        "157"
    ));

    // A floor the JSON counts fall short of is red, as in the text formats.
    let mut one = IndexMap::new();
    one.insert("DOMAINS".to_string(), 158);
    let short = compared("json", &json_path, Some(&inputs), &one);
    assert_eq!(member(&short, &["verdict"]).as_str(), Some("red"));
    assert!(is(
        member(&short, &["reasons"]),
        "[\"floor DOMAINS=158: expected populated 157, actual populated 157\"]"
    ));

    // `--sift s` prints the bare score, which the JSON format carries as sift_score alone.
    let score_corpus = golden().join("GRCh38-sift-score");
    let score_inputs =
        input_records(&score_corpus.join("variants.vcf")).expect("the input is read");
    let score_floors: IndexMap<String, i64> =
        [("SIFT".to_string(), 84), ("PolyPhen".to_string(), 58)]
            .into_iter()
            .collect();
    let score_report = compared(
        "json",
        &score_corpus.join("expected/json.jsonl.gz"),
        Some(&score_inputs),
        &score_floors,
    );
    assert_eq!(member(&score_report, &["verdict"]).as_str(), Some("green"));
    assert!(is(
        member(&score_report, &["columns", "SIFT", "expected_populated"]),
        "84"
    ));
    assert!(is(
        member(
            &score_report,
            &["columns", "sift_score", "expected_populated"]
        ),
        "84"
    ));
    assert!(member(&score_report, &["columns"])
        .get("sift_prediction")
        .is_none());

    // Through the command line, the floors of a fixture row on its JSON format.
    let dir = scratch("suite-json-floors");
    let out = dir.join("json.json");
    let ran = run(
        &dir,
        &[
            "--format",
            "json",
            "--expected",
            path_text(&json_path),
            "--actual",
            path_text(&json_path),
            "--input",
            path_text(&protein_features.join("variants.vcf")),
            "--floors",
            "DOMAINS=157,SIFT=219,PolyPhen=163,ENSP=451,Existing_variation=1",
            "--masks",
            "order_only",
            "--report",
            path_text(&out),
        ],
    );
    assert_eq!(ran.code, 0);
    let written = report(&out);
    assert!(is(
        member(&written, &["floors_not_applicable"]),
        "[\"Existing_variation\"]"
    ));
    assert!(is(
        member(&written, &["masks", "unused"]),
        "[\"order_only\"]"
    ));
}

// The recoder format: two variant_recoder outputs, the reference's keys and values in hash
// order and with the RefSeq entries of its merged database, the candidate's sorted and from
// an Ensembl cache.
const RECODER_EXPECTED: &str = r#"[
 {"T": {
  "vcf_string": ["21-25585733-C-T"],
  "hgvsg": ["NC_000021.9:g.25585733C>T"],
  "hgvsc": ["ENST00000352957.9:c.991G>A", "NM_017446.4:c.991G>A", "ENST00000307301.12:c.*18G>A"],
  "hgvsp": ["NP_059142.3:p.Ala331Thr", "ENSP00000284967.7:p.Ala331Thr", "LRG_1p1:p.Ala331Thr"],
  "ga4gh_vrs": [{
   "state": {"sequence": "T", "type": "SequenceState"},
   "location": {
    "interval": {"end": 25585733, "start": "25585732", "type": "SimpleInterval"},
    "sequence_id": "refseq:NC_000021.9",
    "type": "SequenceLocation"},
   "type": "Allele"}],
  "input": "rs142513484",
  "spdi": ["NC_000021.9:25585732:C:T"],
  "mane_select": [
   {"hgvsc": "ENST00000352957.9:c.991G>A", "hgvsg": "NC_000021.9:g.25585733C>T", "hgvsp": "ENSP00000284967.7:p.Ala331Thr"},
   {"hgvsg": "NC_000021.9:g.25585733C>T", "hgvsc": "NM_017446.4:c.991G>A", "hgvsp": "NP_059142.3:p.Ala331Thr"}],
  "var_synonyms": ["ClinVar: VCV000829541,RCV003317418", "OMIM: 606870.0005"],
  "id": ["rs142513484"]}},
 null,
 {"warnings": ["No variant found with ID 'rs_nope'"]}
]"#;
const RECODER_ACTUAL: &str = r#"[
 {"T": {
  "ga4gh_vrs": [{
   "location": {
    "interval": {"end": 25585733, "start": "25585732", "type": "SimpleInterval"},
    "sequence_id": "refseq:NC_000021.9",
    "type": "SequenceLocation"},
   "state": {"sequence": "T", "type": "SequenceState"},
   "type": "Allele"}],
  "hgvsc": ["ENST00000307301.12:c.*18G>A", "ENST00000352957.9:c.991G>A"],
  "hgvsg": ["NC_000021.9:g.25585733C>T"],
  "hgvsp": ["ENSP00000284967.7:p.Ala331Thr"],
  "id": ["rs142513484"],
  "input": "rs142513484",
  "mane_select": [
   {"hgvsc": "ENST00000352957.9:c.991G>A", "hgvsg": "NC_000021.9:g.25585733C>T", "hgvsp": "ENSP00000284967.7:p.Ala331Thr"}],
  "spdi": ["NC_000021.9:25585732:C:T"],
  "var_synonyms": ["ClinVar: RCV003317418,VCV000829541", "OMIM: 606870.0005"],
  "vcf_string": ["21-25585733-C-T"]}},
 null,
 {"warnings": ["No variant found with ID 'rs_nope'"]}
]"#;
const RECODER_MASKS: &str = "refseq_transcript_entries,lrg_transcript_entries";

/// One step into a JSON value: a member's name or an item's index.
#[derive(Clone, Copy)]
enum Step {
    Name(&'static str),
    At(usize),
}
use Step::{At, Name};

fn node<'a>(value: &'a mut Value, path: &[Step]) -> &'a mut Value {
    path.iter().fold(value, |value, step| match (value, step) {
        (Value::Dict(members), Name(name)) => members.get_mut(*name).expect("the member"),
        (Value::List(items), At(index)) => &mut items[*index],
        _ => panic!("the step does not fit the value"),
    })
}

fn push(value: &mut Value, path: &[Step], item: &str) {
    match node(value, path) {
        Value::List(items) => items.push(json(item)),
        _ => panic!("the value is not a list"),
    }
}

fn set(value: &mut Value, path: &[Step], name: &str, member: Value) {
    match node(value, path) {
        Value::Dict(members) => {
            members.insert(name.to_string(), member);
        }
        _ => panic!("the value is not an object"),
    }
}

const FIRST_ALLELE: [Step; 2] = [At(0), Name("T")];

/// The report of the recoder comparison of `expected` and `actual` under `masks`.
fn compare_recoder(dir: &Path, expected: &Value, actual: &Value, masks: &str) -> Value {
    let plant = |name: &str, objects: &Value| {
        let layout = DumpOptions {
            indent: Some(3),
            sort_keys: false,
            ensure_ascii: true,
        };
        std::fs::write(dir.join(name), pyjson::dumps(objects, &layout)).expect("a file is planted");
    };
    plant("expected.json", expected);
    plant("actual.json", actual);
    run(
        dir,
        &[
            "--format",
            "recoder",
            "--expected",
            "expected.json",
            "--actual",
            "actual.json",
            "--masks",
            masks,
            "--report",
            "report.json",
        ],
    );
    report(&dir.join("report.json"))
}

fn red_with(report: &Value, reasons: &[&str]) {
    assert_eq!(member(report, &["verdict"]).as_str(), Some("red"));
    let given = strings(member(report, &["reasons"]));
    for reason in reasons {
        assert!(given.contains(reason), "{reason} is not among {given:?}");
    }
}

fn green(report: &Value) {
    assert_eq!(
        member(report, &["verdict"]).as_str(),
        Some("green"),
        "{}",
        compact(member(report, &["reasons"]))
    );
}

#[test]
fn recoder_format_compares_each_alleles_keys_as_sets_under_the_masks() {
    let dir = scratch("suite-recoder-sets");
    let report = compare_recoder(
        &dir,
        &json(RECODER_EXPECTED),
        &json(RECODER_ACTUAL),
        RECODER_MASKS,
    );
    green(&report);
    assert_eq!(member(&report, &["format"]).as_str(), Some("recoder"));
    assert!(is(
        member(&report, &["masks", "applied"]),
        "{\"lrg_transcript_entries\": 1, \"refseq_transcript_entries\": 3}"
    ));
    assert!(is(member(&report, &["masks", "unused"]), "[]"));
}

#[test]
fn recoder_format_names_a_differing_key_and_a_missing_object() {
    let dir = scratch("suite-recoder-differing");
    let reference = json(RECODER_EXPECTED);
    let mut actual = json(RECODER_ACTUAL);
    set(
        &mut actual,
        &FIRST_ALLELE,
        "hgvsg",
        json("[\"NC_000021.9:g.25585733C>A\"]"),
    );
    push(
        &mut actual,
        &[At(0), Name("T"), Name("hgvsc")],
        "\"ENST00000999999.1:c.1A>G\"",
    );
    let report = compare_recoder(&dir, &reference, &actual, RECODER_MASKS);
    red_with(
        &report,
        &[
            "field hgvsg: 1 differing entries",
            "field hgvsc: 1 differing entries",
        ],
    );

    let objects = |value: &Value| value.as_list().expect("the objects").to_vec();
    let first_alone = Value::List(objects(&json(RECODER_ACTUAL))[..1].to_vec());
    let report = compare_recoder(&dir, &reference, &first_alone, RECODER_MASKS);
    assert_eq!(member(&report, &["verdict"]).as_str(), Some("red"));
    let reasons = strings(member(&report, &["reasons"]));
    assert!(
        reasons
            .iter()
            .any(|reason| reason.starts_with("missing in the candidate: 2 entries")),
        "{reasons:?}"
    );

    let mut swapped = objects(&reference);
    swapped.swap(0, 1);
    let report = compare_recoder(&dir, &reference, &Value::List(swapped), RECODER_MASKS);
    assert_eq!(member(&report, &["verdict"]).as_str(), Some("red"));
    let reasons = strings(member(&report, &["reasons"]));
    assert!(
        reasons
            .iter()
            .any(|reason| reason.starts_with("order differs")),
        "{reasons:?}"
    );

    // The alleles of one object are hash keys in variant_recoder's output: any order is green.
    let allele = member(&reference, &[]).as_list().expect("the objects")[0]
        .get("T")
        .expect("the allele")
        .clone();
    let pair = |first: &str, second: &str| {
        let mut alleles = IndexMap::new();
        alleles.insert(first.to_string(), allele.clone());
        alleles.insert(second.to_string(), allele.clone());
        Value::List(vec![Value::Dict(alleles)])
    };
    let report = compare_recoder(&dir, &pair("G", "T"), &pair("T", "G"), RECODER_MASKS);
    green(&report);
}

#[test]
fn recoder_format_lrg_class_drops_the_genomic_entries_of_the_lrg_region() {
    let dir = scratch("suite-recoder-lrg");
    let mut reference = json(RECODER_EXPECTED);
    push(
        &mut reference,
        &[At(0), Name("T"), Name("hgvsg")],
        "\"LRG_1:g.12345C>T\"",
    );
    push(
        &mut reference,
        &[At(0), Name("T"), Name("spdi")],
        "\"LRG_1:12344:C:T\"",
    );
    push(
        &mut reference,
        &[At(0), Name("T"), Name("vcf_string")],
        "\"LRG_1-12345-C-T\"",
    );
    let report = compare_recoder(&dir, &reference, &json(RECODER_ACTUAL), RECODER_MASKS);
    green(&report);
    assert!(is(
        member(&report, &["masks", "applied", "lrg_transcript_entries"]),
        "4"
    ));
    let report = compare_recoder(
        &dir,
        &reference,
        &json(RECODER_ACTUAL),
        "refseq_transcript_entries",
    );
    red_with(
        &report,
        &[
            "field hgvsg: 1 differing entries",
            "field vcf_string: 1 differing entries",
        ],
    );
}

#[test]
fn recoder_format_without_the_masks_sees_the_refseq_entries() {
    let dir = scratch("suite-recoder-no-masks");
    let report = compare_recoder(&dir, &json(RECODER_EXPECTED), &json(RECODER_ACTUAL), "");
    red_with(
        &report,
        &[
            "field hgvsc: 1 differing entries",
            "field hgvsp: 1 differing entries",
            "field mane_select: 1 differing entries",
        ],
    );
}

#[test]
fn recoder_format_archive_dbsnp_synonyms_class_drops_the_merged_identifiers() {
    let dir = scratch("suite-recoder-archive");
    let archive = "\"Archive dbSNP: rs60296503,rs17181160\"";
    let with_archive = format!("{RECODER_MASKS},archive_dbsnp_synonyms");
    let mut reference = json(RECODER_EXPECTED);
    push(
        &mut reference,
        &[At(0), Name("T"), Name("var_synonyms")],
        archive,
    );
    let report = compare_recoder(&dir, &reference, &json(RECODER_ACTUAL), &with_archive);
    green(&report);
    assert!(is(
        member(&report, &["masks", "applied", "archive_dbsnp_synonyms"]),
        "1"
    ));
    let report = compare_recoder(&dir, &reference, &json(RECODER_ACTUAL), RECODER_MASKS);
    red_with(&report, &["field var_synonyms: 1 differing entries"]);
    set(
        &mut reference,
        &FIRST_ALLELE,
        "var_synonyms",
        json(&format!("[{archive}]")),
    );
    let mut actual = json(RECODER_ACTUAL);
    match node(&mut actual, &FIRST_ALLELE) {
        Value::Dict(members) => {
            members.shift_remove("var_synonyms");
        }
        _ => panic!("the allele is an object"),
    }
    let report = compare_recoder(&dir, &reference, &actual, &with_archive);
    green(&report);
}

#[test]
fn recoder_format_identifier_indel_class_drops_the_transcript_notations() {
    let dir = scratch("suite-recoder-indel");
    let deletion = r#"{
     "input": "rs2057997486",
     "hgvsg": ["NC_000021.9:g.34886862del"],
     "hgvsc": ["ENST00000300305.7:c.332del"],
     "hgvsp": ["ENSP00000300305.3:p.Leu112CysfsTer10"],
     "mane_select": [{"hgvsc": "ENST00000675419.1:c.332del", "hgvsg": "NC_000021.9:g.34886862del", "hgvsp": "ENSP00000501943.1:p.Leu112CysfsTer10"}]}"#;
    let under_the_deleted_allele = |allele: &Value| {
        let mut alleles = IndexMap::new();
        alleles.insert("-".to_string(), allele.clone());
        Value::List(vec![Value::Dict(alleles)])
    };
    let mut shifted = json(deletion);
    set(
        &mut shifted,
        &[],
        "hgvsc",
        json("[\"ENST00000300305.7:c.334del\"]"),
    );
    set(
        &mut shifted,
        &[Name("mane_select"), At(0)],
        "hgvsc",
        json("\"ENST00000675419.1:c.334del\""),
    );
    let masks = format!("{RECODER_MASKS},identifier_indel_transcript_notations");
    let reference = under_the_deleted_allele(&json(deletion));
    let report = compare_recoder(
        &dir,
        &reference,
        &under_the_deleted_allele(&shifted),
        &masks,
    );
    green(&report);
    assert!(is(
        member(
            &report,
            &["masks", "applied", "identifier_indel_transcript_notations"]
        ),
        "6"
    ));
    let report = compare_recoder(
        &dir,
        &reference,
        &under_the_deleted_allele(&shifted),
        RECODER_MASKS,
    );
    red_with(
        &report,
        &[
            "field hgvsc: 1 differing entries",
            "field mane_select: 1 differing entries",
        ],
    );
    let notation = json("\"NC_000021.9:g.34886862del\"");
    let mut by_notation = json(deletion);
    set(&mut by_notation, &[], "input", notation.clone());
    set(&mut shifted, &[], "input", notation);
    let report = compare_recoder(
        &dir,
        &under_the_deleted_allele(&by_notation),
        &under_the_deleted_allele(&shifted),
        &masks,
    );
    red_with(&report, &["field hgvsc: 1 differing entries"]);
    let mut substitution = json(RECODER_ACTUAL);
    set(
        &mut substitution,
        &FIRST_ALLELE,
        "hgvsc",
        json("[\"ENST00000352957.9:c.991G>T\"]"),
    );
    let report = compare_recoder(&dir, &json(RECODER_EXPECTED), &substitution, &masks);
    red_with(&report, &["field hgvsc: 1 differing entries"]);

    let indel = |label: &str, allele: &str, fields: &str| {
        recoder_identifier_indel(label, allele, &json(fields))
    };
    assert!(indel("rs1", "-", "{}"));
    assert!(indel("rs1", "T", "{\"hgvsg\": [\"NC_000021.9:g.100dup\"]}"));
    assert!(indel("rs1", "T", "{\"spdi\": [\"NC_000021.9:99:A:AT\"]}"));
    assert!(indel("rs1", "T", "{\"vcf_string\": [\"21-99-A-AT\"]}"));
    assert!(!indel(
        "rs1",
        "T",
        "{\"hgvsg\": [\"NC_000021.9:g.100A>T\"], \"spdi\": [\"NC_000021.9:99:A:T\"], \"vcf_string\": [\"21-100-A-T\"]}"
    ));
    assert!(!indel(
        "rs1",
        "T",
        "{\"hgvsg\": \"NC_000021.9:g.100dup\", \"spdi\": 5, \"vcf_string\": \"21-99-A-AT\"}"
    ));
    assert!(!indel(
        "rs1",
        "T",
        "{\"spdi\": [null, 7, \"NC_000021.9:99:A\"], \"vcf_string\": [[\"21-99-A-AT\"]]}"
    ));
    assert!(!indel("NC_000021.9:g.100dup", "-", "{}"));
    assert!(!indel("21\t100\trs1\tA\tAT\t.\t.\t.", "-", "{}"));
}

#[test]
fn recoder_format_compares_ga4gh_vrs_start_as_written() {
    let dir = scratch("suite-recoder-vrs");
    let mut actual = json(RECODER_ACTUAL);
    set(
        &mut actual,
        &[
            At(0),
            Name("T"),
            Name("ga4gh_vrs"),
            At(0),
            Name("location"),
            Name("interval"),
        ],
        "start",
        json("25585732"),
    );
    let report = compare_recoder(&dir, &json(RECODER_EXPECTED), &actual, RECODER_MASKS);
    red_with(&report, &["field ga4gh_vrs: 1 differing entries"]);
}

/// One transcript row of record `record`.
fn planted_row(record: usize, feature: &str, extra: Option<&str>) -> String {
    let own = format!("IMPACT=MODIFIER;STRAND=1;DOMAINS=db:a{record},db:b{record}");
    default_row(
        &format!("rec{record}"),
        &format!("21:{}", 1000 + record),
        &"ACGT"[record % 4..record % 4 + 1],
        feature,
        "intron_variant",
        extra.unwrap_or(&own),
    )
}

/// Two transcript rows for each of eight records.
fn planted_rows() -> Vec<String> {
    (0..8)
        .flat_map(|record| ["ENST1", "ENST2"].map(|feature| planted_row(record, feature, None)))
        .collect()
}

#[test]
fn an_output_whose_path_is_spelt_as_the_help_option_is_read_as_an_output() {
    let dir = scratch("suite-help-spelt-paths");
    let rows = planted_rows();
    let mut other = rows.clone();
    other[6] = planted_row(3, "ENST1", Some("IMPACT=LOW;STRAND=1;DOMAINS=db:a3,db:b3"));
    std::fs::write(dir.join("-h"), default_text(&rows)).expect("a file is planted");
    std::fs::write(dir.join("--help"), default_text(&other)).expect("a file is planted");
    let args = [
        "--format",
        "default",
        "--expected=-h",
        "--actual=--help",
        "--report",
        "whole.json",
    ];
    let ran = run(&dir, &args);
    let lines =
        "default: red\n  field Extra: 1 differing entries\n  field IMPACT: 1 differing entries\n";
    assert_eq!((ran.code, ran.out().as_str()), (1, lines));
    assert!(dir.join("whole.json").is_file());
}
