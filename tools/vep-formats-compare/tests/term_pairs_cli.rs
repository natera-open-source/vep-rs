// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The stop codon class a manifest also documents, the reference bases read from the
//! input's records and two rows the stream hands to the comparison, each through the
//! binary on planted files.

mod compare_common;
mod support;
mod term_pairs_common;

use std::path::Path;

use support::args_of;
use support::compact;
use support::member;
use support::report;
use support::run;
use support::scratch;
use support::Ran;
use term_pairs_common::*;
use vep_formats_compare::consts::WITH_REFERENCE_SETS;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;

/// The fabricated rows of the manifest-reading tests: record `i` is one deletion on one
/// transcript, written in every text format with the consequence set its side gives it.
/// The reference writes `KEPT_SET`; the candidate writes `LOST_SET` in a run without a
/// reference sequence.
const KEPT_SET: &str = "3_prime_UTR_variant,stop_retained_variant";
const LOST_SET: &str = "3_prime_UTR_variant,stop_lost";
const OTHER_SET: &str = "3_prime_UTR_variant,coding_sequence_variant";
const WITH_REFERENCE: &str = "vep_rs_consequence_sets_with_reference_sequence";
const REMOVED_CLASS: &str = "stop_codon_removed_without_stop_codon_test";
const REFERENCE_MASKS: [&str; 2] = ["--masks", REMOVED_CLASS];

fn impact(set: &str) -> &'static str {
    match set {
        KEPT_SET => "LOW",
        LOST_SET => "HIGH",
        _ => "MODIFIER",
    }
}

fn json(text: &str) -> Value {
    pyjson::loads(text).unwrap_or_else(|error| panic!("{text} is not JSON: {error}"))
}

/// Whether `value` is the value `expected` spells, the order of an object's members aside.
fn is(value: &Value, expected: &str) -> bool {
    pyjson::py_eq(value, &json(expected))
}

fn path_text(path: &Path) -> &str {
    path.to_str().expect("the path is UTF-8")
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

/// The fabricated record `ordinal` as a line of the Ensembl default input format.
fn fabricated_ensembl_line(ordinal: usize, separator: &str) -> String {
    let position = 100 * (ordinal + 1);
    [
        "21".to_string(),
        position.to_string(),
        (position + 1).to_string(),
        "AT/-".to_string(),
        "+".to_string(),
    ]
    .join(separator)
}

/// The `divergences` item of the fabricated record `ordinal` under the class `class`, its
/// record ordinals `records`, with `more` members after its own.
fn fabricated_item_of(ordinal: usize, class: &str, records: &str, more: &str) -> String {
    let (_, location, feature) = fabricated_record(ordinal);
    format!(
        "{{\"location\": {}, \"allele\": \"-\", \"feature\": {}, \"feature_type\": \"Transcript\", \"record_indices\": {records}, \"vep_consequence_set\": {}, \"vep_rs_consequence_sets\": [{}], \"expected_divergence\": {}{more}}}",
        quoted(&location),
        quoted(&feature),
        quoted(KEPT_SET),
        quoted(LOST_SET),
        quoted(class)
    )
}

fn fabricated_item(ordinal: usize, class: &str, more: &str) -> String {
    fabricated_item_of(ordinal, class, &format!("[{ordinal}]"), more)
}

/// A manifest of the `divergences` items `items`.
fn divergences(items: &[String]) -> String {
    format!("{{\"divergences\": [{}]}}", items.join(", "))
}

/// The member that states an item's sets of a run with a reference sequence, `value` its
/// JSON text.
fn with_reference(value: &str) -> String {
    format!(", {}: {value}", quoted(WITH_REFERENCE))
}

/// A manifest whose one item carries the reference's own set as its with-reference sets.
fn reference_manifest() -> String {
    divergences(&[fabricated_item(
        0,
        REMOVED_CLASS,
        &with_reference(&quoted_list(&[KEPT_SET])),
    )])
}

/// One row per fabricated record in the format `fmt`, record `i` carrying `sets[i]`;
/// `echoed` is the input line each JSON object echoes where it is not the record's VCF
/// line.
fn fabricated_output(fmt: &str, sets: &[&str], echoed: Option<&[String]>) -> String {
    let records: Vec<_> = (0..sets.len()).map(fabricated_record).collect();
    let rows = records.iter().zip(sets);
    let mut text = String::new();
    match fmt {
        "default" => {
            text.push_str(compare_common::DEFAULT_HEADER);
            text.push('\n');
            for (i, ((_, location, feature), terms)) in rows.enumerate() {
                text.push_str(&compare_common::default_row(
                    &format!("del{i}"),
                    location,
                    "-",
                    feature,
                    terms,
                    &format!("IMPACT={}", impact(terms)),
                ));
                text.push('\n');
            }
        }
        "tab" => {
            text.push_str("#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIMPACT\n");
            for (i, ((_, location, feature), terms)) in rows.enumerate() {
                text.push_str(&format!(
                    "del{i}\t{location}\t-\tENSG1\t{feature}\tTranscript\t{terms}\t{}\n",
                    impact(terms)
                ));
            }
        }
        "vcf" => {
            text.push_str(&csq_header(
                "Allele|Consequence|IMPACT|Feature_type|Feature",
            ));
            for ((line, _, feature), terms) in rows {
                text.push_str(&format!(
                    "{}CSQ=-|{}|{}|Transcript|{feature}\n",
                    open_info(line),
                    terms.replace(',', "&"),
                    impact(terms)
                ));
            }
        }
        "json" => {
            for (i, ((line, _, feature), terms)) in rows.enumerate() {
                let term_list: Vec<&str> = terms.split(',').collect();
                let input = echoed.map_or(line, |echoed| &echoed[i]);
                text.push_str(&format!(
                    "{{\"input\": {}, \"id\": \"del{i}\", \"allele_string\": \"AT/-\", \"most_severe_consequence\": {}, \"transcript_consequences\": [{{\"variant_allele\": \"-\", \"transcript_id\": {}, \"consequence_terms\": {}, \"impact\": {}}}]}}\n",
                    quoted(input),
                    quoted(term_list[0]),
                    quoted(feature),
                    quoted_list(&term_list),
                    quoted(impact(terms))
                ));
            }
        }
        other => panic!("no fabricated output in the format {other}"),
    }
    text
}

/// What a fabricated pair is compared beside: the input where it is not the records' VCF
/// lines, what the JSON objects then echo of it, and the input's file name.
#[derive(Clone, Copy)]
struct Beside<'a> {
    input_lines: Option<&'a [String]>,
    echoed: Option<&'a [String]>,
    input_name: &'a str,
}

const RECORD_LINES: Beside = Beside {
    input_lines: None,
    echoed: None,
    input_name: "variants.vcf",
};

/// The exit code, the report where one was written and the run of the comparator on the
/// fabricated pair in the format `fmt`: the reference carrying the sets `expected` and the
/// candidate `actual`, read beside `manifest` under `options`.
fn compare_fabricated_beside(
    dir: &Path,
    fmt: &str,
    expected: &[&str],
    actual: &[&str],
    manifest: Option<&str>,
    options: &[&str],
    beside: Beside,
) -> (i32, Option<Value>, Ran) {
    std::fs::create_dir_all(dir).expect("the directory is made");
    let plant = |name: &str, text: String| {
        let path = dir.join(name);
        std::fs::write(&path, text).expect("a file is planted");
        path
    };
    let reference = plant("vep.out", fabricated_output(fmt, expected, beside.echoed));
    let candidate = plant(
        "candidate.out",
        fabricated_output(fmt, actual, beside.echoed),
    );
    let lines: String = match beside.input_lines {
        Some(lines) => lines.iter().map(|line| format!("{line}\n")).collect(),
        None => (0..expected.len())
            .map(|i| format!("{}\n", fabricated_record(i).0))
            .collect(),
    };
    let input = plant(beside.input_name, lines);
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

/// `compare_fabricated_beside` the records' own VCF lines, with the report it must write.
fn compare_fabricated(
    dir: &Path,
    fmt: &str,
    expected: &[&str],
    actual: &[&str],
    manifest: &str,
    options: &[&str],
) -> (i32, Value) {
    let (status, report, ran) = compare_fabricated_beside(
        dir,
        fmt,
        expected,
        actual,
        Some(manifest),
        options,
        RECORD_LINES,
    );
    let report = report.unwrap_or_else(|| panic!("no report: {}", ran.err()));
    (status, report)
}

fn options(more: &[&'static str]) -> Vec<&'static str> {
    REFERENCE_MASKS.iter().chain(more).copied().collect()
}

fn reasons(report: &Value) -> String {
    compact(member(report, &["reasons"]))
}

fn applied(report: &Value) -> &Value {
    member(report, &["masks", "applied"])
}

fn skipped(report: &Value) -> &Value {
    member(report, &["counts", "documented_divergences_skipped"])
}

fn removed(count: i64) -> String {
    format!("{{{}: {count}}}", quoted(REMOVED_CLASS))
}

const NOT_OBSERVED: &str = "[\"documented divergences not observed: 1\"]";

/// The red reason of one entry whose Consequence differs, as the format names the field.
fn consequence_reason(fmt: &str) -> String {
    let field = if fmt == "json" {
        "consequence_terms"
    } else {
        "Consequence"
    };
    format!("field {field}: 1 differing entries")
}

/// Crossing 1: under `--reference-sequence absent` an item that carries the with-reference
/// key documents its `vep_rs_consequence_sets`, and a candidate that wrote them is green
/// with the item counted.
fn a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item(
    fmt: &str,
) {
    let dir = scratch(&format!("term-pairs-crossing-1-{fmt}"));
    let (status, report) = compare_fabricated(
        &dir,
        fmt,
        &[KEPT_SET],
        &[LOST_SET],
        &reference_manifest(),
        &options(&["--reference-sequence", "absent"]),
    );
    assert_eq!(status, 0, "{}", reasons(&report));
    assert_eq!(member(&report, &["verdict"]).as_str(), Some("green"));
    assert!(is(applied(&report), &removed(1)));
    assert!(is(skipped(&report), "1"));
    assert_eq!(
        member(&report, &["reference_sequence"]).as_str(),
        Some("absent")
    );
}

/// Crossing 2: under `--reference-sequence present` an item whose with-reference sets are
/// the reference's own documents nothing: the row is compared like any other, no difference
/// is seen, and the item is neither skipped as documented nor reported as not observed.
fn a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly(fmt: &str) {
    let dir = scratch(&format!("term-pairs-crossing-2-{fmt}"));
    let (status, report) = compare_fabricated(
        &dir,
        fmt,
        &[KEPT_SET],
        &[KEPT_SET],
        &reference_manifest(),
        &options(&["--reference-sequence", "present"]),
    );
    assert_eq!(status, 0, "{}", reasons(&report));
    assert_eq!(member(&report, &["verdict"]).as_str(), Some("green"));
    assert!(is(member(&report, &["stale_divergences"]), "[]"));
    assert!(is(skipped(&report), "0"));
    let compared = if matches!(fmt, "default" | "tab") {
        "1"
    } else {
        "2"
    };
    assert!(is(member(&report, &["counts", "compared"]), compared));
    assert_eq!(
        member(&report, &["reference_sequence"]).as_str(),
        Some("present")
    );
}

/// Crossing 3: a run without a reference sequence that wrote the reference's set leaves the
/// item not observed.
fn a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red(fmt: &str) {
    let dir = scratch(&format!("term-pairs-crossing-3-{fmt}"));
    let (status, report) = compare_fabricated(
        &dir,
        fmt,
        &[KEPT_SET],
        &[KEPT_SET],
        &reference_manifest(),
        &options(&["--reference-sequence", "absent"]),
    );
    assert_eq!((status, reasons(&report).as_str()), (1, NOT_OBSERVED));
}

/// Crossing 4: a run with a reference sequence that wrote the sets of a run without one is
/// a plain difference of the Consequence. The item's class is a built-in row class too,
/// whose predicate admits the pair; it takes no row of a key the manifest states as equal.
fn a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red(fmt: &str) {
    let dir = scratch(&format!("term-pairs-crossing-4-{fmt}"));
    let present = options(&["--reference-sequence", "present"]);
    let (status, report) = compare_fabricated(
        &dir,
        fmt,
        &[KEPT_SET],
        &[LOST_SET],
        &reference_manifest(),
        &present,
    );
    assert_eq!(status, 1);
    let stated: Vec<&str> = member(&report, &["reasons"])
        .as_list()
        .expect("the reasons")
        .iter()
        .filter_map(Value::as_str)
        .collect();
    assert!(
        stated.contains(&consequence_reason(fmt).as_str()),
        "{stated:?}"
    );
    assert!(is(applied(&report), "{}"));
    assert!(is(member(&report, &["stale_divergences"]), "[]"));
    // The same pair with no item on the key is the class's own: the predicate admits it.
    let (status, report) = compare_fabricated(
        &dir.join("unstated"),
        fmt,
        &[KEPT_SET],
        &[LOST_SET],
        &divergences(&[]),
        &present,
    );
    assert_eq!(status, 0, "{}", reasons(&report));
    assert!(is(applied(&report), &removed(1)));
}

/// A call that states no setting reads an item's `vep_rs_consequence_sets` whether or not
/// the item carries the with-reference key, and an item without the key reads the same under
/// both settings: a candidate that wrote the documented set is green with the item counted,
/// one that wrote the reference's set leaves it not observed.
fn an_item_without_the_key_and_a_call_without_the_setting_read_alike(fmt: &str) {
    let dir = scratch(&format!("term-pairs-no-key-no-setting-{fmt}"));
    let plain = divergences(&[fabricated_item(0, REMOVED_CLASS, "")]);
    let observed = |name: String, manifest: &str, options: &[&str], context: &str| {
        let (status, report) = compare_fabricated(
            &dir.join(name),
            fmt,
            &[KEPT_SET],
            &[LOST_SET],
            manifest,
            options,
        );
        assert_eq!(status, 0, "{context}: {}", reasons(&report));
        assert!(is(applied(&report), &removed(1)), "{context}");
    };
    let unobserved = |name: String, manifest: &str, options: &[&str], context: &str| {
        let (status, report) = compare_fabricated(
            &dir.join(name),
            fmt,
            &[KEPT_SET],
            &[KEPT_SET],
            manifest,
            options,
        );
        assert_eq!(
            (status, reasons(&report).as_str()),
            (1, NOT_OBSERVED),
            "{context}"
        );
    };
    for (n, manifest) in [reference_manifest(), plain.clone()].iter().enumerate() {
        let context = format!("reading {n}");
        observed(format!("observed{n}"), manifest, &REFERENCE_MASKS, &context);
        unobserved(
            format!("unobserved{n}"),
            manifest,
            &REFERENCE_MASKS,
            &context,
        );
    }
    for (n, setting) in ["present", "absent"].into_iter().enumerate() {
        let options = options(&["--reference-sequence", setting]);
        observed(format!("told{n}"), &plain, &options, setting);
        unobserved(format!("told-unobserved{n}"), &plain, &options, setting);
    }
}

/// One test a text format for each of the five readings of the with-reference key.
macro_rules! per_format {
    ($($(#[$waits:meta])? $name:ident: $case:ident($fmt:literal);)*) => {
        $(
            #[test]
            $(#[$waits])?
            fn $name() {
                $case($fmt);
            }
        )*
    };
}

per_format! {
    a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item_in_the_default_format: a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item("default");
    a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item_in_the_tab_format: a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item("tab");
    a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item_in_the_vcf_format: a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item("vcf");
    a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item_in_the_json_format: a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item("json");
    a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly_in_the_default_format: a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly("default");
    a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly_in_the_tab_format: a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly("tab");
    a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly_in_the_vcf_format: a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly("vcf");
    a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly_in_the_json_format: a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly("json");
    a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red_in_the_default_format: a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red("default");
    a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red_in_the_tab_format: a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red("tab");
    a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red_in_the_vcf_format: a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red("vcf");
    a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red_in_the_json_format: a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red("json");
    a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red_in_the_default_format: a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red("default");
    a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red_in_the_tab_format: a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red("tab");
    a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red_in_the_vcf_format: a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red("vcf");
    a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red_in_the_json_format: a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red("json");
    an_item_without_the_key_and_a_call_without_the_setting_read_alike_in_the_default_format: an_item_without_the_key_and_a_call_without_the_setting_read_alike("default");
    an_item_without_the_key_and_a_call_without_the_setting_read_alike_in_the_tab_format: an_item_without_the_key_and_a_call_without_the_setting_read_alike("tab");
    an_item_without_the_key_and_a_call_without_the_setting_read_alike_in_the_vcf_format: an_item_without_the_key_and_a_call_without_the_setting_read_alike("vcf");
    an_item_without_the_key_and_a_call_without_the_setting_read_alike_in_the_json_format: an_item_without_the_key_and_a_call_without_the_setting_read_alike("json");
}

const TEXT_FORMATS: [&str; 4] = ["default", "tab", "vcf", "json"];

/// The applied count of a class is the entries it sets aside. Of two items of one class,
/// the one whose with-reference sets are the reference's own sets nothing aside in a run
/// with a reference sequence: the class is applied to one entry there and to two in a run
/// without one, while both items stay available.
#[test]
fn an_item_stated_equal_with_a_reference_sequence_is_not_counted_as_applied() {
    let dir = scratch("term-pairs-stated-equal");
    let manifest = divergences(&[
        fabricated_item(0, REMOVED_CLASS, &with_reference(&quoted_list(&[KEPT_SET]))),
        fabricated_item(1, REMOVED_CLASS, ""),
    ]);
    let readings: [(Option<&'static str>, [&str; 2]); 3] = [
        (Some("present"), [KEPT_SET, LOST_SET]),
        (Some("absent"), [LOST_SET, LOST_SET]),
        (None, [LOST_SET, LOST_SET]),
    ];
    for (setting, actual) in readings {
        let options = match setting {
            Some(setting) => options(&["--reference-sequence", setting]),
            None => options(&[]),
        };
        for fmt in TEXT_FORMATS {
            let name = format!("{}-{fmt}", setting.unwrap_or("None"));
            let (status, report) = compare_fabricated(
                &dir.join(&name),
                fmt,
                &[KEPT_SET, KEPT_SET],
                &actual,
                &manifest,
                &options,
            );
            assert_eq!(status, 0, "{name}: {}", reasons(&report));
            assert!(
                is(member(&report, &["masks", "available"]), &removed(2)),
                "{name}"
            );
            assert!(is(member(&report, &["masks", "unused"]), "[]"), "{name}");
            let counted = if setting == Some("present") { "1" } else { "2" };
            assert!(
                is(member(applied(&report), &[REMOVED_CLASS]), counted),
                "{name}"
            );
        }
    }
}

/// The with-reference key is a non-empty sorted list of distinct consequence sets, other
/// than the item's own list, that names the reference's set alone or not at all; any other
/// shape is refused under either setting and under none.
#[test]
fn a_with_reference_key_of_another_shape_is_refused() {
    assert_eq!(WITH_REFERENCE_SETS, WITH_REFERENCE);
    let dir = scratch("term-pairs-with-reference-shape");
    let not_a_list = "is not a non-empty list";
    let not_a_set = "holds an item that is not a sorted, comma-joined set of terms";
    let cases = [
        (quoted(KEPT_SET), not_a_list),
        ("[]".to_string(), not_a_list),
        (
            format!("{{\"vep_rs_consequence_sets\": [{}]}}", quoted(KEPT_SET)),
            not_a_list,
        ),
        (format!("[{}, 1]", quoted(KEPT_SET)), not_a_set),
        (
            quoted_list(&["stop_retained_variant,3_prime_UTR_variant"]),
            not_a_set,
        ),
        (
            quoted_list(&[OTHER_SET, OTHER_SET]),
            "is not a sorted list of distinct sets",
        ),
        (
            quoted_list(&[LOST_SET]),
            "repeats vep_rs_consequence_sets and states no second outcome",
        ),
        (
            quoted_list(&[OTHER_SET, KEPT_SET]),
            "names Ensembl VEP's set beside another",
        ),
    ];
    for (case, (value, reason)) in cases.iter().enumerate() {
        let manifest = divergences(&[fabricated_item(0, REMOVED_CLASS, &with_reference(value))]);
        let settings: [&[&'static str]; 3] = [
            &["--reference-sequence", "present"],
            &["--reference-sequence", "absent"],
            &[],
        ];
        for (n, setting) in settings.into_iter().enumerate() {
            let (status, report, ran) = compare_fabricated_beside(
                &dir.join(format!("{case}-{n}")),
                "default",
                &[KEPT_SET],
                &[LOST_SET],
                Some(&manifest),
                &options(setting),
                RECORD_LINES,
            );
            assert_eq!(status, 2, "{value}");
            assert!(report.is_none(), "{value}");
            assert_eq!(
                ran.err(),
                format!(
                    "ERROR: [compare_vep_formats] the manifest's item on 21:100-101 - ENST1: {WITH_REFERENCE} {reason}\n"
                ),
                "{value}"
            );
        }
    }
}

/// An input no run of the manifest names reads the manifest's own lists, whose record
/// ordinals count the records of the corpus input. In the formats that name a record by its
/// ordinal the item is then looked up by no key: it is not counted, and its row counts once,
/// under the class's own predicate, as many entries as the default format applies.
fn an_item_of_an_input_no_run_names_is_not_looked_up_by_ordinal(fmt: &str) {
    let dir = scratch(&format!("term-pairs-input-no-run-names-{fmt}"));
    // The item's ordinal 1 is its record's in the corpus input; the compared input holds
    // that record first.
    let manifest = |records: &str| {
        format!(
            "{{\"divergences\": [{}], \"flag_sets\": [{{\"name\": \"twin\", \"input\": \"twin.txt\"}}]}}",
            fabricated_item_of(0, REMOVED_CLASS, records, "")
        )
    };
    let lines = [fabricated_record(0).0];
    let one = Beside {
        input_lines: Some(&lines),
        echoed: None,
        input_name: "guess.txt",
    };
    let compared = |name: &str, fmt: &str, records: &str, beside: Beside| {
        let (status, report, ran) = compare_fabricated_beside(
            &dir.join(name),
            fmt,
            &[KEPT_SET],
            &[LOST_SET],
            Some(&manifest(records)),
            &REFERENCE_MASKS,
            beside,
        );
        let report = report.unwrap_or_else(|| panic!("no report: {}", ran.err()));
        assert_eq!(status, 0, "{name}: {}", reasons(&report));
        report
    };
    let located = compared("default", "default", "[1]", one);
    assert!(is(applied(&located), &removed(1)));
    assert!(is(skipped(&located), "1"));
    let report = compared(fmt, fmt, "[1]", one);
    assert!(pyjson::py_eq(applied(&report), applied(&located)));
    assert!(is(skipped(&report), "0"));
    assert!(is(member(&report, &["masks", "unused"]), "[]"));
    // An input a flag set names is numbered by the lists that set reads, and keeps the
    // lookup.
    let named = Beside {
        input_name: "twin.txt",
        ..one
    };
    let report = compared("twin", fmt, "[0]", named);
    assert!(is(applied(&report), &removed(1)));
    assert!(is(skipped(&report), "1"));
}

#[test]
fn an_item_of_an_input_no_run_names_is_not_looked_up_by_ordinal_in_the_vcf_format() {
    an_item_of_an_input_no_run_names_is_not_looked_up_by_ordinal("vcf");
}

#[test]
fn an_item_of_an_input_no_run_names_is_not_looked_up_by_ordinal_in_the_json_format() {
    an_item_of_an_input_no_run_names_is_not_looked_up_by_ordinal("json");
}

/// The VCF format of an input that is not a VCF, and the JSON format where the objects
/// echo another text than the input's lines, name no Location and no input record, so no
/// manifest item meets a row: the item is not counted as applied and the row counts once,
/// under the class's own predicate. The class is then applied to as many entries as in the
/// default format, where the item meets its row, and it is used either way, so a pair that
/// names the class without a count stays green.
fn an_item_no_entry_can_meet_is_not_counted_and_its_row_counts_once(
    fmt: &str,
    echoed: Option<&str>,
) {
    let dir = scratch(&format!("term-pairs-item-no-entry-meets-{fmt}"));
    let manifest = divergences(&[fabricated_item(0, REMOVED_CLASS, "")]);
    let lines = [fabricated_ensembl_line(0, "\t")];
    let echo = echoed.map(|separator| [fabricated_ensembl_line(0, separator)]);
    let ensembl_lines = Beside {
        input_lines: Some(&lines),
        ..RECORD_LINES
    };
    let echoing = Beside {
        echoed: echo.as_ref().map(|echo| &echo[..]),
        ..ensembl_lines
    };
    let compared = |name: &str, fmt: &str, manifest: &str, masks: &[&str], beside: Beside| {
        let (status, report, ran) = compare_fabricated_beside(
            &dir.join(name),
            fmt,
            &[KEPT_SET],
            &[LOST_SET],
            Some(manifest),
            masks,
            beside,
        );
        let report = report.unwrap_or_else(|| panic!("no report: {}", ran.err()));
        (status, report)
    };
    let (status, located) = compared(
        "default",
        "default",
        &manifest,
        &REFERENCE_MASKS,
        ensembl_lines,
    );
    assert_eq!(status, 0, "{}", reasons(&located));
    assert!(is(applied(&located), &removed(1)));
    assert!(is(skipped(&located), "1"));
    let (status, report) = compared(fmt, fmt, &manifest, &REFERENCE_MASKS, echoing);
    assert_eq!(status, 0, "{}", reasons(&report));
    assert_eq!(member(&report, &["verdict"]).as_str(), Some("green"));
    assert!(is(skipped(&report), "0"));
    assert!(pyjson::py_eq(applied(&report), applied(&located)));
    assert!(is(member(&report, &["masks", "unused"]), "[]"));
    // A class with no predicate of its own is applied to no entry there, and the row
    // differs.
    let unmet = divergences(&[fabricated_item(0, "some_class", "")]);
    let (status, report) = compared("unmet", fmt, &unmet, &["--masks", "some_class"], echoing);
    assert_eq!(status, 1);
    assert!(is(applied(&report), "{}"));
    assert!(is(
        member(&report, &["masks", "unused"]),
        "[\"some_class\"]"
    ));
    // An input whose lines the entries name keeps the count as it is read.
    let (status, report) = compared("vcf-input", fmt, &manifest, &REFERENCE_MASKS, RECORD_LINES);
    assert_eq!(status, 0, "{}", reasons(&report));
    assert!(is(applied(&report), &removed(1)));
    assert!(is(skipped(&report), "1"));
}

#[test]
fn an_item_no_entry_can_meet_is_not_counted_and_its_row_counts_once_in_the_vcf_format() {
    an_item_no_entry_can_meet_is_not_counted_and_its_row_counts_once("vcf", None);
}

#[test]
fn an_item_no_entry_can_meet_is_not_counted_and_its_row_counts_once_in_the_json_format() {
    an_item_no_entry_can_meet_is_not_counted_and_its_row_counts_once("json", Some(" "));
}

/// The report of the default-format comparison of the reference's and the candidate's rows
/// under the class of the splice term past the record's span; `records` are the data lines
/// of the input VCF, which `--input` names where there is one.
fn splice_span_compared(directory: &Path, rows: &[(String, String)], records: &[String]) -> Value {
    rows_compared_under(SPLICE_SPAN_CLASS, directory, rows, records)
}

/// `splice_span_compared` under the class `class`.
fn rows_compared_under(
    class: &str,
    directory: &Path,
    rows: &[(String, String)],
    records: &[String],
) -> Value {
    std::fs::create_dir_all(directory).expect("the directory is made");
    let plant = |name: &str, text: String| {
        let path = directory.join(name);
        std::fs::write(&path, text).expect("a file is planted");
        path
    };
    let side = |pick: fn(&(String, String)) -> &String| {
        let rows: Vec<&String> = rows.iter().map(pick).collect();
        compare_common::default_text(&rows)
    };
    let reference = plant("reference.txt", side(|row| &row.0));
    let candidate = plant("candidate.txt", side(|row| &row.1));
    let written = directory.join("report.json");
    let mut argv = args_of(&[
        "--format",
        "default",
        "--expected",
        path_text(&reference),
        "--actual",
        path_text(&candidate),
        "--report",
        path_text(&written),
        "--masks",
        class,
    ]);
    if !records.is_empty() {
        let lines: String = records.iter().map(|record| format!("{record}\n")).collect();
        let input = plant(
            "input.vcf",
            format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{lines}"),
        );
        argv.extend(args_of(&["--input", path_text(&input)]));
    }
    let ran = run(directory, &argv);
    assert!(matches!(ran.code, 0 | 1), "{}", ran.err());
    report(&written)
}

fn verdict(report: &Value) -> &str {
    member(report, &["verdict"]).as_str().expect("a verdict")
}

fn spans(count: i64) -> String {
    format!("{{{}: {count}}}", quoted(SPLICE_SPAN_CLASS))
}

/// Under an identifier a default row spells no alleles, and under the class the comparison
/// reads them by position from the records of --input: two positions inside an exon, a span
/// across an exon's end, a span in an intron and two replaced bases are set aside, each
/// counted once with the IMPACT its terms decide, where a pure insertion at the Location of
/// the two replaced bases stays a difference. Without --input the same rows are read without
/// their records and stay differences, and three positions inside an exon are set aside
/// either way. A record spelt with a chr prefix is the row's; a record on another
/// chromosome, at another position or with another alternate allele and two records that
/// read one row differently state nothing. An input the comparison cannot read is its error
/// where a row asks for a record, and the input is read under no other class.
#[test]
fn splice_term_past_the_record_span_reads_a_row_under_an_identifier_from_the_input_records() {
    let dir = scratch("term-pairs-records-of-the-input");
    let mut records: Vec<SpliceSpanRecord> = splice_span_records()
        .into_iter()
        .map(|(_, record)| record)
        .collect();
    records.push(two_replaced());
    let named: Vec<(String, SpliceSpanRecord)> = records
        .into_iter()
        .enumerate()
        .map(|(number, record)| (format!("rec{}", number + 2), record))
        .collect();
    let rows: Vec<(String, String)> = named
        .iter()
        .map(|(name, record)| record.rows(Some(name)))
        .collect();
    let lines: Vec<String> = named
        .iter()
        .map(|(name, record)| record.input_line(name))
        .collect();
    let report = splice_span_compared(&dir.join("given"), &rows, &lines);
    assert_eq!(verdict(&report), "green");
    assert!(is(applied(&report), &spans(4)));
    assert!(is(
        member(&report, &["columns", "Consequence", "masked"]),
        &spans(4)
    ));
    assert!(is(
        member(&report, &["columns", "Extra", "masked"]),
        &spans(1)
    ));
    let report = splice_span_compared(&dir.join("without"), &rows, &[]);
    assert_eq!(verdict(&report), "red");
    assert!(is(applied(&report), "{}"));
    assert!(is(
        member(&report, &["columns", "Consequence", "differing"]),
        "4"
    ));
    assert!(is(member(&report, &["columns", "Extra", "differing"]), "1"));
    let plain = SpanColumns::default();
    let three_positions = [(
        splice_span_row(DONOR_PAST_THE_SPAN, &plain),
        splice_span_row(REGION_IN_ITS_PLACE, &plain),
    )];
    let its_record = ["7\t100997\trec1\tGACG\tGACGTA\t.\t.\t.".to_string()];
    for records in [&[][..], &its_record] {
        let report = splice_span_compared(
            &dir.join(format!("three_positions_{}", records.len())),
            &three_positions,
            records,
        );
        assert_eq!(verdict(&report), "green", "{records:?}");
        assert!(is(applied(&report), &spans(1)), "{records:?}");
    }
    let inserted = inserted_between_two();
    let report = splice_span_compared(
        &dir.join("inserted"),
        &[inserted.rows(Some("rec6"))],
        &[inserted.input_line("rec6")],
    );
    assert_eq!(verdict(&report), "red");
    assert!(is(applied(&report), "{}"));
    assert!(is(
        member(&report, &["columns", "Consequence", "differing"]),
        "1"
    ));
    let exon = named[0].1.clone();
    let line = |record: &SpliceSpanRecord, part: usize, value: &str| {
        let mut record = record.clone();
        record.line[part] = value.to_string();
        record
    };
    let spelt_with_a_prefix = line(&exon, 0, "chr7");
    let report = splice_span_compared(
        &dir.join("prefix"),
        &[exon.rows(Some("rec2"))],
        &[spelt_with_a_prefix.input_line("rec2")],
    );
    assert_eq!(verdict(&report), "green");
    assert!(is(applied(&report), &spans(1)));
    let stating_nothing = [
        (
            "another chromosome",
            exon.clone(),
            vec![line(&exon, 0, "8")],
        ),
        (
            "another position",
            exon.clone(),
            vec![line(&exon, 1, "100990")],
        ),
        (
            "another alternate allele",
            exon.clone(),
            vec![line(&exon, 3, "ACGTTC")],
        ),
        (
            "two records that read the row differently",
            two_replaced(),
            vec![two_replaced(), inserted.clone()],
        ),
    ];
    for (shape, record, in_the_input) in stating_nothing {
        let lines: Vec<String> = in_the_input
            .iter()
            .map(|other| other.input_line("rec2"))
            .collect();
        let report = splice_span_compared(
            &dir.join(shape.replace(' ', "_")),
            &[record.rows(Some("rec2"))],
            &lines,
        );
        assert_eq!(verdict(&report), "red", "{shape}");
        assert!(is(applied(&report), "{}"), "{shape}");
    }
    // A Location that names no span of whole positions asks for no record.
    let mut no_span = exon.clone();
    no_span.location = "7:10099\u{b2}-101000".to_string();
    let report = splice_span_compared(
        &dir.join("no_span"),
        &[no_span.rows(Some("rec2"))],
        &[exon.input_line("rec2")],
    );
    assert_eq!(verdict(&report), "red");
    assert!(is(applied(&report), "{}"));
    let given = dir.join("given");
    let exit_code = |actual: &str, masks: &str| {
        run(
            &dir,
            &args_of(&[
                "--format",
                "default",
                "--expected",
                path_text(&given.join("reference.txt")),
                "--input",
                path_text(&dir.join("none.vcf")),
                "--actual",
                path_text(&given.join(actual)),
                "--masks",
                masks,
            ]),
        )
        .code
    };
    assert_eq!(exit_code("candidate.txt", SPLICE_SPAN_CLASS), 2);
    assert_eq!(exit_code("reference.txt", SPLICE_SPAN_CLASS), 0);
    assert_eq!(exit_code("candidate.txt", "order_only"), 1);
}

const INSERTION_RANGE_CLASS: &str = "stop_codon_overlap_insertion_range";

/// Under the insertion range class a row at two neighbouring positions is an insertion's
/// where its record's alleles place the allele between them. A first column that spells a
/// padded insertion places it with and without --input. Under an identifier the records of
/// --input place an insertion of two bases, which stays a difference without them, and show
/// two bases replaced by three, which are set aside without them alone. An insertion before
/// its record's first base is placed by no record: its row is read by its allele, and stays a
/// difference beside a record that replaces the same two bases by the same allele.
#[test]
fn stop_codon_overlap_insertion_range_reads_an_insertion_from_the_records_alleles() {
    let dir = scratch("term-pairs-insertion-range-records-of-the-input");
    let pair = |name: &str, allele: &str| {
        let row = |term: &str| {
            let terms = format!("inframe_insertion,{term}");
            stop_codon_row_at(
                allele,
                &terms,
                "7054-7055",
                "6977-6978",
                "19:13208046-13208047",
            )
            .replacen("1993671", name, 1)
        };
        [(row("coding_sequence_variant"), row("stop_retained_variant"))]
    };
    let line = |position: &str, name: &str, reference_allele: &str, alternate: &str| {
        format!("19\t{position}\t{name}\t{reference_allele}\t{alternate}\t.\t.\t.")
    };
    let set_aside = |tag: &str, rows: &[(String, String)], records: &[String]| {
        let report = rows_compared_under(INSERTION_RANGE_CLASS, &dir.join(tag), rows, records);
        let once = format!("{{{}: 1}}", quoted(INSERTION_RANGE_CLASS));
        match verdict(&report) {
            "green" => {
                assert!(is(applied(&report), &once), "{tag}");
                true
            }
            "red" => {
                assert!(is(applied(&report), "{}"), "{tag}");
                false
            }
            other => panic!("{tag}: the verdict {other}"),
        }
    };

    let padded = pair("19_13208045_TA/TATGA", "TGA");
    assert!(set_aside("padded", &padded, &[]));
    let its_record = [line("13208045", ".", "TA", "TATGA")];
    assert!(set_aside("padded_record", &padded, &its_record));

    let two_bases = pair("1993671", "CG");
    assert!(!set_aside("two_bases", &two_bases, &[]));
    for (tag, position, reference_allele, alternate) in [
        ("two_bases_record", "13208046", "A", "ACG"),
        ("two_bases_padded_record", "13208045", "TA", "TACG"),
    ] {
        let its_record = [line(position, "1993671", reference_allele, alternate)];
        assert!(set_aside(tag, &two_bases, &its_record), "{tag}");
    }

    let three_bases = pair("1993671", "CGT");
    assert!(set_aside("replaced", &three_bases, &[]));
    let replacing = [line("13208046", "1993671", "AC", "CGT")];
    assert!(!set_aside("replaced_record", &three_bases, &replacing));

    let unplaced = line("13208047", "1993671", "A", "CGTA");
    assert!(set_aside(
        "unplaced_record",
        &three_bases,
        std::slice::from_ref(&unplaced)
    ));
    let beside = [line("13208046", "1993672", "AC", "CGT"), unplaced];
    assert!(!set_aside("unplaced_beside", &three_bases, &beside));
}

/// The comparison of the default-format outputs of `reference` and `candidate` under
/// `options`: the report where one was written and the run.
fn compare_rows(
    dir: &Path,
    name: &str,
    reference: &[String],
    candidate: &[String],
    options: &[&str],
) -> (Option<Value>, Ran) {
    let plant = |file: String, rows: &[String]| {
        let path = dir.join(file);
        std::fs::write(&path, compare_common::default_text(rows)).expect("a file is planted");
        path
    };
    let expected = plant(format!("expected.{name}.txt"), reference);
    let actual = plant(format!("actual.{name}.txt"), candidate);
    let written = dir.join(format!("whole.{name}.json"));
    let mut argv = args_of(&[
        "--format",
        "default",
        "--expected",
        path_text(&expected),
        "--actual",
        path_text(&actual),
        "--report",
        path_text(&written),
    ]);
    argv.extend(args_of(options));
    let ran = run(dir, &argv);
    let report = written.is_file().then(|| report(&written));
    (report, ran)
}

/// A record's label is its Location, a space and its Uploaded_variation, so a Location
/// holding a space can share a label with a row of another Location: the candidate-only row
/// of `frequency_filter_zero_value` and of `merged_same_stable_id_collapse` is judged by
/// the labels of the reference's rows, a row at another Location among them. The first
/// class is applied under a reference release before 116.3, which corrects its defect.
#[test]
fn a_candidate_only_row_is_judged_by_the_labels_of_the_reference_where_a_location_holds_a_space() {
    let dir = scratch("term-pairs-location-with-a-space");
    let row = compare_common::default_row;
    let reference = [row(
        "-",
        "21:5 x",
        "A",
        "ENST1",
        "intron_variant",
        "IMPACT=MODIFIER",
    )];
    let beside = |extra: &str| {
        vec![
            reference[0].clone(),
            row("x", "21:5", "A", "ENST2", "intron_variant", extra),
        ]
    };
    let (_, ran) = compare_rows(
        &dir,
        "zero",
        &reference,
        &beside("IMPACT=MODIFIER;FREQS=pop:A:0"),
        &[
            "--masks",
            "frequency_filter_zero_value",
            "--reference-release",
            "116.2",
        ],
    );
    assert_eq!(
        (ran.code, ran.out().as_str(), ran.err().as_str()),
        (
            1,
            "default: red\n  unexpected in the candidate: 1 entries\n",
            ""
        )
    );
    let (_, ran) = compare_rows(
        &dir,
        "collapsed",
        &reference,
        &beside("IMPACT=MODIFIER;SOURCE=RefSeq"),
        &["--masks", "merged_same_stable_id_collapse"],
    );
    assert_eq!(
        (ran.out().as_str(), ran.code, ran.err().as_str()),
        ("default: green\n", 0, "")
    );
}

/// One transcript row of record `record`, at a Location of the record's own.
fn planted_row(record: usize, feature: &str, extra: Option<&str>) -> String {
    let own = format!("IMPACT=MODIFIER;STRAND=1;DOMAINS=db:a{record},db:b{record}");
    compare_common::default_row(
        &format!("rec{record}"),
        &format!("21:{}", 1000 + record),
        &"ACGT"[record % 4..record % 4 + 1],
        feature,
        "intron_variant",
        extra.unwrap_or(&own),
    )
}

/// A term-pair class reads the SHIFT_LENGTH of every Transcript row both sides write. A
/// superscript digit, which is a digit and no number, and a number of more digits than the
/// script converts end the comparison with the comparator's error; a SHIFT_LENGTH of 640
/// digits is a number and its row is compared.
#[test]
fn a_shift_length_that_is_no_number_ends_the_comparison_with_the_comparators_error() {
    let dir = scratch("term-pairs-shift-length-no-number");
    let rows: Vec<String> = (0..8)
        .flat_map(|record| ["ENST1", "ENST2"].map(|feature| planted_row(record, feature, None)))
        .collect();
    let shifted = |value: &str| {
        let mut shifted = rows.clone();
        shifted[0] = planted_row(
            0,
            "ENST1",
            Some(&format!("IMPACT=MODIFIER;SHIFT_LENGTH={value}")),
        );
        shifted
    };
    let long = "7".repeat(5000);
    for (name, value, klass, shown) in [
        (
            "superscript",
            "\u{b2}",
            "mature_mirna_unshifted_span",
            "\u{b2}".to_string(),
        ),
        (
            "long",
            long.as_str(),
            "stop_codon_unshifted_span",
            format!("{}...", "7".repeat(24)),
        ),
    ] {
        let rows = shifted(value);
        let (report, ran) = compare_rows(&dir, name, &rows, &rows, &["--masks", klass]);
        assert!(report.is_none(), "{name}");
        assert_eq!(
            (ran.out().as_str(), ran.code),
            ("", 2),
            "{name}: {}",
            ran.err()
        );
        assert_eq!(
            ran.err(),
            format!("ERROR: [compare_vep_formats] SHIFT_LENGTH '{shown}' is not an integer\n"),
            "{name}"
        );
    }
    let rows = shifted(&"7".repeat(640));
    let (_, ran) = compare_rows(
        &dir,
        "read",
        &rows,
        &rows,
        &["--masks", "stop_codon_unshifted_span"],
    );
    assert_eq!(
        (ran.out().as_str(), ran.code, ran.err().as_str()),
        ("default: green\n", 0, "")
    );
}
