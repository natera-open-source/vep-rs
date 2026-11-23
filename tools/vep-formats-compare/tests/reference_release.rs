// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The reference release of a comparison, through the binary on planted files: a class
//! whose defect Ensembl VEP corrects is applied under an earlier release and is not applied
//! under the release that corrects it, where the rows and fields it set aside are
//! differences and the report names the class under `masks.corrected`.
//!
//! `tests/data/reference_release/<case>.json` is the report of a case's command under a
//! release before the one that corrects the case's class.

mod compare_common;
mod support;

use std::ffi::OsString;
use std::path::Path;

use compare_common::default_row;
use compare_common::default_text;
use support::args_of;
use support::run;
use support::run_in;
use support::scratch;
use support::Ran;
use vep_formats_compare::classes;
use vep_formats_compare::classes::DEFAULT_REFERENCE_RELEASE;
use vep_formats_compare::classes::REFERENCE_RELEASES;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::remainder::WHOLE_VARIABLE;
use vep_formats_compare::report::report_text;

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/reference_release");
const USAGE: &str = include_str!("../data/usage.txt");
const ZERO_VALUE: &str = classes::FREQUENCY_FILTER_ZERO_VALUE;
const KEY_ORDER: &str = classes::FREQUENCY_FILTER_ALLELE_KEY_ORDER;

/// The release that corrects the defects of both classes, and the releases before it.
const CORRECTING: &str = "116.3";
const EARLIER: [&str; 2] = ["116.2", "115.2"];

/// The two records of the JSON cases, as the lines of `variants.vcf`.
const RECORDS: [&str; 2] = [
    "21\t5131441\t.\tT\tC\t.\t.\t.",
    "21\t29025626\t.\tA\tG\t.\t.\t.",
];

/// One comparison of planted files under one class Ensembl VEP 116.3 corrects.
struct Case {
    /// The name of the case and of its stored report.
    name: &'static str,
    class: &'static str,
    format: &'static str,
    expected: &'static str,
    actual: &'static str,
    /// The files the class reads beside the two outputs.
    options: &'static [&'static str],
    /// What the class sets aside under a release before 116.3.
    applied: i64,
    /// The manifest's items of the class.
    available: i64,
    /// The record entries counted on a record the manifest states the reference dropped.
    on_skipped_records: i64,
    /// The reasons of the red verdict where the class is not applied, and the member of the
    /// report that counts the rows or fields the class set aside as differences.
    reasons: &'static [&'static str],
    difference: (&'static [&'static str], i64),
}

const CASES: [Case; 5] = [
    // A FREQS the reference reads NA and the candidate 0.
    Case {
        name: "zero_value_field",
        class: ZERO_VALUE,
        format: "default",
        expected: "field.expected.txt",
        actual: "field.actual.txt",
        options: &[],
        applied: 1,
        available: 0,
        on_skipped_records: 0,
        reasons: &[
            "field Extra: 1 differing entries",
            "field FREQS: 1 differing entries",
        ],
        difference: (&["columns", "FREQS", "differing"], 1),
    },
    // The rows of a record the candidate alone keeps, each on a FREQS of zeros.
    Case {
        name: "zero_value_rows",
        class: ZERO_VALUE,
        format: "default",
        expected: "rows.expected.txt",
        actual: "rows.actual.txt",
        options: &[],
        applied: 2,
        available: 0,
        on_skipped_records: 0,
        reasons: &["unexpected in the candidate: 2 entries"],
        difference: (&["one_sided", "actual_only"], 2),
    },
    // A JSON record the candidate alone keeps, which a manifest documents under the class.
    Case {
        name: "zero_value_manifest",
        class: ZERO_VALUE,
        format: "json",
        expected: "kept.expected.jsonl",
        actual: "kept.actual.jsonl",
        options: &[
            "--input",
            "variants.vcf",
            "--manifest",
            "zero_value.manifest.json",
        ],
        applied: 1,
        available: 1,
        on_skipped_records: 1,
        reasons: &["unexpected in the candidate: 2 entries"],
        difference: (&["one_sided", "actual_only"], 2),
    },
    Case {
        name: "key_order_manifest",
        class: KEY_ORDER,
        format: "json",
        expected: "kept.expected.jsonl",
        actual: "kept.actual.jsonl",
        options: &[
            "--input",
            "variants.vcf",
            "--manifest",
            "key_order.manifest.json",
        ],
        applied: 1,
        available: 1,
        on_skipped_records: 1,
        reasons: &["unexpected in the candidate: 2 entries"],
        difference: (&["one_sided", "actual_only"], 2),
    },
    // The default-format rows of a record the candidate alone keeps, each a tuple of a
    // manifest under the class.
    Case {
        name: "key_order_rows",
        class: KEY_ORDER,
        format: "default",
        expected: "rows.expected.txt",
        actual: "rows.actual.txt",
        options: &["--manifest", "key_order_rows.manifest.json"],
        applied: 2,
        available: 2,
        on_skipped_records: 0,
        reasons: &["unexpected in the candidate: 2 entries"],
        difference: (&["one_sided", "actual_only"], 2),
    },
];

fn quoted(text: &str) -> String {
    pyjson::dumps_str(text, true)
}

fn freqs_row(name: &str, location: &str, allele: &str, feature: &str, frequency: &str) -> String {
    default_row(
        name,
        location,
        allele,
        feature,
        "intron_variant",
        &format!("IMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE_AFR:{allele}:{frequency}"),
    )
}

/// The JSON object of the record `RECORDS[record]`: one transcript consequence, and a
/// frequency of 0 on its allele.
fn json_object(record: usize, allele: &str) -> String {
    format!(
        "{{\"input\": {}, \"id\": \"r{record}\", \"most_severe_consequence\": \"intron_variant\", \"colocated_variants\": [{{\"frequencies\": {{{}: {{\"gnomade_afr\": 0}}}}}}], \"transcript_consequences\": [{{\"variant_allele\": {}, \"transcript_id\": \"ENST{record}\", \"consequence_terms\": [\"intron_variant\"]}}]}}\n",
        quoted(RECORDS[record]),
        quoted(allele),
        quoted(allele)
    )
}

/// A manifest whose candidate-only tuples are `tuples`, each a feature of the second
/// record with its record ordinals as a JSON list, under `class`.
fn manifest(class: &str, tuples: &[(&str, &str)]) -> String {
    let items: Vec<String> = tuples
        .iter()
        .map(|(feature, records)| {
            format!(
                "{{\"location\": \"21:29025626\", \"allele\": \"G\", \"feature\": {}, \"record_indices\": {records}, \"expected_divergence\": {}}}",
                quoted(feature),
                quoted(class)
            )
        })
        .collect();
    format!("{{\"vep_rs_only_tuples\": [{}]}}\n", items.join(", "))
}

/// A directory of the test `name` with the files of every case.
fn planted(name: &str) -> std::path::PathBuf {
    let dir = scratch(&format!("reference-release-{name}"));
    let plant = |file: &str, text: String| {
        std::fs::write(dir.join(file), text).expect("a file is planted");
    };
    let read_as_na = default_text(&[freqs_row("a", "21:5131441", "C", "ENST1", "NA")]);
    plant(
        "field.actual.txt",
        read_as_na.replace("GNOMADE_AFR:C:NA", "GNOMADE_AFR:C:0"),
    );
    plant("field.expected.txt", read_as_na);
    let kept_by_both = freqs_row("a", "21:5131441", "C", "ENST1", "0");
    plant("rows.expected.txt", default_text(&[kept_by_both.as_str()]));
    plant(
        "rows.actual.txt",
        default_text(&[
            kept_by_both,
            freqs_row("b", "21:29025626", "G", "ENST2", "0"),
            freqs_row("b", "21:29025626", "G", "ENST3", "0"),
        ]),
    );
    plant("kept.expected.jsonl", json_object(0, "C"));
    plant(
        "kept.actual.jsonl",
        format!("{}{}", json_object(0, "C"), json_object(1, "G")),
    );
    plant(
        "variants.vcf",
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}\n{}\n",
            RECORDS[0], RECORDS[1]
        ),
    );
    plant(
        "zero_value.manifest.json",
        manifest(ZERO_VALUE, &[("ENST1", "[1]")]),
    );
    plant(
        "key_order.manifest.json",
        manifest(KEY_ORDER, &[("ENST1", "[1]")]),
    );
    plant(
        "key_order_rows.manifest.json",
        manifest(KEY_ORDER, &[("ENST2", "[]"), ("ENST3", "[]")]),
    );
    plant(
        "warnings.txt",
        "WARNING: line 3 skipped (21 5 x): invalid\n".to_string(),
    );
    dir
}

/// One run of the binary and the report it wrote, as its bytes.
struct Run {
    ran: Ran,
    report: Option<Vec<u8>>,
}

impl Run {
    /// The exit code and the two standard streams.
    fn ended(&self) -> (i32, String, String) {
        (self.ran.code, self.ran.out(), self.ran.err())
    }

    fn bytes(&self) -> &[u8] {
        self.report.as_deref().expect("the run wrote a report")
    }

    fn value(&self) -> Value {
        let text = std::str::from_utf8(self.bytes()).expect("a report is text");
        pyjson::loads(text).expect("a report is JSON")
    }
}

/// The comparison of `case` in `dir` under `masks` and `more` options, with `env` set.
fn compared_in(
    dir: &Path,
    case: &Case,
    masks: Option<&str>,
    more: &[&str],
    env: &[(OsString, OsString)],
) -> Run {
    let report = dir.join("report.json");
    let _ = std::fs::remove_file(&report);
    let mut words = vec![
        "--format",
        case.format,
        "--expected",
        case.expected,
        "--actual",
        case.actual,
        "--report",
        "report.json",
    ];
    words.extend(case.options);
    if let Some(masks) = masks {
        words.extend(["--masks", masks]);
    }
    words.extend(more);
    let ran = run_in(dir, env, &args_of(&words));
    Run {
        ran,
        report: std::fs::read(&report).ok(),
    }
}

fn compared(dir: &Path, case: &Case, masks: Option<&str>, more: &[&str]) -> Run {
    compared_in(dir, case, masks, more, &[])
}

/// The member of `report` at `path`.
fn at<'a>(report: &'a Value, path: &[&str]) -> &'a Value {
    support::member(report, path)
}

/// Whether `value` is the value the JSON text `expected` spells.
fn is(value: &Value, expected: &str) -> bool {
    pyjson::py_eq(value, &pyjson::loads(expected).expect("the text is JSON"))
}

fn count_of(class: &str, count: i64) -> String {
    if count == 0 {
        "{}".to_string()
    } else {
        format!("{{{}: {count}}}", quoted(class))
    }
}

/// `report` without the member `masks.corrected`, as the text of a report file.
fn without_corrected(report: &Value) -> String {
    let mut members = report.as_dict().expect("a report is an object").clone();
    let mut masks = members["masks"]
        .as_dict()
        .expect("masks is an object")
        .clone();
    assert!(masks.shift_remove("corrected").is_some());
    members.insert("masks".to_string(), Value::Dict(masks));
    report_text(&Value::Dict(members))
}

/// Under a release before the one that corrects it, a class sets aside what its rule and
/// the manifest give it: the verdict is green, the class is requested and applied, and the
/// report is the stored one to the byte.
#[test]
fn a_corrected_class_is_applied_under_an_earlier_reference_release() {
    let dir = planted("earlier");
    for case in &CASES {
        let stored = std::fs::read(format!("{DATA}/{}.json", case.name))
            .unwrap_or_else(|error| panic!("{}: {error}", case.name));
        for release in EARLIER {
            let run = compared(
                &dir,
                case,
                Some(case.class),
                &["--reference-release", release],
            );
            assert_eq!(
                run.ended(),
                (0, format!("{}: green\n", case.format), String::new()),
                "{} under {release}",
                case.name
            );
            let report = run.value();
            let masks = at(&report, &["masks"]);
            assert!(is(
                at(masks, &["requested"]),
                &format!("[{}]", quoted(case.class))
            ));
            assert!(
                is(at(masks, &["applied"]), &count_of(case.class, case.applied)),
                "{} under {release}",
                case.name
            );
            assert!(is(at(masks, &["unused"]), "[]"));
            assert!(is(
                at(masks, &["available"]),
                &count_of(case.class, case.available)
            ));
            assert!(masks.get("corrected").is_none());
            assert!(is(
                at(
                    &report,
                    &["counts", "record_entries_on_reference_skipped_records"]
                ),
                &case.on_skipped_records.to_string()
            ));
            assert!(
                run.bytes() == stored,
                "{} under {release} differs from its stored report",
                case.name
            );
        }
    }
}

/// Under the release that corrects it, named and by default, the class is applied nowhere:
/// the rows and fields it set aside are differences, the verdict is red, a manifest's items
/// of the class are counted under `masks.available` as those of a class the command line
/// does not name, and the class stands under `masks.corrected` and under none of
/// `requested`, `applied` and `unused`. That member aside, the run is the run of the same
/// command line without the class.
#[test]
fn a_corrected_class_is_not_applied_under_the_release_that_corrects_it() {
    assert_eq!(DEFAULT_REFERENCE_RELEASE, CORRECTING);
    let dir = planted("correcting");
    for case in &CASES {
        let unnamed = compared(&dir, case, None, &[]);
        for release in [&["--reference-release", CORRECTING][..], &[]] {
            let run = compared(&dir, case, Some(case.class), release);
            let reasons: String = case
                .reasons
                .iter()
                .map(|reason| format!("  {reason}\n"))
                .collect();
            let ended = (1, format!("{}: red\n{reasons}", case.format), String::new());
            assert_eq!(run.ended(), ended, "{} under {release:?}", case.name);
            let report = run.value();
            assert!(is(at(&report, &["verdict"]), "\"red\""));
            let (member, count) = case.difference;
            assert!(
                is(at(&report, member), &count.to_string()),
                "{} under {release:?}",
                case.name
            );
            assert!(is(
                at(
                    &report,
                    &["counts", "record_entries_on_reference_skipped_records"]
                ),
                "0"
            ));
            let masks = at(&report, &["masks"]);
            assert!(is(at(masks, &["requested"]), "[]"));
            assert!(is(at(masks, &["applied"]), "{}"));
            assert!(is(at(masks, &["unused"]), "[]"));
            assert!(is(
                at(masks, &["available"]),
                &count_of(case.class, case.available)
            ));
            assert!(is(
                at(masks, &["corrected"]),
                &format!("{{{}: {}}}", quoted(case.class), quoted(CORRECTING))
            ));
            let names: Vec<&str> = masks
                .as_dict()
                .expect("masks is an object")
                .keys()
                .map(String::as_str)
                .collect();
            assert_eq!(
                names,
                ["requested", "applied", "unused", "available", "corrected"]
            );
            assert_eq!(unnamed.ended(), ended, "{}", case.name);
            assert!(
                without_corrected(&report).as_bytes() == unnamed.bytes(),
                "{} under {release:?} differs from the run that does not name the class",
                case.name
            );
        }
    }
}

/// A report has no `corrected` member where the command line names no class the reference
/// release corrects: without `--masks` and under another class the report is the same
/// under every release and by default. A manifest's items of the unnamed class are counted
/// as available, as those of any class the command line does not name.
#[test]
fn a_report_has_no_corrected_member_where_no_corrected_class_is_named() {
    let dir = planted("unnamed");
    for case in &CASES {
        let mut reports: Vec<Vec<u8>> = Vec::new();
        for masks in [None, Some("order_only")] {
            for release in REFERENCE_RELEASES {
                let run = compared(&dir, case, masks, &["--reference-release", release]);
                assert_eq!(run.ran.code, 1, "{}", case.name);
                let report = run.value();
                let held = at(&report, &["masks"]);
                assert!(held.get("corrected").is_none(), "{}", case.name);
                assert!(is(
                    at(held, &["available"]),
                    &count_of(case.class, case.available)
                ));
                reports.push(run.bytes().to_vec());
            }
            let run = compared(&dir, case, masks, &[]);
            reports.push(run.bytes().to_vec());
        }
        let (unmasked, order_only) = reports.split_at(REFERENCE_RELEASES.len() + 1);
        assert!(unmasked.iter().all(|report| report == &unmasked[0]));
        assert!(order_only.iter().all(|report| report == &order_only[0]));
        assert_ne!(unmasked[0], order_only[0]);
    }
}

/// A class the reference release corrects leaves the named set alone: the classes named
/// beside it are requested, applied and unused as without it, and every corrected class
/// named stands under `masks.corrected`, in the order of their names.
#[test]
fn the_classes_named_beside_a_corrected_class_are_applied_as_without_it() {
    let dir = planted("beside");
    let case = &CASES[0];
    let named = format!("order_only, {ZERO_VALUE} ,{KEY_ORDER},");
    let run = compared(&dir, case, Some(&named), &[]);
    let report = run.value();
    let masks = at(&report, &["masks"]);
    assert!(is(at(masks, &["requested"]), "[\"order_only\"]"));
    assert!(is(at(masks, &["unused"]), "[\"order_only\"]"));
    assert!(is(
        at(masks, &["corrected"]),
        &format!(
            "{{{}: \"116.3\", {}: \"116.3\"}}",
            quoted(KEY_ORDER),
            quoted(ZERO_VALUE)
        )
    ));
    let corrected: Vec<&str> = at(masks, &["corrected"])
        .as_dict()
        .expect("corrected is an object")
        .keys()
        .map(String::as_str)
        .collect();
    assert_eq!(corrected, [KEY_ORDER, ZERO_VALUE]);
    let beside = compared(&dir, case, Some("order_only"), &[]);
    assert_eq!(run.ended(), beside.ended());
    assert!(without_corrected(&report).as_bytes() == beside.bytes());
    let earlier = compared(&dir, case, Some(&named), &["--reference-release", "116.2"]);
    let report = earlier.value();
    let masks = at(&report, &["masks"]);
    assert!(is(
        at(masks, &["requested"]),
        &format!(
            "[{}, {}, \"order_only\"]",
            quoted(KEY_ORDER),
            quoted(ZERO_VALUE)
        )
    ));
    assert!(masks.get("corrected").is_none());
    assert_eq!(earlier.ran.code, 0);
}

/// The stream of record groups and the comparison of the whole outputs write the same
/// report, standard output and exit code under each release.
#[test]
fn both_comparison_paths_agree_under_each_reference_release() {
    let dir = planted("both-paths");
    let whole = [(OsString::from(WHOLE_VARIABLE), OsString::from("1"))];
    for case in &CASES {
        for release in REFERENCE_RELEASES {
            let options = ["--reference-release", release];
            let streamed = compared(&dir, case, Some(case.class), &options);
            let compared_whole = compared_in(&dir, case, Some(case.class), &options, &whole);
            assert_eq!(
                streamed.ended(),
                compared_whole.ended(),
                "{} under {release}",
                case.name
            );
            assert!(
                streamed.bytes() == compared_whole.bytes(),
                "{} under {release}",
                case.name
            );
            let corrects = *release == CORRECTING;
            assert_eq!(streamed.ran.code, i32::from(corrects));
            assert_eq!(
                at(&streamed.value(), &["masks"]).get("corrected").is_some(),
                corrects
            );
        }
    }
}

/// The option is read in the three spellings every option takes, and a later occurrence
/// replaces an earlier one.
#[test]
fn the_reference_release_is_read_in_every_spelling_of_an_option() {
    let dir = planted("spellings");
    let case = &CASES[1];
    let spelt = compared(
        &dir,
        case,
        Some(case.class),
        &["--reference-release", "116.2"],
    );
    assert_eq!(spelt.ran.code, 0);
    for spelling in [
        &["--reference-release=116.2"][..],
        &["--reference-r", "116.2"],
        &["--reference-rel=116.2"],
        &[
            "--reference-release",
            "116.3",
            "--reference-release",
            "116.2",
        ],
    ] {
        let run = compared(&dir, case, Some(case.class), spelling);
        assert_eq!(run.ended(), spelt.ended(), "{spelling:?}");
        assert!(run.bytes() == spelt.bytes(), "{spelling:?}");
    }
    let last = compared(
        &dir,
        case,
        Some(case.class),
        &["--reference-release", "116.2", "--reference-release=116.3"],
    );
    assert_eq!(last.ran.code, 1);
}

/// A value that is no release of the option's choices is refused as any invalid choice is:
/// the usage block and the parser's message on the standard error, exit code 2 and no
/// report, in every mode.
#[test]
fn a_reference_release_outside_the_choices_is_refused() {
    let dir = planted("refused");
    let case = &CASES[0];
    let refused = |value: &str| {
        format!(
            "{USAGE}vep-formats-compare: error: argument --reference-release: invalid choice: '{value}' (choose from '115.2', '116.2', '116.3')\n"
        )
    };
    for value in ["116.4", "116", "116.3.0", "116.30", "release/116.3", ""] {
        for spelling in [
            vec!["--reference-release", value],
            vec!["--reference-r", value],
        ] {
            let run = compared(&dir, case, Some(case.class), &spelling);
            assert_eq!(
                run.ended(),
                (2, String::new(), refused(value)),
                "{spelling:?}"
            );
            assert!(run.report.is_none(), "{spelling:?}");
        }
    }
    let equals = compared(&dir, case, None, &["--reference-release=117.0"]);
    assert_eq!(equals.ended(), (2, String::new(), refused("117.0")));
    let lacking = compared(&dir, case, None, &["--reference-release"]);
    assert_eq!(
        lacking.ended(),
        (
            2,
            String::new(),
            format!(
                "{USAGE}vep-formats-compare: error: argument --reference-release: expected one argument\n"
            )
        )
    );
    for mode in ["--warnings", "--self-parity"] {
        let ran = run(&dir, &args_of(&[mode, "--reference-release", "116.4"]));
        assert_eq!(
            (ran.code, ran.out(), ran.err()),
            (2, String::new(), refused("116.4")),
            "{mode}"
        );
    }
}

/// The warnings mode takes the option beside `--masks`: neither class reads a warnings
/// file, so the report, which names the classes applied alone, is the same under every
/// release, with the class named and without.
#[test]
fn the_warnings_mode_takes_the_reference_release() {
    let dir = planted("warnings");
    let warnings = |more: &[&str]| -> (Ran, Vec<u8>) {
        let report = dir.join("warnings.json");
        let _ = std::fs::remove_file(&report);
        let mut words = vec![
            "--warnings",
            "--expected",
            "warnings.txt",
            "--actual",
            "warnings.txt",
            "--report",
            "warnings.json",
        ];
        words.extend(more);
        let ran = run(&dir, &args_of(&words));
        (ran, std::fs::read(&report).expect("the run wrote a report"))
    };
    let (plain, plain_report) = warnings(&[]);
    assert_eq!(
        (plain.code, plain.out(), plain.err()),
        (0, "warnings: green\n".to_string(), String::new())
    );
    let named = format!("{ZERO_VALUE},{KEY_ORDER}");
    for release in REFERENCE_RELEASES {
        for masks in [&[][..], &["--masks", named.as_str()]] {
            let mut more = vec!["--reference-release", release];
            more.extend(masks);
            let (ran, report) = warnings(&more);
            assert_eq!(
                (ran.code, ran.out(), ran.err()),
                (plain.code, plain.out(), plain.err()),
                "{more:?}"
            );
            assert!(report == plain_report, "{more:?}");
        }
    }
}
