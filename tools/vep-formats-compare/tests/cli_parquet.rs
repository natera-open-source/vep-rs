// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The Parquet cases of the command line, each run through the binary and held to its
//! report, standard output and exit code: the reference is a tab output and the candidate
//! its rows planted as a Parquet file of the writer's layout.

mod parquet_writer;
mod support;

use parquet_writer::plant_tab_rows;
use std::path::Path;
use support::{args_of, compact, member, report, run, run_in, scratch, Ran};
use vep_formats_compare::pyjson::{self, Value};
use vep_formats_compare::remainder::WHOLE_VARIABLE;

const KEPT_SET: &str = "3_prime_UTR_variant,stop_retained_variant";
const LOST_SET: &str = "3_prime_UTR_variant,stop_lost";
const OTHER_SET: &str = "3_prime_UTR_variant,coding_sequence_variant";
const WITH_REFERENCE: &str = "vep_rs_consequence_sets_with_reference_sequence";
const REMOVED_CLASS: &str = "stop_codon_removed_without_stop_codon_test";
const REFERENCE_MASKS: [&str; 2] = ["--masks", REMOVED_CLASS];
const RUNS_MASKS: [&str; 2] = ["--masks", "some_class"];
/// The line a comparison beside `runs_manifest` ends on where it names `some_class` and no
/// run: three of the manifest's runs annotate the compared input.
const UNNAMED_RUNS: &str = "ERROR: [compare_vep_formats] --run is not given, and the manifest \
    lists 3 runs over variants.vcf (default, other, plain) that document some_class, which \
    --masks names: name the pair's run\n";
const NOT_OBSERVED: &str = "[\"documented divergences not observed: 1\"]";
const CONSEQUENCE_DIFFERS: &str = "field Consequence: 1 differing entries";
const COLUMNS: &str =
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIMPACT";

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

fn quoted(text: &str) -> String {
    pyjson::dumps_str(text, true)
}

fn path_text(path: &Path) -> &str {
    path.to_str().expect("the path is UTF-8")
}

/// The VCF line, the Location and the transcript of the fabricated record `ordinal`: one
/// deletion on one transcript.
fn fabricated_record(ordinal: usize) -> (String, String, String) {
    let position = 100 * (ordinal + 1);
    (
        format!("21\t{}\tdel{ordinal}\tCAT\tC\t.\t.\t.", position - 1),
        format!("21:{position}-{}", position + 1),
        format!("ENST{}", ordinal + 1),
    )
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

fn divergences(items: &[String]) -> String {
    format!("{{\"divergences\": [{}]}}", items.join(", "))
}

/// The member that states an item's sets of a run with a reference sequence.
fn with_reference(set: &str) -> String {
    format!(", {}: [{}]", quoted(WITH_REFERENCE), quoted(set))
}

/// A manifest whose one item carries the reference's own set as its with-reference sets.
fn reference_manifest() -> String {
    divergences(&[fabricated_item(0, REMOVED_CLASS, &with_reference(KEPT_SET))])
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

/// One tab row per fabricated record, record `i` carrying `sets[i]`.
fn fabricated_rows(sets: &[&str]) -> Vec<String> {
    sets.iter()
        .enumerate()
        .map(|(i, terms)| {
            let (_, location, feature) = fabricated_record(i);
            format!(
                "del{i}\t{location}\t-\tENSG1\t{feature}\tTranscript\t{terms}\t{}",
                impact(terms)
            )
        })
        .collect()
}

fn plant(dir: &Path, name: &str, text: &str) -> std::path::PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).expect("a file is planted");
    path
}

/// The run of the comparator on the fabricated pair and the path of its report: the
/// reference a tab output carrying the sets `expected`, the candidate the rows carrying
/// `actual` written as Parquet, read beside `manifest` under `options`.
fn run_fabricated(
    dir: &Path,
    expected: &[&str],
    actual: &[&str],
    manifest: &str,
    options: &[&str],
) -> (Ran, std::path::PathBuf) {
    std::fs::create_dir_all(dir).expect("the directory is made");
    let reference = plant(
        dir,
        "vep.out",
        &format!("{COLUMNS}\n{}\n", fabricated_rows(expected).join("\n")),
    );
    let candidate = dir.join("candidate.parquet");
    let rows = fabricated_rows(actual);
    let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
    plant_tab_rows(&candidate, COLUMNS, &rows);
    let lines: String = (0..expected.len())
        .map(|i| format!("{}\n", fabricated_record(i).0))
        .collect();
    let input = plant(dir, "variants.vcf", &lines);
    let manifest = plant(dir, "manifest.json", manifest);
    let written = dir.join("report.json");
    let mut argv = args_of(&[
        "--format",
        "parquet",
        "--expected",
        path_text(&reference),
        "--actual",
        path_text(&candidate),
        "--report",
        path_text(&written),
        "--input",
        path_text(&input),
        "--manifest",
        path_text(&manifest),
    ]);
    argv.extend(args_of(options));
    // The comparison of the whole outputs first: the run through the stream of record
    // groups writes the same report, standard output, standard error and exit code. The
    // second run sets the variable to a value other than `1`, which asks for the stream
    // whatever the environment of the tests holds.
    let whole = run_in(dir, &[(WHOLE_VARIABLE.into(), "1".into())], &argv);
    let whole_report = std::fs::read(&written).ok();
    let _ = std::fs::remove_file(&written);
    let ran = run_in(dir, &[(WHOLE_VARIABLE.into(), "0".into())], &argv);
    assert_eq!(
        (ran.code, ran.out(), ran.err(), std::fs::read(&written).ok()),
        (whole.code, whole.out(), whole.err(), whole_report)
    );
    (ran, written)
}

/// The exit code and the report of `run_fabricated`, which must write one.
fn compare_fabricated(
    dir: &Path,
    expected: &[&str],
    actual: &[&str],
    manifest: &str,
    options: &[&str],
) -> (i32, Value) {
    let (ran, written) = run_fabricated(dir, expected, actual, manifest, options);
    assert!(written.is_file(), "no report: {}", ran.err());
    (ran.code, report(&written))
}

fn with<'a>(options: &[&'a str], more: &[&'a str]) -> Vec<&'a str> {
    options.iter().chain(more).copied().collect()
}

fn reasons(report: &Value) -> String {
    compact(member(report, &["reasons"]))
}

fn applied(report: &Value) -> &Value {
    member(report, &["masks", "applied"])
}

fn class_count(class: &str, count: i64) -> String {
    format!("{{{}: {count}}}", quoted(class))
}

#[test]
fn a_named_run_reads_its_own_items_alone_in_the_parquet_format() {
    let dir = scratch("parquet-named-run-alone");
    let manifest = runs_manifest();
    let own_pair: (&[&str], &[&str]) = (&[KEPT_SET, KEPT_SET], &[LOST_SET, KEPT_SET]);
    for name in ["default", "plain"] {
        let options = with(&RUNS_MASKS, &["--run", name]);
        let (status, report) =
            compare_fabricated(&dir.join(name), own_pair.0, own_pair.1, &manifest, &options);
        assert_eq!(status, 0, "{name}: {}", reasons(&report));
        assert_eq!(member(&report, &["verdict"]).as_str(), Some("green"));
        assert!(is(applied(&report), &class_count("some_class", 1)));
        assert!(is(
            member(&report, &["masks", "available"]),
            &class_count("some_class", 1)
        ));
        assert_eq!(member(&report, &["manifest_run"]).as_str(), Some(name));
    }
    // Unnamed, the pair would be held to the item of the run `other` as well.
    let (ran, written) = run_fabricated(
        &dir.join("unnamed"),
        own_pair.0,
        own_pair.1,
        &manifest,
        &RUNS_MASKS,
    );
    assert_eq!((ran.code, ran.err().as_str()), (2, UNNAMED_RUNS));
    assert!(!written.is_file());
}

#[test]
fn a_named_run_still_reports_its_own_item_that_is_not_observed_in_the_parquet_format() {
    let dir = scratch("parquet-named-run-not-observed");
    let manifest = runs_manifest();
    for (name, record) in [("default", 0), ("other", 1)] {
        let options = with(&RUNS_MASKS, &["--run", name]);
        let agreeing = [KEPT_SET, KEPT_SET];
        let (status, report) =
            compare_fabricated(&dir.join(name), &agreeing, &agreeing, &manifest, &options);
        assert_eq!(
            (status, reasons(&report).as_str()),
            (1, NOT_OBSERVED),
            "{name}"
        );
        let stale = member(&report, &["stale_divergences"])
            .as_list()
            .expect("the stale items");
        assert_eq!(stale.len(), 1);
        let item = stale[0].as_str().expect("a stale item's text");
        assert!(item.contains(&format!("ENST{}", record + 1)), "{item}");
    }
}

#[test]
fn a_difference_only_another_run_documents_is_not_set_aside_in_the_parquet_format() {
    let dir = scratch("parquet-other-run-difference");
    let manifest = runs_manifest();
    let both_differ: (&[&str], &[&str]) = (&[KEPT_SET, KEPT_SET], &[LOST_SET, LOST_SET]);
    let options = with(&RUNS_MASKS, &["--run", "default"]);
    let (status, report) = compare_fabricated(
        &dir.join("own"),
        both_differ.0,
        both_differ.1,
        &manifest,
        &options,
    );
    assert_eq!(status, 1);
    assert!(
        reasons(&report).contains(CONSEQUENCE_DIFFERS),
        "{}",
        reasons(&report)
    );
    assert!(is(applied(&report), &class_count("some_class", 1)));
    assert!(is(member(&report, &["stale_divergences"]), "[]"));
    // Unnamed, the item of the run `other` would set the second difference aside.
    let (ran, written) = run_fabricated(
        &dir.join("unnamed"),
        both_differ.0,
        both_differ.1,
        &manifest,
        &RUNS_MASKS,
    );
    assert_eq!((ran.code, ran.err().as_str()), (2, UNNAMED_RUNS));
    assert!(!written.is_file());
}

/// Crossing 1: under `--reference-sequence absent` an item that carries the with-reference
/// key documents its `vep_rs_consequence_sets`, and a candidate that wrote them is green
/// with the item counted.
#[test]
fn a_run_without_a_reference_sequence_that_writes_the_sets_without_one_observes_the_item_in_the_parquet_format(
) {
    let dir = scratch("parquet-crossing-1");
    let options = with(&REFERENCE_MASKS, &["--reference-sequence", "absent"]);
    let (status, report) = compare_fabricated(
        &dir,
        &[KEPT_SET],
        &[LOST_SET],
        &reference_manifest(),
        &options,
    );
    assert_eq!(status, 0, "{}", reasons(&report));
    assert_eq!(member(&report, &["verdict"]).as_str(), Some("green"));
    assert!(is(applied(&report), &class_count(REMOVED_CLASS, 1)));
    assert!(is(
        member(&report, &["counts", "documented_divergences_skipped"]),
        "1"
    ));
    assert_eq!(
        member(&report, &["reference_sequence"]).as_str(),
        Some("absent")
    );
}

/// Crossing 2: under `--reference-sequence present` an item whose with-reference sets are
/// the reference's own documents nothing: the row is compared like any other, no difference
/// is seen, and the item is neither skipped as documented nor reported as not observed.
#[test]
fn a_run_with_a_reference_sequence_that_writes_the_stated_sets_is_compared_exactly_in_the_parquet_format(
) {
    let dir = scratch("parquet-crossing-2");
    let options = with(&REFERENCE_MASKS, &["--reference-sequence", "present"]);
    let (status, report) = compare_fabricated(
        &dir,
        &[KEPT_SET],
        &[KEPT_SET],
        &reference_manifest(),
        &options,
    );
    assert_eq!(status, 0, "{}", reasons(&report));
    assert_eq!(member(&report, &["verdict"]).as_str(), Some("green"));
    assert!(is(member(&report, &["stale_divergences"]), "[]"));
    assert!(is(
        member(&report, &["counts", "documented_divergences_skipped"]),
        "0"
    ));
    assert!(is(member(&report, &["counts", "compared"]), "1"));
    assert_eq!(
        member(&report, &["reference_sequence"]).as_str(),
        Some("present")
    );
}

/// Crossing 3: a run without a reference sequence that wrote the reference's set leaves the
/// item not observed.
#[test]
fn a_run_without_a_reference_sequence_that_writes_the_with_reference_sets_is_red_in_the_parquet_format(
) {
    let dir = scratch("parquet-crossing-3");
    let options = with(&REFERENCE_MASKS, &["--reference-sequence", "absent"]);
    let (status, report) = compare_fabricated(
        &dir,
        &[KEPT_SET],
        &[KEPT_SET],
        &reference_manifest(),
        &options,
    );
    assert_eq!((status, reasons(&report).as_str()), (1, NOT_OBSERVED));
}

/// Crossing 4: a run with a reference sequence that wrote the sets of a run without one is
/// a plain difference of the Consequence. The item's class is a built-in row class too,
/// whose predicate admits the pair; it takes no row of a key the manifest states as equal.
#[test]
fn a_run_with_a_reference_sequence_that_writes_the_sets_without_one_is_red_in_the_parquet_format() {
    let dir = scratch("parquet-crossing-4");
    let present = with(&REFERENCE_MASKS, &["--reference-sequence", "present"]);
    let (status, report) = compare_fabricated(
        &dir,
        &[KEPT_SET],
        &[LOST_SET],
        &reference_manifest(),
        &present,
    );
    assert_eq!(status, 1);
    assert!(
        reasons(&report).contains(CONSEQUENCE_DIFFERS),
        "{}",
        reasons(&report)
    );
    assert!(is(applied(&report), "{}"));
    assert!(is(member(&report, &["stale_divergences"]), "[]"));
    // The same pair with no item on the key is the class's own: the predicate admits it.
    let (status, report) = compare_fabricated(
        &dir.join("unstated"),
        &[KEPT_SET],
        &[LOST_SET],
        &divergences(&[]),
        &present,
    );
    assert_eq!(status, 0, "{}", reasons(&report));
    assert!(is(applied(&report), &class_count(REMOVED_CLASS, 1)));
}

/// A call that states no setting reads an item's `vep_rs_consequence_sets` whether or not
/// the item carries the with-reference key, and an item without the key reads the same under
/// both settings: a candidate that wrote the documented set is green with the item counted,
/// one that wrote the reference's set leaves it not observed.
#[test]
fn an_item_without_the_key_and_a_call_without_the_setting_read_alike_in_the_parquet_format() {
    let dir = scratch("parquet-no-key-no-setting");
    let plain = divergences(&[fabricated_item(0, REMOVED_CLASS, "")]);
    let observed = |name: String, manifest: &str, options: &[&str]| {
        let (status, report) = compare_fabricated(
            &dir.join(&name),
            &[KEPT_SET],
            &[LOST_SET],
            manifest,
            options,
        );
        assert_eq!(status, 0, "{name}: {}", reasons(&report));
        assert!(
            is(applied(&report), &class_count(REMOVED_CLASS, 1)),
            "{name}"
        );
    };
    let unobserved = |name: String, manifest: &str, options: &[&str]| {
        let (status, report) = compare_fabricated(
            &dir.join(&name),
            &[KEPT_SET],
            &[KEPT_SET],
            manifest,
            options,
        );
        assert_eq!(
            (status, reasons(&report).as_str()),
            (1, NOT_OBSERVED),
            "{name}"
        );
    };
    for (n, manifest) in [reference_manifest(), plain.clone()].iter().enumerate() {
        observed(format!("observed{n}"), manifest, &REFERENCE_MASKS);
        unobserved(format!("unobserved{n}"), manifest, &REFERENCE_MASKS);
    }
    for (n, setting) in ["present", "absent"].into_iter().enumerate() {
        let options = with(&REFERENCE_MASKS, &["--reference-sequence", setting]);
        observed(format!("told{n}"), &plain, &options);
        unobserved(format!("told-unobserved{n}"), &plain, &options);
    }
}

/// The Parquet leg of two cases that walk every format: the applied count of a class is
/// the entries it sets aside, and with-reference sets other than the reference's are the
/// item's documented sets in a run with a reference sequence.
#[test]
fn an_item_with_reference_sets_is_counted_and_held_to_them_in_the_parquet_format() {
    let dir = scratch("parquet-with-reference-sets");
    let manifest = divergences(&[
        fabricated_item(0, REMOVED_CLASS, &with_reference(KEPT_SET)),
        fabricated_item(1, REMOVED_CLASS, ""),
    ]);
    let settings: [(Option<&str>, [&str; 2], i64); 3] = [
        (Some("present"), [KEPT_SET, LOST_SET], 1),
        (Some("absent"), [LOST_SET, LOST_SET], 2),
        (None, [LOST_SET, LOST_SET], 2),
    ];
    for (setting, actual, count) in settings {
        let options = match setting {
            Some(setting) => with(&REFERENCE_MASKS, &["--reference-sequence", setting]),
            None => REFERENCE_MASKS.to_vec(),
        };
        let name = setting.unwrap_or("untold");
        let (status, report) = compare_fabricated(
            &dir.join(name),
            &[KEPT_SET, KEPT_SET],
            &actual,
            &manifest,
            &options,
        );
        assert_eq!(status, 0, "{name}: {}", reasons(&report));
        assert!(is(
            member(&report, &["masks", "available"]),
            &class_count(REMOVED_CLASS, 2)
        ));
        assert!(is(member(&report, &["masks", "unused"]), "[]"));
        assert!(
            is(applied(&report), &class_count(REMOVED_CLASS, count)),
            "{name}"
        );
    }
    let manifest = divergences(&[fabricated_item(0, "some_class", &with_reference(OTHER_SET))]);
    let options = with(&RUNS_MASKS, &["--reference-sequence", "present"]);
    let (status, report) = compare_fabricated(
        &dir.join("observed"),
        &[KEPT_SET],
        &[OTHER_SET],
        &manifest,
        &options,
    );
    assert_eq!(status, 0, "{}", reasons(&report));
    assert!(is(applied(&report), &class_count("some_class", 1)));
    let (status, report) = compare_fabricated(
        &dir.join("crossed"),
        &[KEPT_SET],
        &[LOST_SET],
        &manifest,
        &options,
    );
    assert_eq!((status, reasons(&report).as_str()), (1, NOT_OBSERVED));
}

/// Mixed records, an unrelated transcript-only record and a motif-only record: the rows of
/// the reference, or of the candidate, which lacks the second record and writes another
/// term on every transcript.
fn regulatory_rows(candidate: bool) -> Vec<String> {
    let rows = [
        (
            0,
            "ENSR1",
            "RegulatoryFeature",
            "regulatory_region_variant",
            "-",
        ),
        (
            0,
            "ENSM1",
            "MotifFeature",
            "TF_binding_site_variant",
            "matrix1",
        ),
        (0, "ENST1", "Transcript", "missense_variant", "-"),
        (1, "ENST2", "Transcript", "intron_variant", "-"),
        (
            2,
            "ENSM2",
            "MotifFeature",
            "TF_binding_site_variant",
            "matrix1",
        ),
    ];
    rows.iter()
        .filter(|(record, ..)| !(candidate && *record == 1))
        .map(|(record, feature, kind, term, matrix)| {
            let term = if candidate && *kind == "Transcript" {
                "stop_gained"
            } else {
                term
            };
            format!("1:{}\tT\t{feature}\t{kind}\t{term}\t{matrix}", 100 + record)
        })
        .collect()
}

#[test]
fn cli_keep_only_regulatory_rows_in_the_parquet_format() {
    let dir = scratch("parquet-keep-only");
    let columns = "#Location\tAllele\tFeature\tFeature_type\tConsequence\tMOTIF_NAME";
    let inputs: String = (0..3)
        .map(|i| format!("1\t{}\tv{i}\tA\tT\t.\t.\t.\n", 100 + i))
        .collect();
    let input = plant(&dir, "input.vcf", &inputs);
    let expected = plant(
        &dir,
        "expected.txt",
        &format!("{columns}\n{}\n", regulatory_rows(false).join("\n")),
    );
    let actual = dir.join("actual.parquet");
    let planted = |rows: &[String]| {
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        plant_tab_rows(&actual, columns, &rows);
    };
    let original = regulatory_rows(true);
    planted(&original);
    let written = dir.join("report.json");
    let compared = |keep_only: bool| {
        let mut argv = args_of(&[
            "--format",
            "parquet",
            "--expected",
            path_text(&expected),
            "--actual",
            path_text(&actual),
            "--input",
            path_text(&input),
            "--report",
            path_text(&written),
        ]);
        if keep_only {
            argv.extend(args_of(&["--keep-only", "regulatory_feature_rows"]));
        }
        let ran = run(&dir, &argv);
        assert!(written.is_file(), "no report: {}", ran.err());
        (ran.code, report(&written))
    };
    assert_eq!(compared(false).0, 1);
    let (status, report) = compared(true);
    assert_eq!(status, 0, "{}", reasons(&report));
    assert_eq!(
        member(&report, &["keep_only"]).as_str(),
        Some("regulatory_feature_rows")
    );
    for count in ["compared", "expected_entries", "actual_entries"] {
        assert!(is(member(&report, &["counts", count]), "3"), "{count}");
    }
    assert!(is(
        member(&report, &["one_sided"]),
        "{\"expected_only\": 0, \"actual_only\": 0}"
    ));
    assert!(is(member(&report, &["masks", "requested"]), "[]"));
    // A motif-field disagreement is reported under the selector.
    let renamed: Vec<String> = original
        .iter()
        .map(|row| row.replace("matrix1", "matrix2"))
        .collect();
    planted(&renamed);
    let (status, report) = compared(true);
    assert_eq!(status, 1);
    assert!(is(
        member(&report, &["columns", "MOTIF_NAME", "differing"]),
        "2"
    ));
    // One-sided regulatory feature identities must not be masked out.
    let other: Vec<String> = original
        .iter()
        .map(|row| row.replace("ENSR1", "ENSR2"))
        .collect();
    planted(&other);
    let (status, report) = compared(true);
    assert_eq!(status, 1);
    assert!(is(
        member(&report, &["one_sided"]),
        "{\"expected_only\": 1, \"actual_only\": 1}"
    ));
}
