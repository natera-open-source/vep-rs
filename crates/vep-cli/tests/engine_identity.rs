// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Engine identity: every output format names the vep-rs version beside the
//! Ensembl VEP header lines it reproduces (`## vep-rs <version>` in the default
//! and tab headers, `##vep-rs=<version>` in the VCF header, `vep_rs_version` in
//! every JSON object), and the golden-corpus comparators set aside exactly those
//! lines, so an output that lacks one, repeats one or names another version is
//! reported like any other header difference.

mod common;

use common::{
    compare_column_descriptions, compare_json_identity, compare_tab_headers, compare_vcf_headers,
    corpora, default_run, read_expected, run_vep, runs,
};

/// The version the binary under test was built from: every crate of the
/// workspace carries the release as its version, this test crate included.
const VERSION: &str = env!("CARGO_PKG_VERSION");

fn header_lines(text: &str) -> Vec<&str> {
    text.lines().filter(|l| l.starts_with('#')).collect()
}

#[test]
fn default_tab_vcf_json_carry_vep_rs_version() {
    let corpus = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-release116")
        .expect("the GRCh38-release116 corpus is committed");
    let run = default_run(&corpus);
    let tmp = tempfile::tempdir().unwrap();

    let tab_identity = format!("## vep-rs {VERSION}");
    for format in ["default", "tab"] {
        let text = run_vep(&run, format, tmp.path(), &[]);
        let lines = header_lines(&text);
        assert!(
            lines[0].starts_with("## ENSEMBL VARIANT EFFECT PREDICTOR v"),
            "{format}: the Ensembl VEP banner stays the first line"
        );
        assert_eq!(
            lines[1], tab_identity,
            "{format}: the identity line is second"
        );
        assert_eq!(
            lines.iter().filter(|l| **l == tab_identity).count(),
            1,
            "{format}: one identity line"
        );
    }

    let vcf = run_vep(&run, "vcf", tmp.path(), &[]);
    let lines = header_lines(&vcf);
    let vep_at = lines
        .iter()
        .position(|l| l.starts_with("##VEP="))
        .expect("the ##VEP= line");
    assert_eq!(lines[vep_at + 1], format!("##vep-rs={VERSION}"));
    assert_eq!(
        lines.iter().filter(|l| l.starts_with("##vep-rs=")).count(),
        1
    );

    let json = run_vep(&run, "json", tmp.path(), &[]);
    let objects: Vec<serde_json::Value> = json
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    // One object per record the run writes: every record without a reference
    // warning, and of the warned ones those vep-rs annotates where Ensembl VEP
    // dropped them (the input line names the record).
    let records = corpus.manifest["records"].as_array().unwrap();
    let unwarned = records
        .iter()
        .filter(|r| r["reference_warning"].is_null())
        .count();
    assert!(
        (unwarned..=records.len()).contains(&objects.len()),
        "one object per written record: {} objects for {} records, {unwarned} unwarned",
        objects.len(),
        records.len()
    );
    for obj in &objects {
        assert!(
            records.iter().any(|r| r["vcf"] == obj["input"]),
            "the object names a record of the corpus: {obj}"
        );
        assert_eq!(obj["vep_rs_version"], VERSION, "{obj}");
    }
}

/// The comparators accept exactly one identity line of exactly this version
/// and report everything else: a missing line, a repeated one, another
/// version, and an object without the key.
#[test]
fn comparators_set_aside_exactly_the_identity_lines() {
    let corpus = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-release116")
        .expect("the GRCh38-release116 corpus is committed");
    let run = default_run(&corpus);
    let tmp = tempfile::tempdir().unwrap();

    let tab_identity = format!("## vep-rs {VERSION}");
    let expected = read_expected(&corpus.dir, "tab.txt");
    let actual = run_vep(&run, "tab", tmp.path(), &[]);
    assert_eq!(compare_tab_headers(&expected, &actual), "");
    assert_eq!(
        compare_tab_headers(
            &expected,
            &actual.replacen(&tab_identity, "## vep-rs 0.0.0", 1)
        ),
        format!("vep-rs writes no {tab_identity:?} header line")
    );
    assert_eq!(
        compare_tab_headers(&expected, &actual.replacen(&tab_identity, "", 1)),
        format!("vep-rs writes no {tab_identity:?} header line")
    );
    let doubled = actual.replacen(&tab_identity, &format!("{tab_identity}\n{tab_identity}"), 1);
    assert_eq!(
        compare_tab_headers(&expected, &doubled),
        format!("header line 2: expected prefix \"## Output produced at \", got {tab_identity:?}"),
        "a repeated identity line is an extra line"
    );

    let vcf_identity = format!("##vep-rs={VERSION}");
    let expected = read_expected(&corpus.dir, "vcf.vcf");
    let actual = run_vep(&run, "vcf", tmp.path(), &[]);
    assert_eq!(compare_vcf_headers(&expected, &actual), "");
    assert_eq!(
        compare_vcf_headers(
            &expected,
            &actual.replacen(&vcf_identity, "##vep-rs=0.0.0", 1)
        ),
        format!("vep-rs writes no {vcf_identity:?} header line")
    );
    let doubled = actual.replacen(&vcf_identity, &format!("{vcf_identity}\n{vcf_identity}"), 1);
    assert!(
        compare_vcf_headers(&expected, &doubled).starts_with("header line count:"),
        "a repeated identity line is an extra line"
    );

    let actual = run_vep(&run, "json", tmp.path(), &[]);
    assert_eq!(compare_json_identity(&actual), "");
    let version_key = format!("\"vep_rs_version\":\"{VERSION}\"");
    assert!(
        actual.contains(&version_key),
        "{version_key} is in every object"
    );
    let other = actual.replacen(&version_key, "\"vep_rs_version\":\"0.0.0\"", 1);
    assert_eq!(
        compare_json_identity(&other),
        format!("JSON object 0: vep_rs_version is Some(\"0.0.0\"), expected {VERSION:?}")
    );
    let without = actual.replacen(&format!(",{version_key}"), "", 1);
    assert_eq!(
        compare_json_identity(&without),
        format!("JSON object 0: vep_rs_version is None, expected {VERSION:?}")
    );
}

/// The `vrs-db` flag set of the GRCh38-identifiers corpus has a database
/// reference, whose tab header compares by its column descriptions and column
/// line alone (`compare_column_descriptions`): silent on the committed pair, and
/// a report for a missing identity line, a changed description, a changed
/// column line and a dropped column line.
#[test]
fn column_descriptions_comparator_reports_each_difference() {
    let corpus = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-identifiers")
        .expect("the GRCh38-identifiers corpus is committed");
    let run = runs(&corpus)
        .into_iter()
        .find(|r| r.name == "vrs-db")
        .expect("the vrs-db flag set is in the manifest");
    assert_eq!(run.perl_mode, "database");
    let tmp = tempfile::tempdir().unwrap();

    let expected = run.expected("tab.txt");
    let actual = run_vep(&run, "tab", tmp.path(), &[]);
    assert_eq!(compare_column_descriptions(&expected, &actual), "");

    let tab_identity = format!("## vep-rs {VERSION}");
    assert_eq!(
        compare_column_descriptions(&expected, &actual.replacen(&tab_identity, "", 1)),
        format!("vep-rs writes no {tab_identity:?} header line")
    );

    let strand_description = "## STRAND : Strand of the feature (1/-1)";
    assert!(
        actual.contains(strand_description),
        "the tab header describes STRAND"
    );
    let changed_description = actual.replacen(strand_description, "## STRAND : Strand", 1);
    let report = compare_column_descriptions(&expected, &changed_description);
    assert!(
        report.starts_with("header line ")
            && report.contains(&format!(
                "expected {strand_description:?}, got \"## STRAND : Strand\""
            )),
        "{report}"
    );

    let changed_columns = actual.replacen("\tGA4GH_VRS", "\tGA4GH_VRX", 1);
    let report = compare_column_descriptions(&expected, &changed_columns);
    assert!(
        report.starts_with("header line ") && report.contains("GA4GH_VRX"),
        "{report}"
    );

    let without_columns: String = actual
        .lines()
        .filter(|l| !l.starts_with("#Uploaded_variation"))
        .map(|l| format!("{l}\n"))
        .collect();
    let report = compare_column_descriptions(&expected, &without_columns);
    assert!(
        report.starts_with("header column lines: expected "),
        "{report}"
    );
}
