// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison of the `pair_*` fixture: two outputs of five keys each that share four,
//! with one field and two Extra keys that differ and one key the reference side repeats.

use super::*;
use report::Json;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

/// The report over the pair with `--examples 5` and the sort intermediates under the
/// default temporary directory.
fn pair_report(keep_only: Option<KeepOnly>) -> Json {
    let args = Args {
        perl: fixture("pair_reference.txt"),
        rust: fixture("pair_scored.txt"),
        out: PathBuf::new(),
        examples: 5,
        tmp_dir: None,
        keep_only,
    };
    compare(&args).unwrap().1
}

fn at<'a>(json: &'a Json, path: &[&str]) -> &'a Json {
    path.iter().fold(json, |node, key| match node {
        Json::Object(items) => items
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value)
            .unwrap_or_else(|| panic!("no key {key:?}")),
        other => panic!("{key:?} asked of {other:?}"),
    })
}

fn object<const N: usize>(items: [(&str, Json); N]) -> Json {
    Json::Object(
        items
            .into_iter()
            .map(|(name, value)| (name.to_string(), value))
            .collect(),
    )
}

fn mismatch(reference: &str, scored: &str, count: u64) -> Json {
    object([
        ("reference", Json::Str(reference.into())),
        ("scored", Json::Str(scored.into())),
        ("count", Json::Int(count)),
    ])
}

#[test]
fn key_totals_and_duplicates() {
    let report = pair_report(None);
    assert_eq!(
        at(&report, &["totals"]),
        &object([
            ("reference_keys", Json::Int(5)),
            ("scored_keys", Json::Int(5)),
            ("intersection", Json::Int(4)),
            // 1:500 is on the reference side alone and 1:600 on the scored side alone.
            ("reference_only", Json::Int(1)),
            ("scored_only", Json::Int(1)),
            // The second 1:500 row spells the same consequence set in another order.
            ("reference_duplicate_keys", Json::Int(1)),
            ("scored_duplicate_keys", Json::Int(0)),
        ])
    );
}

#[test]
fn plain_fields_compare_only_on_intersection() {
    let report = pair_report(None);
    let uploaded = at(&report, &["fields", "Uploaded_variation"]);
    assert_eq!(at(uploaded, &["compared"]), &Json::Int(4));
    assert_eq!(at(uploaded, &["differ"]), &Json::Int(1));
    assert_eq!(
        at(uploaded, &["top_mismatches"]),
        &Json::Array(vec![mismatch("rs1", "1_100_C/T", 1)])
    );
    for name in [
        "Gene",
        "cDNA_position",
        "CDS_position",
        "Protein_position",
        "Amino_acids",
        "Codons",
        "Existing_variation",
    ] {
        assert_eq!(
            at(&report, &["fields", name]),
            &object([
                ("compared", Json::Int(4)),
                ("equal", Json::Int(4)),
                ("differ", Json::Int(0)),
                ("agreement", Json::Float(1.0)),
                ("distinct_mismatch_shapes", Json::Int(0)),
                ("top_mismatches", Json::Array(vec![])),
            ]),
            "{name}"
        );
    }
}

#[test]
fn extra_keys_compared_individually_with_absence_counted() {
    let report = pair_report(None);
    let extra = at(&report, &["extra"]);
    assert_eq!(at(extra, &["IMPACT", "equal"]), &Json::Int(4));
    assert_eq!(at(extra, &["IMPACT", "compared"]), &Json::Int(4));
    assert_eq!(at(extra, &["STRAND", "equal"]), &Json::Int(4));
    // DISTANCE is on one shared tuple on both sides and differs there.
    assert_eq!(
        at(extra, &["DISTANCE"]),
        &object([
            ("compared", Json::Int(1)),
            ("equal", Json::Int(0)),
            ("differ", Json::Int(1)),
            ("agreement", Json::Float(0.0)),
            ("distinct_mismatch_shapes", Json::Int(1)),
            (
                "top_mismatches",
                Json::Array(vec![mismatch("1500", "1499", 1)])
            ),
        ])
    );
    // FLAGS is on the scored side alone, so the reference value is the absent marker.
    let absent = String::from_utf8_lossy(merge::ABSENT);
    assert_eq!(absent, "<absent>");
    assert_eq!(
        at(extra, &["FLAGS", "top_mismatches"]),
        &Json::Array(vec![mismatch(&absent, "cds_end_NF", 1)])
    );
}

#[test]
fn keep_only_regulatory_rows_on_transcript_only_inputs_is_empty() {
    let report = pair_report(Some(KeepOnly::RegulatoryFeatureRows));
    let Json::Object(totals) = at(&report, &["totals"]) else {
        panic!("totals is an object");
    };
    assert_eq!(totals.len(), 7);
    for (name, count) in totals {
        assert_eq!(count, &Json::Int(0), "{name}");
    }
    assert_eq!(at(&report, &["fields", "Gene", "agreement"]), &Json::Null);
}
