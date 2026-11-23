// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `merged_same_dbid_collapse` over fabricated discordant rows and a collector
//! filled by hand with the transcripts a replay found dropped.

use super::*;
use crate::rules::table::{DIVERGENCE_CLASSES_116_2, MERGED_SAME_DBID_COLLAPSE};
use crate::TempDir;

const FILE: &str = "records.txt";

fn row(source: &str, location: &str, allele: &str, feature: &str, feature_type: &str) -> Row {
    let consequence_set = if feature == "-" {
        "intergenic_variant"
    } else {
        "non_coding_transcript_exon_variant"
    };
    Row {
        file_name: FILE.to_string(),
        source: source.to_string(),
        location: location.to_string(),
        allele: allele.to_string(),
        feature: feature.to_string(),
        feature_type: feature_type.to_string(),
        consequence_set: consequence_set.to_string(),
        fields: Vec::new(),
    }
}

/// A row vep-rs alone writes on a transcript.
fn extra(location: &str, feature: &str) -> Row {
    row("extra_in_rust", location, "G", feature, "Transcript")
}

/// The reference's intergenic row.
fn intergenic(location: &str) -> Row {
    row("missing_in_rust", location, "G", "-", "-")
}

/// A collector whose replay dropped each stable id at its Location.
fn scope(dropped: &[(&str, &str)]) -> ReferenceRecordScope {
    let mut scope = ReferenceRecordScope::default();
    scope.record_shared_dbid_drops(FILE, HashMap::new());
    for (location, stable_id) in dropped {
        let one = HashMap::from([(location.to_string(), HashSet::from([stable_id.to_string()]))]);
        scope.record_shared_dbid_drops(FILE, one);
    }
    scope
}

/// The statistics of a pass as (count, excluded, records, perl, rust), and the rows it
/// set aside.
type Outcome = ((u64, u64, u64, u64, u64), Vec<usize>);

/// A pass beside no output of the reference.
fn pass(rows: &[Row], scope: &ReferenceRecordScope) -> Outcome {
    pass_beside(rows, scope, None)
}

/// A pass beside the reference's outputs under `reference_dir`.
fn pass_beside(
    rows: &[Row],
    scope: &ReferenceRecordScope,
    reference_dir: Option<&Path>,
) -> Outcome {
    let mut masked = HashSet::new();
    let stats = apply_shared_dbid_rule(rows, &mut masked, scope, reference_dir)
        .expect("the outputs are read")
        .expect("a cache is named");
    let mut set_aside: Vec<usize> = masked.into_iter().collect();
    set_aside.sort_unstable();
    (
        (
            stats.count,
            stats.excluded,
            stats.records,
            stats.perl,
            stats.rust,
        ),
        set_aside,
    )
}

#[test]
fn the_rows_of_a_dropped_transcript_and_the_intergenic_row_beside_them_are_set_aside() {
    let rows = [
        extra("contig:101500", "ENSEXT00000000001"),
        intergenic("contig:3201000"),
        extra("contig:3201000", "XR_000002.1"),
    ];
    let scope = scope(&[
        ("contig:101500", "ENSEXT00000000001"),
        ("contig:3201000", "XR_000002.1"),
    ]);
    assert_eq!(pass(&rows, &scope), ((3, 3, 2, 1, 2), vec![0, 1, 2]));
}

#[test]
fn a_comparison_that_names_no_cache_applies_no_such_class() {
    let rows = [extra("contig:101500", "ENSEXT00000000001")];
    let mut masked = HashSet::new();
    let unnamed = ReferenceRecordScope::default();
    let absent = Path::new("absent");
    let stats = apply_shared_dbid_rule(&rows, &mut masked, &unnamed, Some(absent));
    assert!(stats.expect("no output is opened").is_none());
    assert!(masked.is_empty());
}

#[test]
fn a_row_of_a_transcript_the_replay_did_not_find_dropped_stays_open() {
    let rows = [
        extra("contig:501000", "ENSEXT00000000009"),
        extra("contig:501000", "ENSEXT00000000001"),
        extra("contig:101500", "ENSEXT00000000001"),
    ];
    // The transcript is dropped at one Location and at no other.
    let scope = scope(&[("contig:101500", "ENSEXT00000000001")]);
    assert_eq!(pass(&rows, &scope), ((1, 1, 1, 0, 1), vec![2]));
    let other_file = Row {
        file_name: "other.txt".to_string(),
        ..extra("contig:101500", "ENSEXT00000000001")
    };
    assert_eq!(pass(&[other_file], &scope), ((0, 0, 0, 0, 0), vec![]));
}

#[test]
fn a_row_that_is_no_lone_transcript_row_of_vep_rs_stays_open() {
    let scope = scope(&[("contig:101500", "ENSEXT00000000001")]);
    // The reference wrote a row of the transcript at the Location, under another allele.
    let rows = [
        extra("contig:101500", "ENSEXT00000000001"),
        row(
            "missing_in_rust",
            "contig:101500",
            "C",
            "ENSEXT00000000001",
            "Transcript",
        ),
    ];
    assert_eq!(pass(&rows, &scope), ((0, 0, 0, 0, 0), vec![]));
    let regulatory = row(
        "extra_in_rust",
        "contig:101500",
        "G",
        "ENSEXT00000000001",
        "RegulatoryFeature",
    );
    assert_eq!(pass(&[regulatory], &scope), ((0, 0, 0, 0, 0), vec![]));
    let reference_only = row(
        "missing_in_rust",
        "contig:101500",
        "G",
        "ENSEXT00000000001",
        "Transcript",
    );
    assert_eq!(pass(&[reference_only], &scope), ((0, 0, 0, 0, 0), vec![]));
}

#[test]
fn a_row_another_rule_took_or_of_a_record_the_reference_skipped_is_left_to_it() {
    let rows = [extra("contig:101500", "ENSEXT00000000001")];
    let mut dropped = scope(&[("contig:101500", "ENSEXT00000000001")]);
    let mut masked = HashSet::from([0]);
    let stats = apply_shared_dbid_rule(&rows, &mut masked, &dropped, None)
        .expect("no output is opened")
        .expect("a cache is named");
    assert_eq!((stats.count, stats.excluded, stats.records), (0, 0, 0));
    dropped
        .absent
        .insert((FILE.to_string(), "contig:101500".to_string()));
    assert_eq!(pass(&rows, &dropped), ((0, 0, 0, 0, 0), vec![]));
    dropped.absent.clear();
    dropped
        .gvcf_absent
        .insert((FILE.to_string(), "contig:101500".to_string()));
    assert_eq!(pass(&rows, &dropped), ((0, 0, 0, 0, 0), vec![]));
}

#[test]
fn the_intergenic_row_stays_open_beside_a_row_of_vep_rs_that_is_not_of_the_class() {
    let scope = scope(&[("contig:3201000", "XR_000002.1")]);
    // vep-rs wrote a second transcript's row there, which no dropped transcript explains.
    let rows = [
        intergenic("contig:3201000"),
        extra("contig:3201000", "XR_000002.1"),
        extra("contig:3201000", "ENSEXT00000000009"),
    ];
    assert_eq!(pass(&rows, &scope), ((1, 1, 1, 0, 1), vec![1]));
    // An intergenic row with no row of vep-rs at its Location and allele.
    let rows = [
        intergenic("contig:3201000"),
        row(
            "extra_in_rust",
            "contig:3201000",
            "C",
            "XR_000002.1",
            "Transcript",
        ),
    ];
    assert_eq!(pass(&rows, &scope), ((1, 1, 1, 0, 1), vec![1]));
    // A row of the reference that is not its intergenic row.
    let rows = [
        row(
            "missing_in_rust",
            "contig:3201000",
            "G",
            "ENSR1",
            "RegulatoryFeature",
        ),
        extra("contig:3201000", "XR_000002.1"),
    ];
    assert_eq!(pass(&rows, &scope), ((1, 1, 1, 0, 1), vec![1]));
}

/// One Transcript row of the default format as an output of the reference holds it.
fn output_row(location: &str, allele: &str, feature: &str) -> String {
    format!(
        "name\t{location}\t{allele}\tGENE\t{feature}\tTranscript\t\
         non_coding_transcript_exon_variant\t7\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
    )
}

#[test]
fn a_row_of_the_transcript_that_both_engines_wrote_at_the_location_keeps_the_rows_open() {
    let scratch = TempDir::create(&std::env::temp_dir(), "shared_dbid_test_").unwrap();
    let rows = [
        intergenic("contig:3201000"),
        extra("contig:3201000", "XR_000002.1"),
    ];
    let scope = scope(&[("contig:3201000", "XR_000002.1")]);
    let beside = |output: String| {
        std::fs::write(scratch.path.join(FILE), output).unwrap();
        pass_beside(&rows, &scope, Some(&scratch.path))
    };
    // The reference wrote the transcript's row there under another allele, as vep-rs did.
    let matched = output_row("contig:3201000", "C", "XR_000002.1");
    assert_eq!(beside(matched), ((0, 0, 0, 0, 0), vec![]));
    // Another transcript's row there, and the transcript's row at another Location.
    let others = output_row("contig:3201000", "C", "ENSEXT00000000009")
        + &output_row("contig:3201001", "G", "XR_000002.1");
    assert_eq!(beside(others), ((2, 2, 1, 1, 1), vec![0, 1]));
    // An output is opened where a row of a dropped transcript asks, and nowhere else.
    let absent = scratch.path.join("absent");
    let mut masked = HashSet::new();
    assert!(apply_shared_dbid_rule(&rows, &mut masked, &scope, Some(&absent)).is_err());
    let unasked = [extra("contig:501000", "ENSEXT00000000009")];
    assert_eq!(
        pass_beside(&unasked, &scope, Some(&absent)),
        ((0, 0, 0, 0, 0), vec![])
    );
}

/// A reference output the class asks for and cannot read ends the pass over
/// `discordant.tsv` on a failure that names that output.
#[test]
fn an_output_the_class_cannot_read_is_named_in_the_failure_of_the_pass() {
    use crate::input::CacheContigNames;
    use crate::rules::filter::{filter_intended_divergences, Collector};
    let scratch = TempDir::create(&std::env::temp_dir(), "shared_dbid_test_").unwrap();
    let discordant = scratch.path.join("discordant.tsv");
    std::fs::write(
        &discordant,
        format!(
            "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n\
             {FILE}\textra_in_rust\tcontig:3201000\tG\tXR_000002.1\tTranscript\t\
             non_coding_transcript_exon_variant\n"
        ),
    )
    .unwrap();
    let scope = scope(&[("contig:3201000", "XR_000002.1")]);
    let absent = scratch.path.join("absent");
    let pass = filter_intended_divergences(
        &discordant,
        &scratch.path.join("discordant_open.tsv"),
        Collector::Release116(&scope, &CacheContigNames::Unknown),
        &[],
        Some((&absent, &absent)),
        None,
    );
    let Err(unread) = pass else {
        panic!("the pass read an output that is not there");
    };
    let output = absent.join(FILE);
    let missing = File::open(&output).unwrap_err();
    assert_eq!(
        unread.to_string(),
        format!("{}: {missing}", output.display())
    );
}

/// The taxonomy of the reference registry among the documentation pages, by a path that
/// holds from the root of every crate under `tools/` that compiles this module.
const TAXONOMY_PAGE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../docs/concordance-provenance/116.2-divergence-classes.md"
);

#[test]
fn the_class_states_a_definition_cites_its_source_and_stands_on_the_taxonomy_page() {
    let class = MERGED_SAME_DBID_COLLAPSE;
    assert!(class.excludes);
    assert_eq!(class.taxonomy_class, class.bucket);
    assert!(class.definition.ends_with('.'));
    assert!(class.definition.split_whitespace().count() >= 12);
    assert!(class
        .perl_citation
        .contains("AnnotationType/Transcript.pm:283-287"));
    assert!(DIVERGENCE_CLASSES_116_2
        .iter()
        .all(|held| held.bucket != class.bucket && held.taxonomy_class != class.taxonomy_class));
    let page =
        std::fs::read_to_string(TAXONOMY_PAGE).unwrap_or_else(|e| panic!("{TAXONOMY_PAGE}: {e}"));
    assert!(page.contains(&format!("`{}`", class.bucket)));
    let heading = "The classes that read the converted cache";
    assert!(
        page.contains(&format!(". {heading}\n")),
        "the page has no section {heading:?}"
    );
}
