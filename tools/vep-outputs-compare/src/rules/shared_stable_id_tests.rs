// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `merged_same_stable_id_collapse` over fabricated discordant rows and a
//! collector filled by hand with what a replay found of the transcripts that share an id.

use super::*;
use crate::rules::scope::SharedStableId;
use crate::rules::table::{DIVERGENCE_CLASSES_116_2, MERGED_SAME_STABLE_ID_COLLAPSE};
use crate::TempDir;

const FILE: &str = "records.txt";

fn row(
    source: &str,
    location: &str,
    allele: &str,
    feature: &str,
    feature_type: &str,
    consequence_set: &str,
) -> Row {
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
fn extra(location: &str, feature: &str, consequence_set: &str) -> Row {
    row(
        "extra_in_rust",
        location,
        "G",
        feature,
        "Transcript",
        consequence_set,
    )
}

/// The reference's intergenic row.
fn intergenic(location: &str) -> Row {
    row(
        "missing_in_rust",
        location,
        "G",
        "-",
        "-",
        "intergenic_variant",
    )
}

/// A collector whose replay found, at each Location, transcripts of the stable id dropped
/// in reach of the record, and the kept one in reach or not.
fn scope(seen: &[(&str, &str, u64, bool)]) -> ReferenceRecordScope {
    let mut scope = ReferenceRecordScope::default();
    scope.record_shared_stable_id_drops(FILE, HashMap::new());
    for (location, stable_id, dropped, kept_in_reach) in seen {
        let seen = SharedStableId {
            dropped: *dropped,
            kept_in_reach: *kept_in_reach,
        };
        let one = HashMap::from([(
            location.to_string(),
            HashMap::from([(stable_id.to_string(), seen)]),
        )]);
        scope.record_shared_stable_id_drops(FILE, one);
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
    let stats = apply_shared_stable_id_rule(rows, &mut masked, scope, reference_dir)
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
fn the_rows_of_dropped_transcripts_and_the_intergenic_row_beside_them_are_set_aside() {
    let rows = [
        intergenic("contig:996672"),
        extra("contig:996672", "XM_900001.1", "stop_gained"),
        extra("contig:996672", "XM_900002.1", "stop_gained"),
        extra("contig:4919323", "XM_900004.1", "intron_variant"),
    ];
    let scope = scope(&[
        ("contig:996672", "XM_900001.1", 1, false),
        ("contig:996672", "XM_900002.1", 1, false),
        ("contig:4919323", "XM_900004.1", 1, false),
    ]);
    assert_eq!(pass(&rows, &scope), ((4, 4, 2, 1, 3), vec![0, 1, 2, 3]));
}

#[test]
fn a_comparison_that_names_no_cache_applies_no_such_class() {
    let rows = [extra("contig:996672", "XM_900001.1", "stop_gained")];
    let mut masked = HashSet::new();
    let unnamed = ReferenceRecordScope::default();
    let absent = Path::new("absent");
    let stats = apply_shared_stable_id_rule(&rows, &mut masked, &unnamed, Some(absent));
    assert!(stats.expect("no output is opened").is_none());
    assert!(masked.is_empty());
}

#[test]
fn a_row_of_an_id_the_replay_did_not_find_dropped_there_stays_open() {
    let rows = [
        extra("contig:501000", "XM_900009.1", "intron_variant"),
        extra("contig:501000", "XM_900001.1", "intron_variant"),
        extra("contig:996672", "XM_900001.1", "stop_gained"),
    ];
    // A transcript of the id is dropped at one Location and at no other.
    let scope = scope(&[("contig:996672", "XM_900001.1", 1, false)]);
    assert_eq!(pass(&rows, &scope), ((1, 1, 1, 0, 1), vec![2]));
    let other_file = Row {
        file_name: "other.txt".to_string(),
        ..extra("contig:996672", "XM_900001.1", "stop_gained")
    };
    assert_eq!(pass(&[other_file], &scope), ((0, 0, 0, 0, 0), vec![]));
}

#[test]
fn more_rows_of_an_id_than_transcripts_dropped_in_reach_stay_open() {
    let rows = [
        extra("contig:2999291", "XM_900003.1", "missense_variant"),
        extra("contig:2999291", "XM_900003.1", "synonymous_variant"),
    ];
    let one = scope(&[("contig:2999291", "XM_900003.1", 1, true)]);
    assert_eq!(pass(&rows, &one), ((0, 0, 0, 0, 0), vec![]));
    let two = scope(&[("contig:2999291", "XM_900003.1", 2, true)]);
    assert_eq!(pass(&rows, &two), ((2, 2, 1, 0, 2), vec![0, 1]));
    // The rows of another allele are counted apart.
    let rows = [
        extra("contig:2999291", "XM_900003.1", "missense_variant"),
        row(
            "extra_in_rust",
            "contig:2999291",
            "C",
            "XM_900003.1",
            "Transcript",
            "synonymous_variant",
        ),
    ];
    assert_eq!(pass(&rows, &one), ((2, 2, 1, 0, 2), vec![0, 1]));
}

#[test]
fn a_row_that_is_no_lone_transcript_row_of_vep_rs_stays_open() {
    let scope = scope(&[("contig:2999291", "XM_900003.1", 1, true)]);
    // The reference wrote a row of the id at the Location that vep-rs lacks.
    let rows = [
        extra("contig:2999291", "XM_900003.1", "missense_variant"),
        row(
            "missing_in_rust",
            "contig:2999291",
            "C",
            "XM_900003.1",
            "Transcript",
            "non_coding_transcript_exon_variant",
        ),
    ];
    assert_eq!(pass(&rows, &scope), ((0, 0, 0, 0, 0), vec![]));
    let regulatory = row(
        "extra_in_rust",
        "contig:2999291",
        "G",
        "XM_900003.1",
        "RegulatoryFeature",
        "regulatory_region_variant",
    );
    assert_eq!(pass(&[regulatory], &scope), ((0, 0, 0, 0, 0), vec![]));
    let reference_only = row(
        "missing_in_rust",
        "contig:2999291",
        "G",
        "XM_900003.1",
        "Transcript",
        "missense_variant",
    );
    assert_eq!(pass(&[reference_only], &scope), ((0, 0, 0, 0, 0), vec![]));
}

#[test]
fn a_row_another_rule_took_or_of_a_record_the_reference_skipped_is_left_to_it() {
    let rows = [extra("contig:996672", "XM_900001.1", "stop_gained")];
    let mut dropped = scope(&[("contig:996672", "XM_900001.1", 1, false)]);
    let mut masked = HashSet::from([0]);
    let stats = apply_shared_stable_id_rule(&rows, &mut masked, &dropped, None)
        .expect("no output is opened")
        .expect("a cache is named");
    assert_eq!((stats.count, stats.excluded, stats.records), (0, 0, 0));
    dropped
        .absent
        .insert((FILE.to_string(), "contig:996672".to_string()));
    assert_eq!(pass(&rows, &dropped), ((0, 0, 0, 0, 0), vec![]));
    dropped.absent.clear();
    dropped
        .gvcf_absent
        .insert((FILE.to_string(), "contig:996672".to_string()));
    assert_eq!(pass(&rows, &dropped), ((0, 0, 0, 0, 0), vec![]));
}

#[test]
fn the_intergenic_row_stays_open_beside_a_row_of_vep_rs_that_is_not_of_the_class() {
    let scope = scope(&[("contig:996672", "XM_900001.1", 1, false)]);
    // vep-rs wrote a second transcript's row there, which no dropped transcript explains.
    let rows = [
        intergenic("contig:996672"),
        extra("contig:996672", "XM_900001.1", "stop_gained"),
        extra("contig:996672", "XM_900009.1", "intron_variant"),
    ];
    assert_eq!(pass(&rows, &scope), ((1, 1, 1, 0, 1), vec![1]));
    // An intergenic row with no row of vep-rs at its Location and allele.
    let rows = [
        intergenic("contig:996672"),
        row(
            "extra_in_rust",
            "contig:996672",
            "C",
            "XM_900001.1",
            "Transcript",
            "stop_gained",
        ),
    ];
    assert_eq!(pass(&rows, &scope), ((1, 1, 1, 0, 1), vec![1]));
    // A row of the reference that is not its intergenic row.
    let rows = [
        row(
            "missing_in_rust",
            "contig:996672",
            "G",
            "ENSR1",
            "RegulatoryFeature",
            "regulatory_region_variant",
        ),
        extra("contig:996672", "XM_900001.1", "stop_gained"),
    ];
    assert_eq!(pass(&rows, &scope), ((1, 1, 1, 0, 1), vec![1]));
}

/// A record of a merged cache that the two passes together left with no transcript: one of
/// vep-rs's rows there is of a transcript the `dbID` pass dropped, which that class took,
/// and the other of a transcript this pass dropped. The intergenic row goes with this
/// class, whose pass ran last; a row another rule took for any other reason holds it open.
#[test]
fn the_intergenic_row_is_set_aside_beside_rows_of_both_passes() {
    let rows = [
        intergenic("contig:996672"),
        extra("contig:996672", "ENSEXT00000000001", "intron_variant"),
        extra("contig:996672", "XM_900001.1", "stop_gained"),
    ];
    let mut both = scope(&[("contig:996672", "XM_900001.1", 1, false)]);
    let of_the_dbid_pass = HashMap::from([(
        "contig:996672".to_string(),
        HashSet::from(["ENSEXT00000000001".to_string()]),
    )]);
    both.record_shared_dbid_drops(FILE, of_the_dbid_pass);
    let taken_before = |scope: &ReferenceRecordScope| {
        let mut masked = HashSet::from([1]);
        let stats = apply_shared_stable_id_rule(&rows, &mut masked, scope, None)
            .expect("no output is opened")
            .expect("a cache is named");
        let mut set_aside: Vec<usize> = masked.into_iter().collect();
        set_aside.sort_unstable();
        ((stats.count, stats.perl, stats.rust), set_aside)
    };
    assert_eq!(taken_before(&both), ((2, 1, 1), vec![0, 1, 2]));
    let stable_id_alone = scope(&[("contig:996672", "XM_900001.1", 1, false)]);
    assert_eq!(taken_before(&stable_id_alone), ((1, 0, 1), vec![1, 2]));
}

/// The pass over `discordant.tsv` applies the `dbID` class before this one: at a record the
/// two passes together left with no transcript this class then finds the other's row taken,
/// and takes the intergenic row.
#[test]
fn the_pass_over_the_discordant_rows_applies_the_dbid_class_first() {
    use crate::input::CacheContigNames;
    use crate::rules::filter::{filter_intended_divergences, Collector};
    use crate::rules::table::MERGED_SAME_DBID_COLLAPSE;
    let scratch = TempDir::create(&std::env::temp_dir(), "shared_stable_id_test_").unwrap();
    let discordant = scratch.path.join("discordant.tsv");
    let open = scratch.path.join("discordant_open.tsv");
    let columns = "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set\n";
    std::fs::write(
        &discordant,
        format!(
            "{columns}{FILE}\tmissing_in_rust\tcontig:996672\tG\t-\t-\tintergenic_variant\n\
             {FILE}\textra_in_rust\tcontig:996672\tG\tENSEXT00000000001\tTranscript\t\
             intron_variant\n\
             {FILE}\textra_in_rust\tcontig:996672\tG\tXM_900001.1\tTranscript\tstop_gained\n"
        ),
    )
    .unwrap();
    let mut both = scope(&[("contig:996672", "XM_900001.1", 1, false)]);
    let of_the_dbid_pass = HashMap::from([(
        "contig:996672".to_string(),
        HashSet::from(["ENSEXT00000000001".to_string()]),
    )]);
    both.record_shared_dbid_drops(FILE, of_the_dbid_pass);
    let pass = filter_intended_divergences(
        &discordant,
        &open,
        Collector::Release116(&both, &CacheContigNames::Unknown),
        &[],
        None,
        None,
    )
    .expect("no output is read");
    assert_eq!((pass.excluded_perl, pass.excluded_rust), (1, 2));
    let stated = |bucket: &str| {
        let class = pass.buckets.iter().find(|(name, _)| name == bucket);
        class.map(|(_, stats)| stats.clone())
    };
    let statistics = |count, perl, rust| {
        Some(vec![
            ("count", count),
            ("excluded", count),
            ("records", 1),
            ("perl", perl),
            ("rust", rust),
        ])
    };
    assert_eq!(
        stated(MERGED_SAME_DBID_COLLAPSE.bucket),
        statistics(1, 0, 1)
    );
    assert_eq!(
        stated(MERGED_SAME_STABLE_ID_COLLAPSE.bucket),
        statistics(2, 1, 1)
    );
    assert_eq!(std::fs::read_to_string(&open).unwrap(), columns);
}

/// One Transcript row of the default format as an output of the reference holds it.
fn output_row(location: &str, allele: &str, feature: &str) -> String {
    format!(
        "name\t{location}\t{allele}\tGENE\t{feature}\tTranscript\t\
         non_coding_transcript_exon_variant\t7\t-\t-\t-\t-\t-\tIMPACT=MODIFIER\n"
    )
}

/// The reference wrote a row of the id at the Location and vep-rs wrote it too, so it is
/// among no discordant rows: it is the kept transcript's where the replay found that one in
/// reach of the record, and vep-rs's further row is then the dropped one's. Where the kept
/// transcript is out of reach the reference's row is of no transcript the replay knows.
#[test]
fn a_row_of_the_id_that_both_engines_wrote_is_the_kept_transcripts_where_it_is_in_reach() {
    let scratch = TempDir::create(&std::env::temp_dir(), "shared_stable_id_test_").unwrap();
    let rows = [extra("contig:2999291", "XM_900003.1", "missense_variant")];
    let matched = output_row("contig:2999291", "G", "XM_900003.1");
    std::fs::write(scratch.path.join(FILE), matched).unwrap();
    let in_reach = scope(&[("contig:2999291", "XM_900003.1", 1, true)]);
    assert_eq!(
        pass_beside(&rows, &in_reach, Some(&scratch.path)),
        ((1, 1, 1, 0, 1), vec![0])
    );
    let out_of_reach = scope(&[("contig:2999291", "XM_900003.1", 1, false)]);
    assert_eq!(
        pass_beside(&rows, &out_of_reach, Some(&scratch.path)),
        ((0, 0, 0, 0, 0), vec![])
    );
    // Another transcript's row there, and the id's row at another Location, say nothing.
    let others = output_row("contig:2999291", "G", "XM_900009.1")
        + &output_row("contig:2999292", "G", "XM_900003.1");
    std::fs::write(scratch.path.join(FILE), others).unwrap();
    assert_eq!(
        pass_beside(&rows, &out_of_reach, Some(&scratch.path)),
        ((1, 1, 1, 0, 1), vec![0])
    );
    // An output is opened where a row of a dropped transcript asks, and nowhere else.
    let absent = scratch.path.join("absent");
    let mut masked = HashSet::new();
    assert!(apply_shared_stable_id_rule(&rows, &mut masked, &in_reach, Some(&absent)).is_err());
    let unasked = [extra("contig:501000", "XM_900009.1", "intron_variant")];
    assert_eq!(
        pass_beside(&unasked, &in_reach, Some(&absent)),
        ((0, 0, 0, 0, 0), vec![])
    );
}

/// Two records of one Location share the key of their rows, so the collector keeps of an
/// id what both saw: the fewer transcripts dropped, and the kept one in reach only where
/// it is for both.
#[test]
fn what_two_records_of_one_location_see_of_an_id_is_kept_where_it_holds_for_both() {
    let scope = scope(&[
        ("contig:2999291", "XM_900003.1", 2, true),
        ("contig:2999291", "XM_900003.1", 1, false),
        ("contig:2999291", "XM_900005.1", 1, true),
    ]);
    let at = &scope.shared_stable_id_drops.as_ref().unwrap()
        [&(FILE.to_string(), "contig:2999291".to_string())];
    let seen = |dropped, kept_in_reach| SharedStableId {
        dropped,
        kept_in_reach,
    };
    assert_eq!(at["XM_900003.1"], seen(1, false));
    assert_eq!(at["XM_900005.1"], seen(1, true));
}

/// The taxonomy of the reference registry among the documentation pages, by a path that
/// holds from the root of every crate under `tools/` that compiles this module.
const TAXONOMY_PAGE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../docs/concordance-provenance/116.2-divergence-classes.md"
);

#[test]
fn the_class_states_a_definition_cites_its_source_and_stands_on_the_taxonomy_page() {
    let class = MERGED_SAME_STABLE_ID_COLLAPSE;
    assert!(class.excludes);
    assert_eq!(class.taxonomy_class, class.bucket);
    assert!(class.definition.ends_with('.'));
    assert!(class.definition.split_whitespace().count() >= 12);
    assert!(class
        .perl_citation
        .contains("AnnotationType/Transcript.pm:275, 303-311, 314-343"));
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
