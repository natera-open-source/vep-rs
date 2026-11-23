// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The one-sided pass of each reference release over fabricated discordant rows: the
//! transcript-selection class of the 115.2 reference under its two guards, and the record
//! classes of the 116.2 reference over a collector filled by hand.

use super::*;
use crate::input::KeyPrep;
use crate::rules::row_pair_rules;
use crate::rules::scope::{Breakends, Group};
use crate::TempDir;

const COLUMNS: &str = "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set";

/// The span of a deletion above the size past which Ensembl VEP marks a record skipped.
const GIANT: &str = "18:62909529-80261528";
/// A span of ten kilobases, below that size.
const SMALL: &str = "18:62909529-62919528";

/// What one pass leaves: the rows set aside per side, every class's statistics in report
/// order and the rows that stay open.
struct Pass {
    excluded: (u64, u64),
    buckets: Vec<(String, BucketStats)>,
    open_rows: Vec<String>,
}

impl Pass {
    fn class(&self, name: &str) -> Option<&BucketStats> {
        self.buckets
            .iter()
            .find(|(bucket, _)| bucket == name)
            .map(|(_, stats)| stats)
    }

    /// The rows the class recognised.
    fn count(&self, name: &str) -> u64 {
        self.class(name).expect("the class")[0].1
    }

    /// One column of every open row.
    fn open_column(&self, column: usize) -> Vec<&str> {
        self.open_rows
            .iter()
            .map(|row| row.split('\t').nth(column).expect("the column"))
            .collect()
    }
}

/// Writes `rows` as a discordance report and runs the pass over it, no output directory
/// given.
fn pass_over(rows: &[String], collector: Collector, row_rules: &[&ExclusionRule]) -> Pass {
    let scratch = TempDir::create(&std::env::temp_dir(), "one_sided_pass_test_").unwrap();
    let discordant = scratch.path.join("discordant.tsv");
    let open = scratch.path.join("discordant_open.tsv");
    std::fs::write(&discordant, format!("{COLUMNS}\n{}", rows.concat())).unwrap();
    let exclusions =
        filter_intended_divergences(&discordant, &open, collector, row_rules, None, None).unwrap();
    let open_rows = std::fs::read_to_string(&open)
        .unwrap()
        .lines()
        .skip(1)
        .map(str::to_string)
        .collect();
    Pass {
        excluded: (exclusions.excluded_perl, exclusions.excluded_rust),
        buckets: exclusions.buckets,
        open_rows,
    }
}

/// A consequence set as the key of a tuple carries it.
fn normalized(consequence: &str) -> String {
    let mut set = Vec::new();
    KeyPrep::default().normalized_set(consequence, &mut set);
    String::from_utf8(set).unwrap()
}

fn names(items: &[&str]) -> HashSet<String> {
    items.iter().map(|item| item.to_string()).collect()
}

/// One discordant row of a deletion on a transcript.
fn deletion_row(source: &str, location: &str, feature: &str, consequence: &str) -> String {
    let set = normalized(consequence);
    format!("f\t{source}\t{location}\tdeletion\t{feature}\tTranscript\t{set}\n")
}

/// One discordant row of a deletion that removes its transcript.
fn ablation_row(source: &str, location: &str, feature: &str) -> String {
    deletion_row(source, location, feature, "transcript_ablation")
}

/// The 115.2 collector holding one deletion record with the transcripts each engine named.
fn selection_scope(location: &str, perl: &[&str], rust: &[&str]) -> ScopeFeatures {
    let record = (location.to_string(), "deletion".to_string());
    ScopeFeatures {
        records: HashMap::from([(record, (names(perl), names(rust)))]),
    }
}

fn selection_pass(rows: &[String], scope: &ScopeFeatures) -> Pass {
    pass_over(rows, Collector::Release115(scope), &[])
}

fn selection(count: u64, excluded: u64, perl: u64, rust: u64) -> BucketStats {
    vec![
        ("count", count),
        ("excluded", excluded),
        ("perl", perl),
        ("rust", rust),
    ]
}

/// One discordant row under the 116.2 reference; a row on no feature is an intergenic one.
fn row(source: &str, location: &str, allele: &str, feature: &str, consequence: &str) -> String {
    let feature_type = if feature == "-" { "-" } else { "Transcript" };
    let set = normalized(consequence);
    format!("f\t{source}\t{location}\t{allele}\t{feature}\t{feature_type}\t{set}\n")
}

fn groups(locations: &[&str]) -> HashSet<Group> {
    locations
        .iter()
        .map(|location| ("f".to_string(), location.to_string()))
        .collect()
}

/// The 116.2 collector holding one breakend Location: the alleles each side wrote there and
/// the mate chromosomes its bracket alleles name.
fn breakend_scope(
    location: &str,
    perl: &[&str],
    rust: &[&str],
    mates: &[&str],
) -> ReferenceRecordScope {
    let mut scope = ReferenceRecordScope::default();
    scope.breakends.insert(
        ("f".to_string(), location.to_string()),
        Breakends {
            perl_alleles: names(perl),
            rust_alleles: names(rust),
            mates: names(mates),
        },
    );
    scope
}

fn record_pass(rows: &[String], scope: &ReferenceRecordScope) -> Pass {
    pass_over(
        rows,
        Collector::Release116(scope, &CacheContigNames::Human),
        &row_pair_rules(),
    )
}

fn one_sided(count: u64, excluded: u64, records: u64, perl: u64, rust: u64) -> BucketStats {
    vec![
        ("count", count),
        ("excluded", excluded),
        ("records", records),
        ("perl", perl),
        ("rust", rust),
    ]
}

#[test]
fn surplus_transcripts_of_the_scored_engine_on_a_giant_record_are_set_aside() {
    let scope = selection_scope(
        GIANT,
        &["T1", "T2", "T3", "T4"],
        &["T1", "T2", "T3", "T4", "T5", "T6", "T7"],
    );
    let rows = ["T5", "T6", "T7"].map(|feature| ablation_row("extra_in_rust", GIANT, feature));
    let pass = selection_pass(&rows, &scope);
    assert_eq!(pass.excluded, (0, 3));
    assert_eq!(
        pass.class(SV_TRANSCRIPT_SELECTION_BUCKET),
        Some(&selection(3, 3, 0, 3))
    );
    assert!(pass.open_rows.is_empty(), "{:?}", pass.open_rows);
}

#[test]
fn a_pair_on_a_transcript_both_engines_named_stays_charged() {
    let scope = selection_scope(GIANT, &["T1", "T2"], &["T1", "T2", "T3"]);
    let rows = [
        deletion_row("missing_in_rust", GIANT, "T1", "intron_variant"),
        deletion_row("extra_in_rust", GIANT, "T1", "intron_variant,stop_lost"),
        ablation_row("extra_in_rust", GIANT, "T3"),
    ];
    let pass = selection_pass(&rows, &scope);
    assert_eq!(pass.excluded, (0, 1));
    assert_eq!(pass.open_rows.len(), 2);
    assert_eq!(
        pass.class(SV_TRANSCRIPT_SELECTION_BUCKET),
        Some(&selection(3, 1, 0, 1)),
        "all three rows sit on the record in scope; one is the class"
    );
}

#[test]
fn a_record_the_scored_engine_named_fewer_transcripts_for_stays_charged() {
    let scope = selection_scope(GIANT, &["T1", "T2", "T3", "T4", "T5"], &["T1", "T2"]);
    let rows = ["T3", "T4", "T5"].map(|feature| ablation_row("missing_in_rust", GIANT, feature));
    let pass = selection_pass(&rows, &scope);
    assert_eq!(pass.excluded, (0, 0));
    assert_eq!(pass.open_rows.len(), 3);
    assert_eq!(
        pass.class(SV_TRANSCRIPT_SELECTION_BUCKET),
        Some(&selection(3, 0, 0, 0))
    );
}

#[test]
fn a_reference_only_transcript_is_set_aside_when_the_scored_engine_named_as_many() {
    let scope = selection_scope(GIANT, &["T1", "T2", "T9"], &["T1", "T2", "T3"]);
    let surplus = ablation_row("extra_in_rust", GIANT, "T3");
    let rows = [
        ablation_row("missing_in_rust", GIANT, "T9"),
        surplus.clone(),
    ];
    let pass = selection_pass(&rows, &scope);
    assert_eq!(
        pass.excluded,
        (1, 0),
        "three transcripts against three: the reference-only one is the class, and the scored \
         engine did not name more, so its own stays"
    );
    assert_eq!(pass.open_rows, [surplus.trim_end()]);
}

#[test]
fn a_record_the_reference_named_no_transcript_for_is_outside_the_class() {
    let scope = selection_scope(GIANT, &[], &["T1", "T2"]);
    let rows = ["T1", "T2"].map(|feature| ablation_row("extra_in_rust", GIANT, feature));
    let pass = selection_pass(&rows, &scope);
    assert_eq!(pass.excluded, (0, 0));
    assert_eq!(pass.open_rows.len(), 2);
}

#[test]
fn a_record_below_the_size_limit_never_enters_the_class() {
    let rows = ["T5", "T6"].map(|feature| ablation_row("extra_in_rust", SMALL, feature));
    let pass = selection_pass(&rows, &ScopeFeatures::default());
    assert_eq!(pass.excluded, (0, 0));
    assert_eq!(pass.class(SV_TRANSCRIPT_SELECTION_BUCKET), None);
    assert_eq!(pass.open_rows.len(), 2);
}

#[test]
fn the_115_2_pass_reports_the_one_class_of_transcript_selection_and_sets_its_rows_aside() {
    assert_eq!(SV_TRANSCRIPT_SELECTION_BUCKET, "sv_transcript_selection");
    let scope = selection_scope(GIANT, &["T1"], &["T1", "T2"]);
    let pass = selection_pass(&[ablation_row("extra_in_rust", GIANT, "T2")], &scope);
    let one_sided_classes: Vec<&str> = pass
        .buckets
        .iter()
        .map(|(bucket, _)| bucket.as_str())
        .filter(|bucket| EXCLUSION_REGISTRY.iter().all(|rule| rule.bucket != *bucket))
        .collect();
    assert_eq!(one_sided_classes, [SV_TRANSCRIPT_SELECTION_BUCKET]);
    assert_eq!(pass.excluded, (0, 1));
}

#[test]
fn both_releases_put_a_pair_to_the_same_pair_rules() {
    let pair = |location: &str, allele: &str, feature: &str, perl: &str, rust: &str| {
        [
            row("missing_in_rust", location, allele, feature, perl),
            row("extra_in_rust", location, allele, feature, rust),
        ]
    };
    let rows = [
        pair(
            "21:100",
            "A",
            "ENST1",
            "missense_variant",
            "missense_variant,splice_region_variant",
        ),
        pair(
            "21:300",
            "T",
            "ENST3",
            "intron_variant",
            "intron_variant,splice_polypyrimidine_tract_variant",
        ),
        pair(
            "21:500",
            "G",
            "ENST5",
            "5_prime_UTR_variant,start_lost,start_retained_variant",
            "5_prime_UTR_variant,start_retained_variant",
        ),
    ]
    .concat();
    let pair_classes = [
        ("intronic_overcall_swap", 1, 0),
        ("ppt_overcall_swap", 0, 0),
        ("frameshift_splice_cascade_swap", 0, 0),
        ("splice_lastwrite_swap", 1, 1),
        ("start_cooccurrence_swap", 1, 1),
    ]
    .map(|(bucket, count, excluded)| {
        (
            bucket.to_string(),
            vec![("count", count), ("excluded", excluded)],
        )
    });
    let earlier = selection_pass(&rows, &ScopeFeatures::default());
    assert_eq!(earlier.buckets, pair_classes);
    assert_eq!(earlier.excluded, (2, 2));
    let later = record_pass(&rows, &ReferenceRecordScope::default());
    assert_eq!(later.buckets[..pair_classes.len()], pair_classes);
    assert_eq!(later.excluded, (2, 2));
    assert_eq!(later.open_rows, earlier.open_rows);
}

#[test]
fn every_scored_row_of_a_skipped_record_is_set_aside() {
    let rows = [
        row("extra_in_rust", "21:800-900", "CPX", "C3", "intron_variant"),
        row(
            "extra_in_rust",
            "21:800-900",
            "CPX",
            "C4",
            "coding_sequence_variant",
        ),
        row(
            "extra_in_rust",
            "21:30000001",
            "deletion",
            "TX_DEL",
            "feature_truncation",
        ),
        row(
            "extra_in_rust",
            "21:600-700",
            "<NON_REF>",
            "-",
            "intergenic_variant",
        ),
        row(
            "missing_in_rust",
            "21:400",
            "G",
            "TX_DOWN",
            "downstream_gene_variant",
        ),
    ];
    let mut scope = ReferenceRecordScope::default();
    scope.absent = groups(&["21:800-900", "21:30000001", "21:600-700"]);
    let pass = record_pass(&rows, &scope);
    assert_eq!(pass.excluded, (0, 4));
    assert_eq!(
        pass.class("reference_skipped_record"),
        Some(&one_sided(4, 4, 3, 0, 4))
    );
    assert_eq!(pass.open_rows.len(), 1, "the reference-only row stays");
}

#[test]
fn a_record_the_reference_kept_is_not_set_aside_under_the_116_2_reference() {
    let rows = [row(
        "extra_in_rust",
        "21:18457513-32484691",
        "N.",
        "T3",
        "transcript_ablation",
    )];
    let pass = record_pass(&rows, &ReferenceRecordScope::default());
    assert_eq!(pass.excluded, (0, 0));
    assert_eq!(pass.open_rows.len(), 1);
    assert_eq!(
        pass.class("reference_skipped_record"),
        Some(&one_sided(0, 0, 0, 0, 0))
    );
    assert_eq!(pass.class(SV_TRANSCRIPT_SELECTION_BUCKET), None);
}

#[test]
fn own_end_rows_the_reference_lost_are_counted_and_stay_charged() {
    let location = "21:26037668";
    let scope = breakend_scope(
        location,
        &["N[22:26037668["],
        &["N[22:26037668[", "chromosome_breakpoint"],
        &["22"],
    );
    let rows = [
        row(
            "extra_in_rust",
            location,
            "chromosome_breakpoint",
            "O1",
            "upstream_gene_variant",
        ),
        row(
            "extra_in_rust",
            location,
            "chromosome_breakpoint",
            "O2",
            "intron_variant",
        ),
    ];
    let pass = record_pass(&rows, &scope);
    assert_eq!(pass.excluded, (0, 0));
    assert_eq!(
        pass.class("bnd_own_chromosome_rows_lost"),
        Some(&one_sided(2, 0, 1, 0, 0))
    );
    assert_eq!(pass.open_rows.len(), 2);
}

#[test]
fn own_end_rows_beside_a_reference_own_end_row_are_not_the_lost_row_shape() {
    let location = "21:26037668";
    let alleles = ["N[22:26037668[", "N."];
    let scope = breakend_scope(location, &alleles, &alleles, &["22"]);
    let rows = [row("extra_in_rust", location, "N.", "O2", "intron_variant")];
    let pass = record_pass(&rows, &scope);
    assert_eq!(pass.count("bnd_own_chromosome_rows_lost"), 0);
}

#[test]
fn a_mate_on_the_records_own_chromosome_is_not_the_lost_row_shape() {
    let location = "21:24808769-35740566";
    let scope = breakend_scope(
        location,
        &["N[21:35740567["],
        &["N[21:35740567[", "N."],
        &["21"],
    );
    let rows = [row("extra_in_rust", location, "N.", "O2", "intron_variant")];
    let pass = record_pass(&rows, &scope);
    assert_eq!(pass.count("bnd_own_chromosome_rows_lost"), 0);
}

#[test]
fn mate_rows_under_a_mate_named_by_an_accession_are_counted() {
    let location = "21:26037668";
    let mate = "N[NC_000022.11:26037668[";
    let scope = breakend_scope(location, &["N."], &["N.", mate], &["NC_000022.11"]);
    let rows = [row("extra_in_rust", location, mate, "M1", "intron_variant")];
    let pass = record_pass(&rows, &scope);
    assert_eq!(
        pass.class("bnd_synonym_mate_rows_lost"),
        Some(&one_sided(1, 0, 1, 0, 0))
    );
    assert_eq!(pass.open_rows.len(), 1);
}

#[test]
fn a_mate_named_as_the_cache_names_it_is_not_the_synonym_shape() {
    let location = "21:26037668";
    let mate = "N[chr22:26037668[";
    let scope = breakend_scope(location, &["N."], &["N.", mate], &["22"]);
    let rows = [row("extra_in_rust", location, mate, "M1", "intron_variant")];
    let pass = record_pass(&rows, &scope);
    assert_eq!(pass.count("bnd_synonym_mate_rows_lost"), 0);
}

#[test]
fn a_mate_is_told_from_a_synonym_by_the_names_the_run_gives_and_by_none_without_them() {
    let location = "II:500";
    let mate = "N[Mito:700[";
    let scope = breakend_scope(location, &["N."], &["N.", mate], &["Mito"]);
    let rows = [row("extra_in_rust", location, mate, "T1", "intron_variant")];
    let class = |cache_names: &CacheContigNames| {
        let collector = Collector::Release116(&scope, cache_names);
        let pass = pass_over(&rows, collector, &row_pair_rules());
        assert_eq!(pass.open_rows.len(), 1);
        pass.class("bnd_synonym_mate_rows_lost").cloned()
    };
    let not_evaluated = |rows: u64| {
        let mut statistics = one_sided(0, 0, 0, 0, 0);
        statistics.push((NOT_EVALUATED, rows));
        statistics
    };
    assert_eq!(
        class(&CacheContigNames::Named(names(&["1", "II", "Mito"]))),
        Some(one_sided(0, 0, 0, 0, 0))
    );
    assert_eq!(
        class(&CacheContigNames::Named(names(&["1", "II", "MT"]))),
        Some(one_sided(1, 0, 1, 0, 0))
    );
    assert_eq!(
        class(&CacheContigNames::Human),
        Some(one_sided(1, 0, 1, 0, 0))
    );
    assert_eq!(class(&CacheContigNames::Unknown), Some(not_evaluated(1)));
}

#[test]
fn a_mitochondrial_mate_written_m_is_a_synonym_of_the_name_the_cache_files_it_under() {
    let location = "21:26037668";
    let counted = |mate: &str, chromosome: &str| {
        let scope = breakend_scope(location, &["N."], &["N.", mate], &[chromosome]);
        let rows = [row(
            "extra_in_rust",
            location,
            mate,
            "MT_TX",
            "upstream_gene_variant",
        )];
        record_pass(&rows, &scope).count("bnd_synonym_mate_rows_lost")
    };
    assert_eq!(counted("N[chrM:100[", "M"), 1);
    assert_eq!(counted("N[MT:100[", "MT"), 0);
}

#[test]
fn the_reference_rows_of_a_record_whose_mate_overflows_the_buffer_are_counted() {
    let rows = [
        row(
            "missing_in_rust",
            "21:30425651",
            "N.",
            "-",
            "intergenic_variant",
        ),
        row(
            "missing_in_rust",
            "21:30425651",
            "N[21:2147478648[",
            "-",
            "intergenic_variant",
        ),
        row(
            "extra_in_rust",
            "21:30425651",
            "N.",
            "T2",
            "feature_truncation,coding_sequence_variant",
        ),
        row(
            "missing_in_rust",
            "21:30425651-30425899",
            "duplication",
            "-",
            "intergenic_variant",
        ),
    ];
    let mut scope = ReferenceRecordScope::default();
    scope.mate_overflow = groups(&["21:30425651"]);
    let pass = record_pass(&rows, &scope);
    assert_eq!(pass.excluded, (0, 0));
    assert_eq!(pass.open_rows.len(), 4);
    assert_eq!(
        pass.class("bnd_mate_32bit_buffer_intergenic"),
        Some(&one_sided(2, 0, 1, 0, 0))
    );
}

#[test]
fn a_dot_breakend_pair_and_a_reference_row_at_a_coordinate_that_is_no_integer_are_set_aside() {
    let span = "21:30425650-30425651";
    let point = "21:30425651";
    let breakpoint = "chromosome_breakpoint";
    let rows = [
        row(
            "missing_in_rust",
            span,
            breakpoint,
            "T2",
            "coding_sequence_variant",
        ),
        row(
            "extra_in_rust",
            point,
            breakpoint,
            "T2",
            "coding_sequence_variant",
        ),
        row(
            "missing_in_rust",
            span,
            breakpoint,
            "T1",
            "upstream_gene_variant",
        ),
        row(
            "missing_in_rust",
            "21:30425651-30425949.5",
            "deletion",
            "T2",
            "feature_truncation,coding_sequence_variant",
        ),
        row(
            "missing_in_rust",
            "21:30425651-3.0425949e7",
            "deletion",
            "T1",
            "upstream_gene_variant",
        ),
        row("extra_in_rust", point, breakpoint, "T3", "intron_variant"),
    ];
    let pass = record_pass(&rows, &ReferenceRecordScope::default());
    assert_eq!(pass.excluded, (3, 1));
    assert_eq!(
        pass.class("dot_breakend_without_svtype_pair"),
        Some(&one_sided(2, 2, 2, 1, 1))
    );
    assert_eq!(
        pass.class("non_integer_end_annotated"),
        Some(&one_sided(2, 2, 2, 2, 0))
    );
    let mut open_features = pass.open_column(4);
    open_features.sort_unstable();
    assert_eq!(
        open_features,
        ["T1", "T3"],
        "the reference's row on a feature the scored engine did not name stays charged"
    );

    let two_sets_on_one_feature = [
        row(
            "missing_in_rust",
            span,
            breakpoint,
            "T2",
            "coding_sequence_variant",
        ),
        row(
            "extra_in_rust",
            point,
            breakpoint,
            "T2",
            "coding_sequence_variant",
        ),
        row(
            "extra_in_rust",
            point,
            breakpoint,
            "T2",
            "coding_sequence_variant,feature_truncation",
        ),
    ];
    let pass = record_pass(&two_sets_on_one_feature, &ReferenceRecordScope::default());
    assert_eq!(pass.excluded, (1, 2));
    assert!(pass.open_rows.is_empty(), "{:?}", pass.open_rows);
    assert_eq!(
        pass.class("dot_breakend_without_svtype_pair"),
        Some(&one_sided(3, 3, 2, 1, 2)),
        "the reference's span row is counted once"
    );
}

#[test]
fn a_reference_intergenic_row_against_scored_transcript_rows_is_counted() {
    let rows = [
        row("missing_in_rust", "21:100", "A", "-", "intergenic_variant"),
        row(
            "extra_in_rust",
            "21:100",
            "A",
            "TX_FIRST",
            "missense_variant",
        ),
        row(
            "extra_in_rust",
            "21:100",
            "A",
            "TX_SECOND",
            "intron_variant",
        ),
        row("missing_in_rust", "21:200", "T", "-", "intergenic_variant"),
        row("extra_in_rust", "21:200", "T", "-", "intergenic_variant"),
    ];
    let pass = record_pass(&rows, &ReferenceRecordScope::default());
    assert_eq!(pass.excluded, (0, 0));
    assert_eq!(
        pass.class("fasta_named_slice_intergenic"),
        Some(&one_sided(1, 0, 1, 0, 0))
    );
    assert_eq!(pass.open_rows.len(), 5);
}

#[test]
fn a_reference_only_location_with_a_negative_first_coordinate_is_the_numified_pos_class() {
    let locations = [
        "21:-4-295",
        "21:-2147483647--2147483348",
        "21:-5",
        "21:1-30425949",
    ];
    let rows = [
        row(
            "missing_in_rust",
            locations[0],
            "deletion",
            "-",
            "intergenic_variant",
        ),
        row(
            "missing_in_rust",
            locations[1],
            "deletion",
            "-",
            "intergenic_variant",
        ),
        row(
            "missing_in_rust",
            locations[2],
            "G",
            "-",
            "intergenic_variant",
        ),
        row(
            "missing_in_rust",
            locations[3],
            "deletion",
            "ENST00000262354",
            "transcript_ablation",
        ),
    ];
    let mut scope = ReferenceRecordScope::default();
    scope.reference_only = groups(&locations);
    let pass = record_pass(&rows, &scope);
    assert_eq!(pass.excluded, (3, 0));
    assert_eq!(
        pass.class("non_numeric_pos_numified"),
        Some(&one_sided(3, 3, 3, 3, 0))
    );
    assert_eq!(
        pass.count("non_integer_end_annotated"),
        0,
        "a negative coordinate is an integer"
    );
    assert_eq!(
        pass.open_column(2),
        ["21:1-30425949"],
        "a start of 1 is not told from a record at position 0"
    );
}

#[test]
fn the_dot_breakend_pair_covers_the_own_and_the_mate_allele_of_a_mate_named_in_info() {
    let span = "21:30425650-30425651";
    let point = "21:30425651";
    let mut rows = Vec::new();
    for allele in ["N.", "N[21:30425700["] {
        for (feature, consequence) in [
            ("ENST00000334068", "upstream_gene_variant"),
            (
                "ENST00000390690",
                "feature_truncation,coding_sequence_variant",
            ),
        ] {
            rows.push(row("missing_in_rust", span, allele, feature, consequence));
            rows.push(row("extra_in_rust", point, allele, feature, consequence));
        }
    }
    let sequence_allele = "ACG";
    for (source, location) in [("missing_in_rust", span), ("extra_in_rust", point)] {
        rows.push(row(
            source,
            location,
            sequence_allele,
            "ENST00000390690",
            "inframe_insertion",
        ));
    }
    let pass = record_pass(&rows, &ReferenceRecordScope::default());
    assert_eq!(pass.excluded, (4, 4));
    assert_eq!(
        pass.class("dot_breakend_without_svtype_pair"),
        Some(&one_sided(8, 8, 2, 4, 4))
    );
    assert_eq!(
        pass.open_column(3),
        [sequence_allele, sequence_allele],
        "a point under a sequence allele pairs with nothing"
    );
}

#[test]
fn the_dot_breakend_pair_takes_the_point_at_the_spans_start_or_the_base_after_it() {
    let point = "21:30425651";
    let breakpoint = "chromosome_breakpoint";
    let features = [
        ("ENST00000334068", "upstream_gene_variant"),
        ("ENST00000390690", "coding_sequence_variant"),
    ];
    for span in [
        "21:30425651-30425652",
        "21:30425651-30425653",
        "21:30425650-30425652",
        "21:30425650-30425653",
    ] {
        let rows: Vec<String> = [("missing_in_rust", span), ("extra_in_rust", point)]
            .iter()
            .flat_map(|(source, location)| {
                features.map(|(feature, consequence)| {
                    row(source, location, breakpoint, feature, consequence)
                })
            })
            .collect();
        let pass = record_pass(&rows, &ReferenceRecordScope::default());
        assert_eq!(pass.excluded, (2, 2), "{span}");
        assert!(pass.open_rows.is_empty(), "{span}: {:?}", pass.open_rows);
        assert_eq!(
            pass.class("dot_breakend_without_svtype_pair"),
            Some(&one_sided(4, 4, 2, 2, 2)),
            "{span}"
        );
    }
}

#[test]
fn a_point_two_bases_past_a_spans_start_is_no_dot_breakend_pair() {
    let rows = [
        row(
            "missing_in_rust",
            "21:30425649-30425651",
            "chromosome_breakpoint",
            "ENST00000390690",
            "coding_sequence_variant",
        ),
        row(
            "extra_in_rust",
            "21:30425651",
            "chromosome_breakpoint",
            "ENST00000390690",
            "coding_sequence_variant",
        ),
    ];
    let pass = record_pass(&rows, &ReferenceRecordScope::default());
    assert_eq!(pass.excluded, (0, 0));
    assert_eq!(pass.open_rows.len(), 2);
    assert_eq!(pass.count("dot_breakend_without_svtype_pair"), 0);
}

#[test]
fn a_mate_prefix_the_reference_does_not_resolve_is_the_synonym_shape() {
    let location = "21:30425651";
    let unresolved = "N[Chr21:30425700[";
    let scope = breakend_scope(location, &["N."], &["N.", unresolved], &["21"]);
    let rows = [
        row(
            "extra_in_rust",
            location,
            unresolved,
            "ENST00000334068",
            "upstream_gene_variant",
        ),
        row(
            "extra_in_rust",
            location,
            unresolved,
            "ENST00000390690",
            "feature_truncation,coding_sequence_variant",
        ),
    ];
    let pass = record_pass(&rows, &scope);
    assert_eq!(
        pass.class("bnd_synonym_mate_rows_lost"),
        Some(&one_sided(2, 0, 1, 0, 0))
    );
    assert_eq!(pass.open_rows.len(), 2);

    let resolved = "N[chr21:30425700[";
    let alleles = ["N.", resolved];
    let scope = breakend_scope(location, &alleles, &alleles, &["21"]);
    let rows = [row(
        "extra_in_rust",
        location,
        resolved,
        "ENST00000418755",
        "downstream_gene_variant",
    )];
    let pass = record_pass(&rows, &scope);
    assert_eq!(pass.count("bnd_synonym_mate_rows_lost"), 0);
}
