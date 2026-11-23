// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The per-consequence-class counts on values written for each test: the two derivations of
//! one term and their guards, the per-term counts of a key file and of a discordance report,
//! their assembly and pooling, the floor and its threshold, the two CSVs, the totals cache and
//! the adjusted F1.

use super::*;
use crate::input::CANONICAL_CONTIGS;
use crate::sort::{KeySorter, DEFAULT_MEMORY_BOUND};
use crate::TempDir;
use std::path::PathBuf;
use std::time::{Duration, UNIX_EPOCH};

fn counts(term: &str, vep_tuples: u64, vep_only: u64, engine_only: u64) -> PerClassCounts {
    PerClassCounts::new(term, vep_tuples, vep_only, engine_only).unwrap()
}

fn counter(entries: &[(&str, u64)]) -> Counter {
    let mut counter = Counter::default();
    for (term, n) in entries {
        counter.add(term, *n);
    }
    counter
}

/// The terms of a counter with their counts, by term.
fn counted(counter: &Counter) -> Vec<(&str, u64)> {
    let mut entries: Vec<(&str, u64)> = counter
        .entries
        .iter()
        .map(|(term, n)| (term.as_str(), *n))
        .collect();
    entries.sort_unstable();
    entries
}

fn terms(rows: &[PerClassCounts]) -> Vec<&str> {
    rows.iter().map(|row| row.term.as_str()).collect()
}

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "classes_test_").unwrap()
}

/// Writes `content` to `name` under `scratch` and returns the path.
fn written(scratch: &TempDir, name: &str, content: &str) -> PathBuf {
    let path = scratch.path.join(name);
    std::fs::write(&path, content).unwrap();
    path
}

/// The message of the refusal `result` holds.
fn refusal<T>(result: io::Result<T>) -> String {
    match result {
        Ok(_) => panic!("accepted"),
        Err(error) => error.to_string(),
    }
}

/// A CSV file as its rows, each a map from column name to cell.
fn csv_rows(path: &Path) -> Vec<HashMap<String, String>> {
    let (header, rows) = csv::read_table(path, b',').unwrap();
    rows.into_iter()
        .map(|row| header.iter().cloned().zip(row).collect())
        .collect()
}

#[test]
fn the_intersection_is_the_vep_total_less_the_vep_only_tuples() {
    assert_eq!(counts("missense_variant", 100, 7, 3).intersection(), 93);
}

#[test]
fn the_engine_total_is_the_intersection_plus_the_engine_only_tuples() {
    assert_eq!(counts("missense_variant", 100, 7, 3).engine_tuples(), 96);
}

#[test]
fn precision_and_recall_divide_by_their_own_totals() {
    let c = counts("t", 100, 7, 3);
    assert_eq!(c.precision(), 93.0 / 96.0);
    assert_eq!(c.recall(), 93.0 / 100.0);
}

#[test]
fn f1_is_twice_the_intersection_over_the_two_totals_to_the_last_bit() {
    let c = counts("t", 100, 7, 3);
    assert_eq!(c.f1(), (2 * 93) as f64 / (100 + 96) as f64);
}

#[test]
fn a_term_without_a_discordant_tuple_reads_one() {
    assert_eq!(counts("t", 50, 0, 0).f1_str(), "1.000000");
}

#[test]
fn a_term_the_reference_never_writes_reads_zero() {
    let c = counts("invented_term", 0, 0, 4);
    assert_eq!(c.intersection(), 0);
    assert_eq!(c.engine_tuples(), 4);
    assert_eq!(c.recall(), 0.0);
    assert_eq!(c.f1_str(), "0.000000");
}

#[test]
fn a_term_neither_side_writes_reads_zero_and_divides_nothing() {
    let c = counts("t", 0, 0, 0);
    assert_eq!(c.f1(), 0.0);
    assert_eq!(c.precision(), 0.0);
}

#[test]
fn more_vep_only_tuples_than_vep_tuples_are_refused_by_term() {
    let message = PerClassCounts::new("synonymous_variant", 1, 5, 0).unwrap_err();
    assert!(message.contains("synonymous_variant"), "{message}");
    assert!(message.contains("5 VEP-only"), "{message}");
}

const BY_DATASET_HEADER: &str =
    "dataset,assembly,engine,sweep_id,binary_md5,term,vep_tuples,vep_only,engine_only\n";

#[test]
fn a_negative_count_is_refused_where_a_per_class_csv_is_read() {
    let scratch = scratch();
    for row in [
        "s01,GRCh37,vep-rs,s,b,t,10,-1,0",
        "s01,GRCh37,vep-rs,s,b,t,10,0,-2",
    ] {
        let path = written(
            &scratch,
            "negative.csv",
            &format!("{BY_DATASET_HEADER}{row}\n"),
        );
        assert!(read_per_class_by_dataset_csv(&path).is_err(), "{row}");
    }
}

#[test]
fn a_term_whose_every_tuple_is_vep_only_is_taken_and_reads_zero() {
    let c = counts("t", 9, 9, 0);
    assert_eq!(c.intersection(), 0);
    assert_eq!(c.f1_str(), "0.000000");
}

#[test]
fn a_tuple_of_several_terms_counts_for_each_of_them() {
    let scratch = scratch();
    let path = written(
        &scratch,
        "keys.sorted.txt",
        "21:100\tA\tENST1\tTranscript\tintron_variant,splice_region_variant\n",
    );
    assert_eq!(
        counted(&count_terms_in_sorted_keys(&path).unwrap()),
        [("intron_variant", 1), ("splice_region_variant", 1)]
    );
}

#[test]
fn a_term_on_two_tuples_counts_twice() {
    let scratch = scratch();
    let path = written(
        &scratch,
        "keys.sorted.txt",
        "21:100\tA\tENST1\tTranscript\tmissense_variant\n\
         21:200\tC\tENST1\tTranscript\tmissense_variant\n",
    );
    let totals = count_terms_in_sorted_keys(&path).unwrap();
    assert_eq!(totals.get("missense_variant"), 2);
}

#[test]
fn a_blank_line_of_a_key_file_counts_nothing() {
    let scratch = scratch();
    let path = written(
        &scratch,
        "keys.txt",
        "21:100\tA\tENST1\tTranscript\tmissense_variant\n\n",
    );
    let totals = count_terms_in_sorted_keys(&path).unwrap();
    assert_eq!(totals.get("missense_variant"), 1);
}

#[test]
fn a_key_line_of_four_fields_is_refused_and_not_read_as_a_term() {
    let scratch = scratch();
    let path = written(
        &scratch,
        "keys.sorted.txt",
        "21:100\tA\tENST1\tmissense_variant\n",
    );
    let message = refusal(count_terms_in_sorted_keys(&path));
    assert!(
        message.ends_with("expected 5 tab-separated key fields, got 4"),
        "{message}"
    );
}

#[test]
fn an_empty_key_file_counts_no_term() {
    let scratch = scratch();
    let path = written(&scratch, "empty.txt", "");
    let totals = count_terms_in_sorted_keys(&path).unwrap();
    assert!(totals.entries.is_empty());
}

#[test]
fn a_tuple_written_twice_counts_once_through_the_sort_the_comparison_scores() {
    let scratch = scratch();
    let keys = [
        "21:100\tA\tENST1\tTranscript\tmissense_variant",
        "21:100\tA\tENST1\tTranscript\tmissense_variant",
        "21:200\tC\tENST1\tTranscript\tmissense_variant",
    ];
    let sorted = scratch.path.join("sorted.txt");
    let mut sorter = KeySorter::new(&scratch.path, "sorted.txt", DEFAULT_MEMORY_BOUND);
    for key in keys {
        sorter.push(key.as_bytes()).unwrap();
    }
    sorter.finish(&sorted).unwrap();
    let totals = count_terms_in_sorted_keys(&sorted).unwrap();
    assert_eq!(totals.get("missense_variant"), 2);
}

const DISCORDANT_COLUMNS: &str =
    "file_name\tsource\tlocation\tallele\tfeature\tfeature_type\tconsequence_set";

/// A discordance report of one row a (source, consequence set), its lines ended as the
/// `csv` module's writer ends them.
fn discordance(scratch: &TempDir, rows: &[(&str, &str)]) -> PathBuf {
    let mut report = format!("{DISCORDANT_COLUMNS}\r\n");
    for (i, (source, csq)) in rows.iter().enumerate() {
        report.push_str(&format!(
            "f.txt\t{source}\t21:{i}\tA\tENST1\tTranscript\t{csq}\r\n"
        ));
    }
    written(scratch, "discordant.tsv", &report)
}

#[test]
fn a_discordance_report_is_counted_a_side_by_its_source_column() {
    let scratch = scratch();
    let path = discordance(
        &scratch,
        &[
            ("missing_in_rust", "intron_variant"),
            ("extra_in_rust", "intron_variant,splice_region_variant"),
        ],
    );
    let (vep_only, engine_only) = count_terms_in_discordance(&path).unwrap();
    assert_eq!(counted(&vep_only), [("intron_variant", 1)]);
    assert_eq!(
        counted(&engine_only),
        [("intron_variant", 1), ("splice_region_variant", 1)]
    );
}

#[test]
fn a_source_that_names_neither_side_is_refused_by_name() {
    let scratch = scratch();
    let path = discordance(&scratch, &[("both_ways", "intron_variant")]);
    let message = refusal(count_terms_in_discordance(&path));
    assert!(message.contains("both_ways"), "{message}");
}

#[test]
fn a_discordance_report_without_its_consequence_column_is_refused_by_column() {
    let scratch = scratch();
    let path = written(&scratch, "d.tsv", "file_name\tsource\n");
    let message = refusal(count_terms_in_discordance(&path));
    assert!(message.contains("consequence_set"), "{message}");
}

#[test]
fn the_terms_of_a_report_row_count_the_same_in_any_order() {
    let scratch = scratch();
    let path = discordance(
        &scratch,
        &[("extra_in_rust", "splice_region_variant,intron_variant")],
    );
    let (_, engine_only) = count_terms_in_discordance(&path).unwrap();
    assert_eq!(engine_only.get("intron_variant"), 1);
    assert_eq!(engine_only.get("splice_region_variant"), 1);
}

#[test]
fn a_discordance_report_of_a_header_alone_counts_nothing_on_either_side() {
    let scratch = scratch();
    let path = discordance(&scratch, &[]);
    let (vep_only, engine_only) = count_terms_in_discordance(&path).unwrap();
    assert!(vep_only.entries.is_empty());
    assert!(engine_only.entries.is_empty());
}

#[test]
fn a_term_on_any_of_the_three_counts_has_a_row() {
    let rows = derive_per_class_counts(
        &counter(&[("a", 10), ("b", 5)]),
        &counter(&[("a", 1)]),
        &counter(&[("c", 2)]),
    )
    .unwrap();
    let mut seen = terms(&rows);
    seen.sort_unstable();
    assert_eq!(seen, ["a", "b", "c"]);
}

#[test]
fn a_term_the_engine_alone_writes_keeps_its_row_with_no_vep_tuples() {
    let rows =
        derive_per_class_counts(&counter(&[("a", 10)]), &counter(&[]), &counter(&[("c", 2)]))
            .unwrap();
    let c = rows.iter().find(|row| row.term == "c").expect("the row");
    assert_eq!(c.vep_tuples, 0);
    assert_eq!(c.engine_only, 2);
}

#[test]
fn the_rows_descend_by_vep_tuples_and_ascend_by_term_within_a_tie() {
    let rows = derive_per_class_counts(
        &counter(&[("low", 1), ("high", 99), ("tie_b", 5), ("tie_a", 5)]),
        &counter(&[]),
        &counter(&[]),
    )
    .unwrap();
    assert_eq!(terms(&rows), ["high", "tie_a", "tie_b", "low"]);
}

#[test]
fn three_empty_counts_give_no_row() {
    let rows = derive_per_class_counts(&counter(&[]), &counter(&[]), &counter(&[])).unwrap();
    assert!(rows.is_empty());
}

#[test]
fn pooling_sums_each_of_the_three_counts_of_a_term() {
    let pooled = pool_per_class_counts(&[counts("a", 100, 2, 1), counts("a", 50, 1, 4)]).unwrap();
    assert_eq!(pooled, [counts("a", 150, 3, 5)]);
}

#[test]
fn the_pooled_derived_counts_are_the_sums_of_the_derived_counts() {
    let parts = [
        counts("a", 100, 2, 1),
        counts("a", 50, 1, 4),
        counts("a", 7, 0, 0),
    ];
    let pooled = &pool_per_class_counts(&parts).unwrap()[0];
    let summed = |derived: fn(&PerClassCounts) -> u64| parts.iter().map(derived).sum::<u64>();
    assert_eq!(pooled.intersection(), summed(PerClassCounts::intersection));
    assert_eq!(
        pooled.engine_tuples(),
        summed(PerClassCounts::engine_tuples)
    );
}

#[test]
fn the_pooled_f1_weighs_each_dataset_by_its_tuples() {
    let pooled =
        pool_per_class_counts(&[counts("a", 1_000_000, 0, 0), counts("a", 2, 2, 0)]).unwrap();
    assert_eq!(pooled[0].f1_str(), "0.999999");
}

#[test]
fn the_pooled_rows_descend_by_vep_tuples() {
    let pooled =
        pool_per_class_counts(&[counts("small", 1, 0, 0), counts("big", 9, 0, 0)]).unwrap();
    assert_eq!(terms(&pooled), ["big", "small"]);
}

#[test]
fn the_floor_is_the_lowest_f1_among_the_terms_at_or_above_the_threshold() {
    let rows = [
        counts("good", 10_000, 1, 0),
        counts("worst", 10_000, 100, 0),
    ];
    let floor = per_class_floor(&rows, 1000);
    assert_eq!(floor.term.as_deref(), Some("worst"));
    assert_eq!(floor.n_eligible, 2);
    assert_eq!(floor.n_below_threshold, 0);
}

#[test]
fn a_term_below_the_threshold_is_counted_and_is_never_the_floor() {
    let rows = [counts("solid", 10_000, 1, 0), counts("rare", 3, 1, 0)];
    let floor = per_class_floor(&rows, 1000);
    assert_eq!(floor.term.as_deref(), Some("solid"));
    assert_eq!(floor.n_below_threshold, 1);
}

#[test]
fn rows_all_below_the_threshold_give_no_floor() {
    let floor = per_class_floor(&[counts("rare", 3, 1, 0)], 1000);
    assert_eq!(floor.f1, None);
    assert_eq!(floor.term, None);
    assert_eq!(floor.n_eligible, 0);
    assert_eq!(floor.n_below_threshold, 1);
}

#[test]
fn terms_at_one_f1_give_the_floor_to_the_first_by_name() {
    let rows = [counts("zeta", 10_000, 5, 0), counts("alpha", 10_000, 5, 0)];
    assert_eq!(per_class_floor(&rows, 1000).term.as_deref(), Some("alpha"));
}

#[test]
fn a_term_at_the_threshold_itself_is_eligible() {
    let floor = per_class_floor(&[counts("edge", 1000, 1, 0)], 1000);
    assert_eq!(floor.n_eligible, 1);
}

#[test]
fn the_threshold_is_one_thousand_tuples() {
    assert_eq!(PER_CLASS_MIN_TUPLES, 1000);
    let floor = per_class_floor(&[counts("t", 999, 1, 0)], PER_CLASS_MIN_TUPLES);
    assert_eq!(floor.n_eligible, 0);
}

#[test]
fn no_row_gives_no_floor() {
    let floor = per_class_floor(&[], PER_CLASS_MIN_TUPLES);
    assert_eq!(floor.f1, None);
    assert_eq!(floor.n_below_threshold, 0);
}

const SWEEP: &str = "sweep-20260101T000000Z";

fn two_terms() -> Vec<PerClassCounts> {
    vec![
        counts("missense_variant", 100, 3, 1),
        counts("intron_variant", 40, 0, 2),
    ]
}

fn pins() -> HashMap<String, (String, String)> {
    [("vep-rs", "deadbeef"), ("fastvep", "v0.0.0")]
        .iter()
        .map(|(engine, binary)| (engine.to_string(), (SWEEP.to_string(), binary.to_string())))
        .collect()
}

fn pin(binary: &str) -> (String, String) {
    (SWEEP.to_string(), binary.to_string())
}

const FIRST_DATASET: DatasetLabels = DatasetLabels {
    dataset: "s01",
    assembly: "GRCh37",
    engine: "vep-rs",
    sweep_id: SWEEP,
    binary_md5: "deadbeef",
};

const SECOND_DATASET: DatasetLabels = DatasetLabels {
    dataset: "s02",
    assembly: "GRCh38",
    engine: "fastvep",
    sweep_id: SWEEP,
    binary_md5: "v0.0.0",
};

/// The pooled CSV of `engines`, each with the two terms, and what is given of `adjusted`.
fn pooled_csv(
    scratch: &TempDir,
    engines: &[&str],
    adjusted: Option<(&str, &Exclusions)>,
) -> PathBuf {
    let path = scratch.path.join("pooled.csv");
    let pooled: Vec<(String, Vec<PerClassCounts>)> = engines
        .iter()
        .map(|engine| (engine.to_string(), two_terms()))
        .collect();
    write_per_class_pooled_csv(&path, &pooled, &pins(), adjusted).unwrap();
    path
}

#[test]
fn a_per_dataset_csv_gives_back_the_three_counts_and_the_pin_of_its_engine() {
    let scratch = scratch();
    let path = scratch.path.join("by_dataset.csv");
    write_per_class_by_dataset_csv(&path, &two_terms(), &FIRST_DATASET, false).unwrap();
    let (back, pins) = read_per_class_by_dataset_csv(&path).unwrap();
    let [(engine, mut rows)] = <[_; 1]>::try_from(back).expect("one engine");
    assert_eq!(engine, "vep-rs");
    assert_eq!(
        pins,
        HashMap::from([("vep-rs".to_string(), pin("deadbeef"))])
    );
    let mut expected = two_terms();
    rows.sort_by(|a, b| a.term.cmp(&b.term));
    expected.sort_by(|a, b| a.term.cmp(&b.term));
    assert_eq!(rows, expected);
}

#[test]
fn an_append_adds_rows_under_the_one_header() {
    let scratch = scratch();
    let path = scratch.path.join("by_dataset.csv");
    write_per_class_by_dataset_csv(&path, &two_terms(), &FIRST_DATASET, false).unwrap();
    write_per_class_by_dataset_csv(&path, &two_terms(), &SECOND_DATASET, true).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    assert_eq!(text.matches("dataset,assembly,engine").count(), 1);
    let (back, pins) = read_per_class_by_dataset_csv(&path).unwrap();
    let mut engines: Vec<&str> = back.iter().map(|(engine, _)| engine.as_str()).collect();
    engines.sort_unstable();
    assert_eq!(engines, ["fastvep", "vep-rs"]);
    assert_eq!(pins["fastvep"], pin("v0.0.0"));
}

#[test]
fn an_append_under_another_header_is_refused() {
    let scratch = scratch();
    let path = written(
        &scratch,
        "by_dataset.csv",
        "term,f1\nmissense_variant,1.0\n",
    );
    let appended = write_per_class_by_dataset_csv(&path, &two_terms(), &FIRST_DATASET, true);
    let message = refusal(appended);
    assert!(message.contains("cannot append"), "{message}");
}

#[test]
fn an_append_to_a_file_that_is_not_there_writes_it_with_its_header_and_its_directory() {
    let scratch = scratch();
    let path = scratch.path.join("nested").join("by_dataset.csv");
    write_per_class_by_dataset_csv(&path, &two_terms(), &FIRST_DATASET, true).unwrap();
    assert!(std::fs::read_to_string(&path)
        .unwrap()
        .starts_with("dataset,"));
}

#[test]
fn the_pooled_csv_names_the_eight_columns_its_readers_select() {
    let scratch = scratch();
    let path = pooled_csv(&scratch, &["vep-rs"], None);
    let (header, _) = csv::read_table(&path, b',').unwrap();
    for column in [
        "engine",
        "term",
        "vep_tuples",
        "engine_tuples",
        "intersection",
        "f1",
        "sweep_id",
        "binary_md5",
    ] {
        assert!(header.iter().any(|name| name == column), "{column}");
    }
}

#[test]
fn the_f1_of_a_pooled_row_is_the_half_up_rounding_of_the_row_s_own_integers() {
    let scratch = scratch();
    let path = pooled_csv(&scratch, &["vep-rs"], None);
    // `Decimal(str(2 * i / (v + e))).quantize(Decimal("0.000001"), rounding=ROUND_HALF_UP)`
    // for the (v, e, i) of each row.
    let rounded = [((100, 98, 97), "0.979798"), ((40, 42, 40), "0.975610")];
    let rows = csv_rows(&path);
    assert_eq!(rows.len(), rounded.len());
    for (row, (integers, f1)) in rows.iter().zip(rounded) {
        let cell = |column: &str| row[column].parse::<u64>().unwrap();
        assert_eq!(
            (
                cell("vep_tuples"),
                cell("engine_tuples"),
                cell("intersection")
            ),
            integers
        );
        assert_eq!(row["f1"], f1);
    }
}

#[test]
fn the_pooled_csv_writes_vep_rs_before_any_other_engine() {
    let scratch = scratch();
    let path = pooled_csv(&scratch, &["fastvep", "vep-rs"], None);
    let rows = csv_rows(&path);
    assert_eq!(rows.first().unwrap()["engine"], "vep-rs");
    assert_eq!(rows.last().unwrap()["engine"], "fastvep");
}

#[test]
fn a_per_class_csv_without_a_count_column_is_refused_by_column() {
    let scratch = scratch();
    let path = written(
        &scratch,
        "bad.csv",
        "engine,term,vep_tuples,sweep_id,binary_md5\nvep-rs,t,1,s,b\n",
    );
    let message = refusal(read_per_class_by_dataset_csv(&path));
    assert!(message.contains("engine_only"), "{message}");
}

#[test]
fn a_per_class_csv_row_of_more_vep_only_than_vep_tuples_is_refused() {
    let scratch = scratch();
    let path = written(
        &scratch,
        "bad.csv",
        &format!("{BY_DATASET_HEADER}s01,GRCh37,vep-rs,s,b,t,1,5,0\n"),
    );
    let message = refusal(read_per_class_by_dataset_csv(&path));
    assert!(message.contains("5 VEP-only"), "{message}");
}

#[test]
fn a_per_class_csv_with_one_engine_under_two_pins_is_refused() {
    let scratch = scratch();
    let path = written(
        &scratch,
        "mixed.csv",
        &format!(
            "{BY_DATASET_HEADER}s01,GRCh37,vep-rs,sweep-A,aaaa,t,10,1,0\n\
             s02,GRCh38,vep-rs,sweep-B,bbbb,t,10,1,0\n"
        ),
    );
    let message = refusal(read_per_class_by_dataset_csv(&path));
    assert!(
        message.contains("more than one (sweep_id, binary_md5)"),
        "{message}"
    );
}

#[test]
fn every_pooled_row_carries_the_pin_of_its_engine() {
    let scratch = scratch();
    let path = pooled_csv(&scratch, &["fastvep", "vep-rs"], None);
    let mut seen: Vec<(String, String, String)> = csv_rows(&path)
        .iter()
        .map(|row| {
            (
                row["engine"].clone(),
                row["sweep_id"].clone(),
                row["binary_md5"].clone(),
            )
        })
        .collect();
    seen.sort_unstable();
    seen.dedup();
    let pinned =
        |engine: &str, binary: &str| (engine.to_string(), SWEEP.to_string(), binary.to_string());
    assert_eq!(
        seen,
        [pinned("fastvep", "v0.0.0"), pinned("vep-rs", "deadbeef")]
    );
}

/// A reference output of three rows, one on a contig outside the human set.
fn reference_output(scratch: &TempDir) -> PathBuf {
    written(
        scratch,
        "output.txt",
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\n\
         v\t21:100\tA\tG\tENST1\tTranscript\tmissense_variant\n\
         v\t21:200\tT\tG\tENST2\tTranscript\tintron_variant,splice_region_variant\n\
         v\tKI270728.1:5\tC\tG\tENST3\tTranscript\tintron_variant\n",
    )
}

#[test]
fn the_fingerprint_reads_the_name_size_and_modification_time_and_the_contig_filter() {
    let scratch = scratch();
    let output = reference_output(&scratch);
    let files = [output.clone()];
    let human: Vec<String> = CANONICAL_CONTIGS.iter().map(|c| c.to_string()).collect();
    let fingerprint = vep_totals_fingerprint(&files, Some(&human)).unwrap();
    let meta = std::fs::metadata(&output).unwrap();
    let modified = meta.modified().unwrap();
    let nanos = modified.duration_since(UNIX_EPOCH).unwrap().as_nanos();
    assert!(
        fingerprint.contains(&format!("output.txt:{}:{nanos}", meta.len())),
        "{fingerprint}"
    );
    assert_ne!(fingerprint, vep_totals_fingerprint(&files, None).unwrap());
    let touched = File::options().write(true).open(&output).unwrap();
    touched
        .set_modified(modified + Duration::from_secs(1))
        .unwrap();
    assert_ne!(
        fingerprint,
        vep_totals_fingerprint(&files, Some(&human)).unwrap()
    );
}

#[test]
fn a_cache_that_is_absent_or_of_another_fingerprint_or_unreadable_gives_no_totals() {
    let scratch = scratch();
    let cache = scratch.path.join("cache").join("s01.csv");
    assert!(read_vep_totals_cache(&cache, "fp").is_none());
    write_vep_totals_cache(&cache, "fp", &counter(&[("missense_variant", 3)])).unwrap();
    let cached = read_vep_totals_cache(&cache, "fp").expect("the totals written");
    assert_eq!(counted(&cached), [("missense_variant", 3)]);
    assert!(read_vep_totals_cache(&cache, "other").is_none());
    std::fs::write(
        &cache,
        "# vep_totals_fingerprint=fp\nterm,vep_tuples\nx,notanumber\n",
    )
    .unwrap();
    assert!(read_vep_totals_cache(&cache, "fp").is_none());
    assert!(!cache.with_file_name("s01.csv.tmp").exists());
}

#[test]
fn the_adjusted_f1_takes_the_excluded_pairs_from_both_sides() {
    // The intersection is 97, with three VEP-only tuples and one engine-only tuple.
    let row = &two_terms()[0];
    assert_eq!(adjusted_f1_str(row, (0, 0)), Ok(row.f1_str()));
    // `2 * 97 / (2 * 97 + 1 + 0)` under the half-up rounding.
    assert_eq!(adjusted_f1_str(row, (2, 1)).as_deref(), Ok("0.994872"));
    assert_eq!(adjusted_f1_str(row, (3, 1)).as_deref(), Ok("1.000000"));
}

#[test]
fn an_exclusion_above_the_one_sided_count_it_is_taken_from_is_refused() {
    assert!(adjusted_f1_str(&two_terms()[0], (4, 0)).is_err());
}

#[test]
fn the_per_term_exclusions_are_read_from_their_two_columns_and_need_both() {
    let scratch = scratch();
    let good = written(
        &scratch,
        "per_term.tsv",
        "term\tvep_only_raw\tengine_only_raw\tvep_only_excluded\tengine_only_excluded\
         \tvep_only_open\tengine_only_open\n\
         missense_variant\t3\t1\t2\t1\t1\t0\n",
    );
    assert_eq!(
        read_per_term_exclusions(&good).unwrap(),
        Exclusions::from([("missense_variant".to_string(), (2, 1))])
    );
    let bad = written(
        &scratch,
        "bad.tsv",
        "term\tvep_only_raw\nmissense_variant\t3\n",
    );
    let message = refusal(read_per_term_exclusions(&bad));
    assert!(message.contains("engine_only_excluded"), "{message}");
}

#[test]
fn the_pooled_csv_fills_the_adjusted_columns_for_the_one_engine_named() {
    let scratch = scratch();
    let exclusions = Exclusions::from([("missense_variant".to_string(), (2, 1))]);
    let path = pooled_csv(
        &scratch,
        &["vep-rs", "fastvep"],
        Some(("vep-rs", &exclusions)),
    );
    let rows = csv_rows(&path);
    let row = |engine: &str, term: &str| {
        let found = rows
            .iter()
            .find(|row| row["engine"] == engine && row["term"] == term);
        found.expect("the row")
    };
    let [missense, intron] = <[_; 2]>::try_from(two_terms()).unwrap();
    let adjusted = row("vep-rs", "missense_variant");
    assert_eq!(
        Ok(&adjusted["adj_f1"]),
        adjusted_f1_str(&missense, (2, 1)).as_ref()
    );
    assert_eq!(adjusted["vep_only_excluded"], "2");
    assert_eq!(adjusted["engine_only_excluded"], "1");
    // A term the per-term file does not name excludes nothing and reads its raw F1.
    let unnamed = row("vep-rs", "intron_variant");
    assert_eq!(unnamed["adj_f1"], intron.f1_str());
    assert_eq!(unnamed["vep_only_excluded"], "0");
    for term in ["missense_variant", "intron_variant"] {
        for column in ["adj_f1", "vep_only_excluded", "engine_only_excluded"] {
            assert_eq!(row("fastvep", term)[column], "", "{term} {column}");
        }
    }
    // The adjusted value follows from the row's own integers, as `f1` does:
    // `2 * 97 / (2 * 97 + (3 - 2) + (1 - 1))` under the half-up rounding.
    let cell = |column: &str| adjusted[column].parse::<u64>().unwrap();
    assert_eq!(
        (
            cell("intersection"),
            cell("vep_only") - cell("vep_only_excluded"),
            cell("engine_only") - cell("engine_only_excluded"),
        ),
        (97, 1, 0)
    );
    assert_eq!(adjusted["adj_f1"], "0.994872");
}

#[test]
fn the_pooled_csv_leaves_the_adjusted_f1_empty_without_exclusions() {
    let scratch = scratch();
    let path = pooled_csv(&scratch, &["vep-rs"], None);
    let rows = csv_rows(&path);
    assert_eq!(rows.len(), 2);
    assert!(rows.iter().all(|row| row["adj_f1"].is_empty()));
}
