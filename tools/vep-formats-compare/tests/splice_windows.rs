// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The splice windows of a cache's transcripts and the class that reads them.

mod compare_common;
mod intron_ends_common;
mod support;

use compare_common::count_map;
use compare_common::counts;
use compare_common::default_row;
use compare_common::default_text;
use compare_common::documented_with;
use compare_common::entry;
use compare_common::json;
use compare_common::verdict;
use intron_ends_common::one_shard_cache;
use intron_ends_common::option_help;
use std::ffi::OsString;
use std::path::Path;
use std::path::PathBuf;
use vep_formats_compare::cli;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::splice_windows::SpliceWindows;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::allele_is_indel;
use vep_formats_compare::shared::uploaded_alleles;

const MAXENTSCAN_CLASS: &str = "maxentscan_overlapping_windows";

/// Chromosome 1 transcripts of one exon layout: an exon of three bases, 201 to 203, between an
/// intron from 101 to 200 and one from 204 to 300. On ENST1 (forward) the first intron's
/// acceptor window is 181 to 203 and the second's donor window 201 to 209; on ENST2 (reverse)
/// the second intron comes first, its acceptor window 201 to 223, and the first intron's donor
/// window is 195 to 203. ENST3 has no such exon: one intron from 101 to 300.
fn maxentscan_transcripts() -> Value {
    let exons = [
        r#"{"start": 1, "end": 100}"#,
        r#"{"start": 201, "end": 203}"#,
        r#"{"start": 301, "end": 400}"#,
    ];
    json(&format!(
        r#"[{{"stable_id": "ENST1", "strand": 1, "exons": [{all}]}}, {{"stable_id": "ENST2", "strand": -1, "exons": [{all}]}}, {{"stable_id": "ENST3", "strand": 1, "exons": [{}, {}]}}]"#,
        exons[0],
        exons[2],
        all = exons.join(", ")
    ))
}

/// A donor window's nine bases and an acceptor window's twenty-three.
fn window_texts() -> (String, String) {
    ("AC".repeat(4) + "G", "AC".repeat(11) + "G")
}

/// One substitution's row as the reference and as the candidate write it: the three scores
/// of two different windows and, under `verbose`, the two window texts, the reference's a
/// donor's and the candidate's an acceptor's.
fn maxentscan_texts(location: &str, feature: &str, verbose: bool) -> (String, String) {
    let (donor, acceptor) = window_texts();
    let side = |scores: [&str; 3], window: &str| {
        let mut extra = "IMPACT=LOW;STRAND=1".to_string();
        for (name, score) in ["MaxEntScan_ref", "MaxEntScan_alt", "MaxEntScan_diff"]
            .iter()
            .zip(scores)
        {
            extra.push_str(&format!(";{name}={score}"));
        }
        if verbose {
            extra.push_str(&format!(
                ";MaxEntScan_ref_seq={window};MaxEntScan_alt_seq={}T",
                &window[..window.len() - 1]
            ));
        }
        default_text(&[default_row(
            "a",
            location,
            "G",
            feature,
            "splice_region_variant",
            &extra,
        )])
    };
    (
        side(["8.100", "2.300", "5.800"], &donor),
        side(["6.400", "1.100", "5.300"], &acceptor),
    )
}

fn windows_of(cache: &Path, input_lines: Option<&[String]>) -> SpliceWindows {
    SpliceWindows::new(cache, input_lines).expect("the cache has its transcripts")
}

fn side_with(windows: SpliceWindows) -> Side {
    Side {
        splice_windows: Some(windows),
        ..Side::default()
    }
}

fn compare_maxentscan(reference: &str, candidate: &str, side: &Side) -> (Value, Documented) {
    let mut documented = documented_with(&[MAXENTSCAN_CLASS]);
    let report = compare_entries(
        &parse_default(reference).unwrap(),
        &parse_default(candidate).unwrap(),
        &mut documented,
        &CompareOptions::new(side),
    )
    .expect("the entries compare");
    (report, documented)
}

fn row_at(location: &str, feature: &str) -> Entry {
    parse_default(&maxentscan_texts(location, feature, false).0)
        .unwrap()
        .remove(0)
}

/// The rows the class takes and the rows of the same shape it leaves: a record in one of the
/// two windows or in neither, on a transcript without the exon or one the cache lacks, over
/// several bases, on another chromosome.
const IN_TWO_WINDOWS: [&str; 3] = ["ENST1", "ENST2", "ENST2.4"];
const OFF_THE_PREDICATE: [(&str, &str); 9] = [
    ("1:199", "ENST1"),
    ("1:205", "ENST1"),
    ("1:150", "ENST1"),
    ("1:199", "ENST2"),
    ("1:205", "ENST2"),
    ("1:202", "ENST3"),
    ("1:202", "ENST9"),
    ("1:201-203", "ENST1"),
    ("2:202", "ENST1"),
];

fn input_lines() -> Vec<String> {
    [
        "1\t202\t.\tA\tG\t.\t.\t.",
        "1\t199\t.\tA\tG\t.\t.\t.",
        "1\t202\t.\tAC\tA\t.\t.\t.",
    ]
    .iter()
    .map(|line| line.to_string())
    .collect()
}

fn entries_by_ordinal() -> Vec<Entry> {
    (0..4)
        .map(|ordinal| {
            entry(
                &format!("{ordinal}|G|ENST1|Transcript"),
                Some(ordinal),
                None,
                "G",
                "ENST1",
                "splice_region_variant",
                &[],
            )
        })
        .collect()
}

fn text(path: &Path) -> String {
    path.to_str().expect("a path of text").to_string()
}

fn cache_under(name: &str) -> (PathBuf, PathBuf) {
    let directory = support::scratch(name);
    let cache = one_shard_cache(&directory.join("cache"), &maxentscan_transcripts());
    (directory, cache)
}

#[test]
fn the_class_takes_a_row_the_cache_places_in_two_windows_and_no_other() {
    let (directory, cache) = cache_under("maxentscan-two-windows");
    let with_windows = side_with(windows_of(&cache, None));
    let without = Side::default();
    for feature in IN_TWO_WINDOWS {
        for verbose in [false, true] {
            let (reference, candidate) = maxentscan_texts("1:202", feature, verbose);
            let (report, documented) = compare_maxentscan(&reference, &candidate, &with_windows);
            assert_eq!(verdict(&report), "green", "{feature} {verbose}");
            assert_eq!(
                counts(&documented.applied),
                count_map(&[(MAXENTSCAN_CLASS, 1)]),
                "{feature} {verbose}"
            );
            let (report, documented) = compare_maxentscan(&reference, &candidate, &without);
            assert_eq!(verdict(&report), "red", "{feature} {verbose}");
            assert!(
                counts(&documented.applied).is_empty(),
                "{feature} {verbose}"
            );
        }
    }
    for (location, feature) in OFF_THE_PREDICATE {
        for verbose in [false, true] {
            let (reference, candidate) = maxentscan_texts(location, feature, verbose);
            let (report, documented) = compare_maxentscan(&reference, &candidate, &with_windows);
            assert_eq!(verdict(&report), "red", "{location} {feature} {verbose}");
            assert!(
                counts(&documented.applied).is_empty(),
                "{location} {feature} {verbose}"
            );
        }
    }
    let held: Vec<bool> = [
        ("1", 200),
        ("1", 201),
        ("1", 203),
        ("1", 204),
        ("chr1", 202),
    ]
    .iter()
    .map(|(name, position)| {
        windows_of(&cache, None)
            .holds(&row_at(&format!("{name}:{position}"), "ENST1"))
            .unwrap()
    })
    .collect();
    assert_eq!(held, [false, true, true, false, true]);
    // A window text is taken where the reference writes a donor's nine bases and the candidate
    // an acceptor's twenty-three, and in no other shape.
    let (donor, acceptor) = window_texts();
    let (reference, candidate) = maxentscan_texts("1:202", "ENST1", true);
    let reversed: String = donor.chars().rev().collect();
    for other in [
        candidate.replace(&acceptor, &reversed),
        candidate.replace(&acceptor, &format!("{acceptor}A")),
        candidate.replace(&format!(";MaxEntScan_ref_seq={acceptor}"), ""),
        candidate.replace(&acceptor, &acceptor.to_lowercase()),
    ] {
        assert_ne!(other, candidate);
        assert_eq!(
            verdict(&compare_maxentscan(&reference, &other, &with_windows).0),
            "red"
        );
    }
    assert_eq!(
        verdict(&compare_maxentscan(&candidate, &reference, &with_windows).0),
        "red"
    );
    // A VCF or JSON entry is placed by its input line: one base of reference at the position.
    let placed = windows_of(&cache, Some(&input_lines()));
    let entries = entries_by_ordinal();
    let held: Vec<bool> = entries.iter().map(|e| placed.holds(e).unwrap()).collect();
    assert_eq!(held, [true, false, false, false]);
    assert!(!with_windows
        .splice_windows
        .as_ref()
        .unwrap()
        .holds(&entries[0])
        .unwrap());
    // Through the command line: without --json-cache the class takes no row.
    let (reference, candidate) = maxentscan_texts("1:202", "ENST1", true);
    let reference_path = directory.join("reference.txt");
    let candidate_path = directory.join("candidate.txt");
    let report_path = directory.join("report.json");
    std::fs::write(&reference_path, reference).unwrap();
    std::fs::write(&candidate_path, candidate).unwrap();
    let argv = |more: &[&str]| -> Vec<OsString> {
        [
            "--format",
            "default",
            "--expected",
            &text(&reference_path),
            "--actual",
            &text(&candidate_path),
            "--masks",
            MAXENTSCAN_CLASS,
            "--report",
            &text(&report_path),
        ]
        .iter()
        .copied()
        .chain(more.iter().copied())
        .map(OsString::from)
        .collect()
    };
    assert_eq!(support::run(&directory, &argv(&[])).code, 1);
    assert_eq!(
        support::compact(support::member(
            &support::report(&report_path),
            &["verdict"]
        )),
        r#""red""#
    );
    let ran = support::run(&directory, &argv(&["--json-cache", &text(&cache)]));
    assert_eq!(ran.code, 0, "{}", ran.err());
    assert_eq!(
        support::compact(support::member(
            &support::report(&report_path),
            &["masks", "applied"]
        )),
        format!(r#"{{"{MAXENTSCAN_CLASS}": 1}}"#)
    );
    let none = text(&directory.join("none"));
    let ran = support::run(&directory, &argv(&["--json-cache", &none]));
    assert_eq!(ran.code, 2);
    assert!(ran.err().contains(&format!(
        "--json-cache: {none} has no transcripts/ directory"
    )));
    assert!(
        option_help(cli::build_parser().format_help(), "--json-cache").contains(MAXENTSCAN_CLASS)
    );
}

#[test]
fn the_windows_hold_a_one_base_record_in_an_acceptor_and_the_next_donor_window() {
    let (_, cache) = cache_under("splice-windows-hold");
    let windows = windows_of(&cache, None);
    for feature in IN_TWO_WINDOWS {
        assert!(
            windows.holds(&row_at("1:202", feature)).unwrap(),
            "{feature}"
        );
    }
    for (location, feature) in OFF_THE_PREDICATE {
        assert!(
            !windows.holds(&row_at(location, feature)).unwrap(),
            "{location} {feature}"
        );
    }
    let held: Vec<bool> = [
        ("1", 200),
        ("1", 201),
        ("1", 203),
        ("1", 204),
        ("chr1", 202),
    ]
    .iter()
    .map(|(name, position)| {
        windows_of(&cache, None)
            .holds(&row_at(&format!("{name}:{position}"), "ENST1"))
            .unwrap()
    })
    .collect();
    assert_eq!(held, [false, true, true, false, true]);
    let placed = windows_of(&cache, Some(&input_lines()));
    let entries = entries_by_ordinal();
    let held: Vec<bool> = entries.iter().map(|e| placed.holds(e).unwrap()).collect();
    assert_eq!(held, [true, false, false, false]);
    assert_eq!(placed.place(&entries[0]), Some(("1".to_string(), 202)));
    assert_eq!(placed.place(&entries[2]), None);
    assert!(!windows.holds(&entries[0]).unwrap());
    let refused = SpliceWindows::new(&cache.join("none"), None)
        .unwrap_err()
        .to_string();
    assert_eq!(
        refused,
        format!(
            "--json-cache: {} has no transcripts/ directory",
            text(&cache.join("none"))
        )
    );
}

#[test]
fn allele_shapes_and_the_class_of_the_overlapping_windows() {
    let spelt = |reference: &str, alternates: &[&str]| {
        Some((
            reference.to_string(),
            alternates
                .iter()
                .map(|allele| allele.to_string())
                .collect::<Vec<_>>(),
        ))
    };
    assert_eq!(
        uploaded_alleles("21_10542414_A/T/AG"),
        spelt("A", &["T", "AG"])
    );
    assert_eq!(uploaded_alleles("rs123"), None);
    assert_eq!(
        uploaded_alleles("chr1_KI270706v1_random_100_A/AAC"),
        spelt("A", &["AAC"])
    );
    assert_eq!(uploaded_alleles("del_padded"), None);
    assert_eq!(allele_is_indel("A", &["AT"], "T", false), Some(true));
    assert_eq!(allele_is_indel("A", &["G"], "G", false), Some(false));
    assert_eq!(allele_is_indel("AT", &["A"], "-", false), Some(true));
    assert_eq!(allele_is_indel("A", &["<DEL>"], "deletion", false), None);
    assert_eq!(
        allele_is_indel("G", &["[21:25881361[G"], "[21:25881361[G", false),
        None
    );
    // A lone alternate keeps its shape however the allele is shortened; an allele string of
    // the JSON format is already chopped, so its parts are read as they are.
    assert_eq!(
        allele_is_indel("TTTTT", &["TTTTTT"], "T", false),
        Some(true)
    );
    assert_eq!(allele_is_indel("TCT", &["TCTC"], "C", false), Some(true));
    assert_eq!(allele_is_indel("T", &["TT", "-"], "TT", true), Some(true));
    assert_eq!(allele_is_indel("T", &["TT", "-"], "-", true), Some(true));
    // The shared leading base is dropped only when every alternate shares it.
    assert_eq!(allele_is_indel("A", &["AT", "G"], "AT", false), Some(true));
    assert_eq!(allele_is_indel("A", &["AT", "G"], "G", false), Some(false));
    assert_eq!(allele_is_indel("A", &["AT", "G"], "T", false), None);
    assert_eq!(allele_is_indel("A", &["AT", "T"], "AT", false), Some(true));
    assert_eq!(allele_is_indel("A", &["AT", "T"], "T", false), Some(false));
    assert_eq!(allele_is_indel("A", &["AT", "*"], "T", false), Some(true));
    // An allele spelt with a `*` is exempt from the chop whatever else it carries.
    assert_eq!(
        allele_is_indel("A", &["AG", "<DEL:*>"], "G", false),
        Some(true)
    );
    assert_eq!(
        allele_is_indel("-", &["G", "<DEL:*>"], "G", true),
        Some(true)
    );
    assert_eq!(
        allele_is_indel("A", &["T", "<DEL:*>"], "T", false),
        Some(false)
    );
    assert_eq!(
        allele_is_indel("AC", &["ACT", "AG"], "CT", false),
        Some(true)
    );
    assert_eq!(
        allele_is_indel("AC", &["ACT", "AG"], "G", false),
        Some(false)
    );
    let (_, cache) = cache_under("maxentscan-allele-shapes");
    let with_windows = side_with(windows_of(&cache, None));
    let mes = default_text(&[default_row(
        "a",
        "1:202",
        "G",
        "ENST1",
        "splice_region_variant",
        "IMPACT=LOW;STRAND=1;MaxEntScan_ref=8.1;MaxEntScan_alt=2.3;MaxEntScan_diff=5.8",
    )]);
    let other_window = mes
        .replace("8.1", "6.4")
        .replace("2.3", "1.1")
        .replace("5.8", "5.3");
    let (report, documented) = compare_maxentscan(&mes, &other_window, &with_windows);
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(MAXENTSCAN_CLASS, 1)])
    );
    let without_a_score = mes.replace("MaxEntScan_alt=2.3;", "");
    assert_eq!(
        verdict(&compare_maxentscan(&mes, &without_a_score, &with_windows).0),
        "red"
    );
    let no_number = mes.replace("8.1", "abc");
    assert_eq!(
        verdict(&compare_maxentscan(&mes, &no_number, &with_windows).0),
        "red"
    );
}
