// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The settling step on threads: one pair of outputs, and the reference against a tuple
//! file of the engine's values, scored on 1, 2, 5 and 16 threads, the sorted files cut into
//! ranges of one Location, of a few and into one range, under classes that set values
//! aside, classes that leave them and a class that refuses. Every run writes the same
//! report and the same residual files, and a refused run ends in the refusal of the first
//! Location that refuses.

mod support;

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use support::{names, refusal, row_with, Planted, PlantedRules, Run, Scratch};
use vep_plugin_pairs::rules::Refusal;
use vep_plugin_pairs::score::report::dumps;

const THREADS: [i64; 4] = [1, 2, 5, 16];
/// The bytes of a sorted file a range holds: a figure every line is longer than, so that a
/// range is one Location; one a few Locations fill; and the program's own, which leaves
/// these outputs in one range.
const RANGE_BYTES: [Option<u64>; 3] = [Some(1), Some(400), None];

/// Forty-eight Locations on three chromosomes, among them Locations that open another
/// (`1:10`, `1:100`), in an order that is not the sorted one.
fn locations() -> Vec<String> {
    let positions = [
        5, 10, 100, 1000, 11, 110, 2, 20, 200, 25, 250, 3, 30, 31, 300, 7,
    ];
    ["2", "1", "10"]
        .iter()
        .flat_map(|chromosome| {
            positions
                .iter()
                .map(move |position| format!("{chromosome}:{position}"))
        })
        .collect()
}

/// The rows of one side, `0` the reference's and `1` the engine's. Every Location has a row
/// on `ENST1` under `T`, every other one a row on `ENST2`, every third a row under `G`, and
/// every tenth two rows under one key. The two sides write the same `PLANT_label`; the
/// other fields differ on some rows in a way one planted class accepts and on others in a
/// way none does. The reference's `MaxEntScan_diff` reads `refuse` at the Locations of
/// `refusing`.
fn rows(side: usize, refusing: &[&str]) -> Vec<String> {
    let mut rows = Vec::new();
    for (at, location) in locations().iter().enumerate() {
        let alleles: &[&str] = if at % 3 == 0 { &["T", "G"] } else { &["T"] };
        let features: &[&str] = if at % 2 == 0 {
            &["ENST1", "ENST2"]
        } else {
            &["ENST1"]
        };
        for (a, allele) in alleles.iter().enumerate() {
            for (f, feature) in features.iter().enumerate() {
                let n = at * 4 + a * 2 + f;
                let reference = side == 0;
                let mut extra = vec!["IMPACT=MODIFIER".to_string()];
                // A value no class reads, left where the engine writes another.
                extra.push(match (n % 4, reference) {
                    (0, false) => format!("PLANT_score=0.9{n}"),
                    _ => format!("PLANT_score=0.{n}"),
                });
                extra.push("PLANT_label=benign".to_string());
                match (n % 5, reference) {
                    (0, false) => extra.push("LOEUF=0".to_string()),
                    (1, false) => extra.push("LOEUF=7".to_string()),
                    _ => {}
                }
                if refusing.contains(&location.as_str()) {
                    extra.push(format!(
                        "MaxEntScan_diff={}",
                        if reference { "refuse" } else { "other" }
                    ));
                } else if n % 7 == 3 {
                    extra.push(format!(
                        "MaxEntScan_diff={}",
                        if reference { "x" } else { "y" }
                    ));
                } else if n % 3 == 1 {
                    extra.push(format!(
                        "MaxEntScan_diff={}",
                        if reference { "1.5" } else { "2.5" }
                    ));
                }
                if n % 6 == 2 && reference {
                    extra.push("GWAS_pmid=17".to_string());
                }
                if n % 8 == 4 {
                    extra.push(format!(
                        "NearestExonJB={}",
                        if reference { "a,b,c" } else { "c,a,b" }
                    ));
                }
                let cdna = if n % 9 == 4 { "1" } else { "-" };
                if n % 9 == 4 || n % 9 == 5 {
                    extra.push(format!(
                        "5UTR_consequence={}",
                        if reference { "x" } else { "y" }
                    ));
                }
                rows.push(row_with(
                    location,
                    allele,
                    feature,
                    &extra.join(";"),
                    cdna,
                    "v",
                ));
                if at % 10 == 9 {
                    let value = if reference { "3.5" } else { "4.5" };
                    rows.push(row_with(
                        location,
                        allele,
                        feature,
                        &format!("IMPACT=MODIFIER;MaxEntScan_diff={value}"),
                        "-",
                        "w",
                    ));
                }
            }
        }
    }
    if side == 1 {
        rows.reverse();
    }
    rows
}

/// The report's text and the residual files, by name, of one run that goes through.
type Written = (String, BTreeMap<String, Vec<u8>>);

fn written(
    planted: &Planted,
    scratch: &Scratch,
    run: &Run,
    threads: i64,
    range_bytes: Option<u64>,
) -> Written {
    let residuals = scratch.path(&format!("residuals-{threads}-{range_bytes:?}"));
    let mut run = run.clone().residuals(&residuals).threads(threads);
    run.range_bytes = range_bytes;
    let (reference, engine) = (scratch.path("reference.txt"), scratch.path("engine.txt"));
    let report = planted.score(scratch, &reference, &[("engine", engine.as_str())], &run);
    let files = names(&residuals)
        .into_iter()
        .map(|name| {
            let bytes = fs::read(Path::new(&residuals).join(&name)).unwrap();
            (name, bytes)
        })
        .collect();
    fs::remove_dir_all(&residuals).unwrap();
    assert!(scratch.names("work").is_empty());
    (dumps(&report), files)
}

/// Scores the outputs under `run` on one thread in one range, then on every thread count
/// and range size, and returns what the first run wrote once every run wrote the same.
fn same_on_every_thread_count(planted: &Planted, scratch: &Scratch, run: &Run) -> Written {
    let first = written(planted, scratch, run, 1, None);
    for threads in THREADS {
        for range_bytes in RANGE_BYTES {
            let other = written(planted, scratch, run, threads, range_bytes);
            assert!(
                first.0 == other.0,
                "the report on {threads} threads in ranges of {range_bytes:?}"
            );
            assert!(
                first.1 == other.1,
                "the residual files on {threads} threads in ranges of {range_bytes:?}"
            );
        }
    }
    first
}

fn lines(files: &BTreeMap<String, Vec<u8>>, name: &str) -> usize {
    files[name].iter().filter(|&&byte| byte == b'\n').count()
}

#[test]
fn classes_that_set_values_aside_score_the_same_on_any_thread_count() {
    let scratch = Scratch::new("settle-threads-aside");
    scratch.output("reference.txt", &rows(0, &[]));
    scratch.output("engine.txt", &rows(1, &[]));
    let planted = Planted::new();
    let (report, files) = same_on_every_thread_count(&planted, &scratch, &Run::new());
    // Every kind of outcome is among the scores: values set aside by a class that reads
    // rows, by one that reads none, and values no class takes.
    for class in [
        "planted_engine_zero",
        "planted_two_numbers",
        "planted_reference_value",
        "planted_row_cdna_one",
        "order_only",
    ] {
        let aside = report.lines().any(|line| {
            line.trim().starts_with(&format!("\"{class}\": ")) && !line.contains(": 0")
        });
        assert!(aside, "{class} set nothing aside");
    }
    for (name, least) in [
        ("Planted.engine.feature.tsv", 20),
        ("Planted.engine.allele.tsv", 20),
        ("LOEUF.engine.feature.tsv", 5),
        ("MaxEntScan.engine.feature.tsv", 5),
        ("UTRAnnotator.engine.feature.tsv", 5),
    ] {
        assert!(lines(&files, name) >= least, "{name}");
    }
    assert_eq!(lines(&files, "GWAS.engine.feature.tsv"), 0);
    assert_eq!(lines(&files, "NearestExonJB.engine.feature.tsv"), 0);
}

#[test]
fn classes_that_leave_values_score_the_same_on_any_thread_count() {
    let scratch = Scratch::new("settle-threads-left");
    scratch.output("reference.txt", &rows(0, &[]));
    scratch.output("engine.txt", &rows(1, &[]));
    let planted = Planted::new();
    // The classes that are not applied accept their values and leave them.
    let run = Run::new().applied(&["planted_engine_zero", "order_only"]);
    let (report, files) = same_on_every_thread_count(&planted, &scratch, &run);
    assert!(report.contains("\"not_applied\": {\n"));
    assert!(lines(&files, "GWAS.engine.feature.tsv") >= 5);
    assert!(lines(&files, "GWAS.engine.allele.tsv") >= 5);
    assert_eq!(lines(&files, "NearestExonJB.engine.feature.tsv"), 0);
    let all = same_on_every_thread_count(&planted, &scratch, &Run::new());
    assert!(
        lines(&files, "MaxEntScan.engine.feature.tsv")
            > lines(&all.1, "MaxEntScan.engine.feature.tsv")
    );
}

/// The plugin values of default-format rows as the lines of a tuple file.
fn tuple_lines(rows: &[String]) -> String {
    let mut lines = String::new();
    for row in rows {
        let columns: Vec<&str> = row.trim_end_matches('\n').split('\t').collect();
        let (location, allele, feature) = (columns[1], columns[2], columns[4]);
        for pair in columns[13].split(';') {
            let (name, value) = pair.split_once('=').unwrap();
            if name != "IMPACT" {
                lines.push_str(&format!(
                    "{location}\t{allele}\t{feature}\t{name}\t{value}\n"
                ));
            }
        }
    }
    lines
}

#[test]
fn tuple_file_engine_scores_the_same_on_any_thread_count() {
    let scratch = Scratch::new("settle-threads-tuples");
    scratch.output("reference.txt", &rows(0, &[]));
    scratch.file("engine.txt", tuple_lines(&rows(1, &[])));
    let planted = Planted::new();
    let (report, files) = same_on_every_thread_count(&planted, &scratch, &Run::new());
    assert!(report.contains("\"format\": \"tuples\""));
    for class in ["planted_engine_zero", "planted_two_numbers", "order_only"] {
        let aside = report.lines().any(|line| {
            line.trim().starts_with(&format!("\"{class}\": ")) && !line.contains(": 0")
        });
        assert!(aside, "{class} set nothing aside");
    }
    assert!(lines(&files, "Planted.engine.feature.tsv") >= 20);
    assert!(lines(&files, "MaxEntScan.engine.feature.tsv") >= 5);
}

#[test]
fn class_that_refuses_ends_the_run_at_the_first_location_on_any_thread_count() {
    let scratch = Scratch::new("settle-threads-refused");
    // In the order of the sorted files `10:200` comes before `1:300` and `2:25`.
    let refusing = ["2:25", "1:300", "10:200"];
    scratch.output("reference.txt", &rows(0, &refusing));
    scratch.output("engine.txt", &rows(1, &refusing));
    let planted = Planted {
        rules: PlantedRules::planted().refusing_class("planted_two_numbers", &[], |consult| {
            match consult.expected {
                "refuse" => Err(Refusal(format!(
                    "planted refusal at {}",
                    consult.e.location.as_deref().unwrap_or("no Location")
                ))),
                _ => Ok(false),
            }
        }),
        ..Planted::new()
    };
    let (reference, engine) = (scratch.path("reference.txt"), scratch.path("engine.txt"));
    for threads in THREADS {
        for range_bytes in RANGE_BYTES {
            let residuals = scratch.path("residuals");
            let mut run = Run::new().residuals(&residuals).threads(threads);
            run.range_bytes = range_bytes;
            let refused = refusal(planted.try_score(
                &scratch,
                &reference,
                &[("engine", engine.as_str())],
                &run,
            ));
            assert_eq!(
                refused, "planted refusal at 10:200",
                "{threads} threads in ranges of {range_bytes:?}"
            );
            assert!(names(&residuals).is_empty(), "{threads} threads");
            assert!(scratch.names("work").is_empty(), "{threads} threads");
        }
    }
}
