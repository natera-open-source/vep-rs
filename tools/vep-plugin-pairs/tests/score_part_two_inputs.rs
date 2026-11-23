// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The inputs of the `score` program under a planted registry and planted classes: a gzip
//! input and its digest, rows and values that share a key, the two formats and the inputs
//! that are refused, the order of the keys, drawn outputs held to multiset arithmetic, and
//! the summary of a run with its class lines.

mod support;

use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::os::unix::process::CommandExt;
use std::process::{Command, Stdio};
use support::{
    command, counts, entry, feature, gzip, reference_rows, refusal, residual_lines, row, row_with,
    sha256, Planted, Run, Scratch, CLASS_TABLE, HEADER, PLANTED,
};
use vep_plugin_pairs::score::reader::Reader;

/// A generator of draws that gives the same sequence in every run.
struct Draws(u64);

impl Draws {
    fn new(seed: u64) -> Draws {
        Draws(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1)
    }

    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 >> 24
    }

    fn below(&mut self, bound: usize) -> usize {
        (self.next() % bound as u64) as usize
    }

    fn choice<'a>(&mut self, of: &[&'a str]) -> &'a str {
        of[self.below(of.len())]
    }

    /// `bytes` drawn bytes as hexadecimal digits.
    fn hex(&mut self, bytes: usize) -> String {
        (0..bytes)
            .map(|_| format!("{:02x}", self.below(256)))
            .collect()
    }
}

#[test]
fn gzip_input_and_the_digest_of_the_file_as_stored() {
    let scratch = Scratch::new("gzip");
    let planted = Planted::new();
    let rows = reference_rows();
    let plain = scratch.output("reference.txt", &rows);
    let stored = gzip(&fs::read(&plain).unwrap());
    let packed = scratch.file("reference.txt.gz", &stored);
    let report = planted.score(&scratch, &packed, &[("engine", &plain)], &Run::new());
    assert_eq!(report["reference"]["sha256"], json!(sha256(&stored)));
    assert_eq!(
        report["engines"][0]["sha256"],
        json!(sha256(&fs::read(&plain).unwrap()))
    );
    assert_eq!(report["reference"]["rows"], rows.len());
    assert_eq!(report["engines"][0]["rows"], rows.len());
    assert_eq!(counts(feature(&report, "Planted")), [4; 6]);
}

#[test]
fn reader_digest_covers_the_file_whenever_it_is_taken() {
    for packed in [false, true] {
        let scratch = Scratch::new("reader-digest");
        let mut draws = Draws::new(116);
        let data: String = (0..40_000)
            .map(|n| format!("1:{n}\tT\tENST1\tPLANT_label\t{}\n", draws.hex(50)))
            .collect();
        let stored = match packed {
            true => gzip(data.as_bytes()),
            false => data.clone().into_bytes(),
        };
        assert!(stored.len() > 1 << 20, "{packed}");
        let path = scratch.file("input", &stored);
        let mut reader = Reader::open(&path).unwrap();
        let mut line = Vec::new();
        assert!(reader.line(&mut line).unwrap(), "{packed}");
        let first = &data.as_bytes()[..=data.find('\n').unwrap()];
        assert_eq!(line, first, "{packed}");
        assert_eq!(reader.sha256().unwrap(), sha256(&stored), "{packed}");
    }
}

#[test]
fn rows_that_share_a_key_are_never_collapsed() {
    let scratch = Scratch::new("shared-rows");
    let planted = Planted::new();
    let rows = reference_rows();
    let twice = scratch.output("twice.txt", &[rows.clone(), rows.clone()].concat());
    let once = scratch.output("once.txt", &rows);
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &twice,
        &[("once", &once), ("twice", &twice)],
        &Run::new().residuals(&residuals),
    );
    assert_eq!(report["reference"]["rows"], 2 * rows.len());
    let single = entry(&report, "Planted", "once", "feature");
    assert_eq!(counts(single), [8, 4, 4, 8, 4, 4]);
    assert_eq!(single["residual_rows"], 4);
    assert_eq!(single["shared_keys"], 3);
    assert_eq!(
        residual_lines(&residuals, "Planted", "once", "feature"),
        [
            "reference\t1:100\tT\tENST1\tPLANT_label\tbenign",
            "reference\t1:100\tT\tENST1\tPLANT_score\t0.5",
            "reference\t1:100\tT\tENST2\tPLANT_score\t0.5",
            "reference\t1:200-201\t-\tENST1\tPLANT_score\t0.25",
        ]
    );
    let double = entry(&report, "Planted", "twice", "feature");
    assert_eq!(counts(double), [8; 6]);
    assert_eq!(double["residual_rows"], 0);
    assert_eq!(double["shared_keys"], 3);
    for ident in ["once", "twice"] {
        let allele = entry(&report, "Planted", ident, "allele");
        assert_eq!(counts(allele), [3; 6], "{ident}");
        assert_eq!(allele["shared_keys"], 0, "{ident}");
    }
    let report = planted.score(&scratch, &once, &[("engine", &once)], &Run::new());
    let alone = feature(&report, "Planted");
    assert_eq!(counts(alone), [4; 6]);
    assert_eq!(alone["shared_keys"], 0);
}

#[test]
fn values_under_a_shared_key_pair_as_often_as_each_side_writes_them() {
    let scratch = Scratch::new("shared-values");
    let planted = Planted::new();
    let of = |uploaded: &str, location: &str, allele: &str, feature: &str, extra: &str| {
        row_with(location, allele, feature, extra, "-", uploaded)
    };
    let reference = scratch.output(
        "reference.txt",
        &[
            of("v1", "1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=a"),
            of("v2", "1:100", "T", "ENST1", "PLANT_score=0.7;PLANT_label=a"),
            of("v3", "1:100", "T", "ENST1", "PLANT_score=0.7"),
            of("v1", "1:100", "T", "ENST2", "PLANT_score=1.0"),
            of("v2", "1:100", "T", "ENST2", "PLANT_score=1.0"),
            row("1:900", "G", "ENST1", "PLANT_score=1"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            of("v2", "1:100", "T", "ENST1", "PLANT_score=0.7;PLANT_label=a"),
            of("v1", "1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=a"),
            of("v3", "1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=a"),
            of("v1", "1:100", "T", "ENST2", "PLANT_score=1"),
            of("v2", "1:100", "T", "ENST2", "PLANT_score=1"),
            row("1:900", "G", "ENST1", "PLANT_score=1"),
        ],
    );
    let listed = scratch.tuples(
        "engine.tsv",
        &[
            ["1:100", "T", "ENST1", "PLANT_score", "0.7"],
            ["1:100", "T", "ENST1", "PLANT_label", "a"],
            ["1:100", "T", "ENST1", "PLANT_score", "0.5"],
            ["1:100", "T", "ENST1", "PLANT_label", "a"],
            ["1:100", "T", "ENST1", "PLANT_score", "0.5"],
            ["1:100", "T", "ENST1", "PLANT_label", "a"],
            ["1:100", "T", "ENST2", "PLANT_score", "1"],
            ["1:100", "T", "ENST2", "PLANT_score", "1"],
            ["1:900", "G", "ENST1", "PLANT_score", "1"],
        ],
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine), ("listed", &listed)],
        &Run::new().residuals(&residuals).no_table(),
    );
    for ident in ["engine", "listed"] {
        let item = entry(&report, "Planted", ident, "feature");
        assert_eq!(counts(item), [8, 9, 7, 8, 9, 7], "{ident}");
        assert_eq!(item["text_only_differences"], 2, "{ident}");
        assert_eq!(item["residual_rows"], 3, "{ident}");
        assert_eq!(item["shared_keys"], 2, "{ident}");
        assert_eq!(
            residual_lines(&residuals, "Planted", ident, "feature"),
            [
                "engine\t1:100\tT\tENST1\tPLANT_label\ta",
                "engine\t1:100\tT\tENST1\tPLANT_score\t0.5",
                "reference\t1:100\tT\tENST1\tPLANT_score\t0.7",
            ],
            "{ident}"
        );
        let allele = entry(&report, "Planted", ident, "allele");
        assert_eq!(counts(allele), [5; 6], "{ident}");
        assert_eq!(allele["text_only_differences"], 1, "{ident}");
        assert_eq!(allele["shared_keys"], 1, "{ident}");
        assert!(
            residual_lines(&residuals, "Planted", ident, "allele").is_empty(),
            "{ident}"
        );
    }
}

#[test]
fn percent_encoded_value_equals_its_decoded_form_in_a_tuple_file() {
    let scratch = Scratch::new("percent");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "T", "ENST1", "PLANT_label=a%3Bb%3Dc")],
    );
    let engine = scratch.tuples(
        "engine.tsv",
        &[["1:100", "T", "ENST1", "PLANT_label", "a;b=c"]],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    assert_eq!(counts(feature(&report, "Planted")), [1; 6]);
}

#[test]
fn location_coordinates_compare_as_numbers() {
    let scratch = Scratch::new("coordinates");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("21:100-0200", "T", "ENST1", "PLANT_score=1"),
            row("21:²", "T", "ENST1", "PLANT_score=1"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("21:100-200", "T", "ENST1", "PLANT_score=1"),
            row("21:²", "T", "ENST1", "PLANT_score=1"),
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    assert_eq!(counts(feature(&report, "Planted")), [2; 6]);
}

#[test]
fn headerless_default_output_and_format_detection() {
    let scratch = Scratch::new("headerless");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let bare = scratch.file("engine.txt", rows.concat());
    let report = planted.score(&scratch, &reference, &[("engine", &bare)], &Run::new());
    assert_eq!(report["engines"][0]["format"], "default");
    assert_eq!(counts(feature(&report, "Planted")), [4; 6]);
    let headed = scratch.file(
        "headed.tsv",
        format!("{HEADER}1:100\tT\tENST1\tPLANT_score\t0.5\n"),
    );
    let refused =
        refusal(planted.try_score(&scratch, &reference, &[("engine", &headed)], &Run::new()));
    assert!(refused.contains("line 2 has 5 columns"), "{refused}");
}

#[test]
fn malformed_inputs_are_refused() {
    let scratch = Scratch::new("malformed");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let refused = |reference: &str, engines: &[(&str, &str)]| {
        refusal(planted.try_score(&scratch, reference, engines, &Run::new()))
    };
    let short = scratch.file(
        "short.tsv",
        "1:100\tT\tENST1\tPLANT_score\t0.5\n1:100\tT\tENST1\tPLANT_score\n",
    );
    let message = refused(&reference, &[("engine", &short)]);
    assert!(message.contains("line 2 has 4 columns"), "{message}");
    let listed = scratch.tuples("t.tsv", &[["1:100", "T", "ENST1", "PLANT_score", "1"]]);
    let message = refused(&listed, &[("engine", &reference)]);
    assert!(
        message.ends_with("line 1 has 5 columns; a default-format row has 14"),
        "{message}"
    );
    let empty = scratch.output("empty.txt", &[]);
    let message = refused(&empty, &[("engine", &reference)]);
    assert!(message.contains("no data row"), "{message}");
    let message = refused(&reference, &[]);
    assert!(message.contains("no engine"), "{message}");
    let odd = scratch.file(
        "odd.txt",
        format!("{HEADER}{}#Uploaded_variation\tLocation\n", rows[0]),
    );
    let message = refused(&reference, &[("engine", &odd)]);
    assert!(
        message.contains("line 3, the column header, names 2 columns, not 14"),
        "{message}"
    );
    let stored = gzip(&fs::read(&reference).unwrap());
    let cut = scratch.file("cut.txt.gz", &stored[..stored.len() - 12]);
    let message = refused(&cut, &[("engine", &reference)]);
    assert!(
        message.contains("cut.txt.gz: cannot read the gzip stream"),
        "{message}"
    );
    assert!(scratch.exists("work") && scratch.names("work").is_empty());
}

#[test]
fn keys_pair_in_the_order_of_their_bytes_whatever_the_callers_locale() {
    let scratch = Scratch::new("byte-order");
    let planted = Planted::new();
    let mut rows = Vec::new();
    for location in [
        "1:117329",
        "10:110755",
        "1:1173290",
        "1:117329-117330",
        "X:5",
        "1_KI270706v1_random:9",
    ] {
        for allele in ["-", "A", "a", "AT"] {
            for feature in ["-", "ENST1", "ENST10", "ENST1.2"] {
                rows.push(row(
                    location,
                    allele,
                    feature,
                    "PLANT_score=1;PLANT_label=x",
                ));
            }
        }
    }
    let reference = scratch.output("reference.txt", &rows);
    rows.reverse();
    let engine = scratch.output("engine.txt", &rows);
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    assert_eq!(counts(feature(&report, "Planted")), [192; 6]);
    // The built program under a locale whose collation is not the bytes' order.
    let no_class = scratch.file("no-class.tsv", "");
    let report_path = scratch.path("report.json");
    let done = command(
        &[
            "--reference",
            &reference,
            "--engine",
            &format!("engine={engine}"),
            "--plugins",
            "",
            "--custom-source",
            "Planted=PLANT_",
            "--class-table",
            &no_class,
            "--work-dir",
            &scratch.path("work"),
            "--report",
            &report_path,
        ],
        &[("LC_ALL", "en_US.UTF-8")],
    );
    assert!(
        done.status.success(),
        "{}",
        String::from_utf8_lossy(&done.stderr)
    );
    let report: Value = serde_json::from_slice(&fs::read(&report_path).unwrap()).unwrap();
    assert_eq!(counts(feature(&report, "Planted")), [192; 6]);
}

type Tuple = [String; 5];

/// Multiset arithmetic over (Location, Allele, Feature, field, value) tuples whose only
/// numerically equal texts are `0.1` and `0.10`, each tuple once where the Feature is
/// dropped: the two sides' counts, the pairs, the keys under which a side holds a field
/// more than once, and the unpaired tuples as residual lines in the tuples' order.
fn expected_counts(
    reference: &[Tuple],
    engine: &[Tuple],
    drop_feature: bool,
) -> (i64, i64, i64, i64, Vec<String>) {
    let grouped = |side: &[Tuple]| {
        let mut grouped: BTreeMap<Vec<String>, BTreeMap<String, i64>> = BTreeMap::new();
        for [location, allele, feature, name, value] in side {
            let mut key = vec![location.clone(), allele.clone()];
            if !drop_feature {
                key.push(feature.clone());
            }
            key.push(name.clone());
            let held = grouped.entry(key).or_default().entry(value.clone());
            let held = held.or_default();
            *held = if drop_feature { 1 } else { *held + 1 };
        }
        grouped
    };
    let sides = [grouped(reference), grouped(engine)];
    let keys: BTreeSet<&Vec<String>> = sides[0].keys().chain(sides[1].keys()).collect();
    let none = BTreeMap::new();
    let mut paired = 0;
    let mut left: Vec<(String, &str)> = Vec::new();
    for key in keys {
        let ref_values = sides[0].get(key).unwrap_or(&none);
        let eng_values = sides[1].get(key).unwrap_or(&none);
        let values: BTreeSet<&String> = ref_values.keys().chain(eng_values.keys()).collect();
        let mut ref_left: BTreeMap<&str, i64> = BTreeMap::new();
        let mut eng_left: BTreeMap<&str, i64> = BTreeMap::new();
        for value in values {
            let held = |values: &BTreeMap<String, i64>| values.get(value).copied().unwrap_or(0);
            let both = held(ref_values).min(held(eng_values));
            paired += both;
            ref_left.insert(value, held(ref_values) - both);
            eng_left.insert(value, held(eng_values) - both);
        }
        for (one, other) in [("0.1", "0.10"), ("0.10", "0.1")] {
            let (Some(ours), Some(theirs)) = (ref_left.get_mut(one), eng_left.get_mut(other))
            else {
                continue;
            };
            let as_numbers = (*ours).min(*theirs);
            paired += as_numbers;
            *ours -= as_numbers;
            *theirs -= as_numbers;
        }
        for (side, held) in [("reference", &ref_left), ("engine", &eng_left)] {
            for (value, times) in held {
                for _ in 0..*times {
                    left.push((format!("{}\t{value}", key.join("\t")), side));
                }
            }
        }
    }
    left.sort_by(|one, other| one.0.as_bytes().cmp(other.0.as_bytes()));
    let lines = left
        .iter()
        .map(|(columns, side)| format!("{side}\t{columns}"))
        .collect();
    let total = |values: &BTreeMap<String, i64>| values.values().sum::<i64>();
    let shared: BTreeSet<&[String]> = sides
        .iter()
        .flatten()
        .filter(|(_, values)| total(values) > 1)
        .map(|(key, _)| &key[..key.len() - 1])
        .collect();
    let totals = sides
        .each_ref()
        .map(|side| side.values().map(total).sum::<i64>());
    (totals[0], totals[1], paired, shared.len() as i64, lines)
}

#[test]
fn drawn_outputs_agree_with_multiset_arithmetic() {
    let scratch = Scratch::new("drawn");
    let planted = Planted::new();
    let mut draws = Draws::new(116);
    let locations = ["1:5", "1:50", "1:5-6", "10:5", "2:5", "MT:5"];
    let alleles = ["A", "-", "AT", "T"];
    let features = ["ENST1", "ENST10", "ENST2", "-"];
    let values = ["a", "b", "0.1", "0.10", "7"];
    let mut side = || -> (Vec<String>, Vec<Tuple>) {
        let (mut rows, mut written) = (Vec::new(), Vec::new());
        for _ in 0..400 {
            let key = [
                draws.choice(&locations),
                draws.choice(&alleles),
                draws.choice(&features),
            ];
            let mut pairs = Vec::new();
            for name in ["PLANT_score", "PLANT_label"] {
                if draws.below(10) < 7 {
                    pairs.push((name, draws.choice(&values)));
                }
            }
            let extra: Vec<String> = pairs
                .iter()
                .map(|(name, value)| format!("{name}={value}"))
                .collect();
            let extra = match extra.is_empty() {
                true => "-".to_string(),
                false => extra.join(";"),
            };
            rows.push(row(key[0], key[1], key[2], &extra));
            for (name, value) in pairs {
                written.push([key[0], key[1], key[2], name, value].map(str::to_string));
            }
        }
        (rows, written)
    };
    let (reference_rows, reference_tuples) = side();
    let (engine_rows, engine_tuples) = side();
    let reference = scratch.output("reference.txt", &reference_rows);
    let as_output = scratch.output("engine.txt", &engine_rows);
    let listed: Vec<[&str; 5]> = engine_tuples
        .iter()
        .map(|tuple| std::array::from_fn(|column| tuple[column].as_str()))
        .collect();
    let as_tuples = scratch.tuples("engine.tsv", &listed);
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("output", &as_output), ("tuples", &as_tuples)],
        &Run::new().residuals(&residuals),
    );
    for (projection, drop_feature) in [("feature", false), ("allele", true)] {
        let (perl, rust, paired, shared, lines) =
            expected_counts(&reference_tuples, &engine_tuples, drop_feature);
        assert!(paired < perl.min(rust) && shared > 0, "{projection}");
        assert_eq!(lines.len() as i64, perl + rust - 2 * paired, "{projection}");
        for engine in ["output", "tuples"] {
            let item = entry(&report, "Planted", engine, projection);
            assert_eq!(
                counts(item),
                [perl, rust, paired, perl, rust, paired],
                "{engine} {projection}"
            );
            assert_eq!(item["residual_rows"], lines.len(), "{engine} {projection}");
            assert_eq!(item["shared_keys"], shared, "{engine} {projection}");
            assert_eq!(
                residual_lines(&residuals, "Planted", engine, projection),
                lines,
                "{engine} {projection}"
            );
        }
    }
}

#[test]
fn summary_prints_each_class_count_and_marks_none_and_above_population() {
    let scratch = Scratch::new("summary");
    let planted = Planted::new();
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let mut changed = rows.clone();
    changed[1] = rows[1].replace("MaxEntScan_diff=1.5", "MaxEntScan_diff=-2.25");
    changed[2] = rows[2].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;LOEUF=0");
    let engine = scratch.output("engine.txt", &changed);
    let table = scratch.file(
        "classes.tsv",
        CLASS_TABLE.replace("GWAS_*\tmet\t-\n", "GWAS_*\n"),
    );
    let report_path = scratch.path("report.json");
    let plugins: Vec<&str> = PLANTED.iter().map(|(name, _)| *name).collect();
    let (code, lines, error) = planted.run_main(&[
        "--reference",
        &reference,
        "--engine",
        &format!("engine={engine}"),
        "--report",
        &report_path,
        "--work-dir",
        &scratch.path("work"),
        "--sort-memory",
        "16M",
        "--class-table",
        &table,
        "--classes",
        "planted_two_numbers, planted_reference_value,order_only",
        "--class-population",
        "planted_two_numbers=0",
        "--class-population",
        "planted_reference_value=5",
        "--plugins",
        &plugins.join(", "),
        "--engine-binary",
        &format!("engine={}", "ab".repeat(32)),
        "--engine-commit",
        "engine=0123abc",
    ]);
    assert_eq!((code, error.as_str()), (0, ""));
    let printed = |line: &str| lines.iter().any(|held| held == line);
    assert!(printed(
        "keys of the reference that no source claims: Custom (1), Custom_SIG (1)"
    ));
    assert!(printed(
        "keys of the engine engine that no source claims: Custom (1), Custom_SIG (1)"
    ));
    assert!(printed(
        "MaxEntScan\tengine\tfeature\tclass planted_two_numbers set aside 1 (above its reachable population of 0) [met, -]"
    ));
    assert!(printed(
        "GWAS\tengine\tfeature\tclass planted_reference_value set aside 0 (none)"
    ));
    assert!(printed(
        "GWAS\tengine\tallele\tclass planted_reference_value set aside 0 (none)"
    ));
    assert!(printed(
        "LOEUF\tengine\tfeature\tclass planted_engine_zero is not applied and would set aside 1: raw counts only [met, -]"
    ));
    assert!(printed(
        "LOEUF\tengine\tfeature\tperl=1 rust=2 intersection=1 text_only_differences=0 residual_rows=1 shared_keys=0"
    ));
    assert!(
        !lines
            .iter()
            .any(|line| line.contains("Planted\t") && line.contains("class")),
        "{lines:?}"
    );
    let report: Value = serde_json::from_slice(&fs::read(&report_path).unwrap()).unwrap();
    let populated: Vec<&Value> = report["classes"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|item| item.get("population").is_some())
        .collect();
    assert_eq!(
        json!(populated),
        json!([
            {
                "class": "planted_two_numbers",
                "plugins": {"MaxEntScan": ["MaxEntScan_*"]},
                "applied": true,
                "further_columns": ["met", "-"],
                "population": 0,
            },
            {
                "class": "planted_reference_value",
                "plugins": {"GWAS": ["GWAS_*"]},
                "applied": true,
                "further_columns": [],
                "population": 5,
            },
        ])
    );
    assert_eq!(
        feature(&report, "MaxEntScan")["masked"],
        json!({"planted_two_numbers": 1})
    );
    assert_eq!(report["plugins"], json!(plugins));
    let described = &report["engines"][0];
    assert_eq!(described["id"], "engine");
    assert_eq!(described["binary_sha256"], json!("ab".repeat(32)));
    assert_eq!(described["commit"], "0123abc");
    assert!(report["reference"].get("binary_sha256").is_none());
}

#[test]
fn summary_needs_no_standard_output() {
    let scratch = Scratch::new("no-stdout");
    let reference = scratch.output("reference.txt", &reference_rows());
    let no_class = scratch.file("no-class.tsv", "");
    let report_path = scratch.path("report.json");
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-plugin-pairs"));
    command
        .arg("score")
        .args(["--reference", &reference])
        .args(["--engine", &format!("engine={reference}")])
        .args(["--plugins", "", "--custom-source", "Planted=PLANT_"])
        .args(["--class-table", &no_class])
        .args(["--work-dir", &scratch.path("work")])
        .args(["--report", &report_path])
        .stdin(Stdio::null());
    // SAFETY: `close` may be called between `fork` and `exec`; the program starts with no
    // open standard output.
    unsafe {
        command.pre_exec(|| {
            libc::close(1);
            Ok(())
        });
    }
    let done = command.output().unwrap();
    assert_eq!(
        (
            done.status.code(),
            String::from_utf8_lossy(&done.stderr).as_ref()
        ),
        (Some(0), "")
    );
    assert!(done.stdout.is_empty());
    let report: Value = serde_json::from_slice(&fs::read(&report_path).unwrap()).unwrap();
    assert_eq!(counts(feature(&report, "Planted")), [4; 6]);
}
