// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The passes of the `score` program under a planted registry and planted classes: the
//! same report and residual files on any number of threads, the values left to the rules,
//! the refusal of an input under the path as it was given, and the command line's
//! `--pairing-tool` and `--pairing-threads`.

mod support;

use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs;
use std::process::Command;
use support::{
    command, counts, feature, gzip, on_any_thread_count, reference_rows, refusal, row, row_with,
    sha256, Planted, PlantedRegistry, PlantedRules, Run, Scratch, ERROR, HEADER,
};
use vep_plugin_pairs::registry::Registry;
use vep_plugin_pairs::rules::decode_extra;

/// A sequence of numbers that is the same in every run (SplitMix64 from a fixed seed).
struct Drawn(u64);

impl Drawn {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut mixed = self.0;
        mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        mixed ^ (mixed >> 31)
    }

    fn choice<'a>(&mut self, items: &[&'a str]) -> &'a str {
        items[(self.next() % items.len() as u64) as usize]
    }

    /// True four times in ten.
    fn four_in_ten(&mut self) -> bool {
        self.next() % 10 < 4
    }
}

const DRAWN_LOCATIONS: [&str; 8] = [
    "1:5",
    "1:50",
    "1:5-06",
    "1:5-6",
    "10:5",
    "MT:+5",
    "21:\u{b2}",
    " 2:5 ",
];
const DRAWN_ALLELES: [&str; 4] = ["A", "-", "AT", " T"];
const DRAWN_FEATURES: [&str; 4] = ["ENST1", "ENST10", "ENST2 ", "-"];
/// The keys a drawn row can write, each with the values it is drawn from: fields of the
/// planted plugins, a core key, a field of the custom source, a key no source claims and a
/// key written without a value.
const DRAWN_FIELDS: [(&str, &[&str]); 12] = [
    ("PLANT_score", &["0.1", "0.10", "7", "7.0000001", "-", ""]),
    ("PLANT_label", &["a", "b", "a%3Bb%3Dc", "x\u{1}y", "%3d"]),
    ("LOEUF", &["0", "0.31"]),
    ("MaxEntScan_diff", &["1.5", "-2", "n"]),
    ("GWAS_pmid", &["17", "18"]),
    ("NearestExonJB", &["a,b", "b,a", "a"]),
    ("5UTR_consequence", &["gained", "lost"]),
    ("REVEL", &["0.2", "0.4"]),
    ("IMPACT", &["LOW"]),
    ("Custom_SIG", &["p", "q"]),
    ("Stray", &["1"]),
    ("Bare", &["yes"]),
];

/// `count` drawn rows and the tuples a tuple file states for them: each row a Location, an
/// Allele and a Feature with about four in ten of the keys of `DRAWN_FIELDS`.
fn drawn_side(drawn: &mut Drawn, count: usize) -> (Vec<String>, Vec<[String; 5]>) {
    let (mut rows, mut listed) = (Vec::new(), Vec::new());
    for _ in 0..count {
        let key = [
            drawn.choice(&DRAWN_LOCATIONS),
            drawn.choice(&DRAWN_ALLELES),
            drawn.choice(&DRAWN_FEATURES),
        ];
        let mut pairs: Vec<(&str, &str)> = Vec::new();
        for (name, values) in DRAWN_FIELDS {
            if drawn.four_in_ten() {
                pairs.push((name, drawn.choice(values)));
            }
        }
        let extra: Vec<String> = pairs
            .iter()
            .map(|(name, value)| match *name {
                "Bare" => "Bare".to_string(),
                _ => format!("{name}={value}"),
            })
            .collect();
        let extra = match extra.is_empty() {
            true => "-".to_string(),
            false => extra.join(";"),
        };
        let cdna = drawn.choice(&["1", "2"]);
        let uploaded = drawn.choice(&["v", "1_5_A/T"]);
        rows.push(row_with(key[0], key[1], key[2], &extra, cdna, uploaded));
        for (name, value) in pairs {
            if name != "Bare" {
                listed.push([
                    key[0].to_string(),
                    key[1].to_string(),
                    key[2].to_string(),
                    name.to_string(),
                    decode_extra(value),
                ]);
            }
        }
    }
    (rows, listed)
}

#[test]
fn thread_count_changes_no_byte_of_the_report_and_residual_files_over_drawn_outputs() {
    let scratch = Scratch::new("drawn-threads");
    let mut drawn = Drawn(1162);
    let reference = scratch.output("reference.txt", &drawn_side(&mut drawn, 300).0);
    let engine_rows = drawn_side(&mut drawn, 300).0;
    // A default-format output with a blank line and a comment between its rows, plain and
    // as two gzip members with zero bytes after them.
    let plain_text = format!(
        "{HEADER}{}\r\n##a comment\r\n{}",
        engine_rows[..150].concat(),
        engine_rows[150..].concat()
    );
    let plain = scratch.file("engine.txt", &plain_text);
    let mut packed_bytes = gzip(&plain_text.as_bytes()[..900]);
    packed_bytes.extend(gzip(&plain_text.as_bytes()[900..]));
    packed_bytes.extend([0u8; 40]);
    let packed = scratch.file("engine.txt.gz", &packed_bytes);
    let listed_tuples = drawn_side(&mut drawn, 300).1;
    let lines: Vec<[&str; 5]> = listed_tuples
        .iter()
        .map(|line| [0, 1, 2, 3, 4].map(|at| line[at].as_str()))
        .collect();
    let listed = scratch.tuples("engine.tsv", &lines);
    let engines = [
        ("output", &*plain),
        ("packed", &*packed),
        ("listed", &*listed),
    ];
    let mut consults = Vec::new();
    let mut last = Value::Null;
    for table in [true, false] {
        let planted = Planted::new();
        let run = Run::new().custom(&[("Custom", "Custom")]);
        let run = if table { run } else { run.no_table() };
        let (report, files) = on_any_thread_count(
            &scratch,
            &reference,
            &engines,
            &run,
            &planted.registry,
            &planted.rules,
        );
        let calls: Vec<String> = planted
            .calls()
            .iter()
            .map(|call| format!("{call:?}"))
            .collect();
        let per_run = calls.len() / 2;
        assert_eq!(calls[..per_run], calls[per_run..], "table {table}");
        consults.push(per_run);
        let sources = report["sources"].as_object().unwrap().len();
        assert_eq!(files.len(), 2 * engines.len() * sources, "table {table}");
        assert!(files.values().map(Vec::len).sum::<usize>() > 0);
        assert_eq!(report["engines"][2]["format"], "tuples");
        assert!(
            report["engines"][2]["tuples_outside_every_source"]
                .as_i64()
                .unwrap()
                > 0
        );
        assert_eq!(report["engines"][1]["sha256"], sha256(&packed_bytes));
        assert_eq!(
            report["engines"][0]["unclaimed_keys"],
            report["engines"][1]["unclaimed_keys"]
        );
        let unclaimed: Vec<&String> = report["engines"][0]["unclaimed_keys"]
            .as_object()
            .unwrap()
            .keys()
            .collect();
        assert_eq!(unclaimed, ["Stray"]);
        for item in report["scores"].as_array().unwrap() {
            if item["source"] == "Planted" {
                assert!(item["shared_keys"].as_i64().unwrap() > 0, "{item}");
            }
        }
        last = report;
    }
    assert!(consults[0] > 0, "{consults:?}");
    assert_eq!(consults[1], 0);
    let scores = last["scores"].as_array().unwrap();
    assert!(!scores.is_empty());
    assert!(scores.iter().all(|item| item["masked"] == json!({})));
}

#[test]
fn missing_extra_and_changed_values_are_left_to_the_rules_on_any_thread_count() {
    let scratch = Scratch::new("left-to-the-rules");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=0.5;PLANT_label=kept"),
            row(
                "1:200",
                "T",
                "ENST1",
                "PLANT_score=0.25;PLANT_label=missing",
            ),
            row("1:300", "T", "ENST1", "PLANT_score=1;PLANT_label=before"),
            row("1:400", "T", "ENST1", "PLANT_score=0.1"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row(
                "1:100",
                "T",
                "ENST1",
                "PLANT_score=0.5;PLANT_label=kept;PlantedDynamic_a=extra",
            ),
            row("1:200", "T", "ENST1", "PLANT_score=0.25"),
            row("1:300", "T", "ENST1", "PLANT_score=1;PLANT_label=after"),
            row("1:400", "T", "ENST1", "PLANT_score=0.10"),
        ],
    );
    let (report, files) = on_any_thread_count(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new().no_table(),
        &planted.registry,
        &planted.rules,
    );
    let item = feature(&report, "Planted");
    assert_eq!(counts(item), [7, 6, 5, 7, 6, 5]);
    assert_eq!(item["text_only_differences"], 1);
    assert_eq!(item["residual_rows"], 3);
    assert_eq!(
        String::from_utf8_lossy(&files["Planted.engine.feature.tsv"]),
        "reference\t1:200\tT\tENST1\tPLANT_label\tmissing\n\
         engine\t1:300\tT\tENST1\tPLANT_label\tafter\n\
         reference\t1:300\tT\tENST1\tPLANT_label\tbefore\n"
    );
    assert_eq!(
        counts(feature(&report, "PlantedDynamic")),
        [0, 1, 0, 0, 1, 0]
    );
    assert_eq!(
        files["PlantedDynamic.engine.feature.tsv"],
        b"engine\t1:100\tT\tENST1\tPlantedDynamic_a\textra\n"
    );
    assert_eq!(
        files["PlantedDynamic.engine.allele.tsv"],
        b"engine\t1:100\tT\tPlantedDynamic_a\textra\n"
    );
}

const NARROW: &str = "line 3 has 3 columns; a default-format row has 14, a tuple file line 5";

#[test]
fn refusal_of_an_input_names_the_path_as_given_and_leaves_nothing_behind() {
    let scratch = Scratch::new("input-refused");
    let rows = reference_rows();
    let reference = scratch.output("reference.txt", &rows);
    let narrow_text = format!("{HEADER}{}1:100\tT\tENST1\n", rows[0]);
    let narrow = scratch.file("narrow.txt", &narrow_text);
    let odd = scratch.file(
        "odd.txt",
        format!("{HEADER}{}#Uploaded_variation\tLocation\n", rows[0]),
    );
    // A file name can carry a form feed and a line separator.
    let boundaries = scratch.file("nar\u{c}row\u{2028}.txt", &narrow_text);
    let planted = Planted::new();
    for (engine, message) in [
        (&narrow, NARROW),
        (&odd, "line 3, the column header, names 2 columns, not 14"),
        (&boundaries, NARROW),
    ] {
        let scored = planted.try_score(
            &scratch,
            &reference,
            &[("engine", &reference), ("later", engine)],
            &Run::new(),
        );
        assert_eq!(refusal(scored), format!("{engine}: {message}"));
        assert!(scratch.0.join("work").is_dir());
        assert!(scratch.names("work").is_empty(), "{engine}");
    }

    // The command line names a file in the directory it runs in by its name alone, one
    // that reads as an option included.
    scratch.file("-h", &narrow_text);
    for name in ["nar\u{c}row\u{2028}.txt", "-h"] {
        let done = Command::new(env!("CARGO_BIN_EXE_vep-plugin-pairs"))
            .current_dir(&scratch.0)
            .env("LC_ALL", "C")
            .args(["score", "--reference", "reference.txt"])
            .args(["--engine", "engine=reference.txt"])
            .args(["--engine", &format!("later={name}")])
            .args(["--plugins", "LOEUF", "--work-dir", "work-named"])
            .args(["--report", "report.json"])
            .output()
            .unwrap();
        assert_eq!(done.status.code(), Some(2), "{name:?}");
        assert_eq!(String::from_utf8_lossy(&done.stdout), "", "{name:?}");
        assert_eq!(
            String::from_utf8_lossy(&done.stderr),
            format!("{ERROR}{name}: {NARROW}\n")
        );
        assert!(!scratch.exists("report.json"));
        assert!(scratch.names("work-named").is_empty(), "{name:?}");
    }

    // A registry that refuses a field of the run ends it with the registry's own words.
    let rows_of = PlantedRegistry::planted();
    let refusing = Planted {
        registry: PlantedRegistry::planted().resolving(move |field, plugins, flags| {
            if field == "PLANT_label" {
                return Err("the run's flags make PLANT_label a key of both engines".to_string());
            }
            rows_of.plugin_of_field(field, plugins, flags)
        }),
        rules: PlantedRules::planted(),
    };
    let scored = refusing.try_score(
        &scratch,
        &reference,
        &[("engine", &reference), ("later", &reference)],
        &Run::new(),
    );
    assert_eq!(
        refusal(scored),
        "the run's flags make PLANT_label a key of both engines"
    );
    assert!(scratch.names("work").is_empty());
}

#[test]
fn command_line_takes_a_pairing_tool_and_a_thread_count_and_writes_the_same_bytes() {
    let scratch = Scratch::new("pairing-options");
    let reference = scratch.output(
        "reference.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "Custom=id1;Custom_SIG=p;Other=1",
        )],
    );
    let engine = scratch.output(
        "engine.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "Custom=id1;Custom_SIG=q;Other=1",
        )],
    );
    let tool = env!("CARGO_BIN_EXE_vep-plugin-pairs");
    let mut written = Vec::new();
    for (name, more) in [
        ("own", vec![]),
        (
            "tool",
            vec!["--pairing-tool", tool, "--pairing-threads", "2"],
        ),
    ] {
        let mut line = vec![
            "--reference".to_string(),
            reference.clone(),
            "--engine".to_string(),
            format!("engine={engine}"),
            "--custom-source".to_string(),
            "Custom=Custom".to_string(),
            "--work-dir".to_string(),
            scratch.path("work"),
            "--report".to_string(),
            scratch.path(&format!("{name}.json")),
            "--residuals".to_string(),
            scratch.path(name),
        ];
        line.extend(more.iter().map(|word| word.to_string()));
        let line: Vec<&str> = line.iter().map(String::as_str).collect();
        let done = command(&line, &[]);
        assert_eq!(done.status.code(), Some(0), "{name}");
        assert_eq!(String::from_utf8_lossy(&done.stderr), "", "{name}");
        let files: BTreeMap<String, Vec<u8>> = scratch
            .names(name)
            .into_iter()
            .map(|file| {
                let bytes = fs::read(scratch.0.join(name).join(&file)).unwrap();
                (file, bytes)
            })
            .collect();
        let report = fs::read(scratch.path(&format!("{name}.json"))).unwrap();
        written.push((String::from_utf8(done.stdout).unwrap(), report, files));
        assert!(scratch.0.join("work").is_dir());
        assert!(scratch.names("work").is_empty(), "{name}");
    }
    assert_eq!(written[0], written[1]);
    assert!(written[1]
        .0
        .contains("keys of the reference that no source claims: Other (1)\n"));
    assert_eq!(
        written[1].2["Custom.engine.feature.tsv"],
        b"reference\t1:100\tT\tENST1\tCustom_SIG\tp\nengine\t1:100\tT\tENST1\tCustom_SIG\tq\n"
    );
    let nowhere = scratch.path("no-such-tool");
    let missing = command(
        &[
            "--reference",
            &reference,
            "--engine",
            &format!("engine={engine}"),
            "--custom-source",
            "Custom=Custom",
            "--report",
            &scratch.path("report.json"),
            "--pairing-tool",
            &nowhere,
        ],
        &[],
    );
    assert_eq!(missing.status.code(), Some(2));
    assert_eq!(String::from_utf8_lossy(&missing.stdout), "");
    assert_eq!(
        String::from_utf8_lossy(&missing.stderr),
        format!("{ERROR}not a file: {nowhere}\n")
    );
}
