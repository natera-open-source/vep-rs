// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `score` program on plugin values at a record the reference dropped: a grid of
//! made-up records that says where the engine's values are set aside and where they stay
//! residuals, with the input records and without them, read alike from a plain and a gzip
//! reference and on any number of threads, and held to the key comparison's own reading
//! of the same two outputs.

mod support;

use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use support::{
    feature, feature_residuals, gzip, on_any_thread_count, row_with, Planted, PlantedRegistry,
    PlantedRules, Run, Scratch,
};
use vep_plugin_pairs::rules::PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD;

const CLASS: &str = PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD;
const PLUGINS: [(&str, &[&str]); 2] = [
    ("DosageSensitivity", &["pHaplo", "pTriplo"]),
    ("UTRAnnotator", &["Existing_uORFs", "5UTR_consequence"]),
];

/// A written row: its first column, Location, Allele, Feature and Extra.
type Written = [&'static str; 5];

fn text(row: &Written) -> String {
    let [uploaded, location, allele, feature, extra] = *row;
    row_with(location, allele, feature, extra, "-", uploaded)
}

/// One cell of the grid: the reference's rows, the engine's rows, and whether the
/// engine's values there are set aside with the input records and without them.
struct Cell {
    name: &'static str,
    reference: &'static [Written],
    engine: &'static [Written],
    with_records: bool,
    without: bool,
}

const GRID: [Cell; 11] = [
    Cell {
        name: "a structural record over the size limit",
        reference: &[],
        engine: &[
            [
                "sv_over",
                "21:7000001-27000000",
                "deletion",
                "ENST1",
                "pHaplo=0.9;pTriplo=0.2",
            ],
            [
                "sv_over",
                "21:7000001-27000000",
                "deletion",
                "ENST2",
                "pHaplo=0.5;Existing_uORFs=2",
            ],
        ],
        with_records: true,
        without: true,
    },
    Cell {
        name: "a structural record under the limit",
        reference: &[],
        engine: &[[
            "sv_under",
            "21:30000001-30000900",
            "duplication",
            "ENST1",
            "pHaplo=0.7",
        ]],
        with_records: true,
        without: true,
    },
    Cell {
        name: "a sequence allele beside the symbolic one",
        reference: &[],
        engine: &[
            [
                "sv_mixed",
                "21:37000001-37000002",
                "deletion",
                "ENST1",
                "pHaplo=0.11",
            ],
            [
                "sv_mixed",
                "21:37000001-37000002",
                "T",
                "ENST1",
                "pHaplo=0.12",
            ],
        ],
        with_records: true,
        without: true,
    },
    Cell {
        name: "every row a sequence allele",
        reference: &[],
        engine: &[["snv_only", "21:31000000", "T", "ENST1", "pHaplo=0.6"]],
        with_records: false,
        without: false,
    },
    Cell {
        name: "the allele tandem_repeat",
        reference: &[],
        engine: &[[
            "tr_1",
            "21:32000001-32000050",
            "tandem_repeat",
            "ENST1",
            "pHaplo=0.4",
        ]],
        with_records: false,
        without: false,
    },
    Cell {
        name: "the reference wrote a row at the Location",
        reference: &[[
            "sv_kept",
            "21:33000001-33000500",
            "deletion",
            "ENST1",
            "IMPACT=MODIFIER",
        ]],
        engine: &[[
            "sv_kept",
            "21:33000001-33000500",
            "deletion",
            "ENST1",
            "pHaplo=0.3",
        ]],
        with_records: false,
        without: false,
    },
    Cell {
        name: "a reference row carries the record's name elsewhere",
        reference: &[[
            "sv_moved",
            "21:34000001-34000500",
            "deletion",
            "ENST1",
            "IMPACT=MODIFIER",
        ]],
        engine: &[[
            "sv_moved",
            "21:34000000-34000500",
            "deletion",
            "ENST1",
            "pHaplo=0.8",
        ]],
        with_records: false,
        without: false,
    },
    Cell {
        name: "an engine row at the Location names no record",
        reference: &[],
        engine: &[
            [
                "sv_named",
                "21:35000001-35900000",
                "deletion",
                "ENST1",
                "pHaplo=0.2",
            ],
            [
                ".",
                "21:35000001-35900000",
                "deletion",
                "ENST2",
                "pHaplo=0.25",
            ],
        ],
        with_records: false,
        without: false,
    },
    Cell {
        name: "a field the class does not read",
        reference: &[],
        engine: &[[
            "sv_field",
            "21:36000001-36900000",
            "deletion",
            "ENST1",
            "5UTR_consequence=gained",
        ]],
        with_records: false,
        without: false,
    },
    Cell {
        name: "a gVCF variant site under an identifier",
        reference: &[],
        engine: &[
            ["gv_1", "21:38000000", "T", "ENST1", "pHaplo=0.13"],
            ["gv_1", "21:38000000", "<NON_REF>", "ENST1", "pHaplo=0.14"],
        ],
        with_records: false,
        without: true,
    },
    Cell {
        name: "a gVCF variant site under a generated name",
        reference: &[],
        engine: &[
            [
                "21_39000000_A/T/<NON_REF>",
                "21:39000000",
                "T",
                "ENST1",
                "pHaplo=0.15",
            ],
            [
                "21_39000000_A/T/<NON_REF>",
                "21:39000000",
                "<NON_REF>",
                "ENST1",
                "pHaplo=0.16",
            ],
        ],
        with_records: false,
        without: false,
    },
];
const INPUT: [[&str; 5]; 5] = [
    ["21", "7000000", "sv_over", "N", "<DEL>"],
    ["21", "30000000", "sv_under", "N", "<DUP>"],
    ["21", "37000000", "sv_mixed", "N", "<DEL>,NT"],
    ["21", "38000000", "gv_1", "A", "T,<NON_REF>"],
    ["21", "39000000", ".", "A", "T,<NON_REF>"],
];

/// The two plugins' fields and the class's name under a predicate that takes no pair: the
/// program decides the class without one.
fn planted() -> Planted {
    Planted {
        registry: PlantedRegistry::new(&PLUGINS),
        rules: PlantedRules::none().class(CLASS, &[], |_| false),
    }
}

fn class_table() -> String {
    format!("{CLASS}\tDosageSensitivity;UTRAnnotator\tpHaplo,pTriplo;Existing_uORFs\tmet\t-\n")
}

/// The run of the grid: the class's table and the two plugins.
fn grid_run() -> Run {
    Run::new()
        .table(&class_table())
        .plugins(&["DosageSensitivity", "UTRAnnotator"])
}

/// The grid's reference and engine outputs, each beside a row both sides write, the
/// reference with a structural row the engine lacks, and the input VCF of their run.
fn outputs(scratch: &Scratch) -> (String, String, String) {
    let kept = text(&["snv_1", "21:100", "T", "ENST9", "pHaplo=0.01"]);
    let mut reference = vec![
        kept.clone(),
        text(&[
            "sv_ref",
            "21:40000001-40000500",
            "deletion",
            "ENST1",
            "pHaplo=0.1",
        ]),
    ];
    reference.extend(GRID.iter().flat_map(|cell| cell.reference.iter().map(text)));
    let mut engine = vec![kept];
    engine.extend(GRID.iter().flat_map(|cell| cell.engine.iter().map(text)));
    (
        scratch.output("reference.txt", &reference),
        scratch.output("engine.txt", &engine),
        scratch.vcf("input.vcf", &INPUT),
    )
}

/// The Location, Allele and Feature of each engine row of a cell.
fn rows_of(cell: &Cell) -> BTreeSet<[String; 3]> {
    cell.engine
        .iter()
        .map(|row| [row[1], row[2], row[3]].map(str::to_string))
        .collect()
}

/// Per cell of the grid, whether no value of its engine rows is left in the residual
/// files.
fn cells_set_aside(residuals: &str) -> Vec<bool> {
    let left: BTreeSet<[String; 3]> = PLUGINS
        .iter()
        .flat_map(|(source, _)| feature_residuals(residuals, source))
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            [cells[1], cells[2], cells[3]].map(str::to_string)
        })
        .collect();
    GRID.iter()
        .map(|cell| rows_of(cell).is_disjoint(&left))
        .collect()
}

#[test]
fn plugin_value_on_a_record_the_reference_dropped_is_set_aside_on_the_grid_and_nothing_wider() {
    for with_input in [true, false] {
        let scratch = Scratch::new("dropped-grid");
        let planted = planted();
        let (reference, engine, run_input) = outputs(&scratch);
        let residuals = scratch.path("residuals");
        let run = grid_run().residuals(&residuals);
        let run = if with_input {
            run.run_input(&run_input)
        } else {
            run
        };
        let report = planted.score(&scratch, &reference, &[("engine", &engine)], &run);
        let expected: Vec<bool> = GRID
            .iter()
            .map(|cell| {
                if with_input {
                    cell.with_records
                } else {
                    cell.without
                }
            })
            .collect();
        let names: Vec<&str> = GRID.iter().map(|cell| cell.name).collect();
        assert_eq!(
            cells_set_aside(&residuals),
            expected,
            "input {with_input}: {names:?}"
        );
        let aside = if with_input { 7 } else { 9 };
        assert_eq!(
            feature(&report, "DosageSensitivity")["masked"],
            json!({CLASS: aside - 1}),
            "input {with_input}"
        );
        assert_eq!(
            feature(&report, "UTRAnnotator")["masked"],
            json!({CLASS: 1}),
            "input {with_input}"
        );
        assert!(feature_residuals(&residuals, "DosageSensitivity")
            .iter()
            .any(|line| line == "reference\t21:40000001-40000500\tdeletion\tENST1\tpHaplo\t0.1"));
    }
}

/// The built key comparison: the file `VEP_OUTPUTS_COMPARE` names, else the binary under
/// the cargo target directory or under the tool's own.
fn built_key_comparison() -> Option<PathBuf> {
    const NAME: &str = "vep-outputs-compare";
    if let Some(named) = std::env::var_os("VEP_OUTPUTS_COMPARE").filter(|named| !named.is_empty()) {
        return Some(PathBuf::from(named));
    }
    let mut targets = vec![Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join(NAME)
        .join("target")];
    if let Some(shared) = std::env::var_os("CARGO_TARGET_DIR").filter(|named| !named.is_empty()) {
        targets.insert(0, PathBuf::from(shared));
    }
    targets
        .iter()
        .flat_map(|target| ["release", "debug"].map(|profile| target.join(profile).join(NAME)))
        .find(|binary| binary.is_file())
}

/// The Location, Allele and Feature of every row of a discordant-rows file after its
/// header line.
fn discordant_rows(path: &Path) -> BTreeSet<[String; 3]> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            [cells[2], cells[3], cells[4]].map(str::to_string)
        })
        .collect()
}

#[test]
fn plugin_value_class_takes_the_rows_the_key_comparison_sets_aside_as_a_skipped_record() {
    let Some(key_comparison) = built_key_comparison() else {
        println!(
            "no vep-outputs-compare binary: set VEP_OUTPUTS_COMPARE or build tools/vep-outputs-compare; the grid is not held to the key comparison"
        );
        return;
    };
    for with_input in [true, false] {
        let scratch = Scratch::new("dropped-key-comparison");
        let planted = planted();
        let (reference, engine, run_input) = outputs(&scratch);
        let residuals = scratch.path("residuals");
        let run = grid_run().residuals(&residuals);
        let run = if with_input {
            run.run_input(&run_input)
        } else {
            run
        };
        planted.score(&scratch, &reference, &[("engine", &engine)], &run);
        let scored = cells_set_aside(&residuals);

        for (side, written) in [("perl", &reference), ("rust", &engine)] {
            scratch.file(&format!("{side}/output.txt"), fs::read(written).unwrap());
        }
        let mut line = vec![
            "--perl-dir".to_string(),
            scratch.path("perl"),
            "--rust-dir".to_string(),
            scratch.path("rust"),
            "--report-dir".to_string(),
            scratch.path("report"),
            "--assembly".to_string(),
            "GRCh38".to_string(),
            "--glob".to_string(),
            "*.txt".to_string(),
        ];
        if with_input {
            line.extend(["--input-vcf".to_string(), run_input.clone()]);
        }
        let done = Command::new(&key_comparison)
            .args(&line)
            .env("VEP_OUTPUTS_COMPARE_THREADS", "2")
            .output()
            .unwrap();
        assert_eq!(
            done.status.code(),
            Some(0),
            "{}",
            String::from_utf8_lossy(&done.stderr)
        );
        let summary: Value =
            serde_json::from_slice(&fs::read(scratch.path("report/summary.json")).unwrap())
                .unwrap();
        let excluded = &summary["aggregate"]["adjusted"]["excluded_categories"];
        let one_sided = discordant_rows(Path::new(&scratch.path("report/discordant.tsv")));
        let left_open = discordant_rows(Path::new(&scratch.path("report/discordant_open.tsv")));
        let rows: Vec<BTreeSet<[String; 3]>> = GRID.iter().map(rows_of).collect();
        // A gVCF variant site its generated name or the input shows is set aside too,
        // under its own class.
        let mut other_class = rows[10].clone();
        if with_input {
            other_class.extend(rows[9].iter().cloned());
        }
        let mut compared: Vec<bool> = rows
            .iter()
            .map(|rows| {
                !rows.is_empty()
                    && rows.iter().all(|row| {
                        one_sided.contains(row)
                            && !left_open.contains(row)
                            && !other_class.contains(row)
                    })
            })
            .collect();
        // The row of the field the class does not read is a dropped record's row all the
        // same.
        compared[8] = false;
        assert_eq!(scored, compared, "input {with_input}");
        let taken: usize = rows
            .iter()
            .zip(&compared)
            .filter(|(_, aside)| **aside)
            .map(|(rows, _)| rows.len())
            .sum();
        assert_eq!(
            excluded["reference_skipped_record"]["excluded"],
            taken + 1,
            "input {with_input}"
        );
        assert_eq!(
            excluded["gvcf_variant_site_dropped"]["excluded"],
            other_class.len(),
            "input {with_input}"
        );
    }
}

#[test]
fn plugin_value_on_a_dropped_record_is_left_for_a_tuple_file_engine_and_read_alike_plain_or_gzip() {
    let scratch = Scratch::new("dropped-stored");
    let planted = planted();
    let (reference, engine, run_input) = outputs(&scratch);
    // A tuple file states no rows, so none of its values is at a dropped record's
    // Location.
    let listed = scratch.tuples(
        "engine.tsv",
        &[["21:7000001-27000000", "deletion", "ENST1", "pHaplo", "0.9"]],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &listed)], &grid_run());
    assert_eq!(
        feature(&report, "DosageSensitivity")["masked"],
        json!({CLASS: 0})
    );
    let stored = scratch.file("reference.txt.gz", gzip(&fs::read(&reference).unwrap()));
    let mut reports = Vec::new();
    for given in [&reference, &stored] {
        let (report, files) = on_any_thread_count(
            &scratch,
            given,
            &[("engine", &engine)],
            &grid_run().run_input(&run_input),
            &planted.registry,
            &planted.rules,
        );
        reports.push((
            feature(&report, "DosageSensitivity")["masked"].clone(),
            files,
        ));
    }
    assert_eq!(reports[0], reports[1]);
    assert_eq!(reports[0].0, json!({CLASS: 6}));
}
