// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the `score` program hands a planted class beside a row's two values: the GWAS
//! catalog and the candidate cache, each to the classes of its own plugin alone, the
//! classes the run applies, an empty value a row writes, and the rows of the record that
//! the run's input names for a row under an identifier.

mod support;

use serde_json::json;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::sync::{Arc, Mutex};
use support::{
    counts, entry, feature, feature_residuals, gzip, raw_counts, residual_lines, row, row_with,
    sha256, Call, Planted, PlantedRegistry, PlantedRules, Run, Scratch, ERROR, HEADER, VCF_HEADER,
};
use vep_plugin_pairs::rules::{
    BuiltIn, Entry, Rules, ABSENT, DOWNSTREAM_START_PAST_OWN_CODON,
    DOWNSTREAM_UTR_FRAME_PAST_CODING_END, PLUGIN_HGVS_SHIFTED_CDS_COORDINATES,
};

/// The addresses of the reads a planted class was handed, one per consult that carried one.
type Addresses = Arc<Mutex<Vec<usize>>>;

/// One command line of `score` over a reference and an engine under a class table.
struct CommandLine<'a> {
    scratch: &'a Scratch,
    reference: &'a str,
    engine: &'a str,
    plugins: &'a str,
    table: &'a str,
}

impl CommandLine<'_> {
    /// The exit status, the lines of standard output and standard error of the command
    /// line with `more` after its options; its residual files go to `name` and its own
    /// scratch to `work-<name>`.
    fn run(&self, planted: &Planted, name: &str, more: &[&str]) -> (i32, Vec<String>, String) {
        let engine = format!("engine={}", self.engine);
        let work = self.scratch.path(&format!("work-{name}"));
        let report = self.scratch.path("report.json");
        let residuals = self.scratch.path(name);
        let mut line = vec![
            "--reference",
            self.reference,
            "--engine",
            &engine,
            "--plugins",
            self.plugins,
            "--class-table",
            self.table,
            "--work-dir",
            &work,
            "--sort-memory",
            "16M",
            "--report",
            &report,
            "--residuals",
            &residuals,
            "--pairing-threads",
            "2",
        ];
        line.extend(more);
        planted.run_main(&line)
    }
}

fn of_class<'a>(calls: &'a [Call], klass: &str) -> Vec<&'a Call> {
    calls.iter().filter(|call| call.klass == klass).collect()
}

const CATALOG_TABLE: &str =
    "planted_catalog_place\tGWAS\tGWAS_*\tmet\t-\nplanted_engine_zero\tLOEUF\tLOEUF\tmet\t-\n";
/// A GWAS catalog associations file of one row, at 1:100 under the risk allele T.
const CATALOG: &str =
    "CHR_ID\tCHR_POS\tSTRONGEST SNP-RISK ALLELE\tSNPS\tSTUDY ACCESSION\n1\t100\trs1-T\trs1\tGCST1\n";
const CATALOG_NOTICE: &str =
    "no --gwas-catalog: a class of GWAS that reads the catalog sets nothing aside";

/// The planted registry and classes with `planted_catalog_place`, which accepts the row at
/// the place of `CATALOG` where it is handed a catalog, and the addresses of the catalogs
/// that class was handed.
fn plant_catalog() -> (Planted, Addresses) {
    let handed = Addresses::default();
    let seen = Arc::clone(&handed);
    let rules = PlantedRules::planted().class("planted_catalog_place", &[], move |c| {
        let Some(pairs) = c.gwas_pairs else {
            return false;
        };
        seen.lock()
            .unwrap()
            .push(std::ptr::from_ref(pairs) as usize);
        c.e.location.as_deref() == Some("1:100") && c.e.allele == "T"
    });
    let planted = Planted {
        registry: PlantedRegistry::planted(),
        rules,
    };
    (planted, handed)
}

/// A catalog, and a reference that writes a GWAS field at the catalog's place and at
/// another where the engine writes none, beside a row on which the engine alone writes a
/// LOEUF of 0.
fn catalog_inputs(scratch: &Scratch) -> (String, String, String) {
    let catalog = scratch.file("associations.tsv", CATALOG);
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "GWAS_pmid=17"),
            row("1:200", "T", "ENST1", "GWAS_pmid=18"),
            row("1:300", "T", "ENST1", "IMPACT=LOW"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "IMPACT=LOW"),
            row("1:200", "T", "ENST1", "IMPACT=LOW"),
            row("1:300", "T", "ENST1", "LOEUF=0"),
        ],
    );
    (catalog, reference, engine)
}

#[test]
fn gwas_catalog_is_handed_to_the_classes_of_gwas_alone() {
    let scratch = Scratch::new("catalog-handed");
    let (catalog, reference, engine) = catalog_inputs(&scratch);
    let listed = scratch.tuples("engine.tsv", &[["1:300", "T", "ENST1", "LOEUF", "0"]]);
    let engines = [("engine", engine.as_str()), ("listed", listed.as_str())];
    let run = Run::new().table(CATALOG_TABLE).plugins(&["GWAS", "LOEUF"]);
    let (planted, catalogs) = plant_catalog();
    let residuals = scratch.path("without");
    let without = planted.score(
        &scratch,
        &reference,
        &engines,
        &run.clone().residuals(&residuals),
    );
    assert!(without["gwas_catalog"].is_null());
    assert!(catalogs.lock().unwrap().is_empty());
    let calls = planted.calls();
    assert!(!calls.is_empty());
    assert!(calls
        .iter()
        .all(|call| call.handed.is_empty() && call.declared.is_empty()));
    for (ident, _) in engines {
        let item = entry(&without, "GWAS", ident, "feature");
        assert_eq!(counts(item), [2, 0, 0, 2, 0, 0], "{ident}");
        assert_eq!(
            item["masked"],
            json!({"planted_catalog_place": 0}),
            "{ident}"
        );
        assert_eq!(
            residual_lines(&residuals, "GWAS", ident, "feature").len(),
            2,
            "{ident}"
        );
    }
    let (planted, catalogs) = plant_catalog();
    let residuals = scratch.path("with");
    let report = planted.score(
        &scratch,
        &reference,
        &engines,
        &run.residuals(&residuals).catalog(&catalog),
    );
    assert_eq!(
        report["gwas_catalog"],
        json!({"path": catalog, "sha256": sha256(&fs::read(&catalog).unwrap())})
    );
    for (ident, _) in engines {
        let item = entry(&report, "GWAS", ident, "feature");
        assert_eq!(counts(item), [2, 0, 0, 1, 0, 0], "{ident}");
        assert_eq!(
            item["masked"],
            json!({"planted_catalog_place": 1}),
            "{ident}"
        );
        assert_eq!(
            residual_lines(&residuals, "GWAS", ident, "feature"),
            ["reference\t1:200\tT\tENST1\tGWAS_pmid\t18"],
            "{ident}"
        );
        let by_allele = entry(&report, "GWAS", ident, "allele");
        assert_eq!(counts(by_allele), [2, 0, 0, 1, 0, 0], "{ident}");
        assert_eq!(
            by_allele["masked"],
            json!({"planted_catalog_place": 1}),
            "{ident}"
        );
        let zero = entry(&report, "LOEUF", ident, "feature");
        assert_eq!(counts(zero), [0, 1, 0, 0, 0, 0], "{ident}");
        assert_eq!(zero["masked"], json!({"planted_engine_zero": 1}), "{ident}");
    }
    let calls = planted.calls();
    let of_gwas = of_class(&calls, "planted_catalog_place");
    assert_eq!(of_gwas.len(), 4);
    assert!(of_gwas
        .iter()
        .all(|call| call.handed == ["gwas_pairs"] && call.declared.is_empty()));
    // One catalog, the same at every consult.
    let catalogs = catalogs.lock().unwrap();
    assert_eq!(catalogs.len(), 4);
    assert!(catalogs.iter().all(|address| *address == catalogs[0]));
    let of_zero = of_class(&calls, "planted_engine_zero");
    assert_eq!(of_zero.len(), 2);
    assert!(of_zero
        .iter()
        .all(|call| call.handed.is_empty() && call.declared.is_empty()));
}

#[test]
fn class_of_gwas_that_is_not_applied_is_handed_the_catalog_and_listed_as_not_applied() {
    let scratch = Scratch::new("catalog-not-applied");
    let (catalog, reference, engine) = catalog_inputs(&scratch);
    let listed = scratch.tuples("engine.tsv", &[["1:300", "T", "ENST1", "LOEUF", "0"]]);
    let engines = [("engine", engine.as_str()), ("listed", listed.as_str())];
    let (planted, _) = plant_catalog();
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &engines,
        &Run::new()
            .table(CATALOG_TABLE)
            .applied(&["planted_engine_zero"])
            .catalog(&catalog)
            .residuals(&residuals),
    );
    for (ident, _) in engines {
        for projection in ["feature", "allele"] {
            let item = entry(&report, "GWAS", ident, projection);
            assert_eq!(raw_counts(item), [2, 0, 0], "{ident} {projection}");
            assert!(
                item.as_object()
                    .unwrap()
                    .keys()
                    .all(|key| !key.starts_with("adj_")),
                "{ident} {projection}"
            );
            assert_eq!(item["masked"], json!({}), "{ident} {projection}");
            assert_eq!(
                item["not_applied"],
                json!({"planted_catalog_place": 1}),
                "{ident} {projection}"
            );
            assert_eq!(item["residual_rows"], 2, "{ident} {projection}");
        }
        assert_eq!(
            residual_lines(&residuals, "GWAS", ident, "feature"),
            [
                "reference\t1:100\tT\tENST1\tGWAS_pmid\t17",
                "reference\t1:200\tT\tENST1\tGWAS_pmid\t18",
            ],
            "{ident}"
        );
    }
    let calls = planted.calls();
    let of_gwas = of_class(&calls, "planted_catalog_place");
    assert_eq!(of_gwas.len(), 4);
    assert!(of_gwas.iter().all(|call| call.handed == ["gwas_pairs"]));
}

#[test]
fn summary_says_that_classes_of_gwas_set_nothing_aside_without_the_catalog() {
    let scratch = Scratch::new("catalog-summary");
    let (catalog, reference, engine) = catalog_inputs(&scratch);
    let table = scratch.file("classes.tsv", CATALOG_TABLE);
    let line = CommandLine {
        scratch: &scratch,
        reference: &reference,
        engine: &engine,
        plugins: "GWAS,LOEUF",
        table: &table,
    };
    let run = |name: &str, more: &[&str]| line.run(&plant_catalog().0, name, more);
    let printed = |lines: &[String]| lines.iter().any(|line| line == CATALOG_NOTICE);
    let set_aside = |keys: &str| {
        format!("GWAS\tengine\tfeature\tclass planted_catalog_place set aside {keys} [met, -]")
    };
    let (code, lines, error) = run("without", &[]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert_eq!(lines.last().map(String::as_str), Some(CATALOG_NOTICE));
    assert!(lines.contains(&set_aside("0 (none)")));
    let (code, lines, error) = run("with", &["--gwas-catalog", &catalog]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(!printed(&lines));
    assert!(lines.contains(&set_aside("1")));
    // A run that does not apply the class, or did not load the plugin, is told nothing.
    let (code, lines, error) = run("unapplied", &["--classes", "planted_engine_zero"]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(!printed(&lines));
    let (code, lines, error) = run("unloaded", &["--plugins", "LOEUF"]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(!printed(&lines));
    let nowhere = scratch.path("nowhere.tsv");
    let (code, lines, error) = run("missing", &["--gwas-catalog", &nowhere]);
    assert_eq!((code, lines.len()), (2, 0));
    assert_eq!(error, format!("{ERROR}not a file: {nowhere}\n"));
    assert!(!scratch.exists("missing") && !scratch.exists("work-missing"));
}

const DOWNSTREAM_TABLE: &str = "planted_recomputed_value\tDownstream\tDownstreamProtein\tmet\t-\nplanted_engine_zero\tLOEUF\tLOEUF\tmet\t-\n";
const DOWNSTREAM_NOTICE: &str =
    "no --json-cache: a class of Downstream that reads the candidate cache sets nothing aside";

/// A registry of Downstream and LOEUF and the planted classes with
/// `planted_recomputed_value`, which accepts a pair where it is handed the candidate
/// cache's sequences, and the addresses of the sequences that class was handed.
fn plant_downstream() -> (Planted, Addresses) {
    let handed = Addresses::default();
    let seen = Arc::clone(&handed);
    let rules = PlantedRules::planted().class("planted_recomputed_value", &[], move |c| {
        let Some(sequences) = c.downstream_sequences else {
            return false;
        };
        seen.lock()
            .unwrap()
            .push(std::ptr::from_ref(sequences) as usize);
        true
    });
    let planted = Planted {
        registry: PlantedRegistry::new(&[
            ("Downstream", &["DownstreamProtein"]),
            ("LOEUF", &["LOEUF"]),
        ]),
        rules,
    };
    (planted, handed)
}

/// A reference and an engine that write different values of a Downstream field at one
/// Location, and a candidate cache directory.
fn downstream_inputs(scratch: &Scratch) -> (String, String, String) {
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "-", "ENST1", "DownstreamProtein=KL"),
            row("1:300", "T", "ENST1", "IMPACT=LOW"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "-", "ENST1", "DownstreamProtein=RKL"),
            row("1:300", "T", "ENST1", "LOEUF=0"),
        ],
    );
    fs::create_dir_all(scratch.0.join("cache/transcripts")).unwrap();
    (reference, engine, scratch.path("cache"))
}

#[test]
fn candidate_cache_is_handed_to_the_classes_of_downstream_alone() {
    let scratch = Scratch::new("cache-handed");
    let (reference, engine, cache) = downstream_inputs(&scratch);
    let engines = [("engine", engine.as_str()), ("again", engine.as_str())];
    let run = Run::new()
        .table(DOWNSTREAM_TABLE)
        .plugins(&["Downstream", "LOEUF"]);
    let (planted, sequences) = plant_downstream();
    let without = planted.score(&scratch, &reference, &engines, &run);
    assert!(without["json_cache"].is_null());
    assert!(sequences.lock().unwrap().is_empty());
    let calls = planted.calls();
    assert!(!calls.is_empty());
    assert!(calls
        .iter()
        .all(|call| call.handed.is_empty() && call.declared.is_empty()));
    let item = feature(&without, "Downstream");
    assert_eq!(counts(item), [1, 1, 0, 1, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_recomputed_value": 0}));
    let (planted, sequences) = plant_downstream();
    let report = planted.score(&scratch, &reference, &engines, &run.cache(&cache));
    assert_eq!(report["json_cache"], json!({"path": cache}));
    for (ident, _) in engines {
        let item = entry(&report, "Downstream", ident, "feature");
        assert_eq!(counts(item), [1, 1, 0, 0, 0, 0], "{ident}");
        assert_eq!(
            item["masked"],
            json!({"planted_recomputed_value": 1}),
            "{ident}"
        );
        let zero = entry(&report, "LOEUF", ident, "feature");
        assert_eq!(counts(zero), [0, 1, 0, 0, 0, 0], "{ident}");
        assert_eq!(zero["masked"], json!({"planted_engine_zero": 1}), "{ident}");
    }
    let calls = planted.calls();
    let of_downstream = of_class(&calls, "planted_recomputed_value");
    assert!(!of_downstream.is_empty());
    // The sequences come with the classes the run applies and with nothing else.
    assert!(of_downstream
        .iter()
        .all(|call| call.handed == ["downstream_sequences"] && !call.declared.is_empty()));
    let sequences = sequences.lock().unwrap();
    assert_eq!(sequences.len(), of_downstream.len());
    assert!(sequences.iter().all(|address| *address == sequences[0]));
    let of_zero = of_class(&calls, "planted_engine_zero");
    assert!(!of_zero.is_empty());
    assert!(of_zero
        .iter()
        .all(|call| call.handed.is_empty() && call.declared.is_empty()));
}

/// A Downstream value two causes change is set aside only beside the second cause's class,
/// so a class of the plugin is handed the applied classes: every declared one where the
/// run names none, the named ones otherwise.
#[test]
fn class_of_downstream_is_handed_the_classes_the_run_applies() {
    let scratch = Scratch::new("cache-declared");
    let (reference, engine, cache) = downstream_inputs(&scratch);
    let run = Run::new()
        .table(DOWNSTREAM_TABLE)
        .plugins(&["Downstream", "LOEUF"])
        .cache(&cache);
    let declared_sets = |run: &Run| -> BTreeSet<BTreeSet<String>> {
        let (planted, _) = plant_downstream();
        planted.score(&scratch, &reference, &[("engine", &engine)], run);
        of_class(&planted.calls(), "planted_recomputed_value")
            .iter()
            .map(|call| call.declared.clone())
            .collect()
    };
    let classes = |names: &[&str]| -> BTreeSet<BTreeSet<String>> {
        BTreeSet::from([names.iter().map(|name| name.to_string()).collect()])
    };
    assert_eq!(
        declared_sets(&run),
        classes(&["planted_recomputed_value", "planted_engine_zero"])
    );
    assert_eq!(
        declared_sets(&run.clone().applied(&["planted_recomputed_value"])),
        classes(&["planted_recomputed_value"])
    );
}

#[test]
fn summary_says_that_classes_of_downstream_set_nothing_aside_without_the_candidate_cache() {
    let scratch = Scratch::new("cache-summary");
    let (reference, engine, cache) = downstream_inputs(&scratch);
    let table = scratch.file("classes.tsv", DOWNSTREAM_TABLE);
    let line = CommandLine {
        scratch: &scratch,
        reference: &reference,
        engine: &engine,
        plugins: "Downstream,LOEUF",
        table: &table,
    };
    let run = |name: &str, more: &[&str]| line.run(&plant_downstream().0, name, more);
    let printed = |lines: &[String]| lines.iter().any(|line| line == DOWNSTREAM_NOTICE);
    let (code, lines, error) = run("without", &[]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert_eq!(lines.last().map(String::as_str), Some(DOWNSTREAM_NOTICE));
    let (code, lines, error) = run("with", &["--json-cache", &cache]);
    assert_eq!((code, error.as_str()), (0, ""));
    assert!(!printed(&lines));
    assert!(lines.iter().any(|line| line
        == "Downstream\tengine\tfeature\tclass planted_recomputed_value set aside 1 [met, -]"));
    let (nowhere, root) = (scratch.path("nowhere"), scratch.0.to_str().unwrap());
    let cases = [
        (
            "missing",
            nowhere.as_str(),
            format!("not a directory: {nowhere}"),
        ),
        (
            "refused",
            root,
            format!("--json-cache: {root} has no transcripts/ directory"),
        ),
    ];
    for (name, path, message) in cases {
        let (code, lines, error) = run(name, &["--json-cache", path]);
        assert_eq!((code, lines.len()), (2, 0), "{name}");
        assert_eq!(error, format!("{ERROR}{message}\n"), "{name}");
        assert!(
            !scratch.exists(name) && !scratch.exists(&format!("work-{name}")),
            "{name}"
        );
    }
}

/// The values a class was handed, as the set of its consults' reference and engine values.
fn values_handed(calls: &[Call], classes: &[&str]) -> BTreeSet<(String, String)> {
    calls
        .iter()
        .filter(|call| classes.contains(&call.klass.as_str()))
        .map(|call| (call.expected.clone(), call.actual.clone()))
        .collect()
}

fn pairs(values: &[(&str, &str)]) -> BTreeSet<(String, String)> {
    values
        .iter()
        .map(|(expected, actual)| (expected.to_string(), actual.to_string()))
        .collect()
}

/// A row that writes `DownstreamProtein=` writes no tuple, so the other side's value is
/// left with none to pair. The two classes of Downstream read the empty value as the
/// plugin's own and are handed it as the empty text; a `-`, a row without the key and a
/// tuple file's missing tuple stay absent, and any other class is handed absent for all.
#[test]
fn classes_of_downstream_are_handed_an_empty_value_as_a_value() {
    let (first, second) = (
        DOWNSTREAM_UTR_FRAME_PAST_CODING_END,
        DOWNSTREAM_START_PAST_OWN_CODON,
    );
    let table: String = ["planted_other_class", first, second]
        .iter()
        .map(|klass| format!("{klass}\tDownstream\tDownstreamProtein\tmet\t-\n"))
        .collect();
    let planted = Planted {
        registry: PlantedRegistry::new(&[("Downstream", &["DownstreamProtein"])]),
        rules: PlantedRules::none().classes(&["planted_other_class", first, second], move |c| {
            if [c.expected, c.actual]
                .iter()
                .any(|value| ["-", ABSENT].contains(value))
            {
                return false;
            }
            (c.klass == first && c.expected.is_empty())
                || (c.klass == second && c.actual.is_empty())
        }),
    };
    let scratch = Scratch::new("empty-value");
    let on = |location: &str, extra: &str| row(location, "-", "ENST1", extra);
    let reference = scratch.output(
        "reference.txt",
        &[
            on("1:100", "DownstreamProtein="),
            on("1:200", "DownstreamProtein=KL"),
            on("1:300", "DownstreamProtein=KL"),
            on("1:400", "DownstreamProtein=-"),
            on("1:500", "DownstreamProtein="),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            on("1:100", "DownstreamProtein=K"),
            on("1:200", "DownstreamProtein="),
            on("1:300", "IMPACT=LOW"),
            on("1:400", "DownstreamProtein=R"),
            on("1:500", "DownstreamProtein="),
        ],
    );
    let listed = scratch.tuples(
        "engine.tsv",
        &[
            ["1:100", "-", "ENST1", "DownstreamProtein", "K"],
            ["1:400", "-", "ENST1", "DownstreamProtein", "R"],
        ],
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine), ("listed", &listed)],
        &Run::new()
            .residuals(&residuals)
            .table(&table)
            .plugins(&["Downstream"]),
    );
    let left = [
        "reference\t1:300\t-\tENST1\tDownstreamProtein\tKL",
        "engine\t1:400\t-\tENST1\tDownstreamProtein\tR",
    ];
    let item = feature(&report, "Downstream");
    assert_eq!(counts(item), [2, 2, 0, 1, 1, 0]);
    assert_eq!(
        item["masked"],
        json!({"planted_other_class": 0, first: 1, second: 1})
    );
    assert_eq!(feature_residuals(&residuals, "Downstream"), left);
    // A tuple file states no empty value: its side of the second row stays absent.
    let item = entry(&report, "Downstream", "listed", "feature");
    assert_eq!(counts(item), [2, 2, 0, 2, 1, 0]);
    assert_eq!(
        item["masked"],
        json!({"planted_other_class": 0, first: 1, second: 0})
    );
    assert_eq!(
        residual_lines(&residuals, "Downstream", "listed", "feature"),
        [
            "reference\t1:200\t-\tENST1\tDownstreamProtein\tKL",
            left[0],
            left[1],
        ]
    );
    let calls = planted.calls();
    assert_eq!(
        values_handed(&calls, &[first, second]),
        pairs(&[("", "K"), ("KL", ""), ("KL", ABSENT), (ABSENT, "R")])
    );
    assert_eq!(
        values_handed(&calls, &["planted_other_class"]),
        pairs(&[(ABSENT, "K"), ("KL", ABSENT), (ABSENT, "R")])
    );
}

/// `plugin_hgvs_shifted_cds_coordinates` reads a `DownstreamProtein` as the two classes of
/// Downstream do, so a side whose row writes the field with an empty value is handed to it
/// as the empty text, on either side; a class that reads no transcript is handed absent
/// for the same rows.
#[test]
fn class_of_the_shifted_cds_coordinates_is_handed_an_empty_value_as_a_value() {
    let shifted = PLUGIN_HGVS_SHIFTED_CDS_COORDINATES;
    assert!(BuiltIn.transcript_reading(shifted));
    assert!(![
        DOWNSTREAM_UTR_FRAME_PAST_CODING_END,
        DOWNSTREAM_START_PAST_OWN_CODON
    ]
    .contains(&shifted));
    let planted = Planted {
        registry: PlantedRegistry::new(&[("Downstream", &["DownstreamProtein"])]),
        rules: PlantedRules::none().classes(&["planted_other_class", shifted], move |c| {
            c.klass == shifted
                && (c.expected.is_empty() || c.actual.is_empty())
                && c.expected != ABSENT
                && c.actual != ABSENT
        }),
    };
    let scratch = Scratch::new("empty-value-shifted");
    let on = |location: &str, extra: &str| row(location, "-", "ENST1", extra);
    let reference = scratch.output(
        "reference.txt",
        &[
            on("1:100", "DownstreamProtein="),
            on("1:200", "DownstreamProtein=KL"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            on("1:100", "DownstreamProtein=K"),
            on("1:200", "DownstreamProtein="),
        ],
    );
    let table = format!(
        "planted_other_class\tDownstream\tDownstreamProtein\tmet\t-\n{shifted}\tDownstream\tDownstreamProtein\tmet\t-\n"
    );
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &Run::new()
            .residuals(&residuals)
            .table(&table)
            .plugins(&["Downstream"]),
    );
    let calls = planted.calls();
    let handed: Vec<(&str, &str, &str)> = calls
        .iter()
        .map(|call| {
            (
                call.klass.as_str(),
                call.expected.as_str(),
                call.actual.as_str(),
            )
        })
        .collect();
    assert_eq!(
        handed,
        [
            ("planted_other_class", ABSENT, "K"),
            (shifted, "", "K"),
            ("planted_other_class", "KL", ABSENT),
            (shifted, "KL", ""),
        ]
    );
    let item = feature(&report, "Downstream");
    assert_eq!(
        item["masked"],
        json!({"planted_other_class": 0, shifted: 2})
    );
    assert_eq!(counts(item), [1, 1, 0, 0, 0, 0]);
    assert!(feature_residuals(&residuals, "Downstream").is_empty());
}

/// A row of the record `uploaded` that writes no cDNA position.
fn named(uploaded: &str, location: &str, allele: &str, feature: &str, extra: &str) -> String {
    row_with(location, allele, feature, extra, "-", uploaded)
}

/// The fields of the entry a record's rows end in under an identifier: its input line.
fn line_fields(line: &Entry) -> BTreeMap<&str, &str> {
    line.fields.iter().collect()
}

fn input_line<'a>(
    name: &'a str,
    position: &'a str,
    alternates: &'a str,
) -> BTreeMap<&'a str, &'a str> {
    BTreeMap::from([
        ("CHROM", "1"),
        ("POS", position),
        ("Uploaded_variation", name),
        ("REF", "C"),
        ("ALT", alternates),
    ])
}

#[test]
fn row_at_a_two_base_location_under_an_identifier_is_handed_its_input_line() {
    let scratch = Scratch::new("two-base-location");
    let planted = Planted {
        registry: PlantedRegistry::planted(),
        rules: PlantedRules::none().class("planted_inserted_length", &[], |c| {
            c.record_rows.is_some_and(|rows| {
                rows.iter().any(|other| {
                    (other.fields.get("REF"), other.fields.get("ALT")) == (Some("C"), Some("CT"))
                })
            })
        }),
    };
    let table = "planted_inserted_length\tUTRAnnotator\t5UTR_consequence\n";
    let on = |term: &str| {
        row_with(
            "1:100-101",
            "T",
            "ENST1",
            &format!("5UTR_consequence={term}"),
            "64-65",
            "ins1",
        )
    };
    let reference = scratch.output("reference.txt", &[on("lost")]);
    let engine = scratch.output("engine.txt", &[on("gained")]);
    let run_input = scratch.vcf(
        "input.vcf",
        &[
            ["1", "100", "ins1", "C", "CT"],
            ["1", "100", "other", "C", "G"],
        ],
    );
    let run = Run::new().table(table).plugins(&["UTRAnnotator"]);
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &run);
    let without = feature(&report, "UTRAnnotator");
    assert_eq!(counts(without), [1, 1, 0, 1, 1, 0]);
    assert_eq!(without["masked"], json!({"planted_inserted_length": 0}));
    let calls = planted.calls();
    assert_eq!(calls.len(), 1);
    assert!(calls[0].record_rows.is_none());
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &run.run_input(&run_input),
    );
    let item = feature(&report, "UTRAnnotator");
    assert_eq!(counts(item), [1, 1, 0, 0, 0, 0]);
    assert_eq!(item["masked"], json!({"planted_inserted_length": 1}));
    let calls = planted.calls();
    assert_eq!(calls.len(), 2);
    let [row_entry, line_entry] = calls[1].record_rows.as_deref().unwrap() else {
        panic!("the record's rows are not the row and its input line");
    };
    assert_eq!(
        (
            row_entry.location.as_deref(),
            row_entry.allele.as_str(),
            row_entry.fields.get("cDNA_position"),
        ),
        (Some("1:100-101"), "T", Some("64-65"))
    );
    assert_eq!(line_fields(line_entry), input_line("ins1", "100", "CT"));
}

const NAMED: [(&str, &str); 3] = [("1:600", "rs6"), ("1:800", "rs8"), ("1:900", "rs9")];
const NAMED_RECORDS: [[&str; 5]; 5] = [
    ["1", "600", "rs6", "C", "A,G"],
    ["1", "600", "rs66", "C", "A"],
    ["1", "800", "rs8", "C", "A"],
    ["1", "800", "rs8", "C", "G"],
    ["1", "50", ".", "C", "T"],
];

/// A reference whose rows under three identifiers write a REVEL score under two alleles
/// each, an engine that writes it under the first allele alone, and a tuple file of the
/// engine's values.
fn named_rows(scratch: &Scratch) -> (String, String, String) {
    let other = named("rs66", "1:600", "A", "ENST2", "REVEL=0.7");
    let rows = |later: &str| -> Vec<String> {
        let mut rows: Vec<String> = NAMED
            .iter()
            .flat_map(|(location, name)| {
                [("A", "REVEL=0.7"), ("G", later)]
                    .map(|(allele, extra)| named(name, location, allele, "ENST1", extra))
            })
            .collect();
        rows.push(other.clone());
        rows
    };
    let reference = scratch.output("reference.txt", &rows("REVEL=0.7"));
    let engine = scratch.output("engine.txt", &rows("IMPACT=LOW"));
    let mut paired: Vec<[&str; 5]> = NAMED
        .iter()
        .map(|(location, _)| [*location, "A", "ENST1", "REVEL", "0.7"])
        .collect();
    paired.push(["1:600", "A", "ENST2", "REVEL", "0.7"]);
    (reference, engine, scratch.tuples("engine.tsv", &paired))
}

/// The Allele and the `ALT` field of each entry of a record's rows.
fn alleles_and_alternates(rows: &[Entry]) -> Vec<(&str, Option<&str>)> {
    rows.iter()
        .map(|other| (other.allele.as_str(), other.fields.get("ALT")))
        .collect()
}

/// The rows of a record under an identifier are those of its Location and name, with the
/// input line that carries the identifier after them; an identifier no line carries names
/// no record. Where two lines carry it, the row is handed the line that writes it.
fn input_names_the_record_of_rows_under_an_identifier(zipped: bool) {
    let scratch = Scratch::new(if zipped { "named-gzip" } else { "named-plain" });
    let (reference, engine, listed) = named_rows(&scratch);
    let mut run_input = scratch.vcf("input.vcf", &NAMED_RECORDS);
    if zipped {
        run_input = scratch.file("input.vcf.gz", gzip(&fs::read(&run_input).unwrap()));
    }
    let engines = [("engine", engine.as_str()), ("listed", listed.as_str())];
    let left: Vec<String> = NAMED
        .iter()
        .map(|(location, _)| format!("reference\t{location}\tG\tENST1\tREVEL\t0.7"))
        .collect();
    let planted = Planted::new();
    let residuals = scratch.path("without");
    let without = planted.score(
        &scratch,
        &reference,
        &engines,
        &Run::new().residuals(&residuals),
    );
    assert!(without["run_input"].is_null());
    for (ident, _) in engines {
        let item = entry(&without, "REVEL", ident, "feature");
        assert_eq!(counts(item), [7, 4, 4, 7, 4, 4], "{ident}");
        assert_eq!(
            residual_lines(&residuals, "REVEL", ident, "feature"),
            left,
            "{ident}"
        );
    }
    let calls = planted.calls();
    assert!(!calls.is_empty());
    assert!(calls.iter().all(|call| call.record_rows.is_none()));
    let planted = Planted::new();
    let residuals = scratch.path("with");
    let report = planted.score(
        &scratch,
        &reference,
        &engines,
        &Run::new().residuals(&residuals).run_input(&run_input),
    );
    assert_eq!(
        report["run_input"],
        json!({
            "path": run_input,
            "sha256": sha256(&fs::read(&run_input).unwrap()),
            "records": 5,
        })
    );
    for (ident, _) in engines {
        let item = entry(&report, "REVEL", ident, "feature");
        assert_eq!(counts(item), [7, 4, 4, 6, 4, 4], "{ident}");
        assert_eq!(
            item["masked"],
            json!({"planted_sibling_value": 1}),
            "{ident}"
        );
        assert_eq!(
            residual_lines(&residuals, "REVEL", ident, "feature"),
            left[1..],
            "{ident}"
        );
    }
    let calls = planted.calls();
    assert_eq!(calls.len(), 6);
    for call in &calls {
        let rows = call.record_rows.as_deref();
        match call.e.location.as_deref().unwrap() {
            "1:600" => assert_eq!(
                alleles_and_alternates(rows.unwrap()),
                [("A", None), ("G", None), ("", Some("A,G"))]
            ),
            "1:800" => {
                let rows = rows.unwrap();
                assert_eq!(alleles_and_alternates(rows), [("G", None), ("", Some("G"))]);
                assert_eq!(line_fields(&rows[1]), input_line("rs8", "800", "G"));
            }
            "1:900" => assert!(rows.is_none()),
            other => panic!("a consult at {other}"),
        }
    }
}

#[test]
fn input_names_the_record_of_rows_under_an_identifier_that_one_line_carries() {
    input_names_the_record_of_rows_under_an_identifier(false);
}

#[test]
fn gzip_input_names_the_record_of_rows_under_an_identifier_that_one_line_carries() {
    input_names_the_record_of_rows_under_an_identifier(true);
}

#[test]
fn allele_order_of_a_record_under_an_identifier_comes_from_its_input_line() {
    let scratch = Scratch::new("allele-order");
    let planted = Planted {
        registry: PlantedRegistry::planted(),
        rules: PlantedRules::none().class("planted_later_alternate", &[], |c| {
            c.record_rows.is_some_and(|rows| {
                rows.iter().any(|other| {
                    other
                        .fields
                        .get("ALT")
                        .is_some_and(|alternates| alternates.split(',').next() != Some(&c.e.allele))
                })
            })
        }),
    };
    let table = "planted_later_alternate\tREVEL\tREVEL\n";
    let rows = |score: &str| -> Vec<String> {
        ["G", "T"]
            .iter()
            .map(|allele| named("rs1", "1:100", allele, "ENST1", &format!("REVEL={score}")))
            .collect()
    };
    let reference = scratch.output("reference.txt", &rows("0.7"));
    let engine = scratch.output("engine.txt", &rows("0.9"));
    let run_input = scratch.vcf("input.vcf", &[["1", "100", "rs1", "A", "G,T"]]);
    let run = Run::new().table(table);
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &run);
    let without = feature(&report, "REVEL");
    assert_eq!(counts(without), [2, 2, 0, 2, 2, 0]);
    assert_eq!(without["masked"], json!({"planted_later_alternate": 0}));
    let residuals = scratch.path("residuals");
    let report = planted.score(
        &scratch,
        &reference,
        &[("engine", &engine)],
        &run.residuals(&residuals).run_input(&run_input),
    );
    let item = feature(&report, "REVEL");
    assert_eq!(counts(item), [2, 2, 0, 1, 1, 0]);
    assert_eq!(item["masked"], json!({"planted_later_alternate": 1}));
    assert_eq!(
        feature_residuals(&residuals, "REVEL"),
        [
            "reference\t1:100\tG\tENST1\tREVEL\t0.7",
            "engine\t1:100\tG\tENST1\tREVEL\t0.9",
        ]
    );
}

/// A row is found in the input by the bytes of its identifier: under an identifier that is
/// no UTF-8 the record's rows end in the line that carries those bytes, shown with the
/// replacement character, and a line elsewhere whose identifier is that shown text is no
/// carrier.
#[test]
fn identifier_that_is_not_utf8_names_the_lines_that_carry_its_bytes() {
    let scratch = Scratch::new("identifier-bytes");
    let (name, shown): (&[u8], &str) = (b"rs\xfe", "rs\u{fffd}");
    let written = |file: &str, rows: &[String]| -> String {
        let text = format!("{HEADER}{}", rows.concat());
        let parts: Vec<&[u8]> = text.split("NAME").map(str::as_bytes).collect();
        scratch.file(file, parts.join(name))
    };
    let reference = written(
        "reference.txt",
        &[
            named("NAME", "1:600", "A", "ENST1", "REVEL=0.7"),
            named("NAME", "1:600", "G", "ENST1", "REVEL=0.7"),
        ],
    );
    let engine = written(
        "engine.txt",
        &[
            named("NAME", "1:600", "A", "ENST1", "REVEL=0.7"),
            named("NAME", "1:600", "G", "ENST1", "IMPACT=LOW"),
        ],
    );
    // A record elsewhere whose identifier is the text the bytes of `name` decode to.
    let other = format!("2\t500\t{shown}\tC\tA,G");
    let cases: [(&str, &[&str], i64, Option<&str>); 3] = [
        ("one", &["A,G"], 1, Some("A,G")),
        ("two", &["A", "G"], 0, Some("G")),
        ("none", &[], 0, None),
    ];
    for (label, carriers, aside, handed) in cases {
        let mut content = VCF_HEADER.as_bytes().to_vec();
        for alternates in carriers {
            content.extend_from_slice(b"1\t600\t");
            content.extend_from_slice(name);
            content.extend_from_slice(format!("\tC\t{alternates}\t.\t.\t.\n").as_bytes());
        }
        content.extend_from_slice(format!("{other}\t.\t.\t.\n").as_bytes());
        let run_input = scratch.file(&format!("{label}.vcf"), content);
        let planted = Planted::new();
        let report = planted.score(
            &scratch,
            &reference,
            &[("engine", &engine)],
            &Run::new().run_input(&run_input),
        );
        assert_eq!(
            feature(&report, "REVEL")["masked"],
            json!({"planted_sibling_value": aside}),
            "{label}"
        );
        let calls = planted.calls();
        assert_eq!(calls.len(), 1, "{label}");
        let line = calls[0]
            .record_rows
            .as_deref()
            .map(|rows| line_fields(rows.last().unwrap()));
        assert_eq!(
            line,
            handed.map(|alternates| input_line(shown, "600", alternates)),
            "{label}"
        );
    }
}
