// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The passes as the `score` program runs them and the settling of the files they leave,
//! under a planted registry and planted classes: the rows a run keeps, the bound on the
//! keys no source claims, the refusal of a file the settling step cannot read forward or
//! that lacks the rows of a key, and a scratch that cannot be removed.

mod support;

use serde_json::json;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use support::{
    counts, feature, gzip, names, raw_counts, reference_rows, refusal, row, sha256, Planted,
    PlantedRegistry, PlantedRules, Run, Scratch, CLASS_TABLE, ERROR, HEADER, PLANTED,
};
use vep_plugin_pairs::resolve::Registry as MetFields;
use vep_plugin_pairs::score::api::{score_engine, Classes, Sources, Spread, Tally, RESOLVED_LIMIT};
use vep_plugin_pairs::score::{pass, reads};

/// The passes over `reference` and `engines` as a run of `score` starts them, their files
/// written into `pairs` of the scratch; the message of the refusal where they are refused.
fn passes(
    scratch: &Scratch,
    sources: &Sources,
    reference: &str,
    engines: &[&str],
    options: &pass::Options,
) -> Result<pass::Passed, String> {
    let direct = pass::Direct { sources };
    let met = MetFields::new(&direct);
    let engines: Vec<PathBuf> = engines.iter().map(PathBuf::from).collect();
    pass::run(
        Path::new(reference),
        &engines,
        &scratch.0.join("pairs"),
        options,
        &met,
        &|_| {},
    )
    .map_err(|refusal| refusal.0)
}

/// The tallies of `engine` against `reference` under the class table `table`: the passes
/// run, `tamper` is handed the directory of their files, and the settling step reads what
/// it left there. The message of the refusal where the settling step refuses.
fn settled(
    scratch: &Scratch,
    planted: &Planted,
    table: &str,
    reference: &str,
    engine: &str,
    tamper: impl FnOnce(&Path),
) -> Result<[Tally; 2], String> {
    let plugins = planted.registry.plugins();
    let sources = Sources::new(&plugins, &[], None, &planted.registry).unwrap();
    let table = scratch.file("classes.tsv", table);
    let classes = Classes::new(
        Some(Path::new(&table)),
        None,
        None,
        &[],
        &sources,
        false,
        &planted.rules,
    )
    .unwrap();
    let mut options = pass::Options::new(!classes.declared().is_empty());
    options.threads = 2;
    let passed = passes(scratch, &sources, reference, &[engine], &options)
        .unwrap_or_else(|message| panic!("the passes were refused: {message}"));
    let out = scratch.0.join("pairs");
    tamper(&out);
    let against = pass::extracted(
        "reference",
        reference,
        &passed.inputs[0],
        &out,
        1,
        "reference",
    );
    let written = pass::extracted("engine", engine, &passed.inputs[1], &out, 1, "engine");
    let context = reads::Context::read(&reads::Named::default(), &classes).unwrap();
    score_engine(
        &against,
        &written,
        &pass::counted(&passed.engines[0]),
        &scratch.0.join("input1.verdicts"),
        &sources,
        &classes,
        &planted.rules,
        &context.reads,
        &context.plugin_reads,
        None,
        None,
        Spread::over(options.threads),
    )
    .map_err(|error| error.0)
}

/// Rewrites the file with the lines `kept` gives for the ones it holds.
fn rewrite(path: &Path, kept: impl FnOnce(Vec<&str>) -> Vec<&str>) {
    let text = fs::read_to_string(path).unwrap();
    let lines = kept(text.split_inclusive('\n').collect());
    fs::write(path, lines.concat()).unwrap();
}

#[test]
fn no_rows_are_kept_when_no_class_is_declared() {
    let scratch = Scratch::new("no-rows");
    // The files of the first engine's remainder while the second engine is read: the
    // registry lists them when it is asked for a field the second engine alone writes.
    let listed: Arc<Mutex<Vec<Vec<String>>>> = Arc::default();
    let (seen, work) = (Arc::clone(&listed), scratch.0.join("work"));
    let planted = Planted {
        registry: PlantedRegistry::planted().resolving(move |field, plugins, _| {
            if field == "Late" {
                let files = names(&work)
                    .iter()
                    .flat_map(|run| names(work.join(run).join("pairs")))
                    .filter(|name| name.starts_with("input1."))
                    .collect();
                seen.lock().unwrap().push(files);
            }
            let writers: Vec<&str> = PLANTED
                .iter()
                .filter(|(name, fields)| {
                    fields.contains(&field) && plugins.iter().any(|loaded| loaded == name)
                })
                .map(|(name, _)| *name)
                .collect();
            Ok(match writers.as_slice() {
                [one] => Some(one.to_string()),
                _ => None,
            })
        }),
        rules: PlantedRules::planted(),
    };
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=1"),
            row("1:200", "T", "ENST1", "PLANT_score=2"),
        ],
    );
    let first = scratch.output(
        "first.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=1"),
            row("1:200", "T", "ENST1", "PLANT_score=3"),
        ],
    );
    let second = scratch.output(
        "second.txt",
        &[
            row("1:100", "T", "ENST1", "PLANT_score=1;Late=1"),
            row("1:200", "T", "ENST1", "PLANT_score=2"),
        ],
    );
    let engines = [("first", first.as_str()), ("second", second.as_str())];
    let tuple_files = [
        "input1.allele.engine",
        "input1.allele.reference",
        "input1.feature.engine",
        "input1.feature.reference",
    ];
    planted.score(&scratch, &reference, &engines, &Run::new().no_table());
    let without_a_class = std::mem::take(&mut *listed.lock().unwrap());
    assert!(!without_a_class.is_empty());
    for files in without_a_class {
        assert_eq!(files, tuple_files);
    }
    planted.score(&scratch, &reference, &engines, &Run::new());
    let with_classes = std::mem::take(&mut *listed.lock().unwrap());
    assert!(!with_classes.is_empty());
    for files in with_classes {
        let mut every = tuple_files.to_vec();
        every.extend(["input1.rows.engine", "input1.rows.reference"]);
        assert_eq!(files, every);
    }
}

#[test]
fn distinct_unclaimed_keys_are_bounded() {
    let scratch = Scratch::new("distinct-keys");
    let planted = Planted::new();
    let plugins = planted.registry.plugins();
    let sources = Sources::new(&plugins, &[], None, &planted.registry).unwrap();
    let mut options = pass::Options::new(false);
    options.threads = 2;
    options.distinct_keys = 2;
    let reference = scratch.output(
        "reference.txt",
        &[row("1:100", "T", "ENST1", "One=1;Two=2;PLANT_score=1")],
    );
    let passed = passes(&scratch, &sources, &reference, &[&reference], &options)
        .unwrap_or_else(|message| panic!("the passes were refused: {message}"));
    let read = pass::extracted(
        "reference",
        &reference,
        &passed.inputs[0],
        &scratch.0.join("pairs"),
        1,
        "reference",
    );
    assert_eq!(
        read.unclaimed,
        [("One".to_string(), 1), ("Two".to_string(), 1)]
    );
    let wide = scratch.output(
        "wide.txt",
        &[row(
            "1:100",
            "T",
            "ENST1",
            "One=1;Two=2;Three=3;PLANT_score=1",
        )],
    );
    let listed = scratch.tuples(
        "wide.tsv",
        &["One", "Two", "Three"].map(|name| ["1:100", "T", "ENST1", name, "1"]),
    );
    for (path, name) in [(&wide, "wide.txt"), (&listed, "wide.tsv")] {
        let refused = passes(&scratch, &sources, &reference, &[path], &options)
            .err()
            .unwrap_or_else(|| panic!("{name} was not refused"));
        assert!(
            refused.ends_with(&format!(
                "{name}: more than 2 distinct keys that no source claims"
            )),
            "{refused}"
        );
        assert!(scratch.names("pairs/shards").is_empty(), "{name}");
    }
}

#[test]
fn distinct_unclaimed_keys_of_a_whole_run_are_bounded_at_the_runs_own_bound() {
    let scratch = Scratch::new("distinct-keys-run");
    let planted = Planted::registry(&[("Planted", &["PLANT_score"])]);
    let run = Run::new().no_table();
    let keys = |count: usize| -> Vec<String> { (0..count).map(|at| format!("K{at}")).collect() };
    let wide_row = |count: usize| {
        let mut extra: Vec<String> = keys(count).iter().map(|name| format!("{name}=1")).collect();
        extra.push("PLANT_score=1".to_string());
        row("1:100", "T", "ENST1", &extra.join(";"))
    };
    let reference = scratch.output("reference.txt", &[wide_row(RESOLVED_LIMIT)]);
    let report = planted.score(&scratch, &reference, &[("engine", &reference)], &run);
    let mut at_the_bound = keys(RESOLVED_LIMIT);
    at_the_bound.sort();
    for stated in [&report["reference"], &report["engines"][0]] {
        let listed: Vec<(&str, Option<u64>)> = stated["unclaimed_keys"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(name, times)| (name.as_str(), times.as_u64()))
            .collect();
        assert_eq!(listed.len(), RESOLVED_LIMIT);
        assert!(listed
            .iter()
            .zip(&at_the_bound)
            .all(|(listed, name)| *listed == (name.as_str(), Some(1))));
    }
    assert_eq!(raw_counts(feature(&report, "Planted")), [1, 1, 1]);
    let wide = scratch.output("wide.txt", &[wide_row(RESOLVED_LIMIT + 1)]);
    let one_more = keys(RESOLVED_LIMIT + 1);
    let lines: Vec<[&str; 5]> = one_more
        .iter()
        .map(|name| ["1:100", "T", "ENST1", name.as_str(), "1"])
        .collect();
    let listed = scratch.tuples("wide.tsv", &lines);
    for path in [&wide, &listed] {
        let refused = refusal(planted.try_score(&scratch, &reference, &[("engine", path)], &run));
        assert_eq!(
            refused,
            format!("{path}: more than {RESOLVED_LIMIT} distinct keys that no source claims")
        );
        assert!(scratch.names("work").is_empty(), "{path}");
    }
}

#[test]
fn passes_on_fewer_threads_than_inputs_write_the_same_files_and_digests() {
    let scratch = Scratch::new("few-threads");
    let planted = Planted::new();
    let plugins = planted.registry.plugins();
    let sources = Sources::new(&plugins, &[], None, &planted.registry).unwrap();
    // Every `changed`-th row differs from the reference's.
    let rows = |changed: usize| -> Vec<String> {
        (1..=300)
            .map(|at| {
                let score = if at % changed == 0 { at + 1 } else { at };
                row(
                    &format!("{}:{}", at % 5 + 1, at * 10),
                    "T",
                    &format!("ENST{}", at % 3),
                    &format!("PLANT_score={score};LOEUF={}", at % 11),
                )
            })
            .collect()
    };
    let reference = scratch.output("reference.txt", &rows(301));
    let plain = scratch.output("plain.txt", &rows(7));
    let packed = scratch.file(
        "packed.txt.gz",
        gzip(format!("{HEADER}{}", rows(5).concat()).as_bytes()),
    );
    let listed = scratch.tuples(
        "listed.tsv",
        &[
            ["1:50", "T", "ENST2", "PLANT_score", "5"],
            ["2:10", "T", "ENST1", "PLANT_score", "2"],
            ["3:20", "T", "ENST2", "LOEUF", "2"],
        ],
    );
    let engines = [plain.as_str(), packed.as_str(), listed.as_str()];
    let stored: Vec<String> = [reference.as_str()]
        .iter()
        .chain(&engines)
        .map(|path| sha256(&fs::read(path).unwrap()))
        .collect();
    let mut runs = Vec::new();
    for threads in [1, 2, 6] {
        // Blocks and shards small enough that every pass hands out more than a thread takes.
        let mut options = pass::Options::new(true);
        options.threads = threads;
        options.shards = Some(5);
        options.block_bytes = 256;
        let passed = passes(&scratch, &sources, &reference, &engines, &options)
            .unwrap_or_else(|message| panic!("the passes were refused: {message}"));
        let digests: Vec<String> = passed
            .inputs
            .iter()
            .map(|input| input.sha256.clone())
            .collect();
        assert_eq!(digests, stored, "{threads} threads");
        let lines: Vec<u64> = passed
            .inputs
            .iter()
            .map(|input| input.extracted.lines)
            .collect();
        assert_eq!(lines, [300, 300, 300, 3], "{threads} threads");
        let files: Vec<(String, Vec<u8>)> = scratch
            .names("pairs")
            .into_iter()
            .map(|name| {
                let bytes = fs::read(scratch.0.join("pairs").join(&name)).unwrap();
                (name, bytes)
            })
            .collect();
        assert_eq!(files.len(), 18, "{threads} threads");
        fs::remove_dir_all(scratch.0.join("pairs")).unwrap();
        runs.push(files);
    }
    assert!(runs[0].iter().any(|(_, bytes)| !bytes.is_empty()));
    assert_eq!(runs[0], runs[1]);
    assert_eq!(runs[0], runs[2]);
}

#[test]
fn remainder_file_that_is_not_ascending_ends_the_settling_in_the_error_line() {
    let scratch = Scratch::new("not-ascending");
    let planted = Planted::new();
    let rows = |values: [&str; 2]| {
        [
            row("1:100", "T", "ENST1", &format!("PLANT_score={}", values[0])),
            row("1:200", "T", "ENST1", &format!("PLANT_score={}", values[1])),
        ]
    };
    let reference = scratch.output("reference.txt", &rows(["1", "2"]));
    let engine = scratch.output("engine.txt", &rows(["3", "4"]));
    let in_order = settled(&scratch, &planted, "", &reference, &engine, |_| {});
    assert!(in_order.is_ok(), "{in_order:?}");
    for part in [
        "feature.reference",
        "feature.engine",
        "allele.reference",
        "allele.engine",
    ] {
        let name = format!("input1.{part}");
        let refused = settled(&scratch, &planted, "", &reference, &engine, |out| {
            rewrite(&out.join(&name), |mut lines| {
                assert_eq!(lines.len(), 2, "{name}");
                lines.reverse();
                lines
            })
        });
        assert_eq!(
            refused.err(),
            Some(format!(
                "{name}: line 2 sorts below line 1: the file is not in ascending order"
            ))
        );
    }
}

#[test]
fn value_with_no_row_under_its_key_among_the_sorted_rows_ends_the_settling_in_the_error_line() {
    let scratch = Scratch::new("no-row");
    let planted = Planted::new();
    let reference = scratch.output(
        "reference.txt",
        &[
            row("1:100", "T", "ENST1", "GWAS_pmid=17"),
            row("1:200", "T", "ENST1", "GWAS_pmid=17"),
        ],
    );
    let engine = scratch.output(
        "engine.txt",
        &[
            row("1:100", "T", "ENST1", "GWAS_pmid=17"),
            row("1:200", "T", "ENST1", "GWAS_pmid=18"),
        ],
    );
    let report = planted.score(&scratch, &reference, &[("engine", &engine)], &Run::new());
    let item = feature(&report, "GWAS");
    assert_eq!(counts(item), [2, 2, 1, 1, 1, 1]);
    assert_eq!(item["masked"], json!({"planted_reference_value": 1}));
    for (side, written) in [("reference", &reference), ("engine", &engine)] {
        let refused = settled(
            &scratch,
            &planted,
            CLASS_TABLE,
            &reference,
            &engine,
            |out| {
                rewrite(&out.join(format!("input1.rows.{side}")), |lines| {
                    lines
                        .into_iter()
                        .filter(|line| !line.starts_with("1:200\t"))
                        .collect()
                })
            },
        );
        assert_eq!(
            refused.err(),
            Some(format!(
                "{written}: no row under Location 1:200, Allele T, Feature ENST1 was read from \
                 the sorted rows, where the input writes a value of GWAS_pmid"
            )),
            "{side}"
        );
    }
}

#[test]
fn scratch_that_cannot_be_removed_is_named_and_leaves_no_residual_file() {
    // SAFETY: `geteuid` reads the process's user and cannot fail.
    if unsafe { libc::geteuid() } == 0 {
        // The superuser removes an entry of any directory.
        return;
    }
    let scratch = Scratch::new("kept-scratch");
    let work = scratch.0.join("work");
    // The work directory loses its write permission once the run has made its scratch
    // there: the registry takes it away when it is asked for a field.
    let locked = work.clone();
    let planted = Planted {
        registry: PlantedRegistry::planted().resolving(move |field, _, _| {
            fs::set_permissions(&locked, fs::Permissions::from_mode(0o555)).unwrap();
            Ok(field.starts_with("PLANT_").then(|| "Planted".to_string()))
        }),
        rules: PlantedRules::planted(),
    };
    let reference = scratch.output("reference.txt", &reference_rows());
    let no_class = scratch.file("no-class.tsv", "");
    let (code, lines, error) = planted.run_main(&[
        "--reference",
        &reference,
        "--engine",
        &format!("engine={reference}"),
        "--plugins",
        "Planted",
        "--class-table",
        &no_class,
        "--work-dir",
        &scratch.path("work"),
        "--report",
        &scratch.path("report.json"),
        "--residuals",
        &scratch.path("residuals"),
    ]);
    let kept = names(&work);
    fs::set_permissions(&work, fs::Permissions::from_mode(0o755)).unwrap();
    assert_eq!((code, lines.len()), (2, 0));
    assert_eq!(kept.len(), 1, "{kept:?}");
    assert_eq!(
        error,
        format!(
            "{ERROR}[Errno 13] Permission denied: '{}'\n",
            work.join(&kept[0]).display()
        )
    );
    assert!(!scratch.exists("report.json"));
    assert!(scratch.exists("residuals") && scratch.names("residuals").is_empty());
}
