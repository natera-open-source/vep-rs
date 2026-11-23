// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over every row of `tests/fixtures/cases.tsv` and compares its exit code,
//! stdout, stderr and report byte for byte with the row and the goldens under
//! `tests/fixtures/expected/`, on one thread, on two and on the host's default.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_vep-gtf-tags-compare");
const ROOT: &str = env!("CARGO_MANIFEST_DIR");
const THREADS_ENV: &str = "VEP_GTF_TAGS_COMPARE_THREADS";

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-gtf-tags-compare-golden-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn fixtures() -> PathBuf {
    Path::new(ROOT).join("tests/fixtures")
}

struct Case {
    name: String,
    exit: i32,
    stdout: String,
    stderr: String,
    args: Vec<String>,
}

fn cases() -> Vec<Case> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    table
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            assert_eq!(columns.len(), 5, "{line}");
            Case {
                name: columns[0].to_string(),
                exit: columns[1].parse().unwrap(),
                stdout: columns[2].to_string(),
                stderr: columns[3].to_string(),
                args: columns[4]
                    .split(' ')
                    .filter(|arg| !arg.is_empty())
                    .map(String::from)
                    .collect(),
            }
        })
        .collect()
}

fn run(args: &[String], threads: Option<&str>) -> Output {
    let mut command = Command::new(BIN);
    command.current_dir(fixtures()).args(args);
    match threads {
        Some(count) => command.env(THREADS_ENV, count),
        None => command.env_remove(THREADS_ENV),
    };
    command.output().unwrap()
}

/// What differs between one run of a case and its row and goldens.
fn differences(case: &Case, threads: Option<&str>, tag: &str) -> Vec<String> {
    let scratch = Scratch::new(&format!("{tag}-{}", case.name));
    let report = scratch.0.join("report.json");
    let args: Vec<String> = case
        .args
        .iter()
        .map(|arg| {
            arg.replace("{out}", report.to_str().unwrap())
                .replace("{scratch}", scratch.0.to_str().unwrap())
        })
        .collect();
    let output = run(&args, threads);
    let expected = fixtures().join("expected");
    let golden = |suffix: &str| fs::read(expected.join(format!("{}.{suffix}", case.name))).ok();
    let mut found = Vec::new();
    if output.status.code() != Some(case.exit) {
        found.push(format!(
            "exit {:?}, expected {}",
            output.status.code(),
            case.exit
        ));
    }
    match case.stdout.as_str() {
        "golden" if Some(&output.stdout) != golden("stdout").as_ref() => {
            found.push("stdout differs from its golden".to_string());
        }
        "empty" if !output.stdout.is_empty() => found.push("stdout is not empty".to_string()),
        "golden" | "empty" | "any" => {}
        other => found.push(format!("unknown stdout kind {other}")),
    }
    match case.stderr.as_str() {
        "golden" if Some(&output.stderr) != golden("stderr").as_ref() => {
            found.push("stderr differs from its golden".to_string());
        }
        "empty" if !output.stderr.is_empty() => found.push(format!(
            "stderr is not empty: {}",
            String::from_utf8_lossy(&output.stderr)
        )),
        "any" if output.stderr.is_empty() => found.push("stderr is empty".to_string()),
        "golden" | "empty" | "any" => {}
        other => found.push(format!("unknown stderr kind {other}")),
    }
    match (fs::read(&report).ok(), golden("json")) {
        (Some(written), Some(golden)) if written != golden => {
            found.push("the report differs from its golden".to_string());
        }
        (Some(_), None) => found.push("a report was written".to_string()),
        (None, Some(_)) => found.push("no report was written".to_string()),
        _ => {}
    }
    found
}

fn assert_every_case(threads: Option<&str>, tag: &str) {
    let cases = cases();
    assert!(cases.len() >= 70, "{} cases", cases.len());
    let failed: Vec<String> = cases
        .iter()
        .filter_map(|case| {
            let found = differences(case, threads, tag);
            (!found.is_empty()).then(|| format!("{}: {}", case.name, found.join("; ")))
        })
        .collect();
    assert!(failed.is_empty(), "{}", failed.join("\n"));
}

#[test]
fn every_case_matches_its_goldens_on_one_thread() {
    assert_every_case(Some("1"), "one");
}

#[test]
fn every_case_matches_its_goldens_on_two_threads() {
    assert_every_case(Some("2"), "two");
}

#[test]
fn every_case_matches_its_goldens_on_the_host_s_cores() {
    assert_every_case(None, "host");
}

#[test]
fn every_golden_belongs_to_a_case() {
    let names: BTreeSet<String> = cases().into_iter().map(|case| case.name).collect();
    for entry in fs::read_dir(fixtures().join("expected")).unwrap() {
        let file = entry.unwrap().file_name().to_string_lossy().into_owned();
        let (name, suffix) = file.rsplit_once('.').unwrap();
        assert!(
            names.contains(name) && ["stdout", "stderr", "json"].contains(&suffix),
            "{file} belongs to no case"
        );
    }
}

#[test]
fn a_thread_bound_that_is_no_whole_number_is_refused_before_any_input_is_read() {
    let args: Vec<String> = ["--gtf", "agree/genes.gtf", "--cache", "agree/cache"]
        .map(String::from)
        .to_vec();
    for bound in ["0", "x", "-1", "1.5"] {
        let output = run(&args, Some(bound));
        assert_eq!(output.status.code(), Some(1), "{bound}");
        assert!(output.stdout.is_empty(), "{bound}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!(
                "ERROR: [vep-gtf-tags-compare] {THREADS_ENV}={bound}: expected a whole number of threads, 1 or more\n"
            )
        );
    }
    let output = run(&args, Some(""));
    assert_eq!(output.status.code(), Some(0));
    let output = run(&args, Some(" 3 "));
    assert_eq!(output.status.code(), Some(0));
}

#[test]
fn a_failure_names_the_tool_and_the_file() {
    let args: Vec<String> = [
        "--gtf",
        "agree/genes.gtf",
        "--cache",
        "faults/trailing-comma",
    ]
    .map(String::from)
    .to_vec();
    for threads in ["1", "4"] {
        let output = run(&args, Some(threads));
        assert_eq!(output.status.code(), Some(1));
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.starts_with(
                "ERROR: [vep-gtf-tags-compare] faults/trailing-comma/transcripts/1/a.json: "
            ),
            "{stderr}"
        );
        assert_eq!(stderr.lines().count(), 1, "{stderr}");
    }
}

#[test]
fn a_stable_id_that_is_an_array_or_an_object_ends_the_run() {
    let scratch = Scratch::new("container-id");
    let shard = scratch.0.join("cache/transcripts/1");
    fs::create_dir_all(&shard).unwrap();
    fs::write(shard.join("a.json"), r#"[{"stable_id": ["T1"]}]"#).unwrap();
    let args: Vec<String> = [
        "--gtf",
        "agree/genes.gtf",
        "--cache",
        scratch.0.join("cache").to_str().unwrap(),
    ]
    .map(String::from)
    .to_vec();
    let output = run(&args, Some("1"));
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr)
        .ends_with("a.json: a stable_id that is an array or an object\n"));
}

/// A cache of many shards in one directory whose transcripts repeat across shards with
/// other codes: whatever the thread count, each stable id counts with the codes of its
/// first copy in file-name order.
#[test]
fn the_first_copy_of_a_stable_id_counts_on_every_thread_count() {
    const CODES: [&str; 4] = [
        "gencode_basic",
        "gencode_primary",
        "cds_start_NF",
        "cds_end_NF",
    ];
    let scratch = Scratch::new("first-copy");
    let dir = scratch.0.join("cache/transcripts/1");
    fs::create_dir_all(&dir).unwrap();
    // Shard `s` holds transcripts `20 * s .. 20 * s + 60`, so each id is stored in up to
    // three shards; the copy in shard `s` carries code `s % 4` alone.
    let shards = 60usize;
    for s in 0..shards {
        let body: Vec<String> = (20 * s..20 * s + 60)
            .map(|t| {
                format!(
                    r#"{{"exons":[{{"stable_id":"E{t}"}}],"stable_id":"T{t:05}.{s}","attributes":[{{"code":"{}"}},{{"code":"TSL"}}]}}"#,
                    CODES[s % 4]
                )
            })
            .collect();
        fs::write(
            dir.join(format!("{s:03}.json")),
            format!("[{}]", body.join(",")),
        )
        .unwrap();
    }
    let transcripts = 20 * (shards - 1) + 60;
    // The GTF tags every transcript with the code of the first shard that stores it.
    let tags = ["basic", "gencode_primary", "cds_start_NF", "cds_end_NF"];
    let gtf: String = (0..transcripts)
        .map(|t| {
            let first = (t / 20).saturating_sub(2);
            format!(
                "1\tsrc\ttranscript\t1\t9\t.\t+\t.\ttranscript_id \"T{t:05}\"; tag \"{}\";\n",
                tags[first % 4]
            )
        })
        .collect();
    let gtf_path = scratch.0.join("genes.gtf");
    fs::write(&gtf_path, gtf).unwrap();
    let report = scratch.0.join("report.json");
    let args: Vec<String> = [
        "--gtf",
        gtf_path.to_str().unwrap(),
        "--cache",
        scratch.0.join("cache").to_str().unwrap(),
        "--out",
        report.to_str().unwrap(),
    ]
    .map(String::from)
    .to_vec();
    let grouped = format!("{},{:03}", transcripts / 1000, transcripts % 1000);
    for threads in [Some("1"), Some("2"), Some("7"), Some("32"), None] {
        let _ = fs::remove_file(&report);
        let output = run(&args, threads);
        assert_eq!(
            output.status.code(),
            Some(0),
            "threads {threads:?}: {}",
            String::from_utf8_lossy(&output.stdout)
        );
        assert!(
            String::from_utf8_lossy(&output.stdout).starts_with(&format!(
            "transcripts: gtf {grouped} cache {grouped} shared {grouped} gtf_only 0 cache_only 0\n"
        ))
        );
        let written = fs::read_to_string(&report).unwrap();
        assert!(written.ends_with(" \"exact_agreement\": true\n}\n"));
        assert_eq!(written.matches("\"disagree\": 0,").count(), 4);
    }
}
