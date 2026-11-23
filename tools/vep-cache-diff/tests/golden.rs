// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs the binary over every row of `tests/fixtures/cases.tsv` and compares its exit code,
//! stdout, stderr and report byte for byte with the row and the goldens under
//! `tests/fixtures/expected/`, under every bound of `common::Bounds::ALL`; then over a
//! generated pair of caches, whose counts are known from how it is made.

mod common;

use common::{fixtures, in_json, run, text, Bounds, Scratch};
use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

struct Case {
    name: String,
    exit: i32,
    stderr_golden: bool,
    report_golden: bool,
    baseline: String,
    builder: String,
}

fn cases() -> Vec<Case> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    table
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            assert_eq!(columns.len(), 6, "{line}");
            let flag = |text: &str, yes: &str, no: &str| match text {
                _ if text == yes => true,
                _ if text == no => false,
                _ => panic!("{line}"),
            };
            Case {
                name: columns[0].to_string(),
                exit: columns[1].parse().unwrap(),
                stderr_golden: flag(columns[2], "golden", "empty"),
                report_golden: flag(columns[3], "golden", "none"),
                baseline: columns[4].to_string(),
                builder: columns[5].to_string(),
            }
        })
        .collect()
}

fn expected(name: &str, suffix: &str) -> String {
    let path = fixtures().join("expected").join(format!("{name}.{suffix}"));
    fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

#[test]
fn every_case_writes_its_goldens_under_every_bound() {
    let scratch = Scratch::new("cases");
    let report = scratch.path("report.json");
    let fixtures = fixtures();
    let on_stdout = |golden: String| {
        golden
            .replace("{{FIXTURES}}", fixtures.to_str().unwrap())
            .replace("{{REPORT}}", report.to_str().unwrap())
    };
    let all = cases();
    assert_eq!(all.len(), 25);
    for case in &all {
        for bounds in Bounds::ALL {
            let _ = fs::remove_file(&report);
            let output = run(
                Path::new(&case.baseline),
                Path::new(&case.builder),
                &report,
                bounds,
            );
            let what = format!("{} under {bounds:?}", case.name);
            assert_eq!(
                output.status.code(),
                Some(case.exit),
                "{what}: {}",
                text(&output.stderr)
            );
            assert_eq!(
                text(&output.stdout),
                on_stdout(expected(&case.name, "stdout")),
                "{what}"
            );
            let stderr = if case.stderr_golden {
                on_stdout(expected(&case.name, "stderr"))
            } else {
                String::new()
            };
            assert_eq!(text(&output.stderr), stderr, "{what}");
            if case.report_golden {
                let golden =
                    expected(&case.name, "json").replace("{{FIXTURES}}", &in_json(&fixtures));
                assert_eq!(fs::read_to_string(&report).unwrap(), golden, "{what}");
            } else {
                assert!(!report.exists(), "{what}: a report was written");
            }
        }
    }
}

#[test]
fn every_golden_belongs_to_a_case() {
    let named: BTreeSet<String> = cases()
        .iter()
        .flat_map(|case| {
            let mut files = vec![format!("{}.stdout", case.name)];
            if case.stderr_golden {
                files.push(format!("{}.stderr", case.name));
            }
            if case.report_golden {
                files.push(format!("{}.json", case.name));
            }
            files
        })
        .collect();
    let present: BTreeSet<String> = fs::read_dir(fixtures().join("expected"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(present, named);
}

/// What the generated pair is made to count.
#[derive(Debug, Default, PartialEq)]
struct Counts {
    names: [usize; 2],
    shared: usize,
    matched: usize,
    differed: usize,
}

/// One record of a generated shard.
fn transcript(id: &str, start: usize) -> String {
    format!(
        r#"{{"stable_id": "{id}", "start": {start}, "biotype": "protein_coding", "vefc": {{"sorted_exons": [{{"start": {start}}}], "translateable_seq": "{}"}}}}"#,
        "acgt".repeat(20)
    )
}

/// Writes 60 shards a side, each of 40 transcripts: one in eight is held by one side
/// alone, one in eight differs in its start.
fn write_transcripts(scratch: &Scratch) -> Counts {
    let mut counts = Counts::default();
    for shard in 0..60 {
        let (mut baseline, mut builder) = (Vec::new(), Vec::new());
        for n in 0..40 {
            let id = format!("T{shard:02}_{n:02}");
            let start = shard * 1000 + n;
            match (shard + n) % 8 {
                0 => {
                    baseline.push(transcript(&id, start));
                    counts.names[0] += 1;
                }
                1 => {
                    builder.push(transcript(&id, start));
                    counts.names[1] += 1;
                }
                rest => {
                    let changed = usize::from(rest == 2);
                    baseline.push(transcript(&id, start));
                    builder.push(transcript(&id, start + changed));
                    counts.names[0] += 1;
                    counts.names[1] += 1;
                    counts.shared += 1;
                    counts.matched += 1 - changed;
                    counts.differed += changed;
                }
            }
        }
        let name = format!("transcripts/{}/{shard}-{}.json", shard % 7, shard + 1);
        scratch.write(&format!("a/{name}"), format!("[{}]", baseline.join(", ")));
        // Every third shard of the builder is gzipped.
        let text = format!("[{}]", builder.join(",\n"));
        if shard % 3 == 0 {
            scratch.write_gzipped(&format!("b/{name}.gz"), text);
        } else {
            scratch.write(&format!("b/{name}"), text);
        }
    }
    counts
}

/// Writes five chromosomes a side, each of 4,000 names: a share of them repeated, changed
/// or held by one side alone.
fn write_variations(scratch: &Scratch) -> Counts {
    let columns =
        r#"["chr", "variation_name", "failed", "somatic", "start", "end", "allele_string", "AF"]"#;
    for side in ["a", "b"] {
        let info = format!(r#"{{"variation_cols": {columns}}}"#);
        scratch.write(&format!("{side}/info.json"), info);
    }
    let mut counts = Counts::default();
    for chromosome in 1..=5 {
        let (mut baseline, mut builder) = (String::new(), String::new());
        for i in 0..4000 {
            // A name repeats across chromosomes and counts once in each.
            let name = format!("rs{}", i * 7919 % 100_003);
            let line = |start: usize, af: &str| {
                format!("{chromosome}\t{name}\t.\t.\t{start}\t.\tA/G\t{af}\n")
            };
            match (i + chromosome) % 10 {
                0 => {
                    baseline += &line(i, ".");
                    baseline += &line(i + 1, ".");
                    counts.names[0] += 1;
                }
                1 => {
                    builder += &line(i, ".");
                    counts.names[1] += 1;
                }
                // The last record of a repeated name is the one compared.
                2 => {
                    baseline += &line(i + 5, "G:0.1");
                    baseline += &line(i, "G:0.2");
                    builder.insert_str(0, &line(i + 9, "G:0.9"));
                    builder += &line(i, "G:0.2");
                    counts.names[0] += 1;
                    counts.names[1] += 1;
                    counts.shared += 1;
                    counts.matched += 1;
                }
                rest => {
                    let changed = usize::from(rest == 3);
                    baseline += &line(i, "G:0.5");
                    builder += &line(i + changed, "G:0.5");
                    counts.names[0] += 1;
                    counts.names[1] += 1;
                    counts.shared += 1;
                    counts.matched += 1 - changed;
                    counts.differed += changed;
                }
            }
        }
        let name = format!("variations/{chromosome}/all_vars.gz");
        scratch.write_gzipped(&format!("a/{name}"), baseline);
        scratch.write_gzipped(&format!("b/{name}"), builder);
    }
    counts
}

/// The block of the report under `name`, as its lines.
fn block<'a>(report: &'a str, name: &str) -> Vec<&'a str> {
    report
        .lines()
        .skip_while(|line| *line != format!("  \"{name}\": {{"))
        .skip(1)
        .take_while(|line| line.starts_with("    "))
        .collect()
}

fn opening_lines(counts: &Counts) -> Vec<String> {
    [
        ("total_baseline", counts.names[0]),
        ("total_builder", counts.names[1]),
        ("shared", counts.shared),
        ("only_baseline", counts.names[0] - counts.shared),
        ("only_builder", counts.names[1] - counts.shared),
        ("matched", counts.matched),
        ("differed", counts.differed),
    ]
    .iter()
    .map(|(name, count)| format!("    \"{name}\": {count},"))
    .collect()
}

#[test]
fn a_generated_pair_reports_its_counts_under_every_bound() {
    let scratch = Scratch::new("generated");
    let transcripts = write_transcripts(&scratch);
    let variations = write_variations(&scratch);
    let report = scratch.path("out/report.json");
    let mut first: Option<(Vec<u8>, String)> = None;
    let bounds = Bounds::ALL.into_iter().chain([
        Bounds {
            threads: Some("7"),
            buffer: Some("50000"),
        },
        Bounds {
            threads: Some("64"),
            buffer: Some("2000"),
        },
    ]);
    for bounds in bounds {
        // A record at a time is tried on the fixtures; here it would take thousands of
        // passes over every partition.
        if bounds.buffer == Some("1") {
            continue;
        }
        let _ = fs::remove_file(&report);
        let output = run(&scratch.path("a"), &scratch.path("b"), &report, bounds);
        assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
        assert!(output.stderr.is_empty(), "{}", text(&output.stderr));
        let written = fs::read_to_string(&report).unwrap();
        match &first {
            None => first = Some((output.stdout, written)),
            Some((stdout, report)) => {
                assert_eq!(&output.stdout, stdout, "{bounds:?}");
                assert_eq!(&written, report, "{bounds:?}");
            }
        }
    }
    let (stdout, written) = first.unwrap();
    assert_eq!(
        block(&written, "transcripts")[..7],
        opening_lines(&transcripts)
    );
    assert_eq!(
        block(&written, "variations")[..7],
        opening_lines(&variations)
    );
    assert_eq!(
        block(&written, "variations")[7..10],
        [
            "    \"field_mismatches\": {",
            &format!("      \"start\": {}", variations.differed),
            "    },",
        ]
    );
    let line = |counts: &Counts| {
        // The rate is rounded to six places, then written to four.
        let rate = format!("{:.6}", counts.matched as f64 / counts.shared as f64);
        let rate: f64 = rate.parse().unwrap();
        format!(
            "    {} shared, {} matched, {} differed, rate={rate:.4}\n",
            counts.shared, counts.matched, counts.differed
        )
    };
    let stdout = text(&stdout);
    assert!(stdout.contains(&line(&transcripts)), "{stdout}");
    assert!(stdout.contains(&line(&variations)), "{stdout}");
    // Fifty of the shard names with differences are listed, the first in sorted order.
    let listed: Vec<&str> = written
        .lines()
        .filter(|line| line.starts_with("        \"region\": "))
        .collect();
    assert_eq!(listed.len(), 50);
    assert_eq!(listed[0], "        \"region\": \"0/0-1.json\",");
    assert_eq!(listed[1], "        \"region\": \"0/14-15.json\",");
}

#[test]
fn records_written_out_leave_no_scratch_file_behind() {
    let scratch = Scratch::new("temporary");
    write_variations(&scratch);
    let temporary = scratch.path("tmp");
    fs::create_dir(&temporary).unwrap();
    let run_under = |buffer: Option<&'static str>, tmpdir: &Path| {
        let args: Vec<std::ffi::OsString> = [
            "--baseline-dir".into(),
            scratch.path("a").into_os_string(),
            "--builder-dir".into(),
            scratch.path("b").into_os_string(),
            "--output-report".into(),
            scratch.path("report.json").into_os_string(),
        ]
        .into();
        std::process::Command::new(common::BIN)
            .args(&args)
            .env("TMPDIR", tmpdir)
            .env(common::THREADS_ENV, "2")
            .env_remove(common::BUFFER_ENV)
            .envs(buffer.map(|bytes| (common::BUFFER_ENV, bytes)))
            .output()
            .unwrap()
    };
    let output = run_under(Some("5000"), &temporary);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    assert_eq!(fs::read_dir(&temporary).unwrap().count(), 0);
    // A `TMPDIR` that cannot be written ends the run where records are
    // written out, and nowhere else.
    let nowhere = scratch.path("tmp/nowhere");
    let output = run_under(None, &nowhere);
    assert_eq!(output.status.code(), Some(0), "{}", text(&output.stderr));
    let _ = fs::remove_file(scratch.path("report.json"));
    let output = run_under(Some("5000"), &nowhere);
    assert_eq!(output.status.code(), Some(1));
    let stderr = text(&output.stderr);
    assert!(
        stderr.starts_with("ERROR: [vep-cache-diff] temporary file of the variation comparison: ")
            && stderr.lines().count() == 1,
        "{stderr}"
    );
    assert!(text(&output.stdout).ends_with("  Comparing variations...\n"));
    assert!(!scratch.path("report.json").exists());
}
