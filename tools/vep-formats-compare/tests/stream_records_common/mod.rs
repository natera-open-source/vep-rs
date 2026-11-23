// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! What the tests of the streamed VCF and JSON record groups share: a scratch directory,
//! the plan of a list of options, and a remainder read as the three texts a comparison
//! reads it from, the same at every worker count and block size.

#![allow(dead_code)]

use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};
use vep_formats_compare::stream::{self, Format, Plan};

static SCRATCH: AtomicU32 = AtomicU32::new(0);

pub struct Scratch(pub PathBuf);

impl Scratch {
    pub fn new() -> Scratch {
        let path = std::env::temp_dir().join(format!(
            "vep-formats-compare-stream-records-test-{}-{}",
            std::process::id(),
            SCRATCH.fetch_add(1, Ordering::Relaxed)
        ));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Scratch(path)
    }

    pub fn file(&self, name: &str, bytes: &[u8]) -> PathBuf {
        let path = self.0.join(name);
        fs::write(&path, bytes).unwrap();
        path
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A remainder as the three texts a comparison reads it from: each side's lines, each
/// ended by `\n`, and the counts.
#[derive(Debug, PartialEq)]
pub struct Remainder {
    pub expected: String,
    pub actual: String,
    pub counters: String,
}

pub fn worker_count() -> usize {
    std::thread::available_parallelism().map_or(1, usize::from)
}

/// The lines of a list file, an empty line being the empty value.
fn list_lines(path: &str) -> Vec<String> {
    let text = fs::read_to_string(path).unwrap();
    let mut lines: Vec<String> = text.split('\n').map(str::to_string).collect();
    if lines.last().is_some_and(String::is_empty) {
        lines.pop();
    }
    lines
}

/// The plan of a format, its input and the options `more`: `--force-records FILE` (one
/// ordinal a line), `--force-text FILE`, `--keep-values COLUMN FILE`, `--floor-fields
/// FILE` (a column and its fields a line, separated by tabs) and `--unordered`.
pub fn plan_of(format: Format, input: &Path, more: &[&str], threads: usize) -> Plan {
    let mut plan = Plan {
        format,
        force_locations: Vec::new(),
        force_text: Vec::new(),
        keep_values: None,
        unordered: false,
        threads,
        input: Some(input.to_path_buf()),
        force_records: Vec::new(),
        floor_fields: Vec::new(),
    };
    let mut more = more.iter();
    while let Some(option) = more.next() {
        let mut value = || *more.next().unwrap();
        match *option {
            "--force-records" => {
                let ordinals = list_lines(value());
                plan.force_records = ordinals.iter().map(|n| n.parse().unwrap()).collect();
            }
            "--force-text" => plan.force_text = list_lines(value()),
            "--keep-values" => {
                let column = value().to_string();
                plan.keep_values = Some((column, list_lines(value())));
            }
            "--floor-fields" => {
                let floor = |line: &String| {
                    let mut names = line.split('\t').map(str::to_string);
                    (names.next().unwrap(), names.collect())
                };
                plan.floor_fields = list_lines(value()).iter().map(floor).collect();
            }
            "--unordered" => plan.unordered = true,
            other => panic!("{other} is no option of a plan"),
        }
    }
    plan
}

fn json_string(text: &str) -> String {
    let mut out = String::from("\"");
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn runs_json(runs: &[(u64, u64)]) -> String {
    let runs: Vec<String> = runs.iter().map(|(n, b)| format!("[{n}, {b}]")).collect();
    format!("[{}]", runs.join(", "))
}

/// The counts of a remainder as one JSON text, for a format whose comparison numbers the
/// lines of an output.
pub fn counters_json(format: Format, r: &stream::Remainder) -> String {
    let whole = r.whole.as_deref().map_or("null".to_string(), json_string);
    let side = |lines: u64, handed_over: u64, bytes: u64| {
        format!(
            "{{\"lines\": {lines}, \"lines_handed_over\": {handed_over}, \"bytes_handed_over\": {bytes}}}"
        )
    };
    let block = |name: &str, counts: &vep_formats_compare::counter::Counter| {
        let items: Vec<String> = counts
            .iter()
            .map(|(name, n)| format!("  {}: {n}", json_string(name)))
            .collect();
        if items.is_empty() {
            format!("\"{name}\": {{}}")
        } else {
            format!("\"{name}\": {{\n{}\n }}", items.join(",\n"))
        }
    };
    let mut items = vec![
        format!("\"whole\": {whole}"),
        format!("\"kept_entries\": {}", r.kept_entries),
        block("kept_populated", &r.kept_populated),
    ];
    if format == Format::Json {
        items.push(block("kept_present", &r.kept_present));
        items.push(block("kept_under", &r.kept_under));
    } else {
        assert!(r.kept_present.is_empty() && r.kept_under.is_empty());
    }
    items.extend([
        format!("\"kept_before\": {}", runs_json(&r.kept_before)),
        format!("\"groups\": {}", r.groups),
        format!("\"groups_handed_over\": {}", r.groups_handed_over),
        format!(
            "\"expected\": {}",
            side(
                r.expected_lines,
                r.expected_handed_over,
                r.expected_bytes_handed_over
            )
        ),
        format!(
            "\"actual\": {}",
            side(
                r.actual_lines,
                r.actual_handed_over,
                r.actual_bytes_handed_over
            )
        ),
        format!(
            "\"line_runs\": {{\"expected\": {}, \"actual\": {}}}",
            runs_json(&r.expected_line_runs),
            runs_json(&r.actual_line_runs)
        ),
    ]);
    format!("{{\n {}\n}}\n", items.join(",\n "))
}

fn file_text(lines: &[String]) -> String {
    lines.iter().fold(String::new(), |mut text, line| {
        text.push_str(line);
        text.push('\n');
        text
    })
}

/// The remainder of one pair at `threads` workers, through blocks of `block` bytes
/// where one is given.
pub fn remainder_at(
    format: Format,
    expected: &Path,
    actual: &Path,
    input: &Path,
    more: &[&str],
    threads: usize,
    block: Option<usize>,
) -> stream::Remainder {
    let plan = plan_of(format, input, more, threads);
    let done = match block {
        Some(block) => stream::remainder_in_blocks(expected, actual, &plan, block),
        None => stream::remainder(expected, actual, &plan),
    };
    let done = done.unwrap_or_else(|error| panic!("{error}"));
    assert_eq!(done.restarts, 0);
    assert_eq!(done.candidate_error, None);
    done
}

pub fn texts_of(format: Format, run: &stream::Remainder) -> Remainder {
    Remainder {
        expected: file_text(&run.expected),
        actual: file_text(&run.actual),
        counters: counters_json(format, run),
    }
}

/// The remainder of one pair, the same at one worker, at two and at the default, and
/// through blocks shorter than a line, of a few lines and of the default size.
pub fn remainder_of(
    format: Format,
    expected: &Path,
    actual: &Path,
    input: &Path,
    more: &[&str],
) -> Remainder {
    let mut runs = Vec::new();
    for (threads, block) in [
        (1, None),
        (2, None),
        (worker_count(), None),
        (2, Some(48)),
        (3, Some(700)),
    ] {
        runs.push(remainder_at(
            format, expected, actual, input, more, threads, block,
        ));
    }
    for run in &runs[1..] {
        assert_eq!(run, &runs[0]);
    }
    texts_of(format, &runs[0])
}

pub fn counter(counters: &str, name: &str) -> String {
    let key = format!("\"{name}\": ");
    let at = counters
        .find(&key)
        .unwrap_or_else(|| panic!("{name} in {counters}"))
        + key.len();
    counters[at..]
        .lines()
        .next()
        .unwrap()
        .trim_end_matches(',')
        .to_string()
}

/// The counts of the block `name` of the counts, in the text's order.
pub fn block(counters: &str, name: &str) -> Vec<(String, u64)> {
    let from = counters.find(&format!("\"{name}\": {{")).unwrap();
    let block = &counters[from..];
    let block = &block[block.find('{').unwrap() + 1..block.find('}').unwrap()];
    block
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let (name, n) = line.trim().trim_end_matches(',').rsplit_once(": ").unwrap();
            (name.trim_matches('"').to_string(), n.parse().unwrap())
        })
        .collect()
}

/// The `kept_populated` counts, in the text's order.
pub fn populated(counters: &str) -> Vec<(String, u64)> {
    block(counters, "kept_populated")
}

pub fn pairs(expected: &[(&str, u64)]) -> Vec<(String, u64)> {
    expected
        .iter()
        .map(|(n, c)| ((*n).to_string(), *c))
        .collect()
}

pub fn line_runs(expected: &str, actual: &str) -> String {
    format!("{{\"expected\": {expected}, \"actual\": {actual}}}")
}
