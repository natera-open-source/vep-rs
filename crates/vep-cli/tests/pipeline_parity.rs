// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Pipeline parity: the reader, coordinator and writer stages must produce the
//! same bytes whatever the batch size and thread count, keep passthrough lines
//! in input order, and fail cleanly rather than hang when a stage errors.

mod common;

use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use common::{corpora, Corpus};

/// Output of one `vep` run: exit status, stdout, stderr and the output file's
/// bytes (empty when the run produced none).
struct Run {
    success: bool,
    stderr: String,
    output: Vec<u8>,
}

fn vep(input: &Path, json_cache: &Path, assembly: &str, out: &Path, args: &[&str]) -> Run {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    cmd.arg("-i")
        .arg(input)
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(json_cache)
        .args(["--species", "homo_sapiens", "--assembly", assembly])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .args(args);
    let result = cmd.output().expect("vep binary runs");
    Run {
        success: result.status.success(),
        stderr: String::from_utf8_lossy(&result.stderr).into_owned(),
        output: fs::read(out).unwrap_or_default(),
    }
}

fn run_corpus(corpus: &Corpus, out: &Path, args: &[&str]) -> Run {
    vep(
        &corpus.dir.join("variants.vcf"),
        &corpus.dir.join("json_cache"),
        &corpus.assembly,
        out,
        args,
    )
}

/// The output without the header lines that describe the run rather than the
/// annotation: the `## Output produced at` and `## VEP command-line:` lines of
/// the default and tab formats, the `##VEP-command-line=` line and the
/// `time="..."` attribute of the `##VEP` line of the VCF format.
fn without_run_headers(output: &[u8]) -> Vec<u8> {
    let text = String::from_utf8_lossy(output);
    let mut kept = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        if line.starts_with("## Output produced at")
            || line.starts_with("## VEP command-line:")
            || line.starts_with("##VEP-command-line=")
        {
            continue;
        }
        if line.starts_with("##VEP=") {
            if let Some(start) = line.find(" time=\"") {
                let rest = &line[start + 7..];
                if let Some(end) = rest.find('"') {
                    kept.push_str(&line[..start]);
                    kept.push_str(&rest[end + 1..]);
                    continue;
                }
            }
        }
        kept.push_str(line);
    }
    kept.into_bytes()
}

const FORMATS: [(&str, &[&str]); 4] = [
    ("default", &[]),
    ("tab", &["--tab"]),
    ("vcf", &["--vcf"]),
    ("json", &["--json"]),
];

#[test]
fn every_format_is_byte_identical_across_buffer_sizes() {
    let tmp = tempfile::tempdir().unwrap();
    let corpora = corpora();
    assert!(!corpora.is_empty(), "no corpus found under tests/golden");
    for corpus in &corpora {
        for (name, flags) in FORMATS {
            let mut reference: Option<Vec<u8>> = None;
            for buffer in ["5000", "7", "1"] {
                let out = tmp
                    .path()
                    .join(format!("{}_{name}_{buffer}.out", corpus.assembly));
                let mut args: Vec<&str> = flags.to_vec();
                args.extend(["--buffer_size", buffer, "--fork", "4"]);
                let run = run_corpus(corpus, &out, &args);
                assert!(
                    run.success,
                    "{}/{name} buffer {buffer}: {}",
                    corpus.assembly, run.stderr
                );
                let body = without_run_headers(&run.output);
                assert!(
                    !body.is_empty(),
                    "{}/{name} buffer {buffer}: empty output",
                    corpus.assembly
                );
                match &reference {
                    None => reference = Some(body),
                    Some(expected) => assert_same(
                        expected,
                        &body,
                        &format!(
                            "{}/{name}: --buffer_size {buffer} output differs from --buffer_size 5000",
                            corpus.assembly
                        ),
                    ),
                }
            }
        }
    }
}

/// Asserts two outputs are byte-identical, naming the first differing line.
fn assert_same(expected: &[u8], actual: &[u8], what: &str) {
    if expected == actual {
        return;
    }
    let (e, a) = (
        String::from_utf8_lossy(expected),
        String::from_utf8_lossy(actual),
    );
    let first_diff = e
        .lines()
        .zip(a.lines())
        .enumerate()
        .find(|(_, (x, y))| x != y)
        .map(|(i, (x, y))| format!("line {}:\n  expected: {x}\n  actual:   {y}", i + 1))
        .unwrap_or_else(|| {
            format!(
                "line counts differ: {} vs {}",
                e.lines().count(),
                a.lines().count()
            )
        });
    panic!("{what}\n{first_diff}");
}

#[test]
fn fork_1_matches_fork_4_on_every_format() {
    let tmp = tempfile::tempdir().unwrap();
    for corpus in &corpora() {
        for (name, flags) in FORMATS {
            let mut outputs = Vec::new();
            for fork in ["1", "4"] {
                let out = tmp
                    .path()
                    .join(format!("{}_{name}_fork{fork}.out", corpus.assembly));
                let mut args: Vec<&str> = flags.to_vec();
                args.extend(["--buffer_size", "500", "--fork", fork]);
                let run = run_corpus(corpus, &out, &args);
                assert!(
                    run.success,
                    "{}/{name} fork {fork}: {}",
                    corpus.assembly, run.stderr
                );
                outputs.push(without_run_headers(&run.output));
            }
            assert_same(
                &outputs[0],
                &outputs[1],
                &format!(
                    "{}/{name}: --fork 1 and --fork 4 outputs differ",
                    corpus.assembly
                ),
            );
        }
    }
}

/// The corpus VCF with a non-variant (`ALT=.`) copy of every third record
/// inserted after it, so `.` lines fall on and around batch boundaries at
/// every small buffer size.
fn with_non_variant_lines(corpus: &Corpus, dir: &Path) -> (std::path::PathBuf, Vec<String>) {
    let text = fs::read_to_string(corpus.dir.join("variants.vcf")).unwrap();
    let mut out = String::new();
    let mut expected_order: Vec<String> = Vec::new();
    let mut data_seen = 0usize;
    for line in text.lines() {
        if line.starts_with('#') {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        out.push_str(line);
        out.push('\n');
        expected_order.push(record_key(line));
        data_seen += 1;
        if data_seen.is_multiple_of(3) {
            let fields: Vec<&str> = line.split('\t').collect();
            let mut nv = fields.clone();
            nv[2] = "non_variant";
            nv[4] = ".";
            let nv_line = nv.join("\t");
            expected_order.push(record_key(&nv_line));
            out.push_str(&nv_line);
            out.push('\n');
        }
    }
    let path = dir.join(format!("non_variant_{}.vcf", corpus.assembly));
    fs::write(&path, out).unwrap();
    (path, expected_order)
}

/// `CHROM:POS:ID:REF:ALT` of a VCF data line.
fn record_key(line: &str) -> String {
    let f: Vec<&str> = line.split('\t').collect();
    format!("{}:{}:{}:{}:{}", f[0], f[1], f[2], f[3], f[4])
}

#[test]
fn allow_non_variant_lines_stay_in_input_order_at_every_buffer_size() {
    let tmp = tempfile::tempdir().unwrap();
    for corpus in &corpora() {
        let (input, expected_order) = with_non_variant_lines(corpus, tmp.path());
        let mut reference: Option<Vec<u8>> = None;
        for buffer in ["5000", "7", "3", "2", "1"] {
            let out = tmp
                .path()
                .join(format!("nv_{}_{buffer}.vcf", corpus.assembly));
            let run = vep(
                &input,
                &corpus.dir.join("json_cache"),
                &corpus.assembly,
                &out,
                &[
                    "--vcf",
                    "--allow_non_variant",
                    "--buffer_size",
                    buffer,
                    "--fork",
                    "4",
                ],
            );
            assert!(
                run.success,
                "{} buffer {buffer}: {}",
                corpus.assembly, run.stderr
            );
            let body = without_run_headers(&run.output);
            let text = String::from_utf8_lossy(&body);
            let order: Vec<String> = text
                .lines()
                .filter(|l| !l.starts_with('#'))
                .map(record_key)
                .collect();
            assert_eq!(
                order, expected_order,
                "{} buffer {buffer}: records are not in input order",
                corpus.assembly
            );
            let non_variant_lines = text
                .lines()
                .filter(|l| !l.starts_with('#') && l.split('\t').nth(2) == Some("non_variant"))
                .count();
            let expected_non_variant = expected_order
                .iter()
                .filter(|k| k.contains(":non_variant:"))
                .count();
            assert!(expected_non_variant > 0);
            assert_eq!(non_variant_lines, expected_non_variant);
            match &reference {
                None => reference = Some(body),
                Some(expected) => assert_same(
                    expected,
                    &body,
                    &format!("{}: --buffer_size {buffer} output differs", corpus.assembly),
                ),
            }
        }
    }
}

/// Runs `cmd` and returns its output, killing it after `limit`.
fn output_within(mut cmd: Command, limit: Duration) -> std::process::Output {
    cmd.stdout(Stdio::piped()).stderr(Stdio::piped());
    let mut child = cmd.spawn().expect("vep binary spawns");
    let started = Instant::now();
    loop {
        if let Some(_status) = child.try_wait().expect("try_wait") {
            return child.wait_with_output().expect("wait_with_output");
        }
        if started.elapsed() > limit {
            let _ = child.kill();
            let _ = child.wait();
            panic!("vep did not exit within {limit:?}: the pipeline hung");
        }
        std::thread::sleep(Duration::from_millis(20));
    }
}

#[test]
fn dont_skip_reports_the_malformed_line_in_a_later_batch_and_exits() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    let text = fs::read_to_string(corpus.dir.join("variants.vcf")).unwrap();
    let header: Vec<&str> = text.lines().filter(|l| l.starts_with('#')).collect();
    let data: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(
        data.len() >= 30,
        "corpus has {} records; need 30",
        data.len()
    );
    let mut lines: Vec<String> = header.iter().map(|l| l.to_string()).collect();
    lines.extend(data[..25].iter().map(|l| l.to_string()));
    // Line 26 of the data, in the third batch of 10, with a non-numeric POS.
    let bad_line_number = header.len() + 26;
    lines.push("21\tnot_a_position\t.\tA\tG\t.\tPASS\t.".to_string());
    lines.extend(data[26..30].iter().map(|l| l.to_string()));
    let input = tmp.path().join("malformed.vcf");
    fs::write(&input, lines.join("\n") + "\n").unwrap();
    let out = tmp.path().join("malformed.out");

    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    cmd.arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.dir.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", &corpus.assembly])
        .args(["--force_overwrite", "--no_stats", "--quiet", "--dont_skip"])
        .args(["--buffer_size", "10", "--fork", "4"]);
    let output = output_within(cmd, Duration::from_secs(120));
    assert!(
        !output.status.success(),
        "--dont_skip must fail on a malformed line"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains(&format!("line {bad_line_number}: parse error")),
        "stderr does not name line {bad_line_number}:\n{stderr}"
    );
    assert!(
        stderr.contains("Invalid POS"),
        "stderr lacks the parse error:\n{stderr}"
    );
}

#[cfg(unix)]
#[test]
fn unwritable_output_fails_instead_of_hanging() {
    let corpus = &corpora()[0];
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_vep"));
    cmd.arg("-i")
        .arg(corpus.dir.join("variants.vcf"))
        .arg("-o")
        .arg("/dev/full")
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.dir.join("json_cache"))
        .args(["--species", "homo_sapiens", "--assembly", &corpus.assembly])
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .args(["--buffer_size", "50", "--fork", "4"]);
    let output = output_within(cmd, Duration::from_secs(120));
    assert!(!output.status.success(), "writing to /dev/full must fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("output"),
        "stderr does not name the output failure:\n{stderr}"
    );
}
