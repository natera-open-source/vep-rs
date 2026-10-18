// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Pipeline parity: the reader, coordinator and writer stages must produce the
//! same bytes whatever the batch size and thread count, keep passthrough lines
//! in input order, and fail cleanly rather than hang when a stage errors.

mod common;

use std::collections::BTreeSet;
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
        &corpus.cache_dir(),
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

/// The default-format output with everything `LoFTEE` adds removed: its
/// `## LoF... :` header lines and every `LoF` key of the Extra column; a column
/// left with no key reads `-`.
fn without_plugin_keys(output: &[u8]) -> Vec<u8> {
    let text = String::from_utf8_lossy(output);
    let mut kept = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        if line.starts_with('#') {
            if !line.starts_with("## LoF") {
                kept.push_str(line);
            }
            continue;
        }
        let (body, newline) = match line.strip_suffix('\n') {
            Some(body) => (body, "\n"),
            None => (line, ""),
        };
        let (columns, extra_column) = body.rsplit_once('\t').unwrap_or((body, ""));
        let extra: Vec<&str> = extra_column
            .split(';')
            .filter(|kv| !kv.starts_with("LoF"))
            .collect();
        kept.push_str(columns);
        kept.push('\t');
        if extra.is_empty() {
            kept.push('-');
        } else {
            kept.push_str(&extra.join(";"));
        }
        kept.push_str(newline);
    }
    kept.into_bytes()
}

/// The corpus VCF with a breakend record appended whose two ends lie in the
/// chromosome 21 centromere, megabases from any transcript of the corpus
/// cache: no allele of it reaches a transcript, so its intergenic rows exist
/// only once the record rules have run.
fn with_intergenic_breakend(corpus: &Corpus, dir: &Path) -> std::path::PathBuf {
    let mut text = fs::read_to_string(corpus.dir.join("variants.vcf")).unwrap();
    text.push_str("21\t12000000\tintergenic_bnd\tN\tN[21:12500000[\t.\t.\tSVTYPE=BND\n");
    let path = dir.join("intergenic_breakend.vcf");
    fs::write(&path, text).unwrap();
    path
}

/// The record step runs at the same point in both coordinator paths: the fused
/// path a run without plugins takes and the batched path a plugin forces
/// (`LoFTEE` loads with no data file and adds only its own `LoF` keys). With
/// `--flag_pick` the two paths mark the same row of every record, and with
/// `--most_severe` they squash every record to the same row, once the plugin's
/// keys are removed from the Extra column. The corpus's sequence variants read
/// the same whether the step runs before or after the record rules; the
/// intergenic breakend record has no row before them, so its flagged row (one
/// row under `--most_severe`) shows the step ran after.
#[test]
fn pick_family_runs_at_the_same_point_in_both_coordinator_paths() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-pick-summary")
        .expect("the pick-summary corpus is committed");
    let input = with_intergenic_breakend(&corpus, tmp.path());
    for flag in ["--flag_pick", "--most_severe"] {
        let fused_out = tmp.path().join(format!("fused{flag}.out"));
        let fused = vep(
            &input,
            &corpus.cache_dir(),
            &corpus.assembly,
            &fused_out,
            &[flag, "--buffer_size", "7", "--fork", "4"],
        );
        assert!(fused.success, "{flag} without plugin: {}", fused.stderr);
        let fused_text = String::from_utf8_lossy(&fused.output);
        let breakend_rows: Vec<&str> = fused_text
            .lines()
            .filter(|l| l.starts_with("intergenic_bnd\t"))
            .collect();
        let picked = breakend_rows
            .iter()
            .filter(|l| {
                l.rsplit('\t')
                    .next()
                    .unwrap()
                    .split(';')
                    .any(|kv| kv == "PICK=1")
            })
            .count();
        match flag {
            "--flag_pick" => assert!(
                !breakend_rows.is_empty() && picked == 1,
                "{flag}: the breakend record's rows {breakend_rows:?} carry {picked} PICK flags, not one"
            ),
            _ => assert_eq!(
                breakend_rows.len(),
                1,
                "{flag}: the breakend record keeps {breakend_rows:?}, not one row"
            ),
        }
        let batched_out = tmp.path().join(format!("batched{flag}.out"));
        let batched = vep(
            &input,
            &corpus.cache_dir(),
            &corpus.assembly,
            &batched_out,
            &[
                flag,
                "--plugin",
                "LoFTEE",
                "--buffer_size",
                "7",
                "--fork",
                "4",
            ],
        );
        assert!(batched.success, "{flag} with plugin: {}", batched.stderr);
        assert_same(
            &without_run_headers(&fused.output),
            &without_plugin_keys(&without_run_headers(&batched.output)),
            &format!("{flag}: the batched path differs from the fused path"),
        );
    }
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

/// The records of a corpus the reference dropped before annotation and vep-rs
/// writes nothing for either: `reference_rows` 0 with no `vep_rs_only_tuples`
/// entry naming the record (the entries name the dropped records vep-rs still
/// annotates, a documented divergence).
fn dropped_by_both(corpus: &Corpus) -> BTreeSet<usize> {
    let manifest = &corpus.manifest;
    let named: BTreeSet<usize> = manifest["vep_rs_only_tuples"]
        .as_array()
        .into_iter()
        .flatten()
        .flat_map(|t| t["record_indices"].as_array().into_iter().flatten())
        .filter_map(|r| r.as_u64().map(|r| r as usize))
        .collect();
    manifest["records"]
        .as_array()
        .into_iter()
        .flatten()
        .enumerate()
        .filter(|(i, r)| r["reference_rows"].as_u64() == Some(0) && !named.contains(i))
        .map(|(i, _)| i)
        .collect()
}

/// The corpus VCF with a non-variant (`ALT=.`) copy of every third record
/// inserted after it, so `.` lines fall on and around batch boundaries at
/// every small buffer size. A record both engines drop ([`dropped_by_both`])
/// stays in the input, is expected in no output, and gets no copy.
fn with_non_variant_lines(corpus: &Corpus, dir: &Path) -> (std::path::PathBuf, Vec<String>) {
    let text = fs::read_to_string(corpus.dir.join("variants.vcf")).unwrap();
    let dropped = dropped_by_both(corpus);
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
        let record = data_seen;
        data_seen += 1;
        if dropped.contains(&record) {
            continue;
        }
        expected_order.push(record_key(line));
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
                &corpus.cache_dir(),
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
    // Data record 26, in the third batch of 10, with a non-numeric POS; the
    // number reported is the record's, as Ensembl VEP counts lines
    // (Parser/VCF.pm:228), header lines not counted.
    let bad_line_number = 26;
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
        .arg(corpus.cache_dir())
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
        .arg(corpus.cache_dir())
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
