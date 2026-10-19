// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Pipeline parity: the reader, coordinator and writer stages must produce the
//! same bytes whatever the batch size, thread count, pipeline mode, pipeline
//! depth, task size, writer strategy and decompression thread count, keep
//! passthrough lines in input order, and fail cleanly rather than hang when a
//! stage errors.

mod common;

use std::collections::BTreeSet;
use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use common::{corpora, run_vep, runs, Corpus};

/// Output of one `vep` run: exit status, stderr and the output bytes (the
/// output file's, or standard output's for a run to `STDOUT`; empty when the
/// run produced none).
struct Run {
    success: bool,
    stderr: String,
    output: Vec<u8>,
}

/// The pipeline modes a case runs under: the batched coordinator first, the
/// reference the comparisons read, then streaming, the default.
const MODES: [(&str, &[&str]); 2] = [
    ("batched", &["--pipeline_mode", "batched"]),
    ("streaming", &["--pipeline_mode", "streaming"]),
];

/// The `vep` invocation every case runs: `input` to `out` (a path or `STDOUT`)
/// against `json_cache` and `assembly`, quiet and overwriting, then `args`.
fn vep_command(
    input: &Path,
    json_cache: &Path,
    assembly: &str,
    out: &Path,
    args: &[&str],
) -> Command {
    let mut cmd = common::isolated_vep_anywhere();
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
    cmd
}

/// Runs `vep` with `-o out`, where `out` is a path or `STDOUT`.
fn vep(input: &Path, json_cache: &Path, assembly: &str, out: &Path, args: &[&str]) -> Run {
    run_command(vep_command(input, json_cache, assembly, out, args), out)
}

/// Runs a `vep` command that writes to `out`, a path or `STDOUT`. A run to
/// `STDOUT` writes its warnings file into the working directory, so it runs in
/// a directory of its own rather than the crate's.
fn run_command(mut cmd: Command, out: &Path) -> Run {
    let to_stdout = out.as_os_str() == "STDOUT";
    let scratch = tempfile::tempdir().unwrap();
    if to_stdout {
        cmd.current_dir(scratch.path());
    }
    let result = cmd.output().expect("vep binary runs");
    let output = if to_stdout {
        result.stdout
    } else {
        fs::read(out).unwrap_or_default()
    };
    Run {
        success: result.status.success(),
        stderr: String::from_utf8_lossy(&result.stderr).into_owned(),
        output,
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

const WRITER_STRATEGIES: [&str; 3] = ["buffered", "direct", "vectored"];

/// Every format, in both pipeline modes, at every buffer size, is the bytes of
/// the default run (batched, `--buffer_size 5000`).
#[test]
fn every_format_is_byte_identical_across_buffer_sizes_and_modes() {
    let tmp = tempfile::tempdir().unwrap();
    let corpora = corpora();
    assert!(!corpora.is_empty(), "no corpus found under tests/golden");
    for corpus in &corpora {
        for (name, flags) in FORMATS {
            let mut reference: Option<Vec<u8>> = None;
            for (mode, mode_flags) in MODES {
                for buffer in ["5000", "7", "1"] {
                    let out = tmp
                        .path()
                        .join(format!("{}_{name}_{mode}_{buffer}.out", corpus.name));
                    let mut args: Vec<&str> = flags.to_vec();
                    args.extend(["--buffer_size", buffer, "--fork", "4"]);
                    args.extend(mode_flags);
                    let run = run_corpus(corpus, &out, &args);
                    assert!(
                        run.success,
                        "{}/{name} {mode} buffer {buffer}: {}",
                        corpus.name, run.stderr
                    );
                    let body = without_run_headers(&run.output);
                    assert!(
                        !body.is_empty(),
                        "{}/{name} {mode} buffer {buffer}: empty output",
                        corpus.name
                    );
                    match &reference {
                        None => reference = Some(body),
                        Some(expected) => assert_same(
                            expected,
                            &body,
                            &format!(
                                "{}/{name}: {mode} --buffer_size {buffer} output differs from the batched --buffer_size 5000 output",
                                corpus.name
                            ),
                        ),
                    }
                }
            }
        }
    }
}

/// A run that reads the cache's variation files takes the batched path
/// whatever `--pipeline_mode` says (`pipeline::Coordinator::streaming`: the
/// variation prefetch holds one batch's regions at a time), so asking for
/// `streaming` with `--check_existing` and the frequency and citation columns
/// on writes the batched run's bytes: every corpus whose cache carries a
/// `variations` directory.
#[test]
fn variation_lookups_write_the_batched_bytes_whatever_the_mode_asks() {
    let tmp = tempfile::tempdir().unwrap();
    let corpora: Vec<Corpus> = corpora()
        .into_iter()
        .filter(|c| c.cache_dir().join("variations").is_dir())
        .collect();
    assert!(
        !corpora.is_empty(),
        "no corpus under tests/golden carries variation files"
    );
    for corpus in &corpora {
        let mut reference: Option<Vec<u8>> = None;
        for (mode, mode_flags) in MODES {
            let out = tmp
                .path()
                .join(format!("{}_variation_{mode}.out", corpus.name));
            let mut args = vec![
                "--check_existing",
                "--af",
                "--af_gnomadg",
                "--pubmed",
                "--fork",
                "4",
            ];
            args.extend(mode_flags);
            let run = run_corpus(corpus, &out, &args);
            assert!(run.success, "{}/{mode}: {}", corpus.name, run.stderr);
            let body = without_run_headers(&run.output);
            assert!(!body.is_empty(), "{}/{mode}: empty output", corpus.name);
            match &reference {
                None => reference = Some(body),
                Some(expected) => assert_same(
                    expected,
                    &body,
                    &format!(
                        "{}: {mode} output with variation lookups differs from the batched output",
                        corpus.name
                    ),
                ),
            }
        }
    }
}

/// Every manifest run of every corpus, its own flags and each of its flag sets,
/// writes the same default-format bytes under `streaming` as under `batched`.
#[test]
fn every_manifest_run_is_byte_identical_in_both_modes() {
    let tmp = tempfile::tempdir().unwrap();
    let corpora = corpora();
    assert!(!corpora.is_empty(), "no corpus found under tests/golden");
    let mut compared = 0usize;
    for corpus in &corpora {
        for run in runs(corpus) {
            let batched = run_vep(&run, "default", tmp.path(), &["--pipeline_mode", "batched"]);
            let streaming = run_vep(
                &run,
                "default",
                tmp.path(),
                &["--pipeline_mode", "streaming"],
            );
            assert_same(
                &without_run_headers(batched.as_bytes()),
                &without_run_headers(streaming.as_bytes()),
                &format!("{}: streaming output differs from batched", run.label()),
            );
            compared += 1;
        }
    }
    assert!(
        compared > corpora.len(),
        "no flag set ran beside the corpora's own runs"
    );
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
fn fork_1_matches_fork_4_on_every_format_in_both_modes() {
    let tmp = tempfile::tempdir().unwrap();
    for corpus in &corpora() {
        for (name, flags) in FORMATS {
            let mut outputs = Vec::new();
            for (mode, mode_flags) in MODES {
                for fork in ["1", "4"] {
                    let out = tmp
                        .path()
                        .join(format!("{}_{name}_{mode}_fork{fork}.out", corpus.name));
                    let mut args: Vec<&str> = flags.to_vec();
                    args.extend(["--buffer_size", "500", "--fork", fork]);
                    args.extend(mode_flags);
                    let run = run_corpus(corpus, &out, &args);
                    assert!(
                        run.success,
                        "{}/{name} {mode} fork {fork}: {}",
                        corpus.name, run.stderr
                    );
                    outputs.push((mode, fork, without_run_headers(&run.output)));
                }
            }
            for (mode, fork, output) in &outputs[1..] {
                assert_same(
                    &outputs[0].2,
                    output,
                    &format!(
                        "{}/{name}: {mode} --fork {fork} differs from batched --fork 1",
                        corpus.name
                    ),
                );
            }
        }
    }
}

/// Every writer strategy, at every buffer size, to a file and to standard
/// output, writes the bytes of the default run (buffered, the built-in size, to
/// a file).
#[test]
fn every_writer_strategy_writes_the_default_bytes() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    for (name, flags) in [FORMATS[0], FORMATS[2]] {
        let reference_out = tmp.path().join(format!("writer_{name}_reference.out"));
        let mut args: Vec<&str> = flags.to_vec();
        args.extend(["--buffer_size", "100", "--fork", "4"]);
        let reference = run_corpus(corpus, &reference_out, &args);
        assert!(reference.success, "{name} reference: {}", reference.stderr);
        let reference = without_run_headers(&reference.output);
        for strategy in WRITER_STRATEGIES {
            for buffer in ["0", "4096", "1"] {
                for destination in ["file", "STDOUT"] {
                    let file_out = tmp
                        .path()
                        .join(format!("writer_{name}_{strategy}_{buffer}.out"));
                    let out: &Path = if destination == "file" {
                        file_out.as_path()
                    } else {
                        Path::new("STDOUT")
                    };
                    let mut args: Vec<&str> = flags.to_vec();
                    args.extend(["--buffer_size", "100", "--fork", "4"]);
                    args.extend(["--writer_strategy", strategy]);
                    args.extend(["--output_buffer_size", buffer]);
                    let run = run_corpus(corpus, out, &args);
                    assert!(
                        run.success,
                        "{name} {strategy} {buffer} {destination}: {}",
                        run.stderr
                    );
                    assert_same(
                        &reference,
                        &without_run_headers(&run.output),
                        &format!(
                            "{name}: --writer_strategy {strategy} --output_buffer_size {buffer} to {destination} differs from the default"
                        ),
                    );
                }
            }
        }
    }
}

/// `--compress_output gzip` and `bgzip`, in both modes and under every writer
/// strategy, decompress to the bytes of the plain default run in the default and
/// VCF formats; the compressed bytes themselves follow the encoder's write
/// boundaries and are not compared.
#[test]
fn compressed_output_decompresses_to_the_plain_bytes_in_both_modes() {
    use std::io::Read;

    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    for (name, flags) in [FORMATS[0], FORMATS[2]] {
        let reference_out = tmp.path().join(format!("plain_{name}.out"));
        let mut args: Vec<&str> = flags.to_vec();
        args.extend(["--buffer_size", "100", "--fork", "4"]);
        let reference = run_corpus(corpus, &reference_out, &args);
        assert!(reference.success, "{name} reference: {}", reference.stderr);
        let reference = without_run_headers(&reference.output);
        for (mode, mode_flags) in MODES {
            for codec in ["gzip", "bgzip"] {
                for strategy in WRITER_STRATEGIES {
                    let out = tmp.path().join(format!(
                        "compressed_{name}_{mode}_{codec}_{strategy}.out.gz"
                    ));
                    let mut args: Vec<&str> = flags.to_vec();
                    args.extend(["--buffer_size", "100", "--fork", "4"]);
                    args.extend(mode_flags);
                    args.extend(["--compress_output", codec]);
                    args.extend(["--writer_strategy", strategy]);
                    let run = run_corpus(corpus, &out, &args);
                    assert!(
                        run.success,
                        "{name} {mode} {codec} {strategy}: {}",
                        run.stderr
                    );
                    let mut plain = Vec::new();
                    flate2::read::MultiGzDecoder::new(&run.output[..])
                        .read_to_end(&mut plain)
                        .unwrap_or_else(|e| {
                            panic!("{name} {mode} {codec} {strategy}: not a gzip stream: {e}")
                        });
                    assert_same(
                        &reference,
                        &without_run_headers(&plain),
                        &format!(
                            "{name}: {mode} --compress_output {codec} --writer_strategy {strategy} differs from the plain run once decompressed"
                        ),
                    );
                }
            }
        }
    }
}

/// Whether `duckdb`, which finalizes a Parquet output, is on the PATH.
fn duckdb_available() -> bool {
    Command::new("duckdb")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|status| status.success())
}

/// A Parquet run, in both modes and under every writer strategy, writes the
/// same TSV intermediate (kept with `--keep_intermediate`) and the same rows
/// once DuckDB has converted it; skipped when `duckdb` is not on the PATH.
#[test]
fn parquet_output_is_identical_in_both_modes_and_strategies() {
    if !duckdb_available() {
        eprintln!("skipped: duckdb not on PATH");
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    let mut reference: Option<(Vec<u8>, String)> = None;
    for (mode, mode_flags) in MODES {
        for strategy in WRITER_STRATEGIES {
            let out = tmp
                .path()
                .join(format!("parquet_{mode}_{strategy}.parquet"));
            let mut args = vec![
                "--output_format",
                "parquet",
                "--parquet_shape",
                "flat",
                "--keep_intermediate",
                "--buffer_size",
                "100",
                "--fork",
                "4",
            ];
            args.extend(mode_flags);
            args.extend(["--writer_strategy", strategy]);
            let run = run_corpus(corpus, &out, &args);
            assert!(run.success, "{mode} {strategy}: {}", run.stderr);
            let intermediate = fs::read(format!("{}.tsv.tmp", out.display())).unwrap();
            assert!(
                !intermediate.is_empty(),
                "{mode} {strategy}: empty intermediate"
            );
            let sql = format!(
                "SELECT * FROM read_parquet('{}/**/*.parquet', hive_partitioning = true) ORDER BY ALL",
                out.display()
            );
            let query = Command::new("duckdb")
                .args(["-csv", "-c", &sql])
                .output()
                .expect("duckdb runs");
            assert!(
                query.status.success(),
                "{mode} {strategy}: duckdb failed: {}",
                String::from_utf8_lossy(&query.stderr)
            );
            let rows = String::from_utf8_lossy(&query.stdout).into_owned();
            assert!(
                rows.lines().count() > 1,
                "{mode} {strategy}: no Parquet rows"
            );
            match &reference {
                None => reference = Some((intermediate, rows)),
                Some((expected_intermediate, expected_rows)) => {
                    assert_same(
                        expected_intermediate,
                        &intermediate,
                        &format!("{mode} {strategy}: the Parquet intermediate differs from batched buffered"),
                    );
                    assert_same(
                        expected_rows.as_bytes(),
                        rows.as_bytes(),
                        &format!(
                            "{mode} {strategy}: the Parquet rows differ from batched buffered"
                        ),
                    );
                }
            }
        }
    }
}

/// Every pipeline depth and task size, in both modes, writes the bytes of the
/// default run.
#[test]
fn pipeline_depth_and_task_size_leave_the_output_unchanged() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    for (name, flags) in [FORMATS[0], FORMATS[2]] {
        let reference_out = tmp.path().join(format!("depth_{name}_reference.out"));
        let mut args: Vec<&str> = flags.to_vec();
        args.extend(["--buffer_size", "50", "--fork", "4"]);
        let reference = run_corpus(corpus, &reference_out, &args);
        assert!(reference.success, "{name} reference: {}", reference.stderr);
        let reference = without_run_headers(&reference.output);
        for (mode, mode_flags) in MODES {
            for depth in ["1", "2", "8"] {
                for records in ["1", "7", "64"] {
                    let out = tmp
                        .path()
                        .join(format!("depth_{name}_{mode}_{depth}_{records}.out"));
                    let mut args: Vec<&str> = flags.to_vec();
                    args.extend(["--buffer_size", "50", "--fork", "4"]);
                    args.extend(mode_flags);
                    args.extend(["--pipeline_depth", depth]);
                    args.extend(["--render_task_records", records]);
                    let run = run_corpus(corpus, &out, &args);
                    assert!(
                        run.success,
                        "{name} {mode} depth {depth} records {records}: {}",
                        run.stderr
                    );
                    assert_same(
                        &reference,
                        &without_run_headers(&run.output),
                        &format!(
                            "{name}: {mode} --pipeline_depth {depth} --render_task_records {records} differs from the default"
                        ),
                    );
                }
            }
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

/// Every corpus runs under its own flags and reference, so a record its
/// manifest marks dropped is dropped by the same checks here (a reference
/// allele `--check_ref` rejects is kept without the flag).
#[test]
fn allow_non_variant_lines_stay_in_input_order_at_every_buffer_size_in_both_modes() {
    let tmp = tempfile::tempdir().unwrap();
    for corpus in &corpora() {
        let (input, expected_order) = with_non_variant_lines(corpus, tmp.path());
        let mut reference: Option<Vec<u8>> = None;
        for (mode, mode_flags) in MODES {
            for buffer in ["5000", "7", "3", "2", "1"] {
                let out = tmp
                    .path()
                    .join(format!("nv_{}_{mode}_{buffer}.vcf", corpus.name));
                let mut corpus_run = common::default_run(corpus);
                corpus_run.input = input.display().to_string();
                let mut cmd = common::vep_command(&corpus_run, tmp.path());
                cmd.arg("-o")
                    .arg(&out)
                    .args([
                        "--vcf",
                        "--allow_non_variant",
                        "--buffer_size",
                        buffer,
                        "--fork",
                        "4",
                    ])
                    .args(mode_flags);
                let run = run_command(cmd, &out);
                assert!(
                    run.success,
                    "{} {mode} buffer {buffer}: {}",
                    corpus.name, run.stderr
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
                    "{} {mode} buffer {buffer}: records are not in input order",
                    corpus.name
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
                        &format!(
                            "{}: {mode} --buffer_size {buffer} output differs",
                            corpus.name
                        ),
                    ),
                }
            }
        }
    }
}

/// A copy of `vcf` at `out`, bgzipped with its tabix index beside it, one block
/// closed every `records_per_block` records; `padding` bytes of an extra INFO
/// field widen every record. A plain copy is written beside it too, and its
/// path returned.
fn bgzip_with_index(
    vcf: &Path,
    out: &Path,
    records_per_block: usize,
    padding: usize,
) -> std::path::PathBuf {
    use std::io::Write;

    use noodles::bgzf;
    use noodles::core::Position;
    use noodles::csi::binning_index::index::header::Builder as HeaderBuilder;
    use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;

    let text = fs::read_to_string(vcf).unwrap();
    let plain_path = out.with_extension("");
    let mut plain = fs::File::create(&plain_path).unwrap();
    let mut writer = bgzf::io::Writer::new(fs::File::create(out).unwrap());
    let mut indexer = noodles::tabix::index::Indexer::default();
    indexer.set_header(HeaderBuilder::vcf().build());
    let pad = format!("PAD={}", "x".repeat(padding));
    let mut records = 0usize;
    for line in text.lines() {
        if line.starts_with('#') {
            writeln!(plain, "{line}").unwrap();
            writeln!(writer, "{line}").unwrap();
            continue;
        }
        let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        if padding > 0 {
            fields[7] = if fields[7] == "." {
                pad.clone()
            } else {
                format!("{};{pad}", fields[7])
            };
        }
        let line = fields.join("\t");
        writeln!(plain, "{line}").unwrap();
        let start = writer.virtual_position();
        writeln!(writer, "{line}").unwrap();
        let pos: usize = fields[1].parse().unwrap();
        let end = pos + fields[3].len() - 1;
        indexer
            .add_record(
                &fields[0],
                Position::try_from(pos).unwrap(),
                Position::try_from(end).unwrap(),
                Chunk::new(start, writer.virtual_position()),
            )
            .unwrap();
        records += 1;
        if records.is_multiple_of(records_per_block) {
            writer.flush().unwrap();
        }
    }
    writer.finish().unwrap();
    let mut index_path = out.as_os_str().to_owned();
    index_path.push(".tbi");
    noodles::tabix::fs::write(index_path, &indexer.build()).unwrap();
    plain_path
}

/// A bgzipped, tabix-indexed copy of the corpus annotates to the plain file's
/// bytes under every form of `--decompression_threads`, in both modes, whether
/// its records inflate to a few hundred bytes or past the width that asks for
/// the inflate-bound pool; the probe arm asserts which side of the threshold
/// each copy's shape falls on and the count `auto` resolves there.
#[test]
fn bgzipped_input_matches_the_plain_file_under_every_decompression_thread_count() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    let cpus = vep_cli::config::host_parallelism();
    let bound = vep_cli::config::INFLATE_BOUND_BYTES_PER_RECORD;
    let host_capped = cpus.min(vep_cli::config::DECOMPRESSION_THREADS_AUTO_CAP);
    let inflate_bound = cpus.min(vep_cli::config::INFLATE_BOUND_INFLATERS);
    for (width, padding, auto_threads, inflated_bytes) in [
        ("narrow", 0usize, host_capped, 0..bound),
        ("wide", 8_192, inflate_bound, bound..u64::MAX),
    ] {
        let gz = tmp.path().join(format!("{width}.vcf.gz"));
        let plain = bgzip_with_index(&corpus.dir.join("variants.vcf"), &gz, 40, padding);
        // On a host with no more CPUs than the `0` rule's cap, and on an
        // architecture whose inflate-bound pool is that cap, both widths resolve
        // `auto` to one count, so the inflated width is asserted as well.
        let probe = vep(
            &gz,
            &corpus.cache_dir(),
            &corpus.assembly,
            Path::new("STDOUT"),
            &["--resolved_config", "--decompression_threads", "auto"],
        );
        assert!(probe.success, "{width} probe: {}", probe.stderr);
        let resolved = String::from_utf8_lossy(&probe.output);
        assert!(
            resolved.contains(&format!("\ndecompression_threads\t{auto_threads}\n"))
                && resolved.contains("\ninput_shape\tindexed\n"),
            "{width}: {resolved}"
        );
        let inflated: u64 = resolved
            .lines()
            .find_map(|line| line.strip_prefix("inflated_bytes_per_record\t"))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("{width}: {resolved}"));
        assert!(
            inflated_bytes.contains(&inflated),
            "{width}: {inflated} inflated bytes per record"
        );
        for (name, flags) in [FORMATS[0], FORMATS[2]] {
            let reference_out = tmp.path().join(format!("gz_{width}_{name}_reference.out"));
            let mut args: Vec<&str> = flags.to_vec();
            args.extend(["--buffer_size", "100", "--fork", "4"]);
            let reference = vep(
                &plain,
                &corpus.cache_dir(),
                &corpus.assembly,
                &reference_out,
                &args,
            );
            assert!(
                reference.success,
                "{width}/{name} reference: {}",
                reference.stderr
            );
            let reference = without_run_headers(&reference.output);
            assert!(!reference.is_empty());
            for (mode, mode_flags) in MODES {
                for threads in ["auto", "1", "0"] {
                    let out = tmp
                        .path()
                        .join(format!("gz_{width}_{name}_{mode}_{threads}.out"));
                    let mut args: Vec<&str> = flags.to_vec();
                    args.extend(["--buffer_size", "100", "--fork", "4"]);
                    args.extend(mode_flags);
                    args.extend(["--decompression_threads", threads]);
                    let run = vep(&gz, &corpus.cache_dir(), &corpus.assembly, &out, &args);
                    assert!(
                        run.success,
                        "{width}/{name} {mode} --decompression_threads {threads}: {}",
                        run.stderr
                    );
                    assert_same(
                        &reference,
                        &without_run_headers(&run.output),
                        &format!(
                            "{width}/{name}: {mode} --decompression_threads {threads} on the bgzipped copy differs from the plain file"
                        ),
                    );
                }
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

/// In both modes, a malformed line in a later batch is skipped with Ensembl
/// VEP's warning at the record's own number under `--dont_skip`, the run
/// completes with exit 0, and every other record is written: `--dont_skip`
/// keeps records that fail a check (Parser.pm:189), and a line that yields no
/// record has nothing to keep.
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
    lines.extend(data[..20].iter().map(|l| l.to_string()));
    // Data record 21, the first of the third batch of 10, with a non-numeric
    // POS; the number reported is the record's, as Ensembl VEP counts lines
    // (Parser/VCF.pm:228), header lines not counted.
    let bad_line_number = 21;
    lines.push("21\tnot_a_position\t.\tA\tG\t.\tPASS\t.".to_string());
    lines.extend(data[21..30].iter().map(|l| l.to_string()));
    let input = tmp.path().join("malformed.vcf");
    fs::write(&input, lines.join("\n") + "\n").unwrap();

    for (mode, mode_flags) in MODES {
        let out = tmp.path().join(format!("malformed_{mode}.out"));
        let mut args = vec!["--dont_skip", "--buffer_size", "10", "--fork", "4"];
        args.extend(mode_flags);
        let cmd = vep_command(&input, &corpus.cache_dir(), &corpus.assembly, &out, &args);
        let output = output_within(cmd, Duration::from_secs(120));
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(
            output.status.code(),
            Some(0),
            "{mode}: --dont_skip keeps the run going past a malformed line:\n{stderr}"
        );
        let warnings = fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap();
        assert_eq!(
            warnings.lines().collect::<Vec<_>>(),
            [format!(
                "WARNING: line {bad_line_number} skipped (21 not_a_position . A G . PASS .): Invalid start 'not_a_position' or end '0' coordinate"
            )],
            "{mode}: the malformed line is the one skip, at its record number"
        );
        let written: Vec<String> = fs::read_to_string(&out)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect();
        let expected_ids = |range: std::ops::Range<usize>| {
            data[range]
                .iter()
                .map(|l| l.split('\t').nth(2).unwrap().to_string())
                .filter(|id| id != ".")
        };
        for id in expected_ids(0..20).chain(expected_ids(21..30)) {
            assert!(written.contains(&id), "{mode}: record {id} is written");
        }
    }
}

/// The line numbers of a warnings file's `line N skipped` lines, in file order;
/// every line of the file is one.
fn warned_lines(warnings: &[u8]) -> Vec<u64> {
    String::from_utf8_lossy(warnings)
        .lines()
        .map(|line| {
            line.strip_prefix("WARNING: line ")
                .and_then(|rest| rest.split(' ').next())
                .and_then(|number| number.parse().ok())
                .unwrap_or_else(|| panic!("not a skip line: {line}"))
        })
        .collect()
}

/// The line numbers of a skipped-variants file's `Line N` rows, in file order,
/// after its header line.
fn skipped_rows(skipped: &[u8]) -> Vec<u64> {
    String::from_utf8_lossy(skipped)
        .lines()
        .skip(1)
        .map(|row| {
            row.strip_prefix("Line ")
                .and_then(|rest| rest.split_whitespace().next())
                .and_then(|number| number.parse().ok())
                .unwrap_or_else(|| panic!("not a skipped-variants row: {row}"))
        })
        .collect()
}

/// In both modes the warnings file and the skipped-variants file list the
/// skipped records in input order, and a streaming run lists them the same way
/// every time: the corpus's first 400 records with every third record's POS not
/// a number and every seventh record's ALT empty, in batches of 50 cut into
/// tasks for 8 workers, so the tasks of a batch finish out of input order.
#[test]
fn skipped_records_are_listed_in_input_order_in_both_modes() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    let text = fs::read_to_string(corpus.dir.join("variants.vcf")).unwrap();
    let data: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(
        data.len() >= 400,
        "corpus has {} records; need 400",
        data.len()
    );
    let mut lines: Vec<String> = text
        .lines()
        .filter(|l| l.starts_with('#'))
        .map(str::to_string)
        .collect();
    for (i, line) in data[..400].iter().enumerate() {
        // The record's number as Ensembl VEP counts lines (Parser/VCF.pm:228),
        // header lines not counted.
        let record = i + 1;
        let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        if record % 3 == 0 {
            fields[1] = format!("abc{record}");
        }
        if record % 7 == 0 {
            fields[4] = String::new();
        }
        lines.push(fields.join("\t"));
    }
    let input = tmp.path().join("skips.vcf");
    fs::write(&input, lines.join("\n") + "\n").unwrap();
    let expected_skips: Vec<u64> = (1..=400u64)
        .filter(|record| record % 3 == 0 || record % 7 == 0)
        .collect();

    /// The files one run wrote.
    struct SkipRun {
        name: String,
        output: Vec<u8>,
        warnings: Vec<u8>,
        skipped: Vec<u8>,
    }
    let mut runs: Vec<SkipRun> = Vec::new();
    for (mode, mode_flags) in MODES {
        let repeats = if mode == "streaming" { 3 } else { 1 };
        for repeat in 0..repeats {
            let name = format!("{mode} run {repeat}");
            let out = tmp.path().join(format!("skips_{mode}_{repeat}.out"));
            let warnings = tmp.path().join(format!("w_{mode}_{repeat}.txt"));
            let skipped = tmp.path().join(format!("s_{mode}_{repeat}.txt"));
            let mut args = vec![
                "--buffer_size",
                "50",
                "--fork",
                "8",
                "--warning_file",
                warnings.to_str().unwrap(),
                "--skipped_variants_file",
                skipped.to_str().unwrap(),
            ];
            args.extend(mode_flags);
            let run = vep(&input, &corpus.cache_dir(), &corpus.assembly, &out, &args);
            assert!(run.success, "{name}: {}", run.stderr);
            runs.push(SkipRun {
                name,
                output: without_run_headers(&run.output),
                warnings: fs::read(&warnings).unwrap(),
                skipped: fs::read(&skipped).unwrap(),
            });
        }
    }

    let batched = &runs[0];
    assert!(!batched.output.is_empty(), "batched: empty output");
    assert_eq!(warned_lines(&batched.warnings), expected_skips);
    assert_eq!(skipped_rows(&batched.skipped), expected_skips);
    for run in &runs[1..] {
        assert_same(
            &batched.output,
            &run.output,
            &format!("{}: output differs from batched", run.name),
        );
        assert_same(
            &batched.warnings,
            &run.warnings,
            &format!("{}: warnings file differs from batched", run.name),
        );
        assert_same(
            &batched.skipped,
            &run.skipped,
            &format!("{}: skipped-variants file differs from batched", run.name),
        );
    }
}

/// In both modes, a warnings file that cannot be written ends the run with the
/// same error line and the same output: the batch whose skip met the failure is
/// written and nothing after it. The skipped record (a chromosome the cache
/// lacks) opens the first batch of 10, so a clean run over the nine records
/// beside it is the reference for both.
#[test]
fn unwritable_warnings_file_ends_both_modes_at_the_same_batch() {
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
    lines.extend(data[..9].iter().map(|l| l.to_string()));
    let clean = tmp.path().join("clean.vcf");
    fs::write(&clean, lines.join("\n") + "\n").unwrap();
    let mut lines: Vec<String> = header.iter().map(|l| l.to_string()).collect();
    lines.push("Z\t100\t.\tA\tG\t.\tPASS\t.".to_string());
    lines.extend(data[..30].iter().map(|l| l.to_string()));
    let input = tmp.path().join("skipped_first.vcf");
    fs::write(&input, lines.join("\n") + "\n").unwrap();
    let warnings = tmp.path().join("missing").join("w.txt");

    let reference_out = tmp.path().join("clean.out");
    let reference = vep(
        &clean,
        &corpus.cache_dir(),
        &corpus.assembly,
        &reference_out,
        &["--buffer_size", "10", "--fork", "4"],
    );
    assert!(reference.success, "clean run: {}", reference.stderr);
    let reference = without_run_headers(&reference.output);
    assert!(!reference.is_empty());

    let mut error_lines: Vec<String> = Vec::new();
    for (mode, mode_flags) in MODES {
        let out = tmp.path().join(format!("skipped_first_{mode}.out"));
        let mut args = vec![
            "--buffer_size",
            "10",
            "--fork",
            "4",
            "--warning_file",
            warnings.to_str().unwrap(),
        ];
        args.extend(mode_flags);
        let cmd = vep_command(&input, &corpus.cache_dir(), &corpus.assembly, &out, &args);
        let output = output_within(cmd, Duration::from_secs(120));
        assert!(
            !output.status.success(),
            "{mode}: a warnings file that cannot be written must fail the run"
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        error_lines.push(
            stderr
                .lines()
                .find(|l| l.starts_with("ERROR: [output] Could not write to warnings file"))
                .unwrap_or_else(|| panic!("{mode}: no warnings-file ERROR line:\n{stderr}"))
                .to_string(),
        );
        assert_same(
            &reference,
            &without_run_headers(&fs::read(&out).unwrap()),
            &format!("{mode}: the output is not the batch whose skip met the failure"),
        );
    }
    assert_eq!(error_lines[0], error_lines[1]);
}

/// In both modes and under every writer strategy, an output that cannot be
/// written ends the run with an error naming the output rather than hanging:
/// the header lines, which the direct and vectored strategies put on the
/// output before the pipeline starts, fail as the buffered strategy's write or
/// final flush does. The warnings go to standard error, so the default warnings
/// file beside the unwritable output is not what fails first.
#[cfg(unix)]
#[test]
fn unwritable_output_fails_instead_of_hanging() {
    let corpus = &corpora()[0];
    for (mode, mode_flags) in MODES {
        for strategy in WRITER_STRATEGIES {
            let mut args = vec![
                "--buffer_size",
                "50",
                "--fork",
                "4",
                "--warning_file",
                "STDERR",
            ];
            args.extend(mode_flags);
            args.extend(["--writer_strategy", strategy]);
            let cmd = vep_command(
                &corpus.dir.join("variants.vcf"),
                &corpus.cache_dir(),
                &corpus.assembly,
                Path::new("/dev/full"),
                &args,
            );
            let output = output_within(cmd, Duration::from_secs(120));
            assert!(
                !output.status.success(),
                "{mode} {strategy}: writing to /dev/full must fail"
            );
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(
                stderr.contains("Failed to write output")
                    || stderr.contains("Failed to flush output"),
                "{mode} {strategy}: stderr does not name the output failure:\n{stderr}"
            );
        }
    }
}
