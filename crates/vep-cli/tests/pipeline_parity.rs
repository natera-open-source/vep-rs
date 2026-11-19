// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Pipeline parity: the reader, coordinator and writer stages must produce the
//! same bytes whatever the batch size, thread count, pipeline mode, pipeline
//! depth, task size, writer strategy, decompression thread count and
//! compression thread count, keep passthrough lines in input order, and fail
//! cleanly rather than hang when a stage errors.

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

/// The `vep` invocation with `--verbose` in place of `--quiet`, so the debug
/// lines the run's decisions write reach standard error.
fn vep_verbose(input: &Path, json_cache: &Path, assembly: &str, out: &Path, args: &[&str]) -> Run {
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("-i")
        .arg(input)
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(json_cache)
        .args(["--species", "homo_sapiens", "--assembly", assembly])
        .args(["--force_overwrite", "--no_stats", "--verbose"])
        .args(args);
    run_command(cmd, out)
}

/// The corpus the single-corpus cases run on: the first human corpus by name,
/// whose records sit on one chromosome in position order and number several
/// hundred, the shape the cases that index, truncate or splice the input need;
/// every command here passes `--species homo_sapiens`.
fn single_corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|corpus| corpus.species == common::DEFAULT_SPECIES)
        .expect("a human golden corpus")
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

/// Runs `vep` over the corpus's `variants.vcf` with the arguments a RefSeq or
/// merged corpus needs ahead of `args`: the cache declares its transcript set,
/// which the run must select, and reads the reference for its edited models.
fn run_corpus(corpus: &Corpus, out: &Path, args: &[&str]) -> Run {
    let set_args = corpus.transcript_set_args();
    let mut all: Vec<&str> = set_args.iter().map(String::as_str).collect();
    all.extend_from_slice(args);
    vep(
        &corpus.dir.join("variants.vcf"),
        &corpus.cache_dir(),
        &corpus.assembly,
        out,
        &all,
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

/// The corpus named `name` under `tests/golden`.
fn corpus_named(name: &str) -> Corpus {
    corpora()
        .into_iter()
        .find(|corpus| corpus.name == name)
        .unwrap_or_else(|| panic!("no corpus {name} under tests/golden"))
}

/// A run whose sources are read ahead of its annotation (the variation files,
/// a `--custom` VCF and bigWig, two plugins whose prefetch reads coordinates
/// alone, gnomADc and REVEL, and one whose prefetch reads the consequence rows,
/// ClinPred, whose slot the annotate phase fills) writes the bytes of its
/// batched `--fork 4 --buffer_size 5000` run in the default and VCF formats
/// whatever the mode, fork count, buffer size and pipeline depth, the pipeline
/// depth of 1 included, which holds the read-ahead to one batch.
#[test]
fn read_ahead_sources_write_the_same_bytes_at_every_fork_buffer_and_depth() {
    let tmp = tempfile::tempdir().unwrap();
    let colocated = corpus_named("GRCh38-colocated");
    let conservation = corpus_named("GRCh38-plugin-conservation");
    let gnomadc = corpus_named("GRCh38-plugin-gnomadc");
    let revel = corpus_named("GRCh38-plugin-revel");
    let clinpred = corpus_named("GRCh38-plugin-clinpred-fixture");
    // A `--custom` VCF of the conservation corpus's own records, each carrying
    // one INFO field, beside the corpus's bigWig.
    let custom_vcf = tmp.path().join("custom_sites.vcf.gz");
    bgzip_with_index(&conservation.dir.join("variants.vcf"), &custom_vcf, 40, 1);
    let custom_vcf = format!(
        "file={},short_name=sites,format=vcf,type=exact,fields=PAD",
        custom_vcf.display()
    );
    let bigwig = format!(
        "file={},short_name=gerp,format=bigwig,type=overlap",
        conservation
            .dir
            .join("plugin_data")
            .join("gerp.bw")
            .display()
    );
    let coverage = format!(
        "gnomADc,{}",
        gnomadc
            .dir
            .join("plugin_data")
            .join("gnomad.coverage.tsv.gz")
            .display()
    );
    let revel_file = format!(
        "REVEL,file={}",
        revel.dir.join("plugin_data").join("revel.tsv.gz").display()
    );
    let clinpred_file = format!(
        "ClinPred,file={}",
        clinpred
            .dir
            .join("data")
            .join("ClinPred_hg38_sorted_tabbed_synthetic.tsv.gz")
            .display()
    );
    let arms: [(&str, &Corpus, Vec<&str>); 5] = [
        (
            "variation",
            &colocated,
            vec!["--check_existing", "--af", "--pubmed"],
        ),
        (
            "custom",
            &conservation,
            vec!["--custom", &custom_vcf, "--custom", &bigwig],
        ),
        ("gnomadc", &gnomadc, vec!["--plugin", &coverage]),
        ("revel", &revel, vec!["--plugin", &revel_file]),
        ("clinpred", &clinpred, vec!["--plugin", &clinpred_file]),
    ];
    for (arm, corpus, flags) in &arms {
        for (name, format_flags) in [FORMATS[0], FORMATS[2]] {
            let reference_out = tmp.path().join(format!("ahead_{arm}_{name}_reference.out"));
            let mut args: Vec<&str> = flags.clone();
            args.extend(format_flags);
            args.extend([
                "--pipeline_mode",
                "batched",
                "--fork",
                "4",
                "--buffer_size",
                "5000",
            ]);
            let reference = run_corpus(corpus, &reference_out, &args);
            assert!(
                reference.success,
                "{arm}/{name} reference: {}",
                reference.stderr
            );
            let reference = without_run_headers(&reference.output);
            assert!(
                !reference.is_empty(),
                "{arm}/{name} reference: empty output"
            );
            for (mode, mode_flags) in MODES {
                for fork in ["1", "4"] {
                    for buffer in ["5000", "7"] {
                        for depth in ["1", "3"] {
                            let out = tmp.path().join(format!(
                                "ahead_{arm}_{name}_{mode}_{fork}_{buffer}_{depth}.out"
                            ));
                            let mut args: Vec<&str> = flags.clone();
                            args.extend(format_flags);
                            args.extend(mode_flags);
                            args.extend(["--fork", fork, "--buffer_size", buffer]);
                            args.extend(["--pipeline_depth", depth]);
                            let run = run_corpus(corpus, &out, &args);
                            assert!(
                                run.success,
                                "{arm}/{name} {mode} fork {fork} buffer {buffer} depth {depth}: {}",
                                run.stderr
                            );
                            assert_same(
                                &reference,
                                &without_run_headers(&run.output),
                                &format!(
                                    "{arm}/{name}: {mode} --fork {fork} --buffer_size {buffer} --pipeline_depth {depth} differs from the batched --fork 4 --buffer_size 5000 output"
                                ),
                            );
                        }
                    }
                }
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
/// (`LoFTEE` under the percentile rule loads with no data file and adds only
/// its own `LoF` keys). With
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
                "LoFTEE,use_gerp_end_trunc=0",
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
    let corpus = &single_corpus();
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
    let corpus = &single_corpus();
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

/// `--compress_output bgzip` under every form of `--compression_threads` (the
/// single-threaded encoder, worker pools of two and three, and the `0` rule),
/// in both modes, at fork 1 and 4, is a BGZF stream that decompresses to the
/// bytes of the plain default run in the default and VCF formats. The
/// compressed streams are not compared across runs: the header carries the
/// run's command line, so they differ by the flag's own value.
#[test]
fn bgzip_output_decompresses_to_the_plain_bytes_at_every_compression_thread_count() {
    use std::io::Read;

    let tmp = tempfile::tempdir().unwrap();
    let corpus = &single_corpus();
    for (name, flags) in [FORMATS[0], FORMATS[2]] {
        let reference_out = tmp.path().join(format!("plain_threads_{name}.out"));
        let mut args: Vec<&str> = flags.to_vec();
        args.extend(["--buffer_size", "100", "--fork", "4"]);
        let reference = run_corpus(corpus, &reference_out, &args);
        assert!(reference.success, "{name} reference: {}", reference.stderr);
        let reference = without_run_headers(&reference.output);
        assert!(!reference.is_empty());
        for (mode, mode_flags) in MODES {
            for fork in ["1", "4"] {
                for threads in ["1", "2", "3", "0"] {
                    let out = tmp
                        .path()
                        .join(format!("bgzip_{name}_{mode}_{fork}_{threads}.out.gz"));
                    let mut args: Vec<&str> = flags.to_vec();
                    args.extend(["--buffer_size", "100", "--fork", fork]);
                    args.extend(mode_flags);
                    args.extend(["--compress_output", "bgzip"]);
                    args.extend(["--compression_threads", threads]);
                    let run = run_corpus(corpus, &out, &args);
                    let case = format!(
                        "{name}: {mode} --fork {fork} --compress_output bgzip --compression_threads {threads}"
                    );
                    assert!(run.success, "{case}: {}", run.stderr);
                    assert_eq!(&run.output[..4], &[0x1f, 0x8b, 0x08, 0x04], "{case}");
                    let mut plain = Vec::new();
                    noodles_bgzf::io::Reader::new(&run.output[..])
                        .read_to_end(&mut plain)
                        .unwrap_or_else(|e| panic!("{case}: not a BGZF stream: {e}"));
                    assert_same(
                        &reference,
                        &without_run_headers(&plain),
                        &format!("{case} differs from the plain run once decompressed"),
                    );
                }
            }
        }
    }
}

/// A Parquet run, in both modes and under every writer strategy, writes the
/// same rows in the same order, read back through the cell contract.
#[cfg(feature = "parquet")]
#[test]
fn parquet_output_is_identical_in_both_modes_and_strategies() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &single_corpus();
    let mut reference: Option<Vec<String>> = None;
    for (mode, mode_flags) in MODES {
        for strategy in WRITER_STRATEGIES {
            let out = tmp
                .path()
                .join(format!("parquet_{mode}_{strategy}.parquet"));
            let mut args = vec![
                "--output_format",
                "parquet",
                "--buffer_size",
                "100",
                "--fork",
                "4",
            ];
            args.extend(mode_flags);
            args.extend(["--writer_strategy", strategy]);
            let run = run_corpus(corpus, &out, &args);
            assert!(run.success, "{mode} {strategy}: {}", run.stderr);
            let rows = vep_io::output::parquet::tab_rows(&out).unwrap();
            assert!(!rows.is_empty(), "{mode} {strategy}: no Parquet rows");
            match &reference {
                None => reference = Some(rows),
                Some(expected) => assert_same(
                    expected.join("\n").as_bytes(),
                    rows.join("\n").as_bytes(),
                    &format!("{mode} {strategy}: the Parquet rows differ from batched buffered"),
                ),
            }
        }
    }
}

/// Every pipeline depth and task size, in both modes, writes the bytes of the
/// default run.
#[test]
fn pipeline_depth_and_task_size_leave_the_output_unchanged() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &single_corpus();
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
/// bytes under every form of `--decompression_threads`, the flag absent
/// included, and under `--fork 0`, in both modes, whether its records inflate
/// to a few hundred bytes or past the width that asks for the inflate-bound
/// pool; the probe arms assert which side of the threshold each copy's shape
/// falls on, the count `auto` resolves there, the count the absent flag
/// resolves (the architecture's default) and the annotator count `--fork 0`
/// resolves beside those inflaters.
#[test]
fn bgzipped_input_matches_the_plain_file_under_every_decompression_thread_count() {
    use vep_cli::config::{
        reported_parallelism, resolve_fork_auto, DecompressionThreads, FORK_AUTO_CAP,
    };
    use vep_cli::input_shape::InputShape;

    let tmp = tempfile::tempdir().unwrap();
    let corpus = &single_corpus();
    let cpus = vep_cli::config::host_parallelism();
    let bound = vep_cli::config::INFLATE_BOUND_BYTES_PER_RECORD;
    let host_capped = cpus.min(vep_cli::config::DECOMPRESSION_THREADS_AUTO_CAP);
    let inflate_bound = cpus.min(vep_cli::config::INFLATE_BOUND_INFLATERS);
    let value_of = |resolved: &str, name: &str| -> u64 {
        resolved
            .lines()
            .find_map(|line| line.strip_prefix(&format!("{name}\t")))
            .and_then(|value| value.parse().ok())
            .unwrap_or_else(|| panic!("no {name} count in {resolved}"))
    };
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
        let inflated = value_of(&resolved, "inflated_bytes_per_record");
        assert!(
            inflated_bytes.contains(&inflated),
            "{width}: {inflated} inflated bytes per record"
        );
        // The absent flag is the architecture's default, and `--fork 0` beside
        // it resolves the annotators from the same shape and the inflaters the
        // default leaves it.
        let default_threads = match DecompressionThreads::default() {
            DecompressionThreads::Auto => auto_threads,
            DecompressionThreads::HostCapped => host_capped,
            DecompressionThreads::Fixed(count) => count,
        };
        let probe = vep(
            &gz,
            &corpus.cache_dir(),
            &corpus.assembly,
            Path::new("STDOUT"),
            &["--resolved_config", "--fork", "0"],
        );
        assert!(probe.success, "{width} default probe: {}", probe.stderr);
        let resolved = String::from_utf8_lossy(&probe.output);
        assert_eq!(
            value_of(&resolved, "decompression_threads"),
            default_threads as u64,
            "{width}: the absent flag's count"
        );
        let shape = InputShape::Indexed {
            inflated_bytes_per_record: inflated,
            records_estimate: value_of(&resolved, "records_estimate"),
            file_bytes: fs::metadata(&gz).unwrap().len(),
        };
        let inflaters = if default_threads > 1 {
            default_threads
        } else {
            0
        };
        let fork = value_of(&resolved, "fork") as usize;
        let host_cpus = reported_parallelism().unwrap_or(1);
        let expected = resolve_fork_auto(&shape, inflaters, host_cpus);
        assert_eq!(
            fork, expected,
            "{width}: fork 0 beside {inflaters} inflaters"
        );
        assert!((1..=FORK_AUTO_CAP).contains(&fork), "{width}: fork {fork}");
        let host_rule = vep_cli::config::resolve_fork(0);
        if width == "narrow" {
            assert_eq!(fork, host_rule, "{width}: the host's rule");
        }
        // An inflater pool as wide as the host leaves the annotators the half
        // floor on the wide copy, below the host's rule on any host of two or
        // more CPUs, so the shape and the inflaters reach the rule whatever the
        // host's width; the narrow copy keeps the host's rule.
        let every_cpu = host_cpus.to_string();
        let probe = vep(
            &gz,
            &corpus.cache_dir(),
            &corpus.assembly,
            Path::new("STDOUT"),
            &[
                "--resolved_config",
                "--fork",
                "0",
                "--decompression_threads",
                &every_cpu,
            ],
        );
        assert!(probe.success, "{width} every-CPU probe: {}", probe.stderr);
        let resolved = String::from_utf8_lossy(&probe.output);
        let fork_every_cpu = value_of(&resolved, "fork") as usize;
        if width == "wide" && host_cpus >= 2 {
            assert_eq!(
                fork_every_cpu,
                (host_rule / 2).max(1),
                "{width}: the half floor"
            );
            assert!(
                fork_every_cpu < host_rule,
                "{width}: {fork_every_cpu} below the host's {host_rule}"
            );
        } else {
            assert_eq!(fork_every_cpu, host_rule, "{width}: the host's rule");
        }
        // The run's own pool, sized where the input opens, says what it
        // resolved; the inflaters are every CPU when there are two or more.
        let inflaters_every_cpu = if host_cpus > 1 { host_cpus } else { 0 };
        let pool_line = format!(
            "--fork 0 resolved {fork_every_cpu} annotators beside {inflaters_every_cpu} inflaters from an input of indexed shape"
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
                // (label, the decompression flag or none, the fork value)
                for (label, threads, fork) in [
                    ("auto", Some("auto"), "4"),
                    ("1", Some("1"), "4"),
                    ("0", Some("0"), "4"),
                    ("default", None, "4"),
                    ("fork-0", None, "0"),
                    ("fork-0-every-cpu", Some(every_cpu.as_str()), "0"),
                ] {
                    let out = tmp
                        .path()
                        .join(format!("gz_{width}_{name}_{mode}_{label}.out"));
                    let mut args: Vec<&str> = flags.to_vec();
                    args.extend(["--buffer_size", "100", "--fork", fork]);
                    args.extend(mode_flags);
                    if let Some(threads) = threads {
                        args.extend(["--decompression_threads", threads]);
                    }
                    let every_cpu_arm = label == "fork-0-every-cpu";
                    let run = if every_cpu_arm {
                        vep_verbose(&gz, &corpus.cache_dir(), &corpus.assembly, &out, &args)
                    } else {
                        vep(&gz, &corpus.cache_dir(), &corpus.assembly, &out, &args)
                    };
                    assert!(run.success, "{width}/{name} {mode} {label}: {}", run.stderr);
                    if every_cpu_arm {
                        assert!(
                            run.stderr.contains(&pool_line),
                            "{width}/{name} {mode}: expected {pool_line:?} in {}",
                            run.stderr
                        );
                    }
                    assert_same(
                        &reference,
                        &without_run_headers(&run.output),
                        &format!(
                            "{width}/{name}: {mode} {label} on the bgzipped copy differs from the plain file"
                        ),
                    );
                }
            }
        }
    }
}

/// A bgzipped, indexed copy of the corpus that spans every chromosome of its
/// cache annotates to the plain file's bytes under every `--prewarm_ahead`
/// value (none ahead, one, two, every indexed contig at the start) and under a
/// pinned `--parse_threads`, in both modes: the order the chromosomes' shards
/// are parsed in, and the pool that parses them, change no output.
#[test]
fn prewarm_ahead_and_parse_threads_write_the_plain_bytes_across_chromosomes() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpus_named("GRCh38");
    let gz = tmp.path().join("chromosomes.vcf.gz");
    let plain = bgzip_with_index(&corpus.dir.join("variants.vcf"), &gz, 20, 0);
    let order = vep(
        &gz,
        &corpus.cache_dir(),
        &corpus.assembly,
        Path::new("STDOUT"),
        &["--resolved_config"],
    );
    assert!(order.success, "{}", order.stderr);
    assert!(
        String::from_utf8_lossy(&order.output).contains("\ninput_shape\tindexed\n"),
        "{}",
        String::from_utf8_lossy(&order.output)
    );
    for (name, flags) in [FORMATS[0], FORMATS[2]] {
        let reference_out = tmp.path().join(format!("chromosomes_{name}_reference.out"));
        let mut args: Vec<&str> = flags.to_vec();
        args.extend(["--buffer_size", "50", "--fork", "4"]);
        let reference = vep(
            &plain,
            &corpus.cache_dir(),
            &corpus.assembly,
            &reference_out,
            &args,
        );
        assert!(reference.success, "{name} reference: {}", reference.stderr);
        let reference = without_run_headers(&reference.output);
        assert!(!reference.is_empty());
        for (mode, mode_flags) in MODES {
            for (ahead, parse_threads) in [
                ("0", "0"),
                ("1", "0"),
                ("2", "0"),
                ("all", "0"),
                ("1", "1"),
                ("all", "2"),
            ] {
                let out = tmp.path().join(format!(
                    "chromosomes_{name}_{mode}_{ahead}_{parse_threads}.out"
                ));
                let mut args: Vec<&str> = flags.to_vec();
                args.extend(["--buffer_size", "50", "--fork", "4"]);
                args.extend(mode_flags);
                args.extend(["--prewarm_ahead", ahead, "--parse_threads", parse_threads]);
                let run = vep(&gz, &corpus.cache_dir(), &corpus.assembly, &out, &args);
                assert!(
                    run.success,
                    "{name} {mode} --prewarm_ahead {ahead} --parse_threads {parse_threads}: {}",
                    run.stderr
                );
                assert_same(
                    &reference,
                    &without_run_headers(&run.output),
                    &format!(
                        "{name}: {mode} --prewarm_ahead {ahead} --parse_threads {parse_threads} differs from the plain file"
                    ),
                );
            }
        }
    }
}

/// Every form of `--reader_core` writes the plain file's bytes on the wide
/// bgzipped copy (the inflate-bound shape, where `auto` may reserve a core) and
/// on the plain file, in both modes; under `--verbose` the run says what it
/// decided, and a fixed CPU says it reserved its core exactly when the rule
/// reserves one on this host (the CPU and its core's other CPUs allowed, with a
/// CPU left over), else that no core was reserved.
#[test]
fn reader_core_forms_write_the_same_bytes_and_say_what_they_reserved() {
    use vep_cli::config::{reserve_reader_core, ReaderCore};
    use vep_cli::input_shape::CpuTopology;

    let tmp = tempfile::tempdir().unwrap();
    let corpus = &single_corpus();
    let gz = tmp.path().join("wide.vcf.gz");
    let plain = bgzip_with_index(&corpus.dir.join("variants.vcf"), &gz, 40, 8_192);
    let allowed = vep_cli::affinity::allowed_cpus().unwrap_or_default();
    let topology = CpuTopology::read();
    let fixed_cpu = allowed.first().copied();
    let fixed = fixed_cpu.map(|cpu| cpu.to_string());
    let mut forms = vec!["none", "auto"];
    if let Some(cpu) = fixed.as_deref() {
        forms.push(cpu);
    }
    for (name, flags) in [FORMATS[0], FORMATS[2]] {
        let reference_out = tmp.path().join(format!("core_{name}_reference.out"));
        let mut args: Vec<&str> = flags.to_vec();
        args.extend(["--buffer_size", "100", "--fork", "4"]);
        let reference = vep(
            &plain,
            &corpus.cache_dir(),
            &corpus.assembly,
            &reference_out,
            &args,
        );
        assert!(reference.success, "{name} reference: {}", reference.stderr);
        let reference = without_run_headers(&reference.output);
        for (input, label) in [(&gz, "wide"), (&plain, "plain")] {
            for (mode, mode_flags) in MODES {
                for form in &forms {
                    let out = tmp
                        .path()
                        .join(format!("core_{name}_{label}_{mode}_{form}.out"));
                    let mut args: Vec<&str> = flags.to_vec();
                    args.extend(["--buffer_size", "100", "--fork", "4"]);
                    args.extend(mode_flags);
                    args.extend(["--reader_core", form]);
                    let run =
                        vep_verbose(input, &corpus.cache_dir(), &corpus.assembly, &out, &args);
                    assert!(
                        run.success,
                        "{name}/{label} {mode} --reader_core {form}: {}",
                        run.stderr
                    );
                    assert_same(
                        &reference,
                        &without_run_headers(&run.output),
                        &format!(
                            "{name}/{label}: {mode} --reader_core {form} differs from the plain default"
                        ),
                    );
                    let decided = format!("--reader_core {form}: ");
                    match (*form, fixed_cpu) {
                        ("none", _) => assert!(
                            !run.stderr.contains(&decided),
                            "{name}/{label} {mode}: none decided something: {}",
                            run.stderr
                        ),
                        ("auto", _) => {}
                        (_, Some(cpu)) => {
                            let reserved = reserve_reader_core(
                                ReaderCore::Fixed(cpu),
                                topology.as_ref(),
                                &allowed,
                                4,
                                0,
                                0,
                            );
                            let expected = match reserved {
                                Some(reservation) => {
                                    format!("{decided}reader on CPU {}", reservation.reader_cpu)
                                }
                                None => format!("{decided}no core reserved"),
                            };
                            assert!(
                                run.stderr.contains(&expected),
                                "{name}/{label} {mode}: expected {expected:?} in {}",
                                run.stderr
                            );
                        }
                        (_, None) => unreachable!("a fixed form needs an allowed CPU"),
                    }
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
    let corpus = &single_corpus();
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
    let corpus = &single_corpus();
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

/// The GRCh38 corpus built with `--refseq` alone: its stored warnings log names
/// two records that draw `Transcript-assembly mismatch`.
fn refseq_corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|corpus| corpus.assembly == "GRCh38" && corpus.transcript_set_flags() == ["--refseq"])
        .expect("a GRCh38 RefSeq golden corpus")
}

/// The records a warnings file's `Transcript-assembly mismatch in
/// <Uploaded_variation>` lines name, in file order. A forked Ensembl VEP run
/// opens a child's lines with its process id (`WARNING: <pid> : WARNING: ...`,
/// Runner.pm:565), which is not part of the record's name.
fn mismatch_warnings(warnings: &str) -> Vec<String> {
    warnings
        .lines()
        .filter_map(|line| {
            line.rsplit_once("Transcript-assembly mismatch in ")
                .map(|(_, record)| record.trim().to_string())
        })
        .collect()
}

/// Each line of a warnings file as a label: `skip N` for a skipped record's
/// line, the record's name for a `Transcript-assembly mismatch` line; any other
/// line is unexpected here.
fn warning_labels(warnings: &str) -> Vec<String> {
    warnings
        .lines()
        .map(|line| {
            if let Some(rest) = line.strip_prefix("WARNING: line ") {
                let number = rest.split(' ').next().unwrap_or_default();
                return format!("skip {number}");
            }
            line.rsplit_once("Transcript-assembly mismatch in ")
                .map(|(_, record)| record.trim().to_string())
                .unwrap_or_else(|| panic!("unexpected warning line: {line}"))
        })
        .collect()
}

/// In both modes the warnings file lists the warnings a record draws while it is
/// annotated in input order, and a streaming run lists them the same way every
/// time: Ensembl VEP writes `Transcript-assembly mismatch in
/// <Uploaded_variation>` as it renders each record of the buffer, in buffer
/// order (OutputFactory.pm:300-303, :1718), and the corpus's stored log holds
/// the two records that draw it in that order. At `--fork 4 --buffer_size 1`
/// every record is a batch of its own, so a warning written as it was drawn
/// would land in task-completion order.
#[test]
fn record_warnings_are_listed_in_input_order_in_both_modes() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = refseq_corpus();
    let reference =
        fs::read_to_string(corpus.dir.join("expected").join("default.warnings.log")).unwrap();
    let expected = mismatch_warnings(&reference);
    assert!(expected.len() >= 2, "{}: {reference}", corpus.name);

    let mut batched: Option<Vec<u8>> = None;
    for (mode, mode_flags) in MODES {
        let repeats = if mode == "streaming" { 12 } else { 1 };
        for repeat in 0..repeats {
            let name = format!("{mode} run {repeat}");
            let out = tmp.path().join(format!("corpus_{mode}_{repeat}.out"));
            let mut args = vec!["--fork", "4", "--buffer_size", "1"];
            args.extend(mode_flags);
            let run = run_corpus(&corpus, &out, &args);
            assert!(run.success, "{name}: {}", run.stderr);
            let warnings = fs::read(format!("{}_warnings.txt", out.display())).unwrap();
            assert_eq!(
                mismatch_warnings(&String::from_utf8_lossy(&warnings)),
                expected,
                "{name}"
            );
            match &batched {
                None => batched = Some(warnings),
                Some(reference) => assert_same(
                    reference,
                    &warnings,
                    &format!("{name}: warnings file differs from batched"),
                ),
            }
        }
    }
}

/// A record's warning keeps its place among the skips in both modes. Ensembl
/// VEP's parser warns about a skipped record while the buffer fills
/// (Parser.pm:189, Runner.pm:423) and the renderer warns about a record once
/// the buffer is annotated (Runner.pm:443, OutputFactory.pm:1718): at
/// `--buffer_size 1` every record is a buffer, so the two kinds alternate in
/// input order, and at a buffer size past the input the batch's skip comes
/// first and the records' warnings follow in record order. Forty copies of the
/// corpus's two records that draw the warning, each under a name of its own,
/// the thirtieth line with a POS that is not a number: the one batch is cut
/// into five tasks of eight (`--render_task_records 8`), so the skip sits in
/// the fourth, after the warnings of three tasks that a writer reporting them
/// as each chunk comes due would list first.
#[test]
fn record_warnings_keep_their_place_among_the_skips_in_both_modes() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = refseq_corpus();
    let reference =
        fs::read_to_string(corpus.dir.join("expected").join("default.warnings.log")).unwrap();
    let drawing = mismatch_warnings(&reference);
    assert!(drawing.len() >= 2, "{}: {reference}", corpus.name);
    let text = fs::read_to_string(corpus.dir.join("variants.vcf")).unwrap();
    let record_line = |id: &str| -> Vec<String> {
        text.lines()
            .find(|line| !line.starts_with('#') && line.split('\t').nth(2) == Some(id))
            .unwrap_or_else(|| panic!("{}: no record {id}", corpus.name))
            .split('\t')
            .map(str::to_string)
            .collect()
    };
    let (first, second) = (record_line(&drawing[0]), record_line(&drawing[1]));
    let mut lines: Vec<String> = text
        .lines()
        .filter(|l| l.starts_with('#'))
        .map(str::to_string)
        .collect();
    let mut in_input_order: Vec<String> = Vec::new();
    for record in 1..=40usize {
        let mut fields = if record <= 20 {
            first.clone()
        } else {
            second.clone()
        };
        fields[2] = format!("m{record:02}");
        if record == 30 {
            fields[1] = "abc".to_string();
            in_input_order.push("skip 30".to_string());
        } else {
            in_input_order.push(fields[2].clone());
        }
        lines.push(fields.join("\t"));
    }
    let input = tmp.path().join("mismatches.vcf");
    fs::write(&input, lines.join("\n") + "\n").unwrap();
    let mut skip_first = vec!["skip 30".to_string()];
    skip_first.extend(in_input_order.iter().filter(|l| *l != "skip 30").cloned());
    let set_args = corpus.transcript_set_args();

    for (buffer_size, expected) in [("1", &in_input_order), ("5000", &skip_first)] {
        let mut batched: Option<Vec<u8>> = None;
        for (mode, mode_flags) in MODES {
            let repeats = if mode == "streaming" { 8 } else { 1 };
            for repeat in 0..repeats {
                let name = format!("{mode} run {repeat} at --buffer_size {buffer_size}");
                let out = tmp
                    .path()
                    .join(format!("mismatches_{buffer_size}_{mode}_{repeat}.out"));
                let mut args: Vec<&str> = set_args.iter().map(String::as_str).collect();
                args.extend(["--fork", "4", "--render_task_records", "8"]);
                args.extend(["--buffer_size", buffer_size]);
                args.extend(mode_flags);
                let run = vep(&input, &corpus.cache_dir(), &corpus.assembly, &out, &args);
                assert!(run.success, "{name}: {}", run.stderr);
                let warnings = fs::read(format!("{}_warnings.txt", out.display())).unwrap();
                assert_eq!(
                    &warning_labels(&String::from_utf8_lossy(&warnings)),
                    expected,
                    "{name}"
                );
                match &batched {
                    None => batched = Some(warnings),
                    Some(reference) => assert_same(
                        reference,
                        &warnings,
                        &format!("{name}: warnings file differs from batched"),
                    ),
                }
            }
        }
    }
}

/// In both modes, a warnings file that cannot be written ends the run at the
/// batch whose record drew the first warning, with that batch written and
/// nothing after it, as it does for a skip: the RefSeq corpus's first record to
/// draw `Transcript-assembly mismatch` is the last one written at
/// `--buffer_size 1`, so a clean run over the records through it is the
/// reference for both.
#[test]
fn unwritable_warnings_file_ends_both_modes_at_the_record_that_warned() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = refseq_corpus();
    let reference_log =
        fs::read_to_string(corpus.dir.join("expected").join("default.warnings.log")).unwrap();
    let first_warned = mismatch_warnings(&reference_log)
        .into_iter()
        .next()
        .unwrap_or_else(|| panic!("{}: {reference_log}", corpus.name));
    let text = fs::read_to_string(corpus.dir.join("variants.vcf")).unwrap();
    let header: Vec<&str> = text.lines().filter(|l| l.starts_with('#')).collect();
    let data: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    let through = data
        .iter()
        .position(|line| line.split('\t').nth(2) == Some(first_warned.as_str()))
        .unwrap_or_else(|| panic!("{}: no record {first_warned}", corpus.name));
    let mut lines: Vec<String> = header.iter().map(|l| l.to_string()).collect();
    lines.extend(data[..=through].iter().map(|l| l.to_string()));
    let clean = tmp.path().join("clean.vcf");
    fs::write(&clean, lines.join("\n") + "\n").unwrap();
    let warnings = tmp.path().join("missing").join("w.txt");
    let set_args = corpus.transcript_set_args();
    let set_args: Vec<&str> = set_args.iter().map(String::as_str).collect();

    let reference_out = tmp.path().join("clean.out");
    let mut args = set_args.clone();
    args.extend(["--buffer_size", "1", "--fork", "4"]);
    let reference = vep(
        &clean,
        &corpus.cache_dir(),
        &corpus.assembly,
        &reference_out,
        &args,
    );
    assert!(reference.success, "clean run: {}", reference.stderr);
    let reference = without_run_headers(&reference.output);
    assert!(!reference.is_empty());

    let mut error_lines: Vec<String> = Vec::new();
    for (mode, mode_flags) in MODES {
        let out = tmp.path().join(format!("warned_{mode}.out"));
        let mut args = set_args.clone();
        args.extend([
            "--buffer_size",
            "1",
            "--fork",
            "4",
            "--warning_file",
            warnings.to_str().unwrap(),
        ]);
        args.extend(mode_flags);
        let cmd = vep_command(
            &corpus.dir.join("variants.vcf"),
            &corpus.cache_dir(),
            &corpus.assembly,
            &out,
            &args,
        );
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
            &format!("{mode}: the output is not the records through the one that warned"),
        );
    }
    assert_eq!(error_lines[0], error_lines[1]);
}

/// In both modes, a warnings file that cannot be written ends the run with the
/// same error line and the same output: the batch whose skip met the failure is
/// written and nothing after it. The skipped record (a chromosome the cache
/// lacks) opens the first batch of 10, so a clean run over the nine records
/// beside it is the reference for both.
#[test]
fn unwritable_warnings_file_ends_both_modes_at_the_same_batch() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &single_corpus();
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
/// file beside the unwritable output is not what fails first. The output is
/// Linux's `/dev/full`, which opens and then fails every write; macOS has no
/// such device.
#[cfg(target_os = "linux")]
#[test]
fn unwritable_output_fails_instead_of_hanging() {
    let corpus = &single_corpus();
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

/// The statistics report without the rows that name the run rather than its
/// counts: the command line (which carries the knobs under test), the three
/// times and the output path.
fn without_run_rows(summary: &str) -> String {
    summary
        .lines()
        .filter(|line| {
            ![
                "Command line options\t",
                "Start time\t",
                "End time\t",
                "Run time\t",
                "Output file\t",
            ]
            .iter()
            .any(|prefix| line.starts_with(prefix))
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// Runs `vep --stats_text` on `corpus` with `args` and returns the report.
fn summary_of(corpus: &Corpus, out: &Path, args: &[&str]) -> String {
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("-i")
        .arg(corpus.dir.join("variants.vcf"))
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(corpus.cache_dir())
        .args(["--species", "homo_sapiens", "--assembly", &corpus.assembly])
        .args(["--force_overwrite", "--stats_text", "--quiet"])
        .args(corpus.transcript_set_args())
        .args(args);
    let result = cmd.output().expect("vep binary runs");
    assert!(
        result.status.success(),
        "{} {args:?}: {}",
        corpus.name,
        String::from_utf8_lossy(&result.stderr)
    );
    let summary = out.with_file_name(format!(
        "{}_summary.txt",
        out.file_name().unwrap().to_string_lossy()
    ));
    fs::read_to_string(summary).unwrap()
}

/// The statistics report is the same whatever the thread count, buffer size,
/// pipeline mode, task size and output format: the chunks' counts merge into one
/// run, and the chart rows are ordered by rule, never by the order the chunks
/// finished in.
#[test]
fn summary_text_is_identical_across_fork_counts() {
    let tmp = tempfile::tempdir().unwrap();
    for corpus in &corpora() {
        let out = tmp.path().join(format!("{}_reference.out", corpus.name));
        let reference = without_run_rows(&summary_of(
            corpus,
            &out,
            &[
                "--fork",
                "1",
                "--buffer_size",
                "5000",
                "--pipeline_mode",
                "batched",
            ],
        ));
        let variants: [&[&str]; 3] = [
            &[
                "--fork",
                "4",
                "--buffer_size",
                "7",
                "--pipeline_mode",
                "streaming",
            ],
            &[
                "--fork",
                "4",
                "--buffer_size",
                "500",
                "--pipeline_mode",
                "batched",
            ],
            &[
                "--fork",
                "3",
                "--buffer_size",
                "64",
                "--render_task_records",
                "1",
            ],
        ];
        for (i, args) in variants.iter().enumerate() {
            let out = tmp.path().join(format!("{}_{i}.out", corpus.name));
            let candidate = without_run_rows(&summary_of(corpus, &out, args));
            assert_eq!(
                candidate, reference,
                "{}: {args:?} differs from fork 1 batched",
                corpus.name
            );
        }
    }
    let corpus = corpora()
        .into_iter()
        .find(|corpus| corpus.name == "GRCh38-statistics")
        .expect("the statistics corpus");
    let out = tmp.path().join("statistics_default.out");
    let reference = without_run_rows(&summary_of(&corpus, &out, &["--fork", "1"]));
    for (name, flags) in FORMATS {
        for (mode, mode_flags) in MODES {
            for fork in ["1", "4"] {
                let out = tmp
                    .path()
                    .join(format!("statistics_{name}_{mode}_{fork}.out"));
                let mut args: Vec<&str> = flags.to_vec();
                args.extend(["--fork", fork, "--buffer_size", "50"]);
                args.extend(mode_flags);
                let candidate = without_run_rows(&summary_of(&corpus, &out, &args));
                assert_eq!(
                    candidate, reference,
                    "GRCh38-statistics {name} {mode} --fork {fork} differs from the default run"
                );
            }
        }
    }
}
