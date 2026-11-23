// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--resolved_config` prints the runtime knobs (the thread counts resolved from
//! the input's shape as the run resolves them, the parse-pool, pipeline, writer,
//! reader-core and prewarm knobs as configured) and the input's shape, one
//! `<name>\t<value>` line each, then exits 0 without opening the output or
//! annotating; a named input file that is not there is refused as a run refuses
//! it. `Config::default()` holds what a command line without a flag resolves to,
//! field by field, but for the fields that resolution reads from the host.

mod common;

use std::path::Path;

use vep_cli::config::{
    host_parallelism, reported_parallelism, resolve_fork, DecompressionThreads,
    COMPRESSION_THREADS_AUTO_CAP, DECOMPRESSION_THREADS_AUTO_CAP,
};
use vep_cli::input_shape::physical_cores;
use vep_cli::Config;

use common::{corpora, isolated_vep_anywhere, Corpus};

/// Runs `vep --resolved_config` on `input` against `json_cache`, with `corpus`'s
/// assembly, then `args`: the exit success, standard output and standard error.
fn resolved_config(
    corpus: &Corpus,
    input: &Path,
    json_cache: &Path,
    out: &Path,
    args: &[&str],
) -> (bool, String, String) {
    let output = isolated_vep_anywhere()
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(out)
        .arg("--offline")
        .arg("--json_cache")
        .arg(json_cache)
        .args(["--species", "homo_sapiens", "--assembly", &corpus.assembly])
        .args(["--quiet", "--resolved_config"])
        .args(args)
        .output()
        .expect("vep binary runs");
    (
        output.status.success(),
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// The `available_parallelism` value the probe prints on this host: the count,
/// or `-` where the host reports none.
fn reported_parallelism_line() -> String {
    reported_parallelism().map_or("-".to_string(), |n| n.to_string())
}

/// The `physical_cores` value the probe prints on this host: the count, or `-`
/// where the host has no topology files.
fn physical_cores_line() -> String {
    physical_cores().map_or("-".to_string(), |n| n.to_string())
}

/// The compressed-input corpus: its BGZF twin `variants.vcf.gz` has a tabix
/// index beside it, so it is the one committed input with a shape.
fn indexed_corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-input-compressed")
        .expect("the GRCh38-input-compressed corpus is committed")
}

/// Every knob given on the command line prints as given (the numeric ones
/// resolved), a plain VCF has no shape, and nothing is written to the output
/// path.
#[test]
fn resolved_config_prints_the_knobs_and_exits_without_annotating() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    let out = tmp.path().join("never_written.out");
    let (success, stdout, stderr) = resolved_config(
        corpus,
        &corpus.dir.join("variants.vcf"),
        &corpus.cache_dir(),
        &out,
        &[
            "--fork",
            "3",
            "--buffer_size",
            "77",
            "--decompression_threads",
            "1",
            "--compression_threads",
            "2",
            "--pipeline_mode",
            "batched",
            "--pipeline_depth",
            "4",
            "--render_task_records",
            "9",
            "--writer_strategy",
            "vectored",
            "--output_buffer_size",
            "1234",
            "--reader_core",
            "none",
            "--parse_threads",
            "5",
            "--prewarm_ahead",
            "all",
        ],
    );
    assert!(success, "{stderr}");
    assert_eq!(
        stdout,
        format!(
            "fork\t3\n\
             buffer_size\t77\n\
             decompression_threads\t1\n\
             compression_threads\t2\n\
             parse_threads\t5\n\
             pipeline_mode\tbatched\n\
             pipeline_depth\t4\n\
             render_task_records\t9\n\
             writer_strategy\tvectored\n\
             output_buffer_size\t1234\n\
             reader_core\tnone\n\
             prewarm_ahead\tall\n\
             input_shape\tunknown\n\
             inflated_bytes_per_record\t-\n\
             records_estimate\t-\n\
             available_parallelism\t{}\n\
             physical_cores\t{}\n",
            reported_parallelism_line(),
            physical_cores_line()
        )
    );
    assert!(!out.exists(), "--resolved_config wrote the output file");
}

/// Every knob the command line leaves alone prints its built-in value, and the
/// committed bgzipped, indexed input prints its shape with the thread counts
/// the defaults resolve from it: the corpus's records inflate to 55 bytes on
/// average, under the inflation threshold, so the absent `--decompression_threads`
/// (`auto` on an aarch64 host, `0` elsewhere) is the `0` rule, the host's
/// parallelism capped at `DECOMPRESSION_THREADS_AUTO_CAP`, and the absent
/// `--fork` is the host's rule. The one block holds every data line whole, so
/// the bytes per record are the data lines' mean exactly; the record estimate
/// is the data line count plus the end-of-file block's share of the file.
#[test]
fn resolved_config_prints_the_defaults_and_an_indexed_input_shape() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = indexed_corpus();
    let out = tmp.path().join("never_written.out");
    let (success, stdout, stderr) = resolved_config(
        &corpus,
        &corpus.dir.join("variants.vcf.gz"),
        &corpus.cache_dir(),
        &out,
        &[],
    );
    assert!(success, "{stderr}");
    assert_eq!(
        DecompressionThreads::default(),
        if cfg!(target_arch = "aarch64") {
            DecompressionThreads::Auto
        } else {
            DecompressionThreads::HostCapped
        }
    );
    let lines: Vec<(&str, &str)> = stdout
        .lines()
        .map(|l| {
            l.split_once('\t')
                .unwrap_or_else(|| panic!("no tab in {l:?}"))
        })
        .collect();
    let cpus = host_parallelism();
    let fork = resolve_fork(0).to_string();
    let threads = cpus.min(DECOMPRESSION_THREADS_AUTO_CAP).to_string();
    let names: Vec<&str> = lines.iter().map(|(name, _)| *name).collect();
    assert_eq!(
        names,
        [
            "fork",
            "buffer_size",
            "decompression_threads",
            "compression_threads",
            "parse_threads",
            "pipeline_mode",
            "pipeline_depth",
            "render_task_records",
            "writer_strategy",
            "output_buffer_size",
            "reader_core",
            "prewarm_ahead",
            "input_shape",
            "inflated_bytes_per_record",
            "records_estimate",
            "available_parallelism",
            "physical_cores",
        ]
    );
    let value = |name: &str| lines.iter().find(|(n, _)| *n == name).unwrap().1;
    assert_eq!(value("fork"), fork);
    assert_eq!(value("buffer_size"), "5000");
    assert_eq!(value("decompression_threads"), threads);
    assert_eq!(
        value("compression_threads"),
        cpus.min(COMPRESSION_THREADS_AUTO_CAP).to_string()
    );
    assert_eq!(value("parse_threads"), "0");
    assert_eq!(value("pipeline_mode"), "streaming");
    assert_eq!(value("pipeline_depth"), "2");
    assert_eq!(value("render_task_records"), "0");
    assert_eq!(value("writer_strategy"), "buffered");
    assert_eq!(value("output_buffer_size"), "4194304");
    assert_eq!(value("reader_core"), "auto");
    assert_eq!(value("prewarm_ahead"), "1");
    assert_eq!(value("input_shape"), "indexed");
    let data_lines: Vec<String> = std::fs::read_to_string(corpus.dir.join("variants.vcf"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect();
    let mean = data_lines.iter().map(|l| l.len() as u64).sum::<u64>() / data_lines.len() as u64;
    assert_eq!(mean, 55);
    assert_eq!(value("inflated_bytes_per_record"), mean.to_string());
    let records: u64 = value("records_estimate").parse().unwrap();
    let count = data_lines.len() as u64;
    assert!(
        (count..=count + 2).contains(&records),
        "{records} records estimated against {count}"
    );
    assert_eq!(value("available_parallelism"), reported_parallelism_line());
    assert_eq!(value("physical_cores"), physical_cores_line());
    assert!(!out.exists());
}

/// The knobs print without a cache: the probe reads the configuration and the
/// input alone, so a cache directory that does not exist is never opened.
#[test]
fn resolved_config_needs_no_cache() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    let out = tmp.path().join("never_written.out");
    let (success, stdout, stderr) = resolved_config(
        corpus,
        &corpus.dir.join("variants.vcf"),
        &tmp.path().join("no_such_cache"),
        &out,
        &["--pipeline_depth", "3"],
    );
    assert!(success, "{stderr}");
    assert!(stdout.contains("pipeline_depth\t3\n"), "{stdout}");
    assert!(stdout.contains("input_shape\tunknown\n"), "{stdout}");
    assert!(!out.exists());
}

/// A named input file that is not there is refused with the run's own line
/// (Parser.pm:239) and exit 1, no knob printed, rather than read as an input of
/// unknown shape.
#[test]
fn resolved_config_refuses_a_missing_input_file() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = &corpora()[0];
    let out = tmp.path().join("never_written.out");
    let missing = tmp.path().join("no_such_input.vcf.gz");
    let (success, stdout, stderr) = resolved_config(
        corpus,
        &missing,
        &corpus.cache_dir(),
        &out,
        &["--decompression_threads", "auto"],
    );
    assert!(!success, "{stdout}");
    assert_eq!(stdout, "");
    assert_eq!(
        stderr,
        format!(
            "ERROR: [input] File \"{}\" does not exist\n",
            missing.display()
        )
    );
    assert!(!out.exists());
}

/// `--input_data` is the run's input whatever `-i` names, so the probe reads no
/// shape from the named file and `auto` resolves as the run does: the `0` rule.
#[test]
fn resolved_config_reads_no_shape_under_input_data() {
    let tmp = tempfile::tempdir().unwrap();
    let corpus = indexed_corpus();
    let out = tmp.path().join("never_written.out");
    let (success, stdout, stderr) = resolved_config(
        &corpus,
        &corpus.dir.join("variants.vcf.gz"),
        &corpus.cache_dir(),
        &out,
        &[
            "--input_data",
            "21\t30282169\t.\tC\tCTATT\t.\t.\t.",
            "--decompression_threads",
            "auto",
        ],
    );
    assert!(success, "{stderr}");
    let threads = host_parallelism().min(DECOMPRESSION_THREADS_AUTO_CAP);
    assert!(
        stdout.contains(&format!("decompression_threads\t{threads}\n")),
        "{stdout}"
    );
    assert!(stdout.contains("input_shape\tunknown\n"), "{stdout}");
    assert!(!out.exists());
}

/// The fields of `config` in declaration order, each with the text its pretty
/// `Debug` form prints for it: a field opens on a line indented by one step and
/// holds every deeper line up to the next field.
fn fields_of(config: &Config) -> Vec<(String, String)> {
    let mut fields: Vec<(String, String)> = Vec::new();
    for line in format!("{config:#?}").lines() {
        let Some(indented) = line.strip_prefix("    ") else {
            continue;
        };
        let opens_a_field = indented.split_once(": ").filter(|(name, _)| {
            !name.is_empty() && name.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
        });
        match (opens_a_field, fields.last_mut()) {
            (Some((name, value)), _) => fields.push((name.to_string(), value.to_string())),
            (None, Some((_, value))) => {
                value.push('\n');
                value.push_str(indented);
            }
            (None, None) => panic!("no field opens before {line:?}"),
        }
    }
    fields
}

/// The fields `Config::default()` is not compared in, each with why the
/// resolution of a command line gives a value of the host there.
const READ_FROM_THE_HOST: [(&str, &str); 4] = [
    ("dir", "the `~` of `~/.vep` is expanded under `HOME`"),
    ("dir_cache", "it is `dir`, expanded under `HOME`"),
    (
        "dir_plugins",
        "`<dir>/Plugins` is expanded under `HOME` and each match is tested for a directory",
    ),
    (
        "compression_threads",
        "the count of the `0` rule is lowered to the host's processors",
    ),
];

/// `vep` requires no flag, so the program name alone is the smallest command
/// line the parser accepts. Every field of the default equals the one that
/// command line resolves to, the fields read from the `Debug` text so that
/// none is named here, except the fields of `READ_FROM_THE_HOST`, where the
/// default holds a value that is the same on every machine.
#[test]
fn default_configuration_is_what_a_command_line_without_a_flag_resolves_to() {
    let resolved = fields_of(&Config::from_argv(["vep"]).expect("the program name alone parses"));
    let default = fields_of(&Config::default());
    let names = |fields: &[(String, String)]| -> Vec<String> {
        fields.iter().map(|(name, _)| name.clone()).collect()
    };
    assert_eq!(names(&default), names(&resolved));
    assert!(default.len() > 150, "{} fields read", default.len());
    for (field, why) in READ_FROM_THE_HOST {
        assert!(
            default.iter().any(|(name, _)| name == field),
            "{field} is no field of the configuration ({why})"
        );
    }
    let differing: Vec<String> = default
        .iter()
        .zip(&resolved)
        .filter(|((name, _), _)| !READ_FROM_THE_HOST.iter().any(|(field, _)| field == name))
        .filter(|((_, default), (_, resolved))| default != resolved)
        .map(|((name, default), (_, resolved))| {
            format!(
                "{name}: the default holds {default} and the command line resolves to {resolved}"
            )
        })
        .collect();
    assert!(differing.is_empty(), "{}", differing.join("\n"));

    let config = Config::default();
    assert_eq!(config.dir, "~/.vep");
    assert_eq!(config.dir_cache, "~/.vep");
    assert_eq!(config.dir_plugins.pattern, "~/.vep/Plugins");
    assert!(config.dir_plugins.dirs.is_empty());
    assert!(config.dir_plugins.not_found.is_empty());
    assert_eq!(config.compression_threads, 4);
}
