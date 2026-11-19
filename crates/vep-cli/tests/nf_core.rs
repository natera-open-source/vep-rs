// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `vep` command line the nf-core `ENSEMBLVEP_VEP` module renders, run as
//! the module runs it: in a task directory where the cache is staged under its
//! own name and named as `${PWD}/<name>`, writing `<prefix>.vcf.gz` that the
//! module's next step indexes with tabix. Every `main.nf` line cited here is
//! of `modules/nf-core/ensemblvep/vep/main.nf` in nf-core/modules at revision
//! `7b3452a23773f687ededb4befd037920e5a51341`, written `7b3452a` beside it.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use noodles::bgzf;
use noodles::core::Region;
use noodles::vcf;

use common::{check_html_report, corpora, isolated_vep, read_expected, Corpus};

/// The run-control corpus, whose cache sits in Ensembl VEP's layout.
fn run_control_corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-run-control")
        .expect("the GRCh38-run-control corpus")
}

/// A task directory holding the corpus input and a link to its cache root, as
/// Nextflow stages a process's `path` inputs. The input keeps the header and
/// the records without an INFO `END`: the corpus header declares none, so the
/// tabix indexer here reads the value as a string and refuses it, where
/// htslib's takes the end from the `END=` text.
fn stage_task(corpus: &Corpus, task: &Path) {
    let input: String = fs::read_to_string(corpus.dir.join("variants.vcf"))
        .unwrap()
        .lines()
        .filter(|line| {
            let has_info_end = line
                .split('\t')
                .nth(7)
                .is_some_and(|info| info.split(';').any(|field| field.starts_with("END=")));
            line.starts_with('#') || !has_info_end
        })
        .map(|line| format!("{line}\n"))
        .collect();
    fs::write(task.join("variants.vcf"), input).unwrap();
    std::os::unix::fs::symlink(corpus.cache_root(), task.join("cache_root")).unwrap();
}

/// Runs the module's script in `task` with `ext_args` as `task.ext.args`,
/// `prefix` as `task.ext.prefix` and two CPUs, and returns the path of
/// `<prefix>.vcf.gz`. The argument order is `main.nf:51-60` at `7b3452a` with no
/// FASTA (:43) and the cache branch of :48, whose `--dir_cache` is the staged
/// name under the task directory (:42).
fn run_module_script(corpus: &Corpus, task: &Path, prefix: &str, ext_args: &[&str]) -> PathBuf {
    let output_name = format!("{prefix}.vcf.gz");
    let dir_cache = task.join("cache_root");
    let mut cmd = isolated_vep(task);
    cmd.args(["-i", "variants.vcf", "-o", &output_name])
        .args(ext_args)
        .args(["--compress_output", "bgzip"])
        .args(["--assembly", &corpus.assembly, "--species", &corpus.species])
        .args(["--cache", "--cache_version", &corpus.release, "--dir_cache"])
        .arg(&dir_cache)
        .args(["--fork", "2"])
        .current_dir(task);
    let output = cmd.output().expect("vep binary runs");
    assert!(
        output.status.success(),
        "the module script failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    task.join(output_name)
}

/// A column as the VCF line writes it: the record reader returns a missing
/// value empty, which the line spells `.`.
fn column(value: &str) -> String {
    if value.is_empty() { "." } else { value }.to_string()
}

/// The lines a tabix index of `path` returns for `region`, `path` indexed as
/// `tabix <prefix>.vcf.gz` (`main.nf:49` at `7b3452a`) indexes it.
fn tabix_query(path: &Path, region: &str) -> Vec<String> {
    let index = vcf::fs::index(path).expect("the output indexes as a BGZF VCF");
    let mut reader = vcf::io::Reader::new(bgzf::io::Reader::new(fs::File::open(path).unwrap()));
    let header = reader.read_header().unwrap();
    let region: Region = region.parse().unwrap();
    reader
        .query(&header, &index, &region)
        .unwrap()
        .records()
        .map(|record| {
            let record = record.unwrap();
            let start = record.variant_start().unwrap().unwrap();
            assert!(record.quality_score().is_none());
            [
                record.reference_sequence_name().to_string(),
                start.to_string(),
                column(record.ids().as_ref()),
                record.reference_bases().to_string(),
                column(record.alternate_bases().as_ref()),
                ".".to_string(),
                column(record.filters().as_ref()),
                column(record.info().as_ref()),
            ]
            .join("\t")
        })
        .collect()
}

/// The module's command line under the corpus flags and `--vcf` writes a BGZF
/// VCF that tabix indexes, and the record a region query returns is the line
/// Ensembl VEP 116.2 wrote for it (`expected/vcf.vcf.gz` of the corpus).
#[test]
fn ensemblvep_vep_script_writes_a_tabix_indexed_vcf_of_the_reference_records() {
    let corpus = run_control_corpus();
    let task = tempfile::tempdir().unwrap();
    stage_task(&corpus, task.path());
    let mut ext_args = vec!["--vcf"];
    for flag in corpus.manifest["flags"].as_array().unwrap() {
        ext_args.push(flag.as_str().unwrap());
    }
    let out = run_module_script(&corpus, task.path(), "test", &ext_args);
    let reference: Vec<String> = read_expected(&corpus.dir, "vcf.vcf")
        .lines()
        .filter(|line| line.split('\t').nth(2) == Some("stop_del2_lost_rev"))
        .map(str::to_string)
        .collect();
    assert_eq!(reference.len(), 1, "the reference holds the record once");
    assert_eq!(tabix_query(&out, "21:30348409-30348409"), reference);
}

/// Sarek's default `ext.args` for the module (`--stats_file <prefix>.summary.html`,
/// `--vcf` and its `vep_custom_args`, with no FASTA) and its `ext.prefix`, the
/// input's stem and `_VEP.ann` (`conf/modules/annotate.config:31-46` and
/// `nextflow.config:135,138,145` of nf-core/sarek at revision
/// `786a11bb4f989751d1c4bb87908ffc6096635ff3`), run through the same script to
/// a tabix-indexed VCF, and the statistics report `--stats_file` names is
/// written in place of `<output>_summary.html`: a self-contained HTML report
/// counting the 26 of the 29 staged records that are not skipped before
/// annotation, which the module's `*.html` output (`main.nf:28` at `7b3452a`)
/// collects.
#[test]
fn sarek_default_arguments_run_through_the_module_script() {
    let corpus = run_control_corpus();
    let task = tempfile::tempdir().unwrap();
    stage_task(&corpus, task.path());
    let out = run_module_script(
        &corpus,
        task.path(),
        "variants_VEP.ann",
        &[
            "--stats_file",
            "variants_VEP.ann.summary.html",
            "--vcf",
            "--everything",
            "--filter_common",
            "--per_gene",
            "--total_length",
            "--offline",
            "--format",
            "vcf",
        ],
    );
    let ids: Vec<String> = tabix_query(&out, "21:30348409-30396602")
        .iter()
        .map(|line| line.split('\t').nth(2).unwrap().to_string())
        .collect();
    assert_eq!(ids, ["stop_del2_lost_rev", "stop_del3_lost_fwd"]);
    let html = fs::read_to_string(task.path().join("variants_VEP.ann.summary.html"))
        .expect("the report --stats_file names is written");
    let lacks = check_html_report(&html);
    assert!(lacks.is_empty(), "{lacks}");
    assert!(html.contains("<td>Variants processed</td><td>26</td>"));
    assert!(!task
        .path()
        .join("variants_VEP.ann.vcf.gz_summary.html")
        .exists());
}
