// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The `vep` command line the nf-core `ENSEMBLVEP_VEP` module renders, run as
//! the module runs it: in a task directory where the cache is staged under its
//! own name and named as `${PWD}/<name>`, writing `<prefix>.vcf.gz` that the
//! module's next step indexes with tabix.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::{corpora, isolated_vep, Corpus};
use noodles::bgzf;
use noodles::core::Region;
use noodles::vcf;

/// The run-control corpus, whose cache sits in Ensembl VEP's layout.
fn run_control_corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-run-control")
        .expect("the GRCh38-run-control corpus")
}

/// A task directory holding the corpus input and a link to its cache root, as
/// Nextflow stages a process's `path` inputs. The input keeps the header and
/// the records whose ALT is bases alone: the tabix indexer here computes no end
/// for a symbolic or breakend allele, where htslib's falls back to the REF.
fn stage_task(corpus: &Corpus, task: &Path) {
    let input: String = fs::read_to_string(corpus.dir.join("variants.vcf"))
        .unwrap()
        .lines()
        .filter(|line| {
            line.starts_with('#')
                || line
                    .split('\t')
                    .nth(4)
                    .is_some_and(|alt| !alt.contains(['<', '[', ']', '*', '.']))
        })
        .map(|line| format!("{line}\n"))
        .collect();
    fs::write(task.join("variants.vcf"), input).unwrap();
    std::os::unix::fs::symlink(corpus.cache_root(), task.join("cache_root")).unwrap();
}

/// Runs the module's script in `task` with `ext_args` as `task.ext.args` and
/// two CPUs, and returns the path of `test.vcf.gz`. The argument order is
/// `ensemblvep/vep/main.nf:51-60` with no FASTA and the cache branch of :48.
fn run_module_script(corpus: &Corpus, task: &Path, ext_args: &[&str]) -> PathBuf {
    let dir_cache = task.join("cache_root");
    let mut cmd = isolated_vep(task);
    cmd.args(["-i", "variants.vcf", "-o", "test.vcf.gz"])
        .args(ext_args)
        .args(["--compress_output", "bgzip"])
        .args(["--assembly", &corpus.assembly, "--species", "homo_sapiens"])
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
    task.join("test.vcf.gz")
}

/// The lines a tabix index of `path` returns for `region`, `path` indexed as
/// `tabix test.vcf.gz` (`ensemblvep/vep/main.nf:49`) indexes it.
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
            let filters = match record.filters().as_ref() {
                "" => ".".to_string(),
                written => written.to_string(),
            };
            [
                record.reference_sequence_name().to_string(),
                start.to_string(),
                record.ids().as_ref().to_string(),
                record.reference_bases().to_string(),
                record.alternate_bases().as_ref().to_string(),
                ".".to_string(),
                filters,
                record.info().as_ref().to_string(),
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
    let flags: Vec<String> = corpus.manifest["flags"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f.as_str().unwrap().to_string())
        .collect();
    ext_args.extend(flags.iter().map(String::as_str));
    let out = run_module_script(&corpus, task.path(), &ext_args);
    assert_eq!(
        tabix_query(&out, "21:30348409-30348409"),
        [
            "21\t30348409\tstop_del2_lost_rev\tCTA\tC\t.\t.\tCSQ=-|stop_lost|HIGH|KRTAP23-1|\
          ENSG00000186980|Transcript|ENST00000334160|protein_coding|1/1||||196-197|196-197|66|\
          */X|TAg/g|||-1||HGNC|HGNC:18928|ENSP00000346536.3"
        ]
    );
}

/// Sarek's default `ext.args` for the module (`--stats_file`, `--vcf` and its
/// `vep_custom_args`) run through the same script to a tabix-indexed VCF, and
/// the statistics file it names is not written.
#[test]
fn sarek_default_arguments_run_through_the_module_script() {
    let corpus = run_control_corpus();
    let task = tempfile::tempdir().unwrap();
    stage_task(&corpus, task.path());
    let out = run_module_script(
        &corpus,
        task.path(),
        &[
            "--stats_file",
            "test_VEP.ann.summary.html",
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
    assert!(!task.path().join("test_VEP.ann.summary.html").exists());
}
