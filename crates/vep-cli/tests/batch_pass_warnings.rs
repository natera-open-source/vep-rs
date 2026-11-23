// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Where the warnings file lists a warning that a `--custom` source or a
//! plugin raises while a whole batch is annotated, through the `vep` binary on
//! two golden corpora: in its record's place among the run's warnings, and for
//! no record beyond the batch a run ends on. The `--custom` source is the
//! custom corpus's promoter file, whose `enhancer` line draws the warning of a
//! feature type the source does not read; the plugin is LoFTEE reading a
//! bigWig that names no chromosome of the corpus, which draws its warning of a
//! chromosome the GERP file lacks.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

mod common;

/// The header of every input the tests write.
const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The warning a record with no consequence term draws under `--most_severe`.
const NO_TERM_WARNING: &str = "WARNING: Unable to assign consequence type";

/// A record on chromosome 21 with no alternate allele, which draws
/// [`NO_TERM_WARNING`] under `--most_severe --allow_non_variant`.
const NON_VARIANT_RECORD: &str = "21\t5000000\tnon_variant\tA\t.\t.\t.\t.";

/// A golden corpus of release 116 under `tests/golden/116/`.
fn corpus(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// One run whose annotation of a whole batch raises a warning: the corpus
/// whose cache and reference it reads, the flags that name the source, the
/// corpus record that draws the warning, with its position, and the line the
/// warnings file holds for it.
struct Raising {
    corpus: PathBuf,
    flags: Vec<String>,
    record: String,
    position: u64,
    warning: String,
}

impl Raising {
    /// The binary on `input` against the corpus's cache and reference with the
    /// source's flags, writing to `out`.
    fn run(&self, dir: &Path, input: &Path, out: &Path) -> Command {
        let mut command = common::isolated_vep(dir);
        command
            .arg("-i")
            .arg(input)
            .arg("-o")
            .arg(out)
            .args(["--offline", "--json_cache"])
            .arg(self.corpus.join("json_cache"))
            .arg("--fasta")
            .arg(self.corpus.join("reference.fa.gz"))
            .args(["--force_overwrite", "--no_stats", "--quiet"])
            .args(&self.flags);
        command
    }
}

/// The line of the corpus's `variants.vcf` that `is_the_record` accepts, with
/// its position.
fn corpus_record(corpus: &Path, is_the_record: impl Fn(&[&str]) -> bool) -> (String, u64) {
    let text = fs::read_to_string(corpus.join("variants.vcf")).unwrap();
    let line = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .find(|line| is_the_record(&line.split('\t').collect::<Vec<_>>()))
        .unwrap();
    let position = line.split('\t').nth(1).unwrap().parse().unwrap();
    (line.to_string(), position)
}

/// The custom corpus's promoter file as a `--custom` source: the corpus record
/// at `21:5028997` lies inside the file's `enhancer` line
/// (`21:5028992-5029004`), a feature type the source does not read.
fn custom_source() -> Raising {
    let corpus = corpus("GRCh38-custom");
    let file = corpus.join("custom/promoters.gff.gz");
    let (record, position) = corpus_record(&corpus, |columns| columns[1] == "5028997");
    Raising {
        flags: vec![
            "--custom".to_string(),
            format!(
                "file={},format=gff,gff_type=gencode_promoter",
                file.display()
            ),
        ],
        warning: format!(
            "WARNING: Ignoring 'enhancer' feature_type from {} GFF/GTF file. \
             This feature_type is not supported in VEP.",
            file.display()
        ),
        corpus,
        record,
        position,
    }
}

/// LoFTEE on its corpus, reading as its GERP file a bigWig whose one sequence
/// is no chromosome of the corpus: the corpus record `64605`, a stop gained in
/// a coding transcript, asks the file for an interval of chromosome 21.
fn plugin() -> Raising {
    let corpus = corpus("GRCh38-plugin-loftee");
    let bigwig = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../vep-builtins/tests/data/cons_negative_zero.bw");
    let (record, position) = corpus_record(&corpus, |columns| columns[2] == "64605");
    Raising {
        flags: vec![
            "--plugin".to_string(),
            format!(
                "LoF,human_ancestor_fa:{},gerp_bigwig:{}",
                corpus.join("plugin_data/human_ancestor.fa.gz").display(),
                bigwig.display()
            ),
        ],
        warning: format!(
            "WARNING: LoFTEE: the GERP bigWig {} holds no chromosome 21; 21:30215248-30215764 \
             counts 0 toward GERP_DIST, as Ensembl VEP's plugin counts it, and later intervals \
             on a chromosome the bigWig lacks are not reported in this run",
            bigwig.display()
        ),
        corpus,
        record,
        position,
    }
}

/// The warning the line `quoted`, record `line` of the input, draws for
/// naming chromosome `Z`, which the cache lacks.
fn unknown_chromosome_warning(line: u64, quoted: &str) -> String {
    format!(
        "WARNING: line {line} skipped ({quoted}): Chromosome Z not found in annotation sources \
         or synonyms; chromosome Z does not overlap any features"
    )
}

/// The warnings file of a complete run over `records` holds exactly
/// `expected`, at any worker count and under either `--pipeline_mode` (a run
/// with a plugin or a `--custom` source annotates batch by batch under both).
fn assert_warnings_in_input_order(raising: &Raising, records: &[&str], expected: &[&str]) {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    fs::write(&input, format!("{VCF_HEADER}{}\n", records.join("\n"))).unwrap();
    let out = dir.path().join("out.txt");
    let warnings = dir.path().join("out.txt_warnings.txt");
    let expected = format!("{}\n", expected.join("\n"));
    for mode in ["streaming", "batched"] {
        for fork in ["1", "2", "8"] {
            let _ = fs::remove_file(&warnings);
            let output = raising
                .run(dir.path(), &input, &out)
                .args(["--most_severe", "--allow_non_variant"])
                .args(["--pipeline_mode", mode, "--fork", fork])
                .output()
                .unwrap();
            let run = format!("{mode}, --fork {fork}");
            assert!(
                output.status.success(),
                "{run}: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert_eq!(fs::read_to_string(&warnings).unwrap(), expected, "{run}");
        }
    }
}

/// A warning raised while a whole batch is annotated stands in its record's
/// place among the batch's warnings: after the skip of an earlier line and the
/// warning an earlier record draws once the batch is annotated, and before the
/// warning a later record draws. It stands there too as the last of 2,002
/// records of one batch, which a `--custom` source matches in chunks of 256
/// records on the run's threads.
fn assert_stands_in_its_records_place(raising: &Raising) {
    let skip = unknown_chromosome_warning(1, "Z 1000 unknown_chromosome A G . . .");
    assert_warnings_in_input_order(
        raising,
        &[
            "Z\t1000\tunknown_chromosome\tA\tG\t.\t.\t.",
            NON_VARIANT_RECORD,
            &raising.record,
        ],
        &[&skip, NO_TERM_WARNING, &raising.warning],
    );
    assert_warnings_in_input_order(
        raising,
        &[&raising.record, NON_VARIANT_RECORD],
        &[&raising.warning, NO_TERM_WARNING],
    );
    let before = raising.position - 1_000_000;
    let plain: Vec<String> = (1..=2_000u64)
        .map(|n| format!("21\t{}\tr{n}\tA\tG\t.\t.\t.", before + n * 100))
        .collect();
    let mut records = vec![NON_VARIANT_RECORD];
    records.extend(plain.iter().map(String::as_str));
    records.push(&raising.record);
    assert_warnings_in_input_order(raising, &records, &[NO_TERM_WARNING, &raising.warning]);
}

#[test]
fn a_custom_source_warning_stands_in_its_records_place() {
    assert_stands_in_its_records_place(&custom_source());
}

#[test]
fn a_plugin_warning_stands_in_its_records_place() {
    assert_stands_in_its_records_place(&plugin());
}

/// A run that ends early leaves no warning of a record beyond the batch it
/// ended on. Of 2,400 records written as Parquet in batches of 600, record 30
/// names a chromosome the cache lacks, record 900, in the second batch, is out
/// of order and record 1,790, in the last chunk of the third, is the corpus
/// record that draws the source's warning: every run is refused at record 900
/// and its warnings file holds the skip of record 30 alone, at any worker
/// count.
#[cfg(feature = "parquet")]
fn assert_no_warning_of_a_later_batch(raising: &Raising) {
    let dir = tempfile::tempdir().unwrap();
    let before = raising.position - 2_000_000;
    let mut text = String::from(VCF_HEADER);
    for line in 1..=2_400u64 {
        text.push_str(&match line {
            30 => format!("Z\t1030\tr{line}\tA\tG\t.\t.\t.\n"),
            900 => format!("21\t5\tr{line}\tA\tG\t.\t.\t.\n"),
            1_790 => format!("{}\n", raising.record),
            _ if line < 1_790 => {
                format!("21\t{}\tr{line}\tA\tG\t.\t.\t.\n", before + line * 100)
            }
            _ => format!(
                "21\t{}\tr{line}\tA\tG\t.\t.\t.\n",
                raising.position + 100_000 + line
            ),
        });
    }
    let input = dir.path().join("in.vcf");
    fs::write(&input, text).unwrap();
    let out = dir.path().join("out.parquet");
    let warnings = dir.path().join("out.parquet_warnings.txt");
    let expected = format!(
        "{}\n",
        unknown_chromosome_warning(30, "Z 1030 r30 A G . . .")
    );
    let refusal = format!(
        "ERROR: [parquet] input is not sorted by chromosome and position at line 900 \
         (21:5 after 21:{}); sort the input or pass --no_check_variants_order",
        before + 899 * 100
    );
    for mode in ["streaming", "batched"] {
        for fork in ["1", "2", "8", "32"] {
            for repeat in 0..3 {
                let _ = fs::remove_file(&warnings);
                let output = raising
                    .run(dir.path(), &input, &out)
                    .args(["--output_format", "parquet", "--buffer_size", "600"])
                    .args(["--pipeline_mode", mode, "--fork", fork])
                    .output()
                    .unwrap();
                let run = format!("{mode}, --fork {fork}, run {repeat}");
                assert_eq!(output.status.code(), Some(1), "{run}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stderr).trim_end(),
                    refusal,
                    "{run}"
                );
                assert!(!out.exists(), "{run}: the refused run leaves no output");
                assert_eq!(fs::read_to_string(&warnings).unwrap(), expected, "{run}");
            }
        }
    }
}

#[cfg(feature = "parquet")]
#[test]
fn a_run_that_ends_early_leaves_no_custom_source_warning_of_a_later_batch() {
    assert_no_warning_of_a_later_batch(&custom_source());
}

#[cfg(feature = "parquet")]
#[test]
fn a_run_that_ends_early_leaves_no_plugin_warning_of_a_later_batch() {
    assert_no_warning_of_a_later_batch(&plugin());
}
