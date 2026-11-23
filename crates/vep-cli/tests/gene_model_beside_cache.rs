// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A GFF or GTF transcript model named beside a cache, through the `vep`
//! binary: `--gff`, `--gtf` and `--custom file=<f>,format=gff|gtf`, any number
//! in one run. Ensembl VEP builds one annotation source per file after the
//! cache's (`AnnotationSourceAdaptor.pm:94-105` at release/116.3), each
//! annotates the input with its own transcripts
//! (`AnnotationType/Transcript.pm:92-190`), and every transcript of a file
//! carries the source's short name as `SOURCE`
//! (`AnnotationSource/File/BaseGXF.pm:635`, `OutputFactory.pm:1592`).
//!
//! The runs read the corpus `tests/golden/116/GRCh38-gene-model-beside-cache`:
//! a slice of the Ensembl cache, the GTF of eight chromosome 21 genes whose
//! transcript identifiers the cache holds too, and a GFF3 of the same models
//! under identifiers of their own.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

fn golden(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(path)
}

fn corpus(file: &str) -> String {
    golden("GRCh38-gene-model-beside-cache")
        .join(file)
        .display()
        .to_string()
}

fn fasta() -> String {
    corpus("reference.fa.gz")
}

const GTF: &str = "chr21_genes.gtf.gz";

fn gtf_as(short_name: &str) -> String {
    format!("file={},format=gtf,short_name={short_name}", corpus(GTF))
}

fn gff_as(short_name: &str) -> String {
    format!(
        "file={},format=gff,short_name={short_name}",
        corpus("novel_models.gff3.gz")
    )
}

/// One run's tab output: its `##` header lines, its column names and its
/// rows as column-to-value maps, with the standard error and the warnings
/// file it left.
struct TabRun {
    code: Option<i32>,
    stderr: String,
    header: Vec<String>,
    columns: Vec<String>,
    rows: Vec<BTreeMap<String, String>>,
    warnings: Vec<String>,
}

impl TabRun {
    fn column(&self, name: &str) -> Vec<&str> {
        self.rows.iter().map(|row| row[name].as_str()).collect()
    }

    /// The rows without the columns of `dropped`, each as its values in
    /// column order.
    fn rows_without(&self, dropped: &[&str]) -> Vec<Vec<&str>> {
        self.rows
            .iter()
            .map(|row| {
                self.columns
                    .iter()
                    .filter(|column| !dropped.contains(&column.as_str()))
                    .map(|column| row[column].as_str())
                    .collect()
            })
            .collect()
    }
}

/// Runs the binary over the corpus's records against `cache` with `flags`,
/// writing the tab format under `dir`.
fn run_on(dir: &Path, cache: &Path, flags: &[&str]) -> TabRun {
    let out = dir.join("out.txt");
    let warnings_file = format!("{}_warnings.txt", out.display());
    // A run that warns of nothing writes no warnings file: an earlier run's
    // must not stand in for it.
    for stale in [out.display().to_string(), warnings_file.clone()] {
        let _ = std::fs::remove_file(stale);
    }
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(corpus("variants.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .args([
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--no_stats",
            "--quiet",
            "--tab",
            "--symbol",
        ])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let text = std::fs::read_to_string(&out).unwrap_or_default();
    let mut header = Vec::new();
    let mut columns: Vec<String> = Vec::new();
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.starts_with("##") {
            header.push(line.to_string());
        } else if let Some(names) = line.strip_prefix('#') {
            columns = names.split('\t').map(str::to_string).collect();
        } else {
            rows.push(
                columns
                    .iter()
                    .cloned()
                    .zip(line.split('\t').map(str::to_string))
                    .collect(),
            );
        }
    }
    let warnings = std::fs::read_to_string(warnings_file)
        .unwrap_or_default()
        .lines()
        .map(str::to_string)
        .collect();
    TabRun {
        code: output.status.code(),
        stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
        header,
        columns,
        rows,
        warnings,
    }
}

/// [`run_on`] against the corpus's cache with the FASTA the files' sequence
/// is read from.
fn run(dir: &Path, flags: &[&str]) -> TabRun {
    let fasta = fasta();
    let mut all = vec!["--fasta", fasta.as_str()];
    all.extend(flags);
    let run = run_on(
        dir,
        &golden("GRCh38-gene-model-beside-cache/json_cache"),
        &all,
    );
    assert_eq!(run.code, Some(0), "{flags:?}: {}", run.stderr);
    run
}

/// Two files and the cache in one run: every row of a file's transcript
/// carries that file's short name as `SOURCE`, the cache's rows carry none,
/// and each source is listed with its header line and its own column, empty
/// on every row since a transcript source annotates no record
/// (`OutputFactory.pm:2444-2462`). The `--custom` source precedes the
/// `--gtf` one, whose option set appends its string to the list
/// (Config.pm:538-541).
#[test]
fn two_files_and_a_cache_are_told_apart_by_source() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = corpus(GTF);
    let run = run(dir.path(), &["--gtf", &gtf, "--custom", &gff_as("novel")]);

    let last_three: Vec<&str> = run.columns[run.columns.len() - 3..]
        .iter()
        .map(String::as_str)
        .collect();
    assert_eq!(last_three, ["SOURCE", "novel", GTF]);
    for line in [
        "## SOURCE : Source of transcript".to_string(),
        "## novel : [PATH]/novel_models.gff3.gz".to_string(),
        format!("## {GTF} : [PATH]/{GTF}"),
    ] {
        assert!(run.header.contains(&line), "{line}: {:?}", run.header);
    }
    assert!(run.column("novel").iter().all(|value| *value == "-"));
    assert!(run.column(GTF).iter().all(|value| *value == "-"));

    let mut by_source: BTreeMap<&str, usize> = BTreeMap::new();
    for row in &run.rows {
        let source = row["SOURCE"].as_str();
        *by_source.entry(source).or_default() += 1;
        match source {
            "novel" => assert!(row["Feature"].starts_with("NOVT"), "{row:?}"),
            _ => assert!(!row["Feature"].starts_with("NOVT"), "{row:?}"),
        }
    }
    assert_eq!(
        by_source.keys().copied().collect::<Vec<_>>(),
        ["-", GTF, "novel"]
    );
    // The two files hold the same models, so a record meets the same number
    // of transcripts in each.
    assert_eq!(by_source[GTF], by_source["novel"]);
    assert!(by_source["-"] > 0 && by_source[GTF] > 0);
}

/// `--transcript_filter` reads a file's transcripts as it reads the cache's
/// (`AnnotationSourceAdaptor.pm:399-403` hands a file source the filter):
/// `_source_cache is <short name>` keeps one file's transcripts, and a
/// filter no file transcript passes leaves the rows of a run that names no
/// file.
#[test]
fn transcript_filter_selects_among_the_sources() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = corpus(GTF);
    let both = ["--gtf", gtf.as_str(), "--custom", &gff_as("novel")];
    let with = |filter: &str| {
        let mut flags = both.to_vec();
        flags.extend(["--transcript_filter", filter]);
        run(dir.path(), &flags)
    };

    let novel = with("_source_cache is novel");
    let sources: BTreeSet<&str> = novel
        .rows
        .iter()
        .filter(|row| row["Feature"] != "-")
        .map(|row| row["SOURCE"].as_str())
        .collect();
    assert_eq!(sources.into_iter().collect::<Vec<_>>(), ["novel"]);

    let all = run(dir.path(), &both);
    let coding = with("biotype is protein_coding");
    let coding_sources: BTreeSet<&str> = coding.column("SOURCE").into_iter().collect();
    assert!(coding_sources.contains(GTF) && coding_sources.contains("novel"));
    assert!(coding.rows.len() < all.rows.len());

    let cache_alone = with("not _source_cache");
    let no_file = run(dir.path(), &[]);
    assert_eq!(
        cache_alone.rows_without(&["SOURCE", "novel", GTF]),
        no_file.rows_without(&[])
    );
}

/// `--gencode_basic` keeps a file's tagged transcripts as it keeps the
/// cache's: every transcript the file's rows name is one the cache's rows
/// name too, and untagged transcripts of either source write no row. Ensembl
/// VEP's file sources annotate every transcript under the flag, the
/// documented divergence `file_cache_honours_gencode_basic`.
#[test]
fn gencode_basic_filters_a_files_transcripts_as_the_caches() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = corpus(GTF);
    let features = |run: &TabRun, source: &str| -> BTreeSet<String> {
        run.rows
            .iter()
            .filter(|row| row["SOURCE"] == source && row["Feature"] != "-")
            .map(|row| row["Feature"].clone())
            .collect()
    };
    let all = run(dir.path(), &["--gtf", &gtf]);
    let basic = run(dir.path(), &["--gtf", &gtf, "--gencode_basic"]);
    let from_file = features(&basic, GTF);
    assert!(!from_file.is_empty());
    assert!(from_file.len() < features(&all, GTF).len());
    // The GTF and the cache are one release's, so a transcript carries the
    // tag in both or in neither.
    let from_cache = features(&basic, "-");
    assert!(
        from_file.is_subset(&from_cache),
        "{:?}",
        from_file.difference(&from_cache).collect::<Vec<_>>()
    );
    let untagged: BTreeSet<String> = features(&all, GTF)
        .difference(&from_file)
        .cloned()
        .collect();
    assert!(untagged.is_disjoint(&from_cache));
}

/// The rows of one transcript identifier come cache first, then each file's
/// in argument order, the `--custom` strings before the `--gtf` one. Ensembl
/// VEP keys such rows by a count of the transcripts it built from the file
/// before (`BaseGXF.pm:489`, `BaseVariationFeature.pm:310-321`), so its own
/// order changes with the input's other records and between runs.
#[test]
fn rows_of_one_transcript_come_cache_first_then_each_file_in_argument_order() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = corpus(GTF);
    let run = run(
        dir.path(),
        &[
            "--gtf",
            &gtf,
            "--custom",
            &gtf_as("second"),
            "--custom",
            &gtf_as("first"),
        ],
    );
    let mut order: BTreeMap<(&str, &str, &str), Vec<&str>> = BTreeMap::new();
    for row in &run.rows {
        if row["Feature"] == "-" {
            continue;
        }
        order
            .entry((&row["Uploaded_variation"], &row["Allele"], &row["Feature"]))
            .or_default()
            .push(&row["SOURCE"]);
    }
    let mut in_all_four = 0;
    for (key, sources) in &order {
        match sources.len() {
            // A transcript of the cache the files do not hold.
            1 => assert_eq!(sources, &["-"], "{key:?}"),
            4 => {
                assert_eq!(sources, &["-", "second", "first", GTF], "{key:?}");
                in_all_four += 1;
            }
            // A transcript the cache slice lacks is in the three sources of
            // the file alone.
            3 => assert_eq!(sources, &["second", "first", GTF], "{key:?}"),
            n => panic!("{key:?}: {n} rows: {sources:?}"),
        }
    }
    assert!(in_all_four > 0);
}

/// The record of the documented difference
/// `file_beside_cache_same_id_row_order`, in the one exon of
/// `ENST00000452969`, which the cache and the GTF both hold: the transcript's
/// two rows are written whole, the cache's first, alone and behind two
/// earlier records of a sorted input. Ensembl VEP 116.3 writes the same two
/// rows, the file's first alone and the cache's first behind those records:
/// it sorts them by the identifier and a count of the transcripts its file
/// source built before (`BaseGXF.pm:489`, `BaseVariationFeature.pm:310-321`,
/// `VariationFeature.pm:666`).
#[test]
fn the_cache_row_of_a_shared_transcript_comes_first_whatever_precedes_the_record() {
    const RECORD: &str = "21\t43783348\ttmem97p1_exon_snv\tA\tG\t.\t.\t.";
    const EARLIER: [&str; 2] = [
        "21\t13381749\tankrd30bp1_down_snv\tT\tC\t.\t.\t.",
        "21\t22097116\tlinc00308_up_snv\tA\tG\t.\t.\t.",
    ];
    const CACHE_ROW: &str = "tmem97p1_exon_snv\t21:43783348\tG\tENSG00000234030\tENST00000452969\tTranscript\tnon_coding_transcript_exon_variant\t226\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1";
    let file_row = format!("{CACHE_ROW};SOURCE={GTF}");
    let dir = tempfile::tempdir().unwrap();
    let rows_of_the_transcript = |name: &str, records: &[&str]| -> Vec<String> {
        let input = dir.path().join(format!("{name}.vcf"));
        let mut text =
            String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for record in records {
            text.push_str(record);
            text.push('\n');
        }
        std::fs::write(&input, text).unwrap();
        let out = dir.path().join(format!("{name}.txt"));
        let output = common::isolated_vep(dir.path())
            .arg("-i")
            .arg(&input)
            .arg("-o")
            .arg(&out)
            .arg("--json_cache")
            .arg(golden("GRCh38-gene-model-beside-cache/json_cache"))
            .args(["--assembly", "GRCh38", "--force_overwrite", "--no_stats"])
            .args(["--quiet", "--fasta", &fasta(), "--gtf", &corpus(GTF)])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(0),
            "{name}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        std::fs::read_to_string(&out)
            .unwrap()
            .lines()
            .filter(|line| {
                line.starts_with("tmem97p1_exon_snv\t") && line.contains("\tENST00000452969\t")
            })
            .map(str::to_string)
            .collect()
    };
    let expected = [CACHE_ROW.to_string(), file_row];
    assert_eq!(rows_of_the_transcript("alone", &[RECORD]), expected);
    assert_eq!(
        rows_of_the_transcript("behind", &[EARLIER[0], EARLIER[1], RECORD]),
        expected
    );
}

/// A chromosome only the file covers is annotated from the file alone: the
/// cache here holds the mitochondrion and nothing of chromosome 21, as
/// Ensembl VEP admits a record on any chromosome one of its sources lists
/// (`BaseRunner.pm:444-460`).
#[test]
fn a_chromosome_only_the_file_covers_is_annotated_from_the_file() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = corpus(GTF);
    let fasta = fasta();
    let run = run_on(
        dir.path(),
        &golden("GRCh38-plugin-gnomadmt/json_cache"),
        &["--fasta", &fasta, "--gtf", &gtf],
    );
    assert_eq!(run.code, Some(0), "{}", run.stderr);
    let transcript_rows: Vec<_> = run
        .rows
        .iter()
        .filter(|row| row["Feature"] != "-")
        .collect();
    assert!(!transcript_rows.is_empty());
    assert!(transcript_rows.iter().all(|row| row["SOURCE"] == GTF));
}

/// A model that writes a sequence otherwise than the cache names it reaches
/// the records on it: the corpus's GFF3 with its sequence written `chr21`,
/// beside the cache and the FASTA that name it `21`, gives every row the
/// committed file gives, the model's transcripts under their `SOURCE` among
/// them. Ensembl VEP reaches a source's sequence under the name the source
/// writes, with or without a leading `chr` (`BaseVEP.pm:689-737` at
/// release/116.3).
#[test]
fn a_model_written_with_a_chr_prefix_annotates_the_records_of_the_caches_sequence() {
    use std::io::Read;

    let dir = tempfile::tempdir().unwrap();
    let mut model = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(corpus("novel_models.gff3.gz")).unwrap())
        .read_to_string(&mut model)
        .unwrap();
    let prefixed: String = model
        .lines()
        .map(|line| match line.strip_prefix("21\t") {
            Some(rest) => format!("chr21\t{rest}\n"),
            None => format!("{line}\n"),
        })
        .collect();
    assert!(prefixed.contains("\nchr21\t"));
    let path = dir.path().join("chr_models.gff3");
    std::fs::write(&path, prefixed).unwrap();

    let as_committed = run(dir.path(), &["--custom", &gff_as("novel")]);
    let written_with_chr = run(
        dir.path(),
        &[
            "--custom",
            &format!("file={},format=gff,short_name=novel", path.display()),
        ],
    );
    assert!(as_committed.rows.iter().any(|row| row["SOURCE"] == "novel"));
    assert_eq!(written_with_chr.rows, as_committed.rows);
}

/// A file source needs a sequence to read and never matches exactly: a run
/// without a FASTA is refused naming the FASTA, with Ensembl VEP's text after
/// it (`BaseGXF.pm:105-106`), and a source of `type=exact` with Ensembl VEP's
/// text (`:109`), each before any output.
#[test]
fn a_file_source_needs_a_fasta_and_is_never_exact() {
    let dir = tempfile::tempdir().unwrap();
    let cache = golden("GRCh38-gene-model-beside-cache/json_cache");
    let gtf = corpus(GTF);
    let fasta = fasta();
    let exact = format!("file={gtf},format=gtf,type=exact");
    for (flags, message) in [
        (
            vec!["--gtf", gtf.as_str()],
            "ERROR: [custom] a GFF or GTF transcript source needs a FASTA file (--fasta) for its \
             sequences: vep-rs never connects to an Ensembl database (Ensembl VEP: GXF annotation \
             requires either database access (--database or --cache) or a FASTA file \
             (--fasta))\n",
        ),
        (
            vec!["--fasta", fasta.as_str(), "--custom", exact.as_str()],
            "ERROR: [custom] GXF annotation sources cannot be used as \"exact\" custom \
             annotation type\n",
        ),
    ] {
        let run = run_on(dir.path(), &cache, &flags);
        assert_eq!(run.code, Some(1), "{flags:?}");
        assert_eq!(run.stderr, message, "{flags:?}");
        assert!(run.rows.is_empty() && run.header.is_empty(), "{flags:?}");
    }
}

/// `--show_cache_info` builds the sources before it prints, so a gene-model
/// file there draws the refusals of a source being built: a run without a
/// FASTA, then a file that is not there. With both in place the cache's
/// information prints and the file is not read.
#[test]
fn show_cache_info_refuses_what_a_file_source_refuses_when_built() {
    let dir = tempfile::tempdir().unwrap();
    let cache = golden("GRCh38-gene-model-beside-cache/json_cache");
    let info = |flags: &[&str]| {
        let output = common::isolated_vep(dir.path())
            .arg("--json_cache")
            .arg(&cache)
            .args(["--assembly", "GRCh38", "--quiet", "--show_cache_info"])
            .args(flags)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        (
            output.status.code(),
            String::from_utf8_lossy(&output.stdout).into_owned(),
            String::from_utf8_lossy(&output.stderr).into_owned(),
        )
    };
    let fasta = fasta();
    let gtf = corpus(GTF);

    let (code, stdout, stderr) = info(&["--gtf", "/nonexistent/genes.gtf.gz"]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(stdout.is_empty());
    assert_eq!(
        stderr,
        "ERROR: [custom] a GFF or GTF transcript source needs a FASTA file (--fasta) for its \
         sequences: vep-rs never connects to an Ensembl database (Ensembl VEP: GXF annotation \
         requires either database access (--database or --cache) or a FASTA file (--fasta))\n"
    );

    let (code, stdout, stderr) = info(&["--fasta", &fasta, "--gtf", "/nonexistent/genes.gtf.gz"]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(stdout.is_empty());
    assert!(
        stderr.starts_with(
            "ERROR: [custom] cannot open the --custom source /nonexistent/genes.gtf.gz: "
        ),
        "{stderr}"
    );

    let (code, stdout, stderr) = info(&["--fasta", &fasta, "--gtf", &gtf]);
    assert_eq!(code, Some(0), "{stderr}");
    assert!(!stdout.is_empty());

    // A file the reader refuses prints the information too: no line of it is
    // read.
    let malformed = dir.path().join("malformed.gtf");
    std::fs::write(
        &malformed,
        "21\tsource\tgene\tx\t20\t.\t+\t.\tgene_id \"g\";\n",
    )
    .unwrap();
    let malformed = malformed.display().to_string();
    let (code, stdout, stderr) = info(&["--fasta", &fasta, "--gtf", &malformed]);
    assert_eq!(code, Some(0), "{stderr}");
    assert!(!stdout.is_empty());
}

/// `--nearest` beside a file source draws Ensembl VEP's warning
/// (`BaseGXF.pm:113-115`) and reads the nearest feature from the cache's
/// transcripts alone, as a run without the file does. Each feature type the
/// reader skips is warned of once, in Ensembl VEP's words
/// (`BaseGXF.pm:239-243`).
#[test]
fn nearest_and_skipped_line_types_warn_as_ensembl_vep_does() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = corpus(GTF);
    let with_file = run(dir.path(), &["--nearest", "symbol", "--gtf", &gtf]);
    let mut expected: Vec<String> = ["five_prime_utr", "start_codon", "three_prime_utr"]
        .iter()
        .map(|feature_type| {
            format!(
                "WARNING: Ignoring '{feature_type}' feature_type from {gtf} GFF/GTF file. This \
                 feature_type is not supported in VEP."
            )
        })
        .collect();
    expected.push(
        "WARNING: --nearest is not currently compatible with GXF annotation sources".to_string(),
    );
    let mut warnings = with_file.warnings.clone();
    warnings.sort();
    expected.sort();
    assert_eq!(warnings, expected);

    let without = run(dir.path(), &["--nearest", "symbol"]);
    let nearest = |run: &TabRun| -> BTreeSet<(String, String)> {
        run.rows
            .iter()
            .map(|row| (row["Uploaded_variation"].clone(), row["NEAREST"].clone()))
            .collect()
    };
    assert_eq!(nearest(&with_file), nearest(&without));
    assert!(without.warnings.is_empty());
}
