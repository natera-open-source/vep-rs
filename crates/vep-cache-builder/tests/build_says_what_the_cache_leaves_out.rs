// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A build ends on one warning for each part its cache leaves out and lists
//! the same parts by name under `left_out` in `info.json`; a `vep` run on the
//! cache then prints one line where its command asks for such a part, and
//! writes the rows it writes on the same cache without the key.
//!
//! The builds read the committed models: the GFF3 slice of
//! `tests/fixtures/grch38_chr21/` with its peptides and table dumps, and the
//! GTF of the golden corpus `GRCh38-gtf-builder` with its genome.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21")
}

fn gtf_corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-gtf-builder")
}

/// The last line of a build without the core table dumps.
const NO_CORE_TABLES: &str = "WARN vep_cache_builder: No SYMBOL_SOURCE, HGNC_ID, MANE_SELECT, MANE_PLUS_CLINICAL, APPRIS, miRNA, RefSeq, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, DOMAINS or GENE_PHENO value in the cache of homo_sapiens GRCh38, and no mature_miRNA_variant from it: the build read no core table dumps; --core-dumps reads the release's (docs/cache-setup.md)";

/// The line of a build that read its model as a file source.
const NO_SOURCE: &str = "WARN vep_cache_builder: No SOURCE in the cache of homo_sapiens GRCh38: its transcripts, read from a GTF or GFF3 with no protein FASTA and no core table dumps, name no source, so a run writes no SOURCE field for them (docs/cache-from-gtf-or-gff3.md)";

/// The line of a run under `--everything` on a cache built without the core
/// table dumps.
const EVERYTHING_ASKS_FOR_CORE_TABLE_COLUMNS: &str = "INFO: --everything asks for SYMBOL_SOURCE, HGNC_ID, MANE_SELECT, MANE_PLUS_CLINICAL, APPRIS, miRNA, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, DOMAINS and GENE_PHENO, which this cache does not carry: it was built without the core table dumps (vep-cache-builder --core-dumps), so the columns stay empty and mature_miRNA_variant is never assigned";

/// The line of a run whose `--fields` list names `SOURCE` on a cache built
/// from a GTF read as the transcript model.
const FIELDS_WRITE_SOURCE: &str = "INFO: --fields writes the SOURCE column, which is empty on every transcript of this cache: a cache built from a GTF or GFF3 read as the transcript model names no source";

/// Builds a GRCh38 cache into `out` with `flags` and returns the log, each
/// line from its level on.
fn build(out: &Path, flags: &[(&str, &Path)]) -> Vec<String> {
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"));
    command.args(["--release", "116", "--assembly", "GRCh38"]);
    for (flag, path) in flags {
        command.arg(flag).arg(path);
    }
    let output = command
        .arg("--output-dir")
        .arg(out)
        .output()
        .expect("vep-cache-builder runs");
    let log = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "{log}");
    log.lines()
        .filter_map(|line| line.split_once(' '))
        .map(|(_, logged)| logged.trim_start().to_string())
        .collect()
}

fn info(cache: &Path) -> Value {
    serde_json::from_str(&std::fs::read_to_string(cache.join("info.json")).unwrap()).unwrap()
}

/// Every transcript record of the cache's shards.
fn transcripts(cache: &Path) -> Vec<Value> {
    let mut records = Vec::new();
    for chromosome in std::fs::read_dir(cache.join("transcripts")).unwrap() {
        for shard in std::fs::read_dir(chromosome.unwrap().path()).unwrap() {
            let text = std::fs::read_to_string(shard.unwrap().path()).unwrap();
            records.extend(serde_json::from_str::<Vec<Value>>(&text).unwrap());
        }
    }
    records
}

/// A copy of `cache` as a builder that writes no `left_out` key leaves it:
/// the same files, `info.json` without the key, and no checksum manifest, so
/// the changed metadata is read unchecked.
fn without_the_key(cache: &Path, copy: &Path) {
    fn copy_tree(from: &Path, to: &Path) {
        std::fs::create_dir_all(to).unwrap();
        for entry in std::fs::read_dir(from).unwrap() {
            let path = entry.unwrap().path();
            let target = to.join(path.file_name().unwrap());
            if path.is_dir() {
                copy_tree(&path, &target);
            } else {
                std::fs::copy(&path, &target).unwrap();
            }
        }
    }
    copy_tree(cache, copy);
    std::fs::remove_file(copy.join("CHECKSUMS.sha256")).unwrap();
    let mut record = info(copy);
    assert!(record.as_object_mut().unwrap().remove("left_out").is_some());
    std::fs::write(
        copy.join("info.json"),
        serde_json::to_string_pretty(&record).unwrap(),
    )
    .unwrap();
}

/// The variable that carries the flags of a run to
/// [`annotation_run_of_a_child_process`], one a line.
const RUN_FLAGS: &str = "VEP_TEST_LEFT_OUT_FLAGS";

/// Runs `vep` in-process with the flags the variable [`RUN_FLAGS`] carries,
/// and does nothing without the variable: the run of the child process
/// [`annotate`] starts.
#[test]
fn annotation_run_of_a_child_process() {
    let Ok(flags) = std::env::var(RUN_FLAGS) else {
        return;
    };
    let config = vep_cli::Config::from_argv(flags.lines()).unwrap();
    vep_cli::runner::Runner::new(config).unwrap().run().unwrap();
}

/// Annotates `input` over `cache` with `flags`, as the `vep` binary would,
/// and returns the data rows of the output and the lines the run wrote to
/// standard error. The run is this test binary started again for
/// [`annotation_run_of_a_child_process`] alone.
fn annotate(
    dir: &Path,
    name: &str,
    cache: &Path,
    input: &Path,
    flags: &[&str],
) -> (Vec<String>, Vec<String>) {
    let output = dir.join(format!("{name}.txt"));
    let mut argv = vec![
        "vep".to_string(),
        "--offline".to_string(),
        "--no_stats".to_string(),
        "--force_overwrite".to_string(),
        "--assembly".to_string(),
        "GRCh38".to_string(),
        "-i".to_string(),
        input.display().to_string(),
        "-o".to_string(),
        output.display().to_string(),
        "--json_cache".to_string(),
        cache.display().to_string(),
    ];
    argv.extend(flags.iter().map(|flag| flag.to_string()));
    let run = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "annotation_run_of_a_child_process",
            "--nocapture",
        ])
        .env(RUN_FLAGS, argv.join("\n"))
        .output()
        .expect("the test binary runs");
    let written = String::from_utf8(run.stderr).unwrap();
    assert!(run.status.success(), "{written}");
    let rows = std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect();
    (rows, written.lines().map(str::to_string).collect())
}

/// How often `line` stands among `lines`.
fn count(lines: &[String], line: &str) -> usize {
    lines.iter().filter(|written| *written == line).count()
}

/// A VCF of one substitution inside each of the first `records` transcripts
/// of `cache` that hold a translation, ten bases into the transcript.
fn records_inside_coding_transcripts(cache: &Path, input: &Path, records: usize) {
    let mut positions: Vec<(String, u64)> = transcripts(cache)
        .iter()
        .filter(|transcript| transcript["translation"].is_object())
        .map(|transcript| {
            (
                transcript["stable_id"].as_str().unwrap().to_string(),
                transcript["start"].as_u64().unwrap() + 10,
            )
        })
        .collect();
    positions.sort();
    positions.truncate(records);
    assert_eq!(positions.len(), records);
    positions.sort_by_key(|(_, position)| *position);
    let mut vcf =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (name, position) in positions {
        vcf.push_str(&format!("21\t{position}\tin_{name}\tA\tG\t.\t.\t.\n"));
    }
    std::fs::write(input, vcf).unwrap();
}

/// A build of the GFF3 slice and its peptides without the core table dumps
/// ends on the line that names the thirteen columns and the consequence term,
/// and its record names `core_tables`. A run under `--everything` on the
/// cache prints the line of the twelve columns that flag asks for once, a run
/// under `--quiet` and a run that asks for none of the columns print none,
/// and the rows are those of the same cache without the key. The same build
/// with the dumps ends on neither line, records no name, and its run prints
/// no such line.
#[test]
fn a_build_without_the_core_tables_says_so_and_a_run_names_the_columns_it_asks_for() {
    let dir = tempfile::tempdir().unwrap();
    let fixture = fixture();
    let (gff3, peptides) = (fixture.join("slice.gff3"), fixture.join("peptides.fa"));
    let model: [(&str, &Path); 2] = [("--gff3", &gff3), ("--fasta", &peptides)];
    let cache = dir.path().join("without");
    let log = build(&cache, &model);
    assert_eq!(
        log.last().map(String::as_str),
        Some(NO_CORE_TABLES),
        "{log:?}"
    );
    assert_eq!(count(&log, NO_CORE_TABLES), 1, "{log:?}");
    assert_eq!(info(&cache)["left_out"], serde_json::json!(["core_tables"]));

    let input = dir.path().join("records.vcf");
    records_inside_coding_transcripts(&cache, &input, 3);
    let (rows, written) = annotate(dir.path(), "everything", &cache, &input, &["--everything"]);
    assert_eq!(
        count(&written, EVERYTHING_ASKS_FOR_CORE_TABLE_COLUMNS),
        1,
        "{written:?}"
    );
    assert!(rows.len() >= 3, "{rows:?}");

    let keyless = dir.path().join("keyless");
    without_the_key(&cache, &keyless);
    let (rows_without_the_key, written) =
        annotate(dir.path(), "keyless", &keyless, &input, &["--everything"]);
    assert!(
        !written.iter().any(|line| line.contains("core table dumps")),
        "{written:?}"
    );
    assert_eq!(rows, rows_without_the_key);

    for flags in [
        &["--everything", "--quiet"][..],
        &["--biotype", "--canonical"],
    ] {
        let (_, written) = annotate(dir.path(), "none", &cache, &input, flags);
        assert!(
            !written.iter().any(|line| line.contains("core table dumps")),
            "{flags:?}: {written:?}"
        );
    }

    let with_tables = dir.path().join("with");
    let mut flags = model.to_vec();
    flags.push(("--core-dumps", &fixture));
    let log = build(&with_tables, &flags);
    assert!(
        !log.iter().any(|line| line.contains("in the cache of")),
        "{log:?}"
    );
    assert!(info(&with_tables).get("left_out").is_none());
    let (_, written) = annotate(
        dir.path(),
        "with_tables",
        &with_tables,
        &input,
        &["--everything"],
    );
    assert!(
        !written.iter().any(|line| line.contains("core table dumps")),
        "{written:?}"
    );
}

/// A build of the committed GTF and its genome, the model read as a file
/// source, ends on the line of the core tables and the line that the cache
/// names no source, and its record names both. A run whose `--fields` list
/// names `SOURCE` prints the line of the column once, a run that does not
/// ask for the column prints none, and the rows are those of the same cache
/// without the key. Every CDS line of the GTF names its protein, so the
/// build counts no transcript without one.
#[test]
fn a_build_of_a_gtf_model_says_it_names_no_source_and_a_run_that_writes_the_column_says_so() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = gtf_corpus();
    let cache = dir.path().join("cache");
    let log = build(
        &cache,
        &[
            ("--gtf", &corpus.join("chr21_genes.gtf.gz")),
            ("--genome-fasta", &corpus.join("reference.fa.gz")),
        ],
    );
    assert_eq!(
        log[log.len() - 2..],
        [NO_CORE_TABLES.to_string(), NO_SOURCE.to_string()],
        "{log:?}"
    );
    assert!(
        !log.iter().any(|line| line.contains("protein_id")),
        "{log:?}"
    );
    assert_eq!(
        info(&cache)["left_out"],
        serde_json::json!(["core_tables", "source"])
    );

    let input = corpus.join("variants.vcf");
    let fields = [
        "--tab",
        "--fields",
        "Uploaded_variation,Feature,Consequence,SOURCE",
    ];
    let (rows, written) = annotate(dir.path(), "source", &cache, &input, &fields);
    assert_eq!(count(&written, FIELDS_WRITE_SOURCE), 1, "{written:?}");
    assert!(!rows.is_empty());
    assert!(rows.iter().all(|row| row.ends_with("\t-")), "{rows:?}");

    let keyless = dir.path().join("keyless");
    without_the_key(&cache, &keyless);
    let (rows_without_the_key, written) =
        annotate(dir.path(), "keyless", &keyless, &input, &fields);
    assert!(
        !written.iter().any(|line| line.contains("SOURCE column")),
        "{written:?}"
    );
    assert_eq!(rows, rows_without_the_key);

    let (_, written) = annotate(
        dir.path(),
        "no_source",
        &cache,
        &input,
        &[
            "--tab",
            "--fields",
            "Uploaded_variation,Feature,Consequence",
        ],
    );
    assert!(
        !written.iter().any(|line| line.contains("SOURCE column")),
        "{written:?}"
    );
}

/// A build named one input at a time reads neither an activity file nor a
/// table of the matrices: with the regulatory-feature and the motif-feature
/// GFF3 files it ends on one line that names the activity of the regulatory
/// features and the two cells of the motif features, and its record names
/// both under `left_out`.
#[test]
fn a_build_of_named_inputs_names_both_cells_of_its_regulatory_tree() {
    let dir = tempfile::tempdir().unwrap();
    let fixture = fixture();
    let regulatory =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../vep-cli/tests/fixtures/regulatory-gff");
    let cache = dir.path().join("cache");
    let log = build(
        &cache,
        &[
            ("--gff3", &fixture.join("slice.gff3")),
            ("--fasta", &fixture.join("peptides.fa")),
            ("--core-dumps", &fixture),
            (
                "--regulatory-gff",
                &regulatory.join("regulatory_features.gff3.gz"),
            ),
            ("--motif-gff", &regulatory.join("motif_features.gff3.gz")),
        ],
    );
    assert_eq!(
        log.last().map(String::as_str),
        Some(
            "WARN vep_cache_builder: No cell-type activity of regulatory features, HIGH_INF_POS \
             or MOTIF_SCORE_CHANGE of motif features in the cache of homo_sapiens GRCh38: the \
             build reads them from no file; a conversion of Ensembl VEP's cache carries them \
             (docs/cache-conversion.md)"
        ),
        "{log:?}"
    );
    assert_eq!(
        log.iter()
            .filter(|line| line.contains("in the cache of"))
            .count(),
        1,
        "{log:?}"
    );
    assert_eq!(
        info(&cache)["left_out"],
        serde_json::json!(["cell_type_activity", "motif_scores"])
    );
}

/// The same GTF with `protein_id` taken off every CDS line builds each
/// coding transcript with a translation that names no protein: the build
/// counts them in one line after the line of the source, and with one CDS
/// line of one transcript alone left without it, counts that transcript in
/// the singular.
#[test]
fn a_build_counts_the_transcripts_whose_cds_line_names_no_protein() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = gtf_corpus();
    let mut gtf = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(corpus.join("chr21_genes.gtf.gz")).unwrap(),
    )
    .read_to_string(&mut gtf)
    .unwrap();
    let without_protein_id = |line: &str| -> String {
        let Some(at) = line.find("protein_id \"") else {
            return line.to_string();
        };
        let rest = &line[at..];
        let end = rest.find("; ").map_or(rest.len(), |end| end + 2);
        format!("{}{}", &line[..at], &rest[end..])
    };
    let genome = corpus.join("reference.fa.gz");
    let build_of = |name: &str, lines: String| -> (PathBuf, Vec<String>) {
        let model = dir.path().join(format!("{name}.gtf"));
        std::fs::write(&model, lines).unwrap();
        let cache = dir.path().join(name);
        let log = build(&cache, &[("--gtf", &model), ("--genome-fasta", &genome)]);
        (cache, log)
    };
    let counting = |log: &[String]| -> Vec<String> {
        log.iter()
            .filter(|line| line.contains("No ENSP or HGVSp for "))
            .cloned()
            .collect()
    };
    let unnamed = |cache: &Path| -> usize {
        transcripts(cache)
            .iter()
            .filter(|transcript| {
                transcript["translation"].is_object() && transcript.get("protein").is_none()
            })
            .count()
    };

    let every: String = gtf
        .lines()
        .map(|line| format!("{}\n", without_protein_id(line)))
        .collect();
    assert!(!every.contains("protein_id"));
    let (cache, log) = build_of("every", every);
    let coding = unnamed(&cache);
    assert!(coding > 1, "{coding}");
    let line = format!(
        "WARN vep_cache_builder: No ENSP or HGVSp for {coding} transcripts of the cache of \
         homo_sapiens GRCh38: the CDS line that opens each one's coding sequence names no \
         protein_id (docs/cache-from-gtf-or-gff3.md)"
    );
    assert_eq!(counting(&log), std::slice::from_ref(&line), "{log:?}");
    assert_eq!(
        log[log.len() - 3..],
        [NO_CORE_TABLES.to_string(), NO_SOURCE.to_string(), line],
        "{log:?}"
    );

    // The CDS lines of one transcript alone: the first transcript the file
    // gives a CDS line.
    let first = gtf
        .lines()
        .find(|line| line.split('\t').nth(2) == Some("CDS"))
        .and_then(|line| line.split("transcript_id \"").nth(1))
        .and_then(|rest| rest.split('"').next())
        .expect("a CDS line names its transcript")
        .to_string();
    let of_first = format!("transcript_id \"{first}\"");
    let one: String = gtf
        .lines()
        .map(|line| {
            if line.split('\t').nth(2) == Some("CDS") && line.contains(&of_first) {
                format!("{}\n", without_protein_id(line))
            } else {
                format!("{line}\n")
            }
        })
        .collect();
    let (cache, log) = build_of("one", one);
    assert_eq!(unnamed(&cache), 1);
    assert_eq!(
        counting(&log),
        [
            "WARN vep_cache_builder: No ENSP or HGVSp for 1 transcript of the cache of \
          homo_sapiens GRCh38: the CDS line that opens its coding sequence names no protein_id \
          (docs/cache-from-gtf-or-gff3.md)"
        ],
        "{log:?}"
    );
}
