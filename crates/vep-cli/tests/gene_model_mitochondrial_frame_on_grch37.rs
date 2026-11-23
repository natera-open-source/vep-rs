// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A gene-model file whose mitochondrial models do not count on the FASTA
//! record of a GRCh37 run, read beside a GRCh37 cache through the `vep`
//! binary and through the library.
//!
//! On GRCh37 two mitochondrial sequences are in use, the assembly's of
//! 16,569 bases and the 16,571-base one of the UCSC hg19 naming, and a
//! file's `chrM` models can count on either. A FASTA record of the
//! assembly's length is read as `MT`, and its length does not say which
//! sequence the file counts on: its coding transcripts do. Where none of
//! them opens on a start codon and runs to its last codon without a stop
//! codon on that record, the run leaves the region's transcripts out and
//! says so once.
//!
//! Every input is made up: a model of one gene, a FASTA record that holds
//! the model's coding bases in place or two positions on, one record on
//! `MT`. The cache is the committed slice of `tests/golden/116/GRCh37`,
//! which holds `MT`.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The made-up model: one protein-coding gene of one exon on `chrM`, its
/// coding sequence at positions 11 to 28.
const MODEL: &str = "##gff-version 3\n\
chrM\tsynthetic\tgene\t1\t30\t.\t+\t.\tID=gene:MITOG1;biotype=protein_coding;Name=MITO1;version=1\n\
chrM\tsynthetic\tmRNA\t1\t30\t.\t+\t.\tID=transcript:MITOT1;Parent=gene:MITOG1;biotype=protein_coding;version=1\n\
chrM\tsynthetic\texon\t1\t30\t.\t+\t.\tParent=transcript:MITOT1;Name=MITOE1;rank=1\n\
chrM\tsynthetic\tCDS\t11\t28\t.\t+\t0\tID=CDS:MITOP1;Parent=transcript:MITOT1;protein_id=MITOP1;version=1\n";

/// The bases the model's coding sequence counts on: six codons that open on
/// `ATG` and close on `TAA`, with no stop codon of codon table 2 between.
const CODING_BASES: &str = "ATGCTTTGAATAAAATAA";

/// One record on `MT`, at position 15.
const RECORD: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
MT\t15\tmade_up\tG\tC\t.\t.\t.\n";

/// The number of bases of GRCh37's mitochondrial sequence.
const ASSEMBLY_LENGTH: usize = 16_569;

/// Where a record holds [`CODING_BASES`]: at the positions the model counts
/// on, or two positions on, as a sequence two bases longer before the model
/// holds them.
#[derive(Clone, Copy, PartialEq)]
enum Held {
    InPlace,
    TwoPositionsOn,
}

fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37/json_cache")
}

/// A FASTA under `dir` of one made-up record `name` of the assembly's
/// length, with its index: [`CODING_BASES`] at positions 11 to 28, or at 13
/// to 30, and a repeating four bases everywhere else, on lines of 60. Two
/// positions on, the record holds `GTA TGC TTT GAA TAA AAT` at the model's
/// coding span: no start codon of codon table 2 opens it and its fifth codon
/// is a stop codon.
fn fasta_of(dir: &Path, name: &str, held: Held) -> PathBuf {
    let mut bases: Vec<u8> = b"ACGT"
        .iter()
        .copied()
        .cycle()
        .take(ASSEMBLY_LENGTH)
        .collect();
    let from = match held {
        Held::InPlace => 10,
        Held::TwoPositionsOn => 12,
    };
    bases[from..from + CODING_BASES.len()].copy_from_slice(CODING_BASES.as_bytes());
    let header = format!(">{name}\n");
    let mut text = header.clone();
    for line in bases.chunks(60) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    let path = dir.join(format!("{name}_{from}_synthetic.fa"));
    std::fs::write(&path, text).unwrap();
    std::fs::write(
        format!("{}.fai", path.display()),
        format!("{name}\t{ASSEMBLY_LENGTH}\t{}\t60\t61\n", header.len()),
    )
    .unwrap();
    path
}

/// One run's rows, each as its tab-separated cells, and the lines of the
/// warnings file it left.
struct Run {
    columns: Vec<String>,
    rows: Vec<Vec<String>>,
    warnings: Vec<String>,
}

impl Run {
    /// The rows of the gene-model file's transcripts: those whose `SOURCE`
    /// cell is the file's short name.
    fn of_the_file(&self) -> Vec<&Vec<String>> {
        let source = self
            .columns
            .iter()
            .position(|column| column == "SOURCE")
            .unwrap();
        self.rows
            .iter()
            .filter(|row| row[source] == "model")
            .collect()
    }
}

/// Runs the binary on [`RECORD`] against `cache` on `assembly` with `flags`,
/// writing the tab format under `dir`.
fn run(dir: &Path, cache: &Path, assembly: &str, flags: &[&str]) -> Run {
    let input = dir.join("record.vcf");
    std::fs::write(&input, RECORD).unwrap();
    let out = dir.join("out.txt");
    let warnings_file = format!("{}_warnings.txt", out.display());
    for stale in [out.display().to_string(), warnings_file.clone()] {
        let _ = std::fs::remove_file(stale);
    }
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache)
        .args([
            "--assembly",
            assembly,
            "--force_overwrite",
            "--no_stats",
            "--quiet",
            "--tab",
        ])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{flags:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    let mut columns = Vec::new();
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.starts_with("##") {
            continue;
        }
        let cells = |line: &str| line.split('\t').map(str::to_string).collect::<Vec<_>>();
        match line.strip_prefix('#') {
            Some(names) => columns = cells(names),
            None => rows.push(cells(line)),
        }
    }
    let warnings = std::fs::read_to_string(warnings_file)
        .unwrap_or_default()
        .lines()
        .map(str::to_string)
        .collect();
    Run {
        columns,
        rows,
        warnings,
    }
}

/// The model under `dir` and the `--custom` value that names it.
fn model_under(dir: &Path) -> String {
    let path = dir.join("model.gff3");
    std::fs::write(&path, MODEL).unwrap();
    format!("file={},format=gff,short_name=model", path.display())
}

/// A synonyms table under `dir` that links `MT` and `chrM`, as the table of
/// Ensembl's GRCh37 cache does.
fn synonyms_under(dir: &Path) -> String {
    let path = dir.join("chr_synonyms.txt");
    std::fs::write(&path, "MT\tchrM\n").unwrap();
    path.display().to_string()
}

/// The model's one row on [`RECORD`] over a record that holds the coding
/// bases two positions on, where a run reads the model on it: position 15 is
/// the second base of the second codon the record gives, `TGC` to `TCC`,
/// cysteine to serine.
fn row_two_positions_on() -> Vec<String> {
    [
        "made_up",
        "MT:15",
        "C",
        "MITOG1",
        "MITOT1",
        "Transcript",
        "missense_variant",
        "15",
        "5",
        "2",
        "C/S",
        "tGc/tCc",
        "-",
        "MODERATE",
        "-",
        "1",
        "-",
        "model",
        "-",
    ]
    .map(str::to_string)
    .to_vec()
}

/// Over a record of the assembly's length that holds the model's coding
/// bases two positions on, under either name and with or without a synonyms
/// table that links the two names, the record on `MT` gets no row of the
/// model: every row it gets is one the run without the model writes, and the
/// warnings file holds the one line that names the file, the region, the
/// count and what the test found.
#[test]
fn a_chrm_model_that_does_not_read_on_a_record_of_the_assemblys_length_is_left_out_of_mt() {
    let dir = tempfile::tempdir().unwrap();
    let model = model_under(dir.path());
    let file = dir.path().join("model.gff3");
    let synonyms = synonyms_under(dir.path());
    for name in ["chrM", "MT"] {
        let fasta = fasta_of(dir.path(), name, Held::TwoPositionsOn);
        let fasta = fasta.display().to_string();
        let without_the_model = run(
            dir.path(),
            &cache(),
            "GRCh37",
            &["--fasta", &fasta, "--synonyms", &synonyms],
        );
        assert!(!without_the_model.rows.is_empty());
        let left_out = format!(
            "WARNING: Leaving out 1 transcript of {} on MT: of 1 coding transcript read from the \
             FASTA record {name}, 16,569 bases long, under codon table 2, none opens on a start \
             codon and runs to its last codon without a stop codon, so this run's FASTA does \
             not hold the sequence the file's models of that region count on",
            file.display()
        );
        for flags in [
            vec![
                "--fasta",
                &fasta,
                "--custom",
                &model,
                "--synonyms",
                &synonyms,
            ],
            vec!["--fasta", &fasta, "--custom", &model],
        ] {
            let run = run(dir.path(), &cache(), "GRCh37", &flags);
            let case = format!("{name}, {flags:?}");
            assert_eq!(run.of_the_file(), [] as [&Vec<String>; 0], "{case}");
            assert_eq!(run.rows.len(), without_the_model.rows.len(), "{case}");
            assert_eq!(run.warnings, [left_out.as_str()], "{case}");
        }
    }
}

/// Over a record that holds the coding bases at the positions the model
/// counts on, the same run reads the model: the record on `MT` gets its one
/// row and nothing is warned of.
#[test]
fn a_chrm_model_that_reads_on_a_record_of_the_assemblys_length_annotates_mt() {
    let dir = tempfile::tempdir().unwrap();
    let model = model_under(dir.path());
    let fasta = fasta_of(dir.path(), "chrM", Held::InPlace);
    let fasta = fasta.display().to_string();
    let run = run(
        dir.path(),
        &cache(),
        "GRCh37",
        &["--fasta", &fasta, "--custom", &model],
    );
    assert_eq!(run.warnings, [] as [&str; 0]);
    let rows = run.of_the_file();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0][4], "MITOT1");
}

/// A transcript the file's own line flags `cds_start_NF` is not tested: its
/// first codon need be no start codon. The same model with that tag on its
/// transcript line, over the record that holds the coding bases two
/// positions on, is read on the record as it stands, and nothing is warned
/// of.
#[test]
fn a_transcript_the_file_flags_as_incomplete_at_its_start_is_not_tested() {
    let dir = tempfile::tempdir().unwrap();
    let flagged = MODEL.replace(
        "biotype=protein_coding;version=1\nchrM\tsynthetic\texon",
        "biotype=protein_coding;version=1;tag=cds_start_NF\nchrM\tsynthetic\texon",
    );
    assert_ne!(flagged, MODEL);
    let path = dir.path().join("flagged.gff3");
    std::fs::write(&path, flagged).unwrap();
    let model = format!("file={},format=gff,short_name=model", path.display());
    let fasta = fasta_of(dir.path(), "chrM", Held::TwoPositionsOn);
    let fasta = fasta.display().to_string();
    let run = run(
        dir.path(),
        &cache(),
        "GRCh37",
        &["--fasta", &fasta, "--custom", &model],
    );
    assert_eq!(run.warnings, [] as [&str; 0]);
    assert_eq!(run.of_the_file(), [&row_two_positions_on()]);
}

/// The test is GRCh37's alone. On another assembly the two names are one
/// sequence, and the same model over the same record beside a GRCh38 cache
/// is read on the record as it stands: the record on `MT` gets the model's
/// row, its coding cells the ones the record's bases give, and nothing is
/// warned of.
#[test]
fn on_another_assembly_the_model_is_read_on_the_record_as_it_stands() {
    let dir = tempfile::tempdir().unwrap();
    let model = model_under(dir.path());
    let fasta = fasta_of(dir.path(), "chrM", Held::TwoPositionsOn);
    let fasta = fasta.display().to_string();
    let cache = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-plugin-gnomadmt/json_cache");
    let run = run(
        dir.path(),
        &cache,
        "GRCh38",
        &["--fasta", &fasta, "--custom", &model],
    );
    assert_eq!(run.warnings, [] as [&str; 0]);
    assert_eq!(run.of_the_file(), [&row_two_positions_on()]);
}

/// The library reads the model as the binary does: an annotator built from
/// the binary's flags over the same cache and model gives the record on `MT`
/// no row of the model over a record that holds its coding bases two
/// positions on, and the model's one row over a record that holds them in
/// place.
#[test]
fn the_library_leaves_the_model_out_as_the_binary_does() {
    let dir = tempfile::tempdir().unwrap();
    let model = model_under(dir.path());
    let cache = cache().display().to_string();
    for (held, rows_of_the_model) in [(Held::TwoPositionsOn, 0), (Held::InPlace, 1)] {
        let fasta = fasta_of(dir.path(), "chrM", held).display().to_string();
        let config = vep_cli::Config::from_argv([
            "vep",
            "--fork",
            "1",
            "--json_cache",
            &cache,
            "--assembly",
            "GRCh37",
            "--fasta",
            &fasta,
            "--custom",
            &model,
        ])
        .unwrap();
        let annotator = vep_cli::Annotator::from_config(&config).unwrap();
        let record = RECORD.lines().last().unwrap();
        let mut batch = vep_cli::vcf_parser::parse_vcf_line(record, false).unwrap();
        let _ = annotator.annotate_batch(&mut batch);
        assert_eq!(batch.len(), 1);
        let of_the_model = batch[0]
            .transcript_consequences
            .iter()
            .filter(|row| row.transcript_id() == "MITOT1")
            .count();
        assert_eq!(of_the_model, rows_of_the_model);
        assert!(batch[0].transcript_consequences.len() > of_the_model);
    }
}

/// The page on gene-model files states the test in the terms the warning
/// and the build's refusal use: what is tested, on which record, and what a
/// run and a build do where no coding transcript reads. The cache page's
/// table of models a build does not take as written holds the refusal's row,
/// linked to the heading the test stands under.
#[test]
fn the_gene_model_page_states_the_test_of_a_mitochondrial_region_on_grch37() {
    let page = std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/cache-from-gtf-or-gff3.md"),
    )
    .unwrap();
    for sentence in [
        "Neither the name nor the record's length says which of the two sequences a file's \
         models count on, so on GRCh37 vep-rs also tests a mitochondrial region's coding \
         transcripts on the FASTA record they are cut from.",
        "A coding sequence cut from the sequence its model counts on opens on a start codon of \
         the region's codon table and holds no stop codon before its last codon.",
        "A transcript whose first codon holds a letter other than the four bases is not tested, \
         as where an `N` pads a coding sequence that opens inside a codon.",
        "Nor is one flagged `cds_start_NF`, by the tag of its own line in the file or by the \
         attribute that a cache build's core dumps or the GTF beside its GFF3 give its record.",
        "Where the FASTA gives the coding sequence of at least one of the region's tested \
         transcripts and none of them reads so, the run leaves the region's transcripts out and \
         warns once, naming the file, the region, the record and its length.",
        "A cache build ends with an error that names them and writes no cache.",
    ] {
        assert!(page.contains(sentence), "the page lacks: {sentence}");
    }
    let builds = std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/cache-setup.md"),
    )
    .unwrap();
    let row = "| On GRCh37, a mitochondrial region none of whose tested coding transcripts, one \
               or more, reads as a coding sequence on the genome FASTA record it is cut from \
               ([the test](cache-from-gtf-or-gff3.md#what-the-run-reads)) | refused before the \
               build writes anything: `sequence region <region> of <file> does not fit the genome \
               FASTA: of <n> coding transcripts read from the FASTA record <record>, <length> \
               bases long, under codon table <table>, none opens on a start codon and runs to \
               its last codon without a stop codon`, then the two mitochondrial sequences \
               GRCh37 has in use and the two remedies: a genome FASTA that holds the sequence \
               the model counts on, or a model without the region's lines |";
    assert!(
        builds.lines().any(|line| line == row),
        "the cache page lacks the row of the build's refusal"
    );
    assert!(
        page.contains("\n### What the run reads\n"),
        "the row links a heading the gene-model page lacks"
    );
}
