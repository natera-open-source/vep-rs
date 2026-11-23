// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A gene-model file that writes its mitochondrial sequence `chrM`, read
//! beside a GRCh37 cache through the `vep` binary and through the library.
//!
//! On GRCh37 `MT` names the assembly's mitochondrial sequence, of 16,569
//! bases, and `chrM` also names the one of the UCSC hg19 naming, of 16,571.
//! A file's `chrM` models are read as the cache's `MT` only where the run's
//! FASTA gives a record of the assembly's length for `MT`, or gives no
//! record for `MT` and names no organelle genome. Over any other FASTA the
//! run leaves them out and says so once: over a record of another length
//! their sequences would be cut from a sequence the cache's `MT` is not.
//!
//! Every input is made up: a model of one gene, a FASTA record of the name
//! and the length each case gives, one record on `MT`. The cache is the committed slice of
//! `tests/golden/116/GRCh37`, which holds `MT`.

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

/// The bases the made-up FASTA holds at positions 11 to 28: six codons, the
/// second of them `CTT`.
const CODING_BASES: &str = "ATGCTTTGAATAAAATAA";

/// One record on `MT`: position 15, the second base of the model's second
/// codon, `T` to `C`.
const RECORD: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
MT\t15\tmade_up\tT\tC\t.\t.\t.\n";

/// The number of bases of GRCh37's mitochondrial sequence.
const ASSEMBLY_LENGTH: usize = 16_569;

/// The number of bases of the mitochondrial sequence the UCSC hg19 naming
/// calls `chrM`.
const OTHER_LENGTH: usize = 16_571;

fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37/json_cache")
}

/// A FASTA under `dir` of one made-up record `name` of `length` bases, with
/// its index: [`CODING_BASES`] at positions 11 to 28 and a repeating four
/// bases everywhere else, on lines of 60.
fn fasta_of(dir: &Path, name: &str, length: usize) -> PathBuf {
    let mut bases: Vec<u8> = b"ACGT".iter().copied().cycle().take(length).collect();
    bases[10..28].copy_from_slice(CODING_BASES.as_bytes());
    let header = format!(">{name}\n");
    let mut text = header.clone();
    for line in bases.chunks(60) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    let path = dir.join(format!("{name}_{length}_synthetic.fa"));
    std::fs::write(&path, text).unwrap();
    std::fs::write(
        format!("{}.fai", path.display()),
        format!("{name}\t{length}\t{}\t60\t61\n", header.len()),
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

/// Runs the binary on [`RECORD`] against the GRCh37 cache with `flags`,
/// writing the tab format under `dir`.
fn run(dir: &Path, flags: &[&str]) -> Run {
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
        .arg(cache())
        .args([
            "--assembly",
            "GRCh37",
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

/// The warning of a run that leaves the model's one transcript out.
fn left_out(model: &str) -> String {
    let file = model
        .strip_prefix("file=")
        .and_then(|rest| rest.split(',').next())
        .unwrap();
    format!(
        "WARNING: Leaving out 1 transcript of {file} on chrM: on GRCh37 a gene-model file's M or \
         chrM is read as MT only where the FASTA gives for MT the assembly's mitochondrial \
         sequence, 16,569 bases long, or gives no record for MT and names no organelle genome, \
         and this run's FASTA does neither"
    )
}

/// Over a FASTA whose `chrM` has the assembly's 16,569 bases the file's
/// `chrM` is the cache's `MT`: the record on `MT` gets the model's row, its
/// coding cells read from that record under the mitochondrial code, and
/// nothing is warned of. The row is derived from the made-up inputs:
/// position 15 is base 15 of the one exon and base 5 of the coding sequence,
/// the second base of the second codon, `CTT` to `CCT`, leucine to proline.
#[test]
fn a_chrm_model_over_a_record_of_the_assemblys_length_annotates_mt() {
    let dir = tempfile::tempdir().unwrap();
    let model = model_under(dir.path());
    let synonyms = synonyms_under(dir.path());
    let fasta = fasta_of(dir.path(), "chrM", ASSEMBLY_LENGTH);
    let fasta = fasta.display().to_string();
    for flags in [
        vec!["--fasta", &fasta, "--custom", &model],
        vec![
            "--fasta",
            &fasta,
            "--custom",
            &model,
            "--synonyms",
            &synonyms,
        ],
    ] {
        let run = run(dir.path(), &flags);
        assert_eq!(run.warnings, [] as [&str; 0], "{flags:?}");
        assert_eq!(
            run.columns,
            [
                "Uploaded_variation",
                "Location",
                "Allele",
                "Gene",
                "Feature",
                "Feature_type",
                "Consequence",
                "cDNA_position",
                "CDS_position",
                "Protein_position",
                "Amino_acids",
                "Codons",
                "Existing_variation",
                "IMPACT",
                "DISTANCE",
                "STRAND",
                "FLAGS",
                "SOURCE",
                "model",
            ],
            "{flags:?}"
        );
        let expected: Vec<String> = [
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
            "L/P",
            "cTt/cCt",
            "-",
            "MODERATE",
            "-",
            "1",
            "-",
            "model",
            "-",
        ]
        .map(str::to_string)
        .to_vec();
        assert_eq!(run.of_the_file(), [&expected], "{flags:?}");
    }
}

/// Over a FASTA whose `chrM` has 16,571 bases the file's `chrM` is not the
/// cache's `MT`, with and without a synonyms table that links the two names:
/// the record on `MT` gets no row of the model, every row it gets is one the
/// run without the model writes, and the warnings file holds the one line
/// that names the file, the region and the count.
#[test]
fn a_chrm_model_over_a_record_of_another_length_is_left_out_of_mt() {
    let dir = tempfile::tempdir().unwrap();
    let model = model_under(dir.path());
    let synonyms = synonyms_under(dir.path());
    let fasta = fasta_of(dir.path(), "chrM", OTHER_LENGTH);
    let fasta = fasta.display().to_string();
    let without_the_model = run(dir.path(), &["--fasta", &fasta, "--synonyms", &synonyms]);
    assert!(!without_the_model.rows.is_empty());
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
        let run = run(dir.path(), &flags);
        assert_eq!(run.of_the_file(), [] as [&Vec<String>; 0], "{flags:?}");
        assert_eq!(run.rows.len(), without_the_model.rows.len(), "{flags:?}");
        assert_eq!(run.warnings, [left_out(&model)], "{flags:?}");
    }
}

/// Over a FASTA that gives no record for `MT` and names no organelle genome,
/// one record named `21` as a part of a genome holds, the file's `chrM` is
/// the cache's `MT` whatever the record's length, and nothing is warned of.
/// The record on `MT` gets the model's row, derived from the made-up inputs:
/// position 15 is base 15 of the one exon and base 5 of the coding sequence,
/// in its second codon, and with no bases read for `MT` the row names a
/// change of the coding sequence and states no codon and no residue.
#[test]
fn a_chrm_model_over_a_fasta_with_no_record_for_mt_and_no_organelle_genome_annotates_mt() {
    for length in [ASSEMBLY_LENGTH, OTHER_LENGTH] {
        let dir = tempfile::tempdir().unwrap();
        let model = model_under(dir.path());
        let synonyms = synonyms_under(dir.path());
        let fasta = fasta_of(dir.path(), "21", length);
        let fasta = fasta.display().to_string();
        for flags in [
            vec!["--fasta", &fasta, "--custom", &model],
            vec![
                "--fasta",
                &fasta,
                "--custom",
                &model,
                "--synonyms",
                &synonyms,
            ],
        ] {
            let run = run(dir.path(), &flags);
            assert_eq!(run.warnings, [] as [&str; 0], "{length} bases, {flags:?}");
            let expected: Vec<String> = [
                "made_up",
                "MT:15",
                "C",
                "MITOG1",
                "MITOT1",
                "Transcript",
                "coding_sequence_variant",
                "15",
                "5",
                "2",
                "-",
                "-",
                "-",
                "MODIFIER",
                "-",
                "1",
                "-",
                "model",
                "-",
            ]
            .map(str::to_string)
            .to_vec();
            assert_eq!(run.of_the_file(), [&expected], "{length} bases, {flags:?}");
        }
    }
}

/// Over a FASTA that gives no record for `MT` and names an organelle genome
/// of its own, one record named `Mito`, the file's `chrM` is not the cache's
/// `MT` whatever the record's length: the record on `MT` gets no row of the
/// model and the warnings file holds the one line.
#[test]
fn a_chrm_model_over_a_fasta_that_names_an_organelle_genome_of_its_own_is_left_out_of_mt() {
    for length in [ASSEMBLY_LENGTH, OTHER_LENGTH] {
        let dir = tempfile::tempdir().unwrap();
        let model = model_under(dir.path());
        let fasta = fasta_of(dir.path(), "Mito", length);
        let fasta = fasta.display().to_string();
        let run = run(dir.path(), &["--fasta", &fasta, "--custom", &model]);
        assert_eq!(run.of_the_file(), [] as [&Vec<String>; 0], "{length} bases");
        assert_eq!(run.warnings, [left_out(&model)], "{length} bases");
    }
}

/// On another assembly the two names are one sequence whatever its length:
/// the same model over the same 16,571-base record beside a GRCh38 cache
/// annotates the record on `MT`, and nothing is warned of.
#[test]
fn on_another_assembly_a_chrm_model_annotates_mt_over_a_record_of_any_length() {
    let dir = tempfile::tempdir().unwrap();
    let model = model_under(dir.path());
    let fasta = fasta_of(dir.path(), "chrM", OTHER_LENGTH);
    let input = dir.path().join("record.vcf");
    std::fs::write(&input, RECORD).unwrap();
    let out = dir.path().join("out.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../tests/golden/116/GRCh38-plugin-gnomadmt/json_cache"),
        )
        .args(["--assembly", "GRCh38", "--force_overwrite", "--no_stats"])
        .args(["--quiet", "--tab", "--custom", &model, "--fasta"])
        .arg(&fasta)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    let of_the_model: Vec<&str> = text
        .lines()
        .filter(|line| line.split('\t').nth(4) == Some("MITOT1"))
        .collect();
    assert_eq!(of_the_model.len(), 1, "{text}");
    assert!(
        of_the_model[0].contains("\tmissense_variant\t15\t5\t2\tL/P\tcTt/cCt\t"),
        "{}",
        of_the_model[0]
    );
    assert!(!Path::new(&format!("{}_warnings.txt", out.display())).exists());
}

/// The library reads the model as the binary does: an annotator built from
/// the binary's flags over the same cache, model and synonyms table gives the
/// record on `MT` no row of the model over a `chrM` record of 16,571 bases,
/// and the model's one row over a `chrM` record of the assembly's length.
#[test]
fn the_library_leaves_a_chrm_model_out_of_mt_as_the_binary_does() {
    let dir = tempfile::tempdir().unwrap();
    let model = model_under(dir.path());
    let synonyms = synonyms_under(dir.path());
    let cache = cache().display().to_string();
    for (length, rows_of_the_model) in [(OTHER_LENGTH, 0), (ASSEMBLY_LENGTH, 1)] {
        let fasta = fasta_of(dir.path(), "chrM", length).display().to_string();
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
            "--synonyms",
            &synonyms,
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
        assert_eq!(of_the_model, rows_of_the_model, "{length} bases");
        assert!(batch[0].transcript_consequences.len() > of_the_model);
    }
}

/// The page on gene-model files states the rule in the terms the warning
/// uses: the two names, the two lengths, what a run does and what a build
/// does.
#[test]
fn the_gene_model_page_states_the_rule_of_a_chrm_region_on_grch37() {
    let page = std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/cache-from-gtf-or-gff3.md"),
    )
    .unwrap();
    for sentence in [
        "A mitochondrion the file writes `chrM` or `M` is read as `MT`.",
        "On GRCh37 those two names also name the 16,571-base mitochondrial sequence of the UCSC \
         hg19 naming.",
        "There the file's `chrM` or `M` is read as `MT` only where the FASTA holds the \
         assembly's mitochondrial sequence, 16,569 bases long, or gives no record for `MT` and \
         names no organelle genome.",
        "Over any other record the run leaves the region's transcripts out and warns once, \
         naming the file, the region and the count.",
        "A cache build keeps such a region under its written name and says so.",
    ] {
        assert!(page.contains(sentence), "the page lacks: {sentence}");
    }
    assert!(
        !page.contains("mitochondrial sequence, 16,569 bases long."),
        "the page reads the name as `MT` over the assembly's sequence alone"
    );
}
