// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The build's log where standard error is not a terminal: every line is
//! plain text, `<timestamp> <level> <target>: <message>`, with no colour code,
//! and a line that counts takes each noun by its count: the singular beside
//! one, the plural beside any other.

use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Whether `text` is a UTC timestamp to the microsecond,
/// `YYYY-MM-DDTHH:MM:SS.ffffffZ`.
fn is_timestamp(text: &str) -> bool {
    text.len() == 27
        && text.bytes().enumerate().all(|(at, byte)| match at {
            4 | 7 => byte == b'-',
            10 => byte == b'T',
            13 | 16 => byte == b':',
            19 => byte == b'.',
            26 => byte == b'Z',
            _ => byte.is_ascii_digit(),
        })
}

/// Writes a one-gene GFF3 model and its peptide FASTA into `dir` and returns
/// the build of them into `dir/cache`, with `NO_COLOR` unset.
fn build(dir: &Path) -> Command {
    std::fs::write(
        dir.join("model.gff3"),
        "##gff-version 3\n\
         1\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;biotype=lncRNA\n\
         1\tx\tlnc_RNA\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=lncRNA\n\
         1\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:T1\n",
    )
    .unwrap();
    std::fs::write(
        dir.join("peptides.fa"),
        ">ENSP1 pep transcript:T1 gene:G1\nMK\n",
    )
    .unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"));
    command
        .args(["--release", "116", "--assembly", "GRCh38", "--gff3"])
        .arg(dir.join("model.gff3"))
        .arg("--fasta")
        .arg(dir.join("peptides.fa"))
        .arg("--output-dir")
        .arg(dir.join("cache"))
        .env_remove("NO_COLOR");
    command
}

/// Asserts that `log` is the whole log of [`build`] in `dir`: ten lines, each
/// a timestamp, one space and the literal rest of the line. The build reads
/// no core table dumps, so its last line names what the cache leaves out.
fn assert_plain_log(log: &str, dir: &Path) {
    let at = |name: &str| dir.join(name).display().to_string();
    let expected = [
        format!(
            " INFO vep_cache_builder: Parsing protein FASTA: {}",
            at("peptides.fa")
        ),
        " INFO vep_cache_builder: Loaded 1 protein sequence from FASTA".to_string(),
        format!(
            " INFO vep_cache_builder: Parsing GFF3: {}",
            at("model.gff3")
        ),
        " INFO vep_cache_builder: Parsed 1 gene, 1 transcript".to_string(),
        format!(
            " INFO vep_cache_builder: Building JSON cache in: {}",
            at("cache")
        ),
        " INFO vep_cache_builder::builder: chr 1: 1 transcript in 1 region".to_string(),
        " INFO vep_cache_builder: Cache built: 1 transcript across 1 chromosome, 1 region file"
            .to_string(),
        format!(
            " INFO vep_cache_builder::info: Wrote {}",
            at("cache/info.json")
        ),
        " INFO vep_cache_builder: Checksum manifest written for 2 files".to_string(),
        " WARN vep_cache_builder: No SYMBOL_SOURCE, HGNC_ID, MANE_SELECT, MANE_PLUS_CLINICAL, \
         APPRIS, miRNA, RefSeq, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, DOMAINS or \
         GENE_PHENO value in the cache of homo_sapiens GRCh38, and no mature_miRNA_variant \
         from it: the build read no core table dumps; --core-dumps reads the release's \
         (docs/cache-setup.md)"
            .to_string(),
    ];
    let lines: Vec<&str> = log.lines().collect();
    let rests: Vec<&str> = lines
        .iter()
        .map(|line| {
            let (stamp, rest) = line.split_once(' ').unwrap_or((line, ""));
            assert!(is_timestamp(stamp), "{line:?} opens with no timestamp");
            rest
        })
        .collect();
    assert_eq!(rests, expected, "{log:?}");
}

#[test]
fn a_build_log_sent_to_a_file_is_plain_text() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("build.log");
    let status = build(dir.path())
        .stderr(File::create(&path).unwrap())
        .status()
        .expect("vep-cache-builder runs");
    let log = std::fs::read_to_string(&path).unwrap();
    assert!(status.success(), "{log}");
    assert_plain_log(&log, dir.path());
}

#[test]
fn a_build_log_sent_to_a_pipe_is_plain_text() {
    let dir = tempfile::tempdir().unwrap();
    let output = build(dir.path()).output().expect("vep-cache-builder runs");
    let log = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "{log}");
    assert_plain_log(&log, dir.path());
}

/// The log of a build that succeeds, each line without its timestamp.
fn log_of(command: &mut Command) -> Vec<String> {
    let output = command.output().expect("vep-cache-builder runs");
    let log = String::from_utf8(output.stderr).unwrap();
    assert!(output.status.success(), "{log}");
    log.lines()
        .map(|line| {
            let (stamp, rest) = line.split_once(' ').unwrap_or((line, ""));
            assert!(is_timestamp(stamp), "{line:?} opens with no timestamp");
            rest.to_string()
        })
        .collect()
}

/// The lines of `log` that `target` wrote, each without its level and target.
fn lines_of(log: &[String], target: &str) -> Vec<String> {
    log.iter()
        .filter_map(|line| line.split_once(&format!(" {target}: ")))
        .map(|(_, message)| message.to_string())
        .collect()
}

/// The whole log of the build of `tests/fixtures/grch38_chr21` with its core
/// table dumps: 27 genes of 65 transcripts on one chromosome, 16 peptide
/// records, one codon-table region, no karyotype region, one gene with a
/// significant phenotype, two transcripts the dumps mark as set aside, seven
/// region files and twelve synonym pairs.
#[test]
fn a_build_with_core_dumps_takes_each_noun_by_its_count() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21");
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    let log = log_of(
        Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
            .args(["--release", "116", "--assembly", "GRCh38", "--gff3"])
            .arg(fixture.join("slice.gff3"))
            .arg("--fasta")
            .arg(fixture.join("peptides.fa"))
            .arg("--core-dumps")
            .arg(&fixture)
            .arg("--output-dir")
            .arg(&cache),
    );
    let expected = [
        format!(
            " INFO vep_cache_builder: Reading core table dumps: {}",
            fixture.display()
        ),
        " INFO vep_cache_builder::core_dumps: Core dumps: 65 transcripts, 27 genes, \
         1 codon-table region, 0 karyotype regions"
            .to_string(),
        " INFO vep_cache_builder::core_dumps: Phenotype features: 1 gene with a significant \
         phenotype"
            .to_string(),
        format!(
            " INFO vep_cache_builder: Parsing protein FASTA: {}",
            fixture.join("peptides.fa").display()
        ),
        " INFO vep_cache_builder: Loaded 16 protein sequences from FASTA".to_string(),
        format!(
            " INFO vep_cache_builder: Parsing GFF3: {}",
            fixture.join("slice.gff3").display()
        ),
        " INFO vep_cache_builder: Parsed 27 genes, 65 transcripts".to_string(),
        format!(
            " INFO vep_cache_builder: Building JSON cache in: {}",
            cache.display()
        ),
        " INFO vep_cache_builder::builder: Dropped 2 artifact or readthrough transcripts, as \
         Ensembl VEP's dumper does"
            .to_string(),
        " INFO vep_cache_builder::builder: chr 21: 63 transcripts in 7 regions".to_string(),
        " INFO vep_cache_builder: Cache built: 63 transcripts across 1 chromosome, 7 region files"
            .to_string(),
        format!(
            " INFO vep_cache_builder::info: Wrote {} (12 synonym pairs)",
            cache.join("chr_synonyms.txt").display()
        ),
        format!(
            " INFO vep_cache_builder::info: Wrote {}",
            cache.join("info.json").display()
        ),
        " INFO vep_cache_builder: Checksum manifest written for 9 files".to_string(),
    ];
    assert_eq!(log, expected);
}

/// The two lines that count a regulatory tree. One regulatory record is one
/// feature in one region, no motif feature, no matrix and one shard. Two
/// regulatory records and two motifs of a matrix each, a megabase apart, are
/// four features in two regions and two shards; one motif beside one
/// regulatory record names one matrix.
#[test]
fn the_regulatory_lines_take_each_noun_by_its_count() {
    let dir = tempfile::tempdir().unwrap();
    let regulatory = |start: u32, id: &str| {
        format!(
            "1\tEnsembl\tenhancer\t{start}\t{}\t.\t.\t.\tID={id}\n",
            start + 99
        )
    };
    let motif = |start: u32, id: &str| {
        format!(
            "1\tEnsembl\tTF_binding_site\t{start}\t{}\t.\t+\t.\tID={id};\
             binding_matrix_id=M{start};transcription_factor=TF{start}\n",
            start + 13
        )
    };
    let cases = [
        (
            "one",
            regulatory(1, "R1"),
            None,
            [
                "chr 1: 1 regulatory or motif feature in 1 region",
                "Regulatory tree written: 1 regulatory feature, 0 motif features naming \
                 0 matrices, 1 shard",
            ],
        ),
        (
            "pair",
            regulatory(1, "R1"),
            Some(motif(30, "F1")),
            [
                "chr 1: 2 regulatory and motif features in 1 region",
                "Regulatory tree written: 1 regulatory feature, 1 motif feature naming \
                 1 matrix, 1 shard",
            ],
        ),
        (
            "four",
            regulatory(1, "R1") + &regulatory(1_500_001, "R2"),
            Some(motif(30, "F1") + &motif(1_500_030, "F2")),
            [
                "chr 1: 4 regulatory and motif features in 2 regions",
                "Regulatory tree written: 2 regulatory features, 2 motif features naming \
                 2 matrices, 2 shards",
            ],
        ),
    ];
    for (name, regulatory, motifs, expected) in cases {
        let inputs = dir.path().join(name);
        std::fs::create_dir(&inputs).unwrap();
        let mut command = build(&inputs);
        let regulatory_file = inputs.join("regulatory.gff");
        std::fs::write(&regulatory_file, regulatory).unwrap();
        command.arg("--regulatory-gff").arg(&regulatory_file);
        if let Some(motifs) = motifs {
            let motif_file = inputs.join("motifs.gff");
            std::fs::write(&motif_file, motifs).unwrap();
            command.arg("--motif-gff").arg(&motif_file);
        }
        let counting: Vec<String> =
            lines_of(&log_of(&mut command), "vep_cache_builder::regulatory")
                .into_iter()
                .filter(|message| !message.starts_with("Reading "))
                .collect();
        assert_eq!(counting, expected, "{name}");
    }
}

/// The two warnings on transcripts a model's lines do not build: one whose
/// biotype no rule names and one `protein_coding` transcript with a CDS line
/// and no exon are each counted in the singular, two of either in the plural.
#[test]
fn the_skipped_transcript_warnings_take_the_noun_by_its_count() {
    let dir = tempfile::tempdir().unwrap();
    let no_biotype = |id: u32| {
        let start = id * 10_000;
        format!(
            "1\tx\tgene\t{start}\t{end}\t.\t+\t.\tID=gene:B{id}\n\
             1\tx\ttranscript\t{start}\t{end}\t.\t+\t.\tID=transcript:NB{id};Parent=gene:B{id}\n\
             1\tx\texon\t{start}\t{end}\t.\t+\t.\tParent=transcript:NB{id}\n",
            end = start + 999
        )
    };
    let no_exon = |id: u32| {
        let start = id * 10_000 + 5_000;
        format!(
            "1\tx\tgene\t{start}\t{end}\t.\t+\t.\tID=gene:E{id};biotype=protein_coding\n\
             1\tx\tmRNA\t{start}\t{end}\t.\t+\t.\tID=transcript:NE{id};Parent=gene:E{id};\
             biotype=protein_coding\n\
             1\tx\tCDS\t{start}\t{end}\t.\t+\t0\tParent=transcript:NE{id}\n",
            end = start + 899
        )
    };
    let cases = [
        (
            "one",
            no_biotype(1) + &no_exon(1),
            [
                "Unable to determine the biotype of 1 transcript, skipped: NB1",
                "No exons found for 1 protein_coding transcript, skipped: NE1",
            ],
        ),
        (
            "two",
            no_biotype(1) + &no_exon(1) + &no_biotype(2) + &no_exon(2),
            [
                "Unable to determine the biotype of 2 transcripts, skipped: NB1, NB2",
                "No exons found for 2 protein_coding transcripts, skipped: NE1, NE2",
            ],
        ),
    ];
    for (name, lines, expected) in cases {
        let inputs = dir.path().join(name);
        std::fs::create_dir(&inputs).unwrap();
        let mut command = build(&inputs);
        let model = inputs.join("model.gff3");
        let built = std::fs::read_to_string(&model).unwrap();
        std::fs::write(&model, built + &lines).unwrap();
        let warnings = lines_of(&log_of(&mut command), "vep_cache_builder::gff3");
        assert_eq!(warnings, expected, "{name}");
    }
}

/// A GTF read for its tags beside the GFF3 model is logged by the transcripts
/// its lines tag: one line, one transcript.
#[test]
fn the_gtf_tags_line_counts_one_transcript_in_the_singular() {
    let dir = tempfile::tempdir().unwrap();
    let mut command = build(dir.path());
    let gtf = dir.path().join("tags.gtf");
    std::fs::write(
        &gtf,
        "1\tx\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; \
         tag \"basic\";\n",
    )
    .unwrap();
    command.arg("--gtf").arg(&gtf);
    let log = log_of(&mut command);
    let read: Vec<&String> = log
        .iter()
        .filter(|line| line.contains("Read tags for"))
        .collect();
    assert_eq!(
        read,
        [" INFO vep_cache_builder: Read tags for 1 transcript from GTF"],
        "{log:?}"
    );
}
