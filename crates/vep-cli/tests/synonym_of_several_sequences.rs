// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A record under a name the synonyms file gives to two or more chromosomes the
//! cache holds is not annotated: vep-rs skips it with a warning that names the
//! chromosomes, since the file does not say which of them the name stands for.
//! Ensembl VEP reads such a record on the first of them in the order `keys`
//! returns the name's synonyms (`BaseVEP.pm:711-716`), which changes between
//! runs of one command. A name the file gives to one chromosome is read as that
//! chromosome, and its rows are the corpus's Ensembl VEP rows for the same
//! position.

use std::fs;
use std::path::{Path, PathBuf};

mod common;

fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-synonyms")
}

/// The position of the corpus record `rs575418327`, a one-base substitution on
/// chromosome 21.
const POSITION: &str = "28885369";

/// The Feature and Consequence of each row Ensembl VEP wrote for that record:
/// the corpus's `expected/default.txt.gz`.
const ROWS_ON_21: [(&str, &str); 13] = [
    ("ENST00000303775", "upstream_gene_variant"),
    ("ENST00000351429", "5_prime_UTR_variant"),
    ("ENST00000447985", "downstream_gene_variant"),
    (
        "ENST00000460212",
        "5_prime_UTR_variant,NMD_transcript_variant",
    ),
    ("ENST00000825033", "non_coding_transcript_exon_variant"),
    ("ENST00000825034", "upstream_gene_variant"),
    ("ENST00000873666", "5_prime_UTR_variant"),
    ("ENST00000922530", "5_prime_UTR_variant"),
    (
        "ENST00001043897",
        "5_prime_UTR_variant,NMD_transcript_variant",
    ),
    ("ENST00001095638", "upstream_gene_variant"),
    (
        "ENST00001095639",
        "5_prime_UTR_variant,NMD_transcript_variant",
    ),
    ("ENST00001122419", "upstream_gene_variant"),
    ("ENST00001131612", "5_prime_UTR_variant"),
];

/// The REF and ALT cells of the corpus record at [`POSITION`].
fn alleles() -> (String, String) {
    let records = fs::read_to_string(corpus().join("variants.vcf")).unwrap();
    let line = records
        .lines()
        .find(|line| line.split('\t').nth(1) == Some(POSITION))
        .expect("the corpus record");
    let cells: Vec<&str> = line.split('\t').collect();
    (cells[3].to_string(), cells[4].to_string())
}

/// The output of one run: the Location, Feature and Consequence of each data
/// row, and the warnings file.
struct Run {
    rows: Vec<(String, String, String)>,
    warnings: String,
}

/// Annotates one record named `name` on the corpus's cache, which holds
/// chromosomes 1, 21 and X, under the synonyms file `synonyms`.
fn annotate(dir: &Path, name: &str, synonyms: &str) -> Run {
    let (reference, alternate) = alleles();
    let input = dir.join("input.vcf");
    let record = format!("{name}\t{POSITION}\trecord\t{reference}\t{alternate}\t.\t.\t.");
    fs::write(
        &input,
        format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{record}\n"),
    )
    .unwrap();
    let synonyms_file = dir.join("synonyms.txt");
    fs::write(&synonyms_file, synonyms).unwrap();
    let output = dir.join("output.txt");
    let warnings = dir.join("warnings.txt");
    let run = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["--offline", "--species", "homo_sapiens", "--assembly"])
        .args(["GRCh38", "--no_stats", "--force_overwrite", "--json_cache"])
        .arg(corpus().join("json_cache"))
        .arg("--synonyms")
        .arg(&synonyms_file)
        .arg("--warning_file")
        .arg(&warnings)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let rows = fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            (
                cells[1].to_string(),
                cells[4].to_string(),
                cells[6].to_string(),
            )
        })
        .collect();
    Run {
        rows,
        warnings: fs::read_to_string(&warnings).unwrap_or_default(),
    }
}

/// The rows of chromosome 21 at the record's position, under the Location `name`
/// gives them.
fn rows_on_21(name: &str) -> Vec<(String, String, String)> {
    let location = format!("{name}:{POSITION}");
    let row = |(feature, consequence): &(&str, &str)| {
        (
            location.clone(),
            feature.to_string(),
            consequence.to_string(),
        )
    };
    ROWS_ON_21.iter().map(row).collect()
}

/// `shared_name` is given to chromosomes 21 and X, both in the cache: the record
/// gets no row, in either order of the file's lines, and the warnings file holds
/// the one line that names both chromosomes in the file's order.
#[test]
fn a_name_the_synonyms_give_to_two_sequences_is_read_as_neither() {
    let (reference, alternate) = alleles();
    for (synonyms, listed) in [
        ("21\tshared_name\nX\tshared_name\n", "21, X"),
        ("X\tshared_name\n21\tshared_name\n", "X, 21"),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let run = annotate(dir.path(), "shared_name", synonyms);
        assert_eq!(run.rows, [], "{listed}");
        assert_eq!(
            run.warnings,
            format!(
                "WARNING: line 1 skipped (shared_name {POSITION} record {reference} {alternate} . . .): \
                 Chromosome shared_name is a synonym of 2 chromosomes of the cache ({listed}) and is \
                 read as none of them\n"
            ),
            "{listed}"
        );
    }
}

/// A name the file gives to three chromosomes names the first two and says that
/// more follow; a chromosome the cache lacks is not counted.
#[test]
fn a_name_given_to_three_sequences_names_the_first_two() {
    let dir = tempfile::tempdir().unwrap();
    let synonyms = "X\tshared_name\n22\tshared_name\n1\tshared_name\n21\tshared_name\n";
    let run = annotate(dir.path(), "shared_name", synonyms);
    assert_eq!(run.rows, []);
    assert!(
        run.warnings.ends_with(
            "): Chromosome shared_name is a synonym of 3 chromosomes of the cache (X, 1, ...) and \
             is read as none of them\n"
        ),
        "{}",
        run.warnings
    );
}

/// A record written `M` is renamed `MT` where neither `M` nor `chrM` names a
/// chromosome of the cache and `MT` is a synonym of one (`Parser.pm:794-796`).
/// Where the file gives `MT` to chromosomes 21 and X, the record gets no row and
/// the warning names both under `MT`, the name the record prints.
#[test]
fn an_m_renamed_to_an_mt_given_to_two_sequences_is_read_as_neither() {
    let dir = tempfile::tempdir().unwrap();
    let run = annotate(dir.path(), "M", "21\tMT\nX\tMT\n");
    assert_eq!(run.rows, []);
    assert_eq!(run.warnings.lines().count(), 1, "{}", run.warnings);
    assert!(
        run.warnings.ends_with(
            "): Chromosome MT is a synonym of 2 chromosomes of the cache (21, X) and is read as \
             none of them\n"
        ),
        "{}",
        run.warnings
    );
}

/// A name the file gives to chromosome 21 alone among the cache's chromosomes is
/// read as chromosome 21: beside a chromosome the cache lacks, and beside another
/// name of chromosome X.
#[test]
fn a_name_the_synonyms_give_to_one_sequence_is_read_as_it() {
    for synonyms in [
        "21\tonly_21\n",
        "22\tonly_21\n21\tonly_21\n",
        "X\tonly_x\n21\tonly_21\n",
    ] {
        let dir = tempfile::tempdir().unwrap();
        let run = annotate(dir.path(), "only_21", synonyms);
        assert_eq!(run.rows, rows_on_21("only_21"), "{synonyms:?}");
        assert_eq!(run.warnings, "", "{synonyms:?}");
    }
}
