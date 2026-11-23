// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The documented difference
//! `mitochondrial_record_of_another_sequence_placed_on_mt`, through the `vep`
//! binary: a GRCh37 run over a cache whose synonym table links `MT` to `chrM`,
//! on input records written `chrM` or `M`, with a reference FASTA whose
//! record of that pair of names has another length than the assembly's
//! 16,569 bases. The records are the reproduction's own,
//! `scripts/divergences/reproductions/mitochondrial_record_of_another_sequence_placed_on_mt.toml`,
//! and the FASTA is built as its prepare steps build it.
//!
//! Ensembl VEP names the cache's chromosome for an input record from the
//! record's name alone, a synonym the cache holds first
//! (`AnnotationSource.pm:179`, `BaseVEP.pm:708-716`), and so writes the rows
//! of the transcripts of `MT` at the record's number. Under a reference check
//! it names the record's slice after the FASTA record (`BaseVEP.pm:413` and
//! `:425-433`), finds no transcript on a slice of that name
//! (`AnnotationType/Transcript.pm:141`) and writes one `intergenic_variant`
//! row.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

const REPRODUCTION: &str =
    "scripts/divergences/reproductions/mitochondrial_record_of_another_sequence_placed_on_mt.toml";
const CACHE: &str = "tests/fixtures/bayesdel-mitochondrial-contig/json_cache";
const REFERENCE: &str = "tests/fixtures/mitochondrial-frame-check/reference_synthetic.fa";

/// The protein-coding transcript whose first exon holds the records'
/// positions on the cache's `MT`.
const CODING: &str = "ENST00000361390";

/// The reproduction's records: the identifier, the chromosome as written, the
/// position and the alleles.
const RECORDS: [(&str, &str, u32, &str, &str); 4] = [
    ("named_chrm_other_base", "chrM", 3349, "A", "C"),
    ("named_chrm_same_base", "chrM", 3350, "T", "A"),
    ("named_m_other_base", "M", 3349, "A", "C"),
    ("named_m_same_base", "M", 3350, "T", "A"),
];

/// The skip line of the record on input line `line`, over a FASTA whose
/// record of 16,571 bases is named `held`.
fn skip_line(line: usize, held: &str) -> String {
    let (id, chromosome, position, reference, alternate) = RECORDS[line - 1];
    format!(
        "WARNING: line {line} skipped ({chromosome} {position} {id} {reference} {alternate} . . .): \
         Chromosome {chromosome} is the FASTA record {held} of 16571 bases, another sequence than \
         the cache's MT of 16569 bases, and is read as no chromosome of the cache\n"
    )
}

/// The warnings file of a run whose FASTA names its record of 16,571 bases
/// `held`: one skip line a record, in input order.
fn skip_lines(held: &str) -> String {
    (1..=RECORDS.len())
        .map(|line| skip_line(line, held))
        .collect()
}

/// The row of the record on input line `line` on the coding transcript, in a
/// run that places it on the cache's `MT`: the record lies on the first exon's
/// codon 15, at its first or second base. Ensembl VEP 116.3 writes that row
/// for a record written `chrM`; a record written `M` it writes as one
/// `intergenic_variant` row, and the row here is vep-rs's reading of a kept
/// `M` as `MT`.
fn row_on_mt(line: usize) -> String {
    let (id, chromosome, position, _, alternate) = RECORDS[line - 1];
    let (cdna, amino_acids, codons) = match position {
        3349 => (43, "I/L", "Atc/Ctc"),
        _ => (44, "I/N", "aTc/aAc"),
    };
    format!(
        "{id}\t{chromosome}:{position}\t{alternate}\tENSG00000198888\tENST00000361390\tTranscript\t\
         missense_variant\t{cdna}\t{cdna}\t15\t{amino_acids}\t{codons}\t-\tIMPACT=MODERATE;STRAND=1"
    )
}

/// The one row a run under `--dont_skip` writes for the record on input line
/// `line`: no feature of the cache lies on the record's sequence.
fn row_on_no_feature(line: usize) -> String {
    let (id, chromosome, position, _, alternate) = RECORDS[line - 1];
    format!(
        "{id}\t{chromosome}:{position}\t{alternate}\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\t\
         IMPACT=MODIFIER"
    )
}

fn crate_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The lines of the TOML basic multi-line string or array that opens on the
/// line starting with `opening`, up to the line that is `closing`.
fn block<'a>(text: &'a str, opening: &str, closing: &str) -> Vec<&'a str> {
    text.lines()
        .skip_while(|line| !line.starts_with(opening))
        .skip(1)
        .take_while(|line| *line != closing)
        .collect()
}

fn reproduction() -> String {
    fs::read_to_string(crate_dir().join("../..").join(REPRODUCTION)).unwrap()
}

/// The reproduction's input records as a VCF file in `dir`, held to
/// [`RECORDS`].
fn records(dir: &Path) -> PathBuf {
    let text = reproduction();
    let lines = block(&text, "text = \"\"\"", "\"\"\"");
    let written: Vec<String> = RECORDS
        .iter()
        .map(|(id, chromosome, position, reference, alternate)| {
            format!("{chromosome}\t{position}\t{id}\t{reference}\t{alternate}\t.\t.\t.")
        })
        .collect();
    assert_eq!(lines, written, "the reproduction's records");
    let vcf = dir.join("records.vcf");
    fs::write(
        &vcf,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}\n",
            lines.join("\n")
        ),
    )
    .unwrap();
    vcf
}

/// The made-up reference's bases, its one record's lines joined.
fn reference_bases() -> String {
    fs::read_to_string(crate_dir().join(REFERENCE))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect()
}

/// An indexed FASTA in `dir` of one record `name` holding `bases`, sixty to a
/// line, with the index line `samtools faidx` writes for it.
fn fasta(dir: &Path, file: &str, name: &str, bases: &str) -> PathBuf {
    let mut text = format!(">{name}\n");
    for line in bases.as_bytes().chunks(60) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    let path = dir.join(file);
    fs::write(&path, text).unwrap();
    fs::write(
        dir.join(format!("{file}.fai")),
        format!("{name}\t{}\t{}\t60\t61\n", bases.len(), name.len() + 2),
    )
    .unwrap();
    path
}

/// The file the reproduction prepares, under the record name `name`: the
/// reference with its first two bases put in front once more, so its base at
/// a position is the reference's base two positions back.
fn two_bases_longer(dir: &Path, name: &str) -> PathBuf {
    let reference = reference_bases();
    assert_eq!(reference.len(), 16_569);
    let longer = format!("{}{reference}", &reference[..2]);
    fasta(
        dir,
        &format!("{name}_two_bases_longer_synthetic.fa"),
        name,
        &longer,
    )
}

/// The rows and the warnings file of a `vep` run over the reproduction's
/// records and cache with `fasta` as its reference.
fn run(dir: &Path, label: &str, fasta: &Path, flags: &[&str]) -> (Vec<String>, String) {
    let out = dir.join(format!("{label}.txt"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(records(dir))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(crate_dir().join(CACHE))
        .arg("--fasta")
        .arg(fasta)
        .args([
            "--offline",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh37",
            "--format",
            "vcf",
            "--buffer_size",
            "5000",
            "--force_overwrite",
            "--quiet",
            "--no_stats",
        ])
        .args(flags)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{label}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let rows = fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect();
    let warnings =
        fs::read_to_string(dir.join(format!("{label}.txt_warnings.txt"))).unwrap_or_default();
    (rows, warnings)
}

/// The one row of `rows` for the record `id` on `transcript`.
fn row_of<'a>(rows: &'a [String], id: &str, transcript: &str) -> &'a str {
    let mut found = rows.iter().filter(|row| {
        let cells: Vec<&str> = row.split('\t').collect();
        cells[0] == id && cells[4] == transcript
    });
    let row = found.next().expect("a row of the record on the transcript");
    assert!(found.next().is_none(), "one row of {id} on {transcript}");
    row
}

/// The strings the reproduction records, in the `vep_rs` list that follows
/// the table header `table`, each less its quotes and its closing comma.
fn recorded_for_vep_rs<'a>(text: &'a str, after: &str) -> Vec<&'a str> {
    let from = text.find(after).expect("the run's table");
    block(&text[from..], "vep_rs = [", "]")
        .into_iter()
        .map(|line| line.trim().trim_end_matches(',').trim_matches('"'))
        .collect()
}

/// A record written `chrM` or `M` over a FASTA whose record of that pair has
/// another length than the assembly's lies on that sequence, which the cache
/// does not hold: the run writes no row for it, with or without a reference
/// check or a notation flag, and one skip line that names the FASTA record
/// and both lengths. The warnings file is the one the reproduction records
/// for vep-rs.
#[test]
fn a_record_named_for_another_mitochondrial_sequence_is_kept_off_the_cache() {
    let dir = tempfile::tempdir().unwrap();
    let longer = two_bases_longer(dir.path(), "chrM");
    let every_flag_set: [(&str, &[&str]); 4] = [
        ("plain", &[]),
        ("notations", &["--hgvs", "--hgvsg"]),
        ("lookup_ref", &["--lookup_ref", "--show_ref_allele"]),
        ("check_ref", &["--check_ref"]),
    ];
    for (label, flags) in every_flag_set {
        let (rows, warnings) = run(dir.path(), label, &longer, flags);
        assert_eq!(rows, Vec::<String>::new(), "{label}");
        assert_eq!(warnings, skip_lines("chrM"), "{label}");
    }

    let text = reproduction();
    let recorded: String = recorded_for_vep_rs(&text, "kind = \"warnings\"")
        .into_iter()
        .map(|line| format!("{line}\n"))
        .collect();
    assert_eq!(recorded, skip_lines("chrM"));
    for table in [
        "[expected]",
        "flags = [\"--lookup_ref\"",
        "flags = [\"--check_ref\"]",
    ] {
        let from = text.find(table).expect("the run's table");
        assert!(
            text[from..].lines().any(|line| line == "vep_rs = []"),
            "{table}"
        );
    }
}

/// The same holds over a FASTA that names the record of another length `M`:
/// both names are of its pair, and the skip line names the record as the
/// FASTA does.
#[test]
fn a_record_named_m_of_another_length_keeps_both_names_off_the_cache() {
    let dir = tempfile::tempdir().unwrap();
    let longer = two_bases_longer(dir.path(), "M");
    let (rows, warnings) = run(dir.path(), "plain", &longer, &[]);
    assert_eq!(rows, Vec::<String>::new());
    assert_eq!(warnings, skip_lines("M"));
}

/// Under `--dont_skip` such a record is kept as every record on a chromosome
/// no annotation source holds is: one row on no feature, beside its skip
/// line.
#[test]
fn dont_skip_writes_such_a_record_on_no_feature() {
    let dir = tempfile::tempdir().unwrap();
    let longer = two_bases_longer(dir.path(), "chrM");
    let (rows, warnings) = run(dir.path(), "kept", &longer, &["--dont_skip"]);
    let expected: Vec<String> = (1..=RECORDS.len()).map(row_on_no_feature).collect();
    assert_eq!(rows, expected);
    assert_eq!(warnings, skip_lines("chrM"));
}

/// Beside a `--custom` file that lists the chromosome, under the cache's name
/// or under the record's own, such a record is kept as any record on a
/// chromosome such a file lists is: one row on no feature and no warning,
/// with the file's record at its position under either name.
#[test]
fn a_custom_file_that_lists_the_chromosome_keeps_such_a_record_on_no_feature() {
    let dir = tempfile::tempdir().unwrap();
    let longer = two_bases_longer(dir.path(), "chrM");
    for listed in ["MT", "chrM"] {
        let file = common::indexed_vcf(
            dir.path(),
            &format!("custom_listing_{listed}_synthetic.vcf.gz"),
            &["##fileformat=VCFv4.2"],
            &format!("{listed}\t3349\tmade_up_record\tA\tC\t.\t.\t."),
        );
        let custom = format!(
            "file={},short_name=made_up,format=vcf,type=exact",
            file.display()
        );
        let (rows, warnings) = run(dir.path(), listed, &longer, &["--custom", &custom]);
        let expected: Vec<String> = (1..=RECORDS.len())
            .map(|line| {
                let row = row_on_no_feature(line);
                match RECORDS[line - 1].2 {
                    3349 => format!("{row};made_up=made_up_record"),
                    _ => row,
                }
            })
            .collect();
        assert_eq!(rows, expected, "{listed}");
        assert_eq!(warnings, "", "{listed}");
    }
}

/// A FASTA record of the assembly's length under either name of the pair is
/// the assembly's sequence, and no FASTA at all tells the two sequences
/// apart: the records are placed on the cache's `MT`, where Ensembl VEP 116.3
/// writes the coding transcript's row for the records written `chrM`.
#[test]
fn a_record_of_the_assembly_length_and_no_fasta_leave_the_records_on_mt() {
    let dir = tempfile::tempdir().unwrap();
    let reference = reference_bases();
    for name in ["chrM", "M", "MT"] {
        let file = fasta(
            dir.path(),
            &format!("{name}_of_the_assembly_length_synthetic.fa"),
            name,
            &reference,
        );
        let (rows, warnings) = run(dir.path(), name, &file, &[]);
        assert_eq!(rows.len(), 80, "{name}");
        for line in 1..=RECORDS.len() {
            assert_eq!(
                row_of(&rows, RECORDS[line - 1].0, CODING),
                row_on_mt(line),
                "{name}"
            );
        }
        assert_eq!(warnings, "", "{name}");
    }

    let out = dir.path().join("no_fasta.txt");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(records(dir.path()))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(crate_dir().join(CACHE))
        .args([
            "--offline",
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh37",
            "--format",
            "vcf",
            "--force_overwrite",
            "--quiet",
            "--no_stats",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success());
    let rows: Vec<String> = fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect();
    assert_eq!(rows.len(), 80);
    for line in 1..=RECORDS.len() {
        assert_eq!(row_of(&rows, RECORDS[line - 1].0, CODING), row_on_mt(line));
    }
}
