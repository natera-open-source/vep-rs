// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A GRCh37 run whose reference FASTA holds two mitochondrial records, `chrM`
//! and `chrMT`, through the `vep` binary: the further runs of the documented
//! difference `fasta_synonym_reads_another_mitochondrial_sequence`. The
//! records are those runs' own, read from
//! `scripts/divergences/reproductions/fasta_synonym_reads_another_mitochondrial_sequence.toml`,
//! and each FASTA is built as the record's prepare steps build it, and once
//! more with its two records listed the other way round.
//!
//! vep-rs reads, for the cache's `MT`, the record of the assembly's 16,569
//! bases, `chrMT` before `chrM` where both have them, whichever the file
//! lists first. Ensembl VEP 116.3 reads `chrMT` for the notations, the name
//! with `chr` before it (`Variation/Utils/FastaSequence.pm:384`), and under a
//! reference check names the slice of a record written `MT` after `chrM`, the
//! synonym the cache's table gives, before it tries that spelling
//! (`BaseVEP.pm:413` and `:708-716`); neither compares a length.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

const REPRODUCTION: &str =
    "scripts/divergences/reproductions/fasta_synonym_reads_another_mitochondrial_sequence.toml";
const CACHE: &str = "tests/fixtures/bayesdel-mitochondrial-contig/json_cache";
const REFERENCE: &str = "tests/fixtures/mitochondrial-frame-check/reference_synthetic.fa";

/// The protein-coding transcript whose codon the substitution changes.
const CODING: &str = "ENST00000361390";

/// The number of rows either engine writes for the substitution on the
/// cache's `MT`: one for each transcript within reach.
const ROWS_ON_MT: usize = 20;

/// The cells every row of the substitution on the coding transcript opens
/// with, after the record's name: the row both engines write over the
/// reference of the assembly's length.
fn row(record: &str, chromosome: &str, extra: &str) -> String {
    format!(
        "{record}\t{chromosome}:3349\tC\tENSG00000198888\tENST00000361390\tTranscript\t\
         missense_variant\t43\t43\t15\tI/L\tAtc/Ctc\t-\t{extra}"
    )
}

/// The last column of that row under `--hgvs --hgvsg` where the notations
/// state the made-up reference's base, the HGVSg under the name as written.
fn notations(chromosome: &str) -> String {
    format!(
        "IMPACT=MODERATE;STRAND=1;HGVSc=ENST00000361390.2:c.43A>C;\
         HGVSp=ENSP00000354687.2:p.Ile15Leu;HGVSg={chromosome}:m.3349A>C"
    )
}

fn crate_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn reproduction() -> String {
    fs::read_to_string(crate_dir().join("../..").join(REPRODUCTION)).unwrap()
}

/// The input records of the further run at `place` of the reproduction's
/// `also` list, counted from 1, as a VCF file in `dir`: the lines of the
/// `text` string of that run's `[also.input]` table.
fn records(dir: &Path, place: usize) -> PathBuf {
    let text = reproduction();
    let run = text
        .split("\n[[also]]\n")
        .nth(place)
        .unwrap_or_else(|| panic!("the reproduction holds no further run {place}"));
    let lines: Vec<&str> = run
        .lines()
        .skip_while(|line| *line != "text = \"\"\"")
        .skip(1)
        .take_while(|line| *line != "\"\"\"")
        .collect();
    assert!(!lines.is_empty(), "the records of further run {place}");
    let vcf = dir.join(format!("records_of_run_{place}.vcf"));
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

/// The values the further run at `place` records for vep-rs.
fn recorded_for_vep_rs(place: usize) -> Vec<String> {
    let text = reproduction();
    let run = text.split("\n[[also]]\n").nth(place).unwrap();
    let run = run.split("\n[implementation]\n").next().unwrap();
    let line = run
        .lines()
        .skip_while(|line| *line != "[also.expected]")
        .find(|line| line.starts_with("vep_rs = "))
        .unwrap_or_else(|| panic!("further run {place} records no vep_rs value"));
    line.trim_start_matches("vep_rs = [")
        .trim_end_matches(']')
        .split("\", \"")
        .map(|value| value.trim_matches('"').to_string())
        .collect()
}

/// The made-up reference's bases, its one record's lines joined.
fn reference_bases() -> String {
    let bases: String = fs::read_to_string(crate_dir().join(REFERENCE))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect();
    assert_eq!(bases.len(), 16_569);
    bases
}

/// The reference with its first two bases put in front once more, so its
/// base at a position is the reference's base two positions back.
fn two_bases_longer() -> String {
    let reference = reference_bases();
    format!("{}{reference}", &reference[..2])
}

/// An indexed FASTA in `dir` of the two records `contigs`, in that order,
/// sixty bases to a line, with the index `samtools faidx` writes for it.
fn fasta(dir: &Path, file: &str, contigs: [(&str, &str); 2]) -> PathBuf {
    let (mut text, mut index) = (String::new(), String::new());
    for (name, bases) in contigs {
        text.push_str(&format!(">{name}\n"));
        index.push_str(&format!(
            "{name}\t{}\t{}\t60\t61\n",
            bases.len(),
            text.len()
        ));
        for line in bases.as_bytes().chunks(60) {
            text.push_str(std::str::from_utf8(line).unwrap());
            text.push('\n');
        }
    }
    let path = dir.join(file);
    fs::write(&path, text).unwrap();
    fs::write(dir.join(format!("{file}.fai")), index).unwrap();
    path
}

/// The file the reproduction prepares for its two runs under a reference
/// check, the shape of the UCSC hg19 file set that adds the assembly's
/// sequence: `chrM` two bases longer than the assembly's sequence and
/// `chrMT` the made-up reference, `chrM` listed first as that set lists it,
/// or `chrMT` first.
fn chrm_longer_beside_chrmt(dir: &Path, chrm_first: bool) -> PathBuf {
    let (longer, reference) = (two_bases_longer(), reference_bases());
    if chrm_first {
        let contigs = [("chrM", longer.as_str()), ("chrMT", reference.as_str())];
        fasta(
            dir,
            "chrM_two_bases_longer_beside_chrMT_synthetic.fa",
            contigs,
        )
    } else {
        let contigs = [("chrMT", reference.as_str()), ("chrM", longer.as_str())];
        fasta(
            dir,
            "chrMT_before_a_chrM_two_bases_longer_synthetic.fa",
            contigs,
        )
    }
}

/// The file the reproduction prepares for its run under the notation flags:
/// `chrM` the made-up reference and `chrMT` two bases longer, `chrM` listed
/// first as the reproduction lists it, or `chrMT` first.
fn chrmt_longer_beside_chrm(dir: &Path, chrm_first: bool) -> PathBuf {
    let (longer, reference) = (two_bases_longer(), reference_bases());
    if chrm_first {
        let contigs = [("chrM", reference.as_str()), ("chrMT", longer.as_str())];
        fasta(
            dir,
            "chrMT_two_bases_longer_beside_chrM_synthetic.fa",
            contigs,
        )
    } else {
        let contigs = [("chrMT", longer.as_str()), ("chrM", reference.as_str())];
        fasta(
            dir,
            "chrMT_two_bases_longer_before_chrM_synthetic.fa",
            contigs,
        )
    }
}

/// The rows and the warnings file of a `vep` run over `records` and the
/// reproduction's cache with `fasta` as its reference.
fn run(
    dir: &Path,
    label: &str,
    records: &Path,
    fasta: &Path,
    flags: &[&str],
) -> (Vec<String>, String) {
    let out = dir.join(format!("{label}.txt"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(records)
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

/// The one row of `rows` for `record` on the coding transcript.
fn coding_row<'a>(rows: &'a [String], record: &str) -> &'a str {
    let mut found = rows.iter().filter(|row| {
        let cells: Vec<&str> = row.split('\t').collect();
        cells[0] == record && cells[4] == CODING
    });
    let row = found.next().expect("a row of the record on the transcript");
    assert!(found.next().is_none(), "one row of {record} on {CODING}");
    row
}

/// The number of rows `rows` holds for `record`.
fn rows_of(rows: &[String], record: &str) -> usize {
    rows.iter()
        .filter(|row| row.split('\t').next() == Some(record))
        .count()
}

/// `row` in the cells `columns` names, joined as a reproduction records its
/// values: a column of the default format by its name, any other name as a
/// key of the last column, `-` where the row lacks it.
fn compared(row: &str, columns: &[&str]) -> String {
    let cells: Vec<&str> = row.split('\t').collect();
    let cell = |column: &str| match column {
        "Uploaded_variation" => cells[0],
        "Location" => cells[1],
        "Allele" => cells[2],
        "Feature" => cells[4],
        "Consequence" => cells[6],
        "Codons" => cells[11],
        key => cells[13]
            .split(';')
            .find_map(|item| item.strip_prefix(key)?.strip_prefix('='))
            .unwrap_or("-"),
    };
    columns
        .iter()
        .map(|column| cell(column))
        .collect::<Vec<&str>>()
        .join(" | ")
}

/// Under `--check_ref`, over a FASTA that holds `chrM` two bases longer
/// beside `chrMT` of the assembly's length, a record the cache places on
/// `MT` is compared with the base of `chrMT`, whichever record the file
/// lists first: the three records that state it are annotated on the
/// transcripts of `MT` and the one that states the base of `chrM` is
/// refused with the mismatch warning alone. Ensembl VEP 116.3 compares a
/// record written `MT` or by the accession with the base of `chrM`: it
/// refuses those two, admits the one vep-rs refuses and writes each record
/// it admits as one `intergenic_variant` row.
#[test]
fn a_reference_check_holds_a_record_on_mt_to_chrmt_beside_a_chrm_of_another_length() {
    let dir = tempfile::tempdir().unwrap();
    let records = records(dir.path(), 1);
    for chrm_first in [true, false] {
        let fasta = chrm_longer_beside_chrmt(dir.path(), chrm_first);
        let label = format!("check_ref_{chrm_first}");
        let (rows, warnings) = run(dir.path(), &label, &records, &fasta, &["--check_ref"]);
        for (record, chromosome) in [
            ("stating_the_base_of_chrmt", "MT"),
            ("by_accession_stating_the_base_of_chrmt", "NC_012920.1"),
            ("named_chrmt_stating_its_base", "chrMT"),
        ] {
            assert_eq!(
                coding_row(&rows, record),
                row(record, chromosome, "IMPACT=MODERATE;STRAND=1"),
                "chrM first: {chrm_first}"
            );
            assert_eq!(rows_of(&rows, record), ROWS_ON_MT, "{record}");
        }
        assert_eq!(rows.len(), 3 * ROWS_ON_MT, "chrM first: {chrm_first}");
        assert_eq!(
            warnings,
            "WARNING: line 2 skipped (MT 3349 stating_the_base_of_chrm G C . . .): Specified \
             reference allele G does not match Ensembl reference allele\n",
            "chrM first: {chrm_first}"
        );
        let columns = ["Uploaded_variation", "Location", "Feature", "Consequence"];
        assert_eq!(
            vec![compared(
                coding_row(&rows, "stating_the_base_of_chrmt"),
                &columns
            )],
            recorded_for_vep_rs(1),
            "the value the reproduction records"
        );
    }
}

/// Under `--lookup_ref`, over the same FASTA, a record the cache places on
/// `MT` is given the base of `chrMT` as its reference and annotated on the
/// transcripts of `MT`, whichever record the file lists first. Ensembl VEP
/// 116.3 gives it the base of `chrM` and writes one `intergenic_variant` row.
#[test]
fn a_looked_up_reference_of_a_record_on_mt_is_the_base_of_chrmt() {
    let dir = tempfile::tempdir().unwrap();
    let records = records(dir.path(), 2);
    for chrm_first in [true, false] {
        let fasta = chrm_longer_beside_chrmt(dir.path(), chrm_first);
        let label = format!("lookup_ref_{chrm_first}");
        let flags = ["--lookup_ref", "--show_ref_allele"];
        let (rows, warnings) = run(dir.path(), &label, &records, &fasta, &flags);
        let named = [
            ("stating_the_base_of_chrmt", "MT"),
            ("by_accession_stating_the_base_of_chrmt", "NC_012920.1"),
        ];
        for (record, chromosome) in named {
            assert_eq!(
                coding_row(&rows, record),
                row(record, chromosome, "REF_ALLELE=A;IMPACT=MODERATE;STRAND=1"),
                "chrM first: {chrm_first}"
            );
            assert_eq!(rows_of(&rows, record), ROWS_ON_MT, "{record}");
        }
        assert_eq!(rows.len(), 2 * ROWS_ON_MT, "chrM first: {chrm_first}");
        assert_eq!(warnings, "", "chrM first: {chrm_first}");
        let columns = [
            "Uploaded_variation",
            "Location",
            "Feature",
            "Consequence",
            "REF_ALLELE",
        ];
        let written: Vec<String> = named
            .iter()
            .map(|(record, _)| compared(coding_row(&rows, record), &columns))
            .collect();
        assert_eq!(written, recorded_for_vep_rs(2));
    }
}

/// Under `--hgvs --hgvsg`, over the same FASTA, every record the cache
/// places on `MT` takes the reference base of its notations from `chrMT`,
/// whichever record the file lists first and whichever name the input
/// writes, and its row on the coding transcript is the row Ensembl VEP
/// 116.3 writes over that file, which reads `chrMT` as well. No row differs
/// between the engines there.
#[test]
fn the_notations_of_a_record_on_mt_read_chrmt_whichever_record_the_fasta_lists_first() {
    let dir = tempfile::tempdir().unwrap();
    let records = records(dir.path(), 1);
    for chrm_first in [true, false] {
        let fasta = chrm_longer_beside_chrmt(dir.path(), chrm_first);
        let label = format!("notations_{chrm_first}");
        let (rows, warnings) = run(dir.path(), &label, &records, &fasta, &["--hgvs", "--hgvsg"]);
        for (record, chromosome) in [
            ("stating_the_base_of_chrmt", "MT"),
            ("stating_the_base_of_chrm", "MT"),
            ("by_accession_stating_the_base_of_chrmt", "NC_012920.1"),
            ("named_chrmt_stating_its_base", "chrMT"),
        ] {
            assert_eq!(
                coding_row(&rows, record),
                row(record, chromosome, &notations(chromosome)),
                "chrM first: {chrm_first}"
            );
            assert_eq!(rows_of(&rows, record), ROWS_ON_MT, "{record}");
        }
        assert_eq!(warnings, "", "chrM first: {chrm_first}");
    }
}

/// Over a FASTA that holds `chrMT` two bases longer beside a `chrM` of the
/// assembly's length, a record on `MT` takes the reference base of its
/// notations from `chrM`, the record of 16,569 bases, whichever the file
/// lists first. Ensembl VEP 116.3 reads `chrMT` there and writes `c.43G>C`
/// and `MT:m.3349G>C`, the base that record holds at the position.
#[test]
fn beside_a_chrmt_of_another_length_mt_reads_the_chrm_of_the_assembly_length() {
    let dir = tempfile::tempdir().unwrap();
    let records = records(dir.path(), 3);
    let record = "stating_the_base_of_chrm_of_the_assembly_length";
    for chrm_first in [true, false] {
        let fasta = chrmt_longer_beside_chrm(dir.path(), chrm_first);
        let label = format!("chrmt_longer_{chrm_first}");
        let (rows, warnings) = run(dir.path(), &label, &records, &fasta, &["--hgvs", "--hgvsg"]);
        assert_eq!(
            coding_row(&rows, record),
            row(record, "MT", &notations("MT")),
            "chrM first: {chrm_first}"
        );
        assert_eq!(rows.len(), ROWS_ON_MT, "chrM first: {chrm_first}");
        assert_eq!(warnings, "", "chrM first: {chrm_first}");
        let columns = [
            "Location",
            "Allele",
            "Feature",
            "Consequence",
            "Codons",
            "HGVSc",
            "HGVSg",
        ];
        assert_eq!(
            vec![compared(coding_row(&rows, record), &columns)],
            recorded_for_vep_rs(3)
        );
    }
}
