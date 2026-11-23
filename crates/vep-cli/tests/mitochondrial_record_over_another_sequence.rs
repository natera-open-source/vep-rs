// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The documented difference `fasta_synonym_reads_another_mitochondrial_sequence`,
//! through the `vep` binary: a GRCh37 run whose reference FASTA holds its only
//! mitochondrial record under the name `chrM`, at another length than the
//! assembly's 16,569 bases, over a cache whose synonym table links `MT` to
//! `chrM`. The records are the reproduction's own,
//! `scripts/divergences/reproductions/fasta_synonym_reads_another_mitochondrial_sequence.toml`,
//! and the FASTA is built as its prepare steps build it.
//!
//! Ensembl VEP's reader serves a slice the FASTA lacks from the first synonym
//! of its name the FASTA holds (`Variation/Utils/FastaSequence.pm:378-397`),
//! and under a reference check names a record's slice over the FASTA's
//! records, a synonym first (`BaseVEP.pm:413` and `:708-716`); neither reads a
//! length. It therefore takes the reference bases of a record the cache places
//! on `MT` from the record named `chrM`, at the positions of another sequence.
//! A record written `chrM` or `M` is a record of that other sequence and is
//! skipped, the documented difference
//! `mitochondrial_record_of_another_sequence_placed_on_mt`, whose own test
//! file is `mitochondrial_record_named_for_another_sequence.rs`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

const REPRODUCTION: &str =
    "scripts/divergences/reproductions/fasta_synonym_reads_another_mitochondrial_sequence.toml";
const CACHE: &str = "tests/fixtures/bayesdel-mitochondrial-contig/json_cache";
const REFERENCE: &str = "tests/fixtures/mitochondrial-frame-check/reference_synthetic.fa";

/// The protein-coding transcript whose codon the substitutions change, and
/// the non-coding one the deletion and the insertion lie in.
const CODING: &str = "ENST00000361390";
const NON_CODING: &str = "ENST00000387347";

/// The substitution whose REF is the reference's base at its own position,
/// where the longer file holds another base, under each name the
/// reproduction writes it with that names the assembly's sequence: the
/// record's identifier and its chromosome as written. The cache places all
/// three on `MT`.
const WRITTEN: [(&str, &str); 3] = [
    ("other_base_two_back", "MT"),
    ("other_base_two_back_named_chrmt", "chrMT"),
    ("other_base_two_back_named_by_accession", "NC_012920.1"),
];

/// The row Ensembl VEP 116.3 writes for `record`, written on `chromosome`, on
/// the coding transcript under `--hgvs --hgvsg` over the reference of the
/// assembly's length: both notations state the record's REF, the first base
/// of the reference codon.
fn row_over_the_reference(record: &str, chromosome: &str) -> String {
    format!(
        "{record}\t{chromosome}:3349\tC\tENSG00000198888\tENST00000361390\tTranscript\t\
         missense_variant\t43\t43\t15\tI/L\tAtc/Ctc\t-\tIMPACT=MODERATE;STRAND=1;\
         HGVSc=ENST00000361390.2:c.43A>C;HGVSp=ENSP00000354687.2:p.Ile15Leu;\
         HGVSg={chromosome}:m.3349A>C"
    )
}

/// That row less the two notations that state a reference base: the row of a
/// run whose FASTA holds no sequence for `MT`. Ensembl VEP 116.3 writes
/// `c.43G>C` and `m.3349G>C` there under each of the three names and under
/// `chrM`, the base the longer file holds at position 3,349.
fn row_with_no_reference_base(record: &str, chromosome: &str) -> String {
    format!(
        "{record}\t{chromosome}:3349\tC\tENSG00000198888\tENST00000361390\tTranscript\t\
         missense_variant\t43\t43\t15\tI/L\tAtc/Ctc\t-\tIMPACT=MODERATE;STRAND=1;\
         HGVSp=ENSP00000354687.2:p.Ile15Leu"
    )
}

/// The deletion's and the insertion's rows on the non-coding transcript as
/// Ensembl VEP 116.3 writes them over the reference of the assembly's length.
/// Over the longer file it shifts the deletion two positions on
/// (`n.1286del`, `MT:m.2956del`, `HGVS_OFFSET=2`) and names the insertion a
/// duplication (`n.1289dup`, `MT:m.2959dup`), each from the other sequence's
/// bases.
const DELETION: &str = "deletion_shifted_over_the_longer_file";
const DELETION_OVER_THE_REFERENCE: &str = "deletion_shifted_over_the_longer_file\tMT:2954\t-\t\
    ENSG00000210082\tENST00000387347\tTranscript\tnon_coding_transcript_exon_variant\t1284\t-\t-\t\
    -\t-\t-\tIMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000387347.2:n.1284del;HGVSg=MT:m.2954del";
const INSERTION: &str = "insertion_duplicated_over_the_longer_file";
const INSERTION_OVER_THE_REFERENCE: &str = "insertion_duplicated_over_the_longer_file\t\
    MT:2959-2960\tA\tENSG00000210082\tENST00000387347\tTranscript\t\
    non_coding_transcript_exon_variant\t1289-1290\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;\
    HGVSc=ENST00000387347.2:n.1289_1290insA;HGVSg=MT:m.2959_2960insA";

/// The skip lines of the two records written `chrM` and `M`, which every run
/// over the longer file writes: each is a record of the FASTA's sequence,
/// which the cache does not hold.
const NAMED_CHRM: &str = "WARNING: line 3 skipped (chrM 3349 other_base_two_back_named_chrm A C . . .): Chromosome chrM is the FASTA record chrM of 16571 bases, another sequence than the cache's MT of 16569 bases, and is read as no chromosome of the cache\n";
const NAMED_M: &str = "WARNING: line 7 skipped (M 3349 other_base_two_back_named_m A C . . .): Chromosome M is the FASTA record chrM of 16571 bases, another sequence than the cache's MT of 16569 bases, and is read as no chromosome of the cache\n";

/// The warnings file of the run under `--check_ref`. A record written `MT`
/// or by the accession is skipped as on a contig the FASTA lacks, with the
/// sub-slice warning and then the mismatch warning (`Parser.pm:607-612` and
/// `:624-632`), as is one written `chrMT`, in Ensembl VEP too. A record
/// written `chrM` or `M` is skipped at the chromosome test, before the
/// check, where Ensembl VEP compares it with the FASTA record it names.
const CHECK_REF_WARNINGS: &str = "\
WARNING: line 1 skipped (MT 3349 other_base_two_back A C . . .): Could not fetch sub-slice from MT:3349-3349(1)
WARNING: line 1 skipped (MT 3349 other_base_two_back A C . . .): Specified reference allele A does not match Ensembl reference allele
WARNING: line 2 skipped (MT 3350 same_base_two_back T A . . .): Could not fetch sub-slice from MT:3350-3350(1)
WARNING: line 2 skipped (MT 3350 same_base_two_back T A . . .): Specified reference allele T does not match Ensembl reference allele
WARNING: line 3 skipped (chrM 3349 other_base_two_back_named_chrm A C . . .): Chromosome chrM is the FASTA record chrM of 16571 bases, another sequence than the cache's MT of 16569 bases, and is read as no chromosome of the cache
WARNING: line 4 skipped (chrMT 3349 other_base_two_back_named_chrmt A C...): Could not fetch sub-slice from chrMT:3349-3349(1)
WARNING: line 4 skipped (chrMT 3349 other_base_two_back_named_chrmt A C...): Specified reference allele A does not match Ensembl reference allele
WARNING: line 5 skipped (NC_012920.1 3349 other_base_two_back_named_by_...): Could not fetch sub-slice from NC_012920.1:3349-3349(1)
WARNING: line 5 skipped (NC_012920.1 3349 other_base_two_back_named_by_...): Specified reference allele A does not match Ensembl reference allele
WARNING: line 6 skipped (NC_012920.1 3350 same_base_two_back_named_by_a...): Could not fetch sub-slice from NC_012920.1:3350-3350(1)
WARNING: line 6 skipped (NC_012920.1 3350 same_base_two_back_named_by_a...): Specified reference allele T does not match Ensembl reference allele
WARNING: line 7 skipped (M 3349 other_base_two_back_named_m A C . . .): Chromosome M is the FASTA record chrM of 16571 bases, another sequence than the cache's MT of 16569 bases, and is read as no chromosome of the cache
WARNING: line 8 skipped (MT 2953 deletion_shifted_over_the_longer_file ...): Could not fetch sub-slice from MT:2954-2954(1)
WARNING: line 8 skipped (MT 2953 deletion_shifted_over_the_longer_file ...): Specified reference allele T does not match Ensembl reference allele
";

/// The warnings file of the run under `--lookup_ref` (`Parser.pm:638-643`):
/// the records written `chrM` and `M` are skipped at the chromosome test.
const LOOKUP_REF_WARNINGS: &str = "\
WARNING: line 1 skipped (MT 3349 other_base_two_back A C . . .): Could not fetch reference allele for MT:3349-3349(1)
WARNING: line 2 skipped (MT 3350 same_base_two_back T A . . .): Could not fetch reference allele for MT:3350-3350(1)
WARNING: line 3 skipped (chrM 3349 other_base_two_back_named_chrm A C . . .): Chromosome chrM is the FASTA record chrM of 16571 bases, another sequence than the cache's MT of 16569 bases, and is read as no chromosome of the cache
WARNING: line 4 skipped (chrMT 3349 other_base_two_back_named_chrmt A C...): Could not fetch reference allele for chrMT:3349-3349(1)
WARNING: line 5 skipped (NC_012920.1 3349 other_base_two_back_named_by_...): Could not fetch reference allele for NC_012920.1:3349-3349(1)
WARNING: line 6 skipped (NC_012920.1 3350 same_base_two_back_named_by_a...): Could not fetch reference allele for NC_012920.1:3350-3350(1)
WARNING: line 7 skipped (M 3349 other_base_two_back_named_m A C . . .): Chromosome M is the FASTA record chrM of 16571 bases, another sequence than the cache's MT of 16569 bases, and is read as no chromosome of the cache
WARNING: line 8 skipped (MT 2953 deletion_shifted_over_the_longer_file ...): Could not fetch reference allele for MT:2954-2954(1)
";

/// The warnings files Ensembl VEP 116.3 writes for the same records under a
/// `--synonyms` table that links no mitochondrial name. A bare `M` is then
/// renamed `MT` (`Parser.pm:794-796`) and read under that name, which reaches
/// no record, and the names the cache does not hold are not found.
const UNLINKED_CHECK_REF_WARNINGS: &str = "\
WARNING: line 1 skipped (MT 3349 other_base_two_back A C . . .): Could not fetch sub-slice from MT:3349-3349(1)
WARNING: line 1 skipped (MT 3349 other_base_two_back A C . . .): Specified reference allele A does not match Ensembl reference allele
WARNING: line 2 skipped (MT 3350 same_base_two_back T A . . .): Could not fetch sub-slice from MT:3350-3350(1)
WARNING: line 2 skipped (MT 3350 same_base_two_back T A . . .): Specified reference allele T does not match Ensembl reference allele
WARNING: line 3 skipped (chrM 3349 other_base_two_back_named_chrm A C . . .): Chromosome chrM not found in annotation sources or synonyms; chromosome chrM does not overlap any features
WARNING: line 4 skipped (chrMT 3349 other_base_two_back_named_chrmt A C...): Could not fetch sub-slice from chrMT:3349-3349(1)
WARNING: line 4 skipped (chrMT 3349 other_base_two_back_named_chrmt A C...): Specified reference allele A does not match Ensembl reference allele
WARNING: line 5 skipped (NC_012920.1 3349 other_base_two_back_named_by_...): Chromosome NC_012920.1 not found in annotation sources or synonyms; chromosome NC_012920.1 does not overlap any features
WARNING: line 6 skipped (NC_012920.1 3350 same_base_two_back_named_by_a...): Chromosome NC_012920.1 not found in annotation sources or synonyms; chromosome NC_012920.1 does not overlap any features
WARNING: line 7 skipped (M 3349 other_base_two_back_named_m A C . . .): Could not fetch sub-slice from MT:3349-3349(1)
WARNING: line 7 skipped (M 3349 other_base_two_back_named_m A C . . .): Specified reference allele A does not match Ensembl reference allele
WARNING: line 8 skipped (MT 2953 deletion_shifted_over_the_longer_file ...): Could not fetch sub-slice from MT:2954-2954(1)
WARNING: line 8 skipped (MT 2953 deletion_shifted_over_the_longer_file ...): Specified reference allele T does not match Ensembl reference allele
";
const UNLINKED_LOOKUP_REF_WARNINGS: &str = "\
WARNING: line 1 skipped (MT 3349 other_base_two_back A C . . .): Could not fetch reference allele for MT:3349-3349(1)
WARNING: line 2 skipped (MT 3350 same_base_two_back T A . . .): Could not fetch reference allele for MT:3350-3350(1)
WARNING: line 3 skipped (chrM 3349 other_base_two_back_named_chrm A C . . .): Chromosome chrM not found in annotation sources or synonyms; chromosome chrM does not overlap any features
WARNING: line 4 skipped (chrMT 3349 other_base_two_back_named_chrmt A C...): Could not fetch reference allele for chrMT:3349-3349(1)
WARNING: line 5 skipped (NC_012920.1 3349 other_base_two_back_named_by_...): Chromosome NC_012920.1 not found in annotation sources or synonyms; chromosome NC_012920.1 does not overlap any features
WARNING: line 6 skipped (NC_012920.1 3350 same_base_two_back_named_by_a...): Chromosome NC_012920.1 not found in annotation sources or synonyms; chromosome NC_012920.1 does not overlap any features
WARNING: line 7 skipped (M 3349 other_base_two_back_named_m A C . . .): Could not fetch reference allele for MT:3349-3349(1)
WARNING: line 8 skipped (MT 2953 deletion_shifted_over_the_longer_file ...): Could not fetch reference allele for MT:2954-2954(1)
";

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

/// The reproduction's input records as a VCF file in `dir`.
fn records(dir: &Path) -> PathBuf {
    let text = fs::read_to_string(crate_dir().join("../..").join(REPRODUCTION)).unwrap();
    let lines = block(&text, "text = \"\"\"", "\"\"\"");
    assert_eq!(lines.len(), 9, "the reproduction's records");
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

/// The file the reproduction prepares: the reference with its first two bases
/// put in front once more, named `chrM`, so its base at a position is the
/// reference's base two positions back.
fn chrm_two_bases_longer(dir: &Path) -> PathBuf {
    let reference = reference_bases();
    assert_eq!(reference.len(), 16_569);
    let longer = format!("{}{reference}", &reference[..2]);
    fasta(dir, "chrM_two_bases_longer_synthetic.fa", "chrM", &longer)
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

/// The one row of `rows` for `record` on `transcript`.
fn row_of<'a>(rows: &'a [String], record: &str, transcript: &str) -> &'a str {
    let mut found = rows.iter().filter(|row| {
        let cells: Vec<&str> = row.split('\t').collect();
        cells[0] == record && cells[4] == transcript
    });
    let row = found.next().expect("a row of the record on the transcript");
    assert!(
        found.next().is_none(),
        "one row of {record} on {transcript}"
    );
    row
}

/// The number of rows `rows` holds for `record`.
fn rows_of(rows: &[String], record: &str) -> usize {
    rows.iter()
        .filter(|row| row.split('\t').next() == Some(record))
        .count()
}

/// `row` in the columns the reproduction compares, joined as its recorded
/// values are: Location, Allele, Feature, Consequence, Codons and the two
/// notations, `-` for a notation the row lacks.
fn compared(row: &str) -> String {
    let cells: Vec<&str> = row.split('\t').collect();
    let item = |key: &str| {
        cells[13]
            .split(';')
            .find_map(|item| item.strip_prefix(key)?.strip_prefix('='))
            .unwrap_or("-")
    };
    [
        cells[1],
        cells[2],
        cells[4],
        cells[6],
        cells[11],
        item("HGVSc"),
        item("HGVSg"),
    ]
    .join(" | ")
}

/// `row` less the items of its last column that open with one of `keys`.
fn less(row: &str, keys: &[&str]) -> String {
    row.split(';')
        .filter(|item| !keys.iter().any(|key| item.starts_with(key)))
        .collect::<Vec<&str>>()
        .join(";")
}

/// A record the cache places on `MT` takes no reference base from a record
/// named `chrM` of another length, though the cache's table links the two
/// names and whichever of the assembly's names the input writes. The two
/// records written `chrM` and `M` are skipped. Under `--hgvs --hgvsg` a
/// substitution's row is the row over the reference less the two notations,
/// a deletion's likewise, and an insertion keeps the HGVSc it has over the
/// reference, unshifted; no row carries an HGVSg, and the compared row is the
/// one the reproduction records for vep-rs.
#[test]
fn a_record_on_mt_takes_no_reference_base_from_a_chrm_of_another_length() {
    let dir = tempfile::tempdir().unwrap();
    let longer = chrm_two_bases_longer(dir.path());
    let (rows, warnings) = run(dir.path(), "notations", &longer, &["--hgvs", "--hgvsg"]);

    assert_eq!(rows.len(), 138);
    assert_eq!(rows_of(&rows, "other_base_two_back_named_chrm"), 0);
    assert_eq!(rows_of(&rows, "other_base_two_back_named_m"), 0);
    for (record, chromosome) in WRITTEN {
        let expected = row_with_no_reference_base(record, chromosome);
        assert_eq!(row_of(&rows, record, CODING), expected, "{record}");
        assert_eq!(
            less(
                &row_over_the_reference(record, chromosome),
                &["HGVSc=", "HGVSg="]
            ),
            expected,
            "{record}"
        );
    }
    assert_eq!(
        row_of(&rows, DELETION, NON_CODING),
        less(DELETION_OVER_THE_REFERENCE, &["HGVSc=", "HGVSg="])
    );
    assert_eq!(
        row_of(&rows, INSERTION, NON_CODING),
        less(INSERTION_OVER_THE_REFERENCE, &["HGVSg="])
    );
    for row in &rows {
        assert!(!row.contains("HGVSg="), "{row}");
        assert!(
            !row.contains("HGVSc=") || row.starts_with(INSERTION),
            "{row}"
        );
    }
    assert_eq!(warnings, format!("{NAMED_CHRM}{NAMED_M}"));

    let text = fs::read_to_string(crate_dir().join("../..").join(REPRODUCTION)).unwrap();
    let recorded = text
        .lines()
        .find_map(|line| line.strip_prefix("vep_rs = [\"")?.strip_suffix("\"]"))
        .expect("the reproduction's recorded vep-rs row");
    assert_eq!(compared(row_of(&rows, WRITTEN[0].0, CODING)), recorded);
}

/// Under the two reference checks a record written `MT`, `chrMT` or by the
/// accession is skipped as on a contig the FASTA lacks, whether or not the
/// base the other sequence holds at its position equals its REF. A record
/// written `chrM` or `M` is skipped before either check, as a record of the
/// sequence the FASTA holds under that name. An insertion, which names no
/// reference base, passes both checks.
#[test]
fn the_reference_checks_skip_a_record_on_mt_over_a_chrm_of_another_length() {
    let dir = tempfile::tempdir().unwrap();
    let longer = chrm_two_bases_longer(dir.path());

    let (rows, warnings) = run(dir.path(), "check_ref", &longer, &["--check_ref"]);
    assert_eq!(rows.len(), 19);
    assert_eq!(rows_of(&rows, INSERTION), 19);
    assert_eq!(warnings, CHECK_REF_WARNINGS);

    let (rows, warnings) = run(
        dir.path(),
        "lookup_ref",
        &longer,
        &["--lookup_ref", "--show_ref_allele"],
    );
    assert_eq!(rows.len(), 19);
    assert_eq!(rows_of(&rows, INSERTION), 19);
    assert_eq!(warnings, LOOKUP_REF_WARNINGS);
}

/// A record named `chrM` of the assembly's length is the assembly's sequence
/// under another name and is read for `MT` through the table: the run writes
/// the rows of the run over the reference itself, each record's with the
/// notations Ensembl VEP 116.3 writes over the reference.
#[test]
fn a_chrm_of_the_assembly_length_is_read_for_a_record_on_mt() {
    let dir = tempfile::tempdir().unwrap();
    let renamed = fasta(
        dir.path(),
        "chrM_of_the_assembly_length_synthetic.fa",
        "chrM",
        &reference_bases(),
    );
    let flags = ["--hgvs", "--hgvsg"];
    let (rows, warnings) = run(dir.path(), "renamed", &renamed, &flags);
    let (over_the_reference, _) = run(
        dir.path(),
        "reference",
        &crate_dir().join(REFERENCE),
        &flags,
    );

    assert_eq!(rows.len(), 178);
    assert_eq!(rows, over_the_reference);
    for (record, chromosome) in WRITTEN {
        assert_eq!(
            row_of(&rows, record, CODING),
            row_over_the_reference(record, chromosome),
            "{record}"
        );
    }
    assert_eq!(
        row_of(&rows, "other_base_two_back_named_chrm", CODING),
        row_over_the_reference("other_base_two_back_named_chrm", "chrM")
    );
    assert_eq!(rows_of(&rows, "other_base_two_back_named_m"), 20);
    assert_eq!(
        row_of(&rows, DELETION, NON_CODING),
        DELETION_OVER_THE_REFERENCE
    );
    assert_eq!(
        row_of(&rows, INSERTION, NON_CODING),
        INSERTION_OVER_THE_REFERENCE
    );
    assert_eq!(warnings, "");
}

/// Under a table that links no mitochondrial name a bare `M` is renamed `MT`,
/// and its reference is read under that name, never under the name before
/// the rename: over the longer file it is skipped as every record on `MT`
/// is, with the warnings Ensembl VEP 116.3 writes.
#[test]
fn a_bare_m_renamed_mt_is_read_under_the_name_it_was_given() {
    let dir = tempfile::tempdir().unwrap();
    let longer = chrm_two_bases_longer(dir.path());
    let table = dir.path().join("table_without_a_mitochondrial_line.txt");
    fs::write(&table, "21\tchr21\n").unwrap();
    let table = table.to_str().unwrap();

    let (rows, warnings) = run(
        dir.path(),
        "unlinked_check_ref",
        &longer,
        &["--check_ref", "--synonyms", table],
    );
    assert_eq!(rows.len(), 19);
    assert_eq!(rows_of(&rows, INSERTION), 19);
    assert_eq!(warnings, UNLINKED_CHECK_REF_WARNINGS);

    let (rows, warnings) = run(
        dir.path(),
        "unlinked_lookup_ref",
        &longer,
        &["--lookup_ref", "--show_ref_allele", "--synonyms", table],
    );
    assert_eq!(rows.len(), 19);
    assert_eq!(rows_of(&rows, INSERTION), 19);
    assert_eq!(warnings, UNLINKED_LOOKUP_REF_WARNINGS);
}
