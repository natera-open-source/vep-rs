// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows the `vep` binary writes under `--hgvsg` for the records of the
//! further run of the documented difference `contig_absent_from_fasta`: a
//! substitution, two deletions and two insertions on the mitochondrial
//! chromosome, annotated over a cache that holds it and a FASTA that does not.
//! The input is read from the entry's record under
//! `scripts/divergences/reproductions`, so each line is spelt in one file.
//! The last test holds the record's `code` row to the functions that hold the
//! entry's rows.

use std::path::{Path, PathBuf};
use std::process::{Output, Stdio};

mod common;

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The records of the further run, each with its Location and Allele cells and
/// the HGVSg both engines write for it over a FASTA that holds the chromosome.
const RECORDS: [(&str, &str, &str, &str); 5] = [
    ("mt_substitution", "MT:700", "C", "MT:m.700A>C"),
    ("mt_deletion", "MT:714", "-", "MT:m.714del"),
    ("mt_deletion_in_a_run", "MT:736", "-", "MT:m.737del"),
    ("mt_insertion", "MT:745-746", "C", "MT:m.745_746insC"),
    (
        "mt_insertion_of_the_base_before_it",
        "MT:755-756",
        "G",
        "MT:m.755dup",
    ),
];

/// The transcripts every one of the records gets a row of, each with its term.
const TRANSCRIPTS: [(&str, &str); 13] = [
    ("ENST00000361390", "upstream_gene_variant"),
    ("ENST00000361453", "upstream_gene_variant"),
    ("ENST00000386347", "upstream_gene_variant"),
    ("ENST00000387314", "downstream_gene_variant"),
    ("ENST00000387342", "upstream_gene_variant"),
    ("ENST00000387347", "upstream_gene_variant"),
    ("ENST00000387365", "upstream_gene_variant"),
    ("ENST00000387372", "downstream_gene_variant"),
    ("ENST00000387377", "upstream_gene_variant"),
    ("ENST00000387382", "upstream_gene_variant"),
    ("ENST00000387392", "downstream_gene_variant"),
    ("ENST00000387400", "downstream_gene_variant"),
    ("ENST00000389680", "non_coding_transcript_exon_variant"),
];

/// The notice a run under `--hgvsg` alone prints at its end for the chromosome.
const NO_SEQUENCE_FOR_MT: &str = "INFO: The FASTA holds no sequence for chromosome MT, which the \
    input has records on: a record there gets no HGVSg";

/// A file or directory of the repository, by its path from the root.
fn committed(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(path)
}

/// The corpus the entry's record names: its cache holds chromosome 21 and the
/// mitochondrial chromosome, and so does its FASTA.
fn corpus() -> PathBuf {
    committed("tests/golden/116/GRCh38-input-hgvs")
}

/// The entry's record, by its path from the repository's root.
const RECORD: &str = "scripts/divergences/reproductions/contig_absent_from_fasta.toml";

/// The input lines of the record's further run: the `text` of its
/// `[also.input]` table, a VCF without its header.
fn further_run_lines() -> String {
    let path = committed(RECORD);
    let record = std::fs::read_to_string(&path).unwrap();
    let lines: Vec<&str> = record
        .lines()
        .skip_while(|line| *line != "[also.input]")
        .skip_while(|line| *line != "text = \"\"\"")
        .skip(1)
        .take_while(|line| *line != "\"\"\"")
        .collect();
    let ids: Vec<&str> = lines
        .iter()
        .filter_map(|line| line.split('\t').nth(2))
        .collect();
    let expected: Vec<&str> = RECORDS.iter().map(|record| record.0).collect();
    assert_eq!(ids, expected, "{}", path.display());
    lines.iter().map(|line| format!("{line}\n")).collect()
}

/// A plain FASTA of one sequence of 120 bases named `21`, with its `.fai`: a
/// FASTA that holds no sequence for the mitochondrial chromosome.
fn fasta_without_the_chromosome(dir: &Path) -> PathBuf {
    let path = dir.join("only_21.fa");
    std::fs::write(
        &path,
        format!(">21\n{}\n{}\n", "ACGT".repeat(15), "ACGT".repeat(15)),
    )
    .unwrap();
    std::fs::write(format!("{}.fai", path.display()), "21\t120\t4\t60\t61\n").unwrap();
    path
}

/// One run of the binary over the further run's records.
struct Run {
    output: Output,
    /// The data rows of the output file.
    rows: Vec<String>,
    /// The warnings file, empty when the run wrote none.
    warnings: String,
}

/// Runs the binary under `--hgvsg` on the further run's records over the
/// corpus's cache and `fasta`, writing `<name>.txt` under `dir`; exit status 0
/// is required.
fn annotate(dir: &Path, name: &str, fasta: &Path) -> Run {
    let input = dir.join("further_run.vcf");
    std::fs::write(&input, format!("{VCF_HEADER}{}", further_run_lines())).unwrap();
    let out = dir.join(format!("{name}.txt"));
    let mut command = common::isolated_vep(dir);
    command
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus().join("json_cache"))
        .arg("--fasta")
        .arg(fasta)
        .args(["--species", "homo_sapiens", "--assembly", "GRCh38"])
        .args(["--format", "vcf", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--hgvsg", "--tab"])
        .args([
            "--fields",
            "Uploaded_variation,Location,Allele,Feature,Consequence,HGVSg",
        ]);
    let output = command.stdin(Stdio::null()).output().unwrap();
    assert_eq!(
        output.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let rows = std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect();
    let warnings =
        std::fs::read_to_string(format!("{}_warnings.txt", out.display())).unwrap_or_default();
    Run {
        output,
        rows,
        warnings,
    }
}

/// The 65 rows of the run, a record's thirteen in the order of `TRANSCRIPTS`,
/// with the HGVSg cell `hgvsg` gives for each record.
fn rows_with(hgvsg: impl Fn(&(&str, &str, &str, &str)) -> String) -> Vec<String> {
    RECORDS
        .iter()
        .flat_map(|record| {
            let (id, location, allele, _) = *record;
            let notation = hgvsg(record);
            TRANSCRIPTS.iter().map(move |(transcript, term)| {
                format!("{id}\t{location}\t{allele}\t{transcript}\t{term}\t{notation}")
            })
        })
        .collect()
}

/// Under `--hgvsg` no record on a chromosome the FASTA holds no sequence for
/// gets an HGVSg, whatever its kind: both ends of a record must be bases the
/// FASTA holds, so vep-rs writes no notation it would have to test against
/// bases it did not read. Ensembl VEP 116.3 reads `N` per base there
/// (`Utils/FastaSequence.pm:597-601`) and writes each record's notation from
/// its own span: `MT:m.700N>C`, `MT:m.714del`, `MT:m.736del`,
/// `MT:m.745_746insC` and `MT:m.755_756insG`, as the entry's record holds
/// them. Over the corpus reference, which holds the chromosome, the rows carry
/// the notation both engines write there: the substitution names the base, the
/// deletion in a run stands on the run's last base and the insertion of a copy
/// of the base before it is a duplication. No warning is written.
/// `docs/cache-setup.md` states the difference.
#[test]
fn no_kind_of_record_on_a_contig_the_fasta_lacks_gets_an_hgvsg() {
    let dir = tempfile::tempdir().unwrap();
    let lacking = annotate(
        dir.path(),
        "lacking",
        &fasta_without_the_chromosome(dir.path()),
    );
    assert_eq!(lacking.rows, rows_with(|_| "-".to_string()));
    assert_eq!(lacking.warnings, "");
    let held = annotate(dir.path(), "held", &corpus().join("reference.fa.gz"));
    assert_eq!(held.rows, rows_with(|record| record.3.to_string()));
    assert_eq!(held.warnings, "");
}

/// The run that leaves the notations out says so once at its end: the notice
/// names the chromosome and the notation the run asked for. A run over a
/// FASTA that holds the chromosome prints no such line.
#[test]
fn the_run_names_the_contig_whose_records_got_no_hgvsg() {
    let dir = tempfile::tempdir().unwrap();
    let notices = |run: &Run| -> Vec<String> {
        String::from_utf8_lossy(&run.output.stderr)
            .lines()
            .filter(|line| line.contains("The FASTA holds no sequence for"))
            .map(str::to_string)
            .collect()
    };
    let lacking = annotate(
        dir.path(),
        "lacking",
        &fasta_without_the_chromosome(dir.path()),
    );
    assert_eq!(notices(&lacking), [NO_SEQUENCE_FOR_MT]);
    let held = annotate(dir.path(), "held", &corpus().join("reference.fa.gz"));
    assert_eq!(notices(&held), Vec::<String>::new());
}

/// The items the `code` row of the entry's record names, in the record's
/// order: the strings of the `code` list of its `[implementation]` table.
fn code_row() -> Vec<String> {
    let path = committed(RECORD);
    let record = std::fs::read_to_string(&path).unwrap();
    let list = record
        .split_once("\n[implementation]\n")
        .and_then(|(_, table)| table.split_once("code = ["))
        .and_then(|(_, rest)| rest.split_once(']'))
        .unwrap_or_else(|| panic!("{}: no code row", path.display()))
        .0;
    list.split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect()
}

/// The record names the eight functions that hold what vep-rs writes on a
/// chromosome the FASTA holds no sequence for, where Ensembl VEP reads `N`
/// per base (`Utils/FastaSequence.pm:597-601`). `IndexedFasta::base` and
/// `IndexedFasta::sequence_in_case` read nothing there. Under `--hgvsg`,
/// `hgvs_genomic` writes no notation for a record whose ends are no bases of
/// the FASTA and `flank_sequence` reads no flank whose first base is none:
/// either keeps a record without an HGVSg where the other reads `N`. Under
/// `--hgvs`, `TranscriptSlice::substr` reads no reference bases and
/// `hgvs_variant_notation` writes no notation whose reference bases cannot be
/// read. Under `--lookup_ref` and `--check_ref`, `reference_sequence` reads
/// nothing and `SkipRules::validate` warns and skips a record that names a
/// reference base. With the part of any one of them given Ensembl VEP's
/// reading, a test the record names fails; for the two functions of
/// `hgvsg.rs`, with both.
///
/// `hgvs_left_out_notices` of `runner.rs` is no row of the list: it names the
/// chromosome on standard error at the run's end and writes no row of the
/// output, so a run with it silent writes the same rows.
/// `the_run_names_the_contig_whose_records_got_no_hgvsg` holds its line.
///
/// No FASTA read of `coding.rs` or `consequences.rs` of `vep-effects` holds a
/// row of the record. The tests the record names reach one of them on such a
/// chromosome, the 3' shift's (`shift_indel_3prime`), which stops at a base
/// that is not its allele's as it stops where it reads none.
#[test]
fn the_record_names_the_functions_that_hold_its_rows_on_a_contig_the_fasta_lacks() {
    assert_eq!(
        code_row(),
        [
            "crates/vep-cli/src/validate.rs::SkipRules::validate",
            "crates/vep-cli/src/validate.rs::reference_sequence",
            "crates/vep-effects/src/hgvs.rs::TranscriptSlice::substr",
            "crates/vep-effects/src/hgvs.rs::hgvs_variant_notation",
            "crates/vep-effects/src/hgvsg.rs::flank_sequence",
            "crates/vep-effects/src/hgvsg.rs::hgvs_genomic",
            "crates/vep-fasta/src/lib.rs::IndexedFasta::base",
            "crates/vep-fasta/src/lib.rs::IndexedFasta::sequence_in_case",
        ]
    );
}
