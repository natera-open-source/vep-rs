// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads a deletion that removes whole codons from the
//! first coding base of a transcript that has no 5' UTR.
//!
//! Two transcripts of one gene share a start codon at 21:32298972-32298974:
//! [`WITH_UTR`] holds 117 bases before it and [`WITHOUT_UTR`] opens on it.
//! A record that removes the three bases removes the start codon of both,
//! and the triplet left at the start of the coding sequence is no `ATG`, so
//! both rows are `start_lost` beside `inframe_deletion` (SO:0002012) with
//! IMPACT `HIGH` and, under `--hgvs`, the protein notation of a lost start.
//! Ensembl VEP 116.3 writes that row on the transcript with the 5' UTR and
//! `inframe_deletion` alone, `MODERATE`, with a `del` notation on the other:
//! its test of the edited sequence returns where the transcript has no 5' UTR
//! (VariationEffect.pm:926), the length test, which reads the start as
//! altered (:1069), is set aside for a deletion of whole codons (:870) and
//! the peptide test returns on the empty alternate peptide (:875). That is
//! the documented difference `start_codon_removed_without_five_prime_utr`.
//!
//! The records are made up on the chromosome 21 reference the corpus carries
//! and read on its release-116 GRCh38 cache.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// A corpus whose chromosome 21 cache slice holds both transcripts, with its
/// reference.
fn corpus() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-plugin-loftee-fixture")
}

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The base before the start codon, 21:32298971, and that base itself.
const POSITION: u64 = 32_298_971;
const SHARED_BASE: &str = "C";

/// The start codon and the codon after it, as the reference holds them.
const START_CODON: &str = "ATG";
const SECOND_CODON: &str = "GCC";

const GENE: &str = "ENSG00000170262";

/// A transcript whose coding sequence starts at its cDNA position 118.
const WITH_UTR: (&str, &str, &str) = ("ENST00000303645", "ENST00000303645.10", "ENSP00000306697.5");
/// A transcript whose coding sequence starts at its first base.
const WITHOUT_UTR: (&str, &str, &str) =
    ("ENST00000339944", "ENST00000339944.4", "ENSP00000343661.4");

/// A record that removes `removed` after the shared base.
fn vcf_line(id: &str, removed: &str) -> String {
    format!("21\t{POSITION}\t{id}\t{SHARED_BASE}{removed}\t{SHARED_BASE}\t.\t.\t.\n")
}

/// The data rows `vep --tab --hgvs` writes for `records` on the corpus.
fn annotate(dir: &Path, records: &str) -> Vec<Vec<String>> {
    let input = dir.join("records.vcf");
    std::fs::write(&input, format!("{VCF_HEADER}{records}")).unwrap();
    let out = dir.join("records.tab");
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus().join("json_cache"))
        .arg("--fasta")
        .arg(corpus().join("reference.fa.gz"))
        .args([
            "--tab",
            "--hgvs",
            "--no_stats",
            "--quiet",
            "--force_overwrite",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

/// The row of `transcript` for the record `id`, its cells joined by tabs.
fn written(rows: &[Vec<String>], id: &str, transcript: (&str, &str, &str)) -> String {
    rows.iter()
        .find(|row| row[0] == id && row[4] == transcript.0)
        .unwrap_or_else(|| panic!("{id}: no row of {}", transcript.0))
        .join("\t")
}

/// A whole row of the tab output under `--tab --hgvs`: the location, the
/// terms, the cDNA, CDS and protein positions, the residues and the codons,
/// the IMPACT, the two notations after their identifiers and the
/// HGVS_OFFSET.
fn row(
    id: &str,
    transcript: (&str, &str, &str),
    location: &str,
    consequence: &str,
    positions: (&str, &str, &str),
    alleles: (&str, &str),
    notations: (&str, &str, &str),
) -> String {
    let (cdna_position, cds_position, protein_position) = positions;
    let (amino_acids, codons) = alleles;
    let (coding, protein, offset) = notations;
    [
        id,
        location,
        "-",
        GENE,
        transcript.0,
        "Transcript",
        consequence,
        cdna_position,
        cds_position,
        protein_position,
        amino_acids,
        codons,
        "-",
        "HIGH",
        "-",
        "1",
        "-",
        &format!("{}:{coding}", transcript.1),
        &format!("{}:{protein}", transcript.2),
        offset,
    ]
    .join("\t")
}

/// The start codon removed whole: the row of the transcript without a 5' UTR
/// is the row of the transcript with one, cell for cell but the cDNA
/// position.
#[test]
fn the_start_codon_removed_whole_is_start_lost_on_a_transcript_without_a_five_prime_utr() {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(dir.path(), &vcf_line("start_codon_removed", START_CODON));
    for (transcript, cdna_position) in [(WITH_UTR, "118-120"), (WITHOUT_UTR, "1-3")] {
        assert_eq!(
            written(&rows, "start_codon_removed", transcript),
            row(
                "start_codon_removed",
                transcript,
                "21:32298972-32298974",
                "start_lost,inframe_deletion",
                (cdna_position, "1-3", "1"),
                ("M/-", "ATG/-"),
                ("c.1_3del", "p.Met1?", "-"),
            )
        );
    }
}

/// The first two codons removed whole: the same terms on both transcripts.
/// The notations are those of the removal one base 3', where the HGVS 3' rule
/// places it.
#[test]
fn two_leading_codons_removed_whole_are_start_lost_on_a_transcript_without_a_five_prime_utr() {
    let dir = tempfile::tempdir().unwrap();
    let removed = format!("{START_CODON}{SECOND_CODON}");
    let rows = annotate(dir.path(), &vcf_line("two_codons_removed", &removed));
    for (transcript, cdna_position) in [(WITH_UTR, "118-123"), (WITHOUT_UTR, "1-6")] {
        assert_eq!(
            written(&rows, "two_codons_removed", transcript),
            row(
                "two_codons_removed",
                transcript,
                "21:32298972-32298977",
                "start_lost,inframe_deletion",
                (cdna_position, "1-6", "1-2"),
                ("MA/-", "ATGGCC/-"),
                ("c.2_7del", "p.MetAla1_?2", "1"),
            )
        );
    }
}

/// One base of the start codon removed shifts the frame: `start_lost` on both
/// transcripts in both engines (VariationEffect.pm:870, :1069).
#[test]
fn one_base_of_the_start_codon_removed_is_start_lost_on_both_transcripts() {
    let dir = tempfile::tempdir().unwrap();
    let rows = annotate(
        dir.path(),
        &vcf_line("first_base_removed", &START_CODON[..1]),
    );
    for (transcript, cdna_position) in [(WITH_UTR, "118"), (WITHOUT_UTR, "1")] {
        assert_eq!(
            written(&rows, "first_base_removed", transcript),
            row(
                "first_base_removed",
                transcript,
                "21:32298972",
                "frameshift_variant,start_lost",
                (cdna_position, "1", "1"),
                ("M/X", "Atg/tg"),
                ("c.1del", "p.Met1?", "-"),
            )
        );
    }
}
