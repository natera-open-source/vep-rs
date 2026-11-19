// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A protein FASTA record is its transcript's peptide only when its length
//! is one the transcript's coding sequence translates to: its whole codons,
//! or one fewer, the stop codon's residue left out. A record of another
//! length is set aside and the transcript is built as one with no record:
//! its peptide is translated from the genome when one is given and absent
//! otherwise, the header's version goes unused, and the build's one warning
//! counts the record.

use std::path::Path;
use std::process::Command;

use serde_json::Value;

/// A GFF3 gene with one transcript `T1` on one exon whose CDS line spans
/// bases 4 to 21 of the genome below: six codons, the last a stop codon, so
/// the transcript's peptide has five residues.
const MODEL: &str = "##gff-version 3\n\
21\tx\tgene\t1\t30\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n\
21\tx\tmRNA\t1\t30\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n\
21\tx\texon\t1\t30\t.\t+\t.\tParent=transcript:T1\n\
21\tx\tCDS\t4\t21\t.\t+\t0\tParent=transcript:T1;protein_id=ENSP1\n";

/// The peptide the genome translates the CDS line to.
const TRANSLATED: &str = "MPKFG";

/// Writes a 60-base genome of sequence region `21` with its index and
/// returns its path.
fn write_genome(dir: &Path) -> std::path::PathBuf {
    let fa = dir.join("ref.fa");
    let seq = "GGGATGCCCAAATTTGGGTAGCCCTTTAAAGGGCCCATGAAACCCGGGTTTTAGCCCAAA";
    std::fs::write(&fa, format!(">21\n{seq}\n")).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), "21\t60\t4\t60\t61\n").unwrap();
    fa
}

/// Builds [`MODEL`] with a protein FASTA of one record for `T1`, its
/// sequence `residues` and its header's protein version 7, and the genome
/// when `genome` holds; returns the cached transcript and the build's log.
fn build(residues: &str, genome: bool) -> (Value, String) {
    build_model(MODEL, residues, genome)
}

/// [`build`] on the GFF3 text `model_text`.
fn build_model(model_text: &str, residues: &str, genome: bool) -> (Value, String) {
    let dir = tempfile::tempdir().unwrap();
    let model = dir.path().join("model.gff3");
    std::fs::write(&model, model_text).unwrap();
    let peptides = dir.path().join("peptides.fa");
    std::fs::write(
        &peptides,
        format!(">ENSP1.7 pep transcript:T1\n{residues}\n"),
    )
    .unwrap();
    let out = dir.path().join("cache");
    let mut command = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"));
    command
        .args(["--release", "116", "--assembly", "GRCh38", "--gff3"])
        .arg(&model)
        .arg("--fasta")
        .arg(&peptides)
        .arg("--output-dir")
        .arg(&out);
    if genome {
        command.arg("--genome-fasta").arg(write_genome(dir.path()));
    }
    let output = command.output().expect("vep-cache-builder runs");
    let log = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{log}");

    let shard = std::fs::read_dir(out.join("transcripts/21"))
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    let mut records: Vec<Value> =
        serde_json::from_str(&std::fs::read_to_string(shard).unwrap()).unwrap();
    assert_eq!(records.len(), 1);
    (records.remove(0), log)
}

/// The peptide a cached transcript carries, `None` when it carries none.
fn peptide(transcript: &Value) -> Option<&str> {
    transcript["variation_effect_feature_cache"]["peptide"].as_str()
}

/// The opening of the line a build logs for one record set aside and no
/// transcript without a record.
const ONE_SET_ASIDE: &str = "0 of 1 translated transcripts have no record in the protein FASTA, \
     and 1 records were set aside for their length";

/// A record of five residues, one fewer than the six codons, and one of six
/// are the transcript's peptide with or without a genome, whatever the
/// genome translates to: the cache carries the record's text and the
/// header's version, and the build logs no line about records.
#[test]
fn a_record_of_the_coding_sequences_length_is_the_peptide() {
    for record in ["MPKFU", "MPKFGX"] {
        for genome in [true, false] {
            let (transcript, log) = build(record, genome);
            assert_eq!(peptide(&transcript), Some(record), "genome {genome}");
            assert_eq!(transcript["translation"]["version"], 7, "{record}");
            assert!(
                !log.contains("have no record in the protein FASTA"),
                "{log}"
            );
        }
    }
}

/// With a genome, a record of three residues or of nine is set aside: the
/// peptide is the genome's translation of the CDS line, the translation
/// takes no version from the record's header, and the warning counts one
/// record set aside.
#[test]
fn a_record_of_another_length_leaves_the_genomes_translation() {
    for record in ["MPK", "MPKFGAAAA"] {
        let (transcript, log) = build(record, true);
        assert_eq!(peptide(&transcript), Some(TRANSLATED), "{record}");
        assert!(transcript["translation"]["version"].is_null(), "{record}");
        assert!(log.contains(ONE_SET_ASIDE), "{log}");
    }
}

/// Without a genome, a record of three residues or of nine is set aside and
/// the transcript carries no peptide, as one with no record: no version
/// from the header, and the same warning.
#[test]
fn a_record_of_another_length_leaves_no_peptide_without_a_genome() {
    for record in ["MPK", "MPKFGAAAA"] {
        let (transcript, log) = build(record, false);
        assert_eq!(peptide(&transcript), None, "{record}");
        assert!(transcript["translation"]["version"].is_null(), "{record}");
        assert!(log.contains(ONE_SET_ASIDE), "{log}");
    }
}

/// The coding sequence a record is held to is led by the start exon's phase
/// padding. A transcript whose CDS line opens its exon with phase column 1
/// has ten coding bases behind two bases of padding, four codons: a record
/// of four residues or of three is its peptide, and one of two, which ten
/// bases alone would admit, is set aside.
#[test]
fn the_start_exons_phase_padding_counts_toward_the_coding_sequence() {
    let padded = "##gff-version 3\n\
        21\tx\tgene\t4\t30\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n\
        21\tx\tmRNA\t4\t30\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n\
        21\tx\texon\t4\t30\t.\t+\t.\tParent=transcript:T1\n\
        21\tx\tCDS\t4\t13\t.\t+\t1\tParent=transcript:T1;protein_id=ENSP1\n";
    for record in ["MPKF", "MPK"] {
        let (transcript, log) = build_model(padded, record, false);
        assert_eq!(peptide(&transcript), Some(record));
        assert!(!log.contains("set aside for their length"), "{log}");
    }
    let (transcript, log) = build_model(padded, "MP", false);
    assert_eq!(peptide(&transcript), None);
    assert!(log.contains(ONE_SET_ASIDE), "{log}");
}
