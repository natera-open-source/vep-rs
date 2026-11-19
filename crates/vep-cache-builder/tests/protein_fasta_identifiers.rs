// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A protein FASTA record goes to the transcript whose identifier its
//! header's `transcript:` token is, whole: an identifier with dots of its own
//! takes its record, transcripts that share a first component each take
//! their own, and a token that ends in a version goes to the model of that
//! version.
//!
//! `tests/fixtures/tair10_slice/` holds the GFF3 lines (`slice.gff3`) and the
//! protein FASTA records (`peptides.fa`) of two Arabidopsis thaliana genes,
//! seven transcripts, of Ensembl Plants release 63 (TAIR10). Their transcript
//! and protein identifiers end in a dot and a number that is the
//! identifier's own (`AT1G01020.2`), and the release writes no version.

use std::collections::BTreeMap;
use std::path::Path;
use std::process::Command;

use serde_json::Value;

/// Builds the GFF3 `model` with the protein FASTA `peptides` and no genome,
/// and returns every transcript of the cache by stable identifier, with the
/// builder's log.
fn build(model: &Path, peptides: &Path) -> (BTreeMap<String, Value>, String) {
    let out = tempfile::tempdir().unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "63", "--assembly", "TAIR10"])
        .args(["--species", "arabidopsis_thaliana", "--gff3"])
        .arg(model)
        .arg("--fasta")
        .arg(peptides)
        .arg("--output-dir")
        .arg(out.path())
        .output()
        .expect("vep-cache-builder runs");
    let log = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{log}");

    let mut transcripts = BTreeMap::new();
    for region in std::fs::read_dir(out.path().join("transcripts")).unwrap() {
        for shard in std::fs::read_dir(region.unwrap().path()).unwrap() {
            let text = std::fs::read_to_string(shard.unwrap().path()).unwrap();
            let records: Vec<Value> = serde_json::from_str(&text).unwrap();
            for record in records {
                transcripts.insert(record["stable_id"].as_str().unwrap().to_string(), record);
            }
        }
    }
    (transcripts, log)
}

/// The records of a protein FASTA, each sequence under the text of its
/// header's `transcript:` token.
fn records(peptides: &Path) -> BTreeMap<String, String> {
    let mut records: BTreeMap<String, String> = BTreeMap::new();
    let mut current = None;
    for line in std::fs::read_to_string(peptides).unwrap().lines() {
        if line.starts_with('>') {
            current = line
                .split_whitespace()
                .find_map(|token| token.strip_prefix("transcript:"))
                .map(str::to_string);
        } else if let Some(id) = &current {
            records.entry(id.clone()).or_default().push_str(line.trim());
        }
    }
    records
}

/// The peptide a cached transcript carries, `None` when it carries none.
fn peptide(transcript: &Value) -> Option<&str> {
    transcript["variation_effect_feature_cache"]["peptide"].as_str()
}

/// Every transcript of the Arabidopsis slice carries the peptide of the
/// record that names it, the six transcripts of `AT1G01020` each their own,
/// and the build counts no transcript without a record. Its translation has
/// no version, as Ensembl VEP's cache holds it: the number that ends the
/// header's `AT1G01020.2` is the protein identifier's own.
#[test]
fn a_build_from_dotted_identifiers_uses_its_peptides() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/tair10_slice");
    let peptides = fixture.join("peptides.fa");
    let (transcripts, log) = build(&fixture.join("slice.gff3"), &peptides);
    let records = records(&peptides);
    assert_eq!(records.len(), 7);
    assert_eq!(
        transcripts.keys().collect::<Vec<_>>(),
        records.keys().collect::<Vec<_>>()
    );
    for (id, transcript) in &transcripts {
        assert_eq!(peptide(transcript), Some(records[id].as_str()), "{id}");
        let translation = &transcript["translation"];
        assert_eq!(translation["stable_id"], id.as_str(), "{id}");
        assert!(translation["version"].is_null(), "{id}: {translation}");
    }
    assert!(
        log.contains("Loaded 7 protein sequences from FASTA"),
        "{log}"
    );
    assert!(
        !log.contains("have no record in the protein FASTA"),
        "{log}"
    );
}

/// A GFF3 gene `G1` whose transcript `T1.1` (and `T1.2` when `second`) the
/// transcript line gives `version` when one is named; each codes on one exon
/// for the protein of its own identifier.
fn model(version: Option<u32>, second: bool) -> String {
    let version_attribute = version.map_or(String::new(), |v| format!(";version={v}"));
    let mut text = String::from(
        "##gff-version 3\n\
         21\tx\tgene\t1000\t2000\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n",
    );
    let ids: &[&str] = if second { &["T1.1", "T1.2"] } else { &["T1.1"] };
    for id in ids {
        text += &format!(
            "21\tx\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:{id};Parent=gene:G1;\
             biotype=protein_coding;transcript_id={id}{version_attribute}\n\
             21\tx\texon\t1000\t2000\t.\t+\t.\tParent=transcript:{id}\n\
             21\tx\tCDS\t1100\t1108\t.\t+\t0\tParent=transcript:{id};protein_id=P-{id}\n"
        );
    }
    text
}

/// A cached transcript's identifier, its peptide and its translation's
/// version.
type Row = (String, Option<String>, Option<u64>);

/// Builds `model` with a protein FASTA of `fasta` and returns the [`Row`] of
/// each transcript, with the log.
fn built(model: &str, fasta: &str) -> (Vec<Row>, String) {
    let dir = tempfile::tempdir().unwrap();
    let model_path = dir.path().join("model.gff3");
    std::fs::write(&model_path, model).unwrap();
    let peptides = dir.path().join("peptides.fa");
    std::fs::write(&peptides, fasta).unwrap();
    let (transcripts, log) = build(&model_path, &peptides);
    let rows = transcripts
        .iter()
        .map(|(id, transcript)| {
            (
                id.clone(),
                peptide(transcript).map(str::to_string),
                transcript["translation"]["version"].as_u64(),
            )
        })
        .collect();
    (rows, log)
}

/// Two transcripts whose identifiers share their first component each take
/// the record that names them whole, in either order of the records; a
/// record that names the shared component alone is neither's.
#[test]
fn transcripts_sharing_a_first_component_take_their_own_records() {
    let first = ">P-T1.1 pep transcript:T1.1\nMK\n";
    let second = ">P-T1.2 pep transcript:T1.2\nMKV\n";
    let own = [
        ("T1.1".to_string(), Some("MK".to_string()), None),
        ("T1.2".to_string(), Some("MKV".to_string()), None),
    ];
    for fasta in [format!("{first}{second}"), format!("{second}{first}")] {
        let (rows, log) = built(&model(None, true), &fasta);
        assert_eq!(rows, own, "{log}");
    }

    let (rows, log) = built(&model(None, true), ">ENSP1 pep transcript:T1\nMKVL\n");
    assert_eq!(
        rows,
        [
            ("T1.1".to_string(), None, None),
            ("T1.2".to_string(), None, None)
        ]
    );
    assert!(
        log.contains("2 of 2 translated transcripts have no record in the protein FASTA"),
        "{log}"
    );
}

/// A token that ends in a version goes to the model of that version: the
/// transcript `T1.1` at version 3 takes the record of `T1.1.3`, with the
/// version its header's protein identifier ends in, and the record of
/// `T1.1` written without a version, and the build counts it with a record
/// either way; the record of `T1.1.4`, another version of the transcript, is
/// not its own and the build counts it without one.
#[test]
fn a_versioned_token_goes_to_the_model_of_that_version() {
    let versioned = model(Some(3), false);
    for fasta in [
        ">P-T1.1.2 pep transcript:T1.1.3\nMK\n",
        ">P-T1.1.2 pep transcript:T1.1\nMK\n",
    ] {
        let (rows, log) = built(&versioned, fasta);
        assert_eq!(
            rows,
            [("T1.1".to_string(), Some("MK".to_string()), Some(2))],
            "{log}"
        );
        assert!(
            !log.contains("have no record in the protein FASTA"),
            "{log}"
        );
    }

    let (rows, log) = built(&versioned, ">P-T1.1.2 pep transcript:T1.1.4\nMK\n");
    assert_eq!(rows, [("T1.1".to_string(), None, None)]);
    assert!(
        log.contains("1 of 1 translated transcripts have no record in the protein FASTA"),
        "{log}"
    );
}
