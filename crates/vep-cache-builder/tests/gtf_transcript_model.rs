// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A GTF and a genome FASTA build a cache on their own: the GTF is the
//! transcript model when no GFF3 is given, and with the genome and no protein
//! FASTA each peptide is translated, so the build reads nothing else.

use std::path::Path;
use std::process::Command;

use serde_json::Value;

/// A one-contig genome on `21`, positions 1 to 60.
fn write_genome(dir: &Path) -> std::path::PathBuf {
    let fa = dir.join("ref.fa");
    let seq = "GGGATGCCCAAATTTGGGTAGCCCTTTAAAGGGCCCATGAAACCCGGGTTTTAGCCCAAA";
    std::fs::write(&fa, format!(">21\n{seq}\n")).unwrap();
    std::fs::write(dir.join("ref.fa.fai"), "21\t60\t4\t60\t61\n").unwrap();
    fa
}

/// A gene with one transcript of two exons (1-12 and 21-30), a CDS from 4
/// to 26 and the stop codon at 27-29, on the genome above: the coding
/// sequence `ATGCCCAAAGCCCTT` plus `TAA`.
const GTF: &str = "#!genome-build GRCh38.p14\n#!genebuild-last-updated 2025-11\n\
21\tensembl_havana\tgene\t1\t30\t.\t+\t.\tgene_id \"ENSG00000000001\"; gene_version \"3\"; gene_name \"GENE1\"; gene_biotype \"protein_coding\";\n\
21\tensembl_havana\ttranscript\t1\t30\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; transcript_version \"5\"; gene_name \"GENE1\"; transcript_biotype \"protein_coding\"; tag \"Ensembl_canonical\"; tag \"gencode_basic\"; transcript_support_level \"1\";\n\
21\tensembl_havana\texon\t1\t12\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"1\"; exon_id \"ENSE00000000001\";\n\
21\tensembl_havana\texon\t21\t30\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"2\"; exon_id \"ENSE00000000002\";\n\
21\tensembl_havana\tCDS\t4\t12\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"1\"; protein_id \"ENSP00000000001\";\n\
21\tensembl_havana\tCDS\t21\t26\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"2\"; protein_id \"ENSP00000000001\";\n\
21\tensembl_havana\tstop_codon\t27\t29\t.\t+\t0\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\"; exon_number \"2\"; protein_id \"ENSP00000000001\";\n\
21\tensembl_havana\tthree_prime_utr\t30\t30\t.\t+\t.\tgene_id \"ENSG00000000001\"; transcript_id \"ENST00000000001\";\n";

/// `--gtf` with `--genome-fasta` and neither `--gff3` nor `--fasta` builds the
/// cache from the two files: the transcript carries the coding sequence through
/// the stop codon (`BaseGXF.pm` `_create_transcript`, a stop codon is a CDS
/// segment), the translation's start and end within its exons, the peptide
/// translated from it with the stop dropped (`MPKAL`, core
/// `Transcript.pm` `translate`), the tags Ensembl VEP keeps, protein version 1
/// (`_add_translation`), no gene version (a GTF writes `gene_version`, not
/// `version`), and `info.json` takes the header's assembly and genebuild.
#[test]
fn gtf_and_genome_build_a_cache_without_a_protein_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let gtf = dir.path().join("model.gtf");
    std::fs::write(&gtf, GTF).unwrap();
    let out = dir.path().join("cache");

    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "GRCh38", "--gtf"])
        .arg(&gtf)
        .arg("--genome-fasta")
        .arg(&genome)
        .arg("--output-dir")
        .arg(&out)
        .output()
        .expect("vep-cache-builder runs");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(
        stderr.contains("No protein FASTA: each peptide is translated from the genome FASTA"),
        "{stderr}"
    );
    assert!(
        stderr.contains("Skipped the feature types three_prime_utr (1)"),
        "{stderr}"
    );

    let shard = std::fs::read_to_string(out.join("transcripts/21/1-1000000.json")).unwrap();
    let records: Vec<Value> = serde_json::from_str(&shard).unwrap();
    assert_eq!(records.len(), 1);
    let tr = &records[0];
    assert_eq!(tr["stable_id"], "ENST00000000001");
    assert_eq!(tr["version"], 5);
    assert_eq!(tr["gene_stable_id"], "ENSG00000000001");
    assert!(tr.get("gene_version").is_none(), "{tr}");
    assert_eq!(tr["gene_symbol"], "GENE1");
    assert_eq!(tr["biotype"], "protein_coding");
    assert_eq!(tr["is_canonical"], 1);
    assert_eq!(tr["ccds"], "-");
    assert_eq!(tr["protein"], "ENSP00000000001");
    assert_eq!(tr["translation"]["version"], 1);
    assert_eq!(tr["translation"]["start"], 4);
    assert_eq!(tr["translation"]["end"], 9);
    assert_eq!(tr["cdna_coding_start"], 4);
    assert_eq!(tr["cdna_coding_end"], 21);
    let codes: Vec<&str> = tr["attributes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|a| a["code"].as_str().unwrap())
        .collect();
    assert_eq!(codes, ["gencode_basic", "TSL"]);
    let vefc = &tr["variation_effect_feature_cache"];
    assert_eq!(vefc["translateable_seq"], "ATGCCCAAAGCCCTTTAA");
    assert_eq!(vefc["peptide"], "MPKAL");
    assert_eq!(vefc["five_prime_utr"], "GGG");
    assert_eq!(vefc["three_prime_utr"], "A");
    assert_eq!(vefc["codon_table"], 1);

    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(out.join("info.json")).unwrap()).unwrap();
    assert_eq!(info["source_versions"]["assembly"], "GRCh38.p14");
    assert_eq!(info["source_versions"]["genebuild"], "2025-11");
    assert_eq!(info["contig_order"], serde_json::json!(["21"]));
}

/// The regulatory files bake beside the GTF model as the run-time reader
/// reads them: a promoter's `extended_start` that is not a whole number is
/// ignored with the warning naming the record and the attribute (Ensembl VEP
/// 116.2 compares the text as a number, `RegFeat.pm:363-364`); two distinct
/// records sharing an `ID` are both baked, in file order, with the collision
/// warned about once (`RegFeat.pm:478-497`); a type column that is not UTF-8
/// once its escapes are decoded, or unsupported once decoded, is skipped with
/// the warning naming the decoded type (`RegFeat.pm:312-322`); each element
/// carries its input slot and its place in the file, so the cache's rows
/// follow the file's.
#[test]
fn regulatory_files_bake_with_their_source_and_order_and_the_readers_warnings() {
    let dir = tempfile::tempdir().unwrap();
    let genome = write_genome(dir.path());
    let gtf = dir.path().join("model.gtf");
    std::fs::write(&gtf, GTF).unwrap();
    let regulatory = dir.path().join("regulatory.gff3");
    std::fs::write(
        &regulatory,
        "##gff-version 3\n\
         21\tEnsembl\tenhancer\t5\t20\t.\t.\t.\tID=TWICE\n\
         21\tEnsembl\tpromoter\t10\t30\t.\t.\t.\tID=TWICE;extended_start=x;extended_end=40\n\
         21\tEnsembl\tcaf%FF\t50\t60\t.\t.\t.\tID=NOT_UTF8_TYPE\n\
         21\tEnsembl\tcaf%CA%80\t70\t80\t.\t.\t.\tID=UNSUPPORTED_TYPE\n",
    )
    .unwrap();
    let emars = dir.path().join("EMARs.gff");
    std::fs::write(&emars, "21\tEnsembl\tEMAR\t5\t20\t.\t.\t.\tID=TWICE\n").unwrap();
    let out = dir.path().join("cache");

    let output = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args(["--release", "116", "--assembly", "GRCh38", "--gtf"])
        .arg(&gtf)
        .arg("--genome-fasta")
        .arg(&genome)
        .arg("--regulatory-gff")
        .arg(&regulatory)
        .arg("--emars-gff")
        .arg(&emars)
        .arg("--extended-promoters")
        .arg("--output-dir")
        .arg(&out)
        .output()
        .expect("vep-cache-builder runs");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    let regulatory_name = regulatory.display().to_string();
    assert!(
        stderr.contains(&format!(
            "Ignoring extended_start 'x' of promoter TWICE in {regulatory_name}; not a whole number"
        )),
        "{stderr}"
    );
    assert_eq!(
        stderr
            .matches(&format!(
                "Duplicate regulatory feature identifier 'TWICE' in {regulatory_name}; identifiers must be unique or annotations will be lost"
            ))
            .count(),
        1,
        "{stderr}"
    );
    assert!(
        !stderr.contains(&format!("'TWICE' in {}", emars.display())),
        "an identifier is unique within its file alone: {stderr}"
    );
    // A type that is not UTF-8 once decoded is skipped with the warning naming
    // it as decoded, the replacement character for the byte, as a type that
    // decodes to unsupported text is.
    for decoded_type in ["caf\u{FFFD}", "caf\u{280}"] {
        assert!(
            stderr.contains(&format!(
                "Ignoring '{decoded_type}' feature type from {regulatory_name}; not a supported regulatory feature type"
            )),
            "{stderr}"
        );
    }

    let shard = std::fs::read_to_string(out.join("regulatory/21/1-1000000.json")).unwrap();
    let records: Vec<Value> = serde_json::from_str(&shard).unwrap();
    let baked: Vec<(&str, u64, u64, u64, u64)> = records
        .iter()
        .map(|record| {
            (
                record["biotype"].as_str().unwrap(),
                record["start"].as_u64().unwrap(),
                record["end"].as_u64().unwrap(),
                record["source"].as_u64().unwrap(),
                record["ordinal"].as_u64().unwrap(),
            )
        })
        .collect();
    assert_eq!(
        baked,
        [
            ("enhancer", 5, 20, 0, 0),
            ("promoter", 10, 40, 0, 1),
            ("EMAR", 5, 20, 1, 0),
        ]
    );
}
