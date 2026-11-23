// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A stop codon that polyadenylation completes, read through the `vep` binary
//! on a model built from a gene-model file and on the cache's record of the
//! same model.
//!
//! The coding sequence of the human mitochondrial model `ENST00000361390`
//! ends two bases into a codon on the transcript's last base, `MT:4262`, and
//! those bases followed by an adenine are a stop codon of the vertebrate
//! mitochondrial code. The record carries codon table 2 and no flag in either
//! cache, which shows the model complete, so a substitution of the codon's
//! first base, `MT:4261`, is read on the completed codon: the stop codon is
//! lost, the residue cell holds the stop against the new residue, the codon
//! cell shows the completing adenine as an unchanged base and the positions
//! are the record's own. Ensembl VEP reads the allele as lying in an
//! incomplete terminal codon (`VariationEffect.pm:1601-1610`).
//!
//! The two caches are those of committed corpora: `GRCh38-native-builder`,
//! which `vep-cache-builder` writes from the release's GFF3 file, and `GRCh38`,
//! converted from the Ensembl VEP cache. Both hold the model with the same
//! coding sequence.

use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

/// The substitution of the first base of the stop codon the model's coding
/// sequence ends inside.
const RECORD: &str = "##fileformat=VCFv4.2\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                      MT\t4261\t.\tT\tC\t.\t.\t.\n";

/// The model whose coding sequence ends at `MT:4262`.
const TRANSCRIPT: &str = "ENST00000361390";

/// The cache of the committed corpus `name`.
fn cache_of(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
        .join("json_cache")
}

/// `(Consequence, cDNA_position, CDS_position, Protein_position, Amino_acids,
/// Codons)` of the record's row on the model, annotated against the cache of
/// the corpus `name`.
fn row_on(name: &str, dir: &Path) -> [String; 6] {
    let input = dir.join("record.vcf");
    std::fs::write(&input, RECORD).unwrap();
    let out = dir.join(format!("{name}.out"));
    let output = common::isolated_vep(dir)
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache_of(name))
        .args(["--assembly", "GRCh38", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{name}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    let row: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|columns| columns[4] == TRANSCRIPT)
        .unwrap_or_else(|| panic!("no row of {TRANSCRIPT} on {name}: {text}"));
    [row[6], row[7], row[8], row[9], row[10], row[11]].map(str::to_string)
}

#[test]
fn substitution_in_a_stop_codon_that_adenines_complete_reads_alike_on_a_gene_model_and_the_cache() {
    let dir = tempfile::tempdir().unwrap();
    let from_the_gene_model = row_on("GRCh38-native-builder", dir.path());
    let from_the_cache = row_on("GRCh38", dir.path());
    assert_eq!(from_the_gene_model, from_the_cache);
    assert_eq!(
        from_the_gene_model,
        ["stop_lost", "955", "955", "319", "*/Q", "Taa/Caa"]
    );
}
