// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A start codon's terms on a cached model that shows no start codon, read
//! through the `vep` binary on the cache of a committed corpus.
//!
//! The coding sequence of the human model `ENST00000420031`, on the reverse
//! strand of chromosome 1 of GRCh37, opens at `1:16901033` behind a 5' UTR of
//! 18 bases. Its record carries no `cds_start_NF` and no consensus coding
//! sequence id, its first codon names Ile, which is no start codon of the
//! standard code, and its recorded translation opens with Ile: no start
//! codon stands at the coding start. A substitution of the codon's third
//! base, `1:16901031`, that keeps Ile is a synonymous change and nothing
//! else. Ensembl VEP reads the first three coding positions as the start
//! codon on every transcript without the flag and writes `start_lost` beside
//! `synonymous_variant` (`VariationEffect.pm:975-986`, `:941`).
//!
//! The cache is that of the committed corpus `GRCh37`, converted from the
//! Ensembl VEP cache of release 116.

use std::path::PathBuf;
use std::process::Stdio;

mod common;

/// The substitution of the third base of the model's first codon that keeps
/// the residue.
const RECORD: &str = "##fileformat=VCFv4.2\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                      1\t16901031\tresidue_kept_1_16901031\tA\tT\t.\t.\t.\n";

/// The model whose coding sequence opens on a codon of Ile at `1:16901033`.
const TRANSCRIPT: &str = "ENST00000420031";

/// The cache of the committed corpus `GRCh37`.
fn cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh37/json_cache")
}

#[test]
fn substitution_that_keeps_the_first_residue_of_a_cached_model_with_no_start_codon_is_synonymous_alone(
) {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("record.vcf");
    std::fs::write(&input, RECORD).unwrap();
    let out = dir.path().join("record.out");
    let output = common::isolated_vep(dir.path())
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(cache())
        .args(["--assembly", "GRCh37", "--offline"])
        .args(["--force_overwrite", "--no_stats", "--no_headers", "--quiet"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&out).unwrap();
    let row: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|columns| columns[4] == TRANSCRIPT)
        .unwrap_or_else(|| panic!("no row of {TRANSCRIPT}: {text}"));
    // Consequence, cDNA_position, CDS_position, Protein_position, Amino_acids,
    // Codons and Extra.
    assert_eq!(
        [row[6], row[7], row[8], row[9], row[10], row[11], row[13]],
        [
            "synonymous_variant",
            "21",
            "3",
            "1",
            "I",
            "atT/atA",
            "IMPACT=LOW;STRAND=-1"
        ]
    );
}
