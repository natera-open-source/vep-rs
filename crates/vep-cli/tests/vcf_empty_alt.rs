// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! How the `vep` binary reads a VCF record with an empty ALT between two
//! others, a shape the VCF format does not define: the record's working string
//! keeps the empty part, as Ensembl VEP's does, and the part is the allele `-`
//! on its rows.

mod common;

use std::path::Path;
use std::process::Stdio;

use common::{corpora, Corpus};

fn corpus() -> Corpus {
    corpora()
        .into_iter()
        .find(|c| c.name == "GRCh38-vcf-edge-cases")
        .expect("the edge-cases corpus is committed")
}

/// The data lines of a run over `21 30289183 . C G,,T` with the corpus cache and
/// `extra`.
fn run(dir: &Path, corpus: &Corpus, extra: &[&str]) -> Vec<String> {
    let input = dir.join("empty-alt.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30289183\t.\tC\tG,,T\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.join(format!("empty-alt{}.out", extra.join("")));
    let mut cmd = common::isolated_vep_anywhere();
    cmd.arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--offline");
    cmd.args(corpus.cache_args());
    cmd.args([
        "--species",
        "homo_sapiens",
        "--assembly",
        &corpus.assembly,
        "--force_overwrite",
        "--no_stats",
        "--quiet",
        "--allele_number",
        "--uploaded_allele",
    ]);
    cmd.args(extra);
    let output = cmd.stdin(Stdio::null()).output().unwrap();
    assert!(
        output.status.success(),
        "{extra:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// `21 30289183 . C G,,T`: the first characters `C`, `G`, none and `T` are not one,
/// so Parser/VCF.pm:308 chops nothing and :342 joins the parts as written into the
/// allele string `C/G//T`, which the JSON `allele_string` prints, where :316 would
/// have written `-` for a part the chop emptied. The empty part is the allele `-`
/// to every overlap (VariationFeatureOverlap.pm:462), numbered 2 between `G` and
/// `T` (:483). Ensembl VEP 116.2 writes `"allele_string":"C/G//T"` and three rows
/// on ENST00000416044: `G` (`stop_lost`), `-` and `T` (`stop_retained_variant`),
/// each named `21_30289183_C/G//T` with `UPLOADED_ALLELE=C/G//T`.
#[test]
fn an_empty_alt_between_others_keeps_the_working_string_of_ensembl_vep() {
    let dir = tempfile::tempdir().unwrap();
    let corpus = corpus();
    let json = run(dir.path(), &corpus, &["--json"]);
    assert_eq!(json.len(), 1, "{json:?}");
    let object: serde_json::Value = serde_json::from_str(&json[0]).unwrap();
    assert_eq!(object["allele_string"], "C/G//T");
    let alleles: Vec<(&str, u64)> = object["transcript_consequences"]
        .as_array()
        .unwrap()
        .iter()
        .map(|tc| {
            (
                tc["variant_allele"].as_str().unwrap(),
                tc["allele_num"].as_u64().unwrap(),
            )
        })
        .collect();
    assert_eq!(alleles, [("G", 1), ("-", 2), ("T", 3)]);

    let rows = run(dir.path(), &corpus, &[]);
    let alleles: Vec<(&str, &str, &str)> = rows
        .iter()
        .map(|row| {
            let columns: Vec<&str> = row.split('\t').collect();
            (columns[0], columns[2], columns[6])
        })
        .collect();
    assert_eq!(
        alleles,
        [
            ("21_30289183_C/G//T", "G", "stop_lost"),
            ("21_30289183_C/G//T", "-", "stop_retained_variant"),
            ("21_30289183_C/G//T", "T", "stop_retained_variant"),
        ]
    );
    for (row, number) in rows.iter().zip(1..) {
        assert!(
            row.ends_with(&format!(
                "\tALLELE_NUM={number};UPLOADED_ALLELE=C/G//T;IMPACT={};STRAND=-1",
                if number == 1 { "HIGH" } else { "LOW" }
            )),
            "{row}"
        );
    }
}
