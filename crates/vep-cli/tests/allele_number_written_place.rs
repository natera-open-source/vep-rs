// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `--allele_number` on a record whose allele list repeats an alternate or names
//! the reference, through the `vep` binary: every alternate that writes a row
//! prints the number of its place in the record, in each output format and from
//! a VCF line as from an Ensembl-format line.

use std::path::{Path, PathBuf};

mod common;

/// The GRCh38-vcf-edge-cases corpus. Position 21:30289183 is the middle base of
/// the stop codon of ENST00000416044 (reverse strand), the one transcript of its
/// cache that the records below reach; the reference base is `C`.
fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-vcf-edge-cases")
}

/// `dup_alt` is line 15 of the corpus's `variants.vcf`: its second alternate
/// repeats the first. The first alternate of `after_reference` is the reference
/// base, and `after_star` opens with a `*`, which holds a place and writes no row.
const VCF_RECORDS: &str = "##fileformat=VCFv4.2\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
    21\t30289183\tdup_alt\tC\tG,G,T\t.\t.\t.\n\
    21\t30289183\tafter_reference\tC\tC,T\t.\t.\t.\n\
    21\t30289183\tafter_star\tC\t*,G,G,T\t.\t.\t.\n";

/// The first two records as Ensembl-format lines.
const ENSEMBL_RECORDS: &str = "21\t30289183\t30289183\tC/G/G/T\t+\tens_dup\n\
    21\t30289183\t30289183\tC/C/T\t+\tens_after_reference\n";

/// The lines of the output `vep` writes for `text`, saved as `name`, under the
/// corpus's own flags and `format`, a format flag or none; header lines whose
/// first two characters are `##` are left out.
fn lines(dir: &Path, name: &str, text: &str, format: Option<&str>) -> Vec<String> {
    let input = dir.join(name);
    std::fs::write(&input, text).unwrap();
    let output = dir.join(format!("{name}.{}.out", format.unwrap_or("default")));
    let result = common::isolated_vep_anywhere()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .args([
            "--json_cache",
            corpus().join("json_cache").to_str().unwrap(),
        ])
        .args(["--offline", "--assembly", "GRCh38"])
        .args(["--allele_number", "--uploaded_allele", "--variant_class"])
        .args(format)
        .args(["--force_overwrite", "--no_stats", "--quiet"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{name} {format:?}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::read_to_string(&output)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with("##"))
        .map(str::to_string)
        .collect()
}

/// A default-format row of ENST00000416044 at the position: `allele` with its
/// term, residues, codons and impact, numbered `number` beside the record's
/// allele list and class.
fn row(name: &str, allele: &str, number: u32, uploaded: &str, class: &str) -> String {
    let (term, residues, codons, impact) = match allele {
        "G" => ("stop_lost", "*/S", "tGa/tCa", "HIGH"),
        _ => ("stop_retained_variant", "*", "tGa/tAa", "LOW"),
    };
    format!(
        "{name}\t21:30289183\t{allele}\tENSG00000232263\tENST00000416044\tTranscript\t{term}\t\
         332\t308\t103\t{residues}\t{codons}\t-\t\
         ALLELE_NUM={number};UPLOADED_ALLELE={uploaded};IMPACT={impact};STRAND=-1;\
         VARIANT_CLASS={class}"
    )
}

/// The third alternate of `dup_alt` prints 3, the alternate after the reference
/// base 2, and the alternates after the `*` and the repeat 2 and 4. Ensembl VEP
/// 116.3 drops a repeated allele and an alternate equal to the reference before
/// it counts (VariationFeatureOverlap.pm:463, `:477-490`) and prints 2, 1, and 2
/// and 3 on those rows (`docs/divergences/allele-number-repeated-or-reference-equal-alt.md`);
/// its stored row for `dup_alt` is the second row below with `ALLELE_NUM=2`.
#[test]
fn an_alternate_prints_the_number_of_its_place_in_the_record() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        lines(dir.path(), "in.vcf", VCF_RECORDS, None)[1..],
        [
            row("dup_alt", "G", 1, "C/G/G/T", "SNV"),
            row("dup_alt", "T", 3, "C/G/G/T", "SNV"),
            row("after_reference", "T", 2, "C/C/T", "SNV"),
            row("after_star", "G", 2, "C/*/G/G/T", "sequence_alteration"),
            row("after_star", "T", 4, "C/*/G/G/T", "sequence_alteration"),
        ]
    );
}

/// The tab column, the `CSQ` field and the JSON `allele_num` of the same rows.
#[test]
fn every_output_format_prints_the_place() {
    let dir = tempfile::tempdir().unwrap();
    let tab: Vec<(String, String, String)> =
        lines(dir.path(), "in.vcf", VCF_RECORDS, Some("--tab"))[1..]
            .iter()
            .map(|line| {
                let cells: Vec<&str> = line.split('\t').collect();
                (
                    cells[0].to_string(),
                    cells[2].to_string(),
                    cells[13].to_string(),
                )
            })
            .collect();
    let expected: Vec<(String, String, String)> = [
        ("dup_alt", "G", "1"),
        ("dup_alt", "T", "3"),
        ("after_reference", "T", "2"),
        ("after_star", "G", "2"),
        ("after_star", "T", "4"),
    ]
    .iter()
    .map(|(name, allele, number)| (name.to_string(), allele.to_string(), number.to_string()))
    .collect();
    assert_eq!(tab, expected);

    let vcf = lines(dir.path(), "in.vcf", VCF_RECORDS, Some("--vcf"));
    let fields = |record: usize| -> Vec<(String, String)> {
        let info = vcf[record].split('\t').nth(7).unwrap();
        info.strip_prefix("CSQ=")
            .unwrap()
            .split(',')
            .map(|entry| {
                let cells: Vec<&str> = entry.split('|').collect();
                (cells[0].to_string(), cells[18].to_string())
            })
            .collect()
    };
    let pairs = |rows: &[(&str, &str)]| -> Vec<(String, String)> {
        rows.iter()
            .map(|(allele, number)| (allele.to_string(), number.to_string()))
            .collect()
    };
    assert_eq!(fields(1), pairs(&[("G", "1"), ("T", "3")]));
    assert_eq!(fields(2), pairs(&[("T", "2")]));
    assert_eq!(fields(3), pairs(&[("G", "2"), ("T", "4")]));

    let json: Vec<(String, String, u64)> = lines(dir.path(), "in.vcf", VCF_RECORDS, Some("--json"))
        .iter()
        .flat_map(|line| {
            let record: serde_json::Value = serde_json::from_str(line).unwrap();
            let name = record["id"].as_str().unwrap().to_string();
            record["transcript_consequences"]
                .as_array()
                .unwrap()
                .iter()
                .map(|entry| {
                    (
                        name.clone(),
                        entry["variant_allele"].as_str().unwrap().to_string(),
                        entry["allele_num"].as_u64().unwrap(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let expected: Vec<(String, String, u64)> = [
        ("dup_alt", "G", 1),
        ("dup_alt", "T", 3),
        ("after_reference", "T", 2),
        ("after_star", "G", 2),
        ("after_star", "T", 4),
    ]
    .iter()
    .map(|(name, allele, number)| (name.to_string(), allele.to_string(), *number))
    .collect();
    assert_eq!(json, expected);
}

/// An Ensembl-format line numbers its alleles by their place in its allele
/// string, as the VCF record of the same alleles does.
#[test]
fn an_ensembl_format_line_prints_the_place_in_its_allele_string() {
    let dir = tempfile::tempdir().unwrap();
    assert_eq!(
        lines(dir.path(), "in.txt", ENSEMBL_RECORDS, None)[1..],
        [
            row("ens_dup", "G", 1, "C/G/G/T", "SNV"),
            row("ens_dup", "T", 3, "C/G/G/T", "SNV"),
            row("ens_after_reference", "T", 2, "C/C/T", "SNV"),
        ]
    );
}
