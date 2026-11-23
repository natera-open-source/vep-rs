// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` GFF source in `gff_type=gencode_promoter` mode through the
//! `vep` binary: the `feature_id` and `associated_gene` fields carry the
//! characters a record's `ID` and `gene_id` attributes encode with `%XX`
//! escapes, in every output format, where Ensembl VEP 116.3 prints the
//! escapes as written (`GFF.pm:199-202` over the raw hash of
//! `GFF3.pm:190-198`, the documented `gff_attribute_escapes_read_raw`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Stdio;

mod common;

fn corpus_cache() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116/GRCh38-release116/json_cache")
}

/// Runs the binary over the one-record input on the promoter source `gff`,
/// writing `<dir>/<name>` in the format `format` selects (`""` for the
/// default), and returns the output text.
fn run_overlay(dir: &Path, gff: &Path, name: &str, format: &str) -> String {
    let input = dir.join("in.vcf");
    fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t30282169\tpromoter_site\tC\tT\t.\t.\t.\n",
    )
    .unwrap();
    let out = dir.join(name);
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus_cache())
        .args([
            "--assembly",
            "GRCh38",
            "--force_overwrite",
            "--quiet",
            "--custom",
        ])
        .arg(format!(
            "file={},format=gff,gff_type=gencode_promoter",
            gff.display()
        ))
        .stdin(Stdio::null());
    if !format.is_empty() {
        cmd.arg(format);
    }
    let output = cmd.output().unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(0), "{name}: {stderr}");
    fs::read_to_string(out).unwrap()
}

/// The rows of a text output, its `#` lines left out.
fn rows(text: &str) -> Vec<&str> {
    text.lines().filter(|line| !line.starts_with('#')).collect()
}

/// A window `ID=PROM%2C1;gene_id=ENSG%3B1` over the record `21:30282169 C>T`
/// gives `feature_id` `PROM,1` and `associated_gene` `ENSG;1`: the tab and
/// JSON formats print the strings; the default format's `Extra` column
/// re-spells a `;` inside a value `%3B` (`OutputFactory/VEP_output.pm:118`),
/// so it prints `PROM,1` and `ENSG%3B1`; the VCF format's `CSQ` entry prints
/// `PROM&1` and `ENSG%3B1`, the spelling every `CSQ` value gives a `,` and a
/// `;` (`OutputFactory/VCF.pm:401-404`). Ensembl VEP 116.3 prints `PROM%2C1`
/// and `ENSG%3B1` in all four formats, the raw read of `GFF3.pm:190-198`.
#[test]
fn promoter_identifiers_carry_the_characters_their_escapes_encode() {
    let dir = tempfile::tempdir().unwrap();
    let gff = common::indexed_gff(
        dir.path(),
        "promoters.gff.gz",
        &["21\t.\tpromoter_window\t30282160\t30282180\t.\t+\t.\tID=PROM%2C1;gene_id=ENSG%3B1"],
    );

    let default = run_overlay(dir.path(), &gff, "default.txt", "");
    let default_rows = rows(&default);
    assert!(!default_rows.is_empty(), "{default}");
    for row in &default_rows {
        assert!(
            row.ends_with(
                "GENCODE_PROMOTER=21:30282160-30282180;\
                 GENCODE_PROMOTER_associated_gene=ENSG%3B1;\
                 GENCODE_PROMOTER_feature_id=PROM,1;\
                 GENCODE_PROMOTER_type=promoter_window"
            ),
            "{row}"
        );
    }

    let tab = run_overlay(dir.path(), &gff, "tab.txt", "--tab");
    let header = tab
        .lines()
        .find(|line| line.starts_with("#Uploaded_variation"))
        .unwrap();
    let columns: Vec<&str> = header.split('\t').collect();
    let column = |name: &str| columns.iter().position(|column| *column == name).unwrap();
    let (feature_id, associated_gene) = (
        column("GENCODE_PROMOTER_feature_id"),
        column("GENCODE_PROMOTER_associated_gene"),
    );
    let tab_rows = rows(&tab);
    assert!(!tab_rows.is_empty(), "{tab}");
    for row in &tab_rows {
        let cells: Vec<&str> = row.split('\t').collect();
        assert_eq!(cells[feature_id], "PROM,1", "{row}");
        assert_eq!(cells[associated_gene], "ENSG;1", "{row}");
    }

    let vcf = run_overlay(dir.path(), &gff, "out.vcf", "--vcf");
    let vcf_rows = rows(&vcf);
    assert_eq!(vcf_rows.len(), 1, "{vcf}");
    let info = vcf_rows[0].split('\t').nth(7).unwrap();
    let csq = info
        .split(';')
        .find_map(|item| item.strip_prefix("CSQ="))
        .unwrap();
    for entry in csq.split(',') {
        assert!(
            entry.ends_with("|21:30282160-30282180|promoter_window|PROM&1|ENSG%3B1"),
            "{entry}"
        );
    }

    let json = run_overlay(dir.path(), &gff, "out.json", "--json");
    let json_rows = rows(&json);
    assert_eq!(json_rows.len(), 1, "{json}");
    let record: serde_json::Value = serde_json::from_str(json_rows[0]).unwrap();
    assert_eq!(
        record["custom_annotations"]["GENCODE_PROMOTER"],
        serde_json::json!([{
            "name": "21:30282160-30282180",
            "fields": {
                "type": "promoter_window",
                "feature_id": "PROM,1",
                "associated_gene": "ENSG;1",
            },
        }]),
        "{json}"
    );
}
