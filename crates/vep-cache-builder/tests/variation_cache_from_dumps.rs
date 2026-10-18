// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The variation shards the binary writes with `--include-variations` carry
//! each record in Ensembl VEP's form and, with `--variation-dumps`, the
//! ClinVar allele and somatic columns of the variation database dumps.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/grch38_chr21")
}

/// A BGZF-compressed variation VCF named as the Ensembl FTP names them, so the
/// `--variation-vcf-dir` scan finds it.
fn write_vcf(dir: &Path, records: &str) {
    let path = dir.join("homo_sapiens-chr21.vcf.gz");
    let mut writer = noodles::bgzf::io::Writer::new(std::fs::File::create(path).unwrap());
    let header = "##fileformat=VCFv4.2\n\
         ##INFO=<ID=E_Phenotype_or_Disease,Number=0,Type=Flag,Description=\"Phenotype_or_Disease\">\n\
         ##INFO=<ID=CLIN_pathogenic,Number=0,Type=Flag,Description=\"pathogenic\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    writer.write_all(header.as_bytes()).unwrap();
    writer.write_all(records.as_bytes()).unwrap();
    writer.finish().unwrap();
}

/// The variation database tables of one ClinVar span, 21:101-102: a germline
/// feature of rs1 (risk allele `-`, significance `pathogenic`, reference allele
/// `TG`) and a somatic feature with a classification.
fn write_dumps(dir: &Path) {
    let tables: &[(&str, &str)] = &[
        (
            "source",
            "1\tdbSNP\t156\t\\N\t\\N\t\\N\tgermline\tvariation\n32\tClinVar\t202509\t\\N\t\\N\t\\N\tmixed\tphenotype_feature\n",
        ),
        ("seq_region", "131537\t21\t4\n"),
        (
            "attrib_type",
            "480\tclinvar_clin_sig\tClinVar clinical significance\t\\N\n483\trisk_allele\tRisk allele\t\\N\n577\tsomatic_clin_sig\tClinVar somatic classification\t\\N\n583\tclinvar_ref_allele\tClinVar reference allele\t\\N\n",
        ),
        ("phenotype", "7\t\\N\t\\N\tBreast cancer\t\\N\n"),
        (
            "phenotype_feature",
            "1001\t7\t32\t\\N\tVariation\trs1\t1\t131537\t101\t102\t1\tGermline\n1002\t7\t32\t\\N\tVariation\trs1\t1\t131537\t101\t102\t1\tSomatic\n",
        ),
        (
            "phenotype_feature_attrib",
            "1001\t480\tpathogenic\n1001\t483\t-\n1001\t583\tTG\n1002\t577\tTier I - Strong\n",
        ),
    ];
    for (table, rows) in tables {
        std::fs::write(dir.join(format!("{table}.txt")), rows).unwrap();
    }
}

/// A `TG` deletion written as `ATG/A` at 100 is cached as `TG/-` at 101-102 with
/// the evidence flag read, and the ClinVar columns join on that span; a record on
/// another span carries no ClinVar column.
#[test]
fn variation_shards_carry_ensembl_form_and_the_clinvar_columns() {
    let fixture = fixture_dir();
    let inputs = tempfile::tempdir().unwrap();
    write_vcf(
        inputs.path(),
        "21\t100\trs1\tATG\tA\t.\t.\tE_Phenotype_or_Disease;CLIN_pathogenic\n\
         21\t300\trs2\tA\tG\t.\t.\t.\n",
    );
    write_dumps(inputs.path());
    let out = tempfile::tempdir().unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_vep-cache-builder"))
        .args([
            "--release",
            "116",
            "--assembly",
            "GRCh38",
            "--include-variations",
        ])
        .arg("--gff3")
        .arg(fixture.join("slice.gff3"))
        .arg("--fasta")
        .arg(fixture.join("peptides.fa"))
        .arg("--output-dir")
        .arg(out.path())
        .arg("--variation-vcf-dir")
        .arg(inputs.path())
        .arg("--variation-dumps")
        .arg(inputs.path())
        .status()
        .expect("vep-cache-builder runs");
    assert!(status.success(), "vep-cache-builder exited {status}");

    let shard = out.path().join("variations/21/1-1000000.json");
    let records: Vec<Value> =
        serde_json::from_str(&std::fs::read_to_string(&shard).unwrap()).unwrap();
    assert_eq!(records.len(), 2);

    let deletion = &records[0];
    assert_eq!(deletion["variation_name"], "rs1");
    assert_eq!(deletion["start"], 101);
    assert_eq!(deletion["end"], 102);
    assert_eq!(deletion["allele_string"], "TG/-");
    assert_eq!(deletion["phenotype_or_disease"], 1);
    assert_eq!(deletion["clin_sig_allele"], "-:pathogenic");
    assert_eq!(deletion["clin_sig_ref_allele"], "TG");
    assert_eq!(
        deletion["clinical_impact"],
        "Breast cancer(Tier I - Strong)"
    );

    let snv = &records[1];
    assert_eq!(snv["allele_string"], "A/G");
    assert_eq!(snv["phenotype_or_disease"], 0);
    assert!(snv.get("clin_sig_allele").is_none());
    assert!(snv.get("clinical_impact").is_none());

    let info: Value =
        serde_json::from_str(&std::fs::read_to_string(out.path().join("info.json")).unwrap())
            .unwrap();
    let columns: Vec<&str> = info["variation_cols"]
        .as_array()
        .expect("variation_cols is written with the records")
        .iter()
        .map(|c| c.as_str().unwrap())
        .collect();
    assert_eq!(columns[0], "variation_name");
    for key in deletion.as_object().unwrap().keys() {
        assert!(
            columns.contains(&key.as_str()),
            "{key} is a record key but no column"
        );
    }
}
