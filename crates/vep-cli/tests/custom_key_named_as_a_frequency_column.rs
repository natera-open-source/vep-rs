// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells a frequency column, on the rows the binary
//! writes for one record: line 9 of the custom corpus's input, for whose site
//! the corpus's `custom/sites.vcf.gz` holds a record, read over the cache of
//! the `GRCh38-everything` corpus, which holds the record's co-located variant
//! and its gnomAD genomes frequency. The expected rows are the ones Ensembl
//! VEP 116.3 writes for the same record and flags, with the key listed once
//! where Ensembl VEP lists it twice under the flag (the documented difference
//! `duplicate_output_column_collapsed`).

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The 1-based line of the record in the custom corpus's `variants.vcf`.
const RECORD_LINE: usize = 9;

/// The rows of the default format without `--af_gnomadg`.
const DEFAULT_OFF: &[&str] = &[
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;gnomADg=exact_2;gnomADg_AF=0.004;gnomADg_CLASS=exact;gnomADg_FILTER=PASS",
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
];

/// The column line and the rows of the tab format without `--af_gnomadg`.
const TAB_OFF: &[&str] = &[
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tgnomADg\tgnomADg_AF\tgnomADg_CLASS",
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
];

/// The `CSQ` format and the record's line of the VCF format without `--af_gnomadg`.
const VCF_OFF: &[&str] = &[
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|gnomADg|gnomADg_AF|gnomADg_CLASS",
    "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||||||||||||||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
];

/// The rows of the default format under `--af_gnomadg`.
const DEFAULT_ON: &[&str] = &[
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;gnomADg_AF=0.004,7.994e-06;gnomADg_AFR_AF=0;gnomADg_AMI_AF=0;gnomADg_AMR_AF=0;gnomADg_ASJ_AF=0;gnomADg_EAS_AF=0;gnomADg_FIN_AF=0;gnomADg_MID_AF=0;gnomADg_NFE_AF=1.993e-05;gnomADg_REMAINING_AF=0;gnomADg_SAS_AF=0;gnomADg=exact_2;gnomADg_CLASS=exact;gnomADg_FILTER=PASS",
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER",
];

/// The column line and the rows of the tab format under `--af_gnomadg`, the key
/// listed once, at the flag's place.
const TAB_ON: &[&str] = &[
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tgnomADg_AF\tgnomADg_AFR_AF\tgnomADg_AMI_AF\tgnomADg_AMR_AF\tgnomADg_ASJ_AF\tgnomADg_EAS_AF\tgnomADg_FIN_AF\tgnomADg_MID_AF\tgnomADg_NFE_AF\tgnomADg_REMAINING_AF\tgnomADg_SAS_AF\tCLIN_SIG\tSOMATIC\tPHENO\tgnomADg\tgnomADg_CLASS",
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t0.004,7.994e-06\t0\t0\t0\t0\t0\t0\t0\t1.993e-05\t0\t0\t-\t-\t-\texact_2\texact",
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
];

/// The `CSQ` format and the record's line of the VCF format under `--af_gnomadg`,
/// the key listed once, at the flag's place.
const VCF_ON: &[&str] = &[
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_REMAINING_AF|gnomADg_SAS_AF|CLIN_SIG|SOMATIC|PHENO|gnomADg|gnomADg_CLASS",
    "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||||||||||||||rs1444031311|||||||0.004&7.994e-06|0|0|0|0|0|0|0|1.993e-05|0|0||||exact_2|exact,G|intergenic_variant|MODIFIER|||||||||||||||rs1444031311||||||||||||||||||||||",
];

/// The directory of the committed corpus `name`.
fn corpus(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/116")
        .join(name)
}

/// The record alone, written under `dir`: the header lines of the custom
/// corpus's input and its line [`RECORD_LINE`].
fn record_file(dir: &Path) -> PathBuf {
    let input = fs::read_to_string(corpus("GRCh38-custom").join("variants.vcf")).unwrap();
    let record = input
        .lines()
        .nth(RECORD_LINE - 1)
        .expect("the record's line");
    assert!(
        record.starts_with("21\t5380280\t"),
        "line {RECORD_LINE} of the custom corpus's input is the record"
    );
    let mut text: String = input
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect();
    text.push_str(record);
    text.push('\n');
    let path = dir.join("record.vcf");
    fs::write(&path, text).unwrap();
    path
}

/// The output of a run over the record in `format` (`default`, `tab`, `vcf`
/// or `parquet`) under `flags`, with the custom corpus's site file as a source
/// named `gnomADg` that asks for its `AF` and `CLASS` fields.
fn annotate(dir: &Path, format: &str, flags: &[&str]) -> PathBuf {
    let out = dir.join(format!("{format}.{}.out", flags.len()));
    let source = format!(
        "file={},short_name=gnomADg,format=vcf,type=exact,fields=AF%CLASS",
        corpus("GRCh38-custom")
            .join("custom/sites.vcf.gz")
            .display()
    );
    let mut cmd = common::isolated_vep(dir);
    cmd.arg("-i")
        .arg(record_file(dir))
        .arg("-o")
        .arg(&out)
        .arg("--json_cache")
        .arg(corpus("GRCh38-everything").join("json_cache"))
        .args(["--offline", "--species", "homo_sapiens"])
        .args(["--assembly", "GRCh38", "--no_stats", "--force_overwrite"])
        .args(["--custom", &source])
        .args(flags);
    match format {
        "default" => {}
        "parquet" => {
            cmd.args(["--output_format", "parquet"]);
        }
        other => {
            cmd.arg(format!("--{other}"));
        }
    }
    let run = cmd.output().unwrap();
    assert!(
        run.status.success(),
        "vep {format} {flags:?}:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    out
}

/// Ensembl VEP keys the line hash by the whole string, so the key
/// `gnomADg_AF` holds the source's value (OutputFactory.pm:1300-1304) and,
/// under `--af_gnomadg`, the co-located variant's frequency merged after it
/// (:1245-1248). The default, tab and VCF outputs print that value on the row
/// of the allele the source's record names, and nothing on the other allele's.
#[test]
fn a_custom_key_named_as_a_frequency_column_prints_the_sources_value_then_the_frequency() {
    let tmp = tempfile::tempdir().unwrap();
    for (flags, default, tab, vcf) in [
        (&[][..], DEFAULT_OFF, TAB_OFF, VCF_OFF),
        (&["--af_gnomadg"][..], DEFAULT_ON, TAB_ON, VCF_ON),
    ] {
        let text = fs::read_to_string(annotate(tmp.path(), "default", flags)).unwrap();
        let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(rows, default, "default format {flags:?}");

        let text = fs::read_to_string(annotate(tmp.path(), "tab", flags)).unwrap();
        let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with("##")).collect();
        assert_eq!(rows, tab, "tab format {flags:?}");

        let text = fs::read_to_string(annotate(tmp.path(), "vcf", flags)).unwrap();
        let format = text
            .lines()
            .find_map(|l| l.strip_prefix("##INFO=<ID=CSQ,"))
            .and_then(|l| l.split_once("Format: "))
            .map(|(_, format)| format.trim_end_matches("\">"))
            .expect("the CSQ header line");
        let record = text
            .lines()
            .find(|l| !l.starts_with('#'))
            .expect("the record's line");
        assert_eq!([format, record], vcf, "VCF format {flags:?}");
    }
}

/// The Parquet output stores the key as a string that holds the `--tab` cell:
/// the value on the row of the allele the source's record names, null on the
/// other allele's.
#[cfg(feature = "parquet")]
#[test]
fn a_custom_key_named_as_a_frequency_column_is_stored_in_parquet_as_the_tab_cell() {
    let tmp = tempfile::tempdir().unwrap();
    for (flags, tab, cell) in [
        (&[][..], TAB_OFF, "0.004"),
        (&["--af_gnomadg"][..], TAB_ON, "0.004,7.994e-06"),
    ] {
        let key = tab[0]
            .trim_start_matches('#')
            .split('\t')
            .position(|column| column == "gnomADg_AF")
            .expect("the key's column");
        let out = annotate(tmp.path(), "parquet", flags);
        let mut cells: Vec<String> = vep_io::output::parquet::tab_rows(&out)
            .unwrap()
            .iter()
            .map(|row| {
                row.split('\t')
                    .nth(key)
                    .expect("the key's cell")
                    .to_string()
            })
            .collect();
        cells.sort_unstable();
        assert_eq!(cells, ["-", cell], "{flags:?}");
    }
}
