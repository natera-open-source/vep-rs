// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells `MAX_AF` under `--max_af`, on the rows the
//! binary writes for one record: line 9 of the custom corpus's input, whose
//! known variant in the cache of the `GRCh38-everything` corpus holds a
//! population frequency for the first of its two alleles and for whose site
//! the corpus's `custom/sites.vcf.gz` holds a record on that allele. On a row
//! that has both, the key lists the source's values and then the maximum,
//! and `MAX_AF_POPS` lists the maximum's populations. Ensembl VEP 116.3
//! compares the maximum with the list the source pushed under the key
//! (OutputFactory.pm:1267-1274) and writes the source's values alone and no
//! `MAX_AF_POPS`: the documented difference `custom_key_named_max_af`. Every
//! other cell of the expected rows is the one Ensembl VEP writes, its second
//! listing of the key dropped (`duplicate_output_column_collapsed`).

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The 1-based line of the record in the custom corpus's `variants.vcf`.
const RECORD_LINE: usize = 9;

/// One way a source's key spells `MAX_AF` and reaches a row, with the rows of
/// the default format and the column line and the rows of the tab format.
struct Shape {
    short_name: &'static str,
    /// The source's `type`: `exact` gives its record to the row of the allele
    /// it names, `overlap` to the whole input record and so to both rows.
    kind: &'static str,
    default: [&'static str; 2],
    tab: [&'static str; 3],
}

/// A source named `MAX` that asks for `AF`: the key holds the record's `AF`.
/// The cache line gives the allele `1.993e-05` in `gnomADg_NFE` and `0` in the
/// other populations, and the combined `gnomADg` column never counts
/// (OutputFactory.pm:1252).
const KEY_OF_A_FIELD: Shape = Shape {
    short_name: "MAX",
    kind: "exact",
    default: [
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;MAX_AF=0.004,1.993e-05;MAX_AF_POPS=gnomADg_NFE;MAX=exact_2;MAX_CLASS=exact;MAX_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER",
    ],
    tab: [
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tMAX\tMAX_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t0.004,1.993e-05\tgnomADg_NFE\t-\t-\t-\texact_2\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
    ],
};

/// The `CSQ` format and the record's line of the VCF under [`KEY_OF_A_FIELD`].
const KEY_OF_A_FIELD_VCF: [&str; 2] = [
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|MAX|MAX_CLASS",
    "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||||||||||||||rs1444031311|||||||0.004&1.993e-05|gnomADg_NFE||||exact_2|exact,G|intergenic_variant|MODIFIER|||||||||||||||rs1444031311|||||||||||||",
];

/// A source named `MAX_AF` itself: the key lists the record's name
/// (OutputFactory.pm:1301).
const KEY_OF_THE_NAME: Shape = Shape {
    short_name: "MAX_AF",
    kind: "exact",
    default: [
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;MAX_AF=exact_2,1.993e-05;MAX_AF_POPS=gnomADg_NFE;MAX_AF_AF=0.004;MAX_AF_CLASS=exact;MAX_AF_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER",
    ],
    tab: [
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tMAX_AF_AF\tMAX_AF_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\texact_2,1.993e-05\tgnomADg_NFE\t-\t-\t-\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
    ],
};

/// The source named `MAX` read as `type=overlap`: its record belongs to the
/// whole input record (OutputFactory.pm:1003-1014), so the second allele's
/// row, which has no frequency, holds the source's value alone, as Ensembl
/// VEP writes it.
const KEY_ON_EVERY_ROW: Shape = Shape {
    short_name: "MAX",
    kind: "overlap",
    default: [
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;MAX_AF=0.004,1.993e-05;MAX_AF_POPS=gnomADg_NFE;MAX=exact_2;MAX_CLASS=exact;MAX_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tIMPACT=MODIFIER;MAX_AF=0.004;MAX=exact_2;MAX_CLASS=exact;MAX_FILTER=PASS",
    ],
    tab: [
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tMAX\tMAX_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t0.004,1.993e-05\tgnomADg_NFE\t-\t-\t-\texact_2\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\trs1444031311\tMODIFIER\t-\t-\t-\t-\t0.004\t-\t-\t-\t-\texact_2\texact",
    ],
};

/// The JSON record Ensembl VEP 116.3 writes under [`KEY_OF_A_FIELD`], its
/// `input` key aside: the frequencies stand under `colocated_variants`, the
/// source's record under `custom_annotations`, and no entry carries the key.
const KEY_OF_A_FIELD_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","colocated_variants":[{"allele_string":"A/C/G","end":5380280,"frequencies":{"C":{"gnomadg":7.994e-06,"gnomadg_afr":0,"gnomadg_ami":0,"gnomadg_amr":0,"gnomadg_asj":0,"gnomadg_eas":0,"gnomadg_fin":0,"gnomadg_mid":0,"gnomadg_nfe":1.993e-05,"gnomadg_remaining":0,"gnomadg_sas":0}},"id":"rs1444031311","seq_region_name":"21","start":5380280,"strand":1}],"custom_annotations":{"MAX":[{"allele":"C","fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

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

/// The output of a run over the record in `format` (`default`, `tab`, `vcf`,
/// `json` or `parquet`) under `--max_af`, with the custom corpus's site file
/// as a source of the shape's short name and type that asks for its `AF` and
/// `CLASS` fields.
fn annotate(dir: &Path, format: &str, shape: &Shape) -> PathBuf {
    let out = dir.join(format!("{}.{}.{format}.out", shape.short_name, shape.kind));
    let source = format!(
        "file={},short_name={},format=vcf,type={},fields=AF%CLASS",
        corpus("GRCh38-custom")
            .join("custom/sites.vcf.gz")
            .display(),
        shape.short_name,
        shape.kind
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
        .args(["--custom", &source, "--max_af"]);
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
        "vep {format} {} {}:\n{}",
        shape.short_name,
        shape.kind,
        String::from_utf8_lossy(&run.stderr)
    );
    out
}

/// Under `--max_af` the key lists the source's values, then the row's maximum,
/// and `MAX_AF_POPS` lists its populations; a row without a frequency holds
/// the source's values alone and a row with neither holds nothing.
#[test]
fn a_custom_key_named_max_af_lists_the_sources_values_then_the_maximum() {
    let tmp = tempfile::tempdir().unwrap();
    for shape in [&KEY_OF_A_FIELD, &KEY_OF_THE_NAME, &KEY_ON_EVERY_ROW] {
        let case = format!("short_name={} type={}", shape.short_name, shape.kind);

        let text = fs::read_to_string(annotate(tmp.path(), "default", shape)).unwrap();
        let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(rows, shape.default, "default format, {case}");

        let text = fs::read_to_string(annotate(tmp.path(), "tab", shape)).unwrap();
        let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with("##")).collect();
        assert_eq!(rows, shape.tab, "tab format, {case}");
    }

    let text = fs::read_to_string(annotate(tmp.path(), "vcf", &KEY_OF_A_FIELD)).unwrap();
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
    assert_eq!(vec![format, record], KEY_OF_A_FIELD_VCF, "VCF format");
}

/// The JSON factory writes the frequencies under `colocated_variants` and
/// moves the custom annotations to the record before any row hash is built
/// (OutputFactory/JSON.pm:249-258), so the two never share a key there and
/// the record is the one Ensembl VEP writes.
#[test]
fn a_custom_key_named_max_af_leaves_the_json_record_as_ensembl_vep_writes_it() {
    let tmp = tempfile::tempdir().unwrap();
    let expected: serde_json::Value = serde_json::from_str(KEY_OF_A_FIELD_JSON).unwrap();
    let text = fs::read_to_string(annotate(tmp.path(), "json", &KEY_OF_A_FIELD)).unwrap();
    let mut record: serde_json::Value =
        serde_json::from_str(text.lines().next().expect("the record's line")).unwrap();
    let object = record.as_object_mut().expect("a JSON object");
    object.remove("input");
    object.remove("vep_rs_version");
    assert_eq!(record, expected);
}

/// The Parquet output stores the key as a string that holds the `--tab` cell:
/// the source's value and the maximum on the row that has both, null on the
/// other allele's.
#[cfg(feature = "parquet")]
#[test]
fn a_custom_key_named_max_af_is_stored_in_parquet_as_the_tab_cell() {
    let tmp = tempfile::tempdir().unwrap();
    let key = KEY_OF_A_FIELD.tab[0]
        .trim_start_matches('#')
        .split('\t')
        .position(|column| column == "MAX_AF")
        .expect("the key's column");
    assert_eq!(
        KEY_OF_A_FIELD.tab[1].split('\t').nth(key),
        Some("0.004,1.993e-05")
    );
    let out = annotate(tmp.path(), "parquet", &KEY_OF_A_FIELD);
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
    assert_eq!(cells, ["-", "0.004,1.993e-05"]);
}
