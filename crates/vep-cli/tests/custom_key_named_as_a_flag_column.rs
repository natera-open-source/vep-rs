// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells a flag column other than a frequency column,
//! on the rows the binary writes for one record: line 9 of the custom corpus's
//! input, for whose site the corpus's `custom/sites.vcf.gz` holds a record on
//! one of the two alleles. A source named `VARIANT` that asks for `CLASS`
//! declares `VARIANT_CLASS`, one named `MAX` that asks for `AF` declares
//! `MAX_AF`, and one named `SYMBOL` lists its record names under `SYMBOL`.
//! The expected rows are the ones Ensembl VEP 116.3 writes for the same record
//! and flags; where it lists the key a second time, vep-rs lists it once (the
//! documented difference `duplicate_output_column_collapsed`).

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The 1-based line of the record in the custom corpus's `variants.vcf`.
const RECORD_LINE: usize = 9;

/// What Ensembl VEP 116.3 writes for the record under one short name and one
/// flag set: the rows of the default format, the column line and the rows of
/// the tab format, and the `CSQ` format and the record's line of the VCF.
struct Reference {
    short_name: &'static str,
    /// The key of the source that spells a flag column.
    key: &'static str,
    flags: &'static [&'static str],
    default: [&'static str; 2],
    tab: [&'static str; 3],
    vcf: [&'static str; 2],
}

/// `VARIANT_CLASS` without `--variant_class`.
const VARIANT_CLASS_OFF: Reference = Reference {
    short_name: "VARIANT",
    key: "VARIANT_CLASS",
    flags: &[],
    default: [
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;VARIANT=exact_2;VARIANT_AF=0.004;VARIANT_CLASS=exact;VARIANT_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: [
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tVARIANT\tVARIANT_AF\tVARIANT_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: [
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|VARIANT|VARIANT_AF|VARIANT_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||||||||||||||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `MAX_AF` without `--max_af`.
const MAX_AF_OFF: Reference = Reference {
    short_name: "MAX",
    key: "MAX_AF",
    flags: &[],
    default: [
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;MAX=exact_2;MAX_AF=0.004;MAX_CLASS=exact;MAX_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: [
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tMAX\tMAX_AF\tMAX_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: [
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|MAX|MAX_AF|MAX_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER||||||||||||||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `SYMBOL` without `--symbol`. The `CSQ` format lists `SYMBOL` whatever the
/// flag, so Ensembl VEP's lists it at places 4 and 25.
const SYMBOL_OFF: Reference = Reference {
    short_name: "SYMBOL",
    key: "SYMBOL",
    flags: &[],
    default: [
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;SYMBOL=exact_2;SYMBOL_AF=0.004;SYMBOL_CLASS=exact;SYMBOL_FILTER=PASS",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER",
    ],
    tab: [
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSOURCE\tSYMBOL\tSYMBOL_AF\tSYMBOL_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: [
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|SYMBOL|SYMBOL_AF|SYMBOL_CLASS",
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|exact_2|||||||||||||||||||||exact_2|0.004|exact,G|intergenic_variant|MODIFIER||||||||||||||||||||||||",
    ],
};

/// `SYMBOL` under `--symbol`, on rows that name no transcript: Ensembl VEP's
/// tab output lists the key at columns 18 and 22.
const SYMBOL_ON: Reference = Reference {
    short_name: "SYMBOL",
    key: "SYMBOL",
    flags: &["--symbol"],
    default: SYMBOL_OFF.default,
    tab: [
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tSYMBOL\tSYMBOL_SOURCE\tHGNC_ID\tSOURCE\tSYMBOL\tSYMBOL_AF\tSYMBOL_CLASS",
        "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\texact_2\t-\t-\t-\texact_2\t0.004\texact",
        "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
    ],
    vcf: SYMBOL_OFF.vcf,
};

/// The record Ensembl VEP 116.3 writes in JSON under a source named `SYMBOL`,
/// with and without `--symbol`, its `input` key aside: the source's record
/// stands under `custom_annotations` and no entry carries a `gene_symbol`.
const SYMBOL_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"SYMBOL":[{"allele":"C","fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1}"#;

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
/// or `parquet`) under the reference's flags, with the custom corpus's site
/// file as a source of the reference's short name that asks for its `AF` and
/// `CLASS` fields.
fn annotate(dir: &Path, format: &str, reference: &Reference) -> PathBuf {
    let out = dir.join(format!(
        "{}.{format}.{}.out",
        reference.short_name,
        reference.flags.len()
    ));
    let source = format!(
        "file={},short_name={},format=vcf,type=exact,fields=AF%CLASS",
        corpus("GRCh38-custom")
            .join("custom/sites.vcf.gz")
            .display(),
        reference.short_name
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
        .args(reference.flags);
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
        "vep {format} {} {:?}:\n{}",
        reference.short_name,
        reference.flags,
        String::from_utf8_lossy(&run.stderr)
    );
    out
}

/// `line` cut at `separator` without the cell at `place`, counted from 0.
fn without_cell(line: &str, separator: char, place: usize) -> String {
    let cells: Vec<&str> = line
        .split(separator)
        .enumerate()
        .filter(|(i, _)| *i != place)
        .map(|(_, cell)| cell)
        .collect();
    cells.join(&separator.to_string())
}

/// The place, counted from 0, at which `names` lists `key` a second time.
fn second_listing(names: &str, separator: char, key: &str) -> Option<usize> {
    names
        .trim_start_matches('#')
        .split(separator)
        .enumerate()
        .filter(|(_, name)| *name == key)
        .map(|(i, _)| i)
        .nth(1)
}

/// The reference's tab lines as vep-rs writes them: the key's second column,
/// where Ensembl VEP lists one, left out of the column line and of each row.
fn tab_listed_once(reference: &Reference) -> Vec<String> {
    match second_listing(reference.tab[0], '\t', reference.key) {
        Some(place) => reference
            .tab
            .iter()
            .map(|line| without_cell(line, '\t', place))
            .collect(),
        None => reference.tab.iter().map(|line| line.to_string()).collect(),
    }
}

/// The reference's `CSQ` format and record line as vep-rs writes them: the
/// key's second field, where Ensembl VEP lists one, left out of the format and
/// of each allele's entry.
fn vcf_listed_once(reference: &Reference) -> Vec<String> {
    let [format, record] = reference.vcf;
    let Some(place) = second_listing(format, '|', reference.key) else {
        return vec![format.to_string(), record.to_string()];
    };
    let (columns, entries) = record.split_once("CSQ=").expect("the CSQ key");
    let entries: Vec<String> = entries
        .split(',')
        .map(|entry| without_cell(entry, '|', place))
        .collect();
    vec![
        without_cell(format, '|', place),
        format!("{columns}CSQ={}", entries.join(",")),
    ]
}

/// Ensembl VEP keys the line hash by the whole string, so a key a source
/// declares holds the source's values (OutputFactory.pm:1300-1304) wherever
/// the flag field of the same name sets none: with the flag off, and under
/// the flag on a row the field has no value for (`SYMBOL` on a row without a
/// transcript, :1529-1534). The default, tab and VCF outputs print those
/// values on the row of the allele the source's record names and nothing on
/// the other allele's.
#[test]
fn a_custom_key_named_as_a_flag_column_prints_the_sources_value_where_the_flag_field_has_none() {
    let tmp = tempfile::tempdir().unwrap();
    for reference in [&VARIANT_CLASS_OFF, &MAX_AF_OFF, &SYMBOL_OFF, &SYMBOL_ON] {
        let case = format!("{} {:?}", reference.key, reference.flags);

        let text = fs::read_to_string(annotate(tmp.path(), "default", reference)).unwrap();
        let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(rows, reference.default, "default format, {case}");

        let text = fs::read_to_string(annotate(tmp.path(), "tab", reference)).unwrap();
        let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with("##")).collect();
        assert_eq!(rows, tab_listed_once(reference), "tab format, {case}");

        let text = fs::read_to_string(annotate(tmp.path(), "vcf", reference)).unwrap();
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
        assert_eq!(
            vec![format, record],
            vcf_listed_once(reference),
            "VCF format, {case}"
        );
    }
}

/// The JSON factory moves the custom annotations to the record before any row
/// hash is built (OutputFactory/JSON.pm:249-258), so the JSON output prints
/// the source's record under `custom_annotations` and the key on no entry,
/// under `--symbol` too.
#[test]
fn a_custom_key_named_as_a_flag_column_stays_out_of_the_json_entries() {
    let tmp = tempfile::tempdir().unwrap();
    let expected: serde_json::Value = serde_json::from_str(SYMBOL_JSON).unwrap();
    for reference in [&SYMBOL_OFF, &SYMBOL_ON] {
        let text = fs::read_to_string(annotate(tmp.path(), "json", reference)).unwrap();
        let mut record: serde_json::Value =
            serde_json::from_str(text.lines().next().expect("the record's line")).unwrap();
        let object = record.as_object_mut().expect("a JSON object");
        object.remove("input");
        object.remove("vep_rs_version");
        assert_eq!(record, expected, "{:?}", reference.flags);
    }
}

/// The helpers that drop Ensembl VEP's second listing leave a reference that
/// lists the key once as it is and take one cell out of each line of one that
/// lists it twice.
#[test]
fn the_second_listing_is_dropped_from_the_column_line_and_from_each_row() {
    assert_eq!(tab_listed_once(&VARIANT_CLASS_OFF), VARIANT_CLASS_OFF.tab);
    assert_eq!(vcf_listed_once(&MAX_AF_OFF), MAX_AF_OFF.vcf);

    let tab = tab_listed_once(&SYMBOL_ON);
    assert_eq!(tab[0].split('\t').count(), 23);
    assert_eq!(tab[0].split('\t').filter(|c| *c == "SYMBOL").count(), 1);
    assert_eq!(tab[0].split('\t').nth(17), Some("SYMBOL"));
    assert_eq!(tab[1].split('\t').nth(17), Some("exact_2"));
    assert_eq!(tab[2].split('\t').nth(17), Some("-"));
    for line in &tab {
        assert_eq!(line.split('\t').count(), 23);
    }

    let vcf = vcf_listed_once(&SYMBOL_OFF);
    assert_eq!(vcf[0].split('|').filter(|f| *f == "SYMBOL").count(), 1);
    assert_eq!(vcf[0].split('|').count(), 26);
    let entries = vcf[1].split_once("CSQ=").unwrap().1;
    for (entry, symbol) in entries.split(',').zip(["exact_2", ""]) {
        assert_eq!(entry.split('|').count(), 26);
        assert_eq!(entry.split('|').nth(3), Some(symbol));
    }
}

/// The Parquet output stores the key as a string that holds the `--tab` cell:
/// the source's value on the row of the allele its record names, null on the
/// other allele's.
#[cfg(feature = "parquet")]
#[test]
fn a_custom_key_named_as_a_flag_column_is_stored_in_parquet_as_the_tab_cell() {
    let tmp = tempfile::tempdir().unwrap();
    for (reference, cell) in [
        (&VARIANT_CLASS_OFF, "exact"),
        (&MAX_AF_OFF, "0.004"),
        (&SYMBOL_OFF, "exact_2"),
        (&SYMBOL_ON, "exact_2"),
    ] {
        let tab = tab_listed_once(reference);
        let key = tab[0]
            .trim_start_matches('#')
            .split('\t')
            .position(|column| column == reference.key)
            .expect("the key's column");
        assert_eq!(tab[1].split('\t').nth(key), Some(cell));
        let out = annotate(tmp.path(), "parquet", reference);
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
        assert_eq!(
            cells,
            ["-", cell],
            "{} {:?}",
            reference.key,
            reference.flags
        );
    }
}
