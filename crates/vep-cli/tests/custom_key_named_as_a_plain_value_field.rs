// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! A `--custom` key that spells a flag field Ensembl VEP sets as a plain
//! value before it pushes the source's values under the same key, on the rows
//! the binary writes for one record: line 9 of the custom corpus's input, for
//! whose site the corpus's `custom/sites.vcf.gz` holds a record on one of the
//! two alleles. A source named `VARIANT` that asks for `CLASS` declares
//! `VARIANT_CLASS`, and sources named `IMPACT`, `ALLELE_NUM` and `PICK` list
//! their record names under those keys. Ensembl VEP 116.3 ends its run at the
//! first row that reaches such a push (`OutputFactory.pm:1301`, `:1303`), so
//! the rows held here are derived: the field's value, then the source's
//! values Ensembl VEP pushes after it. Where Ensembl VEP completes, under an
//! `overlap` source and a field the row assigns and in JSON, the rows are the
//! ones it writes, its second listing of the key left out (the documented
//! difference `duplicate_output_column_collapsed`).

use std::fs;
use std::path::{Path, PathBuf};

mod common;

/// The 1-based line of the record in the custom corpus's `variants.vcf`.
const RECORD_LINE: usize = 9;

/// The thirteen columns of the record's two rows, which no case changes.
const ROW_OF_THE_NAMED_ALLELE: &str =
    "rs1444031311\t21:5380280\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-";
const ROW_OF_THE_OTHER_ALLELE: &str =
    "rs1444031311\t21:5380280\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-";

/// One run over the record: the source's short name and type, the flags, the
/// key that spells a flag field and the cell the tab format prints under it on
/// the row of the allele the source's record names and on the other row.
struct Case {
    short_name: &'static str,
    source_type: &'static str,
    flags: &'static [&'static str],
    key: &'static str,
    tab_cells: [&'static str; 2],
}

/// `VARIANT_CLASS` under `--variant_class`: a value of the record's hash
/// (`OutputFactory.pm:967`), which the record of the row's allele is pushed
/// after (`:1426-1435`).
const VARIANT_CLASS_EXACT: Case = Case {
    short_name: "VARIANT",
    source_type: "exact",
    flags: &["--variant_class"],
    key: "VARIANT_CLASS",
    tab_cells: ["SNV,exact", "SNV"],
};

/// `VARIANT_CLASS` under `--variant_class` and an `overlap` source, whose
/// record belongs to the whole input record and is pushed after the class on
/// the record's hash (`:1004-1014`), so every row carries both.
const VARIANT_CLASS_OVERLAP: Case = Case {
    short_name: "VARIANT",
    source_type: "overlap",
    flags: &["--variant_class"],
    key: "VARIANT_CLASS",
    tab_cells: ["SNV,exact", "SNV,exact"],
};

/// `IMPACT`, which every row assigns (`:1344`) before the record of its
/// allele is pushed.
const IMPACT_EXACT: Case = Case {
    short_name: "IMPACT",
    source_type: "exact",
    flags: &[],
    key: "IMPACT",
    tab_cells: ["MODIFIER,exact_2", "MODIFIER"],
};

/// `IMPACT` under an `overlap` source: the row's assignment writes over the
/// name pushed on the record's hash, and Ensembl VEP 116.3 completes with
/// these cells.
const IMPACT_OVERLAP: Case = Case {
    short_name: "IMPACT",
    source_type: "overlap",
    flags: &[],
    key: "IMPACT",
    tab_cells: ["MODIFIER", "MODIFIER"],
};

/// `ALLELE_NUM` under `--allele_number` (`:1354`).
const ALLELE_NUM_EXACT: Case = Case {
    short_name: "ALLELE_NUM",
    source_type: "exact",
    flags: &["--allele_number"],
    key: "ALLELE_NUM",
    tab_cells: ["1,exact_2", "2"],
};

/// `PICK` under `--flag_pick`: set on the picked row alone (`:1363`), which is
/// the row of the allele the source's record names.
const PICK_EXACT: Case = Case {
    short_name: "PICK",
    source_type: "exact",
    flags: &["--flag_pick"],
    key: "PICK",
    tab_cells: ["1,exact_2", "-"],
};

/// `PICK` under `--flag_pick` and an `overlap` source: the picked row's
/// assignment writes over the name and the other row keeps it, the cells
/// Ensembl VEP 116.3 completes with.
const PICK_OVERLAP: Case = Case {
    short_name: "PICK",
    source_type: "overlap",
    flags: &["--flag_pick"],
    key: "PICK",
    tab_cells: ["1", "exact_2"],
};

/// The record Ensembl VEP 116.3 writes in JSON under a source named `VARIANT`
/// with `--variant_class`, its `input` key aside: the class stands alone under
/// `variant_class` and the source's record under `custom_annotations`.
const VARIANT_CLASS_JSON: &str = r#"{"allele_string":"A/C/G","assembly_name":"GRCh38","custom_annotations":{"VARIANT":[{"allele":"C","fields":{"AF":0.004,"CLASS":"exact","FILTER":["PASS"]},"name":"exact_2"}]},"end":5380280,"id":"rs1444031311","intergenic_consequences":[{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"C"},{"consequence_terms":["intergenic_variant"],"impact":"MODIFIER","variant_allele":"G"}],"most_severe_consequence":"intergenic_variant","seq_region_name":"21","start":5380280,"strand":1,"variant_class":"SNV"}"#;

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
/// `json` or `parquet`) under the case's flags, with the custom corpus's site
/// file as a source of the case's short name and type that asks for its `AF`
/// and `CLASS` fields.
fn annotate(dir: &Path, format: &str, case: &Case) -> PathBuf {
    let out = dir.join(format!(
        "{}.{}.{format}.out",
        case.short_name, case.source_type
    ));
    let source = format!(
        "file={},short_name={},format=vcf,type={},fields=AF%CLASS",
        corpus("GRCh38-custom")
            .join("custom/sites.vcf.gz")
            .display(),
        case.short_name,
        case.source_type
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
        .args(case.flags);
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
        "vep {format} {} {} {:?}:\n{}",
        case.short_name,
        case.source_type,
        case.flags,
        String::from_utf8_lossy(&run.stderr)
    );
    out
}

/// The cells the tab output prints under the case's key, row by row, with the
/// number of columns its column line lists under that name.
fn tab_cells_under_the_key(dir: &Path, case: &Case) -> (usize, Vec<String>) {
    let text = fs::read_to_string(annotate(dir, "tab", case)).unwrap();
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let columns: Vec<&str> = lines
        .next()
        .expect("the column line")
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let listed = columns.iter().filter(|c| **c == case.key).count();
    let place = columns
        .iter()
        .position(|c| *c == case.key)
        .expect("the key's column");
    let cells = lines
        .map(|row| {
            row.split('\t')
                .nth(place)
                .expect("the key's cell")
                .to_string()
        })
        .collect();
    (listed, cells)
}

/// Ensembl VEP sets each of these fields as a plain value and pushes the
/// source's values under the same key afterwards (`OutputFactory.pm:1300-1304`,
/// called at `:1004-1014` for the records of the whole input record and at
/// `:1426-1435` for the records of the row's allele), which ends its run. The
/// tab output holds the field's value and the pushed values after it, in one
/// column: a field of the record's hash is followed by every record of the
/// source, a field the row assigns by the records of the row's allele alone.
#[test]
fn a_custom_key_named_as_a_plain_value_field_prints_the_fields_value_then_the_sources() {
    let tmp = tempfile::tempdir().unwrap();
    for case in [
        &VARIANT_CLASS_EXACT,
        &VARIANT_CLASS_OVERLAP,
        &IMPACT_EXACT,
        &ALLELE_NUM_EXACT,
        &PICK_EXACT,
    ] {
        let (listed, cells) = tab_cells_under_the_key(tmp.path(), case);
        let name = format!("{} {} {:?}", case.key, case.source_type, case.flags);
        assert_eq!(listed, 1, "{name}: the key heads one column");
        assert_eq!(cells, case.tab_cells, "{name}");
    }
}

/// Under an `overlap` source a field the row assigns is written over the name
/// pushed on the record's hash, and Ensembl VEP 116.3 completes: `IMPACT` is
/// the row's own on both rows, and `PICK` is `1` on the picked row and the
/// source's name on the other.
#[test]
fn a_field_the_row_assigns_stands_over_the_values_of_the_whole_record() {
    let tmp = tempfile::tempdir().unwrap();
    for case in [&IMPACT_OVERLAP, &PICK_OVERLAP] {
        let (listed, cells) = tab_cells_under_the_key(tmp.path(), case);
        let name = format!("{} {} {:?}", case.key, case.source_type, case.flags);
        assert_eq!(listed, 1, "{name}: the key heads one column");
        assert_eq!(cells, case.tab_cells, "{name}");
    }
}

/// The whole rows of the default format and the `CSQ` entries of the VCF for
/// `VARIANT_CLASS` under `--variant_class` and for `IMPACT`: the default
/// format joins the two values with a comma under the flag field's own place
/// among the Extra keys, and a `CSQ` entry joins them with `&`.
#[test]
fn the_default_and_vcf_outputs_join_the_fields_value_and_the_sources() {
    let tmp = tempfile::tempdir().unwrap();

    let text = fs::read_to_string(annotate(tmp.path(), "default", &VARIANT_CLASS_EXACT)).unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        rows,
        [
            format!(
                "{ROW_OF_THE_NAMED_ALLELE}\tIMPACT=MODIFIER;VARIANT_CLASS=SNV,exact;VARIANT=exact_2;VARIANT_AF=0.004;VARIANT_FILTER=PASS"
            ),
            format!("{ROW_OF_THE_OTHER_ALLELE}\tIMPACT=MODIFIER;VARIANT_CLASS=SNV"),
        ]
    );

    let text = fs::read_to_string(annotate(tmp.path(), "default", &IMPACT_EXACT)).unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        rows,
        [
            format!(
                "{ROW_OF_THE_NAMED_ALLELE}\tIMPACT=MODIFIER,exact_2;IMPACT_AF=0.004;IMPACT_CLASS=exact;IMPACT_FILTER=PASS"
            ),
            format!("{ROW_OF_THE_OTHER_ALLELE}\tIMPACT=MODIFIER"),
        ]
    );

    let vcf_record = |case: &Case| {
        let text = fs::read_to_string(annotate(tmp.path(), "vcf", case)).unwrap();
        text.lines()
            .find(|l| !l.starts_with('#'))
            .expect("the record's line")
            .to_string()
    };
    assert_eq!(
        vcf_record(&VARIANT_CLASS_EXACT),
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER|||||||||||||||||||SNV&exact||||exact_2|0.004,G|intergenic_variant|MODIFIER|||||||||||||||||||SNV|||||"
    );
    assert_eq!(
        vcf_record(&IMPACT_EXACT),
        "21\t5380280\trs1444031311\tA\tC,G\t.\t.\tCSQ=C|intergenic_variant|MODIFIER&exact_2||||||||||||||||||||||0.004|exact,G|intergenic_variant|MODIFIER|||||||||||||||||||||||"
    );
}

/// The JSON factory moves the custom annotations to the record before any
/// hash is built (`OutputFactory/JSON.pm:249-258`), so Ensembl VEP 116.3
/// completes in JSON: the class stands alone under `variant_class` and the
/// source's record under `custom_annotations`.
#[test]
fn the_json_output_keeps_the_fields_value_alone() {
    let tmp = tempfile::tempdir().unwrap();
    let expected: serde_json::Value = serde_json::from_str(VARIANT_CLASS_JSON).unwrap();
    let text = fs::read_to_string(annotate(tmp.path(), "json", &VARIANT_CLASS_EXACT)).unwrap();
    let mut record: serde_json::Value =
        serde_json::from_str(text.lines().next().expect("the record's line")).unwrap();
    let object = record.as_object_mut().expect("a JSON object");
    object.remove("input");
    object.remove("vep_rs_version");
    assert_eq!(record, expected);
}

/// The Parquet output stores the key as a string that holds the `--tab` cell:
/// the field's value and the source's after it on the row of the allele the
/// source's record names, the field's value on the other row.
#[cfg(feature = "parquet")]
#[test]
fn the_parquet_output_stores_the_tab_cell() {
    let tmp = tempfile::tempdir().unwrap();
    for case in [&VARIANT_CLASS_EXACT, &IMPACT_EXACT] {
        let text = fs::read_to_string(annotate(tmp.path(), "tab", case)).unwrap();
        let place = text
            .lines()
            .find(|l| l.starts_with("#Uploaded_variation"))
            .expect("the column line")
            .trim_start_matches('#')
            .split('\t')
            .position(|column| column == case.key)
            .expect("the key's column");
        let out = annotate(tmp.path(), "parquet", case);
        let mut cells: Vec<String> = vep_io::output::parquet::tab_rows(&out)
            .unwrap()
            .iter()
            .map(|row| {
                row.split('\t')
                    .nth(place)
                    .expect("the key's cell")
                    .to_string()
            })
            .collect();
        cells.sort_unstable();
        let mut expected = case.tab_cells.to_vec();
        expected.sort_unstable();
        assert_eq!(cells, expected, "{}", case.key);
    }
}
