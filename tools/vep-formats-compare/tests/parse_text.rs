// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Reading a file, the default and tab parsers and the reference side's text classes, on
//! the committed Ensembl VEP outputs of `tests/golden/116`. The expected values are read
//! off the committed files.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use vep_formats_compare::classes::PC_DESCRIPTION_TEXT;
use vep_formats_compare::consts::EXAMPLES_PER_GROUP;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Error;
use vep_formats_compare::input::{dedup_keys, input_records, read_text, text_lines};
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::text::{parse_default, parse_tab};
use vep_formats_compare::parse::{parse_format, reference_lines_masked, reference_text_masked};
use vep_formats_compare::pyjson::{self, DumpOptions, Value};
use vep_formats_compare::self_parity::SelfParity;

const RECORDS: usize = 36;
const ROWS: usize = 977;

fn golden() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116")
}

fn corpus() -> PathBuf {
    golden().join("GRCh38-release116")
}

/// The committed output of the corpus in the file `name` of its `expected` directory.
fn expected(name: &str) -> String {
    read_text(&corpus().join("expected").join(name)).unwrap()
}

fn data_lines(text: &str) -> Vec<String> {
    text_lines(text)
        .into_iter()
        .filter(|line| !line.starts_with('#'))
        .collect()
}

fn header_lines(text: &str) -> Vec<String> {
    text_lines(text)
        .into_iter()
        .filter(|line| line.starts_with('#'))
        .collect()
}

fn rebuild(text: &str, rows: &[String]) -> String {
    let mut lines = header_lines(text);
    lines.extend_from_slice(rows);
    lines.join("\n") + "\n"
}

/// A `Documented` of no manifest whose mask set names `class`.
fn documented_with(class: &str) -> Documented {
    let masks = BTreeSet::from([class.to_string()]);
    Documented::from_manifest(None, &masks, None, None, None).unwrap()
}

fn applied(documented: &Documented) -> Vec<(&str, i64)> {
    documented.applied.iter().collect()
}

/// The self-parity report of the four committed formats, `tab`, `vcf` and `json` each
/// replaced where given.
fn self_parity(tab: Option<&str>, vcf: Option<&str>, json: Option<&str>) -> Value {
    let default = expected("default.txt.gz");
    let tab = tab.map_or_else(|| expected("tab.txt.gz"), str::to_string);
    let vcf = vcf.map_or_else(|| expected("vcf.vcf.gz"), str::to_string);
    let json = json.map_or_else(|| expected("json.jsonl.gz"), str::to_string);
    SelfParity::new(&default, &tab, &vcf, &json, None, None, EXAMPLES_PER_GROUP)
        .unwrap()
        .run()
}

/// The member of nested objects `names` leads to.
fn at<'a>(report: &'a Value, names: &[&str]) -> &'a Value {
    names.iter().fold(report, |value, name| {
        value
            .get(name)
            .unwrap_or_else(|| panic!("no member {name}"))
    })
}

fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

/// The JSON output with its first object changed by `change`.
fn json_with_first_object(change: impl FnOnce(&mut indexmap::IndexMap<String, Value>)) -> String {
    let mut lines = text_lines(&expected("json.jsonl.gz"));
    let mut object = pyjson::loads(&lines[0]).unwrap().into_dict().unwrap();
    change(&mut object);
    lines[0] = pyjson::dumps(&Value::Dict(object), &DumpOptions::default());
    lines.join("\n") + "\n"
}

#[test]
fn exact_duplicates_of_a_key_are_counted() {
    let mut entries: Vec<Entry> = ["k", "k", "k", "other"]
        .iter()
        .map(|key| Entry {
            key: key.to_string(),
            ..Entry::default()
        })
        .collect();
    dedup_keys(&mut entries);
    let keys: Vec<&str> = entries.iter().map(|e| e.key.as_str()).collect();
    assert_eq!(keys, ["k", "k#2", "k#3", "other"]);
}

#[test]
fn the_default_parser_reads_the_committed_reference() {
    let entries = parse_default(&expected("default.txt.gz")).unwrap();
    assert_eq!(entries.len(), ROWS);
    let first = &entries[0];
    assert_eq!(first.key, "21:11000001|N.|-|-");
    assert_eq!(first.location.as_deref(), Some("21:11000001"));
    assert_eq!(first.consequence, "intergenic_variant");
    assert_eq!(first.fields["Uploaded_variation"], "bnd_intergenic_pair_a");
    assert_eq!(first.fields["Extra"], "IMPACT=MODIFIER");
    assert_eq!(first.fields["IMPACT"], "MODIFIER");
    let ablation = entries
        .iter()
        .find(|e| e.fields["Uploaded_variation"] == "del_10_6_mb" && e.feature == "ENST00000262354")
        .unwrap();
    assert_eq!(ablation.consequence, "transcript_ablation");
    assert_eq!(
        ablation.fields["Extra"],
        "IMPACT=HIGH;STRAND=-1;OverlapBP=23932;OverlapPC=100.00"
    );
    assert_eq!(ablation.fields["OverlapBP"], "23932");
    assert_eq!(ablation.fields["OverlapPC"], "100.00");
    assert_eq!(ablation.fields["STRAND"], "-1");
    let keys: BTreeSet<&str> = entries.iter().map(|e| e.key.as_str()).collect();
    assert_eq!(keys.len(), ROWS);
}

#[test]
fn the_tab_parser_reads_the_committed_reference() {
    let tab_text = expected("tab.txt.gz");
    let tab = parse_tab(&tab_text).unwrap();
    let default = parse_default(&expected("default.txt.gz")).unwrap();
    assert_eq!(tab.len(), ROWS);
    let keys =
        |entries: &[Entry]| -> Vec<String> { entries.iter().map(|e| e.key.clone()).collect() };
    assert_eq!(keys(&tab), keys(&default));
    assert!(!tab[0].fields.contains_key("Extra"));
    let columns = header_lines(&tab_text).pop().unwrap();
    let named: BTreeSet<&str> = columns[1..].split('\t').collect();
    let read: BTreeSet<&str> = tab[0].fields.keys().map(String::as_str).collect();
    assert_eq!(read, named);
    assert!(tab
        .iter()
        .zip(&default)
        .all(|(t, d)| t.fields["IMPACT"] == d.fields["IMPACT"]));
}

#[test]
fn the_vcf_and_json_formats_are_refused_without_the_input() {
    for (fmt, name) in [("vcf", "vcf.vcf.gz"), ("json", "json.jsonl.gz")] {
        let path = corpus().join("expected").join(name);
        assert_eq!(
            parse_format(fmt, &path, None, None),
            Err(Error::Comparator(format!(
                "--format {fmt} needs --input, the input VCF"
            )))
        );
    }
}

#[test]
fn the_vcf_format_is_parsed_with_the_input() {
    let inputs = input_records(&corpus().join("variants.vcf")).unwrap();
    let path = corpus().join("expected").join("vcf.vcf.gz");
    let entries = parse_format("vcf", &path, Some(&inputs), None).unwrap();
    assert_eq!(entries.len(), RECORDS + ROWS);
}

#[test]
fn self_parity_is_red_on_a_mutated_copy() {
    let tab_text = expected("tab.txt.gz");
    let mut tab_rows = data_lines(&tab_text);
    tab_rows[0] = tab_rows[0].replace("\tMODIFIER\t", "\tLOW\t");
    let report = self_parity(Some(&rebuild(&tab_text, &tab_rows)), None, None);
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    assert_eq!(
        at(&report, &["reasons"]),
        &Value::List(vec![text("tab: column IMPACT differs on 1 entries")])
    );
    assert_eq!(
        at(&report, &["pairs", "tab", "differing", "IMPACT", "count"]),
        &Value::Int(1)
    );

    let json = json_with_first_object(|object| match &mut object["intergenic_consequences"] {
        Value::List(consequences) => {
            consequences.pop();
        }
        other => panic!("no list: {other:?}"),
    });
    let report = self_parity(None, None, Some(&json));
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    let one_sided = [
        ("default_only".to_string(), Value::Int(1)),
        ("json_only".to_string(), Value::Int(0)),
    ];
    assert_eq!(
        at(&report, &["pairs", "json", "one_sided"]),
        &Value::Dict(one_sided.into_iter().collect())
    );

    let vcf_text = expected("vcf.vcf.gz");
    let mut vcf_rows = data_lines(&vcf_text);
    vcf_rows[0] = vcf_rows[0].replacen(
        "|intergenic_variant|MODIFIER|",
        "|downstream_gene_variant|MODIFIER|",
        1,
    );
    let report = self_parity(None, Some(&rebuild(&vcf_text, &vcf_rows)), None);
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    assert_eq!(
        at(
            &report,
            &["pairs", "vcf", "differing", "Consequence", "count"]
        ),
        &Value::Int(1)
    );
}

#[test]
fn self_parity_checks_the_most_severe_consequence() {
    let json = json_with_first_object(|object| {
        object.insert("most_severe_consequence".to_string(), text("stop_gained"));
    });
    let report = self_parity(None, None, Some(&json));
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    assert_eq!(
        at(&report, &["record_checks", "problems"]),
        &Value::List(vec![text(
            "bnd_intergenic_pair_a|21:11000001: most_severe_consequence 'stop_gained', the object's terms give 'intergenic_variant'"
        )])
    );
}

/// The two texts of AnnotationSourceAdaptor.pm:395-397 each name the share the other
/// setting of `reciprocal` computes (File.pm:556-564); under the class the reference's
/// `<key>_PC` header line carries the candidate's text for the share the field reports, in
/// the default and tab shape and in the VCF `##INFO` shape, every other header line and
/// every row as written. The reciprocal=1 lines are the committed bed-overlap-stats
/// reference's (its source is `reciprocal=1`); the reciprocal=0 text is the module's other
/// string.
#[test]
fn the_pc_description_class_maps_both_overlap_texts() {
    let custom = golden().join("GRCh38-custom/expected_bed-overlap-stats");
    let reciprocal_default = read_text(&custom.join("default.txt.gz")).unwrap();
    let reciprocal_vcf = read_text(&custom.join("vcf.vcf.gz")).unwrap();
    let reciprocal_text =
        "Percentage of minimum reciprocal overlap between input variant and reference variant";
    let plain_text = "Percentage of input variant covered by reference variant";
    assert!(reciprocal_default.contains(&format!(
        "## regions_PC : {reciprocal_text} from [PATH]/regions.bed.gz"
    )));
    assert!(reciprocal_vcf.contains(&format!(
        "Description=\"{reciprocal_text} from [PATH]/regions.bed.gz\">"
    )));
    let row = "21_5028997_TCA/T\t21:5028998-5028999\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tregions_PC=100.000,50.000\n";
    let text = format!(
        "## regions_PC : {reciprocal_text} from [PATH]/regions.bed.gz\n\
         ## plain_PC : {plain_text} from [PATH]/plain.bed.gz\n\
         ##INFO=<ID=regions_PC,Number=.,Type=String,Description=\"{reciprocal_text} from [PATH]/regions.bed.gz\">\n\
         ##INFO=<ID=plain_PC,Number=.,Type=String,Description=\"{plain_text} from [PATH]/plain.bed.gz\">\n\
         ## regions_min : min data from [PATH]/regions.bed.gz\n\
         ## other : {plain_text}\n\
         {row}"
    );
    let mut documented = documented_with(PC_DESCRIPTION_TEXT);
    let masked = reference_text_masked(&text, &mut documented);
    let wanted = format!(
        "## regions_PC : Percentage of the input variant covered by the annotation record from [PATH]/regions.bed.gz\n\
         ## plain_PC : Percentage of the annotation record covered by the input variant from [PATH]/plain.bed.gz\n\
         ##INFO=<ID=regions_PC,Number=.,Type=String,Description=\"Percentage of the input variant covered by the annotation record from [PATH]/regions.bed.gz\">\n\
         ##INFO=<ID=plain_PC,Number=.,Type=String,Description=\"Percentage of the annotation record covered by the input variant from [PATH]/plain.bed.gz\">\n\
         ## regions_min : min data from [PATH]/regions.bed.gz\n\
         ## other : {plain_text}\n\
         {row}"
    );
    assert_eq!(masked, wanted);
    assert_eq!(applied(&documented), [(PC_DESCRIPTION_TEXT, 4)]);
    // The lines of the text are rewritten as the text is, with the same count.
    let mut lines = text_lines(&text);
    let mut by_line = documented_with(PC_DESCRIPTION_TEXT);
    reference_lines_masked(&mut lines, &mut by_line);
    assert_eq!(lines, text_lines(&wanted));
    assert_eq!(applied(&by_line), [(PC_DESCRIPTION_TEXT, 4)]);
    // The committed reference's PC line is the one line the class rewrites in each text
    // format; the entries read no header line, so the comparison of the reference with
    // itself stays green with the class counted applied once.
    for (fmt, reference) in [("default", &reciprocal_default), ("vcf", &reciprocal_vcf)] {
        let mut documented = documented_with(PC_DESCRIPTION_TEXT);
        let masked = reference_text_masked(reference, &mut documented);
        assert_eq!(applied(&documented), [(PC_DESCRIPTION_TEXT, 1)], "{fmt}");
        assert_eq!(data_lines(&masked), data_lines(reference));
        let named = |text: &str| -> Vec<String> {
            header_lines(text)
                .into_iter()
                .filter(|line| line.contains("regions_PC"))
                .collect()
        };
        let respelt: Vec<String> = named(reference)
            .iter()
            .map(|line| {
                line.replace(
                    reciprocal_text,
                    "Percentage of the input variant covered by the annotation record",
                )
            })
            .collect();
        assert_eq!(named(&masked), respelt);
    }
}
