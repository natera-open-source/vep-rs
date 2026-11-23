// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The comparison of two output files: the classes read from the documented masks, the
//! floors of each format and the regulatory rows kept alone.

mod compare_common;

use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::path::Path;

use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::classes::CONSERVATION_SHORT_NAME_SPELLING;
use vep_formats_compare::classes::ORDER_ONLY;
use vep_formats_compare::compare::FLOOR_NOT_APPLICABLE;
use vep_formats_compare::compare_files::compare_files_through;
use vep_formats_compare::compare_files::FLOOR_UNDECLARED_NOT_APPLICABLE;
use vep_formats_compare::consts::RECORD_LEVEL_COLUMNS;
use vep_formats_compare::fields::populated_counts;
use vep_formats_compare::floors::json_floor_fields;
use vep_formats_compare::floors::FLOOR_CUSTOM_SOURCE_NOT_APPLICABLE;
use vep_formats_compare::input::input_records;
use vep_formats_compare::input::read_text;
use vep_formats_compare::input::text_lines;
use vep_formats_compare::lists::declared_columns;
use vep_formats_compare::lists::populated_under;
use vep_formats_compare::lists::same_items;
use vep_formats_compare::lists::JSON_ARRAY;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::reference_text_masked;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::vcf::csq_format_fields;
use vep_formats_compare::parse::vcf::header_columns;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::regulatory::only_regulatory_rows;
use vep_formats_compare::rules::Side;

/// The report of two files compared whole in `fmt`, four examples a group.
fn files(
    fmt: &str,
    reference: &Path,
    candidate: &Path,
    inputs: Option<&[String]>,
    documented: &mut Documented,
    masks: &BTreeSet<String>,
    floors: &[(&str, i64)],
) -> Value {
    let floors: IndexMap<String, i64> = floors
        .iter()
        .map(|(name, minimum)| (name.to_string(), *minimum))
        .collect();
    compare_files_through(
        None,
        fmt,
        reference,
        candidate,
        inputs,
        documented,
        masks,
        &floors,
        &IndexMap::new(),
        4,
        false,
        None,
        None,
        None,
        &mut Side::default(),
    )
    .expect("the files compare")
}

/// A reference compared with itself in `fmt` under `floors`.
fn compare_pair(
    fmt: &str,
    path: &Path,
    floors: &[(&str, i64)],
    inputs: Option<&[String]>,
) -> Value {
    files(
        fmt,
        path,
        path,
        inputs,
        &mut no_documented(),
        &BTreeSet::new(),
        floors,
    )
}

fn masks_member(requested: &[&str], applied: &[(&str, i64)], unused: &[&str]) -> Value {
    let list = |items: &[&str]| {
        Value::List(
            items
                .iter()
                .map(|item| Value::Str(item.to_string()))
                .collect(),
        )
    };
    let counted: IndexMap<String, Value> = applied
        .iter()
        .map(|(name, held)| (name.to_string(), Value::Int(*held)))
        .collect();
    Value::Dict(
        [
            ("requested", list(requested)),
            ("applied", Value::Dict(counted)),
            ("unused", list(unused)),
            ("available", Value::Dict(IndexMap::new())),
        ]
        .into_iter()
        .map(|(name, value)| (name.to_string(), value))
        .collect(),
    )
}

/// The frequency-filter reference with a second allele's frequency after the first on the
/// reference's side and before it on the candidate's, and the rows carrying the list.
fn two_allele_frequencies(directory: &Path) -> (std::path::PathBuf, std::path::PathBuf, i64) {
    let freqfilter = golden().join("GRCh38-freqfilter");
    let default = read_text(
        &freqfilter
            .join("expected_freqfilter")
            .join("default.txt.gz"),
    )
    .expect("the output reads");
    let one = "FREQS=AF:C:0.0007987";
    assert!(default.contains(one));
    let reference = planted(
        directory,
        "reference.txt",
        &default.replace(one, "FREQS=AF:C:0.0007987,AF:T:0.5"),
    );
    let candidate = planted(
        directory,
        "candidate.txt",
        &default.replace(one, "FREQS=AF:T:0.5,AF:C:0.0007987"),
    );
    let rows = i64::try_from(default.matches(one).count()).expect("a count");
    (reference, candidate, rows)
}

/// FREQS is printed in hash order (OutputFactory.pm:1152-1153), so a two-allele record's
/// entries agree as a multiset under the class; `compare_files` applies the class from
/// the documented masks and reports it applied or unused, with the other formats untouched.
#[test]
fn the_order_only_class_covers_frequencies_and_is_reported_for_two_files() {
    let directory = scratch("order_only_frequencies");
    let (reference, candidate, rows) = two_allele_frequencies(&directory);
    let none = BTreeSet::new();
    let plain = files(
        "default",
        &reference,
        &candidate,
        None,
        &mut no_documented(),
        &none,
        &[],
    );
    assert_eq!(verdict(&plain), "red");
    assert_eq!(number(&plain, &["columns", "FREQS", "differing"]), rows);
    let order_only = masks(&[ORDER_ONLY]);
    let report = files(
        "default",
        &reference,
        &candidate,
        None,
        &mut documented_with(&[ORDER_ONLY]),
        &order_only,
        &[],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(number(&report, &["columns", "FREQS", "order_only"]), rows);
    assert_eq!(number(&report, &["columns", "Extra", "order_only"]), rows);
    assert_eq!(
        at(&report, &["masks"]),
        &masks_member(&[ORDER_ONLY], &[(ORDER_ONLY, 2 * rows)], &[])
    );
    // The class is unused on a candidate that prints every list in the reference's order.
    let same = files(
        "default",
        &reference,
        &reference,
        None,
        &mut documented_with(&[ORDER_ONLY]),
        &order_only,
        &[],
    );
    assert_eq!(verdict(&same), "green");
    assert_eq!(
        at(&same, &["masks"]),
        &masks_member(&[ORDER_ONLY], &[], &[ORDER_ONLY])
    );
    assert_eq!(number(&same, &["columns", "FREQS", "order_only"]), 0);
    // The record-level FREQS has no JSON key, so the JSON format compares no such field.
    assert!(RECORD_LEVEL_COLUMNS.contains(&"FREQS"));
    assert!(same_items("a,b", "b,a", ",") && !same_items("a,b", "a,a,b", ","));
    assert!(same_items("-", "-", ",") && !same_items("-", "a", ","));
    assert!(same_items(
        r#"{"db":"x","name":"1"},{"db":"y","name":"2"}"#,
        r#"{"db":"y","name":"2"},{"db":"x","name":"1"}"#,
        JSON_ARRAY
    ));
    assert!(!same_items("a,b", "b,a", JSON_ARRAY));
}

/// `compare_files` applies and reports the classes of `documented.masks`; its `masks`
/// argument is the set the caller built `documented` from and is not merged into it.
#[test]
fn two_files_are_compared_under_the_documented_masks_alone() {
    let directory = scratch("documented_masks_alone");
    let (reference, candidate, rows) = two_allele_frequencies(&directory);
    let mut documented = documented_with(&[ORDER_ONLY]);
    let report = files(
        "default",
        &reference,
        &candidate,
        None,
        &mut documented,
        &BTreeSet::new(),
        &[],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        at(&report, &["masks"]),
        &masks_member(&[ORDER_ONLY], &[(ORDER_ONLY, 2 * rows)], &[])
    );

    let mut documented = no_documented();
    let report = files(
        "default",
        &reference,
        &candidate,
        None,
        &mut documented,
        &masks(&[ORDER_ONLY]),
        &[],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "FREQS", "differing"]), rows);
    assert_eq!(at(&report, &["masks"]), &masks_member(&[], &[], &[]));
    assert!(documented.masks.is_empty());
}

/// Rows of the custom corpus's phylop-phastcons reference, those carrying the --phyloP 100
/// source's value in the default format (its header names the key phlyoP100way, Config.pm:555),
/// and the record objects carrying it in the JSON format.
const CUSTOM_ROWS: i64 = 1586;
const PHYLOP_ROWS: i64 = 197;
const PHYLOP_JSON_RECORDS: i64 = 21;
/// Rows of the custom corpus's VCF-source references carrying the `sites` source.
const CUSTOM_SITES_ROWS: i64 = 454;

const FORMAT_FILES: [(&str, &str); 4] = [
    ("default", "default.txt"),
    ("tab", "tab.txt"),
    ("vcf", "vcf.vcf"),
    ("json", "json.jsonl"),
];

/// The rewrites of the misspelt key the class counts on the phylop-phastcons reference of
/// `fmt`: the default format carries the key in its header line and on every annotated row,
/// the JSON on every annotated record, the tab and VCF formats in their header line and
/// column list alone.
fn rewrites(fmt: &str) -> i64 {
    match fmt {
        "default" => PHYLOP_ROWS + 1,
        "json" => PHYLOP_JSON_RECORDS,
        _ => 2,
    }
}

/// Ensembl VEP's --phyloP set names its source `phlyoP<N>way` (Config.pm:555) and the
/// candidate `phyloP<N>way`, so the key heads a column only one side has in every format;
/// under the class the reference's key is read as the candidate's and the two agree, the
/// rewrite counted once per occurrence (header lines and rows alike), and without the class
/// the column is red on both sides. The candidate texts are the committed reference with the
/// key spelt as the candidate spells it, since the two engines agree on every value.
#[test]
fn the_conservation_spelling_class_renames_the_misspelt_column() {
    let directory = scratch("conservation_spelling_column");
    let custom = golden().join("GRCh38-custom");
    let class = [CONSERVATION_SHORT_NAME_SPELLING];
    let inputs = input_records(&custom.join("variants.vcf")).expect("the input reads");
    for (fmt, name) in FORMAT_FILES {
        let reference = read_text(
            &custom
                .join("expected_phylop-phastcons")
                .join(format!("{name}.gz")),
        )
        .expect("the output reads");
        assert!(reference.contains("phlyoP100way"));
        // The header names the file; JSON has none.
        assert_eq!(reference.contains("phyloP100way.bw"), fmt != "json");
        let candidate = reference.replace("phlyoP100way", "phyloP100way");
        let reference_path = planted(&directory, &format!("reference.{name}"), &reference);
        let candidate_path = planted(&directory, &format!("candidate.{name}"), &candidate);
        let plain = files(
            fmt,
            &reference_path,
            &candidate_path,
            Some(&inputs),
            &mut no_documented(),
            &BTreeSet::new(),
            &[],
        );
        assert_eq!(verdict(&plain), "red", "{fmt}");
        if fmt == "json" {
            // The JSON format keys the source under the record's custom_annotations object.
            assert_eq!(
                number(&plain, &["columns", "custom_annotations", "differing"]),
                PHYLOP_JSON_RECORDS
            );
        } else {
            // The default format carries the key on the rows the source annotates; a tab
            // column and a CSQ field are declared on every row, empty where it does not.
            let differing = if fmt == "default" {
                PHYLOP_ROWS
            } else {
                CUSTOM_ROWS
            };
            for column in ["phlyoP100way", "phyloP100way"] {
                assert_eq!(number(&plain, &["columns", column, "differing"]), differing);
            }
        }
        let report = files(
            fmt,
            &reference_path,
            &candidate_path,
            Some(&inputs),
            &mut documented_with(&class),
            &masks(&class),
            &[],
        );
        assert_eq!(verdict(&report), "green", "{fmt} {:?}", reasons(&report));
        assert!(at(&report, &["columns"]).get("phlyoP100way").is_none());
        let applied = number(
            &report,
            &["masks", "applied", CONSERVATION_SHORT_NAME_SPELLING],
        );
        let written = i64::try_from(reference.matches("phlyoP100way").count()).expect("a count");
        assert_eq!(applied, written);
        assert_eq!(applied, rewrites(fmt));
        assert!(texts(at(&report, &["masks", "unused"])).is_empty());
        if fmt != "json" {
            assert_eq!(
                number(&report, &["columns", "phyloP100way", "differing"]),
                0
            );
        }
    }
    // The file's name is spelt correctly on both sides and stays as it is; a key of another
    // width renames with its own number.
    let mut documented = documented_with(&class);
    assert_eq!(
        reference_text_masked(
            "## phlyoP100way : [PATH]/hg38.phyloP100way.bw\n",
            &mut documented
        ),
        "## phyloP100way : [PATH]/hg38.phyloP100way.bw\n"
    );
    assert_eq!(
        reference_text_masked("phlyoP20way=1.5;phlyoP470way=-0.2", &mut documented),
        "phyloP20way=1.5;phyloP470way=-0.2"
    );
    assert_eq!(documented.applied.get(CONSERVATION_SHORT_NAME_SPELLING), 3);
}

/// Per default-format column of a corpus reference with a JSON floor mapping, the
/// default format's populated count beside the JSON format's floor count.
fn json_floor_counts_beside_the_default(corpus: &Path) -> BTreeMap<String, (i64, i64)> {
    let inputs = input_records(&corpus.join("variants.vcf")).expect("the input reads");
    let reference_of = |name: &str| -> String {
        read_text(&corpus.join("expected").join(name)).expect("the output reads")
    };
    let default = parse_default(&reference_of("default.txt.gz")).expect("the reference parses");
    let objects =
        parse_json(&reference_of("json.jsonl.gz"), &inputs).expect("the reference parses");
    let default_counts = populated_counts(&default);
    let columns: BTreeSet<&str> = default
        .iter()
        .flat_map(|e| e.fields.keys().map(String::as_str))
        .collect();
    let mut out = BTreeMap::new();
    for column in columns {
        if let Some(counted) = json_floor_fields(column) {
            let counted: Vec<&str> = counted.iter().map(String::as_str).collect();
            out.insert(
                column.to_string(),
                (
                    default_counts.get(column),
                    populated_under(&objects, &counted),
                ),
            );
        }
    }
    out
}

/// On the GRCh38 and GRCh38-everything references of Ensembl VEP 116.3, every default
/// column with a JSON floor mapping counts the same under its JSON fields as in the default
/// format: the position keys survive a `?` start (JSON.pm:340-341 keep the end), the
/// record objects' own `strand` is not counted, a `domains` array left empty by the
/// re-split still marks a populated row, and a present key is a populated column.
#[test]
fn the_json_floor_counts_equal_the_default_counts_on_every_mapped_column() {
    let corpus = golden().join("GRCh38");
    let grch38 = json_floor_counts_beside_the_default(&corpus);
    assert_eq!(grch38.get("cDNA_position"), Some(&(8879, 8879)));
    assert_eq!(grch38.get("CDS_position"), Some(&(5907, 5907)));
    assert_eq!(grch38.get("STRAND"), Some(&(20621, 20621)));
    let mismatched: Vec<_> = grch38.iter().filter(|(_, held)| held.0 != held.1).collect();
    assert!(mismatched.is_empty(), "{mismatched:?}");
    // A floor spelled with the JSON key itself is counted the same way, over the
    // consequence entries alone: the record objects carry their own strand.
    let inputs = input_records(&corpus.join("variants.vcf")).expect("the input reads");
    let json_path = corpus.join("expected").join("json.jsonl.gz");
    let report = compare_pair("json", &json_path, &[("strand", 20621)], Some(&inputs));
    assert_eq!(
        texts(at(&report, &["columns", "strand", "floor_fields"])),
        ["strand"]
    );
    assert_eq!(
        number(&report, &["columns", "strand", "expected_populated"]),
        20621
    );
    assert_eq!(verdict(&report), "green");
    let everything = json_floor_counts_beside_the_default(&golden().join("GRCh38-everything"));
    assert_eq!(everything.len(), 40);
    assert_eq!(everything.get("DOMAINS"), Some(&(2055, 2055)));
    assert_eq!(everything.get("Protein_position"), Some(&(2734, 2734)));
    let mismatched: Vec<_> = everything
        .iter()
        .filter(|(_, held)| held.0 != held.1)
        .collect();
    assert!(mismatched.is_empty(), "{mismatched:?}");
}

/// A `--custom` source's columns have no row count in the JSON format, which carries a
/// record's custom annotations once, on the record object under `custom_annotations`
/// (JSON.pm:249-258), so a floor on `<source>`, `<source>_<field>` or `<source>_FILTER`
/// is not applicable there, with its own note; the default format counts every column, the
/// VCF source's FILTER included, which it carries in Extra alone: the tab header and the CSQ
/// list name the `fields` the option lists, so a floor on `sites_FILTER` is not applicable
/// in those two formats while the other three columns are counted. A column one header
/// declares and the other does not is counted and stays a difference; a column no entry of
/// either side carries in the default or JSON format stays an unmet floor.
#[test]
fn the_floors_on_custom_source_columns_follow_each_formats_header() {
    let directory = scratch("custom_source_floors");
    let custom = golden().join("GRCh38-custom");
    let inputs = input_records(&custom.join("variants.vcf")).expect("the input reads");
    let inputs = Some(inputs.as_slice());
    let sites = ["sites", "sites_AF", "sites_CLASS", "sites_FILTER"];
    let floors: Vec<(&str, i64)> = sites
        .iter()
        .map(|name| (*name, CUSTOM_SITES_ROWS))
        .collect();
    let exact = custom.join("expected_vcf-exact");
    let default = compare_pair("default", &exact.join("default.txt.gz"), &floors, None);
    assert_eq!(verdict(&default), "green");
    assert!(texts(at(&default, &["floors_not_applicable"])).is_empty());
    for name in sites {
        assert_eq!(
            number(&default, &["columns", name, "expected_populated"]),
            CUSTOM_SITES_ROWS
        );
    }
    for (fmt, name) in [("tab", "tab.txt.gz"), ("vcf", "vcf.vcf.gz")] {
        let report = compare_pair(fmt, &exact.join(name), &floors, inputs);
        assert_eq!(verdict(&report), "green", "{fmt} {:?}", reasons(&report));
        assert_eq!(
            texts(at(&report, &["floors_not_applicable"])),
            ["sites_FILTER"]
        );
        assert_eq!(
            at(&report, &["columns", "sites_FILTER", "floor_note"]).as_str(),
            Some(FLOOR_UNDECLARED_NOT_APPLICABLE)
        );
        assert_eq!(
            at(&report, &["columns", "sites_FILTER", "floor_met"]),
            &Value::Null
        );
        for column in ["sites", "sites_AF", "sites_CLASS"] {
            assert_eq!(
                at(&report, &["columns", column, "floor_met"]),
                &Value::Bool(true)
            );
            assert_eq!(
                number(&report, &["columns", column, "expected_populated"]),
                CUSTOM_SITES_ROWS
            );
            assert!(at(&report, &["columns", column])
                .get("floor_fields")
                .is_none());
        }
    }
    let json_report = compare_pair("json", &exact.join("json.jsonl.gz"), &floors, inputs);
    assert_eq!(verdict(&json_report), "green");
    assert_eq!(texts(at(&json_report, &["floors_not_applicable"])), sites);
    for column in ["sites", "sites_CLASS", "sites_FILTER"] {
        assert_eq!(
            at(&json_report, &["columns", column, "floor_note"]).as_str(),
            Some(FLOOR_CUSTOM_SOURCE_NOT_APPLICABLE)
        );
    }
    // An AF column is record-level by name before it is a source's field.
    assert_eq!(
        at(&json_report, &["columns", "sites_AF", "floor_note"]).as_str(),
        Some(FLOOR_NOT_APPLICABLE)
    );
    assert_eq!(
        number(
            &json_report,
            &["columns", "custom_annotations", "differing"]
        ),
        0
    );
    // The BED, GFF and bigWig sources, their statistics and the shortcut's file-name key alike.
    let sources: [(&str, &[&str]); 4] = [
        (
            "expected_bed-overlap-stats",
            &["regions", "regions_PC", "regions_min", "regions_count"],
        ),
        (
            "expected_gff-promoter",
            &[
                "GENCODE_PROMOTER",
                "GENCODE_PROMOTER_type",
                "GENCODE_PROMOTER_associated_gene",
            ],
        ),
        (
            "expected_bigwig-overlap-stats",
            &["cons", "cons_max", "cons_sum"],
        ),
        ("expected_bigwig-shortcut", &["scores.bw"]),
    ];
    for (expected_dir, columns) in sources {
        let floors: Vec<(&str, i64)> = columns.iter().map(|c| (*c, 1)).collect();
        let report = compare_pair(
            "json",
            &custom.join(expected_dir).join("json.jsonl.gz"),
            &floors,
            inputs,
        );
        assert_eq!(
            verdict(&report),
            "green",
            "{expected_dir} {:?}",
            reasons(&report)
        );
        let mut sorted = columns.to_vec();
        sorted.sort_unstable();
        assert_eq!(texts(at(&report, &["floors_not_applicable"])), sorted);
        for column in columns {
            assert_eq!(
                at(&report, &["columns", column, "floor_note"]).as_str(),
                Some(FLOOR_CUSTOM_SOURCE_NOT_APPLICABLE)
            );
        }
    }
    // A column no source and no entry spells is an unmet floor in the JSON format too.
    let nothing = [("nothing", 1)];
    let unmet = compare_pair("json", &exact.join("json.jsonl.gz"), &nothing, inputs);
    assert_eq!(verdict(&unmet), "red");
    assert_eq!(
        reasons(&unmet),
        ["floor nothing=1: expected populated 0, actual populated 0"]
    );
    let unmet = compare_pair("default", &exact.join("default.txt.gz"), &nothing, None);
    assert_eq!(verdict(&unmet), "red");
    // The candidate's tab lacks the sites_CLASS column the reference declares: the floor is
    // counted on both sides and the column is a difference on every row, never not applicable.
    let tab = read_text(&exact.join("tab.txt.gz")).expect("the output reads");
    let mut rows: Vec<String> = Vec::new();
    let mut dropped = 0;
    for line in text_lines(&tab) {
        if line.starts_with("##") {
            rows.push(line);
            continue;
        }
        let mut cells: Vec<&str> = line.split('\t').collect();
        if line.starts_with('#') {
            dropped = cells
                .iter()
                .position(|cell| *cell == "sites_CLASS")
                .expect("the header declares the column");
        }
        cells.remove(dropped);
        rows.push(cells.join("\t"));
    }
    let reference_path = planted(&directory, "reference.tab", &tab);
    let candidate_path = planted(&directory, "candidate.tab", &(rows.join("\n") + "\n"));
    let one_sided = files(
        "tab",
        &reference_path,
        &candidate_path,
        inputs,
        &mut no_documented(),
        &BTreeSet::new(),
        &floors,
    );
    assert_eq!(verdict(&one_sided), "red");
    assert_eq!(
        texts(at(&one_sided, &["floors_not_applicable"])),
        ["sites_FILTER"]
    );
    assert_eq!(
        number(&one_sided, &["columns", "sites_CLASS", "differing"]),
        CUSTOM_ROWS
    );
    assert_eq!(
        at(&one_sided, &["columns", "sites_CLASS", "floor_met"]),
        &Value::Bool(false)
    );
    let reason = format!(
        "floor sites_CLASS={CUSTOM_SITES_ROWS}: expected populated {CUSTOM_SITES_ROWS}, actual populated 0"
    );
    assert!(reasons(&one_sided).contains(&reason.as_str()));
    // The header decides, not the rows: a tab output with its header alone, and a VCF whose
    // records carry no CSQ, still declare the source's columns, so their floors stay unmet.
    let header: Vec<String> = text_lines(&tab)
        .into_iter()
        .filter(|line| line.starts_with('#'))
        .collect();
    let header_only = planted(&directory, "header_only.tab", &(header.join("\n") + "\n"));
    let set_of = |names: Vec<String>| -> BTreeSet<String> { names.into_iter().collect() };
    assert_eq!(
        declared_columns("tab", &header_only).expect("the header reads"),
        set_of(header_columns(&tab))
    );
    let two_floors = [("sites", 1), ("sites_FILTER", 1)];
    let report = compare_pair("tab", &header_only, &two_floors, None);
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        reasons(&report),
        ["floor sites=1: expected populated 0, actual populated 0"]
    );
    assert_eq!(
        texts(at(&report, &["floors_not_applicable"])),
        ["sites_FILTER"]
    );
    let vcf = read_text(&exact.join("vcf.vcf.gz")).expect("the output reads");
    let csq = regex::Regex::new(r";?CSQ=[^\t\n]*").expect("the pattern compiles");
    let stripped = csq.replace_all(&vcf, "").into_owned();
    let without_csq = planted(&directory, "without_csq.vcf", &stripped);
    assert_eq!(text_lines(&stripped).len(), text_lines(&vcf).len());
    let records = stripped
        .split_once("#CHROM")
        .expect("the output has its column line")
        .1;
    assert!(!records.contains("CSQ="));
    let declared = declared_columns("vcf", &without_csq).expect("the header reads");
    assert!(declared.contains("sites_CLASS") && !declared.contains("sites_FILTER"));
    let report = compare_pair("vcf", &without_csq, &two_floors, inputs);
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        reasons(&report),
        ["floor sites=1: expected populated 0, actual populated 0"]
    );
    let default_text = read_text(&exact.join("default.txt.gz")).expect("the output reads");
    assert_eq!(
        declared_columns("tab", &exact.join("default.txt.gz")).expect("the header reads"),
        set_of(header_columns(&default_text))
    );
    assert!(declared_columns("vcf", &exact.join("tab.txt.gz"))
        .expect("the header reads")
        .is_empty());
    // A --keep_csq output carries two CSQ header lines, the kept list first; the entries are
    // read by the last one, and so are the declared columns.
    let keep_csq = golden()
        .join("GRCh38-writers")
        .join("expected_keep-csq")
        .join("vcf.vcf.gz");
    let keep_csq_text = read_text(&keep_csq).expect("the output reads");
    let csq_lines: Vec<String> = text_lines(&keep_csq_text)
        .into_iter()
        .filter(|line| line.starts_with("##INFO=<ID=CSQ"))
        .collect();
    assert_eq!(csq_lines.len(), 2);
    assert_eq!(
        declared_columns("vcf", &keep_csq).expect("the header reads"),
        set_of(csq_format_fields(&csq_lines[1]).expect("the line names its fields"))
    );
    let two_lists = planted(
        &directory,
        "two_lists.vcf",
        &keep_csq_text.replacen(
            &csq_lines[0],
            &csq_lines[0].replace("|SYMBOL_SOURCE|HGNC_ID", ""),
            1,
        ),
    );
    let declared = declared_columns("vcf", &two_lists).expect("the header reads");
    assert!(declared.contains("SYMBOL_SOURCE") && declared.contains("HGNC_ID"));
}

/// A case's floor on the --phyloP source names the column in either spelling: under the
/// class the reference's key is read as `phyloP<N>way`, and a floor spelt `phlyoP<N>way` is
/// read under the renamed key, so the default, tab and VCF formats count the 197 annotated
/// rows on both sides and report no column under the misspelt name; the JSON format carries
/// the source once per record under `custom_annotations`, so its floor is not applicable
/// there and the renamed objects agree. Without the class the misspelt floor counts the
/// reference's rows alone and is red.
#[test]
fn the_conservation_spelling_class_reads_a_floor_under_the_renamed_key() {
    let directory = scratch("conservation_spelling_floor");
    let custom = golden().join("GRCh38-custom");
    let class = [CONSERVATION_SHORT_NAME_SPELLING];
    let inputs = input_records(&custom.join("variants.vcf")).expect("the input reads");
    for (fmt, name) in FORMAT_FILES {
        let reference = read_text(
            &custom
                .join("expected_phylop-phastcons")
                .join(format!("{name}.gz")),
        )
        .expect("the output reads");
        let reference_path = planted(&directory, &format!("reference.{name}"), &reference);
        let candidate_path = planted(
            &directory,
            &format!("candidate.{name}"),
            &reference.replace("phlyoP100way", "phyloP100way"),
        );
        for spelling in ["phlyoP100way", "phyloP100way"] {
            let floors = [(spelling, PHYLOP_ROWS), ("phastCons100way", PHYLOP_ROWS)];
            let report = files(
                fmt,
                &reference_path,
                &candidate_path,
                Some(&inputs),
                &mut documented_with(&class),
                &masks(&class),
                &floors,
            );
            assert_eq!(
                verdict(&report),
                "green",
                "{fmt} {spelling} {:?}",
                reasons(&report)
            );
            assert!(at(&report, &["columns"]).get("phlyoP100way").is_none());
            let mut applied = IndexMap::new();
            applied.insert(
                CONSERVATION_SHORT_NAME_SPELLING.to_string(),
                Value::Int(rewrites(fmt)),
            );
            assert_eq!(at(&report, &["masks", "applied"]), &Value::Dict(applied));
            let column = ["columns", "phyloP100way"];
            if fmt == "json" {
                assert_eq!(
                    texts(at(&report, &["floors_not_applicable"])),
                    ["phastCons100way", "phyloP100way"]
                );
                assert_eq!(
                    at(&report, &["columns", "phyloP100way", "floor_note"]).as_str(),
                    Some(FLOOR_CUSTOM_SOURCE_NOT_APPLICABLE)
                );
                assert_eq!(
                    number(&report, &["columns", "custom_annotations", "differing"]),
                    0
                );
            } else {
                assert!(texts(at(&report, &["floors_not_applicable"])).is_empty());
                let held = at(&report, &column);
                assert_eq!(held.get("floor_met"), Some(&Value::Bool(true)));
                assert_eq!(
                    held.get("expected_populated"),
                    Some(&Value::Int(PHYLOP_ROWS))
                );
                assert_eq!(held.get("actual_populated"), Some(&Value::Int(PHYLOP_ROWS)));
            }
        }
    }
    let plain = files(
        "default",
        &directory.join("reference.default.txt"),
        &directory.join("candidate.default.txt"),
        None,
        &mut no_documented(),
        &BTreeSet::new(),
        &[("phlyoP100way", PHYLOP_ROWS)],
    );
    assert_eq!(verdict(&plain), "red");
    let reason = format!(
        "floor phlyoP100way={PHYLOP_ROWS}: expected populated {PHYLOP_ROWS}, actual populated 0"
    );
    assert!(reasons(&plain).contains(&reason.as_str()));
    assert_eq!(
        number(&plain, &["columns", "phyloP100way", "differing"]),
        PHYLOP_ROWS
    );
}

#[test]
fn the_regulatory_rows_alone_do_not_hide_a_motif_field_one_side_lacks() {
    let inputs = vec!["1\t100\tv0\tA\tT\t.\t.\t.".to_string()];
    let expected = dumped(
        r#"{"input": "1\t100\tv0\tA\tT\t.\t.\t.", "motif_feature_consequences": [
            {"variant_allele": "T", "motif_feature_id": "ENSM1",
             "consequence_terms": ["TF_binding_site_variant"], "high_inf_pos": "N",
             "motif_score_change": 0.25}]}"#,
    );
    let actual = dumped(
        r#"{"input": "1\t100\tv0\tA\tT\t.\t.\t.", "motif_feature_consequences": [
            {"variant_allele": "T", "motif_feature_id": "ENSM1",
             "consequence_terms": ["TF_binding_site_variant"]}]}"#,
    );
    let regulatory = |text: &str| {
        let mut entries = parse_json(text, &inputs).expect("the text parses");
        only_regulatory_rows(&mut entries).expect("the rows are ranked")
    };
    let report = compare(
        &regulatory(&expected),
        &regulatory(&actual),
        &mut no_documented(),
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(
        number(&report, &["columns", "high_inf_pos", "differing"]),
        1
    );
    assert_eq!(
        number(&report, &["columns", "motif_score_change", "differing"]),
        1
    );
}
