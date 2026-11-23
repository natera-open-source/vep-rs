// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The documented divergences beside the parsers, the comparison and the command line: the
//! `regulatory_feature_rows` class and selector, the recoder format's refusals, and the text
//! classes that are not named.

mod support;

use indexmap::IndexMap;
use std::collections::BTreeSet;
use std::path::Path;
use std::path::PathBuf;
use support::member as at;
use support::report as report_at;
use vep_formats_compare::classes;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::compare_files::compare_files;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::input::input_records;
use vep_formats_compare::input::read_text;
use vep_formats_compare::input::text_lines;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::header_key;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::recoder::recoder_value;
use vep_formats_compare::parse::reference_text_masked;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::vcf::csq_format_fields;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::DumpOptions;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::regulatory::without_regulatory_rows;
use vep_formats_compare::rules::Side;

const HEADER_14: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";

/// The rows of the reference output of `expected_phylop-phastcons` that carry the misspelt key.
const PHYLOP_ROWS: i64 = 197;

/// The file `name` of `directory`, holding `text`.
fn write(directory: &Path, name: &str, text: &str) -> PathBuf {
    let path = directory.join(name);
    std::fs::write(&path, text).expect("the file is written");
    path
}

fn json(text: &str) -> Value {
    pyjson::loads(text).expect("the JSON reads")
}

fn text(text: &str) -> Value {
    Value::Str(text.to_string())
}

fn masks(names: &[&str]) -> BTreeSet<String> {
    names.iter().map(|name| name.to_string()).collect()
}

/// A `Documented` of no manifest whose mask set names `classes`.
fn documented_with(classes: &[&str]) -> Documented {
    Documented::from_manifest(None, &masks(classes), None, None, None).expect("no manifest reads")
}

fn data_lines(text: &str) -> Vec<String> {
    let lines = text_lines(text).into_iter();
    lines.filter(|line| !line.starts_with('#')).collect()
}

fn header_lines(text: &str) -> Vec<String> {
    let lines = text_lines(text).into_iter();
    lines.filter(|line| line.starts_with('#')).collect()
}

fn rebuild(text: &str, rows: &[String]) -> String {
    let mut lines = header_lines(text);
    lines.extend(rows.iter().cloned());
    lines.join("\n") + "\n"
}

/// `compare_entries` with every option at its default.
fn compared(expected: &[Entry], actual: &[Entry]) -> Value {
    let side = Side::default();
    let options = CompareOptions::new(&side);
    compare_entries(expected, actual, &mut Documented::new(), &options)
        .expect("the entries compare")
}

/// `compare_files` of two default-format files with the classes `classes` and nothing else.
fn files_compared(expected: &Path, actual: &Path, classes: &[&str]) -> Value {
    compare_files(
        "default",
        expected,
        actual,
        None,
        &mut documented_with(classes),
        &masks(classes),
        &IndexMap::new(),
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

fn keys(entries: &[Entry]) -> Vec<&str> {
    entries.iter().map(|entry| entry.key.as_str()).collect()
}

/// The exit code, the standard output and the standard error of the tool on `arguments`.
fn tool(arguments: &[&str]) -> (i32, String, String) {
    let ran = support::run(Path::new(env!("CARGO_MANIFEST_DIR")), arguments);
    (ran.code, ran.out(), ran.err())
}

#[test]
fn a_recoder_value_is_its_canonical_text_and_a_masked_entry_is_counted() {
    let recoder_masks: Vec<&str> = classes::tables()
        .recoder_masks()
        .iter()
        .map(|(class, _)| *class)
        .collect();
    let mut documented = documented_with(&recoder_masks);
    let mut value = |key: &str, value: &str| recoder_value(key, &json(value), &mut documented);
    assert_eq!(
        value("hgvsc", r#"["b", "a", "b", "NM_1.1:c.1A>G"]"#).as_deref(),
        Some(r#"["a","b"]"#)
    );
    assert_eq!(
        value("var_synonyms", r#"["ClinVar: VCV2,RCV1"]"#).as_deref(),
        Some(r#"["ClinVar: RCV1,VCV2"]"#)
    );
    assert_eq!(value("input", r#""rs1""#).as_deref(), Some(r#""rs1""#));
    assert_eq!(value("strand", "[1, -1]").as_deref(), Some("[-1,1]"));
    assert_eq!(
        documented.applied.iter().collect::<Vec<_>>(),
        [("refseq_transcript_entries", 1)]
    );
}

/// The `regulatory_feature_rows` mask over parsed rows and over two files; its warnings-mode
/// half is in `tests/warnings.rs`.
#[test]
fn the_regulatory_mask_drops_the_references_regulatory_rows() {
    let corpus = support::golden().join("GRCh38-release116");
    let read =
        |name: &str| read_text(&corpus.join("expected").join(name)).expect("the corpus reads");
    let (default, vcf, json_text) = (
        read("default.txt.gz"),
        read("vcf.vcf.gz"),
        read("json.jsonl.gz"),
    );
    let inputs = input_records(&corpus.join("variants.vcf")).expect("the input reads");

    let regulatory_row = "bnd_intergenic_pair_a\t21:11000001\tN.\t-\tENSR00000347123\tRegulatoryFeature\tregulatory_region_variant\
        \t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;VARIANT_CLASS=SNV;BIOTYPE=promoter";
    let motif_row = "bnd_intergenic_pair_a\t21:11000001\tN.\t-\tENSM00000577720\tMotifFeature\tTF_binding_site_variant\
        \t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;MOTIF_NAME=ENSPFM0086;MOTIF_POS=3;HIGH_INF_POS=N;TRANSCRIPTION_FACTORS=EGR1";
    let rows = data_lines(&default);
    let with_rows = |extra: &[&str]| {
        let mut planted = vec![rows[0].clone()];
        planted.extend(extra.iter().map(|row| row.to_string()));
        planted.extend(rows[1..].iter().cloned());
        rebuild(&default, &planted)
    };
    let mut reference =
        parse_default(&with_rows(&[regulatory_row, motif_row])).expect("the rows parse");
    let actual = parse_default(&default).expect("the rows parse");
    assert_eq!(
        at(
            &compared(&reference, &actual),
            &["one_sided", "expected_only"]
        ),
        &Value::Int(2)
    );
    let kept = without_regulatory_rows(&mut reference).expect("the rows filter");
    assert_eq!(keys(&kept), keys(&actual));
    assert_eq!(at(&compared(&kept, &actual), &["verdict"]), &text("green"));

    let mut vcf_rows = data_lines(&vcf);
    let headers = header_lines(&vcf);
    let csq_header = headers
        .iter()
        .find(|line| line.starts_with("##INFO=<ID=CSQ,"))
        .expect("the VCF names its CSQ fields");
    let values = [
        ("Allele", "N."),
        ("Consequence", "regulatory_region_variant"),
        ("IMPACT", "MODIFIER"),
        ("Feature_type", "RegulatoryFeature"),
        ("Feature", "ENSR00000347123"),
        ("BIOTYPE", "promoter"),
    ];
    let chunk: Vec<&str> = csq_format_fields(csq_header)
        .expect("the header names the fields")
        .iter()
        .map(|name| {
            let held = values.iter().find(|(field, _)| *field == name.as_str());
            held.map_or("", |(_, value)| value)
        })
        .collect();
    vcf_rows[0] = vcf_rows[0].replace(";CSQ=", &format!(";CSQ={},", chunk.join("|")));
    let mut vcf_reference = parse_vcf(&rebuild(&vcf, &vcf_rows), &inputs).expect("the VCF parses");
    let vcf_actual = parse_vcf(&vcf, &inputs).expect("the VCF parses");
    assert_eq!(
        at(
            &compared(&vcf_reference, &vcf_actual),
            &["one_sided", "expected_only"]
        ),
        &Value::Int(1)
    );
    let vcf_kept = without_regulatory_rows(&mut vcf_reference).expect("the rows filter");
    assert_eq!(
        at(&compared(&vcf_kept, &vcf_actual), &["verdict"]),
        &text("green")
    );

    let json_lines: Vec<String> = text_lines(&json_text)
        .into_iter()
        .filter(|line| !line.is_empty())
        .collect();
    let with_members = |line: &str, members: &[(&str, Value)]| {
        let mut object = json(line).into_dict().expect("a line is an object");
        for (name, value) in members {
            object.insert(name.to_string(), value.clone());
        }
        pyjson::dumps(&Value::Dict(object), &DumpOptions::default())
    };
    let first = with_members(
        &json_lines[0],
        &[
            (
                "regulatory_feature_consequences",
                json(
                    r#"[{"variant_allele": "N.", "regulatory_feature_id": "ENSR00000347123",
                         "consequence_terms": ["regulatory_region_variant"], "impact": "MODIFIER", "biotype": "promoter"}]"#,
                ),
            ),
            ("most_severe_consequence", text("regulatory_region_variant")),
        ],
    );
    let json_of = |lines: &[String]| lines.join("\n") + "\n";
    let mut planted = vec![first.clone()];
    planted.extend(json_lines[1..].iter().cloned());
    let mut json_reference = parse_json(&json_of(&planted), &inputs).expect("the JSON parses");
    let json_actual = parse_json(&json_text, &inputs).expect("the JSON parses");
    let unmasked = compared(&json_reference, &json_actual);
    assert_eq!(
        at(&unmasked, &["one_sided", "expected_only"]),
        &Value::Int(1)
    );
    assert_eq!(
        at(
            &unmasked,
            &["columns", "most_severe_consequence", "differing"]
        ),
        &Value::Int(1)
    );
    let kept = without_regulatory_rows(&mut json_reference).expect("the rows filter");
    assert_eq!(
        kept[0].fields["most_severe_consequence"],
        "intergenic_variant"
    );
    assert!(!kept[0]
        .fields
        .contains_key("regulatory_feature_consequences"));
    assert_eq!(
        at(&compared(&kept, &json_actual), &["verdict"]),
        &text("green")
    );
    // A record that lost no row keeps the reference's term, so a wrong one stays red.
    let second = with_members(
        &json_lines[1],
        &[("most_severe_consequence", text("stop_gained"))],
    );
    let mut planted = vec![first, second];
    planted.extend(json_lines[2..].iter().cloned());
    let mut wrong = parse_json(&json_of(&planted), &inputs).expect("the JSON parses");
    let kept = without_regulatory_rows(&mut wrong).expect("the rows filter");
    assert_eq!(
        at(
            &compared(&kept, &json_actual),
            &["columns", "most_severe_consequence", "differing"]
        ),
        &Value::Int(1)
    );

    let scratch = support::scratch("manifest-regulatory-rows");
    let expected_path = write(&scratch, "reference.txt", &with_rows(&[regulatory_row]));
    let actual_path = write(&scratch, "actual.txt", &default);
    let report = files_compared(&expected_path, &actual_path, &["regulatory_feature_rows"]);
    assert_eq!(at(&report, &["verdict"]), &text("green"));
    assert_eq!(
        at(&report, &["masks"]),
        &json(
            r#"{"requested": ["regulatory_feature_rows"], "applied": {"regulatory_feature_rows": 1},
                "unused": [], "available": {}}"#
        )
    );
    let report = files_compared(&actual_path, &actual_path, &["regulatory_feature_rows"]);
    assert_eq!(at(&report, &["verdict"]), &text("green"));
    assert_eq!(at(&report, &["masks", "applied"]), &json("{}"));
    assert_eq!(
        at(&report, &["masks", "unused"]),
        &json(r#"["regulatory_feature_rows"]"#)
    );
}

#[test]
fn text_classes_that_are_not_named_leave_the_headers_alone() {
    let custom = support::golden().join("GRCh38-custom");
    let read = |case: &str| {
        read_text(&custom.join(case).join("default.txt.gz")).expect("the corpus reads")
    };
    let reference = read("expected_phylop-phastcons");
    let pc_reference = read("expected_bed-overlap-stats");
    let named: [&[&str]; 3] = [
        &[],
        &[classes::ORDER_ONLY, classes::REGULATORY_FEATURE_ROWS],
        &[classes::PC_DESCRIPTION_TEXT],
    ];
    for names in named {
        let mut documented = documented_with(names);
        assert_eq!(
            reference_text_masked(&reference, &mut documented),
            reference
        );
        assert!(!documented.applied.contains(classes::PC_DESCRIPTION_TEXT));
    }
    let named: [&[&str]; 2] = [&[], &[classes::CONSERVATION_SHORT_NAME_SPELLING]];
    for names in named {
        let mut documented = documented_with(names);
        assert_eq!(
            reference_text_masked(&pc_reference, &mut documented),
            pc_reference
        );
        assert!(!documented
            .applied
            .contains(classes::CONSERVATION_SHORT_NAME_SPELLING));
    }
    let scratch = support::scratch("manifest-text-classes");
    let reference_path = write(&scratch, "reference.txt", &reference);
    let candidate_path = write(
        &scratch,
        "candidate.txt",
        &reference.replace("phlyoP100way", "phyloP100way"),
    );
    let other = [classes::PC_DESCRIPTION_TEXT, classes::ORDER_ONLY];
    let report = files_compared(&reference_path, &candidate_path, &other);
    assert_eq!(at(&report, &["verdict"]), &text("red"));
    assert_eq!(
        at(&report, &["columns", "phlyoP100way", "differing"]),
        &Value::Int(PHYLOP_ROWS)
    );
    assert_eq!(at(&report, &["masks", "applied"]), &json("{}"));
    assert_eq!(
        at(&report, &["masks", "unused"]),
        &json(r#"["order_only", "pc_description_text"]"#)
    );
    assert_eq!(
        header_key("## regions_PC : text").as_deref(),
        Some("regions_PC")
    );
    assert_eq!(
        header_key(r#"##INFO=<ID=regions_PC,Number=.,Type=String,Description="text">"#).as_deref(),
        Some("regions_PC")
    );
    assert_eq!(header_key("## Column descriptions:"), None);
    assert_eq!(header_key("#Uploaded_variation\tLocation"), None);
    assert_eq!(header_key("21\t100\t.\tA\tG"), None);
}

/// One planted row of the regulatory pair: the record, the feature, its type, its term
/// and its motif name.
type PairRow = (
    usize,
    &'static str,
    &'static str,
    &'static str,
    &'static str,
);

/// Mixed records, an unrelated transcript-only record and a motif-only record, as the
/// reference (`actual` false) or the candidate writes them in `fmt`.
fn regulatory_pair_text(fmt: &str, inputs: &[String], actual: bool) -> String {
    let rows: [PairRow; 5] = [
        (
            0,
            "ENSR1",
            "RegulatoryFeature",
            "regulatory_region_variant",
            "-",
        ),
        (
            0,
            "ENSM1",
            "MotifFeature",
            "TF_binding_site_variant",
            "matrix1",
        ),
        (0, "ENST1", "Transcript", "missense_variant", "-"),
        (1, "ENST2", "Transcript", "intron_variant", "-"),
        (
            2,
            "ENSM2",
            "MotifFeature",
            "TF_binding_site_variant",
            "matrix1",
        ),
    ];
    let selected: Vec<PairRow> = rows
        .into_iter()
        .filter(|row| !(actual && row.0 == 1))
        .map(|(record, feature, kind, term, matrix)| {
            let term = if actual && kind == "Transcript" {
                "stop_gained"
            } else {
                term
            };
            (record, feature, kind, term, matrix)
        })
        .collect();
    if fmt == "default" || fmt == "tab" {
        let lines: Vec<String> = selected
            .iter()
            .map(|(record, feature, kind, term, matrix)| {
                format!("1:{}\tT\t{feature}\t{kind}\t{term}\t{matrix}", 100 + record)
            })
            .collect();
        return format!(
            "#Location\tAllele\tFeature\tFeature_type\tConsequence\tMOTIF_NAME\n{}\n",
            lines.join("\n")
        );
    }
    if fmt == "vcf" {
        let mut lines = Vec::new();
        for (record, input) in inputs.iter().enumerate() {
            let chunks: Vec<String> = selected
                .iter()
                .filter(|row| row.0 == record)
                .map(|(_, feature, kind, term, matrix)| {
                    format!("T|{feature}|{kind}|{term}|{matrix}")
                })
                .collect();
            if !chunks.is_empty() {
                let (columns, _) = input.rsplit_once('\t').expect("an input line has columns");
                lines.push(format!("{columns}\tCSQ={}", chunks.join(",")));
            }
        }
        return format!(
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: Allele|Feature|Feature_type|Consequence|MOTIF_NAME\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}\n",
            lines.join("\n")
        );
    }
    let mut objects = Vec::new();
    for (record, input) in inputs.iter().enumerate() {
        let own: Vec<&PairRow> = selected.iter().filter(|row| row.0 == record).collect();
        let Some(first) = own.first() else {
            continue;
        };
        let most_severe = match (record, actual) {
            (0, true) => "stop_gained",
            (0, false) => "missense_variant",
            _ => first.3,
        };
        let mut object: IndexMap<String, Value> = IndexMap::from([
            ("input".to_string(), text(input)),
            ("id".to_string(), text(&format!("v{record}"))),
            ("most_severe_consequence".to_string(), text(most_severe)),
        ]);
        for (_, feature, kind, term, matrix) in own {
            let kind = match *kind {
                "RegulatoryFeature" => "regulatory_feature",
                "MotifFeature" => "motif_feature",
                _ => "transcript",
            };
            let mut item: IndexMap<String, Value> = IndexMap::from([
                ("variant_allele".to_string(), text("T")),
                (format!("{kind}_id"), text(feature)),
                (
                    "consequence_terms".to_string(),
                    Value::List(vec![text(term)]),
                ),
            ]);
            if *matrix != "-" {
                item.insert("motif_name".to_string(), text(matrix));
            }
            let list = object
                .entry(format!("{kind}_consequences"))
                .or_insert_with(|| Value::List(Vec::new()));
            if let Value::List(items) = list {
                items.push(Value::Dict(item));
            }
        }
        objects.push(pyjson::dumps(&Value::Dict(object), &DumpOptions::default()));
    }
    objects.join("\n") + "\n"
}

/// The selector keeps the regulatory rows of both sides, in the format `fmt`.
fn keep_only_compares_the_regulatory_rows(fmt: &str) {
    let scratch = support::scratch(&format!("manifest-keep-only-{fmt}"));
    let inputs: Vec<String> = (0..3)
        .map(|i| format!("1\t{}\tv{i}\tA\tT\t.\t.\t.", 100 + i))
        .collect();
    let input = write(&scratch, "input.vcf", &(inputs.join("\n") + "\n"));
    let expected = write(
        &scratch,
        "expected.txt",
        &regulatory_pair_text(fmt, &inputs, false),
    );
    let original = regulatory_pair_text(fmt, &inputs, true);
    let actual = write(&scratch, "actual.txt", &original);
    let report_path = scratch.join("report.json");
    let paths =
        [&expected, &actual, &input, &report_path].map(|path| path.to_str().expect("a text path"));
    let plain = [
        "--format",
        fmt,
        "--expected",
        paths[0],
        "--actual",
        paths[1],
        "--input",
        paths[2],
        "--report",
        paths[3],
    ];
    let mut selected = plain.to_vec();
    selected.extend(["--keep-only", "regulatory_feature_rows"]);

    assert_eq!(tool(&plain).0, 1);
    assert_eq!(tool(&selected).0, 0);
    let report = report_at(&report_path);
    assert_eq!(
        at(&report, &["keep_only"]),
        &text("regulatory_feature_rows")
    );
    let compared = Value::Int(if fmt == "vcf" || fmt == "json" { 5 } else { 3 });
    for count in ["compared", "expected_entries", "actual_entries"] {
        assert_eq!(at(&report, &["counts", count]), &compared, "{count}");
    }
    assert_eq!(
        at(&report, &["one_sided"]),
        &json(r#"{"expected_only": 0, "actual_only": 0}"#)
    );
    assert_eq!(at(&report, &["masks", "requested"]), &json("[]"));
    // A motif-field disagreement is reported under the selector.
    std::fs::write(&actual, original.replace("matrix1", "matrix2")).expect("the file is written");
    assert_eq!(tool(&selected).0, 1);
    let motif_name = if fmt == "json" {
        "motif_name"
    } else {
        "MOTIF_NAME"
    };
    assert_eq!(
        at(
            &report_at(&report_path),
            &["columns", motif_name, "differing"]
        ),
        &Value::Int(2)
    );
    // One-sided regulatory feature identities must not be masked out.
    std::fs::write(&actual, original.replace("ENSR1", "ENSR2")).expect("the file is written");
    assert_eq!(tool(&selected).0, 1);
    assert_eq!(
        at(&report_at(&report_path), &["one_sided"]),
        &json(r#"{"expected_only": 1, "actual_only": 1}"#)
    );
}

#[test]
fn keep_only_compares_the_regulatory_rows_of_the_default_format() {
    keep_only_compares_the_regulatory_rows("default");
}

#[test]
fn keep_only_compares_the_regulatory_rows_of_the_tab_format() {
    keep_only_compares_the_regulatory_rows("tab");
}

#[test]
fn keep_only_compares_the_regulatory_rows_of_the_vcf_format() {
    keep_only_compares_the_regulatory_rows("vcf");
}

#[test]
fn keep_only_compares_the_regulatory_rows_of_the_json_format() {
    keep_only_compares_the_regulatory_rows("json");
}

#[test]
fn keep_only_is_refused_in_the_warnings_mode() {
    let (code, _, error) = tool(&["--warnings", "--keep-only", "regulatory_feature_rows"]);
    assert_eq!(code, 2);
    assert!(error.contains("--keep-only is available only in --format comparison mode"));
}

#[test]
fn keep_only_is_refused_in_the_self_parity_mode() {
    let (code, _, error) = tool(&["--self-parity", "--keep-only", "regulatory_feature_rows"]);
    assert_eq!(code, 2);
    assert!(error.contains("--keep-only is available only in --format comparison mode"));
}

#[test]
fn keep_only_is_refused_beside_the_mask_that_drops_the_same_rows() {
    let scratch = support::scratch("manifest-keep-only-and-mask");
    let unused = scratch.join("unused");
    let unused = unused.to_str().expect("a text path");
    let (code, _, error) = tool(&[
        "--format",
        "default",
        "--expected",
        unused,
        "--actual",
        unused,
        "--keep-only",
        "regulatory_feature_rows",
        "--masks",
        "regulatory_feature_rows",
    ]);
    assert_eq!(code, 2);
    assert!(error.contains("conflicts with --masks regulatory_feature_rows"));
}

/// A reference recoder output: one allele object, a `null` and a warnings object.
const RECODER_EXPECTED: &str = r#"[{"T": {"vcf_string": ["21-25585733-C-T"], "hgvsg": ["NC_000021.9:g.25585733C>T"], "hgvsc": ["ENST00000352957.9:c.991G>A", "NM_017446.4:c.991G>A", "ENST00000307301.12:c.*18G>A"], "hgvsp": ["NP_059142.3:p.Ala331Thr", "ENSP00000284967.7:p.Ala331Thr", "LRG_1p1:p.Ala331Thr"], "ga4gh_vrs": [{"state": {"sequence": "T", "type": "SequenceState"}, "location": {"interval": {"end": 25585733, "start": "25585732", "type": "SimpleInterval"}, "sequence_id": "refseq:NC_000021.9", "type": "SequenceLocation"}, "type": "Allele"}], "input": "rs142513484", "spdi": ["NC_000021.9:25585732:C:T"], "mane_select": [{"hgvsc": "ENST00000352957.9:c.991G>A", "hgvsg": "NC_000021.9:g.25585733C>T", "hgvsp": "ENSP00000284967.7:p.Ala331Thr"}, {"hgvsg": "NC_000021.9:g.25585733C>T", "hgvsc": "NM_017446.4:c.991G>A", "hgvsp": "NP_059142.3:p.Ala331Thr"}], "var_synonyms": ["ClinVar: VCV000829541,RCV003317418", "OMIM: 606870.0005"], "id": ["rs142513484"]}}, null, {"warnings": ["No variant found with ID 'rs_nope'"]}]"#;

#[test]
fn the_recoder_format_refuses_what_is_not_an_array_of_objects() {
    let scratch = support::scratch("manifest-recoder-refusals");
    let indented = DumpOptions {
        indent: Some(3),
        ..DumpOptions::default()
    };
    let expected = write(
        &scratch,
        "e.json",
        &pyjson::dumps(&json(RECODER_EXPECTED), &indented),
    );
    let expected = expected.to_str().expect("a text path");
    let against = |actual: &Path, more: &[&str]| {
        let mut arguments = vec![
            "--format",
            "recoder",
            "--expected",
            expected,
            "--actual",
            actual.to_str().expect("a text path"),
        ];
        arguments.extend(more);
        let (code, _, error) = tool(&arguments);
        (code, error)
    };
    let (code, error) = against(&write(&scratch, "object.json", r#"{"T": {}}"#), &[]);
    assert_eq!(code, 2);
    assert!(error.contains("the recoder output is not a JSON array"));
    let (code, error) = against(&write(&scratch, "scalar.json", "[1]"), &[]);
    assert_eq!(code, 2);
    assert!(error.contains("recoder object 1 is neither an object nor null"));
    for option in ["--keep-only", "--masks"] {
        let (code, error) = against(Path::new(expected), &[option, "regulatory_feature_rows"]);
        assert_eq!(code, 2, "{option}");
        assert!(error.contains("no regulatory rows"), "{option}");
    }
}

fn default_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    consequence: &str,
    extra: &str,
    feature_type: &str,
) -> String {
    format!("{name}\t{location}\t{allele}\tENSG1\t{feature}\t{feature_type}\t{consequence}\t-\t-\t-\t-\t-\t-\t{extra}")
}

/// Two transcript rows for each of eight records, each record at a Location of its own.
fn planted_rows() -> Vec<String> {
    let mut rows = Vec::new();
    for record in 0..8 {
        for feature in ["ENST1", "ENST2"] {
            rows.push(default_row(
                &format!("rec{record}"),
                &format!("21:{}", 1000 + record),
                &"ACGT"[record % 4..record % 4 + 1],
                feature,
                "intron_variant",
                &format!("IMPACT=MODIFIER;STRAND=1;DOMAINS=db:a{record},db:b{record}"),
                "Transcript",
            ));
        }
    }
    rows
}

/// The report, the standard output and the exit code of a comparison whose reference holds a
/// regulatory row inside a record's group of transcript rows.
#[test]
fn a_group_holding_a_regulatory_reference_row_is_compared_without_it() {
    let scratch = support::scratch("manifest-regulatory-group");
    let rows = planted_rows();
    let regulatory = default_row(
        "rec2",
        "21:1002",
        "G",
        "ENSR1",
        "regulatory_region_variant",
        "IMPACT=MODIFIER",
        "RegulatoryFeature",
    );
    let mut reference = rows[..6].to_vec();
    reference.push(regulatory);
    reference.extend(rows[6..].iter().cloned());
    let default_text = |rows: &[String]| format!("{HEADER_14}\n{}\n", rows.join("\n"));
    let with_row = write(&scratch, "with_row.txt", &default_text(&reference));
    let without_row = write(&scratch, "without_row.txt", &default_text(&rows));
    let report_path = scratch.join("report.json");
    let run = |actual: &Path, more: &[&str]| {
        let mut arguments = vec![
            "--format",
            "default",
            "--expected",
            with_row.to_str().expect("a text path"),
            "--actual",
            actual.to_str().expect("a text path"),
            "--report",
            report_path.to_str().expect("a text path"),
        ];
        arguments.extend(more);
        let (code, out, _) = tool(&arguments);
        (out, code)
    };
    let masked = ["--masks", "regulatory_feature_rows"];
    assert_eq!(
        run(&without_row, &masked),
        ("default: green\n".to_string(), 0)
    );
    assert_eq!(
        at(&report_at(&report_path), &["masks", "applied"]),
        &json(r#"{"regulatory_feature_rows": 1}"#)
    );
    assert_eq!(
        run(&with_row, &masked),
        (
            "default: red\n  unexpected in the candidate: 1 entries\n".to_string(),
            1
        )
    );
    assert_eq!(run(&with_row, &[]), ("default: green\n".to_string(), 0));
}
