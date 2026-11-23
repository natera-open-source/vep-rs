// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The self-parity mode: one engine's outputs of one run, paired across the formats.

mod compare_common;

use std::collections::BTreeSet;
use std::path::Path;

use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::classes::ORDER_ONLY;
use vep_formats_compare::compare::{compare_entries, CompareOptions};
use vep_formats_compare::consts::{EXAMPLES_PER_GROUP, RECORD_LEVEL_COLUMNS};
use vep_formats_compare::entry::Entry;
use vep_formats_compare::floors::json_floor_fields;
use vep_formats_compare::input::{
    canonical_json_text, input_records, numberified_json_text, read_text, text_lines,
};
use vep_formats_compare::lists::order_only_fields;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::{json_compact, parse_json};
use vep_formats_compare::parse::text::{parse_default, parse_tab};
use vep_formats_compare::parse::vcf::{csq_format_fields, header_columns, parse_vcf};
use vep_formats_compare::pyjson::{self, DumpOptions, Value};
use vep_formats_compare::rules::Side;
use vep_formats_compare::self_parity::{
    format_coords, most_severe, perl_numeric_or_text, perl_numify, record_location,
    reverse_complement, uploaded_variation, values_equal, vcf_escape, IntOrStr, SelfParity,
};
use vep_formats_compare::shared::feature_kind;

/// The self-parity report of four texts, the input read from the JSON objects' own lines.
fn parity(default: &str, tab: &str, vcf: &str, json: &str) -> Value {
    parity_of(default, tab, vcf, json, None)
}

fn parity_of(default: &str, tab: &str, vcf: &str, json: &str, inputs: Option<&[String]>) -> Value {
    SelfParity::new(default, tab, vcf, json, None, inputs, EXAMPLES_PER_GROUP)
        .expect("the outputs parse")
        .run()
}

/// The four outputs of one run under `directory`, and the report of their self-parity with
/// the records of `input`.
fn parity_of_run(directory: &Path, input: &Path) -> Value {
    let text = |name: &str| read_text(&directory.join(name)).expect("the output reads");
    let inputs = input_records(input).expect("the input reads");
    parity_of(
        &text("default.txt.gz"),
        &text("tab.txt.gz"),
        &text("vcf.vcf.gz"),
        &text("json.jsonl.gz"),
        Some(&inputs),
    )
}

fn members(line: &str) -> IndexMap<String, Value> {
    json(line).into_dict().expect("an object")
}

fn dumped_object(members: &IndexMap<String, Value>) -> String {
    pyjson::dumps(&Value::Dict(members.clone()), &DumpOptions::default())
}

fn string(value: &str) -> Value {
    Value::Str(value.to_string())
}

/// The names of the members of the object under `path`.
fn names<'a>(report: &'a Value, path: &[&str]) -> Vec<&'a str> {
    at(report, path)
        .as_dict()
        .expect("an object")
        .keys()
        .map(String::as_str)
        .collect()
}

fn joined(lines: &[String]) -> String {
    lines.join("\n") + "\n"
}

const VRS_HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";
const VRS_ROW: &str = "rs1\t21:10538885\tC\tG1\tT1\tTranscript\tintron_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;GA4GH_VRS={};STRAND=1\n";

/// Ensembl VEP writes the `GA4GH_VRS` object with its keys in hash order, so two rows whose
/// objects agree compare equal whatever the order of their keys, in the default format's
/// `Extra` and as a tab column alike; a value that is not JSON is compared as written.
#[test]
fn ga4gh_vrs_compares_as_a_parsed_object() {
    let perl = concat!(
        r#"{"type":"Allele","state":{"sequence":"C","type":"SequenceState"},"#,
        r#""location":{"type":"SequenceLocation","sequence_id":"refseq:NC_000021.9","#,
        r#""interval":{"type":"SimpleInterval","end":10538885,"start":"10538884"}}}"#
    );
    let rust = concat!(
        r#"{"location":{"type":"SequenceLocation","interval":{"type":"SimpleInterval","#,
        r#""start":"10538884","end":10538885},"sequence_id":"refseq:NC_000021.9"},"#,
        r#""state":{"type":"SequenceState","sequence":"C"},"type":"Allele"}"#
    );
    let canonical = |text: &str| canonical_json_text(text).expect("the text is read");
    assert_eq!(canonical(perl), canonical(rust));
    for written in ["21:10538884:A:C", "-", "1e5"] {
        assert_eq!(canonical(written), written);
    }
    let row_of = |value: &str| VRS_ROW.replace("{}", value);
    let first = |text: String| parse_default(&text).expect("the text parses").remove(0);
    let perl_entry = first(format!("{VRS_HEADER}{}", row_of(perl)));
    let rust_entry = first(format!("{VRS_HEADER}{}", row_of(rust)));
    assert_eq!(
        perl_entry.fields["GA4GH_VRS"],
        rust_entry.fields["GA4GH_VRS"]
    );
    assert_eq!(perl_entry.fields["GA4GH_VRS"], canonical(perl));
    assert_eq!(perl_entry.fields["Extra"], rust_entry.fields["Extra"]);
    assert_eq!(perl_entry.fields["STRAND"], "1");
    // A pair whose value is not an object, and a bare key, are rebuilt as written.
    let bare_key_row = VRS_ROW.replace("GA4GH_VRS={};", "GA4GH_VRS;");
    let entries = parse_default(&format!("{VRS_HEADER}{}{bare_key_row}", row_of("-")))
        .expect("the text parses");
    let (dash_entry, bare_key_entry) = (&entries[0], &entries[1]);
    assert_eq!(
        dash_entry.fields["Extra"],
        "IMPACT=MODIFIER;GA4GH_VRS=-;STRAND=1"
    );
    assert_eq!(dash_entry.fields["GA4GH_VRS"], "-");
    assert_eq!(
        bare_key_entry.fields["Extra"],
        "IMPACT=MODIFIER;GA4GH_VRS;STRAND=1"
    );
    assert_eq!(bare_key_entry.fields["GA4GH_VRS"], "");
    let tab_header = "#Uploaded_variation\tLocation\tAllele\tGA4GH_VRS";
    let tab_row = |value: &str| format!("rs1\t21:10538885\tC\t{value}");
    let tab_entry = |value: &str| {
        parse_tab(&format!("{tab_header}\n{}\n", tab_row(value)))
            .expect("the text parses")
            .remove(0)
    };
    assert_eq!(
        tab_entry(perl).fields["GA4GH_VRS"],
        tab_entry(rust).fields["GA4GH_VRS"]
    );
    // Self-parity: the JSON format nests the object with both bounds numbers, where the text
    // formats carry the start as a string; the VCF format is refused with --ga4gh_vrs, so its
    // CSQ declares no such column.
    let nested = rust.replace(r#""start":"10538884""#, r#""start":10538884"#);
    assert_eq!(numberified_json_text(perl), json_compact(&json(&nested)));
    assert_eq!(numberified_json_text("21:10538884:A:C"), "21:10538884:A:C");
    // A superscript digit is no decimal digit; it stays a string.
    assert_eq!(
        numberified_json_text("{\"start\":\"\u{b2}\"}"),
        "{\"start\":\"\u{b2}\"}"
    );
    let csq = "Allele|Consequence|IMPACT|Gene|Feature_type|Feature|STRAND";
    let vcf = format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {csq}\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t10538885\trs1\tA\tC\t.\t.\tCSQ=C|intron_variant|MODIFIER|G1|Transcript|T1|1\n"
    );
    let object = |nested: &str| {
        format!(
            concat!(
                r#"{{"input": "21\t10538885\trs1\tA\tC\t.\t.\t.", "id": "rs1", "seq_region_name": "21", "#,
                r#""start": 10538885, "end": 10538885, "strand": 1, "allele_string": "A/C", "#,
                r#""most_severe_consequence": "intron_variant", "transcript_consequences": [{{"#,
                r#""variant_allele": "C", "gene_id": "G1", "transcript_id": "T1", "#,
                r#""consequence_terms": ["intron_variant"], "impact": "MODIFIER", "strand": 1, "#,
                r#""ga4gh_vrs": {}}}]}}"#,
                "\n"
            ),
            nested
        )
    };
    let default = format!("{VRS_HEADER}{}", row_of(rust));
    let tab = format!(
        "{tab_header}\tGene\tFeature\tFeature_type\tConsequence\tIMPACT\tSTRAND\n{}\tG1\tT1\tTranscript\tintron_variant\tMODIFIER\t1\n",
        tab_row(rust)
    );
    let report = parity(&default, &tab, &vcf, &object(&nested));
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert!(texts(at(&report, &["pairs", "json", "columns_compared"])).contains(&"GA4GH_VRS"));
    let moved = nested.replace(r#""end":10538885"#, r#""end":10538886"#);
    let report = parity(&default, &tab, &vcf, &object(&moved));
    assert_eq!(
        names(&report, &["pairs", "json", "differing"]),
        ["GA4GH_VRS"]
    );
}

fn committed() -> [String; 4] {
    ["default", "tab", "vcf", "json"].map(text_of)
}

#[test]
fn self_parity_is_green_on_the_committed_formats() {
    let [default, tab, vcf, json] = committed();
    let report = parity(&default, &tab, &vcf, &json);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["records"]), RECORDS as i64);
    assert_eq!(number(&report, &["default_entries"]), ROWS);
    for name in ["tab", "vcf", "json"] {
        assert_eq!(number(&report, &["pairs", name, "shared_entries"]), ROWS);
        assert_eq!(
            compact(at(&report, &["pairs", name, "one_sided"])),
            format!(r#"{{"default_only":0,"{name}_only":0}}"#)
        );
        assert!(names(&report, &["pairs", name, "differing"]).is_empty());
    }
    let compared = |name: &str| texts(at(&report, &["pairs", name, "columns_compared"]));
    assert!(compared("json").contains(&"OverlapBP"));
    assert!(!compared("tab").contains(&"OverlapBP"));
    let not_compared = |name: &str| texts(at(&report, &["not_compared", name]));
    assert_eq!(
        not_compared("vcf_implied"),
        [
            "BIOTYPE",
            "EXON",
            "HGNC_ID",
            "INTRON",
            "SYMBOL",
            "SYMBOL_SOURCE"
        ]
    );
    assert_eq!(not_compared("tab_undeclared"), ["OverlapBP", "OverlapPC"]);
    assert_eq!(not_compared("vcf_undeclared"), ["OverlapBP", "OverlapPC"]);
    assert!(not_compared("default_record_level").is_empty());
}

fn compact(value: &Value) -> String {
    json_compact(value)
}

/// Self-parity pairs the rows of a record that share a key with the default format's rows
/// in their order, so the Parquet side is read in file order, not the sorted canonical form:
/// two rows of one record and feature that differ in `IND` alone, written in reverse order,
/// pair with a default format in the same order. The two texts are the ones the adapter
/// gives for such a file: its rows in file order, and sorted.
#[test]
fn self_parity_reads_the_parquet_rows_in_file_order() {
    let header =
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tIND";
    let rows = [
        "rs1\t21:100\tA\tENSG1\tENST1\tTranscript\tintergenic_variant\tsample_b",
        "rs1\t21:100\tA\tENSG1\tENST1\tTranscript\tintergenic_variant\tsample_a",
    ];
    let mut sorted = rows;
    sorted.sort_unstable();
    let default_rows: Vec<String> = ["sample_b", "sample_a"]
        .iter()
        .map(|ind| {
            format!("rs1\t21:100\tA\tENSG1\tENST1\tTranscript\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;IND={ind}")
        })
        .collect();
    let default = default_text(&default_rows);
    let tab = format!("{header}\n{}\n", rows.join("\n"));
    let [_, _, vcf, json] = committed();
    let parquet_pair = |parquet_rows: &[&str]| {
        let parquet_text = format!("{header}\n{}\n", parquet_rows.join("\n"));
        let parity = SelfParity::new(
            &default,
            &tab,
            &vcf,
            &json,
            Some(&parquet_text),
            None,
            EXAMPLES_PER_GROUP,
        )
        .expect("the outputs parse");
        let parquet = parity.parquet.as_deref().expect("a Parquet side");
        let compare = parity._compare_tab_like(&parity.parquet_columns);
        let pair = parity._pair(
            "parquet",
            &SelfParity::_by_key(&parity.default),
            &SelfParity::_by_key(parquet),
            &compare,
        );
        drop(compare);
        pair
    };
    assert!(names(&parquet_pair(&rows), &["differing"]).is_empty());
    assert_eq!(
        number(&parquet_pair(&sorted), &["differing", "IND", "count"]),
        2
    );
}

#[test]
fn self_parity_accepts_a_name_carrying_a_bar() {
    let [default, tab, vcf, json] =
        committed().map(|text| text.replace("bnd_intergenic_pair_a", "bnd|pair_a"));
    let report = parity(&default, &tab, &vcf, &json);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["pairs", "json", "shared_entries"]), ROWS);
}

/// The self-parity report of one synthetic missense record: the default format's Extra
/// carries IMPACT, STRAND and then `extra`, the tab and VCF formats spell the same values, and
/// the JSON transcript consequence carries the members `consequence` spells beside the fixed
/// columns' members.
fn synthetic_self_parity(extra: &[(&str, &str)], consequence: &str) -> Value {
    let fixed = [
        ("Uploaded_variation", "rs1"),
        ("Location", "1:100"),
        ("Allele", "G"),
        ("Gene", "ENSG1"),
        ("Feature", "ENST1"),
        ("Feature_type", "Transcript"),
        ("Consequence", "missense_variant"),
        ("cDNA_position", "10"),
        ("CDS_position", "5"),
        ("Protein_position", "2"),
        ("Amino_acids", "A/G"),
        ("Codons", "gCa/gGa"),
        ("Existing_variation", "rs1"),
    ];
    let mut extra_column = vec![("IMPACT", "MODERATE"), ("STRAND", "1")];
    extra_column.extend_from_slice(extra);
    let tabbed = |items: &[&str]| items.join("\t");
    let fixed_names: Vec<&str> = fixed.iter().map(|(name, _)| *name).collect();
    let fixed_values: Vec<&str> = fixed.iter().map(|(_, value)| *value).collect();
    let extra_names: Vec<&str> = extra_column.iter().map(|(name, _)| *name).collect();
    let extra_values: Vec<&str> = extra_column.iter().map(|(_, value)| *value).collect();
    let pairs: Vec<String> = extra_column
        .iter()
        .map(|(name, value)| format!("{name}={value}"))
        .collect();
    let default = format!(
        "#{}\tExtra\n{}\t{}\n",
        tabbed(&fixed_names),
        tabbed(&fixed_values),
        pairs.join(";")
    );
    let tab = format!(
        "#{}\t{}\n{}\t{}\n",
        tabbed(&fixed_names),
        tabbed(&extra_names),
        tabbed(&fixed_values),
        tabbed(&extra_values)
    );
    let mut csq = vec![
        "Allele",
        "Consequence",
        "IMPACT",
        "Gene",
        "Feature_type",
        "Feature",
        "cDNA_position",
        "CDS_position",
        "Protein_position",
        "Amino_acids",
        "Codons",
        "Existing_variation",
        "STRAND",
    ];
    csq.extend(extra.iter().map(|(name, _)| *name));
    let value_of = |column: &str| {
        fixed
            .iter()
            .chain(&extra_column)
            .rev()
            .find(|(name, _)| *name == column)
            .map_or("", |(_, value)| *value)
    };
    let chunk: Vec<String> = csq
        .iter()
        .map(|column| vcf_escape(value_of(column)))
        .collect();
    let vcf = format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: {}\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\trs1\tA\tG\t.\t.\tCSQ={}\n",
        csq.join("|"),
        chunk.join("|")
    );
    let mut transcript = members(concat!(
        r#"{"variant_allele": "G", "gene_id": "ENSG1", "transcript_id": "ENST1", "#,
        r#""consequence_terms": ["missense_variant"], "impact": "MODERATE", "strand": 1, "#,
        r#""cdna_start": 10, "cdna_end": 10, "cds_start": 5, "cds_end": 5, "#,
        r#""protein_start": 2, "protein_end": 2, "amino_acids": "A/G", "codons": "gCa/gGa"}"#
    ));
    transcript.extend(members(&format!("{{{consequence}}}")));
    let mut object = members(concat!(
        r#"{"input": "1\t100\trs1\tA\tG\t.\t.\t.", "id": "rs1", "seq_region_name": "1", "#,
        r#""start": 100, "end": 100, "strand": 1, "allele_string": "A/G", "#,
        r#""most_severe_consequence": "missense_variant", "colocated_variants": [{"id": "rs1"}]}"#
    ));
    object.insert(
        "transcript_consequences".to_string(),
        Value::List(vec![Value::Dict(transcript)]),
    );
    parity(&default, &tab, &vcf, &(dumped_object(&object) + "\n"))
}

/// A synthetic record beside the committed --sift, --polyphen and --domains corpora: the
/// JSON prediction/score pair and the db/name objects must match the default format's
/// `deleterious(0.01)` and `Pfam:PF1,PANTHER:PTHR1` spellings.
#[test]
fn self_parity_undoes_the_json_sift_polyphen_and_domains_splits() {
    let extra = [
        ("SIFT", "deleterious(0.01)"),
        ("PolyPhen", "probably_damaging(0.95)"),
        ("DOMAINS", "Pfam:PF1,PANTHER:PTHR1"),
    ];
    let consequence = concat!(
        r#""sift_prediction": "deleterious", "sift_score": 0.01, "#,
        r#""polyphen_prediction": "probably_damaging", "polyphen_score": 0.95, "#,
        r#""domains": [{"db": "Pfam", "name": "PF1"}, {"db": "PANTHER", "name": "PTHR1"}]"#
    );
    let report = synthetic_self_parity(&extra, consequence);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    let compared = texts(at(&report, &["pairs", "json", "columns_compared"]));
    for column in ["SIFT", "PolyPhen", "DOMAINS"] {
        assert!(compared.contains(&column), "{column}");
    }
    let changed = format!(
        r#"{consequence}, "sift_score": 0.02, "domains": [{{"db": "Pfam", "name": "PF1"}}]"#
    );
    let report = synthetic_self_parity(&extra, &changed);
    assert_eq!(
        names(&report, &["pairs", "json", "differing"]),
        ["DOMAINS", "SIFT"]
    );
}

/// Ensembl VEP's JSON re-split (`(\w+):(\w+)`, OutputFactory/JSON.pm:363-371) keeps `$1` and
/// `$2` from the list's last successful match, so an entry the pattern does not match repeats
/// the previous object and an opening one is dropped; `numberify` then makes a digit-only name
/// a number. `Gene3D:1.10.510.10,Low_complexity_(Seg):seg,Pfam:PF1` is therefore
/// `{Gene3D 1}, {Gene3D 1}, {Pfam PF1}` and `Low_complexity_(Seg):seg,Pfam:PF1` is
/// `{Pfam PF1}`. A capture of `0` fails `push ... if $1 && $2`, so
/// `Gene3D:0.10.510.10,Low_complexity_(Seg):seg` is empty and
/// `Pfam:PF1,Gene3D:0.10.510.10,Low_complexity_(Seg):seg` is `{Pfam PF1}` alone.
#[test]
fn self_parity_domains_model_repeats_the_last_match_and_numbers_a_digit_name() {
    let domains = |text: &'static str, objects: &str| {
        synthetic_self_parity(&[("DOMAINS", text)], &format!(r#""domains": {objects}"#))
    };
    let green = |report: Value| assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    let differing = |report: Value| {
        assert_eq!(names(&report, &["pairs", "json", "differing"]), ["DOMAINS"]);
    };
    let pfam = r#"[{"db": "Pfam", "name": "PF1"}]"#;
    green(domains(
        "Gene3D:1.10.510.10,Low_complexity_(Seg):seg,Pfam:PF1",
        r#"[{"db": "Gene3D", "name": 1}, {"db": "Gene3D", "name": 1}, {"db": "Pfam", "name": "PF1"}]"#,
    ));
    green(domains("Low_complexity_(Seg):seg,Pfam:PF1", pfam));
    differing(domains(
        "Gene3D:1.10.510.10,Low_complexity_(Seg):seg,Pfam:PF1",
        r#"[{"db": "Gene3D", "name": "1"}, {"db": "Pfam", "name": "PF1"}]"#,
    ));
    green(domains("Gene3D:0.10.510.10,Low_complexity_(Seg):seg", "[]"));
    green(domains(
        "Pfam:PF1,Gene3D:0.10.510.10,Low_complexity_(Seg):seg",
        pfam,
    ));
    differing(domains(
        "Gene3D:0.10.510.10,Low_complexity_(Seg):seg",
        r#"[{"db": "Gene3D", "name": 0}, {"db": "Gene3D", "name": 0}]"#,
    ));
    assert_eq!(perl_numeric_or_text("007"), IntOrStr::Int(7));
    assert_eq!(
        perl_numeric_or_text("PF1"),
        IntOrStr::Str("PF1".to_string())
    );
    // A superscript digit is no ASCII digit: Perl's looks_like_number refuses it.
    assert_eq!(
        perl_numeric_or_text("\u{b2}"),
        IntOrStr::Str("\u{b2}".to_string())
    );
}

#[test]
fn perl_numify_follows_perl_numification() {
    assert_eq!(perl_numify("30425949"), "30425949");
    assert_eq!(perl_numify("+30425949"), "30425949");
    assert_eq!(perl_numify("30_425_949"), "30");
    assert_eq!(perl_numify("."), "0");
    assert_eq!(perl_numify("0x1D04A1D"), "0");
    assert_eq!(perl_numify("30425949.5"), "30425949.5");
    assert_eq!(perl_numify("1e999"), "Inf");
}

#[test]
fn uploaded_variation_and_location_follow_the_output_factory() {
    let named = json(concat!(
        r#"{"id": "rs1", "seq_region_name": "21", "start": 100, "end": 100, "#,
        r#""allele_string": "A/G", "input": "21\t100\trs1\tA\tG\t.\t.\t."}"#
    ));
    assert_eq!(uploaded_variation(&named), "rs1");
    assert_eq!(record_location(&named), "21:100");
    let unnamed = json(concat!(
        r#"{"id": ".", "seq_region_name": "21", "start": 10795228, "end": 10795228, "#,
        r#""allele_string": "C/G/T", "input": "21\t10795228\t.\tC\tG,T\t.\t.\t."}"#
    ));
    assert_eq!(uploaded_variation(&unnamed), "21_10795228_C/G/T");
    let insertion = json(concat!(
        r#"{"id": ".", "seq_region_name": "21", "start": 9489855, "end": 9489854, "#,
        r#""allele_string": "-/AA/AC", "input": "21\t9489854\t.\tT\tTAA,TAC\t.\t.\t."}"#
    ));
    assert_eq!(uploaded_variation(&insertion), "21_9489854_T/TAA/TAC");
    assert_eq!(record_location(&insertion), "21:9489854-9489855");
    assert_eq!(format_coords(&Value::Int(5), &Value::Null), "5-?");
    assert_eq!(format_coords(&Value::Null, &Value::Int(7)), "?-7");
}

#[test]
fn feature_kind_vcf_escape_and_most_severe() {
    assert_eq!(feature_kind("Transcript"), "transcript");
    assert_eq!(feature_kind("RegulatoryFeature"), "regulatory_feature");
    assert_eq!(feature_kind("MotifFeature"), "motif_feature");
    assert_eq!(feature_kind("-"), "intergenic");
    assert_eq!(feature_kind("regulatory_feature"), "regulatory_feature");
    assert_eq!(vcf_escape("a,b;c d|e"), "a&b%3Bc_d&e");
    assert_eq!(reverse_complement("ACGT"), "ACGT");
    assert_eq!(reverse_complement("AAC"), "GTT");
    assert_eq!(reverse_complement("deletion"), "deletion");
    let terms = |items: &[&str]| -> BTreeSet<String> {
        items.iter().map(|term| term.to_string()).collect()
    };
    assert_eq!(
        most_severe(&terms(&["intron_variant", "splice_region_variant", "?"])),
        "splice_region_variant"
    );
    assert_eq!(most_severe(&terms(&["?"])), "?");
    assert!(values_equal("100.00", "100"));
    assert!(!values_equal("1/1", "1"));
}

/// JSON.pm:327 keeps `-` for Allele, USED_REF and GIVEN_REF where every other `-` is
/// dropped; JSON.pm:229-233 writes an object carrying only `input` for a skipped record.
#[test]
fn self_parity_keeps_the_dash_of_used_ref_and_given_ref_and_skips_an_input_only_object() {
    let [default, tab, vcf, json] = committed();
    let mut rows = data_lines(&default);
    rows[0] = rows[0].replace("IMPACT=MODIFIER", "IMPACT=MODIFIER;USED_REF=-;GIVEN_REF=-");
    let default = rebuild(&default, &rows);
    let mut tab_lines = header_lines(&tab);
    let tab_header = tab_lines.pop().expect("a header line") + "\tUSED_REF\tGIVEN_REF";
    tab_lines.push(tab_header);
    tab_lines.extend(data_lines(&tab).into_iter().map(|row| row + "\t-\t-"));
    let tab = joined(&tab_lines);
    let mut json_lines = text_lines(&json);
    let mut object = members(&json_lines[0]);
    let first_consequence =
        |object: &mut IndexMap<String, Value>, change: &dyn Fn(&mut IndexMap<String, Value>)| {
            match &mut object["intergenic_consequences"] {
                Value::List(consequences) => match &mut consequences[0] {
                    Value::Dict(consequence) => change(consequence),
                    other => panic!("no object: {other:?}"),
                },
                other => panic!("no list: {other:?}"),
            }
        };
    first_consequence(&mut object, &|consequence| {
        consequence.insert("used_ref".to_string(), string("-"));
        consequence.insert("given_ref".to_string(), string("-"));
    });
    json_lines[0] = dumped_object(&object);
    let skipped_line = "21\t99\tskipped_by_size\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=99999999";
    let skipped: IndexMap<String, Value> =
        IndexMap::from([("input".to_string(), string(skipped_line))]);
    json_lines.push(dumped_object(&skipped));
    let mut vcf_rows = data_lines(&vcf);
    vcf_rows.push(skipped_line.to_string());
    let vcf = rebuild(&vcf, &vcf_rows);
    let report = parity(&default, &tab, &vcf, &joined(&json_lines));
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    let compared = texts(at(&report, &["pairs", "json", "columns_compared"]));
    assert!(compared.contains(&"USED_REF") && compared.contains(&"GIVEN_REF"));
    assert_eq!(number(&report, &["records"]), RECORDS as i64 + 1);
    first_consequence(&mut object, &|consequence| {
        consequence.shift_remove("given_ref");
    });
    json_lines[0] = dumped_object(&object);
    let report = parity(&default, &tab, &vcf, &joined(&json_lines));
    assert_eq!(
        number(
            &report,
            &["pairs", "json", "differing", "GIVEN_REF", "count"]
        ),
        1
    );
}

/// A structural variant whose QUAL is `0` is echoed with `.` in the JSON `input` and
/// written back with `.` in the VCF output (Parser/VCF.pm:464), so with the input records
/// given the object and the VCF line of record 11 name the same record only through their
/// first five columns; the committed formats, read with such an input, stay in parity.
#[test]
fn self_parity_keys_json_objects_by_the_first_five_columns_of_the_input() {
    let [default, tab, vcf, json] = committed();
    // The lines whole: the test reads and rewrites the sixth column of one.
    let inputs = data_lines(&read_text(&corpus().join("variants.vcf")).expect("the input reads"));
    let mut qual0 = inputs.clone();
    let columns: Vec<&str> = qual0[11].split('\t').collect();
    assert_eq!(columns[2..6], ["multi_ins_ins", "N", "<INS>,<INS>", "."]);
    qual0[11] = qual0[11].replacen("\t<INS>,<INS>\t.\t", "\t<INS>,<INS>\t0\t", 1);
    let echo = text_lines(&json)
        .iter()
        .find(|line| line.contains("\"multi_ins_ins\""))
        .map(|line| members(line)["input"].clone())
        .expect("the record's object");
    assert_eq!(echo.as_str(), Some(inputs[11].as_str()));
    assert_ne!(inputs[11], qual0[11]);
    let report = parity_of(&default, &tab, &vcf, &json, Some(&qual0));
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["records"]), RECORDS as i64);
}

/// Rows of the individual-zyg reference whose ZYG list names two or more samples, so reversing
/// the list changes the text: 1,530 of the 1,534 rows carrying ZYG in the default format, and
/// 103 of the 104 record objects carrying zyg in the JSON format.
const REORDERABLE_ZYG_ROWS: i64 = 1530;
const REORDERABLE_ZYG_OBJECTS: i64 = 103;

/// `text` with the items of every `column=` list of the Extra text reversed.
fn reversed_lists(text: &str, column: &str) -> String {
    let lists = regex::Regex::new(&format!(r"{column}=([^;\t\n]+)")).expect("the pattern compiles");
    lists
        .replace_all(text, |found: &regex::Captures| {
            let items: Vec<&str> = found[1].split(',').rev().collect();
            format!("{column}={}", items.join(","))
        })
        .into_owned()
}

/// The report of two entry lists compared with the lists of `fmt` read as multisets.
fn compare_as_lists(
    expected: &[Entry],
    actual: &[Entry],
    documented: &mut Documented,
    fmt: &str,
) -> Value {
    let side = Side::default();
    let lists = order_only_fields(fmt);
    let options = CompareOptions {
        order_only: Some(&lists),
        ..CompareOptions::new(&side)
    };
    compare_entries(expected, actual, documented, &options).expect("the entries compare")
}

/// `--individual_zyg` prints one `sample:ZYG` entry per named sample in hash order
/// (OutputFactory.pm:983-987), so a candidate naming the samples in another order agrees as
/// a multiset under the class in every format: the default and tab `ZYG` column with the
/// `Extra` text carrying it, the `&`-joined CSQ field, and the record object's `zyg`
/// array, which the JSON format carries once at the top level, so `ZYG` is a record-level
/// column that self-parity sets aside against the consequence entries and that has no JSON
/// floor; a changed zygosity or a sample fewer stays red.
#[test]
fn the_order_only_class_covers_the_zyg_list_in_every_format() {
    let individual_zyg = golden().join("GRCh38-individual-zyg");
    let expected_dir = individual_zyg.join("expected_individual-zyg");
    let text = |name: &str| read_text(&expected_dir.join(name)).expect("the output reads");
    let inputs = input_records(&individual_zyg.join("variants.vcf")).expect("the input reads");
    let default = text("default.txt.gz");
    let default_entries = |text: &str| parse_default(text).expect("the text parses");
    let expected = default_entries(&default);
    let actual = default_entries(&reversed_lists(&default, "ZYG"));
    let plain = compare(&expected, &actual, &mut no_documented());
    assert_eq!(verdict(&plain), "red");
    for column in ["ZYG", "Extra"] {
        assert_eq!(
            number(&plain, &["columns", column, "differing"]),
            REORDERABLE_ZYG_ROWS
        );
    }
    let mut documented = no_documented();
    let accepted = compare_as_lists(&expected, &actual, &mut documented, "default");
    assert_eq!(verdict(&accepted), "green");
    assert_eq!(number(&accepted, &["columns", "ZYG", "differing"]), 0);
    for column in ["ZYG", "Extra"] {
        assert_eq!(
            number(&accepted, &["columns", column, "order_only"]),
            REORDERABLE_ZYG_ROWS
        );
    }
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(ORDER_ONLY, 2 * REORDERABLE_ZYG_ROWS)])
    );
    let example = texts(at(&accepted, &["columns", "ZYG", "order_only_examples"]))[0];
    assert!(example.starts_with(
        "21:9614117-9614200|-|ENST00000842019|Transcript: expected 'S7:HOMREF,S4:HOMREF,"
    ));
    assert!(example.contains(
        "got 'S1:HOMREF,S5:HOMREF,S6:HOMREF,S8:HOMREF,S3:HOMREF,S2:HOMREF,S4:HOMREF,S7:HOMREF'"
    ));

    // The tab format's ZYG column, the same rule.
    let tab = text("tab.txt.gz");
    let zyg_column = header_columns(&tab)
        .iter()
        .position(|column| column == "ZYG")
        .expect("a ZYG column");
    let tab_rows: Vec<String> = data_lines(&tab)
        .iter()
        .map(|row| {
            let mut cells: Vec<String> = row.split('\t').map(str::to_string).collect();
            let items: Vec<&str> = cells[zyg_column].split(',').rev().collect();
            cells[zyg_column] = items.join(",");
            cells.join("\t")
        })
        .collect();
    let tab_entries = |text: &str| parse_tab(text).expect("the text parses");
    let tab_report = compare_as_lists(
        &tab_entries(&tab),
        &tab_entries(&rebuild(&tab, &tab_rows)),
        &mut no_documented(),
        "tab",
    );
    assert_eq!(verdict(&tab_report), "green");
    assert_eq!(
        number(&tab_report, &["columns", "ZYG", "order_only"]),
        REORDERABLE_ZYG_ROWS
    );

    // The VCF format joins the items with & (VCF.pm:391); the sample columns follow INFO.
    let vcf = text("vcf.vcf.gz");
    let csq_header = text_lines(&vcf)
        .into_iter()
        .find(|line| line.starts_with("##INFO=<ID=CSQ"))
        .expect("a CSQ header");
    let zyg_field = csq_format_fields(&csq_header)
        .expect("a Format list")
        .iter()
        .position(|field| field == "ZYG")
        .expect("a ZYG field");
    let flipped: Vec<String> = data_lines(&vcf)
        .iter()
        .map(|row| {
            let mut cols: Vec<String> = row.split('\t').map(str::to_string).collect();
            let info: Vec<String> = cols[7]
                .split(';')
                .map(|kv| match kv.strip_prefix("CSQ=") {
                    None => kv.to_string(),
                    Some(csq) => {
                        let chunks: Vec<String> = csq
                            .split(',')
                            .map(|chunk| {
                                let mut vals: Vec<String> =
                                    chunk.split('|').map(str::to_string).collect();
                                let items: Vec<&str> = vals[zyg_field].split('&').rev().collect();
                                vals[zyg_field] = items.join("&");
                                vals.join("|")
                            })
                            .collect();
                        format!("CSQ={}", chunks.join(","))
                    }
                })
                .collect();
            cols[7] = info.join(";");
            cols.join("\t")
        })
        .collect();
    let vcf_expected = parse_vcf(&vcf, &inputs).expect("the text parses");
    let vcf_actual = parse_vcf(&rebuild(&vcf, &flipped), &inputs).expect("the text parses");
    let vcf_plain = compare(&vcf_expected, &vcf_actual, &mut no_documented());
    assert_eq!(
        number(&vcf_plain, &["columns", "ZYG", "differing"]),
        REORDERABLE_ZYG_ROWS
    );
    let vcf_report = compare_as_lists(&vcf_expected, &vcf_actual, &mut no_documented(), "vcf");
    assert_eq!(verdict(&vcf_report), "green");
    assert_eq!(
        number(&vcf_report, &["columns", "ZYG", "order_only"]),
        REORDERABLE_ZYG_ROWS
    );

    // The JSON format's zyg array, once per record object, its strings compared as a multiset.
    let json_text = text("json.jsonl.gz");
    let objects: Vec<String> = text_lines(&json_text)
        .iter()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let mut object = members(line);
            if let Some(Value::List(zyg)) = object.get_mut("zyg") {
                zyg.reverse();
            }
            dumped_object(&object)
        })
        .collect();
    let reordered_text = joined(&objects);
    let json_expected = parse_json(&json_text, &inputs).expect("the text parses");
    let json_actual = parse_json(&reordered_text, &inputs).expect("the text parses");
    let json_plain = compare(&json_expected, &json_actual, &mut no_documented());
    assert_eq!(
        number(&json_plain, &["columns", "zyg", "differing"]),
        REORDERABLE_ZYG_OBJECTS
    );
    let json_report = compare_as_lists(&json_expected, &json_actual, &mut no_documented(), "json");
    assert_eq!(verdict(&json_report), "green");
    assert_eq!(
        number(&json_report, &["columns", "zyg", "order_only"]),
        REORDERABLE_ZYG_OBJECTS
    );

    // A changed zygosity, or a sample fewer, is not an order-only agreement.
    let full_list =
        "ZYG=S7:HOMREF,S4:HOMREF,S2:HOMREF,S3:HOMREF,S8:HOMREF,S6:HOMREF,S5:HOMREF,S1:HOMREF";
    let rows = data_lines(&default);
    let target = rows
        .iter()
        .position(|row| row.contains(full_list))
        .expect("a row of the full list");
    let with_list = |list: String| {
        let mut rows = rows.clone();
        rows[target] = rows[target].replacen(full_list, &list, 1);
        default_entries(&rebuild(&default, &rows))
    };
    let changed = compare_as_lists(
        &expected,
        &with_list(full_list.replace("S1:HOMREF", "S1:HET")),
        &mut no_documented(),
        "default",
    );
    assert_eq!(verdict(&changed), "red");
    assert_eq!(number(&changed, &["columns", "ZYG", "differing"]), 1);
    assert_eq!(number(&changed, &["columns", "ZYG", "order_only"]), 0);
    let fewer = compare_as_lists(
        &expected,
        &with_list(full_list.replace(",S1:HOMREF", "")),
        &mut no_documented(),
        "default",
    );
    assert_eq!(verdict(&fewer), "red");
    assert_eq!(number(&fewer, &["columns", "ZYG", "differing"]), 1);

    // ZYG is a record-level column: no JSON floor, and the self-parity of one engine's formats
    // checks it once per record, the rows' list against the object's zyg array as multisets,
    // instead of against the JSON consequence entries, which carry no zyg key.
    assert!(RECORD_LEVEL_COLUMNS.contains(&"ZYG"));
    assert_eq!(json_floor_fields("ZYG"), None);
    let report = parity_of(&default, &tab, &vcf, &json_text, Some(&inputs));
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert!(!texts(at(&report, &["not_compared", "default_record_level"])).contains(&"ZYG"));
    for name in ["vcf", "tab"] {
        assert!(texts(at(&report, &["pairs", name, "columns_compared"])).contains(&"ZYG"));
    }
    let reordered = parity_of(&default, &tab, &vcf, &reordered_text, Some(&inputs));
    assert_eq!(verdict(&reordered), "green", "{:?}", reasons(&reordered));
    let changed = reordered_text.replacen("\"S1:HOMREF\"", "\"S1:HET\"", 1);
    let report = parity_of(&default, &tab, &vcf, &changed, Some(&inputs));
    assert_eq!(verdict(&report), "red");
    let problems = texts(at(&report, &["record_checks", "problems"]));
    assert_eq!(problems.len(), 1);
    assert!(problems[0].ends_with(
        "against the objects' zyg ['S1:HET,S2:HOMREF,S3:HOMREF,S4:HOMREF,S5:HOMREF,S6:HOMREF,S7:HOMREF,S8:HOMREF']"
    ));
}

/// Under `--individual` a reference-homozygous individual's object carries its `zyg`
/// and no consequence, and the default format writes no row for it (OutputFactory.pm:511),
/// so the record check holds only the objects that carry a consequence to the rows.
#[test]
fn self_parity_holds_zyg_to_the_rows_only_for_objects_with_a_consequence() {
    let individual = golden().join("GRCh38-individual");
    let report = parity_of_run(
        &individual.join("expected_individual-all"),
        &individual.join("variants.vcf"),
    );
    assert!(texts(at(&report, &["record_checks", "problems"])).is_empty());
}

/// FREQS is a record-level column the JSON format never carries (JSON.pm:493-495), so
/// the self-parity comparison of one engine's formats lists it under the columns not
/// compared against JSON instead of reporting the missing key as a difference.
#[test]
fn self_parity_sets_freqs_aside_as_a_record_level_column() {
    let freqfilter = golden().join("GRCh38-freqfilter");
    let report = parity_of_run(
        &freqfilter.join("expected_freqfilter"),
        &freqfilter.join("variants.vcf"),
    );
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert!(texts(at(&report, &["not_compared", "default_record_level"])).contains(&"FREQS"));
    assert!(texts(at(&report, &["pairs", "vcf", "columns_compared"])).contains(&"FREQS"));
}
