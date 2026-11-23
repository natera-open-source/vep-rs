// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of the tab, VCF and JSON formats: how a row is found, how a value is spelt and
//! what is refused.

use super::*;
use crate::TempDir;
use std::fs;

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-formats-").unwrap()
}

fn file(directory: &TempDir, name: &str, text: &str) -> String {
    let path = directory.path.join(name);
    fs::write(&path, text).unwrap();
    path.to_str().unwrap().to_string()
}

fn records(lines: &[&str]) -> Vec<RunRecord> {
    let record = |line: &&str| RunRecord {
        vcf: line.to_string(),
        manifest_index: None,
        reference_warning: None,
    };
    lines.iter().map(record).collect()
}

fn key(location: &str, allele: &str, feature: &str, consequences: &str) -> FieldKey {
    FieldKey {
        location: location.to_string(),
        allele: allele.to_string(),
        feature: feature.to_string(),
        feature_type: "Transcript".to_string(),
        consequences: consequences.to_string(),
    }
}

fn refusal<T: std::fmt::Debug>(read: Result<T, Failure>, path: &str, what: &str) {
    assert_eq!(read.unwrap_err(), Failure::own(format!("{path} {what}")));
}

const INPUT: [&str; 3] = [
    "21\t100\trs1\tA\tG\t.\tPASS\t.",
    "21\t200\t.\tA\tT\t.\tPASS\t.",
    "21\t200\t.\tA\tT\t.\tPASS\t.",
];

const CSQ_HEADER: &str = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature|PLUGIN_B|PLUGIN_A\">\n";

#[test]
fn a_float_is_written_as_the_golden_harness_writes_it() {
    let cases: [(f64, &str); 17] = [
        (0.0, "0.0"),
        (-0.0, "-0.0"),
        (1.0, "1.0"),
        (-6.0, "-6.0"),
        (0.15, "0.15"),
        (0.1, "0.1"),
        (12.34, "12.34"),
        (0.001234, "0.001234"),
        (0.00001, "0.00001"),
        (0.000001, "1e-6"),
        (1.5e-7, "1.5e-7"),
        (1e15, "1000000000000000.0"),
        (1e16, "1e+16"),
        (1.2345678901234568e20, "1.2345678901234568e+20"),
        (1.5e100, "1.5e+100"),
        (1e-100, "1e-100"),
        (-2.5e-7, "-2.5e-7"),
    ];
    for (value, text) in cases {
        assert_eq!(float_text(value), text, "{value:e}");
    }
}

#[test]
fn a_member_is_a_text_a_joined_list_or_compact_json_with_its_names_in_byte_order() {
    let member = |text: &str| member_text(&loads(text).unwrap());
    assert_eq!(member(r#""a b""#), "a b");
    assert_eq!(member(r#"["x", "y", 3]"#), "x,y,3");
    assert_eq!(member("-6"), "-6");
    assert_eq!(member("0.10"), "0.1");
    assert_eq!(member("18446744073709551615"), "18446744073709551615");
    assert_eq!(member("18446744073709551616"), "1.8446744073709552e+19");
    assert_eq!(member("true"), "true");
    assert_eq!(member("null"), "null");
    assert_eq!(
        member(r#"{"b": 0.50, "SYMBOL": "G\"1", "a": [1, {"z": null, "y": "é\n"}], "A": -6}"#),
        "{\"A\":-6,\"SYMBOL\":\"G\\\"1\",\"a\":[1,{\"y\":\"é\\n\",\"z\":null}],\"b\":0.5}"
    );
    // A list of objects joins their compact texts, as the harness prints `domains`.
    assert_eq!(
        member(r#"[{"name": "n", "db": "d"}, {"db": "e"}]"#),
        r#"{"db":"d","name":"n"},{"db":"e"}"#
    );
}

#[test]
fn a_tab_row_is_found_by_its_location_and_holds_every_column() {
    let scratch = scratch();
    let text = "## a comment\n\
        #Uploaded_variation\tLocation\tAllele\tFeature\tFeature_type\tConsequence\tPLUGIN\n\
        rs1\t21:100\tG\tENST1\tTranscript\tb_variant,a_variant\tfirst\n\
        rs1\t21:100\tG\tENST1\tTranscript\ta_variant,b_variant\tsecond\n\
        rs1\t21:100\tG\tENST2\tTranscript\ta_variant\t-\n";
    let rows = tab_rows(&file(&scratch, "tab.txt", text)).unwrap();
    let first = rows
        .row(&key("21:100", "G", "ENST1", "a_variant,b_variant"), &[])
        .unwrap();
    assert_eq!(first["PLUGIN"], "first", "the first row of a key is read");
    assert_eq!(first["Consequence"], "b_variant,a_variant");
    assert_eq!(first.len(), 7);
    let other = key("21:100", "G", "ENST2", "a_variant");
    assert_eq!(rows.row(&other, &[]).unwrap()["PLUGIN"], "-");
    assert_eq!(
        rows.row(&key("21:100", "G", "ENST3", "a_variant"), &[0]),
        None
    );

    let headless = file(&scratch, "headless.txt", "rs1\t21:100\tG\n");
    refusal(
        tab_rows(&headless),
        &headless,
        "holds a row before the line that names its columns",
    );
    let short = file(&scratch, "short.txt", "#Location\tAllele\n21:100\n");
    refusal(
        tab_rows(&short),
        &short,
        "holds a row of 1 columns under a line that names 2",
    );
    let nowhere = scratch.path.join("nowhere.txt");
    let read = tab_rows(nowhere.to_str().unwrap()).unwrap_err();
    assert!(
        read.line.contains("is not a readable text file ("),
        "{}",
        read.line
    );
}

#[test]
fn a_tab_row_is_found_without_a_feature_type_column_and_holds_its_vrs_object_by_member_name() {
    let scratch = scratch();
    // The columns of a `--fields` list that leaves out Feature_type, and an object whose
    // members stand in another order than that of their names.
    let text = "#Location\tAllele\tFeature\tConsequence\tGA4GH_VRS\n\
        21:100\tG\tENST1\ta_variant\t{\"type\": \"Allele\", \"state\": {\"b\": 1, \"a\": \"x\"}}\n\
        21:200\tT\tENST2\ta_variant\t-\n";
    let rows = tab_rows(&file(&scratch, "fields.txt", text)).unwrap();
    let row = rows
        .row(&key("21:100", "G", "ENST1", "a_variant"), &[])
        .unwrap();
    assert_eq!(
        row["GA4GH_VRS"],
        r#"{"state":{"a":"x","b":1},"type":"Allele"}"#
    );
    let without = rows
        .row(&key("21:200", "T", "ENST2", "a_variant"), &[])
        .unwrap();
    assert_eq!(without["GA4GH_VRS"], "-", "a value that is no JSON stays");
    assert_eq!(
        rows.row(&key("21:100", "G", "ENST2", "a_variant"), &[]),
        None
    );
}

#[test]
fn a_vcf_entry_is_found_by_its_record_and_spelt_as_the_default_format_spells_it() {
    let scratch = scratch();
    let records = records(&INPUT);
    let pairing = Pairing::of(&records, None);
    // The third line is a line of no record, the fourth and fifth are the two identical
    // records in their order, and an earlier CSQ of a line is not the run's.
    let text = format!(
        "##fileformat=VCFv4.2\n{CSQ_HEADER}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t100\trs1\tA\tG\t.\tPASS\tAC=1;CSQ=G|old|ENST1|x|y;CSQ=G|b_variant&a_variant|ENST1|26&G|,G|a_variant||p.Ala1%3D|0.10\n\
         21\t150\trs9\tA\tG\t.\tPASS\tCSQ=G|a_variant|ENST1|lost|lost\n\
         21\t200\t.\tA\tT\t.\tPASS\tCSQ=T|a_variant|ENST2|second|\n\
         21\t200\t.\tA\tT\t.\tPASS\tCSQ=T|a_variant|ENST2|third|\n\
         21\t100\trs1\tA\tG\t.\tPASS\tCSQ=G|a_variant|ENST1|again|\n\
         21\t100\trs1\tA\tG\t.\tPASS\tAC=1\n"
    );
    let rows = vcf_rows(&file(&scratch, "out.vcf", &text), &pairing, "CSQ").unwrap();
    let transcript = key("21:100", "G", "ENST1", "a_variant,b_variant");
    let row = rows.row(&transcript, &[0]).unwrap();
    assert_eq!(
        row["PLUGIN_B"], "26,G",
        "an `&` is the default format's comma"
    );
    assert_eq!(
        row["PLUGIN_A"], "-",
        "an empty value is the default format's `-`"
    );
    assert_eq!(row["Consequence"], "b_variant,a_variant");
    assert_eq!(
        row.keys().collect::<Vec<_>>(),
        ["Allele", "Consequence", "Feature", "PLUGIN_A", "PLUGIN_B"],
        "the fields stand in the byte order of their names"
    );
    let intergenic = rows
        .row(&key("21:100", "G", "-", "a_variant"), &[0])
        .unwrap();
    assert_eq!(
        (
            intergenic["PLUGIN_B"].as_str(),
            intergenic["PLUGIN_A"].as_str()
        ),
        ("p.Ala1%3D", "0.10")
    );
    // A row is its record's: the key of another record, and of no record, finds none.
    assert_eq!(rows.row(&transcript, &[1]), None);
    assert_eq!(rows.row(&transcript, &[]), None);
    let twin = key("21:200", "T", "ENST2", "a_variant");
    assert_eq!(rows.row(&twin, &[1]).unwrap()["PLUGIN_B"], "second");
    assert_eq!(rows.row(&twin, &[2]).unwrap()["PLUGIN_B"], "third");
    assert_eq!(rows.row(&twin, &[2, 1]).unwrap()["PLUGIN_B"], "third");
    assert_eq!(rows.by_record.len(), 4, "a line of no record is not read");

    // Another INFO key is the run's own.
    let other = text.replace("ID=CSQ,", "ID=ANN,").replace("CSQ=", "ANN=");
    let other = file(&scratch, "other.vcf", &other);
    assert_eq!(
        vcf_rows(&other, &pairing, "CSQ").unwrap(),
        FormatRows::default()
    );
    assert_eq!(
        vcf_rows(&other, &pairing, "ANN").unwrap().by_record.len(),
        4
    );
    let flags = |words: &[&str]| {
        words
            .iter()
            .map(|word| word.to_string())
            .collect::<Vec<_>>()
    };
    assert_eq!(info_field(&flags(&["--hgvs"])), "CSQ");
    assert_eq!(
        info_field(&flags(&["--vcf_info_field", "ANN", "--hgvs"])),
        "ANN"
    );
    assert_eq!(info_field(&flags(&["--vcf_info_field=ANN"])), "ANN");

    let narrow = file(&scratch, "narrow.vcf", "21\t100\trs1\tA\tG\t.\tPASS\n");
    refusal(
        vcf_rows(&narrow, &pairing, "CSQ"),
        &narrow,
        "holds a line of fewer than eight columns",
    );
    let unnamed = file(
        &scratch,
        "unnamed.vcf",
        "21\t100\trs1\tA\tG\t.\tPASS\tCSQ=G|a\n",
    );
    refusal(
        vcf_rows(&unnamed, &pairing, "CSQ"),
        &unnamed,
        "holds an entry of 2 fields under a header that names 0 for CSQ",
    );
    let wide = format!("{CSQ_HEADER}21\t100\trs1\tA\tG\t.\tPASS\tCSQ=G|a|b|c|d|e\n");
    let wide = file(&scratch, "wide.vcf", &wide);
    refusal(
        vcf_rows(&wide, &pairing, "CSQ"),
        &wide,
        "holds an entry of 6 fields under a header that names 5 for CSQ",
    );
}

#[test]
fn a_json_consequence_object_is_found_by_its_record_and_holds_its_members() {
    let scratch = scratch();
    let records = records(&INPUT);
    let pairing = Pairing::of(&records, None);
    let text = concat!(
        r#"{"input": "21\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1", "most_severe_consequence": "b_variant", "transcript_consequences": [{"variant_allele": "G", "transcript_id": "ENST1", "consequence_terms": ["b_variant", "a_variant"], "plugin": {"B": "26,G", "A": 0.10}, "strand": 1}, {"variant_allele": "G", "transcript_id": "ENST1", "consequence_terms": ["a_variant", "b_variant"], "strand": -1}], "regulatory_feature_consequences": [{"variant_allele": "G", "regulatory_feature_id": "ENSR1", "consequence_terms": ["r_variant"]}], "motif_feature_consequences": [{"variant_allele": "G", "motif_feature_id": "ENSM1", "consequence_terms": ["m_variant"]}], "intergenic_consequences": [{"variant_allele": "G", "consequence_terms": ["i_variant"]}]}"#,
        "\n\n",
        r#"{"input": "21\t150\trs9\tA\tG", "transcript_consequences": [{"variant_allele": "G", "transcript_id": "ENST1", "consequence_terms": ["a_variant"]}]}"#,
        "\n",
        r#"{"input": "21\t200\t.\tA\tT\t.\tPASS\t.", "transcript_consequences": [{"transcript_id": "ENST2", "consequence_terms": ["a_variant"]}]}"#,
        "\n",
    );
    let rows = json_rows(&file(&scratch, "out.json", text), &pairing).unwrap();
    let row = rows
        .row(&key("21:100", "G", "ENST1", "a_variant,b_variant"), &[0])
        .unwrap();
    assert_eq!(row["plugin"], r#"{"A":0.1,"B":"26,G"}"#);
    assert_eq!(row["strand"], "1", "the first object of a key is read");
    assert_eq!(row["consequence_terms"], "b_variant,a_variant");
    for (feature, terms) in [
        ("ENSR1", "r_variant"),
        ("ENSM1", "m_variant"),
        ("-", "i_variant"),
    ] {
        assert!(
            rows.row(&key("21:100", "G", feature, terms), &[0])
                .is_some(),
            "{feature}"
        );
    }
    // An object without the allele is a row of the allele `-`, and an object of no record
    // is not read.
    assert!(rows
        .row(&key("21:200", "-", "ENST2", "a_variant"), &[1])
        .is_some());
    assert_eq!(rows.by_record.len(), 5);

    let text_line = file(&scratch, "list.json", "[1]\n");
    refusal(
        json_rows(&text_line, &pairing),
        &text_line,
        "line 1 is no JSON object",
    );
    let cut = file(&scratch, "cut.json", "\n{\"input\": \n");
    refusal(json_rows(&cut, &pairing), &cut, "line 2 is no JSON object");
    let item = r#"{"input": "21\t100\trs1\tA\tG", "motif_feature_consequences": ["x"]}"#;
    let item = file(&scratch, "item.json", item);
    refusal(
        json_rows(&item, &pairing),
        &item,
        "line 1 holds under motif_feature_consequences what is no object",
    );
}

#[test]
fn a_multiset_column_agrees_under_its_name_in_any_letter_case() {
    assert!(agree("domains", "a,b", "b,a", true));
    assert!(agree("DOMAINS", "a,b", "b,a", true));
    assert!(!agree("domains", "a,b", "b,a", false));
    assert!(!agree("plugin", "a,b", "b,a", true));
    assert!(agree("plugin", "a,b", "a,b", false));
}

#[test]
fn one_output_of_a_format_without_the_other_is_refused() {
    let both: [Named<'_>; 3] = [
        ("tab", None, None),
        ("vcf", Some("a"), Some("b")),
        ("json", None, None),
    ];
    assert_eq!(unpaired(&both), None);
    let alone: [Named<'_>; 2] = [("tab", Some("a"), Some("b")), ("json", None, Some("b"))];
    assert_eq!(
        unpaired(&alone),
        Some(Failure::own(
            "--vep-json and --vep-rs-json name the two outputs of one format: give both or neither"
        ))
    );
    // A format named by neither output is not read.
    assert!(load(&both[..1], &Pairing::of(&[], None), &[])
        .unwrap()
        .is_empty());
}

#[test]
fn the_outputs_of_a_reencoded_input_are_paired_by_its_lines_and_in_turn() {
    // A run over `variants.vcf` re-encoded in another format: a JSON object names its record
    // by the input's own line, and the lines of a VCF output, which carry no input line,
    // are the records in turn.
    let scratch = scratch();
    let records = records(&INPUT);
    let lines: Vec<String> = [
        "21\t100\t100\tA/G\t+",
        "21\t200\t200\tA/T\t+\tfirst",
        "21\t200\t200\tA/T\t+\tsecond",
    ]
    .map(str::to_string)
    .to_vec();
    let pairing = Pairing::of(&records, Some(&lines));
    let object = |input: &str, feature: &str| {
        let input = input.replace('\t', "\\t");
        format!(
            r#"{{"input": "{input}", "transcript_consequences": [{{"variant_allele": "G", "transcript_id": "{feature}", "consequence_terms": ["a_variant"]}}]}}"#
        ) + "\n"
    };
    // The objects stand in another order than the records, and the record of the manifest's
    // own line is no line of this input.
    let text =
        object(&lines[2], "ENST3") + &object(&lines[0], "ENST1") + &object(INPUT[0], "ENST9");
    let rows = json_rows(&file(&scratch, "out.json", &text), &pairing).unwrap();
    assert!(rows
        .row(&key("21:100", "G", "ENST1", "a_variant"), &[0])
        .is_some());
    assert!(rows
        .row(&key("21:200", "G", "ENST3", "a_variant"), &[1])
        .is_some());
    assert_eq!(rows.by_record.len(), 2);

    let line = |feature: &str| format!("21\t1\t.\tN\tG\t.\t.\tCSQ=G|a_variant|{feature}||\n");
    let text = format!(
        "{CSQ_HEADER}{}{}{}",
        line("ENST1"),
        line("ENST2"),
        line("ENST3")
    );
    let rows = vcf_rows(&file(&scratch, "out.vcf", &text), &pairing, "CSQ").unwrap();
    for (ordinal, feature) in ["ENST1", "ENST2", "ENST3"].into_iter().enumerate() {
        let found = rows.row(&key("21:1", "G", feature, "a_variant"), &[ordinal]);
        assert!(found.is_some(), "{feature}");
    }
    assert_eq!(
        rows.row(&key("21:1", "G", "ENST1", "a_variant"), &[1]),
        None
    );
    let short = format!("{CSQ_HEADER}{}{}", line("ENST1"), line("ENST2"));
    let short = file(&scratch, "short.vcf", &short);
    refusal(
        vcf_rows(&short, &pairing, "CSQ"),
        &short,
        "holds 2 lines for the 3 records of an input that is no VCF, whose lines are the records in turn",
    );
}

#[test]
fn the_two_further_multiset_columns_agree_under_the_order_only_mask_alone() {
    for name in ["NearestExonJB", "nearestexonjb"] {
        assert!(agree(name, "a,b,b", "b,a,b", true), "{name}");
        assert!(!agree(name, "a,b,b", "b,a,a", true), "{name}");
        assert!(!agree(name, "a,b", "b,a", false), "{name}");
        assert!(!agree(name, "a,b", "-", true), "{name}");
    }
    // The groups of the pair-group column in any order, and the pairs of a group in any
    // order.
    assert!(agree("5UTR_annotation", "a=1:b=2,c=3", "c=3,b=2:a=1", true));
    assert!(!agree(
        "5UTR_annotation",
        "a=1:b=2,c=3",
        "a=1,b=2:c=3",
        true
    ));
    assert!(!agree("5UTR_annotation", "a=1:b=2", "b=2:a=1", false));
    assert!(!agree("5utr_annotation", "a=1:b=2", "b=2:a=1", true));
}

#[test]
fn a_value_with_a_tab_or_a_line_end_is_no_cell_of_a_line() {
    assert!(one_cell("a b|c,d"));
    for value in ["a\tb", "a\nb", "a\rb"] {
        assert!(!one_cell(value), "{value:?}");
    }
}
