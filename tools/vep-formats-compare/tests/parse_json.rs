// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The JSON parser and its value forms, on the committed Ensembl VEP outputs of
//! `tests/golden/116/GRCh38-release116` and on planted lines with their expected values.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::{Error, Uncaught};
use vep_formats_compare::input::{input_records, read_text};
use vep_formats_compare::parse::json::*;
use vep_formats_compare::parse::vcf::{parse_vcf, vcf_line_key, InputIndex};
use vep_formats_compare::pyjson::{loads, Value};

const RECORDS: usize = 36;
const ROWS: usize = 977;

const FIRST_OBJECT: &str = r#"{"input": "1\t100\trs1\tG\tGTAT,C\t.\t.\t.\tmore", "start": 101, "end": 100, "allele_string": "-/TAT/C", "seq_region_name": "1", "id": "rs1", "vep_rs_version": "0.1", "most_severe_consequence": null, "colocated_variants": [{"id": "x", "f": 0.25}, {"id": "y"}], "transcript_consequences": [{"variant_allele": "TAT", "transcript_id": "T1", "consequence_terms": ["a", "b"], "n": 1.5e-05, "strand": -1}, {"variant_allele": "TAT", "transcript_id": "T1"}], "intergenic_consequences": [{"variant_allele": "C", "consequence_terms": ["i"], "intergenic_id": "never read"}]}"#;
const SECOND_OBJECT: &str = r#"{"input": "1\t200\t.\tAT\tA\t.\t.\t.", "start": 201, "end": 201, "allele_string": "T/-", "seq_region_name": "1", "regulatory_feature_consequences": [{"variant_allele": "-", "regulatory_feature_id": "R1"}], "motif_feature_consequences": [{"variant_allele": "-", "motif_feature_id": "M1"}], "transcript_consequences": "not a list"}"#;
const THIRD_OBJECT: &str = r#"{"input": 7, "allele_string": "A/T", "start": true, "end": 5, "transcript_consequences": [{"variant_allele": "T"}]}"#;

fn corpus() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/116/GRCh38-release116")
}

fn strings(items: &[&str]) -> Vec<String> {
    items.iter().map(|item| item.to_string()).collect()
}

fn planted_inputs() -> Vec<String> {
    strings(&[
        "1\t100\trs1\tG\tGTAT,C\t.\t.\t.",
        "1\t200\t.\tAT\tA\t.\t.\t.",
    ])
}

fn value(text: &str) -> Value {
    loads(text).unwrap()
}

fn keys(entries: &[Entry]) -> Vec<&str> {
    entries.iter().map(|entry| entry.key.as_str()).collect()
}

fn fields(entry: &Entry) -> Vec<(&str, &str)> {
    entry
        .fields
        .iter()
        .map(|(name, value)| (name.as_str(), value.as_str()))
        .collect()
}

fn comparator_error(result: Result<Vec<Entry>, Error>) -> String {
    match result {
        Err(Error::Comparator(text)) => text,
        Ok(entries) => panic!("{} entries and no error", entries.len()),
    }
}

/// What the run ends with on an input no reading of the comparison expects.
fn uncaught<T>(run: impl FnOnce() -> T + std::panic::UnwindSafe) -> String {
    match std::panic::catch_unwind(run) {
        Ok(_) => panic!("the run did not end"),
        Err(payload) => match payload.downcast_ref::<Uncaught>() {
            Some(Uncaught(what)) => what.clone(),
            None => std::panic::resume_unwind(payload),
        },
    }
}

/// `leaf` inside `depth` lists.
fn lists(depth: usize, leaf: &str) -> String {
    format!("{}{leaf}{}", "[".repeat(depth), "]".repeat(depth))
}

/// An object whose member `a` is `leaf` inside `depth` lists.
fn nested(depth: usize, leaf: &str) -> String {
    format!("{{\"a\": {}}}", lists(depth, leaf))
}

#[test]
fn the_committed_reference_json_gives_a_record_entry_an_object_and_an_entry_a_consequence() {
    let text = read_text(&corpus().join("expected/json.jsonl.gz")).unwrap();
    let inputs = input_records(&corpus().join("variants.vcf")).unwrap();
    let entries = parse_json(&text, &inputs).unwrap();
    let (records, consequences): (Vec<&Entry>, Vec<&Entry>) = entries
        .iter()
        .partition(|entry| entry.key.ends_with("|record"));
    assert_eq!(records.len(), RECORDS);
    assert_eq!(consequences.len(), ROWS);
    let ins = records
        .iter()
        .find(|entry| entry.key == "rec11|record")
        .unwrap();
    assert_eq!(ins.fields["id"], "multi_ins_ins");
    assert_eq!(ins.fields["most_severe_consequence"], "feature_elongation");
    assert_eq!(ins.fields["start"], "30165781");
    assert_eq!(ins.fields["end"], "30165780");
    assert_eq!(ins.fields["allele_string"], "<INS>/<INS>");
    assert!(!ins.fields.contains_key("input"));
    assert!(!ins.fields.contains_key("transcript_consequences"));
    let items: Vec<&&Entry> = consequences
        .iter()
        .filter(|entry| entry.record == Some(11))
        .collect();
    assert_eq!(
        items
            .iter()
            .map(|entry| entry.key.as_str())
            .collect::<Vec<_>>(),
        [
            "rec11|insertion|ENST00000286808|transcript",
            "rec11|insertion|ENST00000286808|transcript#2",
        ]
    );
    assert_eq!(
        items[0].consequence,
        "feature_elongation,3_prime_UTR_variant"
    );
    assert_eq!(items[0].fields["cdna_start"], "1025");
    assert_eq!(items[0].fields["cdna_end"], "1026");
    assert_eq!(items[0].fields["strand"], "-1");
    let vcf_text = read_text(&corpus().join("expected/vcf.vcf.gz")).unwrap();
    let vcf_entries = parse_vcf(&vcf_text, &inputs).unwrap();
    let named =
        |entries: &[Entry], with_a_record: bool| -> BTreeSet<(Option<i64>, String, String)> {
            entries
                .iter()
                .filter(|entry| !with_a_record || entry.record.is_some())
                .map(|entry| (entry.record, entry.allele.clone(), entry.feature.clone()))
                .collect()
        };
    assert_eq!(named(&entries, false), named(&vcf_entries, true));
}

#[test]
fn a_json_value_is_one_comparable_text() {
    assert_eq!(json_scalar(&value("\"x\"")), "x");
    assert_eq!(json_scalar(&value("[\"a\", \"b\"]")), "a,b");
    assert_eq!(json_scalar(&value("1")), "1");
    assert_eq!(json_scalar(&value("true")), "true");
    assert_eq!(json_scalar(&value("null")), "null");
    assert_eq!(json_scalar(&Value::Float(0.00001)), "1e-5");
    assert_eq!(json_scalar(&Value::Float(1e16)), "1e16");
    assert_eq!(json_scalar(&Value::Float(0.5)), "0.5");
    assert_eq!(
        json_scalar(&value("{\"name\": \"PF1\", \"db\": \"Pfam\"}")),
        "{\"db\":\"Pfam\",\"name\":\"PF1\"}"
    );
    assert_eq!(
        json_scalar(&value("[{\"id\": \"rs1\", \"af\": 0.25}]")),
        "{\"af\":0.25,\"id\":\"rs1\"}"
    );
}

#[test]
fn a_float_and_a_value_have_one_compact_text() {
    assert_eq!(_json_float(1e22), "1e22");
    assert_eq!(_json_float(1.5e-7), "1.5e-7");
    assert_eq!(_json_float(2.0), "2.0");
    assert_eq!(_json_float(f64::NAN), "null");
    assert_eq!(_json_float(f64::NEG_INFINITY), "null");
    assert_eq!(
        json_compact(&value("{\"b\": [1, 2.0, null, true], \"a\": \"\u{e9}\"}")),
        "{\"a\":\"\u{e9}\",\"b\":[1,2.0,null,true]}"
    );
}

#[test]
fn an_object_names_its_record_by_the_first_five_columns_of_its_input() {
    let mut inputs = planted_inputs();
    inputs.push(String::new());
    let mut index = InputIndex::new(&inputs, &mut vcf_line_key);
    let first = value(FIRST_OBJECT);
    assert_eq!(
        json_record_label(&mut index, 5, &first),
        (Some(0), "rec0".to_string())
    );
    assert_eq!(
        json_record_label(&mut index, 5, &first),
        (None, "obj5".to_string())
    );
    // An object without a text `input` asks for the empty key, which an empty input line has.
    assert_eq!(
        json_record_label(&mut index, 6, &value(THIRD_OBJECT)),
        (Some(2), "rec2".to_string())
    );
    assert_eq!(
        json_record_label(&mut index, 7, &value("{}")),
        (None, "obj7".to_string())
    );
}

#[test]
fn an_object_without_two_whole_bounds_has_no_span() {
    for object in [
        "{}",
        "{\"start\": 5}",
        "{\"end\": 5}",
        "{\"start\": true, \"end\": 5}",
        "{\"start\": 5, \"end\": 5.0}",
        "{\"start\": \"5\", \"end\": 5}",
        "{\"start\": 5, \"end\": null}",
    ] {
        assert_eq!(json_object_span(&value(object)), None, "{object}");
    }
}

#[test]
fn an_object_with_two_whole_bounds_has_its_start_and_its_location() {
    assert_eq!(
        json_object_span(&value(
            "{\"start\": 101, \"end\": 100, \"seq_region_name\": \"1\"}"
        )),
        Some((101, "1:100-101".to_string()))
    );
    assert_eq!(
        json_object_span(&value(
            "{\"start\": 7, \"end\": 7, \"seq_region_name\": \"X\"}"
        )),
        Some((7, "X:7".to_string()))
    );
}

#[test]
fn a_line_that_is_no_object_is_an_error_that_numbers_it_in_the_whole_output() {
    let lines = strings(&["{\"a\": 1}", "", "[1, 2]"]);
    assert_eq!(
        comparator_error(parse_json_lines(&lines, &[0, 1, 2], &[])),
        "JSON line 3 is not an object"
    );
    let lines = strings(&["{\"a\": 1}", "{\"a\": 1,}"]);
    assert_eq!(
        comparator_error(parse_json_lines(&lines, &[0, 6], &[])),
        "JSON line 7 does not parse: Illegal trailing comma before end of object: line 1 column 8 \
         (char 7)"
    );
}

#[test]
fn an_item_of_a_consequence_list_that_is_no_object_is_an_error() {
    let lines = strings(&["{\"transcript_consequences\": [{}, 5]}"]);
    assert_eq!(
        comparator_error(parse_json_lines(&lines, &[2], &[])),
        "obj2: a transcript_consequences element is not an object"
    );
    let lines =
        strings(&[r#"{"input": "1\t100\trs1\tG\tGTAT,C", "motif_feature_consequences": [null]}"#]);
    assert_eq!(
        comparator_error(parse_json_lines(&lines, &[2], &planted_inputs())),
        "rec0: a motif_feature_consequences element is not an object"
    );
}

#[test]
fn a_later_line_that_is_no_object_is_reported_before_an_item_that_is_no_object() {
    let lines = strings(&["{\"transcript_consequences\": [5]}", "{\"a\": 1}", "{bad"]);
    assert_eq!(
        comparator_error(parse_json_lines(&lines, &[0, 1, 2], &[])),
        "JSON line 3 does not parse: Expecting property name enclosed in double quotes: line 1 \
         column 2 (char 1)"
    );
    let lines = strings(&["{\"transcript_consequences\": [5]}", "[]"]);
    assert_eq!(
        comparator_error(parse_json_lines(&lines, &[0, 1], &[])),
        "JSON line 2 is not an object"
    );
    let lines = strings(&[&nested(600, "1"), "{bad"]);
    assert_eq!(
        comparator_error(parse_json_lines(&lines, &[0, 1], &[])),
        "JSON line 2 does not parse: Expecting property name enclosed in double quotes: line 1 \
         column 2 (char 1)"
    );
}

#[test]
fn a_line_without_an_ordinal_is_an_error() {
    let lines = strings(&["{}", "{}"]);
    assert_eq!(
        comparator_error(parse_json_lines(&lines, &[0], &[])),
        "a line has no ordinal in its output"
    );
}

#[test]
fn half_a_surrogate_pair_in_a_string_is_refused() {
    let lines = strings(&["{\"a\": \"\\ud800\"}"]);
    assert_eq!(
        comparator_error(parse_json_lines(&lines, &[0], &[])),
        "JSON line 1 does not parse: a string holds half a surrogate pair (char 7)"
    );
}

#[test]
fn a_value_nested_deeper_than_the_interpreter_walks_ends_the_run() {
    // In a comparison a whole number inside 497 lists is one frame too many, and so is a
    // text there.
    for leaf in ["1", "\"s\""] {
        let lines = strings(&[&nested(497, leaf)]);
        assert_eq!(
            uncaught(|| parse_json_lines(&lines, &[0], &[])),
            "maximum recursion depth exceeded",
            "{leaf}"
        );
    }
    // The members of a consequence item and the allele string are read at the same depth.
    for object in [
        format!(
            "{{\"transcript_consequences\": [{{\"x\": {}}}]}}",
            lists(497, "1")
        ),
        format!("{{\"allele_string\": {}}}", lists(497, "1")),
    ] {
        assert_eq!(
            uncaught(|| parse_json_lines(&[object], &[0], &[])),
            "maximum recursion depth exceeded"
        );
    }
    let lines = strings(&[&format!("{}{}", "[".repeat(200_000), "]".repeat(200_000))]);
    assert!(
        uncaught(|| parse_json_lines(&lines, &[0], &[])).contains("while decoding a JSON array")
    );
}

#[test]
fn a_value_one_list_short_of_the_limit_is_read() {
    for (leaf, depth) in [("1", 496), ("\"s\"", 496)] {
        let lines = strings(&[&nested(depth, leaf)]);
        let entries = parse_json_lines(&lines, &[0], &[]).unwrap();
        assert_eq!(keys(&entries), ["obj0|record"], "{leaf}");
    }
    // The self-parity mode calls the parser one frame higher: a text inside 497 lists
    // is read there, and inside 498 it is not.
    let text = format!("{}\n", nested(497, "\"s\""));
    let entries = parse_json_under(&text, &[], SELF_PARITY_CALLERS).unwrap();
    assert_eq!(fields(&entries[0]), [("a", "s")]);
    let text = format!("{}\n", nested(498, "\"s\""));
    assert_eq!(
        uncaught(|| parse_json_under(&text, &[], SELF_PARITY_CALLERS)),
        "maximum recursion depth exceeded"
    );
    assert_eq!(
        uncaught(|| parse_json(&format!("{}\n", nested(497, "\"s\"")), &[])),
        "maximum recursion depth exceeded"
    );
}

/// The parser as one of the two modes calls it.
type Parse<'a> = &'a dyn Fn(&str) -> Result<Vec<Entry>, Error>;

/// The Location of an object with two whole bounds is written two frames above the
/// parser, so its `seq_region_name` is refused one list sooner than any other member: 496
/// lists in a comparison and 497 in the self-parity mode, where 497 and 498 are refused
/// without a span. A name inside objects is refused at 495 of them in both modes, and read
/// there without a span.
#[test]
fn a_seq_region_name_is_refused_sooner_where_the_object_has_a_span() {
    let objects = |depth: usize| format!("{}\"s\"{}", "{\"k\": ".repeat(depth), "}".repeat(depth));
    let object = |name: String, start: &str| {
        format!("{{\"seq_region_name\": {name}, \"start\": {start}, \"end\": 100}}\n")
    };
    let comparison = |text: &str| parse_json(text, &[]);
    let self_parity = |text: &str| parse_json_under(text, &[], SELF_PARITY_CALLERS);
    let read: [(Parse, String); 8] = [
        (&comparison, object(lists(495, "\"s\""), "100")),
        (&comparison, object(lists(496, "\"s\""), "100.5")),
        (&comparison, object(objects(494), "100")),
        (&comparison, object(objects(495), "100.5")),
        (&self_parity, object(lists(496, "\"s\""), "100")),
        (&self_parity, object(lists(497, "\"s\""), "100.5")),
        (&self_parity, object(objects(494), "100")),
        (&self_parity, object(objects(495), "100.5")),
    ];
    for (parse, text) in &read {
        let entries = parse(text).unwrap();
        assert_eq!(keys(&entries), ["obj0|record"]);
    }
    let refused: [(Parse, String); 4] = [
        (&comparison, object(lists(496, "\"s\""), "100")),
        (&comparison, object(objects(495), "100")),
        (&self_parity, object(lists(497, "\"s\""), "100")),
        (&self_parity, object(objects(495), "100")),
    ];
    for (parse, text) in &refused {
        let ended = uncaught(std::panic::AssertUnwindSafe(|| parse(text)));
        assert_eq!(ended, "maximum recursion depth exceeded");
    }
}

#[test]
fn the_objects_of_a_json_lines_output_are_numbered_by_their_lines() {
    let objects = json_objects("{\"a\": 1}\n\n  \n{\"b\": 2}\r\n").unwrap();
    assert_eq!(
        objects,
        [(0, value("{\"a\": 1}")), (3, value("{\"b\": 2}"))]
    );
    assert_eq!(json_objects(""), Ok(Vec::new()));
    assert_eq!(
        json_objects("\n7\n"),
        Err(Error::Comparator(
            "JSON line 2 is not an object".to_string()
        ))
    );
}

#[test]
fn an_object_gives_a_record_entry_and_an_entry_an_item_of_its_consequence_lists() {
    let text = format!("{FIRST_OBJECT}\n\n   \n{SECOND_OBJECT}\n{THIRD_OBJECT}\n");
    let entries = parse_json(&text, &planted_inputs()).unwrap();
    assert_eq!(
        keys(&entries),
        [
            "rec0|record",
            "rec0|TAT|T1|transcript",
            "rec0|TAT|T1|transcript#2",
            "rec0|C||intergenic",
            "rec1|record",
            "rec1|-|R1|regulatory_feature",
            "rec1|-|M1|motif_feature",
            "obj4|record",
            "obj4|T||transcript",
        ]
    );
    assert_eq!(
        fields(&entries[0]),
        [
            ("start", "101"),
            ("end", "100"),
            ("allele_string", "-/TAT/C"),
            ("seq_region_name", "1"),
            ("id", "rs1"),
            ("most_severe_consequence", "null"),
            (
                "colocated_variants",
                "{\"f\":0.25,\"id\":\"x\"},{\"id\":\"y\"}"
            ),
        ]
    );
    assert_eq!(
        fields(&entries[1]),
        [
            ("variant_allele", "TAT"),
            ("transcript_id", "T1"),
            ("consequence_terms", "a,b"),
            ("n", "1.5e-5"),
            ("strand", "-1"),
        ]
    );
    assert_eq!(
        fields(&entries[3]),
        [
            ("variant_allele", "C"),
            ("consequence_terms", "i"),
            ("intergenic_id", "never read"),
        ]
    );
    assert_eq!(
        fields(&entries[4]),
        [
            ("start", "201"),
            ("end", "201"),
            ("allele_string", "T/-"),
            ("seq_region_name", "1"),
        ]
    );
    assert_eq!(
        fields(&entries[7]),
        [("allele_string", "A/T"), ("start", "true"), ("end", "5")]
    );
    let read: Vec<_> = entries
        .iter()
        .map(|entry| {
            (
                entry.record,
                entry.allele.as_str(),
                entry.feature.as_str(),
                entry.consequence.as_str(),
                entry.indel,
                entry.insertion,
                entry.reference.as_deref(),
                entry.allele_location.as_deref(),
            )
        })
        .collect();
    let (yes, no) = (Some(true), Some(false));
    assert_eq!(
        read,
        [
            (Some(0), "", "", "", None, None, None, None),
            (
                Some(0),
                "TAT",
                "T1",
                "a,b",
                yes,
                yes,
                Some(""),
                Some("1:100-101")
            ),
            (
                Some(0),
                "TAT",
                "T1",
                "",
                yes,
                yes,
                Some(""),
                Some("1:100-101")
            ),
            (Some(0), "C", "", "i", yes, yes, Some(""), Some("1:100-101")),
            (Some(1), "", "", "", None, None, None, None),
            (Some(1), "-", "R1", "", yes, no, Some("T"), Some("1:201")),
            (Some(1), "-", "M1", "", yes, no, Some("T"), Some("1:201")),
            (None, "", "", "", None, None, None, None),
            (None, "T", "", "", no, no, None, None),
        ]
    );
    assert!(entries
        .iter()
        .all(|entry| entry.line.is_none() && entry.location.is_none()));
}

#[test]
fn an_object_that_matches_no_input_record_is_named_by_its_ordinal_in_the_whole_output() {
    let lines = strings(&[SECOND_OBJECT, "", THIRD_OBJECT]);
    let entries = parse_json_lines(&lines, &[3, 9, 20], &planted_inputs()).unwrap();
    assert_eq!(
        keys(&entries),
        [
            "rec1|record",
            "rec1|-|R1|regulatory_feature",
            "rec1|-|M1|motif_feature",
            "obj20|record",
            "obj20|T||transcript",
        ]
    );
}
