// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `total_length_json_range_end` through the comparison of two JSON outputs: the
//! end key of a position range under `--total_length`, on planted objects and on the stored
//! reference of the committed corpus that runs the flag.

mod compare_common;

use std::collections::BTreeMap;

use compare_common::*;
use vep_formats_compare::classes::TOTAL_LENGTH_JSON_RANGE_END;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::input::input_records;
use vep_formats_compare::input::read_text;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;

const INPUT: &str = "21\t100\tdel\tAT\tA\t.\t.\t.";

/// The entries of a JSON output of one record whose one transcript object carries `members`.
fn entries(members: &str) -> Vec<Entry> {
    let text = dumped(&format!(
        r#"{{"input": "21\t100\tdel\tAT\tA\t.\t.\t.", "id": "del", "allele_string": "T/-", "most_severe_consequence": "frameshift_variant", "transcript_consequences": [
            {{"variant_allele": "-", "transcript_id": "ENST1", "consequence_terms": ["frameshift_variant"], {members}}}]}}"#
    )) + "\n";
    parse_json(&text, &[INPUT.to_string()]).expect("the text parses")
}

/// The report's verdict and what was applied of two entry lists under `classes`.
fn under(
    classes: &[&str],
    expected: &[Entry],
    actual: &[Entry],
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare(expected, actual, &mut documented);
    (verdict(&report).to_string(), counts(&documented.applied))
}

const REFERENCE: &str = r#""cdna_start": 484, "cdna_end": 484, "cds_start": 116, "cds_end": 116, "protein_start": 39, "protein_end": 39"#;
const CANDIDATE: &str = r#""cdna_start": 484, "cdna_end": "487/2674", "cds_start": 116, "cds_end": "119/1602", "protein_start": 39, "protein_end": "40/533""#;

/// The three end keys of one object are taken under the class, which is credited once for
/// the entry, and are differences without it.
#[test]
fn the_end_keys_of_a_range_are_taken_under_the_class() {
    let klass = TOTAL_LENGTH_JSON_RANGE_END;
    let (reference, candidate) = (entries(REFERENCE), entries(CANDIDATE));
    assert_eq!(
        under(&[klass], &reference, &candidate),
        ("green".to_string(), count_map(&[(klass, 1)]))
    );
    assert_eq!(under(&[], &reference, &candidate).0, "red");
    assert_eq!(under(&["order_only"], &reference, &candidate).0, "red");
}

/// A range whose first position is unknown: the reference writes neither key and the
/// candidate the end key alone.
#[test]
fn the_end_key_of_a_range_of_unknown_start_is_taken() {
    let klass = TOTAL_LENGTH_JSON_RANGE_END;
    let reference = entries(r#""impact": "HIGH""#);
    let candidate = entries(r#""impact": "HIGH", "cdna_end": "487/2674""#);
    assert_eq!(
        under(&[klass], &reference, &candidate),
        ("green".to_string(), count_map(&[(klass, 1)]))
    );
}

/// What the class leaves a difference: a start the candidate moves, an end at the start, an
/// end without its total, and an end key the candidate leaves out.
#[test]
fn another_end_or_start_stays_a_difference_under_the_class() {
    let klass = TOTAL_LENGTH_JSON_RANGE_END;
    let reference = entries(r#""cdna_start": 484, "cdna_end": 484"#);
    for candidate in [
        r#""cdna_start": 485, "cdna_end": "487/2674""#,
        r#""cdna_start": 484, "cdna_end": "484/2674""#,
        r#""cdna_start": 484, "cdna_end": 487"#,
        r#""cdna_start": 484"#,
        r#""cdna_end": "487/2674""#,
    ] {
        assert_eq!(
            under(&[klass], &reference, &entries(candidate)).0,
            "red",
            "{candidate}"
        );
    }
}

/// The stored reference of the corpus that runs the flag, against itself with every end
/// key of a range rewritten as the text format prints the range: every rewritten object is
/// taken under the class and nothing else differs. The tab output of the same flag set
/// gives each row's three columns.
#[test]
fn the_stored_total_length_reference_is_green_against_its_ranges_under_the_class() {
    let corpus = golden().join("GRCh37-output-flags");
    let expected_dir = corpus.join("expected_total-length");
    let inputs = input_records(&corpus.join("variants.vcf")).expect("the input reads");
    let reference_text =
        read_text(&expected_dir.join("json.jsonl.gz")).expect("the stored output reads");
    let tab_text = read_text(&expected_dir.join("tab.txt.gz")).expect("the stored output reads");
    // (record id, allele, transcript) -> the three position columns of the tab row.
    let mut columns: BTreeMap<(String, String, String), [String; 3]> = BTreeMap::new();
    for line in data_lines(&tab_text) {
        let cells: Vec<&str> = line.split('\t').collect();
        columns.insert(
            (cells[0].into(), cells[2].into(), cells[4].into()),
            [cells[7].into(), cells[8].into(), cells[9].into()],
        );
    }
    let mut rewritten_objects = 0;
    let mut candidate_lines = Vec::new();
    for line in data_lines(&reference_text) {
        let parsed = pyjson::loads(&line).expect("a line of the stored output is JSON");
        let mut record = parsed.into_dict().expect("an object");
        let id = record["id"].as_str().expect("an id").to_string();
        if let Some(Value::List(rows)) = record.get_mut("transcript_consequences") {
            for row in rows {
                let Value::Dict(row) = row else {
                    panic!("a consequence is no object");
                };
                let member = |name: &str| row[name].as_str().expect("a text").to_string();
                let key = (
                    id.clone(),
                    member("variant_allele"),
                    member("transcript_id"),
                );
                let mut rewritten = false;
                for (at, end) in ["cdna_end", "cds_end", "protein_end"].iter().enumerate() {
                    let column = &columns[&key][at];
                    let Some((first, last)) = column.split_once('-') else {
                        continue;
                    };
                    // A range of two known positions with its total: `484-487/2674`, and
                    // neither `484-?/2674` nor the placeholder `-/2674`.
                    let digits = |text: &str| {
                        !text.is_empty() && text.bytes().all(|byte| byte.is_ascii_digit())
                    };
                    let known = last.split_once('/').is_some_and(|(last, _)| digits(last));
                    if !digits(first) || !known {
                        continue;
                    }
                    row.insert(end.to_string(), Value::Str(last.to_string()));
                    rewritten = true;
                }
                rewritten_objects += usize::from(rewritten);
            }
        }
        candidate_lines.push(pyjson::dumps_compact(&Value::Dict(record)));
    }
    let candidate_text = candidate_lines.join("\n") + "\n";
    let reference = parse_json(&reference_text, &inputs).expect("the reference parses");
    let candidate = parse_json(&candidate_text, &inputs).expect("the candidate parses");
    let klass = TOTAL_LENGTH_JSON_RANGE_END;
    assert_eq!(rewritten_objects, 48);
    assert_eq!(
        under(&[klass], &reference, &candidate),
        ("green".to_string(), count_map(&[(klass, 48)]))
    );
    assert_eq!(under(&[], &reference, &candidate).0, "red");
}
