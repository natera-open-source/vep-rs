// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The alleles of an individual's copy compared as a set.

mod compare_common;

use std::collections::BTreeSet;

use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::classes::INDIVIDUAL_ALLELE_ORDER;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::compare_files::compare_files_through;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::input::input_records;
use vep_formats_compare::input::read_text;
use vep_formats_compare::input::text_lines;
use vep_formats_compare::lists::order_only_fields;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::Side;

/// The hand-built hb_multi record (G, ALT A,T) and its sample S2, genotype 1/2: the copy's allele
/// string is G/A/T and each of its four transcripts has one A row and one T row, in that order.
const TWO_ALLELE_ROWS: i64 = 8;

fn default_entries(text: &str) -> Vec<Entry> {
    parse_default(text).expect("the text parses")
}

fn with_class() -> Documented {
    documented_with(&[INDIVIDUAL_ALLELE_ORDER])
}

/// The rows of hb_multi's copy of the sample S2.
fn s2_rows(rows: &[String]) -> Vec<usize> {
    (0..rows.len())
        .filter(|i| rows[*i].starts_with("hb_multi\t") && rows[*i].contains("IND=S2;"))
        .collect()
}

fn alleles(rows: &[String], targets: &[usize]) -> Vec<String> {
    targets
        .iter()
        .map(|i| rows[*i].split('\t').nth(2).unwrap_or_default().to_string())
        .collect()
}

fn repeated(items: &[&str], times: usize) -> Vec<String> {
    (0..times)
        .flat_map(|_| items.iter().map(|item| item.to_string()))
        .collect()
}

fn is_s2_of_hb_multi(object: &Value) -> bool {
    object.get("id").and_then(Value::as_str) == Some("hb_multi")
        && object.get("ind").and_then(Value::as_str) == Some("S2")
}

fn consequences_mut(object: &mut Value) -> Option<&mut Vec<Value>> {
    match object {
        Value::Dict(members) => match members.get_mut("transcript_consequences") {
            Some(Value::List(items)) => Some(items),
            _ => None,
        },
        _ => None,
    }
}

fn consequence_alleles(object: &Value) -> Vec<String> {
    object
        .get("transcript_consequences")
        .and_then(Value::as_list)
        .unwrap_or_default()
        .iter()
        .map(|item| {
            item.get("variant_allele")
                .and_then(Value::as_str)
                .unwrap_or_default()
                .to_string()
        })
        .collect()
}

fn dumps(object: &Value) -> String {
    pyjson::dumps(object, &pyjson::DumpOptions::default())
}

fn objects_of(text: &str) -> Vec<Value> {
    text_lines(text)
        .iter()
        .filter(|line| !line.is_empty())
        .map(|line| json(line))
        .collect()
}

fn joined(objects: &[String]) -> String {
    objects.join("\n") + "\n"
}

/// Under `--individual` Ensembl VEP lists the alternate alleles of an individual's copy in
/// hash order (Parser/VCF.pm:639), so a candidate listing hb_multi S2's T before its A agrees
/// under the class in every format: each transcript's two rows in the other order, the JSON
/// record's `allele_string`, and under `--allele_number --uploaded_allele` the rows'
/// `ALLELE_NUM` (one permutation of the alternates' numbers over every row of the copy) and
/// `UPLOADED_ALLELE`; each agreement is counted as `individual_allele_order` and the
/// reordered entries under `counts`. The same rows without `IND` stay red, as do a changed
/// allele, a numbering that is not a permutation, a renumbering of one transcript alone, a
/// reference allele numbered other than 0 and a reference row moved after the alternates;
/// `order_only` alone accepts none of it.
#[test]
fn the_individual_allele_order_class_accepts_a_copys_alleles_in_any_order() {
    let directory = scratch("individual_allele_order");
    let individual = golden().join("GRCh38-individual");
    let expected_dir = individual.join("expected_individual-all");
    let inputs = input_records(&individual.join("variants.vcf")).expect("the input reads");
    let default = read_text(&expected_dir.join("default.txt.gz")).expect("the output reads");
    let rows = data_lines(&default);
    let targets = s2_rows(&rows);
    assert_eq!(targets.len(), 8);
    assert_eq!(alleles(&rows, &targets), repeated(&["A", "T"], 4));
    let mut swapped = rows.clone();
    for pair in targets.chunks(2) {
        swapped.swap(pair[0], pair[1]);
    }
    let expected = default_entries(&default);
    let actual = default_entries(&rebuild(&default, &swapped));
    let first_place = format!("order differs at shared position {}", targets[0]);
    let plain = compare(&expected, &actual, &mut no_documented());
    assert_eq!(verdict(&plain), "red");
    assert_eq!(reasons(&plain), [first_place.as_str()]);
    assert_eq!(
        at(&plain, &["order", "expected"]).as_str(),
        Some("21:34716340|A|ENST00000349499|Transcript#2")
    );
    assert_eq!(
        at(&plain, &["order", "actual"]).as_str(),
        Some("21:34716340|T|ENST00000349499|Transcript")
    );
    let side = Side::default();
    let lists = order_only_fields("default");
    let order_only = compare_entries(
        &expected,
        &actual,
        &mut no_documented(),
        &CompareOptions {
            order_only: Some(&lists),
            ..CompareOptions::new(&side)
        },
    )
    .expect("the entries compare");
    assert_eq!(verdict(&order_only), "red");
    let mut documented = with_class();
    let accepted = compare(&expected, &actual, &mut documented);
    assert_eq!(verdict(&accepted), "green");
    assert_eq!(at(&accepted, &["order"]), &Value::Null);
    assert_eq!(
        number(&accepted, &["counts", "individual_entries_reordered"]),
        TWO_ALLELE_ROWS
    );
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(INDIVIDUAL_ALLELE_ORDER, TWO_ALLELE_ROWS)])
    );
    let columns = at(&accepted, &["columns"]).as_dict().expect("an object");
    assert!(columns
        .values()
        .all(|column| column.get("individual_allele_order") == Some(&Value::Int(0))));

    // The same eight rows without IND are a record of no individual: the order stays red.
    let without_individual = |rows: &[String]| -> Vec<Entry> {
        let eight: Vec<String> = targets
            .iter()
            .map(|i| rows[*i].replace("IND=S2;", ""))
            .collect();
        default_entries(&rebuild(&default, &eight))
    };
    let mut documented = with_class();
    let scoped = compare(
        &without_individual(&rows),
        &without_individual(&swapped),
        &mut documented,
    );
    assert_eq!(verdict(&scoped), "red");
    assert_eq!(reasons(&scoped), ["order differs at shared position 0"]);
    assert_eq!(
        number(&scoped, &["counts", "individual_entries_reordered"]),
        0
    );
    assert!(counts(&documented.applied).is_empty());

    // Compared as two files the class is reported applied, and unused on a candidate in the
    // reference's order.
    let reference_path = planted(&directory, "reference.txt", &default);
    let candidate_path = planted(&directory, "candidate.txt", &rebuild(&default, &swapped));
    let files = |candidate: &std::path::Path| -> Value {
        compare_files_through(
            None,
            "default",
            &reference_path,
            candidate,
            None,
            &mut with_class(),
            &masks(&[INDIVIDUAL_ALLELE_ORDER]),
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
    };
    let report = files(&candidate_path);
    assert_eq!(verdict(&report), "green");
    assert_eq!(
        at(&report, &["masks"]),
        &json(&format!(
            r#"{{"requested": ["{INDIVIDUAL_ALLELE_ORDER}"], "applied": {{"{INDIVIDUAL_ALLELE_ORDER}": {TWO_ALLELE_ROWS}}}, "unused": [], "available": {{}}}}"#
        ))
    );
    let same = files(&reference_path);
    assert_eq!(verdict(&same), "green");
    assert_eq!(at(&same, &["masks", "applied"]), &json("{}"));
    assert_eq!(
        texts(at(&same, &["masks", "unused"])),
        [INDIVIDUAL_ALLELE_ORDER]
    );

    // The JSON record object's allele_string (JSON.pm:241) and its consequence entries.
    let json_text = read_text(&expected_dir.join("json.jsonl.gz")).expect("the output reads");
    let mut objects: Vec<String> = Vec::new();
    for mut object in objects_of(&json_text) {
        if is_s2_of_hb_multi(&object) {
            assert_eq!(
                object.get("allele_string").and_then(Value::as_str),
                Some("G/A/T")
            );
            assert_eq!(consequence_alleles(&object), repeated(&["A", "T"], 4));
            if let Value::Dict(members) = &mut object {
                members.insert("allele_string".to_string(), Value::Str("G/T/A".to_string()));
            }
            let consequences = consequences_mut(&mut object).expect("the record has rows");
            for pair in consequences.chunks_mut(2) {
                pair.swap(0, 1);
            }
        }
        objects.push(dumps(&object));
    }
    let json_entries =
        |text: &str| -> Vec<Entry> { parse_json(text, &inputs).expect("the text parses") };
    let json_expected = json_entries(&json_text);
    let json_actual = json_entries(&joined(&objects));
    let plain = compare(&json_expected, &json_actual, &mut no_documented());
    assert_eq!(verdict(&plain), "red");
    assert_eq!(
        number(&plain, &["columns", "allele_string", "differing"]),
        1
    );
    let mut documented = with_class();
    let json_report = compare(&json_expected, &json_actual, &mut documented);
    assert_eq!(verdict(&json_report), "green");
    assert_eq!(
        number(
            &json_report,
            &["columns", "allele_string", "individual_allele_order"]
        ),
        1
    );
    assert_eq!(
        texts(at(
            &json_report,
            &[
                "columns",
                "allele_string",
                "individual_allele_order_examples"
            ]
        )),
        ["obj434|record: expected 'G/A/T', got 'G/T/A'"]
    );
    assert_eq!(
        number(&json_report, &["counts", "individual_entries_reordered"]),
        TWO_ALLELE_ROWS
    );
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(INDIVIDUAL_ALLELE_ORDER, TWO_ALLELE_ROWS + 1)])
    );
    let changed: Vec<String> = objects
        .iter()
        .map(|o| o.replace(r#""allele_string": "G/T/A""#, r#""allele_string": "G/T/C""#))
        .collect();
    let changed_report = compare(
        &json_expected,
        &json_entries(&joined(&changed)),
        &mut with_class(),
    );
    assert_eq!(verdict(&changed_report), "red");
    assert_eq!(
        reasons(&changed_report),
        ["field allele_string: 1 differing entries"]
    );

    // The VCF format's CSQ entries of the copy's line, each transcript's pair in the other order.
    let vcf = read_text(&expected_dir.join("vcf.vcf.gz")).expect("the output reads");
    let flipped: Vec<String> = data_lines(&vcf)
        .iter()
        .map(|row| {
            let mut cols: Vec<String> = row.split('\t').map(str::to_string).collect();
            if cols[2] == "hb_multi" && cols[7].contains("|S2|") {
                let (prefix, csq) = cols[7]
                    .split_once("CSQ=")
                    .expect("the line carries its entries");
                let mut chunks: Vec<&str> = csq.split(',').collect();
                let first: Vec<&str> = chunks
                    .iter()
                    .map(|chunk| chunk.split('|').next().unwrap_or_default())
                    .collect();
                assert_eq!(first, repeated(&["A", "T"], 4));
                for pair in chunks.chunks_mut(2) {
                    pair.swap(0, 1);
                }
                cols[7] = format!("{prefix}CSQ={}", chunks.join(","));
            }
            cols.join("\t")
        })
        .collect();
    let vcf_expected = parse_vcf(&vcf, &inputs).expect("the reference parses");
    let vcf_actual = parse_vcf(&rebuild(&vcf, &flipped), &inputs).expect("the candidate parses");
    let plain = compare(&vcf_expected, &vcf_actual, &mut no_documented());
    assert_eq!(verdict(&plain), "red");
    let mut documented = with_class();
    let vcf_report = compare(&vcf_expected, &vcf_actual, &mut documented);
    assert_eq!(verdict(&vcf_report), "green");
    assert_eq!(
        number(&vcf_report, &["counts", "individual_entries_reordered"]),
        TWO_ALLELE_ROWS
    );
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(INDIVIDUAL_ALLELE_ORDER, TWO_ALLELE_ROWS)])
    );

    // --process_ref_homs --allele_number --uploaded_allele: the copy's rows are G, A and T per
    // transcript, numbered 0, 1 and 2 in the allele string's order (OutputFactory.pm:1354),
    // UPLOADED_ALLELE the string itself (OutputFactory.pm:1360); the other order numbers T 1
    // and A 2 and prints G/T/A on every row of the copy, the G row included.
    let fields_dir = individual.join("expected_individual-refhoms-fields");
    let fields_default = read_text(&fields_dir.join("default.txt.gz")).expect("the output reads");
    let rows = data_lines(&fields_default);
    let targets = s2_rows(&rows);
    let copy_rows = i64::try_from(targets.len()).expect("a count");
    assert_eq!(alleles(&rows, &targets), repeated(&["G", "A", "T"], 4));
    let mut renumbered = rows.clone();
    for i in &targets {
        let row = rows[*i].replace("UPLOADED_ALLELE=G/A/T;", "UPLOADED_ALLELE=G/T/A;");
        renumbered[*i] = if row.contains("\tA\t") {
            row.replace("ALLELE_NUM=1;", "ALLELE_NUM=2;")
        } else {
            row.replace("ALLELE_NUM=2;", "ALLELE_NUM=1;")
        };
    }
    let reference_rows: Vec<usize> = targets.iter().copied().step_by(3).collect();
    for g in &reference_rows {
        renumbered.swap(g + 1, g + 2);
    }
    let fields_expected = default_entries(&fields_default);
    let against = |rows: &[String], documented: &mut Documented| -> Value {
        compare(
            &fields_expected,
            &default_entries(&rebuild(&fields_default, rows)),
            documented,
        )
    };
    let plain = against(&renumbered, &mut no_documented());
    assert_eq!(verdict(&plain), "red");
    assert_eq!(
        number(&plain, &["columns", "ALLELE_NUM", "differing"]),
        TWO_ALLELE_ROWS
    );
    assert_eq!(
        number(&plain, &["columns", "UPLOADED_ALLELE", "differing"]),
        copy_rows
    );
    let mut documented = with_class();
    let fields_report = against(&renumbered, &mut documented);
    assert_eq!(verdict(&fields_report), "green");
    let accepted_under = |column: &str| {
        number(
            &fields_report,
            &["columns", column, "individual_allele_order"],
        )
    };
    assert_eq!(accepted_under("ALLELE_NUM"), TWO_ALLELE_ROWS);
    assert_eq!(accepted_under("UPLOADED_ALLELE"), copy_rows);
    assert_eq!(accepted_under("Extra"), copy_rows);
    assert_eq!(
        number(&fields_report, &["counts", "individual_entries_reordered"]),
        TWO_ALLELE_ROWS
    );
    // The eight reordered rows, their eight ALLELE_NUM values, and UPLOADED_ALLELE with the
    // Extra text carrying it on all twelve rows of the copy.
    assert_eq!(
        counts(&documented.applied),
        count_map(&[(INDIVIDUAL_ALLELE_ORDER, 2 * TWO_ALLELE_ROWS + 2 * copy_rows)])
    );
    let examples = texts(at(
        &fields_report,
        &["columns", "ALLELE_NUM", "individual_allele_order_examples"],
    ));
    assert_eq!(
        examples[..2],
        [
            "21:34716340|A|ENST00000349499|Transcript#2: expected '1', got '2'",
            "21:34716340|T|ENST00000349499|Transcript: expected '2', got '1'",
        ]
    );
    // Numbering both alternate alleles 1 is not a permutation of the copy's numbers.
    let in_copy: BTreeSet<usize> = targets.iter().copied().collect();
    let doubled: Vec<String> = renumbered
        .iter()
        .enumerate()
        .map(|(i, r)| {
            if in_copy.contains(&i) {
                r.replace("ALLELE_NUM=2;", "ALLELE_NUM=1;")
            } else {
                r.clone()
            }
        })
        .collect();
    let doubled_report = against(&doubled, &mut with_class());
    assert_eq!(verdict(&doubled_report), "red");
    assert_eq!(
        reasons(&doubled_report),
        [
            "field ALLELE_NUM: 4 differing entries",
            "field Extra: 4 differing entries"
        ]
    );
    assert_eq!(
        number(
            &doubled_report,
            &["columns", "ALLELE_NUM", "individual_allele_order"]
        ),
        0
    );
    // A renumbering of one transcript alone is no permutation applied to the copy: each allele
    // carries one number across the copy's rows, its index in the copy's one allele string
    // (OutputFactory.pm:1354).
    let mut partial = rows.clone();
    partial[targets[1]] = rows[targets[1]].replace("ALLELE_NUM=1;", "ALLELE_NUM=2;");
    partial[targets[2]] = rows[targets[2]].replace("ALLELE_NUM=2;", "ALLELE_NUM=1;");
    let partial_report = against(&partial, &mut with_class());
    assert_eq!(verdict(&partial_report), "red");
    assert_eq!(
        number(&partial_report, &["columns", "ALLELE_NUM", "differing"]),
        2
    );
    assert_eq!(
        number(
            &partial_report,
            &["columns", "ALLELE_NUM", "individual_allele_order"]
        ),
        0
    );
    // The reference allele keeps 0 and its row keeps its place (Parser/VCF.pm:639 writes the
    // reference first): G numbered 1 with A numbered 0, or the G row after the alternates, is red.
    let mut shifted_zero = rows.clone();
    for g in &reference_rows {
        shifted_zero[*g] = rows[*g].replace("ALLELE_NUM=0;", "ALLELE_NUM=1;");
        shifted_zero[g + 1] = rows[g + 1].replace("ALLELE_NUM=1;", "ALLELE_NUM=0;");
    }
    let zero_report = against(&shifted_zero, &mut with_class());
    assert_eq!(verdict(&zero_report), "red");
    assert_eq!(
        number(&zero_report, &["columns", "ALLELE_NUM", "differing"]),
        TWO_ALLELE_ROWS
    );
    assert_eq!(
        number(
            &zero_report,
            &["columns", "ALLELE_NUM", "individual_allele_order"]
        ),
        0
    );
    let reference_after_alternates = |rows: &[String], groups: &[usize]| -> Vec<String> {
        let mut moved = rows.to_vec();
        for g in groups {
            moved[*g] = rows[g + 1].clone();
            moved[g + 1] = rows[g + 2].clone();
            moved[g + 2] = rows[*g].clone();
        }
        moved
    };
    let last_report = against(
        &reference_after_alternates(&rows, &reference_rows),
        &mut with_class(),
    );
    assert_eq!(verdict(&last_report), "red");
    let first_place = format!("order differs at shared position {}", targets[0]);
    assert_eq!(reasons(&last_report), [first_place.as_str()]);

    // A JSON consequence entry names the reference allele as ref_allele: without allele_num the
    // reference row is still known by it, so moved after the alternates it reads red while the
    // alternates' own order stays free.
    let reordered = |object: &Value, order: [usize; 3]| -> String {
        let mut moved = object.clone();
        let consequences = consequences_mut(&mut moved).expect("the record has rows");
        let held = consequences.clone();
        let picked: Vec<Value> = (0..held.len())
            .step_by(3)
            .flat_map(|i| order.iter().map(move |j| i + j))
            .map(|at| held[at].clone())
            .collect();
        *consequences = picked;
        dumps(&moved)
    };
    let fields_json = read_text(&fields_dir.join("json.jsonl.gz")).expect("the output reads");
    let mut plain_objects: Vec<String> = Vec::new();
    let mut reference_last_objects: Vec<String> = Vec::new();
    let mut alternates_objects: Vec<String> = Vec::new();
    for mut object in objects_of(&fields_json) {
        if let Some(consequences) = consequences_mut(&mut object) {
            for consequence in consequences {
                if let Value::Dict(held) = consequence {
                    held.shift_remove("allele_num");
                }
            }
        }
        let text = dumps(&object);
        plain_objects.push(text.clone());
        if is_s2_of_hb_multi(&object) {
            assert_eq!(consequence_alleles(&object), repeated(&["G", "A", "T"], 4));
            let consequences = object
                .get("transcript_consequences")
                .and_then(Value::as_list)
                .expect("the record has rows");
            assert!(consequences
                .iter()
                .all(|item| item.get("ref_allele").and_then(Value::as_str) == Some("G")));
            reference_last_objects.push(reordered(&object, [1, 2, 0]));
            alternates_objects.push(reordered(&object, [0, 2, 1]));
        } else {
            reference_last_objects.push(text.clone());
            alternates_objects.push(text);
        }
    }
    let json_fields_expected = json_entries(&joined(&plain_objects));
    let reference_last_json = compare(
        &json_fields_expected,
        &json_entries(&joined(&reference_last_objects)),
        &mut with_class(),
    );
    assert_eq!(verdict(&reference_last_json), "red");
    assert!(reasons(&reference_last_json)[0].starts_with("order differs at shared position"));
    let alternates_json = compare(
        &json_fields_expected,
        &json_entries(&joined(&alternates_objects)),
        &mut with_class(),
    );
    assert_eq!(verdict(&alternates_json), "green");
    assert_eq!(
        number(
            &alternates_json,
            &["counts", "individual_entries_reordered"]
        ),
        TWO_ALLELE_ROWS
    );
    // A number missing on one side is no permutation: the row's own value must be its allele's
    // number in the copy.
    let mut one_missing: Vec<String> = Vec::new();
    for mut object in objects_of(&fields_json) {
        if is_s2_of_hb_multi(&object) {
            let consequences = consequences_mut(&mut object).expect("the record has rows");
            let Value::Dict(first) = &mut consequences[0] else {
                panic!("a consequence entry is an object");
            };
            assert_eq!(first.shift_remove("allele_num"), Some(Value::Int(0)));
        }
        one_missing.push(dumps(&object));
    }
    let one_missing_report = compare(
        &json_entries(&fields_json),
        &json_entries(&joined(&one_missing)),
        &mut with_class(),
    );
    assert_eq!(verdict(&one_missing_report), "red");
    assert_eq!(
        reasons(&one_missing_report),
        ["field allele_num: 1 differing entries"]
    );
    assert_eq!(
        number(
            &one_missing_report,
            &["columns", "allele_num", "individual_allele_order"]
        ),
        0
    );
    // An insertion's reference row prints REF_ALLELE=-, its own allele: without a number it is
    // still the reference row, and moved after the alternates it reads red.
    let indel_targets: Vec<usize> = (0..rows.len())
        .filter(|i| rows[*i].starts_with("hb_multi_indel\t") && rows[*i].contains("IND=S4;"))
        .collect();
    assert_eq!(alleles(&rows, &indel_targets)[..3], ["-", "A", "AA"]);
    let allele_number = regex::Regex::new(r"ALLELE_NUM=\d+;").expect("the pattern compiles");
    let in_indel: BTreeSet<usize> = indel_targets.iter().copied().collect();
    let unnumbered: Vec<String> = rows
        .iter()
        .enumerate()
        .map(|(i, r)| {
            if in_indel.contains(&i) {
                allele_number.replace_all(r, "").into_owned()
            } else {
                r.clone()
            }
        })
        .collect();
    let indel_groups: Vec<usize> = indel_targets.iter().copied().step_by(3).collect();
    let indel_last = reference_after_alternates(&unnumbered, &indel_groups);
    let indel_report = compare(
        &default_entries(&rebuild(&fields_default, &unnumbered)),
        &default_entries(&rebuild(&fields_default, &indel_last)),
        &mut with_class(),
    );
    assert_eq!(verdict(&indel_report), "red");
    let first_place = format!("order differs at shared position {}", indel_targets[0]);
    assert_eq!(reasons(&indel_report), [first_place.as_str()]);
}
