// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The cases on which the golden corpus builder's tests hold the comparison to a stated
//! answer, each asked of the tool: `tests/data/builder_cases/` holds one file a test, one
//! JSON object a line, with the call's arguments and the answer the test states. A text that
//! holds a run of 16 or more capital letters, and an output of more than one row, is written
//! as `{"hex": ...}`, its UTF-8 bytes in hexadecimal.

use std::collections::BTreeSet;
use std::path::Path;
use vep_formats_compare::consts::WITH_REFERENCE_SETS;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::Error;
use vep_formats_compare::manifest::with_reference_sets;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::pyjson::{self, Value};
use vep_formats_compare::rules::coding_end_pairs::AtACodingEnd;
use vep_formats_compare::rules::downstream::field_of;
use vep_formats_compare::rules::hgvsp_stop::{
    stop_codon_protein_comparison, stop_named_on_incomplete_coding_sequence,
};
use vep_formats_compare::rules::term_pairs::{read_entry_references, TERM_PAIR_CLASSES};
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};

/// The cases of the file `name`, one a line.
fn cases(name: &str) -> Vec<Value> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/builder_cases")
        .join(format!("{name}.jsonl"));
    std::fs::read_to_string(path)
        .expect("the cases are read")
        .lines()
        .map(|line| pyjson::loads(line).expect("a case is JSON"))
        .collect()
}

/// The text of a member: as written, or the bytes its `hex` member spells.
fn text(case: &Value, member: &str) -> String {
    let value = case
        .get(member)
        .unwrap_or_else(|| panic!("no member {member}"));
    if let Some(text) = value.as_str() {
        return text.to_string();
    }
    let digits = value
        .get("hex")
        .and_then(Value::as_str)
        .expect("a text or its bytes");
    let bytes: Vec<u8> = (0..digits.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&digits[at..at + 2], 16).expect("two hexadecimal digits"))
        .collect();
    String::from_utf8(bytes).expect("the bytes are UTF-8")
}

fn flag(case: &Value, member: &str) -> bool {
    match case.get(member) {
        Some(Value::Bool(flag)) => *flag,
        other => panic!("{member} is no flag: {other:?}"),
    }
}

/// The first entry of a default-format output of a header and one row.
fn first_row(output: &str) -> Entry {
    parse_default(output)
        .expect("the output is read")
        .into_iter()
        .next()
        .expect("the output holds a row")
}

/// Asks the class of each case of the file `name` about the case's field of its two rows,
/// each a default-format output of one row, with the two values the rows hold for the
/// field, and returns the row pairs asked about, the fields taken and the fields refused.
fn asked(name: &str, class: &str) -> (usize, usize, usize) {
    let side = Side::default();
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    let mut pairs = BTreeSet::new();
    let (mut taken, mut refused) = (0, 0);
    for case in cases(name) {
        assert_eq!(text(&case, "class"), class);
        let (reference, candidate) = (text(&case, "reference"), text(&case, "candidate"));
        let (e, a) = (first_row(&reference), first_row(&candidate));
        let field = text(&case, "field");
        let (ev, av) = (field_of(Some(&e), &[&field]), field_of(Some(&a), &[&field]));
        assert_eq!(ev, text(&case, "reference_value"), "{field} of {reference}");
        assert_eq!(av, text(&case, "candidate_value"), "{field} of {candidate}");
        assert_ne!(ev, av);
        let accepts = field_class_accepts(class, &e, &field, &ev, &av, Some(&a), &context)
            .expect("the class answers");
        assert_eq!(
            accepts,
            flag(&case, "accepts"),
            "{field} of {reference} against {candidate}"
        );
        pairs.insert((reference, candidate));
        if accepts {
            taken += 1;
        } else {
            refused += 1;
        }
    }
    (pairs.len(), taken, refused)
}

const CLIPPED: &str = "hgvs_clipped_reference_against_whole_alternate";
const ADDED_BASES: &str = "hgvs_added_bases_written_outside_the_transcript";

#[test]
fn a_malformed_with_reference_key_is_refused_with_its_reason() {
    assert_eq!(
        WITH_REFERENCE_SETS,
        "vep_rs_consequence_sets_with_reference_sequence"
    );
    let (mut stated, mut absent) = (0, 0);
    let mut reasons: Vec<String> = Vec::new();
    for case in cases("with_reference_sets") {
        let item = case.get("item").expect("an item");
        match (with_reference_sets(item), case.get("error")) {
            (Ok(sets), None) => {
                let expected = case.get("sets").expect("the sets");
                match sets {
                    Some(sets) => {
                        let sets: Vec<Value> = sets.into_iter().map(Value::Str).collect();
                        assert_eq!(&Value::List(sets), expected);
                        assert_eq!(item.get(WITH_REFERENCE_SETS), Some(expected));
                        stated += 1;
                    }
                    None => {
                        assert_eq!(expected, &Value::Null);
                        assert_eq!(item.get(WITH_REFERENCE_SETS), None);
                        absent += 1;
                    }
                }
            }
            (Err(Error::Comparator(refusal)), Some(expected)) => {
                assert_eq!(Some(refusal.as_str()), expected.as_str());
                let opening =
                    format!("the manifest's item on 21:200 T ENST2: {WITH_REFERENCE_SETS} ");
                reasons.push(
                    refusal
                        .strip_prefix(&opening)
                        .expect("the refusal names the item and the key")
                        .to_string(),
                );
            }
            (answer, expected) => panic!("{answer:?} where the case states {expected:?}"),
        }
    }
    assert_eq!((stated, absent), (3, 1));
    assert_eq!(
        reasons,
        [
            "is not a non-empty list",
            "is not a non-empty list",
            "is not a non-empty list",
            "holds an item that is not a sorted, comma-joined set of terms",
            "holds an item that is not a sorted, comma-joined set of terms",
            "holds an item that is not a sorted, comma-joined set of terms",
            "is not a sorted list of distinct sets",
            "is not a sorted list of distinct sets",
            "repeats vep_rs_consequence_sets and states no second outcome",
            "names Ensembl VEP's set beside another",
        ]
    );
}

#[test]
fn the_clipped_class_takes_the_offset_fields_the_reference_writes_beside_its_notation() {
    assert_eq!(
        asked("clipped_reference_offsets_beside_the_notation", CLIPPED),
        (10, 10, 9)
    );
}

#[test]
fn the_clipped_class_takes_the_offset_fields_the_candidate_alone_writes() {
    assert_eq!(
        asked("clipped_reference_offsets_on_the_candidates_side", CLIPPED),
        (11, 8, 22)
    );
}

#[test]
fn the_clipped_class_takes_the_end_offset_before_a_forward_strand_first_base() {
    assert_eq!(
        asked("clipped_reference_end_offset_before_a_first_base", CLIPPED),
        (14, 2, 12)
    );
}

#[test]
fn the_added_bases_class_takes_a_notation_and_an_offset_outside_the_transcript() {
    assert_eq!(
        asked("added_bases_outside_the_transcript", ADDED_BASES),
        (29, 19, 36)
    );
}

#[test]
fn the_added_bases_class_holds_the_offset_fields_to_its_notation() {
    assert_eq!(
        asked("added_bases_intron_offset_fields", ADDED_BASES),
        (13, 28, 15)
    );
}

#[test]
fn the_splice_term_class_takes_the_same_pairs_in_each_reading_of_the_record() {
    let class = "splice_term_past_the_record_span";
    let (_, takes) = TERM_PAIR_CLASSES
        .iter()
        .find(|(held, _)| *held == class)
        .expect("the class is a term-pair class");
    let directory = Path::new(env!("CARGO_TARGET_TMPDIR")).join("builder-cases-splice-term");
    std::fs::create_dir_all(&directory).expect("the scratch directory is made");
    let mut taken_by_reading = Vec::new();
    for case in cases("splice_term_readings") {
        let reading = text(&case, "reading");
        let mut expected = parse_default(&text(&case, "reference")).expect("the output is read");
        let mut actual = parse_default(&text(&case, "candidate")).expect("the output is read");
        // The records of rows under an identifier, where the reading gives them.
        if case.get("records") != Some(&Value::Null) {
            let records = directory.join(format!("{reading}.vcf"));
            std::fs::write(&records, text(&case, "records")).expect("the records are written");
            read_entry_references(&mut expected, &mut actual, &records)
                .expect("the records are read");
        }
        let stated = case
            .get("taken")
            .and_then(Value::as_dict)
            .expect("the answers");
        assert_eq!((expected.len(), actual.len(), stated.len()), (14, 14, 14));
        let mut taken = 0;
        for (e, a) in expected.iter().zip(&actual) {
            let location = e.location.as_deref().expect("a Location");
            let answer = takes(e, a, &AtACodingEnd::NO_CACHE).expect("the class answers");
            assert_eq!(
                Some(&Value::Bool(answer)),
                stated.get(location),
                "{reading} {location}"
            );
            taken += usize::from(answer);
        }
        taken_by_reading.push((reading, taken));
    }
    assert_eq!(
        taken_by_reading,
        [
            ("spelt".to_string(), 8),
            ("records".to_string(), 8),
            ("without".to_string(), 5)
        ]
    );
}

/// Asks `class` about each pair of HGVSp values of the file `name`, with the row's terms
/// and FLAGS, and returns the pairs asked about and the pairs taken.
fn stop_pairs(name: &str, class: fn(&str, &str, &str, &str) -> bool) -> (usize, usize) {
    let cases = cases(name);
    let mut taken = 0;
    for case in &cases {
        let (ev, av) = (text(case, "reference_value"), text(case, "candidate_value"));
        let accepts = class(&text(case, "consequence"), &text(case, "flags"), &ev, &av);
        assert_eq!(accepts, flag(case, "accepts"), "{ev} against {av}");
        taken += usize::from(accepts);
    }
    (cases.len(), taken)
}

#[test]
fn the_stop_codon_protein_comparison_class_takes_its_shapes_and_refuses_the_nearest() {
    assert_eq!(
        stop_pairs(
            "stop_codon_protein_comparison",
            |consequence, flags, ev, av| stop_codon_protein_comparison(consequence, flags, ev, av)
                .expect("the class reads the pair")
        ),
        (66, 22)
    );
}

#[test]
fn the_stop_named_on_an_incomplete_coding_sequence_is_taken_where_the_row_shows_it() {
    assert_eq!(
        stop_pairs(
            "stop_named_on_incomplete_coding_sequence",
            stop_named_on_incomplete_coding_sequence
        ),
        (19, 9)
    );
}
