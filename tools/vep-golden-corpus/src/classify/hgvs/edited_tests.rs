// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the edited transcript models: the `_rna_edit` values of a made-up cache
//! in every shape the reader meets, and the two predicates of a field on an edited model.

use super::*;
use crate::TempDir;
use std::fs;

/// A cache directory under a scratch directory with the shards given, each a chromosome
/// and the text of its one shard.
fn cache(tag: &str, shards: &[(&str, &str)]) -> (TempDir, String) {
    let scratch = TempDir::create(&std::env::temp_dir(), &format!("vep-golden-corpus-{tag}-"))
        .expect("a scratch directory");
    let root = scratch.path.join("json_cache");
    for (chromosome, text) in shards {
        let directory = root.join("transcripts").join(chromosome);
        fs::create_dir_all(&directory).unwrap();
        fs::write(directory.join("1-1000000.json"), text).unwrap();
    }
    let path = root.to_str().unwrap().to_string();
    (scratch, path)
}

/// The models as sorted pairs of a stable id and its distances.
fn listed(models: &EditedModels) -> Vec<(&str, Vec<i128>)> {
    let mut listed: Vec<(&str, Vec<i128>)> = models
        .iter()
        .map(|(id, distances)| (id.as_str(), distances.iter().copied().collect()))
        .collect();
    listed.sort();
    listed
}

#[test]
fn edited_models_read_every_shape_of_an_attribute() {
    // An attribute that is no object, a code that is no text or that the attribute lacks,
    // a value that is no text, one field alone and a field that is no integer state no
    // edit; a value of four fields deletes as one of two does; trailing empty fields are
    // dropped and inner ones are not; an edit whose running sum is zero leaves its
    // transcript among the models with no distance; a later record of a stable id replaces
    // an earlier one where it states an edit; and a base outside ASCII counts as one.
    let first = r#"[
 {"stable_id": "T_mixed", "attributes": [5, "x", null, ["a"], {"code": 7, "value": "1 2"}, {"value": "1 2"}, {"code": "_rna_edit_more", "value": "5 7 TT x"}]},
 {"stable_id": "T_trailing", "attributes": [{"code": "_rna_edit", "value": "50 52  "}]},
 {"stable_id": "T_zero", "attributes": [{"code": "_rna_edit", "value": "10 9 "}]},
 {"stable_id": "T_signs", "attributes": [{"code": "_rna_edit", "value": "+5 -3 TT"}]},
 {"stable_id": "T_leading_space", "attributes": [{"code": "_rna_edit", "value": " 5 7"}]},
 {"stable_id": "T_double_space", "attributes": [{"code": "_rna_edit", "value": "5  7"}]},
 {"stable_id": "T_number_value", "attributes": [{"code": "_rna_edit", "value": 5}, {"code": "_rna_edit", "value": [1, 2]}, {"code": "_rna_edit"}]},
 {"stable_id": "T_object", "attributes": {"code": "_rna_edit", "value": "1 5"}},
 {"stable_id": "T_text", "attributes": "_rna_edit"},
 {"stable_id": "T_null", "attributes": null},
 {"stable_id": "T_absent"},
 {"stable_id": "T_sums", "attributes": [{"code": "_rna_edit", "value": "100 99 TT"}, {"code": "_rna_edit", "value": "20 24"}, {"code": "_rna_edit", "value": "50 49 TTTTT"}]},
 {"stable_id": "T_wide", "attributes": [{"code": "_rna_edit", "value": "5 6 EE"}, {"code": "_rna_edit", "value": "5 6 E"}]},
 {"stable_id": 7, "attributes": [{"code": "_rna_edit", "value": "1 4"}]},
 {"stable_id": "T_twice", "attributes": [{"code": "_rna_edit", "value": "1 4"}]},
 {"stable_id": "T_kept", "attributes": [{"code": "_rna_edit", "value": "1 2"}]},
 {"stable_id": "T_one_field", "attributes": [{"code": "_rna_edit", "value": "12"}]},
 {"stable_id": "T_empty_value", "attributes": [{"code": "_rna_edit", "value": ""}]},
 {"stable_id": "T_false", "attributes": false},
 {"stable_id": "T_zero_attributes", "attributes": 0}
]"#
    .replace('E', "\u{e9}");
    let second = r#"{"22": [{"stable_id": "T_twice", "attributes": [{"code": "_rna_edit", "value": "1 9"}]}, {"stable_id": "T_kept", "attributes": []}]}"#;
    let (_scratch, path) = cache("edited-shapes", &[("21", &first), ("22", second)]);
    let models = indel_edited_models(&path).unwrap_or_else(|failure| panic!("{}", failure.line));
    assert_eq!(
        listed(&models),
        [
            ("7", vec![4]),
            ("T_kept", vec![2]),
            ("T_mixed", vec![3]),
            ("T_signs", vec![2]),
            ("T_sums", vec![2, 5]),
            ("T_trailing", vec![3]),
            ("T_twice", vec![9]),
            ("T_wide", vec![1]),
            ("T_zero", vec![]),
        ]
    );
}

#[test]
fn edited_models_refuse_what_the_reader_cannot_hold() {
    // Attributes that are a number, which CPython cannot walk, and a coordinate 64 bits do
    // not hold end the run on a line of the tool's own.
    for (attributes, detail) in [
        ("5", "a number where its attributes stand"),
        ("true", "a number where its attributes stand"),
        ("1.5", "a number where its attributes stand"),
        (
            r#"[{"code": "_rna_edit", "value": "1 9223372036854775808"}]"#,
            "an _rna_edit coordinate 64 bits do not hold",
        ),
    ] {
        let shard = format!(r#"[{{"stable_id": "T1", "attributes": {attributes}}}]"#);
        let (_scratch, path) = cache("edited-refused", &[("21", &shard)]);
        let failure = indel_edited_models(&path).expect_err("a refusal");
        assert_eq!(
            (failure.line.as_str(), failure.code),
            (
                format!("ERROR: [vep-golden-corpus] the cache record T1 holds {detail}").as_str(),
                1
            ),
            "{attributes}"
        );
    }
}

#[test]
fn an_absent_cache_holds_no_edited_model() {
    let (scratch, _path) = cache("edited-absent", &[]);
    let absent = scratch.path.join("absent");
    let models = indel_edited_models(absent.to_str().unwrap());
    assert!(models.is_ok_and(|models| models.is_empty()));
}

#[test]
fn two_notations_an_edit_apart_differ_in_their_coordinates_alone() {
    let distances = BTreeSet::from([3, 19]);
    let apart = |vep: &str, vep_rs: &str| hgvsc_an_edit_apart(&distances, vep, vep_rs);
    assert!(apart("NM_1.1:c.247+601del", "NM_1.1:c.250+601del"));
    assert!(apart("NM_1.1:c.56+4del", "NM_1.1:c.75+4del"));
    assert!(apart("NM_1.1:c.5_9del", "NM_1.1:c.2+1_9del"));
    // Another name, another shape, no distance of the model, and a value without a name.
    assert!(!apart("NM_1.1:c.247+601del", "NM_1.2:c.250+601del"));
    assert!(!apart("NM_1.1:c.247+601del", "NM_1.1:c.250+601dup"));
    assert!(!apart("NM_1.1:c.247+601del", "NM_1.1:c.248+601del"));
    assert!(!apart("c.247+601del", "NM_1.1:c.250+601del"));
    assert!(!apart("NM_1.1:c.247+601del", "-"));
    // An integer 64 bits do not hold is not read.
    assert!(!apart(
        "NM_1.1:c.9223372036854775808del",
        "NM_1.1:c.9223372036854775811del"
    ));
}

#[test]
fn an_offset_left_empty_is_a_non_zero_integer_against_none() {
    assert!(offset_left_empty("", "21"));
    assert!(offset_left_empty("-", "-7"));
    assert!(!offset_left_empty("3", "21"));
    assert!(!offset_left_empty("-", "0"));
    assert!(!offset_left_empty("-", "dup"));
    assert!(!offset_left_empty("-", "-"));
}
