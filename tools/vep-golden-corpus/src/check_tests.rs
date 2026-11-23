// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The readings `check` makes of a manifest's records and of `variants.vcf`.

use super::*;
use crate::json::loads;
use crate::TempDir;
use std::fs;

fn manifest(text: &str) -> Object {
    match loads(text).unwrap() {
        Value::Object(manifest) => manifest,
        _ => panic!("no object"),
    }
}

#[test]
fn a_corpus_holds_a_name_in_a_dataset_on_a_chromosome() {
    let manifest = manifest(
        r#"{"records": [{"dataset": "a", "vcf": "21\t100\trs1\tA\tG"}, {"dataset": "b", "vcf": "22\t100\trs1\tA\tG"}, {"dataset": 5, "vcf": "21\t300\trs3\tA\tG"}, {"vcf": "21\t400\trs4\tA\tG"}, {"dataset": "a", "vcf": "21\t400\trs4\tA\tG"}, {"dataset": "a", "vcf": "21\t500\t.\tAT\tA"}]}"#,
    );
    let records = corpus_records(&manifest).unwrap();
    let by_name = record_name_index(records.iter().map(|record| record.vcf)).unwrap();
    let corpus = Corpus { records, by_name };
    assert_eq!(corpus.holds("a", "rs1", "21"), Ok(true));
    assert_eq!(corpus.holds("b", "rs1", "22"), Ok(true));
    assert_eq!(corpus.holds("a", "rs1", "22"), Ok(false));
    assert_eq!(corpus.holds("b", "rs1", "21"), Ok(false));
    assert_eq!(corpus.holds("a", "rs9", "21"), Ok(false));
    // A dataset that is no text names none; a record without the key is the error when a
    // name reaches it before a record that holds the name.
    assert_eq!(corpus.holds("5", "rs3", "21"), Ok(false));
    assert_eq!(corpus.holds("a", "rs4", "21"), Err(3));
    // A record without an ID is held under each name written for it.
    assert_eq!(corpus.holds("a", "21_500_AT/A", "21"), Ok(true));
    assert_eq!(corpus.holds("a", "21_501_AT/A", "21"), Ok(true));
}

#[test]
fn combinations_are_what_iterating_the_value_yields() {
    let read = |text: &str| combinations(&loads(text).unwrap(), 4);
    let names = |items: &[&str]| Ok(items.iter().map(|item| item.to_string()).collect());
    assert_eq!(
        read(r#"{"combinations": ["b", "a", "b"]}"#),
        names(&["b", "a", "b"])
    );
    assert_eq!(read(r#"{"combinations": []}"#), names(&[]));
    // Numbers, booleans and null are no names; a text gives its characters and an object
    // its keys.
    assert_eq!(
        read(r#"{"combinations": ["a", 1, null, true, 1.5]}"#),
        names(&["a"])
    );
    assert_eq!(read(r#"{"combinations": "ab"}"#), names(&["a", "b"]));
    assert_eq!(
        read(r#"{"combinations": {"x": 1, "y": 2}}"#),
        names(&["x", "y"])
    );
    for text in [
        r#"{}"#,
        r#"{"combinations": 5}"#,
        r#"{"combinations": null}"#,
        r#"{"combinations": [["a"]]}"#,
        r#"{"combinations": [{"a": 1}]}"#,
        r#"["combinations"]"#,
    ] {
        assert_eq!(
            read(text),
            Err(Failure::own(
                "record 4 of the manifest has no combinations that can be read"
            )),
            "{text}"
        );
    }
}

#[test]
fn sample_columns_are_held_to_the_manifests_samples() {
    let scratch = TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-check-").unwrap();
    let corpus = scratch.path.to_str().unwrap();
    let samples = ["S1".to_string(), "S2".to_string()];
    let fixed = FIXED_COLUMNS.join("\t");
    let write = |text: &str| fs::write(scratch.path.join("variants.vcf"), text).unwrap();
    let errors = || sample_column_errors(corpus, &samples, 2).unwrap();
    let record = "21\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\n";
    write(&format!(
        "##fileformat=VCFv4.2\n{fixed}\tFORMAT\tS1\tS2\n{record}"
    ));
    assert_eq!(errors(), Vec::<String>::new());
    write(&format!(
        "{fixed}\n{record}21\t200\trs2\tA\tG\t.\t.\t.\n#note\n"
    ));
    assert_eq!(
        errors(),
        [
            "variants.vcf header names no sample columns, the manifest keeps FORMAT S1 S2",
            "variants.vcf record 21 200 rs2 A G has 8 columns, the header 11",
            "variants.vcf header names no sample columns, the manifest keeps FORMAT S1 S2",
        ]
    );
    write(&format!("{fixed}\tFORMAT\tS2\tS1\n{record}"));
    assert_eq!(
        errors(),
        ["variants.vcf header names FORMAT S2 S1, the manifest keeps FORMAT S1 S2"]
    );
    write(record);
    assert_eq!(errors(), ["variants.vcf has no #CHROM header line"]);
    fs::remove_file(scratch.path.join("variants.vcf")).unwrap();
    let refusal = sample_column_errors(corpus, &samples, 2).unwrap_err();
    assert_eq!(
        refusal.line,
        format!(
            "ERROR: [build_golden_corpus] {corpus}/variants.vcf is not a readable text file ([Errno 2] No such file or directory: '{corpus}/variants.vcf')"
        )
    );
}

#[test]
fn the_errors_found_come_before_a_refusal_of_the_programs_own() {
    // What is written to standard error is read by the tests of the binary; here, which
    // refusal is passed on.
    let program = failure("in.vcf is not a readable text file (reason)");
    assert_eq!(after_errors(&[], program.clone()), program);
    let own = Failure::own("record 0 of the manifest has no dataset");
    assert_eq!(after_errors(&[], own.clone()), own);
}
