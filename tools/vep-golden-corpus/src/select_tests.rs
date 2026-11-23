// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The scan of the inputs, the named records of a command line and the order of a corpus's
//! records.

use super::*;
use crate::json::loads;
use crate::TempDir;
use std::fs;

fn scratch() -> TempDir {
    TempDir::create(&std::env::temp_dir(), "vep-golden-corpus-select-").unwrap()
}

fn file(directory: &TempDir, name: &str, text: &str) -> String {
    let path = directory.path.join(name);
    fs::write(&path, text).unwrap();
    path.to_str().unwrap().to_string()
}

const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn wanted(sid: &str, names: &[&str]) -> HashMap<String, HashSet<String>> {
    let names = names.iter().map(|name| name.to_string()).collect();
    HashMap::from([(sid.to_string(), names)])
}

fn args(include: &[&str], keep_records: bool) -> SelectArgs {
    SelectArgs {
        release: "116".to_string(),
        assembly: "GRCh38".to_string(),
        species: None,
        dataset: Vec::new(),
        out: "corpus".to_string(),
        k: 3,
        seed: 1,
        max_records: 1500,
        multi_allelic_per_dataset: 4,
        max_span: 1_000_000,
        focus_term: Vec::new(),
        focus_k: None,
        include: include.iter().map(|arg| arg.to_string()).collect(),
        keep_records,
        cache: None,
        flank: 5000,
        covers: Vec::new(),
        family_columns: None,
        keep_samples: 0,
        genotype_exemplars: 8,
        require_genotype_class: Vec::new(),
        record_classes: None,
        record_class_exemplars: 2,
        matrix: None,
        drop_reference_statements: false,
        flag: Vec::new(),
        fasta_name: None,
    }
}

fn manifest(text: &str) -> Object {
    match loads(text).unwrap() {
        Value::Object(manifest) => manifest,
        _ => panic!("no object"),
    }
}

#[test]
fn resolve_records_finds_a_record_under_every_name_it_can_carry() {
    // A wanted Uploaded_variation is resolved to its input record by the names
    // `record_names` lists, a record that has an ID by the names written for one without.
    let scratch = scratch();
    let input = file(
        &scratch,
        "in.vcf",
        &format!(
            "{HEADER}21\t30319403\t0\ttca\tt\t.\t.\t.\n21\t30289183\trs1;rs2\tC\tG\t.\t.\t.\n21\t30289184\trs3\tC\tG\t.\t.\t.\n"
        ),
    );
    let names = ["21_30319404_ca/-", "rs1", "21_30289184_C/G"];
    let inputs = [("s".to_string(), input)];
    let resolved = resolve_records(&inputs, &wanted("s", &names), 0, 1, 0, 0).unwrap();
    let found = &resolved.found["s"];
    assert_eq!(found.len(), 3);
    assert_eq!(
        found.get("21", "21_30319404_ca/-"),
        Some(&["21\t30319403\t0\ttca\tt\t.\t.\t.".to_string()][..])
    );
    assert!(found.get("21", "rs1").is_some() && found.get("21", "21_30289184_C/G").is_some());
    assert_eq!(
        resolved.header,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    );
    assert!(resolved.multi.is_empty() && resolved.samples.is_empty());
}

#[test]
fn resolve_records_samples_the_multi_allelic_records_of_each_stratum() {
    let scratch = scratch();
    let mut text = HEADER.to_string();
    for place in 0..40 {
        let position = 100 + place * 10;
        text += &format!("21\t{position}\tsnv{place}\tA\tG,T\t.\t.\tINFO{place}\textra\n");
        text += &format!("21\t{}\tindel{place}\tCA\tC,CGG\t.\t.\t.\n", position + 1);
        text += &format!("21\t{}\tsym{place}\tA\tG,<DEL>\t.\t.\t.\n", position + 2);
        text += &format!("21\t{}\tstar{place}\tA\tG,*\t.\t.\t.\n", position + 3);
        text += &format!("21\t{}\tone{place}\tA\tG\t.\t.\t.\n", position + 4);
    }
    let input = file(&scratch, "in.vcf", &text);
    let inputs = [("a".to_string(), input.clone()), ("b".to_string(), input)];
    let resolved = resolve_records(&inputs, &HashMap::new(), 2, 1, 0, 0).unwrap();
    let ids = |sid: &str, stratum: &str| -> Vec<String> {
        let lines = &resolved.multi[sid][stratum];
        lines
            .iter()
            .map(|line| line.split('\t').nth(2).unwrap().to_string())
            .collect()
    };
    // The records `random.Random(1)` leaves in each sample of two, the two datasets drawn
    // one after the other from the one generator.
    assert_eq!(ids("a", "multi_allelic_snv"), ["snv16", "snv1"]);
    assert_eq!(ids("a", "multi_allelic_indel"), ["indel26", "indel14"]);
    assert_eq!(ids("b", "multi_allelic_snv"), ["snv10", "snv34"]);
    assert_eq!(ids("b", "multi_allelic_indel"), ["indel32", "indel39"]);
    // A record line is its first eight columns.
    assert!(resolved.multi["a"]["multi_allelic_snv"][0].ends_with("\t.\t.\tINFO16"));
    assert!(resolved.found.is_empty());
    // A count below zero draws and keeps nothing; zero reads no input at all.
    let none = resolve_records(&inputs, &HashMap::new(), -1, 1, 0, 0).unwrap();
    assert!(none
        .multi
        .values()
        .all(|strata| strata.values().all(Vec::is_empty)));
    let absent = [("a".to_string(), "/nonexistent/in.vcf".to_string())];
    let unread = resolve_records(&absent, &HashMap::new(), 0, 1, 0, 0).unwrap();
    assert_eq!(unread.header, FIXED_COLUMNS.join("\t"));
    assert!(resolve_records(&absent, &wanted("a", &["rs1"]), 0, 1, 0, 0).is_err());
}

#[test]
fn named_records_are_the_included_and_the_kept_ones_of_the_datasets_given() {
    let existing = manifest(
        r#"{"records": [{"dataset": "s01", "vcf": "21\t100\t.\tAT\tA\t.\t.\t.", "combinations": []}, {"dataset": "s01", "vcf": "21\t200\trs2;x\tC\tT\t.\t.\t.", "combinations": []}], "included_records": {"s01": ["rs7", "rs8"]}}"#,
    );
    let datasets = HashSet::from(["s01"]);
    let record = |dataset: &str, label: &str, names: &[&str]| NamedRecord {
        dataset: dataset.to_string(),
        label: label.to_string(),
        names: names.iter().map(|name| name.to_string()).collect(),
    };
    assert_eq!(
        named_records(&args(&["s01=rs1", "s01=a=b"], false), &existing, &datasets),
        Ok(vec![
            record("s01", "--include s01=rs1", &["rs1"]),
            record("s01", "--include s01=a=b", &["a=b"]),
        ])
    );
    assert_eq!(
        named_records(&args(&["s01=rs1"], true), &existing, &datasets),
        Ok(vec![
            record("s01", "--include s01=rs1", &["rs1"]),
            record(
                "s01",
                "kept record 21 100 . AT A",
                &["21_100_AT/A", "21_101_AT/A"]
            ),
            record("s01", "kept record 21 200 rs2;x C T", &["rs2"]),
            record("s01", "carried include s01=rs7", &["rs7"]),
            record("s01", "carried include s01=rs8", &["rs8"]),
        ])
    );
    let refused = |include: &[&str], keep: bool, existing: &Object, datasets: &[&str]| {
        let datasets: HashSet<&str> = datasets.iter().copied().collect();
        named_records(&args(include, keep), existing, &datasets)
            .unwrap_err()
            .line
    };
    let program = |message: &str| format!("ERROR: [build_golden_corpus] {message}");
    for arg in ["rs1", "=rs1", "s01=", "="] {
        assert_eq!(
            refused(&[arg], false, &existing, &["s01"]),
            program(&format!("--include wants sid=name, got '{arg}'"))
        );
    }
    assert_eq!(
        refused(&["s01=rs1", "s02=rs1"], false, &existing, &["s01"]),
        program("--include s02=rs1 names dataset 's02', which is not among --dataset")
    );
    assert_eq!(
        refused(&[], true, &Object::new(), &["s01"]),
        program("--keep-records wants a manifest to regenerate in corpus")
    );
    assert_eq!(
        refused(&["s02=x"], true, &existing, &["s02"]),
        program("kept record 21 100 . AT A names dataset 's01', which is not among --dataset")
    );
    // A dataset that is no text is among no `--dataset`, and is refused in its turn.
    let odd = manifest(
        r#"{"records": [{"dataset": "s01", "vcf": "21\t100\trs1\tA\tG"}, {"dataset": 7, "vcf": "21\t200\trs2\tC\tT"}, {"dataset": "s09", "vcf": "21\t300\trs3\tC\tT"}]}"#,
    );
    assert_eq!(
        refused(&[], true, &odd, &["s01"]),
        program("kept record 21 200 rs2 C T names dataset 7, which is not among --dataset")
    );
    assert_eq!(
        refused(&["s05=x"], true, &odd, &["s01"]),
        program("--include s05=x names dataset 's05', which is not among --dataset")
    );
    // Records the tool cannot read as records are refused in its own words.
    let own = |text: &str| refused(&[], true, &manifest(text), &["s01"]);
    for text in [
        r#"{"records": {"dataset": "s01"}}"#,
        r#"{"records": [{"dataset": "s01"}]}"#,
        r#"{"records": [{"vcf": "21\t100\trs1\tA\tG"}]}"#,
        r#"{"records": [{"dataset": ["s01"], "vcf": "21\t100\trs1\tA\tG"}]}"#,
        r#"{"records": [{"dataset": "s01", "vcf": "21\t100"}]}"#,
    ] {
        assert!(
            own(text).starts_with("ERROR: [vep-golden-corpus] "),
            "{text}"
        );
    }
    assert_eq!(
        named_records(
            &args(&[], true),
            &manifest(r#"{"records": [], "covers": ["x"]}"#),
            &datasets
        ),
        Ok(Vec::new())
    );
}

#[test]
fn included_records_are_the_sorted_names_of_each_dataset() {
    let carried = manifest(r#"{"s02": ["rs9", "rs1"], "s01": ["rs5"]}"#);
    let includes = [
        "s01=rs2".to_string(),
        "s01=rs10".to_string(),
        "s01=rs2".to_string(),
        "s00=b=c".to_string(),
    ];
    let text = |object: Object| crate::json::dumps_object(&object, &crate::json::Style::PLAIN);
    assert_eq!(
        text(included_by_dataset(&includes, Some(&carried))),
        r#"{"s00": ["b=c"], "s01": ["rs10", "rs2", "rs5"], "s02": ["rs1", "rs9"]}"#
    );
    assert_eq!(
        text(included_by_dataset(&includes[..1], None)),
        r#"{"s01": ["rs2"]}"#
    );
    assert!(included_by_dataset(&[], None).is_empty());
}

#[test]
fn records_are_ordered_by_chromosome_position_and_alleles() {
    let lines = [
        "MT\t5\t.\tA\tG",
        "GL000008.2\t1\t.\tA\tG",
        "Y\t7\t.\tA\tG",
        "X\t7\t.\tA\tG",
        "10\t99\t.\tA\tG",
        "9\t100\t.\tA\tT",
        "9\t100\t.\tA\tG",
        "9\t100\t.\tAA\tC",
        "09\t100\t.\tA\tA",
        "9\t20\t.\tT\tC",
        "chr1\t1\t.\tA\tG",
        "1000000000000000000000\t1\t.\tA\tG",
        "\t1\t.\tA\tG",
    ];
    let mut ordered: Vec<&str> = lines.to_vec();
    ordered.sort_by_key(|line| record_order(line).unwrap());
    assert_eq!(
        ordered,
        [
            "9\t20\t.\tT\tC",
            "09\t100\t.\tA\tA",
            "9\t100\t.\tA\tG",
            "9\t100\t.\tA\tT",
            "9\t100\t.\tAA\tC",
            "10\t99\t.\tA\tG",
            "1000000000000000000000\t1\t.\tA\tG",
            "X\t7\t.\tA\tG",
            "Y\t7\t.\tA\tG",
            "MT\t5\t.\tA\tG",
            // The other names keep the order they came in.
            "GL000008.2\t1\t.\tA\tG",
            "chr1\t1\t.\tA\tG",
            "\t1\t.\tA\tG",
        ]
    );
    assert!(record_order("9\tx\t.\tA\tG").is_err());
}

#[test]
fn a_file_name_is_the_last_part_of_a_path() {
    let names = [
        ("a/b.txt", "b.txt"),
        ("b.txt", "b.txt"),
        ("/abs/dir/s02.vcf.gz", "s02.vcf.gz"),
        (".", ""),
        ("/", ""),
        ("..", ".."),
        ("a/..", ".."),
    ];
    for (path, name) in names {
        assert_eq!(file_name(path), name, "{path}");
    }
    assert_eq!(chromosome("21:100-200"), "21");
    assert_eq!(chromosome("HLA:A:5"), "HLA");
    assert_eq!(chromosome("21"), "21");
}
