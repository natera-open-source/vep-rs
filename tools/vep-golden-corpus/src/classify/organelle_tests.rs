// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the organelle module: the codon tables of a made-up cache, the shape of a
//! peptide description, and what the class asks of a question before it takes it.

use super::*;
use crate::classify::fields::FieldKey;
use crate::TempDir;

/// A cache whose one shard, of chromosome MT, holds the records of the JSON list given.
fn cache_of(tag: &str, records: &str) -> (TempDir, String) {
    let scratch =
        TempDir::create(&std::env::temp_dir(), &format!("golden_organelle_{tag}_")).unwrap();
    let shard = scratch.path.join("transcripts/MT");
    std::fs::create_dir_all(&shard).unwrap();
    std::fs::write(shard.join("1-1000000.json"), records).unwrap();
    let cache_dir = scratch.path.to_str().unwrap().to_string();
    (scratch, cache_dir)
}

/// A record whose effect cache holds the members given.
fn record(stable_id: &str, members: &str) -> String {
    format!(r#"{{"stable_id": "{stable_id}", "variation_effect_feature_cache": {{{members}}}}}"#)
}

#[test]
fn codon_tables_are_the_tables_other_than_the_standard_code() {
    // A table is a number or the digits of one as text; 0, 1, a negative number, a double,
    // a text that is no run of digits and a record without a table name no organelle
    // transcript. Nor does a code that reads every codon as the standard code does: tables
    // 11 and 24, and genetic code 26, which Ensembl VEP reads under its module's last
    // table (the 64 rows of each id in crates/vep-core/tests/data/genetic_codes.tsv). Nor
    // does a table no genetic code is held for, one past 64 bits among them.
    let records = [
        record("number", r#""codon_table": 2"#),
        record("text", r#""codon_table": "5""#),
        record("padded", r#""codon_table": "004""#),
        record("ascidian", r#""codon_table": 13"#),
        record("bacterial", r#""codon_table": 11"#),
        record("strict", r#""codon_table": "24""#),
        record("pachysolen", r#""codon_table": 26"#),
        record("no_code", r#""codon_table": 7"#),
        record("no_code_text", r#""codon_table": "25""#),
        record("large", r#""codon_table": 92233720368547758070"#),
        record("standard", r#""codon_table": 1"#),
        record("standard_text", r#""codon_table": "1""#),
        record("zero", r#""codon_table": 0"#),
        record("negative", r#""codon_table": -2"#),
        record("double", r#""codon_table": 2.0"#),
        record("signed_text", r#""codon_table": "+2""#),
        record("blank_text", r#""codon_table": " 2""#),
        record("empty_text", r#""codon_table": """#),
        record("named", r#""codon_table": "vertebrate mitochondrial""#),
        record("boolean", r#""codon_table": true"#),
        record("null", r#""codon_table": null"#),
        record("absent", r#""peptide": "M""#),
        r#"{"stable_id": "no_effect_cache"}"#.to_string(),
        r#"{"stable_id": "empty_effect_cache", "variation_effect_feature_cache": null}"#
            .to_string(),
        r#"{"stable_id": 7, "variation_effect_feature_cache": {"codon_table": 2}}"#.to_string(),
    ];
    let (_scratch, cache_dir) = cache_of("tables", &format!("[{}]", records.join(", ")));
    let expected: HashMap<String, i64> =
        [("number", 2), ("text", 5), ("padded", 4), ("ascidian", 13)]
            .into_iter()
            .map(|(stable_id, table)| (stable_id.to_string(), table))
            .collect();
    assert_eq!(organelle_codon_tables(&cache_dir).unwrap(), expected);
}

#[test]
fn a_transcript_filed_twice_keeps_its_last_organelle_table() {
    let records = [
        record("ENSTMT1", r#""codon_table": 2"#),
        record("ENSTMT1", r#""codon_table": 4"#),
        record("ENSTMT1", r#""codon_table": 1"#),
    ];
    let (_scratch, cache_dir) = cache_of("twice", &format!("[{}]", records.join(", ")));
    let tables = organelle_codon_tables(&cache_dir).unwrap();
    assert_eq!(tables.len(), 1);
    assert_eq!(tables["ENSTMT1"], 4);
}

#[test]
fn the_reader_refuses_a_cache_record_outside_what_it_holds() {
    let refusal = |records: &str| {
        let (_scratch, cache_dir) = cache_of("refused", records);
        let failure = organelle_codon_tables(&cache_dir).unwrap_err();
        assert_eq!(failure.code, 1);
        failure.line
    };
    assert_eq!(
        refusal(r#"[{"stable_id": "ENSTMT1", "variation_effect_feature_cache": 5}]"#),
        "ERROR: [vep-golden-corpus] the cache record of ENSTMT1 holds a variation_effect_feature_cache that is no object"
    );
    assert_eq!(
        refusal(&format!("[{}]", record("ENSTMT1", "\"codon_table\": \"\u{b2}\""))),
        "ERROR: [vep-golden-corpus] the cache record of ENSTMT1 holds a codon_table with a digit outside ASCII"
    );
    assert_eq!(
        refusal(
            r#"[{"stable_id": {"id": "ENSTMT1"}, "variation_effect_feature_cache": {"codon_table": 2}}]"#
        ),
        "ERROR: [vep-golden-corpus] a cache record holds a stable_id that is a list or an object"
    );
    // The same id beside the standard table is never filed, so nothing is refused.
    let standard = r#"[{"stable_id": {"id": "ENSTMT1"}, "variation_effect_feature_cache": {"codon_table": 1}}]"#;
    let (_scratch, cache_dir) = cache_of("standard", standard);
    assert!(organelle_codon_tables(&cache_dir).unwrap().is_empty());
}

#[test]
fn an_insertion_is_told_from_a_replacement_by_what_stands_before_it() {
    let shapes = |descriptions: &[&str]| -> Vec<&'static str> {
        descriptions
            .iter()
            .map(|description| peptide_description_shape(description))
            .collect()
    };
    assert_eq!(
        shapes(&[
            "Met31_Gln32insMet",
            "Ala4delinsGly",
            "Ala4delinsGlyinsSer",
            "ins",
            "delins",
            "insdup",
            ""
        ]),
        ["ins", "other", "ins", "ins", "other", "dup", "other"]
    );
}

#[test]
fn the_class_asks_for_an_hgvsp_of_one_protein_on_an_organelle_transcript() {
    let models = FieldModels {
        organelle: Models {
            codon_tables: Some([("ENSTMT1".to_string(), 2)].into_iter().collect()),
        },
        ..FieldModels::default()
    };
    let key = |feature: &str| FieldKey {
        location: "MT:3314".to_string(),
        allele: "-".to_string(),
        feature: feature.to_string(),
        feature_type: "Transcript".to_string(),
        consequences: "frameshift_variant".to_string(),
    };
    let asked = |name: &str, feature: &str, vep: &str, vep_rs: &str, models: &FieldModels| {
        let key = key(feature);
        let question = FieldQuestion {
            name,
            key: &key,
            vep,
            vep_rs,
            strand: "1",
            flags: &[],
            perl_row: None,
            rust_row: None,
            minimal_split_record: false,
            genomic_span: None,
        };
        organelle_peptide_paths_translation_table(&question, models)
    };
    let (vep, vep_rs) = ("ENSPMT1.1:p.Met3ArgfsTer12", "ENSPMT1.1:p.Met3Ter");
    let class = Verdict::Class(ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE, None);
    assert_eq!(asked("HGVSp", "ENSTMT1", vep, vep_rs, &models), class);
    // The protein is the text before the last `:p.`.
    assert_eq!(
        asked(
            "HGVSp",
            "ENSTMT1",
            "a:p.b:p.Met3Ter",
            "a:p.b:p.Met3Arg",
            &models
        ),
        class
    );
    for (name, feature, vep, vep_rs) in [
        ("HGVSc", "ENSTMT1", vep, vep_rs),
        ("HGVSp", "ENST1", vep, vep_rs),
        ("HGVSp", "ENSTMT1", vep, "ENSPMT1.2:p.Met3Ter"),
        ("HGVSp", "ENSTMT1", vep, "-"),
        ("HGVSp", "ENSTMT1", "-", vep_rs),
        ("HGVSp", "ENSTMT1", ":p.Met3ArgfsTer12", ":p.Met3Ter"),
        ("HGVSp", "ENSTMT1", "", ""),
    ] {
        assert_eq!(
            asked(name, feature, vep, vep_rs, &models),
            Verdict::Pass,
            "{name} {feature} {vep} {vep_rs}"
        );
    }
    assert_eq!(
        asked("HGVSp", "ENSTMT1", vep, vep_rs, &FieldModels::default()),
        Verdict::Pass
    );
}
