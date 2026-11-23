// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the LOFTEE field class of a GERP walk without the conservation database:
//! the rows it takes, with and without `LoF_info`, and the nearest rows it passes.

use super::super::fields::{
    classify_field_divergence, FieldKey, FieldModels, FieldQuestion, FieldRow, Verdict,
};
use super::*;

const ITEMS: &str =
    "GERP_DIST:705,BP_DIST:702,PERCENTILE:0.155875299760192,DIST_FROM_LAST_EXON:266,50_BP_RULE:PASS";

fn key(consequences: &str) -> FieldKey {
    FieldKey {
        location: "21:37442543".to_string(),
        allele: "T".to_string(),
        feature: "ENST00000290349".to_string(),
        feature_type: "Transcript".to_string(),
        consequences: consequences.to_string(),
    }
}

/// A row read with the plugin's four fields: the ones `written` names with their values,
/// `-` for the others.
fn row(written: &[(&str, &str)]) -> FieldRow {
    let mut row: FieldRow = LOFTEE_FIELDS
        .iter()
        .map(|name| (name.to_string(), "-".to_string()))
        .collect();
    for (name, value) in written {
        row.insert(name.to_string(), value.to_string());
    }
    row
}

/// The verdict on the field `name` of a key of the terms `consequences`, each row's value
/// of the field read from the row.
fn asked(name: &str, consequences: &str, perl_row: &FieldRow, rust_row: &FieldRow) -> Verdict {
    let key = key(consequences);
    let value = |row: &FieldRow| row.get(name).cloned().unwrap_or_else(|| "-".to_string());
    let (vep, vep_rs) = (value(perl_row), value(rust_row));
    let question = FieldQuestion {
        name,
        key: &key,
        vep: &vep,
        vep_rs: &vep_rs,
        strand: "1",
        flags: &[],
        perl_row: Some(perl_row),
        rust_row: Some(rust_row),
        minimal_split_record: false,
        genomic_span: None,
    };
    let verdict = gerp_walk_without_conservation_database(&question, &FieldModels::default());
    // The dispatch names the class where the function takes the field, and no other class
    // takes a field of the plugin on such a row.
    let named = classify_field_divergence(&question, &FieldModels::default());
    match &verdict {
        Verdict::Class(class, second) => assert_eq!(named, Some((*class, *second))),
        _ => assert_eq!(named, None),
    }
    verdict
}

const TAKEN: Verdict = Verdict::Class(LOFTEE_GERP_WALK_WITHOUT_CONSERVATION_DATABASE, None);

#[test]
fn a_field_vep_rs_alone_writes_on_a_row_of_a_finished_walk_is_the_class() {
    let none = row(&[]);
    let called = row(&[("LoF", "HC"), ("LoF_info", ITEMS)]);
    for name in ["LoF", "LoF_info"] {
        assert_eq!(asked(name, "stop_gained", &none, &called), TAKEN, "{name}");
    }
    let filtered = row(&[
        ("LoF", "LC"),
        ("LoF_filter", "END_TRUNC"),
        ("LoF_info", ITEMS),
    ]);
    for name in ["LoF", "LoF_filter", "LoF_info"] {
        assert_eq!(
            asked(name, "frameshift_variant", &none, &filtered),
            TAKEN,
            "{name}"
        );
    }
    let three_items = "GERP_DIST:1498,BP_DIST:1496,PERCENTILE:0.240412979351032";
    let flagged = row(&[
        ("LoF", "HC"),
        ("LoF_flags", "NO_EXON_NUMBER"),
        ("LoF_info", three_items),
    ]);
    let terms = "frameshift_variant,splice_region_variant";
    for name in ["LoF", "LoF_flags", "LoF_info"] {
        assert_eq!(asked(name, terms, &none, &flagged), TAKEN, "{name}");
    }
    // An empty text is no value either.
    let empty = row(&[("LoF", ""), ("LoF_info", "")]);
    assert_eq!(asked("LoF", "stop_gained", &empty, &called), TAKEN);
    // The total is an item wherever it stands in the value, at either separator.
    for items in ["BP_DIST:702,GERP_DIST:705", "BP_DIST:702&GERP_DIST:705"] {
        let later = row(&[("LoF", "HC"), ("LoF_info", items)]);
        assert_eq!(
            asked("LoF_info", "stop_gained", &none, &later),
            TAKEN,
            "{items}"
        );
    }
}

#[test]
fn a_call_written_without_its_items_under_skip_lof_info_is_the_class() {
    let none = row(&[]);
    let called = row(&[("LoF", "HC")]);
    assert_eq!(asked("LoF", "stop_gained", &none, &called), TAKEN);
    let filtered = row(&[("LoF", "LC"), ("LoF_filter", "END_TRUNC")]);
    for name in ["LoF", "LoF_filter"] {
        assert_eq!(
            asked(name, "stop_gained", &none, &filtered),
            TAKEN,
            "{name}"
        );
    }
    let flagged = row(&[("LoF", "HC"), ("LoF_flags", "NO_EXON_NUMBER")]);
    assert_eq!(
        asked("LoF_flags", "frameshift_variant", &none, &flagged),
        TAKEN
    );
}

#[test]
fn the_nearest_rows_pass() {
    let none = row(&[]);
    let called = row(&[("LoF", "HC"), ("LoF_info", ITEMS)]);
    // A row on which Ensembl VEP wrote one of the four.
    for held in [
        ("LoF", "HC"),
        ("LoF_filter", "END_TRUNC"),
        ("LoF_flags", "SINGLE_EXON"),
        ("LoF_info", "PERCENTILE:0.155875299760192"),
    ] {
        for name in LOFTEE_FIELDS {
            assert_eq!(
                asked(name, "stop_gained", &row(&[held]), &called),
                Verdict::Pass,
                "{} {name}",
                held.0
            );
        }
    }
    // A key that holds neither term.
    for terms in [
        "splice_donor_variant",
        "missense_variant",
        "splice_acceptor_variant,intron_variant",
        "stop_gained_variant",
        "",
    ] {
        assert_eq!(
            asked("LoF", terms, &none, &called),
            Verdict::Pass,
            "{terms}"
        );
    }
    // A row on which vep-rs wrote no call, or items without the total.
    let items_alone = row(&[("LoF_info", ITEMS)]);
    assert_eq!(
        asked("LoF_info", "stop_gained", &none, &items_alone),
        Verdict::Pass
    );
    let filter_alone = row(&[("LoF_filter", "END_TRUNC")]);
    assert_eq!(
        asked("LoF_filter", "stop_gained", &none, &filter_alone),
        Verdict::Pass
    );
    let no_total = row(&[("LoF", "HC"), ("LoF_info", "BP_DIST:702,XGERP_DIST:705")]);
    for name in ["LoF", "LoF_info"] {
        assert_eq!(
            asked(name, "stop_gained", &none, &no_total),
            Verdict::Pass,
            "{name}"
        );
    }
    let longer_key = row(&[("LoF", "HC"), ("LoF_info", "GERP_DISTANCE:705")]);
    assert_eq!(
        asked("LoF_info", "stop_gained", &none, &longer_key),
        Verdict::Pass
    );
    // A field of the four that vep-rs leaves unwritten too, the fields Ensembl VEP alone
    // writes, and a field outside the four.
    assert_eq!(
        asked("LoF_filter", "stop_gained", &none, &called),
        Verdict::Pass
    );
    assert_eq!(asked("LoF", "stop_gained", &called, &none), Verdict::Pass);
    let mut with_other = called.clone();
    with_other.insert("NMD".to_string(), "NMD_escaping_variant".to_string());
    assert_eq!(
        asked("NMD", "stop_gained", &none, &with_other),
        Verdict::Pass
    );
}

#[test]
fn a_question_without_both_rows_passes() {
    let key = key("stop_gained");
    let called = row(&[("LoF", "HC"), ("LoF_info", ITEMS)]);
    let none = row(&[]);
    let mut question = FieldQuestion {
        name: "LoF",
        key: &key,
        vep: "-",
        vep_rs: "HC",
        strand: "1",
        flags: &[],
        perl_row: None,
        rust_row: Some(&called),
        minimal_split_record: false,
        genomic_span: None,
    };
    let models = FieldModels::default();
    assert_eq!(
        gerp_walk_without_conservation_database(&question, &models),
        Verdict::Pass
    );
    question.perl_row = Some(&none);
    question.rust_row = None;
    assert_eq!(
        gerp_walk_without_conservation_database(&question, &models),
        Verdict::Pass
    );
    question.rust_row = Some(&called);
    assert_eq!(
        gerp_walk_without_conservation_database(&question, &models),
        TAKEN
    );
    // A value the caller states for Ensembl VEP is no field left unwritten.
    question.vep = "LC";
    assert_eq!(
        gerp_walk_without_conservation_database(&question, &models),
        Verdict::Pass
    );
}
