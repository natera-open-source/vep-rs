// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the allele-number class of a record that repeats an alternate or lists
//! one equal to the reference allele: the rows it takes under either name of the field,
//! and the nearest rows it passes.

use super::super::fields::{
    classify_field_divergence, FieldKey, FieldModels, FieldQuestion, FieldRow, Verdict,
};
use super::*;

fn key(allele: &str) -> FieldKey {
    FieldKey {
        location: "21:30289183".to_string(),
        allele: allele.to_string(),
        feature: "ENST00000416044".to_string(),
        feature_type: "Transcript".to_string(),
        consequences: "stop_retained_variant".to_string(),
    }
}

/// A row read with the field of the written allele list: the list given, or `-` for a row
/// that carries none.
fn row(list: Option<&str>) -> FieldRow {
    let list = list.unwrap_or("-").to_string();
    FieldRow::from([(WRITTEN_ALLELE_LIST_FIELD.to_string(), list)])
}

/// The verdict on the field `name` of a key of `allele` whose two values are `vep` and
/// `vep_rs`, beside the two rows given.
fn asked(
    name: &str,
    allele: &str,
    values: (&str, &str),
    perl_row: Option<&FieldRow>,
    rust_row: Option<&FieldRow>,
) -> Verdict {
    let key = key(allele);
    let question = FieldQuestion {
        name,
        key: &key,
        vep: values.0,
        vep_rs: values.1,
        strand: "-1",
        flags: &[],
        perl_row,
        rust_row,
        minimal_split_record: false,
        genomic_span: None,
    };
    let verdict = allele_number_repeated_or_reference_equal_alt(&question, &FieldModels::default());
    // The dispatch names the class where the function takes the field, and no other class
    // takes an allele number.
    let named = classify_field_divergence(&question, &FieldModels::default());
    match &verdict {
        Verdict::Class(class, second) => assert_eq!(named, Some((*class, *second))),
        _ => assert_eq!(named, None),
    }
    verdict
}

/// `asked` of `ALLELE_NUM` beside two rows that carry `list`.
fn numbered(list: &str, allele: &str, values: (&str, &str)) -> Verdict {
    let row = row(Some(list));
    asked("ALLELE_NUM", allele, values, Some(&row), Some(&row))
}

const TAKEN: Verdict = Verdict::Class(ALLELE_NUMBER_REPEATED_OR_REFERENCE_EQUAL_ALT, None);

#[test]
fn a_number_after_a_dropped_entry_is_the_class_under_either_name_of_the_field() {
    assert_eq!(numbered("C/G/G/T", "T", ("2", "3")), TAKEN);
    assert_eq!(numbered("A/A/C", "C", ("1", "2")), TAKEN);
    assert_eq!(numbered("C/*/G/G/T", "T", ("3", "4")), TAKEN);
    // The JSON format's name of the field, as a run reads it of that format alone.
    let row = row(Some("C/G/G/T"));
    assert_eq!(
        asked("allele_num", "T", ("2", "3"), Some(&row), Some(&row)),
        TAKEN
    );
    assert_eq!(WRITTEN_ALLELE_LIST_FIELD, "UPLOADED_ALLELE");
}

#[test]
fn a_row_of_an_indel_record_is_the_class_by_the_allele_the_parser_spells() {
    // The key of a row of `CA/C/C/CAT` carries `-` or `AT`, the alleles without the first
    // base they share with the reference allele.
    assert_eq!(numbered("CA/C/C/CAT", "AT", ("2", "3")), TAKEN);
    assert_eq!(numbered("CA/C/C/CAT", "-", ("1", "2")), Verdict::Pass);
    // Alternates that share no first base are carried as written.
    assert_eq!(numbered("CA/C/C/GAT", "GAT", ("2", "3")), TAKEN);
    assert_eq!(numbered("CA/C/C/GAT", "AT", ("2", "3")), Verdict::Pass);
    // An allele that is neither spelling of an alternate.
    assert_eq!(numbered("CA/C/C/CAT", "T", ("2", "3")), Verdict::Pass);
}

#[test]
fn a_row_without_the_list_passes() {
    let listed = row(Some("C/G/G/T"));
    let none = row(None);
    let unread = FieldRow::new();
    let pair = ("2", "3");
    for (perl_row, rust_row) in [
        (&none, &none),
        (&listed, &none),
        (&none, &listed),
        (&unread, &unread),
        (&listed, &unread),
    ] {
        assert_eq!(
            asked("ALLELE_NUM", "T", pair, Some(perl_row), Some(rust_row)),
            Verdict::Pass
        );
    }
    // A question without both rows, and two rows whose lists differ.
    assert_eq!(
        asked("ALLELE_NUM", "T", pair, Some(&listed), None),
        Verdict::Pass
    );
    assert_eq!(
        asked("ALLELE_NUM", "T", pair, None, Some(&listed)),
        Verdict::Pass
    );
    let other = row(Some("C/G/T/T"));
    assert_eq!(
        asked("ALLELE_NUM", "T", pair, Some(&listed), Some(&other)),
        Verdict::Pass
    );
}

#[test]
fn numbers_that_do_not_fit_the_list_and_other_fields_pass() {
    // The candidate numbered 3 where the allele's place is 2, a reference number that is
    // not the count, and the two sides exchanged.
    assert_eq!(numbered("A/A/C", "C", ("1", "3")), Verdict::Pass);
    assert_eq!(numbered("C/G/G/T", "T", ("1", "3")), Verdict::Pass);
    assert_eq!(numbered("C/G/G/T", "T", ("3", "2")), Verdict::Pass);
    // The repeated alternate itself, a list that drops no entry and an allele the list
    // does not hold.
    assert_eq!(numbered("C/G/G/T", "G", ("1", "2")), Verdict::Pass);
    assert_eq!(numbered("C/G/T", "T", ("2", "3")), Verdict::Pass);
    assert_eq!(numbered("C/G/G/T", "A", ("2", "3")), Verdict::Pass);
    // Another field of the row with the same two values, and the field in another case.
    let row = row(Some("C/G/G/T"));
    for name in ["UPLOADED_ALLELE", "STRAND", "Allele_Num", "ALLELE_NUMBER"] {
        assert_eq!(
            asked(name, "T", ("2", "3"), Some(&row), Some(&row)),
            Verdict::Pass,
            "{name}"
        );
    }
}
