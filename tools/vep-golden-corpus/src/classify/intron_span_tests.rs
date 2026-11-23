// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of the intron-span module: the coding bases a `CDS_position` cell spans, and
//! what the class asks of a key's two sets and two rows before it takes the key.

use super::*;

const CLASS: Option<&str> = Some(REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE);

/// Ensembl VEP's set on the documented row and vep-rs's, as `classify` holds a set: the
/// terms sorted.
const PERL_SET: &str = "intron_variant,protein_altering_variant,splice_acceptor_variant,\
splice_donor_5th_base_variant,splice_donor_variant,stop_lost";
const RUST_SET: &str = "inframe_deletion,intron_variant,splice_acceptor_variant,\
splice_donor_5th_base_variant,splice_donor_variant";

/// A reference allele of `bases` bases.
fn allele(bases: usize) -> String {
    "ACGT".chars().cycle().take(bases).collect()
}

/// A row with the three cells the class reads.
fn row(bam_edit: &str, cds_position: &str, used_ref: &str) -> FieldRow {
    [
        ("BAM_EDIT", bam_edit),
        ("CDS_position", cds_position),
        ("USED_REF", used_ref),
    ]
    .into_iter()
    .map(|(name, value)| (name.to_string(), value.to_string()))
    .collect()
}

/// A row of an edited model whose 12 coding bases stand under a reference allele of 30.
fn documented_row() -> FieldRow {
    row("OK", "10-21", &allele(30))
}

/// What the class answers for the documented sets on the two rows given.
fn asked(reference: &FieldRow, candidate: &FieldRow) -> Option<&'static str> {
    classify_intron_span_pair(
        "Transcript",
        "-",
        PERL_SET,
        RUST_SET,
        Some(reference),
        Some(candidate),
    )
}

#[test]
fn coding_bases_reads_a_range_of_two_numbers() {
    assert_eq!(coding_bases("210-569"), Some(360));
    assert_eq!(coding_bases("7-8"), Some(2));
    for cell in [
        "210", "7-7", "-", "?-326", "210-?", "569-210", "-210", "210-", "a-b",
    ] {
        assert_eq!(coding_bases(cell), None, "{cell}");
    }
}

#[test]
fn the_documented_row_takes_the_class() {
    let row = documented_row();
    assert_eq!(asked(&row, &row), CLASS);
}

#[test]
fn a_key_of_another_feature_type_or_without_a_row_takes_no_class() {
    let row = documented_row();
    let ask = |feature_type, reference, candidate| {
        classify_intron_span_pair(feature_type, "-", PERL_SET, RUST_SET, reference, candidate)
    };
    assert_eq!(ask("Transcript", Some(&row), Some(&row)), CLASS);
    assert_eq!(ask("RegulatoryFeature", Some(&row), Some(&row)), None);
    assert_eq!(ask("-", Some(&row), Some(&row)), None);
    assert_eq!(ask("Transcript", None, Some(&row)), None);
    assert_eq!(ask("Transcript", Some(&row), None), None);
}

#[test]
fn a_key_of_an_allele_that_is_no_run_of_bases_takes_no_class() {
    let row = documented_row();
    let ask = |allele| {
        classify_intron_span_pair(
            "Transcript",
            allele,
            PERL_SET,
            RUST_SET,
            Some(&row),
            Some(&row),
        )
    };
    // A deletion's allele and the bases a record puts in the span's place take the class.
    for allele in ["-", "ACG"] {
        assert_eq!(ask(allele), CLASS, "{allele}");
    }
    for allele in ["deletion", "<DEL>", "CNV", ""] {
        assert_eq!(ask(allele), None, "{allele}");
    }
}

#[test]
fn a_row_of_a_model_that_is_not_edited_takes_no_class() {
    let edited = documented_row();
    for status in ["-", "FAILED", ""] {
        let other = row(status, "10-21", &allele(30));
        assert_eq!(asked(&other, &edited), None, "{status}");
        assert_eq!(asked(&edited, &other), None, "{status}");
        assert_eq!(asked(&other, &other), None, "{status}");
    }
}

#[test]
fn a_span_with_an_end_off_the_coding_sequence_takes_no_class() {
    for cell in ["?-21", "10-?", "10", "-"] {
        let other = row("OK", cell, &allele(30));
        assert_eq!(asked(&other, &other), None, "{cell}");
    }
}

#[test]
fn a_reference_allele_no_longer_than_the_coding_span_takes_no_class() {
    // Twelve coding bases: an allele of twelve bases or fewer holds no base the coding
    // span lacks, and a row without the field names no allele.
    for used_ref in [allele(12), allele(11), "-".to_string(), String::new()] {
        let other = row("OK", "10-21", &used_ref);
        assert_eq!(asked(&other, &other), None, "{} bases", used_ref.len());
    }
    let longer = row("OK", "10-21", &allele(13));
    assert_eq!(asked(&longer, &longer), CLASS);
}

#[test]
fn rows_that_differ_in_the_span_or_the_allele_take_no_class() {
    let reference = documented_row();
    assert_eq!(asked(&reference, &row("OK", "10-22", &allele(30))), None);
    assert_eq!(asked(&reference, &row("OK", "10-21", &allele(31))), None);
}

#[test]
fn sets_without_the_terms_of_a_whole_intron_take_no_class() {
    let row = documented_row();
    for missing in [
        "splice_acceptor_variant",
        "splice_donor_variant",
        "intron_variant",
    ] {
        let without = |set: &str| {
            let kept: Vec<&str> = set.split(',').filter(|term| *term != missing).collect();
            kept.join(",")
        };
        let answer = classify_intron_span_pair(
            "Transcript",
            "-",
            &without(PERL_SET),
            &without(RUST_SET),
            Some(&row),
            Some(&row),
        );
        assert_eq!(answer, None, "{missing}");
    }
}

#[test]
fn another_pair_of_coding_terms_takes_no_class() {
    let row = documented_row();
    let shared = "intron_variant,splice_acceptor_variant,splice_donor_variant";
    let ask = |perl_terms: &str, rust_terms: &str| {
        classify_intron_span_pair(
            "Transcript",
            "-",
            &format!("{shared},{perl_terms}"),
            &format!("{shared},{rust_terms}"),
            Some(&row),
            Some(&row),
        )
    };
    assert_eq!(
        ask("protein_altering_variant,stop_lost", "inframe_deletion"),
        CLASS
    );
    // The two sides exchanged, one of Ensembl VEP's two terms alone, a further term on
    // either side, and another term on vep-rs's side.
    assert_eq!(
        ask("inframe_deletion", "protein_altering_variant,stop_lost"),
        None
    );
    assert_eq!(ask("stop_lost", "inframe_deletion"), None);
    assert_eq!(ask("protein_altering_variant", "inframe_deletion"), None);
    assert_eq!(
        ask(
            "protein_altering_variant,stop_lost,stop_gained",
            "inframe_deletion"
        ),
        None
    );
    assert_eq!(
        ask(
            "protein_altering_variant,stop_lost",
            "inframe_deletion,stop_retained_variant"
        ),
        None
    );
    assert_eq!(
        ask("protein_altering_variant,stop_lost", "frameshift_variant"),
        None
    );
}
