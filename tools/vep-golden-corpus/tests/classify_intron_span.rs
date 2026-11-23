// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Runs `classify` on made-up outputs whose two engines differ on a key as
//! `refseq_edited_model_intron_span_coding_reference` documents: a deletion with a whole
//! intron of an edited RefSeq model between its two coding ends. The row of the shape is
//! written under the class, and the nearest rows (a model that is not edited, a reference
//! allele of the coding span's own length, an end off the coding sequence) stay unexplained.
//! A key that a set rule of the tuple comparator names on another of its sets keeps that
//! rule's class.

mod classify_support;
mod common;

use classify_support::{counted, divergence, key, obj, record, records, written, Case};
use common::clean;
use vep_golden_corpus::classify::intron_span::REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE as CLASS;
use vep_golden_corpus::constants::UNEXPLAINED;

/// The terms both engines write around the coding terms of a row of the shape.
const BEFORE: &str = "splice_acceptor_variant,splice_donor_variant";
const AFTER: &str = "splice_donor_5th_base_variant,intron_variant";

/// The two sets as `classify` writes a set: the terms sorted.
const VEP_SET: &str = "intron_variant,protein_altering_variant,splice_acceptor_variant,\
splice_donor_5th_base_variant,splice_donor_variant,stop_lost";
const VEP_RS_SET: &str = "inframe_deletion,intron_variant,splice_acceptor_variant,\
splice_donor_5th_base_variant,splice_donor_variant";

/// A reference allele of `bases` bases.
fn allele(bases: usize) -> String {
    "ACGT".chars().cycle().take(bases).collect()
}

/// The Extra column of a row of an edited model whose reference allele is `used_ref`.
fn edited(used_ref: &str) -> String {
    format!("IMPACT=HIGH;STRAND=-1;BAM_EDIT=OK;GIVEN_REF={used_ref};USED_REF={used_ref}")
}

/// A row of a default-format output of the deletion `X:101-130` on the transcript `feature`,
/// with the coding terms between the terms of a whole intron, its `CDS_position` and its
/// Extra column.
fn row(feature: &str, coding_terms: &str, cds_position: &str, extra: &str) -> String {
    let terms = format!("{BEFORE},{coding_terms},{AFTER}");
    let columns = [
        "del",
        "X:101-130",
        "-",
        "GENE1",
        feature,
        "Transcript",
        &terms,
        "40-51",
        cds_position,
        "4-7",
        "-",
        "-",
        "-",
        extra,
    ];
    columns.join("\t") + "\n"
}

/// The record both outputs annotate: a deletion of 30 bases.
fn held() -> String {
    let line = format!("X\t100\tdel\tG{}\tG\t.\t.\t.", allele(30));
    records(&[record(&line)])
}

/// Runs `classify` on one row a side that differ in their coding terms as the class
/// documents, and returns the class of the one entry it writes.
fn class_of(tag: &str, vep: (&str, &str), vep_rs: (&str, &str)) -> String {
    let case = Case::new(tag);
    let held = held();
    case.manifest(&obj(&[&held]));
    case.outputs(
        &row("NM_1.1", "stop_lost,protein_altering_variant", vep.0, vep.1),
        &row("NM_1.1", "inframe_deletion", vep_rs.0, vep_rs.1),
    );
    let ran = case.classify(&[]);
    for class in [CLASS, UNEXPLAINED] {
        let summary = format!("{{'{class}': 1}}");
        if ran == clean(0, &counted((1, &summary), (0, "{}"), "[]")) {
            let expected = written(
                &[divergence(
                    &key("X:101-130", "-", ("NM_1.1", "Transcript"), &[0]),
                    VEP_SET,
                    &[VEP_RS_SET],
                    class,
                )],
                &[],
                (&format!(r#"{{"{class}": 1}}"#), "{}"),
            );
            assert_eq!(case.after(), obj(&[&held, &expected]), "{tag}");
            return class.to_string();
        }
    }
    panic!("{tag}: classify ended {ran:?}");
}

#[test]
fn a_whole_intron_inside_the_coding_span_of_an_edited_model_takes_the_class() {
    // Twelve coding bases under a reference allele of 30: the allele holds the intron.
    let extra = edited(&allele(30));
    assert_eq!(
        class_of("takes", ("10-21", &extra), ("10-21", &extra)),
        CLASS
    );
}

#[test]
fn the_nearest_rows_stay_unexplained() {
    let extra = edited(&allele(30));
    // A model without an edit status, on either side: the row keeps its two alleles.
    let unedited = extra.replace("BAM_EDIT=OK;", "");
    assert_eq!(
        class_of(
            "unedited_reference",
            ("10-21", &unedited),
            ("10-21", &extra)
        ),
        UNEXPLAINED
    );
    assert_eq!(
        class_of(
            "unedited_candidate",
            ("10-21", &extra),
            ("10-21", &unedited)
        ),
        UNEXPLAINED
    );
    // A reference allele of the coding span's own length holds no intron base.
    let short = edited(&allele(12));
    assert_eq!(
        class_of("allele_of_the_span", ("10-21", &short), ("10-21", &short)),
        UNEXPLAINED
    );
    // An end of the span off the coding sequence, and two sides that place the span apart.
    assert_eq!(
        class_of("open_end", ("?-21", &extra), ("?-21", &extra)),
        UNEXPLAINED
    );
    assert_eq!(
        class_of("spans_apart", ("10-21", &extra), ("10-24", &extra)),
        UNEXPLAINED
    );
}

#[test]
fn a_key_that_a_set_rule_names_keeps_that_class() {
    // vep-rs writes two sets on the key: the class's own, which sorts first, and Ensembl
    // VEP's with `splice_region_variant` beside it, a pair the set rules name
    // `splice_family_swap`.
    let extra = edited(&allele(30));
    let vep_terms = "stop_lost,protein_altering_variant";
    let with_region = format!("{vep_terms},splice_region_variant");
    let case = Case::new("set_rule");
    let held = held();
    case.manifest(&obj(&[&held]));
    case.outputs(
        &row("NM_1.1", vep_terms, "10-21", &extra),
        &[
            row("NM_1.1", "inframe_deletion", "10-21", &extra),
            row("NM_1.1", &with_region, "10-21", &extra),
        ]
        .concat(),
    );
    let summary = "{'splice_family_swap': 1}";
    assert_eq!(
        case.classify(&[]),
        clean(0, &counted((1, summary), (0, "{}"), "[]"))
    );
    let second_set = "intron_variant,protein_altering_variant,splice_acceptor_variant,\
splice_donor_5th_base_variant,splice_donor_variant,splice_region_variant,stop_lost";
    let expected = written(
        &[divergence(
            &key("X:101-130", "-", ("NM_1.1", "Transcript"), &[0]),
            VEP_SET,
            &[VEP_RS_SET, second_set],
            "splice_family_swap",
        )],
        &[],
        (r#"{"splice_family_swap": 1}"#, "{}"),
    );
    assert_eq!(case.after(), obj(&[&held, &expected]));
}
