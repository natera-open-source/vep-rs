// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The term-pair classes, the cells that follow them and the FREQS classes through the
//! comparison of two outputs: what a class sets aside, what it counts and what stays a
//! difference beside it.

mod coding_end_common;
mod compare_common;
mod term_pairs_common;

use std::collections::BTreeMap;
use std::collections::BTreeSet;

use compare_common::*;
use indexmap::IndexMap;
use term_pairs_common::*;
use vep_formats_compare::classes;
use vep_formats_compare::cli::masks_under;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::nearest::NearestSites;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::partial_codon_models::PartialCodonModels;
use vep_formats_compare::rules::term_pairs::splice_term_past_the_record_span_pair;
use vep_formats_compare::rules::term_pairs::stop_codon_overlap_genomic_span_pair;
use vep_formats_compare::rules::term_pairs::stop_codon_overlap_insertion_range_pair;
use vep_formats_compare::rules::term_pairs::stop_codon_removed_without_stop_codon_test_pair;
use vep_formats_compare::rules::term_pairs::stop_codon_unshifted_span_pair;
use vep_formats_compare::rules::term_pairs::TERM_PAIR_CLASSES;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::ABSENT;

const POSITIONS_LEFT_OUT: &str = "\t-\t-\t-\t-\t-\t-\t";

/// The verdict of two default-format texts under the built-in `classes`.
fn verdict_under(reference: &str, candidate: &str, classes: &[&str]) -> String {
    verdict(&compare_default(reference, candidate, classes).0).to_string()
}

/// What was applied, as a count map.
fn applied(documented: &Documented) -> BTreeMap<String, i64> {
    counts(&documented.applied)
}

fn one(klass: &str) -> BTreeMap<String, i64> {
    count_map(&[(klass, 1)])
}

fn none() -> BTreeMap<String, i64> {
    BTreeMap::new()
}

/// The classes a column's differences were set aside under, each with its count; empty for
/// a column without any.
fn masked(report: &Value, column: &str) -> BTreeMap<String, i64> {
    let Some(masked) = at(report, &["columns", column]).get("masked") else {
        return BTreeMap::new();
    };
    masked
        .as_dict()
        .expect("the masked counts are an object")
        .iter()
        .map(|(klass, count)| match count {
            Value::Int(count) => (klass.clone(), *count),
            other => panic!("{klass} is counted {other:?}"),
        })
        .collect()
}

fn has_masked(report: &Value, column: &str) -> bool {
    at(report, &["columns", column]).get("masked").is_some()
}

fn differing(report: &Value, column: &str) -> i64 {
    number(report, &["columns", column, "differing"])
}

/// The report of two entry lists under the built-in `classes`, and what was applied.
fn compared(expected: &[Entry], actual: &[Entry], classes: &[&str]) -> (Value, Documented) {
    let mut documented = documented_with(classes);
    let report = compare(expected, actual, &mut documented);
    (report, documented)
}

fn vcf_entries(text: &str, inputs: &[String]) -> Vec<Entry> {
    parse_vcf(text, inputs).expect("the output parses")
}

fn json_entries(text: &str, inputs: &[String]) -> Vec<Entry> {
    parse_json(text, inputs).expect("the output parses")
}

fn tab_entries(text: &str) -> Vec<Entry> {
    parse_tab(text).expect("the output parses")
}

/// One JSON object of the record `input` with one transcript consequence, `entry`.
fn json_object(
    input: &str,
    id: &str,
    allele_string: &str,
    most_severe: &str,
    entry: &str,
) -> String {
    format!(
        "{{\"input\": {}, \"id\": {}, \"allele_string\": {}, \"most_severe_consequence\": {}, \"transcript_consequences\": [{entry}]}}\n",
        quoted(input),
        quoted(id),
        quoted(allele_string),
        quoted(most_severe)
    )
}

const MIRNA_CLASS: &str = classes::MATURE_MIRNA_UNSHIFTED_SPAN;
const MIRNA_INPUT: &str = "21\t8208517\t21:8208517:A:AG\tA\tAG\t.\t.\t.";

/// Under --shift_3prime 1 --shift_length, two Consequence sets equal but for one engine's
/// mature_miRNA_variant against the other's non_coding_transcript_exon_variant, on a row with
/// a positive SHIFT_LENGTH, is the class in either direction, counted once per row and
/// leaving every other field compared; a zero or absent shift, a set differing by another
/// term or by these two beside another, and a differing HGVSc stay differences.
#[test]
fn mature_mirna_unshifted_span_accepts_the_mature_term_read_at_the_unshifted_span() {
    let klass = MIRNA_CLASS;
    let other =
        |text: &str| text.replace("mature_miRNA_variant", "non_coding_transcript_exon_variant");
    let extra =
        "IMPACT=MODIFIER;STRAND=1;HGVSc=ENST00000615959.1:n.48dup;HGVS_OFFSET=3;SHIFT_LENGTH=3";
    let reference = default_text(&[default_row(
        "21_8208517_A/AG",
        "21:8208520-8208521",
        "G",
        "ENST00000615959",
        "mature_miRNA_variant",
        extra,
    )]);
    let candidate = other(&reference);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "Consequence"), one(klass));
    assert_eq!(verdict_under(&reference, &candidate, &[]), "red");
    // The shift may carry the span onto the mature sequence as well as off it.
    assert_eq!(verdict_under(&candidate, &reference, &[klass]), "green");
    let deletion = reference
        .replace("\tG\t", "\t-\t")
        .replace("21:8208520-8208521", "21:8208524-8208543")
        .replace("n.48dup", "n.52_71del")
        .replace(
            "HGVS_OFFSET=3;SHIFT_LENGTH=3",
            "HGVS_OFFSET=4;SHIFT_LENGTH=4",
        );
    assert_eq!(
        verdict_under(&deletion, &other(&deletion), &[klass]),
        "green"
    );
    for variant in ["SHIFT_LENGTH=0", ""] {
        let replaced = reference.replace("SHIFT_LENGTH=3", variant);
        let unshifted = format!("{}\n", replaced.trim_end_matches([';', '\n']));
        assert_eq!(
            verdict_under(&unshifted, &other(&unshifted), &[klass]),
            "red",
            "{variant}"
        );
    }
    for terms in [
        "non_coding_transcript_variant",
        "intron_variant,non_coding_transcript_variant",
        "mature_miRNA_variant,non_coding_transcript_exon_variant",
        "splice_region_variant,non_coding_transcript_exon_variant",
    ] {
        let differing = reference.replace("mature_miRNA_variant", terms);
        assert_eq!(
            verdict_under(&reference, &differing, &[klass]),
            "red",
            "{terms}"
        );
    }
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace("n.48dup", "n.47dup"),
            &[klass]
        ),
        "red"
    );
}

/// The class counts a row once in the VCF and the JSON format too.
#[test]
fn mature_mirna_unshifted_span_counts_a_row_once_in_the_vcf_and_json_formats() {
    let klass = MIRNA_CLASS;
    let inputs = [MIRNA_INPUT.to_string()];
    let vcf = |term: &str| {
        let text = format!(
            "{}{}CSQ=G|{term}|ENST00000615959|Transcript|ENST00000615959.1:n.48dup|3|3\n",
            csq_header("Allele|Consequence|Feature|Feature_type|HGVSc|HGVS_OFFSET|SHIFT_LENGTH"),
            open_info(MIRNA_INPUT)
        );
        vcf_entries(&text, &inputs)
    };
    let (report, documented) = compared(
        &vcf("mature_miRNA_variant"),
        &vcf("non_coding_transcript_exon_variant"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    let obj = |term: &str| {
        let entry = format!(
            "{{\"variant_allele\": \"G\", \"transcript_id\": \"ENST00000615959\", \"consequence_terms\": [{}], \"hgvsc\": \"ENST00000615959.1:n.48dup\", \"hgvs_offset\": 3, \"shift_length\": 3}}",
            quoted(term)
        );
        json_entries(
            &json_object(MIRNA_INPUT, "21:8208517:A:AG", "-/G", term, &entry),
            &inputs,
        )
    };
    let (report, documented) = compared(
        &obj("mature_miRNA_variant"),
        &obj("non_coding_transcript_exon_variant"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
}

const UNSHIFTED_CLASS: &str = classes::STOP_CODON_UNSHIFTED_SPAN;

/// The bases the record of the class's row deletes, each written `N`.
fn unshifted_deleted() -> String {
    "N".repeat(3804)
}

/// The input line of the class's row: the base the line keeps and the bases it deletes.
fn unshifted_input() -> String {
    format!("16\t172870\t38636\tN{}\tN\t.\t.\t.", unshifted_deleted())
}

/// Under --shift_3prime 1 --shift_length, two Consequence sets equal but for the reference's
/// coding_sequence_variant against the candidate's stop_retained_variant, on a row with a
/// positive SHIFT_LENGTH, is the class, counted once per row and leaving every other field
/// compared; the terms swapped, a zero or absent shift, a set differing by another term or by
/// these two beside another, and a differing HGVSc stay differences.
#[test]
fn stop_codon_unshifted_span_accepts_the_coding_term_read_at_the_unshifted_span() {
    let klass = UNSHIFTED_CLASS;
    let other = |text: &str| text.replace("coding_sequence_variant", "stop_retained_variant");
    let extra =
        "IMPACT=HIGH;STRAND=1;HGVSc=ENST00001097586.1:c.-37_*206del;HGVS_OFFSET=513;SHIFT_LENGTH=513";
    let reference = default_text(&[default_row(
        "38636",
        "16:173384-177187",
        "-",
        "ENST00001097586",
        "splice_acceptor_variant,splice_donor_variant,splice_donor_5th_base_variant,coding_sequence_variant,5_prime_UTR_variant,3_prime_UTR_variant,intron_variant",
        extra,
    )]);
    let candidate = other(&reference);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "Consequence"), one(klass));
    assert_eq!(verdict_under(&reference, &candidate, &[]), "red");
    assert_eq!(verdict_under(&candidate, &reference, &[klass]), "red");
    for variant in ["SHIFT_LENGTH=0", ""] {
        let replaced = reference.replace("SHIFT_LENGTH=513", variant);
        let unshifted = format!("{}\n", replaced.trim_end_matches([';', '\n']));
        assert_eq!(
            verdict_under(&unshifted, &other(&unshifted), &[klass]),
            "red",
            "{variant}"
        );
    }
    for terms in [
        "stop_lost,",
        "stop_retained_variant,coding_sequence_variant,",
        "stop_retained_variant,splice_region_variant,",
        "",
    ] {
        let differing = reference.replace("coding_sequence_variant,", terms);
        assert_eq!(
            verdict_under(&reference, &differing, &[klass]),
            "red",
            "{terms}"
        );
    }
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace("c.-37_*206del", "c.-37_*205del"),
            &[klass]
        ),
        "red"
    );
}

/// The class counts a row once in the VCF and the JSON format too.
#[test]
fn stop_codon_unshifted_span_counts_a_row_once_in_the_vcf_and_json_formats() {
    let klass = UNSHIFTED_CLASS;
    let input = unshifted_input();
    let inputs = [input.clone()];
    let vcf = |term: &str| {
        let text = format!(
            "{}{}CSQ=-|splice_acceptor_variant&splice_donor_variant&splice_donor_5th_base_variant&{term}&5_prime_UTR_variant&3_prime_UTR_variant&intron_variant|ENST00001097586|Transcript|ENST00001097586.1:c.-37_*206del|513|513\n",
            csq_header("Allele|Consequence|Feature|Feature_type|HGVSc|HGVS_OFFSET|SHIFT_LENGTH"),
            open_info(&input)
        );
        vcf_entries(&text, &inputs)
    };
    let (report, documented) = compared(
        &vcf("coding_sequence_variant"),
        &vcf("stop_retained_variant"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    let obj = |term: &str| {
        let entry = format!(
            "{{\"variant_allele\": \"-\", \"transcript_id\": \"ENST00001097586\", \"consequence_terms\": [\"splice_acceptor_variant\", \"splice_donor_variant\", \"splice_donor_5th_base_variant\", {}, \"5_prime_UTR_variant\", \"3_prime_UTR_variant\", \"intron_variant\"], \"hgvsc\": \"ENST00001097586.1:c.-37_*206del\", \"hgvs_offset\": 513, \"shift_length\": 513}}",
            quoted(term)
        );
        json_entries(
            &json_object(
                &input,
                "38636",
                &format!("{}/-", unshifted_deleted()),
                "splice_acceptor_variant",
                &entry,
            ),
            &inputs,
        )
    };
    let (report, documented) = compared(
        &obj("coding_sequence_variant"),
        &obj("stop_retained_variant"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
}

const DELETED_SPAN_CLASS: &str = classes::REFSEQ_EDITED_MODEL_DELETED_SPAN_TERM;
const INTERGENIC: &str = "\tintergenic_variant\t";

fn with_term(text: &str, terms: &str) -> String {
    text.replace(INTERGENIC, &format!("\t{terms}\t"))
}

/// A BAM_EDIT OK transcript row of a sequence allele inside the bases the model's alignment
/// deleted: intergenic_variant in the reference, exactly one of coding_sequence_variant,
/// 5_prime_UTR_variant and 3_prime_UTR_variant in the candidate, no cDNA position and one
/// IMPACT on both, counted once per row in the default and the tab format with every other
/// field compared. A row without BAM_EDIT OK, a cDNA position on either side, a differing
/// IMPACT, a structural allele, a regulatory row, two terms or any other term, and a swap of
/// the two engines stay differences.
#[test]
fn refseq_edited_model_deleted_span_term_accepts_the_default_term_inside_a_deleted_edit_span() {
    let klass = DELETED_SPAN_CLASS;
    let extra = "IMPACT=MODIFIER;STRAND=1;GIVEN_REF=C;USED_REF=C;BAM_EDIT=OK";
    let deleted_extra = extra.replace(
        "GIVEN_REF=C;USED_REF=C",
        "GIVEN_REF=CTGCTGCGCCCCCAG;USED_REF=CTGCTGCGCCCCCAG",
    );
    let mapped_extra = extra.replace(
        "GIVEN_REF=C;USED_REF=C",
        "GIVEN_REF=CCTGCTGCGCCCCCA;USED_REF=CCTGCTGCGCCCCCA",
    );
    let row = |name: &str, location: &str, allele: &str, terms: &str, extra: &str| {
        default_row_of(
            name,
            location,
            allele,
            "NM_198689.3",
            terms,
            extra,
            "Transcript",
            "386675",
        )
    };
    let inside = row(
        "rs373191083",
        "21:46020662",
        "T",
        "intergenic_variant",
        extra,
    );
    let whole_span = row(
        "rs782183638",
        "21:46020656-46020670",
        "-",
        "intergenic_variant",
        &deleted_extra,
    );
    let reaching = row(
        "rs782330177",
        "21:46020655-46020669",
        "-",
        "coding_sequence_variant",
        &mapped_extra,
    )
    .replacen(POSITIONS_LEFT_OUT, "\t159-?\t134-?\t45-?\t-\t-\t-\t", 1);
    let reference = default_text(&[inside.as_str(), &whole_span, &reaching]);
    let candidate = with_term(&reference, "coding_sequence_variant");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(klass, 2)]));
    assert_eq!(masked(&report, "Consequence"), count_map(&[(klass, 2)]));
    assert_eq!(number(&report, &["counts", "compared"]), 3);
    assert_eq!(
        texts(at(&report, &["mask_examples", klass]))[0],
        "21:46020662|T|NM_198689.3|Transcript: Consequence expected 'intergenic_variant', got 'coding_sequence_variant'"
    );
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    assert_eq!(differing(&without, "Consequence"), 2);
    assert_eq!(verdict_under(&candidate, &reference, &[klass]), "red");
    for term in ["3_prime_UTR_variant", "5_prime_UTR_variant"] {
        assert_eq!(
            verdict_under(&reference, &with_term(&reference, term), &[klass]),
            "green",
            "{term}"
        );
    }
    for terms in [
        "missense_variant",
        "intron_variant",
        "coding_sequence_variant,3_prime_UTR_variant",
        "coding_sequence_variant,intergenic_variant",
        "-",
    ] {
        assert_eq!(
            verdict_under(&reference, &with_term(&reference, terms), &[klass]),
            "red",
            "{terms}"
        );
    }
    let unedited = reference.replace(";BAM_EDIT=OK", "");
    assert_eq!(
        verdict_under(
            &unedited,
            &with_term(&unedited, "coding_sequence_variant"),
            &[klass]
        ),
        "red"
    );
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace(";BAM_EDIT=OK", ";BAM_EDIT=FAILED"),
            &[klass]
        ),
        "red"
    );
    let positioned =
        default_text(&[inside.replacen(POSITIONS_LEFT_OUT, "\t165\t-\t-\t-\t-\t-\t", 1)]);
    let positioned_candidate = with_term(&positioned, "coding_sequence_variant");
    assert_eq!(
        verdict_under(&positioned, &positioned_candidate, &[klass]),
        "red"
    );
    assert_eq!(
        verdict_under(
            &default_text(&[inside.as_str()]),
            &positioned_candidate,
            &[klass]
        ),
        "red"
    );
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace("IMPACT=MODIFIER", "IMPACT=LOW"),
            &[klass]
        ),
        "red"
    );
    let structural = reference.replace("\t-\t386675", "\tdeletion\t386675");
    assert_eq!(
        verdict_under(
            &structural,
            &with_term(&structural, "coding_sequence_variant"),
            &[klass]
        ),
        "red"
    );
    let regulatory = default_text(&[default_row_of(
        "rs_reg",
        "21:46020662",
        "T",
        "ENSR1",
        "intergenic_variant",
        extra,
        "RegulatoryFeature",
        "-",
    )]);
    assert_eq!(
        verdict_under(
            &regulatory,
            &with_term(&regulatory, "coding_sequence_variant"),
            &[klass]
        ),
        "red"
    );
    // Every other field of the accepted row stays compared.
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace("STRAND=1;GIVEN_REF=C;", "STRAND=-1;GIVEN_REF=C;"),
            &[klass]
        ),
        "red"
    );
    // The tab format carries the term and the edit flag as fields of their own.
    let tab = |term: &str| {
        tab_entries(&format!(
            "{}\nrs373191083\t21:46020662\tT\t386675\tNM_198689.3\tTranscript\t{term}\t-\t-\t-\t-\t-\t-\tMODIFIER\t1\tOK\n",
            DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tBAM_EDIT")
        ))
    };
    let (report, documented) = compared(
        &tab("intergenic_variant"),
        &tab("coding_sequence_variant"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
}

/// The class counts a row once in the VCF and the JSON format, where the JSON record's
/// most_severe_consequence follows the term and is alone left uncompared.
#[test]
fn refseq_edited_model_deleted_span_term_counts_a_row_once_in_the_vcf_and_json_formats() {
    let klass = DELETED_SPAN_CLASS;
    let input = "21\t46020662\trs373191083\tC\tT\t.\t.\t.";
    let inputs = [input.to_string()];
    let vcf = |term: &str| {
        let text = format!(
            "{}{}CSQ=T|{term}|MODIFIER|NM_198689.3|Transcript||OK\n",
            csq_header("Allele|Consequence|IMPACT|Feature|Feature_type|cDNA_position|BAM_EDIT"),
            open_info(input)
        );
        vcf_entries(&text, &inputs)
    };
    let (report, documented) = compared(
        &vcf("intergenic_variant"),
        &vcf("coding_sequence_variant"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    let obj = |term: &str| {
        let entry = format!(
            "{{\"variant_allele\": \"T\", \"transcript_id\": \"NM_198689.3\", \"gene_id\": \"386675\", \"impact\": \"MODIFIER\", \"consequence_terms\": [{}], \"bam_edit\": \"OK\", \"strand\": 1}}",
            quoted(term)
        );
        json_object(input, "rs373191083", "C/T", term, &entry)
    };
    let reference = json_entries(&obj("intergenic_variant"), &inputs);
    let candidate = json_entries(&obj("coding_sequence_variant"), &inputs);
    let (report, documented) = compared(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    // The record's most_severe_consequence alone is left uncompared: its other keys and the
    // entry's stay compared.
    for (written, changed) in [
        ("\"id\": \"rs373191083\"", "\"id\": \"rs1\""),
        ("\"strand\": 1", "\"strand\": -1"),
    ] {
        let text = obj("coding_sequence_variant").replace(written, changed);
        let (report, _) = compared(&reference, &json_entries(&text, &inputs), &[klass]);
        assert_eq!(verdict(&report), "red", "{changed}");
    }
}

const TAIL_CLASS: &str = classes::REFSEQ_EDITED_MODEL_TAIL_TERM;

/// A BAM_EDIT OK transcript row of a sequence allele on one of the final bases of a model
/// whose edits inserted bases: intergenic_variant in the reference, 3_prime_UTR_variant in the
/// candidate, a cDNA position and one IMPACT on both, counted once per row in the default and
/// the tab format with every other field compared. A row without BAM_EDIT OK or without a
/// cDNA position (the deleted-span shape), a differing IMPACT, a structural allele, a
/// regulatory row, any other candidate term and a swap of the two engines stay differences.
#[test]
fn refseq_edited_model_tail_term_accepts_the_default_term_on_a_lengthened_models_final_bases() {
    let klass = TAIL_CLASS;
    let extra = "IMPACT=MODIFIER;STRAND=-1;GIVEN_REF=T;USED_REF=T;BAM_EDIT=OK";
    let tail = default_row_of(
        "rs2082092473",
        "21:46301131",
        "C",
        "NM_001286463.1",
        "intergenic_variant",
        extra,
        "Transcript",
        "54058",
    )
    .replacen(POSITIONS_LEFT_OUT, "\t1988\t-\t-\t-\t-\t-\t", 1);
    let inside = default_row_of(
        "rs_in",
        "21:46301132",
        "A",
        "NM_001286463.1",
        "3_prime_UTR_variant",
        &extra.replace("GIVEN_REF=T;USED_REF=T", "GIVEN_REF=G;USED_REF=G"),
        "Transcript",
        "54058",
    )
    .replacen(POSITIONS_LEFT_OUT, "\t1987\t-\t-\t-\t-\t-\t", 1);
    let reference = default_text(&[tail.as_str(), &inside]);
    let candidate = with_term(&reference, "3_prime_UTR_variant");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "Consequence"), one(klass));
    assert_eq!(number(&report, &["counts", "compared"]), 2);
    assert_eq!(
        texts(at(&report, &["mask_examples", klass])),
        ["21:46301131|C|NM_001286463.1|Transcript: Consequence expected 'intergenic_variant', got '3_prime_UTR_variant'"]
    );
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    assert_eq!(differing(&without, "Consequence"), 1);
    assert_eq!(verdict_under(&candidate, &reference, &[klass]), "red");
    for terms in [
        "5_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "coding_sequence_variant",
        "3_prime_UTR_variant,intergenic_variant",
        "-",
    ] {
        assert_eq!(
            verdict_under(&reference, &with_term(&reference, terms), &[klass]),
            "red",
            "{terms}"
        );
    }
    let unedited = reference.replace(";BAM_EDIT=OK", "");
    assert_eq!(
        verdict_under(
            &unedited,
            &with_term(&unedited, "3_prime_UTR_variant"),
            &[klass]
        ),
        "red"
    );
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace(";BAM_EDIT=OK", ";BAM_EDIT=FAILED"),
            &[klass]
        ),
        "red"
    );
    // Without a cDNA position the row is the deleted-span class's shape, not this one.
    let unpositioned = reference.replace("\t1988\t", "\t-\t");
    let unpositioned_candidate = with_term(&unpositioned, "3_prime_UTR_variant");
    assert_eq!(
        verdict_under(&unpositioned, &unpositioned_candidate, &[klass]),
        "red"
    );
    assert_eq!(
        verdict_under(
            &unpositioned,
            &unpositioned_candidate,
            &[DELETED_SPAN_CLASS]
        ),
        "green"
    );
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace("\t1988\t", "\t-\t"),
            &[klass]
        ),
        "red"
    );
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace(
                "IMPACT=MODIFIER;STRAND=-1;GIVEN_REF=T",
                "IMPACT=LOW;STRAND=-1;GIVEN_REF=T"
            ),
            &[klass]
        ),
        "red"
    );
    let structural = reference.replace("\tC\t54058", "\tdeletion\t54058");
    assert_eq!(
        verdict_under(
            &structural,
            &with_term(&structural, "3_prime_UTR_variant"),
            &[klass]
        ),
        "red"
    );
    let regulatory = default_text(&[default_row_of(
        "rs_reg",
        "21:46301131",
        "C",
        "ENSR1",
        "intergenic_variant",
        extra,
        "RegulatoryFeature",
        "-",
    )
    .replacen(POSITIONS_LEFT_OUT, "\t1988\t-\t-\t-\t-\t-\t", 1)]);
    assert_eq!(
        verdict_under(
            &regulatory,
            &with_term(&regulatory, "3_prime_UTR_variant"),
            &[klass]
        ),
        "red"
    );
    // Every other field of the accepted row stays compared.
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace("STRAND=-1;GIVEN_REF=T;", "STRAND=1;GIVEN_REF=T;"),
            &[klass]
        ),
        "red"
    );
    // The tab format carries the term, the position and the edit flag as fields of their own.
    let tab = |term: &str| {
        tab_entries(&format!(
            "{}\nrs2082092473\t21:46301131\tC\t54058\tNM_001286463.1\tTranscript\t{term}\t1988\t-\t-\t-\t-\t-\tMODIFIER\t-1\tOK\n",
            DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tBAM_EDIT")
        ))
    };
    let (report, documented) = compared(
        &tab("intergenic_variant"),
        &tab("3_prime_UTR_variant"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
}

/// The class counts a row once in the VCF and the JSON format, where the JSON record's
/// most_severe_consequence follows the term and is alone left uncompared.
#[test]
fn refseq_edited_model_tail_term_counts_a_row_once_in_the_vcf_and_json_formats() {
    let klass = TAIL_CLASS;
    let input = "21\t46301131\trs2082092473\tT\tC\t.\t.\t.";
    let inputs = [input.to_string()];
    let vcf = |term: &str| {
        let text = format!(
            "{}{}CSQ=C|{term}|MODIFIER|NM_001286463.1|Transcript|1988|OK\n",
            csq_header("Allele|Consequence|IMPACT|Feature|Feature_type|cDNA_position|BAM_EDIT"),
            open_info(input)
        );
        vcf_entries(&text, &inputs)
    };
    let (report, documented) = compared(
        &vcf("intergenic_variant"),
        &vcf("3_prime_UTR_variant"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    let obj = |term: &str| {
        let entry = format!(
            "{{\"variant_allele\": \"C\", \"transcript_id\": \"NM_001286463.1\", \"gene_id\": \"54058\", \"impact\": \"MODIFIER\", \"consequence_terms\": [{}], \"bam_edit\": \"OK\", \"strand\": -1, \"cdna_start\": 1988, \"cdna_end\": 1988}}",
            quoted(term)
        );
        json_object(input, "rs2082092473", "T/C", term, &entry)
    };
    let reference = json_entries(&obj("intergenic_variant"), &inputs);
    let candidate = json_entries(&obj("3_prime_UTR_variant"), &inputs);
    let (report, documented) = compared(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    // The record's most_severe_consequence alone is left uncompared: its other keys and the
    // entry's stay compared.
    for (written, changed) in [
        ("\"id\": \"rs2082092473\"", "\"id\": \"rs1\""),
        ("\"strand\": -1", "\"strand\": 1"),
    ] {
        let text = obj("3_prime_UTR_variant").replace(written, changed);
        let (report, _) = compared(&reference, &json_entries(&text, &inputs), &[klass]);
        assert_eq!(verdict(&report), "red", "{changed}");
    }
}

const EQUAL_SPAN_CLASS: &str = classes::REFSEQ_EDITED_MODEL_EQUAL_SPAN_INSERTION_TERM;
const EQUAL_SPAN_POSITIONS: &str = "\t186-188\t75-77\t25-26\t-\t-\t-\t";

/// The candidate's rows for the reference's: `missense_variant` and its IMPACT where the
/// reference writes `intergenic_variant` and its own.
fn equal_span_candidate(reference: &str) -> String {
    with_term(reference, "missense_variant").replace("IMPACT=MODIFIER;", "IMPACT=MODERATE;")
}

/// The default-format text of two rows of one edited model: an insertion of three bases
/// whose CDS span holds three, `intergenic_variant` in the reference, and an insertion of two
/// bases at the same point, `frameshift_variant` on both sides.
fn equal_span_reference() -> String {
    let extra = "IMPACT=MODIFIER;STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK";
    let on_the_added_bases = default_row_of(
        "551063",
        "2:73385942-73385943",
        "GGC",
        "NM_015120.4",
        "intergenic_variant",
        extra,
        "Transcript",
        "7840",
    )
    .replacen(POSITIONS_LEFT_OUT, EQUAL_SPAN_POSITIONS, 1);
    let two_bases = default_row_of(
        "1167157",
        "2:73385942-73385943",
        "GG",
        "NM_015120.4",
        "frameshift_variant",
        &extra.replace("IMPACT=MODIFIER", "IMPACT=HIGH"),
        "Transcript",
        "7840",
    )
    .replacen(POSITIONS_LEFT_OUT, EQUAL_SPAN_POSITIONS, 1);
    default_text(&[on_the_added_bases, two_bases])
}

/// A BAM_EDIT OK transcript row of an insertion between two neighbouring positions whose CDS
/// span, the same on both sides, holds as many bases as the allele: intergenic_variant in
/// the reference, missense_variant in the candidate, counted once per row in the default and
/// the tab format, the IMPACT of each side's term set aside with it and every other field
/// compared. A row without BAM_EDIT OK, an allele of another length than its CDS span, a row
/// without a CDS span or with another on the candidate's side, a Location that is not two
/// neighbouring positions, a deletion, a structural allele, a regulatory row, any other term
/// on either side and a swap of the two engines stay differences, and so does an IMPACT that
/// is not its side's term's.
#[test]
fn refseq_edited_model_equal_span_insertion_term_accepts_the_default_term_on_an_insertion_as_long_as_its_span(
) {
    let klass = EQUAL_SPAN_CLASS;
    let reference = equal_span_reference();
    let candidate = equal_span_candidate(&reference);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "Consequence"), one(klass));
    assert_eq!(masked(&report, "IMPACT"), one(klass));
    assert_eq!(masked(&report, "Extra"), one(klass));
    assert_eq!(number(&report, &["counts", "compared"]), 2);
    assert!(texts(at(&report, &["mask_examples", klass])).contains(
        &"2:73385942-73385943|GGC|NM_015120.4|Transcript: Consequence expected 'intergenic_variant', got 'missense_variant'"
    ));
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    assert_eq!(differing(&without, "Consequence"), 1);
    assert_eq!(differing(&without, "IMPACT"), 1);
    assert_eq!(verdict_under(&candidate, &reference, &[klass]), "red");
    // Any other term on the candidate's side, with the IMPACT that term gives.
    for (terms, impact) in [
        ("synonymous_variant", "LOW"),
        ("coding_sequence_variant", "MODIFIER"),
        ("3_prime_UTR_variant", "MODIFIER"),
        ("inframe_insertion", "MODERATE"),
        ("missense_variant,splice_region_variant", "MODERATE"),
        ("intergenic_variant,missense_variant", "MODERATE"),
        ("-", "MODIFIER"),
    ] {
        let other =
            with_term(&reference, terms).replace("IMPACT=MODIFIER;", &format!("IMPACT={impact};"));
        assert_eq!(
            verdict_under(&reference, &other, &[klass]),
            "red",
            "{terms}"
        );
    }
    // Any other term on the reference's side.
    let coding = with_term(&reference, "coding_sequence_variant");
    assert_eq!(verdict_under(&coding, &candidate, &[klass]), "red");
    let unedited = reference.replace(";BAM_EDIT=OK", "");
    assert_eq!(
        verdict_under(&unedited, &equal_span_candidate(&unedited), &[klass]),
        "red"
    );
    // The edit flag is read on the candidate's side too: with another there the class counts
    // no row.
    let (report, documented) = compare_default(
        &reference,
        &candidate.replace(";BAM_EDIT=OK", ";BAM_EDIT=FAILED"),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    // The allele is as long as the CDS span: four bases on three, three bases on two and
    // three on one position are not the shape.
    for (written, changed) in [
        ("\tGGC\t7840", "\tGGCA\t7840"),
        ("\t75-77\t", "\t75-76\t"),
        ("\t75-77\t", "\t75\t"),
        ("\t75-77\t", "\t?-77\t"),
        ("\t75-77\t25-26\t", "\t-\t-\t"),
    ] {
        let other = reference.replacen(written, changed, 1);
        assert_eq!(
            verdict_under(&other, &equal_span_candidate(&other), &[klass]),
            "red",
            "{changed}"
        );
    }
    // The span is one span, the same on both sides: with another on the candidate's side the
    // class counts no row.
    let (report, documented) = compare_default(
        &reference,
        &candidate.replacen("\t75-77\t", "\t76-78\t", 1),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    // Three bases that replace three are no insertion, and neither is a deletion.
    for (written, changed) in [
        ("2:73385942-73385943\tGGC", "2:73385942-73385944\tGGC"),
        ("2:73385942-73385943\tGGC", "2:73385943\tGGC"),
        ("2:73385942-73385943\tGGC", "2:73385942-73385943\t-"),
        ("2:73385942-73385943\tGGC", "2:73385942-73385943\tinsertion"),
    ] {
        let other = reference.replacen(written, changed, 1);
        assert_eq!(
            verdict_under(&other, &equal_span_candidate(&other), &[klass]),
            "red",
            "{changed}"
        );
    }
    // A deletion's allele is one character and no base: on a CDS span of one position it is
    // as long as the span, and the class counts no row.
    let deleted =
        reference
            .replacen("\tGGC\t7840", "\t-\t7840", 1)
            .replacen("\t75-77\t", "\t75\t", 1);
    let (report, documented) = compare_default(&deleted, &equal_span_candidate(&deleted), &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    // Under an identifier, without the records of the input, a two-base allele between two
    // neighbouring positions cannot be told from two replaced bases, and is left a difference.
    let two =
        reference
            .replacen("\tGGC\t7840", "\tGG\t7840", 1)
            .replacen("\t75-77\t", "\t75-76\t", 1);
    assert_eq!(
        verdict_under(&two, &equal_span_candidate(&two), &[klass]),
        "red"
    );
    let regulatory = default_text(&[default_row_of(
        "551063",
        "2:73385942-73385943",
        "GGC",
        "ENSR1",
        "intergenic_variant",
        "IMPACT=MODIFIER;BAM_EDIT=OK",
        "RegulatoryFeature",
        "-",
    )
    .replacen(POSITIONS_LEFT_OUT, EQUAL_SPAN_POSITIONS, 1)]);
    assert_eq!(
        verdict_under(&regulatory, &equal_span_candidate(&regulatory), &[klass]),
        "red"
    );
    // An IMPACT that is not its side's term's stays a difference on a row the class counts.
    for planted in ["IMPACT=LOW;", ""] {
        let (report, documented) = compare_default(
            &reference,
            &candidate.replace("IMPACT=MODERATE;", planted),
            &[klass],
        );
        assert_eq!(verdict(&report), "red", "{planted}");
        assert_eq!(applied(&documented), one(klass), "{planted}");
        assert_eq!(differing(&report, "Extra"), 1, "{planted}");
        assert!(!has_masked(&report, "Extra"), "{planted}");
    }
    // Every other field of the accepted row stays compared.
    for (written, changed) in [
        (
            "STRAND=1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK",
            "STRAND=-1;GIVEN_REF=-;USED_REF=-;BAM_EDIT=OK",
        ),
        ("\t186-188\t75-77\t", "\t187-189\t75-77\t"),
        ("\t75-77\t25-26\t", "\t75-77\t25\t"),
    ] {
        assert_eq!(
            verdict_under(
                &reference,
                &candidate.replacen(written, changed, 1),
                &[klass]
            ),
            "red",
            "{changed}"
        );
    }
    // The tab format carries the term, the positions, the IMPACT and the edit flag as fields
    // of their own.
    let tab = |term: &str, impact: &str| {
        tab_entries(&format!(
            "{}\n551063\t2:73385942-73385943\tGGC\t7840\tNM_015120.4\tTranscript\t{term}\t186-188\t75-77\t25-26\t-\t-\t-\t{impact}\t1\tOK\n",
            DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tBAM_EDIT")
        ))
    };
    let (report, documented) = compared(
        &tab("intergenic_variant", "MODIFIER"),
        &tab("missense_variant", "MODERATE"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "IMPACT"), one(klass));
    let (report, documented) = compared(
        &tab("intergenic_variant", "MODIFIER"),
        &tab("missense_variant", "HIGH"),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), one(klass));
}

/// The class counts a row once in the VCF and the JSON format, where the record's own
/// alleles show the insertion, the IMPACT follows the two terms and the JSON record's
/// most_severe_consequence is alone left uncompared.
#[test]
fn refseq_edited_model_equal_span_insertion_term_counts_a_row_once_in_the_vcf_and_json_formats() {
    let klass = EQUAL_SPAN_CLASS;
    let input = "2\t73385942\t551063\tA\tAGGC\t.\t.\t.";
    let inputs = [input.to_string()];
    let vcf_of = |input: &str, allele: &str, term: &str, impact: &str| {
        let text = format!(
            "{}{}CSQ={allele}|{term}|{impact}|NM_015120.4|Transcript|186-188|75-77|25-26|OK\n",
            csq_header(
                "Allele|Consequence|IMPACT|Feature|Feature_type|cDNA_position|CDS_position|Protein_position|BAM_EDIT"
            ),
            open_info(input)
        );
        vcf_entries(&text, &[input.to_string()])
    };
    let vcf = |term: &str, impact: &str| vcf_of(input, "GGC", term, impact);
    let (report, documented) = compared(
        &vcf("intergenic_variant", "MODIFIER"),
        &vcf("missense_variant", "MODERATE"),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(
        verdict(&compare(
            &vcf("intergenic_variant", "MODIFIER"),
            &vcf("missense_variant", "MODERATE"),
            &mut no_documented()
        )),
        "red"
    );
    // Three bases that replace the record's three are no insertion.
    let replaced = "2\t73385942\t551063\tAAAA\tAGGC\t.\t.\t.";
    let (report, documented) = compared(
        &vcf_of(replaced, "GGC", "intergenic_variant", "MODIFIER"),
        &vcf_of(replaced, "GGC", "missense_variant", "MODERATE"),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    let obj = |term: &str, impact: &str| {
        let entry = format!(
            "{{\"variant_allele\": \"GGC\", \"transcript_id\": \"NM_015120.4\", \"gene_id\": \"7840\", \"impact\": {}, \"consequence_terms\": [{}], \"bam_edit\": \"OK\", \"strand\": 1, \"cdna_start\": 186, \"cdna_end\": 188, \"cds_start\": 75, \"cds_end\": 77, \"protein_start\": 25, \"protein_end\": 26}}",
            quoted(impact),
            quoted(term)
        );
        json_object(input, "551063", "-/GGC", term, &entry)
    };
    let reference = json_entries(&obj("intergenic_variant", "MODIFIER"), &inputs);
    let candidate = json_entries(&obj("missense_variant", "MODERATE"), &inputs);
    let (report, documented) = compared(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    // The record's most_severe_consequence and the entry's impact alone are set aside: the
    // record's other keys and the entry's stay compared.
    for (written, changed) in [
        ("\"id\": \"551063\"", "\"id\": \"rs1\""),
        ("\"strand\": 1", "\"strand\": -1"),
        ("\"cds_end\": 77", "\"cds_end\": 78"),
        ("\"impact\": \"MODERATE\"", "\"impact\": \"LOW\""),
    ] {
        let text = obj("missense_variant", "MODERATE").replace(written, changed);
        let (report, _) = compared(&reference, &json_entries(&text, &inputs), &[klass]);
        assert_eq!(verdict(&report), "red", "{changed}");
    }
}

const CODONS_CLASS: &str = classes::REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS;
const REFERENCE_CODONS: &str = "\tcCAGGCCCC/c\t";
const CANDIDATE_CODONS: &str = "\tcCAGGCCCA/c\t";

fn with_candidate_codons(text: &str) -> String {
    text.replace(REFERENCE_CODONS, CANDIDATE_CODONS)
}

/// The Codons of a BAM_EDIT OK transcript row of a sequence allele whose USED_REF has more
/// bases than its CDS span: the alternate codons equal, the reference codons of one length
/// agreeing at every lower-case base and differing in an upper-case one, the reference's
/// upper-case bases the leading bases of USED_REF on the transcript's strand, counted once
/// per row in the default and the tab format with every other field compared.
#[test]
fn refseq_edited_model_deleted_span_codons_accepts_the_genomic_bases_spliced_over_a_shorter_span() {
    let klass = CODONS_CLASS;
    let extra =
        "IMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA;BAM_EDIT=OK";
    let row = default_row_of(
        "rs1179352929",
        "21:45504519-45504535",
        "-",
        "NM_030582.4",
        "frameshift_variant",
        extra,
        "Transcript",
        "80781",
    )
    .replacen(
        POSITIONS_LEFT_OUT,
        "\t3380-3387\t3371-3378\t1124-1126\tPGP/X\tcCAGGCCCC/c\t-\t",
        1,
    );
    let sibling = default_row_of(
        "rs1179352929",
        "21:45504519-45504535",
        "-",
        "NM_001379500.1",
        "frameshift_variant",
        &extra.replace(";BAM_EDIT=OK", ""),
        "Transcript",
        "80781",
    )
    .replacen(
        POSITIONS_LEFT_OUT,
        "\t2897-2913\t2831-2847\t944-949\tPGPPGP/X\tcCAGGCCCCCCAGGCCCA/c\t-\t",
        1,
    );
    let reference = default_text(&[row.as_str(), &sibling]);
    let candidate = with_candidate_codons(&reference);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "Codons"), one(klass));
    assert_eq!(number(&report, &["counts", "compared"]), 2);
    assert_eq!(
        texts(at(&report, &["mask_examples", klass])),
        ["21:45504519-45504535|-|NM_030582.4|Transcript: Codons expected 'cCAGGCCCC/c', got 'cCAGGCCCA/c'"]
    );
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    assert_eq!(differing(&without, "Codons"), 1);
    // The reference's upper-case bases must be the allele's leading bases, so a swap of the
    // two engines is red, as is a reference whose upper-case bases are neither engine's.
    assert_eq!(verdict_under(&candidate, &reference, &[klass]), "red");
    for codons in ["\tcTTTTTTTT/c\t", "\tcCAGgCCCC/c\t"] {
        assert_eq!(
            verdict_under(
                &reference.replace(REFERENCE_CODONS, codons),
                &candidate,
                &[klass]
            ),
            "red",
            "{codons}"
        );
    }
    // A lower-case base inside the run on both sides passes the agreement check and fails
    // the prefix test.
    assert_eq!(
        verdict_under(
            &reference.replace(REFERENCE_CODONS, "\tcCAGgCCCC/c\t"),
            &candidate.replace(CANDIDATE_CODONS, "\tcCAGgCCCA/c\t"),
            &[klass]
        ),
        "red"
    );
    for codons in [
        "cCAGGCCCC/cc",
        "cCAGGCCC/c",
        "cCAGGCCCCC/c",
        "gCAGGCCCA/c",
        "CCAGGCCCA/c",
        "-",
    ] {
        assert_eq!(
            verdict_under(
                &reference,
                &reference.replace(REFERENCE_CODONS, &format!("\t{codons}\t")),
                &[klass]
            ),
            "red",
            "{codons}"
        );
    }
    // On a reverse-strand model the allele is read as its reverse complement.
    let reverse = reference.replace(
        "STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA",
        "STRAND=-1;GIVEN_REF=TGGGCCTGGGGGGCCTG;USED_REF=TGGGCCTGGGGGGCCTG",
    );
    assert_eq!(
        verdict_under(&reverse, &with_candidate_codons(&reverse), &[klass]),
        "green"
    );
    let forward_allele = reverse.replace("STRAND=-1;", "STRAND=1;");
    assert_eq!(
        verdict_under(
            &forward_allele,
            &with_candidate_codons(&forward_allele),
            &[klass]
        ),
        "red"
    );
    let structural = reference.replace("\t-\t80781\tNM_030582.4", "\tdeletion\t80781\tNM_030582.4");
    assert_eq!(
        verdict_under(&structural, &with_candidate_codons(&structural), &[klass]),
        "red"
    );
    assert_eq!(
        verdict_under(&reference, &candidate.replace(";BAM_EDIT=OK", ""), &[klass]),
        "red"
    );
    assert_eq!(
        verdict_under(
            &reference.replace(";BAM_EDIT=OK", ""),
            &candidate.replace(";BAM_EDIT=OK", ""),
            &[klass]
        ),
        "red"
    );
    // An allele of the span's length is spliced over its own bases: the codons then agree,
    // and a difference is another defect's.
    let same_length = reference.replace(
        "GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA",
        "GIVEN_REF=CAGGCCCA;USED_REF=CAGGCCCA",
    );
    assert_eq!(
        verdict_under(&same_length, &with_candidate_codons(&same_length), &[klass]),
        "red"
    );
    let unspanned = reference.replace("\t3371-3378\t", "\t?-3378\t");
    assert_eq!(
        verdict_under(&unspanned, &with_candidate_codons(&unspanned), &[klass]),
        "red"
    );
    // Every other field of the accepted row stays compared, as do residues that are not the
    // translation of their own side's reference codons.
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replace("\tPGP/X\tcCAGGCCCA/c", "\tPGS/X\tcCAGGCCCA/c"),
            &[klass]
        ),
        "red"
    );
    assert_eq!(
        verdict_under(
            &reference,
            &candidate.replacen(
                "IMPACT=HIGH;STRAND=1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA;BAM_EDIT=OK\n",
                "IMPACT=HIGH;STRAND=-1;GIVEN_REF=CAGGCCCCCCAGGCCCA;USED_REF=CAGGCCCCCCAGGCCCA;BAM_EDIT=OK\n",
                1
            ),
            &[klass]
        ),
        "red"
    );
    // The tab format carries the codons, the span and the alleles as fields of their own.
    let tab_row = "rs1179352929\t21:45504519-45504535\t-\t80781\tNM_030582.4\tTranscript\tframeshift_variant\t3380-3387\t3371-3378\t1124-1126\tPGP/X\t{}\t-\tHIGH\t1\tCAGGCCCCCCAGGCCCA\tOK\n";
    // Under --total_length the positions carry a /total suffix, which the span reads past.
    let totals = tab_row.replace(
        "\t3380-3387\t3371-3378\t1124-1126\t",
        "\t3380-3387/4560\t3371-3378/4551\t1124-1126/1516\t",
    );
    for row in [tab_row, totals.as_str()] {
        let tab = |codons: &str| {
            tab_entries(&format!(
                "{}\n{}",
                DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tUSED_REF\tBAM_EDIT"),
                row.replace("{}", codons)
            ))
        };
        let (report, documented) = compared(&tab("cCAGGCCCC/c"), &tab("cCAGGCCCA/c"), &[klass]);
        assert_eq!(verdict(&report), "green");
        assert_eq!(applied(&documented), one(klass));
    }
}

/// The class counts a row once in the VCF and the JSON format, whose bounds may carry the
/// /total suffix of --total_length.
#[test]
fn refseq_edited_model_deleted_span_codons_counts_a_row_once_in_the_vcf_and_json_formats() {
    let klass = CODONS_CLASS;
    let input = "21\t45504518\trs1179352929\tCCAGGCCCCCCAGGCCCA\tC\t.\t.\t.";
    let inputs = [input.to_string()];
    let vcf = |codons: &str| {
        let text = format!(
            "{}{}CSQ=-|frameshift_variant|HIGH|NM_030582.4|Transcript|3371-3378|{codons}|CAGGCCCCCCAGGCCCA|OK\n",
            csq_header("Allele|Consequence|IMPACT|Feature|Feature_type|CDS_position|Codons|USED_REF|BAM_EDIT"),
            open_info(input)
        );
        vcf_entries(&text, &inputs)
    };
    let (report, documented) = compared(&vcf("cCAGGCCCC/c"), &vcf("cCAGGCCCA/c"), &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    let obj = |codons: &str| {
        let entry = format!(
            "{{\"variant_allele\": \"-\", \"transcript_id\": \"NM_030582.4\", \"gene_id\": \"80781\", \"impact\": \"HIGH\", \"consequence_terms\": [\"frameshift_variant\"], \"bam_edit\": \"OK\", \"strand\": 1, \"cds_start\": 3371, \"cds_end\": 3378, \"codons\": {}, \"amino_acids\": \"PGP/X\", \"used_ref\": \"CAGGCCCCCCAGGCCCA\"}}",
            quoted(codons)
        );
        json_object(
            input,
            "rs1179352929",
            "CAGGCCCCCCAGGCCCA/-",
            "frameshift_variant",
            &entry,
        )
    };
    let reference = json_entries(&obj("cCAGGCCCC/c"), &inputs);
    let candidate = json_entries(&obj("cCAGGCCCA/c"), &inputs);
    let (report, documented) = compared(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    let longer_span = obj("cCAGGCCCA/c").replace("\"cds_end\": 3378", "\"cds_end\": 3387");
    let (report, _) = compared(&reference, &json_entries(&longer_span, &inputs), &[klass]);
    assert_eq!(verdict(&report), "red");
    // Under --total_length the JSON bounds are strings with a /total suffix on each.
    let with_totals = |text: String| {
        text.replace(
            "\"cds_start\": 3371, \"cds_end\": 3378",
            "\"cds_start\": \"3371/4551\", \"cds_end\": \"3378/4551\"",
        )
    };
    let (report, documented) = compared(
        &json_entries(&with_totals(obj("cCAGGCCCC/c")), &inputs),
        &json_entries(&with_totals(obj("cCAGGCCCA/c")), &inputs),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
}

/// The Amino_acids and Codons cells of an eighteen-base deletion over the nine bases a model
/// lacks, three whole codons of the model: the reference's read from the allele's leading
/// bases, the candidate's from the model's own.
const WHOLE_CODONS_REFERENCE_CELLS: &str = "\tGPP/-\tGGCCCCCCA/-\t";
const WHOLE_CODONS_CANDIDATE_CELLS: &str = "\tGPR/-\tGGCCCACGT/-\t";

/// The reference's default rows of that deletion: the row of the model that lacks the nine
/// bases, written with `terms`, and the row of a model both engines write alike.
fn whole_codons_reference(terms: &str) -> String {
    let row = default_row_of(
        "1338712",
        "21:45504521-45504538",
        "-",
        "NM_030582.4",
        terms,
        "IMPACT=MODERATE;STRAND=1;GIVEN_REF=GGCCCCCCAGGCCCACGT;USED_REF=GGCCCCCCAGGCCCACGT;BAM_EDIT=OK",
        "Transcript",
        "80781",
    )
    .replacen(
        POSITIONS_LEFT_OUT,
        "\t3382-3390\t3373-3381\t1125-1127\tGPP/-\tGGCCCCCCA/-\t-\t",
        1,
    );
    let alike = default_row(
        "1338712",
        "21:45504521-45504538",
        "-",
        "NM_001379500.1",
        "inframe_deletion",
        "IMPACT=MODERATE;STRAND=1",
    );
    default_text(&[row, alike])
}

/// On a row whose Codons the class takes, it takes the Amino_acids as well where the
/// alternate halves are equal and each reference half is the translation of its own side's
/// reference codons; the row is counted once and any other residue stays a difference.
#[test]
fn refseq_edited_model_deleted_span_codons_accepts_the_residues_translated_from_the_spliced_codons()
{
    let klass = CODONS_CLASS;
    let reference = whole_codons_reference("inframe_deletion");
    let candidate = reference.replace(WHOLE_CODONS_REFERENCE_CELLS, WHOLE_CODONS_CANDIDATE_CELLS);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(number(&report, &["counts", "compared"]), 2);
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
        assert_eq!(differing(&report, column), 0, "{column}");
        assert_eq!(differing(&without, column), 1, "{column}");
    }
    // The two engines swapped are another difference.
    assert_eq!(verdict_under(&candidate, &reference, &[klass]), "red");
    // Residues that are not the translation of their own side's reference codons, on either
    // side, a differing alternate half and a reference half of another length stay
    // differences beside the Codons the class takes.
    for (reference_residues, candidate_residues) in [
        ("GPP/-", "GPK/-"),
        ("GPL/-", "GPR/-"),
        ("GPP/-", "GPR/X"),
        ("GPP/-", "GP/-"),
    ] {
        let (report, documented) = compare_default(
            &reference.replace("GPP/-", reference_residues),
            &candidate.replace("GPR/-", candidate_residues),
            &[klass],
        );
        assert_eq!(verdict(&report), "red", "{candidate_residues}");
        assert_eq!(applied(&documented), one(klass), "{candidate_residues}");
        assert_eq!(differing(&report, "Amino_acids"), 1, "{candidate_residues}");
        assert_eq!(
            masked(&report, "Codons"),
            one(klass),
            "{candidate_residues}"
        );
    }
    // Beside Codons that agree the class takes no residue.
    let (report, documented) =
        compare_default(&reference, &reference.replace("GPP/-", "GPR/-"), &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    assert_eq!(differing(&report, "Amino_acids"), 1);
}

/// The reference allele a text's first row gives in GIVEN_REF.
fn given_ref_of(text: &str) -> &str {
    let rest = text
        .split_once("GIVEN_REF=")
        .expect("a row with a GIVEN_REF")
        .1;
    rest.split_once(';').expect("a field after GIVEN_REF").0
}

/// `allele` with its first `shift` bases moved to its end: a deletion's bases once a 3'
/// shift has moved it `shift` bases along the forward strand.
fn shifted_bases(allele: &str, shift: usize) -> String {
    format!("{}{}", &allele[shift..], &allele[..shift])
}

/// `allele` read on the other strand.
fn other_strand(allele: &str) -> String {
    allele
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => panic!("{other} is no base"),
        })
        .collect()
}

/// A text whose row gives one allele in GIVEN_REF and USED_REF, with the cells as Ensembl
/// VEP prints them where the protein notation moved the deletion `shift` bases along the
/// forward strand: USED_REF the shifted bases, GIVEN_REF the allele as the record wrote it
/// and HGVS_OFFSET the shift.
fn with_used_ref_shifted(text: &str, shift: usize) -> String {
    let allele = given_ref_of(text);
    text.replace(
        &format!("GIVEN_REF={allele};USED_REF={allele}"),
        &format!(
            "HGVS_OFFSET={shift};GIVEN_REF={allele};USED_REF={}",
            shifted_bases(allele, shift)
        ),
    )
}

/// Ensembl VEP reads the Codons and the Amino_acids before the HGVS fields
/// (OutputFactory.pm:1742-1743, :1760, :1769), and the protein notation of a deletion it
/// shifts writes the shifted bases over the reference allele
/// (TranscriptVariationAllele.pm:1785-1786), which USED_REF then prints
/// (OutputFactory.pm:1785). On such a row the reference's upper-case bases are the leading
/// bases of GIVEN_REF, the allele the splice wrote, and the class takes the two cells as it
/// takes them on the row of a run that writes no notation.
#[test]
fn refseq_edited_model_deleted_span_codons_reads_the_given_allele_beside_a_shifted_used_ref() {
    let klass = CODONS_CLASS;
    let unshifted = whole_codons_reference("inframe_deletion");
    let allele = given_ref_of(&unshifted).to_string();
    let as_candidate =
        |text: &str| text.replace(WHOLE_CODONS_REFERENCE_CELLS, WHOLE_CODONS_CANDIDATE_CELLS);
    for shift in [1, 3] {
        let reference = with_used_ref_shifted(&unshifted, shift);
        assert_ne!(reference, unshifted, "{shift}");
        let candidate = as_candidate(&reference);
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "green", "{shift}");
        assert_eq!(applied(&documented), one(klass), "{shift}");
        assert_eq!(number(&report, &["counts", "compared"]), 2, "{shift}");
        let (without, _) = compare_default(&reference, &candidate, &[]);
        assert_eq!(verdict(&without), "red", "{shift}");
        for column in ["Codons", "Amino_acids"] {
            assert_eq!(masked(&report, column), one(klass), "{column} {shift}");
            assert_eq!(differing(&report, column), 0, "{column} {shift}");
            assert_eq!(differing(&without, column), 1, "{column} {shift}");
        }
        // The two engines swapped are another difference.
        assert_eq!(
            verdict_under(&candidate, &reference, &[klass]),
            "red",
            "{shift}"
        );
    }
    let reference = with_used_ref_shifted(&unshifted, 3);
    let cells = format!("GIVEN_REF={allele};USED_REF={}", shifted_bases(&allele, 3));
    assert!(reference.contains(&cells));
    // Upper-case bases that are the leading bases of neither cell stay a difference: a
    // GIVEN_REF that holds the shifted bases too, a row with no GIVEN_REF, and a GIVEN_REF
    // of another length than USED_REF.
    for changed in [
        format!("GIVEN_REF={0};USED_REF={0}", shifted_bases(&allele, 3)),
        format!("USED_REF={}", shifted_bases(&allele, 3)),
        format!(
            "GIVEN_REF={};USED_REF={}",
            &allele[..allele.len() - 1],
            shifted_bases(&allele, 3)
        ),
    ] {
        let reference = reference.replace(&cells, &changed);
        let (report, documented) = compare_default(&reference, &as_candidate(&reference), &[klass]);
        assert_eq!(verdict(&report), "red", "{changed}");
        assert_eq!(applied(&documented), none(), "{changed}");
        for column in ["Codons", "Amino_acids"] {
            assert_eq!(differing(&report, column), 1, "{column} {changed}");
        }
    }
    // On a reverse-strand model either cell is read as its reverse complement.
    let reverse = reference.replace(
        &format!("STRAND=1;HGVS_OFFSET=3;{cells}"),
        &format!(
            "STRAND=-1;HGVS_OFFSET=-3;GIVEN_REF={};USED_REF={}",
            other_strand(&allele),
            shifted_bases(&other_strand(&allele), 3)
        ),
    );
    assert_ne!(reverse, reference);
    assert_eq!(
        verdict_under(&reverse, &as_candidate(&reverse), &[klass]),
        "green"
    );
    let forward_allele = reverse.replace("STRAND=-1;HGVS_OFFSET", "STRAND=1;HGVS_OFFSET");
    assert_eq!(
        verdict_under(&forward_allele, &as_candidate(&forward_allele), &[klass]),
        "red"
    );
    // The tab format carries the two alleles as columns of their own.
    let tab = |columns: &str, alleles: &str, cells: &str| {
        tab_entries(&format!(
            "{}\n1338712\t21:45504521-45504538\t-\t80781\tNM_030582.4\tTranscript\tinframe_deletion\t3382-3390\t3373-3381\t1125-1127{cells}-\tMODERATE\t1\t{alleles}\tOK\n",
            DEFAULT_HEADER.replace("\tExtra", &format!("\tIMPACT\tSTRAND\t{columns}\tBAM_EDIT")),
        ))
    };
    let shifted = shifted_bases(&allele, 3);
    for (columns, alleles, verdict_of_the_pair) in [
        (
            "GIVEN_REF\tUSED_REF",
            format!("{allele}\t{shifted}"),
            "green",
        ),
        ("USED_REF", shifted.clone(), "red"),
    ] {
        let (report, _) = compared(
            &tab(columns, &alleles, WHOLE_CODONS_REFERENCE_CELLS),
            &tab(columns, &alleles, WHOLE_CODONS_CANDIDATE_CELLS),
            &[klass],
        );
        assert_eq!(verdict(&report), verdict_of_the_pair, "{columns}");
    }
    // The VCF and the JSON format give the allele under the names those formats have for it.
    let input = format!("21\t45504520\t1338712\tA{allele}\tA\t.\t.\t.");
    let inputs = [input.clone()];
    let vcf = |residues: &str, codons: &str| {
        let text = format!(
            "{}{}CSQ=-|inframe_deletion|MODERATE|NM_030582.4|Transcript|3373-3381|{residues}|{codons}|1|{allele}|{shifted}|OK\n",
            csq_header("Allele|Consequence|IMPACT|Feature|Feature_type|CDS_position|Amino_acids|Codons|STRAND|GIVEN_REF|USED_REF|BAM_EDIT"),
            open_info(&input)
        );
        vcf_entries(&text, &inputs)
    };
    let obj = |residues: &str, codons: &str| {
        let entry = format!(
            "{{\"variant_allele\": \"-\", \"transcript_id\": \"NM_030582.4\", \"gene_id\": \"80781\", \"impact\": \"MODERATE\", \"consequence_terms\": [\"inframe_deletion\"], \"bam_edit\": \"OK\", \"strand\": 1, \"cds_start\": 3373, \"cds_end\": 3381, \"codons\": {}, \"amino_acids\": {}, \"given_ref\": {}, \"used_ref\": {}}}",
            quoted(codons),
            quoted(residues),
            quoted(&allele),
            quoted(&shifted)
        );
        let text = json_object(
            &input,
            "1338712",
            &format!("{allele}/-"),
            "inframe_deletion",
            &entry,
        );
        json_entries(&text, &inputs)
    };
    let cells_of = |cells: &'static str| -> (&'static str, &'static str) {
        let mut cells = cells.trim_matches('\t').split('\t');
        (
            cells.next().expect("the residues"),
            cells.next().expect("the codons"),
        )
    };
    let (reference_residues, reference_codons) = cells_of(WHOLE_CODONS_REFERENCE_CELLS);
    let (candidate_residues, candidate_codons) = cells_of(WHOLE_CODONS_CANDIDATE_CELLS);
    for (reference, candidate, columns) in [
        (
            vcf(reference_residues, reference_codons),
            vcf(candidate_residues, candidate_codons),
            ["Codons", "Amino_acids"],
        ),
        (
            obj(reference_residues, reference_codons),
            obj(candidate_residues, candidate_codons),
            ["codons", "amino_acids"],
        ),
    ] {
        let (report, documented) = compared(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "green", "{columns:?}");
        assert_eq!(applied(&documented), one(klass), "{columns:?}");
        for column in columns {
            assert_eq!(masked(&report, column), one(klass), "{column}");
        }
    }
    // The class of a span that holds a whole intron reads its Codons by the same test.
    let intron_span = with_used_ref_shifted(&default_text(&intron_span_reference_rows()), 2);
    let (report, documented) = compare_default(
        &intron_span,
        &intron_span_candidate(&intron_span),
        &[INTRON_SPAN_CLASS],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(INTRON_SPAN_CLASS));
}

/// The --masks help names the two cells the class reads the spliced allele from.
#[test]
fn the_masks_help_names_the_cells_the_deleted_span_class_reads_the_allele_from() {
    let help = include_str!("../data/help.txt");
    let help: Vec<&str> = help.split_whitespace().collect();
    assert!(help.join(" ").contains(
        "the reference's upper-case bases the leading bases, on the transcript's strand, of \
         USED_REF or of a GIVEN_REF of its length;"
    ));
}

const PROTEIN_ALTERING: &str = "\tprotein_altering_variant\t";
const INFRAME_DELETION: &str = "\tinframe_deletion\t";
/// The Amino_acids and Codons cells of a twelve-base deletion over the six bases a model
/// lacks, six coding bases inside three codons of the model.
const INSIDE_CODONS_REFERENCE_CELLS: &str = "\tPPP/P\tcCACCGCCG/ccc\t";
const INSIDE_CODONS_CANDIDATE_CELLS: &str = "\tPPS/P\tcCGCCGTCC/ccc\t";

/// The reference's default row of that deletion: the allele's leading bases stand in the
/// reference codons, which neither begin nor end with the alternate codon and do not trim to
/// whole codons against it.
fn inside_codons_reference() -> String {
    let row = default_row_of(
        "1415666",
        "12:122064784-122064795",
        "-",
        "NM_032790.4",
        "protein_altering_variant",
        "IMPACT=MODERATE;STRAND=1;GIVEN_REF=CACCGCCGCCGT;USED_REF=CACCGCCGCCGT;BAM_EDIT=OK",
        "Transcript",
        "84876",
    )
    .replacen(
        POSITIONS_LEFT_OUT,
        "\t350-355\t137-142\t46-48\tPPP/P\tcCACCGCCG/ccc\t-\t",
        1,
    );
    default_text(&[row])
}

/// The candidate's text of such a reference: the term, the residues and the codons read from
/// the model's own bases.
fn inside_codons_candidate(reference: &str) -> String {
    reference
        .replace(PROTEIN_ALTERING, INFRAME_DELETION)
        .replace(INSIDE_CODONS_REFERENCE_CELLS, INSIDE_CODONS_CANDIDATE_CELLS)
}

/// On a row whose Codons the class takes, it takes a Consequence the reference writes with
/// `protein_altering_variant` where the candidate writes `inframe_deletion`, every other term
/// shared, where the in-frame test's two readings of the codons refuse the reference's and
/// accept the candidate's.
#[test]
fn refseq_edited_model_deleted_span_codons_accepts_the_term_the_spliced_codons_fail() {
    let klass = CODONS_CLASS;
    let reference = inside_codons_reference();
    let candidate = inside_codons_candidate(&reference);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    for column in ["Consequence", "Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
        assert_eq!(differing(&report, column), 0, "{column}");
        assert_eq!(differing(&without, column), 1, "{column}");
    }
    assert_eq!(verdict_under(&candidate, &reference, &[klass]), "red");
    // Any other pair of term sets stays a difference beside the two cells the class takes.
    for (reference_terms, candidate_terms) in [
        (
            "\tprotein_altering_variant,splice_region_variant\t",
            INFRAME_DELETION,
        ),
        ("\tmissense_variant\t", INFRAME_DELETION),
        (
            PROTEIN_ALTERING,
            "\tinframe_deletion,splice_region_variant\t",
        ),
        (PROTEIN_ALTERING, "\tframeshift_variant\t"),
    ] {
        let (report, documented) = compare_default(
            &reference.replace(PROTEIN_ALTERING, reference_terms),
            &candidate.replace(INFRAME_DELETION, candidate_terms),
            &[klass],
        );
        assert_eq!(verdict(&report), "red", "{reference_terms}");
        assert_eq!(applied(&documented), one(klass), "{reference_terms}");
        assert_eq!(differing(&report, "Consequence"), 1, "{reference_terms}");
        for column in ["Codons", "Amino_acids"] {
            assert_eq!(masked(&report, column), one(klass), "{reference_terms}");
        }
    }
    // The term pair beside reference codons the in-frame test accepts is another difference:
    // a deletion of whole codons leaves no alternate codon, which every reference begins with.
    let whole = whole_codons_reference("protein_altering_variant");
    let (report, documented) = compare_default(
        &whole,
        &whole
            .replacen(PROTEIN_ALTERING, INFRAME_DELETION, 1)
            .replace(WHOLE_CODONS_REFERENCE_CELLS, WHOLE_CODONS_CANDIDATE_CELLS),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(differing(&report, "Consequence"), 1);
    // Beside Codons that agree the class takes no term, and every other field of the
    // accepted row stays compared.
    let (report, documented) = compare_default(
        &reference,
        &reference.replace(PROTEIN_ALTERING, INFRAME_DELETION),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    for (written, changed) in [
        ("IMPACT=MODERATE;", "IMPACT=HIGH;"),
        ("\t350-355\t", "\t350-356\t"),
        ("\t46-48\t", "\t46-49\t"),
    ] {
        assert_eq!(
            verdict_under(&reference, &candidate.replace(written, changed), &[klass]),
            "red",
            "{changed}"
        );
    }
    // The tab format carries the fields of the Extra text as columns of their own.
    let tab = |text: &str| {
        let header = DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tUSED_REF\tBAM_EDIT");
        let row = text.lines().nth(1).expect("the row of the class").replace(
            "IMPACT=MODERATE;STRAND=1;GIVEN_REF=CACCGCCGCCGT;USED_REF=CACCGCCGCCGT;BAM_EDIT=OK",
            "MODERATE\t1\tCACCGCCGCCGT\tOK",
        );
        tab_entries(&format!("{header}\n{row}\n"))
    };
    let (report, documented) = compared(&tab(&reference), &tab(&candidate), &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    for column in ["Consequence", "Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
    }
}

/// The class counts such a row once in the VCF and the JSON format, under the names those
/// formats give the three cells, and leaves uncompared the most_severe_consequence of a JSON
/// record that carries the term pair.
#[test]
fn refseq_edited_model_deleted_span_codons_counts_the_three_cells_once_in_the_vcf_and_json_formats()
{
    let klass = CODONS_CLASS;
    let input = "12\t122064783\t1415666\tGCACCGCCGCCGT\tG\t.\t.\t.";
    let inputs = [input.to_string()];
    let vcf = |term: &str, amino_acids: &str, codons: &str| {
        let text = format!(
            "{}{}CSQ=-|{term}|MODERATE|NM_032790.4|Transcript|137-142|{amino_acids}|{codons}|1|CACCGCCGCCGT|OK\n",
            csq_header("Allele|Consequence|IMPACT|Feature|Feature_type|CDS_position|Amino_acids|Codons|STRAND|USED_REF|BAM_EDIT"),
            open_info(input)
        );
        vcf_entries(&text, &inputs)
    };
    let reference = vcf("protein_altering_variant", "PPP/P", "cCACCGCCG/ccc");
    let candidate = vcf("inframe_deletion", "PPS/P", "cCGCCGTCC/ccc");
    let (report, documented) = compared(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    for column in ["Consequence", "Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
    }
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    let obj = |term: &str, amino_acids: &str, codons: &str| {
        let entry = format!(
            "{{\"variant_allele\": \"-\", \"transcript_id\": \"NM_032790.4\", \"gene_id\": \"84876\", \"impact\": \"MODERATE\", \"consequence_terms\": [{}], \"bam_edit\": \"OK\", \"strand\": 1, \"cds_start\": 137, \"cds_end\": 142, \"codons\": {}, \"amino_acids\": {}, \"used_ref\": \"CACCGCCGCCGT\"}}",
            quoted(term),
            quoted(codons),
            quoted(amino_acids)
        );
        json_object(input, "1415666", "CACCGCCGCCGT/-", term, &entry)
    };
    let reference_text = obj("protein_altering_variant", "PPP/P", "cCACCGCCG/ccc");
    let candidate_text = obj("inframe_deletion", "PPS/P", "cCGCCGTCC/ccc");
    let reference = json_entries(&reference_text, &inputs);
    let candidate = json_entries(&candidate_text, &inputs);
    let (report, documented) = compared(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    for column in ["consequence_terms", "codons", "amino_acids"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
    }
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    // The record's other keys and the entry's stay compared.
    for (written, changed) in [
        ("\"id\": \"1415666\"", "\"id\": \"rs1\""),
        ("\"impact\": \"MODERATE\"", "\"impact\": \"HIGH\""),
        ("\"cds_end\": 142", "\"cds_end\": 143"),
    ] {
        let text = candidate_text.replace(written, changed);
        let (report, _) = compared(&reference, &json_entries(&text, &inputs), &[klass]);
        assert_eq!(verdict(&report), "red", "{changed}");
    }
    // A record whose rows hold no such term pair keeps its most_severe_consequence compared.
    let residues_alone = obj("inframe_deletion", "PPP/P", "cCACCGCCG/ccc").replace(
        "\"most_severe_consequence\": \"inframe_deletion\"",
        "\"most_severe_consequence\": \"protein_altering_variant\"",
    );
    let (report, documented) = compared(
        &json_entries(&residues_alone, &inputs),
        &candidate,
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), one(klass));
}

const INTRON_SPAN_CLASS: &str = classes::REFSEQ_EDITED_MODEL_INTRON_SPAN_CODING_REFERENCE;
const INTRON_SPAN_REFERENCE_TERMS: &str =
    "splice_acceptor_variant,splice_donor_variant,stop_lost,protein_altering_variant,intron_variant";
const INTRON_SPAN_CANDIDATE_TERMS: &str =
    "splice_acceptor_variant,splice_donor_variant,inframe_deletion,intron_variant";
/// The four cells of the class in the default format, and the Extra text that holds one.
const INTRON_SPAN_COLUMNS: [&str; 5] = ["Consequence", "Codons", "Amino_acids", "HGVSp", "Extra"];

/// The reference's default rows of a made-up record that deletes the coding bases 5 to 7 of
/// a forward-strand edited model, the six bases of an intron between bases 5 and 6 among the
/// nine it removes: the row of that model, whose reference codons are read from the allele,
/// and the row of a model both engines write alike.
fn intron_span_reference_rows() -> Vec<String> {
    let extra = "IMPACT=HIGH;STRAND=1;GIVEN_REF=AGTAAAGTT;USED_REF=AGTAAAGTT;BAM_EDIT=OK;HGVSc=NM_990005.1:c.5_7del;HGVSp=NP_990005.1:p.His2_Ser3delextTer?";
    let row = default_row_of(
        "rec1",
        "7:1005-1013",
        "-",
        "NM_990005.1",
        INTRON_SPAN_REFERENCE_TERMS,
        extra,
        "Transcript",
        "ENSG1",
    )
    .replacen(
        POSITIONS_LEFT_OUT,
        "\t105-107\t5-7\t2-3\tQ*/P\tcAGTAA/cca\t-\t",
        1,
    );
    let alike = default_row(
        "rec1",
        "7:1005-1013",
        "-",
        "NM_990006.1",
        "intron_variant",
        "IMPACT=MODIFIER;STRAND=1",
    );
    vec![row, alike]
}

/// The candidate's text of a reference text: the four cells of the class's row read from the
/// model's own codons.
fn intron_span_candidate(reference: &str) -> String {
    reference
        .replace(INTRON_SPAN_REFERENCE_TERMS, INTRON_SPAN_CANDIDATE_TERMS)
        .replace("Q*/P", "HS/P")
        .replace("cAGTAA/cca", "cATTCA/cca")
        .replace("extTer?", "")
}

/// On a BAM_EDIT OK transcript row whose CDS span is one range on both sides, whose USED_REF
/// has more bases than the range and whose terms hold the three of an intron inside the span,
/// the class sets aside the Consequence, the Codons, the Amino_acids and the HGVSp, counted
/// once per row in the default and the tab format with every other field compared.
#[test]
fn refseq_edited_model_intron_span_coding_reference_accepts_the_four_cells_read_from_the_allele() {
    let klass = INTRON_SPAN_CLASS;
    let reference = default_text(&intron_span_reference_rows());
    let candidate = intron_span_candidate(&reference);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(number(&report, &["counts", "compared"]), 2);
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    for column in INTRON_SPAN_COLUMNS {
        assert_eq!(masked(&report, column), one(klass), "{column}");
        assert_eq!(differing(&report, column), 0, "{column}");
        assert_eq!(differing(&without, column), 1, "{column}");
    }
    // The two engines swapped are another difference.
    assert_eq!(verdict_under(&candidate, &reference, &[klass]), "red");
    // Every other field of the accepted row stays compared.
    for (written, changed) in [
        ("IMPACT=HIGH;", "IMPACT=MODERATE;"),
        ("c.5_7del", "c.5_8del"),
        ("\t105-107\t", "\t105-108\t"),
        ("\t2-3\t", "\t2-4\t"),
    ] {
        assert_eq!(
            verdict_under(&reference, &candidate.replace(written, changed), &[klass]),
            "red",
            "{changed}"
        );
    }
    // A notation that is not the reference's less its extension stays a difference beside the
    // three cells the class takes.
    let other_notation = candidate.replace("Ser3del", "Ser3delinsPro");
    let (report, documented) = compare_default(&reference, &other_notation, &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), one(klass));
    for column in ["HGVSp", "Extra"] {
        assert_eq!(differing(&report, column), 1, "{column}");
    }
    for column in ["Consequence", "Codons", "Amino_acids"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
    }
    // A row without the shape keeps all four cells as differences: a model without an applied
    // edit, a span without the intron term, a coding span with an open end and an allele of
    // the span's length.
    for (written, changed) in [
        (";BAM_EDIT=OK", ""),
        (",intron_variant\t", "\t"),
        ("\t5-7\t", "\t?-7\t"),
        (
            "GIVEN_REF=AGTAAAGTT;USED_REF=AGTAAAGTT",
            "GIVEN_REF=AGT;USED_REF=AGT",
        ),
    ] {
        let (reference, candidate) = (
            reference.replace(written, changed),
            candidate.replace(written, changed),
        );
        let (report, documented) = compare_default(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "red", "{changed}");
        assert_eq!(applied(&documented), none(), "{changed}");
        for column in INTRON_SPAN_COLUMNS {
            assert_eq!(differing(&report, column), 1, "{changed} {column}");
        }
    }
    // The tab format carries the fields of the Extra text as columns of their own.
    let tab = |text: &str| {
        let header =
            DEFAULT_HEADER.replace("\tExtra", "\tIMPACT\tSTRAND\tUSED_REF\tBAM_EDIT\tHGVSp");
        let row = text
            .lines()
            .nth(1)
            .expect("the row of the class")
            .replace(
                "IMPACT=HIGH;STRAND=1;GIVEN_REF=AGTAAAGTT;USED_REF=AGTAAAGTT;BAM_EDIT=OK;HGVSc=NM_990005.1:c.5_7del;HGVSp=",
                "HIGH\t1\tAGTAAAGTT\tOK\t",
            );
        tab_entries(&format!("{header}\n{row}\n"))
    };
    let (report, documented) = compared(&tab(&reference), &tab(&candidate), &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    for column in ["Consequence", "Codons", "Amino_acids", "HGVSp"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
    }
    assert_eq!(
        verdict(&compare(
            &tab(&reference),
            &tab(&candidate),
            &mut no_documented()
        )),
        "red"
    );
}

/// The class counts a row once in the VCF and the JSON format, under the names those formats
/// give the four cells, and leaves the JSON record's most_severe_consequence compared.
#[test]
fn refseq_edited_model_intron_span_coding_reference_counts_a_row_once_in_the_vcf_and_json_formats()
{
    let klass = INTRON_SPAN_CLASS;
    let input = "7\t1004\trec1\tGAGTAAAGTT\tG\t.\t.\t.";
    let inputs = [input.to_string()];
    let notation = "NP_990005.1:p.His2_Ser3del";
    let vcf = |terms: &str, amino_acids: &str, codons: &str, hgvsp: &str| {
        let text = format!(
            "{}{}CSQ=-|{}|HIGH|NM_990005.1|Transcript|{hgvsp}|5-7|{amino_acids}|{codons}|1|AGTAAAGTT|OK\n",
            csq_header("Allele|Consequence|IMPACT|Feature|Feature_type|HGVSp|CDS_position|Amino_acids|Codons|STRAND|USED_REF|BAM_EDIT"),
            open_info(input),
            terms.replace(',', "&")
        );
        vcf_entries(&text, &inputs)
    };
    let reference = vcf(
        INTRON_SPAN_REFERENCE_TERMS,
        "Q*/P",
        "cAGTAA/cca",
        &format!("{notation}extTer?"),
    );
    let candidate = vcf(INTRON_SPAN_CANDIDATE_TERMS, "HS/P", "cATTCA/cca", notation);
    let (report, documented) = compared(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    for column in ["Consequence", "Codons", "Amino_acids", "HGVSp"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
    }
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    let obj = |terms: &str, amino_acids: &str, codons: &str, hgvsp: &str| {
        let term_list: Vec<&str> = terms.split(',').collect();
        let entry = format!(
            "{{\"variant_allele\": \"-\", \"transcript_id\": \"NM_990005.1\", \"gene_id\": \"ENSG1\", \"impact\": \"HIGH\", \"consequence_terms\": {}, \"bam_edit\": \"OK\", \"strand\": 1, \"cds_start\": 5, \"cds_end\": 7, \"codons\": {}, \"amino_acids\": {}, \"used_ref\": \"AGTAAAGTT\", \"hgvsp\": {}}}",
            quoted_list(&term_list),
            quoted(codons),
            quoted(amino_acids),
            quoted(hgvsp)
        );
        json_object(
            input,
            "rec1",
            "AGTAAAGTT/-",
            "splice_acceptor_variant",
            &entry,
        )
    };
    let reference_text = obj(
        INTRON_SPAN_REFERENCE_TERMS,
        "Q*/P",
        "cAGTAA/cca",
        &format!("{notation}extTer?"),
    );
    let candidate_text = obj(INTRON_SPAN_CANDIDATE_TERMS, "HS/P", "cATTCA/cca", notation);
    let reference = json_entries(&reference_text, &inputs);
    let candidate = json_entries(&candidate_text, &inputs);
    let (report, documented) = compared(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    for column in ["consequence_terms", "codons", "amino_acids", "hgvsp"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
    }
    assert_eq!(
        verdict(&compare(&reference, &candidate, &mut no_documented())),
        "red"
    );
    // The record's own keys and the entry's other keys stay compared, and a coding span of
    // one position is no range.
    for (written, changed) in [
        (
            "\"most_severe_consequence\": \"splice_acceptor_variant\"",
            "\"most_severe_consequence\": \"splice_donor_variant\"",
        ),
        ("\"impact\": \"HIGH\"", "\"impact\": \"MODERATE\""),
        ("\"cds_end\": 7", "\"cds_end\": 8"),
    ] {
        let text = candidate_text.replace(written, changed);
        let (report, _) = compared(&reference, &json_entries(&text, &inputs), &[klass]);
        assert_eq!(verdict(&report), "red", "{changed}");
    }
    let one_position = |text: &str| text.replace("\"cds_end\": 7", "\"cds_end\": 5");
    let (report, documented) = compared(
        &json_entries(&one_position(&reference_text), &inputs),
        &json_entries(&one_position(&candidate_text), &inputs),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
}

/// A comparison that names the class and `refseq_edited_model_deleted_span_codons` counts the
/// Codons of a row with an intron inside its span under the class, which is asked first, and
/// the Codons of a row without the intron terms under the other. The other class named alone
/// takes the Codons of both rows and the Amino_acids beside them and leaves the class's term
/// pair and notation as differences; the class named alone leaves the other's Codons as one.
#[test]
fn refseq_edited_model_intron_span_coding_reference_is_asked_before_the_deleted_span_class() {
    let klass = INTRON_SPAN_CLASS;
    let deleted_span = default_row_of(
        "rec2",
        "7:2008-2013",
        "-",
        "NM_990007.1",
        "frameshift_variant",
        "IMPACT=HIGH;STRAND=1;GIVEN_REF=ACGTAC;USED_REF=ACGTAC;BAM_EDIT=OK",
        "Transcript",
        "ENSG2",
    )
    .replacen(
        POSITIONS_LEFT_OUT,
        "\t108-110\t8-10\t3-4\tTX/X\tcACGT/c\t-\t",
        1,
    );
    let mut rows = intron_span_reference_rows();
    rows.push(deleted_span);
    let reference = default_text(&rows);
    let candidate = intron_span_candidate(&reference).replace("cACGT/c", "cACTA/c");
    let (report, documented) = compare_default(&reference, &candidate, &[klass, CODONS_CLASS]);
    assert_eq!(verdict(&report), "green");
    let each = count_map(&[(klass, 1), (CODONS_CLASS, 1)]);
    assert_eq!(applied(&documented), each);
    assert_eq!(masked(&report, "Codons"), each);
    for column in ["Consequence", "Amino_acids", "HGVSp", "Extra"] {
        assert_eq!(masked(&report, column), one(klass), "{column}");
    }
    // The other class's own rule reads no intron term, so alone it takes the Codons of both
    // rows and the Amino_acids beside them; the term pair with `stop_lost` is not its own.
    let (report, documented) = compare_default(&reference, &candidate, &[CODONS_CLASS]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), count_map(&[(CODONS_CLASS, 2)]));
    assert_eq!(masked(&report, "Amino_acids"), one(CODONS_CLASS));
    for column in ["Codons", "Amino_acids"] {
        assert_eq!(differing(&report, column), 0, "{column}");
    }
    for column in ["Consequence", "HGVSp", "Extra"] {
        assert_eq!(differing(&report, column), 1, "{column}");
    }
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "Codons"), one(klass));
    assert_eq!(differing(&report, "Codons"), 1);
}

const ZERO_CLASS: &str = classes::FREQUENCY_FILTER_ZERO_VALUE;

/// `ZERO_CLASS` as a comparison names it where Ensembl VEP 116.2 wrote the reference
/// output. Ensembl VEP 116.3 corrects the class's defect, and a comparison under 116.3
/// takes the class out of the named set.
fn zero_class_under_116_2() -> &'static str {
    let named = masks(&[ZERO_CLASS]);
    assert_eq!(
        masks_under(named.clone(), "116.2"),
        (named.clone(), BTreeMap::new())
    );
    let corrected = BTreeMap::from([(ZERO_CLASS.to_string(), "116.3".to_string())]);
    assert_eq!(masks_under(named, "116.3"), (BTreeSet::new(), corrected));
    ZERO_CLASS
}

fn freqs_row(
    name: &str,
    location: &str,
    allele: &str,
    feature: &str,
    terms: &str,
    frequency: &str,
) -> String {
    default_row(
        name,
        location,
        allele,
        feature,
        terms,
        &format!("IMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE_AFR:{allele}:{frequency}"),
    )
}

/// A FREQS the reference reads NA and the candidate 0, item for item, is the class.
#[test]
fn frequency_filter_zero_value_accepts_a_zero_read_as_na() {
    let klass = zero_class_under_116_2();
    let reference = default_text(&[freqs_row(
        "a",
        "21:5131441",
        "C",
        "ENST1",
        "intron_variant",
        "NA",
    )]);
    let candidate = reference.replace("GNOMADE_AFR:C:NA", "GNOMADE_AFR:C:0");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(verdict_under(&reference, &candidate, &[]), "red");
    assert_eq!(
        verdict_under(
            &reference,
            &reference.replace("GNOMADE_AFR:C:NA", "GNOMADE_AFR:C:0.01"),
            &[klass]
        ),
        "red"
    );
    let two = reference.replace("FREQS=GNOMADE_AFR:C:NA", "FREQS=GNOMADE_AFR:C:NA,AF:C:0.2");
    assert_eq!(
        verdict_under(
            &two,
            &two.replace("GNOMADE_AFR:C:NA", "GNOMADE_AFR:C:0"),
            &[klass]
        ),
        "green"
    );
    assert_eq!(
        verdict_under(
            &two,
            &two.replace("GNOMADE_AFR:C:NA,AF:C:0.2", "GNOMADE_AFR:C:0,AF:C:0.3"),
            &[klass]
        ),
        "red"
    );
}

/// The rows of a record the candidate alone keeps whose every FREQS item is 0 are the
/// class; a record kept on a non-zero frequency stays unexpected, and so does an extra row
/// on a record the reference wrote a row for, whatever its FREQS says.
#[test]
fn frequency_filter_zero_value_accepts_the_rows_of_a_record_kept_for_a_zero_frequency() {
    let klass = zero_class_under_116_2();
    let a = freqs_row("a", "21:5131441", "C", "ENST1", "intron_variant", "0");
    let kept = default_text(&[
        a.clone(),
        freqs_row("b", "21:29025626", "G", "ENST2", "intron_variant", "0"),
        freqs_row(
            "b",
            "21:29025626",
            "G",
            "ENST3",
            "upstream_gene_variant",
            "0",
        ),
    ]);
    let only_a = default_text(&[a.as_str()]);
    let (report, documented) = compare_default(&only_a, &kept, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(klass, 2)]));
    let (without, _) = compare_default(&only_a, &kept, &[]);
    assert_eq!(number(&without, &["one_sided", "actual_only"]), 2);
    let zero = "ENST3\tTranscript\tupstream_gene_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;FREQS=GNOMADE_AFR:G:0";
    assert!(kept.contains(zero));
    let mixed = kept.replace(zero, &format!("{zero}.3"));
    assert_eq!(verdict_under(&only_a, &mixed, &[klass]), "red");
    let both_a = default_text(&[
        a.clone(),
        freqs_row(
            "a",
            "21:5131441",
            "C",
            "ENST4",
            "upstream_gene_variant",
            "0",
        ),
    ]);
    let (report, documented) = compare_default(&only_a, &both_a, &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["one_sided", "actual_only"]), 1);
    assert_eq!(applied(&documented), none());
}

/// A JSON consequence entry carries no FREQS, so the built-in row rule finds nothing there
/// and the kept record's entries are unexpected; a manifest's documented tuples of the class
/// mark them as counted, not compared.
#[test]
fn frequency_filter_zero_value_reads_a_kept_json_record_from_the_manifest() {
    let klass = zero_class_under_116_2();
    let inputs = [
        "21\t5131441\t.\tT\tC\t.\t.\t.".to_string(),
        "21\t29025626\t.\tA\tG\t.\t.\t.".to_string(),
    ];
    let obj = |i: usize, allele: &str| {
        format!(
            "{{\"input\": {}, \"id\": \"r{i}\", \"most_severe_consequence\": \"intron_variant\", \"colocated_variants\": [{{\"frequencies\": {{{}: {{\"gnomade_afr\": 0}}}}}}], \"transcript_consequences\": [{{\"variant_allele\": {}, \"transcript_id\": \"ENST{i}\", \"consequence_terms\": [\"intron_variant\"]}}]}}\n",
            quoted(&inputs[i]),
            quoted(allele),
            quoted(allele)
        )
    };
    let reference = json_entries(&obj(0, "C"), &inputs);
    let candidate = json_entries(&format!("{}{}", obj(0, "C"), obj(1, "G")), &inputs);
    let (report, documented) = compared(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    assert_eq!(number(&report, &["one_sided", "actual_only"]), 2);
    let manifest = |klass: &str| {
        json(&format!(
            "{{\"vep_rs_only_tuples\": [{{\"location\": \"21:29025626\", \"allele\": \"G\", \"feature\": \"ENST1\", \"record_indices\": [1], \"expected_divergence\": {}}}]}}",
            quoted(klass)
        ))
    };
    let mut documented = from_manifest(&manifest(klass), &[klass]);
    assert!(documented.reference_skipped(1));
    assert!(!documented.reference_skipped(0));
    let report = compare(&reference, &candidate, &mut documented);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(
        number(
            &report,
            &["counts", "record_entries_on_reference_skipped_records"]
        ),
        1
    );
    let documented = from_manifest(&manifest("other"), &["other"]);
    assert!(!documented.reference_skipped(1));
}

/// Under Ensembl VEP 116.3 the named set holds no `ZERO_CLASS`: the value and the rows the
/// class sets aside under 116.2 are differences, and nothing is applied.
#[test]
fn frequency_filter_zero_value_sets_nothing_aside_under_the_release_that_corrects_it() {
    let compared_under = |reference: &str, candidate: &str, release: &str| {
        let (named, _) = masks_under(masks(&[ZERO_CLASS]), release);
        let named: Vec<&str> = named.iter().map(String::as_str).collect();
        compare_default(reference, candidate, &named)
    };
    let read_as_na = default_text(&[freqs_row(
        "a",
        "21:5131441",
        "C",
        "ENST1",
        "intron_variant",
        "NA",
    )]);
    let read_as_zero = read_as_na.replace("GNOMADE_AFR:C:NA", "GNOMADE_AFR:C:0");
    let (report, documented) = compared_under(&read_as_na, &read_as_zero, "116.2");
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(ZERO_CLASS));
    let (report, documented) = compared_under(&read_as_na, &read_as_zero, "116.3");
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    assert_eq!(differing(&report, "FREQS"), 1);
    let a = freqs_row("a", "21:5131441", "C", "ENST1", "intron_variant", "0");
    let only_a = default_text(&[a.as_str()]);
    let kept = default_text(&[
        a.clone(),
        freqs_row("b", "21:29025626", "G", "ENST2", "intron_variant", "0"),
        freqs_row(
            "b",
            "21:29025626",
            "G",
            "ENST3",
            "upstream_gene_variant",
            "0",
        ),
    ]);
    let (report, documented) = compared_under(&only_a, &kept, "116.2");
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), count_map(&[(ZERO_CLASS, 2)]));
    let (report, documented) = compared_under(&only_a, &kept, "116.3");
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    assert_eq!(number(&report, &["one_sided", "actual_only"]), 2);
}

/// Two input records at one position are two records in the text formats too: the
/// Uploaded_variation tells them apart, so a record the reference wrote rows for keeps its
/// extra candidate row a difference while the other record's rows are set aside.
#[test]
fn the_frequency_filter_row_rule_reads_one_input_record_per_label() {
    let klass = zero_class_under_116_2();
    let row = |name: &str, allele: &str, feature: &str, terms: &str| {
        freqs_row(name, "21:100", allele, feature, terms, "0")
    };
    let reference = default_text(&[row("21_100_A/G", "G", "ENST1", "intron_variant")]);
    // The second record at the position is one the candidate alone keeps: its rows are the
    // class.
    let other_record = format!(
        "{reference}{}\n",
        row("21_100_A/T", "T", "ENST1", "intron_variant")
    );
    let (report, documented) = compare_default(&reference, &other_record, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    // An extra row on the record the reference wrote for stays unexpected.
    let same_record = format!(
        "{reference}{}\n",
        row("21_100_A/G", "G", "ENST2", "upstream_gene_variant")
    );
    let (report, documented) = compare_default(&reference, &same_record, &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    assert_eq!(number(&report, &["one_sided", "actual_only"]), 1);
}

const SWAP_CLASS: &str = classes::SPLICE_FAMILY_SWAP;

/// The default-format rows of a multi-base delins at 12:5000103-5000106 on
/// ENST00000000051: the reference's with the Downstream plugin's two fields, the candidate's
/// with `candidate_fields`.
fn downstream_pair(
    reference_terms: &str,
    candidate_terms: &str,
    candidate_fields: &str,
) -> (String, String) {
    let row = |terms: &str, fields: &str| {
        default_text(&[default_row(
            "delins_swap",
            "12:5000103-5000106",
            "ACACACAC",
            "ENST00000000051",
            terms,
            &format!("IMPACT=HIGH;STRAND=-1{fields}"),
        )])
    };
    (
        row(
            reference_terms,
            ";DownstreamProtein=KL;ProteinLengthChange=-120",
        ),
        row(candidate_terms, candidate_fields),
    )
}

fn first_entries(reference: &str, candidate: &str) -> (Entry, Entry) {
    let first = |text: &str| parse_default(text).expect("the text parses").remove(0);
    (first(reference), first(candidate))
}

/// A multi-base delins whose span touches a splice-region window: the candidate's terms are
/// the reference's plus splice_region_variant, so the Downstream plugin's gate writes nothing
/// in the candidate and the reference alone fills the two fields. The field class accepts
/// each on that key; the comparison itself skips the pair whole and counts it once.
#[test]
fn splice_family_swap_accepts_the_downstream_fields_the_reference_alone_writes() {
    let klass = SWAP_CLASS;
    let (reference, candidate) = downstream_pair(
        "frameshift_variant",
        "frameshift_variant,splice_region_variant",
        "",
    );
    let (e, a) = first_entries(&reference, &candidate);
    for (name, value) in [("DownstreamProtein", "KL"), ("ProteinLengthChange", "-120")] {
        assert!(
            class_accepts(klass, &e, name, value, ABSENT, Some(&a)),
            "{name}"
        );
    }
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(verdict_under(&reference, &candidate, &[]), "red");
}

/// A Downstream field the reference alone writes stays a difference where the two
/// Consequence sets agree or differ by another term; on the class's own key a field both
/// fill, one the candidate alone fills, a field of another plugin, and a call without the
/// candidate's row are refused.
#[test]
fn splice_family_swap_refuses_a_downstream_field_on_any_other_pair_of_rows() {
    let klass = SWAP_CLASS;
    for candidate_terms in [
        "frameshift_variant",
        "frameshift_variant,splice_donor_variant",
    ] {
        let (reference, candidate) = downstream_pair("frameshift_variant", candidate_terms, "");
        let (e, a) = first_entries(&reference, &candidate);
        assert!(
            !class_accepts(klass, &e, "DownstreamProtein", "KL", ABSENT, Some(&a)),
            "{candidate_terms}"
        );
    }
    // The plugin writes on a frameshift row none of whose terms names a splice site: a
    // reference row without the frameshift, or with a splice term of its own, wrote nothing.
    for reference_terms in [
        "intron_variant",
        "frameshift_variant,splice_donor_5th_base_variant",
    ] {
        let (reference, candidate) = downstream_pair(
            reference_terms,
            &format!("{reference_terms},splice_region_variant"),
            "",
        );
        let (e, a) = first_entries(&reference, &candidate);
        assert!(
            !class_accepts(klass, &e, "DownstreamProtein", "KL", ABSENT, Some(&a)),
            "{reference_terms}"
        );
    }
    let (reference, candidate) = downstream_pair("frameshift_variant", "frameshift_variant", "");
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), none());
    assert_eq!(differing(&report, "DownstreamProtein"), 1);
    let (reference, candidate) = downstream_pair(
        "frameshift_variant",
        "frameshift_variant,splice_region_variant",
        "",
    );
    let (e, a) = first_entries(&reference, &candidate);
    assert!(!class_accepts(
        klass,
        &e,
        "DownstreamProtein",
        "KL",
        "KR",
        Some(&a)
    ));
    assert!(!class_accepts(
        klass,
        &a,
        "DownstreamProtein",
        ABSENT,
        "KL",
        Some(&e)
    ));
    assert!(!class_accepts(klass, &e, "REVEL", "0.5", ABSENT, Some(&a)));
    assert!(!class_accepts(
        klass,
        &e,
        "DownstreamProtein",
        "KL",
        ABSENT,
        None
    ));
}

/// The texts of a stop codon row with these terms and this IMPACT, none where it is empty.
#[derive(Clone, Copy)]
struct StopText<'a> {
    cdna: &'a str,
    cds: &'a str,
    allele: &'a str,
    location: &'a str,
}

impl StopText<'_> {
    fn text(&self, terms: &str, impact: &str) -> String {
        let extra = if impact.is_empty() {
            "STRAND=-1".to_string()
        } else {
            format!("IMPACT={impact};STRAND=-1")
        };
        stop_codon_text(
            terms,
            self.cdna,
            self.cds,
            &extra,
            self.allele,
            self.location,
        )
    }
}

/// What every class of the following cells shows through the comparison itself: the pair
/// `taken` on the columns of `at` is green, its Consequence set aside under the class and
/// counted once, its Extra set aside too where the two impacts differ; the same pair on the
/// columns of `elsewhere` stays a difference; and an IMPACT the candidate leaves out or one
/// that is not its terms' (`planted`) stays a difference on a row the class still counts.
fn is_a_built_in_class_of_the_following_cells(
    klass: &str,
    at: StopText,
    elsewhere: StopText,
    taken: [(&str, &str); 2],
    elsewhere_taken: [(&str, &str); 2],
    planted: &str,
) {
    let [(reference_terms, reference_impact), (candidate_terms, candidate_impact)] = taken;
    let reference = at.text(reference_terms, reference_impact);
    let candidate = at.text(candidate_terms, candidate_impact);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green", "{klass}");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "Consequence"), one(klass));
    if reference_impact == candidate_impact {
        assert_eq!(differing(&report, "Extra"), 0);
        assert!(!has_masked(&report, "Extra"));
    } else {
        assert_eq!(masked(&report, "Extra"), one(klass));
    }
    assert_eq!(verdict_under(&reference, &candidate, &[]), "red");
    let [(reference_terms, reference_impact), (candidate_terms, candidate_impact)] =
        elsewhere_taken;
    let (report, documented) = compare_default(
        &elsewhere.text(reference_terms, reference_impact),
        &elsewhere.text(candidate_terms, candidate_impact),
        &[klass],
    );
    assert_eq!(verdict(&report), "red", "{klass}");
    assert_eq!(applied(&documented), none());
    for planted in ["", planted] {
        let (report, documented) =
            compare_default(&reference, &at.text(taken[1].0, planted), &[klass]);
        assert_eq!(verdict(&report), "red", "{klass} {planted}");
        assert_eq!(applied(&documented), one(klass), "{planted}");
        assert_eq!(differing(&report, "Extra"), 1, "{planted}");
        assert!(!has_masked(&report, "Extra"), "{planted}");
    }
}

const DELETION_UNTESTED: StopText = StopText {
    cdna: "?-7057",
    cds: "-",
    allele: "-",
    location: STOP_CODON_LOCATION,
};

/// Through the comparison itself the class sets aside the Consequence of a deletion's row
/// that shows the stop codon test is not made, counted once, and with it the IMPACT that is
/// each side's most severe term's.
#[test]
fn stop_retained_without_stop_codon_test_is_a_built_in_class() {
    let pair = [(RETAINED, "LOW"), (UNKNOWN, "MODIFIER")];
    is_a_built_in_class_of_the_following_cells(
        classes::STOP_RETAINED_WITHOUT_STOP_CODON_TEST,
        DELETION_UNTESTED,
        StopText {
            cdna: "7055-7059",
            cds: "6978-6982",
            ..DELETION_UNTESTED
        },
        pair,
        pair,
        "HIGH",
    );
}

/// Through the comparison itself the class sets aside the Consequence of a deletion's row
/// beside an exon's end, both sides holding splice_region_variant, counted once.
#[test]
fn stop_codon_overlap_genomic_span_is_a_built_in_class() {
    let at = StopText {
        cdna: "1705",
        cds: "1630",
        allele: "-",
        location: "3:36993663",
    };
    is_a_built_in_class_of_the_following_cells(
        classes::STOP_CODON_OVERLAP_GENOMIC_SPAN,
        at,
        at,
        [
            (
                "splice_region_variant,coding_sequence_variant,NMD_transcript_variant",
                "LOW",
            ),
            (
                "splice_region_variant,stop_retained_variant,NMD_transcript_variant",
                "LOW",
            ),
        ],
        [
            ("coding_sequence_variant,NMD_transcript_variant", "MODIFIER"),
            ("stop_retained_variant,NMD_transcript_variant", "LOW"),
        ],
        "MODIFIER",
    );
}

/// Through the comparison itself the class sets aside the Consequence of a deletion's row
/// whose stop codon test is made, counted once, and with it the IMPACT that is each side's
/// most severe term's.
#[test]
fn stop_retained_three_base_stop_codon_test_is_a_built_in_class() {
    let at = StopText {
        cdna: "2558-2564",
        cds: "2249-?",
        allele: "-",
        location: "2:190970701-190970707",
    };
    let pair = [(RETAINED, "LOW"), (UNKNOWN, "MODIFIER")];
    is_a_built_in_class_of_the_following_cells(
        classes::STOP_RETAINED_THREE_BASE_STOP_CODON_TEST,
        at,
        StopText {
            cdna: "?-2564",
            cds: "-",
            ..at
        },
        pair,
        pair,
        "HIGH",
    );
}

/// Through the comparison itself the class sets aside the Consequence of a deletion's row
/// without a cDNA position at one end or a CDS start, counted once, and with it the IMPACT
/// that is each side's most severe term's.
#[test]
fn stop_codon_removed_without_stop_codon_test_is_a_built_in_class() {
    let pair = [("stop_retained_variant", "LOW"), ("stop_lost", "HIGH")];
    is_a_built_in_class_of_the_following_cells(
        classes::STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST,
        DELETION_UNTESTED,
        StopText {
            cdna: "7055-7059",
            cds: "6978-6982",
            ..DELETION_UNTESTED
        },
        pair,
        pair,
        "MODERATE",
    );
}

const INSERTION: StopText = StopText {
    cdna: "7054-7055",
    cds: "6977-6978",
    allele: "NNNN",
    location: "19:13208046-13208047",
};

/// Through the comparison itself the class sets aside the Consequence of an insertion's row
/// with both cDNA positions and a CDS start, counted once, and with it the IMPACT that is
/// each side's most severe term's.
#[test]
fn stop_codon_overlap_insertion_range_is_a_built_in_class() {
    let pair = [
        ("stop_gained,frameshift_variant", "HIGH"),
        ("inframe_insertion,stop_retained_variant", "MODERATE"),
    ];
    is_a_built_in_class_of_the_following_cells(
        classes::STOP_CODON_OVERLAP_INSERTION_RANGE,
        INSERTION,
        StopText {
            cds: "?-6978",
            ..INSERTION
        },
        pair,
        pair,
        "LOW",
    );
}

/// `stop_codon_unshifted_span` stands before every class of the stop codon tests, so a row
/// with a positive SHIFT_LENGTH that has its shape and one of theirs is its row, counted
/// under it alone.
#[test]
fn stop_codon_unshifted_span_is_asked_before_the_stop_codon_test_classes() {
    let tables = classes::tables();
    let names: Vec<&str> = TERM_PAIR_CLASSES.iter().map(|(name, _)| *name).collect();
    let index = |klass: &str| {
        names
            .iter()
            .position(|name| *name == klass)
            .expect("a row class")
    };
    let following = tables.term_pair_following_classes();
    assert!(!following.is_empty());
    for klass in following {
        assert!(index(UNSHIFTED_CLASS) < index(klass), "{klass}");
    }
    let terms = |term: &str| format!("splice_region_variant,{term},3_prime_UTR_variant");
    let shifted = "IMPACT=LOW;STRAND=-1;SHIFT_LENGTH=2";
    let text =
        |term: &str| stop_codon_text(&terms(term), "1705", "1630", shifted, "-", "3:36993663");
    let reference = text("coding_sequence_variant");
    let candidate = text("stop_retained_variant");
    let (e, a) = first_entries(&reference, &candidate);
    assert!(holds(stop_codon_unshifted_span_pair, &e, &a));
    assert!(holds(stop_codon_overlap_genomic_span_pair, &e, &a));
    let mut both: Vec<&str> = following.to_vec();
    both.push(UNSHIFTED_CLASS);
    let (report, documented) = compare_default(&reference, &candidate, &both);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(UNSHIFTED_CLASS));
    let klass = classes::STOP_CODON_OVERLAP_GENOMIC_SPAN;
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
}

/// An insertion's row with a positive SHIFT_LENGTH that has the shape of
/// `stop_codon_unshifted_span` as well is that class's where both are named: it is counted
/// under it alone and its IMPACT is compared like any other field, as on every row of that
/// class. Named alone, the class of the insertion's range takes the row and its IMPACT.
#[test]
fn stop_codon_overlap_insertion_range_leaves_a_shifted_row_to_stop_codon_unshifted_span() {
    let (klass, earlier) = (classes::STOP_CODON_OVERLAP_INSERTION_RANGE, UNSHIFTED_CLASS);
    let text = |terms: &str, impact: &str| {
        stop_codon_text(
            terms,
            INSERTION.cdna,
            INSERTION.cds,
            &format!("IMPACT={impact};STRAND=-1;SHIFT_LENGTH=2"),
            INSERTION.allele,
            INSERTION.location,
        )
    };
    let reference = text(UNKNOWN, "MODIFIER");
    let candidate = text(RETAINED, "LOW");
    let (e, a) = first_entries(&reference, &candidate);
    assert!(holds(stop_codon_unshifted_span_pair, &e, &a));
    assert!(holds(stop_codon_overlap_insertion_range_pair, &e, &a));
    let (report, documented) = compare_default(&reference, &candidate, &[klass, earlier]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), one(earlier));
    assert_eq!(masked(&report, "Consequence"), one(earlier));
    assert_eq!(differing(&report, "IMPACT"), 1);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "Extra"), one(klass));
}

/// Three made-up transcripts: the first two start at one site, the third at another.
const SITE_FIRST: &str = "ENST00000000011";
const SITE_SECOND: &str = "ENST00000000012";
const SITE_OTHER: &str = "ENST00000000013";
const NEAREST_OF_ONE_SITE: &str = "ENST00000000011,ENST00000000012";

/// What the candidate's Extra text carries on the insertion's row, and whether the
/// comparison holds the three transcripts' start sites.
#[derive(Clone, Copy)]
struct NearestRow<'a> {
    impact: &'a str,
    nearest: &'a str,
    strand: &'a str,
    sites: bool,
}

impl Default for NearestRow<'_> {
    fn default() -> Self {
        NearestRow {
            impact: "MODERATE",
            nearest: NEAREST_OF_ONE_SITE,
            strand: "-1",
            sites: true,
        }
    }
}

/// The comparison under `stop_codon_overlap_insertion_range` of one insertion's row whose
/// two sets of terms are the class's. The reference's Extra text carries the IMPACT of its
/// terms, a STRAND and the first transcript as its NEAREST; the candidate's carries the
/// `impact`, `strand` and `nearest` of `row`.
fn insertion_row_with_nearest(row: NearestRow) -> (Value, Documented) {
    let text = |terms: &str, extra: &str| {
        stop_codon_text(
            terms,
            INSERTION.cdna,
            INSERTION.cds,
            extra,
            INSERTION.allele,
            INSERTION.location,
        )
    };
    let reference = text(
        "stop_gained,frameshift_variant",
        &format!("IMPACT=HIGH;STRAND=-1;NEAREST={SITE_FIRST}"),
    );
    let candidate = text(
        "inframe_insertion,stop_retained_variant",
        &format!(
            "IMPACT={};STRAND={};NEAREST={}",
            row.impact, row.strand, row.nearest
        ),
    );
    let start_sites: IndexMap<String, BTreeSet<i64>> =
        [(SITE_FIRST, 1000), (SITE_SECOND, 1000), (SITE_OTHER, 5000)]
            .into_iter()
            .map(|(transcript, site)| (transcript.to_string(), BTreeSet::from([site])))
            .collect();
    let side = Side {
        nearest_sites: row.sites.then(|| NearestSites::new(&start_sites)),
        ..Side::default()
    };
    let mut documented = documented_with(&[classes::STOP_CODON_OVERLAP_INSERTION_RANGE]);
    let report = compare_entries(
        &parse_default(&reference).expect("the reference parses"),
        &parse_default(&candidate).expect("the candidate parses"),
        &mut documented,
        &CompareOptions::new(&side),
    )
    .expect("the entries compare");
    (report, documented)
}

fn nearest_accepted(report: &Value, column: &str) -> i64 {
    number(report, &["columns", column, "nearest_accepted"])
}

/// The row is the class's once, its Extra text among its cells and no difference left.
fn assert_extra_text_under_the_insertion_class(report: &Value, documented: &Documented) {
    let klass = classes::STOP_CODON_OVERLAP_INSERTION_RANGE;
    assert_eq!(verdict(report), "green", "{:?}", reasons(report));
    assert_eq!(applied(documented), one(klass));
    assert_eq!(masked(report, "Extra"), one(klass));
    assert_eq!(differing(report, "Extra"), 0);
    assert_eq!(nearest_accepted(report, "Extra"), 0);
}

/// The row is the class's by its Consequence; its Extra text is a difference under no class.
fn assert_extra_text_a_difference(report: &Value, documented: &Documented) {
    assert_eq!(verdict(report), "red");
    assert_eq!(
        applied(documented),
        one(classes::STOP_CODON_OVERLAP_INSERTION_RANGE)
    );
    assert_eq!(differing(report, "Extra"), 1);
    assert!(!has_masked(report, "Extra"));
    assert_eq!(nearest_accepted(report, "Extra"), 0);
}

/// The row of `NearestRow::default()` is the class's with its Extra text.
fn assert_the_row_of_one_start_site_is_the_classes() {
    let (report, documented) = insertion_row_with_nearest(NearestRow::default());
    assert_extra_text_under_the_insertion_class(&report, &documented);
}

/// On a row of the class, an Extra text that differs in the IMPACT of each side's terms and
/// in a NEAREST whose reference transcript is one of the candidate's two at one start site is
/// set aside once under the class; the NEAREST column itself agrees by its start sites.
#[test]
fn stop_codon_overlap_insertion_range_reads_a_nearest_pair_of_one_start_site_in_its_extra_text() {
    let (report, documented) = insertion_row_with_nearest(NearestRow::default());
    assert_extra_text_under_the_insertion_class(&report, &documented);
    assert_eq!(
        masked(&report, "IMPACT"),
        one(classes::STOP_CODON_OVERLAP_INSERTION_RANGE)
    );
    assert_eq!(nearest_accepted(&report, "NEAREST"), 1);
    assert_eq!(differing(&report, "NEAREST"), 0);
    assert!(!has_masked(&report, "NEAREST"));
}

/// Beside the IMPACT the class takes, a NEAREST naming a transcript of another start site
/// keeps the Extra text a difference.
#[test]
fn stop_codon_overlap_insertion_range_keeps_a_nearest_pair_of_two_start_sites_a_difference() {
    assert_the_row_of_one_start_site_is_the_classes();
    let (report, documented) = insertion_row_with_nearest(NearestRow {
        nearest: &format!("{SITE_FIRST},{SITE_OTHER}"),
        ..NearestRow::default()
    });
    assert_extra_text_a_difference(&report, &documented);
    assert_eq!(differing(&report, "NEAREST"), 1);
    let (report, documented) = insertion_row_with_nearest(NearestRow {
        nearest: SITE_SECOND,
        ..NearestRow::default()
    });
    assert_extra_text_a_difference(&report, &documented);
}

/// Beside a NEAREST pair of one start site, an IMPACT that is not its side's terms' keeps
/// the Extra text a difference.
#[test]
fn stop_codon_overlap_insertion_range_keeps_an_impact_that_is_not_its_terms_a_difference_beside_a_nearest_pair(
) {
    assert_the_row_of_one_start_site_is_the_classes();
    let (report, documented) = insertion_row_with_nearest(NearestRow {
        impact: "LOW",
        ..NearestRow::default()
    });
    assert_extra_text_a_difference(&report, &documented);
    assert_eq!(differing(&report, "IMPACT"), 1);
    assert_eq!(nearest_accepted(&report, "NEAREST"), 1);
}

/// Beside the IMPACT the class takes and a NEAREST pair of one start site, a third differing
/// pair that no class takes keeps the Extra text a difference.
#[test]
fn stop_codon_overlap_insertion_range_keeps_a_third_pair_no_class_takes_a_difference_beside_a_nearest_pair(
) {
    assert_the_row_of_one_start_site_is_the_classes();
    let (report, documented) = insertion_row_with_nearest(NearestRow {
        strand: "1",
        ..NearestRow::default()
    });
    assert_extra_text_a_difference(&report, &documented);
    assert_eq!(differing(&report, "STRAND"), 1);
}

/// A comparison that holds no start sites reads no NEAREST pair: the Extra text and the
/// NEAREST column are differences.
#[test]
fn stop_codon_overlap_insertion_range_keeps_a_nearest_pair_a_difference_without_start_sites() {
    assert_the_row_of_one_start_site_is_the_classes();
    let (report, documented) = insertion_row_with_nearest(NearestRow {
        sites: false,
        ..NearestRow::default()
    });
    assert_extra_text_a_difference(&report, &documented);
    assert_eq!(differing(&report, "NEAREST"), 1);
    assert_eq!(nearest_accepted(&report, "NEAREST"), 0);
}

/// An Extra text in which the NEAREST pair alone differs is counted under the start-site
/// rule and under no class, on a row of the class as on any other.
#[test]
fn stop_codon_overlap_insertion_range_leaves_a_text_differing_in_nearest_alone_to_the_start_site_rule(
) {
    assert_the_row_of_one_start_site_is_the_classes();
    let klass = classes::STOP_CODON_OVERLAP_INSERTION_RANGE;
    let (report, documented) = insertion_row_with_nearest(NearestRow {
        impact: "HIGH",
        ..NearestRow::default()
    });
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(nearest_accepted(&report, "Extra"), 1);
    assert_eq!(differing(&report, "Extra"), 0);
    assert!(!has_masked(&report, "Extra"));
    assert_eq!(masked(&report, "Consequence"), one(klass));
    assert!(!has_masked(&report, "IMPACT"));
}

/// The class's sentence, the same in every text that states the class.
const SPLICE_SPAN_SENTENCE: &str = "On a record whose alternate allele is longer than its reference allele, Ensembl VEP writes a splice donor, splice acceptor, donor fifth base, donor region, polypyrimidine tract or intron term from positions past the end of the record's span; where no base of the span lies in a window of that term, vep-rs does not write it, unless it is a term of an insertion written after the span.";

/// The --masks help states the class's sentence.
#[test]
fn the_masks_help_states_the_sentence_of_splice_term_past_the_record_span() {
    let help = include_str!("../data/help.txt");
    let help: Vec<&str> = help.split_whitespace().collect();
    let stated = format!("{SPLICE_SPAN_CLASS} ({SPLICE_SPAN_SENTENCE} ");
    assert!(help.join(" ").contains(&stated));
}

/// Through the comparison itself the class sets aside the Consequence of a row whose
/// record the default format states, counted once, and with it the IMPACT that is each
/// side's most severe term's. Under an identifier and with no input given, a row whose span
/// of three positions lies in the transcript's sequence under a longer allele is set aside
/// too; the same pair on a row whose span crosses an exon's end and on a row of two
/// neighbouring positions, an IMPACT the candidate leaves out, an IMPACT that is not its
/// terms' and another column that differs stay differences.
#[test]
fn splice_term_past_the_record_span_is_a_built_in_class() {
    let klass = SPLICE_SPAN_CLASS;
    for (shape, record) in splice_span_records() {
        let (reference, candidate) = record.parsed("default");
        let (report, documented) = compared(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "green", "{shape}");
        assert_eq!(applied(&documented), one(klass), "{shape}");
        assert_eq!(
            verdict(&compared(&reference, &candidate, &[]).0),
            "red",
            "{shape}"
        );
    }
    let text = |terms: &str, impact: &str, strand: &str, columns: &SpanColumns| {
        let extra = if impact.is_empty() {
            format!("STRAND={strand}")
        } else {
            format!("IMPACT={impact};STRAND={strand}")
        };
        default_text(&[splice_span_row(
            terms,
            &SpanColumns {
                allele: "ACGTAC",
                extra: &extra,
                ..*columns
            },
        )])
    };
    let plain = SpanColumns::default();
    let reference = text(DONOR_PAST_THE_SPAN_IN_FRAME, "HIGH", "1", &plain);
    let candidate = text(REGION_IN_ITS_PLACE_IN_FRAME, "MODERATE", "1", &plain);
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(masked(&report, "Consequence"), one(klass));
    assert_eq!(masked(&report, "Extra"), one(klass));
    assert_eq!(verdict_under(&reference, &candidate, &[]), "red");
    let (report, documented) = compare_default(
        &text(DONOR_PAST_THE_SPAN, "HIGH", "1", &plain),
        &text(REGION_IN_ITS_PLACE, "HIGH", "1", &plain),
        &[klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(differing(&report, "Extra"), 0);
    assert!(!has_masked(&report, "Extra"));
    let across = SpanColumns {
        cdna: "999-?",
        cds: "899-?",
        ..plain
    };
    let two_neighbours = SpanColumns {
        cdna: "999-1000",
        cds: "899-900",
        location: "7:100999-101000",
        ..plain
    };
    for columns in [across, two_neighbours] {
        let (report, documented) = compare_default(
            &text(DONOR_PAST_THE_SPAN_IN_FRAME, "HIGH", "1", &columns),
            &text(REGION_IN_ITS_PLACE_IN_FRAME, "MODERATE", "1", &columns),
            &[klass],
        );
        assert_eq!(verdict(&report), "red", "{columns:?}");
        assert_eq!(applied(&documented), none(), "{columns:?}");
        assert_eq!(differing(&report, "Consequence"), 1, "{columns:?}");
        assert_eq!(differing(&report, "Extra"), 1, "{columns:?}");
    }
    for planted in ["", "LOW"] {
        let (report, documented) = compare_default(
            &reference,
            &text(REGION_IN_ITS_PLACE_IN_FRAME, planted, "1", &plain),
            &[klass],
        );
        assert_eq!(verdict(&report), "red", "{planted}");
        assert_eq!(applied(&documented), one(klass), "{planted}");
        assert_eq!(differing(&report, "Extra"), 1, "{planted}");
        assert!(!has_masked(&report, "Extra"), "{planted}");
    }
    let (report, documented) = compare_default(
        &reference,
        &text(REGION_IN_ITS_PLACE_IN_FRAME, "MODERATE", "-1", &plain),
        &[klass],
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(differing(&report, "STRAND"), 1);
    assert!(!has_masked(&report, "Extra"));
}

/// A row on which the candidate alone holds the tract term reaches the class as any other
/// row does: the comparison has no rule over two sets of terms that counts a pair without
/// taking it, and neither of the two that take one (`EXCLUDING_RULES`) takes this pair. The
/// record's allele opens with its whole reference allele, in an intron after the tract's last
/// base on the reverse strand.
#[test]
fn splice_term_past_the_record_span_takes_a_row_on_which_the_candidate_alone_holds_the_tract_term()
{
    let klass = SPLICE_SPAN_CLASS;
    let (reference_terms, candidate_terms) = (
        "intron_variant",
        "splice_polypyrimidine_tract_variant,intron_variant",
    );
    let record = SpliceSpanRecord::new(
        ["7", "600998", "TCG", "TCGTAC"],
        "7:600999-601000",
        "CGTAC",
        "-",
        "-",
        (reference_terms, candidate_terms),
    )
    .with_impacts("MODIFIER", "LOW");
    for fmt in ["default", "vcf", "json"] {
        let (reference, candidate) = record.parsed(fmt);
        let (report, documented) = compared(&reference, &candidate, &[klass]);
        assert_eq!(verdict(&report), "green", "{fmt}");
        assert_eq!(applied(&documented), one(klass), "{fmt}");
        assert_eq!(
            verdict(&compared(&reference, &candidate, &[]).0),
            "red",
            "{fmt}"
        );
    }
    for rule in vep_formats_compare::shared::EXCLUDING_RULES {
        assert!(
            !rule.matches(
                reference_terms,
                "intron_variant,splice_polypyrimidine_tract_variant",
                "CGTAC"
            ),
            "{}",
            rule.taxonomy_class
        );
    }
}

/// A row whose record the format states is counted once in the VCF and the JSON format
/// too, where the JSON record's most_severe_consequence is left uncompared.
#[test]
fn splice_term_past_the_record_span_counts_a_row_once_in_the_vcf_and_json_formats() {
    let klass = SPLICE_SPAN_CLASS;
    for (shape, record) in splice_span_records() {
        for fmt in ["vcf", "json"] {
            let (reference, candidate) = record.parsed(fmt);
            let (report, documented) = compared(&reference, &candidate, &[klass]);
            assert_eq!(verdict(&report), "green", "{shape} {fmt}");
            assert_eq!(applied(&documented), one(klass), "{shape} {fmt}");
            assert_eq!(
                verdict(&compared(&reference, &candidate, &[]).0),
                "red",
                "{shape} {fmt}"
            );
        }
    }
}

/// A made-up record whose span begins in the 5' UTR and ends on the stop codon's last base
/// at an exon's end, its row without a CDS start: Ensembl VEP writes three intron-window
/// terms from positions past the span and a kept stop codon, vep-rs `splice_region_variant`
/// and a lost stop codon.
fn span_without_a_cds_start() -> SpliceSpanRecord {
    let gone = "splice_donor_variant,splice_donor_region_variant,intron_variant";
    SpliceSpanRecord::new(
        ["7", "500989", "CCAATGGCTTAA", "TCAATGGCTTAAGTA"],
        "7:500989-501000",
        "TCAATGGCTTAAGTA",
        "989-1000",
        "?-9",
        (
            &format!("{gone},stop_retained_variant"),
            "splice_region_variant,stop_lost",
        ),
    )
}

/// One row can differ by the terms of a stop codon class beside the class's own. A record
/// that replaces a short coding sequence, from its 5' UTR to the stop codon's last base at an
/// exon's end, by a longer allele has no CDS start: the reference writes
/// stop_retained_variant for a stop codon test it does not make and the donor, donor region
/// and intron terms from positions past the span, the candidate stop_lost and
/// splice_region_variant. Every stop codon class refuses the pair for its splice terms, and
/// this class takes it where a stop codon rule takes the rest of the difference: the row is
/// counted once, under this class, with the stop codon class named before it or not, and a
/// span across an exon's end without its cDNA end likewise. The two entries keep their terms.
#[test]
fn splice_term_past_the_record_span_takes_a_row_that_differs_by_a_stop_codon_class_as_well() {
    let klass = SPLICE_SPAN_CLASS;
    let stop_class = classes::STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST;
    let pair: Pair = splice_term_past_the_record_span_pair;
    let index = |name: &str| {
        TERM_PAIR_CLASSES
            .iter()
            .position(|(klass, _)| *klass == name)
            .expect("a row class")
    };
    assert!(index(stop_class) < index(klass));
    let gone = "splice_donor_variant,splice_donor_region_variant,intron_variant";
    let without_cds_start = span_without_a_cds_start();
    let across = splice_span_record("a span across an exon's end");
    let without_cdna_end = across.clone().with_terms(
        &across
            .terms
            .0
            .replace("coding_sequence_variant", "stop_retained_variant"),
        &across
            .terms
            .1
            .replace("coding_sequence_variant", "stop_lost"),
    );
    for record in [&without_cds_start, &without_cdna_end] {
        let (e, a) = record.entries(None);
        assert!(holds(pair, &e, &a), "{}", record.cds);
        assert_eq!(
            (&e.consequence, &a.consequence),
            (&record.terms.0, &record.terms.1),
            "{}",
            record.cds
        );
        for (name, predicate) in TERM_PAIR_CLASSES {
            if *name != klass {
                assert!(
                    !holds_without_a_cache(*predicate, &e, &a),
                    "{} {name}",
                    record.cds
                );
            }
        }
        let (reference, candidate) = record.rows(None);
        let (reference, candidate) = (default_text(&[reference]), default_text(&[candidate]));
        for classes in [&[klass][..], &[stop_class, klass]] {
            let (report, documented) = compare_default(&reference, &candidate, classes);
            assert_eq!(verdict(&report), "green", "{} {classes:?}", record.cds);
            assert_eq!(
                applied(&documented),
                one(klass),
                "{} {classes:?}",
                record.cds
            );
            assert_eq!(
                masked(&report, "Consequence"),
                one(klass),
                "{} {classes:?}",
                record.cds
            );
        }
        let (report, documented) = compare_default(&reference, &candidate, &[stop_class]);
        assert_eq!(verdict(&report), "red", "{}", record.cds);
        assert_eq!(applied(&documented), none(), "{}", record.cds);
    }
    let mut with_a_cds_start = without_cds_start.clone();
    with_a_cds_start.cds = "1-9".to_string();
    let refused = [
        ("a CDS start", with_a_cds_start),
        (
            "nothing where the span gives splice_region_variant",
            without_cds_start
                .clone()
                .with_terms(&without_cds_start.terms.0, "stop_lost"),
        ),
        (
            "the stop codon terms the other way round",
            without_cds_start.clone().with_terms(
                &format!("{gone},stop_lost"),
                "splice_region_variant,stop_retained_variant",
            ),
        ),
    ];
    for (shape, record) in refused {
        let (e, a) = record.entries(None);
        assert!(!holds(pair, &e, &a), "{shape}");
    }
    let stop_codon_alone = without_cds_start.clone().with_terms(
        "splice_region_variant,stop_retained_variant",
        "splice_region_variant,stop_lost",
    );
    let (e, a) = stop_codon_alone.entries(None);
    assert!(!holds(pair, &e, &a));
    assert!(holds(
        stop_codon_removed_without_stop_codon_test_pair,
        &e,
        &a
    ));
    let (reference, candidate) = stop_codon_alone.rows(None);
    let (report, documented) = compare_default(
        &default_text(&[reference]),
        &default_text(&[candidate]),
        &[stop_class, klass],
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(stop_class));
}

/// Whether the class `klass` accepts `name` of two rows of one key, each side its terms
/// and its value, `None` being a field the row leaves out; `columns` are those of
/// `splice_span_row`.
fn splice_span_plugin_value_taken(
    klass: &str,
    name: &str,
    reference: (&str, Option<&str>),
    candidate: (&str, Option<&str>),
    columns: &SpanColumns,
) -> bool {
    let entry = |(terms, value): (&str, Option<&str>)| {
        let extra = match value {
            None => "STRAND=1".to_string(),
            Some(value) => format!("{name}={value};STRAND=1"),
        };
        default_entry(&splice_span_row(
            terms,
            &SpanColumns {
                extra: &extra,
                ..*columns
            },
        ))
    };
    let (e, a) = (entry(reference), entry(candidate));
    class_accepts(klass, &e, name, field(&e, name), field(&a, name), Some(&a))
}

const SPLICE_SPAN_TRUNCATING_VALUES: [(&str, &str); 5] = [
    ("NMD", "NMD_escaping_variant"),
    ("LoF", "HC"),
    ("LoF_info", "PERCENTILE:0.97,GERP_DIST:12.5"),
    ("LoF_filter", "END_TRUNC"),
    ("LoF_flags", "SINGLE_EXON"),
];
const SPLICE_SPAN_DOWNSTREAM_VALUES: [(&str, &str); 2] =
    [("DownstreamProtein", "MKV"), ("ProteinLengthChange", "-12")];

/// The class accepts an NMD or LOFTEE field the reference alone writes on a row of
/// `splice_term_past_the_record_span` whose only truncating term is the donor or acceptor
/// term the candidate does not hold, and a Downstream field the candidate alone writes
/// where the reference holds frameshift_variant beside such a term and the candidate's row
/// holds no splice term. It follows the row class's two readings: without the record on a
/// row under an identifier, three positions inside an exon, and with it on a row whose first
/// column spells the alleles, two positions among them, and five whose first base alone the
/// allele replaces before it runs two bases past them. Named beside the row class, it
/// leaves the row's cells to that class; named alone, it takes them and the Consequence
/// stays a difference.
#[test]
fn plugin_value_follows_splice_term_past_the_record_span_takes_a_plugin_field_one_side_alone_writes(
) {
    let klass = SPLICE_SPAN_PLUGIN_VALUE_CLASS;
    let plain = SpanColumns::default();
    let taken = |name: &str,
                 reference: (&str, Option<&str>),
                 candidate: (&str, Option<&str>),
                 columns: &SpanColumns| {
        splice_span_plugin_value_taken(klass, name, reference, candidate, columns)
    };
    for (name, value) in SPLICE_SPAN_TRUNCATING_VALUES {
        for term in [
            "splice_donor_variant",
            "splice_acceptor_variant",
            "splice_donor_variant,intron_variant",
        ] {
            for candidate in [
                "protein_altering_variant,splice_region_variant",
                "protein_altering_variant",
            ] {
                assert!(
                    taken(
                        name,
                        (&format!("{term},protein_altering_variant"), Some(value)),
                        (candidate, None),
                        &plain
                    ),
                    "{name} {term} {candidate}"
                );
            }
        }
    }
    for (name, value) in SPLICE_SPAN_DOWNSTREAM_VALUES {
        for term in [
            "splice_donor_variant",
            "splice_acceptor_variant",
            "splice_donor_variant,splice_donor_region_variant",
            "splice_donor_variant,splice_donor_region_variant,intron_variant",
        ] {
            assert!(
                taken(
                    name,
                    (&format!("{term},frameshift_variant"), None),
                    ("frameshift_variant", Some(value)),
                    &plain
                ),
                "{name} {term}"
            );
        }
    }
    let exon = splice_span_record("two positions inside an exon");
    let spelt = exon.spelt();
    let two_positions = exon.columns(&spelt);
    let five_positions = SpanColumns {
        name: "7_100995_GACGTA/GGCGTATT",
        allele: "GCGTATT",
        location: "7:100996-101000",
        cdna: "996-1000",
        cds: "896-900",
        ..plain
    };
    for (name, value) in SPLICE_SPAN_TRUNCATING_VALUES {
        let reference = (DONOR_PAST_THE_SPAN_IN_FRAME, Some(value));
        let candidate = (REGION_IN_ITS_PLACE_IN_FRAME, None);
        assert!(taken(name, reference, candidate, &two_positions), "{name}");
        let under_an_identifier = SpanColumns {
            name: "rec2",
            ..two_positions
        };
        assert!(
            !taken(name, reference, candidate, &under_an_identifier),
            "{name}"
        );
    }
    for (name, value) in SPLICE_SPAN_DOWNSTREAM_VALUES {
        assert!(
            taken(
                name,
                (DONOR_PAST_THE_SPAN, None),
                ("frameshift_variant", Some(value)),
                &five_positions
            ),
            "{name}"
        );
    }
    let text = |terms: &str, extra: &str| {
        default_text(&[splice_span_row(
            terms,
            &SpanColumns {
                allele: "ACGTAC",
                extra,
                ..plain
            },
        )])
    };
    let reference = text(
        DONOR_PAST_THE_SPAN_IN_FRAME,
        "STRAND=1;NMD=NMD_escaping_variant;LoF=HC",
    );
    let candidate = text(REGION_IN_ITS_PLACE_IN_FRAME, "STRAND=1");
    let (report, documented) = compare_default(&reference, &candidate, &[SPLICE_SPAN_CLASS, klass]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(SPLICE_SPAN_CLASS));
    assert_eq!(masked(&report, "Extra"), one(SPLICE_SPAN_CLASS));
    let (report, documented) = compare_default(&reference, &candidate, &[klass]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), one(klass));
    assert_eq!(differing(&report, "Consequence"), 1);
    assert_eq!(masked(&report, "Extra"), one(klass));
}

/// The same values on a row the row class does not take, a field both sides write with two
/// values, a field on the side whose terms keep the plugin out, a Downstream field beside
/// the candidate's splice_region_variant, an NMD or LOFTEE field beside a truncating term
/// both sides hold, a cell the class does not name (IMPACT, HGVSp) and a row without the
/// candidate's entry stay differences; so does every field on a row whose intron term alone
/// went, which admits no plugin and keeps none out; the stop codon class takes no row of
/// this one.
#[test]
fn plugin_value_follows_splice_term_past_the_record_span_keeps_every_other_cell_a_difference() {
    let klass = SPLICE_SPAN_PLUGIN_VALUE_CLASS;
    let plain = SpanColumns::default();
    let taken = |name: &str, reference: (&str, Option<&str>), candidate: (&str, Option<&str>)| {
        splice_span_plugin_value_taken(klass, name, reference, candidate, &plain)
    };
    let in_frame = (DONOR_PAST_THE_SPAN_IN_FRAME, REGION_IN_ITS_PLACE_IN_FRAME);
    let shifted = (DONOR_PAST_THE_SPAN, "frameshift_variant");
    let inserted = inserted_between_two();
    let inserted_spelt = inserted.spelt();
    let refused_rows = [
        SpanColumns {
            cdna: "999-1000",
            cds: "899-900",
            location: "7:100999-101000",
            ..plain
        },
        SpanColumns {
            cdna: "999-?",
            cds: "899-?",
            ..plain
        },
        SpanColumns {
            allele: "ACG",
            ..plain
        },
        inserted.columns(&inserted_spelt),
    ];
    for columns in &refused_rows {
        for (name, value) in SPLICE_SPAN_TRUNCATING_VALUES {
            assert!(
                !splice_span_plugin_value_taken(
                    klass,
                    name,
                    (in_frame.0, Some(value)),
                    (in_frame.1, None),
                    columns
                ),
                "{name} {columns:?}"
            );
        }
        for (name, value) in SPLICE_SPAN_DOWNSTREAM_VALUES {
            assert!(
                !splice_span_plugin_value_taken(
                    klass,
                    name,
                    (shifted.0, None),
                    (shifted.1, Some(value)),
                    columns
                ),
                "{name} {columns:?}"
            );
        }
    }
    let intron_alone = ("frameshift_variant,intron_variant", "frameshift_variant");
    assert!(holds(
        splice_term_past_the_record_span_pair,
        &default_entry(&splice_span_row(intron_alone.0, &plain)),
        &default_entry(&splice_span_row(intron_alone.1, &plain))
    ));
    for (name, value) in SPLICE_SPAN_TRUNCATING_VALUES
        .into_iter()
        .chain(SPLICE_SPAN_DOWNSTREAM_VALUES)
    {
        assert!(
            !taken(name, (intron_alone.0, Some(value)), (intron_alone.1, None)),
            "{name}"
        );
        assert!(
            !taken(name, (intron_alone.0, None), (intron_alone.1, Some(value))),
            "{name}"
        );
    }
    for (name, value) in SPLICE_SPAN_TRUNCATING_VALUES {
        let other = format!("{value}0");
        let with_stop_lost = format!("{},stop_lost", in_frame.1);
        assert!(
            !taken(name, (in_frame.0, Some(value)), (&with_stop_lost, None)),
            "{name}"
        );
        assert!(
            !taken(name, (in_frame.0, Some(value)), (in_frame.1, Some(&other))),
            "{name}"
        );
        assert!(
            !taken(name, (in_frame.0, None), (in_frame.1, Some(value))),
            "{name}"
        );
        assert!(
            !taken(name, (shifted.0, Some(value)), (shifted.1, None)),
            "{name}"
        );
        assert!(
            !taken(name, (shifted.0, Some(value)), (shifted.1, Some(&other))),
            "{name}"
        );
        assert!(
            !splice_span_plugin_value_taken(
                PLUGIN_VALUE_CLASS,
                name,
                (in_frame.0, Some(value)),
                (in_frame.1, None),
                &plain
            ),
            "{name}"
        );
    }
    for (name, value) in SPLICE_SPAN_DOWNSTREAM_VALUES {
        let other = format!("{value}0");
        let with_stop_gained = format!("{},stop_gained", shifted.1);
        assert!(
            !taken(name, (shifted.0, None), (&with_stop_gained, Some(value))),
            "{name}"
        );
        assert!(
            !taken(name, (shifted.0, Some(&other)), (shifted.1, Some(value))),
            "{name}"
        );
        assert!(
            !taken(name, (shifted.0, Some(value)), (shifted.1, None)),
            "{name}"
        );
        assert!(
            !taken(name, (shifted.0, None), (REGION_IN_ITS_PLACE, Some(value))),
            "{name}"
        );
        assert!(
            !taken(name, (in_frame.0, None), (in_frame.1, Some(value))),
            "{name}"
        );
        assert!(
            !splice_span_plugin_value_taken(
                PLUGIN_VALUE_CLASS,
                name,
                (shifted.0, None),
                (shifted.1, Some(value)),
                &plain
            ),
            "{name}"
        );
    }
    assert!(!taken(
        "IMPACT",
        (in_frame.0, Some("HIGH")),
        (in_frame.1, Some("MODERATE"))
    ));
    assert!(splice_span_plugin_value_taken(
        SPLICE_SPAN_CLASS,
        "IMPACT",
        (in_frame.0, Some("HIGH")),
        (in_frame.1, Some("MODERATE")),
        &plain
    ));
    let protein = "ENSP90000000001.1:p.";
    assert!(!taken(
        "HGVSp",
        (shifted.0, Some(&format!("{protein}Tyr22CysfsTer4"))),
        (shifted.1, Some(&format!("{protein}Tyr22CysfsTer9")))
    ));
    let entry = default_entry(&splice_span_row(
        in_frame.0,
        &SpanColumns {
            extra: "NMD=NMD_escaping_variant;STRAND=1",
            ..plain
        },
    ));
    assert!(!class_accepts(
        klass,
        &entry,
        "NMD",
        "NMD_escaping_variant",
        ABSENT,
        None
    ));
    let stop_codon = |terms: &str, extra: &str| {
        default_entry(&stop_codon_row("-", terms, "1705-?", "1630-?").replace("IMPACT=LOW", extra))
    };
    let kept = stop_codon("stop_retained_variant", "IMPACT=LOW");
    let lost = stop_codon("frameshift_variant", "IMPACT=HIGH;NMD=NMD_escaping_variant");
    assert!(class_accepts(
        PLUGIN_VALUE_CLASS,
        &kept,
        "NMD",
        ABSENT,
        "NMD_escaping_variant",
        Some(&lost)
    ));
    assert!(!class_accepts(
        klass,
        &kept,
        "NMD",
        ABSENT,
        "NMD_escaping_variant",
        Some(&lost)
    ));
}

/// Beside the candidate's cache the rule that takes the rest of such a difference reads
/// the coding ends of the row's transcript. Where a stop codon stands, the lost stop codon
/// against the kept one is the rest, as without a cache, and so is the coding term against
/// the kept stop codon, the pair of the untested stop codon. Where none stands the lost
/// stop codon is a wrong row of vep-rs and the whole difference stays, and the row vep-rs
/// owes there, the coding term in the stop term's place, is taken beside the splice
/// difference: the rest is the pair of `stop_codon_overlap_without_stop_codon`.
#[test]
fn splice_difference_beside_a_pair_of_a_coding_end() {
    use coding_end_common::{closed, open, planted_cache, Transcript};
    let klass = SPLICE_SPAN_CLASS;
    let lost = span_without_a_cds_start();
    let owed = lost.clone().with_terms(
        &lost.terms.0,
        "splice_region_variant,coding_sequence_variant",
    );
    let feature = "ENST90000000001";
    // The span ends on the genomic base 501000: a transcript that ends there with no
    // stop codon, and one whose stop codon lies under the span.
    let none_stands = Transcript {
        id: feature,
        first: 500980,
        ..open()
    };
    let one_stands = Transcript {
        id: feature,
        first: 500977,
        ..closed()
    };
    let compared = |record: &SpliceSpanRecord, transcript: &Transcript, name: &str| {
        let directory = scratch(&format!("splice-beside-a-coding-end-{name}"));
        let cache = planted_cache(&directory, std::slice::from_ref(transcript));
        let side = Side {
            partial_codon_models: Some(PartialCodonModels::new(&cache).expect("the cache")),
            ..Side::default()
        };
        let (reference, candidate) = record.rows(None);
        let expected = parse_default(&default_text(&[reference])).expect("the reference parses");
        let actual = parse_default(&default_text(&[candidate])).expect("the candidate parses");
        let mut documented = documented_with(&[klass]);
        let report = compare_entries(
            &expected,
            &actual,
            &mut documented,
            &CompareOptions::new(&side),
        )
        .expect("the entries compare");
        (verdict(&report).to_string(), applied(&documented))
    };
    let taken = ("green".to_string(), one(klass));
    let left = ("red".to_string(), none());
    assert_eq!(compared(&lost, &one_stands, "lost-stands"), taken);
    assert_eq!(compared(&lost, &none_stands, "lost-none"), left);
    assert_eq!(compared(&owed, &none_stands, "owed-none"), taken);
    assert_eq!(compared(&owed, &one_stands, "owed-stands"), taken);
}
