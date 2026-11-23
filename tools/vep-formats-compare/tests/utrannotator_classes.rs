// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The UTRAnnotator classes of a row's term: the precedence among the plugin's analyses,
//! a structural allele's text, a frameshift at a stop codon's boundary and a lost stop
//! written before the stop codon, read through the comparison of planted rows.

mod compare_common;
mod support;
mod utrannotator_common;

use std::collections::BTreeSet;

use compare_common::*;
use utrannotator_common::*;
use vep_formats_compare::consts::STRUCTURAL_INSERTION_WORD;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::lists::order_only_fields;
use vep_formats_compare::rules::utrannotator::annotation_groups;
use vep_formats_compare::rules::utrannotator::STRUCTURAL_CLASS_WORDS;
use vep_formats_compare::rules::utrannotator::STRUCTURAL_DELETION_WORDS;
use vep_formats_compare::rules::utrannotator::UTRANNOTATOR_DELETION_TERMS;
use vep_formats_compare::rules::utrannotator::UTRANNOTATOR_REFERENCE_UTR_KEYS;
use vep_formats_compare::rules::utrannotator::UTRANNOTATOR_TERMS;

/// A substitution on which two of the plugin's analyses fire.
const UTR_SUBSTITUTION: Row = ("rs1", "21:35748195", "T", "ENST1", "293");

fn utr_lost_start() -> Values {
    values(
        UTR_LOSS,
        "CapDistanceToStart=290:KozakStrength=Moderate:type=uORF",
    )
}

fn utr_lost_stop() -> Values {
    values(
        UTR_STOP_LOSS,
        "AltStop=False:KozakContext=GGGATGC:KozakStrength=Moderate",
    )
}

/// A substitution on which two of the plugin's analyses fire: the reference reports the
/// lost start, the candidate the lost stop, and Ensembl VEP's own runs under two hash seeds
/// write one term each on the row, the second with its annotation's pairs in another order.
/// The row is counted once, its 5UTR_annotation with its 5UTR_consequence, and a row the
/// engines agree on is untouched.
#[test]
fn the_consequence_precedence_accepts_two_terms_the_seeded_runs_write_on_the_row() {
    let (reference, candidate, other_seed) = precedence_texts();
    let sweep = [reference.as_str(), other_seed.as_str()];
    let outcome = compared_swept(&reference, &candidate, Some(&sweep), &[UTR_PRECEDENCE]);
    assert!(outcome.is("green", &[(UTR_PRECEDENCE, 1)]));
    assert!(outcome.masks("5UTR_consequence", &[(UTR_PRECEDENCE, 1)]));
    assert!(outcome.masks("5UTR_annotation", &[(UTR_PRECEDENCE, 1)]));
    let unnamed = compared_swept(&reference, &candidate, Some(&sweep), &[]);
    assert_eq!(unnamed.verdict(), "red");
    // A deletion whose reference reports a start past the deleted bases beside another term of
    // the candidate's is the window class's when both are named, and the precedence class's
    // when it alone is and the seeded runs write both terms.
    let row = at_cdna(UTR_DELETION, "336-340");
    let deletion = one_row(utr_row(row, &gained(337)));
    let deletion_candidate = one_row(utr_row(row, &utr_shifted()));
    let both_terms = [deletion.as_str(), deletion_candidate.as_str()];
    let outcome = compared_swept(
        &deletion,
        &deletion_candidate,
        Some(&both_terms),
        &[UTR_WINDOW, UTR_PRECEDENCE],
    );
    assert!(outcome.is("green", &[(UTR_WINDOW, 1)]));
    let outcome = compared_swept(
        &deletion,
        &deletion_candidate,
        Some(&both_terms),
        &[UTR_PRECEDENCE],
    );
    assert!(outcome.is("green", &[(UTR_PRECEDENCE, 1)]));
}

/// The reference's and the candidate's outputs of the substitution beside a row the engines
/// agree on, and the output of a second seeded run.
fn precedence_texts() -> (String, String, String) {
    let agreeing = utr_row(("rs2", "21:400", "C", "ENST2", "40"), &gained(38));
    let reference = default_text(&[
        utr_row(UTR_SUBSTITUTION, &utr_lost_start()),
        agreeing.clone(),
    ]);
    let candidate = default_text(&[
        utr_row(UTR_SUBSTITUTION, &utr_lost_stop()),
        agreeing.clone(),
    ]);
    let other_seed = default_text(&[
        utr_row(
            UTR_SUBSTITUTION,
            &values(
                UTR_STOP_LOSS,
                "KozakStrength=Moderate:AltStop=False:KozakContext=GGGATGC",
            ),
        ),
        agreeing,
    ]);
    (reference, candidate, other_seed)
}

/// The command line reads one output per seed.
#[test]
fn the_command_line_reads_one_seeded_output_a_sweep_option() {
    let (reference, candidate, other_seed) = precedence_texts();
    let directory = planted_files(
        "utrannotator-precedence",
        &[
            ("reference", &reference),
            ("candidate", &candidate),
            ("seed_a", &reference),
            ("seed_b", &other_seed),
        ],
    );
    assert_eq!(exit_code(&directory, &["--masks", UTR_PRECEDENCE]), 1);
    assert_eq!(
        exit_code(
            &directory,
            &[
                "--masks",
                UTR_PRECEDENCE,
                "--utrannotator-sweep",
                "seed_a.txt",
                "--utrannotator-sweep",
                "seed_b.txt"
            ]
        ),
        0
    );
}

/// The same value pair stays a difference without the seeded runs, on a row where they
/// write the reference's term alone, so one analysis fires, on a row they do not hold, and
/// on the row of another record at the Location. A - on either side and a value outside the
/// five terms stay differences beside runs that write both terms, an annotation those runs
/// do not write beside the candidate's term stays one under an accepted consequence, or with
/// its consequence where the candidate shares it with the reference, and the VCF and JSON
/// entries, which name no Location, match no row of the runs.
#[test]
fn the_consequence_precedence_refuses_a_term_the_seeded_runs_never_write() {
    let reference = one_row(utr_row(UTR_SUBSTITUTION, &utr_lost_start()));
    let candidate = one_row(utr_row(UTR_SUBSTITUTION, &utr_lost_stop()));
    let both_terms = [reference.as_str(), candidate.as_str()];
    let elsewhere = one_row(utr_row(
        ("rs1", "21:35748196", "T", "ENST1", "294"),
        &utr_lost_stop(),
    ));
    let another_record = one_row(utr_row(
        ("rs9", "21:35748195", "T", "ENST1", "293"),
        &utr_lost_stop(),
    ));
    let sweeps: [Option<Vec<&str>>; 5] = [
        None,
        Some(vec![&reference, &reference]),
        Some(vec![&elsewhere]),
        Some(vec![&reference, &another_record]),
        Some(vec![&candidate]),
    ];
    for sweep in &sweeps {
        let outcome = compared_swept(&reference, &candidate, sweep.as_deref(), &[UTR_PRECEDENCE]);
        assert!(outcome.is("red", &[]), "{sweep:?}");
        assert_eq!(outcome.differing("5UTR_consequence"), 1);
        assert_eq!(outcome.differing("5UTR_annotation"), 1);
    }
    let wrong = [
        utr_row(UTR_SUBSTITUTION, &nothing()),
        utr_row(
            UTR_SUBSTITUTION,
            &values("5_prime_UTR_variant", &utr_lost_stop().1),
        ),
        utr_row(UTR_SUBSTITUTION, &utr_shifted()),
    ];
    for wrong in wrong {
        let outcome = compared_swept(
            &reference,
            &one_row(wrong.clone()),
            Some(&both_terms),
            &[UTR_PRECEDENCE],
        );
        assert!(outcome.is("red", &[]), "{wrong}");
    }
    let outcome = compared_swept(
        &one_row(utr_row(UTR_SUBSTITUTION, &nothing())),
        &candidate,
        Some(&both_terms),
        &[UTR_PRECEDENCE],
    );
    assert!(outcome.is("red", &[]));
    let outcome = compared_swept(
        &reference,
        &candidate.replace("AltStop=False", "AltStop=True"),
        Some(&both_terms),
        &[UTR_PRECEDENCE],
    );
    assert!(outcome.is("red", &[(UTR_PRECEDENCE, 1)]));
    assert!(outcome.masks("5UTR_consequence", &[(UTR_PRECEDENCE, 1)]));
    assert_eq!(outcome.differing("5UTR_annotation"), 1);
    // The reference's own annotation beside the candidate's term, as written or with its pairs
    // in another order, is one the runs write beside the reference's term alone.
    let reordered = "type=uORF:CapDistanceToStart=290:KozakStrength=Moderate";
    for shared in [utr_lost_start().1.as_str(), reordered] {
        let outcome = compared_swept(
            &reference,
            &one_row(utr_row(UTR_SUBSTITUTION, &values(UTR_STOP_LOSS, shared))),
            Some(&both_terms),
            &[UTR_PRECEDENCE],
        );
        assert!(outcome.is("red", &[]), "{shared}");
        assert_eq!(outcome.differing("5UTR_consequence"), 1);
    }
    // order_only reads the reordered pairs as one set, and the consequence stays a difference.
    let outcome = compared_in_order(
        &default_entries(&reference),
        &default_entries(&one_row(utr_row(
            UTR_SUBSTITUTION,
            &values(UTR_STOP_LOSS, reordered),
        ))),
        &["order_only", UTR_PRECEDENCE],
        None,
        Some(seed_sweep(&both_terms)),
        Some(&order_only_fields("default")),
    );
    assert!(outcome.is("red", &[("order_only", 1)]));
    let inputs = lines(&["21\t35748195\trs1\tG\tT\t.\t.\t."]);
    let header =
        csq_header("Allele|Consequence|Feature|Feature_type|5UTR_consequence|5UTR_annotation");
    let vcf_reference = format!(
        "{header}{}CSQ=T|5_prime_UTR_variant|ENST1|Transcript|{UTR_LOSS}|type=uORF\n",
        record_columns(&inputs[0])
    );
    let vcf_candidate = vcf_reference
        .replace(UTR_LOSS, UTR_STOP_LOSS)
        .replace("type=uORF", "AltStop=False");
    let outcome = compared_entries(
        &vcf_entries(&vcf_reference, &inputs),
        &vcf_entries(&vcf_candidate, &inputs),
        &[UTR_PRECEDENCE],
        None,
        Some(seed_sweep(&both_terms)),
    );
    assert!(outcome.is("red", &[]));

    let object = |utr: &str, annotation: &str| {
        json_line(
            &inputs[0],
            "rs1",
            "G/T",
            &format!(
                r#""variant_allele": "T", "transcript_id": "ENST1", "consequence_terms": ["5_prime_UTR_variant"], "5utr_consequence": "{utr}", "5utr_annotation": {{"1": {annotation}}}"#
            ),
        )
    };
    let outcome = compared_entries(
        &json_entries(&object(UTR_LOSS, r#"{"type": "uORF"}"#), &inputs),
        &json_entries(&object(UTR_STOP_LOSS, r#"{"AltStop": "False"}"#), &inputs),
        &[UTR_PRECEDENCE],
        None,
        Some(seed_sweep(&both_terms)),
    );
    assert!(outcome.is("red", &[]));
    let e = &default_entries(&reference)[0];
    let a = &default_entries(&candidate)[0];
    let accepted = |a: Option<&Entry>, swept: bool| {
        accepts(
            UTR_PRECEDENCE,
            e,
            "5UTR_consequence",
            UTR_LOSS,
            UTR_STOP_LOSS,
            a,
            None,
            swept.then(|| seed_sweep(&both_terms)),
        )
    };
    assert!(accepted(Some(a), true));
    assert!(!accepted(Some(a), false));
    assert!(!accepted(None, true));
    assert_eq!(
        annotation_groups("b=2:a=1,c=3"),
        annotation_groups("c=3,a=1:b=2")
    );
}

/// The command line refuses a sweep file that is missing or is no default-format output.
#[test]
fn the_command_line_refuses_a_sweep_file_that_is_no_default_format_output() {
    let reference = one_row(utr_row(UTR_SUBSTITUTION, &utr_lost_start()));
    let candidate = one_row(utr_row(UTR_SUBSTITUTION, &utr_lost_stop()));
    let vcf = format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t35748195\trs1\tG\tT\t.\t.\tCSQ=T|5_prime_UTR_variant|ENST1|Transcript|{UTR_STOP_LOSS}|AltStop=False\n"
    );
    let directory = planted_files(
        "utrannotator-sweep-files",
        &[
            ("reference", &reference),
            ("candidate", &candidate),
            ("vcf", &vcf),
            ("empty", ""),
        ],
    );
    let with_sweep = |second: &str| {
        exit_code(
            &directory,
            &[
                "--masks",
                UTR_PRECEDENCE,
                "--utrannotator-sweep",
                "reference.txt",
                "--utrannotator-sweep",
                second,
            ],
        )
    };
    for sweep_file in ["vcf.txt", "empty.txt", "absent.txt"] {
        assert_eq!(with_sweep(sweep_file), 2, "{sweep_file}");
    }
    assert_eq!(with_sweep("candidate.txt"), 0);
}

fn utr_spliced_text() -> Values {
    values(
        UTR_FRAMESHIFT,
        "alt_type=uORF:alt_type_length=39:ref_type=uORF:ref_type_length=51",
    )
}

const UTR_EXISTING: &str = ";Existing_InFrame_oORFs=0;Existing_OutOfFrame_oORFs=0;Existing_uORFs=6";

/// The class words a structural row carries as its Allele.
const CLASS_WORDS: [&str; 17] = [
    "Alu_deletion",
    "Alu_insertion",
    "HERV_deletion",
    "HERV_insertion",
    "LINE1_deletion",
    "LINE1_insertion",
    "SVA_deletion",
    "SVA_insertion",
    "copy_number_variation",
    "deletion",
    "duplication",
    "insertion",
    "inversion",
    "mobile_element_deletion",
    "mobile_element_insertion",
    "tandem_duplication",
    "tandem_repeat",
];

/// A structural record's 5' UTR row under `allele` with the two UTRAnnotator fields and
/// `more` after them in its Extra.
fn structural_utr_row_of(allele: &str, own: &Values, more: &str, name: &str) -> String {
    utr_row(
        (name, "21:10524589-10524628", allele, "ENST1", "?-208"),
        own,
    ) + more
}

fn structural_utr_row(allele: &str, own: &Values) -> String {
    structural_utr_row_of(allele, own, "", "sv")
}

/// The comparison of two outputs under the structural class, with the seeded runs `sweep`
/// where they are given.
fn structural<S: AsRef<str>>(reference: &str, candidate: &str, sweep: Option<&[S]>) -> Outcome {
    compared_swept(reference, candidate, sweep, &[UTR_STRUCTURAL])
}

/// The comparison of two outputs under the structural class with no seeded run.
fn structural_alone(reference: &str, candidate: &str) -> Outcome {
    compared(reference, candidate, &[UTR_STRUCTURAL])
}

/// Beside a structural row the class takes: a sequence allele's row, a breakend's own
/// string, the breakpoint class's word and a word that names no structural class stay
/// differences; so does a structural row the candidate fills, alone or beside the reference,
/// with either field or both, a reference term outside the plugin's five, an annotation
/// beside no term, the plugin's Existing_ counts and another plugin's field, with the seeded
/// runs given or not.
#[test]
fn the_structural_allele_text_refuses_a_row_outside_the_class() {
    let filled = utr_spliced_text();
    let other = values(UTR_LOSS, "CapDistanceToStart=63:type=uORF");
    let taken = one_row(structural_utr_row("deletion", &filled));
    let seeded = [
        taken.clone(),
        one_row(structural_utr_row("deletion", &other)),
    ];
    let sweeps: [Option<&[String]>; 2] = [None, Some(&seeded)];
    let outcome = structural_alone(&taken, &one_row(structural_utr_row("deletion", &nothing())));
    assert!(outcome.is("green", &[(UTR_STRUCTURAL, 1)]));
    for allele in [
        "-",
        "T",
        "ACGTN",
        "N.",
        ".N",
        "N[21:25881361[",
        "]21:25881361]N",
        "chromosome_breakpoint",
        "complex_structural_alteration",
        "Deletion",
        "deletion2",
    ] {
        let reference = one_row(structural_utr_row(allele, &filled));
        let candidate = one_row(structural_utr_row(allele, &nothing()));
        for sweep in sweeps {
            let outcome = structural(&reference, &candidate, sweep);
            assert!(outcome.is("red", &[]), "{allele}");
            assert_eq!(outcome.differing("5UTR_consequence"), 1, "{allele}");
            assert_eq!(outcome.differing("5UTR_annotation"), 1, "{allele}");
        }
    }
    let pairs: [(&str, Values, Values); 7] = [
        ("the candidate alone", nothing(), filled.clone()),
        ("another term on the candidate's row", filled.clone(), other),
        (
            "the same term with another annotation",
            filled.clone(),
            values(&filled.0, "alt_type=uORF:alt_type_length=12"),
        ),
        (
            "the candidate's term beside no annotation",
            filled.clone(),
            values(&filled.0, "-"),
        ),
        (
            "the candidate's annotation beside no term",
            filled.clone(),
            values("-", &filled.1),
        ),
        (
            "a reference term outside the plugin's five",
            values("5_prime_UTR_variant", &filled.1),
            nothing(),
        ),
        (
            "a reference annotation beside no term",
            values("-", &filled.1),
            nothing(),
        ),
    ];
    for (what, reference_fields, candidate_fields) in &pairs {
        let reference = one_row(structural_utr_row("deletion", reference_fields));
        let candidate = one_row(structural_utr_row("deletion", candidate_fields));
        for sweep in sweeps {
            let outcome = structural(&reference, &candidate, sweep);
            assert!(outcome.is("red", &[]), "{what}");
        }
    }
    // The plugin's other fields and another plugin's field differ beside the two the class takes.
    for (more, column) in [(UTR_EXISTING, "Existing_uORFs"), (";LOEUF=0.31", "LOEUF")] {
        for (reference_more, candidate_more) in [(more, ""), ("", more)] {
            let reference = one_row(structural_utr_row_of(
                "deletion",
                &filled,
                reference_more,
                "sv",
            ));
            let candidate = one_row(structural_utr_row_of(
                "deletion",
                &nothing(),
                candidate_more,
                "sv",
            ));
            let outcome = structural_alone(&reference, &candidate);
            assert!(outcome.is("red", &[(UTR_STRUCTURAL, 1)]), "{column}");
            assert_eq!(outcome.differing(column), 1);
            assert!(outcome.masks("5UTR_consequence", &[(UTR_STRUCTURAL, 1)]));
            let same = one_row(structural_utr_row_of(
                "deletion",
                &nothing(),
                reference_more,
                "sv",
            ));
            let outcome = structural_alone(&same, &candidate);
            assert!(outcome.is("red", &[]), "{column}");
        }
    }
    let e = &default_entries(&taken)[0];
    let a = &default_entries(&one_row(structural_utr_row_of(
        "deletion",
        &nothing(),
        UTR_EXISTING,
        "sv",
    )))[0];
    let accepted = |name: &str, expected_value: &str| {
        accepts(
            UTR_STRUCTURAL,
            e,
            name,
            expected_value,
            "-",
            Some(a),
            None,
            None,
        )
    };
    assert!(accepted("5UTR_consequence", &filled.0));
    assert!(!accepted("Existing_uORFs", "6"));
    assert!(!accepted("LOEUF", "0.31"));
}

/// A record with two symbolic alternates writes two rows on a transcript under one Allele
/// word. The reference fills the two fields on both, its annotation's pairs in another order
/// on the second, and the candidate on neither: both rows are set aside. A candidate that
/// holds one row less leaves a row without its pair, which stays a difference.
#[test]
fn the_structural_allele_text_reads_each_row_of_a_record_with_two_alternates() {
    let word = "copy_number_variation";
    let reordered = "ref_type_length=51:ref_type=uORF:alt_type_length=39:alt_type=uORF";
    let row = |own: &Values| structural_utr_row_of(word, own, "", "cnv");
    let reference = default_text(&[
        row(&utr_spliced_text()),
        row(&values(UTR_FRAMESHIFT, reordered)),
    ]);
    let candidate = default_text(&[row(&nothing()), row(&nothing())]);
    let outcome = structural_alone(&reference, &candidate);
    assert!(outcome.is("green", &[(UTR_STRUCTURAL, 2)]));
    assert!(outcome.masks("5UTR_consequence", &[(UTR_STRUCTURAL, 2)]));
    assert!(outcome.masks("5UTR_annotation", &[(UTR_STRUCTURAL, 2)]));
    let one = one_row(row(&nothing()));
    let outcome = structural_alone(&reference, &one);
    assert!(outcome.is("red", &[(UTR_STRUCTURAL, 1)]));
}

/// The VCF format carries the two fields as CSQ members under the row's Allele word and
/// the JSON format as `5utr_consequence` and the `5utr_annotation` object beside
/// `variant_allele`: the class takes both where the reference alone fills them, and a
/// sequence allele's row of either format stays a difference.
#[test]
fn the_structural_allele_text_reads_the_vcf_and_json_spellings() {
    let inputs = lines(&[
        "21\t10524588\tsv\tN\t<DEL>\t.\t.\tEND=10524628",
        "21\t10524700\tsnv\tA\tG\t.\t.\t.",
    ]);
    let header =
        csq_header("Allele|Consequence|Feature|Feature_type|5UTR_consequence|5UTR_annotation");
    let spliced = utr_spliced_text();
    let filled = format!("{}|{}", spliced.0, spliced.1);
    let vcf_reference = format!(
        "{header}{};CSQ=deletion|5_prime_UTR_variant|ENST1|Transcript|{filled}\n{}CSQ=G|5_prime_UTR_variant|ENST1|Transcript|{filled}\n",
        inputs[0],
        record_columns(&inputs[1])
    );
    let vcf_candidate = vcf_reference.replace(&filled, "|");
    let outcome = compared_entries(
        &vcf_entries(&vcf_reference, &inputs),
        &vcf_entries(&vcf_candidate, &inputs),
        &[UTR_STRUCTURAL],
        None,
        None,
    );
    assert!(outcome.is("red", &[(UTR_STRUCTURAL, 1)]));
    for column in ["5UTR_consequence", "5UTR_annotation"] {
        assert!(outcome.masks(column, &[(UTR_STRUCTURAL, 1)]));
        assert_eq!(outcome.differing(column), 1);
    }

    let object = |at: usize, allele_string: &str, allele: &str, fills: bool| {
        let mut entry = format!(
            r#""variant_allele": "{allele}", "transcript_id": "ENST1", "consequence_terms": ["5_prime_UTR_variant"]"#
        );
        if fills {
            entry.push_str(&format!(
                r#", "5utr_consequence": "{UTR_FRAMESHIFT}", "5utr_annotation": {{"1": {{"alt_type": "uORF", "alt_type_length": 39}}}}"#
            ));
        }
        json_line(&inputs[at], &format!("r{at}"), allele_string, &entry)
    };
    let both =
        |fills: bool| object(0, "deletion", "deletion", fills) + &object(1, "A/G", "G", fills);
    let outcome = compared_entries(
        &json_entries(&both(true), &inputs),
        &json_entries(&both(false), &inputs),
        &[UTR_STRUCTURAL],
        None,
        None,
    );
    assert!(outcome.is("red", &[(UTR_STRUCTURAL, 1)]));
    for column in ["5utr_consequence", "5utr_annotation"] {
        assert!(outcome.masks(column, &[(UTR_STRUCTURAL, 1)]));
        assert_eq!(outcome.differing(column), 1);
    }
}

const UTR_SPAN_LOCATION: &str = "21:10524589-10524628";
const UTR_SPAN_CDNA: &str = "169-208";
const UTR_DELETION_WORDS: [&str; 6] = [
    "Alu_deletion",
    "HERV_deletion",
    "LINE1_deletion",
    "SVA_deletion",
    "deletion",
    "mobile_element_deletion",
];
const UTR_DELETION_TERMS: [&str; 4] = [UTR_GAIN, UTR_LOSS, UTR_STOP_LOSS, UTR_FRAMESHIFT];
const SPAN_TERMS: &str = "feature_truncation,5_prime_UTR_variant";

fn utr_span_effect() -> Values {
    values(
        UTR_STOP_LOSS,
        "AltStop=True:AltStopDistanceToCDS=159:FrameWithCDS=inFrame",
    )
}

/// Where a deletion's row stands and what its Consequence names beside its two fields.
#[derive(Clone, Copy)]
struct SpanShape<'a> {
    location: &'a str,
    cdna: &'a str,
    terms: &'a str,
    more: &'a str,
}

impl Default for SpanShape<'_> {
    fn default() -> Self {
        SpanShape {
            location: UTR_SPAN_LOCATION,
            cdna: UTR_SPAN_CDNA,
            terms: SPAN_TERMS,
            more: "",
        }
    }
}

/// A structural record's row whose Location spans 40 bases between two cDNA positions 40
/// apart, under the Consequence terms of a deletion inside one exon of the 5' UTR, with the
/// two UTRAnnotator fields and `more` after them in its Extra.
fn deletion_span_row_of(allele: &str, own: &Values, shape: SpanShape) -> String {
    let row = utr_row(("sv", shape.location, allele, "ENST1", shape.cdna), own) + shape.more;
    assert_eq!(row.matches("\t5_prime_UTR_variant\t").count(), 1);
    row.replace("\t5_prime_UTR_variant\t", &format!("\t{}\t", shape.terms))
}

fn deletion_span_row(allele: &str, own: &Values) -> String {
    deletion_span_row_of(allele, own, SpanShape::default())
}

/// The candidate annotates a symbolic deletion whose stated span lies inside one exon of
/// the 5' UTR as the sequence deletion of that span, where Ensembl VEP reports the effect of
/// a deletion as long as the ALT text or none. Under each class word of a deletion, with
/// each term a deletion can fire, the two fields the candidate fills are set aside, once a
/// row: beside an empty reference, beside another term and beside the same term with another
/// annotation, where the annotation alone differs; with the seeded runs given or not, and
/// under the NMD term as well. A frameshift is the class's where the span is no multiple of
/// three. The row keeps the reference-alone shape: a candidate that fills neither field
/// there is set aside as on every structural row. Equal values are no difference, and the
/// plugin's Existing_ counts and another plugin's field stay differences beside the two
/// fields the class takes.
#[test]
fn the_structural_allele_text_accepts_the_fields_the_candidate_fills_on_a_deletion_inside_one_utr_exon(
) {
    let text_length = utr_spliced_text();
    for word in UTR_DELETION_WORDS {
        for term in UTR_DELETION_TERMS {
            let effect = values(term, "KozakStrength=Weak:type=uORF");
            let candidate = one_row(deletion_span_row(word, &effect));
            let other = if term != text_length.0 {
                text_length.clone()
            } else {
                utr_lost()
            };
            for reference_fields in [nothing(), other] {
                let reference = one_row(deletion_span_row(word, &reference_fields));
                let seeded = [reference.clone()];
                let sweeps: [Option<&[String]>; 2] = [None, Some(&seeded)];
                for sweep in sweeps {
                    let outcome = structural(&reference, &candidate, sweep);
                    assert!(outcome.is("green", &[(UTR_STRUCTURAL, 1)]), "{word} {term}");
                    assert!(
                        outcome.masks("5UTR_consequence", &[(UTR_STRUCTURAL, 1)]),
                        "{word} {term}"
                    );
                    assert!(
                        outcome.masks("5UTR_annotation", &[(UTR_STRUCTURAL, 1)]),
                        "{word} {term}"
                    );
                    let unnamed = compared_swept(&reference, &candidate, sweep, &[]);
                    assert_eq!(unnamed.verdict(), "red", "{word} {term}");
                }
            }
            let same_term = one_row(deletion_span_row(
                word,
                &values(term, "KozakStrength=Strong:type=uORF"),
            ));
            let outcome = structural_alone(&same_term, &candidate);
            assert!(outcome.is("green", &[(UTR_STRUCTURAL, 1)]), "{word} {term}");
            assert!(outcome.masks("5UTR_consequence", &[]));
            assert!(outcome.masks("5UTR_annotation", &[(UTR_STRUCTURAL, 1)]));
        }
    }
    let nmd = SpanShape {
        terms: "feature_truncation,5_prime_UTR_variant,NMD_transcript_variant",
        ..SpanShape::default()
    };
    let outcome = structural_alone(
        &one_row(deletion_span_row_of("deletion", &nothing(), nmd)),
        &one_row(deletion_span_row_of("deletion", &utr_span_effect(), nmd)),
    );
    assert!(outcome.is("green", &[(UTR_STRUCTURAL, 1)]));
    // A one-base span between equal cDNA positions, and a frameshift on a span of 40.
    let one_base = SpanShape {
        location: "21:10524589",
        cdna: "169",
        ..SpanShape::default()
    };
    for shape in [one_base, SpanShape::default()] {
        let outcome = structural_alone(
            &one_row(deletion_span_row_of("deletion", &nothing(), shape)),
            &one_row(deletion_span_row_of(
                "deletion",
                &values(UTR_FRAMESHIFT, "ref_type=uORF"),
                shape,
            )),
        );
        assert!(
            outcome.is("green", &[(UTR_STRUCTURAL, 1)]),
            "{}",
            shape.location
        );
    }
    // The reference-alone shape on the same row.
    let outcome = structural_alone(
        &one_row(deletion_span_row("deletion", &text_length)),
        &one_row(deletion_span_row("deletion", &nothing())),
    );
    assert!(outcome.is("green", &[(UTR_STRUCTURAL, 1)]));
    // Equal values are no difference, and the row's other fields stay differences.
    let taken = one_row(deletion_span_row("deletion", &utr_span_effect()));
    assert!(structural_alone(&taken, &taken).is("green", &[]));
    for (more, column) in [(UTR_EXISTING, "Existing_uORFs"), (";LOEUF=0.31", "LOEUF")] {
        for (reference_more, candidate_more) in [(more, ""), ("", more)] {
            let with_more = |more| SpanShape {
                more,
                ..SpanShape::default()
            };
            let outcome = structural_alone(
                &one_row(deletion_span_row_of(
                    "deletion",
                    &nothing(),
                    with_more(reference_more),
                )),
                &one_row(deletion_span_row_of(
                    "deletion",
                    &utr_span_effect(),
                    with_more(candidate_more),
                )),
            );
            assert!(outcome.is("red", &[(UTR_STRUCTURAL, 1)]), "{column}");
            assert_eq!(outcome.differing(column), 1);
            assert!(outcome.masks("5UTR_consequence", &[(UTR_STRUCTURAL, 1)]));
        }
    }
    let set =
        |words: &[&'static str]| -> BTreeSet<&'static str> { words.iter().copied().collect() };
    assert_eq!(set(STRUCTURAL_DELETION_WORDS), set(&UTR_DELETION_WORDS));
    assert!(set(STRUCTURAL_DELETION_WORDS).is_subset(&set(STRUCTURAL_CLASS_WORDS)));
    assert!(STRUCTURAL_DELETION_WORDS.len() < STRUCTURAL_CLASS_WORDS.len());
    assert_eq!(set(UTRANNOTATOR_DELETION_TERMS), set(&UTR_DELETION_TERMS));
    let other_terms: Vec<&str> = set(UTRANNOTATOR_TERMS)
        .difference(&set(UTRANNOTATOR_DELETION_TERMS))
        .copied()
        .collect();
    assert_eq!(other_terms, ["5_prime_UTR_uORF_stop_codon_gain_variant"]);
    let e = &default_entries(&one_row(deletion_span_row("deletion", &nothing())))[0];
    let a = &default_entries(&one_row(deletion_span_row("deletion", &utr_span_effect())))[0];
    let effect = utr_span_effect();
    for (name, value) in [
        ("5UTR_consequence", &effect.0),
        ("5UTR_annotation", &effect.1),
    ] {
        assert!(
            accepts(UTR_STRUCTURAL, e, name, "-", value, Some(a), None, None),
            "{name}"
        );
        assert!(
            !accepts(UTR_STRUCTURAL, e, name, "-", value, None, None, None),
            "{name}"
        );
    }
    assert!(!accepts(
        UTR_STRUCTURAL,
        e,
        "Existing_uORFs",
        "-",
        "6",
        Some(a),
        None,
        None
    ));
}

/// The class follows the plugin's pair-order, gained-start window and precedence classes,
/// so on a deletion's row the candidate fills a difference one of them explains is counted
/// under that one: two terms the seeded runs both write on the row under the precedence
/// class, with the annotations written beside them, and an annotation whose pairs differ in
/// order alone under the pair-order class.
#[test]
fn the_structural_allele_text_leaves_a_row_one_of_the_plugins_other_classes_explains_to_that_class()
{
    let reference = one_row(deletion_span_row("deletion", &utr_lost()));
    let candidate = one_row(deletion_span_row("deletion", &utr_shifted()));
    let classes = [UTR_PRECEDENCE, UTR_STRUCTURAL];
    let outcome = compared_swept(
        &reference,
        &candidate,
        Some(&[&reference, &candidate][..]),
        &classes,
    );
    assert!(outcome.is("green", &[(UTR_PRECEDENCE, 1)]));
    assert!(outcome.masks("5UTR_consequence", &[(UTR_PRECEDENCE, 1)]));
    let outcome = compared_swept(&reference, &candidate, Some(&[&reference][..]), &classes);
    assert!(outcome.is("green", &[(UTR_STRUCTURAL, 1)]));
    let lost = utr_lost();
    let reversed: Vec<&str> = lost.1.split(':').rev().collect();
    let reordered = one_row(deletion_span_row(
        "deletion",
        &values(&lost.0, &reversed.join(":")),
    ));
    let outcome = compared(&reordered, &reference, &[PAIR_ORDER, UTR_STRUCTURAL]);
    assert!(outcome.is("green", &[(PAIR_ORDER, 1)]));
    let order = plugin_classes();
    let place = |class: &str| order.iter().position(|held| *held == class);
    for earlier in [PAIR_ORDER, UTR_WINDOW, UTR_PRECEDENCE] {
        assert!(place(earlier) < place(UTR_STRUCTURAL), "{earlier}");
    }
}

/// The class words that are no deletion's, and no plain insertion's.
fn words_but(excluded: &[&str]) -> Vec<&'static str> {
    CLASS_WORDS
        .iter()
        .copied()
        .filter(|word| !excluded.contains(word))
        .collect()
}

/// A row the candidate fills under a deletion's class word is the class's only as the row
/// of a symbolic deletion inside one exon of the 5' UTR. It stays a difference under every
/// other class word but `inversion` and `copy_number_variation`, whose rows have shapes of
/// their own, and under a
/// sequence allele; with a term a deletion cannot fire, a term beside no annotation or an
/// annotation beside no term; where the cDNA positions are unknown, open at an end or not as
/// far apart as the Location's span; where the Consequence names a term of a span that
/// leaves the exon or the 5' UTR; beside a reference term outside the plugin's five or a
/// reference annotation beside no term; and in the VCF and JSON formats, whose entries name
/// no Location.
#[test]
fn the_structural_allele_text_refuses_a_row_the_candidate_fills_outside_the_deletion_shape() {
    let filled = utr_span_effect();
    let taken = one_row(deletion_span_row("deletion", &filled));
    let refused = |reference: &str, candidate: &str, what: &str| {
        let seeded = [reference, candidate];
        let sweeps: [Option<&[&str]>; 2] = [None, Some(&seeded)];
        for sweep in sweeps {
            let outcome = structural(reference, candidate, sweep);
            assert!(outcome.is("red", &[]), "{what}");
        }
    };
    let mut shaped = UTR_DELETION_WORDS.to_vec();
    shaped.extend(["inversion", "copy_number_variation"]);
    let mut alleles = words_but(&shaped);
    alleles.extend(["-", "T", "N[21:25881361[", "Deletion", "deletion2"]);
    for allele in alleles {
        for reference_fields in [nothing(), utr_spliced_text()] {
            refused(
                &one_row(deletion_span_row(allele, &reference_fields)),
                &one_row(deletion_span_row(allele, &filled)),
                allele,
            );
        }
    }
    let shape = SpanShape::default();
    let multiple_of_three = SpanShape {
        location: "21:10524589-10524627",
        cdna: "169-207",
        ..shape
    };
    let at_cdna = |cdna| SpanShape { cdna, ..shape };
    let under = |terms| SpanShape { terms, ..shape };
    let rows: [(&str, SpanShape, Values); 14] = [
        (
            "the gained stop, a substitution's term",
            shape,
            values("5_prime_UTR_uORF_stop_codon_gain_variant", &filled.1),
        ),
        (
            "a term outside the plugin's five",
            shape,
            values("5_prime_UTR_variant", &filled.1),
        ),
        (
            "a frameshift on a span that is a multiple of three",
            multiple_of_three,
            values(UTR_FRAMESHIFT, "ref_type=uORF"),
        ),
        (
            "a term beside no annotation",
            shape,
            values(&filled.0, "-"),
        ),
        (
            "an annotation beside no term",
            shape,
            values("-", &filled.1),
        ),
        (
            "cDNA positions open at the start",
            at_cdna("?-208"),
            filled.clone(),
        ),
        (
            "cDNA positions open at the end",
            at_cdna("169-?"),
            filled.clone(),
        ),
        ("no cDNA positions", at_cdna("-"), filled.clone()),
        (
            "cDNA positions closer than the span, two exons",
            at_cdna("169-190"),
            filled.clone(),
        ),
        (
            "cDNA positions further apart than the span",
            at_cdna("169-250"),
            filled.clone(),
        ),
        (
            "a span into the intron",
            under("feature_truncation,5_prime_UTR_variant,intron_variant"),
            filled.clone(),
        ),
        (
            "a span into the coding sequence",
            under("feature_truncation,coding_sequence_variant,5_prime_UTR_variant"),
            filled.clone(),
        ),
        (
            "a span over the start codon",
            under("frameshift_variant,start_lost,feature_truncation,start_retained_variant,5_prime_UTR_variant"),
            filled.clone(),
        ),
        (
            "a row outside the 5' UTR",
            under("feature_truncation,3_prime_UTR_variant"),
            filled.clone(),
        ),
    ];
    for (what, shape, candidate_fields) in &rows {
        refused(
            &one_row(deletion_span_row_of("deletion", &nothing(), *shape)),
            &one_row(deletion_span_row_of("deletion", candidate_fields, *shape)),
            what,
        );
    }
    for (what, reference_fields) in [
        (
            "a reference term outside the plugin's five",
            values("5_prime_UTR_variant", "type=uORF"),
        ),
        (
            "a reference annotation beside no term",
            values("-", "type=uORF"),
        ),
    ] {
        refused(
            &one_row(deletion_span_row("deletion", &reference_fields)),
            &taken,
            what,
        );
    }
    // The VCF and JSON formats name no Location: a filled entry stays a difference there.
    let inputs = lines(&["21\t10524588\tsv\tN\t<DEL>\t.\t.\tEND=10524628"]);
    let header = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|5UTR_consequence|5UTR_annotation",
    );
    let entry = format!(
        "deletion|feature_truncation&5_prime_UTR_variant|ENST1|Transcript|{UTR_SPAN_CDNA}|"
    );
    let vcf = |fields: &str| format!("{header}{};CSQ={entry}{fields}\n", inputs[0]);
    let outcome = compared_entries(
        &vcf_entries(&vcf("|"), &inputs),
        &vcf_entries(&vcf(&format!("{}|{}", filled.0, filled.1)), &inputs),
        &[UTR_STRUCTURAL],
        None,
        None,
    );
    assert!(outcome.is("red", &[]));

    let object = |fills: bool| {
        let mut consequence = r#""variant_allele": "deletion", "transcript_id": "ENST1", "consequence_terms": ["feature_truncation", "5_prime_UTR_variant"], "cdna_start": 169, "cdna_end": 208"#.to_string();
        if fills {
            consequence.push_str(&format!(
                r#", "5utr_consequence": "{}", "5utr_annotation": {{"1": {{"AltStop": "True"}}}}"#,
                filled.0
            ));
        }
        json_line(&inputs[0], "sv", "<DEL>", &consequence)
    };
    let outcome = compared_entries(
        &json_entries(&object(false), &inputs),
        &json_entries(&object(true), &inputs),
        &[UTR_STRUCTURAL],
        None,
        None,
    );
    assert!(outcome.is("red", &[]));
}

/// A span of three bases and one of four, each between cDNA positions as far apart, under
/// the Consequence of a row that is no deletion's.
const UTR_SPAN_OF_THREE: SpanShape<'static> = SpanShape {
    location: "21:14210905-14210907",
    cdna: "49-51",
    terms: "5_prime_UTR_variant",
    more: "",
};
const UTR_SPAN_OF_FOUR: SpanShape<'static> = SpanShape {
    location: "21:10524589-10524592",
    cdna: "152-155",
    terms: "5_prime_UTR_variant",
    more: "",
};

/// What Ensembl VEP writes on the rows of a structural record over each span: the frameshift
/// of a deletion as long as the first ALT text.
const UTR_TEXT_OF_THREE: &str = "alt_type_length=48:KozakContext=GTGATGA:KozakStrength=Moderate:ref_type_length=1215:ref_type=InFrame_oORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=87";
const UTR_TEXT_OF_FOUR: &str = "alt_type_length=39:KozakContext=AGAATGT:KozakStrength=Moderate:ref_type_length=51:ref_type=uORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=207";

/// The comparison of two outputs under the gained-start window and the structural class
/// with the one input line of the record `sv`, whose ALT is `alternates`.
fn structural_with_record(
    reference: &str,
    candidate: &str,
    alternates: &str,
    shape: SpanShape,
) -> Outcome {
    let (start, end) = shape
        .location
        .split_once(':')
        .and_then(|(_, span)| span.split_once('-'))
        .expect("a span");
    let before: i64 = start.parse::<i64>().expect("a position") - 1;
    let inputs = [format!(
        "21\t{before}\tsv\tN\t{alternates}\t.\t.\tEND={end}"
    )];
    compared_with(
        reference,
        candidate,
        &[UTR_WINDOW, UTR_STRUCTURAL],
        Some(&inputs[..]),
    )
}

/// The candidate annotates a symbolic inversion whose stated span lies inside one exon of
/// the 5' UTR as that span replaced by its reverse complement, where Ensembl VEP reports the
/// effect of a deletion as long as the ALT text or none. With each term a replacement by as
/// many bases can fire, the two fields the candidate fills are set aside, once a row, beside
/// an empty reference and beside another term, on a span of three bases or of four and
/// under the NMD term as well. Where the input states the record's alternates, the row is
/// the class's under `<INV>` alone.
#[test]
fn the_structural_allele_text_accepts_a_replacements_term_on_an_inversion_inside_one_utr_exon() {
    let text_length = values(UTR_FRAMESHIFT, UTR_TEXT_OF_THREE);
    let gained_stop = values(
        UTR_STOP_GAIN,
        "Evidence=False:KozakContext=GTGATGA:KozakStrength=Moderate:newSTOPDistanceToCDS=60:ref_StartDistanceToCDS=87:ref_type=InFrame_oORF",
    );
    let row =
        |own: &Values, shape: SpanShape| one_row(deletion_span_row_of("inversion", own, shape));
    let (reference, candidate) = (
        row(&text_length, UTR_SPAN_OF_THREE),
        row(&gained_stop, UTR_SPAN_OF_THREE),
    );
    let outcome = structural_alone(&reference, &candidate);
    assert!(outcome.is("green", &[(UTR_STRUCTURAL, 1)]));
    for name in ["5UTR_consequence", "5UTR_annotation"] {
        assert!(outcome.masks(name, &[(UTR_STRUCTURAL, 1)]), "{name}");
    }
    assert_eq!(compared(&reference, &candidate, &[]).verdict(), "red");
    let nmd = SpanShape {
        terms: "5_prime_UTR_variant,NMD_transcript_variant",
        ..UTR_SPAN_OF_FOUR
    };
    for term in [UTR_GAIN, UTR_LOSS, UTR_STOP_GAIN, UTR_STOP_LOSS] {
        let effect = values(term, "KozakStrength=Weak:type=uORF");
        for reference_fields in [nothing(), text_length.clone()] {
            for shape in [UTR_SPAN_OF_THREE, UTR_SPAN_OF_FOUR, nmd] {
                let outcome =
                    structural_alone(&row(&reference_fields, shape), &row(&effect, shape));
                assert!(
                    outcome.is("green", &[(UTR_STRUCTURAL, 1)]),
                    "{term} {}",
                    shape.location
                );
            }
        }
    }
    for (alternates, taken) in [
        ("<INV>", true),
        ("<INV>,<DEL>", false),
        ("<DEL>,<INV>", false),
        ("<DUP>", false),
    ] {
        let outcome = structural_with_record(&reference, &candidate, alternates, UTR_SPAN_OF_THREE);
        assert!(
            outcome.taken_is(taken, &[(UTR_STRUCTURAL, 1)]),
            "{alternates}"
        );
    }
}

/// An inversion's row the candidate fills stays a difference with the frameshift, which no
/// replacement by as many bases fires, with a term outside the plugin's five, a term beside
/// no annotation or an annotation beside no term; where the cDNA positions are open at an
/// end or not as far apart as the Location's span; where the Consequence names a deletion's
/// term or a term of a span that leaves the exon or the 5' UTR; and beside a reference term
/// outside the plugin's five or a reference annotation beside no term.
#[test]
fn the_structural_allele_text_refuses_an_inversions_row_outside_the_replaced_span_shape() {
    let filled = values(UTR_LOSS, "KozakStrength=Weak:type=uORF");
    let shape = UTR_SPAN_OF_FOUR;
    let at_cdna = |cdna| SpanShape { cdna, ..shape };
    let under = |terms| SpanShape { terms, ..shape };
    let rows: [(&str, SpanShape, Values); 11] = [
        (
            "the frameshift on a span of four",
            shape,
            values(UTR_FRAMESHIFT, "ref_type=uORF"),
        ),
        (
            "the frameshift on a span of three",
            UTR_SPAN_OF_THREE,
            values(UTR_FRAMESHIFT, "ref_type=uORF"),
        ),
        (
            "a term outside the plugin's five",
            shape,
            values("5_prime_UTR_variant", &filled.1),
        ),
        ("a term beside no annotation", shape, values(&filled.0, "-")),
        (
            "an annotation beside no term",
            shape,
            values("-", &filled.1),
        ),
        (
            "cDNA positions open at the start",
            at_cdna("?-155"),
            filled.clone(),
        ),
        (
            "cDNA positions closer than the span, two exons",
            at_cdna("152-153"),
            filled.clone(),
        ),
        ("no cDNA positions", at_cdna("-"), filled.clone()),
        (
            "a deletion's Consequence",
            under("feature_truncation,5_prime_UTR_variant"),
            filled.clone(),
        ),
        (
            "a span into the intron",
            under("5_prime_UTR_variant,intron_variant"),
            filled.clone(),
        ),
        (
            "a span into the coding sequence",
            under("coding_sequence_variant,5_prime_UTR_variant"),
            filled.clone(),
        ),
    ];
    let row =
        |own: &Values, shape: SpanShape| one_row(deletion_span_row_of("inversion", own, shape));
    for (what, shape, candidate_fields) in &rows {
        let outcome = structural_alone(&row(&nothing(), *shape), &row(candidate_fields, *shape));
        assert!(outcome.is("red", &[]), "{what}");
    }
    for (what, reference_fields) in [
        (
            "a reference term outside the plugin's five",
            values("5_prime_UTR_variant", "type=uORF"),
        ),
        (
            "a reference annotation beside no term",
            values("-", "type=uORF"),
        ),
    ] {
        let outcome = structural_alone(&row(&reference_fields, shape), &row(&filled, shape));
        assert!(outcome.is("red", &[]), "{what}");
    }
}

/// A record with a copy-number loss beside a gain writes two rows on a transcript under
/// `copy_number_variation`. Ensembl VEP writes the effect of a deletion as long as the
/// first ALT text on both; the candidate annotates the loss alternate's row as the deletion
/// of the record's span and leaves the other's two fields empty. Both rows are set aside:
/// the loss alternate's with the deletion's term, the frameshift only where the span is no
/// multiple of three, and the other as the reference alone fills it. Where the input states
/// the record's alternates, the filled row is the class's only as the row of an alternate
/// whose own text is a deletion's, in a record of two or more.
#[test]
fn the_structural_allele_text_accepts_a_deletions_term_on_the_loss_row_of_a_copy_number_record() {
    let word = "copy_number_variation";
    let text_length = values(UTR_FRAMESHIFT, UTR_TEXT_OF_FOUR);
    let span_deleted = values(
        UTR_FRAMESHIFT,
        "Evidence=False:KozakContext=AGAATGT:KozakStrength=Moderate:alt_type=uORF:alt_type_length=102:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51",
    );
    let row = |own: &Values| deletion_span_row_of(word, own, UTR_SPAN_OF_FOUR);
    let reference = default_text(&[row(&text_length), row(&text_length)]);
    let candidate = default_text(&[row(&span_deleted), row(&nothing())]);
    let outcome = structural_alone(&reference, &candidate);
    assert!(outcome.is("green", &[(UTR_STRUCTURAL, 2)]));
    assert!(outcome.masks("5UTR_consequence", &[(UTR_STRUCTURAL, 1)]));
    assert!(outcome.masks("5UTR_annotation", &[(UTR_STRUCTURAL, 2)]));
    for (alternates, taken) in [
        ("<CN0>,<CN2>", true),
        ("<DEL>,<DUP>", true),
        ("<CN=0>,<CN=2>", true),
        ("<CN2>,<CN0>", false),
        ("<DUP>,<INV>", false),
        ("<CN0>,<CN2>,<CN3>", false),
    ] {
        let outcome = structural_with_record(&reference, &candidate, alternates, UTR_SPAN_OF_FOUR);
        if taken {
            assert!(outcome.is("green", &[(UTR_STRUCTURAL, 2)]), "{alternates}");
        } else {
            assert!(outcome.is("red", &[(UTR_STRUCTURAL, 1)]), "{alternates}");
            assert_eq!(outcome.differing("5UTR_annotation"), 1, "{alternates}");
        }
    }
    // The gain alternate's row first: the candidate fills the second row.
    let later = default_text(&[row(&nothing()), row(&span_deleted)]);
    for (alternates, taken) in [("<CN2>,<CN0>", true), ("<CN0>,<CN2>", false)] {
        let outcome = structural_with_record(&reference, &later, alternates, UTR_SPAN_OF_FOUR);
        assert_eq!(outcome.green(), taken, "{alternates}");
    }
    // A record of one alternate writes one row, which is none of a loss beside a gain.
    let (one_reference, one_candidate) = (one_row(row(&text_length)), one_row(row(&span_deleted)));
    assert!(structural_alone(&one_reference, &one_candidate).is("green", &[(UTR_STRUCTURAL, 1)]));
    let outcome = structural_with_record(&one_reference, &one_candidate, "<CN0>", UTR_SPAN_OF_FOUR);
    assert!(outcome.is("red", &[]));
    // On a span of three bases a deletion fires any of its terms but the frameshift.
    let of_three = |own: &Values| one_row(deletion_span_row_of(word, own, UTR_SPAN_OF_THREE));
    for (term, taken) in [
        (UTR_LOSS, true),
        (UTR_GAIN, true),
        (UTR_STOP_LOSS, true),
        (UTR_FRAMESHIFT, false),
        (UTR_STOP_GAIN, false),
    ] {
        for reference_fields in [nothing(), values(UTR_FRAMESHIFT, UTR_TEXT_OF_THREE)] {
            if reference_fields.0 == term {
                continue;
            }
            let outcome = structural_alone(
                &of_three(&reference_fields),
                &of_three(&values(term, "KozakStrength=Weak:type=uORF")),
            );
            assert!(outcome.taken_is(taken, &[(UTR_STRUCTURAL, 1)]), "{term}");
        }
    }
}

/// A copy-number row the candidate fills stays a difference where the cDNA positions are
/// not as far apart as the Location's span, so the span crosses exons, or are open at an
/// end or unknown; where the Consequence names a deletion's term or a term of a span that
/// leaves the exon or the 5' UTR; with a term beside no annotation; and beside a reference
/// term outside the plugin's five.
#[test]
fn the_structural_allele_text_refuses_a_copy_number_row_outside_the_deleted_span_shape() {
    let filled = values(UTR_FRAMESHIFT, "KozakStrength=Weak:ref_type=uORF");
    let shape = UTR_SPAN_OF_FOUR;
    let at_cdna = |cdna| SpanShape { cdna, ..shape };
    let under = |terms| SpanShape { terms, ..shape };
    let rows: [(&str, SpanShape, Values); 8] = [
        (
            "a span that crosses exons",
            at_cdna("152-153"),
            filled.clone(),
        ),
        (
            "cDNA positions further apart than the span",
            at_cdna("152-160"),
            filled.clone(),
        ),
        (
            "cDNA positions open at the end",
            at_cdna("152-?"),
            filled.clone(),
        ),
        ("no cDNA positions", at_cdna("-"), filled.clone()),
        (
            "a deletion's Consequence",
            under("feature_truncation,5_prime_UTR_variant"),
            filled.clone(),
        ),
        (
            "a span into the intron",
            under("5_prime_UTR_variant,intron_variant"),
            filled.clone(),
        ),
        (
            "a row outside the 5' UTR",
            under("3_prime_UTR_variant"),
            filled.clone(),
        ),
        ("a term beside no annotation", shape, values(&filled.0, "-")),
    ];
    let row = |own: &Values, shape: SpanShape| {
        one_row(deletion_span_row_of("copy_number_variation", own, shape))
    };
    for (what, shape, candidate_fields) in &rows {
        let outcome = structural_alone(&row(&nothing(), *shape), &row(candidate_fields, *shape));
        assert!(outcome.is("red", &[]), "{what}");
    }
    let outcome = structural_alone(
        &row(&values("5_prime_UTR_variant", "type=uORF"), shape),
        &row(&filled, shape),
    );
    assert!(outcome.is("red", &[]));
}

const UTR_LENGTH_LOCATION: &str = "21:10524589-10524701";
const UTR_REFERENCE_KEYS: &str = "Evidence=False:KozakContext=GGGATGC:KozakStrength=Moderate:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51";
const UTR_REFERENCE_KEYS_AND_TYPE: &str = "Evidence=False:KozakContext=GGGATGC:KozakStrength=Moderate:alt_type=uORF:ref_StartDistanceToCDS=207:ref_type=uORF:ref_type_length=51";
const LENGTH_TERMS: &str = "feature_elongation,5_prime_UTR_variant";

fn utr_length_effect() -> Values {
    values(UTR_FRAMESHIFT, UTR_REFERENCE_KEYS)
}

fn utr_text_effect() -> Values {
    values(
        UTR_FRAMESHIFT,
        "alt_type_length=39:KozakContext=GGGATGC:KozakStrength=Moderate:ref_type_length=51:ref_type=uORF:alt_type=uORF:Evidence=False:ref_StartDistanceToCDS=207",
    )
}

/// A structural record's row whose Location spans the 113 bases of a stated length, under
/// the Consequence terms of an insertion in the 5' UTR, with the two UTRAnnotator fields.
fn insertion_length_row_of(allele: &str, own: &Values, location: &str, terms: &str) -> String {
    let row = utr_row(("sv", location, allele, "ENST1", "169-170"), own);
    assert_eq!(row.matches("\t5_prime_UTR_variant\t").count(), 1);
    row.replace("\t5_prime_UTR_variant\t", &format!("\t{terms}\t"))
}

fn insertion_length_row(allele: &str, own: &Values) -> String {
    insertion_length_row_of(allele, own, UTR_LENGTH_LOCATION, LENGTH_TERMS)
}

/// The candidate gives a plain symbolic insertion whose stated length is no multiple of
/// three the frameshift that length decides, with the keys read from the reference UTR,
/// where Ensembl VEP reports the effect of a deletion as long as the ALT text or none. Beside
/// the reference's frameshift the consequence is equal and is not handed to the class, and
/// the annotation is set aside, in one effect of the six keys and in two effects, one of
/// them with `alt_type=uORF`; beside another of the plugin's terms and beside an empty
/// reference both fields are set aside; the row's other Consequence terms do not matter.
#[test]
fn the_structural_allele_text_accepts_the_frameshift_the_candidate_writes_on_an_insertion_of_a_stated_length(
) {
    let (six, six_and_type) = (UTR_REFERENCE_KEYS, UTR_REFERENCE_KEYS_AND_TYPE);
    let reference = one_row(insertion_length_row("insertion", &utr_text_effect()));
    for annotation in [six.to_string(), format!("{six},{six_and_type}")] {
        let candidate = one_row(insertion_length_row(
            "insertion",
            &values(UTR_FRAMESHIFT, &annotation),
        ));
        let outcome = structural_alone(&reference, &candidate);
        assert!(outcome.is("green", &[(UTR_STRUCTURAL, 1)]), "{annotation}");
        assert!(outcome.masks("5UTR_consequence", &[]));
        assert!(outcome.masks("5UTR_annotation", &[(UTR_STRUCTURAL, 1)]));
        assert_eq!(
            compared(&reference, &candidate, &[]).verdict(),
            "red",
            "{annotation}"
        );
    }
    for (reference_fields, annotation) in [(utr_lost(), six), (nothing(), six_and_type)] {
        let reference = one_row(insertion_length_row("insertion", &reference_fields));
        let candidate = one_row(insertion_length_row(
            "insertion",
            &values(UTR_FRAMESHIFT, annotation),
        ));
        let outcome = structural_alone(&reference, &candidate);
        assert!(
            outcome.is("green", &[(UTR_STRUCTURAL, 1)]),
            "{reference_fields:?}"
        );
        assert!(outcome.masks("5UTR_consequence", &[(UTR_STRUCTURAL, 1)]));
        assert!(outcome.masks("5UTR_annotation", &[(UTR_STRUCTURAL, 1)]));
    }
    for more in ["NMD_transcript_variant", "intron_variant"] {
        let terms = format!("feature_elongation,5_prime_UTR_variant,{more}");
        let row = |own: &Values| {
            one_row(insertion_length_row_of(
                "insertion",
                own,
                UTR_LENGTH_LOCATION,
                &terms,
            ))
        };
        let outcome = structural_alone(&row(&utr_text_effect()), &row(&utr_length_effect()));
        assert!(outcome.is("green", &[(UTR_STRUCTURAL, 1)]), "{more}");
    }
    let e = &default_entries(&one_row(insertion_length_row(
        "insertion",
        &utr_text_effect(),
    )))[0];
    let a = &default_entries(&one_row(insertion_length_row(
        "insertion",
        &utr_length_effect(),
    )))[0];
    let accepted = |a: Option<&Entry>| {
        accepts(
            UTR_STRUCTURAL,
            e,
            "5UTR_annotation",
            &utr_text_effect().1,
            &utr_length_effect().1,
            a,
            None,
            None,
        )
    };
    assert!(accepted(Some(a)));
    assert!(!accepted(None));
    assert_eq!(STRUCTURAL_INSERTION_WORD, "insertion");
    assert!(STRUCTURAL_CLASS_WORDS.contains(&STRUCTURAL_INSERTION_WORD));
    let keys: BTreeSet<&str> = six
        .split(':')
        .map(|pair| pair.split('=').next().unwrap_or(pair))
        .collect();
    let held: BTreeSet<&str> = UTRANNOTATOR_REFERENCE_UTR_KEYS.iter().copied().collect();
    assert_eq!(held, keys);
}

/// A row the candidate fills under `insertion` is the class's only with the frameshift
/// and, in each effect, the six keys read from the reference UTR, with or without
/// `alt_type=uORF`. It stays a difference where an effect names `alt_type_length`,
/// another `alt_type`, five of the six keys or an eighth key; with another of the plugin's
/// terms or with the frameshift beside no annotation; on a span that is a multiple of three;
/// under a mobile element's word, a duplication's, an inversion's and every other class
/// word; where the Consequence names no `5_prime_UTR_variant`; beside a reference term
/// outside the plugin's five or a reference annotation beside no term; and in the VCF and
/// JSON formats, whose entries name no Location.
#[test]
fn the_structural_allele_text_refuses_an_insertions_row_outside_the_stated_length_shape() {
    let (six, filled, text_length) = (UTR_REFERENCE_KEYS, utr_length_effect(), utr_text_effect());
    let refused = |reference_fields: &Values,
                   candidate_fields: &Values,
                   what: &str,
                   allele: &str,
                   location: &str,
                   terms: &str| {
        let reference = one_row(insertion_length_row_of(
            allele,
            reference_fields,
            location,
            terms,
        ));
        let candidate = one_row(insertion_length_row_of(
            allele,
            candidate_fields,
            location,
            terms,
        ));
        let seeded = [reference.as_str(), candidate.as_str()];
        let sweeps: [Option<&[&str]>; 2] = [None, Some(&seeded)];
        for sweep in sweeps {
            let outcome = structural(&reference, &candidate, sweep);
            assert!(outcome.is("red", &[]), "{what}");
        }
    };
    let plain = |reference_fields: &Values, candidate_fields: &Values, what: &str| {
        refused(
            reference_fields,
            candidate_fields,
            what,
            "insertion",
            UTR_LENGTH_LOCATION,
            LENGTH_TERMS,
        )
    };
    let effects: [(&str, String); 5] = [
        (
            "an effect with alt_type_length",
            format!("{UTR_REFERENCE_KEYS_AND_TYPE}:alt_type_length=66"),
        ),
        (
            "an effect with another alt_type",
            six.replace(
                ":ref_StartDistanceToCDS",
                ":alt_type=OutOfFrame_oORF:ref_StartDistanceToCDS",
            ),
        ),
        (
            "an effect with five of the six keys",
            six.replace("Evidence=False:", ""),
        ),
        (
            "an effect with an eighth key",
            format!("{UTR_REFERENCE_KEYS_AND_TYPE}:type=uORF"),
        ),
        (
            "a second effect with alt_type_length",
            format!("{six},{UTR_REFERENCE_KEYS_AND_TYPE}:alt_type_length=66"),
        ),
    ];
    for (what, annotation) in &effects {
        assert!(annotation != six && annotation != UTR_REFERENCE_KEYS_AND_TYPE);
        plain(&text_length, &values(UTR_FRAMESHIFT, annotation), what);
    }
    plain(
        &text_length,
        &values(UTR_STOP_LOSS, six),
        "another of the plugin's terms beside the six keys",
    );
    plain(
        &text_length,
        &values(UTR_FRAMESHIFT, "-"),
        "the frameshift beside no annotation",
    );
    refused(
        &text_length,
        &filled,
        "a span that is a multiple of three",
        "insertion",
        "21:10524589-10524702",
        LENGTH_TERMS,
    );
    for allele in words_but(&["insertion"]) {
        refused(
            &text_length,
            &filled,
            allele,
            allele,
            UTR_LENGTH_LOCATION,
            LENGTH_TERMS,
        );
    }
    refused(
        &text_length,
        &filled,
        "a Consequence without 5_prime_UTR_variant",
        "insertion",
        UTR_LENGTH_LOCATION,
        "feature_elongation,intron_variant",
    );
    for (what, reference_fields) in [
        (
            "a reference term outside the plugin's five",
            values("5_prime_UTR_variant", "type=uORF"),
        ),
        (
            "a reference annotation beside no term",
            values("-", "type=uORF"),
        ),
    ] {
        plain(&reference_fields, &filled, what);
    }
    // The VCF and JSON formats name no Location: a filled entry stays a difference there.
    let inputs = lines(&["21\t10524588\tsv\tN\t<INS>\t.\t.\tSVLEN=113"]);
    let header =
        csq_header("Allele|Consequence|Feature|Feature_type|5UTR_consequence|5UTR_annotation");
    let entry = format!(
        "insertion|feature_elongation&5_prime_UTR_variant|ENST1|Transcript|{UTR_FRAMESHIFT}|"
    );
    let vcf = |annotation: &str| format!("{header}{};CSQ={entry}{annotation}\n", inputs[0]);
    let outcome = compared_entries(
        &vcf_entries(&vcf(&text_length.1), &inputs),
        &vcf_entries(&vcf(six), &inputs),
        &[UTR_STRUCTURAL],
        None,
        None,
    );
    assert!(outcome.is("red", &[]));
    assert_eq!(outcome.differing("5UTR_annotation"), 1);

    let object = |annotation: &str| {
        let effect: Vec<String> = annotation
            .split(':')
            .map(|pair| {
                let (key, value) = pair.split_once('=').expect("a pair");
                format!("{}: {}", quoted(key), quoted(value))
            })
            .collect();
        json_line(
            &inputs[0],
            "sv",
            "<INS>",
            &format!(
                r#""variant_allele": "insertion", "transcript_id": "ENST1", "consequence_terms": ["feature_elongation", "5_prime_UTR_variant"], "5utr_consequence": "{UTR_FRAMESHIFT}", "5utr_annotation": {{"1": {{{}}}}}"#,
                effect.join(", ")
            ),
        )
    };
    let outcome = compared_entries(
        &json_entries(&object(&text_length.1), &inputs),
        &json_entries(&object(six), &inputs),
        &[UTR_STRUCTURAL],
        None,
        None,
    );
    assert!(outcome.is("red", &[]));
    assert_eq!(outcome.differing("5utr_annotation"), 1);
}

/// The comparison of two one-row outputs under the stop-boundary class.
fn boundary(reference: &str, candidate: &str) -> Outcome {
    compared(reference, candidate, &[UTR_BOUNDARY])
}

/// A made-up grid over an insertion's row on which both engines report the frameshift and
/// the candidate lists one upstream ORF more than the reference's one or two. The row's
/// HGVSc is an `ins` between two neighbouring 5' UTR positions, with and without an
/// `HGVS_OFFSET`, which a reverse-strand transcript's row prints negative and which moves
/// the point back by its size either way, a `dup` of one base or of two, an `ins` between
/// positions that are no neighbours, a `delins`, an intron offset, a coding position, a
/// substitution, a deletion, or absent; the ORF the candidate adds is a uORF whose length
/// changes and whose stop codon begins 54 bases before the coding start, the same uORF with
/// its length unchanged, an overlapping ORF there, or a uORF whose stop codon begins one base
/// further from the coding start. The row is the class's exactly where the notation places
/// the base directly after the insertion 54 bases before the coding start and the added ORF
/// is the uORF whose stop codon begins there with a changed length; it is counted once, on
/// its 5UTR_annotation.
#[test]
fn the_frameshift_at_stop_boundary_takes_an_orf_whose_stop_codon_follows_the_insertion() {
    let added: [(&str, String); 4] = [
        ("at the stop", orf_at_stop()),
        ("same length", shifted_orf(87, 36, 36)),
        (
            "not closed",
            shifted_orf_of(87, "NA", "45", "OutOfFrame_oORF", "GCCATGG"),
        ),
        ("stop elsewhere", shifted_orf(88, 36, 45)),
    ];
    let notations: [(Option<&str>, Option<i64>, Option<i64>); 18] = [
        (Some("c.-55_-54insA"), None, Some(54)),
        (Some("c.-55_-54insA"), Some(0), Some(54)),
        (Some("c.-53_-52insA"), Some(2), Some(54)),
        (Some("c.-55dup"), None, Some(54)),
        (Some("c.-56_-55dup"), None, Some(54)),
        (Some("c.-54dup"), Some(1), Some(54)),
        (Some("c.-54_-53insA"), None, Some(53)),
        (Some("c.-55_-54insA"), Some(3), Some(57)),
        (Some("c.-53_-52insA"), Some(-2), Some(54)),
        (Some("c.-54dup"), Some(-1), Some(54)),
        (Some("c.-55_-54insA"), Some(-3), Some(57)),
        (Some("c.-55_-53insA"), None, None),
        (Some("c.-55_-54delinsAA"), None, None),
        (Some("c.-55+1_-55+2insA"), None, None),
        (Some("c.1_2insA"), None, None),
        (Some("c.-55A>G"), None, None),
        (Some("c.-55del"), None, None),
        (None, None, None),
    ];
    let (mut taken, mut refused) = (0, 0);
    for (notation, offset, after) in notations {
        for listed in [vec![orf_through()], vec![orf_through(), orf_later()]] {
            for (kind, group) in &added {
                let row = |groups: &[String]| {
                    one_row(shifted_row_of(
                        groups,
                        notation,
                        offset,
                        UTR_FRAMESHIFT,
                        UTR_BOUNDARY_ROW,
                    ))
                };
                let mut written = listed.clone();
                written.push(group.clone());
                let outcome = boundary(&row(&listed), &row(&written));
                let case = format!("{notation:?} {offset:?} {} {kind}", listed.len());
                if *kind == "at the stop" && after == Some(54) {
                    assert!(outcome.is("green", &[(UTR_BOUNDARY, 1)]), "{case}");
                    assert!(
                        outcome.masks("5UTR_annotation", &[(UTR_BOUNDARY, 1)]),
                        "{case}"
                    );
                    taken += 1;
                } else {
                    assert!(outcome.is("red", &[]), "{case}");
                    refused += 1;
                }
            }
        }
    }
    // 18 notations, 2 reference lists and 4 added ORFs: the 8 notations that place the base
    // after the insertion 54 bases before the coding start, with the one ORF that ends there.
    assert_eq!((taken, refused), (8 * 2, 18 * 2 * 4 - 8 * 2));
}

/// Two ORFs the candidate adds are both held to the test; an ORF both engines list whose
/// stop codon follows the insertion refuses the row, since the reference's test admits there
/// what it leaves out; so do a group the reference alone lists, the two sides exchanged, and
/// any other term on either side. The reference's pairs in another order beside an added ORF
/// are the class's; the same groups in another order with nothing added are the pair-order
/// class's. The class reads no record and no seeded run, and on a row whose two entries
/// report the frameshift it sets aside 5UTR_annotation alone: the consequence, the plugin's
/// counts and another plugin's field on its row stay differences.
#[test]
fn the_frameshift_at_stop_boundary_refuses_the_side_the_term_and_the_field() {
    let through = orf_through();
    let reversed: Vec<&str> = through.split(':').rev().collect();
    let hash_order = reversed.join(":");
    let cases: [(Vec<String>, Vec<String>, bool); 8] = [
        (
            vec![orf_through()],
            vec![orf_through(), orf_at_stop(), orf_at_stop_too()],
            true,
        ),
        (
            vec![orf_through()],
            vec![orf_at_stop(), orf_through()],
            true,
        ),
        (
            vec![hash_order.clone()],
            vec![orf_through(), orf_at_stop()],
            true,
        ),
        (
            vec![orf_through()],
            vec![orf_through(), orf_at_stop(), orf_later()],
            false,
        ),
        (
            vec![orf_at_stop_too()],
            vec![orf_at_stop_too(), orf_at_stop()],
            false,
        ),
        (
            vec![orf_through(), orf_later()],
            vec![orf_through(), orf_at_stop()],
            false,
        ),
        (
            vec![orf_through(), orf_at_stop()],
            vec![orf_through()],
            false,
        ),
        (vec![hash_order.clone()], vec![orf_through()], false),
    ];
    for (listed, written, taken) in &cases {
        let outcome = boundary(
            &one_row(shifted_row(listed)),
            &one_row(shifted_row(written)),
        );
        assert!(
            outcome.taken_is(*taken, &[(UTR_BOUNDARY, 1)]),
            "{listed:?} {written:?}"
        );
    }
    for classes in [[PAIR_ORDER, UTR_BOUNDARY], [UTR_BOUNDARY, PAIR_ORDER]] {
        let outcome = compared(
            &one_row(shifted_row(std::slice::from_ref(&hash_order))),
            &one_row(shifted_row(&[orf_through()])),
            &classes,
        );
        assert!(outcome.is("green", &[(PAIR_ORDER, 1)]), "{classes:?}");
    }
    let reference = shifted_row(&[orf_through()]);
    let candidate = shifted_row(&[orf_through(), orf_at_stop()]);
    for term in [UTR_LOSS, UTR_STOP_LOSS, UTR_GAIN, "-"] {
        for pair in [
            (reference.replace(UTR_FRAMESHIFT, term), candidate.clone()),
            (reference.clone(), candidate.replace(UTR_FRAMESHIFT, term)),
        ] {
            let outcome = boundary(&one_row(pair.0.clone()), &one_row(pair.1));
            assert!(outcome.is("red", &[]), "{term} {}", pair.0 == reference);
        }
    }
    // A record handed or not, its reference allele empty or not, and the seeded runs change nothing.
    let inputs = lines(&[
        "21\t1000\tgrid_ins\tG\tGA\t.\t.\t.",
        "21\t1000\tgrid_ins\tGC\tGA\t.\t.\t.",
    ]);
    let (reference_text, candidate_text) = (one_row(reference.clone()), one_row(candidate.clone()));
    let cases: [(Option<&[String]>, bool); 4] = [
        (Some(&inputs[..1]), false),
        (Some(&inputs[1..]), false),
        (None, false),
        (None, true),
    ];
    for (lines, swept) in cases {
        let outcome = compared_entries(
            &default_entries(&reference_text),
            &default_entries(&candidate_text),
            &[UTR_BOUNDARY],
            lines,
            swept.then(|| seed_sweep(&[&reference_text])),
        );
        assert!(
            outcome.is("green", &[(UTR_BOUNDARY, 1)]),
            "{lines:?} {swept}"
        );
    }
    let e = &default_entries(&reference_text)[0];
    let a = &default_entries(&candidate_text)[0];
    let (listed, written) = (&e.fields["5UTR_annotation"], &a.fields["5UTR_annotation"]);
    let accepted = |name: &str, a: Option<&Entry>| {
        accepts(UTR_BOUNDARY, e, name, listed, written, a, None, None)
    };
    assert!(accepted("5UTR_annotation", Some(a)));
    assert!(accepted("5utr_annotation", Some(a)));
    assert!(!accepted("5UTR_annotation", None));
    for name in [
        "5UTR_consequence",
        "Existing_uORFs",
        "Existing_InFrame_oORFs",
        "Existing_OutOfFrame_oORFs",
        "LOEUF",
        "HGVSc",
        "Extra",
    ] {
        assert!(
            !accepts(UTR_BOUNDARY, e, name, "1", "2", Some(a), None, None),
            "{name}"
        );
    }
    let moved = candidate.replace("Existing_uORFs=1", "Existing_uORFs=2") + ";LOEUF=0.4";
    let outcome = boundary(&one_row(reference + ";LOEUF=0.3"), &one_row(moved));
    assert!(outcome.is("red", &[(UTR_BOUNDARY, 1)]));
    assert_eq!(outcome.differing("Existing_uORFs"), 1);
    assert_eq!(outcome.differing("LOEUF"), 1);
    assert!(outcome.masks("5UTR_annotation", &[(UTR_BOUNDARY, 1)]));
}

/// A made-up identifier for an insertion of four bases between cDNA positions 493 and 494 of
/// a reverse-strand transcript.
const UTR_KEPT_ROW: Row = ("kept_ins", "21:1000-1001", "TTAA", "ENST1", "493-494");

/// A uORF that starts 314 bases before the coding start and keeps its length of 18 bases:
/// its stop codon begins 299 bases before the coding start.
const ORF_KEPT: &str = "Evidence=False:KozakContext=TTAATGG:KozakStrength=Moderate:alt_type=uORF:alt_type_length=18:ref_StartDistanceToCDS=314:ref_type=uORF:ref_type_length=18";

/// A counted row of `UTR_KEPT_ROW` with the two fields `own`, and the HGVSc body
/// `notation` and the `HGVS_OFFSET` where they are given.
fn kept_row(own: &Values, notation: Option<&str>, offset: Option<i64>) -> String {
    match notation {
        Some(notation) => written_row(UTR_KEPT_ROW, own, "-1", notation, offset),
        None => stranded_row(UTR_KEPT_ROW, own, Some("-1")),
    }
}

/// The frameshift of `groups`.
fn shifted_list(groups: &[&str]) -> Values {
    values(UTR_FRAMESHIFT, &groups.join(","))
}

/// A grid over an insertion's row on which the reference lists one upstream ORF more than
/// the candidate, which reports nothing or the frameshift of the ORFs both list. The ORF
/// the reference adds is a uORF that keeps its length and whose stop codon begins 299 bases
/// before the coding start, the same uORF with a changed length, the same ORF opened into
/// the coding sequence, or an overlapping ORF; the row's HGVSc and `HGVS_OFFSET` place the
/// base after the insertion as the record writes it and as the notation names it. The row is
/// the class's exactly where the added ORF keeps its length and its stop codon begins at or
/// after the record's base and no more than one base beyond the notation's; it is counted
/// once, with both fields where the candidate reports nothing and with its
/// 5UTR_annotation alone where both report the frameshift.
#[test]
fn the_frameshift_at_stop_boundary_takes_an_orf_that_keeps_its_length_at_its_stop_codon() {
    let dropped: [(&str, String); 4] = [
        ("kept", ORF_KEPT.to_string()),
        (
            "changed length",
            ORF_KEPT.replace("alt_type_length=18", "alt_type_length=30"),
        ),
        (
            "opened",
            ORF_KEPT.replace("alt_type=uORF", "alt_type=OutOfFrame_oORF"),
        ),
        (
            "overlapping",
            ORF_KEPT.replace("ref_type=uORF", "ref_type=InFrame_oORF"),
        ),
    ];
    // The record's base and the notation's, in bases before the coding start: 300 and 300,
    // 301 and 300, 302 and 300, 299 and 299, 300 and 300; then 301 and 301, 304 and 301,
    // 298 and 298, and none.
    let notations: [(Option<&str>, Option<i64>, bool); 9] = [
        (Some("c.-301_-300insTTAA"), None, true),
        (Some("c.-301_-300insTTAG"), Some(-1), true),
        (Some("c.-301_-300insTTAG"), Some(2), true),
        (Some("c.-300_-299insTTAA"), None, true),
        (Some("c.-301dup"), None, true),
        (Some("c.-302_-301insTTAA"), None, false),
        (Some("c.-302_-301insTTAA"), Some(3), false),
        (Some("c.-299_-298insTTAA"), None, false),
        (None, None, false),
    ];
    let through = orf_through();
    let (mut taken, mut refused) = (0, 0);
    for (notation, offset, in_window) in notations {
        for (kind, group) in &dropped {
            let row = |own: &Values| one_row(kept_row(own, notation, offset));
            let case = format!("{notation:?} {offset:?} {kind}");
            let fits = *kind == "kept" && in_window;
            let outcome = boundary(&row(&shifted_list(&[group.as_str()])), &row(&nothing()));
            assert!(outcome.taken_is(fits, &[(UTR_BOUNDARY, 1)]), "{case}");
            if fits {
                for name in ["5UTR_consequence", "5UTR_annotation"] {
                    assert!(outcome.masks(name, &[(UTR_BOUNDARY, 1)]), "{case}");
                }
            }
            let outcome = boundary(
                &row(&shifted_list(&[through.as_str(), group.as_str()])),
                &row(&shifted_list(&[through.as_str()])),
            );
            assert!(outcome.taken_is(fits, &[(UTR_BOUNDARY, 1)]), "{case}");
            if fits {
                assert!(outcome.masks("5UTR_annotation", &[(UTR_BOUNDARY, 1)]));
                assert!(outcome.masks("5UTR_consequence", &[]), "{case}");
                taken += 1;
            } else {
                refused += 1;
            }
        }
    }
    // 9 notations and 4 added ORFs: the 5 notations whose window holds the stop codon, with
    // the one ORF that keeps its length.
    assert_eq!((taken, refused), (5, 9 * 4 - 5));
}

/// A shorter list on the candidate's side stays a difference where one of two ORFs the
/// reference adds changes its length, where the candidate adds an ORF of its own, where
/// its row carries no counts, another term, or an annotation beside no term, and with the
/// two sides exchanged. Two equal groups of which the candidate lists one are the class's.
/// The consequence is the class's field only where one entry reports nothing.
#[test]
fn the_frameshift_at_stop_boundary_refuses_a_shorter_list_outside_the_kept_length_shape() {
    let notation = Some("c.-301_-300insTTAA");
    let row = |own: &Values| one_row(kept_row(own, notation, None));
    let changed = ORF_KEPT.replace("alt_type_length=18", "alt_type_length=30");
    let through = orf_through();
    let kept = shifted_list(&[ORF_KEPT]);
    let cases: [(&str, Values, Values, bool); 7] = [
        ("one ORF kept", kept.clone(), nothing(), true),
        (
            "two equal groups, one listed",
            shifted_list(&[ORF_KEPT, ORF_KEPT]),
            kept.clone(),
            true,
        ),
        (
            "one of two changes its length",
            shifted_list(&[ORF_KEPT, changed.as_str()]),
            nothing(),
            false,
        ),
        (
            "an ORF on each side alone",
            shifted_list(&[through.as_str(), ORF_KEPT]),
            shifted_list(&[through.as_str(), changed.as_str()]),
            false,
        ),
        ("another term", kept.clone(), utr_lost_start(), false),
        (
            "an annotation beside no term",
            kept.clone(),
            values("-", ORF_KEPT),
            false,
        ),
        ("the sides exchanged", nothing(), kept.clone(), false),
    ];
    for (what, reference, candidate, taken) in &cases {
        let outcome = boundary(&row(reference), &row(candidate));
        assert!(outcome.taken_is(*taken, &[(UTR_BOUNDARY, 1)]), "{what}");
    }
    let uncounted =
        one_row(utr_row(UTR_KEPT_ROW, &nothing()) + ";HGVSc=ENST1.1:c.-301_-300insTTAA");
    assert!(boundary(&row(&kept), &uncounted).is("red", &[]));
    let e = &default_entries(&row(&kept))[0];
    let nothing_written = &default_entries(&row(&nothing()))[0];
    let fewer = &default_entries(&row(&shifted_list(&[ORF_KEPT, ORF_KEPT])))[0];
    let accepted = |e: &Entry, name: &str, expected_value: &str, actual_value: &str, a: &Entry| {
        accepts(
            UTR_BOUNDARY,
            e,
            name,
            expected_value,
            actual_value,
            Some(a),
            None,
            None,
        )
    };
    for name in ["5UTR_consequence", "5utr_consequence", "5UTR_annotation"] {
        assert!(accepted(e, name, "x", "-", nothing_written), "{name}");
    }
    assert!(accepted(fewer, "5UTR_annotation", "x", "y", e));
    assert!(!accepted(fewer, "5UTR_consequence", "x", "y", e));
    assert!(!accepted(e, "Existing_uORFs", "x", "-", nothing_written));
}

/// A made-up identifier for an insertion of four bases between cDNA positions 288 and 289 of
/// a reverse-strand transcript, written `c.-236_-235insCTGT`.
const UTR_EMPTY_ROW: Row = ("empty_ins", "21:1000-1001", "ACAG", "ENST1", "288-289");

/// A uORF that starts 241 bases before the coding start and whose length changes from 9
/// bases: its stop codon begins 235 bases before the coding start.
const ORF_AT_POINT: &str = "Evidence=False:KozakContext=CCTATGT:KozakStrength=Weak:alt_type=uORF:alt_type_length=111:ref_StartDistanceToCDS=241:ref_type=uORF:ref_type_length=9";

/// The reference reports nothing on an insertion written directly before a stop codon and
/// the candidate the frameshift of the uORF that ends there: the row is the class's, with
/// both fields, where every ORF the candidate lists changes its length and has its stop
/// codon's first base directly after the insertion. It stays a difference where the
/// notation places the insertion elsewhere, where an ORF keeps its length or ends
/// elsewhere, and beside a reference annotation that stands beside no term.
#[test]
fn the_frameshift_at_stop_boundary_takes_a_list_the_reference_leaves_empty() {
    let row = |own: &Values, notation: &str| {
        one_row(written_row(UTR_EMPTY_ROW, own, "-1", notation, None))
    };
    let at_point = "c.-236_-235insCTGT";
    let listed = shifted_list(&[ORF_AT_POINT]);
    let outcome = boundary(&row(&nothing(), at_point), &row(&listed, at_point));
    assert!(outcome.is("green", &[(UTR_BOUNDARY, 1)]));
    for name in ["5UTR_consequence", "5UTR_annotation"] {
        assert!(outcome.masks(name, &[(UTR_BOUNDARY, 1)]), "{name}");
    }
    let unnamed = compared(&row(&nothing(), at_point), &row(&listed, at_point), &[]);
    assert_eq!(unnamed.verdict(), "red");
    let same_length = ORF_AT_POINT.replace("alt_type_length=111", "alt_type_length=9");
    let elsewhere = ORF_AT_POINT.replace("ref_type_length=9", "ref_type_length=12");
    let cases: [(&str, Values, Values, &str); 5] = [
        (
            "the insertion one base on",
            nothing(),
            listed.clone(),
            "c.-235_-234insCTGT",
        ),
        (
            "an ORF that keeps its length",
            nothing(),
            shifted_list(&[same_length.as_str()]),
            at_point,
        ),
        (
            "a second ORF that ends elsewhere",
            nothing(),
            shifted_list(&[ORF_AT_POINT, elsewhere.as_str()]),
            at_point,
        ),
        (
            "a reference annotation beside no term",
            values("-", ORF_AT_POINT),
            listed.clone(),
            at_point,
        ),
        (
            "the reference's lost start",
            utr_lost_start(),
            listed.clone(),
            at_point,
        ),
    ];
    for (what, reference, candidate, notation) in &cases {
        let outcome = boundary(&row(reference, notation), &row(candidate, notation));
        assert!(outcome.is("red", &[]), "{what}");
    }
}

/// The reference's and the candidate's outputs of an insertion written after a start
/// codon's second base that opens with the codon's last base, and the row itself.
fn boundary_after_a_start_codon() -> (Row<'static>, String, String) {
    let row: Row = ("grid_ins", "21:1000-1001", "G", "ENST1", "146-147");
    let text = |groups: &[String]| {
        one_row(shifted_row_of(
            groups,
            Some("c.-55_-54insG"),
            None,
            UTR_FRAMESHIFT,
            row,
        ))
    };
    (
        row,
        text(&[orf_through()]),
        text(&[orf_through(), orf_at_stop()]),
    )
}

/// A row two classes hold is counted under the first in the order a row is put to the
/// classes. On an insertion written after a start codon's second base that opens with the
/// codon's last base, the seeded runs write a gained start at the codon's own index beside
/// the reference's frameshift: the candidate's longer list is the remade-in-place class's
/// with those runs, whichever class is named first, and this class's without them or named
/// alone.
#[test]
fn the_frameshift_at_stop_boundary_follows_a_start_codon_left_in_place() {
    let (row, reference, candidate) = boundary_after_a_start_codon();
    let seeded = [reference.clone(), one_row(counted_row(row, &gained(144)))];
    let cases: [(bool, Classes, Applied); 4] = [
        (true, &[UTR_REMADE, UTR_BOUNDARY], &[(UTR_REMADE, 1)]),
        (true, &[UTR_BOUNDARY, UTR_REMADE], &[(UTR_REMADE, 1)]),
        (true, &[UTR_BOUNDARY], &[(UTR_BOUNDARY, 1)]),
        (false, &[UTR_REMADE, UTR_BOUNDARY], &[(UTR_BOUNDARY, 1)]),
    ];
    for (swept, classes, applied) in cases {
        let outcome = compared_swept(
            &reference,
            &candidate,
            swept.then_some(&seeded[..]),
            classes,
        );
        assert!(outcome.is("green", applied), "{swept} {classes:?}");
    }
}

/// The command line names the stop-boundary class like any other.
#[test]
fn the_command_line_names_the_frameshift_at_stop_boundary_class() {
    let (_, reference, candidate) = boundary_after_a_start_codon();
    let directory = planted_files(
        "utrannotator-boundary",
        &[("reference", &reference), ("candidate", &candidate)],
    );
    assert_eq!(exit_code(&directory, &[]), 1);
    assert_eq!(exit_code(&directory, &["--masks", UTR_BOUNDARY]), 0);
}

/// The tab format carries the HGVSc and the `HGVS_OFFSET` as columns, the VCF format as
/// CSQ fields with `&` between the annotation's groups, and the JSON format as `hgvsc`
/// and `hgvs_offset` beside an object of numbered groups. In each the row is the class's
/// where the notation places the base after the insertion at the added ORF's stop codon,
/// and stays a difference where it places it elsewhere or the format carries no HGVSc.
#[test]
fn the_frameshift_at_stop_boundary_reads_the_tab_vcf_and_json_formats() {
    let cases: [(Option<&str>, &str, bool); 5] = [
        (Some("c.-55_-54insA"), "-", true),
        (Some("c.-53_-52insA"), "2", true),
        (Some("c.-53_-52insA"), "-2", true),
        (Some("c.-54_-53insA"), "-", false),
        (None, "-", false),
    ];
    let listed = vec![orf_through()];
    let written = vec![orf_through(), orf_at_stop()];
    let under_class = |reference: &[Entry], candidate: &[Entry]| {
        compared_entries(reference, candidate, &[UTR_BOUNDARY], None, None)
    };
    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tSTRAND\tHGVSc\tHGVS_OFFSET\t5UTR_consequence\t5UTR_annotation\tExisting_uORFs\n";
    let tab = |notation: Option<&str>, offset: &str, groups: &[String]| {
        let hgvsc = notation.map_or("-".to_string(), |notation| format!("ENST1.1:{notation}"));
        format!(
            "{header}grid_ins\t21:1000-1001\tA\tENSG1\tENST1\tTranscript\t5_prime_UTR_variant\t146-147\t1\t{hgvsc}\t{offset}\t{UTR_FRAMESHIFT}\t{}\t1\n",
            groups.join(",")
        )
    };
    for (notation, offset, taken) in cases {
        let outcome = under_class(
            &tab_entries(&tab(notation, offset, &listed)),
            &tab_entries(&tab(notation, offset, &written)),
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_BOUNDARY, 1)]),
            "{notation:?} {offset}"
        );
    }

    let inputs = lines(&["21\t1000\tgrid_ins\tG\tGA\t.\t.\t."]);
    for (notation, offset, taken) in cases {
        let fields = if notation.is_some() {
            "HGVSc|HGVS_OFFSET|"
        } else {
            ""
        };
        let held = notation.map_or(String::new(), |notation| {
            format!(
                "ENST1.1:{notation}|{}|",
                if offset == "-" { "" } else { offset }
            )
        });
        let csq = csq_header(&format!(
            "Allele|Consequence|Feature|Feature_type|cDNA_position|{fields}5UTR_consequence|5UTR_annotation|Existing_uORFs"
        ));
        let line = format!(
            "{csq}{}CSQ=A|5_prime_UTR_variant|ENST1|Transcript|146-147|{held}{UTR_FRAMESHIFT}|{}|1\n",
            record_columns(&inputs[0]),
            orf_through()
        );
        let longer = line.replace(
            &orf_through(),
            &format!("{}&{}", orf_through(), orf_at_stop()),
        );
        let outcome = under_class(&vcf_entries(&line, &inputs), &vcf_entries(&longer, &inputs));
        assert!(
            outcome.taken_is(taken, &[(UTR_BOUNDARY, 1)]),
            "{notation:?} {offset}"
        );
    }

    let object = |notation: Option<&str>, offset: &str, groups: &[String]| {
        let mut entry = format!(
            r#""variant_allele": "A", "transcript_id": "ENST1", "consequence_terms": ["5_prime_UTR_variant"], "cdna_start": 147, "cdna_end": 146, "strand": 1, "existing_uorfs": 1, "5utr_consequence": "{UTR_FRAMESHIFT}", "5utr_annotation": {}"#,
            json_groups(groups)
        );
        if let Some(notation) = notation {
            entry.push_str(&format!(r#", "hgvsc": "ENST1.1:{notation}""#));
        }
        if offset != "-" {
            entry.push_str(&format!(r#", "hgvs_offset": {offset}"#));
        }
        json_line(&inputs[0], "grid_ins", "-/A", &entry)
    };
    for (notation, offset, taken) in cases {
        let outcome = under_class(
            &json_entries(&object(notation, offset, &listed), &inputs),
            &json_entries(&object(notation, offset, &written), &inputs),
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_BOUNDARY, 1)]),
            "{notation:?} {offset}"
        );
    }
}

/// A uORF that runs into the coding sequence once shifted: no stop codon ends it in the
/// 5' UTR.
const ORF_OPENED: &str = "Evidence=False:KozakContext=GCCATGG:KozakStrength=Strong:alt_type=OutOfFrame_oORF:alt_type_length=NA:ref_StartDistanceToCDS=87:ref_type=uORF:ref_type_length=36";

/// The comparison of two one-row outputs under the lost-stop class.
fn before_stop(reference: &str, candidate: &str) -> Outcome {
    compared(reference, candidate, &[UTR_BEFORE_STOP])
}

/// A counted row of the frameshift with `groups` and no HGVSc.
fn listed_row(groups: &[String], row: Row) -> String {
    shifted_row_of(groups, None, None, UTR_FRAMESHIFT, row)
}

/// A made-up grid over an insertion's row on which the reference reports the frameshift
/// and the candidate the lost stop with one group. The reference lists a uORF that starts 87
/// bases before the coding start and is 45 bases long once shifted, alone, beside an
/// overlapping ORF or beside a second uORF; the allele is two bases long or eight; the
/// candidate's group carries that uORF's Kozak pairs with the next stop codon where its
/// shifted length ends it (the start's distance plus the allele's length, less the length,
/// plus three), one base off, or none, or another Kozak context or strength, or the
/// overlapping ORF's pairs. The row is the class's exactly where the group carries the
/// listed uORF's pairs and that stop codon; it is counted once, with both fields. Where the
/// listed uORF runs into the coding sequence once shifted, the class takes the group that
/// names no alternative stop codon and refuses one that names any.
#[test]
fn the_lost_stop_written_before_stop_codon_takes_the_lost_stop_of_a_listed_orf() {
    let (mut taken, mut refused) = (0, 0);
    for allele in ["AT", "ATCCTAGG"] {
        let row: Row = (
            UTR_BEFORE_STOP_ROW.0,
            UTR_BEFORE_STOP_ROW.1,
            allele,
            UTR_BEFORE_STOP_ROW.3,
            UTR_BEFORE_STOP_ROW.4,
        );
        let ends = 87 + allele.len() as i64 - 45 + 3;
        let groups: [(&str, String); 6] = [
            ("the listed uORF", lost_stop_group(Some(ends))),
            ("one base off", lost_stop_group(Some(ends + 1))),
            ("no alternative stop", lost_stop_group(None)),
            (
                "another context",
                lost_stop_group_of(Some(ends), "AAAATGA", "Strong"),
            ),
            (
                "another strength",
                lost_stop_group_of(Some(ends), "GCCATGG", "Weak"),
            ),
            (
                "the overlapping ORF",
                lost_stop_group_of(Some(ends), "CCCATGC", "Strong"),
            ),
        ];
        for listed in [
            vec![orf_at_stop()],
            vec![orf_through(), orf_at_stop()],
            vec![orf_at_stop(), orf_later()],
        ] {
            for (kind, group) in &groups {
                let outcome = before_stop(
                    &one_row(listed_row(&listed, row)),
                    &one_row(lost_stop_row_of(
                        std::slice::from_ref(group),
                        row,
                        UTR_COUNTS,
                    )),
                );
                let case = format!("{allele} {} {kind}", listed.len());
                if *kind == "the listed uORF" {
                    assert!(outcome.is("green", &[(UTR_BEFORE_STOP, 1)]), "{case}");
                    for name in ["5UTR_consequence", "5UTR_annotation"] {
                        assert!(outcome.masks(name, &[(UTR_BEFORE_STOP, 1)]), "{case}");
                    }
                    taken += 1;
                } else {
                    assert!(outcome.is("red", &[]), "{case}");
                    refused += 1;
                }
            }
        }
    }
    // 2 alleles, 3 reference lists and 6 groups, one of them the listed uORF's.
    assert_eq!((taken, refused), (2 * 3, 2 * 3 * 5));
    for (group, accepted) in [
        (lost_stop_group(None), true),
        (lost_stop_group(Some(47)), false),
    ] {
        let outcome = before_stop(
            &one_row(listed_row(&[ORF_OPENED.to_string()], UTR_BEFORE_STOP_ROW)),
            &one_row(lost_stop_row(std::slice::from_ref(&group))),
        );
        assert!(
            outcome.taken_is(accepted, &[(UTR_BEFORE_STOP, 1)]),
            "{group}"
        );
    }
}

/// Two groups of the candidate's are each held to a different listed uORF; two that carry
/// one uORF's pairs, and a group whose pairs two listed uORFs share, refuse the row. So do
/// the two sides exchanged, any other term on either side, cDNA positions that are not an
/// insertion's two flanking bases, a deletion's or a structural allele's row, a record whose
/// reference allele is not empty, spelt in the first column or handed with the row, and a
/// candidate's row without the plugin's counts. The plugin's counts and another plugin's
/// field on the class's row stay differences.
#[test]
fn the_lost_stop_written_before_stop_codon_refuses_the_side_the_kind_and_the_field() {
    let at_stop = lost_stop_group(Some(47));
    let too = lost_stop_group_of(Some(47), "AAAATGG", "Strong");
    let reference = listed_row(&[orf_through(), orf_at_stop()], UTR_BEFORE_STOP_ROW);
    let candidate = lost_stop_row(std::slice::from_ref(&at_stop));
    let cases: [(Vec<String>, Vec<String>, bool); 5] = [
        (
            vec![orf_at_stop(), orf_at_stop_too()],
            vec![at_stop.clone(), too.clone()],
            true,
        ),
        (vec![orf_at_stop(), orf_at_stop_too()], vec![too], true),
        (
            vec![orf_at_stop(), orf_at_stop_too()],
            vec![at_stop.clone(), at_stop.clone()],
            false,
        ),
        (
            vec![orf_at_stop(), shifted_orf(60, 30, 18)],
            vec![at_stop.clone()],
            false,
        ),
        (vec![orf_through()], vec![at_stop.clone()], false),
    ];
    for (listed, lost, taken) in &cases {
        let outcome = before_stop(
            &one_row(listed_row(listed, UTR_BEFORE_STOP_ROW)),
            &one_row(lost_stop_row(lost)),
        );
        assert!(
            outcome.taken_is(*taken, &[(UTR_BEFORE_STOP, 1)]),
            "{} {lost:?}",
            listed.len()
        );
    }
    let exchanged = before_stop(&one_row(candidate.clone()), &one_row(reference.clone()));
    assert!(exchanged.is("red", &[]));
    for term in [UTR_LOSS, UTR_GAIN, UTR_STOP_LOSS, "-"] {
        let other = if term != UTR_STOP_LOSS {
            term
        } else {
            UTR_FRAMESHIFT
        };
        for pair in [
            (reference.replace(UTR_FRAMESHIFT, term), candidate.clone()),
            (reference.clone(), candidate.replace(UTR_STOP_LOSS, other)),
        ] {
            let outcome = before_stop(&one_row(pair.0.clone()), &one_row(pair.1));
            assert!(outcome.is("red", &[]), "{term} {}", pair.0 == reference);
        }
    }
    let stop_at = |cdna: &'static str| -> Row {
        (
            UTR_BEFORE_STOP_ROW.0,
            UTR_BEFORE_STOP_ROW.1,
            UTR_BEFORE_STOP_ROW.2,
            UTR_BEFORE_STOP_ROW.3,
            cdna,
        )
    };
    let rows: [(Row, &str); 7] = [
        (stop_at("146"), UTR_COUNTS),
        (stop_at("146-148"), UTR_COUNTS),
        (stop_at("?-147"), UTR_COUNTS),
        (("del", "21:1000-1001", "-", "ENST1", "146-147"), UTR_COUNTS),
        (
            ("sv", "21:1000-1001", "<INS>", "ENST1", "146-147"),
            UTR_COUNTS,
        ),
        (
            ("21_1000_GC/AT", "21:1000-1001", "AT", "ENST1", "146-147"),
            UTR_COUNTS,
        ),
        (UTR_BEFORE_STOP_ROW, ""),
    ];
    for (row, counts) in rows {
        let outcome = before_stop(
            &one_row(listed_row(&[orf_at_stop()], row)),
            &one_row(lost_stop_row_of(
                std::slice::from_ref(&at_stop),
                row,
                counts,
            )),
        );
        assert!(outcome.is("red", &[]), "{row:?} {counts:?}");
    }
    // Under an identifier the input's line says whether the reference allele is empty; without
    // it a two-base Location is read by its positions and its allele alone.
    let inputs = lines(&[
        "21\t1000\tgrid_ins\tG\tGAT\t.\t.\t.",
        "21\t500\tmnv\tAC\tAT\t.\t.\t.",
    ]);
    let rows: [(Row, bool, bool); 2] = [
        (UTR_BEFORE_STOP_ROW, true, true),
        (("mnv", "21:500-501", "AT", "ENST1", "146-147"), false, true),
    ];
    let listed: Vec<String> = rows
        .iter()
        .map(|(row, _, _)| listed_row(&[orf_at_stop()], *row))
        .collect();
    let lost: Vec<String> = rows
        .iter()
        .map(|(row, _, _)| lost_stop_row_of(std::slice::from_ref(&at_stop), *row, UTR_COUNTS))
        .collect();
    let listed = default_entries(&default_text(&listed));
    let lost = default_entries(&default_text(&lost));
    for (at, (row, handed, alone)) in rows.iter().enumerate() {
        let found = accepts_with_and_without_the_record(
            UTR_BEFORE_STOP,
            &listed,
            &inputs,
            at,
            &lost[at],
            "1",
            "2",
        );
        assert_eq!(found, [*handed, *handed, *alone, *alone], "{row:?}");
    }
    let (reference_text, candidate_text) = (one_row(reference.clone()), one_row(candidate.clone()));
    let given: [Option<&[String]>; 3] = [Some(&inputs[..1]), Some(&inputs[1..]), None];
    for lines in given {
        let outcome = compared_with(&reference_text, &candidate_text, &[UTR_BEFORE_STOP], lines);
        assert!(outcome.is("green", &[(UTR_BEFORE_STOP, 1)]), "{lines:?}");
    }
    let e = &default_entries(&reference_text)[0];
    let a = &default_entries(&candidate_text)[0];
    for name in ["5UTR_consequence", "5UTR_annotation", "5utr_annotation"] {
        assert!(
            accepts(UTR_BEFORE_STOP, e, name, "1", "2", Some(a), None, None),
            "{name}"
        );
        assert!(
            !accepts(UTR_BEFORE_STOP, e, name, "1", "2", None, None, None),
            "{name}"
        );
    }
    for name in [
        "Existing_uORFs",
        "Existing_InFrame_oORFs",
        "Existing_OutOfFrame_oORFs",
        "LOEUF",
        "Consequence",
        "Extra",
    ] {
        assert!(
            !accepts(UTR_BEFORE_STOP, e, name, "1", "2", Some(a), None, None),
            "{name}"
        );
    }
    let moved = candidate.replace("Existing_uORFs=1", "Existing_uORFs=2") + ";LOEUF=0.4";
    let outcome = before_stop(&one_row(reference + ";LOEUF=0.3"), &one_row(moved));
    assert!(outcome.is("red", &[(UTR_BEFORE_STOP, 1)]));
    assert_eq!(outcome.differing("Existing_uORFs"), 1);
    assert_eq!(outcome.differing("LOEUF"), 1);
    assert!(outcome.masks("5UTR_consequence", &[(UTR_BEFORE_STOP, 1)]));
}

/// An insertion's row under a made-up identifier with its transcript's strand, its HGVSc
/// body and its `HGVS_OFFSET`, and what the reference and the candidate write on it.
struct Written {
    row: Row<'static>,
    strand: &'static str,
    notation: &'static str,
    offset: Option<i64>,
    reference: Values,
    candidate: Values,
}

impl Written {
    fn text(&self, own: &Values) -> String {
        written_row(self.row, own, self.strand, self.notation, self.offset)
    }

    /// The reference's and the candidate's one-row outputs.
    fn texts(&self) -> (String, String) {
        (
            one_row(self.text(&self.reference)),
            one_row(self.text(&self.candidate)),
        )
    }
}

/// Three insertions that can be written one base further on, inside a stop codon: the
/// reference reports nothing, the lost start or the frameshift of another ORF, and the
/// candidate the lost stop.
fn written_further_on() -> [Written; 3] {
    [
        Written {
            row: ("further_one", "21:2000-2001", "AGAA", "ENST1", "625-626"),
            strand: "-1",
            notation: "c.-137_-136insTCTT",
            offset: Some(-1),
            reference: nothing(),
            candidate: values(
                UTR_STOP_LOSS,
                "AltStop=True:AltStopDistanceToCDS=81:Evidence=False:FrameWithCDS=inFrame:KozakContext=AGAATGT:KozakStrength=Moderate",
            ),
        },
        Written {
            row: ("further_two", "21:3000-3001", "TCTCT", "ENST1", "361-362"),
            strand: "1",
            notation: "c.-100_-99insCTCTT",
            offset: Some(1),
            reference: values(
                UTR_LOSS,
                "CapDistanceToStart=360:DistanceToCDS=101:DistanceToStop=84:Evidence=False:KozakContext=GCCATGA:KozakStrength=Moderate:type=uORF",
            ),
            candidate: values(
                UTR_STOP_LOSS,
                "AltStop=True:AltStopDistanceToCDS=93:Evidence=False:FrameWithCDS=inFrame:KozakContext=TACATGG:KozakStrength=Moderate",
            ),
        },
        Written {
            row: ("further_three", "21:4000-4001", "AAACA", "ENST1", "67-68"),
            strand: "-1",
            notation: "c.-108_-107insGTTTT",
            offset: Some(-1),
            reference: values(
                UTR_FRAMESHIFT,
                "Evidence=False:KozakContext=ACAATGT:KozakStrength=Moderate:alt_type=uORF:alt_type_length=27:ref_StartDistanceToCDS=127:ref_type=uORF:ref_type_length=36",
            ),
            candidate: values(
                UTR_STOP_LOSS,
                "AltStop=True:AltStopDistanceToCDS=74:Evidence=False:FrameWithCDS=outOfFrame:KozakContext=TTGATGT:KozakStrength=Weak",
            ),
        },
    ]
}

/// The candidate reports the lost stop of an ORF on an insertion written at or before the
/// stop codon's first base, where the reference reads no stop codon: the row is the
/// class's, with both fields, beside a reference that reports nothing, the lost start, the
/// gained start or the frameshift of other ORFs, where the allele's length is no multiple
/// of three and the row's `HGVS_OFFSET` is not zero. Seeded runs that write the
/// reference's own terms on the row change nothing.
#[test]
fn the_lost_stop_written_before_stop_codon_takes_a_lost_stop_the_reference_does_not_read() {
    let rows = written_further_on();
    let text = |candidate: bool| {
        let rows: Vec<String> = rows
            .iter()
            .map(|row| {
                row.text(if candidate {
                    &row.candidate
                } else {
                    &row.reference
                })
            })
            .collect();
        default_text(&rows)
    };
    let (reference, candidate) = (text(false), text(true));
    let own_terms = [reference.as_str()];
    let sweeps: [Option<&[&str]>; 2] = [None, Some(&own_terms)];
    for sweep in sweeps {
        let outcome = compared_swept(&reference, &candidate, sweep, &[UTR_BEFORE_STOP]);
        assert!(outcome.is("green", &[(UTR_BEFORE_STOP, 3)]));
        for name in ["5UTR_consequence", "5UTR_annotation"] {
            assert!(outcome.masks(name, &[(UTR_BEFORE_STOP, 3)]), "{name}");
        }
    }
    let unnamed = compared(&reference, &candidate, &[]);
    assert_eq!(unnamed.verdict(), "red");
    let gained_start = Written {
        reference: gained(624),
        ..written_further_on().into_iter().next().expect("a row")
    };
    let (reference, candidate) = gained_start.texts();
    assert!(before_stop(&reference, &candidate).is("green", &[(UTR_BEFORE_STOP, 1)]));
}

/// A lost stop the candidate alone reports stays a difference where a seeded run writes the
/// lost stop on the row, where the allele's length is a multiple of three, where the row
/// carries no `HGVS_OFFSET` or a zero, beside a reference that reports the gained stop or
/// an annotation beside no term, beside a frameshift whose list carries the Kozak pairs of
/// the candidate's ORF with another stop codon, where a group of the candidate's states no
/// `KozakStrength`, on a candidate's row without the plugin's counts, and under a record
/// whose reference allele is not empty.
#[test]
fn the_lost_stop_written_before_stop_codon_refuses_a_lost_stop_a_run_writes_or_no_later_writing_explains(
) {
    let [first, _, third] = written_further_on();
    let (reference, candidate) = first.texts();
    assert!(before_stop(&reference, &candidate).is("green", &[(UTR_BEFORE_STOP, 1)]));
    let lost_in_a_run = [reference.as_str(), candidate.as_str()];
    let outcome = compared_swept(
        &reference,
        &candidate,
        Some(&lost_in_a_run[..]),
        &[UTR_BEFORE_STOP],
    );
    assert!(outcome.is("red", &[]));
    let rows: [(&str, Written); 8] = [
        (
            "an allele of three bases",
            Written {
                row: ("further_one", "21:2000-2001", "AGA", "ENST1", "625-626"),
                ..written_further_on().into_iter().next().expect("a row")
            },
        ),
        (
            "an allele of six bases",
            Written {
                row: ("further_one", "21:2000-2001", "AGAAGA", "ENST1", "625-626"),
                ..written_further_on().into_iter().next().expect("a row")
            },
        ),
        (
            "no offset",
            Written {
                offset: None,
                ..written_further_on().into_iter().next().expect("a row")
            },
        ),
        (
            "an offset of zero",
            Written {
                offset: Some(0),
                ..written_further_on().into_iter().next().expect("a row")
            },
        ),
        (
            "the reference's gained stop",
            Written {
                reference: values(
                    "5_prime_UTR_uORF_stop_codon_gain_variant",
                    "KozakStrength=Weak:ref_type=uORF",
                ),
                ..written_further_on().into_iter().next().expect("a row")
            },
        ),
        (
            "a reference annotation beside no term",
            Written {
                reference: values("-", "KozakStrength=Weak:ref_type=uORF"),
                ..written_further_on().into_iter().next().expect("a row")
            },
        ),
        (
            "a listed ORF of the candidate's Kozak pairs",
            Written {
                reference: values(
                    UTR_FRAMESHIFT,
                    &third
                        .reference
                        .1
                        .replace("KozakContext=ACAATGT", "KozakContext=TTGATGT")
                        .replace("KozakStrength=Moderate", "KozakStrength=Weak"),
                ),
                ..written_further_on().into_iter().nth(2).expect("a row")
            },
        ),
        (
            "a group of the candidate's without a KozakStrength",
            Written {
                candidate: values(
                    UTR_STOP_LOSS,
                    &format!(
                        "{},AltStop=False:AltStopDistanceToCDS=NA:KozakContext=TTGATGT",
                        first.candidate.1
                    ),
                ),
                ..written_further_on().into_iter().next().expect("a row")
            },
        ),
    ];
    for (what, row) in &rows {
        let (reference, candidate) = row.texts();
        assert!(before_stop(&reference, &candidate).is("red", &[]), "{what}");
    }
    let uncounted = one_row(
        utr_row(first.row, &first.candidate).replacen("STRAND=1;", "STRAND=-1;", 1)
            + ";HGVSc=ENST1.1:c.-137_-136insTCTT;HGVS_OFFSET=-1",
    );
    assert!(before_stop(&reference, &uncounted).is("red", &[]));
    // The input's line says whether the reference allele is empty.
    let inputs = lines(&[
        "21\t2000\tfurther_one\tG\tGAGAA\t.\t.\t.",
        "21\t2000\treplaced\tGC\tAGAA\t.\t.\t.",
    ]);
    let replaced = Written {
        row: ("replaced", "21:2000-2001", "AGAA", "ENST1", "625-626"),
        ..written_further_on().into_iter().next().expect("a row")
    };
    let listed = default_entries(&default_text(&[
        first.text(&first.reference),
        replaced.text(&replaced.reference),
    ]));
    let lost = default_entries(&default_text(&[
        first.text(&first.candidate),
        replaced.text(&replaced.candidate),
    ]));
    for (at, handed) in [(0, true), (1, false)] {
        let found = accepts_with_and_without_the_record(
            UTR_BEFORE_STOP,
            &listed,
            &inputs,
            at,
            &lost[at],
            "-",
            "2",
        );
        assert_eq!(found, [handed, handed, true, true], "{at}");
    }
}

/// Two insertions of six bases written inside a stop codon, whose mutated sequence can be
/// written at the stop codon's first base: the reference reports the lost stop, with the
/// ORF's next stop codon two bases before the base after the insertion (the ORF's own stop
/// codon, moved on) or four before it (a stop codon the inserted bases put in frame), and
/// the candidate nothing.
fn written_inside_stop_codon() -> [Written; 2] {
    [
        Written {
            row: ("whole_one", "21:5000-5001", "TAAGAA", "ENST1", "69-70"),
            strand: "-1",
            notation: "c.-107_-106insTTCTTA",
            offset: None,
            reference: values(
                UTR_STOP_LOSS,
                "AltStop=True:AltStopDistanceToCDS=108:Evidence=False:FrameWithCDS=inFrame:KozakContext=TTGATGT:KozakStrength=Weak",
            ),
            candidate: nothing(),
        },
        Written {
            row: ("whole_two", "21:6000-6001", "ACTAAT", "ENST1", "362-363"),
            strand: "1",
            notation: "c.-100_-99insACTAAT",
            offset: None,
            reference: values(
                UTR_STOP_LOSS,
                "AltStop=True:AltStopDistanceToCDS=103:Evidence=False:FrameWithCDS=outOfFrame:KozakContext=TACATGG:KozakStrength=Moderate",
            ),
            candidate: nothing(),
        },
    ]
}

/// The reference reports as lost a stop codon that stands whole in frame after an
/// insertion whose length is a multiple of three, and the candidate no lost stop: the row
/// is the class's, with both fields. Over a grid of the next stop codon's distance, the row
/// of a six-base allele is the class's where that distance is one, two, four or five more
/// than the base after the insertion, and the row of a three-base allele where it is one
/// or two more; an allele of four or five bases is refused at every distance.
#[test]
fn the_lost_stop_written_before_stop_codon_takes_a_stop_codon_that_stands_whole_in_frame() {
    let rows = written_inside_stop_codon();
    let text = |candidate: bool| {
        let rows: Vec<String> = rows
            .iter()
            .map(|row| {
                row.text(if candidate {
                    &row.candidate
                } else {
                    &row.reference
                })
            })
            .collect();
        default_text(&rows)
    };
    let outcome = before_stop(&text(false), &text(true));
    assert!(outcome.is("green", &[(UTR_BEFORE_STOP, 2)]));
    for name in ["5UTR_consequence", "5UTR_annotation"] {
        assert!(outcome.masks(name, &[(UTR_BEFORE_STOP, 2)]), "{name}");
    }
    // The base after the insertion stands 106 bases before the coding start.
    let alleles: [(&str, &[i64]); 4] = [
        ("TAAGAA", &[107, 108, 110, 111]),
        ("GAA", &[107, 108]),
        ("TAAG", &[]),
        ("TAAGA", &[]),
    ];
    for (allele, fitting) in alleles {
        let mut taken = Vec::new();
        for distance in 103..=114 {
            let row = Written {
                row: ("whole_one", "21:5000-5001", allele, "ENST1", "69-70"),
                reference: values(
                    UTR_STOP_LOSS,
                    &lost_stop_group_of(Some(distance), "TTGATGT", "Weak"),
                ),
                ..written_inside_stop_codon()
                    .into_iter()
                    .next()
                    .expect("a row")
            };
            let (reference, candidate) = row.texts();
            let outcome = before_stop(&reference, &candidate);
            if outcome.green() {
                assert!(
                    outcome.is("green", &[(UTR_BEFORE_STOP, 1)]),
                    "{allele} {distance}"
                );
                taken.push(distance);
            } else {
                assert!(outcome.is("red", &[]), "{allele} {distance}");
            }
        }
        assert_eq!(taken, fitting, "{allele}");
    }
    // Beside an `HGVS_OFFSET` the notation names a base further on, and the distance is
    // counted from the base after the insertion as the record writes it.
    let shifted = Written {
        notation: "c.-106_-105insTTCTTA",
        offset: Some(-1),
        ..written_inside_stop_codon()
            .into_iter()
            .next()
            .expect("a row")
    };
    let (reference, candidate) = shifted.texts();
    assert!(before_stop(&reference, &candidate).is("green", &[(UTR_BEFORE_STOP, 1)]));
}

/// A lost stop the reference alone reports stays a difference where its ORF meets no
/// further stop codon, with or without a distance beside `AltStop=False`, where one of two
/// ORFs meets its next stop codon elsewhere, on a row without an HGVSc, beside a candidate
/// that reports another lost stop or an annotation beside no term or carries no counts, and
/// with the two sides exchanged. A candidate that reports another of the plugin's terms is
/// the class's, and so are two ORFs that both meet their own stop codons.
#[test]
fn the_lost_stop_written_before_stop_codon_refuses_a_lost_stop_no_whole_stop_codon_explains() {
    let [first, _] = written_inside_stop_codon();
    let near = lost_stop_group_of(Some(108), "TTGATGT", "Weak");
    let also_near = lost_stop_group_of(Some(107), "AAAATGG", "Strong");
    let far = lost_stop_group_of(Some(109), "AAAATGG", "Strong");
    let stop_lost = |groups: &[&String]| {
        let groups: Vec<&str> = groups.iter().map(|group| group.as_str()).collect();
        values(UTR_STOP_LOSS, &groups.join(","))
    };
    let cases: [(&str, Values, Values, bool); 9] = [
        ("one ORF", stop_lost(&[&near]), nothing(), true),
        (
            "two ORFs at their own stop codons",
            stop_lost(&[&near, &also_near]),
            nothing(),
            true,
        ),
        (
            "the candidate's frameshift",
            stop_lost(&[&near]),
            frame_shifted(),
            true,
        ),
        (
            "no further stop codon",
            stop_lost(&[&lost_stop_group_of(None, "TTGATGT", "Weak")]),
            nothing(),
            false,
        ),
        (
            "no further stop codon beside a distance",
            stop_lost(&[&near.replace("AltStop=True", "AltStop=False")]),
            nothing(),
            false,
        ),
        (
            "one of two ORFs elsewhere",
            stop_lost(&[&near, &far]),
            nothing(),
            false,
        ),
        (
            "another lost stop",
            stop_lost(&[&near]),
            stop_lost(&[&also_near]),
            false,
        ),
        (
            "an annotation beside no term",
            stop_lost(&[&near]),
            values("-", "KozakStrength=Weak:type=uORF"),
            false,
        ),
        ("the sides exchanged", nothing(), stop_lost(&[&near]), false),
    ];
    for (what, reference, candidate, taken) in &cases {
        let outcome = before_stop(
            &one_row(first.text(reference)),
            &one_row(first.text(candidate)),
        );
        assert!(outcome.taken_is(*taken, &[(UTR_BEFORE_STOP, 1)]), "{what}");
    }
    let reference = stranded_row(first.row, &first.reference, Some("-1"));
    let candidate = stranded_row(first.row, &nothing(), Some("-1"));
    assert!(before_stop(&one_row(reference), &one_row(candidate)).is("red", &[]));
    let (reference, _) = first.texts();
    let uncounted = one_row(
        utr_row(first.row, &nothing()).replacen("STRAND=1;", "STRAND=-1;", 1)
            + ";HGVSc=ENST1.1:c.-107_-106insTTCTTA",
    );
    assert!(before_stop(&reference, &uncounted).is("red", &[]));
}

/// The reference's frameshift of two listed ORFs and the candidate's lost stop of one.
fn before_stop_texts() -> (String, String) {
    (
        one_row(listed_row(
            &[orf_through(), orf_at_stop()],
            UTR_BEFORE_STOP_ROW,
        )),
        one_row(lost_stop_row(&[lost_stop_group(Some(47))])),
    )
}

/// A row two classes hold is counted under the first in the order a row is put to the
/// classes. Where Ensembl VEP's outputs under other hash seeds write both the frameshift and
/// the lost stop on the row, each with its annotation, the row is the precedence class's,
/// whichever class is named first; with runs that write the frameshift alone, without runs,
/// or named alone, it is this class's. A row on which the reference reports the start codon
/// directly after the insertion lost is the lost-start class's and none of this class's.
#[test]
fn the_lost_stop_written_before_stop_codon_follows_the_plugins_other_classes() {
    let (reference, candidate) = before_stop_texts();
    let both_terms = [reference.as_str(), candidate.as_str()];
    let one_term = [reference.as_str()];
    let cases: [(Seeded, Classes, Applied); 5] = [
        (
            Some(&both_terms),
            &[UTR_PRECEDENCE, UTR_BEFORE_STOP],
            &[(UTR_PRECEDENCE, 1)],
        ),
        (
            Some(&both_terms),
            &[UTR_BEFORE_STOP, UTR_PRECEDENCE],
            &[(UTR_PRECEDENCE, 1)],
        ),
        (
            Some(&both_terms),
            &[UTR_BEFORE_STOP],
            &[(UTR_BEFORE_STOP, 1)],
        ),
        (
            Some(&one_term),
            &[UTR_PRECEDENCE, UTR_BEFORE_STOP],
            &[(UTR_BEFORE_STOP, 1)],
        ),
        (
            None,
            &[UTR_PRECEDENCE, UTR_BEFORE_STOP],
            &[(UTR_BEFORE_STOP, 1)],
        ),
    ];
    for (sweep, classes, applied) in cases {
        let outcome = compared_swept(&reference, &candidate, sweep, classes);
        assert!(
            outcome.is("green", applied),
            "{} {classes:?}",
            sweep.is_none()
        );
    }
    let lost_start = one_row(counted_row(UTR_BEFORE_STOP_ROW, &lost_at(146)));
    let outcome = compared(&lost_start, &candidate, &[UTR_LOST_AFTER, UTR_BEFORE_STOP]);
    assert!(outcome.is("green", &[(UTR_LOST_AFTER, 1)]));
    assert!(before_stop(&lost_start, &candidate).is("red", &[]));
}

/// The command line names the lost-stop class like any other.
#[test]
fn the_command_line_names_the_lost_stop_written_before_stop_codon_class() {
    let (reference, candidate) = before_stop_texts();
    let directory = planted_files(
        "utrannotator-before-stop",
        &[("reference", &reference), ("candidate", &candidate)],
    );
    assert_eq!(exit_code(&directory, &[]), 1);
    assert_eq!(exit_code(&directory, &["--masks", UTR_BEFORE_STOP]), 0);
}

/// The tab format carries the two fields as columns, the VCF format as CSQ fields and the
/// JSON format as `5utr_consequence` beside an object of numbered groups. In each the
/// row is the class's where the candidate's group names the listed uORF's next stop codon,
/// and stays a difference where it names another.
#[test]
fn the_lost_stop_written_before_stop_codon_reads_the_tab_vcf_and_json_formats() {
    let under_class = |reference: &[Entry], candidate: &[Entry]| {
        compared_entries(reference, candidate, &[UTR_BEFORE_STOP], None, None)
    };
    let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tSTRAND\t5UTR_consequence\t5UTR_annotation\tExisting_uORFs\n";
    let tab = |term: &str, annotation: &str| {
        format!(
            "{header}grid_ins\t21:1000-1001\tAT\tENSG1\tENST1\tTranscript\t5_prime_UTR_variant\t146-147\t1\t{term}\t{annotation}\t1\n"
        )
    };
    for (distance, taken) in [(47, true), (48, false)] {
        let outcome = under_class(
            &tab_entries(&tab(UTR_FRAMESHIFT, &orf_at_stop())),
            &tab_entries(&tab(UTR_STOP_LOSS, &lost_stop_group(Some(distance)))),
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_BEFORE_STOP, 1)]),
            "{distance}"
        );
    }

    let inputs = lines(&["21\t1000\tgrid_ins\tG\tGAT\t.\t.\t."]);
    let csq = csq_header(
        "Allele|Consequence|Feature|Feature_type|cDNA_position|5UTR_consequence|5UTR_annotation|Existing_uORFs",
    );
    let vcf_text = |term: &str, annotation: &str| {
        format!(
            "{csq}{}CSQ=AT|5_prime_UTR_variant|ENST1|Transcript|146-147|{term}|{annotation}|1\n",
            record_columns(&inputs[0])
        )
    };
    for (distance, taken) in [(47, true), (48, false)] {
        let outcome = under_class(
            &vcf_entries(
                &vcf_text(
                    UTR_FRAMESHIFT,
                    &format!("{}&{}", orf_through(), orf_at_stop()),
                ),
                &inputs,
            ),
            &vcf_entries(
                &vcf_text(UTR_STOP_LOSS, &lost_stop_group(Some(distance))),
                &inputs,
            ),
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_BEFORE_STOP, 1)]),
            "{distance}"
        );
    }

    let object = |term: &str, groups: &[String]| {
        json_line(
            &inputs[0],
            "grid_ins",
            "-/AT",
            &format!(
                r#""variant_allele": "AT", "transcript_id": "ENST1", "consequence_terms": ["5_prime_UTR_variant"], "cdna_start": 147, "cdna_end": 146, "strand": 1, "existing_uorfs": 1, "5utr_consequence": "{term}", "5utr_annotation": {}"#,
                json_groups(groups)
            ),
        )
    };
    for (distance, taken) in [(47, true), (48, false)] {
        let outcome = under_class(
            &json_entries(
                &object(UTR_FRAMESHIFT, &[orf_through(), orf_at_stop()]),
                &inputs,
            ),
            &json_entries(
                &object(UTR_STOP_LOSS, &[lost_stop_group(Some(distance))]),
                &inputs,
            ),
        );
        assert!(
            outcome.taken_is(taken, &[(UTR_BEFORE_STOP, 1)]),
            "{distance}"
        );
    }
}
