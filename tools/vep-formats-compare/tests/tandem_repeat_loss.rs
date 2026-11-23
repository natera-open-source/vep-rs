// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! `cnv_tr_symbolic_loss_named_a_gain`: the Consequence of a transcript row whose Allele is
//! `tandem_repeat`, on which the reference writes `feature_elongation` and the candidate
//! `feature_truncation` in its place or no term in its place, and the IMPACT that follows
//! where the candidate writes none. The rows are made up: three records of one Location,
//! so that the rows of one key stand apart by their order alone, as the rows of several
//! `<CNV:TR>` records at one position do.

mod compare_common;
mod term_pairs_common;

use std::collections::BTreeMap;

use compare_common::*;
use term_pairs_common::*;
use vep_formats_compare::classes;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::term_pairs::symbolic_loss_named_a_gain_pair;
use vep_formats_compare::rules::term_pairs::TERM_PAIR_CLASSES;
use vep_formats_compare::rules::{field_class_accepts, RuleContext, Side};

const CLASS: &str = classes::CNV_TR_SYMBOLIC_LOSS_NAMED_A_GAIN;
const LOCATION: &str = "7:1001-1009";
const ALLELE: &str = "tandem_repeat";
const FEATURE: &str = "ENST90000000001";
const GAIN: &str = "feature_elongation,coding_sequence_variant";
const LOSS: &str = "feature_truncation,coding_sequence_variant";
const EXTRA: &str = "ALLELE_NUM=1;IMPACT=HIGH;STRAND=1;VARIANT_CLASS=tandem_repeat";
/// The candidate's terms on the row of an allele as long as its run, and its Extra: the
/// IMPACT is `coding_sequence_variant`'s (Utils/Constants.pm:1061-1078).
const SAME: &str = "coding_sequence_variant";
const SAME_EXTRA: &str = "ALLELE_NUM=1;IMPACT=MODIFIER;STRAND=1;VARIANT_CLASS=tandem_repeat";

/// What was applied, as a count map.
fn applied(documented: &Documented) -> BTreeMap<String, i64> {
    counts(&documented.applied)
}

fn one(klass: &str) -> BTreeMap<String, i64> {
    count_map(&[(klass, 1)])
}

/// The classes a column's differences were set aside under, each with its count.
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

fn differing(report: &Value, column: &str) -> i64 {
    number(report, &["columns", column, "differing"])
}

/// One default-format row of the record `name` on the made-up transcript.
fn row(name: &str, allele: &str, terms: &str) -> String {
    default_row_of(
        name,
        LOCATION,
        allele,
        FEATURE,
        terms,
        EXTRA,
        "Transcript",
        "ENSG9",
    )
}

/// A default-format output of three records at one Location: a gain, a loss and an allele as
/// long as the run, the loss's row with the terms `loss`.
fn three_records(loss: &str) -> String {
    default_text(&[
        row("tr_gain", ALLELE, GAIN),
        row("tr_loss", ALLELE, loss),
        row("tr_same", ALLELE, GAIN),
    ])
}

/// The output `three_records(GAIN)` with the row of the allele as long as the run, the
/// third, carrying the terms `same` and the Extra `extra`.
fn as_long_as_the_run(same: &str, extra: &str) -> String {
    default_text(&[
        row("tr_gain", ALLELE, GAIN),
        row("tr_loss", ALLELE, GAIN),
        default_row_of(
            "tr_same",
            LOCATION,
            ALLELE,
            FEATURE,
            same,
            extra,
            "Transcript",
            "ENSG9",
        ),
    ])
}

/// The verdict of two default-format texts under the class.
fn verdict_under_the_class(reference: &str, candidate: &str) -> String {
    verdict(&compare_default(reference, candidate, &[CLASS]).0).to_string()
}

/// The reference entry and the candidate entry of the made-up row at `place`, counted from 0.
fn pair_at(reference: &str, candidate: &str, place: usize) -> (Entry, Entry) {
    let of = |text: &str| {
        let entries = vep_formats_compare::parse::text::parse_default(text).expect("it parses");
        entries.into_iter().nth(place).expect("the row")
    };
    (of(reference), of(candidate))
}

/// The reference entry and the candidate entry of one made-up row, the loss's.
fn pair_of(reference: &str, candidate: &str) -> (Entry, Entry) {
    pair_at(reference, candidate, 1)
}

/// The reference entry and the candidate entry of the row of the allele as long as the run.
fn pair_of_the_same(reference: &str, candidate: &str) -> (Entry, Entry) {
    pair_at(reference, candidate, 2)
}

/// The Consequence of the one row on which the candidate writes the truncation is set aside,
/// counted once, and the two other rows of the key are compared; without the class named the
/// pair is a difference.
#[test]
fn the_class_takes_the_consequence_of_a_loss_the_reference_names_a_gain() {
    let (reference, candidate) = (three_records(GAIN), three_records(LOSS));
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    assert_eq!(differing(&without, "Consequence"), 1);
    let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(CLASS));
    assert_eq!(masked(&report, "Consequence"), one(CLASS));
    assert_eq!(number(&report, &["counts", "compared"]), 3);
    assert!(
        texts(at(&report, &["mask_examples", CLASS])).contains(
            &format!(
                "{LOCATION}|{ALLELE}|{FEATURE}|Transcript#2: Consequence expected '{GAIN}', got '{LOSS}'"
            )
            .as_str()
        )
    );
    let (e, a) = pair_of(&reference, &candidate);
    assert!(holds(symbolic_loss_named_a_gain_pair, &e, &a));
    // The tab format reads the same rows.
    let tab = |text: &str| parse_tab(text).expect("the output parses");
    let mut under_the_class = documented_with(&[CLASS]);
    let report = compare(&tab(&reference), &tab(&candidate), &mut under_the_class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&under_the_class), one(CLASS));
    let report = compare(&tab(&reference), &tab(&candidate), &mut no_documented());
    assert_eq!(verdict(&report), "red");
}

/// The two terms the other way round, a row of another Allele, a regulatory row, a pair of
/// lists that differ in another place as well and a candidate that holds both terms stay
/// differences, and so does every other cell of a row the class takes.
#[test]
fn every_other_pair_and_every_other_cell_stay_differences() {
    let (reference, candidate) = (three_records(GAIN), three_records(LOSS));
    // The reference names the truncation and the candidate the elongation.
    assert_eq!(verdict_under_the_class(&candidate, &reference), "red");
    let (e, a) = pair_of(&candidate, &reference);
    assert!(!holds(symbolic_loss_named_a_gain_pair, &e, &a));
    // A row of another Allele: another structural class's, and a deletion's.
    for allele in ["duplication", "deletion", "-"] {
        let other =
            |terms: &str| default_text(&[row("tr_gain", ALLELE, GAIN), row("sv", allele, terms)]);
        let (report, documented) = compare_default(&other(GAIN), &other(LOSS), &[CLASS]);
        assert_eq!(verdict(&report), "red", "{allele}");
        assert_eq!(applied(&documented), BTreeMap::new(), "{allele}");
    }
    // A row of a feature that is no transcript.
    let regulatory = |terms: &str| {
        default_text(&[default_row_of(
            "tr_loss",
            LOCATION,
            ALLELE,
            "ENSR90000000001",
            terms,
            EXTRA,
            "RegulatoryFeature",
            "-",
        )])
    };
    assert_eq!(
        verdict_under_the_class(&regulatory(GAIN), &regulatory(LOSS)),
        "red"
    );
    // Two lists that differ in another place, in their length or in the order of their terms,
    // and a candidate that holds both terms beside a reference that names the elongation twice.
    for (theirs, ours) in [
        (
            "feature_elongation,feature_elongation,coding_sequence_variant",
            "feature_truncation,feature_elongation,coding_sequence_variant",
        ),
        (
            "feature_elongation,coding_sequence_variant,feature_elongation",
            "feature_elongation,coding_sequence_variant,feature_truncation",
        ),
        (GAIN, "feature_truncation,3_prime_UTR_variant"),
        (GAIN, "feature_truncation"),
        (
            GAIN,
            "feature_truncation,coding_sequence_variant,intron_variant",
        ),
        (GAIN, "coding_sequence_variant,feature_truncation"),
        (GAIN, "feature_elongation,feature_truncation"),
        ("feature_elongation", "transcript_ablation"),
        ("transcript_amplification", "feature_truncation"),
        (
            "feature_elongation,feature_elongation",
            "feature_truncation,feature_truncation",
        ),
    ] {
        let (report, documented) =
            compare_default(&three_records(theirs), &three_records(ours), &[CLASS]);
        assert_eq!(verdict(&report), "red", "{theirs} against {ours}");
        assert_eq!(
            applied(&documented),
            BTreeMap::new(),
            "{theirs} against {ours}"
        );
        let (e, a) = pair_of(&three_records(theirs), &three_records(ours));
        assert!(!holds(symbolic_loss_named_a_gain_pair, &e, &a), "{ours}");
    }
    // The class takes the Consequence alone of a row that names the truncation: an IMPACT,
    // a strand or an allele number that differs on the row is a difference beside it.
    for (written, changed, column) in [
        ("IMPACT=HIGH", "IMPACT=MODERATE", "IMPACT"),
        ("STRAND=1", "STRAND=-1", "STRAND"),
        ("ALLELE_NUM=1", "ALLELE_NUM=2", "ALLELE_NUM"),
    ] {
        let lines: Vec<String> = candidate
            .lines()
            .map(|line| {
                if line.starts_with("tr_loss") {
                    line.replace(written, changed)
                } else {
                    line.to_string()
                }
            })
            .collect();
        let changed_row = lines.join("\n") + "\n";
        let (report, documented) = compare_default(&reference, &changed_row, &[CLASS]);
        assert_eq!(verdict(&report), "red", "{column}");
        assert_eq!(applied(&documented), one(CLASS), "{column}");
        assert_eq!(differing(&report, column), 1, "{column}");
    }
    // A row on which the candidate writes the reference's terms is no row of the class.
    let (report, documented) = compare_default(&reference, &reference, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), BTreeMap::new());
}

/// The row of an allele as long as its run: the candidate writes the reference's terms
/// without `feature_elongation` and the IMPACT its own terms give. The Consequence and the
/// IMPACT of that row are set aside, the row counted once, and the two other rows of the
/// key are compared; without the class named both cells are differences.
#[test]
fn the_class_takes_the_consequence_and_the_impact_of_an_allele_as_long_as_its_run() {
    let (reference, candidate) = (three_records(GAIN), as_long_as_the_run(SAME, SAME_EXTRA));
    let (without, _) = compare_default(&reference, &candidate, &[]);
    assert_eq!(verdict(&without), "red");
    assert_eq!(differing(&without, "Consequence"), 1);
    assert_eq!(differing(&without, "IMPACT"), 1);
    let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(CLASS));
    assert_eq!(masked(&report, "Consequence"), one(CLASS));
    assert_eq!(masked(&report, "IMPACT"), one(CLASS));
    assert_eq!(number(&report, &["counts", "compared"]), 3);
    assert!(
        texts(at(&report, &["mask_examples", CLASS])).contains(
            &format!(
                "{LOCATION}|{ALLELE}|{FEATURE}|Transcript#3: Consequence expected '{GAIN}', got '{SAME}'"
            )
            .as_str()
        )
    );
    let (e, a) = pair_of_the_same(&reference, &candidate);
    assert!(holds(symbolic_loss_named_a_gain_pair, &e, &a));
    // The tab format reads the same rows.
    let tab = |text: &str| parse_tab(text).expect("the output parses");
    let mut under_the_class = documented_with(&[CLASS]);
    let report = compare(&tab(&reference), &tab(&candidate), &mut under_the_class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&under_the_class), one(CLASS));
    let report = compare(&tab(&reference), &tab(&candidate), &mut no_documented());
    assert_eq!(verdict(&report), "red");
    // The terms left keep their order, and the IMPACT is that of the first of them by rank:
    // the reference's own where a term above `feature_elongation` stands.
    for (theirs, ours, extra) in [
        (
            "feature_elongation,coding_sequence_variant,intron_variant",
            "coding_sequence_variant,intron_variant",
            SAME_EXTRA,
        ),
        (
            "feature_elongation,3_prime_UTR_variant",
            "3_prime_UTR_variant",
            SAME_EXTRA,
        ),
        (
            "feature_elongation,non_coding_transcript_exon_variant",
            "non_coding_transcript_exon_variant",
            SAME_EXTRA,
        ),
        (
            "start_lost,feature_elongation,coding_sequence_variant,5_prime_UTR_variant",
            "start_lost,coding_sequence_variant,5_prime_UTR_variant",
            EXTRA,
        ),
    ] {
        let reference = as_long_as_the_run(theirs, EXTRA);
        let candidate = as_long_as_the_run(ours, extra);
        let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
        assert_eq!(verdict(&report), "green", "{theirs} against {ours}");
        assert_eq!(applied(&documented), one(CLASS), "{theirs} against {ours}");
        let (report, _) = compare_default(&reference, &candidate, &[]);
        assert_eq!(verdict(&report), "red", "{theirs} against {ours}");
    }
    // A row that carries no IMPACT has none to hold: the terms alone are the shape.
    let bare = "ALLELE_NUM=1;STRAND=1";
    let (reference, candidate) = (
        as_long_as_the_run(GAIN, bare),
        as_long_as_the_run(SAME, bare),
    );
    let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&documented), one(CLASS));
}

/// Beside the row of an allele as long as its run: the candidate holding
/// `feature_elongation` where the reference holds none, a row of another Allele, a
/// regulatory row, a pair of lists that differ in another term, in two terms or in the order
/// of the terms left, and a candidate that holds `feature_truncation` stay differences.
#[test]
fn every_other_pair_beside_an_allele_as_long_as_its_run_stays_a_difference() {
    let (reference, candidate) = (three_records(GAIN), as_long_as_the_run(SAME, SAME_EXTRA));
    // The candidate names the elongation and the reference does not.
    let (report, documented) = compare_default(&candidate, &reference, &[CLASS]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), BTreeMap::new());
    let (e, a) = pair_of_the_same(&candidate, &reference);
    assert!(!holds(symbolic_loss_named_a_gain_pair, &e, &a));
    // A row of another Allele, and a row of a feature that is no transcript.
    let other = |allele: &str, feature: &str, feature_type: &str, terms: &str, extra: &str| {
        default_text(&[
            row("tr_gain", ALLELE, GAIN),
            default_row_of(
                "sv",
                LOCATION,
                allele,
                feature,
                terms,
                extra,
                feature_type,
                "ENSG9",
            ),
        ])
    };
    for (allele, feature, feature_type) in [
        ("duplication", FEATURE, "Transcript"),
        ("deletion", FEATURE, "Transcript"),
        ("-", FEATURE, "Transcript"),
        (ALLELE, "ENSR90000000001", "RegulatoryFeature"),
        (ALLELE, "ENSM90000000001", "MotifFeature"),
    ] {
        let (report, documented) = compare_default(
            &other(allele, feature, feature_type, GAIN, EXTRA),
            &other(allele, feature, feature_type, SAME, SAME_EXTRA),
            &[CLASS],
        );
        assert_eq!(verdict(&report), "red", "{allele} on {feature_type}");
        assert_eq!(
            applied(&documented),
            BTreeMap::new(),
            "{allele} on {feature_type}"
        );
    }
    // Each candidate row carries the IMPACT its own terms give, so that its terms decide.
    for (theirs, ours, extra) in [
        (GAIN, "3_prime_UTR_variant", SAME_EXTRA),
        (GAIN, "coding_sequence_variant,intron_variant", SAME_EXTRA),
        (
            "feature_elongation,coding_sequence_variant,intron_variant",
            "coding_sequence_variant",
            SAME_EXTRA,
        ),
        (
            "feature_elongation,coding_sequence_variant,intron_variant",
            "intron_variant,coding_sequence_variant",
            SAME_EXTRA,
        ),
        (
            "feature_elongation,feature_truncation,coding_sequence_variant",
            "feature_truncation,coding_sequence_variant",
            EXTRA,
        ),
        (
            "feature_elongation,feature_elongation,coding_sequence_variant",
            "coding_sequence_variant",
            SAME_EXTRA,
        ),
        (
            "feature_elongation,feature_elongation,coding_sequence_variant",
            "feature_elongation,coding_sequence_variant",
            EXTRA,
        ),
        ("feature_elongation", "-", SAME_EXTRA),
        (
            "transcript_amplification,coding_sequence_variant",
            "coding_sequence_variant",
            SAME_EXTRA,
        ),
        (
            "feature_truncation,coding_sequence_variant",
            "coding_sequence_variant",
            SAME_EXTRA,
        ),
    ] {
        let reference = as_long_as_the_run(theirs, EXTRA);
        let candidate = as_long_as_the_run(ours, extra);
        let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
        assert_eq!(verdict(&report), "red", "{theirs} against {ours}");
        assert_eq!(
            applied(&documented),
            BTreeMap::new(),
            "{theirs} against {ours}"
        );
        let (e, a) = pair_of_the_same(&reference, &candidate);
        assert!(!holds(symbolic_loss_named_a_gain_pair, &e, &a), "{ours}");
    }
}

/// The IMPACT of the row of an allele as long as its run is held to the terms: a candidate
/// IMPACT that is not its own terms' leaves the row no row of the class, whether or not it
/// equals the reference's; a reference IMPACT that is not its terms' stays a difference
/// beside the Consequence the class takes; an IMPACT that differs on a row whose terms agree
/// stays a difference; and every other cell of a row the class takes stays compared.
#[test]
fn the_impact_of_an_allele_as_long_as_its_run_is_held_to_each_sides_terms() {
    let reference = three_records(GAIN);
    for impact in ["HIGH", "MODERATE", "LOW"] {
        let extra = SAME_EXTRA.replace("IMPACT=MODIFIER", &format!("IMPACT={impact}"));
        let candidate = as_long_as_the_run(SAME, &extra);
        let (report, documented) = compare_default(&reference, &candidate, &[CLASS]);
        assert_eq!(verdict(&report), "red", "{impact}");
        assert_eq!(applied(&documented), BTreeMap::new(), "{impact}");
        assert_eq!(differing(&report, "Consequence"), 1, "{impact}");
        let (e, a) = pair_of_the_same(&reference, &candidate);
        assert!(!holds(symbolic_loss_named_a_gain_pair, &e, &a), "{impact}");
    }
    // The reference's IMPACT is not the one its terms give.
    let moderate = as_long_as_the_run(GAIN, &EXTRA.replace("IMPACT=HIGH", "IMPACT=MODERATE"));
    let candidate = as_long_as_the_run(SAME, SAME_EXTRA);
    let (report, documented) = compare_default(&moderate, &candidate, &[CLASS]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), one(CLASS));
    assert_eq!(masked(&report, "Consequence"), one(CLASS));
    assert_eq!(differing(&report, "IMPACT"), 1);
    // The terms agree and the IMPACT differs.
    let lowered = as_long_as_the_run(GAIN, SAME_EXTRA);
    let (report, documented) = compare_default(&reference, &lowered, &[CLASS]);
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&documented), BTreeMap::new());
    assert_eq!(differing(&report, "IMPACT"), 1);
    // A strand, an allele number or a class that differs on the row is a difference beside
    // the two cells the class takes.
    for (written, changed, column) in [
        ("STRAND=1", "STRAND=-1", "STRAND"),
        ("ALLELE_NUM=1", "ALLELE_NUM=2", "ALLELE_NUM"),
        (
            "VARIANT_CLASS=tandem_repeat",
            "VARIANT_CLASS=duplication",
            "VARIANT_CLASS",
        ),
    ] {
        let changed_row = as_long_as_the_run(SAME, &SAME_EXTRA.replace(written, changed));
        let (report, documented) = compare_default(&reference, &changed_row, &[CLASS]);
        assert_eq!(verdict(&report), "red", "{column}");
        assert_eq!(applied(&documented), one(CLASS), "{column}");
        assert_eq!(differing(&report, column), 1, "{column}");
        assert_eq!(differing(&report, "IMPACT"), 0, "{column}");
    }
}

const INPUT: &str = "7\t1000\ttr_multi\tG\t<CNV:TR>,<CNV:TR>\t.\t.\t.";

/// A VCF output of the record `INPUT`: two alleles on the made-up transcript, the second
/// with the terms `second` and the IMPACT `impact`.
fn vcf_with(second: &str, impact: &str) -> Vec<Entry> {
    let entry = |terms: &str, impact: &str, number: &str| {
        format!(
            "{ALLELE}|{}|{impact}|{FEATURE}|Transcript|{number}",
            terms.replace(',', "&")
        )
    };
    let text = format!(
        "{}{}CSQ={},{}\n",
        csq_header("Allele|Consequence|IMPACT|Feature|Feature_type|ALLELE_NUM"),
        open_info(INPUT),
        entry(GAIN, "HIGH", "1"),
        entry(second, impact, "2"),
    );
    parse_vcf(&text, &[INPUT.to_string()]).expect("the output parses")
}

/// `vcf_with` for a second allele whose IMPACT is HIGH.
fn vcf(second: &str) -> Vec<Entry> {
    vcf_with(second, "HIGH")
}

/// One member of a JSON record's `transcript_consequences` on the made-up transcript.
fn json_entry(terms: &str, impact: &str, number: u8) -> String {
    let terms: Vec<&str> = terms.split(',').collect();
    format!(
        "{{\"variant_allele\": \"{ALLELE}\", \"transcript_id\": \"{FEATURE}\", \"impact\": \"{impact}\", \"allele_num\": {number}, \"strand\": 1, \"consequence_terms\": {}}}",
        quoted_list(&terms)
    )
}

/// A JSON output of the record `INPUT`: two alleles on the made-up transcript, the second
/// with the terms `second`, under the record's `most_severe`.
fn json_text(second: &str, most_severe: &str) -> String {
    format!(
        "{{\"input\": {}, \"id\": \"tr_multi\", \"allele_string\": \"<CNV:TR>/<CNV:TR>\", \"most_severe_consequence\": {}, \"transcript_consequences\": [{}, {}]}}\n",
        quoted(INPUT),
        quoted(most_severe),
        json_entry(GAIN, "HIGH", 1),
        json_entry(second, "HIGH", 2),
    )
}

const ONE_ALLELE: &str = "7\t1000\ttr_same\tG\t<CNV:TR>\t.\t.\t.";

/// A JSON output of the record `ONE_ALLELE`: its one allele on the made-up transcript with
/// the terms `terms` and the IMPACT `impact`, under the record's `most_severe`.
fn one_allele_json(terms: &str, impact: &str, most_severe: &str) -> Vec<Entry> {
    let text = format!(
        "{{\"input\": {}, \"id\": \"tr_same\", \"allele_string\": \"<CNV:TR>\", \"most_severe_consequence\": {}, \"transcript_consequences\": [{}]}}\n",
        quoted(ONE_ALLELE),
        quoted(most_severe),
        json_entry(terms, impact, 1),
    );
    parse_json(&text, &[ONE_ALLELE.to_string()]).expect("the output parses")
}

fn json_entries(text: &str) -> Vec<Entry> {
    parse_json(text, &[INPUT.to_string()]).expect("the output parses")
}

/// In the VCF and the JSON format the row is the second of its record's two rows of one
/// Allele, counted once; the JSON record's `most_severe_consequence` is left uncompared on
/// the record and its other members stay compared.
#[test]
fn the_class_counts_a_row_once_in_the_vcf_and_json_formats() {
    let mut under_the_class = documented_with(&[CLASS]);
    let report = compare(&vcf(GAIN), &vcf(LOSS), &mut under_the_class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&under_the_class), one(CLASS));
    let report = compare(&vcf(GAIN), &vcf(LOSS), &mut no_documented());
    assert_eq!(verdict(&report), "red");
    let report = compare(&vcf(LOSS), &vcf(GAIN), &mut documented_with(&[CLASS]));
    assert_eq!(verdict(&report), "red");

    // The record's most severe term moves with the row where the row's is the record's.
    let reference = json_entries(&json_text(GAIN, "feature_elongation"));
    for most_severe in ["feature_elongation", "feature_truncation"] {
        let candidate = json_entries(&json_text(LOSS, most_severe));
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(&reference, &candidate, &mut under_the_class);
        assert_eq!(verdict(&report), "green", "{most_severe}");
        assert_eq!(applied(&under_the_class), one(CLASS), "{most_severe}");
        let report = compare(&reference, &candidate, &mut no_documented());
        assert_eq!(verdict(&report), "red", "{most_severe}");
    }
    // A record without a row of the class keeps its most severe term compared.
    let report = compare(
        &reference,
        &json_entries(&json_text(GAIN, "feature_truncation")),
        &mut documented_with(&[CLASS]),
    );
    assert_eq!(verdict(&report), "red");
    // The record's other members and the row's other members stay compared.
    for (written, changed) in [
        ("\"id\": \"tr_multi\"", "\"id\": \"tr_other\""),
        ("\"allele_num\": 2", "\"allele_num\": 3"),
        (
            "\"strand\": 1, \"consequence_terms\": [\"feature_truncation\"",
            "\"strand\": -1, \"consequence_terms\": [\"feature_truncation\"",
        ),
    ] {
        let text = json_text(LOSS, "feature_truncation");
        assert!(text.contains(written), "{written}");
        let candidate = json_entries(&text.replace(written, changed));
        let report = compare(&reference, &candidate, &mut documented_with(&[CLASS]));
        assert_eq!(verdict(&report), "red", "{changed}");
    }
}

/// In the VCF and the JSON format the row of an allele as long as its run is counted once
/// with its Consequence and its IMPACT; the JSON record's `most_severe_consequence`, the
/// candidate's own most severe term, is left uncompared on the record.
#[test]
fn the_class_counts_an_allele_as_long_as_its_run_once_in_the_vcf_and_json_formats() {
    let mut under_the_class = documented_with(&[CLASS]);
    let (reference, candidate) = (vcf(GAIN), vcf_with(SAME, "MODIFIER"));
    let report = compare(&reference, &candidate, &mut under_the_class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&under_the_class), one(CLASS));
    assert_eq!(masked(&report, "Consequence"), one(CLASS));
    assert_eq!(masked(&report, "IMPACT"), one(CLASS));
    let report = compare(&reference, &candidate, &mut no_documented());
    assert_eq!(verdict(&report), "red");
    let report = compare(&candidate, &reference, &mut documented_with(&[CLASS]));
    assert_eq!(verdict(&report), "red");
    // An IMPACT that is not the candidate's terms' leaves the row no row of the class.
    for impact in ["HIGH", "LOW"] {
        let mut under_the_class = documented_with(&[CLASS]);
        let report = compare(&reference, &vcf_with(SAME, impact), &mut under_the_class);
        assert_eq!(verdict(&report), "red", "{impact}");
        assert_eq!(applied(&under_the_class), BTreeMap::new(), "{impact}");
    }

    // The record's most severe term moves with its one row.
    let reference = one_allele_json(GAIN, "HIGH", "feature_elongation");
    let candidate = one_allele_json(SAME, "MODIFIER", SAME);
    let mut under_the_class = documented_with(&[CLASS]);
    let report = compare(&reference, &candidate, &mut under_the_class);
    assert_eq!(verdict(&report), "green");
    assert_eq!(applied(&under_the_class), one(CLASS));
    assert_eq!(masked(&report, "consequence_terms"), one(CLASS));
    assert_eq!(masked(&report, "impact"), one(CLASS));
    let report = compare(&reference, &candidate, &mut no_documented());
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "most_severe_consequence"), 1);
    // A record without a row of the class keeps its most severe term compared.
    let report = compare(
        &reference,
        &one_allele_json(GAIN, "HIGH", SAME),
        &mut documented_with(&[CLASS]),
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "most_severe_consequence"), 1);
    // An IMPACT that is not the candidate's terms' leaves the row and the record compared.
    let mut under_the_class = documented_with(&[CLASS]);
    let report = compare(
        &reference,
        &one_allele_json(SAME, "HIGH", SAME),
        &mut under_the_class,
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(applied(&under_the_class), BTreeMap::new());
    assert_eq!(differing(&report, "most_severe_consequence"), 1);
}

/// A manifest documents the key several records share with the reference's own set among
/// the candidate's sets. The row that differs is the manifest's; a row on which the candidate
/// writes the reference's set is compared like any other, so a field that differs on it is
/// a difference.
#[test]
fn a_manifest_item_of_a_shared_key_leaves_the_agreeing_rows_compared() {
    let manifest = json(&format!(
        r#"{{"divergences": [{{"location": "{LOCATION}", "allele": "{ALLELE}",
            "feature": "{FEATURE}", "feature_type": "Transcript", "record_indices": [1],
            "vep_consequence_set": "coding_sequence_variant,feature_elongation",
            "vep_rs_consequence_sets": ["coding_sequence_variant,feature_elongation",
                "coding_sequence_variant,feature_truncation"],
            "expected_divergence": "{CLASS}"}}]}}"#
    ));
    let (reference, candidate) = (three_records(GAIN), three_records(LOSS));
    let report = compare_texts(
        &reference,
        &candidate,
        &mut from_manifest(&manifest, &[CLASS]),
    );
    assert_eq!(verdict(&report), "green");
    assert_eq!(number(&report, &["counts", "compared"]), 2);
    // The first row is the gain's, on which both engines write one set.
    let changed = candidate.replacen("STRAND=1", "STRAND=-1", 1);
    let report = compare_texts(
        &reference,
        &changed,
        &mut from_manifest(&manifest, &[CLASS]),
    );
    assert_eq!(verdict(&report), "red");
    assert_eq!(differing(&report, "STRAND"), 1);
    // A row that writes a set the item does not state is stale.
    let other = three_records("feature_elongation,3_prime_UTR_variant");
    let report = compare_texts(&reference, &other, &mut from_manifest(&manifest, &[CLASS]));
    assert_eq!(verdict(&report), "red");
    assert_eq!(texts(at(&report, &["stale_divergences"])).len(), 1);
    // An item of three sets, the third the reference's without `feature_elongation`: the
    // row of an allele as long as its run is the manifest's as the loss's is, whatever its
    // IMPACT, and under the item of two sets it is stale.
    let three_sets = json(&format!(
        r#"{{"divergences": [{{"location": "{LOCATION}", "allele": "{ALLELE}",
            "feature": "{FEATURE}", "feature_type": "Transcript", "record_indices": [1, 2],
            "vep_consequence_set": "coding_sequence_variant,feature_elongation",
            "vep_rs_consequence_sets": ["coding_sequence_variant",
                "coding_sequence_variant,feature_elongation",
                "coding_sequence_variant,feature_truncation"],
            "expected_divergence": "{CLASS}"}}]}}"#
    ));
    let both = default_text(&[
        row("tr_gain", ALLELE, GAIN),
        row("tr_loss", ALLELE, LOSS),
        default_row_of(
            "tr_same",
            LOCATION,
            ALLELE,
            FEATURE,
            SAME,
            SAME_EXTRA,
            "Transcript",
            "ENSG9",
        ),
    ]);
    let report = compare_texts(&reference, &both, &mut from_manifest(&three_sets, &[CLASS]));
    assert_eq!(verdict(&report), "green");
    assert_eq!(number(&report, &["counts", "compared"]), 1);
    let report = compare_texts(&reference, &both, &mut from_manifest(&manifest, &[CLASS]));
    assert_eq!(verdict(&report), "red");
    assert_eq!(texts(at(&report, &["stale_divergences"])).len(), 1);
}

/// The class is a row class, asked before the class of the splice term past the record's
/// span, and a field class of the Consequence and of the IMPACT, before the two
/// plugin-value classes; it is no class of `TERM_PAIR_FOLLOWING_CLASSES`, whose classes
/// read the plugin fields and the protein notation as well.
#[test]
fn the_class_stands_in_the_tables_of_the_row_classes_and_of_the_field_classes() {
    let tables = classes::tables();
    let rows = tables.term_pair_classes();
    let place = |klass: &str| rows.iter().position(|held| *held == klass);
    let held: Vec<&str> = TERM_PAIR_CLASSES.iter().map(|(name, _)| *name).collect();
    assert_eq!(held, rows);
    let splice_term = place(classes::SPLICE_TERM_PAST_THE_RECORD_SPAN)
        .expect("the class of the splice term past the record's span is a row class");
    assert_eq!(place(CLASS), Some(splice_term - 1));
    assert!(!tables.in_term_pair_following_classes(CLASS));
    assert_eq!(
        tables.field_classes_of(CLASS),
        Some(["Consequence", "consequence_terms", "IMPACT", "impact"].as_slice())
    );
    let fields: Vec<&str> = tables
        .field_classes()
        .iter()
        .map(|(klass, _)| *klass)
        .collect();
    let field_place = |klass: &str| fields.iter().position(|held| *held == klass);
    assert!(field_place(CLASS).is_some());
    assert!(
        field_place(CLASS)
            < field_place(classes::PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN)
    );
    assert_eq!(
        fields[fields.len() - 2..],
        [
            classes::PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN,
            classes::PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR,
        ]
    );
}

/// Asked as a field class, the class answers for the Consequence as its row predicate
/// does, for the IMPACT of a row on which the candidate writes no term in the place of
/// `feature_elongation`, where each side's IMPACT is the one its terms give, and for no
/// other field.
#[test]
fn the_field_class_answers_as_the_row_class() {
    let (reference, candidate) = (three_records(GAIN), three_records(LOSS));
    let (e, a) = pair_of(&reference, &candidate);
    let side = Side::default();
    let context = RuleContext {
        record_rows: None,
        candidate_rows: None,
        declared: &[],
        assembly: None,
        dbnsfp_columns: &[],
        side: &side,
    };
    let accepts = |e: &Entry, a: Option<&Entry>, name: &str, ev: &str, av: &str| {
        field_class_accepts(CLASS, e, name, ev, av, a, &context).expect("the class answers")
    };
    assert!(accepts(&e, Some(&a), "Consequence", GAIN, LOSS));
    assert!(accepts(&e, Some(&a), "consequence_terms", GAIN, LOSS));
    assert!(!accepts(&e, None, "Consequence", GAIN, LOSS));
    assert!(!accepts(&a, Some(&e), "Consequence", LOSS, GAIN));
    assert!(!accepts(&e, Some(&a), "IMPACT", "HIGH", "MODERATE"));
    assert!(!accepts(
        &e,
        Some(&a),
        "most_severe_consequence",
        GAIN,
        LOSS
    ));
    assert!(!accepts(&e, Some(&a), "IMPACT", "HIGH", "MODIFIER"));
    let candidate = as_long_as_the_run(SAME, SAME_EXTRA);
    let (e, a) = pair_of_the_same(&reference, &candidate);
    assert!(accepts(&e, Some(&a), "Consequence", GAIN, SAME));
    assert!(accepts(&e, Some(&a), "consequence_terms", GAIN, SAME));
    assert!(accepts(&e, Some(&a), "IMPACT", "HIGH", "MODIFIER"));
    assert!(accepts(&e, Some(&a), "impact", "HIGH", "MODIFIER"));
    assert!(!accepts(&e, None, "Consequence", GAIN, SAME));
    assert!(!accepts(&e, None, "IMPACT", "HIGH", "MODIFIER"));
    assert!(!accepts(&a, Some(&e), "Consequence", SAME, GAIN));
    assert!(!accepts(&a, Some(&e), "IMPACT", "MODIFIER", "HIGH"));
    assert!(!accepts(&e, Some(&a), "IMPACT", "HIGH", "LOW"));
    assert!(!accepts(&e, Some(&a), "IMPACT", "MODERATE", "MODIFIER"));
    assert!(!accepts(&e, Some(&a), "IMPACT", "HIGH", "HIGH"));
    assert!(!accepts(&e, Some(&a), "STRAND", "1", "-1"));
    assert!(!accepts(
        &e,
        Some(&a),
        "most_severe_consequence",
        "feature_elongation",
        SAME
    ));
}

/// The masks help states the class and the row it takes.
#[test]
fn the_masks_help_states_the_class() {
    let help = include_str!("../data/help.txt");
    let help: Vec<&str> = help.split_whitespace().collect();
    let stated = format!(
        "{CLASS} (the Consequence of a transcript row whose Allele is tandem_repeat on both \
         sides and whose two lists of terms are equal but for the reference's \
         feature_elongation, where the candidate writes feature_truncation or no term; where \
         it writes none, the IMPACT each side's terms give as well; the JSON record's \
         most_severe_consequence follows), "
    );
    assert!(help.join(" ").contains(&stated));
}
