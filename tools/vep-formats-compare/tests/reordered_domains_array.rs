// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The JSON `domains` arrays of one `DOMAINS` list written in two orders, under `order_only`.
//!
//! Ensembl VEP's JSON writer re-splits each entry of the list with `(\w+):(\w+)` and pushes
//! the two captures (JSON.pm:363-371). An entry the pattern does not match pushes the captures
//! of the last entry that matched once more, and nothing where none has matched yet, so the
//! two orders of one list give arrays that hold the same objects, one of them a different
//! number of times. Every record here is made up.

mod compare_common;

use compare_common::*;
use vep_formats_compare::classes::ORDER_ONLY;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::lists::order_only_fields;
use vep_formats_compare::manifest::Documented;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::Side;

const RECORD: &str = "21\t25587000\tmade_up\tG\tA\t.\t.\t.";

const PFAM: &str = r#"{"db": "Pfam", "name": "PF99001"}"#;
const PANTHER: &str = r#"{"db": "PANTHER", "name": "PTHR99001"}"#;
const MAPPING: &str = r#"{"db": "ENSP_mappings", "name": "AF"}"#;
const GENE3D: &str = r#"{"db": "Gene3D", "name": 1}"#;

/// The entries of one JSON object of `RECORD` whose one transcript consequence carries
/// `domains`.
fn entries(domains: &[&str]) -> Vec<Entry> {
    let line = format!(
        "{{\"input\": \"{}\", \"id\": \"made_up\", \"seq_region_name\": \"21\", \
         \"start\": 25587000, \"end\": 25587000, \"allele_string\": \"G/A\", \"strand\": 1, \
         \"most_severe_consequence\": \"missense_variant\", \"transcript_consequences\": \
         [{{\"variant_allele\": \"A\", \"transcript_id\": \"ENST99000000001\", \
         \"gene_id\": \"ENSG99000000001\", \"consequence_terms\": [\"missense_variant\"], \
         \"impact\": \"MODERATE\", \"strand\": 1, \"domains\": [{}]}}]}}\n",
        RECORD.replace('\t', "\\t"),
        domains.join(", ")
    );
    parse_json(&line, &[RECORD.to_string()]).expect("the made-up object parses")
}

/// The report of two arrays under the fields the JSON format compares under `order_only`,
/// with the classes the comparison applied.
fn under_order_only(reference: &[&str], candidate: &[&str]) -> (Value, Documented) {
    let side = Side::default();
    let lists = order_only_fields("json");
    let options = CompareOptions {
        order_only: Some(&lists),
        ..CompareOptions::new(&side)
    };
    let mut documented = no_documented();
    let report = compare_entries(
        &entries(reference),
        &entries(candidate),
        &mut documented,
        &options,
    )
    .expect("the entries compare");
    (report, documented)
}

fn accepted(reference: &[&str], candidate: &[&str]) {
    let plain = compare(
        &entries(reference),
        &entries(candidate),
        &mut no_documented(),
    );
    assert_eq!(verdict(&plain), "red");
    assert_eq!(number(&plain, &["columns", "domains", "differing"]), 1);
    let (report, documented) = under_order_only(reference, candidate);
    assert_eq!(verdict(&report), "green", "{:?}", reasons(&report));
    assert_eq!(number(&report, &["columns", "domains", "differing"]), 0);
    assert_eq!(number(&report, &["columns", "domains", "order_only"]), 1);
    assert_eq!(counts(&documented.applied), count_map(&[(ORDER_ONLY, 1)]));
}

fn differing(reference: &[&str], candidate: &[&str]) {
    let (report, documented) = under_order_only(reference, candidate);
    assert_eq!(verdict(&report), "red");
    assert_eq!(number(&report, &["columns", "domains", "differing"]), 1);
    assert!(counts(&documented.applied).is_empty());
}

/// The list `Coiled-coils_(Ncoils):Coil,Pfam:PF99001,PANTHER:PTHR99001,
/// AFDB-ENSP_mappings:AF-P99001-F1` opens with the entry the pattern does not match, which
/// pushes nothing; the same four entries as `PANTHER:PTHR99001,Pfam:PF99001,
/// Coiled-coils_(Ncoils):Coil,AFDB-ENSP_mappings:AF-P99001-F1` repeat the Pfam object for it.
/// The two arrays have three and four objects.
#[test]
fn an_unmatched_entry_that_stands_first_in_one_order_alone_leaves_arrays_of_two_lengths() {
    accepted(&[PFAM, PANTHER, MAPPING], &[PANTHER, PFAM, PFAM, MAPPING]);
    accepted(&[PANTHER, PFAM, PFAM, MAPPING], &[PFAM, PANTHER, MAPPING]);
}

/// `Pfam:PF99001,Coiled-coils_(Ncoils):Coil,PANTHER:PTHR99001` repeats the Pfam object and
/// `PANTHER:PTHR99001,Coiled-coils_(Ncoils):Coil,Pfam:PF99001` the PANTHER one: two arrays of
/// three objects that repeat different ones. A name of digits alone is a number in both.
#[test]
fn an_unmatched_entry_after_another_entry_in_each_order_repeats_another_object() {
    accepted(&[PFAM, PFAM, PANTHER], &[PANTHER, PANTHER, PFAM]);
    accepted(&[GENE3D, GENE3D, PFAM], &[PFAM, GENE3D]);
}

/// The same objects the same number of times in another order agree as they did.
#[test]
fn the_same_objects_in_another_order_still_agree() {
    accepted(&[PFAM, PANTHER, MAPPING], &[MAPPING, PFAM, PANTHER]);
    accepted(&[PFAM, PFAM, PANTHER], &[PANTHER, PFAM, PFAM]);
}

/// No order of one list gives these pairs. An object that opens two runs on one side comes
/// from two matched entries, so the other side holds it twice at least; an object one side
/// alone holds is an entry the other list lacks; and a list with a matched entry gives no
/// empty array.
#[test]
fn arrays_that_no_one_list_gives_stay_a_difference() {
    differing(&[PFAM, PANTHER, PFAM], &[PFAM, PANTHER]);
    differing(&[PFAM, PANTHER], &[PANTHER, PFAM, PANTHER, PFAM, PANTHER]);
    differing(&[PFAM, PANTHER, MAPPING], &[PFAM, PANTHER, PANTHER]);
    differing(&[PFAM, PFAM], &[PFAM, PANTHER]);
    differing(&[PFAM], &[]);
    differing(&[], &[PFAM, PFAM]);
}
