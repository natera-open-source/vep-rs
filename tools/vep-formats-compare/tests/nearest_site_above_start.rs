// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The class `nearest_start_site_one_base_above` through the comparison of two outputs: a
//! `NEAREST` value on a row one base below a start site, in the default, tab, VCF and JSON
//! formats, and on a row whose first base is a start site, in the default and VCF formats,
//! with the start sites of a planted table.

mod compare_common;

use std::collections::BTreeMap;
use std::collections::BTreeSet;

use compare_common::*;
use indexmap::IndexMap;
use vep_formats_compare::classes::NEAREST_START_SITE_ONE_BASE_ABOVE;
use vep_formats_compare::classes::ORDER_ONLY;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::nearest::NearestSites;
use vep_formats_compare::parse::json::parse_json;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::parse::text::parse_tab;
use vep_formats_compare::parse::vcf::parse_vcf;
use vep_formats_compare::rules::Side;

const KLASS: &str = NEAREST_START_SITE_ONE_BASE_ABOVE;

/// Sites on sequence 21: two transcripts of `GENE_AB` at 1000, `R` of `GENE_R` at 5000 and
/// `F` of `GENE_F` at 9000.
fn side() -> Side {
    let named: IndexMap<String, BTreeSet<i64>> = [
        ("A", 1000),
        ("B", 1000),
        ("GENE_AB", 1000),
        ("R", 5000),
        ("GENE_R", 5000),
        ("F", 9000),
        ("GENE_F", 9000),
    ]
    .iter()
    .map(|(name, site)| (name.to_string(), BTreeSet::from([*site])))
    .collect();
    Side {
        nearest_sites: Some(NearestSites::on_one_sequence("21", &named)),
        ..Side::default()
    }
}

/// The verdict and what was applied of two entry lists under `classes`, with the start
/// sites beside them or without.
fn under(
    classes: &[&str],
    expected: &[Entry],
    actual: &[Entry],
    side: &Side,
) -> (String, BTreeMap<String, i64>) {
    let mut documented = documented_with(classes);
    let report = compare_entries(
        expected,
        actual,
        &mut documented,
        &CompareOptions::new(side),
    )
    .expect("the entries compare");
    (verdict(&report).to_string(), counts(&documented.applied))
}

/// A default-format output of one upstream row at `location` whose `NEAREST` is `nearest`.
fn default_entries(location: &str, nearest: &str) -> Vec<Entry> {
    let text = default_text(&[default_row(
        "record",
        location,
        "G",
        "R",
        "upstream_gene_variant",
        &format!("IMPACT=MODIFIER;DISTANCE=1;STRAND=1;NEAREST={nearest}"),
    )]);
    parse_default(&text).expect("the output parses")
}

/// A record one base below the site at 5000: the reference names the site its searches
/// reach with that site passed over, 1000, and the candidate the site at 5000. The pair is
/// taken under the class with the start sites, once for the row, and is a difference
/// without the class, without the sites and under another class.
#[test]
fn a_value_one_base_below_a_site_is_taken_under_the_class() {
    let side = side();
    let reference = default_entries("21:4999", "GENE_AB");
    let candidate = default_entries("21:4999", "GENE_R");
    assert_eq!(
        under(&[KLASS], &reference, &candidate, &side),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(under(&[], &reference, &candidate, &side).0, "red");
    assert_eq!(under(&[ORDER_ONLY], &reference, &candidate, &side).0, "red");
    assert_eq!(
        under(&[KLASS], &reference, &candidate, &Side::default()).0,
        "red"
    );
    // One transcript of the two at 1000 against the transcript at 5000.
    let reference = default_entries("21:4999", "B");
    let candidate = default_entries("21:4999", "R");
    assert_eq!(
        under(&[KLASS], &reference, &candidate, &side),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
}

/// Sites on sequence 21 two bases and one base apart, one transcript a site, the first
/// three of `GENE_LMH`: `LOW` at 5000, `MID` at 5002, `HIGH` at 5003 and `F` at 9000.
fn close_sites() -> Side {
    let named: IndexMap<String, BTreeSet<i64>> = [
        ("LOW", vec![5000]),
        ("MID", vec![5002]),
        ("HIGH", vec![5003]),
        ("GENE_LMH", vec![5000, 5002, 5003]),
        ("F", vec![9000]),
        ("GENE_F", vec![9000]),
    ]
    .iter()
    .map(|(name, sites)| (name.to_string(), sites.iter().copied().collect()))
    .collect();
    Side {
        nearest_sites: Some(NearestSites::on_one_sequence("21", &named)),
        ..Side::default()
    }
}

/// A record whose first base is the site at 5000 and whose last is the site at 5002: the
/// reference's searches both reach 5002, since the one above the start begins at 5002, and
/// the candidate's reach both sites at distance zero. A record one base longer covers 5002
/// one base below its end and 5003 on it: the reference names the site on the end alone,
/// the candidate that one and the site on the first base. Each pair is taken under the
/// class with the start sites, and is a difference without the class and without the sites;
/// the gene the three sites share is named by both values, which is no difference.
#[test]
fn a_value_on_a_record_whose_first_base_is_a_site_is_taken_under_the_class() {
    let side = close_sites();
    for (location, reference, candidate) in [
        ("21:5000-5002", "MID", "MID,LOW"),
        ("21:5000-5003", "HIGH", "HIGH,LOW"),
        ("21:5000-5004", "HIGH", "LOW"),
    ] {
        let reference = default_entries(location, reference);
        let candidate = default_entries(location, candidate);
        assert_eq!(
            under(&[KLASS], &reference, &candidate, &side),
            ("green".to_string(), count_map(&[(KLASS, 1)])),
            "{location}"
        );
        assert_eq!(
            under(&[], &reference, &candidate, &side).0,
            "red",
            "{location}"
        );
        assert_eq!(
            under(&[KLASS], &reference, &candidate, &Side::default()).0,
            "red",
            "{location}"
        );
        assert_eq!(
            under(&[KLASS], &candidate, &reference, &side).0,
            "red",
            "{location}"
        );
    }
    let (verdict, applied) = under(
        &[KLASS],
        &default_entries("21:5000-5004", "GENE_LMH"),
        &default_entries("21:5000-5004", "GENE_LMH"),
        &side,
    );
    assert_eq!(verdict, "green");
    assert!(applied.is_empty());
    // An insertion between the sites at 5002 and 5003, as a VCF line states it: the entry
    // is placed on the two positions beside the inserted bases.
    let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations \
                  from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|NEAREST\">\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let record = "21\t5002\trecord\tA\tAGG";
    let input = format!("{record}\t.\t.\t.");
    let vcf = |nearest: &str| -> Vec<Entry> {
        let csq = format!("CSQ=GG|upstream_gene_variant|Transcript|F|{nearest}");
        parse_vcf(
            &format!("{header}\n{record}\t.\t.\t{csq}\n"),
            std::slice::from_ref(&input),
        )
        .expect("the output parses")
    };
    assert_eq!(
        under(&[KLASS], &vcf("HIGH"), &vcf("HIGH&MID"), &side),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(under(&[], &vcf("HIGH"), &vcf("HIGH&MID"), &side).0, "red");
    assert_eq!(
        under(&[KLASS], &vcf("HIGH&MID"), &vcf("HIGH"), &side).0,
        "red"
    );
}

/// The class takes no other difference of the field: a record two bases below the site, a
/// reference value that names a site its searches do not reach, a candidate value that
/// names another site than the one above the start, and the two values exchanged.
#[test]
fn another_nearest_difference_stays_a_difference_under_the_class() {
    let side = side();
    for (location, reference, candidate) in [
        ("21:4998", "GENE_AB", "GENE_R"),
        ("21:4999", "GENE_F", "GENE_R"),
        ("21:4999", "GENE_AB", "GENE_F"),
        ("21:4999", "GENE_R", "GENE_AB"),
        ("21:4999", "GENE_AB", "GENE_R,GENE_F"),
        ("22:4999", "GENE_AB", "GENE_R"),
        ("21:4999-5000", "GENE_AB", "GENE_R"),
    ] {
        let (verdict, applied) = under(
            &[KLASS],
            &default_entries(location, reference),
            &default_entries(location, candidate),
            &side,
        );
        assert_eq!(verdict, "red", "{location} {reference} {candidate}");
        assert!(applied.is_empty(), "{location} {reference} {candidate}");
    }
    // Values that agree by their start sites are settled before any class is asked.
    let (verdict, applied) = under(
        &[KLASS],
        &default_entries("21:4999", "A"),
        &default_entries("21:4999", "A,B"),
        &side,
    );
    assert_eq!(verdict, "green");
    assert!(applied.is_empty());
}

/// The tab format carries the value in a column of its own and the JSON format in the
/// record's `nearest` array: the class takes both.
#[test]
fn the_class_takes_the_tab_column_and_the_json_array() {
    let side = side();
    let tab = |nearest: &str| -> Vec<Entry> {
        let header = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t\
                      Consequence\tNEAREST";
        let row =
            format!("record\t21:4999\tG\tENSG1\tR\tTranscript\tupstream_gene_variant\t{nearest}");
        parse_tab(&format!("{header}\n{row}\n")).expect("the output parses")
    };
    assert_eq!(
        under(&[KLASS], &tab("GENE_AB"), &tab("GENE_R"), &side),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(under(&[], &tab("GENE_AB"), &tab("GENE_R"), &side).0, "red");

    let input = "21\t4999\trecord\tA\tG\t.\t.\t.";
    let json = |nearest: &str| -> Vec<Entry> {
        let text = dumped(&format!(
            r#"{{"input": "21\t4999\trecord\tA\tG\t.\t.\t.", "id": "record", "seq_region_name": "21", "start": 4999, "end": 4999, "allele_string": "A/G", "most_severe_consequence": "upstream_gene_variant", "nearest": ["{nearest}"], "transcript_consequences": [
                {{"variant_allele": "G", "transcript_id": "R", "consequence_terms": ["upstream_gene_variant"]}}]}}"#
        )) + "\n";
        parse_json(&text, &[input.to_string()]).expect("the text parses")
    };
    assert_eq!(
        under(&[KLASS], &json("GENE_AB"), &json("GENE_R"), &side),
        ("green".to_string(), count_map(&[(KLASS, 1)]))
    );
    assert_eq!(
        under(&[], &json("GENE_AB"), &json("GENE_R"), &side).0,
        "red"
    );
}

/// The VCF format carries the value in a CSQ field of the record's line, which states no
/// Location: the class places the entry by the Location the line's POS, REF and ALT give, a
/// deletion's one base above its POS.
#[test]
fn the_class_takes_the_csq_field_of_a_vcf_line() {
    let side = side();
    let header = "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations \
                  from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature|NEAREST\">\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    for (record, allele) in [
        ("21\t4999\trecord\tA\tG", "G"),
        ("21\t4998\trecord\tCA\tC", "-"),
    ] {
        let input = format!("{record}\t.\t.\t.");
        let vcf = |nearest: &str| -> Vec<Entry> {
            let csq = format!("CSQ={allele}|upstream_gene_variant|Transcript|R|{nearest}");
            parse_vcf(
                &format!("{header}\n{record}\t.\t.\t{csq}\n"),
                std::slice::from_ref(&input),
            )
            .expect("the output parses")
        };
        assert_eq!(
            under(&[KLASS], &vcf("GENE_AB"), &vcf("GENE_R"), &side),
            ("green".to_string(), count_map(&[(KLASS, 1)])),
            "{record}"
        );
        assert_eq!(
            under(&[], &vcf("GENE_AB"), &vcf("GENE_R"), &side).0,
            "red",
            "{record}"
        );
        assert_eq!(
            under(&[KLASS], &vcf("GENE_R"), &vcf("GENE_AB"), &side).0,
            "red",
            "{record}"
        );
    }
}
