// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The readers the comparison shares with the tuple comparison
//! (`tools/vep-outputs-compare`), on small cases.

use std::collections::BTreeSet;
use std::panic::{self, AssertUnwindSafe};
use std::path::Path;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::error::{Error, Uncaught};
use vep_formats_compare::shared::*;

/// A Transcript entry of the default format at `location` with the two position columns.
fn row(location: &str, allele: &str, consequence: &str, cdna: &str, cds: &str) -> Entry {
    let mut e = Entry {
        key: format!("{location}|{allele}|ENST00000000001|Transcript"),
        location: Some(location.to_string()),
        allele: allele.to_string(),
        consequence: consequence.to_string(),
        ..Entry::default()
    };
    e.fields
        .insert("cDNA_position".to_string(), cdna.to_string());
    e.fields.insert("CDS_position".to_string(), cds.to_string());
    e
}

fn strings(items: &[&str]) -> Vec<String> {
    items.iter().map(|item| item.to_string()).collect()
}

#[test]
fn an_empty_value_is_blank_in_every_spelling() {
    for value in ["", "-", "<absent>"] {
        assert!(blank(value), "{value:?}");
    }
    for value in ["0", " ", "--", "absent"] {
        assert!(!blank(value), "{value:?}");
    }
    assert_eq!(dash(""), "-");
    assert_eq!(dash("A"), "A");
}

#[test]
fn consequence_terms_are_sorted_and_named_once() {
    assert_eq!(consequence_set("b&a,b"), "a,b");
    assert_eq!(consequence_set(",&"), "");
    assert_eq!(consequence_set("intron_variant"), "intron_variant");
}

#[test]
fn a_feature_type_goes_to_its_consequence_list() {
    for (feature_type, list) in [
        ("Transcript", "transcript"),
        ("RegulatoryFeature", "regulatory_feature"),
        ("MotifFeature", "motif_feature"),
        ("regulatory_feature", "regulatory_feature"),
        ("", "intergenic"),
        ("-", "intergenic"),
    ] {
        assert_eq!(feature_kind(feature_type), list, "{feature_type:?}");
    }
}

#[test]
fn a_generated_name_spells_the_position_and_the_alleles() {
    assert_eq!(
        uploaded_alleles("21_10542414_A/T/AG"),
        Some(("A".to_string(), strings(&["T", "AG"])))
    );
    assert_eq!(uploaded_alleles("rs123"), None);
    assert_eq!(
        uploaded_alleles("21_10542414_A/T\n"),
        Some(("A".to_string(), strings(&["T"])))
    );
    assert_eq!(uploaded_alleles("21_10542414_A/T\n\n"), None);
    assert_eq!(uploaded_position("21_10542414_A/T/AG"), Ok(Some(10542414)));
    assert_eq!(uploaded_position("rs123"), Ok(None));
}

#[test]
fn a_record_decides_whether_an_allele_changes_the_length() {
    assert_eq!(allele_is_indel("A", &["AT", "G"], "AT", false), Some(true));
    assert_eq!(allele_is_indel("A", &["AT", "G"], "G", false), Some(false));
    assert_eq!(allele_is_indel("TCT", &["TCTC"], "C", false), Some(true));
    assert_eq!(allele_is_indel("A", &["<DEL>"], "deletion", false), None);
    assert_eq!(allele_is_indel("T", &["TT", "TA"], "T", false), Some(true));
    assert_eq!(allele_is_indel("T", &["TT", "TA"], "T", true), None);
    assert!(reference_allele_is_empty("T", &["TT", "TA"], false));
    assert!(!reference_allele_is_empty("T", &["TT", "TA"], true));
    assert!(!reference_allele_is_empty("TC", &["T", "TCA"], false));
}

#[test]
fn a_row_reads_the_bases_its_allele_replaces_from_the_record() {
    assert_eq!(
        replaced_reference(100, "ACG", &["ATG"], Some("21:101"), "T"),
        Ok(Some("C".to_string()))
    );
    assert_eq!(
        replaced_reference(100, "A", &["ATT"], Some("21:100-101"), "TT"),
        Ok(Some(String::new()))
    );
    assert_eq!(
        replaced_reference(100, "ACG", &["AG"], Some("21:101"), "-"),
        Ok(Some("C".to_string()))
    );
    assert_eq!(
        replaced_reference(100, "ACG", &["ATG"], Some("21:102"), "T"),
        Ok(None)
    );
    assert_eq!(
        replaced_reference(100, "ACG", &["ATG"], None, "T"),
        Ok(None)
    );
    assert_eq!(
        replaced_reference(100, "ACG", &["ATG"], Some("21:101"), "<DEL>"),
        Ok(None)
    );
}

/// A row that a run under `--shift_genomic 1` writes where the 3' rule moved its record's
/// deletion or insertion (OutputFactory.pm:1711-1715): the bases at the Location are the
/// record's run turned by the count of places moved, on along the chromosome or back, and
/// a moved insertion's allele is the added bases turned so.
#[test]
fn a_moved_row_reads_the_bases_at_its_location_from_the_record() {
    let bases = |bases: &str| Ok(Some(bases.to_string()));
    // Five bases deleted behind the base at 100: one place on, and six.
    let deleted = |location: &str| moved_reference(100, "GAAAAC", &["G"], Some(location), "-");
    assert_eq!(deleted("21:102-106"), bases("AAACA"));
    assert_eq!(deleted("21:107-111"), bases("AAACA"));
    // The run where the record writes it is read by position.
    assert_eq!(deleted("21:101-105"), Ok(None));
    assert_eq!(
        replaced_reference(100, "GAAAAC", &["G"], Some("21:101-105"), "-"),
        bases("AAAAC")
    );
    // Another count of positions than the run has bases, and an allele that has bases.
    assert_eq!(deleted("21:103-106"), Ok(None));
    assert_eq!(deleted("21:102-107"), Ok(None));
    assert_eq!(
        moved_reference(100, "GAAAAC", &["G"], Some("21:102-106"), "A"),
        Ok(None)
    );
    // An alternate that is not the reference allele less a run of its bases.
    for alternate in ["T", "GAAAAT", "GAAAACT"] {
        assert_eq!(
            moved_reference(100, "GAAAAC", &[alternate], Some("21:102-106"), "-"),
            Ok(None),
            "{alternate}"
        );
    }
    // The base before the run is not the run's last, so the run is not moved back; and a
    // run that bases of the record follow, which do not repeat it, is not moved on.
    assert_eq!(deleted("21:99-103"), Ok(None));
    assert_eq!(deleted("21:96-100"), Ok(None));
    assert_eq!(
        moved_reference(100, "GAAAACT", &["GCT"], Some("21:107-110"), "-"),
        Ok(None)
    );
    // The base before the run is the run's last, so the record's alleles give the run at
    // 100 and at 101: both places are read by position, and the run two places back from
    // the second is turned by two.
    let behind_its_last =
        |location: &str| moved_reference(100, "CAAAAC", &["C"], Some(location), "-");
    assert_eq!(behind_its_last("21:100-104"), Ok(None));
    assert_eq!(behind_its_last("21:101-105"), Ok(None));
    assert_eq!(behind_its_last("21:99-103"), bases("ACAAA"));
    assert_eq!(behind_its_last("21:102-106"), bases("AAACA"));
    // One base of two that are alike, moved along a run of that base.
    assert_eq!(
        moved_reference(100, "GTT", &["GT"], Some("21:153"), "-"),
        bases("T")
    );
    // Three bases added behind the base at 100, which is the last of them. Ten places on
    // the allele is the added bases turned by one; nine on and three back it is the added
    // bases whole; two back and one back it is turned the other way, and one back the
    // Location is the two positions before the reference allele's second.
    let added =
        |location: &str, allele: &str| moved_reference(100, "A", &["ATCA"], Some(location), allele);
    assert_eq!(added("21:110-111", "CAT"), bases(""));
    assert_eq!(added("21:109-110", "TCA"), bases(""));
    assert_eq!(added("21:97-98", "TCA"), bases(""));
    assert_eq!(added("21:98-99", "CAT"), bases(""));
    assert_eq!(added("21:99-100", "ATC"), bases(""));
    // The added bases not turned as the place asks, another count of bases, the place the
    // record writes them at, which is read by position, and a Location of three positions.
    assert_eq!(added("21:110-111", "TCA"), Ok(None));
    assert_eq!(added("21:110-111", "TA"), Ok(None));
    assert_eq!(added("21:99-100", "TCA"), Ok(None));
    assert_eq!(added("21:100-101", "TCA"), Ok(None));
    assert_eq!(
        replaced_reference(100, "A", &["ATCA"], Some("21:100-101"), "TCA"),
        bases("")
    );
    assert_eq!(added("21:110-112", "CAT"), Ok(None));
    // The base before the added bases is not their last, so they are not moved back; they
    // are moved on. Bases of the record after them, which do not repeat them, stop that.
    let behind_another_base =
        |location: &str, allele: &str| moved_reference(100, "G", &["GTCA"], Some(location), allele);
    assert_eq!(behind_another_base("21:110-111", "CAT"), bases(""));
    assert_eq!(behind_another_base("21:97-98", "TCA"), Ok(None));
    assert_eq!(behind_another_base("21:99-100", "ATC"), Ok(None));
    assert_eq!(
        moved_reference(100, "GC", &["GTCAC"], Some("21:110-111"), "CAT"),
        Ok(None)
    );
    // Each alternate of a record is asked by itself.
    let either =
        |allele: &str| moved_reference(100, "G", &["GT", "GTT"], Some("21:153-154"), allele);
    assert_eq!(either("T"), bases(""));
    assert_eq!(either("TT"), bases(""));
    assert_eq!(either("TTT"), Ok(None));
    // A Location on the reference allele's own positions is read by position or not at
    // all: the base after the one an alternate lacks, the place at which another alternate
    // lacks its base, and two bases one alternate replaces by the allele that another adds
    // after them.
    assert_eq!(
        moved_reference(100, "GACT", &["GCT"], Some("21:102"), "-"),
        Ok(None)
    );
    assert_eq!(
        moved_reference(100, "GACT", &["GCT", "GAT"], Some("21:101"), "-"),
        Ok(None)
    );
    assert_eq!(
        replaced_reference(100, "GACT", &["GCT", "GAT"], Some("21:101"), "-"),
        bases("A")
    );
    assert_eq!(
        moved_reference(100, "GAC", &["GT", "GACT"], Some("21:101-102"), "T"),
        Ok(None)
    );
    assert_eq!(
        replaced_reference(100, "GAC", &["GT", "GACT"], Some("21:101-102"), "T"),
        bases("AC")
    );
    // No Location, a Location that names no position, and a structural allele.
    assert_eq!(moved_reference(100, "GAAAAC", &["G"], None, "-"), Ok(None));
    assert_eq!(
        moved_reference(100, "GAAAAC", &["G"], Some("21:x-y"), "-"),
        Ok(None)
    );
    assert_eq!(
        moved_reference(100, "GAAAAC", &["<DEL>"], Some("21:102-106"), "deletion"),
        Ok(None)
    );
}

/// The records of an input VCF give a row its bases by position, and a row that none
/// reads so, which the 3' rule moved off its record's alleles, by the record's identifier,
/// which is the row's first column: on the row's chromosome alone, from no record of
/// another name, and not where two records of the name read the row differently.
#[test]
fn the_record_a_moved_row_names_gives_its_bases() {
    let directory =
        Path::new(env!("CARGO_TARGET_TMPDIR")).join("vep-formats-compare-shared-moved-rows");
    std::fs::create_dir_all(&directory).expect("the directory is made");
    let input = directory.join("input.vcf");
    let lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "21\t100\tdeleted\tGAAAAC\tG\t.\t.\t.",
        "21\t200\tadded\tG\tGTCA\t.\t.\t.",
        "21\t300\tbehind_its_last\tCAAAAC\tC\t.\t.\t.",
        "21\t400\ttwice\tGAAAAC\tG\t.\t.\t.",
        "21\t420\ttwice\tGACCCC\tG\t.\t.\t.",
        "21\t600\t.\tGAAAAC\tG\t.\t.\t.",
        "21\t700\ttwo_runs\tGACT\tGCT,GAT\t.\t.\t.",
        "21\t750\tshared\tGC\tG\t.\t.\t.",
        "21\t800\tshared\tGA\tG\t.\t.\t.",
        "22\t100\telsewhere\tGAAAAC\tG\t.\t.\t.",
    ];
    std::fs::write(&input, lines.join("\n") + "\n").expect("the input is written");
    let named = |name: &str, location: &str, allele: &str| {
        let mut e = Entry {
            key: format!("{location}|{allele}|ENST00000000001|Transcript"),
            location: Some(location.to_string()),
            allele: allele.to_string(),
            ..Entry::default()
        };
        e.fields
            .insert("Uploaded_variation".to_string(), name.to_string());
        e
    };
    let cases = [
        // Moved on and back, with and without the `chr` prefix, and an insertion.
        (named("deleted", "21:102-106", "-"), Some("AAACA")),
        (named("deleted", "chr21:107-111", "-"), Some("AAACA")),
        (named("behind_its_last", "21:299-303", "-"), Some("ACAAA")),
        (named("added", "21:210-211", "CAT"), Some("")),
        // Where the record writes it, by position.
        (named("deleted", "21:101-105", "-"), Some("AAAAC")),
        // Each of two alternates that lack one base, at its own place, by position.
        (named("two_runs", "21:701", "-"), Some("A")),
        (named("two_runs", "21:702", "-"), Some("C")),
        // Two records of one name: the row one of them reads by position keeps that
        // reading, whatever the other gives it moved, and a row neither reads by position
        // is the moved row of the one whose own bases let it move there.
        (named("shared", "21:801", "-"), Some("A")),
        (named("shared", "21:761", "-"), Some("C")),
        // Another name, a record without an identifier, and another chromosome.
        (named("another", "21:102-106", "-"), None),
        (named(".", "21:602-606", "-"), None),
        (named("deleted", "22:102-106", "-"), None),
        (named("elsewhere", "21:102-106", "-"), None),
        // Two records of one name that both move to the row and whose runs differ.
        (named("twice", "21:450-454", "-"), None),
        // An insertion's allele that is not the added bases turned.
        (named("added", "21:210-211", "TCA"), None),
    ];
    let mut rows: Vec<Entry> = cases.iter().map(|(row, _)| row.clone()).collect();
    read_record_references(&input, &mut rows).expect("the input is read");
    let read: Vec<(&str, Option<&str>)> = rows
        .iter()
        .map(|row| (row.key.as_str(), row.reference.as_deref()))
        .collect();
    let owed: Vec<(&str, Option<&str>)> = cases
        .iter()
        .map(|(row, bases)| (row.key.as_str(), *bases))
        .collect();
    assert_eq!(read, owed);
}

#[test]
fn digits_are_read_as_a_whole_number_or_end_the_run() {
    assert_eq!(integer_of("0012", "POS"), Ok(12));
    assert_eq!(
        integer_of("\u{b2}", "SHIFT_LENGTH"),
        Err(Error::comparator("SHIFT_LENGTH '\u{b2}' is not an integer"))
    );
    let long = "7".repeat(5000);
    assert_eq!(
        integer_of(&long, "SHIFT_LENGTH"),
        Err(Error::comparator(format!(
            "SHIFT_LENGTH '{}...' is not an integer",
            "7".repeat(24)
        )))
    );
}

#[test]
fn a_span_is_written_low_position_first() {
    assert_eq!(sequence_location("21", 5, 4), "21:4-5");
    assert_eq!(sequence_location("21", 5, 5), "21:5");
    assert_eq!(sequence_location("21", 5, 9), "21:5-9");
    assert_eq!(
        canonical_location("21:30425651-030425949"),
        "21:30425651-30425949"
    );
    assert_eq!(
        canonical_location("21:30425651-+30425949"),
        "21:30425651-30425949"
    );
    assert_eq!(
        canonical_location("21:30425651-30425949.5"),
        "21:30425651-30425949.5"
    );
    assert_eq!(canonical_location("no_colon"), "no_colon");
}

#[test]
fn a_list_is_split_as_the_reference_splits_it() {
    assert_eq!(perl_split("G,,T,,", ","), strings(&["G", "", "T"]));
    assert_eq!(perl_split("<DEL>,", ","), strings(&["<DEL>"]));
    assert_eq!(perl_split("", ","), Vec::<String>::new());
}

#[test]
fn a_position_column_is_read_as_the_candidate_reads_it() {
    for (text, position) in [
        ("3e7", Some(30000000)),
        ("30425650.0", Some(30425650)),
        ("30425650.", Some(30425650)),
        ("+30425650", Some(30425650)),
        (" 30425650 ", Some(30425650)),
        (" 72057594037927937 ", Some(72057594037927936)),
        ("4611686018427387903", Some(4611686018427387903)),
        ("4611686018427387904", None),
        (".", None),
        ("abc", None),
        ("-5", None),
        ("0x1D04A1A", None),
        ("30_425_650", None),
        ("\u{ff11}\u{ff12}", None),
        ("inf", None),
        ("1e400", None),
        ("30425650.5", None),
    ] {
        assert_eq!(vep_rs_pos(text), position, "{text:?}");
    }
}

#[test]
fn a_number_past_i64_is_the_comparators_error_or_the_nearest_number() {
    let past = "9223372036854775808";
    assert_eq!(
        integer_of(past, "Location position"),
        Err(Error::comparator(
            "Location position '9223372036854775808' is not an integer"
        ))
    );
    assert_eq!(
        integer_text_of(past, "Location position"),
        Ok(past.to_string())
    );
    assert_eq!(integer_text_of("007", "SHIFT_LENGTH"), Ok("7".to_string()));
    let e = row("21:5-9223372036854775808", "A", "", "-", "-");
    assert_eq!(
        indel_row(&e),
        Err(Error::comparator(
            "Location position '9223372036854775808' is not an integer"
        ))
    );
    let whole = row("21:0-9223372036854775807", "A", "", "-", "-");
    assert_eq!(location_span(&whole), Ok(Some(i64::MAX)));
    assert_eq!(
        location_span(&row("21:5-7", "A", "", "-", "-")),
        Ok(Some(3))
    );
    assert_eq!(location_span(&row("21:7-5", "A", "", "-", "-")), Ok(None));
}

#[test]
fn digits_the_interpreter_does_not_convert_end_the_run_where_nothing_catches_them() {
    let long = "9".repeat(4301);
    let e = row(&format!("21:5-{long}"), "ACGT", "", "-", "-");
    let ended = panic::catch_unwind(AssertUnwindSafe(|| insertion_row(&e)));
    let what = ended.expect_err("the run ends");
    assert_eq!(
        what.downcast_ref::<Uncaught>(),
        Some(&Uncaught(
            "Exceeds the limit (4300 digits) for integer string conversion: value has 4301 \
             digits; use sys.set_int_max_str_digits() to increase the limit"
                .to_string()
        ))
    );
    let digits = "9".repeat(4300);
    assert!(!insertion_row(&row(
        &format!("21:5-{digits}"),
        "ACGT",
        "",
        "-",
        "-"
    )));
}

#[test]
fn a_stop_codon_rule_takes_its_pair_of_rows_and_no_other() {
    let coding = "coding_sequence_variant&splice_region_variant";
    let retained = "stop_retained_variant,splice_region_variant";
    let e = row("21:100-101", "-", coding, "100-101", "90-91");
    let a = row("21:100-101", "-", retained, "100-101", "90-91");
    assert_eq!(stop_codon_overlap_genomic_span_pair(&e, &a), Ok(true));
    assert_eq!(stop_codon_overlap_genomic_span_pair(&a, &e), Ok(false));
    let away = row(
        "21:100-101",
        "-",
        "coding_sequence_variant",
        "100-101",
        "90-91",
    );
    let kept = row(
        "21:100-101",
        "-",
        "stop_retained_variant",
        "100-101",
        "90-91",
    );
    assert_eq!(
        stop_codon_overlap_genomic_span_pair(&away, &kept),
        Ok(false)
    );

    let untested = |consequence: &str| row("21:100-103", "-", consequence, "100-?", "90-93");
    let e = untested("stop_retained_variant,3_prime_UTR_variant");
    assert_eq!(
        stop_retained_without_stop_codon_test_pair(&e, &untested("3_prime_UTR_variant")),
        Ok(true)
    );
    let tested = |consequence: &str| row("21:100-103", "-", consequence, "100-103", "90-93");
    assert_eq!(
        stop_retained_without_stop_codon_test_pair(
            &tested("stop_retained_variant"),
            &tested("coding_sequence_variant")
        ),
        Ok(false)
    );
    assert_eq!(
        stop_retained_three_base_stop_codon_test_pair(
            &tested("stop_retained_variant"),
            &tested("coding_sequence_variant")
        ),
        Ok(true)
    );
    assert_eq!(
        stop_retained_three_base_stop_codon_test_pair(
            &tested("stop_retained_variant,inframe_deletion"),
            &tested("frameshift_variant,stop_gained")
        ),
        Ok(true)
    );
    let without_start = |consequence: &str| row("21:100-103", "-", consequence, "100-103", "?-93");
    assert_eq!(
        stop_codon_removed_without_stop_codon_test_pair(
            &without_start("stop_retained_variant"),
            &without_start("stop_lost")
        ),
        Ok(true)
    );
    assert_eq!(
        stop_codon_removed_without_stop_codon_test_pair(
            &tested("stop_retained_variant"),
            &tested("stop_lost")
        ),
        Ok(false)
    );
    let inserted = |consequence: &str| row("21:100-101", "A", consequence, "100-101", "90-91");
    assert_eq!(
        stop_codon_overlap_insertion_range_pair(
            &inserted("coding_sequence_variant"),
            &inserted("stop_retained_variant")
        ),
        Ok(true)
    );
    assert_eq!(
        stop_codon_overlap_insertion_range_pair(
            &inserted("frameshift_variant"),
            &inserted("stop_retained_variant,inframe_insertion")
        ),
        Ok(true)
    );
    let two_bases = |consequence: &str| row("21:100-101", "AC", consequence, "100-101", "90-91");
    assert_eq!(
        stop_codon_overlap_insertion_range_pair(
            &two_bases("coding_sequence_variant"),
            &two_bases("stop_retained_variant")
        ),
        Ok(false)
    );
}

#[test]
fn a_splice_term_past_the_span_is_taken_inside_one_exon_without_the_record() {
    let inside = |consequence: &str| row("7:1000-1002", "ACGTA", consequence, "500-502", "-");
    let e = inside("splice_donor_variant,3_prime_UTR_variant");
    let a = inside("splice_region_variant,3_prime_UTR_variant");
    assert_eq!(splice_term_past_the_record_span_pair(&e, &a), Ok(true));
    assert_eq!(splice_term_past_the_record_span_pair(&a, &e), Ok(false));
    let across = |consequence: &str| row("7:1000-1002", "ACGTA", consequence, "500-501", "-");
    assert_eq!(
        splice_term_past_the_record_span_pair(
            &across("splice_donor_variant,3_prime_UTR_variant"),
            &across("splice_region_variant,3_prime_UTR_variant")
        ),
        Ok(false)
    );
}

#[test]
fn the_excluding_rules_are_the_start_pair_and_the_covered_splice_region() {
    let classes: Vec<&str> = EXCLUDING_RULES
        .iter()
        .map(|rule| rule.taxonomy_class)
        .collect();
    assert_eq!(classes, ["start_cooccurrence_swap", "splice_family_swap"]);
    let both = "start_lost,start_retained_variant";
    assert!(EXCLUDING_RULES[0].matches(both, "start_retained_variant", "A"));
    assert!(!EXCLUDING_RULES[0].matches(both, "start_retained_variant", "deletion"));
    assert!(EXCLUDING_RULES[0].matches(both, "start_lost", "deletion"));
    assert!(EXCLUDING_RULES[1].matches(
        "intron_variant",
        "intron_variant,splice_region_variant",
        "A"
    ));
    assert!(!EXCLUDING_RULES[1].matches(
        "intron_variant",
        "intron_variant,splice_polypyrimidine_tract_variant",
        "A"
    ));
}

#[test]
fn a_generated_name_spells_the_alleles_and_they_give_each_allele_its_shape() {
    let named = |name: &str| uploaded_alleles(name);
    assert_eq!(
        named("21_10542414_A/T/AG"),
        Some(("A".to_string(), strings(&["T", "AG"])))
    );
    assert_eq!(named("rs123"), None);
    assert_eq!(
        named("chr1_KI270706v1_random_100_A/AAC"),
        Some(("A".to_string(), strings(&["AAC"])))
    );
    assert_eq!(named("del_padded"), None);
    let shape = |reference: &str, alternates: &[&str], allele: &str| {
        allele_is_indel(reference, alternates, allele, false)
    };
    assert_eq!(shape("A", &["AT"], "T"), Some(true));
    assert_eq!(shape("A", &["G"], "G"), Some(false));
    assert_eq!(shape("AT", &["A"], "-"), Some(true));
    assert_eq!(shape("A", &["<DEL>"], "deletion"), None);
    assert_eq!(shape("G", &["[21:25881361[G"], "[21:25881361[G"), None);
    // A lone alternate keeps its shape however Ensembl VEP shortens the allele; the parts of
    // a JSON allele_string are read as they are.
    assert_eq!(shape("TTTTT", &["TTTTTT"], "T"), Some(true));
    assert_eq!(shape("TCT", &["TCTC"], "C"), Some(true));
    assert_eq!(allele_is_indel("T", &["TT", "-"], "TT", true), Some(true));
    assert_eq!(allele_is_indel("T", &["TT", "-"], "-", true), Some(true));
    // The shared first base is dropped only where every alternate shares it.
    assert_eq!(shape("A", &["AT", "G"], "AT"), Some(true));
    assert_eq!(shape("A", &["AT", "G"], "G"), Some(false));
    assert_eq!(shape("A", &["AT", "G"], "T"), None);
    assert_eq!(shape("A", &["AT", "T"], "AT"), Some(true));
    assert_eq!(shape("A", &["AT", "T"], "T"), Some(false));
    assert_eq!(shape("A", &["AT", "*"], "T"), Some(true));
    // An allele spelt with a `*` is exempt from the drop whatever else it carries.
    assert_eq!(shape("A", &["AG", "<DEL:*>"], "G"), Some(true));
    assert_eq!(
        allele_is_indel("-", &["G", "<DEL:*>"], "G", true),
        Some(true)
    );
    assert_eq!(shape("A", &["T", "<DEL:*>"], "T"), Some(false));
    assert_eq!(shape("AC", &["ACT", "AG"], "CT"), Some(true));
    assert_eq!(shape("AC", &["ACT", "AG"], "G"), Some(false));
}

#[test]
fn a_record_without_an_indel_keeps_its_shared_first_base() {
    let shape = |reference: &str, alternates: &[&str], allele: &str| {
        allele_is_indel(reference, alternates, allele, false)
    };
    assert_eq!(shape("AC", &["AT", "AG"], "AT"), Some(false));
    assert_eq!(shape("AC", &["AT", "AG"], "AG"), Some(false));
    assert_eq!(shape("AC", &["AT", "A"], "T"), Some(false));
    assert_eq!(shape("AC", &["AT", "A"], "-"), Some(true));
    assert_eq!(shape("GC", &["GT", "GA", "*"], "T"), Some(false));
}

/// A Transcript entry whose allele stands in for the reference bases `reference`.
fn record_row(reference: &str, location: &str, allele: &str, terms: &str, cdna: &str) -> Entry {
    Entry {
        reference: Some(reference.to_string()),
        ..row(location, allele, terms, cdna, "-")
    }
}

/// The terms of a comma-separated list.
fn named(terms: &str) -> BTreeSet<String> {
    terms
        .split(',')
        .filter(|term| !term.is_empty())
        .map(str::to_string)
        .collect()
}

/// A span of two bases beside one intron, under an allele that opens with them and holds
/// further bases: the reference's terms come from the bases after the span, the candidate's
/// from the insertion between the span's last base and the base after it.
#[test]
fn an_insertion_after_the_span_is_read_by_its_two_neighbours() {
    let three_further = differing_regions("TA", "TAGTG");
    assert_eq!(three_further, [(2, 4)]);
    let rows = |start: i64, end: i64, forward: bool| {
        intron_rows(&three_further, 2, start, end, forward, false)
    };
    let both = |reference: &str, candidate: &str| (named(reference), named(candidate));
    assert_eq!(
        rows(-15, _FAR, false),
        both(
            "intron_variant",
            "intron_variant,splice_polypyrimidine_tract_variant"
        ),
        "after the tract's last base, on the reverse strand"
    );
    assert_eq!(
        rows(-16, _FAR, false),
        both("intron_variant", "intron_variant"),
        "one base past the tract's last base"
    );
    assert_eq!(
        rows(-_FAR, 3, true),
        both(
            "splice_acceptor_variant",
            "intron_variant,splice_polypyrimidine_tract_variant,splice_region_variant"
        ),
        "between the tract's innermost base and the acceptor bases"
    );
    assert_eq!(
        rows(-_FAR, 3, false),
        both(
            "splice_donor_variant",
            "intron_variant,splice_region_variant"
        ),
        "between the third intronic base and the donor bases, on the reverse strand"
    );
    assert_eq!(
        rows(1, _FAR, true),
        both(
            "splice_donor_variant,splice_donor_region_variant,intron_variant",
            "splice_donor_variant"
        ),
        "between the two donor bases"
    );
    assert_eq!(
        rows(1, _FAR, false),
        both(
            "splice_acceptor_variant,intron_variant",
            "splice_acceptor_variant"
        ),
        "between the two acceptor bases, on the reverse strand"
    );
    let two_further = differing_regions("TA", "TAGT");
    assert_eq!(
        intron_rows(&two_further, 2, 2, _FAR, true, false),
        both("splice_donor_variant", "splice_region_variant"),
        "between an exon's last base and the donor bases"
    );
    assert_eq!(
        intron_rows(&differing_regions("T", "TGTG"), 1, -15, _FAR, false, false),
        both(
            "intron_variant,splice_polypyrimidine_tract_variant",
            "intron_variant,splice_polypyrimidine_tract_variant"
        ),
        "one reference base is read as a change of that base"
    );
    assert_eq!(
        intron_rows(
            &differing_regions("TA", "TGGTG"),
            2,
            -15,
            _FAR,
            false,
            false
        ),
        both(
            "intron_variant,splice_polypyrimidine_tract_variant",
            "intron_variant,splice_polypyrimidine_tract_variant"
        ),
        "an allele that changes the last base of its span is read by its span"
    );
}

/// The reference's terms and the candidate's at the three places where the insertion's own
/// neighbours give the candidate a term the bases after the span do not give the reference.
const AFTER_THE_SPAN_PLACES: [(&str, &str, &str); 3] = [
    (
        "after the tract's last base, on the reverse strand",
        "intron_variant",
        "splice_polypyrimidine_tract_variant,intron_variant",
    ),
    (
        "between the tract's innermost base and the acceptor bases",
        "splice_acceptor_variant",
        "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
    ),
    (
        "between the third intronic base and the donor bases, on the reverse strand",
        "splice_donor_variant",
        "splice_region_variant,intron_variant",
    ),
];

/// Whether the class takes the two rows of the allele `allele` at `location`, whose
/// reference bases are `reference`, their terms as given.
fn after_the_span(
    reference: &str,
    location: &str,
    allele: &str,
    cdna: &str,
    terms: (&str, &str),
) -> bool {
    let row = |terms: &str| record_row(reference, location, allele, terms, cdna);
    splice_term_past_the_record_span_pair(&row(terms.0), &row(terms.1)).expect("no refusal")
}

const TWO_POSITIONS: &str = "7:100999-101000";

#[test]
fn a_splice_term_of_an_insertion_after_the_span_is_taken_with_the_record() {
    for (place, reference, candidate) in AFTER_THE_SPAN_PLACES {
        let terms = (reference, candidate);
        assert!(
            after_the_span("TA", TWO_POSITIONS, "TAGTG", "-", terms),
            "{place}"
        );
        assert!(
            after_the_span("TA", TWO_POSITIONS, "TAG", "-", terms),
            "{place}, one further base"
        );
    }
    let (_, reference, candidate) = AFTER_THE_SPAN_PLACES[0];
    assert!(
        after_the_span(
            "TAC",
            "7:100998-101000",
            "TACGTG",
            "-",
            (reference, candidate)
        ),
        "three reference bases"
    );
    assert!(
        after_the_span(
            "TA",
            TWO_POSITIONS,
            "TAGTG",
            "1000-?",
            (
                "splice_donor_variant,splice_donor_region_variant,coding_sequence_variant,intron_variant",
                "splice_donor_variant,coding_sequence_variant"
            )
        ),
        "between the two donor bases"
    );
    assert!(
        after_the_span(
            "TA",
            TWO_POSITIONS,
            "TAGTG",
            "?-1000",
            (
                "splice_acceptor_variant,coding_sequence_variant,intron_variant",
                "splice_acceptor_variant,coding_sequence_variant"
            )
        ),
        "between the two acceptor bases, on the reverse strand"
    );
}

/// The geometry of line 166 of `tests/golden/116/TAIR10-arabidopsis/variants.vcf`, in
/// made-up bases: three reference bases, a first alternate allele of one base and a second
/// that adds two bases to the reference allele. The second's rows stand at the two positions
/// after the shared first base, an exon's last two bases on one transcript and the two donor
/// bases of an intron on another.
#[test]
fn both_rows_of_a_second_alternate_that_opens_with_its_reference_are_taken() {
    let location = "4:10523300-10523301";
    assert_eq!(
        replaced_reference(10523299, "CTA", &["C", "CTAGT"], Some(location), "TAGT"),
        Ok(Some("TA".to_string()))
    );
    let exon = |terms: (&str, &str)| after_the_span("TA", location, "TAGT", "1444-1445", terms);
    assert!(exon((
        "splice_donor_variant,frameshift_variant",
        "frameshift_variant,splice_region_variant"
    )));
    assert!(
        !exon((
            "splice_donor_variant,frameshift_variant",
            "frameshift_variant"
        )),
        "a candidate's row without the insertion's splice_region_variant"
    );
    assert!(after_the_span(
        "TA",
        location,
        "TAGT",
        "-",
        (
            "splice_donor_region_variant,intron_variant",
            "intron_variant,splice_region_variant"
        )
    ));
    let unread = |terms: &str| row(location, "TAGT", terms, "1444-1445", "1299-1300");
    assert_eq!(
        splice_term_past_the_record_span_pair(
            &unread("splice_donor_variant,frameshift_variant"),
            &unread("frameshift_variant,splice_region_variant")
        ),
        Ok(false),
        "two positions without the record"
    );
}

/// The reference's terms and the candidate's for an allele that opens with its one reference
/// base, at the places where the base read as changed and the insertion after it give
/// different terms, with the cDNA position of the base.
const AFTER_ONE_BASE_PLACES: [(&str, &str, (&str, &str)); 10] = [
    (
        "the base before the tract's outermost base",
        "-",
        (
            "intron_variant",
            "splice_polypyrimidine_tract_variant,intron_variant",
        ),
    ),
    (
        "the intron's last base, before the exon",
        "-",
        ("splice_acceptor_variant", "splice_region_variant"),
    ),
    (
        "the exon's third base",
        "3",
        (
            "splice_region_variant,non_coding_transcript_exon_variant",
            "non_coding_transcript_exon_variant",
        ),
    ),
    (
        "the donor's second base",
        "-",
        (
            "splice_donor_variant",
            "splice_region_variant,intron_variant",
        ),
    ),
    (
        "the intron's fifth base",
        "-",
        (
            "splice_donor_5th_base_variant,intron_variant",
            "splice_donor_region_variant,intron_variant",
        ),
    ),
    (
        "the intron's sixth base",
        "-",
        (
            "splice_donor_region_variant,intron_variant",
            "splice_region_variant,intron_variant",
        ),
    ),
    (
        "the intron's eighth base",
        "-",
        ("splice_region_variant,intron_variant", "intron_variant"),
    ),
    (
        "the acceptor's base beside the tract, on the reverse strand",
        "-",
        (
            "splice_acceptor_variant",
            "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
        ),
    ),
    (
        "the eighth base from the intron's 3' end, on the reverse strand",
        "-",
        (
            "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant",
            "splice_polypyrimidine_tract_variant,intron_variant",
        ),
    ),
    (
        "the donor's first base, before the exon on the reverse strand",
        "-",
        ("splice_donor_variant", "splice_region_variant"),
    ),
];

#[test]
fn a_splice_term_of_an_insertion_after_one_base_is_taken_with_the_record() {
    for (place, cdna, terms) in AFTER_ONE_BASE_PLACES {
        assert!(
            after_the_span("T", "7:101000", "TGTG", cdna, terms),
            "{place}"
        );
        assert!(
            after_the_span("T", "7:101000", "TG", cdna, terms),
            "{place}, one further base"
        );
        assert!(
            !after_the_span("T", "7:101000", "GTGT", cdna, terms),
            "{place}, an allele that does not open with the reference base"
        );
        let unread = |terms: &str| row("7:101000", "TGTG", terms, cdna, "-");
        assert_eq!(
            splice_term_past_the_record_span_pair(&unread(terms.0), &unread(terms.1)),
            Ok(false),
            "{place}, without the record"
        );
    }
}

/// A record of two alternate alleles: the second opens with the reference allele and adds
/// three bases, and its row stands at the one position after the shared first base, the
/// donor's second base.
#[test]
fn the_row_of_an_alternate_that_opens_with_one_kept_base_is_taken() {
    let location = "22:20996968";
    assert_eq!(
        replaced_reference(20996967, "CT", &["C", "CTGTG"], Some(location), "TGTG"),
        Ok(Some("T".to_string()))
    );
    let donor = |terms: (&str, &str)| after_the_span("T", location, "TGTG", "-", terms);
    assert!(donor((
        "splice_donor_variant",
        "splice_region_variant,intron_variant"
    )));
    assert!(
        !donor((
            "splice_donor_variant",
            "splice_donor_region_variant,intron_variant"
        )),
        "a candidate's row with a term no place of the insertion gives beside the donor term"
    );
}

#[test]
fn a_term_an_insertion_after_one_base_does_not_explain_stays_a_difference() {
    let refused = [
        (
            "a donor term the candidate alone holds",
            "-",
            "intron_variant",
            "intron_variant,splice_donor_variant",
        ),
        (
            "another term the candidate alone holds",
            "-",
            "splice_acceptor_variant",
            "splice_region_variant,missense_variant",
        ),
        (
            "another term the reference alone holds",
            "-",
            "splice_donor_variant,non_coding_transcript_variant",
            "splice_region_variant,intron_variant",
        ),
        (
            "the intron term the insertion's place gives, missing after the donor's second base",
            "-",
            "splice_donor_variant,intron_variant",
            "splice_region_variant",
        ),
        (
            "splice_region_variant gone from a base inside the transcript's sequence, with the intron term",
            "3",
            "splice_region_variant,intron_variant",
            "intron_variant",
        ),
        (
            "the same terms on both rows",
            "-",
            "splice_region_variant,intron_variant",
            "splice_region_variant,intron_variant",
        ),
    ];
    for (shape, cdna, reference, candidate) in refused {
        assert!(
            !after_the_span("T", "7:101000", "TGTG", cdna, (reference, candidate)),
            "{shape}"
        );
    }
}

/// The kept base is an intron of one base: the insertion after it stands between the two
/// bases of that intron's first site window, and the base read as changed draws nothing.
#[test]
fn the_site_term_of_an_insertion_after_a_one_base_intron_is_taken_with_the_record() {
    let lone = "non_coding_transcript_variant";
    for site in ["splice_acceptor_variant", "splice_donor_variant"] {
        let candidate = format!("{lone},{site}");
        let rows = (lone, candidate.as_str());
        assert!(after_the_span("T", "7:101000", "TGT", "-", rows), "{site}");
        assert!(
            !after_the_span("T", "7:101000", "TGT", "3", rows),
            "{site} on a base of the transcript's sequence"
        );
        assert!(
            !after_the_span("T", "7:101000", "GGT", "-", rows),
            "{site} under an allele that does not keep the base"
        );
    }
    let refused = |reference: &str, candidate: &str| {
        !after_the_span("T", "7:101000", "TGT", "-", (reference, candidate))
    };
    let site = "splice_acceptor_variant";
    let others_of_the_seven = [
        "splice_donor_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
        "intron_variant",
        "splice_region_variant",
    ];
    for term in others_of_the_seven {
        let held = format!("{lone},{term}");
        assert!(
            refused(&held, &format!("{held},{site}")),
            "a site term beside the reference's {term}"
        );
        assert!(
            refused(&held, &format!("{lone},{site}")),
            "a site term in place of the reference's {term}"
        );
        assert!(
            refused(lone, &format!("{held},{site}")),
            "a site term with {term}"
        );
    }
    for beside in ["upstream_gene_variant", "downstream_gene_variant"] {
        assert!(
            refused(beside, &format!("{beside},{site}")),
            "a site term on a row of {beside}"
        );
    }
}

/// The rows of `--terms display`: both site terms print one name, the fifth base, the donor
/// region, the tract and `splice_region_variant` another, the intron term a third.
#[test]
fn under_the_display_names_the_class_reads_the_names_the_seven_terms_print() {
    assert!(
        after_the_span(
            "T",
            "7:101000",
            "TT",
            "3",
            ("FRAMESHIFT_CODING,SPLICE_SITE", "FRAMESHIFT_CODING")
        ),
        "the exon's third base, one base added after it"
    );
    assert!(
        after_the_span(
            "T",
            "7:101000",
            "TGTG",
            "-",
            ("ESSENTIAL_SPLICE_SITE", "SPLICE_SITE,INTRONIC")
        ),
        "the donor's second base, three bases added after it"
    );
    assert!(
        after_the_span(
            "TA",
            TWO_POSITIONS,
            "GTGTG",
            "999-1000",
            (
                "ESSENTIAL_SPLICE_SITE,FRAMESHIFT_CODING,SPLICE_SITE,INTRONIC",
                "FRAMESHIFT_CODING,SPLICE_SITE"
            )
        ),
        "two replaced bases inside an exon, the donor and intron names gone"
    );
    let refused = [
        (
            "the same names on both rows, the fifth base against the donor region",
            "-",
            "SPLICE_SITE,INTRONIC",
            "SPLICE_SITE,INTRONIC",
        ),
        (
            "another term that differs",
            "3",
            "FRAMESHIFT_CODING,SPLICE_SITE",
            "STOP_GAINED",
        ),
        (
            "names no place of the intron gives",
            "-",
            "INTRONIC",
            "ESSENTIAL_SPLICE_SITE",
        ),
        (
            "a Sequence Ontology name beside a display name",
            "-",
            "ESSENTIAL_SPLICE_SITE",
            "splice_region_variant,INTRONIC",
        ),
        (
            "the name of the donor and the acceptor term on the candidate's row alone",
            "-",
            "WITHIN_NON_CODING_GENE",
            "ESSENTIAL_SPLICE_SITE,WITHIN_NON_CODING_GENE",
        ),
    ];
    for (shape, cdna, reference, candidate) in refused {
        assert!(
            !after_the_span("T", "7:101000", "TGTG", cdna, (reference, candidate)),
            "{shape}"
        );
    }
    assert!(
        !after_the_span(
            "TAC",
            "7:100998-101000",
            "GAT",
            "-",
            ("INTRONIC", "SPLICE_SITE,INTRONIC")
        ),
        "three replaced bases under an allele as long, its first and last letter differing"
    );
    assert!(
        !after_the_span(
            "TAC",
            "7:100998-101000",
            "GACTG",
            "-",
            ("INTRONIC", "SPLICE_SITE,INTRONIC")
        ),
        "a longer allele that does not open with its reference, no name gone from the reference's row"
    );
    let unread = |terms: &str| row("7:101000", "TT", terms, "3", "-");
    assert_eq!(
        splice_term_past_the_record_span_pair(
            &unread("FRAMESHIFT_CODING,SPLICE_SITE"),
            &unread("FRAMESHIFT_CODING")
        ),
        Ok(false),
        "without the record"
    );
}

#[test]
fn an_allele_that_does_not_open_with_its_whole_reference_is_no_insertion_after_the_span() {
    for (place, reference, candidate) in AFTER_THE_SPAN_PLACES {
        let terms = (reference, candidate);
        assert!(
            !after_the_span("TA", TWO_POSITIONS, "TGGTG", "-", terms),
            "{place}, the second base differs"
        );
        assert!(
            !after_the_span("TA", TWO_POSITIONS, "GAGTG", "-", terms),
            "{place}, the first base differs"
        );
        let unread = |terms: &str| row(TWO_POSITIONS, "TAGTG", terms, "-", "-");
        assert_eq!(
            splice_term_past_the_record_span_pair(&unread(reference), &unread(candidate)),
            Ok(false),
            "{place}, without the record"
        );
    }
}

#[test]
fn a_term_an_insertion_after_the_span_does_not_explain_stays_a_difference() {
    let refused = [
        (
            "a donor term the candidate alone holds",
            "intron_variant",
            "splice_polypyrimidine_tract_variant,intron_variant,splice_donor_variant",
        ),
        (
            "splice_region_variant the reference alone holds",
            "splice_region_variant,intron_variant",
            "splice_polypyrimidine_tract_variant,intron_variant",
        ),
        (
            "another term the candidate alone holds",
            "splice_acceptor_variant",
            "splice_region_variant,splice_polypyrimidine_tract_variant,intron_variant,missense_variant",
        ),
        (
            "another term the reference alone holds",
            "splice_donor_variant,non_coding_transcript_variant",
            "splice_region_variant,intron_variant",
        ),
        (
            "the intron term the insertion's place gives, missing beside the acceptor bases",
            "splice_acceptor_variant",
            "splice_region_variant,splice_polypyrimidine_tract_variant",
        ),
        (
            "the intron term the insertion's place gives, missing beside the donor bases",
            "splice_donor_variant",
            "splice_region_variant",
        ),
        (
            "the same terms on both rows",
            "intron_variant",
            "intron_variant",
        ),
    ];
    for (shape, reference, candidate) in refused {
        assert!(
            !after_the_span("TA", TWO_POSITIONS, "TAGTG", "-", (reference, candidate)),
            "{shape}"
        );
    }
}
