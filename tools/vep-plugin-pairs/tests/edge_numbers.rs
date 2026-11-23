// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Whole numbers at the ends of `i128` in a row's Location, its HGVSc, its `HGVS_OFFSET`
//! and a cache's exons: each rule gives the answer that arithmetic without a bound gives.

use std::collections::BTreeSet;
use std::fs;
use vep_plugin_pairs::rules::downstream::classes::{codon_span, hgvs_rotated_allele};
use vep_plugin_pairs::rules::downstream::readings::shifted_pair_readings;
use vep_plugin_pairs::rules::downstream::transcript::TranscriptModel;
use vep_plugin_pairs::rules::entry::{parse_default, replaced_reference, Row};
use vep_plugin_pairs::rules::maxentscan::SpliceWindows;
use vep_plugin_pairs::rules::shared::{
    indel_row, insertion_row, location_span, notation_indel_length, record_allele_lengths,
};
use vep_plugin_pairs::rules::{BuiltIn, Consult, Rules};

const MAX: i128 = i128::MAX;
const MIN: i128 = i128::MIN;
const HEADER: &str = "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\n";

/// A transcript row under an identifier, which states no allele of its record.
fn entry(location: &str, allele: &str, feature: &str, consequence: &str, extra: &str) -> Row {
    let line = format!(
        "rs5\t{location}\t{allele}\tG1\t{feature}\tTranscript\t{consequence}\t-\t-\t-\t-\t-\t-\t{extra}\n"
    );
    parse_default(&format!("{HEADER}{line}")).unwrap().remove(0)
}

fn at(location: &str, allele: &str) -> Row {
    entry(location, allele, "ENST1", "missense_variant", "-")
}

#[test]
fn reverse_strand_rotation_counts_an_offset_of_the_smallest_number() {
    let accepts = |offset: i128, notation: &str| {
        let extra = format!("STRAND=-1;HGVS_OFFSET={offset};HGVSc=ENST1.2:{notation}");
        let e = entry("7:100-101", "A", "ENST1", "frameshift_variant", &extra);
        let declared = BTreeSet::new();
        BuiltIn
            .accepts(&Consult {
                klass: "downstream_reverse_strand_rotation",
                e: &e,
                name: "DownstreamProtein",
                expected: "MK",
                actual: "ML",
                a: Some(&e),
                record_rows: None,
                candidate_rows: None,
                declared: &declared,
                assembly: None,
                source: None,
                gwas_pairs: None,
                revel_table: None,
                utrannotator_sweep: None,
                downstream_sequences: None,
                splice_windows: None,
                spliceai_table: None,
                exon_spans: None,
            })
            .unwrap()
    };
    let to_the_largest = |first: i128| format!("c.{first}_{MAX}del");
    for (offset, notation, expected) in [
        (MIN, "c.5_7del".to_string(), true),
        (MIN, "c.5del".to_string(), false),
        (MIN, "c.5_6insACGT".to_string(), false),
        (MIN, to_the_largest(1), true),
        (MIN, to_the_largest(2), true),
        (MIN, to_the_largest(-1), false),
        (MIN, to_the_largest(MIN), false),
        (MIN, format!("c.{MAX}_1del"), false),
        (MAX, "c.5_7del".to_string(), true),
        (MAX, "c.5_6insACGT".to_string(), true),
        (MAX, to_the_largest(1), false),
        (MAX, to_the_largest(2), true),
        (7, "c.5_7del".to_string(), true),
        (7, to_the_largest(2), false),
    ] {
        assert_eq!(accepts(offset, &notation), expected, "{offset} {notation}");
    }
}

#[test]
fn an_allele_is_turned_by_a_shift_of_any_size() {
    for (allele, shift, strand, expected) in [
        ("ACG", MIN, 1, "GAC"),
        ("ACG", MIN, -1, "CGT"),
        ("ACG", MAX, 1, "CGA"),
        ("ACG", MAX, -1, "CGT"),
        ("ACG", 4, 1, "CGA"),
        ("ACG", 5, 1, "GAC"),
        ("ACG", 4, -1, "CGT"),
        ("ACGT", MIN, 1, "ACGT"),
        ("ACGT", MAX, 1, "TACG"),
        ("ACGT", MAX - 1, 1, "GTAC"),
        ("A", MIN, -1, "T"),
        ("", MIN, 1, ""),
    ] {
        assert_eq!(
            hgvs_rotated_allele(allele, shift, strand),
            expected,
            "{allele} {shift} {strand}"
        );
    }
}

#[test]
fn a_location_that_ends_at_the_largest_number_is_read_by_its_span() {
    let whole = format!("7:0-{MAX}");
    for allele in ["A", "-", "ACGT"] {
        assert_eq!(indel_row(&at(&whole, allele)), Ok(true), "{allele}");
        assert!(!insertion_row(&at(&whole, allele)), "{allele}");
        // The span is one more than the largest `i128` and is read as that one.
        assert_eq!(location_span(&at(&whole, allele)), Ok(Some(MAX)));
    }
    assert_eq!(
        record_allele_lengths(&at(&whole, "ACGT"), None),
        Some((MAX, 4))
    );
    let from_one = format!("7:1-{MAX}");
    assert_eq!(location_span(&at(&from_one, "A")), Ok(Some(MAX)));
    assert_eq!(
        record_allele_lengths(&at(&from_one, "-"), None),
        Some((MAX, 0))
    );
    let one_base = format!("7:{MAX}-{MAX}");
    assert_eq!(indel_row(&at(&one_base, "A")), Ok(false));
    assert!(!insertion_row(&at(&one_base, "A")));
    assert_eq!(location_span(&at(&one_base, "A")), Ok(Some(1)));
    assert_eq!(
        record_allele_lengths(&at(&one_base, "A"), None),
        Some((1, 1))
    );
    let last_two = format!("7:{}-{MAX}", MAX - 1);
    assert_eq!(indel_row(&at(&last_two, "A")), Ok(true));
    assert!(insertion_row(&at(&last_two, "A")));
    assert!(!insertion_row(&at(&last_two, "-")));
    assert_eq!(location_span(&at(&last_two, "A")), Ok(Some(2)));
    assert_eq!(record_allele_lengths(&at(&last_two, "A"), None), None);
    assert_eq!(
        record_allele_lengths(&at(&last_two, "-"), None),
        Some((2, 0))
    );
    let turned = format!("7:{MAX}-{}", MAX - 1);
    assert_eq!(indel_row(&at(&turned, "A")), Ok(true));
    assert!(!insertion_row(&at(&turned, "A")));
    assert_eq!(location_span(&at(&turned, "A")), Ok(None));
    assert_eq!(record_allele_lengths(&at(&turned, "A"), None), Some((0, 1)));
    let down_to_zero = format!("7:{MAX}-0");
    assert!(!insertion_row(&at(&down_to_zero, "A")));
    assert_eq!(location_span(&at(&down_to_zero, "A")), Ok(None));
    assert_eq!(
        record_allele_lengths(&at(&down_to_zero, "A"), None),
        Some((1 - MAX, 1))
    );
}

#[test]
fn the_bases_an_allele_stands_in_for_are_read_at_the_largest_position() {
    let read = |position: i128, alternate: &str, location: &str, allele: &str| {
        replaced_reference(position, "ACG", &[alternate], Some(location), allele)
    };
    assert_eq!(
        read(MAX - 1, "ATG", &format!("7:{MAX}"), "T"),
        Ok(Some("C".to_string()))
    );
    assert_eq!(read(0, "ATG", &format!("7:{MAX}-{MAX}"), "T"), Ok(None));
    assert_eq!(read(0, "ATCG", &format!("7:{MAX}-5"), "T"), Ok(None));
    assert_eq!(read(0, "AC", &format!("7:{MAX}"), "-"), Ok(None));
    assert_eq!(read(MAX, "ATG", &format!("7:{MAX}"), "A"), Ok(None));
    assert_eq!(read(5, "ACTG", "7:6-7", "T"), Ok(Some(String::new())));
    assert_eq!(read(5, "ACTG", "7:7-6", "T"), Ok(None));
}

#[test]
fn the_length_an_hgvsc_states_is_read_up_to_the_largest_number() {
    let length = |first: i128, last: i128, kind: &str| {
        notation_indel_length(&format!("ENST1.2:c.{first}_{last}{kind}"))
    };
    assert_eq!(length(5, 7, "del"), Some(3));
    assert_eq!(length(1, MAX, "del"), Some(MAX));
    assert_eq!(length(2, MAX, "del"), Some(MAX - 1));
    assert_eq!(length(1 - MAX, 1, "del"), Some(MAX));
    assert_eq!(length(MAX, 1, "del"), Some(2 - MAX));
    assert_eq!(length(MAX, MAX, "dup"), Some(1));
    assert_eq!(length(MIN, MIN, "dup"), Some(1));
    // A length past the largest `i128` is none.
    for (first, last) in [
        (-1, MAX),
        (0, MAX),
        (MIN, MAX),
        (MIN, -1),
        (MIN, 1),
        (-2, MAX - 1),
    ] {
        assert_eq!(length(first, last, "del"), None, "{first} {last}");
    }
}

#[test]
fn the_codons_of_a_pair_are_read_at_both_ends_of_the_numbers() {
    let third = 56713727820156410577229101238628035243_i128;
    for (pair, expected) in [
        ((MAX, MAX), (third, third)),
        ((MAX - 1, MAX), (third - 1, third)),
        ((MAX - 2, 1), (1, third - 1)),
        ((MIN, MAX), (1 - third, third)),
        ((MIN, MIN), (1 - third, 1 - third)),
        ((-MAX, 2), (1 - third, 1)),
        ((0, 1), (0, 1)),
        ((1, 3), (1, 1)),
        ((4, 6), (2, 2)),
        ((-1, -3), (-1, 0)),
        ((-4, 7), (-1, 3)),
    ] {
        assert_eq!(codon_span(pair), expected, "{pair:?}");
    }
}

#[test]
fn splice_windows_are_placed_beside_exons_at_the_ends_of_the_numbers() {
    let root = std::env::temp_dir().join(format!(
        "vep-plugin-pairs-edge-windows-{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&root);
    let directory = root.join("transcripts").join("7");
    fs::create_dir_all(&directory).unwrap();
    let transcript = |stable_id: &str, strand: i8, exons: [(i128, i128); 3]| {
        let exons: Vec<String> = exons
            .iter()
            .map(|(start, end)| format!("{{\"start\": {start}, \"end\": {end}}}"))
            .collect();
        format!(
            "{{\"stable_id\": \"{stable_id}\", \"strand\": {strand}, \"exons\": [{}]}}",
            exons.join(", ")
        )
    };
    let inside = [(100, 200), (300, 301), (303, 400)];
    let at_the_largest = [(1, 5), (MAX - 1, MAX), (MAX, MAX)];
    let up_to_the_largest = [(1, MAX - 4), (MAX - 2, MAX - 2), (MAX, MAX)];
    let past_the_largest = [(3, MAX), (4, 10), (20, 30)];
    let at_the_smallest = [(MIN, MIN), (MIN, -5), (0, 10)];
    let from_the_smallest = [(MIN, 2), (4, 5), (7, 20)];
    let mut transcripts = Vec::new();
    for (name, exons) in [
        ("INSIDE", inside),
        ("LARGEST", at_the_largest),
        ("UPTO", up_to_the_largest),
        ("PAST", past_the_largest),
        ("SMALLEST", at_the_smallest),
        ("FROM", from_the_smallest),
    ] {
        transcripts.push(transcript(&format!("ENSTF{name}"), 1, exons));
        transcripts.push(transcript(&format!("ENSTR{name}"), -1, exons));
    }
    fs::write(
        directory.join("0-part.json"),
        format!("[{}]", transcripts.join(", ")),
    )
    .unwrap();
    let windows = SpliceWindows::new(&root).unwrap();
    let low: Vec<i128> = (0..25).chain(295..312).collect();
    let high: Vec<i128> = (0..28).map(|below| MAX - below).collect();
    for (name, held) in [
        ("INSIDE", vec![299, 300, 301, 302]),
        ("LARGEST", vec![MAX - 2, MAX - 1, MAX]),
        ("UPTO", vec![MAX - 4, MAX - 3, MAX - 2, MAX - 1, MAX]),
        ("PAST", vec![]),
        ("SMALLEST", vec![]),
        ("FROM", vec![3, 4, 5, 6]),
    ] {
        for strand in ["F", "R"] {
            let feature = format!("ENST{strand}{name}");
            let mut found: Vec<i128> = low
                .iter()
                .chain(&high)
                .copied()
                .filter(|position| {
                    let e = entry(
                        &format!("7:{position}"),
                        "A",
                        &feature,
                        "intron_variant",
                        "-",
                    );
                    windows.holds(&e).unwrap()
                })
                .collect();
            found.sort_unstable();
            assert_eq!(found, held, "{feature}");
        }
    }
    let _ = fs::remove_dir_all(root);
}

#[test]
fn a_transcript_model_maps_coordinates_at_the_ends_of_the_numbers() {
    let model = |strand: i8, pair: [i128; 4], coding: (i128, i128)| TranscriptModel {
        strand,
        exons: vec![(pair[2], pair[3])],
        coding_start: coding.0,
        coding_end: coding.1,
        padding: 0,
        pairs: vec![[pair[0], pair[1], pair[2], pair[3], i128::from(strand)]],
    };
    // cDNA coordinates up to the largest number.
    let top = model(1, [MAX - 17, MAX, 100, 117], (MAX - 17, MAX));
    assert_eq!(top.cds_pair(117, 117), (Some(18), Some(18)));
    assert_eq!(top.cds_pair(116, 117), (Some(17), Some(18)));
    assert_eq!(top.cds_pair(118, 117), (None, None));
    assert_eq!(top.cds_pair(117, 116), (Some(18), Some(17)));
    assert_eq!(top.genomic_span(17, 18), Some((116, 117)));
    assert_eq!(top.genomic_span(18, 17), Some((117, 116)));
    assert_eq!(top.genomic_span(19, 19), None);
    // A pair whose cDNA coordinates leave the numbers after its sixth base.
    let short = model(1, [MAX - 5, MAX, 100, 117], (MAX - 5, MAX));
    assert_eq!(short.cds_pair(105, 105), (Some(6), Some(6)));
    assert_eq!(short.cds_pair(106, 106), (None, None));
    assert_eq!(short.cds_pair(106, 105), (None, None));
    assert_eq!(short.cds_pair(105, 104), (Some(6), Some(5)));
    assert_eq!(short.genomic_span(6, 5), Some((105, 104)));
    assert_eq!(short.genomic_span(7, 7), None);
    let reverse = model(-1, [MAX - 5, MAX, 100, 117], (MAX - 5, MAX));
    assert_eq!(reverse.cds_pair(112, 112), (Some(6), Some(6)));
    assert_eq!(reverse.cds_pair(111, 111), (None, None));
    assert_eq!(reverse.cds_pair(112, 111), (None, None));
    assert_eq!(reverse.cds_pair(113, 112), (Some(6), Some(5)));
    assert_eq!(reverse.genomic_span(6, 5), Some((113, 112)));
    // Genomic bases up to the largest number, on either strand.
    let far = model(1, [1, 18, MAX - 17, MAX], (1, 18));
    assert_eq!(far.cds_pair(MAX, MAX), (Some(18), Some(18)));
    assert_eq!(far.cds_pair(MAX - 1, MAX), (Some(17), Some(18)));
    assert_eq!(far.cds_pair(MAX, MAX - 1), (Some(18), Some(17)));
    assert_eq!(far.genomic_span(18, 18), Some((MAX, MAX)));
    assert_eq!(far.genomic_span(18, 17), Some((MAX, MAX - 1)));
    assert_eq!(far.genomic_span(19, 18), None);
    let far_reverse = model(-1, [1, 18, MAX - 17, MAX], (1, 18));
    assert_eq!(far_reverse.cds_pair(MAX, MAX), (Some(1), Some(1)));
    assert_eq!(far_reverse.cds_pair(MAX, MAX - 1), (Some(2), Some(1)));
    assert_eq!(far_reverse.genomic_span(1, 1), Some((MAX, MAX)));
    assert_eq!(far_reverse.genomic_span(2, 1), Some((MAX, MAX - 1)));
    // Coding bounds at the smallest number hold no base of the pair.
    let low = model(1, [1, 18, 100, 117], (MIN, MIN));
    assert_eq!(low.cds_pair(100, 100), (None, None));
    assert_eq!(low.cds_pair(101, 100), (None, None));
    // The last 51 bases of the second to last exon, beside exons at either end.
    let mut forward = model(1, [1, 50, 1, 50], (1, 61));
    forward.exons = vec![(1, 50), (MAX - 10, MAX)];
    for (genomic_end, escapes) in [
        (MAX - 10, true),
        (50, true),
        (1, true),
        (0, true),
        (MAX - 11, false),
    ] {
        assert_eq!(
            forward.ends_where_nmd_escapes(genomic_end),
            escapes,
            "{genomic_end}"
        );
    }
    let mut backward = model(-1, [1, 101, 100, 200], (1, 112));
    backward.exons = vec![(100, 200), (MIN, MIN + 10)];
    for (genomic_end, escapes) in [
        (MIN + 10, true),
        (100, true),
        (151, true),
        (152, false),
        (99, false),
    ] {
        assert_eq!(
            backward.ends_where_nmd_escapes(genomic_end),
            escapes,
            "{genomic_end}"
        );
    }
}

#[test]
fn an_insertion_moved_to_the_largest_position_is_read_at_its_point() {
    let model = |strand: i8, coding_start: i128| TranscriptModel {
        strand,
        exons: vec![(MAX - 17, MAX)],
        coding_start,
        coding_end: 42,
        padding: 0,
        pairs: vec![[1, 18, MAX - 17, MAX, i128::from(strand)]],
    };
    // The point between the two bases, moved four bases on, lies after the largest position.
    let e = at(&format!("7:{}-{}", MAX - 4, MAX - 3), "A");
    assert_eq!(
        shifted_pair_readings(&e, &model(1, 1), 4, &[(Some(15), Some(14))]),
        [((Some(15), Some(14)), Some((19, 18)), MAX - 4)]
    );
    assert_eq!(
        shifted_pair_readings(&e, &model(-1, -5), 4, &[(Some(11), Some(10))]),
        [((Some(11), Some(10)), Some((7, 6)), MAX - 4)]
    );
    // A span of bases moved past the largest position has no coding base at that end.
    assert_eq!(
        shifted_pair_readings(&e, &model(1, 1), 4, &[(Some(14), Some(15))]),
        [((Some(14), Some(15)), None, MAX - 3)]
    );
    assert_eq!(
        shifted_pair_readings(&e, &model(1, 1), 5, &[(Some(15), Some(14))]),
        [((Some(15), Some(14)), None, MAX - 4)]
    );
}

#[test]
fn a_gained_start_window_is_not_read_over_a_location_of_more_positions_than_the_numbers() {
    let gained = "5UTR_consequence=5_prime_UTR_premature_start_codon_gain_variant;\
                  5UTR_annotation=CapDistanceToStart=4:type=uORF";
    let row = |location: &str, cdna: &str, utr: &str| {
        let line = format!(
            "rs5\t{location}\tACGT\tG1\tENST1\tTranscript\t5_prime_UTR_variant\t{cdna}\t-\t-\t-\t-\t-\tSTRAND=1;{utr}\n"
        );
        parse_default(&format!("{HEADER}{line}")).unwrap().remove(0)
    };
    let accepts = |location: &str, cdna: &str| {
        let e = row(location, cdna, gained);
        let a = row(location, cdna, "5UTR_consequence=-;5UTR_annotation=-");
        let declared = BTreeSet::new();
        BuiltIn
            .accepts(&Consult {
                klass: "utrannotator_gained_start_window",
                e: &e,
                name: "5UTR_consequence",
                expected: e.fields.get_or("5UTR_consequence", "-"),
                actual: a.fields.get_or("5UTR_consequence", "-"),
                a: Some(&a),
                record_rows: None,
                candidate_rows: None,
                declared: &declared,
                assembly: None,
                source: None,
                gwas_pairs: None,
                revel_table: None,
                utrannotator_sweep: None,
                downstream_sequences: None,
                splice_windows: None,
                spliceai_table: None,
                exon_spans: None,
            })
            .unwrap()
    };
    let width = format!("1-{MAX}");
    // The Location names one position more than the cDNA positions do.
    assert!(!accepts(&format!("7:0-{MAX}"), &width));
    assert!(!accepts(&format!("7:0-{MAX}"), &format!("0-{}", MAX - 1)));
    assert!(accepts(&format!("7:1-{MAX}"), &width));
}
