// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Locations, positions and mate names as the one-sided classes read them.

use super::*;

fn word(value: i128) -> Coordinate {
    Coordinate::Word(value)
}

#[test]
fn a_symbolic_type_without_a_term_is_in_the_batch_scope_under_either_spelling() {
    for allele in ["CPX", "CTX", "<CPX>", "<CTX>"] {
        assert!(
            batch_dependent_annotation_scope("21:800-900", allele),
            "{allele}"
        );
    }
    assert!(!batch_dependent_annotation_scope("21:800-900", "deletion"));
}

#[test]
fn a_location_with_a_negative_coordinate_is_read_and_is_an_integer() {
    assert_eq!(
        location_span("21:-4-295"),
        Some(("21", word(-4), word(295)))
    );
    assert_eq!(
        location_span("21:-2147483647--2147483348"),
        Some(("21", word(-2_147_483_647), word(-2_147_483_348)))
    );
    assert_eq!(location_span("21:-5"), Some(("21", word(-5), word(-5))));
    let past_64_bits = (1 << 64) + 1;
    assert_eq!(
        location_span("21:18446744073709551617"),
        Some(("21", word(past_64_bits), word(past_64_bits)))
    );
    assert_eq!(location_span("21:30425651-NaN"), None);
    for location in [
        "21:-4-295",
        "21:-2147483647--2147483348",
        "21:-5",
        "21:30425651-30425949",
        "21:30425651",
    ] {
        assert!(!non_integer_end(location), "{location}");
    }
    for location in [
        "21:30425651-NaN",
        "21:30425949.0",
        "21: 30425949",
        "21:30425949 ",
        "21:3e7",
        "21:30425651-30425950.7",
    ] {
        assert!(non_integer_end(location), "{location}");
    }
}

#[test]
fn a_mate_chromosome_is_read_as_the_reference_reads_it() {
    for (allele, name) in [
        ("N[chr21:30425700[", "21"),
        ("N[Chr21:30425700[", "Chr21"),
        ("N[CHR21:30425700[", "CHR21"),
        ("N[21:30425700[", "21"),
        ("]NC_000022.11:26037668]N", "NC_000022.11"),
        ("N[chrM:100[", "M"),
    ] {
        assert_eq!(reference_mate_chromosome(allele), Some(name), "{allele}");
    }
    assert_eq!(reference_mate_chromosome("N."), None);
    assert_eq!(reference_mate_chromosome("deletion"), None);
}

#[test]
fn a_pos_is_read_in_every_spelling_of_an_exact_integer_the_scored_engine_accepts() {
    for text in [
        "30425650",
        "3.042565e7",
        "30425650.0",
        "+30425650",
        "30425650.",
        " 30425650 ",
        "030425650",
    ] {
        assert_eq!(vep_rs_pos(text), Some(30_425_650), "{text:?}");
    }
    assert_eq!(vep_rs_pos("3e7"), Some(30_000_000));
    assert_eq!(VEP_RS_COORDINATE_CEILING, (1 << 62) - 1);
    assert_eq!(
        vep_rs_pos(&VEP_RS_COORDINATE_CEILING.to_string()),
        Some(VEP_RS_COORDINATE_CEILING)
    );
    let past_a_double = (1i128 << 56) + 1;
    assert_eq!(
        vep_rs_pos(&past_a_double.to_string()),
        Some(past_a_double),
        "an unpadded digit string is read exactly"
    );
    assert_eq!(
        vep_rs_pos(&format!(" {past_a_double} ")),
        Some(1 << 56),
        "a padded one is read as the double it names"
    );
    let past_the_ceiling = (VEP_RS_COORDINATE_CEILING + 1).to_string();
    let padded_ceiling = format!(" {VEP_RS_COORDINATE_CEILING} ");
    for text in [
        ".",
        "",
        "abc",
        "-30425650",
        "30_425_650",
        "３０４２５６５０",
        "0x1D04A1A",
        "inf",
        "nan",
        "1e400",
        "30425650.5",
        "++30425650",
        &past_the_ceiling,
        &padded_ceiling,
    ] {
        assert_eq!(vep_rs_pos(text), None, "{text:?}");
    }
}

#[test]
fn the_points_of_a_dot_breakend_are_its_spans_start_and_the_base_after_it() {
    let points = |start: i128, end: i128| dot_breakend_points(&("21", word(start), word(end)));
    assert_eq!(
        points(30_425_650, 30_425_651),
        ["21:30425650", "21:30425651"]
    );
    assert_eq!(
        points(30_425_651, 30_425_653),
        ["21:30425651", "21:30425652"]
    );
    assert_eq!(
        points(30_425_650, 30_425_653),
        ["21:30425650", "21:30425651"]
    );
}
