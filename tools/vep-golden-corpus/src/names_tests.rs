// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Unit tests of `names`: the names Ensembl VEP writes for a record without an ID, the names
//! a record line can carry, the multi-allelic strata and the index from a name to records.

use super::*;

fn names(line: &str) -> Vec<String> {
    record_names(line).unwrap()
}

#[test]
fn auto_names_follow_vep_rules() {
    // The first name is release 116's, the raw line as CHR_POS_REF/ALTs; the rest are
    // release 115's minimised names.
    assert_eq!(vep_auto_names("21", 100, "A", &["G"]), ["21_100_A/G"]);
    // A bi-allelic indel: raw alleles, trimmed start.
    assert_eq!(
        vep_auto_names("21", 100, "AT", &["A"]),
        ["21_100_AT/A", "21_101_AT/A"]
    );
    assert_eq!(
        vep_auto_names("21", 100, "A", &["AC"]),
        ["21_100_A/AC", "21_101_A/AC"]
    );
    // A multi-allelic SNV: no chop.
    assert_eq!(
        vep_auto_names("21", 100, "A", &["G", "T"]),
        ["21_100_A/G/T"]
    );
    // A multi-allelic indel: raw alleles, chopped start.
    assert_eq!(
        vep_auto_names("21", 100, "TAG", &["TA", "T"]),
        ["21_100_TAG/TA/T", "21_101_TAG/TA/T"]
    );
    // A repeat: chop, then trim the shared base.
    assert_eq!(
        vep_auto_names("21", 100, "ATT", &["AT"]),
        ["21_100_ATT/AT", "21_102_ATT/AT"]
    );
    // Mixed first bases: no chop.
    assert_eq!(
        vep_auto_names("21", 100, "C", &["CT", "G"]),
        ["21_100_C/CT/G"]
    );
    assert_eq!(
        vep_auto_names("21", 100, "T", &["<CN0>"]),
        ["21_100_T/<CN0>", "21_100_<CN0>", "21_101_<CN0>"]
    );
    // Equal length: raw at POS.
    assert_eq!(vep_auto_names("21", 100, "AT", &["TG"]), ["21_100_AT/TG"]);
}

#[test]
fn auto_names_take_owned_alleles_and_the_last_position_of_64_bits() {
    let alts = vec!["A".to_string()];
    assert_eq!(
        vep_auto_names("21", i64::MAX, "AT", &alts),
        ["21_9223372036854775807_AT/A", "21_9223372036854775808_AT/A"]
    );
    assert_eq!(
        vep_auto_names("21", -1, "N", &["N[21:5["]),
        ["21_-1_N/N[21:5[", "21_-1_N[21:5[", "21_0_N[21:5["]
    );
    // A `*` allele does not hold the chop back.
    assert_eq!(
        vep_auto_names("21", 100, "TA", &["T", "*"]),
        ["21_100_TA/T/*", "21_101_TA/T/*"]
    );
}

#[test]
fn record_names_take_the_first_id_element_perl_reads_as_true() {
    // `rs1;rs2` names its record `rs1`, while `.`, `;rs2`, `0;rs2`, an empty ID column and a
    // multi-allelic `0` leave Ensembl VEP to name the record from its line. An empty or `0`
    // element also offers the name built from the parsed alleles in the case written.
    assert_eq!(names("21\t30289183\trs1;rs2\tC\tG\t.\t.\t."), ["rs1"]);
    assert_eq!(
        names("21\t30289183\t.\tC\tG,T\t.\t.\t."),
        ["21_30289183_C/G/T"]
    );
    assert_eq!(
        names("21\t30289183\t;rs2\tC\tG\t.\t.\t."),
        ["21_30289183_C/G"]
    );
    assert_eq!(
        names("21\t30289183\t0;rs2\tC\tT\t.\t.\t."),
        ["21_30289183_C/T"]
    );
    assert_eq!(
        names("21\t30396602\t\tTTGA\tT\t.\t.\t."),
        [
            "21_30396602_TTGA/T",
            "21_30396603_TTGA/T",
            "21_30396603_TGA/-"
        ]
    );
    assert_eq!(
        names("21\t30396602\t0\tTTGA\tT,TTGATC\t.\t.\t."),
        [
            "21_30396602_TTGA/T/TTGATC",
            "21_30396603_TTGA/T/TTGATC",
            "21_30396603_TGA/-/TGATC",
        ]
    );
    assert_eq!(
        names("21\t30319403\t0\ttca\tt\t.\t.\t."),
        ["21_30319403_tca/t", "21_30319404_tca/t", "21_30319404_ca/-"]
    );
    // A shared trailing base is trimmed from the alleles but moves no start; a `.` element
    // never takes the parsed name, and a symbolic allele has none.
    assert_eq!(
        names("21\t100\t\tattg\tatg\t.\t.\t."),
        ["21_100_attg/atg", "21_102_attg/atg", "21_102_t/-"]
    );
    assert_eq!(
        names("21\t30319403\t.\ttca\tt\t.\t.\t."),
        ["21_30319403_tca/t", "21_30319404_tca/t"]
    );
    assert_eq!(
        names("21\t100\t0\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=200"),
        vep_auto_names("21", 100, "N", &["<DEL>"])
    );
}

#[test]
fn record_names_read_the_position_as_int_reads_it_and_five_columns_alone() {
    assert_eq!(names("21\t 0100 \t.\tA\tG"), ["21_100_A/G"]);
    assert_eq!(names("21\t+1_000\t\tA\tG"), ["21_1000_A/G"]);
    assert_eq!(names("21\t100\trs9\tA\tG\textra"), ["rs9"]);
    // A `*` allele of a multi-allelic indel keeps its one character in the parsed name.
    assert_eq!(
        names("21\t100\t0\tTA\tT,*,TAA\t.\t.\t."),
        ["21_100_TA/T/*/TAA", "21_101_TA/T/*/TAA", "21_101_A/-/*/AA"]
    );
    // An insertion parsed to `-` and bases keeps the anchor chop and trims nothing more.
    assert_eq!(
        names("21\t100\t0\tA\tAA\t.\t.\t."),
        ["21_100_A/AA", "21_101_A/AA", "21_101_-/A"]
    );
    for (line, detail) in [
        (
            "21\tx\trs1\tA\tG",
            "the record 21 x rs1 has a POS that is no integer of 64 bits",
        ),
        (
            "21\t9223372036854775808\t.\tA\tG",
            "the record 21 9223372036854775808 . has a POS that is no integer of 64 bits",
        ),
        (
            "21\t100\trs1\tA",
            "the record line '21\\t100\\trs1\\tA' has fewer than five columns",
        ),
        ("", "the record line '' has fewer than five columns"),
    ] {
        assert_eq!(record_names(line), Err(Failure::own(detail)), "{line:?}");
    }
}

#[test]
fn multi_allelic_strata() {
    let stratum = |ref_allele: &str, alts: &[&str]| multi_allelic_stratum(ref_allele, alts);
    assert_eq!(stratum("A", &["G", "T"]), Ok("multi_allelic_snv"));
    // Chopped to A/-/GG: nothing left to trim.
    assert_eq!(stratum("CA", &["C", "CGG"]), Ok("multi_allelic_indel"));
    // Chopped to A/-/AA: a base is shared.
    assert_eq!(
        stratum("CA", &["C", "CAA"]),
        Ok("multi_allelic_indel_trimmable")
    );
    // REF and the long ALT share their last base once the first is chopped.
    assert_eq!(
        stratum("TCACACA", &["TACACACACACACACACA", "T"]),
        Ok("multi_allelic_indel_trimmable")
    );
    // Chopped to TTT/-/TT: a base is shared.
    assert_eq!(
        stratum("CTTT", &["C", "CTT"]),
        Ok("multi_allelic_indel_trimmable")
    );
    // Alleles that do not share their first base are read as written.
    assert_eq!(stratum("CA", &["G", "TGG"]), Ok("multi_allelic_indel"));
    assert_eq!(
        stratum("CA", &["G", "TGA"]),
        Ok("multi_allelic_indel_trimmable")
    );
    assert_eq!(MULTI_ALLELIC_STRATA[1], "multi_allelic_indel");
    assert_eq!(
        stratum("A", &["", "TT"]),
        Err(Failure::own(
            "a record with REF 'A' and ALT ',TT' has an empty allele"
        ))
    );
    // Empty alleles all share the empty first character and are of one length.
    assert_eq!(stratum("", &["", ""]), Ok("multi_allelic_snv"));
}

#[test]
fn a_name_indexes_every_record_line_that_can_carry_it() {
    let lines = [
        "21\t100\trs1\tA\tG",
        "21\t200\t.\tAT\tA",
        "22\t100\trs1;rs2\tC\tT",
    ];
    let index = record_name_index(lines).unwrap();
    assert_eq!(index["rs1"], [0, 2]);
    assert_eq!(index["21_200_AT/A"], [1]);
    assert_eq!(index["21_201_AT/A"], [1]);
    assert_eq!(index.len(), 3);
    assert_eq!(
        record_name_index(["21\t100\trs1\tA\tG", "21\tone\t.\tA\tG"]),
        Err(Failure::own(
            "the record 21 one . has a POS that is no integer of 64 bits"
        ))
    );
    assert!(record_name_index([]).unwrap().is_empty());
}

#[test]
fn texts_are_cut_at_characters() {
    let text = "\u{e9}a\u{1F600}";
    assert_eq!(length(text), 3);
    assert_eq!(first(text), "\u{e9}");
    assert_eq!(rest(text), "a\u{1F600}");
    assert_eq!(last(text), "\u{1F600}");
    assert_eq!(without_last(text), "\u{e9}a");
    for cut in [first, rest, last, without_last] {
        assert_eq!(cut(""), "");
    }
    assert_eq!(or_dash(""), "-");
    assert_eq!(or_dash("A"), "A");
    assert!(share_first(["AT", "A", "*"]));
    assert!(!share_first(["AT", "T"]));
    assert!(!share_first(["*"]));
    assert!(is_symbolic_or_breakend("<DEL>") && is_symbolic_or_breakend("A]1:5]"));
    assert!(!is_symbolic_or_breakend("A") && !is_symbolic_or_breakend(".A"));
}
