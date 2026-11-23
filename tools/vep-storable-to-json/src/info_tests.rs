// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! Expected values are perl 5.34.0's for the patterns and the sort the comments quote.

use super::{cache_version, contig_order, source_type, Info};
use crate::json::Pretty;

// if ($cache_dir =~ m{/(\d+)(?:_[^/]+)?/?$}) { $version = $1 + 0 }
// if ($cache_dir =~ m{(?:^|/)([^/]+)/\d+(?:_[^/]+)?/?$}) { ... /_refseq$/ ... /_merged$/ ... }
#[test]
fn the_version_and_the_transcript_set_are_read_off_the_path_as_the_patterns_read_them() {
    let integer = |value: i128| Some(Pretty::Integer(value));
    for (path, version, set) in [
        (
            &b"homo_sapiens/116_GRCh38"[..],
            integer(116),
            Some(&b"ensembl"[..]),
        ),
        (
            b"/data/homo_sapiens_refseq/116_GRCh38/",
            integer(116),
            Some(b"refseq"),
        ),
        (
            b"x/homo_sapiens_merged/84_GRCh38",
            integer(84),
            Some(b"merged"),
        ),
        (b"116_GRCh38", None, None),
        (b"/116_GRCh38", integer(116), None),
        (b"//116", integer(116), None),
        (b"a/116", integer(116), Some(b"ensembl")),
        (b"a/116/", integer(116), Some(b"ensembl")),
        (b"a/116//", None, None),
        (b"a/116_", None, None),
        (b"a/116_x/y", None, None),
        (b"a/12abc", None, None),
        (b"a/b_refseq/7_x_y", integer(7), Some(b"refseq")),
        (b"plain", None, None),
        (
            b"sets/homo_sapiens/84_GRCh38",
            integer(84),
            Some(b"ensembl"),
        ),
        (b"a_merged/116_GRCh38", integer(116), Some(b"merged")),
        (b"a_refseq_x/116_GRCh38", integer(116), Some(b"ensembl")),
        (b"a/116\n", integer(116), Some(b"ensembl")),
        (b"a_refseq\n/116", integer(116), Some(b"refseq")),
        (b"a/0116", integer(116), Some(b"ensembl")),
        (
            b"a/99999999999999999999",
            Some(Pretty::Double(1e20)),
            Some(b"ensembl"),
        ),
        (b"input/116_GRCh38", integer(116), Some(b"ensembl")),
        (b"./116", integer(116), Some(b"ensembl")),
        (b"a b/1_c d", integer(1), Some(b"ensembl")),
        (b"_refseq/1", integer(1), Some(b"refseq")),
    ] {
        let shown = String::from_utf8_lossy(path);
        assert_eq!(cache_version(path), version, "cache_version of {shown:?}");
        assert_eq!(source_type(path), set, "source_type of {shown:?}");
    }
}

// sort by_contig_order over the sorted directory names 01 1 10 2 21 GL000008.2
// HSCHR1_CTG1 MT X Y chr1 empty_dir
#[test]
fn contigs_order_numbers_first_then_x_y_mt_then_names() {
    let sorted_names: Vec<Vec<u8>> = [
        "01",
        "1",
        "10",
        "2",
        "21",
        "GL000008.2",
        "HSCHR1_CTG1",
        "MT",
        "X",
        "Y",
        "chr1",
        "empty_dir",
    ]
    .iter()
    .map(|name| name.as_bytes().to_vec())
    .collect();
    let ordered: Vec<String> = contig_order(&sorted_names)
        .into_iter()
        .map(|name| String::from_utf8(name).expect("ASCII"))
        .collect();
    // `01` and `1` are the same number and keep the order they came in.
    assert_eq!(
        ordered,
        [
            "01",
            "1",
            "2",
            "10",
            "21",
            "X",
            "Y",
            "MT",
            "GL000008.2",
            "HSCHR1_CTG1",
            "chr1",
            "empty_dir"
        ]
    );
}

#[test]
fn info_lines_are_taken_as_ensembl_vep_takes_them() {
    let info = Info::from_lines(
        b"# comment\tx\nno tab\nspecies\tfirst\nspecies\tsecond\nassembly\t-\nsift\tb\npolyphen\t\n\
source_sift\t0\nsource_polyphen\t2.2.3\nvariation_cols\ta,,b,\nregulatory\t1\ncell_types\tA,B\n\
tabbed\tone\ttwo\nlast\tno newline",
    );
    assert_eq!(info.string(b"species"), Some(&b"second"[..]));
    assert_eq!(info.string(b"assembly"), None);
    assert_eq!(info.string(b"polyphen"), Some(&b""[..]));
    assert_eq!(info.string(b"tabbed"), Some(&b"one\ttwo"[..]));
    assert_eq!(info.string(b"last"), Some(&b"no newline"[..]));
    assert!(!info.top.contains_key(&b"regulatory"[..]));
    assert!(!info.top.contains_key(&b"cell_types"[..]));
    assert_eq!(
        info.top.get(&b"variation_cols"[..]),
        Some(&Pretty::Array(vec![
            Pretty::String(b"a".to_vec()),
            Pretty::String(b"b".to_vec())
        ]))
    );
    // A version that is the string 0 is a false value in Perl.
    assert!(!info.has_source_version(b"sift"));
    assert!(info.has_source_version(b"polyphen"));
    assert!(!info.has_source_version(b"absent"));
}
