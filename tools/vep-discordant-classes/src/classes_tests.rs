// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The pairing held to hand-worked tables, at every thread bound.

use super::*;
use crate::scan::scan;
use crate::text::Text;

const HEADER: &str = "file_name\tsource\tlocation\tallele\tfeature\tconsequence_set\n";

/// The classes of the rows of a table whose rows are `source location allele feature
/// consequence_set`, by name.
fn classes_of(rows: &[&str], threads: usize) -> Vec<&'static str> {
    let mut table = String::from(HEADER);
    for row in rows {
        table.push_str("o\t");
        table.push_str(&row.split(' ').collect::<Vec<_>>().join("\t"));
        table.push('\n');
    }
    let (header, keys) = match scan(Text::new(table.as_bytes()), threads) {
        Ok(read) => read,
        Err(failure) => panic!("{failure:?}"),
    };
    let names = |classes: Vec<Class>| -> Vec<&'static str> {
        classes
            .into_iter()
            .map(|class| class.name().unwrap_or("as read"))
            .collect()
    };
    let as_read = header.class.is_some();
    let classes = names(classify(&keys, as_read, threads));
    // Whichever half of the table is held whole, the classes are the same.
    for reference_whole in [true, false] {
        assert_eq!(
            names(classify_holding(&keys, as_read, threads, reference_whole)),
            classes,
            "the reference-only rows held whole: {reference_whole}"
        );
    }
    classes
}

fn assert_classes(rows: &[&str], expected: &[&str]) {
    for threads in [1, 2, 3, 8] {
        assert_eq!(classes_of(rows, threads), expected, "{threads} threads");
    }
}

#[test]
fn rows_alike_but_for_the_end_of_their_location_are_coordinate_only() {
    assert_classes(
        &[
            "missing_in_rust 1:100-200 A T1 x",
            "extra_in_rust 1:100-250 A T1 x",
        ],
        &["coordinate_only", "coordinate_only"],
    );
}

#[test]
fn rows_alike_but_for_their_consequence_set_are_a_consequence_swap() {
    assert_classes(
        &["only_perl 2:500 A T2 x", "only_rust 2:500 A T2 y"],
        &["consequence_swap", "consequence_swap"],
    );
}

#[test]
fn a_variant_the_other_engine_wrote_under_other_features_alone_is_transcript_selection() {
    assert_classes(
        &[
            "missing_in_rust 3:700 C T3a x",
            "extra_in_rust 3:700 C T3b x",
            "extra_in_rust 3:700 C T3c x",
        ],
        &["transcript_selection"; 3],
    );
}

#[test]
fn a_variant_the_other_engine_did_not_write_is_missing_or_extra() {
    assert_classes(
        &["missing_in_rust 4:900 G T4 x", "extra_in_rust 5:900 T T5 x"],
        &["true_missing", "true_extra"],
    );
}

#[test]
fn a_feature_both_engines_wrote_with_another_end_and_another_set_is_missing_and_extra() {
    assert_classes(
        &[
            "missing_in_rust 6:10-20 A T6 x",
            "extra_in_rust 6:10-30 A T6 y",
        ],
        &["true_missing", "true_extra"],
    );
}

#[test]
fn a_scored_row_at_the_reference_row_s_own_location_is_passed_for_good() {
    // The first reference row passes the scored row at its own Location and takes the
    // next. The second starts behind both, so it never meets the row that was passed,
    // whose Location differs from its own; the second pairing takes it at its Location.
    assert_classes(
        &[
            "missing_in_rust 7:1-5 A T7 x",
            "missing_in_rust 7:1-6 A T7 x",
            "missing_in_rust 7:1-7 A T7 x",
            "extra_in_rust 7:1-5 A T7 x",
            "extra_in_rust 7:1-8 A T7 x",
            "extra_in_rust 7:1-6 A T7 x",
        ],
        &[
            "coordinate_only",
            "consequence_swap",
            "true_missing",
            "true_extra",
            "coordinate_only",
            "consequence_swap",
        ],
    );
}

#[test]
fn the_second_pairing_takes_scored_rows_in_file_order_one_a_reference_row() {
    assert_classes(
        &[
            "extra_in_rust 8:40 A T8 d",
            "missing_in_rust 8:40 A T8 a",
            "missing_in_rust 8:40 A T8 b",
            "extra_in_rust 8:40 A T8 e",
            "missing_in_rust 8:40 A T8 c",
        ],
        &[
            "consequence_swap",
            "consequence_swap",
            "consequence_swap",
            "consequence_swap",
            "true_missing",
        ],
    );
}

#[test]
fn any_source_but_the_two_reference_only_values_is_the_scored_side() {
    assert_classes(
        &[
            "other 9:1-2 A T9 x",
            "only_perl 9:1-3 A T9 x",
            " 9:5 A T9 x",
        ],
        &["coordinate_only", "coordinate_only", "true_extra"],
    );
}

#[test]
fn a_table_of_many_rows_is_classed_alike_at_any_bound() {
    // Blocks of fifty variants, the reference-only rows of a block before its scored-only
    // rows; the shape of a variant follows its number.
    let mut rows: Vec<String> = Vec::new();
    for block in 0..80 {
        let variants = block * 50..(block + 1) * 50;
        for n in variants.clone() {
            match n % 5 {
                0 => rows.push(format!("missing_in_rust 1:{n}-{} A T x", n + 1)),
                1 => rows.push(format!("missing_in_rust 2:{n} A T x")),
                2 => rows.push(format!("missing_in_rust 3:{n} A Ta x")),
                3 => rows.push(format!("missing_in_rust 4:{n} A T x")),
                _ => {}
            }
        }
        for n in variants {
            match n % 5 {
                0 => rows.push(format!("extra_in_rust 1:{n}-{} A T x", n + 2)),
                1 => rows.push(format!("extra_in_rust 2:{n} A T y")),
                2 => rows.push(format!("extra_in_rust 3:{n} A Tb x")),
                4 => rows.push(format!("extra_in_rust 5:{n} A T x")),
                _ => {}
            }
        }
    }
    let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
    let one = classes_of(&rows, 1);
    let count = |name: &str| one.iter().filter(|class| **class == name).count();
    assert_eq!(count("coordinate_only"), 1600);
    assert_eq!(count("consequence_swap"), 1600);
    assert_eq!(count("transcript_selection"), 1600);
    assert_eq!(count("true_missing"), 800);
    assert_eq!(count("true_extra"), 800);
    for threads in [2, 3, 4, 7, 64] {
        assert_eq!(classes_of(&rows, threads), one, "{threads} threads");
    }
}

#[test]
fn rows_that_collide_on_every_key_are_classed_alike_whichever_half_is_held_whole() {
    // Values from small sets, drawn by a fixed sequence: many rows share each key, on
    // both sides, in tables where the reference-only rows are the few and the many.
    for reference_in_ten in [1, 5, 9] {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut draw = |below: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % below
        };
        let rows: Vec<String> = (0..6000)
            .map(|_| {
                let source = if draw(10) < reference_in_ten {
                    "missing_in_rust"
                } else {
                    "extra_in_rust"
                };
                let end = ["", "-15", "-25"][draw(3) as usize];
                format!(
                    "{source} {}:{}{end} A{} T{} c{}",
                    draw(2),
                    10 * (1 + draw(40)),
                    draw(2),
                    draw(30),
                    draw(3)
                )
            })
            .collect();
        let rows: Vec<&str> = rows.iter().map(String::as_str).collect();
        let one = classes_of(&rows, 1);
        for name in Class::NAMED.map(|class| class.name().unwrap_or_default()) {
            assert!(one.contains(&name), "no row is {name}");
        }
        assert_eq!(classes_of(&rows, 4), one);
    }
}
