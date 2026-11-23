// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The intron-end and intron-start insertion classes on one insertion's row in the four
//! formats: each case through the comparison, and the same rows through the two rules
//! alone.

mod compare_common;
mod intron_ends_common;

use std::path::Path;

use compare_common::at;
use compare_common::compare;
use compare_common::count_map;
use compare_common::counts;
use compare_common::default_row;
use compare_common::default_text;
use compare_common::documented_with;
use compare_common::scratch;
use compare_common::texts as text_items;
use compare_common::with_member;
use intron_ends_common::by_the_rule;
use intron_ends_common::compared;
use intron_ends_common::difference;
use intron_ends_common::dumped;
use intron_ends_common::green;
use intron_ends_common::intron_end_cache;
use intron_ends_common::intron_end_transcripts;
use intron_ends_common::placed;
use intron_ends_common::red;
use intron_ends_common::row;
use intron_ends_common::set_aside;
use intron_ends_common::text_rows;
use intron_ends_common::Verdicts;
use intron_ends_common::INTRON_END;
use intron_ends_common::INTRON_END_PROTEIN;
use intron_ends_common::INTRON_START;
use std::collections::BTreeSet;
use std::io::Write;
use vep_formats_compare::compare::compare_entries;
use vep_formats_compare::compare::CompareOptions;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::pyfmt;
use vep_formats_compare::pyjson::Value;
use vep_formats_compare::rules::intron_ends::clipped_insertions;
use vep_formats_compare::rules::intron_ends::coding_transcript_exons;
use vep_formats_compare::rules::intron_ends::forward_readings;
use vep_formats_compare::rules::intron_ends::held_insertion_gaps;
use vep_formats_compare::rules::intron_ends::intron_end_insertion_form;
use vep_formats_compare::rules::intron_ends::lone_insertion_gaps;
use vep_formats_compare::rules::intron_ends::notation_duplicated_positions;
use vep_formats_compare::rules::intron_ends::notation_inserted_bases;
use vep_formats_compare::rules::intron_ends::notation_insertion_gaps;
use vep_formats_compare::rules::intron_ends::transcript_intron_end_gaps;
use vep_formats_compare::rules::intron_ends::transcript_intron_start_gaps;
use vep_formats_compare::rules::intron_ends::turned_back;
use vep_formats_compare::rules::intron_ends::IntronEnds;
use vep_formats_compare::rules::intron_ends::INTRON_END_INSERTION_RE;
use vep_formats_compare::rules::intron_ends::INTRON_START_INSERTION_RE;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::ABSENT;

/// One case read both ways: through the comparison and through the two rules alone.
macro_rules! both_ways {
    ($compared:ident, $alone:ident, $case:ident) => {
        #[test]
        fn $compared() {
            $case(compared, stringify!($compared));
        }

        #[test]
        fn $alone() {
            $case(by_the_rule, stringify!($alone));
        }
    };
}

fn set(positions: &[i64]) -> BTreeSet<i64> {
    positions.iter().copied().collect()
}

fn lines(items: &[&str]) -> Vec<String> {
    items.iter().map(|item| item.to_string()).collect()
}

/// An insertion written between an intron's last base and the first coding base of the next
/// exon that copies the intron's last base, with no HGVS_OFFSET: the row is the class in every
/// format, counted once, on a forward-strand and on a reverse-strand transcript and at the end
/// of an intron of 14 bases; without the class, and under the class without the candidate's
/// cache, the pair stays a difference.
fn a_copy_of_the_intron_last_base(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let rows = [
        row("ENST00000000001", 11199, "G", "G", "c.51-1dup"),
        row("ENST00000000002", 25299, "G", "C", "c.51-1dup").strand(-1),
        row("ENST00000000004", 41313, "G", "G", "c.286-1dup"),
    ];
    for row in &rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            set_aside(INTRON_END)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
        assert_eq!(verdicts(None, &texts, &[INTRON_END]), difference());
    }
}

#[test]
fn a_copy_of_the_intron_last_base_is_set_aside_under_the_intron_end_class() {
    let name = "a_copy_of_the_intron_last_base_is_set_aside_under_the_intron_end_class";
    a_copy_of_the_intron_last_base(compared, name);
    let cache = intron_end_cache(&scratch(&format!("{name}_report")));
    let (reference, candidate, _) =
        row("ENST00000000001", 11199, "G", "G", "c.51-1dup").text("default");
    let side = Side {
        intron_ends: Some(IntronEnds::new(&cache, None).unwrap()),
        ..Side::default()
    };
    let report = compare_entries(
        &parse_default(&reference).unwrap(),
        &parse_default(&candidate).unwrap(),
        &mut documented_with(&[INTRON_END]),
        &CompareOptions::new(&side),
    )
    .unwrap();
    for column in ["HGVSp", "Extra"] {
        assert_eq!(
            dumped(at(&report, &["columns", column, "masked"])),
            format!(r#"{{"{INTRON_END}": 1}}"#)
        );
    }
    let example = format!(
        "1:11199-11200|G|ENST00000000001|Transcript: HGVSp expected {}, got {}",
        pyfmt::repr_str(INTRON_END_PROTEIN),
        pyfmt::repr_str(ABSENT)
    );
    assert!(text_items(at(&report, &["mask_examples", INTRON_END])).contains(&example.as_str()));
}

#[test]
fn the_intron_end_rule_takes_a_copy_of_the_intron_last_base() {
    a_copy_of_the_intron_last_base(
        by_the_rule,
        "the_intron_end_rule_takes_a_copy_of_the_intron_last_base",
    );
}

/// An insertion written in the same gap that copies two or more of the intron's last bases,
/// its first position numbered back from the exon's first coding base or, in the intron's
/// first half, on from the coding base before the intron, is the class in every format on
/// either strand. A span of one position, a duplication that ends before the intron's last
/// base, one numbered back from another base than the exon's first, and one numbered on from
/// another base than the one before it stay differences.
fn a_duplication_ending_on_the_intron_last_base(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let rows = [
        row("ENST00000000001", 11199, "G", "CAG", "c.51-3_51-1dup"),
        row("ENST00000000002", 25299, "G", "CTG", "c.51-3_51-1dup").strand(-1),
        row(
            "ENST00000000004",
            41313,
            "G",
            "TTCTTCAG",
            "c.285+7_286-1dup",
        ),
    ];
    for row in &rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            set_aside(INTRON_END)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    for notation in [
        "c.51-1_51-1dup",
        "c.51-3_51-2dup",
        "c.52-3_51-1dup",
        "c.49+98_51-1dup",
    ] {
        let texts = row("ENST00000000001", 11199, "G", "CAG", notation).texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            difference(),
            "{notation}"
        );
    }
}

both_ways!(
    a_duplication_ending_on_the_intron_last_base_is_set_aside_under_the_intron_end_class,
    the_intron_end_rule_takes_a_duplication_ending_on_the_intron_last_base,
    a_duplication_ending_on_the_intron_last_base
);

/// An insertion written in the same gap that is no copy of the bases before it is the class in
/// every format on either strand when it is one base, or two or more bases that do not end in
/// one of the two pairs that end an intron. An insertion between two other positions stays a
/// difference.
fn an_insertion_that_cannot_end_the_intron(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let rows = [
        row("ENST00000000001", 11199, "G", "T", "c.51-1_51insT"),
        row("ENST00000000001", 11199, "G", "TAT", "c.51-1_51insTAT"),
        row("ENST00000000002", 25299, "G", "ATA", "c.51-1_51insTAT").strand(-1),
        row("ENST00000000003", 31399, "G", "T", "c.51-1_51insT"),
    ];
    for row in &rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            set_aside(INTRON_END)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    for notation in ["c.51-1_52insT", "c.51-2_51-1insT"] {
        let texts = row("ENST00000000001", 11199, "G", "T", notation).texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            difference(),
            "{notation}"
        );
    }
}

both_ways!(
    an_insertion_that_cannot_end_the_intron_is_set_aside_under_the_intron_end_class,
    the_intron_end_rule_takes_an_insertion_that_cannot_end_the_intron,
    an_insertion_that_cannot_end_the_intron
);

/// An insertion written in the same gap that is two or more bases ending in one of the two
/// pairs that end an intron is the class as well, in every format on either strand, with
/// either pair, with the pair alone and at the end of an intron of 14 bases; without the class
/// each row stays a difference. The same bases inserted between two other positions stay a
/// difference. The form is named apart from the insertion of one base and of bases that end in
/// no pair.
fn an_insertion_that_ends_in_an_intron_end_pair(verdicts: Verdicts, name: &str) {
    for (body, form) in [
        ("c.51-1_51insTAG", "pair"),
        ("c.51-1_51insTAC", "pair"),
        ("c.51-1_51insAG", "pair"),
        ("c.51-1_51insTAT", "insertion"),
        ("c.51-1_51insG", "insertion"),
    ] {
        assert_eq!(intron_end_insertion_form(body).as_deref(), Some(form));
    }
    let cache = intron_end_cache(&scratch(name));
    let rows = [
        row("ENST00000000001", 11199, "G", "TAG", "c.51-1_51insTAG"),
        row("ENST00000000001", 11199, "G", "TAC", "c.51-1_51insTAC"),
        row("ENST00000000001", 11199, "G", "AG", "c.51-1_51insAG"),
        row("ENST00000000002", 25299, "G", "CTA", "c.51-1_51insTAG").strand(-1),
        row("ENST00000000002", 25299, "G", "GTA", "c.51-1_51insTAC").strand(-1),
        row("ENST00000000004", 41313, "G", "TAC", "c.286-1_286insTAC"),
    ];
    for row in &rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            set_aside(INTRON_END)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    for notation in ["c.51-1_52insTAG", "c.51-2_51-1insTAC"] {
        let texts = row("ENST00000000001", 11199, "G", "TAG", notation).texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            difference(),
            "{notation}"
        );
    }
}

both_ways!(
    an_insertion_that_ends_in_an_intron_end_pair_is_set_aside_under_the_intron_end_class,
    the_intron_end_rule_takes_an_insertion_that_ends_in_an_intron_end_pair,
    an_insertion_that_ends_in_an_intron_end_pair
);

/// A record that carries another allele beside an insertion holds the insertion with the
/// reference base it shares: its Location is that base's one position and its Allele the base
/// followed by the inserted bases. The row is the class in every format, counted once, for
/// three inserted bases and for one, on a reverse-strand transcript and at the end of an
/// intron of 14 bases. A span beside no gap of the class, an Allele whose base count or bases
/// do not fit the notation and the end of an intron of 13 bases stay differences, and so does
/// a VCF or JSON entry whose input line's ALT does not clip to the insertion. Beside a
/// duplication's notation the held text row stays a difference, and its entry, placed by the
/// ALT that clips to as many bases as the notation duplicates, is the class.
fn an_intron_end_insertion_the_record_holds_with_a_shared_base(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let others = "GATA,C";
    let held = |position: i64, allele: &'static str, notation: &'static str| {
        row("ENST00000000002", position, "G", allele, notation).strand(-1)
    };
    let rows = [
        held(25299, "GATA", "c.51-1_51insTAT")
            .location("1:25299")
            .record_alternates(others),
        held(25299, "GA", "c.51-1_51insT")
            .location("1:25299")
            .record_alternates("GA,C"),
        row("ENST00000000004", 41313, "G", "GTAT", "c.286-1_286insTAT")
            .location("1:41313")
            .record_alternates("GTAT,C"),
    ];
    for row in &rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            set_aside(INTRON_END)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
        assert_eq!(verdicts(None, &texts, &[INTRON_END]), difference());
    }
    let shapes = [
        (
            "a span beside no gap",
            held(25298, "GATA", "c.51-1_51insTAT")
                .location("1:25298")
                .record_alternates(others),
        ),
        (
            "a base count that does not fit",
            held(25299, "GATA", "c.51-1_51insT")
                .location("1:25299")
                .record_alternates(others),
        ),
        (
            "bases of the other strand",
            held(25299, "GTAT", "c.51-1_51insTAT")
                .location("1:25299")
                .record_alternates("GTAT,C"),
        ),
        (
            "an intron of 13 bases",
            row("ENST00000000004", 41112, "G", "GTAT", "c.101-1_101insTAT")
                .location("1:41112")
                .record_alternates("GTAT,C"),
        ),
    ];
    for (shape, row) in &shapes {
        assert_eq!(
            verdicts(Some(&cache), &row.texts(), &[INTRON_END]),
            difference(),
            "{shape}"
        );
    }
    let unclipped = held(25299, "GATA", "c.51-1_51insTAT")
        .location("1:25299")
        .record_reference("C")
        .record_alternates("GATA,G");
    let duplication = held(25299, "GA", "c.51-1dup")
        .location("1:25299")
        .record_alternates("GA,C");
    assert_eq!(
        verdicts(Some(&cache), &duplication.texts(), &[INTRON_END]),
        placed(red(), green(INTRON_END))
    );
    assert_eq!(
        verdicts(Some(&cache), &unclipped.texts(), &[INTRON_END]),
        placed(green(INTRON_END), red())
    );
}

both_ways!(
    an_intron_end_insertion_a_record_holds_with_a_shared_base_is_set_aside,
    the_intron_end_rule_takes_an_insertion_a_record_holds_with_a_shared_base,
    an_intron_end_insertion_the_record_holds_with_a_shared_base
);

#[test]
fn held_insertion_gaps_names_each_clipping_that_leaves_the_inserted_bases() {
    assert_eq!(held_insertion_gaps(100, 1, "GTAT", "TAT", "1"), [100]);
    assert_eq!(held_insertion_gaps(100, 1, "TATG", "TAT", "1"), [99]);
    assert_eq!(held_insertion_gaps(100, 1, "GATA", "TAT", "-1"), [100]);
    assert_eq!(held_insertion_gaps(100, 1, "GATA", "TAT", "-"), [100]);
    assert_eq!(held_insertion_gaps(100, 1, "GTAT", "TAT", "-"), [100]);
    assert!(held_insertion_gaps(100, 1, "GATA", "TAT", "1").is_empty());
    assert!(held_insertion_gaps(100, 1, "GTAT", "TAT", "-1").is_empty());
    assert_eq!(held_insertion_gaps(100, 1, "GGGG", "GGG", "1"), [99, 100]);
    assert_eq!(held_insertion_gaps(100, 2, "CGTAT", "TAT", "1"), [101]);
    let clipped = |pairs: &[(i64, &str)]| -> Vec<(i64, String)> {
        pairs
            .iter()
            .map(|(gap, bases)| (*gap, bases.to_string()))
            .collect()
    };
    assert_eq!(
        clipped_insertions(100, "G", "GGGG"),
        clipped(&[(99, "GGG"), (100, "GGG")])
    );
    assert_eq!(
        clipped_insertions(100, "G", "GTAT"),
        clipped(&[(100, "TAT")])
    );
    assert!(clipped_insertions(100, "C", "GTAT").is_empty());
    assert_eq!(
        clipped_insertions(100, "CG", "CGTAT"),
        clipped(&[(101, "TAT")])
    );
    assert_eq!(clipped_insertions(100, "", "TAT"), clipped(&[(99, "TAT")]));
    assert!(clipped_insertions(100, "CG", "G").is_empty());
    assert!(clipped_insertions(100, "G", "G").is_empty());
    assert_eq!(forward_readings("TAT", "1"), ["TAT"]);
    assert_eq!(forward_readings("TAT", "-1"), ["ATA"]);
    assert_eq!(forward_readings("TAT", "-"), ["TAT", "ATA"]);
    assert!(held_insertion_gaps(100, 1, "GTAT", "T", "1").is_empty());
    assert!(held_insertion_gaps(100, 2, "GTAT", "TAT", "1").is_empty());
    assert!(held_insertion_gaps(100, 1, "TAT", "TAT", "1").is_empty());
    assert!(held_insertion_gaps(100, 1, "GTAT", "", "1").is_empty());
    assert_eq!(
        notation_inserted_bases(&INTRON_END_INSERTION_RE, "c.51-1_51insTAT"),
        "TAT"
    );
    assert_eq!(
        notation_inserted_bases(&INTRON_START_INSERTION_RE, "c.50_50+1insT"),
        "T"
    );
    assert_eq!(
        notation_inserted_bases(&INTRON_END_INSERTION_RE, "c.51-1dup"),
        ""
    );
}

/// A record that carries another allele beside the insertion keeps its own Location in the
/// default and tab formats, and the class places the row by that span: every gap inside and
/// beside it. A minimised Allele is the insertion alone, the bases an `ins` notation inserts
/// on either strand or, for a `dup` notation, as many bases as the positions it duplicates,
/// and the row is the class in every format. A held Allele the row spells on the reverse
/// strand is the class in default and tab. A minimised Allele that is not the notation's
/// bases, a `dup` row whose Allele's base count does not fit, a duplication that starts in the
/// intron's first half and a span beside no gap of the class stay differences there.
fn an_intron_end_text_row_placed_by_its_location_span(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let end_notation = "c.286-1_286insTAT";
    let minimised = [
        row("ENST00000000004", 41313, "G", "TAT", end_notation)
            .location("1:41313")
            .record_alternates("GTAT,C"),
        row("ENST00000000004", 41313, "G", "TAT", end_notation)
            .location("1:41314")
            .record_alternates("GTAT,C"),
        row("ENST00000000002", 25299, "G", "ATA", "c.51-1_51insTAT")
            .strand(-1)
            .location("1:25299")
            .record_alternates("GATA,C"),
        row("ENST00000000001", 11199, "G", "G", "c.51-1dup")
            .location("1:11199")
            .record_alternates("GG,C"),
        row("ENST00000000001", 11199, "G", "CAG", "c.51-3_51-1dup")
            .location("1:11199")
            .record_alternates("GCAG,C"),
        row("ENST00000000001", 11199, "G", "CAG", "c.51-3_51-1dup")
            .location("1:11198-11200")
            .record_alternates("GCAG,C"),
    ];
    for row in &minimised {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            set_aside(INTRON_END)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    let reverse_strand_records = [
        row("ENST00000000004", 41313, "G", "ATAC", end_notation).location("1:41313"),
        row("ENST00000000002", 25299, "G", "TATC", "c.51-1_51insTAT")
            .strand(-1)
            .location("1:25299"),
    ];
    for row in &reverse_strand_records {
        let verdicts = verdicts(Some(&cache), &row.texts(), &[INTRON_END]);
        assert_eq!(text_rows(&verdicts), [&green(INTRON_END); 2]);
    }
    let shapes = [
        (
            "other bases than the notation's",
            row("ENST00000000004", 41313, "G", "CCG", end_notation).location("1:41313"),
        ),
        (
            "three bases beside a copy of one",
            row("ENST00000000001", 11199, "G", "CAG", "c.51-1dup").location("1:11199"),
        ),
        (
            "one base beside a duplication of three",
            row("ENST00000000001", 11199, "G", "G", "c.51-3_51-1dup").location("1:11199"),
        ),
        (
            "a duplication from the intron's first half",
            row(
                "ENST00000000004",
                41313,
                "G",
                "TTCTTCAG",
                "c.285+7_286-1dup",
            )
            .location("1:41313"),
        ),
        (
            "a span below the gap",
            row("ENST00000000004", 41313, "G", "TAT", end_notation).location("1:41312"),
        ),
        (
            "a span above the gap",
            row("ENST00000000004", 41313, "G", "TAT", end_notation).location("1:41315"),
        ),
    ];
    for (shape, row) in &shapes {
        let verdicts = verdicts(Some(&cache), &row.texts(), &[INTRON_END]);
        assert_eq!(text_rows(&verdicts), [&red(); 2], "{shape}");
    }
}

both_ways!(
    an_intron_end_text_row_is_placed_by_its_location_span,
    the_intron_end_rule_places_a_text_row_by_its_location_span,
    an_intron_end_text_row_placed_by_its_location_span
);

#[test]
fn lone_insertion_gaps_names_every_gap_of_the_span() {
    assert_eq!(lone_insertion_gaps(100, 1, "TAT", "TAT", 0, 0), [99, 100]);
    assert_eq!(lone_insertion_gaps(100, 1, "ATA", "TAT", 0, 0), [99, 100]);
    assert_eq!(
        lone_insertion_gaps(100, 3, "TAT", "TAT", 0, 0),
        [99, 100, 101, 102]
    );
    assert!(lone_insertion_gaps(100, 1, "CCG", "TAT", 0, 0).is_empty());
    assert!(lone_insertion_gaps(100, 1, "GTAT", "TAT", 0, 0).is_empty());
    assert_eq!(lone_insertion_gaps(100, 1, "TAT", "TTA", 0, 2), [99, 100]);
    assert_eq!(lone_insertion_gaps(100, 1, "TAT", "TTA", 0, -5), [99, 100]);
    assert_eq!(lone_insertion_gaps(100, 1, "TAT", "TAT", 0, 2), [99, 100]);
    assert_eq!(lone_insertion_gaps(100, 1, "ATA", "TAT", 0, -2), [99, 100]);
    assert!(lone_insertion_gaps(100, 1, "TAT", "TTA", 0, 0).is_empty());
    assert!(lone_insertion_gaps(100, 1, "ATT", "TTA", 0, 2).is_empty());
    assert_eq!(lone_insertion_gaps(100, 1, "CAG", "", 3, 0), [99, 100]);
    assert_eq!(lone_insertion_gaps(100, 2, "G", "", 1, 0), [99, 100, 101]);
    assert!(lone_insertion_gaps(100, 1, "CAG", "", 1, 0).is_empty());
    assert!(lone_insertion_gaps(100, 1, "CAG", "", 0, 0).is_empty());
    assert!(lone_insertion_gaps(100, 1, "CAG", "", 3, 1).is_empty());
    assert!(lone_insertion_gaps(100, 0, "TAT", "TAT", 0, 0).is_empty());
    assert!(lone_insertion_gaps(100, 1, "DELETION", "", 8, 0).is_empty());
    assert!(lone_insertion_gaps(100, 1, "-", "", 1, 0).is_empty());
    let bodies = [
        "c.51-1dup",
        "c.51-3_51-1dup",
        "c.50+7_51-1dup",
        "c.51+3_51-1dup",
        "c.52-3_51-1dup",
        "c.51-1_51insTAT",
        "c.50_50+1insTAT",
    ];
    assert_eq!(
        bodies.map(notation_duplicated_positions),
        [1, 3, 0, 0, 0, 0, 0]
    );
}

/// The shapes beside the intron-end class stay differences, the class applied to none of
/// them. A text-format row is placed by its Location, whose span must lie over or beside a gap
/// of the class, and a VCF or JSON entry by its input line, whose ALT must clip against the
/// REF to the entry's allele in such a gap.
fn the_shapes_nearest_the_intron_end_class(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let copy = || row("ENST00000000001", 11199, "G", "G", "c.51-1dup");
    assert_eq!(
        verdicts(Some(&cache), &copy().texts(), &[INTRON_END]),
        set_aside(INTRON_END)
    );
    let fs_ter9 = "ENSP00000000001.1:p.Gly17ArgfsTer9";
    assert_eq!(fs_ter9, INTRON_END_PROTEIN.replace("fsTer5", "fsTer9"));
    let other = "ENST00000000002.1:c.51-1dup";
    let shapes = [
        (
            "ends in a pair, with an offset",
            row("ENST00000000001", 11199, "G", "TAG", "c.51-1_51insTAG").both("HGVS_OFFSET", "1"),
        ),
        (
            "ends in a pair, an intron of 13 bases",
            row("ENST00000000004", 41112, "G", "TAG", "c.101-1_101insTAG"),
        ),
        (
            "ends in a pair, an intron of 2 bases",
            row("ENST00000000004", 41201, "G", "TAC", "c.188-1_188insTAC"),
        ),
        (
            "written away from the gap",
            row("ENST00000000001", 11198, "A", "G", "c.51-1dup").both("HGVS_OFFSET", "1"),
        ),
        ("moved, written in the gap", copy().both("HGVS_OFFSET", "1")),
        (
            "an offset on the reference alone",
            copy().reference("HGVS_OFFSET", "1"),
        ),
        (
            "an offset on the candidate alone",
            copy().candidate("HGVS_OFFSET", "1"),
        ),
        (
            "an exonic notation",
            row("ENST00000000001", 11249, "G", "T", "c.100_101insT"),
        ),
        (
            "before a 5' UTR position",
            row("ENST00000000003", 31199, "G", "G", "c.-50-1dup"),
        ),
        (
            "before a 3' UTR position",
            row("ENST00000000003", 31599, "G", "G", "c.*51-1dup"),
        ),
        (
            "two notations that differ",
            copy().candidate("HGVSc", "ENST00000000001.1:c.51-2_51-1dup"),
        ),
        ("another transcript's notation", copy().both("HGVSc", other)),
        (
            "an intron of 13 bases",
            row("ENST00000000004", 41112, "G", "G", "c.101-1dup"),
        ),
        (
            "an intron of 2 bases",
            row("ENST00000000004", 41201, "G", "G", "c.188-1dup"),
        ),
        (
            "a gap before a UTR exon base",
            row("ENST00000000003", 31199, "G", "G", "c.51-1dup"),
        ),
        (
            "another transcript's gap",
            row("ENST00000000001", 25299, "G", "G", "c.51-1dup"),
        ),
        (
            "a transcript the cache lacks",
            row("ENST00000000009", 11199, "G", "G", "c.51-1dup"),
        ),
        (
            "another candidate notation",
            copy().candidate("HGVSp", fs_ter9),
        ),
        (
            "no reference notation",
            copy()
                .reference("HGVSp", "")
                .candidate("HGVSp", INTRON_END_PROTEIN),
        ),
    ];
    for (shape, row) in &shapes {
        assert_eq!(
            verdicts(Some(&cache), &row.texts(), &[INTRON_END]),
            difference(),
            "{shape}"
        );
    }
    for location in ["1:11197-11198", "1:11201-11202", "1:11198", "1:11201-11203"] {
        assert_eq!(
            verdicts(
                Some(&cache),
                &copy().location(location).texts(),
                &[INTRON_END]
            ),
            placed(red(), green(INTRON_END)),
            "{location}"
        );
    }
    assert_eq!(
        verdicts(
            Some(&cache),
            &copy().record_reference("GC").texts(),
            &[INTRON_END]
        ),
        placed(green(INTRON_END), red())
    );
}

both_ways!(
    the_shapes_nearest_the_intron_end_class_stay_differences,
    the_intron_end_rule_leaves_the_nearest_shapes,
    the_shapes_nearest_the_intron_end_class
);

/// The transcript `transcript` under the stable id `ENST0000000000<number>`, its exons `by`
/// bases further along the chromosome.
fn moved(transcript: &Value, number: i64, by: i64) -> Value {
    let exons = transcript.get("exons").and_then(Value::as_list).unwrap();
    let whole = |exon: &Value, name: &str| match exon.get(name) {
        Some(Value::Int(position)) => Value::Int(position + by),
        other => panic!("{name} is {other:?}"),
    };
    let exons = exons
        .iter()
        .map(|exon| {
            let exon = with_member(exon, "start", whole(exon, "start"));
            with_member(&exon, "end", whole(&exon, "end"))
        })
        .collect();
    let transcript = with_member(
        transcript,
        "stable_id",
        Value::Str(format!("ENST0000000000{number}")),
    );
    with_member(&transcript, "exons", Value::List(exons))
}

/// `transcript` with its strand, its coding start and its exons' bounds as text.
fn as_text(transcript: &Value) -> Value {
    let text = |value: &Value| match value {
        Value::Int(number) => Value::Str(number.to_string()),
        other => panic!("{other:?} is no whole number"),
    };
    let exons = transcript
        .get("exons")
        .and_then(Value::as_list)
        .unwrap()
        .iter()
        .map(|exon| {
            let start = text(exon.get("start").unwrap());
            let exon = with_member(exon, "start", start);
            let end = text(exon.get("end").unwrap());
            with_member(&exon, "end", end)
        })
        .collect();
    let transcript = with_member(transcript, "strand", Value::Str("-1".to_string()));
    let transcript = with_member(
        &transcript,
        "cdna_coding_start",
        Value::Str("51".to_string()),
    );
    with_member(&transcript, "exons", Value::List(exons))
}

/// `transcript` without its two coding coordinates.
fn without_coding(transcript: &Value) -> Value {
    Value::Dict(
        transcript
            .as_dict()
            .unwrap()
            .iter()
            .filter(|(name, _)| !name.starts_with("cdna_"))
            .map(|(name, value)| (name.clone(), value.clone()))
            .collect(),
    )
}

/// The names of the shards of chromosome `directory` not read so far.
fn unread(ends: &IntronEnds, directory: &str) -> Vec<String> {
    ends._unread.borrow()[directory]
        .iter()
        .map(|(_, _, shard)| shard.file_name().unwrap().to_string_lossy().into_owned())
        .collect()
}

/// The stable ids whose gaps are held, in the order their shards were read.
fn held_ids(ends: &IntronEnds) -> Vec<String> {
    ends._gaps
        .borrow()
        .keys()
        .map(|(_, stable_id)| stable_id.clone())
        .collect()
}

/// Leaves text that is no JSON in every shard of `directory`: a shard read again ends the
/// lookup.
fn spoil(directory: &Path, shards: &[&str]) {
    for shard in shards {
        std::fs::write(directory.join(shard), "read once").unwrap();
    }
}

#[test]
fn intron_ends_reads_the_shard_at_the_row_first_and_no_shard_twice() {
    let [forward, reverse, utr, short] = intron_end_transcripts();
    assert_eq!(transcript_intron_end_gaps(&forward), set(&[11199, 11499]));
    assert_eq!(transcript_intron_end_gaps(&reverse), set(&[25099, 25299]));
    assert_eq!(transcript_intron_end_gaps(&utr), set(&[31399]));
    assert_eq!(transcript_intron_end_gaps(&short), set(&[41313]));
    assert_eq!(
        transcript_intron_end_gaps(&as_text(&reverse)),
        set(&[25099, 25299])
    );
    assert!(transcript_intron_end_gaps(&without_coding(&forward)).is_empty());
    let coding_past_the_exons = with_member(&forward, "cdna_coding_end", Value::Int(301));
    assert!(transcript_intron_end_gaps(&coding_past_the_exons).is_empty());
    let an_exon_without_its_end = with_member(
        &forward,
        "exons",
        compare_common::json(r#"[{"start": 11000}]"#),
    );
    assert!(transcript_intron_end_gaps(&an_exon_without_its_end).is_empty());

    let root = scratch("intron_ends_reads_the_shard_at_the_row_first_and_no_shard_twice");
    let cache = root.join("cache");
    let directory = cache.join("transcripts").join("7");
    std::fs::create_dir_all(&directory).unwrap();
    let first_shard = Value::List(vec![
        Value::Null,
        moved(&forward, 5, 0),
        with_member(&utr, "stable_id", Value::Str("GENE1.1".to_string())),
        with_member(&short, "stable_id", Value::Str("GENE1.2".to_string())),
    ]);
    std::fs::write(directory.join("1-1000000.json"), dumped(&first_shard)).unwrap();
    let second_shard = format!(
        r#"{{"1000001-2000000": [{}, null]}}"#,
        dumped(&moved(&forward, 6, 1_500_000))
    );
    let mut gzip = flate2::write::GzEncoder::new(
        std::fs::File::create(directory.join("1000001-2000000.json.gz")).unwrap(),
        flate2::Compression::default(),
    );
    gzip.write_all(second_shard.as_bytes()).unwrap();
    gzip.finish().unwrap();
    let third_shard = Value::List(vec![moved(&forward, 7, 2_500_000)]);
    std::fs::write(directory.join("2000001-3000000.json"), dumped(&third_shard)).unwrap();
    std::fs::write(directory.join("notes.json"), "not a shard").unwrap();

    let ends = IntronEnds::new(&cache, None).unwrap();
    assert_eq!(
        ends.gaps("7", "ENST00000000006", 1_511_199).unwrap(),
        set(&[1_511_199, 1_511_499])
    );
    assert_eq!(
        unread(&ends, "7"),
        ["1-1000000.json", "2000001-3000000.json"]
    );
    assert_eq!(
        ends.gaps("7", "ENST00000000007.3", 1_200_000).unwrap(),
        set(&[2_511_199, 2_511_499])
    );
    // The shard nearest the position is read before the one further away.
    assert!(unread(&ends, "7").is_empty());
    assert_eq!(
        held_ids(&ends),
        [
            "ENST00000000006",
            "ENST00000000005",
            "GENE1.1",
            "GENE1.2",
            "ENST00000000007"
        ]
    );
    spoil(
        &directory,
        &[
            "1-1000000.json",
            "1000001-2000000.json.gz",
            "2000001-3000000.json",
        ],
    );
    assert_eq!(
        ends.gaps("chr7", "ENST00000000005", 11_199).unwrap(),
        set(&[11_199, 11_499])
    );
    assert_eq!(ends.gaps("7", "GENE1.1", 31_399).unwrap(), set(&[31_399]));
    assert_eq!(ends.gaps("7", "GENE1.2", 41_313).unwrap(), set(&[41_313]));
    assert!(ends
        .gaps("7", "ENST00000000009", 11_199)
        .unwrap()
        .is_empty());
    assert!(ends
        .gaps("8", "ENST00000000005", 11_199)
        .unwrap()
        .is_empty());
    assert_eq!(held_ids(&ends).len(), 5);
    let row = parse_default(&default_text(&[default_row(
        "rec1",
        "chr7:2511199-2511200",
        "G",
        "ENST00000000007.3",
        "frameshift_variant",
        "IMPACT=HIGH",
    )]))
    .unwrap()
    .remove(0);
    assert!(ends.holds(&row, "", 0).unwrap());
    assert_eq!(
        ends.written_gap(&row),
        Some(("chr7".to_string(), 2_511_199))
    );
    std::fs::write(directory.join("3000001-4000000.json"), "not JSON").unwrap();
    let refused = IntronEnds::new(&cache, None)
        .unwrap()
        .gaps("7", "ENST00000000009", 3_500_000)
        .unwrap_err()
        .to_string();
    assert!(
        refused.contains("3000001-4000000.json is not JSON"),
        "{refused}"
    );
}

/// An insertion the 3' rule leaves between a coding exon's last base and the first base of the
/// intron after it, written in that gap or inside the exon and carried to it, with one
/// HGVS_OFFSET on both sides, positive on a forward-strand transcript and negative on a
/// reverse-strand one. The row is the class in every format, counted once, on either strand
/// and before an intron of 13 bases and of 2 bases; without the class, under the class without
/// the candidate's cache and under the intron-end class alone it stays a difference.
fn an_insertion_left_at_the_start_of_an_intron(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let rows = [
        row("ENST00000000001", 11099, "G", "T", "c.50_50+1insT"),
        row("ENST00000000001", 11299, "G", "TAT", "c.150_150+1insTAT"),
        row("ENST00000000002", 25499, "G", "ATA", "c.50_50+1insTAT").strand(-1),
        row("ENST00000000001", 11097, "G", "TAT", "c.50_50+1insTAT").both("HGVS_OFFSET", "2"),
        row("ENST00000000002", 25501, "G", "ATA", "c.50_50+1insTAT")
            .strand(-1)
            .both("HGVS_OFFSET", "-2"),
        row("ENST00000000004", 41099, "G", "T", "c.100_100+1insT"),
        row("ENST00000000004", 41199, "G", "TAT", "c.187_187+1insTAT"),
    ];
    for row in &rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_START]),
            set_aside(INTRON_START)
        );
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END, INTRON_START]),
            set_aside(INTRON_START)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
        assert_eq!(verdicts(None, &texts, &[INTRON_START]), difference());
        assert_eq!(verdicts(Some(&cache), &texts, &[INTRON_END]), difference());
    }
}

#[test]
fn an_insertion_left_at_the_start_of_an_intron_is_set_aside_under_the_intron_start_class() {
    let name =
        "an_insertion_left_at_the_start_of_an_intron_is_set_aside_under_the_intron_start_class";
    an_insertion_left_at_the_start_of_an_intron(compared, name);
    let cache = intron_end_cache(&scratch(&format!("{name}_report")));
    let (reference, candidate, _) =
        row("ENST00000000001", 11099, "G", "T", "c.50_50+1insT").text("default");
    let side = Side {
        intron_ends: Some(IntronEnds::new(&cache, None).unwrap()),
        ..Side::default()
    };
    let report = compare_entries(
        &parse_default(&reference).unwrap(),
        &parse_default(&candidate).unwrap(),
        &mut documented_with(&[INTRON_START]),
        &CompareOptions::new(&side),
    )
    .unwrap();
    for column in ["HGVSp", "Extra"] {
        assert_eq!(
            dumped(at(&report, &["columns", column, "masked"])),
            format!(r#"{{"{INTRON_START}": 1}}"#)
        );
    }
}

#[test]
fn the_intron_start_rule_takes_an_insertion_left_at_the_start_of_an_intron() {
    an_insertion_left_at_the_start_of_an_intron(
        by_the_rule,
        "the_intron_start_rule_takes_an_insertion_left_at_the_start_of_an_intron",
    );
}

/// A run that prints the moved location writes a text row's Location at the gap the 3' rule
/// leaves the insertion in, beside the same HGVS_OFFSET: the row is the class where the gap
/// its Location names is one in which a coding exon ends before an intron, while the VCF and
/// JSON entries, placed by the input line, are the class by the gap moved. A row with an
/// offset none of whose Location's gaps is such a gap, as named or moved, stays a difference,
/// and an entry placed by an input line written at the gap, beside an offset, is held to the
/// moved gap alone.
fn a_moved_record_located_at_the_gap_it_is_left_in(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let carried = |position: i64| {
        row("ENST00000000001", position, "G", "TAT", "c.50_50+1insTAT").both("HGVS_OFFSET", "2")
    };
    let rows = [
        carried(11097).location("1:11099-11100"),
        row("ENST00000000002", 25501, "G", "ATA", "c.50_50+1insTAT")
            .strand(-1)
            .location("1:25499-25500")
            .both("HGVS_OFFSET", "-2"),
    ];
    for row in &rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_START]),
            set_aside(INTRON_START)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
        assert_eq!(verdicts(None, &texts, &[INTRON_START]), difference());
    }
    assert_eq!(
        verdicts(Some(&cache), &carried(11095).texts(), &[INTRON_START]),
        difference()
    );
    assert_eq!(
        verdicts(Some(&cache), &carried(11099).texts(), &[INTRON_START]),
        placed(green(INTRON_START), red())
    );
}

both_ways!(
    a_moved_record_located_at_the_gap_it_is_left_in_is_set_aside,
    the_intron_start_rule_takes_a_moved_record_located_at_the_gap_it_is_left_in,
    a_moved_record_located_at_the_gap_it_is_left_in
);

/// A record that carries another allele beside an insertion holds the insertion with the
/// reference base it shares. On a forward-strand transcript whose exon ends on the shared
/// base the row is the class in every format, counted once, before an intron of 100 bases and
/// of 2 bases. A span beside no gap of the class, an Allele whose base count or bases do not
/// fit the notation and a row that carries an offset stay differences, and so does a VCF or
/// JSON entry whose allele is no ALT of its input line or does not clip to the insertion
/// there.
fn an_intron_start_insertion_the_record_holds_with_a_shared_base(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let others = "GTAT,C";
    let start_notation = "c.50_50+1insTAT";
    let held = |position: i64, allele: &'static str, notation: &'static str| {
        row("ENST00000000001", position, "G", allele, notation)
    };
    let rows = [
        held(11099, "GTAT", start_notation)
            .location("1:11099")
            .record_alternates(others),
        row("ENST00000000004", 41199, "G", "GT", "c.187_187+1insT")
            .location("1:41199")
            .record_alternates("GT,C"),
    ];
    for row in &rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_START]),
            set_aside(INTRON_START)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
        assert_eq!(verdicts(None, &texts, &[INTRON_START]), difference());
    }
    let shapes = [
        (
            "a span beside no gap",
            held(11097, "GTAT", start_notation)
                .location("1:11097")
                .record_alternates(others),
        ),
        (
            "a base count that does not fit",
            held(11099, "GTAT", "c.50_50+1insT")
                .location("1:11099")
                .record_alternates(others),
        ),
        (
            "bases of the other strand",
            held(11099, "GATA", start_notation)
                .location("1:11099")
                .record_alternates("GATA,C"),
        ),
        (
            "an offset",
            held(11099, "GTAT", start_notation)
                .location("1:11099")
                .record_alternates(others)
                .both("HGVS_OFFSET", "2"),
        ),
    ];
    for (shape, row) in &shapes {
        assert_eq!(
            verdicts(Some(&cache), &row.texts(), &[INTRON_START]),
            difference(),
            "{shape}"
        );
    }
    let unclipped = held(11099, "GTAT", start_notation)
        .location("1:11099")
        .record_reference("C")
        .record_alternates("GTAT,G");
    let no_alternate = held(11099, "GTAT", start_notation)
        .location("1:11099")
        .record_alternates("GT,C");
    for row in [unclipped, no_alternate] {
        assert_eq!(
            verdicts(Some(&cache), &row.texts(), &[INTRON_START]),
            placed(green(INTRON_START), red())
        );
    }
}

both_ways!(
    an_intron_start_insertion_a_record_holds_with_a_shared_base_is_set_aside,
    the_intron_start_rule_takes_an_insertion_a_record_holds_with_a_shared_base,
    an_intron_start_insertion_the_record_holds_with_a_shared_base
);

/// The intron-start class places a text row by its Location's span: every gap inside and
/// beside it. A minimised Allele is the insertion alone, the bases the notation inserts on
/// either strand, turned by the offset where the 3' rule moved the record, and the row is the
/// class in every format. A held Allele the row spells on the reverse strand is the class in
/// default and tab. A minimised Allele that is not the notation's bases, a moved Allele that
/// holds the shared base and a span beside no gap of the class stay differences there.
fn an_intron_start_text_row_placed_by_its_location_span(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let start_notation = "c.50_50+1insTAT";
    let minimised = [
        row("ENST00000000001", 11099, "G", "TAT", start_notation)
            .location("1:11099")
            .record_alternates("GTAT,C"),
        row("ENST00000000001", 11099, "G", "TAT", start_notation)
            .location("1:11100")
            .record_alternates("GTAT,C"),
        row("ENST00000000002", 25499, "G", "ATA", start_notation)
            .strand(-1)
            .location("1:25499")
            .record_alternates("GATA,C"),
        row("ENST00000000001", 11099, "G", "TAT", start_notation)
            .location("1:11098-11100")
            .record_alternates("GTAT,C"),
        row("ENST00000000001", 11097, "G", "TAT", "c.50_50+1insTTA")
            .location("1:11097")
            .record_alternates("GTAT,C")
            .both("HGVS_OFFSET", "2"),
    ];
    for row in &minimised {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_START]),
            set_aside(INTRON_START)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    let reverse_strand_records = [
        row("ENST00000000001", 11099, "G", "ATAC", start_notation).location("1:11099"),
        row("ENST00000000002", 25499, "G", "TATC", start_notation)
            .strand(-1)
            .location("1:25499"),
    ];
    for row in &reverse_strand_records {
        let verdicts = verdicts(Some(&cache), &row.texts(), &[INTRON_START]);
        assert_eq!(text_rows(&verdicts), [&green(INTRON_START); 2]);
    }
    let shapes = [
        (
            "other bases than the notation's",
            row("ENST00000000001", 11099, "G", "CCG", start_notation).location("1:11099"),
        ),
        (
            "moved, held with the shared base",
            row("ENST00000000001", 11097, "G", "GTAT", start_notation)
                .location("1:11097")
                .both("HGVS_OFFSET", "2"),
        ),
        (
            "a span below the gap",
            row("ENST00000000001", 11099, "G", "TAT", start_notation).location("1:11098"),
        ),
        (
            "a span above the gap",
            row("ENST00000000001", 11099, "G", "TAT", start_notation).location("1:11101"),
        ),
    ];
    for (shape, row) in &shapes {
        let verdicts = verdicts(Some(&cache), &row.texts(), &[INTRON_START]);
        assert_eq!(text_rows(&verdicts), [&red(); 2], "{shape}");
    }
}

both_ways!(
    an_intron_start_text_row_is_placed_by_its_location_span,
    the_intron_start_rule_places_a_text_row_by_its_location_span,
    an_intron_start_text_row_placed_by_its_location_span
);

/// A run that prints a moved record's allele as the 3' rule moved it writes the bases the
/// row's HGVSc inserts as they stand, beside the offset. A minimised row of a record that
/// carries another allele is the class in default and tab with the record's own Location and
/// with the moved Location; a VCF or JSON entry is placed by its input line, whose ALT inserts
/// the notation's bases turned back by the offset. A span whose gaps are a gap of the class
/// neither as named nor moved stays a difference. A text row whose allele is the notation's
/// bases under neither reading, or holds the shared base beside the offset, stays one while
/// its entry, placed by the input line, is the class. An entry stays one where its input
/// line's ALT inserts the notation's bases turned by another count, or turned beside no
/// offset.
fn a_moved_allele_printed_as_the_rule_moved_it(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let notation = "c.50_50+1insTTA";
    let (others, reverse_others) = ("GTAT,C", "GATA,C");
    let forward = |position: i64, allele: &'static str, notation: &'static str| {
        row("ENST00000000001", position, "G", allele, notation).both("HGVS_OFFSET", "2")
    };
    let reverse = || {
        row("ENST00000000002", 25501, "G", "TAA", notation)
            .strand(-1)
            .both("HGVS_OFFSET", "-2")
    };
    let rows = [
        forward(11097, "TTA", notation)
            .location("1:11097")
            .record_alternates(others),
        forward(11097, "TTA", notation)
            .location("1:11099")
            .record_alternates(others),
        forward(11097, "TTA", notation).record_alternates("GTAT"),
        reverse()
            .location("1:25501")
            .record_alternates(reverse_others),
        reverse()
            .location("1:25499")
            .record_alternates(reverse_others),
        reverse().record_alternates("GATA"),
    ];
    for row in &rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_START]),
            set_aside(INTRON_START)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    let neither_place = forward(11095, "TTA", notation)
        .location("1:11095")
        .record_alternates(others);
    assert_eq!(
        verdicts(Some(&cache), &neither_place.texts(), &[INTRON_START]),
        difference()
    );
    let text_rows_alone = [
        (
            "neither reading of the allele",
            forward(11097, "ATT", notation)
                .location("1:11097")
                .record_alternates(others),
        ),
        (
            "held with the shared base",
            forward(11097, "GTAT", "c.50_50+1insTAT")
                .location("1:11097")
                .record_alternates(others),
        ),
    ];
    for (shape, row) in &text_rows_alone {
        assert_eq!(
            verdicts(Some(&cache), &row.texts(), &[INTRON_START]),
            placed(red(), green(INTRON_START)),
            "{shape}"
        );
    }
    let entries = [
        (
            "turned by another count",
            forward(11097, "TTA", notation).record_alternates("GATT"),
        ),
        (
            "turned beside no offset",
            row("ENST00000000001", 11099, "G", "TTA", notation).record_alternates("GTAT"),
        ),
    ];
    for (shape, row) in &entries {
        assert_eq!(
            verdicts(Some(&cache), &row.texts(), &[INTRON_START]),
            placed(green(INTRON_START), red()),
            "{shape}"
        );
    }
}

both_ways!(
    a_moved_allele_printed_as_the_rule_moved_it_is_set_aside,
    the_intron_start_rule_takes_a_moved_allele_printed_as_the_rule_moved_it,
    a_moved_allele_printed_as_the_rule_moved_it
);

/// The shapes beside the intron-start class stay differences, the class applied to none of
/// them. A text-format row is placed by its Location, whose span must lie over or beside a gap
/// of the class, and a VCF or JSON entry by its input line, whose ALT must clip against the
/// REF to the entry's allele in such a gap.
fn the_shapes_nearest_the_intron_start_class(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let start_notation = "c.50_50+1insTAT";
    let left = || row("ENST00000000001", 11099, "G", "TAT", start_notation);
    let carried = || row("ENST00000000001", 11097, "G", "TAT", start_notation);
    assert_eq!(
        verdicts(Some(&cache), &left().texts(), &[INTRON_START]),
        set_aside(INTRON_START)
    );
    assert_eq!(
        verdicts(
            Some(&cache),
            &carried().both("HGVS_OFFSET", "2").texts(),
            &[INTRON_START]
        ),
        set_aside(INTRON_START)
    );
    let other = "ENST00000000002.1:c.50_50+1insTAT";
    let shapes = [
        (
            "a notation on both sides",
            left().candidate("HGVSp", "ENSP00000000001.1:p.Gly17ArgfsTer9"),
        ),
        (
            "no reference notation",
            left()
                .reference("HGVSp", "")
                .candidate("HGVSp", INTRON_END_PROTEIN),
        ),
        (
            "a duplication",
            row("ENST00000000001", 11099, "G", "G", "c.50dup"),
        ),
        (
            "between two coding positions",
            row("ENST00000000001", 11099, "G", "TAT", "c.49_50insTAT"),
        ),
        (
            "two different numbers",
            row("ENST00000000001", 11099, "G", "TAT", "c.50_51+1insTAT"),
        ),
        (
            "a 5' UTR position",
            row("ENST00000000003", 31099, "G", "TAT", "c.-51_-51+1insTAT"),
        ),
        (
            "a 3' UTR position",
            row("ENST00000000003", 31499, "G", "TAT", "c.*50_*50+1insTAT"),
        ),
        (
            "two notations that differ",
            left().candidate("HGVSc", "ENST00000000001.1:c.50_50+1insT"),
        ),
        ("another transcript's notation", left().both("HGVSc", other)),
        (
            "offsets that differ",
            carried()
                .reference("HGVS_OFFSET", "2")
                .candidate("HGVS_OFFSET", "3"),
        ),
        (
            "an offset on the reference alone",
            carried().reference("HGVS_OFFSET", "2"),
        ),
        (
            "an offset on the candidate alone",
            carried().candidate("HGVS_OFFSET", "2"),
        ),
        ("an offset of zero", left().both("HGVS_OFFSET", "0")),
        (
            "an offset that is no integer",
            left().both("HGVS_OFFSET", "x"),
        ),
        ("inside the exon, no offset", carried()),
        ("carried past the gap", carried().both("HGVS_OFFSET", "4")),
        ("carried the other way", carried().both("HGVS_OFFSET", "-2")),
        (
            "the intron's end",
            row("ENST00000000001", 11199, "G", "TAT", start_notation),
        ),
        (
            "the transcript's last exon",
            row("ENST00000000001", 11599, "G", "TAT", start_notation),
        ),
        (
            "another transcript's gap",
            row("ENST00000000001", 25499, "G", "TAT", start_notation),
        ),
        (
            "an exon end in the 5' UTR",
            row("ENST00000000003", 31099, "G", "TAT", start_notation),
        ),
        (
            "an exon end in the 3' UTR",
            row("ENST00000000003", 31499, "G", "TAT", start_notation),
        ),
        (
            "a transcript the cache lacks",
            row("ENST00000000009", 11099, "G", "TAT", start_notation),
        ),
    ];
    for (shape, row) in &shapes {
        assert_eq!(
            verdicts(Some(&cache), &row.texts(), &[INTRON_START]),
            difference(),
            "{shape}"
        );
    }
    for location in ["1:11097-11098", "1:11101-11102", "1:11098", "1:11101-11103"] {
        assert_eq!(
            verdicts(
                Some(&cache),
                &left().location(location).texts(),
                &[INTRON_START]
            ),
            placed(red(), green(INTRON_START)),
            "{location}"
        );
    }
    assert_eq!(
        verdicts(
            Some(&cache),
            &left().record_reference("CG").texts(),
            &[INTRON_START]
        ),
        placed(green(INTRON_START), red())
    );
}

both_ways!(
    the_shapes_nearest_the_intron_start_class_stay_differences,
    the_intron_start_rule_leaves_the_nearest_shapes,
    the_shapes_nearest_the_intron_start_class
);

/// A record whose reference allele spans two bases and that carries another allele beside the
/// insertion keeps both bases when its alleles are minimised: the text row's Location is their
/// two adjacent positions and its Allele the insertion alone. Either class places the text row
/// by each gap of the span and the entry by the ALT that clips against the REF to exactly the
/// entry's allele or, on a moved record, to that allele turned back by the offset. A Location
/// none of whose three gaps is a gap of the class stays a difference in every format, and so
/// does an entry whose allele no ALT clips to or whose ALT's other bases are not the REF's.
fn a_minimised_allele_of_a_two_base_reference(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let held = "CGTAT,C";
    let (end_notation, start_notation, moved_notation) =
        ("c.286-1_286insTAT", "c.50_50+1insTAT", "c.50_50+1insTTA");
    let end_rows = [
        row("ENST00000000004", 41312, "G", "TAT", end_notation)
            .record_reference("CG")
            .record_alternates(held),
        row("ENST00000000004", 41313, "G", "TAT", end_notation)
            .location("1:41314-41315")
            .record_alternates("GTAT,C"),
        row("ENST00000000002", 25298, "G", "ATA", "c.51-1_51insTAT")
            .strand(-1)
            .record_reference("CG")
            .record_alternates("CGATA,C"),
    ];
    for row in &end_rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_END]),
            set_aside(INTRON_END)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    let two_bases = |position: i64, allele: &'static str, notation: &'static str| {
        row("ENST00000000001", position, "G", allele, notation)
            .record_reference("CG")
            .record_alternates(held)
    };
    let start_rows = [
        two_bases(11098, "TAT", start_notation),
        row("ENST00000000001", 11099, "G", "TAT", start_notation)
            .location("1:11100-11101")
            .record_alternates("GTAT,C"),
        two_bases(11096, "TAT", moved_notation).both("HGVS_OFFSET", "2"),
        two_bases(11096, "TTA", moved_notation).both("HGVS_OFFSET", "2"),
    ];
    for row in &start_rows {
        let texts = row.texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[INTRON_START]),
            set_aside(INTRON_START)
        );
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    let shapes = [
        (
            "the span below the end gap",
            INTRON_END,
            row("ENST00000000004", 41310, "G", "TAT", end_notation),
        ),
        (
            "the span above the end gap",
            INTRON_END,
            row("ENST00000000004", 41315, "G", "TAT", end_notation),
        ),
        (
            "the span below the start gap",
            INTRON_START,
            row("ENST00000000001", 11096, "G", "TAT", start_notation),
        ),
        (
            "the span above the start gap",
            INTRON_START,
            row("ENST00000000001", 11101, "G", "TAT", start_notation),
        ),
    ];
    for (shape, klass, row) in shapes {
        let texts = row.record_reference("CG").record_alternates(held).texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[klass]),
            difference(),
            "{shape}"
        );
    }
    let entries = [
        (
            "no ALT clips to the allele, end",
            INTRON_END,
            row("ENST00000000004", 41312, "G", "TAT", end_notation),
            "CGCCG,C",
        ),
        (
            "the other bases are not the REF's",
            INTRON_END,
            row("ENST00000000004", 41312, "G", "TAT", end_notation),
            "GGTAT,C",
        ),
        (
            "no ALT clips to the allele, start",
            INTRON_START,
            row("ENST00000000001", 11098, "G", "TAT", start_notation),
            "CGCCG,C",
        ),
    ];
    for (shape, klass, row, alternates) in entries {
        let texts = row
            .record_reference("CG")
            .record_alternates(alternates)
            .texts();
        assert_eq!(
            verdicts(Some(&cache), &texts, &[klass]),
            placed(green(klass), red()),
            "{shape}"
        );
    }
    assert_eq!(turned_back("TTA", 2), "TAT");
    assert_eq!(turned_back("TAA", -2), "ATA");
    assert_eq!(turned_back("TAT", 0), "TAT");
    assert_eq!(turned_back("TAT", 3), "TAT");
    assert_eq!(turned_back("", 2), "");
}

both_ways!(
    a_minimised_allele_of_a_two_base_reference_is_placed_by_either_intron_class,
    the_intron_rules_place_a_minimised_allele_of_a_two_base_reference,
    a_minimised_allele_of_a_two_base_reference
);

/// A VCF or JSON entry is placed from its input line by the insertion its HGVSc describes,
/// however the entry spells its allele. A record whose alleles all share their first base
/// gives an entry whose allele is the ALT less that base; a minimised allele may equal another
/// ALT of its line. Either class places both by the ALT that clips against the REF to the
/// bases the notation inserts, or to as many bases as it duplicates, and the row is the class
/// in every format. An entry stays a difference where no ALT of its line clips to the
/// notation's bases, and so does a row whose ALT clips to them in a gap that is no gap of the
/// class.
fn an_entry_placed_by_the_bases_its_notation_inserts(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let (end_notation, start_notation) = ("c.286-1_286insTAT", "c.50_50+1insTAT");
    let shared =
        |row: intron_ends_common::Row| row.record_reference("GC").record_alternates("GTATC,G");
    let doubled =
        |row: intron_ends_common::Row| row.record_reference("GC").record_alternates("G,GGC");
    let rows = [
        (
            INTRON_END,
            shared(row("ENST00000000004", 41313, "G", "TATC", end_notation).location("1:41314")),
        ),
        (
            INTRON_END,
            doubled(
                row("ENST00000000004", 41313, "G", "G", "c.286-1dup").location("1:41313-41314"),
            ),
        ),
        (
            INTRON_START,
            shared(row("ENST00000000001", 11099, "G", "TATC", start_notation).location("1:11100")),
        ),
        (
            INTRON_START,
            doubled(
                row("ENST00000000001", 11099, "G", "G", "c.50_50+1insG").location("1:11099-11100"),
            ),
        ),
    ];
    for (klass, row) in &rows {
        let texts = row.texts();
        assert_eq!(verdicts(Some(&cache), &texts, &[*klass]), set_aside(klass));
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    let unclipped = [
        (
            INTRON_END,
            row("ENST00000000004", 41313, "G", "TAT", end_notation).location("1:41313"),
        ),
        (
            INTRON_START,
            row("ENST00000000001", 11099, "G", "TAT", start_notation).location("1:11099"),
        ),
    ];
    for (klass, row) in unclipped {
        assert_eq!(
            verdicts(
                Some(&cache),
                &row.record_alternates("GCCG,C").texts(),
                &[klass]
            ),
            placed(green(klass), red()),
            "{klass}"
        );
    }
    let beside = [
        (
            INTRON_END,
            row("ENST00000000004", 41311, "G", "TATC", end_notation).location("1:41312"),
        ),
        (
            INTRON_START,
            row("ENST00000000001", 11097, "G", "TATC", start_notation).location("1:11098"),
        ),
    ];
    for (klass, row) in beside {
        assert_eq!(
            verdicts(Some(&cache), &shared(row).texts(), &[klass]),
            difference(),
            "{klass}"
        );
    }
}

both_ways!(
    an_entry_is_placed_by_the_bases_its_notation_inserts,
    the_intron_rules_place_an_entry_by_the_bases_its_notation_inserts,
    an_entry_placed_by_the_bases_its_notation_inserts
);

/// An input line in the Ensembl default format gives its record's alleles as a VCF line does:
/// an insertion written alone, its start one above its end and its reference allele empty, and
/// one held with a shared base beside another allele, on the forward strand and, its alleles
/// reverse-complemented, on the reverse strand. Either class places both the JSON and the VCF
/// entry by the insertion the HGVSc describes. An entry stays a difference where no allele of
/// its line clips to the notation's bases, where the line's gap is no gap of the class and
/// where the line is of neither format; a VCF entry stays one where its line is written at
/// another position than the input line at its ordinal.
fn an_entry_of_an_ensembl_format_input_line(verdicts: Verdicts, name: &str) {
    let cache = intron_end_cache(&scratch(name));
    let (end_notation, start_notation) = ("c.286-1_286insTAT", "c.50_50+1insTAT");
    let padded =
        |row: intron_ends_common::Row| row.record_reference("CG").record_alternates("CGTAT,CC");
    let end = |position: i64, allele: &'static str| {
        row("ENST00000000004", position, "G", allele, end_notation)
    };
    let start = |position: i64, allele: &'static str| {
        row("ENST00000000001", position, "G", allele, start_notation)
    };
    let rows = [
        (
            INTRON_END,
            end(41313, "TAT").input_line("1\t41314\t41313\t-/TAT\t+\trec1"),
        ),
        (
            INTRON_END,
            end(41313, "TAT").input_line("1\t41314\t41313\t-/ATA\t-\trec1"),
        ),
        (
            INTRON_END,
            padded(
                end(41312, "GTAT")
                    .location("1:41313")
                    .input_line("1\t41313\t41313\tG/GTAT/C\t1\trec1"),
            ),
        ),
        (
            INTRON_END,
            padded(
                end(41312, "ATAC")
                    .location("1:41313")
                    .input_line("1\t41313\t41313\tC/ATAC/G\t-1\trec1"),
            ),
        ),
        (
            INTRON_START,
            start(11099, "TAT").input_line("1\t11100\t11099\t-/TAT\t+\trec1"),
        ),
        (
            INTRON_START,
            start(11099, "TAT").input_line("1\t11100\t11099\t-/ATA\t-\trec1"),
        ),
        (
            INTRON_START,
            padded(
                start(11098, "GTAT")
                    .location("1:11099")
                    .input_line("1\t11099\t11099\tG/GTAT/C\t1\trec1"),
            ),
        ),
        (
            INTRON_START,
            padded(
                start(11098, "ATAC")
                    .location("1:11099")
                    .input_line("1\t11099\t11099\tC/ATAC/G\t-1\trec1"),
            ),
        ),
    ];
    for (klass, row) in &rows {
        let texts = row.texts();
        assert_eq!(verdicts(Some(&cache), &texts, &[*klass]), set_aside(klass));
        assert_eq!(verdicts(Some(&cache), &texts, &[]), difference());
    }
    let entries = [
        (
            "no allele clips to the notation's bases",
            "1\t41314\t41313\t-/CCG\t+\trec1",
        ),
        ("a line of four fields", "1\t41314\t41313\t-/TAT"),
        (
            "a fifth field that names no strand",
            "1\t41314\t41313\t-/TAT\t0\trec1",
        ),
        (
            "an allele string of one allele",
            "1\t41314\t41313\tTAT\t+\trec1",
        ),
    ];
    for (shape, line) in entries {
        assert_eq!(
            verdicts(
                Some(&cache),
                &end(41313, "TAT").input_line(line).texts(),
                &[INTRON_END]
            ),
            placed(green(INTRON_END), red()),
            "{shape}"
        );
    }
    let beside = end(41311, "TAT").input_line("1\t41312\t41311\t-/TAT\t+\trec1");
    assert_eq!(
        verdicts(Some(&cache), &beside.texts(), &[INTRON_END]),
        difference()
    );
    for (position, line) in [
        (41300, "1\t41314\t41313\t-/TAT\t+\trec1"),
        (41313, "2\t41314\t41313\t-/TAT\t+\trec1"),
    ] {
        let row = end(position, "TAT")
            .location("1:41313-41314")
            .input_line(line);
        let verdicts = verdicts(Some(&cache), &row.texts(), &[INTRON_END]);
        assert_eq!(verdicts["vcf"], red(), "{line}");
        assert_eq!(text_rows(&verdicts), [&green(INTRON_END); 2]);
    }
}

both_ways!(
    an_entry_of_an_ensembl_format_input_line_is_placed_by_either_intron_class,
    the_intron_rules_place_an_entry_of_an_ensembl_format_input_line,
    an_entry_of_an_ensembl_format_input_line
);

#[test]
fn notation_insertion_gaps_names_the_gaps_an_input_line_gives_the_notation() {
    let gaps = |first: i64,
                reference: &str,
                alternates: &[&str],
                inserted: &str,
                strand: &str,
                duplicated: i64,
                moved: i64| {
        notation_insertion_gaps(
            first,
            reference,
            &lines(alternates),
            inserted,
            strand,
            duplicated,
            moved,
        )
    };
    assert_eq!(gaps(100, "G", &["GTAT", "C"], "TAT", "1", 0, 0), [100]);
    assert_eq!(gaps(100, "G", &["GATA"], "TAT", "-1", 0, 0), [100]);
    assert!(gaps(100, "G", &["GATA"], "TAT", "1", 0, 0).is_empty());
    assert_eq!(gaps(100, "G", &["GATA"], "TAT", "", 0, 0), [100]);
    assert_eq!(gaps(100, "G", &["GTAT"], "TAT", "", 0, 0), [100]);
    assert_eq!(gaps(100, "G", &["GTAT"], "TTA", "1", 0, 2), [100]);
    assert!(gaps(100, "G", &["GTAT"], "TTA", "1", 0, 0).is_empty());
    assert_eq!(gaps(100, "G", &["GATA"], "TTA", "-1", 0, -2), [100]);
    assert!(gaps(100, "G", &["GATT"], "TTA", "1", 0, 2).is_empty());
    assert_eq!(gaps(100, "GC", &["G", "GTATC"], "TAT", "1", 0, 0), [100]);
    assert_eq!(
        gaps(100, "CG", &["C", "CGTAT", "TATCG"], "TAT", "1", 0, 0),
        [101, 99]
    );
    assert_eq!(gaps(100, "GC", &["G", "GGC"], "", "1", 1, 0), [99, 100]);
    assert_eq!(gaps(100, "G", &["GCAG"], "", "1", 3, 0), [99, 100]);
    assert!(gaps(100, "G", &["GCAG"], "", "1", 1, 0).is_empty());
    assert!(gaps(100, "G", &["GCAG"], "", "1", 3, 1).is_empty());
    assert_eq!(gaps(100, "G", &["GCAG"], "", "1", 0, 0), [100]);
    assert_eq!(gaps(101, "", &["CAG"], "", "1", 0, 0), [100]);
    assert!(gaps(100, "CG", &["CGTAT"], "", "1", 0, 0).is_empty());
    assert!(gaps(100, "G", &["TATG"], "", "1", 0, 0).is_empty());
}

#[test]
fn intron_ends_holds_the_gaps_in_which_a_coding_exon_ends_before_an_intron() {
    let [forward, reverse, utr, short] = intron_end_transcripts();
    assert_eq!(transcript_intron_start_gaps(&forward), set(&[11099, 11299]));
    assert_eq!(transcript_intron_start_gaps(&reverse), set(&[25199, 25499]));
    assert_eq!(transcript_intron_start_gaps(&utr), set(&[31299]));
    assert_eq!(
        transcript_intron_start_gaps(&short),
        set(&[41099, 41199, 41299])
    );
    assert_eq!(
        transcript_intron_start_gaps(&as_text(&reverse)),
        set(&[25199, 25499])
    );
    assert!(transcript_intron_start_gaps(&without_coding(&forward)).is_empty());
    assert_eq!(coding_transcript_exons(&without_coding(&forward)), None);
    assert_eq!(
        coding_transcript_exons(&reverse),
        Some((
            false,
            vec![(25000, 25099), (25200, 25299), (25500, 25599)],
            25050,
            25549
        ))
    );
    let abutting = with_member(&forward, "cdna_coding_end", Value::Int(350));
    let abutting = with_member(
        &abutting,
        "exons",
        compare_common::json(
            r#"[{"start": 11000, "end": 11099}, {"start": 11100, "end": 11299}, {"start": 11301, "end": 11599}]"#,
        ),
    );
    assert_eq!(transcript_intron_start_gaps(&abutting), set(&[11299]));

    let root = scratch("intron_ends_holds_the_gaps_in_which_a_coding_exon_ends_before_an_intron");
    let cache = intron_end_cache(&root.join("cache"));
    let ends = IntronEnds::new(&cache, None).unwrap();
    assert_eq!(
        ends.start_gaps("1", "ENST00000000002", 25499).unwrap(),
        set(&[25199, 25499])
    );
    // One read of the shard serves both kinds of gap.
    assert!(unread(&ends, "1").is_empty());
    spoil(&cache.join("transcripts").join("1"), &["1-1000000.json"]);
    assert_eq!(
        ends.gaps("1", "ENST00000000002", 25299).unwrap(),
        set(&[25099, 25299])
    );
    assert_eq!(
        ends.transcript_gaps("chr1", "ENST00000000004.2", 41199)
            .unwrap(),
        (set(&[41313]), set(&[41099, 41199, 41299]))
    );
    assert!(ends
        .start_gaps("1", "ENST00000000009", 11099)
        .unwrap()
        .is_empty());
    let row = parse_default(&default_text(&[default_row(
        "rec1",
        "1:11097-11098",
        "T",
        "ENST00000000001",
        "frameshift_variant",
        "IMPACT=HIGH",
    )]))
    .unwrap()
    .remove(0);
    assert!(ends.holds_start(&row, 2, "").unwrap());
    assert!(!ends.holds_start(&row, 0, "").unwrap());
    assert!(!ends.holds_start(&row, -2, "").unwrap());
    assert!(!ends.holds(&row, "", 0).unwrap());
}

/// The comparison of two default-format texts under no class and no side table sees the
/// planted row: the helpers of these tests write what the parsers read.
#[test]
fn a_planted_row_differs_in_its_protein_notation_alone() {
    let (reference, candidate, _) =
        row("ENST00000000001", 11199, "G", "G", "c.51-1dup").text("default");
    let expected = parse_default(&reference).unwrap();
    let actual = parse_default(&candidate).unwrap();
    assert_eq!(expected.len(), 1);
    assert_eq!(
        expected[0].fields.get("HGVSp").map(String::as_str),
        Some(INTRON_END_PROTEIN)
    );
    assert_eq!(actual[0].fields.get("HGVSp"), None);
    let mut documented = documented_with(&[]);
    let report = compare(&expected, &actual, &mut documented);
    assert_eq!(compare_common::verdict(&report), "red");
    assert_eq!(counts(&documented.applied), count_map(&[]));
}
