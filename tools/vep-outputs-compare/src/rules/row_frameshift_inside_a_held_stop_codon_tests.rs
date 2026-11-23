// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row-aware pair rule of `frameshift_inside_a_held_stop_codon` on the fabricated rows of
//! `tests/fixtures/row_aware_pair_rules/tables/frameshift_inside_a_held_stop_codon_rows.tsv`:
//! the bases a Codons cell holds in its two halves, the rows the rule takes and those it
//! refuses, its place among the rules of the registry, and the pairs of the other row-aware
//! rules, which it leaves to them.

use super::super::{classify, excluded, row_pair_rules};
use super::*;

const CLASS: &str = "frameshift_inside_a_held_stop_codon";
const SETTLED_STOP: &str = "stop_codon_settled_before_unknown_bases";
const PAST_THE_SPAN: &str = "splice_term_past_the_record_span";
const GENETIC_CODE_26: &str = "pachysolen_nuclear_code_read_as_standard";
const REMOVED_START_CODON: &str = "start_codon_removed_without_five_prime_utr";
const KEPT_START_CODON: &str = "start_codon_kept_without_five_prime_utr";

/// The five classes of a coding end, the first row-aware rules of the registry: each takes
/// a pair only beside the reading of a transcript's cache record, and none where no cache is
/// read.
const CODING_END_CLASSES: [&str; 5] = [
    "structural_deletion_stop_lost_without_stop_codon",
    "stop_codon_completed_by_polyadenylation",
    "stop_codon_overlap_without_stop_codon",
    "stop_codon_behind_last_residue_without_stop_codon",
    "start_codon_overlap_without_start_codon",
];

const TABLE: &str = include_str!(
    "../../tests/fixtures/row_aware_pair_rules/tables/frameshift_inside_a_held_stop_codon_rows.tsv"
);

const COLUMNS: [&str; 13] = [
    "shape",
    "first_column",
    "location",
    "allele",
    "cdna",
    "cds",
    "reference_terms",
    "candidate_terms",
    "reference_amino_acids",
    "candidate_amino_acids",
    "reference_codons",
    "candidate_codons",
    "taken",
];

/// One fabricated Transcript row of the default format as a rule reads it: its first column,
/// its Location, its allele, its terms, its two position columns, its Amino_acids and its
/// Codons.
#[derive(Clone, Copy)]
struct Row<'a> {
    name: &'a str,
    location: &'a str,
    allele: &'a str,
    consequence: &'a str,
    cdna: &'a str,
    cds: &'a str,
    amino_acids: &'a str,
    codons: &'a str,
}

impl Row<'_> {
    fn columns(&self) -> [&str; 12] {
        [
            self.name,
            self.location,
            self.allele,
            "ENSG00000000001",
            "ENST00000000001",
            "Transcript",
            self.consequence,
            self.cdna,
            self.cds,
            "-",
            self.amino_acids,
            self.codons,
        ]
    }

    /// The row read to its Codons column.
    fn read(&self) -> PairRow {
        pair_row_with_codons(&self.columns())
    }

    /// The reference's row under the terms `reference` and vep-rs's under `candidate`, both
    /// with the row's own residues and codons.
    fn pair(&self, reference: &str, candidate: &str) -> (PairRow, PairRow) {
        let under = |consequence| Row {
            consequence,
            ..*self
        };
        (under(reference).read(), under(candidate).read())
    }
}

/// One base lost from inside the codon of a stop the peptide holds, under an identifier: the
/// reference residue is the stop, and the frame after the lost base reads an unknown residue.
const LOST_BASE_ROW: Row = Row {
    name: "a_base_lost_inside_the_stop_codon",
    location: "1:136",
    allele: "-",
    consequence: "",
    cdna: "36",
    cds: "30",
    amino_acids: "*/X",
    codons: "tgA/tg",
};

/// The terms of the class's pair on `LOST_BASE_ROW`: the reference's and vep-rs's.
const LOST_BASE_TERMS: (&str, &str) = ("coding_sequence_variant", "frameshift_variant");

/// The pairs of the table: a shape, the reference's row, vep-rs's row and whether the class
/// takes the pair.
fn table_pairs() -> Vec<(&'static str, PairRow, PairRow, bool)> {
    let mut lines = TABLE.lines();
    assert_eq!(lines.next(), Some(&COLUMNS.join("\t")[..]));
    let pair = |line: &'static str| {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), COLUMNS.len(), "{line}");
        let side = |consequence, amino_acids, codons| Row {
            name: cells[1],
            location: cells[2],
            allele: cells[3],
            consequence,
            cdna: cells[4],
            cds: cells[5],
            amino_acids,
            codons,
        };
        let taken = match cells[12] {
            "yes" => true,
            "no" => false,
            other => panic!("neither yes nor no: {other}"),
        };
        (
            cells[0],
            side(cells[6], cells[8], cells[10]).read(),
            side(cells[7], cells[9], cells[11]).read(),
            taken,
        )
    };
    lines.map(pair).collect()
}

/// The buckets of the row-aware rules of the registry that take the reference's row `e` and
/// vep-rs's row `a`, in the registry's order.
fn row_rules_that_take(e: &PairRow, a: &PairRow) -> Vec<&'static str> {
    row_pair_rules()
        .into_iter()
        .filter(|rule| rule.matches_rows(std::slice::from_ref(e), std::slice::from_ref(a)))
        .map(|rule| rule.bucket)
        .collect()
}

/// A half of `-` holds no base and any other half as many as its characters; a cell without
/// a `/` gives no lengths.
#[test]
fn a_codons_cell_gives_the_bases_of_its_two_halves() {
    for (cell, bases) in [
        ("tgA/tg", Some((3, 2))),
        ("tga/tTga", Some((3, 4))),
        ("tGAGCT/t", Some((6, 1))),
        ("tga/tgAATGa", Some((3, 7))),
        ("TGA/TTGA", Some((3, 4))),
        ("tgN/tg", Some((3, 2))),
        ("-/T", Some((0, 1))),
        ("tga/-", Some((3, 0))),
        ("-/-", Some((0, 0))),
        ("tga", None),
        ("-", None),
        ("", None),
    ] {
        assert_eq!(codon_half_bases(cell), bases, "{cell:?}");
    }
}

/// The class takes vep-rs's `frameshift_variant` on a row whose reference residues open
/// with a stop: the reference's row holds `coding_sequence_variant` alone,
/// `inframe_deletion` alone, or `inframe_insertion` with `coding_sequence_variant` in its
/// place, every other term stands on both rows, the two rows hold one Amino_acids cell whose
/// alternate half holds the unknown residue, and one Codons cell whose two halves differ in
/// length by a number that is no multiple of three. A reference row that holds
/// `stop_retained_variant` or `stop_lost` is no row of the class.
#[test]
fn the_class_takes_a_frameshift_inside_a_stop_codon_the_peptide_holds_and_no_other_row() {
    let pairs = table_pairs();
    let wrong: Vec<&str> = pairs
        .iter()
        .filter(|(_, e, a, taken)| frameshift_inside_a_held_stop_codon_pair(e, a) != *taken)
        .map(|(shape, ..)| *shape)
        .collect();
    assert_eq!(wrong, Vec::<&str>::new());
    let taken = pairs.iter().filter(|(.., taken)| *taken).count();
    assert_eq!((taken, pairs.len()), (10, 42));
    let spelt = pairs.iter().filter(|(_, e, ..)| e.indel.is_some());
    assert_eq!(
        spelt.count(),
        1,
        "the row whose first column spells the alleles"
    );
}

/// The two rows the engines write alike are no pair: the same base added one position
/// earlier, which both call `frameshift_variant`, and a base lost inside the closing stop
/// codon, which both call `stop_lost`.
#[test]
fn rows_both_engines_write_alike_are_no_pair_of_the_class() {
    let alike: Vec<(&str, String)> = table_pairs()
        .into_iter()
        .filter(|(_, e, a, _)| e == a)
        .map(|(shape, e, a, taken)| {
            assert!(!taken, "{shape}");
            assert!(!frameshift_inside_a_held_stop_codon_pair(&e, &a), "{shape}");
            (shape, e.consequence)
        })
        .collect();
    assert_eq!(
        alike,
        [
            (
                "the same base added, written before the stop codon",
                "frameshift_variant".to_string()
            ),
            (
                "a base lost inside the closing stop codon",
                "stop_lost".to_string()
            ),
        ]
    );
}

/// No rule over the two consequence sets takes a pair of the table, so its rows reach the
/// row-aware rules; there the class's rule stands twelfth, after the five rules of a coding
/// end and the six other stop codon rules and before the four rules that follow them. It
/// classifies and sets aside, it is the
/// one rule that takes a pair of the class, and a pair it refuses is taken by no rule.
#[test]
fn a_pair_of_the_class_passes_every_other_rule_and_is_set_aside_under_its_own_bucket() {
    let rules = row_pair_rules();
    let buckets: Vec<&str> = rules.iter().map(|rule| rule.bucket).collect();
    let at = buckets
        .iter()
        .position(|bucket| *bucket == CLASS)
        .expect("the class's rule");
    assert!(rules[at].classifies && rules[at].excludes);
    assert_eq!(
        buckets[..at].last(),
        Some(&SETTLED_STOP),
        "the last of the eleven rules the registry asks before the class's"
    );
    assert_eq!(
        buckets[at + 1..],
        [
            PAST_THE_SPAN,
            GENETIC_CODE_26,
            REMOVED_START_CODON,
            KEPT_START_CODON
        ],
        "the four rules the registry asks after the class's"
    );
    assert_eq!((at + 1, buckets.len()), (12, 16));
    for (shape, e, a, taken) in table_pairs() {
        assert_eq!(
            classify(&e.consequence, &a.consequence, &e.allele),
            None,
            "{shape}"
        );
        assert!(
            !excluded(&e.consequence, &a.consequence, &e.allele),
            "{shape}"
        );
        let expected: &[&str] = if taken { &[CLASS] } else { &[] };
        assert_eq!(row_rules_that_take(&e, &a), expected, "{shape}");
    }
}

/// The same pair read to its ninth column and to its twelfth: a row read to its ninth
/// carries neither Amino_acids nor Codons, and the class takes the pair where both rows
/// carry theirs.
#[test]
fn a_row_read_without_its_codons_is_no_row_of_the_class() {
    let (e, a) = LOST_BASE_ROW.pair(LOST_BASE_TERMS.0, LOST_BASE_TERMS.1);
    assert!(frameshift_inside_a_held_stop_codon_pair(&e, &a));
    let nine = |consequence| {
        let row = Row {
            consequence,
            ..LOST_BASE_ROW
        };
        let columns = row.columns();
        let first: &[&str; 9] = columns.first_chunk().expect("nine of twelve columns");
        pair_row(first)
    };
    let (short_e, short_a) = (nine(LOST_BASE_TERMS.0), nine(LOST_BASE_TERMS.1));
    assert_eq!(
        (short_e.amino_acids.as_deref(), short_e.codons.as_deref()),
        (None, None)
    );
    for (reference, candidate) in [(&short_e, &short_a), (&short_e, &a), (&e, &short_a)] {
        assert!(!frameshift_inside_a_held_stop_codon_pair(
            reference, candidate
        ));
        assert_eq!(
            row_rules_that_take(reference, candidate),
            Vec::<&str>::new()
        );
    }
}

/// The rule reads the terms, the Amino_acids and the Codons of the two rows and no other
/// column: the class's pair is taken under another Location, another Allele and other
/// positions on either row.
#[test]
fn the_class_reads_no_column_but_the_terms_the_residues_and_the_codons() {
    let (e, a) = LOST_BASE_ROW.pair(LOST_BASE_TERMS.0, LOST_BASE_TERMS.1);
    let elsewhere = |row: &PairRow| PairRow {
        location: "7:100998-101000".to_string(),
        allele: "GTNCA".to_string(),
        cdna_position: "?-7057".to_string(),
        cds_position: "-".to_string(),
        ..row.clone()
    };
    for (reference, candidate) in [
        (elsewhere(&e), a.clone()),
        (e.clone(), elsewhere(&a)),
        (elsewhere(&e), elsewhere(&a)),
    ] {
        assert!(frameshift_inside_a_held_stop_codon_pair(
            &reference, &candidate
        ));
    }
}

/// A pair of each other row-aware rule that takes one where no cache is read
/// (`CODING_END_CLASSES` take none), with this class's residues and codons on both rows
/// where that rule reads neither: the five stop codon rules differ by
/// `stop_retained_variant`, the settled stop's by `stop_gained` on a reference row without
/// residues, the past-the-span rule's by a splice term, the genetic code 26 rule's by
/// `synonymous_variant` against `missense_variant`, the one after it by `start_lost` beside an
/// in-frame deletion from the first coding base and the last by a start term. Each is taken
/// by its own rule alone, and none is a pair of this class.
#[test]
fn a_pair_of_another_row_aware_rule_is_left_to_that_rule() {
    let deletion = Row {
        name: "1993671",
        ..LOST_BASE_ROW
    };
    let inserted_at_a_stop_codon = Row {
        name: "21_35201927_G/GAAN",
        location: "21:35201927-35201928",
        allele: "AAN",
        cdna: "3048-3049",
        cds: "3047-3048",
        ..LOST_BASE_ROW
    };
    let (unknown, kept) = (
        "splice_region_variant,coding_sequence_variant",
        "splice_region_variant,stop_retained_variant",
    );
    let settled = |consequence, amino_acids, codons| {
        Row {
            name: "1993671",
            location: "21:10567700-10567701",
            allele: "GAATAGNNNN",
            consequence,
            cdna: "775-776",
            cds: "423-424",
            amino_acids,
            codons,
        }
        .read()
    };
    let on_code_26 = |consequence, amino_acids| {
        Row {
            name: "1993671",
            location: "7:1000",
            allele: "A",
            consequence,
            cdna: "114",
            cds: "114",
            amino_acids,
            codons: "ctG/ctA",
        }
        .read()
    };
    let pairs = [
        (
            "stop_codon_overlap_genomic_span",
            Row {
                location: "3:36993663",
                cdna: "1705",
                cds: "1630",
                ..deletion
            }
            .pair(unknown, kept),
        ),
        (
            "stop_retained_without_stop_codon_test",
            inserted_at_a_stop_codon.pair("stop_retained_variant", "frameshift_variant"),
        ),
        (
            "stop_retained_three_base_stop_codon_test",
            Row {
                location: "21:2101-2107",
                cdna: "961-967",
                cds: "264-270",
                ..deletion
            }
            .pair("stop_retained_variant", "frameshift_variant"),
        ),
        (
            "stop_retained_three_base_stop_codon_test",
            Row {
                location: "21:2101-2107",
                cdna: "961-967",
                cds: "264-270",
                ..deletion
            }
            .pair(
                "stop_retained_variant,inframe_deletion",
                "frameshift_variant",
            ),
        ),
        (
            "stop_codon_removed_without_stop_codon_test",
            Row {
                location: "19:13208046-13208050",
                cdna: "?-7057",
                cds: "-",
                ..deletion
            }
            .pair("stop_retained_variant", "stop_lost"),
        ),
        (
            "stop_codon_overlap_insertion_range",
            inserted_at_a_stop_codon.pair(unknown, kept),
        ),
        (
            "stop_codon_overlap_insertion_range",
            Row {
                name: "1993671",
                allele: "TTGAN",
                ..inserted_at_a_stop_codon
            }
            .pair("frameshift_variant", "stop_retained_variant"),
        ),
        (
            SETTLED_STOP,
            (
                settled("frameshift_variant", "-", "-/GAATAGNNNN"),
                settled("stop_gained,frameshift_variant", "-/E*XX", "-/GAATAGNNNN"),
            ),
        ),
        (
            SETTLED_STOP,
            (
                settled("stop_retained_variant", "-", "-/GAATAGNNNN"),
                settled("stop_gained,frameshift_variant", "-/E*XX", "-/GAATAGNNNN"),
            ),
        ),
        (
            PAST_THE_SPAN,
            Row {
                name: "1993671",
                location: "7:100998-101000",
                allele: "GTNCA",
                cdna: "998-1000",
                cds: "898-900",
                ..LOST_BASE_ROW
            }
            .pair(
                "splice_donor_variant,frameshift_variant,intron_variant",
                "frameshift_variant,splice_region_variant",
            ),
        ),
        (
            GENETIC_CODE_26,
            (
                on_code_26("synonymous_variant", "L"),
                on_code_26("missense_variant", "A/L"),
            ),
        ),
        (
            REMOVED_START_CODON,
            Row {
                name: "21_44817965_ACAT/A",
                location: "21:44817966-44817968",
                allele: "-",
                consequence: "",
                cdna: "1-3",
                cds: "1-3",
                amino_acids: "M/-",
                codons: "ATG/-",
            }
            .pair("inframe_deletion", "start_lost,inframe_deletion"),
        ),
        (
            KEPT_START_CODON,
            Row {
                location: "21:30001-30010",
                cdna: "3-12",
                cds: "3-12",
                ..deletion
            }
            .pair(
                "frameshift_variant,start_lost",
                "frameshift_variant,start_retained_variant",
            ),
        ),
    ];
    let others: Vec<&str> = row_pair_rules()
        .iter()
        .map(|rule| rule.bucket)
        .filter(|bucket| *bucket != CLASS && !CODING_END_CLASSES.contains(bucket))
        .collect();
    let mut shown: Vec<&str> = pairs.iter().map(|(bucket, _)| *bucket).collect();
    shown.dedup();
    assert_eq!(
        shown, others,
        "a pair of every other row-aware rule that takes one where no cache is read"
    );
    for (bucket, (e, a)) in pairs {
        let label = format!("{bucket}: {} -> {}", e.consequence, a.consequence);
        assert_eq!(row_rules_that_take(&e, &a), [bucket], "{label}");
        assert!(!frameshift_inside_a_held_stop_codon_pair(&e, &a), "{label}");
    }
}
