// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row-aware pair rule of `stop_codon_settled_before_unknown_bases` on the fabricated rows
//! of `tests/fixtures/row_aware_pair_rules/tables/settled_stop_codon_rows.tsv`: the stop codon
//! a Codons cell settles, the rows the rule takes and those it refuses, its place among the
//! rules of the registry, the pairs of the other row-aware rules, which it leaves to them,
//! and the one pair of its own that an earlier rule takes first.

use super::super::{classify, excluded, row_pair_rules};
use super::*;

const CLASS: &str = "stop_codon_settled_before_unknown_bases";
const KEPT_INSERTION: &str = "stop_codon_overlap_insertion_range";
const HELD_STOP_CODON: &str = "frameshift_inside_a_held_stop_codon";
const PAST_THE_SPAN: &str = "splice_term_past_the_record_span";
const GENETIC_CODE_26: &str = "pachysolen_nuclear_code_read_as_standard";
const REMOVED_START_CODON: &str = "start_codon_removed_without_five_prime_utr";
const KEPT_START_CODON: &str = "start_codon_kept_without_five_prime_utr";

const TABLE: &str =
    include_str!("../../tests/fixtures/row_aware_pair_rules/tables/settled_stop_codon_rows.tsv");

const COLUMNS: [&str; 15] = [
    "shape",
    "first_column",
    "location",
    "allele",
    "feature_type",
    "cdna",
    "cds",
    "candidate_cds",
    "reference_terms",
    "candidate_terms",
    "reference_amino_acids",
    "candidate_amino_acids",
    "reference_codons",
    "candidate_codons",
    "taken",
];

/// One fabricated row of the default format as a rule reads it: its first column, its
/// Location, its allele, its terms, its two position columns, its Amino_acids and its Codons.
#[derive(Clone, Copy)]
struct Row<'a> {
    name: &'a str,
    location: &'a str,
    allele: &'a str,
    feature_type: &'a str,
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
            "ENSG00000274391",
            "ENST00000622113",
            self.feature_type,
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

    /// The reference's row under the terms `reference`, without residues, and vep-rs's
    /// under `candidate`, with the row's own.
    fn pair(&self, reference: &str, candidate: &str) -> (PairRow, PairRow) {
        let reference = Row {
            consequence: reference,
            amino_acids: "-",
            ..*self
        };
        let candidate = Row {
            consequence: candidate,
            ..*self
        };
        (reference.read(), candidate.read())
    }
}

/// Ten bases between two neighbouring positions under an identifier: a sense codon, a stop
/// codon and four unknown bases, with the residues vep-rs writes for them.
const FRAMESHIFT_ROW: Row = Row {
    name: "1993671",
    location: "21:10567700-10567701",
    allele: "GAATAGNNNN",
    feature_type: "Transcript",
    consequence: "",
    cdna: "775-776",
    cds: "423-424",
    amino_acids: "-/E*XX",
    codons: "-/GAATAGNNNN",
};

/// The pairs of the table: a shape, the reference's row, vep-rs's row and whether the class
/// takes the pair.
fn table_pairs() -> Vec<(&'static str, PairRow, PairRow, bool)> {
    let mut lines = TABLE.lines();
    assert_eq!(lines.next(), Some(&COLUMNS.join("\t")[..]));
    let pair = |line: &'static str| {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), COLUMNS.len(), "{line}");
        let side = |consequence, cds, amino_acids, codons| Row {
            name: cells[1],
            location: cells[2],
            allele: cells[3],
            feature_type: cells[4],
            consequence,
            cdna: cells[5],
            cds,
            amino_acids,
            codons,
        };
        let taken = match cells[14] {
            "yes" => true,
            "no" => false,
            other => panic!("neither yes nor no: {other}"),
        };
        (
            cells[0],
            side(cells[8], cells[6], cells[10], cells[12]).read(),
            side(cells[9], cells[7], cells[11], cells[13]).read(),
            taken,
        )
    };
    lines.map(pair).collect()
}

/// The bucket of the first row-aware rule of the registry that takes the reference's row `e`
/// and vep-rs's row `a`, in the registry's order.
fn first_row_rule(e: &PairRow, a: &PairRow) -> Option<&'static str> {
    row_pair_rules()
        .into_iter()
        .find(|rule| rule.matches_rows(std::slice::from_ref(e), std::slice::from_ref(a)))
        .map(|rule| rule.bucket)
}

/// The stop codon is the first of the standard genetic code among the whole codons of the
/// cell's alternate half, read in upper case, before the first codon that holds a letter
/// other than `A`, `C`, `G` and `T`.
#[test]
fn a_codons_cell_settles_the_first_stop_codon_before_its_first_unknown_codon() {
    for (cell, place) in [
        ("-/GAATAGNNNN", Some(1)),
        ("-/GAATAGNNNNCA", Some(1)),
        ("-/TAGN", Some(0)),
        ("-/TGANNN", Some(0)),
        ("-/gaataaRca", Some(1)),
        ("tgg/tAGNNNNgg", Some(0)),
        ("-/GAATAGTAAN", Some(1)),
        ("-/GAAGCTTGAY", Some(2)),
        ("-/GAATAGAAAACA", None),
        ("-/GAAGCTNNNN", None),
        ("-/GTAGNNNNCA", None),
        ("-/GAANNNTAG", None),
        ("-/NNNTAGGAA", None),
        ("-/GAATANGGA", None),
        ("-/TAG/NNNN", None),
        ("-/GAATAG-NNN", None),
        ("GAATAGNNNN", None),
        ("-/", None),
        ("-", None),
        ("", None),
    ] {
        assert_eq!(settled_stop_codon(cell), place, "{cell:?}");
    }
}

/// The class takes a Transcript row whose Codons cell, the same on both rows, settles a stop
/// codon before its unknown bases: the reference's row writes no residues, vep-rs's residues
/// hold their first stop at that codon's place, and the terms differ by a gained stop or by
/// a kept one. A gained stop is vep-rs's `stop_gained` beside `frameshift_variant` on both
/// rows or in place of the reference's `coding_sequence_variant`, with at most one of the
/// three terms of a length of whole codons. Where the reference's row alone holds
/// `stop_retained_variant` on an insertion or deletion and both rows show that the stop
/// codon test is not made, vep-rs's alone holds the gained stop with `frameshift_variant`,
/// one of the three terms or no other term. A kept stop is vep-rs's `inframe_insertion`
/// beside `stop_retained_variant` on both rows, or the two terms in place of the
/// reference's `coding_sequence_variant` on an insertion between two neighbouring positions
/// that carries both cDNA positions and a CDS start. The Allele and the CDS position are
/// the same text on both rows, and the position carries both ends.
#[test]
fn the_class_takes_a_stop_the_codons_settle_before_their_unknown_bases_and_no_other_row() {
    let pairs = table_pairs();
    let wrong: Vec<&str> = pairs
        .iter()
        .filter(|(_, e, a, taken)| stop_codon_settled_before_unknown_bases_pair(e, a) != *taken)
        .map(|(shape, ..)| *shape)
        .collect();
    assert_eq!(wrong, Vec::<&str>::new());
    let taken = pairs.iter().filter(|(.., taken)| *taken).count();
    assert_eq!((taken, pairs.len()), (28, 83));
    let kept = pairs
        .iter()
        .filter(|(_, _, a, taken)| *taken && a.consequence.contains("stop_retained_variant"));
    assert_eq!(kept.count(), 9, "the kept stops the class takes");
    let spelt = pairs.iter().filter(|(_, e, ..)| e.indel.is_some());
    assert_eq!(
        spelt.count(),
        5,
        "the rows whose first column spells the alleles"
    );
}

/// No rule over the two consequence sets takes a pair of the table, so its rows reach the
/// row-aware rules; there the class's rule is the first to take a pair of the class, it
/// classifies and sets aside, and a pair it refuses is taken by no rule.
#[test]
fn a_pair_of_the_class_passes_every_earlier_rule_and_is_set_aside_under_its_own_bucket() {
    let rules = row_pair_rules();
    let buckets: Vec<&str> = rules.iter().map(|rule| rule.bucket).collect();
    let at = buckets
        .iter()
        .position(|bucket| *bucket == CLASS)
        .expect("the class's rule");
    assert!(rules[at].classifies && rules[at].excludes);
    assert_eq!(
        buckets[at + 1..],
        [
            HELD_STOP_CODON,
            PAST_THE_SPAN,
            GENETIC_CODE_26,
            REMOVED_START_CODON,
            KEPT_START_CODON
        ],
        "the five rules the registry asks after the class's"
    );
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
        assert_eq!(first_row_rule(&e, &a), taken.then_some(CLASS), "{shape}");
    }
}

/// The two rows of a pair carry one allele: the reference's row of one allele against
/// vep-rs's row of another is no pair of the class.
#[test]
fn rows_of_two_alleles_are_no_pair_of_the_class() {
    let (e, a) = FRAMESHIFT_ROW.pair("frameshift_variant", "stop_gained,frameshift_variant");
    assert!(stop_codon_settled_before_unknown_bases_pair(&e, &a));
    let other_allele = PairRow {
        allele: "GAATAGNNNNCA".to_string(),
        ..a.clone()
    };
    assert!(!stop_codon_settled_before_unknown_bases_pair(
        &e,
        &other_allele
    ));
    assert!(!stop_codon_settled_before_unknown_bases_pair(
        &other_allele,
        &a
    ));
}

/// The same pair read to its ninth column and to its twelfth: a row read to its ninth
/// carries no Codons, and the class takes the pair where both rows carry theirs.
#[test]
fn a_row_read_without_its_codons_is_no_row_of_the_class() {
    let terms = ("frameshift_variant", "stop_gained,frameshift_variant");
    let (e, a) = FRAMESHIFT_ROW.pair(terms.0, terms.1);
    assert!(stop_codon_settled_before_unknown_bases_pair(&e, &a));
    let nine = |consequence| {
        let row = Row {
            consequence,
            ..FRAMESHIFT_ROW
        };
        let columns = row.columns();
        let first: &[&str; 9] = columns.first_chunk().expect("nine of twelve columns");
        pair_row(first)
    };
    let (short_e, short_a) = (nine(terms.0), nine(terms.1));
    assert_eq!(
        (short_e.amino_acids.as_deref(), short_e.codons.as_deref()),
        (None, None)
    );
    for (reference, candidate) in [(&short_e, &short_a), (&short_e, &a), (&e, &short_a)] {
        assert!(!stop_codon_settled_before_unknown_bases_pair(
            reference, candidate
        ));
        assert_eq!(first_row_rule(reference, candidate), None);
    }
}

/// A pair of each other row-aware rule, with this class's Codons and residues on its rows
/// and an allele that holds an unknown base where that rule takes one: the five stop codon
/// rules differ by `stop_retained_variant`, the rule after this one by `frameshift_variant`
/// on rows whose reference residues open with a stop, the next by a splice term, the next by
/// `synonymous_variant` against `missense_variant`, the one after it by `start_lost` and the
/// last by a start term, and none of these pairs is one of this class.
#[test]
fn a_pair_of_another_row_aware_rule_is_left_to_that_rule() {
    let deletion = Row {
        allele: "-",
        ..FRAMESHIFT_ROW
    };
    let inserted_at_a_stop_codon = Row {
        name: "21_35201927_G/GAAN",
        location: "21:35201927-35201928",
        allele: "AAN",
        cdna: "3048-3049",
        cds: "3047-3048",
        ..FRAMESHIFT_ROW
    };
    let (unknown, kept) = (
        "splice_region_variant,coding_sequence_variant",
        "splice_region_variant,stop_retained_variant",
    );
    let on_code_26 = |consequence, amino_acids| {
        Row {
            location: "7:1000",
            allele: "A",
            consequence,
            cdna: "114",
            cds: "114",
            amino_acids,
            codons: "ctG/ctA",
            ..FRAMESHIFT_ROW
        }
        .read()
    };
    let at_the_first_codon = |consequence| {
        Row {
            location: "21:44817966-44817968",
            consequence,
            cdna: "1-3",
            cds: "1-3",
            amino_acids: "M/-",
            codons: "ATG/-",
            ..deletion
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
            inserted_at_a_stop_codon.pair(kept, unknown),
        ),
        (
            "stop_retained_three_base_stop_codon_test",
            Row {
                location: "21:2101-2107",
                cdna: "961-967",
                cds: "264-270",
                ..deletion
            }
            .pair("stop_retained_variant", "stop_gained,frameshift_variant"),
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
            .pair("stop_gained,frameshift_variant", "stop_retained_variant"),
        ),
        (
            HELD_STOP_CODON,
            (
                Row {
                    consequence: "coding_sequence_variant",
                    amino_acids: "*/X",
                    codons: "tgA/tg",
                    ..deletion
                }
                .read(),
                Row {
                    consequence: "frameshift_variant",
                    amino_acids: "*/X",
                    codons: "tgA/tg",
                    ..deletion
                }
                .read(),
            ),
        ),
        (
            PAST_THE_SPAN,
            Row {
                location: "7:100998-101000",
                allele: "GTNCA",
                cdna: "998-1000",
                cds: "898-900",
                ..FRAMESHIFT_ROW
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
            (
                at_the_first_codon("inframe_deletion"),
                at_the_first_codon("start_lost,inframe_deletion"),
            ),
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
    for (bucket, (e, a)) in pairs {
        let label = format!("{bucket} {}", e.allele);
        assert_eq!(first_row_rule(&e, &a), Some(bucket), "{label}");
        assert!(
            !stop_codon_settled_before_unknown_bases_pair(&e, &a),
            "{label}"
        );
    }
}

/// Seven bases between the last sense codon and the stop codon of a forward-strand
/// transcript, under an identifier: a stop codon and four unknown bases, with the residues
/// vep-rs writes for them.
const KEPT_STOP_ROW: Row = Row {
    name: "1993671",
    location: "21:10605549-10605550",
    allele: "TAGNNNN",
    feature_type: "Transcript",
    consequence: "",
    cdna: "1951-1952",
    cds: "1599-1600",
    amino_acids: "-/*XX",
    codons: "-/TAGNNNN",
};

/// vep-rs's terms on a row whose inserted bases keep the stop codon.
const KEPT: &str = "inframe_insertion,stop_retained_variant";

/// A kept stop in place of the reference's `frameshift_variant`, on an insertion between two
/// neighbouring positions that carries both cDNA positions and a CDS start, is a pair of the
/// class and of `stop_codon_overlap_insertion_range` as well: that rule reads the two sets
/// of terms without the residues and stands before this one, so the registry counts the
/// pair under its bucket.
#[test]
fn a_kept_stop_against_a_frameshift_is_a_pair_of_the_class_an_earlier_rule_takes_first() {
    let shifted = "frameshift_variant";
    let beside_a_shared_term = |terms: &str| format!("{terms},NMD_transcript_variant");
    let rows = [
        KEPT_STOP_ROW,
        Row {
            name: "21_10605549_T/TTAGNNNN",
            ..KEPT_STOP_ROW
        },
        Row {
            cdna: "1951-1952/2305",
            cds: "1599-1600/1605",
            ..KEPT_STOP_ROW
        },
        // Eight bases inside the last sense codon, which the alternate codons spell again
        // before the stop codon.
        Row {
            location: "21:10605548-10605549",
            allele: "TTAGNNNN",
            cdna: "1950-1951",
            cds: "1598-1599",
            amino_acids: "D/D*XX",
            codons: "gat/gaTTAGNNNNt",
            ..KEPT_STOP_ROW
        },
    ];
    for row in rows {
        for (reference, candidate) in [
            (shifted.to_string(), KEPT.to_string()),
            (beside_a_shared_term(shifted), beside_a_shared_term(KEPT)),
        ] {
            let label = format!("{} {} {reference}", row.name, row.location);
            let (e, a) = row.pair(&reference, &candidate);
            assert!(
                stop_codon_settled_before_unknown_bases_pair(&e, &a),
                "{label}"
            );
            assert_eq!(first_row_rule(&e, &a), Some(KEPT_INSERTION), "{label}");
            assert!(
                !stop_codon_settled_before_unknown_bases_pair(&a, &e),
                "{label}"
            );
        }
    }
}

/// A kept stop against the reference's `frameshift_variant` off the class's shape is no pair
/// of the class: the earlier rule takes the pairs that differ from the shape in their terms
/// or their residues alone, and no rule those on a row that is no insertion between two
/// neighbouring positions with the stop codon test's coordinates.
#[test]
fn a_kept_stop_against_a_frameshift_off_the_shape_is_no_pair_of_the_class() {
    let (e, a) = KEPT_STOP_ROW.pair("frameshift_variant", KEPT);
    assert!(stop_codon_settled_before_unknown_bases_pair(&e, &a));
    let off = [
        (
            "a second term the reference alone holds",
            KEPT_STOP_ROW,
            ("stop_gained,frameshift_variant", KEPT),
            Some(KEPT_INSERTION),
        ),
        (
            "no inframe_insertion on vep-rs's row",
            KEPT_STOP_ROW,
            ("frameshift_variant", "stop_retained_variant"),
            Some(KEPT_INSERTION),
        ),
        (
            "no inframe_insertion in place of coding_sequence_variant",
            KEPT_STOP_ROW,
            ("coding_sequence_variant", "stop_retained_variant"),
            Some(KEPT_INSERTION),
        ),
        (
            "residues with their stop one place later",
            Row {
                amino_acids: "-/X*X",
                ..KEPT_STOP_ROW
            },
            ("frameshift_variant", KEPT),
            Some(KEPT_INSERTION),
        ),
        (
            "a Location of three positions",
            Row {
                location: "21:10605549-10605551",
                ..KEPT_STOP_ROW
            },
            ("frameshift_variant", KEPT),
            None,
        ),
        (
            "no cDNA end",
            Row {
                cdna: "1951-?",
                ..KEPT_STOP_ROW
            },
            ("frameshift_variant", KEPT),
            None,
        ),
        (
            "a record whose reference allele keeps a base",
            Row {
                name: "21_10605549_TA/TTAGNNNN",
                ..KEPT_STOP_ROW
            },
            ("frameshift_variant", KEPT),
            None,
        ),
    ];
    for (shape, row, (reference, candidate), bucket) in off {
        let (e, a) = row.pair(reference, candidate);
        assert!(
            !stop_codon_settled_before_unknown_bases_pair(&e, &a),
            "{shape}"
        );
        assert_eq!(first_row_rule(&e, &a), bucket, "{shape}");
    }
}
