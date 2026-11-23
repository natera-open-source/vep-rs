// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row-aware pair rule of `start_codon_removed_without_five_prime_utr` on the
//! fabricated rows of `tests/fixtures/row_aware_pair_rules/tables/removed_start_codon_rows.tsv`:
//! the rows the rule takes and those it refuses, its place among the rules of the registry,
//! the columns a row is read to, and the pairs of the other row-aware rules, which it leaves
//! to them.

use super::super::{classify, excluded, row_pair_rules, ExclusionRule};
use super::*;

const CLASS: &str = "start_codon_removed_without_five_prime_utr";

const TABLE: &str =
    include_str!("../../tests/fixtures/row_aware_pair_rules/tables/removed_start_codon_rows.tsv");

const COLUMNS: [&str; 15] = [
    "shape",
    "first_column",
    "location",
    "allele",
    "feature_type",
    "reference_cdna",
    "candidate_cdna",
    "reference_cds",
    "candidate_cds",
    "reference_terms",
    "candidate_terms",
    "reference_amino_acids",
    "candidate_amino_acids",
    "codons",
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
            "ENSG00000184900",
            "ENST00000397893",
            self.feature_type,
            self.consequence,
            self.cdna,
            self.cds,
            "1",
            self.amino_acids,
            self.codons,
        ]
    }

    /// The row read to its Codons column.
    fn read(&self) -> PairRow {
        pair_row_with_codons(&self.columns())
    }

    /// The reference's row under the terms `reference` and vep-rs's under `candidate`.
    fn pair(&self, reference: &str, candidate: &str) -> (PairRow, PairRow) {
        let side = |consequence| Row {
            consequence,
            ..*self
        };
        (side(reference).read(), side(candidate).read())
    }
}

/// The three bases of the first codon of a transcript with no 5' UTR, removed by a record
/// whose alleles the first column spells.
const FIRST_CODON_ROW: Row = Row {
    name: "21_44817965_ACAT/A",
    location: "21:44817966-44817968",
    allele: "-",
    feature_type: "Transcript",
    consequence: "",
    cdna: "1-3",
    cds: "1-3",
    amino_acids: "M/-",
    codons: "ATG/-",
};

/// The terms of the reference's row of the class and of vep-rs's.
const TERMS: (&str, &str) = ("inframe_deletion", "start_lost,inframe_deletion");

/// The pairs of the table: a shape, the reference's row, vep-rs's row and whether the class
/// takes the pair.
fn table_pairs() -> Vec<(&'static str, PairRow, PairRow, bool)> {
    let mut lines = TABLE.lines();
    assert_eq!(lines.next(), Some(&COLUMNS.join("\t")[..]));
    let pair = |line: &'static str| {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), COLUMNS.len(), "{line}");
        let side = |cdna, cds, consequence, amino_acids| Row {
            name: cells[1],
            location: cells[2],
            allele: cells[3],
            feature_type: cells[4],
            consequence,
            cdna,
            cds,
            amino_acids,
            codons: cells[13],
        };
        let taken = match cells[14] {
            "yes" => true,
            "no" => false,
            other => panic!("neither yes nor no: {other}"),
        };
        (
            cells[0],
            side(cells[5], cells[7], cells[9], cells[11]).read(),
            side(cells[6], cells[8], cells[10], cells[12]).read(),
            taken,
        )
    };
    lines.map(pair).collect()
}

/// The rule the registry holds under the class's name.
fn class_rule() -> &'static ExclusionRule {
    row_pair_rules()
        .into_iter()
        .find(|rule| rule.bucket == CLASS)
        .expect("the class's rule")
}

/// Whether the class's rule takes the reference's row `e` and vep-rs's row `a`.
fn taken(e: &PairRow, a: &PairRow) -> bool {
    class_rule().matches_rows(std::slice::from_ref(e), std::slice::from_ref(a))
}

/// The bucket of the first row-aware rule of the registry that takes the reference's row `e`
/// and vep-rs's row `a`, in the registry's order.
fn first_row_rule(e: &PairRow, a: &PairRow) -> Option<&'static str> {
    row_pair_rules()
        .into_iter()
        .find(|rule| rule.matches_rows(std::slice::from_ref(e), std::slice::from_ref(a)))
        .map(|rule| rule.bucket)
}

/// The class takes vep-rs's `start_lost` on a Transcript row both of whose position columns
/// open at 1 and end on a number, the same on both rows up to a `/total` suffix: the terms of
/// the two rows are otherwise equal, both hold `inframe_deletion`, the reference's holds no
/// `start_retained_variant`, and the alternate half of both rows' Amino_acids is `-`.
#[test]
fn the_class_takes_a_lost_start_beside_an_inframe_deletion_from_the_first_base_and_no_other_row() {
    let pairs = table_pairs();
    let wrong: Vec<&str> = pairs
        .iter()
        .filter(|(_, e, a, is_taken)| taken(e, a) != *is_taken)
        .map(|(shape, ..)| *shape)
        .collect();
    assert_eq!(wrong, Vec::<&str>::new());
    let taken = pairs.iter().filter(|(.., taken)| *taken).count();
    assert_eq!((taken, pairs.len()), (9, 33));
    let spelt = pairs.iter().filter(|(_, e, ..)| e.indel.is_some());
    assert_eq!(
        spelt.count(),
        1,
        "the rows whose first column spells the alleles"
    );
}

/// No rule over the two consequence sets takes a pair of the table, so its rows reach the
/// row-aware rules; there the class's rule, the last but one the registry asks, is the one
/// rule that takes a pair of the class, it classifies and sets aside, and a pair it refuses is
/// taken by no rule.
#[test]
fn a_pair_of_the_class_passes_every_earlier_rule_and_is_set_aside_under_its_own_bucket() {
    let rules = row_pair_rules();
    let buckets: Vec<&str> = rules.iter().map(|rule| rule.bucket).collect();
    assert_eq!(
        buckets[buckets.len() - 2],
        CLASS,
        "the last rule but one of the registry"
    );
    assert_eq!(buckets.len(), 16);
    let rule = class_rule();
    assert!(rule.classifies && rule.excludes);
    for (shape, e, a, is_taken) in table_pairs() {
        assert_eq!(
            classify(&e.consequence, &a.consequence, &e.allele),
            None,
            "{shape}"
        );
        assert!(
            !excluded(&e.consequence, &a.consequence, &e.allele),
            "{shape}"
        );
        assert_eq!(first_row_rule(&e, &a), is_taken.then_some(CLASS), "{shape}");
    }
}

/// The same pair read to its ninth column and to its twelfth: a row read to its ninth
/// carries no Amino_acids, and the class takes the pair where both rows carry theirs, with or
/// without the Codons after them.
#[test]
fn a_row_read_without_its_amino_acids_is_no_row_of_the_class() {
    let (e, a) = FIRST_CODON_ROW.pair(TERMS.0, TERMS.1);
    assert!(taken(&e, &a));
    let nine = |consequence| {
        let row = Row {
            consequence,
            ..FIRST_CODON_ROW
        };
        let columns = row.columns();
        let first: &[&str; 9] = columns.first_chunk().expect("nine of twelve columns");
        pair_row(first)
    };
    let (short_e, short_a) = (nine(TERMS.0), nine(TERMS.1));
    assert_eq!(
        (short_e.amino_acids.as_deref(), short_e.codons.as_deref()),
        (None, None)
    );
    for (reference, candidate) in [(&short_e, &short_a), (&short_e, &a), (&e, &short_a)] {
        assert!(!taken(reference, candidate));
        assert_eq!(first_row_rule(reference, candidate), None);
    }
    let without_codons = |row: &PairRow| PairRow {
        codons: None,
        ..row.clone()
    };
    assert!(taken(&without_codons(&e), &without_codons(&a)));
}

/// The class reads the terms as a set: the order the two rows write them in, and a term
/// written twice, change nothing.
#[test]
fn the_terms_of_the_two_rows_are_read_as_sets() {
    for (reference, candidate) in [
        (
            "NMD_transcript_variant,inframe_deletion",
            "inframe_deletion,NMD_transcript_variant,start_lost",
        ),
        (
            "inframe_deletion,inframe_deletion",
            "start_lost,inframe_deletion",
        ),
        (
            "inframe_deletion&NMD_transcript_variant",
            "start_lost&inframe_deletion&NMD_transcript_variant",
        ),
    ] {
        let (e, a) = FIRST_CODON_ROW.pair(reference, candidate);
        assert!(taken(&e, &a), "{reference} | {candidate}");
    }
}

/// A pair of each other kind of row-aware rule, on rows that carry this class's positions
/// and residues where that rule admits them: the stop codon rules differ by
/// `stop_retained_variant` or `stop_gained`, or by `frameshift_variant` on rows whose
/// reference residues open with a stop, the rule of a splice term by a term the reference
/// alone holds, the rule of genetic code 26 by `synonymous_variant` against
/// `missense_variant` and the rule of a start codon kept by a `start_retained_variant`
/// vep-rs alone holds, and none of their pairs is one of this class.
#[test]
fn a_pair_of_another_row_aware_rule_is_left_to_that_rule() {
    let deletion = Row {
        name: "1993671",
        location: "19:13208046-13208050",
        cdna: "?-7057",
        cds: "-",
        amino_acids: "-",
        codons: "-",
        ..FIRST_CODON_ROW
    };
    let unknown_bases = Row {
        name: "1993671",
        location: "21:10567700-10567701",
        allele: "GAATAGNNNN",
        cdna: "775-776",
        cds: "423-424",
        ..FIRST_CODON_ROW
    };
    let settled = |consequence, amino_acids| {
        Row {
            consequence,
            amino_acids,
            codons: "-/GAATAGNNNN",
            ..unknown_bases
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
            ..FIRST_CODON_ROW
        }
        .read()
    };
    let pairs = [
        (
            "stop_retained_without_stop_codon_test",
            deletion.pair("stop_retained_variant", "coding_sequence_variant"),
        ),
        (
            "stop_codon_removed_without_stop_codon_test",
            deletion.pair("stop_retained_variant", "stop_lost"),
        ),
        (
            "stop_codon_settled_before_unknown_bases",
            (
                settled("frameshift_variant", "-"),
                settled("stop_gained,frameshift_variant", "-/E*XX"),
            ),
        ),
        (
            "frameshift_inside_a_held_stop_codon",
            Row {
                name: "1993671",
                location: "1:136",
                cdna: "36",
                cds: "30",
                amino_acids: "*/X",
                codons: "tgA/tg",
                ..FIRST_CODON_ROW
            }
            .pair("coding_sequence_variant", "frameshift_variant"),
        ),
        (
            "splice_term_past_the_record_span",
            Row {
                name: "1993671",
                location: "7:100998-101000",
                allele: "GTNCA",
                cdna: "998-1000",
                cds: "898-900",
                ..FIRST_CODON_ROW
            }
            .pair(
                "splice_donor_variant,frameshift_variant,intron_variant",
                "frameshift_variant,splice_region_variant",
            ),
        ),
        (
            "pachysolen_nuclear_code_read_as_standard",
            (
                on_code_26("synonymous_variant", "L"),
                on_code_26("missense_variant", "A/L"),
            ),
        ),
        (
            "start_codon_kept_without_five_prime_utr",
            FIRST_CODON_ROW.pair(
                "inframe_deletion",
                "inframe_deletion,start_retained_variant",
            ),
        ),
    ];
    for (bucket, (e, a)) in pairs {
        assert_eq!(first_row_rule(&e, &a), Some(bucket), "{bucket}");
        assert!(!taken(&e, &a), "{bucket}");
    }
}
