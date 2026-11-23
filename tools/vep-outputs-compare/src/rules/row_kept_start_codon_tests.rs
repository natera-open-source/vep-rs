// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row-aware pair rule of `start_codon_kept_without_five_prime_utr` on the fabricated
//! rows of `tests/fixtures/row_aware_pair_rules/tables/kept_start_codon_rows.tsv` and on rows
//! written here: the two forms of the term difference on each kind of insertion's or
//! deletion's row, the rows the rule refuses, its place among the rules of the registry, the
//! columns a row is read to, and the pairs of the other rules, which it leaves to them.

use super::super::{classify, excluded, row_pair_rules, ExclusionRule};
use super::*;

const CLASS: &str = "start_codon_kept_without_five_prime_utr";

const TABLE: &str =
    include_str!("../../tests/fixtures/row_aware_pair_rules/tables/kept_start_codon_rows.tsv");

const COLUMNS: [&str; 12] = [
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
    "taken",
];

/// One fabricated row of the default format as the rule reads it: its first column, its
/// Location, its allele, its Feature_type, its terms and its two position columns.
#[derive(Clone, Copy)]
struct Row<'a> {
    name: &'a str,
    location: &'a str,
    allele: &'a str,
    feature_type: &'a str,
    consequence: &'a str,
    cdna: &'a str,
    cds: &'a str,
}

impl<'a> Row<'a> {
    fn columns(&self) -> [&'a str; 9] {
        [
            self.name,
            self.location,
            self.allele,
            "ENSG90000000002",
            "ENST90000000002",
            self.feature_type,
            self.consequence,
            self.cdna,
            self.cds,
        ]
    }

    /// The row read to its CDS_position column, the ninth.
    fn read(&self) -> PairRow {
        pair_row(&self.columns())
    }

    /// The reference's row under the terms `reference` and vep-rs's under `candidate`.
    fn pair(&self, reference: &str, candidate: &str) -> (PairRow, PairRow) {
        let side = |consequence| Row {
            consequence,
            ..*self
        };
        (side(reference).read(), side(candidate).read())
    }

    /// The row with `span` in both of its position columns.
    fn spanning(self, span: &'a str) -> Row<'a> {
        Row {
            cdna: span,
            cds: span,
            ..self
        }
    }

    /// The row read to its Codons column, the twelfth, with the three cells after its
    /// CDS_position.
    fn read_with(&self, after: [&'a str; 3]) -> PairRow {
        let mut columns = [""; 12];
        columns[..9].copy_from_slice(&self.columns());
        columns[9..].copy_from_slice(&after);
        pair_row_with_codons(&columns)
    }
}

/// Ten bases removed from the third base of a transcript whose coding sequence opens on
/// its first base, under an identifier.
const DELETION_ROW: Row = Row {
    name: "2750113",
    location: "21:30001-30010",
    allele: "-",
    feature_type: "Transcript",
    consequence: "",
    cdna: "3-12",
    cds: "3-12",
};

/// Three bases put between the second and the third base of such a transcript.
const INSERTION_ROW: Row = Row {
    name: "2750118",
    location: "21:30000-30001",
    allele: "GCA",
    feature_type: "Transcript",
    consequence: "",
    cdna: "2-3",
    cds: "2-3",
};

/// The terms of the reference's row and of vep-rs's on a frameshift of the class.
const TERMS: (&str, &str) = (
    "frameshift_variant,start_lost",
    "frameshift_variant,start_retained_variant",
);

/// The pairs of the table: a shape, the reference's row, vep-rs's row and whether the class
/// takes the pair.
fn table_pairs() -> Vec<(&'static str, PairRow, PairRow, bool)> {
    let mut lines = TABLE.lines();
    assert_eq!(lines.next(), Some(&COLUMNS.join("\t")[..]));
    let pair = |line: &'static str| {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), COLUMNS.len(), "{line}");
        let side = |cdna, cds, consequence| Row {
            name: cells[1],
            location: cells[2],
            allele: cells[3],
            feature_type: cells[4],
            consequence,
            cdna,
            cds,
        };
        let taken = match cells[11] {
            "yes" => true,
            "no" => false,
            other => panic!("neither yes nor no: {other}"),
        };
        (
            cells[0],
            side(cells[5], cells[7], cells[9]).read(),
            side(cells[6], cells[8], cells[10]).read(),
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

/// The class takes vep-rs's `start_retained_variant` on a Transcript row whose
/// cDNA_position and CDS_position are one span of one number or two on both rows, up to a
/// `/total` suffix, that opens on coding base 1, 2 or 3. The reference's row holds no
/// `start_retained_variant`, vep-rs's holds no `start_lost`, the two sets are otherwise
/// equal but for the reference's `start_lost`, and they hold `frameshift_variant`,
/// `inframe_insertion` or `inframe_deletion`.
#[test]
fn the_class_takes_a_kept_start_on_a_span_that_opens_on_the_start_codon_and_no_other_row() {
    let pairs = table_pairs();
    let wrong: Vec<&str> = pairs
        .iter()
        .filter(|(_, e, a, is_taken)| taken(e, a) != *is_taken)
        .map(|(shape, ..)| *shape)
        .collect();
    assert_eq!(wrong, Vec::<&str>::new());
    let taken = pairs.iter().filter(|(.., taken)| *taken).count();
    assert_eq!((taken, pairs.len()), (14, 36));
    let spelt = pairs.iter().filter(|(_, e, ..)| e.indel.is_some());
    assert_eq!(
        spelt.count(),
        1,
        "the rows whose first column spells the alleles"
    );
}

/// The two forms of the term difference on each of the three kinds of row. In the first the
/// reference's row holds `start_lost` and vep-rs's `start_retained_variant` in its place; in
/// the second the reference's row holds neither start term and vep-rs's adds
/// `start_retained_variant`.
#[test]
fn both_forms_of_the_term_difference_are_taken_on_each_kind_of_row() {
    let deletion_of_codons = Row {
        location: "21:30001-30003",
        ..DELETION_ROW
    }
    .spanning("3-5");
    for (row, reference, candidate) in [
        (
            DELETION_ROW,
            "frameshift_variant,start_lost",
            "frameshift_variant,start_retained_variant",
        ),
        (
            DELETION_ROW,
            "frameshift_variant",
            "frameshift_variant,start_retained_variant",
        ),
        (
            INSERTION_ROW,
            "inframe_insertion,start_lost",
            "inframe_insertion,start_retained_variant",
        ),
        (
            INSERTION_ROW,
            "inframe_insertion",
            "inframe_insertion,start_retained_variant",
        ),
        (
            deletion_of_codons,
            "start_lost,inframe_deletion",
            "inframe_deletion,start_retained_variant",
        ),
        (
            deletion_of_codons,
            "inframe_deletion",
            "inframe_deletion,start_retained_variant",
        ),
    ] {
        let (e, a) = row.pair(reference, candidate);
        assert!(taken(&e, &a), "{reference} | {candidate}");
        assert_eq!(first_row_rule(&e, &a), Some(CLASS), "{reference}");
    }
}

/// The nearest rows the class refuses, of four kinds: a cDNA position behind a 5' UTR, a
/// span that opens on the fourth coding base, a further term on the row of vep-rs or on the
/// row of the reference, and a row of vep-rs that keeps the reference's `start_lost` beside
/// the retained start.
#[test]
fn a_five_prime_utr_a_later_span_another_term_and_a_kept_lost_start_are_refused() {
    let behind_a_utr = Row {
        cdna: "104-113",
        ..DELETION_ROW
    };
    let later = DELETION_ROW.spanning("4-13");
    for (label, row, reference, candidate) in [
        ("a 5' UTR", behind_a_utr, TERMS.0, TERMS.1),
        ("past coding base 3", later, TERMS.0, TERMS.1),
        (
            "another term on the row of vep-rs",
            DELETION_ROW,
            TERMS.0,
            "frameshift_variant,start_retained_variant,NMD_transcript_variant",
        ),
        (
            "another term on the row of the reference",
            DELETION_ROW,
            "frameshift_variant,start_lost,NMD_transcript_variant",
            TERMS.1,
        ),
        (
            "a lost start on the row of vep-rs",
            DELETION_ROW,
            TERMS.0,
            "frameshift_variant,start_lost,start_retained_variant",
        ),
    ] {
        let (e, a) = row.pair(reference, candidate);
        assert!(!taken(&e, &a), "{label}");
        assert_eq!(first_row_rule(&e, &a), None, "{label}");
    }
    let (e, a) = DELETION_ROW.pair(TERMS.0, TERMS.1);
    assert!(taken(&e, &a), "the same row without any of the four");
}

/// No rule over the two consequence sets takes a pair of the table, so its rows reach the
/// row-aware rules; there the class's rule, the last the registry asks, is the one rule
/// that takes a pair of the class, it classifies and sets aside, and a pair it refuses is
/// taken by no rule.
#[test]
fn a_pair_of_the_class_reaches_the_last_rule_of_the_registry_and_no_other_takes_it() {
    let rules = row_pair_rules();
    let buckets: Vec<&str> = rules.iter().map(|rule| rule.bucket).collect();
    assert_eq!(
        buckets.last(),
        Some(&CLASS),
        "the last rule of the registry"
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

/// A reference row that holds both start terms is the pair of the set rule
/// `start_cooccurrence_swap`: vep-rs's row adds no term, and the class's rule refuses it.
#[test]
fn both_start_terms_on_the_row_of_the_reference_stay_with_the_set_rule() {
    let (e, a) = DELETION_ROW.pair(
        "frameshift_variant,start_lost,start_retained_variant",
        "frameshift_variant,start_retained_variant",
    );
    assert_eq!(
        classify(&e.consequence, &a.consequence, &e.allele),
        Some("start_cooccurrence_swap")
    );
    assert!(excluded(&e.consequence, &a.consequence, &e.allele));
    assert!(!taken(&e, &a));
}

/// The rule reads a row to its CDS_position: the same pair read to its ninth column and to
/// its twelfth, with residues and codons, is taken alike.
#[test]
fn a_row_is_read_to_its_ninth_column() {
    let (e, a) = DELETION_ROW.pair(TERMS.0, TERMS.1);
    assert_eq!(
        (e.amino_acids.as_deref(), e.codons.as_deref()),
        (None, None)
    );
    assert!(taken(&e, &a));
    let twelve = |consequence| {
        let row = Row {
            consequence,
            ..DELETION_ROW
        };
        row.read_with(["1-4", "MAGT/X", "atGGCGGGAACa/aa"])
    };
    let (long_e, long_a) = (twelve(TERMS.0), twelve(TERMS.1));
    assert_eq!(long_e.amino_acids.as_deref(), Some("MAGT/X"));
    for (reference, candidate) in [(&long_e, &long_a), (&e, &long_a), (&long_e, &a)] {
        assert!(taken(reference, candidate));
    }
}

/// The class reads the terms as a set: the order the two rows write them in, a term
/// written twice and the separator of the VCF format change nothing.
#[test]
fn the_terms_of_the_two_rows_are_read_as_sets() {
    for (reference, candidate) in [
        (
            "NMD_transcript_variant,start_lost,frameshift_variant",
            "start_retained_variant,frameshift_variant,NMD_transcript_variant",
        ),
        (
            "frameshift_variant,start_lost,start_lost",
            "frameshift_variant,start_retained_variant",
        ),
        (
            "frameshift_variant&start_lost&NMD_transcript_variant",
            "frameshift_variant&start_retained_variant&NMD_transcript_variant",
        ),
    ] {
        let (e, a) = DELETION_ROW.pair(reference, candidate);
        assert!(taken(&e, &a), "{reference} | {candidate}");
    }
}

/// A pair of each other kind of row-aware rule, on rows that carry this class's positions
/// where that rule admits them: the stop codon rules differ by `stop_retained_variant` or
/// `stop_gained`, the rule of a splice term by a term the reference alone holds, the rule
/// of genetic code 26 by `synonymous_variant` against `missense_variant` and the rule of a
/// start codon removed whole by a `start_lost` vep-rs alone holds, and none of their pairs
/// is one of this class.
#[test]
fn a_pair_of_any_other_row_aware_rule_is_no_pair_of_the_class() {
    let coding = |row: Row<'static>, consequence, amino_acids, codons| {
        Row { consequence, ..row }.read_with(["-", amino_acids, codons])
    };
    let first_codon = Row {
        location: "21:30001-30003",
        ..DELETION_ROW
    }
    .spanning("1-3");
    let without_a_stop_codon_test = Row {
        location: "19:13208046-13208050",
        cdna: "?-7057",
        cds: "-",
        ..DELETION_ROW
    };
    let unknown_bases = Row {
        location: "21:10567700-10567701",
        allele: "GAATAGNNNN",
        ..DELETION_ROW
    }
    .spanning("2-3");
    let on_code_26 = Row {
        location: "7:1000",
        allele: "A",
        ..DELETION_ROW
    }
    .spanning("3");
    let pairs = [
        (
            "stop_retained_without_stop_codon_test",
            without_a_stop_codon_test.pair("stop_retained_variant", "coding_sequence_variant"),
        ),
        (
            "stop_codon_removed_without_stop_codon_test",
            without_a_stop_codon_test.pair("stop_retained_variant", "stop_lost"),
        ),
        (
            "stop_codon_settled_before_unknown_bases",
            (
                coding(unknown_bases, "frameshift_variant", "-", "-/GAATAGNNNN"),
                coding(
                    unknown_bases,
                    "stop_gained,frameshift_variant",
                    "-/E*XX",
                    "-/GAATAGNNNN",
                ),
            ),
        ),
        (
            "splice_term_past_the_record_span",
            Row {
                location: "7:100998-101000",
                allele: "GTNCA",
                cdna: "998-1000",
                cds: "898-900",
                ..DELETION_ROW
            }
            .pair(
                "splice_donor_variant,frameshift_variant,intron_variant",
                "frameshift_variant,splice_region_variant",
            ),
        ),
        (
            "pachysolen_nuclear_code_read_as_standard",
            (
                coding(on_code_26, "synonymous_variant", "L", "ctG/ctA"),
                coding(on_code_26, "missense_variant", "A/L", "ctG/ctA"),
            ),
        ),
        (
            "start_codon_removed_without_five_prime_utr",
            (
                coding(first_codon, "inframe_deletion", "M/-", "ATG/-"),
                coding(first_codon, "start_lost,inframe_deletion", "M/-", "ATG/-"),
            ),
        ),
    ];
    for (bucket, (e, a)) in pairs {
        assert_eq!(first_row_rule(&e, &a), Some(bucket), "{bucket}");
        assert!(!taken(&e, &a), "{bucket}");
    }
}
