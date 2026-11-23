// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The row-aware pair rules on the fabricated rows of
//! `tests/fixtures/row_aware_pair_rules/tables`: what a row carries of its record's alleles,
//! the three readings `splice_term_past_the_record_span` makes of a record, and the verdict
//! of a rule on the reference's and vep-rs's row of one key.

use super::super::row_pair_rules;
use super::*;

/// The rows of one table of `tests/fixtures/row_aware_pair_rules/tables` under its line of
/// column names.
macro_rules! table {
    ($name:literal, $columns:expr) => {
        table_rows(
            include_str!(concat!(
                "../../tests/fixtures/row_aware_pair_rules/tables/",
                $name
            )),
            &$columns,
        )
    };
}

fn table_rows(text: &'static str, columns: &[&str]) -> Vec<Vec<&'static str>> {
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some(&columns.join("\t")[..]));
    let rows: Vec<Vec<&str>> = lines.map(|line| line.split('\t').collect()).collect();
    assert!(rows.iter().all(|row| row.len() == columns.len()));
    rows
}

/// What a rule reads of one fabricated row of the default format: its first column, its
/// Location, its allele, its feature's type, its terms and its two position columns.
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

/// A deletion's row under an identifier, without terms or positions.
const ROW: Row = Row {
    name: "1993671",
    location: "19:13208046-13208050",
    allele: "-",
    feature_type: "Transcript",
    consequence: "",
    cdna: "-",
    cds: "-",
};

impl Row<'_> {
    fn read(&self) -> PairRow {
        pair_row(&[
            self.name,
            self.location,
            self.allele,
            "ENSG1",
            "ENST00000635895",
            self.feature_type,
            self.consequence,
            self.cdna,
            self.cds,
        ])
    }
}

/// Whether the row-aware pair rule the registry holds under `bucket` takes the reference's
/// row `e` and vep-rs's row `a`.
fn taken_by(bucket: &str, e: &PairRow, a: &PairRow) -> bool {
    let rule = row_pair_rules()
        .into_iter()
        .find(|rule| rule.bucket == bucket)
        .expect("a row-aware pair rule of the registry");
    rule.matches_rows(std::slice::from_ref(e), std::slice::from_ref(a))
}

const PAST_THE_SPAN: &str = "splice_term_past_the_record_span";

fn yes(cell: &str) -> bool {
    match cell {
        "yes" => true,
        "no" => false,
        other => panic!("neither yes nor no: {other}"),
    }
}

/// `number` written in the decimal digits of another script than `0` to `9`.
fn other_digits(number: u64) -> String {
    let digit = |digit: char| char::from_u32(0x0660 + digit.to_digit(10).unwrap()).unwrap();
    number.to_string().chars().map(digit).collect()
}

const INSERTION_RANGE_COLUMNS: [&str; 6] = [
    "shape",
    "first_column",
    "location",
    "allele",
    "taken",
    "insertion",
];

/// The reference's `coding_sequence_variant` against vep-rs's `stop_retained_variant` on an
/// allele of three between two neighbouring positions, where the first column spells an empty
/// reference allele and two alternates of other lengths.
#[test]
fn the_insertion_range_class_takes_an_allele_no_alternate_names_beside_an_empty_reference_allele() {
    let row = |consequence| {
        let row = Row {
            name: "19_13208046_-/CG/T",
            location: "19:13208046-13208047",
            allele: "TGA",
            consequence,
            cdna: "7054-7055",
            cds: "6977-6978",
            ..ROW
        };
        row.read()
    };
    assert!(stop_codon_overlap_insertion_range_pair(
        &row("coding_sequence_variant"),
        &row("stop_retained_variant")
    ));
}

/// The reference's `stop_retained_variant` against vep-rs's `coding_sequence_variant` on an
/// allele of one base, where the first column spells a reference allele and two alternates
/// and none of them, the shared first base dropped, is the row's allele. The record's alleles
/// then name no shape for the row, and its Location decides as under an identifier: over
/// three positions the allele changes the length and the class takes the pair, at one
/// position it does not. An alternate that is the row's allele decides before the Location.
#[test]
fn the_three_base_class_reads_an_allele_no_alternate_names_by_its_location() {
    let row = |location, allele, consequence| {
        let row = Row {
            name: "19_13208044_ACGTT/ACA/AG",
            location,
            allele,
            consequence,
            cdna: "7054-7056",
            cds: "6977-6979",
            ..ROW
        };
        row.read()
    };
    let (retained, coding) = ("stop_retained_variant", "coding_sequence_variant");
    let taken = |location, allele| {
        let (e, a) = (
            row(location, allele, retained),
            row(location, allele, coding),
        );
        let class = "stop_retained_three_base_stop_codon_test";
        assert_eq!(
            stop_retained_three_base_stop_codon_test_pair(&e, &a),
            taken_by(class, &e, &a)
        );
        (e.indel, taken_by(class, &e, &a))
    };
    assert_eq!(taken("19:13208046-13208048", "A"), (None, true));
    assert_eq!(taken("19:13208046", "A"), (None, false));
    assert_eq!(taken("19:13208045", "G"), (Some(true), true));
}

#[test]
fn a_row_carries_whether_its_record_s_alleles_make_every_sequence_allele_an_insertion() {
    let rows = table!("insertion_range_rows.tsv", INSERTION_RANGE_COLUMNS);
    assert_eq!(rows.len(), 14);
    for row in rows {
        let read = Row {
            name: row[1],
            location: row[2],
            allele: row[3],
            consequence: "coding_sequence_variant",
            cdna: "7054-7055",
            cds: "6977-6978",
            ..ROW
        };
        let carried = (row[5] != "none").then(|| yes(row[5]));
        assert_eq!(read.read().insertion, carried, "{}", row[0]);
    }
    assert!(reference_allele_is_empty("-", &["T"]));
    assert!(!reference_allele_is_empty("A", &["T", "G"]));
}

/// Every row stands under an identifier, so none carries its record's alleles.
#[test]
fn without_its_record_the_past_the_span_class_takes_a_longer_allele_over_three_positions_of_the_transcript_and_no_other_row(
) {
    let rows = table!(
        "past_the_span_rows.tsv",
        [
            "shape",
            "reference_terms",
            "candidate_terms",
            "allele",
            "location",
            "cdna",
            "cds",
            "candidate_cdna",
            "candidate_cds",
            "feature_type",
            "taken",
        ]
    );
    let mut taken = 0;
    let mut wrong: Vec<&str> = Vec::new();
    for row in &rows {
        let reference = Row {
            location: row[4],
            allele: row[3],
            feature_type: row[9],
            consequence: row[1],
            cdna: row[5],
            cds: row[6],
            ..ROW
        };
        let candidate = Row {
            consequence: row[2],
            cdna: row[7],
            cds: row[8],
            ..reference
        };
        let (e, a) = (reference.read(), candidate.read());
        assert_eq!((&e.reference, &a.reference), (&None, &None), "{}", row[0]);
        taken += usize::from(yes(row[10]));
        if taken_by(PAST_THE_SPAN, &e, &a) != yes(row[10]) {
            wrong.push(row[0]);
        }
    }
    assert_eq!(wrong, Vec::<&str>::new());
    assert_eq!((taken, rows.len()), (13, 34));
}

const REPLACED_REFERENCE_COLUMNS: [&str; 8] = [
    "shape", "pos", "ref", "alts", "location", "allele", "read", "bases",
];

#[test]
fn the_bases_an_allele_stands_in_for_are_read_by_position_from_the_record_s_alleles() {
    let rows = table!("replaced_reference.tsv", REPLACED_REFERENCE_COLUMNS);
    assert_eq!(rows.len(), 13);
    for row in rows {
        let alternates: Vec<&str> = row[3].split(',').collect();
        let position: i128 = row[1].parse().expect("a position");
        assert_eq!(
            replaced_reference(position, row[2], &alternates, row[4], row[5]),
            yes(row[6]).then(|| row[7].to_string()),
            "{}",
            row[0]
        );
    }
    assert_eq!(
        replaced_reference(100998, "ACG", &["GTTCA"], "", "GTTCA"),
        None,
        "a row without a Location"
    );
}

#[test]
fn a_first_column_that_spells_the_alleles_spells_the_record_s_position_before_them() {
    assert_eq!(uploaded_position("7_100998_ACG/GTTCA"), Some(100998));
    assert_eq!(uploaded_position("7_100997_tacg/t/-"), Some(100997));
    assert_eq!(uploaded_position("chrUn_made_up_100998_A/T"), Some(100998));
    for identifier in ["1993671", "rs1993671", "7_100998", "made_up_7"] {
        assert_eq!(uploaded_position(identifier), None, "{identifier}");
    }
}

#[test]
fn a_row_carries_the_bases_its_first_column_spells_and_none_under_an_identifier() {
    let reference = |name: &str, location, allele| {
        let row = Row {
            name,
            location,
            allele,
            consequence: "frameshift_variant",
            cdna: "998-1000",
            cds: "898-900",
            ..ROW
        };
        row.read().reference
    };
    let rows = table!("replaced_reference.tsv", REPLACED_REFERENCE_COLUMNS);
    let read: Vec<&Vec<&str>> = rows.iter().filter(|row| yes(row[6])).collect();
    assert_eq!(read.len(), 8);
    for row in read {
        let spelt = format!("7_{}_{}/{}", row[1], row[2], row[3].replace(',', "/"));
        assert_eq!(
            reference(&spelt, row[4], row[5]).as_deref(),
            Some(row[7]),
            "{}",
            row[0]
        );
        assert_eq!(reference("1993671", row[4], row[5]), None, "{}", row[0]);
    }
    assert_eq!(
        reference("7_100998_ACG/GTTCAT", "7:100998-101000", "GTTCA"),
        None,
        "an alternate that writes no such row"
    );
}

/// The regions a table writes as `first,last` pairs parted by `;`.
fn regions(cell: &str) -> Vec<(i64, i64)> {
    let region = |pair: &str| {
        let (first, last) = pair.split_once(',').expect("a first and a last offset");
        (first.parse().unwrap(), last.parse().unwrap())
    };
    cell.split(';').map(region).collect()
}

fn letters(sequence: &str) -> Vec<char> {
    sequence.chars().collect()
}

#[test]
fn two_alleles_differ_in_runs_of_indices_and_one_base_gives_the_span() {
    let rows = table!(
        "differing_regions.tsv",
        ["shape", "reference", "allele", "regions"]
    );
    assert_eq!(rows.len(), 7);
    for row in rows {
        assert_eq!(
            differing_regions(&letters(row[1]), &letters(row[2])),
            regions(row[3]),
            "{}",
            row[0]
        );
    }
}

/// The terms of a table's cell among the seven a row of `intron_rows` holds.
fn term_set(cell: &str) -> TermSet {
    let terms: BTreeSet<&str> = cell.split(',').filter(|term| !term.is_empty()).collect();
    assert!(terms.iter().all(|term| INTRON_ROW_TERMS.contains(term)));
    intron_row_terms(&terms)
}

/// The shapes of a table on which `intron_rows` gives other terms than the table's. Every
/// offset counts from the span's first base.
fn intron_rows_that_differ(rows: &[Vec<&'static str>]) -> Vec<&'static str> {
    let number = |cell: &str| -> i64 { cell.parse().expect("an offset") };
    let differs = |row: &&Vec<&str>| {
        let written = intron_rows(
            &regions(row[1]),
            number(row[2]),
            number(row[3]),
            number(row[4]),
            yes(row[5]),
            yes(row[6]),
        );
        written != (term_set(row[7]), term_set(row[8]))
    };
    rows.iter().filter(differs).map(|row| row[0]).collect()
}

const INTRON_ROWS_COLUMNS: [&str; 9] = [
    "shape",
    "regions",
    "span",
    "intron_start",
    "intron_end",
    "forward",
    "wide_exons",
    "reference_terms",
    "candidate_terms",
];

#[test]
fn the_regions_give_the_reference_s_terms_and_the_span_vep_rs_s_beside_one_intron() {
    let rows = table!("intron_rows.tsv", INTRON_ROWS_COLUMNS);
    assert_eq!(rows.len(), 18);
    assert_eq!(intron_rows_that_differ(&rows), Vec::<&str>::new());
}

/// `splice_region_variant` follows the order of the regions: the reference's row holds it
/// where the last region read before one meets the intron's first two or last two bases lies
/// in a window of the term, vep-rs's where any region does.
#[test]
fn the_last_region_read_decides_the_reference_s_region_term_and_any_region_vep_rs_s() {
    let rows = table!("intron_rows_region_term.tsv", INTRON_ROWS_COLUMNS);
    assert_eq!(rows.len(), 5);
    assert_eq!(intron_rows_that_differ(&rows), Vec::<&str>::new());
}

/// The places of an intron beside a span of `span` bases that a row's positions allow, for
/// regions that run `reach` bases past the span's end, in their order.
fn places(span: i64, reach: i64, positions: [&str; 4]) -> Vec<(bool, i64, i64)> {
    let mut places = Vec::new();
    intron_placements(span, reach, &positions, |forward, start, end| {
        places.push((forward, start, end));
        false
    });
    places
}

fn sorted(mut places: Vec<(bool, i64, i64)>) -> Vec<(bool, i64, i64)> {
    places.sort_unstable();
    places
}

#[test]
fn two_numbered_cdna_positions_as_far_apart_as_the_span_put_the_intron_after_the_span_or_before_it()
{
    let places = places(3, 2, ["998", "1000", "898", "900"]);
    let after: Vec<(bool, i64, i64)> = places.iter().copied().filter(|p| p.1 > 0).collect();
    let starts: BTreeSet<(bool, i64)> = after.iter().map(|p| (p.0, p.1)).collect();
    let one_to_seven_bases_on: BTreeSet<(bool, i64)> = [true, false]
        .into_iter()
        .flat_map(|forward| (3..10).map(move |start| (forward, start)))
        .collect();
    assert_eq!(starts, one_to_seven_bases_on);
    assert!(after.iter().all(|(_, start, end)| end >= start));
    assert!(after.contains(&(true, 3, 3)), "an intron of one base");
    assert!(
        after.contains(&(false, 9, 10)),
        "a two-base intron seven bases on"
    );
    assert!(
        after.iter().any(|(_, start, end)| end - start > 1_000_000),
        "an intron longer than any window"
    );
    let before: Vec<(bool, i64, i64)> = places.iter().copied().filter(|p| p.1 <= 0).collect();
    let sixteen_bases_or_fewer_that_end_before_the_span: Vec<(bool, i64, i64)> = [true, false]
        .into_iter()
        .flat_map(|forward| {
            (-16..0).flat_map(move |start| (start..0).map(move |end| (forward, start, end)))
        })
        .collect();
    assert_eq!(
        sorted(before.clone()),
        sorted(sixteen_bases_or_fewer_that_end_before_the_span)
    );
    assert!(
        before.contains(&(true, -2, -1)),
        "a two-base intron that ends beside the span"
    );
    assert_eq!(
        self::places(3, 2, ["997", "1000", "897", "900"]),
        [],
        "numbers further apart than the span"
    );
}

#[test]
fn two_numbered_cdna_positions_closer_than_the_span_put_an_intron_of_the_missing_length_inside_it()
{
    assert_eq!(
        sorted(places(5, 2, ["998", "1000", "898", "900"])),
        [(false, 1, 2), (false, 2, 3), (true, 1, 2), (true, 2, 3)],
        "two bases of five, the span's ends kept"
    );
    assert_eq!(
        sorted(places(3, 2, ["999", "1000", "899", "900"])),
        [(false, 1, 1), (true, 1, 1)],
        "the middle base of three"
    );
    assert_eq!(
        places(2, 2, ["1000", "1000", "900", "900"]),
        [],
        "no base inside a span of two"
    );
}

/// The starts, or the ends, of the places on one strand.
fn ends_on(places: &[(bool, i64, i64)], strand: bool, starts: bool) -> BTreeSet<i64> {
    let on_the_strand = places.iter().filter(|place| place.0 == strand);
    on_the_strand
        .map(|place| if starts { place.1 } else { place.2 })
        .collect()
}

#[test]
fn a_half_numbered_cdna_position_puts_one_end_of_the_intron_inside_the_span() {
    let places = self::places(2, 5, ["1000", "?", "900", "?"]);
    assert_eq!(ends_on(&places, true, true), BTreeSet::from([1]));
    assert_eq!(ends_on(&places, false, false), BTreeSet::from([0]));
    assert!(places.iter().all(|&(forward, start, end)| if forward {
        end >= 1
    } else {
        start <= 0
    }));
    let mirrored = self::places(2, 5, ["?", "901", "?", "801"]);
    assert_eq!(ends_on(&mirrored, true, false), BTreeSet::from([0]));
    assert_eq!(ends_on(&mirrored, false, true), BTreeSet::from([1]));
    let three = self::places(3, 5, ["999", "?", "899", "?"]);
    assert_eq!(ends_on(&three, true, true), BTreeSet::from([1, 2]));
}

#[test]
fn an_unnumbered_cdna_position_puts_the_span_inside_the_intron() {
    let places = self::places(2, 3, ["-", "-", "-", "-"]);
    let strands: BTreeSet<bool> = places.iter().map(|place| place.0).collect();
    assert_eq!(strands, BTreeSet::from([true, false]));
    assert!(places.iter().all(|&(_, start, end)| start <= 0 && end >= 1));
    assert!(places.contains(&(true, 0, 1)), "an intron that is the span");
    assert!(
        places.contains(&(false, -17, 42)),
        "the furthest start and end a window reaches"
    );
    assert!(!places.contains(&(true, 0, 43)));
    let far_ends = places
        .iter()
        .filter(|&&(forward, start, end)| forward && start == 0 && end > 1_000_000);
    assert_eq!(far_ends.count(), 1, "a far end, once");
    let reached = self::places(2, 30, ["-", "-", "-", "-"]);
    let reached: BTreeSet<i64> = reached.iter().map(|place| place.2).collect();
    assert!(reached.contains(&(1 + 30 + 18)), "an end the regions reach");
    assert!(!reached.contains(&(1 + 30 + 19)));
}

/// Each pair stands in its own block of 10,000 positions of chromosome 7, on one of two
/// layouts of the block. On the forward strand's an exon ends at position 1000 (cDNA 1000,
/// CDS 900), its intron runs from 1001 to 1100 and the next exon starts at 1101 (cDNA 1001).
/// On the reverse strand's an exon starts at 1000 (cDNA 901, CDS 801) and runs down the
/// block, the intron before it runs from its 3' end at 1001 to its 5' end at 1100, and the
/// exon before that ends at 1101 (cDNA 900). A two-base intron runs from 1001 to 1002 and a
/// twenty-base intron from 1001 to 1020, each in the place of the layout's intron.
///
/// With the record a pair is set aside where an allele longer than its two or more reference
/// bases gives, at one place of the intron, the reference's terms from its regions and
/// vep-rs's from its span, which lies inside an exon, across its end or in the intron.
/// Without the record only three or more positions of the transcript's sequence are set
/// aside, vep-rs's row holding none of the six. A pure insertion, an allele that is not the
/// longer one, a term gone whose window holds a base of the span, a differing term outside
/// the class and differing positions stay differences either way. A pair whose first column
/// spells the alleles has the record on every reading.
#[test]
fn the_past_the_span_class_gives_each_kind_of_record_its_verdict_with_the_record_and_without_it() {
    let rows = table!(
        "record_kinds.tsv",
        [
            "kind",
            "pos",
            "ref",
            "alt",
            "first_column",
            "location",
            "allele",
            "cdna",
            "cds",
            "candidate_cdna",
            "reference_terms",
            "candidate_terms",
            "with_record",
            "without_record",
        ]
    );
    let kinds: BTreeSet<&str> = rows.iter().map(|row| row[0]).collect();
    assert_eq!((rows.len(), kinds.len()), (85, 85));
    let mut wrong: Vec<&str> = Vec::new();
    let mut by_the_input_alone = 0;
    for row in &rows {
        let kind = row[0];
        let reference = Row {
            name: row[4],
            location: row[5],
            allele: row[6],
            consequence: row[10],
            cdna: row[7],
            cds: row[8],
            ..ROW
        };
        let candidate = Row {
            consequence: row[11],
            cdna: row[9],
            ..reference
        };
        let (e, a) = (reference.read(), candidate.read());
        let position: i128 = row[1].parse().expect("a position");
        let bases = replaced_reference(position, row[2], &[row[3]], row[5], &e.allele);
        assert!(bases.is_some(), "{kind}");
        let verdicts = if e.reference.is_none() {
            let without_record = taken_by(PAST_THE_SPAN, &e, &a);
            let with = |row: &PairRow| PairRow {
                reference: bases.clone(),
                ..row.clone()
            };
            (
                taken_by(PAST_THE_SPAN, &with(&e), &with(&a)),
                without_record,
            )
        } else {
            assert_eq!((&e.reference, &a.reference), (&bases, &bases), "{kind}");
            let spelt = taken_by(PAST_THE_SPAN, &e, &a);
            (spelt, spelt)
        };
        let (with_record, without_record) = (yes(row[12]), yes(row[13]));
        if verdicts != (with_record, without_record) {
            wrong.push(kind);
        }
        if with_record && !without_record {
            // A pair the record alone sets aside has a Location of two positions or a cDNA
            // position that is not two numbers.
            by_the_input_alone += 1;
            let digits = row[7].replace('-', "");
            let numbered = !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit());
            assert!(
                positions_named(&e.location) == Some(2) || !numbered,
                "{kind}"
            );
        }
    }
    assert_eq!(wrong, Vec::<&str>::new());
    assert!(by_the_input_alone > 0);
}

/// Two bases across an exon's end replaced by seven, the first column spelling the alleles:
/// the span lacks a cDNA position at its end, and the intron's fifth base lies past it. A row
/// there can differ by the terms of a stop codon class beside those of the past-the-span
/// class, which takes the pair where the stop codon class's rule takes the two rows read with
/// vep-rs's splice and intron terms on both, and its own reading takes them with the
/// reference's other terms on both.
mod a_row_that_differs_by_the_terms_of_a_stop_codon_class_as_well {
    use super::*;

    const REFERENCE: &str =
        "splice_donor_variant,splice_donor_5th_base_variant,stop_retained_variant,intron_variant";
    const CANDIDATE: &str = "splice_donor_variant,stop_lost,intron_variant";

    /// Whether the past-the-span class's rule, and the rule of
    /// `stop_codon_removed_without_stop_codon_test`, takes the pair of rows.
    fn verdicts(reference: &str, candidate: &str) -> (bool, bool) {
        let row = |consequence| {
            let row = Row {
                name: "7_101000_AC/GTTCAGG",
                location: "7:101000-101001",
                allele: "GTTCAGG",
                consequence,
                cdna: "1000-?",
                cds: "900-?",
                ..ROW
            };
            row.read()
        };
        let (e, a) = (row(reference), row(candidate));
        assert_eq!(
            (e.reference.as_deref(), a.reference.as_deref()),
            (Some("AC"), Some("AC"))
        );
        (
            taken_by(PAST_THE_SPAN, &e, &a),
            taken_by("stop_codon_removed_without_stop_codon_test", &e, &a),
        )
    }

    #[test]
    fn is_taken_by_the_past_the_span_class_s_rule() {
        assert_eq!(verdicts(REFERENCE, CANDIDATE), (true, false));
    }

    #[test]
    fn is_taken_by_neither_rule_with_one_further_term_outside_both_classes() {
        let further = "splice_donor_variant,stop_lost,3_prime_UTR_variant,intron_variant";
        assert_eq!(verdicts(REFERENCE, further), (false, false));
    }

    #[test]
    fn is_taken_by_the_stop_codon_class_s_rule_where_it_differs_by_that_class_s_terms_alone() {
        let stop_codon_alone = "splice_donor_variant,stop_retained_variant,intron_variant";
        assert_eq!(verdicts(stop_codon_alone, CANDIDATE), (false, true));
    }
}

/// One row of an intron's feature under an identifier.
fn intron_row(location: &str, allele: &str, cdna: &str) -> PairRow {
    pair_row(&[
        "rs1",
        location,
        allele,
        "g1",
        "T1",
        "Transcript",
        "intron_variant",
        cdna,
        "-",
    ])
}

#[test]
fn the_row_readers_of_a_location_s_ends_take_the_digits_0_to_9_alone() {
    let spans_five = format!("21:{}-{}", other_digits(1), other_digits(5));
    assert!(indel_row(&intron_row("21:1-5", "AC", "-")));
    assert!(!indel_row(&intron_row(&spans_five, "AC", "-")));
    assert!(length_changing_row(&intron_row("21:1-5", "RY", "-")));
    assert!(!length_changing_row(&intron_row(&spans_five, "RY", "-")));
    let spans_two = format!("21:{}-{}", other_digits(1), other_digits(2));
    assert!(insertion_row(&intron_row("21:1-2", "AC", "-")));
    assert!(!insertion_row(&intron_row(&spans_two, "AC", "-")));
}

#[test]
fn the_readers_of_a_row_s_transcript_positions_take_the_digits_0_to_9_alone() {
    let numbered = places(3, 20, ["10", "12", "4", "6"]);
    let unnumbered = places(3, 20, ["?", "?", "4", "6"]);
    assert_ne!(numbered, unnumbered);
    let (ten, twelve) = (other_digits(10), other_digits(12));
    assert_eq!(places(3, 20, [&ten, &twelve, "4", "6"]), unnumbered);
    assert!(longer_allele_inside_one_exon(&intron_row(
        "21:100-102",
        "ACGTA",
        "10-12"
    )));
    assert!(!longer_allele_inside_one_exon(&intron_row(
        "21:100-102",
        "ACGTA",
        &format!("{ten}-{twelve}")
    )));
}

const INSERTION_RANGE: &str = "stop_codon_overlap_insertion_range";
const GENOMIC_SPAN: &str = "stop_codon_overlap_genomic_span";
const UNKNOWN_BESIDE_REGION: &str = "splice_region_variant,coding_sequence_variant";
const KEPT_BESIDE_REGION: &str = "splice_region_variant,stop_retained_variant";

/// Three bases, the last unknown, between an intron's last base and the exon base that ends
/// a stop codon, the first column spelling the record's alleles.
const UNKNOWN_BASE_ROW: Row = Row {
    name: "21_35201927_G/GAAN",
    location: "21:35201927-35201928",
    allele: "AAN",
    feature_type: "Transcript",
    consequence: "",
    cdna: "3048-3049",
    cds: "3047-3048",
};

/// The reference's row of `columns` under the terms `reference` and vep-rs's under
/// `candidate`.
fn rows_under(columns: Row<'_>, reference: &str, candidate: &str) -> (PairRow, PairRow) {
    let read = |consequence: &str| {
        Row {
            consequence,
            ..columns
        }
        .read()
    };
    (read(reference), read(candidate))
}

/// Whether the rule of `stop_codon_overlap_insertion_range` takes `columns` under the two
/// sets of terms, by the predicate and by the registry alike.
fn in_the_insertion_range(columns: Row<'_>, reference: &str, candidate: &str) -> bool {
    let (e, a) = rows_under(columns, reference, candidate);
    let taken = stop_codon_overlap_insertion_range_pair(&e, &a);
    assert_eq!(taken_by(INSERTION_RANGE, &e, &a), taken);
    taken
}

/// Beside `splice_region_variant` the first pair of the class is taken where the allele
/// holds an unknown base, wherever it stands, with the record's alleles and under an
/// identifier; the split codon's rule refuses those rows and takes the same row of known
/// bases, which this one refuses.
#[test]
fn the_insertion_range_class_takes_its_first_pair_beside_splice_region_variant_for_an_unknown_base()
{
    let (reference, candidate) = (UNKNOWN_BESIDE_REGION, KEPT_BESIDE_REGION);
    for columns in [
        UNKNOWN_BASE_ROW,
        Row {
            name: "1993671",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            name: "1993671",
            allele: "ANA",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            name: "21_35201927_G/GNAAAAT",
            allele: "NAAAAT",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            cdna: "3048-3049/3600",
            cds: "3047-?/3048",
            ..UNKNOWN_BASE_ROW
        },
    ] {
        let label = format!("{} {}", columns.name, columns.allele);
        assert!(
            in_the_insertion_range(columns, reference, candidate),
            "{label}"
        );
        let (e, a) = rows_under(columns, reference, candidate);
        assert!(!taken_by(GENOMIC_SPAN, &e, &a), "{label}");
    }
    for (name, allele) in [("21_35201927_G/GAAT", "AAT"), ("1993671", "AAA")] {
        let known = Row {
            name,
            allele,
            ..UNKNOWN_BASE_ROW
        };
        assert!(
            !in_the_insertion_range(known, reference, candidate),
            "{allele}"
        );
        let (e, a) = rows_under(known, reference, candidate);
        assert!(taken_by(GENOMIC_SPAN, &e, &a), "{allele}");
    }
}

/// The sides exchanged are no row of the class: the reference's `stop_retained_variant`
/// against vep-rs's `coding_sequence_variant` on an allele holding an unknown base is the
/// shape of a stop codon test that is not made.
#[test]
fn the_first_pair_beside_splice_region_variant_is_not_taken_with_the_sides_exchanged() {
    let (reference, candidate) = (UNKNOWN_BESIDE_REGION, KEPT_BESIDE_REGION);
    assert!(in_the_insertion_range(
        UNKNOWN_BASE_ROW,
        reference,
        candidate
    ));
    assert!(!in_the_insertion_range(
        UNKNOWN_BASE_ROW,
        candidate,
        reference
    ));
    let (e, a) = rows_under(UNKNOWN_BASE_ROW, candidate, reference);
    assert!(taken_by("stop_retained_without_stop_codon_test", &e, &a));
}

/// A row off the shape is not taken: an allele holding a letter that is neither a base nor
/// `N`, a deletion, a row without a cDNA end or a CDS start, a Location of one position or
/// of three, a record whose reference allele keeps a base, another feature's row, and terms
/// that are not the pair. The rejoined row of an allele `--minimal` split from a record of
/// two alternates, which prints the record's one-base Location and no position, is one of
/// them.
#[test]
fn a_row_beside_splice_region_variant_off_the_shape_is_not_taken() {
    let (reference, candidate) = (UNKNOWN_BESIDE_REGION, KEPT_BESIDE_REGION);
    assert!(in_the_insertion_range(
        UNKNOWN_BASE_ROW,
        reference,
        candidate
    ));
    for columns in [
        Row {
            name: "1993671",
            allele: "AAR",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            name: "1993671",
            allele: "-",
            location: "21:35201928-35201930",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            cdna: "3048-?",
            cds: "3047-?",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            cds: "?-3048",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            name: "1993671",
            location: "21:35201928",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            name: "1993671",
            location: "21:35201927-35201929",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            name: "21_35201927_GA/GAAN",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            feature_type: "RegulatoryFeature",
            ..UNKNOWN_BASE_ROW
        },
        Row {
            name: "21_35201926_AG/AGAAN/AGNAG",
            location: "21:35201927",
            cdna: "-",
            cds: "-",
            ..UNKNOWN_BASE_ROW
        },
    ] {
        let label = format!("{} {} {}", columns.name, columns.location, columns.cdna);
        assert!(
            !in_the_insertion_range(columns, reference, candidate),
            "{label}"
        );
    }
    for (reference, candidate) in [
        (UNKNOWN_BESIDE_REGION, "splice_region_variant,stop_lost"),
        (UNKNOWN_BESIDE_REGION, "stop_retained_variant"),
        (
            UNKNOWN_BESIDE_REGION,
            "inframe_insertion,splice_region_variant,stop_retained_variant",
        ),
        ("splice_region_variant,stop_gained", KEPT_BESIDE_REGION),
        ("splice_region_variant", KEPT_BESIDE_REGION),
    ] {
        assert!(
            !in_the_insertion_range(UNKNOWN_BASE_ROW, reference, candidate),
            "{reference} {candidate}"
        );
    }
}

const GENETIC_CODE_26: &str = "pachysolen_nuclear_code_read_as_standard";

const GENETIC_CODE_26_COLUMNS: [&str; 9] = [
    "shape",
    "reference_terms",
    "candidate_terms",
    "feature_type",
    "reference_amino_acids",
    "candidate_amino_acids",
    "reference_codons",
    "candidate_codons",
    "taken",
];

/// A row on one codon of a coding sequence under an identifier, read to its Codons column.
fn coded_row(feature_type: &str, consequence: &str, amino_acids: &str, codons: &str) -> PairRow {
    pair_row_with_codons(&[
        "1993671",
        "7:1000",
        "A",
        "GENE1",
        "ENST00000000001",
        feature_type,
        consequence,
        "114",
        "114",
        "38",
        amino_acids,
        codons,
    ])
}

/// Each row of the table against its verdict, and under no other rule of the registry.
#[test]
fn the_genetic_code_26_class_takes_the_rows_of_its_table_marked_taken_and_no_other() {
    let rows = table!("genetic_code_26_rows.tsv", GENETIC_CODE_26_COLUMNS);
    assert_eq!(rows.len(), 24);
    for row in rows {
        let e = coded_row(row[3], row[1], row[4], row[6]);
        let a = coded_row(row[3], row[2], row[5], row[7]);
        assert_eq!(taken_by(GENETIC_CODE_26, &e, &a), yes(row[8]), "{}", row[0]);
        assert_eq!(
            pachysolen_nuclear_code_read_as_standard_pair(&e, &a),
            yes(row[8]),
            "{}",
            row[0]
        );
        for rule in row_pair_rules() {
            if rule.bucket != GENETIC_CODE_26 {
                let (e, a) = (std::slice::from_ref(&e), std::slice::from_ref(&a));
                assert!(!rule.matches_rows(e, a), "{}: {}", row[0], rule.bucket);
            }
        }
    }
}

/// A term the reference's row alone holds beside `synonymous_variant` keeps the pair a
/// difference.
#[test]
fn the_genetic_code_26_class_refuses_a_further_term_on_the_reference_row() {
    let candidate = coded_row("Transcript", "missense_variant", "A/L", "ctG/ctA");
    let alone = coded_row("Transcript", "synonymous_variant", "L", "ctG/ctA");
    assert!(taken_by(GENETIC_CODE_26, &alone, &candidate));
    let beside_a_further_term = coded_row(
        "Transcript",
        "splice_region_variant,synonymous_variant",
        "L",
        "ctG/ctA",
    );
    assert!(!taken_by(
        GENETIC_CODE_26,
        &beside_a_further_term,
        &candidate
    ));
}

/// The same row read to its ninth column and to its twelfth: the first carries neither
/// cell, and the class takes the pair of the second reading alone.
#[test]
fn a_row_read_to_its_ninth_column_carries_no_amino_acids_and_no_codons() {
    let columns = |consequence, amino_acids| {
        [
            "1993671",
            "7:1000",
            "A",
            "GENE1",
            "ENST00000000001",
            "Transcript",
            consequence,
            "114",
            "114",
            "38",
            amino_acids,
            "ctG/ctA",
        ]
    };
    let (reference, candidate) = (
        columns("synonymous_variant", "L"),
        columns("missense_variant", "A/L"),
    );
    let nine = |columns: &[&'static str; 12]| {
        let first: &[&str; 9] = columns.first_chunk().expect("nine of twelve columns");
        pair_row(first)
    };
    let (e, a) = (nine(&reference), nine(&candidate));
    assert_eq!(
        (e.amino_acids.as_deref(), e.codons.as_deref()),
        (None, None)
    );
    assert!(!pachysolen_nuclear_code_read_as_standard_pair(&e, &a));
    let (wide_e, wide_a) = (
        pair_row_with_codons(&reference),
        pair_row_with_codons(&candidate),
    );
    assert_eq!(
        wide_e,
        PairRow {
            protein_position: Some("38".to_string()),
            amino_acids: Some("L".to_string()),
            codons: Some("ctG/ctA".to_string()),
            ..e
        }
    );
    assert!(pachysolen_nuclear_code_read_as_standard_pair(
        &wide_e, &wide_a
    ));
}
