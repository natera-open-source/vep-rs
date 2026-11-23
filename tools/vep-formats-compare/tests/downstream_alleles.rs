// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of a record's later alleles, whose printed coordinates an earlier allele's
//! protein notation left: the plugin values read beside them, the three columns the shift
//! moves, the BLOSUM62 score read at the moved span, and the NMD across two sets of terms.

mod downstream_common;

use downstream_common::*;
use std::collections::BTreeMap;
use std::path::Path;
use vep_formats_compare::classes::START_COOCCURRENCE_SWAP;
use vep_formats_compare::entry::Entry;
use vep_formats_compare::parse::text::parse_default;
use vep_formats_compare::rules::downstream::hgvs_rotated_allele;
use vep_formats_compare::rules::downstream::nmd_reads_shifted_pair;
use vep_formats_compare::rules::downstream::shifted_pair_readings;
use vep_formats_compare::rules::downstream::DownstreamSequences;
use vep_formats_compare::rules::downstream::BLOSUM62_MATRIX;
use vep_formats_compare::rules::downstream::BLOSUM62_RESIDUES;
use vep_formats_compare::rules::registry::registry_class;
use vep_formats_compare::rules::Side;
use vep_formats_compare::shared::ABSENT;

/// A made-up record of two alternates on the toy forward transcript: one base, or two,
/// inserted after coding base 4, inside codon 2 and before three more copies of the base,
/// so each allele's notation moves three bases, into codon 3. The NMD record inserts after
/// coding base 99 of the three-exon transcript.
const LATER_ALLELE_RECORD: &str = "21\t1004\tv\tA\tAA,AAA\t.\t.\t.";
const LATER_ALLELE_NMD_RECORD: &str = "21\t5099\tn\tC\tCG,CGG\t.\t.\t.";

/// One allele's row of the made-up record: its coding columns, an HGVS_OFFSET where
/// `offset` gives one and the Downstream plugin's two fields where `downstream` gives the
/// protein.
#[derive(Clone)]
struct AlleleRow {
    allele: String,
    cds: String,
    protein: String,
    amino_acids: String,
    downstream: Option<String>,
    offset: Option<i64>,
    codons: String,
}

/// The cells of one row a case changes.
#[derive(Clone, Copy, Default)]
struct Change {
    cds: Option<&'static str>,
    protein: Option<&'static str>,
    amino_acids: Option<&'static str>,
    downstream: Option<&'static str>,
    offset: Option<Option<i64>>,
    codons: Option<&'static str>,
}

/// What a case changes of the four rows of the record: the first allele's and the later
/// allele's, on each engine's side.
#[derive(Clone, Copy, Default)]
struct Changes {
    reference_first: Change,
    candidate_first: Change,
    reference_later: Change,
    candidate_later: Change,
}

impl AlleleRow {
    fn of(allele: &str, cds: &str, protein: &str, amino_acids: &str) -> AlleleRow {
        AlleleRow {
            allele: allele.to_string(),
            cds: cds.to_string(),
            protein: protein.to_string(),
            amino_acids: amino_acids.to_string(),
            downstream: None,
            offset: Some(3),
            codons: "-".to_string(),
        }
    }

    fn downstream(mut self, downstream: &str) -> Self {
        self.downstream = Some(downstream.to_string());
        self
    }

    fn offset(mut self, offset: Option<i64>) -> Self {
        self.offset = offset;
        self
    }

    fn changed(mut self, change: &Change) -> Self {
        let set = |cell: &mut String, value: Option<&str>| {
            if let Some(value) = value {
                *cell = value.to_string();
            }
        };
        set(&mut self.cds, change.cds);
        set(&mut self.protein, change.protein);
        set(&mut self.amino_acids, change.amino_acids);
        set(&mut self.codons, change.codons);
        if let Some(downstream) = change.downstream {
            self.downstream = Some(downstream.to_string());
        }
        if let Some(offset) = change.offset {
            self.offset = offset;
        }
        self
    }

    /// The default-format frameshift row on the toy forward transcript.
    fn line(&self) -> String {
        let mut extra = vec!["IMPACT=HIGH".to_string(), "STRAND=1".to_string()];
        extra.extend(self.offset.map(|offset| format!("HGVS_OFFSET={offset}")));
        if let Some(downstream) = &self.downstream {
            extra.push(format!("DownstreamProtein={downstream}"));
            extra.push("ProteinLengthChange=6".to_string());
        }
        with_coding_cells(
            &default_row(
                "v",
                "21:1004-1005",
                &self.allele,
                TOY_FORWARD,
                "frameshift_variant",
                &extra.join(";"),
            ),
            &format!(
                "-\t{}\t{}\t{}\t{}\t-",
                self.cds, self.protein, self.amino_acids, self.codons
            ),
        )
    }
}

fn first_allele() -> AlleleRow {
    AlleleRow::of("A", "4-5", "2", "K/KX")
}

/// The reference's and the candidate's rows of the made-up record's two alleles. The first
/// allele's rows print one pair, the later allele's the traced shape: the reference the
/// record's span moved by the offset, the candidate the record's own. `plugin` adds each
/// row's Downstream fields, read at the shifted pair by the reference.
fn later_allele_texts(plugin: bool, changes: &Changes) -> (String, String) {
    let with = |row: AlleleRow, downstream: &str| {
        if plugin {
            row.downstream(downstream)
        } else {
            row
        }
    };
    let first = with(first_allele(), "NSCCPWLRAQIA");
    let own_first = with(first_allele(), "KNSCCPWLRAQIA");
    let later = with(AlleleRow::of("AA", "7-8", "3", "T/KX"), "KLLLSVVKGPNS");
    let own_later = with(AlleleRow::of("AA", "4-5", "2", "K/KX"), "KKLLLSVVKGPNS");
    (
        default_text(&[
            first.changed(&changes.reference_first).line(),
            later.changed(&changes.reference_later).line(),
        ]),
        default_text(&[
            own_first.changed(&changes.candidate_first).line(),
            own_later.changed(&changes.candidate_later).line(),
        ]),
    )
}

const BOTH_RECORDS: [&str; 2] = [LATER_ALLELE_RECORD, LATER_ALLELE_NMD_RECORD];

/// The rows of the NMD record's two alleles on one engine's side: the reference writes no
/// value and the candidate the plugin's.
fn later_nmd_text(
    candidate: bool,
    first_cds: &str,
    first_offset: i64,
    later_cds: &str,
    later_offset: Option<i64>,
) -> String {
    let nmd = candidate.then_some(NMD_VALUE);
    let row = |allele: &str, cds: &str, offset: Option<i64>| {
        nmd_row(
            "21:5099-5100",
            allele,
            THREE_EXONS,
            cds,
            nmd,
            "frameshift_variant",
            offset,
        )
        .replacen("v\t", "n\t", 1)
    };
    default_text(&[
        row("G", first_cds, Some(first_offset)),
        row("GG", later_cds, later_offset),
    ])
}

/// The alleles of a record share one TranscriptVariation, whose CDS pair the first allele's
/// protein notation leaves at the shifted position (TranscriptVariationAllele.pm:844-858,
/// :2404-2410), so a later allele's row prints that pair as its CDS_position
/// (OutputFactory.pm:1746) where the candidate prints the record's own. The plugins read
/// the pair the later allele's own notation leaves, the span moved by its HGVS_OFFSET. The
/// class takes the row's DownstreamProtein, and its NMD, with the candidate's pair as the
/// own pair. The same two pairs on a first allele's row, a reference pair moved by another
/// size, a candidate pair that is not the record's own, a row whose record's alleles no
/// entry states, a record without an offset on any row and a value that is not its rule's
/// reading stay differences.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_reads_the_records_own_pair_on_a_later_alleles_row() {
    let directory = scratch_of("later-allele");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let no_cache = Side::default();
    let classes = [SHIFTED_CDS_COORDINATES, SHARED_SHIFT];
    let both_rows = count_map(&[(SHIFTED_CDS_COORDINATES, 2)]);
    let first_alone = count_map(&[(SHIFTED_CDS_COORDINATES, 1)]);
    let downstream_taken = |changes: &Changes| {
        let (reference, candidate) = later_allele_texts(true, changes);
        let (report, _) = compare_record(&reference, &candidate, &side, &classes, &BOTH_RECORDS);
        masked(&report, "DownstreamProtein")
    };
    let (reference, candidate) = later_allele_texts(true, &Changes::default());
    let (report, documented) =
        compare_record(&reference, &candidate, &side, &classes, &BOTH_RECORDS);
    assert_eq!(masked(&report, "DownstreamProtein"), both_rows);
    assert_eq!(verdict_of(&report), "green");
    assert_eq!(
        counted(&documented.applied),
        count_map(&[(SHIFTED_CDS_COORDINATES, 2), (SHARED_SHIFT, 1)])
    );
    // The first allele's row alone is the class's in each cell the shape does not reach.
    let (no_alleles, _) = compare_record(&reference, &candidate, &side, &classes, &[]);
    assert_eq!(masked(&no_alleles, "DownstreamProtein"), first_alone);
    let (without_cache, _) =
        compare_record(&reference, &candidate, &no_cache, &classes, &BOTH_RECORDS);
    assert!(masked(&without_cache, "DownstreamProtein").is_empty());
    let cds = |cds: &'static str| Change {
        cds: Some(cds),
        ..Change::default()
    };
    let downstream = |downstream: &'static str| Change {
        downstream: Some(downstream),
        ..Change::default()
    };
    let refused = [
        (
            "a pair moved by another size",
            Changes {
                reference_later: cds("6-7"),
                ..Changes::default()
            },
        ),
        (
            "a candidate pair that is not the record's own",
            Changes {
                candidate_later: cds("5-6"),
                ..Changes::default()
            },
        ),
        (
            "a value that is not the rule's reading",
            Changes {
                reference_later: downstream("LLLSVVKGPNS"),
                ..Changes::default()
            },
        ),
        (
            "a value planted on the candidate's side",
            Changes {
                candidate_later: downstream("KKLLLSVVKGPN"),
                ..Changes::default()
            },
        ),
    ];
    for (label, changes) in &refused {
        assert_eq!(downstream_taken(changes), first_alone, "{label}");
    }
    // A later allele's row without an offset reads the pair the first allele's notation
    // left, with the allele as written, which this record's reference value is; without an
    // offset on any row of the record nothing is the class's.
    let no_offset = Change {
        offset: Some(None),
        ..Change::default()
    };
    assert_eq!(
        downstream_taken(&Changes {
            reference_later: no_offset,
            candidate_later: no_offset,
            ..Changes::default()
        }),
        both_rows
    );
    assert!(downstream_taken(&Changes {
        reference_first: no_offset,
        candidate_first: no_offset,
        reference_later: no_offset,
        candidate_later: no_offset,
    })
    .is_empty());
    // The two pairs on the first allele's row: that row stays a difference, the later one
    // is taken.
    assert_eq!(
        downstream_taken(&Changes {
            reference_first: cds("7-8"),
            ..Changes::default()
        }),
        first_alone
    );
    // The NMD of a later allele: the candidate's pair ends at coding base 99, the
    // reference's, three bases on, at 102, past the first 100 coding bases.
    let nmd_rows = |first_cds: (&str, &str), later_cds: (&str, &str)| {
        (
            later_nmd_text(false, first_cds.0, 3, later_cds.0, Some(3)),
            later_nmd_text(true, first_cds.1, 3, later_cds.1, Some(3)),
        )
    };
    let later_nmd_taken = |texts: &(String, String), inputs: &[&str]| {
        let (report, _) = compare_record(&texts.0, &texts.1, &side, &classes, inputs);
        masked(&report, "NMD")
    };
    let own = ("99-100", "99-100");
    let traced = nmd_rows(own, ("102-103", "99-100"));
    assert_eq!(
        later_nmd_taken(&traced, &[LATER_ALLELE_NMD_RECORD]),
        both_rows
    );
    assert_eq!(later_nmd_taken(&traced, &[]), first_alone);
    assert_eq!(
        later_nmd_taken(
            &nmd_rows(own, ("101-102", "99-100")),
            &[LATER_ALLELE_NMD_RECORD]
        ),
        first_alone
    );
    assert_eq!(
        later_nmd_taken(
            &nmd_rows(own, ("102-103", "98-99")),
            &[LATER_ALLELE_NMD_RECORD]
        ),
        first_alone
    );
    assert_eq!(
        later_nmd_taken(
            &nmd_rows(("102-103", "99-100"), own),
            &[LATER_ALLELE_NMD_RECORD]
        ),
        first_alone
    );
}

/// One cell of the fabricated grid of a later allele whose printed CDS pair an earlier
/// allele's notation left: its name, the record's alternates in order, each with the
/// HGVS_OFFSET both engines print on its row, the pair the reference prints on the last
/// allele's row, the reading its DownstreamProtein is (the span moved by so many bases,
/// the allele turned by so many), and whether the class takes the value. The record
/// inserts after coding base 4 of the toy forward transcript; the last alternate is the
/// row asked.
type EarlierShiftCell = (
    &'static str,
    &'static [(&'static str, Option<i64>)],
    &'static str,
    (i64, i64),
    bool,
);

const EARLIER_SHIFT_GRID: [EarlierShiftCell; 12] = [
    (
        "the pair an earlier allele left, the value at the row's own offset",
        &[("A", Some(3)), ("AC", Some(2))],
        "7-8",
        (2, 2),
        true,
    ),
    (
        "the pair moved by the row's own offset",
        &[("A", Some(3)), ("AC", Some(2))],
        "6-7",
        (2, 2),
        true,
    ),
    (
        "an allele without an offset between the two",
        &[("A", Some(3)), ("AG", None), ("AC", Some(2))],
        "7-8",
        (2, 2),
        true,
    ),
    (
        "the pair of the last earlier allele that has an offset",
        &[("A", Some(1)), ("AG", Some(3)), ("AC", Some(2))],
        "7-8",
        (2, 2),
        true,
    ),
    (
        "the pair of the first allele where a later one has another",
        &[("A", Some(1)), ("AG", Some(3)), ("AC", Some(2))],
        "5-6",
        (2, 2),
        false,
    ),
    (
        "a pair moved by neither offset",
        &[("A", Some(3)), ("AC", Some(2))],
        "5-6",
        (2, 2),
        false,
    ),
    (
        "the value read at the printed pair",
        &[("A", Some(3)), ("AC", Some(2))],
        "7-8",
        (3, 0),
        false,
    ),
    (
        "the first allele's row",
        &[("AC", Some(2))],
        "7-8",
        (2, 2),
        false,
    ),
    (
        "no offset on the row, the pair an earlier allele left, the allele as written",
        &[("A", Some(3)), ("AC", None)],
        "7-8",
        (3, 0),
        true,
    ),
    (
        "no offset on the row and none on an earlier allele's",
        &[("A", None), ("AC", None)],
        "7-8",
        (3, 0),
        false,
    ),
    (
        "no offset on the row, another value",
        &[("A", Some(3)), ("AC", None)],
        "7-8",
        (2, 2),
        false,
    ),
    (
        "no offset on the row, the record's own pair printed",
        &[("A", Some(3)), ("AC", None)],
        "4-5",
        (3, 0),
        false,
    ),
];

fn one_entry(row: &str) -> Entry {
    let mut entries = parse_default(&default_text(&[row])).expect("the row parses");
    assert_eq!(entries.len(), 1);
    entries.remove(0)
}

/// The DownstreamProtein of an insertion of `allele` after coding base 4 of the toy
/// forward transcript, read at its span moved by `moved_by` bases with the allele turned
/// as a shift of `turned_by` turns it: Ensembl VEP's rule, or the candidate's.
fn earlier_shift_reading(
    sequences: &DownstreamSequences,
    allele: &str,
    moved_by: i64,
    turned_by: i64,
    ensembl: bool,
) -> String {
    let e = one_entry(&AlleleRow::of(allele, "4-5", "2", "K/KX").line());
    let transcript = sequences
        .transcript(&e)
        .expect("the cache reads")
        .expect("the cache holds the transcript");
    let mut pair = (5, 4);
    if moved_by != 0 {
        let model = transcript
            .model
            .as_ref()
            .expect("the transcript has a model");
        let readings = shifted_pair_readings(&e, model, moved_by, &[(Some(5), Some(4))]);
        assert_eq!(readings.len(), 1);
        pair = readings[0].1.expect("the moved span is coding");
    }
    transcript
        .downstream_value(
            pair.0,
            pair.1,
            &hgvs_rotated_allele(allele, turned_by, 1),
            false,
            ensembl,
        )
        .expect("the pair has a value")
}

/// The reference's and the candidate's rows of a cell's record, and its input line. The
/// rows of the earlier alleles are alike on both sides; the last allele's row prints the
/// cell's pair and value on the reference's side and the record's own pair with the
/// candidate's rule on its own.
fn earlier_shift_texts(
    sequences: &DownstreamSequences,
    cell: &EarlierShiftCell,
) -> (String, String, String) {
    let (_, alternates, printed, (moved_by, turned_by), _) = *cell;
    let (allele, offset) = alternates[alternates.len() - 1];
    let earlier: Vec<String> = alternates[..alternates.len() - 1]
        .iter()
        .map(|(other, shift)| {
            AlleleRow::of(other, "4-5", "2", "K/KX")
                .downstream("KNS")
                .offset(*shift)
                .line()
        })
        .collect();
    let asked_reference = AlleleRow::of(allele, printed, "3", "T/KX")
        .downstream(&earlier_shift_reading(
            sequences, allele, moved_by, turned_by, true,
        ))
        .offset(offset)
        .line();
    let asked_candidate = AlleleRow::of(allele, "4-5", "2", "K/KX")
        .downstream(&earlier_shift_reading(sequences, allele, 0, 0, false))
        .offset(offset)
        .line();
    let alleles: Vec<String> = alternates
        .iter()
        .map(|(other, _)| format!("A{other}"))
        .collect();
    let record = format!("21\t1004\tv\tA\t{}\t.\t.\t.", alleles.join(","));
    let text = |asked: String| {
        let mut rows = earlier.clone();
        rows.push(asked);
        default_text(&rows)
    };
    (text(asked_reference), text(asked_candidate), record)
}

/// Each allele's protein notation moves the CDS pair on the TranscriptVariation the
/// record's alleles share (TranscriptVariationAllele.pm:844-858), and a row's CDS_position
/// is printed before its own notation runs (OutputFactory.pm:1746, :1769). So a later
/// allele's row prints the record's span moved by the offset of the last earlier allele
/// that has one, while the plugin reads the pair at the row's own offset with the allele
/// turned by it, or, on a row whose notation moves nothing, the printed pair with the
/// allele as written. Over the fabricated grid the class takes the DownstreamProtein
/// exactly there. The readings the grid tells apart differ in value.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_reads_the_pair_an_earlier_allele_left_and_nothing_wider() {
    let directory = scratch_of("earlier-shift-grid");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let no_cache = Side::default();
    let sequences = sequences_of(&side);
    let readings: Vec<String> = [(0, 0), (2, 2), (3, 0)]
        .iter()
        .map(|(moved_by, turned_by)| {
            earlier_shift_reading(sequences, "AC", *moved_by, *turned_by, true)
        })
        .collect();
    assert!(readings[0] != readings[1] && readings[0] != readings[2] && readings[1] != readings[2]);
    let once = count_map(&[(SHIFTED_CDS_COORDINATES, 1)]);
    for cell in &EARLIER_SHIFT_GRID {
        let (reference, candidate, record) = earlier_shift_texts(sequences, cell);
        let (report, _) = compare_record(
            &reference,
            &candidate,
            &side,
            &[SHIFTED_CDS_COORDINATES],
            &[record.as_str()],
        );
        let taken = masked(&report, "DownstreamProtein");
        assert_eq!(
            differing(&report, "DownstreamProtein") + taken.values().sum::<i64>(),
            1,
            "{}",
            cell.0
        );
        assert_eq!(taken == once, cell.4, "{}", cell.0);
        // Without the record's alleles in order, and without the cache, no cell is the
        // class's.
        let (no_alleles, _) = compare_record(
            &reference,
            &candidate,
            &side,
            &[SHIFTED_CDS_COORDINATES],
            &[],
        );
        assert!(
            masked(&no_alleles, "DownstreamProtein").is_empty(),
            "{}",
            cell.0
        );
        let (without_cache, _) = compare_record(
            &reference,
            &candidate,
            &no_cache,
            &[SHIFTED_CDS_COORDINATES],
            &[record.as_str()],
        );
        assert!(
            masked(&without_cache, "DownstreamProtein").is_empty(),
            "{}",
            cell.0
        );
    }
}

/// The NMD of a later allele whose own notation moves it three bases while the first
/// allele's moved four: the reference prints the span four bases on and its plugin reads
/// the pair three bases on, past the first 100 coding bases, where the candidate's own pair
/// ends at base 99. The class takes it with the candidate's pair as the own pair. The same
/// row without an offset of its own stays a difference: the class reads an NMD at the row's
/// own offset alone.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_takes_an_nmd_at_the_rows_own_offset_beside_an_earlier_alleles_pair(
) {
    let directory = scratch_of("later-allele-nmd");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let classes = [SHIFTED_CDS_COORDINATES, SHARED_SHIFT];
    let taken = |later_offset: Option<i64>, later_cds: &str| {
        let reference = later_nmd_text(false, "99-100", 4, later_cds, later_offset);
        let candidate = later_nmd_text(true, "99-100", 4, "99-100", later_offset);
        let (report, _) = compare_record(
            &reference,
            &candidate,
            &side,
            &classes,
            &[LATER_ALLELE_NMD_RECORD],
        );
        masked(&report, "NMD")
    };
    let both_rows = count_map(&[(SHIFTED_CDS_COORDINATES, 2)]);
    let first_alone = count_map(&[(SHIFTED_CDS_COORDINATES, 1)]);
    assert_eq!(taken(Some(3), "103-104"), both_rows);
    assert_eq!(taken(Some(3), "102-103"), both_rows);
    assert_eq!(taken(Some(3), "101-102"), first_alone);
    assert_eq!(taken(None, "103-104"), first_alone);
}

/// A later allele's row prints the CDS pair and the translation pair the first allele's
/// protein notation left on the TranscriptVariation the record's alleles share
/// (OutputFactory.pm:1746, :1750; TranscriptVariationAllele.pm:1724-1733, :844-858), and
/// the shared reference allele's peptide read again there (:1780-1783, :1793). The class
/// takes the row's CDS_position, Protein_position and Amino_acids where the reference's
/// pair is the record's span moved by the first allele's HGVS_OFFSET and the candidate's
/// the span itself. A first allele's row, a first allele without an offset or with another,
/// a move of another size, a candidate that does not print the own pair, a column outside
/// the class, a protein position or a peptide that is not its pair's, and a row without the
/// cache or the record's alleles stay differences.
#[test]
fn hgvs_shift_shared_across_alleles_takes_three_columns_of_a_later_alleles_row() {
    let directory = scratch_of("shared-shift");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let no_cache = Side::default();
    let columns = ["CDS_position", "Protein_position", "Amino_acids"];
    let once = count_map(&[(SHARED_SHIFT, 1)]);
    let taken = |read: &Side, inputs: &[&str], changes: &Changes| -> Vec<&str> {
        let (reference, candidate) = later_allele_texts(false, changes);
        let (report, documented) =
            compare_record(&reference, &candidate, read, &[SHARED_SHIFT], inputs);
        let accepted: Vec<&str> = columns
            .iter()
            .copied()
            .filter(|column| masked(&report, column) == once)
            .collect();
        assert_eq!(
            verdict_of(&report) == "green",
            accepted == columns && changes.reference_later.codons.is_none()
        );
        let applied: &[(&str, i64)] = if accepted.is_empty() {
            &[]
        } else {
            &[(SHARED_SHIFT, 1)]
        };
        assert_eq!(counted(&documented.applied), count_map(applied));
        accepted
    };
    let with_cache = |changes: &Changes| taken(&side, &[LATER_ALLELE_RECORD], changes);
    assert_eq!(with_cache(&Changes::default()), columns);
    let moved = Change {
        cds: Some("7-8"),
        protein: Some("3"),
        amino_acids: Some("T/KX"),
        ..Change::default()
    };
    let own = Change {
        cds: Some("4-5"),
        protein: Some("2"),
        amino_acids: Some("K/KX"),
        ..Change::default()
    };
    let offset = |offset: Option<i64>| Change {
        offset: Some(offset),
        ..Change::default()
    };
    let cds = |cds: &'static str| Change {
        cds: Some(cds),
        ..Change::default()
    };
    let refused = [
        (
            "the moved values on a first allele's row",
            Changes {
                reference_first: moved,
                reference_later: own,
                ..Changes::default()
            },
        ),
        (
            "a first allele without an offset",
            Changes {
                reference_first: offset(None),
                candidate_first: offset(None),
                ..Changes::default()
            },
        ),
        (
            "a first allele that moved two bases",
            Changes {
                reference_first: offset(Some(2)),
                candidate_first: offset(Some(2)),
                ..Changes::default()
            },
        ),
        (
            "a move of another size",
            Changes {
                reference_later: cds("6-7"),
                ..Changes::default()
            },
        ),
        (
            "a candidate that does not print the own pair",
            Changes {
                candidate_later: cds("5-6"),
                ..Changes::default()
            },
        ),
    ];
    for (label, changes) in &refused {
        assert!(with_cache(changes).is_empty(), "{label}");
    }
    assert!(taken(&no_cache, &[LATER_ALLELE_RECORD], &Changes::default()).is_empty());
    assert!(taken(&side, &[], &Changes::default()).is_empty());
    // A column outside the class stays a difference beside the three.
    let reference_later = |change: Change| Changes {
        reference_later: change,
        ..Changes::default()
    };
    assert_eq!(
        with_cache(&reference_later(Change {
            codons: Some("aCt/aAAct"),
            ..Change::default()
        })),
        columns
    );
    // A protein position that is not the moved pair's codon leaves the peptide open as
    // well; a peptide with another alternate part, or whose reference part is not the moved
    // codon's residue, stays open alone.
    assert_eq!(
        with_cache(&reference_later(Change {
            protein: Some("4"),
            ..Change::default()
        })),
        ["CDS_position"]
    );
    let amino_acids = |amino_acids: &'static str| Change {
        amino_acids: Some(amino_acids),
        ..Change::default()
    };
    assert_eq!(
        with_cache(&reference_later(amino_acids("T/TX"))),
        columns[..2]
    );
    assert_eq!(
        with_cache(&reference_later(amino_acids("L/KX"))),
        columns[..2]
    );
    assert_eq!(
        with_cache(&Changes {
            candidate_later: amino_acids("L/KX"),
            ..Changes::default()
        }),
        columns[..2]
    );
}

/// The last allele's row of a cell of the BLOSUM62 grid on one engine's side: its CDS
/// position, protein position, Amino_acids and BLOSUM62.
type Blosum62Row = (
    &'static str,
    &'static str,
    &'static str,
    Option<&'static str>,
);

/// One cell of the fabricated grid of the class: a record of the toy forward transcript
/// that inserts after coding base 6, the last base of codon 2. Its name, the record's
/// alternates in order with the HGVS_OFFSET on each one's row, the last allele's row on
/// the reference's side and on the candidate's, and whether the class takes the value. The
/// last alternate is the row asked: three bases between two codons, whose own pair has no
/// reference residue.
type Blosum62Cell = (
    &'static str,
    &'static [(&'static str, Option<i64>)],
    Blosum62Row,
    Blosum62Row,
    bool,
);

const BLOSUM62_OWN: Blosum62Row = ("6-7", "2-3", "-/K", None);
const BLOSUM62_MOVED: Blosum62Row = ("7-8", "3", "T/K", Some("-1"));
const BLOSUM62_GRID: [Blosum62Cell; 10] = [
    (
        "a later allele read one base on, the first allele's offset",
        &[("C", Some(1)), ("AAA", None)],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        true,
    ),
    (
        "a third allele, the first allele's offset",
        &[("C", Some(1)), ("G", None), ("AAA", None)],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        true,
    ),
    (
        "the first allele's row",
        &[("AAA", Some(1))],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        false,
    ),
    (
        "a first allele without an offset",
        &[("C", None), ("AAA", None)],
        BLOSUM62_MOVED,
        BLOSUM62_OWN,
        false,
    ),
    (
        "a value that is not the matrix cell of its pair",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "T/K", Some("5")),
        BLOSUM62_OWN,
        false,
    ),
    (
        "a reference residue that is not the moved codon's",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "S/K", Some("0")),
        BLOSUM62_OWN,
        false,
    ),
    (
        "alternate residues that differ",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "T/R", Some("-1")),
        BLOSUM62_OWN,
        false,
    ),
    (
        "a candidate pair the plugin scores",
        &[("C", Some(1)), ("AAA", None)],
        BLOSUM62_MOVED,
        ("6-7", "2-3", "N/K", None),
        false,
    ),
    (
        "a pair moved by another size",
        &[("C", Some(1)), ("AAA", None)],
        ("8-9", "3", "T/K", Some("-1")),
        BLOSUM62_OWN,
        false,
    ),
    (
        "the value on the candidate's side alone",
        &[("C", Some(1)), ("AAA", None)],
        ("7-8", "3", "T/K", None),
        ("6-7", "2-3", "-/K", Some("-1")),
        false,
    ),
];

fn blosum62_row(allele: &str, row: Blosum62Row, offset: Option<i64>) -> String {
    let (cds, protein, amino_acids, blosum) = row;
    let mut extra = vec!["IMPACT=MODERATE".to_string(), "STRAND=1".to_string()];
    extra.extend(offset.map(|offset| format!("HGVS_OFFSET={offset}")));
    extra.extend(blosum.map(|blosum| format!("BLOSUM62={blosum}")));
    with_coding_cells(
        &default_row(
            "b",
            "21:1006-1007",
            allele,
            TOY_FORWARD,
            "inframe_insertion",
            &extra.join(";"),
        ),
        &format!("-\t{cds}\t{protein}\t{amino_acids}\t-\t-"),
    )
}

/// The reference's and the candidate's rows of a cell's record, alike on the earlier
/// alleles' rows, and its input line.
fn blosum62_texts(cell: &Blosum62Cell) -> (String, String, String) {
    let (_, alternates, reference, candidate, _) = *cell;
    let earlier: Vec<String> = alternates[..alternates.len() - 1]
        .iter()
        .map(|(other, shift)| blosum62_row(other, ("6-7", "2-3", "-/X", None), *shift))
        .collect();
    let (allele, offset) = alternates[alternates.len() - 1];
    let alleles: Vec<String> = alternates
        .iter()
        .map(|(other, _)| format!("A{other}"))
        .collect();
    let record = format!("21\t1006\tb\tA\t{}\t.\t.\t.", alleles.join(","));
    let text = |row: Blosum62Row| {
        let mut rows = earlier.clone();
        rows.push(blosum62_row(allele, row, offset));
        default_text(&rows)
    };
    (text(reference), text(candidate), record)
}

/// Blosum62.pm:110 scores a row whose peptide pair is one residue over one. A later
/// allele's row carries the coordinates the first allele's notation left on the object the
/// record's alleles share, so Ensembl VEP reads the reference residue one codon on and
/// scores an insertion between two codons, whose own pair has no reference residue. Over
/// the fabricated grid the class takes the BLOSUM62 the reference alone writes exactly
/// where the row's CDS and protein positions are the base class's, the reference's pair is
/// the moved codon's residue over the candidate's own alternate residue, the candidate's
/// pair is none the plugin scores, and the value is the matrix cell of the reference's
/// pair. A run without the cache and one without the record's alleles stay differences.
#[test]
fn blosum62_hgvs_shift_shared_across_alleles_takes_the_score_read_at_the_moved_span_and_nothing_wider(
) {
    let directory = scratch_of("blosum62-grid");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let no_cache = Side::default();
    assert_eq!(
        [('T', 'K'), ('S', 'K'), ('T', 'R')].map(|pair| BLOSUM62_MATRIX[&pair]),
        ["-1", "0", "-1"]
    );
    let once = count_map(&[(BLOSUM62_SHARED_SHIFT, 1)]);
    for cell in &BLOSUM62_GRID {
        let (reference, candidate, record) = blosum62_texts(cell);
        let (report, _) = compare_record(
            &reference,
            &candidate,
            &side,
            &[BLOSUM62_SHARED_SHIFT],
            &[record.as_str()],
        );
        let taken = masked(&report, "BLOSUM62");
        assert_eq!(
            differing(&report, "BLOSUM62") + taken.values().sum::<i64>(),
            1,
            "{}",
            cell.0
        );
        assert_eq!(taken == once, cell.4, "{}", cell.0);
        let (without_cache, _) = compare_record(
            &reference,
            &candidate,
            &no_cache,
            &[BLOSUM62_SHARED_SHIFT],
            &[record.as_str()],
        );
        assert!(masked(&without_cache, "BLOSUM62").is_empty(), "{}", cell.0);
        let (no_alleles, _) =
            compare_record(&reference, &candidate, &side, &[BLOSUM62_SHARED_SHIFT], &[]);
        assert!(masked(&no_alleles, "BLOSUM62").is_empty(), "{}", cell.0);
    }
}

/// The comparator's matrix is the plugin's own (Blosum62.pm:50-71, residues in the order
/// of :73), cell for cell the one the plugin of `crates/vep-builtins` holds.
#[test]
fn blosum62_matrix_is_the_one_the_plugin_scores_with() {
    let source = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../crates/vep-builtins/src/plugins/blosum62.rs"),
    )
    .expect("the plugin's module reads");
    let residues = regex::Regex::new(r#"RESIDUES: &\[u8; 20\] = b"([A-Z]{20})";"#)
        .unwrap()
        .captures(&source)
        .expect("the module names its residues")[1]
        .to_string();
    let opening = "const BLOSUM62: [[i8; 20]; 20] = [";
    let start = source.find(opening).expect("the module holds its matrix") + opening.len();
    let end = start + source[start..].find("\n];").expect("the matrix ends");
    let cells: Vec<&str> = regex::Regex::new(r"-?\d+")
        .unwrap()
        .find_iter(&source[start..end])
        .map(|cell| cell.as_str())
        .collect();
    assert_eq!(residues, BLOSUM62_RESIDUES);
    assert_eq!(cells.len(), 400);
    let mut held = BTreeMap::new();
    for (i, first) in residues.chars().enumerate() {
        for (j, second) in residues.chars().enumerate() {
            held.insert((first, second), cells[20 * i + j]);
        }
    }
    assert_eq!(held, *BLOSUM62_MATRIX);
    assert!(held
        .keys()
        .all(|(first, second)| BLOSUM62_MATRIX[&(*first, *second)]
            == BLOSUM62_MATRIX[&(*second, *first)]));
}

/// The two pairs of consequence sets of a row on which Ensembl VEP writes
/// stop_retained_variant for a stop codon test it does not make and the candidate
/// coding_sequence_variant, each set with a term of the NMD plugin's.
const RETAINED_STOP_PAIRS: [(&str, &str); 2] = [
    (
        "intron_variant,splice_acceptor_variant,splice_donor_5th_base_variant,splice_donor_variant,stop_retained_variant",
        "coding_sequence_variant,intron_variant,splice_acceptor_variant,splice_donor_5th_base_variant,splice_donor_variant",
    ),
    (
        "splice_acceptor_variant,stop_retained_variant",
        "coding_sequence_variant,splice_acceptor_variant",
    ),
];

type NmdRow = (&'static str, &'static str, &'static str, &'static str);

const NMD_ACCEPTOR_ROW: NmdRow = ("21:3013-3014", "-", TOY_TWO_EXONS, "?-10");
const NMD_PENULTIMATE_ROW: NmdRow = ("21:5200-5201", "-", THREE_EXONS, "?-111");
const NMD_LATE_ROW: NmdRow = ("21:5100", "-", THREE_EXONS, "100");

/// The row an NMD question is asked of: its cells, its HGVS_OFFSET, the value each engine
/// writes (`None` for no value) and the cDNA_position both print.
#[derive(Clone, Copy)]
struct NmdAsk {
    row: NmdRow,
    offset: Option<i64>,
    vep: Option<&'static str>,
    vep_rs: Option<&'static str>,
    cdna: &'static str,
}

impl Default for NmdAsk {
    fn default() -> Self {
        NmdAsk {
            row: NMD_ACCEPTOR_ROW,
            offset: Some(1),
            vep: Some(NMD_VALUE),
            vep_rs: None,
            cdna: "-",
        }
    }
}

/// Whether `plugin_hgvs_shifted_cds_coordinates` takes the NMD of the asked row that the
/// reference's entry writes under `reference_terms` and the candidate's under
/// `candidate_terms`: the comparison's verdict on the field, held to the rule's own return.
fn nmd_taken(reference_terms: &str, candidate_terms: &str, side: &Side, ask: NmdAsk) -> bool {
    let (location, allele, feature, cds) = ask.row;
    let text = |value: Option<&str>, terms: &str| {
        default_text(&[
            nmd_row(location, allele, feature, cds, value, terms, ask.offset).replacen(
                &format!("\t-\t{cds}\t"),
                &format!("\t{}\t{cds}\t", ask.cdna),
                1,
            ),
        ])
    };
    let reference = text(ask.vep, reference_terms);
    let candidate = text(ask.vep_rs, candidate_terms);
    let (report, documented) =
        compare_downstream(&reference, &candidate, side, &[SHIFTED_CDS_COORDINATES]);
    let once = count_map(&[(SHIFTED_CDS_COORDINATES, 1)]);
    let applied = counted(&documented.applied);
    let taken = applied == once;
    assert_eq!(taken, masked(&report, "NMD") == once);
    assert!(taken || applied.is_empty());
    let (e, a) = (
        one_entry(reference.lines().nth(1).expect("the reference has a row")),
        one_entry(candidate.lines().nth(1).expect("the candidate has a row")),
    );
    let by_rule = side.downstream_sequences.as_ref().is_some_and(|sequences| {
        nmd_reads_shifted_pair(
            &e,
            &a,
            ask.vep.unwrap_or(ABSENT),
            ask.vep_rs.unwrap_or(ABSENT),
            sequences,
            None,
        )
        .expect("the cache reads")
    });
    assert_eq!(taken, by_rule);
    taken
}

/// Ensembl VEP writes stop_retained_variant on a deletion's row whose stop codon test it
/// does not make, where the candidate writes coding_sequence_variant, so the two entries of
/// a row the class reads can differ by that pair of terms. The acceptor deletions of the
/// two-exon and the three-exon transcripts have no CDS start, so the test is made on
/// neither, and under both pairs of sets such rows carry the class takes the value the
/// reference alone writes at the shifted coordinates. The Consequence stays a difference,
/// since the class takes the plugin's field alone, and without the cache it takes nothing.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_takes_an_nmd_across_the_retained_stop_term_pair() {
    let directory = scratch_of("nmd-term-pair");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let no_cache = Side::default();
    for (reference_terms, candidate_terms) in RETAINED_STOP_PAIRS {
        for row in [NMD_ACCEPTOR_ROW, NMD_PENULTIMATE_ROW] {
            let ask = NmdAsk {
                row,
                ..NmdAsk::default()
            };
            assert!(
                nmd_taken(reference_terms, candidate_terms, &side, ask),
                "{reference_terms} {row:?}"
            );
            assert!(
                !nmd_taken(reference_terms, candidate_terms, &no_cache, ask),
                "{reference_terms} {row:?}"
            );
        }
    }
    // Three rows in one comparison: the two acceptor rows under the pair and the deletion
    // of coding base 100 under one set.
    let (reference_terms, candidate_terms) = RETAINED_STOP_PAIRS[0];
    let rows = |acceptor: Option<&str>, terms: &str, late: Option<&str>| {
        let line = |row: NmdRow, nmd: Option<&str>, terms: &str, offset: i64| {
            nmd_row(row.0, row.1, row.2, row.3, nmd, terms, Some(offset))
        };
        default_text(&[
            line(NMD_ACCEPTOR_ROW, acceptor, terms, 1),
            line(NMD_PENULTIMATE_ROW, acceptor, terms, 1),
            line(NMD_LATE_ROW, late, "frameshift_variant", 3),
        ])
    };
    let reference = rows(Some(NMD_VALUE), reference_terms, None);
    let candidate = rows(None, candidate_terms, Some(NMD_VALUE));
    let (report, documented) =
        compare_downstream(&reference, &candidate, &side, &[SHIFTED_CDS_COORDINATES]);
    let all_three = count_map(&[(SHIFTED_CDS_COORDINATES, 3)]);
    assert_eq!(counted(&documented.applied), all_three);
    assert_eq!(masked(&report, "NMD"), all_three);
    assert_eq!(verdict_of(&report), "red");
    assert_eq!(differing(&report, "Consequence"), 2);
    assert_eq!(differing(&report, "NMD"), 0);
}

/// Two differing consequence sets pass the class's NMD test only as the pair
/// `stop_retained_without_stop_codon_test_pair` admits, each set naming one of the
/// plugin's terms. Refused on the acceptor row the class takes under that pair: the two
/// terms the other way round, stop_lost in coding_sequence_variant's place, a further term
/// on either side, the registry's `start_cooccurrence_swap` pair, which the Downstream arm
/// of the class admits, and the pair where a set names none of the plugin's terms. The pair
/// on a row with a cDNA position and a CDS start, where both engines make the stop codon
/// test, is refused where the same row without the cDNA position is taken.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_refuses_an_nmd_across_a_term_pair_the_retained_stop_pair_does_not_admit(
) {
    let directory = scratch_of("nmd-term-pair-refused");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    let (retained, coding) = RETAINED_STOP_PAIRS[1];
    let asked = NmdAsk::default();
    assert!(nmd_taken(retained, coding, &side, asked));
    assert!(nmd_taken(
        retained,
        "frameshift_variant,splice_acceptor_variant",
        &side,
        asked
    ));
    let coding_and_intron = format!("{coding},intron_variant");
    let retained_and_intron = format!("{retained},intron_variant");
    let start_pair = (
        "splice_acceptor_variant,start_lost,start_retained_variant",
        "splice_acceptor_variant,start_retained_variant",
    );
    let refused = [
        ("the two terms the other way round", (coding, retained)),
        (
            "stop_lost in the term's place",
            (retained, "splice_acceptor_variant,stop_lost"),
        ),
        (
            "a further term on the candidate's side",
            (retained, coding_and_intron.as_str()),
        ),
        (
            "a further term on the reference's side",
            (retained_and_intron.as_str(), coding),
        ),
        ("the registry's start pair", start_pair),
        (
            "no term of the plugin's in the reference's set",
            ("stop_retained_variant", "frameshift_variant"),
        ),
        (
            "no term of the plugin's in either set",
            (
                "intron_variant,stop_retained_variant",
                "coding_sequence_variant,intron_variant",
            ),
        ),
    ];
    for (label, (reference_terms, candidate_terms)) in refused {
        assert!(
            !nmd_taken(reference_terms, candidate_terms, &side, asked),
            "{label}"
        );
    }
    assert_eq!(
        registry_class(start_pair.0, start_pair.1, "-", &[START_COOCCURRENCE_SWAP]).as_deref(),
        Some(START_COOCCURRENCE_SWAP)
    );
    // The deletion of the three-exon transcript's coding base 100, whose value the
    // candidate alone writes: without a cDNA position the row shows a stop codon test that
    // is not made, and with one beside its CDS start a test both engines make.
    let late = NmdAsk {
        row: NMD_LATE_ROW,
        offset: Some(3),
        vep: None,
        vep_rs: Some(NMD_VALUE),
        cdna: "-",
    };
    let terms = (
        "frameshift_variant,stop_retained_variant",
        "frameshift_variant",
    );
    assert!(nmd_taken(terms.0, terms.1, &side, late));
    assert!(!nmd_taken(
        terms.0,
        terms.1,
        &side,
        NmdAsk {
            cdna: "300",
            ..late
        }
    ));
}

/// The class's NMD test on two entries of one consequence set is the test it makes across
/// the retained-stop pair, less the pair: on the two acceptor rows the class takes the
/// value under either set of each pair carried by both entries as it does across the pair,
/// a set without one of the plugin's terms is refused on both entries as across the pair,
/// and a row without an HGVS_OFFSET or with a shift that leaves the exon is refused both
/// ways.
#[test]
fn plugin_hgvs_shifted_cds_coordinates_reads_an_nmd_under_one_consequence_set_as_across_the_pair() {
    let directory = scratch_of("nmd-one-set");
    let side = transcripts_of(&downstream_cache(&directory.join("cache")));
    for row in [NMD_ACCEPTOR_ROW, NMD_PENULTIMATE_ROW] {
        let asked = NmdAsk {
            row,
            ..NmdAsk::default()
        };
        for pair in RETAINED_STOP_PAIRS {
            for terms in [pair.0, pair.1] {
                assert!(nmd_taken(terms, terms, &side, asked), "{terms} {row:?}");
            }
            assert!(nmd_taken(pair.0, pair.1, &side, asked), "{pair:?} {row:?}");
            for one_set in [(pair.0, pair.0), (pair.1, pair.1), pair] {
                for offset in [None, Some(300)] {
                    assert!(
                        !nmd_taken(one_set.0, one_set.1, &side, NmdAsk { offset, ..asked }),
                        "{one_set:?} {row:?} {offset:?}"
                    );
                }
            }
        }
    }
    for terms in [
        "intron_variant,stop_retained_variant",
        "coding_sequence_variant,intron_variant",
        "inframe_deletion",
    ] {
        assert!(
            !nmd_taken(terms, terms, &side, NmdAsk::default()),
            "{terms}"
        );
    }
}
