// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The coding terms vep-rs owes on a row at a coding end (`expected_terms`), the cells
//! derived beside them, the classes of pair (`CodingEndClass`) and the place a landed
//! rule holds beside a cache (`landed_rule_may_take`), each on made-up cache records.
//!
//! A record is one exon whose first cDNA base on the forward strand, or last on the
//! reverse, is the genomic base `FIRST`, unless a test gives it mapper pairs of its own:
//! `leader` cDNA bases of 5' UTR, the coding sequence, `trailer` bases of 3' UTR. A row's
//! cells are Ensembl VEP's on the record, written by hand.

use std::collections::BTreeSet;
use vep_formats_compare::coding_ends::terms::{
    amino_acids_on_a_completed_stop_codon, class_of_pair, coding_family,
    codons_on_a_completed_stop_codon, count_alone_differs_at_the_first_residue, expected_terms,
    frameshift_counted_at_the_first_residue, frameshift_the_coding_end_cuts_short,
    landed_rule_may_take, outside_the_coding_family, position_cell,
    protein_notation_on_a_completed_stop_codon, protein_notation_where_no_start_codon_is_shown,
    Cells, CodingEndClass, CodingTerm, LandedRule, NotDerived, RowFacts, WhatWasRead,
    STRUCTURAL_CLASS_WORDS,
};
use vep_formats_compare::coding_ends::{
    coding_ends, coding_ends_with_three_prime_utr, CodingEnds, MapperFacts, RecordFacts,
};

const FIRST: i64 = 1001;

/// Codons of methionine, alanine, glycine, glutamate and serine, and a stop codon of the
/// standard and the vertebrate mitochondrial code.
const MET: &str = "ATG";
const ALA: &str = "GCT";
const GLY: &str = "GGA";
const GLU: &str = "GAA";
const SER: &str = "TCA";
const STOP: &str = "TAA";

type Pair = (i64, i64, i64, i64, i64);

#[derive(Clone)]
struct Made {
    codes: Vec<&'static str>,
    consensus_id: Option<&'static str>,
    codon_table: Option<i64>,
    coding: String,
    peptide: Option<&'static str>,
    strand: i64,
    start_phase: i64,
    leader: i64,
    trailer: i64,
    pairs: Option<Vec<Pair>>,
    span: Option<(i64, i64)>,
    introns: Vec<(i64, i64)>,
    /// The cached 3' UTR; without one, as many cytosines as `trailer`.
    three_prime_utr: Option<&'static str>,
    /// Whether the reader hands the record's 3' UTR to the reading.
    utr_read: bool,
}

/// A record of the forward strand whose coding sequence is `pieces` joined, without a
/// cached peptide, with six bases of 5' UTR and `trailer` bases of 3' UTR.
fn made(pieces: &[&str], trailer: i64) -> Made {
    Made {
        codes: Vec::new(),
        consensus_id: None,
        codon_table: None,
        coding: pieces.concat(),
        peptide: None,
        strand: 1,
        start_phase: 0,
        leader: 6,
        trailer,
        pairs: None,
        span: None,
        introns: Vec::new(),
        three_prime_utr: None,
        utr_read: true,
    }
}

impl Made {
    fn cdna_length(&self) -> i64 {
        self.leader + self.coding.len() as i64 - self.start_phase.max(0) + self.trailer
    }

    /// The genomic base of a cDNA position of the one exon, on past either end.
    fn genomic(&self, cdna: i64) -> i64 {
        match self.strand {
            1 => FIRST + cdna - 1,
            _ => FIRST + self.cdna_length() - cdna,
        }
    }

    /// The Location of the cDNA positions `first` to `last`, lowest base first.
    fn at(&self, first: i64, last: i64) -> (i64, i64) {
        let (one, other) = (self.genomic(first), self.genomic(last));
        (one.min(other), one.max(other))
    }

    fn read(&self) -> CodingEnds {
        let on_the_cdna = self.coding.len() as i64 - self.start_phase.max(0);
        let last = self.cdna_length();
        let one_exon = [(1, last, FIRST, FIRST + last - 1, self.strand)];
        let cytosines = "C".repeat(self.trailer as usize);
        let three_prime_utr = self.three_prime_utr.unwrap_or(&cytosines);
        let facts = RecordFacts {
            attribute_codes: &self.codes,
            consensus_id: self.consensus_id,
            strand: Some(self.strand),
            codon_table: self.codon_table,
            coding: Some(&self.coding),
            peptide: self.peptide,
            mapper: Some(MapperFacts {
                start_phase: self.start_phase,
                coding_start: self.leader + 1,
                coding_end: self.leader + on_the_cdna,
                pairs: self.pairs.as_deref().unwrap_or(&one_exon),
            }),
            span: self.span.or(Some((FIRST, FIRST + last - 1))),
            introns: &self.introns,
        };
        match self.utr_read {
            true => coding_ends_with_three_prime_utr(&facts, Some(three_prime_utr)),
            false => coding_ends(&facts),
        }
    }
}

/// Whole codons whose last names a residue, with a 3' UTR: no stop codon stands.
fn open_whole() -> Made {
    made(&[MET, ALA, GLY], 6)
}

/// `open_whole` with `cds_end_NF`.
fn flagged() -> Made {
    Made {
        codes: vec!["cds_end_NF"],
        ..open_whole()
    }
}

/// Whole codons closed by a stop codon, with a 3' UTR.
fn closed() -> Made {
    made(&[MET, ALA, GLY, STOP], 6)
}

/// A coding sequence that ends `tail` bases into a stop codon on the transcript's last
/// base, under the vertebrate mitochondrial code: polyadenylation completes the codon.
fn shown_tail(tail: usize) -> Made {
    Made {
        codon_table: Some(2),
        ..made(&[MET, ALA, GLY, &STOP[..tail]], 0)
    }
}

/// The same end under the standard code, shown complete by a consensus id.
fn shown_by_an_id(tail: usize) -> Made {
    Made {
        codon_table: None,
        consensus_id: Some("CCDS1.1"),
        ..shown_tail(tail)
    }
}

/// The same end under the standard code without an id: the record does not show the
/// model complete, and no stop codon stands.
fn tail_not_shown() -> Made {
    Made {
        codon_table: None,
        ..shown_tail(2)
    }
}

/// A first codon that is no start codon, closed by a stop codon: no start codon is
/// shown, and a stop codon stands.
fn start_open() -> Made {
    made(&[GLU, ALA, GLY, STOP], 6)
}

/// Neither a start codon nor a stop codon.
fn both_open() -> Made {
    made(&[GLU, ALA, GLY], 6)
}

/// One pair's row and what a comparison read of it.
#[derive(Clone)]
struct Row {
    allele: String,
    location: (i64, i64),
    reference_length: Option<u64>,
    theirs: &'static str,
    cells: Option<[&'static str; 5]>,
    ours: Option<[&'static str; 5]>,
    transcript_row: bool,
    reverse_record: bool,
}

/// A row of `allele` at `location` that stands in for `reference_length` bases, with
/// Ensembl VEP's terms `theirs` and its cells: `cDNA_position`, `CDS_position`,
/// `Protein_position`, `Amino_acids` and `Codons`, which vep-rs's row prints alike.
fn row(
    allele: &str,
    location: (i64, i64),
    reference_length: u64,
    theirs: &'static str,
    cells: [&'static str; 5],
) -> Row {
    Row {
        allele: allele.to_string(),
        location,
        reference_length: Some(reference_length),
        theirs,
        cells: Some(cells),
        ours: None,
        transcript_row: true,
        reverse_record: false,
    }
}

const BLANK: [&str; 5] = ["-", "-", "-", "-", "-"];

fn cells_of(cells: &[&'static str; 5]) -> Cells<'static> {
    Cells {
        cdna_position: cells[0],
        cds_position: cells[1],
        protein_position: cells[2],
        amino_acids: cells[3],
        codons: cells[4],
    }
}

type Reference<'a> = Option<&'a dyn Fn(i64) -> Option<u8>>;

fn with<T>(
    row: &Row,
    reference: Reference<'_>,
    ask: impl FnOnce(&RowFacts<'_>, &WhatWasRead<'_>) -> T,
) -> T {
    let facts = RowFacts {
        transcript_row: row.transcript_row,
        allele: &row.allele,
        location: row.location,
        ensembl_vep_terms: row.theirs,
    };
    let read = WhatWasRead {
        cells: row.cells.as_ref().map(|cells| {
            (
                cells_of(cells),
                cells_of(row.ours.as_ref().unwrap_or(cells)),
            )
        }),
        reference_length: row.reference_length,
        reference_bases: reference,
        record_on_the_reverse_strand: row.reverse_record,
    };
    ask(&facts, &read)
}

fn names(terms: &BTreeSet<CodingTerm>) -> String {
    terms
        .iter()
        .map(|term| term.name())
        .collect::<Vec<_>>()
        .join(",")
}

/// The derived terms of the row on the record, joined by commas in the order vep-rs
/// writes them.
fn owed(record: &Made, row: &Row) -> Result<String, NotDerived> {
    owed_beside(record, row, None)
}

fn owed_beside(record: &Made, row: &Row, reference: Reference<'_>) -> Result<String, NotDerived> {
    let ends = record.read();
    with(row, reference, |facts, read| {
        expected_terms(facts, Some(&ends), read).map(|terms| names(&terms))
    })
}

fn ok(terms: &str) -> Result<String, NotDerived> {
    Ok(terms.to_string())
}

/// The class that takes the pair of the row and vep-rs's terms `ours`.
fn class(record: &Made, row: &Row, ours: &str) -> Option<CodingEndClass> {
    let ends = record.read();
    with(row, None, |facts, read| {
        class_of_pair(facts, ours, Some(&ends), read)
    })
}

/// The reference bases of a genome that holds the record's cDNA on the forward strand
/// from `FIRST` on, cytosines for its 5' UTR, and `beyond` after its last base.
fn genome(record: &Made, beyond: &'static str) -> impl Fn(i64) -> Option<u8> {
    let mut bases = vec![b'C'; record.leader as usize];
    bases.extend_from_slice(record.coding.as_bytes());
    bases.extend_from_slice(beyond.as_bytes());
    move |position: i64| {
        usize::try_from(position - FIRST)
            .ok()
            .and_then(|place| bases.get(place).copied())
    }
}

// The stop end: no stop codon stands.

#[test]
fn a_deletion_in_the_last_whole_codon_where_no_stop_codon_stands_owes_a_frameshift_alone() {
    let deletion = |record: &Made| {
        row(
            "-",
            record.at(14, 14),
            1,
            "frameshift_variant,stop_lost",
            ["14", "8", "3", "G/X", "gGa/ga"],
        )
    };
    for record in [open_whole(), tail_not_shown()] {
        assert_eq!(owed(&record, &deletion(&record)), ok("frameshift_variant"));
        assert_eq!(
            class(&record, &deletion(&record), "frameshift_variant"),
            Some(CodingEndClass::StopCodonOverlap)
        );
    }
    // On the reverse strand the same row lies at the base that strand gives it.
    let reverse = Made {
        strand: -1,
        ..open_whole()
    };
    assert_eq!(
        owed(&reverse, &deletion(&reverse)),
        ok("frameshift_variant")
    );
}

#[test]
fn frameshift_alone_where_the_inserted_bases_hold_a_stop_codon() {
    let record = open_whole();
    let insertion = row(
        "CTTAA",
        record.at(13, 14),
        0,
        "frameshift_variant,stop_lost",
        ["13-14", "7-8", "3", "G/A*X", "gga/gCTTAAga"],
    );
    assert_eq!(
        owed(&record, &insertion),
        ok("stop_gained,frameshift_variant")
    );
    assert_eq!(class(&record, &insertion, "frameshift_variant"), None);
    assert_eq!(
        class(&record, &insertion, "stop_gained,frameshift_variant"),
        Some(CodingEndClass::StopCodonOverlap)
    );
}

#[test]
fn a_stop_codon_of_the_utr_at_the_old_coding_end_keeps_no_stop() {
    let record = open_whole();
    // One base deleted, and one inserted, in the last whole codon: Ensembl VEP reads a
    // stop codon kept and an in-frame term behind it.
    let deletion = row(
        "-",
        record.at(14, 14),
        1,
        "inframe_deletion,stop_retained_variant",
        ["14", "8", "3", "G/X", "gGa/ga"],
    );
    let insertion = row(
        "C",
        record.at(14, 15),
        0,
        "inframe_insertion,stop_retained_variant",
        ["14-15", "8-9", "3", "G/GX", "gga/ggCa"],
    );
    for kept in [&deletion, &insertion] {
        assert_eq!(owed(&record, kept), ok("frameshift_variant"));
        assert_eq!(
            class(&record, kept, "frameshift_variant"),
            Some(CodingEndClass::StopCodonOverlap)
        );
    }
    // The last whole codon deleted before a stop codon of the 3' UTR.
    let whole_codon = row(
        "-",
        record.at(13, 15),
        3,
        "inframe_deletion,stop_retained_variant",
        ["13-15", "7-9", "3", "G/-", "GGA/-"],
    );
    assert_eq!(owed(&record, &whole_codon), ok("inframe_deletion"));
    assert_eq!(
        class(&record, &whole_codon, "inframe_deletion"),
        Some(CodingEndClass::StopCodonOverlap)
    );
}

#[test]
fn a_deletion_from_the_coding_sequence_past_its_end_owes_the_coding_term() {
    let record = open_whole();
    // Into the 3' UTR, and past the transcript's last base.
    let into_the_utr = row(
        "-",
        record.at(14, 17),
        4,
        "stop_lost,3_prime_UTR_variant",
        ["14-17", "8-?", "3-?", "-", "-"],
    );
    let past_the_end = row(
        "-",
        record.at(14, 30),
        17,
        "stop_retained_variant,3_prime_UTR_variant",
        ["14-?", "8-?", "3-?", "-", "-"],
    );
    for span in [&into_the_utr, &past_the_end] {
        assert_eq!(owed(&record, span), ok("coding_sequence_variant"));
        assert_eq!(
            class(&record, span, "coding_sequence_variant,3_prime_UTR_variant"),
            Some(CodingEndClass::StopCodonOverlap)
        );
    }
}

/// A record of two exons, the cDNA positions 1 to 9 on the bases 1001 to 1009 and 10 to
/// 21 on 1101 to 1112, whose coding sequence is `open_whole`'s.
fn two_exons() -> Made {
    Made {
        pairs: Some(vec![(1, 9, 1001, 1009, 1), (10, 21, 1101, 1112, 1)]),
        span: Some((1001, 1112)),
        introns: vec![(1010, 1100)],
        ..open_whole()
    }
}

#[test]
fn deletion_from_the_intron_before_the_last_coding_exon() {
    let record = two_exons();
    // Its 3' end in the last whole codon: the cells lack their starts.
    let into_the_codon = row(
        "-",
        (1050, 1105),
        56,
        "stop_retained_variant,splice_acceptor_variant,intron_variant",
        ["?-14", "?-8", "?-3", "-", "-"],
    );
    assert_eq!(
        owed(&record, &into_the_codon),
        ok("coding_sequence_variant")
    );
    assert_eq!(
        class(
            &record,
            &into_the_codon,
            "coding_sequence_variant,splice_acceptor_variant,intron_variant"
        ),
        Some(CodingEndClass::StopCodonOverlap)
    );
    assert_eq!(
        class(
            &record,
            &into_the_codon,
            "frameshift_variant,coding_sequence_variant,splice_acceptor_variant,intron_variant"
        ),
        None
    );
    // From the intron past the transcript's end every cell is blank, and the span holds
    // the genomic base of a coding position.
    let past_the_end = row(
        "-",
        (1050, 1130),
        81,
        "stop_retained_variant,intron_variant",
        BLANK,
    );
    assert_eq!(owed(&record, &past_the_end), ok("coding_sequence_variant"));
    // Inside the intron it holds none.
    let in_the_intron = row("-", (1050, 1060), 11, "intron_variant", BLANK);
    assert_eq!(owed(&record, &in_the_intron), ok(""));
}

#[test]
fn an_insertion_before_the_last_three_coding_positions_owes_a_frameshift() {
    let record = open_whole();
    let insertion = row(
        "C",
        record.at(12, 13),
        0,
        "inframe_insertion,stop_retained_variant",
        ["12-13", "6-7", "2-3", "-/X", "-/C"],
    );
    assert_eq!(owed(&record, &insertion), ok("frameshift_variant"));
    assert_eq!(
        class(&record, &insertion, "frameshift_variant"),
        Some(CodingEndClass::StopCodonOverlap)
    );
}

#[test]
fn inframe_insertion_alone_where_a_stop_codon_comes_behind_the_last_residue() {
    for record in [open_whole(), flagged(), tail_not_shown()] {
        let insertion = row(
            "CTA",
            record.at(14, 15),
            0,
            "inframe_insertion,stop_retained_variant",
            ["14-15", "8-9", "3", "G/G*", "gga/ggCTAa"],
        );
        assert_eq!(
            owed(&record, &insertion),
            ok("stop_gained,inframe_insertion")
        );
        assert_eq!(class(&record, &insertion, "inframe_insertion"), None);
        assert_eq!(
            class(&record, &insertion, "stop_gained,inframe_insertion"),
            Some(CodingEndClass::StopCodonBehindLastResidue)
        );
    }
    // The substitution of as many bases that brings the stop codon the same way.
    let record = open_whole();
    let substitution = row(
        "TAA",
        record.at(16, 18),
        3,
        "stop_retained_variant,3_prime_UTR_variant",
        ["16-18", "-", "-", "-", "-"],
    );
    assert_eq!(owed(&record, &substitution), ok(""));
}

#[test]
fn stop_term_of_a_stop_codon_inside_the_coding_sequence_dropped() {
    // The codon before the last is a stop codon; three bases deleted over its last two
    // and the last codon's first leave a codon of serine, or a stop codon.
    let lost = made(&[MET, "TAG", ALA], 6);
    let deletion = row(
        "-",
        lost.at(11, 13),
        3,
        "stop_lost,inframe_deletion",
        ["11-13", "5-7", "2-3", "*A/S", "tAGGct/tct"],
    );
    assert_eq!(owed(&lost, &deletion), ok("stop_lost,inframe_deletion"));
    assert_eq!(class(&lost, &deletion, "inframe_deletion"), None);
    let kept = made(&[MET, "TAG", "CAA"], 6);
    let deletion = row(
        "-",
        kept.at(11, 13),
        3,
        "inframe_deletion,stop_retained_variant",
        ["11-13", "5-7", "2-3", "*Q/*", "tAGCaa/taa"],
    );
    assert_eq!(
        owed(&kept, &deletion),
        ok("inframe_deletion,stop_retained_variant")
    );
    assert_eq!(class(&kept, &deletion, "inframe_deletion"), None);
}

#[test]
fn a_substitution_in_the_last_whole_codon_owes_the_terms_of_its_residues() {
    let record = open_whole();
    let substitution = row(
        "C",
        record.at(14, 14),
        1,
        "missense_variant",
        ["14", "8", "3", "G/A", "gGa/gCa"],
    );
    assert_eq!(owed(&record, &substitution), ok("missense_variant"));
    assert_eq!(class(&record, &substitution, "missense_variant"), None);
}

// The stop end: polyadenylation completes the stop codon.

const INCOMPLETE: &str = "incomplete_terminal_codon_variant,coding_sequence_variant";

#[test]
fn a_substitution_in_the_tail_loses_or_keeps_the_completed_stop_codon() {
    // One tail base under the mitochondrial code: any other base there names a residue.
    let record = shown_tail(1);
    let lost = row(
        "C",
        record.at(16, 16),
        1,
        INCOMPLETE,
        ["16", "10", "4", "X", "T/C"],
    );
    assert_eq!(owed(&record, &lost), ok("stop_lost"));
    assert_eq!(
        class(&record, &lost, "stop_lost"),
        Some(CodingEndClass::StopCodonCompleted)
    );
    let cell = |derive: fn(
        &RowFacts<'_>,
        Option<&CodingEnds>,
        &WhatWasRead<'_>,
    ) -> Result<String, NotDerived>,
                record: &Made,
                row: &Row| {
        let ends = record.read();
        with(row, None, |facts, read| derive(facts, Some(&ends), read))
    };
    assert_eq!(
        cell(amino_acids_on_a_completed_stop_codon, &record, &lost),
        ok("*/Q")
    );
    assert_eq!(
        cell(codons_on_a_completed_stop_codon, &record, &lost),
        ok("Taa/Caa")
    );
    assert_eq!(
        cell(protein_notation_on_a_completed_stop_codon, &record, &lost),
        ok("p.Ter4GlnextTer?")
    );
    // Two tail bases: the second base changed leaves a stop codon of the standard code,
    // which is a codon of tryptophan under the mitochondrial one.
    for (record, owes) in [
        (shown_by_an_id(2), "stop_retained_variant"),
        (shown_tail(2), "stop_lost"),
    ] {
        let second = row(
            "G",
            record.at(17, 17),
            1,
            INCOMPLETE,
            ["17", "11", "4", "X", "ta/tG"],
        );
        assert_eq!(owed(&record, &second), ok(owes));
        assert_eq!(
            class(&record, &second, owes),
            Some(CodingEndClass::StopCodonCompleted)
        );
        // On the reverse strand the cell holds the complement.
        let reverse = Made {
            strand: -1,
            ..record.clone()
        };
        let second = row(
            "C",
            reverse.at(17, 17),
            1,
            INCOMPLETE,
            ["17", "11", "4", "X", "ta/tG"],
        );
        assert_eq!(owed(&reverse, &second), ok(owes));
    }
}

#[test]
fn two_base_allele_on_the_reverse_strand() {
    // The two tail bases replaced by a base pair that keeps the stop codon of the
    // standard code: the cell holds it reversed and complemented.
    let record = Made {
        strand: -1,
        ..shown_by_an_id(2)
    };
    let substitution = row(
        "CA",
        record.at(16, 17),
        2,
        INCOMPLETE,
        ["16-17", "10-11", "4", "X", "ta/TG"],
    );
    assert_eq!(owed(&record, &substitution), ok("stop_retained_variant"));
    assert_eq!(class(&record, &substitution, "stop_lost"), None);
    assert_eq!(
        class(&record, &substitution, "stop_retained_variant"),
        Some(CodingEndClass::StopCodonCompleted)
    );
    // A record given on the reverse strand is not derived where the bases are read.
    let given_reversed = Row {
        reverse_record: true,
        ..substitution
    };
    assert_eq!(
        owed(&record, &given_reversed),
        Err(NotDerived::RecordOnTheReverseStrand)
    );
}

#[test]
fn a_substitution_of_several_bases_into_the_tail_loses_the_completed_stop_codon() {
    let record = shown_tail(1);
    let notation = |row: &Row| {
        let ends = record.read();
        with(row, None, |facts, read| {
            protein_notation_on_a_completed_stop_codon(facts, Some(&ends), read)
        })
    };
    // The last residue kept, and changed.
    let kept = row(
        "CC",
        record.at(15, 16),
        2,
        "missense_variant,coding_sequence_variant",
        ["15-16", "9-10", "3-4", "GX", "ggAT/ggCC"],
    );
    assert_eq!(owed(&record, &kept), ok("stop_lost"));
    assert_eq!(notation(&kept), ok("p.Ter4GlnextTer?"));
    let changed = row(
        "ACC",
        record.at(14, 16),
        3,
        "missense_variant,coding_sequence_variant",
        ["14-16", "8-10", "3-4", "GX", "gGAT/gACC"],
    );
    assert_eq!(owed(&record, &changed), ok("stop_lost"));
    assert_eq!(notation(&changed), ok("p.Gly3_Ter4delinsAspGln"));
    assert_eq!(
        class(&record, &changed, "stop_lost"),
        Some(CodingEndClass::StopCodonCompleted)
    );
}

#[test]
fn stop_lost_beside_coding_sequence_variant_on_a_completed_stop_codon() {
    // A deletion in the tail: the adenines that follow name a residue.
    let record = shown_tail(1);
    let deletion = row(
        "-",
        record.at(16, 16),
        1,
        "protein_altering_variant,incomplete_terminal_codon_variant",
        ["16", "10", "4", "X/-", "T/-"],
    );
    assert_eq!(owed(&record, &deletion), ok("stop_lost"));
    assert_eq!(
        class(&record, &deletion, "stop_lost,coding_sequence_variant"),
        None
    );
    assert_eq!(
        class(
            &record,
            &deletion,
            "stop_lost,incomplete_terminal_codon_variant"
        ),
        None
    );
    assert_eq!(
        class(&record, &deletion, "stop_lost"),
        Some(CodingEndClass::StopCodonCompleted)
    );
    let ends = record.read();
    with(&deletion, None, |facts, read| {
        assert_eq!(
            amino_acids_on_a_completed_stop_codon(facts, Some(&ends), read),
            ok("*/X")
        );
        assert_eq!(
            codons_on_a_completed_stop_codon(facts, Some(&ends), read),
            ok("Taa/aa")
        );
        assert_eq!(
            protein_notation_on_a_completed_stop_codon(facts, Some(&ends), read),
            ok("p.Ter4LysextTer?")
        );
    });
    // The second of two tail bases deleted: the adenine after it leaves the codon.
    let record = shown_tail(2);
    let deletion = row(
        "-",
        record.at(17, 17),
        1,
        INCOMPLETE,
        ["17", "11", "4", "X/-", "tA/t"],
    );
    assert_eq!(owed(&record, &deletion), ok("stop_retained_variant"));
}

#[test]
fn stop_codon_inserted_before_the_tail() {
    let record = shown_tail(1);
    let before_the_tail = |bases: &str| {
        row(
            bases,
            record.at(15, 16),
            0,
            "inframe_insertion,incomplete_terminal_codon_variant,coding_sequence_variant",
            ["15-16", "9-10", "3-4", "-/X", "-/N"],
        )
    };
    for (bases, owes) in [
        ("C", "frameshift_variant"),
        ("CCC", "inframe_insertion"),
        ("TAA", "inframe_insertion,stop_retained_variant"),
        ("TAAC", "inframe_insertion,stop_retained_variant"),
    ] {
        let insertion = before_the_tail(bases);
        assert_eq!(owed(&record, &insertion), ok(owes), "{}", bases.len());
        assert_eq!(
            class(&record, &insertion, owes),
            Some(CodingEndClass::StopCodonCompleted),
            "{}",
            bases.len()
        );
    }
    assert_eq!(
        class(&record, &before_the_tail("TAA"), "inframe_insertion"),
        None
    );
    assert_eq!(
        class(&record, &before_the_tail("TAAC"), "frameshift_variant"),
        None
    );
    let ends = record.read();
    with(&before_the_tail("C"), None, |facts, read| {
        assert_eq!(
            protein_notation_on_a_completed_stop_codon(facts, Some(&ends), read),
            ok("p.Ter4LeuextTer?")
        );
    });
}

#[test]
fn an_insertion_between_the_two_tail_bases_loses_or_keeps_the_stop_codon() {
    for (record, owes) in [
        (shown_tail(2), "stop_lost"),
        (shown_by_an_id(2), "stop_retained_variant"),
    ] {
        let insertion = row(
            "G",
            record.at(16, 17),
            0,
            "stop_gained,incomplete_terminal_codon_variant,coding_sequence_variant",
            ["16-17", "10-11", "4", "X/X", "ta/tGa"],
        );
        assert_eq!(owed(&record, &insertion), ok(owes));
        // A multiple of three bases there is a shape left out.
        let in_frame = Row {
            allele: "GGG".to_string(),
            ..insertion
        };
        assert_eq!(owed(&record, &in_frame), Err(NotDerived::ShapeNotCovered));
    }
}

#[test]
fn frameshift_alone_before_a_completed_stop_codon_that_is_kept() {
    let record = shown_tail(1);
    // A deletion in the last whole codon reaches no base of the stop codon.
    let deletion = row(
        "-",
        record.at(14, 14),
        1,
        "frameshift_variant,stop_lost",
        ["14", "8", "3", "G/X", "gGa/ga"],
    );
    assert_eq!(owed(&record, &deletion), ok("frameshift_variant"));
    assert_eq!(
        class(&record, &deletion, "frameshift_variant"),
        Some(CodingEndClass::StopCodonOverlap)
    );
    let ends = record.read();
    with(&deletion, None, |facts, read| {
        assert_eq!(
            protein_notation_on_a_completed_stop_codon(facts, Some(&ends), read),
            ok("p.Gly3AspfsTer?")
        );
    });
    // Four bases inserted in the last whole codon that keep its residue and bring a
    // stop codon to the old stop position.
    for record in [
        shown_tail(1),
        Made {
            strand: -1,
            ..shown_tail(1)
        },
    ] {
        let bases = if record.strand == 1 { "CTAA" } else { "TTAG" };
        let insertion = row(
            bases,
            record.at(14, 15),
            0,
            "frameshift_variant,stop_lost",
            ["14-15", "8-9", "3", "G/GX", "gga/ggNa"],
        );
        assert_eq!(
            owed(&record, &insertion),
            ok("inframe_insertion,stop_retained_variant")
        );
        assert_eq!(class(&record, &insertion, "frameshift_variant"), None);
        // The right row has no class: Ensembl VEP's row shows no incomplete codon read,
        // and the derived set holds a stop term.
        assert_eq!(
            class(
                &record,
                &insertion,
                "inframe_insertion,stop_retained_variant"
            ),
            None
        );
    }
}

#[test]
fn three_base_deletion_into_the_tail() {
    // Two bases over the last whole codon's last and the tail, and three.
    let record = shown_tail(1);
    let two = row(
        "-",
        record.at(15, 16),
        2,
        "frameshift_variant",
        ["15-16", "9-10", "3-4", "GX/X", "ggAT/gg"],
    );
    assert_eq!(owed(&record, &two), ok("frameshift_variant,stop_lost"));
    let three = row(
        "-",
        record.at(14, 16),
        3,
        "stop_lost,inframe_deletion",
        ["14-16", "8-10", "3-4", "GX/X", "gGAT/g"],
    );
    assert_eq!(owed(&record, &three), ok("stop_lost,inframe_deletion"));
    assert_eq!(class(&record, &three, "inframe_deletion"), None);
    // Where the base left with the adenines is a stop codon, no stop is lost.
    let record = Made {
        coding: [MET, ALA, SER, "T"].concat(),
        ..shown_tail(1)
    };
    let three = row(
        "-",
        record.at(14, 16),
        3,
        "stop_lost,inframe_deletion",
        ["14-16", "8-10", "3-4", "SX/X", "tCAT/t"],
    );
    assert_eq!(owed(&record, &three), ok("inframe_deletion"));
    assert_eq!(class(&record, &three, "stop_lost,inframe_deletion"), None);
    assert_eq!(
        class(&record, &three, "inframe_deletion"),
        Some(CodingEndClass::StopCodonCompleted)
    );
}

#[test]
fn deletion_past_the_end_that_can_be_written_off_the_transcript() {
    let record = shown_tail(1);
    let from_the_tail = row(
        "-",
        record.at(16, 18),
        3,
        INCOMPLETE,
        ["16-?", "10-?", "4-?", "-", "-"],
    );
    // The bases after the deletion do not repeat its first: the adenines that follow
    // the coding bases left name a residue.
    let apart = genome(&record, "GGG");
    assert_eq!(
        owed_beside(&record, &from_the_tail, Some(&apart)),
        ok("stop_lost")
    );
    // The base after it repeats the tail base: written one place on, the deletion takes
    // no coding base.
    let repeated = genome(&record, "GGT");
    assert_eq!(
        owed_beside(&record, &from_the_tail, Some(&repeated)),
        ok("stop_retained_variant")
    );
    // The bases before it repeat its last: written two places back it ends inside the
    // transcript.
    let back = genome(&record, "GAC");
    assert_eq!(
        owed_beside(&record, &from_the_tail, Some(&back)),
        Err(NotDerived::PlacementInsideTheTranscript)
    );
    assert_eq!(
        owed(&record, &from_the_tail),
        Err(NotDerived::NoReferenceBases)
    );
    // From the last whole codon past the end, the span Ensembl VEP reads as untested.
    let from_the_codon = row(
        "-",
        record.at(15, 18),
        4,
        "stop_retained_variant",
        ["15-?", "9-?", "3-?", "-", "-"],
    );
    let apart = genome(&record, "CCC");
    assert_eq!(
        owed_beside(&record, &from_the_codon, Some(&apart)),
        ok("stop_lost")
    );
}

#[test]
fn deletion_past_the_end_that_leaves_the_first_tail_base_before_the_adenines() {
    // Two tail bases; the base after the deletion repeats the first and the next does
    // not repeat the second: written one place on, the deletion leaves the first tail
    // base, which two adenines complete.
    let record = shown_tail(2);
    let from_the_tail = row(
        "-",
        record.at(16, 19),
        4,
        INCOMPLETE,
        ["16-?", "10-?", "4-?", "-", "-"],
    );
    let repeated_once = genome(&record, "GGTC");
    assert_eq!(
        owed_beside(&record, &from_the_tail, Some(&repeated_once)),
        ok("stop_retained_variant")
    );
    let apart = genome(&record, "GGCC");
    assert_eq!(
        owed_beside(&record, &from_the_tail, Some(&apart)),
        ok("stop_lost")
    );
    // From the second tail base past the end: a shape left out.
    let from_the_second = row(
        "-",
        record.at(17, 19),
        3,
        INCOMPLETE,
        ["17-?", "11-?", "4-?", "-", "-"],
    );
    assert_eq!(
        owed_beside(&record, &from_the_second, Some(&apart)),
        Err(NotDerived::ShapeNotCovered)
    );
}

#[test]
fn a_substitution_from_the_tail_past_the_end_owes_the_coding_term() {
    let record = shown_tail(1);
    let substitution = row(
        "CC",
        record.at(16, 17),
        2,
        INCOMPLETE,
        ["16-?", "10-?", "4-?", "-", "-"],
    );
    assert_eq!(owed(&record, &substitution), ok("coding_sequence_variant"));
    assert_eq!(
        class(&record, &substitution, "coding_sequence_variant"),
        Some(CodingEndClass::StopCodonCompleted)
    );
}

// The start end.

#[test]
fn a_substitution_in_a_first_codon_that_is_no_start_codon_owes_the_terms_of_its_residues() {
    let record = start_open();
    let notation = |row: &Row| {
        let ends = record.read();
        with(row, None, |facts, read| {
            protein_notation_where_no_start_codon_is_shown(facts, Some(&ends), read, ("", ""))
        })
    };
    for (bases, place, theirs, cells, owes, written) in [
        (
            "G",
            9,
            "start_lost,synonymous_variant",
            ["9", "3", "1", "E", "gaA/gaG"],
            "synonymous_variant",
            "p.Glu1=",
        ),
        (
            "T",
            9,
            "start_lost",
            ["9", "3", "1", "E/D", "gaA/gaT"],
            "missense_variant",
            "p.Glu1Asp",
        ),
        (
            "T",
            9,
            "start_lost,start_retained_variant",
            ["9", "3", "1", "E/D", "gaA/gaT"],
            "missense_variant",
            "p.Glu1Asp",
        ),
        (
            "T",
            7,
            "stop_gained,start_lost",
            ["7", "1", "1", "E/*", "Gaa/Taa"],
            "stop_gained",
            "p.Glu1Ter",
        ),
    ] {
        let substitution = row(bases, record.at(place, place), 1, theirs, cells);
        assert_eq!(owed(&record, &substitution), ok(owes), "{theirs}");
        assert_eq!(notation(&substitution), ok(written), "{theirs}");
        assert_eq!(
            class(&record, &substitution, owes),
            Some(CodingEndClass::StartCodonOverlap),
            "{theirs}"
        );
        assert_eq!(class(&record, &substitution, ""), None, "{theirs}");
    }
}

#[test]
fn a_frameshift_in_a_first_codon_that_is_no_start_codon_owes_the_frameshift_alone() {
    let record = start_open();
    let deletion = row(
        "-",
        record.at(8, 8),
        1,
        "frameshift_variant,start_lost",
        ["8", "2", "1", "E/X", "gAa/ga"],
    );
    assert_eq!(owed(&record, &deletion), ok("frameshift_variant"));
    assert_eq!(
        class(&record, &deletion, "frameshift_variant"),
        Some(CodingEndClass::StartCodonOverlap)
    );
    // A first codon that the start exon's phase pads.
    let padded = Made {
        coding: ["NNT", ALA, GLY, STOP].concat(),
        start_phase: 2,
        ..start_open()
    };
    let deletion = row(
        "-",
        padded.at(7, 7),
        1,
        "frameshift_variant,start_lost",
        ["7", "3", "1", "X/X", "nnT/nn"],
    );
    assert_eq!(owed(&padded, &deletion), ok("frameshift_variant"));
}

#[test]
fn inframe_insertion_in_a_padded_first_codon() {
    let padded = Made {
        coding: ["NNT", ALA, GLY, STOP].concat(),
        start_phase: 2,
        ..start_open()
    };
    let insertion = row(
        "CCC",
        padded.at(7, 8),
        0,
        "start_lost",
        ["7-8", "3-4", "1-2", "X/XP", "nnt/nntCCC"],
    );
    assert_eq!(
        owed(&padded, &insertion),
        ok("inframe_insertion,coding_sequence_variant")
    );
    assert_eq!(class(&padded, &insertion, "inframe_insertion"), None);
}

#[test]
fn an_inframe_allele_at_a_first_codon_that_is_no_start_codon_owes_the_term_of_its_peptides() {
    let record = start_open();
    let deletion = row(
        "-",
        record.at(7, 9),
        3,
        "start_lost,inframe_deletion",
        ["7-9", "1-3", "1", "E/-", "GAA/-"],
    );
    assert_eq!(owed(&record, &deletion), ok("inframe_deletion"));
    let after_the_codon = row(
        "GGA",
        record.at(9, 10),
        0,
        "start_lost",
        ["9-10", "3-4", "1-2", "-/G", "-/GGA"],
    );
    assert_eq!(owed(&record, &after_the_codon), ok("inframe_insertion"));
    assert_eq!(
        class(&record, &after_the_codon, "inframe_insertion"),
        Some(CodingEndClass::StartCodonOverlap)
    );
}

#[test]
fn inframe_insertion_where_a_protein_altering_variant_is_owed() {
    let record = Made {
        leader: 0,
        ..start_open()
    };
    let insertion = row(
        "CTC",
        record.at(1, 2),
        0,
        "start_lost",
        ["1-2", "1-2", "1", "E/AQ", "gaa/gCTCaa"],
    );
    assert_eq!(owed(&record, &insertion), ok("protein_altering_variant"));
    assert_eq!(class(&record, &insertion, "inframe_insertion"), None);
    assert_eq!(
        class(&record, &insertion, "protein_altering_variant"),
        Some(CodingEndClass::StartCodonOverlap)
    );
}

#[test]
fn stop_gained_alone_where_an_inframe_insertion_is_owed_beside_it() {
    let record = start_open();
    let kept = row(
        "ATA",
        record.at(8, 9),
        0,
        "stop_gained,start_lost",
        ["8-9", "2-3", "1", "E/E*", "gaa/gaATAa"],
    );
    assert_eq!(owed(&record, &kept), ok("stop_gained,inframe_insertion"));
    assert_eq!(class(&record, &kept, "stop_gained"), None);
    assert_eq!(
        class(&record, &kept, "stop_gained,inframe_insertion"),
        Some(CodingEndClass::StartCodonOverlap)
    );
    // Where the first residue changes before the stop codon, the protein is altered.
    let changed = row(
        "CTT",
        record.at(7, 8),
        0,
        "stop_gained,start_lost",
        ["7-8", "1-2", "1", "E/A*", "gaa/gCTTaa"],
    );
    assert_eq!(
        owed(&record, &changed),
        ok("stop_gained,protein_altering_variant")
    );
}

#[test]
fn a_span_from_the_utr_into_a_first_codon_that_is_no_start_codon_owes_the_coding_term() {
    let record = start_open();
    let deletion = row(
        "-",
        record.at(5, 8),
        4,
        "start_lost,5_prime_UTR_variant",
        ["5-8", "?-2", "?-1", "-", "-"],
    );
    assert_eq!(owed(&record, &deletion), ok("coding_sequence_variant"));
    assert_eq!(
        class(
            &record,
            &deletion,
            "coding_sequence_variant,5_prime_UTR_variant"
        ),
        Some(CodingEndClass::StartCodonOverlap)
    );
    assert_eq!(class(&record, &deletion, "coding_sequence_variant"), None);
    // An insertion between the last base of the 5' UTR and the first coding base the
    // genome holds, on a model whose first codon the genome lacks.
    let edited = Made {
        pairs: Some(vec![(1, 6, 1001, 1006, 1), (10, 24, 1007, 1021, 1)]),
        span: Some((1001, 1021)),
        ..start_open()
    };
    assert!(edited.read().edited_model);
    let insertion = row(
        "CCCC",
        (1006, 1007),
        0,
        "start_lost,start_retained_variant,5_prime_UTR_variant",
        ["7-9", "?-3", "?-1", "-", "-"],
    );
    assert_eq!(owed(&edited, &insertion), ok("coding_sequence_variant"));
    assert_eq!(
        class(
            &edited,
            &insertion,
            "coding_sequence_variant,5_prime_UTR_variant"
        ),
        Some(CodingEndClass::StartCodonOverlap)
    );
}

// A structural allele.

#[test]
fn a_structural_deletion_over_a_coding_end_where_no_codon_stands_loses_none() {
    let over_the_end = |record: &Made| {
        row(
            "deletion",
            record.at(10, 40),
            31,
            "stop_lost,coding_sequence_variant,3_prime_UTR_variant,feature_truncation",
            BLANK,
        )
    };
    for record in [flagged(), open_whole(), tail_not_shown()] {
        let deletion = over_the_end(&record);
        assert_eq!(owed(&record, &deletion), ok("coding_sequence_variant"));
        assert_eq!(
            class(
                &record,
                &deletion,
                "coding_sequence_variant,3_prime_UTR_variant,feature_truncation"
            ),
            Some(CodingEndClass::StructuralDeletionStopLost)
        );
    }
    assert_eq!(
        owed(&closed(), &over_the_end(&closed())),
        ok("stop_lost,coding_sequence_variant")
    );
    // A copy number loss is read alike.
    let loss = Row {
        allele: "copy_number_loss".to_string(),
        ..over_the_end(&flagged())
    };
    assert_eq!(owed(&flagged(), &loss), ok("coding_sequence_variant"));
    // A start codon that is not shown.
    let record = start_open();
    let over_the_start = row(
        "deletion",
        record.at(-20, 8),
        28,
        "start_lost,coding_sequence_variant,5_prime_UTR_variant",
        BLANK,
    );
    assert_eq!(
        owed(&record, &over_the_start),
        ok("coding_sequence_variant")
    );
    assert_eq!(
        class(
            &record,
            &over_the_start,
            "coding_sequence_variant,5_prime_UTR_variant"
        ),
        Some(CodingEndClass::StartCodonOverlap)
    );
}

#[test]
fn structural_deletion_ending_in_the_last_whole_codon() {
    let record = shown_tail(1);
    let short_of_the_tail = row(
        "deletion",
        record.at(3, 15),
        13,
        "stop_lost,coding_sequence_variant,feature_truncation",
        BLANK,
    );
    assert_eq!(
        owed(&record, &short_of_the_tail),
        ok("coding_sequence_variant")
    );
    assert_eq!(
        class(
            &record,
            &short_of_the_tail,
            "coding_sequence_variant,feature_truncation"
        ),
        Some(CodingEndClass::StructuralDeletionStopLost)
    );
    assert_eq!(
        class(
            &record,
            &short_of_the_tail,
            "stop_lost,coding_sequence_variant,feature_truncation"
        ),
        None
    );
    let over_the_tail = row(
        "deletion",
        record.at(3, 16),
        14,
        "stop_lost,coding_sequence_variant,feature_truncation",
        BLANK,
    );
    assert_eq!(
        owed(&record, &over_the_tail),
        ok("stop_lost,coding_sequence_variant")
    );
}

// Both ends at once.

#[test]
fn deletion_through_both_coding_ends() {
    let record = both_open();
    let deletion = row(
        "-",
        record.at(3, 18),
        16,
        "start_lost,stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant",
        ["3-18", "-", "-", "-", "-"],
    );
    assert_eq!(owed(&record, &deletion), ok("coding_sequence_variant"));
    let outside = "5_prime_UTR_variant,3_prime_UTR_variant";
    for kept in ["start_lost", "stop_retained_variant"] {
        assert_eq!(
            class(
                &record,
                &deletion,
                &format!("{kept},coding_sequence_variant,{outside}")
            ),
            None
        );
    }
    assert_eq!(
        class(
            &record,
            &deletion,
            &format!("coding_sequence_variant,{outside}")
        ),
        Some(CodingEndClass::StopCodonOverlap)
    );
    let structural = row(
        "deletion",
        record.at(3, 18),
        16,
        "stop_lost,start_lost,coding_sequence_variant,feature_truncation",
        BLANK,
    );
    assert_eq!(owed(&record, &structural), ok("coding_sequence_variant"));
    assert_eq!(
        class(
            &record,
            &structural,
            "coding_sequence_variant,feature_truncation"
        ),
        Some(CodingEndClass::StructuralDeletionStopLost)
    );
}

#[test]
fn a_stop_term_of_a_stop_codon_that_stands_is_not_read_from_the_row_under_the_allele() {
    // A deletion from the 5' UTR into the 3' UTR over a stop codon that stands: vep-rs
    // reads the sequence the record gives there, which the cells do not show.
    let record = start_open();
    let deletion = row(
        "-",
        record.at(3, 20),
        18,
        "start_lost,stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant",
        ["3-20", "-", "-", "-", "-"],
    );
    assert_eq!(
        owed(&record, &deletion),
        Err(NotDerived::AtAStandingStopCodon)
    );
    assert_eq!(
        class(
            &record,
            &deletion,
            "stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant"
        ),
        None
    );
    // An insertion within its own length of the stop codon is not clear of it.
    let record = closed();
    let near = row(
        "CC",
        record.at(13, 14),
        0,
        "frameshift_variant",
        ["13-14", "7-8", "3", "G/GX", "gga/gCCga"],
    );
    assert_eq!(owed(&record, &near), Err(NotDerived::AtAStandingStopCodon));
    let clear = row(
        "C",
        record.at(10, 11),
        0,
        "frameshift_variant",
        ["10-11", "4-5", "2", "A/AX", "gct/gCct"],
    );
    assert_eq!(owed(&record, &clear), ok("frameshift_variant"));
}

#[test]
fn a_start_term_of_a_start_codon_that_stands_is_not_read_from_the_row_where_a_term_it_follows_moves(
) {
    // The coding sequence opens the transcript; a deletion from its second base to the
    // base before its last leaves a stop codon of the 3' UTR at the old stop position.
    // Ensembl VEP reads the stop kept and the deletion in frame; where no stop codon
    // stands the frameshift holds, and with it the start term.
    let record = Made {
        leader: 0,
        ..made(&[MET, ALA, GLY], 9)
    };
    let deletion = row(
        "-",
        record.at(2, 8),
        7,
        "inframe_deletion,stop_retained_variant",
        ["2-8", "2-8", "1-3", "MAG/X", "aTGGCTGGa/aa"],
    );
    assert_eq!(
        owed(&record, &deletion),
        Err(NotDerived::AtAStandingStartCodon)
    );
    assert_eq!(class(&record, &deletion, "frameshift_variant"), None);
}

#[test]
fn both_start_terms_on_a_structural_row_where_a_start_codon_stands() {
    let record = open_whole();
    let deletion = row(
        "deletion",
        record.at(3, 18),
        16,
        "stop_lost,start_lost,start_retained_variant,coding_sequence_variant",
        BLANK,
    );
    assert_eq!(
        owed(&record, &deletion),
        Err(NotDerived::BothStartTermsWhereAStartCodonStands)
    );
    assert_eq!(
        class(
            &record,
            &deletion,
            "start_lost,start_retained_variant,coding_sequence_variant"
        ),
        None
    );
}

#[test]
fn no_coding_term_on_a_span_over_the_whole_transcript() {
    let record = both_open();
    let longer = "ACGT".repeat(8);
    let replaced = Row {
        allele: longer,
        ..row("-", record.at(-3, 26), 30, "stop_retained_variant", BLANK)
    };
    assert_eq!(owed(&record, &replaced), ok(""));
    assert_eq!(class(&record, &replaced, "coding_sequence_variant"), None);
}

// What is not derived.

#[test]
fn allele_letter_that_names_no_base_is_not_derived() {
    let record = start_open();
    for letters in ["N", "R", "a", "*", "CNV", "<DEL>", ""] {
        let insertion = Row {
            allele: letters.to_string(),
            ..row(
                "-",
                record.at(7, 8),
                0,
                "start_retained_variant",
                ["7-8", "1-2", "1", "-", "-"],
            )
        };
        assert_eq!(
            owed(&record, &insertion),
            Err(NotDerived::AlleleNamesNoBases),
            "{letters}"
        );
        assert_eq!(class(&record, &insertion, ""), None, "{letters}");
    }
    assert!(STRUCTURAL_CLASS_WORDS.contains(&"deletion"));
    assert!(STRUCTURAL_CLASS_WORDS.contains(&"mobile_element_deletion"));
}

#[test]
fn a_row_that_lacks_what_the_derivation_reads_is_not_derived() {
    let record = open_whole();
    let deletion = row(
        "-",
        record.at(14, 14),
        1,
        "frameshift_variant,stop_lost",
        ["14", "8", "3", "G/X", "gGa/ga"],
    );
    let derived = |change: &dyn Fn(&mut Row)| {
        let mut changed = deletion.clone();
        change(&mut changed);
        owed(&record, &changed)
    };
    assert_eq!(
        derived(&|row| row.transcript_row = false),
        Err(NotDerived::NotATranscriptRow)
    );
    assert_eq!(derived(&|row| row.cells = None), Err(NotDerived::NoCells));
    assert_eq!(
        derived(&|row| row.reference_length = None),
        Err(NotDerived::NoReferenceLength)
    );
    for (cell, name) in [
        (0, "cDNA_position"),
        (1, "CDS_position"),
        (2, "Protein_position"),
        (3, "Amino_acids"),
        (4, "Codons"),
    ] {
        let mut ours = deletion.cells.unwrap();
        ours[cell] = "7";
        assert_eq!(
            derived(&|row| row.ours = Some(ours)),
            Err(NotDerived::CellsDiffer(name)),
            "{name}"
        );
    }
    // A length after the position is cut before the cells are held equal.
    let mut with_lengths = deletion.cells.unwrap();
    with_lengths[..3].copy_from_slice(&["14/21", "8/9", "3/3"]);
    assert_eq!(
        derived(&|row| row.ours = Some(with_lengths)),
        ok("frameshift_variant")
    );
    assert_eq!(
        derived(&|row| row.cells = Some(["14", "eight", "3", "G/X", "gGa/ga"])),
        Err(NotDerived::CellNotRead("CDS_position"))
    );
    // No record answers, or the record is not read.
    with(&deletion, None, |facts, read| {
        assert_eq!(
            expected_terms(facts, None, read),
            Err(NotDerived::NoVerdict)
        );
        assert_eq!(class_of_pair(facts, "frameshift_variant", None, read), None);
    });
    let unread = Made {
        codon_table: Some(7),
        ..open_whole()
    };
    assert_eq!(owed(&unread, &deletion), Err(NotDerived::NoVerdict));
    // Both start terms where a start codon stands.
    let record = closed();
    let both = row(
        "C",
        record.at(8, 8),
        1,
        "start_lost,start_retained_variant",
        ["8", "2", "1", "M/T", "aTg/aCg"],
    );
    assert_eq!(
        owed(&record, &both),
        Err(NotDerived::BothStartTermsWhereAStartCodonStands)
    );
}

// A row whose `Codons` and `Amino_acids` differ as the class
// `alternate_codon_past_the_coding_end` documents.

/// A coding sequence that ends two bases into a codon before a 3' UTR, without
/// `cds_end_NF`: no stop codon stands, and the last codon lacks one base.
fn open_tail() -> Made {
    made(&[MET, ALA, GLY, &STOP[..2]], 6)
}

const KEPT_STOP: &str = "inframe_deletion,stop_retained_variant";

/// The cells of a deletion of the last two bases of the last whole codon and of the two
/// bases behind it, on `open_tail`. Ensembl VEP's alternate codon reads on one base into
/// the 3' UTR (`READ_ON`); vep-rs's ends where the coding sequence ends (`CUT`).
const READ_ON: [&str; 5] = ["14-17", "8-11", "3-4", "GX/X", "gGATA/gc"];
const CUT: [&str; 5] = ["14-17", "8-11", "3-4", "GX/X", "gGATA/g"];

/// That deletion's row on `record`, with Ensembl VEP's cells `theirs` and vep-rs's `ours`.
fn deletion_into_the_tail(
    record: &Made,
    theirs: [&'static str; 5],
    ours: [&'static str; 5],
) -> Row {
    Row {
        ours: Some(ours),
        ..row("-", record.at(14, 17), 4, KEPT_STOP, theirs)
    }
}

#[test]
fn a_row_whose_alternate_codon_vep_rs_ends_at_the_coding_end_is_derived_on_the_cells_of_ensembl_vep(
) {
    let reverse = Made {
        strand: -1,
        ..open_tail()
    };
    for record in [open_tail(), reverse] {
        let deletion = deletion_into_the_tail(&record, READ_ON, CUT);
        assert_eq!(owed(&record, &deletion), ok("frameshift_variant"));
        assert_eq!(
            class(&record, &deletion, "frameshift_variant"),
            Some(CodingEndClass::StopCodonOverlap),
            "strand {}",
            record.strand
        );
        // No stop codon stands, so the rule of a three-base stop codon takes nothing
        // beside the cache.
        assert!(!may(
            LandedRule::StopRetainedThreeBaseStopCodonTest,
            Some(&record),
            &deletion
        ));
    }
    // A substitution of the first of the two bases: the base of the 3' UTR completes a
    // codon of asparagine in Ensembl VEP's row, and vep-rs's names the unknown residue
    // once. Both rows hold the same terms, which are the derived ones.
    let record = open_tail();
    let incomplete = "incomplete_terminal_codon_variant,coding_sequence_variant";
    let substitution = Row {
        ours: Some(["16", "10", "4", "X", "Ta/Aa"]),
        ..row(
            "A",
            record.at(16, 16),
            1,
            incomplete,
            ["16", "10", "4", "X/N", "Ta/Aac"],
        )
    };
    assert_eq!(owed(&record, &substitution), ok(incomplete));
    assert_eq!(class(&record, &substitution, incomplete), None);
}

#[test]
fn another_difference_of_the_codons_or_the_residues_is_not_derived() {
    let record = open_tail();
    let derived = |theirs: [&'static str; 5], ours: [&'static str; 5]| {
        owed(&record, &deletion_into_the_tail(&record, theirs, ours))
    };
    assert_eq!(derived(READ_ON, CUT), ok("frameshift_variant"));
    for (theirs, ours, what) in [
        ("gGATA/gc", "gGATA/a", "a changed base"),
        ("gGATA/gc", "gGATA/gcc", "a longer alternate codon"),
        ("gGATA/gc", "gGAtA/g", "another reference codon"),
        ("gGATA/gC", "gGATA/g", "a base of the allele"),
        ("gGATA/gcc", "gGATA/g", "two bases past the end"),
        ("gGATA/gccc", "gGATA/g", "three bases past the end"),
        ("gGATA/g", "gGATA/gc", "the two rows exchanged"),
        ("gGATA/gc", "-", "no codons"),
    ] {
        let mut reference = READ_ON;
        reference[4] = theirs;
        let mut written = CUT;
        written[4] = ours;
        assert_eq!(
            derived(reference, written),
            Err(NotDerived::CellsDiffer("Codons")),
            "{what}"
        );
    }
    // Residues that differ beside equal codons, and residues that are no translation
    // of vep-rs's alternate codon.
    for (codons, residues, what) in [
        ("gGATA/gc", "GX/G", "equal codons"),
        ("gGATA/g", "GX/G", "another residue"),
        ("gGATA/g", "GX", "one peptide"),
        ("gGATA/g", "AX/X", "another reference peptide"),
    ] {
        let mut written = CUT;
        written[3] = residues;
        written[4] = codons;
        assert_eq!(
            derived(READ_ON, written),
            Err(NotDerived::CellsDiffer("Amino_acids")),
            "{what}"
        );
    }
    // A position cell is held equal as on any row.
    let mut moved = CUT;
    moved[1] = "8-10";
    assert_eq!(
        derived(READ_ON, moved),
        Err(NotDerived::CellsDiffer("CDS_position"))
    );
}

#[test]
fn a_pair_whose_codons_differ_so_is_refused_on_other_terms_a_stop_codon_that_stands_and_without_a_record(
) {
    let record = open_tail();
    let deletion = deletion_into_the_tail(&record, READ_ON, CUT);
    assert_eq!(
        class(&record, &deletion, "frameshift_variant"),
        Some(CodingEndClass::StopCodonOverlap)
    );
    // Ensembl VEP's terms kept, a term dropped, a term added, a term outside the family.
    for ours in [
        KEPT_STOP,
        "inframe_deletion",
        "stop_retained_variant",
        "coding_sequence_variant",
        "frameshift_variant,stop_retained_variant",
        "frameshift_variant,splice_region_variant",
    ] {
        assert_eq!(class(&record, &deletion, ours), None, "{ours}");
    }
    // A stop codon that stands: vep-rs reads the edited sequence there.
    let stands = closed();
    let on_a_stop_codon = deletion_into_the_tail(&stands, READ_ON, CUT);
    assert_eq!(
        owed(&stands, &on_a_stop_codon),
        Err(NotDerived::AtAStandingStopCodon)
    );
    assert_eq!(class(&stands, &on_a_stop_codon, "frameshift_variant"), None);
    // `cds_end_NF`: the terms are derived, and Ensembl VEP reads the flag in its overlap.
    let flagged = Made {
        codes: vec!["cds_end_NF"],
        ..open_tail()
    };
    let on_a_flagged_end = deletion_into_the_tail(&flagged, READ_ON, CUT);
    assert_eq!(owed(&flagged, &on_a_flagged_end), ok("frameshift_variant"));
    assert_eq!(
        class(&flagged, &on_a_flagged_end, "frameshift_variant"),
        None
    );
    // No record answers.
    with(&deletion, None, |facts, read| {
        assert_eq!(
            expected_terms(facts, None, read),
            Err(NotDerived::NoVerdict)
        );
        assert_eq!(class_of_pair(facts, "frameshift_variant", None, read), None);
    });
}

#[test]
fn a_row_whose_record_keeps_too_little_is_not_derived() {
    // A flagged end beside a start codon that stands keeps no coordinates.
    let record = flagged();
    assert!(record.read().model.is_none());
    let past_the_end = row("-", (1010, 1030), 21, "3_prime_UTR_variant", BLANK);
    assert_eq!(owed(&record, &past_the_end), Err(NotDerived::ModelNotKept));
    let edited = Made {
        pairs: Some(vec![(1, 12, 1001, 1012, 1), (16, 24, 1013, 1021, 1)]),
        span: Some((1001, 1021)),
        ..flagged()
    };
    let insertion = row(
        "CTA",
        (1012, 1013),
        0,
        "inframe_insertion,stop_retained_variant",
        ["14-15", "8-9", "3", "G/G*", "gga/ggCTAa"],
    );
    assert_eq!(
        owed(&edited, &insertion),
        Err(NotDerived::InsertionOnAnEditedModel)
    );
    // A genetic code the table module lacks keeps no sequence.
    let record = Made {
        codon_table: Some(26),
        ..shown_by_an_id(2)
    };
    let substitution = row(
        "G",
        record.at(17, 17),
        1,
        INCOMPLETE,
        ["17", "11", "4", "X", "ta/tG"],
    );
    assert_eq!(
        owed(&record, &substitution),
        Err(NotDerived::SequenceNotKept)
    );
}

#[test]
fn a_row_on_a_completed_stop_codon_that_the_derivation_does_not_cover_is_not_derived() {
    let record = shown_tail(1);
    // The position cells differ between the rows.
    let lost = Row {
        ours: Some(["16", "10-11", "4", "*/Q", "Taa/Caa"]),
        ..row(
            "C",
            record.at(16, 16),
            1,
            INCOMPLETE,
            ["16", "10", "4", "X", "T/C"],
        )
    };
    assert_eq!(
        owed(&record, &lost),
        Err(NotDerived::CellsDiffer("CDS_position"))
    );
    // The codons hold protein position 1, or a residue the cached peptide states
    // otherwise than its codon.
    let first = row(
        "C",
        record.at(8, 8),
        1,
        "missense_variant",
        ["8", "2", "1", "M/T", "aTg/aCg"],
    );
    assert_eq!(owed(&record, &first), Err(NotDerived::SequenceEdit));
    let edited = Made {
        peptide: Some("MAU"),
        ..shown_tail(1)
    };
    let third = row(
        "C",
        edited.at(14, 14),
        1,
        "missense_variant",
        ["14", "8", "3", "U/A", "gGa/gCa"],
    );
    assert_eq!(owed(&edited, &third), Err(NotDerived::SequenceEdit));
    let second = row(
        "C",
        edited.at(11, 11),
        1,
        "missense_variant",
        ["11", "5", "2", "A/A", "gCt/gCt"],
    );
    // A residue of the same record that the cached peptide states as its codon does is
    // derived.
    assert_eq!(owed(&edited, &second), ok("synonymous_variant"));
    // A span over the tail that lacks its cDNA start.
    let from_before = row(
        "-",
        record.at(-5, 16),
        22,
        "stop_retained_variant",
        ["?-16", "?-10", "?-4", "-", "-"],
    );
    let apart = genome(&record, "GGG");
    assert_eq!(
        owed_beside(&record, &from_before, Some(&apart)),
        Err(NotDerived::SpanLacksACoordinate)
    );
}

// The classes.

#[test]
fn each_class_refuses_a_pair_on_a_transcript_whose_codon_stands_or_that_carries_the_flag() {
    let deletion = |record: &Made| {
        row(
            "-",
            record.at(14, 14),
            1,
            "frameshift_variant,stop_lost",
            ["14", "8", "3", "G/X", "gGa/ga"],
        )
    };
    // The stop codon overlap: a stop codon that stands, and the flag, alone and beside
    // another attribute code.
    for record in [
        closed(),
        flagged(),
        Made {
            codes: vec!["TSL", "cds_end_NF"],
            ..open_whole()
        },
    ] {
        assert_eq!(
            class(&record, &deletion(&record), "frameshift_variant"),
            None
        );
    }
    // The completed stop codon: whole codons, a tail the record does not show complete,
    // and the flag.
    let lost = |record: &Made| {
        row(
            "C",
            record.at(16, 16),
            1,
            INCOMPLETE,
            ["16", "10", "4", "X", "T/C"],
        )
    };
    for record in [
        closed(),
        Made {
            codon_table: None,
            ..shown_tail(1)
        },
        Made {
            codes: vec!["cds_end_NF"],
            ..shown_tail(1)
        },
    ] {
        assert_eq!(class(&record, &lost(&record), "stop_lost"), None);
    }
    // The stop codon behind the last residue: a stop codon that stands; and the two
    // peptides holding their first stop at one place.
    let behind = |record: &Made| {
        row(
            "CTA",
            record.at(14, 15),
            0,
            "inframe_insertion,stop_retained_variant",
            ["14-15", "8-9", "3", "G/G*", "gga/ggCTAa"],
        )
    };
    assert_eq!(
        class(
            &closed(),
            &behind(&closed()),
            "stop_gained,inframe_insertion"
        ),
        None
    );
    let inner = made(&[MET, "TAG", "CAA"], 6);
    let one_place = row(
        "-",
        inner.at(11, 13),
        3,
        "inframe_deletion,stop_retained_variant",
        ["11-13", "5-7", "2-3", "*Q/*", "tAGCaa/taa"],
    );
    assert_eq!(class(&inner, &one_place, "inframe_deletion"), None);
    // The structural deletion: a stop codon that stands, and a row without the term.
    let over_the_end =
        |record: &Made, theirs: &'static str| row("deletion", record.at(10, 40), 31, theirs, BLANK);
    assert_eq!(
        class(
            &closed(),
            &over_the_end(&closed(), "stop_lost,coding_sequence_variant"),
            "coding_sequence_variant"
        ),
        None
    );
    assert_eq!(
        class(
            &flagged(),
            &over_the_end(&flagged(), "coding_sequence_variant"),
            "coding_sequence_variant"
        ),
        None
    );
    // The start codon overlap: a start codon that stands, by its codon or by a
    // consensus id, and the flag.
    let substitution = |record: &Made| {
        row(
            "T",
            record.at(9, 9),
            1,
            "start_lost",
            ["9", "3", "1", "E/D", "gaA/gaT"],
        )
    };
    for record in [
        closed(),
        Made {
            consensus_id: Some("CCDS1.1"),
            ..start_open()
        },
        Made {
            peptide: Some("MAG"),
            ..start_open()
        },
        Made {
            codes: vec!["cds_start_NF"],
            ..start_open()
        },
    ] {
        assert_eq!(
            class(&record, &substitution(&record), "missense_variant"),
            None
        );
    }
}

#[test]
fn each_class_refuses_a_row_of_vep_rs_that_drops_adds_or_keeps_a_term() {
    let open = open_whole();
    let deletion = row(
        "-",
        open.at(14, 14),
        1,
        "frameshift_variant,stop_lost,splice_region_variant",
        ["14", "8", "3", "G/X", "gGa/ga"],
    );
    let taken = "frameshift_variant,splice_region_variant";
    assert_eq!(
        class(&open, &deletion, taken),
        Some(CodingEndClass::StopCodonOverlap)
    );
    for ours in [
        "splice_region_variant",
        "frameshift_variant",
        "stop_gained,frameshift_variant,splice_region_variant",
        "frameshift_variant,stop_lost,splice_region_variant",
        "frameshift_variant,stop_retained_variant,splice_region_variant",
        "frameshift_variant,splice_region_variant,intron_variant",
    ] {
        assert_eq!(class(&open, &deletion, ours), None, "{ours}");
    }
    let shown = shown_tail(1);
    let lost = row(
        "C",
        shown.at(16, 16),
        1,
        INCOMPLETE,
        ["16", "10", "4", "X", "T/C"],
    );
    for ours in [
        "",
        "stop_retained_variant",
        "stop_lost,incomplete_terminal_codon_variant",
        INCOMPLETE,
    ] {
        assert_eq!(class(&shown, &lost, ours), None, "{ours}");
    }
    let behind = row(
        "CTA",
        open.at(14, 15),
        0,
        "inframe_insertion,stop_retained_variant",
        ["14-15", "8-9", "3", "G/G*", "gga/ggCTAa"],
    );
    for ours in [
        "stop_gained",
        "inframe_insertion",
        "stop_gained,inframe_insertion,stop_retained_variant",
        "inframe_insertion,stop_retained_variant",
    ] {
        assert_eq!(class(&open, &behind, ours), None, "{ours}");
    }
    let end = flagged();
    let structural = row(
        "deletion",
        end.at(10, 40),
        31,
        "stop_lost,coding_sequence_variant,feature_truncation",
        BLANK,
    );
    for ours in [
        "coding_sequence_variant",
        "feature_truncation",
        "stop_retained_variant,coding_sequence_variant,feature_truncation",
        "stop_lost,coding_sequence_variant,feature_truncation",
    ] {
        assert_eq!(class(&end, &structural, ours), None, "{ours}");
    }
    let start = start_open();
    let substitution = row(
        "T",
        start.at(9, 9),
        1,
        "start_lost,start_retained_variant",
        ["9", "3", "1", "E/D", "gaA/gaT"],
    );
    for ours in [
        "",
        "missense_variant,stop_gained",
        "start_retained_variant",
        "start_lost,start_retained_variant",
        "start_retained_variant,missense_variant",
    ] {
        assert_eq!(class(&start, &substitution, ours), None, "{ours}");
    }
}

#[test]
fn a_class_is_named_as_the_class_file_names_it() {
    let names: Vec<&str> = CodingEndClass::IN_ORDER
        .iter()
        .map(|class| class.name())
        .collect();
    assert_eq!(
        names,
        [
            "structural_deletion_stop_lost_without_stop_codon",
            "stop_codon_completed_by_polyadenylation",
            "stop_codon_overlap_without_stop_codon",
            "stop_codon_behind_last_residue_without_stop_codon",
            "start_codon_overlap_without_start_codon",
        ]
    );
    for class in CodingEndClass::IN_ORDER {
        assert_eq!(CodingEndClass::of_name(class.name()), Some(class));
    }
    assert_eq!(CodingEndClass::of_name("start_cooccurrence_swap"), None);
    for rule in LandedRule::ALL {
        assert_eq!(LandedRule::of_name(rule.name()), Some(rule));
    }
    assert_eq!(
        coding_family("stop_lost,intron_variant, stop_gained"),
        BTreeSet::from([CodingTerm::StopGained, CodingTerm::StopLost])
    );
    assert_eq!(
        outside_the_coding_family("stop_lost,intron_variant"),
        BTreeSet::from(["intron_variant"])
    );
    assert_eq!(CodingTerm::ALL.len(), 13);
    assert_eq!(position_cell(Some("8"), Some("8")), "8");
    assert_eq!(position_cell(Some("8"), None), "8-?");
    assert_eq!(position_cell(None, Some("9")), "?-9");
    assert_eq!(position_cell(Some("8"), Some("9")), "8-9");
    assert_eq!(position_cell(None, None), "-");
}

// The landed rules beside a cache.

/// Whether a comparison given the cache may let `rule` take a pair on the row, on the
/// record or with no record at all.
fn may(rule: LandedRule, record: Option<&Made>, row: &Row) -> bool {
    let ends = record.map(Made::read);
    with(row, None, |facts, read| {
        landed_rule_may_take(rule, facts, ends.as_ref(), read)
    })
}

const STOP_RULES: [LandedRule; 6] = [
    LandedRule::StopCodonRemovedWithoutStopCodonTest,
    LandedRule::StopRetainedWithoutStopCodonTest,
    LandedRule::StopRetainedThreeBaseStopCodonTest,
    LandedRule::StopCodonOverlapInsertionRange,
    LandedRule::StopCodonOverlapGenomicSpan,
    LandedRule::StopCodonUnshiftedSpan,
];

#[test]
fn a_landed_rule_takes_its_own_pairs_wherever_the_stop_codon_stands() {
    use LandedRule::*;
    // An insertion that keeps a stop codon beside a splice site: at a stop codon that
    // an intron splits after its first base, and at one that no intron splits.
    let split = Made {
        pairs: Some(vec![(1, 16, 1001, 1016, 1), (17, 24, 1101, 1108, 1)]),
        span: Some((1001, 1108)),
        ..closed()
    };
    let pair = "coding_sequence_variant,inframe_insertion,splice_region_variant";
    let at_the_junction = row(
        "CTAA",
        (1015, 1016),
        0,
        pair,
        ["15-16", "9-10", "3-4", "-/X", "-/N"],
    );
    let whole = closed();
    let unsplit = row(
        "CTAA",
        whole.at(16, 17),
        0,
        pair,
        ["16-17", "10-11", "4", "*/X", "taa/tNaa"],
    );
    // A deletion from the 5' UTR into the 3' UTR, over the whole coding sequence and
    // every base of the stop codon, that Ensembl VEP reads as a stop codon kept.
    let over_the_coding_sequence = row(
        "-",
        whole.at(3, 20),
        18,
        "3_prime_UTR_variant,5_prime_UTR_variant,start_lost,stop_retained_variant",
        ["3-20", "-", "-", "-", "-"],
    );
    // An insertion on the reverse strand before the stop codon by fewer bases than it
    // inserts, which Ensembl VEP's overlap reaches: eight bases four before the codon,
    // and thirty-one bases twenty-seven before it.
    let reverse = Made {
        strand: -1,
        ..made(
            &[MET, ALA, GLY, ALA, GLY, ALA, GLY, ALA, GLY, ALA, GLY, STOP],
            6,
        )
    };
    let four_before = row(
        "CCCCCCCC",
        reverse.at(35, 36),
        0,
        "protein_altering_variant,stop_retained_variant",
        ["35-36", "29-30", "10", "A/X", "gct/gcN"],
    );
    let longer = "ACGT".repeat(8);
    let twenty_seven_before = Row {
        allele: longer[..31].to_string(),
        ..row(
            "-",
            reverse.at(12, 13),
            0,
            "inframe_insertion,stop_retained_variant",
            ["12-13", "6-7", "2-3", "-/X", "-/N"],
        )
    };
    let shapes = [
        (
            "an insertion at a stop codon that an intron splits",
            may(StopCodonOverlapGenomicSpan, Some(&split), &at_the_junction),
        ),
        (
            "an insertion at a stop codon that no intron splits",
            may(StopCodonOverlapGenomicSpan, Some(&whole), &unsplit),
        ),
        (
            "a deletion over the whole coding sequence",
            may(
                StopRetainedWithoutStopCodonTest,
                Some(&whole),
                &over_the_coding_sequence,
            ),
        ),
        (
            "eight bases inserted four before the stop codon on the reverse strand",
            may(
                StopRetainedThreeBaseStopCodonTest,
                Some(&reverse),
                &four_before,
            ),
        ),
        (
            "thirty-one bases inserted twenty-seven before it",
            may(
                StopRetainedThreeBaseStopCodonTest,
                Some(&reverse),
                &twenty_seven_before,
            ),
        ),
    ];
    let refused: Vec<&str> = shapes
        .iter()
        .filter(|(_, taken)| !taken)
        .map(|(shape, _)| *shape)
        .collect();
    assert!(refused.is_empty(), "refused: {refused:?}");
    // Wherever a whole stop codon stands every rule of a stop codon takes, whatever
    // the place of the allele: the place is the rule's own reading of the two rows.
    for row in [
        &at_the_junction,
        &unsplit,
        &over_the_coding_sequence,
        &four_before,
    ] {
        for rule in STOP_RULES {
            for record in [&whole, &split, &reverse] {
                assert!(may(rule, Some(record), row), "{}", rule.name());
            }
        }
    }
    // Where polyadenylation completes the stop codon each takes too, but the rule of
    // one triplet read two ways, since the engines read two different triplets there,
    // and the rule of a stop codon kept for want of a test, which holds only where the
    // genome carries the codon whole.
    for shown in [shown_tail(1), shown_tail(2), shown_by_an_id(2)] {
        for rule in STOP_RULES {
            assert_eq!(
                may(rule, Some(&shown), &over_the_coding_sequence),
                rule != StopRetainedThreeBaseStopCodonTest
                    && rule != StopRetainedWithoutStopCodonTest,
                "{}",
                rule.name()
            );
        }
    }
}

#[test]
fn a_landed_rule_takes_nothing_beside_a_cache_where_its_codon_does_not_stand() {
    use LandedRule::*;
    let whole = closed();
    let rows = [
        row(
            "-",
            whole.at(14, 30),
            17,
            "stop_retained_variant",
            ["14-?", "8-?", "3-?", "-", "-"],
        ),
        row(
            "C",
            whole.at(16, 17),
            0,
            "stop_retained_variant",
            ["16-17", "10-11", "4", "*/X", "taa/tCaa"],
        ),
        row(
            "-",
            whole.at(3, 20),
            18,
            "start_lost,stop_retained_variant",
            ["3-20", "-", "-", "-", "-"],
        ),
    ];
    // No stop codon stands (whole codons whose last names a residue, the flag, a tail
    // no record shows completed), the record is not read, or no record holds the row.
    let unread = Made {
        codon_table: Some(7),
        ..closed()
    };
    for record in [
        Some(open_whole()),
        Some(flagged()),
        Some(tail_not_shown()),
        Some(unread.clone()),
        None,
    ] {
        for rule in STOP_RULES {
            for row in &rows {
                assert!(!may(rule, record.as_ref(), row), "{}", rule.name());
            }
        }
    }
    // The start pair: only where a start codon stands.
    let first = row(
        "C",
        whole.at(8, 8),
        1,
        "start_lost,start_retained_variant",
        ["8", "2", "1", "M/T", "aTg/aCg"],
    );
    assert!(may(StartCooccurrenceSwap, Some(&whole), &first));
    assert!(may(StartCooccurrenceSwap, Some(&open_whole()), &first));
    let flagged_start = Made {
        codes: vec!["cds_start_NF"],
        ..closed()
    };
    for record in [Some(start_open()), Some(flagged_start), Some(unread), None] {
        assert!(!may(StartCooccurrenceSwap, record.as_ref(), &first));
    }
    // The stop named on an incomplete coding sequence reads as it does without the
    // coding ends, but for a completed stop codon, on which vep-rs owes a notation.
    assert!(may(
        HgvspStopNamedOnIncompleteCodingSequence,
        Some(&tail_not_shown()),
        &first
    ));
    assert!(may(
        HgvspStopNamedOnIncompleteCodingSequence,
        Some(&whole),
        &first
    ));
    assert!(may(HgvspStopNamedOnIncompleteCodingSequence, None, &first));
    for shown in [shown_tail(1), shown_by_an_id(2)] {
        assert!(!may(
            HgvspStopNamedOnIncompleteCodingSequence,
            Some(&shown),
            &first
        ));
    }
}

#[test]
fn the_rule_of_a_stop_codon_kept_for_want_of_a_test_takes_only_where_the_genome_carries_the_codon_whole(
) {
    use LandedRule::*;
    // A deletion from inside the coding sequence past the transcript's last base, on a
    // transcript whose stop codon polyadenylation completes and on one whose stop
    // codon's last two bases an edit of the model adds, with no genome base.
    let completed = shown_tail(1);
    let edited = Made {
        pairs: Some(vec![(1, 16, 1001, 1016, 1)]),
        span: Some((1001, 1016)),
        ..made(&[MET, ALA, GLY, STOP], 0)
    };
    assert_eq!(edited.read().last_codon, [Some(1016), None, None]);
    let past_the_end = row(
        "-",
        (1012, 1030),
        19,
        "stop_retained_variant,3_prime_UTR_variant",
        ["12-?", "6-?", "2-?", "-", "-"],
    );
    // The row with no stop term on vep-rs's side has the rule's shape and is none of
    // its pairs there; the row with the lost stop is the removed stop codon's on both.
    let taken: Vec<&str> = [("completed", &completed), ("edited", &edited)]
        .into_iter()
        .filter(|(_, record)| {
            may(
                StopRetainedWithoutStopCodonTest,
                Some(record),
                &past_the_end,
            )
        })
        .map(|(name, _)| name)
        .collect();
    assert!(taken.is_empty(), "taken: {taken:?}");
    for record in [&completed, &edited] {
        assert!(may(
            StopCodonRemovedWithoutStopCodonTest,
            Some(record),
            &past_the_end
        ));
    }
    // A deletion over the whole coding sequence of a transcript whose genome carries
    // the stop codon whole stays the rule's, on either strand and across an intron.
    let whole = closed();
    let over_the_coding_sequence = row(
        "-",
        whole.at(3, 20),
        18,
        "3_prime_UTR_variant,5_prime_UTR_variant,start_lost,stop_retained_variant",
        ["3-20", "-", "-", "-", "-"],
    );
    let split = Made {
        pairs: Some(vec![(1, 16, 1001, 1016, 1), (17, 24, 1101, 1108, 1)]),
        span: Some((1001, 1108)),
        ..closed()
    };
    let reverse = Made {
        strand: -1,
        ..closed()
    };
    for record in [&whole, &split, &reverse] {
        assert!(may(
            StopRetainedWithoutStopCodonTest,
            Some(record),
            &over_the_coding_sequence
        ));
    }
}

// The protein notation of a frameshift at a start end where no start codon is shown.

#[test]
fn a_frameshift_at_a_first_codon_that_is_no_start_codon_is_written_to_the_new_stop() {
    let written = |record: &Made, row: &Row| {
        let ends = record.read();
        with(row, None, |facts, read| {
            protein_notation_where_no_start_codon_is_shown(facts, Some(&ends), read, ("", ""))
        })
    };
    let lost = "frameshift_variant,start_lost";
    // A first codon of glutamate; the frame one base on meets a stop codon in the
    // second codon's place.
    let early_stop = made(&[GLU, "CTA", "AGC", GLY, STOP], 6);
    let deletion = row(
        "-",
        early_stop.at(7, 7),
        1,
        lost,
        ["7", "1", "1", "E/X", "Gaa/aa"],
    );
    assert_eq!(written(&early_stop, &deletion), ok("p.Glu1AsnfsTer2"));
    // One base inserted in the first codon, with a stop codon in the third codon's
    // place of the new frame.
    let third = made(&[GLU, "CTT", "AAG", GLY, STOP], 6);
    let insertion = row(
        "C",
        third.at(7, 8),
        0,
        lost,
        ["7-8", "1-2", "1", "E/AX", "gaa/gCaa"],
    );
    assert_eq!(written(&third, &insertion), ok("p.Glu1AlafsTer3"));
    // The new frame runs on into the 3' UTR: a stop codon there, and none.
    let stop_in_the_utr = Made {
        three_prime_utr: Some("CTAGCC"),
        ..start_open()
    };
    let no_stop = Made {
        three_prime_utr: Some("CCCCCC"),
        ..start_open()
    };
    let deletion = |record: &Made| {
        row(
            "-",
            record.at(8, 8),
            1,
            lost,
            ["8", "2", "1", "E/X", "gAa/ga"],
        )
    };
    assert_eq!(
        written(&stop_in_the_utr, &deletion(&stop_in_the_utr)),
        ok("p.Ala2LeufsTer4")
    );
    assert_eq!(
        written(&no_stop, &deletion(&no_stop)),
        ok("p.Ala2LeufsTer?")
    );
    let insertion = |record: &Made| {
        row(
            "C",
            record.at(7, 8),
            0,
            lost,
            ["7-8", "1-2", "1", "E/AX", "gaa/gCaa"],
        )
    };
    let later_stop = Made {
        three_prime_utr: Some("CCTAAC"),
        ..start_open()
    };
    assert_eq!(
        written(&later_stop, &insertion(&later_stop)),
        ok("p.Glu1AlafsTer6")
    );
    assert_eq!(
        written(&no_stop, &insertion(&no_stop)),
        ok("p.Glu1AlafsTer?")
    );
    // On the reverse strand the cell holds the inserted base's complement.
    let reverse = Made {
        strand: -1,
        ..later_stop.clone()
    };
    let complement = Row {
        allele: "G".to_string(),
        ..insertion(&reverse)
    };
    assert_eq!(written(&reverse, &complement), ok("p.Glu1AlafsTer6"));
    // Not derived: a start codon that is shown, by its codon or by a consensus id; a
    // start that carries the flag; a row on which Ensembl VEP writes no start term; a
    // record whose 3' UTR the reader did not hand over; a record given on the reverse
    // strand.
    let shown = made(&[MET, ALA, GLY, STOP], 6);
    let by_an_id = Made {
        consensus_id: Some("CCDS1.1"),
        ..no_stop.clone()
    };
    let flagged_start = Made {
        codes: vec!["cds_start_NF"],
        ..no_stop.clone()
    };
    let unread_utr = Made {
        utr_read: false,
        ..no_stop.clone()
    };
    for record in [&shown, &by_an_id, &flagged_start, &unread_utr] {
        assert_eq!(
            written(record, &deletion(record)),
            Err(NotDerived::CellNotDerived)
        );
    }
    let elsewhere = Row {
        theirs: "frameshift_variant",
        ..deletion(&no_stop)
    };
    assert_eq!(
        written(&no_stop, &elsewhere),
        Err(NotDerived::CellNotDerived)
    );
    let given_reversed = Row {
        reverse_record: true,
        ..deletion(&no_stop)
    };
    assert_eq!(
        written(&no_stop, &given_reversed),
        Err(NotDerived::CellNotDerived)
    );
}

// A deletion across a coding start that shows no start codon, on which Ensembl VEP reads a
// lost start and the notation is made on the deletion moved 3'.

#[test]
fn a_deletion_across_a_coding_start_that_shows_no_start_codon_is_read_on_its_moved_span() {
    let written = |record: &Made, row: &Row, offsets: (&str, &str)| {
        let ends = record.read();
        with(row, None, |facts, read| {
            protein_notation_where_no_start_codon_is_shown(facts, Some(&ends), read, offsets)
        })
    };
    // Six bases of 5' UTR before a first codon of proline. The last base of the 5' UTR
    // and the first coding base are deleted; moved one base 3' the deletion takes coding
    // positions 1 and 2, and the frame two bases on meets a stop codon in the 3' UTR.
    let record = Made {
        three_prime_utr: Some(STOP_IN_THE_UTR),
        ..made(&FLAGGED_CODONS, 6)
    };
    let across = |record: &Made, last: i64, cds: &'static str| {
        row(
            "-",
            record.at(6, last),
            (last - 5) as u64,
            "start_lost,5_prime_UTR_variant",
            [if last == 7 { "6-7" } else { "6-8" }, cds, cds, "-", "-"],
        )
    };
    let deletion = across(&record, 7, "?-1");
    assert_eq!(
        owed(&record, &deletion),
        ok("coding_sequence_variant"),
        "the terms derived for the row hold no frameshift"
    );
    assert_eq!(
        written(&record, &deletion, ("1", "1")),
        ok("p.Pro1AsnfsTer6")
    );
    // Not derived: no offset, an offset on one row alone, two offsets that differ, a
    // deletion of three bases, which moves no frame, a span the offset moves past the
    // coding sequence, and a transcript that carries the flag.
    for offsets in [("", ""), ("-", "-"), ("1", ""), ("", "1"), ("1", "2")] {
        assert_eq!(
            written(&record, &deletion, offsets),
            Err(NotDerived::CellNotDerived),
            "{offsets:?}"
        );
    }
    assert_eq!(
        written(&record, &across(&record, 8, "?-2"), ("1", "1")),
        Err(NotDerived::CellNotDerived)
    );
    assert_eq!(
        written(&record, &deletion, ("15", "15")),
        Err(NotDerived::CellNotDerived)
    );
    let flagged_start = Made {
        codes: vec!["cds_start_NF"],
        ..record.clone()
    };
    assert_eq!(
        written(&flagged_start, &deletion, ("1", "1")),
        Err(NotDerived::CellNotDerived)
    );
}

// The frameshift named at residue 1 where no start codon stands and neither engine reads a
// lost start.

/// The codons of a made-up coding sequence that opens with proline and is closed by a stop
/// codon. One base on, its frame reads glutamine and a stop codon; two bases on, it reads
/// asparagine, lysine, arginine and isoleucine and runs on into the 3' UTR.
const FLAGGED_CODONS: [&str; 5] = ["CCA", "ATA", "AGC", GLY, STOP];

/// A 3' UTR of six bases that holds a stop codon in the frame two bases on.
const STOP_IN_THE_UTR: &str = concat!("GC", "TAA", "C");

/// `FLAGGED_CODONS` with `cds_start_NF`, no 5' UTR and the 3' UTR `STOP_IN_THE_UTR`.
fn flagged_start() -> Made {
    Made {
        codes: vec!["cds_start_NF"],
        leader: 0,
        three_prime_utr: Some(STOP_IN_THE_UTR),
        ..made(&FLAGGED_CODONS, 6)
    }
}

const NO_OFFSET: (&str, &str) = ("", "");

/// The two notations derived for the row on the record, vep-rs's and Ensembl VEP's, each
/// after `p.`; `ours` is vep-rs's `Consequence` cell and `offsets` the two `HGVS_OFFSET`
/// cells.
fn counted_beside(
    record: Option<&Made>,
    row: &Row,
    ours: &str,
    offsets: (&str, &str),
) -> Result<(String, String), NotDerived> {
    let ends = record.map(Made::read);
    with(row, None, |facts, read| {
        frameshift_counted_at_the_first_residue(facts, ours, ends.as_ref(), read, offsets)
            .map(|counted| (counted.vep_rs_notation(), counted.ensembl_vep_notation()))
    })
}

/// `counted_beside` for a row both engines give the terms of `row.theirs`.
fn counted(record: &Made, row: &Row) -> Result<(String, String), NotDerived> {
    counted_beside(Some(record), row, row.theirs, NO_OFFSET)
}

/// The pair of notations of the frameshift `frameshift`, `<Ref>1<Alt>fsTer`, counted to
/// `count`.
fn to_the_stop(frameshift: &str, count: usize) -> Result<(String, String), NotDerived> {
    Ok((format!("{frameshift}{count}"), format!("{frameshift}?")))
}

const NOT_DERIVED: Result<(String, String), NotDerived> = Err(NotDerived::CellNotDerived);

/// The deletion of the first coding base of `record`, cDNA position `cdna`.
fn first_base_deleted(record: &Made, cdna: i64) -> Row {
    let position: &'static str = if cdna == 1 { "1" } else { "7" };
    row(
        "-",
        record.at(cdna, cdna),
        1,
        FRAMESHIFT,
        [position, "1", "1", "P/X", "Cca/ca"],
    )
}

/// The deletion of the first two coding bases of a record without a 5' UTR.
fn first_two_bases_deleted(record: &Made) -> Row {
    row(
        "-",
        record.at(1, 2),
        2,
        FRAMESHIFT,
        ["1-2", "1-2", "1", "P/X", "CCa/a"],
    )
}

#[test]
fn a_frameshift_at_residue_1_of_a_flagged_start_is_counted_to_the_stop_codon_of_its_frame() {
    let record = flagged_start();
    // The new frame meets a stop codon inside the coding sequence, and in the 3' UTR.
    assert_eq!(
        counted(&record, &first_base_deleted(&record, 1)),
        to_the_stop("Pro1GlnfsTer", 2)
    );
    assert_eq!(
        counted(&record, &first_two_bases_deleted(&record)),
        to_the_stop("Pro1AsnfsTer", 6)
    );
    // The same on the reverse strand, where a deletion's cell holds no base.
    let reverse = Made {
        strand: -1,
        ..record.clone()
    };
    assert_eq!(
        counted(&reverse, &first_base_deleted(&reverse, 1)),
        to_the_stop("Pro1GlnfsTer", 2)
    );
    assert_eq!(
        counted(&reverse, &first_two_bases_deleted(&reverse)),
        to_the_stop("Pro1AsnfsTer", 6)
    );
    // With a 5' UTR the cDNA position counts its bases and the coding position does not.
    let with_a_leader = Made {
        leader: 6,
        ..record.clone()
    };
    assert_eq!(
        counted(&with_a_leader, &first_base_deleted(&with_a_leader, 7)),
        to_the_stop("Pro1GlnfsTer", 2)
    );
    // A transcript that ends at its coding end has no 3' UTR to read: the stop codon
    // inside the coding sequence is counted to, with one flag and with both.
    let ends_at_its_coding_end = Made {
        three_prime_utr: None,
        trailer: 0,
        ..record.clone()
    };
    let flagged_at_both_ends = Made {
        codes: vec!["cds_start_NF", "cds_end_NF"],
        coding: FLAGGED_CODONS[..4].concat(),
        ..ends_at_its_coding_end.clone()
    };
    for record in [&ends_at_its_coding_end, &flagged_at_both_ends] {
        assert_eq!(
            counted(record, &first_base_deleted(record, 1)),
            to_the_stop("Pro1GlnfsTer", 2)
        );
    }
    // The flag decides: a first codon that is the usual start codon, and a consensus id.
    let opens_with_a_start_codon = Made {
        coding: [MET, "ATA", "AGC", GLY, STOP].concat(),
        ..record.clone()
    };
    let two_deleted = Row {
        cells: Some(["1-2", "1-2", "1", "M/X", "ATg/g"]),
        ..first_two_bases_deleted(&opens_with_a_start_codon)
    };
    assert_eq!(
        counted(&opens_with_a_start_codon, &two_deleted),
        to_the_stop("Met1AspfsTer", 6)
    );
    let with_a_consensus_id = Made {
        consensus_id: Some("CCDS1.1"),
        ..record.clone()
    };
    assert_eq!(
        counted(
            &with_a_consensus_id,
            &first_base_deleted(&with_a_consensus_id, 1)
        ),
        to_the_stop("Pro1GlnfsTer", 2)
    );
    // A cell that carries the length of `--total_length` is read up to it.
    let with_lengths = Row {
        cells: Some(["1/21", "1/15", "1/4", "P/X", "Cca/ca"]),
        ..first_base_deleted(&record, 1)
    };
    assert_eq!(
        counted(&record, &with_lengths),
        to_the_stop("Pro1GlnfsTer", 2)
    );
}

#[test]
fn the_reference_length_tells_an_insertion_from_a_replacement_of_two_bases() {
    let record = flagged_start();
    // One base between coding positions 1 and 2, and the same base in the place of both:
    // the two rows print one allele and one `CDS_position`.
    let inserted = row(
        "A",
        record.at(1, 2),
        0,
        FRAMESHIFT,
        ["1-2", "1-2", "1", "P/HX", "cca/cAca"],
    );
    let replaced = row(
        "A",
        record.at(1, 2),
        2,
        FRAMESHIFT,
        ["1-2", "1-2", "1", "P/X", "CCa/Aa"],
    );
    assert_eq!(counted(&record, &inserted), to_the_stop("Pro1HisfsTer", 7));
    assert_eq!(counted(&record, &replaced), to_the_stop("Pro1LysfsTer", 2));
    assert_eq!(
        counted(&record, &first_two_bases_deleted(&record)),
        to_the_stop("Pro1AsnfsTer", 6)
    );
    // On the reverse strand the cell holds the complement of the base the transcript
    // reads.
    let reverse = Made {
        strand: -1,
        ..record.clone()
    };
    let complement = |row: &Row| Row {
        allele: "T".to_string(),
        location: reverse.at(1, 2),
        ..row.clone()
    };
    assert_eq!(
        counted(&reverse, &complement(&inserted)),
        to_the_stop("Pro1HisfsTer", 7)
    );
    assert_eq!(
        counted(&reverse, &complement(&replaced)),
        to_the_stop("Pro1LysfsTer", 2)
    );
    assert_eq!(counted(&reverse, &inserted), to_the_stop("Pro1LeufsTer", 7));
}

/// A deletion of `length` bases that ends on the coding position `printed_end` of a
/// record without a 5' UTR and begins before the record's first base: the cells are open
/// at their start and the terms hold no frameshift.
fn deletion_from_before_the_first_base(record: &Made, length: i64, printed_end: i64) -> Row {
    let cell: &'static str = match printed_end {
        1 => "?-1",
        2 => "?-2",
        _ => "?-3",
    };
    row(
        "-",
        record.at(printed_end + 1 - length, printed_end),
        length as u64,
        "coding_sequence_variant,5_prime_UTR_variant",
        [cell, cell, cell, "-", "-"],
    )
}

#[test]
fn a_deletion_from_before_the_coding_start_is_read_on_the_span_its_offset_moves_it_to() {
    let record = flagged_start();
    let offset =
        |row: &Row, offset: &str| counted_beside(Some(&record), row, row.theirs, (offset, offset));
    // A whole printed span is read as it is printed, whatever the offset.
    let first = first_base_deleted(&record, 1);
    assert_eq!(offset(&first, "1"), to_the_stop("Pro1GlnfsTer", 2));
    assert_eq!(
        counted_beside(Some(&record), &first, FRAMESHIFT, ("1", "")),
        to_the_stop("Pro1GlnfsTer", 2)
    );
    // Two bases that end on coding position 1, moved one base: coding positions 1 and 2.
    let two = deletion_from_before_the_first_base(&record, 2, 1);
    assert_eq!(offset(&two, "1"), to_the_stop("Pro1AsnfsTer", 6));
    // Five bases that end on coding position 2, moved three: coding positions 1 to 5.
    let five = deletion_from_before_the_first_base(&record, 5, 2);
    assert_eq!(offset(&five, "3"), to_the_stop("Pro1LysfsTer", 5));
    // Moved one base further the two bases leave another residue at residue 1.
    assert_eq!(offset(&two, "2"), to_the_stop("Pro1HisfsTer", 6));
    // The reverse strand, and a transcript that carries no flag and shows no start codon.
    let reverse = Made {
        strand: -1,
        ..record.clone()
    };
    let without_the_flag = Made {
        codes: Vec::new(),
        ..record.clone()
    };
    for record in [&reverse, &without_the_flag] {
        let two = deletion_from_before_the_first_base(record, 2, 1);
        assert_eq!(
            counted_beside(Some(record), &two, two.theirs, ("1", "1")),
            to_the_stop("Pro1AsnfsTer", 6)
        );
    }
    // Not derived: no offset, an offset on one row alone, two that differ, an offset that
    // is no whole number above 0, a moved span that begins before the coding sequence or
    // ends past it, a deletion that moves no frame, and an allele with bases over an open
    // start.
    for offsets in [
        ("", ""),
        ("-", "-"),
        ("1", ""),
        ("", "1"),
        ("1", "2"),
        ("0", "0"),
        ("-1", "-1"),
        ("+1", "+1"),
        ("1.0", "1.0"),
    ] {
        assert_eq!(
            counted_beside(Some(&record), &two, two.theirs, offsets),
            NOT_DERIVED,
            "{offsets:?}"
        );
    }
    assert_eq!(offset(&five, "1"), NOT_DERIVED);
    assert_eq!(offset(&two, "20"), NOT_DERIVED);
    let three = deletion_from_before_the_first_base(&record, 3, 2);
    assert_eq!(offset(&three, "1"), NOT_DERIVED);
    let inserted = Row {
        allele: "A".to_string(),
        reference_length: Some(0),
        ..two.clone()
    };
    assert_eq!(offset(&inserted, "1"), NOT_DERIVED);
    let replaced = Row {
        allele: "A".to_string(),
        ..two.clone()
    };
    assert_eq!(offset(&replaced, "1"), NOT_DERIVED);
}

#[test]
fn no_count_is_derived_where_the_frame_meets_no_stop_codon_or_the_first_codon_is_padded() {
    // A 3' UTR without a stop codon in the frame two bases on, a 3' UTR the transcript
    // does not have: both engines owe the unknown count there.
    let no_stop_in_the_utr = Made {
        three_prime_utr: None,
        ..flagged_start()
    };
    let no_utr = Made {
        three_prime_utr: None,
        trailer: 0,
        ..flagged_start()
    };
    for record in [&no_stop_in_the_utr, &no_utr] {
        assert_eq!(
            counted(record, &first_two_bases_deleted(record)),
            NOT_DERIVED
        );
        assert_eq!(
            counted(record, &first_base_deleted(record, 1)),
            to_the_stop("Pro1GlnfsTer", 2)
        );
    }
    // A start phase of 1 and of 2: the first codon holds padding letters on both sides
    // and reads as the unknown residue on both, so the walk passes residue 1.
    for (phase, first_codon, position) in [(1, "NCA", "2"), (2, "NNA", "3")] {
        let padded = Made {
            start_phase: phase,
            coding: [first_codon, "ATA", "AGC", GLY, STOP].concat(),
            ..flagged_start()
        };
        let first_genome_base = row(
            "-",
            padded.at(1, 1),
            1,
            FRAMESHIFT,
            ["1", position, "1", "X/X", "-"],
        );
        assert_eq!(counted(&padded, &first_genome_base), NOT_DERIVED, "{phase}");
    }
}

#[test]
fn a_row_the_count_at_the_first_residue_is_not_derived_for() {
    let record = flagged_start();
    let deletion = first_base_deleted(&record, 1);
    assert_eq!(counted(&record, &deletion), to_the_stop("Pro1GlnfsTer", 2));
    let refused = |record: &Made, row: &Row| counted(record, row).err();
    // The record: none, one that is not read, a start codon that stands by its codon and
    // one a consensus id shows.
    assert_eq!(
        counted_beside(None, &deletion, FRAMESHIFT, NO_OFFSET),
        Err(NotDerived::NoVerdict)
    );
    let unread = Made {
        codon_table: Some(77),
        ..record.clone()
    };
    assert_eq!(refused(&unread, &deletion), Some(NotDerived::NoVerdict));
    let stands = Made {
        codes: Vec::new(),
        coding: [MET, "ATA", "AGC", GLY, STOP].concat(),
        ..record.clone()
    };
    let shown_by_an_id = Made {
        codes: Vec::new(),
        consensus_id: Some("CCDS1.1"),
        ..record.clone()
    };
    for record in [&stands, &shown_by_an_id] {
        assert_eq!(
            refused(record, &first_two_bases_deleted(record)),
            Some(NotDerived::CellNotDerived)
        );
    }
    // The sequence is not kept: a 3' UTR the reader did not hand over, a stop codon that
    // polyadenylation completes, a genetic code the module holds no table for.
    let unread_utr = Made {
        utr_read: false,
        ..record.clone()
    };
    let completed = Made {
        codes: vec!["cds_start_NF"],
        leader: 0,
        codon_table: Some(2),
        ..made(&["CCA", "ATA", "AGC", GLY, &STOP[..1]], 0)
    };
    let another_module = Made {
        codon_table: Some(26),
        ..record.clone()
    };
    for record in [&unread_utr, &completed, &another_module] {
        assert_eq!(
            refused(record, &first_base_deleted(record, 1)),
            Some(NotDerived::SequenceNotKept)
        );
    }
    // The row: no Transcript's, a record given on the reverse strand, no cells, no
    // reference length, an allele that names no bases, a structural allele.
    let altered = |alter: &dyn Fn(&mut Row)| {
        let mut row = deletion.clone();
        alter(&mut row);
        refused(&record, &row)
    };
    assert_eq!(
        altered(&|row| row.transcript_row = false),
        Some(NotDerived::NotATranscriptRow)
    );
    assert_eq!(
        altered(&|row| row.reverse_record = true),
        Some(NotDerived::RecordOnTheReverseStrand)
    );
    assert_eq!(altered(&|row| row.cells = None), Some(NotDerived::NoCells));
    assert_eq!(
        altered(&|row| row.reference_length = None),
        Some(NotDerived::NoReferenceLength)
    );
    assert_eq!(
        altered(&|row| row.allele = "N".to_string()),
        Some(NotDerived::AlleleNamesNoBases)
    );
    assert_eq!(
        altered(&|row| row.allele = "deletion".to_string()),
        Some(NotDerived::CellNotDerived)
    );
    // A lost start on either row.
    assert_eq!(
        altered(&|row| row.theirs = "frameshift_variant,start_lost"),
        Some(NotDerived::CellNotDerived)
    );
    assert_eq!(
        counted_beside(
            Some(&record),
            &deletion,
            "frameshift_variant,start_lost",
            NO_OFFSET
        ),
        NOT_DERIVED
    );
    let lost_on_their_side = Row {
        theirs: "frameshift_variant,start_lost",
        ..deletion.clone()
    };
    assert_eq!(
        counted_beside(Some(&record), &lost_on_their_side, FRAMESHIFT, NO_OFFSET),
        NOT_DERIVED
    );
    // A span that is not read: a blank cell, a cell open at its end or at both, a
    // `CDS_position` that differs between the rows.
    for cds in ["-", "1-?", "?", "?-?"] {
        assert_eq!(
            altered(&|row| row.cells = Some(["1", cds, "1", "P/X", "Cca/ca"])),
            Some(match cds {
                "?-?" => NotDerived::CellNotRead("CDS_position"),
                _ => NotDerived::CellNotDerived,
            }),
            "{cds}"
        );
    }
    assert_eq!(
        altered(&|row| row.ours = Some(["1", "2", "1", "P/X", "Cca/ca"])),
        Some(NotDerived::CellsDiffer("CDS_position"))
    );
    // An allele that moves no frame, a frameshift named at residue 2, one whose walk
    // passes a first residue that stays, and a residue 1 that becomes a stop.
    let three_deleted = row(
        "-",
        record.at(1, 3),
        3,
        "inframe_deletion",
        ["1-3", "1-3", "1", "P/-", "CCA/-"],
    );
    let fourth_base_deleted = row(
        "-",
        record.at(4, 4),
        1,
        FRAMESHIFT,
        ["4", "4", "2", "I/X", "Ata/ta"],
    );
    let third_base_deleted = row(
        "-",
        record.at(3, 3),
        1,
        FRAMESHIFT,
        ["3", "3", "1", "P/X", "ccA/cc"],
    );
    for row in [&three_deleted, &fourth_base_deleted, &third_base_deleted] {
        assert_eq!(refused(&record, row), Some(NotDerived::CellNotDerived));
    }
    let becomes_a_stop = Made {
        coding: ["CTA", "AGC", GLY, STOP].concat(),
        ..record.clone()
    };
    assert_eq!(
        refused(&becomes_a_stop, &first_base_deleted(&becomes_a_stop, 1)),
        Some(NotDerived::CellNotDerived)
    );
    // A first residue outside the twenty: a selenocysteine, the unknown residue, a stop.
    for peptide in ["UISG", "XISG", "*ISG"] {
        let outside = Made {
            peptide: Some(peptide),
            ..record.clone()
        };
        assert_eq!(
            refused(&outside, &first_base_deleted(&outside, 1)),
            Some(NotDerived::CellNotDerived),
            "{peptide}"
        );
    }
}

#[test]
fn the_two_derived_notations_are_held_to_the_two_cells_behind_one_protein() {
    let record = flagged_start();
    let deletion = first_base_deleted(&record, 1);
    let ends = record.read();
    let derived = with(&deletion, None, |facts, read| {
        frameshift_counted_at_the_first_residue(facts, FRAMESHIFT, Some(&ends), read, NO_OFFSET)
    })
    .expect("the row is derived");
    let pairs = [
        // The pair, and the pair in prediction parentheses on both sides.
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP1.1:p.Pro1GlnfsTer2", true),
        (
            "ENSP1.1:p.(Pro1GlnfsTer?)",
            "ENSP1.1:p.(Pro1GlnfsTer2)",
            true,
        ),
        // Parentheses on one side alone, another protein, the sides exchanged.
        (
            "ENSP1.1:p.(Pro1GlnfsTer?)",
            "ENSP1.1:p.Pro1GlnfsTer2",
            false,
        ),
        (
            "ENSP1.1:p.Pro1GlnfsTer?",
            "ENSP1.1:p.(Pro1GlnfsTer2)",
            false,
        ),
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP2.1:p.Pro1GlnfsTer2", false),
        ("ENSP1.1:p.Pro1GlnfsTer2", "ENSP1.1:p.Pro1GlnfsTer?", false),
        // Another count, a count on both sides, the unknown count on both.
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP1.1:p.Pro1GlnfsTer3", false),
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP1.1:p.Pro1GlnfsTer1", false),
        ("ENSP1.1:p.Pro1GlnfsTer2", "ENSP1.1:p.Pro1GlnfsTer2", false),
        ("ENSP1.1:p.Pro1GlnfsTer3", "ENSP1.1:p.Pro1GlnfsTer2", false),
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP1.1:p.Pro1GlnfsTer?", false),
        // Another residue on either side or on both, another number.
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP1.1:p.Pro1HisfsTer2", false),
        ("ENSP1.1:p.Pro1HisfsTer?", "ENSP1.1:p.Pro1GlnfsTer2", false),
        ("ENSP1.1:p.Pro1HisfsTer?", "ENSP1.1:p.Pro1HisfsTer2", false),
        ("ENSP1.1:p.Ala1GlnfsTer?", "ENSP1.1:p.Ala1GlnfsTer2", false),
        ("ENSP1.1:p.Pro2GlnfsTer?", "ENSP1.1:p.Pro2GlnfsTer2", false),
        // Ensembl VEP's lost start, a cell without a notation, a blank cell.
        ("ENSP1.1:p.Pro1?", "ENSP1.1:p.Pro1GlnfsTer2", false),
        ("-", "ENSP1.1:p.Pro1GlnfsTer2", false),
        ("ENSP1.1:p.Pro1GlnfsTer?", "", false),
    ];
    for (theirs, ours, taken) in pairs {
        assert_eq!(
            derived.takes(theirs, ours),
            taken,
            "{theirs} against {ours}"
        );
    }
    // The shape of the pair, read without a record: one frameshift named at residue 1
    // behind one protein, the unknown count against a number.
    let shaped = [
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP1.1:p.Pro1GlnfsTer2", true),
        (
            "ENSP1.1:p.(Pro1GlnfsTer?)",
            "ENSP1.1:p.(Pro1GlnfsTer2)",
            true,
        ),
        ("ENSP1.1:p.Ala1HisfsTer?", "ENSP1.1:p.Ala1HisfsTer310", true),
        ("ENSP1.1:p.Pro1GlnfsTer2", "ENSP1.1:p.Pro1GlnfsTer?", false),
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP1.1:p.Pro1GlnfsTer?", false),
        ("ENSP1.1:p.Pro1GlnfsTer3", "ENSP1.1:p.Pro1GlnfsTer2", false),
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP1.1:p.Pro1HisfsTer2", false),
        ("ENSP1.1:p.Pro2GlnfsTer?", "ENSP1.1:p.Pro2GlnfsTer2", false),
        (
            "ENSP1.1:p.Pro11GlnfsTer?",
            "ENSP1.1:p.Pro11GlnfsTer2",
            false,
        ),
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP2.1:p.Pro1GlnfsTer2", false),
        (
            "ENSP1.1:p.(Pro1GlnfsTer?)",
            "ENSP1.1:p.Pro1GlnfsTer2",
            false,
        ),
        (
            "ENSP1.1:p.Pro1GlnextTer?",
            "ENSP1.1:p.Pro1GlnextTer2",
            false,
        ),
        ("ENSP1.1:p.Pro1?", "ENSP1.1:p.Pro1GlnfsTer2", false),
        ("ENSP1.1:p.Pro1GlnfsTer?", "ENSP1.1:p.Pro1GlnfsTer2x", false),
        ("-", "-", false),
    ];
    for (theirs, ours, shape) in shaped {
        assert_eq!(
            count_alone_differs_at_the_first_residue(theirs, ours),
            shape,
            "{theirs} against {ours}"
        );
    }
}

#[test]
fn no_count_is_derived_from_a_frame_the_sequence_does_not_give() {
    let record = flagged_start();
    let two_deleted = first_two_bases_deleted(&record);
    assert_eq!(
        counted(&record, &two_deleted),
        to_the_stop("Pro1AsnfsTer", 6)
    );
    // A codon of the new frame before its first stop codon holds a letter that is no
    // base, in the coding sequence or in the 3' UTR: that codon may itself be a stop
    // codon. A stop codon met before it is counted to.
    let unread_in_the_coding_sequence = Made {
        coding: ["CCA", "ATA", "AGC", "GNA", STOP].concat(),
        ..record.clone()
    };
    let unread_in_the_utr = Made {
        three_prime_utr: Some(concat!("NC", "TAA", "C")),
        ..record.clone()
    };
    for record in [&unread_in_the_coding_sequence, &unread_in_the_utr] {
        assert_eq!(counted(record, &two_deleted), NOT_DERIVED);
        assert_eq!(
            counted(record, &first_base_deleted(record, 1)),
            to_the_stop("Pro1GlnfsTer", 2)
        );
    }
    // Fewer than three coding bases are left: both engines drop them and read the 3' UTR
    // from its own first base, a frame the transcript does not have.
    let short = Made {
        codes: vec!["cds_start_NF", "cds_end_NF"],
        coding: ["CCA", "ATA"].concat(),
        three_prime_utr: Some(concat!("AAG", "TAA")),
        ..record.clone()
    };
    let four_deleted = row(
        "-",
        short.at(1, 4),
        4,
        FRAMESHIFT,
        ["1-4", "1-4", "1-2", "PI/X", "CCAAta/ta"],
    );
    assert_eq!(counted(&short, &four_deleted), NOT_DERIVED);
    assert_eq!(
        counted(&short, &first_base_deleted(&short, 1)),
        to_the_stop("Pro1GlnfsTer", 2)
    );
    // The record's reference bases are more than the span's positions, as where a
    // deletion takes an intron with two coding bases: the cells do not show the edit.
    let with_an_intron = Row {
        reference_length: Some(102),
        ..two_deleted.clone()
    };
    assert_eq!(counted(&record, &with_an_intron), NOT_DERIVED);
    // An insertion's cell is two neighbouring positions: one position, or two that lie
    // further apart, is a span of bases the record does not replace.
    for cds in ["1", "1-3"] {
        let inserted = row(
            "A",
            record.at(1, 2),
            0,
            FRAMESHIFT,
            ["1-2", cds, "1", "P/HX", "cca/cAca"],
        );
        assert_eq!(counted(&record, &inserted), NOT_DERIVED, "{cds}");
    }
    // Three bases taken away move no frame, on the span as printed and on a moved span.
    let three_deleted = row(
        "-",
        record.at(1, 3),
        3,
        FRAMESHIFT,
        ["1-3", "1-3", "1", "P/-", "CCA/-"],
    );
    assert_eq!(counted(&record, &three_deleted), NOT_DERIVED);
    let three_moved = deletion_from_before_the_first_base(&record, 3, 1);
    assert_eq!(
        counted_beside(Some(&record), &three_moved, three_moved.theirs, ("2", "2")),
        NOT_DERIVED
    );
}

// The protein notation of a frameshift that the coding end cuts short, where no stop codon
// stands.

const LOST_STOP: &str = "frameshift_variant,stop_lost";
const FRAMESHIFT: &str = "frameshift_variant";

/// Whole codons whose last names glutamate on the transcript's last base: no stop codon
/// stands and no 3' UTR follows.
fn open_end() -> Made {
    made(&[MET, ALA, GLY, GLU], 0)
}

/// A coding sequence that opens one base into a codon, with `cds_start_NF` and no 5' UTR,
/// and ends two bases into a codon on the transcript's last base; adenines after those
/// two spell the codon of serine.
fn padded_tail() -> Made {
    Made {
        codes: vec!["cds_start_NF"],
        start_phase: 1,
        leader: 0,
        ..made(&["NGA", ALA, GLY, "TC"], 0)
    }
}

/// The deletion of the first base of the last codon of a record shaped as `open_end`.
fn first_base_of_the_last_codon(record: &Made) -> Row {
    row(
        "-",
        record.at(16, 16),
        1,
        LOST_STOP,
        ["16", "10", "4", "E/X", "Gaa/aa"],
    )
}

/// vep-rs's notation and Ensembl VEP's, each after `p.`.
type Notations = Result<(Option<String>, String), NotDerived>;

/// What `frameshift_the_coding_end_cuts_short` gives for the row with vep-rs's terms
/// `ours`.
fn cut_short(record: Option<&Made>, row: &Row, ours: &str) -> Notations {
    let ends = record.map(Made::read);
    with(row, None, |facts, read| {
        frameshift_the_coding_end_cuts_short(facts, ours, ends.as_ref(), read)
            .map(|cut| (cut.vep_rs_notation(), cut.ensembl_vep_notation()))
    })
}

/// The pair of a row whose frameshift routine names `deletion`: vep-rs's `frameshift` and
/// Ensembl VEP's deletion with an extension of unknown length.
fn notations(frameshift: &str, deletion: &str) -> Notations {
    Ok((Some(frameshift.to_string()), format!("{deletion}extTer?")))
}

#[test]
fn a_frameshift_the_coding_end_cuts_short_is_named_at_the_first_of_its_reference_residues() {
    // One base of the last of whole codons, on either strand.
    for strand in [1, -1] {
        let record = Made {
            strand,
            ..open_end()
        };
        let row = first_base_of_the_last_codon(&record);
        assert_eq!(
            cut_short(Some(&record), &row, FRAMESHIFT),
            notations("Glu4fs", "Glu4del"),
            "strand {strand}"
        );
    }
    // The last two bases of one codon and the first two of the last: two residues, of
    // which the frameshift names the first.
    let record = open_end();
    let two = row(
        "-",
        record.at(14, 17),
        4,
        LOST_STOP,
        ["14-17", "8-11", "3-4", "GE/X", "gGAGAa/ga"],
    );
    assert_eq!(
        cut_short(Some(&record), &two, FRAMESHIFT),
        notations("Gly3fs", "Gly3_Glu4del")
    );
    // The last two bases of the last whole codon and a tail of two, behind a padded first
    // codon, on either strand: the unknown residue of the tail stands on both peptides and
    // is clipped.
    for strand in [1, -1] {
        let record = Made {
            strand,
            ..padded_tail()
        };
        let tail = row(
            "-",
            record.at(7, 10),
            4,
            LOST_STOP,
            ["7-10", "8-11", "3-4", "GX/X", "gGATC/g"],
        );
        assert_eq!(
            cut_short(Some(&record), &tail, FRAMESHIFT),
            notations("Gly3fs", "Gly3del"),
            "strand {strand}"
        );
    }
}

#[test]
fn a_frameshift_notation_is_derived_only_where_the_coding_end_cuts_it_short() {
    let open = open_end();
    let row = first_base_of_the_last_codon(&open);
    let derived =
        |record: Option<&Made>, row: &Row, ours: &str| cut_short(record, row, ours).is_ok();
    let closed = made(&[MET, ALA, GLY, GLU, STOP], 0);
    let flagged = Made {
        codes: vec!["cds_end_NF"],
        ..open_end()
    };
    let before_a_utr = made(&[MET, ALA, GLY, GLU], 6);
    let utr_not_read = Made {
        utr_read: false,
        ..before_a_utr.clone()
    };
    let edited = Made {
        peptide: Some("MAGK"),
        ..open_end()
    };
    let longer = made(&[MET, ALA, GLY, GLU, GLU], 0);
    let cases = [
        (
            "the row on the open end",
            derived(Some(&open), &row, FRAMESHIFT),
        ),
        ("no record", derived(None, &row, FRAMESHIFT)),
        (
            "a stop codon stands",
            derived(Some(&closed), &row, FRAMESHIFT),
        ),
        ("cds_end_NF", derived(Some(&flagged), &row, FRAMESHIFT)),
        (
            "Ensembl VEP keeps the stop",
            derived(
                Some(&open),
                &Row {
                    theirs: "inframe_deletion,stop_retained_variant",
                    ..row.clone()
                },
                FRAMESHIFT,
            ),
        ),
        (
            "Ensembl VEP loses no stop",
            derived(
                Some(&open),
                &Row {
                    theirs: FRAMESHIFT,
                    ..row.clone()
                },
                FRAMESHIFT,
            ),
        ),
        (
            "vep-rs loses the stop",
            derived(Some(&open), &row, LOST_STOP),
        ),
        (
            "vep-rs writes another term",
            derived(Some(&open), &row, "frameshift_variant,stop_gained"),
        ),
        (
            "a 3' UTR the frame reads on into",
            derived(Some(&before_a_utr), &row, FRAMESHIFT),
        ),
        (
            "a 3' UTR that is not read",
            derived(Some(&utr_not_read), &row, FRAMESHIFT),
        ),
        (
            "a codon the sequence still holds whole",
            derived(Some(&longer), &row, FRAMESHIFT),
        ),
        (
            "a cached peptide that names another residue",
            derived(Some(&edited), &row, FRAMESHIFT),
        ),
        (
            "a residue cell that is not the record's",
            derived(
                Some(&open),
                &Row {
                    cells: Some(["16", "10", "4", "D/X", "Gaa/aa"]),
                    ..row.clone()
                },
                FRAMESHIFT,
            ),
        ),
        (
            "a record given on the reverse strand",
            derived(
                Some(&open),
                &Row {
                    reverse_record: true,
                    ..row.clone()
                },
                FRAMESHIFT,
            ),
        ),
        (
            "no cells",
            derived(
                Some(&open),
                &Row {
                    cells: None,
                    ..row.clone()
                },
                FRAMESHIFT,
            ),
        ),
        (
            "no reference length",
            derived(
                Some(&open),
                &Row {
                    reference_length: None,
                    ..row.clone()
                },
                FRAMESHIFT,
            ),
        ),
    ];
    let owed: Vec<(&str, bool)> = cases
        .iter()
        .enumerate()
        .map(|(place, (case, _))| (*case, place == 0))
        .collect();
    assert_eq!(cases.to_vec(), owed);
}

// Numbers at the edge of what 64 bits hold: a cell, a Location and a mapper pair are read
// as any other number, and a sum past that edge as the number it is.

/// The largest whole number 64 bits hold, as a cell writes it.
const LARGEST: &str = "9223372036854775807";

/// A reference sequence that holds no base.
const NO_BASE: &dyn Fn(i64) -> Option<u8> = &|_| None;

/// The reading of a record of the coding sequence `coding` under the vertebrate
/// mitochondrial code, whose span is the bases 5 to 20 and whose mapper states the two
/// coding bounds `bounds` and `pairs`.
fn read_on(
    coding: &str,
    strand: i64,
    start_phase: i64,
    bounds: (i64, i64),
    pairs: &[Pair],
) -> CodingEnds {
    coding_ends(&RecordFacts {
        strand: Some(strand),
        codon_table: Some(2),
        coding: Some(coding),
        mapper: Some(MapperFacts {
            start_phase,
            coding_start: bounds.0,
            coding_end: bounds.1,
            pairs,
        }),
        span: Some((5, 20)),
        ..RecordFacts::default()
    })
}

fn owed_on(ends: &CodingEnds, row: &Row, reference: Reference<'_>) -> Result<String, NotDerived> {
    with(row, reference, |facts, read| {
        expected_terms(facts, Some(ends), read).map(|terms| names(&terms))
    })
}

#[test]
fn a_translation_end_too_large_to_count_in_bases_is_read_to_the_coding_end() {
    let record = open_end();
    let row = Row {
        cells: Some(["16", "10", "4-4000000000000000000", "E/X", "Gaa/aa"]),
        ..first_base_of_the_last_codon(&record)
    };
    assert_eq!(
        cut_short(Some(&record), &row, FRAMESHIFT),
        notations("Glu4fs", "Glu4del")
    );
}

#[test]
fn a_coding_span_longer_than_64_bits_count_is_a_span_of_that_length() {
    let record = open_whole();
    // One position more than the largest whole number, and a start that far past its end.
    let deletion = row(
        "-",
        record.at(14, 14),
        1,
        "frameshift_variant,stop_lost",
        ["14", "0-9223372036854775807", "3", "G/X", "gGa/ga"],
    );
    assert_eq!(owed(&record, &deletion), ok("frameshift_variant"));
    let replacement = row(
        "CC",
        record.at(14, 14),
        1,
        "frameshift_variant,stop_lost",
        ["14", "9223372036854775807-0", "3", "G/X", "gGa/gCCa"],
    );
    assert_eq!(owed(&record, &replacement), ok("frameshift_variant"));
}

#[test]
fn an_insertion_whose_location_opens_at_the_lowest_whole_number_is_clear_of_the_stop_codon() {
    let record = start_open();
    let insertion = row("A", (i64::MIN, 5), 0, "upstream_gene_variant", BLANK);
    assert_eq!(owed(&record, &insertion), ok(""));
}

#[test]
fn a_mapper_pair_at_the_edge_of_64_bits_places_a_row_as_any_other_pair() {
    // One pair over every genomic base: an insertion stands between two neighbours of it.
    let everywhere = Made {
        pairs: Some(vec![(1, 21, i64::MIN, i64::MAX, 1)]),
        ..open_whole()
    };
    let insertion = row("A", (1014, 1015), 0, "downstream_gene_variant", BLANK);
    assert_eq!(owed(&everywhere, &insertion), ok(""));
    // A pair whose coding bases lie past the largest genomic base holds no base of a row.
    let past_the_last_base = Made {
        pairs: Some(vec![(1, 21, i64::MAX - 3, i64::MAX, 1)]),
        ..open_whole()
    };
    let deletion = row("-", (1014, 1014), 1, "downstream_gene_variant", BLANK);
    assert_eq!(owed(&past_the_last_base, &deletion), ok(""));
}

#[test]
fn a_completed_stop_codon_is_read_on_numbers_at_the_edge_of_64_bits() {
    let coding = [GLU, ALA, GLY, "TA"].concat();
    let in_the_tail = |cdna_position: &'static str, cds_position: &'static str| {
        row(
            "C",
            (FIRST + 15, FIRST + 15),
            1,
            INCOMPLETE,
            [cdna_position, cds_position, "4", "X", "Taa/Caa"],
        )
    };
    // A coding end at the largest cDNA position.
    let last_position = read_on(
        &coding,
        1,
        0,
        (i64::MAX - 10, i64::MAX),
        &[(1, i64::MAX, FIRST, FIRST + 16, 1)],
    );
    assert_eq!(
        owed_on(&last_position, &in_the_tail("16", "10"), None),
        ok("stop_lost")
    );
    // A cDNA position that the padding of the start exon's phase carries past that edge.
    let padded = read_on(&coding, 1, 2, (1, 9), &[(1, 9, FIRST, FIRST + 8, 1)]);
    assert_eq!(
        owed_on(&padded, &in_the_tail(LARGEST, "10"), None),
        ok("stop_lost")
    );
    // Two tail bases as far apart as two genomic bases lie.
    let apart = read_on(
        &coding,
        1,
        0,
        (1, 11),
        &[
            (1, 10, i64::MIN, i64::MIN + 9, 1),
            (11, 11, i64::MAX, i64::MAX, 1),
        ],
    );
    let elsewhere = row("-", (5, 6), 2, "downstream_gene_variant", BLANK);
    assert_eq!(owed_on(&apart, &elsewhere, None), ok(""));
    // A translation start and a coding span past what 64 bits hold name no codon.
    let record = shown_tail(2);
    let window = |cds_position: &'static str, protein_position: &'static str| {
        row(
            "C",
            record.at(16, 16),
            1,
            INCOMPLETE,
            ["16", cds_position, protein_position, "X", "Taa/Caa"],
        )
    };
    assert_eq!(owed(&record, &window("10", "4000000000000000000")), ok(""));
    assert_eq!(
        owed(&record, &window("0-9223372036854775807", "4")),
        ok("stop_lost")
    );
}

#[test]
fn a_deletion_past_the_transcript_end_is_read_on_numbers_at_the_edge_of_64_bits() {
    let past_the_end = |location: (i64, i64)| {
        row(
            "-",
            location,
            3,
            INCOMPLETE,
            ["16-?", "10-?", "4-?", "-", "-"],
        )
    };
    // To the largest genomic base: no base follows the deletion.
    let record = shown_tail(1);
    let to_the_last_base = past_the_end((record.genomic(16), i64::MAX));
    assert_eq!(
        owed_beside(&record, &to_the_last_base, Some(NO_BASE)),
        ok("stop_lost")
    );
    // On the reverse strand from the lowest genomic base, where the tail base lies.
    let reverse = read_on(
        &[GLU, ALA, GLY, "T"].concat(),
        -1,
        0,
        (7, 16),
        &[(1, 16, i64::MIN, i64::MIN + 15, -1)],
    );
    assert_eq!(
        owed_on(&reverse, &past_the_end((i64::MIN, 4)), Some(NO_BASE)),
        ok("stop_lost")
    );
}
