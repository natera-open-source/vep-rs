// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The coding sequence a transcript's `translation` row of the core dumps
//! gives it, where the row and the model's CDS records disagree.
//!
//! Ensembl VEP's cache is dumped from the core database, so a cached
//! translation starts and ends at its row's `seq_start` and `seq_end`. A
//! release's GFF3 can place the same coding sequence otherwise: a `CDS` line
//! over a whole exon whose first bases the row, and the file's own
//! `five_prime_UTR` line, leave untranslated. Ensembl VEP's `--gff` and
//! `--gtf` sources read the CDS record's own ends (`BaseGXF.pm:738-756`), so
//! only a build that reads the core dumps holds a transcript to its row.

use std::borrow::Cow;

use crate::core_dumps::ExonPositions;
use crate::gff3::{CdsRecord, ExonRecord};

/// What a build makes of a transcript's translation row beside its CDS
/// records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RowUse {
    /// No row is read, the transcript has no coding sequence, or the row and
    /// the CDS records place the coding sequence alike.
    Agrees,
    /// The row places the coding sequence otherwise, inside the transcript's
    /// one exon: the transcript is built as the row places it.
    Placed,
    /// The row places the coding sequence otherwise and the exon each of its
    /// positions counts in is not known: the transcript is built as its CDS
    /// records place it.
    Unplaced,
}

/// The CDS records a transcript is built from, and what its translation row
/// came to.
///
/// `written` is the translation's start and end as `cds_regions`, the model's
/// records, place them: the first coding base's position within its exon and
/// the last coding base's within its own. Where `row` gives the same two
/// positions, or there is no row or no coding sequence, the records are the
/// model's.
///
/// Where the row differs and `ordered`, the transcript's exons, is one exon
/// that holds both of the row's positions, the coding sequence is that exon's
/// bases from `row.start` to `row.end`: one record, carrying the protein
/// identifier, version and line of the model's first CDS record in transcript
/// order, the record a translation's identifiers are read from. Its phase is
/// that record's where the row leaves the start in place and 0 where the row
/// moves it, a translation opening on a codon's first base.
///
/// A row that differs in a transcript of several exons places nothing: the
/// dumps name the row's two exons by database identifier, which no line of
/// the model carries. Neither does a row whose positions lie outside the one
/// exon, or whose two exon identifiers differ.
pub(crate) fn coding_records<'a>(
    ordered: &[&ExonRecord],
    cds_regions: &'a [CdsRecord],
    strand: i8,
    written: Option<(u64, u64)>,
    row: Option<ExonPositions>,
) -> (Cow<'a, [CdsRecord]>, RowUse) {
    let model = Cow::Borrowed(cds_regions);
    let (Some(written), Some(row)) = (written, row) else {
        return (model, RowUse::Agrees);
    };
    if written == (row.start, row.end) {
        return (model, RowUse::Agrees);
    }
    let [exon] = ordered else {
        return (model, RowUse::Unplaced);
    };
    let length = exon.end - exon.start + 1;
    if !row.in_one_exon || row.start < 1 || row.start > row.end || row.end > length {
        return (model, RowUse::Unplaced);
    }
    let (start, end) = if strand < 0 {
        (exon.end + 1 - row.end, exon.end + 1 - row.start)
    } else {
        (exon.start + row.start - 1, exon.start + row.end - 1)
    };
    let first = if strand < 0 {
        cds_regions.iter().max_by_key(|c| c.start)
    } else {
        cds_regions.iter().min_by_key(|c| c.start)
    };
    let Some(first) = first else {
        return (model, RowUse::Agrees);
    };
    let phase = if row.start == written.0 {
        first.phase
    } else {
        0
    };
    let record = CdsRecord {
        start,
        end,
        phase,
        ..first.clone()
    };
    (Cow::Owned(vec![record]), RowUse::Placed)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn exon(start: u64, end: u64, strand: i8) -> ExonRecord {
        ExonRecord {
            exon_id: None,
            start,
            end,
            strand,
            rank: None,
            ensembl_phase: None,
            ensembl_end_phase: None,
        }
    }

    fn cds(start: u64, end: u64, strand: i8, phase: i8, line: usize) -> CdsRecord {
        CdsRecord {
            start,
            end,
            strand,
            phase,
            protein_id: Some(format!("P{line}")),
            protein_version: Some(line as u32),
            line,
        }
    }

    fn row(start: u64, end: u64) -> Option<ExonPositions> {
        Some(ExonPositions {
            start,
            end,
            in_one_exon: true,
        })
    }

    /// The spans, phases and lines of `records`.
    fn spans(records: &[CdsRecord]) -> Vec<(u64, u64, i8, usize)> {
        records
            .iter()
            .map(|c| (c.start, c.end, c.phase, c.line))
            .collect()
    }

    /// A forward exon at 101 to 160 under a CDS record over all of it: a row
    /// that starts at base 19 gives one record from 119 to 160 at phase 0, with
    /// the model record's identifier, version and line. On the reverse strand
    /// base 19 of the exon is 142, so the record runs from 101 to 142.
    #[test]
    fn a_row_that_starts_inside_the_one_exon_gives_the_coding_sequence_its_span() {
        let forward = exon(101, 160, 1);
        let model = [cds(101, 160, 1, 0, 7)];
        let (records, row_use) = coding_records(&[&forward], &model, 1, Some((1, 60)), row(19, 60));
        assert_eq!(row_use, RowUse::Placed);
        assert_eq!(spans(&records), [(119, 160, 0, 7)]);
        assert_eq!(records[0].protein_id.as_deref(), Some("P7"));
        assert_eq!(records[0].protein_version, Some(7));

        let reverse = exon(101, 160, -1);
        let model = [cds(101, 160, -1, 0, 7)];
        let (records, row_use) =
            coding_records(&[&reverse], &model, -1, Some((1, 60)), row(19, 60));
        assert_eq!(row_use, RowUse::Placed);
        assert_eq!(spans(&records), [(101, 142, 0, 7)]);
    }

    /// A row that ends the coding sequence earlier than the records and
    /// leaves its start in place keeps the phase of the model's first record
    /// in transcript order, and its identifiers: the record of lowest start on
    /// the forward strand, of highest start on the reverse.
    #[test]
    fn a_row_that_moves_the_end_alone_keeps_the_first_records_phase() {
        let forward = exon(101, 160, 1);
        let model = [cds(131, 160, 1, 0, 9), cds(101, 130, 1, 2, 8)];
        let (records, row_use) = coding_records(&[&forward], &model, 1, Some((1, 60)), row(1, 45));
        assert_eq!(row_use, RowUse::Placed);
        assert_eq!(spans(&records), [(101, 145, 2, 8)]);

        let reverse = exon(101, 160, -1);
        let model = [cds(101, 130, -1, 0, 8), cds(131, 160, -1, 1, 9)];
        let (records, row_use) = coding_records(&[&reverse], &model, -1, Some((1, 60)), row(1, 45));
        assert_eq!(row_use, RowUse::Placed);
        assert_eq!(spans(&records), [(116, 160, 1, 9)]);
    }

    /// The model's records are kept, and nothing is counted, where the row
    /// gives the positions the records give, where no row is read and where
    /// the transcript has no coding sequence.
    #[test]
    fn a_row_that_agrees_with_the_records_changes_nothing() {
        let forward = exon(101, 160, 1);
        let model = [cds(107, 148, 1, 0, 7)];
        for (written, positions) in [
            (Some((7, 48)), row(7, 48)),
            (Some((7, 48)), None),
            (None, row(19, 60)),
        ] {
            let (records, row_use) = coding_records(&[&forward], &model, 1, written, positions);
            assert_eq!(row_use, RowUse::Agrees);
            assert_eq!(spans(&records), [(107, 148, 0, 7)]);
        }
    }

    /// A row that differs places nothing in a transcript of two exons, nor
    /// in one exon where its two exon identifiers differ, where it ends past
    /// the exon's 60 bases, or where its end comes before its start: the
    /// records stay the model's and the transcript is counted.
    #[test]
    fn a_row_whose_exons_are_not_known_leaves_the_records_as_the_model_has_them() {
        let first = exon(101, 160, 1);
        let second = exon(201, 260, 1);
        let model = [cds(101, 160, 1, 0, 7)];
        let two_ids = Some(ExonPositions {
            start: 19,
            end: 60,
            in_one_exon: false,
        });
        let cases: [(&[&ExonRecord], Option<ExonPositions>); 4] = [
            (&[&first, &second], row(19, 60)),
            (&[&first], two_ids),
            (&[&first], row(19, 61)),
            (&[&first], row(40, 39)),
        ];
        for (ordered, positions) in cases {
            let (records, row_use) = coding_records(ordered, &model, 1, Some((1, 60)), positions);
            assert_eq!(row_use, RowUse::Unplaced, "{positions:?}");
            assert_eq!(spans(&records), [(101, 160, 0, 7)]);
        }
    }
}
