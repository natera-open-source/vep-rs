// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The REVEL plugin's classes and the score table they read.
//!
//! A table row scores one substitution, from its `aaref` to its `aaalt`, on the transcripts
//! it lists. REVEL.pm:170-199 walk the rows at a record's span and return the first whose
//! alleles match one of the record's alternates, whose transcript ids hold the output row's
//! and whose `aaalt` is the output row's peptide; the `aaref` cell is read and never tested.

use super::entry::{blank, field_of, uploaded_alleles, uploaded_position, Entry, Row};
use super::shared::{record_entry, vcf_alleles_as_spelt};
use super::sibling::AMINO_ACIDS_FIELDS;
use super::tabix::{IndexedTable, TableError, TableRow};
use super::{Refusal, Res};
use crate::pytext;
use std::collections::BTreeSet;
use std::path::Path;
use std::sync::Arc;

/// One row of a REVEL table as REVEL.pm:207-241 read it: the cell of the run's position
/// column as written, the two allele cells (`None` for a cell at or past the row's last
/// filled one, the undefined value `split /\t/` leaves), the two amino-acid cells, the score
/// cell and the transcript ids cell (`None` likewise, and for the seven-column layout, which
/// has no ids).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RevelRow<'a> {
    pub position: &'a str,
    pub reference: Option<&'a str>,
    pub alternate: Option<&'a str>,
    pub reference_amino_acid: &'a str,
    pub alternate_amino_acid: &'a str,
    pub score: Option<&'a str>,
    pub transcript_ids: Option<&'a str>,
}

/// The REVEL table a run's plugin read, read through its tabix index as the plugin reads
/// it: the rows at a span of positions, in file order. The header, the last of the file's
/// leading `#` lines, fixes the layout by its column count (REVEL.pm:127-130, :209-241):
/// seven columns are `chr, hg19_pos, ref, alt, aaref, aaalt, REVEL`; eight and nine are
/// `chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL, Ensembl_transcriptid`. The
/// position column is the one the index was built on, which is the run's assembly's.
pub struct RevelTable {
    table: IndexedTable,
    /// The columns of the reference allele, the alternate allele, the two amino acids, the
    /// score and the transcript ids.
    columns: [Option<usize>; 6],
    no_match: bool,
}

fn refusal(error: TableError) -> Refusal {
    Refusal(format!("--revel-table: {}", error.0))
}

impl RevelTable {
    /// The table at `path` with its tabix index beside it; `no_match` says the run
    /// set the plugin's `no_match` (REVEL.pm:107-110), under which no row's transcript ids
    /// are read. A file that is no indexed REVEL table is refused by name.
    pub fn new(path: &Path, no_match: bool) -> Result<RevelTable, Refusal> {
        let table = IndexedTable::open(path).map_err(refusal)?;
        let columns = match table.header().len() {
            7 => [Some(2), Some(3), Some(4), Some(5), Some(6), None],
            8 | 9 => [Some(3), Some(4), Some(5), Some(6), Some(7), Some(8)],
            _ => {
                return Err(Refusal(format!(
                    "--revel-table: {} is not a REVEL table (a header of 7, 8 or 9 columns)",
                    path.display()
                )))
            }
        };
        Ok(RevelTable {
            table,
            columns,
            no_match,
        })
    }

    /// Whether the run set the plugin's `no_match`.
    pub fn no_match(&self) -> bool {
        self.no_match
    }

    /// The table's rows at the positions `start` to `end` of `chromosome`, in file order:
    /// those REVEL.pm:168 fetches for a variant of that span. The chromosome is looked up
    /// as given, then without and with a `chr` prefix.
    pub fn rows(&self, chromosome: &str, start: i128, end: i128) -> Res<Arc<Vec<TableRow>>> {
        self.table.rows(chromosome, start, end).map_err(refusal)
    }

    /// The cells of a table row under the table's layout.
    pub fn row<'a>(&self, row: &'a TableRow) -> RevelRow<'a> {
        let cell = |at: usize| self.columns[at].and_then(|column| row.cell(column));
        RevelRow {
            position: &row.position,
            reference: cell(0),
            alternate: cell(1),
            reference_amino_acid: cell(2).unwrap_or(""),
            alternate_amino_acid: cell(3).unwrap_or(""),
            score: cell(4),
            transcript_ids: cell(5),
        }
    }
}

/// A reference allele, an alternate allele and the position of the reference's first base.
pub type AlleleKey = (String, String, i128);

/// The (reference, alternate, start) of an allele pair after each trimming order
/// `get_matched_variant_alleles` tries (ensembl-variation Utils/Sequence.pm:1190-1209 and
/// :1250-1258, `trim_sequences` :965-1038 on the forward strand): the shared bases cut from
/// the left and then from the right, and from the right first as well where either allele is
/// longer than one base; an allele left empty is written `-`.
pub fn trimmed_allele_keys(
    reference: &str,
    alternate: &str,
    position: i128,
) -> BTreeSet<AlleleKey> {
    let reference: Vec<char> = reference.chars().collect();
    let alternate: Vec<char> = alternate.chars().collect();
    let orders: &[bool] = if reference.len() > 1 || alternate.len() > 1 {
        &[false, true]
    } else {
        &[false]
    };
    let spelt = |allele: &[char]| {
        if allele.is_empty() {
            "-".to_string()
        } else {
            allele.iter().collect()
        }
    };
    let mut keys = BTreeSet::new();
    for end_first in orders {
        let (mut reference, mut alternate, mut start) = (&reference[..], &alternate[..], position);
        for left in [!end_first, *end_first] {
            while !reference.is_empty() && !alternate.is_empty() {
                if left && reference[0] == alternate[0] {
                    // A position past the largest `i128` wraps: the keys are only told apart,
                    // and no two positions of one record wrap to one.
                    (reference, alternate, start) =
                        (&reference[1..], &alternate[1..], start.wrapping_add(1));
                } else if !left && reference.last() == alternate.last() {
                    reference = &reference[..reference.len() - 1];
                    alternate = &alternate[..alternate.len() - 1];
                } else {
                    break;
                }
            }
        }
        keys.insert((spelt(reference), spelt(alternate), start));
    }
    keys
}

/// The integer of a text `str.isdigit()` takes; `None` for any other text, and for digits
/// `int()` does not read or a number outside `i128`.
fn digits_value(text: &str) -> Option<i128> {
    if pytext::isdigit(text) {
        pytext::int(text)
    } else {
        None
    }
}

/// A record as the plugin's match reads it: its chromosome, the start and end of its span,
/// the allele keys of each of its alternates and the place of a row's own alternate among
/// them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RevelRecord {
    pub chromosome: String,
    pub start: i128,
    pub end: i128,
    pub keys: Vec<BTreeSet<AlleleKey>>,
    pub place: usize,
}

/// The record of the row `e` as the plugin's match reads it (REVEL.pm:168-188). The alleles
/// are read as Ensembl VEP spells them, an indel record's shared first base dropped and its
/// start one base on (`vcf_alleles_as_spelt`): from an identifier-less `Uploaded_variation`,
/// else from the record entry among `record_rows` under `e`'s own label, a VCF line's
/// `CHROM`, `POS`, `REF` and `ALT` or the JSON `seq_region_name`, `start` and
/// `allele_string`. `None` where neither states the record, where a VCF line's start lies
/// past the largest `i128`, and where no alternate, or more than one, spells the row's allele.
pub fn revel_record(e: &Entry, record_rows: Option<&[Row]>) -> Res<Option<RevelRecord>> {
    let name = e.fields.get_or("Uploaded_variation", "");
    let stated: Option<(String, Option<i128>, String, Vec<String>)> =
        match (uploaded_alleles(name), e.location.as_deref()) {
            (Some((reference, alternates)), Some(location)) => Some((
                pytext::rpartition(location, ":").0.to_string(),
                uploaded_position(name)?,
                reference,
                alternates,
            )),
            _ => {
                let record = record_entry(e, record_rows.unwrap_or(&[]));
                let field = |name: &str| record.and_then(|record| record.fields.get(name));
                let line = (field("CHROM"), field("POS"), field("REF"), field("ALT"));
                let object = (
                    field("seq_region_name"),
                    field("start"),
                    field("allele_string"),
                );
                match (line, object) {
                    ((Some(chromosome), Some(position), Some(reference), Some(alternates)), _)
                        if pytext::isdigit(position) =>
                    {
                        let alternates: Vec<&str> = alternates.split(',').collect();
                        let (spelt, alternates) = vcf_alleles_as_spelt(reference, &alternates);
                        let moved = i128::from(spelt != reference);
                        // A start past the largest `i128` states no record.
                        let start =
                            pytext::int(position).and_then(|position| position.checked_add(moved));
                        Some((chromosome.to_string(), start, spelt, alternates))
                    }
                    (_, (Some(chromosome), Some(start), Some(alleles)))
                        if pytext::isdigit(start) =>
                    {
                        let mut alleles = alleles.split('/').map(str::to_string);
                        let reference = alleles.next().unwrap_or_default();
                        Some((
                            chromosome.to_string(),
                            pytext::int(start),
                            reference,
                            alleles.collect(),
                        ))
                    }
                    _ => None,
                }
            }
        };
    let Some((chromosome, Some(start), reference, alternates)) = stated else {
        return Ok(None);
    };
    let named = alternates.iter().filter(|a| **a == e.allele).count();
    let place = alternates.iter().position(|a| *a == e.allele);
    let (false, 1, Some(place)) = (chromosome.is_empty(), named, place) else {
        return Ok(None);
    };
    let keys = alternates
        .iter()
        .map(|alternate| trimmed_allele_keys(&reference, alternate, start))
        .collect();
    let length = reference.trim_matches('-').chars().count() as i128;
    Ok(Some(RevelRecord {
        chromosome,
        start,
        // No table line stands past the largest `i128`, so the span ends there.
        end: start.saturating_add(length - 1),
        keys,
        place,
    }))
}

/// The keys `get_matched_variant_alleles` gives a table row (Utils/Sequence.pm:1219-1232):
/// its alleles trimmed as `trimmed_allele_keys` trims them, with its position cell as text
/// where no base is cut from the left. `None` for a row with no reference or alternate
/// allele cell, which the match dies on (Utils/Sequence.pm:968); no key for a row whose
/// position cell is false to Perl.
pub fn revel_row_keys(row: &RevelRow) -> Option<BTreeSet<(String, String, String)>> {
    let (reference, alternate) = (row.reference?, row.alternate?);
    if matches!(row.position, "" | "0") {
        return Some(BTreeSet::new());
    }
    let spelt = |allele: &str| if allele.is_empty() { "-" } else { allele }.to_string();
    let Some(start) = digits_value(row.position) else {
        return Some(BTreeSet::from([(
            spelt(reference),
            spelt(alternate),
            row.position.to_string(),
        )]));
    };
    Some(
        trimmed_allele_keys(reference, alternate, start)
            .into_iter()
            .map(|(reference, alternate, at)| {
                let position = if at == start {
                    row.position.to_string()
                } else {
                    at.to_string()
                };
                (reference, alternate, position)
            })
            .collect(),
    )
}

/// The tests REVEL.pm:191-197 make of a row whose alleles matched: where its transcript ids
/// cell is true to Perl and the run did not set `no_match`, one of its `;`-separated ids is
/// the row's `Feature`; and its alternate amino acid is the row's.
pub fn revel_row_passes(
    table: &RevelTable,
    row: &RevelRow,
    feature: &str,
    alternate_residue: &str,
) -> bool {
    if let Some(ids) = row.transcript_ids.filter(|ids| !matches!(*ids, "" | "0")) {
        if !table.no_match && !ids.split(';').any(|id| id == feature) {
            return false;
        }
    }
    row.alternate_amino_acid == alternate_residue
}

/// The first of `rows`, the table's rows over a record's span in file order, that carries
/// one of the allele keys `keys` and passes `revel_row_passes`: the row whose score
/// REVEL.pm:170-199 return for an output row. With `reference_residue` the row's reference
/// amino acid is tested as well. `None` where the walk ends without such a row, and where
/// it meets a row with no allele cell first, on which the plugin's match dies.
pub fn revel_walk<'a>(
    table: &RevelTable,
    rows: &'a [TableRow],
    keys: &BTreeSet<AlleleKey>,
    feature: &str,
    alternate_residue: &str,
    reference_residue: Option<&str>,
) -> Option<RevelRow<'a>> {
    let wanted: BTreeSet<(String, String, String)> = keys
        .iter()
        .map(|(reference, alternate, at)| (reference.clone(), alternate.clone(), at.to_string()))
        .collect();
    for row in rows {
        let row = table.row(row);
        let row_keys = revel_row_keys(&row)?;
        if row_keys.is_disjoint(&wanted)
            || !revel_row_passes(table, &row, feature, alternate_residue)
        {
            continue;
        }
        if reference_residue.is_none_or(|residue| row.reference_amino_acid == residue) {
            return Some(row);
        }
    }
    None
}

/// The chromosome and position of a row whose Location is one position and whose allele is
/// one base, read from the row's Location or from the Location its record gives its allele
/// (`Entry.allele_location`); `None` for any other row.
pub fn revel_one_base_site(e: &Entry) -> Option<(String, i128)> {
    let location = e.location.as_deref().or(e.allele_location.as_deref())?;
    if !matches!(e.allele.as_str(), "A" | "C" | "G" | "T") {
        return None;
    }
    let (chromosome, _, position) = pytext::rpartition(location, ":");
    if chromosome.is_empty() {
        return None;
    }
    Some((chromosome.to_string(), digits_value(position)?))
}

/// Whether the reference's `REVEL` on the row `e` is the score of a table row whose
/// reference amino acid is not the row's own, and the engine's is what the walk that tests
/// both cells gives. Both outputs spell the row's `Amino_acids` alike, as one pair of
/// residues.
///
/// Where the record is stated (`revel_record`), the reference's value must be the score of
/// the row its walk reaches over every alternate of the record (`revel_walk`), and that
/// row's reference amino acid must be another than the first residue of `Amino_acids`.
/// Where it is not, a one-base row at one position (`revel_one_base_site`) is read without
/// it: the rows the reference's walk can reach, whatever other alternates the record holds,
/// are the first passing row of the row's own allele and, for each other allele, its first
/// passing row before that one; every one of them that carries the reference's value must
/// name another reference amino acid, and one must carry it.
///
/// The engine's value must be the score of the first row of the row's own allele that
/// passes and names both residues of `Amino_acids`, or be empty where there is none.
/// Nothing is accepted without the table, without the engine's row, on a row whose walk
/// meets a table row with no allele cell, or on a value no such table row explains.
pub fn revel_other_reference_amino_acid_row(
    e: &Entry,
    a: Option<&Entry>,
    expected_value: &str,
    actual_value: &str,
    record_rows: Option<&[Row]>,
    table: Option<&RevelTable>,
) -> Res<bool> {
    let (Some(table), Some(a)) = (table, a) else {
        return Ok(false);
    };
    if blank(expected_value) {
        return Ok(false);
    }
    let amino_acids = field_of(e, &AMINO_ACIDS_FIELDS);
    let (reference_residue, slash, alternate_residue) = pytext::partition(amino_acids, "/");
    if slash.is_empty() || reference_residue.is_empty() || alternate_residue.is_empty() {
        return Ok(false);
    }
    if field_of(a, &AMINO_ACIDS_FIELDS) != amino_acids {
        return Ok(false);
    }
    let rows;
    let reached: Vec<RevelRow>;
    let own: Option<RevelRow>;
    if let Some(record) = revel_record(e, record_rows)? {
        rows = table.rows(&record.chromosome, record.start, record.end)?;
        let every: BTreeSet<AlleleKey> = record.keys.iter().flatten().cloned().collect();
        let walk = |keys: &BTreeSet<AlleleKey>, residue: Option<&str>| {
            revel_walk(table, &rows, keys, &e.feature, alternate_residue, residue)
        };
        reached = walk(&every, None).into_iter().collect();
        own = walk(&record.keys[record.place], Some(reference_residue));
    } else {
        let Some((chromosome, position)) = revel_one_base_site(e) else {
            return Ok(false);
        };
        rows = table.rows(&chromosome, position, position)?;
        let read: Vec<RevelRow> = rows.iter().map(|row| table.row(row)).collect();
        if read.iter().any(|row| revel_row_keys(row).is_none()) {
            return Ok(false);
        }
        let at = position.to_string();
        let passing: Vec<RevelRow> = read
            .into_iter()
            .filter(|row| {
                row.position == at
                    && row.alternate.unwrap_or("").chars().count() == 1
                    && revel_row_passes(table, row, &e.feature, alternate_residue)
            })
            .collect();
        let own_allele = |row: &RevelRow| row.alternate == Some(e.allele.as_str());
        let first_own = passing.iter().position(own_allele).unwrap_or(passing.len());
        let mut first_rows: Vec<RevelRow> = Vec::new();
        for row in passing.iter().take(first_own + 1) {
            if first_rows
                .iter()
                .all(|earlier| row.alternate != earlier.alternate)
            {
                first_rows.push(*row);
            }
        }
        reached = first_rows;
        own = passing
            .iter()
            .find(|row| own_allele(row) && row.reference_amino_acid == reference_residue)
            .copied();
    }
    let mut explaining = reached
        .iter()
        .filter(|row| row.score == Some(expected_value))
        .peekable();
    if explaining.peek().is_none()
        || explaining.any(|row| row.reference_amino_acid == reference_residue)
    {
        return Ok(false);
    }
    let wanted = own.and_then(|row| row.score).unwrap_or("");
    Ok(if blank(wanted) {
        blank(actual_value)
    } else {
        actual_value == wanted
    })
}
