// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The rows of one input record.

use crate::entry::{Entry, Fields};
use crate::error::Result;
use crate::fields::record_label;
use crate::shared::{perl_split, replaced_reference, vep_rs_pos, UPLOADED_ALLELES_RE};
use indexmap::IndexMap;
use std::collections::HashMap;

/// One input line: the place of its name in `lines_by_name` and its place among the lines
/// of that name.
type Line = (usize, usize);

/// The entries of one output grouped by input record (`record_label`): the row group
/// `compare_entries` hands a field class whose shape spans a record's rows
/// (`RECORD_READING_CLASSES`). `input_lines` are the input VCF's data lines, read for the
/// default and tab rows that name their record by an identifier. A line carries the first
/// `;`-joined element of its ID column, the name both engines print for its rows, and no
/// name where that element is a dot, empty or `0`.
///
/// The entries are held by their place in the slice the rows were built from, and an entry
/// is known again by its address in that slice. A line is held as its first five columns,
/// the ones a row is read against.
#[derive(Debug)]
pub struct RecordRows<'a> {
    entries: &'a [Entry],
    /// The places of the entries of each label, in the entries' order.
    pub by_record: IndexMap<String, Vec<usize>>,
    /// The input lines by identifier, each stored with the identifier as its third column.
    pub lines_by_name: IndexMap<String, Vec<Vec<String>>>,
    /// The line of each entry that several lines were asked about, by the entry's place:
    /// asking a label's entries for their lines again at every look-up costs the square of
    /// its entries times the lines.
    _line_by_entry: HashMap<usize, Option<Line>>,
    /// Per label, the places of its entries under each line, in their order.
    _entries_by_line: HashMap<String, HashMap<Line, Vec<usize>>>,
    /// The entries that stand in for a record's own, one a line and Location.
    record_entries: Vec<Entry>,
    record_entry_at: HashMap<(Line, String), usize>,
}

impl<'a> RecordRows<'a> {
    pub fn new(entries: &'a [Entry], input_lines: Option<&[String]>) -> Self {
        let mut by_record: IndexMap<String, Vec<usize>> = IndexMap::new();
        for (at, e) in entries.iter().enumerate() {
            by_record.entry(record_label(e)).or_default().push(at);
        }
        let mut lines_by_name: IndexMap<String, Vec<Vec<String>>> = IndexMap::new();
        for line in input_lines.unwrap_or_default() {
            let mut columns: Vec<&str> = line.splitn(6, '\t').collect();
            if columns.len() < 5 {
                continue;
            }
            columns.truncate(5);
            columns[2] = columns[2].split(';').next().unwrap_or("");
            if !matches!(columns[2], "." | "" | "0") {
                lines_by_name
                    .entry(columns[2].to_string())
                    .or_default()
                    .push(columns.into_iter().map(str::to_string).collect());
            }
        }
        RecordRows {
            entries,
            by_record,
            lines_by_name,
            _line_by_entry: HashMap::new(),
            _entries_by_line: HashMap::new(),
            record_entries: Vec::new(),
            record_entry_at: HashMap::new(),
        }
    }

    /// The entries of `label`, in their order.
    pub fn rows(&self, label: &str) -> Vec<&'a Entry> {
        let entries = self.entries;
        self.by_record
            .get(label)
            .map(|rows| rows.iter().map(|row| &entries[*row]).collect())
            .unwrap_or_default()
    }

    /// The place of `e` among the entries the rows were built from; `None` for any other
    /// entry, an equal one included.
    fn place(&self, e: &Entry) -> Option<usize> {
        let size = std::mem::size_of::<Entry>();
        let offset = (std::ptr::from_ref(e) as usize).checked_sub(self.entries.as_ptr() as usize)?;
        (offset % size == 0 && offset / size < self.entries.len()).then_some(offset / size)
    }

    fn line_at(&mut self, e: &Entry) -> Result<Option<Line>> {
        let name = e
            .fields
            .get("Uploaded_variation")
            .map_or("", String::as_str);
        let Some((name_at, _, lines)) = self.lines_by_name.get_full(name) else {
            return Ok(None);
        };
        if lines.len() < 2 {
            return Ok((!lines.is_empty()).then_some((name_at, 0)));
        }
        let place = self.place(e);
        if let Some(asked) = place.and_then(|place| self._line_by_entry.get(&place)) {
            return Ok(*asked);
        }
        let mut writing = Vec::new();
        for (ordinal, columns) in lines.iter().enumerate() {
            if line_writes_row(columns, e)? {
                writing.push(ordinal);
            }
        }
        let line = match writing[..] {
            [ordinal] => Some((name_at, ordinal)),
            _ => None,
        };
        if let Some(place) = place {
            self._line_by_entry.insert(place, line);
        }
        Ok(line)
    }

    /// The input line of the record a default or tab row under an identifier belongs to:
    /// the one line that carries the identifier, or of several the one that writes the row
    /// at its Location with its allele (`line_writes_row`). `None` for an identifier no
    /// line carries, and where none or more than one of several lines writes the row. Several
    /// lines are asked about an entry of this output once, when its own line is wanted.
    pub fn line_of(&mut self, e: &Entry) -> Result<Option<Vec<String>>> {
        Ok(self
            .line_at(e)?
            .map(|(name_at, ordinal)| self.lines_by_name[name_at][ordinal].clone()))
    }

    fn written_at(&mut self, label: &str, line: Line) -> Result<Vec<usize>> {
        if !self._entries_by_line.contains_key(label) {
            let entries = self.entries;
            let rows = self.by_record.get(label).cloned().unwrap_or_default();
            let mut by_line: HashMap<Line, Vec<usize>> = HashMap::new();
            for row in rows {
                if let Some(written) = self.line_at(&entries[row])? {
                    by_line.entry(written).or_default().push(row);
                }
            }
            self._entries_by_line.insert(label.to_string(), by_line);
        }
        Ok(self._entries_by_line[label]
            .get(&line)
            .cloned()
            .unwrap_or_default())
    }

    /// The entries of `label` whose line (`line_of`) is `line`, in the label's order.
    pub fn _written_by(&mut self, label: &str, line: &[String]) -> Result<Vec<&'a Entry>> {
        let held = line.get(2).and_then(|name| {
            let (name_at, _, lines) = self.lines_by_name.get_full(name)?;
            let ordinal = lines.iter().position(|columns| columns == line)?;
            Some((name_at, ordinal))
        });
        let Some(held) = held else {
            return Ok(Vec::new());
        };
        let entries = self.entries;
        Ok(self
            .written_at(label, held)?
            .into_iter()
            .map(|row| &entries[row])
            .collect())
    }

    /// Whether `e`'s label names one input record: a VCF or JSON entry's ordinal, a
    /// default or tab row whose identifier-less `Uploaded_variation` spells the record's
    /// alleles, or one whose identifier names one line of the input (`line_of`). Under an
    /// identifier two records at one Location share a label, so without the input their
    /// rows do not tell a multi-allelic record from two one-allele records.
    pub fn one_record(&mut self, e: &Entry) -> Result<bool> {
        if e.location.is_none() {
            return Ok(true);
        }
        let uploaded = e
            .fields
            .get("Uploaded_variation")
            .map_or("", String::as_str);
        // `uploaded_alleles(uploaded) is not None`.
        Ok(UPLOADED_ALLELES_RE.is_match(uploaded) || self.line_at(e)?.is_some())
    }

    /// The entries of `e`'s input record, or `None` where the label may name two
    /// records (`one_record`). The record's own entry is among them in the VCF and JSON
    /// formats, and for a default or tab row under an identifier one built from its input
    /// line (`CHROM`, `POS`, `REF`, `ALT`) stands in for it. Where several lines
    /// carry the identifier, the rows are those of the label that the same line writes.
    pub fn of(&mut self, e: &Entry) -> Result<Option<Vec<&Entry>>> {
        if !self.one_record(e)? {
            return Ok(None);
        }
        let label = record_label(e);
        let (Some(location), Some(line)) = (
            e.location.as_deref(),
            match e.location {
                Some(_) => self.line_at(e)?,
                None => None,
            },
        ) else {
            return Ok(Some(self.rows(&label)));
        };
        let rows = if self.lines_by_name[line.0].len() > 1 {
            self.written_at(&label, line)?
        } else {
            self.by_record.get(&label).cloned().unwrap_or_default()
        };
        let key = (line, location.to_string());
        let record_entry = match self.record_entry_at.get(&key) {
            Some(held) => *held,
            None => {
                let columns = &self.lines_by_name[line.0][line.1];
                let fields: Fields = [
                    ("Uploaded_variation", &columns[2]),
                    ("CHROM", &columns[0]),
                    ("POS", &columns[1]),
                    ("REF", &columns[3]),
                    ("ALT", &columns[4]),
                ]
                .into_iter()
                .map(|(name, value)| (name.to_string(), value.clone()))
                .collect();
                self.record_entries.push(Entry {
                    key: format!("{location}|record"),
                    location: Some(location.to_string()),
                    fields,
                    ..Entry::default()
                });
                self.record_entry_at
                    .insert(key, self.record_entries.len() - 1);
                self.record_entries.len() - 1
            }
        };
        let entries = self.entries;
        let mut found: Vec<&Entry> = rows.into_iter().map(|row| &entries[row]).collect();
        found.push(&self.record_entries[record_entry]);
        Ok(Some(found))
    }
}

#[cfg(test)]
thread_local! {
    /// The ALT column of every line `line_writes_row` was asked about on this thread.
    static ASKED: std::cell::RefCell<Vec<String>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Whether the input line `columns` (`CHROM`, `POS`, `ID`, `REF`, `ALT`)
/// writes the default or tab row `e`: it stands on the row's chromosome, a `chr` prefix
/// set aside, and one of its alternates is its reference allele with the row's allele in
/// place of the Location's bases, or between the Location's two positions
/// (`replaced_reference`). Each alternate is asked by itself, so a line two of whose
/// alternates give the row in two ways writes it too.
pub fn line_writes_row(columns: &[String], e: &Entry) -> Result<bool> {
    #[cfg(test)]
    ASKED.with(|asked| asked.borrow_mut().extend(columns.get(4).cloned()));
    let Some(location) = e.location.as_deref() else {
        return Ok(false);
    };
    let column = |at: usize| -> &str {
        match columns.get(at) {
            Some(column) => column,
            None => crate::uncaught!("list index out of range"),
        }
    };
    let unprefixed = |name| str::strip_prefix(name, "chr").unwrap_or(name);
    let Some(position) = vep_rs_pos(column(1)) else {
        return Ok(false);
    };
    if unprefixed(column(0)) != unprefixed(location.split(':').next().unwrap_or("")) {
        return Ok(false);
    }
    for alternate in perl_split(column(4), ",") {
        let replaced = replaced_reference(
            position,
            column(3),
            &[alternate.as_str()],
            Some(location),
            &e.allele,
        )?;
        if replaced.is_some() {
            return Ok(true);
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(location: &str, allele: &str, feature: &str, name: &str) -> Entry {
        Entry {
            key: format!("{location}|{allele}|{feature}|Transcript"),
            location: Some(location.to_string()),
            allele: allele.to_string(),
            feature: feature.to_string(),
            fields: [("Uploaded_variation".to_string(), name.to_string())]
                .into_iter()
                .collect(),
            ..Entry::default()
        }
    }

    #[test]
    fn an_input_line_stands_under_the_first_identifier_of_its_id_column_alone() {
        let names = [
            "rs1;rs2", "rs3", ".;rs4", ";rs5", "0;rs6", "0", "", ".", "rs1",
        ];
        let lines: Vec<String> = names
            .iter()
            .enumerate()
            .map(|(at, name)| format!("21\t{}\t{name}\tA\tT\t.\t.\t.", 100 + at))
            .collect();
        let records = RecordRows::new(&[], Some(&lines));
        let found: Vec<(&str, Vec<&[String]>)> = records
            .lines_by_name
            .iter()
            .map(|(name, found)| (name.as_str(), found.iter().map(|line| &line[..3]).collect()))
            .collect();
        let line = |position: &str, name: &str| {
            vec!["21".to_string(), position.to_string(), name.to_string()]
        };
        let (first, last, other) = (line("100", "rs1"), line("108", "rs1"), line("101", "rs3"));
        assert_eq!(
            found,
            [
                ("rs1", vec![&first[..], &last[..]]),
                ("rs3", vec![&other[..]]),
            ]
        );
    }

    /// What `RecordRows.of` hands for a row, each entry as its allele, its feature and its
    /// `ALT` field.
    type Handed = Option<Vec<(String, String, Option<String>)>>;

    #[test]
    fn each_row_under_several_lines_is_asked_for_its_line_once() {
        let (location, features) = ("21:1000", ["ENST1", "ENST2", "ENST3", "ENST4"]);
        let carried = [
            ("rs1", "C"),
            ("rs1", "G"),
            ("rs1", "T"),
            ("rs2", "C"),
            ("rs2", "G"),
            ("rs2", "G"),
            ("rs3", "T"),
        ];
        let lines: Vec<String> = carried
            .iter()
            .map(|(name, alternate)| format!("21\t1000\t{name}\tA\t{alternate}\t.\t.\t."))
            .collect();
        let written = [
            ("rs1", "C"),
            ("rs2", "G"),
            ("rs1", "G"),
            ("rs3", "C"),
            ("rs2", "C"),
            ("rs1", "T"),
            ("rs2", "T"),
            ("rs3", "T"),
        ];
        let entries: Vec<Entry> = features
            .iter()
            .flat_map(|feature| {
                written
                    .iter()
                    .map(move |(name, allele)| row(location, allele, feature, name))
            })
            .collect();
        let put = || ASKED.with(|asked| std::mem::take(&mut *asked.borrow_mut()));
        put();
        assert!(RecordRows::new(&entries, Some(&lines))
            .of(&entries[1])
            .unwrap()
            .is_none());
        assert_eq!(put(), ["C", "G", "G"]);

        // A scan of the label's rows for those the row's own line writes.
        let scan = RecordRows::new(&entries, Some(&lines));
        let line = |e: &Entry| -> Option<Vec<String>> {
            let carrying = &scan.lines_by_name[&e.fields["Uploaded_variation"]];
            let writing: Vec<&Vec<String>> = carrying
                .iter()
                .filter(|columns| line_writes_row(columns, e).unwrap())
                .collect();
            match (carrying.len(), writing.len()) {
                (1, _) => Some(carrying[0].clone()),
                (_, 1) => Some(writing[0].clone()),
                _ => None,
            }
        };
        let scanned = |e: &Entry| -> Handed {
            let own = line(e)?;
            let mut rows: Vec<(String, String, Option<String>)> = scan
                .rows(&record_label(e))
                .into_iter()
                .filter(|row| line(row).as_ref() == Some(&own))
                .map(|row| (row.allele.clone(), row.feature.clone(), None))
                .collect();
            rows.push((String::new(), String::new(), Some(own[4].clone())));
            Some(rows)
        };
        let expected: Vec<Handed> = entries.iter().map(scanned).collect();
        put();

        let mut records = RecordRows::new(&entries, Some(&lines));
        for _ in 0..2 {
            for (e, expected) in entries.iter().zip(&expected) {
                let handed: Handed = records.of(e).unwrap().map(|rows| {
                    rows.into_iter()
                        .map(|row| {
                            (
                                row.allele.clone(),
                                row.feature.clone(),
                                row.fields.get("ALT").cloned(),
                            )
                        })
                        .collect()
                });
                assert_eq!(
                    &handed,
                    expected,
                    "{:?}",
                    (&e.fields, &e.allele, &e.feature)
                );
            }
        }
        let asked = put();
        let rows_of = |alleles: &[&str]| -> Vec<(String, String, Option<String>)> {
            features
                .iter()
                .flat_map(|feature| {
                    alleles
                        .iter()
                        .map(move |allele| (allele.to_string(), feature.to_string(), None))
                })
                .collect()
        };
        let with_line = |mut rows: Vec<(String, String, Option<String>)>, alternate: &str| {
            rows.push((String::new(), String::new(), Some(alternate.to_string())));
            Some(rows)
        };
        assert_eq!(expected[0], with_line(rows_of(&["C"]), "C"));
        assert!(expected[1].is_none() && expected[6].is_none());
        assert_eq!(expected[3], with_line(rows_of(&["C", "T"]), "T"));
        let under_several_lines = entries
            .iter()
            .filter(|e| matches!(e.fields["Uploaded_variation"].as_str(), "rs1" | "rs2"))
            .count();
        assert_eq!(asked.len(), 3 * under_several_lines);
        assert_eq!(asked.len(), 72);
    }
}
